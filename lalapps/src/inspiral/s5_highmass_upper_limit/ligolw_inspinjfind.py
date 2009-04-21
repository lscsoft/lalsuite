# $Id: ligolw_binjfind.py,v 1.44 2009/03/08 01:40:04 kipp Exp $
#
# Copyright (C) 2006  Kipp C. Cannon
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#


"""
Inspiral injection identification library.  Contains code providing the
capacity to search a list of sngl_inspiral candidates for events
matching entries in a sim_inspiral list of software injections, recording the
matches as inspiral <--> injection coincidences using the standard coincidence
infrastructure.  Also, any pre-recorded inspiral <--> inspiral coincidences are
checked for cases where all the inspiral events in a coincidence match an
injection, and these are recorded as coinc <--> injection coincidences,
again using the standard coincidence infrastructure.
"""


import bisect
import sys


from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw.utils import process as ligolw_process
from pylal import ligolw_thinca
from pylal import llwapp
from pylal import SimInspiralUtils
from pylal import SnglInspiralUtils
from pylal.date import LIGOTimeGPS
from pylal.xlal import tools


lsctables.CoincMapTable.RowType = lsctables.CoincMap = tools.CoincMap


__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__version__ = "$Revision: 1.44 $"[11:-2]
__date__ = "$Date: 2009/03/08 01:40:04 $"[7:-2]


#
# =============================================================================
#
#                                 Speed Hacks
#
# =============================================================================
#


lsctables.LIGOTimeGPS = LIGOTimeGPS


def sngl_inspiral___cmp__(self, other):
	# compare self's end time to the LIGOTimeGPS instance other
	return cmp(self.end_time, other.seconds) or cmp(self.end_time_ns, other.nanoseconds)


lsctables.SnglInspiral.__cmp__ = sngl_inspiral___cmp__


#
# =============================================================================
#
#                              Document Interface
#
# =============================================================================
#


#FIXME This is in ligolw_thinca already
#InspiralCoincDef = lsctables.CoincDef(search = u"inspiral", search_coinc_type = 0, description = u"sngl_inspiral<-->sngl_inspiral coincidences")

InspiralSICoincDef = lsctables.CoincDef(search = u"inspiral", search_coinc_type = 1, description = u"sim_inspiral<-->sngl_inspiral coincidences")
InspiralSCNearCoincDef = lsctables.CoincDef(search = u"inspiral", search_coinc_type = 2, description = u"sim_inspiral<-->coinc_event coincidences (nearby)")

class DocContents(object):
	"""
	A wrapper interface to the XML document.
	"""
	def __init__(self, xmldoc, bbdef, sbdef, scndef, process):
		#
		# store the process row
		#

		self.process = process

		#
		# locate the sngl_inspiral and sim_inspiral tables
		#

		self.snglinspiraltable = table.get_table(xmldoc, lsctables.SnglInspiralTable.tableName)
		self.siminspiraltable = table.get_table(xmldoc, lsctables.SimInspiralTable.tableName)

		#
		# construct the zero-lag time slide needed to cover the
		# instruments listed in all the triggers, then determine
		# its ID (or create it if needed)
		#
		# FIXME:  in the future, the sim_inspiral table should
		# indicate time slide at which the injection was done
		#

		self.tisi_id = llwapp.get_time_slide_id(xmldoc, {}.fromkeys(self.snglinspiraltable.getColumnByName("ifo"), 0.0), create_new = process)

		#
		# get coinc_definer row for sim_inspiral <--> sngl_inspiral
		# coincs; this creates a coinc_definer table if the
		# document doesn't have one
		#

		self.sb_coinc_def_id = llwapp.get_coinc_def_id(xmldoc, sbdef.search, sbdef.search_coinc_type, create_new = True, description = sbdef.description)

		#
		# get coinc_def_id's for sngl_inspiral <--> sngl_inspiral, and
		# the sim_inspiral <--> coinc_event coincs.  set all
		# to None if this document does not contain any sngl_inspiral
		# <--> sngl_inspiral coincs.
		#

		try:
			bb_coinc_def_id = llwapp.get_coinc_def_id(xmldoc, bbdef.search, bbdef.search_coinc_type, create_new = False)
		except KeyError:
			bb_coinc_def_id = None
			self.scn_coinc_def_id = None
		else:
			self.scn_coinc_def_id = llwapp.get_coinc_def_id(xmldoc, scndef.search, scndef.search_coinc_type, create_new = True, description = scndef.description)

		#
		# get coinc table, create one if needed
		#

		try:
			self.coinctable = table.get_table(xmldoc, lsctables.CoincTable.tableName)
		except ValueError:
			self.coinctable = lsctables.New(lsctables.CoincTable)
			xmldoc.childNodes[0].appendChild(self.coinctable)
		self.coinctable.sync_next_id()

		#
		# get coinc_map table, create one if needed
		#

		try:
			self.coincmaptable = table.get_table(xmldoc, lsctables.CoincMapTable.tableName)
		except ValueError:
			self.coincmaptable = lsctables.New(lsctables.CoincMapTable)
			xmldoc.childNodes[0].appendChild(self.coincmaptable)

		#
		# index the document
		#
		# FIXME:  inspiral<-->inspiral coincs should be organized by time
		# slide ID, but since injections are only done at zero lag
		# for now this is ignored.
		#

		# index sngl_inspiral table
		index = {}
		for row in self.snglinspiraltable:
			index[row.event_id] = row
		# find IDs of inspiral<-->inspiral coincs
		self.coincs = {}
		for coinc in self.coinctable:
			if coinc.coinc_def_id == bb_coinc_def_id:
				self.coincs[coinc.coinc_event_id] = []
		# construct event list for each inspiral<-->inspiral coinc
		for row in self.coincmaptable:
			if row.coinc_event_id in self.coincs:
				self.coincs[row.coinc_event_id].append(index[row.event_id])
		del index
		# sort each event list by end time and convert to tuples
		# for speed
		for coinc_event_id in self.coincs.keys():
			events = self.coincs[coinc_event_id]
			events.sort(lambda a, b: cmp(a.end_time, b.end_time) or cmp(a.end_time_ns, b.end_time_ns))
			self.coincs[coinc_event_id] = tuple(events)
		# convert dictionary to a list
		self.coincs = self.coincs.items()

		#
		# FIXME Is this true for inspirals too?
		# sort sngl_inspiral table by end time, and sort the coincs
		# list by the end time of the first (earliest) event in
		# each coinc (recall that the event tuple for each coinc
		# has been time-ordered)
		#

		self.snglinspiraltable.sort(lambda a, b: cmp(a.end_time, b.end_time) or cmp(a.end_time_ns, b.end_time_ns))
		self.coincs.sort(lambda (id_a, a), (id_b, b): cmp(a[0].end_time, b[0].end_time) or cmp(a[0].end_time_ns, b[0].end_time_ns))

		#
		# set the window for inspirals_near_endtime().  this window
		# is the amount of time such that if an injection's end
		# time and a inspiral event's end time differ by more than
		# this it is *impossible* for them to match one another.
		#

 		# FIXME I'll just make the windows one second

                self.inspiral_end_time_window = 9.0
                self.coinc_end_time_window = 9.0


	def inspirals_near_endtime(self, t):
		"""
		Return a list of the inspiral events whose peak times are
		within self.inspiral_end_time_window of t.
		"""
		return self.snglinspiraltable[bisect.bisect_left(self.snglinspiraltable, t - self.inspiral_end_time_window):bisect.bisect_right(self.snglinspiraltable, t + self.inspiral_end_time_window)]

	def coincs_near_endtime(self, t):
		"""
		Return a list of the (coinc_event_id, event list) tuples in
		which at least one inspiral event's end time is within
		self.coinc_end_time_window of t.
		"""
		# FIXME:  this test does not consider the time slide
		# offsets that should be applied to the coinc, but for now
		# injections are done at zero lag so this isn't a problem
		# yet
		return [(coinc_event_id, inspirals) for coinc_event_id, inspirals in self.coincs if (t - self.coinc_end_time_window <= inspirals[-1].get_end()) and (inspirals[0].get_end() <= t + self.coinc_end_time_window)]

	def sort_triggers_by_id(self):
		"""
		Sort the sngl_burst table's rows by ID (tidy-up document
		for output).
		"""
		self.snglinspiraltable.sort(lambda a, b: cmp(a.event_id, b.event_id))

	def new_coinc(self, coinc_def_id):
		"""
		Construct a new coinc_event row attached to the given
		process, and belonging to the set of coincidences defined
		by the given coinc_def_id.
		"""
		coinc = lsctables.Coinc()
		coinc.process_id = self.process.process_id
		coinc.coinc_def_id = coinc_def_id
		coinc.coinc_event_id = self.coinctable.get_next_id()
		coinc.time_slide_id = self.tisi_id
		coinc.set_instruments(None)
		coinc.nevents = 0
		coinc.likelihood = None
		self.coinctable.append(coinc)
		return coinc


#
# =============================================================================
#
#                           Add Process Information
#
# =============================================================================
#


process_program_name = "ligolw_inspinjfind"


def append_process(xmldoc, match_algorithm, comment):
	"""
	Convenience wrapper for adding process metadata to the document.
	"""
	process = llwapp.append_process(xmldoc, program = process_program_name, version = __version__, cvs_repository = u"lscsoft", cvs_entry_time = __date__, comment = comment)

	params = [(u"--match-algorithm", u"lstring", match_algorithm)]
	ligolw_process.append_process_params(xmldoc, process, params)

	return process


#
# =============================================================================
#
#                 Injection <--> Inspiral Event Comparison Tests
#
# =============================================================================
#

def InspiralSnglCompare(sim, inspiral):
        """
	Return False if the peak time of the sim is within 9 seconds of the inspiral event.
        """
	return SnglInspiralUtils.CompareSnglInspiral(sim, inspiral, twindow = LIGOTimeGPS(9))


def NearCoincCompare(sim, inspiral):
	"""
	Return False if the peak time of the sim is within 9 seconds of the inspiral event.
	"""
	return SnglInspiralUtils.CompareSnglInspiral(sim, inspiral, twindow = LIGOTimeGPS(9))
	#return not SimBurstUtils.burst_is_near_injection(sim, burst.start_time, burst.start_time_ns, burst.duration, burst.ifo)


#
# =============================================================================
#
#                 Build sim_burst <--> sngl_burst Coincidences
#
# =============================================================================
#


def find_sngl_inspiral_matches(contents, sim, comparefunc):
	"""
	Scan the inspiral table for triggers matching sim.
	"""
	return [inspiral for inspiral in contents.inspirals_near_endtime(sim.get_end()) if not comparefunc(sim, inspiral)]


def add_sim_inspiral_coinc(contents, sim, inspirals):
	"""
	Create a coinc_event in the coinc table, and add arcs in the
	coinc_event_map table linking the sim_inspiral row and the list of
	sngl_inspiral rows to the new coinc_event row.
	"""
	coinc = contents.new_coinc(contents.sb_coinc_def_id)
	coinc.set_instruments(set(event.ifo for event in inspirals))
	coinc.nevents = len(inspirals)

	coincmap = lsctables.CoincMap()
	coincmap.coinc_event_id = coinc.coinc_event_id
	coincmap.table_name = sim.simulation_id.table_name
	coincmap.event_id = sim.simulation_id
	contents.coincmaptable.append(coincmap)

	for event in inspirals:
		coincmap = lsctables.CoincMap()
		coincmap.coinc_event_id = coinc.coinc_event_id
		coincmap.table_name = event.event_id.table_name
		coincmap.event_id = event.event_id
		contents.coincmaptable.append(coincmap)


#
# =============================================================================
#
#                   Build sim_burst <--> coinc Coincidences
#
# =============================================================================
#


def find_coinc_matches(coincs, sim, comparefunc):
	"""
	Return a list of the coinc_event_ids of the inspiral<-->inspiral coincs
	in which all inspiral events match sim.
	"""
	# FIXME:  this test does not consider the time slide offsets that
	# should be applied to the coinc, but for now injections are done
	# at zero lag so this isn't a problem yet
	return [coinc_event_id for coinc_event_id, inspirals in coincs if True not in [bool(comparefunc(sim, inspiral)) for inspiral in inspirals]]


def find_near_coinc_matches(coincs, sim, comparefunc):
	"""
	Return a list of the coinc_event_ids of the inspiral<-->inspiral coincs
	in which at least one inspiral event matches sim.
	"""
	# FIXME:  this test does not consider the time slide offsets that
	# should be applied to the coinc, but for now injections are done
	# at zero lag so this isn't a problem yet
	return [coinc_event_id for coinc_event_id, inspirals in coincs if False in [bool(comparefunc(sim, inspiral)) for inspiral in inspirals]]


def add_sim_coinc_coinc(contents, sim, coinc_event_ids, coinc_def_id):
	"""
	Create a coinc_event in the coinc table, and add arcs in the
	coinc_event_map table linking the sim_inspiral row and the list of
	coinc_event rows to the new coinc_event row.
	"""
	coinc = contents.new_coinc(coinc_def_id)
	coinc.nevents = len(coinc_event_ids)

	coincmap = lsctables.CoincMap()
	coincmap.coinc_event_id = coinc.coinc_event_id
	coincmap.table_name = sim.simulation_id.table_name
	coincmap.event_id = sim.simulation_id
	contents.coincmaptable.append(coincmap)

	for coinc_event_id in coinc_event_ids:
		coincmap = lsctables.CoincMap()
		coincmap.coinc_event_id = coinc.coinc_event_id
		coincmap.table_name = coinc_event_id.table_name
		coincmap.event_id = coinc_event_id
		contents.coincmaptable.append(coincmap)


#
# =============================================================================
#
#                                 Library API
#
# =============================================================================
#


def ligolw_inspinjfind(xmldoc, process, search, snglcomparefunc, nearcoinccomparefunc, verbose = False):
	#
	# Analyze the document's contents.
	#

	if verbose:
		print >>sys.stderr, "indexing ..."

	bbdef = {"inspiral": ligolw_thinca.InspiralCoincDef}[search]
	sbdef = {"inspiral": InspiralSICoincDef}[search]
	scndef = {"inspiral": InspiralSCNearCoincDef}[search]

	contents = DocContents(xmldoc = xmldoc, bbdef = bbdef, sbdef = sbdef, scndef = scndef, process = process)
	N = len(contents.siminspiraltable)

	#
	# Find sim_burst <--> sngl_inspiral coincidences.
	#

	if verbose:
		print >>sys.stderr, "constructing %s:" % sbdef.description
	for n, sim in enumerate(contents.siminspiraltable):
		if verbose:
			print >>sys.stderr, "\t%.1f%%\r" % (100.0 * n / N),
		inspirals = find_sngl_inspiral_matches(contents, sim, snglcomparefunc)
		if inspirals:
			add_sim_inspiral_coinc(contents, sim, inspirals)
	if verbose:
		print >>sys.stderr, "\t100.0%"

	#
	# Find sim_inspiral <--> coinc_event coincidences.
	#

	if contents.scn_coinc_def_id:
		if verbose:
			print >>sys.stderr, "constructing %s:" % (scndef.description)
		for n, sim in enumerate(contents.siminspiraltable):
			if verbose:
				print >>sys.stderr, "\t%.1f%%\r" % (100.0 * n / N),
			coincs = contents.coincs_near_endtime(sim.get_end())
			coinc_event_ids = find_near_coinc_matches(coincs, sim, nearcoinccomparefunc)
			if coinc_event_ids:
				add_sim_coinc_coinc(contents, sim, coinc_event_ids, contents.scn_coinc_def_id)
		if verbose:
			print >>sys.stderr, "\t100.0%"

	#
	# Restore the original event order.
	#

	if verbose:
		print >>sys.stderr, "finishing ..."
	contents.sort_triggers_by_id()

	#
	# Done.
	#

	return xmldoc
