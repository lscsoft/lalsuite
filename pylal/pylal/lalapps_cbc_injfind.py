# Copyright (C) 2006  Kipp Cannon
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


from glue.ligolw import lsctables
from glue.ligolw.utils import coincs as ligolw_coincs
from glue.ligolw.utils import process as ligolw_process
from pylal import git_version
from pylal import ligolw_thinca
from pylal import ligolw_rinca
from lalburst import timeslides as ligolw_tisi
from pylal import SimInspiralUtils
from pylal import SnglInspiralUtils
from pylal.xlal import tools
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS


lsctables.CoincMapTable.RowType = lsctables.CoincMap = tools.CoincMap


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


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

def sngl_ringdown___cmp__(self, other):
        # compare self's end time to the LIGOTimeGPS instance other
        return cmp(self.start_time, other.seconds) or cmp(self.start_time_ns, other.nanoseconds)


lsctables.SnglRingdown.__cmp__ = sngl_ringdown___cmp__

#
# =============================================================================
#
#                              Document Interface
#
# =============================================================================
#


InspiralSICoincDef = lsctables.CoincDef(search = u"inspiral", search_coinc_type = 1, description = u"sim_inspiral<-->sngl_inspiral coincidences")
InspiralSCNearCoincDef = lsctables.CoincDef(search = u"inspiral", search_coinc_type = 2, description = u"sim_inspiral<-->coinc_event coincidences (nearby)")

RingdownSICoincDef = lsctables.CoincDef(search = u"ringdown", search_coinc_type = 1, description = u"sim_ringdown<-->sngl_ringdown coincidences")
RingdownSCNearCoincDef= lsctables.CoincDef(search = u"ringdown", search_coinc_type = 2, description = u"sim_ringdown<-->coinc_event coincidences (nearby)")

class DocContents(object):
	"""
	A wrapper interface to the XML document.
	"""
	def __init__(self, xmldoc, bbdef, sbdef, scndef, process, sngl_type, sim_type, get_sngl_time):
		#
		# store the process row
		#

		self.process = process

		#
		# locate the sngl_inspiral and sim_inspiral tables
		#

		self.sngltable = sngl_type.get_table(xmldoc)
		try:
			self.simtable = sim_type.get_table(xmldoc)
		except ValueError:
			self.simtable = lsctables.SimInspiralTable.get_table(xmldoc)
			print >>sys.stderr,"No SimRingdownTable, use SimInspiralTable instead!"

		#
		# construct the zero-lag time slide needed to cover the
		# instruments listed in all the triggers, then determine
		# its ID (or create it if needed)
		#
		# FIXME:  in the future, the sim_inspiral table should
		# indicate time slide at which the injection was done
		#

		self.tisi_id = ligolw_tisi.get_time_slide_id(xmldoc, {}.fromkeys(self.sngltable.getColumnByName("ifo"), 0.0), create_new = process)

		#
		# get coinc_definer row for sim_type <--> sngl_type
		# coincs; this creates a coinc_definer table if the
		# document doesn't have one
		#

		self.sb_coinc_def_id = ligolw_coincs.get_coinc_def_id(xmldoc, sbdef.search, sbdef.search_coinc_type, create_new = True, description = sbdef.description)

		#
		# get coinc_def_id's for sngl_type <--> sngl_type, and
		# the sim_type <--> coinc_event coincs.  set all
		# to None if this document does not contain any sngl_type
		# <--> sngl_type coincs.
		#

		try:
			bb_coinc_def_id = ligolw_coincs.get_coinc_def_id(xmldoc, bbdef.search, bbdef.search_coinc_type, create_new = False)
		except KeyError:
			bb_coinc_def_id = None
			self.scn_coinc_def_id = None
		else:
			self.scn_coinc_def_id = ligolw_coincs.get_coinc_def_id(xmldoc, scndef.search, scndef.search_coinc_type, create_new = True, description = scndef.description)

		#
		# get coinc table, create one if needed
		#

		try:
			self.coinctable = lsctables.CoincTable.get_table(xmldoc)
		except ValueError:
			self.coinctable = lsctables.New(lsctables.CoincTable)
			xmldoc.childNodes[0].appendChild(self.coinctable)
		self.coinctable.sync_next_id()

		#
		# get coinc_map table, create one if needed
		#

		try:
			self.coincmaptable = lsctables.CoincMapTable.get_table(xmldoc)
		except ValueError:
			self.coincmaptable = lsctables.New(lsctables.CoincMapTable)
			xmldoc.childNodes[0].appendChild(self.coincmaptable)

		#
		# index the document
		#
		# FIXME:  type<-->type coincs should be organized by time
		# slide ID, but since injections are only done at zero lag
		# for now this is ignored.
		#

		# index sngl_type table
		index = {}
		for row in self.sngltable:
			index[row.event_id] = row
		# find IDs of type<-->type coincs
		self.coincs = {}
		for coinc in self.coinctable:
			if coinc.coinc_def_id == bb_coinc_def_id:
				self.coincs[coinc.coinc_event_id] = []
		# construct event list for each type<-->type coinc
		for row in self.coincmaptable:
			if row.coinc_event_id in self.coincs:
				self.coincs[row.coinc_event_id].append(index[row.event_id])
		del index
		# sort each event list by end/start time and convert to tuples
		# for speed

		for coinc_event_id, events in self.coincs.iteritems():
			events.sort(key=get_sngl_time)
			self.coincs[coinc_event_id]=tuple(events)
		# convert dictionary to a list

		self.coincs = self.coincs.items()

		#
		# FIXME Is this true for inspirals too?
		# sort sngl_type table by end/start time, and sort the coincs
		# list by the end/start time of the first (earliest) event in
		# each coinc (recall that the event tuple for each coinc
		# has been time-ordered)
		#

		self.sngltable.sort(key=get_sngl_time)
		self.coincs.sort(key=lambda(id,a): get_sngl_time(a[0]))

		#
		# set the window for type_near_time().  this window
		# is the amount of time such that if an injection's end
		# time and a inspiral event's end time differ by more than
		# this it is *impossible* for them to match one another.
		#

 		# FIXME I'll just make the windows 1.0 s

                self.search_time_window = 1.0
                self.coinc_time_window = 1.0


	def type_near_time(self, t):
		"""
		Return a list of the inspiral events whose peak times are
		within self.search_time_window of t.
		"""
		return self.sngltable[bisect.bisect_left(self.sngltable, t - self.search_time_window):bisect.bisect_right(self.sngltable, t + self.search_time_window)]

	def coincs_near_time(self, t, get_time):
		"""
		Return a list of the (coinc_event_id, event list) tuples in
		which at least one type event's end time is within
		self.coinc_time_window of t.
		"""
		# FIXME:  this test does not consider the time slide
		# offsets that should be applied to the coinc, but for now
		# injections are done at zero lag so this isn't a problem
		# yet
		return [(coinc_event_id, search_type) for coinc_event_id, search_type in self.coincs if (t - self.coinc_time_window <= get_time(search_type[-1])) and (get_time(search_type[0]) <= t + self.coinc_time_window)]
	
	def sort_triggers_by_id(self):
		"""
		Sort the sngl_table's rows by ID (tidy-up document
		for output).
		"""
		self.sngltable.sort(lambda a, b: cmp(a.event_id, b.event_id))

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


process_program_name = "lalapps_cbc_injfind"


def append_process(xmldoc, match_algorithm, comment):
	"""
	Convenience wrapper for adding process metadata to the document.
	"""
	process = ligolw_process.append_process(xmldoc, program = process_program_name, version = __version__, cvs_repository = u"lscsoft", cvs_entry_time = __date__, comment = comment)

	params = [(u"--match-algorithm", u"lstring", match_algorithm)]
	ligolw_process.append_process_params(xmldoc, process, params)

	return process


#
# =============================================================================
#
#                 Injection <--> Inspiral/Ringdown Event Comparison Tests
#
# =============================================================================
#

def InspiralSnglCompare(sim, inspiral):
        """
	Return False if the peak time of the sim is within 9 seconds of the inspiral event.
        """
	return SnglInspiralUtils.CompareSnglInspiral(sim, inspiral, twindow = LIGOTimeGPS(9))


def InspiralNearCoincCompare(sim, inspiral):
	"""
	Return False if the peak time of the sim is within 9 seconds of the inspiral event.
	"""
	return SnglInspiralUtils.CompareSnglInspiral(sim, inspiral, twindow = LIGOTimeGPS(9))

def cmp_sngl_sim(sim, sngl, get_sim_time, get_sngl_time, twindow = LIGOTimeGPS(9)):
	tdiff = abs(get_sngl_time(sngl) - get_sim_time(sim))
	if tdiff < twindow:
	  return 0
	else:
	  return cmp(get_sngl_time(sngl), get_sim_time(sim))


#
# =============================================================================
#
#                 Build sim_type <--> sngl_type Coincidences
#
# =============================================================================
#


def find_sngl_type_matches(contents, sim, comparefunc, get_sim_time, get_sngl_time) :
	"""
	Scan the type table for triggers matching sim.
	"""
	return [type_table for type_table in contents.type_near_time(get_sim_time(sim)) if not comparefunc(sim, type_table, get_sim_time, get_sngl_time)]


def add_sim_type_coinc(contents, sim, types):
	"""
	Create a coinc_event in the coinc table, and add arcs in the
	coinc_event_map table linking the sim_type row and the list of
	sngl_type rows to the new coinc_event row.
	"""
	coinc = contents.new_coinc(contents.sb_coinc_def_id)
	coinc.set_instruments(set(event.ifo for event in types))
	coinc.nevents = len(types)

	coincmap = lsctables.CoincMap()
	coincmap.coinc_event_id = coinc.coinc_event_id
	coincmap.table_name = sim.simulation_id.table_name
	coincmap.event_id = sim.simulation_id
	contents.coincmaptable.append(coincmap)

	for event in types:
		coincmap = lsctables.CoincMap()
		coincmap.coinc_event_id = coinc.coinc_event_id
		coincmap.table_name = event.event_id.table_name
		coincmap.event_id = event.event_id
		contents.coincmaptable.append(coincmap)


#
# =============================================================================
#
#                   Build sim_type <--> coinc Coincidences
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

def find_near_coinc_matches(coincs, sim, comparefunc, get_sim_time, get_sngl_time):
	"""
	Return a list of the coinc_event_ids of the type<-->type coincs
	in which at least one type event matches sim.
	"""
	# FIXME:  this test does not consider the time slide offsets that
	# should be applied to the coinc, but for now injections are done
	# at zero lag so this isn't a problem yet
	return [coinc_event_id for coinc_event_id, coinc_types in coincs if False in [bool(comparefunc(sim, coinc_type, get_sim_time, get_sngl_time)) for coinc_type in coinc_types]]


def add_sim_coinc_coinc(contents, sim, coinc_event_ids, coinc_def_id):
	"""
	Create a coinc_event in the coinc table, and add arcs in the
	coinc_event_map table linking the sim_type row and the list of
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


def lalapps_cbc_injfind(xmldoc, process, search, snglcomparefunc, nearcoinccomparefunc, verbose = False):
	#
	# Analyze the document's contents.
	#

	if verbose:
		print >>sys.stderr, "indexing ..."

	bbdef = {"inspiral": ligolw_thinca.InspiralCoincDef, "ringdown": ligolw_rinca.RingdownCoincDef}[search]
	sbdef = {"inspiral": InspiralSICoincDef, "ringdown": RingdownSICoincDef}[search]
	scndef = {"inspiral": InspiralSCNearCoincDef, "ringdown": RingdownSCNearCoincDef}[search]
	sngl_type = {"inspiral": lsctables.SnglInspiralTable, "ringdown": lsctables.SnglRingdownTable}[search]
	sim_type = {"inspiral": lsctables.SimInspiralTable, "ringdown": lsctables.SimRingdownTable}[search]
	get_sngl_time = {"inspiral": lsctables.SnglInspiral.get_end, "ringdown": lsctables.SnglRingdown.get_start}[search]

	contents = DocContents(xmldoc = xmldoc, bbdef = bbdef, sbdef = sbdef, scndef = scndef, process = process, sngl_type = sngl_type, sim_type = sim_type, get_sngl_time = get_sngl_time)
	N = len(contents.simtable)

	if isinstance(contents.simtable, lsctables.SimInspiralTable):
		get_sim_time = lsctables.SimInspiral.get_end
	elif isinstance(contents.simtable, lsctables.SimRingdownTable):
		get_sim_time = lsctables.SimRingdown.get_start
	else:
		raise ValueError, "Unknown sim table"

	#
	# Find sim_type <--> sngl_type coincidences.
	#

	if verbose:
		print >>sys.stderr, "constructing %s:" % sbdef.description
	for n, sim in enumerate(contents.simtable):
		if verbose:
			print >>sys.stderr, "\t%.1f%%\r" % (100.0 * n / N),
		types = find_sngl_type_matches(contents, sim, snglcomparefunc, get_sim_time, get_sngl_time)
		if types:
			add_sim_type_coinc(contents, sim, types)
	if verbose:
		print >>sys.stderr, "\t100.0%"

	#
	# Find sim_type <--> coinc_event coincidences.
	#

	if contents.scn_coinc_def_id:
		if verbose:
			print >>sys.stderr, "constructing %s:" % (scndef.description)
		for n, sim in enumerate(contents.simtable):
			if verbose:
				print >>sys.stderr, "\t%.1f%%\r" % (100.0 * n / N),
			coincs = contents.coincs_near_time(get_sim_time(sim), get_sngl_time)
			coinc_event_ids = find_near_coinc_matches(coincs, sim, nearcoinccomparefunc, get_sim_time, get_sngl_time)
			if coinc_event_ids:
				add_sim_coinc_coinc(contents, sim, coinc_event_ids, contents.scn_coinc_def_id)
		if verbose:
			print >>sys.stderr, "\t100.0%"

	#
	# Find sim_type <--> sim_type mappings
	# if both a sim_inspiral and sim_ring table
	# are present.
	#

	if verbose:
		print >>sys.stderr, "Checking for both SimInspiralTable and SimRingdownTable"
	if isinstance(contents.simtable, lsctables.SimRingdownTable):
		try:
			insp_simtable = lsctables.SimInspiralTable.get_table(xmldoc)
		except ValueError:
			print >>sys.stderr,"No SimInspiralTable, only SimRingdownTable present"
			insp_simtable = None
		if insp_simtable is not None:
			if verbose:
				print >> sys.stderr, "found SimInspiralTable, creating maps to SimRingdownTable"
			# create an index of the sim_ringdown geocent_start_times
			sim_ring_time_map = dict([ [(str(row.process_id), row.geocent_start_time), str(row.simulation_id)] for row in contents.simtable])
			# create an index of the sim_ringdown's coinc_event_ids
			sim_ring_ceid_map = dict([ [str(row.event_id), row.coinc_event_id] for row in contents.coincmaptable if row.table_name == "sim_ringdown" ])
			# cycle over the sim_inspiral table, creating the maps
			for this_sim_insp in insp_simtable:
				if (str(this_sim_insp.process_id), this_sim_insp.geocent_end_time) in sim_ring_time_map:
					this_sim_ring_id = sim_ring_time_map[( str(this_sim_insp.process_id), this_sim_insp.geocent_end_time )]
				else:
					continue
				if str(this_sim_ring_id) in sim_ring_ceid_map:
					this_ceid = sim_ring_ceid_map[ str(this_sim_ring_id) ]
				else:
					continue
				# add to the coinc_event_map table
				new_mapping = lsctables.CoincMap()
				new_mapping.table_name = "sim_inspiral"
				new_mapping.event_id = this_sim_insp.simulation_id
				new_mapping.coinc_event_id = this_ceid
				contents.coincmaptable.append(new_mapping)

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
