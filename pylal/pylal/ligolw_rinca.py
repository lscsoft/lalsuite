# Copyright (C) 2008  Kipp Cannon, Drew G. Keppel
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


import bisect
import math
import sys


from glue import iterutils
from glue.ligolw import lsctables
from glue.ligolw.utils import process as ligolw_process
from glue.ligolw.utils import search_summary as ligolw_search_summary
from pylal import git_version
from pylal import llwapp
from pylal import snglcoinc
from pylal.xlal import tools as xlaltools
from pylal.xlal.datatypes import snglringdowntable
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS


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


#
# Use C row classes for memory and speed
#


lsctables.CoincMapTable.RowType = lsctables.CoincMap = xlaltools.CoincMap


#
# Construct a subclass of the C sngl_ringdown row class with the methods
# that are needed
#


class SnglRingdown(snglringdowntable.SnglRingdownTable):
	__slots__ = ()

	def get_start(self):
		return LIGOTimeGPS(self.start_time, self.start_time_ns)

	def set_start(self, gps):
		self.start_time, self.start_time_ns = gps.seconds, gps.nanoseconds

	def __cmp__(self, other):
		# compare self's start time to the LIGOTimeGPS instance
		# other.  allows bisection searches by GPS time to find
		# ranges of triggers quickly
		return cmp(self.start_time, other.seconds) or cmp(self.start_time_ns, other.nanoseconds)


#
# Use C LIGOTimeGPS type
#


lsctables.LIGOTimeGPS = LIGOTimeGPS


#
# Utilities
#

def coinc_ringdown_start(events, offset_vector):
	
	events = sorted(events, lambda a, b: cmp(a.ifo, b.ifo))
	tstart = events[0].get_start() + offset_vector[events[0].ifo]
	return tstart + sum(event.snr * float(event.get_start() + offset_vector[event.ifo] - tstart) for event in events) / sum(event.snr for event in events)

#
# =============================================================================
#
#                           Add Process Information
#
# =============================================================================
#


process_program_name = "ligolw_rinca"


def append_process(xmldoc, comment = None, force = None, ds_sq_threshold = None, save_small_coincs = None, vetoes_name = None, coinc_end_time_segment = None, verbose = None):
	process = ligolw_process.append_process(xmldoc, program = process_program_name, version = __version__, cvs_repository = u"lscsoft", cvs_entry_time = __date__, comment = comment)

	params = [
		(u"--ds-sq-threshold", u"real_8", ds_sq_threshold)
	]
	if comment is not None:
		params += [(u"--comment", u"lstring", comment)]
	if force is not None:
		params += [(u"--force", None, None)]
	if save_small_coincs is not None:
		params += [(u"--save-small-coincs", None, None)]
	if vetoes_name is not None:
		params += [(u"--vetoes-name", u"lstring", vetoes_name)]
	if coinc_end_time_segment is not None:
		params += [(u"--coinc-end-time-segment", u"lstring", coinc_end_time_segment)]
	if verbose is not None:
		params += [(u"--verbose", None, None)]

	ligolw_process.append_process_params(xmldoc, process, params)

	return process


#
# =============================================================================
#
#                          CoincTables Customizations
#
# =============================================================================
#


#
# The sngl_ringdown <--> sngl_ringdown coinc type.
#


RingdownCoincDef = lsctables.CoincDef(search = u"ring", search_coinc_type = 0, description = u"sngl_ringdown<-->sngl_ringdown coincidences")


#
# Custom snglcoinc.CoincTables subclass.
#


class RingdownCoincTables(snglcoinc.CoincTables):
	def __init__(self, xmldoc, vetoes = None, program = u"lalapps_ring"):
		snglcoinc.CoincTables.__init__(self, xmldoc)

		#
		# create a string uniquifier
		#

		self.uniquifier = {}

		#
		# find the coinc_ringdown table or create one if not found
		#

		try:
			self.coinc_ringdown_table = lsctables.table.get_table(xmldoc, lsctables.CoincRingdownTable.tableName)
		except ValueError:
			self.coinc_ringdown_table = lsctables.New(lsctables.CoincRingdownTable)
			xmldoc.childNodes[0].appendChild(self.coinc_ringdown_table)

		#
		# extract the coalesced out segment lists from lalapps_ring
		#

		self.seglists = ligolw_search_summary.segmentlistdict_fromsearchsummary(xmldoc, program = program).coalesce()
		if vetoes is not None:
			self.seglists -= vetoes

	def append_coinc(self, process_id, node, coinc_def_id, events):
		#
		# populate the coinc_event and coinc_event_map tables
		#
		
		time_slide_id = node.time_slide_id
		
		coinc = snglcoinc.CoincTables.append_coinc(self, process_id, time_slide_id, coinc_def_id, events)

		#
		# populate the coinc_ringdown table:
		#
		# - start_time is the SNR-weighted mean of the start times
		# - frequency and quality are the SNR-weighted means of the
		#   frequencies and qualities
		# - snr is root-sum-square of SNRs
		# - false-alarm rate is blank
		#

		coinc_ringdown = self.coinc_ringdown_table.RowType()
		coinc_ringdown.coinc_event_id = coinc.coinc_event_id
		coinc_ringdown.snr = sum(event.snr**2. for event in events)**.5
		coinc_ringdown.false_alarm_rate = None
		# use the time of event[0] as an epoch
		tstart = coinc_ringdown_start(events, node.offset_vector)
		coinc_ringdown.set_start(tstart)
		#tstart = events[0].get_start() + self.time_slide_index[time_slide_id][events[0].ifo]
		#coinc_ringdown.set_start(tstart + sum(event.snr * float(event.get_start() + self.time_slide_index[time_slide_id][event.ifo] - tstart) for event in events) / sum(event.snr for event in events))
		coinc_ringdown.set_ifos(event.ifo for event in events)
		coinc_ringdown.frequency = sum(event.snr * event.frequency for event in events) / sum(event.snr for event in events)
		coinc_ringdown.quality = sum(event.snr * event.quality for event in events) / sum(event.snr for event in events)
		coinc_ringdown.mass = None
		coinc_ringdown.spin = None
		coinc_ringdown.snr_sq = None
		coinc_ringdown.choppedl_snr = None
		coinc_ringdown.eff_coh_snr = None
		coinc_ringdown.null_stat = None
		coinc_ringdown.kappa = None
		coinc_ringdown.snr_ratio = None
		coinc_ringdown.combined_far = None
		self.coinc_ringdown_table.append(coinc_ringdown)

		#
		# record the instruments that were on at the time of the
		# coinc.  note that the start time of the coinc must be
		# unslid to compare with the instrument segment lists
		#

		#tstart = coinc_ringdown.get_start()
		instruments = set([event.ifo for event in events])
		instruments |= set([instrument for instrument, segs in self.seglists.items() if tstart - self.time_slide_index[time_slide_id][instrument] in segs])
		coinc.set_instruments(instruments)

		#
		# save memory by re-using strings
		#

		coinc.instruments = self.uniquifier.setdefault(coinc.instruments, coinc.instruments)
		coinc_ringdown.ifos = self.uniquifier.setdefault(coinc_ringdown.ifos, coinc_ringdown.ifos)

		#
		# done
		#

		return coinc


#
# =============================================================================
#
#                            Event List Management
#
# =============================================================================
#


class RingdownEventList(snglcoinc.EventList):
	"""
	A customization of the EventList class for use with the ringdown
	search.
	"""
	def make_index(self):
		"""
		Sort events by start time so that a bisection search can
		retrieve them.  Note that the bisection search relies on
		the __cmp__() method of the SnglRingdown row class having
		previously been set to compare the event's start time to a
		LIGOTimeGPS.
		"""
		self.sort(lambda a, b: cmp(a.start_time, b.start_time) or cmp(a.start_time_ns, b.start_time_ns))

	def set_dt(self, dt):
		"""
		If an event's start time differs by more than this many
		seconds from the start time of another event then it is
		*impossible* for them to be coincident.
		"""
		# add 1% for safety, and pre-convert to LIGOTimeGPS to
		# avoid doing type conversion in loops
		self.dt = LIGOTimeGPS(dt * 1.01)

	def get_coincs(self, event_a, offset_a, light_travel_time, ds_sq_threshold, comparefunc):
		#
		# event_a's start time with time shift applied
		#

		start = event_a.get_start() + offset_a - self.offset

		#
		# extract the subset of events from this list that pass
		# coincidence with event_a (use bisection searches for the
		# minimum and maximum allowed start times to quickly identify
		# a subset of the full list)
		#

		return [event_b for event_b in self[bisect.bisect_left(self, start - self.dt) : bisect.bisect_right(self, start + self.dt)] if not comparefunc(event_a, offset_a, event_b, self.offset, light_travel_time, ds_sq_threshold)]


#
# =============================================================================
#
#                              Coincidence Tests
#
# =============================================================================
#


def ringdown_max_dt(events, ds_sq_threshold):
	"""
	Given a ds_sq threshold and a list of sngl_ringdown events,
	return the greatest \Delta t that can separate two events and they
	still be considered coincident.
	"""
	# for each instrument present in the event list, compute the
	# largest \Delta t interval for the events from that instrument,
	# and return the sum of the largest two such \Delta t's.

	# FIXME: get these from somewhere else
	LAL_REARTH_SI = 6.378140e6 # m
	LAL_C_SI = 299792458 # m s^-1

	return sum(sorted(max(xlaltools.XLALRingdownTimeError(event, ds_sq_threshold) for event in events if event.ifo == instrument) for instrument in set(event.ifo for event in events))[-2:]) + 2. * LAL_REARTH_SI / LAL_C_SI


def ringdown_coinc_compare(a, offseta, b, offsetb, light_travel_time, ds_sq_threshold):
	"""
	Returns False (a & b are coincident) if they pass the metric
	rinca test.
	"""
	if offseta: a.set_start(a.get_start() + offseta)
	if offsetb: b.set_start(b.get_start() + offsetb)
	try:
		# FIXME:  should it be "<" or "<="?
		coincident = xlaltools.XLAL3DRinca(a, b) <= ds_sq_threshold
	except ValueError:
		# ds_sq test failed to converge == events are not
		# coincident
		coincident = False
	if offseta: a.set_start(a.get_start() - offseta)
	if offsetb: b.set_start(b.get_start() - offsetb)
	return not coincident


#
# =============================================================================
#
#                              Compare Functions
#
# =============================================================================
#


def default_ntuple_comparefunc(events, offset_vector):
	"""
	Default ntuple test function.  Accept all ntuples.
	"""
	return False


#
# =============================================================================
#
#                                 Library API
#
# =============================================================================
#


def replicate_threshold(ds_sq_threshold, instruments):
	"""
	From a single threshold and a list of instruments, return a
	dictionary whose keys are every instrument pair (both orders), and
	whose values are all the same single threshold.

	Example:

	>>> replicate_threshold(6, ["H1", "H2"])
	{("H1", "H2"): 6, ("H2", "H1"): 6}
	"""
	instruments = sorted(instruments)
	thresholds = dict((pair, ds_sq_threshold) for pair in iterutils.choices(instruments, 2))
	instruments.reverse()
	thresholds.update(dict((pair, ds_sq_threshold) for pair in iterutils.choices(instruments, 2)))
	return thresholds


def ligolw_rinca(
	xmldoc,
	process_id,
	EventListType,
	CoincTables,
	coinc_definer_row,
	event_comparefunc,
	thresholds,
	ntuple_comparefunc = lambda events, offset_vector: False,
	small_coincs = False,
	veto_segments = None,
	coinc_end_time_segment = None,
	verbose = False
):
	#
	# prepare the coincidence table interface.
	#

	if verbose:
		print >>sys.stderr, "indexing ..."
	coinc_tables = CoincTables(xmldoc, vetoes = veto_segments)
	coinc_def_id = llwapp.get_coinc_def_id(xmldoc, coinc_definer_row.search, coinc_definer_row.search_coinc_type, create_new = True, description = coinc_definer_row.description)
	sngl_index = dict((row.event_id, row) for row in lsctables.table.get_table(xmldoc, lsctables.SnglRingdownTable.tableName))

	#
	# build the event list accessors, populated with events from those
	# processes that can participate in a coincidence.  apply vetoes by
	# removing events from the lists that fall in vetoed segments
	#

	eventlists = snglcoinc.EventListDict(EventListType, lsctables.SnglRingdownTable.get_table(xmldoc))
	if veto_segments is not None:
		for eventlist in eventlists.values():
			iterutils.inplace_filter((lambda event: event.ifo not in veto_segments or event.get_start() not in veto_segments[event.ifo]), eventlist)

	#
	# set the \Delta t parameter on all the event lists
	#

	max_dt = ringdown_max_dt(lsctables.table.get_table(xmldoc, lsctables.SnglRingdownTable.tableName), thresholds)
	if verbose:
		print >>sys.stderr, "event bisection search window will be %.16g s" % max_dt
	for eventlist in eventlists.values():
		eventlist.set_dt(max_dt)

	#
	# replicate the ds_sq threshold for every possible instrument
	# pair
	#

	thresholds = replicate_threshold(thresholds, set(eventlists))

	#
	# construct offset vector assembly graph
	#

	time_slide_graph = snglcoinc.TimeSlideGraph(coinc_tables.time_slide_index, verbose = verbose)

	#
	# retrieve all coincidences, apply the final n-tuple compare func
	# and record the survivors
	#

	for node, coinc in time_slide_graph.get_coincs(eventlists, event_comparefunc, thresholds, include_small_coincs = small_coincs, verbose = verbose):
		ntuple = tuple(sngl_index[id] for id in coinc)
		if not ntuple_comparefunc(ntuple, node.offset_vector):
			coinc_tables.append_coinc(process_id, node, coinc_def_id, ntuple)

	#
	# remove time offsets from events
	#

	del eventlists.offsetvector

	#
	# done
	#

	return xmldoc
