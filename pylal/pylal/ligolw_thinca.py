# Copyright (C) 2008--2015  Kipp Cannon, Drew G. Keppel
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
from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw.utils import search_summary as ligolw_search_summary
from glue.ligolw.utils import coincs as ligolw_coincs
from glue import offsetvector
import lal
from pylal import git_version
from pylal import snglcoinc
from pylal.xlal import tools as xlaltools
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
from pylal.xlal.datatypes import snglinspiraltable


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
# Construct a subclass of the C sngl_inspiral row class with the methods
# that are needed
#


class SnglInspiral(snglinspiraltable.SnglInspiralTable):
	__slots__ = ()

	spin1 = lsctables.SnglInspiral.spin1
	spin2 = lsctables.SnglInspiral.spin2

	__eq__ = lsctables.SnglInspiral.__eq__

	def __cmp__(self, other):
		# compare self's end time to the LIGOTimeGPS instance
		# other.  allows bisection searches by GPS time to find
		# ranges of triggers quickly
		return cmp(self.end, other)


#
# Use C LIGOTimeGPS type
#


lsctables.LIGOTimeGPS = LIGOTimeGPS


#
# =============================================================================
#
#                          CoincTables Customizations
#
# =============================================================================
#


#
# The sngl_inspiral <--> sngl_inspiral coinc type.
#


InspiralCoincDef = lsctables.CoincDef(search = u"inspiral", search_coinc_type = 0, description = u"sngl_inspiral<-->sngl_inspiral coincidences")


#
# Custom snglcoinc.CoincTables subclass.
#


class InspiralCoincTables(snglcoinc.CoincTables):
	def __init__(self, xmldoc, vetoes = None, program = u"inspiral", likelihood_func = None, likelihood_params_func = None):
		snglcoinc.CoincTables.__init__(self, xmldoc)

		#
		# configure the likelihood ratio evaluator
		#

		if likelihood_func is None and likelihood_params_func is not None or likelihood_func is not None and likelihood_params_func is None:
			raise ValueError("must provide both a likelihood function and a parameter function or neither")
		self.likelihood_func = likelihood_func
		self.likelihood_params_func = likelihood_params_func

		#
		# create a string uniquifier
		#

		self.uniquifier = {}

		#
		# find the coinc_inspiral table or create one if not found
		#

		try:
			self.coinc_inspiral_table = lsctables.CoincInspiralTable.get_table(xmldoc)
		except ValueError:
			self.coinc_inspiral_table = lsctables.New(lsctables.CoincInspiralTable)
			xmldoc.childNodes[0].appendChild(self.coinc_inspiral_table)

		#
		# extract the coalesced out segment lists from the trigger generator
		#

		self.seglists = ligolw_search_summary.segmentlistdict_fromsearchsummary(xmldoc, program = program).coalesce()
		if vetoes is not None:
			self.seglists -= vetoes

	def append_coinc(self, process_id, time_slide_id, coinc_def_id, events):
		#
		# populate the coinc_event and coinc_event_map tables
		#

		coinc = snglcoinc.CoincTables.append_coinc(self, process_id, time_slide_id, coinc_def_id, events)

		#
		# populate the coinc_inspiral table:
		#
		# - end_time is the end time of the first trigger in
		#   alphabetical order by instrument (!?) time-shifted
		#   according to the coinc's offset vector
		# - mass is average of total masses
		# - mchirp is average of mchirps
		# - snr is root-sum-square of SNRs
		# - false-alarm rates are blank
		#

		coinc_inspiral = self.coinc_inspiral_table.RowType()
		coinc_inspiral.coinc_event_id = coinc.coinc_event_id
		coinc_inspiral.mass = sum(event.mass1 + event.mass2 for event in events) / len(events)
		coinc_inspiral.mchirp = sum(event.mchirp for event in events) / len(events)
		coinc_inspiral.snr = math.sqrt(sum(event.snr**2 for event in events))
		coinc_inspiral.false_alarm_rate = None
		coinc_inspiral.combined_far = None
		coinc_inspiral.minimum_duration = None
		offsetvector = self.time_slide_index[time_slide_id]
		coinc_inspiral.end = coinc_inspiral_end_time(events, offsetvector)
		coinc_inspiral.instruments = (event.ifo for event in events)
		self.coinc_inspiral_table.append(coinc_inspiral)

		#
		# record the instruments that were on at the time of the
		# coinc.  note that the start time of the coinc must be
		# unslid to compare with the instrument segment lists
		#

		tstart = coinc_inspiral.end
		instruments = set(event.ifo for event in events)
		instruments |= set(instrument for instrument, segs in self.seglists.items() if tstart - offsetvector[instrument] in segs)
		coinc.set_instruments(instruments)

		#
		# if a likelihood ratio calculator is available, assign a
		# likelihood ratio to the coinc
		#

		if self.likelihood_func is not None:
			coinc.likelihood = self.likelihood_func(self.likelihood_params_func(events, offsetvector))

		#
		# save memory by re-using strings
		#

		coinc.instruments = self.uniquifier.setdefault(coinc.instruments, coinc.instruments)
		coinc_inspiral.ifos = self.uniquifier.setdefault(coinc_inspiral.ifos, coinc_inspiral.ifos)

		#
		# done
		#

		return coinc


#
# Custom function to compute the coinc_inspiral.end_time
#


def coinc_inspiral_end_time(events, offset_vector):
	"""
	Function to compute the end time of an inspiral coincidence.
	@events: a tuple of sngl_inspiral triggers making up a single
	coinc_inspiral trigger
	@offset_vector: a dictionary of offsets to apply to different
	detectors keyed by detector name
	"""
	event = max(events, key = lambda event: event.snr)
	return event.end + offset_vector[event.ifo]


#
# =============================================================================
#
#                            Event List Management
#
# =============================================================================
#


class InspiralEventList(snglcoinc.EventList):
	"""
	A customization of the EventList class for use with the inspiral
	search.
	"""
	def make_index(self):
		"""
		Sort events by end time so that a bisection search can
		retrieve them.  Note that the bisection search relies on
		the __cmp__() method of the SnglInspiral row class having
		previously been set to compare the event's end time to a
		LIGOTimeGPS.
		"""
		self.sort(key = lambda event: event.end)

	def set_dt(self, dt):
		"""
		If an event's end time differs by more than this many
		seconds from the end time of another event then it is
		*impossible* for them to be coincident.
		"""
		# add 1% for safety, and pre-convert to LIGOTimeGPS to
		# avoid doing type conversion in loops
		self.dt = LIGOTimeGPS(dt * 1.01)

	def get_coincs(self, event_a, offset_a, light_travel_time, threshold, comparefunc):
		#
		# event_a's end time, with time shift applied
		#

		end = event_a.end + offset_a - self.offset

		#
		# extract the subset of events from this list that pass
		# coincidence with event_a (use bisection searches for the
		# minimum and maximum allowed end times to quickly identify
		# a subset of the full list)
		#

		return [event_b for event_b in self[bisect.bisect_left(self, end - self.dt) : bisect.bisect_right(self, end + self.dt)] if not comparefunc(event_a, offset_a, event_b, self.offset, light_travel_time, threshold)]


#
# =============================================================================
#
#                              Coincidence Tests
#
# =============================================================================
#


def inspiral_ethinca_max_dt(events, e_thinca_parameter):
	"""
	Given an e-thinca parameter and a list of sngl_inspiral events,
	return the greatest \Delta t that can separate two events and they
	still be considered coincident.
	"""
	# for each instrument present in the event list, compute the
	# largest \Delta t interval for the events from that instrument,
	# and return the sum of the largest two such \Delta t's.
	return sum(sorted(max(xlaltools.XLALSnglInspiralTimeError(event, e_thinca_parameter) for event in events if event.ifo == instrument) for instrument in set(event.ifo for event in events))[-2:]) + 2. * lal.REARTH_SI / lal.C_SI


def inspiral_coinc_compare(a, offseta, b, offsetb, light_travel_time, e_thinca_parameter):
	"""
	Returns False (a & b are coincident) if they pass the ellipsoidal
	thinca test.
	"""
	if offseta: a.end += offseta
	if offsetb: b.end += offsetb
	try:
		# FIXME:  should it be "<" or "<="?
		coincident = xlaltools.XLALCalculateEThincaParameter(a, b) <= e_thinca_parameter
	except ValueError:
		# ethinca test failed to converge == events are not
		# coincident
		coincident = False
	if offseta: a.end -= offseta
	if offsetb: b.end -= offsetb
	return not coincident


def inspiral_coinc_compare_exact(a, offseta, b, offsetb, light_travel_time, delta_t):
	"""
	Returns False (a & b are coincident) if they have the same masses
	and spins (assumed to indicate the templates are identical) and
	their end times are within delta_t of each other.
	"""
	return abs(a.end + offseta - b.end - offsetb) > light_travel_time + delta_t or (a.mass1 != b.mass1) or (a.mass2 != b.mass2) or (a.spin1 != b.spin1).any() or (a.spin2 != b.spin2).any()


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


def replicate_threshold(threshold, instruments):
	"""
	From a single threshold and a list of instruments, return a
	dictionary whose keys are every instrument pair (both orders), and
	whose values are all the same single threshold.

	Example:

	>>> replicate_threshold(6, ["H1", "H2"])
	{("H1", "H2"): 6, ("H2", "H1"): 6}
	"""
	instruments = sorted(instruments)
	thresholds = dict((pair, threshold) for pair in iterutils.choices(instruments, 2))
	instruments.reverse()
	thresholds.update(dict((pair, threshold) for pair in iterutils.choices(instruments, 2)))
	return thresholds


def ligolw_thinca(
	xmldoc,
	process_id,
	coinc_definer_row,
	event_comparefunc,
	thresholds,
	max_dt,
	ntuple_comparefunc = default_ntuple_comparefunc,
	veto_segments = None,
	trigger_program = u"inspiral",
	likelihood_func = None,
	likelihood_params_func = None,
	verbose = False
):
	#
	# prepare the coincidence table interface.
	#

	if verbose:
		print >>sys.stderr, "indexing ..."
	coinc_tables = InspiralCoincTables(xmldoc, vetoes = veto_segments, program = trigger_program, likelihood_func = likelihood_func, likelihood_params_func = likelihood_params_func)
	coinc_def_id = ligolw_coincs.get_coinc_def_id(xmldoc, coinc_definer_row.search, coinc_definer_row.search_coinc_type, create_new = True, description = coinc_definer_row.description)
	sngl_index = dict((row.event_id, row) for row in lsctables.SnglInspiralTable.get_table(xmldoc))

	#
	# build the event list accessors, populated with events from those
	# processes that can participate in a coincidence.  apply vetoes by
	# removing events from the lists that fall in vetoed segments
	#

	eventlists = snglcoinc.EventListDict(InspiralEventList, lsctables.SnglInspiralTable.get_table(xmldoc))
	if veto_segments is not None:
		for eventlist in eventlists.values():
			iterutils.inplace_filter((lambda event: event.ifo not in veto_segments or event.end not in veto_segments[event.ifo]), eventlist)

	#
	# set the \Delta t parameter on all the event lists
	#

	if verbose:
		print >>sys.stderr, "event bisection search window will be %.16g s" % max_dt
	for eventlist in eventlists.values():
		eventlist.set_dt(max_dt)

	#
	# replicate the ethinca parameter for every possible instrument
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

	for node, coinc in time_slide_graph.get_coincs(eventlists, event_comparefunc, thresholds, verbose = verbose):
		coinc = tuple(sngl_index[event_id] for event_id in coinc)
		if not ntuple_comparefunc(coinc, node.offset_vector):
			coinc_tables.append_coinc(process_id, node.time_slide_id, coinc_def_id, coinc)

	#
	# remove time offsets from events
	#

	del eventlists.offsetvector

	#
	# done
	#

	return xmldoc


#
# =============================================================================
#
#                              GraceDB Utilities
#
# =============================================================================
#


#
# Device to extract sngl_inspiral coincs from a source XML document tree.
#


class sngl_inspiral_coincs(object):
	"""
	Dictionary-like device to extract XML document trees containing
	individual sngl_inspiral coincs from a source XML document tree
	containing several.

	An instance of the class is initialized with an XML document tree.
	The coinc event ID of a sngl_inspiral<-->sngl_inspiral coinc in
	the document can then be used like a dictionary key to retrieve a
	newly-constructed XML document containing that coinc by itself.
	The output document trees are complete, self-describing, documents
	with all metadata about the event from the source document
	preserved.

	Example:

	>>> coincs = sngl_inspiral_coincs(xmldoc)
	>>> print coincs.coinc_def_id
	coinc_definer:coinc_def_id:0
	>>> coincs.keys()
	[<glue.ligolw.ilwd.cached_ilwdchar_class object at 0x41a4328>]
	>>> coinc_id = coincs.keys()[0]
	>>> print coinc_id
	coinc_event:coinc_event_id:83763
	>>> coincs[coinc_id].write()
	<?xml version='1.0' encoding='utf-8'?>
	<!DOCTYPE LIGO_LW SYSTEM "http://ldas-sw.ligo.caltech.edu/doc/ligolwAPI/html/ligolw_dtd.txt">
	<LIGO_LW>
		<Table Name="process:table">
			<Column Type="lstring" Name="process:comment"/>
			<Column Type="lstring" Name="process:node"/>
	...

	The XML documents returned from this class share references to the
	row objects in the original document.  Modifications to the row
	objects in the tables returned by this class will affect both the
	original document and all other documents returned by this class.
	However, each retrieval constructs a new document from scratch,
	they are not cached nor re-used, therefore this operation can be
	time consuming if it needs to be performed repeatedly but the table
	objects and document trees can be edited without affecting each
	other.

	If the source document is modified after this class has been
	instantiated, the behaviour is undefined.

	To assist with memory clean-up, it is helpful to invoke the
	.unlink() method on the XML trees returned by this class when they
	are no longer needed.
	"""
	def __init__(self, xmldoc):
		"""
		Initialize an instance of the class.  xmldoc is the source
		XML document tree from which the
		sngl_inspiral<-->sngl_inspiral coincs will be extracted.
		"""
		#
		# find all tables
		#

		self.process_table = lsctables.ProcessTable.get_table(xmldoc)
		self.process_params_table = lsctables.ProcessParamsTable.get_table(xmldoc)
		self.search_summary_table = lsctables.SearchSummaryTable.get_table(xmldoc)
		self.sngl_inspiral_table = lsctables.SnglInspiralTable.get_table(xmldoc)
		self.coinc_def_table = lsctables.CoincDefTable.get_table(xmldoc)
		self.coinc_event_table = lsctables.CoincTable.get_table(xmldoc)
		self.coinc_inspiral_table = lsctables.CoincInspiralTable.get_table(xmldoc)
		self.coinc_event_map_table = lsctables.CoincMapTable.get_table(xmldoc)
		self.time_slide_table = lsctables.TimeSlideTable.get_table(xmldoc)

		#
		# index the process, process params, search_summary,
		# sngl_inspiral and time_slide tables
		#

		self.process_index = dict((row.process_id, row) for row in self.process_table)
		self.process_params_index = {}
		for row in self.process_params_table:
			self.process_params_index.setdefault(row.process_id, []).append(row)
		self.search_summary_index = dict((row.process_id, row) for row in self.search_summary_table)
		self.sngl_inspiral_index = dict((row.event_id, row) for row in self.sngl_inspiral_table)
		self.time_slide_index = {}
		for row in self.time_slide_table:
			self.time_slide_index.setdefault(row.time_slide_id, []).append(row)

		#
		# find the sngl_inspiral<-->sngl_inspiral coincs
		#

		self.coinc_def, = (row for row in self.coinc_def_table if row.search == InspiralCoincDef.search and row.search_coinc_type == InspiralCoincDef.search_coinc_type)
		self.coinc_event_index = dict((row.coinc_event_id, row) for row in self.coinc_event_table if row.coinc_def_id == self.coinc_def.coinc_def_id)
		self.coinc_inspiral_index = dict((row.coinc_event_id, row) for row in self.coinc_inspiral_table)
		assert set(self.coinc_event_index) == set(self.coinc_inspiral_index)
		self.coinc_event_map_index = dict((coinc_event_id, []) for coinc_event_id in self.coinc_event_index)
		for row in self.coinc_event_map_table:
			try:
				self.coinc_event_map_index[row.coinc_event_id].append(row)
			except KeyError:
				continue

	@property
	def coinc_def_id(self):
		"""
		The coinc_def_id of the sngl_inspiral<-->sngl_inspiral
		coincs in the source XML document.
		"""
		return self.coinc_def.coinc_def_id

	def sngl_inspirals(self, coinc_event_id):
		"""
		Return a list of the sngl_inspiral rows that participated
		in the coincidence given by coinc_event_id.
		"""
		return [self.sngl_inspiral_index[row.event_id] for row in self.coinc_event_map_index[coinc_event_id]]

	def offset_vector(self, time_slide_id):
		"""
		Return the offsetvector given by time_slide_id.
		"""
		return offsetvector.offsetvector((row.instrument, row.offset) for row in self.time_slide_index[time_slide_id])

	def __getitem__(self, coinc_event_id):
		"""
		Construct and return an XML document containing the
		sngl_inspiral<-->sngl_inspiral coinc carrying the given
		coinc_event_id.
		"""
		newxmldoc = ligolw.Document()
		newxmldoc.appendChild(ligolw.LIGO_LW())

		# when making these, we can't use table.new_from_template()
		# because we need to ensure we have a Table subclass, not a
		# DBTable subclass
		new_process_table = newxmldoc.childNodes[-1].appendChild(lsctables.New(lsctables.ProcessTable, self.process_table.columnnames))
		new_process_params_table = newxmldoc.childNodes[-1].appendChild(lsctables.New(lsctables.ProcessParamsTable, self.process_params_table.columnnames))
		new_search_summary_table = newxmldoc.childNodes[-1].appendChild(lsctables.New(lsctables.SearchSummaryTable, self.search_summary_table.columnnames))
		new_sngl_inspiral_table = newxmldoc.childNodes[-1].appendChild(lsctables.New(lsctables.SnglInspiralTable, self.sngl_inspiral_table.columnnames))
		new_coinc_def_table = newxmldoc.childNodes[-1].appendChild(lsctables.New(lsctables.CoincDefTable, self.coinc_def_table.columnnames))
		new_coinc_event_table = newxmldoc.childNodes[-1].appendChild(lsctables.New(lsctables.CoincTable, self.coinc_event_table.columnnames))
		new_coinc_inspiral_table = newxmldoc.childNodes[-1].appendChild(lsctables.New(lsctables.CoincInspiralTable, self.coinc_inspiral_table.columnnames))
		new_coinc_event_map_table = newxmldoc.childNodes[-1].appendChild(lsctables.New(lsctables.CoincMapTable, self.coinc_event_map_table.columnnames))
		new_time_slide_table = newxmldoc.childNodes[-1].appendChild(lsctables.New(lsctables.TimeSlideTable, self.time_slide_table.columnnames))

		new_coinc_def_table.append(self.coinc_def)
		coinc_event = self.coinc_event_index[coinc_event_id]
		new_coinc_event_table.append(coinc_event)
		new_coinc_inspiral_table.append(self.coinc_inspiral_index[coinc_event_id])
		map(new_coinc_event_map_table.append, self.coinc_event_map_index[coinc_event_id])
		map(new_time_slide_table.append, self.time_slide_index[coinc_event.time_slide_id])
		for row in new_coinc_event_map_table:
			new_sngl_inspiral_table.append(self.sngl_inspiral_index[row.event_id])

		for process_id in set(new_sngl_inspiral_table.getColumnByName("process_id")) | set(new_coinc_event_table.getColumnByName("process_id")) | set(new_time_slide_table.getColumnByName("process_id")):
			# process row is required
			new_process_table.append(self.process_index[process_id])
			try:
				map(new_process_params_table.append, self.process_params_index[process_id])
			except KeyError:
				# process_params rows are optional
				pass
			try:
				new_search_summary_table.append(self.search_summary_index[process_id])
			except KeyError:
				# search_summary rows are optional
				pass

		return newxmldoc

	def __iter__(self):
		"""
		Iterate over the coinc_event_id's in the source document.
		"""
		return iter(self.coinc_event_index)

	def __nonzero__(self):
		return bool(self.coinc_event_index)

	def keys(self):
		"""
		A list of the coinc_event_id's of the
		sngl_inspiral<-->sngl_inspiral coincs available in the
		source XML document.
		"""
		return self.coinc_event_index.keys()

	def items(self):
		"""
		Yield a sequence of (coinc_event_id, XML tree) tuples, one
		for each sngl_inspiral<-->sngl_inspiral coinc in the source
		document.

		NOTE:  to allow this to work more easily with very large
		documents, instead of returning the complete sequence as a
		pre-constructed list this method is implemented as a
		generator.
		"""
		for coinc_event_id in self:
			yield (coinc_event_id, self[coinc_event_id])

	iteritems = items
