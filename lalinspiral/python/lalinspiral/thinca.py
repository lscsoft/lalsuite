# Copyright (C) 2008--2017  Kipp Cannon, Drew G. Keppel
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


from __future__ import print_function
from bisect import bisect_left
import itertools
import math
import sys
import time


from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue import offsetvector
import lal
from lalburst import snglcoinc


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
from .git_version import date as __date__
from .git_version import version as __version__


#
# =============================================================================
#
#                                 Speed Hacks
#
# =============================================================================
#


#
# Construct a subclass of the sngl_inspiral row class with the methods that
# are needed
#


class SnglInspiral(lsctables.SnglInspiral):
	__slots__ = ()
	def __cmp__(self, other):
		# compare self's end time to the LIGOTimeGPS instance
		# other.  allows bisection searches by GPS time to find
		# ranges of triggers quickly
		return cmp(self.end, other)

	#
	# simulate mtotal, eta, and mchirp from mass1 and mass2.  this (a)
	# works around documents that have incorrect values in those three
	# columns (yes, yes people do that) and (b) allows us to process
	# documents that don't have the columns at all
	#

	@property
	def mtotal(self):
		return self.mass1 + self.mass2

	@property
	def eta(self):
		return self.mass1 * self.mass2 / self.mtotal**2.

	@property
	def mchirp(self):
		return self.mtotal * self.eta**0.6


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
# Definition of coinc_inspiral.end_time
#


def coinc_inspiral_end_time(events, offset_vector):
	"""
	Compute the end time of an inspiral coincidence.  events is an
	iterable of sngl_inspiral triggers, offset_vector is a dictionary
	mapping instrument to offset.

	This function returns the time shifted end time of the trigger with
	the highest SNR.  The end time reported by this function gets used
	for things like plot titles, alert messages, and so on.  It is not
	meant to be an accurate estimate of the time at which the
	gravitational wave passed through the geocentre, or any other such
	thing.

	This end time is also used to parallelize thinca by allowing a
	single lock stretch to be split across several jobs without missing
	or double counting any coincs.  This is achieved by using a
	definition that is guaranteed to return a bit-identical "end time"
	for a given set of triggers.  Guaranteeing that allows thinca to
	clip coincs to a sequence of contiguous segments and know that no
	coinc will missed or double counted.
	"""
	event = max(events, key = lambda event: event.snr)
	return event.end + offset_vector[event.ifo]


#
# Custom snglcoinc.CoincTables subclass.
#


class InspiralCoincTables(snglcoinc.CoincTables):
	def __init__(self, xmldoc, coinc_definer_row):
		super(InspiralCoincTables, self).__init__(xmldoc, coinc_definer_row)

		#
		# find the coinc_inspiral table or create one if not found
		#

		try:
			self.coinc_inspiral_table = lsctables.CoincInspiralTable.get_table(xmldoc)
		except ValueError:
			self.coinc_inspiral_table = lsctables.New(lsctables.CoincInspiralTable)
			xmldoc.childNodes[0].appendChild(self.coinc_inspiral_table)


	def coinc_rows(self, process_id, time_slide_id, events, seglists = None):
		coinc, coincmaps = super(InspiralCoincTables, self).coinc_rows(process_id, time_slide_id, events)

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

		offsetvector = self.time_slide_index[time_slide_id]
		end = coinc_inspiral_end_time(events, offsetvector)
		coinc_inspiral = self.coinc_inspiral_table.RowType(
			coinc_event_id = coinc.coinc_event_id,	# = None
			mass = sum(event.mass1 + event.mass2 for event in events) / len(events),
			mchirp = sum(event.mchirp for event in events) / len(events),
			snr = math.sqrt(sum(event.snr**2. for event in events)),
			false_alarm_rate = None,
			combined_far = None,
			minimum_duration = None,
			end = end,
			instruments = (event.ifo for event in events)
		)

		#
		# record the instruments that were on at the time of the
		# coinc.  instruments that provide triggers are, by
		# definition, on.  note that the end time of the coinc
		# must be unslid to compare with the instrument segment
		# lists
		#

		coinc.insts = set(event.ifo for event in events)
		if seglists is not None:
			coinc.insts |= set(instrument for instrument, segs in seglists.items() if end - offsetvector[instrument] in segs)

		#
		# done
		#

		return coinc, coincmaps, coinc_inspiral


	@staticmethod
	def ntuple_comparefunc(events, offset_vector):
		"""
		Default ntuple test function.  Accept all ntuples.
		"""
		return False


	def append_coinc(self, coinc_event, coinc_event_maps, coinc_inspiral):
		coinc_event = super(InspiralCoincTables, self).append_coinc(coinc_event, coinc_event_maps)
		coinc_inspiral.coinc_event_id = coinc_event.coinc_event_id
		self.coinc_inspiral_table.append(coinc_inspiral)
		return coinc_event


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
	@staticmethod
	def template(event):
		"""
		Returns an immutable hashable object (it can be used as a
		dictionary key) uniquely identifying the template that
		produced the given event.
		"""
		return event.mass1, event.mass2, event.spin1x, event.spin1y, event.spin1z, event.spin2x, event.spin2y, event.spin2z

	def make_index(self):
		"""
		Sort events into bins according to their template so that a
		dictionary look-up can retrieve all triggers from a given
		template.  Then sort bins by end time so that a bisection
		search can retrieve triggers from a template within a
		window of time.  Note that the bisection search relies on
		the __cmp__() method of the SnglInspiral row class having
		previously been set to compare the event's end time to a
		LIGOTimeGPS.
		"""
		self.index = {}
		for event in self:
			self.index.setdefault(self.template(event), []).append(event)
		for events in self.index.values():
			events.sort(key = lambda event: event.end)

	def get_coincs(self, event_a, offset_a, light_travel_time, delta_t):
		#
		# event_a's end time, shifted to be with respect to end
		# times in this list.
		#

		end = event_a.end + offset_a

		#
		# the coincidence window
		#

		coincidence_window = light_travel_time + delta_t

		#
		# extract the subset of events from this list that pass
		# coincidence with event_a.  use a bisection search for the
		# minimum allowed end time and a brute-force scan for the
		# maximum allowed end time.  because the number of events
		# in the coincidence window is generally quite small, the
		# brute-force scan has a lower expected operation count
		# than a second bisection search to find the upper bound in
		# the sequence
		#

		try:
			events = self.index[self.template(event_a)]
		except KeyError:
			# that template didn't produce any events in this
			# instrument.  this is rare given the SNR thresholds
			# typical today so trapping the exception is more
			# efficient than testing
			return ()
		stop = end + coincidence_window
		return tuple(itertools.takewhile(lambda event: event.end <= stop, events[bisect_left(events, end - coincidence_window):]))


#
# =============================================================================
#
#                                 Library API
#
# =============================================================================
#


def ligolw_thinca(
	xmldoc,
	process_id,
	delta_t,
	ntuple_comparefunc = InspiralCoincTables.ntuple_comparefunc,
	seglists = None,
	veto_segments = None,
	likelihood_func = None,
	fapfar = None,
	min_instruments = 2,
	min_log_L = None,
	coinc_definer_row = InspiralCoincDef,
	verbose = False
):
	#
	# prepare the coincidence table interface.
	#

	if verbose:
		print("indexing ...", file=sys.stderr)
	coinc_tables = InspiralCoincTables(xmldoc, coinc_definer_row)

	#
	# build the event list accessors.  apply vetoes by excluding events
	# from the lists that fall in vetoed segments
	#

	sngl_inspiral_table = lsctables.SnglInspiralTable.get_table(xmldoc)
	if veto_segments is not None:
		sngl_inspiral_table = (event for event in sngl_inspiral_table if event.ifo not in veto_segments or event.end not in veto_segments[event.ifo])
		if seglists is not None:
			# don't do in-place
			seglists = seglists - veto_segments

	#
	# construct offset vector assembly graph
	#

	time_slide_graph = snglcoinc.TimeSlideGraph(coinc_tables.time_slide_index, min_instruments = min_instruments, verbose = verbose)

	#
	# retrieve all coincidences, apply the final n-tuple compare func
	# and record the survivors
	#

	gps_time_now = float(lal.UTCToGPS(time.gmtime()))
	for node, events in time_slide_graph.get_coincs(
		snglcoinc.EventListDict(InspiralEventList, sngl_inspiral_table, instruments = set(coinc_tables.time_slide_table.getColumnByName("instrument"))),
		delta_t,
		verbose = verbose
	):
		if not ntuple_comparefunc(events, node.offset_vector):
			coinc, coincmaps, coinc_inspiral = coinc_tables.coinc_rows(process_id, node.time_slide_id, events, seglists = seglists)
			if likelihood_func is not None:
				coinc.likelihood = likelihood_func(events, node.offset_vector)
				if fapfar is not None:
					# FIXME:  add proper columns to
					# store these values in
					coinc_inspiral.combined_far = fapfar.far_from_rank(coinc.likelihood)
					coinc_inspiral.false_alarm_rate = fapfar.fap_from_rank(coinc.likelihood)
			# if min_log_L is None, this test always passes,
			# regardless of the value of .likelihood, be it
			# None, some number, -inf or even nan.
			if coinc.likelihood >= min_log_L:
				# set latency.
				# NOTE:  this is nonsense unless running
				# live.
				# FIXME: add a proper column for this
				coinc_inspiral.minimum_duration = gps_time_now - float(coinc_inspiral.end)
				# finally, append coinc to tables
				coinc_tables.append_coinc(coinc, coincmaps, coinc_inspiral)

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
	>>> print(coincs.coinc_def_id)
	coinc_definer:coinc_def_id:0
	>>> coincs.keys()
	[<glue.ligolw.ilwd.cached_ilwdchar_class object at 0x41a4328>]
	>>> coinc_id = coincs.keys()[0]
	>>> print(coinc_id)
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
		self.sngl_inspiral_table = lsctables.SnglInspiralTable.get_table(xmldoc)
		self.coinc_def_table = lsctables.CoincDefTable.get_table(xmldoc)
		self.coinc_event_table = lsctables.CoincTable.get_table(xmldoc)
		self.coinc_inspiral_table = lsctables.CoincInspiralTable.get_table(xmldoc)
		self.coinc_event_map_table = lsctables.CoincMapTable.get_table(xmldoc)
		self.time_slide_table = lsctables.TimeSlideTable.get_table(xmldoc)

		#
		# index the process, process params, sngl_inspiral and
		# time_slide tables
		#

		self.process_index = dict((row.process_id, row) for row in self.process_table)
		self.process_params_index = {}
		for row in self.process_params_table:
			self.process_params_index.setdefault(row.process_id, []).append(row)
		self.sngl_inspiral_index = dict((row.event_id, row) for row in self.sngl_inspiral_table)
		self.time_slide_index = {}
		for row in self.time_slide_table:
			self.time_slide_index.setdefault(row.time_slide_id, []).append(row)
		self.zero_lag_time_slide_ids = frozenset(time_slide_id for time_slide_id, offset_vector in self.time_slide_table.as_dict().items() if not any(offset_vector.values()))

		#
		# find the sngl_inspiral<-->sngl_inspiral coincs
		#

		self.coinc_def, = (row for row in self.coinc_def_table if row.search == InspiralCoincDef.search and row.search_coinc_type == InspiralCoincDef.search_coinc_type)
		coinc_event_map_ids = frozenset(row.coinc_event_id for row in self.coinc_event_map_table)
		self.coinc_event_index = dict((row.coinc_event_id, row) for row in self.coinc_event_table if row.coinc_def_id == self.coinc_def.coinc_def_id and row.coinc_event_id in coinc_event_map_ids)
		self.coinc_inspiral_index = dict((row.coinc_event_id, row) for row in self.coinc_inspiral_table if row.coinc_event_id in coinc_event_map_ids)
		assert frozenset(self.coinc_event_index) == frozenset(self.coinc_inspiral_index)
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

	def single_sngl_inspirals(self):
		"""
		Generator returns a sequence of the sngl_inspiral table
		rows that formed zero-lag single-instrument "coincs".

		This is only meaningful if the coincidence engine was run
		with min_instruments = 1, otherwise this sequence will be
		empty by construction.  Also, if there was no zero-lag time
		slide included in the time slide graph then this sequence
		will be empty.

		This method is used by codes that want lists of
		non-coincident triggers for background models even if
		min_instruments has been set below 2.

		The constraint that they be "zero-lag" singles might at
		first seem nonsensical but is included to exclude triggers
		that form genuine coincidences at zero-lag but are present
		only as single-detector candidates in one or more time
		slides.
		"""
		for coinc_event_id, coinc_event in self.coinc_event_index.items():
			if coinc_event.time_slide_id in self.zero_lag_time_slide_ids and coinc_event.nevents < 2:
				row, = self.coinc_event_map_index[coinc_event_id]
				yield self.sngl_inspiral_index[row.event_id]

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
		ligolw_elem = newxmldoc.appendChild(ligolw.LIGO_LW())

		# when making these, we can't use .copy() method of Table
		# instances because we need to ensure we have a Table
		# subclass, not a DBTable subclass
		new_process_table = ligolw_elem.appendChild(lsctables.New(lsctables.ProcessTable, self.process_table.columnnames))
		new_process_params_table = ligolw_elem.appendChild(lsctables.New(lsctables.ProcessParamsTable, self.process_params_table.columnnames))
		new_sngl_inspiral_table = ligolw_elem.appendChild(lsctables.New(lsctables.SnglInspiralTable, self.sngl_inspiral_table.columnnames))
		new_coinc_def_table = ligolw_elem.appendChild(lsctables.New(lsctables.CoincDefTable, self.coinc_def_table.columnnames))
		new_coinc_event_table = ligolw_elem.appendChild(lsctables.New(lsctables.CoincTable, self.coinc_event_table.columnnames))
		new_coinc_inspiral_table = ligolw_elem.appendChild(lsctables.New(lsctables.CoincInspiralTable, self.coinc_inspiral_table.columnnames))
		new_coinc_event_map_table = ligolw_elem.appendChild(lsctables.New(lsctables.CoincMapTable, self.coinc_event_map_table.columnnames))
		new_time_slide_table = ligolw_elem.appendChild(lsctables.New(lsctables.TimeSlideTable, self.time_slide_table.columnnames))

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
