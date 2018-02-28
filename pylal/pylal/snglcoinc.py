# Copyright (C) 2006--2014  Kipp Cannon, Drew G. Keppel, Jolien Creighton
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
Generic coincidence engine for use with time-based event lists in LIGO
Light Weight XML documents.
"""


import bisect
try:
	from fpconst import NaN, NegInf, PosInf
except ImportError:
	# fpconst is not part of the standard library and might not be
	# available
	NaN = float("nan")
	NegInf = float("-inf")
	PosInf = float("+inf")
import itertools
import math
import numpy
import random
from scipy.constants import c as speed_of_light
import scipy.optimize
import sys
import warnings


from glue import iterutils
from glue import offsetvector
from glue import segmentsUtils
from glue.ligolw import ligolw
from glue.ligolw import array as ligolw_array
from glue.ligolw import param as ligolw_param
from glue.ligolw import lsctables
from glue.text_progress_bar import ProgressBar
from pylal import git_version
from pylal import inject
from pylal import rate


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                             Event List Interface
#
# =============================================================================
#


class EventList(list):
	"""
	A parent class for managing a list of events:  applying time
	offsets, and retrieving subsets of the list selected by time
	interval.  To be useful, this class must be subclassed with
	overrides provided for certain methods.  The only methods that
	*must* be overridden in a subclass are the _add_offset() and
	get_coincs() methods.  The make_index() method can be overridden if
	needed.  None of the other methods inherited from the list parent
	class need to be overridden, indeed they probably should not be
	unless you know what you're doing.
	"""
	def __init__(self, instrument):
		# the offset that should be added to the times of events in
		# this list when comparing to the times of other events.
		# used to implement time-shifted coincidence tests
		self.offset = lsctables.LIGOTimeGPS(0)

		# the name of the instrument from which the events in this
		# list have been taken
		self.instrument = instrument

	def make_index(self):
		"""
		Provided to allow for search-specific look-up tables or
		other indexes to be constructed for use in increasing the
		speed of the get_coincs() method.  This will be called
		after all events have been added to the list, and again if
		the list is ever modified, and before get_coincs() is ever
		called.
		"""
		pass

	def set_offset(self, offset):
		"""
		Set an offset on the times of all events in the list.
		"""
		# cast offset to LIGOTimeGPS to avoid repeated conversion
		# when applying the offset to each event.
		self.offset = lsctables.LIGOTimeGPS(offset)

	def get_coincs(self, event_a, offset_a, light_travel_time, threshold, comparefunc):
		"""
		Return a list of the events from this list that are
		coincident with event_a.

		offset_a is the time shift to be added to the time of
		event_a before comparing to the times of events in this
		list.  The offset attribute of this object will contain the
		time shift to be added to the times of the events in this
		list before comparing to event_a.  That is, the times of
		arrival of the events in this list should have (self.offset
		- offset_a) added to them before comparing to the time of
		arrival of event_a.  Or, equivalently, the time of arrival
		of event_a should have (offset_a - self.offset) added to it
		before comparing to the times of arrival of the events in
		this list.  This behaviour is to support the construction
		of time shifted coincidences.

		Because it is frequently needed by implementations of this
		method, the distance in light seconds between the two
		instruments is provided as the light_travel_time parameter.

		The threshold argument will be the thresholds appropriate
		for "instrument_a, instrument_b", in that order, where
		instrument_a is the instrument for event_a, and
		instrument_b is the instrument for the events in this
		EventList.

		comparefunc is the function to use to compare events in
		this list to event_a.
		"""
		raise NotImplementedError


class EventListDict(dict):
	"""
	A dictionary of EventList objects, indexed by instrument,
	initialized from an XML trigger table and a list of process IDs
	whose events should be included.
	"""
	def __new__(cls, *args, **kwargs):
		# wrapper to shield dict.__new__() from our arguments.
		return dict.__new__(cls)

	def __init__(self, EventListType, event_table, process_ids = None):
		"""
		Initialize a newly-created instance.  EventListType is a
		subclass of EventList (the subclass itself, not an instance
		of the subclass).  event_table is a list of events (e.g.,
		an instance of a glue.ligolw.table.Table subclass).  If the
		optional process_ids arguments is not None, then it is
		assumed to be a list or set or other thing implementing the
		"in" operator which is used to define the set of
		process_ids whose events should be considered in the
		coincidence analysis, otherwise all events are considered.
		"""
		for event in event_table:
			if (process_ids is None) or (event.process_id in process_ids):
				# FIXME:  only works when the instrument
				# name is in the "ifo" column.  true for
				# inspirals, bursts and ringdowns
				if event.ifo not in self:
					self[event.ifo] = EventListType(event.ifo)
				self[event.ifo].append(event)
		for l in self.values():
			l.make_index()

	@property
	def offsetvector(self):
		"""
		offsetvector of the offsets carried by the event lists.
		When assigning to this property, any event list whose
		instrument is not in the dictionary of offsets is not
		modified, and KeyError is raised if the offset dictionary
		contains an instrument that is not in this dictionary of
		event lists.  As a short-cut to setting all offsets to 0,
		the attribute can be deleted.
		"""
		return offsetvector.offsetvector((instrument, eventlist.offset) for instrument, eventlist in self.items())

	@offsetvector.setter
	def offsetvector(self, offsetvector):
		for instrument, offset in offsetvector.items():
			self[instrument].set_offset(offset)

	@offsetvector.deleter
	def offsetvector(self):
		for eventlist in self.values():
			eventlist.set_offset(0)


#
# =============================================================================
#
#                         Double Coincidence Iterator
#
# =============================================================================
#


def get_doubles(eventlists, comparefunc, instruments, thresholds, verbose = False):
	"""
	Given an instance of an EventListDict, an event comparison
	function, an iterable (e.g., a list) of instruments, and a
	dictionary mapping instrument pair to threshold data for use by the
	event comparison function, generate a sequence of tuples of
	mutually coincident events.

	The signature of the comparison function should be

	>>> comparefunc(event1, offset1, event2, offset2, light_travel_time, threshold_data)

	where event1 and event2 are two objects drawn from the event lists
	(of different instruments), offset1 and offset2 are the time shifts
	that should be added to the arrival times of event1 and event2
	respectively, light_travel_time is the distance in light seconds
	between the instruments from which event1 and event2 have been
	drawn, and threshold_data is the value contained in the thresholds
	dictionary for that pair of instruments.  The return value should
	be 0 (False) if the events are coincident, and non-zero otherwise
	(the behaviour of the comparison function is like a subtraction
	operator, returning 0 when the two events are "the same").

	The thresholds dictionary should look like

	>>> {("H1", "L1"): 10.0, ("L1", "H1"): -10.0}

	i.e., the keys are tuples of instrument pairs and the values
	specify the "threshold data" for that instrument pair.  The
	threshold data itself is an arbitrary Python object.  Floats are
	shown in the example above, but any Python object can be provided
	and will be passed to the comparefunc().  Note that it is assumed
	that order matters in the comparison function and so the thresholds
	dictionary must provide a threshold for the instruments in both
	orders.

	Each tuple returned by this generator will contain exactly two
	events, one from each of the two instruments in the instruments
	sequence.

	NOTE:  the instruments sequence must contain exactly two
	instruments.

	NOTE:  the order of the events in each tuple returned by this
	function is arbitrary, in particular it does not necessarily match
	the order of the instruments sequence.
	"""
	# retrieve the event lists for the requested instrument combination

	instruments = tuple(instruments)
	assert len(instruments) == 2
	for instrument in instruments:
		assert eventlists[instrument].instrument == instrument
	eventlista, eventlistb = [eventlists[instrument] for instrument in instruments]

	# insure eventlist a is the shorter of the two event lists;  record
	# the length of the shortest

	if len(eventlista) > len(eventlistb):
		eventlista, eventlistb = eventlistb, eventlista
	length = len(eventlista)

	# extract the thresholds and pre-compute the light travel time

	try:
		threshold_data = thresholds[(eventlista.instrument, eventlistb.instrument)]
	except KeyError as e:
		raise KeyError("no coincidence thresholds provided for instrument pair %s, %s" % e.args[0])
	light_travel_time = inject.light_travel_time(eventlista.instrument, eventlistb.instrument)

	# for each event in the shortest list

	for n, eventa in enumerate(eventlista):
		if verbose and not (n % 2000):
			print >>sys.stderr, "\t%.1f%%\r" % (100.0 * n / length),

		# iterate over events from the other list that are
		# coincident with the event, and return the pairs

		for eventb in eventlistb.get_coincs(eventa, eventlista.offset, light_travel_time, threshold_data, comparefunc):
			yield (eventa, eventb)
	if verbose:
		print >>sys.stderr, "\t100.0%"

	# done


#
# =============================================================================
#
#                               Time Slide Graph
#
# =============================================================================
#


class TimeSlideGraphNode(object):
	def __init__(self, offset_vector, time_slide_id = None):
		self.time_slide_id = time_slide_id
		self.offset_vector = offset_vector
		self.deltas = frozenset(offset_vector.deltas.items())
		self.components = None
		self.coincs = None
		self.unused_coincs = set()

	def name(self):
		return self.offset_vector.__str__(compact = True)

	def get_coincs(self, eventlists, event_comparefunc, thresholds, verbose = False):
		#
		# has this node already been visited?  if so, return the
		# answer we already know
		#

		if self.coincs is not None:
			if verbose:
				print >>sys.stderr, "\treusing %s" % str(self.offset_vector)
			return self.coincs

		#
		# is this a leaf node?  construct the coincs explicitly
		#

		if self.components is None:
			if verbose:
				print >>sys.stderr, "\tconstructing %s ..." % str(self.offset_vector)
			#
			# can we do it?
			#

			assert len(self.offset_vector) == 2
			avail_instruments = set(eventlists)
			offset_instruments = set(self.offset_vector)
			if not offset_instruments.issubset(avail_instruments):
				if verbose:
					print >>sys.stderr, "\twarning: do not have data for instrument(s) %s ... assuming 0 coincs" % ", ".join(offset_instruments - avail_instruments)
				self.coincs = tuple()
				return self.coincs

			#
			# apply offsets to events
			#

			if verbose:
				print >>sys.stderr, "\tapplying offsets ..."
			eventlists.offsetvector = self.offset_vector

			#
			# search for and record coincidences.  coincs is a
			# sorted tuple of event ID pairs, where each pair
			# of IDs is, itself, ordered alphabetically by
			# instrument name
			#

			if verbose:
				print >>sys.stderr, "\tsearching ..."
			# FIXME:  assumes the instrument column is named
			# "ifo".  works for inspirals, bursts, and
			# ring-downs.  note that the event order in each
			# tuple returned by get_doubles() is arbitrary so
			# we need to sort each tuple by instrument name
			# explicitly
			self.coincs = tuple(sorted((a.event_id, b.event_id) if a.ifo <= b.ifo else (b.event_id, a.event_id) for (a, b) in get_doubles(eventlists, event_comparefunc, offset_instruments, thresholds, verbose = verbose)))
			return self.coincs

		#
		# is this a head node, or some other node that magically
		# has only one component?  copy coincs from component
		#

		if len(self.components) == 1:
			if verbose:
				print >>sys.stderr, "\tgetting coincs from %s ..." % str(self.components[0].offset_vector)
			self.coincs = self.components[0].get_coincs(eventlists, event_comparefunc, thresholds, verbose = verbose)
			self.unused_coincs = self.components[0].unused_coincs

			#
			# done.  unlink the graph as we go to release
			# memory
			#

			self.components = None
			return self.coincs

		#
		# len(self.components) == 2 is impossible
		#

		assert len(self.components) > 2

		#
		# this is a regular node in the graph.  use coincidence
		# synthesis algorithm to populate its coincs
		#

		self.coincs = []

		# all coincs with n-1 instruments from the component time
		# slides are potentially unused.  they all go in, we'll
		# remove things from this set as we use them
		# NOTE:  this function call is the recursion into the
		# components to ensure they are initialized, it must be
		# executed before any of what follows
		for component in self.components:
			self.unused_coincs |= set(component.get_coincs(eventlists, event_comparefunc, thresholds, verbose = verbose))
		# of the (< n-1)-instrument coincs that were not used in
		# forming the (n-1)-instrument coincs, any that remained
		# unused after forming two compontents cannot have been
		# used by any other components, they definitely won't be
		# used to construct our n-instrument coincs, and so they go
		# into our unused pile
		for componenta, componentb in iterutils.choices(self.components, 2):
			self.unused_coincs |= componenta.unused_coincs & componentb.unused_coincs

		if verbose:
			print >>sys.stderr, "\tassembling %s ..." % str(self.offset_vector)
		# magic:  we can form all n-instrument coincs by knowing
		# just three sets of the (n-1)-instrument coincs no matter
		# what n is (n > 2).  note that we pass verbose=False
		# because we've already called the .get_coincs() methods
		# above, these are no-ops to retrieve the answers again
		allcoincs0 = self.components[0].get_coincs(eventlists, event_comparefunc, thresholds, verbose = False)
		allcoincs1 = self.components[1].get_coincs(eventlists, event_comparefunc, thresholds, verbose = False)
		allcoincs2 = self.components[-1].get_coincs(eventlists, event_comparefunc, thresholds, verbose = False)
		# for each coinc in list 0
		length = len(allcoincs0)
		for n, coinc0 in enumerate(allcoincs0):
			if verbose and not (n % 200):
				print >>sys.stderr, "\t%.1f%%\r" % (100.0 * n / length),
			# find all the coincs in list 1 whose first (n-2)
			# event IDs are the same as the first (n-2) event
			# IDs in coinc0.  note that they are guaranteed to
			# be arranged together in the list of coincs and
			# can be identified with two bisection searches
			# note:  cannot use bisect_right() because we're
			# only comparing against the first (n-2) of (n-1)
			# things in each tuple, we need to use bisect_left
			# after incrementing the last of the (n-2) things
			# by one to obtain the correct range of indexes
			coincs1 = allcoincs1[bisect.bisect_left(allcoincs1, coinc0[:-1]):bisect.bisect_left(allcoincs1, coinc0[:-2] + (coinc0[-2] + 1,))]
			# find all the coincs in list 2 whose first (n-2)
			# event IDs are the same as the last (n-2) event
			# IDs in coinc0.  note that they are guaranteed to
			# be arranged together in the list and can be
			# identified with two bisection searches
			coincs2 = allcoincs2[bisect.bisect_left(allcoincs2, coinc0[1:]):bisect.bisect_left(allcoincs2, coinc0[1:-1] + (coinc0[-1] + 1,))]
			# for each coinc extracted from list 1 above search
			# for a coinc extracted from list 2 above whose
			# first (n-2) event IDs are the last (n-2) event
			# IDs in coinc 0 and whose last event ID is the
			# last event ID in coinc 1.  when found, the first
			# ID from coinc 0 prepended to the (n-1) coinc IDs
			# from coinc 2 forms an n-instrument coinc.  how
			# this works is as follows:  coinc 0 and coinc 1,
			# both (n-1)-instrument coincs, together identify a
			# unique potential n-instrument coinc.  coinc 2's
			# role is to confirm the coincidence by showing
			# that the event from the instrument in coinc 1
			# that isn't found in coinc 0 is coincident with
			# all the other events that are in coinc 1.  if the
			# coincidence holds then that combination of event
			# IDs must be found in the coincs2 list, because we
			# assume the coincs2 list is complete  the
			# bisection search above to extract the coincs2
			# list could be skipped, but by starting with a
			# shorter list the bisection searches inside the
			# following loop are faster.
			for coinc1 in coincs1:
				i = bisect.bisect_left(coincs2, coinc0[1:] + coinc1[-1:])
				if i < len(coincs2) and coincs2[i] == coinc0[1:] + coinc1[-1:]:
					new_coinc = coinc0[:1] + coincs2[i]
					# break the new coinc into
					# (n-1)-instrument components and
					# remove them from the unused list
					# because we just used them, then
					# record the coinc and move on
					self.unused_coincs -= set(iterutils.choices(new_coinc, len(new_coinc) - 1))
					self.coincs.append(new_coinc)
		if verbose:
			print >>sys.stderr, "\t100.0%"
		# sort the coincs we just constructed by the component
		# event IDs and convert to a tuple for speed
		self.coincs.sort()
		self.coincs = tuple(self.coincs)

		#
		# done.  we won't be back here again so unlink the graph as
		# we go to release memory
		#

		self.components = None
		return self.coincs


class TimeSlideGraph(object):
	def __init__(self, offset_vector_dict, verbose = False):
		#
		# validate input
		#

		if min(len(offset_vector) for offset_vector in offset_vector_dict.values()) < 2:
			raise ValueError("offset vectors must have at least two instruments")

		#
		# populate the graph head nodes.  these represent the
		# target offset vectors requested by the calling code.
		#

		if verbose:
			print >>sys.stderr, "constructing coincidence assembly graph for %d target offset vectors ..." % len(offset_vector_dict)
		self.head = tuple(TimeSlideGraphNode(offset_vector, time_slide_id) for time_slide_id, offset_vector in sorted(offset_vector_dict.items()))

		#
		# populate the graph generations.  generations[n] is a
		# tuple of the nodes in the graph representing all unique
		# n-instrument offset vectors to be constructed as part of
		# the analysis (including normalized forms of any
		# n-instrument target offset vectors).
		#

		self.generations = {}
		n = max(len(offset_vector) for offset_vector in offset_vector_dict.values())
		self.generations[n] = tuple(TimeSlideGraphNode(offset_vector) for offset_vector in offsetvector.component_offsetvectors((node.offset_vector for node in self.head if len(node.offset_vector) == n), n))
		for n in range(n, 2, -1):	# [n, n-1, ..., 3]
			#
			# collect all offset vectors of length n that we
			# need to be able to construct
			#

			offset_vectors = [node.offset_vector for node in self.head if len(node.offset_vector) == n] + [node.offset_vector for node in self.generations[n]]

			#
			# determine the smallest set of offset vectors of
			# length n-1 required to construct the length-n
			# offset vectors, build a graph node for each of
			# the vectors of length n-1, and record the nodes
			# as the n-1'st generation
			#

			self.generations[n - 1] = tuple(TimeSlideGraphNode(offset_vector) for offset_vector in offsetvector.component_offsetvectors(offset_vectors, n - 1))

		#
		# link each n-instrument node to the n-1 instrument nodes
		# from which it will be constructed.  NOTE:  the components
		# are sorted according to the alphabetically-sorted tuples
		# of instrument names involved in each component;  this is
		# a critical part of the coincidence synthesis algorithm
		#

		for node in self.head:
			#
			# the offset vector of a head node should be found
			# directly in its generation, and it must be unique
			# or there's a bug above.  despite this, we still
			# go to the trouble of sorting to make it easy to
			# keep this code in sync with the code for other
			# graph nodes below, but the assert makes sure the
			# result contains just one entry
			#

			node.components = tuple(sorted((component for component in self.generations[len(node.offset_vector)] if node.deltas == component.deltas), key = lambda x: sorted(x.offset_vector)))
			assert len(node.components) == 1

		for n, nodes in self.generations.items():
			assert n >= 2	# failure indicates bug in code that constructed generations
			if n == 2:
				# leaf nodes have no components
				continue
			for node in nodes:
				component_deltas = set(frozenset(offset_vector.deltas.items()) for offset_vector in offsetvector.component_offsetvectors([node.offset_vector], n - 1))
				node.components = tuple(sorted((component for component in self.generations[n - 1] if component.deltas in component_deltas), key = lambda x: sorted(x.offset_vector)))

		#
		# done
		#

		if verbose:
			print >>sys.stderr, "graph contains:"
			for n in sorted(self.generations):
				print >>sys.stderr,"\t%d %d-insrument offset vectors (%s)" % (len(self.generations[n]), n, ((n == 2) and "to be constructed directly" or "to be constructed indirectly"))
			print >>sys.stderr, "\t%d offset vectors total" % sum(len(self.generations[n]) for n in self.generations)


	def get_coincs(self, eventlists, event_comparefunc, thresholds, include_small_coincs = True, verbose = False):
		if verbose:
			print >>sys.stderr, "constructing coincs for target offset vectors ..."
		for n, node in enumerate(self.head, start = 1):
			if verbose:
				print >>sys.stderr, "%d/%d: %s" % (n, len(self.head), str(node.offset_vector))
			if include_small_coincs:
				# note that unused_coincs must be retrieved
				# after the call to .get_coincs() because
				# the former is computed as a side effect
				# of the latter
				iterator = itertools.chain(node.get_coincs(eventlists, event_comparefunc, thresholds, verbose), node.unused_coincs)
			else:
				iterator = node.get_coincs(eventlists, event_comparefunc, thresholds, verbose)
			for coinc in iterator:
				yield node, coinc


	def write(self, fileobj):
		"""
		Write a DOT graph representation of the time slide graph to
		fileobj.
		"""
		print >>fileobj, "digraph \"Time Slides\" {"
		for node in itertools.chain(*self.generations.values()):
			print >>fileobj, "\t\"%s\" [shape=box];" % node.name()
			if node.components is not None:
				for component in node.components:
					print >>fileobj, "\t\"%s\" -> \"%s\";" % (component.name(), node.name())
		for node in self.head:
			print >>fileobj, "\t\"%s\" [shape=ellipse];" % node.name()
			for component in node.components:
				print >>fileobj, "\t\"%s\" -> \"%s\";" % (component.name(), node.name())
		print >>fileobj, "}"


#
# =============================================================================
#
#                              Document Interface
#
# =============================================================================
#


class CoincTables(object):
	"""
	A convenience interface to the XML document's coincidence tables,
	allowing for easy addition of coincidence events.
	"""
	def __init__(self, xmldoc):
		# find the coinc table or create one if not found
		try:
			self.coinctable = lsctables.CoincTable.get_table(xmldoc)
		except ValueError:
			self.coinctable = lsctables.New(lsctables.CoincTable)
			xmldoc.childNodes[0].appendChild(self.coinctable)
		self.coinctable.sync_next_id()

		# find the coinc_map table or create one if not found
		try:
			self.coincmaptable = lsctables.CoincMapTable.get_table(xmldoc)
		except ValueError:
			self.coincmaptable = lsctables.New(lsctables.CoincMapTable)
			xmldoc.childNodes[0].appendChild(self.coincmaptable)

		# find the time_slide table
		self.time_slide_table = lsctables.TimeSlideTable.get_table(xmldoc)
		self.time_slide_index = self.time_slide_table.as_dict()

		# cast all offsets to LIGOTimeGPS for reversable arithmetic
		# FIXME:  I believe the arithmetic in the time slide graph
		# construction can be cleaned up so that this isn't
		# required.  when that is fixed, remove this
		self.time_slide_index = dict((time_slide_id, type(offset_vector)((instrument, lsctables.LIGOTimeGPS(offset)) for instrument, offset in offset_vector.items())) for time_slide_id, offset_vector in self.time_slide_index.items())

	def append_coinc(self, process_id, time_slide_id, coinc_def_id, events):
		"""
		Takes a process ID, a time slide ID, and a list of events,
		and adds the events as a new coincidence to the coinc_event
		and coinc_map tables.

		Subclasses that wish to override this method should first
		chain to this method to construct and initialize the
		coinc_event and coinc_event_map rows.  When subclassing
		this method, if the time shifts that were applied to the
		events in constructing the coincidence are required to
		compute additional metadata, they can be retrieved from
		self.time_slide_index using the time_slide_id.
		"""
		# so we can iterate over it more than once incase we've
		# been given a generator expression.
		events = tuple(events)

		coinc = self.coinctable.RowType()
		coinc.process_id = process_id
		coinc.coinc_def_id = coinc_def_id
		coinc.coinc_event_id = self.coinctable.get_next_id()
		coinc.time_slide_id = time_slide_id
		coinc.set_instruments(None)
		coinc.nevents = len(events)
		coinc.likelihood = None
		self.coinctable.append(coinc)
		for event in events:
			coincmap = self.coincmaptable.RowType()
			coincmap.coinc_event_id = coinc.coinc_event_id
			coincmap.table_name = event.event_id.table_name
			coincmap.event_id = event.event_id
			self.coincmaptable.append(coincmap)
		return coinc


#
# =============================================================================
#
#                       Time-slideless Coinc Synthesizer
#
# =============================================================================
#


class CoincSynthesizer(object):
	"""
	Class to collect the information required to predict the rate at
	which different instrument combinations will participate in
	background coincidences, and compute various probabilities and
	rates related to the problem of doing so.
	"""

	def __init__(self, eventlists = None, segmentlists = None, delta_t = None, min_instruments = 2, abundance_rel_accuracy = 1e-4):
		"""
		eventlists is either a dictionary mapping instrument name
		to a list of the events (arbitrary objects) seen in that
		instrument or a dictionary mapping instrument name to a
		total count of events seen in that instrument (some
		features will not be available if a dictionary of counts is
		provided).  segmentlists is a glue.segments.segmentlistdict
		object describing the observation segments for each of the
		instruments.  A segment list must be provided for (at
		least) each instrument in eventlists.  delta_t is a time
		window in seconds, the light travel time between instrument
		pairs is added to this internally to set the maximum
		allowed coincidence window between a pair of instruments.
		min_instruments sets the minimum number of instruments that
		must participate in a coincidence (default is 2).

		abundance_rel_accuracy sets the fractional error tolerated
		in the Monte Carlo integrator used to estimate the relative
		abundances of the different kinds of coincs.

		Example:

		>>> from glue.segments import *
		>>> eventlists = {"H1": [0, 1, 2, 3], "L1": [10, 11, 12, 13], "V1": [20, 21, 22, 23]}
		>>> seglists = segmentlistdict({"H1": segmentlist([segment(0, 30)]), "L1": segmentlist([segment(10, 50)]), "V1": segmentlist([segment(20, 70)])})
		>>> coinc_synth = CoincSynthesizer(eventlists, seglists, 0.001)
		>>> coinc_synth.mu
		{'V1': 0.08, 'H1': 0.13333333333333333, 'L1': 0.1}
		>>> coinc_synth.tau
		{frozenset(['V1', 'H1']): 0.028287979933844225, frozenset(['H1', 'L1']): 0.011012846152223924, frozenset(['V1', 'L1']): 0.027448341016726496}
		>>> coinc_synth.rates
		{frozenset(['V1', 'H1']): 0.0006034769052553435, frozenset(['V1', 'H1', 'L1']): 1.1793108172576082e-06, frozenset(['H1', 'L1']): 0.000293675897392638, frozenset(['V1', 'L1']): 0.00043917345626762395}
		>>> coinc_synth.P_live
		{frozenset(['V1', 'H1']): 0.0, frozenset(['V1', 'H1', 'L1']): 0.25, frozenset(['H1', 'L1']): 0.25, frozenset(['V1', 'L1']): 0.5}
		>>>
		>>>
		>>> coinc_synth = CoincSynthesizer(eventlists, seglists, 0.001, min_instruments = 1)
		>>> coinc_synth.rates
		{frozenset(['V1']): 0.08, frozenset(['H1']): 0.13333333333333333, frozenset(['V1', 'H1']): 0.0006034769052553435, frozenset(['L1']): 0.1, frozenset(['V1', 'L1']): 0.00043917345626762395, frozenset(['V1', 'H1', 'L1']): 1.179508868912594e-06, frozenset(['H1', 'L1']): 0.000293675897392638}
		"""
		self.eventlists = eventlists if eventlists is not None else dict.fromkeys(segmentlists, 0) if segmentlists is not None else {}
		self.segmentlists = segmentlists if segmentlists is not None else segmentsUtils.segments.segmentlistdict()
		if set(self.eventlists) > set(self.segmentlists):
			raise ValueError("require a segmentlist for each event list")
		self.delta_t = delta_t
		if min_instruments < 1:
			raise ValueError("min_instruments must be >= 1")
		self.min_instruments = min_instruments
		self.abundance_rel_accuracy = abundance_rel_accuracy

		self.verbose = False	# turn on for diagnostics


	def reset(self):
		"""
		Reset all internally-cached data.  This method must be
		invoked if the .eventlists, .segmentlists or .delta_t
		attributes (or their contents) are modified.  This class
		relies heavily on pre-computed quantities that are derived
		from the input parameters and cached;  invoking this method
		forces the recalculation of cached data (the next time it's
		needed).  Until this method is invoked, derived data like
		coincidence window sizes and mean event rates might reflect
		the previous state of this class.
		"""
		try:
			del self._P_live
		except AttributeError:
			pass
		try:
			del self._mu
		except AttributeError:
			pass
		try:
			del self._tau
		except AttributeError:
			pass
		try:
			del self._rates
		except AttributeError:
			pass


	@property
	def all_instrument_combos(self):
		"""
		A tuple of all possible instrument combinations (as
		frozensets).
		"""
		all_instruments = tuple(self.eventlists)
		return tuple(frozenset(instruments) for n in range(self.min_instruments, len(all_instruments) + 1) for instruments in iterutils.choices(all_instruments, n))


	@property
	def P_live(self):
		"""
		Dictionary mapping instrument combination (as a frozenset)
		to fraction of the total time in which the minimum required
		number of instruments were on during which precisely that
		combination of instruments (and no other instruments) are
		on.  E.g., P_live[frozenset(("H1", "L1"))] gives the
		probability that precisely H1 and L1 are the only
		instruments operating.
		"""
		try:
			return self._P_live
		except AttributeError:
			livetime = float(abs(segmentsUtils.vote(self.segmentlists.values(), self.min_instruments)))
			all_instruments = set(self.segmentlists)
			self._P_live = dict((instruments, float(abs(self.segmentlists.intersection(instruments) - self.segmentlists.union(all_instruments - instruments))) / livetime) for instruments in self.all_instrument_combos)
			# check normalization
			total = sum(sorted(self._P_live.values()))
			assert abs(1.0 - total) < 1e-14
			for key in self._P_live:
				self._P_live[key] /= total
			# done
			return self._P_live


	@property
	def mu(self):
		"""
		Dictionary mapping instrument name to mean event rate in
		Hz.  This is a reference to an internally-cached
		dictionary.  Modifications will be retained, or an
		externally supplied dictionary can be assigned to this
		attribute to override it entirely.
		"""
		try:
			return self._mu
		except AttributeError:
			try:
				# try treating eventlists values as lists
				# and measure their lengths
				counts = dict((instrument, len(events)) for instrument, events in self.eventlists.items())
			except TypeError:
				# failed.  assume they're scalars giving
				# the number of events directly
				counts = self.eventlists
			self._mu = dict((instrument, count / float(abs(self.segmentlists[instrument]))) for instrument, count in counts.items())
			return self._mu


	@mu.setter
	def mu(self, val):
		self._mu = val
		# force re-computation of coincidence rates
		try:
			del self._rates
		except AttributeError:
			pass


	@property
	def tau(self):
		"""
		Dictionary mapping pair of instrument names (as a
		frozenset) to coincidence window in seconds.  This is a
		reference to an internally-cached dictionary.
		Modifications will be retained, or an externally supplied
		dictionary can be assigned to this attribute to override it
		entirely.
		"""
		try:
			return self._tau
		except AttributeError:
			self._tau = dict((frozenset(ab), self.delta_t + inject.light_travel_time(*ab)) for ab in iterutils.choices(tuple(self.eventlists), 2))
			return self._tau


	@tau.setter
	def tau(self, val):
		self._tau = val
		# force re-computation of coincidence rates
		try:
			del self._rates
		except AttributeError:
			pass


	@property
	def rates(self):
		"""
		Dictionary mapping instrument combination (as a frozenset)
		to mean rate in Hz at which that combination of instruments
		can be found in a coincidence under the assumption that all
		instruments are on and able to participate in coincidences.
		Corrections (e.g., based on the contents of the .P_live
		attribute) are required if that assumption does not hold.
		Note the difference between the contents of this dictionary
		and the rates of various kinds of coincidences.  For
		example, the rate for frozenset(("H1", "L1")) is the rate,
		in Hz, at which that combination of instruments
		participates in coincidences, not the rate of H1,L1
		doubles.  This is a reference to a cached internal
		dictionary.  Modifications will be retained until the
		cached data is regenerated (after .reset() is invoked or
		the .tau or .mu attributes are assigned to).
		"""
		try:
			return self._rates
		except AttributeError:
			all_instruments = set(self.mu)
			self._rates = {}
			for instruments in self.all_instrument_combos:
		# choose the instrument whose TOA forms the "epoch" of the
		# coinc.  to improve the convergence rate this should be
		# the instrument with the smallest Cartesian product of
		# coincidence windows with other instruments (so that
		# coincidence with this instrument provides the tightest
		# prior constraint on the time differences between the
		# other instruments).
				key = instruments
				anchor = min(instruments, key = lambda a: sum(math.log(self.tau[frozenset((a, b))]) for b in instruments - set([a])))
				instruments = tuple(instruments - set([anchor]))
		# compute \mu_{1} * \mu_{2} ... \mu_{N} * 2 * \tau_{12} * 2
		# * \tau_{13} ... 2 * \tau_{1N}.  this is the rate at which
		# events from instrument 1 are coincident with events from
		# all of instruments 2...N.  later, we will multiply this
		# by the probability that events from instruments 2...N
		# known to be coincident with an event from instrument 1
		# are themselves mutually coincident
				rate = self.mu[anchor]
		# the factor of 2 is because to be coincident the time
		# difference can be anywhere in [-tau, +tau], so the size
		# of the coincidence window is 2 tau
				for instrument in instruments:
					rate *= self.mu[instrument] * 2 * self.tau[frozenset((anchor, instrument))]
				if self.verbose:
					print >>sys.stderr, "%s uncorrected mean event rate = %g Hz" % (",".join(sorted(key)), rate)

		# if there are more than two instruments, correct for the
		# probability of full N-way coincidence by computing the
		# volume of the allowed parameter space by stone throwing.
		# FIXME:  it might be practical to solve this with some
		# sort of computational geometry library and convex hull
		# volume calculator.
		# FIXME:  in any case, these correction factors depend only
		# on the coincidence windows and can be computed and saved
		# when self.tau is updated allowing the rates, here, to be
		# computed very quickly if only the single-instrument
		# trigger rates have changed and not the coincidence
		# windows.
				if len(instruments) > 1:
		# for each instrument 2...N, the interval within which an
		# event is coincident with instrument 1
					windows = tuple((-self.tau[frozenset((anchor, instrument))], +self.tau[frozenset((anchor, instrument))]) for instrument in instruments)
		# pre-assemble a sequence of instrument index pairs and the
		# maximum allowed \Delta t between them to avoid doing the
		# work associated with assembling the sequence inside a
		# loop
					ijseq = tuple((i, j, self.tau[frozenset((instruments[i], instruments[j]))]) for (i, j) in iterutils.choices(range(len(instruments)), 2))
		# compute the numerator and denominator of the fraction of
		# events coincident with the anchor instrument that are
		# also mutually coincident.  this is done by picking a
		# vector of allowed \Delta ts and testing them against the
		# coincidence windows.  the loop's exit criterion is
		# arrived at as follows.  after d trials, the number of
		# successful outcomes is a binomially-distributed RV with
		# variance = d p (1 - p) <= d/4 where p is the probability
		# of a successful outcome.  we quit when the ratio of the
		# bound on the standard deviation of the number of
		# successful outcomes (d/4) to the actual number of
		# successful outcomes (n) falls below rel accuracy:
		# \sqrt{d/4} / n < rel accuracy.  note that if the true
		# probability is 0, so that n=0 identically, then the loop
		# will never terminate; from the nature of the problem we
		# know 0<p<1 so the loop will, eventually, terminate.  note
		# that if instead of using the upper bound on the variance,
		# we replace p with the estimate of p at the current
		# iteration (=n/d) and use that to estimate the variance
		# the loop can be shown to require many fewer iterations to
		# meet the desired accuracy, but that choice creates a
		# rather strong bias that, to overcome, requires some extra
		# hacks to force the loop to run for additional iterations.
		# the approach used here is much simpler.
					math_sqrt = math.sqrt
					random_uniform = random.uniform
					epsilon = self.abundance_rel_accuracy
					n, d = 0, 0
					while math_sqrt(d) >= epsilon * n:
						dt = tuple(random_uniform(*window) for window in windows)
						if all(abs(dt[i] - dt[j]) <= maxdt for i, j, maxdt in ijseq):
		# instead of adding 1 here and multiplying n by 2 in the
		# loop exit test, we increment n by 2 and then fix it
		# afterwards.
							n += 2
						d += 1
		# fix n (see above)
					n //= 2

					rate *= float(n) / float(d)
					if self.verbose:
						print >>sys.stderr, "	multi-instrument correction factor = %g" % (float(n)/float(d))
						print >>sys.stderr, "	%s mean event rate = %g Hz" % (",".join(sorted(key)), rate)

				self._rates[key] = rate
				if self.verbose:
					print >>sys.stderr, "%s mean event rate = %g Hz" % (",".join(sorted(key)), rate)

		# self._rates now contains the mean rate at which each
		# combination of instruments can be found in a coincidence
		# during the times when at least those instruments are
		# available to form coincidences.  Note:  the rate, e.g.,
		# for the combo "H1,L1" is the sum of the rate of "H1,L1"
		# doubles as well as the rate of "H1,L1,V1" triples and all
		# other higher-order coincidences in which H1 and L1
		# participate.

			# done
			return self._rates


	@property
	def mean_coinc_rate(self):
		"""
		Dictionary mapping instrument combo (as a frozenset) to the
		mean rate at which coincidences involving precisely that
		combination of instruments occur, averaged over times when
		at least the minimum required number of instruments are
		operating --- the mean rate during times when coincidences
		are possible, not the mean rate over all time.  The result
		is not cached.
		"""
		coinc_rate = dict.fromkeys(self.rates, 0.0)
		# iterate over probabilities in order for better numerical
		# accuracy
		for on_instruments, P_on_instruments in sorted(self.P_live.items(), key = lambda (ignored, P): P):
			# short cut
			if not P_on_instruments:
				continue

			# rates for instrument combinations that are
			# possible given the instruments that are on
			allowed_rates = dict((participating_instruments, rate) for participating_instruments, rate in self.rates.items() if participating_instruments <= on_instruments)

			# subtract from each rate the rate at which that
			# combination of instruments is found in (allowed)
			# higher-order coincs.  after this, allowed_rates
			# maps instrument combo to rate of coincs involving
			# exactly that combo given the instruments that are
			# on
			for key in sorted(allowed_rates, key = lambda x: len(x), reverse = True):
				allowed_rates[key] -= sum(sorted(rate for otherkey, rate in allowed_rates.items() if key < otherkey))

			for combo, rate in allowed_rates.items():
				assert rate >= 0.
				coinc_rate[combo] += P_on_instruments * rate
		return coinc_rate


	@property
	def P_instrument_combo(self):
		"""
		A dictionary mapping instrument combination (as a
		frozenset) to the probability that a background coincidence
		involves precisely that combination of instruments.  This
		is derived from the live times and the mean rates at which
		the different instrument combinations participate in
		coincidences.  The result is not cached.
		"""
		# convert rates to relative abundances
		mean_rates = self.mean_coinc_rate	# calculate once
		total_rate = sum(sorted(mean_rates.values()))
		P = dict((key, rate / total_rate) for key, rate in mean_rates.items())
		# make sure normalization is good
		total = sum(sorted(P.values()))
		assert abs(1.0 - total) < 1e-14
		for key in P:
			P[key] /= total
		return P


	@property
	def mean_coinc_count(self):
		"""
		A dictionary mapping instrument combination (as a
		frozenset) to the total number of coincidences involving
		precisely that combination of instruments expected from the
		background.  The result is not cached.
		"""
		T = float(abs(segmentsUtils.vote(self.segmentlists.values(), self.min_instruments)))
		return dict((instruments, rate * T) for instruments, rate in self.mean_coinc_rate.items())


	def instrument_combos(self):
		"""
		Generator that yields random instrument combinations (as
		frozensets) in relative abundances that match the expected
		relative abundances of background coincidences given the
		live times, mean single-instrument event rates, and
		coincidence windows.

		Example:

		>>> from glue.segments import *
		>>> eventlists = {"H1": [0, 1, 2, 3], "L1": [10, 11, 12, 13], "V1": [20, 21, 22, 23]}
		>>> seglists = segmentlistdict({"H1": segmentlist([segment(0, 30)]), "L1": segmentlist([segment(10, 50)]), "V1": segmentlist([segment(20, 70)])})
		>>> coinc_synth = CoincSynthesizer(eventlists, seglists, 0.001)
		>>> combos = coinc_synth.instrument_combos()
		>>> combos.next()	# returns a frozenset of instruments
		"""
		#
		# retrieve sorted tuple of (probability mass, instrument
		# combo) pairs.  remove instrument combos whose probability
		# mass is 0.  if no combos remain then we can't form
		# coincidences
		#

		P = tuple(sorted([mass, instruments] for instruments, mass in self.P_instrument_combo.items() if mass != 0))
		if not P:
			return

		#
		# replace the probability masses with cummulative probabilities
		#

		for i in range(1, len(P)):
			P[i][0] += P[i - 1][0]

		#
		# normalize (should be already, just be certain)
		#

		assert abs(P[-1][0] - 1.0) < 1e-14
		for i in range(len(P)):
			P[i][0] /= P[-1][0]
		assert P[-1][0] == 1.0

		#
		# generate random instrument combos
		#

		while 1:	# 1 is immutable, so faster than True
			yield P[bisect.bisect_left(P, [random.uniform(0.0, 1.0)])][1]


	def coincs(self, timefunc, allow_zero_lag = False):
		"""
		Generator to yield time shifted coincident event tuples
		without the use of explicit time shift vectors.  This
		generator can only be used if the eventlists dictionary
		with which this object was initialized contained lists of
		event objects and not merely counts of events.

		timefunc is a function for computing the "time" of an
		event, its signature should be

			t = timefunc(event)

		This function will be applied to the event objects
		contained in self.eventlists.

		If allow_zero_lag is False (the default), then only event tuples
		with no genuine zero-lag coincidences are returned, that is
		only tuples in which no event pairs would be considered to
		be coincident without time shifts applied.  Note that
		single-instrument "coincidences", if allowed, are *not*
		considered to be zero-lag coincidences.

		Example:

		>>> from glue.segments import *
		>>> eventlists = {"H1": [0, 1, 2, 3], "L1": [10, 11, 12, 13], "V1": [20, 21, 22, 23]}
		>>> seglists = segmentlistdict({"H1": segmentlist([segment(0, 30)]), "L1": segmentlist([segment(10, 50)]), "V1": segmentlist([segment(20, 70)])})
		>>> coinc_synth = CoincSynthesizer(eventlists, seglists, 0.001)
		>>> coincs = coinc_synth.coincs((lambda x: 0), allow_zero_lag = True)
		>>> coincs.next()	# returns a tuple of events
		"""
		for instruments in self.instrument_combos():
			# randomly selected events from those instruments
			instruments = tuple(instruments)
			events = tuple(random.choice(self.eventlists[instrument]) for instrument in instruments)

			# test for a genuine zero-lag coincidence among them
			if not allow_zero_lag and any(abs(ta - tb) < self.tau[frozenset((instrumenta, instrumentb))] for (instrumenta, ta), (instrumentb, tb) in iterutils.choices(zip(instruments, (timefunc(event) for event in events)), 2)):
				continue

			# return acceptable event tuples
			yield events


	def plausible_toas(self, instruments):
		"""
		Generator that yields dictionaries of random event
		time-of-arrivals for the instruments in instruments such
		that the time-of-arrivals are mutually coincident given the
		maximum allowed inter-instrument \Delta t's.

		Example:

		>>> tau = {frozenset(['V1', 'H1']): 0.028287979933844225, frozenset(['H1', 'L1']): 0.011012846152223924, frozenset(['V1', 'L1']): 0.027448341016726496}
		>>> instruments = set(("H1", "L1", "V1"))
		>>> coinc_synth = CoincSynthesizer()
		>>> coinc_synth.tau = tau	# override
		>>> toas = coinc_synth.plausible_toas(instruments)
		>>> toas.next()
		>>> toas.next()
		"""
		# this algorithm is documented in slideless_coinc_generator_rates()
		instruments = tuple(instruments)
		anchor, instruments = instruments[0], instruments[1:]
		windows = tuple((-self.tau[frozenset((anchor, instrument))], +self.tau[frozenset((anchor, instrument))]) for instrument in instruments)
		ijseq = tuple((i, j, self.tau[frozenset((instruments[i], instruments[j]))]) for (i, j) in iterutils.choices(range(len(instruments)), 2))
		while True:
			dt = tuple(random.uniform(*window) for window in windows)
			if all(abs(dt[i] - dt[j]) <= maxdt for i, j, maxdt in ijseq):
				yield dict([(anchor, 0.0)] + zip(instruments, dt))


#
# =============================================================================
#
#                                Triangulation
#
# =============================================================================
#


class TOATriangulator(object):
	"""
	Time-of-arrival triangulator.  See section 6.6.4 of
	"Gravitational-Wave Physics and Astronomy" by Creighton and
	Anderson.

	An instance of this class is a function-like object that accepts a
	tuple of event arival times and returns a tuple providing
	information derived by solving for the maximum-likelihood source
	location assuming Gaussian-distributed timing errors.
	"""
	def __init__(self, rs, sigmas, v = speed_of_light):
		"""
		Create and initialize a triangulator object.

		rs is a sequence of location 3-vectors, sigmas is a
		sequence of the timing uncertainties for those locations.
		Both sequences must be in the same order --- the first
		sigma in the sequence is interpreted as belonging to the
		first location 3-vector --- and, of course, they must be
		the same length.

		v is the speed at which the wave carrying the signals
		travels.  The rs 3-vectors carry units of distance, the
		sigmas carry units of time, v carries units of
		distance/time.  What units are used for the three is
		arbitrary, but they must be mutually consistent.  The
		default value for v in c, the speed of light, in
		metres/second, therefore the location 3-vectors should be
		given in metres and the sigmas should be given in seconds
		unless a value for v is provided with different units.

		Example:

		>>> from numpy import array
		>>> triangulator = TOATriangulator([
			array([-2161414.92636, -3834695.17889, 4600350.22664]),
			array([  -74276.0447238, -5496283.71971  ,  3224257.01744  ]),
			array([ 4546374.099   ,   842989.697626,  4378576.96241 ])
		], [
			0.005,
			0.005,
			0.005
		])

		This creates a TOATriangulator instance configured for the
		LIGO Hanford, LIGO Livingston and Virgo antennas with 5 ms
		time-of-arrival uncertainties at each location.

		Note:  rs and sigmas may be iterated over multiple times.
		"""
		assert len(rs) == len(sigmas)
		assert len(rs) >= 2

		self.rs = numpy.vstack(rs)
		self.sigmas = numpy.array(sigmas)
		self.v = v

		# sigma^-2 -weighted mean of locations
		rbar = sum(self.rs / self.sigmas[:,numpy.newaxis]**2) / sum(1 / self.sigmas**2)

		# the ith row is r - \bar{r} for the ith location
		self.R = self.rs - rbar

		# ith row is \sigma_i^-2 (r_i - \bar{r}) / c
		M = self.R / (self.v * self.sigmas[:,numpy.newaxis]**2)

		if len(rs) >= 3:
			self.U, self.S, self.VT = numpy.linalg.svd(M)

			# if the smallest singular value is less than 10^-8 * the
			# largest singular value, assume the network is degenerate
			self.singular = abs(self.S.min() / self.S.max()) < 1e-8
		else:
			# len(rs) == 2
			self.max_dt = numpy.dot(self.rs[1] - self.rs[0], self.rs[1] - self.rs[0])**.5 / self.v

	def __call__(self, ts):
		"""
		Triangulate the direction to the source of a signal based
		on a tuple of times when the signal was observed.  ts is a
		sequence of signal arrival times.  One arrival time must be
		provided for each of the observation locations provided
		when the instance was created, and the units of the arrival
		times must be the same as the units used for the sequence
		of sigmas.

		The return value is a tuple of information derived by
		solving for the maximum-likelihood source location assuming
		Gaussian-distributed timing errors.  The return value is

			(n, toa, chi2 / DOF, dt)

		where n is a unit 3-vector pointing from the co-ordinate
		origin towards the source of the signal, toa is the
		time-of-arrival of the signal at the co-ordinate origin,
		chi2 / DOF is the \chi^{2} per degree-of-freedom from to
		the arrival time residuals, and dt is the root-sum-square
		of the arrival time residuals.

		Example:

		>>> n, toa, chi2_per_dof, dt = triangulator([
			794546669.429688,
			794546669.41333,
			794546669.431885
		])
		>>> n
		array([ 0.28747132, -0.37035214,  0.88328904])
		>>> toa
		794546669.40874898
		>>> chi2_per_dof
		2.7407579727907194
		>>> dt
		0.01433725384999875
		"""
		assert len(ts) == len(self.sigmas)

		# change of t co-ordinate to avoid LIGOTimeGPS overflow
		t0 = min(ts)
		ts = numpy.array([float(t - t0) for t in ts])

		# sigma^-2 -weighted mean of arrival times
		tbar = sum(ts / self.sigmas**2) / sum(1 / self.sigmas**2)
		# the i-th element is ts - tbar for the i-th location
		tau = ts - tbar

		if len(self.rs) >= 3:
			tau_prime = numpy.dot(self.U.T, tau)[:3]

			if self.singular:
				l = 0.0
				np = tau_prime / self.S
				try:
					np[2] = math.sqrt(1.0 - np[0]**2 - np[1]**2)
				except ValueError:
					np[2] = 0.0
					np /= math.sqrt(numpy.dot(np, np))
			else:
				def n_prime(l, Stauprime = self.S * tau_prime, S2 = self.S * self.S):
					return Stauprime / (S2 + l)
				def secular_equation(l):
					np = n_prime(l)
					return numpy.dot(np, np) - 1

				# values of l that make the denominator of n'(l) 0
				lsing = -self.S * self.S
				# least negative of them is used as lower bound for
				# bisection search root finder (elements of S are
				# ordered from greatest to least, so the last
				# element of lsing is the least negative)
				l_lo = lsing[-1]

				# find a suitable upper bound for the root finder
				# FIXME:  in Jolien's original code l_hi was
				# hard-coded to 1 but we can't figure out why the
				# root must be <= 1, so I put this loop to be safe
				# but at some point it would be good to figure out
				# if 1.0 can be used because it would allow this
				# loop to be skipped
				l_hi = 1.0
				while secular_equation(l_lo) / secular_equation(l_hi) > 0:
					l_lo, l_hi = l_hi, l_hi * 2

				# solve for l
				l = scipy.optimize.brentq(secular_equation, l_lo, l_hi)

				# compute n'
				np = n_prime(l)

			# compute n from n'
			n = numpy.dot(self.VT.T, np)

			# safety check the nomalization of the result
			assert abs(numpy.dot(n, n) - 1.0) < 1e-8

			# arrival time at origin
			toa = sum((ts - numpy.dot(self.rs, n) / self.v) / self.sigmas**2) / sum(1 / self.sigmas**2)

			# chi^{2}
			chi2 = sum(((numpy.dot(self.R, n) / self.v - tau) / self.sigmas)**2)

			# root-sum-square timing residual
			dt = ts - toa - numpy.dot(self.rs, n) / self.v
			dt = math.sqrt(numpy.dot(dt, dt))
		else:
			# len(rs) == 2
			# FIXME:  fill in n and toa (is chi2 right?)
			n = numpy.zeros((3,), dtype = "double")
			toa = 0.0
			dt = max(abs(ts[1] - ts[0]) - self.max_dt, 0)
			chi2 = dt**2 / sum(self.sigmas**2)

		# done
		return n, t0 + toa, chi2 / len(self.sigmas), dt


#
# =============================================================================
#
#                     Coincidence Parameter Distributions
#
# =============================================================================
#


#
# A binning for instrument combinations
#
# FIXME:  we decided that the coherent and null stream naming convention
# would look like
#
# H1H2:LSC-STRAIN_HPLUS, H1H2:LSC-STRAIN_HNULL
#
# and so on.  i.e., the +, x and null streams from a coherent network would
# be different channels from a single instrument whose name would be the
# mash-up of the names of the instruments in the network.  that is
# inconsisntent with the "H1H2+", "H1H2-" shown here, so this needs to be
# fixed but I don't know how.  maybe it'll go away before it needs to be
# fixed.
#


def InstrumentBins(names = ("E0", "E1", "E2", "E3", "G1", "H1", "H2", "H1H2+", "H1H2-", "L1", "V1")):
	"""
	Example:

	>>> x = InstrumentBins()
	>>> x[frozenset(("H1", "L1"))]
	55
	>>> x.centres()[55]
	frozenset(['H1', 'L1'])
	"""
	return rate.HashableBins(frozenset(combo) for n in range(len(names) + 1) for combo in iterutils.choices(names, n))


#
# A class for measuring parameter distributions
#


class CoincParamsDistributions(object):
	"""
	A class for histograming the parameters of coincidences (or of
	single events).  It is assumed there is a fixed, pre-determined,
	set of parameters that one wishes to histogram, and that each
	parameter has a name.  To use this, it must be sub-classed and the
	derived class must provide dictionaries of binnings and functions
	for performing the kernel density estimation transform to obtain
	PDFs from histograms of counts.  The binnings is a dictionary
	mapping parameter names to rate.NDBins instances describing the
	binning to be used for each paramter.  The pdf_from_rates_func
	dictionary maps parameter names to functions to smooth and
	normalize bin count data into PDFs.  As a special case, a default
	function is provided and will be used for any parameters whose
	names do not appear in the pdf_from_rates_func dictionary.  The
	default function looks for a smoothing filter in the filters
	dictionary, applies it if found, then invokes the .to_pdf() method
	of the binned array object.  Subclasses must also provide a
	.coinc_params() static method that will transform a list of
	single-instrument events into a dictionary mapping paramter name to
	parameter value.

	This class maintains three sets of histograms, one set for noise
	(or "background") events, one set for signal (or "injection")
	events and one set for observed (or "zero lag") events.  The bin
	counts are floating point values (not integers).
	"""
	#
	# sub-classes may override the following
	#

	ligo_lw_name_suffix = u"pylal_snglcoinc_coincparamsdistributions"

	#
	# Default content handler for loading CoincParamsDistributions
	# objects from XML documents
	#

	class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
		pass
	ligolw_array.use_in(LIGOLWContentHandler)
	lsctables.use_in(LIGOLWContentHandler)
	ligolw_param.use_in(LIGOLWContentHandler)

	#
	# sub-classes must override the following
	#

	binnings = {}

	pdf_from_rates_func = {}
	filters = {}

	@staticmethod
	def coinc_params(*args, **kwargs):
		"""
		Given a sequence of single-instrument events (rows from an
		event table) that form a coincidence, compute and return a
		dictionary mapping parameter name to parameter values,
		suitable for being passed to one of the .add_*() methods.
		This function may return None.
		"""
		raise NotImplementedError("subclass must implement .coinc_params() method")

	#
	# begin implementation
	#

	def __init__(self, process_id = None):
		if not self.binnings:
			raise NotImplementedError("subclass must provide dictionary of binnings")
		self.zero_lag_rates = dict((param, rate.BinnedArray(binning)) for param, binning in self.binnings.items())
		self.background_rates = dict((param, rate.BinnedArray(binning)) for param, binning in self.binnings.items())
		self.injection_rates = dict((param, rate.BinnedArray(binning)) for param, binning in self.binnings.items())
		self.zero_lag_pdf = {}
		self.background_pdf = {}
		self.injection_pdf = {}
		self.zero_lag_lnpdf_interp = {}
		self.background_lnpdf_interp = {}
		self.injection_lnpdf_interp = {}
		self.process_id = process_id

	def _rebuild_interpolators(self, keys = None):
		"""
		Initialize the interp dictionaries from the discretely
		sampled PDF data.  For internal use only.
		"""
		self.zero_lag_lnpdf_interp.clear()
		self.background_lnpdf_interp.clear()
		self.injection_lnpdf_interp.clear()
		# if a specific set of keys wasn't given, do them all
		if keys is None:
			keys = set(self.zero_lag_pdf)
		# build interpolators for the requested keys
		def mkinterp(binnedarray):
			with numpy.errstate(invalid = "ignore"):
				assert not (binnedarray.array < 0.).any()
			binnedarray = binnedarray.copy()
			with numpy.errstate(divide = "ignore"):
				binnedarray.array = numpy.log(binnedarray.array)
			return rate.InterpBinnedArray(binnedarray, fill_value = NegInf)
		for key, binnedarray in self.zero_lag_pdf.items():
			if key in keys:
				self.zero_lag_lnpdf_interp[key] = mkinterp(binnedarray)
		for key, binnedarray in self.background_pdf.items():
			if key in keys:
				self.background_lnpdf_interp[key] = mkinterp(binnedarray)
		for key, binnedarray in self.injection_pdf.items():
			if key in keys:
				self.injection_lnpdf_interp[key] = mkinterp(binnedarray)

	@staticmethod
	def addbinnedarrays(rate_target_dict, rate_source_dict, pdf_target_dict, pdf_source_dict):
		"""
		For internal use.
		"""
		weight_target = {}
		weight_source = {}
		for name, binnedarray in rate_source_dict.items():
			if name in rate_target_dict:
				weight_target[name] = rate_target_dict[name].array.sum()
				weight_source[name] = rate_source_dict[name].array.sum()
				rate_target_dict[name] += binnedarray
			else:
				rate_target_dict[name] = binnedarray.copy()
		for name, binnedarray in pdf_source_dict.items():
			if name in pdf_target_dict:
				binnedarray = binnedarray.copy()
				binnedarray.array *= weight_source[name]
				pdf_target_dict[name].array *= weight_target[name]
				pdf_target_dict[name] += binnedarray
				pdf_target_dict[name].array /= weight_source[name] + weight_target[name]
			else:
				pdf_target_dict[name] = binnedarray.copy()

	def __iadd__(self, other):
		if type(other) != type(self):
			raise TypeError(other)

		self.addbinnedarrays(self.zero_lag_rates, other.zero_lag_rates, self.zero_lag_pdf, other.zero_lag_pdf)
		self.addbinnedarrays(self.background_rates, other.background_rates, self.background_pdf, other.background_pdf)
		self.addbinnedarrays(self.injection_rates, other.injection_rates, self.injection_pdf, other.injection_pdf)

		#
		# rebuild interpolators
		#

		self._rebuild_interpolators()

		#
		# done
		#

		return self

	def copy(self):
		new = type(self)(process_id = self.process_id)
		new += self
		return new

	def add_zero_lag(self, param_dict, weight = 1.0):
		"""
		Increment a bin in one or more of the observed data (or
		"zero lag") histograms by weight (default 1).  The names of
		the histograms to increment, and the parameters identifying
		the bin in each histogram, are given by the param_dict
		dictionary.
		"""
		for param, value in param_dict.items():
			try:
				self.zero_lag_rates[param][value] += weight
			except IndexError:
				# param value out of range
				pass

	def add_background(self, param_dict, weight = 1.0):
		"""
		Increment a bin in one or more of the noise (or
		"background") histograms by weight (default 1).  The names
		of the histograms to increment, and the parameters
		identifying the bin in each histogram, are given by the
		param_dict dictionary.
		"""
		for param, value in param_dict.items():
			try:
				self.background_rates[param][value] += weight
			except IndexError:
				# param value out of range
				pass

	def add_injection(self, param_dict, weight = 1.0):
		"""
		Increment a bin in one or more of the signal (or
		"injection") histograms by weight (default 1).  The names
		of the histograms to increment, and the parameters
		identifying the bin in each histogram, are given by the
		param_dict dictionary.
		"""
		for param, value in param_dict.items():
			try:
				self.injection_rates[param][value] += weight
			except IndexError:
				# param value out of range
				pass

	def default_pdf_from_rates(self, key, pdf_dict):
		"""
		For internal use by the CoincParamsDistributions class.
		"""
		binnedarray = pdf_dict[key]
		if key in self.filters:
			rate.filter_array(binnedarray.array, self.filters[key])
		binnedarray.to_pdf()

	def finish(self, verbose = False):
		"""
		Populate the discrete PDF dictionaries from the contents of
		the rates dictionaries, and then the PDF interpolator
		dictionaries from the discrete PDFs.  The raw bin counts
		from the rates dictionaries are copied verbatim, smoothed
		using the dictionary of filters carried by this class
		instance, and converted to normalized PDFs using the bin
		volumes.  Finally the dictionary of PDF interpolators is
		populated from the discretely sampled PDF data.
		"""
		#
		# convert raw bin counts into normalized PDFs
		#

		self.zero_lag_pdf.clear()
		self.background_pdf.clear()
		self.injection_pdf.clear()
		progressbar = ProgressBar(text = "Computing Parameter PDFs", max = len(self.zero_lag_rates) + len(self.background_rates) + len(self.injection_rates)) if verbose else None
		for key, (msg, rates_dict, pdf_dict) in itertools.chain(
				zip(self.zero_lag_rates, itertools.repeat(("zero lag", self.zero_lag_rates, self.zero_lag_pdf))),
				zip(self.background_rates, itertools.repeat(("background", self.background_rates, self.background_pdf))),
				zip(self.injection_rates, itertools.repeat(("injections", self.injection_rates, self.injection_pdf)))
		):
			assert numpy.isfinite(rates_dict[key].array).all() and (rates_dict[key].array >= 0).all(), "%s %s counts are not valid" % (key, msg)
			pdf_dict[key] = rates_dict[key].copy()
			try:
				pdf_from_rates_func = self.pdf_from_rates_func[key]
			except KeyError:
				pdf_from_rates_func = self.default_pdf_from_rates
			if pdf_from_rates_func is not None:
				pdf_from_rates_func(key, pdf_dict)
			if progressbar is not None:
				progressbar.increment()

		#
		# rebuild interpolators
		#

		self._rebuild_interpolators()

	def lnP_noise(self, params):
		"""
		From a parameter value dictionary as returned by
		self.coinc_params(), compute and return the natural
		logarithm of the noise probability density at that point in
		parameter space.

		The .finish() method must have been invoked before this
		method does meaningful things.  No attempt is made to
		ensure the .finish() method has been invoked nor, if it has
		been invoked, that no manipulations have occured that might
		require it to be re-invoked (e.g., the contents of the
		parameter distributions have been modified and require
		re-normalization).

		This default implementation assumes the individual PDFs
		containined in the noise dictionary are for
		statistically-independent random variables, and computes
		and returns the logarithm of their product.  Sub-classes
		that require more sophisticated calculations can override
		this method.
		"""
		__getitem__ = self.background_lnpdf_interp.__getitem__
		return sum(__getitem__(name)(*value) for name, value in params.items())

	def lnP_signal(self, params):
		"""
		From a parameter value dictionary as returned by
		self.coinc_params(), compute and return the natural
		logarithm of the signal probability density at that point
		in parameter space.

		The .finish() method must have been invoked before this
		method does meaningful things.  No attempt is made to
		ensure the .finish() method has been invoked nor, if it has
		been invoked, that no manipulations have occured that might
		require it to be re-invoked (e.g., the contents of the
		parameter distributions have been modified and require
		re-normalization).

		This default implementation assumes the individual PDFs
		containined in the signal dictionary are for
		statistically-independent random variables, and computes
		and returns the logarithm of their product.  Sub-classes
		that require more sophisticated calculations can override
		this method.
		"""
		__getitem__ = self.injection_lnpdf_interp.__getitem__
		return sum(__getitem__(name)(*value) for name, value in params.items())

	def get_xml_root(self, xml, name):
		"""
		Sub-classes can use this in their overrides of the
		.from_xml() method to find the root element of the XML
		serialization.
		"""
		name = u"%s:%s" % (name, self.ligo_lw_name_suffix)
		xml = [elem for elem in xml.getElementsByTagName(ligolw.LIGO_LW.tagName) if elem.hasAttribute(u"Name") and elem.Name == name]
		if len(xml) != 1:
			raise ValueError("XML tree must contain exactly one %s element named %s" % (ligolw.LIGO_LW.tagName, name))
		return xml[0]

	@classmethod
	def from_xml(cls, xml, name, **kwargs):
		"""
		In the XML document tree rooted at xml, search for the
		serialized CoincParamsDistributions object named name, and
		deserialize it.  The return value is a two-element tuple.
		The first element is the deserialized
		CoincParamsDistributions object, the second is the process
		ID recorded when it was written to XML.
		"""
		# find the root element of the XML serialization
		xml = self.get_xml_root(xml, name)

		# retrieve the process ID
		process_id = ligolw_param.get_pyvalue(xml, u"process_id")

		# create an instance
		self = cls(process_id = process_id, **kwargs)

		# reconstruct the BinnedArray objects
		def reconstruct(xml, prefix, target_dict):
			for name in [elem.Name.split(u":")[1] for elem in xml.childNodes if elem.Name.startswith(u"%s:" % prefix)]:
				target_dict[str(name)] = rate.BinnedArray.from_xml(xml, u"%s:%s" % (prefix, name))
		reconstruct(xml, u"zero_lag", self.zero_lag_rates)
		reconstruct(xml, u"zero_lag_pdf", self.zero_lag_pdf)
		reconstruct(xml, u"background", self.background_rates)
		reconstruct(xml, u"background_pdf", self.background_pdf)
		reconstruct(xml, u"injection", self.injection_rates)
		reconstruct(xml, u"injection_pdf", self.injection_pdf)

		#
		# rebuild interpolators
		#

		self._rebuild_interpolators()

		#
		# done
		#

		return self

	def to_xml(self, name):
		"""
		Serialize this CoincParamsDistributions object to an XML
		fragment and return the root element of the resulting XML
		tree.  The .process_id attribute of process will be
		recorded in the serialized XML, and the object will be
		given the name name.
		"""
		xml = ligolw.LIGO_LW({u"Name": u"%s:%s" % (name, self.ligo_lw_name_suffix)})
		# FIXME: remove try/except when we can rely on new-enough
		# glue to provide .from_pyvalue() class method
		try:
			xml.appendChild(ligolw_param.Param.from_pyvalue(u"process_id", self.process_id))
		except AttributeError:
			xml.appendChild(ligolw_param.from_pyvalue(u"process_id", self.process_id))
		def store(xml, prefix, source_dict):
			for name, binnedarray in sorted(source_dict.items()):
				xml.appendChild(binnedarray.to_xml(u"%s:%s" % (prefix, name)))
		store(xml, u"zero_lag", self.zero_lag_rates)
		store(xml, u"zero_lag_pdf", self.zero_lag_pdf)
		store(xml, u"background", self.background_rates)
		store(xml, u"background_pdf", self.background_pdf)
		store(xml, u"injection", self.injection_rates)
		store(xml, u"injection_pdf", self.injection_pdf)

		return xml


#
# Likelihood Ratio
#


# starting from Bayes' theorem:
#
# P(coinc is a g.w. | its parameters)
#     P(those parameters | a coinc known to be a g.w.) * P(coinc is g.w.)
#   = -------------------------------------------------------------------
#                                P(parameters)
#
#     P(those parameters | a coinc known to be a g.w.) * P(coinc is g.w.)
#   = -------------------------------------------------------------------
#     P(noise params) * P(coinc is not g.w.) + P(inj params) * P(coinc is g.w.)
#
#                       P(inj params) * P(coinc is g.w.)
#   = -------------------------------------------------------------------
#     P(noise params) * [1 - P(coinc is g.w.)] + P(inj params) * P(coinc is g.w.)
#
#                        P(inj params) * P(coinc is g.w.)
#   = ----------------------------------------------------------------------
#     P(noise params) + [P(inj params) - P(noise params)] * P(coinc is g.w.)
#
# this last form above is used below to compute the LHS
#
#          [P(inj params) / P(noise params)] * P(coinc is g.w.)
#   = --------------------------------------------------------------
#     1 + [[P(inj params) / P(noise params)] - 1] * P(coinc is g.w.)
#
#          Lambda * P(coinc is g.w.)                       P(inj params)
#   = -----------------------------------  where Lambda = ---------------
#     1 + (Lambda - 1) * P(coinc is g.w.)                 P(noise params)
#
# Differentiating w.r.t. Lambda shows the derivative is always positive, so
# thresholding on Lambda is equivalent to thresholding on P(coinc is a g.w.
# | its parameters).  The limits:  Lambda=0 --> P(coinc is a g.w. | its
# parameters)=0, Lambda=+inf --> P(coinc is a g.w. | its parameters)=1.  We
# interpret Lambda=0/0 to mean P(coinc is a g.w. | its parameters)=0 since
# although it can't be noise it's definitely not a g.w..  We do not protect
# against NaNs in the Lambda = +inf/+inf case.


class LnLikelihoodRatio(object):
	"""
	Class for computing signal hypothesis / noise hypothesis likelihood
	ratios from the measurements in a
	snglcoinc.CoincParamsDistributions instance.
	"""
	def __init__(self, coinc_param_distributions):
		self.lnP_noise = coinc_param_distributions.lnP_noise
		self.lnP_signal = coinc_param_distributions.lnP_signal

	def __call__(self, *args, **kwargs):
		"""
		Return the natural logarithm of the likelihood ratio for
		the hypothesis that the list of events are the result of a
		gravitational wave.  The likelihood ratio is the ratio
		P(inj params) / P(noise params).  The probability that the
		events are the result of a gravitiational wave is a
		monotonically increasing function of the likelihood ratio,
		so ranking events from "most like a gravitational wave" to
		"least like a gravitational wave" can be performed by
		calculating the (logarithm of the) likelihood ratios, which
		has the advantage of not requiring a prior probability to
		be provided (knowing how many gravitational waves you've
		actually detected).

		The arguments are passed verbatim to the .lnP_noise and
		.lnP_signal() methods of the
		snglcoinc.CoincParamsDistributions instance with which this
		object is associated.
		"""
		lnP_noise = self.lnP_noise(*args, **kwargs)
		lnP_signal = self.lnP_signal(*args, **kwargs)
		if math.isinf(lnP_noise) and math.isinf(lnP_signal):
			# need to handle a special case
			if lnP_noise < 0. and lnP_signal < 0.:
				# both probabilities are 0.  "correct"
				# answer is -inf, because if a candidate is
				# in a region of parameter space where the
				# probability of a signal occuring is 0
				# then there is no way it is a signal.
				# there is also, aparently, no way it's a
				# noise event, which is puzzling, but
				# that's irrelevant because we are supposed
				# to be computing something that is a
				# monotonically increasing function of the
				# probability that a candidate is a signal,
				# which is 0 in this part of the parameter
				# space.
				return NegInf
			# all remaining cases are handled correctly by the
			# expression that follows, but one still deserves a
			# warning
			if lnP_noise > 0. and lnP_signal > 0.:
				# both probabilities are +inf.  no correct
				# answer.
				warnings.warn("inf/inf encountered")
		return  lnP_signal - lnP_noise

	def samples(self, random_params_seq, sampler_coinc_params = None, **kwargs):
		"""
		Generator that yields an unending sequence of 3-element
		tuples.  Each tuple's elements are a value of the natural
		logarithm of the likelihood rato, the natural logarithm of
		the probability density of that likelihood ratio in the
		signal population, the natural logarithm of the probability
		density of that likelihood ratio in the noise population.

		random_params_seq should be a sequence (or generator) that
		yielding 2-element tuples whose first element is a choice
		of parameter values and whose second element is the natural
		logarithm of the probability density from which the
		parameters have been drawn evaluated at the parameters.

		The parameter values yielded by the random_params_seq are
		passed as the first argument, verbatim, to the .lnP_noise()
		an .lnP_signal() methods of the CoincParamsDistributions
		object with which this object is associated, followed by
		any (optional) key-word arguments.
		"""
		if sampler_coinc_params is None:
			lnP_noise_func = self.lnP_noise
			lnP_signal_func = self.lnP_signal
		else:
			lnP_noise_func = sampler_coinc_params.lnP_noise
			lnP_signal_func = sampler_coinc_params.lnP_signal
		isinf = math.isinf
		for params, lnP_params in random_params_seq:
			lnP_noise = lnP_noise_func(params, **kwargs)
			lnP_signal = lnP_signal_func(params, **kwargs)
			yield self(params, **kwargs), lnP_signal - lnP_params, lnP_noise - lnP_params
