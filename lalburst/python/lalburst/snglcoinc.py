# Copyright (C) 2006--2017  Kipp Cannon, Drew G. Keppel, Jolien Creighton
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


"""
Generic coincidence engine for use with time-based event lists in LIGO
Light Weight XML documents.
"""


from bisect import bisect_left
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
import scipy.optimize
from scipy import spatial
import sys
import warnings


import lal
from lal import rate


from glue import offsetvector
from glue import segmentsUtils
from glue.ligolw import ligolw
from glue.ligolw import array as ligolw_array
from glue.ligolw import param as ligolw_param
from glue.ligolw import lsctables
from glue.text_progress_bar import ProgressBar


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
from .git_version import date as __date__
from .git_version import version as __version__


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
	def __init__(self):
		# the offset that should be added to the times of events in
		# this list when comparing to the times of other events.
		# used to implement time-shifted coincidence tests
		self.offset = lal.LIGOTimeGPS(0)

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
		self.offset = lal.LIGOTimeGPS(offset)

	def get_coincs(self, event_a, offset_a, light_travel_time, threshold):
		"""
		Return a sequence of the events from this list that are
		coincident with event_a.  The object returned by this
		method must support being passed to bool() to determine if
		the sequence is empty.

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

	def __init__(self, EventListType, events, instruments):
		"""
		Initialize a newly-created instance.  EventListType is a
		subclass of EventList (the subclass itself, not an instance
		of the subclass).  events is an iterable of events (e.g.,
		an instance of a glue.ligolw.table.Table subclass).
		instruments is an iterable of instrument names.

		For each instrument in instruments, an instance of
		EventListType will be created, and the event objects in
		events whose .ifo attribute equals that instrument will be
		.append()'ed to that list.  If an event has a .ifo value
		that is not equal to one of the instruments, KeyError is
		raised.  It is not an error for there to be no events for a
		given instrument.

		NOTE:  the event objects in events must have a .ifo
		attribute.

		NOTE:  both the events and instruments iterables will only
		be iterated over once, so generator expressions are
		acceptable.
		"""
		for instrument in instruments:
			self[instrument] = EventListType()
		self.idx = {}
		for event in events:
			self[event.ifo].append(event)
			self.idx[id(event)] = event
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


def light_travel_time(instrument1, instrument2):
	"""
	Compute and return the time required for light to travel through
	free space the distance separating the two instruments.  The inputs
	are two instrument prefixes (e.g., "H1"), and the result is
	returned in seconds.  Note how this differs from LAL's
	XLALLightTravelTime() function, which takes two detector objects as
	input, and returns the time truncated to integer nanoseconds.
	"""
	dx = lal.cached_detector_by_prefix[instrument1].location - lal.cached_detector_by_prefix[instrument2].location
	return math.sqrt((dx * dx).sum()) / lal.C_SI


def get_doubles(eventlists, instruments, thresholds, unused):
	"""
	Given an instance of an EventListDict, an iterable (e.g., a list)
	of instruments, and a dictionary mapping instrument pair to
	threshold data for use by the event comparison test defined by the
	EventListDict (or None if this feature is not used), generate a
	sequence of tuples of Python IDs of mutually coincident events, and
	populate a set (unused) of 1-element tuples of the Python IDs of
	the events that did not participate in coincidences.

	If not set to None, the thresholds dictionary should look like

	{("H1", "L1"): 10.0, ("L1", "H1"): -10.0}

	i.e., the keys are tuples of instrument pairs and the values
	specify the "threshold data" for that instrument pair.  The
	threshold data itself is an arbitrary Python object.  Floats are
	shown in the example above, but any Python object can be provided
	and will be passed to the .get_coincs() method of an EventList
	object in the EventListDict.  Note that it is assumed that order
	matters in the selection of the threshold object function and so
	the thresholds dictionary must provide a threshold for the
	instruments in both orders.

	Each tuple returned by this generator will contain exactly two
	Python IDs, one from each of the two instruments in the instruments
	sequence.

	NOTE:  the instruments sequence must contain exactly two
	instruments;  it may be a generator, it will be iterated over only
	once.

	NOTE:  the "unused" parameter passed to this function must be a set
	or set-like object.  It will be cleared by invoking .clear(), then
	populated by invoking .update(), .add(), and .remove().

	NOTE:  the order of the IDs in each tuple returned by this function
	matches the order of the instruments sequence.
	"""
	# retrieve the event lists for the requested instrument combination

	instruments = tuple(instruments)
	assert len(instruments) == 2, "instruments must be an iterable of exactly two names, not %d" % len(instruments)
	eventlista, eventlistb = eventlists[instruments[0]], eventlists[instruments[1]]

	# choose the shorter of the two lists for the outer loop

	if len(eventlistb) < len(eventlista):
		eventlista, eventlistb = eventlistb, eventlista
		instruments = instruments[1], instruments[0]
		unswap = lambda a, b: (b, a)
	else:
		unswap = lambda a, b: (a, b)

	# extract the thresholds and pre-compute the light travel time.
	# need to do this after swapping the event lists (if they need to
	# be swapped).

	try:
		threshold_data = thresholds[instruments] if thresholds is not None else None
	except KeyError as e:
		raise KeyError("no coincidence thresholds provided for instrument pair %s, %s" % e.args[0])
	dt = light_travel_time(*instruments)

	# populate the unused set with all IDs from list B

	unused.clear()
	unused.update((id(event),) for event in eventlistb)

	# for each event in list A, iterate over events from the other list
	# that are coincident with the event, and return the pairs.  if
	# nothing is coincident with it add its ID to the set of unused
	# IDs, otherwise remove the IDs of the things that are coincident
	# with it from the set.

	offset_a = eventlista.offset
	eventlistb_get_coincs = eventlistb.get_coincs
	for eventa in eventlista:
		eventa_id = id(eventa)
		matches = eventlistb_get_coincs(eventa, offset_a, dt, threshold_data)
		if matches:
			for eventb in matches:
				eventb_id = id(eventb)
				unused.discard((eventb_id,))
				yield unswap(eventa_id, eventb_id)
		else:
			unused.add((eventa_id,))

	# done


#
# =============================================================================
#
#                               Time Slide Graph
#
# =============================================================================
#


class TimeSlideGraphNode(object):
	def __init__(self, offset_vector, time_slide_id = None, keep_unused = True):
		self.time_slide_id = time_slide_id
		self.offset_vector = offset_vector
		self.deltas = frozenset(offset_vector.deltas.items())
		self.components = None
		self.coincs = None
		self.unused_coincs = set()
		self.keep_unused = keep_unused

	@property
	def name(self):
		return self.offset_vector.__str__(compact = True)

	def get_coincs(self, eventlists, thresholds, verbose = False):
		#
		# has this node already been visited?  if so, return the
		# answer we already know
		#

		if self.coincs is not None:
			if verbose:
				print("\treusing %s" % str(self.offset_vector), file=sys.stderr)
			return self.coincs

		#
		# is this a leaf node?  construct the coincs explicitly
		#

		if self.components is None:
			if verbose:
				print("\tconstructing %s ..." % str(self.offset_vector), file=sys.stderr)

			#
			# sanity check input
			#

			assert len(self.offset_vector) == 2, "broken graph:  node with no components has %d-component offset vector, must be 2" % len(self.offset_vector)
			assert set(self.offset_vector) <= set(eventlists), "no event list for instrument(s) %s" % ", ".join(sorted(set(self.offset_vector) - set(eventlists)))

			#
			# apply offsets to events
			#

			eventlists.offsetvector = self.offset_vector

			#
			# search for and record coincidences.  coincs is a
			# sorted tuple of event ID pairs, where each pair
			# of IDs is, itself, ordered alphabetically by
			# instrument name.  note that the event order in
			# each tuple returned by get_doubles() is set by
			# the order of the instrument names passed to it,
			# which we make be alphabetical
			#

			self.coincs = tuple(sorted(get_doubles(eventlists, sorted(self.offset_vector), thresholds, self.unused_coincs)))
			if not self.keep_unused:
				self.unused_coincs.clear()
			return self.coincs

		#
		# is this a head node, or some other node that magically
		# has only one component?  copy coincs from component
		#

		if len(self.components) == 1:
			if verbose:
				print("\tcopying from %s ..." % str(self.components[0].offset_vector), file=sys.stderr)
			self.coincs = self.components[0].get_coincs(eventlists, thresholds, verbose = verbose)
			if self.keep_unused:
				# don't copy reference, always do in-place
				# manipulation of our own .unused_coincs so
				# that calling codes that might hold a
				# reference to it get the correct result
				self.unused_coincs.update(self.components[0].unused_coincs)

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
		# executed before any of what follows, and in particular it
		# must be done whether or not we'll be keeping the
		# collection of return values
		for component in self.components:
			self.unused_coincs.update(component.get_coincs(eventlists, thresholds, verbose = verbose))
		# of the (< n-1)-instrument coincs that were not used in
		# forming the (n-1)-instrument coincs, any that remained
		# unused after forming two compontents cannot have been
		# used by any other components, they definitely won't be
		# used to construct our n-instrument coincs, and so they go
		# into our unused pile
		if self.keep_unused:
			for componenta, componentb in itertools.combinations(self.components, 2):
				self.unused_coincs |= componenta.unused_coincs & componentb.unused_coincs
		else:
			self.unused_coincs.clear()

		if verbose:
			print("\tassembling %s ..." % str(self.offset_vector), file=sys.stderr)
		# magic:  we can form all n-instrument coincs by knowing
		# just three sets of the (n-1)-instrument coincs no matter
		# what n is (n > 2).  note that we pass verbose=False
		# because we've already called the .get_coincs() methods
		# above, these are no-ops to retrieve the answers again
		allcoincs0 = self.components[0].get_coincs(eventlists, thresholds, verbose = False)
		allcoincs1 = self.components[1].get_coincs(eventlists, thresholds, verbose = False)
		allcoincs2 = self.components[-1].get_coincs(eventlists, thresholds, verbose = False)
		# for each coinc in list 0
		for coinc0 in allcoincs0:
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
			coincs1 = allcoincs1[bisect_left(allcoincs1, coinc0[:-1]):bisect_left(allcoincs1, coinc0[:-2] + (coinc0[-2] + 1,))]
			# find all the coincs in list 2 whose first (n-2)
			# event IDs are the same as the last (n-2) event
			# IDs in coinc0.  note that they are guaranteed to
			# be arranged together in the list and can be
			# identified with two bisection searches
			coincs2 = allcoincs2[bisect_left(allcoincs2, coinc0[1:]):bisect_left(allcoincs2, coinc0[1:-1] + (coinc0[-1] + 1,))]
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
			# assume the coincs2 list is complete the
			# bisection search above to extract the coincs2
			# list could be skipped, but by starting with a
			# shorter list the bisection searches inside the
			# following loop are faster.
			for coinc1 in coincs1:
				confirmation = coinc0[1:] + coinc1[-1:]
				i = bisect_left(coincs2, confirmation)
				if i < len(coincs2) and coincs2[i] == confirmation:
					new_coinc = coinc0[:1] + confirmation
					# break the new coinc into
					# (n-1)-instrument components and
					# remove them from the unused list
					# because we just used them, then
					# record the coinc and move on
					self.unused_coincs.difference_update(itertools.combinations(new_coinc, len(new_coinc) - 1))
					self.coincs.append(new_coinc)
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
	def __init__(self, offset_vector_dict, min_instruments = 2, verbose = False):
		#
		# safety check input
		#

		if min_instruments < 1:
			raise ValueError("min_instruments must be >= 1: %d" % min_instruments)
		offset_vector = min(offset_vector_dict.values(), key = lambda x: len(x))
		if len(offset_vector) < 2:
			raise ValueError("encountered offset vector with fewer than 2 instruments: %s", str(offset_vector))
		if len(offset_vector) < min_instruments:
			# this test is part of the logic that ensures we
			# will only extract coincs that meet the
			# min_instruments criterion
			raise ValueError("encountered offset vector smaller than min_instruments (%d): %s", (min_instruments, str(offset_vector)))

		#
		# plays no role in the coincidence engine.  used for early
		# memory clean-up, to remove coincs from intermediate data
		# products that the calling code will not retain.  also as
		# a convenience for the calling code, implementing the
		# ubiquitous "minimum instruments" cut here in the
		# .get_coincs() method
		#

		self.min_instruments = min_instruments

		#
		# populate the graph head nodes.  these represent the
		# target offset vectors requested by the calling code.
		#

		if verbose:
			print("constructing coincidence assembly graph for %d target offset vectors ..." % len(offset_vector_dict), file=sys.stderr)
		self.head = tuple(TimeSlideGraphNode(offset_vector, time_slide_id = time_slide_id, keep_unused = len(offset_vector) > min_instruments) for time_slide_id, offset_vector in sorted(offset_vector_dict.items()))

		#
		# populate the graph generations.  generations[n] is a
		# tuple of the nodes in the graph representing all unique
		# n-instrument offset vectors to be constructed as part of
		# the analysis (including normalized forms of any
		# n-instrument target offset vectors).
		#

		self.generations = {}
		n = max(len(offset_vector) for offset_vector in offset_vector_dict.values())
		self.generations[n] = tuple(TimeSlideGraphNode(offset_vector, keep_unused = len(offset_vector) > min_instruments) for offset_vector in offsetvector.component_offsetvectors((node.offset_vector for node in self.head if len(node.offset_vector) == n), n))
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

			self.generations[n - 1] = tuple(TimeSlideGraphNode(offset_vector, keep_unused = len(offset_vector) > min_instruments) for offset_vector in offsetvector.component_offsetvectors(offset_vectors, n - 1))

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
			print("graph contains:", file=sys.stderr)
			for n in sorted(self.generations):
				print("\t%d %d-insrument offset vectors (%s)" % (len(self.generations[n]), n, ("to be constructed directly" if n == 2 else "to be constructed indirectly")), file=sys.stderr)
			print("\t%d offset vectors total" % sum(len(self.generations[n]) for n in self.generations), file=sys.stderr)


	def get_coincs(self, eventlists, thresholds, verbose = False):
		if verbose:
			print("constructing coincs for target offset vectors ...", file=sys.stderr)
		# don't do attribute look-ups in the loop
		sngl_index = eventlists.idx
		min_instruments = self.min_instruments
		for n, node in enumerate(self.head, start = 1):
			if verbose:
				print("%d/%d: %s" % (n, len(self.head), str(node.offset_vector)), file=sys.stderr)
			# the contents of .unused_coincs must be iterated
			# over after the call to .get_coincs() because the
			# former is populated as a side effect of the
			# latter, but .unused_coincs is populated in-place
			# so the following is not sensitive to the order of
			# evaluation of arguments (the .unused_coincs
			# attribute look-up can occur before or after
			# .get_coincs() is called as long as its contents
			# are not iterated over until after).
			for coinc in itertools.chain(node.get_coincs(eventlists, thresholds, verbose), node.unused_coincs):
				# we don't need to check that coinc
				# contains at least min_instruments events
				# because those that don't meet the
				# criteria are excluded during coinc
				# construction.
				yield node, tuple(sngl_index[event_id] for event_id in coinc)


	def write(self, fileobj):
		"""
		Write a DOT graph representation of the time slide graph to
		fileobj.
		"""
		print("digraph \"Time Slides\" {", file=fileobj)
		for node in itertools.chain(*self.generations.values()):
			print("\t\"%s\" [shape=box];" % node.name, file=fileobj)
			if node.components is not None:
				for component in node.components:
					print("\t\"%s\" -> \"%s\";" % (component.name, node.name), file=fileobj)
		for node in self.head:
			print("\t\"%s\" [shape=ellipse];" % node.name, file=fileobj)
			for component in node.components:
				print("\t\"%s\" -> \"%s\";" % (component.name, node.name), file=fileobj)
		print("}", file=fileobj)


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
		self.time_slide_index = dict((time_slide_id, type(offset_vector)((instrument, lal.LIGOTimeGPS(offset)) for instrument, offset in offset_vector.items())) for time_slide_id, offset_vector in self.time_slide_index.items())

	def coinc_rows(self, process_id, time_slide_id, coinc_def_id, events):
		"""
		From a process ID, a time slide ID, and a sequence of
		events (generator expressions are OK), constructs and
		initializes a coinc_event table row object and a sequence
		of coinc_event_map table row objects describing the
		coincident event.  The return value is the coinc_event row
		and a sequence of the coinc_event_map rows.

		The coinc_event is *not* assigned a coinc_event_id by this
		method.  It is expected that will be done in
		.append_coinc().  This allows sub-classes to defer the
		question of whether or not to include the coincidence in
		the search results without consuming additional IDs.

		The coinc_event row's .instruments and .likelihood
		attributes are initialized to null values.  The calling
		code should populate as needed.

		When subclassing this method, if the time shifts that were
		applied to the events in constructing the coincidence are
		required to compute additional metadata, they can be
		retrieved from self.time_slide_index using the
		time_slide_id.
		"""
		coinc = self.coinctable.RowType()
		coinc.process_id = process_id
		coinc.coinc_def_id = coinc_def_id
		coinc.coinc_event_id = None
		coinc.time_slide_id = time_slide_id
		coinc.insts = None
		coinc.likelihood = None

		coincmaps = [self.coincmaptable.RowType(
			coinc_event_id = None,
			table_name = event.event_id.table_name,
			event_id = event.event_id
		) for event in events]

		coinc.nevents = len(coincmaps)

		return coinc, coincmaps

	def append_coinc(self, coinc_event_row, coinc_event_map_rows):
		"""
		Appends the coinc_event row object and coinc_event_map row
		objects to the coinc_event and coinc_event_map tables
		respectively after assigning a coinc_event_id to the
		coincidence.  Returns the coinc_event row object.
		"""
		coinc_event_row.coinc_event_id = self.coinctable.get_next_id()
		self.coinctable.append(coinc_event_row)
		for row in coinc_event_map_rows:
			row.coinc_event_id = coinc_event_row.coinc_event_id
			self.coincmaptable.append(row)
		return coinc_event_row


#
# =============================================================================
#
#                       Time-slideless Coinc Synthesizer
#
# =============================================================================
#


class CoincRates(object):
	def __init__(self, instruments, delta_t, min_instruments = 2, abundance_rel_accuracy = 1e-4):
		"""
		Model for coincidences among noise events collected from
		several antennas.  Independent Poisson processes are
		assumed.  The coincidence window for each pair of
		instruments is (delta_t + light travel time between that
		pair).  A coincidence is assumed to require full N-way
		coincidence and require at least min_instruments to
		participate.

		Initial configuration requires some relatively expensive
		pre-calculation of internal quantities, but once
		initialized coincidence rates can be computed quickly from
		observed event rates.  Several other related quantities can
		be computed.

		Example:

		>>> coincrates = CoincRates(("H1", "L1", "V1"), 0.005, 2)
		"""
		self.instruments = frozenset(instruments)
		self.delta_t = delta_t
		self.min_instruments = min_instruments
		if self.min_instruments > len(self.instruments):
			raise ValueError("require len(instruments) >= min_instruments")
		if self.delta_t < 0.:
			raise ValueError("require delta_t >= 0.")
		if abundance_rel_accuracy <= 0.:
			raise ValueError("require abundance_rel_accuracy > 0.")

		# dictionary mapping pair of instrument names (as a
		# frozenset) to coincidence window in seconds = delta_t +
		# light travel time
		self.tau = dict((frozenset(ab), self.delta_t + light_travel_time(*ab)) for ab in itertools.combinations(tuple(self.instruments), 2))

		# for instruments {1, ..., N}, with rates \\mu_{1}, ...,
		# \\mu_{N}, the rate of coincidences is
		#
		# 	\\propto \\prod_{i} \\mu_{i}.
		#
		# the proportionality constant depends only on the
		# coincidence windows.  the following computes a dictionary
		# of these proportionality constants keyed by instrument
		# set.

		# initialize the proportionality constants
		self.rate_factors = dict.fromkeys(self.all_instrument_combos, 1.0)

		# fast-path for gstlal-inspiral pipeline:  hard-coded
		# result for H,L,V network, 5 ms coincidence window, 1 or 2
		# minimum instruments required.  computing using qhull's
		# half-plane intersection code on a machine where that's
		# available
		if self.instruments == set(("H1", "L1", "V1")) and self.delta_t == 0.005:
			if self.min_instruments == 1:
				self.rate_factors = {
					frozenset(["H1"]): 1.0,
					frozenset(["L1"]): 1.0,
					frozenset(["V1"]): 1.0,
					frozenset(["H1", "L1"]): 0.030025692304447849,
					frozenset(["H1", "V1"]): 0.064575959867688451,
					frozenset(["L1", "V1"]): 0.062896682033452986,
					frozenset(["H1", "L1", "V1"]): 0.0016876366183778862
				}
				return
			elif self.min_instruments == 2:
				self.rate_factors = {
					frozenset(["H1", "L1"]): 0.030025692304447849,
					frozenset(["H1", "V1"]): 0.064575959867688451,
					frozenset(["L1", "V1"]): 0.062896682033452986,
					frozenset(["H1", "L1", "V1"]): 0.0016876366183778862
				}
				return

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
		# the computation of a coincidence rate starts by computing
		# \mu_{1} * \mu_{2} ... \mu_{N} * 2 * \tau_{12} * 2 *
		# \tau_{13} ... 2 * \tau_{1N}.  this is the rate at which
		# events from instrument 1 are coincident with events from
		# all of instruments 2...N.  the factors of 2 are because
		# to be coincident the time difference can be anywhere in
		# [-tau, +tau], so the size of the coincidence window is 2
		# tau.  removing the factor of
		#
		#	\prod_{i} \mu_{i}
		#
		# leaves
		#
		#	\prod_{i} 2 \tau_{1i}.
		#
		# in the N-1 dimensional space defined by the time
		# differences between each instrument and the anchor
		# instrument, the coincidence windows between instruments
		# define pairs of half-space boundaries.  for the
		# coincidence windows between each instrument and the
		# anchor instrument these boundaries are perpendicular to
		# co-ordinate axes, while for other pairs of instruments
		# the coincidence windows correspond to planes angled at 45
		# degrees in various orientations.  altogether they define
		# a convex polyhedron containing the origin.
		#
		# the product of taus, above, is the volume of the
		# rectangular polyhedron defined by the anchor instrument
		# constraints alone.  it's aparent that the final quantity
		# we seek here is the volume of the convex polyhedron
		# defined by all of the constraints.  this can be computed
		# using the qhull library's half-space intersection
		# implementation, but the Python interface is not available
		# in scipy versions below 0.19, therefore we give up in
		# frustration and do the following:  we first compute the
		# volume of the rectangular polyhedron defined by the
		# anchor instrument constraints, and then use
		# stone-throwing to estimate the ratio of the desired
		# volume to that volume and multiply by that factor.
			for instrument in instruments:
				self.rate_factors[key] *= 2. * self.tau[frozenset((anchor, instrument))]
		# compute the ratio of the desired volume to that volume by
		# stone throwing.
			if len(instruments) > 1:
		# for each instrument 2...N, the interval within which an
		# event is coincident with instrument 1
				windows = tuple((-self.tau[frozenset((anchor, instrument))], +self.tau[frozenset((anchor, instrument))]) for instrument in instruments)
		# pre-assemble a sequence of instrument index pairs and the
		# maximum allowed \Delta t between them to avoid doing the
		# work associated with assembling the sequence inside a
		# loop
				ijseq = tuple((i, j, self.tau[frozenset((instruments[i], instruments[j]))]) for (i, j) in itertools.combinations(range(len(instruments)), 2))
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
		# successful outcomes (\sqrt{d/4}) to the actual number of
		# successful outcomes (n) falls below rel accuracy:
		# \sqrt{d/4} / n < rel accuracy, or
		#
		# \sqrt{d} < 2 * rel accuracy * n
		#
		# note that if the true probability is 0, so that n=0
		# identically, then the loop will never terminate; from the
		# nature of the problem we know 0<p<1 so the loop will,
		# eventually, terminate.  note that if instead of using the
		# upper bound on the variance, we replace p with the
		# estimate of p at the current iteration (=n/d) and use
		# that to estimate the variance the loop can be shown to
		# require many fewer iterations to meet the desired
		# accuracy, but that choice creates a rather strong bias
		# that, to overcome, requires some extra hacks to force the
		# loop to run for additional iterations.  the approach used
		# here is much simpler.
				math_sqrt = math.sqrt
				random_uniform = random.uniform
				two_epsilon = 2. * abundance_rel_accuracy
				n, d = 0, 0
				while math_sqrt(d) >= two_epsilon * n:
					dt = tuple(random_uniform(*window) for window in windows)
					if all(abs(dt[i] - dt[j]) <= maxdt for i, j, maxdt in ijseq):
						n += 1
					d += 1
				self.rate_factors[key] *= float(n) / float(d)

		# done computing rate_factors

		# FIXME:  commented-out implementation that requires scipy
		# >= 0.19.  saving it for later.

		# the half-space instersection code assumes constraints of
		# the form
		#
		#	A x + b <= 0,
		#
		# where A has size (n constraints x m dimensions), where
		# for N instruments m = N-1.  each coincidence window
		# between an instrument, i, and the anchor imposes two
		# constraints of the form
		#
		#	+/-t_{i} - \tau_{1i} <= 0
		#
		# for a total of 2*(N-1) constraints.  each coincidence
		# window between a pair of (non-anchor) instruments imposes
		# two constraints of the form
		#
		#	+/-(t_{i} - t_{j}) - \tau_{ij} <= 0
		#
		# for a total (N-1)*(N-2) constraints.  altogether there
		# are
		#
		#	n = (N-1)^2 + (N-1) = m * (m + 1)
		#
		# constraints
		#	if len(instruments) > 1:
		#		dimensions = len(instruments)	# anchor not included
		#		halfspaces = numpy.zeros((dimensions * (dimensions + 1), dimensions + 1), dtype = "double")
		#		# anchor constraints
		#		for i, instrument in enumerate(instruments):
		#			j = i
		#			i *= 2
		#			halfspaces[i, j] = +1.
		#			halfspaces[i + 1, j] = -1.
		#			halfspaces[i, -1] = halfspaces[i + 1, -1] = -self.tau[frozenset((anchor, instrument))]
		#		# non-anchor constraints
		#		for i, ((j1, a), (j2, b)) in enumerate(itertools.combinations(enumerate(instruments), 2), dimensions):
		#			i *= 2
		#			halfspaces[i, j1] = +1.
		#			halfspaces[i, j2] = -1.
		#			halfspaces[i + 1, j1] = -1.
		#			halfspaces[i + 1, j2] = +1.
		#			halfspaces[i, -1] = halfspaces[i + 1, -1] = -self.tau[frozenset((a, b))]
		#		# the origin is in the interior
		#		interior = numpy.zeros((len(instruments),), dtype = "double")
		#		# compute volume
		#		self.rate_factors[key] = spatial.ConvexHull(spatial.HalfspaceIntersection(halfspaces, interior).intersections).volume


	@property
	def all_instrument_combos(self):
		"""
		A tuple of all possible instrument combinations (as
		frozensets).

		Example:

		>>> coincrates = CoincRates(("H1", "L1", "V1"), 0.005, 1)
		>>> coincrates.all_instrument_combos
		(frozenset(['V1']), frozenset(['H1']), frozenset(['L1']), frozenset(['V1', 'H1']), frozenset(['V1', 'L1']), frozenset(['H1', 'L1']), frozenset(['V1', 'H1', 'L1']))
		>>> coincrates = CoincRates(("H1", "L1", "V1"), 0.005, 2)
		>>> coincrates.all_instrument_combos
		(frozenset(['V1', 'H1']), frozenset(['V1', 'L1']), frozenset(['H1', 'L1']), frozenset(['V1', 'H1', 'L1']))
		"""
		all_instruments = tuple(self.instruments)
		return tuple(frozenset(instruments) for n in range(self.min_instruments, len(all_instruments) + 1) for instruments in itertools.combinations(all_instruments, n))


	def coinc_rates(self, **rates):
		"""
		Given the event rates for a collection of instruments,
		compute the rates at which N-way coincidences occur among
		them where N >= min_instruments.  The return value is a
		dictionary whose keys are frozensets of instruments and
		whose values are the rate of coincidences for that set.

		NOTE:  the computed rates are the rates at which
		coincidences among at least those instruments occur, not
		the rate at which coincidences among exactly those
		instruments occur.  e.g., considering the H1, L1, V1
		network, for the pair H1, L1 the return value is the sum of
		the rate at which those two instruments form double
		coincidences and also the rate at which they participate in
		H1, L1, V1 triple coincidences.

		See also .strict_coinc_rates().

		Example:

		>>> coincrates = CoincRates(("H1", "L1", "V1"), 0.005, 2)
		>>> coincrates.coinc_rates(H1 = 0.001, L1 = 0.002, V1 = 0.003)
		{frozenset(['V1', 'H1']): 1.9372787960306537e-07, frozenset(['V1', 'H1', 'L1']): 1.0125819710267318e-11, frozenset(['H1', 'L1']): 6.00513846088957e-08, frozenset(['V1', 'L1']): 3.77380092200718e-07}
		>>> coincrates.coinc_rates(H1 = 0.001, L1 = 0.002, V1 = 0.002)
		{frozenset(['V1', 'H1']): 1.291519197353769e-07, frozenset(['V1', 'H1', 'L1']): 6.750546473511545e-12, frozenset(['H1', 'L1']): 6.00513846088957e-08, frozenset(['V1', 'L1']): 2.5158672813381197e-07}
		>>> coincrates.coinc_rates(H1 = 0.001, L1 = 0.002, V1 = 0.001)
		{frozenset(['V1', 'H1']): 6.457595986768845e-08, frozenset(['V1', 'H1', 'L1']): 3.3752732367557724e-12, frozenset(['H1', 'L1']): 6.00513846088957e-08, frozenset(['V1', 'L1']): 1.2579336406690598e-07}
		"""
		if set(rates) != self.instruments:
			raise ValueError("require event rates for %s, got rates for %s" % (", ".join(sorted(self.instruments)), ", ".join(sorted(rates))))
		if any(rate < 0. for rate in rates.values()):
			# they don't have to be non-negative for this
			# method to succede, but negative values are
			# nonsensical and other things will certainly fail.
			# it's best to catch and report the problem as
			# early in the code path as it can be identified,
			# so that when the problem is encountered it's
			# easier to identify the cause
			raise ValueError("rates must be >= 0")
		if max(rates.values()) * max(self.tau.values()) >= 1.:
			raise ValueError("events per coincidence window must be << 1: rates = %s, max window = %g" % (rates, max(self.tau.values())))

		# compute \mu_{1} * \mu_{2} ... \mu_{N} * FACTOR where
		# FACTOR is the previously-computed proportionality
		# constant from self.rate_factors

		coinc_rates = dict(self.rate_factors)
		for instruments in coinc_rates:
			for instrument in instruments:
				coinc_rates[instruments] *= rates[instrument]
		return coinc_rates


	def strict_coinc_rates(self, **rates):
		"""
		Given the event rates for a collection of instruments,
		compute the rates at which strict N-way coincidences occur
		among them where N >= min_instruments.  The return value is
		a dictionary whose keys are frozensets of instruments and
		whose values are the rate of coincidences for that set.

		NOTE:  the computed rates are the rates at which
		coincidences occur among exactly those instrument
		combinations, excluding the rate at which each combination
		participates in higher-order coincs.  e.g., considering the
		H1, L1, V1 network, for the pair H1, L1 the return value is
		the rate at which H1, L1 doubles occur, not including the
		rate at which the H1, L1 pair participates in H1, L1, V1
		triples.

		See also .coinc_rates().

		Example:

		>>> coincrates = CoincRates(("H1", "L1", "V1"), 0.005, 2)
		>>> coincrates.strict_coinc_rates(H1 = 0.001, L1 = 0.002, V1 = 0.003)
		{frozenset(['V1', 'H1']): 1.937177537833551e-07, frozenset(['V1', 'H1', 'L1']): 1.0125819710267318e-11, frozenset(['H1', 'L1']): 6.004125878918543e-08, frozenset(['V1', 'L1']): 3.7736996638100773e-07}
		>>> coincrates.strict_coinc_rates(H1 = 0.001, L1 = 0.002, V1 = 0.002)
		{frozenset(['V1', 'H1']): 1.2914516918890337e-07, frozenset(['V1', 'H1', 'L1']): 6.750546473511545e-12, frozenset(['H1', 'L1']): 6.004463406242219e-08, frozenset(['V1', 'L1']): 2.5157997758733847e-07}
		>>> coincrates.strict_coinc_rates(H1 = 0.001, L1 = 0.002, V1 = 0.001)
		{frozenset(['V1', 'H1']): 6.457258459445168e-08, frozenset(['V1', 'H1', 'L1']): 3.3752732367557724e-12, frozenset(['H1', 'L1']): 6.004800933565894e-08, frozenset(['V1', 'L1']): 1.2578998879366924e-07}
		"""
		# initialize from the plain coinc rates
		strict_coinc_rates = self.coinc_rates(**rates)
		# iterating over the instrument combos from the combo with
		# the most instruments to the combo with the least
		# instruments ...
		for instruments in sorted(strict_coinc_rates, reverse = True, key = lambda instruments: len(instruments)):
			# ... subtract from its rate the rate at which
			# combos containing it occur (which, because we're
			# moving from biggest combo to smallest combo, have
			# already had the rates of higher order combos
			# containing themselves subtracted)
			for key, rate in strict_coinc_rates.items():
				if instruments < key:
					strict_coinc_rates[instruments] -= rate
		# make sure this didn't produce any negative rates
		assert all(rate >= 0. for rate in strict_coinc_rates.values()), "encountered negative rate: %s" % strict_coinc_rates
		return strict_coinc_rates


	def marginalized_strict_coinc_counts(self, seglists, **rates):
		"""
		A dictionary mapping instrument combination (as a
		frozenset) to the total number of coincidences involving
		precisely that combination of instruments expected from the
		background.
		"""
		if set(seglists) != self.instruments:
			raise ValueError("require segmentlists for %s, got %s" % (", ".join(sorted(self.instruments)), ", ".join(sorted(seglists))))

		# time when exactly a given set of instruments are on
		livetime = dict((on_instruments, float(abs(seglists.intersection(on_instruments) - seglists.union(self.instruments - on_instruments)))) for on_instruments in self.all_instrument_combos)

		coinc_count = dict.fromkeys(livetime, 0.0)
		for on_instruments, T in livetime.items():
			for coinc_instruments, rate in self.strict_coinc_rates(**dict((instrument, (rate if instrument in on_instruments else 0.)) for instrument, rate in rates.items())).items():
				coinc_count[coinc_instruments] += T * rate

		return coinc_count


	def lnP_instruments(self, **rates):
		"""
		Given the event rates for a collection of instruments,
		compute the natural logarithm of the probability that a
		coincidence is found to involve exactly a given set of
		instruments.  This is equivalent to the ratios of the
		values in the dictionary returned by .strict_coinc_rates()
		to their sum.

		Raises ZeroDivisionError if all coincidence rates are 0.

		See also .strict_coinc_rates().

		Example:

		>>> coincrates = CoincRates(("H1", "L1", "V1"), 0.005, 2)
		>>> coincrates.lnP_instruments(H1 = 0.001, L1 = 0.002, V1 = 0.003)
		{frozenset(['V1', 'H1']): -1.181124067253893, frozenset(['V1', 'H1', 'L1']): -11.040192999777876, frozenset(['H1', 'L1']): -2.352494317162074, frozenset(['V1', 'L1']): -0.5143002401188091}
		"""
		strict_coinc_rates = self.strict_coinc_rates(**rates)
		total_rate = sum(strict_coinc_rates.values())
		if total_rate == 0.:
			raise ZeroDivisionError("all rates are 0")
		P_instruments = dict((instruments, rate / total_rate) for instruments, rate in strict_coinc_rates.items())
		norm = sum(sorted(P_instruments.values()))
		# safety check:  result should be nearly exactly normalized
		assert abs(1.0 - norm) < 1e-14
		return dict((instruments, (math.log(P / norm) if P else NegInf)) for instruments, P in P_instruments.items())


	def random_instruments(self, **rates):
		"""
		Generator that, given the event rates for a collection of
		instruments, yields a sequence of two-element tuples each
		containing a randomly-selected frozen set of instruments
		and the natural logarithm of the ratio of the rate at which
		that combination of instruments forms coincidences to the
		rate at which it is being yielded by this generator.

		See also .lnP_instruments().

		Example:

		>>> coincrates = CoincRates(("H1", "L1", "V1"), 0.005, 2)
		>>> x = iter(coincrates.random_instruments(H1 = 0.001, L1 = 0.002, V1 = 0.003))
		>>> x.next()	# doctest: +SKIP
		(frozenset(['H1', 'L1']), -3.738788683913535)
		"""
		# guaranteed non-empty
		lnP_instruments = self.lnP_instruments(**rates)
		lnN = math.log(len(lnP_instruments))
		results = tuple((instruments, lnP - lnN) for instruments, lnP in lnP_instruments.items())
		choice = random.choice
		while 1:
			yield choice(results)


	def plausible_toas(self, instruments):
		"""
		Generator that yields dictionaries of random event
		time-of-arrival offsets for the given instruments such that
		the time-of-arrivals are mutually coincident given the
		maximum allowed inter-instrument \\Delta t's.  The values
		returned are offsets, and would need to be added to some
		common time to yield absolute arrival times.

		Example:

		>>> coincrates = CoincRates(("H1", "L1", "V1"), 0.005, 2)
		>>> x = iter(coincrates.plausible_toas(("H1", "L1")))
		>>> x.next()	# doctest: +SKIP
		{'H1': 0.0, 'L1': -0.010229226372297711}
		"""
		instruments = tuple(instruments)
		if set(instruments) > self.instruments:
			raise ValueError("not configured for %s" % ", ".join(sorted(set(instruments) - self.instruments)))
		if len(instruments) < self.min_instruments:
			raise ValueError("require at least %d instruments, got %d" % (self.min_instruments, len(instruments)))
		anchor, instruments = instruments[0], instruments[1:]
		anchor_offset = ((anchor, 0.0),)	 # don't build inside loop
		uniform = random.uniform
		windows = tuple((instrument, -self.tau[frozenset((anchor, instrument))], +self.tau[frozenset((anchor, instrument))]) for instrument in instruments)
		ijseq = tuple((i, j, self.tau[frozenset((instruments[i], instruments[j]))]) for (i, j) in itertools.combinations(range(len(instruments)), 2))
		while 1:
			dt = tuple((instrument, uniform(lo, hi)) for instrument, lo, hi in windows)
			if all(abs(dt[i][1] - dt[j][1]) <= maxdt for i, j, maxdt in ijseq):
				yield dict(anchor_offset + dt)


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
	def __init__(self, rs, sigmas, v = lal.C_SI):
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
		...	array([-2161414.92636, -3834695.17889, 4600350.22664]),
		...	array([  -74276.0447238, -5496283.71971  ,  3224257.01744  ]),
		...	array([ 4546374.099   ,   842989.697626,  4378576.96241 ])
		... ], [
		...	0.005,
		...	0.005,
		...	0.005
		... ])
		...
		>>>

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
		chi2 / DOF is the \\chi^{2} per degree-of-freedom from to
		the arrival time residuals, and dt is the root-sum-square
		of the arrival time residuals.

		Example:

		>>> from numpy import array
		>>> triangulator = TOATriangulator([
		...	array([-2161414.92636, -3834695.17889, 4600350.22664]),
		...	array([  -74276.0447238, -5496283.71971  ,  3224257.01744  ]),
		...	array([ 4546374.099   ,   842989.697626,  4378576.96241 ])
		... ], [
		...	0.005,
		...	0.005,
		...	0.005
		... ])
		...
		>>> n, toa, chi2_per_dof, dt = triangulator([
		...	794546669.429688,
		...	794546669.41333,
		...	794546669.431885
		... ])
		...
		>>> n
		array([ 0.28747132, -0.37035214,  0.88328904])
		>>> print(toa)
		794546669.409
		>>> print(chi2_per_dof)
		2.74075797279
		>>> print(dt)
		0.01433725385
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

				# values of l that make the denominator of
				# n'(l) 0
				lsing = -self.S * self.S
				# least negative of them is used as lower
				# bound for bisection search root finder
				# (elements of S are ordered from greatest
				# to least, so the last element of lsing is
				# the least negative)
				l_lo = lsing[-1]

				# find a suitable upper bound for the root
				# finder FIXME:  in Jolien's original code
				# l_hi was hard-coded to 1 but we can't
				# figure out why the root must be <= 1, so
				# I put this loop to be safe but at some
				# point it would be good to figure out if
				# 1.0 can be used because it would allow
				# this loop to be skipped
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
	return rate.HashableBins(frozenset(combo) for n in range(len(names) + 1) for combo in itertools.combinations(names, n))


#
# Base class for parameter distribution densities for use in log likelihood
# ratio ranking statistics
#


class LnLRDensity(object):
	"""
	Base class for parameter distribution densities for use in log
	likelihood ratio ranking statistics.  Generally several instances
	of (subclasses of) this will be grouped together to construct a log
	likelihood ratio class for use as a ranking statistic in a
	trigger-based search.  For example, as a minimum one would expect
	one instance for the numerator and another for the denominator, but
	additional ones might be included in a practical ranking statistic
	implementation, for example a third might be used for storing a
	histogram of the candidates observed in a search.

	Typically, the ranking statistic implementation will provide a
	function to transform a candidate to the arguments to use with the
	.__call__() implementation, and so in this way a LnLRDensity object
	is generally only meaningful in the context of the ranking
	statistic class for which it has been constructed.
	"""
	def __call__(self, *args, **kwargs):
		"""
		Evaluate.  Return the natural logarithm of the density
		evaluated at the given parameters.
		"""
		raise NotImplementedError

	def __iadd__(self, other):
		"""
		Marginalize the two densities.
		"""
		raise NotImplementedError

	def increment(self, *args, **kwargs):
		"""
		Increment the counts defining this density at the given
		parameters.
		"""
		raise NotImplementedError

	def copy(self):
		"""
		Return a duplicate copy of this object.
		"""
		raise NotImplementedError

	def finish(self):
		"""
		Ensure all internal densities are normalized, and
		initialize interpolator objects as needed for smooth
		evaluation.  Must be invoked before .__call__() will yield
		sensible results.

		NOTE:  for some implementations this operation will
		irreversibly alter the contents of the counts array, for
		example often this operation will involve the convolution
		of the counts with a density estimation kernel.  If it is
		necessary to preserve a pristine copy of the counts data,
		use the .copy() method to obtain a copy of the data, first,
		and then .finish() the copy.
		"""
		raise NotImplementedError

	def samples(self):
		"""
		Generator returning a sequence of parameter values drawn
		from the distribution density.  Some subclasses might
		choose not to implement this, and those that do might
		choose to use an MCMC-style sample generator and so the
		samples should not be assumed to be statistically
		independent.
		"""
		raise NotImplementedError

	def to_xml(self, name):
		"""
		Serialize to an XML fragment and return the root element of
		the resulting XML tree.

		Subclasses must chain to this method, then customize the
		return value as needed.
		"""
		return ligolw.LIGO_LW({u"Name": u"%s:lnlrdensity" % name})

	@classmethod
	def get_xml_root(cls, xml, name):
		"""
		Sub-classes can use this in their overrides of the
		.from_xml() method to find the root element of the XML
		serialization.
		"""
		name = u"%s:lnlrdensity" % name
		xml = [elem for elem in xml.getElementsByTagName(ligolw.LIGO_LW.tagName) if elem.hasAttribute(u"Name") and elem.Name == name]
		if len(xml) != 1:
			raise ValueError("XML tree must contain exactly one %s element named %s" % (ligolw.LIGO_LW.tagName, name))
		return xml[0]

	@classmethod
	def from_xml(cls, xml, name):
		"""
		In the XML document tree rooted at xml, search for the
		serialized LnLRDensity object named name, and deserialize
		it.  The return value is the deserialized LnLRDensity
		object.
		"""
		# Generally implementations should start with something
		# like this:
		#xml = cls.get_xml_root(xml, name)
		#self = cls()
		#return self
		raise NotImplementedError


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

	@classmethod
	def get_xml_root(cls, xml, name):
		"""
		Sub-classes can use this in their overrides of the
		.from_xml() method to find the root element of the XML
		serialization.
		"""
		name = u"%s:%s" % (name, cls.ligo_lw_name_suffix)
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
		xml = cls.get_xml_root(xml, name)

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
		xml.appendChild(ligolw_param.Param.from_pyvalue(u"process_id", self.process_id))
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
#     P(those parameters | coinc is g.w.) * P(coinc is g.w.)
#   = ------------------------------------------------------
#                         P(parameters)
#
#               P(those parameters | coinc is g.w.) * P(coinc is g.w.)
#   = -------------------------------------------------------------------------
#     P(noise params) * P(coinc is not g.w.) + P(inj params) * P(coinc is g.w.)
#
#                        P(inj params) * P(coinc is g.w.)
#   = ---------------------------------------------------------------------------
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


class LnLikelihoodRatioMixin(object):
	"""
	Mixin class to provide the standard log likelihood ratio methods.
	Intended to be added to the parent classes of a ranking statistic
	class defining .numerator and .denominator attributes that are both
	instances of (subclasses of) the LnLRDensity class.  The ranking
	statistic class will then acquire a .__call__() method allowing it
	to be used as a log likelihood ratio function, and also a
	.ln_lr_samples() method providing importance-weighted sampling of
	the log likelihood ratio distribution in the signal and noise
	(numerator and denominator) populations.
	"""
	def __call__(self, *args, **kwargs):
		"""
		Return the natural logarithm of the likelihood ratio for
		the given parameters.  The likelihood ratio is P(params |
		signal) / P(params | noise).  The probability that the
		events are the result of a gravitiational wave is a
		monotonically increasing function of the likelihood ratio,
		so ranking events from "most like a gravitational wave" to
		"least like a gravitational wave" can be performed by
		calculating the (logarithm of the) likelihood ratios.

		The arguments are passed verbatim to the .__call__()
		methods of the .numerator and .denominator attributes of
		self.
		"""
		lnP_signal = self.numerator(*args, **kwargs)
		lnP_noise = self.denominator(*args, **kwargs)
		if math.isinf(lnP_noise) and math.isinf(lnP_signal):
			# need to handle a special case
			if lnP_noise < 0. and lnP_signal < 0.:
				# both probabilities are 0.  "correct"
				# answer is -inf, because if a candidate is
				# in a region of parameter space where the
				# probability of a signal occuring is 0
				# then it is not a signal.  is it also,
				# aparently, not noise, which is curious
				# but irrelevant because we are seeking a
				# result that is a monotonically increasing
				# function of the probability that a
				# candidate is a signal, which is
				# impossible in this part of the parameter
				# space.
				return NegInf
			# all remaining cases are handled correctly by the
			# expression that follows, but one still deserves a
			# warning
			if lnP_noise > 0. and lnP_signal > 0.:
				# both probabilities are +inf.  no correct
				# answer.  NaN will be returned in this
				# case, and it helps to have a record in
				# the log of why that happened.
				warnings.warn("inf/inf encountered")
		return  lnP_signal - lnP_noise

	def ln_lr_samples(self, random_params_seq, sampler_coinc_params = None):
		"""
		Generator that yields an unending sequence of 3-element
		tuples.  Each tuple's elements are a value of the natural
		logarithm of the likelihood rato, the natural logarithm of
		the relative frequency of occurance of that likelihood
		ratio in the signal population corrected for the relative
		frequency at which the sampler is yielding that value, and
		the natural logarithm of the relative frequency of
		occurance of that likelihood ratio in the noise population
		similarly corrected for the relative frequency at which the
		sampler is yielding that value.  The intention is for the
		return values to be added to histograms using the given
		probability densities as weights, i.e., the two relative
		frequencies give the number of times one should consider
		this one draw of log likelihood ratio to have occured in
		the two populations.

		random_params_seq is a sequence (generator is OK) yielding
		3-element tuples whose first two elements provide the *args
		and **kwargs values passed to the numerator and denominator
		density functions, and whose thrid element is the natural
		logarithm of the probability density from which the
		parameters have been drawn evaluated at the parameters.

		On each iteration, the *args and **kwargs values yielded by
		random_params_seq is passed to our own .__call__() method
		to evalute the log likelihood ratio at that choice of
		parameter values.  If sampler_coinc_params is None the
		parameters are also passed to the .__call__() mehods of the
		.numerator and .denominator attributes of self to obtain
		the signal and noise population densities at those
		parameters.  If sample_coinc_params is not None then,
		instead, the parameters are passed to the .__call__()
		methods of its .numerator and .denominator attributes.

		If histograming the results as described above, the effect
		is to draw paramter values from the signal and noise
		populations defined by sampler_coinc_params' PDFs but with
		log likelihood ratios evaluted using our own PDFs.
		"""
		if sampler_coinc_params is None:
			lnP_signal_func = self.numerator
			lnP_noise_func = self.denominator
		else:
			lnP_signal_func = sampler_coinc_params.numerator
			lnP_noise_func = sampler_coinc_params.denominator
		for args, kwargs, lnP_params in random_params_seq:
			lnP_signal = lnP_signal_func(*args, **kwargs)
			lnP_noise = lnP_noise_func(*args, **kwargs)
			yield self(*args, **kwargs), lnP_signal - lnP_params, lnP_noise - lnP_params
