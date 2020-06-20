# Copyright (C) 2006--2018  Kipp Cannon, Drew G. Keppel, Jolien Creighton
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
Generic time interval coincidence engine.
"""


from bisect import bisect_left
import collections
try:
	from fpconst import NegInf
except ImportError:
	# fpconst is not part of the standard library and might not be
	# available
	NegInf = float("-inf")
import itertools
import math
import numpy
import random
import scipy.optimize
from scipy import spatial
import sys
try:
	# python 2
	from UserDict import DictMixin
	class UserDict(object, DictMixin):
		pass
except ImportError:
	# python 3
	from collections import UserDict
import warnings


from ligo.lw import ligolw
from ligo.lw import lsctables
from ligo.lw.utils import coincs as ligolw_coincs
from glue.text_progress_bar import ProgressBar
import lal
from ligo.segments import NegInfinity
from . import offsetvector


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
from .git_version import date as __date__
from .git_version import version as __version__


#
# =============================================================================
#
#                    Streaming Double Coincidence Generator
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


class singlesqueue(object):
	"""
	Buffer a partially time-ordered stream of events, sort them into
	time order, maintain a record of the time up-to which the list of
	events is complete, and provide events in time order to coincidence
	tests.

	Searches must define a class subclassed from this that defines an
	override of the .event_time() method.  See coincgen_doubles for
	information on what to do with the subclass.
	"""
	@staticmethod
	def event_time(event):
		"""
		Override with a method to extract the "time" of the given
		event.  The "time" object that is returned should be a
		lal.LIGOTimeGPS object.  The object must support arithmetic
		and comparison with python float objects and
		lal.LIGOTimeGPS objects, and support comparison with
		ligo.segments.PosInfinity and ligo.segments.NegInfinity,
		and it must have a .__dict__ or other mechanism allowing a
		an additional attribute named .event to be set on the
		object.  float is not suitable.
		"""
		raise NotImplementedError


	def queueentry_from_event(self, event):
		"""
		For internal use.
		"""
		entry = self.event_time(event)
		entry.event = event
		return entry

	def __init__(self, coinc_window):
		# using .event_time() to define the times of events, this
		# is the minimum time that must separate events in this
		# queue from other queues for them to be certain to fail
		# coincidence.  it is not necessary for this to be exactly
		# the coincidence window, it is sufficient to bound it from
		# above, but performance is improved by making this close
		# to, or ideally exactly, the coincidence window.  NOTE:
		# this must include light travel time, and any additional
		# window required to allow for time shift offsets.
		self.coinc_window = coinc_window
		# using .event_time() to define the times of events, the
		# list of events is complete upto but not including
		# .t_valid.  we can't use fpconst.NegInf because
		# LIGOTimeGPS objects refuse to be compared to it.  the
		# NegInfinity object from the segments library, however, is
		# compatible with LIGOTimeGPS.
		self.t_complete = NegInfinity
		# queues of events.  the queues' contents are time-ordered,
		# they contain the events upto .t_complete and from
		# .t_complete on, respectively.
		self.complete = []
		self.incomplete = []
		# id() --> event mapping for the contents of queue.  sets
		# will be used to track the status of events, e.g. which
		# have and haven't been used to form candidates.  we don't
		# require the event objects to be hashable and suitable for
		# inclusion in sets, instead we put their Python IDs into
		# the sets, but then we need a way to turn an ID back into
		# an event, and that's what this provides
		self.index = {}

	@property
	def age(self):
		"""
		Using .event_time() to define the times of events, the time
		of the oldest event in the queue or self.t_complete if the
		queue is empty.
		"""
		return self.complete[0] if self.complete else self.t_complete

	@property
	def t_coinc_complete(self):
		"""
		Using .event_time() to define the times of events, events
		whose times are less than this are sufficiently far behind
		.t_complete that all coincidences that can be formed with
		those events can now be formed (assuming all queues
		contributing events to the coincidences share the same
		.t_complete).  Events with times equal to or later than
		this are too close to .t_complete to guarantee that all
		coincidences in which they will participate can now be
		formed.  After processing the current contents of the
		queue, it could be flushed up to this time without loss of
		information.
		"""
		return self.t_complete - self.coinc_window

	def push(self, events, t_complete):
		"""
		Add events to the queue.  Mark the queue complete up to
		t_complete, meaning you promise no further events will be
		added to the queue earlier than t_complete.

		NOTE:  events will be iterated over multiple times.  It may
		not be a generator.
		"""
		if t_complete < self.t_complete:
			raise ValueError("t_complete has gone backwards:  last was %s, new is %s" % (self.t_complete, t_complete))

		if events:
			# add the new events to the ID index
			self.index.update((id(event), event) for event in events)

			# construct queue entries from the new events and
			# insort with current "incomplete" queue
			self.incomplete += map(self.queueentry_from_event, events)
			self.incomplete.sort()
			if self.incomplete[0] < self.t_complete:
				raise ValueError("t_complete violation: earliest event is %s, previous t_complete was %s" % (self.incomplete[0], self.t_complete))

		# move events preceding the new t_complete from the
		# incomplete to the complete queue.
		i = bisect_left(self.incomplete, t_complete)
		self.complete += self.incomplete[:i]
		del self.incomplete[:i]

		# update the marker labelling time up to which the event
		# list is complete
		self.t_complete = t_complete

	def pull(self, t):
		"""
		Returns a pair of tuples.  The first contains the events
		from the queue upto t;  the second contains any additional
		events in the queue that can form coincidences with events
		upto t, i.e., the events from the queue at or following t
		up to (not including) t + coinc_window.  The events
		returned in the first of the two tuples are removed from
		the queue.  The times of events are defined by
		.event_time(), and both tuples are time-ordered.

		If t is None, then the queue is flushed.  All remaining
		events are pulled from the queue, .t_complete is reset to
		-inf and all other internal state is reset.  After calling
		.pull() with t = None, the object's state is equivalent to
		its initial state and it is ready to process a new stream
		of events.
		"""
		if t is None:
			events = tuple(entry.event for entry in self.complete) + tuple(entry.event for entry in self.incomplete)
			del self.complete[:]
			del self.incomplete[:]
			self.index.clear()
			self.t_complete = NegInfinity
			return events, ()

		if t > self.t_coinc_complete:
			raise ValueError("pull to %g exceeds time to which queue is coinc-complete, %g" % (t, self.t_coinc_complete))

		# these events will never be used again.  remove them from
		# the queue, and remove their IDs from the index
		i = bisect_left(self.complete, t)
		events = tuple(entry.event for entry in self.complete[:i])
		del self.complete[:i]
		for event in events:
			self.index.pop(id(event))
		# collect other events that might be used again but that
		# can form coincidences with things up to t
		t += self.coinc_window
		if t <= self.t_complete:
			other_events = tuple(entry.event for entry in self.complete[:bisect_left(self.complete, t)])
		else:
			other_events = tuple(entry.event for entry in self.complete) + tuple(entry.event for entry in self.incomplete[:bisect_left(self.incomplete, t)])
		return events, other_events


class multidict(UserDict):
	"""
	Read-only dictionary view into a collection of dictionaries.

	Example:

	>>> x = {"a": 10., "b": 100.}
	>>> y = {"c": 1000., "d": 10000.}
	>>> z = multidict(x, y)
	>>> z["a"]
	10.0
	>>> z["c"]
	1000.0
	>>> "d" in z
	True
	>>> "e" in z
	False
	"""
	def __init__(self, *dicts):
		self.dicts = dicts
		if not self.dicts:
			raise ValueError("len(dicts) must be > 0")

	def __getitem__(self, x):
		for d in self.dicts:
			if x in d:
				return d[x]
		raise KeyError(x)

	def __nonzero__(self):
		return any(self.dicts)

	def __contains__(self, x):
		return any(x in d for d in self.dicts)

	def __len__(self):
		return sum(len(d) for d in self.dicts)

	def __iter__(self):
		return itertools.chain(*(iter(d) for d in self.dicts))

	def iteritems(self):
		return itertools.chain(*(d.iteritems() for d in self.dicts))

	def keys(self):
		return list(self)


class coincgen_doubles(object):
	"""
	Using a pair of singlesqueue objects, constructs pairs of
	coincident events from streams of partially time-ordered events.

	Searches must subclass this.  The .singlesqueue class attribute
	must be set to the search-specific singlesqueue implementation to
	be used.  The internally-defined .get_coincs class must be
	overridden with an implementation that provides the required
	.__init__() and .__call__() methods.
	"""
	# subclass must override this attribute with a reference to an
	# implementation of singlesqueue that implements the .event_time()
	# method.
	singlesqueue = singlesqueue

	class get_coincs(object):
		"""
		Override with a class whose .__init__() accepts a tuple of
		events, and whose .__call__ behaves as described below.
		The events will be provided to .__init__() in time order.

		It is not required that the implementation be subclassed
		from this.  This placeholder implementation is merely
		provided to document the required interface.
		"""
		def __init__(self, events):
			"""
			Prepare to search a collection of events for
			coincidences with other single events.  events is a
			time-ordered iterable of events.  It is recommended
			that any additional indexing required to improve
			search performance be performed in this method.
			"""
			raise NotImplementedError

		def __call__(self, event_a, offset_a, coinc_window):
			"""
			Return a sequence of the events from those passed
			to .__init__() that are coincident with event_a.
			The sequence need not be in time order.  The object
			returned by this method must be iterable and
			support being passed to bool() to test if it is
			empty.

			offset_a is the time shift to be added to the time
			of event_a before comparing to the times of events
			passed to .__init__().  This behaviour is to
			support the construction of time shifted
			coincidences.

			coinc_window is the maximum time, in seconds,
			separating events from the shifted time of event_a.
			This is the value passed to coincgen_doubles()'s
			.__init__() method with the light travel time
			between this event list's detector and the detector
			of event_a added.
			"""
			raise NotImplementedError

	def __init__(self, offset_vector, coinc_window):
		"""
		offset_vector must be a two-instrument offset vector.  This
		sets which two instruments' events are to be processed by
		this object, and the offsets that should be applied to
		their events when searching for coincidences.  coinc_window
		is the coincidence window in seconds, *not including* the
		light travel time between the instruments.  For covenience,
		the light travel time will be added internally by this
		object.
		"""
		if len(offset_vector) != 2:
			raise ValueError("offset_vector must name exactly 2 instruments (got %d)" % len(offset_vector))
		if coinc_window < 0.:
			raise ValueError("coinc_window must be non-negative (got %g)" % coinc_window)
		self.offset_vector = offset_vector
		# add the light travel time to coinc_window
		self.coinc_window = coinc_window + light_travel_time(*offset_vector)
		# FIXME:  the latency can be reduced by teaching
		# singlesqueue about the asymmetry of the offsetvector
		self.queues = dict((instrument, self.singlesqueue(coinc_window + abs(offset_vector))) for instrument in offset_vector)
		# view into the id() --> event indexes of the queues
		self.index = multidict(*(queue.index for queue in self.queues.values()))
		# Python id()s of events currently in the queues that have
		# been reported in coincidences
		self.used = set()

	@property
	def age(self):
		"""
		The earliest of the internal queues' .age.
		"""
		return min(queue.age for queue in self.queues.values())

	@property
	def t_complete(self):
		"""
		The earliest of the internal queues' .t_complete.
		"""
		return min(queue.t_complete for queue in self.queues.values())

	@property
	def t_coinc_complete(self):
		"""
		The earliest of the internal queues' .t_coinc_complete.
		"""
		return min(queue.t_coinc_complete for queue in self.queues.values())

	def push(self, instrument, events, t_complete):
		"""
		Push new events from some instrument into the internal
		queues.  The iterable of events need not be time-ordered.
		With self.singlesqueue.event_time() defining the time of
		events, t_complete sets the time up-to which the collection
		of events is known to be complete.  That is, in the future
		no new events will be pushed whose times are earlier than
		t_complete.
		"""
		self.queues[instrument].push(events, t_complete)

	def doublesgen(self, eventsa, offset_a, eventsb):
		"""
		For internal use only.
		"""
		# choose the longer of the two lists for the outer loop.
		# yes, it's counter-interuitive, it's because of the high
		# cost of indexing events for the inner loop.  in any case
		# we must still return coincident pairs with the events
		# ordered as supplied so if the event lists are swapped we
		# need to unswap the pairs that we return.
		# FIXME:  hopefully later improvements to the scaling will
		# change the logic here, and we'll need to switch to using
		# the shorter list for the outer loop after-all.

		if len(eventsb) > len(eventsa):
			eventsa, eventsb = eventsb, eventsa
			offset_a = -offset_a
			unswap = lambda a, b: (b, a)
		else:
			unswap = lambda a, b: (a, b)

		# for each event in list A, iterate over events from the
		# other list that are coincident with the event, and return
		# the pairs.  while doing this, collect the events that are
		# used in coincidences in a set for later logic.
		coinc_window = self.coinc_window
		used = set()
		used_add = used.add
		queueb_get_coincs = self.get_coincs(eventsb)
		for eventa in eventsa:
			matches = queueb_get_coincs(eventa, offset_a, coinc_window)
			if matches:
				eventa_id = id(eventa)
				used_add(eventa_id)
				for eventb in matches:
					eventb_id = id(eventb)
					used_add(eventb_id)
					yield unswap(eventa_id, eventb_id)
		self.used |= used

	def pull(self, t, flushed, flushed_unused):
		"""
		Generate a sequence of 2-element tuples of Python IDs of
		coincident events up to t, i.e., requiring the time of at
		least one of the events in each coincident pair to precede
		t.  The internal queues are flushed to t.  t cannot be
		greater than self.t_coinc_complete, or it can be the
		special value None (see below).

		The order of the IDs in each tuple is alphabetical by
		instrument.

		The Python IDs of all events flushed from the queues are
		placed in the flushed set.  The IDs of the flushed events
		that were never reported in a coincident pair, either now
		or in a previous call to .pull(), are also placed in the
		flushed_unused set.

		If t is None, then all events are flushed from the queues,
		and all remaining coincidences are constructed.  The object
		is reset to its initial state and can be re-used for a new
		event stream.

		NOTE:  the flushed and flushed_unused parameters passed to
		this function must be sets or set-like objects.  They must
		support Python set manipulation methods, like .clear(),
		.update(), and so on.
		"""
		# get the instrument names in alphabetical order, and
		# compute the time offset of the one wrt the other.

		instrumenta, instrumentb = sorted(self.offset_vector)
		offset_a = self.offset_vector[instrumenta] - self.offset_vector[instrumentb]

		# retrieve the event lists.  the .pull() methods will
		# safety check t, we don't have to do it here.

		flusheda, fornexttimea = self.queues[instrumenta].pull(t)
		flushedb, fornexttimeb = self.queues[instrumentb].pull(t)

		# construct coincident pairs.  NOTE: because we require any
		# coincidence that we report to contain an event preceding
		# t, the loop is done as a two-step process.  firstly, we
		# iterate over the "flushed" sequence for list A but search
		# for coincidences with any event from the flushed and "for
		# next time" sequences for list B.  lastly, we iterate over
		# the "for next time" sequence for list A but search for
		# coincidences only in the "flushed" sequence for list B.

		# FIXME:  when switching to Python 3 use
		#yield from self.doublesgen(flusheda, offset_a, flushedb + fornexttimeb)
		#yield from self.doublesgen(fornexttimea, offset_a, flushedb)
		for pair in self.doublesgen(flusheda, offset_a, flushedb + fornexttimeb):
			yield pair
		for pair in self.doublesgen(fornexttimea, offset_a, flushedb):
			yield pair

		# populate the flushed and flushed_unused sets

		flushed.clear()
		flushed.update(id(event) for event in flusheda)
		flushed.update(id(event) for event in flushedb)
		flushed_unused.clear()
		flushed_unused |= flushed - self.used

		# remove the events being flushed from the used set.  if
		# we've been flushed, there better not be anything left

		self.used -= flushed
		assert t is not None or not self.used


#
# =============================================================================
#
#                               Time Slide Graph
#
# =============================================================================
#


class TimeSlideGraphNode(object):
	def __init__(self, coincgen_doubles_type, offset_vector, coinc_window, min_instruments, time_slide_id = None):
		#
		# initialize
		#

		# time_slide_id non-None only in head nodes
		self.time_slide_id = time_slide_id
		self.offset_vector = offset_vector
		# not used here, but helpful to calling codes to avoid
		# repeating this test inside loops
		self.is_zero_lag = not any(offset_vector.values())
		# keep_partial is part of the logic that ensures we only
		# return coincs that meet the min_instruments criterion
		self.keep_partial = len(offset_vector) > min_instruments
		if len(offset_vector) > 2:
			self.components = tuple(TimeSlideGraphNode(coincgen_doubles_type, offset_vector, coinc_window, min_instruments) for offset_vector in offsetvector.component_offsetvectors([offset_vector], len(offset_vector) - 1))
			# view into the id() --> event indexes of the
			# nodes
			self.index = multidict(*(node.index for node in self.components))
		elif len(offset_vector) == 2:
			self.components = (coincgen_doubles_type(offset_vector, coinc_window),)
			self.index = self.components[0].index
		elif len(offset_vector) == 1:
			# use a singles queue directly, no coincidence
			self.components = (coincgen_doubles_type.singlesqueue(coinc_window),)
			self.index = self.components[0].index
		else:
			raise ValueError("offset_vector cannot be empty")

	@property
	def age(self):
		"""
		The earliest of the component nodes' .age.
		"""
		return min(node.age for node in self.components)

	@property
	def t_complete(self):
		"""
		The earliest of the component nodes' .t_complete.
		"""
		return min(node.t_complete for node in self.components)

	@property
	def t_coinc_complete(self):
		"""
		The earliest of the component nodes' .t_coinc_complete.
		"""
		return min(node.t_coinc_complete for node in self.components)

	def push(self, instrument, events, t_complete):
		"""
		Push new events from some instrument into the internal
		queues.  The iterable of events need not be time-ordered.
		With self.singlesqueue.event_time() defining the time of
		events, t_complete sets the time up-to which the collection
		of events is known to be complete.  That is, in the future
		no new events will be pushed whose times are earlier than
		t_complete.
		"""
		if len(self.offset_vector) == 1:
			# push directly into the singles queue
			assert (instrument,) == self.offset_vector.keys()
			self.components[0].push(events, t_complete)
		else:
			for node in self.components:
				if instrument in node.offset_vector:
					node.push(instrument, events, t_complete)

	def pull(self, t, verbose = False):
		"""
		Using the events contained in the internal queues upto t,
		construct and return all coincidences they participate in.
		Three objects are returned:  a tuple of the (n=N)-way
		coincidences, where N is the number of instruments in this
		node's offset vector, a tuple of the
		(min_instruments<=n<N)-way coincidences, which will be
		empty if N < min_instruments, and a set of the Python IDs
		of events that have been flushed from all internal queues
		by this operation.

		For the two tuples, the coincidences are reported as tuples
		of the Python id()'s of the events that form coincidences.
		The .index look-up table can be used to map these id()'s
		back to the original events.  To do that, however, a copy
		of the contents of the .index mapping must be made before
		calling this method because this method will result in some
		of the events in question being removed from the queues,
		and thus also from the .index mapping.
		"""
		#
		# NOTE:  it is essential that the same t be used for all
		# nodes when walking the graph.  the code assumes this to
		# be true to determine which events have been and have not
		# been used to construct coincidences:  we assume that
		# events that are being flushed from one instrument's queue
		# in one part of the graph must also be getting flushed
		# from all other queues for that same instrument elsewhere
		# in the graph, and so if they have not been used to form a
		# coincidence by now then they won't ever be.
		#
		# It is not possible to encode a safety check to ensure
		# this invariant, it is left to developers not to screw
		# with the graph traversal in that way.
		#

		#
		# is this a "no-op" node?  return 1-detector "coincs" and a
		# null set of "partial coincs".  NOTE:  the "no-op" node
		# concept exists to enable single-detector searches, i.e.,
		# searches where there is not, in fact, any coincidence
		# analysis to perform.  it's here so that calling codes
		# don't need to implement special cases themselves to
		# handle this.
		#

		if len(self.offset_vector) == 1:
			flushed = self.components[0].pull(t)[0]
			return tuple((eventid,) for eventid in flushed), set(), set(flushed)

		#
		# is this a leaf node?  construct the coincs explicitly
		#

		if len(self.offset_vector) == 2:
			if verbose:
				print("\tconstructing %s ..." % str(self.offset_vector), file=sys.stderr)

			#
			# search for and record coincidences.  coincs is a
			# sorted tuple of event ID pairs, where each pair
			# of IDs is, itself, ordered alphabetically by
			# instrument name.  here, for the two-instrument
			# case, the "partial coincs" are any singles that
			# failed to form coincidences, which corresponds to
			# the "flushed_unused" parameter of the
			# coincgen_doubles' .pull() method, but we convert
			# the set into a set of 1-element tuples so they
			# take the form of 1-detector coincs
			#

			flushed = set()
			partial_coincs = set()
			coincs = tuple(sorted(self.components[0].pull(t, flushed, partial_coincs)))
			return coincs, (set((eventid,) for eventid in partial_coincs) if self.keep_partial else set()), flushed

		#
		# this is a regular node in the graph.  use coincidence
		# synthesis algorithm to populate its coincs
		#

		#
		# double-check the consistency of our structure
		#

		assert len(self.components) > 2

		# first collect all coincs and partial coincs from the
		# component nodes in the graph
		component_coincs_and_partial_coincs_and_flushed = tuple(component.pull(t, verbose = verbose) for component in self.components)
		component_coincs = tuple(elem[0] for elem in component_coincs_and_partial_coincs_and_flushed)
		flushed = set.union(*(elem[2] for elem in component_coincs_and_partial_coincs_and_flushed))

		if self.keep_partial:
			# any coinc with n-1 instruments from the component
			# time slides might fail to form an n instrument
			# coincidence.  we put them all into our
			# partial_coincs set now, and we'll remove them
			# from this set as we discover which particiapte in
			# n instrument coincidences
			partial_coincs = set(itertools.chain(*component_coincs))

			# the (< n-1)-instrument partial coincs are more
			# complicated.  for example, H+L doubles are sought
			# out when forming H+L+V triples and also when
			# forming H+K+L triples, and if we are now trying
			# to form H+K+L+V quadruples we will have been
			# given two sets of "unused" H+L doubles, and there
			# is no reason for them to agree.  it can be shown
			# that if a (< n-1)-instrument coinc went unused in
			# forming any two of our (n-1)-instrument
			# components, then it cannot have been used to form
			# any of our (n-1)-instrument components.  we find
			# the unused (< n-1)-instrument partial coincs from
			# the union of all pair-wise intersections of the
			# (< n-1)-instrument partial coincs from our
			# components.
			for partial_coincsa, partial_coincsb in itertools.combinations((elem[1] for elem in component_coincs_and_partial_coincs_and_flushed), 2):
				partial_coincs |= partial_coincsa & partial_coincsb
		else:
			partial_coincs = set()
		del component_coincs_and_partial_coincs_and_flushed

		if verbose:
			print("\tassembling %s ..." % str(self.offset_vector), file=sys.stderr)
		# magic:  we can form all n-instrument coincs by knowing
		# just three sets of the (n-1)-instrument coincs no matter
		# what n is (n > 2).
		coincs = []
		component_coincs0 = component_coincs[0]
		component_coincs1 = component_coincs[1]
		component_coincs2 = component_coincs[-1]
		# for each coinc in list 0
		for coinc0 in component_coincs0:
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
			coincs1 = component_coincs1[bisect_left(component_coincs1, coinc0[:-1]):bisect_left(component_coincs1, coinc0[:-2] + (coinc0[-2] + 1,))]
			# find all the coincs in list 2 whose first (n-2)
			# event IDs are the same as the last (n-2) event
			# IDs in coinc0.  note that they are guaranteed to
			# be arranged together in the list and can be
			# identified with two bisection searches
			coincs2 = component_coincs2[bisect_left(component_coincs2, coinc0[1:]):bisect_left(component_coincs2, coinc0[1:-1] + (coinc0[-1] + 1,))]
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
					partial_coincs.difference_update(itertools.combinations(new_coinc, len(new_coinc) - 1))
					coincs.append(new_coinc)
		# sort the coincs we just constructed by the component
		# event IDs and convert to a tuple for speed
		coincs.sort()
		coincs = tuple(coincs)

		#
		# done.
		#

		return coincs, partial_coincs, flushed


class TimeSlideGraph(object):
	def __init__(self, coincgen_doubles_type, offset_vector_dict, coinc_window, min_instruments = 2, verbose = False):
		#
		# populate the graph head nodes.  these represent the
		# target offset vectors requested by the calling code.
		#

		if min_instruments < 1:
			raise ValueError("require min_instruments >= 1 (%d)" % min_instruments)
		if min(len(offset_vector) for offset_vector in offset_vector_dict.values()) < min_instruments:
			# this test is part of the logic that ensures we
			# will only extract coincs that meet the
			# min_instruments criterion
			raise ValueError("encountered offset vector (%s) smaller than min_instruments (%d)", (str(min(offset_vector_dict.values(), key = lambda offset_vector: len(offset_vector))), min_instruments))

		if verbose:
			print("constructing coincidence assembly graph for %d offset vectors ..." % len(offset_vector_dict), file=sys.stderr)
		self.head = tuple(
			TimeSlideGraphNode(
				coincgen_doubles_type, offset_vector, coinc_window, min_instruments, time_slide_id = time_slide_id
			) for time_slide_id, offset_vector in sorted(offset_vector_dict.items())
		)
		self.index = multidict(*(node.index for node in self.head))
		# the set of the Python id()'s of the events contained in
		# the internal queues that have formed coincident
		# candidates (including single-detector coincidences if
		# min_instruments = 1).  this is used to allow .pull() to
		# determine which of the events being flushed from the
		# internal queues has never been reported in a candidate.
		# calling codes can use this information to identify noise
		# samples for use in defining background models.
		self.used_ids = set()

		# the set of the Python id()'s of the events contained in
		# the internal queues that have been reported in coincident
		# candidates (including single-detector coincidences if
		# min_instruments = 1).  this is used to allow .pull() to
		# report which of the events it has found in coincidences
		# is being reported in a coincidence for the first time.
		# this can be used by calling code to be informed of which
		# events should be moved to an output document for storage.
		self.reported_ids = set()

		#
		# done
		#

		if verbose:
			def walk(node):
				# return the number of leaf and non-leaf
				# nodes in the graph rooted on this node
				return numpy.array((1, 0)) if len(node.components) == 1 else numpy.array((0, 1)) + sum(walk(node) for node in node.components)
			print("graph contains %d fundamental nodes, %d higher-order nodes" % tuple(sum(walk(node) for node in self.head)), file=sys.stderr)


	@property
	def age(self):
		"""
		The earliest of the graph's head nodes' .age.
		"""
		return min(node.age for node in self.head)


	@property
	def t_coinc_complete(self):
		"""
		The earliest of the graph's head nodes' .t_coinc_complete.
		"""
		return min(node.t_coinc_complete for node in self.head)


	def push(self, instrument, events, t_complete):
		"""
		Push new events from some instrument into the internal
		queues.  The iterable of events need not be time-ordered.
		With self.singlesqueue.event_time() defining the time of
		events, t_complete sets the time up-to which the collection
		of events is known to be complete.  That is, in the future
		no new events will be pushed whose times are earlier than
		t_complete.

		Returns True if the graph state has changed in a way that
		might allow new candidates to be constructed, False if it
		is not possible for new candidates to be constructed.
		"""
		t_before = self.t_coinc_complete
		for node in self.head:
			if instrument in node.offset_vector:
				node.push(instrument, events, t_complete)
		return self.t_coinc_complete != t_before


	def pull(self, newly_reported = None, flushed = None, flushed_unused = None, flush = False, coinc_sieve = None, event_collector = None, verbose = False):
		if verbose:
			print("constructing coincs for target offset vectors ...", file=sys.stderr)

		# NOTE:  we can only rely on the "flushed" and so on event
		# sets from the graph head nodes to be compatible with each
		# other if all nodes are .pull()ed to the same t.
		# otherwise what has been flushed from one node in the
		# graph might not yet be flushed from another.  we do use a
		# single time, but that means the latency is always the
		# worst-case scenario:  whichever graph node requires
		# triggers to be held the longest for the coincidence test
		# determines the latency of all graph nodes.  online we
		# only do a zero-lag analysis, and offline we don't care
		# about latency, so for now this is OK, but be aware of
		# this limitation if you are tempted to mess about with
		# this code's logic.
		t = None if flush else self.t_coinc_complete

		# flatten ID index for faster performance in loop.  NOTE:
		# this must be done before calling .pull() on the graph
		# because that operation will flush events from the
		# internal index, leaving it incomplete for our needs.  we
		# do it outside the loop over nodes because we'll need it
		# when the loop terminates
		index = dict(self.index.iteritems())

		# default coinc sieve
		if coinc_sieve is None:
			coinc_sieve = lambda events, offset_vector: False

		# default event_collector:
		if event_collector is None:
			event_collector_push = lambda event_ids, offset_vector: None
		else:
			event_collector_push = event_collector.push

		used_ids = set()
		newly_reported_ids = set()
		flushed_ids = set()
		for n, node in enumerate(self.head, start = 1):
			if verbose:
				print("%d/%d: %s" % (n, len(self.head), str(node.offset_vector)), file=sys.stderr)
			# we don't need to check that the coincs or partial
			# coincs contain at least min_instruments events
			# because those that don't meet the criteria are
			# excluded during coinc construction.
			coincs, partial_coincs, _flushed_ids = node.pull(t, verbose)
			flushed_ids |= _flushed_ids
			for coinc in itertools.chain(coincs, partial_coincs):
				used_ids.update(coinc)
				event_collector_push(coinc, node.offset_vector)
				# use the index to convert Python IDs back
				# to event objects
				events = tuple(index[event_id] for event_id in coinc)
				# apply the coinc sieve test.  this lever
				# allows calling code to reduce the event
				# rate going into output documents
				if not coinc_sieve(events, node.offset_vector):
					newly_reported_ids.update(coinc)
					yield node, events
		self.used_ids |= used_ids
		if newly_reported_ids:
			newly_reported_ids -= self.reported_ids
			self.reported_ids |= newly_reported_ids
		if flushed_ids:
			flushed_unused_ids = flushed_ids - self.used_ids
			self.used_ids -= flushed_ids
			self.reported_ids -= flushed_ids
		else:
			flushed_unused_ids = set()

		# if we've been flushed then there can't be any events left
		# in the queues (the or in parentheses is a boolean
		# operation, not a set operation)
		assert not flush or not (self.used_ids or self.reported_ids)

		# use the index to populate newly_used, flushed, and
		# flushed_unused with event objects
		if newly_reported is not None:
			newly_reported[:] = (index[event_id] for event_id in newly_reported_ids)
		if flushed is not None:
			flushed[:] = (index[event_id] for event_id in flushed_ids)
		if flushed_unused is not None:
			flushed_unused[:] = (index[event_id] for event_id in flushed_unused_ids)


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
	def __init__(self, xmldoc, coinc_definer_row):
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

		# look-up the coinc_def_id, creating a new one if required
		self.coinc_def_id = ligolw_coincs.get_coinc_def_id(xmldoc, coinc_definer_row.search, coinc_definer_row.search_coinc_type, create_new = True, description = coinc_definer_row.description)

		# find the time_slide table
		self.time_slide_table = lsctables.TimeSlideTable.get_table(xmldoc)
		self.time_slide_index = self.time_slide_table.as_dict()

	def coinc_rows(self, process_id, time_slide_id, events, table_name):
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
		coincmaps = [self.coincmaptable.RowType(
			coinc_event_id = None,
			table_name = table_name,
			event_id = event.event_id
		) for event in events]
		assert coincmaps, "coincs must contain >= 1 event"

		coinc = self.coinctable.RowType(
			process_id = process_id,
			coinc_def_id = self.coinc_def_id,
			coinc_event_id = None,
			time_slide_id = time_slide_id,
			insts = None,
			nevents = len(coincmaps),
			likelihood = None
		)

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

		self.rate_factors = {}
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
		# implementation.
		#
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
			if not instruments:
				# one-instrument case, no-op
				self.rate_factors[key] = 1.
			elif len(instruments) == 1:
				# two-instrument (1-D) case, don't use qhull
				self.rate_factors[key] = 2. * self.tau[frozenset((anchor, instruments[0]))]
			else:
				# three- and more instrument (2-D and
				# higher) case
				dimensions = len(instruments)	# anchor not included
				halfspaces = numpy.zeros((dimensions * (dimensions + 1), dimensions + 1), dtype = "double")
				# anchor constraints
				for i, instrument in enumerate(instruments):
					j = i
					i *= 2
					halfspaces[i, j] = +1.
					halfspaces[i + 1, j] = -1.
					halfspaces[i, -1] = halfspaces[i + 1, -1] = -self.tau[frozenset((anchor, instrument))]
				# non-anchor constraints
				for i, ((j1, a), (j2, b)) in enumerate(itertools.combinations(enumerate(instruments), 2), dimensions):
					i *= 2
					halfspaces[i, j1] = +1.
					halfspaces[i, j2] = -1.
					halfspaces[i + 1, j1] = -1.
					halfspaces[i + 1, j2] = +1.
					halfspaces[i, -1] = halfspaces[i + 1, -1] = -self.tau[frozenset((a, b))]
				# the origin is in the interior
				interior = numpy.zeros((len(instruments),), dtype = "double")
				# compute volume
				try:
					self.rate_factors[key] = spatial.ConvexHull(spatial.HalfspaceIntersection(halfspaces, interior).intersections).volume
				except AttributeError:
					# fall-through to old version
					pass
				else:
					# it worked, continue
					continue

				# old stone-throwing version in case qhull
				# is not available.  FIXME:  remove when we
				# are sure it's not needed.  for each
				# instrument 2...N, the interval within
				# which an event is coincident with
				# instrument 1
				windows = tuple((-self.tau[frozenset((anchor, instrument))], +self.tau[frozenset((anchor, instrument))]) for instrument in instruments)
				# pre-assemble a sequence of instrument
				# index pairs and the maximum allowed
				# \Delta t between them to avoid doing the
				# work associated with assembling the
				# sequence inside a loop
				ijseq = tuple((i, j, self.tau[frozenset((instruments[i], instruments[j]))]) for (i, j) in itertools.combinations(range(len(instruments)), 2))
				# compute the numerator and denominator of
				# the fraction of events coincident with
				# the anchor instrument that are also
				# mutually coincident.  this is done by
				# picking a vector of allowed \Delta ts and
				# testing them against the coincidence
				# windows.  the loop's exit criterion is
				# arrived at as follows.  after d trials,
				# the number of successful outcomes is a
				# binomially-distributed RV with variance =
				# d p (1 - p) <= d/4 where p is the
				# probability of a successful outcome.  we
				# quit when the ratio of the bound on the
				# standard deviation of the number of
				# successful outcomes (\sqrt{d/4}) to the
				# actual number of successful outcomes (n)
				# falls below rel accuracy: \sqrt{d/4} / n
				# < rel accuracy, or
				#
				# \sqrt{d} < 2 * rel accuracy * n
				#
				# note that if the true probability is 0,
				# so that n=0 identically, then the loop
				# will never terminate; from the nature of
				# the problem we know 0<p<1 so the loop
				# will, eventually, terminate.
				math_sqrt = math.sqrt
				random_uniform = random.uniform
				two_epsilon = 2. * abundance_rel_accuracy
				n, d = 0, 0
				while math_sqrt(d) >= two_epsilon * n:
					dt = tuple(random_uniform(*window) for window in windows)
					if all(abs(dt[i] - dt[j]) <= maxdt for i, j, maxdt in ijseq):
						n += 1
					d += 1
				self.rate_factors[key] = float(n) / float(d)
				for instrument in instruments:
					self.rate_factors[key] *= 2. * self.tau[frozenset((anchor, instrument))]

		# done computing rate_factors


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
		if self.tau and max(rates.values()) * max(self.tau.values()) >= 1.:
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
		M = self.R / (self.v * self.sigmas[:,numpy.newaxis])

		self.U, self.S, self.VT = numpy.linalg.svd(M)

		if len(rs) >= 3:
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
		>>> from numpy import testing
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
		array([[-0.45605637,  0.75800934,  0.46629865],
		       [-0.45605637,  0.75800934,  0.46629865]])
		>>> testing.assert_approx_equal(toa, 794546669.4269662)
		>>> testing.assert_approx_equal(chi2_per_dof, 0.47941941158371465)
		>>> testing.assert_approx_equal(dt, 0.005996370224459011)

		NOTE: if len(rs) >= 4, n is a 1x3 array.
		      if len(rs) == 3, n is a 2x3 array.
		      if len(rs) == 2, n is None.
		NOTE: n is not the source direction but the propagation
		      direction of GW.
		      Therefore, if you want source direction, you have to
		      multiply -1.
		NOTE: n is represented by earth fixed coordinate, not
		      celestial coordinate.
		      Up to your purpose, you should transform \\phi -> RA.
		      To do it, you can use dir2coord.
		"""
		assert len(ts) == len(self.sigmas)

		# change of t co-ordinate to avoid LIGOTimeGPS overflow
		t0 = min(ts)
		ts = numpy.array([float(t - t0) for t in ts])

		# sigma^-2 -weighted mean of arrival times
		tbar = sum(ts / self.sigmas**2) / sum(1 / self.sigmas**2)
		# the i-th element is ts - tbar for the i-th location
		tau = (ts - tbar) / self.sigmas

		tau_prime = numpy.dot(self.U.T, tau)[:3]

		if len(self.rs) >= 3:
			if self.singular:
				# len(rs) == 3
				np = numpy.array([tau_prime / self.S, tau_prime / self.S])
				try:
					np[0][2] =  math.sqrt(1.0 - np[0][0]**2 - np[0][1]**2)
					np[1][2] = -math.sqrt(1.0 - np[1][0]**2 - np[1][1]**2)
				except ValueError:
					# two point is mergered, n_0 = n_1
					np.T[2] = 0.0
					np /= math.sqrt(numpy.dot(np[0], np[0]))

				# compute n from n'
				n = numpy.array([numpy.zeros(3), numpy.zeros(3)])
				n = numpy.dot(self.VT.T, np.T).T

				# safety check the nomalization of the result
				assert abs(numpy.dot(n[0], n[0]) - 1.0) < 1e-8
				assert abs(numpy.dot(n[1], n[1]) - 1.0) < 1e-8

				# arrival time at origin
				toa = sum((ts - numpy.dot(self.rs, n[0]) / self.v) / self.sigmas**2) / sum(1 / self.sigmas**2)

				# chi^{2}
				chi2 = sum((numpy.dot(self.R, n[0]) / (self.v * self.sigmas) - tau)**2)

				# root-sum-square timing residual
				dt = ts - toa - numpy.dot(self.rs, n[0]) / self.v
				dt = math.sqrt(numpy.dot(dt, dt))
			else:
				# len(rs) >= 4
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
				chi2 = sum((numpy.dot(self.R, n) / (self.v * self.sigmas) - tau)**2)

				# root-sum-square timing residual
				dt = ts - toa - numpy.dot(self.rs, n) / self.v
				dt = math.sqrt(numpy.dot(dt, dt))
		else:
			# len(rs) == 2
			# s_1 == s_2 == 0

			# compute an n'
			xp = tau_prime[0] / self.S[0]
			try:
				np = numpy.array([xp, 0, math.sqrt(1-xp**2)])
			except ValueError:
				# two point is mergered, n_0 = n_1
				np = numpy.array([xp, 0, 0])

			# compute n from n'
			n = numpy.dot(self.VT.T, np)

			# safety check the nomalization of the result
			assert abs(numpy.dot(n, n) - 1.0) < 1e-8

			# arrival time at origin
			toa = sum((ts - numpy.dot(self.rs, n) / self.v) / self.sigmas**2) / sum(1 / self.sigmas**2)

			# chi^{2}
			chi2 = sum((numpy.dot(self.R, n) / (self.v * self.sigmas) - tau)**2)

			# root-sum-square timing residual
			dt = ts - toa - numpy.dot(self.rs, n) / self.v
			dt = math.sqrt(numpy.dot(dt, dt))

			# set n None
			n = None

		# done
		return n, t0 + toa, chi2 / len(self.sigmas), dt

	@staticmethod
	def dir2coord(n, gps):
		"""
		This transforms from propagation direction vector to right
		ascension and declination source co-ordinates.

		The input is the propagation vector, n, in Earth fixed
		co-ordinates (x axis through Greenwich merdian and equator,
		z axis through North Pole), and the time.  n does not need
		to be a unit vector.  The return value is the (ra, dec)
		pair, in radians.  NOTE:  right ascension is not
		necessarily in [0, 2\\pi).
		"""
		# safety checks
		if len(n) != 3:
			raise ValueError("n must be a 3-vector")

		# normalize n
		n /= math.sqrt(numpy.dot(n, n))

		# transform from propagation to source direction
		n = -n

		# compute ra, dec
		RA = numpy.arctan2(n[1], n[0]) + lal.GreenwichMeanSiderealTime(gps)
		DEC = numpy.arcsin(n[2])

		# done
		return RA, DEC


#
# =============================================================================
#
#                     Coincidence Parameter Distributions
#
# =============================================================================
#


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

		NOTE:  it is possible for sub-classes to override this
		method, and chain to it if they wish.  There is no
		requirement that this method evaluate the ratio .numerator
		/ .denominator, for example it would not invalidate the
		output of .ln_lr_samples() if the computation of the return
		value includes some kind of cuts, or other non-trivial
		logic.
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

	def ln_lr_samples(self, random_params_seq, signal_noise_pdfs = None):
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
		first element of each tuple to be added to histograms using
		the two relative frequencies as weights, i.e., the two
		relative frequencies give the number of times one should
		consider this one draw of log likelihood ratio to have
		occured in the two populations.

		random_params_seq is a sequence (generator is OK) yielding
		3-element tuples whose first two elements provide the *args
		and **kwargs values passed to the numerator and denominator
		density functions, and whose thrid element is the natural
		logarithm of the probability density from which the
		parameters have been drawn evaluated at the parameters.

		On each iteration, the *args and **kwargs values yielded by
		random_params_seq are passed to our own .__call__() method
		to evalute the log likelihood ratio at that choice of
		parameter values.  If signal_noise_pdfs is None the
		parameters are also passed to the .__call__() mehods of our
		own .numerator and .denominator attributes to obtain the
		signal and noise population densities at those parameters.
		If signal_noise_pdfs is not None then, instead, the
		parameters are passed to the .__call__() methods of its
		.numerator and .denominator attributes to obtain those
		densities.

		If histograming the results as described above, the effect
		is to draw paramter values from the signal and noise
		populations defined by signal_noise_pdfs' PDFs but with log
		likelihood ratios evaluated using our own PDFs.
		"""
		if signal_noise_pdfs is None:
			lnP_signal_func = self.numerator
			lnP_noise_func = self.denominator
		else:
			lnP_signal_func = signal_noise_pdfs.numerator
			lnP_noise_func = signal_noise_pdfs.denominator
		for args, kwargs, lnP_params in random_params_seq:
			lnP_signal = lnP_signal_func(*args, **kwargs)
			lnP_noise = lnP_noise_func(*args, **kwargs)
			yield self(*args, **kwargs), lnP_signal - lnP_params, lnP_noise - lnP_params
