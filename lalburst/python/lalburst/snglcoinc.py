# Copyright (C) 2006--2021  Kipp Cannon, Drew G. Keppel, Jolien Creighton
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
Generic time interval coincidence engine.
"""


from bisect import bisect_left
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
from collections import ChainMap, Counter
import warnings


from ligo.lw import ligolw
from ligo.lw import lsctables
from ligo.lw.utils import coincs as ligolw_coincs
from ligo import segments
import lal
from . import offsetvector


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
from .git_version import date as __date__
from .git_version import version as __version__


#
# =============================================================================
#
#                                  Utilities
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


#
# =============================================================================
#
#                    Streaming Double Coincidence Generator
#
# =============================================================================
#


#
# Algorithm notes:
#
# The ingredients are shown in the diagram below:
#
#      #1 -----X----------X--X------------------X---X------)----->
#
#           #2 ------X--------------X------X-X------------X)---X-X---->
#
# #3 ----X------X---------x--------------x)----------------->
# ----------------------------------------------------------------------> t
#
#
# "X" = an event.  Each event has a unique discrete time associated with
# it.  Physically the event might represent a signal lasting many tens
# of minutes, but regardless of that we require a "time" to be assigned to
# that signal, and that time will form the basis of the coincidence test.
#
# "#1", etc. = a detector that provides a sequence of events.  Each
# detector has a time offset applied to its events
#
# ")" = t_complete for that detector.  Calling code guarantees that the
# event sequence for that detector is complete up to (not necessarily
# including) this time.  Events might already have been collected at or
# following this time, but no further events will be collected preceding
# that time.  In the diagram above the event lists are shown displaced
# according to their respective offsets, and the t_completes shifted
# accordingly.
#
# For each pair of detectors there is a coincidence window which is assumed
# to equal some fixed time interval common to all detector pairs + the
# light travel time between that pair.  Two events are coincident if their
# times differ by no more than the coincidence window for their respective
# detector pair.  three or more events are coincident if all pairs of
# events within the n-tuple are coincident.
#
# Looking at the events again:
#
#      #1 -----X----------X--X------------------X---X------)----->
#
#           #2 ------X--------------X------X-X------------X)---X-X---->
#
# #3 ----X------X---------x--------------x)----------------->
# ----------------------------------------------------------------------> t
#                                         ^                ^
#                                         B                A
#
# The event list for detector #2 is complete up to
#
# t = t_complete(#2) + offset(#2) = A,
#
# and if we were only concerned about single-detector candidates from
# detector #2 then we could scan up to that time and report everything that
# we find and know that we had found all such candidates without omissions
# or double-counting any.  However, that last event in detector #2, for
# example, just before that detector's t_complete, might later be found to
# be coincident with something from detector #3, whose event list at this
# time is still only complete up to
#
# t = t_complete(#3) + offset(#3) = B
#
# Likewise, #1+#2 doubles can only be constructed up to approximately t =
# B because they might actually be found to form triples once the event
# list for detector #3 is filled in.
#
# We need to determine which candidates or potential candidates we have
# sufficiently complete information about to determine if they form a
# coincident n-tuple;  and, after that, which events can be disposed of
# because all possible n-tuple candidates that can be constructed from them
# have been constructed (or, equivalently, which events to retain because
# we have not yet been able to decide all the candidates they might
# participate in);  and, additionally, we need to do this in a way that
# upon repeated applications of this process no candidate or part thereof
# is reported more than once and no candidate that can be constructed from
# the (complete) event sequences is missed.
#
# The algorithm we use for achieving this is the following.  We need the
# following ingredients:
#
# - an exactly reproducible method by which to assign a "time" to an
# n-tuple candidate so that two different codes (possibly using different
# arithmetic optimiziations, or the same code doing the calculation a
# second time) will agree whether a given candidate is before or after (or
# at) some given time boundary;
#
# - given that definition, a method by which to determine the latest such
# time up to which all possible candidates can be constructed;
#
# - a method by which to determine which events can no longer participate
# in candidates on subsequent invocations of the algorithm.
#
# In particular, we need these to act together so that the decision as to
# whether or not a candidate is before or after the time bound we choose is
# stable against the addition of new events to the queues (the bound must
# precede the times up to which event lists are complete by a sufficent
# margin), and we will need to retain events long enough to make this
# definition stable against removal of events from the queues at least
# until the risk of mistaking subsets of the n-tuple as new candidates has
# passed.  It is desirable for the algorithm to achieve the lowest possible
# latency, meaning always reporting as many candidates as can be decided
# given the information available, but performance is the primary
# consideration so we also require an algorithm that can make decisions
# about candidates very quickly without using much memory.  For example,
# comparing every potential candidate against a list of everything that has
# previously been reported is impractical.
#
# We will define the "time" of an n-tuple candidate to be the earliest time
# of the events in the candidate with their respective offsets applied.
# For n-tuple candidates whoses "times" precede the time up to which all
# detectors' event lists are complete this is stable against addition of
# new events to the detector queues, even for those candidates that are
# modified by the addition of extra coincident events as the event
# sequences are advanced.
#
# Looking at our event sequences from above again, we now identify 6 times
# of interest:
#
#      #1 -----O---|----|-X--X--|----|----|-|---X---X------)----->
#                  |    |       |    |    | |
#           #2 ----|-X--|-------|---X|----|X|X------------X)---X-X---->
#                  |    |       |    |    | |
# #3 ----O------O--|----|-x-----|----|---x)-|--------------->
# -----------------|----[=======|====)----|-|---------------------------> t
#                  A    B       C    D    E F
#
# E is the earliest of the time-shifted t_compeletes, in this case it
# corresponds to #3's time-shifted t_complete.  D is the closest
# time-shifted time an event in another detector can have to E before we
# cannot say if it forms a coincidence with events in #3 that we have not
# yet collected.  Since the #3's event list is complete upto but not
# including E, we can collect all n-tuples whose time (time of earliest
# candidate) is upto but not including D.  B was the boundary upto
# which we collected candidates in the previous iteration.  Since at that
# time we did not allow candidates whose time was exactly B to be
# collected, we now include them.  Therefore the range of candidate times
# to be collected is [B, D) (closed from below, open from above).  In order
# to know that an n-tuple's earliest coincident member is not prior to B
# (and therefore it would have already been reported as part of a larger
# coincidence involving that event) we must retain and include for
# consideration events as early as A.  A precedes B by the largest
# coincidence window between any pair of detectors.  In this iteration, any
# events that had preceded A have been disposed of (and are indicated now
# as an "O").  After this iteration is complete, we will have reported
# exactly all candidates who canonical times fall in [B, D) and they will
# not be reported again, but for the same reason we retained events between
# A and B we must not discard all of them yet.  We must retain all the
# events whose time-shifted times are C or later where C precedes D by the
# largest coincidence window between any pair of detectors.  Since we are
# only forming coincidences involving at least one event preceding D, we do
# not need to include in the coincidence analysis at this iteration any of
# the events whose time-sifted times are F or later, where F antecedes
# D by the largest coincidence window between any pair of detectors.
#
# The procedure is to (i) compute these boundary times, (ii) collect all
# the events from all detectors upto but not including F, (iii) perform a
# full coincidence analysis, (iv) report as newly-identified any candidates
# whose canonical times fall in the interval [B, D), ignoring the rest, (v)
# then discard all events up to C, (vi) record D to be used as B next time,
# (vii) and return.
#
# There is an additional complication.  One of the problems we need to
# solve is to provide the set of single-detector events that participate in
# the n-tuple candidates, across all time-shift offset vectors.  For
# example, we need a list of the detector #1 events that participate in any
# n-tuple candidate, but we need each one to be listed only once,
# regardless of how many different offset vectors were found to result in
# it participating in an n-tuple candidate.
#
# We solve this problem the following way.  We maintain a set of the events
# that have participated in n-tuple coincidences.  The first time an event
# is found to participate in an n-tuple candidate it is added to the set
# and reported as a newly used single-detector event to be inserted into
# the calling code's output document.  Events that participate in n-tuple
# coincidences but are found to already be in that set are not reported to
# the calling code (the n-tuple coincidence is, but not the single-detector
# event).  To prevent the set from growing without bound, causing both memory
# and performance problems, events are removed from the set when they are
# removed from the last queue that contains them.
#
# How to compute the times.  For each time slide,
#
# - D is the minimum across all instruments of
#
# D = min(t_complete + offset - max(coincidence window with other detectors))
#
# Note that this is slightly different than the algorithm described above,
# where we explained that D was computed relative to E, the minimum of the
# time-shifted completes.  In fact we compute the boundary uniquely for
# each detector and take the minimum of those.  This is simpler, but also
# is necessary in the event that two detectors' t_completes are very close
# and the later one has a larger maximum coincidence window than the
# earlier one, larger by an amount greater than the difference between the
# t_completes.
#
# - B is this value computed in the previous iteration.
#
# - F is computed as a unique value for each detector as
#
# F = D + max(coincidence window with other detectors)
#
# For simplicity, the description above implied F was common to all
# detectors, but we can reduce the computational cost by trimming it
# as tightly as we can for each detector.
#
# - likewise, A is computed as a unique value for each detector as
#
# A = B - max(coincidence window with other detectors)
#
# - and C, also, is computed as a unique value for each detector as
#
# C = D - max(coincidence window with other detectors)
#
# For each detector, all events in the interval [A, F) for that detector
# are collected and fed into the coincidence constructor.  When finished,
# all events upto but not including C are discarded and reported to the
# top-level of the coincidence engine as being no-longer in storage, which
# is information it uses to manage the set of events that have participated
# in n-tuple coincidences.
#
#
# Some pipelines require a list of the events known to not form
# coincidences.  This information might be used, for example, to assist in
# construting noise models for false-alarm probability estimation.
# Therefore, in addition to constructing n-tuple coincs we also construct a
# stream of non-coincident single-detector events (which might also be
# reported as n-tuple candidates if n is 1).  This is achieved by
# maintainingg a set of all events that participate in (n>=2)-tuple coincs,
# and then at the point where events are flushed from the queues any that
# are not found in that set are reported in the "non-coincident event"
# stream.  Any that are found in the set as they are being flushed are at
# that time removed from the set to keep it from growing.
#
# Implementation Notes:
#
# The snglesqueue, coincgen_doubles, and the TimeSlideGraphNode object
# export APIs that are partially compatible with each other.  This is done
# so that they can be used interchangably within the time slide graph to
# implement coincidence generators of different orders without the need to
# wrap them in conditionals to select the correct code path.  This is why,
# for example, the concept of the coincidence window and the distinction
# between t_complete and t_coinc_complete have been pushed down into the
# snglesqueue object even though they seem to be nonsensical there.
#


class singlesqueue(object):
	"""
	Queue a stream of partially time ordered events:  as new events are
	added, sort them into the queue in time order;  maintain a record
	of the time up-to which the list of events is complete;  and
	provide events sorted in time order to coincidence tests.

	Searches must define a class subclassed from this that defines an
	override of the .event_time() method.  See coincgen_doubles for
	information on what to do with the subclass.

	Implementation notes:

	What is actually inserted into the queue is the time of the event
	as defined by the .event_time() method provided by the subclass.
	This allows sorts and bisection searches to be performed quickly,
	without the need for much additional code.  When we need to
	retrieve an event, we need some way of turning the time of
	the event (i.e., the thing that is actually in the queue) back into
	that event, and this is done by adding an attribute to the time
	object that points back to the original event.  This means that
	builtin Python numeric types cannot be used for the times of the
	events because new attributes cannot be attached to them.
	Subclasses of the builtin numeric types can be used as time
	objects, as can be LAL's LIGOTimeGPS.
	"""
	@staticmethod
	def event_time(event):
		"""
		Override with a method to return the "time" of the given
		event.  The "time" object that is returned is typically a
		lal.LIGOTimeGPS object.  If some other type is used, it
		must support arithmetic and comparison with python float
		objects and lal.LIGOTimeGPS objects, and support comparison
		with ligo.segments.PosInfinity and
		ligo.segments.NegInfinity, and it must have a .__dict__ or
		other mechanism allowing an additional attribute named
		.event to be set on the object.  The builtin float type is
		not suitable, but a subclass of float would be.
		"""
		raise NotImplementedError

	def queueentry_from_event(self, event):
		"""
		For internal use.  Construct the object that is to be
		inserted into the queue from an event.
		"""
		# retrieve the time and attach a reference back to the
		# original event object
		entry = self.event_time(event)
		entry.event = event
		return entry

	def __init__(self, offset, max_coinc_window):
		"""
		Initialize a new, empty, singles queue.  offset is the time
		that will be added to events in this queue when forming
		coincidences with events from other queues, and
		max_coinc_window is an upper bound on the maximum time that
		can separate events in this queue from events in other
		queues and they still be able to form coincidences --- it
		is not that maximum time, necessarily, but the true maximum
		time cannot be larger than max_coinc_window.

		Note that max_coinc_window is not defining the coincidence
		test, that is defined by the get_coincs machinery within
		the coincgen_doubles class.  max_coinc_window is used by
		the streaming logic to determine where in the event
		sequences complete n-tuple candidates can be constructed.
		Although max_coinc_window need only bound the time interval
		described above, the latency of the coincidence engine is
		reduced and the number of events stored in the queue (and
		the computational overhead they cause) is reduced if the
		number is as small as possible.
		"""
		self.offset = offset
		if max_coinc_window < 0.:
			raise ValueError("max_coinc_window < 0 (%g)" % max_coinc_window)
		self.max_coinc_window = max_coinc_window
		# using .event_time() to define the times of events, the
		# list of events is complete upto but not including
		# .t_complete.  we can't use fpconst.NegInf because
		# LIGOTimeGPS objects refuse to be compared to it.  the
		# NegInfinity object from the segments library, however, is
		# compatible with both LIGOTimeGPS and with native Python
		# numeric types.
		self.t_complete = segments.NegInfinity
		# queue of events.  the queue's contents are time-ordered.
		self.queue = []
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
		queue is empty.  The time scale for .age is that of the
		events themselves, not their shifted times (the offset
		parameter is not applied).

		This is not used by the coincidence engine.  This
		information is provided as a convenience for calling code.
		For example, some performance improvements might be
		realized if segment lists or other information are clipped
		to the range of times spanned by the contents of the
		coincidence engine queues, and this allows calling codes to
		see what interval that is, in physical event times not
		time-shifted event times.
		"""
		return self.queue[0] if self.queue else self.t_complete

	@property
	def t_coinc_complete(self):
		"""
		The time up to which the time-shifted events in other
		detectors' queues can be considered complete for the
		purpose of forming n-tuple coincidences with the
		time-shifted events in this queue.  This precedes the
		(time-shifted) time up to which this event list is complete
		by the maximum coincidence window that can separate an
		event in this queue from events in others and they still be
		considered coincident in time.

		The earliest such time across all detectors is the time up
		to which the n-tuple candidate list can be completely
		constructed, assuming the "time" of an n-tuple candidate is
		defined to be the earliest time-shifted time of the events
		in the n-tuple.
		"""
		return self.t_complete + self.offset - self.max_coinc_window

	def push(self, events, t_complete):
		"""
		Add new events to the queue.  Mark the queue complete up to
		t_complete, meaning you promise no further events will be
		added to the queue earlier than t_complete.  The time scale
		for t_complete is that of the events themselves, not their
		shifted times (the offset parameter is not applied).  The
		value of t_complete is stored for later retrieval.

		NOTE:  t_complete is not permitted to decrease.

		NOTE:  the events sequence will be iterated over multiple
		times.  It may not be a generator.
		"""
		if t_complete < self.t_complete:
			raise ValueError("t_complete has gone backwards:  last was %s, new is %s" % (self.t_complete, t_complete))

		if events:
			# add the new events to the ID index
			self.index.update((id(event), event) for event in events)

			# construct queue entries from the new events and
			# insort with current "incomplete" queue
			self.queue.extend(map(self.queueentry_from_event, events))
			self.queue.sort()

		# update the marker labelling time up to which the event
		# list is complete
		self.t_complete = t_complete

	def pull(self, t):
		"""
		t is the time up to which the calling code is in the
		process of constructing all available n-tuple coincident
		candidates.  The return value is a tuple of all events from
		the queue whose time-shifted end times are up to (not
		including) t + max_coinc_window, which are all the events
		from this queue that could form a coincidence with an event
		up to (not including) t.

		t is defined with respect to the time-shifted event times.
		The contents of the returned tuple is time-ordered as
		defined by .event_time().

		Assuming that t cannot go backwards, any of the reported
		events that cannot be used again are removed from the
		internal queue.

		If t is None, then the queue is completely flushed.  All
		remaining events are pulled from the queue and reported in
		the tuple, .t_complete is reset to -inf and all other
		internal state is reset.  After calling .pull() with t =
		None, the object's state is equivalent to its initial state
		and it is ready to process a new stream of events.
		"""
		# the input parameter t is the time D in the algorithm
		# description above

		# are we being flushed?  if so, return everything we've got
		# and reset out state to .__init__() equivalent
		if t is None:
			events = tuple(entry.event for entry in self.queue)
			del self.queue[:]
			self.index.clear()
			self.t_complete = segments.NegInfinity
			return events

		# unshift back to our events' time scale
		t -= self.offset

		# safety check
		if t + self.max_coinc_window > self.t_complete:
			raise ValueError("pull to %g fails to precede time to which queue is complete, (%g + %g), by the required margin of %g s" % (t + self.offset, self.t_complete, self.offset, self.max_coinc_window))

		# these events will never be used again.  remove them from
		# the queue, and remove their IDs from the index.  this is
		# calculating time C for this event list in the description
		# above.
		i = bisect_left(self.queue, t - self.max_coinc_window)
		flushed_events = tuple(entry.event for entry in self.queue[:i])
		del self.queue[:i]
		for event in flushed_events:
			self.index.pop(id(event))
		# collect other events that might be used again but that
		# can form coincidences with things in other detectors up
		# to t.  this is computing time F for this event list in
		# the description above.
		return flushed_events + tuple(entry.event for entry in self.queue[:bisect_left(self.queue, t + self.max_coinc_window)])


class coincgen_doubles(object):
	"""
	Using a pair of singlesqueue objects, constructs pairs of
	coincident events from streams of partially time-ordered events.

	Searches must subclass this.  The .singlesqueue class attribute
	must be set to the search-specific singlesqueue implementation to
	be used, i.e., a subclass of singlesqueue with an appropriate
	.event_time() override.  The internally-defined .get_coincs class
	must be overridden with an implementation that provides the
	required .__init__() and .__call__() methods.
	"""
	# subclass must override this attribute with a reference to an
	# implementation of singlesqueue that implements the .event_time()
	# method.
	singlesqueue = singlesqueue

	class get_coincs(object):
		"""
		This class defines the coincidence test.  An instance is
		initialized with a sequence of events, and is a callable
		object.  When the instance is called with a single event
		from some other instrument, a time offset to apply to that
		event (relative to the events in this list) and a time
		coincidence window, the return value must be a (possibly
		empty) sequence of the initial events that are coincident
		with that given event.

		The sequence of events with which the instance is
		initialized is passed in time order.

		It is not required that the implementation be subclassed
		from this.  This placeholder implementation is merely
		provided to document the required interface.

		A minimal example:

		class get_coincs(object):
			def __init__(self, events):
				self.events = events
			def __call__(self, event_a, offset_a, coinc_window):
				return [event_b for event_b in self.events if abs(event_a.time + offset_a - event_b.time) < coinc_window]

		This is performance-critical code and a naive
		implementation such as the one above will likely be found
		to be inadequate.  Expect to implement this code in C for
		best results.
		"""
		# FIXME:  move offset_a and coinc_window parameters to
		# __init__() method
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
			separating coincident events from the shifted time of
			event_a.  Here, this is the interval that defines
			the coincidence test between the two detectors, not
			the bound on all such intervals used by the
			singlesqueue object to determine the interval for
			which n-tuple candidates can be constructed.
			"""
			raise NotImplementedError

	def __init__(self, offset_vector, coinc_windows):
		"""
		offset_vector must be a two-instrument offset vector.  This
		sets which two instruments' events are to be processed by
		this object, and the offsets that should be applied to
		their events when searching for coincidences.
		coinc_windows is a dictionary mapping pairs of instruments
		(as sets) to the largest time that can separate events from
		those pairs of instruments (after applying the required
		time offsets) and they still be considered coincident.
		"""
		if len(offset_vector) != 2:
			raise ValueError("offset_vector must name exactly 2 instruments (got %d)" % len(offset_vector))
		self.offset_vector = offset_vector
		# the coincidence window for our pair of detectors
		self.coinc_window = coinc_windows[frozenset(offset_vector)]
		if self.coinc_window < 0.:
			raise ValueError("coinc_window < 0 (%g)" % coinc_window)
		# the largest coincidence window between any pair of
		# detectors in the network
		max_coinc_window = max(coinc_windows.values())
		if max_coinc_window < 0.:
			raise ValueError("max(coinc_window) < 0 (%g)" % max_coinc_window)
		# initialize the singlesqueue objects, one for each of our
		# detectors
		# FIXME: the max_coinc_window is the largest of all
		# coincidence windows, but as described above in the notes
		# about this algorithm we only need it to be the largest of
		# the windows that include the specific instrument.  fixing
		# this will reduce latency slightly (tens of milliseconds)
		# and keep the in-ram event lists smaller and reduce
		# coincidence overhead
		self.queues = dict((instrument, self.singlesqueue(offset, max_coinc_window)) for instrument, offset in offset_vector.items())
		# view into the id() --> event indexes of the queues
		self.index = ChainMap(*(queue.index for queue in self.queues.values()))

		# pre-sort the instruments into alphabetical order
		self.instrumenta, self.instrumentb = sorted(self.queues)
		# pre-compute the offset of the 1st wrt to the 2nd
		self.offset_a = self.queues[self.instrumenta].offset - self.queues[self.instrumentb].offset

	@property
	def age(self):
		"""
		The earliest of the internal queues' .age.
		"""
		return min(queue.age for queue in self.queues.values())

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
		t_complete.  t_complete is not permitted to decrease.
		"""
		self.queues[instrument].push(events, t_complete)

	def doublesgen(self, eventsa, offset_a, eventsb, singles_ids):
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

		# initialize the singles_ids set

		singles_ids.clear()
		singles_ids.update(id(event) for event in eventsb)
		singles_ids_add = singles_ids.add
		singles_ids_discard = singles_ids.discard

		# for each event in list A, iterate over events from the
		# other list that are coincident with the event, and return
		# the pairs.  while doing this, collect the IDs of events
		# that are used in coincidences in a set for later logic.
		coinc_window = self.coinc_window
		queueb_get_coincs = self.get_coincs(eventsb)
		for eventa in eventsa:
			matches = queueb_get_coincs(eventa, offset_a, coinc_window)
			eventa = id(eventa)
			if matches:
				for eventb in matches:
					eventb = id(eventb)
					singles_ids_discard(eventb)
					yield unswap(eventa, eventb)
			else:
				singles_ids_add(eventa)

	def pull(self, t, singles_ids):
		"""
		Generate a sequence of 2-element tuples of Python IDs of
		coincident event pairs.  Calling code is assumed to be in
		the process of constructing all possible n-tuple
		coincidences such that at least one time-shifted event in
		the n-tuple is before t, and so we construct and return all
		pairs required to decide that set, meaning we will
		construct and return pairs that might contain only events
		after t so that the calling code can test them for
		coincidence with events from other instruments that precede
		t.  The internal queues are flushed of events that cannot
		be used again assuming that t never goes backwards.  If t
		is the special value None then all coincident pairs that
		can be constructed from the internal queues are constructed
		and returned, the queues are emptied and returned to their
		initial state.

		The order of the IDs in each tuple is alphabetical by
		instrument.

		The Python IDs of events that are not reported in a
		coincident pair are placed in the singles_ids set.

		NOTE:  the singles_ids parameter passed to this function
		must a set or set-like object.  It must support Python set
		manipulation methods, like .clear(), .update(), and so on.
		"""
		# retrieve the event lists to be considered for coincidence
		# using .pull() on the singlesqueues.  these return
		# sequences of event objects, not their python IDs.  the
		# .pull() methods will safety check t against their
		# time-shifted t_completes, we don't have to do it here.
		# from the event sequences, we construct coincident pairs.
		# the pairs are the python IDs of the coincident events,
		# not the events themselves.  while doing this, the
		# singles_ids set is populated with python IDs of events
		# that do not form pairs.

		return sorted(self.doublesgen(self.queues[self.instrumenta].pull(t), self.offset_a, self.queues[self.instrumentb].pull(t), singles_ids))


#
# =============================================================================
#
#                               Time Slide Graph
#
# =============================================================================
#


class TimeSlideGraphNode(object):
	def __init__(self, coincgen_doubles_type, offset_vector, coinc_windows, min_instruments, time_slide_id = None):
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
		# construct this node's children
		if len(offset_vector) > 2:
			self.components = tuple(TimeSlideGraphNode(coincgen_doubles_type, component_offset_vector, coinc_windows, min_instruments) for component_offset_vector in self.component_offsetvectors(offset_vector, len(offset_vector) - 1))
			# view into the id() --> event indexes of the
			# nodes.  we assume they are all ChainMap
			# instances, and, for performance reasons, flatten
			# the hierarchy one level.  if that was also done
			# at the lower levels, ensures that, globally, the
			# ChainMap hierarchy doesn't grow in depth as the
			# graph is constructed.  it's only ever a ChainMap
			# of dicts, not a ChainMap of ChainMaps of ...
			self.index = ChainMap(*(index for node in self.components for index in node.index.maps))
		elif len(offset_vector) == 2:
			self.components = (coincgen_doubles_type(offset_vector, coinc_windows),)
			self.index = self.components[0].index
		elif len(offset_vector) == 1:
			# use a singles queue directly, no coincidence.
			# note that this configuration is not how singles
			# are generated in a normal analysis.  this is a
			# special case to support analyses processing only
			# a single detector (possibly being run for
			# diagnostic or performance testing purposes).
			# normal multi-detector analyses get singles from
			# the unused events returned by the doubles
			# generator in the len==2 code path.
			assert not coinc_windows
			self.components = (coincgen_doubles_type.singlesqueue(0.),)
			self.index = self.components[0].index
		else:
			raise ValueError("offset_vector cannot be empty")
		# the time up to which we have been pulled already
		self.previous_t = segments.NegInfinity
		# local reference to event_time() function.  required for
		# cutting candidate sets to time boundaries.
		self.event_time = coincgen_doubles_type.singlesqueue.event_time

	@staticmethod
	def component_offsetvectors(offset_vector, n):
		"""
		Equivalent to offsetvector.component_offsetvectors() except
		only one input offsetvector is considered, not a sequence
		of them, and the output offsetvector objects preserve the
		absolute offsets of each instrument, not only the relative
		offsets between instruments.
		"""
		assert 1 <= n <= len(offset_vector)
		for selected_instruments in itertools.combinations(offset_vector, n):
			yield offsetvector.offsetvector((instrument, offset_vector[instrument]) for instrument in selected_instruments)

	@property
	def age(self):
		"""
		The earliest of the component nodes' .age.
		"""
		return min(node.age for node in self.components)

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

	def pull(self, t):
		"""
		Using the events contained in the internal queues construct
		and return all coincidences they participate in.  It is
		assumed calling code is in the process of constructing
		n-tuple candidates such that each candidate contains at
		least one event preceding t, and so we construct and return
		all candidates required to decide that set, which might
		include candidates all of whose events come after t.

		Two objects are returned:  a tuple of the (n=N)-way
		coincidences, where N is the number of instruments in this
		node's offset vector, and a set of the
		(min_instruments<=n<N)-way coincidences, which will be
		empty if N < min_instruments.

		Since the (min_instruments<=n<N)-tuple coincidences cannot
		go on to form (n>N)-tuple coincidences, any that fail to
		contain at least one event preceding t could, at this stage,
		be discarded, however testing for that would require paying
		the price of the ID-->event look-up on all events in all
		candidates, which we will have to perform at each level of
		the graph, including the final stage before reporting
		candidates, so it is more efficient to leave the invalid
		candidates in the stream at this stage and cull them in a
		single pass at the end.

		The coincidences are reported as tuples of the Python IDs
		of the events that form coincidences.  The .index look-up
		table can be used to map these IDs back to the original
		events.  To do that, however, a copy of the contents of the
		.index mapping must be made before calling this method
		because this method will result in some of the events in
		question being removed from the queues, and thus also from
		the .index mapping.
		"""
		#
		# no-op unless time has advanced.  .previous_t is never
		# None, so flushing cannot be mistaken for a no-op
		#

		if t == self.previous_t:
			return tuple(), set()

		#
		# record the t to which we are being pulled for use in
		# defining the candidate interval in the next interation.
		# if we are being flushed, then reset to -inf.  this is
		# time B in the algorithm description above
		#

		self.previous_t = t if t is not None else segments.NegInfinity

		#
		# is this a "no-op" node?  return 1-detector "coincs" and a
		# null set of "partial coincs".  NOTE:  the "no-op" node
		# concept exists to enable single-detector searches, i.e.,
		# searches where there is not, in fact, any coincidence
		# analysis to perform.  it is not here to enable
		# single-detector candidate reporting in multi-detector
		# searches.  typically this streaming coincidence engine is
		# wrapped in additional code that updates the contents of
		# an output document based on the results, so by adding
		# this feature to the "coincidence engine", calling codes
		# don't need to support special configurations for
		# one-detector anlyses, they can set up their normal event
		# handling code with a one-detector offset vector, push
		# events in and collect candidates out.
		#

		if len(self.offset_vector) == 1:
			return tuple((id(event),) for event in self.components[0].pull(t)[0]), set()

		#
		# is this a leaf node?  construct the coincs explicitly
		#

		if len(self.offset_vector) == 2:
			#
			# search for and record coincidences.  coincs is a
			# sorted tuple of event ID pairs, where each pair
			# of IDs is, itself, ordered alphabetically by
			# instrument name.  here, for the two-instrument
			# case, the "partial coincs" are any singles that
			# failed to form coincidences, but we convert the
			# set into a set of 1-element tuples so they take
			# the form of 1-detector coincs
			#

			singles_ids = set()
			coincs = tuple(self.components[0].pull(t, singles_ids))
			return coincs, (set((eventid,) for eventid in singles_ids) if self.keep_partial else set())

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
		component_coincs_and_partial_coincs = tuple(component.pull(t) for component in self.components)
		component_coincs = tuple(elem[0] for elem in component_coincs_and_partial_coincs)

		if self.keep_partial:
			# any coinc with n-1 instruments from the component
			# time slides might fail to form an n instrument
			# coincidence.  we put them all into our
			# partial_coincs set now, and we'll remove them
			# from this set as we discover which particiapte in
			# n instrument coincidences
			partial_coincs = set(itertools.chain(*component_coincs))

			# the (< n-1)-instrument partial coincs are more
			# complicated.  for example, H+L doubles are
			# constructed as an intermediate step when forming
			# H+L+V triples and also when forming H+K+L
			# triples.  The H+L doubles that failed to form a
			# coincidence with V are reported as unused in the
			# former case, while the H+L doubles that failed to
			# form a coincidence with K are reported as unused
			# in the latter case, and there is no reason why
			# these two sets should be the same:  some H+L
			# doubles that were not used in forming H+L+V
			# triples might have been used to form H+L+K
			# triples.  if we are now using those two sets of
			# triples (among others) to form H+K+L+V
			# quadruples, we need to figure out which unused
			# H+L doubles have, in fact, really not been used
			# at this stage.  it can be shown that a (<
			# n-1)-instrument coinc will go unused in forming
			# any two of our (n-1)-instrument components if and
			# only if it is is not used to form any of our
			# (n-1)-instrument components.  therefore, we
			# obtain the unused (< n-1)-instrument partial
			# coincs from the union of all pair-wise
			# intersections of the (< n-1)-instrument partial
			# coincs from our components.  since each partial
			# coinc can only be listed at most once in the
			# result from each component, and because we must
			# compute all possible pair-wise intersections of
			# the sets, this is equivalent to finding the
			# partial coincs that appear two or more times in
			# the concatenated sequence.
			partial_coincs.update(coinc for coinc, count in Counter(itertools.chain(*(elem[1] for elem in component_coincs_and_partial_coincs))).items() if count >= 2)
		else:
			partial_coincs = set()
		del component_coincs_and_partial_coincs

		# magic:  we can form all n-instrument coincs by knowing
		# just three sets of the (n-1)-instrument coincs no matter
		# what n is (n > 2).
		coincs = []
		component_coincs0 = component_coincs[0]
		component_coincs1 = component_coincs[1]
		component_coincs2 = component_coincs[-1]
		del component_coincs
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

		return coincs, partial_coincs


class TimeSlideGraph(object):
	def __init__(self, coincgen_doubles_type, offset_vector_dict, coinc_window, min_instruments = 2, verbose = False):
		#
		# some initial input safety checks
		#

		if min_instruments < 1:
			raise ValueError("require min_instruments >= 1 (%d)" % min_instruments)
		if min(len(offset_vector) for offset_vector in offset_vector_dict.values()) < min_instruments:
			# this test is part of the logic that ensures we
			# will only extract coincs that meet the
			# min_instruments criterion, i.e., the condition
			# being checked here must be true or we will get
			# incorrect results but this is the only place it
			# is checked so be careful.  don't delete this
			# check.
			raise ValueError("encountered offset vector (%s) smaller than min_instruments (%d)", (str(min(offset_vector_dict.values(), key = lambda offset_vector: len(offset_vector))), min_instruments))

		#
		# for each pair of instruments compute the coincidence
		# window by adding the common coinc_window parameter to
		# the light travel time unique to each pair
		#

		coinc_windows = dict((frozenset(pair), coinc_window + light_travel_time(*pair)) for pair in itertools.combinations(set(instrument for offset_vector in offset_vector_dict.values() for instrument in offset_vector), 2))
		if verbose:
			print("coincidence windows:\n\t%s" % ",\n\t".join("%s = %g s" % ("+".join(sorted(pair)), dt) for pair, dt in coinc_windows.items()))

		#
		# populate the graph head nodes.  these represent the
		# target offset vectors requested by the calling code.
		#

		if verbose:
			print("constructing coincidence assembly graph for %d offset vectors ..." % len(offset_vector_dict), file=sys.stderr)
		self.head = tuple(
			TimeSlideGraphNode(
				coincgen_doubles_type,
				offset_vector,
				coinc_windows,
				min_instruments,
				time_slide_id = time_slide_id
			) for time_slide_id, offset_vector in sorted(offset_vector_dict.items())
		)

		# view into the id() --> event indexes of the graph nodes
		self.index = ChainMap(*(node.index for node in self.head))

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


	def push(self, instrument, events, t_complete):
		"""
		Push new events from some instrument into the internal
		queues.  The iterable of events need not be time-ordered.
		With self.singlesqueue.event_time() defining the time of
		events, t_complete sets the time up-to which the event
		stream from this instrument is known to be complete.  That
		is, in the future no new events will be pushed from this
		instrument whose times are earlier than t_complete.
		t_complete is measured with respect to the unshifted times
		of the events.

		Returns True if the graph's .t_coinc_complete has changed,
		indicating that it might be possible to form new
		candidates;  False if the addition of these events does not
		yet allow new candidates to be formed.  If the return value
		is False, .pull() is guaranteed to return an empty
		candidate list unless the graph is flushed of all
		candidates.
		"""
		t_changed = False
		for node in self.head:
			if instrument in node.offset_vector:
				t_before = node.t_coinc_complete
				node.push(instrument, events, t_complete)
				t_changed |= t_before != node.t_coinc_complete
		return t_changed


	def pull(self, newly_reported = None, flushed = None, flushed_unused = None, flush = False, coinc_sieve = None, event_collector = None):
		"""
		Extract newly available coincident candidates from the
		internal queues.  If flush if True, then all remaining
		candidates are extracted and the internal queues are reset,
		otherwise only those candidates that can be constructed
		given the times up to which the individual event queues are
		known to be complete and given the offset vectors being
		considered will be reported.

		This method returns a generator whose elements consist of
		(node, events) pairs, wherein events is a tuple of
		single-detector event objects comprising a single
		coincident n-tuple, and node is the TimeSlideGraphNode
		object from which the coincidence was retrieved.  The node
		object can be consulted to retrieve the offset vector
		applied to the event times to form the coincidence, and
		other information about the state of the analysis.  If the
		optional coinc_sieve() function is supplied (see below)
		then only those n-tuples that pass that function's test are
		included in the sequence.

		Two output streams are available.  One stream is the
		generator returned to the calling code explained above,
		which consists only of n-tuples that pass the optional
		coinc_sieve() test.  The other stream exits this code via
		calls to event_collector() (see below), and consists of all
		coincident n-tuples formed by the coincidence engine,
		regardless of whether they pass the optional coinc_sieve()
		test or not.  The n-tuples passed to event_collector()
		contain the Python IDs of the events, not the event objects
		themselves.  If coinc_sieve() is not supplied, both streams
		contain the same n-tuples, one of tuples of event objects,
		the other of tuples of the Python IDs of those events.

		event_collector(), if supplied, must be an object with a
		.push() method with the signature

			event_collector.push(event_ids, offset_vector)

		The return value is ignored.  It will be called once for
		every coincident n-tuple that is formed by the coincidence
		engine, regardless of whether or not that n-tuple is
		ultimately reported to the calling code.  event_ids will be
		a tuple of Python IDs of the single-detector events in the
		coincidence.

		coinc_sieve(), if supplied, must be a callable object with
		the signature

			coinc_sieve(events, offset_vector)

		If the return value is False the event tuple is reported to
		the calling code, otherwise it is discarded.  The default
		is to report all n-tuples.  This feature can be used to
		remove obviously uninteresting candidates from the reported
		candidate stream, to reduce the memory requirements of the
		calling code and the disk requirements of output documents.

		The optional parameters newly_reported, flushed, and
		flushed_unused may be used to pass lists to this code,
		which will be populated with additional information to
		assist the calling code.  If lists are supplied, their
		contents will be erased and populated with event objects as
		follows:

		newly_reported:  list of single-detector events that have
		participated in n-tuple candidates reported to the calling
		code for the first time this invocation of .pull().  Note
		that only those events that particate in n-tuples that
		survive the optional coinc_sieve() test are included.  This
		is meant to allow the calling code to populate an event
		table with the list of just the events needed to describe
		the surviving n-tuple candidates.

		flushed:  list of single-detector events that have been
		removed from the internal queues because they are now
		sufficiently old that all n-tuples that they can
		participate in have been formed and reported.  This can be
		used to assist the calling code with house keeping tasks,
		for example emptying its event_collector implementation of
		unnecessary internal state.

		flushed_unused:  list of the single-detector events in the
		flushed list that never participated in an n-tuple
		candidate.  Note that events are considered to have
		participated in an n-tuple candidate if the coincidence
		engine reported them in one, regardless of whether or not
		the optional coinc_sieve() test ultimately rejected the
		candidate before it was reported to the calling code.

		NOTE:  the intent of the coinc_sieve() feature is not to
		alter the definition of the coincidence test, but merely
		reduce the data volume going into the output document.
		Therefore, the coinc_sieve() feature affects the events
		appearing in the newly_reported list because that list is
		expected to be used to construct the output document, but
		does not affect the events appearing in the flushed_unused
		list because that list is meant to convey information about
		which events failed to form coincidenes, for example for
		the purpose of informing false-alarm probability noise
		models.
		"""
		# flatten ID index for faster performance in loop.  NOTE:
		# this is also done to freeze the contents of the index.
		# calling .pull() on the graph nodes flushes events from
		# the internal queues which removes them from the index,
		# but we need the index to turn IDs back into events.
		index = dict(self.index.items())

		# default coinc sieve
		if coinc_sieve is None:
			coinc_sieve = lambda events, offset_vector: False

		# default event_collector:
		if event_collector is None:
			event_collector_push = lambda event_ids, offset_vector: None
		else:
			event_collector_push = event_collector.push

		newly_reported_ids = set()
		# initialize the set of flushed events to the complete
		# contents of the internal queues
		flushed_ids = set(index)
		# avoid attribute look-ups in loops
		newly_reported_update = newly_reported_ids.update
		used_update = self.used_ids.update
		id_to_event = index.__getitem__
		for node in self.head:
			# avoid attribute look-ups loops
			event_time = node.event_time
			offset_vector = node.offset_vector

			if flush:
				t = None
				candidate_seg = segments.segment(node.previous_t, segments.PosInfinity)
			else:
				t = node.t_coinc_complete
				candidate_seg = segments.segment(node.previous_t, t)

			# we don't need to check that the coincs or partial
			# coincs contain at least min_instruments events
			# because those that don't meet the criteria are
			# excluded during coinc construction.
			for event_ids in itertools.chain(*node.pull(t)):
				# use the index to convert Python IDs back
				# to event objects
				events = tuple(map(id_to_event, event_ids))
				# check the candidate's time, and skip
				# those that don't meet the
				# t_coinc_complete constraint
				if min(event_time(event) + offset_vector[event.ifo] for event in events) not in candidate_seg:
					continue
				# update the used IDs state and inform the
				# event collector of a candidate
				used_update(event_ids)
				event_collector_push(event_ids, offset_vector)
				# apply the coinc sieve test.  if the
				# calling code has some mechanism to decide
				# which n-tuples are worth saving for later
				# analysis, it can be applied here (instead
				# of by the calling code itself) so that
				# the machinery that is tracking which
				# events need to be recorded in the output
				# document can see that some will not be
				# needed.
				if not coinc_sieve(events, offset_vector):
					newly_reported_update(event_ids)
					yield node, events
		# finish computing the set of newly reported event ids
		if newly_reported_ids:
			newly_reported_ids -= self.reported_ids
			self.reported_ids |= newly_reported_ids
		# finish computing the set of flushed events by removing
		# any from the set that still remain in the queues
		flushed_ids -= set(self.index)
		# compute the set of flushed events that were never used
		# for a candidate, and use flushed ids to clean up other
		# sets to avoid them growing without bound
		if flushed_ids:
			flushed_unused_ids = flushed_ids - self.used_ids
			self.used_ids -= flushed_ids
			self.reported_ids -= flushed_ids
		else:
			flushed_unused_ids = set()

		# if we've been flushed then there can't be any events left
		# in the queues (note that the "or" in parentheses is a
		# boolean operation, not a set operation)
		assert not flush or not (self.used_ids or self.reported_ids)

		# use the index to populate newly_used, flushed, and
		# flushed_unused with event objects
		if newly_reported is not None:
			newly_reported[:] = map(id_to_event, newly_reported_ids)
		if flushed is not None:
			flushed[:] = map(id_to_event, flushed_ids)
		if flushed_unused is not None:
			flushed_unused[:] = map(id_to_event, flushed_unused_ids)


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
#                        Poisson Model for Coincidences
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
#                         Ranking Statistic Components
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
		return ligolw.LIGO_LW({"Name": "%s:lnlrdensity" % name})

	@classmethod
	def get_xml_root(cls, xml, name):
		"""
		Sub-classes can use this in their overrides of the
		.from_xml() method to find the root element of the XML
		serialization.
		"""
		name = "%s:lnlrdensity" % name
		xml = [elem for elem in xml.getElementsByTagName(ligolw.LIGO_LW.tagName) if elem.hasAttribute("Name") and elem.Name == name]
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


class LnLikelihoodRatioMixin(object):
	"""
	Mixin to assist in implementing a log likelihood ratio ranking
	statistic class.  The ranking statistic class that inherits from
	this must:  (i) define a callable .numerator attribute that returns
	ln P(*args, **kwargs | signal);  (ii) define a callable
	.denominator attribute that returns ln P(*args, **kwargs | noise).

	Inheriting from this will:

	1.  Add a .__call__() method that returns the natural logarithm of
	the likelihood ratio

	ln P(*args, **kwargs | signal) - ln P(*args, **kwargs | noise)

	The implementation handles various special cases sensibly, such as
	when either or both of the logarithms of the numerator and
	denominator diverge.

	2.  Add a .ln_lr_samples() method that makes use of the .numerator
	and .denominator attributes, together with the .__call__() method
	to transform a sequence of (*args, **kwargs) into a sequence of log
	likelihood ratios and their respective relative frequencies.  This
	can be used to construct histograms of P(ln L | signal) and P(ln L
	| noise).  These distributions are required for, for example,
	signal rate estimation and false-alarm probability estimation.

	Why is the likelihood ratio useful?  Starting from Bayes' theorem,
	and using the fact that "signal" and "noise" are the only two
	choices:

	                   P(data | signal) * P(signal)
	P(signal | data) = ----------------------------
	                              P(data)

	                  P(data | signal) * P(signal)
	  = ---------------------------------------------------------
	    P(data | noise) * P(noise) + P(data | signal) * P(signal)

	            [P(data | signal) / P(data | noise)] * P(signal)
	  = ----------------------------------------------------------------
	    1 - P(signal) + [P(data | signal) / P(data | noise)] * P(signal)

	                        Lambda * P(signal)
	P(signal | data) = ----------------------------
	                   1 + (Lambda - 1) * P(signal)

	               P(data | signal)
	where Lambda = ----------------
	               P(data | noise)

	Differentiating P(signal | data) w.r.t. Lambda shows the derivative
	is always positive, so the probability that a candidate is the
	result of a gravitiational wave is a monotonically increasing
	function of the likelihood ratio.  Ranking events from "most likely
	to be a genuine signal" to "least likely to be a genuine signal"
	can be performed by sorting candidates by likelihood ratio.  Or, if
	one wanted to set a threshold on P(signal | data) to select a
	subset of candidates such that the candidates selected have a given
	purity (the fraction of them that are real is fixed), that is
	equivalent to selecting candidates using a likelihood ratio
	threshold.

	These are expressions of the Neyman-Pearson lemma which tells us
	that thresholding on Lambda creates a detector that extremizes the
	detection efficiency at fixed false-alarm rate.
	"""
	def __call__(self, *args, **kwargs):
		"""
		Return the natural logarithm of the likelihood ratio for
		the given parameters,

		ln P(*args, **kwargs | signal) - ln P(*args, **kwargs | noise)

		The arguments are passed verbatim to the .__call__()
		methods of the .numerator and .denominator attributes of
		self and the return value is computed from the results.

		NOTE:  sub-classes may override this method, possibly
		chaining to it if they wish.  The .ln_lr_samples()
		mechanism does not require this method to return exactly
		the natural logarithm of the .numerator/.denominator ratio.
		The .ln_lr_samples() mechanism does not assume the
		numerator, denominator and ranking statistic are related to
		each other as the latter being the ratio of the former two,
		it evaluates all three separately.  For this reason, the
		.__call__() method that implements the ranking statistic is
		free to include other logic, such as hierarchical cuts or
		bail-outs that are not stricly equivalent to the ratio of
		the numerator and denominator.

		Special cases:

		.numerator/.denominator=0/0 is mapped to ln Lambda = -inf,
		meaning P(signal | data) = 0.  Although this condition is
		nonsensical because the data is claimed to be inconsistent
		with both noise and signal --- the only two things it can
		be --- our task here is to compute something that is
		monotonic in the probability the data is the result of a
		signal, and since this data cannot be from a signal the
		correct answer is -inf.

		.numerator/.denominator = +inf/+inf is mapped to ln Lambda
		= NaN.  This is sufficiently nonsensical that there is no
		correct interpretation.  A warning will be displayed when
		this is encountered.
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
		Generator that transforms a sequence of candidate parameter
		samples into a sequence of log likelihood ratio samples.

		random_params_seq is a sequence (generator is OK) yielding
		3-element tuples whose first two elements provide the *args
		and **kwargs values to be passed to the .numerator and
		.denominator functions, and whose third element is the
		natural logarithm of the probability density from which the
		(*args, **kwargs) parameters have been drawn evaluated at
		those parameters.

		The output of this generator is a sequence of 3-element
		tuples, each of whose elements are:

		1.  a value of the natural logarithm of the likelihood
		ratio,

		2.  the natural logarithm of the relative frequency of
		occurance of that likelihood ratio in the signal population
		corrected for the relative frequency at which the
		random_params_seq sampler is causing that value to be
		returned, and

		3.  the natural logarithm of the relative frequency of
		occurance of that likelihood ratio in the noise population
		similarly corrected for the relative frequency at which the
		random_params_seq sampler is causing that value to be
		returned.

		The intention is for the first element of each tuple to be
		added to histograms using the two relative frequencies as
		weights, i.e., the two relative frequencies give the number
		of times one should consider this one draw of log
		likelihood ratio to have occured in the two populations.

		On each iteration, the *args and **kwargs values yielded by
		random_params_seq are passed to our own .__call__() method
		to evalute the log likelihood ratio at that choice of
		parameter values.  The parameters are also passed to the
		.__call__() mehods of our own .numerator and .denominator
		attributes to obtain the signal and noise population
		densities at those parameters.

		If signal_noise_pdfs is not None then, instead of using our
		own .numerator and .denominator attributes, the parameters
		are passed to the .__call__() methods of its .numerator and
		.denominator attributes to obtain those densities.  This
		allows the distribution of ranking statistic values
		obtained from alternative signal and noise populations to
		be modelled.  This is sometimes useful for diagnostic
		purposes.

		Normalizations:

		Within the context of the intended application, it is
		sufficient for all of the probabilities involved (the
		.numerator and .denominator probability densities, and the
		probability density supplied by the random_params_seq
		geneator) to be correct up to unknown normalization
		constants, i.e., the natural logarithms of the probabilties
		to be correct up to unknown additive constants.  That is
		why the two probability densities yielded by each iteration
		of this generator are described as relative frequencies:
		the ratios among their values are meaningful, but not their
		absolute values.

		If all of the supplied probabilities are, in fact, properly
		normalized, then the relative frequencies returned by this
		generator are, also, correctly normalized probability
		densities.
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
