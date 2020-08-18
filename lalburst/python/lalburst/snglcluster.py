# Copyright (C) 2009,2010,2016,2017  Kipp Cannon
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


import math
import sys


from glue.text_progress_bar import ProgressBar
from lal import iterutils
from ligo import segments


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
from .git_version import date as __date__
from .git_version import version as __version__


#
# =============================================================================
#
#                               Clustering Tools
#
# =============================================================================
#


def smallest_enclosing_seg(a, b):
	"""
	Return the smallest segment that contains both a and b.
	"""
	return segments.segment(min(a[0], b[0]), max(a[1], b[1]))


def weighted_average_seg(seg1, weight1, seg2, weight2):
	"""
	Return the segment whose start and ends are the weighted arithmetic
	means of the start and ends of the two input segments, using the
	two weights given.

	The arithmetic means are computed in a way that is safe for
	LIGOTimeGPS objects.
	"""
	return segments.segment(seg1[0] + weight2 * float(seg2[0] - seg1[0]) / (weight1 + weight2), seg1[1] + weight2 * float(seg2[1] - seg1[1]) / (weight1 + weight2))


#
# =============================================================================
#
#                               Clustering Loop
#
# =============================================================================
#


def cluster_events(events, testfunc, clusterfunc, sortkeyfunc = None, bailoutfunc = None, verbose = False):
	"""
	Cluster the events in an event list.  testfunc will be passed a
	pair of events in random order, and must return 0 (or False) if
	they should be clustered.  clusterfunc will be passed a pair of
	events in random order, and must return an event that is the
	"cluster" of the two.  clusterfunc is free to return a new events,
	or modify one or the other of its parameters in place and return
	it.

	If sortkeyfunc and bailoutfunc are both not None (if one is
	provided the other must be as well), the events will be sorted into
	"increasing" order using sortkeyfunc as a sort key operator, and
	then only pairs of events for which bailoutfunc returns 0 (or
	False) will be considered for clustering.

	The return value is True if the events in the event list were
	modified, and False if they were not (although their order might
	have changed).
	"""
	# changed indicates if the event list has changed
	changed = False
	while True:
		if verbose:
			progress = ProgressBar("clustering %d events" % len(events), max = len(events))
			progress.show()
		else:
			progress = None

		if sortkeyfunc is not None:
			events.sort(key = sortkeyfunc)

		# outer_did_cluster indicates if the event list changes on
		# this pass
		outer_did_cluster = False
		i = 0
		while i < len(events):
			if progress is not None:
				progress.update(i)
			if events[i] is None:
				# this event has been clustered away
				i += 1
				continue
			# inner_did_cluster indicates if events[i] has
			# changed
			inner_did_cluster = False
			for j, event_j in enumerate(events[i + 1:], 1):
				if event_j is not None:
					if not testfunc(events[i], event_j):
						events[i] = clusterfunc(events[i], event_j)
						events[i + j] = None
						inner_did_cluster = True
					elif (sortkeyfunc is not None) and bailoutfunc(events[i], event_j):
						break
			if inner_did_cluster:
				outer_did_cluster = True
				# don't advance until events[i] stops
				# changing
			else:
				i += 1
		del progress
		# repeat until we do a pass without the listing changing
		if not outer_did_cluster:
			break
		iterutils.inplace_filter(lambda event: event is not None, events)
		changed = True
	return changed
