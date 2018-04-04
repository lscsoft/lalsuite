# Copyright (C) 2009  Kipp Cannon
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


from glue import segments
from glue import iterutils
from pylal import git_version


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


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


def cluster_events(events, testfunc, clusterfunc, sortfunc = None, bailoutfunc = None, verbose = False):
	"""
	Cluster the events in an event list.  testfunc will be passed a
	pair of events in random order, and must return 0 (or False) if
	they should be clustered.  clusterfunc will be passed a pair of
	events in random order, and must return an event that is the
	"cluster" of the two.  clusterfunc is free to return a new events,
	or modify one or the other of its parameters in place and return
	it.

	If sortfunc and bailoutfunc are both not None (if one is provided
	the other must be as well), the events will be sorted into
	"increasing" order using sortfunc as a comparison operator, and
	then only pairs of events for which bailoutfunc returns 0 (or
	False) will be considered for clustering.

	The return value is True if the events in the event list were
	modified, and False if they were not (although their order might
	have changed).
	"""
	changed = False
	while True:
		if verbose:
			print >>sys.stderr, "clustering pass:"

		if sortfunc is not None:
			if verbose:
				print >>sys.stderr, "\tsorting ..."
			events.sort(sortfunc)

		outer_did_cluster = False
		i = 0
		while i < len(events):
			if events[i] is None:
				i += 1
				continue
			if verbose and not (i % 13):
				print >>sys.stderr, "\t%d / %d%s\r" % (i + 1, len(events), " " * (int(math.floor(math.log10(len(events) or 1))) + 1)),
			inner_did_cluster = False
			for j in xrange(i + 1, len(events)):
				if events[j] is not None:
					if not testfunc(events[i], events[j]):
						events[i] = clusterfunc(events[i], events[j])
						events[j] = None
						inner_did_cluster = True
					elif (sortfunc is not None) and bailoutfunc(events[i], events[j]):
						break
			if inner_did_cluster:
				outer_did_cluster = True
			else:
				i += 1
		if verbose:
			print >>sys.stderr, "\t%d / %d%s" % (len(events), len(events), " " * (int(math.floor(math.log10(len(events) or 1))) + 1))
		if not outer_did_cluster:
			if verbose:
				print >>sys.stderr, "\tno change"
			break
		iterutils.inplace_filter(lambda event: event is not None, events)
		changed = True
	return changed
