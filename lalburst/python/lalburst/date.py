# Copyright (C) 2006,2007,2009,2015-2017  Kipp Cannon
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
This module provides a few time-related utilities that have been used in
burst searches in the past.
"""


import math


from glue import segments
import lal


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
from .git_version import date as __date__
from .git_version import version as __version__


#
# =============================================================================
#
#                              Function Wrappers
#
# =============================================================================
#


def XLALTimeDelayFromEarthCenter(pos, ra, dec, gps):
	return lal.ArrivalTimeDiff(pos, (0.0, 0.0, 0.0), ra, dec, gps)


#
# =============================================================================
#
#                                Plotting Tools
#
# =============================================================================
#


def utc_midnight(gps):
	"""
	Truncate a LIGOTimeGPS to UTC midnight.
	"""
	# convert to UTC (as list so we can edit it)
	tm = list(lal.GPSToUTC(int(gps)))

	# truncate to midnight
	tm[3] = 0       # hours
	tm[4] = 0       # minutes
	tm[5] = 0       # seconds

	# convert back to LIGOTimeGPS
	return lal.LIGOTimeGPS(lal.UTCToGPS(tuple(tm)))


def UTCMidnights(start, end):
	"""
	Iterator for generating LIGOTimeGPS objects for UTC midnights.
	"""
	# 86402 == more seconds than there are in 1 day, but less than
	# there are in 2 days.  Can 1 day have more than 1 leap second?  If
	# so, this constant should be increased.
	midnight = utc_midnight(start)
	if midnight < start:
		midnight = utc_midnight(midnight + 86402)
	while midnight < end:
		yield midnight
		midnight = utc_midnight(midnight + 86402)


def gmst_0h(gps):
	"""
	Truncate a LIGOTimeGPS to Greenwich mean sidereal 0 rad.
	"""
	gmst = lal.GreenwichMeanSiderealTime(gps)
	residual = gmst % (2.0 * math.pi)
	if residual:
		gmst -= residual
	return lal.GreenwichMeanSiderealTimeToGPS(gmst)


def GMST_0hs(start, end):
	"""
	Iterator for generating LIGOTimeGPS objects for Greenwich Mean
	Sidereal 0h.
	"""
	midnight = gmst_0h(start)
	if midnight < start:
		midnight = gmst_0h(midnight + 86402)
	while midnight < end:
		yield midnight
		midnight = gmst_0h(midnight + 86402)


#
# =============================================================================
#
#                                Segment Lists
#
# =============================================================================
#


def gmst_days(gps_start, gps_stop):
	"""
	Generates a segmentlist whose segments are the Greenwich Mean
	Sidereal days spanning the given range of GPS times.  Input and
	output times are all GPS seconds.
	"""
	# construct an iterator for computing sequential 0h GMST in
	# LIGOTimeGPS.  86402 == same constant seen above in in
	# UTCMidnights.
	gmst_0hs = GMST_0hs(gmst_0h(gps_start), gps_stop + 86402)

	# initialize a segmentlist with the first sidereal day
	l = segments.segmentlist([segments.segment(gmst_0hs.next(), gmst_0hs.next())])

	# append each subsequent sideral day as another segment
	while l[-1][1] < gps_stop:
		l.append(segments.segment(l[-1][1], gmst_0hs.next()))

	# done
	return l
