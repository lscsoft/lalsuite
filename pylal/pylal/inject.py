# Copyright (C) 2006--2011,2013  Kipp Cannon
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
This module provides a few convenience wrappers around bits of LAL's inject
package.
"""


import math


from lal import C_SI as _LAL_C_SI
from lal import CachedDetectors as _lalCachedDetectors
from pylal import git_version


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                                Look-up Tables
#
# =============================================================================
#


cached_detector_by_prefix = dict((cd.frDetector.prefix, cd) for cd in _lalCachedDetectors)
# make sure there were no duplicates
assert len(cached_detector_by_prefix) == len(_lalCachedDetectors)


cached_detector_by_name = dict((cd.frDetector.name, cd) for cd in _lalCachedDetectors)
# make sure there were no duplicates
assert len(cached_detector_by_name) == len(_lalCachedDetectors)

# FIXME:  backwards compatibility.  remove when no longer used
cached_detector = cached_detector_by_name


name_to_prefix = dict((name, detector.frDetector.prefix) for name, detector in cached_detector_by_name.items())
prefix_to_name = dict((prefix, name) for name, prefix in name_to_prefix.items())


# FIXME:  this is a hack to allow inject.light_travel_time(), which is used
# internally by snglcoinc.py, to work with the H1H2 coherent and null
# detectors.  the use of light_travel_time() internally by snglcoinc.py
# might be inappropriate and should be re-considered.  in the meantime,
# this will get the sub-solar mass search working.
prefix_to_name["H1H2"] = "LHO_4k"


#
# =============================================================================
#
#                              Function Wrappers
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
	dx = cached_detector_by_prefix[instrument1].location - cached_detector_by_prefix[instrument2].location
	return math.sqrt((dx * dx).sum()) / _LAL_C_SI
