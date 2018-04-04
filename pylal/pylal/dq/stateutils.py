#!/usr/bin/env python

# Copyright (C) 2011 Duncan Macleod
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
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

"""
This module provides some utilities for manipulating data from a bit-mask state vector channel.
"""

from numpy import arange
from glue import segments

# =============================================================================
# Convert a data quality bit mask into segments
# =============================================================================

def _bits(i, n=8):
    """
    Convert integer bit mask into binary bits. Returns a list of 0s or 1s
    from 2^0 up to 2^n.

    Example:
    
    >>> _bits(295, n=8)
    [1, 1, 1, 0, 0, 1, 0, 0]
    """
    return [(0, 1)[int(i)>>j & 1] for j in xrange(n)]

def tosegmentlistdict(timeseries, bitmask):
    """
    Returns a glue.segments.segmentlistdict of active segments for each bit
    in a bit-masked state vector TimeSeries.
    """

    bits,flags = zip(*sorted(bitmask.items(), key=lambda (k,v): k))

    segdict = segments.segmentlistdict()
    for flag in flags:
        segdict[flag] = segments.segmentlist()

    # convert DQ bits into segments
    tarray = arange(timeseries.data.length) * float(timeseries.deltaT) +\
             float(timeseries.epoch)
    for t,d in zip(tarray.astype(float), timeseries.data.data):
        binary = _bits(d, bits[-1]+1)
        seg    = segments.segment(t, t-timeseries.deltaT)
        for bit,flag in zip(bits, flags):
            if binary[bit] == 1:
                segdict[flag].append(seg)

    segdict.coalesce()
    return segdict
