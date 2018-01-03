# Copyright (C) 2012 Duncan M. Macleod
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

"""Methods to calculate and manipulate trigger rates
"""

import numpy

from lal import (CreateREAL8TimeSeries, LIGOTimeGPS, lalHertzUnit)
from laldetchar import (git_version, triggers)

__author__ = "Duncan M. Macleod <duncan.macleod@ligo.org>"
__version__ = git_version.id
__date__ = git_version.date


def rate(table, stride, start=None, end=None):
    """@returns a TimeSeries of rate over time for all triggers in the
    given LIGO_LW table.
    """
    # get time
    tarray = triggers.get_time_column(table).astype(float)
    tarray.sort()
    # get limits
    if not start:
        start = tarray[0]
    if not end:
        end = tarray[-1]
    start = float(start)
    end = float(end)
    duration = end - start
    # contruct time bins
    stride = float(stride)
    duration = stride * round(duration/stride)
    bins = numpy.linspace(start, start+duration, num = duration//stride)
    # calculate rate
    rate = CreateREAL8TimeSeries("Rate (Hz)", LIGOTimeGPS(start), 0,
                                 stride, lalHertzUnit, bins.size-1)
    hist,_ = numpy.histogram(tarray, bins=bins)
    rate.data.data = hist.astype(numpy.float64)/stride
    return rate


def rate_per_bin(table, stride, column, bins, start=None, end=None):
    """@returns a list of TimeSeries representing the rate of events
    in each bin for the given LIGO_LW table
    """
    # get time
    tarray = triggers.get_time_column(table).astype(float)
    tarray.sort()
    # get limits
    if not start:
        start = tarray[0]
    if not end:
        end = tarray[-1]
    start = float(start)
    end = float(end)
    duration = end - start
    # contruct time bins
    stride = float(stride)
    duration = stride * round(duration/stride)
    bins = numpy.linspace(start, start+duration, num = duration//stride)
    # calculate rate per bin
    carray = triggers.get_column(str(column))
    out = []
    for bin_l,bin_r in bins:
        in_bin = (bin_l <= carray) & (carray < bin_r)
        rate = CreateREAL8TimeSeries("Rate (Hz) [%s <= %s < %s]"
                                     % (bin_l, column, bin_r),
                                     LIGOTimeGPS(start), 0,
                                     stride, lalHertzUnit, bins.size-1)
        hist,_ = numpy.histogram(tarray, bins=bins)
        rate.data.data = (hist / stride).astype(numpy.float64)
        out.append(rate)
    return out
