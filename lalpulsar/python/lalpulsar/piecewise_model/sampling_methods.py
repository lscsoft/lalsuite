# Copyright (C) 2019--2023 Benjamin Grace
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

## \file
## \ingroup lalpulsar_python_piecewise_model
"""
Sample the general torque equation using various methods.
"""

import logging

import numpy as np

from . import basis_functions as bf
from . import gte_and_other_methods as gom
from . import errors


# Returns sample points where a 'ppint' number of points are evenly spaced between each knot
def pointsbyknot(ppint):
    ints = len(bf.knotslist) - 1
    points = []

    for i in range(ints):
        ps = bf.p(i)
        pe = bf.p(i + 1)
        for j in range(ppint):
            point = ps + j / ppint * (pe - ps)
            points.append(point)

    return points


# Returns sample points that correspond with values of frequency that are evenly spaced (number of points = ppint *
# ints)
def pointsbyfrequency(ppint, f0, ngte, kgte):
    points = []
    dur = bf.knotslist[-1] - bf.knotslist[0]
    ints = len(bf.knotslist) - 1

    freqps = np.linspace(
        gom.gte(0, f0, ngte, kgte), gom.gte(dur, f0, ngte, kgte), ppint * ints
    )

    for freq in freqps:
        points.append(gom.gteinv(freq, f0, ngte, kgte) + bf.knotslist[0])

    if points[-1] > bf.knotslist[-1]:
        points[-1] = bf.knotslist[-1]

    return points


# Returns sample points that are randomly chosen. Number of sample points = ppint * ints
def pointsrandom(ppint):
    points = []
    ints = len(bf.knotslist) - 1
    dur = bf.knotslist[-1] - bf.knotslist[0]

    for i in range(ppint * ints):
        point = np.random.uniform(0, dur)
        points.append(point)

    return points


normalpoints = False  # Points chosen to be evenly spaced between knots
freqpoints = True  # Points chosen in time to correspond with points evenly spaced in frequency by the GTE (runs into problems for long signal durations)

# Returns an integer value corresponding to the interval that point lies in
def thisint(point):
    if point < bf.knotslist[0] or point > bf.knotslist[-1]:
        logging.debug(
            "Sample point value not correctly chosen. Point and knot extrema are: "
            + str([point, bf.knotslist[0], bf.knotslist[-1]])
        )
        raise errors.PointNotWithinKnotBoundaries

    ints = len(bf.knotslist) - 1

    for i in range(ints + 1):
        ps = bf.p(i)

        if point < ps:
            return i - 1

    if point == bf.knotslist[-1]:
        return ints - 1


# Returns an integer list with numbers corresponding to which intervals do not contain sample points
def samplepointcheck(points, checkfirstintervals=True):

    if normalpoints:
        return []
    ints = len(bf.knotslist) - 1
    intscontainingnopoint = [i for i in range(ints)]

    for point in points:
        j = thisint(point)

        if j in intscontainingnopoint:
            intscontainingnopoint.remove(j)

    if not checkfirstintervals:
        if len(intscontainingnopoint) != 0:
            for j in range(ints - 1):
                if j in intscontainingnopoint:
                    intscontainingnopoint.remove(j)

    return intscontainingnopoint


# Returns our sample points depending upon the boolean values above
def samplepoints(ppint, f0, ngte, kgte):
    if normalpoints:
        points = pointsbyknot(ppint)
    elif freqpoints:
        points = pointsbyfrequency(ppint, f0, ngte, kgte)
    else:
        points = pointsrandom(ppint)

    sadints = samplepointcheck(points)

    if len(sadints) != 0:
        logging.error("These intervals contain no sample point: " + str(sadints))
        raise errors.SegmentContainsNoSamplePoints
    return points


# As the pointsbyfrequency method but now instead we choose sample points between the times corresponding with the
# knotnumbers knotnuma and knotnumb. Can be used for calculating the MOLS model of the GTE between two specific times
# instead of the default times 0 to the signal duration
def pointsbyfrequencywithinknots(knotnuma, knotnumb, ppint, f0, ngte, kgte):
    points = []

    ts = bf.knotslist[knotnuma] - bf.knotslist[0]
    te = bf.knotslist[knotnumb] - bf.knotslist[0]

    ints = knotnumb - knotnuma

    freqps = np.linspace(
        gom.gte(ts, f0, ngte, kgte), gom.gte(te, f0, ngte, kgte), ppint * ints
    )

    for freq in freqps:
        points.append(gom.gteinv(freq, f0, ngte, kgte) + bf.knotslist[0])

    if points[0] < ts + bf.knotslist[0]:
        points[0] = ts + bf.knotslist[0]
    if points[-1] > te + bf.knotslist[0]:
        points[-1] = te + bf.knotslist[0]

    return points


# As the the samplepoints method but now instead only generates sample points between the two given knots. As the above
# method, can be used for calculating the MOLS model between two specific knots instead of the default times 0 to the
# signal duration.
def samplepointswithinknots(knotnuma, knotnumb, ppint, f0, ngte, kgte):
    points = pointsbyfrequencywithinknots(knotnuma, knotnumb, ppint, f0, ngte, kgte)

    sadints = samplepointcheck(points, checkfirstintervals=False)

    if len(sadints) != 0:
        for thisint in sadints:
            if knotnuma + 1 <= thisint and thisint < knotnumb:
                logging.error(
                    "These intervals contain no sample point: " + str(sadints)
                )
                raise errors.SegmentContainsNoSamplePoints

    return points
