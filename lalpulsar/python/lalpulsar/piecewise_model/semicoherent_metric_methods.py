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
Functions for building the semi-coherent metric for the piecewise model.
"""

import logging

import numpy as np
from scipy import integrate

from . import basis_functions as bf
from . import errors

# In this module we build the methods that generate the semi-coherent metric for our piecewise model. Each piecewise
# segment being a semi-coherent segment.

# This method checks whether an integral on the lth segment (starting from a segment number 0) will have a non-zero
# value if we integrate a basis function attached to the parameter with coordinates (i, borc) (ommitting the s value
# as this does not effect whether this integral is non-zero or not).
#
# Further on in this document we loop over all parameters for our signal model and conduct integrals over the basis
# functions for these parameters on each interval. This list of parameters we loop over have coordinates (i, k),
# where we by default have borc = 1. This is done because using the coordinates (i, borc, k) has redundancies. For
# example, the parameter coordinates (1, 0, k) is the same as (0, 1, k). To get rid of these redundancies we
# implement the new coordinates (i, k) with borc by default set as 1. To then check whether a parameter's basis
# function has a non-zero integral over the lth interval, if i == l then we know that that parameter is present on
# that interval. The other possibility is that i = l - 1 but we have borc = 1. We then check these two cases only,
# if either one mis True then we know that the associated integral will be non-zero.
#
#
# One might also note that if i = l + 1 and borc = 0 then that parameter will also appear on the lth interval and
# hence lead to a non-zero integral. This is true, however later in this module we have by default set borc = 1,
# thus if i = l + 1 by default borc = 1 and that parameter will not appear on the lth interval. Rather unintuitively,
# if borc = 1 by default then the first parameter actually appears on the -1th segment (despite this segment not
# being 'real').
def integralzerocheck(i, borc, l):
    if i == l or (i + 1 == l and borc == 1):
        return True
    else:
        return False


# The same as the integralzerocheck method except if we are integrating the product of two basis functions on intervals
# i1 and i2. Again included only for computational efficiency
def doubleintegralzerocheck(i1, i2, borc1, borc2, l):
    if integralzerocheck(i1, borc1, l) and integralzerocheck(i2, borc2, l):
        return True
    else:
        return False


# Integrates a basis function on the lth interval corresponding to coordinates l, borc and s and then multiplied by
# 2pi. For the remainder of this file l will refer to the interval we are working on, whereas a coordinate i will
# refer to the coordinate of the basis function or parameter
def integratebasisfunctiononl(pe, l, borc, s, coeffs):
    basisfunc = lambda t: bf.basisfunctionvalue(t, l, borc, s, coeffs)

    ps = bf.p(l)

    phasederiv = 2 * np.pi * integrate.quad(basisfunc, ps, pe, epsabs=0)[0]

    return phasederiv


# Calculates the time average for the part of our phase that belongs to the lth interval for a parameter with
# coordinates i, borc and k
def singlephaseaverageonl(l, i, borc, s, coeffs):
    if integralzerocheck(i, borc, l):
        if i == l:
            phasederiv = lambda t: integratebasisfunctiononl(t, l, borc, s, coeffs)
        else:
            phasederiv = lambda t: integratebasisfunctiononl(t, l, borc - 1, s, coeffs)
    else:
        phasederiv = 0

    # Changed ints to large value because changing of how we generate knots
    ps = bf.p(l)

    # Case put in for when we generate knots and the l + 1th knot is currently being generated and doesn't exist yet
    if l + 1 >= len(bf.knotslist):
        pe = bf.knotslist[-1]
    else:
        pe = bf.p(l + 1)

    return 1 / (pe - ps) * integrate.quad(phasederiv, ps, pe)[0]


# Calculates the product of two averaged phase functions for the given parameter coordinates on the lth interval
def singlephaseaveragesquaredonl(l, i1, i2, borc1, borc2, s1, s2, coeffs):
    if doubleintegralzerocheck(i1, i2, borc1, borc2, l):
        if i1 == l:
            parta = singlephaseaverageonl(l, l, borc1, s1, coeffs)
        else:
            parta = singlephaseaverageonl(l, l, borc1 - 1, s1, coeffs)

        if i2 == l:
            partb = singlephaseaverageonl(l, l, borc2, s2, coeffs)
        else:
            partb = singlephaseaverageonl(l, l, borc2 - 1, s2, coeffs)

        return parta * partb

    else:
        return 0


# Calculates the average of the product of two phase functions for the given parameter coordinates on the lth interval
def doublephaseaverageonl(l, i1, i2, borc1, borc2, s1, s2, coeffs):
    if doubleintegralzerocheck(i1, i2, borc1, borc2, l):
        if i1 == l:
            parta = lambda t: integratebasisfunctiononl(t, l, borc1, s1, coeffs)
        else:
            parta = lambda t: integratebasisfunctiononl(t, l, borc1 - 1, s1, coeffs)

        if i2 == l:
            partb = lambda t: integratebasisfunctiononl(t, l, borc2, s2, coeffs)
        else:
            partb = lambda t: integratebasisfunctiononl(t, l, borc2 - 1, s2, coeffs)

        # Changed ints to large value because changing of how we generate knots
        ps = bf.p(l)

        # Case put in for when we generate knots and the l + 1th knot is currently being generated and doesn't exist yet
        if l + 1 >= len(bf.knotslist):
            pe = bf.knotslist[-1]
        else:
            pe = bf.p(l + 1)

        doubleaverage = lambda x: parta(x) * partb(x)

        return 1 / (pe - ps) * integrate.quad(doubleaverage, ps, pe)[0]

    else:
        return 0


# Calculates the metric element corresponding to the parameter with the given coordinates on the lth interval
def metricelementonl(l, i1, i2, borc1, borc2, s1, s2, coeffs):
    singleaveragedsquared = singlephaseaveragesquaredonl(
        l, i1, i2, borc1, borc2, s1, s2, coeffs
    )
    doubleaverage = doublephaseaverageonl(l, i1, i2, borc1, borc2, s1, s2, coeffs)

    return doubleaverage - singleaveragedsquared


# Converts from the coordinates (i, k) to (i, borc, k)
def twocoordstothreecoords(i, s):
    if i == -1:
        return [0, 0, s]

    else:
        return [i, 1, s]


# Returns the list of coordinates that correspond to all parameters that will go into the computation of the metric
# without repetition
def coordslist(s):
    coords = []
    ints = len(bf.knotslist) - 1
    for i in range(-1, ints):
        for thiss in range(s):
            coords.append([i, thiss])

    return coords


# Calculates the metric associated with the lth interval. The appropriate coordinates for all parameters and
# coefficients for all basis functions must be given as input. Returns matrix with dimensions equal to the full
# semicoherent metric. To use the metric with the appropriate dimensions for the given interval, use the method
# metriconlcorrectdims
def metriconl(l, coords, coeffs):
    thismet = []

    newcoords = []

    for c in coords:
        newcoords.append(twocoordstothreecoords(c[0], c[1]))

    for nc in newcoords:
        thisrow = []
        for nnc in newcoords:
            metelem = metricelementonl(
                l, nc[0], nnc[0], nc[1], nnc[1], nc[2], nnc[2], coeffs
            )

            thisrow.append(metelem)

        thismet.append(thisrow)

    return np.array(thismet)


# Calculates the semicoherent metric for the given parameters
def metric(s):
    coeffs = bf.allcoeffs(s)
    coords = coordslist(s)

    metsonints = []
    ints = len(bf.knotslist) - 1
    for l in range(ints):
        metsonints.append(metriconl(l, coords, coeffs))

    thismet = np.zeros(np.shape(metsonints[0]))

    for met in metsonints:
        thismet += met

    return 1 / ints * thismet


# Returns the metric component for a single segment of length Tdata, symbolically computed in Mathematica.
# As the value of the metric should be independent of the start time, we do not need the knot values p0, p1,
# but only the length of the segment, Tdata. While it may appear that some of the float denominators are not
# exact (e.g. 1.89189e7), this is not the case and they are the exact value as calculated by Mathematica
def PreCompSingleSegMetric(Tdata, s):
    p1 = float(Tdata)
    Pi = np.pi

    if s == 3:
        firstrow = [
            (5975 * p1**2 * Pi**2) / 63063.0,
            (13859 * p1**3 * Pi**2) / 630630.0,
            (475 * p1**4 * Pi**2) / 252252.0,
            (9071 * p1**2 * Pi**2) / 126126.0,
            (-23333 * p1**3 * Pi**2) / 1.26126e6,
            (4259 * p1**4 * Pi**2) / 2.52252e6,
        ]
        secondrow = [
            (13859 * p1**3 * Pi**2) / 630630.0,
            (2764 * p1**4 * Pi**2) / 525525.0,
            (321 * p1**5 * Pi**2) / 700700.0,
            (23333 * p1**3 * Pi**2) / 1.26126e6,
            (-29713 * p1**4 * Pi**2) / 6.3063e6,
            (8077 * p1**5 * Pi**2) / 1.89189e7,
        ]
        thirdrow = [
            (475 * p1**4 * Pi**2) / 252252.0,
            (321 * p1**5 * Pi**2) / 700700.0,
            (127 * p1**6 * Pi**2) / 3.15315e6,
            (4259 * p1**4 * Pi**2) / 2.52252e6,
            (-8077 * p1**5 * Pi**2) / 1.89189e7,
            (233 * p1**6 * Pi**2) / 6.054048e6,
        ]
        fourthrow = [
            (9071 * p1**2 * Pi**2) / 126126.0,
            (23333 * p1**3 * Pi**2) / 1.26126e6,
            (4259 * p1**4 * Pi**2) / 2.52252e6,
            (5975 * p1**2 * Pi**2) / 63063.0,
            (-13859 * p1**3 * Pi**2) / 630630.0,
            (475 * p1**4 * Pi**2) / 252252.0,
        ]
        fifthrow = [
            (-23333 * p1**3 * Pi**2) / 1.26126e6,
            (-29713 * p1**4 * Pi**2) / 6.3063e6,
            (-8077 * p1**5 * Pi**2) / 1.89189e7,
            (-13859 * p1**3 * Pi**2) / 630630.0,
            (2764 * p1**4 * Pi**2) / 525525.0,
            (-321 * p1**5 * Pi**2) / 700700.0,
        ]
        sixthrow = [
            (4259 * p1**4 * Pi**2) / 2.52252e6,
            (8077 * p1**5 * Pi**2) / 1.89189e7,
            (233 * p1**6 * Pi**2) / 6.054048e6,
            (475 * p1**4 * Pi**2) / 252252.0,
            (-321 * p1**5 * Pi**2) / 700700.0,
            (127 * p1**6 * Pi**2) / 3.15315e6,
        ]

        metric = [firstrow, secondrow, thirdrow, fourthrow, fifthrow, sixthrow]

        return metric
    elif s == 2:
        firstrow = [
            (583 * p1**2 * Pi**2) / 6300.0,
            (223 * p1**3 * Pi**2) / 12600.0,
            (467 * p1**2 * Pi**2) / 6300.0,
            (-197 * p1**3 * Pi**2) / 12600.0,
        ]
        secondrow = [
            (223 * p1**3 * Pi**2) / 12600.0,
            (11 * p1**4 * Pi**2) / 3150.0,
            (197 * p1**3 * Pi**2) / 12600.0,
            (-41 * p1**4 * Pi**2) / 12600.0,
        ]
        thirdrow = [
            (467 * p1**2 * Pi**2) / 6300.0,
            (197 * p1**3 * Pi**2) / 12600.0,
            (583 * p1**2 * Pi**2) / 6300.0,
            (-223 * p1**3 * Pi**2) / 12600.0,
        ]
        fourthrow = [
            (-197 * p1**3 * Pi**2) / 12600.0,
            (-41 * p1**4 * Pi**2) / 12600.0,
            (-223 * p1**3 * Pi**2) / 12600.0,
            (11 * p1**4 * Pi**2) / 3150.0,
        ]

        return [firstrow, secondrow, thirdrow, fourthrow]

    else:
        logging.debug(
            "Given value of S has not yet been accounted for. Value of S given: "
            + str(s)
        )
        raise errors.ValueOfSNotAccountedFor()


def PreCompMetric(s):
    if s != 3 and s != 2:
        logging.debug(
            "Given value of S has not yet been accounted for. Value of S given: "
            + str(s)
        )
        raise errors.ValueOfSNotAccountedFor()

    segs = len(bf.knotslist) - 1

    segmentmetrics = []

    for seg in range(segs):
        Tseg = bf.knotslist[seg + 1] - bf.knotslist[seg]

        thismetric = PreCompSingleSegMetric(Tseg, s)

        segmentmetrics.append(thismetric)

    metric = np.zeros((s * (segs + 1), s * (segs + 1)))

    for i, thismetric in enumerate(segmentmetrics):

        prependzeros = [0 for j in range(i * s)]
        appendzeros = [0 for j in range((segs - 1 - i) * s)]
        zerorow = [0 for j in range((segs + 1) * s)]

        submetric = []

        prependzerorows = i * s
        appendzerorows = (segs - 1 - i) * s

        for j in range(prependzerorows):
            submetric.append(zerorow)

        for row in thismetric:
            thisrow = []

            thisrow = thisrow + prependzeros
            thisrow = thisrow + row
            thisrow = thisrow + appendzeros

            submetric.append(thisrow)

        for j in range(appendzerorows):
            submetric.append(zerorow)

        metric += np.array(submetric)

    return 1 / segs * metric
