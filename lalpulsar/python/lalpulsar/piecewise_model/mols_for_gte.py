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
Method Of Least Squares (MOLS) for matching a piecewise model to the general
torque equation (GTE).
"""

import logging

import numpy as np
from scipy import integrate

from . import basis_functions as bf
from . import gte_and_other_methods as gom
from . import errors
from . import sampling_methods as sm


# Our b vector for MOLS
def bvec(points, f0, ngte, kgte):
    return [gom.gte(t - bf.knotslist[0], f0, ngte, kgte) for t in points]


# Constructs an individual row of the a matrix for the sample point 'point' and the coefficients of our basis functions
# 'coeffs'
def rowofa(point, coeffs, s):
    j = sm.thisint(point)
    ints = len(bf.knotslist) - 1
    firstzeros = [0] * j * s
    lastzeros = [0] * (ints - 1 - j) * s

    coeffsb = []
    coeffsc = []

    for thiss in range(s):
        coeffb = bf.basisfunctionvalue(point, j, 0, thiss, coeffs)
        coeffc = bf.basisfunctionvalue(point, j, 1, thiss, coeffs)

        coeffsb.append(coeffb)
        coeffsc.append(coeffc)

    return firstzeros + coeffsb + coeffsc + lastzeros


# Builds our a matrix
def a(points, coeffs, s):
    amat = []

    for point in points:
        amat.append(rowofa(point, coeffs, s))

    return amat


# Returns a diagonal matrix containing the diagonal elements of mat to power -1/2
def pmat(mat):
    return np.diag(bf.listpseudoinv(np.diagonal(mat), 1 / 2))


# Returns mat multiplied by its transpose
def ata(mat):
    return np.matmul(np.transpose(mat), mat)


# Returns the solutions for our LSTSQ problem
def sols(coeffs, ppint, s, f0, ngte, kgte, conditioning=True):
    while True:
        try:
            points = sm.samplepoints(ppint, f0, ngte, kgte)
            break
        except errors.SegmentContainsNoSamplePoints:
            logging.debug("Reattempting MOLS fitting with greater sampling")
            ppint *= 2

    b = bvec(points, f0, ngte, kgte)
    amat = a(points, coeffs, s)

    atamat = ata(amat)

    if conditioning:
        pm = pmat(atamat)
    else:
        pm = np.identity(len(atamat))

    lhs = np.matmul(pm, np.matmul(atamat, pm))
    rhs = np.matmul(np.matmul(pm, np.transpose(amat)), b)

    logging.debug(
        "Condition number for matrix A  is: " + "{:.2E}".format(np.linalg.cond(amat))
    )
    logging.debug(
        "Condition number for matrix A' is: " + "{:.2E}".format(np.linalg.cond(lhs))
    )

    try:
        params = np.matmul(pm, np.linalg.solve(lhs, rhs))
    except np.linalg.LinAlgError as error:
        logging.debug(error)
        logging.debug(
            "Error in calculating MOLS parameters, using Python LSTSQ method instead"
        )
        params = np.matmul(pm, np.linalg.lstsq(lhs, rhs, rcond=-1)[0])

    return params


# Partitions the 1D list params into a 3D list for our parameters. Extract parameter by [ int ][ B or C ][ k ]
def solsbyint(params, s):
    ints = int(len(params) / s) - 1

    partedsols = np.zeros((ints, 2, s))

    for i in range(ints):
        solsb = np.array(params[i * s : (i + 1) * s])
        solsc = np.array(params[(i + 1) * s : (i + 2) * s])

        partedsols[i][0] = solsb
        partedsols[i][1] = solsc

    return partedsols


def solsbetweenknots(
    knotnuma, knotnumb, coeffs, ppint, s, f0, ngte, kgte, conditioning=True
):

    while True:
        try:
            points = sm.samplepointswithinknots(
                knotnuma, knotnumb, ppint, f0, ngte, kgte
            )
            break
        except errors.SegmentContainsNoSamplePoints:
            logging.debug("Reattempting MOLS fitting with greater sampling")
            ppint *= 2

    b = bvec(points, f0, ngte, kgte)
    amat = a(points, coeffs, s)

    for i in range(knotnuma * s):
        [row.pop(0) for row in amat]

    atamat = ata(amat)

    if conditioning:
        pm = pmat(atamat)
    else:
        pm = np.identity(len(atamat))

    lhs = np.matmul(pm, np.matmul(atamat, pm))
    rhs = np.matmul(np.matmul(pm, np.transpose(amat)), b)

    logging.debug(
        "Condition number for matrix A  is: " + "{:.2E}".format(np.linalg.cond(amat))
    )
    logging.debug(
        "Condition number for matrix A' is: " + "{:.2E}".format(np.linalg.cond(lhs))
    )

    try:
        params = np.matmul(pm, np.linalg.solve(lhs, rhs))
    except np.linalg.LinAlgError as error:
        logging.debug(error)
        logging.debug(
            "Error in calculating MOLS parameters, using Python LSTSQ method instead"
        )
        params = np.matmul(pm, np.linalg.lstsq(lhs, rhs)[0])

    zerobuff = [0] * (s * knotnuma)
    paramszerobuff = zerobuff + list(params)

    return paramszerobuff


points = []

# Calculates model value at given point for basis function coeffs 'coeffs' and parameter values 'params'
def modelvalueatpoint(point, coeffs, params, ignoreintcheck=False, singleseg=False):

    s = len(params[0][0])

    if ignoreintcheck:
        try:
            j = sm.thisint(point)
        except errors.PointNotWithinKnotBoundaries:

            t = point - bf.knotslist[0]

            if t < 0:
                f0 = params[0][0][0]
                f1 = params[0][0][1]
                f2 = 0

                if s > 2:
                    f2 = params[0][0][2]

                return f0 + f1 * t + 1 / 2 * f2 * t**2
            elif t > bf.knotslist[-1] - bf.knotslist[0]:
                f0 = params[-1][1][0]
                f1 = params[-1][1][1]
                f2 = 0

                if s > 2:
                    f2 = params[-1][1][2]

                return f0 + f1 * t + 1 / 2 * f2 * t**2

            return 0
    else:
        j = sm.thisint(point)

    # if singleseg:
    #    j = 0

    points.append(point)

    s = np.shape(params)[2]

    val = 0

    for thiss in range(s):
        basisfuncs = bf.basisfunctionvalue(point, j, 0, thiss, coeffs)
        basisfunce = bf.basisfunctionvalue(point, j, 1, thiss, coeffs)

        if singleseg:
            val += params[0][0][thiss] * basisfuncs
            val += params[0][1][thiss] * basisfunce
        else:
            val += params[j][0][thiss] * basisfuncs
            val += params[j][1][thiss] * basisfunce

    return val


def phase(point, coeffs, params, ignoreintcheck=False):

    model = lambda t: modelvalueatpoint(
        t, coeffs, params, ignoreintcheck=ignoreintcheck
    )
    # logging.debug("Integral bounds: " + str([bf.knotslist[0], point]))
    phasemodel = 2 * np.pi * integrate.quad(model, bf.knotslist[0], point, epsabs=0)[0]

    return phasemodel


# Returns the value of the error of our model at a given point
def errorvalueatpoint(point, coeffs, params, f0, ngte, kgte):
    error = np.abs(
        modelvalueatpoint(point, coeffs, params)
        - gom.gte(point - bf.knotslist[0], f0, ngte, kgte)
    )

    return error
