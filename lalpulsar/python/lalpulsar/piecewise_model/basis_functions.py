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
Build the methods required for creating the basis functions of our piecewise model.
"""

import logging

import numpy as np

# It is worth noting here that initially conditioning of the matrices used here
# was undertaken. Initially this was necessary as when plotting the basis
# functions all that was returned was numerical fuzz, however many months after
# having written these methods, I cannot recover this behaviour and the
# conditioning used here seems to increase the condition number of the matrices
# built here. For now the conditioning methods in this module have been
# commented out until a later time, however for now without conditioning the
# results produced appear to be returned accurately despite the
# ill-conditionedness of the problem.

knotslist = [0.0]

# Returns the value of the ith knot. If we are using the methods in the estimating_knots module, this method should
# simply extract the ith element from the knotslist variable above.
def p(i):

    return knotslist[i]


# Our basis function
def u(t, i):
    return (t - p(i)) / (p(i + 1) - p(i))


# Builds the W matrix whose inverse is the coefficients of our basis function
def w(i, s):
    try:
        from sympy import Symbol, diff, lambdify
    except:
        raise ImportError("Cannot import sympy")

    t = Symbol("t")

    ps = p(i)
    pe = p(i + 1)

    matrixupper = []
    matrixlower = []

    thisrowderivs = [u(t, i) ** m for m in range(2 * s)]

    for row in range(s):

        if row != 0:
            for j, elem in enumerate(thisrowderivs):
                thisrowderivs[j] = diff(elem, t, 1)

        thisrows = []
        thisrowe = []

        for elem in thisrowderivs:
            thiselem = lambdify(t, elem, "numpy")

            thisrows.append(thiselem(ps))
            thisrowe.append(thiselem(pe))

        matrixupper.append(thisrows)
        matrixlower.append(thisrowe)

    return np.array(matrixupper + matrixlower)


# Returns the pseudo inverse of a list to the given power
def listpseudoinv(lst, pwr=1):
    invlist = []

    for elem in lst:
        if elem != 0:
            invlist.append(elem**-pwr)
        else:
            invlist.append(0)

    return invlist


# Builds a diagonal matrix with elements equal to the inverse of the diagonal elements of the given matrix
def d(mat):
    return np.diag(listpseudoinv(np.diagonal(mat)))


# Returns the coefficients of our basis function in a 3D list. Reference elements by [ B or C ][ k ][ specific coeff ]
def basiscoeffs(i, s, conditioning=True):
    wmat = w(i, s)

    if conditioning:
        dmat = d(wmat)
        dwmat = np.matmul(dmat, wmat)

    else:
        dmat = np.identity(len(wmat))
        dwmat = wmat

    try:
        coeffs = np.transpose(np.linalg.solve(dwmat, dmat))
    except np.linalg.LinAlgError as error:
        logging.error(error)
        logging.error(
            "Error in calculating basis functions. Using Python LSTSQ method instead"
        )
        coeffs = np.linalg.lstsq(dwmat, dmat)[0]

    blist = coeffs[0:s]
    clist = coeffs[s : 2 * s]

    return np.array([blist, clist])


# Returns the coefficients of all basis functions in a 4D list. Reference elements by [ int ][ B or C][ k ][ specific
# coeff ]
def allcoeffs(s):
    coeffs = [basiscoeffs(i, s) for i in range(len(knotslist) - 1)]

    return np.array(coeffs)


# allcoeffs(3, 10, 10000)

# Returns the value of a specified basis function given the 4D list of coefficients coeffs.
def basisfunctionvalue(t, i, borc, s, coeffs):
    val = 0

    if t < p(i) or p(i + 1) < t:
        return 0

    thesecoeffs = coeffs[i][borc][s]

    for m, coeff in enumerate(thesecoeffs):
        val += coeff * u(t, i) ** m

    return val
