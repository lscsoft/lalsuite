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
Implementations of the general torque equation and other useful functions.
"""

import logging

import numpy as np


# The general torque equation we have formulated
def gte(t, f0, ngte, kgte):
    return f0 * (1 + (ngte - 1) * f0 ** (ngte - 1) * kgte * t) ** (1 / (1 - ngte))


# Inverse of the GTE
def gteinv(freq, f0, ngte, kgte):
    return ((freq / f0) ** (1 - ngte) - 1) / ((ngte - 1) * f0 ** (ngte - 1) * kgte)


# Returns value of kgte for a purely GW (ngte=5) source. Typical values: Izz \in [1, 3] * 10^38, ellip \in [1e-4, 1e-8].
# Izz default value is the same as that from the long transient search paper. Defaults give kgte=1.7182e-20
def kforGWsource(Izz=4.34e38, ellip=1e-4, radius=1e4):
    c = 299792458
    G = 6.6743e-11

    numerator = 32 * G * Izz * np.pi**4 * ellip**2
    denominator = 5 * c**5

    return numerator / denominator  # units: s^3


# The matrix which transforms from the PW parameters to tref parameters
def ParamTransformationMatrix(tstart, tend, reftime, s):

    dt = tend - tstart
    dref = reftime - tstart

    if s == 3:
        matrix = [
            [
                1
                - (6 * dref**5) / dt**5
                + (15 * dref**4) / dt**4
                - (10 * dref**3) / dt**3,
                dref
                - (3 * dref**5) / dt**4
                + (8 * dref**4) / dt**3
                - (6 * dref**3) / dt**2,
                dref**2 / 2.0
                - dref**5 / (2.0 * dt**3)
                + (3 * dref**4) / (2.0 * dt**2)
                - (3 * dref**3) / (2.0 * dt),
                (6 * dref**5) / dt**5
                - (15 * dref**4) / dt**4
                + (10 * dref**3) / dt**3,
                (-3 * dref**5) / dt**4
                + (7 * dref**4) / dt**3
                - (4 * dref**3) / dt**2,
                dref**5 / (2.0 * dt**3)
                - dref**4 / dt**2
                + dref**3 / (2.0 * dt),
            ],
            [
                (-30 * dref**4) / dt**5
                + (60 * dref**3) / dt**4
                - (30 * dref**2) / dt**3,
                1
                - (15 * dref**4) / dt**4
                + (32 * dref**3) / dt**3
                - (18 * dref**2) / dt**2,
                dref
                - (5 * dref**4) / (2.0 * dt**3)
                + (6 * dref**3) / dt**2
                - (9 * dref**2) / (2.0 * dt),
                (30 * dref**4) / dt**5
                - (60 * dref**3) / dt**4
                + (30 * dref**2) / dt**3,
                (-15 * dref**4) / dt**4
                + (28 * dref**3) / dt**3
                - (12 * dref**2) / dt**2,
                (5 * dref**4) / (2.0 * dt**3)
                - (4 * dref**3) / dt**2
                + (3 * dref**2) / (2.0 * dt),
            ],
            [
                (-120 * dref**3) / dt**5
                + (180 * dref**2) / dt**4
                - (60 * dref) / dt**3,
                (-60 * dref**3) / dt**4
                + (96 * dref**2) / dt**3
                - (36 * dref) / dt**2,
                1
                - (10 * dref**3) / dt**3
                + (18 * dref**2) / dt**2
                - (9 * dref) / dt,
                (120 * dref**3) / dt**5
                - (180 * dref**2) / dt**4
                + (60 * dref) / dt**3,
                (-60 * dref**3) / dt**4
                + (84 * dref**2) / dt**3
                - (24 * dref) / dt**2,
                (10 * dref**3) / dt**3
                - (12 * dref**2) / dt**2
                + (3 * dref) / dt,
            ],
            [
                (-360 * dref**2) / dt**5 + (360 * dref) / dt**4 - 60 / dt**3,
                (-180 * dref**2) / dt**4 + (192 * dref) / dt**3 - 36 / dt**2,
                (-30 * dref**2) / dt**3 + (36 * dref) / dt**2 - 9 / dt,
                (360 * dref**2) / dt**5 - (360 * dref) / dt**4 + 60 / dt**3,
                (-180 * dref**2) / dt**4 + (168 * dref) / dt**3 - 24 / dt**2,
                (30 * dref**2) / dt**3 - (24 * dref) / dt**2 + 3 / dt,
            ],
            [
                (-720 * dref) / dt**5 + 360 / dt**4,
                (-360 * dref) / dt**4 + 192 / dt**3,
                (-60 * dref) / dt**3 + 36 / dt**2,
                (720 * dref) / dt**5 - 360 / dt**4,
                (-360 * dref) / dt**4 + 168 / dt**3,
                (60 * dref) / dt**3 - 24 / dt**2,
            ],
            [
                -720 / dt**5,
                -360 / dt**4,
                -60 / dt**3,
                720 / dt**5,
                -360 / dt**4,
                60 / dt**3,
            ],
        ]

    elif s == 2:

        matrix = [
            [
                1 - (3 * dref**2) / dt**2 - (2 * dref**3) / dt**3,
                dref + (2 * dref**2) / dt + dref**3 / dt**2,
                (3 * dref**2) / dt**2 + (2 * dref**3) / dt**3,
                dref**2 / dt + dref**3 / dt**2,
            ],
            [
                (-6 * dref) / dt**2 - (6 * dref**2) / dt**3,
                1 + (4 * dref) / dt + (3 * dref**2) / dt**2,
                (6 * dref) / dt**2 + (6 * dref**2) / dt**3,
                (2 * dref) / dt + (3 * dref**2) / dt**2,
            ],
            [
                -6 / dt**2 - (12 * dref) / dt**3,
                4 / dt + (6 * dref) / dt**2,
                6 / dt**2 + (12 * dref) / dt**3,
                2 / dt + (6 * dref) / dt**2,
            ],
            [-12 / dt**3, 6 / dt**2, 12 / dt**3, 6 / dt**2],
        ]

    return matrix


# An alternative method to transoform our PW parameters to Doppler parameters in a taylor expansion with (t - tref) terms.
# Like the above methods, this is ill-conditioned, however less so than the above. Again similar to the above methods,
# this method works best if the start time (bf.knotslist[0]) of your data is subtracted from all time elements, p0, p1
# and tref. If you are using the conditioning in this method, it slows down dramatically. Good luck
def PWParamstoTrefParams(pwparams, p0, p1, tref, s):

    matrix = ParamTransformationMatrix(p0, p1, tref, s)

    logging.debug(len(matrix))
    logging.debug(len(pwparams))
    logging.debug(pwparams)
    # If no conditioning is required
    return np.matmul(matrix, pwparams)
