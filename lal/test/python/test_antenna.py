#!/usr/bin/env python
#
# Copyright (C) 2018  Matthew Pitkin
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

"""
Test code for antenna.py.
"""

import sys
import numpy as np
from lal import DAYSID_SI
from lal.antenna import AntennaResponse


def test_antenna():
    """
    Test that the default, LAL and lookup table implementations match.
    """

    # set values
    psi = np.random.uniform(0., 2.*np.pi)
    ra = np.random.uniform(0., 2.*np.pi)
    dec = np.random.uniform(-np.pi/2., np.pi/2.)

    # set detectors
    detectors = ['H1', 'L1', 'V1']

    # set times over one sidereal day
    t0 = 1234567890.
    times = np.linspace(t0, t0 + DAYSID_SI, 1000)

    # test different detectors
    for det in detectors:
        # default antenna pattern (all tensor, vector and scalar modes)
        A = AntennaResponse(det, ra=ra, dec=dec, psi=psi, times=times,
                            scalar=True, vector=True)

        # LAL-wrapped antenna pattern
        B = AntennaResponse(det, ra=ra, dec=dec, psi=psi, times=times,
                            use_lal=True, scalar=True, vector=True)

        # look-up table
        C = AntennaResponse(det, ra=ra, dec=dec, psi=psi, times=times,
                            lookup=True, scalar=True, vector=True)

        for key in ['plus', 'cross', 'x', 'y', 'b', 'l']:
            # check values are the same (to within a small amount)
            if np.any(np.abs(A.response[key] - B.response[key]).astype('float16') != 0.):
                return False

            if np.any(np.abs(A.response[key] - C.response[key]) > 1e-3):
                return False

    return True


if test_antenna():
    sys.exit(0)
else:
    sys.exit(1)
