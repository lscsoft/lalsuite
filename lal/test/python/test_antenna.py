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
from numpy.testing import assert_allclose

from lal import DAYSID_SI
from lal.antenna import AntennaResponse

import pytest

# set values
PSI = np.random.uniform(0., 2.*np.pi)
RA = np.random.uniform(0., 2.*np.pi)
DEC = np.random.uniform(-np.pi/2., np.pi/2.)

# set times over one sidereal day
T0 = 1234567890.
TIMES = np.linspace(T0, T0 + DAYSID_SI, 1000)


@pytest.mark.parametrize("detector", ("H1", "L1", "V1"))
def test_antenna(detector):
    """
    Test that the default, LAL and lookup table implementations match.
    """
    # default antenna pattern (all tensor, vector and scalar modes)
    A = AntennaResponse(detector, ra=RA, dec=DEC, psi=PSI, times=TIMES,
                        use_lal=False, scalar=True, vector=True)

    # LAL-wrapped antenna pattern
    B = AntennaResponse(detector, ra=RA, dec=DEC, psi=PSI, times=TIMES,
                        use_lal=True, scalar=True, vector=True)

    # look-up table
    C = AntennaResponse(detector, ra=RA, dec=DEC, psi=PSI, times=TIMES,
                        lookup=True, scalar=True, vector=True)

    for key in ['plus', 'cross', 'x', 'y', 'b', 'l']:
        assert_allclose(A.response[key], B.response[key])
        assert_allclose(A.response[key], C.response[key], atol=1e-3)


if __name__ == '__main__':
    args = sys.argv[1:] or ["-v", "-rs", "--junit-xml=junit-antenna.xml"]
    sys.exit(pytest.main(args=[__file__] + args))
