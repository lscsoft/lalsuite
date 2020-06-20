# Test program for functions mapping between Equatorial and Detector-based coordinate systems
# Copyright (C) 2016 John Veitch <john.veitch@ligo.org>
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

import sys

import numpy as np
from numpy.testing import assert_allclose

import pytest

import lal
import lalinference as li

NTEST = 1000

LHO = lal.CachedDetectors[lal.LALDetectorIndexLHODIFF]
LLO = lal.CachedDetectors[lal.LALDetectorIndexLLODIFF]

RAS = np.random.uniform(low=0, high=lal.TWOPI, size=NTEST)
DECS = np.random.uniform(low=-lal.PI/2.0, high=lal.PI/2.0, size=NTEST)
TIMES = np.random.uniform(low=0, high=lal.DAYSID_SI, size=NTEST)


def test_invertable():
    res = np.empty((3, NTEST))
    for i, (ra, dec, time) in enumerate(zip(RAS, DECS, TIMES)):
        forward = li.EquatorialToDetFrame(LHO, LLO, time, ra, dec)
        res[:, i] = li.DetFrameToEquatorial(LHO, LLO, *forward)
    for a, b, tol in zip((TIMES, RAS, DECS), res, (1e-6, 1e-5, 1e-5)):
        assert_allclose(a, b, atol=tol)


if __name__ == '__main__':
    args = sys.argv[1:] or ["-v", "-rs", "--junit-xml=junit-detframe.xml"]
    sys.exit(pytest.main(args=[__file__] + args))
