#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2020 Sebastian Khan
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http: //www.gnu.org/licenses/>.

"""Test PhenomD utility functions
"""

import sys
import pytest
import lal
import lalsimulation

# -- test functions ---------------------

FDamp0_test_data = [
    (100,1,2.3570523154029924e-21),
    (110,2,1.4260166508188107e-21),
    (120,3,1.1313851113934362e-21),
    (130,4,9.958546032577642e-22)
]

@pytest.mark.parametrize("mtot, distance, expected", FDamp0_test_data)
def test_SimPhenomUtilsFDamp0(mtot, distance, expected):
    """
    mtot: total mass in solar masses
    distance: distance to source in Mpc
    expected: expected value from lalsimulation.SimPhenomUtilsFDamp0
    """
    value = lalsimulation.SimPhenomUtilsFDamp0(mtot, distance*1e6*lal.PC_SI)

    assert value == expected, (
        "lalsimulation.SimPhenomUtilsFDamp0 test failed"
    )

# -- run the tests ------------------------------

if __name__ == '__main__':
    args = sys.argv[1:] or ["-v", "-rs", "--junit-xml=junit-phenomD_utils.xml"]
    sys.exit(pytest.main(args=[__file__] + args))
