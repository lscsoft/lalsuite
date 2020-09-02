#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2020 Nathan K. Johnson-McDaniel
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

"""
This script checks that all implemented waveform
models have entries in:

XLALSimInspiralGetSpinSupportFromApproximant()
XLALSimInspiralGetSpinFreqFromApproximant()
XLALSimInspiralApproximantAcceptTestGRParams()
"""

import sys
import pytest
import lalsimulation as lalsim

# Function to get the name of a function
def get_name(f):
    return f.__name__

# Set of functions to test
func_list = [
    lalsim.SimInspiralGetSpinSupportFromApproximant,
    lalsim.SimInspiralGetSpinFreqFromApproximant,
    lalsim.SimInspiralApproximantAcceptTestGRParams
]

# Loop over all implemented TD or FD approximants

IMPLEMENTED = [
    k for k in range(lalsim.NumApproximants) if
    lalsim.SimInspiralImplementedTDApproximants(k) or
    lalsim.SimInspiralImplementedFDApproximants(k)
]

@pytest.mark.parametrize("func", func_list, ids=get_name)
@pytest.mark.parametrize("k", IMPLEMENTED, ids=lalsim.GetStringFromApproximant)
def test_wf_property_lists(k, func):
    func(k)

if __name__ == '__main__':
    args = sys.argv[1:] or ["-v", "-rs", "--junit-xml=junit-wf_property_lists.xml"]
    sys.exit(pytest.main(args=[__file__] + args))
