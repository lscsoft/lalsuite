#
# Copyright (C) 2015-2016  Leo Singer
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
"""
Run doctests for selected LALInference modules.
"""

# If we are missing Astropy and Healpy, then return exit code 77
# to mark the entire test as skipped.
import pkg_resources
try:
    pkg_resources.require(['astropy', 'healpy >= 1.9'])
except pkg_resources.ResolutionError as e:
    print(e)
    raise SystemExit(77)

import sys
import doctest
import numpy as np
import lalinference.bayestar.sky_map
import lalinference.bayestar.distance
import lalinference.bayestar.filter
import lalinference.io.fits
import lalinference.io.hdf5
import lalinference.healpix_tree
import lalinference.bayestar.timing
import lalinference.bayestar.postprocess

modules = [
    lalinference.bayestar.distance,
    lalinference.bayestar.filter,
    lalinference.bayestar.sky_map,
    lalinference.io.fits,
    lalinference.io.hdf5,
    lalinference.healpix_tree,
    lalinference.bayestar.timing,
    lalinference.bayestar.postprocess,
]


finder = doctest.DocTestFinder()
runner = doctest.DocTestRunner()
tests = []

for module in modules:
    # Find doctests in the module
    tests += finder.find(module)
    # Find doctests in Numpy external C ufuncs
    for ufunc in set(_ for _ in module.__dict__.values() if isinstance(_, np.ufunc)):
        tests += finder.find(ufunc, module=module)

total_failures = 0
for test in tests:
    failures, tests = runner.run(test)
    total_failures += failures

if total_failures > 0:
    sys.exit(1)
