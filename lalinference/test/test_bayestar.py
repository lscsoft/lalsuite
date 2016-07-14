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
import sys
import doctest
import lalinference.bayestar.sky_map
import lalinference.bayestar.distance
import lalinference.bayestar.filter
import lalinference.fits
import lalinference.bayestar.timing

print('Running C unit tests.')
total_failures = lalinference.bayestar.sky_map.test()

print('Running Python unit tests.')
modules = [
    lalinference.bayestar.distance,
    lalinference.bayestar.filter,
    lalinference.fits,
    lalinference.bayestar.timing,
]
for module in modules:
    failures, tests = doctest.testmod(module)
    total_failures += failures

if total_failures > 0:
    sys.exit(1)
