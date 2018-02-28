#!/usr/bin/env python
# Copyright (C) 2010  Leo Singer
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
Benchmark script for pylal.low_latency_clustering.
"""
__author__ = "Leo Singer <leo.singer@ligo.org>"

from timeit import timeit


setup = """
from pylal.low_latency_clustering import clustered
import numpy as np

count_trials = 100000
keys = np.cumsum(np.random.random_integers(1, 10, count_trials))
vals = np.random.rand(count_trials)
indata = zip(keys, vals)
key = lambda x: x[0]
val = lambda x: x[1]
window = %g
"""

stmt = """
for item in clustered(indata, key, val, window):
	pass
"""


windows = range(10, 100, 10) + range(100, 1000, 100) + range(1000, 5000, 1000)
times = []

for window in range(10, 100, 10) + range(100, 1000, 100) + range(1000, 5000, 1000):
	print window
	times.append(timeit(stmt, setup % window, number=1) / 100000.)


import pylab
pylab.plot(windows, times)
pylab.xlabel('window duration')
pylab.ylabel('run time')
pylab.savefig('bench.pdf')
