#
# Copyright (C) 2014-2016  Leo Singer
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
"""Construct a LIGO-LW XML power spectral density file for a network of
detectors by evaluating a model power noise sensitivity curve."""

import argparse
import inspect
import lal
import lalsimulation
from lalinference.bayestar import command

# Get names of PSD functions.
psd_name_prefix = 'SimNoisePSD'
psd_names = sorted(
    name[len(psd_name_prefix):]
    for name, func in inspect.getmembers(lalsimulation)
    if name.startswith(psd_name_prefix) and callable(func)
    and (
        '(double f) -> double' in func.__doc__ or
        '(REAL8FrequencySeries psd, double flow) -> int' in func.__doc__))

parser = command.ArgumentParser()
parser.add_argument(
    '-o', '--output', metavar='OUT.xml[.gz]', type=argparse.FileType('w'),
    default='-', help='Name of output file [default: stdout]')
parser.add_argument(
    '--df', metavar='Hz', type=float, default=1.0,
    help='Frequency step size [default: %(default)s]')
parser.add_argument(
    '--f-max', metavar='Hz', type=float, default=2048.0,
    help='Maximum frequency [default: %(default)s]')

# Add options for individual detectors
detectors = []
for detector in lal.CachedDetectors:
    name = detector.frDetector.name
    prefix = detector.frDetector.prefix
    detectors.append(prefix)
    parser.add_argument('--' + prefix, choices=psd_names, metavar='func',
        help='PSD function for {0} detector [optional]'.format(name))
    parser.add_argument('--' + prefix + '-scale', type=float, default=1.0,
        help='Scale range for {0} detector [default: %(default)s]'.format(name))

# Add list of vaild PSD functions.
parser.description += '''

The following options are supported for all detectors:

'''
for psd_name in psd_names:
    parser.description += '  ' + psd_name + '\n'

opts = parser.parse_args()


import glue.ligolw.utils
import lal.series
import os.path
import numpy as np
from lalinference.bayestar.timing import vectorize_swig_psd_func

# Add basic options.

psds = {}

n = int(opts.f_max // opts.df)
f = np.arange(n) * opts.df

for detector in detectors:
    psd_name = getattr(opts, detector)
    scale = 1 / np.square(getattr(opts, detector + '_scale'))
    if psd_name is None:
        continue
    func = getattr(lalsimulation, psd_name_prefix + psd_name)
    series = lal.CreateREAL8FrequencySeries(psd_name, 0, 0, opts.df, lal.SecondUnit, n)
    if '(double f) -> double' in func.__doc__:
        series.data.data = vectorize_swig_psd_func(psd_name_prefix + psd_name)(f)
    else:
        func(series, 0.0)

        # Find indices of first and last nonzero samples.
        nonzero = np.flatnonzero(series.data.data)
        first_nonzero = nonzero[0]
        last_nonzero = nonzero[-1]

        # Truncate
        series = lal.CutREAL8FrequencySeries(series, first_nonzero, last_nonzero - first_nonzero + 1)
        series.f0 = first_nonzero * series.deltaF

        series.name = psd_name
    series.data.data *= scale
    psds[detector] = series

glue.ligolw.utils.write_fileobj(
    lal.series.make_psd_xmldoc(psds), opts.output,
    gz=(os.path.splitext(opts.output.name)[-1]==".gz"))
