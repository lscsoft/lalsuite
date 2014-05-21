#!@PYTHON@
#
# Copyright (C) 2014  Leo Singer
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

import glue.ligolw.utils
import glue.text_progress_bar
import lal
import lal.series
import lalsimulation
import inspect
import optparse
import os.path
from lalinference.bayestar import command

# Get names of PSD functions.
psd_name_prefix = 'SimNoisePSD'
psd_names = sorted(
    name[len(psd_name_prefix):]
    for name, func in inspect.getmembers(lalsimulation)
    if name.startswith(psd_name_prefix) and callable(func)
    and '(double f) -> double' in func.__doc__)

# Create command line option parser.
parser = optparse.OptionParser(
    formatter=command.NewlinePreservingHelpFormatter(),
    description=__doc__,
    usage='%prog --H1=func ... [options] [-o FILENAME.xml[.gz]]')

# Add list of vaild PSD functions.
parser.description += '''
The following options are supported for all detectors:

'''
for psd_name in psd_names:
    parser.description += '  ' + psd_name + '\n'

# Add basic options.
parser.add_option(
    '-o', '--output', metavar='FILENAME.xml[.gz]', default='/dev/stdout',
    help='Name of optionally gzip-compressed output file [default: %default]')
parser.add_option(
    '--df', metavar='Hz', type=float, default=1,
    help='Frequency step size [default: %default]')
parser.add_option(
    '--f-max', metavar='Hz', type=float, default=2048,
    help='Maximum frequency [default: %default]')

detectors = []

# Add options for individual detectors
for detector in lal.CachedDetectors:
    name = detector.frDetector.name
    prefix = detector.frDetector.prefix
    detectors.append(prefix)
    parser.add_option('--' + prefix, choices=psd_names, metavar='func',
        help='PSD function for {0} detector [optional]'.format(name))

# Parse command line.
opts, args = parser.parse_args()
if args:
    parser.error('Did not expect any positional command line arguments')

psds = {}

unit = lal.Unit()
unit = lal.UnitInvert(unit, lal.HertzUnit)
n = int(opts.f_max // opts.df)
epoch = lal.LIGOTimeGPS()

progress = glue.text_progress_bar.ProgressBar()

for detector in detectors:
    psd_name = getattr(opts, detector)
    if psd_name is None:
        continue
    psd_func = getattr(lalsimulation, psd_name_prefix + psd_name)
    series = lal.CreateREAL8FrequencySeries(None, epoch, 0, opts.df, unit, n)
    fmt = '%s (%%d / %d)' % (detector, n)
    for i in progress.iterate(range(1, n), format=fmt):
        f = i * opts.df
        series.data.data[i] = psd_func(f)
    series.data.data[0] = series.data.data[1]
    psds[detector] = series

progress.update(-1, 'writing ' + opts.output)
glue.ligolw.utils.write_filename(
    lal.series.make_psd_xmldoc(psds), opts.output,
    gz=(os.path.splitext(opts.output)[-1]==".gz"))
