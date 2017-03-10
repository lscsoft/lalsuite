# Copyright (C) 2016  Nickolas Fotopoulos, Stephen Privitera, Ian Harry
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
%prog --nbanks N [options] input_bank.hdf5

Split a bank into smaller chirp mass regions based on the expected
computational scale for generating a template bank with the chosen
chirp mass boundaries.  The edges are taken from the input parameters
of the ligolw_cbc_sbank job that produced them unless overridden.
"""

import sys
import numpy
from bisect import bisect_left, bisect_right
from operator import attrgetter
from optparse import OptionParser
import h5py

def parse_command_line():
    parser = OptionParser(usage=__doc__)
    parser.add_option("-o", "--output-file", default=sys.stdout, help="Output list of N + 1 mchirp boundaries to this file (default: stdout)")
    parser.add_option("-n", "--nbanks", metavar="N", type = "int", help = "Set the number of subbanks to split the input bank.")
    parser.add_option("-w", "--template-weight", metavar="METHOD", default="equal", help = "Choose the method by which the chirp mass boundaries are chosen. By default, all templates are given equal weight, leading to subbanks of equal size. This choice is appropriate when using a metric for the match calculation. When a metric is not used, it is better to weight the templates by their duration, leading to finer splits at low mass where the templates are longer.", choices=("equal", "duration"))
    parser.add_option("-a", "--mchirp-min", metavar="MC_MIN", type="float", help="Override the lowest mchirp boundary to be MC_MIN")
    parser.add_option("-b", "--mchirp-max", metavar="MC_MAX", type="float", help="Override the lowest mchirp boundary to be MC_MAX")
    parser.add_option("-v", "--verbose", action="store_true", default=False, help="Be verbose.")
    options, filenames = parser.parse_args()

    if options.nbanks is None:
        parser.error("--nbanks is a required option")

    if options.nbanks <= 0:
        parser.error("--nbanks argument must be positive")

    if len(filenames) == 0:
        parser.error("must provide input file")

    if options.mchirp_min and options.mchirp_max and options.mchirp_max <= options.mchirp_min:
        parser.error("--mchirp-max argument must be greater than --mchirp-min")

    return options, filenames

options, filenames = parse_command_line()

# read input document
chirp_masses = []
for filename in filenames:
    hdf_fp = h5py.File(filename, 'r')
    mass1 = hdf_fp['mass1'][:]
    mass2 = hdf_fp['mass2'][:]
    mtot = mass1 + mass2
    eta = (mass1 * mass2) / (mtot*mtot)
    mchirp = mtot * eta**(3./5.)
    chirp_masses += list(mchirp)

if len(chirp_masses) < options.nbanks:
    raise ValueError("len(templates) < opts.nbanks; coarse bank too coarse")

# sort by mchirp
chirp_masses.sort()

# down-select if --mchirp-min or --mchirp-max are restrictive
if options.mchirp_min or options.mchirp_max:
    mchirps = [t.mchirp for t in templates]
    if options.mchirp_min:
        low = bisect_right(mchirps, options.mchirp_min)
    else:
        low = 0
    if options.mchirp_max:
        high = bisect_left(mchirps, options.mchirp_max)
    else:
        high = len(mchirps)
    chirp_masses = chirp_masses[low:high]


#
# Define boundaries based on templates in coarse bank. The step size
# is determined by a rough estimate of the computational scale for the
# bank generation between the given mchirp boundaries.
#
if options.template_weight == "equal":
    # More appropriate for metric match.
    template_weights = [1.0] * len(chirp_masses)
elif options.template_weight == "duration":
    raise NotImplementedError()
    # More appropriate for brute match (or bank sim). Split template
    # banks up to have equal amounts of sum( deltaF * deltaT ), which
    # should be ~ proportional to the computational scale for
    # generating this sub-bank.
    # template_weights = [row.tau0 * row.f_final for row in templates]

template_weights = numpy.array(template_weights).cumsum()
stride = numpy.searchsorted(template_weights, numpy.linspace(0, template_weights[-1], options.nbanks+1))

# We output the boundaries excluding the largest and smallest chirp
# mass values. The idea is that if the "coarse" bank already contains
# strict boundary points, then including extremal boundary points in
# the split boundary file will lead to a job with a zero-measure
# parameter space. This job gets stuck in a while 1 loop (this
# behavior is a bug, which should be caught by the code and lead to a
# crash, but regardless we want to avoid such cases here).
boundaries = [chirp_masses[k] for k in stride[1:-1]] # exclude high mass and low mass end points


# adjust based on user preferences
if options.mchirp_min:
    boundaries[0] = options.mchirp_min

# output
with open(options.output_file, "w") as outfile:
    print >>outfile, "\n".join(map(str, boundaries))
