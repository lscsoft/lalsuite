# Copyright (C) 2012  Nickolas Fotopoulos, Stephen Privitera
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
%prog --nbanks N [options] input_bank.xml

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
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils

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

    if len(filenames) != 1:
        parser.error("must provide input file")

    if options.mchirp_min and options.mchirp_max and options.mchirp_max <= options.mchirp_min:
        parser.error("--mchirp-max argument must be greater than --mchirp-min")

    return options, filenames[0]

options, filename = parse_command_line()

# read input document
xmldoc = utils.load_filename(filename, verbose=options.verbose)
templates = table.get_table(xmldoc, lsctables.SnglInspiralTable.tableName)
if len(templates) < options.nbanks:
    raise ValueError("len(templates) < opts.nbanks; coarse bank too coarse")

# sort by mchirp
templates.sort(key=attrgetter("mchirp"))

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
    templates = templates[low:high]

# define boundaries based on templates in coarse bank
if options.template_weight == "equal":
    template_weights = [1.0 for row in templates]
elif options.template_weight == "duration":
    template_weights = [row.tau0 * row.f_final for row in templates]
template_weights = numpy.array(template_weights).cumsum()

# split template banks up to have equal amounts of sum( deltaF *
# deltaT ), which should be ~ proportional to the computational scale
# for generating this sub-bank
stride = numpy.searchsorted(template_weights, numpy.linspace(0, template_weights[-1], options.nbanks))
boundaries = [templates[k].mchirp for k in stride[:-1]] # exclude high mass end point

# adjust based on user preferences
if options.mchirp_min:
    boundaries[0] = options.mchirp_min

# output
with open(options.output_file, "w") as outfile:
    print >>outfile, "\n".join(map(str, boundaries))
