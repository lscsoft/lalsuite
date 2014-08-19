#!/usr/bin/env python
# script to combine toplists from lalapps_pulsar_crosscorr_v2

# Copyright (C) 2014 John Whelan

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with with program; see the file COPYING. If not, write to the
# Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA

from __future__ import division
from argparse import ArgumentParser
import numpy as np

freq_col = 0
tp_col = 1
argp_col = 2
asini_col = 3
ecc_col = 4
period_col = 5
estSens_col = 6
evSquared_col = 7
rho_col = 8

def read_and_sort_toplist(filename,max_cands=-1,dropped_rho=float("-inf")):
    data = np.loadtxt(filename)
    rho = data[:,rho_col]
    sorted_inds = rho.argsort()[::-1]
    # If the toplist is longer than the maximum number of candidates,
    # truncate it and record the highest SNR value we discarded
    if max_cands > 0 and max_cands < len(sorted_inds):
        if dropped_rho < rho[sorted_inds[max_cands]]:
            dropped_rho = rho[sorted_inds[max_cands]]
        sorted_inds = sorted_inds[:max_cands]

    return data[sorted_inds,:], dropped_rho

parser = ArgumentParser()
parser.add_argument("--input-toplist-files", action="store", nargs="+",
                    required=True,
                    help='A space-separated list of toplist files to combine')
parser.add_argument("--output-toplist-file", action="store", required=True,
                    help='Filename for output toplist')
parser.add_argument("--max-cands-per-toplist", action="store", type=int,
                    default=-1,
                    help='Maximum number of candidates to keep from each toplist')

args = parser.parse_args()

dropped_rho = float("-inf")

outfile = open(args.output_toplist_file,'w')
for filename in args.input_toplist_files:
    (data,
     dropped_rho) = read_and_sort_toplist(filename=filename,
                                          max_cands=args.max_cands_per_toplist,
                                          dropped_rho=dropped_rho)
    for line in data:
        outfile.write("%.10f %.10f %.10g %.10f %.10g %.5f %.10f %.10g %.10g\n"
                      % tuple(line))
print "Highest discarded candidate SNR was %10.g" % dropped_rho
outfile.close()
