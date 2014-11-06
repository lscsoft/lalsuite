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

def read_and_sort_toplist(filename,min_snr=float("inf"),
                          min_cands=-1,dropped_rho=float("-inf")):
    data = np.loadtxt(filename)
    rho = data[:,rho_col]
    sorted_inds = rho.argsort()[::-1]
    # If the toplist is longer than the maximum number of candidates,
    # truncate it and record the highest SNR value we discarded
    # unless that would bring us below the minimum snr to keep
    if min_cands > 0 and min_cands < len(sorted_inds):
        rho_to_drop = rho[sorted_inds[min_cands]]
        if rho_to_drop > min_snr:
            # count how many candidates are over the minimum
            keep_cands = np.sum(rho > min_snr)
            rho_to_drop = rho[sorted_inds[keep_cands]]
        else:
            keep_cands = min_cands
        if dropped_rho < rho[sorted_inds[keep_cands]]:
            dropped_rho = rho[sorted_inds[keep_cands]]
        sorted_inds = sorted_inds[:keep_cands]

    return data[sorted_inds,:], dropped_rho

parser = ArgumentParser()
parser.add_argument("--input-toplist-files", action="store", nargs="+",
                    required=True,
                    help='A space-separated list of toplist files to combine')
parser.add_argument("--output-toplist-file", action="store", required=True,
                    help='Filename for output toplist')
parser.add_argument("--min-cands-per-toplist", action="store", type=int,
                    default=-1,
                    help='Keep at least his many candidates from each toplist')
parser.add_argument("--min-snr-to-keep", action="store", type=float,
                    default=float("inf"),
                    help='Keep all candidates with at least this SNR')

args = parser.parse_args()

dropped_rho = float("-inf")

outfile = open(args.output_toplist_file,'w')
datatuple = tuple()
for filename in args.input_toplist_files:
    (newdata, dropped_rho) = \
        read_and_sort_toplist(filename=filename,
                              min_snr=args.min_snr_to_keep,
                              min_cands=args.min_cands_per_toplist,
                              dropped_rho=dropped_rho)
    datatuple += (newdata,)

data = np.concatenate(datatuple)
sorted_inds = data[:,rho_col].argsort()[::-1]
data = data[sorted_inds]

for line in data:
        outfile.write("%.10f %.10f %.10g %.10f %.10g %.5f %.10g %.10g %.10g\n"
                      % tuple(line))
print "Highest discarded candidate SNR was %f" % dropped_rho
outfile.close()
