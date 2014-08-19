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

def read_and_sort_toplist(filename,max_cands=-1,dropped_rho=0):
    data = np.loadtxt(filename,usecols=(freq_col,tp_col,asini_col,rho_col))

    my_freq_col = 0
    my_tp_col = 1
    my_asini_col = 2
    my_rho_col = 3

    rho = data[:,my_rho_col]
    sorted_inds = rho.argsort()[::-1]
    # If the toplist is longer than the maximum number of candidates,
    # truncate it and record the highest SNR value we discarded
    if max_cands > 0 and max_cands < len(sorted_inds):
        if dropped_rho < rho[sorted_inds[max_cands]]:
            dropped_rho = rho[sorted_inds[max_cands]]
        sorted_inds = sorted_inds[:max_cands]

    rho = rho[sorted_inds]
    freq = data[sorted_inds,my_freq_col]
    tp = data[sorted_inds,my_tp_col]
    asini = data[sorted_inds,my_asini_col]

    return freq, tp, asini, rho, dropped_rho

parser = ArgumentParser()
parser.add_argument("--input-toplist-files", action="store", nargs="*"
                    help='A space-separated list of toplist files to combine')
parser.add_argument("--output-toplist-files", action="store",
                    help='Filename for output toplist')
parser.add_argument("--max-cands-per-toplist", action="store", type=int,
                    help='Maximum number of candidates to keep from each toplist')

args = parser.parse_args()

data = parse_data(args.data_file)
