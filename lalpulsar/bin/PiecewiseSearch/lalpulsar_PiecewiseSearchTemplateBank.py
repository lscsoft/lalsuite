# Copyright (C) 2019--2023 Benjamin Grace
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

import argparse as ap
import cProfile
import io
import logging
import pstats

import matplotlib.pyplot as plt
import numpy as np

import lalpulsar.piecewise_model.basis_functions as bf
import lalpulsar.piecewise_model.class_definitions as cd
import lalpulsar.piecewise_model.estimating_knots as ek
import lalpulsar.piecewise_model.gte_and_other_methods as gom
import lalpulsar.piecewise_model.tbank_estimates as tbe

# Initialise profiler
pr = cProfile.Profile()
pr.enable()

parser = ap.ArgumentParser()

parser.add_argument(
    "--tbankcode",
    type=str,
    help="Specifies the tbank object to use for a search",
    required=True,
)
parser.add_argument(
    "--j", type=int, help="The job array number on OzSTAR", required=True
)
parser.add_argument(
    "--SizeOrTemps",
    type=str,
    help="Use a capital S or T to specify whether you want to calculate the size of the template bank, or the templates themselves",
    required=True,
)

# Non-parameter space/non-tbank object optional arguments
parser.add_argument("--tstart", type=int, help="Start time of SFTs", default=900000000)
parser.add_argument(
    "--stats",
    help="Whether to write the statistics of the tiling to a file or not",
    action="store_true",
)

# Optional parameter space arguments
parser.add_argument(
    "--s", type=int, help="Set the s parameter", nargs="?", const=1, default=-1
)
parser.add_argument(
    "--fmin", type=float, help="Set the fmin parameter", nargs="?", const=1, default=-1
)
parser.add_argument(
    "--fmax", type=float, help="Set the fmax parameter", nargs="?", const=1, default=-1
)
parser.add_argument(
    "--nmin", type=float, help="Set the nmin parameter", nargs="?", const=1, default=-1
)
parser.add_argument(
    "--nmax", type=float, help="Set the nmax parameter", nargs="?", const=1, default=-1
)
parser.add_argument(
    "--kmin", type=float, help="Set the kmin parameter", nargs="?", const=1, default=-1
)
parser.add_argument(
    "--kmax", type=float, help="Set the kmax parameter", nargs="?", const=1, default=-1
)
parser.add_argument(
    "--dur",
    type=float,
    help="The maximum duration we want the knot algorithm to calculate up to",
    nargs="?",
    const=1,
    default=-1,
)
parser.add_argument(
    "--knots", type=float, help="Use user defined knots", nargs="+", default=[0]
)
parser.add_argument(
    "--mismatch",
    type=float,
    help="Set the mismatch parameter",
    nargs="?",
    const=1,
    default=-1,
)

args = parser.parse_args()

# Required arguments
tbankcode = args.tbankcode
j = args.j
SizeOrTemps = args.SizeOrTemps

# Optional arguments
tstart = args.tstart
stats = args.stats

# Logging
logging.basicConfig(
    filename="TBankCommandLineLog_" + str(j) + ".log",
    level=logging.DEBUG,
    filemode="w",
    format="%(asctime)s %(message)s",
)
log = logging.getLogger()
log.setLevel(logging.DEBUG)

# Constructing tbank object and setting parameters
tbank = cd.TBank()

if tbankcode == "GW170817":
    tbank.SetDefaultBNSR()

if tbankcode == "1987A":
    tbank.SetDefault1987A()

if tbankcode == "Small":
    tbank.SetSmallTestCase()

s = args.s
fmin = args.fmin
fmax = args.fmax
nmin = args.nmin
nmax = args.nmax
kmin = args.kmin
kmax = args.kmax
dur = args.dur
knots = args.knots
mismatch = args.mismatch

# If these parameters are not the set default value, then we make the corresponding change to the tbank object
if s != -1:
    tbank.s = s
if fmin != -1:
    tbank.fmin = fmin
if fmax != -1:
    tbank.fmax = fmax
if nmin != -1:
    tbank.nmin = nmin
if nmax != -1:
    tbank.nmax = nmax
if kmin != -1:
    tbank.kmin = kmin
if kmax != -1:
    tbank.kmax = kmax
if dur != -1:
    tbank.dur = dur
if mismatch != -1:
    tbank.maxmismatch = mismatch

# Checking if we are using user defined knots or if we need to recalculate the knots in case any changed parameters affect knot choice
if knots != [0]:
    bf.knotslist = knots
    tbank.dur = knots[-1] - knots[0]
else:
    ek.allidealisedknots(
        tbank.s, tbank.dur, 40, tbank.fmax, tbank.nmax, tbank.kmin, tbank.maxmismatch
    )

# Adjusting knot start time
if bf.knotslist[0] != tstart:
    for i, knot in enumerate(bf.knotslist):
        bf.knotslist[i] = knot + tstart

logging.info(tbank.toString())

logging.info("Knots are: " + str(bf.knotslist))
logging.info("Knots: %s", str(bf.knotslist))

tbank.knots = bf.knotslist

logging.info(bf.knotslist)

if SizeOrTemps == "S":
    temps = tbe.PWTBankSizeWithObject(tbank, stats=stats)
    logging.info("Number of temps found: " + str(temps))
elif SizeOrTemps == "T":
    temps = tbe.PWTBankWithObject(tbank)
    logging.info("Number of temps found: " + str(len(temps)))

    midway = int(np.ceil(len(temps) / 2))

    first = temps[0]
    midway = temps[midway]
    last = temps[-1]

    logging.info(
        "First:  "
        + str(first[0 : (2 * tbank.s)])
        + " ... "
        + str(first[-(2 * tbank.s + 1) : -1])
    )
    logging.info(
        "Midway: "
        + str(midway[0 : (2 * tbank.s)])
        + " ... "
        + str(midway[-(2 * tbank.s + 1) : -1])
    )
    logging.info(
        "Last:   "
        + str(last[0 : (2 * tbank.s)])
        + " ... "
        + str(last[-(2 * tbank.s + 1) : -1])
    )

    gom.PlotPWModel(first, show=False, label="First", linewidth=4)
    gom.PlotPWModel(midway, show=False, label="Midway", linewidth=3)
    gom.PlotPWModel(last, show=False, label="Last", linewidth=2)
    plt.legend()
    plt.show()

# Profiling
pr.disable()
s = io.StringIO()
sortby = pstats.SortKey.CUMULATIVE
ps = pstats.Stats(pr, stream=s).sort_stats("cumtime")
ps.print_stats()

with open("TBankCommandLineProfile_" + str(j) + ".txt", "w+") as f:
    f.write(s.getvalue())
