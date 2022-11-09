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

import lalpulsar as lp
import lal

import lalpulsar.piecewise_model.BasisFunctions as bf
import lalpulsar.piecewise_model.ClassDefinitions as cd
import lalpulsar.piecewise_model.PWFStat as pwf
import lalpulsar.piecewise_model.PWModelSimulations as pwsim
import lalpulsar.piecewise_model.GTEandOtherMethods as gom
import lalpulsar.piecewise_model.EstimatingKnots as ek
import lalpulsar.piecewise_model.SemicoherentMetricMethods as scmm
import lalpulsar.piecewise_model.TBankEstimates as tbe

import argparse as ap
import numpy as np
import os
import sys
import logging
import cProfile
import pstats
import io
import signal

# Allow liblalpulsar C code to be interrupted with Ctrl+C
signal.signal(signal.SIGINT, signal.SIG_DFL)

parser = ap.ArgumentParser()

# Required parameters
parser.add_argument("--tbankcode",      type=str,       help="Specifies the tbank object to use for a search", required=True)
parser.add_argument("--j",              type=int,       help="The job array number on OzSTAR", required=True)

# Parameters used for specifying signals to inject in SFTs as well as other parameters for F-Statistic calculations. Do tstart and tref need to be ints?
parser.add_argument("--tstart",         type=int,       help="Start time of SFTs", default=900000000)
parser.add_argument("--trefsegfrac",    type=float,     help="What percentage tref occurs on each segment. E.g. if 1/2, then tref will be halfway through each segment", default=0.)
parser.add_argument("--dirname",        type=str,       help="The directory name to use for finding where the SFTs are if already generated", default=".")
parser.add_argument("--noiseratio",     type=float,     help="What fraction of h0 the noise in the SFTs will be. E.g. If set to 100, then the noise will be at h0 / 100", default=100)
parser.add_argument("--tempsperfile",   type=int,       help="The number of templates with the greatest fstats to record", default=1000)
parser.add_argument("--rtnsum",                         help="Whether to return the 2F stat averaged over all segments", action='store_true')
parser.add_argument("--SFTFiles",                       help="Whether to write SFT files or have them stored in memory", action='store_true')

# Parameters to be used when creating the necessary elements to use lp.ComputeFStat
parser.add_argument("--Tsft",           type=int,       help="The length of each SFT to be built", default=5)
parser.add_argument("--h0",             type=float,     help="Strain of the injected signal", default=1e-24)
parser.add_argument("--cosi",           type=float,     help="Default injection cos-iota", default=0.123)
parser.add_argument("--psi",            type=float,     help="Default injection psi", default=2.345)
parser.add_argument("--phi0",           type=float,     help="Default injection phi0", default=3.210)
parser.add_argument("--dt_wf",          type=float,     help="Sampling time for waveform generation", default=10)
parser.add_argument("--Alpha",          type=float,     help="Right Ascension to sky position of source", default=6.12)
parser.add_argument("--Delta",          type=float,     help="Declination to sky position of source", default=1.02)
parser.add_argument("--detector",       type=str,       help="Detector being used. Options H1, V1, etc.", default='H1')
parser.add_argument("--assume_sqrtSh",  type=float,     help="Assumed noise level for calculating 2F", default=1e-23)
parser.add_argument("--dfreq",          type=float,     help="Frequency spacing for ComputeFstat()", default=0)
parser.add_argument("--sourceDeltaT",   type=float,     help="sourceDeltaT option for ComputeFstat()", default=0.1)

# Parameters used in defining the parameter space
parser.add_argument("--s",              type=int,       help="Set the s parameter",                                                     nargs='?', const=1, default=-1)
parser.add_argument("--fmin",           type=float,     help="Set the fmin parameter",                                                  nargs='?', const=1, default=-1)
parser.add_argument("--fmax",           type=float,     help="Set the fmax parameter",                                                  nargs='?', const=1, default=-1)
parser.add_argument("--fmaxtrue",       type=float,     help="Set the fmaxtrue parameter",                                              nargs='?', const=1, default=-1)
parser.add_argument("--nmin",           type=float,     help="Set the nmin parameter",                                                  nargs='?', const=1, default=-1)
parser.add_argument("--nmax",           type=float,     help="Set the nmax parameter",                                                  nargs='?', const=1, default=-1)
parser.add_argument("--nmin0",          type=float,     help="Set the nmin0 parameter",                                                 nargs='?', const=1, default=-1)
parser.add_argument("--nmax0",          type=float,     help="Set the nmax parameter",                                                  nargs='?', const=1, default=-1)
parser.add_argument("--ntol",           type=float,     help="The time over which the braking index is allowed to change by 1%",        nargs='?', const=1, default=-1)
parser.add_argument("--taumin",         type=float,     help="Set the taumin parameter",                                                nargs='?', const=1, default=-1)
parser.add_argument("--taumax",         type=float,     help="Set the taumax parameter",                                                nargs='?', const=1, default=-1)
parser.add_argument("--ktol",           type=float,     help="The time over which the k value is allowed to change by 1%",              nargs='?', const=1, default=-1)
parser.add_argument("--dur",            type=float,     help="The maximum duration the knot algorithm to calculate up to",              nargs='?', const=1, default=-1)
parser.add_argument("--knots",          type=float,     help="Use user defined knots",                                                  nargs='+',          default=[0])
parser.add_argument("--maxmismatch",    type=float,     help="Set the mismatch parameter",                                              nargs='?', const=1, default=-1)
parser.add_argument("--maxtemps",       type=int,       help="The maximum number of templates to calculate using this tbank",                       default=-1)

# Parameters controlling outputs
parser.add_argument("--outbasedir",     type=str,       help="Output base directory",                                  default='.')
parser.add_argument("--logfile",                        help="If true, log to file; if false, log to standard output", action='store_true')
parser.add_argument("--loglevel",                       help="Level at which to log messages",                         default=logging.INFO)
parser.add_argument("--profile",                        help="Profile the code",                                       action='store_true')

args = parser.parse_args()
basedirectory = args.outbasedir

ek.setknotarchivepath(basedirectory)

# For logging
logging_format = '%(asctime)s : %(levelname)s : %(message)s'
if args.logfile:
        logging_filename = os.path.join(basedirectory, f'PiecewiseSearchLog_{j}.log')
        logging.basicConfig(filename=logging_filename, filemode='w', level=args.loglevel, format=logging_format)
        logging.info(f"Logging to {logging_filename}")
else:
        logging.basicConfig(level=args.loglevel, format=logging_format)
        logging.info(f"Logging to standard output")
logging.captureWarnings(True)
if logging.getLogger().isEnabledFor(logging.INFO):
        lp.globalvar.LatticeTilingProgressLogLevel = lal.LOG_NORMAL

# For profiling
if args.profile:
        logging.info("Making the profiler")
        pr = cProfile.Profile()
        pr.enable()

logging.debug("Doing the arg stuff")

# Required arguments
tbankcode     = args.tbankcode
j             = args.j

# Optional arguments
tstart        = args.tstart
trefsegfrac   = args.trefsegfrac
noiseratio    = args.noiseratio
rtnsum        = args.rtnsum
tempsperfile  = args.tempsperfile
SFTFiles      = args.SFTFiles

# Padding flags
padding_flags_bbox = args.flags_bbox
padding_flags_int  = args.flags_int
reset              = -1
if args.reset: reset = 1

# Setting the FInput object
finputdata               = cd.FInput()
finputdata.Tsft          = args.Tsft
finputdata.h0            = args.h0
finputdata.cosi          = args.cosi
finputdata.psi           = args.psi
finputdata.phi0          = args.phi0
finputdata.dt_wf         = args.dt_wf
finputdata.Alpha         = args.Alpha
finputdata.Delta         = args.Delta
finputdata.detector      = args.detector
finputdata.assume_sqrtSh = args.assume_sqrtSh
finputdata.dfreq         = args.dfreq
finputdata.sourceDeltaT  = args.sourceDeltaT

h0 = finputdata.h0

# Setting the template bank object and changing user specified parameters
tbank = cd.TBank()

if tbankcode == "GW170817":
        logging.info(f"Starting with default parameter space: {tbankcode}")
        tbank.SetDefaultBNSR()

if tbankcode == "1987A":
        logging.info(f"Starting with default parameter space: {tbankcode}")
        tbank.SetDefault1987A()

if tbankcode == "CasA":
        logging.info(f"Starting with default parameter space: {tbankcode}")
        tbank.SetDefaultCasA()

# Override any tbank arguments set on the command line
tbank.SetTBankParams(args)

# Checking if we are using user defined knots or if we need to recalculate the knots in case any changed parameters affect knot choice
if tbank.knots != [0]:
        logging.info("Setting knots from command line")
        bf.knotslist = args.knots
else:
        logging.info("Setting knots from algorithm")
        tbank.SetKnotsByAlg()

tbank.knots = bf.knotslist
tbank.dur = tbank.knots[-1]

# Adjusting knot start time
if bf.knotslist[0] != tstart:
        for i, knot in enumerate(bf.knotslist):
                bf.knotslist[i] = knot + tstart

logging.info("Knots: %s", str(bf.knotslist))

# Building the name of the directory where the SFTs will be saved. Directory name is based off the tbank object. In this way, any searches completed
# using the same tbank will have all SFTs in the same directory for neatness. Within the tbankdirectory, sub folders are used to contain the SFTs
# and results of each different search
tbankdirectory = os.path.join(basedirectory, tbank.toString())
if not os.path.isdir(tbankdirectory):
        os.mkdir(tbankdirectory)
logging.info(f"TBank directory is {tbankdirectory}")

# The name of the file that 2F and templates are stored to is only changed by whether the 2F values are averaged or only summed
sumstr = ""
if rtnsum:
        sumstr = "_Sum"

tempsfilename = "Temps" + sumstr + ".txt"
mismatchtempname = "Temps" + sumstr + "_Mismatch.txt"

# Build the SFTs and returns the random injected signal parameters. In the case that SFTFiles is false (the default) the SFTs are created in memory when
# the fstatinput method is called in PWFStat, however the randomly generate signal parameters are generated here

if SFTFiles:
        signalparams, SFTdirectory = pwsim.buildSFTs(tbank.dur, tbank, tstart=tstart, trefsegfrac=trefsegfrac, Tsft=finputdata.Tsft, parentdirectory=tbankdirectory, h0=h0, noiseratio=noiseratio)
else:
        logging.debug("Creating tiling lattice")
        tiling = lp.CreateLatticeTiling(tbank.s * len(bf.knotslist))

        logging.debug("Building metric")
        metric = scmm.PreCompMetric(tbank.s)

        logging.debug("Setting Bounds")
        tbe.setbounds(tiling, tbank)

        logging.debug("Setting tiling lattice and metric")
        lp.SetTilingLatticeAndMetric(tiling, lp.TILING_LATTICE_ANSTAR, metric, tbank.maxmismatch)

        logging.debug("Creating random signal params")
        randparams = lal.CreateRandomParams(0)
        signalparams = lal.gsl_matrix(tbank.s * len(bf.knotslist), 1)
        lp.RandomLatticeTilingPoints(tiling, 0, randparams, signalparams);                                      # Chance 'scale' parameter to be non-zero (within range [-1, 0])

        signalparams = np.transpose(signalparams.data)[0]

        logging.info("Random Signal Params are: %s", str(signalparams))

        f0 = signalparams[0]
        f1 = signalparams[1]
        f2 = signalparams[2]

        SFTdirectory = os.path.join(tbankdirectory, f"SFTs_h0-{h0:.2e}_f0-{f0:.3f}_f1-{f1:.2e}_f2-{f2:.2e}_dur-{tbank.dur}_tstart-{tstart}")

        if not os.path.isdir(SFTdirectory):
                os.mkdir(SFTdirectory)

# Create the text file where 2F and template data will be stored. The first line is added to the file which contains the column titles for all data.
with open(os.path.join(SFTdirectory, tempsfilename), "w") as reader:

        commentline = "{:20s}".format("#FStat") + "     " + str("{:20s}".format("Mismatch")) + "     "

        knotnum = -1

        for i in range(len(signalparams)):

                derivnum = i % tbank.s

                if derivnum == 0:
                        knotnum += 1

                columntitle = "f_" + str(knotnum) + "_" + str(i % 3)
                commentline = commentline + " " + "{:^20s}".format(columntitle)

        reader.write(commentline + "\n")

# Create the text file where 2F and template data will be stored but ordered by mismatch. The first line is added which contains the column titles for all data.
with open(os.path.join(SFTdirectory, mismatchtempname), "w") as reader:

        commentline = "{:20s}".format("#Mismatch") + "     " + str("{:20s}".format("FStat")) + "     "

        knotnum = -1

        for i in range(len(signalparams)):

                derivnum = i % tbank.s

                if derivnum == 0:
                        knotnum += 1

                columntitle = "f_" + str(knotnum) + "_" + str(i % 3)
                commentline = commentline + " " + "{:^20s}".format(columntitle)

        reader.write(commentline + "\n")

logging.debug("Setting Antenna Pattern")
# Set the maximum condition number for the antenna pattern above its default (get a warning using the piecewise model if not manually changed to be larger)
lp.SetAntennaPatternMaxCond(10 ** 5)

# The minimum and maximum frequencies needed to load in from SFTs to cover all frequencies any template may cover
SFTfmin = max(gom.gte(tbank.dur, tbank.fmin, tbank.nmax, gom.kwhichresultsingivenhalflife(tbank.taumin, tbank.fmax, tbank.nmin)), 50)
SFTfmax = tbank.fmax + 50
logging.info(f"SFTfmin/fmax: [{SFTfmin}, {SFTfmax}]")

# Build parameter space and begin calculated 2Fs
logging.info(f"SFT files is: {SFTFiles}")
pwf.semifstatcatalogue(SFTfmin, SFTfmax, tbank, finputdata, signalparams, tempsfilename, directory=SFTdirectory, trefsegfrac=trefsegfrac, rtnsum=rtnsum, SFTFiles=SFTFiles, tempsperfile=tempsperfile)

# For logging and profiling
if args.profile:
        pr.disable()
        s = io.StringIO()
        sortby = pstats.SortKey.CUMULATIVE
        ps = pstats.Stats(pr, stream=s).sort_stats('cumtime')
        ps.print_stats()
        with open(os.path.join(basedirectory, f'PiecewiseSearchProfile_{j}.txt'), 'w+') as f:
                f.write(s.getvalue())
