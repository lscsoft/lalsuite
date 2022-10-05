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
import lalpulsar.piecewise_model.MOLSforGTE as mols

import argparse as ap
import numpy as np
import os
import logging
import cProfile, pstats, io

print("Making the profiler")

pr = cProfile.Profile()
pr.enable()

# Change ints I am using as Bools to bools
# Change  no signal names by including H_0 parameters

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
parser.add_argument("--Tsft",           type=int,       help="The length of each SFT to be built", default=60)
parser.add_argument("--h0",             type=float,     help="Strain of the injected signal", default=1e-24)
parser.add_argument("--cosi",           type=float,     help="", default=0.123)
parser.add_argument("--psi",            type=float,     help="", default=2.345)
parser.add_argument("--phi0",           type=float,     help="", default=3.210)
parser.add_argument("--dt_wf",          type=float,     help="", default=10)
parser.add_argument("--Alpha",          type=float,     help="Right Ascension to sky position of source", default=6.12)
parser.add_argument("--Delta",          type=float,     help="Declination to sky position of source", default=1.02)
parser.add_argument("--detector",       type=str,       help="Detector being used. Options H1, V1, uhhhhh", default='H1')
parser.add_argument("--assume_sqrtSh",  type=float,     help="Assumed noise level for calculating 2F", default=1e-23)
parser.add_argument("--dfreq",          type=float,     help="", default=0)
parser.add_argument("--sourceDeltaT",   type=float,     help="", default=2)

# Parameter space padding flags
parser.add_argument("--flags_bbox",     type=float,     help="List of what fraction of bounding box to use for padding flags",          nargs='+',          default=[0])
parser.add_argument("--flags_int",      type=int,       help="List of what number of templates to use as padding flags",                nargs='+',          default=[0])
parser.add_argument("--reset",                          help="Whether or not to use the reseting methods in PiecewiseModel. If Either flags_bbox or flags_int are non-zero, this will automatically be set to true", action='store_true')

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

args = parser.parse_args()

# Required arguments
tbankcode     = args.tbankcode
j             = args.j

print("Doing the arg stuff")

# Optional arguments
tstart        = args.tstart
trefsegfrac   = args.trefsegfrac
directoryname = args.dirname
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

print("Making the logger")
# For logging
logging.basicConfig(filename='PWFCommandLineLog_' + str(j) + '.log', level=logging.DEBUG, filemode='w', format='%(asctime)s %(message)s')
log = logging.getLogger()
log.setLevel(logging.DEBUG)

# Setting the template bank object and changing user specified parameters
tbank = cd.TBank()

if tbankcode == "GW170817":
        tbank.SetDefaultBNSR()

if tbankcode == "1987A":
        tbank.SetDefault1987A()

if tbankcode == "Small":
        tbank.SetSmallTestCase()
print(tbank.dur)

s            = args.s
fmin         = args.fmin
fmax         = args.fmax
fmaxtrue     = args.fmaxtrue
nmin         = args.nmin
nmax         = args.nmax
nmin0        = args.nmin0
nmax0        = args.nmax0
ntol         = 0.01 / args.ntol
taumin       = args.taumin
taumax       = args.taumax
ktol         = 0.01 / args.ktol
dur          = args.dur
knots        = args.knots
max_mismatch = args.maxmismatch
maxtemps     = args.maxtemps

# Changing the appropriate parameters
if s                != -1: tbank.s        = s
if fmin             != -1: tbank.fmin     = fmin
if fmax             != -1: tbank.fmax     = fmax
if fmaxtrue         != -1: tbank.fmaxtrue = fmaxtrue
if nmin             != -1: tbank.nmin     = nmin
if nmax             != -1: tbank.nmax     = nmax
if nmin0            != -1: tbank.nmin0    = nmin0
if nmax0            != -1: tbank.nmax0    = nmax0
if ntol             != -0.01: tbank.ntol  = ntol
if taumin           != -1: tbank.taumin   = taumin
if taumax           != -1: tbank.taumax   = taumax
if ktol             != -0.01: tbank.ktol  = ktol
if dur              != -1: tbank.dur      = dur
if max_mismatch     != -1: tbank.mismatch = max_mismatch
if maxtemps         != -1: tbank.maxtemps = maxtemps

print("Setting knots in PWFCommnadLine")

# Checking if we are using user defined knots or if we need to recalculate the knots in case any changed parameters affect knot choice
if knots != [0]:
        bf.knotslist = knots
        tbank.dur = knots[-1]
else:
        ek.allidealisedknots(tbank.s, tbank.dur, 40, tbank.fmaxtrue, tbank.nmax, tbank.taumin, tbank.mismatch)

# Adjusting knot start time
if bf.knotslist[0] != tstart:
        for i, knot in enumerate(bf.knotslist):
                bf.knotslist[i] = knot + tstart

print("Knots are: " + str(bf.knotslist))
logging.info("Knots: %s", str(bf.knotslist))

tbank.knots = bf.knotslist

print("Setting padding flag array in PWFCommandLine")
# Creating appropriate padding flags lists if they are not provided
if padding_flags_bbox == [0]:
        padding_flags_bbox = [0] * tbank.s * len(bf.knotslist)
if padding_flags_int == [0]:
        padding_flags_int = [0] * tbank.s * len(bf.knotslist)

# Building the name of the directory where the SFTs will be saved. Directory name is based off the tbank object. In this way, any searches completed
# using the same tbank will have all SFTs in the same directory for neatness. Within the tbankdirectory, sub folders are used to contain the SFTs
# and results of each different search
tbankstr = tbank.toString()
splitstr = tbankstr.split()

tbankdirectory = str(tbank.s) + "_"

for elem in splitstr[1:]:
        tbankdirectory = tbankdirectory + str("{:.3f}".format(float(elem[:-1] ))) + "_"

tbankdirectory = tbankdirectory[:-1]

if not os.path.isdir(tbankdirectory):
        os.mkdir(tbankdirectory)

# The name of the file that 2F and templates are stored to is only changed by whether the 2F values are averaged or only summed
sumstr = ""
if rtnsum:
        sumstr = "_Sum"

tempsfilename = "Temps" + sumstr + ".txt"
mismatchtempname = "Temps" + sumstr + "_Mismatch.txt"

# Build the SFTs and returns the random injected signal parameters. In the case that SFTFiles is false (the default) the SFTs are created in memory when
# the fstatinput method is called in PWFStat, however the randomly generate signal parameters are generated here

if SFTFiles:
        signalparams = pwsim.buildSFTs(tbank.dur, tbank, tstart=tstart, trefsegfrac=trefsegfrac, Tsft=finputdata.Tsft, parentdirectory=tbankdirectory, h0=h0, noiseratio=noiseratio)
else:
        print("Creating tiling lattice")
        tiling = lp.CreateLatticeTiling(tbank.s * len(bf.knotslist))

        print("Building metric")
        metric = scmm.PreCompMetric(tbank.s)

        print("Setting Bounds")
        tbe.setbounds(tiling, padding_flags_bbox, padding_flags_int, tbank, reset)

        print("Setting tiling lattice and metric")
        print(padding_flags_int)
        lp.SetTilingLatticeAndMetric(tiling, lp.TILING_LATTICE_ANSTAR, metric, tbank.mismatch)

        print("Creating random signal params")
        randparams = lal.CreateRandomParams(0)
        signalparams = lal.gsl_matrix(tbank.s * len(bf.knotslist), 1)
        lp.RandomLatticeTilingPoints(tiling, 0, randparams, signalparams);                                      # Chance 'scale' parameter to be non-zero (within range [-1, 0])

        signalparams = np.transpose(signalparams.data)[0]

print("Random Signal params are: " + str(signalparams))
logging.info("Random Signal Params are: %s", str(signalparams))

h0str = "{:.2E}".format(h0)
f0 = "_{:.3f}".format(signalparams[0]) + "_"

if tbank.s >= 2:
        f1 = "{:.2E}".format(signalparams[1]) +"_"
else:
        f1 = ""
if tbank.s == 3:
        f2 = "{:.2E}".format(signalparams[2]) + "_"
else:
        f2 = ""

# The directory where the SFTs will be saved to. This is the same directory name as in pwsim, so should not be altered unless also altered there (or SFTs are written
# in memory and not to files in which case it shouldn't matter if the name is altered)
directoryname = h0str + f0 + f1 + f2 + str(tbank.dur) + "_" + str(tstart)

if not os.path.isdir(tbankdirectory + "/" + directoryname):
                os.mkdir(tbankdirectory + "/" + directoryname)

# Create the text file where 2F and template data will be stored. The first line is added to the file which contains the column titles for all data.
with open(tbankdirectory + "/" + directoryname + "/" + tempsfilename, "w") as reader:

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
with open(tbankdirectory + "/" + directoryname + "/" + mismatchtempname, "w") as reader:

        commentline = "{:20s}".format("#Mismatch") + "     " + str("{:20s}".format("FStat")) + "     "

        knotnum = -1

        for i in range(len(signalparams)):

                derivnum = i % tbank.s

                if derivnum == 0:
                        knotnum += 1

                columntitle = "f_" + str(knotnum) + "_" + str(i % 3)
                commentline = commentline + " " + "{:^20s}".format(columntitle)

        reader.write(commentline + "\n")

print("Setting Antenna Pattern")
# Set the maximum condition number for the antenna pattern above its default (get a warning using the piecewise model if not manually changed to be larger)
lp.SetAntennaPatternMaxCond(10 ** 5)

# The minimum and maximum frequencies needed to load in from SFTs to cover all frequencies any template may cover
SFTfmin = gom.gte(tbank.dur, tbank.fmin, tbank.nmax, gom.kwhichresultsingivenhalflife(tbank.taumin, tbank.fmax, tbank.nmax)) - 10
SFTfmax = tbank.fmax + 50

if SFTfmin <= 0:
        SFTfmin = 1

# Build parameter space and begin calculated 2Fs
print("SFTfmin/fmax: " + str([SFTfmin, SFTfmax]))
print("SFT files is: " + str(SFTFiles))
pwf.semifstatcatalogue(SFTfmin, SFTfmax, tbank, finputdata, signalparams, tempsfilename, reset, padding_flags_bbox=padding_flags_bbox, padding_flags_int=padding_flags_int, directory=tbankdirectory + "/" + directoryname, trefsegfrac=trefsegfrac, rtnsum=rtnsum, SFTFiles=SFTFiles, tempsperfile=tempsperfile)

# For logging and profiling
pr.disable()
s = io.StringIO()
sortby = pstats.SortKey.CUMULATIVE
ps = pstats.Stats(pr, stream=s).sort_stats('cumtime')
ps.print_stats()

with open("PWFCommandLineProfile.txt", 'w+') as f:
        f.write(s.getvalue())
