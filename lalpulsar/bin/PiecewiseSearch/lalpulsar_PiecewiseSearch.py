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
import os
import pstats
import random
import signal
import sys
import time

import matplotlib.pyplot as plt
import numpy as np

import lal
import lalpulsar as lp
import lalpulsar.piecewise_model.BasisFunctions as bf
import lalpulsar.piecewise_model.ClassDefinitions as cd
import lalpulsar.piecewise_model.EstimatingKnots as ek
import lalpulsar.piecewise_model.GTEandOtherMethods as gom
import lalpulsar.piecewise_model.PWFStat as pwf
import lalpulsar.piecewise_model.PWModelSimulations as pwsim
import lalpulsar.piecewise_model.SemicoherentMetricMethods as scmm
import lalpulsar.piecewise_model.TBankEstimates as tbe

start_time = time.time()

job_num = os.getenv("SLURM_ARRAY_TASK_ID")
print("Job_array_num_is_" + str(job_num))

# Allow liblalpulsar C code to be interrupted with Ctrl+C
signal.signal(signal.SIGINT, signal.SIG_DFL)

parser = ap.ArgumentParser()

# Thoughts on this notebook
# It is beginning to get quite convoluted. Is it worth writing a new python module where a lot of this set up is written as functions? It may make this notebook more succinct
# and readable

# Required parameters
parser.add_argument(
    "--tbankcode",
    type=str,
    help="Specifies the tbank object to use for a search",
    required=True,
)

# Random seed, typically chosen as the job array number of OzSTAR. To ensure unique random seeds, j < 10,000
parser.add_argument(
    "--j",
    type=int,
    help="The job array number on OzSTAR, to ensure unique random seeds, j < 10,000",
    default=0,
)

# Parameters used for specifying signals to inject in SFTs as well as other parameters for F-Statistic calculations
parser.add_argument("--tstart", type=int, help="Start time of SFTs", default=None)
parser.add_argument(
    "--trefsegfrac",
    type=float,
    help="What percentage tref occurs on each segment. E.g. if 1/2, then tref will be halfway through each segment",
    default=0.0,
)
parser.add_argument(
    "--noise_sqrt_Sh",
    type=float,
    help="Noise level to be added to data",
    nargs="+",
    default=[],
)
parser.add_argument(
    "--tempsperfile",
    type=int,
    help="The number of templates with the greatest fstats to record",
    default=1000,
)
parser.add_argument(
    "--rtnsum",
    help="Whether to return the 2F stat averaged over all segments",
    action="store_true",
)
parser.add_argument(
    "--SFTFiles",
    help="Whether to write SFT files or have them stored in memory",
    action="store_true",
)
parser.add_argument(
    "--UseLocal",
    help="If selected all files will be written to local storage and deleted once job finished except return values of pwf.semifstatcatalogue",
    action="store_true",
)
parser.add_argument(
    "--BaseSeed",
    type=int,
    help="BaseSeed * 10e7 + job number * 10e3 + iterations is the seed for generating random signal params",
    default=0,
)

# Parameters to be used when creating the necessary elements to use lp.ComputeFStat
parser.add_argument(
    "--Tsft", type=int, help="The length of each SFT to be built", default=10
)
parser.add_argument("--h0", type=float, help="Strain of the injected signal", default=0)
parser.add_argument(
    "--cosi", type=float, help="Default injection cos-iota", default=None
)
parser.add_argument("--psi", type=float, help="Default injection psi", default=None)
parser.add_argument("--phi0", type=float, help="Default injection phi0", default=None)
parser.add_argument(
    "--dt_wf", type=float, help="Sampling time for waveform generation", default=10
)
parser.add_argument(
    "--Alpha",
    type=float,
    help="Right Ascension to sky position of source",
    default=None,
)
parser.add_argument(
    "--Delta", type=float, help="Declination to sky position of source", default=None
)
parser.add_argument(
    "--detectors",
    type=str,
    help="Set detector being used. Options H1, V1, etc.",
    nargs="+",
    default=[],
)
parser.add_argument(
    "--dfreq", type=float, help="Frequency spacing for ComputeFstat()", default=0
)
parser.add_argument(
    "--sourceDeltaT",
    type=float,
    help="sourceDeltaT option for ComputeFstat()",
    default=0.1,
)

# Parameters used in defining the parameter space
parser.add_argument(
    "--s", type=int, help="Set the s parameter", nargs="?", const=1, default=None
)
parser.add_argument(
    "--fmin",
    type=float,
    help="Set the fmin parameter",
    nargs="?",
    const=1,
    default=None,
)
parser.add_argument(
    "--fmax",
    type=float,
    help="Set the fmax parameter",
    nargs="?",
    const=1,
    default=None,
)
parser.add_argument(
    "--nmin",
    type=float,
    help="Set the nmin parameter",
    nargs="?",
    const=1,
    default=None,
)
parser.add_argument(
    "--nmax",
    type=float,
    help="Set the nmax parameter",
    nargs="?",
    const=1,
    default=None,
)
parser.add_argument(
    "--Izz",
    type=float,
    help="Inertia of the NS. Used to calculate default k values",
    nargs="?",
    const=1,
    default=None,
)
parser.add_argument(
    "--ellip",
    type=float,
    help="Ellipticity of NS. Used to calculate default k values",
    nargs="?",
    const=1,
    default=None,
)
parser.add_argument(
    "--radius",
    type=float,
    help="Radius of the  NS. Used to calculate default k values",
    nargs="?",
    const=1,
    default=None,
)
parser.add_argument(
    "--kmin",
    type=float,
    help="Set the kmin parameter. Overrides Izz, ellip and radius values",
    nargs="?",
    const=1,
    default=None,
)
parser.add_argument(
    "--kmax",
    type=float,
    help="Set the kmax parameter. Overrides Izz, ellip and radius values",
    nargs="?",
    const=1,
    default=None,
)
parser.add_argument(
    "--dur",
    type=float,
    help="The maximum duration the knot algorithm to calculate up to",
    nargs="?",
    const=1,
    default=None,
)
parser.add_argument(
    "--knots", type=float, help="Use user defined knots", nargs="+", default=[0]
)
parser.add_argument(
    "--maxmismatch",
    type=float,
    help="Set the mismatch parameter",
    nargs="?",
    const=1,
    default=None,
)
parser.add_argument(
    "--knotnum",
    type=int,
    help="Use a knotnum number of knots generated by algorithm",
    nargs="?",
    const=1,
    default=None,
)
parser.add_argument(
    "--knot_alg_dur",
    type=float,
    help="Uses knots generated by the algorithm extending to this duration",
    nargs="?",
    const=1,
    default=None,
)
parser.add_argument(
    "--flags_bbox",
    type=float,
    help="Fractional bounding box padding to use on each dimension",
    nargs="+",
    default=[],
)
parser.add_argument(
    "--flags_intbox",
    type=int,
    help="Integer points to use for padding on each dimension",
    nargs="+",
    default=[],
)

# Different modes for the code to run
parser.add_argument(
    "--det_prob",
    help="If selected, does an iterations number of searches, with signal h0 varying from min_h0 to max_h0. Also includes h0 = 0 for computing 2F* simultaneously",
    action="store_true",
)
parser.add_argument(
    "--threshold_2F",
    help="Does an iterations number of searches with no injected signal. Useful for calculating 2F* without calculating det_prob",
    action="store_true",
)
parser.add_argument(
    "--min_h0",
    type=float,
    help="Minimum h0 for calculating detection probabilities, det_prob must be true",
    default=1e-30,
)
parser.add_argument(
    "--max_h0",
    type=float,
    help="Maximum h0 for calculating detection probabilities, det_prob must be true",
    default=1e-10,
)
parser.add_argument(
    "--template_count",
    help="Returns the template count for the parameter space. Does not compute anything else",
    action="store_true",
)
parser.add_argument(
    "--mis_hist",
    help="If selected with Fstat_mismatch, only returns the lowest mis_match",
    action="store_true",
)
parser.add_argument(
    "--fake_data",
    help="If selected, carries out a search on fake data with an injected signal",
    action="store_true",
)
parser.add_argument(
    "--inject_data",
    help="If selected, injects a randomly generated signal into the SFTs located at the path set by . Can be selected alongside other modes",
    action="store_true",
)

# Other options available
parser.add_argument(
    "--Fstat_mismatch",
    help="If selected, determines mismatches using the 2F definition, otherwise uses the metric",
    action="store_true",
)
parser.add_argument(
    "--signal_as_template",
    help="If selected, the random injected signal will coincide with a template",
    action="store_true",
)
parser.add_argument(
    "--iterations",
    type=int,
    help="Number of loops for the code to perform. For unique random seeds, iterations < 1000",
    default=1,
)
parser.add_argument(
    "--freq_bands",
    type=float,
    help="Defines frequency bands for search to cover. Frequency bands chosen as fmin=freq_bands[j], fmax=freq_bands[j+1]",
    nargs="+",
    default=[],
)
parser.add_argument(
    "--fstat_hist",
    help="If selected, returns data on the distribution of 2Fs as a histogram",
    action="store_true",
)
parser.add_argument(
    "--data_psd",
    help="Use simulated data with a noise level equal to the PSD of data",
    action="store_true",
)

# Parameters controlling outputs
parser.add_argument("--outbasedir", type=str, help="Output base directory", default=".")
parser.add_argument(
    "--logfile",
    help="If true, log to file; if false, log to standard output",
    action="store_true",
)
parser.add_argument(
    "--loglevel", help="Level at which to log messages", default=logging.INFO
)
parser.add_argument("--profile", help="Profile the code", action="store_true")
parser.add_argument(
    "--noise_path",
    type=str,
    help="Path to noise data. Default option harded coded in tbank object",
    default=None,
)
parser.add_argument(
    "--sft_path",
    type=str,
    help="Path to sfts. Default option hard coded in tbank object",
    default=None,
)

# Hacky way to know the absolute path to the noise curve data by specifying which machine you are on as either 'PC' or 'OzSTAR'
parser.add_argument(
    "--machine",
    type=str,
    help="Specify the machine you are working on to find the path to noise curves",
    default="OzSTAR",
)

# -------------------------------------------------------------------------------------------------------------------------------------------------------
# Logging, profiling and setting the directory for the knots list
# -------------------------------------------------------------------------------------------------------------------------------------------------------

args = parser.parse_args()
outbasedirectory = args.outbasedir

ek.setknotarchivepath(outbasedirectory)

# For logging
logging_format = "%(asctime)s : %(levelname)s : %(message)s"
if args.logfile:
    logging_filename = os.path.join(outbasedirectory, f"PiecewiseSearchLog_{j}.log")
    logging.basicConfig(
        filename=logging_filename,
        filemode="w",
        level=args.loglevel,
        format=logging_format,
    )
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

# -------------------------------------------------------------------------------------------------------------------------------------------------------
# Required arguments
# -------------------------------------------------------------------------------------------------------------------------------------------------------

tbankcode = args.tbankcode

# -------------------------------------------------------------------------------------------------------------------------------------------------------
# Optional arguments and Random Seed
# -------------------------------------------------------------------------------------------------------------------------------------------------------

j = args.j
trefsegfrac = args.trefsegfrac
rtnsum = args.rtnsum
tempsperfile = args.tempsperfile
SFTFiles = args.SFTFiles
UseLocal = args.UseLocal
BaseSeed = args.BaseSeed
iterations = args.iterations

Random_Seed = BaseSeed * 10000000 + j * 1000
Random_Seed_Str = "{:0>9d}".format(Random_Seed)

# -------------------------------------------------------------------------------------------------------------------------------------------------------
# Setting the mode the code is run in, influencing its output
# -------------------------------------------------------------------------------------------------------------------------------------------------------
mode = "search"

if args.det_prob:
    mode = "det_prob"
elif args.threshold_2F:
    mode = "threshold_2F"
elif args.template_count:
    mode = "template_count"
elif args.mis_hist:
    mode = "mis_hist"
elif args.fake_data:
    mode = "fake_data"

logging.info(f"Search mode: {mode}")

# -------------------------------------------------------------------------------------------------------------------------------------------------------
# Setting the template bank object and changing user specified parameters
# -------------------------------------------------------------------------------------------------------------------------------------------------------
tbank = cd.TBank()

if tbankcode == "GW170817":
    logging.info(f"Starting with default parameter space: {tbankcode}")
    tbank.SetDefaultGW170817()

elif tbankcode == "GW190425":
    logging.info(f"Starting with default parameter space: {tbankcode}")
    tbank.SetDefaultGW190425()

elif tbankcode == "1987A":
    logging.info(f"Starting with default parameter space: {tbankcode}")
    tbank.SetDefault1987A()

elif tbankcode == "CasA":
    logging.info(f"Starting with default parameter space: {tbankcode}")
    tbank.SetDefaultCasA()
else:
    print("No tbank selected! Put an error here!")

# If using custom frequency bands
if args.freq_bands:
    args.fmin = args.freq_bands[j]
    args.fmax = args.freq_bands[j + 1]

# Override any tbank arguments set on the command line
tbank.SetTBankParams(args)

# -------------------------------------------------------------------------------------------------------------------------------------------------------
# Checking for user defined knots, or if knots need to be recalculated based on any changed default parameters
# -------------------------------------------------------------------------------------------------------------------------------------------------------
if args.knots != [0]:
    logging.info("Setting knots from command line")
    bf.knotslist = args.knots
elif args.knot_alg_dur:
    logging.info("Setting knots from algorithm")
    tbank.dur = args.knot_alg_dur
    tbank.SetKnotsByAlg()
else:
    logging.info("Using default knots")
    bf.knotslist = tbank.knots

tbank.knots = bf.knotslist
tbank.dur = tbank.knots[-1] - tbank.knots[0]

# Adjusting knot start time and duration
if bf.knotslist[0] == 0:
    for i, knot in enumerate(bf.knotslist):
        bf.knotslist[i] = knot + tbank.tstart

logging.info("Knots: %s", str(bf.knotslist))

# -------------------------------------------------------------------------------------------------------------------------------------------------------
# Padding flag customisation
# -------------------------------------------------------------------------------------------------------------------------------------------------------
if args.flags_bbox == []:
    tbank.flags_bbox = [-1] * tbank.s * len(bf.knotslist)

if args.flags_intbox == []:
    tbank.flags_intbox = [-1] * tbank.s * len(bf.knotslist)

# -------------------------------------------------------------------------------------------------------------------------------------------------------
# Setting detector specifics for noise curve data
# -------------------------------------------------------------------------------------------------------------------------------------------------------

detectors = tbank.detectors

logging.info(f"Detectors: {detectors}")
logging.info(f"Arg Detectors: {args.detectors}")
logging.info(f"TBank Detectors: {tbank.detectors}")

num_det = len(detectors)

detector_str_list = lal.CreateStringVector(detectors[0])

# I'm sure there's a better way to initialise the detector LALStringVector
for detector in detectors[1:]:
    lal.AppendString2Vector(detector_str_list, detector)

machine = args.machine

if machine == "PC":
    noise_path = tbank.pc_noise_path
    psd_path = tbank.pc_psd_path
elif machine == "OzSTAR":
    noise_path = tbank.ozstar_noise_path
    psd_path = tbank.ozstar_psd_path

noise_curve_data = []

if args.noise_sqrt_Sh == [] and mode != "search":

    for detector in detectors:
        path_to_data = noise_path + detector + "asd.txt"

        if args.data_psd:
            path_to_data = psd_path + detector + ".txt"

        min_freq_strain = -1
        max_freq_strain = -1

        if mode == "template_count":
            break

        with open(path_to_data) as noise_curve:
            for line in noise_curve:

                this_line = line.split()
                frequency = float(this_line[0])
                strain = float(this_line[1])

                if min_freq_strain == -1 and frequency > args.fmin:
                    min_freq_strain = strain

                if max_freq_strain == -1 and frequency > args.fmax:
                    max_freq_strain = strain
                    break

        this_det_noise = (min_freq_strain + max_freq_strain) / 2
        noise_curve_data.append(this_det_noise)

logging.info("Noise curve data for given frequency band: %s", str(noise_curve_data))
logging.info("User specified noise used (if given): %s", str(args.noise_sqrt_Sh))

# -------------------------------------------------------------------------------------------------------------------------------------------------------
# Setting the FInput object
# -------------------------------------------------------------------------------------------------------------------------------------------------------

# Note that the parameters cosi, psi and phi0 are randomised in a for loop below.

finputdata = cd.FInput()
finputdata.Tsft = args.Tsft
finputdata.tstart = tbank.tstart
finputdata.trefsegfrac = args.trefsegfrac
finputdata.h0 = args.h0
finputdata.cosi = args.cosi
finputdata.psi = args.psi
finputdata.phi0 = args.phi0
finputdata.dt_wf = args.dt_wf
finputdata.Alpha = tbank.Alpha
finputdata.Delta = tbank.Delta
finputdata.detectors = detector_str_list
finputdata.noise_sqrt_Sh = (
    args.noise_sqrt_Sh if args.noise_sqrt_Sh != [] else noise_curve_data
)
finputdata.dfreq = args.dfreq
finputdata.sourceDeltaT = args.sourceDeltaT
finputdata.inject_data = args.inject_data

# -------------------------------------------------------------------------------------------------------------------------------------------------------
# Building the name of the directory where the SFTs will be saved. Directory name is based off the tbank object. In this way, any searches completed
# using the same tbank will have all SFTs in the same directory for neatness. Within the tbankdirectory, sub folders are used to contain the SFTs
# and results of each different search. If UseLocal, then all files are written to OzSTAR local disks, which are deleted when jobs are finished.
# -------------------------------------------------------------------------------------------------------------------------------------------------------

append_string = "_BS-" + str(args.BaseSeed) + "_j-" + str(j)

if UseLocal:
    local_disk_directory = os.getenv("JOBFS")
    tbankdirectory = (
        os.path.join(local_disk_directory, tbank.toString()) + append_string
    )
else:
    tbankdirectory = os.path.join(outbasedirectory, tbank.toString()) + append_string

if not os.path.isdir(tbankdirectory):
    os.mkdir(tbankdirectory)

logging.info(f"TBank directory is {tbankdirectory}")

# -------------------------------------------------------------------------------------------------------------------------------------------------------
# Building strain list for injected signals, only necessary if --det_prob option is chosen
# -------------------------------------------------------------------------------------------------------------------------------------------------------
minimum_h0 = args.min_h0
maximum_h0 = args.max_h0

power_list = np.linspace(np.log10(minimum_h0), np.log10(maximum_h0), iterations)

h0_list = [0] + [10**i for i in power_list]

plus_one = 0
if args.det_prob:
    plus_one = 1

# -------------------------------------------------------------------------------------------------------------------------------------------------------
# Begin search set up. Some options may benefit from running through multiple independent searches, hence all search setup is done in a for loop
# -------------------------------------------------------------------------------------------------------------------------------------------------------
for i in range(iterations + plus_one):
    logging.info("Performing iteration: %s", str(i))

    # Change strain
    if args.det_prob:
        finputdata.h0 = h0_list[i]

    if args.threshold_2F:
        finputdata.h0 = 0

    logging.debug("Signal strength: %s", str(finputdata.h0))

    # Random seed for this iteration
    this_seed = Random_Seed + i

    # The name of the file that 2F and templates are stored to is only changed by whether the 2F values are averaged or only summed
    sumstr = ""
    if rtnsum:
        sumstr = "_Sum"

    iteration_string = "_i-" + str(i)

    tempsfilename = "Temps_For_" + tbank.toString() + iteration_string + sumstr + ".txt"
    mismatchtempname = (
        "Mismatch_Temps_For_" + tbank.toString() + iteration_string + sumstr + ".txt"
    )

    # Randomise the cosi, phi0 and psi parameters if we are not doing a targetted search
    if mode != "search" or args.inject_data:

        seed = random.seed(this_seed)

        finputdata.phi0 = random.uniform(0, 2 * np.pi)
        finputdata.psi = random.uniform(0, 2 * np.pi)
        finputdata.cosi = random.uniform(-1, 1)

        logging.info("phi0: %s", str(finputdata.phi0))
        logging.info("psi: %s", str(finputdata.psi))
        logging.info("cosi: %s", str(finputdata.cosi))

    # This chain of else-if statements determines the path to the SFTs being used and the signal parameters that will be/have been injected into them. The first option
    # is if we are generating our own fake data with gaussian noise by building SFT files. Using the SFT files is advantageous as we can inject any type of signal we wish,
    # depending on what signal model is defined in pwsim. Signal parameters for the first option are randomly generated inside pwsim. The second option randomly generates
    # signal parameters inside the parameter space. This option is used if we do not want to write SFT files, but rather store SFTs in memory. This option should be used for
    # all configurations which are not a search on experimental data OR we are doing a search with experimental data but we want to inject a signal into that data. If the case
    # is the latter, the random signal parameters are generated, but an enclosed if statement resets the SFT directory to the appropriate path for the experimental data.
    # The final option is for if we are doing a search, in which case we use the SFT directory to the data SFTs and set the injected signal to all zeros.

    if SFTFiles:
        signalparams, SFTdirectory = pwsim.buildSFTs(
            tbank.dur,
            tbank,
            finputdata,
            trefsegfrac=trefsegfrac,
            parentdirectory=tbankdirectory,
            rand_seed=this_seed,
        )

    elif mode != "search" or (mode == "search" and args.inject_data):
        logging.debug("Creating tiling lattice")
        tiling = lp.CreateLatticeTiling(tbank.s * len(bf.knotslist))

        logging.debug("Building metric")
        metric = scmm.PreCompMetric(tbank.s)

        logging.debug("Setting Bounds")

        tbe.setbounds(tiling, tbank)

        logging.debug("Setting tiling lattice and metric")
        lp.SetTilingLatticeAndMetric(
            tiling, lp.TILING_LATTICE_ANSTAR, metric, tbank.maxmismatch
        )

        logging.debug("Creating random signal params")

        randparams = lal.CreateRandomParams(this_seed)
        signalparams = lal.gsl_matrix(tbank.s * len(bf.knotslist), 1)
        lp.RandomLatticeTilingPoints(tiling, 0, randparams, signalparams)

        # -------------------------------------------------------------------------------------------------------------------------------------------------------
        # If we want a signal to coincide exactly with a template from the template bank. A random point in the parameter space is generated, then the nearest
        # template found to that point is chosen to be the injected signal
        # -------------------------------------------------------------------------------------------------------------------------------------------------------
        if args.signal_as_template:
            locator = lp.CreateLatticeTilingLocator(tiling)

            signal_params_vec = lal.gsl_vector(tbank.s * len(bf.knotslist))

            # Done as a for loop because setting signal_params_vec.data = signalparams.data[0] doesn't seem to work
            for i, elem in enumerate(signalparams.data[0]):
                signal_params_vec.data[i] = elem

            nearest_point = lal.gsl_vector(tbank.s * len(bf.knotslist))
            nearest_index = lal.CreateUINT8Vector(tbank.s * len(bf.knotslist))

            lp.NearestLatticeTilingPoint(
                locator, signal_params_vec, nearest_point, nearest_index
            )

            signalparams = np.transpose(nearest_point.data)

        else:
            signalparams = np.transpose(signalparams.data)[0]

        logging.info("Random Signal Params are: %s", str(signalparams))

        # Details for the SFT directory

        f0 = signalparams[0]
        f1 = signalparams[1]
        f2 = signalparams[2]

        SFTdirectory = os.path.join(
            tbankdirectory,
            f"SFTs_h0-{finputdata.h0:.2e}_f0-{f0:.3f}_f1-{f1:.2e}_f2-{f2:.2e}_dur-{tbank.dur}_tstart-{tbank.tstart}/",
        )

        if finputdata.inject_data:

            if machine == "PC":
                SFTdirectory = tbank.pc_sft_path
            elif machine == "OzSTAR":
                SFTdirectory = tbank.ozstar_sft_path

    elif mode == "search" and not args.inject_data:

        if machine == "PC":
            SFTdirectory = tbank.pc_sft_path
        elif machine == "OzSTAR":
            SFTdirectory = tbank.ozstar_sft_path

        signalparams = [-1] * (tbank.s * len(bf.knotslist))

    logging.info("SFT's for this iteration have path: %s", SFTdirectory)

    # Create the text file where 2F and template data will be stored. The first line is added to the file which contains the column titles for all data. We also write
    # the signal parameters as the first data entry to this file, however we skip that step if we are carrying out a search
    with open(os.path.join(tbankdirectory, tempsfilename), "w") as reader:

        commentline = "{:20s}".format("#FStat") + "     "

        for detector_string in detectors:
            commentline += "{:20s}".format(detector_string) + "     "

        commentline += str("{:20s}".format("Mismatch")) + "     "

        knotnum = -1

        for i in range(len(signalparams)):

            if mode == "search":
                break

            derivnum = i % tbank.s

            if derivnum == 0:
                knotnum += 1

            columntitle = "f_" + str(knotnum) + "_" + str(derivnum)
            commentline += " " + "{:^20s}".format(columntitle)

        reader.write(commentline + "\n")

    # Create the text file where 2F and template data will be stored but ordered by mismatch. The first line is added which contains the column titles for all data. We skip
    # this if we are doing a search as we can't calculate the mismatch
    if mode != "search":
        with open(os.path.join(tbankdirectory, mismatchtempname), "w") as reader:

            commentline = (
                "{:20s}".format("#Mismatch")
                + "     "
                + str("{:20s}".format("FStat"))
                + "     "
            )

            for detector_string in detectors:
                commentline += "{:20s}".format(detector_string) + "     "

            knotnum = -1

            for i in range(len(signalparams)):

                derivnum = i % tbank.s

                if derivnum == 0:
                    knotnum += 1

                columntitle = "f_" + str(knotnum) + "_" + str(derivnum)
                commentline += " " + "{:^20s}".format(columntitle)

            reader.write(commentline + "\n")

    # Set the maximum condition number for the antenna pattern above its default (get a warning using the piecewise model if not manually changed to be larger)
    logging.debug("Setting Antenna Pattern")
    lp.SetAntennaPatternMaxCond(10**5)

    # The minimum and maximum frequencies needed to load in from SFTs to cover all frequencies any template may cover
    SFTfpad = 10  # 1800 / tbank.Tsft + tbank.dur / 86400 + 5
    SFTfmin = gom.gte(tbank.dur, tbank.fmin, tbank.nmax, tbank.kmax) - SFTfpad - 10
    SFTfmax = gom.gte(0, tbank.fmax, tbank.nmin, tbank.kmin) + SFTfpad

    logging.info(f"SFTfmin/fmax: [{SFTfmin}, {SFTfmax}]")
    logging.info(f"SFT files is: {SFTFiles}")

    Fstat_mismatch = args.Fstat_mismatch
    logging.info(f"2F mismatch is: {Fstat_mismatch}")

    # Output from the search, which is carried out in PWFStat.py. Output changes depending on run mode. All search results are written to files within this method
    if mode == "search" and not args.fstat_hist:
        (
            lowest_mismatch_metric,
            lowest_mismatch_Fstat,
            fstat_max,
            template_count,
        ) = pwf.pw_fstat_search_catalogue(
            SFTfmin,
            SFTfmax,
            tbank,
            finputdata,
            tempsfilename,
            tbankdirectory=tbankdirectory,
            SFTdirectory=SFTdirectory,
            trefsegfrac=trefsegfrac,
            rtnsum=rtnsum,
            tempsperfile=tempsperfile,
        )
    else:
        (
            lowest_mismatch_metric,
            lowest_mismatch_Fstat,
            fstat_max,
            template_count,
        ) = pwf.semifstatcatalogue(
            SFTfmin,
            SFTfmax,
            tbank,
            finputdata,
            signalparams,
            tempsfilename,
            out_directory=tbankdirectory,
            SFTdirectory=SFTdirectory,
            trefsegfrac=trefsegfrac,
            rtnsum=rtnsum,
            SFTFiles=SFTFiles,
            tempsperfile=tempsperfile,
            Fstat_mismatch=args.Fstat_mismatch,
            fstat_hist=args.fstat_hist,
            mode=mode,
        )

    # For each non-zero return value above, we create a directory and write the data into a file based on Random_Seed (not the same random seed used inside the four loop)

    # Save the template count of the search and the time taken to conduct the search (if template_count option selected, the time taken is not indicative of
    # how long the search takes)
    if template_count != 0:
        temp_count_directory = "Template_Count_And_Timing_For_" + tbank.toString()
        temp_count_file = "Template_Count_And_Timing_" + Random_Seed_Str + ".txt"

        if not os.path.isdir(temp_count_directory):
            os.mkdir(temp_count_directory)

        # Writing the total time taken and total templates counted
        time_taken = time.time() - start_time

        with open(temp_count_directory + "/" + temp_count_file, "a+") as reader:
            line = str(time_taken) + "        " + str(template_count) + "\n"

            reader.write(line)

    # If we are in search mode, the lowest mismatches aren't calculated, so we can skip this step. Mismatches written in two columns. First column is metric mismatches, second is Fstat mismatches.
    # We also don't save the largest 2Fs in search mode, as they should be saved with our templates separately.
    if mode != "search":

        # Save the highest found fstat
        if fstat_max != 0:
            Two_F_directory = "Largest_2Fs_For_" + tbank.toString()
            Two_F_file = "Largest_2F_" + Random_Seed_Str + ".txt"

            if not os.path.isdir(Two_F_directory):
                os.mkdir(Two_F_directory)

            with open(Two_F_directory + "/" + Two_F_file, "a+") as reader:
                line = str(finputdata.h0) + "        " + str(fstat_max) + "\n"

                reader.write(line)

        # Save the lowest mismatches
        if lowest_mismatch_metric != 0:
            mismatch_directory = "Mismatches_For_" + tbank.toString()
            mismatch_file = "Mismatch_" + Random_Seed_Str + ".txt"

            if not os.path.isdir(mismatch_directory):
                os.mkdir(mismatch_directory)

            with open(mismatch_directory + "/" + mismatch_file, "a+") as reader:
                line = (
                    str(lowest_mismatch_metric)
                    + "        "
                    + str(lowest_mismatch_Fstat)
                    + "\n"
                )
                reader.write(line)

logging.info("Done")

# For logging and profiling
if args.profile:
    pr.disable()
    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats("cumtime")
    ps.print_stats()
    with open(
        os.path.join(basedirectory, f"PiecewiseSearchProfile_{j}.txt"), "w+"
    ) as f:
        f.write(s.getvalue())
