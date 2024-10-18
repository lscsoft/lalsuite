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

## \file
## \ingroup lalpulsar_python_piecewise_model
"""
Compute the F-statistic for a piecewise model template.
"""

import copy
import heapq as hq
import logging
import time

import numpy as np
from tqdm import trange

import lal
import lalpulsar as lp

from . import basis_functions as bf
from . import gte_and_other_methods as gom
from . import semicoherent_metric_methods as scmm
from . import tbank_estimates as tbe

# Returns our FstatInput object, currently uses SFTs built in pwsim
# fmax (and fmin) should be chosen with the same reasoning outlined in the pwsim notbook. fmax = f0 + 10^-4 * f0 + 58/TSFT,
# fmin is the same except fmin = f0 - .... AND fmin should be lower than the minimum value any given template we might use
# over the time period of the relevant SFTs!
def fstatinputfromSFTFiles(
    SFTfmin,
    SFTfmax,
    finputdata,
    timeint=[-1, -2],
    sft_path=".",
    detectors=[],
    inject_params=[],
):

    Fstat_opt_args = lp.FstatOptionalArgs(lp.FstatOptionalArgsDefaults)
    Fstat_opt_args.FstatMethod = lp.FMETHOD_DEMOD_BEST
    Fstat_opt_args.sourceDeltaT = finputdata.sourceDeltaT

    logging.info(f"SFT path: {sft_path}")

    all_minus_one = True

    for elem in inject_params:
        if elem != -1:
            all_minus_one = False
            break
    logging.debug(inject_params)
    logging.debug(all_minus_one)
    logging.debug(not inject_params == [] and not all_minus_one)

    # This code injects a signal, if inject_params is non-empty and all elements are not -1, then we inject the signal present in inject_params
    if inject_params != [] and not all_minus_one:
        Tdata = timeint[-1] - timeint[0]

        h0 = finputdata.h0
        cosi = finputdata.cosi
        psi = finputdata.psi
        phi0 = finputdata.phi0
        tref = timeint[0] + Tdata * finputdata.trefsegfrac
        Alpha = finputdata.Alpha
        Delta = finputdata.Delta

        # create injection parameters
        Fstat_signal = lp.CreatePulsarParamsVector(1)
        Fstat_signal.data[0].Amp.aPlus = 0.5 * h0 * (1.0 + cosi * cosi)
        Fstat_signal.data[0].Amp.aCross = h0 * cosi
        Fstat_signal.data[0].Amp.psi = psi
        Fstat_signal.data[0].Amp.phi0 = phi0
        Fstat_signal.data[0].Doppler.refTime = tref
        Fstat_signal.data[0].Doppler.Alpha = Alpha
        Fstat_signal.data[0].Doppler.Delta = Delta

        logging.info(f"Injecting parameters: {inject_params} with h0: {h0}")

        # Set doppler params. Done as a for loop in case inject_params changes length (it shouldn't, but just in case)
        for i in range(len(inject_params)):
            Fstat_signal.data[0].Doppler.fkdot[i] = inject_params[i]

        Fstat_opt_args.injectSources = Fstat_signal

    dfreq = finputdata.dfreq
    emphdata = lp.InitBarycenter("earth00-40-DE405.dat.gz", "sun00-40-DE405.dat.gz")

    constraints = lp.SFTConstraints()

    constraints.minStartTime = lal.LIGOTimeGPS(timeint[0])
    constraints.maxStartTime = lal.LIGOTimeGPS(timeint[1])

    if len(detectors) == 0:
        catalogue = lp.SFTdataFind(sft_path + "*.sft", constraints)

        sft_path_master = sft_path + "*.sft"

    elif len(detectors) == 1:

        det1_path = sft_path + detectors[0] + "/*.sft"

        catalogue = lp.SFTdataFind(det1_path, constraints)

        sft_path_master = det1_path

    elif len(detectors) == 2:

        det1_path = sft_path + detectors[0] + "/*.sft"
        det2_path = sft_path + detectors[1] + "/*.sft"

        sft_path_master = det1_path + ";" + det2_path

        catalogue = lp.SFTdataFind(sft_path_master, constraints)

    elif len(detectors) == 3:

        det1_path = sft_path + detectors[0] + "/*.sft"
        det2_path = sft_path + detectors[1] + "/*.sft"
        det2_path = sft_path + detectors[3] + "/*.sft"

        sft_path_master = det1_path + ";" + det2_path + ";" + det3_path

        catalogue = lp.SFTdataFind(sft_path_master, constraints)

    logging.debug("Heya!")

    logging.info(f"Path to SFT's being used: {sft_path_master}")

    fstat = lp.CreateFstatInput(
        catalogue, SFTfmin, SFTfmax, dfreq, emphdata, Fstat_opt_args
    )

    return fstat


# Create fake sft catalog
def SFTCatalog(tstart, Tdata, finputdata):

    num_det = finputdata.detectors.length
    Tsft = finputdata.Tsft

    sft_ts = lp.MakeMultiTimestamps(tstart, Tdata, Tsft, 0, num_det)
    sft_catalog = lp.MultiAddToFakeSFTCatalog(None, finputdata.detectors, sft_ts)

    return sft_catalog


# Builds SFTs in memory rather than files
def fstatinput(
    SFTfmin,
    SFTfmax,
    sft_catalog,
    finputdata,
    inject_params,
    timeint=[-1, -2],
    trefsegfrac=0.0,
):

    Tdata = timeint[1] - timeint[0]
    dfreq = finputdata.dfreq
    # create default F-statistic optional arguments
    Fstat_opt_args = lp.FstatOptionalArgs(lp.FstatOptionalArgsDefaults)
    Fstat_opt_args.FstatMethod = lp.FMETHOD_DEMOD_BEST
    Fstat_opt_args.sourceDeltaT = finputdata.sourceDeltaT

    h0 = finputdata.h0
    cosi = finputdata.cosi
    psi = finputdata.psi
    phi0 = finputdata.phi0
    tref = timeint[0] + Tdata * trefsegfrac
    Alpha = finputdata.Alpha
    Delta = finputdata.Delta

    # create injection parameters
    Fstat_signal = lp.CreatePulsarParamsVector(1)
    Fstat_signal.data[0].Amp.aPlus = 0.5 * h0 * (1.0 + cosi * cosi)
    Fstat_signal.data[0].Amp.aCross = h0 * cosi
    Fstat_signal.data[0].Amp.psi = psi
    Fstat_signal.data[0].Amp.phi0 = phi0
    Fstat_signal.data[0].Doppler.refTime = tref
    Fstat_signal.data[0].Doppler.Alpha = Alpha
    Fstat_signal.data[0].Doppler.Delta = Delta

    logging.debug(f"Injecting parameters: {inject_params} with h0: {h0}")

    # Set doppler params. Done as a for loop in case inject_params changes length (it shouldn't, but just in case)
    for i in range(len(inject_params)):
        Fstat_signal.data[0].Doppler.fkdot[i] = inject_params[i]

    Fstat_opt_args.injectSources = Fstat_signal

    Fstat_injection_noise = lp.MultiNoiseFloor()
    Fstat_injection_noise.length = finputdata.detectors.length
    # Fstat_injection_noise.sqrtSn = finputdata.noise_sqrt_Sh

    for i, noise in enumerate(finputdata.noise_sqrt_Sh):
        Fstat_injection_noise.sqrtSn[i] = noise

    Fstat_opt_args.injectSqrtSX = Fstat_injection_noise

    ephemerides = lp.InitBarycenter("earth00-19-DE405.dat.gz", "sun00-19-DE405.dat.gz")

    # create F-statistic input data
    Fstat_input = lp.CreateFstatInput(
        sft_catalog, SFTfmin, SFTfmax, dfreq, ephemerides, Fstat_opt_args
    )

    return Fstat_input


# Returns our fstat results for a given template (using sft's from fstatin).
def computefstat(
    template,
    finputdata,
    fstatinputarray,
    transformmatrices,
    trefsegfrac=0.0,
    rtnsum=False,
):

    s = int(len(transformmatrices[0]) / 2)

    segs = int(len(template) / s - 1)

    dopparams = lp.PulsarDopplerParams()
    dopparams.Alpha = finputdata.Alpha
    dopparams.Delta = finputdata.Delta

    counter = 0

    fstat_res_array = []

    for seg in range(segs):
        counter += 1
        segtemp = template[seg * s : (seg + 2) * s]

        tref = (
            bf.knotslist[seg]
            + (bf.knotslist[seg + 1] - bf.knotslist[seg]) * trefsegfrac
        )
        dopparams.refTime = tref

        p0 = bf.knotslist[seg]
        p1 = bf.knotslist[seg + 1]

        fstatin = fstatinputarray[seg]

        tstart = bf.knotslist[0]
        dopplerparams = np.matmul(transformmatrices[seg], segtemp)

        # If you want to use conditioning and conditioning matrices, the below line is what should be used. The addition of the appropriate conditioning matrix will
        # be required also. The conditioning matrix that should be used is also given (for appropriate choices of p0, p1 and tref).
        # doplerparams = gom.PWParamstoTrefParamsPreComputed(segtemp, transformmatrices[seg], condmatrices[seg])
        # condmatrix = gom.ParamTransformationConditioningMatrix(p0, p1, tref)

        dopparams.fkdot[0 : len(dopplerparams)] = dopplerparams

        numFreqBins = 1
        whatToCompute = lp.FSTATQ_2F_PER_DET + lp.FSTATQ_2F
        fstatresults = lp.FstatResults()
        fstatresults = lp.ComputeFstat(
            fstatresults, fstatin, dopparams, numFreqBins, whatToCompute
        )

        fstat_res_array.append(fstatresults)

    return fstat_res_array


def fstat_res_heap_to_file(filename, fstat_res_heap, fstats_first=True):
    with open(filename, "a+") as reader:

        for elem in fstat_res_heap:
            if fstats_first:
                mismatch = elem[1]
            else:
                mismatch = elem[0]

            template = elem[3]
            fstat_res_array = elem[4]

            fstats_on_segs = [fstat_res.twoF[0] for fstat_res in fstat_res_array]
            fstat = sum(fstats_on_segs) / len(fstats_on_segs)

            # 2D list. Each sublist corresponds to a given detector. Each element of the sublists are the 2F's on each piecewise segment
            twoF_on_segs_on_detectors = []

            numdets = fstat_res_array[0].numDetectors

            for i in range(numdets):
                two_F_on_segs = [
                    fstat_res.twoFPerDet(i)[0] for fstat_res in fstat_res_array
                ]
                twoF_on_segs_on_detectors.append(two_F_on_segs)

            two_F_per_det = []

            for detector_2Fs in twoF_on_segs_on_detectors:
                this_det_2F = sum(detector_2Fs) / len(detector_2Fs)

                two_F_per_det.append(this_det_2F)

            # Note, the addition of each templateline string starting with a space makes that string have a length of 16, not 15
            templateline = " ".join("{:20.12E}".format(x) for x in template)

            fstat_line = "{:20.10f}".format(fstat) + "     "

            for twoF in two_F_per_det:
                fstat_line += "{:20.10f}".format(twoF) + "     "

            mismatch_line = "{:20.10f}".format(mismatch) + "     "

            if fstats_first:
                line = fstat_line + mismatch_line + templateline + "\n"
            else:
                line = mismatch_line + fstat_line + templateline + "\n"

            reader.write(line)


def add_fstat_to_dic(fstat, dic):

    base = 0

    if fstat < 50:
        base = 1
    elif 50 <= fstat < 500:
        base = 10
    else:
        base = 100

    lower_round = int(fstat - fstat % base)

    if lower_round in dic:
        dic[lower_round] += 1
    else:
        dic[lower_round] = 1


def dic_to_file(filename, dic):
    with open(filename, "a+") as reader:

        for k, v in dic.items():
            line = str(k) + ", " + str(v) + "\n"

            reader.write(line)


# For a template bank associated with the parameter space associated with tbank, finds the 'tbank.maxtemps' number of templates with the greatest FStat.
# SFTfmin and SFTfmax are the frequency range of SFTs we want to load in, tbank is the parameter space parameters we use. filename is the name of the file
# where the F-statistics and their corresponding templates will be saved. tbank.maxtemps is the top 2F's and templates to store and saved. directory is the
# path to where filename will be saved. trefsegfrac is how far along each segment the reference time will be defined. E.g. If set to 0, tref will be at the
# start of each segment, if trefsegfrac = 0.5, then tref will occur half way through each segment and so on. signalparams are the signal parameters of the signal
# injected into the SFTs and rtnsum tells us whether to save the summed 2F across all segments or the average. rtnsum = False means the average is saved,
# True and the sum will be saved
def semifstatcatalogue(
    SFTfmin,
    SFTfmax,
    tbank,
    finputdata,
    signalparams,
    filename,
    out_directory=".",
    SFTdirectory=".",
    trefsegfrac=0.0,
    rtnsum=False,
    SFTFiles=False,
    tempsperfile=1000,
    Fstat_mismatch=False,
    fstat_hist=False,
    mode="search",
):

    segs = len(bf.knotslist) - 1
    logging.info("Knots: %s", str(bf.knotslist))

    finalknot = len(bf.knotslist)
    s = tbank.s

    # The below lines are required for finding the nearest lattice point to the signal parameters

    # Create LatticeTiling object
    tiling = lp.CreateLatticeTiling(s * finalknot)

    logging.debug("Doing metric")
    metric = scmm.PreCompMetric(s)

    logging.debug("Doing Bounds")
    tbe.setbounds(tiling, tbank)

    logging.debug("Setting tiling lattice and metric")

    # Set metric, mismatch and lattice type
    lp.SetTilingLatticeAndMetric(
        tiling, lp.TILING_LATTICE_ANSTAR, metric, tbank.maxmismatch
    )

    signalpoint = lal.gsl_vector(s * finalknot)
    signalpoint.data = signalparams

    nearestpoint = lal.gsl_vector(s * finalknot)

    nearestpoint.data = signalpoint.data
    nearesttemp = nearestpoint.data

    diffvec = []

    for i in range(len(nearesttemp)):
        diffvec.append(nearesttemp[i] - signalparams[i])

    # The nearest point and its mismatch are calculated. An iterator is now constructed, the Finputs built and 2F is calculated for each template
    # derived from tbank

    # Create Iterator
    iterator = lp.CreateLatticeTilingIterator(tiling, s * finalknot)

    for i in range(s * finalknot):
        stats = lp.LatticeTilingStatistics(tiling, i)
        logging.info(f"Total   points in dim={i}: {stats.total_points}")
        logging.info(f"Minimum points in dim={i}: {stats.min_points}")
        logging.info(f"Maximum points in dim={i}: {stats.max_points}")

    fstatinputarray = []
    transformmatrices = []

    logging.info("Setting up Fstatistic input and transformation matrices")

    for seg in range(segs):
        logging.debug(f"Doing seg: {seg}")
        p0 = bf.knotslist[seg]
        p1 = bf.knotslist[seg + 1]
        tref = p0 + (p1 - p0) * trefsegfrac
        tstart = bf.knotslist[0]

        # Note that inject_params are only injected into data IF h0 != 0 otherwise an 'empty' signal is injected

        signalparamsseg = signalparams[seg * s : (seg + 2) * s]
        inject_params = gom.PWParamstoTrefParams(
            signalparamsseg, p0 - tstart, p1 - tstart, tref - tstart, s
        )

        transformmatrices.append(gom.ParamTransformationMatrix(p0, p1, tref, s))

        if SFTFiles or finputdata.inject_data:
            finput = fstatinputfromSFTFiles(
                SFTfmin,
                SFTfmax,
                finputdata,
                timeint=[p0, p1],
                sft_path=SFTdirectory,
                detectors=finputdata.detectors.data,
                inject_params=inject_params,
            )
        else:
            sft_catalog = SFTCatalog(p0, p1 - p0, finputdata)
            finput = fstatinput(
                SFTfmin,
                SFTfmax,
                sft_catalog,
                finputdata,
                inject_params,
                timeint=[p0, p1],
                trefsegfrac=trefsegfrac,
            )
        fstatinputarray.append(finput)

    logging.info("Finputs done")

    # 2F for the signal parameters and nearest template to them. They are then written to a file. The first two lines of data in these files are then
    # the signal parameters and nearest point and are therefore not a part of the search. The "NearestLatticeTilingPoint" method currently selects
    # a point well outside of the parameter space, so to prevent errors we have just set the nearestfstat value to -1.
    sigfstat = computefstat(
        signalparams,
        finputdata,
        fstatinputarray,
        transformmatrices,
        trefsegfrac=trefsegfrac,
        rtnsum=rtnsum,
    )

    all_sig_twoFs = [elem.twoF[0] for elem in sigfstat]
    sigtwoF = sum(all_sig_twoFs) / (len(all_sig_twoFs))

    sig_fstat_heap = [[sigtwoF, -1, -1, signalparams, sigfstat]]

    fstat_res_heap_to_file(
        out_directory + "/" + filename, sig_fstat_heap, fstats_first=True
    )

    mismatchfilename = "Mismatch_" + filename[:-4] + ".txt"
    # fstattemptofile(directory + "/" + mismatchfilename, signalparams, 0, sigfstat)
    # fstattemptofile(directory + "/" + mismatchfilename, nearestpoint.data, nearestmismatch, nearestfstat)
    fstat_res_heap_to_file(
        out_directory + "/" + mismatchfilename, sig_fstat_heap, fstats_first=False
    )

    # fin is iterated through by the NextLatticeTilingPoint, when fin == 0, the parameter space associated with tbank has been completed tiled.
    fin = -1

    # The temp heap, which will contain the highest 2F's and their associated templates. The elements of the tempheap have form (fstat, counter, template)
    template = lal.gsl_vector(s * finalknot)

    starttime = time.time()

    logging.info("Counting templates ...")
    tempcount = lp.TotalLatticeTilingPoints(iterator)
    logging.info(f"... finished: number of templates = {tempcount}")

    if mode == "template_count":
        return 0, 0, 0, tempcount

    lowest_mismatch_metric = 1000000000000000000000000
    lowest_mismatch_fstat = 1000000000000000000000000

    tqdm_file = open(filename[:-4] + "_tqdm_Output.txt", "w+")

    # If we are creating a mismatch histogram and are using the metric definition of the mismatch, we do not need to calculate any 2F's, so for efficiency
    # we do not calculate them and return the lowest found mismatch from the generated template bank
    if mode == "mis_hist" and not Fstat_mismatch:

        for counter in trange(tempcount, file=tqdm_file):
            fin = lp.NextLatticeTilingPoint(iterator, template)
            thistemp = copy.copy(template.data)

            tempdiffvec = [
                signalparams[i] - thistemp[i] for i in range(len(signalparams))
            ]

            tempmismatch = np.dot(tempdiffvec, np.dot(metric, tempdiffvec))

            if tempmismatch < lowest_mismatch_metric:
                lowest_mismatch_metric = tempmismatch

        return lowest_mismatch_metric, 0, 0, tempcount

    fstat_res_heap = []
    mismatch_fstat_res_heap = []

    fstatmax = 0
    lowest_mismatch = 100000000000

    fstat_hist_dic = {}

    # Iterate through the template bank and calculate 2F's for each as well as their mismatches. The templates with the greatest 2F's (and lowest mismatches)
    # are saved to a heap, which is then written to a file once the entire template bank has been iterated through
    for counter in trange(tempcount, file=tqdm_file):
        fin = lp.NextLatticeTilingPoint(iterator, template)
        thistemp = copy.copy(template.data)

        fstat = computefstat(
            thistemp,
            finputdata,
            fstatinputarray,
            transformmatrices,
            trefsegfrac=trefsegfrac,
            rtnsum=rtnsum,
        )

        all_twoFs = [res.twoF[0] for res in fstat]

        twoF = sum(all_twoFs) / len(all_twoFs)

        if twoF > fstatmax:
            fstatmax = twoF

        tempdiffvec = [signalparams[i] - thistemp[i] for i in range(len(signalparams))]

        temp_mismatch_metric = np.dot(tempdiffvec, np.dot(metric, tempdiffvec))
        temp_mismatch_fstat = (sigtwoF - twoF) / (sigtwoF - 4)

        if temp_mismatch_metric < lowest_mismatch_metric:
            lowest_mismatch_metric = temp_mismatch_metric

        if temp_mismatch_fstat < lowest_mismatch_fstat:
            lowest_mismatch_fstat = temp_mismatch_fstat

        if Fstat_mismatch:
            this_mismatch = temp_mismatch_fstat
        else:
            this_mismatch = temp_mismatch_metric

        if this_mismatch < lowest_mismatch:
            lowest_mismatch = this_mismatch

        # We include a counter when using heappush and heappushpop in the case that we comes acros the case where two f-statistics are equal
        # and hence will then try and find the next 'greatest' template, which will throw an error, given that the templates are arrays.
        if counter < tempsperfile:
            hq.heappush(fstat_res_heap, (twoF, this_mismatch, counter, thistemp, fstat))
            hq.heappush(
                mismatch_fstat_res_heap,
                (-this_mismatch, twoF, counter, thistemp, fstat),
            )
        else:
            hq.heappushpop(
                fstat_res_heap, (twoF, this_mismatch, counter, thistemp, fstat)
            )
            hq.heappushpop(
                mismatch_fstat_res_heap,
                (-this_mismatch, twoF, counter, thistemp, fstat),
            )

        if fstat_hist:
            add_fstat_to_dic(twoF, fstat_hist_dic)

    if fstat_hist:
        fstat_hist_file_name = out_directory + "/" + "Fstat_hist_" + filename
        dic_to_file(fstat_hist_file_name, fstat_hist_dic)

    tqdm_file.close()

    for i, elem in enumerate(mismatch_fstat_res_heap):
        new_mismatch = -1 * elem[0]
        mismatch_fstat_res_heap[i] = (new_mismatch, elem[1], elem[2], elem[3], elem[4])

    # The mismatch heap is sorted to have lowest mismatches first, while the fstat heap is sorted to have highest element first
    mismatch_fstat_res_heap.sort()
    fstat_res_heap.sort(reverse=True)

    # Check that end of template bank has been reached
    fin = lp.NextLatticeTilingPoint(iterator, template)
    if fin != 0:
        raise RuntimeError("end of template bank has not been reached")

    # General stats
    logging.info(f"Length of tempheap: {len(fstat_res_heap)}")
    logging.info(f"Templates per file is: {tempsperfile}")
    logging.info(f"Templates counted: {counter}")
    logging.info(f"Time Elapsed: {time.time() - starttime}")
    logging.info(f"Templates per second: {counter / (time.time() - starttime)}")

    # Minimum and maximum 2F's and mismatches
    logging.info(f"Sig f-stat: {sigtwoF}")
    logging.info("Running maximum f-stat: %s", str(fstatmax))
    logging.info("Heap    maximum f-stat: %s", str(fstat_res_heap[0][0]))
    logging.info(f"Running lowest mismatch: {lowest_mismatch}")
    logging.info(f"Heap lowest mismatch: {mismatch_fstat_res_heap[0][0]}")

    # Write both heaps to files
    fstat_res_heap_to_file(
        out_directory + "/" + filename, fstat_res_heap, fstats_first=True
    )
    fstat_res_heap_to_file(
        out_directory + "/" + mismatchfilename,
        mismatch_fstat_res_heap,
        fstats_first=False,
    )

    return lowest_mismatch_metric, lowest_mismatch_fstat, fstatmax, tempcount


# As the above, but this is what should be used for doing a search. It has all the superfluous code used for testing removed for fewer over heads and greater efficiency.
def pw_fstat_search_catalogue(
    SFTfmin,
    SFTfmax,
    tbank,
    finputdata,
    filename,
    tbankdirectory=".",
    SFTdirectory=".",
    trefsegfrac=0.0,
    rtnsum=False,
    tempsperfile=1000,
):

    segs = len(bf.knotslist) - 1
    logging.info("Knots: %s", str(bf.knotslist))

    finalknot = len(bf.knotslist)
    s = tbank.s

    # The below lines are required for finding the nearest lattice point to the signal parameters

    # Create LatticeTiling object
    tiling = lp.CreateLatticeTiling(s * finalknot)

    logging.debug("Doing metric")
    metric = scmm.PreCompMetric(s)

    logging.debug("Doing Bounds")
    tbe.setbounds(tiling, tbank)

    logging.debug("Setting tiling lattice and metric")

    # Set metric, mismatch and lattice type
    lp.SetTilingLatticeAndMetric(
        tiling, lp.TILING_LATTICE_ANSTAR, metric, tbank.maxmismatch
    )

    # Create Iterator
    iterator = lp.CreateLatticeTilingIterator(tiling, s * finalknot)

    for i in range(s * finalknot):
        stats = lp.LatticeTilingStatistics(tiling, i)
        logging.info(f"Total   points in dim={i}: {stats.total_points}")
        logging.info(f"Minimum points in dim={i}: {stats.min_points}")
        logging.info(f"Maximum points in dim={i}: {stats.max_points}")

    fstatinputarray = []
    transformmatrices = []

    logging.info("Setting up Fstatistic input and transformation matrices")

    for seg in range(segs):
        logging.debug(f"Doing seg: {seg}")
        p0 = bf.knotslist[seg]
        p1 = bf.knotslist[seg + 1]
        tref = p0 + (p1 - p0) * trefsegfrac
        tstart = bf.knotslist[0]

        transformmatrices.append(gom.ParamTransformationMatrix(p0, p1, tref, s))

        # We assume that SFTs are follow a directory structure "GW170817/Detector_X/SFT_X.sft". To extract all sfts for our FStatInput we need to have an sft_path with
        # the pattern Source/*/*.sft, to get all detectors and all of their SFTS. The SFTdirectory path is currently hardcoded in lalpulsar_PiecewiseSearch.py to get us
        # to the "GW170817" directory level. So when we come to remove this hardcoded path and allow it to be a user input, that user input also necessarily needs to
        # get to that level directory and should have the same structure. That is, the structure of "GW170817/Detector_X/SFT_X.sft".
        sft_path = SFTdirectory

        finput = fstatinputfromSFTFiles(
            SFTfmin,
            SFTfmax,
            finputdata,
            timeint=[p0, p1],
            sft_path=sft_path,
            detectors=tbank.detectors,
        )

        fstatinputarray.append(finput)

    logging.info("Finputs done")

    # fin is iterated through by the NextLatticeTilingPoint, when fin == 0, the parameter space associated with tbank has been completed tiled.
    fin = -1

    fstatmax = 0

    # The temp heap, which will contain the highest 2F's and their associated templates. The elements of the heap have form (fstat, mismatch, counter, template, FStatResults())
    template = lal.gsl_vector(s * finalknot)
    fstat_res_heap = []
    starttime = time.time()

    logging.info("Counting templates ...")
    tempcount = lp.TotalLatticeTilingPoints(iterator)
    logging.info(f"... finished: number of templates = {tempcount}")

    tqdm_file = open(filename[:-4] + "_tqdm_Output.txt", "w+")

    # Iterate through the template bank and calculate 2F's for each as well as their mismatches. The templates with the greatest 2F's (and lowest mismatches)
    # are saved to a heap, which is then written to a file once the entire template bank has been iterated through
    for counter in trange(tempcount, file=tqdm_file):
        fin = lp.NextLatticeTilingPoint(iterator, template)
        thistemp = copy.copy(template.data)

        fstat = computefstat(
            thistemp,
            finputdata,
            fstatinputarray,
            transformmatrices,
            trefsegfrac=trefsegfrac,
            rtnsum=rtnsum,
        )

        all_twoFs = [res.twoF[0] for res in fstat]

        twoF = sum(all_twoFs) / len(all_twoFs)

        if twoF > fstatmax:
            fstatmax = twoF

        # We include a counter when using heappush and heappushpop in the case that we comes acros the case where two f-statistics are equal
        # and hence will then try and find the next 'greatest' template, which will throw an error, given that the templates are arrays.
        if counter <= tempsperfile:
            hq.heappush(fstat_res_heap, (twoF, 0, counter, thistemp, fstat))
        else:
            hq.heappushpop(fstat_res_heap, (twoF, 0, counter, thistemp, fstat))

    tqdm_file.close()

    # Check that end of template bank has been reached
    fin = lp.NextLatticeTilingPoint(iterator, template)
    if fin != 0:
        raise RuntimeError("end of template bank has not been reached")

    # General stats
    logging.info(f"Length of tempheap: {len(fstat_res_heap)}")
    logging.info(f"Templates per file is: {tempsperfile}")
    logging.info(f"Templates counted: {counter}")
    logging.info(f"Time Elapsed: {time.time() - starttime}")
    logging.info(f"Templates per second: {counter / (time.time() - starttime)}")

    # Minimum and maximum 2F's and mismatches
    logging.info("Running maximum f-stat: %s", str(fstatmax))
    logging.info("Heap    maximum f-stat: %s", str(fstat_res_heap[-1][0]))

    fstat_res_heap.sort(reverse=True)

    # Write heap to file
    fstat_res_heap_to_file(tbankdirectory + "/" + filename, fstat_res_heap)

    return 0, 0, fstatmax, tempcount
