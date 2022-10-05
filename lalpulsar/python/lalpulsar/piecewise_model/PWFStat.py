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

from . import SemicoherentMetricMethods as scmm
from . import GTEandOtherMethods as gom
from . import BasisFunctions as bf
from . import TBankEstimates as tbe
from . import TemplatePlottingMethods as tpm
from . import TemplateChecks as tc

import numpy as np
import heapq as hq
import time
import logging
import copy

import matplotlib.pyplot as plt

log = logging.getLogger()
log.setLevel(logging.DEBUG)

# Returns our FstatInput object, currently uses SFTs built in pwsim
# fmax (and fmin) should be chosen with the same reasoning outlined in the pwsim notbook. fmax = f0 + 10^-4 * f0 + 58/TSFT,
# fmin is the same except fmin = f0 - .... AND fmin should be lower than the minimum value any given template we might use
# over the time period of the relevant SFTs!

def fstatinputfromSFTFiles(SFTfmin, SFTfmax, finputdata, timeint=[-1, -2], directory=""):

        Fstat_opt_args = lp.FstatOptionalArgs(lp.FstatOptionalArgsDefaults)
        Fstat_opt_args.FstatMethod = lp.FMETHOD_DEMOD_BEST
        Fstat_opt_args.sourceDeltaT = finputdata.sourceDeltaT

        dfreq = finputdata.dfreq
        emphdata = lp.InitBarycenter("earth00-40-DE405.dat.gz", "sun00-40-DE405.dat.gz")

        constraints = lp.SFTConstraints()

        constraints.minStartTime = lal.LIGOTimeGPS(timeint[0])
        constraints.maxStartTime = lal.LIGOTimeGPS(timeint[1])
        catalogue = lp.SFTdataFind(directory + "/*.sft", constraints)
        fstat = lp.CreateFstatInput(catalogue, SFTfmin, SFTfmax, dfreq, emphdata, Fstat_opt_args)
        return fstat

# Create fake sft catalog
def SFTCatalog(tstart, Tdata, Tsft, detector='H1'):

        sft_ts = lp.MakeTimestamps(tstart, Tdata, Tsft, 0)
        sft_catalog = lp.AddToFakeSFTCatalog(None, detector, sft_ts)

        return sft_catalog

# Builds SFTs in memory rather than files
def fstatinput(SFTfmin, SFTfmax, sft_catalog, finputdata, dopplerparams, timeint=[-1, -2], trefsegfrac=0.):

        Tdata = timeint[1] - timeint[0]
        dfreq = finputdata.dfreq

        # create default F-statistic optional arguments
        Fstat_opt_args = lp.FstatOptionalArgs(lp.FstatOptionalArgsDefaults)
        Fstat_opt_args.FstatMethod = lp.FMETHOD_DEMOD_BEST
        Fstat_opt_args.sourceDeltaT = finputdata.sourceDeltaT

        h0            = finputdata.h0
        cosi          = finputdata.cosi
        psi           = finputdata.psi
        phi0          = finputdata.phi0
        tref          = timeint[0] + Tdata * trefsegfrac
        Alpha         = finputdata.Alpha
        Delta         = finputdata.Delta
        assume_sqrtSh = finputdata.h0

        # create injection parameters
        Fstat_signal = lp.CreatePulsarParamsVector(1);
        Fstat_signal.data[0].Amp.aPlus       = 0.5 * h0 * (1.0 + cosi * cosi)
        Fstat_signal.data[0].Amp.aCross      = h0 * cosi
        Fstat_signal.data[0].Amp.psi         = psi
        Fstat_signal.data[0].Amp.phi0        = phi0
        Fstat_signal.data[0].Doppler.refTime = tref
        Fstat_signal.data[0].Doppler.Alpha   = Alpha
        Fstat_signal.data[0].Doppler.Delta   = Delta
        print(dopplerparams)
        # Set doppler params. Done as a for loop in case dopplerparams changes length (it shouldn't, but just in case)
        for i in range(len(dopplerparams)):
                Fstat_signal.data[0].Doppler.fkdot[i] = dopplerparams[i]

        Fstat_opt_args.injectSources = Fstat_signal
        Fstat_assume_noise           = lp.MultiNoiseFloor()
        Fstat_assume_noise.length    = 1
        Fstat_assume_noise.sqrtSn[0] = assume_sqrtSh
        Fstat_opt_args.assumeSqrtSX  = Fstat_assume_noise

        ephemerides = lp.InitBarycenter('earth00-19-DE405.dat.gz', 'sun00-19-DE405.dat.gz')

        # create F-statistic input data
        Fstat_input = lp.CreateFstatInput(sft_catalog, SFTfmin, SFTfmax, dfreq, ephemerides, Fstat_opt_args)

        return Fstat_input

# Returns our fstat results for a given template (using sft's from fstatin).
def computefstat(template, finputdata, fstatinputarray, transformmatrices, trefsegfrac=0., rtnsum=False):

        s = int(len(transformmatrices[0]) / 2)

        segs = int(len(template)/s - 1)

        dopparams       = lp.PulsarDopplerParams()
        dopparams.Alpha = finputdata.Alpha
        dopparams.Delta = finputdata.Delta

        fs = []
        counter = 0

        for seg in range(segs):
                counter += 1
                segtemp = template[seg * s:(seg + 2) * s]

                tref = bf.knotslist[seg] + (bf.knotslist[seg + 1] - bf.knotslist[seg]) * trefsegfrac
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

                dopparams.fkdot[0:len(dopplerparams)] = dopplerparams

                numFreqBins   = 1
                whatToCompute = lp.FSTATQ_2F_PER_DET + lp.FSTATQ_2F
                fstatresults  = lp.FstatResults()
                fstatresults  = lp.ComputeFstat(fstatresults, fstatin, dopparams, numFreqBins, whatToCompute)

                fs.append(fstatresults.twoF[0])

        fstat = sum(fs)/len(fs)

        if rtnsum:
                return sum(fs)
        else:
                return fstat

# Will append a single fstatistic and template to a file
def fstattemptofile(filename, template, vara, varb):
        with open(filename, "a+") as reader:

                # Note, the addition of each templateline string starting with a space makes that string have a length of 16, not 15
                templateline = " ".join("{:20.12E}".format(x) for x in template)
                line = "{:20.10f}".format(vara) + "     " + "{:20.4f}".format(varb) + "     " + templateline + "\n"

                reader.write(line)

# Will write a python heap to a file name. The heap should have form (Vara, Varb, template). Ascending refers to whether we write the heap elements from the lowest value of Vara (ascending = 1),
# or the highest value of Vara (ascending = -1). This method assumes that the heap has elements ordered from lowest Vara to highest Vara.
def writetempheaptofile(tempheap, filename, ascending=-1):

        mismatchheap = []
        mismatchfilename = filename[:-4] + "_Mismatch.txt"
        counter = 0

        for ft in tempheap[::ascending]:
                vara = ft[0]
                varb = ft[1]
                temp = ft[3]

                fstattemptofile(filename, temp, vara, varb)

# For a template bank associated with the parameter space associated with tbank, finds the 'tbank.maxtemps' number of templates with the greatest FStat.
# SFTfmin and SFTfmax are the frequency range of SFTs we want to load in, tbank is the parameter space parameters we use. filename is the name of the file
# where the F-statistics and their corresponding templates will be saved. tbank.maxtemps is the top 2F's and templates to store and saved. directory is the
# path to where filename will be saved. trefsegfrac is how far along each segment the reference time will be defined. E.g. If set to 0, tref will be at the
# start of each segment, if trefsegfrac = 0.5, then tref will occur half way through each segment and so on. signalparams are the signal parameters of the signal
# injected into the SFTs and rtnsum tells us whether to save the summed 2F across all segments or the average. rtnsum = False means the average is saved,
# True and the sum will be saved
def semifstatcatalogue(SFTfmin, SFTfmax, tbank, finputdata, signalparams, filename, reset, padding_flags_bbox=[], padding_flags_int=[], directory="", trefsegfrac=0., rtnsum=False, SFTFiles=False, tempsperfile=1000):
        print("Beginning semifstatcatalogue")
        if padding_flags_bbox == []:
                padding_flags_bbox = [0] * tbank.s * len(bf.knotslist)
        if padding_flags_int == []:
                padding_flags_int  = [0] * tbank.s * len(bf.knotslist)

        max_mismatch = tbank.mismatch

        #gom.PlotPWModel(signalparams, show=True, label="", linewidth=2)

        segs = len(bf.knotslist) - 1
        logging.info("Knots: %s", str(bf.knotslist))

        finalknot = len(bf.knotslist)
        s = tbank.s
        print(bf.knotslist)

        # The below lines are required for finding the nearest lattice point to the signal parameters

        # Create LatticeTiling object
        tiling = lp.CreateLatticeTiling(s * finalknot)

        print("Doing metric")
        metric = scmm.PreCompMetric(s)

        print("Doing Bounds")
        # Set Bounds and reset knots (as knots are reset when setting bounds)
        tbe.setbounds(tiling, padding_flags_bbox, padding_flags_int, tbank, reset)

        print("Setting tiling lattice and metric")
        # Set metric, mismatch and lattice type
        lp.SetTilingLatticeAndMetric(tiling, lp.TILING_LATTICE_ANSTAR, metric, max_mismatch)

        signalpoint      = lal.gsl_vector(s * finalknot)
        signalpoint.data = signalparams

        nearestpoint = lal.gsl_vector(s * finalknot)
        UINT8vec     = lal.CreateUINT8Vector(s * finalknot)

        # The below was used to see which point is the closest in the template bank, however, it is very slow. It has been commented out
        # for now, and nearesttemp and nearestpoint.data have just been set to be the signal parameters
        """
        locator      = lp.CreateLatticeTilingLocator(tiling)
        lp.NearestLatticeTilingPoint(locator, signalpoint, nearestpoint, UINT8vec)
        nearesttemp = nearestpoint.data
        """

        nearestpoint.data = signalpoint.data
        nearesttemp       = nearestpoint.data

        diffvec = []

        for i in range(len(nearesttemp)):
                                diffvec.append(nearesttemp[i] - signalparams[i])

        nearestmismatch = np.dot(diffvec, np.dot(metric, diffvec))

        # The nearest point and its mismatch are calculated. An iterator is now constructed, the Finputs built and 2F is calculated for each template
        # derived from tbank

        # Create Iterator
        iterator = lp.CreateLatticeTilingIterator(tiling, s * finalknot)

        fstatinputarray = []
        transformmatrices = []
        condmatrices = []
        print("Doing Finputs matrices")
        for seg in range(segs):
                print("Doing seg: " + str(seg))
                p0     = bf.knotslist[seg]
                p1     = bf.knotslist[seg + 1]
                tref   = p0 + (p1 - p0) * trefsegfrac
                tstart = bf.knotslist[0]

                signalparamsseg = signalparams[seg * s: (seg + 2) * s]
                dopplerparams   = gom.PWParamstoTrefParams(signalparamsseg, p0 - tstart, p1 - tstart, tref - tstart, s)

                transformmatrices.append(gom.ParamTransformationMatrix(p0, p1, tref, s))

                if SFTFiles:
                        finput = fstatinputfromSFTFiles(SFTfmin, SFTfmax, finputdata, timeint=[p0, p1], directory="")
                else:
                        sft_catalog = SFTCatalog(p0, p1 - p0, finputdata.Tsft, detector=finputdata.detector)
                        finput = fstatinput(SFTfmin, SFTfmax, sft_catalog, finputdata, dopplerparams, timeint=[p0, p1], trefsegfrac=trefsegfrac)

                fstatinputarray.append(finput)

        logging.info("Finputs done")

        # 2F for the signal parameters and nearest template to them. They are then written to a file. The first two lines of data in these files are then
        # the signal parameters and nearest point and are therefore not a part of the search. The "NearestLatticeTilingPoint" method currently selects
        # a point well outside of the parameter space, so to prevent errors we have just set the nearestfstat value to -1.
        sigfstat     = computefstat(signalparams, finputdata, fstatinputarray, transformmatrices, trefsegfrac=trefsegfrac, rtnsum=rtnsum)
        nearestfstat = -1 #computefstat(nearestpoint.data, finputdata, fstatinputarray, transformmatrices, trefsegfrac=trefsegfrac, rtnsum=rtnsum)

        fstattemptofile(directory + "/" + filename, signalparams, sigfstat, -1)
        fstattemptofile(directory + "/" + filename, nearestpoint.data, nearestfstat, nearestmismatch)

        mismatchfilename = filename[:-4] + "_Mismatch.txt"
        fstattemptofile(directory + "/" + mismatchfilename, signalparams, -1, sigfstat)
        fstattemptofile(directory + "/" + mismatchfilename, nearestpoint.data, nearestmismatch, nearestfstat)

        # fin is iterated through by the NextLatticeTilingPoint, when fin == 0, the parameter space associated with tbank has been completed tiled.
        # The parameter maxtempsreached indicated whether the tempheap has exceed the specified maximum number of templates to store (specified in
        # the tbank value tbank.maxtemps. The counter is used to print out how many templates have been counted, a print out message is made using
        # counter every thousand templates
        fin = -1

        fstatmax = 0
        maxtempsreached = False
        counter = 0

        # The temp heap, which will contain the highest 2F's and their associated templates. The elements of the tempheap have form (fstat, counter, template)
        template = lal.gsl_vector(s * finalknot)
        tempheap = []
        mismatchheap = []
        starttime = time.time()

        templates = []

        valid_temps   = 0
        invalid_temps = 0
        valid_stats = np.array([0, 0, 0, 0, 0])
        #freqp0 = []
        #f1dotp0 = []
        print("max temps is: " + str(tbank.maxtemps))
        while counter <= tbank.maxtemps:
                fin = lp.NextLatticeTilingPoint(iterator, template)
                thistemp = copy.copy(template.data)

                """
                this_valid_stats = tc.is_valid_temp(thistemp, tbank)
                valid_stats += this_valid_stats[1:]

                if this_valid_stats[0]:
                        valid_temps += 1
                else:
                        invalid_temps += 1
                """
                #freqp0.append(thistemp[0])
                #f1dotp0.append(thistemp[1])

                templates.append(thistemp)

                fstat = computefstat(thistemp, finputdata, fstatinputarray, transformmatrices, trefsegfrac=trefsegfrac, rtnsum=rtnsum)

                if fstat > fstatmax:
                        fstatmax = fstat

                tempdiffvec  = [signalparams[i] - thistemp[i] for i in range(len(signalparams))]
                tempmismatch = np.dot(tempdiffvec, np.dot(metric, tempdiffvec))
                #tempmismatch = (sigfstat - fstat) / (sigfstat - 4)

                # We include a counter when using heappush and heappushpop in the case that we comes acros the case where two f-statistics are equal
                # and hence will then try and find the next 'greatest' template, which will throw an error, given that the templates are arrays.
                if counter < tempsperfile:
                        hq.heappush(tempheap,     (fstat, tempmismatch, counter, thistemp))
                        hq.heappush(mismatchheap, (-tempmismatch, fstat, counter, thistemp))
                else:
                        hq.heappushpop(tempheap,     (fstat, tempmismatch, counter, thistemp))
                        hq.heappushpop(mismatchheap, (-tempmismatch, fstat, counter, thistemp))

                counter += 1

                if counter % 10000 == 0:
                        print("Counter is at: " + str(counter))
                        logging.info("Counter is at: %s", str(counter))

                if fin == 0:
                        break

        signal_valid_stats = tc.is_valid_temp(signalparams, tbank)
        print("Signal valid stats: " + str(signal_valid_stats))

        print("Valid and Invalid temps are: " + str([valid_temps, invalid_temps]))
        print("Valid stats are: " + str(valid_stats))

        print("Lengt of tempheap: " + str(len(tempheap)))
        print("Temps per file is: " + str(tbank.maxtemps))
        print("Templates counted: " + str(counter))
        print("Time Elapsed: " + str(time.time() - starttime))
        print("Templates per second: " + str(counter / (time.time() - starttime)))

        logging.info("Last template: %s", str(template.data))

        logging.info("Maximum f-stat: %s", str(fstatmax))

        tempheap.sort()
        mismatchheap.sort()

        print("Sig f-stat: " + str(sigfstat))
        print("Max f-stat: " + str(tempheap[-1][0]))
        print("Lowest mismatch: " + str(-mismatchheap[-1][0]))
        print("Mismatch of loudest: " + str(tempheap[-1][1]))

        logging.info("Heaps sorted")

        # Write both heaps to files
        writetempheaptofile(tempheap,     directory + "/" + filename,         ascending=-1)
        writetempheaptofile(mismatchheap, directory + "/" + mismatchfilename, ascending=-1)

        # Some model plotting

        # Firstly, we plot the loudest, lowest mismatch and signal templates together
        highest2F = tempheap[-1][-1]
        lowestmu = mismatchheap[0][-1]

        gom.PlotPWModel(highest2F,    show=False, label="Loudest",         linewidth=4)
        gom.PlotPWModel(lowestmu,     show=False, label="Lowest mismatch", linewidth=3)
        gom.PlotPWModel(signalparams, show=False, label="Signal",          linewidth=2)
        plt.legend()
        plt.show()

        # Secondly, we plot the template bank for two parameters relative to given knots (the dependency on knot
        # is contained within the freqindex and f1dotindex, although these parameters may not actually reference
        # a frequency and first spindown parameter).

        kmin = gom.kwhichresultsingivenhalflife(tbank.taumax, tbank.fmax, tbank.nmax)
        kmax = gom.kwhichresultsingivenhalflife(tbank.taumin, tbank.fmax, tbank.nmax)

        f0s = np.linspace(tbank.fmin, tbank.fmax, 100)

        f1upper = [-kmin * f0 ** tbank.nmin for f0 in f0s]
        f1lower = [-kmax * f0 ** tbank.nmax for f0 in f0s]

        for seg in range(1):

                freqp0  = []
                f1dotp0 = []

                freqindex = 0
                f1dotindex = 1

                for template in templates:
                        freqp0.append(template[freqindex])
                        f1dotp0.append(template[f1dotindex])

                templates.append(signalparams)

                crosssectiontemp = highest2F #templates[int(np.floor(len(templates) / 2))]

                tpm.plottemplates(templates, crosssectiontemp, freqindex, f1dotindex, metric, tbank.mismatch, show=False)

                signalf0f1   = [signalparams[freqindex], signalparams[f1dotindex]]
                highestf0f1  = [highest2F[freqindex], highest2F[f1dotindex]]
                lowestmuf0f1 = [lowestmu[freqindex], lowestmu[f1dotindex]]

                #plt.plot(f0s, f1upper, label="F1 Upper Bound", color="blue")
                #plt.plot(f0s, f1lower, label="F1 Lower Bound", color="blue")
                plt.scatter(freqp0, f1dotp0, label="Templates", s=2.5, color="blue")
                plt.scatter([signalf0f1[0]], [signalf0f1[1]], label="Signal", color="green")
                plt.scatter([highestf0f1[0]], [highestf0f1[1]], label="Loudest", s=40, color="red")
                plt.scatter([lowestmuf0f1[0]], [lowestmuf0f1[1]], label="Lowest mismatch", s=30, color="orange")
                plt.legend()
                plt.xlabel("Frequency")
                plt.ylabel("F1dot")
                plt.title("Knot: p_" + str(seg))
                plt.show()
