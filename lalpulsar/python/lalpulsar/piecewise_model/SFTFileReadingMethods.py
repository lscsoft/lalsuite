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

from . import PWModelSimulations as pwsim
from . import TExpSimulations as tsim
from . import MOLSforGTE as mols
from . import BasisFunctions as bf
from . import ClassDefinitions as cd

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2
from math import factorial
from random import randrange

def SFTVector(fmin, fmax, buildsfts=None, filepattern="*.sft", constraints=None):
        if buildsfts:
                if buildsfts == "tsim":
                        tsim.runthisfile()
                elif buildsfts == "pwsim":
                        pwsim.runthisfile()

        catalogue = lp.SFTdataFind(filepattern, constraints)

        sftvector = lp.LoadSFTs(catalogue, fmin, fmax)

        return sftvector

def SFTVectorData(vector, realimag="m"):

        times = []
        mags  = []

        sft0 = vector.data[0]

        freqs = np.linspace(sft0.f0, sft0.f0 + sft0.deltaF * sft0.data.length, sft0.data.length)

        for sft in vector.data:
                GPStime = sft.epoch
                times.append(GPStime.gpsSeconds)

                sftmags = sft.data.data
                thismags = []

                for sftcoeffs in sftmags:
                        real = sftcoeffs.real
                        imag = sftcoeffs.imag

                        if realimag == "m":
                                thismags.append((real ** 2 + imag ** 2) ** 0.5)
                        elif realimag == "i":
                                thismags.append(imag)
                        elif realimag == "r":
                                thismags.append(real)

                mags.append(thismags)

        return [times, freqs, mags]

def DetectorPosVel(gpstime):
        posvec = lp.PosVel3D_t()
        velvec = lp.PosVel3D_t()

        laldetector = lal.CachedDetectors[5] # H1

        motiontype = lp.DETMOTION_ORBIT

        emphdata = lp.InitBarycenter("earth00-40-DE200.dat.gz", "sun00-40-DE200.dat.gz")

        lp.DetectorPosVel(posvec, velvec, gpstime, laldetector, emphdata, motiontype)

        pos = posvec.pos * lal.C_SI # I believe this is the appropriate conversion
        vel = velvec.vel * lal.C_SI

        return [pos, vel]

def plotsft(bigsftdata, plottemps=[], tbank=None, log=False, labels=[], tay=False, dopp=False):
        times = bigsftdata[0]
        freqs = bigsftdata[1]

        if log:
                maglist = np.matrix.transpose(np.log10(bigsftdata[2]))[::-1]
        else:
                maglist = np.matrix.transpose(np.array(bigsftdata[2]))[::-1]

        if tbank:
                bf.knotslist = tbank.knots

        if plottemps != []:
                for temp in plottemps:
                        res = 200
                        x = np.linspace(times[0], times[-1], res)

                        y = []

                        # This is a bit tricky. There is no way to know what the knots list is from the SFTs, so a tbank object with its knots should be supplied
                        bf.knotslist = tbank.knots

                        s = tbank.s
                        print(s)
                        print(temp)
                        if not tay:
                                coeffs = bf.allcoeffs(s)
                                params = mols.solsbyint(temp, s)

                        counter = 0
                        doppfactor = 1

                        for time in x:
                                modelvalue = 0

                                if tay:
                                        for i, param in enumerate(temp):
                                                modelvalue += (1 / factorial(i)) * param * (time - bf.knotslist[0]) ** i
                                else:
                                        modelvalue = mols.modelvalueatpoint(time, coeffs, params, ignoreintcheck=True)

                                # Doppler correction, slows down plotting significantly
                                if dopp:
                                        if (counter % (res // 10)) == 0:

                                                cosa = np.cos(pwsim.Alpha)
                                                sina = np.sin(pwsim.Alpha)
                                                cosd = np.cos(pwsim.Delta)
                                                sind = np.sin(pwsim.Delta);

                                                sourceskypos = [cosd * cosa, cosd * sina, sind]

                                                pos, vel = DetectorPosVel(time)

                                                doppfactor = (1 + np.dot(sourceskypos, vel) / lal.C_SI)
                                        counter += 1

                                y.append(doppfactor * modelvalue)

                        plt.plot(x, y)

        plt.imshow(maglist, cmap='viridis', extent=[times[0], times[-1], freqs[0], freqs[-1]], aspect='auto')
        plt.colorbar()
        plt.legend(labels)

        plt.show()

def plotsftwithtempsfromfile(path, knots, SFTfmin, SFTfmax, log=False, tay=False, dopp=False):
        signal = []
        nearest = []
        largest2F = []
        randomtemp = []
        lowestmismatchtemp = []

        randomtempnum = randrange(4, 10e5)

        if tay:
                paramtogoto = 2 + 6
        else:
                paramtogoto = 2 + 3 * len(knots)

        with open(path + "Temps_Sum.txt") as fstatfile:

                for linenum, line in enumerate(fstatfile):
                        if linenum == 0:
                                continue

                        if linenum == 1:
                                signalsplitstr = line.split()
                                signal = [float(elem) for elem in signalsplitstr[2:paramtogoto]]

                        if linenum == 2:
                                nearestsplitstr = line.split()
                                nearest = [float(elem) for elem in nearestsplitstr[2:paramtogoto]]

                        if linenum == 3:
                                largest2Fsplitstr = line.split()
                                largest2F = [float(elem) for elem in largest2Fsplitstr[2:paramtogoto]]

                        if linenum == randomtempnum:
                                randomtempsplitstr = line.split()
                                randomtemp = [float(elem) for elem in randomtempsplitstr[2:paramtogoto]]

        with open(path + "Temps_Sum_Mismatch.txt") as mismatchfile:
                for linenum, line in enumerate(mismatchfile):
                        if linenum == 0 or linenum == 1 or linenum == 2:
                                continue

                        elif linenum == 3:
                                lowestmismatchtempsplitstr = line.split()
                                lowestmismatchtemp = [float(elem) for elem in lowestmismatchtempsplitstr[2:paramtogoto]]

                        elif linenum > 3:
                                break

        tbank = cd.TBank()
        tbank.SetDefaultBNSR()
        tbank.knots = knots

        if tay:
                tbank.s = len(signal)
        else:
                tbank.s = len(signal) / len(knots)

        filepattern = path + "*.sft"

        vec = SFTVector(SFTfmin, SFTfmax, filepattern=filepattern)
        vecdata = SFTVectorData(vec)

        print("Signal: " + str(signal))
        print("Nearest: " + str(nearest))
        print("Largest 2F: " + str(largest2F))
        print("Random: " + str(randomtemp))
        print("Lowest Mismatch: " + str(lowestmismatchtemp))
        print()

        if len(randomtemp) != 0:
                plotsft(vecdata, plottemps=[signal, nearest, largest2F, lowestmismatchtemp, randomtemp], tbank=tbank, log=log, labels=["Signal", "Nearest Template", "Largest 2F", "Lowest Mismatch", "Random Temp"], tay=tay, dopp=dopp)
        else:
                plotsft(vecdata, plottemps=[signal, nearest, largest2F, lowestmismatchtemp], tbank=tbank, log=log, labels=["Signal", "Nearest Template", "Largest 2F", "Lowest Mismatch"], tay=tay, dopp=dopp)

def fstathistogram(path, bins, chi=False, log=False):
        fstats = []

        knots = 0

        with open(path) as fstatfile:
                for line in fstatfile:

                        if line.startswith("#"):
                                continue

                        linestr = line.split()

                        if float(linestr[0]) == 100000000000000 or float(linestr[0]) == 10000.0000000000:
                                knots = (len(linestr) - 2) / 3
                                continue

                        fstats.append(float(linestr[0]))

        if chi:
                x = np.arange(0, np.max(fstats), 0.1)
                plt.plot(x, chi2.pdf(x, df=4 * (knots - 1)))
        if log:
                plt.yscale('log')

        plt.hist(fstats, bins=bins, density=True)
        plt.xlabel("2F")
        plt.ylabel("Normalised Frequency")
        plt.show()

def paramhistogram(path, dim, bins, log=False):
        params = []

        with open(path) as fstatfile:
                for line in fstatfile:

                        if line.startswith("#"):
                                continue

                        linestr = line.split()

                        params.append([float(linestr[2 + dim]), float(linestr[2 + dim + 1])])

        if log:
                plt.yscale('log')

        # Will plot parameter values as dot points

        #plt.plot([0, len(params)], [params[0], params[0]], color='red', label="Signal")
        #plt.plot([0, len(params)], [params[1], params[1]], color='green', label="Nearest Point")
        #plt.plot([0, len(params)], [params[2], params[2]], color='orange', label="Highest 2F Template")
        #params.sort()
        #plt.plot(params)

        # Plots a line between the value of one parameter and the next parameter (scaled to the same order of magnitude) for each set of parameters in
        # params (if the params.append line in the for loop above appends a list of two adjacent parameters)

        # Scales everything to approximately the same order of magnitude
        firstparamlog = int(np.log10(np.abs(params[0][0])))
        secondparamlog = int(np.log10(np.abs(params[0][1])))

        for twoparams in params:
                plt.plot([np.abs(twoparams[0]) / (10 ** firstparamlog), np.abs(twoparams[1]) / (10 ** secondparamlog)], color='blue')

        plt.plot([np.abs(params[0][0]) / (10 ** firstparamlog), np.abs(params[0][1]) / (10 ** secondparamlog)], color='red', label="Signal")
        plt.plot([np.abs(params[1][0]) / (10 ** firstparamlog), np.abs(params[1][1]) / (10 ** secondparamlog)], color='green', label="Nearest Point")
        plt.plot([np.abs(params[2][0]) / (10 ** firstparamlog), np.abs(params[2][1]) / (10 ** secondparamlog)], color='orange', label="Highest 2F Template")

        # Will plot parameter values in a histogram
        """
        plt.axvline(params[0], color='red', label="Signal")
        plt.axvline(params[1], color='green', label="Nearest Point")
        plt.axvline(params[2], color='orange', label="Highest 2F Template")

        plt.hist(params, bins=bins, density=True)
        plt.xlabel("2F")
        plt.ylabel("Normalised Frequency")
        """

        plt.legend()
        plt.show()
