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

from . import BasisFunctions as bf
from . import GTEandOtherMethods as gom
from . import EstimatingKnots as ek

import numpy as np

# Converts a template given as a list of parameters into a template made up of a list of knot templates
def templatetoknottemplates(temp, s):
        template = []
        knottemplate = []

        for i, param in enumerate(temp):
                if len(knottemplate) == s - 1:
                        knottemplate.append(param)
                        template.append(knottemplate)
                        knottemplate = []
                else:
                        knottemplate.append(param)

        return template

# Checks to see that every value in values lies within the defined global minimum and maximum range and are within the appropriate tolerance of one another
def checkrangesandtolerances(values, globalmin, globalmax, tolerance):

        outofrange = 0
        outoftolerance = 0

        for i, val in enumerate(values):
                if not 0.999 * globalmin <= val <= globalmax * 1.001:
                        print("Failed global ranges: " + str(val))
                        outofrange += 1

                if i == 0:
                        continue

                Dp = bf.knotslist[i] - bf.knotslist[i - 1]

                lower_bound = values[i - 1] * (1 - Dp * tolerance)
                upper_bound = values[i - 1]

                if not 0.999 * lower_bound <= val <= upper_bound * 1.001:
                        #print(values)
                        print("Failed tolerance: " + str([lower_bound, val, upper_bound]))
                        outoftolerance += 1

        return outofrange + outoftolerance

# Takes a template (list of parameters) and tbank object and returns a list [nfails, kfails, gtefails] indicating which criteria and how many times this template fails those criteria
def is_valid_temp(temp, tbank):

        template = templatetoknottemplates(temp, tbank.s)

        fmin     = tbank.fmin
        fmax     = tbank.fmax
        fmaxtrue = tbank.fmaxtrue
        nmin     = tbank.nmin
        nmax     = tbank.nmax
        nmin0    = tbank.nmin0
        nmax0    = tbank.nmax0
        ntol     = tbank.ntol
        taumin   = tbank.taumin
        taumax   = tbank.taumax
        ktol     = tbank.ktol

        kmin = gom.kwhichresultsingivenhalflife(taumax, fmaxtrue, nmax)
        kmax = gom.kwhichresultsingivenhalflife(taumin, fmaxtrue, nmax)

        ns = []
        ks = []

        for i, knottemp in enumerate(template):

                n = -1

                if tbank.s == 3:
                        n = knottemp[2] * knottemp[0] / (knottemp[1] ** 2)

                elif tbank.s == 2:

                        if i == 0:
                                continue

                        f0prev = template[i - 1][0]
                        f1prev = template[i - 1][1]

                        f0 = knottemp[0]
                        f1 = knottemp[1]

                        segment_length = bf.knotslist[i] - bf.knotslist[i - 1]

                        n = 1 + (f0prev * f1 - f0 * f1prev) / (f1prev * f1 * segment_length);

                        if i == 1:
                                ns.append(n)

                try:
                        k = - knottemp[1] / (knottemp[0] ** n)
                except ZeroDivisionError:
                        k = 10 ** 6
                except OverflowError:
                        k = 0

                ns.append(n)
                ks.append(k)
        """
        print([nmin, nmax])
        print([kmin, kmax])
        print(ns)
        print(ks)
        print()
        """

        nfails = 0
        kfails = 0
        f0fails = 0
        f1fails = 0
        f2fails = 0

        if not 0.9999999 * nmin0 <= ns[0] <= nmax0 * 1.0000001:
                print("First knot braking index fail: " + str(ns[0]))
                nfails += 1

        nfails += checkrangesandtolerances(ns, nmin, nmax, ntol)
        kfails += checkrangesandtolerances(ks, kmin, kmax, ktol)

        for i, knottemp in enumerate(template):
                if i == 0:
                        if not fmin <= knottemp[0] <= fmax:
                                #print("First knot frequency fail: " + str([fmin, knottemp[0], fmax]))
                                f0fails += 1
                        continue

                f0n1 = template[i - 1][0]
                Dp = bf.knotslist[i] - bf.knotslist[i - 1]

                ntolmin = ns[i - 1] * (1 - Dp * ntol)
                ntolmax = ns[i - 1]
                ktolmin = ks[i - 1] * (1 - Dp * ktol)
                ktolmax = ks[i - 1]

                try:
                        if not 0.999 * gom.gtederivs(Dp, f0n1, ntolmax, ktolmax, 0) <= knottemp[0] <= gom.gtederivs(Dp, f0n1, ntolmin, ktolmin, 0) * 1.001:
                                #print("f0 fail")
                                print("F0 failed by: " + str([knottemp[0], gom.gtederivs(Dp, f0n1, ntolmax, ktolmax, 0) - knottemp[0], knottemp[0] - gom.gtederivs(Dp, f0n1, ntolmin, ktolmin, 0)]))
                                f0fails += 1
                except OverflowError:
                        f0fails += 1

                try:
                        if not 1.001 * gom.gtederivs(Dp, f0n1, ntolmax, ktolmax, 1) <= knottemp[1] <= gom.gtederivs(Dp, f0n1, ntolmin, ktolmin, 1) * 0.999:
                                #print("f1 fail")
                                print("F1 failed by: " + str([knottemp[1], gom.gtederivs(Dp, f0n1, ntolmax, ktolmax, 1) - knottemp[1], knottemp[1] - gom.gtederivs(Dp, f0n1, ntolmin, ktolmin, 1)]))
                                f1fails += 1
                except OverflowError:
                        f1fails += 1

                if tbank.s == 3:
                        try:
                                if not 0.999 * gom.gtederivs(Dp, f0n1, ntolmin, ktolmin, 2) <= knottemp[2] <= gom.gtederivs(Dp, f0n1, ntolmax, ktolmax, 2) * 1.001:
                                        #print("f2 fail")
                                        print("F2 failed by: " + str([knottemp[2], gom.gtederivs(Dp, f0n1, ntolmax, ktolmax, 2) - knottemp[2], knottemp[2] - gom.gtederivs(Dp, f0n1, ntolmin, ktolmin, 2)]))
                                        f2fails += 1
                        except OverflowError:
                                f2fails += 1

        verdict = True

        if nfails + kfails + f0fails + f1fails + f2fails != 0:
                verdict = False

        return [verdict, nfails, kfails, f0fails, f1fails, f2fails]

def checkfstatfile(path, tbank, maxtemps = 10):

        templates = []

        tempsfailed = []

        with open(path) as fstatfile:
                for line in fstatfile:

                        if line.startswith("#"):
                                continue

                        linestr = line.split()[2:]

                        temp = []

                        for elem in linestr:
                                temp.append(float(elem))

                        templates.append(temp)

        for i, template in enumerate(templates):
                if i > maxtemps:
                        break
                print()
                print(i)
                fails = isvalidtemp(template, tbank)

                if fails != [0, 0, 0]:
                        tempsfailed.append([i, fails])

        return tempsfailed

# Builds a valid template using the lp.PiecewiseParameterBounds method. Builds a template by sequentially calculating the bounds on a dimension using the values
# of the parameters prior. The value of the parameter for the current dimension is then chosen to be halfway between these bounds. This is repeated until an entire
# template has been consturcted.
def auto_valid_template(tbank, build_knots=[], reset=1):

        template = []

        kmin = gom.kwhichresultsingivenhalflife(tbank.taumax, tbank.fmax, tbank.nmax)
        kmax = gom.kwhichresultsingivenhalflife(tbank.taumin, tbank.fmax, tbank.nmax)

        if build_knots == []:
                ek.allidealisedknots(tbank.s, tbank.dur, 40, tbank.fmaxtrue, tbank.nmax, tbank.taumin, tbank.mismatch)
        else:
                bf.knotslist = build_knots

        for dim in range(tbank.s * len(bf.knotslist)):

                print("At dim: " + str(dim))
                print(template)

                if dim == 0:
                        template.append(1/2 * (tbank.fmax + tbank.fmin))
                        continue

                segment = int(np.floor(dim / tbank.s))

                seg_length = bf.knotslist[segment] - bf.knotslist[segment - 1]
                print(segment)
                print(seg_length)

                lower = lp.PiecewiseParameterBounds(dim, template, -1, tbank.fmin, tbank.fmax, tbank.nmin, tbank.nmax, tbank.ntol, kmin, kmax, tbank.ktol, seg_length, reset)
                upper = lp.PiecewiseParameterBounds(dim, template,  1, tbank.fmin, tbank.fmax, tbank.nmin, tbank.nmax, tbank.ntol, kmin, kmax, tbank.ktol, seg_length, reset)

                print([lower, upper])

                template.append(1/2 * (upper + lower))

        return template
