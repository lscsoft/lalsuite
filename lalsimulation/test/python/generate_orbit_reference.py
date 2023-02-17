#!/usr/bin/env python
# Copyright (C) 2021 Riccardo Sturani
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http: //www.gnu.org/licenses/>.

""" Code writing the values for regression tests for the SimInspiralSpinTaylorOrbitalDriver() function of
    lalsimulation/src/LALSimInspiralSpinTaylor.c.
"""

#import numpy as np

import lal
import lalsimulation as lalsim

DEFAULT_FILE = 'reviewed_evolveorbit.ini'
NEW_DATA_STR = '######### NEW DATASET #############\n'

lal_pars=lal.CreateDict()
phaseO = -1
spinO  = -1
tidalO = 0
phiref = 0.
deltaT = 0.0125
m1     = 22.
m2     = 22.
fstart = 30.
fref   = 30.
s1x    = -0.28360049
s1y    = 0.74726744
s1z    = 0.03481657
s2x    = -0.01902506
s2y    = 0.79400483
s2z    = 0.09233356
lnhatx = 0.
lnhaty = 0.
lnhatz = 1.
e1x    = 1.
e1y    = 0.
e1z    = 0.
approxs=["SpinTaylorT1","SpinTaylorT4","SpinTaylorT5"]

with open(DEFAULT_FILE, "w") as outfile:

    for approx in approxs:
        print(NEW_DATA_STR, file=outfile)
        print("[parameters]", file=outfile)
        print("approximant = %s" % approx, file=outfile)
        print("phaseO = %d" % phaseO, file=outfile)
        print("spinO = %d" % spinO, file=outfile)
        print("tidalO = %d" % tidalO, file=outfile)
        print("phiref = %f" % phiref, file=outfile)
        print("deltaT = %f" % deltaT, file=outfile)
        print("m1 = %f" % m1, file=outfile)
        print("m2 = %f" % m2, file=outfile)
        print("fstart = %f" % fstart, file=outfile)
        print("fref = %f" % fref, file=outfile)
        print("s1x = %10.8f" % s1x, file=outfile)
        print("s1y = %10.8f" % s1y, file=outfile)
        print("s1z = %10.8f" % s1z, file=outfile)
        print("s2x = %10.8f" % s2x, file=outfile)
        print("s2y = %10.8f" % s2y, file=outfile)
        print("s2z = %10.8f" % s2z, file=outfile)
        print("lnhatx = %f" % lnhatx, file=outfile)
        print("lnhaty = %f" % lnhaty, file=outfile)
        print("lnhatz = %f" % lnhatz, file=outfile)
        print("e1x = %f" % e1x, file=outfile)
        print("e1y = %f" % e1y, file=outfile)
        print("e1z = %f\n" % e1z, file=outfile)

        _, _, S1x, S1y, S1z, S2x, S2y, S2z, _, _, _, _, _, _ = lalsim.SimInspiralSpinTaylorOrbitalDriver(phiref, deltaT, m1*lal.MSUN_SI, m2*lal.MSUN_SI, fstart, fref, s1x, s1y, s1z, s2x, s2y, s2z, lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lal_pars, lalsim.GetApproximantFromString(approx))
        print("[spin-data]", file=outfile)
        print("spin1x =", end='', file=outfile)
        for x in S1x.data.data:
            print(" %.16e" % x, end='', file=outfile)
        print("\nspin1y =", end='', file=outfile)
        for x in S1y.data.data:
            print(" %.16e" % x, end='', file=outfile)
        print("\nspin1z =", end='', file=outfile)
        for x in S1z.data.data:
            print(" %.16e" % x, end='', file=outfile)
        print("\nspin2x =", end='', file=outfile)
        for x in S2x.data.data:
            print(" %.16e" % x, end='', file=outfile)
        print("\nspin2y =", end='', file=outfile)
        for x in S2y.data.data:
            print(" %.16e" % x, end='', file=outfile)
        print("\nspin2z =", end='', file=outfile)
        for x in S2z.data.data:
            print(" %.16e" % x, end='', file=outfile)
        print("\n", file=outfile)

outfile.close()
