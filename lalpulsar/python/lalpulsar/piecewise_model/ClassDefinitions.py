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

from . import BasisFunctions as bf
from . import EstimatingKnots as ek
import numpy as np

class TBank:
        s          = None
        fmin       = None
        fmax       = None
        fmaxtrue   = None
        nmin       = None
        nmax       = None
        nmin0      = None
        nmax0      = None
        ntol       = None
        taumin     = None
        taumax     = None
        ktol       = None
        knots      = []
        mismatch   = None
        dur        = None
        maxtemps   = None
        flagoption = 0

        temps = 0
        currenttemplate = 0

        def SetTBankParams(self, s, fmin, fmax, fmaxtrue, nmin, nmax, nmin0, nmax0, ntol, taumin, taumax, ktol, knots, mismatch, dur, maxtemps=1000000):
                self.s        = s
                self.fmin     = fmin
                self.fmax     = fmax
                self.fmaxtrue = fmaxtrue
                self.nmin     = nmin
                self.nmax     = nmax
                self.nmin0    = nmin0
                self.nmax0    = nmax0
                self.ntol     = ntol
                self.taumin   = taumin
                self.taumax   = taumax
                self.ktol     = ktol
                self.knots    = knots
                self.mismatch = mismatch
                self.dur      = dur
                self.maxtemps = maxtemps

                bf.knotslist = knots

        def SetDefaultBNSR(self, maxtemps=1000000):
                self.s        = 3
                self.fmin     = 999   #Abitrary minimum
                self.fmax     = 1000
                self.fmaxtrue = self.fmax
                self.nmin     = 2
                self.nmax     = 5
                self.nmin0    = self.nmin
                self.nmax0    = self.nmax
                self.ntol     = 0.01 / (1/3 * 24 * 3600 * 365) # Four months, arbitrary
                self.taumin   = 3 * 3600
                self.taumax   = 6 * 3600
                self.ktol     = 0.01 / (1/3 * 24 * 3600 * 365) # Four months, arbitrary
                self.mismatch = 0.2
                self.dur      = 6 * 3600
                self.maxtemps = maxtemps

                ek.allidealisedknots(self.s, self.dur, 40, self.fmaxtrue, self.nmax, self.taumin, mu=self.mismatch)
                knotnum = ek.getknotnum(self.s,self.dur, self.fmaxtrue, self.nmax, self.taumin, mu=self.mismatch)
                knots = bf.knotslist[:knotnum + 1]

                # Knots for this object are:
                # [0, 132.275, 367.4173537605237, 733.6451197773617, 1416.1006518110376, 2551.2430055715613, 4091.8124206136563, 6172.951250697846, 9835.228910866226, 13497.506571034606,
                # 18325.35364624508, 28367.585631795242]

                self.knots = knots
                bf.knotslist = knots

                self.maxtemps = maxtemps

        def SetDefault1987A(self, maxtemps=1000000):
                self.s        = 3
                self.fmin     = 400 # Arbitrary minimum
                self.fmax     = 500
                self.fmaxtrue = self.fmax
                self.nmin     = 2
                self.nmax     = 7
                self.nmin0    = self.nmin
                self.nmax0    = self.nmax
                self.ntol     = 0.01 / (10 * 365 * 24 * 3600) # Ten years, arbitrary
                self.taumin   = 40  * 365 * 24 * 3600       # Forty years
                self.taumax   = 400 * 365 * 24 * 3600       # Four hundred years, arbitrary
                self.ktol     = 0.01 / (10 * 365 * 24 * 3600) # Ten years, arbitrary
                self.mismatch = 0.2
                self.dur      = 20 * 86400 # Twenty days, arbitrary
                self.maxtemps = maxtemps

                knots = np.linspace(0, self.dur, 30)

                self.knots   = knots
                bf.knotslist = knots

                # ek.allidealisedknots(self.s, dur, 40, self.fmaxtrue, self.nmax, self.taumin, self.mismatch)
                # knotnum = ek.getknotnum(self.s, dur, self.fmaxtrue, self.nmax, self.taumin, self.mismatch)
                # self.knots = bf.knotslist[:knotnum + 1]

                # Knots built from the algorithm are:
                # [0, 2883584.55, 6545862.21016838]

                self.maxtemps = maxtemps

        def SetDefaultCassA(self, maxtemps=1000000):
                self.s        = 3
                self.fmin     = 100 # Arbitrary minimum
                self.fmax     = 300
                self.fmaxtrue = self.fmax
                self.nmin     = 2
                self.nmax     = 7
                self.nmin0    = self.nmin
                self.nmax0    = self.nmax
                self.ntol     = 0.01 / (10 * 365 * 24 * 3600) # Ten years, arbitrary
                self.taumin   = 300  * 365 * 24 * 3600      # Three hundred years
                self.taumax   = 3000 * 365 * 24 * 3600      # Three thousand years, arbitrary
                self.ktol     = 0.01 / (10 * 365 * 24 * 3600) # Ten years, arbitrary
                self.mismatch = 0.2
                self.dur      = 100 * 86400 # 100 days, arbitrary
                self.maxtemps = maxtemps

                ek.allidealisedknots(self.s, self.dur, 40, self.fmaxtrue, self.nmax, self.taumin, self.mismatch)
                knotnum = ek.getknotnum(self.s, self.dur, self.fmaxtrue, self.nmax, self.taumin, self.mismatch)
                knots = bf.knotslist[:knotnum + 1]

                self.knots = knots
                bf.knotslist = knots

                # Knots for this object are:
                # [0, 11534336.55]

                self.maxtemps = maxtemps

        # Approx 2.5 * 10^5 templates
        def SetSmallTestCase(self, maxtemps=1000000):
                self.s        = 3
                self.fmin     = 9 # Arbitrary minimum
                self.fmax     = 11
                self.fmaxtrue = self.fmax
                self.nmin     = 2
                self.nmax     = 2.5
                self.nmin0    = self.nmin
                self.nmax0    = self.nmax
                self.ntol     = 0.01 / (10 * 365 * 24 * 3600) # Ten years, arbitrary
                self.taumin   = 365 * 86400                 # Quarter of a day (Six hours)
                self.taumax   = 10 * 365 * 86400            # Half a day (12 hours)
                self.ktol     = 0.01 / (10 * 365 * 24 * 3600) # Ten years, arbitrary
                self.knots    = np.array([0, 86400])
                self.dur      = 86400
                self.mismatch = 0.4
                self.maxtemps = maxtemps

                bf.knotslist = self.knots

        # Approx 1.6 * 10^5 templates
        def SetSmallMultiSegCase(self, maxtemps=1000000):
                self.s        = 3
                self.fmin     = 9 # Arbitrary minimum
                self.fmax     = 11
                self.fmaxtrue = self.fmax
                self.nmin     = 2
                self.nmax     = 2.5
                self.nmin0    = self.nmin
                self.nmax0    = self.nmax
                self.ntol     = 0.01 / (10 * 365 * 24 * 3600) # Ten years, arbitrary
                self.taumin   = 365 * 86400                 # Quarter of a day (Six hours)
                self.taumax   = 10 * 365 * 86400            # Half a day (12 hours)
                self.ktol     = 0.01 / (10 * 365 * 24 * 3600) # Ten years, arbitrary
                self.knots    = np.array([0, 43200, 86400])
                self.dur      = 86400
                self.mismatch = 0.4
                self.maxtemps = maxtemps

                bf.knotslist = self.knots

        # Sets the knots of a TBank object to be those calculated by our knot algorithm. A duration or knot number must first be set
        def SetKnotsByAlg(self, knotnum=None):
                s        = self.s
                steps    = 40   # These are the number of steps used to find the maximum error in the bi-section method for calculating knots. A value of at least 30 usually gives good enough results
                f0       = self.fmaxtrue
                nmax     = self.nmax
                tau      = self.taumin
                mismatch = self.mismatch
                dur      = self.dur

                self.knots = ek.allidealisedknots(s, dur, steps, f0, nmax, tau, mismatch, knotnum=knotnum)
                bf.knotslist = self.knots

        # Returns a string containing all parameters for this particular TBank, excluding the knots
        def toString(self):
                s        = self.s
                fmin     = self.fmin
                fmax     = self.fmax
                fmaxtrue = self.fmaxtrue
                nmin     = self.nmin
                nmax     = self.nmax
                nmin0    = self.nmin0
                nmax0    = self.nmax0
                ntol     = self.ntol
                taumin   = self.taumin
                taumax   = self.taumax
                ktol     = self.ktol
                mismatch = self.mismatch
                dur      = self.dur
                maxtemps = self.maxtemps

                return str([s, fmin, fmax, fmaxtrue, nmin, nmax, nmin0, nmax0, ntol, taumin, taumax, ktol, mismatch, dur, maxtemps])

class FInput:
        Tsft          = 60
        h0            = 1e-24
        cosi          = 0.123
        psi           = 2.345
        phi0          = 3.210

        dt_wf         = 10
        Alpha         = 6.12 #3.45
        Delta         = 1.02 #-0.42
        detector      = 'H1'
        assume_sqrtSh = 1e-23

        dfreq         = 0.1 / (6 * 3600) # Arbitrary, default BNSR signal duration

        sourceDeltaT  = 2
