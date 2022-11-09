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
from . import GTEandOtherMethods as gom

import numpy as np

HOUR = 3600
DAY = 24 * HOUR
YEAR = 365.25 * DAY

class TBank:
        s          = None
        fmin       = None
        fmax       = None
        nmin       = None
        nmax       = None
        kmin       = None
        kmax       = None
        mismatch   = None
        dur        = None
        knotnum     = None

        flagoption = 0
        temps = 0
        currenttemplate = 0

        def SetTBankParams(self, args):
                if args.s             is not None: self.s             = args.s
                if args.fmin          is not None: self.fmin          = args.fmin
                if args.fmax          is not None: self.fmax          = args.fmax
                if args.nmin          is not None: self.nmin          = args.nmin
                if args.nmax          is not None: self.nmax          = args.nmax
                if args.kmin          is not None: self.kmin          = args.kmin
                if args.kmax          is not None: self.kmax          = args.kmax
                if args.dur           is not None: self.dur           = args.dur
                if args.maxmismatch   is not None: self.maxmismatch   = args.maxmismatch

        def SetDefaultBNSR(self):
                self.s             = 2
                self.fmin          = 999
                self.fmax          = 1000
                self.nmin          = 3
                self.nmax          = 5
                self.kmin          = gom.kforGWsource() * 0.1
                self.kmax          = gom.kforGWsource()

                self.maxmismatch   = 0.2
                self.Tsft          = 10

                self.dur           = 1800

        def SetDefault1987A(self):
                self.s             = 2
                self.fmin          = 49.99
                self.fmax          = 50
                self.nmin          = 3
                self.nmax          = 5
                self.kmin          = gom.kforGWsource() * 0.1
                self.kmax          = gom.kforGWsource()

                self.maxmismatch   = 0.2
                self.Tsft          = 1800

                self.dur           = 20 * DAY

        def SetDefaultCasA(self):
                self.s             = 2
                self.fmin          = 29.99
                self.fmax          = 30
                self.nmin          = 3
                self.nmax          = 5
                self.kmin          = gom.kforGWsource() * 0.1
                self.kmax          = gom.kforGWsource()

                self.maxmismatch   = 0.2
                self.Tsft          = 1800

                self.dur           = 20 * DAY

        # Sets the knots of a TBank object to be those calculated by our knot algorithm. A duration or knot number must first be set
        def SetKnotsByAlg(self):
                s        = self.s
                steps    = 40   # These are the number of steps used to find the maximum error in the bi-section method for calculating knots. A value of at least 30 usually gives good enough results
                f0       = self.fmax
                nmax     = self.nmax
                kmax     = self.kmax
                mismatch = self.maxmismatch
                dur      = self.dur
                knotnum  = self.knotnum

                self.knots = ek.allidealisedknots(s, dur, steps, f0, nmax, kmax, mismatch, knotnum=knotnum)
                if self.knotnum:
                        self.knots = self.knots[0:self.knotnum]
                        self.dur = self.knots[-1] - self.knots[0]
                        self.knotnum = len(self.knots)

                bf.knotslist = self.knots

        # Returns a string containing all parameters for this particular TBank, excluding the knots
        def toString(self):
                string = 'TBank'

                for k, v in vars(self).items():
                        if k == 'knots':
                                continue
                        if v is not None:
                                string += f"_{k}-{v:.3g}"

                return string

class FInput:
        pass
