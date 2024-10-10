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
Common classes for piecewise model module.
"""

from . import basis_functions as bf
from . import estimating_knots as ek
from . import gte_and_other_methods as gom

HOUR = 3600
DAY = 24 * HOUR


class TBank:
    name = None
    s = None
    fmin = None
    fmax = None
    nmin = None
    nmax = None
    Izz = None
    ellip = None
    radius = None
    kmin = None
    kmax = None
    mismatch = None
    dur = None
    knotnum = None
    knots = None

    flags_bbox = []
    flags_intbox = []

    tstart = None

    cosi = None
    psi = None
    phi0 = None
    Alpha = None
    Delta = None

    detectors = []

    noise_path = None
    sft_path = None
    psd_path = None

    def SetTBankParams(self, args):

        if args.s:
            self.s = args.s
        if args.fmin:
            self.fmin = args.fmin
        if args.fmax:
            self.fmax = args.fmax
        if args.nmin:
            self.nmin = args.nmin
        if args.nmax:
            self.nmax = args.nmax
        if args.Izz:
            self.Izz = args.Izz
        if args.ellip:
            self.ellip = args.ellip
        if args.radius:
            self.radius = args.radius

        self.kmin = (
            args.kmin
            if args.kmin is not None
            else gom.kforGWsource(Izz=self.Izz, ellip=self.ellip, radius=self.radius)
            * 0.1
        )
        self.kmax = (
            args.kmax
            if args.kmax is not None
            else gom.kforGWsource(Izz=self.Izz, ellip=self.ellip, radius=self.radius)
        )

        if args.dur:
            self.dur = args.dur
        if args.maxmismatch:
            self.maxmismatch = args.maxmismatch
        if args.knotnum:
            self.knotnum = args.knotnum

        if args.flags_bbox is not []:
            self.flags_bbox = args.flags_bbox
        if args.flags_intbox is not []:
            self.flags_intbox = args.flags_intbox

        if args.tstart:
            self.tstart = args.tstart

        if args.cosi:
            self.cosi = args.cosi
        if args.psi:
            self.psi = args.psi
        if args.phi0:
            self.phi0 = args.phi0
        if args.Alpha:
            self.Alpha = args.Alpha
        if args.Delta:
            self.Delta = args.Delta

        if not (args.detectors == []):
            self.detectors = args.detectors

        if args.noise_path:
            self.noise_path = args.noise_path
        if args.sft_path:
            self.sft_path = args.sft_path
        if args.psd_path:
            self.psd_path = args.psd_path

    def SetDefaultGW170817(self):
        self.name = "GW170817"

        # Physical Parameters
        self.s = 2
        self.fmin = 999
        self.fmax = 1000
        self.nmin = 3
        self.nmax = 5
        self.Izz = 1e38  # The long transient paper uses a value of 4.34e38
        self.ellip = 1e-4
        self.radius = 1e4

        self.kmin = (
            gom.kforGWsource(Izz=self.Izz, ellip=self.ellip, radius=self.radius) * 0.1
        )
        self.kmax = gom.kforGWsource(Izz=self.Izz, ellip=self.ellip, radius=self.radius)

        # Relevant to search data
        self.maxmismatch = 0.2
        self.Tsft = 10

        # Timing
        self.dur = 1800
        self.knotnum = 0
        self.knots = [0, 1800]
        self.tstart = 1187008882

        # Sky Positions and orientations
        self.cosi = None  # 0.882948           # iota = 28 Degrees
        self.psi = None  # 0
        self.phi0 = None  # 0
        self.Alpha = 3.446
        self.Delta = -0.408

        self.detectors = ["H1", "L1"]

    def SetDefaultGW190425(self):
        self.name = "GW190425"

        # Physical Parameters
        self.s = 2
        self.fmin = 999
        self.fmax = 1000
        self.nmin = 3
        self.nmax = 5
        self.Izz = 1e38  # The long transient paper uses a value of 4.34e38
        self.ellip = 1e-4
        self.radius = 1e4

        self.kmin = (
            gom.kforGWsource(Izz=self.Izz, ellip=self.ellip, radius=self.radius) * 0.1
        )
        self.kmax = gom.kforGWsource(Izz=self.Izz, ellip=self.ellip, radius=self.radius)

        # Relevant to search data
        self.maxmismatch = 0.2
        self.Tsft = 10

        # Timing
        self.dur = 1800
        self.knotnum = 0
        self.knots = [0, 1800]
        self.tstart = 1240215503

        # Sky Positions and orientations
        self.cosi = None  # 0
        self.psi = None  # 0
        self.phi0 = None  # 0
        self.Alpha = 4.218
        self.Delta = 0.401

        self.detectors = ["L1", "V1"]

    def SetDefault1987A(self):
        self.name = "1987A"

        self.s = 2
        self.fmin = 49.99
        self.fmax = 50
        self.nmin = 3
        self.nmax = 5
        self.kmin = gom.kforGWsource() * 0.1
        self.kmax = gom.kforGWsource()

        self.maxmismatch = 0.2
        self.Tsft = 1800

        self.dur = 20 * DAY
        self.knotnum = 0
        self.knots = [0, self.dur]

    def SetDefaultCasA(self):
        self.name = "CasA"

        self.s = 2
        self.fmin = 29.99
        self.fmax = 30
        self.nmin = 3
        self.nmax = 5
        self.kmin = gom.kforGWsource() * 0.1
        self.kmax = gom.kforGWsource()

        self.maxmismatch = 0.2
        self.Tsft = 1800

        self.dur = 20 * DAY
        self.knotnum = 0
        self.knots = [0, self.dur]

    # Sets the knots of a TBank object to be those calculated by our knot algorithm. A duration or knot number must first be set
    def SetKnotsByAlg(self):
        s = self.s
        steps = 40  # These are the number of steps used to find the maximum error in the bi-section method for calculating knots. A value of at least 30 usually gives good enough results
        f0 = self.fmax
        nmax = self.nmax
        kmax = self.kmax
        mismatch = self.maxmismatch
        dur = self.dur
        knotnum = self.knotnum

        self.knots = ek.allidealisedknots(
            s, dur, steps, f0, nmax, kmax, mismatch, knotnum=knotnum
        )

        if self.knotnum != 0:
            self.knots = self.knots[0 : self.knotnum + 1]
            self.dur = self.knots[-1] - self.knots[0]
            self.knotnum = len(self.knots)

        bf.knotslist = self.knots

    # Returns a string containing all parameters for this particular TBank, excluding the knots
    def toString(self):
        string = self.name

        tbank = TBank()

        if string == "GW170817":
            tbank.SetDefaultGW170817()
        elif string == "GW190425":
            tbank.SetDefaultGW190425()

        tbank.fmin = -1
        tbank.fmax = -1

        original_dict = vars(self)
        default_dict = vars(tbank)

        different_variables = []

        skip_variables = [
            "name",
            "knots",
            "flags_bbox",
            "flags_intbox",
            "noise_path",
            "sft_path",
            "psd_path",
        ]

        for key, value in original_dict.items():

            if key in skip_variables:
                continue

            if value is None:
                continue

            default_variable = default_dict[key]

            if value != default_variable:

                if key == "detectors":

                    value_string = ""

                    for det in value:
                        value_string += det
                        different_variables.append(value_string)
                else:
                    value_string = key + "-" + str(value)
                    different_variables.append(value_string)

        for diff_var in different_variables:
            string += "_" + diff_var

        return string


class FInput:
    pass
