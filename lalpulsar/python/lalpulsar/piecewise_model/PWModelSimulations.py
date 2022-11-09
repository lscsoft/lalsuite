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

import lal
import lalpulsar as lp
from lalpulsar import simulateCW

from . import BasisFunctions as bf
from . import MOLSforGTE as mols
from . import SemicoherentMetricMethods as scmm
from . import TBankEstimates as tbe

import numpy as np
import os
import logging

#dts = []

def waveform(h0, cosi, params, tref):
        s = len(params[0][0])
        coeffs = bf.allcoeffs(s)

        def wf(dt):

                ap = h0 * (1.0 + cosi ** 2) / 2.0
                ax = h0 * cosi

                dphi = mols.phase(dt + tref, coeffs, params, ignoreintcheck=True)

                #freq  = params[0][0][0]
                #f1dot = params[0][0][1]
                #f2dot = params[0][0][2]
                #dphi = lal.TWOPI * (freq * dt + f1dot * 0.5 * dt**2 + f2dot * 1/6 * dt**3)

                return dphi, ap, ax
        return wf

# Deciding to just keep these variables here, not in the buildSFTs method, have no reason to change them just yet.
cosi     = 0.123
psi      = 2.345
phi0     = 3.210

dt_wf    = 10
Alpha    = 6.12 #3.45
Delta    = 1.02 #-0.42
detector = 'H1'

# Builds and saves SFTs to files where a signal of strain h0 is randomly chosen and injected. Tdata is the length of data
# to be used. tbank contains all the parameters which defines the parameter space. tstart is the start time of the data
# trefsegfrac tells us how far through the data tref is defined. E.g. If set to 0, tref occurs at the start of the data,
# if 0.5 tref occurs half way through the data and so on. Tsft is the length of SFTs to be built. parentdirectory is the
# path leading to the directory where the SFTs will be saved. h0 is the strain of the injected signal and noiseratio is
# what fraction weaker the noise is to the signal. E.g. if set to 100, then the noise will be at 1/100 of h0.
def buildSFTs(Tdata, tbank, tstart=900000000, trefsegfrac=0., Tsft=60, parentdirectory=".", h0=1e-24, noiseratio=100):

        # Double check the start time of the knots agrees with tstart and set them to be equal if they don't
        if bf.knotslist[0] != tstart:
                for i in range(len(bf.knotslist)):
                        bf.knotslist[i] = bf.knotslist[i] + tstart

        logging.info("Knots: " + str(bf.knotslist))

        # Set Bounds, metric, mismatch and lattice type
        finalknot = len(bf.knotslist)
        s = tbank.s
        mismatch = tbank.maxmismatch

        tiling = lp.CreateLatticeTiling(s * finalknot)

        metric = scmm.metric(s)
        print("Knots before bounds set: " + str(bf.knotslist))
        tbe.setbounds(tiling, tbank)
        print("Knots after bounds set: " + str(bf.knotslist))
        lp.SetTilingLatticeAndMetric(tiling, lp.TILING_LATTICE_ANSTAR, metric, mismatch)

        # Create random parameters to be used as an injected signal
        randparams = lal.CreateRandomParams(0)
        signalparams = lal.gsl_matrix(s * finalknot, 1)
        lp.RandomLatticeTilingPoints(tiling, 0, randparams, signalparams);

        params = mols.solsbyint(np.transpose(signalparams.data)[0], tbank.s)

        logging.info("Randomly generated Signal Params: %s", str(np.transpose(signalparams.data)[0]))

        # Waveform and simulateCW object
        tref = bf.knotslist[0] + Tdata * trefsegfrac
        wf = waveform(h0, cosi, params, tref)
        print("Tref in pwsim is: " + str(tref))
        S = simulateCW.CWSimulator(tref, bf.knotslist[0], Tdata, wf, dt_wf, phi0, psi, Alpha, Delta, detector, tref_at_det=True)

        # Create directory to save SFTs into using strain and signal parameters as directory name
        f0 = params[0][0][0]
        f1 = params[0][0][1]
        f2 = params[0][0][2]
        directoryname = os.path.join(parentdirectory, f"SFTs_h0-{h0:.2e}_f0-{f0:.3f}_f1-{f1:.2e}_f2-{f2:.2e}_dur-{tbank.dur}_tstart-{tstart}")

        if not os.path.isdir(directoryname):
                os.mkdir(directoryname)

        # Write SFT files

        # Rule of thumb given by Karl for fmax -> fmax = f0 + 10^-4 * f0 + (50 + 8)/Tsft
        # The 10^-4 comes from an approximation doppler modulation, 50 + 8 is the number of bins extra either side we should have
        for file, i, N in S.write_sft_files(fmax=tbank.fmax + 20, Tsft=Tsft, noise_sqrt_Sh=h0 / noiseratio, comment="PWModPhase", out_dir=directoryname):
                print('Generated SFT file %s (%i of %i)' % (file, i+1, N))

        # Return signal parameters so they can be known outside of this method

        #return [float(elem) for elem in paramsplitstr]
        return np.transpose(signalparams.data)[0], directoryname
