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
Functions used in simulation studies of the piecewise model.
"""

import logging
import os

import numpy as np

import lal
import lalpulsar as lp
from lalpulsar import simulateCW

from . import basis_functions as bf
from . import mols_for_gte as mols
from . import semicoherent_metric_methods as scmm
from . import tbank_estimates as tbe


def waveform(h0, params, cosi, tref):

    s = len(params[0][0])
    coeffs = bf.allcoeffs(s)
    f0 = params[0][0][0]

    def wf(dt):

        freq = mols.modelvalueatpoint(dt + tref, coeffs, params, ignoreintcheck=True)

        h_t = h0 * (freq / f0) ** 2

        ap = h_t * (1.0 + cosi**2) / 2.0
        ax = h_t * cosi

        dphi = mols.phase(dt + tref, coeffs, params, ignoreintcheck=True)

        return dphi, ap, ax

    return wf


# Builds and saves SFTs to files where a signal of strain h0 is randomly chosen and injected. Tdata is the length of data
# to be used. tbank contains all the parameters which defines the parameter space. tstart is the start time of the data
# trefsegfrac tells us how far through the data tref is defined. E.g. If set to 0, tref occurs at the start of the data,
# if 0.5 tref occurs half way through the data and so on. Tsft is the length of SFTs to be built. parentdirectory is the
# path leading to the directory where the SFTs will be saved. h0 is the strain of the injected signal and noiseratio is
# what fraction weaker the noise is to the signal. E.g. if set to 100, then the noise will be at 1/100 of h0.
def buildSFTs(
    Tdata, tbank, finputdata, trefsegfrac=0.0, parentdirectory=".", rand_seed=0
):

    # Double check the start time of the knots agrees with tstart and set them to be equal if they don't
    if bf.knotslist[0] != finputdata.tstart:
        for i in range(len(bf.knotslist)):
            bf.knotslist[i] = bf.knotslist[i] + finputdata.tstart

    logging.info("Knots: " + str(bf.knotslist))

    # Set Bounds, metric, mismatch and lattice type
    finalknot = len(bf.knotslist)
    s = tbank.s
    mismatch = tbank.maxmismatch

    tiling = lp.CreateLatticeTiling(s * finalknot)

    metric = scmm.PreCompMetric(s)

    logging.info("Knots before bounds set: " + str(bf.knotslist))
    tbe.setbounds(tiling, tbank)

    logging.info("Knots after bounds set: " + str(bf.knotslist))
    lp.SetTilingLatticeAndMetric(tiling, lp.TILING_LATTICE_ANSTAR, metric, mismatch)

    # Create random parameters to be used as an injected signal
    randparams = lal.CreateRandomParams(rand_seed)
    signalparams = lal.gsl_matrix(s * finalknot, 1)
    lp.RandomLatticeTilingPoints(tiling, 0, randparams, signalparams)

    f0 = signalparams.data[0][0]
    f1 = signalparams.data[1][0]
    f2 = signalparams.data[2][0]

    logging.info(
        "Randomly generated Signal Params: %s", str(np.transpose(signalparams.data)[0])
    )

    # Waveform and simulateCW object
    tref = bf.knotslist[0] + Tdata * trefsegfrac
    params = mols.solsbyint(np.transpose(signalparams.data)[0], s)
    wf = waveform(finputdata.h0, params, finputdata.cosi, tref)

    for i, detector in enumerate(finputdata.detectors.data):

        S = simulateCW.CWSimulator(
            tref,
            bf.knotslist[0],
            Tdata,
            wf,
            finputdata.dt_wf,
            finputdata.phi0,
            finputdata.psi,
            finputdata.Alpha,
            finputdata.Delta,
            detector,
            tref_at_det=True,
        )

        # Create directory to save SFTs into using strain and signal parameters as directory name
        # directoryname = os.path.join(parentdirectory, f"SFTs_h0-{finputdata.h0:.2e}_f0-{f0:.3f}_f1-{f1:.2e}_f2-{f2:.2e}_dur-{tbank.dur}_tstart-{finputdata.tstart}/")

        path_to_detector_sfts = os.path.join(
            parentdirectory,
            f"SFTs_h0-{finputdata.h0:.2e}_f0-{f0:.3f}_f1-{f1:.2e}_f2-{f2:.2e}_dur-{tbank.dur}_tstart-{finputdata.tstart}/",
        )
        directoryname = os.path.join(path_to_detector_sfts, detector + "/")

        if not os.path.isdir(path_to_detector_sfts):
            os.mkdir(path_to_detector_sfts)

        if not os.path.isdir(directoryname):
            os.mkdir(directoryname)

        # Write SFT files

        # Rule of thumb given by Karl for fmax -> fmax = f0 + 10^-4 * f0 + (50 + 8)/Tsft
        # The 10^-4 comes from an approximation doppler modulation, 50 + 8 is the number of bins extra either side we should have
        for file, i, N in S.write_sft_files(
            fmax=tbank.fmax + 20,
            Tsft=finputdata.Tsft,
            noise_sqrt_Sh=finputdata.noise_sqrt_Sh[i],
            comment="PWModPhase",
            out_dir=directoryname,
        ):
            logging.info("Generated SFT file %s (%i of %i)" % (file, i + 1, N))

    # Return signal parameters so they can be known outside of this method

    return np.transpose(signalparams.data)[0], path_to_detector_sfts
