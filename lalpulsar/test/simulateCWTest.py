# Copyright (C) 2017 Karl Wette
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

import os
import sys
import shutil
import numpy as np

import lal
import lalpulsar
from lalpulsar import simulateCW as simCW

# load Earth and Sun ephemerides
earth_ephem_file = os.path.join(os.environ['LAL_TEST_PKGDATADIR'], 'earth00-19-DE405.dat.gz')
sun_ephem_file = os.path.join(os.environ['LAL_TEST_PKGDATADIR'], 'sun00-19-DE405.dat.gz')
ephemerides = lalpulsar.InitBarycenter(earth_ephem_file, sun_ephem_file)

# amplitude parameters
h0 = 1e-24
assume_sqrtSh = 1e-23
cosi = 0.123
psi = 2.345
phi0 = 3.210

# phase parameters
freq  = 10.0
fcmin =  5.0
fcmax = 15.0
f1dot = -1.35e-8
alpha = 6.12
delta = 1.02

# detector parameters
detector = 'V1'

# data parameters
tref   = 900043200
tstart = 900000000
Tdata  = 86400

# SFT parameters
fmax = 32.0
Tsft = 1800

# F-statistic parameters
dfreq = 0.1 / Tdata
Nfs = 111
fsmin = freq - (Nfs - 1) * 0.5 * dfreq

# compute F-statistic of reference spindown signal
def compute_Fstat_spindown_reference():

    # create fake SFT catalog
    sft_ts = lalpulsar.MakeTimestamps(tstart, Tdata, Tsft, 0)
    sft_catalog = lalpulsar.AddToFakeSFTCatalog(None, detector, sft_ts)

    # create default F-statistic optional arguments
    Fstat_opt_args = lalpulsar.FstatOptionalArgs(lalpulsar.FstatOptionalArgsDefaults)
    Fstat_opt_args.FstatMethod = lalpulsar.FMETHOD_DEMOD_BEST

    # create injection parameters
    Fstat_signal = lalpulsar.CreatePulsarParamsVector(1);
    Fstat_signal.data[0].Amp.h0 = h0
    Fstat_signal.data[0].Amp.cosi = cosi
    Fstat_signal.data[0].Amp.psi = psi
    Fstat_signal.data[0].Amp.phi0 = phi0
    Fstat_signal.data[0].Doppler.refTime = tref
    Fstat_signal.data[0].Doppler.Alpha = alpha
    Fstat_signal.data[0].Doppler.Delta = delta
    Fstat_signal.data[0].Doppler.fkdot[0] = freq
    Fstat_signal.data[0].Doppler.fkdot[1] = f1dot
    Fstat_opt_args.injectSources = Fstat_signal
    Fstat_assume_noise = lalpulsar.MultiNoiseFloor()
    Fstat_assume_noise.length = 1
    Fstat_assume_noise.sqrtSn[0] = assume_sqrtSh
    Fstat_opt_args.assumeSqrtSX = Fstat_assume_noise

    # create F-statistic input data
    Fstat_input = lalpulsar.CreateFstatInput(sft_catalog, fcmin, fcmax, dfreq, ephemerides, Fstat_opt_args)
    Fstat_res = 0
    doppler = lalpulsar.PulsarDopplerParams(Fstat_signal.data[0].Doppler)
    doppler.fkdot[0] = fsmin

    # search SFTs using F-statistic
    Fstat_res = lalpulsar.ComputeFstat(Fstat_res, Fstat_input, doppler, Nfs, lalpulsar.FSTATQ_2F)

    return Fstat_res.twoF

# compute F-statistic of spindown signal simulated by simulateCW
def compute_Fstat_spindown_simulateCW():

    # waveform: simple spindown model
    def waveform(h0, cosi, freq, f1dot):
        def wf(dt):
            dphi = lal.TWOPI * (freq * dt + f1dot * 0.5 * dt**2)
            ap = h0 * (1.0 + cosi**2) / 2.0
            ax = h0 * cosi
            return dphi, ap, ax
        return wf

    # waveform parameters
    dt_wf = 5

    # create simulator
    S = simCW.CWSimulator(tref, tstart, Tdata, waveform(h0, cosi, freq, f1dot), dt_wf, phi0, psi, alpha, delta, detector, earth_ephem_file=earth_ephem_file, sun_ephem_file=sun_ephem_file)

    # write SFTs
    for path, i, N in S.write_sft_files(fmax, Tsft, 'simulateCWTest'):
        pass

    # load SFTs from catalog
    sft_catalog = lalpulsar.SFTdataFind('V-1_V1_1800SFT_simCW_simulateCWTest-*.sft', None)

    # create default F-statistic optional arguments
    Fstat_opt_args = lalpulsar.FstatOptionalArgs(lalpulsar.FstatOptionalArgsDefaults)
    Fstat_opt_args.FstatMethod = lalpulsar.FMETHOD_DEMOD_BEST
    Fstat_assume_noise = lalpulsar.MultiNoiseFloor()
    Fstat_assume_noise.length = 1
    Fstat_assume_noise.sqrtSn[0] = assume_sqrtSh
    Fstat_opt_args.assumeSqrtSX = Fstat_assume_noise

    # create F-statistic input data
    Fstat_input = lalpulsar.CreateFstatInput(sft_catalog, fcmin, fcmax, dfreq, ephemerides, Fstat_opt_args)
    Fstat_res = 0
    doppler = lalpulsar.PulsarDopplerParams()
    doppler.refTime = tref
    doppler.Alpha = alpha
    doppler.Delta = delta
    doppler.fkdot[0] = fsmin
    doppler.fkdot[1] = f1dot

    # search SFTs using F-statistic
    Fstat_res = lalpulsar.ComputeFstat(Fstat_res, Fstat_input, doppler, Nfs, lalpulsar.FSTATQ_2F)

    return Fstat_res.twoF

# compute F-statistic of reference spindown signal
Fstat_ref = compute_Fstat_spindown_reference()
sys.stdout.write('F-statistic of reference spindown signal:\n   %s\n\n' % (', '.join(['%0.4f'] * len(Fstat_ref)) % tuple(Fstat_ref)))

# compute F-statistic of spindown signal from simulateCW
Fstat_simCW = compute_Fstat_spindown_simulateCW()
sys.stdout.write('F-statistic of spindown signal from simulateCW:\n   %s\n\n' % ('  '.join(['%0.4f'] * len(Fstat_simCW)) % tuple(Fstat_simCW)))

# compute root-mean-square error between F-statistics
Fstat_rmserr = np.sqrt(np.mean((1 - (Fstat_simCW / Fstat_ref))**2))
Fstat_rmserr_OK = (Fstat_rmserr < 1e-3)
sys.stdout.write('Root-mean-square error between F-statistics:\n   %0.2e => %s!\n\n' % (Fstat_rmserr, 'OK' if Fstat_rmserr_OK else 'ERROR'))

sys.exit(0 if Fstat_rmserr_OK else 1)
