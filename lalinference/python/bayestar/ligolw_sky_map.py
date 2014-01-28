#
# Copyright (C) 2013  Leo Singer
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
#
"""
Convenience function to produce a sky map from LIGO-LW rows.
"""
from __future__ import division
__author__ = "Leo Singer <leo.singer@ligo.org>"


import time
import numpy as np
from . import filter
from . import postprocess
from . import timing
from . import sky_map
import lal, lalsimulation


def ligolw_sky_map(sngl_inspirals, approximant, amplitude_order, phase_order, f_low, min_distance=None, max_distance=None, prior_distance_power=None, method="toa_phoa_snr", reference_frequency=None, psds=None, nside=-1, chain_dump=None):
    """Convenience function to produce a sky map from LIGO-LW rows. Note that
    min_distance and max_distance should be in Mpc."""

    ifos = [sngl_inspiral.ifo for sngl_inspiral in sngl_inspirals]

    # Extract masses from the table.
    mass1s = np.asarray([sngl_inspiral.mass1 for sngl_inspiral in sngl_inspirals])
    mass2s = np.asarray([sngl_inspiral.mass2 for sngl_inspiral in sngl_inspirals])

    # Extract SNRs from table.
    # FIXME: should get complex SNR, but MBTAOnline events don't populate the
    # coa_phase column, and we are not using coa_phase yet.
    snrs = np.asarray([sngl_inspiral.snr for sngl_inspiral in sngl_inspirals])

    # Extract TOAs from table.
    toas_ns = np.asarray([sngl_inspiral.get_end().ns()
        for sngl_inspiral in sngl_inspirals], dtype=np.int64)

    # Optionally apply reference frequency shift.
    if reference_frequency is not None:
        toas_ns -= [int(round(1e9 * lalsimulation.
            SimInspiralTaylorF2ReducedSpinChirpTime(
            reference_frequency,
            m1 * lal.LAL_MSUN_SI,
            m2 * lal.LAL_MSUN_SI,
            0, 4))) for m1, m2 in zip(mass1s, mass2s)]

    # Find average Greenwich mean sidereal time of event.
    mean_toa_ns = sum(toas_ns) // len(toas_ns)
    epoch = lal.LIGOTimeGPS(0, long(mean_toa_ns))
    gmst = lal.GreenwichMeanSiderealTime(epoch)

    # Convert TOAs from nanoseconds to seconds, subtracting off mean TOA
    # to keep numbers small.
    toas = 1e-9 * (toas_ns - mean_toa_ns)

    # Retrieve phases on arrival from table.
    phoas = np.asarray([sngl_inspiral.coa_phase for sngl_inspiral in sngl_inspirals])

    # Power spectra for each detector.
    if psds is None:
        psds = [timing.get_noise_psd_func(ifo) for ifo in ifos]

    # Signal models for each detector.
    signal_models = [timing.SignalModel(mass1, mass2, psd, f_low, approximant, amplitude_order, phase_order)
        for mass1, mass2, psd in zip(mass1s, mass2s, psds)]

    # Get SNR=1 horizon distances for each detector.
    horizons = [signal_model.get_horizon_distance()
        for signal_model in signal_models]

    # Estimate TOA uncertainty (squared) using CRB or BRB evaluated at MEASURED
    # values of the SNRs.
    w_toas = [1/np.square(signal_model.get_toa_uncert(np.abs(snr)))
        for signal_model, snr in zip(signal_models, snrs)]

    w1s = [signal_model.get_sn_moment(1) for signal_model in signal_models]
    w2s = [signal_model.get_sn_moment(2) for signal_model in signal_models]

    # Look up physical parameters for detector.
    detectors = [lalsimulation.DetectorPrefixToLALDetector(str(ifo))
        for ifo in ifos]
    responses = [det.response for det in detectors]
    locations = [det.location for det in detectors]

    # If minimum distance is not specified, then default to 0 Mpc.
    if min_distance is None:
        min_distance = 0

    # If maximum distance is not specified, then default to the SNR=4
    # horizon distance of the most sensitive detector.
    if max_distance is None:
        max_distance = max(horizons) / 4

    # If prior_distance_power is not specified, then default to 2
    # (p(r) ~ r^2, uniform in volume).
    if prior_distance_power is None:
        prior_distance_power = 2

    # Raise an exception if 0 Mpc is the minimum effective distance and the prior
    # is of the form r**k for k<0
    if min_distance == 0 and prior_distance_power < 0:
        raise ValueError(("Prior is a power law r^k with k={}, "
            + "undefined at min_distance=0").format(prior_distance_power))

    # Time and run sky localization.
    start_time = time.time()
    if method == "toa":
        prob = sky_map.toa(gmst, toas, w_toas, locations, nside=nside)
    elif method == "toa_snr":
        prob = sky_map.toa_snr(gmst, toas, snrs, w_toas, responses, locations, horizons, min_distance, max_distance, prior_distance_power, nside=nside)
    elif method == "toa_phoa_snr":
        prob = sky_map.toa_phoa_snr(gmst, toas, phoas, snrs, w_toas, w1s, w2s, responses, locations, horizons, min_distance, max_distance, prior_distance_power, nside=nside)
    elif method == "toa_snr_mcmc":
        import emcee

        ntemps = 20
        nwalkers = 100
        ndim = 5
        sampler = emcee.PTSampler(
            ntemps=ntemps,
            nwalkers=nwalkers,
            dim=ndim,
            logl=(lambda args: sky_map.log_posterior_toa_snr(*args,
                gmst=gmst,
                toas=toas,
                snrs=snrs,
                w_toas=w_toas,
                responses=responses,
                locations=locations,
                horizons=horizons,
                prior_distance_power=prior_distance_power)),
            logp=(lambda (ra, sin_dec, distance, u, twopsi):
                1 if 0 <= ra < 2*np.pi
                and -1 <= sin_dec <= 1
                and min_distance <= distance <= max_distance
                and 0 <= u <= 1
                and 0 <= twopsi < 2*np.pi
                else -np.inf))
        p0 = np.random.uniform(
            [0, -1, min_distance, 0, 0],
            [2*np.pi, 1, max_distance, 1, 2*np.pi], (ntemps, nwalkers, ndim))
        sampler.run_mcmc(p0, 1000)
        if chain_dump is not None:
            np.save(chain_dump, sampler.chain)
        ra, sin_dec, _, _, _ = np.concatenate(sampler.chain[0, :, 100:]).T
        theta = np.arccos(sin_dec)
        phi = ra
        prob = postprocess.adaptive_healpix_histogram(theta, phi)
    else:
        raise ValueError("Unrecognized method: %s" % method)
    end_time = time.time()

    # Find elapsed run time.
    elapsed_time = end_time - start_time

    # Done!
    return prob, epoch, elapsed_time


def gracedb_sky_map(coinc_file, psd_file, waveform, f_low, min_distance=None, max_distance=None, prior_distance_power=None, reference_frequency=None, nside=-1):
    # LIGO-LW XML imports.
    from glue.ligolw import table as ligolw_table
    from glue.ligolw import utils as ligolw_utils
    from glue.ligolw import lsctables
    import lal.series

    # Determine approximant, amplitude order, and phase order from command line arguments.
    approximant, amplitude_order, phase_order = timing.get_approximant_and_orders_from_string(waveform)

    # Read input file.
    xmldoc, _ = ligolw_utils.load_fileobj(coinc_file)

    # Locate the tables that we need.
    coinc_inspiral_table = ligolw_table.get_table(xmldoc,
        lsctables.CoincInspiralTable.tableName)
    coinc_map_table = ligolw_table.get_table(xmldoc,
        lsctables.CoincMapTable.tableName)
    sngl_inspiral_table = ligolw_table.get_table(xmldoc,
        lsctables.SnglInspiralTable.tableName)

    # Locate the sngl_inspiral rows that we need.
    coinc_inspiral = coinc_inspiral_table[0]
    coinc_event_id = coinc_inspiral.coinc_event_id
    event_ids = [coinc_map.event_id for coinc_map in coinc_map_table
        if coinc_map.coinc_event_id == coinc_event_id]
    sngl_inspirals = [(sngl_inspiral for sngl_inspiral in sngl_inspiral_table
        if sngl_inspiral.event_id == event_id).next() for event_id in event_ids]
    instruments = set(sngl_inspiral.ifo for sngl_inspiral in sngl_inspirals)

    # Read PSDs.
    if psd_file is None:
        psds = None
    else:
        xmldoc, _ = ligolw_utils.load_fileobj(psd_file)
        psds = lal.series.read_psd_xmldoc(xmldoc)

        # Rearrange PSDs into the same order as the sngl_inspirals.
        psds = [psds[sngl_inspiral.ifo] for sngl_inspiral in sngl_inspirals]

        # Interpolate PSDs.
        psds = [timing.InterpolatedPSD(filter.abscissa(psd), psd.data.data) for psd in psds]

    # TOA+SNR sky localization
    prob, epoch, elapsed_time = ligolw_sky_map(sngl_inspirals, approximant,
        amplitude_order, phase_order, f_low,
        min_distance, max_distance, prior_distance_power,
        reference_frequency=reference_frequency, nside=nside, psds=psds)

    return prob, epoch, elapsed_time, instruments
