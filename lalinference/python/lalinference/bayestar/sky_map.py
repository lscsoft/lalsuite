#
# Copyright (C) 2013-2017  Leo Singer
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


import itertools
import logging
import time
import numpy as np
import healpy as hp
from astropy.table import Table
from .decorator import with_numpy_random_seed
from . import distance
from . import ligolw
from . import filter
from . import postprocess
from . import timing
try:
    from . import _moc as moc
except ImportError:
    raise ImportError(
        'Could not import the lalinference.bayestar._moc Python C '
        'extension module. This probably means that LALInfernece was built '
        'without HEALPix support. Please install CHEALPix '
        '(https://sourceforge.net/projects/healpix/files/Healpix_3.30/'
        'chealpix-3.30.0.tar.gz), rebuild LALInference, and try again.')
try:
    from . import _sky_map
except ImportError:
    raise ImportError(
        'Could not import the lalinference.bayestar._sky_map Python C '
        'extension module. This probably means that LALInfernece was built '
        'without HEALPix support. Please install CHEALPix '
        '(https://sourceforge.net/projects/healpix/files/Healpix_3.30/'
        'chealpix-3.30.0.tar.gz), rebuild LALInference, and try again.')
import lal, lalsimulation
from glue.ligolw import table as ligolw_table
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import lsctables

log = logging.getLogger('BAYESTAR')


def toa_phoa_snr_log_prior(
        params, min_distance, max_distance, prior_distance_power, max_abs_t):
    ra, sin_dec, distance, u, twopsi, t = params
    return (
        prior_distance_power * np.log(distance)
        if 0 <= ra < 2*np.pi
        and -1 <= sin_dec <= 1
        and min_distance <= distance <= max_distance
        and -1 <= u <= 1
        and 0 <= twopsi < 2*np.pi
        and -max_abs_t <= t <= max_abs_t
        else -np.inf)


@with_numpy_random_seed
def emcee_sky_map(
        logl, loglargs, logp, logpargs, xmin, xmax,
        nside=-1, kde=False, chain_dump=None, max_horizon=1.0):
    # Set up sampler
    import emcee
    ntemps = 20
    nwalkers = 100
    nburnin = 1000
    nthin = 10
    niter = 10000 + nburnin
    ndim = len(xmin)
    sampler = emcee.PTSampler(
        ntemps=ntemps, nwalkers=nwalkers, dim=ndim, logl=logl, logp=logp,
        loglargs=loglargs, logpargs=logpargs)

    # Draw initial state from multivariate uniform distribution
    p0 = np.random.uniform(xmin, xmax, (ntemps, nwalkers, ndim))

    # Collect samples. The .copy() is important because PTSampler.sample()
    # reuses p on every iteration.
    chain = np.vstack([
        p[0, :, :].copy() for p, _, _
        in itertools.islice(
            sampler.sample(p0, iterations=niter, storechain=False),
            nburnin, niter, nthin
        )])

    # Extract polar coordinates. For all likelihoodds, the first two parameters
    # are ra, sin(dec).
    theta = np.arccos(chain[:, 1])
    phi = chain[:, 0]
    dist = chain[:, 2]

    # Do adaptive histogram binning if the user has not selected the KDE,
    # or if the user has selected the KDE but we need to guess the resolution.
    if nside == -1 or not kde:
        samples_per_bin = int(np.ceil(0.005 * len(chain)))
        prob = postprocess.adaptive_healpix_histogram(
            theta, phi, samples_per_bin, nside=nside, nest=True)

        # Determine what nside what actually used.
        nside = hp.npix2nside(len(prob))
        # Go one level finer, because the KDE will find more detail.
        nside *= 2

    if kde:
        from sky_area.sky_area_clustering import Clustered3DKDEPosterior
        ra = phi
        dec = 0.5 * np.pi - theta
        pts = np.column_stack((ra, dec, dist))
        # Pass a random subset of 1000 points to the KDE, to save time.
        pts = np.random.permutation(pts)[:1000, :]
        prob = Clustered3DKDEPosterior(pts).as_healpix(nside)

    # Optionally save posterior sample chain to file.
    # Read back in with np.load().
    if chain_dump:
        # Undo numerical conditioning of distances; convert back to Mpc
        chain[:, 2] *= max_horizon
        names = 'ra sin_dec distance cos_inclination twopsi time'.split()[:ndim]
        np.save(chain_dump, np.rec.fromrecords(chain, names=names))

    # Done!
    return prob


def ligolw_sky_map(
        sngl_inspirals, waveform, f_low,
        min_distance=None, max_distance=None, prior_distance_power=None,
        method="toa_phoa_snr", psds=None, nside=-1, chain_dump=None,
        phase_convention='antifindchirp', snr_series=None,
        enable_snr_series=False):
    """Convenience function to produce a sky map from LIGO-LW rows. Note that
    min_distance and max_distance should be in Mpc.

    Returns a 'NESTED' ordering HEALPix image as a Numpy array.
    """

    # Ensure that sngl_inspiral is either a single template or a list of
    # identical templates
    for key in 'mass1 mass2 spin1x spin1y spin1z spin2x spin2y spin2z'.split():
        if hasattr(sngl_inspirals[0], key):
            value = getattr(sngl_inspirals[0], key)
            if any(value != getattr(_, key) for _ in sngl_inspirals):
                raise ValueError(
                    '{0} field is not the same for all detectors'.format(key))

    ifos = [sngl_inspiral.ifo for sngl_inspiral in sngl_inspirals]

    # Extract SNRs from table.
    snrs = np.asarray([sngl_inspiral.snr
        for sngl_inspiral in sngl_inspirals])

    # Look up physical parameters for detector.
    detectors = [lalsimulation.DetectorPrefixToLALDetector(str(ifo))
        for ifo in ifos]
    responses = np.asarray([det.response for det in detectors])
    locations = np.asarray([det.location for det in detectors])

    # Power spectra for each detector.
    if psds is None:
        psds = [timing.get_noise_psd_func(ifo) for ifo in ifos]

    log.debug('calculating templates')
    H = filter.sngl_inspiral_psd(sngl_inspirals[0], waveform, f_min=f_low)

    log.debug('calculating noise PSDs')
    HS = [filter.signal_psd_series(H, S) for S in psds]

    # Signal models for each detector.
    log.debug('calculating Fisher matrix elements')
    signal_models = [timing.SignalModel(_) for _ in HS]

    # Get SNR=1 horizon distances for each detector.
    horizons = np.asarray([signal_model.get_horizon_distance()
        for signal_model in signal_models])

    weights = [1/np.square(signal_model.get_crb_toa_uncert(snr))
        for signal_model, snr in zip(signal_models, snrs)]

    # Center detector array.
    locations -= np.average(locations, weights=weights, axis=0)

    if enable_snr_series:
        log.warn('Enabling input of SNR time series. This feature is UNREVIEWED.')
    else:
        snr_series = None

    if snr_series is None:
        log.warn("No SNR time series found, so we are creating a zero-noise "
                 "SNR time series from the whitened template's autocorrelation "
                 "sequence. The sky localization uncertainty may be "
                 "underestimated.")

        # Maximum barycentered arrival time error:
        # |distance from array barycenter to furthest detector| / c + 5 ms.
        # For LHO+LLO, this is 15.0 ms.
        # For an arbitrary terrestrial detector network, the maximum is 26.3 ms.
        max_abs_t = np.max(
            np.sqrt(np.sum(np.square(locations / lal.C_SI), axis=1))) + 0.005

        acors, sample_rates = zip(
            *[filter.autocorrelation(_, max_abs_t) for _ in HS])
        sample_rate = sample_rates[0]
        deltaT = 1 / sample_rate
        nsamples = len(acors[0])
        assert all(sample_rate == _ for _ in sample_rates)
        assert all(nsamples == len(_) for _ in acors)
        nsamples = nsamples * 2 - 1

        snr_series = []
        for acor, sngl in zip(acors, sngl_inspirals):
            series = lal.CreateCOMPLEX8TimeSeries(
                'fake SNR', 0, 0, deltaT, lal.StrainUnit, nsamples)
            series.epoch = sngl.end - 0.5 * (nsamples - 1) * deltaT
            acor = np.concatenate((np.conj(acor[:0:-1]), acor))
            if phase_convention.lower() == 'antifindchirp':
                # The matched filter phase convention does NOT affect the
                # template autocorrelation sequence; however it DOES affect
                # the maximum-likelihood phase estimate AND the SNR time series.
                # So if we are going to apply the anti-findchirp phase
                # correction later, we'll have to apply a complex conjugate to
                # the autocorrelation sequence to cancel it here.
                acor = np.conj(acor)
            series.data.data = sngl.snr * filter.exp_i(sngl.coa_phase) * acor
            snr_series.append(series)

    # Ensure that all of the SNR time series have the same sample rate.
    # FIXME: for now, the Python wrapper expects all of the SNR time sries to
    # also be the same length.
    deltaT = snr_series[0].deltaT
    sample_rate = 1 / deltaT
    if any(deltaT != series.deltaT for series in snr_series):
        raise ValueError('BAYESTAR does not yet support SNR time series with mixed sample rates')

    # FIXME: for now, the Python wrapper expects all of the SNR time sries to
    # also be the same length.
    nsamples = len(snr_series[0].data.data)
    if any(nsamples != len(series.data.data) for series in snr_series):
        raise ValueError('BAYESTAR does not yet support SNR time series of mixed lengths')

    # Perform sanity checks that the middle sample of the SNR time series match
    # the sngl_inspiral records.
    for sngl_inspiral, series in zip(sngl_inspirals, snr_series):
        if np.abs(0.5 * (nsamples - 1) * series.deltaT
                  + float(series.epoch - sngl_inspiral.end)) >= 0.5 * deltaT:
            raise ValueError('BAYESTAR expects the SNR time series to be centered on the sngl_inspiral end times')

    # Extract the TOAs in GPS nanoseconds from the SNR time series, assuming
    # that the trigger happened in the middle.
    toas_ns = [series.epoch.ns() + 1e9 * 0.5 * (len(series.data.data) - 1)
               * series.deltaT for series in snr_series]

    # Collect all of the SNR series in one array.
    snr_series = np.vstack([series.data.data for series in snr_series])

    # Fudge factor for excess estimation error in gstlal_inspiral.
    fudge = 0.83
    snr_series *= fudge

    # If using 'findchirp' phase convention rather than gstlal/mbta,
    # then flip signs of phases.
    if phase_convention.lower() == 'antifindchirp':
        log.warn('Using anti-FINDCHIRP phase convention; inverting phases. '
                 'This is currently the default and it is appropriate for '
                 'gstlal and MBTA but not pycbc as of observing run 1 ("O1"). '
                 'The default setting is likely to change in the future.')
        snr_series = np.conj(snr_series)

    # Center times of arrival and compute GMST at mean arrival time.
    # Pre-center in integer nanoseconds to preserve precision of
    # initial datatype.
    epoch = sum(toas_ns) // len(toas_ns)
    toas = 1e-9 * (np.asarray(toas_ns) - epoch)
    mean_toa = np.average(toas, weights=weights, axis=0)
    toas -= mean_toa
    epoch += int(np.round(1e9 * mean_toa))
    epoch = lal.LIGOTimeGPS(0, int(epoch))
    gmst = lal.GreenwichMeanSiderealTime(epoch)

    # Translate SNR time series back to time of first sample.
    toas -= 0.5 * (nsamples - 1) * deltaT

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

    # Rescale distances to horizon distance of most sensitive detector.
    max_horizon = np.max(horizons)
    horizons /= max_horizon
    min_distance /= max_horizon
    max_distance /= max_horizon

    # Use KDE for density estimation?
    if method.endswith('_kde'):
        method = method[:-4]
        kde = True
    else:
        kde = False

    # Time and run sky localization.
    log.debug('starting computationally-intensive section')
    start_time = lal.GPSTimeNow()
    if method == "toa_phoa_snr":
        skymap = Table(_sky_map.toa_phoa_snr(
            min_distance, max_distance, prior_distance_power, gmst, sample_rate,
            toas, snr_series, responses, locations, horizons))
    elif method == "toa_phoa_snr_mcmc":
        # prob = emcee_sky_map(
        #     logl=_sky_map.log_likelihood_toa_phoa_snr,
        #     loglargs=(gmst, sample_rate, toas, snr_series, responses, locations,
        #         horizons),
        #     logp=toa_phoa_snr_log_prior,
        #     logpargs=(min_distance, max_distance, prior_distance_power,
        #         max_abs_t),
        #     xmin=[0, -1, min_distance, -1, 0, 0],
        #     xmax=[2*np.pi, 1, max_distance, 1, 2*np.pi, 2 * max_abs_t],
        #     nside=nside, kde=kde, chain_dump=chain_dump,
        #     max_horizon=max_horizon)
        raise NotImplemented('Working on porting this to new MOC format')
    else:
        raise ValueError("Unrecognized method: %s" % method)

    # Convert distance moments to parameters
    distmean = skymap.columns.pop('DISTMEAN')
    diststd = skymap.columns.pop('DISTSTD')
    skymap['DISTMU'], skymap['DISTSIGMA'], skymap['DISTNORM'] = \
        distance.moments_to_parameters(distmean, diststd)

    # Add marginal distance moments
    good = np.isfinite(distmean) & np.isfinite(diststd)
    prob = (moc.uniq2pixarea(skymap['UNIQ']) * skymap['PROBDENSITY'])[good]
    distmean = distmean[good]
    diststd = diststd[good]
    rbar = (prob * distmean).sum()
    r2bar = (prob * (np.square(diststd) + np.square(distmean))).sum()
    skymap.meta['distmean'] = rbar
    skymap.meta['diststd'] = np.sqrt(r2bar - np.square(rbar))

    # Rescale
    rescale = max_horizon * fudge
    skymap['DISTMU'] *= rescale
    skymap['DISTSIGMA'] *= rescale
    skymap.meta['distmean'] *= rescale
    skymap.meta['diststd'] *= rescale
    skymap['DISTNORM'] /= np.square(rescale)

    end_time = lal.GPSTimeNow()
    log.debug('finished computationally-intensive section')

    # Fill in metadata and return.
    skymap.meta['creator'] = 'BAYESTAR'
    skymap.meta['origin'] = 'LIGO/Virgo'
    skymap.meta['gps_time'] = float(epoch)
    skymap.meta['runtime'] = float(end_time - start_time)
    skymap.meta['instruments'] = {sngl_inspiral.ifo for sngl_inspiral in sngl_inspirals}
    skymap.meta['gps_creation_time'] = end_time

    return skymap


def gracedb_sky_map(
        coinc_file, psd_file, waveform, f_low, min_distance=None,
        max_distance=None, prior_distance_power=None,
        method="toa_phoa_snr", nside=-1, chain_dump=None,
        f_high_truncate=1.0, enable_snr_series=False):
    # Read input file.
    log.debug('reading coinc file')
    xmldoc, _ = ligolw_utils.load_fileobj(
        coinc_file, contenthandler=ligolw.LSCTablesAndSeriesContentHandler)

    # Locate the tables that we need.
    coinc_inspiral_table = ligolw_table.get_table(xmldoc,
        lsctables.CoincInspiralTable.tableName)
    coinc_map_table = ligolw_table.get_table(xmldoc,
        lsctables.CoincMapTable.tableName)
    sngl_inspiral_table = ligolw_table.get_table(xmldoc,
        lsctables.SnglInspiralTable.tableName)

    # Attempt to determine phase convention from process table.
    try:
        process_table = ligolw_table.get_table(xmldoc,
            lsctables.ProcessTable.tableName)
        process, = process_table
        process_name = process.program.lower()
    except ValueError:
        process_name = ''
    if 'pycbc' in process_name:
        phase_convention = 'findchirp'
    else:
        phase_convention = 'antifindchirp'

    # Locate the sngl_inspiral rows that we need.
    coinc_inspiral = coinc_inspiral_table[0]
    coinc_event_id = coinc_inspiral.coinc_event_id
    event_ids = [coinc_map.event_id for coinc_map in coinc_map_table
        if coinc_map.coinc_event_id == coinc_event_id]
    sngl_inspirals = [next((sngl_inspiral for sngl_inspiral in sngl_inspiral_table
        if sngl_inspiral.event_id == event_id)) for event_id in event_ids]

    # Try to load complex SNR time series.
    snrs = ligolw.snr_series_by_sngl_inspiral_id_for_xmldoc(xmldoc)
    try:
        snrs = [snrs[sngl.event_id] for sngl in sngl_inspirals]
    except KeyError:
        snrs = None

    # Read PSDs.
    log.debug('reading PSDs time series')
    xmldoc, _ = ligolw_utils.load_fileobj(
        psd_file, contenthandler=lal.series.PSDContentHandler)
    psds = lal.series.read_psd_xmldoc(xmldoc, root_name=None)

    # Rearrange PSDs into the same order as the sngl_inspirals.
    psds = [psds[sngl_inspiral.ifo] for sngl_inspiral in sngl_inspirals]

    # Interpolate PSDs.
    log.debug('constructing PSD interpolants')
    psds = [timing.InterpolatedPSD(filter.abscissa(psd), psd.data.data,
            f_high_truncate=f_high_truncate)
        for psd in psds]

    # Run sky localization
    return ligolw_sky_map(sngl_inspirals, waveform, f_low,
        min_distance, max_distance, prior_distance_power, method=method,
        nside=nside, psds=psds, phase_convention=phase_convention,
        chain_dump=chain_dump, snr_series=snrs,
        enable_snr_series=enable_snr_series)


def rasterize(skymap):
    skymap = Table(moc.rasterize(skymap), meta=skymap.meta)
    skymap.rename_column('PROBDENSITY', 'PROB')
    skymap['PROB'] *= 4 * np.pi / len(skymap)
    return skymap


def test():
    """Run BAYESTAR C unit tests.
    >>> test()
    0
    """
    return int(_sky_map.test())
