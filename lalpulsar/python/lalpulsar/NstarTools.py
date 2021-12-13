# Copyright (C) 2017 Gregory Ashton
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
"""

Tools to estimate Nstar and related functions

"""
from __future__ import division, absolute_import, print_function

import pprint
import itertools
import logging
import numpy as np

try:
    import lal
except ImportError:
    raise ImportError("SWIG wrappings of LAL cannot be imported")

try:
    import lalpulsar
except ImportError:
    raise ImportError("SWIG wrappings of LALPulsar cannot be imported")


def _extract_data_from_prior(prior):
    """ Calculate the input data from the prior

    Parameters
    ----------
    prior: dict
        For each key in 'Alpha', 'Delta', 'F0', 'F1', 'F2', either a scalar
        value (not searched over), or a pair of the upper and lower limit.

    Returns
    -------
    p : ndarray
        Matrix with columns being the edges of the uniform bounding box
    spindowns : int
        The number of spindowns
    sky : bool
        If true, search includes the sky position
    fiducial_freq : float
        Fidicual frequency

    """
    keys = ['Alpha', 'Delta', 'F0', 'F1', 'F2']
    spindown_keys = keys[3:]
    sky_keys = keys[:2]
    lims = []
    lims_keys = []
    scalar_keys = []
    unfound_keys = []
    lims_idxs = []
    for i, key in enumerate(keys):
        array_like = isinstance(prior[key], (list, tuple, np.ndarray))
        if key not in keys:
            raise ValueError(
                'Key {}={} not understood'.format(key, prior[key]))
        elif array_like and len(prior[key]) == 2:
            lims.append([prior[key][0], prior[key][1]])
            lims_keys.append(key)
            lims_idxs.append(i)
        elif np.isscalar(prior[key]) or len(prior[key] == 1):
            scalar_keys.append(key)
            if key in sky_keys:
                lims.append([prior[key], 0])
        else:
            unfound_keys.append(key)
    logging.info(
        'Computing Nstar over dimensions {} while {} are scalar and {} unused'
        .format(lims_keys, scalar_keys, unfound_keys))
    lims = np.array(lims)
    lims_keys = np.array(lims_keys)
    base = lims[:, 0]
    p = [base]
    for i in lims_idxs:
        basex = base.copy()
        basex[i] = lims[i, 1]
        p.append(basex)
    spindowns = np.sum([np.sum(lims_keys == k) for k in spindown_keys])
    sky = any([key in lims_keys for key in sky_keys])
    fiducial_freq = np.max(prior['F0'])

    return np.array(p).T, spindowns, sky, fiducial_freq


def _convert_array_to_gsl_matrix(array):
    gsl_matrix = lal.gsl_matrix(*array.shape)
    gsl_matrix.data = array
    return gsl_matrix


def get_Nstar_estimate(
        nsegs, tref, minStartTime, maxStartTime, prior, detector_names,
        earth_ephem_file="earth00-40-DE405.dat.gz",
        sun_ephem_file="sun00-40-DE405.dat.gz"):
    """ Returns N* estimated from the super-sky metric

    Nstar is the approximate number of unit-mismatch templates, see
    https://dcc.ligo.org/P1700455 for further details.

    Parameters
    ----------
    nsegs : int
        Number of semi-coherent segments
    tref : int
        Reference time in GPS seconds
    minStartTime, maxStartTime : int
        Minimum and maximum SFT timestamps
    prior: dict
        For each key in 'Alpha', 'Delta', 'F0', 'F1', 'F2', either a scalar
        value (not searched over), or a pair of the upper and lower limit.
    detector_names : list
        List of detectors to average over, e.g. ['H1']

    Returns
    -------
    Nstar: int
        The estimated approximate number of templates to cover the prior
        parameter space at a mismatch of unity, assuming the normalised
        thickness is unity.

    Example
    --------

    An example for a directed search where the Nstar can be estimated from the
    metric directly. This estimate is used to define the frequency and
    spin-down uncertainty. The calculated estimate using `get_Nstar_estimate`
    agrees with this input value.

    >>> from lalpulsar import NstarTools
    >>> nsegs = 1
    >>> minStartTime = 1000000000
    >>> duration = 10 * 86400
    >>> maxStartTime = minStartTime + duration
    >>> tref = minStartTime + .5*duration
    >>> detector_names = ['H1']
    >>> F0 = 30
    >>> F1 = -1e-10
    >>> F2 = 0
    >>> Alpha = 0.5
    >>> Delta = 1.5
    >>> Nstar = 1e3
    >>> DeltaF0 = np.sqrt(Nstar) * np.sqrt(3)/(np.pi*duration)
    >>> DeltaF1 = np.sqrt(Nstar) * np.sqrt(180)/(np.pi*duration**2)
    >>> prior = {'F0': [F0-DeltaF0/2., F0+DeltaF0/2],
    >>>          'F1': [F1-DeltaF1/2., F1+DeltaF1/2],
    >>>          'F2': F2,
    >>>          'Alpha': Alpha
    >>>          'Delta': Delta}
    >>> print NstarTools.get_Nstar_estimate(
    >>>     nsegs, tref, minStartTime, maxStartTime, prior, detector_names)
    1000.00000009

    Note
    ----

    To see detailed information about each set of dimensions Nstar, add the
    following before calling get_Nstar_estimate()
    >>> import logging
    >>> logging.basicConfig(level=logging.DEBUG)

    """
    detector_names = list(detector_names)

    in_phys, spindowns, sky, fiducial_freq = _extract_data_from_prior(prior)
    out_rssky = np.zeros(in_phys.shape)

    # Convert to Alpha, Delta, Fn, Fn-1, Fn-2,.. ordering
    in_phys[2:] = in_phys[2:][::-1]

    in_phys = _convert_array_to_gsl_matrix(in_phys)
    out_rssky = _convert_array_to_gsl_matrix(out_rssky)

    tboundaries = np.linspace(minStartTime, maxStartTime, nsegs+1)

    ref_time = lal.LIGOTimeGPS(tref)
    segments = lal.SegListCreate()
    for j in range(len(tboundaries)-1):
        seg = lal.SegCreate(lal.LIGOTimeGPS(tboundaries[j]),
                            lal.LIGOTimeGPS(tboundaries[j+1]),
                            j)
        lal.SegListAppend(segments, seg)
    detNames = lal.CreateStringVector(*detector_names)
    detectors = lalpulsar.MultiLALDetector()
    lalpulsar.ParseMultiLALDetector(detectors, detNames)
    detector_weights = None
    detector_motion = (lalpulsar.DETMOTION_SPIN
                       + lalpulsar.DETMOTION_ORBIT)
    ephemeris = lalpulsar.InitBarycenter(earth_ephem_file, sun_ephem_file)
    try:
        SSkyMetric = lalpulsar.ComputeSuperskyMetrics(
            lalpulsar.SUPERSKY_METRIC_TYPE, spindowns, ref_time, segments,
            fiducial_freq, detectors, detector_weights, detector_motion,
            ephemeris)
    except RuntimeError as e:
        logging.warning('Encountered run-time error {}'.format(e))
        raise RuntimeError("Calculation of the SSkyMetric failed")

    if sky:
        i = 0
    else:
        i = 2

    lalpulsar.ConvertPhysicalToSuperskyPoints(
        out_rssky, in_phys, SSkyMetric.semi_rssky_transf)

    d = out_rssky.data
    g = SSkyMetric.semi_rssky_metric.data

    g[2:, 2:] = g[2:, 2:][::-1, ::-1]  # Convert to Alpha, Delta, F0, F1, F2

    # Remove sky if required
    g = g[i:, i:]
    d = d[i:, i:]

    # Calculate parallelepiped of differences
    parallelepiped = (d[:, 1:].T - d[:, 0]).T

    labels = np.array(['Alpha', 'Delta', 'F0', 'F1', 'F2'])[i:]

    # Compute Nstar over all allowed combinations then take maximum
    Nstars = {}
    for j in range(1, len(g)+1):
        for idx in itertools.combinations(np.arange(len(g)), j):
            idx = list(idx)
            dV = np.abs(np.linalg.det(parallelepiped[idx][:, idx]))
            sqrtdetG = np.sqrt(np.abs(np.linalg.det(g[idx][:, idx])))
            Nstars[''.join(labels[idx])] = sqrtdetG * dV
    logging.debug(pprint.pformat(Nstars))
    return np.max(Nstars.values())
