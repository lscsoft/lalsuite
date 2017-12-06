#
# Copyright (C) 2013-2016  Leo Singer
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
Postprocessing utilities for HEALPix sky maps
"""
from __future__ import division
__author__ = "Leo Singer <leo.singer@ligo.org>"


import numpy as np
import healpy as hp
import collections
import itertools
import lal
import lalsimulation
from ..healpix_tree import *


def flood_fill(nside, ipix, m, nest=False):
    """Stack-based flood fill algorithm in HEALPix coordinates.
    Based on <http://en.wikipedia.org/w/index.php?title=Flood_fill&oldid=566525693#Alternative_implementations>.
    """
    # Initialize stack with starting pixel index.
    stack = [ipix]
    while stack:
        # Pop last pixel off of the stack.
        ipix = stack.pop()
        # Is this pixel in need of filling?
        if m[ipix]:
            # Fill in this pixel.
            m[ipix] = False
            # Find the pixels neighbors.
            neighbors = hp.get_all_neighbours(nside, ipix, nest=nest)
            # All pixels have up to 8 neighbors. If a pixel has less than 8
            # neighbors, then some entries of the array are set to -1. We
            # have to skip those.
            neighbors = neighbors[neighbors != -1]
            # Push neighboring pixels onto the stack.
            stack.extend(neighbors)


def count_modes(m, nest=False):
    """Count the number of modes in a binary HEALPix image by repeatedly
    applying the flood-fill algorithm.

    WARNING: The input array is clobbered in the process."""
    npix = len(m)
    nside = hp.npix2nside(npix)
    for nmodes in range(npix):
        nonzeroipix = np.flatnonzero(m)
        if len(nonzeroipix):
            flood_fill(nside, nonzeroipix[0], m, nest=nest)
        else:
            break
    return nmodes


def indicator(n, i):
    """Create a binary array of length n that is True for every index that is in
    i and False for every other index. Named after the indicator function."""
    m = np.zeros(n, dtype=np.bool)
    np.put(m, i, True)
    return m


def cos_angle_distance(theta0, phi0, theta1, phi1):
    """Cosine of angular separation in radians between two points on the unit sphere."""
    cos_angle_distance = (np.cos(phi1 - phi0) * np.sin(theta0) * np.sin(theta1)
        + np.cos(theta0) * np.cos(theta1))
    return np.clip(cos_angle_distance, -1, 1)

def angle_distance(theta0, phi0, theta1, phi1):
    """Angular separation in radians between two points on the unit sphere."""
    return np.arccos(cos_angle_distance(theta0, phi0, theta1, phi1))


# Class to hold return value of find_injection method
FoundInjection = collections.namedtuple('FoundInjection',
    'searched_area searched_prob offset searched_modes contour_areas area_probs contour_modes')


def find_injection(sky_map, true_ra, true_dec, contours=(), areas=(), modes=False, nest=False):
    """
    Given a sky map and the true right ascension and declination (in radians),
    find the smallest area in deg^2 that would have to be searched to find the
    source, the smallest posterior mass, and the angular offset in degrees from
    the true location to the maximum (mode) of the posterior. Optionally, also
    compute the areas of and numbers of modes within the smallest contours
    containing a given total probability.
    """

    # Compute the HEALPix lateral resolution parameter for this sky map.
    npix = len(sky_map)
    nside = hp.npix2nside(npix)
    deg2perpix = hp.nside2pixarea(nside, degrees=True)

    # Convert from ra, dec to conventional spherical polar coordinates.
    true_theta = 0.5 * np.pi - true_dec
    true_phi = true_ra

    # Find the HEALPix pixel index of the mode of the posterior and of the
    # true sky location.
    mode_pix = np.argmax(sky_map)
    true_pix = hp.ang2pix(nside, true_theta, true_phi, nest=nest)

    # Compute spherical polar coordinates of true location.
    mode_theta, mode_phi = hp.pix2ang(nside, mode_pix, nest=nest)

    # Sort the pixels in the sky map by descending posterior probability and
    # form the cumulative sum.  Record the total value.
    indices = np.argsort(sky_map)[::-1]
    cum_sky_map = np.cumsum(sky_map[indices])

    # Find the index of the true location in the cumulative distribution.
    idx = next((i for i, pix in enumerate(indices) if pix == true_pix))

    # Find the smallest area that would have to be searched to find
    # the true location. Note that 1 is added to the index because we want
    # the **length** of the array up to and including the idx'th element,
    # not the index itself.
    searched_area = (idx + 1) * deg2perpix

    # Find the smallest posterior mass that would have to be searched to find
    # the true location.
    searched_prob = cum_sky_map[idx]

    # Get the total number of pixels that lie inside each contour.
    ipix = np.searchsorted(cum_sky_map, contours)

    # For each of the given confidence levels, compute the area of the
    # smallest region containing that probability.
    contour_areas = (deg2perpix * (ipix + 1)).tolist()

    # For each listed area, find the probability contained within the
    # smallest credible region of that area.
    area_probs = cum_sky_map[
        np.round(np.asarray(areas) / deg2perpix).astype(np.intp)].tolist()

    # Find the angular offset between the mode and true locations.
    offset = np.rad2deg(angle_distance(true_theta, true_phi,
        mode_theta, mode_phi))

    if modes:
        # Count up the number of modes in each of the given contours.
        searched_modes = count_modes(indicator(npix, indices[:idx+1]), nest=nest)
        contour_modes = [count_modes(indicator(npix, indices[:i+1]), nest=nest)
            for i in ipix]
    else:
        searched_modes = None
        contour_modes = None

    # Done.
    return FoundInjection(
        searched_area, searched_prob, offset, searched_modes, contour_areas,
        area_probs, contour_modes)


def _norm(vertices):
    return np.sqrt(np.sum(np.square(vertices), -1))


def _adjacent_triangle_areas(vertices):
    return 0.5 * _norm(np.cross(
        np.roll(vertices, -1, axis=0) - vertices,
        np.roll(vertices, +1, axis=0) - vertices))


def _simplify(vertices, min_area):
    """Visvalingam's algorithm (see http://bost.ocks.org/mike/simplify/)
    for linear rings on a sphere. This is a naive, slow implementation."""
    area = _adjacent_triangle_areas(vertices)

    while True:
        i_min_area = np.argmin(area)
        if area[i_min_area] > min_area:
            break

        vertices = np.delete(vertices, i_min_area, axis=0)
        area = np.delete(area, i_min_area)
        new_area = _adjacent_triangle_areas(vertices)
        area = np.maximum(area, new_area)

    return vertices


def _vec2radec(vertices, degrees=False):
    theta, phi = hp.vec2ang(np.asarray(vertices))
    ret = np.column_stack((phi % (2 * np.pi), 0.5 * np.pi - theta))
    if degrees:
        ret = np.rad2deg(ret)
    return ret


def contour(m, levels, nest=False, degrees=False, simplify=True):
    import pkg_resources
    try:
        pkg_resources.require('healpy >= 1.9.0')
    except:
        raise RuntimeError('This function requires healpy >= 1.9.0.')
    try:
        import networkx as nx
    except:
        raise RuntimeError('This function requires the networkx package.')

    # Determine HEALPix resolution
    npix = len(m)
    nside = hp.npix2nside(npix)
    min_area = 0.4 * hp.nside2pixarea(nside)

    # Compute faces, vertices, and neighbors.
    # vertices is an N X 3 array of the distinct vertices of the HEALPix faces.
    # faces is an npix X 4 array mapping HEALPix faces to their vertices.
    # neighbors is an npix X 4 array mapping faces to their nearest neighbors.
    faces = np.ascontiguousarray(
            np.rollaxis(hp.boundaries(nside, np.arange(npix), nest=nest), 2, 1))
    dtype = faces.dtype
    faces = faces.view(np.dtype((np.void, dtype.itemsize * 3)))
    vertices, faces = np.unique(faces.ravel(), return_inverse=True)
    faces = faces.reshape(-1, 4)
    vertices = vertices.view(dtype).reshape(-1, 3)
    neighbors = hp.get_all_neighbours(nside, np.arange(npix), nest=nest)[::2].T

    # Loop over the requested contours.
    paths = []
    for level in levels:

        # Find credible region
        indicator = (m >= level)

        # Construct a graph of the eges of the contour.
        graph = nx.Graph()
        face_pairs = set()
        for ipix1, ipix2 in enumerate(neighbors):
            for ipix2 in ipix2:
                # Determine if we have already considered this pair of faces.
                new_face_pair = frozenset((ipix1, ipix2))
                if new_face_pair in face_pairs: continue
                face_pairs.add(new_face_pair)

                # Determine if this pair of faces are on a boundary of the
                # credible level.
                if indicator[ipix1] == indicator[ipix2]: continue

                # Add all common edges of this pair of faces.
                i1 = np.concatenate((faces[ipix1], [faces[ipix1][0]]))
                i2 = np.concatenate((faces[ipix2], [faces[ipix2][0]]))
                edges1 = frozenset(frozenset(_) for _ in zip(i1[:-1], i1[1:]))
                edges2 = frozenset(frozenset(_) for _ in zip(i2[:-1], i2[1:]))
                for edge in edges1 & edges2:
                    graph.add_edge(*edge)
        graph = nx.freeze(graph)

        # Record a closed path for each cycle in the graph.
        cycles = [
            np.take(vertices, cycle, axis=0) for cycle in nx.cycle_basis(graph)]

        # Simplify paths if requested
        if simplify:
            cycles = [_simplify(cycle, min_area) for cycle in cycles]
            cycles = [cycle for cycle in cycles if len(cycle) > 2]

        # Convert to lists
        cycles = [
            _vec2radec(cycle, degrees=degrees).tolist() for cycle in cycles]

        # Add to output paths
        paths.append([cycle + [cycle[0]] for cycle in cycles])

    return paths


def find_greedy_credible_levels(p, ranking=None):
    """Find the greedy credible levels of a (possibly multi-dimensional) array.

    Parameters
    ----------

    p : np.ndarray
        The input array, typically a HEALPix image.

    ranking : np.ndarray, optional
        The array to rank in order to determine the greedy order.
        The default is `p` itself.

    Returns
    -------

    cls : np.ndarray
        An array with the same shape as `p`, with values ranging from `0`
        to `p.sum()`, representing the greedy credible level to which each
        entry in the array belongs.
    """
    p = np.asarray(p)
    pflat = p.ravel()
    if ranking is None:
        ranking = pflat
    else:
        ranking = np.ravel(ranking)
    i = np.flipud(np.argsort(ranking))
    cs = np.cumsum(pflat[i])
    cls = np.empty_like(pflat)
    cls[i] = cs
    return cls.reshape(p.shape)


def get_detector_pair_axis(ifo1, ifo2, gmst):
    """Find the sky position where the line between two detectors pierces the
    celestial sphere.

    Parameters
    ----------

    ifo1 : str or `~lal.Detector` or `~np.ndarray`
        The first detector; either the name of the detector (e.g. `'H1'`), or a
        `lal.Detector` object (e.g., as returned by
        `lalsimulation.DetectorPrefixToLALDetector('H1')` or the geocentric
        Cartesian position of the detection in meters.

    ifo2 : str or `~lal.Detector` or `~np.ndarray`
        The second detector; same as described above.

    gmst : float
        The Greenwich mean sidereal time in radians, as returned by
        `lal.GreenwichMeanSiderealTime`.

    Returns
    -------

    pole_ra : float
        The right ascension in radians at which a ray from `ifo1` to `ifo2`
        would pierce the celestial sphere.

    pole_dec : float
        The declination in radians at which a ray from `ifo1` to `ifo2` would
        pierce the celestial sphere.

    light_travel_time : float
        The light travel time from `ifo1` to `ifo2` in seconds.
    """

    # Get location of detectors if ifo1, ifo2 are LAL detector structs
    try:
        ifo1 = lalsimulation.DetectorPrefixToLALDetector(ifo1)
    except TypeError:
        pass
    try:
        ifo1 = ifo1.location
    except AttributeError:
        pass
    try:
        ifo2 = lalsimulation.DetectorPrefixToLALDetector(ifo2)
    except TypeError:
        pass
    try:
        ifo2 = ifo2.location
    except AttributeError:
        pass

    n = ifo2 - ifo1
    light_travel_time = np.sqrt(np.sum(np.square(n))) / lal.C_SI
    (theta,), (phi,) = hp.vec2ang(n)
    pole_ra = (gmst + phi) % (2 * np.pi)
    pole_dec = 0.5 * np.pi - theta
    return pole_ra, pole_dec, light_travel_time


def rotate_map_to_axis(m, ra, dec, nest=False, method='direct'):
    """Rotate a sky map to place a given line of sight on the +z axis.

    Parameters
    ----------

    m : np.ndarray
        The input HEALPix array.

    ra : float
        Right ascension of axis in radians.

        To specify the axis in geocentric coordinates, supply ra=(lon + gmst),
        where lon is the geocentric longitude and gmst is the Greenwich mean
        sidereal time in radians.

    dec : float
        Declination of axis in radians.

        To specify the axis in geocentric coordinates, supply dec=lat,
        where lat is the geocentric latitude in radians.

    nest : bool, default=False
        Indicates whether the input sky map is in nested rather than
        ring-indexed HEALPix coordinates (default: ring).

    method : 'direct' or 'fft'
        Select whether to use spherical harmonic transformation ('fft') or
        direct coordinate transformation ('direct')

    Returns
    -------

    m_rotated : np.ndarray
        The rotated HEALPix array.
    """
    npix = len(m)
    nside = hp.npix2nside(npix)

    theta = 0.5 * np.pi - dec
    phi = ra

    if method == 'fft':
        if nest:
            m = hp.reorder(m, n2r=True)
        alm = hp.map2alm(m)
        hp.rotate_alm(alm, -phi, -theta, 0.0)
        ret = hp.alm2map(alm, nside, verbose=False)
        if nest:
            ret = hp.reorder(ret, r2n=True)
    elif method == 'direct':
        R = hp.Rotator(
            rot=np.asarray([0, theta, -phi]),
            deg=False, inv=False, eulertype='Y')
        theta, phi = hp.pix2ang(nside, np.arange(npix), nest=nest)
        ipix = hp.ang2pix(nside, *R(theta, phi), nest=nest)
        ret = m[ipix]
    else:
        raise ValueError('Unrecognized method: {0}'.format(method))

    return ret


def polar_profile(m, nest=False):
    """Obtain the marginalized polar profile of sky map.

    Parameters
    ----------

    m : np.ndarray
        The input HEALPix array.

    nest : bool, default=False
        Indicates whether the input sky map is in nested rather than
        ring-indexed HEALPix coordinates (default: ring).

    Returns
    -------

    theta : np.ndarray
        The polar angles (i.e., the colatitudes) of the isolatitude rings.

    m_int : np.ndarray
        The normalized probability density, such that `np.trapz(m_int, theta)`
        is approximately `np.sum(m)`.
    """
    npix = len(m)
    nside = hp.npix2nside(npix)
    nrings, = hp.pix2ring(nside, np.asarray([npix]))
    startpix, ringpix, costheta, sintheta, _ = hp.ringinfo(
        nside, np.arange(1, nrings))

    if nest:
        m = hp.reorder(m, n2r=True)

    theta = np.arccos(costheta)
    m_int = np.asarray([m[i:i+j].sum() * stheta * 0.5 * npix / j
        for i, j, stheta in zip(startpix, ringpix, sintheta)])

    return theta, m_int


def smooth_ud_grade(m, nside, nest=False):
    """Resample a sky map to a new resolution using bilinear interpolation.

    Parameters
    ----------

    m : np.ndarray
        The input HEALPix array.

    nest : bool, default=False
        Indicates whether the input sky map is in nested rather than
        ring-indexed HEALPix coordinates (default: ring).

    Returns
    -------

    new_m : np.ndarray
        The resampled HEALPix array. The sum of `m` is approximately preserved.
    """
    npix = hp.nside2npix(nside)
    theta, phi = hp.pix2ang(nside, np.arange(npix), nest=nest)
    new_m = hp.get_interp_val(m, theta, phi, nest=nest)
    return new_m * len(m) / len(new_m)
