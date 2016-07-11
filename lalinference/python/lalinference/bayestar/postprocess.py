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


# Maximum 64-bit HEALPix resolution.
HEALPIX_MACHINE_ORDER = 29
HEALPIX_MACHINE_NSIDE = hp.order2nside(HEALPIX_MACHINE_ORDER)


def nside2order(nside):
    """Convert lateral HEALPix resolution to order.
    FIXME: see https://github.com/healpy/healpy/issues/163"""
    order = np.log2(nside)
    int_order = int(order)
    if order != int_order:
        raise ValueError('not a valid value for nside: {0}'.format(nside))
    return int_order


def order2nside(order):
    return 1 << order


class HEALPixTree(object):
    """Data structure used internally by the function
    adaptive_healpix_histogram()."""

    def __init__(self, samples, max_samples_per_pixel, max_order, order=0, needs_sort=True):
        if needs_sort:
            samples = np.sort(samples)
        if len(samples) >= max_samples_per_pixel and order < max_order:
            # All nodes have 4 children, except for the root node, which has 12.
            nchildren = 12 if order == 0 else 4
            self.samples = None
            self.children = [HEALPixTree([], max_samples_per_pixel,
                max_order, order=order + 1) for i in range(nchildren)]
            for ipix, samples in itertools.groupby(samples,
                    self.key_for_order(order)):
                self.children[np.uint64(ipix % nchildren)] = HEALPixTree(
                    list(samples), max_samples_per_pixel, max_order,
                    order=order + 1, needs_sort=False)
        else:
            # There are few enough samples that we can make this cell a leaf.
            self.samples = list(samples)
            self.children = None

    @staticmethod
    def key_for_order(order):
        """Create a function that downsamples full-resolution pixel indices."""
        return lambda ipix: ipix >> np.uint64(2 * (HEALPIX_MACHINE_ORDER - order))

    @property
    def order(self):
        """Return the maximum HEALPix order required to represent this tree,
        which is the same as the tree depth."""
        if self.children is None:
            return 0
        else:
            return 1 + max(child.order for child in self.children)

    def _visit(self, order, full_order, ipix):
        if self.children is None:
            nside = 1 << order
            full_nside = 1 << order
            ipix0 = ipix << 2 * (full_order - order)
            ipix1 = (ipix + 1) << 2 * (full_order - order)
            yield nside, full_nside, ipix, ipix0, ipix1, self.samples
        else:
            for i, child in enumerate(self.children):
                # FIXME: Replace with `yield from` in Python 3
                for _ in child._visit(order + 1, full_order, (ipix << 2) + i):
                    yield _

    def visit(self):
        """Evaluate a function on each leaf node of the HEALPix tree.

        Yields
        -------------------
        nside : int
            The HEALPix resolution of the node.

        full_nside : int
            The HEALPix resolution of the deepest node in the tree.

        ipix : int
            The nested HEALPix index of the node.

        ipix0 : int
            The start index of the range of pixels spanned by the node at the
            resolution `full_nside`.

        ipix1 : int
            The end index of the range of pixels spanned by the node at the
            resolution `full_nside`.

        samples : list
            The list of samples contained in the node.
        """
        order = self.order
        for ipix, child in enumerate(self.children):
            # FIXME: Replace with `yield from` in Python 3
            for _ in child._visit(0, order, ipix):
                yield _

    @property
    def flat_bitmap(self):
        """Return flattened HEALPix representation."""
        m = np.empty(hp.nside2npix(hp.order2nside(self.order)))
        for nside, full_nside, ipix, ipix0, ipix1, samples in self.visit():
            m[ipix0:ipix1] = len(samples) / hp.nside2pixarea(nside)
        return m


def adaptive_healpix_histogram(theta, phi, max_samples_per_pixel, nside=-1, max_nside=-1, nest=False):
    """Adaptively histogram the posterior samples represented by the
    (theta, phi) points using a recursively subdivided HEALPix tree. Nodes are
    subdivided until each leaf contains no more than max_samples_per_pixel
    samples. Finally, the tree is flattened to a fixed-resolution HEALPix image
    with a resolution appropriate for the depth of the tree. If nside is
    specified, the result is resampled to another desired HEALPix resolution."""
    # Calculate pixel index of every sample, at the maximum 64-bit resolution.
    #
    # At this resolution, each pixel is only 0.2 mas across; we'll use the
    # 64-bit pixel indices as a proxy for the true sample coordinates so that
    # we don't have to do any trigonometry (aside from the initial hp.ang2pix
    # call).
    #
    # FIXME: Cast to uint64 needed because Healpy returns signed indices.
    ipix = hp.ang2pix(HEALPIX_MACHINE_NSIDE, theta, phi, nest=True).astype(np.uint64)

    # Build tree structure.
    if nside == -1 and max_nside == -1:
        max_order = HEALPIX_MACHINE_ORDER
    elif nside == -1:
        max_order = nside2order(max_nside)
    elif max_nside == -1:
        max_order = nside2order(nside)
    else:
        max_order = nside2order(min(nside, max_nside))
    tree = HEALPixTree(ipix, max_samples_per_pixel, max_order)

    # Compute a flattened bitmap representation of the tree.
    p = tree.flat_bitmap

    # If requested, resample the tree to the output resolution.
    if nside != -1:
        p = hp.ud_grade(p, nside, order_in='NESTED', order_out='NESTED')

    # Normalize.
    p /= np.sum(p)

    if not nest:
        p = hp.reorder(p, n2r=True)

    # Done!
    return p


def _interpolate_level(m):
    """Recursive multi-resolution interpolation. Modifies `m` in place."""
    # Determine resolution.
    npix = len(m)

    if npix > 12:
        # Determine which pixels comprise multi-pixel tiles.
        ipix = np.flatnonzero(
            (m[0::4] == m[1::4]) &
            (m[0::4] == m[2::4]) &
            (m[0::4] == m[3::4]))

        if len(ipix):
            ipix = (4 * ipix +
                np.expand_dims(np.arange(4, dtype=np.intp), 1)).T.ravel()

            nside = hp.npix2nside(npix)

            # Downsample.
            m_lores = hp.ud_grade(
                m, nside // 2, order_in='NESTED', order_out='NESTED')

            # Interpolate recursively.
            _interpolate_level(m_lores)

            # Record interpolated multi-pixel tiles.
            m[ipix] = hp.get_interp_val(
                m_lores, *hp.pix2ang(nside, ipix, nest=True), nest=True)


def interpolate_nested(m, nest=False):
    """
    Apply bilinear interpolation to a multiresolution HEALPix map, assuming
    that runs of pixels containing identical values are nodes of the tree. This
    smooths out the stair-step effect that may be noticeable in contour plots.

    Here is how it works. Consider a coarse tile surrounded by base tiles, like
    this:

                +---+---+
                |   |   |
                +-------+
                |   |   |
        +---+---+---+---+---+---+
        |   |   |       |   |   |
        +-------+       +-------+
        |   |   |       |   |   |
        +---+---+---+---+---+---+
                |   |   |
                +-------+
                |   |   |
                +---+---+

    The value within the central coarse tile is computed by downsampling the
    sky map (averaging the fine tiles), upsampling again (with bilinear
    interpolation), and then finally copying the interpolated values within the
    coarse tile back to the full-resolution sky map. This process is applied
    recursively at all successive HEALPix resolutions.

    Note that this method suffers from a minor discontinuity artifact at the
    edges of regions of coarse tiles, because it temporarily treats the
    bordering fine tiles as constant. However, this artifact seems to have only
    a minor effect on generating contour plots.

    Parameters
    ----------

    m: `~numpy.ndarray`
        a HEALPix array

    nest: bool, default: False
        Whether the input array is stored in the `NESTED` indexing scheme (True)
        or the `RING` indexing scheme (False).

    """
    # Convert to nest indexing if necessary, and make sure that we are working
    # on a copy.
    if nest:
        m = m.copy()
    else:
        m = hp.reorder(m, r2n=True)

    _interpolate_level(m)

    # Convert to back ring indexing if necessary
    if not nest:
        m = hp.reorder(m, n2r=True)

    # Done!
    return m


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


def find_greedy_credible_levels(p):
    """Find the greedy credible levels of a (possibly multi-dimensional) array.

    Parameters
    ----------

    p : np.ndarray
        The input array, typically a HEALPix image.

    Returns
    -------

    cls : np.ndarray
        An array with the same shape as `p`, with values ranging from `0`
        to `p.sum()`, representing the greedy credible level to which each
        entry in the array belongs.
    """
    pflat = p.ravel()
    i = np.flipud(np.argsort(pflat))
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
