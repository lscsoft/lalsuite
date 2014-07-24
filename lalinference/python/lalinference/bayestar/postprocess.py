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
Postprocessing utilities for HEALPix sky maps
"""
from __future__ import division
__author__ = "Leo Singer <leo.singer@ligo.org>"


import numpy as np
import healpy as hp
import collections
import itertools


_max_order = 29


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


class _HEALPixNode(object):
    """Data structure used internally by the function
    adaptive_healpix_histogram()."""

    def __init__(self, samples, max_samples_per_pixel, max_order, order=0, needs_sort=True):
        if needs_sort:
            samples = np.sort(samples)
        if len(samples) > max_samples_per_pixel and order < max_order:
            # All nodes have 4 children, except for the root node, which has 12.
            nchildren = 12 if order == 0 else 4
            self.samples = None
            self.children = [_HEALPixNode([], max_samples_per_pixel,
                max_order, order=order + 1) for i in range(nchildren)]
            for ipix, samples in itertools.groupby(samples,
                    self.key_for_order(order)):
                self.children[np.uint64(ipix % nchildren)] = _HEALPixNode(
                    list(samples), max_samples_per_pixel, max_order,
                    order=order + 1, needs_sort=False)
        else:
            # There are few enough samples that we can make this cell a leaf.
            self.samples = list(samples)
            self.children = None

    @staticmethod
    def key_for_order(order):
        """Create a function that downsamples full-resolution pixel indices."""
        return lambda ipix: ipix >> np.uint64(2 * (_max_order - order))

    @property
    def order(self):
        """Return the maximum HEALPix order required to represent this tree,
        which is the same as the tree depth."""
        if self.children is None:
            return 0
        else:
            return 1 + max(child.order for child in self.children)

    def _flat_bitmap(self, order, full_order, ipix, m):
        if self.children is None:
            nside = 1 << order
            ipix0 = ipix << 2 * (full_order - order)
            ipix1 = (ipix + 1) << 2 * (full_order - order)
            m[ipix0:ipix1] = len(self.samples) / hp.nside2pixarea(nside)
        else:
            for i, child in enumerate(self.children):
                child._flat_bitmap(order + 1, full_order, (ipix << 2) + i, m)

    @property
    def flat_bitmap(self):
        """Return flattened HEALPix representation."""
        order = self.order
        nside = 1 << order
        npix = hp.nside2npix(nside)
        m = np.empty(npix)
        for ipix, child in enumerate(self.children):
            child._flat_bitmap(0, order, ipix, m)
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
    ipix = hp.ang2pix(1 << _max_order, theta, phi, nest=True).astype(np.uint64)

    # Build tree structure.
    if nside == -1 and max_nside == -1:
        max_order = _max_order
    elif nside == -1:
        max_order = nside2order(max_nside)
    elif max_nside == -1:
        max_order = nside2order(nside)
    else:
        max_order = nside2order(min(nside, max_nside))
    tree = _HEALPixNode(ipix, max_samples_per_pixel, max_order)

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
    for nmodes in xrange(npix):
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
    'searched_area searched_prob offset searched_modes contour_areas contour_modes')


def find_injection(sky_map, true_ra, true_dec, contours=(), modes=False, nest=False):
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
    idx = (i for i, pix in enumerate(indices) if pix == true_pix).next()

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
    return FoundInjection(searched_area, searched_prob, offset, searched_modes, contour_areas, contour_modes)
