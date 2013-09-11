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


def flood_fill(nside, ipix, m):
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
            neighbors = hp.get_all_neighbours(nside, ipix)
            # All pixels have up to 8 neighbors. If a pixel has less than 8
            # neighbors, then some entries of the array are set to -1. We
            # have to skip those.
            neighbors = neighbors[neighbors != -1]
            # Push neighboring pixels onto the stack.
            stack.extend(neighbors)


def count_modes(m):
    """Count the number of modes in a binary HEALPix image by repeatedly
    applying the flood-fill algorithm.

    WARNING: The input array is clobbered in the process."""
    npix = len(m)
    nside = hp.npix2nside(npix)
    for nmodes in xrange(npix):
        nonzeroipix = np.flatnonzero(m)
        if len(nonzeroipix):
            flood_fill(nside, nonzeroipix[0], m)
        else:
            break
    return nmodes


def indicator(n, i):
    """Create a binary array of length n that is True for every index that is in
    i and False for every other index. Named after the indicator function."""
    m = np.zeros(n, dtype=np.bool)
    np.put(m, i, True)
    return m


def angle_distance(theta0, phi0, theta1, phi1):
    """Angular separation in radians between two points on the unit sphere."""
    cos_angle_distance = (np.cos(phi1 - phi0) * np.sin(theta0) * np.sin(theta1)
        + np.cos(theta0) * np.cos(theta1))
    if cos_angle_distance > 1:
        return 0.
    elif cos_angle_distance < -1:
        return np.pi
    else:
        return np.arccos(cos_angle_distance)


# Class to hold return value of find_injection method
FoundInjection = collections.namedtuple('FoundInjection',
    'searched_area searched_prob offset searched_modes contour_areas contour_modes')


def find_injection(sky_map, true_ra, true_dec, contours=()):
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
    true_pix = hp.ang2pix(nside, true_theta, true_phi)

    # Compute spherical polar coordinates of true location.
    mode_theta, mode_phi = hp.pix2ang(nside, mode_pix)

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

    # Count up the number of modes in each of the given contours.
    searched_modes = count_modes(indicator(npix, indices[:idx+1]))
    contour_modes = [count_modes(indicator(npix, indices[:i+1]))
        for i in ipix]

    # Done.
    return FoundInjection(searched_area, searched_prob, offset, searched_modes, contour_areas, contour_modes)
