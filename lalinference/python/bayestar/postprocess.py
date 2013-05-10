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


def angle_distance(theta0, phi0, theta1, phi1):
    """Angular separation in radians between two points on the unit sphere."""
    cos_angle_distance = np.cos(phi1 - phi0) * np.sin(theta0) * np.sin(theta1) + np.cos(theta0) * np.cos(theta1)
    if cos_angle_distance > 1:
        return 0.
    elif cos_angle_distance < -1:
        return np.pi
    else:
        return np.arccos(cos_angle_distance)


def find_injection(sky_map, true_ra, true_dec):
    """
    Given a sky map and the true right ascension and declination (in radians),
    find the smallest area in deg^2 that would have to be searched to find the
    source, the smallest posterior mass, and the angular offset in degrees from
    the true location to the maximum (mode) of the posterior.
    """

    # Compute the HEALPix lateral resolution parameter for this sky map.
    npix = len(sky_map)
    nside = hp.npix2nside(npix)

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

    # Find the smallest area that would have to be searched to find the true
    # location.
    searched_area = (idx + 1) * hp.nside2pixarea(nside, degrees=True)

    # Find the smallest posterior mass that would have to be searched to find
    # the true location.
    searched_prob = cum_sky_map[idx]

    # Permute the cumulative distribution so that it is indexed the same way
    # as the original sky map.
    cum_sky_map[indices] = cum_sky_map

    # Find the angular offset between the mode and true locations.
    offset = np.rad2deg(angle_distance(true_theta, true_phi, mode_theta, mode_phi))

    # Done.
    return searched_area, searched_prob, offset

