# Copyright (C) 2010  Nickolas Fotopoulos
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
Utilities to perform operations on (polar, azimuthal) vectors.
"""

from __future__ import division

import numpy as np
np.seterr(all="raise")

import lal

__author__ = "Nickolas Fotopoulos <nvf@gravity.phys.uwm.edu>"


#
# Rotation utilities
#

def rotate_euler(sph_coords, alpha, beta, gamma):
    """
    Take an Nx2 array of N (theta, phi) vectors on the unit sphere
    (that is, (polar, azimuthal) angles in radians) and apply
    rotations through the Euler angles alpha, beta, and gamma
    in radians, using the ZXZ convention.
    """
    c = np.cos
    s = np.sin

    # Define rotation matrix
    R = np.array(
        [[c(alpha) * c(gamma) - c(beta) * s(alpha) * s(gamma),
          c(gamma) * s(alpha) + c(alpha) * c(beta) * s(gamma),
          s(beta) * s(gamma)],
         [-c(beta) * c(gamma) * s(alpha) - c(alpha) * s(gamma),
          c(alpha) * c(beta) * c(gamma) - s(alpha) * s(gamma),
          c(gamma) * s(beta)],
         [s(alpha) * s(beta), -c(alpha) * s(beta), c(beta)]], dtype=float)

    # Convert to intermediate cartesian representation (N x 3)
    cart_orig = np.empty(shape=(len(sph_coords), 3), dtype=float)
    cart_orig[:, 0] = c(sph_coords[:, 1]) # * s(sph_coords[:, 0])
    cart_orig[:, 1] = s(sph_coords[:, 1]) # * s(sph_coords[:, 0])
    cart_orig[:, 0:2] *= s(sph_coords[:, 0])[:, None]
    cart_orig[:, 2] = c(sph_coords[:, 0])

    # Rotate by x_i = R_{ij} * x_j
    # NB: dot contracts last dim of A with second-to-last dim of B
    cart_new = np.dot(cart_orig, R)

    # Extract new spherical coordinates
    sph_new = np.empty(shape=(len(sph_coords), 2), dtype=float)
    sph_new[:, 0] = np.arccos(cart_new[:, 2])
    sph_new[:, 1] = np.arctan2(cart_new[:, 1], cart_new[:, 0])

    return sph_new

def new_z_to_euler(new_z):
    """
    From the new Z axis expressed in (polar, azimuthal) angles of the
    initial coordinate system, return the (alpha, beta) Euler angles
    that rotate the old Z axis (0, 0) to the new Z axis.
    """
    return (lal.PI_2 + new_z[1]) % (2 * lal.PI), new_z[0]

def rotate_about_axis(x, axis, ang):
    """
    Return the result of rotating x about axis by ang (in radians).
    axis must be a unit vector.
    """
    if abs(np.dot(axis, axis) - 1) >= 1e-6:
        raise ValueError, "axis must be a unit vector"
    if len(x) != 3:
        raise ValueError, "x must be three-dimensional"
    if len(axis) != 3:
        raise ValueError, "axis must be three-dimensional"
    cosa = np.cos(ang)
    sina = np.sin(ang)

    # Rodrigues' rotation formula
    R = np.array(
        [[cosa+axis[0]*axis[0]*(1.0-cosa), 
          axis[0]*axis[1]*(1.0-cosa)-axis[2]*sina,
          axis[0]*axis[2]*(1.0-cosa)+axis[1]*sina],
         [axis[1]*axis[0]*(1.0-cosa)+axis[2]*sina,
          cosa+axis[1]*axis[1]*(1.0-cosa),
          axis[1]*axis[2]*(1.0-cosa)-axis[0]*sina],
         [axis[2]*axis[0]*(1.0-cosa)-axis[1]*sina,
          axis[2]*axis[1]*(1.0-cosa)+axis[0]*sina,
          cosa+axis[2]*axis[2]*(1.0-cosa)]])
    return np.dot(R, x)
#
# Utilities to find the angle between two points
#

def _abs_diff(c):
    """
    For some angular difference c = |a - b| in radians, find the
    magnitude of the difference, taking into account the wrap-around at 2*pi.
    """
    c = abs(c) % (2 * lal.PI)
    # XXX: numpy 1.3.0 introduces fmin, which is more elegant
    return np.where(c < lal.PI, c, 2 * lal.PI - c)

def _haversine(angle):
    return np.sin(angle / 2)**2

def _archaversine(h):
    """
    Compute the inverse of the _haversine function, using clip as protection
    for the antipodes.
    """
    h = np.clip(h, 0., 1.)
    return 2 * np.arcsin(np.sqrt(h))

def angle_between_points(a, b):
    """
    Find the angle in radians between a and b, each expressed in
    (polar, azimuthal) angles in radians. If a and b are Nx2 arrays
    of vectors, then return the N pairwise angles between them.

    This formula is the law of _haversines, which is derivable from the
    spherical law of cosines, but is more numerically stable about a == b.
    This technique is slightly unstable for antipodal a and b.
    """
    s = np.sin
    dtheta, dphi = (a - b).T
    h = _haversine(dtheta) + s(a[..., 0]) * s(b[..., 0]) * _haversine(dphi)
    result = _abs_diff(_archaversine(h))
    if result.shape == (1,):
        return result[0]
    else:
        return result

#
# Implement the Fisher distribution
#

def fisher_rvs(mu, kappa, size=1):
    """
    Return a random (polar, azimuthal) angle drawn from the Fisher
    distribution. Assume that the concentration parameter (kappa) is large
    so that we can use a Rayleigh distribution about the north pole and
    rotate it to be centered at the (polar, azimuthal) coordinate mu.

    Assume kappa = 1 / sigma**2

    References:
      * http://en.wikipedia.org/wiki/Von_Mises-Fisher_distribution
      * http://arxiv.org/pdf/0902.0737v1 (states the Rayleigh limit)
    """
    rayleigh_rv = \
        np.array((np.random.rayleigh(scale=1. / np.sqrt(kappa), size=size),
                  np.random.uniform(low=0, high=2*lal.PI, size=size)))\
                .reshape((2, size)).T  # guarantee 2D and transpose
    a, b = new_z_to_euler(mu)
    return rotate_euler(rayleigh_rv, a, b, 0)

def fisher_pdf(theta, kappa):
    """
    Return the PDF of theta, the opening angle of X with respect to mu where
    X is Fisher-distributed about mu. See fisher_rvs for the definition of mu.
    """
    return kappa / (2 * np.sinh(kappa)) * np.exp(kappa * np.cos(theta))\
        * np.sin(theta)

def fisher_cdf(theta, kappa):
    return 0.5 * (np.exp(kappa) - np.exp(kappa * np.cos(theta))) / np.sinh(kappa)
