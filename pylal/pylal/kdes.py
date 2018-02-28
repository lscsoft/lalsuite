# -*- coding: utf-8 -*-
#
#       kdes.py: KDE utilities for density estimation in unusual
#       topologies.
#
#       Copyright 2012
#       Will M. Farr <will.farr@ligo.org>
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.


import numpy as np
from scipy.stats import gaussian_kde

class BoundedKDE(gaussian_kde):
    """Density estimation using a KDE on bounded domains.

    Bounds can be any combination of low or high (if no bound, set to
    ``float('inf')`` or ``float('-inf')``), and can be periodic or
    non-periodic.  Cannot handle topologies that have
    multi-dimensional periodicities; will only handle topologies that
    are direct products of (arbitrary numbers of) R, [0,1], and S1.

    :param pts:
        ``(Ndim, Npts)`` shaped array of points (as in :class:`gaussian_kde`).

    :param low: 
        Lower bounds; if ``None``, assume no lower bounds.

    :param high:
        Upper bounds; if ``None``, assume no upper bounds.

    :param periodic:
        Boolean array giving periodicity in each dimension; if
        ``None`` assume no dimension is periodic.

    :param bw_method: (optional)
        Bandwidth estimation method (see :class:`gaussian_kde`)."""
    def __init__(self, pts, low=None, high=None, periodic=None, bw_method=None):

        if bw_method is None:
            bw_method='scott'

        if low is None:
            self._low = np.zeros(pts.shape[0])
            self._low[:] = float('-inf')
        else:
            self._low = low

        if high is None:
            self._high = np.zeros(pts.shape[0])
            self._high[:] = float('inf')
        else:
            self._high = high

        if periodic is None:
            self._periodic = np.zeros(pts.shape[0], dtype=np.bool)
        else:
            self._periodic = periodic

        super(BoundedKDE, self).__init__(pts, bw_method=bw_method)

    def evaluate(self, pts):
        """Evaluate the KDE at the given points."""

        pts_orig=pts
        pts=np.copy(pts_orig)

        den=super(BoundedKDE, self).evaluate(pts)

        for i, (low,high,period) in enumerate(zip(self._low, self._high, self._periodic)):
            if period:
                P = high-low
                
                pts[i, :] += P
                den += super(BoundedKDE, self).evaluate(pts)

                pts[i,:] -= 2.0*P
                den += super(BoundedKDE, self).evaluate(pts)

                pts[i,:]=pts_orig[i,:]

            else:
                if not (low == float('-inf')):
                    pts[i,:] = 2.0*low - pts[i,:]
                    den += super(BoundedKDE, self).evaluate(pts)
                    pts[i,:] = pts_orig[i,:]

                if not (high == float('inf')):
                    pts[i,:] = 2.0*high - pts[i,:]
                    den += super(BoundedKDE, self).evaluate(pts)
                    pts[i,:] = pts_orig[i,:]

        return den

    __call__ = evaluate

    def quantile(self, pt):
        """Quantile of ``pt``, evaluated by a greedy algorithm.

        :param pt:
            The point at which the quantile value is to be computed.

        The quantile of ``pt`` is the fraction of points used to
        construct the KDE that have a lower KDE density than ``pt``."""

        return float(np.count_nonzero(self.evaluate(self.dataset) < self.evaluate(pt)))/float(self.n)
