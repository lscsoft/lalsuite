# Copyright (C) 2024 Karl Wette
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

## \defgroup lalpulsar_py_metric_utils Parameter Space Metric Utilities
## \ingroup lalpulsar_python
"""
Utilities for working with parameter space metrics.
"""

import numpy as np


def mismatch_ellipse(g, mu_max, tol=0.001):
    """Return points plotting the mismatch ellipse of a metric.

    :param g: Metric.
    :param mu_max: Maximum mismatch.
    :param tol: Tolerance in generating points to plot.
    """
    assert g.shape == (2, 2)
    assert mu_max > 0
    assert tol > 0

    # starting point
    zz = [0]
    x = np.cos(zz[-1])
    y = np.sin(zz[-1])
    rr = [
        np.sqrt(
            mu_max / (g[0, 0] * x**2 + (g[0, 1] + g[1, 0]) * x * y + g[1, 1] * y**2)
        )
    ]
    dz = 2 * np.pi / 1000

    # until ellipse is complete
    while zz[-1] < 2 * np.pi:

        # next point
        z = zz[-1] + dz
        x = np.cos(z)
        y = np.sin(z)
        r = np.sqrt(
            mu_max / (g[0, 0] * x**2 + (g[0, 1] + g[1, 0]) * x * y + g[1, 1] * y**2)
        )

        # decide if next point is well-spaced from previous point
        e = abs(r - rr[-1]) / rr[-1]
        if e < 0.1 * tol and 2 * np.pi / dz > 1000:
            dz *= 1.5
        elif e > tol:
            dz /= 3
        else:

            # store new point
            zz.append(z)
            rr.append(r)

    # finish back at the starting point
    zz[-1] = zz[0]
    rr[-1] = rr[0]

    xx = rr * np.cos(zz)
    yy = rr * np.sin(zz)

    return xx, yy
