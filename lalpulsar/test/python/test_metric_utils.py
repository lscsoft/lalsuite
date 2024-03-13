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

import sys
import pytest

import numpy as np
from numpy.testing import assert_allclose


def test_mismatch_ellipse():
    from lalpulsar.metric_utils import mismatch_ellipse

    mu_max = 0.1

    with pytest.raises(AssertionError):
        g = np.ones((3, 2))
        xx, yy = mismatch_ellipse(g, mu_max)

    g = np.array([[7, 1], [1, 1]])

    with pytest.raises(AssertionError):
        xx, yy = mismatch_ellipse(g, 0.0)

    xx, yy = mismatch_ellipse(g, mu_max)

    assert len(xx) == len(yy)

    mu = []
    for i in range(len(xx)):
        xy = np.array([xx[i], yy[i]])
        mu.append(xy @ g @ xy)

    assert len(mu) == len(xx)

    assert_allclose(mu, mu_max, rtol=1e-5)


if __name__ == "__main__":
    args = sys.argv[1:] or ["-v", "-rs", "--junit-xml=junit-metric_utils.xml"]
    sys.exit(pytest.main(args=[__file__] + args))
