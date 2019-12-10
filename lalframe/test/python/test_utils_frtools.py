# -*- coding: utf-8 -*-
#
# Copyright (C) 2019 Duncan Macleod
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

"""Tests for lalframe.utils.frtools
"""

import sys
import os
try:
    from pathlib import Path
except ImportError as exc:  # probably macports
    import warnings
    warnings.warn(str(exc))
    sys.exit(77)

import pytest

from lalframe import utils

# find test GWF file
TEST_PATH = Path(os.getenv(
    "LAL_TEST_SRCDIR",
    Path(__file__).parent,
)).absolute().parent
TEST_GWF = TEST_PATH / "F-TEST-600000060-60.gwf"

# test result
TEST_CHANNELS = [
    "H1:LSC-AS_Q",
]


def test_get_channels():
    assert utils.get_channels(TEST_GWF) == TEST_CHANNELS


if __name__ == '__main__':
    if "-v" not in " ".join(sys.argv[1:]):  # default to verbose
        sys.argv.append("-v")
    sys.exit(pytest.main(args=[__file__] + sys.argv[1:]))
