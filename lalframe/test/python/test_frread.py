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

"""Tests for lalframe.frread
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

import lal
from lal.utils.cache import CacheEntry

from lalframe import frread

try:
    from glue.lal import Cache
except ImportError:
    Cache = None
else:
    Cache.entry_class = CacheEntry

# find test GWF file
TEST_PATH = Path(os.getenv(
    "LAL_TEST_SRCDIR",
    Path(__file__).parent,
)).absolute().parent
TEST_GWF = TEST_PATH / "F-TEST-600000060-60.gwf"

# parametrize sources
SOURCES = [
    str(TEST_GWF),
    lal.CacheGlob(str(TEST_GWF.parent), TEST_GWF.name),
]
if Cache is not None:
    SOURCES.append(Cache([CacheEntry.from_T050017(str(TEST_GWF))]))


@pytest.mark.parametrize("source", SOURCES)
def test_read_timeseries(source):
    ts = frread.read_timeseries(
        source,
        "H1:LSC-AS_Q",
    )
    assert ts.name == "H1:LSC-AS_Q"
    assert ts.data.length == 16384 * 60
    assert ts.deltaT == 6.103515625e-05


@pytest.mark.parametrize("instart, induration, outstart, outduration", [
    (None, None, 600000060, 60),
    (600000061, None, 600000061, 59),
    (None, 30, 600000060, 30),
    (600000061, 1, 600000061, 1),
])
def test_read_timeseries_start_duration(
        instart,
        induration,
        outstart,
        outduration,
):
    ts = frread.read_timeseries(
        str(TEST_GWF),
        "H1:LSC-AS_Q",
        start=instart,
        duration=induration,
    )
    assert ts.epoch == outstart
    assert ts.data.length * ts.deltaT == outduration


@pytest.mark.parametrize("inputs, message", [
    (("does-not-exist.gwf", "channel"),
     "Internal function call failed: I/O error"),
    ((str(TEST_GWF), "bad-channel"), "Wrong name"),
])
def test_read_timeseries_error(inputs, message):
    with pytest.raises(RuntimeError) as exc:
        frread.read_timeseries(*inputs)
    assert str(exc.value) == message


if __name__ == '__main__':
    args = sys.argv[1:] or ["-v", "-rs", "--junit-xml=junit-frread.xml"]
    sys.exit(pytest.main(args=[__file__] + args))
