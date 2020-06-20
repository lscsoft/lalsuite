#!/usr/bin/env python
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

import datetime
import sys
try:
    from unittest import mock
except ImportError:  # python < 3
    import mock

import pytest

from freezegun import freeze_time

from lal import (LIGOTimeGPS, gpstime)


@mock.patch("lal.gpstime._gps_time_now", return_value=LIGOTimeGPS(100))
def test_gps_time_now(_):
    assert gpstime.gps_time_now() == 100.


@pytest.mark.parametrize('date, gps', [
    (datetime.datetime(2000, 1, 1, 0, 0, 0), 630720013),
    (datetime.date(2000, 1, 2), 630806413),
])
def test_utc_to_gps(date, gps):
    out = gpstime.utc_to_gps(date)
    assert isinstance(out, LIGOTimeGPS)
    assert out == LIGOTimeGPS(gps)


@pytest.mark.parametrize('date, gps', [
    (630720013, datetime.datetime(2000, 1, 1, 0, 0, 0), ),
])
def test_gps_to_utc(gps, date):
    assert gpstime.gps_to_utc(
        630720013,
    ) == datetime.datetime(2000, 1, 1, 0, 0, 0)


@mock.patch("lal.gpstime._gps_time_now", return_value=LIGOTimeGPS(100))
def test_utc_time_now(_):
    assert gpstime.utc_time_now() == datetime.datetime(
        1980, 1, 6, 0, 1, 40,
    )


@freeze_time("2015-09-14 09:50:45.391")
@mock.patch("lal.gpstime._gps_time_now", return_value=LIGOTimeGPS(1126259462))
@pytest.mark.parametrize('text, gps', [
    (None, 1126259462),
    ("now", 1126259462),
    ("today", 1126224017),
    ("tomorrow", 1126310417),
    ("yesterday", 1126137617),
    ("Sep 14 2015 09:50:45.391", LIGOTimeGPS(1126259462, 391000000)),
])
def test_str_to_gps(_, text, gps):
    out = gpstime.str_to_gps(text)
    print(type(out))
    assert isinstance(out, LIGOTimeGPS)
    assert out == LIGOTimeGPS(gps)


@pytest.mark.parametrize('gps, text', [
    (LIGOTimeGPS(1126259462, 391000000),
     "September 14 2015, 09:50:45.391000 UTC"),
    (1126259462, "September 14 2015, 09:50:45 UTC"),
])
def test_gps_to_str(gps, text):
    assert gpstime.gps_to_str(gps) == text


def test_gps_to_str_form():
    assert gpstime.gps_to_str(12345, form="%Y%m%d") == "19800106"


@mock.patch("lal.gpstime._gps_time_now", return_value=LIGOTimeGPS(1126259462))
@pytest.mark.parametrize('in_, result', [
    ("now", LIGOTimeGPS(1126259462)),
    (LIGOTimeGPS(1126259462, 391000000),
     "September 14 2015, 09:50:45.391000 UTC"),
])
def test_tconvert(_, in_, result):
    out = gpstime.tconvert(in_)
    assert type(out) == type(result)
    assert out == result


if __name__ == '__main__':
    args = sys.argv[1:] or ["-v", "-rs", "--junit-xml=junit-gpstime.xml"]
    sys.exit(pytest.main(args=[__file__] + args))
