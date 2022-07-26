# Copyright (C) 2021  Cardiff University
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
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

## \defgroup lal_py_gpstime GPSTime
## \ingroup lal_python
"""Utilties for calculating and modifying GPS times using the LAL date
package
"""
#
# ### Synopsis ###
#
# ~~~
# from lal import gpstime
# ~~~
# This module wraps the LAL \ref Date_h module into Python providing
# functions to convert from the `datetime.datetime` objects from the
# Python standard library into LIGOTimeGPS numbers, and vice-versa.
# See in particular ::gps_to_utc and ::utc_to_gps.
#
# This module also provides a Python implementation of the `tconvert'
# module, allowing conversion from strings of dates and times into
# LIGOTimeGPS, and vice-versa. See in particular ::gps_to_str and
# ::str_to_gps.
#
# \author Duncan Macleod <duncan.macleod@ligo.org>
#@{

import datetime as _datetime
from dateutil.parser import parse as str_to_utc

from .lal import (
    LIGOTimeGPS,
    GPSTimeNow as _gps_time_now,
    GPSToUTC as _gps_to_utc,
    UTCToGPS as _utc_to_gps,
)
from . import git_version

__author__ = "Duncan Macleod <duncan.macleod@ligo.org>"
__version__ = git_version.verbose_msg
__date__ = git_version.date

TIME_ZONES = {
    "PST": -8*3600,
    "PDT": -7*3600,
    "CST": -6*3600,
    "CDT": -5*3600,
    "EST": -5*3600,
    "EDT": -4*3600,
    "GMT": 0,
    "UTC": 0,
    "BST": 1*3600,
    "CET": 1*3600,
    "CEST": 2*3600,
}
GPS_EPOCH = _datetime.datetime(1980, 1, 6, 0, 0, 0)
LAL_GPS_MAX = _datetime.datetime(2048, 1, 24, 3, 13, 55)


def _check_utc(utc_time):
    """Check whether this UTC datetime will convert to a valid LIGOTimeGPS
    """
    if utc_time < GPS_EPOCH:
        raise ValueError("Given UTC time is before the start of the GPS era.")
    if utc_time > LAL_GPS_MAX:
        raise ValueError("Given UTC time is too far in the future, "
                         "LAL will SegmentationFault.")
    return


def gps_time_now():
    """Get the current time in GPS seconds

    @returns a LIGOTimeGPS
    """
    return _gps_time_now()


def utc_to_gps(utc_time):
    """Convert the given `datetime.datetime` into a GPS time

    @returns a LIGOTimeGPS
    """
    if not isinstance(utc_time, _datetime.datetime):
        utc_time = _datetime.datetime.combine(utc_time, _datetime.time())
    _check_utc(utc_time)
    return LIGOTimeGPS(_utc_to_gps(utc_time.utctimetuple()))


def gps_to_utc(gps):
    """Convert a GPS time into a `datetime.datetime`

    @returns a Python `datetime.datetime` object in UTC
    """
    gps = LIGOTimeGPS(gps)
    dt = _datetime.datetime(
        *_gps_to_utc(gps.gpsSeconds)[:7]
    ).replace(microsecond=0)  # force microseconds to 0
    if gps.gpsNanoSeconds:
        return dt.replace(
            microsecond=int(round(gps.gpsNanoSeconds * 1e-3)),
        )
    return dt


def utc_time_now():
    """Get the current date and time in UTC

    @returns a Python `datetime.datetime` object in UTC
    """
    return gps_to_utc(_gps_time_now())


def str_to_gps(time_string=None):
    r"""Converts a date/time string into a GPS time.

    The following special words are permitted:
        - "now"
        - "today"
        - "yesterday"
        - "tomorrow"

    Example:
    \code
    >>> gpstime.str_to_gps("September 14 2011, 01:46:25")
    1000000000.000000000
    \endcode

    @returns a LIGOTimeGPS
    """
    if not time_string or time_string.lower() == "now":
        return _gps_time_now()
    elif time_string == "today":
        date = _datetime.date.today()
        return utc_to_gps(_datetime.datetime.combine(date, _datetime.time()))
    elif time_string == "tomorrow":
        today = _datetime.datetime.combine(_datetime.date.today(),
                                           _datetime.time())
        tomorrow = today + _datetime.timedelta(days=1)
        return utc_to_gps(tomorrow)
    elif time_string == "yesterday":
        today = _datetime.datetime.combine(_datetime.date.today(),
                                           _datetime.time())
        yesterday = today - _datetime.timedelta(days=1)
        return utc_to_gps(yesterday)
    # otherwise parse the string as a date/time
    utc = str_to_utc(time_string, tzinfos=TIME_ZONES)
    micro = utc.microsecond
    gps = utc_to_gps(utc)
    return gps + micro / 1000000.0


def gps_to_str(gps, form=None):
    r"""
    Convert a LIGOTimeGPS time object into a string.
    The output format can be given explicitly, but will default
    as shown in the example.

    Example:

    \code
    >>> gps_to_str(1000000000)
    'September 14 2011, 01:46:25 UTC'
    \endcode

    @returns a string with the given format.
    """
    gps = LIGOTimeGPS(gps)
    utc = gps_to_utc(gps)
    if gps.gpsNanoSeconds and not form:
        form = "%B %d %Y, %H:%M:%S.%f UTC"
    elif not form:
        form = "%B %d %Y, %H:%M:%S UTC"
    return utc.strftime(form)


def tconvert(arg=None, form=None):
    r"""Convert date/time strings to and from GPS times.
    If no argument is given, the current GPS time is returned.

    The following special words are permitted:
        - "now"
        - "today"
        - "yesterday"
        - "tomorrow"

    Example:

    \code
    >>> tconvert()
    1048275013.000000000
    >>> tconvert("January 6 1980 00:00:00")
    0.000000000
    >>> tconvert(1000000000)
    'September 14 2011, 01:46:25 UTC'
    \endcode

    @returns the LIGOTimeGPS of the given time string, OR, string
    representing the given GPS time
    """
    try:
        float(arg)
    except ValueError:
        return str_to_gps(arg)
    else:
        return gps_to_str(arg, form=form)


##@}


if __name__ == "__main__":
    now = tconvert()
    now_utc = gps_to_str(now)
    print("The date/time now is %s (%d)" % (now_utc, now))
