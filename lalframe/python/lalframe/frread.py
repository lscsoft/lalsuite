# Copyright (C) 2013 Duncan Macleod
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with with program; see the file COPYING. If not, write to the
# Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
# MA  02110-1301  USA

"""Wrappings of the LALFrame input methods for reading GWF files

This module provides the read_timeseries function, the primary method
by which users can load TimeSeries data from a variety of sources,
including:
    - a Gravitational-Wave Frame file (.gwf file extension)
    - a LAL-format cache file (.lcf or .cache file extension)
    - a LAL-format cache object (either XLALCache() or glue.lal.Cache [1]

[1] https://www.lsc-group.phys.uwm.edu/daswg/projects/glue/doc/glue.lal.Cache-class.html
"""

import os
import re

from six import string_types

from lal import (utils as lalutils, lal)

try:
    from glue import lal as gcache
except ImportError:
    _HAS_GLUE = False
else:
    _HAS_GLUE = True

from . import (lalframe, git_version)
__author__ = "Duncan Macleod <duncan.macleod@ligo.org>"
__version__ = git_version.id
__date__ = git_version.date

__all__ = ['read_timeseries']


def read_timeseries(source, channel, start=None, duration=None,
                    datatype=None, verbose=False):
    r"""Read a TimeSeries of channel data from a source.

    Acceptable sources are:
        - a .gwf-format framefile (string ending in '.gwf')
        - a LAL-format cache file (string ending in '.lcf' or '.cache')
        - a lal.Cache object (either from SWIG-LAL or GLUE)

    @param source
        input source, see above for details
    @param channel
        string name of channel, e.g. 'L1:LDAS-STRAIN', or a list of channel
        names
    @param start
        LIGOTimeGPS start time for output TimeSeries
    @param duration
        float duration (seconds) for output TimeSeries
    @param datatype
        datatype, either an integer from the LALTYPECODE, a string
        matchine the corresponding type, or a numpy dtype
    @param verbose
        print verbose output, default: False

    @returns a TimeSeries of the imported data

    Example 1, reading from a frame file:

    \code
    >>> out = read_timeseries('L-R-1061499968-32.gwf', 'L1:PSL-ISS_PDB_OUT_DQ')
    >>> print(type(out))
    <type 'REAL4TimeSeries'>
    >>> print(out.name, float(out.epoch), out.deltaT)
    ('L1:PSL-ISS_PDB_OUT_DQ', 1061499968.0, 3.0517578125e-05)
    \endcode

    Example 2, reading from a cache:

    \code
    >>> import lal
    >>> cache = lal.CacheGlob('/scratch/ER4/L0/L1/L-R-10614', 'L-R-1061499968*')
    >>> out = frread.read_timeseries(cache, 'L1:PSL-ISS_PDB_OUT_DQ')
    >>> print(out.name, float(out.epoch), out.deltaT)
    ('L1:PSL-ISS_PDB_OUT_DQ', 1061499968.0, 3.0517578125e-05)
    \endcode

    Example 3, restricting data input:

    \code
    >>> out = read_timeseries('L-R-1061499968-32.gwf', 'L1:PSL-ISS_PDB_OUT_DQ',
                              start=1061499970, duration=10)
    >>> print(out.name, float(out.epoch), out.data.length)
    ('L1:PSL-ISS_PDB_OUT_DQ', 1061499970.0, 327680)
    \endcode

    Example 4, specifying data type:

    \code
    >>> out = read_timeseries('L-R-1061499968-32.gwf',
                              'L1:PSL-ODC_CHANNEL_OUT_DQ')
    >>> print(type(out), out.data.data[:4])
    (<type 'REAL4TimeSeries'>,
     array([ 4259839.,  4259839.,  4259839.,  4259839.], dtype=float32))
    >>> out = read_timeseries('L-R-1061499968-32.gwf',
                              'L1:PSL-ODC_CHANNEL_OUT_DQ', datatype='int8')
    >>> print(type(out), out.data.data[:4])
    (<type 'INT8TimeSeries'>, array([4259839, 4259839, 4259839, 4259839]))
    \endcode

    """
    # parse channels
    if isinstance(channel, string_types):
        channels = [channel]
    else:
        channels = list(channel)

    # read cache file
    if (isinstance(source, string_types) and
          re.search(r'(.lcf|.cache)\Z', source)):
        source = lal.CacheImport(os.path.expanduser(source))
    # convert GLUE cache file
    if _HAS_GLUE and isinstance(source, gcache.Cache):
        source = lalutils.lalcache_from_gluecache(source)

    # read from single frame
    if isinstance(source, string_types) and source.endswith('.gwf'):
        out = _ts_from_frame_file(source, channels, start=start,
                                  duration=duration, datatype=datatype,
                                  verbose=verbose)
    # read from XLALCache
    elif isinstance(source, lal.Cache):
        out = _ts_from_cache(source, channels, start=start,
                             duration=duration, datatype=datatype,
                             verbose=verbose)
    # otherwise barf
    else:
        raise ValueError("Cannot interpret source '%s'." % source)

    # return
    if isinstance(channel, string_types):
        return out[0]
    else:
        return out


def _ts_from_cache(cache, channels, start=None, duration=None, datatype=None,
                   verbose=False):
    """Read a TimeSeries of channel data from a LAL Cache object

    @param cache
        XLALCache() containing list of GWF file paths
    @param channels
        list of channel names
    @param start
        LIGOTimeGPS start time for output TimeSeries
    @param duration
        float duration (seconds) for output TimeSeries
    @param datatype
        datatype, either an integer from the LALTYPECODE, a string
        matchine the corresponding type, or a numpy dtype
    @param verbose UNDOCUMENTED

    @returns a TimeSeries of the imported data
    """
    # open the cache into a stream
    stream = lalframe.FrCacheOpen(cache)
    # read the stream
    return _ts_from_stream(stream, channels, start=start, duration=duration,
                           datatype=datatype, verbose=verbose)


def _ts_from_frame_file(framefile, channels, start=None, duration=None,
                        datatype=None, verbose=False):
    """Read a TimeSeries of channel data from a GWF-format framefile

    @param framefile
        path to GWF-format framefile to read
    @param channels
        list of channel names
    @param start
        LIGOTimeGPS start time for output TimeSeries
    @param duration
        float duration (seconds) for output TimeSeries
    @param datatype
        datatype, either an integer from the LALTYPECODE, a string
        matchine the corresponding type, or a numpy dtype
    @param verbose
        print verbose output, default: False

    @returns a TimeSeries of the imported data
    """
    # open the file into a stream
    framefile = os.path.abspath(framefile)
    stream = lalframe.FrStreamOpen('', framefile)
    # read the stream
    return _ts_from_stream(stream, channels, start=start, duration=duration,
                           datatype=datatype, verbose=verbose)


def _ts_from_stream(stream, channels, start=None, duration=None, datatype=None,
                    verbose=False):
    """Read a TimeSeries of channel data from an open FrStream

    @param stream
        XLALFrStream() of data from which to read
    @param channels
        list of channel names
    @param start
        LIGOTimeGPS start time for output TimeSeries
    @param duration
        float duration (seconds) for output TimeSeries
    @param datatype
        datatype, either an integer from the LALTYPECODE, a string
        matchine the corresponding type, or a numpy dtype
    @param verbose
        print verbose output, default: False

    @returns a TimeSeries of the imported data
    """
    # set verbosity
    lalframe.FrSetMode(verbose and lalframe.FR_STREAM_VERBOSE_MODE or
                       lalframe.FR_STREAM_DEFAULT_MODE, stream)
    # determine default start time and duration
    epoch = lal.LIGOTimeGPS(stream.epoch)
    if start is None:
        start = epoch
    if not duration:
        startoffset = float(start - epoch)
        duration = float(get_stream_duration(stream)) - startoffset

    out = []
    for channel in channels:
        out.append(read_channel_from_stream(stream, channel, start, duration,
                                            datatype=datatype))
        lalframe.FrStreamSeek(stream, epoch)
    return out


def read_channel_from_stream(stream, channel, start, duration, datatype=None):
    """Read the TimeSeries of a single channel from an open stream
    """
    # get series type
    frdatatype = lalframe.FrStreamGetTimeSeriesType(channel, stream)
    if datatype is None:
        datatype = frdatatype
    else:
        datatype = lalutils.get_lal_type(datatype)

    # read original data
    read = getattr(lalframe, 'FrStreamRead%sTimeSeries'
                             % lalutils.get_lal_type_str(frdatatype))
    origin = read(stream, channel, start, duration, 0)
    # format to output data-type if required
    if datatype == frdatatype:
        return origin
    if datatype != frdatatype:
        create = lalutils.func_factory(
            'create', '%stimeseries' % lalutils.get_lal_type_str(datatype))
        series = create(channel, start, origin.f0, origin.deltaT,
                        origin.sampleUnits, origin.data.length)
        series.data.data = origin.data.data.astype(
                               lalutils.get_numpy_type(datatype))
        return series


def get_stream_length(stream, channel):
    """Find the number of samples represented in a frame stream

    @param stream
        XLALFrStream() of data to measure
    @param channel
        string name of channel to measure

    @returns the integer length of the data for this channel
    """
    epoch = lal.LIGOTimeGPS(stream.epoch.gpsSeconds,
                            stream.epoch.gpsNanoSeconds)
    # loop over each file in the stream cache and query its vector length
    nfile = stream.cache.length
    length = 0
    for i in range(nfile):
        for j in range(lalframe.FrFileQueryNFrame(stream.file)):
            length += lalframe.FrFileQueryChanVectorLength(stream.file,
                                                           channel,0)
            lalframe.FrStreamNext(stream)
    # rewind the stream and return
    lalframe.FrStreamSeek(stream, epoch)
    return length


def get_stream_duration(stream):
    """Find the duration of time stored in a frame stream

    @param stream
        XLALFrStream() of data to measure

    @returns the float duration (seconds) of the data for this channel
    """
    epoch = lal.LIGOTimeGPS(stream.epoch.gpsSeconds,
                            stream.epoch.gpsNanoSeconds)
    # loop over each file in the stream cache and query its duration
    nfile = stream.cache.length
    duration = 0
    for i in range(nfile):
        for j in range(lalframe.FrFileQueryNFrame(stream.file)):
            duration += lalframe.FrFileQueryDt(stream.file, 0)
            lalframe.FrStreamNext(stream)
    # rewind stream and return
    lalframe.FrStreamSeek(stream, epoch)
    return duration
