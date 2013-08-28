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
# Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA

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
import tempfile
import re

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
    """Read a TimeSeries of channel data from a source.

    Acceptable sources are:
        - a .gwf-format framefile (string ending in '.gwf')
        - a LAL-format cache file (string ending in '.lcf' or '.cache')
        - a lal.Cache object (either from SWIG-LAL or GLUE)

    @param source
        input source, see above for details
    @param channel
        string name of channel, e.g. 'L1:LDAS-STRAIN'
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
    # read from single frame
    if isinstance(source, basestring) and source.endswith('.gwf'):
        return ts_from_frame_file(source, channel, start=start,
                                  duration=duration, datatype=datatype,
                                  verbose=verbose)
    # read from cache file
    elif (isinstance(source, basestring) and
          re.search('(.lcf|.cache)\Z', source)):
        cache = lal.CacheImport(os.path.expanduser(source))
        return ts_from_cache(cache, channel, start=start,
                             duration=duration, datatype=datatype,
                             verbose=verbose)
    # read from XLALCache
    elif isinstance(source, lal.Cache):
        return ts_from_cache(source, channel, start=start,
                             duration=duration, datatype=datatype,
                             verbose=verbose)
    # read from GLUE cache
    elif _HAS_GLUE and isinstance(source, gcache.Cache):
        cache = lalutils.lalcache_from_gluecache(source)
        return ts_from_cache(cache, channel, start=start,
                             duration=duration, datatype=datatype,
                             verbose=verbose)
    # otherwise barf
    else:
        raise ValueError("Cannot interpret source '%s'." % source)


def ts_from_cache(cache, channel, start=None, duration=None, datatype=None,
                  verbose=False):
    """Read a TimeSeries of channel data from a LAL Cache object

    @param cache
        XLALCAche() containing list of GWF file paths
    @param channel
        string name of channel, e.g. 'L1:LDAS-STRAIN'
    @param start
        LIGOTimeGPS start time for output TimeSeries
    @param duration
        float duration (seconds) for output TimeSeries
    @param datatype
        datatype, either an integer from the LALTYPECODE, a string
        matchine the corresponding type, or a numpy dtype

    @returns a TimeSeries of the imported data
    """
    # open the cache into a stream
    stream = lalframe.FrCacheOpen(cache)
    # read the stream
    return ts_from_stream(stream, channel, start=start, duration=duration,
                       datatype=datatype, verbose=verbose)


def ts_from_frame_file(framefile, channel, start=None, duration=None,
                       datatype=None, verbose=False):
    """Read a TimeSeries of channel data from a GWF-format framefile

    @param framefile
        path to GWF-format framefile to read
    @param channel
        string name of channel, e.g. 'L1:LDAS-STRAIN'
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
    stream = lalframe.FrOpen('', framefile)
    # read the stream
    return ts_from_stream(stream, channel, start=start, duration=duration,
                       datatype=datatype, verbose=verbose)


def ts_from_stream(stream, channel, start=None, duration=None, datatype=None,
                verbose=False):
    """Read a TimeSeries of channel data from an open FrStream

    @param stream
        XLALFrStream() of data from which to read
    @param channel
        string name of channel, e.g. 'L1:LDAS-STRAIN'
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
    lalframe.FrSetMode(verbose and lalframe.LAL_FR_STREAM_VERBOSE_MODE or
                       lalframe.LAL_FR_STREAM_DEFAULT_MODE, stream)
    # determine default start time and duration
    epoch = stream.epoch
    if start is None:
        start = epoch
    if not duration:
        startoffset = float(start - epoch)
        duration = float(get_stream_duration(stream)) - startoffset

    # get series type
    frdatatype = lalframe.FrStreamGetTimeSeriesType(channel, stream)
    if datatype is None:
        datatype = frdatatype
    else:
        datatype = _parse_datatype(datatype)

    # read original data
    read = getattr(lalframe, 'FrStreamRead%sTimeSeries'
                             % lalutils.LAL_TYPE_STR[frdatatype])
    origin = read(stream, channel, start, duration, 0)
    # format to output data-type if required
    if datatype == frdatatype:
        series = origin
    if datatype != frdatatype:
        create = lalutils.func_factory('create', '%stimeseries'
                                              % lalutils.LAL_TYPE_STR[datatype])
        series = create(channel, start, origin.f0, origin.deltaT,
                        origin.sampleUnits, origin.data.length)
        series.data.data = origin.data.data.astype(
                               lalutils.NUMPY_TYPE_FROM_LAL[datatype])
    # rewind the stream and return
    lalframe.FrStreamSeek(stream, epoch)
    return series


def get_stream_length(stream, channel):
    """Find the number of samples represented in a frame stream

    @param stream
        XLALFrStream() of data to measure
    @param channel
        string name of channel to measure

    @returns the integer length of the data for this channel
    """
    epoch = lal.LIGOTimeGPS(stream.epoch)
    # loop over each file in the stream cache and query its vector length
    nfile = stream.cache.length
    length = 0
    for i in range(nfile):
        length += lalframe.FrStreamGetVectorLength(stream, channel)
        lalframe.FrStreamNext(stream)
    # rewind the stream and return
    lalframe.FrStreamSeek(stream, epoch)
    return length


def get_stream_duration(stream):
    """Find the duration of time stored in a frame stream

    @param stream
        XLALFrStream() of data to measure
    @param channel
        string name of channel to measure

    @returns the float duration (seconds) of the data for this channel
    """
    epoch = lal.LIGOTimeGPS(stream.epoch)
    # loop over each file in the stream cache and query its duration
    nfile = stream.cache.length
    duration = 0
    for i in range(nfile):
        duration += lalframe.FrFileQueryDt(stream.file, 0)
        lalframe.FrStreamNext(stream)
    # rewind stream and return
    lalframe.FrStreamSeek(stream, epoch)
    return duration


def _parse_datatype(datatype):
    """Internal helper to format an abitrary datatype reference into
    an integer element of the LALTYPECODE enum
    """
    # return already formatted int
    if isinstance(datatype, int) and lalutils.LAL_TYPE_STR.has_key(datatype):
        return datatype
    # convert from string, e.g. 'real8'
    elif (isinstance(datatype, basestring) and
         lalutils.LAL_TYPE_FROM_STR.has_key(datatype.upper())):
        return lalutils.LAL_TYPE_FROM_STR[datatype.upper()]
    # convert from numpy dtype
    elif lalutils.LAL_TYPE_FROM_NUMPY.has_key(datatype):
        return lalutils.LAL_TYPE_FROM_NUMPY[datatype]
    # otherwise barf
    raise ValueError("Cannot interpret datatype '%s'" % datatype)
