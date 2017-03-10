# Copyright (C) 2013 Duncan Macleod
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

"""Methods to generate LAL Time- and FrequencySeries objects in python
"""

import numpy
import re

try:
    from .. import lal
except ImportError:
    raise ImportError("The SWIG-wrappings of LAL cannot be imported.")

from .. import git_version
__author__ = "Duncan Macleod <duncan.macleod@ligo.org>"
__version__ = git_version.id
__date__ = git_version.date

# -----------------------------------------------------------------------------
# utility constants

# map LAL type codes to strings (for function factory)
LAL_TYPE_STR = {
    lal.I2_TYPE_CODE: 'INT2',
    lal.I4_TYPE_CODE: 'INT4',
    lal.I8_TYPE_CODE: 'INT8',
    lal.U2_TYPE_CODE: 'UINT2',
    lal.U4_TYPE_CODE: 'UINT4',
    lal.U8_TYPE_CODE: 'UINT8',
    lal.S_TYPE_CODE: 'REAL4',
    lal.D_TYPE_CODE: 'REAL8',
    lal.C_TYPE_CODE: 'COMPLEX8',
    lal.Z_TYPE_CODE: 'COMPLEX16',
}
LAL_TYPE_FROM_STR = dict((v, k) for k, v in LAL_TYPE_STR.items())
LAL_TYPE_STR_REGEX = re.compile(
    '(?P<dtype>(%s))' % ('|'.join(LAL_TYPE_FROM_STR.keys())), re.I)

# map numpy dtypes to LAL type codes
LAL_TYPE_FROM_NUMPY = {
    numpy.int16: lal.I2_TYPE_CODE,
    numpy.int32: lal.I4_TYPE_CODE,
    numpy.int64: lal.I8_TYPE_CODE,
    numpy.uint16: lal.U2_TYPE_CODE,
    numpy.uint32: lal.U4_TYPE_CODE,
    numpy.uint64: lal.U8_TYPE_CODE,
    numpy.float32: lal.S_TYPE_CODE,
    numpy.float64: lal.D_TYPE_CODE,
    numpy.complex64: lal.C_TYPE_CODE,
    numpy.complex128: lal.Z_TYPE_CODE,
}
NUMPY_TYPE_FROM_LAL = dict((v, k) for k, v in LAL_TYPE_FROM_NUMPY.items())


# structure definers
SERIES_OPERATIONS = ['create', 'destroy', 'cut', 'resize', 'shrink', 'add']
SERIES_TYPES = ['Time', 'Frequency']
STRUCT_TYPES = ['Sequence', 'Vector']

SERIES_REGEX = re.compile(
    '%s(?P<stype>(%s))Series\Z'
    % (LAL_TYPE_STR_REGEX.pattern, '|'.join(SERIES_TYPES)), re.I)
ARRAY_REGEX = re.compile(
    '%sArray(?:(?P<dir>(L|V))?)' % LAL_TYPE_STR_REGEX.pattern, re.I)
STRUCT_REGEX = re.compile(
    '%s(?P<struct>(%s))\Z'
    % (LAL_TYPE_STR_REGEX.pattern, '|'.join(STRUCT_TYPES)), re.I)


# -----------------------------------------------------------------------------
# utility methods

def func_factory(operation, dtype):
    """Returns the LAL function to perform the given operation for the
    relevant data type.

    Example::

        >>> create = func_factory('create', 'real8timeseries')
        >>> create
        lal.CreateREAL8TimeSeries
        >>> ts = create(name, epoch, f0, deltaT, sampleUnits, length)
        >>> func_factory('resize', ts)
        lal.ResizeREAL8TimeSeries
    """
    # verify operation
    try:
        SERIES_OPERATIONS.index(operation.lower())
    except ValueError as e:
        e.args("Operation '%s' not understood for LAL series. "
               "Please select one of: %s"
               % (operation, ", ".join(SERIES_OPERATIONS)),)
        raise e
    # verify data type
    struct = get_struct_name(dtype)
    return getattr(lal, ''.join([operation.title(), struct]))


def get_struct_name(series):
    """Format a structure name into the understood type for LAL

    Example::

        >>> get_struct_name('real8timeseries')
        'REAL8TimeSeries'
    """
    # get name of object
    if isinstance(series, basestring):
        typestr = series
    else:
        typestr = type(series).__name__

    # attempt to match as a series type
    try:
        match = SERIES_REGEX.match(typestr).groupdict()
    except AttributeError:
        pass
    else:
        return '%s%sSeries' % (match['dtype'].upper(), match['stype'].title())

    # attempt to match as an array (with optional dimension)
    try:
        match = ARRAY_REGEX.match(typestr).groupdict()
    except AttributeError:
        pass
    else:
        return '%sArray%s' % (match['dtype'].upper(),
                              match['dir'] and match['dir'].upper() or '')

    # attempt to match as a structure
    try:
        match = STRUCT_REGEX.match(typestr).groupdict()
    except AttributeError:
        raise ValueError(
            "Input %s cannot be parsed into LAL struct name" % series)
    else:
        return '%s%s' % (match['dtype'].upper(), match['struct'].title())


def get_series_type(series):
    """Find the LAL type enum for this series

    @param series
        a LAL series object (e.g. REAL8TimeSeries)

    @returns the LAL type enum (integer) for the series
    """
    try:
        match = TYPE_REGEX.match(type(series).__name__).groupdict()
    except AttributeError:
        raise ValueError("Data type for series type %r unknown."
                         % type(series).__name__)
    else:
        return get_lal_type(match['dtype'].upper())


def get_lal_type_str(datatype):
    """Return the LAL type str for the given `datatype`

    @param datatype
        a dtype representation, normally a string, or a python/numpy type
        object

    @returns the LAL type str for the given datatype

    Example::

        >>> get_lal_type_str('uint32')
        'UINT4'
        >>> get_lal_type_str(float)
        'REAL8'
    """
    return LAL_TYPE_STR[get_lal_type(datatype)]


def get_lal_type(datatype):
    """Return the LAL type enum for the given `datatype`

    @param datatype
        a dtype representation, normally a string, or a python/numpy type
        object

    @returns the LAL type enum (integer) for the given datatype

    Example::

        >>> get_lal_type('uint32')
        34
        >>> get_lal_type(float)
        11
    """
    # parse a LAL type enum
    try:
        LAL_TYPE_STR[datatype]
    except KeyError:
        pass
    else:
        return datatype
    # map a LAL type string
    try:
        return LAL_TYPE_FROM_STR[datatype]
    except KeyError:
        pass
    # test again for 'real4' or 'real8' (lower-case)
    #     can't do this with others because they match numpy names
    if re.match('real(4|8)\Z', str(datatype), re.I):
        return LAL_TYPE_FROM_STR[datatype.upper()]
    # format as a numpy data type and parse
    try:
        dtype = numpy.dtype(datatype).type
    except TypeError:
        pass
    else:
        try:
            return LAL_TYPE_FROM_NUMPY[dtype]
        except KeyError as e:
            e.args = ('LAL has no support for numpy.%s' % dtype.__name__,)
            raise
    raise ValueError("Cannot interpret datatype %r" % datatype)


def get_numpy_type(datatype):
    """Return the numpy type for the given `datatype`

    @param datatype
        a dtype representation, normally a LAL type enum (int),
        a LAL type string, or a python/numpy type object

    @returns the numpy type corresponding to the given datatype

    Example::

        >>> get_numpy_type(float)
        numpy.float64
        >>> get_numpy_type('REAL8')
        numpy.float64
        >>> get_numpy_type(11)
        numpy.float64
    """
    try:
        return NUMPY_TYPE_FROM_LAL[get_lal_type(datatype)]
    except KeyError as e:
        e.args('numpy has no support for %s'
               % get_lal_type(str(get_lal_type_str(datatype))),)
        raise


def duplicate(series):
    """
    Duplicate a TimeSeries or FrequencySeries.

    Arguments:

        series : [ TimeSeries | FrequencySeries ]
            input series to duplicate
    """
    create = func_factory('create', series)
    stype = series.__class__.__name__
    if stype.endswith('FrequencySeries'):
        out = create(series.name, series.epoch, series.f0, series.deltaF,
                     series.sampleUnits, series.data.length)
    elif stype.endswith('TimeSeries'):
        out = create(series.name, series.epoch, series.f0, series.deltaT,
                     series.sampleUnits, series.data.length)
    else:
        raise NotImplementedError("A duplicator for the %s has not been "
                                  "implemented: here's your chance!" % series)
    out.data.data = series.data.data
    return out
