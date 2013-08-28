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

# map LAL type codes to strings (for function factory)
LAL_TYPE_STR = {lal.LAL_I2_TYPE_CODE: 'INT2',
                lal.LAL_I4_TYPE_CODE: 'INT4',
                lal.LAL_I8_TYPE_CODE: 'INT8',
                lal.LAL_U2_TYPE_CODE: 'UINT2',
                lal.LAL_U4_TYPE_CODE: 'UINT4',
                lal.LAL_U8_TYPE_CODE: 'UINT8',
                lal.LAL_S_TYPE_CODE:  'REAL4',
                lal.LAL_D_TYPE_CODE:  'REAL8',
                lal.LAL_C_TYPE_CODE:  'COMPLEX8',
                lal.LAL_Z_TYPE_CODE:  'COMPLEX16'}

# map strings to LAL type codes
LAL_TYPE_FROM_STR = dict((v,k) for k, v in LAL_TYPE_STR.iteritems())

# map numpy dtypes to LAL type codes
LAL_TYPE_FROM_NUMPY = {numpy.int16:lal.LAL_I2_TYPE_CODE,
                       numpy.int32:lal.LAL_I4_TYPE_CODE,
                       numpy.int64:lal.LAL_I8_TYPE_CODE,
                       numpy.uint16:lal.LAL_U2_TYPE_CODE,
                       numpy.uint32:lal.LAL_U4_TYPE_CODE,
                       numpy.uint64:lal.LAL_U8_TYPE_CODE,
                       numpy.float32:lal.LAL_S_TYPE_CODE,
                       numpy.float64:lal.LAL_D_TYPE_CODE,
                       numpy.complex64:lal.LAL_C_TYPE_CODE,
                       numpy.complex128:lal.LAL_Z_TYPE_CODE}

# map LAL type codes to numpy dtypes
NUMPY_TYPE_FROM_LAL = dict((v,k) for k, v in LAL_TYPE_FROM_NUMPY.iteritems())


VALID_OPERATIONS = ['create', 'destroy', 'cut', 'resize', 'shrink', 'add']
_REGEX_TYPE = re.compile('(%s)' % '|'.join(LAL_TYPE_FROM_STR.keys()), re.I)
def func_factory(operation, dtype, series=None):
    """Returns the LAL function to perform the given operation for the
    relevant data type.

    Example::

        >>> create = func_factory(create, 'real8timeseries')
        >>> create
        lal.CreateREAL8TimeSeries
        >>> ts = create(name, epoch, f0, deltaT, sampleUnits, length)
        >>> func_factory('resize', ts)
        lal.ResizeREAL8TimeSeries
    """
    # verify operation
    try:
        VALID_OPERATIONS.index(operation)
    except ValueError:
        raise ValueError("Operation '%s' not understood for LAL series. "
                         "Please select one of: %s"
                         % (operation, ", ".join(VALID_OPERATIONS)))
    # verify data type
    struct = _struct_name(dtype)
    return getattr(lal, ''.join([operation.title(), struct]))


_STRUCTS = ['Sequence', 'Vector']
_SERIES_TYPES = ['Time', 'Frequency']
_SERIES_REGEX = re.compile('\A(?P<dtype>(%s))(?P<stype>(%s))Series\Z'
                            % ('|'.join(LAL_TYPE_FROM_STR.keys()),
                               '|'.join(_SERIES_TYPES)), re.I)
_ARRAY_REGEX = re.compile('\A(?P<dtype>(%s))Array(?:(?P<dir>(L|V))?)'
                          % ('|'.join(LAL_TYPE_FROM_STR.keys())), re.I)
_STRUCT_REGEX = re.compile('\A(?P<dtype>(%s))(?P<struct>(%s))\Z'
                            % ('|'.join(LAL_TYPE_FROM_STR.keys()),
                               '|'.join(_STRUCTS)), re.I)
def _struct_name(series):
    """Format a structure name into the understood type for LAL

    Example::

        >>> _struct_name('real8timeseries')
        'REAL8TimeSeries'
    """
    if isinstance(series, basestring):
        typestr = series
    else:
        typestr = series.__class__.__name__

    match = _SERIES_REGEX.match(typestr)
    if match:
        match = match.groupdict()
        return '%s%sSeries' % (match['dtype'].upper(), match['stype'].title())
    match = _ARRAY_REGEX.match(typestr)
    if match:
        match = match.groupdict()
        return '%sArray%s' % (match['dtype'].upper(),
                              match['dir'] and match['dir'].upper() or '')
    match = _STRUCT_REGEX
    try:
        match = match.groupdict()
    except AttributeError:
        raise ValueError("Input '%s' cannot be parsed into LAL struct name"
                         % series)
    else:
        return '%s%s' % (match['dtype'].upper(), match['struct'].title())


_TYPE_REGEX = re.compile('\A(?P<dtype>(%s))'
                         % ('|'.join(LAL_TYPE_FROM_STR.keys())), re.I)
def dtype(series):
    """Find the type string name for this series

    Example::

        >> dtype(timeseries)
        'REAL8'
    """
    match = _TYPE_REGEX.match(series.__class__.__name__)
    try:
        match = match.groupdict()
    except AttributeError:
        raise ValueError("Data type for series type '%s' unknown."
                         % series.__class__.__name__)
    return match['dtype'].upper()


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
        out  = create(series.name, series.epoch, series.f0, series.deltaF,\
                      series.sampleUnits, series.data.length)
    elif stype.endswith('TimeSeries'):
        out  = create(series.name, series.epoch, series.f0, series.deltaT,\
                      series.sampleUnits, series.data.length)
    else:
        raise NotImplementedError("A duplicator for the %s has not been "
                                  "implemented: here's your chance!")
    out.data.data = series.data.data
    return out
