# Copyright (C) 2012 Duncan M. Macleod
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

"""
This module provides a bunch of user-friendly wrappers to the HDF5
data format, for reading and writing SWIG-bound TimeSeries and
FrequencySeries objects from LAL. Issues should be raised on redmine:
https://bugs.ligo.org/redmine/projects/lalsuite
"""

# =============================================================================
# Preamble
# =============================================================================

from __future__ import division
import numpy
import h5py
import lal

from pylal import seriesutils
from pylal import git_version

# set metadata
__author__ = "Duncan M. Macleod <duncan.macleod@ligo.org>"
__version__ = git_version.id
__date__ = git_version.date

# =============================================================================
# Useful global variables
# =============================================================================

# type code dict
_numpy_lal_typemap = {numpy.int16: lal.LAL_I2_TYPE_CODE,\
                      numpy.int32: lal.LAL_I4_TYPE_CODE,\
                      numpy.int64: lal.LAL_I8_TYPE_CODE,\
                      numpy.uint16: lal.LAL_U2_TYPE_CODE,\
                      numpy.uint32: lal.LAL_U4_TYPE_CODE,\
                      numpy.uint64: lal.LAL_U8_TYPE_CODE,\
                      numpy.float32: lal.LAL_S_TYPE_CODE,\
                      numpy.float64: lal.LAL_D_TYPE_CODE,\
                      numpy.complex64: lal.LAL_C_TYPE_CODE,\
                      numpy.complex128: lal.LAL_Z_TYPE_CODE}
_lal_numpy_typemap = dict((v,k) for k, v in _numpy_lal_typemap.items())

# =============================================================================
# Data reading
# =============================================================================

def readTimeSeries(h5file, name, group="", start=None, duration=None,\
                   datatype=None):
    """
    Read a 1-D array from an HDF5 file into a LAL TimeSeries.

    Returns SWIG-bound LAL TimeSeries.

    Arguments:

        h5file : [ h5py.File | str ]
            open HDF5 file object, or path to HDF5 file on disk.
        name : str
            name of data object in HDF5 file relative in it's group.

    Keyword arguments:

        group : str
            name of HDF5 Group containing data object required.
        start : LIGOTimeGPS
            GPS start time of data requested, defaults to first data point.
        duration : float
            length of data (seconds) requested, default to all data.
        datatype : int
            LAL typecode for output datatype, defaults to type of data found.
    """
    own = False
    if not isinstance(h5file, h5py.File):
       h5file = h5py.File(h5file, "r")
       own = True
    
    # read data
    dataset, metadata = readArray(h5file, name, group=group)    
    if own:
        h5file.close()

    # parse metadata
    epoch = lal.LIGOTimeGPS(metadata.pop("epoch", 0))
    f0 = float(metadata.pop("f0", 0))
    deltaT = float(metadata.pop("dx", 0))
    sampleUnits = lal.lalDimensionlessUnit # FIXME

    # get series type
    if datatype is None:
        datatype = _numpy_lal_typemap[dataset.dtype.type]
    TYPECODE = seriesutils._typestr[datatype]
    numpytype = _lal_numpy_typemap[datatype]

    # cut data to size
    if duration is None:
        length = dataset.size
    else:
        length = int(duration / deltaT)
    if start is None:
        startidx = 0
    else:
        start = float(epoch) + (float(start)-float(epoch))//deltaT*deltaT
        if start < epoch:
            start = epoch
        startidx = int(max(0, (float(start)-float(epoch))/deltaT))
        epoch = lal.LIGOTimeGPS(start)
    endidx = min(dataset.size, startidx + length)
    length = endidx-startidx

    # build series
    createTimeSeries = getattr(lal, "Create%sTimeSeries" % TYPECODE)
    series = createTimeSeries(name, epoch, f0, deltaT, sampleUnits,\
                              length)
    series.data.data = dataset[startidx:endidx].astype(numpytype)

    return series

def readFrequencySeries(h5file, name, group=None, fmin=None, fmax=None,\
                        datatype=None):
    """
    Read FrequencySeries object from HDF5 file.

    Returns SWIG-bound LAL FrequencySeries.

    Arguments:

        h5file : [ h5py.File | str ]
            open HDF5 file object, or path to HDF5 file on disk.
        name : str
            name of data object in HDF5 file relative in it's group.

    Keyword arguments:

        group : str
            name of HDF5 Group containing data object required.
        fmin : float
            lower frequency bound on data returned
        fmax : float
            upper frequency bound on data returned
        datatype : int
            LAL typecode for output datatype, defaults to type of data found.
    """
    own = False
    if not isinstance(h5file, h5py.File):
       h5file = h5py.File(h5file, "r")
       own = True

    # read data
    dataset, metadata =  readArray(h5file, name, group=group)
    if own:
        h5file.close()

    # parse metadata
    epoch = lal.LIGOTimeGPS(metadata.pop("epoch", 0))
    f0 = float(metadata.pop("f0", 0))
    deltaF = float(metadata.pop("dx", 0))
    sampleUnits = lal.lalDimensionlessUnit # FIXME

    # get series type
    if not datatype:
        datatype = _numpy_lal_typemap[dataset.dtype.type]
    TYPECODE = seriesutils._typestr[datatype]
    numpytype = _lal_numpy_typemap[datatype]

    # cut data to size
    if fmin is None:
        fmin = f0
    else:
        fmin = f0 + (fmin-f0)//deltaF*deltaF
    startidx = int(max(0, (fmin-f0)/deltaF))
    if fmax is None:
        fmax = f0 + dataset.size * deltaF
    else:
        fmax = f0 + (fmax-f0)//deltaF*deltaF
    endidx = int(min(dataset.size, (fmax-fmin)/deltaF))
    length = endidx-startidx
    
    # build series
    createFrequencySeries = getattr(lal, "Create%sFrequencySeries" % TYPECODE)
    series = createFrequencySeries(name, epoch, f0, deltaF, sampleUnits,\
                                   length)
    series.data.data = dataset[startidx:endidx].astype(numpytype)

    return series

def readVectorSequence(h5file, name, group=None, start=None, duration=None,\
                       fmin=None, fmax=None, datatype=None):
    """
    Read VectorSequence object from HDF5 file.

    Returns SWIG-bound LAL VectorSequence, LIGOTimeGPS epoch, float deltaT,
    float f0, float deltaF.

    Arguments:

        h5file : [ h5py.File | str ]
            open HDF5 file object, or path to HDF5 file on disk.
        name : str
            name of data object in HDF5 file relative in it's group.

    Keyword arguments:

        group : str
            name of HDF5 Group containing data object required.
        start : LIGOTimeGPS
            GPS start time of data requested, defaults to first data point.
        duration : float
            length of data (seconds) requested, default to all data.
        fmin : float
            lower frequency bound on data returned
        fmax : float
            upper frequency bound on data returned
        datatype : int
            LAL typecode for output datatype, defaults to type of data found.
    """
    own = False
    if not isinstance(h5file, h5py.File):
       h5file = h5py.File(h5file, "r")
       own = True

    # read data
    dataset, metadata =  readArray(h5file, name, group=group)
    if own:
        h5file.close()

    # parse metadata
    epoch = lal.LIGOTimeGPS(metadata.pop("epoch", 0))
    deltaT = float(metadata.pop("dx", 0))
    f0 = float(metadata.pop("f0", 0))
    deltaF = float(metadata.pop("dy", 0))

    # get series type
    if not datatype:
        datatype = _numpy_lal_typemap[dataset.dtype.type]
    TYPECODE = seriesutils._typestr[datatype]
    numpytype = _lal_numpy_typemap[datatype]

    # cut data to size
    if duration is None:
        length = dataset.shape[0]
    else:
        length = int(duration / deltaT)
    if start is None:
        x0 = 0
    else:
        start = float(epoch) + (float(start)-float(epoch))//deltaT*deltaT
        if start < epoch:
            start = epoch
        x0 = int(max(0, (float(start)-float(epoch))/deltaT))
        epoch = lal.LIGOTimeGPS(start)
    length = min(dataset.size, x0 + length) - x0
    if fmin is None:
        fmin = f0
    else:
        fmin = f0 + (fmin-f0)//deltaF*deltaF
    y0 = int(max(0, (fmin-f0)/deltaF))
    if fmax is None:
        fmax = f0 + dataset.shape[1] * deltaF
    else:
        fmax = f0 + ((fmax-f0)//deltaF)*deltaF
    vectorLength = int(min(dataset.size, (fmax-fmin)/deltaF)) - y0

    # build series
    createVectorSequence = getattr(lal, "Create%sVectorSequence" % TYPECODE)
    sequence = createVectorSequence(length, vectorLength)

    for i,x in enumerate(numpy.arange(length)+x0):
        # get correct array
        sequence.data[i] = dataset[x,:][y0:y0+vectorLength].astype(numpytype)
        
    return sequence, epoch, deltaT, f0, deltaF

    series.data.data = dataset[startidx:endidx].astype(numpytype)

    return series

def readArray(h5file, name, group=None):
    """
    Read numpy.ndarray from h5py.File object h5file.

    Returns numpy.ndarray and metadata dict. 

    Arguments:

        h5file : [ h5py.File | str ]
            open HDF5 file object, or path to HDF5 file on disk.
        name : str
            name of data object in HDF5 file relative in it's group.

    Keyword arguments:

        group : str
            name of HDF5 Group containing data object required.
    """
    # format key
    if group.endswith(name):
        group = group[:len(name)]
    key = "%s/%s" % (group, name)

    # get dataset
    dataset = h5file[key]

    return dataset[...], dict(dataset.attrs)

# =============================================================================
# Data writing
# =============================================================================

def writeTimeSeries(h5file, series, group=None, **kwargs):
    """
    Write TimeSeries object to the given h5py.File object h5file.

    No return.

    Arguments:

        h5file : h5py.File
            open HDF5 file object
        series : TimeSeries
            SWIG-bound LAL TimeSeries object

    Keyword arguments:

        group : str
            name of HDF5 Group to write to, group is generated if required.

    All other keyword arguments are passed to
    h5py.highlevel.Group.create_dataset.
    """
    h5group = _create_groups(h5file, group)
    seriesToDataset(h5group, series, **kwargs)
    return

def writeFrequencySeries(h5file, series, group=None, **kwargs):
    """
    Write FrequencySeries object to the given h5py.File object h5file.

    No return.

    Arguments:

        h5file : h5py.File
            open HDF5 file object
        series : FrequencySeries
            SWIG-bound LAL FrequencySeries object

    Keyword arguments:

        group : str
            name of HDF5 Group to write to, group is generated if required.

    All other keyword arguments are passed to
    h5py.highlevel.Group.create_dataset.
    """
    h5group = _create_groups(h5file, group)
    seriesToDataset(h5group, series, **kwargs)
    return

def seriesToDataset(h5group, series, **kwargs):
    """
    Write {Time,Frequency}Series object to a new data set in the
    given h5py.Group object h5group.

    Returns h5py.highlevel.Dataset.

    Arguments:

        h5group : h5py.highlevel.Group
            Group object inside an HDF5 file
        series : [ TimeSeries | FrequencySeries ]
            LAL Time or FrequencySeries object

    All kwargs are passed to h5py.highlevel.Group.create_dataset.
    """
    # extract metadata
    metadata = dict()
    metadata["name"] = str(series.name)
    metadata["epoch"] = float(series.epoch)
    metadata["f0"] = float(series.f0)
    metadata["dx"] = float(hasattr(series, "deltaT") and series.deltaT\
                     or series.deltaF)
    # FIXME metadata["sampleUnits"] = series.sampleUnits
    
    # format data
    data = series.data.data
    
    # create dataset
    h5dataset = arrayToDataset(h5group, metadata["name"], data, metadata,\
                                   **kwargs) 
    
    return h5dataset

def arrayToDataset(h5group, name, array, metadata={}, **kwargs):
    """
    Write the given numpy array object to a new data set in the given
    h5py.highlevel.Group object h5group, including each element in the
    metadata dict.

    Returns h5py.highlevel.Dataset.

    Arguments:

        h5group : h5py.highlevel.Group
            Group object inside an HDF5 file
        name : str
            descriptive name for the new dataset.
        array: numpy.ndarray
            data object to write to dataset
        
    Keyword arguments:

        metadata: dict 
            dict containing metadata to write as attributes in the
            dataset, default={} (empty)

    All kwargs are passed to h5py.highlevel.Group.create_dataset.
    """
    # create dataset
    h5dataset = h5group.create_dataset(name, data=array, **kwargs)
    
    # create metadata
    metadata = dict(metadata)
    for key,value in metadata.iteritems():
        h5dataset.attrs.create(key, value)

    return h5dataset

def _create_groups(h5file, group):
    """
    Helper to return a group from an h5py.File object, creating all parents
    and the group itself if needed.

    Returns a new h5py.highlevel.Group, except in the special case of
    group=None, in which case the original h5py.File is returned.

    Argument:

        h5file : [ h5py.File | str ]
            open HDF5 file object, or path to HDF5 file on disk.
        group : str
            name of Group to create. If a composite group is given
            (e.g. "/data/spectra") the full tree will be generated
            as appropriate.
    """
    if group is None:
        return h5file
    group = group.lstrip("/")
    if group in h5file.keys():
        h5group = h5file[group]
    else:
        groups = group.split("/")
        h5group = h5file
        for g in groups:
            try:
                h5group = h5group.create_group(g)
            except ValueError:
                pass
    return h5group
