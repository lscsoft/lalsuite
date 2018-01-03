# Copyright (C) 2012 Duncan Macleod
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

"""This module provides methods to compute the Rayleigh spectrum of a
data set
"""

import numpy

try:
    from .. import lal
except ImportError:
    raise ImportError("The SWIG-wrappings of LAL cannot be imported.")

from .. import git_version
__author__ = "Duncan Macleod <duncan.macleod@ligo.org>"
__version__ = git_version.id
__date__ = git_version.date

__all__ = ['rayleigh']

## \addtogroup lal_py_spectrum
#@{


def rayleigh(series, segment_length, segment_overlap, window=None, plan=None,
             unit=lal.lalDimensionlessUnit):
    """Compute the Rayleigh spectrum of the data time-series

    @param series
        input TimeSeries
    @param segment_length
        number of samples for single Fourier transform
    @param segment_overlap
        number of samples between successive Fourier transforms
    @param window
        time-domain window to apply to the data (default: 24-point Kaiser)
    @param plan
        FFT memory plan for Fourier calculations
    @param unit
        output unit for FrequencySeries

    @returns a FrequencySeries with the same data-type as the input
    """
    # check data length
    N = series.data.length
    numsegs = 1 + int((N - segment_length) / segment_overlap)
    required_N = (numsegs - 1) * segment_overlap + segment_length
    if N != required_N:
        warnings.warn("Data array is the wrong size for the correct number "
                      "of averages given the input parameters. The trailing %d "
                      "samples will not be used in this calculation."
                      % (N - required_N))
        resize = func_factory('resize', timeseries)
        timeseries = resize(timeseries, 0, required_N)

    # generate window
    destroywindow = not window
    destroyplan   = not plan
    if not window:
        func   = getattr(lal, "CreateKaiser%sWindow" % TYPESTR)
        window = func(seglen, 24)
    # generate FFT plan
    if not plan:
        func = getattr(lal, "CreateForward%sFFTPlan" % TYPESTR)
        plan = func(seglen, 1)

    # generate output spectrum
    stype = dtype(timeseries)
    f0 = (1/series.deltaT) * (1/seglen)
    deltaF = (1/seglen)
    create = func_factory('create', '%sfrequencyseries' % stype)
    spectrum = create(series.name, series.epoch, f0, deltaF, lal.lalStrainUnit,
                      seglen//2+1)

    # compute Rayleigh spectrum average
    cut = func_factory('cut', series)
    data = numpy.ndarray((numsegs, segment_length//2+1))
    for i in range(numsegs):
        segseries = cut(series, i*segment_overlap,
                                i*segment_overlap + segment_length)
        data[i,:] = spectrum.bartlett(segseries, segment_length).data.data
    std = data.std(axis=0)
    mean = data.mean(axis=0)
    spectrum.data.data = std / mean
    return spectrum

##@}
