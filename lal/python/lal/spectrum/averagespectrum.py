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

"""This module provides functions to calculated the average spectrum (PSD)
of a given data set, using wrappers of the LAL FFT package
"""

import re
import warnings

try:
    from .. import lal
except ImportError:
    raise ImportError("The SWIG-wrappings of LAL cannot be imported.")

from .. import git_version
__author__ = "Duncan Macleod <duncan.macleod@ligo.org>"
__version__ = git_version.id
__date__ = git_version.date

from .. import utils

## \addtogroup lal_py_spectrum
#@{

def median_mean(series, segment_length, segment_overlap, window=None,
                plan=None, unit=lal.lalStrainUnit):
    """Computes the spectrum of the timeseries using the
    median-mean average method.

    For more details see XLALREAL8AverageSpectrumMedianMean().

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
    return _psd('medianmean', series, segment_length, segment_overlap,
                window=window, plan=plan, unit=unit)


def welch(series, segment_length, segment_overlap, window=None,
          plan=None, unit=lal.lalStrainUnit):
    """Computes the spectrum of the timeseries using the Welch
    average method.

    For more details see XLALREAL8AverageSpectrumWelch().

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
    return _psd('welch', series, segment_length, segment_overlap,
                window=window, plan=plan, unit=unit)


def median(series, segment_length, segment_overlap, window=None,
                plan=None, unit=lal.lalStrainUnit):
    """Computes the spectrum of the timeseries using the median
    average method.

    For more details see XLALREAL8AverageSpectrumMedian().

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
    return _psd('median', series, segment_length, segment_overlap,
                window=window, plan=plan, unit=unit)


def bartlett(series, segment_length, window=None,
             plan=None, unit=lal.lalStrainUnit):
    """Computes the spectrum of the timeseries using the Bartlett
    average method.

    For more details see XLALREAL8AverageSpectrumWelch().

    @param series
        input TimeSeries
    @param segment_length
        number of samples for single Fourier transform
    @param window
        time-domain window to apply to the data (default: 24-point Kaiser)
    @param plan
        FFT memory plan for Fourier calculations
    @param unit
        output unit for FrequencySeries

    @returns a FrequencySeries with the same data-type as the input
    """
    return _psd('welch', series, segment_length, 0,
                window=window, plan=plan, unit=unit)

##@}

def _psd(method, series, segment_length, segment_overlap, window=None,
         plan=None, unit=lal.lalStrainUnit):
    """Computes the spectrum of the timeseries using the
    median average method.
    """
    method = method.lower()
    stype = utils.dtype(series)
    segment_length = int(segment_length)
    segment_overlap = int(segment_overlap)

    # check data length
    N = series.data.length
    numsegs = 1 + int((N - segment_length) / segment_overlap)
    if re.match('median(.*)mean', method) and numsegs % 2 == 1:
        numsegs -= 1
    required_N = int((numsegs - 1) * segment_overlap + segment_length)
    if N != required_N:
        warnings.warn("Data array is the wrong size for the correct number "
                      "of averages given the input parameters. The trailing %d "
                      "samples will not be used in this calculation."
                      % (N - required_N))
        resize = utils.func_factory('resize', series)
        series = utils.duplicate(series)
        timeseries = resize(series, 0, required_N)

    # generate window
    destroywindow = not window
    destroyplan   = not plan
    if not window:
        func   = getattr(lal, "CreateKaiser%sWindow" % stype)
        window = func(segment_length, 24)
    # generate FFT plan
    if not plan:
        func = getattr(lal, "CreateForward%sFFTPlan" % stype)
        plan = func(segment_length, 1)

    # generate output spectrum
    f0 = (1/series.deltaT) * (1/segment_length)
    deltaF = (1/segment_length)
    create = utils.func_factory('create', '%sfrequencyseries' % stype)
    spectrum = create(series.name, series.epoch, f0, deltaF, lal.lalStrainUnit,
                      segment_length//2+1)

    # calculate medianmean spectrum
    if re.match('median?mean\Z', method, re.I):
        average_spectrum = getattr(lal, "%sAverageSpectrumMedianMean" % stype)
    elif re.match('median\Z', method, re.I):
        average_spectrum = getattr(lal, "%sAverageSpectrumMedian" % stype)
    elif re.match('welch\Z', method, re.I):
        average_spectrum = getattr(lal, "%sAverageSpectrumWelch" % stype)
    else:
        raise NotImplementedError("Sorry, only 'median' and 'medianmean' "+\
                                  "and 'welch' average methods are available.")
    average_spectrum(spectrum, series, segment_length, segment_overlap,
                     window, plan)
    return spectrum
