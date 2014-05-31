#
# Copyright (C) 2013  Leo Singer
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
#
from __future__ import division
"""
Functions related to matched filtering.
"""
__author__ = "Leo Singer <leo.singer@ligo.org>"


# General imports
import numpy as np
import math
from scipy import optimize

# LAL imports
import lal
import lalsimulation

# My own imports
from .decorator import memoized


# Useful sample units
unitInverseHertz = lal.Unit()
unitInverseSqrtHertz = lal.Unit()
lal.UnitInvert(unitInverseHertz, lal.lalHertzUnit)
lal.UnitSqrt(unitInverseSqrtHertz, unitInverseHertz)


# Memoize FFT plans
CreateForwardCOMPLEX16FFTPlan = memoized(lal.CreateForwardCOMPLEX16FFTPlan)
CreateForwardREAL8FFTPlan = memoized(lal.CreateForwardREAL8FFTPlan)
CreateReverseCOMPLEX16FFTPlan = memoized(lal.CreateReverseCOMPLEX16FFTPlan)
CreateReverseREAL8FFTPlan = memoized(lal.CreateReverseREAL8FFTPlan)


def ceil_pow_2(n):
    """Return the least integer power of 2 that is greater than or equal to n.

    >>> ceil_pow_2(128.0)
    128.0
    >>> ceil_pow_2(0.125)
    0.125
    >>> ceil_pow_2(129.0)
    256.0
    >>> ceil_pow_2(0.126)
    0.25
    >>> ceil_pow_2(1.0)
    1.0
    """
    # frexp splits floats into mantissa and exponent, ldexp does the opposite.
    # For positive numbers, mantissa is in [0.5, 1.).
    mantissa, exponent = math.frexp(n)
    return math.ldexp(
        1 if mantissa >= 0 else float('nan'),
        exponent - 1 if mantissa == 0.5 else exponent
    )


def fftfilt(b, x):
    """Apply the FIR filter with coefficients b to the signal x, as if the filter's
    state is initially all zeros. The output has the same length as x."""

    # Zero-pad by at least (len(b) - 1).
    nfft = int(ceil_pow_2(len(x) + len(b) - 1))

    # Create FFT plans.
    forwardplan = CreateForwardCOMPLEX16FFTPlan(nfft, 0)
    reverseplan = CreateReverseCOMPLEX16FFTPlan(nfft, 0)

    # Create temporary workspaces.
    workspace1 = lal.CreateCOMPLEX16Vector(nfft)
    workspace2 = lal.CreateCOMPLEX16Vector(nfft)
    workspace3 = lal.CreateCOMPLEX16Vector(nfft)

    workspace1.data[:len(x)] = x
    workspace1.data[len(x):] = 0
    lal.COMPLEX16VectorFFT(workspace2, workspace1, forwardplan)
    workspace1.data[:len(b)] = b
    workspace1.data[len(b):] = 0
    lal.COMPLEX16VectorFFT(workspace3, workspace1, forwardplan)
    workspace2.data *= workspace3.data
    lal.COMPLEX16VectorFFT(workspace1, workspace2, reverseplan)

    # Return result with zero-padding stripped.
    return workspace1.data[:len(x)] / nfft


def abscissa(series):
    """Produce the independent variable for a lal TimeSeries or FrequencySeries."""
    try:
        delta = series.deltaT
        x0 = float(series.epoch)
    except AttributeError:
        delta = series.deltaF
        x0 = series.f0
    return x0 + delta * np.arange(len(series.data.data))


def colored_noise(epoch, duration, sample_rate, psd):
    """Generate a REAL8TimeSeries containing duration seconds of colored
    Gaussian noise at the given sample rate, with the start time given by epoch.
    psd should be an instance of REAL8FrequencySeries containing a discretely
    sample power spectrum with f0=0, deltaF=1/duration, and a length of
    ((duration * sample_rate) // 2 + 1) samples.
    """
    data_length = duration * sample_rate
    plan = CreateReverseREAL8FFTPlan(data_length, 0)
    x = lal.CreateREAL8TimeSeries(None, lal.LIGOTimeGPS(0), 0, 0,
        lal.lalDimensionlessUnit, data_length)
    xf = lal.CreateCOMPLEX16FrequencySeries(None, epoch, 0, 1 / duration,
        lal.lalDimensionlessUnit, data_length // 2 + 1)
    white_noise = (np.random.randn(len(xf.data.data)) + np.random.randn(len(xf.data.data)) * 1j)

    # On line 1288 of lal's AverageSpectrum.c, in the code comments for
    # XLALWhitenCOMPLEX8FrequencySeries, it says that according to the LAL
    # conventions a whitened frequency series should consist of bins whose
    # real and imaginary parts each have a variance of 1/2.
    white_noise /= np.sqrt(2)

    # The factor of sqrt(2 * psd.deltaF) comes from the value of 'norm' on
    # line 1362 of AverageSpectrum.c.
    xf.data.data = white_noise * np.sqrt(psd.data.data / (2 * psd.deltaF))

    # Detrend the data: no DC component.
    xf.data.data[0] = 0

    # Return to time domain.
    lal.REAL8FreqTimeFFT(x, xf, plan)

    # Copy over metadata.
    x.epoch = epoch
    x.sampleUnits = lal.lalStrainUnit

    # Done.
    return x


def add_quadrature_phase(rseries):
    rseries_len = len(rseries.data.data)
    cseries = lal.CreateCOMPLEX16FrequencySeries(rseries.name, rseries.epoch,
        rseries.f0 - rseries.deltaF * (rseries_len - 1), rseries.deltaF,
        rseries.sampleUnits, 2 * (rseries_len - 1))
    cseries.data.data[:rseries_len] = 0
    cseries.data.data[rseries_len:] = 2 * rseries.data.data[1:-1]
    return cseries


def matched_filter_real_fd(template, psd):
    fdfilter = lal.CreateCOMPLEX16FrequencySeries(template.name, template.epoch,
        template.f0, template.deltaF, template.sampleUnits, len(template.data.data))
    fdfilter.data.data = template.data.data
    fdfilter = lal.WhitenCOMPLEX16FrequencySeries(fdfilter, psd)
    fdfilter.data.data /= np.sqrt(np.sum(np.abs(fdfilter.data.data)**2))
    fdfilter = lal.WhitenCOMPLEX16FrequencySeries(fdfilter, psd)
    return fdfilter


def matched_filter_spa(template, psd):
    """Create a complex matched filter kernel from a stationary phase approximation
    template and a PSD."""
    fdfilter = matched_filter_real_fd(template, psd)
    fdfilter2 = add_quadrature_phase(fdfilter)
    tdfilter = lal.CreateCOMPLEX16TimeSeries(None, lal.LIGOTimeGPS(0), 0, 0,
        lal.lalDimensionlessUnit, len(fdfilter2.data.data))
    plan = CreateReverseCOMPLEX16FFTPlan(len(fdfilter2.data.data), 0)
    lal.COMPLEX16FreqTimeFFT(tdfilter, fdfilter2, plan)
    return tdfilter


def matched_filter_real(template, psd):
    fdfilter = matched_filter_real_fd(template, psd)
    tdfilter = lal.CreateREAL8TimeSeries(None, lal.LIGOTimeGPS(0), 0, 0,
        lal.lalDimensionlessUnit, 2 * (len(fdfilter.data.data) - 1))
    plan = CreateReverseREAL8FFTPlan(len(tdfilter.data.data), 0)
    lal.REAL8FreqTimeFFT(tdfilter, fdfilter, plan)
    return tdfilter


def generate_template(mass1, mass2, S, f_low, sample_rate, template_duration, approximant, amplitude_order, phase_order):
    template_length = sample_rate * template_duration
    if approximant == lalsimulation.TaylorF2:
        zf, _ = lalsimulation.SimInspiralChooseFDWaveform(0,
            1 / template_duration,
            mass1 * lal.LAL_MSUN_SI, mass2 * lal.LAL_MSUN_SI,
            0, 0, 0, 0, 0, 0, f_low, 0, 1e6 * lal.LAL_PC_SI,
            0, 0, 0, None, None, amplitude_order, phase_order, approximant)
        lal.ResizeCOMPLEX16FrequencySeries(zf, 0, template_length // 2 + 1)

        # Generate over-whitened template
        psd = lal.CreateREAL8FrequencySeries(None, zf.epoch, zf.f0, zf.deltaF,
            lal.lalDimensionlessUnit, len(zf.data.data))
        psd.data.data = S(abscissa(psd))
        zW = matched_filter_spa(zf, psd)
    elif approximant == lalsimulation.TaylorT4:
        hplus, hcross = lalsimulation.SimInspiralChooseTDWaveform(
            0, 1 / sample_rate,
            mass1 * lal.LAL_MSUN_SI, mass2 * lal.LAL_MSUN_SI,
            0, 0, 0, 0, 0, 0,
            f_low, f_low,
            1e6 * lal.LAL_PC_SI,
            0, 0, 0,
            None, None,
            amplitude_order,
            phase_order,
            approximant)

        ht = lal.CreateREAL8TimeSeries(None, lal.LIGOTimeGPS(-template_duration), hplus.f0, hplus.deltaT, hplus.sampleUnits, template_length)
        hf = lal.CreateCOMPLEX16FrequencySeries(None, lal.LIGOTimeGPS(0), 0, 0, lal.lalDimensionlessUnit, template_length // 2 + 1)
        plan = CreateForwardREAL8FFTPlan(template_length, 0)

        ht.data.data[:-len(hplus.data.data)] = 0
        ht.data.data[-len(hplus.data.data):] = hplus.data.data
        lal.REAL8TimeFreqFFT(hf, ht, plan)

        psd = lal.CreateREAL8FrequencySeries(None, hf.epoch, hf.f0, hf.deltaF,
            lal.lalDimensionlessUnit, len(hf.data.data))
        psd.data.data = S(abscissa(psd))

        zWreal = matched_filter_real(hf, psd)

        ht.data.data[:-len(hcross.data.data)] = 0
        ht.data.data[-len(hcross.data.data):] = hcross.data.data

        lal.REAL8TimeFreqFFT(hf, ht, plan)
        zWimag = matched_filter_real(hf, psd)

        zW = lal.CreateCOMPLEX16TimeSeries(None, zWreal.epoch, zWreal.f0, zWreal.deltaT, zWreal.sampleUnits, len(zWreal.data.data))
        zW.data.data = zWreal.data.data + zWimag.data.data * 1j
    else:
        raise ValueError("unrecognized approximant")
    return zW.data.data[::-1].conj() * np.sqrt(2) * template_duration / sample_rate / 2


def abs2(y):
    """Return the |z|^2 for a complex number z."""
    return np.square(y.real) + np.square(y.imag)


#
# Lanczos interpolation
#


def lanczos(t, a):
    """The Lanczos kernel."""
    return np.sinc(t) * np.sinc(t / a)


def lanczos_interpolant(t, y):
    """An interpolant constructed by convolution of the Lanczos kernel with
    a set of discrete samples at unit intervals."""
    a = len(y) // 2
    return sum(lanczos(t - i + a, a) * yi for i, yi in enumerate(y))


def lanczos_interpolant_utility_func(t, y):
    """Utility function for Lanczos interpolation."""
    return -abs2(lanczos_interpolant(t, y))


def interpolate_max_lanczos(imax, y, window_length):
    """Find the time and maximum absolute value of a time series by Lanczos interpolation."""
    yi = y[imax-window_length:imax+window_length+1]
    tmax = optimize.fminbound(lanczos_interpolant_utility_func,
        -1., 1., (yi,), xtol=1e-5)
    tmax = np.asscalar(tmax)
    ymax = np.asscalar(lanczos_interpolant(tmax, yi))
    return imax + tmax, ymax


#
# Catmull-Rom spline interpolation
#


def poly_catmull_rom(y):
    return np.poly1d([
        -0.5 * y[0] + 1.5 * y[1] - 1.5 * y[2] + 0.5 * y[3],
        y[0] - 2.5 * y[1] + 2 * y[2] - 0.5 * y[3],
        -0.5 * y[0] + 0.5 * y[2],
        y[1]
    ])


def interpolate_max_catmull_rom_even(y):

    # Construct Catmull-Rom interpolating polynomials for real and imaginary parts
    poly_re = poly_catmull_rom(y.real)
    poly_im = poly_catmull_rom(y.imag)

    # Find the roots of d(|y|^2)/dt as approximated
    roots = (poly_re * poly_re.deriv() + poly_im * poly_im.deriv()).r

    # Find which of the two matched interior points has a greater magnitude
    t_max = 0.
    y_max = y[1]
    y_max_abs2 = abs2(y_max)

    new_t_max = 1.
    new_y_max = y[2]
    new_y_max_abs2 = abs2(new_y_max)

    if new_y_max_abs2 > y_max_abs2:
        t_max = new_t_max
        y_max = new_y_max
        y_max_abs2 = new_y_max_abs2

    # Find any real root in (0, 1) that has a magnitude greater than the greatest endpoint
    for root in roots:
        if np.isreal(root) and 0 < root < 1:
            new_t_max = root
            new_y_max = poly_re(new_t_max) + poly_im(new_t_max) * 1j
            new_y_max_abs2 = abs2(new_y_max)
            if new_y_max_abs2 > y_max_abs2:
                t_max = new_t_max
                y_max = new_y_max
                y_max_abs2 = new_y_max_abs2

    # Done
    return t_max, y_max


def interpolate_max_catmull_rom(imax, y, window_length):
    t_max, y_max = interpolate_max_catmull_rom_even(y[imax - 2:imax + 2])
    y_max_abs2 = abs2(y_max)
    t_max = t_max - 1

    new_t_max, new_y_max = interpolate_max_catmull_rom_even(y[imax - 1:imax + 3])
    new_y_max_abs2 = abs2(new_y_max)

    if new_y_max_abs2 > y_max_abs2:
        t_max = new_t_max
        y_max = new_y_max
        y_max_abs2 = new_y_max_abs2

    return imax + t_max, y_max


#
# Quadratic fit
#


def interpolate_max_quadratic_fit(imax, y, window_length):
    """Quadratic fit to absolute value of y. Note that this one does not alter the value at the maximum."""

    poly = np.polyfit(
        np.arange(-window_length, window_length + 1.),
        np.abs(y[imax - window_length:imax + window_length + 1]),
        2)

    # Find which of the two matched interior points has a greater magnitude
    t_max = -1.
    y_max = y[imax - 1]
    y_max_abs2 = abs2(y_max)

    new_t_max = 1.
    new_y_max = y[imax + 1]
    new_y_max_abs2 = abs2(new_y_max)

    if new_y_max_abs2 > y_max_abs2:
        t_max = new_t_max
        y_max = new_y_max
        y_max_abs2 = new_y_max_abs2

    # Determine if the global extremum of the polynomial is a local maximum in (-1, 1)
    A, B, C = poly
    new_t_max = -0.5 * B / A
    new_y_max_abs2 = np.square(np.polyval(poly, new_t_max))
    if -1 < new_t_max < 1 and new_y_max_abs2 > y_max_abs2:
        t_max = new_t_max

    return imax + t_max, y[imax]


#
# Nearest neighbor interpolation
#


def interpolate_max_nearest_neighbor(imax, y, window_length):
    """Trivial, nearest-neighbor interpolation"""
    return imax, y[imax]


#
# Set default interpolation scheme
#


def interpolate_max(imax, y, window_length, method='catmull-rom'):
    if method == 'catmull-rom':
        func = interpolate_max_catmull_rom
    elif method == 'lanczos':
        func = interpolate_max_lanczos
    elif method == 'nearest-neighbor':
        func = interpolate_max_nearest_neighbor
    elif method == 'quadratic-fit':
        func = interpolate_max_quadratic_fit
    else:
        raise ValueError('unrecognized method: %s' % method)
    return func(imax, y, window_length)
