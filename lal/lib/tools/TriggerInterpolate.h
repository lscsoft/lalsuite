
/*
 * Interpolation of time and SNR of triggers
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA
 *
 * Copyright (C) 2012-2020 Leo Singer
 */

#ifndef _TRIGGERINTERPOLATE_H
#define _TRIGGERINTERPOLATE_H

/* SWIG bindings are not correct for these functions. */
#ifndef SWIG

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif


/**
 * \defgroup TriggerInterpolation Trigger Interpolation
 * \ingroup lal_tools
 *
 * This module implements several algorithms for performing sub-sample
 * interpolation of peaks in the output of a matched filter. The interface for
 * each algorithm is the same. To illustrate, here is a description of the
 * interface for the 'QuadraticFit' algorithm.
 *
 * In this example, the interpolation will take into account a given sample,
 * the 5 samples before, and the 5 samples after, for a total of 11 samples.
 *
 * To apply the interpolant to a time series, call
 * \c XLALTriggerInterpolateQuadraticFit. It is requried that |y[0]| >= |y[i]|
 * for all i from -\c window to +\c window. Suppose that you have a complex
 * array \c y:
 *
 * \code{.c}
 * COMPLEX16 y[] = {...};
 * \endcode
 *
 * and suppose that the maximum of |y[i]| occurs at i = 8.
 * To interpolate the time and complex value at which the maximum occurs, do
 * the following:
 *
 * \code{.c}
 * double tmax;
 * COMPLEX16 ymax;
 * int result = XLALTriggerInterpolateQuadraticFit(interp, &tmax, &ymax, &y[8], 5);
 * \endcode
 *
 * Upon success, the return value is zero, the interpolated index of the
 * maximum sample is 8 + \c *tmax, and the interpolated value is in \c ymax.
 * Upon failure, the return value is nonzero and neither \c *tmax nor \c *ymax
 * are modified.
 *
 * The recommended algorithms are \c XLALTriggerInterpolateCubicSplineAmpPhase.
 * and \c XLALTriggerInterpolateQuadraticFit. The
 * \c XLALTriggerInterpolateCubicSpline and \c XLALTriggerInterpolateLanczos
 * methods should not be used because they tend to exhibit oscillations between
 * samples and may perform no better than the nearest-neighbor method; these
 * two functions are preserved only for historical interest. The
 * \c XLALTriggerInterpolateCubicSplineAmpPhase method is preferred over the
 * \c XLALTriggerInterpolateQuadraticFit method because the former provides a
 * continuous and smooth model for the signal at all times, whereas the latter
 * does not.
 *
 * \{
 */


/**
 * Catmull-Rom cubic spline interpolation on amplitude phase
 *
 * Find the Catmull-Rom cubic spline interpolants that pass through the
 * amplitude and phase of the data.
 *
 * \warning It is an error if the window size is less than 2.
 */
int XLALTriggerInterpolateCubicSplineAmpPhase(
    double *tmax,
    double _Complex *ymax,
    const double _Complex *y,
    unsigned int window);


/**
 * Catmull-Rom cubic spline interpolation
 *
 * Find the Catmull-Rom cubic spline interpolants that pass through the real and
 * imaginary parts of the data. Find the interpolated time by solving for the
 * maximum of the sum of the squares of the two spline interpolants.
 *
 * \warning It is an error if the window size is less than 2.
 */
int XLALTriggerInterpolateCubicSpline(
    double *tmax,
    double _Complex *ymax,
    const double _Complex *y,
    unsigned int window);


/**
 * Lanczos interpolation
 *
 * Convolve the data with a Lanczos reconstruction kernel and find the maximum
 * of the absolute value using a one-dimensional numerical method.
 *
 * \warning It is an error if the window size is less than 1.
 */
int XLALTriggerInterpolateLanczos(
    double *tmax,
    double _Complex *ymax,
    const double _Complex *y,
    unsigned int window);


/**
 * Nearest-neighbor interpolation
 *
 * Perform trivial, nearest-neighbor interpolation.
 * Always sets <tt>*tmax = 0, *ymax = data[0]</tt>.
 */
int XLALTriggerInterpolateNearestNeighbor(
    double *tmax,
    double _Complex *ymax,
    const double _Complex *y,
    unsigned int window);


/**
 * Quadratic fit
 *
 * Fit a quadratic function to the absolute value of the data. Report the sample
 * index of the vertex if the vertex is a maximum. Note that this is not
 * actually an interpolant because it is not guaranteed to agree with the input
 * data points, and that it always sets *ymax = data[0].
 *
 * \warning It is an error if the window size is less than 1.
 */
int XLALTriggerInterpolateQuadraticFit(
    double *tmax,
    double _Complex *ymax,
    const double _Complex *y,
    unsigned int window);


/** \} */


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
} /* extern "C" */
#endif

#endif /* ifndef SWIG */

#endif /* _TRIGGERINTERPOLATE_h */
