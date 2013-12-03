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
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA  02111-1307  USA
 *
 * Copyright (C) 2012 Leo Singer
 */

#include <lal/LALAtomicDatatypes.h>

#ifndef _TRIGGERINTERPOLATION_H
#define _TRIGGERINTERPOLATION_H

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif


/**
 * This module implements several algorithms for performing sub-sample
 * interpolation of peaks in the output of a matched filter. The interface for
 * each algorithm is the same. To illustrate, here is a description of the
 * interface for the 'Lanczos' algorithm.
 *
 * There is an opaque data type called \c LanczosTriggerInterpolant that holds
 * workspaces for performing the interpolation. To create one, call the
 * function \c XLALCreateLanczosTriggerInterpolant as follows:
 *
 * \code{.c}
 * LanczosTriggerInterpolant *interp = XLALCreateLanczosTriggerInterpolant(5);
 * \endcode
 *
 * The '5' is the interpolation window. In this example, the interpolation will
 * take into account a given sample, the 5 samples before, and the 5 samples
 * after, for a total of 11 samples.
 *
 * To apply the interpolant to a time series, call
 * \c XLALApplyLanczosTriggerInterpolant. It is requried that |y[0]| >= |y[i]|
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
 * int result = XLALApplyLanczosTriggerInterpolant(interp, &tmax, &ymax, &y[8]);
 * \endcode
 *
 * Upon success, the return value of \c XLALApplyLanczosTriggerInterpolant is
 * zero, the interpolated index of the maximum sample is 8 + \c *tmax, and the
 * interpolated value is in \c ymax. Upon failure, the return value is nonzero
 * and neither \c *tmax nor \c *ymax are modified.
 *
 * When you are done, release all of the workspace resources associated with
 * the interpolant with:
 *
 * \code{.c}
 * XLALDestroyLanczosTriggerInterpolant(interp);
 * \endcode
 *
 * \addtogroup TriggerInterpolation
 * \{
 */


/**
 * Catmull-Rom cubic spline interpolation
 *
 * Find the Catmull-Rom cubic spline interpolants that pass through the real and
 * imaginary parts of the data. Find the interpolated time by solving for the
 * maximum of the sum of the squares of the two spline interpolants.
 *
 * \warning It is an error to request a Catmull-Rom workspace with a window size
 * not equal to 2.
 *
 * \addtogroup CubicSplineTriggerInterpolant
 * \{
 */

struct tagCubicSplineTriggerInterpolant;
typedef struct tagCubicSplineTriggerInterpolant CubicSplineTriggerInterpolant;

/**
 * Constructor.
 *
 * Allocate a new interpolation workspace and return a pointer to
 * it, or NULL if an error occurs.
 *
 * \warning It is an error to request a Catmull-Rom workspace with a window size
 * not equal to 2.
 */
CubicSplineTriggerInterpolant *XLALCreateCubicSplineTriggerInterpolant(unsigned int window);

/**
 * Destructor.
 *
 * If the workspace is non-NULL, deallocate it.
 */
void XLALDestroyCubicSplineTriggerInterpolant(CubicSplineTriggerInterpolant *);

/**
 * Perform interpolation using the allocated workspace.
 *
 * Perform interpolation around the matched-filter output data pointed to by
 * \c y. There should exist \c window samples before and \c window samples
 * after this pointer.
 *
 * On success, set \c *tmax to the interpolated time, \c *ymax to the
 * complex interpolated signal, and return 0. On failure, return a non-zero GSL
 * error code.
 */
int XLALCOMPLEX16ApplyCubicSplineTriggerInterpolant(
    CubicSplineTriggerInterpolant *interp,
    double *tmax,
    COMPLEX16 *ymax,
    const COMPLEX16 *y);

int XLALCOMPLEX8ApplyCubicSplineTriggerInterpolant(
    CubicSplineTriggerInterpolant *interp,
    double *tmax,
    COMPLEX8 *ymax,
    const COMPLEX8 *y);

int XLALREAL8ApplyCubicSplineTriggerInterpolant(
    CubicSplineTriggerInterpolant *interp,
    double *tmax,
    REAL8 *ymax,
    const REAL8 *y);

int XLALREAL4ApplyCubicSplineTriggerInterpolant(
    CubicSplineTriggerInterpolant *interp,
    double *tmax,
    REAL4 *ymax,
    const REAL4 *y);

/* \} */


/**
 * Lanczos interpolation
 *
 * Convolve the data with a Lanczos reconstruction kernel and find the maximum
 * of the absolute value using a one-dimensional numerical method.
 *
 * \addtogroup LanczosTriggerInterpolant
 * \{
 */

struct tagLanczosTriggerInterpolant;
typedef struct tagLanczosTriggerInterpolant LanczosTriggerInterpolant;

/**
 * Constructor.
 *
 * Allocate a new interpolation workspace and return a pointer to
 * it, or NULL if an error occurs.
 */
LanczosTriggerInterpolant *XLALCreateLanczosTriggerInterpolant(unsigned int window);

/**
 * Destructor.
 *
 * If the workspace is non-NULL, deallocate it.
 */
void XLALDestroyLanczosTriggerInterpolant(LanczosTriggerInterpolant *);

/**
 * Perform interpolation using the allocated workspace.
 *
 * Perform interpolation around the matched-filter output data pointed to by
 * \c y. There should exist \c window samples before and \c window samples
 * after this pointer.
 *
 * On success, set \c *tmax to the interpolated time, \c *ymax to the
 * complex interpolated signal, and return 0. On failure, return a non-zero GSL
 * error code.
 */
int XLALCOMPLEX16ApplyLanczosTriggerInterpolant(
    LanczosTriggerInterpolant *interp,
    double *tmax,
    COMPLEX16 *ymax,
    const COMPLEX16 *y);

int XLALCOMPLEX8ApplyLanczosTriggerInterpolant(
    LanczosTriggerInterpolant *interp,
    double *tmax,
    COMPLEX8 *ymax,
    const COMPLEX8 *y);

int XLALREAL8ApplyLanczosTriggerInterpolant(
    LanczosTriggerInterpolant *interp,
    double *tmax,
    REAL8 *ymax,
    const REAL8 *y);

int XLALREAL4ApplyLanczosTriggerInterpolant(
    LanczosTriggerInterpolant *interp,
    double *tmax,
    REAL4 *ymax,
    const REAL4 *y);

/* \} */


/**
 * Nearest-neighbor interpolation
 *
 * Perform trivial, nearest-neighbor interpolation.
 * Always sets <tt>*tmax = 0, *ymax = data[0]</tt>.
 *
 * \warning It is an error to request a nearest-neighbor workspace with a window
 * size not equal to 0.
 *
 * \addtogroup NearestNeighborTriggerInterpolant
 * \{
 */

struct tagNearestNeighborTriggerInterpolant;
typedef struct tagNearestNeighborTriggerInterpolant NearestNeighborTriggerInterpolant;

/**
 * Constructor.
 *
 * Allocate a new interpolation workspace and return a pointer to
 * it, or NULL if an error occurs.
 *
 * \warning It is an error to request a nearest-neighbor workspace with a window
 * size not equal to 0.
 */
NearestNeighborTriggerInterpolant *XLALCreateNearestNeighborTriggerInterpolant(unsigned int window);

/**
 * Destructor.
 *
 * If the workspace is non-NULL, deallocate it.
 */
void XLALDestroyNearestNeighborTriggerInterpolant(NearestNeighborTriggerInterpolant *);

/**
 * Perform interpolation using the allocated workspace.
 *
 * Perform interpolation around the matched-filter output data pointed to by
 * \c y. There should exist \c window samples before and \c window samples
 * after this pointer.
 *
 * On success, set \c *tmax to the interpolated time, \c *ymax to the
 * complex interpolated signal, and return 0. On failure, return a non-zero GSL
 * error code.
 */
int XLALCOMPLEX16ApplyNearestNeighborTriggerInterpolant(
    NearestNeighborTriggerInterpolant *interp,
    double *tmax,
    COMPLEX16 *ymax,
    const COMPLEX16 *y);

int XLALCOMPLEX8ApplyNearestNeighborTriggerInterpolant(
    NearestNeighborTriggerInterpolant *interp,
    double *tmax,
    COMPLEX8 *ymax,
    const COMPLEX8 *y);

int XLALREAL8ApplyNearestNeighborTriggerInterpolant(
    NearestNeighborTriggerInterpolant *interp,
    double *tmax,
    REAL8 *ymax,
    const REAL8 *y);

int XLALREAL4ApplyNearestNeighborTriggerInterpolant(
    NearestNeighborTriggerInterpolant *interp,
    double *tmax,
    REAL4 *ymax,
    const REAL4 *y);

/* \} */


/**
 * Quadratic fit
 *
 * Fit a quadratic function to the absolute value of the data. Report the sample
 * index of the vertex if the vertex is a maximum. Note that this is not
 * actually an interpolant because it is not guaranteed to agree with the input
 * data points, and that it always sets *ymax = data[0].
 *
 * \addtogroup QuadraticFitTriggerInterpolant
 * \{
 */

struct tagQuadraticFitTriggerInterpolant;
typedef struct tagQuadraticFitTriggerInterpolant QuadraticFitTriggerInterpolant;

/**
 * Constructor.
 *
 * Allocate a new interpolation workspace and return a pointer to
 * it, or NULL if an error occurs.
 */
QuadraticFitTriggerInterpolant *XLALCreateQuadraticFitTriggerInterpolant(unsigned int window);

/**
 * Destructor.
 *
 * If the workspace is non-NULL, deallocate it.
 */
void XLALDestroyQuadraticFitTriggerInterpolant(QuadraticFitTriggerInterpolant *);

/**
 * Perform interpolation using the allocated workspace.
 *
 * Perform interpolation around the matched-filter output data pointed to by
 * \c y. There should exist \c window samples before and \c window samples
 * after this pointer.
 *
 * On success, set \c *tmax to the interpolated time, \c *ymax to the
 * complex interpolated signal, and return 0. On failure, return a non-zero GSL
 * error code.
 */
int XLALCOMPLEX16ApplyQuadraticFitTriggerInterpolant(
    QuadraticFitTriggerInterpolant *interp,
    double *tmax,
    COMPLEX16 *ymax,
    const COMPLEX16 *y);

int XLALCOMPLEX8ApplyQuadraticFitTriggerInterpolant(
    QuadraticFitTriggerInterpolant *interp,
    double *tmax,
    COMPLEX8 *ymax,
    const COMPLEX8 *y);

int XLALREAL8ApplyQuadraticFitTriggerInterpolant(
    QuadraticFitTriggerInterpolant *interp,
    double *tmax,
    REAL8 *ymax,
    const REAL8 *y);

int XLALREAL4ApplyQuadraticFitTriggerInterpolant(
    QuadraticFitTriggerInterpolant *interp,
    double *tmax,
    REAL4 *ymax,
    const REAL4 *y);

/* \} */


/* \} */


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
} /* extern "C" */
#endif

#endif /* _TRIGGERINTERPOLATION_h */
