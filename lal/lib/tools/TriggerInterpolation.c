/*
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

#include <complex.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <lal/TriggerInterpolate.h>
#include <lal/TriggerInterpolation.h>


#define APPLY_FUNC(NAME, TYPE) \
int XLAL ## TYPE ## Apply ## NAME ## TriggerInterpolant( \
    NAME ## TriggerInterpolant *interp, \
    double *tmax, \
    TYPE *ymax, \
    const TYPE *y) \
{ \
    int result; \
    double complex ymax_full; \
    double complex data_full[2 * interp->window + 1]; \
    double complex *const y_full = &data_full[interp->window]; \
    for (int i = -(int)interp->window; i <= (int)interp->window; i ++) \
        y_full[i] = y[i]; \
    result = XLALTriggerInterpolate ## NAME(tmax, &ymax_full, y_full, interp->window); \
    if (result == GSL_SUCCESS) \
        *ymax = ymax_full; \
    return result; \
}



#define LEGACY_API(NAME) \
struct tag ## NAME ## TriggerInterpolant { \
    unsigned int window; \
}; \
\
NAME ## TriggerInterpolant *XLALCreate ## NAME ## TriggerInterpolant(unsigned int window) \
{ \
    NAME ## TriggerInterpolant *interp = malloc(sizeof(NAME ## TriggerInterpolant)); \
    if (!interp) \
        GSL_ERROR_NULL("Failed to allocate interpolant", GSL_ENOMEM); \
    interp->window = window; \
    return interp; \
} \
\
void XLALDestroy ## NAME ## TriggerInterpolant(NAME ## TriggerInterpolant *interp) \
{ \
    free(interp); \
} \
\
APPLY_FUNC(NAME, COMPLEX16) \
APPLY_FUNC(NAME, COMPLEX8) \
APPLY_FUNC(NAME, REAL8) \
APPLY_FUNC(NAME, REAL4)

LEGACY_API(CubicSpline)
LEGACY_API(CubicSplineAmpPhase)
LEGACY_API(Lanczos)
LEGACY_API(NearestNeighbor)
LEGACY_API(QuadraticFit)
