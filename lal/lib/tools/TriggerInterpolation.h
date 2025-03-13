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

#ifndef _TRIGGERINTERPOLATION_H
#define _TRIGGERINTERPOLATION_H

#include <lal/LALAtomicDatatypes.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

#ifdef SWIG
#define TRIGGER_INTERPOLATION_DEPRECATED
#else
#define TRIGGER_INTERPOLATION_DEPRECATED __attribute__ ((deprecated("TriggerInterpolation.h is deprecated. Use TriggerInterpolate.h instead.")));
#endif

#define LEGACY_API(NAME) \
struct tag ## NAME ## TriggerInterpolant; \
typedef struct tag ## NAME ## TriggerInterpolant NAME ## TriggerInterpolant; \
NAME ## TriggerInterpolant *XLALCreate ## NAME ## TriggerInterpolant(unsigned int window) TRIGGER_INTERPOLATION_DEPRECATED; \
void XLALDestroy ## NAME ## TriggerInterpolant(NAME ## TriggerInterpolant *) TRIGGER_INTERPOLATION_DEPRECATED; \
int XLALCOMPLEX16Apply ## NAME ## TriggerInterpolant(NAME ## TriggerInterpolant *, double *tmax, COMPLEX16 *ymax, const COMPLEX16 *y) TRIGGER_INTERPOLATION_DEPRECATED; \
int XLALCOMPLEX8Apply ## NAME ## TriggerInterpolant(NAME ## TriggerInterpolant *, double *tmax, COMPLEX8 *ymax, const COMPLEX8 *y) TRIGGER_INTERPOLATION_DEPRECATED; \
int XLALREAL8Apply ## NAME ## TriggerInterpolant(NAME ## TriggerInterpolant *, double *tmax, REAL8 *ymax, const REAL8 *y) TRIGGER_INTERPOLATION_DEPRECATED; \
int XLALREAL4Apply ## NAME ## TriggerInterpolant(NAME ## TriggerInterpolant *, double *tmax, REAL4 *ymax, const REAL4 *y) TRIGGER_INTERPOLATION_DEPRECATED;

LEGACY_API(CubicSpline)
LEGACY_API(CubicSplineAmpPhase)
LEGACY_API(Lanczos)
LEGACY_API(NearestNeighbor)
LEGACY_API(QuadraticFit)

#undef LEGACY_API
#undef TRIGGER_INTERPOLATION_DEPRECATED

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
} /* extern "C" */
#endif

#endif /* _TRIGGERINTERPOLATION_h */
