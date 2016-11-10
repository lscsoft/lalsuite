/*
 * Copyright (C) 2013 Evan Ochsner and Will M. Farr
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

#ifndef _LALSIMINSPIRALWAVEFORMCACHE_H
#define _LALSIMINSPIRALWAVEFORMCACHE_H

#include <lal/LALSimInspiral.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * @defgroup LALSimInspiralWaveformCache_h Header LALSimInspiralWaveformCache.h
 * @ingroup lalsimulation_inspiral
 *
 * @brief Routines for saving previously-computed waveforms for reuse.
 *
 * @{
 */

/**
 * Stores previously-computed waveforms and parameters to take
 * advantage of approximant- and parameter-specific opportunities for
 * accelerating waveform computation.
 */
typedef struct
tagLALSimInspiralWaveformCacheOld {
    REAL8TimeSeries *hplus;
    REAL8TimeSeries *hcross;
    COMPLEX16FrequencySeries *hptilde;
    COMPLEX16FrequencySeries *hctilde;
    REAL8 phiRef;
    REAL8 deltaTF;
    REAL8 m1;
    REAL8 m2;
    REAL8 S1x;
    REAL8 S1y;
    REAL8 S1z;
    REAL8 S2x;
    REAL8 S2y;
    REAL8 S2z;
    REAL8 f_min;
    REAL8 f_ref;
    REAL8 f_max;
    REAL8 r;
    REAL8 i;
    REAL8 lambda1;
    REAL8 lambda2;
    LALSimInspiralWaveformFlags *waveFlags;
    LALSimInspiralTestGRParam *nonGRparams; /* Non-NULL pointers here are not allowed b/c it's impossible to know which fields are present */
    int amplitudeO;
    int phaseO;
    Approximant approximant;
    REAL8Sequence *frequencies;
} LALSimInspiralWaveformCacheOld;

typedef struct
tagLALSimInspiralWaveformCache {
    REAL8TimeSeries *hplus;
    REAL8TimeSeries *hcross;
    COMPLEX16FrequencySeries *hptilde;
    COMPLEX16FrequencySeries *hctilde;
    REAL8 phiRef;
    REAL8 deltaTF;
    REAL8 m1;
    REAL8 m2;
    REAL8 S1x;
    REAL8 S1y;
    REAL8 S1z;
    REAL8 S2x;
    REAL8 S2y;
    REAL8 S2z;
    REAL8 f_min;
    REAL8 f_ref;
    REAL8 f_max;
    REAL8 r;
    REAL8 i;
    LALDict *LALpars;
    Approximant approximant;
    REAL8Sequence *frequencies;
} LALSimInspiralWaveformCache;

/** @} */

LALSimInspiralWaveformCache *XLALCreateSimInspiralWaveformCache(void);

void XLALDestroySimInspiralWaveformCache(LALSimInspiralWaveformCache *cache);

int XLALSimInspiralChooseTDWaveformFromCache(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 f_min, REAL8 f_ref, REAL8 r, REAL8 i, LALDict *LALpars, Approximant approximant, LALSimInspiralWaveformCache *cache);

int XLALSimInspiralChooseFDWaveformFromCache(COMPLEX16FrequencySeries **hptilde, COMPLEX16FrequencySeries **hctilde, REAL8 phiRef, REAL8 deltaF, REAL8 m1, REAL8 m2, REAL8 S1x, REAL8 S1y, REAL8 S1z, REAL8 S2x, REAL8 S2y, REAL8 S2z, REAL8 f_min, REAL8 f_max, REAL8 f_ref, REAL8 r, REAL8 i, LALDict *LALpars, Approximant approximant, LALSimInspiralWaveformCache *cache, REAL8Sequence *frequencies);

int XLALSimInspiralChooseFDWaveformSequence(COMPLEX16FrequencySeries **hptilde, COMPLEX16FrequencySeries **hctilde, REAL8 phiRef, REAL8 m1, REAL8 m2, REAL8 S1x, REAL8 S1y, REAL8 S1z, REAL8 S2x, REAL8 S2y, REAL8 S2z, REAL8 f_ref, REAL8 r, REAL8 i, LALDict *LALpars, Approximant approximant, REAL8Sequence *frequencies);

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMINSPIRAL_H */
