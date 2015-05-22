/*
 * Copyright (C) 2013 Evan Ochsner and Will M. Farr
 *   2014 Salvatore Vitale
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

#ifndef _LALSIMBURSTWAVEFORMCACHE_H
#define _LALSIMBURSTWAVEFORMCACHE_H

#include <lal/LALSimBurst.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * Stores previously-computed waveforms and parameters to take
 * advantage of approximant- and parameter-specific opportunities for
 * accelerating waveform computation.
 */
typedef struct
tagLALSimBurstWaveformCache {
    REAL8TimeSeries *hplus;
    REAL8TimeSeries *hcross;
    COMPLEX16FrequencySeries *hptilde;
    COMPLEX16FrequencySeries *hctilde;
    REAL8 deltaT;
    REAL8 deltaF;
    REAL8 f0;
    REAL8 q,tau;
    REAL8 f_min;
    REAL8 f_max;
    REAL8 hrss;
    REAL8 polar_angle;
    REAL8 polar_ecc;
    LALSimBurstExtraParam *extraParams;
    BurstApproximant approximant;
} LALSimBurstWaveformCache;


LALSimBurstWaveformCache *XLALCreateSimBurstWaveformCache(void);

void XLALDestroySimBurstWaveformCache(LALSimBurstWaveformCache *cache);
int XLALSimBurstChooseTDWaveformFromCache(
        REAL8TimeSeries **hplus,                /**< +-polarization waveform */
        REAL8TimeSeries **hcross,               /**< x-polarization waveform */
        REAL8 deltaT,                           /**< sampling interval (s) */
        REAL8 f0,                               /**< central frequency [Hz] */
        REAL8 q,                               /**< quality */
        REAL8 tau,                              /**< duration */
        REAL8 f_min,                            /**< starting GW frequency (Hz) */
        REAL8 f_max,                            /**< max GW frequency (Hz) */
        REAL8 hrss,                             /**< hrss */
        REAL8 polar_angle,                      /**< Together with polar_ecc, controls ratio plus/cross polarization (rad) */
        REAL8 polar_ecc,                        /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
        LALSimBurstExtraParam *extraParams, /**< Linked list of extra Parameters (includes alpha and phase). Pass in NULL (or None in python) to neglect */
        BurstApproximant approximant,           /**< Burst approximant to use for waveform production */
        LALSimBurstWaveformCache *cache      /**< waveform cache structure; use NULL for no caching */
        );

int XLALSimBurstChooseFDWaveformFromCache(
        COMPLEX16FrequencySeries **hptilde,     /**< +-polarization waveform */
        COMPLEX16FrequencySeries **hctilde,     /**< x-polarization waveform */
        REAL8 deltaF,                           /**< sampling interval (Hz) */
        REAL8 deltaT,                           /**< time step corresponding to consec */
        REAL8 f0,                               /**< central frequency (Hz) */
        REAL8 q,                                /**< Q (==sqrt(2) \pi f0 tau ) [dless]*/
        REAL8 tau,                              /**< Duration [s] */
        REAL8 f_min,                            /**< starting GW frequency (Hz) */
        REAL8 f_max,                            /**< ending GW frequency (Hz) (0 for Nyquist) */
        REAL8 hrss,                             /**< hrss [strain] */
        REAL8 polar_angle,                      /**< Polar_ellipse_angle as defined in the burst table. Together with polar_ellipse_eccentricity below will fix the ratio of + vs x aplitude. Some WFs uses a single parameter alpha for this. Alpha is passed through extraParams*/
        REAL8 polar_ecc,                        /**< See above */
        LALSimBurstExtraParam *extraParams, /**< Linked list of extra burst parameters. Pass in NULL (or None in python) to neglect these */
        BurstApproximant approximant ,                /**< Burst approximant  */
        LALSimBurstWaveformCache *cache      /**< waveform cache structure; use NULL for no caching */
        );
#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMINSPIRAL_H */
