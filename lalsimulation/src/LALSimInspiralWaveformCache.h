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

/** Stores previously-computed waveforms and parameters to take
    advantage of approximant- and parameter-specific opportunities for
    accelerating waveform computation. */
typedef struct tagLALSimInspiralWaveformCache LALSimInspiralWaveformCache;

LALSimInspiralWaveformCache *XLALCreateSimInspiralWaveformCache(void);

void XLALDestroySimInspiralWaveformCache(LALSimInspiralWaveformCache *cache);

int XLALSimInspiralChooseTDWaveformFromCache(
    REAL8TimeSeries **hplus,    /**< +-polarization waveform */
    REAL8TimeSeries **hcross,   /**< x-polarization waveform */
    REAL8 phiRef,               /**< reference orbital phase (rad) */
    REAL8 deltaT,               /**< sampling interval (s) */
    REAL8 m1,                   /**< mass of companion 1 (kg) */
    REAL8 m2,                   /**< mass of companion 2 (kg) */
    REAL8 s1x,                  /**< x-component of the dimensionless spin of object 1 */
    REAL8 s1y,                  /**< y-component of the dimensionless spin of object 1 */
    REAL8 s1z,                  /**< z-component of the dimensionless spin of object 1 */
    REAL8 s2x,                  /**< x-component of the dimensionless spin of object 2 */
    REAL8 s2y,                  /**< y-component of the dimensionless spin of object 2 */
    REAL8 s2z,                  /**< z-component of the dimensionless spin of object 2 */
    REAL8 f_min,                /**< starting GW frequency (Hz) */
    REAL8 f_ref,                /**< reference GW frequency (Hz) */
    REAL8 r,                    /**< distance of source (m) */
    REAL8 i,                    /**< inclination of source (rad) */
    REAL8 lambda1,              /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
    REAL8 lambda2,              /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
    LALSimInspiralWaveformFlags *waveFlags, /**< Set of flags to control special behavior of some waveform families. Pass in NULL (or None in python) for default flags */
    LALSimInspiralTestGRParam *nonGRparams, /**< Linked list of non-GR parameters. Pass in NULL (or None in python) for standard GR waveforms */
    int amplitudeO,             /**< twice post-Newtonian amplitude order */
    int phaseO,                 /**< twice post-Newtonian phase order */
    Approximant approximant,    /**< post-Newtonian approximant to use for waveform production */
    LALSimInspiralWaveformCache *cache  /**< waveform cache structure; use NULL for no caching */
    );

int XLALSimInspiralChooseFDWaveformFromCache(
    COMPLEX16FrequencySeries **hptilde,         /**< FD plus polarization */
    COMPLEX16FrequencySeries **hctilde,         /**< FD cross polarization */
    REAL8 phiRef,                               /**< reference orbital phase (rad) */
    REAL8 deltaF,                               /**< sampling interval (Hz) */
    REAL8 m1,                                   /**< mass of companion 1 (kg) */
    REAL8 m2,                                   /**< mass of companion 2 (kg) */
    REAL8 S1x,                                  /**< x-component of the dimensionless spin of object 1 */
    REAL8 S1y,                                  /**< y-component of the dimensionless spin of object 1 */
    REAL8 S1z,                                  /**< z-component of the dimensionless spin of object 1 */
    REAL8 S2x,                                  /**< x-component of the dimensionless spin of object 2 */
    REAL8 S2y,                                  /**< y-component of the dimensionless spin of object 2 */
    REAL8 S2z,                                  /**< z-component of the dimensionless spin of object 2 */
    REAL8 f_min,                                /**< starting GW frequency (Hz) */
    REAL8 f_max,                                /**< ending GW frequency (Hz) */
    REAL8 r,                                    /**< distance of source (m) */
    REAL8 i,                                    /**< inclination of source (rad) */
    REAL8 lambda1,                              /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
    REAL8 lambda2,                              /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
    LALSimInspiralWaveformFlags *waveFlags,     /**< Set of flags to control special behavior of some waveform families. Pass in NULL (or None in python) for default flags */
    LALSimInspiralTestGRParam *nonGRparams, 	/**< Linked list of non-GR parameters. Pass in NULL (or None in python) for standard GR waveforms */
    int amplitudeO,                             /**< twice post-Newtonian amplitude order */
    int phaseO,                                 /**< twice post-Newtonian order */
    Approximant approximant,                    /**< post-Newtonian approximant to use for waveform production */
    LALSimInspiralWaveformCache *cache         /**< waveform cache structure; use NULL for no caching */
    );

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMINSPIRAL_H */
