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

#include <math.h>
#include <LALSimInspiralWaveformCache.h>
#include <lal/LALSimInspiral.h>
#include <lal/FrequencySeries.h>

/**
 * Bitmask enumerating which parameters have changed, to determine
 * if the requested waveform can be transformed from a cached waveform
 * or if it must be generated from scratch.
 */
typedef enum {
    NO_DIFFERENCE = 0,
    INTRINSIC = 1,
    DISTANCE = 2,
    PHI_REF = 4,
    INCLINATION = 8
} CacheVariableDiffersBitmask;

static CacheVariableDiffersBitmask CacheArgsDifferenceBitmask(
        LALSimInspiralWaveformCache *cache,
        REAL8 phiRef,
        REAL8 deltaTF,
        REAL8 m1,
        REAL8 m2,
        REAL8 S1x, REAL8 S1y, REAL8 S1z,
        REAL8 S2x, REAL8 S2y, REAL8 S2z,
        REAL8 f_min, REAL8 f_ref_max,
        REAL8 r,
        REAL8 i,
        REAL8 lambda1,
        REAL8 lambda2,
        LALSimInspiralWaveformFlags *waveFlags,
        LALSimInspiralTestGRParam *nonGRparams,
        int amplitudeO,
        int phaseO,
        Approximant approximant);

static int StoreTDHCache(LALSimInspiralWaveformCache *cache,
        REAL8TimeSeries *hplus,
        REAL8TimeSeries *hcross,
        REAL8 phiRef,
        REAL8 deltaT,
        REAL8 m1, REAL8 m2,
        REAL8 S1x, REAL8 S1y, REAL8 S1z,
        REAL8 S2x, REAL8 S2y, REAL8 S2z,
        REAL8 f_min, REAL8 f_ref_max,
        REAL8 r,
        REAL8 i,
        REAL8 lambda1, REAL8 lambda2,
        LALSimInspiralWaveformFlags *waveFlags,
        LALSimInspiralTestGRParam *nonGRparams,
        int amplitudeO,
        int phaseO,
        Approximant approximant);

static int StoreFDHCache(LALSimInspiralWaveformCache *cache,
        COMPLEX16FrequencySeries *hptilde,
        COMPLEX16FrequencySeries *hctilde,
        REAL8 phiRef,
        REAL8 deltaT,
        REAL8 m1, REAL8 m2,
        REAL8 S1x, REAL8 S1y, REAL8 S1z,
        REAL8 S2x, REAL8 S2y, REAL8 S2z,
        REAL8 f_min, REAL8 f_ref_max,
        REAL8 r,
        REAL8 i,
        REAL8 lambda1, REAL8 lambda2,
        LALSimInspiralWaveformFlags *waveFlags,
        LALSimInspiralTestGRParam *nonGRparams,
        int amplitudeO,
        int phaseO,
        Approximant approximant);

/**
 * Chooses between different approximants when requesting a waveform to be generated
 * Returns the waveform in the time domain.
 * The parameters passed must be in SI units.
 *
 * This version allows caching of waveforms. The most recently generated
 * waveform and its parameters are stored. If the next call requests a waveform
 * that can be obtained by a simple transformation, then it is done.
 * This bypasses the waveform generation and speeds up the code.
 */
int XLALSimInspiralChooseTDWaveformFromCache(
        REAL8TimeSeries **hplus,                /**< +-polarization waveform */
        REAL8TimeSeries **hcross,               /**< x-polarization waveform */
        REAL8 phiRef,                           /**< reference orbital phase (rad) */
        REAL8 deltaT,                           /**< sampling interval (s) */
        REAL8 m1,                               /**< mass of companion 1 (kg) */
        REAL8 m2,                               /**< mass of companion 2 (kg) */
        REAL8 S1x,                              /**< x-component of the dimensionless spin of object 1 */
        REAL8 S1y,                              /**< y-component of the dimensionless spin of object 1 */
        REAL8 S1z,                              /**< z-component of the dimensionless spin of object 1 */
        REAL8 S2x,                              /**< x-component of the dimensionless spin of object 2 */
        REAL8 S2y,                              /**< y-component of the dimensionless spin of object 2 */
        REAL8 S2z,                              /**< z-component of the dimensionless spin of object 2 */
        REAL8 f_min,                            /**< starting GW frequency (Hz) */
        REAL8 f_ref,                            /**< reference GW frequency (Hz) */
        REAL8 r,                                /**< distance of source (m) */
        REAL8 i,                                /**< inclination of source (rad) */
        REAL8 lambda1,                          /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
        REAL8 lambda2,                          /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
        LALSimInspiralWaveformFlags *waveFlags, /**< Set of flags to control special behavior of some waveform families. Pass in NULL (or None in python) for default flags */
        LALSimInspiralTestGRParam *nonGRparams, /**< Linked list of non-GR parameters. Pass in NULL (or None in python) for standard GR waveforms */
        int amplitudeO,                         /**< twice post-Newtonian amplitude order */
        int phaseO,                             /**< twice post-Newtonian order */
        Approximant approximant,                /**< post-Newtonian approximant to use for waveform production */
        LALSimInspiralWaveformCache *cache      /**< waveform cache structure; use NULL for no caching */
        )
{
    int status;
    size_t j;
    REAL8 phasediff, dist_ratio, incl_ratio_plus, incl_ratio_cross;
    REAL8 cosrot, sinrot;
    CacheVariableDiffersBitmask changedParams;

    // If nonGRparams are not NULL, don't even try to cache.
    if ( nonGRparams != NULL || (!cache) )
        return XLALSimInspiralChooseTDWaveform(hplus, hcross, phiRef, deltaT,
                m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref, r, i,
                lambda1, lambda2, waveFlags, nonGRparams, amplitudeO, phaseO,
                approximant);

    // Check which parameters have changed
    changedParams = CacheArgsDifferenceBitmask(cache, phiRef, deltaT,
            m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref, r, i,
            lambda1, lambda2, waveFlags, nonGRparams, amplitudeO,
            phaseO, approximant);

    // No parameters have changed! Copy the cached polarizations
    if( changedParams == NO_DIFFERENCE ) {
        *hplus = XLALCutREAL8TimeSeries(cache->hplus, 0,
                cache->hplus->data->length);
        if (*hplus == NULL) return XLAL_ENOMEM;
        *hcross = XLALCutREAL8TimeSeries(cache->hcross, 0,
                cache->hcross->data->length);
        if (*hcross == NULL) {
            XLALDestroyREAL8TimeSeries(*hplus);
            *hplus = NULL;
            return XLAL_ENOMEM;
        }

        return XLAL_SUCCESS;
    }

    // Intrinsic parameters have changed. We must generate a new waveform
    if( (changedParams & INTRINSIC) != 0 ) {
        status = XLALSimInspiralChooseTDWaveform(hplus, hcross, phiRef, deltaT,
                m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref, r, i,
                lambda1, lambda2, waveFlags, nonGRparams, amplitudeO, phaseO,
                approximant);
        if (status == XLAL_FAILURE) return status;

        // FIXME: Need to add hlms, dynamic variables, etc. in cache
        return StoreTDHCache(cache, *hplus, *hcross, phiRef, deltaT, m1, m2,
            S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref, r, i, lambda1, lambda2,
            waveFlags, nonGRparams, amplitudeO, phaseO, approximant);
    }

    // case 1: Precessing waveforms
    if( approximant == SpinTaylorT4 || approximant == SpinTaylorT2 ) {
        // If polarizations are not cached we must generate a fresh waveform
        // FIXME: Will need to check hlms and/or dynamical variables as well
        if( cache->hplus == NULL || cache->hcross == NULL) {
            status = XLALSimInspiralChooseTDWaveform(hplus, hcross, phiRef,
                    deltaT, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref,
                    r, i, lambda1, lambda2, waveFlags, nonGRparams,
                    amplitudeO, phaseO, approximant);
            if (status == XLAL_FAILURE) return status;

            // FIXME: Need to add hlms, dynamic variables, etc. in cache
            return StoreTDHCache(cache, *hplus, *hcross, phiRef, deltaT, m1, m2,
                    S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref, r, i,
                    lambda1, lambda2, waveFlags, nonGRparams,
                    amplitudeO, phaseO, approximant);
        }

        if( changedParams & INCLINATION ) {
            // FIXME: For now just treat as intrinsic parameter.
            // Will come back and put in transformation
            status = XLALSimInspiralChooseTDWaveform(hplus, hcross, phiRef,
                    deltaT, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref,
                    r, i, lambda1, lambda2, waveFlags, nonGRparams, amplitudeO,
                    phaseO, approximant);
            if (status == XLAL_FAILURE) return status;

            // FIXME: Need to add hlms, dynamic variables, etc. in cache
            return StoreTDHCache(cache, *hplus, *hcross, phiRef, deltaT, m1, m2,
                    S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref, r, i,
                    lambda1, lambda2, waveFlags, nonGRparams, amplitudeO,
                    phaseO, approximant);
        }
        if( changedParams & PHI_REF ) {
            // FIXME: For now just treat as intrinsic parameter.
            // Will come back and put in transformation
            status = XLALSimInspiralChooseTDWaveform(hplus, hcross, phiRef,
                    deltaT, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref,
                    r, i, lambda1, lambda2, waveFlags, nonGRparams, amplitudeO,
                    phaseO, approximant);
            if (status == XLAL_FAILURE) return status;

            // FIXME: Need to add hlms, dynamic variables, etc. in cache
            return StoreTDHCache(cache, *hplus, *hcross, phiRef, deltaT, m1, m2,
                    S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref, r, i,
                    lambda1, lambda2, waveFlags, nonGRparams, amplitudeO,
                    phaseO, approximant);
        }
        if( (changedParams & DISTANCE) != 0 ) {
            // Return rescaled copy of cached polarizations
            dist_ratio = cache->r / r;
            *hplus = XLALCreateREAL8TimeSeries(cache->hplus->name,
                    &(cache->hplus->epoch), cache->hplus->f0,
                    cache->hplus->deltaT, &(cache->hplus->sampleUnits),
                    cache->hplus->data->length);
            if (*hplus == NULL) return XLAL_ENOMEM;

            *hcross = XLALCreateREAL8TimeSeries(cache->hcross->name,
                    &(cache->hcross->epoch), cache->hcross->f0,
                    cache->hcross->deltaT, &(cache->hcross->sampleUnits),
                    cache->hcross->data->length);
            if (*hcross == NULL) {
                XLALDestroyREAL8TimeSeries(*hplus);
                *hplus = NULL;
                return XLAL_ENOMEM;
            }

            for (j = 0; j < cache->hplus->data->length; j++) {
                (*hplus)->data->data[j] = cache->hplus->data->data[j]
                        * dist_ratio;
                (*hcross)->data->data[j] = cache->hcross->data->data[j]
                        * dist_ratio;
            }
        }

        return XLAL_SUCCESS;
    }

    // case 2: Non-precessing, ampO = 0
    else if( amplitudeO==0 && (approximant==TaylorT1 || approximant==TaylorT2
                || approximant==TaylorT3 || approximant==TaylorT4
                || approximant==EOBNRv2 || approximant==SEOBNRv1) ) {
        // If polarizations are not cached we must generate a fresh waveform
        if( cache->hplus == NULL || cache->hcross == NULL) {
            status = XLALSimInspiralChooseTDWaveform(hplus, hcross, phiRef,
                    deltaT, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref,
                    r, i, lambda1, lambda2, waveFlags, nonGRparams,
                    amplitudeO, phaseO, approximant);
            if (status == XLAL_FAILURE) return status;

            // FIXME: Need to add hlms, dynamic variables, etc. in cache
            return StoreTDHCache(cache, *hplus, *hcross, phiRef, deltaT, m1, m2,
                    S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref, r, i,
                    lambda1, lambda2, waveFlags, nonGRparams,
                    amplitudeO, phaseO, approximant);
        }

        // Set transformation coefficients for identity transformation.
        // We'll adjust them depending on which extrinsic parameters changed.
        dist_ratio = incl_ratio_plus = incl_ratio_cross = cosrot = 1.;
        phasediff = sinrot = 0.;

        if( changedParams & PHI_REF ) {
            // Only 2nd harmonic present, so {h+,hx} rotates by 2*deltaphiRef
            phasediff = 2.*(phiRef - cache->phiRef);
            cosrot = cos(phasediff);
            sinrot = sin(phasediff);
        }
        if( changedParams & INCLINATION) {
            // Rescale h+, hx by ratio of new/old inclination dependence
            incl_ratio_plus = (1.0 + cos(i)*cos(i))
                    / (1.0 + cos(cache->i)*cos(cache->i));
            incl_ratio_cross = cos(i) / cos(cache->i);
        }
        if( changedParams & DISTANCE ) {
            // Rescale h+, hx by ratio of (1/new_dist)/(1/old_dist) = old/new
            dist_ratio = cache->r / r;
        }

        // Create the output polarizations
        *hplus = XLALCreateREAL8TimeSeries(cache->hplus->name,
                &(cache->hplus->epoch), cache->hplus->f0,
                cache->hplus->deltaT, &(cache->hplus->sampleUnits),
                cache->hplus->data->length);
        if (*hplus == NULL) return XLAL_ENOMEM;
        *hcross = XLALCreateREAL8TimeSeries(cache->hcross->name,
                &(cache->hcross->epoch), cache->hcross->f0,
                cache->hcross->deltaT, &(cache->hcross->sampleUnits),
                cache->hcross->data->length);
        if (*hcross == NULL) {
            XLALDestroyREAL8TimeSeries(*hplus);
            *hplus = NULL;
            return XLAL_ENOMEM;
        }

        // Get new polarizations by transforming the old
        incl_ratio_plus *= dist_ratio;
        incl_ratio_cross *= dist_ratio;
        // FIXME: Do changing phiRef and inclination commute?!?!
        for (j = 0; j < cache->hplus->data->length; j++) {
            (*hplus)->data->data[j] = incl_ratio_plus
                    * (cosrot*cache->hplus->data->data[j]
                    - sinrot*cache->hcross->data->data[j]);
            (*hcross)->data->data[j] = incl_ratio_cross
                    * (sinrot*cache->hplus->data->data[j]
                    + cosrot*cache->hcross->data->data[j]);
        }

        return XLAL_SUCCESS;
    }
    // case 3: Non-precessing, ampO > 0
    // FIXME: EOBNRv2HM actually ignores ampO. If it's given with ampO==0,
    // it will fall to the catch-all and not be cached.
    else if( (amplitudeO==-1 || amplitudeO>0) && (approximant==TaylorT1
                || approximant==TaylorT2 || approximant==TaylorT3
                || approximant==TaylorT4 || approximant==EOBNRv2HM) ) {
        // If polarizations are not cached we must generate a fresh waveform
        // FIXME: Add in check that hlms non-NULL
        if( cache->hplus == NULL || cache->hcross == NULL) {
            // FIXME: This will change to a code-path: inputs->hlms->{h+,hx}
            status = XLALSimInspiralChooseTDWaveform(hplus, hcross, phiRef,
                    deltaT, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref,
                    r, i, lambda1, lambda2, waveFlags, nonGRparams,
                    amplitudeO, phaseO, approximant);
            if (status == XLAL_FAILURE) return status;

            // FIXME: Need to add hlms, dynamic variables, etc. in cache
            return StoreTDHCache(cache, *hplus, *hcross, phiRef, deltaT, m1, m2,
                    S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref, r, i,
                    lambda1, lambda2, waveFlags, nonGRparams,
                    amplitudeO, phaseO, approximant);
        }

        if( changedParams & INCLINATION) {
            // FIXME: For now just treat as intrinsic parameter.
            // Will come back and put in transformation
            status = XLALSimInspiralChooseTDWaveform(hplus, hcross, phiRef,
                    deltaT, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref,
                    r, i, lambda1, lambda2, waveFlags, nonGRparams, amplitudeO,
                    phaseO, approximant);
            if (status == XLAL_FAILURE) return status;

            // FIXME: Need to add hlms, dynamic variables, etc. in cache
            return StoreTDHCache(cache, *hplus, *hcross, phiRef, deltaT, m1, m2,
                    S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref, r, i,
                    lambda1, lambda2, waveFlags, nonGRparams, amplitudeO,
                    phaseO, approximant);

        }
        if( changedParams & PHI_REF ) {
            // FIXME: For now just treat as intrinsic parameter.
            // Will come back and put in transformation
            status = XLALSimInspiralChooseTDWaveform(hplus, hcross, phiRef,
                    deltaT, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref,
                    r, i, lambda1, lambda2, waveFlags, nonGRparams, amplitudeO,
                    phaseO, approximant);
            if (status == XLAL_FAILURE) return status;

            // FIXME: Need to add hlms, dynamic variables, etc. in cache
            return StoreTDHCache(cache, *hplus, *hcross, phiRef, deltaT, m1, m2,
                    S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref, r, i,
                    lambda1, lambda2, waveFlags, nonGRparams, amplitudeO,
                    phaseO, approximant);

        }
        if( changedParams & DISTANCE ) {
            // Return rescaled copy of cached polarizations
            dist_ratio = cache->r / r;
            *hplus = XLALCreateREAL8TimeSeries(cache->hplus->name,
                    &(cache->hplus->epoch), cache->hplus->f0,
                    cache->hplus->deltaT, &(cache->hplus->sampleUnits),
                    cache->hplus->data->length);
            if (*hplus == NULL) return XLAL_ENOMEM;

            *hcross = XLALCreateREAL8TimeSeries(cache->hcross->name,
                    &(cache->hcross->epoch), cache->hcross->f0,
                    cache->hcross->deltaT, &(cache->hcross->sampleUnits),
                    cache->hcross->data->length);
            if (*hcross == NULL) {
                XLALDestroyREAL8TimeSeries(*hplus);
                *hplus = NULL;
                return XLAL_ENOMEM;
            }

            for (j = 0; j < cache->hplus->data->length; j++) {
                (*hplus)->data->data[j] = cache->hplus->data->data[j]
                        * dist_ratio;
                (*hcross)->data->data[j] = cache->hcross->data->data[j]
                        * dist_ratio;
            }
        }

        return XLAL_SUCCESS;
    }

    // Catch-all. Unsure what to do, don't try to cache.
    // Basically, you requested a waveform type which is not setup for caching
    // b/c of lack of interest or it's unclear what/how to cache for that model
    else {
        return XLALSimInspiralChooseTDWaveform(hplus, hcross, phiRef, deltaT,
                m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref, r, i,
                lambda1, lambda2, waveFlags, nonGRparams, amplitudeO, phaseO,
                approximant);
    }
}

/**
 * Chooses between different approximants when requesting a waveform to be generated
 * Returns the waveform in the frequency domain.
 * The parameters passed must be in SI units.
 *
 * This version allows caching of waveforms. The most recently generated
 * waveform and its parameters are stored. If the next call requests a waveform
 * that can be obtained by a simple transformation, then it is done.
 * This bypasses the waveform generation and speeds up the code.
 */
int XLALSimInspiralChooseFDWaveformFromCache(
        COMPLEX16FrequencySeries **hptilde,     /**< +-polarization waveform */
        COMPLEX16FrequencySeries **hctilde,     /**< x-polarization waveform */
        REAL8 phiRef,                           /**< reference orbital phase (rad) */
        REAL8 deltaF,                           /**< sampling interval (Hz) */
        REAL8 m1,                               /**< mass of companion 1 (kg) */
        REAL8 m2,                               /**< mass of companion 2 (kg) */
        REAL8 S1x,                              /**< x-component of the dimensionless spin of object 1 */
        REAL8 S1y,                              /**< y-component of the dimensionless spin of object 1 */
        REAL8 S1z,                              /**< z-component of the dimensionless spin of object 1 */
        REAL8 S2x,                              /**< x-component of the dimensionless spin of object 2 */
        REAL8 S2y,                              /**< y-component of the dimensionless spin of object 2 */
        REAL8 S2z,                              /**< z-component of the dimensionless spin of object 2 */
        REAL8 f_min,                            /**< starting GW frequency (Hz) */
        REAL8 f_max,                            /**< ending GW frequency (Hz) */
        REAL8 r,                                /**< distance of source (m) */
        REAL8 i,                                /**< inclination of source (rad) */
        REAL8 lambda1,                          /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
        REAL8 lambda2,                          /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
        LALSimInspiralWaveformFlags *waveFlags, /**< Set of flags to control special behavior of some waveform families. Pass in NULL (or None in python) for default flags */
        LALSimInspiralTestGRParam *nonGRparams, /**< Linked list of non-GR parameters. Pass in NULL (or None in python) for standard GR waveforms */
        int amplitudeO,                         /**< twice post-Newtonian amplitude order */
        int phaseO,                             /**< twice post-Newtonian order */
        Approximant approximant,                /**< post-Newtonian approximant to use for waveform production */
        LALSimInspiralWaveformCache *cache      /**< waveform cache structure */
        )
{
    int status;
    size_t j;
    REAL8 dist_ratio, incl_ratio_plus, incl_ratio_cross, phase_diff;
    COMPLEX16 exp_dphi;
    CacheVariableDiffersBitmask changedParams;

    // If nonGRparams are not NULL, don't even try to cache.
    if ( nonGRparams != NULL || (!cache) )
        return XLALSimInspiralChooseFDWaveform(hptilde, hctilde, phiRef, deltaF,
                m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_max, r, i,
                lambda1, lambda2, waveFlags, nonGRparams, amplitudeO, phaseO,
                approximant);

    // Check which parameters have changed
    changedParams = CacheArgsDifferenceBitmask(cache, phiRef, deltaF,
            m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_max, r, i,
            lambda1, lambda2, waveFlags, nonGRparams, amplitudeO,
            phaseO, approximant);

    // No parameters have changed! Copy the cached polarizations
    if( changedParams == NO_DIFFERENCE ) {
        *hptilde = XLALCutCOMPLEX16FrequencySeries(cache->hptilde, 0,
                cache->hptilde->data->length);
        if (*hptilde == NULL) return XLAL_ENOMEM;
        *hctilde = XLALCutCOMPLEX16FrequencySeries(cache->hctilde, 0,
                cache->hctilde->data->length);
        if (*hctilde == NULL) {
            XLALDestroyCOMPLEX16FrequencySeries(*hptilde);
            *hptilde = NULL;
            return XLAL_ENOMEM;
        }

        return XLAL_SUCCESS;
    }

    // Intrinsic parameters have changed. We must generate a new waveform
    if( (changedParams & INTRINSIC) != 0 ) {
        status = XLALSimInspiralChooseFDWaveform(hptilde, hctilde, phiRef,
                deltaF, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_max,
                r, i, lambda1, lambda2, waveFlags, nonGRparams, amplitudeO,
                phaseO, approximant);
        if (status == XLAL_FAILURE) return status;

        return StoreFDHCache(cache, *hptilde, *hctilde, phiRef, deltaF, m1, m2,
            S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_max, r, i, lambda1, lambda2,
            waveFlags, nonGRparams, amplitudeO, phaseO, approximant);
    }

    // case 1: Non-precessing, 2nd harmonic only
    if( approximant == TaylorF2 || approximant == TaylorF2RedSpin
                || approximant == TaylorF2RedSpinTidal
                || approximant == IMRPhenomA || approximant == IMRPhenomB
                || approximant == IMRPhenomC ) {
        // If polarizations are not cached we must generate a fresh waveform
        // FIXME: Will need to check hlms and/or dynamical variables as well
        if( cache->hptilde == NULL || cache->hctilde == NULL) {
            status = XLALSimInspiralChooseFDWaveform(hptilde, hctilde, phiRef,
                    deltaF, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_max,
                    r, i, lambda1, lambda2, waveFlags, nonGRparams, amplitudeO,
                    phaseO, approximant);
            if (status == XLAL_FAILURE) return status;

            return StoreFDHCache(cache, *hptilde, *hctilde, phiRef, deltaF,
                    m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_max, r, i,
                    lambda1, lambda2, waveFlags, nonGRparams, amplitudeO,
                    phaseO, approximant);
        }

        // Set transformation coefficients for identity transformation.
        // We'll adjust them depending on which extrinsic parameters changed.
        dist_ratio = incl_ratio_plus = incl_ratio_cross = 1.;
        phase_diff = 0.;
        exp_dphi = 1.;

        if( changedParams & PHI_REF ) {
            // Only 2nd harmonic present, so {h+,hx} \propto e^(2 i phiRef)
            phase_diff = 2.*(phiRef - cache->phiRef);
            exp_dphi = cpolar(1., phase_diff);
        }
        if( changedParams & INCLINATION) {
            // Rescale h+, hx by ratio of new/old inclination dependence
            incl_ratio_plus = (1.0 + cos(i)*cos(i))
                    / (1.0 + cos(cache->i)*cos(cache->i));
            incl_ratio_cross = cos(i) / cos(cache->i);
        }
        if( changedParams & DISTANCE ) {
            // Rescale h+, hx by ratio of (1/new_dist)/(1/old_dist) = old/new
            dist_ratio = cache->r / r;
        }

        // Create the output polarizations
        *hptilde = XLALCreateCOMPLEX16FrequencySeries(cache->hptilde->name,
                &(cache->hptilde->epoch), cache->hptilde->f0,
                cache->hptilde->deltaF, &(cache->hptilde->sampleUnits),
                cache->hptilde->data->length);
        if (*hptilde == NULL) return XLAL_ENOMEM;

        *hctilde = XLALCreateCOMPLEX16FrequencySeries(cache->hctilde->name,
                &(cache->hctilde->epoch), cache->hctilde->f0,
                cache->hctilde->deltaF, &(cache->hctilde->sampleUnits),
                cache->hctilde->data->length);
        if (*hctilde == NULL) {
            XLALDestroyCOMPLEX16FrequencySeries(*hptilde);
            *hptilde = NULL;
            return XLAL_ENOMEM;
        }

        // Get new polarizations by transforming the old
        incl_ratio_plus *= dist_ratio;
        incl_ratio_cross *= dist_ratio;
        for (j = 0; j < cache->hptilde->data->length; j++) {
            (*hptilde)->data->data[j] = exp_dphi * incl_ratio_plus
                    * cache->hptilde->data->data[j];
            (*hctilde)->data->data[j] = exp_dphi * incl_ratio_cross
                    * cache->hctilde->data->data[j];
        }

        return XLAL_SUCCESS;
    }

    // case 2: Precessing
    /*else if( approximant == SpinTaylorF2 ) {

    }*/

    // Catch-all. Unsure what to do, don't try to cache.
    // Basically, you requested a waveform type which is not setup for caching
    // b/c of lack of interest or it's unclear what/how to cache for that model
    else {
        return XLALSimInspiralChooseFDWaveform(hptilde, hctilde, phiRef, deltaF,
                m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_max, r, i,
                lambda1, lambda2, waveFlags, nonGRparams, amplitudeO, phaseO,
                approximant);
    }

}

/**
 * Construct and initialize a waveform cache.  Caches are used to
 * avoid re-computation of waveforms that differ only by simple
 * scaling relations in extrinsic parameters.
 */
LALSimInspiralWaveformCache *XLALCreateSimInspiralWaveformCache()
{
    LALSimInspiralWaveformCache *cache = XLALCalloc(1,
            sizeof(LALSimInspiralWaveformCache));

    return cache;
}

/**
 * Destroy a waveform cache.
 */
void XLALDestroySimInspiralWaveformCache(LALSimInspiralWaveformCache *cache)
{
    if (cache != NULL) {
        XLALDestroyREAL8TimeSeries(cache->hplus);
        XLALDestroyREAL8TimeSeries(cache->hcross);
        XLALDestroyCOMPLEX16FrequencySeries(cache->hptilde);
        XLALDestroyCOMPLEX16FrequencySeries(cache->hctilde);

        XLALFree(cache);
    }
}

/**
 * Function to compare the requested arguments to those stored in the cache,
 * returns a bitmask which determines if a cached waveform can be recycled.
 */
static CacheVariableDiffersBitmask CacheArgsDifferenceBitmask(
        LALSimInspiralWaveformCache *cache,
        REAL8 phiRef,
        REAL8 deltaTF,
        REAL8 m1,
        REAL8 m2,
        REAL8 S1x, REAL8 S1y, REAL8 S1z,
        REAL8 S2x, REAL8 S2y, REAL8 S2z,
        REAL8 f_min, REAL8 f_ref_max,
        REAL8 r,
        REAL8 i,
        REAL8 lambda1,
        REAL8 lambda2,
        LALSimInspiralWaveformFlags *waveFlags,
        LALSimInspiralTestGRParam *nonGRparams,
        int amplitudeO,
        int phaseO,
        Approximant approximant
        )
{
    CacheVariableDiffersBitmask difference = NO_DIFFERENCE;

    if (cache == NULL) return INTRINSIC;

    if ( !XLALSimInspiralWaveformFlagsEqual(waveFlags, cache->waveFlags) )
        return INTRINSIC;

    if ( deltaTF != cache->deltaTF) return INTRINSIC;
    if ( m1 != cache->m1) return INTRINSIC;
    if ( m2 != cache->m2) return INTRINSIC;
    if ( S1x != cache->S1x) return INTRINSIC;
    if ( S1y != cache->S1y) return INTRINSIC;
    if ( S1z != cache->S1z) return INTRINSIC;
    if ( S2x != cache->S2x) return INTRINSIC;
    if ( S2y != cache->S2y) return INTRINSIC;
    if ( S2z != cache->S2z) return INTRINSIC;
    if ( f_min != cache->f_min) return INTRINSIC;
    if ( f_ref_max != cache->f_ref_max) return INTRINSIC;
    if ( lambda1 != cache->lambda1) return INTRINSIC;
    if ( lambda2 != cache->lambda2) return INTRINSIC;
    if ( nonGRparams != cache->nonGRparams) return INTRINSIC;
    if ( amplitudeO != cache->amplitudeO) return INTRINSIC;
    if ( phaseO != cache->phaseO) return INTRINSIC;
    if ( approximant != cache->approximant) return INTRINSIC;

    if (r != cache->r) difference = difference | DISTANCE;
    if (phiRef != cache->phiRef) difference = difference | PHI_REF;
    if (i != cache->i) difference = difference | INCLINATION;

    return difference;
}

/** Store the output TD hplus and hcross in the cache. */
static int StoreTDHCache(LALSimInspiralWaveformCache *cache,
        REAL8TimeSeries *hplus,
        REAL8TimeSeries *hcross,
        REAL8 phiRef,
        REAL8 deltaT,
        REAL8 m1, REAL8 m2,
        REAL8 S1x, REAL8 S1y, REAL8 S1z,
        REAL8 S2x, REAL8 S2y, REAL8 S2z,
        REAL8 f_min, REAL8 f_ref_max,
        REAL8 r,
        REAL8 i,
        REAL8 lambda1, REAL8 lambda2,
        LALSimInspiralWaveformFlags *waveFlags,
        LALSimInspiralTestGRParam *nonGRparams,
        int amplitudeO,
        int phaseO,
        Approximant approximant
        )
{
    /* Clear any frequency-domain data. */
    if (cache->hptilde != NULL) {
        XLALDestroyCOMPLEX16FrequencySeries(cache->hptilde);
        cache->hptilde = NULL;
    }

    if (cache->hctilde != NULL) {
        XLALDestroyCOMPLEX16FrequencySeries(cache->hctilde);
        cache->hctilde = NULL;
    }

    /* Store params in cache */
    cache->phiRef = phiRef;
    cache->deltaTF = deltaT;
    cache->m1 = m1;
    cache->m2 = m2;
    cache->S1x = S1x;
    cache->S1y = S1y;
    cache->S1z = S1z;
    cache->S2x = S2x;
    cache->S2y = S2y;
    cache->S2z = S2z;
    cache->f_min = f_min;
    cache->f_ref_max = f_ref_max;
    cache->r = r;
    cache->i = i;
    cache->lambda1 = lambda1;
    cache->lambda2 = lambda2;
    cache->waveFlags = waveFlags;
    cache->nonGRparams = nonGRparams;
    cache->amplitudeO = amplitudeO;
    cache->phaseO = phaseO;
    cache->approximant = approximant;

    // Copy over the waveforms
    // NB: XLALCut... creates a new Series object and copies data and metadata
    XLALDestroyREAL8TimeSeries(cache->hplus);
    XLALDestroyREAL8TimeSeries(cache->hcross);
    cache->hplus = XLALCutREAL8TimeSeries(hplus, 0, hplus->data->length);
    if (cache->hplus == NULL) return XLAL_ENOMEM;
    cache->hcross = XLALCutREAL8TimeSeries(hcross, 0, hcross->data->length);
    if (cache->hcross == NULL) {
        XLALDestroyREAL8TimeSeries(cache->hplus);
        cache->hplus = NULL;
        return XLAL_ENOMEM;
    }

    return XLAL_SUCCESS;
}

/** Store the output FD hptilde and hctilde in cache. */
static int StoreFDHCache(LALSimInspiralWaveformCache *cache,
        COMPLEX16FrequencySeries *hptilde,
        COMPLEX16FrequencySeries *hctilde,
        REAL8 phiRef,
        REAL8 deltaT,
        REAL8 m1, REAL8 m2,
        REAL8 S1x, REAL8 S1y, REAL8 S1z,
        REAL8 S2x, REAL8 S2y, REAL8 S2z,
        REAL8 f_min, REAL8 f_ref_max,
        REAL8 r,
        REAL8 i,
        REAL8 lambda1, REAL8 lambda2,
        LALSimInspiralWaveformFlags *waveFlags,
        LALSimInspiralTestGRParam *nonGRparams,
        int amplitudeO,
        int phaseO,
        Approximant approximant
        )
{
    /* Clear any time-domain data. */
    if (cache->hplus != NULL) {
        XLALDestroyREAL8TimeSeries(cache->hplus);
        cache->hplus = NULL;
    }

    if (cache->hcross != NULL) {
        XLALDestroyREAL8TimeSeries(cache->hcross);
        cache->hcross = NULL;
    }

    /* Store params in cache */
    cache->phiRef = phiRef;
    cache->deltaTF = deltaT;
    cache->m1 = m1;
    cache->m2 = m2;
    cache->S1x = S1x;
    cache->S1y = S1y;
    cache->S1z = S1z;
    cache->S2x = S2x;
    cache->S2y = S2y;
    cache->S2z = S2z;
    cache->f_min = f_min;
    cache->f_ref_max = f_ref_max;
    cache->r = r;
    cache->i = i;
    cache->lambda1 = lambda1;
    cache->lambda2 = lambda2;
    cache->waveFlags = waveFlags;
    cache->nonGRparams = nonGRparams;
    cache->amplitudeO = amplitudeO;
    cache->phaseO = phaseO;
    cache->approximant = approximant;

    // Copy over the waveforms
    // NB: XLALCut... creates a new Series object and copies data and metadata
    XLALDestroyCOMPLEX16FrequencySeries(cache->hptilde);
    XLALDestroyCOMPLEX16FrequencySeries(cache->hctilde);
    cache->hptilde = XLALCutCOMPLEX16FrequencySeries(hptilde, 0,
            hptilde->data->length);
    if (cache->hptilde == NULL) return XLAL_ENOMEM;
    cache->hctilde = XLALCutCOMPLEX16FrequencySeries(hctilde, 0,
            hctilde->data->length);
    if (cache->hctilde == NULL) {
        XLALDestroyCOMPLEX16FrequencySeries(cache->hptilde);
        cache->hptilde = NULL;
        return XLAL_ENOMEM;
    }

    return XLAL_SUCCESS;
}
