/*
 * Copyright (C) 2013 Evan Ochsner and Will M. Farr, 
 *  2014 Salvatore Vitale
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
#include <LALSimBurstWaveformCache.h>
#include <LALSimBurstExtraParams.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
/**
 * Bitmask enumerating which parameters have changed, to determine
 * if the requested waveform can be transformed from a cached waveform
 * or if it must be generated from scratch.
 */
typedef enum {
    NO_DIFFERENCE = 0,
    INTRINSIC = 1,
    HRSS = 2,
    ALPHA = 4
} CacheVariableDiffersBitmask;

static CacheVariableDiffersBitmask CacheArgsDifferenceBitmask(
        LALSimBurstWaveformCache *cache,
        REAL8 deltaT,
        REAL8 deltaF,
        REAL8 f0,
        REAL8 q,
        REAL8 tau,
        REAL8 f_min,
        REAL8 f_max,
        REAL8 hrss,
        REAL8 polar_angle,
        REAL8 polar_ecc,
        LALSimBurstExtraParam *extraParams,
        BurstApproximant approximant
);

static int StoreTDHCache(
        LALSimBurstWaveformCache *cache,
        REAL8TimeSeries *hplus,
        REAL8TimeSeries *hcross,
        REAL8 deltaT,
        REAL8 f0,
        REAL8 q, REAL8 tau,
        REAL8 f_min, REAL8 f_max,
        REAL8 hrss,
        REAL8 polar_angle,
        REAL8 polar_ecc,
        LALSimBurstExtraParam *extraParams,
        BurstApproximant approximant);

static int StoreFDHCache(
        LALSimBurstWaveformCache *cache,
        COMPLEX16FrequencySeries *hptilde,
        COMPLEX16FrequencySeries *hctilde,
        REAL8 deltaF, 
        REAL8 deltaT,
        REAL8 f0,
        REAL8 q, REAL8 tau,
        REAL8 f_min, REAL8 f_max,
        REAL8 hrss,
        REAL8 polar_angle,
        REAL8 polar_ecc,
        LALSimBurstExtraParam *extraParams,
        BurstApproximant approximant);

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
        LALSimBurstExtraParam *extraParams, /**< Linked list of extra Parameters (includes alpha). Pass in NULL (or None in python) to neglect */
        BurstApproximant approximant,           /**< Burst approximant to use for waveform production */
        LALSimBurstWaveformCache *cache      /**< waveform cache structure; use NULL for no caching */
        )
{
    int status;
    size_t j;
    REAL8 hrss_ratio, alpha_ratio_plus,alpha_ratio_cross;
    CacheVariableDiffersBitmask changedParams;
    REAL8 deltaF=0.0; // UNUSED
    if ( (!cache) )
        return XLALSimBurstChooseTDWaveform(hplus, hcross,deltaT,f0,q,tau,f_min,f_max,hrss,polar_angle,polar_ecc,extraParams,approximant);

    // Check which parameters have changed
    changedParams = CacheArgsDifferenceBitmask(cache, deltaT,deltaF,f0,q, tau,f_min,f_max, hrss,polar_angle, polar_ecc,extraParams,approximant);

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
        status = XLALSimBurstChooseTDWaveform(hplus, hcross,deltaT,f0,q,tau,f_min,f_max,hrss,polar_angle,polar_ecc,extraParams,approximant);
        if (status == XLAL_FAILURE) return status;

        // FIXME: Need to add hlms, dynamic variables, etc. in cache
        return StoreTDHCache(cache, *hplus, *hcross, deltaT,f0,q, tau,f_min,f_max, hrss,polar_angle, polar_ecc,extraParams,approximant);
    }

    // If polarizations are not cached we must generate a fresh waveform
    if( cache->hplus == NULL || cache->hcross == NULL) {
        status = XLALSimBurstChooseTDWaveform(hplus, hcross,deltaT,f0,q,tau,f_min,f_max,hrss,polar_angle,polar_ecc,extraParams,approximant);
        if (status == XLAL_FAILURE) return status;

        // FIXME: Need to add hlms, dynamic variables, etc. in cache
        return StoreTDHCache(cache, *hplus, *hcross, deltaT,f0,q, tau,f_min,f_max, hrss,polar_angle, polar_ecc,extraParams,approximant);
    }

    // Set transformation coefficients for identity transformation.
    // We'll adjust them depending on which extrinsic parameters changed.
    hrss_ratio = alpha_ratio_plus = alpha_ratio_cross=1.0;

    if( changedParams & ALPHA) {
        // Rescale h+, hx by ratio of new/old alpha dependence
        REAL8 alpha=XLALSimBurstGetExtraParam(extraParams,"alpha");
        REAL8 cached_alpha=XLALSimBurstGetExtraParam(cache->extraParams,"alpha");
        alpha_ratio_plus = cos(alpha)/cos(cached_alpha);
        alpha_ratio_cross = sin(alpha)/sin(cached_alpha);
    }
    if( changedParams & HRSS ) {
        // Rescale h+, hx by ratio of new_hrss/old_hrss
        hrss_ratio = hrss/cache->hrss;
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
    alpha_ratio_plus *= hrss_ratio;
    alpha_ratio_cross *= hrss_ratio;
    // FIXME: Do changing phiRef and inclination commute?!?!
    for (j = 0; j < cache->hplus->data->length; j++) {
        (*hplus)->data->data[j] = alpha_ratio_plus
                * (cache->hplus->data->data[j]);
        (*hcross)->data->data[j] = alpha_ratio_cross
                *cache->hcross->data->data[j];
    }

    return XLAL_SUCCESS;

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
        BurstApproximant approximant,                 /**< Burst approximant  */
        LALSimBurstWaveformCache *cache      /**< waveform cache structure; use NULL for no caching */
        )
{
    int status;
    size_t j;
    REAL8 hrss_ratio, alpha_ratio_plus,alpha_ratio_cross;
    COMPLEX16 exp_dphi;
    CacheVariableDiffersBitmask changedParams;
  if ((!cache) ){
     
        return XLALSimBurstChooseFDWaveform(hptilde, hctilde,deltaF,deltaT,f0,q, tau,f_min,f_max, hrss,polar_angle, polar_ecc,extraParams,approximant);
     }

    // Check which parameters have changed
    changedParams = CacheArgsDifferenceBitmask(cache, deltaT,deltaF,f0,q, tau,f_min,f_max, hrss,polar_angle, polar_ecc,extraParams,approximant);

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
        status = XLALSimBurstChooseFDWaveform(hptilde, hctilde, deltaF,deltaT,f0,q, tau,f_min,f_max, hrss,polar_angle, polar_ecc,extraParams,approximant);
        if (status == XLAL_FAILURE) return status;

        return StoreFDHCache(cache, *hptilde, *hctilde, deltaF,deltaT,f0,q, tau,f_min,f_max, hrss,polar_angle, polar_ecc,extraParams,approximant);
    }

    // If polarizations are not cached we must generate a fresh waveform
    // FIXME: Will need to check hlms and/or dynamical variables as well
    if( cache->hptilde == NULL || cache->hctilde == NULL) {
        status = XLALSimBurstChooseFDWaveform(hptilde, hctilde, deltaF,deltaT,f0,q, tau,f_min,f_max, hrss,polar_angle, polar_ecc,extraParams,approximant);
        if (status == XLAL_FAILURE) return status;

        return StoreFDHCache(cache, *hptilde, *hctilde, deltaF,deltaT,f0,q, tau,f_min,f_max, hrss,polar_angle, polar_ecc,extraParams,approximant);
    }

    // Set transformation coefficients for identity transformation.
    // We'll adjust them depending on which extrinsic parameters changed.
    hrss_ratio = alpha_ratio_plus = alpha_ratio_cross = 1.;
    exp_dphi = 1.;


    if( changedParams & ALPHA) {
        // Rescale h+, hx by ratio of new/old alpha dependence
        REAL8 alpha=XLALSimBurstGetExtraParam(extraParams,"alpha");
        REAL8 cached_alpha=XLALSimBurstGetExtraParam(cache->extraParams,"alpha");
        alpha_ratio_plus = cos(alpha)/cos(cached_alpha);
        alpha_ratio_cross = sin(alpha)/sin(cached_alpha);
    }
    if( changedParams & HRSS ) {
        // Rescale h+, hx by ratio of new_hrss/old_hrss
        hrss_ratio = hrss/cache->hrss;
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
    alpha_ratio_plus *= hrss_ratio;
    alpha_ratio_cross *= hrss_ratio;
    for (j = 0; j < cache->hptilde->data->length; j++) {
        (*hptilde)->data->data[j] = exp_dphi * alpha_ratio_plus
                * cache->hptilde->data->data[j];
        (*hctilde)->data->data[j] = exp_dphi * alpha_ratio_cross
                * cache->hctilde->data->data[j];
    }

  return XLAL_SUCCESS;
  
}

/**
 * Construct and initialize a waveform cache.  Caches are used to
 * avoid re-computation of waveforms that differ only by simple
 * scaling relations in extrinsic parameters.
 */
LALSimBurstWaveformCache *XLALCreateSimBurstWaveformCache()
{
    LALSimBurstWaveformCache *cache = XLALCalloc(1,
            sizeof(LALSimBurstWaveformCache));

    return cache;
}

/**
 * Destroy a waveform cache.
 */
void XLALDestroySimBurstWaveformCache(LALSimBurstWaveformCache *cache)
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
    LALSimBurstWaveformCache *cache ,     /**< waveform cache structure; use NULL for no caching */
    REAL8 deltaT,
    REAL8 deltaF,
    REAL8 f0,                               /**< central frequency (Hz) */
    REAL8 q,                                /**< Q (==sqrt(2) \pi f0 tau ) [dless]*/
    REAL8 tau,                              /**< Duration [s] */
    REAL8 f_min,                            /**< starting GW frequency (Hz) */
    REAL8 f_max,                            /**< ending GW frequency (Hz) (0 for Nyquist) */
    REAL8 hrss,                             /**< hrss [strain] */
    REAL8 polar_angle,                      /**< Polar_ellipse_angle as defined in the burst table. Together with polar_ellipse_eccentricity below will fix the ratio of + vs x aplitude. Some WFs uses a single parameter alpha for this. Alpha is passed through extraParams*/
    REAL8 polar_ecc,                        /**< See above */
    LALSimBurstExtraParam *extraParams, /**< Linked list of non-GR parameters. Pass in NULL (or None in python) to neglect these */
    BurstApproximant approximant                 /**< Burst approximant  */
    )
{
    CacheVariableDiffersBitmask difference = NO_DIFFERENCE;

    if (cache == NULL) return INTRINSIC;

    if ( deltaT != cache->deltaT) return INTRINSIC;
    if ( deltaF != cache->deltaF) return INTRINSIC;
    if ( f0 != cache->f0) return INTRINSIC;
    if ( q != cache->q) return INTRINSIC;
    if ( tau != cache->tau) return INTRINSIC;
    if ( f_min != cache->f_min) return INTRINSIC;
    if ( f_max != cache->f_max) return INTRINSIC;
    if ( polar_angle != cache->polar_angle) return INTRINSIC;
    if ( polar_ecc != cache->polar_ecc) return INTRINSIC;
    if ( approximant != cache->approximant) return INTRINSIC;
    if ( polar_angle != cache->polar_angle) return INTRINSIC;
    if ( polar_ecc != cache->polar_ecc) return INTRINSIC;
    
    if (hrss != cache->hrss) difference = difference | HRSS;
    
    if (extraParams!=NULL){
      if (XLALSimBurstExtraParamExists(extraParams,"alpha")){ 
        REAL8 alpha=XLALSimBurstGetExtraParam(extraParams,"alpha");
        REAL8 cached_alpha=1.0;
        if ( XLALSimBurstExtraParamExists(cache->extraParams,"alpha"))
          cached_alpha=XLALSimBurstGetExtraParam(cache->extraParams,"alpha");
        else 
          return INTRINSIC;
        if (alpha!= cached_alpha) difference = difference | ALPHA;
        }
    }
    return difference;
}

/** Store the output TD hplus and hcross in the cache. */
static int StoreTDHCache(
        LALSimBurstWaveformCache *cache,
        REAL8TimeSeries *hplus,
        REAL8TimeSeries *hcross,
        REAL8 deltaT,                           /**< time step corresponding to consec */
        REAL8 f0,                               /**< central frequency (Hz) */
        REAL8 q,                                /**< Q (==sqrt(2) \pi f0 tau ) [dless]*/
        REAL8 tau,                              /**< Duration [s] */
        REAL8 f_min,                            /**< starting GW frequency (Hz) */
        REAL8 f_max,                            /**< ending GW frequency (Hz) (0 for Nyquist) */
        REAL8 hrss,                             /**< hrss [strain] */
        REAL8 polar_angle,                      /**< Polar_ellipse_angle as defined in the burst table. Together with polar_ellipse_eccentricity below will fix the ratio of + vs x aplitude. Some WFs uses a single parameter alpha for this. Alpha is passed through extraParams*/
        REAL8 polar_ecc,                        /**< See above */
        LALSimBurstExtraParam *extraParams, /**< Linked list of non-GR parameters. Pass in NULL (or None in python) to neglect these */
        BurstApproximant approximant                 /**< Burst approximant  */
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
    cache->deltaT = deltaT;
    cache->f0 = f0;
    cache->q = q;
    cache->tau = tau;
    cache->f_min = f_min;
    cache->f_max = f_max;
    cache->hrss = hrss;
    cache->polar_angle =polar_angle;
    cache->polar_ecc = polar_ecc;
    if (extraParams==NULL)
      cache->extraParams=NULL;
    else if (cache->extraParams==NULL){
      /* Initialize to something that won't make the ratio of sin alphas to blow up */
      cache->extraParams=XLALSimBurstCreateExtraParam("alpha",1.0);
    }
    else{
      XLALSimBurstSetExtraParam(cache->extraParams,"alpha",XLALSimBurstGetExtraParam(extraParams,"alpha"));
    }
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
static int StoreFDHCache(LALSimBurstWaveformCache *cache,
        COMPLEX16FrequencySeries *hptilde,
        COMPLEX16FrequencySeries *hctilde,
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
        BurstApproximant approximant                 /**< Burst approximant  */
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
    cache->deltaT = deltaT;
    cache->deltaF = deltaF;
    cache->f0 = f0;
    cache->q = q;
    cache->tau = tau;
    cache->f_min = f_min;
    cache->f_max = f_max;
    cache->hrss = hrss;
    cache->polar_angle =polar_angle;
    cache->polar_ecc = polar_ecc;
    if (extraParams==NULL)
      cache->extraParams=NULL;
    else if (cache->extraParams==NULL){
      /* Initialize to something that won't make the ratio of sin alphas to blow up */
      cache->extraParams=XLALSimBurstCreateExtraParam("alpha",1.0);
    }
    else{
      XLALSimBurstSetExtraParam(cache->extraParams,"alpha",XLALSimBurstGetExtraParam(extraParams,"alpha"));
    }
    
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
