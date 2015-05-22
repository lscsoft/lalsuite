/*
 * Copyright (C) 2008 J. Creighton, K. Cannon
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with with program; see the file COPYING. If not, write to the Free
 * Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 */


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */

#ifndef _LALSIMBURST_H
#define _LALSIMBURST_H

#include <gsl/gsl_rng.h>
#include <lal/LALDatatypes.h>
#include <lal/LALSimBurstExtraParams.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/*
 * ============================================================================
 *
 *                            Function Prototypes
 *
 * ============================================================================
 */


int XLALGenerateImpulseBurst(
	REAL8TimeSeries **hplus,
	REAL8TimeSeries **hcross,
	REAL8 hpeak,
	REAL8 delta_t
);


int XLALGenerateBandAndTimeLimitedWhiteNoiseBurst(
	REAL8TimeSeries **hplus,
	REAL8TimeSeries **hcross,
	REAL8 duration,
	REAL8 frequency,
	REAL8 bandwidth,
	REAL8 int_hdot_squared,
	REAL8 delta_t,
	gsl_rng *rng
);


int XLALGenerateStringCusp(
	REAL8TimeSeries **hplus,
	REAL8TimeSeries **hcross,
	REAL8 amplitude,
	REAL8 f_high,
	REAL8 delta_t
);


int XLALSimBurstSineGaussian(
	REAL8TimeSeries **hplus,
	REAL8TimeSeries **hcross,
	REAL8 Q,
	REAL8 centre_frequency,
	REAL8 hrss,
	REAL8 eccentricity,
	REAL8 polarization,
	REAL8 delta_t
);

int XLALSimBurstGaussian(
	REAL8TimeSeries **hplus,
	REAL8TimeSeries **hcross,
	REAL8 duration,
	REAL8 hrss,
	REAL8 eccentricity,
	REAL8 polarization,
	REAL8 delta_t
);

int XLALSimBurstImg(
	REAL8TimeSeries **hplus,
	REAL8TimeSeries **hcross, 
	REAL8Array *image,
	double dt,
	double df,
	double fstart,
	double hrss,
	double deltaT,
	gsl_rng *rng
);


int XLALSimUnicorn(
	REAL8TimeSeries **hplus,
	REAL8TimeSeries **hcross,
	double f_min,
	double f_max,
	double V,
	double hrss,
	double deltaT,
	gsl_rng *rng
);


int XLALSimBurstSineGaussianF(
	COMPLEX16FrequencySeries **hplus,
	COMPLEX16FrequencySeries **hcross,
	REAL8 Q,
	REAL8 centre_frequency,
	REAL8 hrss,
  REAL8 alpha,
	REAL8 deltaF,
    REAL8 deltaT
);

int XLALSimBurstGaussianF(
	COMPLEX16FrequencySeries **hplus,
	COMPLEX16FrequencySeries **hcross,
	REAL8 duration,
	REAL8 hrss,
	REAL8 alpha,
	REAL8 deltaF,
  REAL8 deltaT
);

int XLALSimBurstSineGaussianFFast(
	COMPLEX16FrequencySeries **hplus,
	COMPLEX16FrequencySeries **hcross,
	REAL8 Q,
	REAL8 centre_frequency,
	REAL8 hrss,
	REAL8 alpha,
	REAL8 deltaF,
  REAL8 deltaT
);

int XLALSimBurstDampedSinusoid(
        REAL8TimeSeries **hplus,
        REAL8TimeSeries **hcross,
        REAL8 Q,
        REAL8 centre_frequency,
        REAL8 hrss,
        REAL8 eccentricity,
        REAL8 polarization,
        REAL8 delta_t
);

int XLALSimBurstDampedSinusoidF(
        COMPLEX16FrequencySeries **hplus,
        COMPLEX16FrequencySeries **hcross,
        REAL8 Q,
        REAL8 centre_frequency,
        REAL8 hrss,
        REAL8 alpha,
        REAL8 deltaF,
    REAL8 deltaT
);

REAL8 XLALMeasureHPeak(const REAL8TimeSeries *);
REAL8 XLALMeasureIntS1S2DT(const REAL8TimeSeries *, const REAL8TimeSeries *);
REAL8 XLALMeasureHrss(const REAL8TimeSeries *, const REAL8TimeSeries *);
REAL8 XLALMeasureIntHDotSquaredDT(const COMPLEX16FrequencySeries *);
REAL8 XLALMeasureEoverRsquared(REAL8TimeSeries *, REAL8TimeSeries *);


/** Enum that specifies the PN approximant to be used in computing the waveform.
*/
typedef enum {
   SineGaussianF,
   SineGaussian,
   GaussianF,
   Gaussian,
   RingdownF,
   DampedSinusoidF,
   DampedSinusoid,
   NumBurstApproximants	/**< Number of elements in enum, useful for checking bounds */
 } BurstApproximant;

int XLALSimBurstImplementedTDApproximants( 
BurstApproximant approximant /**< Burst approximant (see enum in LALSimBurst.h) */
    );
int XLALSimBurstImplementedFDApproximants( 
BurstApproximant approximant /**< Burst approximant (see enum in LALSimBurst.h) */
    );    
    /** Enumeration to specify time or frequency domain */
int XLALSimBurstChooseFDWaveform(
    COMPLEX16FrequencySeries **hptilde,     /**< FD plus polarization */
    COMPLEX16FrequencySeries **hctilde,     /**< FD cross polarization */
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
    LALSimBurstExtraParam *extraParams, /**< Linked list of non-GR parameters. Pass in NULL (or None in python) to neglect these */
    BurstApproximant approximant                 /**< Burst approximant  */
    );
    
int XLALSimBurstChooseTDWaveform(
    REAL8TimeSeries **hplus,                    /**< +-polarization waveform */
    REAL8TimeSeries **hcross,                   /**< x-polarization waveform */
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
    );
/** 
 * XLAL function to determine burst approximant from a string.  The string need not 
 * match exactly, only contain a member of the BurstApproximant enum.
 */
int XLALGetBurstApproximantFromString(const CHAR *inString);
char* XLALGetStringFromBurstApproximant(BurstApproximant approximant);
int XLALCheckBurstApproximantFromString(const CHAR *inString);


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif
#endif
