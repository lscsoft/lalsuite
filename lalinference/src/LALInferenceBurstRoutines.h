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

#ifndef _LALINFERENCEBURSTROUTINES_H
#define _LALINFERENCEBURSTROUTINES_H

#include <gsl/gsl_rng.h>
#include <lal/LALDatatypes.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif



/** 
 * XLAL function to determine burst approximant from a string.  The string need not 
 * match exactly, only contain a member of the BurstApproximant enum.
 */
int XLALGetBurstApproximantFromString(const CHAR *inString);
int XLALCheckBurstApproximantFromString(const CHAR *inString);


/** Enum that specifies the PN approximant to be used in computing the waveform.
*/
typedef enum {
   SineGaussianFFast,
   SineGaussianF,
   SineGaussian,
   GaussianF,
   Gaussian,
   RingdownF,
   DampedSinusoidF,
   DampedSinusoid,
   NumBurstApproximants	/**< Number of elements in enum, useful for checking bounds */
 } BurstApproximant;

char* XLALGetStringFromBurstApproximant(BurstApproximant approximant);

int XLALSimBurstImplementedTDApproximants( 
BurstApproximant approximant /**< Burst approximant (see enum in LALInferenceBurst.h) */
    );
int XLALSimBurstImplementedFDApproximants( 
BurstApproximant approximant /**< Burst approximant (see enum in LALInferenceBurst.h) */
    );    
    /** Enumeration to specify time or frequency domain */

int XLALInferenceBurstSineGaussian(
	REAL8TimeSeries **hplus,
	REAL8TimeSeries **hcross,
	REAL8 Q,
	REAL8 centre_frequency,
	REAL8 hrss,
	REAL8 eccentricity,
	REAL8 phase,
	REAL8 delta_t // 1 over srate
);

int XLALInferenceBurstSineGaussianF(
	COMPLEX16FrequencySeries **hplus,
	COMPLEX16FrequencySeries **hcross,
	REAL8 Q,
	REAL8 centre_frequency,
	REAL8 hrss,
	REAL8 eccentricity,
        REAL8 phase,
	REAL8 deltaF,
    REAL8 deltaT
);

int XLALInferenceBurstSineGaussianFFast(
  COMPLEX16FrequencySeries **hplus,
  COMPLEX16FrequencySeries **hcross,
  REAL8 Q,
  REAL8 centre_frequency,
  REAL8 hrss,
  REAL8 eccentricity,
  REAL8 phase,  
  REAL8 deltaF,
  REAL8 deltaT
);

int XLALInferenceBurstGaussianF(
	COMPLEX16FrequencySeries **hplus,
	COMPLEX16FrequencySeries **hcross,
	REAL8 duration,
	REAL8 hrss,
	REAL8 alpha,
	REAL8 deltaF,
  REAL8 deltaT
);

int XLALInferenceBurstDampedSinusoid(
	REAL8TimeSeries **hplus,
	REAL8TimeSeries **hcross,
	REAL8 Q,
	REAL8 centre_frequency,
	REAL8 hrss,
	REAL8 eccentricity,
	REAL8 phase,
	REAL8 delta_t // 1 over srate
);

int XLALInferenceBurstDampedSinusoidF(
  COMPLEX16FrequencySeries **hplus,
  COMPLEX16FrequencySeries **hcross,
  REAL8 Q,
  REAL8 centre_frequency,
  REAL8 hrss,
  REAL8 eccentricity,
  REAL8 phase,  
  REAL8 deltaF,
  REAL8 deltaT
);

int XLALInferenceBurstGaussian(
	REAL8TimeSeries **hplus,
	REAL8TimeSeries **hcross,
	REAL8 duration,
	REAL8 hrss,
	REAL8 eccentricity,
	REAL8 polarization,
	REAL8 delta_t
);

/* ======================= Extra params ============================= */


/**
 * Linked list node for the testing GR parameters
 */
typedef struct tagLALSimBurstExtraParamData
{
    char name[32]; 	/**< Name of the variable */
    double value;  	/**< Value of the variable */
} LALSimBurstExtraParamData;

/**
 * Linked list of any number of parameters for testing GR
 */
typedef struct tagLALSimBurstExtraParam
{
    struct tagLALSimBurstExtraParamData *data; /**< Current variable */
    struct tagLALSimBurstExtraParam *next; /**< The next variable in linked list */
}  LALSimBurstExtraParam;

/**
 * Function that creates the head node of the extra burst parameters linked list.
 * It is initialized with a single parameter with given name and value
 */
 
#ifdef SWIG   // SWIG interface directives
SWIGLAL(INOUT_STRUCTS(LALSimBurstExtraParam**, parameter));
#endif 

LALSimBurstExtraParam *XLALSimBurstCreateExtraParam(
        const char *name, /**< Name of first parameter in new linked list */
        double value 	 /**< Value of first parameter in new linked list */
        );

/**
 * Function that adds a parameter to the extra burst parameters linked list. If the
 * parameter already exists, it throws an error.
 */
int XLALSimBurstAddExtraParam(
        LALSimBurstExtraParam **parameter, /**< Pointer to the head node of the linked list of parameters */
        const char *name, 			/**< Parameter name */
        const double value 			/**< Parameter value */
        );

/**
 * Function that sets the value of the desired parameter in the burst extra
 * parameters linked list to 'value'.  Throws an error if the parameter
 * is not found
 */
int XLALSimBurstSetExtraParam(
        LALSimBurstExtraParam *parameter, /**< Linked list to be modified */
        const char *name, 		/**< Name of parameter to be modified */
        const double value 		/**< New value for parameter */
        );

/**
 * Function that returns the value of the desired parameters in the
 * burst extra parameters linked list.  Aborts if the parameter is not found
 */
double XLALSimBurstGetExtraParam(
        const LALSimBurstExtraParam *parameter, /**< Linked list to retrieve from */
        const char *name 	/**< Name of parameter to be retrieved */
        );

/**
 * Function that checks whether the requested parameter exists within the
 * burst extra parameters linked list.  Returns true (1) or false (0) accordingly
 */
int XLALSimBurstExtraParamExists(
        const LALSimBurstExtraParam *parameter, /**< Linked list to check */
        const char *name 		/**< Parameter name to check for */
        );

/** Function that prints the whole burst extra parameters linked list */
int XLALSimBurstPrintExtraParam(
        FILE *fp, 			/** FILE pointer to write to */
        LALSimBurstExtraParam *parameter /**< Linked list to print */
        );

/** Function that destroys the whole extra burst parameters linked list */
void XLALSimBurstDestroyExtraParam(
        LALSimBurstExtraParam *parameter /**< Linked list to destroy */
        );


/* ======================= WF cache params ============================= */


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
        REAL8 q,                                /**< Q (==sqrt(2) \f$\pi\f$ f0 tau ) [dless]*/
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

/* ======================= EVERYTHING ELSE  ===========================*/

int XLALSimBurstChooseFDWaveform(
    COMPLEX16FrequencySeries **hptilde,     /**< FD plus polarization */
    COMPLEX16FrequencySeries **hctilde,     /**< FD cross polarization */
    REAL8 deltaF,                           /**< sampling interval (Hz) */
    REAL8 deltaT,                           /**< time step corresponding to consec */
    REAL8 f0,                               /**< central frequency (Hz) */
    REAL8 q,                                /**< Q (==sqrt(2) \f$\pi\f$ f0 tau ) [dless]*/
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
    REAL8 q,                                /**< Q (==sqrt(2) \f$\pi\f$ f0 tau ) [dless]*/
    REAL8 tau,                              /**< Duration [s] */
    REAL8 f_min,                            /**< starting GW frequency (Hz) */
    REAL8 f_max,                            /**< ending GW frequency (Hz) (0 for Nyquist) */
    REAL8 hrss,                             /**< hrss [strain] */
    REAL8 polar_angle,                      /**< Polar_ellipse_angle as defined in the burst table. Together with polar_ellipse_eccentricity below will fix the ratio of + vs x aplitude. Some WFs uses a single parameter alpha for this. Alpha is passed through extraParams*/
    REAL8 polar_ecc,                        /**< See above */
    LALSimBurstExtraParam *extraParams, /**< Linked list of non-GR parameters. Pass in NULL (or None in python) to neglect these */
    BurstApproximant approximant                 /**< Burst approximant  */
    );

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif
#endif


