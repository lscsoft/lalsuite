/*
 * Copyright (C) 2019 Cecilio García Quirós, Sascha Husa
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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */
 /* This code applies the multibanding technique described in arXiv:2001.10897 to the model IMRPhenomXHM described in arXiv:2001.10914. */

#ifndef _LALSIM_IMR_PHENOMXHM_MULTIBAND_H
#define _LALSIM_IMR_PHENOMXHM_MULTIBAND_H

#ifdef __cplusplus
extern "C" {
  #endif
  #ifdef __GNUC__
  #define UNUSED __attribute__ ((unused))
  #else
  #define UNUSED
  #endif

#include "LALSimIMRPhenomX_internals.h"

/*** Core functions to call the waveform with multibanding ***/
int IMRPhenomXHMMultiBandOneMode(
    COMPLEX16FrequencySeries **htildelm, /**< [out] FD waveform **/
    IMRPhenomXWaveformStruct *pWF,       /**< Structure of 22 mode **/
    UINT4 ell,                           /**< First index (l,m) mode **/
    UINT4 emm,                           /**< Second incex (l,m) mode **/
    LALDict *lalParams                   /**< LAL dictionary **/
);

int IMRPhenomXHMMultiBandOneModeMixing(
    COMPLEX16FrequencySeries **htildelm, /**< [out] FD waveform */
    COMPLEX16FrequencySeries *htilde22,  /**< Recycle the 22 mode if previously computed **/
    IMRPhenomXWaveformStruct *pWF,       /**< Structure of 22 mode **/
    UINT4 ell,                           /**< First index (l,m) mode **/
    UINT4 emm,                           /**< Second incex (l,m) mode **/
    LALDict *lalParams                   /**< LAL dictionary **/
);

/*** Set up frequency array and initialize amplitude/phase frequency series ***/
static int SetupWFArraysReal(
  REAL8Sequence **freqs,           /**<[out] Frequency array to evaluate model **/
  REAL8FrequencySeries **amphase,  /**<[out] Initialize amplitude or phase with the length of freqs **/
  REAL8Sequence *freqs_In,         /**< Input frequency array or fmin, fmax **/
  IMRPhenomXWaveformStruct *pWF,   /**< Structure of the 22 mode **/
  LIGOTimeGPS ligotimegps_zero     /**< Needed to initialize amphase **/
);

/*** Functions to compute coarse amplitude and phase ***/

static int IMRPhenomXHM_Amplitude(
  REAL8FrequencySeries **amplm,           /**<[out] amplitude of hlm mode **/
  REAL8Sequence *freqs_In,                /**< Frequency array to evaluate model or fmin, fmax  **/
  IMRPhenomXWaveformStruct *pWF,          /**< Structure of the 22 mode **/
  IMRPhenomXAmpCoefficients *pAmp22,      /**< Amplitude coefficients 22 */
  IMRPhenomXPhaseCoefficients *pPhase22,  /**< Phase coefficients 22 */
  IMRPhenomXHMWaveformStruct *pWFHM,      /**< waveform parameters lm mode */
  IMRPhenomXHMAmpCoefficients *pAmp,      /**< Amplitude coefficients lm */
  IMRPhenomXHMPhaseCoefficients *pPhase   /**< Phase coefficients 22 */
);

static int IMRPhenomXHM_AmplitudeMixing(
  REAL8FrequencySeries **amplm,           /**<[out] amplitude of hlm mode **/
  REAL8Sequence *freqs_In,                /**< Frequency array to evaluate model or fmin, fmax  **/
  IMRPhenomXWaveformStruct *pWF,          /**< Structure of the 22 mode **/
  IMRPhenomXHMWaveformStruct *pWFHM,      /**< waveform parameters lm mode */
  IMRPhenomXHMAmpCoefficients *pAmp,      /**< Amplitude coefficients lm */
  IMRPhenomXHMPhaseCoefficients *pPhase   /**< Phase coefficients 22 */
);

static int IMRPhenomXHM_Phase(
  REAL8FrequencySeries **phaselm,         /**<[out] phase of hlm mode **/
  REAL8Sequence *freqs_In,                /**< Frequency array to evaluate model or fmin, fmax  **/
  IMRPhenomXWaveformStruct *pWF,          /**< Structure of the 22 mode **/
  IMRPhenomXAmpCoefficients *pAmp22,      /**< Amplitude coefficients 22 */
  IMRPhenomXPhaseCoefficients *pPhase22,  /**< Phase coefficients 22 */
  IMRPhenomXHMWaveformStruct *pWFHM,      /**< waveform parameters lm mode */
  IMRPhenomXHMAmpCoefficients *pAmp,      /**< Amplitude coefficients lm */
  IMRPhenomXHMPhaseCoefficients *pPhase   /**< Phase coefficients 22 */
);

static int IMRPhenomXHM_PhaseMixing(
  REAL8FrequencySeries **phaselm,         /**<[out] phase of hlm mode **/
  REAL8Sequence *freqs_In,                /**< Frequency array to evaluate model or fmin, fmax  **/
  IMRPhenomXWaveformStruct *pWF,          /**< Structure of the 22 mode **/
  IMRPhenomXHMWaveformStruct *pWFHM,      /**< waveform parameters lm mode */
  IMRPhenomXHMPhaseCoefficients *pPhase   /**< Phase coefficients 22 */
);

/* Log function */
static double logbase(double base, double x);


/**** Multibanding grids ****/

typedef struct tagIMRPhenomXMultiBandingGridStruct
{
        /* Debug flag */
        INT4 debug;

        /* Model Version Parameters */
        INT4  nIntervals;           //
        INT4  Length;               //
        INT4  intdfRatio;

        REAL8 xStart;               //
        REAL8 xEndRequested;        //
        REAL8 xEndFrom_xStart_dx;   //
        REAL8 xMax;                 //
        REAL8 deltax;               //

}IMRPhenomXMultiBandingGridStruct;

/* Equally spaced grid */
IMRPhenomXMultiBandingGridStruct XLALSimIMRPhenomXGridComp(
  REAL8 fSTART,  /**< Starting frequency of the uniform bin **/
  REAL8 fEND,    /**< Ending frequency of the uniform bin **/
  REAL8 mydf     /**< Frequency spacing of the bin **/
);

/* Non-uniform spaced grid */
INT4 XLALSimIMRPhenomXMultibandingGrid(
  REAL8 fstartIn,                             /**< Minimun frequency in NR unit s**/
  REAL8 fend,                                 /**< End of inspiral frequency bins **/
  REAL8 MfLorentzianEnd,                      /**< Determines the last frequency bin **/
  REAL8 Mfmax,                                /**< Maximun frequency in NR units **/
  REAL8 evaldMf,                              /**< Spacing of the uniform frequency grid (NR units) **/
  REAL8 dfpower,                              /**< decaying frequency power to estimate frequency spacing **/
  REAL8 dfcoefficient,                        /**< multiplying factor to the estimate of the frequency spacing **/
  IMRPhenomXMultiBandingGridStruct *allGrids, /**<[out] list of non-uniform frequency bins**/
  REAL8 dfmerger,                             /**<[out] Spacing merger bin**/
  REAL8 dfringdown                            /**<[out] Spacing ringdown bin**/
);

/**** Interpolating functions ****/
static int interpolateAmplitude(
  double *fineAmp,      /**<[out] amplitude in the fine uniform grid **/
  double coarsefreqs[], /**< non-uniform frequency array**/
  double coarseAmp[],   /**< amplitude in the non-uniform frequency array **/
  double finefreqs[],   /**< uniform fine frequency grid**/
  int lengthCoarse,     /**< length of non-uniform freq array **/
  int lengthFine,       /**< length of uniform fine freq array **/
  int ampinterpolorder     /**< order of the gsl interpolation **/
);

static int interpolateAmplitudeMixing(
  double *fineAmp,      /**<[out] amplitude in the fine uniform grid **/
  double *fineAmpSS,      /**<[out] spheroidal amplitude in the fine uniform grid **/
  double coarsefreqs[], /**< non-uniform frequency array**/
  double coarsefreqsSS[], /**< non-uniform frequency array spheroidal **/
  double coarseAmp[],   /**< amplitude in the non-uniform frequency array **/
  double coarseAmpSS[],  /**< spheroidal amplitude in the non-uniform frequency array **/
  double finefreqs[],   /**< uniform fine frequency grid**/
  int lengthCoarse,     /**< length of non-uniform freq array **/
  int lengthCoarseSS,     /**< length of non-uniform freq array Sphroidal **/
  int lengthFine,        /**< length of uniform fine freq array **/
  int sphericalfinecount,     /**< length of spherical fine grid **/
  int sphericalfinecountMax,     /**< length of spherical fine grid **/
  int ampinterpolorder  /**< order of interpolation **/
);

static double deltaF_mergerBin(REAL8 fdamp, REAL8 alpha4, REAL8 abserror);
static double deltaF_ringdownBin(REAL8 fdamp, REAL8 alpha4, REAL8 LAMBDA, REAL8 abserror);

  #ifdef __cplusplus
}
#endif

#endif /* LALSimIMRPhenomXHM_multiband_h */
