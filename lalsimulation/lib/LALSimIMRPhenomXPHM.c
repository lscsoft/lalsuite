/*
* Copyright (C) 2019 Cecilio García Quirós, Geraint Pratten
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

#include <lal/LALSimIMR.h>
#include <lal/SphericalHarmonics.h>
#include <lal/Sequence.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/LALConstants.h>
#include <lal/FrequencySeries.h>
#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/XLALError.h>

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef _OPENMP
#define omp ignore
#endif

#define L_MAX 4

#ifndef PHENOMXHMDEBUG
#define DEBUG 0
#define PHENOMXDEBUG 0
#else
#define DEBUG 1 //print debugging info
#define PHENOMXDEBUG 1
#endif

#include "LALSimIMRPhenomXPHM.h"


/* Generic routine for twisting up higher multipole models */
static int IMRPhenomXPHMTwistUp(
  const REAL8 Mf,                           /**< Frequency (Hz) */
  const COMPLEX16 hHM,                      /**< Underlying aligned-spin IMRPhenomXAS waveform*/
  IMRPhenomXWaveformStruct *pWF,            /**< IMRPhenomX Waveform Struct */
  IMRPhenomXPrecessionStruct *pPrec,        /**< IMRPhenomXP Precession Struct */
  INT4  ell,                                /**< l index of the (l,m) mode */
  INT4  emmprime,                           /**< m index of the (l,m) mode */
  COMPLEX16 *hp,                            /**< [out] h_+ polarization \f$\tilde h_+\f$ */
  COMPLEX16 *hc                             /**< [out] h_x polarization \f$\tilde h_x\f$ */
);

/*
  Core twisting up routine for one single mode.
  Twist the waveform in the precessing L-frame to the inertial J-frame for one frequency point.
  This function will be inside a loop of frequencies insid a loop over mprime >0 up to l.
*/
static int IMRPhenomXPHMTwistUpOneMode(
  const REAL8 Mf,                          /**< Frequency (Mf geometric units) */
  const COMPLEX16 hlmprime,                /**< Underlying aligned-spin IMRPhenomXHM waveform. The loop is with mprime positive, but the mode has to be the negative one for positive frequencies.*/
  IMRPhenomXWaveformStruct *pWF,           /**< IMRPhenomX Waveform Struct */
  IMRPhenomXPrecessionStruct *pPrec,       /**< IMRPhenomXP Precession Struct */
  UINT4  l,                                /**< l index of the (l,m) (non-)precessing mode */
  UINT4  mprime,                           /**< second index of the (l,mprime) non-precessing mode  */
  INT4   m,                                /**< second index of the (l,m) precessing mode */
  COMPLEX16Sequence *hlminertial           /**< [out] hlm for one frequency in the inertial frame (precessing mode)  */
);


//This is a wrapper function for adding higher modes to the ModeArray
LALDict *IMRPhenomXPHM_setup_mode_array(LALDict *lalParams);


static int IMRPhenomXPHM_hplushcross(
  COMPLEX16FrequencySeries **hptilde,  /**< [out] Frequency domain h+ GW strain */
  COMPLEX16FrequencySeries **hctilde,  /**< [out] Frequency domain hx GW strain */
  REAL8Sequence *freqs_In,             /**< Frequency array to evaluate the model. (fmin, fmax) for equally spaced grids. */
  IMRPhenomXWaveformStruct *pWF,       /**< IMRPhenomX Waveform Struct  */
  IMRPhenomXPrecessionStruct *pPrec,   /**< IMRPhenomXP Precession Struct  */
  LALDict *lalParams                   /**< LAL Dictionary Structure    */
);

static int IMRPhenomXPHM_hplushcross_from_modes(
  COMPLEX16FrequencySeries **hptilde,  /**< [out] Frequency domain h+ GW strain */
  COMPLEX16FrequencySeries **hctilde,  /**< [out] Frequency domain hx GW strain */
  REAL8Sequence *freqs_In,             /**< Frequency array to evaluate the model. (fmin, fmax) for equally spaced grids. */
  IMRPhenomXWaveformStruct *pWF,       /**< IMRPhenomX Waveform Struct  */
  IMRPhenomXPrecessionStruct *pPrec,   /**< IMRPhenomXP Precession Struct  */
  LALDict *lalParams                   /**< LAL Dictionary Structure    */
);

static int IMRPhenomXPHM_OneMode(
  COMPLEX16FrequencySeries **hlmpos,    /**< [out] Frequency domain hlm GW strain inertial frame positive frequencies */
  COMPLEX16FrequencySeries **hlmneg,    /**< [out] Frequency domain hlm GW strain inertial frame negative frequencies */
  REAL8Sequence *freqs_In,              /**< Input frequency grid        */
  IMRPhenomXWaveformStruct *pWF,        /**< IMRPhenomX Waveform Struct  */
  IMRPhenomXPrecessionStruct *pPrec,    /**< IMRPhenomXP Precession Struct  */
  UINT4 ell,                            /**< l index of the (l,m) precessing mode */
  INT4  m,                              /**< m index of the (l,m) precessing mode */
  LALDict *lalParams                    /**< LAL Dictionary Structure    */
);

/* Return the 3 Post-Newtonian Euler angles evaluated in (2*pi*Mf/mprime) */
static int Get_alpha_beta_epsilon(
  REAL8 *alpha,                       /**< [out] Azimuthal angle of L w.r.t J */
  REAL8 *cBetah,                      /**< [out] Cosine of polar angle between L and J */
  REAL8 *sBetah,                      /**< [out] Sine of polar angle between L and J */
  REAL8 *epsilon,                     /**< [out] Minus the third Euler angle (-gamma) describing L w.r.t J, fixed by minimal rotation condition */
  INT4 mprime,                        /**< Second index of the non-precesssing mode (l, mprime) */
  REAL8 Mf,                           /**< Frequency geometric units */
  IMRPhenomXPrecessionStruct *pPrec,  /**< IMRPhenomXP Precessing structure*/
  IMRPhenomXWaveformStruct *pWF       /**< IMRPhenomX Waveform structure*/
);

/* Return the offset at reference frequency for alpha and epsilon Euler angles for a particular non-precessing mode.
  alpha_offset_mprime = alpha(2*pi*MfRef/mprime) - alpha0. Used for Pv2 and Pv3 angles. */
static double Get_alpha_epsilon_offset(
  REAL8 *alpha_offset_mprime,          /**< [out] Offset alpha angle at reference frequency */
  REAL8 *epsilon_offset_mprime,        /**< [out] Offset epsilon angle at reference frequency */
  INT4 mprime,                         /**< Second index of the non-precesssing mode (l, mprime) */
  IMRPhenomXPrecessionStruct *pPrec    /**< IMRPhenomXP Precessing structure*/
);



/**
 * @addtogroup LALSimIMRPhenomX_c
 * @{
 *
 * @name Routines for IMRPhenomXPHM
 * @{
 *
 */


/*********************************************/
/*                                           */
/*      MULTIMODE PRECESSING FUNCTIONS       */
/*                                           */
/*********************************************/

/** Returns hptilde and hctilde of the multimode precessing waveform for positive frequencies in an equally spaced grid.

This is the default function used when calling ChooseFDWaveform.

It computes the non-precessing modes just once and do the twisting up according to eqs. 3.5-3.7 in the Precessing paper.

It is just a wrapper of the internal function that actually carries out the calculation: IMRPhenomXPHM_hplushcross.

*/
int XLALSimIMRPhenomXPHM(
  COMPLEX16FrequencySeries **hptilde, /**< [out] Frequency-domain waveform h+ */
  COMPLEX16FrequencySeries **hctilde, /**< [out] Frequency-domain waveform hx */
  REAL8 m1_SI,                        /**< mass of companion 1 (kg) */
  REAL8 m2_SI,                        /**< mass of companion 2 (kg) */
  REAL8 chi1x,                        /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi1y,                        /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi1z,                        /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2x,                        /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2y,                        /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2z,                        /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  REAL8 distance,                     /**< Distance of source (m) */
  REAL8 inclination,                  /**< inclination of source (rad) */
  REAL8 phiRef,                       /**< Orbital phase (rad) at reference frequency */
  REAL8 f_min,                        /**< Starting GW frequency (Hz) */
  REAL8 f_max,                        /**< Ending GW frequency (Hz); Defaults to Mf = 0.3 if no f_max is specified. */
  REAL8 deltaF,                       /**< Sampling frequency (Hz). To use non-uniform frequency grid, set deltaF <= 0. */
  REAL8 fRef_In,                      /**< Reference frequency (Hz) */
  LALDict *lalParams                  /**< LAL Dictionary struct */
)
{
  /* Variable to check correct calls to functions. */
  UINT4 status;

  /* Check if m1 > m2, swap the bodies otherwise. */
  status = XLALIMRPhenomXPCheckMassesAndSpins(&m1_SI,&m2_SI,&chi1x,&chi1y,&chi1z,&chi2x,&chi2y,&chi2z);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: XLALIMRPhenomXPCheckMassesAndSpins failed.\n");

  #if DEBUG == 1
  printf("fRef_In : %e\n",fRef_In);
  printf("m1_SI   : %e\n",m1_SI);
  printf("m2_SI   : %e\n",m2_SI);
  printf("chi1L   : %e\n",chi1z);
  printf("chi2L   : %e\n",chi2z);
  printf("phiRef  : %e\n",phiRef);
  printf("Prec V. : %d\n\n",XLALSimInspiralWaveformParamsLookupPhenomXPrecVersion(lalParams));
  printf("Performing sanity checks...\n");
  #endif

  /* Perform initial sanity checks */
  XLAL_CHECK(NULL != hptilde, XLAL_EFAULT, "Error: hptilde already defined.                        \n");
  XLAL_CHECK(NULL != hctilde, XLAL_EFAULT, "Error: hctilde already defined.                        \n");
  XLAL_CHECK(fRef_In  >= 0, XLAL_EFUNC,    "Error: fRef_In must be positive or set to 0 to ignore. \n");
  XLAL_CHECK(deltaF   >  0, XLAL_EFUNC,    "Error: deltaF must be positive and greater than 0.     \n");
  XLAL_CHECK(m1_SI    >  0, XLAL_EFUNC,    "Error: m1 must be positive and greater than 0.         \n");
  XLAL_CHECK(m2_SI    >  0, XLAL_EFUNC,    "Error: m2 must be positive and greater than 0.         \n");
  XLAL_CHECK(f_min    >  0, XLAL_EFUNC,    "Error: f_min must be positive and greater than 0.      \n");
  XLAL_CHECK(f_max    >= 0, XLAL_EFUNC,    "Error: f_max must be non-negative.                     \n");
  XLAL_CHECK(distance >  0, XLAL_EFUNC,    "Error: Distance must be positive and greater than 0.   \n");

  /*
  Perform a basic sanity check on the region of the parameter space in which model is evaluated.
  Behaviour is as follows, consistent with the choices for IMRPhenomXAS/IMRPhenomXHM
    - For mass ratios <= 20.0 and spins <= 0.99: no warning messages.
    - For 1000 > mass ratio > 20 and spins <= 0.99: print a warning message that we are extrapolating outside of *NR* calibration domain.
    - For mass ratios > 1000: throw a hard error that model is not valid.
    - For spins > 0.99: throw a warning that we are extrapolating the model to extremal

  */
  REAL8 mass_ratio;
  if(m1_SI > m2_SI)
  {
    mass_ratio = m1_SI / m2_SI;
  }
  else
  {
    mass_ratio = m2_SI / m1_SI;
  }
  if(mass_ratio > 20.0  ) { XLAL_PRINT_WARNING("Warning: Extrapolating outside of Numerical Relativity calibration domain. NNLO angles may become pathological at large mass ratios.\n"); }
  if(mass_ratio > 1000. && fabs(mass_ratio - 1000) > 1e-12) { XLAL_ERROR(XLAL_EDOM, "ERROR: Model not valid at mass ratios beyond 1000.\n"); } // The 1e-12 is to avoid rounding errors
  if(fabs(chi1z) > 0.99 || fabs(chi2z) > 0.99) { XLAL_PRINT_WARNING("Warning: Extrapolating to extremal spins, model is not trusted.\n"); }

  /* Check that the modes chosen are available for the model */
  XLAL_CHECK(check_input_mode_array(lalParams) == XLAL_SUCCESS, XLAL_EFAULT, "Not available mode chosen.\n");

  /* If no reference frequency is given, set it to the starting gravitational wave frequency. */
  REAL8 fRef = (fRef_In == 0.0) ? f_min : fRef_In;
  
  /* Use an auxiliar laldict to not overwrite the input argument */
  LALDict *lalParams_aux;
  /* setup mode array */
  if (lalParams == NULL)
  {
      lalParams_aux = XLALCreateDict();
  }
  else{
      lalParams_aux = XLALDictDuplicate(lalParams);
  }
  
  #if DEBUG == 1
  printf("\n\n **** Initializing waveform struct... **** \n\n");
  #endif

  /* Initialize the useful powers of LAL_PI */
  status = IMRPhenomX_Initialize_Powers(&powers_of_lalpi, LAL_PI);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initialize useful powers of LAL_PI.\n");

  /* Initialize IMRPhenomX Waveform struct and check that it initialized correctly */
  IMRPhenomXWaveformStruct *pWF;
  pWF    = XLALMalloc(sizeof(IMRPhenomXWaveformStruct));
  status = IMRPhenomXSetWaveformVariables(pWF, m1_SI, m2_SI, chi1z, chi2z, deltaF, fRef, phiRef, f_min, f_max, distance, inclination, lalParams_aux, PHENOMXDEBUG);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetWaveformVariables failed.\n");

  /*
      Create a REAL8 frequency series.
      Use fLow, fHigh, deltaF to compute frequency sequence. Only pass the boundaries (fMin, fMax).
  */
  REAL8Sequence *freqs = XLALCreateREAL8Sequence(2);
  freqs->data[0] = pWF->fMin;
  freqs->data[1] = pWF->f_max_prime;

  #if DEBUG == 1
  printf("\n\n **** Initializing precession struct... **** \n\n");
  #endif

  /* Initialize IMRPhenomX Precession struct and check that it generated successfully */
  IMRPhenomXPrecessionStruct *pPrec;
  pPrec  = XLALMalloc(sizeof(IMRPhenomXPrecessionStruct));

  status = IMRPhenomXGetAndSetPrecessionVariables(
             pWF,
             pPrec,
             m1_SI,
             m2_SI,
             chi1x,
             chi1y,
             chi1z,
             chi2x,
             chi2y,
             chi2z,
             lalParams_aux,
             PHENOMXDEBUG
           );
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetPrecessionVariables failed.\n");


  #if DEBUG == 1
  printf("\n\n **** Calling IMRPhenomXPHM_hplushcross... **** \n\n");
  #endif
  
  /* We now call the core IMRPhenomXPHM waveform generator */
  status = IMRPhenomXPHM_hplushcross(hptilde, hctilde, freqs, pWF, pPrec, lalParams_aux);
  XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "IMRPhenomXPHM_hplushcross failed to generate IMRPhenomXHM waveform.\n");

  #if DEBUG == 1
  printf("\n\n **** Call to IMRPhenomXPHM_hplus_hcross complete. **** \n\n");
  #endif


  /* Resize hptilde, hctilde */
  REAL8 lastfreq;
  if (pWF->f_max_prime < pWF->fMax)
  {
    /* The user has requested a higher f_max than Mf = fCut.
    Resize the frequency series to fill with zeros beyond the cutoff frequency. */
    lastfreq = pWF->fMax;
    XLAL_PRINT_WARNING("The input f_max = %.2f Hz is larger than the internal cutoff of Mf=0.3 (%.2f Hz). Array will be filled with zeroes between these two frequencies.\n", pWF->fMax, pWF->f_max_prime);
  }
  else{  // We have to look for a power of 2 anyway.
    lastfreq = pWF->f_max_prime;
  }
  // We want to have the length be a power of 2 + 1
  size_t n_full = NextPow2(lastfreq / deltaF) + 1;
  size_t n = (*hptilde)->data->length;

  /* Resize the COMPLEX16 frequency series */
  *hptilde = XLALResizeCOMPLEX16FrequencySeries(*hptilde, 0, n_full);
  XLAL_CHECK (*hptilde, XLAL_ENOMEM, "Failed to resize h_+ COMPLEX16FrequencySeries of length %zu (for internal fCut=%f) to new length %zu (for user-requested f_max=%f).", n, pWF->fCut, n_full, pWF->fMax );

  /* Resize the COMPLEX16 frequency series */
  *hctilde = XLALResizeCOMPLEX16FrequencySeries(*hctilde, 0, n_full);
  XLAL_CHECK (*hctilde, XLAL_ENOMEM, "Failed to resize h_x COMPLEX16FrequencySeries of length %zu (for internal fCut=%f) to new length %zu (for user-requested f_max=%f).", n, pWF->fCut, n_full, pWF->fMax );


  /* Free memory */
  LALFree(pWF);
  LALFree(pPrec);
  XLALDestroyREAL8Sequence(freqs);
  XLALDestroyDict(lalParams_aux);
  
  return XLAL_SUCCESS;
}

/** Returns hptilde and hctilde of the multimode precessing waveform for positive frequencies in an equally space grid.

This function is equivalent to XLALSimIMRPhenomXPHM, gives the same result but it is computed by first calling
all the precessing indiviudal modes in the J-frame and then sum them all with Ylm(thetaJN, 0) since the J-frame is aligned
such that the line of sight N is in the x-z plane of the J-frame. See appendix C and in particular eq. C8 in Precessing paper.

This function is slower since it has to compute the non-precessing modes again for every precessing mode.

It is just a wrapper of the internal function that actually carries out the calculation: IMRPhenomXPHM_hplushcross_from_modes.

*/
int XLALSimIMRPhenomXPHMFromModes(
  COMPLEX16FrequencySeries **hptilde, /**< [out] Frequency-domain waveform h+ */
  COMPLEX16FrequencySeries **hctilde, /**< [out] Frequency-domain waveform hx */
  REAL8 m1_SI,                        /**< mass of companion 1 (kg) */
  REAL8 m2_SI,                        /**< mass of companion 2 (kg) */
  REAL8 chi1x,                        /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi1y,                        /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi1z,                        /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2x,                        /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2y,                        /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2z,                        /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  REAL8 distance,                     /**< Distance of source (m) */
  REAL8 inclination,                  /**< inclination of source (rad) */
  REAL8 phiRef,                       /**< Orbital phase (rad) at reference frequency */
  REAL8 f_min,                        /**< Starting GW frequency (Hz) */
  REAL8 f_max,                        /**< Ending GW frequency (Hz); Defaults to Mf = 0.3 if no f_max is specified. */
  REAL8 deltaF,                       /**< Sampling frequency (Hz). To use non-uniform frequency grid, set deltaF <= 0. */
  REAL8 fRef_In,                      /**< Reference frequency (Hz) */
  LALDict *lalParams                  /**< LAL Dictionary struct */
)
{
  /* Variable to check correct calls to functions. */
  UINT4 status;

  /* Check if m1 > m2, swap the bodies otherwise. */
  status = XLALIMRPhenomXPCheckMassesAndSpins(&m1_SI,&m2_SI,&chi1x,&chi1y,&chi1z,&chi2x,&chi2y,&chi2z);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: XLALIMRPhenomXPCheckMassesAndSpins failed.\n");

  #if DEBUG == 1
  printf("fRef_In : %e\n",fRef_In);
  printf("m1_SI   : %e\n",m1_SI);
  printf("m2_SI   : %e\n",m2_SI);
  printf("chi1L   : %e\n",chi1z);
  printf("chi2L   : %e\n",chi2z);
  printf("phiRef  : %e\n",phiRef);
  printf("Prec V. : %d\n\n",XLALSimInspiralWaveformParamsLookupPhenomXPrecVersion(lalParams));
  printf("Performing sanity checks...\n");
  #endif

  /* Perform initial sanity checks */
  XLAL_CHECK(NULL != hptilde, XLAL_EFAULT, "Error: hptilde already defined.                        \n");
  XLAL_CHECK(NULL != hctilde, XLAL_EFAULT, "Error: hctilde already defined.                        \n");
  XLAL_CHECK(fRef_In  >= 0, XLAL_EFUNC,    "Error: fRef_In must be positive or set to 0 to ignore. \n");
  XLAL_CHECK(deltaF   >  0, XLAL_EFUNC,    "Error: deltaF must be positive and greater than 0.     \n");
  XLAL_CHECK(m1_SI    >  0, XLAL_EFUNC,    "Error: m1 must be positive and greater than 0.         \n");
  XLAL_CHECK(m2_SI    >  0, XLAL_EFUNC,    "Error: m2 must be positive and greater than 0.         \n");
  XLAL_CHECK(f_min    >  0, XLAL_EFUNC,    "Error: f_min must be positive and greater than 0.      \n");
  XLAL_CHECK(f_max    >= 0, XLAL_EFUNC,    "Error: f_max must be non-negative.                     \n");
  XLAL_CHECK(distance >  0, XLAL_EFUNC,    "Error: Distance must be positive and greater than 0.   \n");

  /*
  Perform a basic sanity check on the region of the parameter space in which model is evaluated.
  Behaviour is as follows, consistent with the choices for IMRPhenomXAS/IMRPhenomXHM
    - For mass ratios <= 20.0 and spins <= 0.99: no warning messages.
    - For 1000 > mass ratio > 20 and spins <= 0.99: print a warning message that we are extrapolating outside of *NR* calibration domain.
    - For mass ratios > 1000: throw a hard error that model is not valid.
    - For spins > 0.99: throw a warning that we are extrapolating the model to extremal

  */
  REAL8 mass_ratio;
  if(m1_SI > m2_SI)
  {
    mass_ratio = m1_SI / m2_SI;
  }
  else
  {
    mass_ratio = m2_SI / m1_SI;
  }
  if(mass_ratio > 20.0  ) { XLAL_PRINT_WARNING("Warning: Extrapolating outside of Numerical Relativity calibration domain. NNLO angles may become pathological at large mass ratios.\n"); }
  if(mass_ratio > 1000. && fabs(mass_ratio - 1000) > 1e-12) { XLAL_ERROR(XLAL_EDOM, "ERROR: Model not valid at mass ratios beyond 1000."); } // The 1e-12 is to avoid rounding errors
  if(fabs(chi1z) > 0.99 || fabs(chi2z) > 0.99) { XLAL_PRINT_WARNING("Warning: Extrapolating to extremal spins, model is not trusted."); }

  /* Check that the modes chosen are available for the model */
  XLAL_CHECK(check_input_mode_array(lalParams) == XLAL_SUCCESS, XLAL_EFAULT, "Not available mode chosen.\n");

  /* If no reference frequency is given, set it to the starting gravitational wave frequency */
  REAL8 fRef = (fRef_In == 0.0) ? f_min : fRef_In;
  
  /* Use an auxiliar laldict to not overwrite the input argument */
  LALDict *lalParams_aux;
  /* setup mode array */
  if (lalParams == NULL)
  {
      lalParams_aux = XLALCreateDict();
  }
  else{
      lalParams_aux = XLALDictDuplicate(lalParams);
  }
  
  #if DEBUG == 1
  printf("\n\n **** Initializing waveform struct... **** \n\n");
  #endif

  /* Initialize the useful powers of LAL_PI */
  status = IMRPhenomX_Initialize_Powers(&powers_of_lalpi, LAL_PI);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initialize useful powers of LAL_PI.");

  /* Initialize IMR PhenomX Waveform struct and check that it initialized correctly */
  IMRPhenomXWaveformStruct *pWF;
  pWF    = XLALMalloc(sizeof(IMRPhenomXWaveformStruct));
  status = IMRPhenomXSetWaveformVariables(pWF, m1_SI, m2_SI, chi1z, chi2z, deltaF, fRef, phiRef, f_min, f_max, distance, inclination, lalParams_aux, PHENOMXDEBUG);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetWaveformVariables failed.\n");

  /*
      Create a REAL8 frequency series.
      Use fLow, fHigh, deltaF to compute frequency sequence. Only pass the boundaries (fMin, fMax).
  */
  REAL8Sequence *freqs = XLALCreateREAL8Sequence(2);
  freqs->data[0] = pWF->fMin;
  freqs->data[1] = pWF->f_max_prime;

  #if DEBUG == 1
  printf("\n\n **** Initializing precession struct... **** \n\n");
  #endif

  /* Initialize IMRPhenomX Precession struct and check that it generated successfully. */
  IMRPhenomXPrecessionStruct *pPrec;
  pPrec  = XLALMalloc(sizeof(IMRPhenomXPrecessionStruct));

  status = IMRPhenomXGetAndSetPrecessionVariables(
             pWF,
             pPrec,
             m1_SI,
             m2_SI,
             chi1x,
             chi1y,
             chi1z,
             chi2x,
             chi2y,
             chi2z,
             lalParams_aux,
             PHENOMXDEBUG
           );
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetPrecessionVariables failed.\n");


  #if DEBUG == 1
  printf("\n\n **** Calling IMRPhenomXPHM_hplushcross_from_modes... **** \n\n");
  #endif

  /* We now call the core IMRPhenomXPHMFromModes waveform generator */
  status = IMRPhenomXPHM_hplushcross_from_modes(hptilde, hctilde, freqs, pWF, pPrec, lalParams_aux);
  XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "IMRPhenomXPHM_hplushcross_from_modes failed to generate IMRPhenomXPHM waveform.");

  #if DEBUG == 1
  printf("\n\n **** Call to IMRPhenomXPHM_from _modes complete. **** \n\n");
  #endif

  /* Resize hptilde, hctilde */
  REAL8 lastfreq;
  if (pWF->f_max_prime < pWF->fMax)
  {
    /* The user has requested a higher f_max than Mf = fCut.
    Resize the frequency series to fill with zeros beyond the cutoff frequency. */
    lastfreq = pWF->fMax;
    XLAL_PRINT_WARNING("The input f_max = %.2f Hz is larger than the internal cutoff of Mf=0.3 (%.2f Hz). Array will be filled with zeroes between these two frequencies.\n", pWF->fMax, pWF->f_max_prime);
  }
  else{  // We have to look for a power of 2 anyway.
    lastfreq = pWF->f_max_prime;
  }
  // We want to have the length be a power of 2 + 1
  size_t n_full = NextPow2(lastfreq / deltaF) + 1;
  size_t n = (*hptilde)->data->length;

  /* Resize the COMPLEX16 frequency series */
  *hptilde = XLALResizeCOMPLEX16FrequencySeries(*hptilde, 0, n_full);
  XLAL_CHECK (*hptilde, XLAL_ENOMEM, "Failed to resize h_+ COMPLEX16FrequencySeries of length %zu (for internal fCut=%f) to new length %zu (for user-requested f_max=%f).", n, pWF->fCut, n_full, pWF->fMax );

  /* Resize the COMPLEX16 frequency series */
  *hctilde = XLALResizeCOMPLEX16FrequencySeries(*hctilde, 0, n_full);
  XLAL_CHECK (*hctilde, XLAL_ENOMEM, "Failed to resize h_x COMPLEX16FrequencySeries of length %zu (for internal fCut=%f) to new length %zu (for user-requested f_max=%f).", n, pWF->fCut, n_full, pWF->fMax );


  /* Free memory */
  LALFree(pWF);
  LALFree(pPrec);
  XLALDestroyREAL8Sequence(freqs);
  XLALDestroyDict(lalParams_aux);

  return XLAL_SUCCESS;
}


/**
 * Returns hptilde and hctilde as a complex frequency series with entries exactly at the frequencies specified in
 * the REAL8Sequence *freqs (which can be unequally spaced). No zeros are added. Assumes positive frequencies.
 * This is the function used when calling ChooseFDWaveformSequence.
 * It is a wrapper that calls the function IMRPhenomXPHM_hplushcross or IMRPhenomXPHM_hplushcross_from_modes
 * and we can change which one is used by modifying the option 'UseModes' in the LAL dictionary.
 * No multibanding is used since this technique is only for equal-spaced frequency grids.
 */
 int XLALSimIMRPhenomXPHMFrequencySequence(
    COMPLEX16FrequencySeries **hptilde, /**< [out] Frequency-domain waveform h+  */
    COMPLEX16FrequencySeries **hctilde, /**< [out] Frequency-domain waveform hx  */
    REAL8Sequence *freqs,               /**< Input Frequency series [Hz]         */
    REAL8 m1_SI,                        /**< mass of companion 1 (kg) */
    REAL8 m2_SI,                        /**< mass of companion 2 (kg) */
    REAL8 chi1x,                        /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 chi1y,                        /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 chi1z,                        /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 chi2x,                        /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    REAL8 chi2y,                        /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    REAL8 chi2z,                        /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    REAL8 distance,                     /**< Distance of source (m) */
    REAL8 inclination,                  /**< inclination of source (rad) */
    REAL8 phiRef,                       /**< Orbital phase (rad) at reference frequency */
    REAL8 fRef_In,                      /**< Reference frequency (Hz) */
    LALDict *lalParams                  /**< LAL Dictionary struct */
)
{
  /* Variable to check correct calls to functions. */
  INT4 status;

  /* Check if m1 > m2, swap the bodies otherwise. */
  status = XLALIMRPhenomXPCheckMassesAndSpins(&m1_SI,&m2_SI,&chi1x,&chi1y,&chi1z,&chi2x,&chi2y,&chi2z);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: XLALIMRPhenomXPCheckMassesAndSpins failed.\n");

  /* Perform initial sanity checks */
  XLAL_CHECK(NULL != hptilde, XLAL_EFAULT, "Error: hptilde already defined.                        \n");
  XLAL_CHECK(NULL != hctilde, XLAL_EFAULT, "Error: hctilde already defined.                        \n");
  XLAL_CHECK(NULL != freqs,   XLAL_EFAULT, "Error: Input freq array must be defined.               \n");
  XLAL_CHECK(fRef_In  >= 0, XLAL_EFUNC,    "Error: fRef must be positive and greater than 0.       \n");
  XLAL_CHECK(m1_SI    >  0, XLAL_EFUNC,    "Error: m1 must be positive and greater than 0.         \n");
  XLAL_CHECK(m2_SI    >  0, XLAL_EFUNC,    "Error: m2 must be positive and greater than 0.         \n");
  XLAL_CHECK(distance >  0, XLAL_EFUNC,    "Error: Distance must be positive and greater than 0.   \n");

  /*
  Perform a basic sanity check on the region of the parameter space in which model is evaluated.
  Behaviour is as follows, consistent with the choices for IMRPhenomXAS/IMRPhenomXHM
    - For mass ratios <= 20.0 and spins <= 0.99: no warning messages.
    - For 1000 > mass ratio > 20 and spins <= 0.99: print a warning message that we are extrapolating outside of *NR* calibration domain.
    - For mass ratios > 1000: throw a hard error that model is not valid.
    - For spins > 0.99: throw a warning that we are extrapolating the model to extremal

  */
  REAL8 mass_ratio;
  if(m1_SI > m2_SI)
  {
    mass_ratio = m1_SI / m2_SI;
  }
  else
  {
    mass_ratio = m2_SI / m1_SI;
  }
  if(mass_ratio > 20.0  ) { XLAL_PRINT_INFO("Warning: Extrapolating outside of Numerical Relativity calibration domain."); }
  if(mass_ratio > 1000. && fabs(mass_ratio - 1000) > 1e-12) { XLAL_ERROR(XLAL_EDOM, "ERROR: Model not valid at mass ratios beyond 1000."); } // The 1e-12 is to avoid rounding errors
  if(fabs(chi1z) > 0.99 || fabs(chi2z) > 0.99) { XLAL_PRINT_INFO("Warning: Extrapolating to extremal spins, model is not trusted."); }

  /* Check that the modes chosen are available for the model */
  XLAL_CHECK(check_input_mode_array(lalParams) == XLAL_SUCCESS, XLAL_EFAULT, "Not available mode chosen.\n");

  /* If no reference frequency is given, set it to the starting gravitational wave frequency */
  REAL8 fRef = (fRef_In == 0.0) ? freqs->data[0] : fRef_In; //It is giving valgrind error, but it is not needed. f_ref = f_min in WaveformCache.c and SimInspiral.c.

  /* Get minimum and maximum frequencies. */
  REAL8 f_min_In  = freqs->data[0];
  REAL8 f_max_In  = freqs->data[freqs->length - 1];

  /*
    Passing deltaF = 0 implies that freqs is a frequency grid with non-uniform spacing.
    The function waveform then start at lowest given frequency.
    The Multibanding has to be switched off since it is built only for equally spaced frequency grid.
  */
  /* Use an auxiliar laldict to not overwrite the input argument */
  LALDict *lalParams_aux;
  /* setup mode array */
  if (lalParams == NULL)
  {
      lalParams_aux = XLALCreateDict();
  }
  else{
      lalParams_aux = XLALDictDuplicate(lalParams);
  }
  if(XLALSimInspiralWaveformParamsLookupPhenomXHMThresholdMband(lalParams_aux)!=0 || XLALSimInspiralWaveformParamsLookupPhenomXPHMThresholdMband(lalParams_aux)!=0)
  {
    XLAL_PRINT_WARNING("Warning: Function is aimed for non-uniform frequency grid, switching off Multibanding.");
    XLALSimInspiralWaveformParamsInsertPhenomXHMThresholdMband(lalParams_aux, 0);
    XLALSimInspiralWaveformParamsInsertPhenomXPHMThresholdMband(lalParams_aux, 0);
  }

  status = IMRPhenomX_Initialize_Powers(&powers_of_lalpi, LAL_PI);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initialize useful powers of LAL_PI.");

  /* Initialize IMRPhenomX waveform struct and perform sanity check. */
  IMRPhenomXWaveformStruct *pWF;
  pWF = XLALMalloc(sizeof(IMRPhenomXWaveformStruct));
  status = IMRPhenomXSetWaveformVariables(pWF, m1_SI, m2_SI, chi1z, chi2z, 0.0, fRef, phiRef, f_min_In, f_max_In, distance, inclination, lalParams_aux, DEBUG);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetWaveformVariables failed.\n");

  /* Initialize IMR PhenomX Precession struct and check that it generated successfully */
  IMRPhenomXPrecessionStruct *pPrec;
  pPrec  = XLALMalloc(sizeof(IMRPhenomXPrecessionStruct));

  status = IMRPhenomXGetAndSetPrecessionVariables(
           pWF,
           pPrec,
           m1_SI,
           m2_SI,
           chi1x,
           chi1y,
           chi1z,
           chi2x,
           chi2y,
           chi2z,
           lalParams_aux,
           PHENOMXDEBUG
         );
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetPrecessionVariables failed.\n");


   /* Now call the core IMRPhenomXPHM waveform generator */
   /* Choose between the two options for computing the polarizations. The default is 0. */
   if(XLALSimInspiralWaveformParamsLookupPhenomXPHMUseModes(lalParams_aux) == 0)
   {
     status = IMRPhenomXPHM_hplushcross(hptilde, hctilde, freqs, pWF, pPrec, lalParams_aux);
   }
   else
   {
     status = IMRPhenomXPHM_hplushcross_from_modes(hptilde, hctilde, freqs, pWF, pPrec, lalParams_aux);
   }
   XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "IMRPhenomXPHM_hplushcross failed to generate IMRPhenomXPHM waveform.");

   /* Free memory */
   LALFree(pPrec);
   LALFree(pWF);
   XLALDestroyDict(lalParams_aux);

   return XLAL_SUCCESS;
 }
 /** @} **
 * @} **/


 /**
  Core function of XLALSimIMRPhenomXPHM and XLALSimIMRPhenomXPHMFrequencySequence.
  Returns hptilde, hctilde for positive frequencies.
  The default non-precessing modes are 2|2|, 2|1|, 3|3|, 3|2| and 4|4|.
  It returns also the contribution of the corresponding negatives modes.
  It can be evaulated in a non-uniform frequency grid. Assumes positive frequencies.
 */
 static int IMRPhenomXPHM_hplushcross(
   COMPLEX16FrequencySeries **hptilde,  /**< [out] Frequency domain h+ GW strain */
   COMPLEX16FrequencySeries **hctilde,  /**< [out] Frequency domain hx GW strain */
   REAL8Sequence *freqs_In,             /**< Frequency array to evaluate the model. (fmin, fmax) for equally spaced grids. */
   IMRPhenomXWaveformStruct *pWF,       /**< IMRPhenomX Waveform Struct  */
   IMRPhenomXPrecessionStruct *pPrec,   /**< IMRPhenomXP Precession Struct  */
   LALDict *lalParams                   /**< LAL Dictionary Structure    */
 )
 {
   /* Set LIGOTimeGPS */
   LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0,0}

   REAL8 deltaF = pWF->deltaF;

   lalParams = IMRPhenomXPHM_setup_mode_array(lalParams);
   LALValue *ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(lalParams);

   /* At this point ModeArray should contain the list of modes
   and therefore if NULL then something is wrong and abort. */
   if (ModeArray == NULL)
   {
     XLAL_ERROR(XLAL_EDOM, "ModeArray is NULL when it shouldn't be. Aborting.\n");
   }

   INT4 status = 0; //Variable to check correct functions calls.

   /* Build the frequency array and initialize hptilde to the length of freqs. */
   REAL8Sequence *freqs;
   UINT4 offset = SetupWFArrays(&freqs, hptilde, freqs_In, pWF, ligotimegps_zero);

   /* Initialize hctilde according to hptilde. */
   size_t npts = (*hptilde)->data->length;
   *hctilde = XLALCreateCOMPLEX16FrequencySeries("hctilde: FD waveform", &(*hptilde)->epoch, (*hptilde)->f0, pWF->deltaF, &lalStrainUnit, npts);
   XLAL_CHECK (*hctilde, XLAL_ENOMEM, "Failed to allocated waveform COMPLEX16FrequencySeries of length %zu.", npts);
   memset((*hctilde)->data->data, 0, npts * sizeof(COMPLEX16));
   XLALUnitMultiply(&((*hctilde)->sampleUnits), &((*hctilde)->sampleUnits), &lalSecondUnit);

   /* Object to store the non-precessing 22 mode waveform and to be recycled when calling the 32 mode in multibanding. */
   COMPLEX16FrequencySeries *htilde22 = NULL;

   /*
      Take input/default value for the threshold of the Multibanding for the hlms modes.
      If = 0 then do not use Multibanding. Default value defined in XLALSimInspiralWaveformParams.c.
      If the input freqs_In is non-uniform the Multibanding has been already switche off.
   */
   REAL8 thresholdMB  = XLALSimInspiralWaveformParamsLookupPhenomXHMThresholdMband(lalParams);

   /* Initialize the power of pi for the HM internal functions. */
   status = IMRPhenomX_Initialize_Powers(&powers_of_lalpiHM, LAL_PI);
   XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Failed to initialize useful powers of LAL_PI.");


   SphHarmFrequencySeries **hlms = XLALMalloc(sizeof(SphHarmFrequencySeries));
   *hlms = NULL;
   if (XLALSimInspiralWaveformParamsLookupPhenomXPHMTwistPhenomHM(lalParams)==1)
   {
     /* evaluate all hlm modes */
      status = XLALSimIMRPhenomHMGethlmModes(
         hlms,
         freqs,
         pWF->m1_SI,
         pWF->m2_SI,
         pPrec->chi1x,
         pPrec->chi1y,
         pWF->chi1L,
         pPrec->chi2x,
         pPrec->chi2y,
         pWF->chi2L,
         pWF->phi0,
         //pWF->deltaF,
         0,
         pWF->fRef,
         lalParams);
     XLAL_CHECK(XLAL_SUCCESS == status,
                XLAL_EFUNC, "XLALSimIMRPhenomHMGethlmModes failed");
   }




   /***** Loop over non-precessing modes ******/
   for (UINT4 ell = 2; ell <= L_MAX; ell++)
   {
     for (UINT4 emmprime = 1; emmprime <= ell; emmprime++)
     {
       /* Loop over only positive mprime is intentional.
          The single mode function returns the negative mode h_l-mprime, and the positive
          is added automatically in during the twisting up in IMRPhenomXPHMTwistUp.
          First check if (l,m) mode is 'activated' in the ModeArray.
          If activated then generate the mode, else skip this mode.
       */
       if (XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, emmprime) != 1)
       { /* skip mode */
         continue;
       } /* else: generate mode */

       #if DEBUG == 1
       printf("\n\n*********************************\n*Non-precessing Mode %i%i \n******************************\n",ell, emmprime);
       // Save the hlm mode into a file
       FILE *fileangle;
       char fileSpec[40];

       if(pPrec->MBandPrecVersion == 0)
       {
           sprintf(fileSpec, "angles_hphc_%i%i.dat", ell, emmprime);
       }
       else
       {
           sprintf(fileSpec, "angles_hphc_MB_%i%i.dat", ell, emmprime);
       }
       printf("\nOutput angle file: %s\r\n", fileSpec);
       fileangle = fopen(fileSpec,"w");

       fprintf(fileangle,"# q = %.16e m1 = %.16e m2 = %.16e chi1 = %.16e chi2 = %.16e lm = %i%i Mtot = %.16e distance = %.16e\n", pWF->q, pWF->m1, pWF->m2, pWF->chi1L, pWF->chi2L, ell, emmprime, pWF->Mtot, pWF->distance/LAL_PC_SI/1e6);
       fprintf(fileangle,"#fHz   cexp_i_alpha(re im)   cexp_i_epsilon(re im)    cexp_i_betah(re im)\n");

       fclose(fileangle);
       #endif

       /* Variable to store the strain of only one (negative) mode: h_l-mprime */
       COMPLEX16FrequencySeries *htildelm = NULL;

       if (XLALSimInspiralWaveformParamsLookupPhenomXPHMTwistPhenomHM(lalParams)==1)
       {
          INT4 minus1l = 1;
          if(ell % 2 !=0) minus1l = -1;
          COMPLEX16FrequencySeries *htildelmPhenomHM = NULL;
          /* Initialize the htilde frequency series */
          htildelm = XLALCreateCOMPLEX16FrequencySeries("htildelm: FD waveform", &ligotimegps_zero, 0, pWF->deltaF, &lalStrainUnit, npts);
          /* Check that frequency series generated okay */
          XLAL_CHECK(htildelm,XLAL_ENOMEM,"Failed to allocate COMPLEX16FrequencySeries of length %zu for f_max = %f, deltaF = %g.\n", npts, freqs_In->data[freqs_In->length - 1], pWF->deltaF);
          memset((htildelm)->data->data, 0, npts * sizeof(COMPLEX16));
          XLALUnitMultiply(&((htildelm)->sampleUnits), &((htildelm)->sampleUnits), &lalSecondUnit);

          htildelmPhenomHM = XLALSphHarmFrequencySeriesGetMode(*hlms, ell, emmprime);
          for(UINT4 idx = 0; idx < freqs->length; idx++)
          {
            htildelm->data->data[idx+offset] = minus1l * htildelmPhenomHM->data->data[idx] * pWF->amp0;
          }
          //XLALDestroyCOMPLEX16FrequencySeries(htildelmPhenomHM);
       }
       else
       {
         /* Compute non-precessing mode */
         if (thresholdMB == 0){  // No multibanding
           if(ell == 2 && emmprime == 2)
           {
              status = IMRPhenomXASGenerateFD(&htildelm, freqs, pWF, lalParams);
              XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "IMRPhenomXASGenerateFD failed to generate IMRPhenomXHM waveform.");
           }
           else
           {
             status = IMRPhenomXHMGenerateFDOneMode(&htildelm, freqs, pWF, ell, emmprime, lalParams);
             XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "IMRPhenomXHMGenerateFDOneMode failed to generate IMRPhenomXHM waveform.");
           }
         }
         else{               // With multibanding
           if(ell==3 && emmprime==2){  // mode with mode-mixing
              status = IMRPhenomXHMMultiBandOneModeMixing(&htildelm, htilde22, pWF, ell, emmprime, lalParams);
              XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "IMRPhenomXHMMultiBandOneModeMixing failed to generate IMRPhenomXHM waveform.");
           }
           else{                  // modes without mode-mixing including 22 mode
              status = IMRPhenomXHMMultiBandOneMode(&htildelm, pWF, ell, emmprime, lalParams);
              XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "IMRPhenomXHMMultiBandOneMode failed to generate IMRPhenomXHM waveform.");
           }

           /* IMRPhenomXHMMultiBandOneMode* functions set pWF->deltaF=0 internally, we put it back here. */
           pWF->deltaF = deltaF;

           /* If the 22 and 32 modes are active, we recycle the 22 mode for the mixing in the 32 and it is passed to IMRPhenomXHMMultiBandOneModeMixing.
              The 22 mode is always computed first than the 32, we store the 22 mode in the variable htilde22. */
           if(ell==2 && emmprime==2 && XLALSimInspiralModeArrayIsModeActive(ModeArray, 3, 2)==1){
             htilde22 = XLALCreateCOMPLEX16FrequencySeries("hptilde: FD waveform", &(ligotimegps_zero), 0.0, pWF->deltaF, &lalStrainUnit, htildelm->data->length);
             for(UINT4 idx = 0; idx < htildelm->data->length; idx++){
               htilde22->data->data[idx] = htildelm->data->data[idx];
             }
           }
         }
      }

      if (!(htildelm)){ XLAL_ERROR(XLAL_EFUNC); }

      /*
         For very special cases of deltaF, it can happen that building htildelm with 'freqs_In' or with 'freqs' gives different lengths.
         In case that happens we resize here to the correct length. We could also have called GenerateFD passing freqs_In,
         but in that ways we would be computing the uniform frequency array twice.
         Alsom doing as here we cover the multibanding and PhenomHM cases.
      */
      if(htildelm->data->length != npts)
      {
        htildelm = XLALResizeCOMPLEX16FrequencySeries(htildelm, 0, npts);
        XLAL_CHECK (htildelm, XLAL_ENOMEM, "Failed to resize hlm COMPLEX16FrequencySeries" );
      }

       /* htildelm is recomputed every time in the loop. Check that it always comes out with the same length */
       XLAL_CHECK (    ((*hptilde)->data->length==htildelm->data->length)
                    && ((*hctilde)->data->length==htildelm->data->length),
                    XLAL_EBADLEN,
                    "Inconsistent lengths between frequency series htildelm (%d), hptilde (%d) and hctilde (%d).",
                    htildelm->data->length, (*hptilde)->data->length, (*hctilde)->data->length
                  );

        /* Skip twisting-up if the non-precessing mode is zero. */
        if((pWF->q == 1) && (pWF->chi1L == pWF->chi2L) && (emmprime % 2 != 0))
        {
           XLALDestroyCOMPLEX16FrequencySeries(htildelm);
          continue;
        }

       /*
                                TWISTING UP
            Transform modes from the precessing L-frame to inertial J-frame.
       */


       /* Variable to store the non-precessing waveform in one frequency point. */
       COMPLEX16 hlmcoprec;

       /* No Multibanding for the angles. */
       if(pPrec->MBandPrecVersion == 0)
       {
         #if DEBUG == 1
         printf("\n****************************************************************\n");
         printf("\n*              NOT USING MBAND FOR ANGLES %i                *\n", offset);
         printf("\n****************************************************************\n");
         #endif

         for (UINT4 idx = 0; idx < freqs->length; idx++)
         {
           double Mf             = pWF->M_sec * freqs->data[idx];
           hlmcoprec             = htildelm->data->data[idx + offset];  /* Co-precessing waveform for one freq point */
           COMPLEX16 hplus       = 0.0;  /* h_+ */
           COMPLEX16 hcross      = 0.0;  /* h_x */

           IMRPhenomXPHMTwistUp(Mf, hlmcoprec, pWF, pPrec, ell, emmprime, &hplus, &hcross);

           (*hptilde)->data->data[idx + offset] += hplus;
           (*hctilde)->data->data[idx + offset] += hcross;
         }
       }
       else
       {
         /*
            Multibanding for the angles.

            - In this first release we use the same coarse grid that is used for computing the non-precessing modes.
            - This grid is discussed in section II-A of arXiv:2001.10897. See also section D of Precessing paper.
            - This grid is computed with the function XLALSimIMRPhenomXMultibandingVersion defined in LALSimIMRPhenomXHM_multiband.c.
            - The version of the coarse grid will be changed with the option 'MBandPrecVersion' defined in LALSimInspiralWaveformParams.c.
            - Currently there is only one version available and the option value for that is 0, which is the default value.
         */

         #if DEBUG == 1
         printf("\n****************************************************************\n");
         printf("\n*                 USING MBAND FOR ANGLES                       *\n");
         printf("\n****************************************************************\n");
         #endif



         /* Compute non-uniform coarse frequency grid as 1D array */
         REAL8Sequence *coarseFreqs;
         XLALSimIMRPhenomXPHMMultibandingGrid(&coarseFreqs, ell, emmprime, pWF, lalParams);

         UINT4 lenCoarseArray = coarseFreqs->length;

         /* Euler angles */
         REAL8 alpha        = 0.0;
         REAL8 epsilon      = 0.0;

         REAL8 cBetah       = 0.0;
         REAL8 sBetah       = 0.0;

         /* Variables to store the Euler angles in the coarse frequency grid. */
         REAL8 *valpha      = (REAL8*)XLALMalloc(lenCoarseArray * sizeof(REAL8));
         REAL8 *vepsilon    = (REAL8*)XLALMalloc(lenCoarseArray * sizeof(REAL8));
         REAL8 *vbetah      = (REAL8*)XLALMalloc(lenCoarseArray * sizeof(REAL8));

         switch(pPrec->IMRPhenomXPrecVersion)
         {
           case 101:
           case 102:
           case 103:
           case 104:
           {
             /* Use NNLO PN Euler angles */
             /* Evaluate angles in coarse freq grid */
             for(UINT4 j=0; j<lenCoarseArray; j++)
             {
               REAL8 Mf = coarseFreqs->data[j];

               /* This function already add the offsets to the angles. */
               Get_alpha_beta_epsilon(&alpha, &cBetah, &sBetah, &epsilon, emmprime, Mf, pPrec, pWF);

               valpha[j]   = alpha;
               vepsilon[j] = epsilon;
               vbetah[j]   = acos(cBetah);
             }
             break;
           }
           case 220:
           case 221:
           case 222:
           case 223:
           case 224:
           {
             /* Use MSA Euler angles. */
             /* Evaluate angles in coarse freq grid */
             for(UINT4 j=0; j<lenCoarseArray; j++)
             {
               /* Get Euler angles. */
               REAL8 Mf = coarseFreqs->data[j];
               const REAL8 v        = cbrt (LAL_PI * Mf * (2.0 / emmprime) );
               const vector vangles = IMRPhenomX_Return_phi_zeta_costhetaL_MSA(v,pWF,pPrec);
               REAL8 cos_beta  = 0.0;

               /* Get the offset for the Euler angles alpha and epsilon. */
               REAL8 alpha_offset_mprime = 0, epsilon_offset_mprime = 0;
               Get_alpha_epsilon_offset(&alpha_offset_mprime, &epsilon_offset_mprime, emmprime, pPrec);

               valpha[j]   = vangles.x - alpha_offset_mprime;
               vepsilon[j] = vangles.y - epsilon_offset_mprime;
               cos_beta    = vangles.z;
               
               status = IMRPhenomXWignerdCoefficients_cosbeta(&cBetah, &sBetah, cos_beta);
               XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "Call to IMRPhenomXWignerdCoefficients_cosbeta failed.");

               vbetah[j]   = acos(cBetah);
             }
             break;
           }
           default:
           {
             XLAL_ERROR(XLAL_EINVAL,"Error: IMRPhenomXPrecVersion not recognized. Recommended default is 223.\n");
             break;
           }
         }

         /*
            We have the three Euler angles evaluated in the coarse frequency grid.
            Now we have to carry out the iterative linear interpolation for the complex exponential of each Euler angle. This follows the procedure of eq. 2.32 in arXiv:2001.10897..
            The result will be three arrays of complex exponential evaluated in the finefreqs.
         */
         UINT4 fine_count = 0, ratio;
         REAL8 Omega_alpha, Omega_epsilon, Omega_betah, Qalpha, Qepsilon, Qbetah;
         REAL8 Mfhere, Mfnext, evaldMf;
         Mfnext = coarseFreqs->data[0];
         evaldMf = XLALSimIMRPhenomXUtilsHztoMf(pWF->deltaF, pWF->Mtot);

         /*
            Number of points where the waveform will be computed.
            It is the same for all the modes and could be computed outside the loop, it is here for clarity since it is not used anywhere else.
         */
         size_t iStop  = (size_t) (pWF->f_max_prime / pWF->deltaF) + 1 - offset;

         UINT4 length_fine_grid = iStop + 3; // This is just to reserve memory, add 3 points of buffer.

         COMPLEX16 *cexp_i_alpha   = (COMPLEX16*)XLALMalloc(length_fine_grid * sizeof(COMPLEX16));
         COMPLEX16 *cexp_i_epsilon = (COMPLEX16*)XLALMalloc(length_fine_grid * sizeof(COMPLEX16));
         COMPLEX16 *cexp_i_betah   = (COMPLEX16*)XLALMalloc(length_fine_grid * sizeof(COMPLEX16));

         #if DEBUG == 1
         printf("\n\nLENGTHS fine grid estimate, coarseFreqs->length = %i %i\n", length_fine_grid, lenCoarseArray);
         printf("fine_count, htildelm->length, offset = %i %i %i\n", fine_count, htildelm->data->length, offset);
         #endif

         /* Loop over the coarse freq points */
         for(UINT4 j = 0; j<lenCoarseArray-1 && fine_count < iStop; j++)
         {
           Mfhere = Mfnext;
           Mfnext = coarseFreqs->data[j+1];

           Omega_alpha   = (valpha[j + 1]   - valpha[j])  /(Mfnext - Mfhere);
           Omega_epsilon = (vepsilon[j + 1] - vepsilon[j])/(Mfnext - Mfhere);
           Omega_betah   = (vbetah[j + 1]   - vbetah[j])  /(Mfnext - Mfhere);

           cexp_i_alpha[fine_count]   = cexp(I*valpha[j]);
           cexp_i_epsilon[fine_count] = cexp(I*vepsilon[j]);
           cexp_i_betah[fine_count]   = cexp(I*vbetah[j]);

           Qalpha   = cexp(I*evaldMf*Omega_alpha);
           Qepsilon = cexp(I*evaldMf*Omega_epsilon);
           Qbetah   = cexp(I*evaldMf*Omega_betah);

           fine_count++;

           REAL8 dratio = (Mfnext-Mfhere)/evaldMf;
           UINT4 ceil_ratio  = ceil(dratio);
           UINT4 floor_ratio = floor(dratio);

           /* Make sure the rounding is done correctly. */
           if(fabs(dratio-ceil_ratio) < fabs(dratio-floor_ratio))
           {
             ratio = ceil_ratio;
           }
           else
           {
             ratio = floor_ratio;
           }

          /* Compute complex exponential in fine points between two coarse points */
          /* This loop carry out the eq. 2.32 in arXiv:2001.10897 */
           for(UINT4 kk = 1; kk < ratio && fine_count < iStop; kk++){
             cexp_i_alpha[fine_count]   = Qalpha*cexp_i_alpha[fine_count-1];
             cexp_i_epsilon[fine_count] = Qepsilon*cexp_i_epsilon[fine_count-1];
             cexp_i_betah[fine_count]   = Qbetah*cexp_i_betah[fine_count-1];
             fine_count++;
           }
         }// Loop over coarse grid

         /*
          Now we have the complex exponentials of the three Euler angles alpha, beta, epsilon evaluated in the fine frequency grid.
          Next step is do the twisting up with these.
         */

        #if DEBUG == 1
        printf("fine_count, htildelm->length, offset = %i %i %i\n", fine_count, htildelm->data->length, offset);
        #endif

         /************** TWISTING UP in the fine grid *****************/
         for (UINT4 idx = 0; idx < fine_count; idx++)
         {
           double Mf   = pWF->M_sec * (idx + offset)*pWF->deltaF;

           hlmcoprec   = htildelm->data->data[idx + offset];  /* Co-precessing waveform */

           COMPLEX16 hplus       = 0.0;  /* h_+ */
           COMPLEX16 hcross      = 0.0;  /* h_x */

           pPrec->cexp_i_alpha   = cexp_i_alpha[idx];
           pPrec->cexp_i_epsilon = cexp_i_epsilon[idx];
           pPrec->cexp_i_betah   = cexp_i_betah[idx];

           IMRPhenomXPHMTwistUp(Mf, hlmcoprec, pWF, pPrec, ell, emmprime, &hplus, &hcross);

           (*hptilde)->data->data[idx + offset] += hplus;
           (*hctilde)->data->data[idx + offset] += hcross;
         }

         XLALDestroyREAL8Sequence(coarseFreqs);
         LALFree(valpha);
         LALFree(vepsilon);
         LALFree(vbetah);
         LALFree(cexp_i_alpha);
         LALFree(cexp_i_epsilon);
         LALFree(cexp_i_betah);
       }// End of Multibanding-specific.

     XLALDestroyCOMPLEX16FrequencySeries(htildelm);

   }//Loop over emmprime
 }//Loop over ell

    XLALDestroySphHarmFrequencySeries(*hlms);
    XLALFree(hlms);
 /*
      Loop over h+ and hx to rotate waveform by 2 \zeta.
      See discussion in Appendix C: Frame Transformation and Polarization Basis.
      The formula for \zeta is given by eq. C26.
*/
  if(fabs(pPrec->zeta_polarization) > 0)
  {
    COMPLEX16 PhPpolp, PhPpolc;
    REAL8 cosPolFac, sinPolFac;

    cosPolFac = cos(2.0 * pPrec->zeta_polarization);
    sinPolFac = sin(2.0 * pPrec->zeta_polarization);

    for (UINT4 i = offset; i < (*hptilde)->data->length; i++)
    {
        PhPpolp = (*hptilde)->data->data[i];
        PhPpolc = (*hctilde)->data->data[i];

        (*hptilde)->data->data[i] = cosPolFac * PhPpolp + sinPolFac * PhPpolc;
        (*hctilde)->data->data[i] = cosPolFac * PhPpolc - sinPolFac * PhPpolp;
    }
  }

  /* Free memory */
  XLALDestroyCOMPLEX16FrequencySeries(htilde22);
  XLALDestroyValue(ModeArray);
  XLALDestroyREAL8Sequence(freqs);

  #if DEBUG == 1
  printf("\n******Leaving IMRPhenomXPHM_hplushcross*****\n");
  #endif

  return XLAL_SUCCESS;
}



/*
  Core function of XLALSimIMRPhenomXPHMFromModes and XLALSimIMRPhenomXPHMFrequencySequence.
  Returns hptilde, hctilde for positive frequencies.
  The default non-precessing modes twisted up are 2|2|, 2|1|, 3|3|, 3|2| and 4|4|.
  It returns also the contribution of the corresponding negatives modes.
  It returns the same result than IMRPhenomXPHM_hplushcross but here it calls the individual precessing modes
  and then sum them all. It is therefor slower and does not include Multibanding for the angles.
  It can be evaulated in a non-uniform frequency grid. Assume positive frequencies.
*/
static int IMRPhenomXPHM_hplushcross_from_modes(
  COMPLEX16FrequencySeries **hptilde,  /**< [out] Frequency domain h+ GW strain */
  COMPLEX16FrequencySeries **hctilde,  /**< [out] Frequency domain hx GW strain */
  REAL8Sequence *freqs_In,             /**< Frequency array to evaluate the model. (fmin, fmax) for equally spaced grids. */
  IMRPhenomXWaveformStruct *pWF,       /**< IMRPhenomX Waveform Struct  */
  IMRPhenomXPrecessionStruct *pPrec,   /**< IMRPhenomXP Precession Struct  */
  LALDict *lalParams                   /**< LAL Dictionary Structure    */
)
{
  /* Set LIGOTimeGPS */
  LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0,0}

  lalParams = IMRPhenomXPHM_setup_mode_array(lalParams);
  LALValue *ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(lalParams);

  /* At this point ModeArray should contain the list of modes
  and therefore if NULL then something is wrong and abort. */
  if (ModeArray == NULL)
  {
    XLAL_ERROR(XLAL_EDOM, "ModeArray is NULL when it shouldn't be. Aborting.\n");
  }

  INT4 status = 0; //Variable to check correc functions calls


  /* Build the frequency array and initialize hctilde to the length of freqs. */
  REAL8Sequence *freqs;
  UINT4 offset = SetupWFArrays(&freqs, hptilde, freqs_In, pWF, ligotimegps_zero);

  /* Initialize hctilde according to hptilde. */
  size_t npts = (*hptilde)->data->length;
  *hctilde = XLALCreateCOMPLEX16FrequencySeries("hctilde: FD waveform", &ligotimegps_zero, (*hptilde)->f0, pWF->deltaF, &lalStrainUnit, npts);
  XLAL_CHECK (*hctilde, XLAL_ENOMEM, "Failed to allocated waveform COMPLEX16FrequencySeries of length %zu.", npts);
  memset((*hctilde)->data->data, 0, npts * sizeof(COMPLEX16));
  XLALUnitMultiply(&((*hctilde)->sampleUnits), &((*hctilde)->sampleUnits), &lalSecondUnit);

  /* Initialize useful powers of pi for the higher modes internal code. */
  status = IMRPhenomX_Initialize_Powers(&powers_of_lalpiHM, LAL_PI);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initialize useful powers of LAL_PI.");


  /* Loop over precessing modes */
  for(UINT4 ell = 2; ell <= 4; ell++)
  {
    for(INT4 emm = -1*ell; emm <= (INT4)ell; emm++)
    {
      #if DEBUG == 1
      printf("\n*****************************************************************\n");
      printf("                     Precessing mode (%i%i)                          ", ell, emm);
      printf("*******************************************************************\n");
      #endif

      COMPLEX16FrequencySeries *hlmpos = NULL;
      COMPLEX16FrequencySeries *hlmneg = NULL;

      /* We now call one single precessing mode.  */
      status = IMRPhenomXPHM_OneMode(&hlmpos, &hlmneg, freqs, pWF, pPrec, ell, emm, lalParams);
      XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "IMRPhenomXPHM_OneMode failed to generate IMRPhenomXHM waveform.");

      if (!(hlmpos)){ XLAL_ERROR(XLAL_EFUNC);}
      if (!(hlmneg)){ XLAL_ERROR(XLAL_EFUNC);}


      /* hlmpos and hlmneg are recomputed every time in the loop. Check that they always come out with the same length. */
      XLAL_CHECK ( ((*hptilde)->data->length == hlmpos->data->length) && ((*hptilde)->data->length  == hlmpos->data->length)
                   && (hlmpos->data->length  == hlmneg->data->length),
                   XLAL_EBADLEN,
                   "Inconsistent lengths between frequency series hlmpos (%d), hlmneg (%d), hptilde (%d) and hctilde (%d).",
                   hlmpos->data->length, hlmneg->data->length, (*hptilde)->data->length, (*hctilde)->data->length
                 );

      /*
         The precessing modes hlm that we just computed are in the J-frame.
         For computing hplus and hcross we have to sum them all with Ylm(thetaJN, 0) since the J-frame is aligned
         such that the line of sight N is in the x-z plane of the J-frame.
         See appendix C of Precessing Paper, in particular eq. C8.
      */
      COMPLEX16 Ylm = XLALSpinWeightedSphericalHarmonic(pPrec->thetaJN, 0, -2, ell, emm);
      COMPLEX16 Ylmstar = conj(Ylm);

      for(UINT4 i = offset; i < (hlmpos)->data->length; i++)
      {
        (*hptilde)->data->data[i] +=   0.5*(hlmpos->data->data[i] * Ylm + conj(hlmneg->data->data[i]) * Ylmstar);
        (*hctilde)->data->data[i] += I*0.5*(hlmpos->data->data[i] * Ylm - conj(hlmneg->data->data[i]) * Ylmstar);
      }

      XLALDestroyCOMPLEX16FrequencySeries(hlmpos);
      XLALDestroyCOMPLEX16FrequencySeries(hlmneg);
    }
  }// End loop over precessing modes

  /*
       Loop over h+ and hx to rotate waveform by 2 \zeta.
       See discussion in Appendix C: Frame Transformation and Polarization Basis.
       The formula for \zeta is given by eq. C24.
 */
  if(fabs(pPrec->zeta_polarization) > 0)
  {
    COMPLEX16 PhPpolp, PhPpolc;
    REAL8 cosPolFac, sinPolFac;

    cosPolFac = cos(2.0 * pPrec->zeta_polarization);
    sinPolFac = sin(2.0 * pPrec->zeta_polarization);

    for (UINT4 i = offset; i < (*hptilde)->data->length; i++)
    {
        PhPpolp = (*hptilde)->data->data[i];
        PhPpolc = (*hctilde)->data->data[i];

        (*hptilde)->data->data[i] = cosPolFac * PhPpolp + sinPolFac * PhPpolc;
        (*hctilde)->data->data[i] = cosPolFac * PhPpolc - sinPolFac * PhPpolp;
    }
  }

  /* Free memory */
  XLALDestroyValue(ModeArray);
  XLALDestroyREAL8Sequence(freqs);

  #if DEBUG == 1
  printf("\n******Leaving IMRPhenomXPHM_hplushcross_from_modes*****\n");
  #endif

  return XLAL_SUCCESS;
}


/*
  Core twisting up routine to get hptilde and hctilde.
  Twist one h_lmprime waveform in the precessing L-frame to the inertial J-frame for one frequency point
  as described in section III of Precessing paper.
  The explicit formula used implemented in this function correspond to eqs. E18, E19 in Precessing paper.
  This function is used inside a loop over frequencies inside a loop over mprime >0 up to l.
*/
static int IMRPhenomXPHMTwistUp(
  const REAL8 Mf,                          /**< Frequency (Hz) */
  const COMPLEX16 hlmprime,                /**< Underlying aligned-spin IMRPhenomXHM waveform. The loop is with mprime positive, but the mode has to be the negative one for positive frequencies.*/
  IMRPhenomXWaveformStruct *pWF,           /**< IMRPhenomX Waveform Struct */
  IMRPhenomXPrecessionStruct *pPrec,       /**< IMRPhenomXP Precession Struct */
  INT4  l,                                 /**< First index of the non-precessing (l,mprime) mode */
  INT4  mprime,                            /**< Second index of the non-precessing (l,mprime) mode */
  COMPLEX16 *hp,                           /**< [out] h_+ polarization \f$\tilde h_+\f$ */
  COMPLEX16 *hc                            /**< [out] h_x polarization \f$\tilde h_x\f$ */
)
{
  XLAL_CHECK(hp  != NULL, XLAL_EFAULT);
  XLAL_CHECK(hc  != NULL, XLAL_EFAULT);

  /* Euler angles */
  double alpha       = 0.0;
  double epsilon     = 0.0;

  double cBetah      = 0.0;
  double sBetah      = 0.0;

  COMPLEX16 cexp_i_alpha, cexp_i_epsilon = 1;

  if(pPrec->MBandPrecVersion == 0) /* No multibanding for angles */
  {
    switch(pPrec->IMRPhenomXPrecVersion)
    {
      case 101:    /* Post-Newtonian Euler angles. Single spin approximantion. See sections IV-B and IV-C in Precessing paper. */
      case 102:    /* The different number 10i means different PN order. */
      case 103:
      case 104:
      {
        Get_alpha_beta_epsilon(&alpha, &cBetah, &sBetah, &epsilon, mprime, Mf, pPrec, pWF);
        break;
      }
      case 220:    /* Use MSA angles. See section IV-D in Precessing paper. */
      case 221:
      case 222:
      case 223:
      case 224:
      {
        /* Get Euler angles. */
        const double v        = cbrt (LAL_PI * Mf * (2.0 / mprime) );
        const vector vangles  = IMRPhenomX_Return_phi_zeta_costhetaL_MSA(v,pWF,pPrec);
        double cos_beta       = 0.0;

        /* Get the offset for the Euler angles alpha and epsilon. */
        REAL8 alpha_offset_mprime = 0, epsilon_offset_mprime = 0;
        Get_alpha_epsilon_offset(&alpha_offset_mprime, &epsilon_offset_mprime, mprime, pPrec);

        alpha       = vangles.x - alpha_offset_mprime;
        epsilon     = vangles.y - epsilon_offset_mprime;
        cos_beta    = vangles.z;
        
        INT4 status = 0;
        status = IMRPhenomXWignerdCoefficients_cosbeta(&cBetah, &sBetah, cos_beta);
        XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "Call to IMRPhenomXWignerdCoefficients_cosbeta failed.");

        break;
      }
      default:
      {
        XLAL_ERROR(XLAL_EINVAL,"Error. IMRPhenomXPrecVersion not recognized. Recommended default is 223.\n");
        break;
      }
    }

    cexp_i_alpha   = cexp(+I*alpha);

    #if DEBUG == 1
    // Used for writing angles to file in debug mode.
    cexp_i_epsilon = cexp(+I*epsilon);
    pPrec->cexp_i_betah = cBetah + I*sBetah;
    #endif

   } // End of no multibanding
   else{ /*  For Multibanding  */
     cexp_i_alpha   = pPrec->cexp_i_alpha;
     cexp_i_epsilon = pPrec->cexp_i_epsilon;
     cBetah = (pPrec->cexp_i_betah + 1./pPrec->cexp_i_betah)*0.5;
     sBetah = (pPrec->cexp_i_betah - 1./pPrec->cexp_i_betah)*0.5/I;
   } // End of Multibanding-specific

  /* Useful powers of the Wigner coefficients */
  REAL8 cBetah2 = cBetah * cBetah;
  REAL8 cBetah3 = cBetah * cBetah2;
  REAL8 cBetah4 = cBetah * cBetah3;
  REAL8 cBetah5 = cBetah * cBetah4;
  REAL8 cBetah6 = cBetah * cBetah5;
  REAL8 cBetah7 = cBetah * cBetah6;
  REAL8 cBetah8 = cBetah * cBetah7;

  REAL8 sBetah2 = sBetah * sBetah;
  REAL8 sBetah3 = sBetah * sBetah2;
  REAL8 sBetah4 = sBetah * sBetah3;
  REAL8 sBetah5 = sBetah * sBetah4;
  REAL8 sBetah6 = sBetah * sBetah5;
  REAL8 sBetah7 = sBetah * sBetah6;
  REAL8 sBetah8 = sBetah * sBetah7;

  /*
      The following expressions for the Wigner-d coefficients correspond to those in appendix A of the Precessing paper.
      They are the same expressions used in IMRPhenomXPHMTwistUpOneMode.
  */

  COMPLEX16 hp_sum  = 0;
  COMPLEX16 hc_sum  = 0;

  /* Sum over l = 2 modes */
  if (l == 2 && mprime == 2){
    /* Sum up contributions to \tilde{h}_+ and \tilde{h}_x */
    /* Precompute powers of e^{i m alpha} */
    COMPLEX16 cexp_2i_alpha     = cexp_i_alpha  * cexp_i_alpha;

    COMPLEX16 cexp_mi_alpha     = 1.0 / cexp_i_alpha;
    COMPLEX16 cexp_m2i_alpha    = cexp_mi_alpha * cexp_mi_alpha;

    COMPLEX16 cexp_im_alpha_l2[5]  = {cexp_m2i_alpha, cexp_mi_alpha, 1.0, cexp_i_alpha, cexp_2i_alpha};

    COMPLEX16 Y2mA[5] = {pPrec->Y2m2, pPrec->Y2m1, pPrec->Y20, pPrec->Y21, pPrec->Y22};

    //                  d^2_{-2,2}    d^2_{-1,2}             d^2_{0,2}                   d^2_{1,2}       d^2_{2,2}
    COMPLEX16 d22[5]   = {sBetah4, 2.0*cBetah*sBetah3, pPrec->sqrt6*sBetah2*cBetah2, 2.0*cBetah3*sBetah, cBetah4};
    //                  d^2_{-2,-2}  d^2_{-1,-2}  d^2_{0,-2}  d^2_{1,-2}  d^2_{2,-2}
    COMPLEX16 d2m2[5]  = {d22[4],    -d22[3],      d22[2],     -d22[1],     d22[0]}; /* Exploit symmetry d^2_{-m,-2} = (-1)^m d^2_{-m,2}. See eq. A2 of Precessing paper */

    for(int m=-2; m<=2; m++)
    {
      /* Transfer functions, see eqs. 3.5-3.7 in Precessing paper */
      COMPLEX16 A2m2emm  = cexp_im_alpha_l2[-m+2] * d2m2[m+2]  * Y2mA[m+2];
      COMPLEX16 A22emmstar = cexp_im_alpha_l2[m+2] * d22[m+2] *  conj(Y2mA[m+2]);
      hp_sum +=    A2m2emm + A22emmstar;
      hc_sum += I*(A2m2emm - A22emmstar);
    }
  }

  if (l == 2 && mprime == 1){
    /* Sum up contributions to \tilde{h}_+ and \tilde{h}_x */
    /* Precompute powers of e^{i m alpha} */
    COMPLEX16 cexp_2i_alpha     = cexp_i_alpha  * cexp_i_alpha;

    COMPLEX16 cexp_mi_alpha     = 1.0 / cexp_i_alpha;
    COMPLEX16 cexp_m2i_alpha    = cexp_mi_alpha * cexp_mi_alpha;

    COMPLEX16 cexp_im_alpha_l2[5]  = {cexp_m2i_alpha, cexp_mi_alpha, 1.0, cexp_i_alpha, cexp_2i_alpha};

    COMPLEX16 Y2mA[5] = {pPrec->Y2m2, pPrec->Y2m1, pPrec->Y20, pPrec->Y21, pPrec->Y22};

    //                    d^2_{-2,1}          d^2_{-1,1}                       d^2_{0,1}                                       d^2_{1,1}                        d^2_{2,1}
    COMPLEX16 d21[5]   = {2.0*cBetah*sBetah3, 3.0*cBetah2*sBetah2 - sBetah4, pPrec->sqrt6*(cBetah3*sBetah - cBetah*sBetah3), cBetah2*(cBetah2 - 3.0*sBetah2), -2.0*cBetah3*sBetah};
    //                  d^2_{-2,-1}  d^2_{-1,-1}  d^2_{0,-1}  d^2_{1,-1}  d^2_{2,-1}
    COMPLEX16 d2m1[5]  = {-d21[4],   d21[3],     -d21[2],    d21[1],     -d21[0]}; /* Exploit symmetry d^2_{-m,-1} = -(-1)^m d^2_{m,1}.  See eq. A2 of Precessing paper.  */


    for(int m=-2; m<=2; m++)
    {
      /* Transfer functions, see eqs. 3.5-3.7 in Precessing paper. */
      COMPLEX16 A2m1emm  = cexp_im_alpha_l2[-m+2]  * d2m1[m+2]  * Y2mA[m+2];
      COMPLEX16 A21emmstar = cexp_im_alpha_l2[m+2] * d21[m+2] *  conj(Y2mA[m+2]);
      hp_sum +=    A2m1emm + A21emmstar;
      hc_sum += I*(A2m1emm - A21emmstar);
     }

  }

  /* Sum over l = 3 modes */
  if (l == 3 && mprime == 3){
    /* Sum up contributions to \tilde{h}_+ and \tilde{h}_x */
    /* Precompute powers of e^{i m alpha} */
    COMPLEX16 cexp_2i_alpha     = cexp_i_alpha  * cexp_i_alpha;
    COMPLEX16 cexp_3i_alpha     = cexp_i_alpha  * cexp_2i_alpha;

    COMPLEX16 cexp_mi_alpha     = 1.0 / cexp_i_alpha;
    COMPLEX16 cexp_m2i_alpha    = cexp_mi_alpha * cexp_mi_alpha;
    COMPLEX16 cexp_m3i_alpha    = cexp_mi_alpha * cexp_m2i_alpha;

    COMPLEX16 cexp_im_alpha_l3[7]  = {cexp_m3i_alpha, cexp_m2i_alpha, cexp_mi_alpha, 1.0, cexp_i_alpha, cexp_2i_alpha, cexp_3i_alpha};

    COMPLEX16 Y3mA[7] = {pPrec->Y3m3, pPrec->Y3m2, pPrec->Y3m1, pPrec->Y30, pPrec->Y31, pPrec->Y32, pPrec->Y33};

    //                  d^3_{-3,3}    d^3_{-2,3}                 d^3_{-1,3}                    d^3_{0,3}                         d^3_{1,3}                      d^3_{2,3}                   d^3_{3,3}
    COMPLEX16 d33[7]   = {sBetah6, pPrec->sqrt6*cBetah*sBetah5, pPrec->sqrt15*cBetah2*sBetah4, 2.0*pPrec->sqrt5*cBetah3*sBetah3, pPrec->sqrt15*cBetah4*sBetah2, pPrec->sqrt6*cBetah5*sBetah, cBetah6};
    //                  d^3_{-3,-3}  d^3_{-2,-3}  d^3_{-1,-3}  d^3_{0,-3}  d^3_{1,-3}  d^3_{2,-3}  d^3_{3,-3}
    COMPLEX16 d3m3[7]  = {d33[6],    -d33[5],     d33[4],      -d33[3],    d33[2],     -d33[1],    d33[0]}; /* Exploit symmetry d^3_{-m,-3} = -(-1)^m d^3_{m,3}. See eq. A2 of Precessing paper. */

    for(int m=-3; m<=3; m++)
    {
      /* Transfer functions */
      COMPLEX16 A3m3emm  = cexp_im_alpha_l3[-m+3]  * d3m3[m+3]  * Y3mA[m+3];
      COMPLEX16 A33emmstar   = cexp_im_alpha_l3[m+3] * d33[m+3] *  conj(Y3mA[m+3]);
      hp_sum +=    A3m3emm - A33emmstar;
      hc_sum += I*(A3m3emm + A33emmstar);
    }
  }

  /* Sum over l = 3 modes */
  if (l == 3 && mprime == 2){
    /* Sum up contributions to \tilde{h}_+ and \tilde{h}_x */
    /* Precompute powers of e^{i m alpha} */
    COMPLEX16 cexp_2i_alpha     = cexp_i_alpha  * cexp_i_alpha;
    COMPLEX16 cexp_3i_alpha     = cexp_i_alpha  * cexp_2i_alpha;

    COMPLEX16 cexp_mi_alpha     = 1.0 / cexp_i_alpha;
    COMPLEX16 cexp_m2i_alpha    = cexp_mi_alpha * cexp_mi_alpha;
    COMPLEX16 cexp_m3i_alpha    = cexp_mi_alpha * cexp_m2i_alpha;

    COMPLEX16 cexp_im_alpha_l3[7]  = {cexp_m3i_alpha, cexp_m2i_alpha, cexp_mi_alpha, 1.0, cexp_i_alpha, cexp_2i_alpha, cexp_3i_alpha};

    COMPLEX16 Y3mA[7] = {pPrec->Y3m3, pPrec->Y3m2, pPrec->Y3m1, pPrec->Y30, pPrec->Y31, pPrec->Y32, pPrec->Y33};

    //                    d^3_{-3,2}                     d^3_{-2,2}                     d^3_{-1,21}                                            d^3_{0,2}                                         d^3_{1,2}                                              d^3_{2,2}                        d^3_{3,2}
    COMPLEX16 d32[7]   = {pPrec->sqrt6*cBetah*sBetah5, sBetah4*(5.0*cBetah2 - sBetah2), pPrec->sqrt10*sBetah3*(2.0*cBetah3 - cBetah*sBetah2), pPrec->sqrt30*cBetah2*(cBetah2 - sBetah2)*sBetah2, pPrec->sqrt10*cBetah3*(cBetah2*sBetah - 2.0*sBetah3), cBetah4*(cBetah2 - 5.0*sBetah2), -1.0*pPrec->sqrt6*cBetah5*sBetah};
    //                  d^3_{-3,-2}  d^3_{-2,-2}  d^3_{-1,-2}  d^3_{0,-2}  d^3_{1,-2}  d^3_{2,-2}  d^3_{3,-2}
    COMPLEX16 d3m2[7]  = {-d32[6],   d32[5],      -d32[4],     d32[3],     -d32[2],    d32[1],     -d32[0]}; /* Exploit symmetry d^3_{-m,-2} = (-1)^m d^3_{m,2}. See eq. A2 of Precessing paper.  */


    for(int m=-3; m<=3; m++)
    {
      /* Transfer functions, see eqs. 3.5-3.7 in Precessing paper. */
      COMPLEX16 A3m2emm  =  cexp_im_alpha_l3[-m+3]  * d3m2[m+3] * Y3mA[m+3];
      COMPLEX16 A32emmstar  =  cexp_im_alpha_l3[m+3] * d32[m+3] *  conj(Y3mA[m+3]);
      hp_sum +=    A3m2emm - A32emmstar;
      hc_sum += I*(A3m2emm + A32emmstar);
    }

  }

  /* Sum over l = 4 modes */
  if (l == 4 && mprime == 4){
    /* Sum up contributions to \tilde{h}_+ and \tilde{h}_x */
    /* Precompute powers of e^{i m alpha} */
    COMPLEX16 cexp_2i_alpha     = cexp_i_alpha  * cexp_i_alpha;
    COMPLEX16 cexp_3i_alpha     = cexp_i_alpha  * cexp_2i_alpha;
    COMPLEX16 cexp_4i_alpha     = cexp_i_alpha  * cexp_3i_alpha;

    COMPLEX16 cexp_mi_alpha     = 1.0 / cexp_i_alpha;
    COMPLEX16 cexp_m2i_alpha    = cexp_mi_alpha * cexp_mi_alpha;
    COMPLEX16 cexp_m3i_alpha    = cexp_mi_alpha * cexp_m2i_alpha;
    COMPLEX16 cexp_m4i_alpha    = cexp_mi_alpha * cexp_m3i_alpha;

    COMPLEX16 cexp_im_alpha_l4[9]  = {cexp_m4i_alpha, cexp_m3i_alpha, cexp_m2i_alpha, cexp_mi_alpha, 1.0, cexp_i_alpha, cexp_2i_alpha, cexp_3i_alpha, cexp_4i_alpha};

    COMPLEX16 Y4mA[9] = {pPrec->Y4m4, pPrec->Y4m3, pPrec->Y4m2, pPrec->Y4m1, pPrec->Y40, pPrec->Y41, pPrec->Y42, pPrec->Y43, pPrec->Y44};

    //                    d^4_{-4,4}         d^4_{-3,4}               d^4_{-2,4}                         d^4_{-1,4}                          d^4_{0,4}                       d^4_{1,4}                           d^4_{2,4}                          d^4_{3,4}                  d^4_{44}
    COMPLEX16 d44[9]   = {sBetah8, 2.0*pPrec->sqrt2*cBetah*sBetah7, 2.0*pPrec->sqrt7*cBetah2*sBetah6, 2.0*pPrec->sqrt14*cBetah3*sBetah5, pPrec->sqrt70*cBetah4*sBetah4, 2.0*pPrec->sqrt14*cBetah5*sBetah3, 2.0*pPrec->sqrt7*cBetah6*sBetah2, 2.0*pPrec->sqrt2*cBetah7*sBetah, cBetah8};
      //                  d^4_{4,-4}  d^4_{-3,-4}  d^4_{-2,-4}  d^4_{-1,-4}  d^4_{0,-4}  d^4_{1,-4}  d^4_{2,-4}  d^4_{3,-4}  d^4_{4-4}
    COMPLEX16 d4m4[9]  = {d44[8],     -d44[7],     d44[6],      -d44[5],     d44[4],    -d44[3],     d44[2],     -d44[1],    d44[0]}; /* Exploit symmetry d^4_{-m,-4} = (-1)^m d^4_{m,4}. See eq. A2 of Precessing paper.  */

    for(int m=-4; m<=4; m++)
    {
      /* Transfer functions, see eqs. 3.5-3.7 in Precessing paper. */
      COMPLEX16 A4m4emm  =  cexp_im_alpha_l4[-m+4]  * d4m4[m+4]  * Y4mA[m+4];
      COMPLEX16 A44emmstar  =  cexp_im_alpha_l4[m+4] * d44[m+4] *  conj(Y4mA[m+4]);
      hp_sum +=    A4m4emm + A44emmstar;
      hc_sum += I*(A4m4emm - A44emmstar);
    }
  }

  /* Sum over l = 4 modes. This is only used when twisting PhenomHM. */
  if (l == 4 && mprime == 3){
    /* Sum up contributions to \tilde{h}_+ and \tilde{h}_x */
    /* Precompute powers of e^{i m alpha} */
    COMPLEX16 cexp_2i_alpha     = cexp_i_alpha  * cexp_i_alpha;
    COMPLEX16 cexp_3i_alpha     = cexp_i_alpha  * cexp_2i_alpha;
    COMPLEX16 cexp_4i_alpha     = cexp_i_alpha  * cexp_3i_alpha;

    COMPLEX16 cexp_mi_alpha     = 1.0 / cexp_i_alpha;
    COMPLEX16 cexp_m2i_alpha    = cexp_mi_alpha * cexp_mi_alpha;
    COMPLEX16 cexp_m3i_alpha    = cexp_mi_alpha * cexp_m2i_alpha;
    COMPLEX16 cexp_m4i_alpha    = cexp_mi_alpha * cexp_m3i_alpha;

    COMPLEX16 cexp_im_alpha_l4[9]  = {cexp_m4i_alpha, cexp_m3i_alpha, cexp_m2i_alpha, cexp_mi_alpha, 1.0, cexp_i_alpha, cexp_2i_alpha, cexp_3i_alpha, cexp_4i_alpha};

    COMPLEX16 Y4mA[9] = {pPrec->Y4m4, pPrec->Y4m3, pPrec->Y4m2, pPrec->Y4m1, pPrec->Y40, pPrec->Y41, pPrec->Y42, pPrec->Y43, pPrec->Y44};

    //                    d^4_{-4,3}         d^4_{-3,3}               d^4_{-2,3}                         d^4_{-1,3}                          d^4_{0,3}                       d^4_{1,3}                           d^4_{2,3}                          d^4_{3,3}                  d^4_{43}
    COMPLEX16 d43[9]   = {2*pPrec->sqrt2*cBetah*sBetah7, 7*cBetah2*sBetah6-sBetah8, pPrec->sqrt14*(3*cBetah3*sBetah5-cBetah*sBetah7), pPrec->sqrt7*(5*cBetah4*sBetah4-3*cBetah2*sBetah6), 2*5.916079783099616*(cBetah5*sBetah3-cBetah3*sBetah5), pPrec->sqrt7*(3*cBetah6*sBetah2-5*cBetah4*sBetah4), pPrec->sqrt14*(cBetah7*sBetah-3*cBetah5*sBetah3), cBetah8-7*cBetah6*sBetah2, -2.*pPrec->sqrt2*cBetah7*sBetah};
      //                  d^4_{4,-3}  d^4_{-3,-3}  d^4_{-2,-3}  d^4_{-1,-3}  d^4_{0,-3}  d^4_{1,-3}  d^4_{2,-3}  d^4_{3,-3}  d^4_{4-3}
    COMPLEX16 d4m3[9]  = {-d43[8],     d43[7],     -d43[6],      d43[5],     -d43[4],    d43[3],     -d43[2],     d43[1],    -d43[0]}; /* Exploit symmetry d^4_{-m,-3} = -(-1)^m d^4_{m,3}. See eq. A2 of Precessing paper.  */

    for(int m=-4; m<=4; m++)
    {
      /* Transfer functions, see eqs. 3.5-3.7 in Precessing paper. */
      COMPLEX16 A4m3emm  =  cexp_im_alpha_l4[-m+4]  * d4m3[m+4]  * Y4mA[m+4];
      COMPLEX16 A43emmstar  =  cexp_im_alpha_l4[m+4] * d43[m+4] *  conj(Y4mA[m+4]);
      hp_sum +=    A4m3emm + A43emmstar;
      hc_sum += I*(A4m3emm - A43emmstar);
    }
  }

  COMPLEX16 eps_phase_hP_lmprime;

  if(pPrec->MBandPrecVersion == 0) // No multibanding
  {
     eps_phase_hP_lmprime = cexp(-1.*mprime*I*epsilon) * hlmprime / 2.0;
  }
  else{          // With multibanding
    COMPLEX16 exp_imprime_epsilon = cexp_i_epsilon;
    for(INT4 i=1; i<mprime; i++)
    {
      exp_imprime_epsilon *= cexp_i_epsilon;
    }
    eps_phase_hP_lmprime = 1./exp_imprime_epsilon * hlmprime / 2.0;
  }


  /* Return h_+ and h_x */
  *hp += eps_phase_hP_lmprime * hp_sum;
  *hc += eps_phase_hP_lmprime * hc_sum;

  #if DEBUG == 1
  /* Save angles in output file.  */
  FILE *fileangle;
  char fileSpec[40];

  if(pPrec->MBandPrecVersion == 0)
  {
      sprintf(fileSpec, "angles_hphc_%i%i.dat", l, mprime);
  }
  else
  {
      sprintf(fileSpec, "angles_hphc_MB_%i%i.dat", l, mprime);
  }
  fileangle = fopen(fileSpec,"a");

  fprintf(fileangle, "%.16e  %.16e  %.16e  %.16e  %.16e  %.16e  %.16e\n",  XLALSimIMRPhenomXUtilsMftoHz(Mf, pWF->Mtot), creal(cexp_i_alpha), cimag(cexp_i_alpha), creal(cexp_i_epsilon), cimag(cexp_i_epsilon), creal(pPrec->cexp_i_betah), cimag(pPrec->cexp_i_betah));
  fclose(fileangle);
  #endif

  return XLAL_SUCCESS;
}

/* @} */
/* @} */

/** @addtogroup LALSimIMRPhenomX_c
* @{
* @name Routines for IMRPhenomXPHM
* @{
* @author Cecilio García Quirós, Geraint Pratten
*
* @brief C code for IMRPhenomXPHM phenomenological waveform model.
*
* This is a frequency domain precessing model based on the twisting-up of the algined spin model with higher modes IMRPhenomXHM.
* See G.Pratten et al for details. Any studies that use this waveform model should include
* a reference to this paper.
*
* @note DCC link to the paper: https://dcc.ligo.org/LIGO-P2000039. This paper will be refered in the code as the "Precessing paper".
*
* Waveform flags:
*
* All the flags for IMRPhenomXP apply here plus the following ones:
*
*   TwistPhenomHM: option to twist-up the AS model PhenomHM instead of PhenomXHM. It is only available for the polarizations, not for individual modes.
*       - 0: (DEFAULT) twist-up PhenomXHM
*       - 1: twist-up PhenomHM
*
*   UseModes: Determine how the polarizations hp, hc are computed.
*       - 0: (DEFAULT) Compute the non-precessing modes once and do the twistin up as in eq. 3.5-3.7 in the Precessing paper.
*       - 1: Compute first the individual precessing modes in the inertial J-frame and sum them to get the polarizations.
*
*   ModesL0Frame: Determine in which frame the individual precessing modes are returned.
*       - 0: inertial J-frame (DEFAULT).
*       - 1: inertial L0-frame (only working near the aligned spin limit).
*
*   PrecModes: Determine which indiviual modes are returned, the non-precessing or the precessing.
*       - 0: (DEFAULT) Return the precessing individual mode in the J-frame.
*       - 1: Return the non-precessing individual mode before the twisting-up with the modified final spin.
*
* Multibanding flags:
*
*   PrecThresholdMband: Determines the accuracy and speed of the Multibanding algorithm for the Euler angles. The higher the threshold the faster is the algorithm but also less accurate.
*        - 0.001 (DEFAULT)
*        - 0: Switch off the multibanding.
*
*   MBandPrecVersion: Determines the algorithm to build the non-uniform frequency grid for the Euler angles.
*        - 0: (DEFAULT) Not use multibanding.  Activated to 1 when PrecThresholdMband is non-zero.
*        - 1: Use the same grid that for the non-precessing modes. Activated when PrecThresholdMband is non-zero.
**/



/*********************************************/
/*                                           */
/*      SINGLE-MODE PRECESSING FUNCTIONS     */
/*                                           */
/*********************************************/

/**
    Function to compute one hlm precessing mode in an uniform frequency grid.
    By default the mode is given in the inertial J-frame. It can be transform to the L0-frame with the option "L0Frame" which
    currently only works for cases near AS limit.
    Returns two frequency series, one for the positive frequencies and other for the negative frequencies since, as opposite to the
    aligned spin case, in the precessing case all the modes have support in the whole frequency regime.
    This is a wrapper of the internal core function that actually does the calculation IMRPhenomXPHM_OneMode.
*/
int XLALSimIMRPhenomXPHMOneMode(
  COMPLEX16FrequencySeries **hlmpos,      /**< [out] Frequency-domain waveform hlm inertial frame positive frequencies */
  COMPLEX16FrequencySeries **hlmneg,      /**< [out] Frequency-domain waveform hlm inertial frame negative frequencies */
  const UINT4 l,                          /**< First index of the (l,m) precessing mode */
  const INT4  m,                          /**< Second index of the (l,m) precessing mode */
  REAL8 m1_SI,                            /**< mass of companion 1 (kg) */
  REAL8 m2_SI,                            /**< mass of companion 2 (kg) */
  REAL8 chi1x,                            /**< x-component of the dimensionless spin of object 1  w.r.t. Lhat = (0,0,1) */
  REAL8 chi1y,                            /**< y-component of the dimensionless spin of object 1  w.r.t. Lhat = (0,0,1) */
  REAL8 chi1z,                            /**< z-component of the dimensionless spin of object 1  w.r.t. Lhat = (0,0,1) */
  REAL8 chi2x,                            /**< x-component of the dimensionless spin of object 2  w.r.t. Lhat = (0,0,1) */
  REAL8 chi2y,                            /**< y-component of the dimensionless spin of object 2  w.r.t. Lhat = (0,0,1) */
  REAL8 chi2z,                            /**< z-component of the dimensionless spin of object 2  w.r.t. Lhat = (0,0,1) */
  const REAL8 distance,                   /**< distance of source (m) */
  const REAL8 inclination,                /**< inclination of source (rad) */
  const REAL8 phiRef,                     /**< reference orbital phase (rad) */
  const REAL8 deltaF,                     /**< Sampling frequency (Hz) */
  const REAL8 f_min,                      /**< Starting GW frequency (Hz) */
  const REAL8 f_max,                      /**< End frequency; 0 defaults to ringdown cutoff freq */
  const REAL8 fRef_In,                    /**< Reference frequency */
  LALDict *lalParams                      /**<LAL Dictionary */
)
{
  /* Variable to check correct calls to functions. */
  INT4 status;

  /* Check if m1 > m2, swap the bodies otherwise. */
  status = XLALIMRPhenomXPCheckMassesAndSpins(&m1_SI,&m2_SI,&chi1x,&chi1y,&chi1z,&chi2x,&chi2y,&chi2z);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: XLALIMRPhenomXPCheckMassesAndSpins failed.\n");

  #if DEBUG == 1
  printf("fRef_In : %e\n",fRef_In);
  printf("m1_SI   : %e\n",m1_SI);
  printf("m2_SI   : %e\n",m2_SI);
  printf("chi1_l   : %e\n",chi1z);
  printf("chi2_l   : %e\n",chi2z);
  printf("phiRef  : %e\n",phiRef);
  printf("Prec V. : %d\n\n",XLALSimInspiralWaveformParamsLookupPhenomXPrecVersion(lalParams));
  printf("Performing sanity checks...\n");
  #endif

  /* Perform initial sanity checks */
  XLAL_CHECK(NULL != hlmpos, XLAL_EFAULT,  "Error: hlmpos already defined.                         \n");
  XLAL_CHECK(NULL != hlmneg, XLAL_EFAULT,  "Error: hlmneg already defined.                         \n");
  XLAL_CHECK(fRef_In  >= 0, XLAL_EFUNC,    "Error: fRef_In must be positive or set to 0 to ignore. \n");
  XLAL_CHECK(deltaF   >  0, XLAL_EFUNC,    "Error: deltaF must be positive and greater than 0.     \n");
  XLAL_CHECK(m1_SI    >  0, XLAL_EFUNC,    "Error: m1 must be positive and greater than 0.         \n");
  XLAL_CHECK(m2_SI    >  0, XLAL_EFUNC,    "Error: m2 must be positive and greater than 0.         \n");
  XLAL_CHECK(f_min    >  0, XLAL_EFUNC,    "Error: f_min must be positive and greater than 0.      \n");
  XLAL_CHECK(f_max    >= 0, XLAL_EFUNC,    "Error: f_max must be non-negative.                     \n");
  XLAL_CHECK(distance >  0, XLAL_EFUNC,    "Error: Distance must be positive and greater than 0.   \n");

  /*
  Perform a basic sanity check on the region of the parameter space in which model is evaluated.
  Behaviour is as follows, consistent with the choices for IMRPhenomXAS/IMRPhenomXHM
    - For mass ratios <= 20.0 and spins <= 0.99: no warning messages.
    - For 1000 > mass ratio > 20 and spins <= 0.99: print a warning message that we are extrapolating outside of *NR* calibration domain.
    - For mass ratios > 1000: throw a hard error that model is not valid.
    - For spins > 0.99: throw a warning that we are extrapolating the model to extremal

  */
  REAL8 mass_ratio;
  if(m1_SI > m2_SI)
  {
    mass_ratio = m1_SI / m2_SI;
  }
  else
  {
    mass_ratio = m2_SI / m1_SI;
  }
  if(mass_ratio > 20.0  ) { XLAL_PRINT_WARNING("Warning: Extrapolating outside of Numerical Relativity calibration domain. NNLO angles may become pathological at large mass ratios.\n"); }
  if(mass_ratio > 1000. && fabs(mass_ratio - 1000) > 1e-12) { XLAL_ERROR(XLAL_EDOM, "ERROR: Model not valid at mass ratios beyond 1000.\n"); } // The 1e-12 is to avoid rounding errors
  if(fabs(chi1z) > 0.99 || fabs(chi2z) > 0.99) { XLAL_PRINT_WARNING("Warning: Extrapolating to extremal spins, model is not trusted.\n"); }

  /* If no reference frequency is given, set it to the starting gravitational wave frequency. */
  REAL8 fRef = (fRef_In == 0.0) ? f_min : fRef_In;

  /* Use an auxiliar laldict to not overwrite the input argument */
  LALDict *lalParams_aux;
  /* setup mode array */
  if (lalParams == NULL)
  {
      lalParams_aux = XLALCreateDict();
  }
  else{
      lalParams_aux = XLALDictDuplicate(lalParams);
  }
  /* Check that the modes chosen are available for the model */
  XLAL_CHECK(check_input_mode_array(lalParams_aux) == XLAL_SUCCESS, XLAL_EFAULT, "Not available mode chosen.\n");

  #if DEBUG == 1
  printf("\n\n **** Initializing waveform struct... **** \n\n");
  #endif

  /* Initialize the useful powers of LAL_PI */
  status = IMRPhenomX_Initialize_Powers(&powers_of_lalpi, LAL_PI);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initialize useful powers of LAL_PI.");

  /* Initialize IMR PhenomX Waveform struct and check that it initialized correctly. */
  /* We pass inclination 0 since for the individual modes is not relevant. */
  IMRPhenomXWaveformStruct *pWF;
  pWF    = XLALMalloc(sizeof(IMRPhenomXWaveformStruct));
  status = IMRPhenomXSetWaveformVariables(pWF, m1_SI, m2_SI, chi1z, chi2z, deltaF, fRef, phiRef, f_min, f_max, distance, inclination, lalParams_aux, PHENOMXDEBUG);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetWaveformVariables failed.\n");

  /*
      Create a REAL8 frequency series.
      Use fLow, fHigh, deltaF to compute frequency sequence. Only pass the boundaries (fMin, fMax).
  */
  REAL8Sequence *freqs = XLALCreateREAL8Sequence(2);
  freqs->data[0] = pWF->fMin;
  freqs->data[1] = pWF->f_max_prime;


  #if DEBUG == 1
  printf("\n\n **** Initializing precession struct... **** \n\n");
  #endif

  /* Initialize IMRPhenomX Precession struct and check that it generated successfully. */
  IMRPhenomXPrecessionStruct *pPrec;
  pPrec  = XLALMalloc(sizeof(IMRPhenomXPrecessionStruct));

  status = IMRPhenomXGetAndSetPrecessionVariables(
             pWF,
             pPrec,
             m1_SI,
             m2_SI,
             chi1x,
             chi1y,
             chi1z,
             chi2x,
             chi2y,
             chi2z,
             lalParams_aux,
             PHENOMXDEBUG
           );
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetPrecessionVariables failed.\n");

  /* Ensure recovering AS limit when modes are in the L0 frame. */
  if(XLALSimInspiralWaveformParamsLookupPhenomXPHMModesL0Frame(lalParams_aux)==1)
  {
    XLAL_PRINT_WARNING("The L0Frame option only works near the AS limit, it should not be used otherwise.");
    switch(XLALSimInspiralWaveformParamsLookupPhenomXPConvention(lalParams_aux))
    {
      case 0:
      case 5:
        //pWF->phi0 = pPrec->phi0_aligned;
        break;
      case 1:
      case 6:
      case 7:
      {
        //pWF->phi0 = pPrec->epsilon0 - pPrec->alpha0 + phiRef;
        pWF->phi0 = phiRef;
        break;
      }
    }
  }


  #if DEBUG == 1
  printf("\n\n **** Calling IMRPhenomXPHM_OneMode... **** \n\n");
  #endif

  /* We now call the core IMRPhenomXPHM_OneMode waveform generator */
  status = IMRPhenomXPHM_OneMode(hlmpos, hlmneg, freqs, pWF, pPrec, l, m, lalParams_aux);
  XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "IMRPhenomXPHM_OneMode failed to generate IMRPhenomXHM waveform.");

  /* Tranform modes to L0-frame if requested. It only works for (near) AS cases. */
  if(XLALSimInspiralWaveformParamsLookupPhenomXPHMModesL0Frame(lalParams_aux)==1)
  {
    switch(XLALSimInspiralWaveformParamsLookupPhenomXPConvention(lalParams_aux))
    {
      case 0:
      case 5:
        //pWF->phi0 = pPrec->phi0_aligned;
        break;
      case 1:
      case 6:
      case 7:
      {
        COMPLEX16 shiftpos = cexp( abs(m)*I*( pPrec->epsilon0 - pPrec->alpha0) );
        COMPLEX16 shiftneg = 1./shiftpos;

        for(UINT4 i = 0; i<(*hlmpos)->data->length; i++)
        {
          (*hlmpos)->data->data[i] *= shiftpos;
          (*hlmneg)->data->data[i] *= shiftneg;
        }
        break;
      }
    }
  }

  #if DEBUG == 1
  printf("\n\n **** Call to IMRPhenomXPHM_OneMode complete. **** \n\n");
  #endif

  /* Resize hptilde, hctilde */
  REAL8 lastfreq;
  if (pWF->f_max_prime < pWF->fMax)
  {
    /* The user has requested a higher f_max than Mf = fCut.
    Resize the frequency series to fill with zeros beyond the cutoff frequency. */
    lastfreq = pWF->fMax;
    XLAL_PRINT_WARNING("The input f_max = %.2f Hz is larger than the internal cutoff of Mf=0.3 (%.2f Hz). Array will be filled with zeroes between these two frequencies.\n", pWF->fMax, pWF->f_max_prime);
  }
  else{  // We have to look for a power of 2 anyway.
    lastfreq = pWF->f_max_prime;
  }
  // We want to have the length be a power of 2 + 1
  size_t n_full = NextPow2(lastfreq / deltaF) + 1;
  size_t n = (*hlmpos)->data->length;

  /* Resize the COMPLEX16 frequency series */
  *hlmpos = XLALResizeCOMPLEX16FrequencySeries(*hlmpos, 0, n_full);
  XLAL_CHECK (*hlmpos, XLAL_ENOMEM, "Failed to resize hlmpos COMPLEX16FrequencySeries of length %zu (for internal fCut=%f) to new length %zu (for user-requested f_max=%f).", n, pWF->fCut, n_full, pWF->fMax );

  /* Resize the COMPLEX16 frequency series */
  *hlmneg = XLALResizeCOMPLEX16FrequencySeries(*hlmneg, 0, n_full);
  XLAL_CHECK (*hlmneg, XLAL_ENOMEM, "Failed to resize hlmneg COMPLEX16FrequencySeries of length %zu (for internal fCut=%f) to new length %zu (for user-requested f_max=%f).", n, pWF->fCut, n_full, pWF->fMax );


  /* Free memory */
  LALFree(pWF);
  LALFree(pPrec);
  XLALDestroyREAL8Sequence(freqs);
  XLALDestroyDict(lalParams_aux);

  return XLAL_SUCCESS;
}
/** @}
* @} **/


/**
    Core funciton to compute an individual precessing mode hlm in the inertial J-frame.
    Returns two frequency series, one for the positive frequencies and other for the negative frequencies.
    It can be evaluated in a non-uniform frequency grid through the argument REAL8Seuqnce *freqs_In. This is in fact done when Calling
    XLALSimIMRPhenomXPHMFrequencySequence with the option of 'UseModes' activated.
*/
static int IMRPhenomXPHM_OneMode(
  COMPLEX16FrequencySeries **hlmpos,    /**< [out] Frequency domain hlm GW strain inertial frame positive frequencies */
  COMPLEX16FrequencySeries **hlmneg,    /**< [out] Frequency domain hlm GW strain inertial frame negative frequencies */
  REAL8Sequence *freqs_In,              /**< Input frequency grid        */
  IMRPhenomXWaveformStruct *pWF,        /**< IMRPhenomX Waveform Struct  */
  IMRPhenomXPrecessionStruct *pPrec,    /**< IMRPhenomXP Precession Struct  */
  UINT4 ell,                            /**< l index of the (l,m) precessing mode */
  INT4  m,                              /**< m index of the (l,m) precessing mode */
  LALDict *lalParams                    /**< LAL Dictionary Structure    */
)
{
  /* Set LIGOTimeGPS */
  LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0,0}

  REAL8 deltaF = pWF->deltaF;

  lalParams = IMRPhenomXPHM_setup_mode_array(lalParams);
  LALValue *ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(lalParams);

  /* At this point ModeArray should contain the list of modes
  and therefore if NULL then something is wrong and abort. */
  if (ModeArray == NULL)
  {
    XLAL_ERROR(XLAL_EDOM, "ModeArray is NULL when it shouldn't be. Aborting.\n");
  }
  
  /* Check that the co-precessing ModeArray has at least one ell mode. If not, twisting-up is not possible. */
  bool mode_arrays_consistent = false;
  INT4 emm = -(INT4)ell;
  while (mode_arrays_consistent == false && emm<=(INT4)ell){
    if (XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, emm) ==1){
      mode_arrays_consistent = true;
    }
    emm++;
  }
  if(mode_arrays_consistent == false){
    XLAL_ERROR(XLAL_EDOM, "ModeArrays are not consistent. The (%i,%i) mode in the inertial J-frame requires at least one mode with l=%i in the ModeArray (L-frame) option.\n", ell, emm-1, ell);
  }

  INT4 status = 0; //Variable to check correct functions calls.

  /* Build the frequency array and initialize hctilde to the length of freqs. */
  REAL8Sequence *freqs;
  UINT4 offset = SetupWFArrays(&freqs, hlmpos, freqs_In, pWF, ligotimegps_zero);

  /* Initialize hlmneg according to hlmpos. */
  size_t npts = (*hlmpos)->data->length;
  *hlmneg = XLALCreateCOMPLEX16FrequencySeries("hlmneg: FD waveform", &(*hlmpos)->epoch, (*hlmpos)->f0, pWF->deltaF, &lalStrainUnit, npts);
  XLAL_CHECK (*hlmneg, XLAL_ENOMEM, "Failed to allocated waveform COMPLEX16FrequencySeries of length %zu.", npts);
  memset((*hlmneg)->data->data, 0, npts * sizeof(COMPLEX16));
  XLALUnitMultiply(&((*hlmneg)->sampleUnits), &((*hlmneg)->sampleUnits), &lalSecondUnit);

  /* Variable to store the strain of only one (negative) mode: h_l-mprime */
  COMPLEX16FrequencySeries *htilde22 = NULL;

  /*
     Take input/default value for the threshold of the Multibanding for the hlms modes.
     If = 0 then do not use Multibanding. Default value defined in XLALSimInspiralWaveformParams.c.
     If the input freqs_In is non-uniform the Multibanding has been already switched off.
  */
  REAL8 thresholdMB  = XLALSimInspiralWaveformParamsLookupPhenomXHMThresholdMband(lalParams);

  /* Initialize the power of pi for the HM internal functions. */
  status = IMRPhenomX_Initialize_Powers(&powers_of_lalpiHM, LAL_PI);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initialize useful powers of LAL_PI.");

  UINT4 n_coprec_modes = 0;

  /***** Loop over non-precessing modes ******/
  for (UINT4 emmprime = 1; emmprime <= ell; emmprime++)
  {
    /* Loop over only positive mprime is intentional.
       The single mode function returns the negative mode h_l-mprime, and the positive
       is added automatically in during the twisting up in IMRPhenomXPHMTwistUpOneMode.
       First check if (l,m) mode is 'activated' in the ModeArray.
       If activated then generate the mode, else skip this mode.
    */
    if (XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, emmprime) != 1)
    { /* skip mode */
      continue;
    } /* else: generate mode */

    if(XLALSimInspiralWaveformParamsLookupPhenomXPHMPrecModes(lalParams) == 1 && (INT4)emmprime!=-m)
    {
      continue;
    }

    n_coprec_modes++;

    #if DEBUG == 1
    printf("\n*************************************************\n Non-precessing Mode %i%i\n************************************",ell, emmprime);
    #endif

    /* Variable to store the strain of only one (negative) mode: h_l-mprime */
    COMPLEX16FrequencySeries *htildelm = NULL;

    /* Compute non-precessing mode */
    if (thresholdMB == 0){  // No multibanding
      if(ell == 2 && emmprime == 2)
      {
         status = IMRPhenomXASGenerateFD(&htildelm, freqs_In, pWF, lalParams);
         XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "IMRPhenomXASGenerateFD failed to generate IMRPhenomXHM waveform.");
      }
      else
      {
        status = IMRPhenomXHMGenerateFDOneMode(&htildelm, freqs_In, pWF, ell, emmprime, lalParams);
        XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "IMRPhenomXHMGenerateFDOneMode failed to generate IMRPhenomXHM waveform.");
      }
    }
    else{               // With multibanding
      if(ell==3 && emmprime==2){  // mode with mode-mixing
         status = IMRPhenomXHMMultiBandOneModeMixing(&htildelm, htilde22, pWF, ell, emmprime, lalParams);
         XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "IMRPhenomXHMMultiBandOneModeMixing failed to generate IMRPhenomXHM waveform.");
      }
      else{                  // modes without mode-mixing including 22 mode
         status = IMRPhenomXHMMultiBandOneMode(&htildelm, pWF, ell, emmprime, lalParams);
         XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "IMRPhenomXHMMultiBandOneMode failed to generate IMRPhenomXHM waveform.");
      }

      /* IMRPhenomXHMMultiBandOneMode* functions set pWF->deltaF=0 internally, we put it back here. */
      pWF->deltaF = deltaF;

      /* If the 22 and 32 modes are active, we recycle the 22 mode for the mixing in the 32 and it is passed to IMRPhenomXHMMultiBandOneModeMixing.
         The 22 mode is always computed first than the 32, we store the 22 mode in the variable htilde22. */
      if(ell==2 && emmprime==2 && XLALSimInspiralModeArrayIsModeActive(ModeArray, 3, 2)==1){
        htilde22 = XLALCreateCOMPLEX16FrequencySeries("hptilde: FD waveform", &(ligotimegps_zero), 0.0, pWF->deltaF, &lalStrainUnit, htildelm->data->length);
        for(UINT4 idx = 0; idx < htildelm->data->length; idx++){
          htilde22->data->data[idx] = htildelm->data->data[idx];
        }
      }
    }

    if (!(htildelm)){ XLAL_ERROR(XLAL_EFUNC);}


    /* htildelm is recomputed every time in the loop. Check that it always comes out with the same length */
    XLAL_CHECK (    ((*hlmpos)->data->length==htildelm->data->length)
                 && ((*hlmneg)->data->length==htildelm->data->length),
                 XLAL_EBADLEN,
                 "Inconsistent lengths between frequency series htildelm (%d), hlmpos (%d) and hlmneg (%d).",
                 htildelm->data->length, (*hlmpos)->data->length, (*hlmneg)->data->length
               );

     /* Skip twisting-up if the non-precessing mode is zero. */
     if((pWF->q == 1) && (pWF->chi1L == pWF->chi2L) && (emmprime % 2 != 0))
     {
        XLALDestroyCOMPLEX16FrequencySeries(htildelm);
       continue;
     }
     /*
                              TWISTING UP
          Transform modes from the precessing L-frame to inertial J-frame.
     */

     /* Variable to store the non-precessing waveform in one frequency point. */
     COMPLEX16 hlmcoprec;

     if(XLALSimInspiralWaveformParamsLookupPhenomXPHMPrecModes(lalParams) == 1)
     {
       for (UINT4 idx = 0; idx < freqs->length; idx++)
       {
         hlmcoprec  = htildelm->data->data[idx + offset];  /* Co-precessing waveform */
         if(m < 0) (*hlmpos)->data->data[idx + offset] = hlmcoprec;     // Positive frequencies. Freqs do 0, df, 2df, ...., fmax
         if(m > 0) (*hlmneg)->data->data[idx + offset] = hlmcoprec;     // Negative frequencies. Freqs do 0, -df, -2df, ...., -fmax
       }
     }
     else{
       /* Loop over frequencies. Only where waveform is non zero. */
       for (UINT4 idx = 0; idx < freqs->length; idx++)
       {
          REAL8 Mf = pWF->M_sec*freqs->data[idx];
          hlmcoprec  = htildelm->data->data[idx + offset];  /* Co-precessing waveform */
          COMPLEX16Sequence *hlm;
          hlm = XLALCreateCOMPLEX16Sequence(2);
          IMRPhenomXPHMTwistUpOneMode(Mf, hlmcoprec, pWF, pPrec, ell, emmprime, m, hlm);
          (*hlmpos)->data->data[idx + offset] += hlm->data[0];     // Positive frequencies. Freqs do 0, df, 2df, ...., fmax
          (*hlmneg)->data->data[idx + offset] += hlm->data[1];     // Negative frequencies. Freqs do 0, -df, -2df, ...., -fmax
          XLALDestroyCOMPLEX16Sequence(hlm);
        }
     }
     
     XLALDestroyCOMPLEX16FrequencySeries(htildelm);

  }// End Loop over emmprime

  if(n_coprec_modes == 0)
  {
    XLAL_PRINT_ERROR("For computing the mode (%i,%i) in the inertial J-frame we need at least one l=%i mode activated in the co-precessing L-frame. \nConsider activate some l=%i modes in L-frame with the ModeArray option of the LAL dictionary. \nWe filled the (%i,%i) mode with zeroes." , ell, m, ell, ell, ell, m);
  }

  XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "IMRPhenomXPHM_OneMode failed to generate IMRPhenomXPHM waveform.");

/* Free memory */
XLALDestroyREAL8Sequence(freqs);
XLALDestroyCOMPLEX16FrequencySeries(htilde22);
XLALDestroyValue(ModeArray);


#if DEBUG == 1
printf("\n******Leaving IMRPhenomXPHM_OneMode*****\n");
#endif

return XLAL_SUCCESS;
}


/*
  Core twisting up routine to get one single precessing mode.
  Twist the waveform in the precessing L-frame to the inertial J-frame for one frequency point.
  This function will be inside a loop of frequencies insid a loop over the non-precessing modes.
  It carries out the operation specified in eqs. E3-E4 in the Precessing paper.
*/
static int IMRPhenomXPHMTwistUpOneMode(
  const REAL8 Mf,                          /**< Frequency (Mf geometric units) */
  const COMPLEX16 hlmprime,                /**< Underlying aligned-spin IMRPhenomXHM waveform. The loop is with mprime positive, but the mode has to be the negative one for positive frequencies.*/
  IMRPhenomXWaveformStruct *pWF,           /**< IMRPhenomX Waveform Struct */
  IMRPhenomXPrecessionStruct *pPrec,       /**< IMRPhenomXP Precession Struct */
  UINT4  l,                                /**< l index of the (l,m) (non-)precessing mode */
  UINT4  mprime,                           /**< second index of the (l,mprime) non-precessing mode  */
  INT4   m,                                /**< second index of the (l,m) precessing mode */
  COMPLEX16Sequence *hlminertial           /**< [out] hlm in the inertial J-frame for one frequency point, precessing waveform  */
)
{
  XLAL_CHECK(hlminertial  != NULL, XLAL_EFAULT);

  /* Euler angles */
  double alpha       = 0.0;
  double epsilon     = 0.0;

  double cBetah      = 0.0;
  double sBetah      = 0.0;


  switch(pPrec->IMRPhenomXPrecVersion)
  {
    case 101:        /* Post-Newtonian Euler angles. Single spin approximantion. See sections IV-B and IV-C in Precessing paper. */
    case 102:        /* The different number 10i means different PN order. */
    case 103:
    case 104:
    {
      /* NNLO PN Euler Angles */
      Get_alpha_beta_epsilon(&alpha, &cBetah, &sBetah, &epsilon, mprime, Mf, pPrec, pWF);
      break;
    }
    case 220:
    case 221:
    case 222:
    case 223:
    case 224:
    {
      /* ~~~~~ Euler Angles from Chatziioannou et al, PRD 95, 104004, (2017)  ~~~~~ */
      const double v            = cbrt(LAL_PI * Mf * (2.0/mprime) );
      const vector vangles      = IMRPhenomX_Return_phi_zeta_costhetaL_MSA(v,pWF,pPrec);
      double cos_beta           = 0.0;

      REAL8 alpha_offset_mprime = 0, epsilon_offset_mprime = 0;
      Get_alpha_epsilon_offset(&alpha_offset_mprime, &epsilon_offset_mprime, mprime, pPrec);

      alpha    = vangles.x - alpha_offset_mprime;
      epsilon  = vangles.y - epsilon_offset_mprime;
      cos_beta = vangles.z;
      
      INT4 status = 0;
      status = IMRPhenomXWignerdCoefficients_cosbeta(&cBetah, &sBetah, cos_beta);
      XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "Call to IMRPhenomXWignerdCoefficients_cosbeta failed.");

      break;
    }
    default:
    {
      XLAL_ERROR(XLAL_EINVAL,"Error. IMRPhenomXPrecVersion not recognized. Recommended default is 223.\n");
      break;
    }
  }


  /* Useful powers of the Wigner coefficients */
  REAL8 cBetah2 = cBetah * cBetah;
  REAL8 cBetah3 = cBetah * cBetah2;
  REAL8 cBetah4 = cBetah * cBetah3;
  REAL8 cBetah5 = cBetah * cBetah4;
  REAL8 cBetah6 = cBetah * cBetah5;
  REAL8 cBetah7 = cBetah * cBetah6;
  REAL8 cBetah8 = cBetah * cBetah7;

  REAL8 sBetah2 = sBetah * sBetah;
  REAL8 sBetah3 = sBetah * sBetah2;
  REAL8 sBetah4 = sBetah * sBetah3;
  REAL8 sBetah5 = sBetah * sBetah4;
  REAL8 sBetah6 = sBetah * sBetah5;
  REAL8 sBetah7 = sBetah * sBetah6;
  REAL8 sBetah8 = sBetah * sBetah7;

  /*
      The following expressions for the Wigner-d coefficients correspond to those in appendix A of the Precessing paper.
      They are the same expressions used in IMRPhenomXPHMTwistUp.
  */

  COMPLEX16 hlm = 0, hlmneg = 0;
  INT4 minus1l = 1;

  /* Sum over l = 2 modes */
  if (l == 2 && mprime == 2){
    //                  d^2_{-2,2}    d^2_{-1,2}             d^2_{0,2}                   d^2_{1,2}       d^2_{2,2}
    COMPLEX16 d22[5]   = {sBetah4, 2.0*cBetah*sBetah3, pPrec->sqrt6*sBetah2*cBetah2, 2.0*cBetah3*sBetah, cBetah4};
    //                  d^2_{-2,-2}  d^2_{-1,-2}  d^2_{0,-2}  d^2_{1,-2}    d^2_{2,-2}
    COMPLEX16 d2m2[5]  = {d22[4],    -d22[3],      d22[2],     -d22[1],     d22[0]}; /* Exploit symmetry d^2_{-m,-2} = (-1)^m d^2_{-m,2}. See eq. A2 of Precessing paper */

    /* See eqs. E3-E4 in Precessing paper. */
    COMPLEX16 cexp_im_alpha = cexp(-1.*I*m*alpha);
    hlm += cexp_im_alpha * d2m2[m+2];
    hlmneg += cexp_im_alpha * d22[m+2];
  }

  /* Case (l, mprime) = (2, 1) */
  if (l == 2 && mprime == 1){
    //                    d^2_{-2,1}          d^2_{-1,1}                       d^2_{0,1}                                       d^2_{1,1}                        d^2_{2,1}
    COMPLEX16 d21[5]   = {2.0*cBetah*sBetah3, 3.0*sBetah2*cBetah2 - sBetah4, pPrec->sqrt6*(cBetah3*sBetah - cBetah*sBetah3), cBetah2*(cBetah2 - 3.0*sBetah2), -2.0*cBetah3*sBetah};
    //                  d^2_{-2,-1}  d^2_{-1,-1}  d^2_{0,-1}  d^2_{1,-1}  d^2_{2,-1}
    COMPLEX16 d2m1[5]  = {-d21[4],   d21[3],     -d21[2],    d21[1],     -d21[0]}; /* Exploit symmetry d^2_{-m,-1} = -(-1)^m d^2_{m,1}.  See eq. A2 of Precessing paper.  */

    /* See eqs. E3-E4 in Precessing paper. */
    COMPLEX16 cexp_im_alpha = cexp(-1.*I*m*alpha);
    hlm += cexp_im_alpha * d2m1[m+2];
    hlmneg += cexp_im_alpha * d21[m+2];
  }

  /* Case (l, mprime) = (3, 3) */
  if (l == 3 && mprime == 3){
    minus1l = -1;
    //                  d^3_{-3,3}    d^3_{-2,3}                 d^3_{-1,3}                    d^3_{0,3}                         d^3_{1,3}                      d^3_{2,3}                   d^3_{3,3}
    COMPLEX16 d33[7]   = {sBetah6, pPrec->sqrt6*cBetah*sBetah5, pPrec->sqrt15*cBetah2*sBetah4, 2.0*pPrec->sqrt5*cBetah3*sBetah3, pPrec->sqrt15*cBetah4*sBetah2, pPrec->sqrt6*cBetah5*sBetah, cBetah6};
    //                  d^3_{-3,-3}  d^3_{-2,-3}  d^3_{-1,-3}  d^3_{0,-3}  d^3_{1,-3}  d^3_{2,-3}  d^3_{3,-3}
    COMPLEX16 d3m3[7]  = {d33[6],    -d33[5],     d33[4],      -d33[3],    d33[2],     -d33[1],    d33[0]}; /* Exploit symmetry d^3_{-m,-3} = -(-1)^m d^3_{m,3}. See eq. A2 of Precessing paper. */

    /* See eqs. E3-E4 in Precessing paper. */
    COMPLEX16 cexp_im_alpha = cexp(-1.*I*m*alpha);
    hlm += cexp_im_alpha * d3m3[m+3];
    hlmneg += cexp_im_alpha * d33[m+3];
  }

  /* Case (l, mprime) = (3, 2) */
  if (l == 3 && mprime == 2){
    minus1l = -1;
    //                    d^3_{-3,2}                     d^3_{-2,2}                     d^3_{-1,21}                                            d^3_{0,2}                                         d^3_{1,2}                                              d^3_{2,2}                        d^3_{3,2}
    COMPLEX16 d32[7]   = {pPrec->sqrt6*cBetah*sBetah5, sBetah4*(5.0*cBetah2 - sBetah2), pPrec->sqrt10*sBetah3*(2.0*cBetah3 - cBetah*sBetah2), pPrec->sqrt30*cBetah2*(cBetah2 - sBetah2)*sBetah2, pPrec->sqrt10*cBetah3*(cBetah2*sBetah - 2.0*sBetah3), cBetah4*(cBetah2 - 5.0*sBetah2), -1.0*pPrec->sqrt6*cBetah5*sBetah};
    //                  d^3_{-3,-2}  d^3_{-2,-2}  d^3_{-1,-2}  d^3_{0,-2}  d^3_{1,-2}  d^3_{2,-2}  d^3_{3,-2}
    COMPLEX16 d3m2[7]  = {-d32[6],   d32[5],     -d32[4],      d32[3],     -d32[2],    d32[1],    -d32[0]}; /* Exploit symmetry d^3_{-m,-2} = (-1)^m d^3_{m,2}. See eq. A2 of Precessing paper.  */

    /* See eqs. E3-E4 in Precessing paper. */
    COMPLEX16 cexp_im_alpha = cexp(-1.*I*m*alpha);
    hlm += cexp_im_alpha * d3m2[m+3];
    hlmneg += cexp_im_alpha * d32[m+3];
  }

  /* Case (l, mprime) = (4, 4) */
  if (l == 4 && mprime == 4){
    //                    d^4_{-4,4}         d^4_{-3,4}               d^4_{-2,4}                         d^4_{-1,4}                          d^4_{0,4}                       d^4_{1,4}                           d^4_{2,4}                          d^4_{3,4}                  d^4_{44}
    COMPLEX16 d44[9]   = {sBetah8, 2.0*pPrec->sqrt2*cBetah*sBetah7, 2.0*pPrec->sqrt7*cBetah2*sBetah6, 2.0*pPrec->sqrt14*cBetah3*sBetah5, pPrec->sqrt70*cBetah4*sBetah4, 2.0*pPrec->sqrt14*cBetah5*sBetah3, 2.0*pPrec->sqrt7*cBetah6*sBetah2, 2.0*pPrec->sqrt2*cBetah7*sBetah, cBetah8};
    //                  d^4_{4,-4}  d^4_{-3,-4}  d^4_{-2,-4}  d^4_{-1,-4}  d^4_{0,-4}  d^4_{1,-4}  d^4_{2,-4}  d^4_{3,-4}  d^4_{4-4}
    COMPLEX16 d4m4[9]  = {d44[8],   -d44[7],      d44[6],     -d44[5],     d44[4],    -d44[3],     d44[2],    -d44[1],     d44[0]}; /* Exploit symmetry d^4_{-m,-4} = (-1)^m d^4_{m,4}. See eq. A2 of Precessing paper.  */

    /* See eqs. E3-E4 in Precessing paper. */
    COMPLEX16 cexp_im_alpha = cexp(-1.*I*m*alpha);
    hlm += cexp_im_alpha * d4m4[m+4];
    hlmneg += cexp_im_alpha * d44[m+4];
  }

  /* See eqs. E3-E4 in Precessing paper. */
  COMPLEX16 exp_imprime_epsilon = cexp(mprime*I*epsilon);
  COMPLEX16 eps_phase_hP_lmprime = 1./exp_imprime_epsilon * hlmprime;
  COMPLEX16 eps_phase_hP_lmprime_neg = exp_imprime_epsilon * minus1l * conj(hlmprime);

 /* Return h_lminertail */
 (hlminertial)->data[0] = eps_phase_hP_lmprime * hlm;
 (hlminertial)->data[1] = eps_phase_hP_lmprime_neg * hlmneg;


  return XLAL_SUCCESS;
}


/** Function to obtain a SphHarmFrequencySeries with the individual modes h_lm.
    By default it returns all the modes available in the model, positive and negatives.
    With the mode array option in the LAL dictionary, the user can specify a custom mode array.
    The modes are computed in the inertial J-frame, so the mode array option does not refers to
    the modes in the co-precessing frame conversely to the functions for the polarizations XLALSimIMRPhenomXPHM.
    This function is to be used by ChooseFDModes.
*/
int XLALSimIMRPhenomXPHMModes(
      SphHarmFrequencySeries **hlms,              /**< [out] list with single modes h_lm in the J-frame */
  	  REAL8 m1_SI,                                /**< mass of companion 1 (kg) */
      REAL8 m2_SI,                                /**< mass of companion 2 (kg) */
      REAL8 S1x,                                  /**< x-component of the dimensionless spin of object 1  w.r.t. Lhat = (0,0,1) */
      REAL8 S1y,                                  /**< y-component of the dimensionless spin of object 1  w.r.t. Lhat = (0,0,1) */
      REAL8 S1z,                                  /**< z-component of the dimensionless spin of object 1  w.r.t. Lhat = (0,0,1) */
      REAL8 S2x,                                  /**< x-component of the dimensionless spin of object 2  w.r.t. Lhat = (0,0,1) */
      REAL8 S2y,                                  /**< y-component of the dimensionless spin of object 2  w.r.t. Lhat = (0,0,1) */
      REAL8 S2z,                                  /**< z-component of the dimensionless spin of object 2  w.r.t. Lhat = (0,0,1) */
      REAL8 deltaF,                               /**< frequency spacing (Hz) */
  		REAL8 f_min,                                /**< starting GW frequency (Hz) */
  		REAL8 f_max,                                /**< ending GW frequency (Hz) */
      REAL8 f_ref,                                /**< reference GW frequency (Hz) */
      REAL8 phiRef,                               /**< phase shift at reference frequency */
      REAL8 distance,                             /**< distance of source (m) */
      REAL8 inclination,                          /**< inclination of source (rad) */
  		LALDict *lalParams                          /**< LAL dictionary with extra options */
)
{
    LALValue *ModeArrayJframe = NULL;  // Modes in the precessing J-frame. Specified through the new LAL dictionary option "ModeArrayJframe"
    
    /* Use an auxiliar laldict to not overwrite the input argument */
    LALDict *lalParams_aux;
    /* setup mode array */
    if (lalParams == NULL)
    {
        lalParams_aux = XLALCreateDict();
    }
    else{
        lalParams_aux = XLALDictDuplicate(lalParams);
    }
    
    /* Check that the co-precessing modes chosen are available for the model */
    XLAL_CHECK(check_input_mode_array(lalParams_aux) == XLAL_SUCCESS, XLAL_EFAULT, "Not available co-precessing mode chosen.\n");
    
    
    /* Read mode array from LAL dictionary */
    ModeArrayJframe = XLALSimInspiralWaveformParamsLookupModeArrayJframe(lalParams_aux);
        
    /* If input LAL dictionary does not have mode array, use all the modes available for XPHM (l<=4)  */
    if(ModeArrayJframe == NULL)
    {
      ModeArrayJframe = XLALSimInspiralModeArrayActivateAllModesAtL(XLALSimInspiralCreateModeArray(), 2);
      ModeArrayJframe = XLALSimInspiralModeArrayActivateAllModesAtL(ModeArrayJframe, 3);
      ModeArrayJframe = XLALSimInspiralModeArrayActivateAllModesAtL(ModeArrayJframe, 4);
    }
    else{
      /* Check that the modes chosen are available for the model */
      XLAL_CHECK(check_input_mode_array_Jframe(ModeArrayJframe) == XLAL_SUCCESS, XLAL_EFAULT, "Not available mode chosen. l must be lower than %i\n", L_MAX);
    }

  
    INT4 length = 0;
    /***** Loop over modes ******/
    for (UINT4 ell = 2; ell <= LAL_SIM_L_MAX_MODE_ARRAY; ell++)
    {
      for (INT4 emm = -(INT4)ell; emm <= (INT4)ell; emm++)
      {
        if(XLALSimInspiralModeArrayIsModeActive(ModeArrayJframe, ell, emm) !=1)
        {
          /* Skip mode if user did not specified it. */
          continue;
        }
        //Variable to store the strain of only one (positive/negative) mode: h_lm
        COMPLEX16FrequencySeries *hlmpos = NULL;
        COMPLEX16FrequencySeries *hlmneg = NULL;
        COMPLEX16FrequencySeries *hlmall = NULL;

        /* Compute precessing single mode */
        XLALSimIMRPhenomXPHMOneMode(&hlmpos, &hlmneg, ell, emm, m1_SI, m2_SI, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, deltaF, f_min, f_max, f_ref, lalParams_aux);

        if (!(hlmpos) || !hlmneg){ XLAL_ERROR(XLAL_EFUNC);}

        length = hlmpos->data->length-1;

        hlmall = XLALCreateCOMPLEX16FrequencySeries("hlmall: precessing FD mode",  &(hlmpos->epoch), hlmpos->f0, hlmpos->deltaF, &(hlmpos->sampleUnits), 2*length+1);

        for(INT4 i=0; i<=length; i++)
        {
          hlmall->data->data[i+length] = hlmpos->data->data[i];
          hlmall->data->data[i] = hlmneg->data->data[length-i];
        }


        // Add single mode to list
        *hlms = XLALSphHarmFrequencySeriesAddMode(*hlms, hlmall, ell, emm);

        // Free memory
        XLALDestroyCOMPLEX16FrequencySeries(hlmpos);
        XLALDestroyCOMPLEX16FrequencySeries(hlmneg);
        XLALDestroyCOMPLEX16FrequencySeries(hlmall);
      }
    } /* End loop over modes */


    /* Add frequency array to SphHarmFrequencySeries */
    REAL8Sequence *freqs = XLALCreateREAL8Sequence(2*length+1);
    for (INT4 i = -length; i<=length; i++)
    {
      freqs->data[i+length] = i*deltaF;
    }
    XLALSphHarmFrequencySeriesSetFData(*hlms, freqs);

    /* Free memory */
    XLALDestroyDict(lalParams_aux);
    XLALDestroyValue(ModeArrayJframe);
    

    return XLAL_SUCCESS;

}


/*********************************************/
/*                                           */
/*           AUXILIARY FUNCTIONS             */
/*                                           */
/*********************************************/


/** Wrapper function for setup ModeArray of modes in the precessing frame to be twisted up. */
LALDict *IMRPhenomXPHM_setup_mode_array(LALDict *lalParams)
{
  /* setup ModeArray */
  INT4 lalParams_In = 0;
  if (lalParams == NULL)
  {
    lalParams_In = 1;
    lalParams = XLALCreateDict();
  }
  LALValue *ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(lalParams);

  /* If the mode array is empty, populate using a default choice of modes */
  if (ModeArray == NULL)
  {
    /* Default behaviour */
    XLAL_PRINT_INFO("Using default non-precessing modes for IMRPhenomXPHM: 2|2|, 2|1|, 3|3|, 3|2|, 4|4|.\n");
    ModeArray = XLALSimInspiralCreateModeArray();

    /* IMRPhenomXHM has the following calibrated modes. 22 mode taken from IMRPhenomXAS */
    XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 2);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 1);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 3, 3);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 3, 2);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 4, 4);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 2, -2);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 2, -1);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 3, -3);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 3, -2);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 4, -4);
    XLALSimInspiralWaveformParamsInsertModeArray(lalParams, ModeArray);
  }
  else {XLAL_PRINT_INFO("Using custom non-precessing modes for PhenomXPHM.\n"); }

  XLALDestroyValue(ModeArray);
  if(lalParams_In == 1)
  {
    XLALDestroyDict(lalParams);
  }

  return lalParams;
}


/* Function to check if the input mode array in the J-frame contains unsupported modes */
INT4 check_input_mode_array_Jframe(LALValue *ModeArrayJframe){
  /* Check if the input array has a too high l. */
  for(INT4 ell=2; ell<=LAL_SIM_L_MAX_MODE_ARRAY; ell++)
  {
    for(INT4 emm=0; emm<=ell; emm++)
    {
      if(XLALSimInspiralModeArrayIsModeActive(ModeArrayJframe, ell, emm) == 1 && ell>L_MAX){
        XLALDestroyValue(ModeArrayJframe);
        return XLAL_FAILURE;
      }
    }
  }
  return XLAL_SUCCESS;
}

/*
  Return the offset at reference frequency for alpha and epsilon Euler angles for a particular non-precessing mode.
  E.g.: alpha_offset_mprime = alpha(2*pi*MfRef/mprime) - alpha0. Used for Pv2 and Pv3 angles.
*/
static double Get_alpha_epsilon_offset(
  REAL8 *alpha_offset_mprime,          /**< [out] Offset alpha angle at reference frequency */
  REAL8 *epsilon_offset_mprime,        /**< [out] Offset epsilon angle at reference frequency */
  INT4 mprime,                         /**< Second index of the non-precesssing mode (l, mprime) */
  IMRPhenomXPrecessionStruct *pPrec    /**< IMRPhenomXP Precessing structure*/
)
{
  /* The angles are evaluated at frequency 2*pi*MfRef/mprime so the offset depend on mprime.  */
  switch(mprime)
  {
     case 1:
     {
       *alpha_offset_mprime   = pPrec->alpha_offset_1;
       *epsilon_offset_mprime = pPrec->epsilon_offset_1;
       break;
     }
     case 2:
     {
       *alpha_offset_mprime   = pPrec->alpha_offset;   // These variable was already used in the XP code,
       *epsilon_offset_mprime = pPrec->epsilon_offset; // that is why we did not add the _2 here.
       break;
     }
     case 3:
     {
       *alpha_offset_mprime   = pPrec->alpha_offset_3;
       *epsilon_offset_mprime = pPrec->epsilon_offset_3;
       break;
     }
     case 4:
     {
       *alpha_offset_mprime   = pPrec->alpha_offset_4;
       *epsilon_offset_mprime = pPrec->epsilon_offset_4;
       break;
     }
     default:
     {
       XLAL_ERROR(XLAL_EINVAL,"Error: mprime not supported, check available non-precessing modes.\n");
       break;
     }
  }

  return XLAL_SUCCESS;
}


/* Return the 3 Post-Newtonian Euler angles evaluated at frequency (2*pi*Mf/mprime). */
/* See equations C19/C20 and G7/G8 in the Precessing paper for the expressions of alpha/epsilon. beta is computed accoring to eq. 4.4. */
static int Get_alpha_beta_epsilon(
  REAL8 *alpha,                       /**< [out] Azimuthal angle of L w.r.t J */
  REAL8 *cBetah,                      /**< [out] Cosine of polar angle between L and J */
  REAL8 *sBetah,                      /**< [out] Sine of polar angle between L and J */
  REAL8 *epsilon,                     /**< [out] Minus the third Euler angle (-gamma) describing L w.r.t J, fixed by minimal rotation condition */
  INT4 mprime,                        /**< Second index of the non-precesssing mode (l, mprime) */
  REAL8 Mf,                           /**< Frequency geometric units */
  IMRPhenomXPrecessionStruct *pPrec,  /**< IMRPhenomXP Precessing structure*/
  IMRPhenomXWaveformStruct *pWF       /**< IMRPhenomX Waveform structure*/
)
{
  double omega       = LAL_PI * Mf *2./mprime;
  double logomega    = log(omega);
  double omega_cbrt  = cbrt(omega);
  double omega_cbrt2 = omega_cbrt * omega_cbrt;

  REAL8 alpha_offset_mprime = 0., epsilon_offset_mprime = 0.;

  Get_alpha_epsilon_offset(&alpha_offset_mprime, &epsilon_offset_mprime, mprime, pPrec);

  *alpha = (
      pPrec->alpha1  / omega
    + pPrec->alpha2  / omega_cbrt2
    + pPrec->alpha3  / omega_cbrt
    + pPrec->alpha4L * logomega
    + pPrec->alpha5  * omega_cbrt
    - alpha_offset_mprime
  );

  *epsilon =  (
        pPrec->epsilon1  / omega
      + pPrec->epsilon2  / omega_cbrt2
      + pPrec->epsilon3  / omega_cbrt
      + pPrec->epsilon4L * logomega
      + pPrec->epsilon5  * omega_cbrt
      - epsilon_offset_mprime
  );

  INT4 status = 0;
  status = IMRPhenomXWignerdCoefficients(cBetah, sBetah, omega_cbrt, pWF, pPrec);
  XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "Call to IMRPhenomXWignerdCoefficients failed.");

  return status;
}


/*
   Non-uniform/coarse frequency grid for the Multibanding of the angles.

   - In this first release we use the same coarse grid that is used for computing the non-precessing modes.
   - This grid is discussed in section II of arXiv:2001.10897. See also section D of Precessing paper.
   - This grid is computed with the function XLALSimIMRPhenomXMultibandingVersion defined in LALSimIMRPhenomXHM_multiband.c.
   - The version of the coarse grid will be changed with the option 'MBandPrecVersion' defined in LALSimInspiralWaveformParams.c.
   - Currently there is only one version available and the option value for that is 0, which is the default value.
*/
INT4 XLALSimIMRPhenomXPHMMultibandingGrid(
  REAL8Sequence **coarseFreqs,      /**<[out] Non-uniform coarse frequency grid (1D array) */
  UINT4 ell,                        /**< First index non-precessing mode */
  UINT4 emmprime,                   /**< Second index non-precessing mode */
  IMRPhenomXWaveformStruct *pWF,    /**< IMRPhenomX Waveform Struct*/
  LALDict *lalParams)               /**< LAL dictionary */
{
  /* This function is basically a copy of the first part of IMRPhenomXHMMultiBandOneMode and IMRPhenomXHMMultiBandOneModeMixing. */

 /* Create non-uniform grid for each mode. */
 REAL8 thresholdMB  = XLALSimInspiralWaveformParamsLookupPhenomXPHMThresholdMband(lalParams);

 /* Compute the coarse frequency array. It is stored in a list of grids. */
 size_t iStart = (size_t) (pWF->fMin / pWF->deltaF);

 /* Final grid spacing, adimensional (NR) units */
 REAL8 evaldMf = XLALSimIMRPhenomXUtilsHztoMf(pWF->deltaF, pWF->Mtot);

 /* Variable for the Multibanding criteria. See eq. 2.8-2.9 of arXiv:2001.10897. */
 REAL8 dfpower = 11./6.;
 REAL8 dfcoefficient = 8. * sqrt(3./5.) * LAL_PI * powers_of_lalpi.m_one_sixth * sqrt(2.)*cbrt(2) /(cbrt(emmprime)*emmprime) * sqrt(thresholdMB * pWF->eta);

 /* Variables for the coarse frequency grid */
 REAL8 Mfmin = XLALSimIMRPhenomXUtilsHztoMf(iStart*pWF->deltaF, pWF->Mtot);
 REAL8 Mfmax = XLALSimIMRPhenomXUtilsHztoMf(pWF->f_max_prime, pWF->Mtot);
 REAL8 MfMECO, MfLorentzianEnd;
 REAL8 dfmerger = 0., dfringdown = 0.;
 UINT4 lengthallGrids = 20;

 IMRPhenomXMultiBandingGridStruct *allGrids = (IMRPhenomXMultiBandingGridStruct*)XLALMalloc(lengthallGrids * sizeof(IMRPhenomXMultiBandingGridStruct));
 IMRPhenomXHMWaveformStruct *pWFHM = (IMRPhenomXHMWaveformStruct *) XLALMalloc(sizeof(IMRPhenomXHMWaveformStruct));

 //populate coefficients of 22 mode, to rotate to spherical
 IMRPhenomXAmpCoefficients   *pAmp22   = (IMRPhenomXAmpCoefficients *) XLALMalloc(sizeof(IMRPhenomXAmpCoefficients));
 IMRPhenomXPhaseCoefficients *pPhase22 = (IMRPhenomXPhaseCoefficients *) XLALMalloc(sizeof(IMRPhenomXPhaseCoefficients));
 IMRPhenomXGetPhaseCoefficients(pWF, pPhase22);

 /* Allocate and initialize the PhenomXHM lm amplitude coefficients struct */
 IMRPhenomXHMAmpCoefficients *pAmp = (IMRPhenomXHMAmpCoefficients*)XLALMalloc(sizeof(IMRPhenomXHMAmpCoefficients));
 IMRPhenomXHMPhaseCoefficients *pPhase = (IMRPhenomXHMPhaseCoefficients*)XLALMalloc(sizeof(IMRPhenomXHMPhaseCoefficients));

 if(ell == 2 && emmprime ==2){
   MfMECO = pWF->fMECO;
   #if DEBUG == 1
   printf("\nfRING = %e\n",pWF->fRING);
   printf("fDAMP = %e\n",pWF->fDAMP);
   printf("alphaL22 = %.16e", pPhase22->cLovfda/pWF->eta);
   #endif
   MfLorentzianEnd = pWF->fRING + 2*pWF->fDAMP;
   IMRPhenomXGetAmplitudeCoefficients(pWF, pAmp22);
   dfmerger = deltaF_mergerBin(pWF->fDAMP, pPhase22->cLovfda/pWF->eta, thresholdMB);
   dfringdown = deltaF_ringdownBin(pWF->fDAMP, pPhase22->cLovfda/pWF->eta, pAmp22->gamma2/(pAmp22->gamma3*pWF->fDAMP),thresholdMB);
 }
 else{
   // allocate qnm struct
   QNMFits *qnms = (QNMFits *) XLALMalloc(sizeof(QNMFits));
   IMRPhenomXHM_Initialize_QNMs(qnms);
   // Populate pWFHM
   IMRPhenomXHM_SetHMWaveformVariables(ell, emmprime, pWFHM, pWF, qnms, lalParams);
   LALFree(qnms);

   /* Allocate and initialize the PhenomXHM lm phase and amp coefficients struct */
   IMRPhenomXHM_FillAmpFitsArray(pAmp);
   IMRPhenomXHM_FillPhaseFitsArray(pPhase);

   if(pWFHM->MixingOn == 1){
     /* Get coefficients for Amplitude and phase */
     GetSpheroidalCoefficients(pPhase, pPhase22, pWFHM, pWF);
     IMRPhenomXGetAmplitudeCoefficients(pWF, pAmp22);
   }

   /* Get coefficients for Amplitude and phase */
   IMRPhenomXHM_GetAmplitudeCoefficients(pAmp, pPhase, pAmp22, pPhase22, pWFHM, pWF);
   IMRPhenomXHM_GetPhaseCoefficients(pAmp, pPhase, pAmp22, pPhase22, pWFHM, pWF,lalParams);

   MfMECO = pWFHM->fMECOlm;
   MfLorentzianEnd = pWFHM->fRING + 2*pWFHM->fDAMP;
   #if DEBUG == 1
   printf("\nfRING = %e\n",pWFHM->fRING);
   printf("fDAMP = %e\n",pWFHM->fDAMP);
   printf("alphaL = %.16e", pPhase->alphaL);
   #endif
   dfmerger = deltaF_mergerBin(pWFHM->fDAMP, pPhase->alphaL, thresholdMB);
   dfringdown = deltaF_ringdownBin(pWFHM->fDAMP, pPhase->alphaL, pAmp->lambda/(pAmp->sigma*pWFHM->fDAMP), thresholdMB);
 }
 LALFree(pWFHM);
 LALFree(pAmp);
 LALFree(pPhase);
 LALFree(pAmp22);
 LALFree(pPhase22);

 UINT4 nGridsUsed = XLALSimIMRPhenomXMultibandingGrid(Mfmin, MfMECO, MfLorentzianEnd, Mfmax, evaldMf, dfpower, dfcoefficient, allGrids, dfmerger, dfringdown);
 if (allGrids == NULL)
 {
   #if DEBUG == 1
   printf("Malloc of allGrids failed!\n");
   printf("nGridsUsed = %i\n", nGridsUsed);
   #endif
   return -1;
 }

 /* Number of fine frequencies per coarse interval in every coarse grid */
 /* Actual number of subgrids to be used. We allocated more than needed. */
 UINT4 actualnumberofGrids = 0;
 /* Length of coarse frequency array */
 UINT4 lenCoarseArray = 0;

 /* Transform the coarse frequency array to 1D array. */
 // Take only the subgrids needed
 for(UINT4 kk = 0; kk < nGridsUsed; kk++){
   lenCoarseArray = lenCoarseArray + allGrids[kk].Length;
   actualnumberofGrids++;

   #if DEBUG == 1
   printf("\nkk = %i\n",kk);
   printf("xStart: %.16e\n", allGrids[kk].xStart);
   printf("xEnd: %.16e\n", allGrids[kk].xEndRequested);
   printf("Length: %i\n", allGrids[kk].Length);
   printf("deltax: %.16e\n", allGrids[kk].deltax);
   printf("evaldMf: %.16e\n", evaldMf);
   printf("xMax: %.16e\n", allGrids[kk].xMax);
   printf("Last grid.xMax = %.16f\n", allGrids[actualnumberofGrids-1].xMax);
   #endif

   if(allGrids[kk].xMax + evaldMf >= Mfmax){
     break;
   }
 }

 // Add extra points to the coarse grid if the last freq is lower than Mfmax
 while(allGrids[actualnumberofGrids-1].xMax < Mfmax){
   allGrids[actualnumberofGrids-1].xMax =   allGrids[actualnumberofGrids-1].xMax +  allGrids[actualnumberofGrids-1].deltax;
   allGrids[actualnumberofGrids-1].Length =  allGrids[actualnumberofGrids-1].Length + 1;
   lenCoarseArray++;
 }

 #if DEBUG == 1
 printf("\nactualnumberofGrids = %i\n", actualnumberofGrids);
 printf("lenCoarseArray = %i\n", lenCoarseArray);
 printf("Last grid.xMax = %.16f", allGrids[actualnumberofGrids-1].xMax);
 #endif

 // Transform coarse frequency array to 1D vector
 /* Allocate memory for frequency array and terminate if this fails */
 (*coarseFreqs) = XLALCreateREAL8Sequence(lenCoarseArray);
 if (!(*coarseFreqs)) { XLAL_ERROR(XLAL_EFUNC, "Frequency array allocation failed.");}

 UINT4 lencoarseFreqs = 0;

 for(UINT4 kk = 0; kk < actualnumberofGrids; kk++){
   for(INT4 ll = 0; ll < allGrids[kk].Length; ll++){
     (*coarseFreqs)->data[lencoarseFreqs] = (allGrids[kk].xStart + allGrids[kk].deltax*ll);
     lencoarseFreqs++;
   }
 }

 /* End of coarse frequency array. */

 #if DEBUG == 1
 printf("\n******** Coarse frequencies array done ********* \n");
 printf("\nlencoarseFreqs, coarse[0], coarse[-1], Mfmax = %i %.16e %.16e %.16e\n",lencoarseFreqs, XLALSimIMRPhenomXUtilsMftoHz((*coarseFreqs)->data[0],pWF->Mtot), XLALSimIMRPhenomXUtilsMftoHz((*coarseFreqs)->data[lencoarseFreqs-1], pWF->Mtot), XLALSimIMRPhenomXUtilsMftoHz(Mfmax,pWF->Mtot));
 #endif


 LALFree(allGrids);
 return actualnumberofGrids;
}

/* @} */
/* @} */
