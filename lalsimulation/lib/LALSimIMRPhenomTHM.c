/*
 *  Copyright (C) 2020 Hector Estelles
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

/* Standard LAL */
#include <lal/Sequence.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/XLALError.h>

/* LAL datatypes and constants */
#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>

/* Time series, frequency series and spherical harmonics */
#include <lal/TimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/SphericalHarmonics.h>
#include <lal/LALSimSphHarmMode.h>
#include <lal/FrequencySeries.h>

/* LALSimulation */
#include <lal/LALSimIMR.h>
#include <lal/LALSimInspiral.h>

/* Standard C */
#include <math.h>
#include <complex.h>
#include <stdbool.h>

/* Link IMRPhenomTHM routines */
#include "LALSimIMRPhenomTHM_internals.h"
#include "LALSimIMRPhenomTHM_internals.c"
#include "LALSimIMRPhenomTHM.h"

/* *********** ALIGNED SPIN IMR MULTIMODE PHENOMENOLOGICAL TIME DOMAIN WAVEFORM MODEL: IMRPhenomTHM *********** */

/* EXTERNAL ROUTINES */

/**
 * @addtogroup LALSimIMRPhenomT_c
 * @brief Routines to produce IMRPhenomT-family of phenomenological inspiral-merger-ringdown waveforms.
 * These are time-domain models for compact binaries at comparable and extreme mass ratios, tuned to numerical-relativity simulations.
 *
 * - IMRPhenomT: model for the 22 mode of non-precessing binaries. Model concept together with precession was published in https://arxiv.org/abs/2004.08302. 
 * Although an updated version together with the higher mode extension is under the DCC: https://dcc.ligo.org/LIGO-P2000524.  
 * - IMRPhenomTHM: model with subdominant harmonics for non-precessing binaries. DCC link: https://dcc.ligo.org/LIGO-P2000524.
 * - IMRPhenomTP: model for precessing binaries. Twisted-up version of IMRPhenomT. Original concept: https://arxiv.org/abs/2004.08302.
 * - IMRPhenomTPHM: model for precessing binaries including higher harmonics. Twisted-up version of IMRPhenomTHM.
 *
 *
 * @review IMRPhenomT(HM) has been review by Maria Haney, Jacob Lange, Jonathan Thompson and Eleanor Hamilton. Review wiki: https://git.ligo.org/waveforms/reviews/phenomt/-/wikis/home.
 * IMRPhenomTP(HM) are being reviewed by the same team under the same wiki page.
 *
 *
 */

 /**
  * @addtogroup LALSimIMRPhenomT_c
  * @{
  *
  * @name Routines for IMRPhenomT
  * @{
  *
  * @author Héctor Estellés
  *
  * @brief C code for the IMRPhenomT phenomenological waveform model
  *
  *
  *
  */
  
/**
 * Routine to compute time domain polarisations for IMRPhenomT model. It returns two real time series for the plus and the cross polarisations,
 * constructed with the modes (l,m)=(2,2), (2,-2).
 */

int XLALSimIMRPhenomT(
  REAL8TimeSeries **hp, /**< [out] TD waveform for plus polarisation */
  REAL8TimeSeries **hc, /**< [out] TD waveform for cross polarisation */
  REAL8 m1_SI,      /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,      /**< Mass of companion 2 (kg) */
  REAL8 chi1L,      /**< Dimensionless aligned spin of companion 1 */
  REAL8 chi2L,      /**< Dimensionless aligned spin of companion 2 */
  REAL8 distance,   /**< Luminosity distance (m) */
  REAL8 inclination,  /**< inclination of source (rad) */
  REAL8 deltaT,     /**< sampling interval (s) */
  REAL8 fmin,     /**< starting GW frequency (Hz) */
  REAL8 fRef,     /**< reference GW frequency (Hz) */
  REAL8 phiRef,     /**< reference orbital phase (rad) */
  LALDict *lalParams  /**< LAL dictionary containing accessory parameters */
  )
{

  INT4 status;

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
  LALValue *ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(lalParams_aux);
  
  /* Check if an mode array is NULL and then populate it with the (2,2) and (2,-2) modes. */
  if(ModeArray==NULL)
  {
    ModeArray = XLALSimInspiralCreateModeArray();
    /* Activate (2,2) and (2,-2) modes*/
    XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 2);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 2, -2);
    /* Insert ModeArray into lalParams */
    XLALSimInspiralWaveformParamsInsertModeArray(lalParams_aux, ModeArray);
    XLALDestroyValue(ModeArray);
  }

  UINT4 only22 = 1; /* Internal flag for mode array check */

  /* Compute modes dominant modes */
  SphHarmTimeSeries *hlms = NULL;
  status = LALSimIMRPhenomTHM_Modes(&hlms, m1_SI, m2_SI, chi1L, chi2L, distance, deltaT, fmin, fRef, phiRef, lalParams_aux, only22);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: Internal function LALSimIMRPhenomTHM_Modes has failed producing the modes.");

  /* Obtain length and epoch from modes (they are the same for all the modes) */
  INT4 length = hlms->mode->data->length;
  LIGOTimeGPS epoch = hlms->mode->epoch;

  /* Auxiliary time series for setting the polarisations to zero */
  REAL8TimeSeries *hplus = XLALCreateREAL8TimeSeries ("hplus", &epoch, 0.0, deltaT, &lalStrainUnit,length);
  REAL8TimeSeries *hcross = XLALCreateREAL8TimeSeries ("hcross", &epoch, 0.0, deltaT, &lalStrainUnit,length);
  memset( hplus->data->data, 0, hplus->data->length * sizeof(REAL8) );
  memset( hcross->data->data, 0, hcross->data->length * sizeof(REAL8) );

  /* Add the modes to the polarisations using lalsim routines.
  Negative modes are explicitely passed, instead of using the symmetry flag */

  SphHarmTimeSeries *hlms_temp = hlms;
  while ( hlms_temp )
  {
      XLALSimAddMode(hplus, hcross, hlms_temp->mode, inclination, LAL_PI/2. - phiRef, hlms_temp->l, hlms_temp->m, 0);
      hlms_temp = hlms_temp->next;
  }

  /* Point the output pointers to the relevant time series */
  (*hp) = hplus;
  (*hc) = hcross;

  /* Destroy intermediate time series */
  XLALDestroySphHarmTimeSeries(hlms);
  XLALDestroySphHarmTimeSeries(hlms_temp);

  /* Destroy lalParams_aux. */
  XLALDestroyDict(lalParams_aux);
  

  return status;
}

/** @} */
/**
* @name Routines for IMRPhenomTHM
* @{
*
* @author Héctor Estellés
*
* @brief C code for the IMRPhenomTHM phenomenological waveform model
*
* IMRPhenomTHM is a phenomenological model in the time domain for aligned-spin binary black hole coalescences, calibrated to Numerical Relativity simulations for
* comparable mass ratio and Teukolsky waveforms for extreme mass ratio. The model produces IMR waveforms for the Spin-Weighted Spherical Harmonic modes (l,m)=(2,2), (2,1), (3,3), (4,4) and (5,5),
* obtaining the corresponding negative m modes by symmetry.
*
* For selecting a particular list of modes to be returned or to be employed in the polarisations construction, the user can follow the usual procedure:
* - Create a mode array object with lalsimulation.SimInspiralCreateModeArray
* - Activate the desired modes with lalsim.SimInspiralModeArrayActivateMode
* - Insert the mode array into a LAL dictionary with lalsim.SimInspiralWaveformParamsInsertModeArray
* - Pass the LAL ditionary to ChooseTDWaveform or ChooseTDModes.
* For a user specified mode array, only the implemented modes in the model will be computed.
*
* User option for selecting (2,2) phase and frequency reconstruction through LAL parameter PhenomTHMInspiralVersion. Default (0) will reconstruct with
* only 1 inspiral region modelled by TaylorT3 with the merger time parameter tt0 set to 0 and additional higher order terms.
* Non-default will provide an early inspiral region for imensionless PN time parameter theta<0.33, constructed with pure TaylorT3 (without higher orders) with tt0 calibrated
* across parameter space. Both options maintained for code historical reasons, but not clear benefit from non-default option was found.
*
* Model has been calibrated to 531 BBH non-precessing NR simulations from the
* last release of the SXS Catalog (https://iopscience.iop.org/article/10.1088/1361-6382/ab34e2), additional BAM NR simulations at q=4, q=8 and q=18, and numerical Teukolsky waveforms placed at q=200 and q=1000. Calibration procedure has followed the 
* hierarchical data-driven fitting approach (Xisco Jimenez-Forteza et al https://arxiv.org/abs/1611.00332) using the symmetric mass ratio eta, dimensionless effective spin Shat=(m1^2*chi1+m2^2*chi2)/(m1^2+m2^2)
* and spin difference dchi=chi1-chi2. Supplementary material for the fits of the various collocation points
* and phenomenological coefficients is available at https://git.ligo.org/waveforms/reviews/phenomt/-/tree/master/SupplementaryMaterial/Fits3DPhenomTHM.
*
*/

/**
 * Routine to compute time domain polarisations for IMRPhenomTHM model. It returns two real time series for the plus and the cross polarisations,
 * constructed with the modes specified in lalParams that correspond to included modes in the model: (l,m)=(2,2), (2,1), (3,3), (4,4) and (5,5) and
 * the corresponding negative m modes.
 */

int XLALSimIMRPhenomTHM(
  REAL8TimeSeries **hp, /**< [out] TD waveform for plus polarisation */
  REAL8TimeSeries **hc, /**< [out] TD waveform for cross polarisation */
  REAL8 m1_SI,			/**< Mass of companion 1 (kg) */
  REAL8 m2_SI,			/**< Mass of companion 2 (kg) */
  REAL8 chi1L,			/**< Dimensionless aligned spin of companion 1 */
  REAL8 chi2L,			/**< Dimensionless aligned spin of companion 2 */
  REAL8 distance,		/**< Luminosity distance (m) */
  REAL8 inclination,	/**< inclination of source (rad) */
  REAL8 deltaT,			/**< sampling interval (s) */
  REAL8 fmin,			/**< starting GW frequency (Hz) */
  REAL8 fRef,			/**< reference GW frequency (Hz) */
  REAL8 phiRef,			/**< reference orbital phase (rad) */
  LALDict *lalParams 	/**< LAL dictionary containing accessory parameters */
  )
{

  int status;

  /* Sanity checks pointers. */
  XLAL_CHECK(NULL != hp, XLAL_EFAULT);
  XLAL_CHECK(NULL != hc, XLAL_EFAULT);
  XLAL_CHECK(*hp == NULL, XLAL_EFAULT);
  XLAL_CHECK(*hc == NULL, XLAL_EFAULT);

  UINT4 only22 = 0; /* Internal flag for mode array check */

  /* Compute modes (if a subset of modes is desired, it should be passed through lalParams) */
  SphHarmTimeSeries *hlms = NULL;
  status = LALSimIMRPhenomTHM_Modes(&hlms, m1_SI, m2_SI, chi1L, chi2L, distance, deltaT, fmin, fRef, phiRef, lalParams, only22);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: Internal function LALSimIMRPhenomTHM_Modes has failed producing the modes.");

  /* Obtain length and epoch from modes (they are the same for all the modes) */
  INT4 length = hlms->mode->data->length;
  LIGOTimeGPS epoch = hlms->mode->epoch;

  /* Auxiliary time series for setting the polarisations to zero */
  REAL8TimeSeries *hplus = XLALCreateREAL8TimeSeries ("hplus", &epoch, 0.0, deltaT, &lalStrainUnit,length);
  REAL8TimeSeries *hcross = XLALCreateREAL8TimeSeries ("hcross", &epoch, 0.0, deltaT, &lalStrainUnit,length);
  memset( hplus->data->data, 0, hplus->data->length * sizeof(REAL8) );
  memset( hcross->data->data, 0, hcross->data->length * sizeof(REAL8) );

  /* Add the modes to the polarisations using lalsim routines.
  Negative modes are explicitely passed, instead of using the symmetry flag */

  SphHarmTimeSeries *hlms_temp = hlms;
  while ( hlms_temp )
  {
      XLALSimAddMode(hplus, hcross, hlms_temp->mode, inclination, LAL_PI/2. - phiRef, hlms_temp->l, hlms_temp->m, 0);
      hlms_temp = hlms_temp->next;
  }

  /* Point the output pointers to the relevant time series */
  (*hp) = hplus;
  (*hc) = hcross;

  /* Destroy intermediate time series */
  XLALDestroySphHarmTimeSeries(hlms);
  XLALDestroySphHarmTimeSeries(hlms_temp);

  return status;
}

/**
 * Routine to compute time domain Spin-Weighted Spherical Harmonic modes for IMRPhenomTHM model. It returns a SphHarmTimeSeries structure, with complex time series for each mode generated.
 */

SphHarmTimeSeries *XLALSimIMRPhenomTHM_Modes(
  REAL8 m1_SI,								/**< Mass of companion 1 (kg) */
  REAL8 m2_SI,								/**< Mass of companion 2 (kg) */
  REAL8 chi1L,								/**< Dimensionless aligned spin of companion 1 */
  REAL8 chi2L,								/**< Dimensionless aligned spin of companion 2 */
  REAL8 distance,							/**< Luminosity distance (m) */
  REAL8 deltaT,								/**< sampling interval (s) */
  REAL8 fmin,								/**< starting GW frequency (Hz) */
  REAL8 fRef,								/**< reference GW frequency (Hz) */
  REAL8 phiRef,								/**< reference orbital phase (rad) */
  LALDict *lalParams 						/**< LAL dictionary containing accessory parameters */
  )
{

  INT4 status;

  UINT4 only22 = 0; /* Internal flag for mode array check */

  /* Compute modes (if a subset of modes is desired, it should be passed through lalParams */
  SphHarmTimeSeries *hlms = NULL;
  status = LALSimIMRPhenomTHM_Modes(&hlms, m1_SI, m2_SI, chi1L, chi2L, distance, deltaT, fmin, fRef, phiRef, lalParams, only22);
  XLAL_CHECK_NULL(status != XLAL_FAILURE, XLAL_EFUNC, "Error: Internal function LALSimIMRPhenomTHM_Modes has failed producing the modes.");

  return hlms;
}

/** @} */
/** @} */

/*
 * Internal routine to compute time domain Spin-Weighted Spherical Harmonic modes for IMRPhenomTHM model. It returns a SphHarmTimeSeries structure, with complex time series for each mode generated.
 */

int LALSimIMRPhenomTHM_Modes(
  SphHarmTimeSeries **hlms,   /**< [out] Time domain modes */
  REAL8 m1_SI,                /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,                /**< Mass of companion 2 (kg) */
  REAL8 chi1L,                /**< Dimensionless aligned spin of companion 1 */
  REAL8 chi2L,                /**< Dimensionless aligned spin of companion 2 */
  REAL8 distance,             /**< Luminosity distance (m) */
  REAL8 deltaT,               /**< sampling interval (s) */
  REAL8 fmin,               /**< starting GW frequency (Hz) */
  REAL8 fRef,               /**< reference GW frequency (Hz) */
  REAL8 phiRef,               /**< reference orbital phase (rad) */
  LALDict *lalParams,            /**< LAL dictionary containing accessory parameters */
  UINT4 only22              /** Internal flag for mode array check */
  )
{

  /* Pointers sanity check */
  XLAL_CHECK(NULL != hlms, XLAL_EFAULT);
  XLAL_CHECK(*hlms == NULL, XLAL_EFAULT);

  /* Sanity checks */
  if(fRef  <  0.0) { XLAL_ERROR(XLAL_EDOM, "fRef_In must be positive or set to 0 to ignore.\n");  }
  if(deltaT   <= 0.0) { XLAL_ERROR(XLAL_EDOM, "deltaF must be positive.\n");                         }
  if(m1_SI    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "m1 must be positive.\n");                             }
  if(m2_SI    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "m2 must be positive.\n");                             }
  if(fmin    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "f_min must be positive.\n");                          }
  if(distance <= 0.0) { XLAL_ERROR(XLAL_EDOM, "Distance must be positive and greater than 0.\n");    }

  /* Check that the modes chosen are available for the model. Included ones are (2,2), (2,1), (3,3), (4,4) and (5,5) and corresponding negative m ones.
  If only22=1, then checks that only (2,2) (2,-2) are included for dominant mode approximant */

  if(only22==0)
  {
    XLAL_CHECK(check_input_mode_array_THM(lalParams) == XLAL_SUCCESS, XLAL_EFAULT, "Not available mode chosen.\n");
  }
  else if(only22==1)
  {
    XLAL_CHECK(check_input_mode_array_22_THM(lalParams) == XLAL_SUCCESS, XLAL_EFAULT, "Not available mode chosen.\n");
  }

  /*
    Perform a basic sanity check on the region of the parameter space in which model is evaluated. Behaviour is as follows:
      - For mass ratios <= 20.0 and spins <= 0.99: no warning messages.
      - For 200 > mass ratio > 20 and spins <= 0.9: print a warning message that model can be pathological.
      - For mass ratios > 200: throw a hard error that model is not valid.
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
  if(mass_ratio > 20.0 && chi1L < 0.9 && m2_SI/LAL_MSUN_SI >= 0.5  ) { XLAL_PRINT_INFO("Warning: Extrapolating outside of Numerical Relativity calibration domain."); }
  if(mass_ratio > 20.0 && (chi1L >= 0.9 || m2_SI/LAL_MSUN_SI < 0.5) ) { XLAL_PRINT_INFO("Warning: Model can be pathological at these parameters."); }
  if(mass_ratio > 200. && fabs(mass_ratio - 200) > 1e-12) { XLAL_ERROR(XLAL_EDOM, "ERROR: Model not valid at mass ratios beyond 200."); } // The 1e-12 is to avoid rounding errors
  if(fabs(chi1L) > 0.99 || fabs(chi2L) > 0.99) { XLAL_PRINT_INFO("Warning: Extrapolating to extremal spins, model is not trusted."); }

  INT4 status;

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
  lalParams_aux = IMRPhenomTHM_setup_mode_array(lalParams_aux);
  LALValue *ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(lalParams_aux);

  /* Initialise waveform struct and phase/frequency struct for the 22, needed for any mode */

  IMRPhenomTWaveformStruct *pWF;
  pWF    = XLALMalloc(sizeof(IMRPhenomTWaveformStruct));
  status = IMRPhenomTSetWaveformVariables(pWF, m1_SI, m2_SI, chi1L, chi2L, distance, deltaT, fmin, fRef, phiRef, lalParams_aux);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: Internal function IMRPhenomTSetWaveformVariables has failed.");


  IMRPhenomTPhase22Struct *pPhase;
  pPhase = XLALMalloc(sizeof(IMRPhenomTPhase22Struct));
  status   = IMRPhenomTSetPhase22Coefficients(pPhase, pWF);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: Internal function IMRPhenomTSetPhase22Coefficients has failed.");

  /* Length of the required sequence and length of the different inspiral regions of the 22 phase. */
  size_t length = pPhase->wflength;
  size_t length_insp_early = pPhase->wflength_insp_early;
  size_t length_insp_late = pPhase->wflength_insp_late;
  size_t length_insp = length_insp_early + length_insp_late;

  /*Initialize time */
  REAL8 t, thetabar, theta, w22, ph22;
  REAL8 factheta = pow(5.0,1./8);

  /* Initialize REAL8 sequences for storing the phase, frequency and PN expansion parameter x of the 22 */
  REAL8Sequence *phi22 = NULL;
  REAL8Sequence *xorb = NULL;
  xorb = XLALCreateREAL8Sequence(length);
  phi22 = XLALCreateREAL8Sequence(length);
  

  /* Compute 22 phase and x=(0.5*omega_22)^(2/3), needed for all modes. 
  Depending on the inspiral region, theta is defined with a fitted t0 parameter or with t0=0. */

  if(pWF->inspVersion!=0) // If reconstruction is non-default, this will compute frequency, phase and PN expansion parameter x from TaylorT3 with fitted t0 for the early inspiral region-
  {
    for(UINT4 jdx = 0; jdx < length_insp_early; jdx++)
    {
      t = pPhase->tmin + jdx*pWF->dtM;

      thetabar = pow(pWF->eta*(pPhase->tt0-t),-1./8);
    
      theta = factheta*thetabar;

      w22 = IMRPhenomTomega22(t, theta, pWF, pPhase);
      xorb->data[jdx] = pow(0.5*w22,2./3);

      ph22 = IMRPhenomTPhase22(t, thetabar, pWF, pPhase);
    
      phi22->data[jdx] = ph22;
    }

    for(UINT4 jdx = length_insp_early; jdx < length_insp; jdx++) // For times later than the early-late inspiral boundary, it computes phase, frequency and x with the extended TaylorT3 (with tt0=0)
    {
      t = pPhase->tmin + jdx*pWF->dtM;

      thetabar = pow(-pWF->eta*t,-1./8);
    
      theta = factheta*thetabar;

      w22 = IMRPhenomTomega22(t, theta, pWF, pPhase);
      xorb->data[jdx] = pow(0.5*w22,2./3);

      ph22 = IMRPhenomTPhase22(t, thetabar, pWF, pPhase);
    
      phi22->data[jdx] = ph22;
    }
  }

  else // If default reconstruction, only one inspiral region with extended TaylorT3 (with tt0=0) is employed up to the inspiral-merger boundary time
  {
      for(UINT4 jdx = 0; jdx < length_insp; jdx++)
      {
        t = pPhase->tmin + jdx*pWF->dtM;

        thetabar = pow(-pWF->eta*t,-1./8);
    
        theta = factheta*thetabar;

        w22 = IMRPhenomTomega22(t, theta, pWF, pPhase);
        xorb->data[jdx] = pow(0.5*w22,2./3);

        ph22 = IMRPhenomTPhase22(t, thetabar, pWF, pPhase);
    
        phi22->data[jdx] = ph22;
      }

  }

  /* During the 22 phase merger region, phase and x are also computed because the higher modes end the inspiral at a fixed time,
  which can occur later than the end of the 22 inspiral. */
  for(UINT4 jdx = length_insp; jdx < length; jdx++)
  {
    t = pPhase->tmin + jdx*pWF->dtM;

    w22 = IMRPhenomTomega22(t, 0.0, pWF, pPhase);
    xorb->data[jdx] = pow(0.5*w22,2./3);

    ph22 = IMRPhenomTPhase22(t, 0.0, pWF, pPhase);
    
    phi22->data[jdx] = ph22;
  }

  /***** Loop over modes ******/

  INT4 posMode, negMode;

  for (UINT4 ell = 2; ell <= 5; ell++)
  {
    for (UINT4 emm = 1; emm <= ell; emm++)
    { /* Loop over only positive m is intentional.*/

      posMode = XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, emm);
      negMode = XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, -emm);
      if ( posMode != 1 && negMode != 1)
      { /* skip mode */
        continue;
      } /* else: generate mode */

      COMPLEX16TimeSeries *tmp_mode = NULL;

      status = LALSimIMRPhenomTHM_OneMode(&tmp_mode, pWF, pPhase, phi22, xorb, ell, emm);

      if(posMode==1) /* Modes are generated only for positive m. If positive m is requested, simply add to the SphHarmTimeSeries structure */
      {
        *hlms = XLALSphHarmTimeSeriesAddMode(*hlms, tmp_mode, ell, emm);
      }

      if(negMode==1) /* If negative m is requested, transform from positive m mode using symmetry property for aligned spin systems */
      {
        for(UINT4 jdx = 0; jdx < length; jdx++)
        {
          ((*tmp_mode).data)->data[jdx] = pow(-1,ell)*conj(((*tmp_mode).data)->data[jdx]);
        }
        *hlms = XLALSphHarmTimeSeriesAddMode(*hlms, tmp_mode, ell, -emm);
      }

      XLALDestroyCOMPLEX16TimeSeries(tmp_mode); /* Destroy temporary structure once the desired mode is added */
    }

  }

  /*Free structures and destroy sequences and dict */
  XLALDestroyValue(ModeArray);
  LALFree(pPhase);
  LALFree(pWF);

  XLALDestroyREAL8Sequence(phi22);

  XLALDestroyREAL8Sequence(xorb);

  XLALDestroyDict(lalParams_aux);

  return status;
}

/* Internal function for generating one mode. See section II of PhenomTHM paper for details on mode construction: https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf */

int LALSimIMRPhenomTHM_OneMode(
	COMPLEX16TimeSeries **hlm,			/**< [out] Time domain waveform of the requested (l,m) mode */
	IMRPhenomTWaveformStruct *pWF,		/**< Waveform structure */
	IMRPhenomTPhase22Struct *pPhase,	/**< 22 phase and frequency structure */
  REAL8Sequence *phase22,       /**< Values of the 22 phase for the waveform time array */
	REAL8Sequence *xorb,				/**< Values of the 22 frequency for the waveform time array */
	UINT4 ell,							/**< l value of the requested mode */
	UINT4 emm							/**< m value of the requested mode */
	)
{

  /* Sanity checks were performed in XLALSimIMRPhenomTHM_Modes, from which this function is called. */
  UINT4 status;

  /*Length of the required sequence. Read it from 22 phase structure. */
  size_t length = pPhase->wflength;

  /*Declare time, amplitude and phase */
  REAL8 t;

  COMPLEX16 amplm;
  UNUSED REAL8 philm;

  /* Set LIGOTimeGPS */
  LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0,0}
  XLALGPSAdd(&ligotimegps_zero, pPhase->tminSec); /* By model construction, this implies that the 22 amplitude peaks exactly at t=0 */

  /* Initialize time series for mode */
  *hlm = XLALCreateCOMPLEX16TimeSeries("hlm: TD waveform",&ligotimegps_zero,0.0,pWF->deltaT,&lalStrainUnit,length);

  /* Initialize variables for 22 phase and frequency */
  REAL8 x = 0.0;
  UNUSED REAL8 phi22 = 0.0;
  UNUSED COMPLEX16 expPhi;
  COMPLEX16 wf;

  /* phiref0 is the value of the 22 phase at the reference frequency, which will be substracted in order to impose the chosen orbital reference phase */
  REAL8 thetabarRef = pow(-pPhase->tRef*pWF->eta,-1./8);
  REAL8 phiref0 = IMRPhenomTPhase22(pPhase->tRef, thetabarRef, pWF, pPhase);

  IMRPhenomTHMAmpStruct *pAmplm;
  pAmplm = XLALMalloc(sizeof(IMRPhenomTHMAmpStruct));
  status   = IMRPhenomTSetHMAmplitudeCoefficients(ell,emm, pAmplm, pPhase, pWF);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomTSetHMAmplitudeCoefficients failed for %d,%d.\n",ell,emm);

	/* For l=2, m=2 mode, phase already computed. Amplitude is computed in absolute value since we don't want the complex phase contribution in the 22 (frequency is already calibrated in the inspiral) */
  if(ell==2 && emm==2)
	{
     /* Loop over time */
    /* Mode construction, see equation 2 of THM paper https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf */
    for(UINT4 jdx = 0; jdx < length; jdx++)
    {
        t = pPhase->tmin + jdx*pWF->dtM;
        x = (xorb)->data[jdx];
        amplm = pWF->ampfac*cabs(IMRPhenomTHMAmp(t, x, pAmplm));

        phi22 = (phase22)->data[jdx] - phiref0;
        wf = amplm*cexp(-I*phi22);

        ((*hlm)->data->data)[jdx] = wf;
    }
	}

  else if(emm%2 != 0 && pWF->delta<1E-10 && fabs(pWF->chi1L-pWF->chi2L)<1E-10) // This is for not computing odd modes in the equal BH limit. Instead, they are set to 0 directly.
  {

    memset( (*hlm)->data->data, 0.0, (*hlm)->data->length * sizeof(COMPLEX16) );
        
  }
	
  else // Higher modes
	{
		/* Initialize phase struct for desired mode */

		  IMRPhenomTHMPhaseStruct *pPhaselm;
    	pPhaselm = XLALMalloc(sizeof(IMRPhenomTHMPhaseStruct));
    	status   = IMRPhenomTSetHMPhaseCoefficients(ell,emm, pPhaselm, pPhase, pAmplm, pWF);
    	XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomTSetHMPhaseCoefficients failed for %d,%d.\n",ell,emm);

      /* Phase offset, as explained in eq 13 of PhenomTHM paper (https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf). */
      
      REAL8 phoff = 0.0; // Constant phase offset employed for rotating PN amplitudes so much of the content is real.
      if(ell==2 && emm==1)
      {
        phoff = 0.5*LAL_PI;
      }
      else if(ell==3 && emm==3)
      {
        phoff = -0.5*LAL_PI;
      }
      else if(ell==5 && emm==5)
      {
        phoff = 0.5*LAL_PI;
      }
      else if(ell==4 && emm==4)
      {
        phoff = LAL_PI;
      }

      /* Loop over time */
      for(UINT4 jdx = 0; jdx < length; jdx++)
      {

        t = pPhase->tmin + jdx*pWF->dtM;
        x = (xorb)->data[jdx]; // 22 frequency, needed for evaluating inspiral amplitude
        amplm = pWF->ampfac*IMRPhenomTHMAmp(t, x, pAmplm);
        
        phi22 = (phase22)->data[jdx]; // 22 phase, needed for rescaling lm inspiral phase
        
        philm = IMRPhenomTHMPhase(t, phi22, pPhaselm,pAmplm) - (emm/2.0)*phiref0 - phoff; // Phase construction
        wf = amplm*cexp(-I*philm); /* Mode construction. Note than amplitude is now complex during the inspiral, which will contribute to the phase of the mode */

        ((*hlm)->data->data)[jdx] = wf;

      }
      LALFree(pPhaselm);

	}

  LALFree(pAmplm);

  return status;
}

/* Wrapper function for adding higher modes to the ModeArray, almost identical to IMRPhenomXHM_setup_mode_array */

LALDict *IMRPhenomTHM_setup_mode_array(LALDict *lalParams)
{
  LALValue *ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(lalParams);

  /* If the mode array is empty, populate using a default choice of modes */
  if (ModeArray == NULL)
  {
    /* Default behaviour */
    XLAL_PRINT_INFO("Using default modes for IMRPhenomTHM.\n");
    ModeArray = XLALSimInspiralCreateModeArray();


    /* IMRPhenomTHM has the following calibrated modes.*/
    XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 2);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 1);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 3, 3);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 4, 4);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 5, 5);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 2, -2);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 2, -1);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 3, -3);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 4, -4);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 5, -5);

    XLALSimInspiralWaveformParamsInsertModeArray(lalParams, ModeArray);
  }
  else {XLAL_PRINT_INFO("Using custom modes for PhenomTHM.\n"); }
  XLALDestroyValue(ModeArray);

  return lalParams;
}

/* Function for checking if input mode array has supported modes by the model. Ported from IMRPhenomXHM, with modifications in the supported modes. */
INT4 check_input_mode_array_THM(LALDict *lalParams)
{
	UINT4 flagTrue = 0;
	
  if(lalParams == NULL) return XLAL_SUCCESS;
  
  LALValue *ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(lalParams);
  
  if(ModeArray!=NULL)
  {
    INT4 larray[5] = {2, 2, 3, 4, 5};
    INT4 marray[5] = {2, 1, 3, 4, 5};
    
    for(INT4 ell=2; ell<=LAL_SIM_L_MAX_MODE_ARRAY; ell++)
		{
			for(INT4 emm=0; emm<=ell; emm++)
			{
				if(XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, emm)==1 || XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, -1*emm)==1)
				{
					for(UINT4 idx=0; idx<5; idx++)
					{
			      if(ell==larray[idx] && abs(emm)==marray[idx])
						{
							flagTrue = 1;
						}
					}
					// If flagTrue !=1 means that the input mode array has a mode that is not supported by the model.
					if(flagTrue!=1){
						XLALPrintError ("Mode (%d,%d) is not available by the model.\n", ell, emm);
						XLALDestroyValue(ModeArray);
						return XLAL_FAILURE;
					}
					flagTrue = 0;					
				}				
			}//End loop over emm
		}//End loop over ell
  }//End of if block
	
  XLALDestroyValue(ModeArray);
	
  return XLAL_SUCCESS;
}

/* Function for checking if input mode array has supported modes by the model. Ported from IMRPhenomXHM, with modifications in the supported modes. In this cases, for IMRPhenomT, only 22 and 2-2 modes are supported. */
INT4 check_input_mode_array_22_THM(LALDict *lalParams)
{
  UINT4 flagTrue = 0;
  
  if(lalParams == NULL) return XLAL_SUCCESS;
  
  LALValue *ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(lalParams);
  
  if(ModeArray!=NULL)
  { 
    for(INT4 ell=2; ell<=LAL_SIM_L_MAX_MODE_ARRAY; ell++)
    {
      for(INT4 emm=0; emm<=ell; emm++)
      {
        if(XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, emm)==1 || XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, -1*emm)==1)
        {
          
          if(ell==2 && abs(emm)==2)
          {
            flagTrue = 1;
          }
          
          // If flagTrue !=1 means that the input mode array has a mode that is not supported by the model.
          if(flagTrue!=1){
            XLALPrintError ("Mode (%d,%d) is not available by the model.\n", ell, emm);
            XLALDestroyValue(ModeArray);
            return XLAL_FAILURE;
          }
          flagTrue = 0;         
        }       
      }//End loop over emm
    }//End loop over ell
  }//End of if block
  
  XLALDestroyValue(ModeArray);
  
  return XLAL_SUCCESS;
}
