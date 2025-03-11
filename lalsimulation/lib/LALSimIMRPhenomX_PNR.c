/*
 * Copyright (C) 2022 Cardiff University
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

#ifdef __cplusplus
extern "C"
{
#endif

#include <lal/SphericalHarmonics.h>
#include "LALSimIMR.h"

#include "LALSimIMRPhenomX.h"
#include "LALSimIMRPhenomXPHM.h"
#include "LALSimIMRPhenomX_internals.h"
#include "LALSimIMRPhenomXUtilities.h"
#include "LALSimIMRPhenomX_precession.h"

#include "LALSimIMRPhenomX_PNR.h"
#include "LALSimIMRPhenomX_PNR_internals.h"
#include "LALSimIMRPhenomX_PNR_internals.c"

#ifndef _OPENMP
#define omp ignore
#endif

#ifndef PHENOMXHMDEBUG
#define DEBUG 0
#else
#define DEBUG 1
#endif

  /**
   * @addtogroup LALSimIMRPhenomX_c
   * @{
   *
   * @name Routines for PNR angles
   * @{
   *
   * @author Eleanor Hamilton, Sebastian Khan, Jonathan E. Thompson
   *
   * @brief C code for routines to implement PNR angles in IMRPhenomXP/IMRPhenomXPHM.
   *
   * This is a C-code impementation of the PNR angle model
   * outlined in arXiv:2107.08876.
   *
   * Available flags:
   *  PhenomXPNRUseTunedAngles:
   *    - 0 : Disable PNR angles. Use the default NNLO/MSA precession angles in IMRPhenomXP/HM without tuning.
   *    - 1 : Enable PNR angles.
   */

  /**
   * This is an external wrapper to generate the (2,2) PNR angles,
   * following the prescription outlined in arXiv:2107.08876,
   * given the standard inputs given to generate FD waveforms.
   */
  int XLALSimIMRPhenomX_PNR_GeneratePNRAngles(
      REAL8Sequence **alphaPNR, /**< [out] Alpha precession angle (rad) */
      REAL8Sequence **betaPNR,  /**< [out] Beta precession angle (rad) */
      REAL8Sequence **gammaPNR, /**< [out] Gamma precession angle (rad) */
      REAL8Sequence **freqs,    /**< [out] Frequency array (Hz) */
      REAL8 *alphaPNR_ref,      /**< [out] reference value of alpha (rad) */
      REAL8 *gammaPNR_ref,      /**< [out] reference value of gamma (rad) */
      REAL8 m1_SI,              /**< mass of companion 1 (kg) */
      REAL8 m2_SI,              /**< mass of companion 2 (kg) */
      REAL8 chi1x,              /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
      REAL8 chi1y,              /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
      REAL8 chi1z,              /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
      REAL8 chi2x,              /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
      REAL8 chi2y,              /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
      REAL8 chi2z,              /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
      REAL8 inclination,        /**< Angle between orbital angular momentum and line-of-sight vector at reference frequency (rad) */
      REAL8 deltaF,             /**< Frequency spacing (Hz) */
      REAL8 f_min,              /**< Starting GW frequency (Hz) */
      REAL8 f_max,              /**< Ending GW frequency (Hz) */
      REAL8 fRef_In,            /**< Reference frequency (Hz) */
      LALDict *lalParams        /**< LAL Dictionary struct */
  )
  {
    /* Simple check on masses and spins */
    UINT4 status = 0;
    status = XLALIMRPhenomXPCheckMassesAndSpins(&m1_SI, &m2_SI, &chi1x, &chi1y, &chi1z, &chi2x, &chi2y, &chi2z);
    XLAL_CHECK(
        XLAL_SUCCESS == status,
        XLAL_EFUNC,
        "Error: XLALIMRPhenomXPCheckMassesAndSpins failed in XLALSimIMRPhenomX_PNR_GeneratePNRAngles.\n");

    /* Ensure we have a dictionary */
    LALDict *lalParams_aux;
    if (lalParams == NULL)
    {
      lalParams_aux = XLALCreateDict();
    }
    else
    {
      lalParams_aux = XLALDictDuplicate(lalParams);
    }

    /* Make sure we're calling the tuned angles */
    int UsePNR = XLALSimInspiralWaveformParamsLookupPhenomXPNRUseTunedAngles(lalParams_aux);
    if (!UsePNR)
    {
      UsePNR = 1;
      XLALSimInspiralWaveformParamsInsertPhenomXPNRUseTunedAngles(lalParams_aux, UsePNR);
    }

    /* Map fRef to the start frequency if need be, then make sure it's within the frequency range */
    REAL8 fRef = (fRef_In == 0.0) ? f_min : fRef_In;
    XLAL_CHECK(
        ((f_min <= fRef) && (fRef <= f_max)),
        XLAL_EFUNC,
        "Error: fRef needs to be within the specified minimum and maximum frequency values!\n");

    /* Check that the passed deltaF is non-zero */
    XLAL_CHECK(
        deltaF > 0,
        XLAL_EFUNC,
        "Error: deltaF needs to be a positive number!\n");

    /* Generate a uniformly sampled frequency grid of spacing deltaF. */
    /* Frequencies will be set using only the lower and upper bounds that we passed */
    size_t iStart = (size_t)(f_min / deltaF);
    size_t iStop = (size_t)(f_max / deltaF) + 1;

    XLAL_CHECK(
        (iStart <= iStop),
        XLAL_EDOM,
        "Error: the starting frequency index is greater than the stopping index! Please ensure that f_min <= f_max.\n");

    /* Allocate memory for frequency array and terminate if this fails */
    *freqs = XLALCreateREAL8Sequence(iStop - iStart);
    if (!(*freqs))
    {
      XLAL_ERROR(XLAL_EFUNC, "Frequency array allocation failed.");
    }
    *alphaPNR = XLALCreateREAL8Sequence(iStop - iStart);
    if (!(*alphaPNR))
    {
      XLAL_ERROR(XLAL_EFUNC, "Alpha array allocation failed.");
    }
    *betaPNR = XLALCreateREAL8Sequence(iStop - iStart);
    if (!(*betaPNR))
    {
      XLAL_ERROR(XLAL_EFUNC, "Beta array allocation failed.");
    }
    *gammaPNR = XLALCreateREAL8Sequence(iStop - iStart);
    if (!(*gammaPNR))
    {
      XLAL_ERROR(XLAL_EFUNC, "Gamma array allocation failed.");
    }

    /* Populate frequency array */
    for (UINT4 i = iStart; i < iStop; i++)
    {
      (*freqs)->data[i - iStart] = i * deltaF;
    }

    /* Specify arbitrary parameters for the angle generation */
    REAL8 distance = 1.0;
    REAL8 phiRef = 0.0;

    /* Use the true minimum and maximum frequency values */
    REAL8 f_min_eval = (*freqs)->data[0];
    REAL8 f_max_eval = (*freqs)->data[(*freqs)->length - 1];

    /* Initialize PhenomX Waveform struct and check that it initialized correctly */
    IMRPhenomXWaveformStruct *pWF;
    pWF = XLALMalloc(sizeof(IMRPhenomXWaveformStruct));
    status = IMRPhenomXSetWaveformVariables(pWF, m1_SI, m2_SI, chi1z, chi2z, deltaF, fRef, phiRef, f_min_eval, f_max_eval, distance, inclination, lalParams_aux, DEBUG);
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetWaveformVariables failed.\n");

    /* Initialize PhenomX Precession struct and check that it generated successfully */
    IMRPhenomXPrecessionStruct *pPrec;
    pPrec = XLALMalloc(sizeof(IMRPhenomXPrecessionStruct));
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
        DEBUG);
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetPrecessionVariables failed.\n");

    /* See if PNR angles were turned off in pPrec check */
    UsePNR = pPrec->IMRPhenomXPNRUseTunedAngles;

    if(UsePNR){
      /* Generate the tuned precession angles */
      status = IMRPhenomX_PNR_GeneratePNRAngles(*alphaPNR, *betaPNR, *gammaPNR, *freqs, pWF, pPrec, lalParams_aux);
      XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomX_PNR_GeneratePNRAngles failed.\n");

      /* Pass through the reference values of alpha and gamma
      * NOTE: gamma = - epsilon */
      *alphaPNR_ref = pPrec->alpha_offset;
      *gammaPNR_ref = -(pPrec->epsilon_offset);
    }
    else{
      /* Generate PN angles */
      if(XLALSimInspiralWaveformParamsLookupPhenomXPrecVersion(lalParams_aux) > 200 ){
        /* MSA angles generated using built-in function */

        REAL8Sequence *cosBeta = XLALCreateREAL8Sequence(iStop - iStart);
        if (!cosBeta)
        {
          XLAL_ERROR(XLAL_EFUNC, "cosBeta array allocation failed.");
        }

        status = XLALSimIMRPhenomXPMSAAngles(
          alphaPNR, gammaPNR, &cosBeta, *freqs,
          m1_SI,m2_SI,
          chi1x,chi1y,chi1z,
          chi2x,chi2y,chi2z,
          inclination,
          fRef_In, 2, lalParams_aux
        );

        for (UINT4 i = 0; i < (*freqs)->length; i++)
        {
          (*betaPNR)->data[i] = acos(cosBeta->data[i]);
        }

        XLALDestroyREAL8Sequence(cosBeta);

        /* Reference values already applied in function */
        *alphaPNR_ref = 0.0;
        *gammaPNR_ref = 0.0;
      }
      else{
        /* NNLO angles */

        REAL8Sequence *cosBeta = XLALCreateREAL8Sequence(iStop - iStart);
        if (!cosBeta)
        {
          XLAL_ERROR(XLAL_EFUNC, "cosBeta array allocation failed.");
        }

        status = XLALSimIMRPhenomXPPNAngles(
          alphaPNR, gammaPNR, &cosBeta, *freqs,
          m1_SI,m2_SI,
          chi1x,chi1y,chi1z,
          chi2x,chi2y,chi2z,
          inclination,
          fRef_In, 2, lalParams_aux
        );

        for (UINT4 i = 0; i < (*freqs)->length; i++)
        {
          (*betaPNR)->data[i] = acos(cosBeta->data[i]);
        }

        XLALDestroyREAL8Sequence(cosBeta);

        /* Pass through the reference values of alpha and gamma
        * NOTE: gamma = - epsilon */
        *alphaPNR_ref = pPrec->alpha_offset;
        *gammaPNR_ref = -(pPrec->epsilon_offset);
      }
    }

    /* Clean up memory allocation */
    LALFree(pPrec);
    LALFree(pWF);
    XLALDestroyDict(lalParams_aux);

    return XLAL_SUCCESS;
  }

  /**
   * This is an external wrapper to generate the (l,m) PNR angles,
   * following the prescriptions outlined in arXiv:2107.08876
   * and arXiv:##### FIXME: add reference,
   * given the standard inputs given to generate FD waveforms.
   */
  int XLALSimIMRPhenomX_PNR_GeneratePNRAnglesHM(
      REAL8Sequence **alphaPNR, /**< [out] Alpha precession angle (rad) */
      REAL8Sequence **betaPNR,  /**< [out] Beta precession angle (rad) */
      REAL8Sequence **gammaPNR, /**< [out] Gamma precession angle (rad) */
      REAL8Sequence **freqs,    /**< [out] Frequency array (Hz) */
      REAL8 *alphaPNR_ref,      /**< [out] Reference value of alpha (rad) */
      REAL8 *gammaPNR_ref,      /**< [out] Reference value of gamma (rad) */
      REAL8 m1_SI,              /**< mass of companion 1 (kg) */
      REAL8 m2_SI,              /**< mass of companion 2 (kg) */
      REAL8 chi1x,              /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
      REAL8 chi1y,              /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
      REAL8 chi1z,              /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
      REAL8 chi2x,              /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
      REAL8 chi2y,              /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
      REAL8 chi2z,              /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
      REAL8 inclination,        /**< Angle between orbital angular momentum and line-of-sight vector at reference frequency (rad) */
      REAL8 deltaF,             /**< Frequency spacing (Hz) */
      REAL8 f_min,              /**< Starting GW frequency (Hz) */
      REAL8 f_max,              /**< Ending GW frequency (Hz); Defaults to Mf = 0.3 if no f_max is specified. */
      REAL8 fRef_In,            /**< Reference frequency (Hz) */
      INT4 ell,                 /**< Orbital index (int) */
      INT4 emmprime,            /**< Azimuthal index (int) */
      LALDict *lalParams        /**< LAL Dictionary struct */
  )
  {
    /* Simple check on masses and spins */
    UINT4 status = 0;
    status = XLALIMRPhenomXPCheckMassesAndSpins(&m1_SI, &m2_SI, &chi1x, &chi1y, &chi1z, &chi2x, &chi2y, &chi2z);
    XLAL_CHECK(
        XLAL_SUCCESS == status,
        XLAL_EFUNC,
        "Error: XLALIMRPhenomXPCheckMassesAndSpins failed in XLALSimIMRPhenomX_PNR_GeneratePNRAnglesHM.\n");

    /* Ensure we have a dictionary */
    LALDict *lalParams_aux;
    if (lalParams == NULL)
    {
      lalParams_aux = XLALCreateDict();
    }
    else
    {
      lalParams_aux = XLALDictDuplicate(lalParams);
    }

    /* Make sure we're calling the tuned angles */
    int UsePNR = XLALSimInspiralWaveformParamsLookupPhenomXPNRUseTunedAngles(lalParams_aux);
    if (!UsePNR)
    {
      UsePNR = 1;
      XLALSimInspiralWaveformParamsInsertPhenomXPNRUseTunedAngles(lalParams_aux, UsePNR);
    }

    /* make sure the HM multipoles are activated
     * we need both (2,2) and (l,m) */
    LALValue *ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(lalParams_aux);

    /* If the mode array is empty, populate using a default choice of modes */
    if (ModeArray == NULL)
    {
      ModeArray = XLALSimInspiralCreateModeArray();
    }
    XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 2);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 2, -2);
    XLALSimInspiralModeArrayActivateMode(ModeArray, ell, emmprime);
    XLALSimInspiralModeArrayActivateMode(ModeArray, ell, -emmprime);

    XLALSimInspiralWaveformParamsInsertModeArray(lalParams_aux, ModeArray);
    XLALDestroyValue(ModeArray);

    /* Map fRef to the start frequency if need be, then make sure it's within the frequency range */
    REAL8 fRef = (fRef_In == 0.0) ? f_min : fRef_In;
    XLAL_CHECK(
        ((f_min <= fRef) && (fRef <= f_max)),
        XLAL_EFUNC,
        "Error: fRef needs to be within the specified minimum and maximum frequency values!\n");

    /* Check that the passed deltaF is non-zero */
    XLAL_CHECK(
        deltaF > 0,
        XLAL_EFUNC,
        "Error: deltaF needs to be a positive number!\n");

    /* Generate a uniformly sampled frequency grid of spacing deltaF. */
    /* Frequencies will be set using only the lower and upper bounds that we passed */
    size_t iStart = (size_t)(f_min / deltaF);
    size_t iStop = (size_t)(f_max / deltaF) + 1;

    XLAL_CHECK(
        (iStart <= iStop),
        XLAL_EDOM,
        "Error: the starting frequency index is greater than the stopping index! Please ensure that f_min <= f_max.\n");

    /* Allocate memory for frequency array and terminate if this fails */
    *freqs = XLALCreateREAL8Sequence(iStop - iStart);
    if (!(*freqs))
    {
      XLAL_ERROR(XLAL_EFUNC, "Frequency array allocation failed.");
    }
    *alphaPNR = XLALCreateREAL8Sequence(iStop - iStart);
    if (!(*alphaPNR))
    {
      XLAL_ERROR(XLAL_EFUNC, "Alpha array allocation failed.");
    }
    *betaPNR = XLALCreateREAL8Sequence(iStop - iStart);
    if (!(*betaPNR))
    {
      XLAL_ERROR(XLAL_EFUNC, "Beta array allocation failed.");
    }
    *gammaPNR = XLALCreateREAL8Sequence(iStop - iStart);
    if (!(*gammaPNR))
    {
      XLAL_ERROR(XLAL_EFUNC, "Gamma array allocation failed.");
    }

    /* Populate frequency array */
    for (UINT4 i = iStart; i < iStop; i++)
    {
      (*freqs)->data[i - iStart] = i * deltaF;
    }

    /* Specify arbitrary parameters for the angle generation */
    REAL8 distance = 1.0;
    REAL8 phiRef = 0.0;

    /* Use the true minimum and maximum frequency values */
    REAL8 f_min_eval = (*freqs)->data[0];
    REAL8 f_max_eval = (*freqs)->data[(*freqs)->length - 1];

    /* Initialize PhenomX Waveform struct and check that it initialized correctly */
    IMRPhenomXWaveformStruct *pWF;
    pWF = XLALMalloc(sizeof(IMRPhenomXWaveformStruct));
    status = IMRPhenomXSetWaveformVariables(pWF, m1_SI, m2_SI, chi1z, chi2z, deltaF, fRef, phiRef, f_min_eval, f_max_eval, distance, inclination, lalParams_aux, DEBUG);
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetWaveformVariables failed.\n");

    /* Initialize PhenomX Precession struct and check that it generated successfully */
    IMRPhenomXPrecessionStruct *pPrec;
    pPrec = XLALMalloc(sizeof(IMRPhenomXPrecessionStruct));
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
        DEBUG);
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetPrecessionVariables failed.\n");

    /* See if PNR angles were turned off in pPrec check */
    UsePNR = pPrec->IMRPhenomXPNRUseTunedAngles;

    if(UsePNR){
      /* First generate (2,2)-angle interpolants */
      IMRPhenomX_PNR_angle_spline *hm_angle_spline = (IMRPhenomX_PNR_angle_spline *)XLALMalloc(sizeof(IMRPhenomX_PNR_angle_spline));
      status = IMRPhenomX_PNR_GeneratePNRAngleInterpolants(hm_angle_spline, pWF, pPrec, lalParams_aux);
      XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomX_PNR_GeneratePNRAngleInterpolants failed.\n");

      REAL8 M = pWF->Mtot;

      /* Interpolate the angles */
      /* If the function is called with (l,m)=(2,2) then just interpolate onto the frequency grid*/
      if ((ell == 2) && (emmprime == 2))
      {
        for (size_t i = 0; i < (*freqs)->length; i++)
        {
          (*alphaPNR)->data[i] = gsl_spline_eval(hm_angle_spline->alpha_spline, (*freqs)->data[i], hm_angle_spline->alpha_acc);
          (*betaPNR)->data[i] = gsl_spline_eval(hm_angle_spline->beta_spline, (*freqs)->data[i], hm_angle_spline->beta_acc);
          (*gammaPNR)->data[i] = gsl_spline_eval(hm_angle_spline->gamma_spline, (*freqs)->data[i], hm_angle_spline->gamma_acc);
        }
      }
      else
      { /* Called with l!=2 or m!=2, so we map the (2,2) angles onto a rescaled frequency grid */

        /* Collect the (2,2) and (l,m) ringdown frequencies for use in the frequency map */
        REAL8 Mf_RD_22 = pWF->fRING;
        REAL8 Mf_RD_lm = IMRPhenomXHM_GenerateRingdownFrequency(ell, emmprime, pWF);

        /* Generate interpolation transition frequencies */
        REAL8 Mf_high = 0.0;
        REAL8 Mf_low = 0.0;

        status = IMRPhenomX_PNR_LinearFrequencyMapTransitionFrequencies(&Mf_low, &Mf_high, emmprime, Mf_RD_22, Mf_RD_lm, pPrec);
        XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomX_PNR_LinearFrequencyMapTransitionFrequencies failed.\n");

        UINT4 toggleInspiralScaling = pPrec->PNRInspiralScaling;

  #if DEBUG == 1
        // Save angles into a file
        FILE *fileangle;
        char fileSpecII[40];
        sprintf(fileSpecII, "interpolation_frequencies.dat");

        fileangle = fopen(fileSpecII, "w");

        fprintf(fileangle, "Mf_RD_22: %.16e\n", Mf_RD_22);
        fprintf(fileangle, "Mf_RD_lm: %.16e\n", Mf_RD_lm);

        fprintf(fileangle, "Mf_low: %.16e\n", Mf_low);
        fprintf(fileangle, "Mf_high: %.16e\n", Mf_high);

        fclose(fileangle);
  #endif
        for (size_t i = 0; i < (*freqs)->length; i++)
        {
          REAL8 Mf = XLALSimIMRPhenomXUtilsHztoMf((*freqs)->data[i], M);
          /* Calculate the mapped frequency */
          REAL8 Mf_mapped = IMRPhenomX_PNR_LinearFrequencyMap(Mf, ell, emmprime, Mf_low, Mf_high, Mf_RD_22, Mf_RD_lm, toggleInspiralScaling);
          REAL8 f_mapped = XLALSimIMRPhenomXUtilsMftoHz(Mf_mapped, M);

          (*alphaPNR)->data[i] = gsl_spline_eval(hm_angle_spline->alpha_spline, f_mapped, hm_angle_spline->alpha_acc);
          (*betaPNR)->data[i] = gsl_spline_eval(hm_angle_spline->beta_spline, f_mapped, hm_angle_spline->beta_acc);
          (*gammaPNR)->data[i] = gsl_spline_eval(hm_angle_spline->gamma_spline, f_mapped, hm_angle_spline->gamma_acc);
        }
      }

      /* Here we assign the reference values of alpha and gamma to their values in the precession struct */
      /* NOTE: the contribution from pPrec->alpha0 is assigned in IMRPhenomX_PNR_RemapThetaJSF */
      pPrec->alpha_offset = gsl_spline_eval(hm_angle_spline->alpha_spline, pWF->fRef, hm_angle_spline->alpha_acc);
      /* NOTE: the sign is flipped between gamma and epsilon */
      pPrec->epsilon_offset = -gsl_spline_eval(hm_angle_spline->gamma_spline, pWF->fRef, hm_angle_spline->gamma_acc) - pPrec->epsilon0; // note the sign difference between gamma and epsilon

      /* Remap the J-frame sky location to use beta instead of ThetaJN */
      REAL8 betaPNR_ref = gsl_spline_eval(hm_angle_spline->beta_spline, pWF->fRef, hm_angle_spline->beta_acc);
      status = IMRPhenomX_PNR_RemapThetaJSF(betaPNR_ref, pWF, pPrec, lalParams_aux);
      XLAL_CHECK(
          XLAL_SUCCESS == status,
          XLAL_EFUNC,
          "Error: IMRPhenomX_PNR_RemapThetaJSF failed in IMRPhenomX_PNR_GeneratePNRAnglesHM.");

      /* Pass through the reference values of alpha and gamma */
      *alphaPNR_ref = pPrec->alpha_offset;
      *gammaPNR_ref = -(pPrec->epsilon_offset);

      gsl_spline_free(hm_angle_spline->alpha_spline);
      gsl_spline_free(hm_angle_spline->beta_spline);
      gsl_spline_free(hm_angle_spline->gamma_spline);

      gsl_interp_accel_free(hm_angle_spline->alpha_acc);
      gsl_interp_accel_free(hm_angle_spline->beta_acc);
      gsl_interp_accel_free(hm_angle_spline->gamma_acc);

      LALFree(hm_angle_spline);
    }
    else{
      /* Generate PN angles */
      if(XLALSimInspiralWaveformParamsLookupPhenomXPrecVersion(lalParams_aux) > 200 ){
        /* MSA angles generated using built-in function */

        REAL8Sequence *cosBeta = XLALCreateREAL8Sequence(iStop - iStart);
        if (!cosBeta)
        {
          XLAL_ERROR(XLAL_EFUNC, "cosBeta array allocation failed.");
        }

        status = XLALSimIMRPhenomXPMSAAngles(
          alphaPNR, gammaPNR, &cosBeta, *freqs,
          m1_SI,m2_SI,
          chi1x,chi1y,chi1z,
          chi2x,chi2y,chi2z,
          inclination,
          fRef_In, emmprime, lalParams_aux
        );

        for (UINT4 i = 0; i < (*freqs)->length; i++)
        {
          (*betaPNR)->data[i] = acos(cosBeta->data[i]);
        }

        XLALDestroyREAL8Sequence(cosBeta);

        /* Reference values already applied in function */
        *alphaPNR_ref = 0.0;
        *gammaPNR_ref = 0.0;
      }
      else{
        /* NNLO angles */

        REAL8Sequence *cosBeta = XLALCreateREAL8Sequence(iStop - iStart);
        if (!cosBeta)
        {
          XLAL_ERROR(XLAL_EFUNC, "cosBeta array allocation failed.");
        }

        status = XLALSimIMRPhenomXPPNAngles(
          alphaPNR, gammaPNR, &cosBeta, *freqs,
          m1_SI,m2_SI,
          chi1x,chi1y,chi1z,
          chi2x,chi2y,chi2z,
          inclination,
          fRef_In, emmprime, lalParams_aux
        );

        for (UINT4 i = 0; i < (*freqs)->length; i++)
        {
          (*betaPNR)->data[i] = acos(cosBeta->data[i]);
        }

        XLALDestroyREAL8Sequence(cosBeta);

        /* Pass through the reference values of alpha and gamma
        * NOTE: gamma = - epsilon */
        *alphaPNR_ref = pPrec->alpha_offset;
        *gammaPNR_ref = -(pPrec->epsilon_offset);
      }
    }

    /* Clean up memory allocation */
    LALFree(pPrec);
    LALFree(pWF);
    XLALDestroyDict(lalParams_aux);

    return XLAL_SUCCESS;
  }

  /**
   * Generate the tuned precession angles outlined in arXiv:2107.08876.
   * The general flow is as follows:
   *
   * - First check if deltaF > 0:
   *
   * --> if it is, then we generate angles sampled on a uniform frequency grid
   *     and then compute the values of the angles at fRef using linear interpolation.
   *
   * --> otherwise we expect non-uniform frequencies, so we generate interpolants
   *     and then evaluate them at fRef.
   *
   * - The reference values of alpha and gamma are stored in the pPrec struct for use later,
   *   and then we re-map the J-frame sky position using the reference value of beta.
   *
   * - NOTE: alpha0 is recomputed in IMRPhenomX_PNR_RemapThetaJSF, hence it is not used here
   *   unlike epsilon0. Instead it is applied inside IMRPhenomX_PNR_RemapThetaJSF.
   *
   */
  int IMRPhenomX_PNR_GeneratePNRAngles(
      REAL8Sequence *alphaPNR,           /**< alpha precession angle (rad) */
      REAL8Sequence *betaPNR,            /**< beta precession angle (rad) */
      REAL8Sequence *gammaPNR,           /**< gamma precession angle (rad) */
      const REAL8Sequence *freqs,        /**< input frequency array (Hz) */
      IMRPhenomXWaveformStruct *pWF,     /**< waveform struct */
      IMRPhenomXPrecessionStruct *pPrec, /**< precession struct */
      LALDict *lalParams                 /**< LAL dictionary struct */
  )
  {
    /* Make sure the incoming pointers lead to something initialized */
    XLAL_CHECK(alphaPNR != NULL, XLAL_EFAULT);
    XLAL_CHECK(betaPNR != NULL, XLAL_EFAULT);
    XLAL_CHECK(gammaPNR != NULL, XLAL_EFAULT);
    XLAL_CHECK(freqs != NULL, XLAL_EFAULT);
    XLAL_CHECK(pWF != NULL, XLAL_EFAULT);
    XLAL_CHECK(pPrec != NULL, XLAL_EFAULT);
    XLAL_CHECK(lalParams != NULL, XLAL_EFAULT);

    REAL8 f_ref = pWF->fRef;
    REAL8 deltaF = pWF->deltaF;

    /* Make sure we're supposed to be here */
    int UsePNR = pPrec->IMRPhenomXPNRUseTunedAngles;
    XLAL_CHECK(
        UsePNR,
        XLAL_EFUNC,
        "Error: PNR angles called without being activated!\n");

    INT4 status;

    /* Check for uniform frequency series */
    if (deltaF > 0.0)
    {

      /* Generate PNR structs */
      IMRPhenomXWaveformStruct *pWF_SingleSpin = NULL;
      IMRPhenomXPrecessionStruct *pPrec_SingleSpin = NULL;
      IMRPhenomX_PNR_alpha_parameters *alphaParams = NULL;
      IMRPhenomX_PNR_beta_parameters *betaParams = NULL;

      status = IMRPhenomX_PNR_PopulateStructs(
        &pWF_SingleSpin,
        &pPrec_SingleSpin,
        &alphaParams,
        &betaParams,
        pWF,
        pPrec,
        lalParams);
      XLAL_CHECK(
        XLAL_SUCCESS == status,
        XLAL_EFUNC,
        "Error: IMRPhenomX_PNR_PopulateStructs failed!\n");

      /* generate angles */
      status = IMRPhenomX_PNR_GeneratePNRAngles_UniformFrequencies(
          alphaPNR,
          betaPNR,
          gammaPNR,
          freqs,
          pWF_SingleSpin,
          pPrec_SingleSpin,
          alphaParams,
          betaParams,
          pWF,
          pPrec,
          lalParams);
      XLAL_CHECK(
          XLAL_SUCCESS == status,
          XLAL_EFUNC,
          "Error: IMRPhenomX_PNR_GeneratePNRAngles_UniformFrequencies failed in IMRPhenomX_PNR_GeneratePNRAngles.");

      /* Here we assign the reference values of alpha and gamma to their values in the precession struct.
       * For the uniformly-sampled arrays, we linearly interpolate the angles between the
       * next-to-nearest frequencies to get the reference values. */

      /* NOTE: the contribution from pPrec->alpha0 is assigned in IMRPhenomX_PNR_RemapThetaJSF */
      pPrec->alpha_offset = IMRPhenomX_PNR_AngleAtFRef(alphaPNR, f_ref, freqs, deltaF);
      /* NOTE: the sign is flipped between gamma and epsilon */
      pPrec->epsilon_offset = -IMRPhenomX_PNR_AngleAtFRef(gammaPNR, f_ref, freqs, deltaF) - pPrec->epsilon0;

      /* Remap the J-frame sky location to use beta instead of ThetaJN */
      REAL8 betaPNR_ref = IMRPhenomX_PNR_AngleAtFRef(betaPNR, f_ref, freqs, deltaF);
      status = IMRPhenomX_PNR_RemapThetaJSF(betaPNR_ref, pWF, pPrec, lalParams);
      XLAL_CHECK(
          XLAL_SUCCESS == status,
          XLAL_EFUNC,
          "Error: IMRPhenomX_PNR_RemapThetaJSF failed in IMRPhenomX_PNR_GeneratePNRAngles.");

      /* clean this up */
      IMRPhenomX_PNR_FreeStructs(
        &pWF_SingleSpin,
        &pPrec_SingleSpin,
        &alphaParams,
        &betaParams);
    }
    else
    { /* Non-uniform frequency array, so we generate interpolants to evaluate */
      IMRPhenomX_PNR_angle_spline *angle_spline = (IMRPhenomX_PNR_angle_spline *)XLALMalloc(sizeof(IMRPhenomX_PNR_angle_spline));
      /* Generate the angle interpolants */
      status = IMRPhenomX_PNR_GeneratePNRAngleInterpolants(
          angle_spline,
          pWF, pPrec,
          lalParams);
      XLAL_CHECK(
          XLAL_SUCCESS == status,
          XLAL_EFUNC,
          "Error: IMRPhenomX_PNR_GeneratePNRAngleInterpolants failed in IMRPhenomX_PNR_GeneratePNRAngles");

      /* Fill the angle arrays with the values of the angles at the frequncy points */
      REAL8 fval = 0;
      for (UINT4 i = 0; i < freqs->length; i++)
      {
        fval = freqs->data[i];
        alphaPNR->data[i] = gsl_spline_eval(angle_spline->alpha_spline, fval, angle_spline->alpha_acc);
        betaPNR->data[i] = gsl_spline_eval(angle_spline->beta_spline, fval, angle_spline->beta_acc);
        gammaPNR->data[i] = gsl_spline_eval(angle_spline->gamma_spline, fval, angle_spline->gamma_acc);
      }

      /* Here we assign the reference values of alpha and gamma to their values in the precession struct */
      /* NOTE: the contribution from pPrec->alpha0 is assigned in IMRPhenomX_PNR_RemapThetaJSF */
      pPrec->alpha_offset = gsl_spline_eval(angle_spline->alpha_spline, f_ref, angle_spline->alpha_acc);
      /* NOTE: the sign is flipped between gamma and epsilon */
      pPrec->epsilon_offset = -gsl_spline_eval(angle_spline->gamma_spline, f_ref, angle_spline->gamma_acc) - pPrec->epsilon0;

      /* Remap the J-frame sky location to use beta instead of ThetaJN */
      REAL8 betaPNR_ref = gsl_spline_eval(angle_spline->beta_spline, f_ref, angle_spline->beta_acc);
      status = IMRPhenomX_PNR_RemapThetaJSF(betaPNR_ref, pWF, pPrec, lalParams);
      XLAL_CHECK(
          XLAL_SUCCESS == status,
          XLAL_EFUNC,
          "Error: IMRPhenomX_PNR_RemapThetaJSF failed in IMRPhenomX_PNR_GeneratePNRAngles.");

      /* Free up memory */
      gsl_spline_free(angle_spline->alpha_spline);
      gsl_spline_free(angle_spline->beta_spline);
      gsl_spline_free(angle_spline->gamma_spline);

      gsl_interp_accel_free(angle_spline->alpha_acc);
      gsl_interp_accel_free(angle_spline->beta_acc);
      gsl_interp_accel_free(angle_spline->gamma_acc);

      LALFree(angle_spline);
    }

    return XLAL_SUCCESS;
  }

  /**
   * Generate the tuned precession angles outlined in arXiv:2107.08876
   * specifically on a uniform frequency grid.
   *
   * - First, we populate the required alpha and beta parameter structs,
   *   along with the two-spin structs if needed.
   *
   * - Next we store the two connection frequencies possibly needed for the
   *   HM angle frequency mapping.
   *
   * - Check if we should be attaching the MR contributions to beta: this is
   *   determined by calling the function IMRPhenomX_PNR_AttachMRBeta.
   *
   * --> If yes, loop through the frequencies and generate full alpha and beta
   *
   * --> If no, loop through the frequencies and generate full alpha but only
   *     the inspiral beta
   *
   * - Finally, if an initialized struct is passed for gamma, compute it.
   *
   */
  int IMRPhenomX_PNR_GeneratePNRAngles_UniformFrequencies(
      REAL8Sequence *alphaPNR,
      REAL8Sequence *betaPNR,
      REAL8Sequence *gammaPNR,
      const REAL8Sequence *freqs,
      IMRPhenomXWaveformStruct *pWF_SingleSpin,
      IMRPhenomXPrecessionStruct *pPrec_SingleSpin,
      IMRPhenomX_PNR_alpha_parameters *alphaParams,
      IMRPhenomX_PNR_beta_parameters *betaParams,
      IMRPhenomXWaveformStruct *pWF,
      IMRPhenomXPrecessionStruct *pPrec,
      LALDict *lalParams)
  {
    /* Make sure the incoming pointers lead to something initialized */
    /* Gamma could be a NULL pointer, so we don't check here */
    XLAL_CHECK(alphaPNR != NULL, XLAL_EFAULT);
    XLAL_CHECK(betaPNR != NULL, XLAL_EFAULT);
    XLAL_CHECK(freqs != NULL, XLAL_EFAULT);
    XLAL_CHECK(pWF != NULL, XLAL_EFAULT);
    XLAL_CHECK(pPrec != NULL, XLAL_EFAULT);
    XLAL_CHECK(lalParams != NULL, XLAL_EFAULT);

    REAL8 M = pWF->Mtot;

    /* Make sure we're supposed to be here */
    int UsePNR = pPrec->IMRPhenomXPNRUseTunedAngles;
    XLAL_CHECK(
        UsePNR,
        XLAL_EFUNC,
        "Error: PNR angles called without being activated!\n");

    INT4 status;

   #if DEBUG == 1
    IMRPhenomX_PNR_AngleParameterDebugPrint(alphaParams, betaParams);
   #endif

    REAL8 Mf = 0.0;
    /* generate PNR angles */
    REAL8 q = pWF->q;
    REAL8 chi = pPrec->chi_singleSpin;
    /* inside calibration region */
    if ((q <= pPrec->PNR_q_window_lower) && (chi <= pPrec->PNR_chi_window_lower))
    {
      /* First check to see if we attach the MR tuning to beta */
      UINT4 attach_MR_beta = IMRPhenomX_PNR_AttachMRBeta(betaParams);
      if (attach_MR_beta) /* yes we do! */
	{
	  for (size_t i = 0; i < freqs->length; i++)
	    {
	      Mf = XLALSimIMRPhenomXUtilsHztoMf(freqs->data[i], M);
	      /* generate alpha and beta with MR tuning */
	      alphaPNR->data[i] = IMRPhenomX_PNR_GeneratePNRAlphaAtMf(Mf, alphaParams, pWF, pPrec);
	      betaPNR->data[i] = IMRPhenomX_PNR_GeneratePNRBetaAtMf(Mf, betaParams, pWF, pPrec, pWF_SingleSpin, pPrec_SingleSpin);
	    }
	}
      else /* don't attach MR tuning to beta */
	{
	  for (size_t i = 0; i < freqs->length; i++)
	    {
	      Mf = XLALSimIMRPhenomXUtilsHztoMf(freqs->data[i], M);

	      /* generate alpha, generate beta with no MR tuning */
	      alphaPNR->data[i] = IMRPhenomX_PNR_GeneratePNRAlphaAtMf(Mf, alphaParams, pWF, pPrec);
	      betaPNR->data[i] = IMRPhenomX_PNR_GeneratePNRBetaNoMR(Mf, pWF, pPrec);
	    }
	}
    }
    /* inside transition region */
    else if ((q <= pPrec->PNR_q_window_upper) && (chi <= pPrec->PNR_chi_window_upper))
      {
      /* First check to see if we attach the MR tuning to beta */
      UINT4 attach_MR_beta = IMRPhenomX_PNR_AttachMRBeta(betaParams);
      if (attach_MR_beta) /* yes we do! */
        {
          for (size_t i = 0; i < freqs->length; i++)
            {
              Mf = XLALSimIMRPhenomXUtilsHztoMf(freqs->data[i], M);
              /* generate alpha and beta with MR tuning */
              alphaPNR->data[i] = IMRPhenomX_PNR_GenerateMergedPNRAlphaAtMf(Mf, alphaParams, pWF, pPrec);
	      betaPNR->data[i] = IMRPhenomX_PNR_GenerateMergedPNRBetaAtMf(Mf, betaParams, pWF, pPrec, pWF_SingleSpin, pPrec_SingleSpin);
            }
        }
      else /* don't attach MR tuning to beta */
        {
          for (size_t i = 0; i < freqs->length; i++)
            {
              Mf = XLALSimIMRPhenomXUtilsHztoMf(freqs->data[i], M);

              /* generate alpha, generate beta with no MR tuning */
              alphaPNR->data[i] = IMRPhenomX_PNR_GenerateMergedPNRAlphaAtMf(Mf, alphaParams, pWF, pPrec);
              betaPNR->data[i] = IMRPhenomX_PNR_GeneratePNRBetaNoMR(Mf, pWF, pPrec);
            }
        }

      }
    /* fully in outside calibration region */
    else{
      for (size_t i = 0; i < freqs->length; i++)
            {
              Mf = XLALSimIMRPhenomXUtilsHztoMf(freqs->data[i], M);

              /* generate MSA alpha, generate beta with no MR tuning */
              alphaPNR->data[i] = IMRPhenomX_PNR_GetPNAlphaAtFreq(Mf, pWF, pPrec);
              betaPNR->data[i] = IMRPhenomX_PNR_GeneratePNRBetaNoMR(Mf, pWF, pPrec);
            }
    }

    /* If gamma is a NULL pointer, we assume that we don't want to
     * compute it. This allows us to skip computing gamma for the
     * interpolation code to reduce redundant steps */
    if (gammaPNR != NULL)
      {
	/* generate gamma */
	status = IMRPhenomX_PNR_GeneratePNRGamma(gammaPNR, freqs, alphaPNR, betaPNR);
	XLAL_CHECK(
		   XLAL_SUCCESS == status,
		   XLAL_EFUNC,
		   "Error: IMRPhenomX_PNR_GeneratePNRGamma failed");
      }

    return XLAL_SUCCESS;
  }

    /**
   * Internal helper function to generate PNR alpha for the
   * antisymmetric waveform.
   */
  int IMRPhenomX_PNR_GeneratePNRAlphaForAntisymmetry(
      REAL8Sequence *alphaPNR,
      const REAL8Sequence *freqs,
      IMRPhenomXWaveformStruct *pWF,
      IMRPhenomXPrecessionStruct *pPrec,
      LALDict *lalParams)
  {
    /* Make sure the incoming pointers lead to something initialized */
    /* Gamma could be a NULL pointer, so we don't check here */
    XLAL_CHECK(alphaPNR != NULL, XLAL_EFAULT);
    XLAL_CHECK(freqs != NULL, XLAL_EFAULT);
    XLAL_CHECK(pWF != NULL, XLAL_EFAULT);
    XLAL_CHECK(pPrec != NULL, XLAL_EFAULT);
    XLAL_CHECK(lalParams != NULL, XLAL_EFAULT);

    REAL8 M = pWF->Mtot;

    INT4 status;

    /* generate alpha parameters */
    IMRPhenomX_PNR_alpha_parameters *alphaParams = XLALMalloc(sizeof(IMRPhenomX_PNR_alpha_parameters));
    status = IMRPhenomX_PNR_precompute_alpha_coefficients(alphaParams, pWF, pPrec);
    XLAL_CHECK(
        XLAL_SUCCESS == status,
        XLAL_EFUNC,
        "Error: IMRPhenomX_PNR_precompute_alpha_coefficients failed.\n");

    status = IMRPhenomX_PNR_alpha_connection_parameters(alphaParams, pWF, pPrec);
    XLAL_CHECK(
        XLAL_SUCCESS == status,
        XLAL_EFUNC,
        "Error: IMRPhenomX_PNR_alpha_connection_parameters failed.\n");

    REAL8 Mf = 0.0;
    /* generate PNR angles */
    REAL8 q = pWF->q;
    REAL8 chi = pPrec->chi_singleSpin;
    /* inside calibration region */
    if ((q <= pPrec->PNR_q_window_lower) && (chi <= pPrec->PNR_chi_window_lower))
    {
      for (size_t i = 0; i < freqs->length; i++)
      {
        Mf = XLALSimIMRPhenomXUtilsHztoMf(freqs->data[i], M);

        /* generate alpha */
        alphaPNR->data[i] = IMRPhenomX_PNR_GeneratePNRAlphaAtMf(Mf, alphaParams, pWF, pPrec);
      }
    }
    /* inside transition region */
    else if ((q <= pPrec->PNR_q_window_upper) && (chi <= pPrec->PNR_chi_window_upper))
    {
      for (size_t i = 0; i < freqs->length; i++)
      {
        Mf = XLALSimIMRPhenomXUtilsHztoMf(freqs->data[i], M);

        /* generate alpha */
        alphaPNR->data[i] = IMRPhenomX_PNR_GenerateMergedPNRAlphaAtMf(Mf, alphaParams, pWF, pPrec);
      }
    }
    /* fully in outside calibration region */
    else
    {
      for (size_t i = 0; i < freqs->length; i++)
      {
        Mf = XLALSimIMRPhenomXUtilsHztoMf(freqs->data[i], M);

        /* generate MSA or NNLO alpha */
        alphaPNR->data[i] = IMRPhenomX_PNR_GetPNAlphaAtFreq(Mf, pWF, pPrec);
      }
    }

    LALFree(alphaParams);

    return XLAL_SUCCESS;
  }

  /**
   * Generate the tuned precession angles outlined in arXiv:2107.08876
   * specifically populate the angle interpolant struct.
   *
   * - First, we check for an activated mode array.
   *
   * --> If it exists, then we're likely computing HM angles with these interpolants
   *     and we need to extend the frequency range to capture the HM frequency
   *     map. Use the activated modes to determine this.
   *
   * --> Otherwise we use the default f_min and f_max.
   *
   * - Get the required deltaF from IMRPhenomX_PNR_HMInterpolationDeltaF.
   *
   * - Compute uniform alpha and beta
   *
   * - Populate interpolants for alpha and beta, then compute gamma using these interpolants
   *
   */
  int IMRPhenomX_PNR_GeneratePNRAngleInterpolants(
      IMRPhenomX_PNR_angle_spline *angle_spline,
      IMRPhenomXWaveformStruct *pWF,
      IMRPhenomXPrecessionStruct *pPrec,
      LALDict *lalparams)
  {
    /* check for initialization */
    XLAL_CHECK(angle_spline != NULL, XLAL_EFAULT);
    XLAL_CHECK(pWF != NULL, XLAL_EFAULT);
    XLAL_CHECK(pPrec != NULL, XLAL_EFAULT);
    XLAL_CHECK(lalparams != NULL, XLAL_EFAULT);

    /* Generate PNR structs */
    IMRPhenomXWaveformStruct *pWF_SingleSpin = NULL;
    IMRPhenomXPrecessionStruct *pPrec_SingleSpin = NULL;
    IMRPhenomX_PNR_alpha_parameters *alphaParams = NULL;
    IMRPhenomX_PNR_beta_parameters *betaParams = NULL;

    UINT4 status = 0;
    status = IMRPhenomX_PNR_PopulateStructs(
        &pWF_SingleSpin,
        &pPrec_SingleSpin,
        &alphaParams,
        &betaParams,
        pWF,
        pPrec,
        lalparams);
    XLAL_CHECK(
        XLAL_SUCCESS == status,
        XLAL_EFUNC,
        "Error: IMRPhenomX_PNR_PopulateStructs failed!\n");

    /* Store connection frequencies */
    pPrec->PNR_HM_Mfhigh = betaParams->Mf_beta_lower;
    pPrec->PNR_HM_Mflow = alphaParams->Mf_alpha_lower;

    /* Check for connection frequencies potentially being poorly defined for frequency interpolation */
    /* Check that the lower and upper connection frequencies are well behaved relative to fCut, fRING, and each other.
     * fCutDef is Mf = 0.3 is except for large chiEff. */
    if((pPrec->PNR_HM_Mfhigh > pWF->fCutDef) || (pPrec->PNR_HM_Mfhigh < 0.1 * pWF->fRING)){
      pPrec->PNR_HM_Mfhigh = pWF->fRING;
    }
    if((pPrec->PNR_HM_Mflow > pWF->fCutDef) || (pPrec->PNR_HM_Mfhigh < pPrec->PNR_HM_Mflow)){
      pPrec->PNR_HM_Mflow = pPrec->PNR_HM_Mfhigh / 2.0;
    }

    REAL8 f_min_22 = pWF->fMin;
    REAL8 f_max_22 = pWF->f_max_prime;

    if(pWF->deltaF != 0){
      /* need to account for finite frequency spacing in fMin */
      size_t iStart = (size_t) (f_min_22 / pWF->deltaF);

      REAL8 fMinFromDeltaF = iStart*pWF->deltaF;
      f_min_22 = (fMinFromDeltaF < f_min_22) ? fMinFromDeltaF : f_min_22;
    }

    /* Grab the mode array to find the minimum and maximum
     * frequency values that need to be interpolated. */
    LALValue *ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(lalparams);

    if (ModeArray != NULL)
    {
      /* Assume we possibly have higher multipoles. */
      for (UINT4 ell = 2; ell <= 4; ell++)
      {
        for (UINT4 emmprime = 1; emmprime <= ell; emmprime++)
        {
          /* First check for (2,2) and multibanding */
          if(XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, emmprime) && ((ell == 2)&&(emmprime == 2)))
          {
            if(pPrec->MBandPrecVersion != 0)
            {
              /* For (2,2) with angle MB, find max (2,2) frequency */
              REAL8Sequence *coarseFreqs;
              XLALSimIMRPhenomXPHMMultibandingGrid(&coarseFreqs, ell, emmprime, pWF, lalparams);

              REAL8 MfCoarseFMax = coarseFreqs->data[coarseFreqs->length-1];
              REAL8 coarseFMax = XLALSimIMRPhenomXUtilsMftoHz(MfCoarseFMax, pWF->Mtot);
              if(coarseFMax > f_max_22){
                f_max_22 = coarseFMax;
              }

              XLALDestroyREAL8Sequence(coarseFreqs);
            }
          }
          else if(XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, emmprime) && !((ell == 2)&&(emmprime == 2)))
          {
            if((pWF->q == 1) && (pWF->chi1L == pWF->chi2L) && (emmprime % 2 != 0))
            {
              continue;
            }
            REAL8 Mf_RD_22 = pWF->fRING;
            REAL8 Mf_RD_lm = IMRPhenomXHM_GenerateRingdownFrequency(ell, emmprime, pWF);

            REAL8 Mf_low = 0.0;
            REAL8 Mf_high = 0.0;

            IMRPhenomX_PNR_LinearFrequencyMapTransitionFrequencies(&Mf_low, &Mf_high, emmprime, Mf_RD_22, Mf_RD_lm, pPrec);

            UINT4 toggleInspiralScaling = pPrec->PNRInspiralScaling;

            if(pPrec->MBandPrecVersion != 0)
            {
              /* For (l,m) with angle MB, find max mapped (2,2) frequency */
              REAL8Sequence *coarseFreqs;
              XLALSimIMRPhenomXPHMMultibandingGrid(&coarseFreqs, ell, emmprime, pWF, lalparams);

              REAL8 MfCoarseFMin = coarseFreqs->data[0];
              REAL8 MfCoarseFMax = coarseFreqs->data[coarseFreqs->length-1];

              REAL8 MfMinMapped = IMRPhenomX_PNR_LinearFrequencyMap(MfCoarseFMin, ell, emmprime, Mf_low, Mf_high, Mf_RD_22, Mf_RD_lm, toggleInspiralScaling);
              REAL8 MfMaxMapped = IMRPhenomX_PNR_LinearFrequencyMap(MfCoarseFMax, ell, emmprime, Mf_low, Mf_high, Mf_RD_22, Mf_RD_lm, toggleInspiralScaling);

              REAL8 coarseFMin = XLALSimIMRPhenomXUtilsMftoHz(MfMinMapped, pWF->Mtot);
              REAL8 coarseFMax = XLALSimIMRPhenomXUtilsMftoHz(MfMaxMapped, pWF->Mtot);

              if(coarseFMin < f_min_22)
              {
                f_min_22 = coarseFMin;
              }

              if(coarseFMax > f_max_22)
              {
                f_max_22 = coarseFMax;
              }

              XLALDestroyREAL8Sequence(coarseFreqs);
            }
            else
            {
              REAL8 MfFMin = XLALSimIMRPhenomXUtilsHztoMf(f_min_22,pWF->Mtot);
              REAL8 MfFMax = XLALSimIMRPhenomXUtilsHztoMf(f_max_22,pWF->Mtot);

              REAL8 MfMinMapped = IMRPhenomX_PNR_LinearFrequencyMap(MfFMin, ell, emmprime, Mf_low, Mf_high, Mf_RD_22, Mf_RD_lm, toggleInspiralScaling);
              REAL8 MfMaxMapped = IMRPhenomX_PNR_LinearFrequencyMap(MfFMax, ell, emmprime, Mf_low, Mf_high, Mf_RD_22, Mf_RD_lm, toggleInspiralScaling);

              REAL8 FMin = XLALSimIMRPhenomXUtilsMftoHz(MfMinMapped, pWF->Mtot);
              REAL8 FMax = XLALSimIMRPhenomXUtilsMftoHz(MfMaxMapped, pWF->Mtot);

              if(FMin < f_min_22)
              {
                f_min_22 = FMin;
              }

              if(FMax > f_max_22)
              {
                f_max_22 = FMax;
              }
            }
          }
        }
      }
    } /* Otherwise f_min and f_max are already set to the (2,2) values */

    REAL8 f_min = f_min_22;
    REAL8 f_max = f_max_22;

    /* Get appropriate frequency spacing */
    REAL8 deltaF_int = IMRPhenomX_PNR_HMInterpolationDeltaF(f_min, pWF, pPrec);

    /* guarantee additional space to ensure no extrapolation and to reduce effects
     * of natural boundary conditions */
    f_min = (f_min - 2.0 * deltaF_int < 0) ? f_min / 2.0 : f_min - 2.0 * deltaF_int;
    f_max += 2.0 * deltaF_int;

    /* Frequencies will be set using the lower and upper bounds and computed deltaF */
    size_t iStart = (size_t)(f_min / deltaF_int);
    size_t iStop = (size_t)(f_max / deltaF_int) + 1;

#if DEBUG == 1
    /* Save interpolation parameters into a file */
    FILE *fileangle;
    char fileSpecII[40];
    sprintf(fileSpecII, "interpolation_parameters.dat");

    fileangle = fopen(fileSpecII, "w");

    fprintf(fileangle, "f_min_22: %.16e\n", f_min_22);
    fprintf(fileangle, "f_max_22: %.16e\n", f_max_22);
    fprintf(fileangle, "f_min: %.16e\n", f_min);
    fprintf(fileangle, "f_max: %.16e\n", f_max);
    fprintf(fileangle, "fCut: %.16e\n", pWF->fCut);
    fprintf(fileangle, "f_min from df: %.16e\n", iStart * deltaF_int);
    fprintf(fileangle, "f_max from df: %.16e\n", iStop * deltaF_int);
    fprintf(fileangle, "deltaF: %.16e\n\n", deltaF_int);

    fclose(fileangle);
#endif

    XLAL_CHECK((iStart <= iStop), XLAL_EDOM,
               "minimum freq index %zu and maximum freq index %zu do not fulfill 0<=ind_min<=ind_max.\n", iStart, iStop);

    /* Allocate memory for frequency array and angles, and terminate if this fails */
    REAL8Sequence *alphaPNR_22 = NULL;
    REAL8Sequence *betaPNR_22 = NULL;
    REAL8Sequence *gammaPNR_22 = NULL;
    REAL8Sequence *freqs_22 = NULL;

    freqs_22 = XLALCreateREAL8Sequence(iStop - iStart);
    if (!freqs_22)
    {
      XLAL_ERROR(XLAL_EFUNC, "freqs_22 array allocation failed.");
    }
    alphaPNR_22 = XLALCreateREAL8Sequence(iStop - iStart);
    if (!alphaPNR_22)
    {
      XLAL_ERROR(XLAL_EFUNC, "alphaPNR_22 array allocation failed.");
    }
    betaPNR_22 = XLALCreateREAL8Sequence(iStop - iStart);
    if (!betaPNR_22)
    {
      XLAL_ERROR(XLAL_EFUNC, "gammaPNR_22 array allocation failed.");
    }
    gammaPNR_22 = XLALCreateREAL8Sequence(iStop - iStart);
    if (!gammaPNR_22)
    {
      XLAL_ERROR(XLAL_EFUNC, "gammaPNR_22 array allocation failed.");
    }

    /* Populate frequency array */
    for (UINT4 i = iStart; i < iStop; i++)
    {
      freqs_22->data[i - iStart] = i * deltaF_int;
    }

    /* Generate the (2,2) angles on this uniform frequency grid.
     * Pass gamma as a NULL pointer, since we can save a step by
     * re-using the gsl splines constructed to compute gamma later on */
    status = IMRPhenomX_PNR_GeneratePNRAngles_UniformFrequencies(
      alphaPNR_22,
      betaPNR_22,
      NULL,
      freqs_22,
      pWF_SingleSpin,
      pPrec_SingleSpin,
      alphaParams,
      betaParams,
      pWF,
      pPrec,
      lalparams);
    XLAL_CHECK(
        XLAL_SUCCESS == status,
        XLAL_EFUNC,
        "Error: IMRPhenomX_PNR_GeneratePNRAngles_UniformFrequencies failed in IMRPhenomX_PNR_GeneratePNRAngleInterpolants.\n");

    /* Construct the splines for interpolation and populate the alpha and beta splines */
    UINT4 flen = freqs_22->length;

    gsl_spline *a_spline = gsl_spline_alloc(gsl_interp_cspline, flen);
    gsl_spline *b_spline = gsl_spline_alloc(gsl_interp_cspline, flen);
    gsl_spline *g_spline = gsl_spline_alloc(gsl_interp_cspline, flen);

    gsl_interp_accel *a_acc = gsl_interp_accel_alloc();
    gsl_interp_accel *b_acc = gsl_interp_accel_alloc();
    gsl_interp_accel *g_acc = gsl_interp_accel_alloc();

    gsl_spline_init(a_spline, freqs_22->data, alphaPNR_22->data, flen);
    gsl_spline_init(b_spline, freqs_22->data, betaPNR_22->data, flen);

    angle_spline->alpha_spline = a_spline;
    angle_spline->beta_spline = b_spline;

    angle_spline->alpha_acc = a_acc;
    angle_spline->beta_acc = b_acc;

    /* Now use the alpha and beta splines to compute gamma */
    status = IMRPhenomX_PNR_GeneratePNRGamma_FromInterpolants(
        gammaPNR_22,
        freqs_22,
        angle_spline);
    XLAL_CHECK(
        XLAL_SUCCESS == status,
        XLAL_EFUNC,
        "Error: IMRPhenomX_PNR_GeneratePNRGamma_FromInterpolants failed in IMRPhenomX_PNR_GeneratePNRAngleInterpolants.\n");

    /* Reset the gsl accelerators for alpha and beta
     * Not sure if this is helpful =\_(~~)_/= */
    gsl_interp_accel_reset(angle_spline->alpha_acc);
    gsl_interp_accel_reset(angle_spline->beta_acc);

    /* Populate the gamma spline */
    gsl_spline_init(g_spline, freqs_22->data, gammaPNR_22->data, flen);

    angle_spline->gamma_spline = g_spline;
    angle_spline->gamma_acc = g_acc;

    /* Clean up memory allocation
     * In particular, the splines copy over the data
     * so we can get rid of the 22 angle arrays */
    XLALDestroyValue(ModeArray);

    XLALDestroyREAL8Sequence(alphaPNR_22);
    XLALDestroyREAL8Sequence(betaPNR_22);
    XLALDestroyREAL8Sequence(gammaPNR_22);
    XLALDestroyREAL8Sequence(freqs_22);

        /* clean this up */
    IMRPhenomX_PNR_FreeStructs(
        &pWF_SingleSpin,
        &pPrec_SingleSpin,
        &alphaParams,
        &betaParams);

    return XLAL_SUCCESS;
  }

  /**
   * External wrapper for IMRPhenomX_PNR_GenerateEffectiveRingdownFreq
   *
   * Computes the effective ringdown frequency for a given (l,m) multipole
   * following the prescription given in arXiv:##### FIXME: add citation
   *
   */
  REAL8 XLALSimIMRPhenomX_PNR_GenerateEffectiveRingdownFreq(
      REAL8 m1_SI,       /**< mass of companion 1 (kg) */
      REAL8 m2_SI,       /**< mass of companion 2 (kg) */
      REAL8 chi1x,       /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
      REAL8 chi1y,       /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
      REAL8 chi1z,       /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
      REAL8 chi2x,       /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
      REAL8 chi2y,       /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
      REAL8 chi2z,       /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
      REAL8 f_min,       /**< Starting GW frequency (Hz) */
      REAL8 f_max,       /**< Ending GW frequency (Hz) */
      REAL8 fRef,        /**< Reference frequency (Hz) */
      UINT4 ell,         /**< Orbital index */
      UINT4 emmprime,    /**< azimuthal index */
      LALDict *lalParams /**< LAL Dictionary struct */
  )
  {
    /* Ensure we have a dictionary */
    LALDict *lalParams_aux;
    if (lalParams == NULL)
    {
      lalParams_aux = XLALCreateDict();
    }
    else
    {
      lalParams_aux = XLALDictDuplicate(lalParams);
    }

    UINT4 status = 0;
    status = XLALIMRPhenomXPCheckMassesAndSpins(&m1_SI, &m2_SI, &chi1x, &chi1y, &chi1z, &chi2x, &chi2y, &chi2z);
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: XLALIMRPhenomXPCheckMassesAndSpins failed.\n");

    IMRPhenomXWaveformStruct *pWF;
    pWF = XLALMalloc(sizeof(IMRPhenomXWaveformStruct));
    status = IMRPhenomXSetWaveformVariables(pWF, m1_SI, m2_SI, chi1z, chi2z, 0.0, fRef, 0.0, f_min, f_max, 1.0, 0.0, lalParams_aux, 0);
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetWaveformVariables failed.\n");

    /* Initialize IMR PhenomX Precession struct and check that it generated successfully */
    IMRPhenomXPrecessionStruct *pPrec;
    pPrec = XLALMalloc(sizeof(IMRPhenomXPrecessionStruct));
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
        DEBUG);
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetPrecessionVariables failed.\n");

    /* get effective ringdown frequency */
    REAL8 effRD = IMRPhenomX_PNR_GenerateEffectiveRingdownFreq(pWF, ell, emmprime, lalParams_aux);

    /* clean up memory allocation */
    LALFree(pWF);
    LALFree(pPrec);
    XLALDestroyDict(lalParams_aux);

    return effRD;
  }

  /**
   * External wrapper for IMRPhenomX_PNR_GenerateRingdownPNRBeta.
   *
   * Generate PNR beta as described in Eq. 60 or arXiv:2107.08876 and
   * evaluate at the upper connection frequency f_f described in Sec. 8C.
   */
  REAL8 XLALSimIMRPhenomX_PNR_GenerateRingdownPNRBeta(
      REAL8 m1_SI,       /**< mass of companion 1 (kg) */
      REAL8 m2_SI,       /**< mass of companion 2 (kg) */
      REAL8 chi1x,       /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
      REAL8 chi1y,       /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
      REAL8 chi1z,       /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
      REAL8 chi2x,       /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
      REAL8 chi2y,       /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
      REAL8 chi2z,       /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
      REAL8 deltaF,      /**< Sampling Frequency (Hz) */
      REAL8 f_min,       /**< Starting GW frequency (Hz) */
      REAL8 f_max,       /**< Ending GW frequency (Hz) */
      REAL8 fRef_In,     /**< Reference frequency (Hz) */
      LALDict *lalParams /**< LAL Dictionary struct */
  )
  {
    /* Ensure we have a dictionary */
    LALDict *lalParams_aux;
    if (lalParams == NULL)
    {
      lalParams_aux = XLALCreateDict();
    }
    else
    {
      lalParams_aux = XLALDictDuplicate(lalParams);
    }

    int UsePNR = XLALSimInspiralWaveformParamsLookupPhenomXPNRUseTunedAngles(lalParams_aux);
    if (!UsePNR)
    {
      UsePNR = 1;
      XLALSimInspiralWaveformParamsInsertPhenomXPNRUseTunedAngles(lalParams_aux, UsePNR);
    }

    REAL8 fRef = (fRef_In == 0.0) ? f_min : fRef_In;

    UINT4 status;

    REAL8 distance = 1.0;
    REAL8 inclination = 0.0;
    REAL8 phiRef = 0.0;

    /* Initialize IMR PhenomX Waveform struct and check that it initialized correctly */
    IMRPhenomXWaveformStruct *pWF;
    pWF = XLALMalloc(sizeof(IMRPhenomXWaveformStruct));
    status = IMRPhenomXSetWaveformVariables(pWF, m1_SI, m2_SI, chi1z, chi2z, deltaF, fRef, phiRef, f_min, f_max, distance, inclination, lalParams_aux, DEBUG);
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetWaveformVariables failed.\n");

    /* Initialize IMR PhenomX Precession struct and check that it generated successfully */
    IMRPhenomXPrecessionStruct *pPrec;
    pPrec = XLALMalloc(sizeof(IMRPhenomXPrecessionStruct));
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
        DEBUG);
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetPrecessionVariables failed.\n");

    /* get the RD beta */
    REAL8 finalBeta = IMRPhenomX_PNR_GenerateRingdownPNRBeta(pWF, pPrec);

    /* free up memory allocation */
    LALFree(pPrec);
    LALFree(pWF);
    XLALDestroyDict(lalParams_aux);

    return finalBeta;
  }

  /* Function for testing deltaF spacing for interpolation */
  int XLALSimIMRPhenomX_PNR_InterpolationDeltaF(
      INT4 *twospin, /**< UNDOCUMENTED */
      INT4 *precversion, /**< UNDOCUMENTED */
      REAL8 *deltaF, /**< UNDOCUMENTED */
      REAL8 m1_SI, /**< mass of companion 1 (kg) */
      REAL8 m2_SI, /**< mass of companion 2 (kg) */
      REAL8 chi1x, /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
      REAL8 chi1y, /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
      REAL8 chi1z, /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
      REAL8 chi2x, /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
      REAL8 chi2y, /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
      REAL8 chi2z, /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
      REAL8 f_min, /**< Starting GW frequency (Hz) */
      REAL8 f_max, /**< Ending GW frequency (Hz); Defaults to Mf = 0.3 if no f_max is specified. */
      REAL8 fRef,  /**< Reference frequency (Hz) */
      LALDict *lalParams /**< LAL Dictionary struct */)
  {
    /* Ensure we have a dictionary */
    LALDict *lalParams_aux;
    if (lalParams == NULL)
    {
      lalParams_aux = XLALCreateDict();
    }
    else
    {
      lalParams_aux = XLALDictDuplicate(lalParams);
    }

    UINT4 status = 0;
    status = XLALIMRPhenomXPCheckMassesAndSpins(&m1_SI, &m2_SI, &chi1x, &chi1y, &chi1z, &chi2x, &chi2y, &chi2z);
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: XLALIMRPhenomXPCheckMassesAndSpins failed.\n");

    IMRPhenomXWaveformStruct *pWF;
    pWF = XLALMalloc(sizeof(IMRPhenomXWaveformStruct));
    status = IMRPhenomXSetWaveformVariables(pWF, m1_SI, m2_SI, chi1z, chi2z, 0.0, fRef, 0.0, f_min, f_max, 1.0, 0.0, lalParams_aux, 0);
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetWaveformVariables failed.\n");

    /* Initialize IMR PhenomX Precession struct and check that it generated successfully */
    IMRPhenomXPrecessionStruct *pPrec;
    pPrec = XLALMalloc(sizeof(IMRPhenomXPrecessionStruct));

    status = IMRPhenomXGetAndSetPrecessionVariables( // needed to adjust ringdown frequency in pWF
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
        DEBUG);
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetPrecessionVariables failed.\n");

    *deltaF = IMRPhenomX_PNR_HMInterpolationDeltaF(f_min, pWF, pPrec);
    *twospin = IMRPhenomX_PNR_CheckTwoSpin(pPrec);
    *precversion = pPrec->IMRPhenomXPrecVersion;

    LALFree(pWF);
    LALFree(pPrec);
    XLALDestroyDict(lalParams_aux);

    return XLAL_SUCCESS;
  }

  /** @} */
  /** @} */

#ifdef __cplusplus
}
#endif
