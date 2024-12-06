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

  /**
   * \author Eleanor Hamilton, Sebastian Khan, Jonathan E. Thompson
   *
   */

#include "LALSimIMRPhenomXUtilities.h"
#include "LALSimIMRPhenomUtils.h"
#include "LALSimIMRPhenomX_PNR_coefficients.h"
#include "LALSimIMRPhenomX_PNR_internals.h"
#include "LALSimIMRPhenomX_PNR_coefficients.c"
#include "LALSimIMRPhenomX_PNR_alpha.c"
#include "LALSimIMRPhenomX_PNR_beta.c"

#ifndef _OPENMP
#define omp ignore
#endif

#ifndef PHENOMXHMDEBUG
#define DEBUG 0
#else
#define DEBUG 1
#endif

  /**
   * This function computes the required single-spin
   * quantities used to parameterize the MR tuned functions from arXiv:2107.08876.
   *
   * We place these quantities in the already-existing precession struct
   * to avoid extensive code modifications.
   */
  INT4 IMRPhenomX_PNR_GetAndSetPNRVariables(
      IMRPhenomXWaveformStruct *pWF,    /**< PhenomX waveform struct */
      IMRPhenomXPrecessionStruct *pPrec /**< PhenomX precession struct */
  )
  {
    /* check for uninitialized structs */
    XLAL_CHECK(pWF != NULL, XLAL_EFAULT);
    XLAL_CHECK(pPrec != NULL, XLAL_EFAULT);

    /* get needed quantities */
    REAL8 m1 = pWF->m1 * pWF->Mtot;
    REAL8 m2 = pWF->m2 * pWF->Mtot;
    REAL8 q = pWF->q;
    REAL8 chi1x = pPrec->chi1x;
    REAL8 chi1y = pPrec->chi1y;
    REAL8 chi2x = pPrec->chi2x;
    REAL8 chi2y = pPrec->chi2y;

    REAL8 chieff = pWF->chiEff;
    REAL8 chipar = pWF->Mtot * chieff / m1;
    REAL8 chiperp = 0.0;
    REAL8 costheta = 0.0;
    REAL8 chiperp_antisymmetric = 0.0;
    REAL8 theta_antisymmetric = 0.0;

    /* compute effective in-plane spin contribution from Eq. 17 of arXiv:2107.08876 */
    if (q <= 1.5)
    {
      REAL8 chis = sqrt((m1 * m1 * chi1x + m2 * m2 * chi2x) * (m1 * m1 * chi1x + m2 * m2 * chi2x) + (m1 * m1 * chi1y + m2 * m2 * chi2y) * (m1 * m1 * chi1y + m2 * m2 * chi2y)) / (m1 * m1);
      chiperp = sin((q - 1.0) * LAL_PI) * sin((q - 1.0) * LAL_PI) * pPrec->chi_p + cos((q - 1.0) * LAL_PI) * cos((q - 1.0) * LAL_PI) * chis;

      REAL8 antisymmetric_chis = sqrt((m1 * m1 * chi1x - m2 * m2 * chi2x) * (m1 * m1 * chi1x - m2 * m2 * chi2x) + (m1 * m1 * chi1y - m2 * m2 * chi2y) * (m1 * m1 * chi1y - m2 * m2 * chi2y)) / (m1 * m1);
      chiperp_antisymmetric = sin((q - 1.0) * LAL_PI) * sin((q - 1.0) * LAL_PI) * pPrec->chi_p + cos((q - 1.0) * LAL_PI) * cos((q - 1.0) * LAL_PI) * antisymmetric_chis;
    }
    else
    {
      chiperp = pPrec->chi_p;
      chiperp_antisymmetric = pPrec->chi_p;
    }

    /* get the total magnitude, Eq. 18 of arXiv:2107.08876 */
    REAL8 chi_mag = sqrt(chipar * chipar + chiperp * chiperp);
    pPrec->chi_singleSpin = chi_mag;

    REAL8 chi_mag_antisymmetric = sqrt(chipar * chipar + chiperp_antisymmetric * chiperp_antisymmetric);
    pPrec->chi_singleSpin_antisymmetric = chi_mag_antisymmetric;

    /* get the opening angle of the single spin, Eq. 19 of arXiv:2107.08876 */
    if (chi_mag >= 1.0e-6)
    {
      costheta = chipar / chi_mag;
    }
    else
    {
      costheta = 0.;
    }
    if (chi_mag_antisymmetric >= 1.0e-6)
    {
      theta_antisymmetric = acos(chipar / chi_mag_antisymmetric);
    }
    else
    {
      theta_antisymmetric = 0.0;
    }

    pPrec->costheta_singleSpin = costheta;
    pPrec->theta_antisymmetric = theta_antisymmetric; 

    /* compute an approximate final spin using single-spin mapping FIXME: add documentation */
    REAL8 chi1L = chi_mag * costheta;
    REAL8 chi2L = 0.0;

    REAL8 Xfparr = XLALSimIMRPhenomXFinalSpin2017(pWF->eta, chi1L, chi2L);

    /* rescale Xfperp to use the final total mass of 1 */
    REAL8 qfactor = q / (1.0 + q);
    REAL8 Xfperp = qfactor * qfactor * chi_mag * sqrt(1.0 - costheta * costheta);
    REAL8 xf = sqrt(Xfparr * Xfparr + Xfperp * Xfperp);
    if (xf > 1.0e-6)
    {
      pPrec->costheta_final_singleSpin = Xfparr / xf;
    }
    else
    {
      pPrec->costheta_final_singleSpin = 0.;
    }

    /* Initialize frequency values */
    pPrec->PNR_HM_Mflow = 0.0;
    pPrec->PNR_HM_Mfhigh = 0.0;

    /* set angle window boundaries */
    pPrec->PNR_q_window_lower = 8.5;
    pPrec->PNR_q_window_upper = 12.0;
    pPrec->PNR_chi_window_lower = 0.85;
    pPrec->PNR_chi_window_upper = 1.2;

    /* set inspiral scaling flag for HM frequency map */
    pPrec->PNRInspiralScaling = 0;
    if ((q > pPrec->PNR_q_window_upper) || (pPrec->chi_singleSpin > pPrec->PNR_chi_window_upper))
    {
      pPrec->PNRInspiralScaling = 1;
    }

    #if DEBUG == 1
      printf("\nSetting PNR-related single-spin quantities:\n");
      printf("chi_singleSpin                     : %e\n", pPrec->chi_singleSpin);
      printf("chi_singleSpin_antisymmetric       : %e\n", pPrec->chi_singleSpin_antisymmetric);
      printf("costheta_singleSpin                : %e\n", pPrec->costheta_singleSpin);
      /* printf("theta_singleSpin_antisymmetric  : %e\n", pPrec->theta_singleSpin_antisymmetric); */
      printf("theta_antisymmetric  : %e\n", pPrec->theta_antisymmetric); 
      printf("costheta_final_singleSpin          : %e\n\n", pPrec->costheta_final_singleSpin);
    #endif

    return XLAL_SUCCESS;
  }

  /**
   * These two functions create and populate the required
   * parameter structs for PNR.
   *
   * - First, we check to see if we need to do the two-spin mapping;
   *   if so, the approximate single-spin PhenomX structs are computed.
   *
   * - Next we compute the structs for alpha and beta.
   *
   * - Finally, the second function is a wrapper for cleaning up memory.
   */
  INT4 IMRPhenomX_PNR_PopulateStructs(
      IMRPhenomXWaveformStruct **pWF_SingleSpin,     /**< PhenomX waveform struct with single spin parameters */
      IMRPhenomXPrecessionStruct **pPrec_SingleSpin, /**< PhenomX precession struct with single spin parameters */
      IMRPhenomX_PNR_alpha_parameters **alphaParams, /**< alpha paramter struct */
      IMRPhenomX_PNR_beta_parameters **betaParams,   /**< beta paramter struct */
      IMRPhenomXWaveformStruct *pWF,                 /**< PhenomX waveform struct */
      IMRPhenomXPrecessionStruct *pPrec,             /**< PhenomX precession struct */
      LALDict *lalParams                             /**< LAL dictionary struct */
  )
  {
    /* check we should be here */
    INT4 UsePNR = XLALSimInspiralWaveformParamsLookupPhenomXPNRUseTunedAngles(lalParams);
    XLAL_CHECK(
        UsePNR,
        XLAL_EFUNC,
        "Error: IMRPhenomX_PNR_PopulateStructs called without PNR angle activated!\n");

    INT4 status = 0;

    /* check for two-spin effects, if so generate associated single-spin structs */
    if (IMRPhenomX_PNR_CheckTwoSpin(pPrec))
    {
      REAL8 chi = pPrec->chi_singleSpin;
      REAL8 costheta = pPrec->costheta_singleSpin;
      *pWF_SingleSpin = XLALMalloc(sizeof(IMRPhenomXWaveformStruct));
      status = IMRPhenomXSetWaveformVariables(
          *pWF_SingleSpin,
          pWF->m1_SI,
          pWF->m2_SI,
          pWF->chi1L,
          pWF->chi2L,
          pWF->deltaF,
          pWF->fRef,
          pWF->phi0,
          pWF->fMin,
          pWF->fMax,
          pWF->distance,
          pWF->inclination,
          lalParams,
          DEBUG);
      XLAL_CHECK(
          XLAL_SUCCESS == status,
          XLAL_EFUNC,
          "Error: IMRPhenomXSetWaveformVariables failed.\n");

      /* the only needed changes for the single-spin mapping in pWF
       * are to set the condition PNR_SINGLE_SPIN to overcome the
       * Kerr limit checks for pPrec, and to update chiEff */
      (*pWF_SingleSpin)->PNR_SINGLE_SPIN = 1;
      (*pWF_SingleSpin)->chiEff = XLALSimIMRPhenomXchiEff(pWF->eta, chi * costheta, 0.0);

      *pPrec_SingleSpin = XLALMalloc(sizeof(IMRPhenomXPrecessionStruct));
      status = IMRPhenomXGetAndSetPrecessionVariables(
          *pWF_SingleSpin,
          *pPrec_SingleSpin,
          pWF->m1_SI,
          pWF->m2_SI,
          chi * sin(acos(costheta)),
          0.0,
          chi * costheta,
          0.0,
          0.0,
          0.0,
          lalParams,
          DEBUG);
      XLAL_CHECK(
          XLAL_SUCCESS == status,
          XLAL_EFUNC,
          "Error: IMRPhenomXGetAndSetPrecessionVariables failed.\n");
    }

    /* generate alpha parameters */
    *alphaParams = XLALMalloc(sizeof(IMRPhenomX_PNR_alpha_parameters));
    status = IMRPhenomX_PNR_precompute_alpha_coefficients(*alphaParams, pWF, pPrec);
    XLAL_CHECK(
        XLAL_SUCCESS == status,
        XLAL_EFUNC,
        "Error: IMRPhenomX_PNR_precompute_alpha_coefficients failed.\n");

    status = IMRPhenomX_PNR_alpha_connection_parameters(*alphaParams, pWF, pPrec);
    XLAL_CHECK(
        XLAL_SUCCESS == status,
        XLAL_EFUNC,
        "Error: IMRPhenomX_PNR_alpha_connection_parameters failed.\n");

    /* generate beta parameters */
    *betaParams = XLALMalloc(sizeof(IMRPhenomX_PNR_beta_parameters));
    status = IMRPhenomX_PNR_precompute_beta_coefficients(*betaParams, pWF, pPrec);
    XLAL_CHECK(
        XLAL_SUCCESS == status,
        XLAL_EFUNC,
        "Error: IMRPhenomX_PNR_precompute_beta_coefficients failed.\n");

    status = IMRPhenomX_PNR_beta_connection_parameters(*betaParams, pWF, pPrec, *pWF_SingleSpin, *pPrec_SingleSpin);
    XLAL_CHECK(
        XLAL_SUCCESS == status,
        XLAL_EFUNC,
        "Error: IMRPhenomX_PNR_beta_connection_parameters failed.\n");

#if DEBUG == 1
    IMRPhenomX_PNR_AngleParameterDebugPrint(*alphaParams, *betaParams);
#endif

    return XLAL_SUCCESS;
  }

  void IMRPhenomX_PNR_FreeStructs(
      IMRPhenomXWaveformStruct **pWF_SingleSpin,
      IMRPhenomXPrecessionStruct **pPrec_SingleSpin,
      IMRPhenomX_PNR_alpha_parameters **alphaParams,
      IMRPhenomX_PNR_beta_parameters **betaParams)
  {
    if (*pWF_SingleSpin != NULL)
    {
      LALFree(*pWF_SingleSpin);
    }
    if (*pPrec_SingleSpin != NULL)
    {
      if((*pPrec_SingleSpin)->pWF22AS)
      {
        LALFree((*pPrec_SingleSpin)->pWF22AS);
      }
      LALFree(*pPrec_SingleSpin);
    }
    if (*alphaParams != NULL)
    {
      LALFree(*alphaParams);
    }
    if (*betaParams != NULL)
    {
      LALFree(*betaParams);
    }
  }

  /**
   * This function computes the frequency integral
   * outlined in Eq. A7 of arXiv:2107.08876.
   *
   * Given frequency arrays of alpha and beta, compute
   * spline interpolants and then use Boole's rule to
   * numerically integrate.
   */
  INT4 IMRPhenomX_PNR_GeneratePNRGamma(
      REAL8Sequence *gamma,       /**< gamma frequency series (rad) */
      const REAL8Sequence *freqs, /**< input frequencies (Hz) */
      const REAL8Sequence *alpha, /**< alpha frequency series (rad) */
      const REAL8Sequence *beta   /**< beta frequency series (rad) */
  )
  {
    /* Check that all is well with the inputs */
    XLAL_CHECK(gamma != NULL, XLAL_EFAULT);
    XLAL_CHECK(freqs != NULL, XLAL_EFAULT);
    XLAL_CHECK(alpha != NULL, XLAL_EFAULT);
    XLAL_CHECK(beta != NULL, XLAL_EFAULT);

    /* Initialise the struct containing splines of alpha and beta */
    IMRPhenomX_PNR_angle_spline ab_splines;

    /* Construct splines of alpha and beta*/
    gsl_spline *alpha_spline = gsl_spline_alloc(gsl_interp_cspline, freqs->length);
    gsl_spline *beta_spline = gsl_spline_alloc(gsl_interp_cspline, freqs->length);

    gsl_interp_accel *alpha_acc = gsl_interp_accel_alloc();
    gsl_interp_accel *beta_acc = gsl_interp_accel_alloc();

    gsl_spline_init(alpha_spline, freqs->data, alpha->data, freqs->length);
    gsl_spline_init(beta_spline, freqs->data, beta->data, freqs->length);

    /* Populate combined angle spline struct */
    ab_splines.alpha_spline = alpha_spline;
    ab_splines.alpha_acc = alpha_acc;
    ab_splines.beta_spline = beta_spline;
    ab_splines.beta_acc = beta_acc;

    /* Initialize gamma to zero. Will be corrected with the reference value later. */
    gamma->data[0] = 0.0;

    /* Numerically integrate \gamma' = -\alpha' \cos(\beta) to get gamma
     * Eq. A7 of arXiv:2107.08876
     * Integration is done using Boole's Rule. */
    for (size_t i = 1; i < freqs->length; i++)
    {
      REAL8 t1 = freqs->data[i - 1];
      REAL8 t2 = freqs->data[i];
      /* Boole's Rule with stepsize h = deltaF / 4 */
      gamma->data[i] = (gamma->data[i - 1] + (1.0 / 90.0) * (t2 - t1) * (7.0 * IMRPhenomX_PNR_alphadot_cosbeta(t1, &ab_splines) + 32.0 * IMRPhenomX_PNR_alphadot_cosbeta((t1 + 3.0 * t2) / 4.0, &ab_splines) + 12.0 * IMRPhenomX_PNR_alphadot_cosbeta((t1 + t2) / 2.0, &ab_splines) + 32.0 * IMRPhenomX_PNR_alphadot_cosbeta((3.0 * t1 + t2) / 4.0, &ab_splines) + 7.0 * IMRPhenomX_PNR_alphadot_cosbeta(t2, &ab_splines)));
    }

    /* Free up memory allocation */
    gsl_spline_free(alpha_spline);
    gsl_spline_free(beta_spline);

    gsl_interp_accel_free(alpha_acc);
    gsl_interp_accel_free(beta_acc);

    return XLAL_SUCCESS;
  }

  /**
   * This function computes the frequency integral
   * outlined in Eq. A7 of arXiv:2107.08876.
   *
   * Given spline interpolants of alpha and beta, compute
   * use Boole's rule to numerically integrate.
   */
  int IMRPhenomX_PNR_GeneratePNRGamma_FromInterpolants(
      REAL8Sequence *gamma,                  /**< frequency series for gamma (rad) */
      const REAL8Sequence *freqs,            /**< input frequencies (Hz) */
      IMRPhenomX_PNR_angle_spline *ab_splines /**< pnr angle spline struct */
  )
  {
    /* Check that everything is initialized */
    XLAL_CHECK(gamma != NULL, XLAL_EFAULT);
    XLAL_CHECK(freqs != NULL, XLAL_EFAULT);
    XLAL_CHECK(ab_splines != NULL, XLAL_EFAULT);

    /* Initialize gamma to zero. Will be corrected with the reference value later. */
    gamma->data[0] = 0.0;

    /* Numerically integrate \gamma' = -\alpha' \cos(\beta) to get gamma
     * Eq. A7 of arXiv:2107.08876
     * Integration is done using Boole's Rule. */
    for (size_t i = 1; i < freqs->length; i++)
    {
      REAL8 t1 = freqs->data[i - 1];
      REAL8 t2 = freqs->data[i];
      /* Boole's Rule with stepsize h = deltaF / 4 */
      gamma->data[i] = (gamma->data[i - 1] + (1.0 / 90.0) * (t2 - t1) * (7.0 * IMRPhenomX_PNR_alphadot_cosbeta(t1, ab_splines) + 32.0 * IMRPhenomX_PNR_alphadot_cosbeta((t1 + 3.0 * t2) / 4.0, ab_splines) + 12.0 * IMRPhenomX_PNR_alphadot_cosbeta((t1 + t2) / 2.0, ab_splines) + 32.0 * IMRPhenomX_PNR_alphadot_cosbeta((3.0 * t1 + t2) / 4.0, ab_splines) + 7.0 * IMRPhenomX_PNR_alphadot_cosbeta(t2, ab_splines)));
    }

    return XLAL_SUCCESS;
  }

  /**
   * Wrapper function for computing the integrand in
   * Eq. A7 of arXiv:2107.08876
   */
  REAL8 IMRPhenomX_PNR_alphadot_cosbeta(
      REAL8 f,                            /**< frequency (Hz) */
      IMRPhenomX_PNR_angle_spline *params /**< pnr angle interpolant struct */
  )
  {
    REAL8 alphadot = gsl_spline_eval_deriv(params->alpha_spline, f, params->alpha_acc);
    REAL8 beta = gsl_spline_eval(params->beta_spline, f, params->beta_acc);

    return -1.0 * alphadot * cos(beta);
  }

  /**
   * Computes a linear frequency map for the HM PNR angles
   * based on the linear frequency mapping used originally in
   * PhenomHM (Eq. 5 of arXiv:1708.00404)
   *
   * The specifics are given in Eq. #### of FIXME: add documentation
   */
  REAL8 IMRPhenomX_PNR_LinearFrequencyMap(
      REAL8 Mf,       /**< geometric evaluation frequency */
      REAL8 ell,      /**< polar index */
      REAL8 emm,      /**< azimuthal index */
      REAL8 Mf_lower, /**< lower geometric transition frequency */
      REAL8 Mf_upper, /**< upper geometric transition frequency */
      REAL8 Mf_RD_22, /**< (2,2) geometric ringdown frequency */
      REAL8 Mf_RD_lm,  /**< (l,m) geometric ringdown frequency */
      UINT4 INSPIRAL  /**< flag to toggle inspiral scaling only */
  )
  {
    /* initial break if mapping (2,2) */
    if ((ell == 2) && (emm == 2))
    {
      return Mf;
    }

    /* below the lower transition, or outside calibration region, using SPA scaling */
    if ((Mf <= Mf_lower) || INSPIRAL)
    {
      return 2.0 * Mf / emm;
    }

    /* check that the lower and upper transition frequencies make sense */
    XLAL_CHECK(
        Mf_lower <= Mf_upper,
        XLAL_EFUNC,
        "Error: Low-frequency cutoff is greater than high-frequency cutoff in IMRPhenomX_PNR_LinearFrequencyMap.\n");

    /* above the upper transition, shift by the difference of the RD frequencies */
    if (Mf > Mf_upper)
    {
      return Mf - (Mf_RD_lm - Mf_RD_22);
    }

    /* construct a linear transition between Mf_lower and Mf_upper scalings */
    REAL8 A = IMRPhenomX_PNR_LinearFrequencySlope(emm, Mf_lower, Mf_upper, Mf_RD_22, Mf_RD_lm);
    REAL8 B = 2.0 * Mf_lower / emm;

    return A * (Mf - Mf_lower) + B;
  }

  /**
   * Computes the slope of the linear frequency interpolation
   * mapping used for the HM PNR angles.
   *
   * See Eq. ### of FIXME: add citation
   */
  REAL8 IMRPhenomX_PNR_LinearFrequencySlope(
      REAL8 emm,      /**< azimuthal index */
      REAL8 Mf_lower, /**< lower geometric transition frequency */
      REAL8 Mf_upper, /**< upper geometric transition frequency */
      REAL8 Mf_RD_22, /**< (2,2) geometric ringdown frequency */
      REAL8 Mf_RD_lm  /**< (l,m) geometric ringdown frequency */
  )
  {
    /* check that the lower and upper transition frequencies make sense */
    XLAL_CHECK(
        Mf_lower <= Mf_upper,
        XLAL_EFUNC,
        "Error: Low-frequency cutoff is greater than high-frequency cutoff in IMRPhenomX_PNR_LinearFrequencySlope.\n");

    /* compute the slope */
    return (Mf_upper - (Mf_RD_lm - Mf_RD_22) - 2.0 * Mf_lower / emm) / (Mf_upper - Mf_lower);
  }

  /**
   * Compute the transition frequencies for the HM PNR angle
   * mapping.
   *
   * See Eq. ### of FIXME: add documentation
   */
  INT4 IMRPhenomX_PNR_LinearFrequencyMapTransitionFrequencies(
      REAL8 *Mf_low,                    /**< lower transition frequency */
      REAL8 *Mf_high,                   /**< upper transition frequency */
      REAL8 emmprime,                   /**< azimuthal index */
      REAL8 Mf_RD_22,                   /**< (2,2) geometric ringdown frequency */
      REAL8 Mf_RD_lm,                   /**< (l,m) geometric ringdown frequency */
      IMRPhenomXPrecessionStruct *pPrec /**< PhenomX precession struct */
  )
  {
    /* Check that everything is initialized */
    XLAL_CHECK(Mf_low != NULL, XLAL_EFAULT);
    XLAL_CHECK(Mf_high != NULL, XLAL_EFAULT);
    XLAL_CHECK(pPrec != NULL, XLAL_EFAULT);

    /* specify choices for connection frequencies */
    *Mf_low = 0.65 * pPrec->PNR_HM_Mflow * emmprime / 2.0;
    *Mf_high = 1.1 * (pPrec->PNR_HM_Mfhigh + (Mf_RD_lm - Mf_RD_22));

    /* check for cases with (2,1) multipole where Mf_high might misbehave */
    if ((*Mf_high < 0.0)||((Mf_RD_lm - Mf_RD_22 < 0.0)&&(*Mf_high < pPrec->PNR_HM_Mfhigh / 2.0))){
      *Mf_high = pPrec->PNR_HM_Mfhigh;
    }

    return XLAL_SUCCESS;
  }

  /**
   * Here we compute an appropriate deltaF to be used when generating
   * the (2,2) angle interpolants and mapping them to the HMs.
   *
   * FIXME: add documentation
   */
  REAL8 IMRPhenomX_PNR_HMInterpolationDeltaF(
      REAL8 f_min,                      /**< minimum starting frequency (Hz) */
      IMRPhenomXWaveformStruct *pWF,    /**< PhenomX waveform struct */
      IMRPhenomXPrecessionStruct *pPrec /**< PhenomX precession struct */
  )
  {
    /* First compute the deltaF required for a single-spin system. This assumes error
     * scaling roughly like deltaF**4 in the cubic spline. */

    REAL8 error_tolerance = pPrec->IMRPhenomXPNRInterpTolerance;

    /* aligned-spin limit means we don't need high resolution on the prec angles (which should be zero) */
    if((pPrec->chi1_perp==0.0)&&(pPrec->chi2_perp==0.0))
    {
      if(pWF->deltaF != 0){
        return pWF->deltaF;
      }

      return 0.1;
    }

    /* Compute approximate scaling from leading-order alpha contribution, a la
     * - Eq. 5.6 of arXiv:2004.06503
     * - Eq. 2.8 of arXiv:2001.10897 */
    REAL8 Mf = XLALSimPhenomUtilsHztoMf(f_min, pWF->Mtot);
    REAL8 eta_term = sqrt(1.0 - 4.0 * pWF->eta);
    REAL8 numerator = 3.0 * LAL_PI * Mf * Mf * Mf * Mf * Mf * error_tolerance * (1.0 + eta_term);
    REAL8 denominator = 7.0 + 13.0 * eta_term;
    REAL8 constant = 4.0 * sqrt(2.0 / 5.0);

    REAL8 deltaF_alpha = XLALSimPhenomUtilsMftoHz(constant * sqrt(sqrt(numerator / denominator)), pWF->Mtot);

    /* Check for two-spin oscillations!
     * We approximate the frequency sampling required
     * by first approximating the frequency of oscillation
     * in the magnitude of the total spin computed by the MSA expansion */
    if (IMRPhenomX_PNR_CheckTwoSpin(pPrec))
    {
      /* we have a well-behaved MSA expansion, so we grab precomputed terms */
      REAL8 g0 = pPrec->g0;
      REAL8 deltam = pPrec->delta_qq;
      REAL8 psi1 = pPrec->psi1;
      REAL8 psi2 = pPrec->psi2;

      REAL8 v0 = cbrt(LAL_PI * Mf);
      REAL8 v02 = v0 * v0;
      REAL8 v06 = v02 * v02 * v02;

      /* this is the frequency derivative of psi, Eq. 51 of arXiv:1703.03967
       * approximate frequency of the two-spin oscillations */
      REAL8 dpsi = g0 * deltam * LAL_PI / (4.0 * v06) * (3.0 + 2.0 * psi1 * v0 + psi2 * v02);
      /* we can approximate the period by the inverse of dpsi */
      REAL8 dpsiInv = fabs(1.0 / dpsi);

      /* when oscillations in beta bring its value close to zero, alpha can undergo
       * sudden shifts by factors of pi. We need to predict this to provide enough sampling
       * for good interpolation.
       *
       * To do this, we first make sure that we have a large-enough mass-ratio to provide noticeable drift
       * in the in-plane spin orientations.
       *
       * Next we compute the minimum and maximum values that beta can take at the starting frequency
       * of the interpolation. If betaMin is very close to zero AND it is significantly smaller than betaMax,
       * then we have large oscillations and the potential for sharp jumps in alpha. We then divide the
       * approximate two-spin oscillation period by 4 to increase sampling resolution. */

      REAL8 L_fmin = pWF->Mtot * pWF->Mtot * XLALSimIMRPhenomXLPNAnsatz(v0, pWF->eta / v0, pPrec->L0, pPrec->L1, pPrec->L2, pPrec->L3, pPrec->L4, pPrec->L5, pPrec->L6, pPrec->L7, pPrec->L8, pPrec->L8L);

      REAL8 betaMin = atan2(fabs(pPrec->S1_perp - pPrec->S2_perp), L_fmin + pPrec->SL);
      REAL8 betaMax = atan2(pPrec->S1_perp + pPrec->S2_perp, L_fmin + pPrec->SL);

#if DEBUG == 1
        // Save angles into a file
        FILE *fileangle;
        char fileSpecIII[40];
        sprintf(fileSpecIII, "PNR_HM_interpolation_data.dat");
        fileangle = fopen(fileSpecIII, "w");

        fprintf(fileangle, "dpsi = %.16e \n", dpsi);
        fprintf(fileangle, "pPrec->S1_perp = %.16e \n", pPrec->S1_perp);
        fprintf(fileangle, "pPrec->S2_perp = %.16e \n", pPrec->S2_perp);
        fprintf(fileangle, "pPrec->SL = %.16e \n", pPrec->SL);
        fprintf(fileangle, "pPrec->LRef = %.16e \n\n", pPrec->LRef);
        fprintf(fileangle, "L_fmin = %.16e \n\n", L_fmin);
        fprintf(fileangle, "betaMin = %.16e \n", betaMin);
        fprintf(fileangle, "betaMax = %.16e \n", betaMax);

        fclose(fileangle);
#endif

      if ((betaMin < 0.01) && (betaMin / betaMax < 0.55))
      {
        dpsiInv /= 4;
      }


      /* sample 4 points per oscillaion */
      REAL8 deltaF_twospin = XLALSimPhenomUtilsMftoHz(dpsiInv / 4.0, pWF->Mtot);

      /* check that deltaF_twospin is smaller than the spacing required by alpha */
      if ((deltaF_twospin < deltaF_alpha) && !isnan(dpsi))
      {
        /* put a hard limit of 0.01 on deltaF to avoid memory saturation */
        return (deltaF_twospin < 1e-2) ? 1e-2 : deltaF_twospin;
      }
    }
    /* put a hard limit of 0.01 on deltaF to avoid memory saturation */
    return (deltaF_alpha < 1e-2) ? 1e-2 : deltaF_alpha;
  }

  /**
   * This function quickly checks to see if we expect two-spin effects
   * to be present in the inspiral precession angles.
   */
  INT4 IMRPhenomX_PNR_CheckTwoSpin(
      IMRPhenomXPrecessionStruct *pPrec /**< PhenomX precession struct */
  )
  {
    /* Ensure we have:
     * - non-zero spin on the primary
     * - non-trivial spin on the secondary
     * - activated the MSA angles
     */
    if ((pPrec->chi1_norm != 0.0) && (pPrec->chi2_norm >= 1.0e-3) && (pPrec->IMRPhenomXPrecVersion >= 200))
    {
      return 1;
    }

    return 0;
  }

  /**
   * Evaluates a function at two points and interpolates between them.
   */
  REAL8 IMRPhenomX_PNR_AngleAtFRef(
      const REAL8Sequence *angle, /**< input angle array (rad) */
      const REAL8 f_ref,          /**< reference frequency value (Hz) */
      const REAL8Sequence *freqs, /**< evaluation frequency */
      const REAL8 deltaF          /**< evaluation frequency */
  )
  {
    /* this function assumes a uniform frequency sampling */
    /* Check that everything is initialized */
    XLAL_CHECK(angle != NULL, XLAL_EFAULT);
    XLAL_CHECK(freqs != NULL, XLAL_EFAULT);

    /* check that f_ref is within range */
    REAL8 f_min = freqs->data[0];
    REAL8 f_max = freqs->data[freqs->length - 1];

    XLAL_CHECK((f_min <= f_ref) && (f_ref <= f_max),
               XLAL_EFUNC,
               "Error: f_ref does not fall within the evaluated frequencies of the angle in IMRPhenomX_PNR_AngleAtFRef.\n");

    size_t idx_fmin = (size_t)(f_min / deltaF);
    size_t idx_eval = (f_ref == f_min) ? 0 : (size_t)(f_ref / deltaF) - idx_fmin;
    size_t idx_eval_p1 = idx_eval + 1;

    return IMRPhenomX_PNR_LinearInterpolate(
        angle->data[idx_eval], angle->data[idx_eval_p1],
        freqs->data[idx_eval], freqs->data[idx_eval_p1],
        f_ref);
  }

  /**
   * Evaluates a function at two points and interpolates between them.
   */
  REAL8 IMRPhenomX_PNR_LinearInterpolate(
      REAL8 a0,   /**< function evaluated at f0 */
      REAL8 a1,   /**< function evaluated at f1 */
      REAL8 f0,   /**< lower frequency value */
      REAL8 f1,   /**< upper frequency value */
      REAL8 feval /**< evaluation frequency */
  )
  {
    XLAL_CHECK(
        (f0 <= feval) && (feval <= f1),
        XLAL_EFUNC,
        "Error: the evaluation frequency is not between the two sampling frequencies in IMRPhenomX_PNR_LinearInterpolate.\n");

    /* evaluate linear interpolation between f0 and f1 */
    return a0 + (feval - f0) * (a1 - a0) / (f1 - f0);
  }

  /**
   * This code recomputes the skymapped locations in the J-frame
   * using the new value of beta computed from the model. This beta
   * corresponds to the orientation of the maximal emission direction
   * relative to J, as opposed to the orientation of L.
   *
   * This code is mostly copied from LAlSimIMRPhenomX_precession.c with slight modifications.
   */
  INT4 IMRPhenomX_PNR_RemapThetaJSF(
      REAL8 beta_ref,                    /**< reference opening angle (rad) */
      IMRPhenomXWaveformStruct *pWF,     /**< PhenomX waveform struct */
      IMRPhenomXPrecessionStruct *pPrec, /**< PhenomX precession struct */
      LALDict *lalParams                 /**< LAL dictionary struct */
  )
  {
    /* Check that everything is initialized */
    XLAL_CHECK(pWF != NULL, XLAL_EFAULT);
    XLAL_CHECK(pPrec != NULL, XLAL_EFAULT);
    XLAL_CHECK(lalParams != NULL, XLAL_EFAULT);

    /* first check to see if nothing needs doing */
    if ((pPrec->J0 < 1e-10))
    {
      return XLAL_SUCCESS;
    }
    REAL8 chi_in_plane = sqrt(pPrec->chi1x*pPrec->chi1x+pPrec->chi1y*pPrec->chi1y+pPrec->chi2x*pPrec->chi2x+pPrec->chi2y*pPrec->chi2y);
    if (chi_in_plane < 1e-3)
    {
      pPrec->alpha_offset -= pPrec->alpha0;
      return XLAL_SUCCESS;
    }

    /* replace thetaJN with beta at f_Ref */
    pPrec->thetaJ_Sf = beta_ref;

    /* re-run all the following code to get thetaJN, correct polarizations, and angle offsets */
    /* annoying to duplicate code, but it's the easiest way to modify thetaJN */
    const double phiRef = pWF->phiRef_In;

    INT4 convention = XLALSimInspiralWaveformParamsLookupPhenomXPConvention(lalParams);
    if (!(convention == 0 || convention == 1 || convention == 5 || convention == 6 || convention == 7))
    {
      XLAL_ERROR(XLAL_EINVAL, "Error: IMRPhenomXPConvention not recognized. Requires version 0, 1, 5, 6 or 7.\n");
    }

    /*
    Here we follow the same prescription as in IMRPhenomPv2:

    Now rotate from SF to J frame to compute alpha0, the azimuthal angle of LN, as well as
    thetaJ, the angle between J and N.

    The J frame is defined by imposing that J points in the z-direction and the line of sight N is in the xz-plane
    (with positive projection along x).

    The components of any vector in the (new) J-frame can be obtained by rotation from the (old) source frame (SF).
    This is done by multiplying by: RZ[-kappa].RY[-thetaJ].RZ[-phiJ]

    Note that kappa is determined by rotating N with RY[-thetaJ].RZ[-phiJ], which brings J to the z-axis, and
    taking the opposite of the azimuthal angle of the rotated N.
    */

    /* Determine kappa via rotations, as above */
    pPrec->Nx_Sf = sin(pWF->inclination) * cos((LAL_PI / 2.0) - phiRef);
    pPrec->Ny_Sf = sin(pWF->inclination) * sin((LAL_PI / 2.0) - phiRef);
    pPrec->Nz_Sf = cos(pWF->inclination);

    REAL8 tmp_x = pPrec->Nx_Sf;
    REAL8 tmp_y = pPrec->Ny_Sf;
    REAL8 tmp_z = pPrec->Nz_Sf;

    IMRPhenomX_rotate_z(-pPrec->phiJ_Sf, &tmp_x, &tmp_y, &tmp_z);
    IMRPhenomX_rotate_y(-pPrec->thetaJ_Sf, &tmp_x, &tmp_y, &tmp_z);

    /* Note difference in overall - sign w.r.t PhenomPv2 code */
    pPrec->kappa = XLALSimIMRPhenomXatan2tol(tmp_y, tmp_x, MAX_TOL_ATAN);

    /* Now determine alpha0 by rotating LN. In the source frame, LN = {0,0,1} */
    tmp_x = 0.0;
    tmp_y = 0.0;
    tmp_z = 1.0;
    IMRPhenomX_rotate_z(-pPrec->phiJ_Sf, &tmp_x, &tmp_y, &tmp_z);
    IMRPhenomX_rotate_y(-pPrec->thetaJ_Sf, &tmp_x, &tmp_y, &tmp_z);
    IMRPhenomX_rotate_z(-pPrec->kappa, &tmp_x, &tmp_y, &tmp_z);

    if (fabs(tmp_x) < MAX_TOL_ATAN && fabs(tmp_y) < MAX_TOL_ATAN)
    {
/* This is the aligned spin case */
#if DEBUG == 1
      printf("\nAligned-spin case.\n");
#endif

      switch (convention)
      {
      case 0:
      case 5:
      {
        pPrec->alpha0 = LAL_PI;
        break;
      }
      case 1:
      case 6:
      case 7:
      {
        pPrec->alpha0 = LAL_PI - pPrec->kappa;
        break;
      }
      }
    }
    else
    {
      switch (convention)
      {
      case 0:
      case 5:
      {
        pPrec->alpha0 = atan2(tmp_y, tmp_x);
        break;
      }
      case 1:
      case 6:
      case 7:
      {
        pPrec->alpha0 = LAL_PI - pPrec->kappa;
        break;
      }
      }
    }

    /* update alpha_offset with alpha0 */

    pPrec->alpha_offset -= pPrec->alpha0;

    /* remove convention options for thetaJN using PNR, so that PNR beta is used instead of theta_JL */

    /* Now determine thetaJN by rotating N */
    tmp_x = pPrec->Nx_Sf;
    tmp_y = pPrec->Ny_Sf;
    tmp_z = pPrec->Nz_Sf;
    IMRPhenomX_rotate_z(-pPrec->phiJ_Sf, &tmp_x, &tmp_y, &tmp_z);
    IMRPhenomX_rotate_y(-pPrec->thetaJ_Sf, &tmp_x, &tmp_y, &tmp_z);
    IMRPhenomX_rotate_z(-pPrec->kappa, &tmp_x, &tmp_y, &tmp_z);

    /* We don't need the y-component but we will store it anyway */
    pPrec->Nx_Jf = tmp_x;
    pPrec->Ny_Jf = tmp_y;
    pPrec->Nz_Jf = tmp_z;

    /* This is a unit vector, so no normalization */
    pPrec->thetaJN = acos(pPrec->Nz_Jf);

    /*
      Define the polarizations used. This follows the conventions adopted for IMRPhenomPv2.

      The IMRPhenomP polarizations are defined following the conventions in Arun et al (arXiv:0810.5336),
      i.e. projecting the metric onto the P, Q, N triad defining where: P = (N x J) / |N x J|.

      However, the triad X,Y,N used in LAL (the "waveframe") follows the definition in the
      NR Injection Infrastructure (Schmidt et al, arXiv:1703.01076).

      The triads differ from each other by a rotation around N by an angle \zeta. We therefore need to rotate
      the polarizations by an angle 2 \zeta.
    */
    pPrec->Xx_Sf = -cos(pWF->inclination) * sin(phiRef);
    pPrec->Xy_Sf = -cos(pWF->inclination) * cos(phiRef);
    pPrec->Xz_Sf = +sin(pWF->inclination);

    tmp_x = pPrec->Xx_Sf;
    tmp_y = pPrec->Xy_Sf;
    tmp_z = pPrec->Xz_Sf;

    IMRPhenomX_rotate_z(-pPrec->phiJ_Sf, &tmp_x, &tmp_y, &tmp_z);
    IMRPhenomX_rotate_y(-pPrec->thetaJ_Sf, &tmp_x, &tmp_y, &tmp_z);
    IMRPhenomX_rotate_z(-pPrec->kappa, &tmp_x, &tmp_y, &tmp_z);

    /*
      The components tmp_i are now the components of X in the J frame.

      We now need the polar angle of this vector in the P, Q basis of Arun et al:

          P = (N x J) / |NxJ|

      Note, that we put N in the (pos x)z half plane of the J frame
    */

    switch (convention)
    {
    case 0:
    case 5:
    {
      /* Get polar angle of X vector in J frame in the P,Q basis of Arun et al */
      pPrec->PArunx_Jf = +0.0;
      pPrec->PAruny_Jf = -1.0;
      pPrec->PArunz_Jf = +0.0;

      /* Q = (N x P) by construction */
      pPrec->QArunx_Jf = pPrec->Nz_Jf;
      pPrec->QAruny_Jf = 0.0;
      pPrec->QArunz_Jf = -pPrec->Nx_Jf;
      break;
    }
    case 1:
    case 6:
    case 7:
    {
      /* Get polar angle of X vector in J frame in the P,Q basis of Arun et al */
      pPrec->PArunx_Jf = pPrec->Nz_Jf;
      pPrec->PAruny_Jf = 0;
      pPrec->PArunz_Jf = -pPrec->Nx_Jf;

      /* Q = (N x P) by construction */
      pPrec->QArunx_Jf = 0;
      pPrec->QAruny_Jf = 1;
      pPrec->QArunz_Jf = 0;
      break;
    }
    }

    // (X . P)
    pPrec->XdotPArun = (tmp_x * pPrec->PArunx_Jf) + (tmp_y * pPrec->PAruny_Jf) + (tmp_z * pPrec->PArunz_Jf);

    // (X . Q)
    pPrec->XdotQArun = (tmp_x * pPrec->QArunx_Jf) + (tmp_y * pPrec->QAruny_Jf) + (tmp_z * pPrec->QArunz_Jf);

    /* Now get the angle zeta */
    pPrec->zeta_polarization = atan2(pPrec->XdotQArun, pPrec->XdotPArun);

    const REAL8 ytheta = pPrec->thetaJN;
    const REAL8 yphi = 0.0;
    pPrec->Y2m2 = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 2, -2);
    pPrec->Y2m1 = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 2, -1);
    pPrec->Y20 = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 2, 0);
    pPrec->Y21 = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 2, 1);
    pPrec->Y22 = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 2, 2);
    pPrec->Y3m3 = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 3, -3);
    pPrec->Y3m2 = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 3, -2);
    pPrec->Y3m1 = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 3, -1);
    pPrec->Y30 = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 3, 0);
    pPrec->Y31 = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 3, 1);
    pPrec->Y32 = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 3, 2);
    pPrec->Y33 = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 3, 3);
    pPrec->Y4m4 = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 4, -4);
    pPrec->Y4m3 = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 4, -3);
    pPrec->Y4m2 = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 4, -2);
    pPrec->Y4m1 = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 4, -1);
    pPrec->Y40 = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 4, 0);
    pPrec->Y41 = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 4, 1);
    pPrec->Y42 = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 4, 2);
    pPrec->Y43 = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 4, 3);
    pPrec->Y44 = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 4, 4);

    return XLAL_SUCCESS;
  }

  /**
   * This code recomputes the skymapped locations in the J-frame
   * using the new value of beta computed from the model. This beta
   * corresponds to the orientation of the maximal emission direction
   * relative to J, as opposed to the orientation of L.
   *
   * This code is mostly copied from LAlSimIMRPhenomX_precession.c with slight modifications.
   */
  REAL8 IMRPhenomX_PNR_GenerateEffectiveRingdownFreq(
      IMRPhenomXWaveformStruct *pWF, /**< PhenomX waveform struct */
      UINT4 ell,                     /**< polar index */
      UINT4 emmprime,                /**< azimuthal index */
      LALDict *lalParams             /**< LAL Dictionary struct */
  )
  {
    /* Check that everything is initialized */
    XLAL_CHECK(pWF != NULL, XLAL_EFAULT);
    XLAL_CHECK(lalParams != NULL, XLAL_EFAULT);

    REAL8 effRD = 0.0;

    /* if (2,2) return normal RD frequency */
    if ((ell == 2) & (emmprime == 2))
    {
      effRD = pWF->fRING;
    }
    else
    {
      /* compute QNM and pWFHM structs */
      QNMFits *qnms = (QNMFits *)XLALMalloc(sizeof(QNMFits));
      IMRPhenomXHM_Initialize_QNMs(qnms);
      IMRPhenomXHMWaveformStruct *pWFHM = (IMRPhenomXHMWaveformStruct *)XLALMalloc(sizeof(IMRPhenomXHMWaveformStruct));
      IMRPhenomXHM_SetHMWaveformVariables(ell, emmprime, pWFHM, pWF, qnms, lalParams);

      /* grab effective RD frequency */
      effRD = pWFHM->fRING;

      /* clean up memory allocation */
      LALFree(pWFHM);
      LALFree(qnms);
    }

    return effRD;
  }

  /**
   * Print various parameters in the angle structs.
   */
  void IMRPhenomX_PNR_AngleParameterDebugPrint(
      IMRPhenomX_PNR_alpha_parameters *alphaParams, /**< alpha parameter struct */
      IMRPhenomX_PNR_beta_parameters *betaParams    /**< beta parameter struct */
  )
  {
    // Save stuff into a file
    FILE *fileparams;
    char fileSpec[40];
    sprintf(fileSpec, "pnr_angle_parameters.dat");

    fileparams = fopen(fileSpec, "w");

    fprintf(fileparams, "~~~~~~ Alpha Coefficients ~~~~~~\n");
    fprintf(fileparams, "A1 = %.16e\n", alphaParams->A1);
    fprintf(fileparams, "A2 = %.16e\n", alphaParams->A2);
    fprintf(fileparams, "A3 = %.16e\n", alphaParams->A3);
    fprintf(fileparams, "A4 = %.16e\n\n", alphaParams->A4);

    fprintf(fileparams, "~~~~~~ Beta Coefficients ~~~~~~\n");
    fprintf(fileparams, "B0 = %.16e\n", betaParams->B0);
    fprintf(fileparams, "B1 = %.16e\n", betaParams->B1);
    fprintf(fileparams, "B2 = %.16e\n", betaParams->B2);
    fprintf(fileparams, "B3 = %.16e\n", betaParams->B3);
    fprintf(fileparams, "B4 = %.16e\n", betaParams->B4);
    fprintf(fileparams, "B5 = %.16e\n\n", betaParams->B5);

    fprintf(fileparams, "~~~~~~ Connection Values Alpha ~~~~~~\n");
    fprintf(fileparams, "Mf_alpha_lower = %.16e\n", alphaParams->Mf_alpha_lower);
    fprintf(fileparams, "Mf_alpha_upper = %.16e\n", alphaParams->Mf_alpha_upper);
    fprintf(fileparams, "alpha_lower = %.16e\n", alphaParams->alpha_lower);
    fprintf(fileparams, "alpha_upper = %.16e\n", alphaParams->alpha_upper);
    fprintf(fileparams, "derivative_alpha_lower = %.16e\n", alphaParams->derivative_alpha_lower);
    fprintf(fileparams, "derivative_alpha_upper = %.16e\n", alphaParams->derivative_alpha_upper);
    fprintf(fileparams, "alpha_interp_0 = %.16e\n", alphaParams->alpha_interp_0);
    fprintf(fileparams, "alpha_interp_1 = %.16e\n", alphaParams->alpha_interp_1);
    fprintf(fileparams, "alpha_interp_2 = %.16e\n", alphaParams->alpha_interp_2);
    fprintf(fileparams, "alpha_interp_3 = %.16e\n\n", alphaParams->alpha_interp_3);

    fprintf(fileparams, "~~~~~~ Connection Values Beta ~~~~~~\n");
    fprintf(fileparams, "Mf_beta_lower = %.16e\n", betaParams->Mf_beta_lower);
    fprintf(fileparams, "Mf_beta_upper = %.16e\n", betaParams->Mf_beta_upper);
    fprintf(fileparams, "beta_lower = %.16e\n", betaParams->beta_lower);
    fprintf(fileparams, "beta_upper = %.16e\n", betaParams->beta_upper);
    fprintf(fileparams, "derivative_beta_lower = %.16e\n", betaParams->derivative_beta_lower);
    fprintf(fileparams, "derivative_beta_upper = %.16e\n", betaParams->derivative_beta_upper);
    fprintf(fileparams, "beta_rescale_1 = %.16e\n", betaParams->beta_rescale_1);
    fprintf(fileparams, "beta_rescale_2 = %.16e\n\n", betaParams->beta_rescale_2);

    fclose(fileparams);
  }


  /* Compute window function to ensure smooth transition to PN expression for angles outside calibration. */
  REAL8 IMRPhenomX_PNR_AnglesWindow(
				    IMRPhenomXWaveformStruct *pWF,
				    IMRPhenomXPrecessionStruct *pPrec
				    ){

    double window_q_boundary   = 8.5;
    double window_chi_boundary = 0.85;

    double window_width_q   = 3.5;
    double window_width_chi = 0.35;

    double window_q_argument   = (pWF->q - window_q_boundary ) / window_width_q;
    double window_chi_argument = (pPrec->chi_singleSpin - window_chi_boundary) / window_width_chi;

    // Window mass ratio
    double window_q_value;
    if ( (window_q_argument>0) && (window_q_argument<=1) ){
      window_q_value = 0.5*cos( window_q_argument     * LAL_PI ) + 0.5;
    } else if ( window_q_argument>1 ) {
      window_q_value = 0;
    } else {
      window_q_value = 1;
    } 

    // Window spin magnitude
    double window_chi_value;
    if ( (window_chi_argument>0) && (window_chi_argument<=1) ){
      window_chi_value = 0.5*cos( window_chi_argument     * LAL_PI ) + 0.5;
    } else if ( window_chi_argument>1  ) {
      window_chi_value = 0;
    } else {
      window_chi_value = 1;
    }
    
    /* the final window is the product of the individiual windows */
    double angles_window  = window_q_value * window_chi_value;

    return angles_window;

  }
    
  
  /* Compute window function which controls use of PNR coprecessing deviations. NOTE that it only makes sense to use this function inside of LALSimIMRPhenomX_precession.c */
  REAL8 IMRPhenomX_PNR_CoprecWindow(
    IMRPhenomXWaveformStruct *pWF
  ){

    // NOTE that it only makes sense to use this function inside of LALSimIMRPhenomX_precession.c

    // Set window parameters based on selected model version
    double window_q_boundary;
    double width_window_q;
      
    // Location of boundaries
    window_q_boundary     = 10.0; // LEGACY VALUE --->  8.0;
    // Width of the transition AFTER the boundary location
    width_window_q        = 10.0; // LEGACY VALUE --->  0.50;
      
    // // Location of boundaries
    // double window_theta_boundary = 150.0*LAL_PI/180.0;
    // double window_a1_boundary    = 0.8; 
    
    // // Width of the transition AFTER the boundary location
    // double width_window_theta = 0.50;
    // double width_window_a1    = 0.02;

    double window_q_argument     = (pWF->q   - window_q_boundary     )/ width_window_q;
    // double window_theta_argument = (pWF->theta_LS - window_theta_boundary )/ width_window_theta;
    // double window_a1_argument    = (pWF->a1       - window_a1_boundary    )/ width_window_a1;

    // When the arguments are <=0, then the model is on
    // when they are between 0 and 1, the model is turning off
    // and when they are greater than 1, the model is off

    // For q
    double window_q_value;
    if ( (window_q_argument>0) && (window_q_argument<=1) ){
      window_q_value = 0.5*cos( window_q_argument     * LAL_PI ) + 0.5;
    } else if ( window_q_argument>1  ) {
      window_q_value = 0;
    } else {
      window_q_value = 1;
    } //

    // // For theta
    // double window_theta_value;
    // if ( (window_theta_argument>0) && (window_theta_argument<=1) ){
    //   window_theta_value = 0.5*cos( window_theta_argument     * LAL_PI ) + 0.5;
    // } else if ( window_theta_argument>1  ) {
    //   window_theta_value = 0;
    // } else {
    //   window_theta_value = 1;
    // } //

    // // For a1
    // double window_a1_value;
    // if ( (window_a1_argument>0) && (window_a1_argument<=1) ){
    //   window_a1_value = 0.5*cos( window_a1_argument     * LAL_PI ) + 0.5;
    // } else if ( window_a1_argument>1  ) {
    //   window_a1_value = 0;
    // } else {
    //   window_a1_value = 1;
    // } //

    // the NET window will be the product of the individual windows so that any one parameter may turn the model off. We're only interested in this if PNR is desired via the pflag option.
    double pnr_window = window_q_value;
    // double pnr_window = window_q_value * window_theta_value * window_a1_value;

    //
    return pnr_window;

  }

  /* Compute XAS phase and phase derivative a reference frequency "f_inspiral_align" */
  INT4 IMRPhenomX_PNR_SetPhaseAlignmentParams(
      IMRPhenomXWaveformStruct *pWF,
      IMRPhenomXPrecessionStruct *pPrec
  ){


    /*
    Copy the current state of pWF to pPrec for use in PNR+XPHM, where we will similarly want to enforce phase alignment with XHM during inspiral.

    The code immediately below is very typical of annoying C language code:
    to copy the structure, one first must allocate memory for the vessel to
    hold the copy. Then one must copy the struct into that allocated
    momory. While there are "correct" versions of this code that do not
    require use of malloc, these versions essentially copy the pointer, so
    when pWF is changed, so is pWF22AS. We do not want that, so the use of
    malloc is essential.
    */
    IMRPhenomXWaveformStruct *pWF22AS = malloc(sizeof(IMRPhenomXWaveformStruct));
    memcpy( pWF22AS, pWF, sizeof(IMRPhenomXWaveformStruct) );
    pPrec->pWF22AS = pWF22AS;
    /* NOTE that the copy must be deleted upstream in order to avoid a memory leak */

    // /* Dev printing */
    // printf("--)                 pWF->fRING = %f\n",pWF->fRING);

    /* Define alignment frequency in fM (aka Mf). This is the
    frequency at which PNR coprecessing phase and phase
    derivaive will be aligned with corresponding XAS and XHM
    values.  */
    pWF->f_inspiral_align = 0.004;

    // Define an int to hold status values
    UINT4 status;

    // NOTE that just below we compute the non-precessing phase parameters
    // BEFORE any changes are made to pWF -- SO the pWF input must not
    // contain any changes due to precession.
    IMRPhenomXPhaseCoefficients *pPhaseAS;
    pPhaseAS = XLALMalloc(sizeof(IMRPhenomXPhaseCoefficients));
    status   = IMRPhenomXGetPhaseCoefficients(pWF,pPhaseAS);
    // Check for error
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXGetPhaseCoefficients failed.\n");
    // Clean up

    /*
    Below we use IMRPhenomX_FullPhase_22 to somultaneously compute
    the XAS phase and phase derivative at the point of interest.
    */

    /**********************************************************/
    // Initialize holders for the phase and phase derivative
    double phase, dphase;
    // Define the values inside of IMRPhenomX_FullPhase_22
    IMRPhenomX_FullPhase_22(&phase,&dphase,pWF->f_inspiral_align,pPhaseAS,pWF);
    // Store the phase and phase derivative for later use
    pWF->XAS_phase_at_f_inspiral_align = phase;
    pWF->XAS_dphase_at_f_inspiral_align = dphase;//full_dphase_value;
    /**********************************************************/

    /*
    Now, once all other model changes have been made, but before the
    final phase is output in IMRPhenomXASGenerateFD, we want to force
    the PNR CoPrecessing phase and phase derivative to be pWF->XAS_phase_at_f_inspiral_align and pWF->XAS_dphase_at_f_inspiral_align, respectively. This effort
    is facilitated by IMRPhenomX_PNR_EnforceXASPhaseAlignment below.
    */

    // // Printing for development
    // printf("##>> XAS_phase_at_f_inspiral_align = %f\n",pWF->XAS_phase_at_f_inspiral_align);
    // printf("##>> XAS_dphase_at_f_inspiral_align = %f\n",pWF->XAS_dphase_at_f_inspiral_align);
    

    LALFree(pPhaseAS);

    //
    return status;

  }



  /* Compute  XHM phase and phase derivative a reference frequency "f_inspiral_align" */
  INT4 IMRPhenomXHM_PNR_SetPhaseAlignmentParams(
      INT4 ell,
      INT4 emm,
      IMRPhenomXWaveformStruct *pWF,
      IMRPhenomXPrecessionStruct *pPrec,
      LALDict *lalParams
  ){

    // Define an int to hold status values
    UINT4 status;

    //
    double f_inspiral_align = pWF->f_inspiral_align;

    //
    double phase,dphase;
    status = IMRPhenomXHM_Phase_for_Initialization(&phase,&dphase,f_inspiral_align,ell,emm,pPrec->pWF22AS,lalParams);

    //
    pWF->XHM_phase_at_f_inspiral_align  = phase;
    pWF->XHM_dphase_at_f_inspiral_align = dphase;

    // // Dev printing
    // printf("(lal)>> XHM_phase_at_f_inspiral_align  = %1.8f\n",phase);
    // printf("(lal)>> XHM_dphase_at_f_inspiral_align = %1.8f\n",dphase);

    //
    return status;

  }



  /* Align the PNR CoPrec phase and phase derivative at
  "f_inspiral_align" by changing the effective value of
  phifRef and linb */
  void IMRPhenomX_PNR_EnforceXASPhaseAlignment(
      double* linb,
      IMRPhenomXWaveformStruct *pWF,
      IMRPhenomXPhaseCoefficients *pPhase
  ){

    /**********************************************************/
    // Initialize holders for the phase and phase derivative
    double phi_at_f_inspiral_align;
    double dphi_at_f_inspiral_align;
    // Define the values inside of IMRPhenomX_FullPhase_22
    IMRPhenomX_FullPhase_22(&phi_at_f_inspiral_align,&dphi_at_f_inspiral_align,pWF->f_inspiral_align,pPhase,pWF);
    /**********************************************************/

    /* If dphi is the phase derivative at some frequency,
       then we want

       dphi = dphi - dphi_at_f_inspiral_align +
       pWF->XAS_dphase_at_f_inspiral_align

       So the constant of interest is ...
    */
    double shift_in_linb = - dphi_at_f_inspiral_align +
       pWF->XAS_dphase_at_f_inspiral_align;

    // So the new value of linb is
    *linb += shift_in_linb;

    /* NOTE that phifRef need not be changed due to its dependence on linb.
    In other words, given the new value of linb, we leave the computation of
    phifRef to an external routine. */

    // /* If phi is the phase at some frequency, then we want
    //    phi = phi - phi_at_f_inspiral_align +
    //    pWF->XAS_phase_at_f_inspiral_align
    //    So the constant of interest is ...
    // */
    // double shift_in_phifRef = - phi_at_f_inspiral_align +
    //    pWF->XAS_phase_at_f_inspiral_align;
    // // So the new value of phifRef is
    // *phifRef += shift_in_phifRef;

  }



  /* Align the PNR HM CoPrec phase and phase derivative at
  "f_inspiral_align" by providing the needed phase and time shifts*/
  void IMRPhenomXHM_PNR_EnforceXHMPhaseAlignment(
      double* lina,
      double* linb,
      INT4 ell,
      INT4 emm,
      IMRPhenomXWaveformStruct *pWF,
      LALDict *lalParams
  ){

    /**********************************************************/
    // Initialize holders for the phase and phase derivative
    double phi_at_f_inspiral_align;
    double dphi_at_f_inspiral_align;
    // Define the values
    IMRPhenomXHM_Phase_for_Initialization(&phi_at_f_inspiral_align,&dphi_at_f_inspiral_align,pWF->f_inspiral_align,ell,emm,pWF,lalParams);
    /**********************************************************/

    // printf("**  phi_at_f_inspiral_align = %f\n",phi_at_f_inspiral_align);
    // printf("** dphi_at_f_inspiral_align = %f\n",dphi_at_f_inspiral_align);

    /* If phi is the phase at some frequency, then we want
       phi = phi - phi_at_f_inspiral_align +
       pWF->XHM_phase_at_f_inspiral_align
       So the constant of interest is ...

       NOTE that at this point in the program flow,
       pWF->XHM_phase_at_f_inspiral_align should be the
       correct value for the (ell,emm) multipole moment.
       See
    */
    double shift_in_phase = - phi_at_f_inspiral_align +
       pWF->XHM_phase_at_f_inspiral_align;

    /* So the value of lina is below. But note that we are
    not done yet -- to get the correct answer, we must take
    into account the effect of adding a non-zero phase
    derivative shift. For this see below ... */
    *lina = shift_in_phase;

    /* If dphi is the phase derivative at some frequency,
       then we want

       dphi = dphi - dphi_at_f_inspiral_align +
       pWF->XHM_dphase_at_f_inspiral_align

       So the constant of interest is ...

       NOTE that at this point in the program flow,
       pWF->XHM_phase_at_f_inspiral_align should be the
       correct value for the (ell,emm) multipole moment.
       See
    */
    double shift_in_linb = - dphi_at_f_inspiral_align +
       pWF->XHM_dphase_at_f_inspiral_align;

    // So the value of linb is
    *linb = shift_in_linb;

    /* With the value of linb in hand, we now want to note that
    given the phase, we will add a line:

    phase += lina + linb*Mf .

    when Mf=f_inspiral_align, then we want alignment to hold.
    For this to hapen, lina must include the negative of
    linb*Mf at f_inspiral_align:

    So the new value of lina is */
    *lina = shift_in_phase - shift_in_linb*pWF->f_inspiral_align;

  }

/* Function to get and or store coprec params
into pWF and pPrec */
INT4 IMRPhenomX_PNR_GetAndSetCoPrecParams(
    IMRPhenomXWaveformStruct *pWF,
    IMRPhenomXPrecessionStruct *pPrec,
    LALDict *lalParams
){


  //
  INT4 status = 0;

  // Get toggle for outputting coprecesing model from LAL dictionary
  pPrec->IMRPhenomXReturnCoPrec = XLALSimInspiralWaveformParamsLookupPhenomXReturnCoPrec(lalParams);

  // Get toggle for PNR coprecessing tuning
  INT4 PNRUseTunedCoprec = XLALSimInspiralWaveformParamsLookupPhenomXPNRUseTunedCoprec(lalParams);
  pWF->IMRPhenomXPNRUseTunedCoprec   = PNRUseTunedCoprec;
  pPrec->IMRPhenomXPNRUseTunedCoprec = PNRUseTunedCoprec;
  // Same as above but for 33
  pPrec->IMRPhenomXPNRUseTunedCoprec33 = XLALSimInspiralWaveformParamsLookupPhenomXPNRUseTunedCoprec33(lalParams) * PNRUseTunedCoprec;
  pWF->IMRPhenomXPNRUseTunedCoprec33 = pPrec->IMRPhenomXPNRUseTunedCoprec33;
  
  // Throw error if preferred value of PNRUseTunedCoprec33 is not found
  if ( pPrec->IMRPhenomXPNRUseTunedCoprec33 ) {
    XLAL_ERROR(XLAL_EFUNC,"Error: Coprecessing tuning for l=|m|=3 must be off.\n");
  }

  // Get toggle for enforced use of non-precessing spin as is required during tuning of PNR's coprecessing model
  INT4 PNRUseInputCoprecDeviations = XLALSimInspiralWaveformParamsLookupPhenomXPNRUseInputCoprecDeviations(lalParams);
  pPrec->IMRPhenomXPNRUseInputCoprecDeviations = PNRUseInputCoprecDeviations;

  // Get toggle for forcing inspiral phase and phase derivative alignment with XHM/AS
  INT4 PNRForceXHMAlignment = XLALSimInspiralWaveformParamsLookupPhenomXPNRForceXHMAlignment(lalParams);
  pWF->IMRPhenomXPNRForceXHMAlignment = PNRForceXHMAlignment;
  
  // Throw error if preferred value of PNRForceXHMAlignment is not found
  if ( PNRForceXHMAlignment ) {
    XLAL_ERROR(XLAL_EFUNC,"Error: PNRForceXHMAlignment must be off.\n");
  }

#if DEBUG == 1
  if ( PNRForceXHMAlignment ){
    printf("lal>> Dev toggle for XHM alignment ON in LALSimIMRPhenomX_precession.c.\n");
  }
#endif

  /*-~-~-~-~-~-~-~-~-~-~-~-~-~*
  Validate PNR usage options
  *-~-~-~-~-~-~-~-~-~-~-~-~-~*/
  if (  PNRUseInputCoprecDeviations && PNRUseTunedCoprec ) {
    //
    XLAL_ERROR(XLAL_EDOM,"Error: PNRUseTunedCoprec and PNRUseInputCoprecDeviations must not be enabled simultaneously.\n");
  }

  // Define high-level toggle for whether to apply deviations. NOTE that this is imposed at the definition of PNR_DEV_PARAMETER, rather than in a series of IF-ELSE conditions.
  INT4 APPLY_PNR_DEVIATIONS = PNRUseTunedCoprec || PNRUseInputCoprecDeviations;
  pWF->APPLY_PNR_DEVIATIONS = APPLY_PNR_DEVIATIONS;

  // If applying PNR deviations, then we want to be able to refer to some non-PNR waveform properties. For that, we must compute the struct for when PNR is off (and specifically, XAS is wanted).
  if ( APPLY_PNR_DEVIATIONS && PNRForceXHMAlignment ) {

    /*<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.<-.
    Compute and store the value of the XAS phase derivative at pPrec->f_inspiral_align
      - Note that this routine also copies the current state of pWF to pPrec for use in PNR+XPHM, where we will similarly want to enforce phase alignment with XHM during inspiral.
    ->.->.->.->.->.->.->.->.->.->.->.->.->.->.->.->.->.->.->.->.->.->.->.->.->.->.*/
    IMRPhenomX_PNR_SetPhaseAlignmentParams(pWF,pPrec);

  }

  /*-~---~---~---~---~---~---~---~---~---~---~---~---~--*/
  /*  Define single spin parameters for fit evaluation  */
  /*-~---~---~---~---~---~---~---~---~---~---~---~---~--*/

  //
  REAL8 a1 = pPrec->chi_singleSpin;
  pWF->a1 = a1;
  REAL8 cos_theta = pPrec->costheta_singleSpin;

  //
  double theta_LS = acos( cos_theta );
  pWF->theta_LS = theta_LS;


  // Use external function to compute window of tuning deviations. The value of the window will only differ from unity if PNRUseTunedCoprec is equipotent to True. NOTE that correct evaluation of this function requires that e,g, pWF->a1 and pWF->theta_LS be defined above.
  double pnr_window = 0.0; /* Making the defualt to be zero here, meaning that by default tunings will be off regardless of physical case, or other option flags.*/
  if (PNRUseTunedCoprec) {
    // Store for output in related XLAL function
    pnr_window = IMRPhenomX_PNR_CoprecWindow(pWF);
  }
  pWF->pnr_window = pnr_window;


  /* Store XCP deviation parameter: NOTE that we only want to apply the window if PNR is being used, not e.g. if we are calibrating the related coprecessing model */
  pWF->PNR_DEV_PARAMETER = a1 * sin( pWF->theta_LS ) * APPLY_PNR_DEVIATIONS;
  if ( PNRUseTunedCoprec ){
    pWF->PNR_DEV_PARAMETER = pnr_window * (pWF->PNR_DEV_PARAMETER);
    // NOTE that PNR_DEV_PARAMETER for l=m=3 is derived from (and directly proportional to) the one defined just above.
  }

  /* Store deviations to be used in PhenomXCP (PNRUseInputCoprecDeviations) */
  // Get them from the laldict (also used as a way to get default values)
  // For information about how deviations are applied, see code chunk immediately below.
  /* NOTE the following for the code just below:
      - all default values are zero
      - we could toggle the code chunk with PNRUseInputCoprecDeviations, but doing so would be non-orthogonal to the comment above about default values.
      - In any case, the user must set PNRUseInputCoprecDeviations=True, AND manually set the deviations using the LALDict interface.
  */
  pWF->MU1   = XLALSimInspiralWaveformParamsLookupPhenomXCPMU1(lalParams);
  pWF->MU2   = XLALSimInspiralWaveformParamsLookupPhenomXCPMU2(lalParams);
  pWF->MU3   = XLALSimInspiralWaveformParamsLookupPhenomXCPMU3(lalParams);
  pWF->MU4   = XLALSimInspiralWaveformParamsLookupPhenomXCPMU4(lalParams);
  pWF->NU0   = XLALSimInspiralWaveformParamsLookupPhenomXCPNU0(lalParams);
  pWF->NU4   = XLALSimInspiralWaveformParamsLookupPhenomXCPNU4(lalParams);
  pWF->NU5   = XLALSimInspiralWaveformParamsLookupPhenomXCPNU5(lalParams);
  pWF->NU6   = XLALSimInspiralWaveformParamsLookupPhenomXCPNU6(lalParams);
  pWF->ZETA1 = XLALSimInspiralWaveformParamsLookupPhenomXCPZETA1(lalParams);
  pWF->ZETA2 = XLALSimInspiralWaveformParamsLookupPhenomXCPZETA2(lalParams);

  //
  #if DEBUG == 1
    printf("** >>>>>>>>>>>> PhenomXCP Model domain >>>>>>>>>>> **\n");
    printf("theta : %f\n",theta_LS*180.0/LAL_PI);
    printf("eta   : %f\n",pWF->eta);
    printf("a1    : %f\n",a1);
    printf("** >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> **\n\n");
  #endif

  //
  if( PNRUseTunedCoprec )
  {

    /* ------------------------------------------------------ >>
     Get them from the stored model fits that define PhenomXCP
     within PhenomXPNR. NOTE that most but not all
     modifications take place in LALSimIMRPhenomX_internals.c.
     For example, fRING and fDAMP are modified in this file.
     NOTE that each tuned parameter requires pWF->PNR_DEV_PARAMETER
     to be unchanged from the value used during tuning e.g. a1*sin(theta)
    << ------------------------------------------------------ */

    // Flatten mass-ratio dependence to limit extrapolation artifacts outside of calibration region
    double coprec_eta = ( pWF->eta >= 0.09876 ) ? pWF->eta : 0.09876;
    // Flatten spin dependence to limit extrapolation artifacts outside of calibration region
    double coprec_a1 = pWF->a1;
    coprec_a1  = ( coprec_a1  <= 0.8     ) ? coprec_a1  : 0.8;
	  coprec_a1  = ( coprec_a1  >= 0.2     ) ? coprec_a1  : 0.2;


    /* MU1 modifies pAmp->v1RD */
    pWF->MU1     = XLALSimIMRPhenomXCP_MU1_l2m2(   theta_LS, coprec_eta, coprec_a1 );

    // NOTE that the function for MU2 is not defined in the model
    /* MU2 would modify pAmp->gamma2 */

    /* MU2  */
    pWF->MU2     = XLALSimIMRPhenomXCP_MU2_l2m2(   theta_LS, coprec_eta, coprec_a1 );

    /* MU3 modifies pAmp->gamma3 */
    pWF->MU3     = XLALSimIMRPhenomXCP_MU3_l2m2(   theta_LS, coprec_eta, coprec_a1 );

    /* MU4 modifies V2 for the intermediate amplitude
    for the DEFAULT value of IMRPhenomXIntermediateAmpVersion
    use in IMRPhenomXPHM */
    // pWF->MU4     = XLALSimIMRPhenomXCP_MU4_l2m2(   theta_LS, coprec_eta, coprec_a1 );

    /* NU0 modifies the output of IMRPhenomX_TimeShift_22() */
    pWF->NU0     = XLALSimIMRPhenomXCP_NU0_l2m2(   theta_LS, coprec_eta, coprec_a1 );

    /* NU4 modifies pPhase->cL */
    pWF->NU4     = XLALSimIMRPhenomXCP_NU4_l2m2(   theta_LS, coprec_eta, coprec_a1 );

    /* NU5 modifies pWF->fRING [EXTRAP-PASS-TRUE] */
    pWF->NU5     = XLALSimIMRPhenomXCP_NU5_l2m2(   theta_LS, coprec_eta, coprec_a1 );

    /* NU6 modifies pWF->fDAMP [EXTRAP-PASS-TRUE] */
    pWF->NU6     = XLALSimIMRPhenomXCP_NU6_l2m2(   theta_LS, coprec_eta, coprec_a1 );

    /* ZETA1 modifies pPhase->b4 */
    pWF->ZETA1   = XLALSimIMRPhenomXCP_ZETA1_l2m2( theta_LS, coprec_eta, coprec_a1 );

    /* ZETA2 modifies pPhase->b1  */
    pWF->ZETA2   = XLALSimIMRPhenomXCP_ZETA2_l2m2( theta_LS, coprec_eta, coprec_a1 );

  }

  //
  pWF->NU0 = 0;

  //
  #if DEBUG == 1
    printf("** >>>>> PhenomXCP Model Parameters >>>>> **\n");
    printf("MU1 : %f\n",pWF->MU1);
    printf("MU2 : %f\n",pWF->MU2);
    printf("MU3 : %f\n",pWF->MU3);
    printf("MU4 : %f\n",pWF->MU4);
    printf("NU4 : %f\n",pWF->NU4);
    printf("NU5 : %f\n",pWF->NU5);
    printf("NU6 : %f\n",pWF->NU6);
    printf("ZETA1 : %f\n",pWF->ZETA1);
    printf("ZETA2 : %f\n",pWF->ZETA2);
    printf("** >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> **\n\n");
  #endif

  //
  return status;

}





#ifdef __cplusplus
}
#endif
