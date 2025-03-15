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

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/XLALError.h>

#include "LALSimIMRPhenomX_PNR_alpha.h"

#ifndef _OPENMP
#define omp ignore
#endif

#ifndef PHENOMXHMDEBUG
#define DEBUG 0
#else
#define DEBUG 1
#endif

  /**
   * This function evaluates Eq. 59 of arXiv:2107.08876.
   */
  REAL8 IMRPhenomX_PNR_GeneratePNRAlphaAtMf(
      REAL8 Mf,                                           /**< geometric frequency */
      const IMRPhenomX_PNR_alpha_parameters *alphaParams, /**< Alpha parameter struct */
      IMRPhenomXWaveformStruct *pWF,                      /**< PhenomX Waveform struct */
      IMRPhenomXPrecessionStruct *pPrec                   /**< PhenomX Precession struct */
  )
  {
    /* Alpha transition frequencies */
    REAL8 Mf_alpha_lower = alphaParams->Mf_alpha_lower;
    REAL8 Mf_alpha_upper = alphaParams->Mf_alpha_upper;

    /* Continuity offset between intermediate and MR region */
    REAL8 alpha_MR_offset = alphaParams->alpha_MR_offset;

    if (Mf <= Mf_alpha_lower)
    {
      /* Below low-frequency cutoff, use either NNLO or MSA alpha depending on precession flag */
      return IMRPhenomX_PNR_GetPNAlphaAtFreq(Mf, pWF, pPrec);
    }

    if (Mf >= Mf_alpha_upper)
    {
      /* above high-frequency cutoff, use MR expression for alpha shifted for continuity */
      return IMRPhenomX_PNR_MR_alpha_expression(Mf, alphaParams) + alpha_MR_offset;
    }

    /* evaluate interpolating function for alpha between Mf_alpha_lower and Mf_alpha_upper */
    return IMRPhenomX_PNR_intermediate_alpha_expression(Mf, alphaParams);
  }

  /**
   * This function generates the blended PNR and PN expressions for alpha for the transition region of parameter space.
   */
  REAL8 IMRPhenomX_PNR_GenerateMergedPNRAlphaAtMf(
      REAL8 Mf,                                           /**< geometric frequency */
      const IMRPhenomX_PNR_alpha_parameters *alphaParams, /**< Alpha parameter struct */
      IMRPhenomXWaveformStruct *pWF,                      /**< PhenomX Waveform struct */
      IMRPhenomXPrecessionStruct *pPrec                   /**< PhenomX Precession struct */
  )
  {
    /* evaluate blending window */
    double pnr_window = IMRPhenomX_PNR_AnglesWindow(pWF, pPrec);
    double msa_window = 1 - pnr_window;

    /* Alpha transition frequencies */
    REAL8 Mf_alpha_lower = alphaParams->Mf_alpha_lower;
    REAL8 Mf_alpha_upper = alphaParams->Mf_alpha_upper;

    /* Continuity offset between intermediate and MR region */
    REAL8 alpha_MR_offset = alphaParams->alpha_MR_offset;

    if (Mf <= Mf_alpha_lower)
    {
      /* Below low-frequency cutoff, use either NNLO or MSA alpha depending on precession flag */
      return IMRPhenomX_PNR_GetPNAlphaAtFreq(Mf, pWF, pPrec);
    }

    if (Mf >= Mf_alpha_upper)
    {
      /* above high-frequency cutoff, use MR expression for alpha shifted for continuity */
      REAL8 alpha_PN = IMRPhenomX_PNR_GetPNAlphaAtFreq(Mf, pWF, pPrec);
      REAL8 alpha_MR = IMRPhenomX_PNR_MR_alpha_expression(Mf, alphaParams) + alpha_MR_offset;
      return pnr_window * alpha_MR + msa_window * alpha_PN;
    }

    /* evaluate interpolating function for alpha between Mf_alpha_lower and Mf_alpha_upper */
    REAL8 alpha_PN = IMRPhenomX_PNR_GetPNAlphaAtFreq(Mf, pWF, pPrec);
    REAL8 alpha_IM = IMRPhenomX_PNR_intermediate_alpha_expression(Mf, alphaParams);
    return pnr_window * alpha_IM + msa_window * alpha_PN;
  }


  /**
   * Here we write a wrapper function to produce either MSA or NNLO
   * alpha for use in IMRPhenomX_PNR_GeneratePNRAlphaAtMf. This function will
   * toggle which to use based on IMRPhenomXPrecVersion.
   */
  REAL8 IMRPhenomX_PNR_GetPNAlphaAtFreq(
      REAL8 Mf,                         /**< geometric frequency */
      IMRPhenomXWaveformStruct *pWF,    /**< PhenomX waveform struct */
      IMRPhenomXPrecessionStruct *pPrec /**< PhenomX precession struct */
  )
  {
    /* Define PN expansion parameter v = (pi M f)^(1/3)*/
    const double omega = LAL_PI * Mf;
    const double omega_cbrt = cbrt(omega);

    REAL8 alpha;

    /* Toggle between NNLO and MSA alpha */
    switch (pPrec->IMRPhenomXPrecVersion)
    {
    /* ~~~~~ Use NNLO PN Euler Angles - Appendix G of arXiv:2004.06503 and https://dcc.ligo.org/LIGO-T1500602 ~~~~~ */
    case 101:
    case 102:
    case 103:
    case 104:
    {
      const double omega_cbrt2 = omega_cbrt * omega_cbrt;
      const double logomega = log(omega);

      alpha = IMRPhenomX_PN_Euler_alpha_NNLO(pPrec, omega, omega_cbrt2, omega_cbrt, logomega);

      break;
    }
    case 220:
    case 221:
    case 222:
    case 223:
    case 224:
    {
      vector vangles = {0., 0., 0.};

      /* ~~~~~ Euler Angles from Chatziioannou et al, PRD 95, 104004, (2017), arXiv:1703.03967 ~~~~~ */
      vangles = IMRPhenomX_Return_phi_zeta_costhetaL_MSA(omega_cbrt, pWF, pPrec);

      alpha = vangles.x;

      break;
    }
    case 330:
    {
      alpha = alpha_SpinTaylor_IMR(Mf,pWF,pPrec);
      if(isnan(alpha) || isinf(alpha)) XLAL_ERROR(XLAL_EDOM, "Error in %s: alpha_SpinTaylor_IMR returned invalid value.\n",__func__);
      break;
    }

    default:
    {
      XLAL_ERROR(XLAL_EINVAL, "Error: IMRPhenomXPrecessionVersion not recognized in IMRPhenomX_PNR_GetPNAlphaAtFreq.\n");
      break;
    }
    }

    return alpha;
  }

  /**
   * This function evaluates the coefficients outlined in
   * Sec 7A of arXiv:2107.08876 for alpha.
   */
  int IMRPhenomX_PNR_precompute_alpha_coefficients(
      IMRPhenomX_PNR_alpha_parameters *alphaParams, /**< Alpha parameter struct */
      IMRPhenomXWaveformStruct *pWF,                /**< PhenomX waveform struct */
      IMRPhenomXPrecessionStruct *pPrec             /**< PhenomX precession struct */
  )
  {
    /* safety first */
    XLAL_CHECK(alphaParams != NULL, XLAL_EFAULT);
    XLAL_CHECK(pWF != NULL, XLAL_EFAULT);
    XLAL_CHECK(pPrec != NULL, XLAL_EFAULT);

    /* Grab the single-spin MR parameters */
    REAL8 eta = pWF->eta;
    if( pPrec->IMRPhenomXPrecVersion==330 ){
      eta = ( pWF->eta >= 0.09 ) ? pPrec->eta : 0.09;
    }
    else{
      eta = pWF->eta;
    }
    REAL8 chiboundary = 0.80 - 0.20 * exp( -pow((pWF->q - 6.0)/1.5, 8) );
    REAL8 chi;
    if( pPrec->IMRPhenomXPrecVersion==330 ){
      chi = ( pPrec->chi_singleSpin <= chiboundary ) ? pPrec->chi_singleSpin : chiboundary;
    }
    else{
      chi = pPrec->chi_singleSpin;
    }
    REAL8 costheta = pPrec->costheta_singleSpin;

    /* Evaluate coefficients as described in Sec. 8D in arXiv:2107.08876 */
    REAL8 A2 = IMRPhenomX_PNR_alpha_A2_coefficient(eta, chi, costheta);
    if (A2 > 0.0)
    {
      A2 = 0.0;
    }
    REAL8 A3 = fabs(IMRPhenomX_PNR_alpha_A3_coefficient(eta, chi, costheta));
    if (A3 < 0.00001)
    {
      A3 = 0.00001;
    }
    REAL8 maxA2 = LAL_PI*LAL_PI*sqrt(A3);
    if (A2 < -maxA2)
      {
        A2 = -maxA2;
      }

    alphaParams->A1 = fabs(IMRPhenomX_PNR_alpha_A1_coefficient(eta, chi, costheta));
    alphaParams->A2 = A2;
    alphaParams->A3 = A3;
    alphaParams->A4 = IMRPhenomX_PNR_alpha_A4_coefficient(eta, chi, costheta);

    return XLAL_SUCCESS;
  }

  /**
   * The following four functions are used to compute the C1 interpolation function
   * used to connect the inspiral and MR expressions for alpha. They are evaluated at
   * the two transition frequency values Mf1 = 2 A4 / 7 and Mf2 = A4 / 3.
   * See Sec. 8A of arXiv:2107.08876.
   *
   * Specifically, the next four functions are Eq. 52 in arXiv:2107.08876.
   */
  REAL8 IMRPhenomX_PNR_alpha_interpolation_0(
      REAL8 Mf1,     /**< lower connection frequency in geometric units */
      REAL8 Mf2,     /**< upper connection frequency in geometric units */
      REAL8 alpha1,  /**< value of alpha at Mf1 */
      REAL8 alpha2,  /**< value of alpha at Mf2 */
      REAL8 dalpha1, /**< derivative of alpha at Mf1 */
      REAL8 dalpha2  /**< derivative of alpha at Mf2 */
  )
  {

    REAL8 D = (Mf2 - Mf1) * (Mf2 - Mf1) * (Mf2 - Mf1);
    REAL8 N = 2.0 * (Mf1 * alpha1 - Mf2 * alpha2) - (Mf1 - Mf2) * ((Mf1 * dalpha1 + Mf2 * dalpha2) + (alpha1 + alpha2));

    return N / D;
  }

  REAL8 IMRPhenomX_PNR_alpha_interpolation_1(
      REAL8 Mf1,     /**< lower connection frequency in geometric units */
      REAL8 Mf2,     /**< upper connection frequency in geometric units */
      REAL8 alpha1,  /**< value of alpha at Mf1 */
      REAL8 alpha2,  /**< value of alpha at Mf2 */
      REAL8 dalpha1, /**< derivative of alpha at Mf1 */
      REAL8 dalpha2  /**< derivative of alpha at Mf2 */
  )
  {

    REAL8 D = (Mf2 - Mf1) * (Mf2 - Mf1) * (Mf2 - Mf1);
    REAL8 N = 3.0 * (Mf1 + Mf2) * (Mf2 * alpha2 - Mf1 * alpha1) + (Mf1 - Mf2) * ((Mf1 + 2.0 * Mf2) * (Mf1 * dalpha1 + alpha1) + (2.0 * Mf1 + Mf2) * (Mf2 * dalpha2 + alpha2));

    return N / D;
  }

  REAL8 IMRPhenomX_PNR_alpha_interpolation_2(
      REAL8 Mf1,     /**< lower connection frequency in geometric units */
      REAL8 Mf2,     /**< upper connection frequency in geometric units */
      REAL8 alpha1,  /**< value of alpha at Mf1 */
      REAL8 alpha2,  /**< value of alpha at Mf2 */
      REAL8 dalpha1, /**< derivative of alpha at Mf1 */
      REAL8 dalpha2  /**< derivative of alpha at Mf2 */
  )
  {

    REAL8 D = (Mf2 - Mf1) * (Mf2 - Mf1) * (Mf2 - Mf1);
    REAL8 N = 6.0 * Mf1 * Mf2 * (Mf1 * alpha1 - Mf2 * alpha2) - (Mf1 - Mf2) * (Mf2 * (2.0 * Mf1 + Mf2) * (Mf1 * dalpha1 + alpha1) + Mf1 * (Mf1 + 2.0 * Mf2) * (Mf2 * dalpha2 + alpha2));

    return N / D;
  }

  REAL8 IMRPhenomX_PNR_alpha_interpolation_3(
      REAL8 Mf1,     /**< lower connection frequency in geometric units */
      REAL8 Mf2,     /**< upper connection frequency in geometric units */
      REAL8 alpha1,  /**< value of alpha at Mf1 */
      REAL8 alpha2,  /**< value of alpha at Mf2 */
      REAL8 dalpha1, /**< derivative of alpha at Mf1 */
      REAL8 dalpha2  /**< derivative of alpha at Mf2 */
  )
  {

    REAL8 D = (Mf2 - Mf1) * (Mf2 - Mf1) * (Mf2 - Mf1);
    REAL8 N = Mf1 * Mf2 * Mf2 * (Mf2 - 3.0 * Mf1) * alpha1 - Mf1 * Mf1 * Mf2 * (Mf1 - 3.0 * Mf2) * alpha2 + Mf1 * Mf2 * (Mf1 - Mf2) * (Mf2 * (Mf1 * dalpha1 + alpha1) + Mf1 * (Mf2 * dalpha2 + alpha2));

    return N / D;
  }

  /**
   * This function evaluates Eq. 51 of arXiv:2107.08876
   */
  REAL8 IMRPhenomX_PNR_intermediate_alpha_expression(
      REAL8 Mf,                                          /**< frequency in geometric units */
      const IMRPhenomX_PNR_alpha_parameters *alphaParams /**< Alpha parameter struct */
  )
  {
    /* check we have data */
    XLAL_CHECK(alphaParams != NULL, XLAL_EFAULT);

    /* get interpolation coefficients */
    REAL8 a0 = alphaParams->alpha_interp_0;
    REAL8 a1 = alphaParams->alpha_interp_1;
    REAL8 a2 = alphaParams->alpha_interp_2;
    REAL8 a3 = alphaParams->alpha_interp_3;

    /* evaluate */
    return a0 * Mf * Mf + a1 * Mf + a2 + a3 / Mf;
  }

  /**
   * This function evaluates Eq. 48 of arXiv:2107.08876
   */
  REAL8 IMRPhenomX_PNR_MR_alpha_expression(
      REAL8 Mf,                                          /**< frequency in geometric units */
      const IMRPhenomX_PNR_alpha_parameters *alphaParams /**< Alpha parameter struct */
  )
  {
    /* check we have data */
    XLAL_CHECK(alphaParams != NULL, XLAL_EFAULT);

    /* get MR coefficients */
    REAL8 A1 = alphaParams->A1;
    REAL8 A2 = alphaParams->A2;
    REAL8 A3 = alphaParams->A3;
    REAL8 A4 = alphaParams->A4;

    /* evaluate */
    return -(A1 / Mf + (A2 * sqrt(A3)) / (A3 + (Mf - A4) * (Mf - A4)));
  }

  /**
   * This function evaluates the connection frequencies Mf1 and Mf2 detailed in Sec. 8A
   * of arXiv:2107.08876. It then computes the values of alpha and its derivatives
   * at these frequencies to compute the coefficients used for the intermediate
   * interpolation function Eq. 51.
   *
   * The derivatives are computed using simple centered finite differences.
   */
  int IMRPhenomX_PNR_alpha_connection_parameters(
      IMRPhenomX_PNR_alpha_parameters *alphaParams, /**< Alpha parameter struct */
      IMRPhenomXWaveformStruct *pWF,                /**< PhenomX waveform struct */
      IMRPhenomXPrecessionStruct *pPrec             /**< PhenomX precession struct */
  )
  {
    /* safety first */
    XLAL_CHECK(alphaParams != NULL, XLAL_EFAULT);
    XLAL_CHECK(pWF != NULL, XLAL_EFAULT);
    XLAL_CHECK(pPrec != NULL, XLAL_EFAULT);

    /* define frequency spacing over which to calculate derivatives */
    /* here set to be 0.0005Mf */
    REAL8 dMf = 0.0005;

    /* define connection frequencies for alpha */
    REAL8 Mf_alpha_upper = alphaParams->A4 / 3.0;
    REAL8 Mf_alpha_lower = (3.0 / 3.5) * Mf_alpha_upper;

    /* catch cases where Mf_alpha_lower < 2 Hz and just use MSA alpha */
    if (Mf_alpha_upper < XLALSimIMRPhenomXUtilsHztoMf(2.0, pWF->Mtot))
    {
      Mf_alpha_lower = 100.0;
      Mf_alpha_upper = 100.0;
    }

    /* catch cases where Mf_alpha_lower lies below the starting frequency of the integration of the PN spin-precession eqs */
    if (( pPrec->IMRPhenomXPrecVersion==330)&&(Mf_alpha_lower-2.*dMf < pPrec->Mfmin_integration) )
    {
      Mf_alpha_lower = 100.0;
      Mf_alpha_upper = 100.0;
    }

    alphaParams->Mf_alpha_lower = Mf_alpha_lower;
    alphaParams->Mf_alpha_upper = Mf_alpha_upper;

    /* evaluate expressions for alpha at and around the connection frequencies */
    REAL8 a1 = IMRPhenomX_PNR_GetPNAlphaAtFreq(Mf_alpha_lower - dMf, pWF, pPrec);
    REAL8 alpha_lower = IMRPhenomX_PNR_GetPNAlphaAtFreq(Mf_alpha_lower, pWF, pPrec);
    REAL8 a3 = IMRPhenomX_PNR_GetPNAlphaAtFreq(Mf_alpha_lower + dMf, pWF, pPrec);

    REAL8 a4 = IMRPhenomX_PNR_MR_alpha_expression(Mf_alpha_upper - dMf, alphaParams);
    REAL8 alpha_upper = IMRPhenomX_PNR_GetPNAlphaAtFreq(Mf_alpha_upper, pWF, pPrec);
    REAL8 a6 = IMRPhenomX_PNR_MR_alpha_expression(Mf_alpha_upper + dMf, alphaParams);

    /* calculate derivative of alpha at connection frequencies using finite difference central difference method */
    REAL8 derivative_alpha_lower = (a3 - a1) / (2.0 * dMf);
    REAL8 derivative_alpha_upper = (a6 - a4) / (2.0 * dMf);

    /* evaluate the interpolation coefficients */
    REAL8 alpha_interp_0 = IMRPhenomX_PNR_alpha_interpolation_0(Mf_alpha_lower, Mf_alpha_upper, alpha_lower, alpha_upper, derivative_alpha_lower, derivative_alpha_upper);
    REAL8 alpha_interp_1 = IMRPhenomX_PNR_alpha_interpolation_1(Mf_alpha_lower, Mf_alpha_upper, alpha_lower, alpha_upper, derivative_alpha_lower, derivative_alpha_upper);
    REAL8 alpha_interp_2 = IMRPhenomX_PNR_alpha_interpolation_2(Mf_alpha_lower, Mf_alpha_upper, alpha_lower, alpha_upper, derivative_alpha_lower, derivative_alpha_upper);
    REAL8 alpha_interp_3 = IMRPhenomX_PNR_alpha_interpolation_3(Mf_alpha_lower, Mf_alpha_upper, alpha_lower, alpha_upper, derivative_alpha_lower, derivative_alpha_upper);

    alpha_interp_0 = isnan(alpha_interp_0) ? 0.0 : alpha_interp_0;
    alpha_interp_1 = isnan(alpha_interp_1) ? 0.0 : alpha_interp_1;
    alpha_interp_2 = isnan(alpha_interp_2) ? 0.0 : alpha_interp_2;
    alpha_interp_3 = isnan(alpha_interp_3) ? 0.0 : alpha_interp_3;

    /* calculate the offset required for the MR contributions to alpha to continuously
     * connect to the inspiral. */
    REAL8 MR_alpha_at_Mf_upper = IMRPhenomX_PNR_MR_alpha_expression(Mf_alpha_upper, alphaParams);
    REAL8 alpha_MR_offset = alpha_upper - MR_alpha_at_Mf_upper;

    /* save alpha values into the struct */
    alphaParams->alpha_lower = alpha_lower;
    alphaParams->alpha_upper = alpha_upper;
    alphaParams->derivative_alpha_lower = derivative_alpha_lower;
    alphaParams->derivative_alpha_upper = derivative_alpha_upper;
    alphaParams->alpha_interp_0 = alpha_interp_0;
    alphaParams->alpha_interp_1 = alpha_interp_1;
    alphaParams->alpha_interp_2 = alpha_interp_2;
    alphaParams->alpha_interp_3 = alpha_interp_3;
    alphaParams->alpha_MR_offset = alpha_MR_offset;

    return XLAL_SUCCESS;
  }

#ifdef __cplusplus
}
#endif
