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
   * \author Eleanor Hamilton, Sebastian Khan, Jonathan E. Thompson, Marta Colleoni
   *
   */

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/XLALError.h>

#include "LALSimIMRPhenomX_PNR_beta.h"

#ifndef _OPENMP
#define omp ignore
#endif

#ifndef PHENOMXHMDEBUG
#define DEBUG 0
#else
#define DEBUG 1
#endif

  /**
   * This function evaluates Eqs. 60 and 61 of arXiv:2107.08876.
   */
  REAL8 IMRPhenomX_PNR_GeneratePNRBetaAtMf(
      REAL8 Mf,                                         /**< geometric frequency */
      const IMRPhenomX_PNR_beta_parameters *betaParams, /**< beta parameter struct */
      IMRPhenomXWaveformStruct *pWF,                    /**< PhenomX waveform struct */
      IMRPhenomXPrecessionStruct *pPrec,                /**< PhenomX precession struct */
      IMRPhenomXWaveformStruct *pWF_SingleSpin,         /**< PhenomX waveform struct with approximate single spin */
      IMRPhenomXPrecessionStruct *pPrec_SingleSpin      /**< PhenomX waveform struct with approximate single spin */
  )
  {
    /* get beta connection frequencies */
    REAL8 Mf_beta_lower = betaParams->Mf_beta_lower;
    REAL8 Mf_beta_upper = betaParams->Mf_beta_upper;

    REAL8 beta_lower = betaParams->beta_lower;
    REAL8 beta_upper = betaParams->beta_upper;

    if (Mf <= Mf_beta_lower)
    { /* Below Mf_beta_lower, we use a rescaled form of the PN beta.
       * In the case of a two-spin system, the PN beta used may have the two-spin
       * oscillations tapered away to connect nicely at Mf_beta_lower */

       if ((pPrec->IMRPhenomXPrecVersion==223) && (beta_lower < 0.01 * beta_upper))
       {
          /* Catch cases in the almost anti-aligned spin limit where beta drops
           * to almost zero before reaching Mf_beta_lower */
          REAL8 MR_beta = IMRPhenomX_PNR_MR_beta_expression(Mf_beta_lower, betaParams);
          return IMRPhenomX_PNR_arctan_window(MR_beta);
       }

      REAL8 betaPN = IMRPhenomX_PNR_GetPNBetaAtFreq(Mf, betaParams, pWF, pPrec, pWF_SingleSpin, pPrec_SingleSpin);
      REAL8 beta_waveform = IMRPhenomX_PNR_PNWaveformBetaWrapper(Mf, betaPN, pWF, pPrec);
      REAL8 rescaled_beta = beta_waveform * IMRPhenomX_PNR_rescale_beta_expression(Mf, betaParams);
      return IMRPhenomX_PNR_arctan_window(rescaled_beta);
    }

    if (Mf >= Mf_beta_upper)
    { /* above Mf_beta_upper, attach the final value of beta*/
      REAL8 final_beta = IMRPhenomX_PNR_MR_beta_expression(Mf_beta_upper, betaParams);
      return IMRPhenomX_PNR_arctan_window(final_beta);
    }

    /* in between we evaluate the MR Ansatz */
    REAL8 MR_beta = IMRPhenomX_PNR_MR_beta_expression(Mf, betaParams);
    return IMRPhenomX_PNR_arctan_window(MR_beta);
}


  /**
   * This function evaluates only the rescaled inspiral beta given in
   * Eq. 41 of arXiv:2107.08876, without attaching the MR model or tapering
   * the two-spin oscillations.
   */
  REAL8 IMRPhenomX_PNR_GeneratePNRBetaNoMR(
      REAL8 Mf,                         /**< geometric frequency */
      const IMRPhenomX_PNR_beta_parameters *betaParams, /**< beta parameter struct */
      IMRPhenomXWaveformStruct *pWF,    /**< PhenomX wavefrom struct */
      IMRPhenomXPrecessionStruct *pPrec /**< PhenomX precession struct */
  )
  {

    /* generate scaled PN beta without MR contributions */
    REAL8 waveform_beta;
    REAL8 pn_beta = IMRPhenomX_PNR_GetPNBetaAtFreq_fulltwospin(Mf, pWF, pPrec, betaParams);
    REAL8 pn_waveform_beta = IMRPhenomX_PNR_PNWaveformBetaWrapper(Mf, pn_beta, pWF, pPrec);

    waveform_beta = pn_waveform_beta;

    return IMRPhenomX_PNR_arctan_window(waveform_beta);
  }

  /**
   *  This function generates beta with the tuned angles and PN expressions blended during merger-ringdown
   */
  REAL8 IMRPhenomX_PNR_GenerateMergedPNRBetaAtMf(
      REAL8 Mf,                                         /**< geometric frequency */
      const IMRPhenomX_PNR_beta_parameters *betaParams, /**< beta parameter struct */
      IMRPhenomXWaveformStruct *pWF,                    /**< PhenomX waveform struct */
      IMRPhenomXPrecessionStruct *pPrec,                /**< PhenomX precession struct */
      IMRPhenomXWaveformStruct *pWF_SingleSpin,         /**< PhenomX waveform struct with approximate single spin */
      IMRPhenomXPrecessionStruct *pPrec_SingleSpin      /**< PhenomX waveform struct with approximate single spin */
						  )
  {
    /* evaluate blending window */
    double pnr_window = IMRPhenomX_PNR_AnglesWindow(pWF, pPrec);
    double msa_window = 1-pnr_window;

    /* get beta connection frequencies */
    REAL8 Mf_beta_lower = betaParams->Mf_beta_lower;
    REAL8 Mf_beta_upper = betaParams->Mf_beta_upper;
    if (Mf <= Mf_beta_lower)
      { /* Below Mf_beta_lower, we use a rescaled form of the PN beta.
	 * In the case of a two-spin system, the PN beta used may have the two-spin
	 * oscillations tapered away to connect nicely at Mf_beta_lower */
	REAL8 pnr_beta, pn_beta;
	if (pPrec->IMRPhenomXPrecVersion!=330){
	  REAL8 betaPN = IMRPhenomX_PNR_GetPNBetaAtFreq(Mf, betaParams, pWF, pPrec, pWF_SingleSpin, pPrec_SingleSpin);
	  REAL8 beta_waveform = IMRPhenomX_PNR_PNWaveformBetaWrapper(Mf, betaPN, pWF, pPrec);
	  REAL8 rescaled_beta = beta_waveform * IMRPhenomX_PNR_rescale_beta_expression(Mf, betaParams);
	  pnr_beta = rescaled_beta;
	  pn_beta = beta_waveform;
	}
	else{
	  pPrec->UseMRbeta = 1;
	  REAL8 betaPN1 = IMRPhenomX_PNR_GetPNBetaAtFreq(Mf, betaParams, pWF, pPrec, pWF_SingleSpin, pPrec_SingleSpin);
	  REAL8 beta_waveform1 = IMRPhenomX_PNR_PNWaveformBetaWrapper(Mf, betaPN1, pWF, pPrec);
	  REAL8 rescaled_beta1 = beta_waveform1 * IMRPhenomX_PNR_rescale_beta_expression(Mf, betaParams);
	  pPrec->UseMRbeta = 0;
	  REAL8 betaPN2 = IMRPhenomX_PNR_GetPNBetaAtFreq(Mf, betaParams, pWF, pPrec, pWF_SingleSpin, pPrec_SingleSpin);
	  REAL8 beta_waveform2 = IMRPhenomX_PNR_PNWaveformBetaWrapper(Mf, betaPN2, pWF, pPrec);
	  pnr_beta = rescaled_beta1;
	  pn_beta = beta_waveform2;
	}
	return IMRPhenomX_PNR_arctan_window(pnr_window * pnr_beta + msa_window * pn_beta);
      }

    if (Mf >= Mf_beta_upper)
      { /* above Mf_beta_upper, attach the final value of beta*/
	pPrec->UseMRbeta = 0;
	REAL8 betaPN = IMRPhenomX_PNR_GetPNBetaAtFreq(Mf, betaParams, pWF, pPrec, pWF_SingleSpin, pPrec_SingleSpin);
	REAL8 beta_waveform = IMRPhenomX_PNR_PNWaveformBetaWrapper(Mf, betaPN, pWF, pPrec);
	REAL8 final_beta = IMRPhenomX_PNR_MR_beta_expression(Mf_beta_upper, betaParams);
	return IMRPhenomX_PNR_arctan_window(pnr_window * final_beta + msa_window * beta_waveform);
      }

    /* in between we evaluate the MR Ansatz */
    pPrec->UseMRbeta = 0;
    REAL8 betaPN = IMRPhenomX_PNR_GetPNBetaAtFreq(Mf, betaParams, pWF, pPrec, pWF_SingleSpin, pPrec_SingleSpin);
    REAL8 beta_waveform = IMRPhenomX_PNR_PNWaveformBetaWrapper(Mf, betaPN, pWF, pPrec);
    REAL8 MR_beta = IMRPhenomX_PNR_MR_beta_expression(Mf, betaParams);
    return IMRPhenomX_PNR_arctan_window(pnr_window * MR_beta + msa_window * beta_waveform);
  }
  /**
   * We evaluate beta at the final Mf_beta_upper connection frequency
   * to approximate the final value of beta during ringdown. This
   * is required to analytically approximate the effective ringdown frequency. FIXME: add citation
   */

  REAL8 IMRPhenomX_PNR_GenerateRingdownPNRBeta(
      IMRPhenomXWaveformStruct *pWF,     /**< PhenomX waveform struct */
      IMRPhenomXPrecessionStruct *pPrec) /**< PhenomX precession struct */
  {
    /* may we continue? */
    XLAL_CHECK(pWF != NULL, XLAL_EFAULT);
    XLAL_CHECK(pPrec != NULL, XLAL_EFAULT);

    /* get effective single spin parameters */
    REAL8 eta = pWF->eta;
    REAL8 chi = pPrec->chi_singleSpin;
    REAL8 costheta = pPrec->costheta_singleSpin;

    /* approximate orientation of final spin */
    REAL8 costhetaf = pPrec->costheta_final_singleSpin;

    REAL8 betafinal = IMRPhenomX_PNR_arctan_window( acos(costhetaf) - IMRPhenomX_PNR_beta_Bf_coefficient(eta, chi, costheta) );

    return betafinal;
  }

  /**
   * A wrapper to produce either the NNLO or MSA beta depending on the IMRPhenomXPrecVersion.
   *
   * Should the MSA angle be called, we taper any potential oscillations induced by a time-varying
   * total spin magnitude so that we return an effective single-spin value for beta at the
   * lower connection frequency. This is described in Sec. 6C of arXiv:2107.08876
   */

  REAL8 IMRPhenomX_PNR_GetPNBetaAtFreq(
      REAL8 Mf,                                         /**< geometric frequency */
      const IMRPhenomX_PNR_beta_parameters *betaParams, /**< beta parameter struct */
      IMRPhenomXWaveformStruct *pWF,                    /**< PhenomX waveform struct */
      IMRPhenomXPrecessionStruct *pPrec,                /**< PhenomX precession struct */
      IMRPhenomXWaveformStruct *pWF_SingleSpin,         /**< PhenomX waveform struct with approximate single spin */
      IMRPhenomXPrecessionStruct *pPrec_SingleSpin      /**< PhenomX waveform struct with approximate single spin */
				       )
  {

    REAL8 beta;

    /* get PN expansion parameter v = (pi M f)^(1/3) */
    const double omega = LAL_PI * Mf;
    const double omega_cbrt = cbrt(omega);

    switch (pPrec->IMRPhenomXPrecVersion)
      {
	/* ~~~~~ Use NNLO PN Euler Angles - Appendix G of arXiv:2004.06503 and https://dcc.ligo.org/LIGO-T1500602 ~~~~~ */
      case 101:
      case 102:
      case 103:
      case 104:
	{
	  REAL8 L = XLALSimIMRPhenomXLPNAnsatz(omega_cbrt, pWF->eta / omega_cbrt, pPrec->L0, pPrec->L1, pPrec->L2, pPrec->L3, pPrec->L4, pPrec->L5, pPrec->L6, pPrec->L7, pPrec->L8, pPrec->L8L);

	  /*
	    Comment from IMRPhenomX_precession.c:
	    We ignore the sign of L + SL below:
	    s := Sp / (L + SL)
	  */
	  REAL8 s = pPrec->Sperp / (L + pPrec->SL);
	  REAL8 s2 = s * s;
	  beta = acos(copysign(1.0, L + pPrec->SL) / sqrt(1.0 + s2));

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
	  /* First we grab the full MSA angle with possible two-spin oscillations */
	  vangles = IMRPhenomX_Return_phi_zeta_costhetaL_MSA(omega_cbrt, pWF, pPrec);
	  REAL8 beta_full = acos(vangles.z);

	  /* lower connection frequency as target for taper */
	  REAL8 Mf_beta_lower = betaParams->Mf_beta_lower;

	  /* check if we have meaningful two-spin contributions */
	  if (IMRPhenomX_PNR_CheckTwoSpin(pPrec))
	    { /* we do */
	      /* compute an effective single-spin beta that averages through oscillations */
	      const vector vangles_SingleSpin = IMRPhenomX_Return_phi_zeta_costhetaL_MSA(omega_cbrt, pWF_SingleSpin, pPrec_SingleSpin);
	      REAL8 beta_SingleSpin = acos(vangles_SingleSpin.z);

	      /* if we are below the connection frequency, taper the oscillations */
	      if (Mf <= Mf_beta_lower)
		{
		  /* tapering is described in Eq. 45 of arXiv:2107.08876 */
		  REAL8 oscillations = beta_full - beta_SingleSpin;
		  REAL8 envelope = cos(2.0 * LAL_PI * Mf / (4.0 * Mf_beta_lower)) * cos(2.0 * LAL_PI * Mf / (4.0 * Mf_beta_lower));
		  beta = (beta_SingleSpin + oscillations * envelope);
		}
	      else
		{ /* otherwise just return single-spin beta */
		  beta = beta_SingleSpin;
		}
	    }
	  else
	    { /* no oscillations, just return full beta */
	      beta = beta_full;
	    }

	  break;
	}
      case 330:
	{

	  REAL8 beta_full= beta_SpinTaylor_IMR(Mf,pWF,pPrec,betaParams);
    if(isnan(beta_full) || isinf(beta_full)) XLAL_ERROR(XLAL_EINVAL, "Error in %s: beta_SpinTaylor_IMR returned invalid value.\n",__func__);

	  /* lower connection frequency as target for taper */
	  REAL8 Mf_beta_lower = betaParams->Mf_beta_lower;


	  /* check if we have meaningful two-spin contributions */
	  if (IMRPhenomX_PNR_CheckTwoSpin(pPrec))
	    { /* we do */
	      /* compute an effective single-spin beta that averages through oscillations */
	      const vector vangles_SingleSpin = IMRPhenomX_Return_phi_zeta_costhetaL_MSA(omega_cbrt, pWF_SingleSpin, pPrec_SingleSpin);
	      REAL8 beta_SingleSpin = acos(vangles_SingleSpin.z);

	      /* if we are below the connection frequency, taper the oscillations */
	      if (Mf <= Mf_beta_lower)
		{
		  /* tapering is described in Eq. 45 of arXiv:2107.08876 */
		  REAL8 oscillations = beta_full - beta_SingleSpin;
		  REAL8 envelope = cos(2.0 * LAL_PI * Mf / (4.0 * Mf_beta_lower)) * cos(2.0 * LAL_PI * Mf / (4.0 * Mf_beta_lower));
		  beta = (beta_SingleSpin + oscillations * envelope);
		}
	      else
		{ /* otherwise just return single-spin beta */
		  beta = beta_SingleSpin;
		}
	    }
	  else
	    { /* no oscillations, just return full beta */
	      beta = beta_full;
	    }

	  break;
	}
      default:
	{
	  XLAL_ERROR(XLAL_EINVAL, "Error: IMRPhenomXPrecessionVersion not recognized in IMRPhenomX_PNR_GetPNBetaAtFreq.\n");
	  break;
	}
      }

     return beta;
   }

   /**
    * A wrapper to produce either the NNLO or MSA beta depending on the IMRPhenomXPrecVersion.
    *
    * This version does not modify the two-spin oscillations.
    */

   REAL8 IMRPhenomX_PNR_GetPNBetaAtFreq_fulltwospin(
       REAL8 Mf,                          /**< geometric frequency */
       IMRPhenomXWaveformStruct *pWF,     /**< PhenomX waveform struct */
       IMRPhenomXPrecessionStruct *pPrec, /**< PhenomX precession struct */
       const IMRPhenomX_PNR_beta_parameters *betaParams /**< beta parameter struct */
						    )
   {
     REAL8 beta;

     /* get PN expansion parameter v = (pi M f)^(1/3) */
     const double omega = LAL_PI * Mf;
     const double omega_cbrt = cbrt(omega);

     switch (pPrec->IMRPhenomXPrecVersion)
       {
	 /* ~~~~~ Use NNLO PN Euler Angles - Appendix G of arXiv:2004.06503 and https://dcc.ligo.org/LIGO-T1500602 ~~~~~ */
       case 101:
       case 102:
       case 103:
       case 104:
	 {

	   const REAL8 L = XLALSimIMRPhenomXLPNAnsatz(omega_cbrt, pWF->eta / omega_cbrt, pPrec->L0, pPrec->L1, pPrec->L2, pPrec->L3, pPrec->L4, pPrec->L5, pPrec->L6, pPrec->L7, pPrec->L8, pPrec->L8L);

	   /*  Comment from IMRPhenomX_precession.c:
	       We ignore the sign of L + SL below:
	       s := Sp / (L + SL)
	   */
	   REAL8 s = pPrec->Sperp / (L + pPrec->SL);
	   REAL8 s2 = s * s;
	   beta = acos(copysign(1.0, L + pPrec->SL) / sqrt(1.0 + s2));

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

	   beta = acos(vangles.z);

	   break;
	 }
       case 330:
	 {

	   beta= beta_SpinTaylor_IMR(Mf,pWF,pPrec,betaParams);
     if(isnan(beta) || isinf(beta)) XLAL_ERROR(XLAL_EINVAL, "Error in %s: beta_SpinTaylor_IMR returned invalid value.\n",__func__);

	   break;

	 }
       default:
	 {
	   XLAL_ERROR(XLAL_EINVAL, "Error: IMRPhenomXPrecessionVersion not recognized in IMRPhenomX_PNR_GetPNBetaAtFreq_fulltwospin.\n");
	   break;
	 }
       }

     return beta;
   }

  /**
   * A wrapper to generate the "waveform" PN beta from Eq. 41 of arXiv:2107.08876.
   *
   * The wrapper goes through the trouble of computing the frequency-dependent Sperp
   * given by Eq. 47 of arXiv:2107.08876.
   */

  REAL8 IMRPhenomX_PNR_PNWaveformBetaWrapper(
      REAL8 Mf,                         /**< geometric frequency */
      REAL8 pn_beta,                    /**< MSA or NNLO beta */
      IMRPhenomXWaveformStruct *pWF,    /**< PhenomX waveform struct */
      IMRPhenomXPrecessionStruct *pPrec /**< PhenomX precession struct */
  )
  {
    /* grab needed parameters */
    REAL8 M = pWF->Mtot;
    REAL8 m1 = pWF->m1 * M;
    REAL8 m2 = pWF->m2 * M;
    REAL8 J0 = pPrec->J0;
    REAL8 L0 = pPrec->LRef;
    REAL8 chi_eff = pWF->chiEff;

    /* compute derived terms */
    REAL8 chi_parr = (m1 + m2) * chi_eff / m1;

    /* get orbital angular momentum */
    REAL8 v = cbrt(LAL_PI * Mf);
    REAL8 L_norm = pWF->eta / v;
    REAL8 L_3PN = M*M*XLALSimIMRPhenomXLPNAnsatz(v, L_norm, pPrec->L0, pPrec->L1, pPrec->L2, pPrec->L3, pPrec->L4, pPrec->L5, pPrec->L6, pPrec->L7, pPrec->L8, pPrec->L8L);

    /* compute spin and costheta in Eqs. 18-19 of arXiv:2107.08876 at given Mf */
    REAL8 chi_temp = IMRPhenomX_PNR_chi_calc(m1, L_3PN, J0, L0, chi_parr, pn_beta);
    REAL8 costheta_temp = chi_parr / chi_temp;

    /* evaluate Eq. 41 of arXiv:2107.08876 */
    return IMRPhenomX_PNR_PNWaveformBeta(Mf, pn_beta, m1, m2, chi_temp, costheta_temp);
  }

  /**
   * The magnitude of the effective total spin is computed from
   * the total and orbital angular momenta, J0 and L0 resp., along with the opening
   * angle, beta, between them.
   *
   * This procedure is outlined in Eqs. 47 and 18 of arXiv:2107.08876.
   */

  REAL8 IMRPhenomX_PNR_chi_calc(
      REAL8 m1,       /**< mass of primary (Msun) */
      REAL8 L,        /**< magnitude of L and Mf */
      REAL8 J0,       /**< initial magnitude of J at Mf_ref */
      REAL8 L0,       /**< initial magnitude of L at Mf_ref */
      REAL8 chi_parr, /**< combined spin parallel to L0 */
      REAL8 beta      /**< PN opening angle, either MSA or NNLO */
  )
  {
    /* compute frequency-dependent Sperp and scale by mass */
    REAL8 S_perp = (J0 - (L0 - L)) * sin(beta);
    REAL8 chi_p = S_perp / (m1 * m1);

    /* find magnitude */
    REAL8 chi = sqrt(chi_parr * chi_parr + chi_p * chi_p);

    return chi;
  }

  /**
   * The "waveform" PN beta from Eq. 41 of arXiv:2107.08876.
   *
   * This function maps the dynamics iota to a version of beta
   * more closely resembling the angle associated with the optimal
   * emission direction described in Sec. 6B of arXiv:2107.08876.
   */

  REAL8 IMRPhenomX_PNR_PNWaveformBeta(
      REAL8 Mf,      /**< geometric frequency */
      REAL8 iota,    /**< dynamics precession cone opening angle */
      REAL8 m1,      /**< mass of primary (scaled to total mass 1) */
      REAL8 m2,      /**< mass of secondary (scaled to total mass 1) */
      REAL8 chi,     /**< effective single spin magnitude */
      REAL8 costheta /**< effective single spin polar angle (rad) */
  )
  {

    /* calculate combinations of masses */
    REAL8 M = m1 + m2;
    REAL8 nu = (m1 * m2) / (M * M);
    REAL8 delta = (m1 - m2) / M;

    /* get polar angle */
    REAL8 theta = acos(costheta);

    /* calculate invariant velocity */
    REAL8 w_orb = LAL_PI * Mf;
    REAL8 v = pow(w_orb, 1.0 / 3.0);
    REAL8 v2 = v * v;
    REAL8 v3 = v * v2;

    /* calculate needed trig functions */
    REAL8 cos_iota = cos(iota);
    REAL8 sin_iota = sin(iota);

    REAL8 cos_half_iota = cos(iota / 2.0);
    REAL8 sin_half_iota = sin(iota / 2.0);

    REAL8 cos_theta = cos(theta);
    REAL8 sin_theta = sin(theta);

    /* compute numerator expansion of h21/h22: N0 + N2 * v**2 + N3 * v**3*/
    REAL8 N0 = 84.0 * sin_iota;
    REAL8 N2 = 2.0 * (55.0 * nu - 107.0) * sin_iota;
    REAL8 N3 = -7.0 * (5.0 * nu + 6.0 * delta + 6.0) * chi * (2.0 * cos_iota - 1.0) * sin_theta + 56.0 * (3.0 * LAL_PI - (1.0 + delta - nu) * chi * cos_theta) * sin_iota;

    REAL8 N = (N0 + N2 * v2 + N3 * v3) / cos_half_iota;

    /* compute denominator expansion of h21/h22: D0 + D2 * v**2 + D3 * v**3*/
    REAL8 D0 = 84.0 * cos_half_iota;
    REAL8 D2 = 2.0 * (55.0 * nu - 107.0) * cos_half_iota;
    REAL8 D3 = 14.0 * 4.0 * (3.0 * LAL_PI + (nu - 1.0 - delta) * chi * cos_theta) * cos_half_iota + 14.0 * (6.0 + 6.0 * delta + 5.0 * nu) * chi * sin_theta * sin_half_iota;

    REAL8 D = D0 + D2 * v2 + D3 * v3;

    /* evaluate Eq. 41 of arXiv:2107.08876. NOTE: there's a typo in arXiv:2107.08876
     * whereby the factor of 2 in the denominator of Eq. 39 was dropped. */
    return 2.0 * XLALSimIMRPhenomXatan2tol(N, 2.0 * D, MAX_TOL_ATAN);
  }

  /**
   * This function evaluates the Ansatz coefficients of beta outlined in
   * Eq. 49 of arXiv:2107.08876.
   *
   * See the discussion in Sec. 8D of arXiv:2107.08876 for an explanation of
   * the condition on B4.
   *
   * The definition of B0 has since changed and now depends on the angle of
   * the final spin; see discussion in the technical document. FIXME: add citation
   */

  int IMRPhenomX_PNR_precompute_beta_coefficients(
      IMRPhenomX_PNR_beta_parameters *betaParams, /**< beta parameter struct */
      IMRPhenomXWaveformStruct *pWF,              /**< PhenomX waveform struct */
      IMRPhenomXPrecessionStruct *pPrec           /**< PhenomX precession struct */
  )
  {

    /* get effective single spin parameters */
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

    /* approximate orientation of final spin */
    REAL8 costhetaf = pPrec->costheta_final_singleSpin;

    /* ensure B4 is sufficiently large */
    REAL8 B4 = IMRPhenomX_PNR_beta_B4_coefficient(eta, chi, costheta);
    if (B4 <= 175.0)
    {
      B4 = 175.0;
    }

    betaParams->B0 = acos(costhetaf) - IMRPhenomX_PNR_beta_B0_coefficient(eta, chi, costheta);
    betaParams->B1 = IMRPhenomX_PNR_beta_B1_coefficient(eta, chi, costheta);
    betaParams->B2 = IMRPhenomX_PNR_beta_B2_coefficient(eta, chi, costheta);
    betaParams->B3 = betaParams->B2 * IMRPhenomX_PNR_beta_B3_coefficient(eta, chi, costheta);
    betaParams->B4 = B4;
    betaParams->B5 = IMRPhenomX_PNR_beta_B5_coefficient(eta, chi, costheta);

    return XLAL_SUCCESS;
  }

  /**
   * These three functions produce the inspiral rescaling
   * of beta described in Sec. 8B of arXiv:2107.08876.
   *
   * - IMRPhenomX_PNR_beta_rescaling_1 computes b1 in Eq. 54
   * - IMRPhenomX_PNR_beta_rescaling_2 computes b2 in Eq. 55
   * - IMRPhenomX_PNR_rescale_beta_expression combines the results
   *   to produce Eq. 53.
   */

  REAL8 IMRPhenomX_PNR_beta_rescaling_1(
      REAL8 Mf,     /**< geometric frequency */
      REAL8 beta1,  /**< PN beta evaluated at Mf */
      REAL8 beta2,  /**< MR beta evaluated at Mf */
      REAL8 dbeta1, /**< derivative of PN beta at Mf */
      REAL8 dbeta2  /**< derivative of MR beta at Mf */
  )
  {

    REAL8 D = beta1 * beta1 * Mf;
    REAL8 N = -2.0 * beta1 * (beta2 - beta1) + (beta1 * dbeta2 - beta2 * dbeta1) * Mf;

    return -N / D;
  }

  REAL8 IMRPhenomX_PNR_beta_rescaling_2(
      REAL8 Mf,     /**< geometric frequency */
      REAL8 beta1,  /**< PN beta evaluated at Mf */
      REAL8 beta2,  /**< MR beta evaluated at Mf */
      REAL8 dbeta1, /**< derivative of PN beta at Mf */
      REAL8 dbeta2  /**< derivative of MR beta at Mf */
  )
  {

    REAL8 D = (beta1 * Mf) * (beta1 * Mf);
    REAL8 N = beta1 * (beta2 - beta1) - (beta1 * dbeta2 - beta2 * dbeta1) * Mf;

    return -N / D;
  }

  REAL8 IMRPhenomX_PNR_rescale_beta_expression(
      REAL8 Mf,                                        /**< geometric frequency */
      const IMRPhenomX_PNR_beta_parameters *betaParams /**< beta parameter struct */
  )
  {

    REAL8 b1 = betaParams->beta_rescale_1;
    REAL8 b2 = betaParams->beta_rescale_2;

    return 1.0 + b1 * Mf + b2 * Mf * Mf;
  }

  /**
   * These four functions produce the MR Ansatz
   * of beta described in Sec. 7A of arXiv:2107.08876.
   *
   * - IMRPhenomX_PNR_MR_beta_expression computes the MR Ansatz
   *   detailed in Eq. 49
   * - IMRPhenomX_PNR_MR_dbeta_expression computes the first derivative of Eq. 49
   * - IMRPhenomX_PNR_MR_ddbeta_expression computes the second derivative of Eq. 49
   * - IMRPhenomX_PNR_MR_dddbeta_expression computes the third derivative of Eq. 49
   */
  REAL8 IMRPhenomX_PNR_MR_beta_expression(
      REAL8 Mf,                                        /**< geometric frequency */
      const IMRPhenomX_PNR_beta_parameters *betaParams /**< beta parameter struct */
  )
  {

    REAL8 B0 = betaParams->B0;
    REAL8 B1 = betaParams->B1;
    REAL8 B2 = betaParams->B2;
    REAL8 B3 = betaParams->B3;
    REAL8 B4 = betaParams->B4;
    REAL8 B5 = betaParams->B5;

    return B0 + (B1 + B2 * Mf + B3 * Mf * Mf) / (1.0 + B4 * (Mf + B5) * (Mf + B5));
  }

  /**
   *  expression for first derivative of beta in merger-ringdown regime
   */
  REAL8 IMRPhenomX_PNR_MR_dbeta_expression(
      REAL8 Mf,                                        /**< geometric frequency */
      const IMRPhenomX_PNR_beta_parameters *betaParams /**< beta parameter struct */
  )
  {

    REAL8 B1 = betaParams->B1;
    REAL8 B2 = betaParams->B2;
    REAL8 B3 = betaParams->B3;
    REAL8 B4 = betaParams->B4;
    REAL8 B5 = betaParams->B5;

    return ((2.0 * B3 * B4 * B5 - B2 * B4) * Mf * Mf + (2.0 * B3 - 2.0 * B1 * B4 + 2.0 * B3 * B4 * B5 * B5) * Mf + (B2 - 2.0 * B1 * B4 * B5 + B2 * B4 * B5 * B5)) / pow((1.0 + B4 * (Mf + B5) * (Mf + B5)), 2);
  }

  /**
   *  expression for second derivative of beta in merger-ringdown regime
   */
  REAL8 IMRPhenomX_PNR_MR_ddbeta_expression(
      REAL8 Mf,                                        /**< geometric frequency */
      const IMRPhenomX_PNR_beta_parameters *betaParams /**< beta parameter struct */
  )
  {

    REAL8 B1 = betaParams->B1;
    REAL8 B2 = betaParams->B2;
    REAL8 B3 = betaParams->B3;
    REAL8 B4 = betaParams->B4;
    REAL8 B5 = betaParams->B5;

    REAL8 a = B2 * B4 * B4 - 2.0 * B3 * B4 * B4 * B5;
    REAL8 b = -3.0 * B3 * B4 + 3.0 * B1 * B4 * B4 - 3.0 * B3 * B4 * B4 * B5 * B5;
    REAL8 c = -3.0 * B2 * B4 + 6.0 * B1 * B4 * B4 * B5 - 3.0 * B2 * B4 * B4 * B5 * B5;
    REAL8 d = B3 - B1 * B4 - 2.0 * B2 * B4 * B5 + 2.0 * B3 * B4 * B5 * B5 + 3.0 * B1 * B4 * B4 * B5 * B5 - 2.0 * B2 * B4 * B4 * B5 * B5 * B5 + B3 * B4 * B4 * B5 * B5 * B5 * B5;

    return 2.0 * (a * Mf * Mf * Mf + b * Mf * Mf + c * Mf + d) / pow((1.0 + B4 * (B5 + Mf) * (B5 + Mf)), 3);
  }

  /**
   *  expression for third derivative of beta in merger-ringdown regime
   */
  REAL8 IMRPhenomX_PNR_MR_dddbeta_expression(
      REAL8 Mf,                                        /**< geometric frequency */
      const IMRPhenomX_PNR_beta_parameters *betaParams /**< beta parameter struct */
  )
  {

    REAL8 B1 = betaParams->B1;
    REAL8 B2 = betaParams->B2;
    REAL8 B3 = betaParams->B3;
    REAL8 B4 = betaParams->B4;
    REAL8 B5 = betaParams->B5;

    REAL8 a = -B2 * B4 * B4 + 2.0 * B3 * B4 * B4 * B5;
    REAL8 b = 4.0 * B3 * B4 - 4.0 * B1 * B4 * B4 + 4.0 * B3 * B4 * B4 * B5 * B5;
    REAL8 c = 6.0 * B2 * B4 - 12.0 * B1 * B4 * B4 * B5 + 6.0 * B2 * B4 * B4 * B5 * B5;
    REAL8 d = -4.0 * B3 + 4.0 * B1 * B4 + 8.0 * B2 * B4 * B5 - 8.0 * B3 * B4 * B5 * B5 - 12.0 * B1 * B4 * B4 * B5 * B5 + 8.0 * B2 * B4 * B4 * B5 * B5 * B5 - 4.0 * B3 * B4 * B4 * B5 * B5 * B5 * B5;
    REAL8 e = -B2 - 2.0 * B3 * B5 + 4.0 * B1 * B4 * B5 + 2.0 * B2 * B4 * B5 * B5 - 4.0 * B3 * B4 * B5 * B5 * B5 - 4.0 * B1 * B4 * B4 * B5 * B5 * B5 + 3.0 * B2 * B4 * B4 * B5 * B5 * B5 * B5 - 2.0 * B3 * B4 * B4 * B5 * B5 * B5 * B5 * B5;

    return 6.0 * B4 * (a * Mf * Mf * Mf * Mf + b * Mf * Mf * Mf + c * Mf * Mf + d * Mf + e) / pow((1.0 + B4 * (B5 + Mf) * (B5 + Mf)), 4);
  }

  /**
   * Here we work through the construction of the connection frequency
   * for beta, outlined in Sec. 8B of arXiv:2107.08876, along with
   * discussion in Sec. 8D.
   *
   * In particular, this function performs the following tasks in order:
   * - compute the inflection frequency required to get beta_inf
   * - get derivative of beta at that frequency
   * - get the extremal frequencies of beta
   * - choose the lower and upper connection frequencies
   *
   */
  int IMRPhenomX_PNR_BetaConnectionFrequencies(
      IMRPhenomX_PNR_beta_parameters *betaParams /**< beta parameter struct */
  )
  {
    /* check for initialization */
    XLAL_CHECK(betaParams != NULL, XLAL_EFAULT);

    REAL8 B1 = betaParams->B1;
    REAL8 B2 = betaParams->B2;
    REAL8 B3 = betaParams->B3;
    REAL8 B4 = betaParams->B4;
    REAL8 B5 = betaParams->B5;

    REAL8 Mf_inflection;
    REAL8 root_term;
    REAL8 Mf_plus;
    REAL8 Mf_minus;
    REAL8 ddbeta_Mf_plus;
    REAL8 Mf_at_minimum;
    REAL8 Mf_at_maximum;
    REAL8 Mf_low;
    REAL8 Mf_IM;
    REAL8 Mf_MR;
    REAL8 ddbeta;
    REAL8 dddbeta;
    REAL8 dbeta_inflection;
    REAL8 chosen_dbeta;
    REAL8 d1;
    REAL8 d2;
    REAL8 delta;

    /* calculate inflection frequency to get beta_inf */
    Mf_inflection = IMRPhenomX_PNR_single_inflection_point(betaParams);

    /* calculate the derivative of beta at the inflection, dbeta_inf*/
    dbeta_inflection = IMRPhenomX_PNR_MR_dbeta_expression(Mf_inflection, betaParams);

    /* next we seek to compute Mf_max used in Eq. 57 of arXiv:2107.08876 */
    /* start by calculating the two extrema of beta FIXME: add documentation */
    root_term = B4 * (B2 - 2.0 * B3 * B5) * (B2 - 2.0 * B1 * B4 * B5 + B2 * B4 * B5 * B5) + (B3 - B1 * B4 + B3 * B4 * B5 * B5) * (B3 - B1 * B4 + B3 * B4 * B5 * B5);

    Mf_plus = ((B3 - B1 * B4 + B3 * B4 * B5 * B5) + sqrt(root_term)) / (B4 * (B2 - 2.0 * B3 * B5));
    Mf_minus = ((B3 - B1 * B4 + B3 * B4 * B5 * B5) - sqrt(root_term)) / (B4 * (B2 - 2.0 * B3 * B5));

    /* classify extrema as either maximum or minimum depending on sign of second derivative */
    ddbeta_Mf_plus = IMRPhenomX_PNR_MR_ddbeta_expression(Mf_plus, betaParams);

    if (ddbeta_Mf_plus > 0.0)
    {
      /* smiley face, Mf_plus is the minimum */
      Mf_at_minimum = Mf_plus;
      Mf_at_maximum = Mf_minus;
    }
    else
    {
      /* frowny face, Mf_plus is the maximum */
      Mf_at_minimum = Mf_minus;
      Mf_at_maximum = Mf_plus;
    }

    /* calculate the second and third derivatives at the maximum frequency */
    ddbeta = IMRPhenomX_PNR_MR_ddbeta_expression(Mf_at_maximum, betaParams);
    dddbeta = IMRPhenomX_PNR_MR_dddbeta_expression(Mf_at_maximum, betaParams);

    /* ensure that the sign of dbeta_inflection is maintained in Eq. 56 of arXiv:2107.08876,
     * even though dbeta_inflection is squared */
    REAL8 sign = 0.0;
    if (dbeta_inflection > 0.0)
    {
      sign = +1.0;
    }
    else
    {
      sign = -1.0;
    }
    /* compute Eq. 56 of arXiv:2107.08876 */
    chosen_dbeta = sign * ((dbeta_inflection / 100.0) * (dbeta_inflection / 100.0)) * 25.0;

    /* compute the two possible square roots for Eq. 57 of arXiv:2107.08876 */
    d1 = 1.0 / dddbeta * (-ddbeta + sqrt(ddbeta * ddbeta + 2.0 * dddbeta * chosen_dbeta));
    d2 = 1.0 / dddbeta * (-ddbeta - sqrt(ddbeta * ddbeta + 2.0 * dddbeta * chosen_dbeta));

    /* select the most appropriate case. The logic is as follows:
     * - if both d1 and d2 are positive, select the smallest
     * - if d2 < 0 and d1 > 0, choose d1
     * - if d1 < 0, choose d2 */
    if (d1 > 0.0)
    {
      if (d2 > 0.0)
      {
        /* choose the smallest positive solution */
        delta = (d1 < d2) ? d1 : d2;
      }
      else
      {
        delta = d1;
      }
    }
    else
    {
      delta = d2;
    }

    /* specify the connection frequency based on Eq. 58 of arXiv:2107.08876
     * for cases where the turnover is not present within the fitting region.
     * We will decide whether to use this or not below. */
    if (Mf_inflection >= 0.06)
    {
      Mf_low = Mf_inflection - 0.03;
    }
    else
    {
      Mf_low = 3.0 * Mf_inflection / 5.0;
    }

    /* quantify the shape of beta for this case: see Fig. 10 of arXiv:2107.08876 for
     * visualization of the three different cases.
     *
     * We start by checking to see if the minimum happens at a higher frequency to the maximum
     * and if the chosen inflection point is also at higher frequency than the maximum. In this case,
     * we are in the first two panels of Fig. 10. */
    if ((Mf_at_minimum > Mf_at_maximum) || (Mf_inflection > Mf_at_maximum))
    {
      /* If the maximum occurs at a higher frequency than in Eq. 58, pick Eq. 57 as the connection frequency.
       * Otherwise use Eq. 58. */
      if (Mf_at_maximum >= Mf_low)
      {
        Mf_IM = Mf_at_maximum + delta;
      }
      else
      {
        Mf_IM = Mf_low;
      }
    }
    else /* Here we are in the right-most panel of Fig. 10 */
    {
      /* ensure that we have enough room to start below the minimum frequency */
      if (Mf_at_minimum > 0.06)
      {
        Mf_IM = Mf_at_minimum - 0.03;
      }
      else
      {
        Mf_IM = 3.0 * Mf_at_minimum / 5.0;
      }
    }

    /* Specify the upper connection frequency where we either transition to a constant
     * value for beta at the minimum to enforce zero slope, or we are in the first panel
     * of Fig. 10 where beta will asymptote to some value.
     * In the latter case, we simply set the upper connection frequency to be very large */
    if (Mf_at_minimum > Mf_inflection)
    {
      Mf_MR = Mf_at_minimum;
    }
    else
    {
      Mf_MR = 100.;
    }

    /* failsafe if the lower connection frequency is negative, then we won't be using the MR fit */
    if ((Mf_IM < 0.0) || isnan(Mf_IM))
    {
      betaParams->Mf_beta_lower = 100.0;
      betaParams->Mf_beta_upper = 100.0;
    }
    else
    {
      betaParams->Mf_beta_lower = Mf_IM;
      betaParams->Mf_beta_upper = Mf_MR;
    }

    return XLAL_SUCCESS;
  }

  /**
   * Compute the three roots of a depressed cubic given by
   * Eq. 68 of arXiv:2107.08876.
   *
   */
  COMPLEX16 *IMRPhenomX_PNR_three_inflection_points(
      const IMRPhenomX_PNR_beta_parameters *betaParams /**< beta parameter struct */
  )
  {

    REAL8 B1 = betaParams->B1;
    REAL8 B2 = betaParams->B2;
    REAL8 B3 = betaParams->B3;
    REAL8 B4 = betaParams->B4;
    REAL8 B5 = betaParams->B5;

    COMPLEX16 a;
    COMPLEX16 b;
    COMPLEX16 c;
    COMPLEX16 d;
    COMPLEX16 r;
    COMPLEX16 s;
    COMPLEX16 phi;
    static COMPLEX16 f[3];

    /* cubic ax^3 + bx^2 + cx + d FIXME: add documentation */
    a = 2 * (B2 * B4 * B4 - 2 * B3 * B4 * B4 * B5);
    b = 6 * (-B3 * B4 + B1 * B4 * B4 - B3 * B4 * B4 * B5 * B5);
    c = 6 * (-B2 * B4 + 2 * B1 * B4 * B4 * B5 - B2 * B4 * B4 * B5 * B5);
    d = 2 * (B3 - B1 * B4 - 2 * B2 * B4 * B5 + 2 * B3 * B4 * B5 * B5 + 3 * B1 * B4 * B4 * B5 * B5 - 2 * B2 * B4 * B4 * B5 * B5 * B5 + B3 * B4 * B4 * B5 * B5 * B5 * B5);

    /* convert to depressed cubic t^3 + rt +s, t = x + b/3a */
    r = (3 * a * c - b * b) / (3 * a * a);
    s = (2 * b * b * b - 9 * a * b * c + 27 * a * a * d) / (27 * a * a * a);

    phi = acos(((3 * s) / (2 * r)) * csqrt(-3 / r));

    for (int n = 0; n < 3; n++)
    {
      f[n] = 2 * csqrt(-r / 3) * cos((phi - 2 * n * LAL_PI) / 3) - b / (3 * a);
    }

    return f;
  }

  /**
   * Compute the two roots of a depressed cubic given by
   * Eq. 68 of arXiv:2107.08876 when the cubic term vanishes.
   *
   */
  COMPLEX16 *IMRPhenomX_PNR_two_inflection_points(
      const IMRPhenomX_PNR_beta_parameters *betaParams /**< beta parameter struct */
  )
  {

    REAL8 B1 = betaParams->B1;
    REAL8 B2 = betaParams->B2;
    REAL8 B3 = betaParams->B3;
    REAL8 B4 = betaParams->B4;
    REAL8 B5 = betaParams->B5;

    REAL8 b;
    REAL8 c;
    REAL8 d;

    static COMPLEX16 f[2];

    /* cubic ax^3 + bx^2 + cx + d, here a=0 */
    b = 6 * (-B3 * B4 + B1 * B4 * B4 - B3 * B4 * B4 * B5 * B5);
    c = 6 * (-B2 * B4 + 2 * B1 * B4 * B4 * B5 - B2 * B4 * B4 * B5 * B5);
    d = 2 * (B3 - B1 * B4 - 2 * B2 * B4 * B5 + 2 * B3 * B4 * B5 * B5 + 3 * B1 * B4 * B4 * B5 * B5 - 2 * B2 * B4 * B4 * B5 * B5 * B5 + B3 * B4 * B4 * B5 * B5 * B5 * B5);

    f[0] = (1 / (2 * b)) * (-c + csqrt(c * c - 4 * b * d));
    f[1] = (1 / (2 * b)) * (-c - csqrt(c * c - 4 * b * d));

    return f;
  }

  /**
   * Compute the roots of a depressed cubic given by
   * Eq. 68 of arXiv:2107.08876, then select the correct
   * root to use for the inflection frequency based on
   * Eq. 69 and discussion in Sec. 8D.
   *
   */
  REAL8 IMRPhenomX_PNR_single_inflection_point(
      const IMRPhenomX_PNR_beta_parameters *betaParams /**< beta parameter struct */
  )
  {

    REAL8 B1 = betaParams->B1;
    REAL8 B2 = betaParams->B2;
    REAL8 B3 = betaParams->B3;
    REAL8 B4 = betaParams->B4;
    REAL8 B5 = betaParams->B5;

    // cubic terms
    REAL8 a;
    REAL8 b;
    REAL8 c;
    REAL8 d;

    int i;

    REAL8 grad;
    REAL8 f_inflection = 0.0;

    /* cubic ax^3 + bx^2 + cx + d */
    a = 2.0 * (B2 * B4 * B4 - 2.0 * B3 * B4 * B4 * B5);
    b = 6.0 * (-B3 * B4 + B1 * B4 * B4 - B3 * B4 * B4 * B5 * B5);
    c = 6.0 * (-B2 * B4 + 2.0 * B1 * B4 * B4 * B5 - B2 * B4 * B4 * B5 * B5);
    d = 2.0 * (B3 - B1 * B4 - 2.0 * B2 * B4 * B5 + 2.0 * B3 * B4 * B5 * B5 + 3.0 * B1 * B4 * B4 * B5 * B5 - 2.0 * B2 * B4 * B4 * B5 * B5 * B5 + B3 * B4 * B4 * B5 * B5 * B5 * B5);

    /* treat cases where co-efficients of cubic term is zero (or less than 10^-10 to account for numerical error) */
    if (fabs(a) < 1e-10)
    {
      /* co-efficient of quadratic term is also zero */
      if (fabs(b) < 2e-10)
      {
        /* calculate single inflection point */
        f_inflection = -d / c;
      }
      else
      {
        COMPLEX16 *f_inf;

        /* calculate 2 inflection points */
        f_inf = IMRPhenomX_PNR_two_inflection_points(betaParams);

        /* choose inflection point with negative gradient */
        for (i = 0; i < 2; i++)
        {
          grad = IMRPhenomX_PNR_MR_dbeta_expression(f_inf[i], betaParams);
          if (grad < 0.0)
          {
            f_inflection = f_inf[i];
          }
        }
      }
    }
    /* treat cases with 3 roots */
    else
    {
      COMPLEX16 *f_inf;
      REAL8 f_temp = 0;
      REAL8 f_IM = 0;
      int w = 0; // initialise counter

      /* calculate all 3 inflection points */
      f_inf = IMRPhenomX_PNR_three_inflection_points(betaParams);

      /* check for substantial imaginary component and select real root if this component exists */
      for (i = 0; i < 3; i++)
      {
        f_IM = cimag(f_inf[i]);
        if (f_IM < 1e-10)
        {
          w += 1;
          f_temp = creal(f_inf[i]);
        }
      }
      /* case where we have just one real root */
      if (w == 1)
      {
        f_inflection = f_temp;
      }
      /* case where we have all 3 real roots */
      else
      {
        if (a < 0.0)
        {
          /* take the central root */
          f_inflection = creal(f_inf[1]);
        }
        else
        {
          /* enforce Eq. 69 */
          if (b / (3.0 * a) > B5 / 2.0 - (2141.0 / 90988.0)) // FIXME: add documentation
          {
            f_inflection = creal(f_inf[0]);
          }
          else
          {
            f_inflection = creal(f_inf[2]);
          }
        }
      }
    }

    return f_inflection;
  }

  /**
   * This function combines several functions together to
   * compute the connection frequencies and the inspiral rescaling.
   *
   */
  int IMRPhenomX_PNR_beta_connection_parameters(
      IMRPhenomX_PNR_beta_parameters *betaParams,  /**< beta parameter struct */
      IMRPhenomXWaveformStruct *pWF,               /**< PhenomX waveform struct */
      IMRPhenomXPrecessionStruct *pPrec,           /**< PhenomX precession struct */
      IMRPhenomXWaveformStruct *pWF_SingleSpin,    /**< PhenomX waveform struct with single-spin mapping */
      IMRPhenomXPrecessionStruct *pPrec_SingleSpin /**< PhenomX precession struct with single-spin mapping */
  )
  {
    /* may we continue? Single-spin structs might be NULL */
    XLAL_CHECK(betaParams != NULL, XLAL_EFAULT);
    XLAL_CHECK(pWF != NULL, XLAL_EFAULT);
    XLAL_CHECK(pPrec != NULL, XLAL_EFAULT);

    /* define frequency spacing over which to calculate derivatives */
    /* here set to be 0.0005Mf */
    REAL8 dMf = 0.0005;

    /* generate beta connection frequencies*/
    IMRPhenomX_PNR_BetaConnectionFrequencies(betaParams);

    if((pPrec->IMRPhenomXPrecVersion==330)&&(betaParams->Mf_beta_lower-2.*dMf<pPrec->Mfmin_integration))
      betaParams->Mf_beta_lower = 100.;

    /* evaluate expressions for beta at lower connection frequency */
    REAL8 Mf_beta_lower = betaParams->Mf_beta_lower;

    /* get the values of PN and MR beta at Mf_beta_lower, as well as around that
     * frequency for the centered FD derivatives */
    REAL8 beta_temp = IMRPhenomX_PNR_GetPNBetaAtFreq(Mf_beta_lower - dMf, betaParams, pWF, pPrec, pWF_SingleSpin, pPrec_SingleSpin);
    REAL8 b1 = IMRPhenomX_PNR_PNWaveformBetaWrapper(Mf_beta_lower - dMf, beta_temp, pWF, pPrec);

    beta_temp = IMRPhenomX_PNR_GetPNBetaAtFreq(Mf_beta_lower, betaParams, pWF, pPrec, pWF_SingleSpin, pPrec_SingleSpin);
    REAL8 beta_lower = IMRPhenomX_PNR_PNWaveformBetaWrapper(Mf_beta_lower, beta_temp, pWF, pPrec);

    beta_temp = IMRPhenomX_PNR_GetPNBetaAtFreq(Mf_beta_lower + dMf, betaParams, pWF, pPrec, pWF_SingleSpin, pPrec_SingleSpin);
    REAL8 b3 = IMRPhenomX_PNR_PNWaveformBetaWrapper(Mf_beta_lower + dMf, beta_temp, pWF, pPrec);

    REAL8 b4 = IMRPhenomX_PNR_MR_beta_expression(Mf_beta_lower - dMf, betaParams);
    REAL8 beta_upper = IMRPhenomX_PNR_MR_beta_expression(Mf_beta_lower, betaParams);
    REAL8 b6 = IMRPhenomX_PNR_MR_beta_expression(Mf_beta_lower + dMf, betaParams);

    /* calculate derivative of alpha at connection frequencies using central finite difference method */
    REAL8 derivative_beta_lower = (b3 - b1) / (2.0 * dMf);
    REAL8 derivative_beta_upper = (b6 - b4) / (2.0 * dMf);

    /* compute the inspiral rescaling */
    REAL8 beta_rescale_1 = IMRPhenomX_PNR_beta_rescaling_1(Mf_beta_lower, beta_lower, beta_upper, derivative_beta_lower, derivative_beta_upper);
    REAL8 beta_rescale_2 = IMRPhenomX_PNR_beta_rescaling_2(Mf_beta_lower, beta_lower, beta_upper, derivative_beta_lower, derivative_beta_upper);

    beta_rescale_1 = isnan(beta_rescale_1) ? 0.0 : beta_rescale_1;
    beta_rescale_2 = isnan(beta_rescale_2) ? 0.0 : beta_rescale_2;

    /* save beta values to the parameter dict */
    betaParams->beta_lower = beta_lower;
    betaParams->beta_upper = beta_upper;
    betaParams->derivative_beta_lower = derivative_beta_lower;
    betaParams->derivative_beta_upper = derivative_beta_upper;
    betaParams->beta_rescale_1 = beta_rescale_1;
    betaParams->beta_rescale_2 = beta_rescale_2;

    return XLAL_SUCCESS;
  }

  /**
   * Function to calculate beta and first derivative at connection frequency between SpinTaylor and PNR angles
   */

  int IMRPhenomX_ST_PNR_beta_connection_parameters(
      IMRPhenomX_PNR_beta_parameters *betaParams,  /**< beta parameter struct */
      IMRPhenomXWaveformStruct *pWF,               /**< PhenomX waveform struct */
      IMRPhenomXPrecessionStruct *pPrec,           /**< PhenomX precession struct */
      IMRPhenomXWaveformStruct *pWF_SingleSpin,       /**< PhenomX waveform struct with max-spin mapping */
      IMRPhenomXPrecessionStruct *pPrec_SingleSpin    /**< PhenomX precession struct with max-spin mapping */

  )
  {
    /* may we continue? Single-spin structs might be NULL */
    XLAL_CHECK(betaParams != NULL, XLAL_EFAULT);
    XLAL_CHECK(pPrec != NULL, XLAL_EFAULT);

    /* define frequency spacing over which to calculate derivatives */
    /* here set to be 0.0005Mf */
    REAL8 dMf = 0.0005;

    /* generate beta connection frequencies*/
    IMRPhenomX_PNR_BetaConnectionFrequencies(betaParams);

    /* frequencies marking beginning and end of interpolation region */
    REAL8 ftrans = ( pPrec->ftrans_MRD <= 0.9*betaParams->Mf_beta_lower ) ? pPrec->ftrans_MRD : 0.9*betaParams->Mf_beta_lower;
    REAL8 MfA = ( ftrans-dMf  < pPrec->Mfmin_integration + 2.*dMf ? pPrec->Mfmin_integration + 2*dMf : ftrans - dMf );
    REAL8 MfB = betaParams->Mf_beta_lower;

    /* evaulate SpinTaylor angles around lower connection frequency */
    REAL8 cosbeta=0.;
    int status=GSL_SUCCESS;

    status=gsl_spline_eval_e(pPrec->cosbeta_spline, MfA - dMf, pPrec->cosbeta_acc,&cosbeta);
    XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "%s: error in beta interpolation for q=%.12f, Mtot=%.12f and chi1=[%.12f, %.12f,%.12f], chi2=[%.12f, %.12f,%.12f]\n",__func__,pWF->q,pWF->M, pPrec->chi1x,pPrec->chi1y,pPrec->chi1z,pPrec->chi2x,pPrec->chi2y,pPrec->chi2z);
    REAL8 bA1 = acos(MAX(-1, MIN(1, cosbeta)));

    status=gsl_spline_eval_e(pPrec->cosbeta_spline, MfA, pPrec->cosbeta_acc,&cosbeta);
    XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "%s: error in beta interpolation for q=%.12f, Mtot=%.12f and chi1=[%.12f, %.12f,%.12f], chi2=[%.12f, %.12f,%.12f]\n",__func__,pWF->q,pWF->M, pPrec->chi1x,pPrec->chi1y,pPrec->chi1z,pPrec->chi2x,pPrec->chi2y,pPrec->chi2z);
    REAL8 bA = acos(MAX(-1, MIN(1, cosbeta)));

    status=gsl_spline_eval_e(pPrec->cosbeta_spline, MfA + dMf, pPrec->cosbeta_acc,&cosbeta);
    XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "%s: error in beta interpolation for q=%.12f, Mtot=%.12f and chi1=[%.12f, %.12f,%.12f], chi2=[%.12f, %.12f,%.12f]\n",__func__,pWF->q,pWF->M, pPrec->chi1x,pPrec->chi1y,pPrec->chi1z,pPrec->chi2x,pPrec->chi2y,pPrec->chi2z);
    REAL8 bA2 = acos(MAX(-1, MIN(1, cosbeta)));

    /* evaluate single-spin MSA angles at upper connection frequency */
    const double omega = LAL_PI * MfB;
    const double domega = LAL_PI * dMf;

    REAL8 bB1, bB, bB2;
    if (IMRPhenomX_PNR_CheckTwoSpin(pPrec))
      {
	const vector vangles_SingleSpin1 = IMRPhenomX_Return_phi_zeta_costhetaL_MSA(cbrt(omega - domega), pWF_SingleSpin, pPrec_SingleSpin);
	bB1 = acos(vangles_SingleSpin1.z);

	const vector vangles_SingleSpin2 = IMRPhenomX_Return_phi_zeta_costhetaL_MSA(cbrt(omega), pWF_SingleSpin, pPrec_SingleSpin);
	bB = acos(vangles_SingleSpin2.z);

	const vector vangles_SingleSpin3 = IMRPhenomX_Return_phi_zeta_costhetaL_MSA(cbrt(omega + domega), pWF_SingleSpin, pPrec_SingleSpin);
	bB2 = acos(vangles_SingleSpin3.z);
      }
    else
      {
	int user_version = pPrec->IMRPhenomXPrecVersion;
        pPrec->IMRPhenomXPrecVersion=223;
        IMRPhenomX_Initialize_MSA_System(pWF,pPrec,pPrec->ExpansionOrder);
        pPrec->IMRPhenomXPrecVersion=user_version;

        const vector vangles_SingleSpin1 = IMRPhenomX_Return_phi_zeta_costhetaL_MSA(cbrt(omega - domega), pWF, pPrec);
        bB1 = acos(vangles_SingleSpin1.z);

        const vector vangles_SingleSpin2 = IMRPhenomX_Return_phi_zeta_costhetaL_MSA(cbrt(omega), pWF, pPrec);
        bB = acos(vangles_SingleSpin2.z);

        const vector vangles_SingleSpin3 = IMRPhenomX_Return_phi_zeta_costhetaL_MSA(cbrt(omega + domega), pWF, pPrec);
        bB2 = acos(vangles_SingleSpin3.z);
      }

    /* compute numerical gradient */
    REAL8 dbA = (bA2-bA1)/(2.0*dMf);
    REAL8 dbB = (bB2-bB1)/(2.0*dMf);

    /* compute and store connection co-efficients */
    REAL8 b0 = IMRPhenomX_ST_PNR_beta_coeffs_0( MfA, MfB, bA, bA, dbA, dbB );
    REAL8 b1 = IMRPhenomX_ST_PNR_beta_coeffs_1( MfA, MfB, bA, bA, dbA, dbB );
    REAL8 b2 = IMRPhenomX_ST_PNR_beta_coeffs_2( MfA, MfB, bA, bA, dbA, dbB );
    REAL8 b3 = IMRPhenomX_ST_PNR_beta_coeffs_3( MfA, MfB, bA, bA, dbA, dbB );

    betaParams->beta_interp_0 = b0;
    betaParams->beta_interp_1 =	b1;
    betaParams->beta_interp_2 =	b2;
    betaParams->beta_interp_3 =	b3;

    /* Now parameters for rescaling */

    REAL8 bBnew = IMRPhenomX_PNR_CheckTwoSpin(pPrec) ? bB : bA;
    // MSA value if two spin because already tapered to this point
    // SpinTaylor otherwise to avoid introducing kink in beta where ST and MSA apparoximants disagree

    REAL8 bC1 = IMRPhenomX_PNR_PNWaveformBetaWrapper(MfB - dMf, bB1, pWF, pPrec);
    REAL8 bC  = IMRPhenomX_PNR_PNWaveformBetaWrapper(MfB, bBnew, pWF, pPrec);
    REAL8 bC2 = IMRPhenomX_PNR_PNWaveformBetaWrapper(MfB + dMf, bB2, pWF, pPrec);

    REAL8 bD1 = IMRPhenomX_PNR_MR_beta_expression(MfB - dMf, betaParams);
    REAL8 bD  = IMRPhenomX_PNR_MR_beta_expression(MfB, betaParams);
    REAL8 bD2 = IMRPhenomX_PNR_MR_beta_expression(MfB + dMf, betaParams);

    /* calculate derivative of beta at connection frequencies using central finite difference method */
    REAL8 dbC = (bC2 - bC1) / (2.0 * dMf);
    REAL8 dbD = (bD2 - bD1) / (2.0 * dMf);

    /* compute the inspiral rescaling */
    REAL8 beta_rescale_1 = IMRPhenomX_PNR_beta_rescaling_1(MfB, bC, bD, dbC, dbD);
    REAL8 beta_rescale_2 = IMRPhenomX_PNR_beta_rescaling_2(MfB, bC, bD, dbC, dbD);

    beta_rescale_1 = isnan(beta_rescale_1) ? 0.0 : beta_rescale_1;
    beta_rescale_2 = isnan(beta_rescale_2) ? 0.0 : beta_rescale_2;

    /* save beta values to the parameter dict */
    betaParams->beta_lower = bC;
    betaParams->beta_upper = bD;
    betaParams->derivative_beta_lower = dbC;
    betaParams->derivative_beta_upper = dbD;
    betaParams->beta_rescale_1 = beta_rescale_1;
    betaParams->beta_rescale_2 = beta_rescale_2;

    //printf("%f, %f \n", beta_rescale_1, beta_rescale_2);
    //printf("%f, %f, %f, %f, %f, %f \n", MfA, MfB, bA, bB, dbA, dbB);
    //printf("%f, %f, %f, %f\n", betaParams->beta_interp_0, betaParams->beta_interp_1, betaParams->beta_interp_2, betaParams->beta_interp_3);
    //printf("-------\n");

    return XLAL_SUCCESS;
  }

  /**
   * Series of functions to calculation co-efficients in beta interpolation function
   * These expressions come from enforcing C(1) connectivity to the ST and PNR expressions
   * See discussion below Eq. (19) in https://dcc.ligo.org/LIGO-T2400368
   */

  REAL8 IMRPhenomX_ST_PNR_beta_coeffs_0(
      REAL8 Mf1,     /**< geometric frequency at lower connection point */
      REAL8 Mf2,     /**< geometric frequency at lower connection point */
      REAL8 beta1,  /**< PN beta evaluated at Mf1 */
      REAL8 beta2,  /**< MR beta evaluated at Mf2 */
      REAL8 dbeta1, /**< derivative of PN beta at Mf1 */
      REAL8 dbeta2  /**< derivative of MR beta at Mf2 */
  )
  {
    REAL8 Mf1_2 = Mf1*Mf1;
    REAL8 Mf1_3 = Mf1_2*Mf1;

    REAL8 Mf2_2 = Mf2*Mf2;
    REAL8 Mf2_3 = Mf2_2*Mf2;

    REAL8 N = -beta2*Mf1_3 + 3*beta2*Mf1_2*Mf2 + dbeta2*Mf1_3*Mf2 - 3*beta1*Mf1*Mf2_2 + dbeta1*Mf1_2*Mf2_2 - dbeta2*Mf1_2*Mf2_2 + beta1*Mf2_3 - dbeta1*Mf1*Mf2_3;
    REAL8 D = (Mf1-Mf2)*(Mf1-Mf2)*(Mf1-Mf2);

    return -N/D;
  }

  REAL8 IMRPhenomX_ST_PNR_beta_coeffs_1(
      REAL8 Mf1,     /**< geometric frequency at lower connection point */
      REAL8 Mf2,     /**< geometric frequency at lower connection point */
      REAL8 beta1,  /**< PN beta evaluated at Mf1 */
      REAL8 beta2,  /**< MR beta evaluated at Mf2 */
      REAL8 dbeta1, /**< derivative of PN beta at Mf1 */
      REAL8 dbeta2  /**< derivative of MR beta at Mf2 */
  )
  {
    REAL8 Mf1_2 = Mf1*Mf1;
    REAL8 Mf1_3 = Mf1_2*Mf1;

    REAL8 Mf2_2 = Mf2*Mf2;
    REAL8 Mf2_3 = Mf2_2*Mf2;

    REAL8 N = -dbeta2*Mf1_3 + 6*beta1*Mf1*Mf2 - 6*beta2*Mf1*Mf2 - 2*dbeta1*Mf1_2*Mf2 - dbeta2*Mf1_2*Mf2 + dbeta1*Mf1*Mf2_2 + 2*dbeta2*Mf1*Mf2_2 + dbeta1*Mf2_3;
    REAL8 D = (Mf1-Mf2)*(Mf1-Mf2)*(Mf1-Mf2);

    return -N/D;
    }

  REAL8 IMRPhenomX_ST_PNR_beta_coeffs_2(
      REAL8 Mf1,     /**< geometric frequency at lower connection point */
      REAL8 Mf2,     /**< geometric frequency at lower connection point */
      REAL8 beta1,  /**< PN beta evaluated at Mf1 */
      REAL8 beta2,  /**< MR beta evaluated at Mf2 */
      REAL8 dbeta1, /**< derivative of PN beta at Mf1 */
      REAL8 dbeta2  /**< derivative of MR beta at Mf2 */
    )
    {
      REAL8 Mf1_2 = Mf1*Mf1;
      REAL8 Mf2_2 = Mf2*Mf2;

      REAL8 N = -3*(beta1 - beta2)*(Mf1 + Mf2) + (dbeta1 + 2*dbeta2)*Mf1_2 + (dbeta1 - dbeta2)*Mf1*Mf2 - (2*dbeta1 + dbeta2)*Mf2_2;
      REAL8 D = (Mf1-Mf2)*(Mf1-Mf2)*(Mf1-Mf2);

      return -N/D;
    }

    REAL8 IMRPhenomX_ST_PNR_beta_coeffs_3(
        REAL8 Mf1,     /**< geometric frequency at lower connection point */
        REAL8 Mf2,     /**< geometric frequency at lower connection point */
        REAL8 beta1,  /**< PN beta evaluated at Mf1 */
	REAL8 beta2,  /**< MR beta evaluated at Mf2 */
        REAL8 dbeta1, /**< derivative of PN beta at Mf1 */
	REAL8 dbeta2  /**< derivative of MR beta at Mf2 */
    )
    {

      REAL8 N = 2*(beta1-beta2) - (dbeta1+dbeta2)*(Mf1-Mf2);
      REAL8 D = (Mf1-Mf2)*(Mf1-Mf2)*(Mf1-Mf2);

      return -N/D;
    }

  /**
   * Utility function to compute the arctan windowing described
   * in Eq. 62 of arXiv:2107.08876.
   *
   */
  REAL8 IMRPhenomX_PNR_arctan_window(
      REAL8 beta /**< beta angle (rad) */
  )
  {
    /* specify the width as described in Sec. 8C of arXiv:2107.08876*/
    REAL8 window_border = 0.01;

    /* if beta is close to zero or PI, window it */
    if ((beta <= window_border) || (beta >= LAL_PI - window_border))
    {
      /* precompute some things, then evaluate the arctan */
      REAL8 p = 0.002;
      REAL8 PI_by_2 = 1.570796326794897;
      REAL8 PI_by_2_1mp = 1.569378278348018;
      REAL8 PI_by_2_1oq = 7.308338225719002e97;
      REAL8 sign = (beta < PI_by_2) ? -1.0 : 1.0;

      return sign * PI_by_2_1mp * pow(atan2(pow(beta - PI_by_2, 1.0 / p), PI_by_2_1oq), p) + PI_by_2;
    }

    return beta;
  }

  /**
   * Determine whether to attach the MR contributions to beta.
   * - If the connection frequency is too low, < 0.009
   * - If the connection frequency is 100.0
   * - if the final Ansatz value of beta is negative at the
   *   lower connection frequency
   *
   */
  UINT4 IMRPhenomX_PNR_AttachMRBeta(
      const IMRPhenomX_PNR_beta_parameters *betaParams /**< beta parameter struct */
  )
  {
    /* check for initialization */
    XLAL_CHECK(betaParams != NULL, XLAL_EFAULT);

    if ((betaParams->Mf_beta_lower < 0.009) || (betaParams->Mf_beta_lower == 100.0) || (betaParams->beta_upper <= 0.0) || (betaParams->beta_upper > 5.0 * (betaParams->B0 + 0.1)))
    {
      return 0;
    }

    return 1;
  }

#ifdef __cplusplus
}
#endif
