/*
 * Copyright (C) 2018 Geraint Pratten
 *
 * 	This code builds on:
 * 		LALSimIMRPhenomD_internals.c
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


/**
 * \author Geraint Pratten
 *
 * Internal functions for IMRPhenomX, arXiv:2001.11412
 */

/* IMRPhenomX Header Files */
#include "LALSimIMRPhenomXUtilities.h"
#include "LALSimIMRPhenomX_internals.h"

#include "LALSimIMRPhenomX_inspiral.c"
#include "LALSimIMRPhenomX_intermediate.c"
#include "LALSimIMRPhenomX_ringdown.c"
#include "LALSimIMRPhenomX_qnm.c"
#include "LALSimIMRPhenomX_PNR_deviations.c"

/* LAL Header Files */
#include <lal/LALSimIMR.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>

/* GSL Header Files */
#include <gsl/gsl_linalg.h>

/* Link PN coefficients (needed for NRTidal spin-induced quadrupole terms) */
#include "LALSimInspiralPNCoefficients.c"

/* This struct is used to pre-cache useful powers of frequency, avoiding numerous expensive operations */
int IMRPhenomX_Initialize_Powers(IMRPhenomX_UsefulPowers *p, REAL8 number)
{
	XLAL_CHECK(0 != p, XLAL_EFAULT, "p is NULL");
	XLAL_CHECK(number >= 0, XLAL_EDOM, "number must be non-negative");

	double sixth      = pow(number, 1.0 / 6.0);
	double m_sixth    = 1.0 / sixth;

	p->one_sixth      = sixth;
	p->m_one_sixth    = m_sixth;

	p->m_one          = 1.0 / number;
	p->itself         = number;

	p->one_third      = sixth * sixth;
	p->m_one_third    = 1.0 / (p->one_third);

	p->two_thirds     = p->one_third * p->one_third;
	p->m_two_thirds   = 1.0 / p->two_thirds;

	p->four_thirds    = p->two_thirds * p->two_thirds;
	p->m_four_thirds  = 1.0 / p->four_thirds;

	p->five_thirds    = p->four_thirds * p->one_third;
	p->m_five_thirds  = 1.0 / p->five_thirds;

	p->seven_thirds   = p->four_thirds * number;
	p->m_seven_thirds = 1.0 / p->seven_thirds;

	p->eight_thirds   = p->seven_thirds * p->one_third;
	p->m_eight_thirds = 1.0 / p->eight_thirds;

    p->ten_thirds   = p->seven_thirds * number;
    p->m_ten_thirds = 1.0 / p->ten_thirds;

	p->two            = number   * number;
	p->three          = p->two   * number;
	p->four           = p->three * number;
	p->five           = p->four  * number;

	p->m_two          = 1.0 / p->two;
	p->m_three        = 1.0 / p->three;
	p->m_four         = 1.0 / p->four;
	p->m_five         = 1.0 / p->five;

	p->seven_sixths   = p->one_sixth   * p->itself;
	p->m_seven_sixths = p->m_one_sixth * p->m_one;

	p->log            = log(number);
	p->sqrt           = p->one_sixth*p->one_sixth*p->one_sixth;

	return XLAL_SUCCESS;
}

/* A stripped down version of IMRPhenomX_Initialize_Powers for main production loop */
int IMRPhenomX_Initialize_Powers_Light(IMRPhenomX_UsefulPowers *p, REAL8 number)
{
	XLAL_CHECK(0 != p, XLAL_EFAULT, "p is NULL");
	XLAL_CHECK(number >= 0, XLAL_EDOM, "number must be non-negative");

	double sixth      = pow(number, 1.0 / 6.0);
	double m_sixth    = 1.0 / sixth;

	p->one_sixth      = sixth;
	p->m_one_sixth    = m_sixth;

	p->m_one          = 1.0 / number;
	p->itself         = number;

	p->one_third      = sixth * sixth;
	p->two_thirds     = p->one_third * p->one_third;

	p->m_two          = p->m_one   * p->m_one;
	p->m_three        = p->m_two   * p->m_one;

	p->seven_sixths   = p->one_sixth   * p->itself;
	p->m_seven_sixths = p->m_one_sixth * p->m_one;

	p->log            = log(number);

	return XLAL_SUCCESS;
}


int IMRPhenomXSetWaveformVariables(
	IMRPhenomXWaveformStruct *wf,
  const REAL8 m1_SI,
  const REAL8 m2_SI,
  const REAL8 chi1L_In,
  const REAL8 chi2L_In,
  const REAL8 deltaF,
  const REAL8 fRef,
  const REAL8 phi0,
  const REAL8 f_min,
  const REAL8 f_max,
  const REAL8 distance,
  const REAL8 inclination,
  LALDict *LALParams,
	UNUSED const UINT4 debug
)
{
	/* Place the LALparams into pWF */
	wf->LALparams = LALParams;

	/* Copy model version to struct */
	wf->IMRPhenomXInspiralPhaseVersion      = XLALSimInspiralWaveformParamsLookupPhenomXInspiralPhaseVersion(LALParams);
	wf->IMRPhenomXIntermediatePhaseVersion  = XLALSimInspiralWaveformParamsLookupPhenomXIntermediatePhaseVersion(LALParams);
	wf->IMRPhenomXRingdownPhaseVersion      = XLALSimInspiralWaveformParamsLookupPhenomXRingdownPhaseVersion(LALParams);

	wf->IMRPhenomXInspiralAmpVersion        = XLALSimInspiralWaveformParamsLookupPhenomXInspiralAmpVersion(LALParams);
	wf->IMRPhenomXIntermediateAmpVersion    = XLALSimInspiralWaveformParamsLookupPhenomXIntermediateAmpVersion(LALParams);
	wf->IMRPhenomXRingdownAmpVersion        = XLALSimInspiralWaveformParamsLookupPhenomXRingdownAmpVersion(LALParams);

	wf->IMRPhenomXPNRUseTunedCoprec = XLALSimInspiralWaveformParamsLookupPhenomXPNRUseTunedCoprec(LALParams);
	// NOTE that the line below means that 33 tuning can only be on IFF 22 tuning is on
	wf->IMRPhenomXPNRUseTunedCoprec33 = XLALSimInspiralWaveformParamsLookupPhenomXPNRUseTunedCoprec(LALParams) * (wf->IMRPhenomXPNRUseTunedCoprec);

    wf->PhenomXOnlyReturnPhase = XLALSimInspiralWaveformParamsLookupPhenomXOnlyReturnPhase(LALParams);

	// Get toggle for forcing inspiral phase and phase derivative alignment with XHM/AS
	INT4 PNRForceXHMAlignment = XLALSimInspiralWaveformParamsLookupPhenomXPNRForceXHMAlignment(LALParams);
	wf->IMRPhenomXPNRForceXHMAlignment = PNRForceXHMAlignment;

	wf->debug = PHENOMXDEBUG;

	if(PHENOMXDEBUG)
	{
		printf("\n");
		printf("Inspiral Amp Version		: %d\n",wf->IMRPhenomXInspiralAmpVersion);
		printf("Intermediate Amp Version	: %d\n",wf->IMRPhenomXIntermediateAmpVersion);
		printf("Ringdown Amp Version		: %d\n",wf->IMRPhenomXRingdownAmpVersion);
		printf("\n");
		printf("Inspiral Phase Version		: %d\n",wf->IMRPhenomXInspiralPhaseVersion);
		printf("Intermediate Phase Version	: %d\n",wf->IMRPhenomXIntermediatePhaseVersion);
		printf("Ringdown Phase Version		: %d\n",wf->IMRPhenomXRingdownPhaseVersion);
		printf("\n");
	}

	/*
	First perform a sanity check to make sure that a legitimate waveform version has been called. If not, fail immediately.
	Every time a new version for one of the regions is added, user must add case below.
	*/
	int InspPhaseFlag = wf->IMRPhenomXInspiralPhaseVersion;

	/* Phase : Inspiral */
	switch( InspPhaseFlag )
	{
		// Canonical TaylorF2 up to 3.5PN with 4 pseudo-PN coefficients
		case 104:
		{
			break;
		}
		// Canonical TaylorF2 up to 3.5PN with 5 pseudo-PN coefficients
		case 105:
		{
			break;
		}
		// Extended TaylorF2 up to 4.5PN with 4 pseudo-PN coefficients
		case 114:
		{
			break;
		}
		// Extended TaylorF2 up to 4.5PN with 5 pseudo-PN coefficients
		case 115:
		{
			break;
		}
		// We must pass a recognised inspiral phase version. Otherwise fail.
		default:
		{
			XLAL_ERROR(XLAL_EINVAL, "Error in IMRPhenomX_SetWaveformVariables: IMRPhenomXInspiralPhaseVersion not recognized. Recommended flag is 104.\n");
		}
	}

	int IntPhaseFlag = wf->IMRPhenomXIntermediatePhaseVersion;
	/* Phase : Intermediate */
	switch( IntPhaseFlag )
	{
		// 4 intermediate coefficients
		case 104:
		{
			break;
		}
		// 5 intermediate coefficients
		case 105:
		{
			break;
		}
		default:
		{
			XLAL_ERROR(XLAL_EINVAL,"Error in IMRPhenomX_SetWaveformVariables: IMRPhenomXIntermediatePhaseVersion not recognized. Recommended flag is 105.\n");
		}
	}

	int RDPhaseFlag = wf->IMRPhenomXRingdownPhaseVersion;
	/* Phase : Ringdown */
	switch( RDPhaseFlag )
	{
		// 5 ringdown coefficients
		case 105:
		{
			break;
		}
		default:
		{
			XLAL_ERROR(XLAL_EINVAL,"Error in IMRPhenomX_SetWaveformVariables: IMRPhenomXRingdownPhaseVersion not recognized. Recommended flag is 105.\n");
		}
	}

	int InsAmpFlag = wf->IMRPhenomXInspiralAmpVersion;
	/* Amplitude : Inspiral */
	switch( InsAmpFlag )
	{
		// Canonical TaylorF2 with 4 pseudo-PN coefficients
		case 103:
		{
			break;
		}
		default:
		{
			XLAL_ERROR(XLAL_EINVAL,"Error in IMRPhenomX_SetWaveformVariables: IMRPhenomXInspiralAmpVersion not recognized. Recommended flag is 103.\n");
		}
	}

	int IntAmpFlag = wf->IMRPhenomXIntermediateAmpVersion;
	/* Amplitude : Intermediate */
	switch( IntAmpFlag )
	{
		// Intermediate region constructed from 4th order polynomial
		case 1043:
		{
			break;
		}
		// Intermediate region constructed from a 4th order polynomial
		case 104:
		{
			break;
		}
		// Intermediate region constructed from a 5th order polynomial
		case 105:
		{
			break;
		}
		default:
		{
			XLAL_ERROR(XLAL_EINVAL,"Error in IMRPhenomX_SetWaveformVariables: IMRPhenomXIntermediateAmpVersion not recognized. Recommended flag is 104.\n");
		}
	}

	int RDAmpFlag = wf->IMRPhenomXRingdownAmpVersion;
	/* Amplitude : Ringdown */
	switch( RDAmpFlag )
	{
		case 103:
		{
			break;
		}
		default:
		{
			XLAL_ERROR(XLAL_EINVAL,"Error in IMRPhenomX_SetWaveformVariables: IMRPhenomXRingdownAmpVersion not recognized. Recommended flag is 103.\n");
		}
	}

	/* Rescale to mass in solar masses */
	REAL8 m1_In      = m1_SI / LAL_MSUN_SI; // Mass 1 in solar masses
	REAL8 m2_In      = m2_SI / LAL_MSUN_SI; // Mass 2 in solar masses

	/* Set matter parameters */
	REAL8 lambda1_In = 0, lambda2_In = 0, quadparam1_In = 1, quadparam2_In = 1; // Tidal deformabilities and quadrupole spin-induced deformation parameters

	if(XLALSimInspiralWaveformParamsLookupPhenomXTidalFlag(LALParams)!=0)
	{

      lambda1_In = XLALSimInspiralWaveformParamsLookupTidalLambda1(LALParams);
      lambda2_In = XLALSimInspiralWaveformParamsLookupTidalLambda2(LALParams);
      if( lambda1_In < 0 || lambda2_In < 0 )
                XLAL_ERROR(XLAL_EFUNC, "lambda1 = %f, lambda2 = %f. Both should be greater than zero for NRTidalv2 and NRTidalv3", lambda1_In, lambda2_In);

      int retcode;

      retcode = XLALSimInspiralSetQuadMonParamsFromLambdas(LALParams);
	  XLAL_CHECK(retcode == XLAL_SUCCESS, XLAL_EFUNC, "Failed to set quadparams from Universal relation.\n");

	  quadparam1_In = 1. + XLALSimInspiralWaveformParamsLookupdQuadMon1(LALParams);
	  quadparam2_In = 1. + XLALSimInspiralWaveformParamsLookupdQuadMon2(LALParams);
	}

	REAL8 m1, m2, chi1L, chi2L, lambda1, lambda2, quadparam1, quadparam2;

	/* Check if m1 >= m2, if not then swap masses/spins/lambdas/quadparams */
	if(m1_In >= m2_In)
	{
		chi1L      = chi1L_In;
		chi2L      = chi2L_In;
		m1         = m1_In;
		m2         = m2_In;
		lambda1    = lambda1_In;
		lambda2    = lambda2_In;
		quadparam1 = quadparam1_In;
		quadparam2 = quadparam2_In;
	}
	else
	{
		XLAL_PRINT_WARNING("Warning: m1 < m2, swapping the masses, spins, tidal deformabilities, and spin-induced quadrupole parameters.\n");
		chi1L      = chi2L_In;
		chi2L      = chi1L_In;
		m1         = m2_In;
		m2         = m1_In;
		lambda1    = lambda2_In;
		lambda2    = lambda1_In;
		quadparam1 = quadparam2_In;
		quadparam2 = quadparam1_In;
	}

	if(chi1L > 1.0)
	{
		IMRPhenomX_InternalNudge(chi1L,1.0,1e-6);
	}
	if(chi2L > 1.0)
	{
		IMRPhenomX_InternalNudge(chi2L,1.0,1e-6);
	}
	if(chi1L < -1.0)
	{
		IMRPhenomX_InternalNudge(chi1L,-1.0,1e-6);
	}
	if(chi2L < -1.0)
	{
		IMRPhenomX_InternalNudge(chi2L,-1.0,1e-6);
	}
	/* If spins are still unphysical after checking for small round-off errors, fail. */
	if( chi1L > 1.0 || chi1L < -1.0 || chi2L > 1.0 || chi2L < -1.0)
	{
		//XLALPrintError("Unphyiscal spins: must be in the range [-1,1].\n");
		XLAL_ERROR(XLAL_EDOM, "Unphysical spins requested: must obey the Kerr bound [-1,1].\n");
	}

	// Symmetric mass ratio
	REAL8 delta = fabs((m1 - m2) / (m1+m2));
	REAL8 eta   = fabs(0.25 * (1.0 - delta*delta) ); // use fabs to prevent negative sign due to roundoff
	REAL8 q   = ((m1 > m2) ? (m1 / m2) : (m2 / m1));

	/* If eta > 0.25, e.g. roundoff, then set to 0.25 */
	if(eta > 0.25)
	{
		eta = 0.25;
	}
	if(eta > 0.25 || eta < 0.0)
	{
		XLAL_ERROR(XLAL_EDOM,"Unphysical eta: must be between 0 and 0.25\n");
	}
	if(eta == 0.25) q = 1.;
	//if(eta < MAX_ALLOWED_ETA)
	//{
	//	XLAL_PRINT_WARNING("Warning: The model is not calibrated against numerical-relativity mass-ratios greater than 18.\n");
	//}

	/* Check that masses are correctly ordered. Swap if necessary.
	return_code = IMRPhenomX_CheckMassOrdering(m1,m2,chi1L,chi2L);
	XLAL_CHECK(XLAL_SUCCESS == return_code, XLAL_EFUNC, "Failure: IMRPhenomX_CheckMassOrdering");
	*/



	/* Masses definitions. Note that m1 and m2 are the component masses in solar masses. */
	wf->m1_SI     = m1 * LAL_MSUN_SI;																						// Mass 1 (larger) in SI units
	wf->m2_SI     = m2 * LAL_MSUN_SI;																						// Mass 2 (smaller) in SI units
	wf->q         = q;																								    			// Mass ratio >= 1
	wf->eta       = eta;                                                        // Symmetric mass ratio
	wf->Mtot_SI   = wf->m1_SI + wf->m2_SI;                                      // Total mass in SI units
	wf->Mtot      = m1 + m2;																										// Total mass in solar masss
	wf->m1        = m1 / (wf->Mtot);                         						        // Mass 1 (larger), dimensionless (i.e. m1 \in [0,1])
	wf->m2        = m2 / (wf->Mtot);						                                // Mass 2 (smaller), dimensionless
	wf->M_sec     = wf->Mtot * LAL_MTSUN_SI;                                    // Conversion factor Hz to geometric frequency, total mass in seconds
	wf->delta     = delta;                                                      // PN symmetry coefficient

	wf->eta2 = eta * eta;
	wf->eta3 = eta * wf->eta2;

	if(wf->debug)
	{
		printf("m1 (Msun) = %.6f\n",m1);
		printf("m2 (Msun) = %.6f\n",m2);
		printf("m1_SI     = %.6f\n",wf->m1_SI);
		printf("m2_SI     = %.6f\n",wf->m2_SI);
		printf("m1        = %.6f\n",wf->m1);
		printf("m2        = %.6f\n",wf->m2);
	}

	if(wf->q > 1000.0)
	{
		XLAL_PRINT_WARNING("Warning: The model is not supported for mass ratios > 1000.\n");
	}

	/* Spins */
	wf->chi1L     = chi1L;
	wf->chi2L     = chi2L;

	wf->chi1L2L   = chi1L * chi2L;

	/* Useful powers of spin */
	wf->chi1L2    = chi1L*chi1L;
	wf->chi1L3    = chi1L*chi1L*chi1L;

	wf->chi2L2    = chi2L*chi2L;
	wf->chi2L3    = chi2L*chi2L*chi2L;

	/* Spin parameterisations */
	wf->chiEff    = XLALSimIMRPhenomXchiEff(eta,chi1L,chi2L);
	wf->chiPNHat  = XLALSimIMRPhenomXchiPNHat(eta,chi1L,chi2L);
	wf->STotR     = XLALSimIMRPhenomXSTotR(eta,chi1L,chi2L);
	wf->dchi      = XLALSimIMRPhenomXdchi(chi1L,chi2L);
	wf->dchi_half = wf->dchi*0.5;

	wf->SigmaL    = (wf->chi2L * wf->m2) - (wf->chi1L * wf->m1); 										// SigmaL = (M/m2)*(S2.L) - (M/m2)*(S1.L)
	wf->SL        = wf->chi1L * (wf->m1 * wf->m1) + wf->chi2L * (wf->m2 * wf->m2);  // SL = S1.L + S2.L

	if(wf->debug)
	{
		printf("eta     : %.16e\n",wf->eta);
		printf("q       : %.16e\n",wf->q);
		printf("chi1L   : %.16e\n",wf->chi1L);
		printf("chi2L   : %.16e\n",wf->chi2L);
		printf("chi_eff : %.16e\n",wf->chiEff);
		printf("chi_hat : %.16e\n",wf->chiPNHat);
		printf("STotR   : %.16e\n",wf->STotR);
	}

	/* Matter parameters */
	wf->lambda1    = lambda1;
	wf->lambda2    = lambda2;
	wf->quadparam1 = quadparam1;
	wf->quadparam2 = quadparam2;
    wf->kappa2T=XLALSimNRTunedTidesComputeKappa2T(wf->m1_SI, wf->m2_SI, wf->lambda1, wf->lambda2);
    wf->fmerger    = XLALSimNRTunedTidesMergerFrequency(wf->Mtot, wf->kappa2T, wf->q);

	/* If no reference frequency is passed, set it to the starting GW frequency */
	wf->fRef      = (fRef == 0.0) ? f_min : fRef;
	wf->phiRef_In = phi0;
	wf->phi0      = phi0;							// Orbital phase at reference frequency (as passed from lalsimulation)
	wf->beta      = LAL_PI*0.5 - phi0; 				// Azimuthal angle of binary at reference frequency
	wf->phifRef   = 0.0;							// This is calculated later

	/* Geometric reference frequency */
	wf->MfRef     = XLALSimIMRPhenomXUtilsHztoMf(wf->fRef,wf->Mtot);
	wf->piM       = LAL_PI * wf->M_sec;
	wf->v_ref     = cbrt(wf->piM * wf->fRef);

	wf->deltaF    = deltaF;
    wf->deltaMF   = XLALSimIMRPhenomXUtilsHztoMf(wf->deltaF,wf->Mtot);

	/* Define the default end of the waveform as: 0.3 Mf. This value is chosen such that the 44 mode also shows the ringdown part. */
	wf->fCutDef   = 0.3;
	// If chieff is very high the ringdown of the 44 is almost cut it out when using 0.3, so we increase a little bit the cut of freq up 0.33.
	if(wf->chiEff > 0.99) wf->fCutDef = 0.33;

	/* Minimum and maximum frequency */
	wf->fMin      = f_min;
	wf->fMax      = f_max;
    wf->MfMax     = XLALSimIMRPhenomXUtilsHztoMf(wf->fMax,wf->Mtot);

	/* Convert fCut to physical cut-off frequency */
	wf->fCut      = wf->fCutDef / wf->M_sec;

	/* Sanity check that minimum start frequency is less than cut-off frequency */
	if (wf->fCut <= wf->fMin)
	{
		XLALPrintError("(fCut = %g Hz) <= f_min = %g\n", wf->fCut, wf->fMin);
	}

	if(wf->debug)
	{
		printf("fRef : %.6f\n",wf->fRef);
		printf("phi0 : %.6f\n",wf->phi0);
		printf("fCut : %.6f\n",wf->fCut);
		printf("fMin : %.6f\n",wf->fMin);
		printf("fMax : %.6f\n",wf->fMax);
	}

	/* By default f_max_prime is f_max. If fCut < fMax, then use fCut, i.e. waveform up to fMax will be zeros */
	wf->f_max_prime   = wf->fMax;
	wf->f_max_prime   = wf->fMax ? wf->fMax : wf->fCut;
	wf->f_max_prime   = (wf->f_max_prime > wf->fCut) ? wf->fCut : wf->f_max_prime;

	if (wf->f_max_prime <= wf->fMin)
	{
		XLALPrintError("f_max <= f_min\n");
	}

	if(wf->debug)
	{
		printf("fMin        = %.6f\n",wf->fMin);
		printf("fMax        = %.6f\n",wf->fMax);
		printf("f_max_prime = %.6f\n",wf->f_max_prime);
	}

	/* Final Mass and Spin */
	// NOTE: These are only default values
	wf->Mfinal    = XLALSimIMRPhenomXFinalMass2017(wf->eta,wf->chi1L,wf->chi2L);
	wf->afinal    = XLALSimIMRPhenomXFinalSpin2017(wf->eta,wf->chi1L,wf->chi2L);

	/* (500) Set default values of physically specific final spin parameters for use with PNR/XCP */
	wf->afinal_nonprec = wf->afinal;     // NOTE: This is only a default value; see LALSimIMRPhenomX_precession.c
	wf->afinal_prec    = wf->afinal;     // NOTE: This is only a default value; see LALSimIMRPhenomX_precession.c

	/* Ringdown and damping frequency of final BH */
#if QNMfits == 1
	wf->fRING     = evaluate_QNMfit_fring22(wf->afinal) / (wf->Mfinal);
	wf->fDAMP     = evaluate_QNMfit_fdamp22(wf->afinal) / (wf->Mfinal);
#else
	wf->fRING     = interpolateQNMData_fring_22(wf->afinal) / (wf->Mfinal);
	wf->fDAMP     = interpolateQNMData_fdamp_22(wf->afinal) / (wf->Mfinal);
#endif


	if(wf->debug)
	{
		printf("Mf  = %.6f\n",wf->Mfinal);
		printf("af  = %.6f\n",wf->afinal);
		printf("frd = %.6f\n",wf->fRING);
		printf("fda = %.6f\n",wf->fDAMP);
	}

	if(wf->Mfinal > 1.0)
	{
		XLAL_ERROR(XLAL_EDOM, "IMRPhenomX_FinalMass2018: Final mass > 1.0 not physical.");
	}
	if(fabs(wf->afinal) > 1.0)
	{
		XLAL_ERROR(XLAL_EDOM, "IMRPhenomX_FinalSpin2018: Final spin > 1.0 is not physical.");
	}

	/* Fit to the hybrid minimum energy circular orbit (MECO), Cabero et al, Phys.Rev. D95 (2017) */
	wf->fMECO       = XLALSimIMRPhenomXfMECO(wf->eta,wf->chi1L,wf->chi2L);

	/* Innermost stable circular orbit (ISCO), e.g. Ori et al, Phys.Rev. D62 (2000) 124022 */
	wf->fISCO       = XLALSimIMRPhenomXfISCO(wf->afinal);

	if(wf->debug)
	{
		printf("fMECO = %.6f\n",wf->fMECO);
		printf("fISCO = %.6f\n",wf->fISCO);
	}

	if(wf->fMECO > wf->fISCO)
	{
		/* If MECO > fISCO, throw an error - this may be the case for very corner cases in the parameter space (e.g. q ~ 1000, chi ~ 0.999) */
		XLAL_ERROR(XLAL_EDOM, "Error: f_MECO cannot be greater than f_ISCO.");
	}

	/* Distance and inclination */
	wf->distance    = distance;
	wf->inclination = inclination;

	/* Amplitude normalization */
	wf->amp0        = wf->Mtot * LAL_MRSUN_SI * wf->Mtot * LAL_MTSUN_SI / wf->distance;
	wf->ampNorm     = sqrt(2.0/3.0) * sqrt(wf->eta) * powers_of_lalpi.m_one_sixth;

	if(wf->debug)
	{
		printf("\n\neta     = %.6f\n",wf->eta);
		printf("\nampNorm = %e\n",wf->ampNorm);
		printf("amp0 : %e",wf->amp0);
	}

	/* Phase normalization */
	wf->dphase0     = 5.0 / (128.0 * pow(LAL_PI,5.0/3.0));

	if(wf->debug)
	{
		printf("\n\n **** Sanity checks complete. Waveform struct has been initialized. **** \n\n");
	}

	/* Set nonprecessing value of select precession quantities (PNRUseTunedCoprec)*/
	wf->chiTot_perp = 0.0;
	wf->chi_p = 0.0;
	wf->theta_LS = 0.0;
	wf->a1 = 0.0;
	wf->PNR_DEV_PARAMETER = 0.0;
	wf->PNR_SINGLE_SPIN = 0;
	wf->MU1 = 0;
	wf->MU2 = 0;
	wf->MU3 = 0;
	wf->MU4 = 0;
	wf->NU0 = 0;
	wf->NU4 = 0;
	wf->NU5 = 0;
	wf->NU6 = 0;
	wf->ZETA1 = 0;
	wf->ZETA2 = 0;
	wf->fRINGEffShiftDividedByEmm = 0;

	wf->f_inspiral_align = 0.0;
	wf->XAS_dphase_at_f_inspiral_align = 0.0;
	wf->XAS_phase_at_f_inspiral_align = 0.0;
	wf->XHM_dphase_at_f_inspiral_align = 0.0;
	wf->XHM_phase_at_f_inspiral_align = 0.0;

	wf->betaRD = 0.0;
	wf->fRING22_prec = 0.0;
	wf->fRINGCP = 0.0;
	wf->pnr_window = 0.0;

	wf->APPLY_PNR_DEVIATIONS = 0;

	return XLAL_SUCCESS;
}

int IMRPhenomXGetAmplitudeCoefficients(
	IMRPhenomXWaveformStruct *pWF,
	IMRPhenomXAmpCoefficients *pAmp
)
{
	INT4 debug = PHENOMXDEBUG;

	REAL8 V1, V2, V3, V4;
	REAL8 F1, F2, F3, F4;
	REAL8 d1, d4;

	/* Declare local spin variables for convenience */
	REAL8 chi1L         = pWF->chi1L;
	REAL8 chi2L         = pWF->chi2L;
	REAL8 chi1L2        = pWF->chi1L2;
	REAL8 chi1L3        = pWF->chi1L3;
	REAL8 chi2L2        = pWF->chi2L2;

	// Local PN symmetry parameter
	REAL8 delta         = pWF->delta;

	// Local powers of the symmetric mass ratio
	REAL8 eta           = pWF->eta;
	REAL8 eta2          = pWF->eta2;
	REAL8 eta3          = pWF->eta3;

	if(debug)
	{
		printf("chi1L  = %.6f\n",chi1L);
		printf("chi1L2 = %.6f\n",chi1L2);
		printf("chi1L3 = %.6f\n",chi1L3);
		printf("chi2L2 = %.6f\n",chi2L2);
		printf("delta  = %.6f\n",delta);
		printf("eta2   = %.6f\n",eta2);
		printf("eta3   = %.6f\n",eta3);
	}

	// Phenomenological ringdown coefficients: note that \gamma2 = \lambda in arXiv:2001.11412 and \gamma3 = \sigma in arXiv:2001.11412
	pAmp->gamma2        = IMRPhenomX_Ringdown_Amp_22_gamma2(pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXRingdownAmpVersion);
	pAmp->gamma3        = IMRPhenomX_Ringdown_Amp_22_gamma3(pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXRingdownAmpVersion);

	/* Apply modifications in the style of PhenomDCP: https://arxiv.org/pdf/2107.08876.pdf (500)*/
	// pAmp->gamma2 = pAmp->gamma2  +  ( pWF->PNR_DEV_PARAMETER * pWF->MU2 );
	pAmp->gamma3 = pAmp->gamma3  +  ( pWF->PNR_DEV_PARAMETER * pWF->MU2 );

	/* Get peak ringdown frequency, Eq. 5.14 in arXiv:2001.11412: Abs[fring + fdamp * gamma3 * (Sqrt[1 - gamma2^2] - 1)/gamma2 ] */
	pAmp->fAmpRDMin     = IMRPhenomX_Ringdown_Amp_22_PeakFrequency(pAmp->gamma2,pAmp->gamma3,pWF->fRING,pWF->fDAMP,pWF->IMRPhenomXRingdownAmpVersion);

	// Value of v1RD at fAmpRDMin is used to calcualte gamma1 \propto the amplitude of the deformed Lorentzian
	pAmp->v1RD          = IMRPhenomX_Ringdown_Amp_22_v1(pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXRingdownAmpVersion);
	F1                  = pAmp->fAmpRDMin;

	/* Apply modifications in the style of PhenomDCP: https://arxiv.org/pdf/2107.08876.pdf -- here we modify the ringdown amplitude (500)*/
	pAmp->v1RD = pAmp->v1RD  +  ( pWF->PNR_DEV_PARAMETER * pWF->MU1 );

	// Solving linear equation: ansatzRD(F1) == v1
	pAmp->gamma1 = ( pAmp->v1RD / (pWF->fDAMP * pAmp->gamma3) ) * (F1*F1 - 2.0*F1*pWF->fRING + pWF->fRING*pWF->fRING + pWF->fDAMP*pWF->fDAMP*pAmp->gamma3*pAmp->gamma3)
													* exp( ( (F1 - pWF->fRING) * pAmp->gamma2 ) / (pWF->fDAMP * pAmp->gamma3) );

	/* Pre-cache these constants for the ringdown amplitude */
	pAmp->gammaR   = pAmp->gamma2 / (pWF->fDAMP * pAmp->gamma3);
	pAmp->gammaD2  = (pAmp->gamma3 * pWF->fDAMP) * (pAmp->gamma3 * pWF->fDAMP);
	pAmp->gammaD13 = pWF->fDAMP * pAmp->gamma1 * pAmp->gamma3;

	if(debug)
	{
		printf("V1     = %.6f\n",pAmp->v1RD);
		printf("gamma1 = %.6f\n",pAmp->gamma1);
		printf("gamma2 = %.6f\n",pAmp->gamma2);
		printf("gamma3 = %.6f\n",pAmp->gamma3);
		printf("fmax   = %.6f\n",pAmp->fAmpRDMin);
	}

	/* Inspiral-Intermediate transition frequencies, see Eq. 5.7 in arXiv:2001.11412 */
	pAmp->fAmpInsMin    = 0.0026;
	pAmp->fAmpInsMax    = pWF->fMECO + 0.25*(pWF->fISCO - pWF->fMECO);

	/* Inspiral-Intermediate matching frequency, Eq. 5.16 in arXiv:2001.11412 */
	pAmp->fAmpMatchIN   = pAmp->fAmpInsMax;

	/* Intermediate-Ringdown matching frequency */
	//pAmp->fAmpMatchIM   = pAmp->fAmpRDMin;

	/* End of intermerdiate region, Eq. 5.16 in arXiv:2001.11412 */
	pAmp->fAmpIntMax    = pAmp->fAmpRDMin;

	/* TaylorF2 PN Amplitude Coefficients from usual PN re-expansion of Hlm's */
	pAmp->pnInitial     = 1.0;
	pAmp->pnOneThird    = 0.0;
	pAmp->pnTwoThirds   = ( (-969 + 1804*eta)/672. ) * powers_of_lalpi.two_thirds;
	pAmp->pnThreeThirds = ( (81*(chi1L + chi2L) + 81*chi1L*delta - 81*chi2L*delta - 44*(chi1L + chi2L)*eta)/48. ) * powers_of_lalpi.itself;
	pAmp->pnFourThirds  = ( (-27312085 - 10287648*chi1L2*(1 + delta)
													+ 24*(428652*chi2L2*(-1 + delta) + (-1975055 + 10584*(81*chi1L2 - 94*chi1L*chi2L + 81*chi2L2))*eta
													+ 1473794*eta2))/8.128512e6 ) * powers_of_lalpi.one_third * powers_of_lalpi.itself;
	pAmp->pnFiveThirds  = ( (-6048*chi1L3*(-1 - delta + (3 + delta)*eta) + chi2L*(-((287213 + 6048*chi2L2)*(-1 + delta))
													+ 4*(-93414 + 1512*chi2L2*(-3 + delta) + 2083*delta)*eta - 35632*eta2)
													+ chi1L*(287213*(1 + delta) - 4*eta*(93414 + 2083*delta + 8908*eta))
													+ 42840*(-1 + 4*eta)*LAL_PI)/32256. ) * powers_of_lalpi.five_thirds;
	pAmp->pnSixThirds   = ( (-1242641879927 + 12.0*(28.0*(-3248849057.0
													+ 11088*(163199*chi1L2 - 266498*chi1L*chi2L + 163199*chi2L2))*eta2
													+ 27026893936*eta3 - 116424*(147117*(-(chi2L2*(-1.0 + delta)) + chi1L2*(1.0 + delta))
													+ 60928*(chi1L + chi2L + chi1L*delta - chi2L*delta)*LAL_PI)
													+ eta*(545384828789.0 - 77616*(638642*chi1L*chi2L + chi1L2*(-158633 + 282718*delta)
													- chi2L2*(158633.0 + 282718.0*delta) - 107520.0*(chi1L + chi2L)*LAL_PI
													+ 275520*powers_of_lalpi.two))))/6.0085960704e10 ) * powers_of_lalpi.two;


	/* These are the same order as the canonical pseudo-PN terms but we can add higher order PN information as and when available here */
	pAmp->pnSevenThirds = 0.0;
	pAmp->pnEightThirds = 0.0;
	pAmp->pnNineThirds  = 0.0;

	/* Generate Pseudo-PN Coefficients for Amplitude. */
	switch( pWF->IMRPhenomXInspiralAmpVersion)
	{
		case 103:
		{
			pAmp->CollocationValuesAmpIns[0] = 0.0; // This is not currently used but may be invoked in future recalibration. Leave here, its harmless.
			pAmp->CollocationValuesAmpIns[1] = IMRPhenomX_Inspiral_Amp_22_v2(pWF->eta,pWF->chiPNHat,pWF->dchi,pWF->delta,pWF->IMRPhenomXInspiralAmpVersion);
			pAmp->CollocationValuesAmpIns[2] = IMRPhenomX_Inspiral_Amp_22_v3(pWF->eta,pWF->chiPNHat,pWF->dchi,pWF->delta,pWF->IMRPhenomXInspiralAmpVersion);
			pAmp->CollocationValuesAmpIns[3] = IMRPhenomX_Inspiral_Amp_22_v4(pWF->eta,pWF->chiPNHat,pWF->dchi,pWF->delta,pWF->IMRPhenomXInspiralAmpVersion);

			pAmp->CollocationPointsAmpIns[0] = 0.25 * pAmp->fAmpMatchIN;
			pAmp->CollocationPointsAmpIns[1] = 0.50 * pAmp->fAmpMatchIN;
			pAmp->CollocationPointsAmpIns[2] = 0.75 * pAmp->fAmpMatchIN;
			pAmp->CollocationPointsAmpIns[3] = 1.00 * pAmp->fAmpMatchIN;
			break;
		}
		default:
		{
			XLAL_ERROR(XLAL_EINVAL, "Error in IMRPhenomXGetAmplitudeCoefficients: IMRPhenomXInspiralAmpVersion is not valid.\n");
		}
	}

	// Set collocation points and values
	V1 = pAmp->CollocationValuesAmpIns[1];
	V2 = pAmp->CollocationValuesAmpIns[2];
	V3 = pAmp->CollocationValuesAmpIns[3];
	F1 = pAmp->CollocationPointsAmpIns[1];
	F2 = pAmp->CollocationPointsAmpIns[2];
	F3 = pAmp->CollocationPointsAmpIns[3];

	if(debug)
	{
		printf("\nAmplitude pseudo PN coeffcients\n");
		printf("fAmpMatchIN = %.6f\n",pAmp->fAmpMatchIN);
		printf("V1   = %.6f\n",V1);
		printf("V2   = %.6f\n",V2);
		printf("V3   = %.6f\n",V3);
		printf("F1   = %e\n",F1);
		printf("F2   = %e\n",F2);
		printf("F3   = %e\n",F3);
	}

	/* Using the collocation points and their values, we solve a linear system of equations to recover the pseudo-PN coefficients */
	pAmp->rho1 = IMRPhenomX_Inspiral_Amp_22_rho1(V1,V2,V3,F1,F2,F3,pWF->IMRPhenomXInspiralAmpVersion);
	pAmp->rho2 = IMRPhenomX_Inspiral_Amp_22_rho2(V1,V2,V3,F1,F2,F3,pWF->IMRPhenomXInspiralAmpVersion);
	pAmp->rho3 = IMRPhenomX_Inspiral_Amp_22_rho3(V1,V2,V3,F1,F2,F3,pWF->IMRPhenomXInspiralAmpVersion);

	/* Set TaylorF2 pseudo-PN coefficients */
	pAmp->pnSevenThirds = pAmp->rho1;
	pAmp->pnEightThirds = pAmp->rho2;
	pAmp->pnNineThirds  = pAmp->rho3;

	if(debug)
	{
		printf("\nTaylorF2 PN Amplitude Coefficients\n");
		printf("pnTwoThirds   = %.6f\n",pAmp->pnTwoThirds);
		printf("pnThreeThirds = %.6f\n",pAmp->pnThreeThirds);
		printf("pnFourThirds  = %.6f\n",pAmp->pnFourThirds);
		printf("pnFiveThirds  = %.6f\n",pAmp->pnFiveThirds);
		printf("pnSixThirds   = %.6f\n",pAmp->pnSixThirds);
		printf("\n");
		printf("powers_of_lalpi.itself = %.6f\n",powers_of_lalpi.itself);
		printf("powers_of_lalpi.four_thirds = %.6f\n",powers_of_lalpi.four_thirds);
		printf("powers_of_lalpi.five_thirds = %.6f\n",powers_of_lalpi.five_thirds);
		printf("\n");
		printf("Pseudo-PN Amplitude Coefficients (Agrees with MMA).\n");
		printf("alpha1 = %.6f\n",pAmp->pnSevenThirds);
		printf("alpha2 = %.6f\n",pAmp->pnEightThirds);
		printf("alpha3  = %.6f\n",pAmp->pnNineThirds);
	}

	/* Now we can reconstruct the intermediate region, which depends on the derivative of the inspiral and ringdown fits at the boundaries */
	F1     = pAmp->fAmpMatchIN;
	F4     = pAmp->fAmpRDMin;

	// Initialise useful powers for the boundary points
	IMRPhenomX_UsefulPowers powers_of_F1;
	IMRPhenomX_Initialize_Powers(&powers_of_F1,F1);

	IMRPhenomX_UsefulPowers powers_of_F4;
	IMRPhenomX_Initialize_Powers(&powers_of_F4,F4);

	/* Calculate derivative of amplitude on boundary points */
	d1     = IMRPhenomX_Inspiral_Amp_22_DAnsatz(F1,pWF,pAmp);
	d4     = IMRPhenomX_Ringdown_Amp_22_DAnsatz(F4,pWF,pAmp);

	double inspF1 = IMRPhenomX_Inspiral_Amp_22_Ansatz(F1,&powers_of_F1,pWF,pAmp);
	double rdF4   = IMRPhenomX_Ringdown_Amp_22_Ansatz(F4,pWF,pAmp);

	/*
		Use d1 and d4 calculated above to get the derivative of the amplitude on the boundaries:

		D[ f^{7/6} Ah^{-1}[f] , f] = ( (7/6) * f^{1/6} / Ah[f] ) - ( f^{7/6} Ah'[f] ) / (Ah[f]^2)

		where d1 = Ah'[F1], inspF1 = Ah[F1] etc.

	*/
	d1     = ((7.0/6.0) * powers_of_F1.one_sixth / inspF1) - ( powers_of_F1.seven_sixths * d1 / (inspF1*inspF1) );
	d4     = ((7.0/6.0) * powers_of_F4.one_sixth / rdF4)   - ( powers_of_F4.seven_sixths * d4 / (rdF4*rdF4) );

	if(debug)
	{
		printf("d1 = %.6f\n",d1);
		printf("d4 = %.6f\n",d4);
	}

	/* Pass inverse of points to reconstruction function... */
	switch ( pWF->IMRPhenomXIntermediateAmpVersion )
	{
		case 1043:
		{
			// Use a 5th order polynomial in intermediate - great agreement to calibration points but poor extrapolation
			F2     = F1 + (1.0/3.0) * (F4-F1);
			F3     = F1 + (2.0/3.0) * (F4-F1);

			V1     = powers_of_F1.m_seven_sixths * inspF1;
			V2     = IMRPhenomX_Intermediate_Amp_22_v2(pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXIntermediateAmpVersion);
			V3     = IMRPhenomX_Intermediate_Amp_22_v3(pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXIntermediateAmpVersion);
			V4     = powers_of_F4.m_seven_sixths * rdF4;

			V1 		 = 1.0 / V1;
			V2 		 = 1.0 / V2;
			V3 		 = 1.0 / V3;
			V4 		 = 1.0 / V4;
			break;
		}
		case 104:
		{
			// Use a 4th order polynomial in intermediate - good extrapolation, recommended default fit
			F2     = F1 + (1.0/2.0) * (F4-F1);
			F3     = 0.0;

			V1     = powers_of_F1.m_seven_sixths * inspF1;
			V2     = IMRPhenomX_Intermediate_Amp_22_vA(pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXIntermediateAmpVersion);
			V4     = powers_of_F4.m_seven_sixths * rdF4;

			V1 		 = 1.0 / V1;
			V2 		 = 1.0 / V2;
			V3 		 = 0.0;
			V4 		 = 1.0 / V4;

			break;
		}
		case 105:
		{
			// Use a 5th order polynomial in intermediate - great agreement to calibration points but poor extrapolation
			F2     = F1 + (1.0/3.0) * (F4-F1);
			F3     = F1 + (2.0/3.0) * (F4-F1);

			V1     = powers_of_F1.m_seven_sixths * inspF1;
			V2     = IMRPhenomX_Intermediate_Amp_22_v2(pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXIntermediateAmpVersion);
			V3     = IMRPhenomX_Intermediate_Amp_22_v3(pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXIntermediateAmpVersion);
			V4     = powers_of_F4.m_seven_sixths * rdF4;

			V1 		 = 1.0 / V1;
			V2 		 = 1.0 / V2;
			V3 		 = 1.0 / V3;
			V4 		 = 1.0 / V4;
			break;
		}
		default:
		{
			XLAL_ERROR(XLAL_EINVAL, "Error: IMRPhenomXIntermediateAmpVersion is not valid.\n");
		}
	}

	/* Let's try directly modifying the (inverse) amplitude values used in the linear system with collocation points (500). See table I of https://arxiv.org/pdf/2001.11412.pdf. NOTE that the 4th order polynomial fit is used by default (500)*/
	V2 = V2 + ( pWF->PNR_DEV_PARAMETER * pWF->MU3 );
	// V3 = V3 + ( pWF->PNR_DEV_PARAMETER * pWF->MU4 ); // NOTE that V3 is not used by default in PhenomX, therefore V3 and its deviations have no effect

	if(debug)
	{
		printf("\nIntermediate Region: \n");
		printf("F1 = %.6f\n",F1);
		printf("F2 = %.6f\n",F2);
		printf("F3 = %.6f\n",F3);
		printf("F4 = %.6f\n",F4);
		printf("\n");
		printf("Insp@F1 = %.6f\n",IMRPhenomX_Inspiral_Amp_22_Ansatz(F1,&powers_of_F1,pWF,pAmp));
		printf("d1 = %.6f\n",d1);
		printf("d4 = %.6f\n",d4);
		printf("V1 = %.6f\n",V1);
		printf("V2 = %.6f\n",V2);
		printf("V3 = %.6f\n",V3);
		printf("V4 = %.6f\n",V4);
	}

	/*
	Reconstruct the phenomenological coefficients for the intermediate ansatz:
	*/
	pAmp->delta0 = IMRPhenomX_Intermediate_Amp_22_delta0(d1,d4,V1,V2,V3,V4,F1,F2,F3,F4,pWF->IMRPhenomXIntermediateAmpVersion);
	pAmp->delta1 = IMRPhenomX_Intermediate_Amp_22_delta1(d1,d4,V1,V2,V3,V4,F1,F2,F3,F4,pWF->IMRPhenomXIntermediateAmpVersion);
	pAmp->delta2 = IMRPhenomX_Intermediate_Amp_22_delta2(d1,d4,V1,V2,V3,V4,F1,F2,F3,F4,pWF->IMRPhenomXIntermediateAmpVersion);
	pAmp->delta3 = IMRPhenomX_Intermediate_Amp_22_delta3(d1,d4,V1,V2,V3,V4,F1,F2,F3,F4,pWF->IMRPhenomXIntermediateAmpVersion);
	pAmp->delta4 = IMRPhenomX_Intermediate_Amp_22_delta4(d1,d4,V1,V2,V3,V4,F1,F2,F3,F4,pWF->IMRPhenomXIntermediateAmpVersion);
	pAmp->delta5 = IMRPhenomX_Intermediate_Amp_22_delta5(d1,d4,V1,V2,V3,V4,F1,F2,F3,F4,pWF->IMRPhenomXIntermediateAmpVersion);

	if(debug)
	{
  	    printf("delta0 = %.6f\n",pAmp->delta0);
		printf("delta1 = %.6f\n",pAmp->delta1);
		printf("delta2 = %.6f\n",pAmp->delta2);
		printf("delta3 = %.6f\n",pAmp->delta3);
		printf("delta4 = %.6f\n",pAmp->delta4);
		printf("delta5 = %.6f\n",pAmp->delta5);
	}

	return XLAL_SUCCESS;
}

/*
 * Function to populate the IMRPhenomXPhaseCoefficients struct:
*/
int IMRPhenomXGetPhaseCoefficients(
	IMRPhenomXWaveformStruct *pWF,
	IMRPhenomXPhaseCoefficients *pPhase
)
{
	/* Get LALparams */
	LALDict *LALparams    = pWF->LALparams;
	const INT4 debug      = PHENOMXDEBUG;

	/* GSL objects for solving system of equations via LU decomposition */
	gsl_vector *b, *x;
	gsl_matrix *A;
	gsl_permutation *p;
	int s;

	REAL8 deltax;
	REAL8 xmin;
	REAL8 fi;

	REAL8 gpoints4[4]     = {0.0, 1.0/4.0, 3.0/4.0, 1.0};
	REAL8 gpoints5[5]     = {0.0, 1.0/2 - 1.0/(2*sqrt(2.0)), 1.0/2.0, 1.0/2 + 1.0/(2.0*sqrt(2.0)), 1.0};

	// Matching regions

	/* This is Eq. 5.11 in the paper */
	double fIMmatch = 0.6 * (0.5 * pWF->fRING + pWF->fISCO);

	/* This is the MECO frequency */
	double fINmatch = pWF->fMECO;

	/* This is Eq. 5.10 in the paper */
	double deltaf   = (fIMmatch - fINmatch) * 0.03;

	// Transition frequency is just below the MECO frequency and just above the RD fitting region

	/* These are defined in Eq. 7.7 and the text just below, f_H = fPhaseMatchIM and f_L = fPhaseMatchIN */
	pPhase->fPhaseMatchIN  = fINmatch - 1.0*deltaf;
	pPhase->fPhaseMatchIM  = fIMmatch + 0.5*deltaf;

	/* Defined in Eq. 7.4, this is f_L */
	pPhase->fPhaseInsMin  = 0.0026;

	/* Defined in Eq. 7.4, this is f_H */
	pPhase->fPhaseInsMax  = 1.020 * pWF->fMECO;

	/* Defined in Eq. 7.12, this is f_L */
	pPhase->fPhaseRDMin   = fIMmatch;

	/* Defined in Eq. 7.12, this is f_L */
	pPhase->fPhaseRDMax   = pWF->fRING + 1.25*pWF->fDAMP;

	pPhase->phiNorm    		= -(3.0 * powers_of_lalpi.m_five_thirds) / (128.0);

	/* For convenience, define some variables here */
	REAL8 chi1L           = pWF->chi1L;
	REAL8 chi2L           = pWF->chi2L;

	REAL8 chi1L2L         = chi1L * chi2L;

	REAL8 chi1L2          = pWF->chi1L * pWF->chi1L;
	REAL8 chi1L3          = pWF->chi1L * chi1L2;

	REAL8 chi2L2          = pWF->chi2L * pWF->chi2L;
	REAL8 chi2L3          = pWF->chi2L * chi2L2;

	REAL8 eta             = pWF->eta;
	REAL8 eta2            = eta*eta;
	REAL8 eta3            = eta*eta2;

	REAL8 delta           = pWF->delta;

	/* Pre-initialize all phenomenological coefficients */
	pPhase->a0 = 0.0;
	pPhase->a1 = 0.0;
	pPhase->a2 = 0.0;
	pPhase->a3 = 0.0;
	pPhase->a4 = 0.0;

	pPhase->b0 = 0.0;
	pPhase->b1 = 0.0;
	pPhase->b2 = 0.0;
	pPhase->b3 = 0.0;
	pPhase->b4 = 0.0;

	pPhase->c0 = 0.0;
	pPhase->c1 = 0.0;
	pPhase->c2 = 0.0;
	pPhase->c3 = 0.0;
	pPhase->c4 = 0.0;
	pPhase->cL = 0.0;

    pPhase->c2PN_tidal   = 0.;
    pPhase->c3PN_tidal   = 0.;
    pPhase->c3p5PN_tidal = 0.;

	pPhase->sigma0 = 0.0;
	pPhase->sigma1 = 0.0;
	pPhase->sigma2 = 0.0;
	pPhase->sigma3 = 0.0;
	pPhase->sigma4 = 0.0;
	pPhase->sigma5 = 0.0;


	/*
		The general strategy is to initialize a linear system of equations:

		A.x = b

		- A is a matrix with the coefficients of the ansatz evaluated at the collocation nodes.
		- b is a vector of the value of the collocation points
		- x is the solution vector (i.e. the coefficients) that we must solve for.

		We choose to do this using a standard LU decomposition.
	*/

	/* Generate list of collocation points */
	/*
	The Gauss-Chebyshev Points are given by:
	GCPoints[n] = Table[Cos[i * pi/n] , {i, 0, n}]

	SamplePoints[xmin,xmax,n] = {
		pWF->delta = xmax - xmin;
		gpoints = 0.5 * (1 + GCPoints[n-1]);

		return {xmin + pWF->delta * gpoints}
	}

	gpoints4 = [1.0, 3.0/4, 1.0/4, 0.0];
	gpoints5 = [1.0, 1.0/2 + 1.0/(2.0*sqrt(2.0)), 1.0/2, 1.0/2 - 1.0/(2*sqrt(2.0)), 0.]

	*/

	/*
	Ringdown phase collocation points:
	Here we only use the first N+1 points in the array where N = the
	number of pseudo PN terms.

	The size of the array is controlled by: N_MAX_COLLOCATION_POINTS_PHASE_RD

	Default is to use 5 collocation points.
	*/
	deltax = pPhase->fPhaseRDMax - pPhase->fPhaseRDMin;
	xmin   = pPhase->fPhaseRDMin;
	int i;

	double phaseRD; // This is used in intermediate phase reconstruction.

	if(debug)
	{
		printf("\n");
		printf("Solving system of equations for RD phase...\n");
	}

	// Initialize collocation points
	for(i = 0; i < 5; i++)
	{
		pPhase->CollocationPointsPhaseRD[i] = gpoints5[i] * deltax + xmin;
	}
	// Collocation point 4 is set to the ringdown frequency ~ dip in Lorentzian
	pPhase->CollocationPointsPhaseRD[3] = pWF->fRING;

	if(debug)
	{
		printf("Rigndown collocation points : \n");
		printf("F1 : %.6f\n",pPhase->CollocationPointsPhaseRD[0]);
		printf("F2 : %.6f\n",pPhase->CollocationPointsPhaseRD[1]);
		printf("F3 : %.6f\n",pPhase->CollocationPointsPhaseRD[2]);
		printf("F4 : %.6f\n",pPhase->CollocationPointsPhaseRD[3]);
		printf("F5 : %.6f\n",pPhase->CollocationPointsPhaseRD[4]);
	}

	switch(pWF->IMRPhenomXRingdownPhaseVersion)
	{
		case 105:
		{
			pPhase->NCollocationPointsRD = 5;
			break;
		}
		default:
		{
			XLAL_ERROR(XLAL_EINVAL, "Error: IMRPhenomXRingdownPhaseVersion is not valid.\n");
		}
	}

	if(debug)
	{
		printf("NCollRD = %d\n",pPhase->NCollocationPointsRD);
	}

	// Eq. 7.13 in arXiv:2001.11412
	double RDv4 = IMRPhenomX_Ringdown_Phase_22_v4(pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXRingdownPhaseVersion);

	/* These are the calibrated collocation points, as per Eq. 7.13 */
	pPhase->CollocationValuesPhaseRD[0] = IMRPhenomX_Ringdown_Phase_22_d12(pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXRingdownPhaseVersion);
	pPhase->CollocationValuesPhaseRD[1] = IMRPhenomX_Ringdown_Phase_22_d24(pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXRingdownPhaseVersion);
	pPhase->CollocationValuesPhaseRD[2] = IMRPhenomX_Ringdown_Phase_22_d34( pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXRingdownPhaseVersion);
	pPhase->CollocationValuesPhaseRD[3] = RDv4;
	pPhase->CollocationValuesPhaseRD[4] = IMRPhenomX_Ringdown_Phase_22_d54(pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXRingdownPhaseVersion);

	/* v_j = d_{j4} + v4 */
	pPhase->CollocationValuesPhaseRD[4] = pPhase->CollocationValuesPhaseRD[4] + pPhase->CollocationValuesPhaseRD[3]; // v5 = d54  + v4
	pPhase->CollocationValuesPhaseRD[2] = pPhase->CollocationValuesPhaseRD[2] + pPhase->CollocationValuesPhaseRD[3]; // v3 = d34  + v4
	pPhase->CollocationValuesPhaseRD[1] = pPhase->CollocationValuesPhaseRD[1] + pPhase->CollocationValuesPhaseRD[3]; // v2 = d24  + v4
	pPhase->CollocationValuesPhaseRD[0] = pPhase->CollocationValuesPhaseRD[0] + pPhase->CollocationValuesPhaseRD[1]; // v1 = d12  + v2

	// Debugging information. Leave for convenience later on.
	if(debug)
	{
		printf("\n");
		printf("Ringdown Collocation Points: \n");
		printf("v1 : %.6f\n",pPhase->CollocationValuesPhaseRD[0]);
		printf("v2 : %.6f\n",pPhase->CollocationValuesPhaseRD[1]);
		printf("v3 : %.6f\n",pPhase->CollocationValuesPhaseRD[2]);
		printf("v4 : %.6f\n",pPhase->CollocationValuesPhaseRD[3]);
		printf("v5 : %.6f\n",pPhase->CollocationValuesPhaseRD[4]);
		printf("\n");
	}

	phaseRD = pPhase->CollocationValuesPhaseRD[0];

	p = gsl_permutation_alloc(pPhase->NCollocationPointsRD);
	b = gsl_vector_alloc(pPhase->NCollocationPointsRD);
	x = gsl_vector_alloc(pPhase->NCollocationPointsRD);
	A = gsl_matrix_alloc(pPhase->NCollocationPointsRD,pPhase->NCollocationPointsRD);

	/*
	Populate the b vector
	*/
	gsl_vector_set(b,0,pPhase->CollocationValuesPhaseRD[0]);
	gsl_vector_set(b,1,pPhase->CollocationValuesPhaseRD[1]);
	gsl_vector_set(b,2,pPhase->CollocationValuesPhaseRD[2]);
	gsl_vector_set(b,3,pPhase->CollocationValuesPhaseRD[3]);
	gsl_vector_set(b,4,pPhase->CollocationValuesPhaseRD[4]);

	/*
			Eq. 7.12 in arXiv:2001.11412

			ansatzRD(f) = a_0 + a_1 f^(-1/3) + a_2 f^(-2) + a_3 f^(-3) + a_4 f^(-4) + ( aRD ) / ( (f_damp^2 + (f - f_ring)^2 ) )

			Canonical ansatz sets a_3 to 0.
	*/

	/*
		We now set up and solve a linear system of equations.
		First we populate the matrix A_{ij}
	*/

	/* Note that ff0 is always 1 */
	REAL8 ff, invff, ff0, ff1, ff2, ff3, ff4;

	/* A_{0,i} */
	ff    = pPhase->CollocationPointsPhaseRD[0];
	invff = 1.0 / ff;
	ff1   = cbrt(invff);   // f^{-1/3}
	ff2   = invff * invff; // f^{-2}
	ff3   = ff2 * ff2;		 // f^{-4}
	ff4   = -(pWF->dphase0) / (pWF->fDAMP*pWF->fDAMP + (ff - pWF->fRING)*(ff - pWF->fRING));
	gsl_matrix_set(A,0,0,1.0); // Constant
	gsl_matrix_set(A,0,1,ff1); // f^(-1/3) term
	gsl_matrix_set(A,0,2,ff2); // f^(-2) term
	gsl_matrix_set(A,0,3,ff3); // f^(-4) term
	gsl_matrix_set(A,0,4,ff4); // Lorentzian term

	if(debug)
	{
		printf("For row 0: a0 + a1 %.6f + a2 %.6f + a4 %.6f + aRD %.6f\n",ff1,ff2,ff3,ff4);
	}

	/* A_{1,i} */
	ff    = pPhase->CollocationPointsPhaseRD[1];
	invff = 1.0 / ff;
	ff1   = cbrt(invff);
	ff2   = invff * invff;
	ff3   = ff2 * ff2;
	ff4   = -(pWF->dphase0) / (pWF->fDAMP*pWF->fDAMP + (ff - pWF->fRING)*(ff - pWF->fRING));
	gsl_matrix_set(A,1,0,1.0);
	gsl_matrix_set(A,1,1,ff1);
	gsl_matrix_set(A,1,2,ff2);
	gsl_matrix_set(A,1,3,ff3);
	gsl_matrix_set(A,1,4,ff4);

	/* A_{2,i} */
	ff    =  pPhase->CollocationPointsPhaseRD[2];
	invff = 1.0 / ff;
	ff1   = cbrt(invff);
	ff2   = invff * invff;
	ff3   = ff2 * ff2;
	ff4   = -(pWF->dphase0) / (pWF->fDAMP*pWF->fDAMP + (ff - pWF->fRING)*(ff - pWF->fRING));
	gsl_matrix_set(A,2,0,1.0);
	gsl_matrix_set(A,2,1,ff1);
	gsl_matrix_set(A,2,2,ff2);
	gsl_matrix_set(A,2,3,ff3);
	gsl_matrix_set(A,2,4,ff4);

	/* A_{3,i} */
	ff    = pPhase->CollocationPointsPhaseRD[3];
	invff = 1.0 / ff;
	ff1   = cbrt(invff);
	ff2   = invff * invff;
	ff3   = ff2 * ff2;
	ff4   = -(pWF->dphase0) / (pWF->fDAMP*pWF->fDAMP + (ff - pWF->fRING)*(ff - pWF->fRING));
	gsl_matrix_set(A,3,0,1.0);
	gsl_matrix_set(A,3,1,ff1);
	gsl_matrix_set(A,3,2,ff2);
	gsl_matrix_set(A,3,3,ff3);
	gsl_matrix_set(A,3,4,ff4);

	/* A_{4,i} */
	ff    = pPhase->CollocationPointsPhaseRD[4];
	invff = 1.0 / ff;
	ff1   = cbrt(invff);
	ff2   = invff * invff;
	ff3   = ff2 * ff2;
	ff4   = -(pWF->dphase0) / (pWF->fDAMP*pWF->fDAMP + (ff - pWF->fRING)*(ff - pWF->fRING));
	gsl_matrix_set(A,4,0,1.0);
	gsl_matrix_set(A,4,1,ff1);
	gsl_matrix_set(A,4,2,ff2);
	gsl_matrix_set(A,4,3,ff3);
	gsl_matrix_set(A,4,4,ff4);

	/* We now solve the system A x = b via an LU decomposition */
	gsl_linalg_LU_decomp(A,p,&s);
	gsl_linalg_LU_solve(A,p,b,x);

	pPhase->c0  = gsl_vector_get(x,0); // x[0]; 	// a0
	pPhase->c1  = gsl_vector_get(x,1); // x[1];		// a1
	pPhase->c2  = gsl_vector_get(x,2); // x[2]; 	// a2
	pPhase->c4  = gsl_vector_get(x,3); // x[3]; 	// a4
	pPhase->cRD = gsl_vector_get(x,4);
	pPhase->cL  = -(pWF->dphase0 * pPhase->cRD); // ~ x[4] // cL = - a_{RD} * dphase0

	/* Apply NR tuning for precessing cases (500) */
	pPhase->cL = pPhase->cL + ( pWF->PNR_DEV_PARAMETER * pWF->NU4 );
	// pPhase->c0 = pPhase->c0 + ( pWF->PNR_DEV_PARAMETER * pWF->NU0 );

	if(debug)
	{
		printf("\n");
		printf("Ringdown Coefficients: \n");
		printf("c0  : %.6f\n",pPhase->c0);
		printf("c1  : %.6f\n",pPhase->c1);
		printf("c2  : %.6f\n",pPhase->c2);
		printf("c4  : %e\n",pPhase->c4);
		printf("cRD : %.6f\n",gsl_vector_get(x,4));
		printf("d0  : %.6f\n",pWF->dphase0);
		printf("cL  : %e\n",pPhase->cL);
		printf("\n");

		printf("Freeing arrays...\n");
	}

	/* Tidy up in preparation for next GSL solve ... */
	gsl_vector_free(b);
	gsl_vector_free(x);
	gsl_matrix_free(A);
	gsl_permutation_free(p);

	/*
	Inspiral phase collocation points:
	Here we only use the first N+1 points in the array where N = the
	number of pseudo PN terms. E.g. for 4 pseudo-PN terms, we will
	need 4 collocation points. An ansatz with n free coefficients
	needs n pieces of information in order to constrain the ansatz.

	The size of the array is controlled by: N_MAX_COLLOCATION_POINTS_PHASE_INS

	Default is to use 4 pseudo-PN coefficients and hence 4 collocation points.

	GC points as per Eq. 7.4 and 7.5, where f_L = pPhase->fPhaseInsMin and f_H = pPhase->fPhaseInsMax
	*/
	deltax      = pPhase->fPhaseInsMax - pPhase->fPhaseInsMin;
	xmin        = pPhase->fPhaseInsMin;

	/*
			Set number of pseudo-PN coefficients:
				- If you add a new PN inspiral approximant, update with new version here.
	*/
	switch(pWF->IMRPhenomXInspiralPhaseVersion)
	{
		case 104:
		{
			pPhase->NPseudoPN = 4;
			pPhase->NCollocationPointsPhaseIns = 4;
			break;
		}
		case 105:
		{
			pPhase->NPseudoPN = 5;
			pPhase->NCollocationPointsPhaseIns = 5;
			break;
		}
		case 114:
		{
			pPhase->NPseudoPN = 4;
			pPhase->NCollocationPointsPhaseIns = 4;
			break;
		}
		case 115:
		{
			pPhase->NPseudoPN = 5;
			pPhase->NCollocationPointsPhaseIns = 5;
			break;
		}
		default:
		{
			XLAL_ERROR_REAL8(XLAL_EINVAL, "Error: IMRPhenomXInspiralPhaseVersion is not valid.\n");
		}
	}

	if(debug)
	{
		printf("\n");
		printf("NPseudoPN : %d\n",pPhase->NPseudoPN);
		printf("NColl : %d\n",pPhase->NCollocationPointsPhaseIns);
		printf("\n");
	}

	p = gsl_permutation_alloc(pPhase->NCollocationPointsPhaseIns);
	b = gsl_vector_alloc(pPhase->NCollocationPointsPhaseIns);
	x = gsl_vector_alloc(pPhase->NCollocationPointsPhaseIns);
	A = gsl_matrix_alloc(pPhase->NCollocationPointsPhaseIns,pPhase->NCollocationPointsPhaseIns);

	/*
	If we are using 4 pseudo-PN coefficients, call the routines below.
	The inspiral phase version is still passed to the individual functions.
	*/
	if(pPhase->NPseudoPN == 4)
	{
		// By default all models implemented use the following GC points.
		// If a new model is calibrated with different choice of collocation points, edit this.
		for(i = 0; i < pPhase->NCollocationPointsPhaseIns; i++)
		{
			fi = gpoints4[i] * deltax + xmin;
			pPhase->CollocationPointsPhaseIns[i] = fi;
		}

		// Calculate the value of the differences between the ith and 3rd collocation points at the GC nodes
		pPhase->CollocationValuesPhaseIns[0] = IMRPhenomX_Inspiral_Phase_22_d13(pWF->eta,pWF->chiPNHat,pWF->dchi,pWF->delta,pWF->IMRPhenomXInspiralPhaseVersion);
		pPhase->CollocationValuesPhaseIns[1] = IMRPhenomX_Inspiral_Phase_22_d23(pWF->eta,pWF->chiPNHat,pWF->dchi,pWF->delta,pWF->IMRPhenomXInspiralPhaseVersion);
		pPhase->CollocationValuesPhaseIns[2] = IMRPhenomX_Inspiral_Phase_22_v3( pWF->eta,pWF->chiPNHat,pWF->dchi,pWF->delta,pWF->IMRPhenomXInspiralPhaseVersion);
		pPhase->CollocationValuesPhaseIns[3] = IMRPhenomX_Inspiral_Phase_22_d43(pWF->eta,pWF->chiPNHat,pWF->dchi,pWF->delta,pWF->IMRPhenomXInspiralPhaseVersion);

		// Calculate the value of the collocation points at GC nodes via: v_i = d_i3 + v3
		pPhase->CollocationValuesPhaseIns[0] = pPhase->CollocationValuesPhaseIns[0] + pPhase->CollocationValuesPhaseIns[2];
		pPhase->CollocationValuesPhaseIns[1] = pPhase->CollocationValuesPhaseIns[1] + pPhase->CollocationValuesPhaseIns[2];
		pPhase->CollocationValuesPhaseIns[3] = pPhase->CollocationValuesPhaseIns[3] + pPhase->CollocationValuesPhaseIns[2];

		if(debug)
		{
			printf("\n");
			printf("Inspiral Collocation Points and Values:\n");
			printf("F1 : %.6f\n",pPhase->CollocationPointsPhaseIns[0]);
			printf("F2 : %.6f\n",pPhase->CollocationPointsPhaseIns[1]);
			printf("F3 : %.6f\n",pPhase->CollocationPointsPhaseIns[2]);
			printf("F4 : %.6f\n",pPhase->CollocationPointsPhaseIns[3]);
			printf("\n");
			printf("V1 : %.6f\n",pPhase->CollocationValuesPhaseIns[0]);
			printf("V2 : %.6f\n",pPhase->CollocationValuesPhaseIns[1]);
			printf("V3 : %.6f\n",pPhase->CollocationValuesPhaseIns[2]);
			printf("V4 : %.6f\n",pPhase->CollocationValuesPhaseIns[3]);
			printf("\n");
		}

		gsl_vector_set(b,0,pPhase->CollocationValuesPhaseIns[0]);
		gsl_vector_set(b,1,pPhase->CollocationValuesPhaseIns[1]);
		gsl_vector_set(b,2,pPhase->CollocationValuesPhaseIns[2]);
		gsl_vector_set(b,3,pPhase->CollocationValuesPhaseIns[3]);

		/* A_{0,i} */
		ff  = pPhase->CollocationPointsPhaseIns[0];
		ff1 = cbrt(ff);
		ff2 = ff1 * ff1;
		ff3 = ff;
		ff0 = 1.0;
		gsl_matrix_set(A,0,0,1.0);
		gsl_matrix_set(A,0,1,ff1);
		gsl_matrix_set(A,0,2,ff2);
		gsl_matrix_set(A,0,3,ff3);

		/* A_{1,i} */
		ff  = pPhase->CollocationPointsPhaseIns[1];
		ff1 = cbrt(ff);
		ff2 = ff1 * ff1;
		ff3 = ff;
		ff0 = 1.0;
		gsl_matrix_set(A,1,0,1.0);
		gsl_matrix_set(A,1,1,ff1);
		gsl_matrix_set(A,1,2,ff2);
		gsl_matrix_set(A,1,3,ff3);

		/* A_{2,i} */
		ff  = pPhase->CollocationPointsPhaseIns[2];
		ff1 = cbrt(ff);
		ff2 = ff1 * ff1;
		ff3 = ff;
		ff0 = 1.0;
		gsl_matrix_set(A,2,0,1.0);
		gsl_matrix_set(A,2,1,ff1);
		gsl_matrix_set(A,2,2,ff2);
		gsl_matrix_set(A,2,3,ff3);

		/* A_{3,i} */
		ff  = pPhase->CollocationPointsPhaseIns[3];
		ff1 = cbrt(ff);
		ff2 = ff1 * ff1;
		ff3 = ff;
		ff0 = 1.0;
		gsl_matrix_set(A,3,0,1.0);
		gsl_matrix_set(A,3,1,ff1);
		gsl_matrix_set(A,3,2,ff2);
		gsl_matrix_set(A,3,3,ff3);

		/* We now solve the system A x = b via an LU decomposition */
		gsl_linalg_LU_decomp(A,p,&s);
		gsl_linalg_LU_solve(A,p,b,x);

		/* Set inspiral phenomenological coefficients from solution to A x = b */
		pPhase->a0 = gsl_vector_get(x,0); // x[0]; // alpha_0
		pPhase->a1 = gsl_vector_get(x,1); // x[1]; // alpha_1
		pPhase->a2 = gsl_vector_get(x,2); // x[2]; // alpha_2
		pPhase->a3 = gsl_vector_get(x,3); // x[3]; // alpha_3
		pPhase->a4 = 0.0;

		/*
				PSEUDO PN TERMS WORK:
					- 104 works.
					- 105 not tested.
					- 114 not tested.
					- 115 not tested.
		*/
		if(debug)
		{
			printf("\n");
			printf("3pPN\n");
			printf("Inspiral Pseudo-PN Coefficients:\n");
			printf("a0 : %.6f\n",pPhase->a0);
			printf("a1 : %.6f\n",pPhase->a1);
			printf("a2 : %.6f\n",pPhase->a2);
			printf("a3 : %.6f\n",pPhase->a3);
			printf("a4 : %.6f\n",pPhase->a4);
			printf("\n");
		}

		/* Tidy up in preparation for next GSL solve ... */
		gsl_vector_free(b);
		gsl_vector_free(x);
		gsl_matrix_free(A);
		gsl_permutation_free(p);

	}
	else if(pPhase->NPseudoPN == 5)
	{
		// Using 5 pseudo-PN coefficients so set 5 collocation points
		for(i = 0; i < 5; i++)
		{
			fi = gpoints5[i] * deltax + xmin;
			pPhase->CollocationPointsPhaseIns[i] = fi;
		}
		pPhase->CollocationValuesPhaseIns[0] = IMRPhenomX_Inspiral_Phase_22_d13(pWF->eta,pWF->chiPNHat,pWF->dchi,pWF->delta,pWF->IMRPhenomXInspiralPhaseVersion);
		pPhase->CollocationValuesPhaseIns[1] = IMRPhenomX_Inspiral_Phase_22_d23(pWF->eta,pWF->chiPNHat,pWF->dchi,pWF->delta,pWF->IMRPhenomXInspiralPhaseVersion);
		pPhase->CollocationValuesPhaseIns[2] = IMRPhenomX_Inspiral_Phase_22_v3( pWF->eta,pWF->chiPNHat,pWF->dchi,pWF->delta,pWF->IMRPhenomXInspiralPhaseVersion);
		pPhase->CollocationValuesPhaseIns[3] = IMRPhenomX_Inspiral_Phase_22_d43(pWF->eta,pWF->chiPNHat,pWF->dchi,pWF->delta,pWF->IMRPhenomXInspiralPhaseVersion);
		pPhase->CollocationValuesPhaseIns[4] = IMRPhenomX_Inspiral_Phase_22_d53(pWF->eta,pWF->chiPNHat,pWF->dchi,pWF->delta,pWF->IMRPhenomXInspiralPhaseVersion);

		/* v_j = d_j3 + v_3 */
		pPhase->CollocationValuesPhaseIns[0] = pPhase->CollocationValuesPhaseIns[0] + pPhase->CollocationValuesPhaseIns[2];
		pPhase->CollocationValuesPhaseIns[1] = pPhase->CollocationValuesPhaseIns[1] + pPhase->CollocationValuesPhaseIns[2];
		pPhase->CollocationValuesPhaseIns[3] = pPhase->CollocationValuesPhaseIns[3] + pPhase->CollocationValuesPhaseIns[2];
		pPhase->CollocationValuesPhaseIns[4] = pPhase->CollocationValuesPhaseIns[4] + pPhase->CollocationValuesPhaseIns[2];

		if(debug)
		{
			printf("\n");
			printf("Inspiral Collocation Points and Values:\n");
			printf("F1 : %.6f\n",pPhase->CollocationPointsPhaseIns[0]);
			printf("F2 : %.6f\n",pPhase->CollocationPointsPhaseIns[1]);
			printf("F3 : %.6f\n",pPhase->CollocationPointsPhaseIns[2]);
			printf("F4 : %.6f\n",pPhase->CollocationPointsPhaseIns[3]);
			printf("F5 : %.6f\n",pPhase->CollocationPointsPhaseIns[4]);
			printf("\n");
			printf("V1 : %.6f\n",pPhase->CollocationValuesPhaseIns[0]);
			printf("V2 : %.6f\n",pPhase->CollocationValuesPhaseIns[1]);
			printf("V3 : %.6f\n",pPhase->CollocationValuesPhaseIns[2]);
			printf("V4 : %.6f\n",pPhase->CollocationValuesPhaseIns[3]);
			printf("V5 : %.6f\n",pPhase->CollocationValuesPhaseIns[4]);
			printf("\n");
		}

		gsl_vector_set(b,0,pPhase->CollocationValuesPhaseIns[0]);
		gsl_vector_set(b,1,pPhase->CollocationValuesPhaseIns[1]);
		gsl_vector_set(b,2,pPhase->CollocationValuesPhaseIns[2]);
		gsl_vector_set(b,3,pPhase->CollocationValuesPhaseIns[3]);
		gsl_vector_set(b,4,pPhase->CollocationValuesPhaseIns[4]);

		/* A_{0,i} */
		ff  = pPhase->CollocationPointsPhaseIns[0];
		ff1 = cbrt(ff);
		ff2 = ff1 * ff1;
		ff3 = ff;
		ff4 = ff * ff1;
		gsl_matrix_set(A,0,0,1.0);
		gsl_matrix_set(A,0,1,ff1);
		gsl_matrix_set(A,0,2,ff2);
		gsl_matrix_set(A,0,3,ff3);
		gsl_matrix_set(A,0,4,ff4);

		/* A_{1,i} */
		ff  = pPhase->CollocationPointsPhaseIns[1];
		ff1 = cbrt(ff);
		ff2 = ff1 * ff1;
		ff3 = ff;
		ff4 = ff * ff1;
		gsl_matrix_set(A,1,0,1.0);
		gsl_matrix_set(A,1,1,ff1);
		gsl_matrix_set(A,1,2,ff2);
		gsl_matrix_set(A,1,3,ff3);
		gsl_matrix_set(A,1,4,ff4);

		/* A_{2,i} */
		ff  = pPhase->CollocationPointsPhaseIns[2];
		ff1 = cbrt(ff);
		ff2 = ff1 * ff1;
		ff3 = ff;
		ff4 = ff * ff1;
		gsl_matrix_set(A,2,0,1.0);
		gsl_matrix_set(A,2,1,ff1);
		gsl_matrix_set(A,2,2,ff2);
		gsl_matrix_set(A,2,3,ff3);
		gsl_matrix_set(A,2,4,ff4);

		/* A_{3,i} */
		ff  = pPhase->CollocationPointsPhaseIns[3];
		ff1 = cbrt(ff);
		ff2 = ff1 * ff1;
		ff3 = ff;
		ff4 = ff * ff1;
		gsl_matrix_set(A,3,0,1.0);
		gsl_matrix_set(A,3,1,ff1);
		gsl_matrix_set(A,3,2,ff2);
		gsl_matrix_set(A,3,3,ff3);
		gsl_matrix_set(A,3,4,ff4);

		/* A_{4,i} */
		ff  = pPhase->CollocationPointsPhaseIns[4];
		ff1 = cbrt(ff);
		ff2 = ff1 * ff1;
		ff3 = ff;
		ff4 = ff * ff1;
		gsl_matrix_set(A,4,0,1.0);
		gsl_matrix_set(A,4,1,ff1);
		gsl_matrix_set(A,4,2,ff2);
		gsl_matrix_set(A,4,3,ff3);
		gsl_matrix_set(A,4,4,ff4);

		/* We now solve the system A x = b via an LU decomposition */
		gsl_linalg_LU_decomp(A,p,&s);
		gsl_linalg_LU_solve(A,p,b,x);

		/* Set inspiral phenomenological coefficients from solution to A x = b */
		pPhase->a0 = gsl_vector_get(x,0); // x[0];
		pPhase->a1 = gsl_vector_get(x,1); // x[1];
		pPhase->a2 = gsl_vector_get(x,2); // x[2];
		pPhase->a3 = gsl_vector_get(x,3); // x[3];
		pPhase->a4 = gsl_vector_get(x,4); // x[4];

		if(debug)
		{
			printf("\n");
			printf("4pPN\n");
			printf("Inspiral Pseudo-PN Coefficients:\n");
			printf("a0 : %.6f\n",pPhase->a0);
			printf("a1 : %.6f\n",pPhase->a1);
			printf("a2 : %.6f\n",pPhase->a2);
			printf("a3 : %.6f\n",pPhase->a3);
			printf("a4 : %.6f\n",pPhase->a4);
			printf("\n");
		}

		/* Tidy up in preparation for next GSL solve ... */
		gsl_vector_free(b);
		gsl_vector_free(x);
		gsl_matrix_free(A);
		gsl_permutation_free(p);
	}
	else
	{
		XLALPrintError("Error in ComputeIMRPhenomXWaveformVariables: NPseudoPN requested is not valid.\n");
	}

	/* The pseudo-PN coefficients are normalized such that: (dphase0 / eta) * f^{8/3} * a_j */
	/* So we must re-scale these terms by an extra factor of f^{-8/3} in the PN phasing */
	pPhase->sigma1 = (-5.0/3.0) * pPhase->a0;
	pPhase->sigma2 = (-5.0/4.0) * pPhase->a1;
	pPhase->sigma3 = (-5.0/5.0) * pPhase->a2;
	pPhase->sigma4 = (-5.0/6.0) * pPhase->a3;
	pPhase->sigma5 = (-5.0/7.0) * pPhase->a4;

	/* Initialize TaylorF2 PN coefficients  */
	pPhase->dphi0  = 0.0;
	pPhase->dphi1  = 0.0;
	pPhase->dphi2  = 0.0;
	pPhase->dphi3  = 0.0;
	pPhase->dphi4  = 0.0;
	pPhase->dphi5  = 0.0;
	pPhase->dphi6  = 0.0;
	pPhase->dphi7  = 0.0;
	pPhase->dphi8  = 0.0;
	pPhase->dphi9  = 0.0;
	pPhase->dphi10 = 0.0;
	pPhase->dphi11 = 0.0;
	pPhase->dphi12 = 0.0;

	pPhase->dphi5L = 0.0;
	pPhase->dphi6L = 0.0;
	pPhase->dphi8L = 0.0;
	pPhase->dphi9L = 0.0;

	pPhase->phi0   = 0.0;
	pPhase->phi1   = 0.0;
	pPhase->phi2   = 0.0;
	pPhase->phi3   = 0.0;
	pPhase->phi4   = 0.0;
	pPhase->phi5   = 0.0;
	pPhase->phi6   = 0.0;
	pPhase->phi7   = 0.0;
	pPhase->phi8   = 0.0;
	pPhase->phi9   = 0.0;
	pPhase->phi10  = 0.0;
	pPhase->phi11  = 0.0;
	pPhase->phi12  = 0.0;

	pPhase->phi5L  = 0.0;
	pPhase->phi6L  = 0.0;
	pPhase->phi8L  = 0.0;
	pPhase->phi9L  = 0.0;

	/* **** TaylorF2 PN Coefficients: Phase **** */

	/*
			- These are the PN coefficients normalised by: 3 / (128 * eta * [pi M f]^{5/3} ).
			- We add in powers of (M f)^{N/3} later but add powers of pi^{N/3} here
			- The log terms are *always* in terms of log(v), so we multiply by log(v) when summing PN phasing series.
			- We *do not* overwrite the PN phasing series with pseudo-PN terms. These are added separately.

			PN terms can be found in:
				- Marsat et al, CQG, 32, 085008, (2015)
				- Bohe et al, CQG, 32, 195010, (2015)
				- Bernard et al, PRD, 95, 044026, (2017)
				- Bernard et al, PRD, 93, 084037, (2016)
				- Damour et al, PRD, 89, 064058, (2014)
				- Damour et al, PRD, 95, 084005, (2017)
				- Bernard et al, PRD, 96, 104043, (2017)
				- Marchand et al, PRD, 97, 044023, (2018)
				- Marchand et al, CQG, 33, 244003, (2016)
				- Nagar et al, PRD, 99, 044007, (2019)
				- Messina et al, PRD, 97, 084016, (2018)
	*/

	/* Split into non-spinning and spin-dependent coefficients */
	UNUSED REAL8 phi0S = 0.0, phi1S = 0.0, phi2S = 0.0;
	REAL8 phi0NS  = 0.0,  phi1NS = 0.0,  phi2NS = 0.0;
	REAL8 phi3NS  = 0.0,  phi3S  = 0.0,  phi4NS = 0.0,  phi4S   = 0.0,  phi5NS  = 0.0,  phi5S  = 0.0;
	REAL8 phi5LNS = 0.0,  phi5LS = 0.0,  phi6NS = 0.0,  phi6S   = 0.0,  phi6LNS = 0.0,  phi6LS = 0.0;
	REAL8 phi7NS  = 0.0,  phi7S  = 0.0,  phi8NS = 0.0,  phi8S   = 0.0,  phi8LNS = 0.0;
	REAL8 phi8LS  = 0.0,  phi9NS = 0.0,  phi9S  = 0.0,  phi9LNS = 0.0,  phi9LS  = 0.0;


	/* Analytically known PN coefficients */
	/* Newtonian */
	phi0NS         = 1.0;

	/* ~~ 0.5 PN ~~ */
	phi1NS         = 0.0;

	/* ~~ 1.0 PN ~~ */
	/* 1.0PN, Non-Spinning */
	phi2NS         = (3715/756. + (55*eta)/9.) * powers_of_lalpi.two_thirds;

	/* ~~ 1.5 PN ~~ */
	/* 1.5PN, Non-Spinning */
	phi3NS         = -16.0 * powers_of_lalpi.two;
	/* 1.5PN, Spin-Orbit */
	phi3S          = ( (113*(chi1L + chi2L + chi1L*delta - chi2L*delta) - 76*(chi1L + chi2L)*eta)/6. ) * powers_of_lalpi.itself;

	/* ~~ 2.0 PN ~~ */
	/* 2.0PN, Non-Spinning */
	phi4NS         = ( 15293365/508032. + (27145*eta)/504. + (3085*eta2)/72. ) * powers_of_lalpi.four_thirds;
	/* 2.0PN, Spin-Spin */
	phi4S          = ( (-5*(81*chi1L2*(1 + delta - 2*eta) + 316*chi1L2L*eta - 81*chi2L2*(-1 + delta + 2*eta)))/16. ) * powers_of_lalpi.four_thirds;

	/* ~~ 2.5 PN ~~ */
	phi5NS         = 0.0;
	phi5S          = 0.0;

	/* ~~ 2.5 PN, Log Term ~~ */
	/* 2.5PN, Non-Spinning */
	phi5LNS        = ( (5*(46374 - 6552*eta)*LAL_PI)/4536. ) * powers_of_lalpi.five_thirds;
	/* 2.5PN, Spin-Orbit */
	phi5LS         = ( (-732985*(chi1L + chi2L + chi1L*delta - chi2L*delta) - 560*(-1213*(chi1L + chi2L)
												+ 63*(chi1L - chi2L)*delta)*eta + 85680*(chi1L + chi2L)*eta2)/4536. ) * powers_of_lalpi.five_thirds;

	/* ~~ 3.0 PN ~~ */
	/* 3.0 PN, Non-Spinning */
	phi6NS         = ( 11583231236531/4.69421568e9 - (5*eta*(3147553127 + 588*eta*(-45633 + 102260*eta)))/3.048192e6 - (6848*LAL_GAMMA)/21.
						- (640*powers_of_lalpi.two)/3. + (2255*eta*powers_of_lalpi.two)/12. - (13696*log(2))/21. - (6848*powers_of_lalpi.log)/63. ) * powers_of_lalpi.two;
	/* 3.0 PN, Spin-Orbit */
	phi6S          = ( (5*(227*(chi1L + chi2L + chi1L*delta - chi2L*delta) - 156*(chi1L + chi2L)*eta)*LAL_PI)/3. ) * powers_of_lalpi.two;
	/* 3.0 PN, Spin-Spin */
	phi6S         += ( (5*(20*chi1L2L*eta*(11763 + 12488*eta) + 7*chi2L2*(-15103*(-1 + delta) + 2*(-21683 + 6580*delta)*eta - 9808*eta2) -
							7*chi1L2*(-15103*(1 + delta) + 2*(21683 + 6580*delta)*eta + 9808*eta2)))/4032. ) * powers_of_lalpi.two;

	/* ~~ 3.0 PN, Log Term ~~ */
	phi6LNS        = (-6848/63.) * powers_of_lalpi.two;
	phi6LS         = 0.0;

	/* ~~ 3.5 PN ~~ */
	/* 3.5 PN, Non-Spinning */
	phi7NS         = ( (5*(15419335 + 168*(75703 - 29618*eta)*eta)*LAL_PI)/254016. ) * powers_of_lalpi.seven_thirds;
	/* 3.5 PN, Spin-Orbit */
	phi7S          = ( (5*(-5030016755*(chi1L + chi2L + chi1L*delta - chi2L*delta) + 4*(2113331119*(chi1L + chi2L) + 675484362*(chi1L - chi2L)*delta)*eta - 1008*(208433*(chi1L + chi2L) + 25011*(chi1L - chi2L)*delta)*eta2 + 90514368*(chi1L + chi2L)*eta3))/6.096384e6 ) * powers_of_lalpi.seven_thirds;
	/* 3.5 PN, Spin-Spin */
	phi7S         += ( -5*(57*chi1L2*(1 + delta - 2*eta) + 220*chi1L2L*eta - 57*chi2L2*(-1 + delta + 2*eta))*LAL_PI ) * powers_of_lalpi.seven_thirds;
	/* 3.5 PN, Cubic-in-Spin */
	phi7S         += ( (14585*(-(chi2L3*(-1 + delta)) + chi1L3*(1 + delta)) - 5*(chi2L3*(8819 - 2985*delta) + 8439*chi1L*chi2L2*(-1 + delta) - 8439*chi1L2*chi2L*(1 + delta) + chi1L3*(8819 + 2985*delta))*eta + 40*(chi1L + chi2L)*(17*chi1L2 - 14*chi1L2L + 17*chi2L2)*eta2)/48. ) * powers_of_lalpi.seven_thirds;

		/* ~~ 4.0 PN ~~ */
	/* 4.0 PN, Non-Spinning */
	phi8NS         = 0.0;
	/* 4.0 PN, Spin-Orbit */
	phi8S          = ( (-5*(1263141*(chi1L + chi2L + chi1L*delta - chi2L*delta) - 2*(794075*(chi1L + chi2L) + 178533*(chi1L - chi2L)*delta)*eta + 94344*(chi1L + chi2L)*eta2)*LAL_PI*(-1 + powers_of_lalpi.log))/9072. ) * powers_of_lalpi.eight_thirds;

	/* ~~ 4.0 PN, Log Term ~~ */
	/* 4.0 PN, log term, Non-Spinning */
	phi8LNS        = 0.0;
	/* 4.0 PN, log term, Spin-Orbit */
	phi8LS         = ((-5*(1263141*(chi1L + chi2L + chi1L*delta - chi2L*delta) - 2*(794075*(chi1L + chi2L) + 178533*(chi1L - chi2L)*delta)*eta
							+ 94344*(chi1L + chi2L)*eta2)*LAL_PI)/9072.) * powers_of_lalpi.eight_thirds;

	/* ~~ 4.5 PN ~~ */
	phi9NS         = 0.0;
	phi9S          = 0.0;

		/* ~~ 4.5 PN, Log Term ~~ */
	phi9LNS        = 0.0;
	phi9LS         = 0.0;

	/* This version of TaylorF2 contains an additional 4.5PN tail term and a LO-SS tail term at 3.5PN */
	if(pWF->IMRPhenomXInspiralPhaseVersion == 114 || pWF->IMRPhenomXInspiralPhaseVersion == 115)
	{
			/* 3.5PN, Leading Order Spin-Spin Tail Term */
			phi7S         += ( (5*(65*chi1L2*(1 + delta - 2*eta) + 252*chi1L2L*eta - 65*chi2L2*(-1 + delta + 2*eta))*LAL_PI)/4. ) * powers_of_lalpi.seven_thirds;

			/* 4.5PN, Tail Term */
			phi9NS        += ( (5*(-256 + 451*eta)*powers_of_lalpi.three)/6. + (LAL_PI*(105344279473163 + 700*eta*(-298583452147 + 96*eta*(99645337 + 14453257*eta)) -
																						12246091038720*LAL_GAMMA - 24492182077440*log(2.0)))/1.877686272e10 - (13696*LAL_PI*powers_of_lalpi.log)/63. ) * powers_of_lalpi.three;

			/* 4.5PN, Log Term */
			phi9LNS       += (  (-13696*LAL_PI)/63.0  ) * powers_of_lalpi.three;
	}

	/* 0.0 PN */
	pPhase->phi0   = phi0NS;

	/* 0.5 PN */
	pPhase->phi1   = phi1NS;

	/* 1.0 PN */
	pPhase->phi2   = phi2NS;

	/* 1.5 PN */
	pPhase->phi3   = phi3NS + phi3S;

	/* 2.0 PN */
	pPhase->phi4   = phi4NS + phi4S;

	/* 2.5 PN */
	pPhase->phi5   = phi5NS + phi5S;

	/* 2.5 PN, Log Terms */
	pPhase->phi5L  = phi5LNS + phi5LS;

	/* 3.0 PN */
	pPhase->phi6   = phi6NS + phi6S;

	/* 3.0 PN, Log Term */
	pPhase->phi6L  = phi6LNS + phi6LS;

	/* 3.5PN */
	pPhase->phi7   = phi7NS + phi7S;

	/* 4.0PN */
	pPhase->phi8   = phi8NS + phi8S;

	/* 4.0 PN, Log Terms */
	pPhase->phi8L  = phi8LNS + phi8LS;

	/* 4.5 PN */
	pPhase->phi9   = phi9NS + phi9S;

	/* 4.5 PN, Log Terms */
	pPhase->phi9L  = phi9LNS + phi9LS;


	if(debug)
	{
		printf("TaylorF2 PN Coefficients: \n");
		printf("phi0   : %.6f\n",pPhase->phi0);
		printf("phi1   : %.6f\n",pPhase->phi1);
		printf("phi2   : %.6f\n",pPhase->phi2);
		printf("phi3   : %.6f\n",pPhase->phi3);
		printf("phi4   : %.6f\n",pPhase->phi4);
		printf("phi5   : %.6f\n",pPhase->phi5);
		printf("phi6   : %.6f\n",pPhase->phi6);
		printf("phi7   : %.6f\n",pPhase->phi7);
		printf("phi8   : %.6f\n",pPhase->phi8);

		printf("phi5L  : %.6f\n",pPhase->phi5L);
		printf("phi6L  : %.6f\n",pPhase->phi6L);
		printf("phi8L  : %.6f\n",pPhase->phi8L);

		printf("phi8P  : %.6f\n",pPhase->sigma1);
		printf("phi9P  : %.6f\n",pPhase->sigma2);
		printf("phi10P : %.6f\n",pPhase->sigma3);
		printf("phi11P : %.6f\n",pPhase->sigma4);
		printf("phi12P : %.6f\n",pPhase->sigma5);
	}

	pPhase->phi_initial = - LAL_PI_4;

	/* **** TaylorF2 PN Coefficients: Normalized Phase Derivative **** */
	pPhase->dphi0  = pPhase->phi0;
	pPhase->dphi1  = 4.0 / 5.0 * pPhase->phi1;
	pPhase->dphi2  = 3.0 / 5.0 * pPhase->phi2;
	pPhase->dphi3  = 2.0 / 5.0 * pPhase->phi3;
	pPhase->dphi4  = 1.0 / 5.0 * pPhase->phi4;
	pPhase->dphi5  = -3.0 / 5.0 * pPhase->phi5L;
	pPhase->dphi6  = -1.0 / 5.0 * pPhase->phi6 - 3.0 / 5.0 * pPhase->phi6L;
	pPhase->dphi6L = -1.0 / 5.0 * pPhase->phi6L;
	pPhase->dphi7  = -2.0 / 5.0 * pPhase->phi7;
	pPhase->dphi8  = -3.0 / 5.0 * pPhase->phi8 - 3.0 / 5.0 * pPhase->phi8L;
	pPhase->dphi8L = -3.0 / 5.0 * pPhase->phi8L;
	pPhase->dphi9  = -4.0 / 5.0 * pPhase->phi9 - 3.0 / 5.0 * pPhase->phi9L;
	pPhase->dphi9L = -3.0 / 5.0 * pPhase->phi9L;

	if(debug)
	{
		printf("\nTaylorF2 PN Derivative Coefficients\n");
		printf("dphi0  : %.6f\n",pPhase->dphi0);
		printf("dphi1  : %.6f\n",pPhase->dphi1);
		printf("dphi2  : %.6f\n",pPhase->dphi2);
		printf("dphi3  : %.6f\n",pPhase->dphi3);
		printf("dphi4  : %.6f\n",pPhase->dphi4);
		printf("dphi5  : %.6f\n",pPhase->dphi5);
		printf("dphi6  : %.6f\n",pPhase->dphi6);
		printf("dphi7  : %.6f\n",pPhase->dphi7);
		printf("dphi8  : %.6f\n",pPhase->dphi8);
		printf("dphi9  : %.6f\n",pPhase->dphi9);
		printf("\n");
		printf("dphi6L : %.6f\n",pPhase->dphi6L);
		printf("dphi8L : %.6f\n",pPhase->dphi8L);
		printf("dphi9L : %.6f\n",pPhase->dphi9L);
	}

	/*
			Calculate phase at fmatchIN. This will be used as the collocation point for the intermediate fit.
			In practice, the transition point is just below the MECO frequency.
	*/
	if(debug)
	{
		printf("\nTransition frequency for ins to int : %.6f\n",pPhase->fPhaseMatchIN);
	}

	IMRPhenomX_UsefulPowers powers_of_fmatchIN;
	IMRPhenomX_Initialize_Powers(&powers_of_fmatchIN,pPhase->fPhaseMatchIN);

	double phaseIN;
	phaseIN  = pPhase->dphi0; 																	// f^{0/3}
	phaseIN += pPhase->dphi1 	* powers_of_fmatchIN.one_third; 								// f^{1/3}
	phaseIN += pPhase->dphi2 	* powers_of_fmatchIN.two_thirds; 								// f^{2/3}
	phaseIN += pPhase->dphi3 	* powers_of_fmatchIN.itself; 									// f^{3/3}
	phaseIN += pPhase->dphi4 	* powers_of_fmatchIN.four_thirds; 								// f^{4/3}
	phaseIN += pPhase->dphi5 	* powers_of_fmatchIN.five_thirds; 								// f^{5/3}
	phaseIN += pPhase->dphi6  	* powers_of_fmatchIN.two;										// f^{6/3}
	phaseIN += pPhase->dphi6L 	* powers_of_fmatchIN.two * powers_of_fmatchIN.log;				// f^{6/3}, Log[f]
	phaseIN += pPhase->dphi7  	* powers_of_fmatchIN.seven_thirds;								// f^{7/3}
	phaseIN += pPhase->dphi8  	* powers_of_fmatchIN.eight_thirds;								// f^{8/3}
	phaseIN += pPhase->dphi8L 	* powers_of_fmatchIN.eight_thirds * powers_of_fmatchIN.log;		// f^{8/3}
	phaseIN += pPhase->dphi9  	* powers_of_fmatchIN.three;										// f^{9/3}
	phaseIN += pPhase->dphi9L 	* powers_of_fmatchIN.three * powers_of_fmatchIN.log;			// f^{9/3}

	// Add pseudo-PN Coefficient
	phaseIN += ( 		pPhase->a0 * powers_of_fmatchIN.eight_thirds
								+ pPhase->a1 * powers_of_fmatchIN.three
								+ pPhase->a2 * powers_of_fmatchIN.eight_thirds * powers_of_fmatchIN.two_thirds
								+ pPhase->a3 * powers_of_fmatchIN.eight_thirds * powers_of_fmatchIN.itself
								+ pPhase->a4 * powers_of_fmatchIN.eight_thirds * powers_of_fmatchIN.four_thirds
							);

	phaseIN  = phaseIN * powers_of_fmatchIN.m_eight_thirds * pWF->dphase0;

	/*
	Intermediate phase collocation points:
	Here we only use the first N points in the array where N = the
	number of intermediate collocation points.

	The size of the array is controlled by: N_MAX_COLLOCATION_POINTS_PHASE_INT

	Default is to use 5 collocation points.

	See. Eq. 7.7 and 7.8 where f_H = pPhase->fPhaseMatchIM and f_L = pPhase->fPhaseMatchIN
	*/
	deltax      = pPhase->fPhaseMatchIM - pPhase->fPhaseMatchIN;
	xmin        = pPhase->fPhaseMatchIN;

	switch(pWF->IMRPhenomXIntermediatePhaseVersion)
	{
		case 104:
		{
			// Fourth order polynomial ansatz
			pPhase->NCollocationPointsInt = 4;
			break;
		}
		case 105:
		{
			// Fifth order polynomial ansatz
			pPhase->NCollocationPointsInt = 5;
			break;
		}
		default:
		{
			XLAL_ERROR(XLAL_EINVAL, "Error: IMRPhenomXIntermediatePhaseVersion is not valid.\n");
		}
	}

	if(debug)
	{
	printf("\nNColPointsInt : %d\n",pPhase->NCollocationPointsInt);
	}

	p = gsl_permutation_alloc(pPhase->NCollocationPointsInt);
	b = gsl_vector_alloc(pPhase->NCollocationPointsInt);
	x = gsl_vector_alloc(pPhase->NCollocationPointsInt);
	A = gsl_matrix_alloc(pPhase->NCollocationPointsInt,pPhase->NCollocationPointsInt);

	// Canonical intermediate model using 4 collocation points
	if(pWF->IMRPhenomXIntermediatePhaseVersion == 104)
	{
		// Using 4 collocation points in intermediate region
		for(i = 0; i < 4; i++)
		{
			fi = gpoints4[i] * deltax + xmin;

			pPhase->CollocationPointsPhaseInt[i] = fi;
		}

		// v2IM - v4RD. Using v4RD helps condition the fits with v4RD being very a robust fit.
		double v2IMmRDv4 = IMRPhenomX_Intermediate_Phase_22_v2mRDv4(pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXIntermediatePhaseVersion);

		// v3IM - v4RD. Using v4RD helps condition the fits with v4RD being very a robust fit.
		double v3IMmRDv4 = IMRPhenomX_Intermediate_Phase_22_v3mRDv4(pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXIntermediatePhaseVersion);

		// Direct fit to the collocation point at F2. We will take a weighted average of the direct and conditioned fit.
		double v2IM      = IMRPhenomX_Intermediate_Phase_22_v2(pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXIntermediatePhaseVersion);

		/* Evaluate collocation points */
		pPhase->CollocationValuesPhaseInt[0] = phaseIN;

		// Take a weighted average for these points? Can help condition the fit.
		pPhase->CollocationValuesPhaseInt[1] = 0.75*(v2IMmRDv4 + RDv4) + 0.25*v2IM;

		// Use just v2 - v4RD to reconstruct the fit?
		//pPhase->CollocationValuesPhaseInt[1] = v2IMmRDv4 + RDv4);

		pPhase->CollocationValuesPhaseInt[2] = v3IMmRDv4 + RDv4;

		pPhase->CollocationValuesPhaseInt[3] = phaseRD;

		/* A_{0,i} */
		ff  = pPhase->CollocationPointsPhaseInt[0];
		ff1 = pWF->fRING / ff;
		ff2 = ff1 * ff1;
		ff3 = ff1 * ff2;
		ff0 = (4 * pPhase->cL) / (4.0*pWF->fDAMP*pWF->fDAMP + (ff - pWF->fRING)*(ff - pWF->fRING));
		gsl_matrix_set(A,0,0,1.0);
		gsl_matrix_set(A,0,1,ff1);
		gsl_matrix_set(A,0,2,ff2);
		gsl_matrix_set(A,0,3,ff3);
		gsl_vector_set(b,0,pPhase->CollocationValuesPhaseInt[0] - ff0);

		/* A_{1,i} */
		ff  = pPhase->CollocationPointsPhaseInt[1];
		ff1 = 1.0 / (ff / pWF->fRING);
		ff2 = ff1 * ff1;
		ff3 = ff1 * ff2;
		ff0 = (4 * pPhase->cL) / (4.0*pWF->fDAMP*pWF->fDAMP + (ff - pWF->fRING)*(ff - pWF->fRING));
		gsl_matrix_set(A,1,0,1);
		gsl_matrix_set(A,1,1,ff1);
		gsl_matrix_set(A,1,2,ff2);
		gsl_matrix_set(A,1,3,ff3);
		gsl_vector_set(b,1,pPhase->CollocationValuesPhaseInt[1] - ff0);

		/* A_{2,i} */
		ff  = pPhase->CollocationPointsPhaseInt[2];
		ff1 = 1.0 / (ff / pWF->fRING);
		ff2 = ff1 * ff1;
		ff3 = ff1 * ff2;
		ff0 = (4 * pPhase->cL) / (4.0*pWF->fDAMP*pWF->fDAMP + (ff - pWF->fRING)*(ff - pWF->fRING));
		gsl_matrix_set(A,2,0,1);
		gsl_matrix_set(A,2,1,ff1);
		gsl_matrix_set(A,2,2,ff2);
		gsl_matrix_set(A,2,3,ff3);
		gsl_vector_set(b,2,pPhase->CollocationValuesPhaseInt[2] - ff0);

		/* A_{3,i} */
		ff  = pPhase->CollocationPointsPhaseInt[3];
		ff1 = 1.0 / (ff / pWF->fRING);
		ff2 = ff1 * ff1;
		ff3 = ff1 * ff2;
		ff0 = (4 * pPhase->cL) / (4.0*pWF->fDAMP*pWF->fDAMP + (ff - pWF->fRING)*(ff - pWF->fRING));
		gsl_matrix_set(A,3,0,1);
		gsl_matrix_set(A,3,1,ff1);
		gsl_matrix_set(A,3,2,ff2);
		gsl_matrix_set(A,3,3,ff3);
		gsl_vector_set(b,3,pPhase->CollocationValuesPhaseInt[3] - ff0);

		/* We now solve the system A x = b via an LU decomposition */
		gsl_linalg_LU_decomp(A,p,&s);
		gsl_linalg_LU_solve(A,p,b,x);

		/* Set intermediate phenomenological coefficients from solution to A x = b */
		pPhase->b0 = gsl_vector_get(x,0);                                                        // x[0] // Constant
		pPhase->b1 = gsl_vector_get(x,1) * pWF->fRING;                                           // x[1] // f^{-1}
		pPhase->b2 = gsl_vector_get(x,2) * pWF->fRING * pWF->fRING;                              // x[2] // f^{-2}
		//pPhase->b3 = 0.0;
		pPhase->b4 = gsl_vector_get(x,3) * pWF->fRING * pWF->fRING * pWF->fRING * pWF->fRING;    // x[3] // f^{-4}

		/* Tidy up in preparation for next GSL solve ... */
		gsl_vector_free(b);
		gsl_vector_free(x);
		gsl_matrix_free(A);
		gsl_permutation_free(p);
	}
	// Canonical intermediate model using 5 collocation points
	else if(pWF->IMRPhenomXIntermediatePhaseVersion == 105)
	{
		// Using 5 collocation points in intermediate region
		for(i = 0; i < 5; i++)
		{
			fi = gpoints5[i] * deltax + xmin;

			pPhase->CollocationPointsPhaseInt[i] = fi;
		}

		/* Evaluate collocation points */

		/* The first and last collocation points for the intermediate region are set from the inspiral fit and ringdown respectively */
		pPhase->CollocationValuesPhaseInt[0] = phaseIN;
		pPhase->CollocationValuesPhaseInt[4] = phaseRD;

		// v2IM - v4RD. Using v4RD helps condition the fits with v4RD being very a robust fit.
		double v2IMmRDv4 = IMRPhenomX_Intermediate_Phase_22_v2mRDv4(pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXIntermediatePhaseVersion);

		// v3IM - v4RD. Using v4RD helps condition the fits with v4RD being very a robust fit.
		double v3IMmRDv4 = IMRPhenomX_Intermediate_Phase_22_v3mRDv4(pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXIntermediatePhaseVersion);

		// Direct fit to the collocation point at F2. We will take a weighted average of the direct and conditioned fit.
		double v2IM      = IMRPhenomX_Intermediate_Phase_22_v2(pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXIntermediatePhaseVersion);

		// Take a weighted average for these points. Helps condition the fit.
		pPhase->CollocationValuesPhaseInt[1] = 0.75*(v2IMmRDv4 + RDv4) + 0.25*v2IM;

		pPhase->CollocationValuesPhaseInt[2] = v3IMmRDv4 + RDv4;
		pPhase->CollocationValuesPhaseInt[3] = IMRPhenomX_Intermediate_Phase_22_d43(pWF->eta,pWF->STotR,pWF->dchi,pWF->delta,pWF->IMRPhenomXIntermediatePhaseVersion);


		// Collocation points: v4 = d43 + v3
		pPhase->CollocationValuesPhaseInt[3] = pPhase->CollocationValuesPhaseInt[3] + pPhase->CollocationValuesPhaseInt[2];

		/* A_{0,i} */
		ff  = pPhase->CollocationPointsPhaseInt[0];
		ff1 = pWF->fRING / ff;
		ff2 = ff1 * ff1;
		ff3 = ff1 * ff2;
		ff4 = ff2 * ff2;
		ff0 = (4.0 * pPhase->cL) / ((2.0*pWF->fDAMP)*(2.0*pWF->fDAMP) + (ff - pWF->fRING)*(ff - pWF->fRING));
		gsl_matrix_set(A,0,0,1.0);
		gsl_matrix_set(A,0,1,ff1);
		gsl_matrix_set(A,0,2,ff2);
		gsl_matrix_set(A,0,3,ff3);
		gsl_matrix_set(A,0,4,ff4);
		gsl_vector_set(b,0,pPhase->CollocationValuesPhaseInt[0] - ff0);

		if(debug)
		{
		printf("For row 0: a0 + a1 %.6f + a2 %.6f + a3 %.6f + a4 %.6f = %.6f , ff0 = %.6f, ff = %.6f\n",ff1,ff2,ff3,ff4,pPhase->CollocationValuesPhaseInt[0] - ff0,ff0,ff);
		}

		/* A_{1,i} */
		ff  = pPhase->CollocationPointsPhaseInt[1];
		ff1 = pWF->fRING / ff;
		ff2 = ff1 * ff1;
		ff3 = ff1 * ff2;
		ff4 = ff2 * ff2;
		ff0 = (4 * pPhase->cL) / (4.0*pWF->fDAMP*pWF->fDAMP + (ff - pWF->fRING)*(ff - pWF->fRING));
		gsl_matrix_set(A,1,0,1.0);
		gsl_matrix_set(A,1,1,ff1);
		gsl_matrix_set(A,1,2,ff2);
		gsl_matrix_set(A,1,3,ff3);
		gsl_matrix_set(A,1,4,ff4);
		gsl_vector_set(b,1,pPhase->CollocationValuesPhaseInt[1] - ff0);

		if(debug)
		{
		printf("For row 1: a0 + a1 %.6f + a2 %.6f + a3 %.6f + a4 %.6f = %.6f\n",ff1,ff2,ff3,ff4,pPhase->CollocationValuesPhaseInt[1] - ff0);
		}

		/* A_{2,i} */
		ff  = pPhase->CollocationPointsPhaseInt[2];
		ff1 = pWF->fRING / ff;
		ff2 = ff1 * ff1;
		ff3 = ff1 * ff2;
		ff4 = ff2 * ff2;
		ff0 = (4 * pPhase->cL) / (4.0*pWF->fDAMP*pWF->fDAMP + (ff - pWF->fRING)*(ff - pWF->fRING));
		gsl_matrix_set(A,2,0,1.0);
		gsl_matrix_set(A,2,1,ff1);
		gsl_matrix_set(A,2,2,ff2);
		gsl_matrix_set(A,2,3,ff3);
		gsl_matrix_set(A,2,4,ff4);
		gsl_vector_set(b,2,pPhase->CollocationValuesPhaseInt[2] - ff0);

		if(debug)
		{
		printf("For row 2: a0 + a1 %.6f + a2 %.6f + a3 %.6f + a4 %.6f = %.6f\n",ff1,ff2,ff3,ff4,pPhase->CollocationValuesPhaseInt[2] - ff0);
		}

		/* A_{3,i} */
		ff  = pPhase->CollocationPointsPhaseInt[3];
		ff1 = pWF->fRING / ff;
		ff2 = ff1 * ff1;
		ff3 = ff1 * ff2;
		ff4 = ff2 * ff2;
		ff0 = (4 * pPhase->cL) / (4.0*pWF->fDAMP*pWF->fDAMP + (ff - pWF->fRING)*(ff - pWF->fRING));
		gsl_matrix_set(A,3,0,1.0);
		gsl_matrix_set(A,3,1,ff1);
		gsl_matrix_set(A,3,2,ff2);
		gsl_matrix_set(A,3,3,ff3);
		gsl_matrix_set(A,3,4,ff4);
		gsl_vector_set(b,3,pPhase->CollocationValuesPhaseInt[3] - ff0);

		if(debug)
		{
		printf("For row 3: a0 + a1 %.6f + a2 %.6f + a3 %.6f + a4 %.6f = %.6f\n",ff1,ff2,ff3,ff4,pPhase->CollocationValuesPhaseInt[3] - ff0);
		}

		/* A_{4,i} */
		ff  = pPhase->CollocationPointsPhaseInt[4];
		ff1 = pWF->fRING / ff;
		ff2 = ff1 * ff1;
		ff3 = ff1 * ff2;
		ff4 = ff2 * ff2;
		ff0 = (4 * pPhase->cL) / (4.0*pWF->fDAMP*pWF->fDAMP + (ff - pWF->fRING)*(ff - pWF->fRING));
		gsl_matrix_set(A,4,0,1.0);
		gsl_matrix_set(A,4,1,ff1);
		gsl_matrix_set(A,4,2,ff2);
		gsl_matrix_set(A,4,3,ff3);
		gsl_matrix_set(A,4,4,ff4);
		gsl_vector_set(b,4,pPhase->CollocationValuesPhaseInt[4] - ff0);

		if(debug)
		{
		printf("For row 4: a0 + a1 %.6f + a2 %.6f + a3 %.6f + a4 %.6f = %.6f\n",ff1,ff2,ff3,ff4,pPhase->CollocationValuesPhaseInt[4] - ff0);
		}

		/* We now solve the system A x = b via an LU decomposition */
		gsl_linalg_LU_decomp(A,p,&s);
		gsl_linalg_LU_solve(A,p,b,x);

		/* Set intermediate phenomenological coefficients from solution to A x = b */
		pPhase->b0 = gsl_vector_get(x,0);                                                      // x[0] // Const.
		pPhase->b1 = gsl_vector_get(x,1) * pWF->fRING;                                         // x[1] // f^{-1}
		pPhase->b2 = gsl_vector_get(x,2) * pWF->fRING * pWF->fRING;                            // x[2] // f^{-2}
		pPhase->b3 = gsl_vector_get(x,3) * pWF->fRING * pWF->fRING * pWF->fRING;               // x[3] // f^{-3}
		pPhase->b4 = gsl_vector_get(x,4) * pWF->fRING * pWF->fRING * pWF->fRING * pWF->fRING;  // x[4] // f^{-4}

		if(debug)
		{
		printf("\n");
		printf("Intermediate Collocation Points and Values:\n");
		printf("F1 : %.7f\n",pPhase->CollocationPointsPhaseInt[0]);
		printf("F2 : %.7f\n",pPhase->CollocationPointsPhaseInt[1]);
		printf("F3 : %.7f\n",pPhase->CollocationPointsPhaseInt[2]);
		printf("F4 : %.7f\n",pPhase->CollocationPointsPhaseInt[3]);
		printf("F5 : %.7f\n",pPhase->CollocationPointsPhaseInt[4]);
		printf("\n");
		printf("V's agree with Mathematica...\n");
		printf("V1 : %.7f\n",pPhase->CollocationValuesPhaseInt[0]);
		printf("V2 : %.7f\n",pPhase->CollocationValuesPhaseInt[1]);
		printf("V3 : %.7f\n",pPhase->CollocationValuesPhaseInt[2]);
		printf("V4 : %.7f\n",pPhase->CollocationValuesPhaseInt[3]);
		printf("V5 : %.7f\n",pPhase->CollocationValuesPhaseInt[4]);
		printf("\n");
		printf("g0 : %.7f\n",gsl_vector_get(x,0));
		printf("g1 : %.7f\n",gsl_vector_get(x,1));
		printf("g2 : %.7f\n",gsl_vector_get(x,2));
		printf("g3 : %.7f\n",gsl_vector_get(x,3));
		printf("g4 : %.7f\n",gsl_vector_get(x,4));
		printf("\n");
		printf("b0 : %.7f\n",pPhase->b0);
		printf("b1 : %.7f\n",pPhase->b1);
		printf("b2 : %.7f\n",pPhase->b2);
		printf("b3 : %.7f\n",pPhase->b3);
		printf("b4 : %.7f\n",pPhase->b4);
		printf("\n");
		}

		/* Tidy up */
		gsl_vector_free(b);
		gsl_vector_free(x);
		gsl_matrix_free(A);
		gsl_permutation_free(p);
	}
	else
	{
		XLALPrintError("Error in ComputeIMRPhenomXWaveformVariables: IMRPhenomXIntermediatePhaseVersion is not valid.\n");
	}

	/* Ringdown coefficients */
	REAL8 nonGR_dc1   = XLALSimInspiralWaveformParamsLookupNonGRDC1(LALparams);
	REAL8 nonGR_dc2   = XLALSimInspiralWaveformParamsLookupNonGRDC2(LALparams);
	REAL8 nonGR_dc4   = XLALSimInspiralWaveformParamsLookupNonGRDC4(LALparams);
	REAL8 nonGR_dcl   = XLALSimInspiralWaveformParamsLookupNonGRDCL(LALparams);

	/* Intermediate coefficients */
	REAL8 nonGR_db1   = XLALSimInspiralWaveformParamsLookupNonGRDB1(LALparams);
	REAL8 nonGR_db2   = XLALSimInspiralWaveformParamsLookupNonGRDB2(LALparams);
	REAL8 nonGR_db3   = XLALSimInspiralWaveformParamsLookupNonGRDB3(LALparams);
	REAL8 nonGR_db4   = XLALSimInspiralWaveformParamsLookupNonGRDB4(LALparams);

	/* Inspiral coefficients */
	REAL8 dchi_minus2 = XLALSimInspiralWaveformParamsLookupNonGRDChiMinus2(LALparams);
	REAL8 dchi_minus1 = XLALSimInspiralWaveformParamsLookupNonGRDChiMinus1(LALparams);
	REAL8 dchi0       = XLALSimInspiralWaveformParamsLookupNonGRDChi0(LALparams);
	REAL8 dchi1       = XLALSimInspiralWaveformParamsLookupNonGRDChi1(LALparams);
	REAL8 dchi2       = XLALSimInspiralWaveformParamsLookupNonGRDChi2(LALparams);
	REAL8 dchi3       = XLALSimInspiralWaveformParamsLookupNonGRDChi3(LALparams);
	REAL8 dchi4       = XLALSimInspiralWaveformParamsLookupNonGRDChi4(LALparams);
	REAL8 dchi5       = XLALSimInspiralWaveformParamsLookupNonGRDChi5(LALparams);
	REAL8 dchi5L      = XLALSimInspiralWaveformParamsLookupNonGRDChi5L(LALparams);
	REAL8 dchi6       = XLALSimInspiralWaveformParamsLookupNonGRDChi6(LALparams);
	REAL8 dchi6L      = XLALSimInspiralWaveformParamsLookupNonGRDChi6L(LALparams);
	REAL8 dchi7       = XLALSimInspiralWaveformParamsLookupNonGRDChi7(LALparams);

	/* Can include these terms in the future as desired... */
	REAL8 dchi8       = 0.0;
	REAL8 dchi8L      = 0.0;
	REAL8 dchi9       = 0.0;
	REAL8 dchi9L      = 0.0;

	/* ~~~~ RINGDOWN ~~~~ */
	pPhase->cLGR  = pPhase->cL; // Store GR value for reference
	pPhase->c1   *= (1.0 + nonGR_dc1);
	pPhase->c2   *= (1.0 + nonGR_dc2);
	pPhase->c4   *= (1.0 + nonGR_dc4);
	pPhase->cL   *= (1.0 + nonGR_dcl);

	/* Set pre-cached variables */
	pPhase->c4ov3   = pPhase->c4 / 3.0;
	pPhase->cLovfda = pPhase->cL / pWF->fDAMP;

	/* Apply NR tuning for precessing cases (500) */
	pPhase->b1 = pPhase->b1  +  ( pWF->PNR_DEV_PARAMETER * pWF->ZETA2 );
	pPhase->b4 = pPhase->b4  +  ( pWF->PNR_DEV_PARAMETER * pWF->ZETA1 );

	/* ~~~~ INTERMEDIATE ~~~~ */
	if(pWF->IMRPhenomXIntermediatePhaseVersion == 104)
	{
		pPhase->b1 *= (1.0 + nonGR_db1);
		pPhase->b2 *= (1.0 + nonGR_db2);
		pPhase->b4 *= (1.0 + nonGR_db4);
	}
	else if(pWF->IMRPhenomXIntermediatePhaseVersion == 105)
	{
		pPhase->b1 *= (1.0 + nonGR_db1);
		pPhase->b2 *= (1.0 + nonGR_db2);
		pPhase->b3 *= (1.0 + nonGR_db3);
		pPhase->b4 *= (1.0 + nonGR_db4);
	}
	else
	{
		XLALPrintError("Error in ComputeIMRPhenomXWaveformVariables: IMRPhenomXIntermediatePhaseVersion is not valid.\n");
	}

	/* ~~~~ INSPIRAL ~~~~ */
	/* Initialize -1PN coefficient*/
	pPhase->phi_minus2   = 0.0;
	pPhase->dphi_minus2  = 0.0;

	pPhase->phi_minus1   = 0.0;
	pPhase->dphi_minus1  = 0.0;

	/*
		If tgr_parameterization = 1, deform complete PN coefficient. This is an FTA-like parameterization.
		If tgr_parameterization = 0, only deform non-spinning coefficient. This is the original TIGER-like implementation.
	*/
	int tgr_parameterization = 0;
	tgr_parameterization     = XLALSimInspiralWaveformParamsLookupNonGRParameterization(LALparams);

	if(tgr_parameterization == 1)
	{
			/* -1.0 PN: This vanishes in GR, so is parameterized as an absolute deviation */
			pPhase->phi_minus2 = dchi_minus2 / powers_of_lalpi.two_thirds;

			/* -0.5 PN: This vanishes in GR, so is parameterized as an absolute deviation */
			pPhase->phi_minus1 = dchi_minus1 / powers_of_lalpi.one_third;

			/* 0.0 PN */
			pPhase->phi0       = (phi0NS + phi0S)*(1.0 + dchi0);

			/* 0.5 PN: This vanishes in GR, so is parameterized as an absolute deviation */
			pPhase->phi1       = dchi1 * powers_of_lalpi.one_third;

			/* 1.0 PN */
			pPhase->phi2       = (phi2NS + phi2S)*(1.0 + dchi2);

			/* 1.5 PN */
			pPhase->phi3       = (phi3NS + phi3S)*(1.0 + dchi3);

			/* 2.0 PN */
			pPhase->phi4       = (phi4NS + phi4S)*(1.0 + dchi4);

			/* 2.5 PN */
			pPhase->phi5       = (phi5NS + phi5S)*(1.0 + dchi5);

			/* 2.5 PN, Log Terms */
			pPhase->phi5L      = (phi5LNS + phi5LS)*(1.0 + dchi5L);

			/* 3.0 PN */
			pPhase->phi6       = (phi6NS + phi6S)*(1.0 + dchi6);

			/* 3.0 PN, Log Term */
			pPhase->phi6L      = (phi6LNS + phi6LS)*(1.0 + dchi6L);

			/* 3.5PN */
			pPhase->phi7       = (phi7NS + phi7S)*(1.0 + dchi7);

			/* 4.0PN */
			pPhase->phi8       = (phi8NS + phi8S)*(1.0 + dchi8);

			/* 4.0 PN, Log Terms */
			pPhase->phi8L      = (phi8LNS + phi8LS)*(1.0 + dchi8L);

			/* 4.0 PN */
			pPhase->phi9       = (phi9NS + phi9S)*(1.0 + dchi9);

			/* 4.0 PN, Log Terms */
			pPhase->phi9L      = (phi9LNS + phi9LS)*(1.0 + dchi9L);
	}
	else if(tgr_parameterization == 0)
	{
			/* -1.0 PN: This vanishes in GR, so is parameterized as an absolute deviation */
			pPhase->phi_minus2 = dchi_minus2 / powers_of_lalpi.two_thirds;

			/* -0.5 PN: This vanishes in GR, so is parameterized as an absolute deviation */
			pPhase->phi_minus1 = dchi_minus1 / powers_of_lalpi.one_third;

			/* 0.0 PN */
			pPhase->phi0       = phi0NS*(1.0 + dchi0) + phi0S;

			/* 0.5 PN: This vanishes in GR, so is parameterized as an absolute deviation */
			pPhase->phi1       = dchi1 * powers_of_lalpi.one_third;

			/* 1.0 PN */
			pPhase->phi2       = phi2NS*(1.0 + dchi2) + phi2S;

			/* 1.5 PN */
			pPhase->phi3       = phi3NS*(1.0 + dchi3)+ phi3S;

			/* 2.0 PN */
			pPhase->phi4       = phi4NS*(1.0 + dchi4) + phi4S;

			/* 2.5 PN */
			pPhase->phi5       = phi5NS*(1.0 + dchi5) + phi5S;

			/* 2.5 PN, Log Terms */
			pPhase->phi5L      = phi5LNS*(1.0 + dchi5L) + phi5LS;

			/* 3.0 PN */
			pPhase->phi6       = phi6NS*(1.0 + dchi6) + phi6S;

			/* 3.0 PN, Log Term */
			pPhase->phi6L      = phi6LNS*(1.0 + dchi6L) + phi6LS;

			/* 3.5PN */
			pPhase->phi7       = phi7NS*(1.0 + dchi7) + phi7S;

			/* 4.0PN */
			pPhase->phi8       = phi8NS*(1.0 + dchi8) + phi8S;

			/* 4.0 PN, Log Terms */
			pPhase->phi8L      = phi8LNS*(1.0 + dchi8L) + phi8LS;

			/* 4.0 PN */
			pPhase->phi9       = phi9NS*(1.0 + dchi9) + phi9S;

			/* 4.0 PN, Log Terms */
			pPhase->phi9L      = phi9LNS*(1.0 + dchi9L) + phi9LS;
	}
	else
	{
			XLALPrintError("Error in IMRPhenomXGetPhaseCoefficients: TGR Parameterizataion is not valid.\n");
	}

	/* Recalculate phase derivatives including TGR corrections */
	pPhase->dphi_minus2 = +(7.0 / 5.0) * pPhase->phi_minus2;
	pPhase->dphi_minus1 = +(6.0 / 5.0) * pPhase->phi_minus1;
	pPhase->dphi0       = +(5.0 / 5.0) * pPhase->phi0;
	pPhase->dphi1       = +(4.0 / 5.0) * pPhase->phi1;
	pPhase->dphi2       = +(3.0 / 5.0) * pPhase->phi2;
	pPhase->dphi3       = +(2.0 / 5.0) * pPhase->phi3;
	pPhase->dphi4       = +(1.0 / 5.0) * pPhase->phi4;
	pPhase->dphi5       = -(3.0 / 5.0) * pPhase->phi5L;
	pPhase->dphi6       = -(1.0 / 5.0) * pPhase->phi6 - (3.0 / 5.0) * pPhase->phi6L;
	pPhase->dphi6L      = -(1.0 / 5.0) * pPhase->phi6L;
	pPhase->dphi7       = -(2.0 / 5.0) * pPhase->phi7;
	pPhase->dphi8       = -(3.0 / 5.0) * pPhase->phi8 - (3.0 / 5.0) * pPhase->phi8L;
	pPhase->dphi8L      = -(3.0 / 5.0) * pPhase->phi8L;
	pPhase->dphi9       = -(4.0 / 5.0) * pPhase->phi9 - (3.0 / 5.0) * pPhase->phi9L;
	pPhase->dphi9L      = -(3.0 / 5.0) * pPhase->phi9L;

	/* Initialize connection coefficients */
	pPhase->C1Int = 0;
	pPhase->C2Int = 0;
	pPhase->C1MRD = 0;
	pPhase->C2MRD = 0;

	/* END OF ROUTINE; RETURN STRUCT */
	return XLAL_SUCCESS;
}


void IMRPhenomX_Phase_22_ConnectionCoefficients(IMRPhenomXWaveformStruct *pWF, IMRPhenomXPhaseCoefficients *pPhase)
{
	int debug = pWF->debug;

	double fIns = pPhase->fPhaseMatchIN;
	double fInt = pPhase->fPhaseMatchIM;

	/*
			Assume an ansatz of the form:

			phi_Inspiral (f) = phi_Intermediate (f) + C1 + C2 * f

			where transition frequency is fIns

			phi_Inspiral (fIns) = phi_Intermediate (fIns) + C1 + C2 * fIns
			phi_Inspiral'(fIns) = phi_Intermediate'(fIns) + C2

			Solving for C1 and C2

			C2 = phi_Inspiral'(fIns) - phi_Intermediate'(fIns)
			C1 = phi_Inspiral (fIns) - phi_Intermediate (fIns) - C2 * fIns

	*/
	IMRPhenomX_UsefulPowers powers_of_fIns;
	IMRPhenomX_Initialize_Powers(&powers_of_fIns,fIns);

	double DPhiIns = IMRPhenomX_Inspiral_Phase_22_Ansatz(fIns,&powers_of_fIns,pPhase);
	double DPhiInt = IMRPhenomX_Intermediate_Phase_22_Ansatz(fIns,&powers_of_fIns,pWF,pPhase);

	pPhase->C2Int  = DPhiIns - DPhiInt;

	double phiIN = IMRPhenomX_Inspiral_Phase_22_AnsatzInt(fIns,&powers_of_fIns,pPhase);
	double phiIM = IMRPhenomX_Intermediate_Phase_22_AnsatzInt(fIns,&powers_of_fIns,pWF,pPhase);

	if(debug)
	{
	printf("\n");
	printf("dphiIM = %.6f and dphiIN = %.6f\n",DPhiInt,DPhiIns);
	printf("phiIN(fIns)  : %.7f\n",phiIN);
	printf("phiIM(fIns)  : %.7f\n",phiIM);
	printf("fIns         : %.7f\n",fIns);
	printf("C2           : %.7f\n",pPhase->C2Int);
	printf("\n");
	}

	pPhase->C1Int = phiIN - phiIM - (pPhase->C2Int * fIns);

	/*
			Assume an ansatz of the form:

			phi_Intermediate (f)    = phi_Ringdown (f) + C1 + C2 * f

			where transition frequency is fIM

			phi_Intermediate (fIM) = phi_Ringdown (fRD) + C1 + C2 * fIM
			phi_Intermediate'(fIM) = phi_Ringdown'(fRD) + C2

			Solving for C1 and C2

			C2 = phi_Inspiral'(fIM) - phi_Intermediate'(fIM)
			C1 = phi_Inspiral (fIM) - phi_Intermediate (fIM) - C2 * fIM

	*/
	IMRPhenomX_UsefulPowers powers_of_fInt;
	IMRPhenomX_Initialize_Powers(&powers_of_fInt,fInt);

	double phiIMC         = IMRPhenomX_Intermediate_Phase_22_AnsatzInt(fInt,&powers_of_fInt,pWF,pPhase) + pPhase->C1Int + pPhase->C2Int*fInt;
	double phiRD          = IMRPhenomX_Ringdown_Phase_22_AnsatzInt(fInt,&powers_of_fInt,pWF,pPhase);
	double DPhiIntC       = IMRPhenomX_Intermediate_Phase_22_Ansatz(fInt,&powers_of_fInt,pWF,pPhase) + pPhase->C2Int;
	double DPhiRD         = IMRPhenomX_Ringdown_Phase_22_Ansatz(fInt,&powers_of_fInt,pWF,pPhase);

	pPhase->C2MRD = DPhiIntC - DPhiRD;
	pPhase->C1MRD = phiIMC - phiRD - pPhase->C2MRD*fInt;

	if(debug)
	{
	printf("\n");
	printf("phiIMC(fInt) : %.7f\n",phiIMC);
	printf("phiRD(fInt)  : %.7f\n",phiRD);
	printf("fInt         : %.7f\n",fInt);
	printf("C2           : %.7f\n",pPhase->C2Int);
	printf("\n");
	}

	if(debug)
	{
	printf("dphiIM = %.6f and dphiRD = %.6f\n",DPhiIntC,DPhiRD);
	printf("\nContinuity Coefficients\n");
	printf("C1Int : %.6f\n",pPhase->C1Int);
	printf("C2Int : %.6f\n",pPhase->C2Int);
	printf("C1MRD : %.6f\n",pPhase->C1MRD);
	printf("C2MRD : %.6f\n",pPhase->C2MRD);
	}

  return;
}

double IMRPhenomX_TimeShift_22(IMRPhenomXPhaseCoefficients *pPhase, IMRPhenomXWaveformStruct *pWF){

    REAL8 linb, tshift, dphi22Ref,frefFit;

    // here we align the model to the hybrids, for which psi4 peaks 500M before the end of the waveform
    // linb is a parameter-space fit of dphi22(fring22-fdamp22), evaluated on the calibration dataset
    linb    = XLALSimIMRPhenomXLinb(pWF->eta, pWF->STotR, pWF->dchi, pWF->delta);
    frefFit = pWF->fRING-pWF->fDAMP;
    IMRPhenomX_UsefulPowers powers_of_frefFit;
    IMRPhenomX_Initialize_Powers(&powers_of_frefFit,frefFit);
    dphi22Ref = 1.0 / (pWF->eta)*IMRPhenomX_dPhase_22(frefFit, &powers_of_frefFit, pPhase, pWF);
    // here we correct the time-alignment of the waveform by first aligning the peak of psi4, and then adding a correction to align the peak of strain instead
    REAL8 psi4tostrain=XLALSimIMRPhenomXPsi4ToStrain(pWF->eta, pWF->STotR, pWF->dchi);
    tshift = linb-dphi22Ref -2.*LAL_PI*(500+psi4tostrain);

	// Apply PNR deviation (500)
	tshift = tshift + ( pWF->PNR_DEV_PARAMETER * pWF->NU0 );

    //phX phase will read phi22=1/eta*IMRPhenomX_Phase_22+tshift f, modulo a residual phase-shift
    return(tshift);

}



/*
 * ********** ********** ********** ********** ********** ********** ********** ********** ********** **********
 * This function computes the IMRPhenomX phase given a phase coefficients struct pPhase.
 * ********** ********** ********** ********** ********** ********** ********** ********** ********** **********
 */
double IMRPhenomX_Phase_22(double ff, IMRPhenomX_UsefulPowers *powers_of_f, IMRPhenomXPhaseCoefficients *pPhase, IMRPhenomXWaveformStruct *pWF)
{
	// Inspiral region, f < fPhaseInsMax
  if (!IMRPhenomX_StepFuncBool(ff, pPhase->fPhaseMatchIN))
  {
	  double PhiIns = IMRPhenomX_Inspiral_Phase_22_AnsatzInt(ff, powers_of_f, pPhase);
	  return PhiIns;
  }

	// Ringdown region, f > fPhaseIntMax
  if (IMRPhenomX_StepFuncBool(ff, pPhase->fPhaseMatchIM))
  {
	  double PhiMRD = IMRPhenomX_Ringdown_Phase_22_AnsatzInt(ff, powers_of_f, pWF, pPhase)
                            + pPhase->C1MRD + (pPhase->C2MRD * ff);
	  return PhiMRD;
  }

  //	Intermediate region, fPhaseInsMax < f < fPhaseIntMax
  double PhiInt = IMRPhenomX_Intermediate_Phase_22_AnsatzInt(ff, powers_of_f, pWF, pPhase)
                            + pPhase->C1Int + (pPhase->C2Int * ff);

  return PhiInt;
}

/*
 * ********** ********** ********** ********** ********** ********** ********** ********** ********** **********
 * This function computes the IMRPhenomX phase derivative given a phase coefficients struct pPhase.
 * ********** ********** ********** ********** ********** ********** ********** ********** ********** **********
 */
double IMRPhenomX_dPhase_22(double ff, IMRPhenomX_UsefulPowers *powers_of_f, IMRPhenomXPhaseCoefficients *pPhase, IMRPhenomXWaveformStruct *pWF)
{
	// Inspiral region, f < fPhaseInsMax
  if (!IMRPhenomX_StepFuncBool(ff, pPhase->fPhaseMatchIN))
  {
	  double dPhiIns = IMRPhenomX_Inspiral_Phase_22_Ansatz(ff, powers_of_f, pPhase);
	  return dPhiIns;
  }

	// Ringdown region, f > fPhaseIntMax
  if (IMRPhenomX_StepFuncBool(ff, pPhase->fPhaseMatchIM))
  {
	  double dPhiMRD = IMRPhenomX_Ringdown_Phase_22_Ansatz(ff, powers_of_f, pWF, pPhase)
                            + (pPhase->C2MRD);
	  return dPhiMRD;
  }

  //	Intermediate region, fPhaseInsMax < f < fPhaseIntMax
  double dPhiInt = IMRPhenomX_Intermediate_Phase_22_Ansatz(ff, powers_of_f, pWF, pPhase)
                            + (pPhase->C2Int);

  return dPhiInt;
}

/*
 * ********** ********** ********** ********** ********** ********** ********** ********** ********** **********
 * This function computes the IMRPhenomX amplitude given an amplitude coefficients struct pAmp.
 * ********** ********** ********** ********** ********** ********** ********** ********** ********** **********
 */
double IMRPhenomX_Amplitude_22(double ff, IMRPhenomX_UsefulPowers *powers_of_f, IMRPhenomXAmpCoefficients *pAmp, IMRPhenomXWaveformStruct *pWF) {
  //double f_seven_sixths = powers_of_f->seven_sixths;
  //double AmpPreFac      = pWF->ampNorm / f_seven_sixths; // Use this if we want return normalized amplitudes
	double AmpPreFac      = pWF->ampNorm / powers_of_f->seven_sixths;

  // Inspiral region, f < fAmpMatchIN
  if (!IMRPhenomX_StepFuncBool(ff, pAmp->fAmpMatchIN))
  {
    double AmpIns =  AmpPreFac * IMRPhenomX_Inspiral_Amp_22_Ansatz(ff, powers_of_f, pWF, pAmp);
	  return AmpIns;
  }

	// Ringdown region, f > fAmpRDMin
  if (IMRPhenomX_StepFuncBool(ff, pAmp->fAmpRDMin))
  {
    double AmpMRD = AmpPreFac * IMRPhenomX_Ringdown_Amp_22_Ansatz(ff, pWF, pAmp);
    return AmpMRD;
  }

  // Intermediate region, fAmpMatchIN < f < fAmpRDMin
  double AmpInt = AmpPreFac * IMRPhenomX_Intermediate_Amp_22_Ansatz(ff, powers_of_f, pWF, pAmp);
  return AmpInt;
}


/* Function to check if the input mode array contains unsupported modes */
INT4 check_input_mode_array(LALDict *lalParams)
{
	UINT4 flagTrue = 0;

  if(lalParams == NULL) return XLAL_SUCCESS;

  LALValue *ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(lalParams);

  if(ModeArray!=NULL)
  {
    INT4 larray[5] = {2, 2, 3, 3, 4};
    INT4 marray[5] = {2, 1, 3, 2, 4};

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

/* Function to compute full model phase. This function is designed to be used in in initialization routines, and not for evaluating the phase at many frequencies. */
INT4 IMRPhenomX_FullPhase_22(double *phase, double *dphase, double Mf, IMRPhenomXPhaseCoefficients *pPhase, IMRPhenomXWaveformStruct *pWF){


    /*
    Function to compute full XAS phase at a single frequency point.
	See IMRPhenomXASGenerateFD for reference.
    */

	/*--(*)--(*)--(*)--(*)--(*)--(*)--(*)--(*)--(*)--*/
	/*            Define useful powers               */
	/*--(*)--(*)--(*)--(*)--(*)--(*)--(*)--(*)--(*)--*/

	// Status indicator
	INT4 status;

	// Get useful powers of Mf
	IMRPhenomX_UsefulPowers powers_of_Mf;
    status = IMRPhenomX_Initialize_Powers(&powers_of_Mf, Mf);
	XLAL_CHECK(XLAL_SUCCESS == status, status, "IMRPhenomX_Initialize_Powers failed for Mf.\n");

	/* Initialize a struct containing useful powers of Mf at fRef */
	IMRPhenomX_UsefulPowers powers_of_MfRef;
	status = IMRPhenomX_Initialize_Powers(&powers_of_MfRef,pWF->MfRef);
	XLAL_CHECK(XLAL_SUCCESS == status, status, "IMRPhenomX_Initialize_Powers failed for MfRef.\n");

	/*--(*)--(*)--(*)--(*)--(*)--(*)--(*)--(*)--(*)--*/
	/*           Define needed constants             */
	/*--(*)--(*)--(*)--(*)--(*)--(*)--(*)--(*)--(*)--*/

	/* 1/eta is used to re-scale the pre-phase quantity */
	REAL8 inveta    = (1.0 / pWF->eta);

	/* We keep this phase shift to ease comparison with
	original phase routines */
	REAL8 lina = 0;

	/* Get phase connection coefficients
	and store them to pWF. This step is to make
	sure that teh coefficients are up-to-date */
	IMRPhenomX_Phase_22_ConnectionCoefficients(pWF,pPhase);

	/* Compute the timeshift that PhenomXAS uses to align waveforms
	with the hybrids used to make their model */
	double linb=IMRPhenomX_TimeShift_22(pPhase, pWF);

	/* Calculate phase at reference frequency: phifRef = 2.0*phi0 + LAL_PI_4 + PhenomXPhase(fRef) */
	double phifRef = -(inveta * IMRPhenomX_Phase_22(pWF->MfRef, &powers_of_MfRef, pPhase, pWF) + linb*pWF->MfRef + lina) + 2.0*pWF->phi0 + LAL_PI_4;

	/* ~~+~~+~~+~~+~~+~~+~~+~~+~~+~~+~~+~~+~~+~~+~~+~~+~~+~~+~~+~~+
	Note that we do not store the value of phifRef to pWF as is done in
	IMRPhenomXASGenerateFD. We choose to not do so in order to avoid
	potential confusion (e.g. if this function is called within a
	workflow that assumes the value defined in IMRPhenomXASGenerateFD).
	Note that this concern may not be valid.
	~~+~~+~~+~~+~~+~~+~~+~~+~~+~~+~~+~~+~~+~~+~~+~~+~~+~~+~~+~~+~~+ */

	/*--(*)--(*)--(*)--(*)--(*)--(*)--(*)--(*)--(*)--*/
	/*        Compute the full model phase           */
	/*--(*)--(*)--(*)--(*)--(*)--(*)--(*)--(*)--(*)--*/

	/* Use previously made function to compute what we call
	here the pre-phase, becuase it's not actually a phase! */
	double pre_phase = IMRPhenomX_Phase_22(Mf,&powers_of_Mf,pPhase,pWF);

	/* Given the pre-phase, we need to scale and shift according to the
	XAS construction */
	*phase   = pre_phase * inveta;
	*phase  += linb*Mf + lina + phifRef;

	/* Repeat the excercise above for the phase derivative:
	"dphase" is (d/df)phase at Mf */
	double pre_dphase = IMRPhenomX_dPhase_22(Mf,&powers_of_Mf,pPhase,pWF);
    *dphase  = pre_dphase * inveta;
    *dphase += linb;

	//
	return status;
}


NRTidal_version_type IMRPhenomX_SetTidalVersion(LALDict *lalParams){

    int tidal_version=XLALSimInspiralWaveformParamsLookupPhenomXTidalFlag(lalParams);

    NRTidal_version_type version;

    switch(tidal_version){
          case 0:
              version=NoNRT_V;
              break;
          case 1:
              version=NRTidal_V;
              break;
          case 2:
              version=NRTidalv2_V;
              break;
		  case 3:
              version=NRTidalv3_V;
			  break;
          default:
              {
                  XLAL_ERROR(XLAL_EINVAL, "Error: Tidal version not recognized. Only NRTidal, NRTidalv2, NRTidalv3, and NoNRT are allowed, and NRTidal is not implemented completely in IMRPhenomX*.\n");
              }
    }
    return(version);
}


void IMRPhenomXGetTidalPhaseCoefficients(
    IMRPhenomXWaveformStruct *pWF,
    IMRPhenomXPhaseCoefficients *pPhase,
    NRTidal_version_type NRTidal_version){


    REAL8 quadparam1 = pWF->quadparam1;
    REAL8 quadparam2 = pWF->quadparam2;
    /* declare HO 3.5PN spin-spin and spin-cubed terms added separately in Pv2_NRTidalv2 */
    REAL8 SS_3p5PN = 0., SSS_3p5PN = 0.;

    /* New variables needed for the NRTidalv2 model */
    REAL8 X_A = pWF->m1; // Already scaled by Mtot
    REAL8 X_B = pWF->m2; // Ibid.

    REAL8 chi1L_sq = pWF->chi1L2;
    REAL8 chi2L_sq = pWF->chi2L2;

     pPhase->c2PN_tidal=((XLALSimInspiralTaylorF2Phasing_4PNQM2SOCoeff(X_A) + XLALSimInspiralTaylorF2Phasing_4PNQM2SCoeff(X_A)) * (quadparam1 - 1.) * chi1L_sq
       + (XLALSimInspiralTaylorF2Phasing_4PNQM2SOCoeff(X_B) + XLALSimInspiralTaylorF2Phasing_4PNQM2SCoeff(X_B)) * (quadparam2 - 1.) * chi2L_sq);
     pPhase->c3PN_tidal=(XLALSimInspiralTaylorF2Phasing_6PNQM2SCoeff(X_A) * (quadparam1 - 1.) * chi1L_sq
     + XLALSimInspiralTaylorF2Phasing_6PNQM2SCoeff(X_B) * (quadparam2 - 1.) * chi2L_sq);

      if (NRTidal_version == NRTidalv2_V || NRTidal_version == NRTidalv3_V) {
     /* Get the PN SS-tail and SSS terms */
       XLALSimInspiralGetHOSpinTerms(&SS_3p5PN, &SSS_3p5PN, X_A, X_B, pWF->chi1L, pWF->chi2L, quadparam1, quadparam2);
       pPhase->c3p5PN_tidal=(SS_3p5PN + SSS_3p5PN);
      }

}

REAL8 IMRPhenomX_TidalPhase(IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXWaveformStruct *pWF, IMRPhenomXPhaseCoefficients *pPhase, NRTidal_version_type NRTidal_version){

    REAL8 X_A = pWF->m1; // Already scaled by Mtot
    REAL8 X_B = pWF->m2;

    REAL8 Mf=powers_of_Mf->itself;

    REAL8 phaseTidal=0.;

    REAL8 pfaN=3./(128.*X_A*X_B);

    REAL8 c2pn=pPhase->c2PN_tidal, c3pn=pPhase->c3PN_tidal, c3p5pn=pPhase->c3p5PN_tidal;

    /* 2PN terms */
    phaseTidal += pfaN * c2pn* powers_of_lalpi.m_one_third * powers_of_Mf->m_one_third;

    /* 3PN terms */
    phaseTidal += pfaN * c3pn* powers_of_lalpi.one_third * powers_of_Mf->one_third;

    if (NRTidal_version == NRTidalv2_V)
    {
        REAL8 NRTidalv2_coeffs[9];

        int errcode;
        errcode = XLALSimNRTunedTidesSetFDTidalPhase_v2_Coeffs(NRTidalv2_coeffs);
        XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "Setting NRTidalv2 coefficients failed.\n");

        const REAL8 cNewt   = NRTidalv2_coeffs[0];
        const REAL8 n1      = NRTidalv2_coeffs[1];
        const REAL8 n3over2 = NRTidalv2_coeffs[2];
        const REAL8 n2      = NRTidalv2_coeffs[3];
        const REAL8 n5over2 = NRTidalv2_coeffs[4];
        const REAL8 n3      = NRTidalv2_coeffs[5];
        const REAL8 d1      = NRTidalv2_coeffs[6];
        const REAL8 d3over2 = NRTidalv2_coeffs[7];
        const REAL8 d2      = NRTidalv2_coeffs[8];

        REAL8 kappa2T = pWF->kappa2T;

        REAL8 NRphase=-((cNewt*kappa2T*powers_of_Mf->five_thirds*powers_of_lalpi.five_thirds*(1 + powers_of_Mf->two_thirds*n1*powers_of_lalpi.two_thirds + Mf*n3over2*LAL_PI + powers_of_Mf->four_thirds*n2*powers_of_lalpi.four_thirds + powers_of_Mf->five_thirds*n5over2*powers_of_lalpi.five_thirds + pow(Mf,2)*n3*powers_of_lalpi.two))/((1 + d1*powers_of_Mf->two_thirds*powers_of_lalpi.two_thirds + d3over2*Mf*LAL_PI + d2*powers_of_Mf->four_thirds*powers_of_lalpi.four_thirds)*X_A*X_B));

        phaseTidal+=NRphase;

        /* Get the PN SS-tail and SSS terms */
        phaseTidal += pfaN * c3p5pn * powers_of_lalpi.two_thirds * powers_of_Mf->two_thirds;

    }
	else if (NRTidal_version == NRTidalv3_V)
    {
         /* local variable declaration */
		REAL8 Xa = X_A;

		REAL8 PN_coeffs[10];
		XLALSimNRTunedTidesSetFDTidalPhase_PN_Coeffs(PN_coeffs, Xa);

		REAL8 lambda1 = pWF->lambda1;
		REAL8 lambda2 = pWF->lambda2;

		REAL8 mtot = pWF->Mtot;

		REAL8 NRTidalv3_coeffs[20];
    	XLALSimNRTunedTidesSetFDTidalPhase_v3_Coeffs(NRTidalv3_coeffs, Xa, mtot, lambda1, lambda2, PN_coeffs);

		REAL8 s1 = NRTidalv3_coeffs[0];
		REAL8 s2 = NRTidalv3_coeffs[1];
		REAL8 s3 = NRTidalv3_coeffs[2];


		REAL8 s2s3 = s2*s3;
		REAL8 s2Mf = -s2*Mf*LAL_PI*2.0;
		REAL8 exps2s3 = cosh(s2s3) + sinh(s2s3);
		REAL8 exps2Mf = cosh(s2Mf) + sinh(s2Mf);

		/* Rewriting Eq. (27) of arXiv:2311.07456 */
  		REAL8 dynk2barfunc = 1.0 + ((s1) - 1)*(1.0/(1.0 + exps2Mf*exps2s3)) - ((s1-1.0)/(1.0 + exps2s3)) - 2.0*(Mf*LAL_PI)*((s1) - 1)*s2*exps2s3/((1.0 + exps2s3)*(1.0 + exps2s3));

		REAL8 kappaA = NRTidalv3_coeffs[4];
		REAL8 kappaB = NRTidalv3_coeffs[5];

		REAL8 dynkappaA = kappaA*dynk2barfunc;
		REAL8 dynkappaB = kappaB*dynk2barfunc;

		/* Pade Coefficients, Table II of arXiv:2311.07456 */
		REAL8 n_5over2A = NRTidalv3_coeffs[6];
		REAL8 n_3A = NRTidalv3_coeffs[7];
		REAL8 d_1A = NRTidalv3_coeffs[8];

		REAL8 n_5over2B = NRTidalv3_coeffs[9];
		REAL8 n_3B = NRTidalv3_coeffs[10];
		REAL8 d_1B = NRTidalv3_coeffs[11];

		/* 7.5PN Coefficients */
		REAL8 c_NewtA = PN_coeffs[0];
		REAL8 c_1A = PN_coeffs[1];
		REAL8 c_3over2A = PN_coeffs[2];
		REAL8 c_2A = PN_coeffs[3];
		REAL8 c_5over2A = PN_coeffs[4];

		REAL8 c_NewtB = PN_coeffs[5];
		REAL8 c_1B = PN_coeffs[6];
		REAL8 c_3over2B = PN_coeffs[7];
		REAL8 c_2B = PN_coeffs[8];
		REAL8 c_5over2B = PN_coeffs[9];

		/* Pade Coefficients constrained with PN, see Eq. (33) in arXiv:2311.07456 */
		REAL8 n_1A = NRTidalv3_coeffs[12];
		REAL8 n_3over2A = NRTidalv3_coeffs[13];
		REAL8 n_2A = NRTidalv3_coeffs[14];
		REAL8 d_3over2A = NRTidalv3_coeffs[15];

		REAL8 n_1B = NRTidalv3_coeffs[16];
		REAL8 n_3over2B = NRTidalv3_coeffs[17];
		REAL8 n_2B = NRTidalv3_coeffs[18];
		REAL8 d_3over2B = NRTidalv3_coeffs[19];

		/* Rewriting Eq. (30) and (32) in arXiv:2311.07456 in terms of Mf instead of x */
        REAL8 NRphasetermA=-((c_NewtA*dynkappaA*powers_of_Mf->five_thirds*powers_of_lalpi.five_thirds*(1 + powers_of_Mf->two_thirds*n_1A*powers_of_lalpi.two_thirds + Mf*n_3over2A*LAL_PI + powers_of_Mf->four_thirds*n_2A*powers_of_lalpi.four_thirds + powers_of_Mf->five_thirds*n_5over2A*powers_of_lalpi.five_thirds + pow(Mf,2)*n_3A*powers_of_lalpi.two))/((1 + d_1A*powers_of_Mf->two_thirds*powers_of_lalpi.two_thirds + d_3over2A*Mf*LAL_PI)));
        REAL8 NRphasetermB=-((c_NewtB*dynkappaB*powers_of_Mf->five_thirds*powers_of_lalpi.five_thirds*(1 + powers_of_Mf->two_thirds*n_1B*powers_of_lalpi.two_thirds + Mf*n_3over2B*LAL_PI + powers_of_Mf->four_thirds*n_2B*powers_of_lalpi.four_thirds + powers_of_Mf->five_thirds*n_5over2B*powers_of_lalpi.five_thirds + pow(Mf,2)*n_3B*powers_of_lalpi.two))/((1 + d_1B*powers_of_Mf->two_thirds*powers_of_lalpi.two_thirds + d_3over2B*Mf*LAL_PI)));

		REAL8 NRphaseNRT = NRphasetermA + NRphasetermB;

		REAL8 PNtidalphaseA = -((c_NewtA*kappaA*powers_of_Mf->five_thirds*powers_of_lalpi.five_thirds*(1 + powers_of_Mf->two_thirds*c_1A*powers_of_lalpi.two_thirds + Mf*c_3over2A*LAL_PI + powers_of_Mf->four_thirds*c_2A*powers_of_lalpi.four_thirds + powers_of_Mf->five_thirds*c_5over2A*powers_of_lalpi.five_thirds )));
		REAL8 PNtidalphaseB = -((c_NewtB*kappaB*powers_of_Mf->five_thirds*powers_of_lalpi.five_thirds*(1 + powers_of_Mf->two_thirds*c_1B*powers_of_lalpi.two_thirds + Mf*c_3over2B*LAL_PI + powers_of_Mf->four_thirds*c_2B*powers_of_lalpi.four_thirds + powers_of_Mf->five_thirds*c_5over2B*powers_of_lalpi.five_thirds )));

		REAL8 PNtidalphase = PNtidalphaseA + PNtidalphaseB;

		REAL8 Mfmerger = XLALSimNRTunedTidesMergerFrequency_v3(pWF->Mtot, pWF->lambda1, pWF->lambda2, pWF->q, pWF->chi1L, pWF->chi2L) * (pWF->Mtot * LAL_MTSUN_SI) ;

		REAL8 Mftaperstart = 1.15*Mfmerger;
		REAL8 Mftaperend = 1.35*Mfmerger;

		/* Eq. (46) in arXiv:2311.07456 */
		REAL8 exp_arg;
		REAL8 exp_fn;
		REAL8 plancktaperfn;

		if (Mf <= Mftaperstart){
		plancktaperfn = 0.0;
		}
		else if (Mf >= Mftaperend){
		plancktaperfn = 1.0;
		}
		else{
		exp_arg = (Mftaperend - Mftaperstart)/(Mf - Mftaperstart) + (Mftaperend - Mftaperstart)/(Mf - Mftaperend);
		exp_fn = cosh(exp_arg) + sinh(exp_arg);
		plancktaperfn = 1/(exp_fn + 1.0);
		}

		/* Eq. (45) in arXiv:2311.07456 */
		REAL8 NRphase = NRphaseNRT*(1-  plancktaperfn) + PNtidalphase*plancktaperfn;

		phaseTidal+=NRphase;

        /* Get the PN SS-tail and SSS terms */
        phaseTidal += pfaN * c3p5pn * powers_of_lalpi.two_thirds * powers_of_Mf->two_thirds;
    }
    else
    {
      XLAL_ERROR( XLAL_EINVAL, "Error in IMRPhenomX_TidalPhase: Unsupported NRTidal_version. This function currently only supports NRTidalv2 and NRTidalv3.\n");
    }

    return(phaseTidal);


}

REAL8 IMRPhenomX_TidalPhaseDerivative(IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXWaveformStruct *pWF, IMRPhenomXPhaseCoefficients *pPhase, NRTidal_version_type NRTidal_version){

    REAL8 X_A = pWF->m1; // Already scaled by Mtot
    REAL8 X_B = pWF->m2;

    REAL8 Mf_ten_thirds=powers_of_Mf->eight_thirds*powers_of_Mf->two_thirds;
    REAL8 pi_ten_thirds=powers_of_lalpi.eight_thirds*powers_of_lalpi.two_thirds;
    REAL8 Mf=powers_of_Mf->itself;

    REAL8 c2pn=pPhase->c2PN_tidal, c3pn=pPhase->c3PN_tidal, c3p5pn=pPhase->c3p5PN_tidal;

    REAL8 NRTuned_dphase=0.;

    REAL8 pfaN=3./(128.*X_A*X_B);

    REAL8 threePN_dphase=(pfaN*(-c2pn + c3pn*powers_of_Mf->two_thirds*powers_of_lalpi.two_thirds))/(3.*powers_of_Mf->four_thirds*powers_of_lalpi.one_third);

    REAL8 dphase=threePN_dphase;


    if(NRTidal_version==NRTidalv2_V){

        dphase+=(2.*c3p5pn*pfaN*powers_of_lalpi.two_thirds)/(3.*powers_of_Mf->one_third);

        REAL8 NRTidalv2_coeffs[9];

        int errcode;
        errcode = XLALSimNRTunedTidesSetFDTidalPhase_v2_Coeffs(NRTidalv2_coeffs);
        XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "Setting NRTidalv2 coefficients failed.\n");

        const REAL8 cNewt   = NRTidalv2_coeffs[0];
        const REAL8 n1      = NRTidalv2_coeffs[1];
        const REAL8 n3over2 = NRTidalv2_coeffs[2];
        const REAL8 n2      = NRTidalv2_coeffs[3];
        const REAL8 n5over2 = NRTidalv2_coeffs[4];
        const REAL8 n3      = NRTidalv2_coeffs[5];
        const REAL8 d1      = NRTidalv2_coeffs[6];
        const REAL8 d3over2 = NRTidalv2_coeffs[7];
        const REAL8 d2      = NRTidalv2_coeffs[8];

        REAL8 kappa2T = pWF->kappa2T;

        NRTuned_dphase=((cNewt*kappa2T*powers_of_Mf->two_thirds*powers_of_lalpi.five_thirds*(-5 - powers_of_Mf->two_thirds*(3*d1 + 7*n1)*powers_of_lalpi.two_thirds - 2*Mf*(d3over2 + 4*n3over2)*LAL_PI - powers_of_Mf->four_thirds*(d2 + 5*d1*n1 + 9*n2)*powers_of_lalpi.four_thirds - 2*powers_of_Mf->five_thirds*(2*d3over2*n1 + 3*d1*n3over2 + 5*n5over2)*powers_of_lalpi.five_thirds - pow(Mf,2)*(3*d2*n1 + 7*d1*n2 + 11*n3 + 5*d3over2*n3over2)*powers_of_lalpi.two - 2*powers_of_Mf->seven_thirds*(3*d3over2*n2 + 2*d2*n3over2 + 4*d1*n5over2)*powers_of_lalpi.seven_thirds - powers_of_Mf->eight_thirds*(5*d2*n2 + 9*d1*n3 + 7*d3over2*n5over2)*powers_of_lalpi.eight_thirds - 2*pow(Mf,3)*(4*d3over2*n3 + 3*d2*n5over2)*powers_of_lalpi.three - 7*d2*Mf_ten_thirds*n3*pi_ten_thirds))/(3.*pow(1 + d1*powers_of_Mf->two_thirds*powers_of_lalpi.two_thirds + d3over2*Mf*LAL_PI + d2*powers_of_Mf->four_thirds*powers_of_lalpi.four_thirds,2)*X_A*X_B));

    }
	else if(NRTidal_version==NRTidalv3_V){

        dphase+=(2.*c3p5pn*pfaN*powers_of_lalpi.two_thirds)/(3.*powers_of_Mf->one_third);

         /* local variable declaration */
		REAL8 Xa = X_A;

		REAL8 PN_coeffs[10];
		XLALSimNRTunedTidesSetFDTidalPhase_PN_Coeffs(PN_coeffs, Xa);

		REAL8 lambda1 = pWF->lambda1;
		REAL8 lambda2 = pWF->lambda2;

		REAL8 mtot = pWF->Mtot;

		REAL8 NRTidalv3_coeffs[20];
    	XLALSimNRTunedTidesSetFDTidalPhase_v3_Coeffs(NRTidalv3_coeffs, Xa, mtot, lambda1, lambda2, PN_coeffs);

		REAL8 s1 = NRTidalv3_coeffs[0];
		REAL8 s2 = NRTidalv3_coeffs[1];
		REAL8 s3 = NRTidalv3_coeffs[2];

		REAL8 s2s3 = s2*s3;
		REAL8 s2Mf = -s2*Mf*LAL_PI*2.0;
		REAL8 exps2s3 = cosh(s2s3) + sinh(s2s3);
		REAL8 exps2Mf = cosh(s2Mf) + sinh(s2Mf);

		/* Rewriting Eq. (27) of arXiv:2311.07456 */
  		REAL8 dynk2barfunc = 1.0 + ((s1) - 1)*(1.0/(1.0 + exps2Mf*exps2s3)) - ((s1-1.0)/(1.0 + exps2s3)) - 2.0*(Mf*LAL_PI)*((s1) - 1)*s2*exps2s3/((1.0 + exps2s3)*(1.0 + exps2s3));

		REAL8 kappaA = NRTidalv3_coeffs[4];
		REAL8 kappaB = NRTidalv3_coeffs[5];

		REAL8 dynkappaA = kappaA*dynk2barfunc;
		REAL8 dynkappaB = kappaB*dynk2barfunc;

		/* Pade Coefficients, Table II of arXiv:2311.07456 */
		REAL8 n_5over2A = NRTidalv3_coeffs[6];
		REAL8 n_3A = NRTidalv3_coeffs[7];
		REAL8 d_1A = NRTidalv3_coeffs[8];

		REAL8 n_5over2B = NRTidalv3_coeffs[9];
		REAL8 n_3B = NRTidalv3_coeffs[10];
		REAL8 d_1B = NRTidalv3_coeffs[11];

		/* 7.5PN Coefficients */
		REAL8 c_NewtA = PN_coeffs[0];
		REAL8 c_1A = PN_coeffs[1];
		REAL8 c_3over2A = PN_coeffs[2];
		REAL8 c_2A = PN_coeffs[3];
		REAL8 c_5over2A = PN_coeffs[4];

		REAL8 c_NewtB = PN_coeffs[5];
		REAL8 c_1B = PN_coeffs[6];
		REAL8 c_3over2B = PN_coeffs[7];
		REAL8 c_2B = PN_coeffs[8];
		REAL8 c_5over2B = PN_coeffs[9];

		/* Pade Coefficients constrained with PN, see Eq. (33) in arXiv:2311.07456 */
		REAL8 n_1A = NRTidalv3_coeffs[12];
		REAL8 n_3over2A = NRTidalv3_coeffs[13];
		REAL8 n_2A = NRTidalv3_coeffs[14];
		REAL8 d_3over2A = NRTidalv3_coeffs[15];

		REAL8 n_1B = NRTidalv3_coeffs[16];
		REAL8 n_3over2B = NRTidalv3_coeffs[17];
		REAL8 n_2B = NRTidalv3_coeffs[18];
		REAL8 d_3over2B = NRTidalv3_coeffs[19];

		/* Rewriting Eq. (30) and (32) in arXiv:2311.07456 in terms of Mf instead of x */
		REAL8 NRphasetermA=-((c_NewtA*dynkappaA*powers_of_Mf->five_thirds*powers_of_lalpi.five_thirds*(1 + powers_of_Mf->two_thirds*n_1A*powers_of_lalpi.two_thirds + Mf*n_3over2A*LAL_PI + powers_of_Mf->four_thirds*n_2A*powers_of_lalpi.four_thirds + powers_of_Mf->five_thirds*n_5over2A*powers_of_lalpi.five_thirds + pow(Mf,2)*n_3A*powers_of_lalpi.two))/((1 + d_1A*powers_of_Mf->two_thirds*powers_of_lalpi.two_thirds + d_3over2A*Mf*LAL_PI )));
        REAL8 NRphasetermB=-((c_NewtB*dynkappaB*powers_of_Mf->five_thirds*powers_of_lalpi.five_thirds*(1 + powers_of_Mf->two_thirds*n_1B*powers_of_lalpi.two_thirds + Mf*n_3over2B*LAL_PI + powers_of_Mf->four_thirds*n_2B*powers_of_lalpi.four_thirds + powers_of_Mf->five_thirds*n_5over2B*powers_of_lalpi.five_thirds + pow(Mf,2)*n_3B*powers_of_lalpi.two))/((1 + d_1B*powers_of_Mf->two_thirds*powers_of_lalpi.two_thirds + d_3over2B*Mf*LAL_PI )));

        REAL8 NRphaseNRT = NRphasetermA + NRphasetermB;

		REAL8 PNtidalphaseA = -((c_NewtA*kappaA*powers_of_Mf->five_thirds*powers_of_lalpi.five_thirds*(1 + powers_of_Mf->two_thirds*c_1A*powers_of_lalpi.two_thirds + Mf*c_3over2A*LAL_PI + powers_of_Mf->four_thirds*c_2A*powers_of_lalpi.four_thirds + powers_of_Mf->five_thirds*c_5over2A*powers_of_lalpi.five_thirds )));
		REAL8 PNtidalphaseB = -((c_NewtB*kappaB*powers_of_Mf->five_thirds*powers_of_lalpi.five_thirds*(1 + powers_of_Mf->two_thirds*c_1B*powers_of_lalpi.two_thirds + Mf*c_3over2B*LAL_PI + powers_of_Mf->four_thirds*c_2B*powers_of_lalpi.four_thirds + powers_of_Mf->five_thirds*c_5over2B*powers_of_lalpi.five_thirds )));

		REAL8 PNtidalphase = PNtidalphaseA + PNtidalphaseB;

		//DERIVATIVES

		REAL8 s2s3b = -s2*s3;
		REAL8 s2Mfb = s2*Mf*LAL_PI*2.0;
		REAL8 exps2s3b = cosh(s2s3b) + sinh(s2s3b);
		REAL8 exps2Mfb = cosh(s2Mfb) + sinh(s2Mfb);

		REAL8 dynk2barfunc_deriv = (s1-1.0)*(2.0*LAL_PI)*s2*exps2Mfb*exps2s3b*(1/((1 + exps2Mfb*exps2s3b)*(1 + exps2Mfb*exps2s3b))) -2.0*LAL_PI*((s1) - 1.0)*s2*exps2s3/((1.0 + exps2s3)*(1.0 + exps2s3));

		REAL8 NRTuned_dphaseA1=-((c_NewtA*kappaA*dynk2barfunc_deriv*powers_of_Mf->five_thirds*powers_of_lalpi.five_thirds*(1 + powers_of_Mf->two_thirds*n_1A*powers_of_lalpi.two_thirds + Mf*n_3over2A*LAL_PI + powers_of_Mf->four_thirds*n_2A*powers_of_lalpi.four_thirds + powers_of_Mf->five_thirds*n_5over2A*powers_of_lalpi.five_thirds + pow(Mf,2)*n_3A*powers_of_lalpi.two))/((1 + d_1A*powers_of_Mf->two_thirds*powers_of_lalpi.two_thirds + d_3over2A*Mf*LAL_PI)));
        REAL8 NRTuned_dphaseB1=-((c_NewtB*kappaB*dynk2barfunc_deriv*powers_of_Mf->five_thirds*powers_of_lalpi.five_thirds*(1 + powers_of_Mf->two_thirds*n_1B*powers_of_lalpi.two_thirds + Mf*n_3over2B*LAL_PI + powers_of_Mf->four_thirds*n_2B*powers_of_lalpi.four_thirds + powers_of_Mf->five_thirds*n_5over2B*powers_of_lalpi.five_thirds + pow(Mf,2)*n_3B*powers_of_lalpi.two))/((1 + d_1B*powers_of_Mf->two_thirds*powers_of_lalpi.two_thirds + d_3over2B*Mf*LAL_PI)));

		REAL8 NRTuned_dphaseA2=((c_NewtA*kappaA*dynk2barfunc*powers_of_Mf->two_thirds*powers_of_lalpi.five_thirds*(-5 - powers_of_Mf->two_thirds*(3*d_1A + 7*n_1A)*powers_of_lalpi.two_thirds - 2*Mf*(d_3over2A + 4*n_3over2A)*LAL_PI - powers_of_Mf->four_thirds*(5*d_1A*n_1A + 9*n_2A)*powers_of_lalpi.four_thirds - 2*powers_of_Mf->five_thirds*(2*d_3over2A*n_1A + 3*d_1A*n_3over2A + 5*n_5over2A)*powers_of_lalpi.five_thirds - pow(Mf,2)*(7*d_1A*n_2A + 11*n_3A + 5*d_3over2A*n_3over2A)*powers_of_lalpi.two - 2*powers_of_Mf->seven_thirds*(3*d_3over2A*n_2A + 4*d_1A*n_5over2A)*powers_of_lalpi.seven_thirds - powers_of_Mf->eight_thirds*(9*d_1A*n_3A + 7*d_3over2A*n_5over2A)*powers_of_lalpi.eight_thirds - 2*pow(Mf,3)*(4*d_3over2A*n_3A)*powers_of_lalpi.three))/(3.*pow(1 + d_1A*powers_of_Mf->two_thirds*powers_of_lalpi.two_thirds + d_3over2A*Mf*LAL_PI, 2)));
        REAL8 NRTuned_dphaseB2=((c_NewtB*kappaB*dynk2barfunc*powers_of_Mf->two_thirds*powers_of_lalpi.five_thirds*(-5 - powers_of_Mf->two_thirds*(3*d_1B + 7*n_1B)*powers_of_lalpi.two_thirds - 2*Mf*(d_3over2B + 4*n_3over2B)*LAL_PI - powers_of_Mf->four_thirds*(5*d_1B*n_1B + 9*n_2B)*powers_of_lalpi.four_thirds - 2*powers_of_Mf->five_thirds*(2*d_3over2B*n_1B + 3*d_1B*n_3over2B + 5*n_5over2B)*powers_of_lalpi.five_thirds - pow(Mf,2)*(7*d_1B*n_2B + 11*n_3B + 5*d_3over2B*n_3over2B)*powers_of_lalpi.two - 2*powers_of_Mf->seven_thirds*(3*d_3over2B*n_2B + 4*d_1B*n_5over2B)*powers_of_lalpi.seven_thirds - powers_of_Mf->eight_thirds*(9*d_1B*n_3B + 7*d_3over2B*n_5over2B)*powers_of_lalpi.eight_thirds - 2*pow(Mf,3)*(4*d_3over2B*n_3B)*powers_of_lalpi.three))/(3.*pow(1 + d_1B*powers_of_Mf->two_thirds*powers_of_lalpi.two_thirds + d_3over2B*Mf*LAL_PI, 2)));

		REAL8 factorA = -c_NewtA*kappaA;
		REAL8 factorB = -c_NewtB*kappaB;

		REAL8 NRTuned_dphaseNRT = NRTuned_dphaseA1 + NRTuned_dphaseA2 + NRTuned_dphaseB1 + NRTuned_dphaseB2;

		REAL8 PNpolyA = (1 + powers_of_Mf->two_thirds*c_1A*powers_of_lalpi.two_thirds + Mf*c_3over2A*LAL_PI + powers_of_Mf->four_thirds*c_2A*powers_of_lalpi.four_thirds + powers_of_Mf->five_thirds*c_5over2A*powers_of_lalpi.five_thirds);
		REAL8 PNpolyB = (1 + powers_of_Mf->two_thirds*c_1B*powers_of_lalpi.two_thirds + Mf*c_3over2B*LAL_PI + powers_of_Mf->four_thirds*c_2B*powers_of_lalpi.four_thirds + powers_of_Mf->five_thirds*c_5over2B*powers_of_lalpi.five_thirds);


		REAL8 dPNpolyA = (2./3.)*pow(Mf, -1.0/3.0)*c_1A*powers_of_lalpi.two_thirds + c_3over2A*LAL_PI + (4/3)*pow(Mf, 1.0/3.0)*c_2A*powers_of_lalpi.four_thirds + (5./3.)*powers_of_Mf->two_thirds*c_5over2A*powers_of_lalpi.five_thirds;
		REAL8 dPNpolyB = (2./3.)*pow(Mf, -1.0/3.0)*c_1B*powers_of_lalpi.two_thirds + c_3over2B*LAL_PI + (4/3)*pow(Mf, 1.0/3.0)*c_2B*powers_of_lalpi.four_thirds + (5./3.)*powers_of_Mf->two_thirds*c_5over2B*powers_of_lalpi.five_thirds;

		REAL8 PNTuned_dphaseA = factorA*((5./3.)*powers_of_Mf->two_thirds*powers_of_lalpi.five_thirds*PNpolyA + powers_of_Mf->five_thirds*powers_of_lalpi.five_thirds*dPNpolyA);
		REAL8 PNTuned_dphaseB = factorB*((5./3.)*powers_of_Mf->two_thirds*powers_of_lalpi.five_thirds*PNpolyB + powers_of_Mf->five_thirds*powers_of_lalpi.five_thirds*dPNpolyB);

		REAL8 PNTuned_dphase = PNTuned_dphaseA + PNTuned_dphaseB;

		REAL8 Mfmerger = XLALSimNRTunedTidesMergerFrequency_v3(pWF->Mtot, pWF->lambda1, pWF->lambda2, pWF->q, pWF->chi1L, pWF->chi2L) * (pWF->Mtot * LAL_MTSUN_SI) ;

		REAL8 Mftaperstart = 1.15*Mfmerger;
		REAL8 Mftaperend = 1.35*Mfmerger;


		/* Eq. (46) in arXiv:2311.07456 */
		REAL8 exp_arg;
		REAL8 exp_fn;
		REAL8 plancktaperfn;
		REAL8 dfactortaper = 0.0;
		REAL8 dPlancktaper = 0.0;

		if (Mf <= Mftaperstart){
		plancktaperfn = 0.0;
		}
		else if (Mf >= Mftaperend){
		plancktaperfn = 1.0;
		}
		else{
		exp_arg = (Mftaperend - Mftaperstart)/(Mf - Mftaperstart) + (Mftaperend - Mftaperstart)/(Mf - Mftaperend);
		exp_fn = cosh(exp_arg) + sinh(exp_arg);
		plancktaperfn = 1/(exp_fn + 1.0);
		dfactortaper = -(exp_fn)*(-(Mftaperend - Mftaperstart)/((Mf - Mftaperstart) * (Mf - Mftaperstart)) - (Mftaperend - Mftaperstart)/((Mf - Mftaperend) * (Mf - Mftaperend)));
		dPlancktaper = dfactortaper*(plancktaperfn * plancktaperfn);
		}

		NRTuned_dphase = NRTuned_dphaseNRT*(1.0 - plancktaperfn) + PNTuned_dphase*plancktaperfn - (NRphaseNRT - PNtidalphase)*dPlancktaper;
    }
    else
    {
      XLAL_ERROR( XLAL_EINVAL, "Error in IMRPhenomX_TidalPhaseDerivative: Unsupported NRTidal_version. This function currently only supports NRTidalv2 and NRTidalv3.\n");
    }

      dphase+=NRTuned_dphase;

    return(dphase);

}
