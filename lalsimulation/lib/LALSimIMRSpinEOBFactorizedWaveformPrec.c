/**
 * \author Craig Robinson, Yi Pan, Prayush Kumar, Stas Babak, Andrea Taracchini
 *
 * \brief Function to compute the factorized waveform as used in the SEOBNRv1 model.
 * Waveform expressions are given by
 * Taracchini et al. ( PRD 86, 024011 (2012), arXiv 1202.0790 ).
 * All equation numbers in this file refer to equations of this paper,
 * unless otherwise specified.
 * Coefficients of the so-called "deltalm" terms are given by
 * Damour et al. PRD 79, 064004 (2009) and Pan et al. PRD 83, 064003 (2011),
 * henceforth DIN and PBFRT.
 *
 * Functions to compute the factorized waveform for the purpose of computing the
 * flux, i.e. returning only the absolute value of the multipoles. The tail term
 * Tlm is used in its resummed form, given by Eq. (42) of Damour, Nagar and
 * Bernuzzi, PRD 87 084035 (2013), called DNB here.
 *
 */

#ifndef _LALSIMIMRSPINPRECEOBFACTORIZEDWAVEFORM_C
#define _LALSIMIMRSPINPRECEOBFACTORIZEDWAVEFORM_C
#include <gsl/gsl_deriv.h>
#include <complex.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_linalg.h>

#include "LALSimIMREOBNRv2.h"
#include "LALSimIMRSpinEOB.h"

#include "LALSimIMREOBNewtonianMultipole.c"
#include "LALSimIMRSpinEOBAuxFuncs.c"
#include "LALSimIMREOBNQCCorrection.c"
#include "LALSimIMRSpinEOBHamiltonian.c"
#include "LALSimIMRSpinEOBHamiltonianPrec.c"
#include "LALSimIMRSpinEOBFactorizedWaveform.c"
#define v4Pwave 451
/* #include "LALSimIMRSpinEOBFactorizedWaveform.c" */


/*------------------------------------------------------------------------------------------
 *
 *          Prototypes of functions defined in this code.
 *
 *------------------------------------------------------------------------------------------
 */
static INT4
XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform(
					       COMPLEX16 * restrict hlm,	/**< OUTPUT, hlm waveforms */
					       REAL8Vector * restrict values,	/**< dyanmical variables: \f$(r,\phi,p_r,p_\phi)\f$ */
					  REAL8Vector * restrict cartvalues,	/**< dyanmical variables */
					       const REAL8 v,	/**< velocity */
					       const REAL8 Hreal,	/**< real Hamiltonian */
					       const INT4 l,	/**< l mode index */
					       const INT4 m,	/**< m mode index */
					     SpinEOBParams * restrict params	/**< Spin EOB parameters */
);

static INT4
XLALSimIMRSpinEOBFluxGetPrecSpinFactorizedWaveform(
						   COMPLEX16 * restrict hlmTab,	/**< OUTPUT, hlm waveforms */
					      REAL8Vector * restrict values,	/**< dyanmical variables: \f$(r,\phi,p_r,p_\phi)\f$ */
					  REAL8Vector * restrict cartvalues,	/**< dyanmical variables */
						   const REAL8 v,	/**< velocity */
						   const REAL8 Hreal,	/**< real Hamiltonian */
						   const INT4 lMax,	/**< maximum l mode to compute, compute 0 < m <= lMax */
					     SpinEOBParams * restrict params	/**< Spin EOB parameters */
);

UNUSED static INT4 XLALSimIMRSpinEOBGetAmplitudeResidualPrec (
	COMPLEX16 * restrict rholmpwrl,
	const REAL8 v,
	const REAL8 Hreal,
	const INT4 modeL,
	const INT4 modeM,
	SpinEOBParams * restrict params
);

UNUSED static INT4 XLALSimIMREOBCalcCalibCoefficientHigherModesPrec (
               SpinEOBParams * restrict UNUSED params, /**<< Output */
               const UINT4 modeL, /**<< Mode index L */
              const UINT4 modeM, /**<< Mode index M */
              SEOBdynamics *seobdynamics /** Dynamics array*/,
               const REAL8 timeorb, /**<< Time of the peak of the orbital frequency */
               const REAL8 m1, /**<< Component mass 1 */
               const REAL8 m2, /**<< Component mass 2 */
               const REAL8 UNUSED deltaT  /**<< Sampling interval */
						 );





/*------------------------------------------------------------------------------------------
 *
 *          Defintions of functions.
 *
 *------------------------------------------------------------------------------------------
 */


/** ######################################################### */
/** FOR PRECESSING EOB
 * This function calculates hlm mode factorized-resummed waveform
 * for given dynamical variables. This is optimized for flux calculation,
 * by ignoring complex arguments and keeping only absolute values.
 * Changes:
 * (i) Complex Argument of Tlm not exponentiated.
 * (ii) exp(i deltalm) set to 1.
 * Eq. 17 and the entire Appendix of PRD 86, 024011 (2012) + changes
 * described in ￼the section "Factorized waveforms" of https://dcc.ligo.org/T1400476
 */
static INT4
XLALSimIMRSpinEOBFluxGetPrecSpinFactorizedWaveform(
						   COMPLEX16 * restrict hlmTab,	/**< OUTPUT, hlm waveforms */
					      REAL8Vector * restrict values,	/**< dyanmical variables: \f$(r,\phi,p_r,p_\phi)\f$ */
					  REAL8Vector * restrict cartvalues,	/**< dyanmical variables */
						   const REAL8 v,	/**< velocity */
						   const REAL8 Hreal,	/**< real Hamiltonian */
						   const INT4 lMax,	/**< maximum l mode to compute, compute 0 < m <= lMax */
					     SpinEOBParams * restrict params	/**< Spin EOB parameters */
)
{
    int		debugPK = 0;
	const	REAL8 vPhiKepler = params->alignedSpins ?
					XLALSimIMRSpinAlignedEOBNonKeplerCoeff(values->data, params) :
					XLALSimIMRSpinPrecEOBNonKeplerCoeff(cartvalues->data, params);
	if (XLAL_IS_REAL8_FAIL_NAN(vPhiKepler)) {
		XLAL_ERROR(XLAL_EFUNC);
	}

	for (INT4 l = 2; l <= lMax; l++)
	{
		for (INT4 m = 1; m <= l; m++)
		{
			COMPLEX16 *hlm = &hlmTab[l*(lMax+1)+m];


			/* Status of function calls */
			INT4		status;
			INT4		i;

			REAL8		eta;
			REAL8 UNUSED	r , pp, Omega, v2, /* vh, vh3, */ k, hathatk, eulerlogxabs;
			//pr
				REAL8 UNUSED rcrossp_x, rcrossp_y, rcrossp_z;
			REAL8		Slm     , rholm, rholmPwrl;
			REAL8		auxflm = 0.0;
			REAL8		hathatksq4, hathatk4pi, Tlmprefac, Tlmprodfac;
			REAL8		Tlm;
			COMPLEX16	hNewton;
			gsl_sf_result	z2;

			/* Non-Keplerian velocity */
			REAL8		vPhi    , vPhi2;

			/* Pre-computed coefficients */
			FacWaveformCoeffs *hCoeffs = params->eobParams->hCoeffs;

			if (abs(m) > (INT4) l) {
				XLAL_ERROR(XLAL_EINVAL);
			}
			if (m == 0) {
				XLAL_ERROR(XLAL_EINVAL);
			}
			eta = params->eobParams->eta;

			/* Check our eta was sensible */
			if (eta > 0.25) {
				XLALPrintError("XLAL Error - %s: Eta seems to be > 0.25 - this isn't allowed!\n", __func__);
				XLAL_ERROR(XLAL_EINVAL);
			}
			/*
			 * else if ( eta == 0.25 && m % 2 ) { // If m is odd and dM = 0, hLM
			 * will be zero memset( hlm, 0, sizeof( COMPLEX16 ) ); return
			 * XLAL_SUCCESS; }
			 */

			//r = sqrt(values->data[0] * values->data[0] + values->data[1] * values->data[1] + values->data[2] * values->data[2]);
			//pr = values->data[2];
			r = values->data[0];
			pp = values->data[3];

			rcrossp_x = cartvalues->data[1] * cartvalues->data[5] - cartvalues->data[2] * cartvalues->data[4];
			rcrossp_y = cartvalues->data[2] * cartvalues->data[3] - cartvalues->data[0] * cartvalues->data[5];
			rcrossp_z = cartvalues->data[0] * cartvalues->data[4] - cartvalues->data[1] * cartvalues->data[3];

			//pp = sqrt(rcrossp_x * rcrossp_x + rcrossp_y * rcrossp_y + rcrossp_z * rcrossp_z);

			v2 = v * v;
			Omega = v2 * v;
			//vh3 = Hreal * Omega;
			//vh = cbrt(vh3);
			eulerlogxabs = LAL_GAMMA + log(2.0 * (REAL8) m * v);

			/* Calculate the non-Keplerian velocity */
			vPhi = vPhiKepler;

			vPhi = r * cbrt(vPhi);

			if (debugPK)
				XLAL_PRINT_INFO("In XLALSimIMRSpinEOBFluxCalculateNewtonianMultipole, getting rW = %.12e\n",
				       vPhi);
			vPhi *= Omega;
			vPhi2 = vPhi * vPhi;

			/*
			 * Calculate the newtonian multipole, 1st term in Eq. 17, given by
			 * Eq. A1
			 */
			//debugPK
			if (debugPK) {
//				XLAL_PRINT_INFO("\nValues inside XLALSimIMRSpinEOBFluxGetPrecSpinFactorizedWaveform:\n");
//				for (i = 0; i < 14; i++)
//					XLAL_PRINT_INFO("values[%d] = %.12e\n", i, values->data[i]);

				XLAL_PRINT_INFO("Calculating hNewton, with v = %.12e, vPhi = %.12e, r = %.12e, Phi = %.12e, l = %d, m = %d\n",
				       v, vPhi, r, values->data[1], (UINT4) l, (UINT4) m);
			}
			status = XLALSimIMRSpinEOBFluxCalculateNewtonianMultipole(&hNewton,
				vPhi2, r, values->data[1], (UINT4) l, m, params->eobParams);
			if (status == XLAL_FAILURE) {
				XLAL_ERROR(XLAL_EFUNC);
			}
			/* Calculate the source term, 2nd term in Eq. 17, given by Eq. A5 */
			if (((l + m) % 2) == 0) {
				Slm = (Hreal * Hreal - 1.) / (2. * eta) + 1.;
			} else {
				Slm = v * pp;
				//Slm = v * sqrt(rcrossp_x * rcrossp_x + rcrossp_y * rcrossp_y + rcrossp_z * rcrossp_z);
			}
			if (debugPK)
				XLAL_PRINT_INFO("In XLALSimIMRSpinEOBFluxGetSpinFactorizedWaveform: Hreal = %e, Slm = %e, eta = %e\n", Hreal, Slm, eta);

			/*
			 * Calculate the absolute value of the Tail term, 3rd term in Eq. 17,
			 * given by Eq. A6, and Eq. (42) of
			 * http://arxiv.org/pdf/1212.4357.pdf
			 */
			k = m * Omega;
			hathatk = Hreal * k;
			hathatksq4 = 4. * hathatk * hathatk;
			hathatk4pi = 4. * LAL_PI * hathatk;
			/*
			 * gsl_sf_result lnr1, arg1; XLAL_CALLGSL( status =
			 * gsl_sf_lngamma_complex_e( l+1.0, -2.0*hathatk, &lnr1, &arg1 ) );
			 * if (status != GSL_SUCCESS) { XLALPrintError("XLAL Error - %s:
			 * Error in GSL function\n", __func__ ); XLAL_ERROR( XLAL_EFUNC ); }
			 */
			XLAL_CALLGSL(status = gsl_sf_fact_e(l, &z2));
			if (status != GSL_SUCCESS) {
				XLALPrintError("XLAL Error - %s: Error in GSL function\n", __func__);
				XLAL_ERROR(XLAL_EFUNC);
			}
			/*
			 * COMPLEX16 Tlmold; Tlmold = cexp( ( lnr1.val + LAL_PI * hathatk ) +
			 * I * ( arg1.val + 2.0 * hathatk * log(4.0*k/sqrt(LAL_E)) ) );
			 * Tlmold /= z2.val;
			 */
			/* Calculating the prefactor of Tlm, outside the multiple product */
			Tlmprefac = sqrt(hathatk4pi / (1. - exp(-hathatk4pi))) / z2.val;

			/* Calculating the multiple product factor */
			for (Tlmprodfac = 1., i = 1; i <= l; i++) {
				Tlmprodfac *= (hathatksq4 + (REAL8) i * i);
			}

			Tlm = Tlmprefac * sqrt(Tlmprodfac);

			/* Calculate the residue phase and amplitude terms */
			/*
			 * deltalm is the 4th term in Eq. 17, delta 22 given by Eq. A15,
			 * others
			 */
			/* rholm is the 5th term in Eq. 17, given by Eqs. A8 - A14 */
			/*
			 * auxflm is a special part of the 5th term in Eq. 17, given by Eq.
			 * A15
			 */
			/*
			 * Actual values of the coefficients are defined in the next function
			 * of this file
			 */
			switch (l) {
			case 2:
				switch (abs(m)) {
				case 2:
					rholm = 1. + v2 * (hCoeffs->rho22v2 + v * (hCoeffs->rho22v3
								     + v * (hCoeffs->rho22v4
					     + v * (hCoeffs->rho22v5 + v * (hCoeffs->rho22v6
									    + hCoeffs->rho22v6l * eulerlogxabs + v * (hCoeffs->rho22v7
														      + v * (hCoeffs->rho22v8 + hCoeffs->rho22v8l * eulerlogxabs
															     + (hCoeffs->rho22v10 + hCoeffs->rho22v10l * eulerlogxabs) * v2)))))));
					//FIXME
						if (debugPK){
                            XLAL_PRINT_INFO("PK:: rho22v2 = %.12e, rho22v3 = %.12e, rho22v4 = %.12e,\n rho22v5 = %.16e, rho22v6 = %.16e, rho22v6LOG = %.16e, \n rho22v7 = %.12e, rho22v8 = %.16e, rho22v8LOG = %.16e, \n rho22v10 = %.16e, rho22v10LOG = %.16e\n, rho22v6 = %.12e, rho22v8 = %.12e, rho22v10 = %.12e\n",
						       hCoeffs->rho22v2, hCoeffs->rho22v3, hCoeffs->rho22v4,
						       hCoeffs->rho22v5, hCoeffs->rho22v6, hCoeffs->rho22v6l,
						       hCoeffs->rho22v7, hCoeffs->rho22v8, hCoeffs->rho22v8l,
						       hCoeffs->rho22v10, hCoeffs->rho22v10l,
						       hCoeffs->rho22v6 + hCoeffs->rho22v6l * eulerlogxabs,
						       hCoeffs->rho22v8 + hCoeffs->rho22v8l * eulerlogxabs,
						       hCoeffs->rho22v10 + hCoeffs->rho22v10l * eulerlogxabs);}
					break;
				case 1:
					{
						rholm = 1. + v * (hCoeffs->rho21v1
								  + v * (hCoeffs->rho21v2 + v * (hCoeffs->rho21v3 + v * (hCoeffs->rho21v4
															 + v * (hCoeffs->rho21v5 + v * (hCoeffs->rho21v6 + hCoeffs->rho21v6l * eulerlogxabs
																			+ v * (hCoeffs->rho21v7 + hCoeffs->rho21v7l * eulerlogxabs
																			       + v * (hCoeffs->rho21v8 + hCoeffs->rho21v8l * eulerlogxabs
																				      + (hCoeffs->rho21v10 + hCoeffs->rho21v10l * eulerlogxabs) * v2))))))));
						auxflm = v * hCoeffs->f21v1 + v2 * v * hCoeffs->f21v3;
					}
					break;
				default:
					XLAL_ERROR(XLAL_EINVAL);
					break;
				}
				break;
			case 3:
				switch (m) {
				case 3:
					rholm = 1. + v2 * (hCoeffs->rho33v2 + v * (hCoeffs->rho33v3 + v * (hCoeffs->rho33v4
													   + v * (hCoeffs->rho33v5 + v * (hCoeffs->rho33v6 + hCoeffs->rho33v6l * eulerlogxabs
																	  + v * (hCoeffs->rho33v7 + (hCoeffs->rho33v8 + hCoeffs->rho33v8l * eulerlogxabs) * v))))));
					auxflm = v * v2 * hCoeffs->f33v3;
					break;
				case 2:
					rholm = 1. + v * (hCoeffs->rho32v
							  + v * (hCoeffs->rho32v2 + v * (hCoeffs->rho32v3 + v * (hCoeffs->rho32v4 + v * (hCoeffs->rho32v5
																	 + v * (hCoeffs->rho32v6 + hCoeffs->rho32v6l * eulerlogxabs
																		+ (hCoeffs->rho32v8 + hCoeffs->rho32v8l * eulerlogxabs) * v2))))));
					break;
				case 1:
					rholm = 1. + v2 * (hCoeffs->rho31v2 + v * (hCoeffs->rho31v3 + v * (hCoeffs->rho31v4
													   + v * (hCoeffs->rho31v5 + v * (hCoeffs->rho31v6 + hCoeffs->rho31v6l * eulerlogxabs
																	  + v * (hCoeffs->rho31v7 + (hCoeffs->rho31v8 + hCoeffs->rho31v8l * eulerlogxabs) * v))))));
					auxflm = v * v2 * hCoeffs->f31v3;
					break;
				default:
					XLAL_ERROR(XLAL_EINVAL);
					break;
				}
				break;
			case 4:
				switch (m) {
				case 4:

					rholm = 1. + v2 * (hCoeffs->rho44v2
					     + v * (hCoeffs->rho44v3 + v * (hCoeffs->rho44v4
						 + v * (hCoeffs->rho44v5 + (hCoeffs->rho44v6
						+ hCoeffs->rho44v6l * eulerlogxabs) * v))));
					break;
				case 3:
					rholm = 1. + v * (hCoeffs->rho43v
							  + v * (hCoeffs->rho43v2
					    + v2 * (hCoeffs->rho43v4 + v * (hCoeffs->rho43v5
									    + (hCoeffs->rho43v6 + hCoeffs->rho43v6l * eulerlogxabs) * v))));
					auxflm = v * hCoeffs->f43v;
					break;
				case 2:
					rholm = 1. + v2 * (hCoeffs->rho42v2
							   + v * (hCoeffs->rho42v3 + v * (hCoeffs->rho42v4 + v * (hCoeffs->rho42v5
														  + (hCoeffs->rho42v6 + hCoeffs->rho42v6l * eulerlogxabs) * v))));
					break;
				case 1:
					rholm = 1. + v * (hCoeffs->rho41v
							  + v * (hCoeffs->rho41v2
					    + v2 * (hCoeffs->rho41v4 + v * (hCoeffs->rho41v5
									    + (hCoeffs->rho41v6 + hCoeffs->rho41v6l * eulerlogxabs) * v))));
					auxflm = v * hCoeffs->f41v;
					break;
				default:
					XLAL_ERROR(XLAL_EINVAL);
					break;
				}
				break;
			case 5:
				switch (m) {
				case 5:
					rholm = 1. + v2 * (hCoeffs->rho55v2
					     + v * (hCoeffs->rho55v3 + v * (hCoeffs->rho55v4
					 + v * (hCoeffs->rho55v5 + hCoeffs->rho55v6 * v))));
					break;
				case 4:
					rholm = 1. + v2 * (hCoeffs->rho54v2 + v * (hCoeffs->rho54v3
								   + hCoeffs->rho54v4 * v));
					break;
				case 3:
					rholm = 1. + v2 * (hCoeffs->rho53v2
							   + v * (hCoeffs->rho53v3 + v * (hCoeffs->rho53v4 + hCoeffs->rho53v5 * v)));
					break;
				case 2:
					rholm = 1. + v2 * (hCoeffs->rho52v2 + v * (hCoeffs->rho52v3
								   + hCoeffs->rho52v4 * v));
					break;
				case 1:
					rholm = 1. + v2 * (hCoeffs->rho51v2
							   + v * (hCoeffs->rho51v3 + v * (hCoeffs->rho51v4 + hCoeffs->rho51v5 * v)));
					break;
				default:
					XLAL_ERROR(XLAL_EINVAL);
					break;
				}
				break;
			case 6:
				switch (m) {
				case 6:
					rholm = 1. + v2 * (hCoeffs->rho66v2 + v * (hCoeffs->rho66v3
								   + hCoeffs->rho66v4 * v));
					break;
				case 5:
					rholm = 1. + v2 * (hCoeffs->rho65v2 + hCoeffs->rho65v3 * v);
					break;
				case 4:
					rholm = 1. + v2 * (hCoeffs->rho64v2 + v * (hCoeffs->rho64v3
								   + hCoeffs->rho64v4 * v));
					break;
				case 3:
					rholm = 1. + v2 * (hCoeffs->rho63v2 + hCoeffs->rho63v3 * v);
					break;
				case 2:
					rholm = 1. + v2 * (hCoeffs->rho62v2 + v * (hCoeffs->rho62v3
								   + hCoeffs->rho62v4 * v));
					break;
				case 1:
					rholm = 1. + v2 * (hCoeffs->rho61v2 + hCoeffs->rho61v3 * v);
					break;
				default:
					XLAL_ERROR(XLAL_EINVAL);
					break;
				}
				break;
			case 7:
				switch (m) {
				case 7:
					rholm = 1. + v2 * (hCoeffs->rho77v2 + hCoeffs->rho77v3 * v);
					break;
				case 6:
					rholm = 1. + hCoeffs->rho76v2 * v2;
					break;
				case 5:
					rholm = 1. + v2 * (hCoeffs->rho75v2 + hCoeffs->rho75v3 * v);
					break;
				case 4:
					rholm = 1. + hCoeffs->rho74v2 * v2;
					break;
				case 3:
					rholm = 1. + v2 * (hCoeffs->rho73v2 + hCoeffs->rho73v3 * v);
					break;
				case 2:
					rholm = 1. + hCoeffs->rho72v2 * v2;
					break;
				case 1:
					rholm = 1. + v2 * (hCoeffs->rho71v2 + hCoeffs->rho71v3 * v);
					break;
				default:
					XLAL_ERROR(XLAL_EINVAL);
					break;
				}
				break;
			case 8:
				switch (m) {
				case 8:
					rholm = 1. + hCoeffs->rho88v2 * v2;
					break;
				case 7:
					rholm = 1. + hCoeffs->rho87v2 * v2;
					break;
				case 6:
					rholm = 1. + hCoeffs->rho86v2 * v2;
					break;
				case 5:
					rholm = 1. + hCoeffs->rho85v2 * v2;
					break;
				case 4:
					rholm = 1. + hCoeffs->rho84v2 * v2;
					break;
				case 3:
					rholm = 1. + hCoeffs->rho83v2 * v2;
					break;
				case 2:
					rholm = 1. + hCoeffs->rho82v2 * v2;
					break;
				case 1:
					rholm = 1. + hCoeffs->rho81v2 * v2;
					break;
				default:
					XLAL_ERROR(XLAL_EINVAL);
					break;
				}
				break;
			default:
				XLAL_ERROR(XLAL_EINVAL);
				break;
			}

			//debugPK
				if (debugPK)
				XLAL_PRINT_INFO("rho_%d_%d = %.12e \n", l, m, rholm);
			/* Raise rholm to the lth power */
			rholmPwrl = 1.0;
			i = l;
			while (i--) {
				rholmPwrl *= rholm;
			}
			/*
			 * In the equal-mass odd m case, there is no contribution from
			 * nonspin terms,  and the only contribution comes from the auxflm
			 * term that is proportional to chiA (asymmetric spins). In this
			 * case, we must ignore the nonspin terms directly, since the leading
			 * term defined by CalculateThisMultipolePrefix in
			 * LALSimIMREOBNewtonianMultipole.c is not zero (see comments there).
			 */
			if (eta == 0.25 && m % 2) {
				rholmPwrl = auxflm;
			} else {
				rholmPwrl += auxflm;
			}

			if (r > 0.0 && debugPK) {
				XLAL_PRINT_INFO("YP::dynamics variables in waveform: %i, %i, r = %.12e, v = %.12e\n", l, m, r, v);
				XLAL_PRINT_INFO("rholm^l = %.16e, Tlm = %.16e + i %.16e, \nSlm = %.16e, hNewton = %.16e + i %.16e, delta = %.16e\n", rholmPwrl, Tlm, 0.0, Slm, creal(hNewton), cimag(hNewton), 0.0);
			}
			/* Put all factors in Eq. 17 together */
			*hlm = Tlm * Slm * rholmPwrl;
			*hlm *= hNewton;
			/*
			 * if (r > 8.5) { XLAL_PRINT_INFO("YP::FullWave: %.16e,%.16e,
			 * %.16e\n",hlm->re,hlm->im,sqrt(hlm->re*hlm->re+hlm->im*hlm->im)); }
			 */
		}
	}
	return XLAL_SUCCESS;
}

/**
 * This function calculates hlm mode factorized-resummed waveform
 * for given dynamical variables.
 * Eq. 17 and the entire Appendix of PRD 86, 024011 (2012) + changes
 * described in ￼the section "Factorized waveforms" of https://dcc.ligo.org/T1400476
 */
static INT4
XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform(
					       COMPLEX16 * restrict hlm,	/**< OUTPUT, hlm waveforms */
					       REAL8Vector * restrict values,	/**< dyanmical variables: \f$(r,\phi,p_r,p_\phi)\f$ */
					  REAL8Vector * restrict cartvalues,	/**< dyanmical variables */
					       const REAL8 v,	/**< velocity */
					       const REAL8 Hreal,	/**< real Hamiltonian */
					       const INT4 l,	/**< l mode index */
					       const INT4 m,	/**< m mode index */
					     SpinEOBParams * restrict params	/**< Spin EOB parameters */
)
{
	int		debugPK = 0;
	/* Status of function calls */
	INT4		status;
	INT4		i;

	REAL8		eta;
	REAL8 UNUSED	r , pp, Omega, v2, Omegav2, vh, vh3, k, hathatk, eulerlogxabs;
	//pr
		REAL8 UNUSED rcrossp_x, rcrossp_y, rcrossp_z;
	REAL8		Slm     , deltalm;
	COMPLEX16		auxflm = 0.0;
	REAL8		hathatksq4, hathatk4pi, Tlmprefac, Tlmprodfac;
	COMPLEX16	Tlm, rholmPwrl,rholm;
	COMPLEX16	hNewton;
	gsl_sf_result	z2;

	/* Non-Keplerian velocity */
	REAL8		vPhi    , vPhi2;

	/* Pre-computed coefficients */

	FacWaveformCoeffs *hCoeffs = params->eobParams->hCoeffs;

	if (abs(m) > (INT4) l) {
		XLAL_ERROR(XLAL_EINVAL);
	}
	if (m == 0) {
		XLAL_ERROR(XLAL_EINVAL);
	}
	eta = params->eobParams->eta;

	/* Check our eta was sensible */
	if (eta > 0.25) {
		XLALPrintError("XLAL Error - %s: Eta seems to be > 0.25 - this isn't allowed!\n", __func__);
		XLAL_ERROR(XLAL_EINVAL);
	}
	/*
	 * else if ( eta == 0.25 && m % 2 ) { // If m is odd and dM = 0, hLM
	 * will be zero memset( hlm, 0, sizeof( COMPLEX16 ) ); return
	 * XLAL_SUCCESS; }
	 */

	r = values->data[0];
	pp = values->data[3];

	rcrossp_x = cartvalues->data[1] * cartvalues->data[5] - cartvalues->data[2] * cartvalues->data[4];
	rcrossp_y = cartvalues->data[2] * cartvalues->data[3] - cartvalues->data[0] * cartvalues->data[5];
	rcrossp_z = cartvalues->data[0] * cartvalues->data[4] - cartvalues->data[1] * cartvalues->data[3];

	//pp = sqrt(rcrossp_x * rcrossp_x + rcrossp_y * rcrossp_y + rcrossp_z * rcrossp_z);

	v2 = v * v;
	Omega = v2 * v;
    Omegav2 = Omega*v2;
	vh3 = Hreal * Omega;
	vh = cbrt(vh3);
	eulerlogxabs = LAL_GAMMA + log(2.0 * (REAL8) m * v);

	FacWaveformCoeffs oldCoeffs = *(params->eobParams->hCoeffs); //RC: the function XLALSimIMRSpinPrecEOBNonKeplerCoeff is calculating again the coefficients hCoeffs, for the omega. These are different because
	// for the dynamics we are using different coefficients for the 21 mode. I store the old coefficients and the I recover them ofter vPhi is computed.

	/* Calculate the non-Keplerian velocity */
	if (params->alignedSpins) {

		vPhi = XLALSimIMRSpinAlignedEOBNonKeplerCoeff(values->data, params);

		if (XLAL_IS_REAL8_FAIL_NAN(vPhi)) {
			XLAL_ERROR(XLAL_EFUNC);
		}
		vPhi = r * cbrt(vPhi);

		if (debugPK)
			XLAL_PRINT_INFO("In XLALSimIMRSpinEOBFluxCalculateNewtonianMultipole, getting rW = %.12e\n",
			       vPhi);
		vPhi *= Omega;
		vPhi2 = vPhi * vPhi;
	} else {
		vPhi = XLALSimIMRSpinPrecEOBNonKeplerCoeff(cartvalues->data, params);
		if (XLAL_IS_REAL8_FAIL_NAN(vPhi)) {
			XLAL_ERROR(XLAL_EFUNC);
		}
		vPhi = r * cbrt(vPhi);

		if (debugPK)
			XLAL_PRINT_INFO("In XLALSimIMRSpinEOBFluxCalculateNewtonianMultipole, getting rW = %.12e\n",
			       vPhi);
		vPhi *= Omega;
		vPhi2 = vPhi * vPhi;
	}
	*hCoeffs = oldCoeffs; //RC: Here I recover the old coefficients
	/*
	 * Calculate the newtonian multipole, 1st term in Eq. 17, given by
	 * Eq. A1
	 */
	if (debugPK) {
		XLAL_PRINT_INFO("\nValues inside XLALSimIMRSpinEOBFluxGetSpinFactorizedWaveform:\n");
		for (i = 0; i < 11; i++)
			XLAL_PRINT_INFO("values[%d] = %.12e\n", i, cartvalues->data[i]);

		XLAL_PRINT_INFO("Calculating hNewton, with v = %.12e, vPhi = %.12e, r = %.12e, Phi = %.12e, l = %d, m = %d\n",
		       v, vPhi, r, values->data[1], (UINT4) l, (UINT4) m);
	}
	status = XLALSimIMRSpinEOBCalculateNewtonianMultipole(&hNewton,
							      vPhi2, r, cartvalues->data[12] + cartvalues->data[13], (UINT4) l, m, params->eobParams);
							  
	if (status == XLAL_FAILURE) {
		XLAL_ERROR(XLAL_EFUNC);
	}
	/* Calculate the source term, 2nd term in Eq. 17, given by Eq. A5, Hreal is given by Eq.5 and Heff is in Eq.2 */
	if (((l + m) % 2) == 0) {
		Slm = (Hreal * Hreal - 1.) / (2. * eta) + 1.;
	} else {
		Slm = v * pp;
		//Slm = v * sqrt(rcrossp_x * rcrossp_x + rcrossp_y * rcrossp_y + rcrossp_z * rcrossp_z);
	}
	if (debugPK)
		XLAL_PRINT_INFO("In XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform: Hreal = %e, Slm = %e, eta = %e\n", Hreal, Slm, eta);

	/*
	 * Calculate the Tail term, 3rd term in Eq. 17,
	 * given by Eq. A6, and Eq. (42) of
	 * http://arxiv.org/pdf/1212.4357.pdf (or PRD 87 084035 (2013))
	 */
	k = m * Omega;
	hathatk = Hreal * k;

	gsl_sf_result	lnr1, arg1;
	XLAL_CALLGSL(status = gsl_sf_lngamma_complex_e(l + 1.0, -2.0 * hathatk, &lnr1, &arg1));
	if (status != GSL_SUCCESS) {
		XLALPrintError("XLAL Error - %s: Error in GSL function\n", __func__);
		XLAL_ERROR(XLAL_EFUNC);
	}
	XLAL_CALLGSL(status = gsl_sf_fact_e(l, &z2));
	if (status != GSL_SUCCESS) {
		XLALPrintError("XLAL Error - %s: Error in GSL function\n", __func__);
		XLAL_ERROR(XLAL_EFUNC);
	}
	Tlm = cexp((lnr1.val + LAL_PI * hathatk) + I * (
		    arg1.val + 2.0 * hathatk * log(4.0 * k / sqrt(LAL_E))));
	Tlm /= z2.val;


	if (debugPK){
        hathatksq4 = 4. * hathatk * hathatk;
        hathatk4pi = 4. * LAL_PI * hathatk;
        /* Calculating the prefactor of Tlm, outside the multiple product */
        Tlmprefac = sqrt(hathatk4pi / (1. - exp(-hathatk4pi))) / z2.val;

        /* Calculating the multiple product factor */
        for (Tlmprodfac = 1., i = 1; i <= l; i++)
            Tlmprodfac *= (hathatksq4 + (REAL8) i * i);

        REAL8		Tlmold;
        Tlmold = Tlmprefac * sqrt(Tlmprodfac);
		XLAL_PRINT_INFO("Tlm = %e + i%e, |Tlm| = %.16e (should be %.16e)\n", creal(Tlm), cimag(Tlm), cabs(Tlm), Tlmold);
    }

	/* Calculate the residue phase and amplitude terms */
	/*
	 * deltalm is the 4th term in Eq. 17, delta 22 given by Eq. A15,
	 * others
	 */
	/* rholm is the 5th term in Eq. 17, given by Eqs. A8 - A14 */
	/*
	 * auxflm is a special part of the 5th term in Eq. 17, given by Eq.
	 * A15
	 */
	/*
	 * Actual values of the coefficients are defined in the another function
	 * see file LALSimIMRSpinEOBFactorizedWaveformCoefficientsPrec.c
	 */
	switch (l) {
	case 2:
		switch (abs(m)) {
		case 2:
			deltalm = vh3 * (hCoeffs->delta22vh3 + vh3 * (hCoeffs->delta22vh6
				    + vh * vh * (hCoeffs->delta22vh9 * vh)))
				+ Omega*(hCoeffs->delta22v5 * v2 + Omega*(hCoeffs->delta22v6  + hCoeffs->delta22v8 *v2));
			rholm = 1. + v2 * (hCoeffs->rho22v2 + v * (hCoeffs->rho22v3
						     + v * (hCoeffs->rho22v4
			     + v * (hCoeffs->rho22v5 + v * (hCoeffs->rho22v6
							    + hCoeffs->rho22v6l * eulerlogxabs + v * (hCoeffs->rho22v7
												      + v * (hCoeffs->rho22v8 + hCoeffs->rho22v8l * eulerlogxabs
													     + (hCoeffs->rho22v10 + hCoeffs->rho22v10l * eulerlogxabs) * v2)))))));
			//FIXME
				if (debugPK){
                    XLAL_PRINT_INFO("PK:: rho22v2 = %.12e, rho22v3 = %.12e, rho22v4 = %.12e,\n rho22v5 = %.16e, rho22v6 = %.16e, rho22v6LOG = %.16e, \n rho22v7 = %.12e, rho22v8 = %.16e, rho22v8LOG = %.16e, \n rho22v10 = %.16e, rho22v10LOG = %.16e\n, rho22v6 = %.12e, rho22v8 = %.12e, rho22v10 = %.12e\n",
				       hCoeffs->rho22v2, hCoeffs->rho22v3, hCoeffs->rho22v4,
				       hCoeffs->rho22v5, hCoeffs->rho22v6, hCoeffs->rho22v6l,
				       hCoeffs->rho22v7, hCoeffs->rho22v8, hCoeffs->rho22v8l,
				       hCoeffs->rho22v10, hCoeffs->rho22v10l,
				       hCoeffs->rho22v6 + hCoeffs->rho22v6l * eulerlogxabs,
				       hCoeffs->rho22v8 + hCoeffs->rho22v8l * eulerlogxabs,
				       hCoeffs->rho22v10 + hCoeffs->rho22v10l * eulerlogxabs);}
			break;
		case 1:
			{
				deltalm = vh3 * (hCoeffs->delta21vh3 + vh3 * (hCoeffs->delta21vh6
									      + vh * (hCoeffs->delta21vh7 + (hCoeffs->delta21vh9) * vh * vh)))
					+ Omegav2*(hCoeffs->delta21v5  + hCoeffs->delta21v7 * v2);
				rholm = 1. + v * (hCoeffs->rho21v1
						  + v * (hCoeffs->rho21v2 + v * (hCoeffs->rho21v3 + v * (hCoeffs->rho21v4
													 + v * (hCoeffs->rho21v5 + v * (hCoeffs->rho21v6 + hCoeffs->rho21v6l * eulerlogxabs
																	+ v * (hCoeffs->rho21v7 + hCoeffs->rho21v7l * eulerlogxabs
																	       + v * (hCoeffs->rho21v8 + hCoeffs->rho21v8l * eulerlogxabs
																		      + (hCoeffs->rho21v10 + hCoeffs->rho21v10l * eulerlogxabs) * v2))))))));
				auxflm = v * (hCoeffs->f21v1 + v2 * (hCoeffs->f21v3 + v * hCoeffs->f21v4 + v2 * (hCoeffs->f21v5 + v * hCoeffs->f21v6 + v2 * hCoeffs->f21v7c)));
			}
			break;
		default:
			XLAL_ERROR(XLAL_EINVAL);
			break;
		}
		break;
	case 3:
		switch (m) {
		case 3:
			deltalm = vh3 * (hCoeffs->delta33vh3 + vh3 * (hCoeffs->delta33vh6 + hCoeffs->delta33vh9 * vh3))
				+ hCoeffs->delta33v5 * v * v2 * v2 + hCoeffs->delta33v7 * v2 * v2 * v2 * v;	
			//R.C: delta33v7 is set to 0, whoever is adding it here as a coefficient is evil, TODO: double check that is zero and then remove it
			//RC: This terms are in Eq.A6 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
			rholm = 1. + v2 * (hCoeffs->rho33v2 + v * (hCoeffs->rho33v3 + v * (hCoeffs->rho33v4 + v * (hCoeffs->rho33v5 + v * (hCoeffs->rho33v6 +
                	hCoeffs->rho33v6l * eulerlogxabs + v * (hCoeffs->rho33v7 + v * (hCoeffs->rho33v8 + hCoeffs->rho33v8l * eulerlogxabs +
                    v2*(hCoeffs->rho33v10 + hCoeffs->rho33v10l*eulerlogxabs))))))));
			//RC: This terms are in Eq.A10 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
        	auxflm = v * (v2 * (hCoeffs->f33v3 + v * (hCoeffs->f33v4 + v * (hCoeffs->f33v5  + v * hCoeffs->f33v6)))) + _Complex_I * vh3 * vh3 * hCoeffs->f33vh6;
			break;
		case 2:
			deltalm = vh3 * (hCoeffs->delta32vh3 + vh * (hCoeffs->delta32vh4 + vh * vh * (hCoeffs->delta32vh6
					     + hCoeffs->delta32vh9 * vh3)));
			rholm = 1. + v * (hCoeffs->rho32v
					  + v * (hCoeffs->rho32v2 + v * (hCoeffs->rho32v3 + v * (hCoeffs->rho32v4 + v * (hCoeffs->rho32v5
															 + v * (hCoeffs->rho32v6 + hCoeffs->rho32v6l * eulerlogxabs
																+ (hCoeffs->rho32v8 + hCoeffs->rho32v8l * eulerlogxabs) * v2))))));
			break;
		case 1:
			deltalm = vh3 * (hCoeffs->delta31vh3 + vh3 * (hCoeffs->delta31vh6
								      + vh * (hCoeffs->delta31vh7 + hCoeffs->delta31vh9 * vh * vh)))
				+ hCoeffs->delta31v5 * v * v2 * v2;
			rholm = 1. + v2 * (hCoeffs->rho31v2 + v * (hCoeffs->rho31v3 + v * (hCoeffs->rho31v4
											   + v * (hCoeffs->rho31v5 + v * (hCoeffs->rho31v6 + hCoeffs->rho31v6l * eulerlogxabs
															  + v * (hCoeffs->rho31v7 + (hCoeffs->rho31v8 + hCoeffs->rho31v8l * eulerlogxabs) * v))))));
			auxflm = v * v2 * hCoeffs->f31v3;
			break;
		default:
			XLAL_ERROR(XLAL_EINVAL);
			break;
		}
		break;
	case 4:
		switch (m) {
		case 4:
			//RC: This terms are in Eq.A15 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
      		deltalm = vh3 * (hCoeffs->delta44vh3 + vh3 * (hCoeffs->delta44vh6 + vh3 * hCoeffs->delta44vh9))
        			+ hCoeffs->delta44v5 * v2 * v2 * v;
			//RC: This terms are in Eq.A8 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
      		rholm = 1. + v2 * (hCoeffs->rho44v2
               		+ v * (hCoeffs->rho44v3 + v * (hCoeffs->rho44v4 + v * (hCoeffs->rho44v5 + v * (hCoeffs->rho44v6 + hCoeffs->rho44v6l *
                	eulerlogxabs + v2 *( hCoeffs->rho44v8 + hCoeffs->rho44v8l * eulerlogxabs +v2 * (hCoeffs->rho44v10 + hCoeffs->rho44v10l * eulerlogxabs) ) )))));
			break;
		case 3:
			deltalm = vh3 * (hCoeffs->delta43vh3 + vh * (hCoeffs->delta43vh4
					  + hCoeffs->delta43vh6 * vh * vh));
			rholm = 1. + v * (hCoeffs->rho43v
					  + v * (hCoeffs->rho43v2
			    + v2 * (hCoeffs->rho43v4 + v * (hCoeffs->rho43v5
							    + (hCoeffs->rho43v6 + hCoeffs->rho43v6l * eulerlogxabs) * v))));
			auxflm = v * hCoeffs->f43v;
			break;
		case 2:
			deltalm = vh3 * (hCoeffs->delta42vh3 + hCoeffs->delta42vh6 * vh3);
			rholm = 1. + v2 * (hCoeffs->rho42v2
					   + v * (hCoeffs->rho42v3 + v * (hCoeffs->rho42v4 + v * (hCoeffs->rho42v5
												  + (hCoeffs->rho42v6 + hCoeffs->rho42v6l * eulerlogxabs) * v))));
			break;
		case 1:
			deltalm = vh3 * (hCoeffs->delta41vh3 + vh * (hCoeffs->delta41vh4
					  + hCoeffs->delta41vh6 * vh * vh));
			rholm = 1. + v * (hCoeffs->rho41v
					  + v * (hCoeffs->rho41v2
			    + v2 * (hCoeffs->rho41v4 + v * (hCoeffs->rho41v5
							    + (hCoeffs->rho41v6 + hCoeffs->rho41v6l * eulerlogxabs) * v))));
			auxflm = v * hCoeffs->f41v;
			break;
		default:
			XLAL_ERROR(XLAL_EINVAL);
			break;
		}
		break;
	case 5:
		switch (m) {
		case 5:
			//RC: This terms are in Eq.A16 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
    		deltalm = vh3 *(hCoeffs->delta55vh3 +vh3*(hCoeffs->delta55vh6 +vh3 *(hCoeffs->delta55vh9)))
       				+ hCoeffs->delta55v5 * v2 * v2 * v;
			//RC: This terms are in Eq.A9 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
    		rholm = 1. + v2 * (hCoeffs->rho55v2 + v * (hCoeffs->rho55v3 + v * (hCoeffs->rho55v4 + v * (hCoeffs->rho55v5 +
            		v * (hCoeffs->rho55v6 + hCoeffs->rho55v6l * eulerlogxabs +
          			v2 * (hCoeffs->rho55v8 + hCoeffs->rho55v8l * eulerlogxabs +
        			v2 * (hCoeffs->rho55v10 + hCoeffs->rho55v10l * eulerlogxabs )))))));
		//RC: This terms are in Eq.A12 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
    		auxflm = v2 * v * (hCoeffs->f55v3 + v * (hCoeffs->f55v4 + v * (hCoeffs->f55v5c) ));
			break;
		case 4:
			deltalm = vh3 * (hCoeffs->delta54vh3 + hCoeffs->delta54vh4 * vh);
			rholm = 1. + v2 * (hCoeffs->rho54v2 + v * (hCoeffs->rho54v3
						   + hCoeffs->rho54v4 * v));
			break;
		case 3:
			deltalm = hCoeffs->delta53vh3 * vh3;
			rholm = 1. + v2 * (hCoeffs->rho53v2
					   + v * (hCoeffs->rho53v3 + v * (hCoeffs->rho53v4 + hCoeffs->rho53v5 * v)));
			break;
		case 2:
			deltalm = vh3 * (hCoeffs->delta52vh3 + hCoeffs->delta52vh4 * vh);
			rholm = 1. + v2 * (hCoeffs->rho52v2 + v * (hCoeffs->rho52v3
						   + hCoeffs->rho52v4 * v));
			break;
		case 1:
			deltalm = hCoeffs->delta51vh3 * vh3;
			rholm = 1. + v2 * (hCoeffs->rho51v2
					   + v * (hCoeffs->rho51v3 + v * (hCoeffs->rho51v4 + hCoeffs->rho51v5 * v)));
			break;
		default:
			XLAL_ERROR(XLAL_EINVAL);
			break;
		}
		break;
	case 6:
		switch (m) {
		case 6:
			deltalm = hCoeffs->delta66vh3 * vh3;
			rholm = 1. + v2 * (hCoeffs->rho66v2 + v * (hCoeffs->rho66v3
						   + hCoeffs->rho66v4 * v));
			break;
		case 5:
			deltalm = hCoeffs->delta65vh3 * vh3;
			rholm = 1. + v2 * (hCoeffs->rho65v2 + hCoeffs->rho65v3 * v);
			break;
		case 4:
			deltalm = hCoeffs->delta64vh3 * vh3;
			rholm = 1. + v2 * (hCoeffs->rho64v2 + v * (hCoeffs->rho64v3
						   + hCoeffs->rho64v4 * v));
			break;
		case 3:
			deltalm = hCoeffs->delta63vh3 * vh3;
			rholm = 1. + v2 * (hCoeffs->rho63v2 + hCoeffs->rho63v3 * v);
			break;
		case 2:
			deltalm = hCoeffs->delta62vh3 * vh3;
			rholm = 1. + v2 * (hCoeffs->rho62v2 + v * (hCoeffs->rho62v3
						   + hCoeffs->rho62v4 * v));
			break;
		case 1:
			deltalm = hCoeffs->delta61vh3 * vh3;
			rholm = 1. + v2 * (hCoeffs->rho61v2 + hCoeffs->rho61v3 * v);
			break;
		default:
			XLAL_ERROR(XLAL_EINVAL);
			break;
		}
		break;
	case 7:
		switch (m) {
		case 7:
			deltalm = hCoeffs->delta77vh3 * vh3;
			rholm = 1. + v2 * (hCoeffs->rho77v2 + hCoeffs->rho77v3 * v);
			break;
		case 6:
			deltalm = 0.0;
			rholm = 1. + hCoeffs->rho76v2 * v2;
			break;
		case 5:
			deltalm = hCoeffs->delta75vh3 * vh3;
			rholm = 1. + v2 * (hCoeffs->rho75v2 + hCoeffs->rho75v3 * v);
			break;
		case 4:
			deltalm = 0.0;
			rholm = 1. + hCoeffs->rho74v2 * v2;
			break;
		case 3:
			deltalm = hCoeffs->delta73vh3 * vh3;
			rholm = 1. + v2 * (hCoeffs->rho73v2 + hCoeffs->rho73v3 * v);
			break;
		case 2:
			deltalm = 0.0;
			rholm = 1. + hCoeffs->rho72v2 * v2;
			break;
		case 1:
			deltalm = hCoeffs->delta71vh3 * vh3;
			rholm = 1. + v2 * (hCoeffs->rho71v2 + hCoeffs->rho71v3 * v);
			break;
		default:
			XLAL_ERROR(XLAL_EINVAL);
			break;
		}
		break;
	case 8:
		switch (m) {
		case 8:
			deltalm = 0.0;
			rholm = 1. + hCoeffs->rho88v2 * v2;
			break;
		case 7:
			deltalm = 0.0;
			rholm = 1. + hCoeffs->rho87v2 * v2;
			break;
		case 6:
			deltalm = 0.0;
			rholm = 1. + hCoeffs->rho86v2 * v2;
			break;
		case 5:
			deltalm = 0.0;
			rholm = 1. + hCoeffs->rho85v2 * v2;
			break;
		case 4:
			deltalm = 0.0;
			rholm = 1. + hCoeffs->rho84v2 * v2;
			break;
		case 3:
			deltalm = 0.0;
			rholm = 1. + hCoeffs->rho83v2 * v2;
			break;
		case 2:
			deltalm = 0.0;
			rholm = 1. + hCoeffs->rho82v2 * v2;
			break;
		case 1:
			deltalm = 0.0;
			rholm = 1. + hCoeffs->rho81v2 * v2;
			break;
		default:
			XLAL_ERROR(XLAL_EINVAL);
			break;
		}
		break;
	default:
		XLAL_ERROR(XLAL_EINVAL);
		break;
	}

	//debugPK
		if (debugPK)
		XLAL_PRINT_INFO("rho_%d_%d = %.12e + %.12e \n", l, m, creal(rholm),cimag(rholm));
	if (debugPK)
		XLAL_PRINT_INFO("exp(delta_%d_%d) = %.16e + i%.16e (abs=%e)\n", l, m, creal(cexp(I * deltalm)),
		       cimag(cexp(I * deltalm)), cabs(cexp(I * deltalm)));
	/* Raise rholm to the lth power */
	rholmPwrl = 1.0;
	i = l;
	while (i--) {
		rholmPwrl *= rholm;
	}
	/*
	 * In the equal-mass odd m case, there is no contribution from
	 * nonspin terms,  and the only contribution comes from the auxflm
	 * term that is proportional to chiA (asymmetric spins). In this
	 * case, we must ignore the nonspin terms directly, since the leading
	 * term defined by CalculateThisMultipolePrefix in
	 * LALSimIMREOBNewtonianMultipole.c is not zero (see comments there).
	 */
	if (eta == 0.25 && m % 2) {
		rholmPwrl = auxflm;
	} else {
		rholmPwrl += auxflm;
	}
	// if (m==1) {
	// 	printf("f21v1 = %.16f f21v3 = %.16f f21v4 = %.16f f21v5 = %.16f\n", hCoeffs->f21v1, hCoeffs->f21v3, hCoeffs->f21v4, hCoeffs->f21v5);
	// }

	if (r > 0.0 && debugPK) {
		XLAL_PRINT_INFO("YP::dynamics variables in waveform: %i, %i, r = %.12e, v = %.12e\n", l, m, r, v);
		XLAL_PRINT_INFO("rholm^l = %.16e + %.16e, Tlm = %.16e + i %.16e, \nSlm = %.16e, hNewton = %.16e + i %.16e, delta = %.16e\n", creal(rholmPwrl), cimag(rholmPwrl), creal(Tlm), cimag(Tlm), Slm, creal(hNewton), cimag(hNewton), 0.0);
	}
	/* Put all factors in Eq. 17 together */
	*hlm = Tlm * cexp(I * deltalm) * Slm * rholmPwrl;
	*hlm *= hNewton;
	if (r > 8.5 && debugPK) {
		XLAL_PRINT_INFO("YP::FullWave: %.16e,%.16e, %.16e\n", creal(*hlm), cimag(*hlm), cabs(*hlm));
	}
	return XLAL_SUCCESS;
}

/**This function calculate the residue amplitude terms */
/* rholm is the 4th term in Eq. 4.5 of https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701] given by Eqs. A6 - A9 */
/* auxflm is a special part of the 4th term in Eq. 4.5 of https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701], given by Eq. A10-A12 */
UNUSED static INT4 XLALSimIMRSpinEOBGetAmplitudeResidualPrec (COMPLEX16 * restrict rholmpwrl, const REAL8 v, const UNUSED REAL8 Hreal, const INT4 modeL, const INT4 modeM, SpinEOBParams * restrict params)
  {
    INT4 i = 0;
    REAL8 eta;
    REAL8 eulerlogxabs;
    REAL8 rholm;
    REAL8 v2 = v * v;
    COMPLEX16 auxflm = 0.0;
    COMPLEX16 rholmPwrl;
    FacWaveformCoeffs *hCoeffs = params->eobParams->hCoeffs;
    eulerlogxabs = LAL_GAMMA + log (2.0 * (REAL8) modeM * v);

    if (abs (modeM) > (INT4) modeL)
      {
        XLAL_ERROR (XLAL_EINVAL);
      }
      if (modeM == 0)
      {
        XLAL_ERROR (XLAL_EINVAL);
      }
      eta = params->eobParams->eta;

      /* Check our eta was sensible */
      if (eta > 0.25 && eta < 0.25 +1e-4) {
        eta = 0.25;
      }
      if (eta > 0.25)
      {
        XLALPrintError
	       ("XLAL Error - %s: Eta seems to be > 0.25 - this isn't allowed!\n",
	        __func__);
          XLAL_ERROR (XLAL_EINVAL);
        }
  switch (modeL)
  {
          case 2:
      switch (abs (modeM))
	{
	case 2:
	  rholm =
	    1. + v2 * (hCoeffs->rho22v2 +
		       v * (hCoeffs->rho22v3 +
			    v * (hCoeffs->rho22v4 +
				 v * (hCoeffs->rho22v5 +
				      v * (hCoeffs->rho22v6 +
					   hCoeffs->rho22v6l * eulerlogxabs +
					   v * (hCoeffs->rho22v7 +
						v * (hCoeffs->rho22v8 +
						     hCoeffs->rho22v8l *
						     eulerlogxabs +
						     (hCoeffs->rho22v10 +
						      hCoeffs->rho22v10l *
						      eulerlogxabs) *
						     v2)))))));
	  break;
	case 1:
	  {
	    rholm =
	      1. + v * (hCoeffs->rho21v1 +
			v * (hCoeffs->rho21v2 +
			     v * (hCoeffs->rho21v3 +
				  v * (hCoeffs->rho21v4 +
				       v * (hCoeffs->rho21v5 +
					    v * (hCoeffs->rho21v6 +
						 hCoeffs->rho21v6l *
						 eulerlogxabs +
						 v * (hCoeffs->rho21v7 +
						      hCoeffs->rho21v7l *
						      eulerlogxabs +
						      v * (hCoeffs->rho21v8 +
							   hCoeffs->rho21v8l *
							   eulerlogxabs +
							   (hCoeffs->
							    rho21v10 +
							    hCoeffs->
							    rho21v10l *
							    eulerlogxabs) *
							   v2))))))));
                   auxflm = v * (hCoeffs->f21v1 + v2 * (hCoeffs->f21v3 + v * hCoeffs->f21v4 + v2 * (hCoeffs->f21v5 + v * hCoeffs->f21v6)));
               }
	  break;
	default:
	  XLAL_ERROR (XLAL_EINVAL);
	  break;
	}
      break;
    case 3:
      switch (modeM)
	{
	case 3:
          rholm =
              1. + v2 * (hCoeffs->rho33v2 +
                         v * (hCoeffs->rho33v3 +
                              v * (hCoeffs->rho33v4 +
                                   v * (hCoeffs->rho33v5 +
                                        v * (hCoeffs->rho33v6 +
                                             hCoeffs->rho33v6l * eulerlogxabs +
                                             v * (hCoeffs->rho33v7 +
                                                  v * (hCoeffs->rho33v8 +
                                                   hCoeffs->rho33v8l *
                                                   eulerlogxabs)))))));
          auxflm = v * v2 * hCoeffs->f33v3;
  break;
	case 2:
	  rholm =
	    1. + v * (hCoeffs->rho32v +
		      v * (hCoeffs->rho32v2 +
			   v * (hCoeffs->rho32v3 +
				v * (hCoeffs->rho32v4 +
				     v * (hCoeffs->rho32v5 +
					  v * (hCoeffs->rho32v6 +
					       hCoeffs->rho32v6l *
					       eulerlogxabs +
					       (hCoeffs->rho32v8 +
						hCoeffs->rho32v8l *
						eulerlogxabs) * v2))))));
	  break;
	case 1:
	  rholm =
	    1. + v2 * (hCoeffs->rho31v2 +
		       v * (hCoeffs->rho31v3 +
			    v * (hCoeffs->rho31v4 +
				 v * (hCoeffs->rho31v5 +
				      v * (hCoeffs->rho31v6 +
					   hCoeffs->rho31v6l * eulerlogxabs +
					   v * (hCoeffs->rho31v7 +
						(hCoeffs->rho31v8 +
						 hCoeffs->rho31v8l *
						 eulerlogxabs) * v))))));
	  auxflm = v * v2 * hCoeffs->f31v3;
	  break;
	default:
	  XLAL_ERROR (XLAL_EINVAL);
	  break;
	}
      break;
    case 4:
      switch (modeM)
	{
    case 4:
      rholm = 1. + v2 * (hCoeffs->rho44v2
               + v * (hCoeffs->rho44v3 + v * (hCoeffs->rho44v4
                      +
                      v *
                      (hCoeffs->
                       rho44v5 +
                       (hCoeffs->
                        rho44v6 +
                        hCoeffs->
                        rho44v6l *
                        eulerlogxabs) *
                       v))));

  	  break;
	case 3:
	  rholm =
	    1. + v * (hCoeffs->rho43v +
		      v * (hCoeffs->rho43v2 +
			   v2 * (hCoeffs->rho43v4 +
				 v * (hCoeffs->rho43v5 +
				      (hCoeffs->rho43v6 +
				       hCoeffs->rho43v6l * eulerlogxabs) *
				      v))));
	  auxflm = v * hCoeffs->f43v;
	  break;
	case 2:
	  rholm = 1. + v2 * (hCoeffs->rho42v2
			     + v * (hCoeffs->rho42v3 +
				    v * (hCoeffs->rho42v4 +
					 v * (hCoeffs->rho42v5 +
					      (hCoeffs->rho42v6 +
					       hCoeffs->rho42v6l *
					       eulerlogxabs) * v))));
	  break;
	case 1:
	  rholm =
	    1. + v * (hCoeffs->rho41v +
		      v * (hCoeffs->rho41v2 +
			   v2 * (hCoeffs->rho41v4 +
				 v * (hCoeffs->rho41v5 +
				      (hCoeffs->rho41v6 +
				       hCoeffs->rho41v6l * eulerlogxabs) *
				      v))));
	  auxflm = v * hCoeffs->f41v;
	  break;
	default:
	  XLAL_ERROR (XLAL_EINVAL);
	  break;
	}
      break;
    case 5:
      switch (modeM)
	{
    case 5:
      //RC: This terms are in Eq.A9 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
      rholm =
        1. + v2 * (hCoeffs->rho55v2 +
             v * (hCoeffs->rho55v3 +
            v * (hCoeffs->rho55v4 +
           v * (hCoeffs->rho55v5 +
              v * (hCoeffs->rho55v6 + hCoeffs->rho55v6l * eulerlogxabs +
            v2 * (hCoeffs->rho55v8 + hCoeffs->rho55v8l * eulerlogxabs +
          v2 * (hCoeffs->rho55v10 + hCoeffs->rho55v10l * eulerlogxabs )))))));
      //RC: This terms are in Eq.A12 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
      auxflm = v2 * v * (hCoeffs->f55v3 + v * (hCoeffs->f55v4));

  	  break;
	case 4:
	  rholm = 1. + v2 * (hCoeffs->rho54v2 + v * (hCoeffs->rho54v3
						     + hCoeffs->rho54v4 * v));
	  break;
	case 3:
	  rholm = 1. + v2 * (hCoeffs->rho53v2
			     + v * (hCoeffs->rho53v3 +
				    v * (hCoeffs->rho53v4 +
					 hCoeffs->rho53v5 * v)));
	  break;
	case 2:
	  rholm = 1. + v2 * (hCoeffs->rho52v2 + v * (hCoeffs->rho52v3
						     + hCoeffs->rho52v4 * v));
	  break;
	case 1:
	  rholm = 1. + v2 * (hCoeffs->rho51v2
			     + v * (hCoeffs->rho51v3 +
				    v * (hCoeffs->rho51v4 +
					 hCoeffs->rho51v5 * v)));
	  break;
	default:
	  XLAL_ERROR (XLAL_EINVAL);
	  break;
	}
      break;
    case 6:
      switch (modeM)
	{
	case 6:
	  rholm = 1. + v2 * (hCoeffs->rho66v2 + v * (hCoeffs->rho66v3
						     + hCoeffs->rho66v4 * v));
	  break;
	case 5:
	  rholm = 1. + v2 * (hCoeffs->rho65v2 + hCoeffs->rho65v3 * v);
	  break;
	case 4:
	  rholm = 1. + v2 * (hCoeffs->rho64v2 + v * (hCoeffs->rho64v3
						     + hCoeffs->rho64v4 * v));
	  break;
	case 3:
	  rholm = 1. + v2 * (hCoeffs->rho63v2 + hCoeffs->rho63v3 * v);
	  break;
	case 2:
	  rholm = 1. + v2 * (hCoeffs->rho62v2 + v * (hCoeffs->rho62v3
						     + hCoeffs->rho62v4 * v));
	  break;
	case 1:
	  rholm = 1. + v2 * (hCoeffs->rho61v2 + hCoeffs->rho61v3 * v);
	  break;
	default:
	  XLAL_ERROR (XLAL_EINVAL);
	  break;
	}
      break;
    case 7:
      switch (modeM)
	{
	case 7:
	  rholm = 1. + v2 * (hCoeffs->rho77v2 + hCoeffs->rho77v3 * v);
	  break;
	case 6:
	  rholm = 1. + hCoeffs->rho76v2 * v2;
	  break;
	case 5:
	  rholm = 1. + v2 * (hCoeffs->rho75v2 + hCoeffs->rho75v3 * v);
	  break;
	case 4:
	  rholm = 1. + hCoeffs->rho74v2 * v2;
	  break;
	case 3:
	  rholm = 1. + v2 * (hCoeffs->rho73v2 + hCoeffs->rho73v3 * v);
	  break;
	case 2:
	  rholm = 1. + hCoeffs->rho72v2 * v2;
	  break;
	case 1:
	  rholm = 1. + v2 * (hCoeffs->rho71v2 + hCoeffs->rho71v3 * v);
	  break;
	default:
	  XLAL_ERROR (XLAL_EINVAL);
	  break;
	}
      break;
    case 8:
      switch (modeM)
	{
	case 8:
	  rholm = 1. + hCoeffs->rho88v2 * v2;
	  break;
	case 7:
	  rholm = 1. + hCoeffs->rho87v2 * v2;
	  break;
	case 6:
	  rholm = 1. + hCoeffs->rho86v2 * v2;
	  break;
	case 5:
	  rholm = 1. + hCoeffs->rho85v2 * v2;
	  break;
	case 4:
	  rholm = 1. + hCoeffs->rho84v2 * v2;
	  break;
	case 3:
	  rholm = 1. + hCoeffs->rho83v2 * v2;
	  break;
	case 2:
	  rholm = 1. + hCoeffs->rho82v2 * v2;
	  break;
	case 1:
	  rholm = 1. + hCoeffs->rho81v2 * v2;
	  break;
	default:
	  XLAL_ERROR (XLAL_EINVAL);
	  break;
	}
      break;
    default:
      XLAL_ERROR (XLAL_EINVAL);
      break;
    }

  /* Raise rholm to the lth power */
  rholmPwrl = 1.0;
  i = modeL;
  while (i--)
    {
      rholmPwrl *= rholm;
    }
  /* In the equal-mass odd m case, there is no contribution from nonspin terms,
   * and the only contribution comes from the auxflm term that is proportional to chiA (asymmetric spins).
   * In this case, we must ignore the nonspin terms directly, since the leading term defined by
   * CalculateThisMultipolePrefix in LALSimIMREOBNewtonianMultipole.c is not zero (see comments there).
   */
  if (eta == 0.25 && modeM % 2)
    {
      rholmPwrl = auxflm;
    }
  else
    {
      rholmPwrl += auxflm;
    }
    *rholmpwrl = rholmPwrl;

return XLAL_SUCCESS;
}

/**
* This function calculate the calibration parameter for the amplitude of
* the factorized-resummed waveforms in the case of the 21 and the 55 mode
*/
UNUSED static INT4 XLALSimIMREOBCalcCalibCoefficientHigherModesPrec (
               SpinEOBParams * restrict UNUSED params, /**Output **/
               const UINT4 modeL, /*<< Mode index L */
               const UINT4 modeM, /*<< Mode index M */
               SEOBdynamics *seobdynamics, /*<< Dynamics vector */
               const REAL8 timeorb, /*<< Time of the peak of the orbital frequency */
               const REAL8 m1, /**Component mass 1 */
               const REAL8 m2, /**Component mass 2 */
               const REAL8 UNUSED deltaT  /**<< Sampling interval */
             )
               {
                 /* Status of function calls */
                 UINT4 i;
                 UINT4 debugRC = 0;

                 /** Physical quantities */
                 REAL8Vector *timeVec = NULL, *hLMdivrholmVec = NULL;
								 REAL8Vector polarDynamics, values;
								 REAL8 valuesdata[14] = {0.};
					       REAL8 polarDynamicsdata[4] = {0.};
                 polarDynamics.length = 4;
								 values.length = 14;
              	 polarDynamics.data = polarDynamicsdata;
								 values.data = valuesdata;
                 REAL8 omegaAttachmentPoint, hLMdivrholmAttachmentPoint, rholmNRAttachmentPoint,rholmBeforeCalibAttachmentPoint;
                 REAL8 v;
                 COMPLEX16Vector *hLMVec = NULL, *rholmpwrlVec = NULL;
                 REAL8Vector *rholmpwrlVecReal = NULL;

								 REAL8Vector orbOmegaVec, HamVec;
								 UINT4 retLen = seobdynamics->length;

								 /** Find the vaulues of the final spins */
								 REAL8 mtot = m1+m2;
								 REAL8 eta = params->eobParams->eta;
								 REAL8 s1dotZ = 0, s2dotZ = 0, chi1dotZ = 0, chi2dotZ = 0;
								 REAL8 chiS = 0;
								 REAL8 chiA = 0;
								 REAL8 tplspin = 0;
								 REAL8 spin1z_omegaPeak = params->spin1z_omegaPeak;
								 REAL8 spin2z_omegaPeak = params->spin2z_omegaPeak;
							   REAL8 chiS_omegaPeak = 0.5*(spin1z_omegaPeak+ spin2z_omegaPeak);
							   REAL8 chiA_omegaPeak = 0.5*(spin1z_omegaPeak-spin2z_omegaPeak);

                /* Create dynamical arrays */

                 orbOmegaVec.length = HamVec.length = retLen;

                 hLMVec = XLALCreateCOMPLEX16Vector (retLen);
                 rholmpwrlVec = XLALCreateCOMPLEX16Vector (retLen);
                 timeVec = XLALCreateREAL8Vector (retLen);
                 hLMdivrholmVec = XLALCreateREAL8Vector (retLen);
                 rholmpwrlVecReal = XLALCreateREAL8Vector (retLen);
                 if (!hLMVec || !rholmpwrlVec|| !timeVec|| !rholmpwrlVec || !rholmpwrlVecReal)
                     XLAL_ERROR (XLAL_ENOMEM);

                 orbOmegaVec.data = seobdynamics->omegaVec;
                 HamVec.data = seobdynamics->hamVec;


                 /* Stuff for interpolating function */
                 gsl_spline *spline = NULL;
                 gsl_interp_accel *acc = NULL;

                 /* The calibration parameter is only used for 21 and 55 modes, if you try to use this function for other modes, you get an error */

                 if (!((modeL == 2 && modeM == 1) || (modeL == 5 && modeM == 5)))
                   {
                     XLALPrintError ("Mode %d,%d is not supported by this function.\n", modeL,
                        modeM);
                     XLAL_ERROR (XLAL_EINVAL);
                   }
                 /* Populate time vector as necessary */
                 for (i = 0; i < timeVec->length; i++)
                   {
                     timeVec->data[i] = i * deltaT;
                   }
                 /**Initializing stuff for interpolation */
                 spline = gsl_spline_alloc (gsl_interp_cspline, orbOmegaVec.length);
                 acc = gsl_interp_accel_alloc ();
                 /* Calculation of the frequency at the attachment point */
								 REAL8 timewavePeak = timeorb-XLALSimIMREOBGetNRSpinPeakDeltaTv4 (modeL, modeM, m1, m2, spin1z_omegaPeak, spin2z_omegaPeak);
                 gsl_spline_init (spline, timeVec->data, orbOmegaVec.data, orbOmegaVec.length);
                 gsl_interp_accel_reset (acc);
                 omegaAttachmentPoint = gsl_spline_eval (spline, timewavePeak, acc);

                 /** Calculation and interpolation at the matching point of rho_lm^l + f_lm */
                for(i=0; i<orbOmegaVec.length; i++){
										 for (UINT4 j=0; j<14; j++) {
												 values.data[j] = seobdynamics->array->data[i + (j+1)*retLen];
										 }
										 s1dotZ = seobdynamics->s1dotZVec[i];
								     s2dotZ = seobdynamics->s2dotZVec[i];
										 polarDynamics.data[0] = seobdynamics->polarrVec[i];
										 polarDynamics.data[1] = seobdynamics->polarphiVec[i];
										 polarDynamics.data[2] = seobdynamics->polarprVec[i];
										 polarDynamics.data[3] = seobdynamics->polarpphiVec[i];
								     chi1dotZ = s1dotZ * mtot*mtot / (m1*m1);
								     chi2dotZ = s2dotZ * mtot*mtot / (m2*m2);
										 chiS = 0.5*(chi1dotZ+chi2dotZ);
										 chiA = 0.5*(chi1dotZ-chi2dotZ);
										 tplspin = (1.-2.*eta) * chiS + (m1 - m2)/(m1 + m2) * chiA;
                     v = cbrt (orbOmegaVec.data[i]);
										 if ( XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients( params->eobParams->hCoeffs, m1, m2, eta, tplspin, chiS, chiA, v4Pwave ) == XLAL_FAILURE )
								     {
								       XLALPrintError("XLAL Error - %s: failure in XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients at step %d of the loop.\n", __func__, i);
								       XLAL_ERROR( XLAL_EFUNC );
								     }
										 if ( XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform( &(hLMVec->data[i]), &polarDynamics, &values, v, HamVec.data[i], modeL, modeM, params ) == XLAL_FAILURE )
								     {
								       XLALPrintError("XLAL Error - %s: failure in XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform at step %d of the loop.\n", __func__, i);
								       XLAL_ERROR( XLAL_EDOM );
								     }
                     if (XLALSimIMRSpinEOBGetAmplitudeResidualPrec (&(rholmpwrlVec->data[i]), v, HamVec.data[i], modeL, modeM, params) == XLAL_FAILURE)
                     //RC: for the 21 and 55 mode rholmpwrlVec is always real. This is not true for the 33 mode. For this reason, in order to make this function general, we use a complex variable for it.
                     {
                         /* TODO: Clean-up */
                         XLAL_ERROR (XLAL_EFUNC);
                     }
                     rholmpwrlVecReal->data[i] = (REAL8)creal(rholmpwrlVec->data[i]);
                     hLMdivrholmVec->data[i] = ((REAL8)cabs(hLMVec->data[i]))/fabs(rholmpwrlVecReal->data[i]);
                    }
                    gsl_spline_init (spline, timeVec->data, hLMdivrholmVec->data, hLMdivrholmVec->length);
                    gsl_interp_accel_reset (acc);
                    hLMdivrholmAttachmentPoint = gsl_spline_eval (spline, timewavePeak, acc);
					REAL8 nra = 0;
                    nra = XLALSimIMREOBGetNRSpinPeakAmplitudeV4 (modeL, modeM, m1, m2,chiS_omegaPeak, chiA_omegaPeak);
                    if((fabs(nra/eta)< 3e-2) && ((modeL == 2) && (modeM == 1))){
                      //R.C.: safeguard to avoid the 21 mode to go to 0
                      nra = GSL_SIGN(nra)*eta*3e-2;
                    }
										if((fabs(nra/eta)< 1e-4) && ((modeL == 5) && (modeM == 5))){
                      //R.C.: safeguard to avoid the 55 mode to go to 0
                      nra = GSL_SIGN(nra)*eta*1e-4;
                    }
                    rholmNRAttachmentPoint = nra/hLMdivrholmAttachmentPoint;
                    gsl_spline_init (spline, timeVec->data, rholmpwrlVecReal->data, rholmpwrlVecReal->length);
                    gsl_interp_accel_reset (acc);
                    /* rho_lm^l + f_lm before the calibration parameter is set */
                    rholmBeforeCalibAttachmentPoint = gsl_spline_eval (spline, timewavePeak, acc);
                    if(debugRC==1){
                      FILE *timeomegapeak = NULL;
                      timeomegapeak = fopen ("timeomegapeak.dat", "w");
                      fprintf(timeomegapeak, "%.16f\n", timewavePeak);
                      fclose(timeomegapeak);
                    }
                    if ((modeL == 2) && (modeM == 1))
                      {
                      /* Here we compute ((rho_lm^l + f_lm + CalPar*omega^7/3)_NR - (rho_lm^l + f_lm)_EOB)/omega^7/3 to get CalPar.
                      The factor rholmpwrlVecReal->data[0])/cabs(rholmpwrlVecReal->data[0]) is used to know the sign of the function (rho_lm^l + f_lm + CalPar*omega^7/3)_NR which is computed as absolute value */
                      params->cal21 = (rholmNRAttachmentPoint - rholmBeforeCalibAttachmentPoint)/(pow(omegaAttachmentPoint,7./3.));
											//printf("cal21 = %.16f\n", params->cal21);
                      }
					if ((modeL == 5) && (modeM == 5))
                      {
                      /* Here we compute ((rho_lm^l + f_lm + CalPar*omega^7/3)_NR - (rho_lm^l + f_lm)_EOB)/omega^7/3 to get CalPar.
                      The factor rholmpwrlVecReal->data[0])/cabs(rholmpwrlVecReal->data[0]) is used to know the sign of the function (rho_lm^l + f_lm + CalPar*omega^7/3)_NR which is computed as absolute value */
                      params->cal55 = (rholmNRAttachmentPoint - rholmBeforeCalibAttachmentPoint)/(pow(omegaAttachmentPoint,5./3.));
											//printf("params->cal55 = %.16f\n",params->cal55);
                      }


                    gsl_spline_free (spline);
                    gsl_interp_accel_free (acc);
                    XLALDestroyREAL8Vector (hLMdivrholmVec);
                    XLALDestroyCOMPLEX16Vector (hLMVec);
                    XLALDestroyREAL8Vector (timeVec);
                    XLALDestroyCOMPLEX16Vector (rholmpwrlVec);
                    XLALDestroyREAL8Vector (rholmpwrlVecReal);



										return XLAL_SUCCESS;}

#endif				/* _LALSIMIMRSPINPRECEOBFACTORIZEDWAVEFORM */
