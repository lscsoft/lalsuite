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

#ifndef _LALSIMIMRSPINPRECEOBFACTORIZEDWAVEFORM_V3OPT_C
#define _LALSIMIMRSPINPRECEOBFACTORIZEDWAVEFORM_V3OPT_C
#include <gsl/gsl_deriv.h>
#include <complex.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

#include "LALSimIMREOBNRv2.h"
#include "LALSimIMRSpinEOB.h"

#include "LALSimIMREOBNewtonianMultipole.c"
#include "LALSimIMRSpinEOBAuxFuncs.c"
#include "LALSimIMREOBNQCCorrection.c"
#include "LALSimIMRSpinEOBHamiltonian.c"
#include "LALSimIMRSpinEOBHamiltonianPrec.c"

#include "LALSimIMRSpinPrecEOBHcapExactDerivative.c"
#include "LALSpinPrecHcapRvecDerivative_v3opt.c"
/* #include "LALSimIMRSpinEOBFactorizedWaveform.c" */
/*------------------------------------------------------------------------------------------
 *
 *          Prototypes of functions defined in this code.
 *
 *------------------------------------------------------------------------------------------
 */
static INT4
XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform_exact(
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
XLALSimIMRSpinEOBFluxGetPrecSpinFactorizedWaveform_exact(
						   COMPLEX16 * restrict hlmTab,	/**< OUTPUT, hlm waveforms */
					      REAL8Vector * restrict values,	/**< dyanmical variables: \f$(r,\phi,p_r,p_\phi)\f$ */
					  REAL8Vector * restrict cartvalues,	/**< dyanmical variables */
						   const REAL8 v,	/**< velocity */
						   const REAL8 Hreal,	/**< real Hamiltonian */
						   const INT4 lMax,	/**< maximum l mode to compute, compute 0 < m <= lMax */
					     SpinEOBParams * restrict params	/**< Spin EOB parameters */
);


static REAL8
XLALSimIMRSpinPrecEOBNonKeplerCoeff_exact(
                      const REAL8           values[],   /**<< Dynamical variables */
                      SpinEOBParams         *funcParams /**<< EOB parameters */
                      );

static REAL8 XLALSimIMRSpinPrecEOBCalcOmega_exact(
                      const REAL8           values[],   /**<< Dynamical variables */
                      SpinEOBParams         *funcParams /**<< EOB parameters */
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
XLALSimIMRSpinEOBFluxGetPrecSpinFactorizedWaveform_exact(
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
					XLALSimIMRSpinPrecEOBNonKeplerCoeff_exact(cartvalues->data, params);
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
UNUSED
static INT4
XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform_exact(
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
	REAL8		Slm     , deltalm, rholm, rholmPwrl;
	REAL8		auxflm = 0.0;
	REAL8		hathatksq4, hathatk4pi, Tlmprefac, Tlmprodfac;
	COMPLEX16	Tlm;
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
		vPhi = XLALSimIMRSpinPrecEOBNonKeplerCoeff_exact(cartvalues->data, params);


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
				auxflm = v * hCoeffs->f21v1  + v2 * v * hCoeffs->f21v3;
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
			rholm = 1. + v2 * (hCoeffs->rho33v2 + v * (hCoeffs->rho33v3 + v * (hCoeffs->rho33v4
											   + v * (hCoeffs->rho33v5 + v * (hCoeffs->rho33v6 + hCoeffs->rho33v6l * eulerlogxabs
															  + v * (hCoeffs->rho33v7 + (hCoeffs->rho33v8 + hCoeffs->rho33v8l * eulerlogxabs) * v))))));
			auxflm = v * v2 * hCoeffs->f33v3;
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

			deltalm = vh3 * (hCoeffs->delta44vh3 + hCoeffs->delta44vh6 * vh3)
				+ hCoeffs->delta44v5 * v2 * v2 * v;
			rholm = 1. + v2 * (hCoeffs->rho44v2
			     + v * (hCoeffs->rho44v3 + v * (hCoeffs->rho44v4
				 + v * (hCoeffs->rho44v5 + (hCoeffs->rho44v6
				+ hCoeffs->rho44v6l * eulerlogxabs) * v))));
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
			deltalm = hCoeffs->delta55vh3 * vh3 + hCoeffs->delta55v5 * v2 * v2 * v;
			rholm = 1. + v2 * (hCoeffs->rho55v2
			     + v * (hCoeffs->rho55v3 + v * (hCoeffs->rho55v4
			 + v * (hCoeffs->rho55v5 + hCoeffs->rho55v6 * v))));
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
		XLAL_PRINT_INFO("rho_%d_%d = %.12e \n", l, m, rholm);
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

	if (r > 0.0 && debugPK) {
		XLAL_PRINT_INFO("YP::dynamics variables in waveform: %i, %i, r = %.12e, v = %.12e\n", l, m, r, v);
		XLAL_PRINT_INFO("rholm^l = %.16e, Tlm = %.16e + i %.16e, \nSlm = %.16e, hNewton = %.16e + i %.16e, delta = %.16e\n", rholmPwrl, creal(Tlm), cimag(Tlm), Slm, creal(hNewton), cimag(hNewton), 0.0);
	}
	/* Put all factors in Eq. 17 together */
	*hlm = Tlm * cexp(I * deltalm) * Slm * rholmPwrl;
	*hlm *= hNewton;
	if (r > 8.5 && debugPK) {
		XLAL_PRINT_INFO("YP::FullWave: %.16e,%.16e, %.16e\n", creal(*hlm), cimag(*hlm), cabs(*hlm));
	}
	return XLAL_SUCCESS;
}


/*
 * Function to calculate the value of omega for the PRECESSING EOB waveform.
 * Needs the dynamics in Cartesian coordinates.
 *
 * First, the frame is rotated so that L is along the y-axis.
 * this rotation includes the spins.
 *
 * Second, \f$\vec{r}\f$ and \f$\vec{p}\f$ are converted to polar coordinates
 * (and not the spins). As L is along the y-axis, \f$\theta\f$ (defined as the
 * angle between L and the y-axis) is 0, which is a cyclic coordinate now and
 * that fixes
 * \f$p_\theta = 0\f$.
 *
 * Third, \f$p_r\f$ is set to 0.
 *
 * Fourth, the polar \f$(r,\phi,p_r=0,p_\phi)\f$ and the Cartesian spin vectors
 * are used to compute the derivative
 * \f$\partial Hreal/\partial p_\phi |p_r=0\f$.
 */

static REAL8 XLALSimIMRSpinPrecEOBCalcOmega_exact(
                      const REAL8           values[],   /**<< Dynamical variables */
                      SpinEOBParams         *funcParams /**<< EOB parameters */
                      )
{
  int debugPK = 1;
  if (debugPK){
    for(int i =0; i < 12; i++)
      if( isnan(values[i]) ) {
        XLAL_PRINT_INFO("XLALSimIMRSpinPrecEOBCalcOmega_exact::values %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f\n", values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9], values[10], values[11]);
          XLALPrintError( "XLAL Error - %s: nan in input values  \n", __func__);
          XLAL_ERROR( XLAL_EINVAL );
      }
  }

  /* ********************************************************************* */
  /* ************ Memory Allocation ************************************** */
  /* ********************************************************************* */
  static const REAL8 STEP_SIZE = 1.0e-4;
  REAL8 tmpvar = 0;

  HcapDerivParams params;

  /* Cartesian values for calculating the Hamiltonian */
  REAL8 cartValues[14] = {0.}, dvalues[14] = {0.};
  REAL8 cartvalues[14] = {0.}, polarvalues[6] = {0.}; /* The rotated cartesian/polar values */
  REAL8 polarRPcartSvalues[14] = {0.};
  memcpy( cartValues, values, 14 * sizeof(REAL8) );

  INT4 i, j;

  REAL8 rvec[3]  = {0.,0,0}, pvec[3]  = {0.,0,0};
  REAL8 s1vec[3] = {0.,0,0}, s2vec[3] = {0.,0,0};

  REAL8 rdotvec[3] = {0.,0,0};
  REAL8 rvecprime[3] = {0.,0,0}, pvecprime[3] = {0.,0,0},
        s1vecprime[3]= {0.,0,0}, s2vecprime[3]= {0.,0,0};
  REAL8 rvectmp[3]   = {0.,0,0}, pvectmp[3] = {0.,0,0},
        s1vectmp[3]  = {0.,0,0}, s2vectmp[3]= {0.,0,0};
  REAL8 LNhatprime[3]= {0.,0,0}, LNhatTmp[3]= {0.,0,0};
  REAL8 rcrossrdot[3] = {0.,0,0};

  REAL8 Rot1[3][3] ={{0.}}; // Rotation matrix for prevention of blowing up
  REAL8 Rot2[3][3] ={{0.}} ;
  REAL8 LNhat[3] = {0.,0,0};

  REAL8        Xhat[3] = {1, 0, 0};
  UNUSED REAL8 Yhat[3] = {0, 1, 0};
  UNUSED REAL8 Zhat[3] = {0, 0, 1};

  REAL8 Xprime[3] = {0.,0,0}, Yprime[3] = {0.,0,0}, Zprime[3] = {0.,0,0};

  gsl_function F;
  INT4         gslStatus;

  REAL8 omega;

  /* The error in a derivative as measured by GSL */
  REAL8 absErr;

  /* ********************************************************************* */
  /* ************ Main Logic begins ************************************ */
  /* ********************************************************************* */

  /* Copy over the coordinates and spins */
  memcpy( rvec,  values,   3*sizeof(REAL8) );
  memcpy( pvec,  values+3, 3*sizeof(REAL8) );
  memcpy( s1vec, values+6, 3*sizeof(REAL8) );
  memcpy( s2vec, values+9, 3*sizeof(REAL8) );

  /* Calculate rDot = \f$\partial Hreal / \partial p_r\f$ */
  memset( dvalues, 0, 14 * sizeof(REAL8) );
  if( XLALSpinPrecHcapRvecDerivative_exact( 0, values, dvalues,
                                  (void*) funcParams) == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }
  memcpy( rdotvec, dvalues, 3*sizeof(REAL8) );

  if (debugPK){
    for(int ii =0; ii < 12; ii++)
      if( isnan(dvalues[ii]) ) {
        XLAL_PRINT_INFO("XLALSimIMRSpinPrecEOBCalcOmega::dvalues %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f\n", dvalues[0], dvalues[1], dvalues[2], dvalues[3], dvalues[4], dvalues[5], dvalues[6], dvalues[7], dvalues[8], dvalues[9], dvalues[10], dvalues[11]);
          XLALPrintError( "XLAL Error - %s: nan in dvalues \n", __func__);
          XLAL_ERROR( XLAL_EINVAL );
      }
  }

  /* Calculate LN by computing r cross rDot and normalizing */
  cross_product( rvec, rdotvec, rcrossrdot );
  REAL8 rcrossrdotNorm = sqrt(inner_product( rcrossrdot, rcrossrdot ));
  for( i = 0; i < 3; i++ )
    rcrossrdot[i] /= rcrossrdotNorm;
  memcpy( LNhat, rcrossrdot, 3 * sizeof(REAL8) );


  /* ********************************************************************* */
  /* First, the frame is rotated so that L is along the y-axis. */
  /* this rotation includes the spins. */
  /* ********************************************************************* */

  /* For now, set first rotation matrix to identity.  If LNhat and Xhat are
     too closely aligned, rotate LNhat by pi/4 in the x,y plane.  Then
     LNhat determines Xprime, Yprime is the normalized cross product of
     Xprime and Xhat, and Zprime is the normalized cross product of Xprime
     and Zprime.
  */
  if( inner_product(LNhat, Xhat) < 0.9 )
  {
	  Rot1[0][0] = 1.; Rot1[0][1] = 0; Rot1[0][2] = 0;
	  Rot1[1][0] = 0.; Rot1[1][1] = 1; Rot1[1][2] = 0;
	  Rot1[2][0] = 0.; Rot1[2][1] = 0; Rot1[2][2] = 1;

	  memcpy(Xprime, LNhat, 3 * sizeof(REAL8));
  }
  else
  {
	  Rot1[0][0] = 1./sqrt(2); Rot1[0][1] = -1/sqrt(2); Rot1[0][2] = 0;
	  Rot1[1][0] = 1./sqrt(2); Rot1[1][1] = 1./sqrt(2); Rot1[1][2] = 0;
	  Rot1[2][0] = 0.;         Rot1[2][1] = 0;          Rot1[2][2] = 1;
	  LNhatTmp[0] = LNhatTmp[1] = LNhatTmp[2] = 0.;

	  for(i=0; i<3; i++)
	    for(j=0; j<3; j++)
	      LNhatTmp[i] += Rot1[i][j]*LNhat[j];

	  memcpy(Xprime, LNhatTmp, 3*sizeof(REAL8));
  }

  cross_product( Xprime, Xhat, Yprime );
  tmpvar = sqrt(inner_product(Yprime, Yprime));

  for( i=0; i<3; i++) Yprime[i] /= tmpvar;

  cross_product(Xprime, Yprime, Zprime);
  tmpvar = sqrt(inner_product(Zprime, Zprime));
  for( i=0; i<3; i++) Zprime[i] /= tmpvar;

  Rot2[0][0] = Xprime[0]; Rot2[0][1] = Xprime[1]; Rot2[0][2] = Xprime[2];
  Rot2[1][0] = Yprime[0]; Rot2[1][1] = Yprime[1]; Rot2[1][2] = Yprime[2];
  Rot2[2][0] = Zprime[0]; Rot2[2][1] = Zprime[1]; Rot2[2][2] = Zprime[2];

  memset( rvectmp,    0, 3 * sizeof(REAL8) );
  memset( pvectmp,    0, 3 * sizeof(REAL8) );
  memset( s1vectmp,   0, 3 * sizeof(REAL8) );
  memset( s2vectmp,   0, 3 * sizeof(REAL8) );
  memset( rvecprime,  0, 3 * sizeof(REAL8) );
  memset( pvecprime,  0, 3 * sizeof(REAL8) );
  memset( s1vecprime, 0, 3 * sizeof(REAL8) );
  memset( s2vecprime, 0, 3 * sizeof(REAL8) );
  memset( LNhatprime, 0, 3 * sizeof(REAL8) );
  memset( LNhatTmp,   0, 3 * sizeof(REAL8) );

  /* Perform the actual rotation */
  for (i=0; i<3; i++)
    for(j=0; j<3; j++)
      {
        rvectmp[i]  += Rot1[i][j]*rvec[j];
        pvectmp[i]  += Rot1[i][j]*pvec[j];
        s1vectmp[i] += Rot1[i][j]*s1vec[j];
        s2vectmp[i] += Rot1[i][j]*s2vec[j];
        LNhatTmp[i] += Rot1[i][j]*LNhat[j];
      }
  for (i=0; i<3; i++)
    for(j=0; j<3; j++)
      {
        rvecprime[i]  += Rot2[i][j]*rvectmp[j];
        pvecprime[i]  += Rot2[i][j]*pvectmp[j];
        s1vecprime[i] += Rot2[i][j]*s1vectmp[j];
        s2vecprime[i] += Rot2[i][j]*s2vectmp[j];
        LNhatprime[i] += Rot2[i][j]*LNhatTmp[j];
      }

  memcpy(cartvalues,   rvecprime,  3*sizeof(REAL8));
  memcpy(cartvalues+3, pvecprime,  3*sizeof(REAL8));
  memcpy(cartvalues+6, s1vecprime, 3*sizeof(REAL8));
  memcpy(cartvalues+9, s2vecprime, 3*sizeof(REAL8));

  /* ********************************************************************* */
  /* Second, \f$\vec{r}\f$ and \f$\vec{p}\f$ are converted to polar
   * coordinates (and not the spins).
   * As L is along the y-axis, \f$\theta\f$ (defined as the angle between
   * L and the y-axis) is 0, which is a cyclic coordinate now and that fixes
   * \f$p_\theta = 0\f$. */
  /* ********************************************************************* */

  /** The polarvalues, respectively, are
   * \f${r, \theta, \phi, p_r, p_\theta, p_\phi}\f$.
   * Below we transform \f$\vec{r}\f$ and \f$\vec{p}\f$ to the six
   * corresponding polar coordinates. */
  polarvalues[0] = sqrt(inner_product(rvecprime,rvecprime));
  polarvalues[1] = acos(rvecprime[0] / polarvalues[0]);
  polarvalues[2] = atan2(-rvecprime[1], rvecprime[2]);

  /* FIX p_r = 0 */
  polarvalues[3] = 0;

  REAL8 rvecprime_x_xhat[3] = {0.}, rvecprime_x_xhat_x_rvecprime[3] = {0.};
  cross_product(rvecprime, Xhat, rvecprime_x_xhat);
  cross_product(rvecprime_x_xhat, rvecprime, rvecprime_x_xhat_x_rvecprime);

  polarvalues[4] = -inner_product(rvecprime_x_xhat_x_rvecprime, pvecprime)
                              / polarvalues[0] / sin(polarvalues[1]);
  polarvalues[5] = -inner_product(rvecprime_x_xhat, pvecprime);


  /* ********************************************************************* */
  /* Finally, Differentiate Hamiltonian w.r.t. p_\phi, keeping p_r = 0 */
  /* ********************************************************************* */

  /* Populate the vector specifying the dynamical variables in mixed frames */
  memcpy( polarRPcartSvalues, cartvalues, 12*sizeof(REAL8));
  memcpy( polarRPcartSvalues, polarvalues, 6*sizeof(REAL8));

  /* Set up pointers for GSL */
  params.values  = polarRPcartSvalues;
  params.params  = funcParams;

  F.function = &GSLSpinPrecHamiltonianWrapperFordHdpphi;
  F.params   = &params;

  /* Now calculate omega. In the chosen co-ordinate system, */
  /* we need dH/dpphi to calculate this, i.e. varyParam = 5 */
  params.varyParam = 5;
  XLAL_CALLGSL( gslStatus = gsl_deriv_central( &F, polarvalues[5],
                  STEP_SIZE, &omega, &absErr ) );

  if ( gslStatus != GSL_SUCCESS )
  {
    XLALPrintError( "XLAL Error - %s: Failure in GSL function\n", __func__ );
    XLAL_ERROR_REAL8( XLAL_EFUNC );
  }

  return omega;
}

/**
 * Function to calculate the non-Keplerian coefficient for the PRECESSING EOB
 * model.
 *
 * radius \f$r\f$ times the cuberoot of the returned number is \f$r_\Omega\f$
 * defined in Eq. A2, i.e. the function returns
 * \f$(r_{\Omega} / r)^3\f$
 *     = \f$1/(r^3 (\partial Hreal/\partial p_\phi |p_r=0)^2)\f$.
 */
static REAL8
XLALSimIMRSpinPrecEOBNonKeplerCoeff_exact(
                      const REAL8           values[],   /**<< Dynamical variables */
                      SpinEOBParams         *funcParams /**<< EOB parameters */
                      )
{
  int debugPK = 1;
  if (debugPK){
    for(int i =0; i < 12; i++)
      if( isnan(values[i]) ) {
        XLAL_PRINT_INFO("XLALSimIMRSpinPrecEOBNonKeplerCoeff_exact::values %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f\n",
        values[0], values[1], values[2], values[3], values[4], values[5],
        values[6], values[7], values[8], values[9], values[10], values[11]);
          XLALPrintError( "XLAL Error - %s: nan in values  \n", __func__);
          XLAL_ERROR( XLAL_EINVAL );
      }
  }

  REAL8 omegaCirc = 0;
  REAL8 tmpValues[14]= {0.};
  REAL8 r3;

  /* We need to find the values of omega assuming pr = 0 */
  memcpy( tmpValues, values, sizeof(tmpValues) );
  omegaCirc = XLALSimIMRSpinPrecEOBCalcOmega_exact( tmpValues, funcParams );

  if ( XLAL_IS_REAL8_FAIL_NAN( omegaCirc ) )
  {
    XLAL_ERROR_REAL8( XLAL_EFUNC );
  }

  /* inner_product operates on the first three entries of the vector
   * passed to it, and therefore inner_product(values,values) = r^2 */
  r3 = pow(inner_product(values, values), 3./2.);
  return 1.0/(omegaCirc*omegaCirc*r3);
}

#endif				/* _LALSIMIMRSPINPRECEOBFACTORIZEDWAVEFORM_V3OPT */
