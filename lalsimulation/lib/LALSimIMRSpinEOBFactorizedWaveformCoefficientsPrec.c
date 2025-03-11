/*
 * \author Stas Babak, Prayush Kumar, Andrea Taracchini
 */
#ifndef _LALSIMIMRSPINPRECEOBFACTORIZEDWAVEFORMCOEFFICIENTS_C
#define _LALSIMIMRSPINPRECEOBFACTORIZEDWAVEFORMCOEFFICIENTS_C

#define v4Pwave 451


#include <complex.h>
#include <lal/LALSimInspiral.h>
//#include <lal/LALSimIMR.h>

#include "LALSimIMREOBNRv2.h"
#include "LALSimIMRSpinEOB.h"

/**
 *
 * Waveform expressions are given by
 * Taracchini et al. ( PRD 86, 024011 (2012), arXiv 1202.0790 ).
 * All equation numbers in this file refer to equations of this paper,
 * unless otherwise specified.
 * Coefficients of the so-called "deltalm" terms are given by
 * Damour et al. PRD 79, 064004 (2009) [arXiv:0811.2069] and Pan et al. PRD 83, 064003 (2011) [arXiv:1006.0431],
 * henceforth DIN and PBFRT.
 *
 * Functions to compute the factorized waveform for the purpose of computing the
 * flux, i.e. returning only the absolute value of the multipoles. The tail term
 * Tlm is used in its resummed form, given by Eq. (42) of Damour, Nagar and
 * Bernuzzi, PRD 87 084035 (2013) [arxiv:1212.4357], called DNB here.
 *
 */


/*--------------------------------------------------------------*/
/**
 * Spin Factors
 */

static int
XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(
					   FacWaveformCoeffs * const coeffs,	/**< OUTPUT, pre-computed waveform coefficients */
						 const REAL8 m1,	/**< mass 1 */
						 const REAL8 m2,	/**< mass 2 */
						 const REAL8 eta,	/**< symmetric mass ratio */
						 const REAL8 a,	/**< Kerr spin parameter for test-particle terms */
						 const REAL8 chiS,	/**< (chi1+chi2)/2 */
						 const REAL8 chiA,	/**< (chi1-chi2)/2 */
					   const UINT4 SpinAlignedEOBversion	/**< 1 for SEOBNRv1; 2 for SEOBNRv2; 3 for a small special treatment for SEOBNRv3*/
);


/****************************************************
 *	Definition of Functions
 *	**********************************************/


/**
 * This function calculates coefficients for hlm mode factorized-resummed waveform.
 * The coefficients are pre-computed and stored in the SpinEOBParams structure.
 * Eq. 17 and the entire Appendix of PRD 86, 024011 (2012) + changes
 * described in ￼the section "Factorized waveforms" of https://dcc.ligo.org/T1400476
 */

static int
XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(
					   FacWaveformCoeffs * const coeffs,	/**< OUTPUT, pre-computed waveform coefficients */
						 const REAL8 m1,	/**< mass 1 */
						 const REAL8 m2,	/**< mass 2 */
						 const REAL8 eta,	/**< symmetric mass ratio */
						 const REAL8 tmpa,	/**< Kerr spin parameter for test-particle terms */
						 const REAL8 chiS,	/**< (chi1+chi2)/2 */
						 const REAL8 chiA,	/**< (chi1-chi2)/2 */
					   UINT4 SpinAlignedEOBversion	/**< 1 for SEOBNRv1; 2 for SEOBNRv2; 4 for the coefficients
						 in the flux of v4P and v4Pwave for the coefficients in the waveform of v4P */
						 )
{
	int		debugPK = 0;
	if (debugPK) {
		XLAL_PRINT_INFO("In XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients: Renewing hlm coefficients.\n");
			XLAL_PRINT_INFO("PK:: chiS = %.12e, chiA = %.12e, a = %.12e (UNUSED), EOBVERSION = %d\n",
				   chiS, chiA, tmpa, SpinAlignedEOBversion);
	}
	REAL8		a = tmpa;
	INT4 waveform = 0;
	REAL8		eta2 = eta * eta;
	REAL8		eta3 = eta2 * eta;

	if (SpinAlignedEOBversion == v4Pwave) {
		SpinAlignedEOBversion = 4;
		waveform = 1;
	}
	 //RC: Here I am taking care of the fact that the flux used to calibrate SEOBNRv4 is missing some terms
	 //in the 21 mode that I have added later in SEOBNRv4HM. For this reason when computing the Coefficients
	 //for the 21 mode one has to specify if they are needed for the flux or for the waveform.

	REAL8		dM	  , dM2, chiA2, chiS2;
	//dM3;
	REAL8		aDelta  , a2, a3;

	/* Combination which appears a lot */
	REAL8		m1Plus3eta, m1Plus3eta2;
	REAL8 UNUSED m1Plus3eta3;

	dM2 = 1. - 4. * eta;
	chiA2 = chiA*chiA;
	chiS2 = chiS*chiS;
	REAL8 chiS3 = chiS2*chiS;
	REAL8 chiA3 = chiA2*chiA;


	//XLAL_PRINT_INFO("****************************** a = %e *********************************\n", a);

	/* Check that deltaM has a reasonable value */
	if (dM2 < 0) {
		XLALPrintError("XLAL Error - %s: eta seems to be > 0.25 - this isn't allowed!\n", __func__);
		XLAL_ERROR(XLAL_EINVAL);
	}
	dM = sqrt(dM2);
	if (m1 < m2) {
		dM = -dM;
	}
	//dM3 = dM2 * dM;

	aDelta = 0.;
	//a value in delta_lm is 0 in both SEOBNRv1 and SEOBNRv2
	a2 = a * a;
	a3 = a2 * a;

	m1Plus3eta = -1. + 3. * eta;
	m1Plus3eta2 = m1Plus3eta * m1Plus3eta;
	m1Plus3eta3 = m1Plus3eta * m1Plus3eta2;

	/* Initialize all coefficients to zero */
	/* This is important, as we will not set some if dM is zero */
	memset(coeffs, 0, sizeof(FacWaveformCoeffs));


	/*
	 * l = 2, Eqs. A8a and A8b for rho, Eq. A15a for f, Eqs. 20 and 21 of
	 * DIN and Eqs. 27a and 27b of PBFRT for delta as well as eqns 28-29 of PBFRT
	 */

	coeffs->delta22vh3 = 7. / 3.;
	coeffs->delta22vh6 = (-4. * aDelta) / 3. + (428. * LAL_PI) / 105.;
	if (SpinAlignedEOBversion == 4)
		{
			coeffs->delta22vh6 =
	-4. / 3. * (dM * chiA + chiS * (1 - 2 * eta)) +
	(428. * LAL_PI) / 105.;
		}
	coeffs->delta22v8 = (20. * aDelta) / 63.;
	coeffs->delta22vh9 = -2203. / 81. + (1712. * LAL_PI * LAL_PI) / 315.;
	coeffs->delta22v5 = -24. * eta;
	coeffs->delta22v6 = 0.0;
	if (SpinAlignedEOBversion == 2 && chiS + chiA * dM / (1. - 2. * eta) > 0.) {
		coeffs->delta22v6 = -540. * eta * (chiS + chiA * dM / (1. - 2. * eta));
//		double chi = (chiS + chiA * dM / (1. - 2. * eta));
//		coeffs->delta22v6 = eta*(1./4.*(1. - 1080.*chi - sqrt((1. - 1080.*chi)*(1. - 1080*chi) + 8.*(13.5 +270.*chi +13.5*chi*chi))));
		}
   if (SpinAlignedEOBversion == 3 /*&& chiS + chiA * dM / (1. - 2. * eta) > 0.*/) {
//				coeffs->delta22v6 = -540. * eta * (chiS + chiA * dM / (1. - 2. * eta));
	   double chi = (chiS + chiA * dM / (1. - 2. * eta));
	   coeffs->delta22v6 = eta*(1./4.*(1. - 1080.*chi - sqrt((1. - 1080.*chi)*(1. - 1080*chi) + 8.*(13.5 +270.*chi +13.5*chi*chi))));
   }
	coeffs->rho22v2 = -43. / 42. + (55. * eta) / 84.;
	coeffs->rho22v3 = (-2. * (chiS + chiA * dM - chiS * eta)) / 3.;
	switch (SpinAlignedEOBversion) {
	case 1:
		coeffs->rho22v4 = -20555. / 10584. + 0.5 * a2
			- (33025. * eta) / 21168. + (19583. * eta2) / 42336.;
		break;
	case 2:
		coeffs->rho22v4 = -20555. / 10584. + 0.5 * (chiS + chiA * dM) * (chiS + chiA * dM)
			- (33025. * eta) / 21168. + (19583. * eta2) / 42336.;
		break;
	case 3:
		coeffs->rho22v4 = -20555. / 10584. + 0.5 * (chiS + chiA * dM) * (chiS + chiA * dM)
		- (33025. * eta) / 21168. + (19583. * eta2) / 42336.;
		break;
	case 4:
		coeffs->rho22v4 =
			-20555. / 10584. + 0.5 * (chiS + chiA * dM) * (chiS + chiA * dM) -
			(33025. * eta) / 21168. + (19583. * eta2) / 42336.;
		break;
	default:
		XLALPrintError("XLAL Error - %s: wrong SpinAlignedEOBversion value, must be 1, 2 or 3!\n", __func__);
		XLAL_ERROR(XLAL_EINVAL);
		break;
	}
	coeffs->rho22v5 = (-34. * a) / 21.;
if (SpinAlignedEOBversion == 4)
			{
				coeffs->rho22v5 =
		(-34. / 21. + 49. * eta / 18. + 209. * eta2 / 126.) * chiS +
		(-34. / 21. - 19. * eta / 42.) * dM * chiA;
			}
	coeffs->rho22v6 = 1556919113. / 122245200. + (89. * a2) / 252. - (48993925. * eta) / 9779616.
		- (6292061. * eta2) / 3259872. + (10620745. * eta3) / 39118464.
		+ (41. * eta * LAL_PI * LAL_PI) / 192.;
	coeffs->rho22v6l = -428. / 105.;
	coeffs->rho22v7 = (18733. * a) / 15876. + a * a2 / 3.;
	if (SpinAlignedEOBversion == 4)
			{
				coeffs->rho22v7 =
		a3 / 3. + chiA * dM * (18733. / 15876. + (50140. * eta) / 3969. +
							 (97865. * eta2) / 63504.) +
		chiS * (18733. / 15876. + (74749. * eta) / 5292. -
			(245717. * eta2) / 63504. + (50803. * eta3) / 63504.);
			}
	switch (SpinAlignedEOBversion) {
	case 1:
		coeffs->rho22v8 = eta * (-5.6 - 117.6 * eta) - 387216563023. / 160190110080. + (18353. * a2) / 21168. - a2 * a2 / 8.;
		break;
	case 2:
		coeffs->rho22v8 = -387216563023. / 160190110080. + (18353. * a2) / 21168. - a2 * a2 / 8.;
		break;
	case 3:
		coeffs->rho22v8 = -387216563023. / 160190110080. + (18353. * a2) / 21168. - a2 * a2 / 8.;
		break;
	case 4:
		coeffs->rho22v8 =
			-387216563023. / 160190110080. + (18353. * a2) / 21168. -
			a2 * a2 / 8.;
		break;
	default:
		XLALPrintError("XLAL Error - %s: wrong SpinAlignedEOBversion value, must be 1, 2 or 3!\n", __func__);
		XLAL_ERROR(XLAL_EINVAL);
		break;
	}
	coeffs->rho22v8l = 9202. / 2205.;
	coeffs->rho22v10 = -16094530514677. / 533967033600.;
	coeffs->rho22v10l = 439877. / 55566.;

	if (debugPK) {
		XLAL_PRINT_INFO("\nPK:: dM, eta, chiS, chiA while renewing hlm coeffs: %e, %e, %e, %e\n",
			   dM, eta, chiS, chiA);
		XLAL_PRINT_INFO("PK:: Renewed rho-lm coeffs: v2 = %.16e, v3 = %.16e, v4 = %.16e, v5 = %.16e\nv6 = %.16e, v6l = %.16e v7 = %.16e v8 = %.16e, v8l = %.16e v10 = %.16e v10l = %.16e\n",
			   coeffs->rho22v2, coeffs->rho22v3, coeffs->rho22v4, coeffs->rho22v5,
			   coeffs->rho22v6, coeffs->rho22v6l, coeffs->rho22v7, coeffs->rho22v8,
			 coeffs->rho22v8l, coeffs->rho22v10, coeffs->rho22v10l);
	}
	if(waveform == 1){
		a = 0;
		a2 = 0;
		a3 = 0;
	}
	coeffs->delta21vh3 = 2. / 3.;
	coeffs->delta21vh6 = (-17. * aDelta) / 35. + (107. * LAL_PI) / 105.;
	coeffs->delta21vh7 = (3. * aDelta * aDelta) / 140.;
	coeffs->delta21vh9 = -272. / 81. + (214. * LAL_PI * LAL_PI) / 315.;
	coeffs->delta21v5 = -493. * eta / 42.;
	coeffs->delta21v7 = 0.0;
	if (dM2) {

		//coeffs->rho21v1 = (-3. * (chiS + chiA / dM)) / (4.);
		coeffs->rho21v1 = 0.0;
		//coeffs->rho21v2 = -59. / 56 - (9. * chiAPlusChiSdM * chiAPlusChiSdM) / (32. * dM2) + (23. * eta) / 84.;
		switch (SpinAlignedEOBversion) {
			case 1:
				coeffs->rho21v2 = -59. / 56 + (23. * eta) / 84. - 9. / 32. * a2;
				coeffs->rho21v3 = 1177. / 672. * a - 27. / 128. * a3;
				break;
			case 2:
				coeffs->rho21v2 = -59. / 56 + (23. * eta) / 84.;
				coeffs->rho21v3 = 0.0;
				break;
			case 3:
				coeffs->rho21v2 = -59. / 56 + (23. * eta) / 84.;
				coeffs->rho21v3 = 0.0;
				break;
			case 4:
				coeffs->rho21v2 = -59. / 56 + (23. * eta) / 84.;
				coeffs->rho21v3 = 0.0;
				break;
			default:
				XLALPrintError("XLAL Error - %s: wrong SpinAlignedEOBversion value, must be 1 or 2!\n", __func__);
				XLAL_ERROR(XLAL_EINVAL);
				break;
		}
		/*
		 * coeffs->rho21v3   = (-567.*chiA*chiA*chiA -
		 * 1701.*chiA*chiA*chiS*dM + chiA*(-4708. + 1701.*chiS*chiS -
		 * 2648.*eta)*(-1. + 4.*eta) + chiS* dM3 *(4708. -
		 * 567.*chiS*chiS + 1816.*eta))/(2688.*dM3);
		 */
		coeffs->rho21v4 = -47009. / 56448. - (865. * a2) / 1792. - (405. * a2 * a2) / 2048. - (10993. * eta) / 14112.
			+ (617. * eta2) / 4704.;
		coeffs->rho21v5 = (-98635. * a) / 75264. + (2031. * a * a2) / 7168. - (1701. * a2 * a3) / 8192.;
		coeffs->rho21v6 = 7613184941. / 2607897600. + (9032393. * a2) / 1806336. + (3897. * a2 * a2) / 16384.
			- (15309. * a3 * a3) / 65536.;
		coeffs->rho21v6l = -107. / 105.;
		coeffs->rho21v7 = (-3859374457. * a) / 1159065600. - (55169. * a3) / 16384.
			+ (18603. * a2 * a3) / 65536. - (72171. * a2 * a2 * a3) / 262144.;
		coeffs->rho21v7l = 107. * a / 140.;
		coeffs->rho21v8 = -1168617463883. / 911303737344.;
		coeffs->rho21v8l = 6313. / 5880.;
		coeffs->rho21v10 = -63735873771463. / 16569158860800.;
		coeffs->rho21v10l = 5029963. / 5927040.;

		coeffs->f21v1 = (-3. * (chiS + chiA / dM)) / (2.);
		switch (SpinAlignedEOBversion) {
		case 1:
			coeffs->f21v3 = 0.0;
			break;
		case 2:
			coeffs->f21v3 = (chiS * dM * (427. + 79. * eta) + chiA * (147. + 280. * dM * dM + 1251. * eta)) / 84. / dM;
			break;
		case 3:
			coeffs->f21v3 = (chiS * dM * (427. + 79. * eta) + chiA * (147. + 280. * dM * dM + 1251. * eta)) / 84. / dM;
			break;
		case 4:
			coeffs->f21v3 =
			(chiS * dM * (427. + 79. * eta) +
				chiA * (147. + 280. * dM * dM + 1251. * eta)) / 84. / dM;
			/*RC: New terms for SEOBNRv4HM, they are put to zero if use_hm == 0 */
			if (waveform == 0)
			{
				coeffs->f21v4 = 0.0;
				coeffs->f21v5 = 0.0;
				coeffs->f21v6 = 0.0;
				coeffs->f21v7c = 0;
			}
			else{
				coeffs->f21v4 = (-3.-2.*eta)*chiA2 + (-3.+ eta/2.)*chiS2 + (-6.+21.*eta/2.)*chiS*chiA/dM;
				coeffs->f21v5 = (3./4.-3.*eta)*chiA3/dM + (-81./16. +1709.*eta/1008. + 613.*eta2/1008.+(9./4.-3*eta)*chiA2)*chiS + 3./4.*chiS3
				+ (-81./16. - 703.*eta2/112. + 8797.*eta/1008.+(9./4. - 6.*eta)*chiS2)*chiA/dM;
				coeffs->f21v6 = (4163./252.-9287.*eta/1008. - 85.*eta2/112.)*chiA2 + (4163./252. - 2633.*eta/1008. + 461.*eta2/1008.)*chiS2 + (4163./126.-1636.*eta/21. + 1088.*eta2/63.)*chiS*chiA/dM;
			}
			/* End new terms for SEOBNRv4HM */
			break;
		default:
			XLALPrintError("XLAL Error - %s: wrong SpinAlignedEOBversion value, must be 1, 2 or 3!\n", __func__);
			XLAL_ERROR(XLAL_EINVAL);
			break;
		}
	} else {
		coeffs->f21v1 = -3. * chiA / 2.; // Odd modes in dm->0 limit, note that there is a typo in Taracchini et.al. discussion after A7 (even->odd)
		switch (SpinAlignedEOBversion) {
		case 1:
			coeffs->f21v3 = 0.0;
			break;
		case 2:
			coeffs->f21v3 = (chiS * dM * (427. + 79. * eta) + chiA * (147. + 280. * dM * dM + 1251. * eta)) / 84.;
			break;
		case 3:
			coeffs->f21v3 = (chiS * dM * (427. + 79. * eta) + chiA * (147. + 280. * dM * dM + 1251. * eta)) / 84.;
			break;
		case 4:
			coeffs->f21v3 =
			(chiS * dM * (427. + 79. * eta) +
				chiA * (147. + 280. * dM * dM + 1251. * eta)) / 84.;
			/* New terms for SEOBNRv4HM, they are put to zero if use_hm == 0 */
			if (waveform == 0)
			{
				coeffs->f21v4 = 0.0;
				coeffs->f21v5 = 0.0;
				coeffs->f21v6 = 0.0;
			}
			else{
				coeffs->f21v4 = (-6+21*eta/2.)*chiS*chiA;
				coeffs->f21v5 = (3./4.-3.*eta)*chiA3  + (-81./16. - 703.*eta2/112. + 8797.*eta/1008. + (9./4. - 6.*eta)*chiS2)*chiA;
				coeffs->f21v6 = (4163./126.-1636.*eta/21. + 1088.*eta2/63.)*chiS*chiA;
				// printf("f21v1 = %.16f\n f21v3 = %.16f\n f21v4 = %.16f\n f21v5 = %.16f\n f21v6 = %.16f\n", coeffs->f21v1, coeffs->f21v3, coeffs->f21v4, coeffs->f21v5, coeffs->f21v6);
			}
			/* End new terms for SEOBNRv4HM */
			break;
		default:
			XLALPrintError("XLAL Error - %s: wrong SpinAlignedEOBversion value, must be 1 or 2!\n", __func__);
			XLAL_ERROR(XLAL_EINVAL);
			break;
		}
	}

	/*
	 * l = 3, Eqs. A9a - A9c for rho, Eqs. A15b and A15c for f, Eqs. 22 -
	 * 24 of DIN and Eqs. 27c - 27e of PBFRT for delta
	 */
	coeffs->delta33vh3 = 13. / 10.;
	coeffs->delta33vh6 = (-81. * aDelta) / 20. + (39. * LAL_PI) / 7.;
	coeffs->delta33vh9 = -227827. / 3000. + (78. * LAL_PI * LAL_PI) / 7.;
	coeffs->delta33v5 = -80897. * eta / 2430.;
	if (dM2) {
		coeffs->rho33v2 = -7. / 6. + (2. * eta) / 3.;
		//coeffs->rho33v3 = (chiS * dM * (-4. + 5. * eta) + chiA * (-4. + 19. * eta)) / (6. * dM);
		coeffs->rho33v3 = 0.0;
		coeffs->rho33v4 = -6719. / 3960. + a2 / 2. - (1861. * eta) / 990. + (149. * eta2) / 330.;
		coeffs->rho33v5 = (-4. * a) / 3.;
		coeffs->rho33v6 = 3203101567. / 227026800. + (5. * a2) / 36.;
		coeffs->rho33v6l = -26. / 7.;
		coeffs->rho33v7 = (5297. * a) / 2970. + a * a2 / 3.;
		coeffs->rho33v8 = -57566572157. / 8562153600.;
		coeffs->rho33v8l = 13. / 3.;
		if(waveform == 1){
			//RC: This terms are in Eq.A6 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
			coeffs->rho33v6 = 3203101567. / 227026800. + (5. * a2) / 36. + (-129509./25740. + 41./192. * LAL_PI*LAL_PI)*eta - 274621./154440.*eta2+ 12011./46332.*eta3;
			coeffs->rho33v10 = -903823148417327./30566888352000.;
			coeffs->rho33v10l = 87347./13860.;
		}

		coeffs->f33v3 = (chiS * dM * (-4. + 5. * eta) + chiA * (-4. + 19. * eta)) / (2. * dM);
		coeffs->f33v4 = 0;
		coeffs->f33v5 = 0;
		coeffs->f33v6 = 0;
		coeffs->f33vh6 = 0;
		if(waveform == 1){
			//RC: This terms are in Eq.A10 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
			coeffs->f33v4 = (3./2. * chiS2 * dM + (3. - 12 * eta) * chiA * chiS + dM * (3./2. -6. * eta) * chiA2)/(dM);
			coeffs->f33v5 = (dM * (241./30. * eta2 + 11./20. * eta + 2./3.) * chiS + (407./30. * eta2 - 593./60. * eta + 2./3.)* chiA)/(dM);
			coeffs->f33v6 = (dM * (6. * eta2 -27. / 2. * eta - 7./ 4.) * chiS2 + (44. * eta2 - 1. * eta - 7./2.) * chiA * chiS + dM * (-12 * eta2 + 11./2. * eta - 7./4.) * chiA2)/dM;
			coeffs->f33vh6 = (dM * (593. / 108. * eta - 81./20.) * chiS + (7339./540. * eta - 81./20.) * chiA)/(dM);
		}
	}
	else {
		coeffs->f33v3 = chiA * 3. / 8.;
		if(waveform == 1){
			//RC: This terms are in Eq.A10 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
			coeffs->f33v4 = ((3. - 12 * eta) * chiA * chiS);
			coeffs->f33v5 = ((407./30. * eta2 - 593./60. * eta + 2./3.)* chiA);
			coeffs->f33v6 = ((44. * eta2 - 1. * eta - 7./2.) * chiA * chiS);
			coeffs->f33vh6 = ((7339./540. * eta - 81./20.) * chiA);
		}
	}

	coeffs->delta32vh3 = (10. + 33. * eta) / (-15. * m1Plus3eta);
	coeffs->delta32vh4 = 4. * aDelta;
	coeffs->delta32vh6 = (-136. * aDelta) / 45. + (52. * LAL_PI) / 21.;
	coeffs->delta32vh9 = -9112. / 405. + (208. * LAL_PI * LAL_PI) / 63.;

	coeffs->rho32v = (4. * chiS * eta) / (-3. * m1Plus3eta);
	/** TODO The term proportional to eta^2 a^2 in coeffs->rho32v2 is wrong, but it was used in the calibration of SEOBNRv2 */
	coeffs->rho32v2 = (-4. * a2 * eta2) / (9. * m1Plus3eta2) + (328. - 1115. * eta
						 + 320. * eta2) / (270. * m1Plus3eta);
	if (SpinAlignedEOBversion == 4) {
		coeffs->rho32v2 = (328. - 1115. * eta +
			320. * eta2) / (270. * m1Plus3eta);
	}
	//coeffs->rho32v3 = (2. * (45. * a * m1Plus3eta3
	//	 - a * eta * (328. - 2099. * eta + 5. * (733. + 20. * a2) * eta2
	//		  - 960. * eta3))) / (405. * m1Plus3eta3);
	coeffs->rho32v3 = 2. / 9. * a;
	coeffs->rho32v4 = a2 / 3. + (-1444528.
			   + 8050045. * eta - 4725605. * eta2 - 20338960. * eta3
			   + 3085640. * eta2 * eta2) / (1603800. * m1Plus3eta2);
	coeffs->rho32v5 = (-2788. * a) / 1215.;
	coeffs->rho32v6 = 5849948554. / 940355325. + (488. * a2) / 405.;
	coeffs->rho32v6l = -104. / 63.;
	coeffs->rho32v8 = -10607269449358. / 3072140846775.;
	coeffs->rho32v8l = 17056. / 8505.;

	if (dM2) {
		coeffs->delta31vh3 = 13. / 30.;
		coeffs->delta31vh6 = (61. * aDelta) / 20. + (13. * LAL_PI) / 21.;
		coeffs->delta31vh7 = (-24. * aDelta * aDelta) / 5.;
		coeffs->delta31vh9 = -227827. / 81000. + (26. * LAL_PI * LAL_PI) / 63.;
		coeffs->delta31v5 = -17. * eta / 10.;

		coeffs->rho31v2 = -13. / 18. - (2. * eta) / 9.;
		//coeffs->rho31v3 = (chiA * (-4. + 11. * eta) + chiS * dM * (-4. + 13. * eta)) / (6. * dM);
		coeffs->rho31v3 = 0.0;
		coeffs->rho31v4 = 101. / 7128.
			- (5. * a2) / 6. - (1685. * eta) / 1782. - (829. * eta2) / 1782.;
		coeffs->rho31v5 = (4. * a) / 9.;
		coeffs->rho31v6 = 11706720301. / 6129723600. - (49. * a2) / 108.;
		coeffs->rho31v6l = -26. / 63.;
		coeffs->rho31v7 = (-2579. * a) / 5346. + a * a2 / 9.;
		coeffs->rho31v8 = 2606097992581. / 4854741091200.;
		coeffs->rho31v8l = 169. / 567.;

		coeffs->f31v3 = (chiA * (-4. + 11. * eta) + chiS * dM * (-4. + 13. * eta)) / (2. * dM);
	} else {
		coeffs->f31v3 = -chiA * 5. / 8.;
	}

	/*
	 * l = 4, Eqs. A10a - A10d for rho, Eq. A15d for f Eqs. 25 - 28 of
	 * DIN and Eqs. 27f - 27i of PBFRT for delta
	 */

	coeffs->delta44vh3 = (112. + 219. * eta) / (-120. * m1Plus3eta);
	coeffs->delta44vh6 = (-464. * aDelta) / 75. + (25136. * LAL_PI) / 3465.;
	coeffs->delta44vh9 = 0.;
  	if(waveform == 1){
		//RC: This terms are in Eq.A15 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
  		coeffs->delta44vh9 = -55144./375. + 201088.*LAL_PI*LAL_PI/10395.;
  	}

	coeffs->rho44v2 = (1614. - 5870. * eta + 2625. * eta2) / (1320. * m1Plus3eta);
	coeffs->rho44v3 = (chiA * (10. - 39. * eta) * dM + chiS * (10. - 41. * eta
					+ 42. * eta2)) / (15. * m1Plus3eta);
	coeffs->rho44v4 = a2 / 2. + (-511573572.
		+ 2338945704. * eta - 313857376. * eta2 - 6733146000. * eta3
		  + 1252563795. * eta2 * eta2) / (317116800. * m1Plus3eta2);
	coeffs->rho44v5 = (-69. * a) / 55.;
	coeffs->rho44v8 = 0.;
 	coeffs->rho44v8l = 0.;
 	coeffs->rho44v10 = 0.;
 	coeffs->rho44v10l = 0;
  	if(waveform == 1){
	  //RC: This terms are in Eq.A8 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
	  coeffs->rho44v4 =
		(-511573572. + 2338945704. * eta - 313857376. * eta2 -
		 6733146000. * eta3 +
		 1252563795. * eta2 * eta2) / (317116800. * m1Plus3eta2)
		+ chiS2/2. + dM*chiS*chiA + dM2*chiA2/2.;
	  coeffs->rho44v5 = chiA*dM*(-8280. + 42716.*eta - 57990.*eta2 + 8955*eta3)/(6600.*m1Plus3eta2)
		+ chiS*(-8280. + 66284.*eta-176418.*eta2+128085.*eta3 + 88650*eta2*eta2)/(6600.*m1Plus3eta2);
	  coeffs->rho44v8 = -172066910136202271./19426955708160000.;
	  coeffs->rho44v8l = 845198./190575.;
	  coeffs->rho44v10 = - 17154485653213713419357./568432724020761600000.;
	  coeffs->rho44v10l = 22324502267./3815311500;
	}
	coeffs->rho44v6 = 16600939332793. / 1098809712000. + (217. * a2) / 3960.;
	coeffs->rho44v6l = -12568. / 3465.;

	if (dM2) {
		coeffs->delta43vh3 = (486. + 4961. * eta) / (810. * (1. - 2. * eta));
		coeffs->delta43vh4 = (11. * aDelta) / 4.;
		coeffs->delta43vh6 = 1571. * LAL_PI / 385.;

		//coeffs->rho43v = (5. * (chiA - chiS * dM) * eta) / (8. * dM * (-1. + 2. * eta));
		coeffs->rho43v = 0.0;
		coeffs->rho43v2 = (222. - 547. * eta + 160. * eta2) / (176. * (-1. + 2. * eta));
		coeffs->rho43v4 = -6894273. / 7047040. + (3. * a2) / 8.;
		coeffs->rho43v5 = (-12113. * a) / 6160.;
		coeffs->rho43v6 = 1664224207351. / 195343948800.;
		coeffs->rho43v6l = -1571. / 770.;

		coeffs->f43v = (5. * (chiA - chiS * dM) * eta) / (2. * dM * (-1. + 2. * eta));
	} else {
		coeffs->f43v = -5. * chiA / 4.;
	}

	coeffs->delta42vh3 = (7. * (1. + 6. * eta)) / (-15. * m1Plus3eta);
	coeffs->delta42vh6 = (212. * aDelta) / 75. + (6284. * LAL_PI) / 3465.;

	coeffs->rho42v2 = (1146. - 3530. * eta + 285. * eta2) / (1320. * m1Plus3eta);
	coeffs->rho42v3 = (chiA * (10. - 21. * eta) * dM + chiS * (10. - 59. * eta
					+ 78. * eta2)) / (15. * m1Plus3eta);
	coeffs->rho42v4 = a2 / 2. + (-114859044. + 295834536. * eta + 1204388696. * eta2 - 3047981160. * eta3
		   - 379526805. * eta2 * eta2) / (317116800. * m1Plus3eta2);
	coeffs->rho42v5 = (-7. * a) / 110.;
	coeffs->rho42v6 = 848238724511. / 219761942400. + (2323. * a2) / 3960.;
	coeffs->rho42v6l = -3142. / 3465.;

	if (dM2) {
		coeffs->delta41vh3 = (2. + 507. * eta) / (10. * (1. - 2. * eta));
		coeffs->delta41vh4 = (11. * aDelta) / 12.;
		coeffs->delta41vh6 = 1571. * LAL_PI / 3465.;

		//coeffs->rho41v = (5. * (chiA - chiS * dM) * eta) / (8. * dM * (-1. + 2. * eta));
		coeffs->rho41v = 0.0;
		coeffs->rho41v2 = (602. - 1385. * eta + 288. * eta2) / (528. * (-1. + 2. * eta));
		coeffs->rho41v4 = -7775491. / 21141120. + (3. * a2) / 8.;
		coeffs->rho41v5 = (-20033. * a) / 55440. - (5 * a * a2) / 6.;
		coeffs->rho41v6 = 1227423222031. / 1758095539200.;
		coeffs->rho41v6l = -1571. / 6930.;

		coeffs->f41v = (5. * (chiA - chiS * dM) * eta) / (2. * dM * (-1. + 2. * eta));
	} else {
		coeffs->f41v = -5. * chiA / 4.;
	}

	/*
	 * l = 5, Eqs. A11a - A11e for rho, Eq. 29 of DIN and Eqs. E1a and
	 * E1b of PBFRT for delta
	 */
	coeffs->delta55vh3 =
		(96875. + 857528. * eta) / (131250. * (1 - 2 * eta));
	coeffs->delta55vh6 = 0;
	coeffs->delta55vh9 = 0;
 	if(waveform == 1){
		//RC: This terms are in Eq.A16 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
		coeffs->delta55vh6 = 3865./429.*LAL_PI;
		coeffs->delta55vh9 = (-7686949127. + 954500400.*LAL_PI*LAL_PI)/31783752.;
 	}
	if (dM2) {

		coeffs->rho55v2 = (487. - 1298. * eta + 512. * eta2) / (390. * (-1. + 2. * eta));
		coeffs->rho55v3 = (-2. * a) / 3.;
		coeffs->rho55v4 = -3353747. / 2129400. + a2 / 2.;
		coeffs->rho55v5 = -241. * a / 195.;
		coeffs->rho55v6 = 0.;
		coeffs->rho55v6l = 0.;
		coeffs->rho55v8 = 0.;
		coeffs->rho55v8l = 0.;
		coeffs->rho55v10 = 0.;
		coeffs->rho55v10l = 0.;
		coeffs->f55v3 = 0.;
		coeffs->f55v4 = 0.;
		coeffs->f55v5c = 0;
		if(waveform == 1){
			//RC: This terms are in Eq.A9 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
			coeffs->rho55v6 = 190606537999247./11957879934000.;
			coeffs->rho55v6l = - 1546./429.;
			coeffs->rho55v8 = - 1213641959949291437./118143853747920000.;
			coeffs->rho55v8l = 376451./83655.;
			coeffs->rho55v10 = -150082616449726042201261./4837990810977324000000.;
			coeffs->rho55v10l = 2592446431./456756300.;

			coeffs->f55v3 = chiA/dM *(10./(3.*(-1.+2.*eta)) - 70.*eta/(3.*(-1.+2.*eta)) + 110.*eta2/(3.*(-1.+2.*eta)) ) +
				chiS*(10./(3.*(-1.+2.*eta)) -10.*eta/(-1.+2.*eta) + 10*eta2/(-1.+2.*eta));
			coeffs->f55v4 = chiS2*(-5./(2.*(-1.+2.*eta)) + 5.*eta/(-1.+2.*eta)) +
				chiA*chiS/dM *(-5./(-1.+2.*eta) + 30.*eta/(-1.+2.*eta) - 40.*eta2/(-1.+2.*eta)) +
				chiA2*(-5./(2.*(-1.+2.*eta)) + 15.*eta/(-1.+2.*eta) - 20.*eta2/(-1.+2.*eta));
			coeffs->f55v5c = 0; //RC: this is the calibration parameter which is initially set to 0.
	  }
	}
	else{
		coeffs->f55v3 = 0;
		coeffs->f55v4 = 0;
		coeffs->f55v5c = 0;
		if(waveform == 1){
			//RC: This terms are in Eq.A12 in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.98.084028 [arXiv:1803.10701]
			coeffs->f55v3 = chiA *(10./(3.*(-1.+2.*eta)) - 70.*eta/(3.*(-1.+2.*eta)) + 110.*eta2/(3.*(-1.+2.*eta)) );
			coeffs->f55v4 = chiA*chiS *(-5./(-1.+2.*eta) + 30.*eta/(-1.+2.*eta) - 40.*eta2/(-1.+2.*eta));
			coeffs->f55v5c = 0;
		}
	}
	coeffs->delta54vh3 = 8. / 15.;
	coeffs->delta54vh4 = 12. * aDelta / 5.;

	coeffs->rho54v2 = (-17448. + 96019. * eta - 127610. * eta2
		  + 33320. * eta3) / (13650. * (1. - 5. * eta + 5. * eta2));
	coeffs->rho54v3 = (-2. * a) / 15.;
	coeffs->rho54v4 = -16213384. / 15526875. + (2. * a2) / 5.;

	if (dM2) {
		coeffs->delta53vh3 = 31. / 70.;

		coeffs->rho53v2 = (375. - 850. * eta + 176. * eta2) / (390. * (-1. + 2. * eta));
		coeffs->rho53v3 = (-2. * a) / 3.;
		coeffs->rho53v4 = -410833. / 709800. + a2 / 2.;
		coeffs->rho53v5 = -103. * a / 325.;
	}
	coeffs->delta52vh3 = 4. / 15.;
	coeffs->delta52vh4 = 6. * aDelta / 5.;

	coeffs->rho52v2 = (-15828. + 84679. * eta - 104930. * eta2
		  + 21980. * eta3) / (13650. * (1. - 5. * eta + 5. * eta2));
	coeffs->rho52v3 = (-2. * a) / 15.;
	coeffs->rho52v4 = -7187914. / 15526875. + (2. * a2) / 5.;

	if (dM2) {
		coeffs->delta51vh3 = 31. / 210.;

		coeffs->rho51v2 = (319. - 626. * eta + 8. * eta2) / (390. * (-1. + 2. * eta));
		coeffs->rho51v3 = (-2. * a) / 3.;
		coeffs->rho51v4 = -31877. / 304200. + a2 / 2.;
		coeffs->rho51v5 = 139. * a / 975.;
	}
	/*
	 * l = 6, Eqs. A12a - A12f for rho, Eqs. E1c and E1d of PBFRT for
	 * delta
	 */

	coeffs->delta66vh3 = 43. / 70.;

	coeffs->rho66v2 = (-106. + 602. * eta - 861. * eta2
			   + 273. * eta3) / (84. * (1. - 5. * eta + 5. * eta2));
	coeffs->rho66v3 = (-2. * a) / 3.;
	coeffs->rho66v4 = -1025435. / 659736. + a2 / 2.;

	if (dM2) {
		coeffs->delta65vh3 = 10. / 21.;

		coeffs->rho65v2 = (-185. + 838. * eta - 910. * eta2
				+ 220. * eta3) / (144. * (dM2 + 3. * eta2));
		coeffs->rho65v3 = -2. * a / 9.;
	}
	coeffs->delta64vh3 = 43. / 105.;

	coeffs->rho64v2 = (-86. + 462. * eta - 581. * eta2
			   + 133. * eta3) / (84. * (1. - 5. * eta + 5. * eta2));
	coeffs->rho64v3 = (-2. * a) / 3.;
	coeffs->rho64v4 = -476887. / 659736. + a2 / 2.;

	if (dM2) {
		coeffs->delta63vh3 = 2. / 7.;

		coeffs->rho63v2 = (-169. + 742. * eta - 750. * eta2
				+ 156. * eta3) / (144. * (dM2 + 3. * eta2));
		coeffs->rho63v3 = -2. * a / 9.;
	}
	coeffs->delta62vh3 = 43. / 210.;

	coeffs->rho62v2 = (-74. + 378. * eta - 413. * eta2
			+ 49. * eta3) / (84. * (1. - 5. * eta + 5. * eta2));
	coeffs->rho62v3 = (-2. * a) / 3.;
	coeffs->rho62v4 = -817991. / 3298680. + a2 / 2.;

	if (dM2) {
		coeffs->delta61vh3 = 2. / 21.;

		coeffs->rho61v2 = (-161. + 694. * eta - 670. * eta2
				+ 124. * eta3) / (144. * (dM2 + 3. * eta2));
		coeffs->rho61v3 = -2. * a / 9.;
	}
	/*
	 * l = 7, Eqs. A13a - A13g for rho, Eqs. E1e and E1f of PBFRT for
	 * delta
	 */
	if (dM2) {
		coeffs->delta77vh3 = 19. / 36.;

		coeffs->rho77v2 = (-906. + 4246. * eta - 4963. * eta2
				   + 1380. * eta3) / (714. * (dM2 + 3. * eta2));
		coeffs->rho77v3 = -2. * a / 3.;
	}
	coeffs->rho76v2 = (2144. - 16185. * eta + 37828. * eta2 - 29351. * eta3
		 + 6104. * eta2 * eta2) / (1666. * (-1 + 7 * eta - 14 * eta2
							+ 7 * eta3));

	if (dM2) {
		coeffs->delta75vh3 = 95. / 252.;

		coeffs->rho75v2 = (-762. + 3382. * eta - 3523. * eta2
				+ 804. * eta3) / (714. * (dM2 + 3. * eta2));
		coeffs->rho75v3 = -2. * a / 3.;
	}
	coeffs->rho74v2 = (17756. - 131805. * eta + 298872. * eta2 - 217959. * eta3
		+ 41076. * eta2 * eta2) / (14994. * (-1. + 7. * eta - 14. * eta2
						 + 7. * eta3));

	if (dM2) {
		coeffs->delta73vh3 = 19. / 84.;

		coeffs->rho73v2 = (-666. + 2806. * eta - 2563. * eta2
				+ 420. * eta3) / (714. * (dM2 + 3. * eta2));
		coeffs->rho73v3 = -2. * a / 3.;
	}
	coeffs->rho72v2 = (16832. - 123489. * eta + 273924. * eta2 - 190239. * eta3
		+ 32760. * eta2 * eta2) / (14994. * (-1. + 7. * eta - 14. * eta2
						 + 7. * eta3));

	if (dM2) {
		coeffs->delta71vh3 = 19. / 252.;

		coeffs->rho71v2 = (-618. + 2518. * eta - 2083. * eta2
				+ 228. * eta3) / (714. * (dM2 + 3. * eta2));
		coeffs->rho71v3 = -2. * a / 3.;
	}
	/* l = 8, Eqs. A14a - A14h */

	coeffs->rho88v2 = (3482. - 26778. * eta + 64659. * eta2 - 53445. * eta3
		 + 12243. * eta2 * eta2) / (2736. * (-1. + 7. * eta - 14. * eta2
						 + 7. * eta3));

	if (dM2) {
		coeffs->rho87v2 = (23478. - 154099. * eta + 309498. * eta2 - 207550. * eta3
		+ 38920 * eta2 * eta2) / (18240. * (-1 + 6 * eta - 10 * eta2
							+ 4 * eta3));
	}
	coeffs->rho86v2 = (1002. - 7498. * eta + 17269. * eta2 - 13055. * eta3
		   + 2653. * eta2 * eta2) / (912. * (-1. + 7. * eta - 14. * eta2
						 + 7. * eta3));

	if (dM2) {
		coeffs->rho85v2 = (4350. - 28055. * eta + 54642. * eta2 - 34598. * eta3
				   + 6056. * eta2 * eta2) / (3648. * (-1. + 6. * eta - 10. * eta2
								  + 4. * eta3));
	}
	coeffs->rho84v2 = (2666. - 19434. * eta + 42627. * eta2 - 28965. * eta3
		  + 4899. * eta2 * eta2) / (2736. * (-1. + 7. * eta - 14. * eta2
						 + 7. * eta3));

	if (dM2) {
		coeffs->rho83v2 = (20598. - 131059. * eta + 249018. * eta2 - 149950. * eta3
				   + 24520. * eta2 * eta2) / (18240. * (-1. + 6. * eta - 10. * eta2
								  + 4. * eta3));
	}
	coeffs->rho82v2 = (2462. - 17598. * eta + 37119. * eta2 - 22845. * eta3
		  + 3063. * eta2 * eta2) / (2736. * (-1. + 7. * eta - 14. * eta2
						 + 7. * eta3));

	if (dM2) {
		coeffs->rho81v2 = (20022. - 126451. * eta + 236922. * eta2 - 138430. * eta3
				   + 21640. * eta2 * eta2) / (18240. * (-1. + 6. * eta - 10. * eta2
								  + 4. * eta3));
	}
	/* All relevant coefficients should be set, so we return */

	return XLAL_SUCCESS;
}

#endif
