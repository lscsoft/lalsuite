/*
 * Copyright (C) 2010 Craig Robinson, Yi Pan
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option)
 * any later version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * with program; see the file COPYING. If not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */


/**
 * \author Craig Robinson, Yi Pan
 *
 * \brief In newer versions of the EOBNR approximant, we
 * do not have an analytic expression for the derivative of the waveform.
 * As such, it is necessary to calculate the derivatives numerically. This
 * function provides the means to do just that.
 *
 */

#ifndef _LALSIMIMRSPINEOBHCAPNUMERICALDERIVATIVE_C
#define _LALSIMIMRSPINEOBHCAPNUMERICALDERIVATIVE_C

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/* #include "LALSimIMRSpinEOBHcapNumericalDerivative.h" */

#include <unistd.h>
#include <math.h>
#include <gsl/gsl_deriv.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

#include "LALSimIMRSpinEOB.h"

#include "LALSimIMRSpinEOBAuxFuncs.c"
#include "LALSimIMRSpinEOBHamiltonian.c"
#include "LALSimIMRSpinEOBFactorizedFlux.c"
#include "LALSimIMRSpinEOBFactorizedWaveformCoefficients.c"
#include "LALSimIMREOBFactorizedWaveform.c"

// int		UsePrec = 1;

/*------------------------------------------------------------------------------------------
 *
 *          Prototypes of functions defined in this code.
 *
 *------------------------------------------------------------------------------------------
 */

static UNUSED REAL8 GSLSpinHamiltonianWrapper(double x, void *params);

static double	GSLSpinHamiltonianWrapperV2(double x, void *params);

static int 
XLALSpinHcapNumericalDerivative(
				double UNUSED t,	/**<< UNUSED */
				const REAL8 values[],	/**<< Dynamical variables */
				REAL8 dvalues[],	/**<< Time derivatives of variables (returned) */
				void *funcParams	/**<< EOB parameters */
);

	static UNUSED REAL8 XLALSpinHcapNumDerivWRTParam(
			    		const		INT4	paramIdx,	/**<< Index of the parameters */
			    		const		REAL8	values[],	/**<< Dynamical variables */
				   		SpinEOBParams * funcParams	/**<< EOB Parameters */
);

	static int	XLALSpinHcapNumericalDerivativeNoFlux(
					   		double	UNUSED	t,	/**<< UNUSED */
			    		const		REAL8	values[],	/**<< Dynamical variables */
				   		REAL8		dvalues[],	/**<< Time derivatives of variables (returned) */
				     		void         *funcParams	/**<< EOB parameters */
);



/*------------------------------------------------------------------------------------------
 *
 *          Defintions of functions.
 *
 *------------------------------------------------------------------------------------------
 */

/**
 * Function to calculate numerical derivatives of the spin EOB Hamiltonian,
 * which correspond to time derivatives of the dynamical variables in conservative dynamcis.
 * All derivatives, including those on two terms of the orbital phase, are returned together.
 * The derivatives are combined with energy flux to give right hand side of the ODEs
 * of a generic spin EOB model, as decribed in Eqs. 21, 22, 26 and 27 of
 * Pan et al. PRD 81, 084041 (2010)
 * This function is not used by the spin-aligned SEOBNRv1 model.
 */
	static int	XLALSpinHcapNumericalDerivative(
					     		double	UNUSED	t,	/**<< UNUSED */
			      		const		REAL8	values[],	/**<< Dynamical variables */
				     		REAL8		dvalues[],	/**<< Time derivatives of variables (returned) */
				       		void         *funcParams	/**<< EOB parameters */
)
{
	int		debugPK = 0;
	static const REAL8 STEP_SIZE = 1.0e-4;

	static const INT4 lMax = 8;

	HcapDerivParams	params;

	/* Since we take numerical derivatives wrt dynamical variables */
	/* but we want them wrt time, we use this temporary vector in  */
	/* the conversion */
	REAL8		tmpDValues[14];

	REAL8		H;
	//Hamiltonian
		REAL8 flux;

	gsl_function	F;
	INT4		gslStatus;
	UINT4		SpinAlignedEOBversion;

	UINT4		i       , j, k, l;

	REAL8Vector	rVec, pVec;
	REAL8		rData    [3], pData[3];

	/* We need r, phi, pr, pPhi to calculate the flux */
	REAL8		r;
	REAL8Vector	polarDynamics, cartDynamics;
	REAL8		polData  [4];

	REAL8		mass1   , mass2, eta;
	REAL8 UNUSED	rrTerm2, pDotS1, pDotS2;
	REAL8Vector	s1 , s2, s1norm, s2norm, sKerr, sStar;
	REAL8		s1Data   [3], s2Data[3], s1DataNorm[3], s2DataNorm[3];
	REAL8		sKerrData[3], sStarData[3];
	REAL8 /* magS1, magS2, */ chiS, chiA, a, tplspin;
	REAL8 UNUSED	s1dotL, s2dotL;
	REAL8 UNUSED	rcrossrDot[3], rcrossrDotMag, s1dotLN, s2dotLN;


	/* Orbital angular momentum */
	REAL8		Lx      , Ly, Lz, magL;
	REAL8		Lhatx   , Lhaty, Lhatz;
	REAL8		dLx     , dLy, dLz;
	REAL8		dLhatx  , dLhaty, dMagL;

	REAL8		alphadotcosi;

	REAL8		rCrossV_x, rCrossV_y, rCrossV_z, omega;

	/* The error in a derivative as measured by GSL */
	REAL8		absErr;

	REAL8		tmpP     [3], rMag, rMag2, prT;
	REAL8		u       , u2, u3, u4, u5, w2, a2;
	REAL8		D       , m1PlusetaKK, bulk, logTerms, deltaU, deltaT, deltaR;
	REAL8 UNUSED	eobD_r, deltaU_u, deltaU_r, deltaT_r;
	REAL8		dcsi    , csi;

	REAL8		tmpValues[12];
	REAL8		Tmatrix  [3][3], invTmatrix[3][3], dTijdXk[3][3][3];
	REAL8		tmpPdotT1[3], tmpPdotT2[3], tmpPdotT3[3];
	//3 terms of Eq.A5

	/* NQC coefficients container */
		EOBNonQCCoeffs * nqcCoeffs = NULL;

	/* Set up pointers for GSL */
	params.values = values;
	params.params = (SpinEOBParams *) funcParams;
	nqcCoeffs = params.params->nqcCoeffs;

	F.function = &GSLSpinHamiltonianWrapper;
	F.params = &params;

	mass1 = params.params->eobParams->m1;
	mass2 = params.params->eobParams->m2;
	eta = params.params->eobParams->eta;
	SpinAlignedEOBversion = params.params->seobCoeffs->SpinAlignedEOBversion;
	SpinEOBHCoeffs *coeffs = (SpinEOBHCoeffs *) params.params->seobCoeffs;

	/*
	 * For precessing binaries, the effective spin of the Kerr background
	 * evolves with time. The coefficients used to compute the
	 * Hamiltonian depend on the Kerr spin, and hence need to be updated
	 * for the current spin values
	 */

	/*
	 * Set the position/momenta vectors to point to the appropriate
	 * things
	 */
	rVec.length = pVec.length = 3;
	rVec.data = rData;
	pVec.data = pData;
	memcpy(rData, values, sizeof(rData));
	memcpy(pData, values + 3, sizeof(pData));

	/*
	 * We need to re-calculate the parameters at each step as precessing
	 * spins will not be constant
	 */
	/* TODO: Modify so that only spin terms get re-calculated */

	/*
	 * We cannot point to the values vector directly as it leads to a
	 * warning
	 */
	s1.length = s2.length = s1norm.length = s2norm.length = 3;
	s1.data = s1Data;
	s2.data = s2Data;
	s1norm.data = s1DataNorm;
	s2norm.data = s2DataNorm;

	memcpy(s1Data, values + 6, 3 * sizeof(REAL8));
	memcpy(s2Data, values + 9, 3 * sizeof(REAL8));
	memcpy(s1DataNorm, values + 6, 3 * sizeof(REAL8));
	memcpy(s2DataNorm, values + 9, 3 * sizeof(REAL8));

	for (i = 0; i < 3; i++) {
		s1Data[i] *= (mass1 + mass2) * (mass1 + mass2);
		s2Data[i] *= (mass1 + mass2) * (mass1 + mass2);
	}

	sKerr.length = 3;
	sKerr.data = sKerrData;
	XLALSimIMRSpinEOBCalculateSigmaKerr(&sKerr, mass1, mass2, &s1, &s2);

	sStar.length = 3;
	sStar.data = sStarData;
	XLALSimIMRSpinEOBCalculateSigmaStar(&sStar, mass1, mass2, &s1, &s2);

	a = sqrt(sKerr.data[0] * sKerr.data[0] + sKerr.data[1] * sKerr.data[1]
		 + sKerr.data[2] * sKerr.data[2]);

	///*set the tortoise flag to 2 * /
		//INT4 oldTortoise = params.params->tortoise;
	//params.params->tortoise = 2;

	/* Convert momenta to p */
	rMag = sqrt(rData[0] * rData[0] + rData[1] * rData[1] + rData[2] * rData[2]);
	prT = pData[0] * (rData[0] / rMag) + pData[1] * (rData[1] / rMag)
		+ pData[2] * (rData[2] / rMag);

	rMag2 = rMag * rMag;
	u = 1. / rMag;
	u2 = u * u;
	u3 = u2 * u;
	u4 = u2 * u2;
	u5 = u4 * u;
	a2 = a * a;
	w2 = rMag2 + a2;
	/* Eq. 5.83 of BB1, inverse */
	D = 1. + log(1. + 6. * eta * u2 + 2. * (26. - 3. * eta) * eta * u3);
	eobD_r = (u2 / (D * D)) * (12. * eta * u + 6. * (26. - 3. * eta) * eta * u2) / (1.
			+ 6. * eta * u2 + 2. * (26. - 3. * eta) * eta * u3);
	m1PlusetaKK = -1. + eta * coeffs->KK;
	/* Eq. 5.75 of BB1 */
	bulk = 1. / (m1PlusetaKK * m1PlusetaKK) + (2. * u) / m1PlusetaKK + a2 * u2;
	/* Eq. 5.73 of BB1 */
	logTerms = 1. + eta * coeffs->k0 + eta * log(1. + coeffs->k1 * u
		       + coeffs->k2 * u2 + coeffs->k3 * u3 + coeffs->k4 * u4
			     + coeffs->k5 * u5 + coeffs->k5l * u5 * log(u));
	/* Eq. 5.73 of BB1 */
	deltaU = bulk * logTerms;
	/* Eq. 5.71 of BB1 */
	deltaT = rMag2 * deltaU;
	/* ddeltaU/du */
	deltaU_u = 2. * (1. / m1PlusetaKK + a2 * u) * logTerms +
		bulk * (eta * (coeffs->k1 + u * (2. * coeffs->k2 + u * (3. * coeffs->k3
									+ u * (4. * coeffs->k4 + 5. * (coeffs->k5 + coeffs->k5l * log(u)) * u)))))
		/ (1. + coeffs->k1 * u + coeffs->k2 * u2 + coeffs->k3 * u3
	      + coeffs->k4 * u4 + (coeffs->k5 + coeffs->k5l * log(u)) * u5);
	deltaU_r = -u2 * deltaU_u;
	/* Eq. 5.38 of BB1 */
	deltaR = deltaT * D;
	if (params.params->tortoise)
		csi = sqrt(deltaT * deltaR) / w2;	/* Eq. 28 of Pan et al.
							 * PRD 81, 084041 (2010) */
	else
		csi = 1.0;

	for (i = 0; i < 3; i++) {
		tmpP[i] = pData[i] - (rData[i] / rMag) * prT * (csi - 1.) / csi;
	}

	if (debugPK) {
		printf("csi = %.12e\n", csi);
		for (i = 0; i < 3; i++)
			printf("p,p*: %.12e\t%.12e\n", pData[i], tmpP[i]);
	}
	/*
	 * Calculate the T-matrix, required to convert P from tortoise to
	 * non-tortoise coordinates, and/or vice-versa. This is given
	 * explicitly in Eq. A3 of 0912.3466
	 */
	for (i = 0; i < 3; i++)
		for (j = 0; j <= i; j++) {
			Tmatrix[i][j] = Tmatrix[j][i] = (rData[i] * rData[j] / rMag2)
				* (csi - 1.);

			invTmatrix[i][j] = invTmatrix[j][i] =
				-(csi - 1) / csi * (rData[i] * rData[j] / rMag2);

			if (i == j) {
				Tmatrix[i][j]++;
				invTmatrix[i][j]++;
			}
		}

	dcsi = csi * (2. / rMag + deltaU_r / deltaU) + csi * csi * csi
		/ (2. * rMag2 * rMag2 * deltaU * deltaU) * (rMag * (-4. * w2) / D - eobD_r * (w2 * w2));

	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			for (k = 0; k < 3; k++) {
				dTijdXk[i][j][k] =
					(rData[i] * XLALKronecker(j, k) + XLALKronecker(i, k) * rData[j])
					* (csi - 1.) / rMag2
					+ rData[i] * rData[j] * rData[k] / rMag2 / rMag * (-2. / rMag * (csi - 1.) + dcsi);
			}

	//Print out the T - matrix for comparison
		if (debugPK) {
			printf("\nT-Matrix:\n");
			for (i = 0; i < 3; i++)
				printf("%le\t%le\t%le\n", Tmatrix[i][0], Tmatrix[i][1], Tmatrix[i][2]);

			for (i = 0; i < 3; i++) {
				printf("dT[%d][j]dX[k]:\n", i);
				for (j = 0; j < 3; j++)
					printf("%.12e\t%.12e\t%.12e\n", dTijdXk[i][j][0],
					dTijdXk[i][j][1], dTijdXk[i][j][2]);
				printf("\n");
			}
		}
	/* Now calculate derivatives w.r.t. each parameter */
	for (i = 0; i < 12; i++) {
		params.varyParam = i;
		if (i >= 6 && i < 9) {
			params.params->seobCoeffs->updateHCoeffs = 1;
			XLAL_CALLGSL(gslStatus = gsl_deriv_central(&F, values[i],
			STEP_SIZE * mass1 * mass1, &tmpDValues[i], &absErr));
		} else if (i >= 9) {
			params.params->seobCoeffs->updateHCoeffs = 1;
			XLAL_CALLGSL(gslStatus = gsl_deriv_central(&F, values[i],
			STEP_SIZE * mass2 * mass2, &tmpDValues[i], &absErr));
		} else if (i < 3) {
			params.params->tortoise = 2;
			memcpy(tmpValues, params.values, sizeof(tmpValues));
			tmpValues[3] = tmpP[0];
			tmpValues[4] = tmpP[1];
			tmpValues[5] = tmpP[2];
			params.values = tmpValues;

			XLAL_CALLGSL(gslStatus = gsl_deriv_central(&F, values[i],
				       STEP_SIZE, &tmpDValues[i], &absErr));

			params.values = values;
			params.params->tortoise = 1;
		} else {
			params.params->seobCoeffs->updateHCoeffs = 1;
			XLAL_CALLGSL(gslStatus = gsl_deriv_central(&F, values[i],
				       STEP_SIZE, &tmpDValues[i], &absErr));
		}
		if (gslStatus != GSL_SUCCESS) {
			XLALPrintError("XLAL Error %s - Failure in GSL function\n", __func__);
			XLAL_ERROR(XLAL_EFUNC);
		}
	}

	/* Now make the conversion */
	/* rVectorDot */
	for (i = 0; i < 3; i++)
		for (j = 0, dvalues[i] = 0.; j < 3; j++)
			dvalues[i] += tmpDValues[j + 3] * Tmatrix[i][j];

	//dvalues[0] = tmpDValues[3] * Tmatrix[0][0] + tmpDValues[4] * Tmatrix[0][1]
		// +tmpDValues[5] * Tmatrix[0][2];
	//dvalues[1] = tmpDValues[3] * Tmatrix[1][0] + tmpDValues[4] * Tmatrix[1][1]
		// +tmpDValues[5] * Tmatrix[1][2];
	//dvalues[2] = tmpDValues[3] * Tmatrix[2][0] + tmpDValues[4] * Tmatrix[2][1]
		// +tmpDValues[5] * Tmatrix[2][2];

	/* Calculate the orbital angular momentum */
	Lx = values[1] * values[5] - values[2] * values[4];
	Ly = values[2] * values[3] - values[0] * values[5];
	Lz = values[0] * values[4] - values[1] * values[3];

	magL = sqrt(Lx * Lx + Ly * Ly + Lz * Lz);

	Lhatx = Lx / magL;
	Lhaty = Ly / magL;
	Lhatz = Lz / magL;

	/* Calculate the polar data */
	polarDynamics.length = 4;
	polarDynamics.data = polData;

	r = polData[0] = sqrt(values[0] * values[0] + values[1] * values[1]
			      + values[2] * values[2]);
	polData[1] = 0;
	polData[2] = (values[0] * values[3] + values[1] * values[4]
		      + values[2] * values[5]) / polData[0];
	polData[3] = magL;


	/* Compute \vec{S_i} \dot \vec{L}	 */
	s1dotL = (s1Data[0] * Lhatx + s1Data[1] * Lhaty + s1Data[2] * Lhatz)
		/ (mass1 * mass1);
	s2dotL = (s2Data[0] * Lhatx + s2Data[1] * Lhaty + s2Data[2] * Lhatz)
		/ (mass2 * mass2);

	/* Now calculate rCrossRdot and omega */
	rCrossV_x = values[1] * dvalues[2] - values[2] * dvalues[1];
	rCrossV_y = values[2] * dvalues[0] - values[0] * dvalues[2];
	rCrossV_z = values[0] * dvalues[1] - values[1] * dvalues[0];

	omega = sqrt(rCrossV_x * rCrossV_x + rCrossV_y * rCrossV_y + rCrossV_z * rCrossV_z) / (r * r);

	/*
	 * Compute \vec{L_N} = \vec{r} \times \.{\vec{r}}, \vec{S_i} \dot
	 * \vec{L_N} and chiS and chiA
	 */
	/*
	 * THIS IS NOT QUITE RIGHT. ITS rcrossr*Dot INSTEAD, as it is
	 * \partial H/\partial p_*
	 */
	/*
	 * rcrossrDot[0] = values[1]*tmpDValues[5] - values[2]*tmpDValues[4];
	 * rcrossrDot[1] = values[2]*tmpDValues[3] - values[0]*tmpDValues[5];
	 * rcrossrDot[2] = values[0]*tmpDValues[4] - values[1]*tmpDValues[3];
	 * rcrossrDotMag = sqrt( rcrossrDot[0]*rcrossrDot[0] +
	 * rcrossrDot[1]*rcrossrDot[1]	+ rcrossrDot[2]*rcrossrDot[2] );
	 * 
	 * rcrossrDot[0] /= rcrossrDotMag; rcrossrDot[1] /= rcrossrDotMag;
	 * rcrossrDot[2] /= rcrossrDotMag;
	 */

	s1dotLN = (s1Data[0] * rCrossV_x + s1Data[1] * rCrossV_y + s1Data[2] * rCrossV_z) /
		(r * r * omega * mass1 * mass1);
	s2dotLN = (s2Data[0] * rCrossV_x + s2Data[1] * rCrossV_y + s2Data[2] * rCrossV_z) /
		(r * r * omega * mass2 * mass2);

	chiS = 0.5 * (s1dotLN + s2dotLN);
	chiA = 0.5 * (s1dotLN - s2dotLN);

	if (debugPK) {
		printf("chiS = %.12e, chiA = %.12e\n", chiS, chiA);
		fflush(NULL);
	}
	REAL8		sscaling1 = 1;
	//(mass1 + mass2) * (mass1 + mass2) / (mass1 * mass1);
	REAL8		sscaling2 = 1;
	//(mass1 + mass2) * (mass1 + mass2) / (mass2 * mass2);

	if (debugPK) {
		printf("Computing derivatives for values\n%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n\n",
		    (double)values[0], (double)values[1], (double)values[2],
		    (double)values[3], (double)values[4], (double)values[5],
		(double)sscaling1 * values[6], (double)sscaling1 * values[7],
		(double)sscaling1 * values[8], (double)sscaling2 * values[9],
		       (double)sscaling2 * values[10], (double)sscaling2 * values[11]);
		printf("tmpDvalues\n%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t\n",
		       (double)tmpDValues[0], (double)tmpDValues[1], (double)tmpDValues[2],
		       (double)tmpDValues[3], (double)tmpDValues[4], (double)tmpDValues[5],
		       (double)tmpDValues[6], (double)tmpDValues[7], (double)tmpDValues[8],
		       (double)tmpDValues[9], (double)tmpDValues[10], (double)tmpDValues[11]);
	}
	/*
	 * Compute the test-particle limit spin of the deformed-Kerr
	 * background
	 */
	/* TODO: Check this is actually the way it works in LAL */
	switch (SpinAlignedEOBversion) {
	case 1:
		tplspin = 0.0;
		break;
	case 2:
		tplspin = (1. - 2. * eta) * chiS + (mass1 - mass2) / (mass1 + mass2) * chiA;
		break;
	default:
		XLALPrintError("XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
		XLAL_ERROR(XLAL_EINVAL);
		break;
	}

	for (i = 0; i < 3; i++) {
		params.params->s1Vec->data[i] = s1norm.data[i];
		params.params->s2Vec->data[i] = s2norm.data[i];
		params.params->sigmaStar->data[i] = sStar.data[i];
		params.params->sigmaKerr->data[i] = sKerr.data[i];
	}

	//params.params->s1Vec = &s1norm;
	//params.params->s2Vec = &s2norm;
	//params.params->sigmaStar = &sStar;
	//params.params->sigmaKerr = &sKerr;
	params.params->a = a;

	XLALSimIMREOBCalcSpinFacWaveformCoefficients(
	      params.params->eobParams->hCoeffs, mass1, mass2, eta, tplspin,
					 chiS, chiA, SpinAlignedEOBversion);
	XLALSimIMRCalculateSpinEOBHCoeffs(params.params->seobCoeffs, eta, a,
					  SpinAlignedEOBversion);

	H = XLALSimIMRSpinEOBHamiltonian(eta, &rVec, &pVec, &s1norm, &s2norm,
	&sKerr, &sStar, params.params->tortoise, params.params->seobCoeffs);
	H = H * (mass1 + mass2);

	/* Now we have the ingredients to compute the flux */
	memcpy(tmpValues, values, 12 * sizeof(REAL8));
	cartDynamics.data = tmpValues;
	if (debugPK) {
		printf("params.params->a = %.12e, %.12e\n", a, params.params->a);
		fflush(NULL);
	}
	flux = XLALInspiralPrecSpinFactorizedFlux(&polarDynamics, &cartDynamics,
						  nqcCoeffs, omega, params.params, H / (mass1 + mass2), lMax, SpinAlignedEOBversion);

	/*
	 * Looking at the non-spinning model, I think we need to divide the
	 * flux by eta
	 */
	flux = flux / eta;

	pDotS1 = pData[0] * s1Data[0] + pData[1] * s1Data[1] + pData[2] * s1Data[2];
	pDotS2 = pData[0] * s2Data[0] + pData[1] * s2Data[1] + pData[2] * s2Data[2];
	rrTerm2 = 8. / 15. * eta * eta * pow(omega, 8. / 3.) / (magL * magL * r) * ((61. + 48. * mass2 / mass1) * pDotS1 + (61. + 48. * mass1 / mass2) * pDotS2);

	if (debugPK) {
		printf("omega = %.12e \n flux = %.12e \n Lmag = %.12e\n", omega, flux, magL);
		printf("rrForce = %.12e %.12e %.12e\n", -flux * values[3] / (omega * magL), -flux * values[4] / (omega * magL), -flux * values[5] / (omega * magL));
	}
	/* Now pDot */
	/* Compute the first and second terms in eq. A5 of 0912.3466 */
	for (i = 0; i < 3; i++) {
		for (j = 0, tmpPdotT1[i] = 0.; j < 3; j++)
			tmpPdotT1[i] += -tmpDValues[j] * Tmatrix[i][j];

		tmpPdotT2[i] = -flux * values[i + 3] / (omega * magL);
	}

	/* Compute the third term in eq. A5 */
	REAL8		tmpPdotT3T11[3][3][3], tmpPdotT3T12[3][3], tmpPdotT3T2[3];

	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			for (l = 0; l < 3; l++)
				for (k = 0, tmpPdotT3T11[i][j][l] = 0.; k < 3; k++)
					tmpPdotT3T11[i][j][l] += dTijdXk[i][k][j] * invTmatrix[k][l];

	if (debugPK) {
		for (i = 0; i < 3; i++)
			for (j = 0; j < 1; j++)
				for (l = 0; l < 3; l++) {
					double		sum = 0;
					for (k = 0; k < 3; k++)
						sum += dTijdXk[i][k][j] * invTmatrix[k][l];
					printf("\n sum[%d][%d][%d] = %.12e", i, j, l, sum);

				}
		printf("\n\n Printing dTdX * Tinverse:\n");
		for (l = 0; l < 3; l++) {
			for (i = 0; i < 3; i++)
				for (j = 0; j < 3; j++) {
					double		sum = 0;
					for (k = 0; k < 3; k++) {
						sum += dTijdXk[i][k][l] * invTmatrix[k][j];
						printf("\n sum[%d][%d][%d] = %.12e", l, i, j, sum);
					}
				}
		}
	}
	if (debugPK)
		printf("\npData: {%.12e, %.12e, %.12e}\n", pData[0], pData[1], pData[2]);
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			for (k = 0, tmpPdotT3T12[i][j] = 0.; k < 3; k++)
				tmpPdotT3T12[i][j] += tmpPdotT3T11[i][j][k] * pData[k];

	for (i = 0; i < 3; i++)
		for (j = 0, tmpPdotT3T2[i] = 0.; j < 3; j++)
			tmpPdotT3T2[i] += tmpDValues[j + 3] * Tmatrix[i][j];

	for (i = 0; i < 3; i++)
		for (j = 0, tmpPdotT3[i] = 0.; j < 3; j++)
			tmpPdotT3[i] += tmpPdotT3T12[i][j] * tmpPdotT3T2[j];

	/* Add them to obtain pDot */
	for (i = 0; i < 3; i++)
		dvalues[i + 3] = tmpPdotT1[i] + tmpPdotT2[i] + tmpPdotT3[i];

	if (debugPK) {
		printf("\ntmpPdotT3 = ");
		for (i = 0; i < 3; i++)
			printf("%.12e ", tmpPdotT3[i]);
		printf("\n");
	}
	//dvalues[3] = -tmpDValues[0] - flux * values[3] / (omega * magL) + rrTerm2 * Lx;
	//dvalues[4] = -tmpDValues[1] - flux * values[4] / (omega * magL) + rrTerm2 * Ly;
	//dvalues[5] = -tmpDValues[2] - flux * values[5] / (omega * magL) + rrTerm2 * Lz;

	/* spin1 */
	if (debugPK) {
		printf("Raw spin1 derivatives = %.12e %.12e %.12e\n", tmpDValues[6], tmpDValues[7], tmpDValues[8]);
		printf("Raw spin2 derivatives = %.12e %.12e %.12e\n", tmpDValues[9], tmpDValues[10], tmpDValues[11]);
	}
	dvalues[6] = eta * (tmpDValues[7] * values[8] - tmpDValues[8] * values[7]);
	dvalues[7] = eta * (tmpDValues[8] * values[6] - tmpDValues[6] * values[8]);
	dvalues[8] = eta * (tmpDValues[6] * values[7] - tmpDValues[7] * values[6]);

	/* spin2 */
	dvalues[9] = eta * (tmpDValues[10] * values[11] - tmpDValues[11] * values[10]);
	dvalues[10] = eta * (tmpDValues[11] * values[9] - tmpDValues[9] * values[11]);
	dvalues[11] = eta * (tmpDValues[9] * values[10] - tmpDValues[10] * values[9]);

	/* phase and precessing bit */
	dLx = dvalues[1] * values[5] - dvalues[2] * values[4]
		+ values[1] * dvalues[5] - values[2] * dvalues[4];

	dLy = dvalues[2] * values[3] - dvalues[0] * values[5]
		+ values[2] * dvalues[3] - values[0] * dvalues[5];

	dLz = dvalues[0] * values[4] - dvalues[1] * values[3]
		+ values[0] * dvalues[4] - values[1] * dvalues[3];

	dMagL = (Lx * dLx + Ly * dLy + Lz * dLz) / magL;

	dLhatx = (dLx * magL - Lx * dMagL) / (magL * magL);
	dLhaty = (dLy * magL - Ly * dMagL) / (magL * magL);

	/*
	 * Finn Chernoff convention is used here. TODO: implement the
	 * geometric precessing convention
	 */
	if (Lhatx == 0.0 && Lhaty == 0.0) {
		alphadotcosi = 0.0;
	} else {
		alphadotcosi = Lhatz * (Lhatx * dLhaty - Lhaty * dLhatx) / (Lhatx * Lhatx + Lhaty * Lhaty);
	}

	dvalues[12] = omega - alphadotcosi;
	dvalues[13] = alphadotcosi;

	if (debugPK) {
    printf("\nIn XLALSpinHcapNumericalDerivative:\n");
		/* Print out all mass parameters */
		printf("m1 = %12.12lf, m2 = %12.12lf, eta = %12.12lf\n", 
          (double)mass1, (double)mass2, (double)eta);
		/* Print out all spin parameters */
		printf("spin1 = {%12.12lf,%12.12lf,%12.12lf}, spin2 = {%12.12lf,%12.12lf,%12.12lf}\n",
            (double)s1.data[0], (double)s1.data[1], (double)s1.data[2],
            (double)s2.data[0], (double)s2.data[1], (double)s2.data[2]);
		printf("sigmaStar = {%12.12lf,%12.12lf,%12.12lf}, sigmaKerr = {%12.12lf,%12.12lf,%12.12lf}\n",
		       (double)sStar.data[0], (double)sStar.data[1],
		       (double)sStar.data[2], (double)sKerr.data[0],
		       (double)sKerr.data[1], (double)sKerr.data[2]);
		printf("L = {%12.12lf,%12.12lf,%12.12lf}, |L| = %12.12lf\n", 
        (double)Lx, (double)Ly, (double)Lz, (double)magL);
		printf("dLdt = {%12.12lf,%12.12lf,%12.12lf}, d|L|dt = %12.12lf\n", 
        (double)dLx, (double)dLy, (double)dLz, (double)dMagL);
		printf("Polar coordinates = {%12.12lf, %12.12lf, %12.12lf, %12.12lf}\n",
		       (double)polData[0], (double)polData[1], (double)polData[2],
           (double)polData[3]);

		printf("Cartesian coordinates: {%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf}\n",
		       (double)values[0], (double)values[1], (double)values[2],
           (double)values[3], (double)values[4], (double)values[5], 
           (double)values[6], (double)values[7], (double)values[8], 
           (double)values[9], (double)values[10], (double)values[11],
		       (double)values[12], (double)values[13]);
		printf("Cartesian derivatives: {%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf}\n",
		       (double)dvalues[0], (double)dvalues[1], (double)dvalues[2],
           (double)dvalues[3], (double)dvalues[4], (double)dvalues[5],
           (double)dvalues[6], (double)dvalues[7], (double)dvalues[8],
           (double)dvalues[9], (double)dvalues[10], (double)dvalues[11],
           (double)dvalues[12], (double)dvalues[13]);

     printf("Hamiltonian = %12.12lf, Flux = %12.12lf, Omega = %12.12lf\n", 
              H/ (mass1 + mass2), eta*flux, omega);
		fflush(NULL);
	}
	
    for(i=0; i<14; i++){
        if(isnan(dvalues[i])){
            dvalues[i] = 1.;
          }
    }



  for(i=0; i<14; i++)
    if(dvalues[i] > 1e3){
      printf("\nIn XLALSpinHcapNumericalDerivative:\n");
      printf("Derivatives have blown up!\n");
      for(j=0; j<14; printf("dvalues[%d] = %3.12f\n", j, dvalues[j]), j++);
      printf("Flux = %3.12f\n\n", flux);
      break;
    }
      
  return XLAL_SUCCESS;
}

/**
 * Function to calculate numerical derivatives of the spin EOB Hamiltonian,
 * which correspond to time derivatives of the dynamical variables in conservative dynamcis.
 * All derivatives, including those on two terms of the orbital phase, are returned together.
 * The derivatives are combined with energy flux to give right hand side of the ODEs
 * of a generic spin EOB model, as decribed in Eqs. 21, 22, 26 and 27 of
 * Pan et al. PRD 81, 084041 (2010)
 * This function is not used by the spin-aligned SEOBNRv1 model.
 */
static int 
XLALSpinHcapNumericalDerivativeNoFlux(
				      double UNUSED t,	/**<< UNUSED */
				      const REAL8 values[],	/**<< Dynamical variables */
				      REAL8 dvalues[],	/**<< Time derivatives of variables (returned) */
				      void *funcParams	/**<< EOB parameters */
)
{
	int		debugPK = 0;
	static const REAL8 STEP_SIZE = 1.0e-4;

	static const INT4 lMax = 8;

	HcapDerivParams	params;

	/* Since we take numerical derivatives wrt dynamical variables */
	/* but we want them wrt time, we use this temporary vector in  */
	/* the conversion */
	REAL8		tmpDValues[14];

	REAL8		H;
	//Hamiltonian
		REAL8 flux;

	gsl_function	F;
	INT4		gslStatus;
	UINT4		SpinAlignedEOBversion;

	UINT4		i       , j, k, l;

	REAL8Vector	rVec, pVec;
	REAL8		rData    [3], pData[3];

	/* We need r, phi, pr, pPhi to calculate the flux */
	REAL8		r;
	REAL8Vector	polarDynamics, cartDynamics;
	REAL8		polData  [4];

	REAL8		mass1   , mass2, eta;
	REAL8 UNUSED	rrTerm2, pDotS1, pDotS2;
	REAL8Vector	s1 , s2, s1norm, s2norm, sKerr, sStar;
	REAL8		s1Data   [3], s2Data[3], s1DataNorm[3], s2DataNorm[3];
	REAL8		sKerrData[3], sStarData[3];
	REAL8 /* magS1, magS2, */ chiS, chiA, a, tplspin;
	REAL8 UNUSED	s1dotL, s2dotL;
	REAL8 UNUSED	rcrossrDot[3], rcrossrDotMag, s1dotLN, s2dotLN;


	/* Orbital angular momentum */
	REAL8		Lx      , Ly, Lz, magL;
	REAL8		Lhatx   , Lhaty, Lhatz;
	REAL8		dLx     , dLy, dLz;
	REAL8		dLhatx  , dLhaty, dMagL;

	REAL8		alphadotcosi;

	REAL8		rCrossV_x, rCrossV_y, rCrossV_z, omega;

	/* The error in a derivative as measured by GSL */
	REAL8		absErr;

	REAL8		tmpP     [3], rMag, rMag2, prT;
	REAL8		u       , u2, u3, u4, u5, w2, a2;
	REAL8		D       , m1PlusetaKK, bulk, logTerms, deltaU, deltaT, deltaR;
	REAL8 UNUSED	eobD_r, deltaU_u, deltaU_r, deltaT_r;
	REAL8		dcsi    , csi;

	REAL8		tmpValues[12];
	REAL8		Tmatrix  [3][3], invTmatrix[3][3], dTijdXk[3][3][3];
	REAL8		tmpPdotT1[3], tmpPdotT2[3], tmpPdotT3[3];
	//3 terms of Eq.A5

	/* NQC coefficients container */
		EOBNonQCCoeffs * nqcCoeffs = NULL;

	/* Set up pointers for GSL */
	params.values = values;
	params.params = (SpinEOBParams *) funcParams;
	nqcCoeffs = params.params->nqcCoeffs;

	F.function = &GSLSpinHamiltonianWrapper;
	F.params = &params;

	mass1 = params.params->eobParams->m1;
	mass2 = params.params->eobParams->m2;
	eta = params.params->eobParams->eta;
	SpinAlignedEOBversion = params.params->seobCoeffs->SpinAlignedEOBversion;
	SpinEOBHCoeffs *coeffs = (SpinEOBHCoeffs *) params.params->seobCoeffs;

	/*
	 * For precessing binaries, the effective spin of the Kerr background
	 * evolves with time. The coefficients used to compute the
	 * Hamiltonian depend on the Kerr spin, and hence need to be updated
	 * for the current spin values
	 */

	/*
	 * Set the position/momenta vectors to point to the appropriate
	 * things
	 */
	rVec.length = pVec.length = 3;
	rVec.data = rData;
	pVec.data = pData;
	memcpy(rData, values, sizeof(rData));
	memcpy(pData, values + 3, sizeof(pData));

	/*
	 * We need to re-calculate the parameters at each step as precessing
	 * spins will not be constant
	 */
	/* TODO: Modify so that only spin terms get re-calculated */

	/*
	 * We cannot point to the values vector directly as it leads to a
	 * warning
	 */
	s1.length = s2.length = s1norm.length = s2norm.length = 3;
	s1.data = s1Data;
	s2.data = s2Data;
	s1norm.data = s1DataNorm;
	s2norm.data = s2DataNorm;

	memcpy(s1Data, values + 6, 3 * sizeof(REAL8));
	memcpy(s2Data, values + 9, 3 * sizeof(REAL8));
	memcpy(s1DataNorm, values + 6, 3 * sizeof(REAL8));
	memcpy(s2DataNorm, values + 9, 3 * sizeof(REAL8));

	for (i = 0; i < 3; i++) {
		s1Data[i] *= (mass1 + mass2) * (mass1 + mass2);
		s2Data[i] *= (mass1 + mass2) * (mass1 + mass2);
	}

	sKerr.length = 3;
	sKerr.data = sKerrData;
	XLALSimIMRSpinEOBCalculateSigmaKerr(&sKerr, mass1, mass2, &s1, &s2);

	sStar.length = 3;
	sStar.data = sStarData;
	XLALSimIMRSpinEOBCalculateSigmaStar(&sStar, mass1, mass2, &s1, &s2);

	a = sqrt(sKerr.data[0] * sKerr.data[0] + sKerr.data[1] * sKerr.data[1]
		 + sKerr.data[2] * sKerr.data[2]);

	///*set the tortoise flag to 2 * /
		//INT4 oldTortoise = params.params->tortoise;
	//params.params->tortoise = 2;

	/* Convert momenta to p */
	rMag = sqrt(rData[0] * rData[0] + rData[1] * rData[1] + rData[2] * rData[2]);
	prT = pData[0] * (rData[0] / rMag) + pData[1] * (rData[1] / rMag)
		+ pData[2] * (rData[2] / rMag);

	rMag2 = rMag * rMag;
	u = 1. / rMag;
	u2 = u * u;
	u3 = u2 * u;
	u4 = u2 * u2;
	u5 = u4 * u;
	a2 = a * a;
	w2 = rMag2 + a2;
	/* Eq. 5.83 of BB1, inverse */
	D = 1. + log(1. + 6. * eta * u2 + 2. * (26. - 3. * eta) * eta * u3);
	eobD_r = (u2 / (D * D)) * (12. * eta * u + 6. * (26. - 3. * eta) * eta * u2) / (1.
			+ 6. * eta * u2 + 2. * (26. - 3. * eta) * eta * u3);
	m1PlusetaKK = -1. + eta * coeffs->KK;
	/* Eq. 5.75 of BB1 */
	bulk = 1. / (m1PlusetaKK * m1PlusetaKK) + (2. * u) / m1PlusetaKK + a2 * u2;
	/* Eq. 5.73 of BB1 */
	logTerms = 1. + eta * coeffs->k0 + eta * log(1. + coeffs->k1 * u
		       + coeffs->k2 * u2 + coeffs->k3 * u3 + coeffs->k4 * u4
			     + coeffs->k5 * u5 + coeffs->k5l * u5 * log(u));
	/* Eq. 5.73 of BB1 */
	deltaU = bulk * logTerms;
	/* Eq. 5.71 of BB1 */
	deltaT = rMag2 * deltaU;
	/* ddeltaU/du */
	deltaU_u = 2. * (1. / m1PlusetaKK + a2 * u) * logTerms +
		bulk * (eta * (coeffs->k1 + u * (2. * coeffs->k2 + u * (3. * coeffs->k3
									+ u * (4. * coeffs->k4 + 5. * (coeffs->k5 + coeffs->k5l * log(u)) * u)))))
		/ (1. + coeffs->k1 * u + coeffs->k2 * u2 + coeffs->k3 * u3
	      + coeffs->k4 * u4 + (coeffs->k5 + coeffs->k5l * log(u)) * u5);
	deltaU_r = -u2 * deltaU_u;
	/* Eq. 5.38 of BB1 */
	deltaR = deltaT * D;
	if (params.params->tortoise)
		csi = sqrt(deltaT * deltaR) / w2;	/* Eq. 28 of Pan et al.
							 * PRD 81, 084041 (2010) */
	else
		csi = 1.0;

	for (i = 0; i < 3; i++) {
		tmpP[i] = pData[i] - (rData[i] / rMag) * prT * (csi - 1.) / csi;
	}

	if (debugPK) {
		printf("csi = %.12e\n", csi);
		for (i = 0; i < 3; i++)
			printf("p,p*: %.12e\t%.12e\n", pData[i], tmpP[i]);
	}
	/*
	 * Calculate the T-matrix, required to convert P from tortoise to
	 * non-tortoise coordinates, and/or vice-versa. This is given
	 * explicitly in Eq. A3 of 0912.3466
	 */
	for (i = 0; i < 3; i++)
		for (j = 0; j <= i; j++) {
			Tmatrix[i][j] = Tmatrix[j][i] = (rData[i] * rData[j] / rMag2)
				* (csi - 1.);

			invTmatrix[i][j] = invTmatrix[j][i] =
				-(csi - 1) / csi * (rData[i] * rData[j] / rMag2);

			if (i == j) {
				Tmatrix[i][j]++;
				invTmatrix[i][j]++;
			}
		}

	dcsi = csi * (2. / rMag + deltaU_r / deltaU) + csi * csi * csi
		/ (2. * rMag2 * rMag2 * deltaU * deltaU) * (rMag * (-4. * w2) / D - eobD_r * (w2 * w2));

	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			for (k = 0; k < 3; k++) {
				dTijdXk[i][j][k] =
					(rData[i] * XLALKronecker(j, k) + XLALKronecker(i, k) * rData[j])
					* (csi - 1.) / rMag2
					+ rData[i] * rData[j] * rData[k] / rMag2 / rMag * (-2. / rMag * (csi - 1.) + dcsi);
			}

	//Print out the T - matrix for comparison
		if (debugPK) {
			printf("\nT-Matrix:\n");
			for (i = 0; i < 3; i++)
				printf("%le\t%le\t%le\n", Tmatrix[i][0], Tmatrix[i][1], Tmatrix[i][2]);

			for (i = 0; i < 3; i++) {
				printf("dT[%d][j]dX[k]:\n", i);
				for (j = 0; j < 3; j++)
					printf("%.12e\t%.12e\t%.12e\n", dTijdXk[i][j][0],
					dTijdXk[i][j][1], dTijdXk[i][j][2]);
				printf("\n");
			}
		}
	/* Now calculate derivatives w.r.t. each parameter */
	for (i = 0; i < 12; i++) {
		params.varyParam = i;
		if (i >= 6 && i < 9) {
			params.params->seobCoeffs->updateHCoeffs = 1;
			XLAL_CALLGSL(gslStatus = gsl_deriv_central(&F, values[i],
			STEP_SIZE * mass1 * mass1, &tmpDValues[i], &absErr));
		} else if (i >= 9) {
			params.params->seobCoeffs->updateHCoeffs = 1;
			XLAL_CALLGSL(gslStatus = gsl_deriv_central(&F, values[i],
			STEP_SIZE * mass2 * mass2, &tmpDValues[i], &absErr));
		} else if (i < 3) {
			params.params->tortoise = 2;
			memcpy(tmpValues, params.values, sizeof(tmpValues));
			tmpValues[3] = tmpP[0];
			tmpValues[4] = tmpP[1];
			tmpValues[5] = tmpP[2];
			params.values = tmpValues;

			XLAL_CALLGSL(gslStatus = gsl_deriv_central(&F, values[i],
				       STEP_SIZE, &tmpDValues[i], &absErr));

			params.values = values;
			params.params->tortoise = 1;
		} else {
			params.params->seobCoeffs->updateHCoeffs = 1;
			XLAL_CALLGSL(gslStatus = gsl_deriv_central(&F, values[i],
				       STEP_SIZE, &tmpDValues[i], &absErr));
		}
		if (gslStatus != GSL_SUCCESS) {
			XLALPrintError("XLAL Error %s - Failure in GSL function\n", __func__);
			XLAL_ERROR(XLAL_EFUNC);
		}
	}

	/* Now make the conversion */
	/* rVectorDot */
	for (i = 0; i < 3; i++)
		for (j = 0, dvalues[i] = 0.; j < 3; j++)
			dvalues[i] += tmpDValues[j + 3] * Tmatrix[i][j];

	//dvalues[0] = tmpDValues[3] * Tmatrix[0][0] + tmpDValues[4] * Tmatrix[0][1]
		// +tmpDValues[5] * Tmatrix[0][2];
	//dvalues[1] = tmpDValues[3] * Tmatrix[1][0] + tmpDValues[4] * Tmatrix[1][1]
		// +tmpDValues[5] * Tmatrix[1][2];
	//dvalues[2] = tmpDValues[3] * Tmatrix[2][0] + tmpDValues[4] * Tmatrix[2][1]
		// +tmpDValues[5] * Tmatrix[2][2];

	/* Calculate the orbital angular momentum */
	Lx = values[1] * values[5] - values[2] * values[4];
	Ly = values[2] * values[3] - values[0] * values[5];
	Lz = values[0] * values[4] - values[1] * values[3];

	magL = sqrt(Lx * Lx + Ly * Ly + Lz * Lz);

	Lhatx = Lx / magL;
	Lhaty = Ly / magL;
	Lhatz = Lz / magL;

	/* Calculate the polar data */
	polarDynamics.length = 4;
	polarDynamics.data = polData;

	r = polData[0] = sqrt(values[0] * values[0] + values[1] * values[1]
			      + values[2] * values[2]);
	polData[1] = 0;
	polData[2] = (values[0] * values[3] + values[1] * values[4]
		      + values[2] * values[5]) / polData[0];
	polData[3] = magL;


	/* Compute \vec{S_i} \dot \vec{L}	 */
	s1dotL = (s1Data[0] * Lhatx + s1Data[1] * Lhaty + s1Data[2] * Lhatz)
		/ (mass1 * mass1);
	s2dotL = (s2Data[0] * Lhatx + s2Data[1] * Lhaty + s2Data[2] * Lhatz)
		/ (mass2 * mass2);

	/* Now calculate rCrossRdot and omega */
	rCrossV_x = values[1] * dvalues[2] - values[2] * dvalues[1];
	rCrossV_y = values[2] * dvalues[0] - values[0] * dvalues[2];
	rCrossV_z = values[0] * dvalues[1] - values[1] * dvalues[0];

	omega = sqrt(rCrossV_x * rCrossV_x + rCrossV_y * rCrossV_y + rCrossV_z * rCrossV_z) / (r * r);

	/*
	 * Compute \vec{L_N} = \vec{r} \times \.{\vec{r}}, \vec{S_i} \dot
	 * \vec{L_N} and chiS and chiA
	 */
	/*
	 * THIS IS NOT QUITE RIGHT. ITS rcrossr*Dot INSTEAD, as it is
	 * \partial H/\partial p_*
	 */
	/*
	 * rcrossrDot[0] = values[1]*tmpDValues[5] - values[2]*tmpDValues[4];
	 * rcrossrDot[1] = values[2]*tmpDValues[3] - values[0]*tmpDValues[5];
	 * rcrossrDot[2] = values[0]*tmpDValues[4] - values[1]*tmpDValues[3];
	 * rcrossrDotMag = sqrt( rcrossrDot[0]*rcrossrDot[0] +
	 * rcrossrDot[1]*rcrossrDot[1]	+ rcrossrDot[2]*rcrossrDot[2] );
	 * 
	 * rcrossrDot[0] /= rcrossrDotMag; rcrossrDot[1] /= rcrossrDotMag;
	 * rcrossrDot[2] /= rcrossrDotMag;
	 */

	s1dotLN = (s1Data[0] * rCrossV_x + s1Data[1] * rCrossV_y + s1Data[2] * rCrossV_z) /
		(r * r * omega * mass1 * mass1);
	s2dotLN = (s2Data[0] * rCrossV_x + s2Data[1] * rCrossV_y + s2Data[2] * rCrossV_z) /
		(r * r * omega * mass2 * mass2);

	chiS = 0.5 * (s1dotLN + s2dotLN);
	chiA = 0.5 * (s1dotLN - s2dotLN);

	if (debugPK) {
		printf("chiS = %.12e, chiA = %.12e\n", chiS, chiA);
		fflush(NULL);
	}
	REAL8		sscaling1 = 1;
	//(mass1 + mass2) * (mass1 + mass2) / (mass1 * mass1);
	REAL8		sscaling2 = 1;
	//(mass1 + mass2) * (mass1 + mass2) / (mass2 * mass2);

	if (debugPK) {
		printf("Computing derivatives for values\n%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n\n",
		    (double)values[0], (double)values[1], (double)values[2],
		    (double)values[3], (double)values[4], (double)values[5],
		(double)sscaling1 * values[6], (double)sscaling1 * values[7],
		(double)sscaling1 * values[8], (double)sscaling2 * values[9],
		       (double)sscaling2 * values[10], (double)sscaling2 * values[11]);
		printf("tmpDvalues\n%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t\n",
		       (double)tmpDValues[0], (double)tmpDValues[1], (double)tmpDValues[2],
		       (double)tmpDValues[3], (double)tmpDValues[4], (double)tmpDValues[5],
		       (double)tmpDValues[6], (double)tmpDValues[7], (double)tmpDValues[8],
		       (double)tmpDValues[9], (double)tmpDValues[10], (double)tmpDValues[11]);
	}
	/*
	 * Compute the test-particle limit spin of the deformed-Kerr
	 * background
	 */
	/* TODO: Check this is actually the way it works in LAL */
	switch (SpinAlignedEOBversion) {
	case 1:
		tplspin = 0.0;
		break;
	case 2:
		tplspin = (1. - 2. * eta) * chiS + (mass1 - mass2) / (mass1 + mass2) * chiA;
		break;
	default:
		XLALPrintError("XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
		XLAL_ERROR(XLAL_EINVAL);
		break;
	}

	for (i = 0; i < 3; i++) {
		params.params->s1Vec->data[i] = s1norm.data[i];
		params.params->s2Vec->data[i] = s2norm.data[i];
		params.params->sigmaStar->data[i] = sStar.data[i];
		params.params->sigmaKerr->data[i] = sKerr.data[i];
	}

	//params.params->s1Vec = &s1norm;
	//params.params->s2Vec = &s2norm;
	//params.params->sigmaStar = &sStar;
	//params.params->sigmaKerr = &sKerr;
	params.params->a = a;

	XLALSimIMREOBCalcSpinFacWaveformCoefficients(
	      params.params->eobParams->hCoeffs, mass1, mass2, eta, tplspin,
					 chiS, chiA, SpinAlignedEOBversion);
	XLALSimIMRCalculateSpinEOBHCoeffs(params.params->seobCoeffs, eta, a,
					  SpinAlignedEOBversion);

	H = XLALSimIMRSpinEOBHamiltonian(eta, &rVec, &pVec, &s1norm, &s2norm,
	&sKerr, &sStar, params.params->tortoise, params.params->seobCoeffs);
	H = H * (mass1 + mass2);

	/* Now we have the ingredients to compute the flux */
	memcpy(tmpValues, values, 12 * sizeof(REAL8));
	cartDynamics.data = tmpValues;
	if (debugPK) {
		printf("params.params->a = %.12e, %.12e\n", a, params.params->a);
		fflush(NULL);
	}
	flux = 0 * XLALInspiralPrecSpinFactorizedFlux(&polarDynamics, &cartDynamics,
						      nqcCoeffs, omega, params.params, H / (mass1 + mass2), lMax, SpinAlignedEOBversion);

	/*
	 * Looking at the non-spinning model, I think we need to divide the
	 * flux by eta
	 */
	flux = flux / eta;

	pDotS1 = pData[0] * s1Data[0] + pData[1] * s1Data[1] + pData[2] * s1Data[2];
	pDotS2 = pData[0] * s2Data[0] + pData[1] * s2Data[1] + pData[2] * s2Data[2];
	rrTerm2 = 8. / 15. * eta * eta * pow(omega, 8. / 3.) / (magL * magL * r) * ((61. + 48. * mass2 / mass1) * pDotS1 + (61. + 48. * mass1 / mass2) * pDotS2);

	if (debugPK) {
		printf("omega = %.12e \n flux = %.12e \n Lmag = %.12e\n", omega, flux, magL);
		printf("rrForce = %.12e %.12e %.12e\n", -flux * values[3] / (omega * magL), -flux * values[4] / (omega * magL), -flux * values[5] / (omega * magL));
	}
	/* Now pDot */
	/* Compute the first and second terms in eq. A5 of 0912.3466 */
	for (i = 0; i < 3; i++) {
		for (j = 0, tmpPdotT1[i] = 0.; j < 3; j++)
			tmpPdotT1[i] += -tmpDValues[j] * Tmatrix[i][j];

		tmpPdotT2[i] = -flux * values[i + 3] / (omega * magL);
	}

	/* Compute the third term in eq. A5 */
	REAL8		tmpPdotT3T11[3][3][3], tmpPdotT3T12[3][3], tmpPdotT3T2[3];

	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			for (l = 0; l < 3; l++)
				for (k = 0, tmpPdotT3T11[i][j][l] = 0.; k < 3; k++)
					tmpPdotT3T11[i][j][l] += dTijdXk[i][k][j] * invTmatrix[k][l];

	if (debugPK) {
		for (i = 0; i < 3; i++)
			for (j = 0; j < 1; j++)
				for (l = 0; l < 3; l++) {
					double		sum = 0;
					for (k = 0; k < 3; k++)
						sum += dTijdXk[i][k][j] * invTmatrix[k][l];
					printf("\n sum[%d][%d][%d] = %.12e", i, j, l, sum);

				}
		printf("\n\n Printing dTdX * Tinverse:\n");
		for (l = 0; l < 3; l++) {
			for (i = 0; i < 3; i++)
				for (j = 0; j < 3; j++) {
					double		sum = 0;
					for (k = 0; k < 3; k++) {
						sum += dTijdXk[i][k][l] * invTmatrix[k][j];
						printf("\n sum[%d][%d][%d] = %.12e", l, i, j, sum);
					}
				}
		}
	}
	if (debugPK)
		printf("\npData: {%.12e, %.12e, %.12e}\n", pData[0], pData[1], pData[2]);
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			for (k = 0, tmpPdotT3T12[i][j] = 0.; k < 3; k++)
				tmpPdotT3T12[i][j] += tmpPdotT3T11[i][j][k] * pData[k];

	for (i = 0; i < 3; i++)
		for (j = 0, tmpPdotT3T2[i] = 0.; j < 3; j++)
			tmpPdotT3T2[i] += tmpDValues[j + 3] * Tmatrix[i][j];

	for (i = 0; i < 3; i++)
		for (j = 0, tmpPdotT3[i] = 0.; j < 3; j++)
			tmpPdotT3[i] += tmpPdotT3T12[i][j] * tmpPdotT3T2[j];

	/* Add them to obtain pDot */
	for (i = 0; i < 3; i++)
		dvalues[i + 3] = tmpPdotT1[i] + tmpPdotT2[i] + tmpPdotT3[i];

	if (debugPK) {
		printf("\ntmpPdotT3 = ");
		for (i = 0; i < 3; i++)
			printf("%.12e ", tmpPdotT3[i]);
		printf("\n");
	}
	//dvalues[3] = -tmpDValues[0] - flux * values[3] / (omega * magL) + rrTerm2 * Lx;
	//dvalues[4] = -tmpDValues[1] - flux * values[4] / (omega * magL) + rrTerm2 * Ly;
	//dvalues[5] = -tmpDValues[2] - flux * values[5] / (omega * magL) + rrTerm2 * Lz;

	/* spin1 */
	if (debugPK) {
		printf("Raw spin1 derivatives = %.12e %.12e %.12e\n", tmpDValues[6], tmpDValues[7], tmpDValues[8]);
		printf("Raw spin2 derivatives = %.12e %.12e %.12e\n", tmpDValues[9], tmpDValues[10], tmpDValues[11]);
	}
	dvalues[6] = eta * (tmpDValues[7] * values[8] - tmpDValues[8] * values[7]);
	dvalues[7] = eta * (tmpDValues[8] * values[6] - tmpDValues[6] * values[8]);
	dvalues[8] = eta * (tmpDValues[6] * values[7] - tmpDValues[7] * values[6]);

	/* spin2 */
	dvalues[9] = eta * (tmpDValues[10] * values[11] - tmpDValues[11] * values[10]);
	dvalues[10] = eta * (tmpDValues[11] * values[9] - tmpDValues[9] * values[11]);
	dvalues[11] = eta * (tmpDValues[9] * values[10] - tmpDValues[10] * values[9]);

	/* phase and precessing bit */
	dLx = dvalues[1] * values[5] - dvalues[2] * values[4]
		+ values[1] * dvalues[5] - values[2] * dvalues[4];

	dLy = dvalues[2] * values[3] - dvalues[0] * values[5]
		+ values[2] * dvalues[3] - values[0] * dvalues[5];

	dLz = dvalues[0] * values[4] - dvalues[1] * values[3]
		+ values[0] * dvalues[4] - values[1] * dvalues[3];

	dMagL = (Lx * dLx + Ly * dLy + Lz * dLz) / magL;

	dLhatx = (dLx * magL - Lx * dMagL) / (magL * magL);
	dLhaty = (dLy * magL - Ly * dMagL) / (magL * magL);

	/*
	 * Finn Chernoff convention is used here. TODO: implement the
	 * geometric precessing convention
	 */
	if (Lhatx == 0.0 && Lhaty == 0.0) {
		alphadotcosi = 0.0;
	} else {
		alphadotcosi = Lhatz * (Lhatx * dLhaty - Lhaty * dLhatx) / (Lhatx * Lhatx + Lhaty * Lhaty);
	}

	dvalues[12] = omega - alphadotcosi;
	dvalues[13] = alphadotcosi;

	if (debugPK) {
		printf(" r = %e %e %e (mag = %e)\n", values[0], values[1], values[2], sqrt(values[0] * values[0] + values[1] * values[1] + values[2] * values[2]));
		printf(" p = %e %e %e (mag = %e)\n", values[3], values[4], values[5], sqrt(values[3] * values[3] + values[4] * values[4] + values[5] * values[5]));
		printf("Derivatives:\n");
		for (i = 0; i < 12; i++) {
			printf("%.12e\n", dvalues[i]);
		}
		printf("\n");
	}
	if (debugPK) {
#if 0
		/* Print out all mass parameters */
		printf("\nIn XLALSpinHcapNumericalDerivative:\n");
		printf("m1 = %12.12lf, m2 = %12.12lf, eta = %12.12lf\n", (double)mass1,
		       (double)mass2, (double)eta);
		/* Print out all spin parameters */
		printf("spin1 = {%12.12lf,%12.12lf,%12.12lf}, spin2 = {%12.12lf,%12.12lf,%12.12lf}\n",
		 (double)s1.data[0], (double)s1.data[1], (double)s1.data[2],
		(double)s2.data[0], (double)s2.data[1], (double)s2.data[2]);
		printf("sigmaStar = {%12.12lf,%12.12lf,%12.12lf}, sigmaKerr = {%12.12lf,%12.12lf,%12.12lf}\n",
		       (double)sStar.data[0], (double)sStar.data[1],
		       (double)sStar.data[2], (double)sKerr.data[0],
		       (double)sKerr.data[1], (double)sKerr.data[2]);
		printf("L = {%12.12lf,%12.12lf,%12.12lf}, |L| = %12.12lf\n", (double)Lx, (double)Ly,
		       (double)Lz, (double)magL);
		printf("dLdt = {%12.12lf,%12.12lf,%12.12lf}, d|L|dt = %12.12lf\n", (double)dLx,
		       (double)dLy, (double)dLz, (double)dMagL);
		printf("Polar coordinates = {%12.12lf, %12.12lf, %12.12lf, %12.12lf}\n",
		       (double)polData[0], (double)polData[1], (double)polData[2], (double)polData[3]);

		printf("Cartesian coordinates: {%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf}\n",
		       (double)values[0], (double)values[1], (double)values[2], (double)values[3],
		       (double)values[4], (double)values[5], (double)values[6], (double)values[7],
		       (double)values[8], (double)values[9], (double)values[10], (double)values[11],
		       (double)values[12], (double)values[13]);
		printf("Cartesian derivatives: {%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf}\n",
		       (double)dvalues[0], (double)dvalues[1], (double)dvalues[2], (double)dvalues[3],
		       (double)dvalues[4], (double)dvalues[5], (double)dvalues[6], (double)dvalues[7],
		       (double)dvalues[8], (double)dvalues[9], (double)dvalues[10], (double)dvalues[11],
		       (double)dvalues[12], (double)dvalues[13]);

		printf("Hamiltonian = %12.12lf, Flux = %12.12lf, Omega = %12.12lf\n", H, flux, omega);
#endif
#if 1
		if (debugPK)
			printf("%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n\n",
			       (double)values[0], (double)values[1], (double)values[2],
			       (double)values[3], (double)values[4], (double)values[5],
			       (double)sscaling1 * values[6], (double)sscaling1 * values[7],
			       (double)sscaling1 * values[8], (double)sscaling2 * values[9],
			       (double)sscaling2 * values[10], (double)sscaling2 * values[11],
			       (double)dvalues[0], (double)dvalues[1], (double)dvalues[2],
			       (double)dvalues[3], (double)dvalues[4], (double)dvalues[5],
			       (double)sscaling1 * dvalues[6], (double)sscaling1 * dvalues[7],
			       (double)sscaling1 * dvalues[8], (double)sscaling2 * dvalues[9],
			       (double)sscaling2 * dvalues[10], (double)sscaling2 * dvalues[11],
			       (double)polData[0], (double)polData[1], (double)polData[2],
			       (double)polData[3], H / (mass1 + mass2), flux * eta, omega);
#endif
		fflush(NULL);
	}
	return XLAL_SUCCESS;
}

/**
 * Calculate the derivative of the Hamiltonian w.r.t. a specific parameter
 * Used by generic spin EOB model, including initial conditions solver.
 */
static REAL8 
XLALSpinHcapNumDerivWRTParam(
			     const INT4 paramIdx,	/**<< Index of the parameters */
			     const REAL8 values[],	/**<< Dynamical variables */
			     SpinEOBParams * funcParams	/**<< EOB Parameters */
)
{
	static const REAL8 STEP_SIZE = 1.0e-3;

	HcapDerivParams	params;

	REAL8		result;

	gsl_function	F;
	INT4		gslStatus;

	REAL8		mass1   , mass2;

	/* The error in a derivative as measured by GSL */
	REAL8		absErr;

	/* Set up pointers for GSL */
	params.values = values;
	params.params = funcParams;

	F.function = &GSLSpinHamiltonianWrapperV2;
	F.params = &params;
	params.varyParam = paramIdx;

	mass1 = params.params->eobParams->m1;
	mass2 = params.params->eobParams->m2;

	/* Now calculate derivatives w.r.t. the required parameter */
	if (paramIdx >= 6 && paramIdx < 9) {
		XLAL_CALLGSL(gslStatus = gsl_deriv_central(&F, values[paramIdx],
			      STEP_SIZE * mass1 * mass1, &result, &absErr));
	} else if (paramIdx >= 9) {
		XLAL_CALLGSL(gslStatus = gsl_deriv_central(&F, values[paramIdx],
			      STEP_SIZE * mass2 * mass2, &result, &absErr));
	} else {
		XLAL_CALLGSL(gslStatus = gsl_deriv_central(&F, values[paramIdx],
					      STEP_SIZE, &result, &absErr));
	}
	if (gslStatus != GSL_SUCCESS) {
		XLALPrintError("XLAL Error %s - Failure in GSL function\n", __func__);
		XLAL_ERROR_REAL8(XLAL_EFUNC);
	}
	//printf("Abserr = %e\n", absErr);

	return result;
}


/**
 * Wrapper for GSL to call the Hamiltonian function
 */
/* static */
REAL8 
GSLSpinHamiltonianWrapper(double x, void *params)
{
	int		debugPK = 0;
	HcapDerivParams *dParams = (HcapDerivParams *) params;

	EOBParams      *eobParams = (EOBParams *) dParams->params->eobParams;
	SpinEOBHCoeffs UNUSED *coeffs = (SpinEOBHCoeffs *) dParams->params->seobCoeffs;

	REAL8		tmpVec   [12];
	REAL8		s1normData[3], s2normData[3], sKerrData[3], sStarData[3];

	/*
	 * These are the vectors which will be used in the call to the
	 * Hamiltonian
	 */
	REAL8Vector	r  , p, spin1, spin2, spin1norm, spin2norm;
	REAL8Vector	sigmaKerr, sigmaStar;

	INT4		i;
	REAL8		a;
	REAL8		m1 = eobParams->m1;
	REAL8		m2 = eobParams->m2;
	REAL8 UNUSED	mT2 = (m1 + m2) * (m1 + m2);
	REAL8 UNUSED	eta = m1 * m2 / mT2;

	INT4		oldTortoise = dParams->params->tortoise;
	/* Use a temporary vector to avoid corrupting the main function */
	memcpy(tmpVec, dParams->values, sizeof(tmpVec));

	/* Set the relevant entry in the vector to the correct value */
	tmpVec[dParams->varyParam] = x;

	/* Set the LAL-style vectors to point to the appropriate things */
	r.length = p.length = spin1.length = spin2.length = spin1norm.length = spin2norm.length = 3;
	sigmaKerr.length = sigmaStar.length = 3;
	r.data = tmpVec;
	p.data = tmpVec + 3;

	spin1.data = tmpVec + 6;
	spin2.data = tmpVec + 9;
	spin1norm.data = s1normData;
	spin2norm.data = s2normData;
	sigmaKerr.data = sKerrData;
	sigmaStar.data = sStarData;

	memcpy(s1normData, tmpVec + 6, 3 * sizeof(REAL8));
	memcpy(s2normData, tmpVec + 9, 3 * sizeof(REAL8));

	/*
	 * To compute the SigmaKerr and SigmaStar, we need the non-normalized
	 * spin values, i.e. S_i. The spins being evolved are S_i/M^2.
	 */
	for (i = 0; i < 3; i++) {
		spin1.data[i] *= mT2;
		spin2.data[i] *= mT2;

		//s1normData[i] /= mT2;
		//s2normData[i] /= mT2;
	}

	/* Calculate various spin parameters */
	XLALSimIMRSpinEOBCalculateSigmaKerr(&sigmaKerr, eobParams->m1,
					    eobParams->m2, &spin1, &spin2);
	XLALSimIMRSpinEOBCalculateSigmaStar(&sigmaStar, eobParams->m1,
					    eobParams->m2, &spin1, &spin2);
	a = sqrt(sigmaKerr.data[0] * sigmaKerr.data[0]
		 + sigmaKerr.data[1] * sigmaKerr.data[1]
		 + sigmaKerr.data[2] * sigmaKerr.data[2]);
	//printf("a = %e\n", a);
	//printf("aStar = %e\n", sqrt(sigmaStar.data[0] * sigmaStar.data[0] + sigmaStar.data[1] * sigmaStar.data[1] + sigmaStar.data[2] * sigmaStar.data[2]));
	if (isnan(a) ) {
          XLALPrintError( "XLAL Error - %s: a = nan  \n", __func__);
          XLAL_ERROR( XLAL_EINVAL );
	}
	if (isnan(a) && debugPK) {
		printf("a is nan!!\n");
	}
	//XLALSimIMRCalculateSpinEOBHCoeffs(dParams->params->seobCoeffs, eobParams->eta, a);
	/*
	 * If computing the derivative w.r.t. the position vector, we need to
	 * pass p and not p* to the Hamiltonian. This is so because we want
	 * to hold p constant as we compute dH/dx. The way to do this is set
	 * the tortoise flag = 2, and convert the momenta being evolved, i.e.
	 * p*, to p and pass that as input
	 */
	//TODO
#if 0
		if (dParams->varyParam < 3) {
		/* set the tortoise flag to 2 */
		dParams->params->tortoise = 2;

		/* Convert momenta to p */
		REAL8		tmpP     [3];
		REAL8		rMag = sqrt(r.data[0] * r.data[0] + r.data[1] * r.data[1]
				   + r.data[2] * r.data[2]);
		REAL8		prT = p.data[0] * (r.data[0] / rMag) + p.data[1] * (r.data[1] / rMag)
		+ p.data[2] * (r.data[2] / rMag);

		REAL8		u       , u2, u3, u4, u5, w2, a2;
		REAL8		csi;

		u = 1. / rMag;
		u2 = u * u;
		u3 = u2 * u;
		u4 = u2 * u2;
		u5 = u4 * u;
		a2 = a * a;
		w2 = rMag * rMag + a2;
		/* Eq. 5.83 of BB1, inverse */
		REAL8		D = 1. + log(1. + 6. * eta * u2 + 2. * (26. - 3. * eta) * eta * u3);
		REAL8		m1PlusetaKK = -1. + eta * coeffs->KK;
		/* Eq. 5.75 of BB1 */
		REAL8		bulk = 1. / (m1PlusetaKK * m1PlusetaKK) + (2. * u) / m1PlusetaKK + a2 * u2;
		/* Eq. 5.73 of BB1 */
		REAL8		logTerms = 1. + eta * coeffs->k0 + eta * log(1. + coeffs->k1 * u
		       + coeffs->k2 * u2 + coeffs->k3 * u3 + coeffs->k4 * u4
			     + coeffs->k5 * u5 + coeffs->k5l * u5 * log(u));
		/* Eq. 5.73 of BB1 */
		REAL8		deltaU = bulk * logTerms;
		/* Eq. 5.71 of BB1 */
		REAL8		deltaT = rMag * rMag * deltaU;
		/* Eq. 5.38 of BB1 */
		REAL8		deltaR = deltaT * D;
		if (oldTortoise)
			csi = sqrt(deltaT * deltaR) / w2;	/* Eq. 28 of Pan et al.
								 * PRD 81, 084041 (2010) */
		else
			csi = 1.0;

		for (i = 0; i < 3; i++) {
			tmpP[i] = p.data[i] - (r.data[i] / rMag) * prT * (csi - 1.) / csi;
		}
		memcpy(p.data, tmpP, sizeof(tmpP));
	}
#endif
	//printf("Hamiltonian = %e\n", XLALSimIMRSpinEOBHamiltonian(eobParams->eta, &r, &p, &sigmaKerr, &sigmaStar, dParams->params->seobCoeffs));
	REAL8		SpinEOBH = XLALSimIMRSpinEOBHamiltonian(eobParams->eta, &r, &p, &spin1norm, &spin2norm, &sigmaKerr, &sigmaStar, dParams->params->tortoise, dParams->params->seobCoeffs) / eobParams->eta;

	if (dParams->varyParam < 3)
		dParams->params->tortoise = oldTortoise;
	return SpinEOBH;
}


/*
 * Wrapper for GSL to call the Hamiltonian function
 */
static double 
GSLSpinHamiltonianWrapperV2(double x, void *params)
{
	bool		UsePrec = false;
	HcapDerivParams *dParams = (HcapDerivParams *) params;

	EOBParams      *eobParams = dParams->params->eobParams;

	REAL8		tmpVec   [12];
	REAL8		s1normData[3], s2normData[3], sKerrData[3], sStarData[3];

	/*
	 * These are the vectors which will be used in the call to the
	 * Hamiltonian
	 */
	REAL8Vector	r  , p, spin1, spin2, spin1norm, spin2norm;
	REAL8Vector	sigmaKerr, sigmaStar;

	int		i;
	REAL8		a;
	REAL8		m1 = eobParams->m1;
	REAL8		m2 = eobParams->m2;
	REAL8		mT2 = (m1 + m2) * (m1 + m2);

	/* Use a temporary vector to avoid corrupting the main function */
	memcpy(tmpVec, dParams->values, sizeof(tmpVec));

	/* Set the relevant entry in the vector to the correct value */
	tmpVec[dParams->varyParam] = x;

	/* Set the LAL-style vectors to point to the appropriate things */
	r.length = p.length = spin1.length = spin2.length = spin1norm.length = spin2norm.length = 3;
	sigmaKerr.length = sigmaStar.length = 3;
	r.data = tmpVec;
	p.data = tmpVec + 3;
	spin1.data = tmpVec + 6;
	spin2.data = tmpVec + 9;
	spin1norm.data = s1normData;
	spin2norm.data = s2normData;
	sigmaKerr.data = sKerrData;
	sigmaStar.data = sStarData;

	memcpy(s1normData, tmpVec + 6, 3 * sizeof(REAL8));
	memcpy(s2normData, tmpVec + 9, 3 * sizeof(REAL8));

	for (i = 0; i < 3; i++) {
		s1normData[i] /= mT2;
		s2normData[i] /= mT2;
	}

	/* Calculate various spin parameters */
	XLALSimIMRSpinEOBCalculateSigmaKerr(&sigmaKerr, eobParams->m1,
					    eobParams->m2, &spin1, &spin2);
	XLALSimIMRSpinEOBCalculateSigmaStar(&sigmaStar, eobParams->m1,
					    eobParams->m2, &spin1, &spin2);
	a = sqrt(sigmaKerr.data[0] * sigmaKerr.data[0] + sigmaKerr.data[1] * sigmaKerr.data[1]
		 + sigmaKerr.data[2] * sigmaKerr.data[2]);
	//printf("a = %e\n", a);
	//printf("aStar = %e\n", sqrt(sigmaStar.data[0] * sigmaStar.data[0] + sigmaStar.data[1] * sigmaStar.data[1] + sigmaStar.data[2] * sigmaStar.data[2]));
	if (isnan(a)) {
		  printf("a is nan!!\n");
          XLALPrintError( "XLAL Error - %s: a = nan  \n", __func__);
          XLAL_ERROR( XLAL_EINVAL );
	}
	//XLALSimIMRCalculateSpinEOBHCoeffs(dParams->params->seobCoeffs, eobParams->eta, a);
	if (UsePrec) {
		/* Set up structures and calculate necessary PN parameters */
		/*
		 * Due to precession, these need to get calculated in every
		 * step
		 */
		/* TODO: Only calculate non-spinning parts once */
		memset(dParams->params->seobCoeffs, 0, sizeof(SpinEOBHCoeffs));

		/*
		 * Update the Z component of the Kerr background spin
		 * Pre-compute the Hamiltonian coefficients
		 */
		REAL8Vector    *delsigmaKerr = dParams->params->sigmaKerr;
		dParams->params->sigmaKerr = &sigmaKerr;
		dParams->params->a = a;
		//tmpsigmaKerr->data[2];
		if (XLALSimIMRCalculateSpinEOBHCoeffs(dParams->params->seobCoeffs,
						      eobParams->eta, a,
						      dParams->params->seobCoeffs->SpinAlignedEOBversion) == XLAL_FAILURE) {
			XLALDestroyREAL8Vector(dParams->params->sigmaKerr);
			XLAL_ERROR(XLAL_EFUNC);
		}
		/* Release the old memory */
		XLALDestroyREAL8Vector(delsigmaKerr);
	}
	//printf("Hamiltonian = %e\n", XLALSimIMRSpinEOBHamiltonian(eobParams->eta, &r, &p, &sigmaKerr, &sigmaStar, dParams->params->seobCoeffs));
	return XLALSimIMRSpinEOBHamiltonian(eobParams->eta, &r, &p, &spin1norm, &spin2norm, &sigmaKerr, &sigmaStar, dParams->params->tortoise, dParams->params->seobCoeffs) / eobParams->eta;
}


#endif				/* _LALSIMIMRSPINEOBHCAPNUMERICALDERIVATIVE_C */
