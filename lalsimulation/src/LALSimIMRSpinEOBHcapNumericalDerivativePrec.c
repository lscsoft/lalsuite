/*
 *  Copyright (C) 2014 Craig Robinson, Enrico Barausse, Yi Pan, Prayush Kumar, 
 *  Stanislav Babak, Andrea Taracchini
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
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */
/**
 * \author Craig Robinson, Yi Pan, Stas Babak, Prayush Kumar, Andrea Taracchini
 * \brief Functions to compute the RHSs of the SEOBNRv3 ODEs
 */

#ifndef _LALSIMIMRSPINPRECEOBHCAPNUMERICALDERIVATIVE_C
#define _LALSIMIMRSPINPRECEOBHCAPNUMERICALDERIVATIVE_C

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/* #include "LALSimIMRSpinEOBHcapNumericalDerivative.h" */

#include <math.h>
#include <gsl/gsl_deriv.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

#include "LALSimIMRSpinEOB.h"

#include "LALSimIMRSpinEOBAuxFuncs.c"
#include "LALSimIMRSpinEOBHamiltonian.c"
#include "LALSimIMRSpinEOBHamiltonianPrec.c"
#include "LALSimIMRSpinEOBFactorizedFluxPrec.c"
#include "LALSimIMRSpinEOBFactorizedWaveformCoefficientsPrec.c"
#include "LALSimIMREOBFactorizedWaveform.c"

// int		UsePrec = 1;

/*------------------------------------------------------------------------------------------
 *
 *          Prototypes of functions defined in this code.
 *
 *------------------------------------------------------------------------------------------
 */
static REAL8 GSLSpinPrecHamiltonianWrapper(double x, void *params);

static int
XLALSpinPrecHcapNumericalDerivative(
				double UNUSED t,	/**<< UNUSED */
				const REAL8 values[],	/**<< Dynamical variables */
				REAL8 dvalues[],	/**<< Time derivatives of variables (returned) */
				void *funcParams	/**<< EOB parameters */
);

static REAL8 XLALSpinPrecHcapNumDerivWRTParam(
					const		INT4	paramIdx,	/**<< Index of the parameters */
					const		REAL8	values[],	/**<< Dynamical variables */
						SpinEOBParams * funcParams	/**<< EOB Parameters */
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
 * All derivatives are returned together.
 * The derivatives are combined with energy flux to give right hand side of the ODEs
 * of a generic spin EOB model, as decribed in Eqs. A4, A5, 26 and 27 of
 * Pan et al. PRD 81, 084041 (2010)
 * Later on we use BB1, that is Barausse and Buonanno PRD 81, 084024 (2010)
 * This function is not used by the spin-aligned model.
 */
	static int	XLALSpinPrecHcapNumericalDerivative(
							double	UNUSED	t,	/**<< UNUSED */
					const		REAL8	values[],	/**<< Dynamical variables */
						REAL8		dvalues[],	/**<< Time derivatives of variables (returned) */
						void         *funcParams	/**<< EOB parameters */
)
{
	int		debugPK = 0;
	static const REAL8 STEP_SIZE = 2.0e-4;

    /** lMax: l index up to which h_{lm} modes are included in the computation of the GW enegy flux: see Eq. in 13 in PRD 86,  024011 (2012) */
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
	UNUSED REAL8 	rrTerm2, pDotS1, pDotS2;
	REAL8Vector	s1 , s2, s1norm, s2norm, sKerr, sStar;
	REAL8		s1Data   [3], s2Data[3], s1DataNorm[3], s2DataNorm[3];
	REAL8		sKerrData[3], sStarData[3];
	REAL8 chiS, chiA, a, tplspin;
	REAL8 	s1dotLN, s2dotLN;


	/* Orbital angular momentum */
	REAL8		Lx      , Ly, Lz, magL;
	REAL8		Lhatx   , Lhaty, Lhatz;
	REAL8		dLx     , dLy, dLz;
	REAL8		dLhatx  , dLhaty, dMagL;

	REAL8		alphadotcosi;

	REAL8		rCrossV_x, rCrossV_y, rCrossV_z, omega;

	/* The error in a derivative as measured by GSL */
	REAL8		absErr;

	REAL8		tmpP[3], rMag, rMag2, prT;
	REAL8		u, u2, u3, u4, u5, w2, a2;
	REAL8		D, m1PlusetaKK, bulk, logTerms, deltaU, deltaT, deltaR;
    REAL8  	eobD_r, deltaU_u, deltaU_r;
	REAL8		dcsi, csi;

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

	F.function = &GSLSpinPrecHamiltonianWrapper;
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
    /* Here pvec is the reduced tortoise p^* vector of Pan et al. PRD 81, 084041 (2010) */
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
		s1.data[i] *= (mass1 + mass2) * (mass1 + mass2);
		s2.data[i] *= (mass1 + mass2) * (mass1 + mass2);
	}

	sKerr.length = 3;
	sKerr.data = sKerrData;
	XLALSimIMRSpinEOBCalculateSigmaKerr(&sKerr, mass1, mass2, &s1, &s2);

	sStar.length = 3;
	sStar.data = sStarData;
	XLALSimIMRSpinEOBCalculateSigmaStar(&sStar, mass1, mass2, &s1, &s2);

	a = sqrt(sKerr.data[0] * sKerr.data[0] + sKerr.data[1] * sKerr.data[1]
		 + sKerr.data[2] * sKerr.data[2]);

	/* Convert momenta to p Eq. A3 of PRD 81, 084041 (2010) */
	rMag = sqrt(rVec.data[0] * rVec.data[0] + rVec.data[1] * rVec.data[1] + rVec.data[2] * rVec.data[2]);
    /* This is p^*.r/|r| */
	prT = pVec.data[0] * (rVec.data[0] / rMag) + pVec.data[1] * (rVec.data[1] / rMag)
		+ pVec.data[2] * (rVec.data[2] / rMag);

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
     /* d(Eq. 5.83 of BB1)/dr */
	eobD_r = (u2 / (D * D)) * (12. * eta * u + 6. * (26. - 3. * eta) * eta * u2) / (1.
			+ 6. * eta * u2 + 2. * (26. - 3. * eta) * eta * u3);
	m1PlusetaKK = -1. + eta * coeffs->KK;
	/* Eq. 5.75 of BB1 */
    /* a as returned by XLALSimIMRSpinEOBCalculateSigmaKerr is S/M^2 so that a is unitless, i.e. 1/M^2 is absorbed in a2. Also, u = M/|r| is unitless */
	bulk = 1. / (m1PlusetaKK * m1PlusetaKK) + (2. * u) / m1PlusetaKK + a2 * u2;
	/* Eq. 5.73 of BB1 */
    /* The 4PN term with coefficients k5 and k5l are defined in the SEOBNRv2 review document https://dcc.ligo.org/T1400476 */
	logTerms = 1. + eta * coeffs->k0 + eta * log(1. + coeffs->k1 * u
		       + coeffs->k2 * u2 + coeffs->k3 * u3 + coeffs->k4 * u4
			     + coeffs->k5 * u5 + coeffs->k5l * u5 * log(u));
	/* Eq. 5.73 of BB1 */
	deltaU = bulk * logTerms;
	/* Eq. 5.71 of BB1 */
	deltaT = rMag2 * deltaU;
	/* ddeltaU/du */
    /* The log(u) is treated as a constant when taking the derivative wrt u */
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

    /* This is A3 of PRD 81, 084041 (2010) explicitly evaluated */
	for (i = 0; i < 3; i++) {
		tmpP[i] = pVec.data[i] - (rVec.data[i] / rMag) * prT * (csi - 1.) / csi;
	}

	if (debugPK) {
		XLAL_PRINT_INFO("csi = %.12e\n", csi);
		for (i = 0; i < 3; i++)
			XLAL_PRINT_INFO("p,p*: %.12e\t%.12e\n", pData[i], tmpP[i]);
	}
	/*
	 * Calculate the T-matrix, required to convert P from tortoise to
	 * non-tortoise coordinates, and/or vice-versa. This is given
	 * explicitly in Eq. A3 of 0912.3466
	 */
	for (i = 0; i < 3; i++)
		for (j = 0; j <= i; j++) {
			Tmatrix[i][j] = Tmatrix[j][i] = (rVec.data[i] * rVec.data[j] / rMag2)
				* (csi - 1.);

			invTmatrix[i][j] = invTmatrix[j][i] =
				-(csi - 1) / csi * (rVec.data[i] * rVec.data[j] / rMag2);

			if (i == j) {
				Tmatrix[i][j]++;
				invTmatrix[i][j]++;
			}
		}

     /* This is dcsi/dr: this is needed in the last term of Eq. A5 of PRD 81, 084041 (2010) */
	dcsi = csi * (2. / rMag + deltaU_r / deltaU) + csi * csi * csi
		/ (2. * rMag2 * rMag2 * deltaU * deltaU) * (rMag * (-4. * w2) / D - eobD_r * (w2 * w2));

	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			for (k = 0; k < 3; k++) {
				dTijdXk[i][j][k] =
					(rVec.data[i] * KRONECKER_DELTA(j, k) + KRONECKER_DELTA(i, k) * rVec.data[j])
					* (csi - 1.) / rMag2
					+ rVec.data[i] * rVec.data[j] * rVec.data[k] / rMag2 / rMag * (-2. / rMag * (csi - 1.) + dcsi);
			}

	//Print out the T - matrix for comparison
		if (debugPK) {
			XLAL_PRINT_INFO("\nT-Matrix:\n");
			for (i = 0; i < 3; i++)
				XLAL_PRINT_INFO("%le\t%le\t%le\n", Tmatrix[i][0], Tmatrix[i][1], Tmatrix[i][2]);

			for (i = 0; i < 3; i++) {
				XLAL_PRINT_INFO("dT[%d][j]dX[k]:\n", i);
				for (j = 0; j < 3; j++)
					XLAL_PRINT_INFO("%.12e\t%.12e\t%.12e\n", dTijdXk[i][j][0],
					dTijdXk[i][j][1], dTijdXk[i][j][2]);
				XLAL_PRINT_INFO("\n");
			}
		}
    INT4   updateHCoeffsOld =  params.params->seobCoeffs->updateHCoeffs;
	/* Now calculate derivatives w.r.t. each parameter */
	for (i = 0; i < 3; i++) {
		params.varyParam = i;
        params.params->tortoise = 2;
        memcpy(tmpValues, params.values, sizeof(tmpValues));
        tmpValues[3] = tmpP[0];
        tmpValues[4] = tmpP[1];
        tmpValues[5] = tmpP[2];
        params.values = tmpValues;
        /* We need derivatives of H wrt to P (and not P^*) */
        /* Note that in the 1st term on the last line of Eq. A5 of PRD 81, 084041 (2010) one needs
         * dH/dX @ fixed P, not P^*, hence the need for what follows  */
        XLAL_CALLGSL(gslStatus = gsl_deriv_central(&F, values[i],
        STEP_SIZE, &tmpDValues[i], &absErr));
        params.values = values;
        params.params->tortoise = 1;
        if (gslStatus != GSL_SUCCESS) {
            XLALPrintError("XLAL Error %s - Failure in GSL function\n", __func__);
            XLAL_ERROR(XLAL_EFUNC);
        }
    }
    for (i = 3; i < 6; i++) {
        params.varyParam = i;
        params.params->seobCoeffs->updateHCoeffs = 1;
        XLAL_CALLGSL(gslStatus = gsl_deriv_central(&F, values[i],
        STEP_SIZE, &tmpDValues[i], &absErr));
        if (gslStatus != GSL_SUCCESS) {
            XLALPrintError("XLAL Error %s - Failure in GSL function\n", __func__);
            XLAL_ERROR(XLAL_EFUNC);
        }
    }
    for (i = 6; i < 9; i++) {
        params.varyParam = i;
        params.params->seobCoeffs->updateHCoeffs = 1;
        XLAL_CALLGSL(gslStatus = gsl_deriv_central(&F, values[i],
        STEP_SIZE * mass1 * mass1, &tmpDValues[i], &absErr));
        if (gslStatus != GSL_SUCCESS) {
            XLALPrintError("XLAL Error %s - Failure in GSL function\n", __func__);
            XLAL_ERROR(XLAL_EFUNC);
        }
    }
    for (i = 9; i < 12; i++) {
        params.varyParam = i;
        params.params->seobCoeffs->updateHCoeffs = 1;
        XLAL_CALLGSL(gslStatus = gsl_deriv_central(&F, values[i],
        STEP_SIZE * mass2 * mass2, &tmpDValues[i], &absErr));
        if (gslStatus != GSL_SUCCESS) {
            XLALPrintError("XLAL Error %s - Failure in GSL function\n", __func__);
            XLAL_ERROR(XLAL_EFUNC);
        }
    }
    params.params->seobCoeffs->updateHCoeffs = updateHCoeffsOld;

	/* Now make the conversion */
	/* rVectorDot */
    /* Eq. A4 of PRD 81, 084041 (2010).  Note that dvalues[i] = \dot{X^i} but tmpDValues[j] = dH/dvalues[j] */
	for (i = 0; i < 3; i++)
		for (j = 0, dvalues[i] = 0.; j < 3; j++)
			dvalues[i] += tmpDValues[j + 3] * Tmatrix[i][j];

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

	r = polarDynamics.data[0] = sqrt(values[0] * values[0] + values[1] * values[1]
			      + values[2] * values[2]);
	polarDynamics.data[1] = 0;
    /* Note that poldata[2] and poldata[3] differ from v2 in which one is normalized by polData[0]. Behaviour is reverted. */
	polarDynamics.data[2] = (values[0] * values[3] + values[1] * values[4]
		      + values[2] * values[5]) / polData[0];
	polarDynamics.data[3] = magL;


	/* Now calculate rCrossRdot and omega */
	rCrossV_x = values[1] * dvalues[2] - values[2] * dvalues[1];
	rCrossV_y = values[2] * dvalues[0] - values[0] * dvalues[2];
	rCrossV_z = values[0] * dvalues[1] - values[1] * dvalues[0];

	omega = sqrt(rCrossV_x * rCrossV_x + rCrossV_y * rCrossV_y + rCrossV_z * rCrossV_z) / (r * r);

	/*
	 * Compute \vec{L_N} = \vec{r} \times \.{\vec{r}}, \vec{S_i} \dot
	 * \vec{L_N} and chiS and chiA
	 */
    
    /* Eq. 16 of PRD 89, 084006 (2014): it's S_{1,2}/m_{1,2}^2.LNhat */
	s1dotLN = (s1.data[0] * rCrossV_x + s1.data[1] * rCrossV_y + s1.data[2] * rCrossV_z) /
		(r * r * omega * mass1 * mass1);
	s2dotLN = (s2.data[0] * rCrossV_x + s2.data[1] * rCrossV_y + s2.data[2] * rCrossV_z) /
		(r * r * omega * mass2 * mass2);

	chiS = 0.5 * (s1dotLN + s2dotLN);
	chiA = 0.5 * (s1dotLN - s2dotLN);

	if (debugPK) {
		XLAL_PRINT_INFO("chiS = %.12e, chiA = %.12e\n", chiS, chiA);
		fflush(NULL);
	}
	if (debugPK) {
		XLAL_PRINT_INFO("Computing derivatives for values\n%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n\n",
		    (double)values[0], (double)values[1], (double)values[2],
		    (double)values[3], (double)values[4], (double)values[5],
            (double)values[6], (double)values[7], (double)values[8],
            (double)values[9], (double)values[10], (double)values[11]);
		XLAL_PRINT_INFO("tmpDvalues\n%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t\n",
		       (double)tmpDValues[0], (double)tmpDValues[1], (double)tmpDValues[2],
		       (double)tmpDValues[3], (double)tmpDValues[4], (double)tmpDValues[5],
		       (double)tmpDValues[6], (double)tmpDValues[7], (double)tmpDValues[8],
		       (double)tmpDValues[9], (double)tmpDValues[10], (double)tmpDValues[11]);
	}
	/*
	 * Compute the test-particle limit spin of the deformed-Kerr
	 * background
	 */
	switch (SpinAlignedEOBversion) {
	case 1:
        /* See below Eq. 17 of PRD 86, 041011 (2012) */
		tplspin = 0.0;
		break;
	case 2:
        /* See below Eq. 4 of PRD 89, 061502(R) (2014)*/
		tplspin = (1. - 2. * eta) * chiS + (mass1 - mass2) / (mass1 + mass2) * chiA;
		break;
	default:
		XLALPrintError("XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
		XLAL_ERROR(XLAL_EINVAL);
		break;
	}


    memcpy(params.params->s1Vec->data, s1norm.data, 3*sizeof(*params.params->s1Vec->data));
    memcpy(params.params->s2Vec->data, s2norm.data, 3*sizeof(*params.params->s2Vec->data));
    memcpy(params.params->sigmaStar->data, sStar.data, 3*sizeof(*params.params->sigmaStar->data));
    memcpy(params.params->sigmaKerr->data, sKerr.data, 3*sizeof(*params.params->sigmaKerr->data));


	params.params->a = a;
    if (params.params->alignedSpins==1) {
        XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(
	      params.params->eobParams->hCoeffs, mass1, mass2, eta, tplspin,
					 chiS, chiA, SpinAlignedEOBversion);
    }
    else {
        XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(
                                                     params.params->eobParams->hCoeffs, mass1, mass2, eta, tplspin,
                                                     chiS, chiA, 3);
    }
	XLALSimIMRCalculateSpinPrecEOBHCoeffs(params.params->seobCoeffs, eta, a,
					  SpinAlignedEOBversion);

	H = XLALSimIMRSpinPrecEOBHamiltonian(eta, &rVec, &pVec, &s1norm, &s2norm,
	&sKerr, &sStar, params.params->tortoise, params.params->seobCoeffs);
	H = H * (mass1 + mass2);

	/* Now we have the ingredients to compute the flux */
	memcpy(tmpValues, values, 12 * sizeof(REAL8));
	cartDynamics.data = tmpValues;
	if (debugPK) {
		XLAL_PRINT_INFO("params.params->a = %.12e, %.12e\n", a, params.params->a);
		fflush(NULL);
	}
    /* Eq. 13 of PRD 86, 024011 (2012) */
    if ( params.params->ignoreflux == 0 ) {
        flux = XLALInspiralPrecSpinFactorizedFlux(&polarDynamics, &cartDynamics,
						  nqcCoeffs, omega, params.params, H / (mass1 + mass2), lMax, SpinAlignedEOBversion);
    }
    else if ( params.params->ignoreflux  == 1) {
            flux = 0.;
    }
    else {
        XLAL_PRINT_INFO("Wrong ignorflux option in XLALSpinPrecHcapNumericalDerivative!\n");
        XLAL_ERROR(XLAL_EFUNC);
    }
	/*
	 * Looking at consistency with the non-spinning model, we have to divide the
	 * flux by eta
	 */
	flux = flux / eta;

	pDotS1 = pData[0] * s1.data[0] + pVec.data[1] * s1.data[1] + pVec.data[2] * s1.data[2];
	pDotS2 = pVec.data[0] * s2.data[0] + pVec.data[1] * s2.data[1] + pVec.data[2] * s2.data[2];
	rrTerm2 = 8. / 15. * eta * eta * pow(omega, 8. / 3.) / (magL * magL * r) * ((61. + 48. * mass2 / mass1) * pDotS1 + (61. + 48. * mass1 / mass2) * pDotS2);

	if (debugPK) {
		XLAL_PRINT_INFO("omega = %.12e \n flux = %.12e \n Lmag = %.12e\n", omega, flux, magL);
		XLAL_PRINT_INFO("rrForce = %.12e %.12e %.12e\n", -flux * values[3] / (omega * magL), -flux * values[4] / (omega * magL), -flux * values[5] / (omega * magL));
	}
	/* Now pDot */
	/* Compute the first and second terms in eq. A5 of PRD 81, 084041 (2010) */
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
					XLAL_PRINT_INFO("\n sum[%d][%d][%d] = %.12e", i, j, l, sum);

				}
		XLAL_PRINT_INFO("\n\n Printing dTdX * Tinverse:\n");
		for (l = 0; l < 3; l++) {
			for (i = 0; i < 3; i++)
				for (j = 0; j < 3; j++) {
					double		sum = 0;
					for (k = 0; k < 3; k++) {
						sum += dTijdXk[i][k][l] * invTmatrix[k][j];
						XLAL_PRINT_INFO("\n sum[%d][%d][%d] = %.12e", l, i, j, sum);
					}
				}
		}
	}
	if (debugPK)
		XLAL_PRINT_INFO("\npData: {%.12e, %.12e, %.12e}\n", pData[0], pData[1], pData[2]);
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
		XLAL_PRINT_INFO("\ntmpPdotT3 = ");
		for (i = 0; i < 3; i++)
			XLAL_PRINT_INFO("%.12e ", tmpPdotT3[i]);
		XLAL_PRINT_INFO("\n");
	}

    /* Eqs. 11c-11d of PRD 89, 084006 (2014) */
	/* spin1 */
	if (debugPK) {
		XLAL_PRINT_INFO("Raw spin1 derivatives = %.12e %.12e %.12e\n", tmpDValues[6], tmpDValues[7], tmpDValues[8]);
		XLAL_PRINT_INFO("Raw spin2 derivatives = %.12e %.12e %.12e\n", tmpDValues[9], tmpDValues[10], tmpDValues[11]);
	}
     /* The factor eta is there becasue Hreal is normalized by eta */
	dvalues[6] = eta * (tmpDValues[7] * values[8] - tmpDValues[8] * values[7]);
	dvalues[7] = eta * (tmpDValues[8] * values[6] - tmpDValues[6] * values[8]);
	dvalues[8] = eta * (tmpDValues[6] * values[7] - tmpDValues[7] * values[6]);

	/* spin2 */
    /* The factor eta is there becasue Hreal is normalized by eta */
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
	 * Finn Chernoff convention is used here. 
     */
    /* Eqs. 19-20 of PRD 89, 084006 (2014) */
	if (Lhatx == 0.0 && Lhaty == 0.0) {
		alphadotcosi = 0.0;
	} else {
		alphadotcosi = Lhatz * (Lhatx * dLhaty - Lhaty * dLhatx) / (Lhatx * Lhatx + Lhaty * Lhaty);
	}
    /* These are ODEs for the phase that enters the h_{lm}: see Eq. 3.11 of PRD 79, 104023 (2009) */
	dvalues[12] = omega - alphadotcosi;
	dvalues[13] = alphadotcosi;

	if (debugPK) {
    XLAL_PRINT_INFO("\nIn XLALSpinPrecHcapNumericalDerivative:\n");
		/* Print out all mass parameters */
		XLAL_PRINT_INFO("m1 = %12.12lf, m2 = %12.12lf, eta = %12.12lf\n",
          (double)mass1, (double)mass2, (double)eta);
		/* Print out all spin parameters */
		XLAL_PRINT_INFO("spin1 = {%12.12lf,%12.12lf,%12.12lf}, spin2 = {%12.12lf,%12.12lf,%12.12lf}\n",
            (double)s1.data[0], (double)s1.data[1], (double)s1.data[2],
            (double)s2.data[0], (double)s2.data[1], (double)s2.data[2]);
		XLAL_PRINT_INFO("sigmaStar = {%12.12lf,%12.12lf,%12.12lf}, sigmaKerr = {%12.12lf,%12.12lf,%12.12lf}\n",
		       (double)sStar.data[0], (double)sStar.data[1],
		       (double)sStar.data[2], (double)sKerr.data[0],
		       (double)sKerr.data[1], (double)sKerr.data[2]);
		XLAL_PRINT_INFO("L = {%12.12lf,%12.12lf,%12.12lf}, |L| = %12.12lf\n",
        (double)Lx, (double)Ly, (double)Lz, (double)magL);
		XLAL_PRINT_INFO("dLdt = {%12.12lf,%12.12lf,%12.12lf}, d|L|dt = %12.12lf\n",
        (double)dLx, (double)dLy, (double)dLz, (double)dMagL);
		XLAL_PRINT_INFO("Polar coordinates = {%12.12lf, %12.12lf, %12.12lf, %12.12lf}\n",
		       (double)polData[0], (double)polData[1], (double)polData[2],
           (double)polData[3]);

		XLAL_PRINT_INFO("Cartesian coordinates: {%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf}\n",
		       (double)values[0], (double)values[1], (double)values[2],
           (double)values[3], (double)values[4], (double)values[5],
           (double)values[6], (double)values[7], (double)values[8],
           (double)values[9], (double)values[10], (double)values[11],
		       (double)values[12], (double)values[13]);
		XLAL_PRINT_INFO("Cartesian derivatives: {%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf,%12.12lf}\n",
		       (double)dvalues[0], (double)dvalues[1], (double)dvalues[2],
           (double)dvalues[3], (double)dvalues[4], (double)dvalues[5],
           (double)dvalues[6], (double)dvalues[7], (double)dvalues[8],
           (double)dvalues[9], (double)dvalues[10], (double)dvalues[11],
           (double)dvalues[12], (double)dvalues[13]);

     XLAL_PRINT_INFO("Hamiltonian = %12.12lf, Flux = %12.12lf, Omega = %12.12lf\n",
              H/ (mass1 + mass2), eta*flux, omega);
		fflush(NULL);
	}

    if(debugPK){
  for(i=0; i<14; i++)
    if(dvalues[i] > 1e3){
      XLAL_PRINT_INFO("\nIn XLALSpinPrecHcapNumericalDerivative:\n");
      XLAL_PRINT_INFO("Derivatives have blown up!\n");
      for(j=0; j<14; XLAL_PRINT_INFO("dvalues[%d] = %3.12f\n", j, dvalues[j]), j++);
      XLAL_PRINT_INFO("Flux = %3.12f\n\n", flux);
      break;
    }
    }
  return XLAL_SUCCESS;
}


/**
 * Calculate the derivative of the Hamiltonian w.r.t. a specific parameter
 * Used by generic spin EOB model, including initial conditions solver.
 */
static REAL8
XLALSpinPrecHcapNumDerivWRTParam(
			     const INT4 paramIdx,	/**<< Index of the parameters */
			     const REAL8 values[],	/**<< Dynamical variables */
			     SpinEOBParams * funcParams	/**<< EOB Parameters */
)
{
	static const REAL8 STEP_SIZE = 2.0e-3;

	HcapDerivParams	params;

	REAL8		result;

	gsl_function	F;
	INT4		gslStatus;

	REAL8		mass1   , mass2;

	/* The error in a derivative as measured by GSL */
	REAL8		absErr;

	/* Set up pointers for GSL */
	params.params = funcParams;
    

	F.function = &GSLSpinPrecHamiltonianWrapper;
	F.params = &params;
	params.varyParam = paramIdx;

	mass1 = params.params->eobParams->m1;
	mass2 = params.params->eobParams->m2;
    
    REAL8 mT2 = (mass1 + mass2) * (mass1 + mass2);
    REAL8 tmpValues[14];
    for (UINT4 i = 0; i < 3; i++) {
        tmpValues[i] = values[i];
        tmpValues[i + 3] = values[i + 3];
        tmpValues[i + 6] = values[i + 6]/mT2;
        tmpValues[i + 9] = values[i + 9]/mT2;
    }
    tmpValues[12] = values[12];
    tmpValues[13] = values[13];
    params.values = tmpValues;


	/* Now calculate derivatives w.r.t. the required parameter */
    if (paramIdx >=0 && paramIdx < 6) {
        XLAL_CALLGSL(gslStatus = gsl_deriv_central(&F, values[paramIdx],
                                                   STEP_SIZE, &result, &absErr));
    }
    else if (paramIdx >= 6 && paramIdx < 9) {
		XLAL_CALLGSL(gslStatus = gsl_deriv_central(&F, values[paramIdx],
			      STEP_SIZE * mass1 * mass1, &result, &absErr));
	}
    else if (paramIdx >= 9 && paramIdx < 12) {
		XLAL_CALLGSL(gslStatus = gsl_deriv_central(&F, values[paramIdx],
			      STEP_SIZE * mass2 * mass2, &result, &absErr));
	}
    else {
        XLAL_PRINT_INFO("Requested partial derivative is not availabble\n");
        XLAL_ERROR( XLAL_EINVAL );
    }
	if (gslStatus != GSL_SUCCESS) {
		XLALPrintError("XLAL Error %s - Failure in GSL function\n", __func__);
		XLAL_ERROR_REAL8(XLAL_EFUNC);
	}
	//XLAL_PRINT_INFO("Abserr = %e\n", absErr);

	return result;
}


/**
 * Wrapper for GSL to call the Hamiltonian function
 */
static REAL8
GSLSpinPrecHamiltonianWrapper(double x, void *params)
{
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

	memcpy(spin1norm.data, tmpVec + 6, 3 * sizeof(REAL8));
	memcpy(spin2norm.data, tmpVec + 9, 3 * sizeof(REAL8));

	/*
	 * To compute the SigmaKerr and SigmaStar, we need the non-normalized
	 * spin values, i.e. S_i. The spins being evolved are S_i/M^2.
	 */
	for (i = 0; i < 3; i++) {
		spin1.data[i] *= mT2;
		spin2.data[i] *= mT2;
	}

	/* Calculate various spin parameters */
    /* Note that XLALSimIMRSpinEOBCalculateSigmaKerr and XLALSimIMRSpinEOBCalculateSigmaStar return a unitless quantity */
	XLALSimIMRSpinEOBCalculateSigmaKerr(&sigmaKerr, eobParams->m1,
					    eobParams->m2, &spin1, &spin2);
	XLALSimIMRSpinEOBCalculateSigmaStar(&sigmaStar, eobParams->m1,
					    eobParams->m2, &spin1, &spin2);
	a = sqrt(sigmaKerr.data[0] * sigmaKerr.data[0]
		 + sigmaKerr.data[1] * sigmaKerr.data[1]
		 + sigmaKerr.data[2] * sigmaKerr.data[2]);

	if (isnan(a) ) {
          XLALPrintError( "XLAL Error - %s: a = nan  \n", __func__);
          XLAL_ERROR( XLAL_EINVAL );
	}

	REAL8		SpinEOBH = XLALSimIMRSpinPrecEOBHamiltonian(eobParams->eta, &r, &p, &spin1norm, &spin2norm, &sigmaKerr, &sigmaStar, dParams->params->tortoise, dParams->params->seobCoeffs) / eobParams->eta;

	if (dParams->varyParam < 3)
		dParams->params->tortoise = oldTortoise;
	return SpinEOBH;
}

#endif				/* _LALSIMIMRSPINPRECEOBHCAPNUMERICALDERIVATIVE_C */
