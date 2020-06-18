#include <lal/LALSimIMR.h>
#include <lal/LALSimInspiral.h>
#include <lal/Date.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/LALAdaptiveRungeKuttaIntegrator.h>
#include <lal/SphericalHarmonics.h>
#include <lal/LALSimSphHarmMode.h>
#include <LALSimInspiralWaveformFlags.h>
#include <lal/LALDict.h>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_deriv.h>

#include "LALSimInspiralEOBPostAdiabatic.h"
#include "LALSimIMREOBNRv2.h"
#include "LALSimInspiralPrecess.h"
#include "LALSimBlackHoleRingdown.h"
#include "LALSimIMRSpinEOB.h"

#include "LALSimIMRSpinEOBHamiltonian.c"
#include "LALSimIMRSpinEOBFactorizedFlux.c"
#include "LALSimIMREOBNewtonianMultipole.c"
// #include "LALSimIMREOBNQCTables.c"

int
XLALSimInspiralEOBPACalculateRadialGrid(
	REAL8Vector *rVec,
	LALDict *LALParams)
{
	const REAL8 rInitial = XLALDictLookupREAL8Value(LALParams, "rInitial");
	const UINT4 rSize = XLALDictLookupUINT4Value(LALParams, "rSize");
	const REAL8 dr = XLALDictLookupREAL8Value(LALParams, "dr");

	UINT4 i;

	for (i = 0; i < rSize; i++)
	{
		rVec->data[i] = rInitial - i*dr;
	}

	return XLAL_SUCCESS;
}

double
XLALSimInspiralEOBPostAdiabaticdprstarFunc(
	REAL8 prstar_sol,
	void *params)
{
	struct PostAdiabaticRootSolveParams *prstarParams = (struct PostAdiabaticRootSolveParams *)params;

	REAL8 r = prstarParams->r;
	REAL8 prstar = prstar_sol;
	REAL8 pphi = prstarParams->pphi;
	REAL8 omega = prstarParams->omega;
	REAL8 dprstarBydpr = prstarParams->csi;
	SpinEOBParams *seobParams = prstarParams->seobParams;
	EOBNonQCCoeffs *nqcCoeffs = prstarParams->nqcCoeffs;

	SpinEOBHCoeffs *seobCoeffs = seobParams->seobCoeffs;

	REAL8 dpphiBydr;
	dpphiBydr = prstarParams->dpphiBydr;

	LALDict *LALParams = prstarParams->LALParams;

	REAL8 partialHBypartialprstar;
	partialHBypartialprstar = XLALSimInspiralEOBPAHamiltonianPartialDerivativeprstar(
																					1.e-5,
																					r,
																					prstar,
																					pphi,
																					seobCoeffs,
																					LALParams
																				);
	REAL8 H;
	H = XLALSimInspiralEOBPAHamiltonianWrapper(
							r,
							prstar,
							pphi,
							seobCoeffs,
							LALParams
						);

	REAL8 flux;
    flux = XLALSimInspiralEOBPAFluxWrapper(
				r,
				prstar,
				pphi,
				omega,
				H,
				seobParams,
				nqcCoeffs,
				LALParams
			);

	REAL8 result;

	result = dpphiBydr*partialHBypartialprstar*dprstarBydpr - flux;

	return result;
}

double 
XLALSimInspiralEOBPostAdiabaticdpphiFunc(
	REAL8 pphi_sol,
	void *params)
{
	struct PostAdiabaticRootSolveParams *pphiParams = (struct PostAdiabaticRootSolveParams *)params;

	REAL8 r = pphiParams->r;
	REAL8 prstar = pphiParams->prstar;
	REAL8 pphi = pphi_sol;
	SpinEOBHCoeffs *seobCoeffs = pphiParams->seobParams->seobCoeffs;
	LALDict *LALParams = pphiParams->LALParams;

	const REAL8 dr = XLALDictLookupREAL8Value(LALParams, "dr");

	REAL8 partialHBypartialprstar;
	partialHBypartialprstar = XLALSimInspiralEOBPAHamiltonianPartialDerivativeprstar(
									1.e-5,
									r,
									prstar,
									pphi,
									seobCoeffs,
									LALParams
								);

	REAL8 dprstarBydpr;
	dprstarBydpr = pphiParams->csi;

	REAL8 partialHBypartialr;
	partialHBypartialr = XLALSimInspiralEOBPAHamiltonianDerivative(
								dr,
								r,
								prstar,
								pphi,
								seobCoeffs,
								LALParams
							);

	REAL8 dprstarBydr;
	dprstarBydr = pphiParams->dprstarBydr;
	
	REAL8 result;
    result = partialHBypartialprstar*dprstarBydpr - partialHBypartialr/dprstarBydr;

	return result;
}

double
XLALSimInspiralEOBPostAdiabaticj0Func(
	REAL8 j0_sol,
	void *params)
{
	struct PostAdiabaticRootSolveParams *j0Params = (struct PostAdiabaticRootSolveParams *)params;

	SpinEOBHCoeffs *seobCoeffs = j0Params->seobParams->seobCoeffs;

	const REAL8 dr = XLALDictLookupREAL8Value(j0Params->LALParams, "dr");

	REAL8 r = j0Params->r;
   	REAL8 prstar = j0Params->prstar;
    LALDict *LALParams = j0Params->LALParams;

    REAL8 result;
    result = XLALSimInspiralEOBPAHamiltonianDerivative(dr, r, prstar, j0_sol, seobCoeffs, LALParams);

	return result;
}

int
XLALSimInspiralEOBPostAdiabaticRootFinder(
	struct PostAdiabaticRoot *result,
	double (*Func)(REAL8, void *),
	struct PostAdiabaticRootSolveParams *params,
	REAL8 x_lower,
	REAL8 x_upper,
	REAL8 absTol,
	REAL8 relTol)
{
	INT4 maxIters = 1000;
	INT4 iters;
	INT4 status;
	REAL8 x;

	iters = 0;
	x = 0.0;

	gsl_function F;

	F.function = Func;
	F.params = params;

    const gsl_root_fsolver_type *solver_type;
    // solver_type = gsl_root_fsolver_falsepos;
    solver_type = gsl_root_fsolver_brent;

	gsl_root_fsolver *solver;
	solver = gsl_root_fsolver_alloc(solver_type);

	REAL8 F_lower = Func(x_lower, params);
	REAL8 F_upper = Func(x_upper, params);

	if (F_lower*F_upper >= 0.0)
	{
		x_lower = -5.e-1;
		x_upper = 5.e-1;
	}

	gsl_root_fsolver_set(solver, &F, x_lower, x_upper);

	status = GSL_CONTINUE;

	do
	{
		iters++;

        /* iterate one step of the solver */
        status = gsl_root_fsolver_iterate(solver);

        if (status != GSL_SUCCESS)
            break;

        /* get the solver's current best solution and bounds */
        x = gsl_root_fsolver_root(solver);
        x_lower = gsl_root_fsolver_x_lower(solver);
        x_upper = gsl_root_fsolver_x_upper(solver);

        /* Check to see if the solution is within 0.001 */
        status = gsl_root_test_interval(x_lower, x_upper, absTol, relTol);
    }
    while (iters <= maxIters && status == GSL_CONTINUE);

    result->root = x;
    result->status = status;
    result->nIter = iters;

    return XLAL_SUCCESS;
}

REAL8Vector
XLALReverseREAL8Vector(
	REAL8Vector *Vec)
{
	UINT4 vecLength;
	vecLength = Vec->length;

	REAL8Vector *reverseVec = XLALCreateREAL8Vector(vecLength);
	memset(reverseVec->data, 0, reverseVec->length * sizeof(REAL8));

	UINT4 i;

	for (i = 0; i < vecLength; i++)
	{
		reverseVec->data[i] = Vec->data[vecLength-i-1];
	}

	return *reverseVec;
}

REAL8Vector
XLALOffsetREAL8Vector(
	REAL8Vector *Vec,
	REAL8 offset)
{
	UINT4 vecLength;
	vecLength = Vec->length;

	REAL8Vector *offsetVec = XLALCreateREAL8Vector(vecLength);
	memset(offsetVec->data, 0, offsetVec->length * sizeof(REAL8));

	UINT4 i;

	for (i = 0; i < vecLength; i++)
	{
		offsetVec->data[i] = Vec->data[i] + offset;
	}

	return *offsetVec;
}

int
XLALSimInspiralEOBPACalculateAdiabaticDynamics(
	REAL8Vector *rVec,
	REAL8Vector *phiVec,
	REAL8Vector *prstarVec,
	REAL8Vector *pphiVec,
	REAL8Vector *pphi0Vec,
	REAL8Vector *dpphiBydrVec,
	REAL8Vector *dpphiBydr0Vec,
	REAL8Vector *HVec,
	REAL8Vector *csiVec,
	REAL8Vector *omegaVec,
	SpinEOBParams *seobParams,
	LALDict *LALParams
)
{
	const UINT4 rSize = XLALDictLookupUINT4Value(LALParams, "rSize");
	const REAL8 nu = XLALDictLookupREAL8Value(LALParams, "nu");
	const REAL8 sigmaKerr = XLALDictLookupREAL8Value(LALParams, "aK");
	const REAL8 dr = XLALDictLookupREAL8Value(LALParams, "dr");

	UINT4 i;

	REAL8 DeltaT;
	REAL8 DeltaR;

	for (i = 0; i < rSize; i++)
	{
		REAL8 Newtonianj0;
		Newtonianj0 = XLALSimInspiralEOBPACalculateNewtonianj0(rVec->data[i]);

		struct PostAdiabaticRootSolveParams pphiParams;

        pphiParams.r = rVec->data[i];
        pphiParams.prstar = prstarVec->data[i];
        pphiParams.seobParams = seobParams;
	    pphiParams.LALParams = LALParams;

	    REAL8 pphi0_lower;
	    REAL8 pphi0_upper;

	    pphi0_lower = 0.1 * Newtonianj0;
		pphi0_upper = 1.9 * Newtonianj0;

		struct PostAdiabaticRoot pphiRoot;

	    XLALSimInspiralEOBPostAdiabaticRootFinder(
	    	&pphiRoot,
			XLALSimInspiralEOBPostAdiabaticj0Func,
			&pphiParams,
			pphi0_lower,
			pphi0_upper,
			1.e-17,
			1.e-20
		);

		pphiVec->data[i] = pphiRoot.root;
        pphi0Vec->data[i] = pphiRoot.root;

		HVec->data[i] = XLALSimInspiralEOBPAHamiltonianWrapper(
							rVec->data[i],
							prstarVec->data[i],
							pphiVec->data[i],
							seobParams->seobCoeffs,
							LALParams
						);

		DeltaT = XLALSimIMRSpinEOBHamiltonianDeltaT(seobParams->seobCoeffs, rVec->data[i], nu, sigmaKerr);
		DeltaR = XLALSimIMRSpinEOBHamiltonianDeltaR(seobParams->seobCoeffs, rVec->data[i], nu, sigmaKerr);

		csiVec->data[i] = sqrt(DeltaT*DeltaR) / (rVec->data[i]*rVec->data[i]+sigmaKerr*sigmaKerr);

		REAL8 polarDynamics[4];
		
		polarDynamics[0] = rVec->data[i];
		polarDynamics[1] = phiVec->data[i];
		polarDynamics[2] = prstarVec->data[i];
		polarDynamics[3] = pphiVec->data[i];

		omegaVec->data[i] = XLALSimIMRSpinAlignedEOBCalcOmega(
								polarDynamics,
								seobParams,
								dr
							);
	}

	REAL8Vector *rReverseVec = XLALCreateREAL8Vector(rSize);
	memset(rReverseVec->data, 0, rReverseVec->length * sizeof(REAL8));

	REAL8Vector *pphiReverseVec = XLALCreateREAL8Vector(rSize);
	memset(pphiReverseVec->data, 0, pphiReverseVec->length * sizeof(REAL8));

	REAL8Vector *dpphiBydrReverseVec = XLALCreateREAL8Vector(rSize);
	memset(dpphiBydrReverseVec->data, 0, dpphiBydrReverseVec->length * sizeof(REAL8));

	*rReverseVec = XLALReverseREAL8Vector(rVec);
    *pphiReverseVec = XLALReverseREAL8Vector(pphiVec);

    *dpphiBydrReverseVec = XLALFDDerivative1Order2(rReverseVec, pphiReverseVec);

    *dpphiBydrVec = XLALReverseREAL8Vector(dpphiBydrReverseVec);
    *dpphiBydr0Vec = XLALReverseREAL8Vector(dpphiBydrReverseVec);

    XLALDestroyREAL8Vector(rReverseVec);
    XLALDestroyREAL8Vector(pphiReverseVec);
    XLALDestroyREAL8Vector(dpphiBydrReverseVec);

	return XLAL_SUCCESS;
}

REAL8
XLALSimInspiralEOBPAHamiltonianDerivative(
	REAL8 h,
	REAL8 r,
	REAL8 prstar,
	REAL8 pphi,
	SpinEOBHCoeffs *seobCoeffs,
	LALDict *LALParams
)
{
	REAL8 coeffs[9] = {1./280., -4./105., 1./5., -4./5., 0, 4./5., -1./5., 4./105., -1./280.};

	REAL8 Hderivative;
	Hderivative = 0.;

	INT4 i;

	for (i = -4; i <= 4; i++)
	{
		Hderivative += coeffs[i+4] * XLALSimInspiralEOBPAHamiltonianWrapper(r + i*h, prstar, pphi, seobCoeffs, LALParams);
	}

	Hderivative /= h;

	return Hderivative;
}

REAL8
XLALSimInspiralEOBPAHamiltonianPartialDerivativeprstar(
	REAL8 h,
	REAL8 r,
	REAL8 prstar,
	REAL8 pphi,
	SpinEOBHCoeffs *seobCoeffs,
	LALDict *LALParams
)
{
	REAL8 coeffs[9] = {1./280., -4./105., 1./5., -4./5., 0, 4./5., -1./5., 4./105., -1./280.};

	REAL8 Hderivative;
	Hderivative = 0.;

	INT4 i;

	for (i = -4; i <= 4; i++)
	{
		Hderivative += coeffs[i+4] * XLALSimInspiralEOBPAHamiltonianWrapper(r, prstar + i*h, pphi, seobCoeffs, LALParams);
	}

	Hderivative /= h;

	return Hderivative;
}

REAL8
XLALSimInspiralEOBPAHamiltonianWrapper(
	REAL8 r,
	REAL8 prstar,
	REAL8 pphi,
	SpinEOBHCoeffs *seobCoeffs,
	LALDict *LALParams
)
{
	const REAL8 nu = XLALDictLookupREAL8Value(LALParams, "nu");
	const REAL8 a1 = XLALDictLookupREAL8Value(LALParams, "a1");
	const REAL8 a2 = XLALDictLookupREAL8Value(LALParams, "a2");
	const REAL8 aK = XLALDictLookupREAL8Value(LALParams, "aK");
	const REAL8 Sstar = XLALDictLookupREAL8Value(LALParams, "Sstar");

	REAL8 H;

	REAL8Vector *xCartVec = XLALCreateREAL8Vector(3);
	memset(xCartVec->data, 0, xCartVec->length * sizeof(REAL8));

	REAL8Vector *pCartVec = XLALCreateREAL8Vector(3);
	memset(pCartVec->data, 0, pCartVec->length * sizeof(REAL8));

	REAL8Vector *a1CartVec = XLALCreateREAL8Vector(3);
	memset(a1CartVec->data, 0, a1CartVec->length * sizeof(REAL8));

    REAL8Vector *a2CartVec = XLALCreateREAL8Vector(3);
    memset(a2CartVec->data, 0, a2CartVec->length * sizeof(REAL8));

    REAL8Vector *aKCartVec = XLALCreateREAL8Vector(3);
    memset(aKCartVec->data, 0, aKCartVec->length * sizeof(REAL8));

    REAL8Vector *SstarCartVec = XLALCreateREAL8Vector(3);
    memset(SstarCartVec->data, 0, SstarCartVec->length * sizeof(REAL8));

    xCartVec->data[0] = r;

    pCartVec->data[0] = prstar;
    pCartVec->data[1] = pphi / r;

    a1CartVec->data[2] = a1;

    a2CartVec->data[2] = a2;

    aKCartVec->data[2] = aK;

    SstarCartVec->data[2] = Sstar;

    // tortoise flag:
    // 0 - unstarred coordinates pr, pphi
    // 1 - starred coordinates prstar, pphistar = pphi
    UINT4 tortoiseFlag;
    tortoiseFlag = 1;

    H = XLALSimIMRSpinEOBHamiltonian(nu, xCartVec, pCartVec, a1CartVec, a2CartVec, aKCartVec, SstarCartVec, tortoiseFlag, seobCoeffs);

    return H;
}

REAL8
XLALSimInspiralEOBPAFluxWrapper(
	REAL8 r,
	REAL8 prstar,
	REAL8 pphi,
	REAL8 omega,
	REAL8 H,
	SpinEOBParams *seobParams,
	EOBNonQCCoeffs *nqcCoeffs,
	UNUSED LALDict *LALParams
)
{
	/* polarDynamics contains r, phi, pr, pphi */
	REAL8Vector *polarDynamics = XLALCreateREAL8Vector(4);
	memset(polarDynamics->data, 0, polarDynamics->length * sizeof(REAL8));

	polarDynamics->data[0] = r;
	polarDynamics->data[1] = 0.0;
	polarDynamics->data[2] = prstar;
	polarDynamics->data[3] = pphi;

	const UINT4 lMax = 8;

	const UINT4 SpinAlignedEOBversion = 4;
	
 	REAL8 Flux;

    Flux = XLALInspiralSpinFactorizedFlux(polarDynamics, nqcCoeffs, omega, seobParams, H, lMax, SpinAlignedEOBversion);

    return Flux;
}

int
XLALSimInspiralEOBPostAdiabatic(
	UNUSED REAL8Vector **dynamics,
	/**<< OUTPUT, real part of the modes */
	UNUSED REAL8 deltaT,
	/**<< sampling time step */
	UNUSED const REAL8 m1SI,
	/**<< mass-1 in SI unit */
	UNUSED const REAL8 m2SI,
	/**<< mass-2 in SI unit */
	const REAL8 spin1z,
	/**<< z-component of spin-1, dimensionless */
	const REAL8 spin2z,
	/**<< z-component of spin-2, dimensionless */
	UINT4 SpinAlignedEOBversion,
	/**<< 1 for SEOBNRv1, 2 for SEOBNRv2, 4 for SEOBNRv4, 201 for SEOBNRv2T, 401 for SEOBNRv4T, 41 for SEOBNRv4HM */
	SpinEOBParams *seobParams,
	EOBNonQCCoeffs *nqcCoeffs
)
{
	if (SpinAlignedEOBversion != 4)
	{
		XLALPrintError("XLALSimInspiralEOBPostAdiabatic can only be used with SpinAlignedEOBversion = 4.");
		return XLAL_FAILURE;
	}

	// input variables
	UNUSED REAL8 M;
	UNUSED REAL8 q;

	UNUSED REAL8 chi1;
	UNUSED REAL8 chi2;

	UNUSED INT4 useSpins;
	UNUSED INT4 useTidal;

	UNUSED REAL8 fInit;

	UNUSED UINT4 rSize;

	UNUSED REAL8 dr;

	UNUSED UINT4 PA_order;

	// swap these for readings
	M = 1.;
	q = 1.;

	chi1 = spin1z;
	chi2 = spin2z;

	useSpins = 1;
	useTidal = 0;

	fInit = 0.00059105892307;

	const char centrifugalRadius[] = "LO";
	const char useFlm[] = "SSLO";

	// dr = 1.9774396742060579e-01;

	PA_order = 8;
 
	// set by the code
	// REAL8 m1; // unused
	// REAL8 m2; // unused
	REAL8 nu;
	REAL8 X1;
	REAL8 X2;

	REAL8 a1;
	REAL8 a2;
	REAL8 aK;
	REAL8 aK2;

	REAL8 S1;
	REAL8 S2;
	REAL8 S;
	REAL8 Sstar;

	REAL8 c3;

	REAL8 C_Q1;
	REAL8 C_Q2;

	// REAL8 time_units_factor; // unused

	// REAL8 f0; // unused
	REAL8 r0;
	REAL8 rMin;
	REAL8 z3;

	// set the above variables

	nu = XLALSimInspiralEOBPACalculateSymmetricMassRatio(q);
	X1 = XLALSimInspiralEOBPostAdiabaticX1(nu);
	X2 = XLALSimInspiralEOBPostAdiabaticX2(nu);

	a1 = X1 * chi1;
	a2 = X2 * chi2;
	aK = a1 + a2;
	aK2 = aK * aK;

	S1 = X1 * X1 * chi1;
	S2 = X2 * X2 * chi2;
	S = S1 + S2;
	Sstar = X2 * a1 + X1 * a2;

	c3 = XLALSimInspiralEOBPostAdiabaticFitGlobalc3(nu, a1, a2);

	// Self-spin coefficients
	C_Q1 = 1.;
	C_Q2 = 1.;

	// time_units_factor = XLALSimInspiralEOBPostAdiabaticTimeUnitsFactor(M);

	// f0 = fInit / time_units_factor; // unused
	r0 = XLALSimInspiralEOBPostAdiabaticDynr0Kepler(fInit); // should be (f0)
	
	// rMin = XLALSimInspiralEOBPostAdiabaticFinalRadius(q, a1, a2);
	rMin = 9.;

	z3 = XLALSimInspiralEOBPostAdiabaticz3(nu);

	rSize = 300;

	dr = (r0-rMin) / (rSize-1);

	LALDict *LALparams = XLALCreateDict();

	XLALDictInsertREAL8Value(LALparams, "M", M);
	XLALDictInsertREAL8Value(LALparams, "q", q);

	XLALDictInsertINT4Value(LALparams, "useSpins", useSpins);
	XLALDictInsertINT4Value(LALparams, "useTidal", useTidal);

	XLALDictInsertREAL8Value(LALparams, "chi1", chi1);
	XLALDictInsertREAL8Value(LALparams, "chi2", chi2);

	XLALDictInsertREAL8Value(LALparams, "rInitial", r0);
	XLALDictInsertREAL8Value(LALparams, "rFinal", rMin);
	XLALDictInsertUINT4Value(LALparams, "rSize", rSize);
	XLALDictInsertREAL8Value(LALparams, "dr", dr);

	XLALDictInsertStringValue(LALparams, "centrifugalRadius", centrifugalRadius);
	XLALDictInsertStringValue(LALparams, "useFlm", useFlm);

	XLALDictInsertREAL8Value(LALparams, "nu", nu);

	XLALDictInsertREAL8Value(LALparams, "X1", X1);
	XLALDictInsertREAL8Value(LALparams, "X2", X2);

	XLALDictInsertREAL8Value(LALparams, "a1", a1);
	XLALDictInsertREAL8Value(LALparams, "a2", a2);
	XLALDictInsertREAL8Value(LALparams, "aK", aK);
	XLALDictInsertREAL8Value(LALparams, "aK2", aK2);

	XLALDictInsertREAL8Value(LALparams, "S1", S1);
	XLALDictInsertREAL8Value(LALparams, "S2", S2);
	XLALDictInsertREAL8Value(LALparams, "S", S);
	XLALDictInsertREAL8Value(LALparams, "Sstar", Sstar);

	XLALDictInsertREAL8Value(LALparams, "c3", c3);
	XLALDictInsertREAL8Value(LALparams, "z3", z3);

	XLALDictInsertREAL8Value(LALparams, "C_Q1", C_Q1);
	XLALDictInsertREAL8Value(LALparams, "C_Q2", C_Q2);

	REAL8Vector *rVec = XLALCreateREAL8Vector(rSize);
	memset(rVec->data, 0, rVec->length * sizeof(REAL8));

	REAL8Vector *rReverseVec = XLALCreateREAL8Vector(rSize);
	memset(rReverseVec->data, 0, rReverseVec->length * sizeof(REAL8));

	REAL8Vector *tVec = XLALCreateREAL8Vector(rSize);
    memset(tVec->data, 0, tVec->length * sizeof(REAL8));

    REAL8Vector *tReverseVec = XLALCreateREAL8Vector(rSize);
    memset(tReverseVec->data, 0, tReverseVec->length * sizeof(REAL8));

    REAL8Vector *dtBydrVec = XLALCreateREAL8Vector(rSize);
    memset(dtBydrVec->data, 0, dtBydrVec->length * sizeof(REAL8));

    REAL8Vector *dtBydrReverseVec = XLALCreateREAL8Vector(rSize);
    memset(dtBydrReverseVec->data, 0, dtBydrReverseVec->length * sizeof(REAL8));

	REAL8Vector *phiVec = XLALCreateREAL8Vector(rSize);
	memset(phiVec->data, 0, phiVec->length * sizeof(REAL8));

	REAL8Vector *phiReverseVec = XLALCreateREAL8Vector(rSize);
    memset(phiReverseVec->data, 0, phiReverseVec->length * sizeof(REAL8));

    REAL8Vector *dphiBydrVec = XLALCreateREAL8Vector(rSize);
    memset(dphiBydrVec->data, 0, dphiBydrVec->length * sizeof(REAL8));

    REAL8Vector *dphiBydrReverseVec = XLALCreateREAL8Vector(rSize);
    memset(dphiBydrReverseVec->data, 0, dphiBydrReverseVec->length * sizeof(REAL8));

	REAL8Vector *HVec = XLALCreateREAL8Vector(rSize);
	memset(HVec->data, 0, HVec->length * sizeof(REAL8));

	REAL8Vector *csiVec = XLALCreateREAL8Vector(rSize);
	memset(csiVec->data, 0, csiVec->length * sizeof(REAL8));

	REAL8Vector *omegaVec = XLALCreateREAL8Vector(rSize);
	memset(omegaVec->data, 0, omegaVec->length * sizeof(REAL8));

	REAL8Vector *fluxVec = XLALCreateREAL8Vector(rSize);
	memset(fluxVec->data, 0, fluxVec->length * sizeof(REAL8));

    REAL8Vector *pphiVec = XLALCreateREAL8Vector(rSize);
    memset(pphiVec->data, 0, pphiVec->length * sizeof(REAL8));

    REAL8Vector *pphi0Vec = XLALCreateREAL8Vector(rSize);
    memset(pphi0Vec->data, 0, pphi0Vec->length * sizeof(REAL8));

    REAL8Vector *pphiReverseVec = XLALCreateREAL8Vector(rSize);
	memset(pphiReverseVec->data, 0, pphiReverseVec->length * sizeof(REAL8));

    REAL8Vector *prstarVec = XLALCreateREAL8Vector(rSize);
    memset(prstarVec->data, 0, prstarVec->length * sizeof(REAL8));

    REAL8Vector *dprstarBydrVec = XLALCreateREAL8Vector(rSize);
    memset(dprstarBydrVec->data, 0, dprstarBydrVec->length * sizeof(REAL8));

    REAL8Vector *dpphiBydrVec = XLALCreateREAL8Vector(rSize);
    memset(dpphiBydrVec->data, 0, dpphiBydrVec->length * sizeof(REAL8));

    REAL8Vector *dpphiBydr0Vec = XLALCreateREAL8Vector(rSize);
    memset(dpphiBydr0Vec->data, 0, dpphiBydr0Vec->length * sizeof(REAL8));

    REAL8Vector *dpphiBydrReverseVec = XLALCreateREAL8Vector(rSize);
	memset(dpphiBydrReverseVec->data, 0, dpphiBydrReverseVec->length * sizeof(REAL8));

    REAL8Vector *prstarReverseVec = XLALCreateREAL8Vector(rSize);
    memset(prstarReverseVec->data, 0, prstarReverseVec->length * sizeof(REAL8));

    REAL8Vector *dprstarBydrReverseVec = XLALCreateREAL8Vector(rSize);
    memset(dprstarBydrReverseVec->data, 0, dprstarBydrReverseVec->length * sizeof(REAL8));

	// REAL8 r;

	UINT4 i;

 	// compute r grid
 	XLALSimInspiralEOBPACalculateRadialGrid(
 		rVec,
 		LALparams
 	);

 	XLALSimInspiralEOBPACalculateAdiabaticDynamics(
		rVec,
		phiVec,
		prstarVec,
		pphiVec,
		pphi0Vec,
		dpphiBydrVec,
		dpphiBydr0Vec,
		HVec,
		csiVec,
		omegaVec,
		seobParams,
		LALparams
	);

	*rReverseVec = XLALReverseREAL8Vector(rVec);

    UINT4 n;
    UINT4 parity;

    UNUSED REAL8 dHBydprstar;
    REAL8 partialHBypartialprstar;

    // UINT4 count = 1;

    for (n = 1; n <= PA_order; n++)
    {
    	// Separating even and odd orders
    	if (n%2 == 0)
    	{
    		parity = 0;
    	}
    	else
    	{
    		parity = 1;
    	}

    	for (i = 0; i < rSize; i++)
    	{
    		if (parity)
    		{
    			// Odd PA orders: corrections to prstar only

                struct PostAdiabaticRootSolveParams prstarParams;
                prstarParams.r = rVec->data[i];
                prstarParams.phi = 0.;
			    prstarParams.pphi = pphiVec->data[i];
			    prstarParams.dpphiBydr = dpphiBydrVec->data[i];
			    prstarParams.omega = omegaVec->data[i];
			    prstarParams.csi = csiVec->data[i];
			    prstarParams.LALParams = LALparams;
			    prstarParams.seobParams = seobParams;
			    prstarParams.nqcCoeffs = nqcCoeffs;

    			REAL8 x_lower;
			    REAL8 x_upper;

			    if (prstarVec->data[i] > 0.)
			    {
					x_lower = 0.9 * prstarVec->data[i];
					x_upper = 1.1 * prstarVec->data[i];
				}
				else
				{
					x_upper = 0.9 * prstarVec->data[i];
					x_lower = 1.1 * prstarVec->data[i];
				}

				struct PostAdiabaticRoot prstarRoot;

			    XLALSimInspiralEOBPostAdiabaticRootFinder(
			    	&prstarRoot,
					XLALSimInspiralEOBPostAdiabaticdprstarFunc,
					&prstarParams,
					x_lower,
					x_upper,
					1.e-8,
					1.e-8
				);

    			prstarVec->data[i] = prstarRoot.root;

    			// compute Omega again
    		}
    		else
    		{
    			// Even PA orders: corrections to pphi only

    			struct PostAdiabaticRoot pphiRoot;

    			struct PostAdiabaticRootSolveParams pphiParams;
    			pphiParams.r = rVec->data[i];
			    pphiParams.csi = csiVec->data[i];
			    pphiParams.prstar = prstarVec->data[i];
			    pphiParams.dr = prstarVec->data[i];
			    pphiParams.dprstarBydr = dprstarBydrVec->data[i];

			    REAL8 x_lower = 0.9 * pphiVec->data[i];
				REAL8 x_upper = 1.1 * pphiVec->data[i];

			    XLALSimInspiralEOBPostAdiabaticRootFinder(
			    	&pphiRoot,
					XLALSimInspiralEOBPostAdiabaticdpphiFunc,
					&pphiParams,
					x_lower,
					x_upper,
					1.e-17,
					1.e-20
				);

				pphiVec->data[i] = pphiRoot.root;
    		}

    		REAL8 polarDynamics[4];
		
			polarDynamics[0] = rVec->data[i];
			polarDynamics[1] = phiVec->data[i];
			polarDynamics[2] = prstarVec->data[i];
			polarDynamics[3] = pphiVec->data[i];

			omegaVec->data[i] = XLALSimIMRSpinAlignedEOBCalcOmega(
									polarDynamics,
									seobParams,
									dr
								);

			partialHBypartialprstar = XLALSimInspiralEOBPAHamiltonianPartialDerivativeprstar(
																					1.e-5,
																					rVec->data[i],
																					prstarVec->data[i],
																					pphiVec->data[i],
																					seobParams->seobCoeffs,
																					LALparams
																				);
	        dtBydrVec->data[i] = 1. / partialHBypartialprstar;

	        dphiBydrVec->data[i] = omegaVec->data[i] * dtBydrVec->data[i];
    	}

    	if (parity)
		{
		    *prstarReverseVec = XLALReverseREAL8Vector(prstarVec);
		    *dprstarBydrReverseVec = XLALFDDerivative1Order8(rReverseVec, prstarReverseVec);
		    *dprstarBydrVec = XLALReverseREAL8Vector(dprstarBydrReverseVec);
		}
		else
		{
			*pphiReverseVec = XLALReverseREAL8Vector(pphiVec);
			*dpphiBydrReverseVec = XLALFDDerivative1Order8(rReverseVec, pphiReverseVec);
			*dpphiBydrVec = XLALReverseREAL8Vector(dpphiBydrReverseVec);
		}
    }

    *dtBydrReverseVec = XLALReverseREAL8Vector(dtBydrVec);
    *tReverseVec = XLALCumulativeIntegral3(rReverseVec, dtBydrReverseVec);
    *tVec = XLALReverseREAL8Vector(tReverseVec);

    // *tVec = XLALOffsetREAL8Vector(tVec, -tVec->data[0]);

    *dphiBydrReverseVec = XLALReverseREAL8Vector(dphiBydrVec);
    *phiReverseVec = XLALCumulativeIntegral3(rReverseVec, dphiBydrReverseVec);
    *phiVec = XLALReverseREAL8Vector(phiReverseVec);

    // *phiVec = XLALOffsetREAL8Vector(phiVec, -phiVec->data[0]);

    for (i = 0; i < rSize; i++)
    {
        printf("%.18e, %.18e, %.18e\n", rVec->data[i], tVec->data[i], phiVec->data[i]);
    }

    exit(0);

    // Output ready
    // tVec, rVec, phiVec, pphiVec, prstarVec, pphi0Vec, rcVec, AVec, dpphiBydrVec

    return XLAL_SUCCESS;
}