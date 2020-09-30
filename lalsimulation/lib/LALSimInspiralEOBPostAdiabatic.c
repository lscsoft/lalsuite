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
#include "LALSimIMREOBNewtonianMultipole.c"
#include "LALSimIMRSpinEOBFactorizedFlux_PA.c"


#define PA_AD_THRS 8.0e-3

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
	REAL8 dprstar = prstarParams->dprstar;
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
									dprstar,
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
				seobParams,
				nqcCoeffs,
				LALParams
			);

	REAL8 result;

	result = dpphiBydr*partialHBypartialprstar*dprstarBydpr - flux/omega;

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
	REAL8 dprstar = pphiParams->dprstar;
	REAL8 pphi = pphi_sol;
	REAL8 omega = pphiParams->omega;
	SpinEOBHCoeffs *seobCoeffs = pphiParams->seobParams->seobCoeffs;
	LALDict *LALParams = pphiParams->LALParams;

	const REAL8 dr = XLALDictLookupREAL8Value(LALParams, "dr");

	REAL8 partialHBypartialprstar;
	partialHBypartialprstar = XLALSimInspiralEOBPAHamiltonianPartialDerivativeprstar(
									dprstar,
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

	REAL8 flux;
    flux = XLALSimInspiralEOBPAFluxWrapper(
				r,
				prstar,
				pphi,
				omega,
				pphiParams->seobParams,
				pphiParams->nqcCoeffs,
				LALParams
			);
	
	REAL8 result;
    result = partialHBypartialprstar*dprstarBydr + partialHBypartialr - (prstar/pphi)*flux/(omega*dprstarBydpr);

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

    REAL8 partialHBypartialr;
    partialHBypartialr = XLALSimInspiralEOBPAHamiltonianDerivative(dr, r, prstar, j0_sol, seobCoeffs, LALParams);

    // REAL8 result;
    // result = partialHBypartialr;

	return partialHBypartialr;
}

int
XLALSimInspiralEOBPostAdiabaticRootFinder(
	struct PostAdiabaticRoot *result,
	double (*Func)(REAL8, void *),
	struct PostAdiabaticRootSolveParams *params,
	REAL8 x_lower,
	REAL8 x_upper,
	REAL8 absTol,
	REAL8 relTol,
	INT2 parity
)
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
    solver_type = gsl_root_fsolver_falsepos;
    // solver_type = gsl_root_fsolver_brent;

	gsl_root_fsolver *solver;
	solver = gsl_root_fsolver_alloc(solver_type);

	REAL8 F_lower = Func(x_lower, params);
	REAL8 F_upper = Func(x_upper, params);

	if (parity)
	{
		if (F_lower*F_upper >= 0.0)
		{
			x_lower = -0.5;
			x_upper = -1.e-16;
		}

		F_lower = Func(x_lower, params);
		F_upper = Func(x_upper, params);

		if (F_lower*F_upper >= 0.0)
		{
			REAL8 y;

			for (y = -1.0; y <= 1.0; y+=0.05)
			{
				printf("%.18e %.18e\n", y, Func(y, params));
			}

			printf("Derivatives have the wrong sign.\n");

			XLAL_ERROR(XLAL_EFUNC);
		}
	}
	else
	{
		while (F_lower*F_upper >= 0.0)
		{
			x_lower *= 0.9;
			x_upper *= 1.1;

			F_lower = Func(x_lower, params);
			F_upper = Func(x_upper, params);
		}
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

    if (status != GSL_SUCCESS)
    {
		printf("Root finding status: %d\n", status);
		XLAL_ERROR(XLAL_EFUNC);
	}

	gsl_root_fsolver_free (solver);
    return XLAL_SUCCESS;
}

UNUSED REAL8 XLALSimInspiralEOBPACalculateAdibaticParameter(
	REAL8 r,
	REAL8 prstar,
	REAL8 pphi,
	SpinEOBHCoeffs * seobCoeffs,
	REAL8 csi,
	LALDict *LALParams,
	REAL8 omega,
	REAL8 domegadr
)
{
	// Compute the adiabatic parameter \dot{\Omega}/2\Omega^{2}
	// where \Omega is the *orbital* frequency

	// We compute dOmega/dt = dOmega/dr * dr/dt
	
	//dr/dt = dHdpr
	REAL8 partialHBypartialprstar;
	partialHBypartialprstar = XLALSimInspiralEOBPAHamiltonianPartialDerivativeprstar(
													1e-4,
													r,
													prstar,
													pphi,
													seobCoeffs,
													LALParams
												);
	REAL8 drdt = partialHBypartialprstar*csi;
	return drdt*domegadr/(2*omega*omega);
}

int
XLALReverseREAL8Vector(
		       REAL8Vector *Vec,
		       REAL8Vector *reverseVec
)
{
	UINT4 vecLength;
	vecLength = Vec->length;

	UINT4 i;

	for (i = 0; i < vecLength; i++)
	{
		reverseVec->data[i] = Vec->data[vecLength-i-1];
	}

	return XLAL_SUCCESS;
}

int
XLALOffsetREAL8Vector(
	REAL8Vector *Vec,
	REAL8 offset,
	REAL8Vector *offsetVec
)
{
	UINT4 vecLength;
	vecLength = Vec->length;

	UINT4 i;

	for (i = 0; i < vecLength; i++)
	{
		offsetVec->data[i] = Vec->data[i] + offset;
	}

	return XLAL_SUCCESS;
}

int
XLALRescaleREAL8Vector(
	REAL8Vector *Vec,
	REAL8 factor,
	REAL8Vector *rescaledVec
)
{
	UINT4 vecLength;
	vecLength = Vec->length;

	UINT4 i;

	for (i = 0; i < vecLength; i++)
	{
		rescaledVec->data[i] = factor * Vec->data[i];
	}

	return XLAL_SUCCESS;
}

int
XLALSimInspiralEOBPACalculateAdiabaticDynamics(
	REAL8Vector *rVec,
	REAL8Vector *phiVec,
	REAL8Vector *prstarVec,
	REAL8Vector *pphiVec,
	REAL8Vector *pphi0Vec,
	REAL8Vector *dpphiBydrVec,
	UNUSED REAL8Vector *dpphiBydr0Vec,
	REAL8Vector *csiVec,
	REAL8Vector *omegaVec,
	SpinEOBParams *seobParams,
	UNUSED EOBNonQCCoeffs *nqcCoeffs,
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
	struct PostAdiabaticRootSolveParams pphiParams;
	for (i = 0; i < rSize; i++)
	{
		REAL8 Newtonianj0;
		Newtonianj0 = XLALSimInspiralEOBPACalculateNewtonianj0(rVec->data[i]);

		

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
			1.e-14,
			1.e-16,
			0
		);

		pphiVec->data[i] = pphiRoot.root;
        pphi0Vec->data[i] = pphiRoot.root;

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

	XLALReverseREAL8Vector(rVec,rReverseVec);
	XLALReverseREAL8Vector(pphiVec,pphiReverseVec);
	XLALFDDerivative1Order8(rReverseVec, pphiReverseVec, dpphiBydrReverseVec);
	XLALReverseREAL8Vector(dpphiBydrReverseVec, dpphiBydrVec);

    XLALDestroyREAL8Vector(rReverseVec);
    XLALDestroyREAL8Vector(pphiReverseVec);
    XLALDestroyREAL8Vector(dpphiBydrReverseVec);

	return XLAL_SUCCESS;
}

int
XLALSimInspiralEOBPACalculatePostAdiabaticDynamics(
	REAL8Vector *rVec,
	REAL8Vector *phiVec,
	REAL8Vector *dphiBydrVec,
	REAL8Vector *prstarVec,
	REAL8Vector *dprstarBydrVec,
	REAL8Vector *pphiVec,
	REAL8Vector *dpphiBydrVec,
	REAL8Vector *dtBydrVec,
	REAL8Vector *csiVec,
	REAL8Vector *omegaVec,
	SpinEOBParams *seobParams,
	EOBNonQCCoeffs *nqcCoeffs,
	LALDict *LALParams
)
{
	const UINT4 PAOrder = XLALDictLookupUINT4Value(LALParams, "PAOrder");
	const UINT4 rSize = XLALDictLookupUINT4Value(LALParams, "rSize");
	const REAL8 dr = XLALDictLookupREAL8Value(LALParams, "dr");

	UINT4 n;
	UINT4 parity;
	UINT4 i;

	REAL8Vector *fluxVec = XLALCreateREAL8Vector(rSize);
	memset(fluxVec->data, 0, fluxVec->length * sizeof(REAL8));

	REAL8Vector *rReverseVec = XLALCreateREAL8Vector(rSize);
	memset(rReverseVec->data, 0, rReverseVec->length * sizeof(REAL8));
	XLALReverseREAL8Vector(rVec,rReverseVec);

	REAL8Vector *pphiReverseVec = XLALCreateREAL8Vector(rSize);
	memset(pphiReverseVec->data, 0, pphiReverseVec->length * sizeof(REAL8));

	REAL8Vector *dpphiBydrReverseVec = XLALCreateREAL8Vector(rSize);
	memset(dpphiBydrReverseVec->data, 0, dpphiBydrReverseVec->length * sizeof(REAL8));

	REAL8Vector *prstarReverseVec = XLALCreateREAL8Vector(rSize);
	memset(prstarReverseVec->data, 0, prstarReverseVec->length * sizeof(REAL8));

	REAL8Vector *dprstarBydrReverseVec = XLALCreateREAL8Vector(rSize);
	memset(dprstarBydrReverseVec->data, 0, dprstarBydrReverseVec->length * sizeof(REAL8));

	// Test
	REAL8Vector *dprstarVec = XLALCreateREAL8Vector(rSize);
	memset(dprstarVec->data, 0, dprstarVec->length * sizeof(REAL8));

	for (i = 0; i < rSize; i++)
	{
		dprstarVec->data[i] = 1.e-4;
	}
	// Test ends

	REAL8 polarDynamics[4];

	for (n = 1; n <= PAOrder; n++)
	{
		// printf("Calculating %d PA order\n", n);

		parity = n%2;

		if (n > 1)
		{
			XLALSimInspiralEOBPAMeanValueOrder8(prstarVec, dprstarVec);

			// for (i = 0; i < rSize; i++)
	  //   	{
	  //   		printf("dprstar = %.18e\n", dprstarVec->data[i]);
	  //   	}

	  //   	exit(0);
		}

		if (parity)
		{
			for (i = 0; i < rSize; i++)
	    	{
                struct PostAdiabaticRootSolveParams prstarParams;
                prstarParams.r = rVec->data[i];
                prstarParams.dprstar = dprstarVec->data[i];
                prstarParams.phi = 0.;
			    prstarParams.pphi = pphiVec->data[i];
			    prstarParams.dpphiBydr = dpphiBydrVec->data[i];
			    prstarParams.omega = omegaVec->data[i];
			    prstarParams.csi = csiVec->data[i];
			    prstarParams.LALParams = LALParams;
			    prstarParams.seobParams = seobParams;
			    prstarParams.nqcCoeffs = nqcCoeffs;
    			REAL8 x_lower;
			    REAL8 x_upper;

				x_upper = 0.8 * prstarVec->data[i];
				x_lower = 1.2 * prstarVec->data[i];

				struct PostAdiabaticRoot prstarRoot;

			    XLALSimInspiralEOBPostAdiabaticRootFinder(
			    	&prstarRoot,
					XLALSimInspiralEOBPostAdiabaticdprstarFunc,
					&prstarParams,
					x_lower,
					x_upper,
					1.e-14,
					1.e-16,
					parity
				);

    			prstarVec->data[i] = prstarRoot.root;

    			// printf("%d %.18e %.18e\n", i, prstarVec->data[i], dpphiBydrVec->data[i]);

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

			XLALReverseREAL8Vector(prstarVec,prstarReverseVec);
			memset(dprstarBydrReverseVec->data, 0, dprstarBydrReverseVec->length * sizeof(REAL8)); // Otherwise we are in trouble
			XLALFDDerivative1Order8(rReverseVec, prstarReverseVec,dprstarBydrReverseVec);
			XLALReverseREAL8Vector(dprstarBydrReverseVec,dprstarBydrVec);
		}
		else
		{
			for (i = 0; i < rSize; i++)
	    	{
                struct PostAdiabaticRoot pphiRoot;

    			struct PostAdiabaticRootSolveParams pphiParams;
    			pphiParams.r = rVec->data[i];
			    pphiParams.csi = csiVec->data[i];
			    pphiParams.prstar = prstarVec->data[i];
			    pphiParams.dprstar = dprstarVec->data[i];
			    pphiParams.omega = omegaVec->data[i];
			    pphiParams.dprstarBydr = dprstarBydrVec->data[i];

			    REAL8 x_lower = 0.9 * pphiVec->data[i];
				REAL8 x_upper = 1.1 * pphiVec->data[i];

			    XLALSimInspiralEOBPostAdiabaticRootFinder(
			    	&pphiRoot,
					XLALSimInspiralEOBPostAdiabaticdpphiFunc,
					&pphiParams,
					x_lower,
					x_upper,
					1.e-14,
					1.e-16,
					parity
				);

				pphiVec->data[i] = pphiRoot.root;
		
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

			XLALReverseREAL8Vector(pphiVec,pphiReverseVec);
			memset(dpphiBydrReverseVec->data, 0, dpphiBydrReverseVec->length * sizeof(REAL8)); // Otherwise we are in trouble
			XLALFDDerivative1Order8(rReverseVec, pphiReverseVec,dpphiBydrReverseVec);
			XLALReverseREAL8Vector(dpphiBydrReverseVec,dpphiBydrVec);
		}
	}

	// Loop over and compute the derivatives we need to do quadratures later
	REAL8 partialHBypartialprstar;

	for(i = 0; i < rSize; i++)
	{
		partialHBypartialprstar = XLALSimInspiralEOBPAHamiltonianPartialDerivativeprstar(
													1.e-4,
													rVec->data[i],
													prstarVec->data[i],
													pphiVec->data[i],
													seobParams->seobCoeffs,
													LALParams
												);
		dtBydrVec->data[i] = 1. / (partialHBypartialprstar*csiVec->data[i]);
	    dphiBydrVec->data[i] = omegaVec->data[i] * dtBydrVec->data[i];
	}

	XLALDestroyREAL8Vector(fluxVec);
	XLALDestroyREAL8Vector(rReverseVec);
    XLALDestroyREAL8Vector(pphiReverseVec);
    XLALDestroyREAL8Vector(dpphiBydrReverseVec);
    XLALDestroyREAL8Vector(prstarReverseVec);
    XLALDestroyREAL8Vector(dprstarBydrReverseVec);
    XLALDestroyREAL8Vector(dprstarVec);

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
		if (i != 0)
		{
			Hderivative += coeffs[i+4] * XLALSimInspiralEOBPAHamiltonianWrapper(r + i*h, prstar, pphi, seobCoeffs, LALParams);
		}
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
		if (i!=0)
		{
			Hderivative += coeffs[i+4] * XLALSimInspiralEOBPAHamiltonianWrapper(r, prstar + i*h, pphi, seobCoeffs, LALParams);
		}
	}

	Hderivative /= h;

	return Hderivative;
}

int
XLALSimInspiralEOBPAMeanValueOrder8(
	REAL8Vector *inputVec,
	REAL8Vector *meanVec
)
{
	UINT4 vecLength = inputVec->length;

	UINT4 i;
	UINT4 j;

	for (i = 0; i < vecLength; i++)
	{
		if (i == 0)
		{
			for (j = 0; j < 8; j++)
			{
				meanVec->data[i] += fabs(inputVec->data[j+1] - inputVec->data[j]);
			}
		}
		else if (i == 1)
		{
			for (j = 0; j < 8; j++)
			{
				meanVec->data[i] += fabs(inputVec->data[j+1] - inputVec->data[j]);
			}
		}
		else if (i == 2)
		{
			for (j = 0; j < 8; j++)
			{
				meanVec->data[i] += fabs(inputVec->data[j+1] - inputVec->data[j]);
			}
		}
		else if (i == 3)
		{
			for (j = 0; j < 8; j++)
			{
				meanVec->data[i] += fabs(inputVec->data[j+1] - inputVec->data[j]);
			}
		}
		else if (i == vecLength-4)
		{
			for (j = 0; j < 8; j++)
			{
				meanVec->data[i] += fabs(inputVec->data[i+j-5] - inputVec->data[i+j-4]);
			}
		}
		else if (i == vecLength-3)
		{
			for (j = 0; j < 8; j++)
			{
				meanVec->data[i] += fabs(inputVec->data[i+j-6] - inputVec->data[i+j-5]);
			}
		}
		else if (i == vecLength-2)
		{
			for (j = 0; j < 8; j++)
			{
				meanVec->data[i] += fabs(inputVec->data[i+j-7] - inputVec->data[i+j-6]);
			}
		}
		else if (i == vecLength-1)
		{
			for (j = 0; j < 8; j++)
			{
				meanVec->data[i] += fabs(inputVec->data[i+j-8] - inputVec->data[i+j-7]);
			}
		}
		else
		{
			for (j = 0; j < 8; j++)
			{
				meanVec->data[i] += fabs(inputVec->data[i+j-4] - inputVec->data[i+j-3]);
			}
		}

		meanVec->data[i] /= 8;
	}

	return XLAL_SUCCESS;
}

REAL8Vector
XLALSimInspiralEOBPAHamiltonianPartialDerivativeprstarBetter(
	REAL8Vector *rVec,
    REAL8Vector *prstarVec,
    REAL8Vector *pphiVec,
    SpinEOBHCoeffs *seobCoeffs,
    LALDict *LALParams
)
{
	REAL8 eightOrderCoeffs[9][9] = {
		{-761./280., 8., -14., 56./3., -35./2., 56./5., -14./3., 8./7., -1./8.},
		{-1./8., -223./140., 7./2., -7./2., 35./12., -7./4., 7./10., -1./6., 1./56.},
		{1./56., -2./7., -19./20., 2., -5./4., 2./3., -1./4., 2./35., -1./168.},
		{-1./168., 1./14., -1./2., -9./20., 5./4., -1./2., 1./6., -1./28., 1./280.},
		{1./280., -4./105., 1./5., -4./5., 0, 4./5., -1./5., 4./105., -1./280.},
		{-1./280., 1./28., -1./6., 1./2., -5./4., 9./20., 1./2., -1./14., 1./168.},
		{1./168., -2./35., 1./4., -2./3., 5./4., -2., 19./20., 2./7., -1./56.},
		{-1./56., 1./6., -7./10., 7./4., -35./12., 7./2., -7./2., 223./140., 1./8.},
		{1./8., -8./7., 14./3., -56./5., 35./2., -56./3., 14., -8., 761./280.}
	};

	const UINT4 rSize = XLALDictLookupUINT4Value(LALParams, "rSize");

	REAL8Vector *meanprstarVec = XLALCreateREAL8Vector(rSize);
	memset(meanprstarVec->data, 0, meanprstarVec->length * sizeof(REAL8));
	XLALSimInspiralEOBPAMeanValueOrder8(prstarVec, meanprstarVec);

	REAL8Vector *partialHByPartialprstarVec = XLALCreateREAL8Vector(rSize);
	memset(partialHByPartialprstarVec->data, 0, partialHByPartialprstarVec->length * sizeof(REAL8));

	UINT4 i;
	INT4 j;

	REAL8 instantaneousH;

	for (i = 0; i < rSize; i++)
	{
		if (i == 0)
		{
			for (j = 0; j <= 8; j++)
			{
				instantaneousH = XLALSimInspiralEOBPAHamiltonianWrapper(
										rVec->data[i],
										prstarVec->data[i] + j*meanprstarVec->data[i],
										pphiVec->data[i],
										seobCoeffs,
										LALParams
									);

				partialHByPartialprstarVec->data[i] += eightOrderCoeffs[0][j] * instantaneousH;
			}
		}
		else if (i == 1)
		{
			for (j = 0; j <= 8; j++)
			{
				instantaneousH = XLALSimInspiralEOBPAHamiltonianWrapper(
										rVec->data[i],
										prstarVec->data[i] + (j-1)*meanprstarVec->data[i],
										pphiVec->data[i],
										seobCoeffs,
										LALParams
									);

				partialHByPartialprstarVec->data[i] += eightOrderCoeffs[1][j] * instantaneousH;
			}
		}
		else if (i == 2)
		{
			for (j = 0; j <= 8; j++)
			{
				instantaneousH = XLALSimInspiralEOBPAHamiltonianWrapper(
										rVec->data[i],
										prstarVec->data[i] + (j-2)*meanprstarVec->data[i],
										pphiVec->data[i],
										seobCoeffs,
										LALParams
									);

				partialHByPartialprstarVec->data[i] += eightOrderCoeffs[2][j] * instantaneousH;
			}
		}
		else if (i == 3)
		{
			for (j = 0; j <= 8; j++)
			{
				instantaneousH = XLALSimInspiralEOBPAHamiltonianWrapper(
										rVec->data[i],
										prstarVec->data[i] + (j-3)*meanprstarVec->data[i],
										pphiVec->data[i],
										seobCoeffs,
										LALParams
									);

				partialHByPartialprstarVec->data[i] += eightOrderCoeffs[3][j] * instantaneousH;
			}
		}
		else if (i == rSize-4)
		{
			for (j = 0; j <= 8; j++)
			{
				instantaneousH = XLALSimInspiralEOBPAHamiltonianWrapper(
										rVec->data[i],
										prstarVec->data[i] + (j-5)*meanprstarVec->data[i],
										pphiVec->data[i],
										seobCoeffs,
										LALParams
									);

				partialHByPartialprstarVec->data[i] += eightOrderCoeffs[5][j] * instantaneousH;
			}
		}
		else if (i == rSize-3)
		{
			for (j = 0; j <= 8; j++)
			{
				instantaneousH = XLALSimInspiralEOBPAHamiltonianWrapper(
										rVec->data[i],
										prstarVec->data[i] + (j-6)*meanprstarVec->data[i],
										pphiVec->data[i],
										seobCoeffs,
										LALParams
									);

				partialHByPartialprstarVec->data[i] += eightOrderCoeffs[6][j] * instantaneousH;
			}
		}
		else if (i == rSize-2)
		{
			for (j = 0; j <= 8; j++)
			{
				instantaneousH = XLALSimInspiralEOBPAHamiltonianWrapper(
										rVec->data[i],
										prstarVec->data[i] + (j-7)*meanprstarVec->data[i],
										pphiVec->data[i],
										seobCoeffs,
										LALParams
									);

				partialHByPartialprstarVec->data[i] += eightOrderCoeffs[7][j] * instantaneousH;
			}
		}
		else if (i == rSize-1)
		{
			for (j = 0; j <= 8; j++)
			{
				instantaneousH = XLALSimInspiralEOBPAHamiltonianWrapper(
										rVec->data[i],
										prstarVec->data[i] + (j-8)*meanprstarVec->data[i],
										pphiVec->data[i],
										seobCoeffs,
										LALParams
									);

				partialHByPartialprstarVec->data[i] += eightOrderCoeffs[8][j] * instantaneousH;
			}
		}
		else
		{
			for (j = 0; j <= 8; j++)
			{
				instantaneousH = XLALSimInspiralEOBPAHamiltonianWrapper(
										rVec->data[i],
										prstarVec->data[i] + (j-4)*meanprstarVec->data[i],
										pphiVec->data[i],
										seobCoeffs,
										LALParams
									);

				partialHByPartialprstarVec->data[i] += eightOrderCoeffs[4][j] * instantaneousH;
			}
		}

		partialHByPartialprstarVec->data[i] /= meanprstarVec->data[i];
	}

	return *partialHByPartialprstarVec;
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

	REAL8Vector xCartVec,pCartVec;
	xCartVec.length = pCartVec.length = 3;
	REAL8 xCart[3] = {r, 0., 0.};
	REAL8 pCart[3] = {prstar, pphi/r, 0.};
	xCartVec.data = xCart;
	pCartVec.data = pCart;

	REAL8Vector a1CartVec, a2CartVec, aKCartVec, SstarCartVec;
	a1CartVec.length = a2CartVec.length = aKCartVec.length = SstarCartVec.length = 3;
	REAL8 spin1[3] = {0., 0., a1};
	REAL8 spin2[3] = {0., 0., a2};
	REAL8 aKV[3] = {0., 0., aK};
	REAL8 SstarV[3] = {0., 0., Sstar};
	a1CartVec.data = spin1;
	a2CartVec.data = spin2;
	aKCartVec.data = aKV;
	SstarCartVec.data = SstarV;
    
    // tortoise flag:
    // 0 - unstarred coordinates pr, pphi
    // 1 - starred coordinates prstar, pphistar = pphi
    UINT4 tortoiseFlag;
    tortoiseFlag = 1;

    H = XLALSimIMRSpinEOBHamiltonian(
    		nu,
    		&xCartVec,
    		&pCartVec,
    		&a1CartVec,
    		&a2CartVec,
    		&aKCartVec,
    		&SstarCartVec,
    		tortoiseFlag,
    		seobCoeffs
    	);

    if (H != H)
    {
    	printf("r = %.18e\n", r);
    	printf("prstar = %.18e\n", prstar);
    	printf("pphi = %.18e\n", pphi);
    	printf("pphi/r = %.18e\n", pphi/r);
    	printf("a1 = %.18e\n", a1);
    	printf("a2 = %.18e\n", a2);
    	printf("aK = %.18e\n", aK);
    	printf("Sstar = %.18e\n\n", Sstar);
    }

	H /= nu;

    return H;
}

REAL8
XLALSimInspiralEOBPAFluxWrapper(
	REAL8 r,
	REAL8 prstar,
	REAL8 pphi,
	REAL8 omega,
	SpinEOBParams *seobParams,
	EOBNonQCCoeffs *nqcCoeffs,
	LALDict *LALParams
)
{
	const REAL8 nu = XLALDictLookupREAL8Value(LALParams, "nu");

	REAL8 H;
	H = XLALSimInspiralEOBPAHamiltonianWrapper(
			r,
			prstar,
			pphi,
			seobParams->seobCoeffs,
			LALParams
		);
	H *= nu;

	/* polarDynamics contains r, phi, pr, pphi */
	REAL8Vector polarDynamics;
	polarDynamics.length = 4;
	REAL8 pol[4] = {r,0.,prstar,pphi};
	polarDynamics.data = pol;
	
	const UINT4 lMax = 8;

	const UINT4 SpinAlignedEOBversion = 4;
	
 	REAL8 Flux;

    Flux = XLALInspiralSpinFactorizedFlux_PA(
    			&polarDynamics,
    			nqcCoeffs,
    			omega,
    			seobParams,
    			H,
    			lMax,
    			SpinAlignedEOBversion
    		);

    Flux /= nu;
    Flux *= -1.;

    return Flux;
}

int
XLALSimInspiralEOBPostAdiabatic(
	REAL8Array **dynamics,
	/**<< OUTPUT, real part of the modes */
	const REAL8 m1,
	/**<< mass-1 */
	const REAL8 m2,
	/**<< mass-2 */
	const REAL8 spin1z,
	/**<< z-component of spin-1, dimensionless */
	const REAL8 spin2z,
	/**<< z-component of spin-2, dimensionless */
	const REAL8Vector initVals,
	UINT4 SpinAlignedEOBversion,
	/**<< 1 for SEOBNRv1, 2 for SEOBNRv2, 4 for SEOBNRv4, 201 for SEOBNRv2T, 401 for SEOBNRv4T, 41 for SEOBNRv4HM */
	SpinEOBParams *seobParams,
	EOBNonQCCoeffs *nqcCoeffs,
	UNUSED LALDict *PAParams
)
{
	if (SpinAlignedEOBversion != 4)
	{
		XLALPrintError("XLALSimInspiralEOBPostAdiabatic can only be used with SpinAlignedEOBversion = 4.");
		return XLAL_FAILURE;
	}

	REAL8 M;
	M = m1 + m2;

	REAL8 q;
	q = XLALSimInspiralEOBPACalculateMassRatio(m1, m2);

	REAL8 nu;
	nu = XLALSimInspiralEOBPACalculateSymmetricMassRatio(q);

	REAL8 chi1;
	chi1 = spin1z;

	REAL8 chi2;
	chi2 = spin2z;

	const UINT4 PAOrder = XLALDictLookupUINT4Value(PAParams, "PAOrder");

	REAL8 X1;
	X1 = XLALSimInspiralEOBPACalculateX1(nu);

	REAL8 X2;
	X2 = XLALSimInspiralEOBPACalculateX2(nu);

	REAL8 a1;
	a1 = XLALSimInspiralEOBPACalculatea(X1, chi1);

	REAL8 a2;
	a2 = XLALSimInspiralEOBPACalculatea(X2, chi2);

	REAL8 aK;
	aK = a1 + a2;

	REAL8 S1;
	S1 = XLALSimInspiralEOBPACalculateS(X1, chi1);

	REAL8 S2;
	S2 = XLALSimInspiralEOBPACalculateS(X2, chi2);

	REAL8 S;
	S = S1 + S2;

	REAL8 Sstar;
	Sstar = XLALSimInspiralEOBPACalculateSstar(X1, X2, chi1, chi2);

	REAL8 rInitial;
	rInitial = initVals.data[0];

	REAL8 rFinal;
	REAL8 rFinalPrefactor = XLALDictLookupREAL8Value(PAParams, "rFinal");
	rFinal = rFinalPrefactor * XLALSimInspiralEOBPostAdiabaticFinalRadius(q, a1, a2);

	if (rInitial <= rFinal)
	{
		XLAL_ERROR(XLAL_EFUNC);
	}

	UINT4 rSize;
	rSize = XLALDictLookupUINT4Value(PAParams, "rSize");
	
	REAL8 dr;
	dr = XLALSimInspiralEOBPACalculatedr(rInitial, rFinal, rSize);

	LALDict *LALparams = XLALCreateDict();

	XLALDictInsertREAL8Value(LALparams, "M", M);
	XLALDictInsertREAL8Value(LALparams, "q", q);
	XLALDictInsertREAL8Value(LALparams, "nu", nu);

	XLALDictInsertREAL8Value(LALparams, "chi1", chi1);
	XLALDictInsertREAL8Value(LALparams, "chi2", chi2);

	XLALDictInsertUINT4Value(LALparams, "PAOrder", PAOrder);

	XLALDictInsertREAL8Value(LALparams, "X1", X1);
	XLALDictInsertREAL8Value(LALparams, "X2", X2);

	XLALDictInsertREAL8Value(LALparams, "a1", a1);
	XLALDictInsertREAL8Value(LALparams, "a2", a2);
	XLALDictInsertREAL8Value(LALparams, "aK", aK);

	XLALDictInsertREAL8Value(LALparams, "S1", S1);
	XLALDictInsertREAL8Value(LALparams, "S2", S2);
	XLALDictInsertREAL8Value(LALparams, "S", S);
	XLALDictInsertREAL8Value(LALparams, "Sstar", Sstar);

	XLALDictInsertREAL8Value(LALparams, "rInitial", rInitial);
	XLALDictInsertREAL8Value(LALparams, "rFinal", rFinal);
	XLALDictInsertUINT4Value(LALparams, "rSize", rSize);
	XLALDictInsertREAL8Value(LALparams, "dr", dr);

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

	REAL8Vector *csiVec = XLALCreateREAL8Vector(rSize);
	memset(csiVec->data, 0, csiVec->length * sizeof(REAL8));

	REAL8Vector *omegaVec = XLALCreateREAL8Vector(rSize);
	memset(omegaVec->data, 0, omegaVec->length * sizeof(REAL8));

	REAL8Vector *omegaReverseVec = XLALCreateREAL8Vector(rSize);
	memset(omegaReverseVec->data, 0, omegaReverseVec->length * sizeof(REAL8));

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

	REAL8Vector *domegadrVec = XLALCreateREAL8Vector(rSize);
	memset(domegadrVec->data, 0, domegadrVec->length * sizeof(REAL8));

	REAL8Vector *domegadrReverseVec = XLALCreateREAL8Vector(rSize);
	memset(domegadrReverseVec->data, 0, domegadrReverseVec->length * sizeof(REAL8));

	REAL8Vector *adiabatic_param_Vec = XLALCreateREAL8Vector(rSize);
	memset(adiabatic_param_Vec->data, 0, adiabatic_param_Vec->length * sizeof(REAL8));

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
		csiVec,
		omegaVec,
		seobParams,
		nqcCoeffs,
		LALparams
	);

 	XLALReverseREAL8Vector(rVec,rReverseVec);

	if (PAOrder > 0)
	{
		XLALSimInspiralEOBPACalculatePostAdiabaticDynamics(
			rVec,
			phiVec,
			dphiBydrVec,
			prstarVec,
			dprstarBydrVec,
			pphiVec,
			dpphiBydrVec,
			dtBydrVec,
			csiVec,
			omegaVec,
			seobParams,
			nqcCoeffs,
			LALparams
		);
	}

	XLALCumulativeIntegral3(rVec, dtBydrVec, tVec);

	XLALCumulativeIntegral3(rVec, dphiBydrVec, phiVec);

    UINT4 i,j;
	XLALReverseREAL8Vector(omegaVec, omegaReverseVec);
	XLALFDDerivative1Order8(rReverseVec, omegaReverseVec,domegadrReverseVec);
	XLALReverseREAL8Vector(domegadrReverseVec, domegadrVec);
    FILE *out = fopen ("pa_dyn.dat", "w");
	// Figure out where we are going to stop
	UNUSED REAL8 adiabatic_param = 0.0;
	UNUSED REAL8 r,prstar,pphi,csi,omega,domegadr;

	UINT4 idx_stop = 0;

	// for (j=0; j<rSize; j++)
	// {
	// 	r = rVec->data[j];
	// 	prstar = prstarVec->data[j];
	// 	pphi = pphiVec->data[j];
	// 	csi = csiVec->data[j];
	// 	omega = omegaVec->data[j];
	// 	domegadr = domegadrVec->data[j];
	// 	adiabatic_param = XLALSimInspiralEOBPACalculateAdibaticParameter(r,prstar,pphi,seobParams->seobCoeffs,csi,LALparams, omega, domegadr);
	// 	adiabatic_param_Vec->data[j]=adiabatic_param;
	// 	//printf("r=%.17f, omega=%.17f,domegadr=%.17f, adibatic_param = %.17f\n",r,omega,domegadr,adiabatic_param_Vec->data[j]);
	// 	if(idx_stop==0 && adiabatic_param>PA_AD_THRS){
	// 		printf("r=%.18e Q = %.18e\n",r,adiabatic_param);
	// 		idx_stop = j;
	// 	}
	// }

	REAL8 rSwitch;
	REAL8 rSwitchPrefactor = XLALDictLookupREAL8Value(PAParams, "rSwitch");
	rSwitch = rSwitchPrefactor * XLALSimInspiralEOBPostAdiabaticFinalRadius(q, a1, a2);

	if (rInitial <= rSwitch)
	{
		XLAL_ERROR(XLAL_EFUNC);
	}

	for (j=0; j<rSize; j++)
	{
		r = rVec->data[j];

		if (idx_stop == 0 && rSwitch > r)
		{
			idx_stop = j-1;
			break;
		}
	}

	if (idx_stop == 0)
	{
		idx_stop = rSize - 1;
	}

	for (i = 0; i < rSize; i++)
    {
        fprintf(
        	out,
        	"%.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e\n",
        	tVec->data[i],
        	rVec->data[i],
            phiVec->data[i],
            prstarVec->data[i],
            pphiVec->data[i],
            dtBydrVec->data[i],
            domegadrVec->data[i],
            adiabatic_param_Vec->data[i]
        );
    }

    fclose (out);

  	UINT4 outSize = idx_stop;

	*dynamics = XLALCreateREAL8ArrayL(2, 5, outSize);

	for (i = 0; i < outSize; i++)
    {
    	(*dynamics)->data[i] = tVec->data[i];
    	(*dynamics)->data[outSize + i] = rVec->data[i];
    	(*dynamics)->data[2*outSize + i] = phiVec->data[i];
    	(*dynamics)->data[3*outSize + i] = prstarVec->data[i];
    	(*dynamics)->data[4*outSize + i] = pphiVec->data[i];
    } 

	XLALDestroyREAL8Vector(tVec);
	XLALDestroyREAL8Vector(tReverseVec);
	XLALDestroyREAL8Vector(rVec);
	XLALDestroyREAL8Vector(rReverseVec);
	XLALDestroyREAL8Vector(phiVec);
	XLALDestroyREAL8Vector(phiReverseVec);
	XLALDestroyREAL8Vector(prstarVec);
	XLALDestroyREAL8Vector(pphiVec);
	XLALDestroyREAL8Vector(pphiReverseVec);
	XLALDestroyREAL8Vector(pphi0Vec);
	XLALDestroyREAL8Vector(omegaVec);
	XLALDestroyREAL8Vector(fluxVec);
	XLALDestroyREAL8Vector(csiVec);
	XLALDestroyREAL8Vector(prstarReverseVec);
	XLALDestroyREAL8Vector(dprstarBydrVec );
	XLALDestroyREAL8Vector(dpphiBydrVec);
	XLALDestroyREAL8Vector(dpphiBydr0Vec);
	XLALDestroyREAL8Vector(dpphiBydrReverseVec);
	XLALDestroyREAL8Vector(dprstarBydrReverseVec);
	XLALDestroyREAL8Vector(dphiBydrVec);
	XLALDestroyREAL8Vector(dphiBydrReverseVec);
	XLALDestroyREAL8Vector(dtBydrVec);
	XLALDestroyREAL8Vector(dtBydrReverseVec);
	XLALDestroyREAL8Vector(omegaReverseVec);
	XLALDestroyREAL8Vector(domegadrVec);
	XLALDestroyREAL8Vector(domegadrReverseVec);
	XLALDestroyREAL8Vector(adiabatic_param_Vec);

	XLALDestroyDict(LALparams);

	// XLALDestroyREAL8Array(outputDynamics);
	
    return XLAL_SUCCESS;
}
