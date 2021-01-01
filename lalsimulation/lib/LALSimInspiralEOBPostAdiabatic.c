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

#include "LALSimIMREOBNewtonianMultipole.c"

#include "LALSimIMRSpinEOBHamiltonian.c"
#include "LALSimIMRSpinEOBHamiltonianOptimized.c"

#include "LALSimIMRSpinEOBFactorizedFlux_PA.c"
#include "LALSimIMRSpinEOBFactorizedFluxOptimized.c"


#define PA_AD_THRS 8.0e-3

#define ROOT_SOLVER_ABS_TOL 1.0e-10
#define ROOT_SOLVER_REL_TOL 1.0e-8

int
XLALSimInspiralEOBPACalculateRadialGrid(
	REAL8Vector *rVec,
	LALDict *LALParams
)
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
XLALSimInspiralEOBPostAdiabaticj0Func(
	REAL8 j0_sol,
	void *params
)
{
	struct PostAdiabaticRootSolveParams *j0Params = (struct PostAdiabaticRootSolveParams *)params;

	const REAL8 dr = XLALDictLookupREAL8Value(j0Params->LALParams, "dr");

	REAL8 r = j0Params->r;
   	REAL8 prstar = j0Params->prstar;
    LALDict *LALParams = j0Params->LALParams;

    // RHS of eq. (2) of the document
    REAL8 partialHBypartialr = XLALSimInspiralEOBPAPartialHByPartialr(
		dr,
		r,
		prstar,
		j0_sol,
		j0Params->seobParams,
		LALParams
	);

	return partialHBypartialr;
}

double
XLALSimInspiralEOBPostAdiabaticdprstarFunc(
	REAL8 prstar_sol,
	void *params
)
{
	struct PostAdiabaticRootSolveParams *prstarParams = (struct PostAdiabaticRootSolveParams *)params;

	REAL8 r = prstarParams->r;
	REAL8 prstar = prstar_sol;
	REAL8 dprstar = prstarParams->dprstar;
	REAL8 pphi = prstarParams->pphi;
	REAL8 dprstarBydpr = prstarParams->csi;
	SpinEOBParams *seobParams = prstarParams->seobParams;
	EOBNonQCCoeffs *nqcCoeffs = prstarParams->nqcCoeffs;

	REAL8 dpphiBydr = prstarParams->dpphiBydr;

	LALDict *LALParams = prstarParams->LALParams;

	const UINT2 analyticFlag = XLALDictLookupUINT2Value(LALParams, "analyticFlag");

	REAL8 pol[4] = {r, 0., prstar, pphi};
	REAL8 omega;
	REAL8 flux;

	REAL8 result;

	if (analyticFlag == 0)
	{
		REAL8 partialHByPartialprstar = XLALSimInspiralEOBPAPartialHByPartialprstar(
			dprstar,
			r,
			prstar,
			pphi,
			seobParams,
			LALParams
		);

		omega = XLALSimIMRSpinAlignedEOBPACalculateOmega(
			pol,
			dprstar,
			seobParams,
			LALParams
		);

	    flux = XLALSimInspiralEOBPAFluxWrapper(
			r,
			prstar,
			pphi,
			omega,
			seobParams,
			nqcCoeffs,
			LALParams
		);

		result = dpphiBydr*partialHByPartialprstar*dprstarBydpr - flux/omega;
	}
	else
	{
	    REAL8 cartCoordinates[6] = {r, 0., 0., prstar, pphi/r, 0.};

        REAL8 dHBydpx = XLALSpinHcapExactDerivWRTParam(
            3,
            cartCoordinates,
            seobParams
        );

		omega = XLALSimIMRSpinAlignedEOBPACalculateOmega(
			pol,
			dprstar,
			seobParams,
			LALParams
		);

        flux = XLALSimInspiralEOBPAFluxWrapper(
			r,
			prstar,
			pphi,
			omega,
			seobParams,
			nqcCoeffs,
			LALParams
		);

        result = dpphiBydr*dHBydpx*dprstarBydpr - flux/omega;
	}

	return result;
}

double 
XLALSimInspiralEOBPostAdiabaticdpphiFunc(
	REAL8 pphi_sol,
	void *params
)
{
	struct PostAdiabaticRootSolveParams *pphiParams = (struct PostAdiabaticRootSolveParams *)params;

	REAL8 r = pphiParams->r;
	REAL8 prstar = pphiParams->prstar;
	REAL8 dprstar = pphiParams->dprstar;
	REAL8 pphi = pphi_sol;
	REAL8 dprstarBydr = pphiParams->dprstarBydr;
	REAL8 dprstarBydpr = pphiParams->csi;
	// SpinEOBParams *seobParams = pphiParams->seobParams;
	LALDict *LALParams = pphiParams->LALParams;

	const REAL8 dr = XLALDictLookupREAL8Value(LALParams, "dr");
	const UINT2 analyticFlag = XLALDictLookupUINT2Value(LALParams, "analyticFlag");

	REAL8 pol[4] = {r, 0., prstar, pphi};
	REAL8 omega;
	REAL8 flux;

	REAL8 result;

	if (analyticFlag == 0)
	{
		REAL8 partialHByPartialprstar = XLALSimInspiralEOBPAPartialHByPartialprstar(
			dprstar,
			r,
			prstar,
			pphi,
			pphiParams->seobParams,
			LALParams
		);

		REAL8 partialHByPartialr = XLALSimInspiralEOBPAPartialHByPartialr(
			dr,
			r,
			prstar,
			pphi,
			pphiParams->seobParams,
			LALParams
		);

		omega = XLALSimIMRSpinAlignedEOBPACalculateOmega(
			pol,
			dprstar,
			pphiParams->seobParams,
			LALParams
		);

	    flux = XLALSimInspiralEOBPAFluxWrapper(
			r,
			prstar,
			pphi,
			omega,
			pphiParams->seobParams,
			pphiParams->nqcCoeffs,
			LALParams
		);

		result = dprstarBydr*partialHByPartialprstar + partialHByPartialr - (prstar/pphi)*flux/(omega*dprstarBydpr);
	}
	else
	{
	    REAL8 cartCoordinates[6] = {r, 0., 0., prstar, pphi/r, 0.};

	    REAL8 dHBydx = XLALSpinHcapExactDerivWRTParam(
            0,
            cartCoordinates,
            pphiParams->seobParams
        );

        REAL8 dHBydpx = XLALSpinHcapExactDerivWRTParam(
            3,
            cartCoordinates,
            pphiParams->seobParams
        );

        REAL8 dHBydpy = XLALSpinHcapExactDerivWRTParam(
            4,
            cartCoordinates,
            pphiParams->seobParams
        );

        omega = XLALSimIMRSpinAlignedEOBPACalculateOmega(
			pol,
			dprstar,
			pphiParams->seobParams,
			LALParams
		);

        flux = XLALSimInspiralEOBPAFluxWrapper(
			r,
			prstar,
			pphi,
			omega,
			pphiParams->seobParams,
			pphiParams->nqcCoeffs,
			LALParams
		);

	    result = dprstarBydr*dHBydpx + (dHBydx-dHBydpy*pphi/(r*r)) - (prstar/pphi)*(flux/omega)/dprstarBydpr;
	}

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
			XLAL_ERROR(XLAL_EFUNC, "Derivatives have the wrong sign.");
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

        /* Check to see if the solution is within the tolerance */
        status = gsl_root_test_interval(x_lower, x_upper, absTol, relTol);
    }
    while (iters <= maxIters && status == GSL_CONTINUE);

    result->root = x;
    result->status = status;
    result->nIter = iters;

    if (status != GSL_SUCCESS)
    {
		XLALPrintError("Root finding status: %d\n", status);
		XLAL_ERROR(XLAL_EFUNC);
	}

	gsl_root_fsolver_free (solver);

    return XLAL_SUCCESS;
}

UNUSED REAL8 XLALSimInspiralEOBPACalculateAdibaticParameter(
	REAL8 r,
	REAL8 prstar,
	REAL8 pphi,
	SpinEOBParams *seobParams,
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
	partialHBypartialprstar = XLALSimInspiralEOBPAPartialHByPartialprstar(
		1e-4,
		r,
		prstar,
		pphi,
		seobParams,
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
	struct PostAdiabaticRootSolveParams j0Params;

	for (i = 0; i < rSize; i++)
	{
		DeltaT = XLALSimIMRSpinEOBHamiltonianDeltaT(seobParams->seobCoeffs, rVec->data[i], nu, sigmaKerr);
		DeltaR = XLALSimIMRSpinEOBHamiltonianDeltaR(seobParams->seobCoeffs, rVec->data[i], nu, sigmaKerr);

		csiVec->data[i] = sqrt(DeltaT*DeltaR) / (rVec->data[i]*rVec->data[i]+sigmaKerr*sigmaKerr);

		REAL8 Newtonianj0;
		Newtonianj0 = XLALSimInspiralEOBPACalculateNewtonianj0(rVec->data[i]);

        j0Params.r = rVec->data[i];
        j0Params.prstar = prstarVec->data[i];
        j0Params.seobParams = seobParams;
	    j0Params.LALParams = LALParams;

	    REAL8 pphi0_lower;
	    REAL8 pphi0_upper;

	    pphi0_lower = 0.1 * Newtonianj0;
		pphi0_upper = 1.9 * Newtonianj0;

		struct PostAdiabaticRoot pphiRoot;

	    if (XLALSimInspiralEOBPostAdiabaticRootFinder(
		    	&pphiRoot,
				XLALSimInspiralEOBPostAdiabaticj0Func,
				&j0Params,
				pphi0_lower,
				pphi0_upper,
				ROOT_SOLVER_ABS_TOL,
				ROOT_SOLVER_REL_TOL,
				0
			) != XLAL_SUCCESS)
	    {
	      	XLAL_ERROR(XLAL_EFUNC, "Root solver failed!");
	    }

		pphiVec->data[i] = pphiRoot.root;
        pphi0Vec->data[i] = pphiRoot.root;

		REAL8 polarDynamics[4];
		
		polarDynamics[0] = rVec->data[i];
		polarDynamics[1] = phiVec->data[i];
		polarDynamics[2] = prstarVec->data[i];
		polarDynamics[3] = pphiVec->data[i];

		omegaVec->data[i] = XLALSimIMRSpinAlignedEOBPACalculateOmega(
			polarDynamics,
			dr,
			seobParams,
			LALParams
		);
	}

	REAL8Vector *rReverseVec = XLALCreateREAL8Vector(rSize);
	memset(rReverseVec->data, 0, rReverseVec->length * sizeof(REAL8));

	REAL8Vector *pphiReverseVec = XLALCreateREAL8Vector(rSize);
	memset(pphiReverseVec->data, 0, pphiReverseVec->length * sizeof(REAL8));

	REAL8Vector *dpphiBydrReverseVec = XLALCreateREAL8Vector(rSize);
	memset(dpphiBydrReverseVec->data, 0, dpphiBydrReverseVec->length * sizeof(REAL8));

	XLALReverseREAL8Vector(rVec, rReverseVec);
	XLALReverseREAL8Vector(pphiVec, pphiReverseVec);
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
	const UINT2 analyticFlag = XLALDictLookupUINT2Value(LALParams, "analyticFlag");

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

	REAL8Vector *dprstarVec = XLALCreateREAL8Vector(rSize);
	memset(dprstarVec->data, 0, dprstarVec->length * sizeof(REAL8));

	UINT4 i;

	if (analyticFlag == 0)
	{
		for (i = 0; i < rSize; i++)
		{
			dprstarVec->data[i] = 1.e-4;
		}
	}

	REAL8 polarDynamics[4];

	UINT4 n;
	UINT4 parity = 0;

	for (n = 1; n <= PAOrder; n++)
	{
		parity = n%2;

		if ((n > 1) && (analyticFlag == 0))
		{
			XLALSimInspiralEOBPAMeanValueOrder8(prstarVec, dprstarVec);
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
			    prstarParams.csi = csiVec->data[i];
			    prstarParams.LALParams = LALParams;
			    prstarParams.seobParams = seobParams;
			    prstarParams.nqcCoeffs = nqcCoeffs;
    			REAL8 x_lower;
			    REAL8 x_upper;

				x_upper = 0.8 * prstarVec->data[i];
				x_lower = 1.2 * prstarVec->data[i];

				struct PostAdiabaticRoot prstarRoot;

				if (XLALSimInspiralEOBPostAdiabaticRootFinder(
				    	&prstarRoot,
						XLALSimInspiralEOBPostAdiabaticdprstarFunc,
						&prstarParams,
						x_lower,
						x_upper,
						ROOT_SOLVER_ABS_TOL,
						ROOT_SOLVER_REL_TOL,
						parity
					) != XLAL_SUCCESS)
				{
					XLAL_ERROR(XLAL_EFUNC, "Root finder failed!");
				}

    			prstarVec->data[i] = prstarRoot.root;

				polarDynamics[0] = rVec->data[i];
				polarDynamics[1] = phiVec->data[i];
				polarDynamics[2] = prstarVec->data[i];
				polarDynamics[3] = pphiVec->data[i];

				omegaVec->data[i] = XLALSimIMRSpinAlignedEOBPACalculateOmega(
					polarDynamics,
					dr,
					seobParams,
					LALParams
				);	
	    	}

			XLALReverseREAL8Vector(prstarVec, prstarReverseVec);
			memset(dprstarBydrReverseVec->data, 0, dprstarBydrReverseVec->length * sizeof(REAL8));
			XLALFDDerivative1Order8(rReverseVec, prstarReverseVec, dprstarBydrReverseVec);
			XLALReverseREAL8Vector(dprstarBydrReverseVec, dprstarBydrVec);
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
			    pphiParams.dprstarBydr = dprstarBydrVec->data[i];
			    pphiParams.LALParams = LALParams;
			    pphiParams.seobParams = seobParams;
			    pphiParams.nqcCoeffs = nqcCoeffs;

			    REAL8 x_lower = 0.9 * pphiVec->data[i];
				REAL8 x_upper = 1.1 * pphiVec->data[i];

				if (XLALSimInspiralEOBPostAdiabaticRootFinder(
				    	&pphiRoot,
						XLALSimInspiralEOBPostAdiabaticdpphiFunc,
						&pphiParams,
						x_lower,
						x_upper,
						ROOT_SOLVER_ABS_TOL,
						ROOT_SOLVER_REL_TOL,
						parity
					) != XLAL_SUCCESS)
				{
					XLAL_ERROR(XLAL_EFUNC, "Root finder failed!");
				}

				pphiVec->data[i] = pphiRoot.root;
		
				polarDynamics[0] = rVec->data[i];
				polarDynamics[1] = phiVec->data[i];
				polarDynamics[2] = prstarVec->data[i];
				polarDynamics[3] = pphiVec->data[i];

				omegaVec->data[i] = XLALSimIMRSpinAlignedEOBPACalculateOmega(
					polarDynamics,
					dr,
					seobParams,
					LALParams
				);
	    	}

			XLALReverseREAL8Vector(pphiVec, pphiReverseVec);
			memset(dpphiBydrReverseVec->data, 0, dpphiBydrReverseVec->length * sizeof(REAL8));
			XLALFDDerivative1Order8(rReverseVec, pphiReverseVec, dpphiBydrReverseVec);
			XLALReverseREAL8Vector(dpphiBydrReverseVec, dpphiBydrVec);
		}
	}

	if ((parity) && (analyticFlag == 0))
	{
		XLALSimInspiralEOBPAMeanValueOrder8(prstarVec, dprstarVec);
	}

	// Loop over and compute the derivatives we need to do quadratures later
	REAL8 partialHBypartialprstar;

	for(i = 0; i < rSize; i++)
	{
		partialHBypartialprstar = XLALSimInspiralEOBPAPartialHByPartialprstar(
			dprstarVec->data[i],
			rVec->data[i],
			prstarVec->data[i],
			pphiVec->data[i],
			seobParams,
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
XLALSimInspiralEOBPAPartialHByPartialr(
	REAL8 h,
	REAL8 r,
	REAL8 prstar,
	REAL8 pphi,
	SpinEOBParams *seobParams,
	LALDict *LALParams
)
{
	const UINT2 analyticFlag = XLALDictLookupUINT2Value(LALParams, "analyticFlag");

	REAL8 partialHByPartialr = 0.;

	if (analyticFlag == 0)
	{
		REAL8 coeffs[9] = {1./280., -4./105., 1./5., -4./5., 0, 4./5., -1./5., 4./105., -1./280.};

		INT4 i;

		for (i = -4; i <= 4; i++)
		{
			if (i != 0)
			{
				partialHByPartialr += coeffs[i+4] * XLALSimInspiralEOBPAHamiltonianWrapper(
					r + i*h,
					prstar,
					pphi,
					seobParams->seobCoeffs,
					LALParams
				);
			}
		}

		partialHByPartialr /= h;
	}
	else
	{
		REAL8 cartCoordinates[6] = {r, 0., 0., prstar, pphi/r, 0.};

	    REAL8 dHBydx = XLALSpinHcapExactDerivWRTParam(
            0,
            cartCoordinates,
            seobParams
        );

        REAL8 dHBydpy = XLALSpinHcapExactDerivWRTParam(
            4,
            cartCoordinates,
            seobParams
        );

        partialHByPartialr = dHBydx - dHBydpy*pphi/(r*r);
	}

	return partialHByPartialr;
}

REAL8
XLALSimInspiralEOBPAPartialHByPartialprstar(
	REAL8 h,
	REAL8 r,
	REAL8 prstar,
	REAL8 pphi,
	SpinEOBParams *seobParams,
	LALDict *LALParams
)
{
	REAL8 coeffs[9] = {1./280., -4./105., 1./5., -4./5., 0, 4./5., -1./5., 4./105., -1./280.};

	const UINT2 analyticFlag = XLALDictLookupUINT2Value(LALParams, "analyticFlag");

	REAL8 partialHByPartialprstar = 0.;

	if (analyticFlag == 0)
	{
		INT4 i;
	
		for (i = -4; i <= 4; i++)
		{
			if (i != 0)
			{
				partialHByPartialprstar += coeffs[i+4] * XLALSimInspiralEOBPAHamiltonianWrapper(
					r,
					prstar + i*h,
					pphi,
					seobParams->seobCoeffs,
					LALParams
				);
			}
		}
	
		partialHByPartialprstar /= h;
	}
	else
	{
		REAL8 cartCoordinates[6] = {r, 0., 0., prstar, pphi/r, 0.};

	    REAL8 dHBydpx = XLALSpinHcapExactDerivWRTParam(
            3,
            cartCoordinates,
            seobParams
        );

		partialHByPartialprstar = dHBydpx;
	}

	return partialHByPartialprstar;
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

		meanVec->data[i] /= 8.0;
	}

	return XLAL_SUCCESS;
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
	const UINT2 analyticFlag = XLALDictLookupUINT2Value(LALParams, "analyticFlag");

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

    if (analyticFlag == 0)
    {
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
	}
	else
	{
		H = XLALSimIMRSpinEOBHamiltonianOptimized(
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
	const UINT2 analyticFlag = XLALDictLookupUINT2Value(LALParams, "analyticFlag");

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
	REAL8 pol[4] = {r, 0., prstar, pphi};
	polarDynamics.data = pol;
	
	const UINT4 lMax = 8;

	const UINT4 SpinAlignedEOBversion = 4;
	
 	REAL8 Flux;

    if (analyticFlag == 0)
    {
	    Flux = XLALInspiralSpinFactorizedFlux_PA(
	    			&polarDynamics,
	    			nqcCoeffs,
	    			omega,
	    			seobParams,
	    			H,
	    			lMax,
	    			SpinAlignedEOBversion
	    		);
	}
	else
	{
		Flux = XLALInspiralSpinFactorizedFluxOptimized(
					&polarDynamics,
	    			nqcCoeffs,
	    			omega,
	    			seobParams,
	    			H,
	    			lMax,
	    			SpinAlignedEOBversion
	    		);
	}

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
	/**<< initial values (r, phi, pr, pphi) for the post-adiabatic routine */
	UINT4 SpinAlignedEOBversion,
	/**<< 1 for SEOBNRv1, 2 for SEOBNRv2, 4 for SEOBNRv4, 201 for SEOBNRv2T, 401 for SEOBNRv4T, 41 for SEOBNRv4HM */
	SpinEOBParams *seobParams,
	/**<< SEOB parameters */
	EOBNonQCCoeffs *nqcCoeffs,
	/**<< NQC coefficients */
	LALDict *PAParams
	/**<< pointer to a dictionary containing parameters for the post-adiabatic routine */
)
{
	if (SpinAlignedEOBversion != 4)
	{
		XLAL_ERROR(XLAL_EFUNC, "XLALSimInspiralEOBPostAdiabatic can only be used with SpinAlignedEOBversion = 4.\n");
	}

	REAL8 q = XLALSimInspiralEOBPACalculateMassRatio(m1, m2);
	REAL8 nu = XLALSimInspiralEOBPACalculateSymmetricMassRatio(q);

	REAL8 chi1 = spin1z;
	REAL8 chi2 = spin2z;

	const UINT4 PAOrder = XLALDictLookupUINT4Value(PAParams, "PAOrder");

	const UINT2 analyticFlag = XLALDictLookupUINT2Value(PAParams, "analyticFlag");

	REAL8 X1 = XLALSimInspiralEOBPACalculateX1(nu);
	REAL8 X2 = XLALSimInspiralEOBPACalculateX2(nu);

	REAL8 a1 = XLALSimInspiralEOBPACalculatea(X1, chi1);
	REAL8 a2 = XLALSimInspiralEOBPACalculatea(X2, chi2);

	REAL8 aK = a1 + a2;

	REAL8 Sstar = XLALSimInspiralEOBPACalculateSstar(X1, X2, chi1, chi2);

	REAL8 rInitial = initVals.data[0];

	// REAL8 rFinalPrefactor = XLALDictLookupREAL8Value(PAParams, "rFinal");
	REAL8 rFinalPrefactor = 1.8;
	REAL8 rFinal = rFinalPrefactor * XLALSimInspiralEOBPostAdiabaticFinalRadius(q, a1, a2);

	if (rInitial <= rFinal)
	{
		XLAL_ERROR(XLAL_ERANGE);
	}

	// REAL8 rSwitchPrefactor = XLALDictLookupREAL8Value(PAParams, "rSwitch");
	REAL8 rSwitchPrefactor = 1.8;
	REAL8 rSwitch = rSwitchPrefactor * XLALSimInspiralEOBPostAdiabaticFinalRadius(q, a1, a2);

	if (rInitial <= rSwitch)
	{
		XLAL_ERROR(XLAL_ERANGE);
	}

	// UINT4 rSize = 100;
	// UINT4 rSize = XLALDictLookupUINT4Value(PAParams, "rSize");
	
	// REAL8 dr = XLALSimInspiralEOBPACalculatedr(rInitial, rFinal, rSize);

	REAL8 dr = 0.3;
	
	UINT4 rSize = (int) ceil((rInitial-rFinal)/dr);
	
	if (rSize <= 4)
	{
		XLAL_ERROR(XLAL_ERANGE);
	}
	else if (rSize < 10)
	{
		rSize = 10;
	}
	
	dr = XLALSimInspiralEOBPACalculatedr(rInitial, rFinal, rSize);

	LALDict *LALparams = XLALCreateDict();

	XLALDictInsertREAL8Value(LALparams, "q", q);
	XLALDictInsertREAL8Value(LALparams, "nu", nu);

	XLALDictInsertREAL8Value(LALparams, "chi1", chi1);
	XLALDictInsertREAL8Value(LALparams, "chi2", chi2);

	XLALDictInsertUINT4Value(LALparams, "PAOrder", PAOrder);
	XLALDictInsertUINT2Value(LALparams, "analyticFlag", analyticFlag);

	XLALDictInsertREAL8Value(LALparams, "X1", X1);
	XLALDictInsertREAL8Value(LALparams, "X2", X2);

	XLALDictInsertREAL8Value(LALparams, "a1", a1);
	XLALDictInsertREAL8Value(LALparams, "a2", a2);
	XLALDictInsertREAL8Value(LALparams, "aK", aK);

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

 	if (XLALSimInspiralEOBPACalculateAdiabaticDynamics(
			rVec,
			phiVec,
			prstarVec,
			pphiVec,
			pphi0Vec,
			dpphiBydrVec,
			csiVec,
			omegaVec,
			seobParams,
			LALparams
		) != XLAL_SUCCESS)
	{
		XLAL_ERROR(XLAL_EFUNC, "XLALSimInspiralEOBPACalculateAdiabaticDynamics failed.");
	}

 	XLALReverseREAL8Vector(rVec, rReverseVec);

	if (PAOrder > 0)
	{
		if (XLALSimInspiralEOBPACalculatePostAdiabaticDynamics(
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
			) != XLAL_SUCCESS)
		{
			XLAL_ERROR(XLAL_EFUNC, "XLALSimInspiralEOBPACalculatePostAdiabaticDynamics failed.");
		}
	}

	XLALCumulativeIntegral3(rVec, dtBydrVec, tVec);

	XLALCumulativeIntegral3(rVec, dphiBydrVec, phiVec);

    UINT4 i, j;
	XLALReverseREAL8Vector(omegaVec, omegaReverseVec);
	XLALFDDerivative1Order8(rReverseVec, omegaReverseVec, domegadrReverseVec);
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
	// 	adiabatic_param = XLALSimInspiralEOBPACalculateAdibaticParameter(r,prstar,pphi,seobParams,csi,LALparams, omega, domegadr);
	// 	adiabatic_param_Vec->data[j]=adiabatic_param;
	// 	//printf("r=%.17f, omega=%.17f,domegadr=%.17f, adibatic_param = %.17f\n",r,omega,domegadr,adiabatic_param_Vec->data[j]);
	// 	if(idx_stop==0 && adiabatic_param>PA_AD_THRS){
	// 		printf("r=%.18e Q = %.18e\n",r,adiabatic_param);
	// 		idx_stop = j;
	// 	}
	// }

	for (j=0; j<rSize; j++)
	{
		r = rVec->data[j];

		if (idx_stop == 0 && rSwitch >= r)
		{
			idx_stop = j;
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

	// This needs to be revisited
	if (outSize < 4)
	{
		XLAL_ERROR(XLAL_EFUNC, "PA inspiral doesn't have enough points for interpolation.");
	}
	// please fix before the end
	
    return XLAL_SUCCESS;
}
