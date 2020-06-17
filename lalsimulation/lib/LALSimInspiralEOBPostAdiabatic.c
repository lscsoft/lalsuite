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



double
XLALSimInspiralEOBPACalculateSymmetricMassRatio(
	const REAL8 q
	/**< mass ratio */)
{
	REAL8 nu = 0.;

	if (q > 0.)
	{
		nu = q / ((q+1.)*(q+1.));
	}

	return nu;
}

double
XLALSimInspiralEOBPostAdiabaticX1(
	const REAL8 nu
	/**< symmetric mass ratio */)
{
	if ( (nu < 0.) || (nu > 0.25) )
	{
		XLALPrintError("XLAL Error - %s: Symmetric mass ratio is 0 <= nu <= 1/4.\n", __func__);
		XLAL_ERROR(XLAL_EINVAL);
	}

	REAL8 X1 = 0.5 * (1.+sqrt(1.-4.*nu));

	return X1;
}


// Potentially remove or change at a later stage
double
XLALSimInspiralEOBPostAdiabaticX2(
	const REAL8 nu
	/**< symmetric mass ratio */)
{
	REAL8 X1 = XLALSimInspiralEOBPostAdiabaticX1(nu);

	REAL8 X2 = 1. - X1;

	return X2;
}

double
XLALSimInspiralEOBPostAdiabaticTimeUnitsFactor(
	REAL8 M
	/**< Mass in SI units */)
{
	REAL8 time_units_factor = 1. / (M*LAL_MTSUN_SI);

	return time_units_factor;
}

double
XLALSimInspiralEOBPostAdiabaticDynr0Kepler(
	REAL8 f0
	/**< Mass in SI units */)
{
	const REAL8 omg_orb0 = LAL_PI * f0; // = 2 * Pi * (f0/2)

	REAL8 r0 = pow(omg_orb0, -2./3.);

	return r0;
}

REAL8
XLALSimInspiralEOBPostAdiabaticTotalSpin(
	REAL8 q,
	REAL8 a1,
	REAL8 a2)
{
	REAL8 invQ;
	REAL8 onePlusInvQ;
	REAL8 aTotal;

	invQ = 1. / q;
	onePlusInvQ = 1. + invQ;

	aTotal = (a1 + a2*invQ*invQ) / (onePlusInvQ * onePlusInvQ);

	return aTotal;
}

REAL8
XLALSimInspiralEOBPostAdiabaticFinalRadius(
	REAL8 q,
	REAL8 a1,
	REAL8 a2)
{
	REAL8 aTotal;

	REAL8 rISCO;
	REAL8 finalRadiusPrefactor;
	REAL8 rFinal;

	aTotal = XLALSimInspiralEOBPostAdiabaticTotalSpin(q, a1, a2);

	rISCO = XLALSimRadiusKerrISCO(aTotal);
	finalRadiusPrefactor = 1.35;

	rFinal = finalRadiusPrefactor * rISCO;

	return rFinal;
}

/* 
 * Fit of c3, TEOBResumS paper Nagar et al. (2018) 
 * c3 = 0 with tides
 */
double
XLALSimInspiralEOBPostAdiabaticFitGlobalc3(
	REAL8 nu,
	REAL8 a1,
	REAL8 a2)
{  
    const REAL8 nu2 = nu * nu;
    const REAL8 nu3 = nu2 * nu;
    const REAL8 X12 = sqrt(1. - 4.*nu);
    const REAL8 a12 = a1 + a2;
    
    /* Equal-mass, equal-spin coefficients */
    const REAL8 c0 = 43.371638;
    const REAL8 n1 = -1.174839;
    const REAL8 n2 =  0.354064;
    const REAL8 d1 = -0.151961;
	
    const REAL8 c3_eq = c0 * (1.+n1*a12+n2*a12*a12) / (1.+d1*a12);
  
    /* --- */
    const REAL8 cnu = 929.579;
    const REAL8 cnu2 = -9178.87;
    const REAL8 cnu3 = 23632.3;
    const REAL8 ca1_a2 = -104.891;
  
    const REAL8 c3_uneq = cnu*a12*nu*X12 + cnu2*a12*nu2*X12 + cnu3*a12*nu3*X12 + ca1_a2*(a1-a2)*nu2;
  
    const REAL8 c3 = c3_eq + c3_uneq;

    return c3;
}

double
XLALSimInspiralEOBPostAdiabaticz3(
	const REAL8 nu
	/**< symmetric mass ratio */)
{
	REAL8 z3 = 2. * nu * (4.-3.*nu);

	return z3;
}

int
XLALSimInspiralEOBPACalculateRadialGrid(
	REAL8Vector *rVec,
	LALDict *LALParams)
{
	const REAL8 rInitial = XLALDictLookupREAL8Value(LALParams, "rInitial");
	const UINT4 rSize = XLALDictLookupREAL8Value(LALParams, "rSize");
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
	SpinEOBParams *seobParams = prstarParams->seobParams;
	EOBNonQCCoeffs *nqcCoeffs = prstarParams->nqcCoeffs;

	SpinEOBHCoeffs *seobCoeffs = seobParams->seobCoeffs;

	REAL8 dpphiBydr;
	dpphiBydr = prstarParams->dpphiBydr;

	LALDict *LALParams = prstarParams->LALParams;

	REAL8 partialHBypartialprstar;
	partialHBypartialprstar = XLALSimInspiralEOBPAHamiltonianPartialDerivativeprstar(
																					1.,
																					prstarParams->r,
																					prstar_sol,
																					prstarParams->pphi,
																					seobCoeffs,
																					LALParams
																				);
	REAL8 H;
	H = XLALSimInspiralEOBPAHamiltonianWrapper(
							prstarParams->r,
							prstar_sol,
							prstarParams->pphi,
							seobCoeffs,
							LALParams
						);

	REAL8 flux;
    flux = XLALSimInspiralEOBPAFluxWrapper(
				prstarParams->r,
				prstar_sol,
				prstarParams->pphi,
				prstarParams->omega,
				H,
				seobParams,
				nqcCoeffs,
				LALParams
			);

	REAL8 result;

	result = dpphiBydr*partialHBypartialprstar*prstarParams->csi - flux;

	return result;
}

double 
XLALSimInspiralEOBPostAdiabaticdpphiFunc(
	REAL8 pphi_sol,
	void *params)
{
	struct PostAdiabaticRootSolveParams *pphiParams = (struct PostAdiabaticRootSolveParams *)params;

	SpinEOBHCoeffs *seobCoeffs = pphiParams->seobParams->seobCoeffs;

	REAL8 r = pphiParams->r;
	REAL8 prstar = pphiParams->prstar;
	REAL8 pphi = pphi_sol;
	const REAL8 dr = XLALDictLookupREAL8Value(pphiParams->LALParams, "dr");
	LALDict *LALParams = pphiParams->LALParams;

	REAL8 partialHBypartialprstar;
	partialHBypartialprstar = XLALSimInspiralEOBPAHamiltonianPartialDerivativeprstar(
									1.,
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

REAL8
XLALSimInspiralEOBPANewtonianj0(
	REAL8 r)
{
	REAL8 Newtonianj0;
	Newtonianj0 = sqrt(r);

	return Newtonianj0;
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

REAL8Vector
XLALPostAdiabaticSplineDerivative(
	REAL8Vector *VecX,
	REAL8Vector *VecY)
{
	UINT4 vecLength;
	vecLength = VecX->length;

	REAL8Vector *splineDerivative = XLALCreateREAL8Vector(vecLength);
	memset(splineDerivative->data, 0, splineDerivative->length * sizeof(REAL8));

    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_akima, vecLength);

    gsl_spline_init(spline, VecX->data, VecY->data, vecLength);

    UINT4 i;

    for (i = 0; i < vecLength; i++)
    {
        splineDerivative->data[i] = gsl_spline_eval_deriv(spline, VecX->data[i], acc);
    }

    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);

    return *splineDerivative;
}

REAL8Vector
XLALFDDerivative1Order4(
	REAL8Vector *XVec,
	REAL8Vector *YVec)
{
	REAL8 fourthOrderCoeffs[5][5] = {
		{-25./12., 4., -3., 4./3., -1./4.},
		{-1./4., -5./6., 3./2., -1./2., 1./12.},
		{1./12., -2./3., 0, 2./3., -1./12.},
		{-1./12., 1./2., -3./2., 5./6., 1./4.},
		{1./4., -4./3., 3., -4., 25./12.}
	};


	UINT4 vecLength;
	vecLength = XVec->length;

	REAL8Vector *derivativeVec = XLALCreateREAL8Vector(vecLength);
	memset(derivativeVec->data, 0, derivativeVec->length * sizeof(REAL8));

	REAL8 h;
	h = fabs(XVec->data[0] - XVec->data[1]);

	UINT4 i;
	UINT4 j;

	for (i = 0; i <= vecLength-1; i++)
	{
		if (i == 0)
		{
			for (j = 0; j <= 4; j++)
			{
				derivativeVec->data[i] += fourthOrderCoeffs[0][j]*YVec->data[j];
			}
		}
		else if (i == 1)
		{
			for (j = 0; j <= 4; j++)
			{
				derivativeVec->data[i] += fourthOrderCoeffs[1][j]*YVec->data[j];
			}
		}
		else if (i == vecLength-2)
		{
			for (j = 0; j <= 4; j++)
			{
				derivativeVec->data[i] += fourthOrderCoeffs[3][j]*YVec->data[i+j-3];
			}
		}
		else if (i == vecLength-1)
		{
			for (j = 0; j <= 4; j++)
			{
				derivativeVec->data[i] += fourthOrderCoeffs[4][j]*YVec->data[i+j-4];
			}
		}
		else
		{
			for (j = 0; j <= 4; j++)
			{
				derivativeVec->data[i] += fourthOrderCoeffs[2][j]*YVec->data[i+j-2];
			}
		}

		derivativeVec->data[i] /= h;
	}
	
    return *derivativeVec;
}

REAL8Vector
XLALFDDerivative1Order2(
	REAL8Vector *XVec,
	REAL8Vector *YVec)
{
	REAL8 secondOrderCoeffs[3][3] = {
		{-3./2., 2, -1./2.},
		{-1./2., 0, 1./2.},
		{1./2., -2., 3./2.}
	};

	UINT4 vecLength;
	vecLength = XVec->length;

	REAL8Vector *derivativeVec = XLALCreateREAL8Vector(vecLength);
	memset(derivativeVec->data, 0, derivativeVec->length * sizeof(REAL8));

	REAL8 h;
	h = fabs(XVec->data[0] - XVec->data[1]);

	UINT4 i;
	UINT4 j;

	for (i = 0; i <= vecLength-1; i++)
	{
		if (i == 0)
		{
			for (j = 0; j <= 2; j++)
			{
				derivativeVec->data[i] += secondOrderCoeffs[0][j] * YVec->data[j];
			}
		}
		else if (i == vecLength-1)
		{
			for (j = 0; j <= 2; j++)
			{
				derivativeVec->data[i] += secondOrderCoeffs[2][j] * YVec->data[i+j-2];
			}
		}
		else
		{
			for (j = 0; j <= 2; j++)
			{
				derivativeVec->data[i] += secondOrderCoeffs[1][j] * YVec->data[i+j-1];
			}
		}

		derivativeVec->data[i] /= h;
	}
	
    return *derivativeVec;
}

REAL8Vector
XLALFDDerivative1Order6(
	REAL8Vector *XVec,
	REAL8Vector *YVec)
{
	REAL8 sixthOrderCoeffs[7][7] = {
		{-49./20., 6., -15./2., 20./3., -15./4., 6./5., -1./6.},
		{-1./6., -77./60., 5./2., -5./3., 5./6., -1./4., 1./30.},
		{1./30., -2./5., -7./12., 4./3., -1./2., 2./15., -1./60.},
		{-1./60., 3./20., -3./4., 0, 3./4., -3./20., 1./60.},
		{1./60., -2./15., 1./2., -4./3., 7./12., 2./5., -1./30.},
		{-1./30., 1./4., -5./6., 5./3., -5./2., 77./60., 1./6.},
		{1./6., -6./5., 15./4., -20./3., 15./2., -6., 49./20.}
	};

	UINT4 vecLength;
	vecLength = XVec->length;

	REAL8Vector *derivativeVec = XLALCreateREAL8Vector(vecLength);
	memset(derivativeVec->data, 0, derivativeVec->length * sizeof(REAL8));

	REAL8 h;
	h = fabs(XVec->data[0] - XVec->data[1]);

	UINT4 i;
	UINT4 j;

	for (i = 0; i <= vecLength-1; i++)
	{
		if (i == 0)
		{
			for (j = 0; j <= 6; j++)
			{
				derivativeVec->data[i] += sixthOrderCoeffs[0][j] * YVec->data[j];
			}
		}
		else if (i == 1)
		{
			for (j = 0; j <= 6; j++)
			{
				derivativeVec->data[i] += sixthOrderCoeffs[1][j] * YVec->data[j];
			}
		}
		else if (i == 2)
		{
			for (j = 0; j <= 6; j++)
			{
				derivativeVec->data[i] += sixthOrderCoeffs[2][j] * YVec->data[j];
			}
		}
		else if (i == vecLength-3)
		{
			for (j = 0; j <= 6; j++)
			{
				derivativeVec->data[i] += sixthOrderCoeffs[4][j] * YVec->data[i+j-4];
			}
		}
		else if (i == vecLength-2)
		{
			for (j = 0; j <= 6; j++)
			{
				derivativeVec->data[i] += sixthOrderCoeffs[5][j] * YVec->data[i+j-5];
			}
		}
		else if (i == vecLength-1)
		{
			for (j = 0; j <= 6; j++)
			{
				derivativeVec->data[i] += sixthOrderCoeffs[6][j] * YVec->data[i+j-6];
			}
		}
		else
		{
			for (j = 0; j <= 6; j++)
			{
				derivativeVec->data[i] += sixthOrderCoeffs[3][j] * YVec->data[i+j-3];
			}
		}

		derivativeVec->data[i] /= h;
	}
	
    return *derivativeVec;
}

REAL8Vector
XLALFDDerivative1Order8(
	REAL8Vector *XVec,
	REAL8Vector *YVec)
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

	UINT4 vecLength;
	vecLength = XVec->length;

	REAL8Vector *derivativeVec = XLALCreateREAL8Vector(vecLength);
	memset(derivativeVec->data, 0, derivativeVec->length * sizeof(REAL8));

	REAL8 h;
	h = fabs(XVec->data[0] - XVec->data[1]);

	UINT4 i;
	UINT4 j;

	for (i = 0; i <= vecLength-1; i++)
	{
		if (i == 0)
		{
			for (j = 0; j <= 8; j++)
			{
				derivativeVec->data[i] += eightOrderCoeffs[0][j] * YVec->data[j];
			}
		}
		else if (i == 1)
		{
			for (j = 0; j <= 8; j++)
			{
				derivativeVec->data[i] += eightOrderCoeffs[1][j] * YVec->data[j];
			}
		}
		else if (i == 2)
		{
			for (j = 0; j <= 8; j++)
			{
				derivativeVec->data[i] += eightOrderCoeffs[2][j] * YVec->data[j];
			}
		}
		else if (i == 3)
		{
			for (j = 0; j <= 8; j++)
			{
				derivativeVec->data[i] += eightOrderCoeffs[3][j] * YVec->data[j];
			}
		}
		else if (i == vecLength-4)
		{
			for (j = 0; j <= 8; j++)
			{
				derivativeVec->data[i] += eightOrderCoeffs[5][j] * YVec->data[i+j-5];
			}
		}
		else if (i == vecLength-3)
		{
			for (j = 0; j <= 8; j++)
			{
				derivativeVec->data[i] += eightOrderCoeffs[6][j] * YVec->data[i+j-6];
			}
		}
		else if (i == vecLength-2)
		{
			for (j = 0; j <= 8; j++)
			{
				derivativeVec->data[i] += eightOrderCoeffs[7][j] * YVec->data[i+j-7];
			}
		}
		else if (i == vecLength-1)
		{
			for (j = 0; j <= 8; j++)
			{
				derivativeVec->data[i] += eightOrderCoeffs[8][j] * YVec->data[i+j-8];
			}
		}
		else
		{
			for (j = 0; j <= 8; j++)
			{
				derivativeVec->data[i] += eightOrderCoeffs[4][j] * YVec->data[i+j-4];
			}
		}

		derivativeVec->data[i] /= h;
	}
	
    return *derivativeVec;
}

REAL8Vector
XLALCumulativeIntegral3(
	REAL8Vector *XVec,
	REAL8Vector *YVec)
{
	UINT4 vecLength;
	vecLength = XVec->length;

	REAL8Vector *XVecExt = XLALCreateREAL8Vector(vecLength+2);
	REAL8Vector *YVecExt = XLALCreateREAL8Vector(vecLength+2);

	memset(XVecExt->data, 0, XVecExt->length * sizeof(REAL8));
	memset(YVecExt->data, 0, YVecExt->length * sizeof(REAL8));

	REAL8 *X0, *X1, *X2, *X3;
	REAL8 *Y0, *Y1, *Y2, *Y3;

	REAL8 a, b, c, d, e, h, g, z;
	const REAL8 oo12 = 0.08333333333333333;

	UINT4 i;

	for (i=1; i < vecLength+1; i++)
	{
		XVecExt->data[i] = XVec->data[i-1];
		YVecExt->data[i] = YVec->data[i-1];
	}

	XVecExt->data[0] = XVec->data[3];
	XVecExt->data[vecLength+1] = XVec->data[vecLength-4];

	YVecExt->data[0] = YVec->data[3];
	YVecExt->data[vecLength+1] = YVec->data[vecLength-4];

	X0 = &XVecExt->data[0];
	X1 = &XVecExt->data[1];
	X2 = &XVecExt->data[2];
	X3 = &XVecExt->data[3];

	Y0 = &YVecExt->data[0];
	Y1 = &YVecExt->data[1];
	Y2 = &YVecExt->data[2];
	Y3 = &YVecExt->data[3];

	REAL8Vector *integralVec = XLALCreateREAL8Vector(vecLength);
	memset(integralVec->data, 0, integralVec->length * sizeof(REAL8));

	for (i=0; i < vecLength-1; i++)
	{
		a = X1[i] - X0[i];
		b = X2[i] - X1[i];
		c = X3[i] - X2[i];
		d = Y1[i] - Y0[i];
		e = Y2[i] - Y1[i];
		h = Y3[i] - Y2[i];
		g = 0.5 * (Y1[i]+Y2[i]);
		z = b*g + oo12*b*b*(c*b*(2*c+b)*(c+b)*d-a*c*(c-a)*(2*c+2*a+3*b)*e-a*b*(2*a+b)*(a+b)*h)/(a*c*(a+b)*(c+b)*(c+a+b));
		integralVec->data[i+1] = integralVec->data[i] + z;
	}

	return *integralVec;
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
		printf("XLALSimInspiralEOBPostAdiabatic can only be used with SpinAlignedEOBversion = 4.");
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

	REAL8Vector *HVec = XLALCreateREAL8Vector(rSize);
	memset(HVec->data, 0, HVec->length * sizeof(REAL8));

	REAL8Vector *csiVec = XLALCreateREAL8Vector(rSize);
	memset(csiVec->data, 0, csiVec->length * sizeof(REAL8));

	REAL8Vector *omegaVec = XLALCreateREAL8Vector(rSize);
	memset(omegaVec->data, 0, omegaVec->length * sizeof(REAL8));

	REAL8Vector *fluxVec = XLALCreateREAL8Vector(rSize);
	memset(fluxVec->data, 0, fluxVec->length * sizeof(REAL8));

	REAL8Vector *AVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *BVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *dAVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *d2AVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *dBVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *rcVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *drcBydrVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *d2rBydr2Vec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *dGBydrVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *dGBydprstarVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *dGBydprstarbyprstarVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *G0Vec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *dGBydr0Vec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *sqrtAByBVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *uc2Vec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *ducBydrVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *dAuc2BydrVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *pphiVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *pphi0Vec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *prstarVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *dprstarBydrVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *HeffVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *HeffOrbVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *dHeffBydpphiVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *EVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *OmgVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *OmgOrbVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *dtBydrVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *dphiBydrVec = XLALCreateREAL8Vector(rSize);
    

    REAL8Vector *rReverseVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *pphiReverseVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *dpphiBydrReverseVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *dpphiBydrVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *dpphiBydr0Vec = XLALCreateREAL8Vector(rSize);

    REAL8Vector *prstarReverseVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *dprstarBydrReverseVec = XLALCreateREAL8Vector(rSize);

    memset(AVec->data, 0, AVec->length * sizeof(REAL8));
    memset(BVec->data, 0, BVec->length * sizeof(REAL8));
    memset(dAVec->data, 0, dAVec->length * sizeof(REAL8));
    memset(d2AVec->data, 0, d2AVec->length * sizeof(REAL8));
    memset(dBVec->data, 0, dBVec->length * sizeof(REAL8));
    memset(rcVec->data, 0, rcVec->length * sizeof(REAL8));
    memset(drcBydrVec->data, 0, drcBydrVec->length * sizeof(REAL8));
    memset(d2rBydr2Vec->data, 0, d2rBydr2Vec->length * sizeof(REAL8));
    memset(dGBydrVec->data, 0, dGBydrVec->length * sizeof(REAL8));
    memset(dGBydprstarVec->data, 0, dGBydprstarVec->length * sizeof(REAL8));
    memset(dGBydprstarbyprstarVec->data, 0, dGBydprstarbyprstarVec->length * sizeof(REAL8));
    memset(G0Vec->data, 0, G0Vec->length * sizeof(REAL8));
    memset(dGBydr0Vec->data, 0, dGBydr0Vec->length * sizeof(REAL8));
    memset(sqrtAByBVec->data, 0, sqrtAByBVec->length * sizeof(REAL8));
    memset(uc2Vec->data, 0, uc2Vec->length * sizeof(REAL8));
    memset(ducBydrVec->data, 0, ducBydrVec->length * sizeof(REAL8));
    memset(dAuc2BydrVec->data, 0, dAuc2BydrVec->length * sizeof(REAL8));
    memset(pphiVec->data, 0, pphiVec->length * sizeof(REAL8));
    memset(pphi0Vec->data, 0, pphi0Vec->length * sizeof(REAL8));
    memset(prstarVec->data, 0, prstarVec->length * sizeof(REAL8));
    memset(dprstarBydrVec->data, 0, dprstarBydrVec->length * sizeof(REAL8));
    memset(HeffVec->data, 0, HeffVec->length * sizeof(REAL8));
    memset(HeffOrbVec->data, 0, HeffOrbVec->length * sizeof(REAL8));
    memset(dHeffBydpphiVec->data, 0, dHeffBydpphiVec->length * sizeof(REAL8));
    memset(EVec->data, 0, EVec->length * sizeof(REAL8));
    memset(OmgVec->data, 0, OmgVec->length * sizeof(REAL8));
    memset(OmgOrbVec->data, 0, OmgOrbVec->length * sizeof(REAL8));
    memset(dtBydrVec->data, 0, dtBydrVec->length * sizeof(REAL8));
    memset(dphiBydrVec->data, 0, dphiBydrVec->length * sizeof(REAL8));
    
    memset(rReverseVec->data, 0, rReverseVec->length * sizeof(REAL8));
    memset(pphiReverseVec->data, 0, pphiReverseVec->length * sizeof(REAL8));
    memset(dpphiBydrReverseVec->data, 0, dpphiBydrReverseVec->length * sizeof(REAL8));
    memset(dpphiBydrVec->data, 0, dpphiBydrVec->length * sizeof(REAL8));
    memset(dpphiBydr0Vec->data, 0, dpphiBydr0Vec->length * sizeof(REAL8));

    memset(prstarReverseVec->data, 0, prstarReverseVec->length * sizeof(REAL8));
    memset(dprstarBydrReverseVec->data, 0, dprstarBydrReverseVec->length * sizeof(REAL8));

	// REAL8 r;

	UINT4 i;

	SpinEOBHCoeffs *seobCoeffs = seobParams->seobCoeffs;

 	// compute r grid
 	XLALSimInspiralEOBPACalculateRadialGrid(
 		rVec,
 		LALparams
 	);

	// compute SEOBNRv4 Hamiltonian
	for (i = 0; i < rSize; i++)
	{
		HVec->data[i] = XLALSimInspiralEOBPAHamiltonianWrapper(
							rVec->data[i],
							prstarVec->data[i],
							pphiVec->data[i],
							seobCoeffs,
							LALparams
						);
		// printf("%.18e\n", HVec->data[i]);
	}

	// compute j0
	for (i = 0; i < rSize; i++)
	{
		REAL8 Newtonianj0;
		Newtonianj0 = XLALSimInspiralEOBPANewtonianj0(rVec->data[i]);

		struct PostAdiabaticRoot j0Root;

        struct PostAdiabaticRootSolveParams j0Params;
        j0Params.r = rVec->data[i];
        j0Params.prstar = prstarVec->data[i];
        j0Params.seobParams = seobParams;
	    j0Params.LALParams = LALparams;

	    REAL8 pphi0_lower;
	    REAL8 pphi0_upper;

	    pphi0_lower = 0.1 * Newtonianj0;
		pphi0_upper = 1.9 * Newtonianj0;

	    XLALSimInspiralEOBPostAdiabaticRootFinder(
	    	&j0Root,
			XLALSimInspiralEOBPostAdiabaticj0Func,
			&j0Params,
			pphi0_lower,
			pphi0_upper,
			1.e-17,
			1.e-20
		);

		pphiVec->data[i] = j0Root.root;
        pphi0Vec->data[i] = j0Root.root;
	}

	*rReverseVec = XLALReverseREAL8Vector(rVec);
    *pphiReverseVec = XLALReverseREAL8Vector(pphiVec);

    *dpphiBydrReverseVec = XLALFDDerivative1Order2(rReverseVec, pphiReverseVec);

    *dpphiBydrVec = XLALReverseREAL8Vector(dpphiBydrReverseVec);
    *dpphiBydr0Vec = XLALReverseREAL8Vector(dpphiBydrReverseVec);

	REAL8 DeltaT;
	REAL8 DeltaR;

	// compute csiVec
	for (i = 0; i < rSize; i++)
	{
		DeltaT = XLALSimIMRSpinEOBHamiltonianDeltaT(seobCoeffs, rVec->data[i], nu, aK);
		DeltaR = XLALSimIMRSpinEOBHamiltonianDeltaR(seobCoeffs, rVec->data[i], nu, aK);

		csiVec->data[i] = sqrt(DeltaT*DeltaR) / (rVec->data[i]*rVec->data[i]+aK*aK);
	}

	// compute omegaVec
	for (i = 0; i < rSize; i++)
	{
		/* polarDynamics contains r, phi, pr, pphi */
		REAL8 polarDynamics[4];
		
		polarDynamics[0] = rVec->data[i];
		polarDynamics[1] = 0.0;
		polarDynamics[2] = prstarVec->data[i];
		polarDynamics[3] = pphiVec->data[i];

		omegaVec->data[i] = XLALSimIMRSpinAlignedEOBCalcOmega(
								polarDynamics,
								seobParams,
								dr
							);
		printf("%.18e\n",omegaVec->data[i]);
	}

    UINT4 n;
    UINT4 parity;

    UNUSED REAL8 dHBydprstar;

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

    		// 	UNUSED REAL8 Flux;
    		// 	Flux = XLALSimInspiralEOBPAFluxWrapper(
						// 	rVec->data[i],
						// 	0.0,
						// 	pphiVec->data[i],
						// 	omegaVec->data[i],
						// 	HVec->data[i],
						// 	seobParams,
						// 	nqcCoeffs,
						// 	LALparams
						// );
    		// 	printf("%.18e\n", Flux);

    			struct PostAdiabaticRoot prstarRoot;

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
    	}

    	if (parity)
		{
		    *prstarReverseVec = XLALReverseREAL8Vector(prstarVec);
		    *dprstarBydrReverseVec = XLALFDDerivative1Order2(rReverseVec, prstarReverseVec);
		    *dprstarBydrVec = XLALReverseREAL8Vector(dprstarBydrReverseVec);
		}
		else
		{
			*pphiReverseVec = XLALReverseREAL8Vector(pphiVec);
			*dpphiBydrReverseVec = XLALFDDerivative1Order2(rReverseVec, pphiReverseVec);
			*dpphiBydrVec = XLALReverseREAL8Vector(dpphiBydrReverseVec);
		}
    }

	exit(0);

	// for (i = 0; i < rSize; i++)
	// {
	// 	r = r0 - i*dr;

	// 	rVec->data[i] = r;

	// 	XLALSimInspiralEOBPostAdiabaticMetricS(&AVec->data[i], &BVec->data[i], &dAVec->data[i], &d2AVec->data[i], &dBVec->data[i], r, LALparams);

	// 	XLALSimInspiralEOBPostAdiabaticGetCentrifugalRadiusLO(&rcVec->data[i], &drcBydrVec->data[i], &d2rBydr2Vec->data[i], r, LALparams);

	// 	REAL8 ggm[14];
	// 	XLALSimInspiralEOBPostAdiabaticsGSDynamics(ggm, r, rcVec->data[i], drcBydrVec->data[i], 0.0, LALparams);

	//     REAL8 G;
	//     G = ggm[2]*S + ggm[3]*Sstar;

 //        dGBydrVec->data[i] = ggm[6]*S + ggm[7]*Sstar;
 //        dGBydprstarVec->data[i] = ggm[4]*S + ggm[5]*Sstar;
 //        dGBydprstarbyprstarVec->data[i] = ggm[10]*S + ggm[11]*Sstar;
        
 //        /* Circular values for flux */
 //        G0Vec->data[i] = G;
 //        dGBydr0Vec->data[i] = dGBydrVec->data[i];

 //        /* Auxiliary variables */
 //        REAL8 uc;

 //        sqrtAByBVec->data[i] = sqrt(AVec->data[i]/BVec->data[i]);
 //        uc = 1. / rcVec->data[i];
 //        uc2Vec->data[i] = uc * uc;
 //        ducBydrVec->data[i] = -uc2Vec->data[i] * drcBydrVec->data[i];
 //        dAuc2BydrVec->data[i] = uc2Vec->data[i] * (dAVec->data[i]-2*AVec->data[i]*uc*drcBydrVec->data[i]);

 //        REAL8 aCoeff;
 //        REAL8 bCoeff;
 //        REAL8 cCoeff;

 //        aCoeff = dAuc2BydrVec->data[i]*dAuc2BydrVec->data[i] - 4*AVec->data[i]*uc2Vec->data[i]*dGBydrVec->data[i]*dGBydrVec->data[i];
 //        bCoeff = 2*dAVec->data[i]*dAuc2BydrVec->data[i] - 4*AVec->data[i]*dGBydrVec->data[i]*dGBydrVec->data[i];
 //        cCoeff = dAVec->data[i] * dAVec->data[i];

 //        REAL8 root1;
 //        REAL8 root2;

 //        gsl_poly_solve_quadratic(aCoeff, bCoeff, cCoeff, &root1, &root2);

 //        REAL8 j02;

 //        // dGdr sign determines choice of solution
 //        if (((dGBydr0Vec->data[i] > 0.) && (aCoeff > 0.)) || ((dGBydr0Vec->data[i] < 0.) && (aCoeff < 0.)))
 //        {
 //            j02 = root2;
 //        }
 //        else
 //        {
 //            j02 = root1;
 //        }

 //        // circular orbital angular momentum
 //        pphiVec->data[i] = sqrt(j02);
 //        pphi0Vec->data[i] = sqrt(j02);

 //        REAL8 H;
 //        REAL8 dHeff_dr;
 //        REAL8 dHeff_dprstar;
 //        REAL8 d2Heff_dprstar20;

 //        XLALSimInspiralEOBPostAdiabaticHamiltonianS(&H, &HeffVec->data[i], &HeffOrbVec->data[i], &dHeff_dr, &dHeff_dprstar, &dHeffBydpphiVec->data[i], &d2Heff_dprstar20, r, rcVec->data[i], drcBydrVec->data[i], pphiVec->data[i], prstarVec->data[i], S, Sstar, AVec->data[i], dAVec->data[i], LALparams);

 //        // printf("%.18e %.18e %.18e\n", rVec->data[i], H, SpinAlignedH);

 //        EVec->data[i] = nu * H;

 //        // Circular orbital frequency 
 //        OmgVec->data[i] = dHeffBydpphiVec->data[i] / EVec->data[i];

 //        /* Circular real orbital frequency */
 //        OmgOrbVec->data[i] = (pphiVec->data[i]*AVec->data[i]*uc2Vec->data[i]) / (EVec->data[i]*HeffOrbVec->data[i]);
 //    }

    // *rReverseVec = XLALReverseREAL8Vector(rVec);
    // *pphiReverseVec = XLALReverseREAL8Vector(pphiVec);

    // *dpphiBydrReverseVec = XLALFDDerivative1Order2(rReverseVec, pphiReverseVec);

    // *dpphiBydrVec = XLALReverseREAL8Vector(dpphiBydrReverseVec);
    // *dpphiBydr0Vec = XLALReverseREAL8Vector(dpphiBydrReverseVec); // Figure out how to copy REAL8Vector here instead of reversing it again

    // UINT4 n;
    // UINT4 parity;

    // UINT4 count = 1;

//     for (n = 1; n <= PA_order; n++)
//     {
//     	// Separating even and odd orders
//     	if (n%2 == 0)
//     	{
//     		parity = 0;
//     	}
//     	else
//     	{
//     		parity = 1;
//     	}

//     	for (i = 0; i < rSize; i++)
//     	{
//     		r = rVec->data[i];

//     		if (parity)
//     		{
//     			// Odd PA orders: corrections to prstar only
//     			REAL8 Heff_orb_f;
//     			REAL8 Heff_f;
//     			REAL8 E_f;
//     			REAL8 psi;
//     			REAL8 r_omg2;
//     			REAL8 v_phi;
//     			REAL8 x;
//     			REAL8 jhat;
//     			REAL8 prstar_fake;
//     			REAL8 ddotr_fake;
//     			REAL8 Fphi;

//     			// Variables for which Kepler's law is still valid
//     			Heff_orb_f = sqrt(AVec->data[i] * (1.+pphiVec->data[i]*pphiVec->data[i]*uc2Vec->data[i]));
//                 Heff_f = G0Vec->data[i]*pphiVec->data[i] + Heff_orb_f;
//                 E_f = sqrt(1 + 2*nu*(Heff_f-1));

//                 psi = (ducBydrVec->data[i]+dGBydr0Vec->data[i]*rcVec->data[i]*sqrt(AVec->data[i]/(pphiVec->data[i]*pphiVec->data[i]) + AVec->data[i]*uc2Vec->data[i])/AVec->data[i]) / (-0.5*dAVec->data[i]);
//                 r_omg2 = pow((((1./sqrt(rcVec->data[i]*rcVec->data[i]*rcVec->data[i]*psi))+G0Vec->data[i])/(E_f)), -2./3.);
//                 v_phi = r_omg2 * OmgVec->data[i];
//                 x = v_phi * v_phi;
//                 jhat = pphiVec->data[i] / (r_omg2*v_phi);

//                 // Add save to REAL8Array for these variables

//                 ddotr_fake  = 0.0;
//                 prstar_fake = 0.0;

//                 XLALSimInspiralEOBPostAdiabaticFluxS(
//                 	&Fphi,
//                 	x,
//                 	OmgVec->data[i],
//                 	r_omg2,
//                 	EVec->data[i],
//                 	HeffVec->data[i],
//                 	jhat,
//                 	r,
//                 	prstar_fake,
//                 	ddotr_fake,
//                 	LALparams);

//                 fluxVec->data[i] = Fphi;

//                 struct PostAdiabaticRoot prstarRoot;

//                 struct PostAdiabaticRootSolveParams prstarParams;
//                 prstarParams.r = rVec->data[i];
// 			    prstarParams.rc = rcVec->data[i];
// 			    prstarParams.drcBydr = drcBydrVec->data[i];
// 			    prstarParams.uc2 = uc2Vec->data[i];
// 			    prstarParams.ducBydr = ducBydrVec->data[i];
// 			    prstarParams.pphi = pphiVec->data[i];
// 			    prstarParams.dpphiBydr = dpphiBydrVec->data[i];
// 			    prstarParams.A = AVec->data[i];
// 			    prstarParams.B = BVec->data[i];
// 			    prstarParams.dA = dAVec->data[i];
// 			    prstarParams.LALParams = LALparams;

// 			    REAL8 x_lower;
// 			    REAL8 x_upper;

// 			    if (prstarVec->data[i] > 0.)
// 			    {
// 					x_lower = 0.9 * prstarVec->data[i];
// 					x_upper = 1.1 * prstarVec->data[i];
// 				}
// 				else
// 				{
// 					x_upper = 0.9 * prstarVec->data[i];
// 					x_lower = 1.1 * prstarVec->data[i];
// 				}

// 			    XLALSimInspiralEOBPostAdiabaticRootFinder(
// 			    	&prstarRoot,
// 					XLALSimInspiralEOBPostAdiabaticdpphiFunc,
// 					&prstarParams,
// 					x_lower,
// 					x_upper,
// 					1.e-17,
// 					1.e-20
// 				);

// 				// printf("%.18e\n", prstarRoot.root);

// 				// printf("%d, %.18e, %.18e, %.18e, %.18e, %.18e, %.18e, %.18e, %.18e\n", count, rVec->data[i], rcVec->data[i], drcBydrVec->data[i], prstarVec->data[i], dprstarBydrVec->data[i], pphiVec->data[i], dpphiBydrVec->data[i], prstarRoot.root);
// 				// count++;

// 				prstarVec->data[i] = prstarRoot.root;
				

//                 // REAL8 dHeff_dprstarbyprstar, dr_dtbyprstar;

// 				// dHeff_dprstarbyprstar = pphiVec->data[i]*dGBydprstarbyprstarVec->data[i] + (1+2*z3*AVec->data[i]*uc2Vec->data[i]*prstarVec->data[i]*prstarVec->data[i])/HeffOrbVec->data[i];
// 	   //          dr_dtbyprstar = sqrtAByBVec->data[i]/(EVec->data[i]) * dHeff_dprstarbyprstar;
// 	   //          prstarVec->data[i] = Fphi / dpphiBydrVec->data[i] / dr_dtbyprstar;

// 	            REAL8 ggm[14];
// 				XLALSimInspiralEOBPostAdiabaticsGSDynamics(ggm, r, rcVec->data[i], drcBydrVec->data[i], prstarVec->data[i], LALparams);

//                 dGBydrVec->data[i] = ggm[6]*S + ggm[7]*Sstar;
//                 dGBydprstarVec->data[i] = ggm[4]*S + ggm[5]*Sstar;
//                 dGBydprstarbyprstarVec->data[i] = ggm[10]*S + ggm[11]*Sstar;
//     		}
//     		else
//     		{
//     			// Even PA orders: corrections to pphi only
//     			struct PostAdiabaticRoot pphiRoot;

//     			struct PostAdiabaticRootSolveParams pphiParams;
//     			pphiParams.r = rVec->data[i];
// 			    pphiParams.rc = rcVec->data[i];
// 			    pphiParams.drcBydr = drcBydrVec->data[i];
//                 pphiParams.dAuc2Bydr = dAuc2BydrVec->data[i];
// 			    pphiParams.dprstarBydr = dprstarBydrVec->data[i];
// 			    pphiParams.dA = dAVec->data[i];
// 			   	pphiParams.prstar = prstarVec->data[i];
// 			   	pphiParams.A = AVec->data[i];
// 			    pphiParams.uc2 = uc2Vec->data[i];

// 			    REAL8 x_lower = 0.9 * pphiVec->data[i];
// 				REAL8 x_upper = 1.1 * pphiVec->data[i];

// 			    XLALSimInspiralEOBPostAdiabaticRootFinder(
// 			    	&pphiRoot,
// 					XLALSimInspiralEOBPostAdiabaticdprstarFunc,
// 					&pphiParams,
// 					x_lower,
// 					x_upper,
// 					1.e-17,
// 					1.e-20
// 				);

// 				pphiVec->data[i] = pphiRoot.root;
// 				// printf("%.18e\n", pphiVec->data[i]);

// 				// printf("%d, %.18e, %.18e, %.18e, %.18e, %.18e, %.18e, %.18e, %.18e\n", count, rVec->data[i], rcVec->data[i], drcBydrVec->data[i], prstarVec->data[i], dprstarBydrVec->data[i], pphiVec->data[i], dpphiBydrVec->data[i], pphiVec->data[i]);
//     // 			count++;
//     		}

//     		REAL8 H;
// 			REAL8 dHeff_dr;
// 			REAL8 dHeff_dprstar;
// 			REAL8 d2Heff_dprstar20;

// 			XLALSimInspiralEOBPostAdiabaticHamiltonianS(&H, &HeffVec->data[i], &HeffOrbVec->data[i], &dHeff_dr, &dHeff_dprstar, &dHeffBydpphiVec->data[i], &d2Heff_dprstar20, r, rcVec->data[i], drcBydrVec->data[i], pphiVec->data[i], prstarVec->data[i], S, Sstar, AVec->data[i], dAVec->data[i], LALparams);

// 	        EVec->data[i] = nu * H;

// 	        /* Circular orbital frequency */
// 	        OmgVec->data[i] = dHeffBydpphiVec->data[i] / EVec->data[i];
	        
// 	        /* Circular real orbital frequency */
// 	        OmgOrbVec->data[i] = (pphiVec->data[i]*AVec->data[i]*uc2Vec->data[i]) / (EVec->data[i]*HeffOrbVec->data[i]);

// 	        dtBydrVec->data[i] = EVec->data[i] / (sqrtAByBVec->data[i]*dHeff_dprstar); // dt_dr = 1/dr_dt 
// 	        dphiBydrVec->data[i] = OmgVec->data[i] * dtBydrVec->data[i];
//     	}

// 		if (parity)
// 		{
// 		    *prstarReverseVec = XLALReverseREAL8Vector(prstarVec);

// 		    *dprstarBydrReverseVec = XLALFDDerivative1Order2(rReverseVec, prstarReverseVec);

// 		    *dprstarBydrVec = XLALReverseREAL8Vector(dprstarBydrReverseVec);
// 		}
// 		else
// 		{
// 			*pphiReverseVec = XLALReverseREAL8Vector(pphiVec);

// 			*dpphiBydrReverseVec = XLALFDDerivative1Order2(rReverseVec, pphiReverseVec);

// 			*dpphiBydrVec = XLALReverseREAL8Vector(dpphiBydrReverseVec);
// 		}
//     }

//     REAL8Vector *dtBydrReverseVec = XLALCreateREAL8Vector(rSize);
//     REAL8Vector *tReverseVec = XLALCreateREAL8Vector(rSize);
//     REAL8Vector *tVec = XLALCreateREAL8Vector(rSize);
//     REAL8Vector *dphiBydrReverseVec = XLALCreateREAL8Vector(rSize);
//     REAL8Vector *phiVec = XLALCreateREAL8Vector(rSize);
//     REAL8Vector *phiReverseVec = XLALCreateREAL8Vector(rSize);

//     memset(dtBydrReverseVec->data, 0, dtBydrReverseVec->length * sizeof(REAL8));
//     memset(tReverseVec->data, 0, tReverseVec->length * sizeof(REAL8));
//     memset(tVec->data, 0, tReverseVec->length * sizeof(REAL8));
//     memset(dphiBydrReverseVec->data, 0, dphiBydrReverseVec->length * sizeof(REAL8));
//     memset(phiVec->data, 0, phiVec->length * sizeof(REAL8));
//     memset(phiReverseVec->data, 0, phiReverseVec->length * sizeof(REAL8));

//     *dtBydrReverseVec = XLALReverseREAL8Vector(dtBydrVec);
//     *tReverseVec = XLALCumulativeIntegral3(rReverseVec, dtBydrReverseVec);
//     *tVec = XLALReverseREAL8Vector(tReverseVec);

//     *tVec = XLALOffsetREAL8Vector(tVec, -tVec->data[0]);

//     *dphiBydrReverseVec = XLALReverseREAL8Vector(dphiBydrVec);
//     *phiReverseVec = XLALCumulativeIntegral3(rReverseVec, dphiBydrReverseVec);
//     *phiVec = XLALReverseREAL8Vector(phiReverseVec);

//     *phiVec = XLALOffsetREAL8Vector(phiVec, -phiVec->data[0]);

//     for (i = 0; i < rSize; i++)
//     {
//         printf("%.18e, %.18e, %.18e\n", rVec->data[i], tVec->data[i], phiVec->data[i]);
//     }

//     exit(0);

//     // Output ready
//     // tVec, rVec, phiVec, pphiVec, prstarVec, pphi0Vec, rcVec, AVec, dpphiBydrVec

    return XLAL_SUCCESS;
}