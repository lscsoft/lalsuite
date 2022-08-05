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

#include "LALSimIMRSpinEOBHamiltonian.c"
#include "LALSimIMRSpinEOBHamiltonianOptimized.c"

#include "LALSimBlackHoleRingdown.h"

/**
 * Function which calculates the mass ratio q from the masses m1 and m2
 */
REAL8
XLALSimInspiralEOBPACalculateMassRatio(
	const REAL8 m1,
    /**< Mass of the primary */
	const REAL8 m2
    /**< Mass of the secondary */
)
{
	REAL8 q;

	q = m1 / m2;

	if (q > 1.)
	{
		q = m2 / m1;
	}

	return q;
}

/**
 * Function which calculates the symmetric mass ratio nu from the mass 
 * ratio q
 */
REAL8
XLALSimInspiralEOBPACalculateSymmetricMassRatio(
	const REAL8 q
	/**< Mass ratio */
)
{
	REAL8 nu = 0.;

	if (q > 0.)
	{
		nu = q / ((q+1.)*(q+1.));
	}

	return nu;
}

/**
 * Function which calculates the parameter X1 from the symmetric mass 
 * ratio nu
 */
REAL8
XLALSimInspiralEOBPACalculateX1(
	const REAL8 nu
	/**< Symmetric mass ratio */
)
{
	if ( (nu < 0.) || (nu > 0.25) )
	{
		XLALPrintError("XLAL Error - %s: Symmetric mass ratio is 0 <= nu <= 1/4.\n", __func__);
		XLAL_ERROR(XLAL_EINVAL);
	}

	REAL8 X1 = 0.5 * (1.+sqrt(1.-4.*nu));

	return X1;
}

/**
 * Function which calculates the parameter X2 from the symmetric mass 
 * ratio nu
 */
REAL8
XLALSimInspiralEOBPACalculateX2(
	const REAL8 nu
	/**< Symmetric mass ratio */
)
{
	REAL8 X1 = XLALSimInspiralEOBPACalculateX1(nu);

	REAL8 X2 = 1. - X1;

	return X2;
}

/**
 * Function which calculates the spin parameter a
 */
REAL8
XLALSimInspiralEOBPACalculatea(
	REAL8 X,
    /**< Parameter X for the binary component */
	REAL8 chi
    /**< Spin of the binary component */
)
{
	REAL8 a;
	a = X * X * chi;

	return a;
}

/**
 * Function which calculates the spin parameter S
 */
// REAL8
// XLALSimInspiralEOBPACalculateS(
// 	REAL8 X,
// 	REAL8 chi
// )
// {
// 	REAL8 S;
// 	S = X * X * chi;

// 	return S;
// }

/**
 * Function which calculates the spin parameter Sstar (S*)
 */
REAL8
XLALSimInspiralEOBPACalculateSstar(
	REAL8 X1,
    /**< Parameter X1 */
	REAL8 X2,
    /**< Parameter X2 */
	REAL8 chi1,
    /**< Spin of the primary component */
	REAL8 chi2
    /**< Spin of the secondary component */
)
{
	REAL8 Sstar;
	Sstar = X1 * X2 * (chi1+chi2);

	return Sstar;
}

/**
 * Function which calculates the frequency Omega
 */
REAL8
XLALSimIMRSpinAlignedEOBPACalculateOmega(
    REAL8 polarDynamics[],
    /**< The polar coordinates of a point along the binary inspiral */
    REAL8 dr,
    /**< The spacing of the radial grid */
    SpinEOBParams *seobParams,
    /**< Struct of additional parameters */
    LALDict *LALParams
    /**< Pointer to a dictionary containing additional */
)
{
	const UINT2 analyticFlag = XLALDictLookupUINT2Value(LALParams, "analyticFlag");

	REAL8 omega;

	if (analyticFlag == 0)
	{
		omega = XLALSimIMRSpinAlignedEOBCalcOmega(
			polarDynamics,
			seobParams,
			dr
		);
	}
	else
	{
		omega = XLALSimIMRSpinAlignedEOBCalcOmegaOptimized(
			polarDynamics,
			seobParams
		);
	}

	return omega;
}

/**
 * Function which calculates the time units factor
 */
// double
// XLALSimInspiralEOBPostAdiabaticTimeUnitsFactor(
// 	REAL8 M
// 	/**< Total mass in SI units */)
// {
// 	REAL8 time_units_factor = 1. / (M*LAL_MTSUN_SI);

// 	return time_units_factor;
// }

/**
 * Function which calculates the time units factor
 */
// double
// XLALSimInspiralEOBPostAdiabaticDynr0Kepler(
// 	REAL8 f0
// 	/**< Mass in SI units */)
// {
// 	const REAL8 omg_orb0 = LAL_PI * f0; // = 2 * Pi * (f0/2)

// 	REAL8 r0 = pow(omg_orb0, -2./3.);

// 	return r0;
// }

/**
 * Function which calculates the total spin
 */
// REAL8
// XLALSimInspiralEOBPostAdiabaticTotalSpin(
// 	REAL8 q,
// 	REAL8 a1,
// 	REAL8 a2)
// {
// 	REAL8 invQ;
// 	REAL8 onePlusInvQ;
// 	REAL8 aTotal;

// 	invQ = 1. / q;
// 	onePlusInvQ = 1. + invQ;

// 	aTotal = (a1 + a2*invQ*invQ) / (onePlusInvQ * onePlusInvQ);

// 	return aTotal;
// }

/**
 * Function which calculates the final radius at which the post-adiabatic
 * routine stops
 */
// REAL8
// XLALSimInspiralEOBPostAdiabaticFinalRadius(
// 	REAL8 q,
// 	REAL8 a1,
// 	REAL8 a2)
// {
// 	REAL8 aTotal;

// 	REAL8 rISCO;
// 	REAL8 finalRadiusPrefactor;
// 	REAL8 rFinal;

// 	aTotal = XLALSimInspiralEOBPostAdiabaticTotalSpin(q, a1, a2);
	
// 	rISCO = XLALSimRadiusKerrISCO(aTotal);
// 	finalRadiusPrefactor = 1.0;

// 	rFinal = finalRadiusPrefactor * rISCO;

// 	return rFinal;
// }

/**
 * Function which calculates the final radius at which the post-adiabatic
 * routine stops
 */
REAL8
XLALSimInspiralEOBPostAdiabaticFinalRadiusAlternative(
    REAL8 a
    /**< Spin parameter a */
)
{

  REAL8 rISCO;
  REAL8 finalRadiusPrefactor;
  REAL8 rFinal;


  rISCO = XLALSimRadiusKerrISCO(a);
  finalRadiusPrefactor = 1.0;

  rFinal = finalRadiusPrefactor * rISCO;

  return rFinal;
}

/**
 * Function which calculates the spacing of the radial grid
 */
REAL8
XLALSimInspiralEOBPACalculatedr(
	REAL8 rStart,
    /**< The starting radius */
	REAL8 rFinal,
    /**< The final radius */
	UINT4 rSize
    /**< The number of points along the post-adiabatic inspiral */
)
{
	REAL8 dr;

	dr = (rStart-rFinal) / (rSize-1);

	return dr;
}

/**
 * Function which calculates the circular angular momentum j0
 */
REAL8
XLALSimInspiralEOBPACalculateNewtonianj0(
	REAL8 r
    /**< Value of the radius */
)
{
	REAL8 Newtonianj0;
	Newtonianj0 = sqrt(r);

	return Newtonianj0;
}
