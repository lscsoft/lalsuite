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

#include "LALSimBlackHoleRingdown.h"

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

REAL8
XLALSimInspiralEOBPACalculateNewtonianj0(
	REAL8 r)
{
	REAL8 Newtonianj0;
	Newtonianj0 = sqrt(r);

	return Newtonianj0;
}