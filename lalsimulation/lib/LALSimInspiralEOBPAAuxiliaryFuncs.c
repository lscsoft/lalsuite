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

REAL8
XLALSimInspiralEOBPACalculateMassRatio(
	const REAL8 m1,
	const REAL8 m2
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

REAL8
XLALSimInspiralEOBPACalculateSymmetricMassRatio(
	const REAL8 q
	/**< mass ratio */
)
{
	REAL8 nu = 0.;

	if (q > 0.)
	{
		nu = q / ((q+1.)*(q+1.));
	}

	return nu;
}

REAL8
XLALSimInspiralEOBPACalculateX1(
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
REAL8
XLALSimInspiralEOBPACalculateX2(
	const REAL8 nu
	/**< symmetric mass ratio */)
{
	REAL8 X1 = XLALSimInspiralEOBPACalculateX1(nu);

	REAL8 X2 = 1. - X1;

	return X2;
}

REAL8
XLALSimInspiralEOBPACalculatea(
	REAL8 X,
	REAL8 chi
)
{
	REAL8 a;
	a = X * X * chi;

	return a;
}

REAL8
XLALSimInspiralEOBPACalculateS(
	REAL8 X,
	REAL8 chi
)
{
	REAL8 S;
	S = X * X * chi;

	return S;
}

REAL8
XLALSimInspiralEOBPACalculateSstar(
	REAL8 X1,
	REAL8 X2,
	REAL8 chi1,
	REAL8 chi2
)
{
	REAL8 Sstar;
	Sstar = X1 * X2 * (chi1+chi2);

	return Sstar;
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
	finalRadiusPrefactor = 1.6;

	rFinal = finalRadiusPrefactor * rISCO;

	return rFinal;
}

REAL8
XLALSimInspiralEOBPACalculatedr(
	REAL8 rStart,
	REAL8 rFinal,
	UINT4 rSize
)
{
	REAL8 dr;

	dr = (rStart-rFinal) / (rSize-1);

	return dr;
}

REAL8
XLALSimInspiralEOBPACalculateNewtonianj0(
	REAL8 r
)
{
	REAL8 Newtonianj0;
	Newtonianj0 = sqrt(r);

	return Newtonianj0;
}