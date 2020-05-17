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

#define KMAX 35 /** Multipolar linear index, max value */
#define PMTERMS_eps 1

static const REAL8 CNlm[35] = {
						    0.17777777777777778, 6.4, 
						    0.0007936507936507937, 0.5079365079365079, 8.678571428571429, 
						    2.2675736961451248e-05, 0.008062484252960444, 1.0414285714285714, 14.447971781305114, 
						    5.010421677088344e-08, 0.0006384836014465644, 0.03106534090909091, 1.9614216236438458, 25.688197074915823, 
						    8.898517797026696e-10, 4.464920289836115e-06, 0.003830520129221428, 0.08778390483440988, 3.584607999290539, 46.98226573426573, 
						    1.0695333890657086e-12, 2.3656367312710969e-07, 4.972309783123969e-05, 0.013643024456211269, 0.21542115380351798, 6.472046810332524, 87.1329124642076, 
						    1.2233225038333268e-14, 9.27700678929842e-10, 4.468579053461083e-06, 0.0002675102834551229, 0.03813282951314997, 0.4894825318738884, 11.627213264559293, 162.79300083906728
};

const INT4 LINDEX[KMAX] = {
					    2,2,
					    3,3,3,
					    4,4,4,4,
					    5,5,5,5,5,
					    6,6,6,6,6,6,
					    7,7,7,7,7,7,7,
					    8,8,8,8,8,8,8,8
};

const INT4 MINDEX[KMAX] = {
					    1,2,
					    1,2,3,
					    1,2,3,4,
					    1,2,3,4,5,
					    1,2,3,4,5,6,
					    1,2,3,4,5,6,7,
					    1,2,3,4,5,6,7,8
};

/* Factorials evaluated for the tail term */
static const REAL8 f14[] = {
							 1.,         1.,          2.,
						     6.,         24.,         120.,
						     720.,       5040.,       40320.,
						     362880.,    3628800.,    39916800.,
						     479001600., 6227020800., 87178291200.
};

/* Include all the static function files we need */
// #include "LALSimIMREOBFactorizedWaveform.c"
// #include "LALSimIMREOBNewtonianMultipole.c"
// #include "LALSimIMREOBNQCCorrection.c"
// #include "LALSimIMRSpinEOBInitialConditions.c"
// #include "LALSimIMRSpinEOBAuxFuncs.c"
// #include "LALSimIMRSpinAlignedEOBHcapDerivative.c"
// #include "LALSimIMRSpinEOBHamiltonian.c"
// #include "LALSimIMRSpinEOBFactorizedWaveform.c"
// #include "LALSimIMRSpinEOBFactorizedFlux.c"
/* OPTIMIZED */
// #include "LALSimIMRSpinEOBHamiltonianOptimized.c"
// #include "LALSimIMRSpinEOBComputeAmpPhasefromEOMSoln.c"
// #include "LALSimIMRSpinAlignedEOBGSLOptimizedInterpolation.c"
// #include "LALSimIMRSpinAlignedEOBHcapDerivativeOptimized.c"
/* END OPTIMIZED */

double
XLALSimInspiralEOBPostAdiabaticSymmetricMassRatio(
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

/*
 * EOB Metric potentials A(r), B(r), and their derivatives, spin version
 */
int
XLALSimInspiralEOBPostAdiabaticMetricS(
	REAL8 *A,
	REAL8 *B,
	REAL8 *dA,
	REAL8 *d2A,
	REAL8 *dB,
	REAL8 r,
	LALDict *LALParams)
{
    UNUSED REAL8 nu = XLALDictLookupREAL8Value(LALParams, "nu");
    INT4 useTidal = XLALDictLookupINT4Value(LALParams, "useTidal");

    REAL8 u = 1. / r;
    REAL8 u2 = u * u;
    REAL8 u3 = u2 * u;
    REAL8 u4 = u2 * u2;
     
    REAL8 rc, drcBydr, d2rcBydr2;
    XLALSimInspiralEOBPostAdiabaticGetCentrifugalRadius(&rc, &drcBydr, &d2rcBydr2, r, LALParams);

    /* A potential and derivative with respect to u */
    REAL8 Aorb, dAorbBydu, d2AorbBydu2;
    XLALSimInspiralEOBPostAdiabaticMetricA5PNlog(&Aorb, &dAorbBydu, &d2AorbBydu2, rc, LALParams);
    
    /* Add here tides if needed */
    if (useTidal)
    {
    	// Include tides code if needed
        ;
    }
    
    /* A potential and derivative with respect to r */  
    REAL8 uc  = 1. / rc;
    REAL8 uc2 = uc * uc;
    REAL8 uc3 = uc2 * uc;
    REAL8 uc4 = uc2 * uc2;
    
    REAL8 dAorb  = -dAorbBydu * uc2;
    REAL8 d2Aorb = 2.*dAorbBydu*uc3 + d2AorbBydu2*uc4;
    
    /* Correct A for spin */
    REAL8 AKerr_Multipole = (1.+2.*uc) / (1.+2.*u);
    REAL8 fss = 1.;
         
    *A = Aorb * AKerr_Multipole * fss;    
    *dA = dAorb*drcBydr*(1.+2.*uc)/(1.+2.*u) - 2.*Aorb*drcBydr*uc2/(1.+2.*u) + 2.*Aorb*(1.+2.*uc)*u2/((1.+2.*u)*(1.+2.*u));
    *d2A = d2Aorb*(1.+2.*uc)/(1.+2.*u) + 4.*dAorb*( u2*(1.+2.*uc)/((1.+2.*u)*(1.+2.*u)) - uc2/(1.+2.*u)*drcBydr) + Aorb*(-4.*u3*(1.+2.*uc)/((1.+2.*u)*(1.+2.*u)) + 8.*u4*(1.+2.*uc)/((1.+2.*u)*(1.+2.*u)*(1.+2.*u))+4.*uc3*(1.+2.*u)*drcBydr*drcBydr - 2.*uc2/(1.+2.*u)*d2rcBydr2);
    
    // /* D potential and derivative with respect to r */
    REAL8 Dp = 1.0 + 6.*nu*uc2 - 2.*(3.0*nu-26.0)*nu*uc3;
    REAL8 D  = 1./Dp;
    REAL8 dD = 6.*uc2*(2.*nu*uc-(3.*nu-26.)*nu*uc2)*D*D;
    
    // /* B potential and derivative with respect to r */
    *B = r * r * uc2 * D / (*A);
    *dB = (dD*(*A) - D*(*dA))/((*A)*(*A));

    return XLAL_SUCCESS;
}

int
XLALSimInspiralEOBPostAdiabaticGetCentrifugalRadius(
	REAL8 *rc,
	REAL8 *drcBydr,
	REAL8 *d2rcBydr2,
	REAL8 r,
	LALDict *LALParams)
{
	const char *centrifugalRadius = XLALDictLookupStringValue(LALParams, "centrifugalRadius");
	// const char centrifugalRadius[] = "LO";

	if (strcmp(centrifugalRadius, "LO") == 0)
	{
		XLALSimInspiralEOBPostAdiabaticGetCentrifugalRadiusLO(rc, drcBydr, d2rcBydr2, r, LALParams);
	}
	else if (strcmp(centrifugalRadius, "NLO") == 0)
	{
		XLALPrintError("The option NLO is not available yet.\n");
	}
	else if (strcmp(centrifugalRadius, "NNLO") == 0)
	{
		XLALPrintError("The option NNLO is not available yet.\n");
	}
	else if (strcmp(centrifugalRadius, "NNLOS4") == 0)
	{
		XLALPrintError("The option NNLOS4 is not available yet.\n");
	}
	else if (strcmp(centrifugalRadius, "NOSPIN") == 0)
	{
		XLALPrintError("The option NOSPIN is not available yet.\n");
	}
	else if (strcmp(centrifugalRadius, "NOTIDES") == 0)
	{
		XLALPrintError("The option NOTIDES is not available yet.\n");
	}
	else
	{
		XLALPrintError("Wrong option for centrifugalRadius parameter, use one of \"LO\", \"NLO\", \"NNLO\", \"NNLOS4\", \"NOSPIN\", and \"NOTIDES\".\n");
	}

	return XLAL_SUCCESS;
}

int
XLALSimInspiralEOBPostAdiabaticGetCentrifugalRadiusLO(
	REAL8 *rc,
	REAL8 *drcBydr,
	REAL8 *d2rcBydr2,
	REAL8 r,
	LALDict *LALParams)
{
	REAL8 u = 1. / r;
	REAL8 u3 = u * u * u;

	const INT4 useTidal = XLALDictLookupINT4Value(LALParams, "useTidal");

	if (useTidal)
	{
		// Include tides code here
		;
	}
	else
	{
		REAL8 nu = XLALDictLookupREAL8Value(LALParams, "nu");
		REAL8 a1 = XLALDictLookupREAL8Value(LALParams, "a1");
		REAL8 a2 = XLALDictLookupREAL8Value(LALParams, "a2");
		REAL8 aK2 = XLALDictLookupREAL8Value(LALParams, "aK2");

		REAL8 X12 = sqrt(1. - 4.*nu);   
	    REAL8 cssNLO = (- a2*a2*(1.25 + 1.25*X12 + 0.5*nu) - a1*a1*(1.25 - 1.25*X12 + 0.5*nu) + a1*a2*(-2.+nu));
	    REAL8 rc2 = r*r + aK2*(1. + 2.*u) + u*cssNLO;
	    REAL8 divrc = 1.0/sqrt(rc2);

	    *rc = sqrt(rc2);
	    *drcBydr = r * divrc * (1-(aK2 + 0.5*cssNLO)*u3);	
	    *d2rcBydr2 = divrc * (1.-(*drcBydr)*r*divrc*(1.-(aK2+0.5*cssNLO)*u3)+(2.*aK2 + cssNLO)*u3);
	}

	return XLAL_SUCCESS;
}

int
XLALSimInspiralEOBPostAdiabaticMetricA5PNlog(
	REAL8 *Aorb,
	REAL8 *dAorbBydu,
	REAL8 *d2AorbBydu2,
	REAL8 r,
	LALDict *LALParams)
{
	REAL8 nu = XLALDictLookupREAL8Value(LALParams, "nu");

	REAL8 nu2 = nu * nu;
    const REAL8 PI2 = LAL_PI * LAL_PI;
    const REAL8 PI4 = PI2 * PI2;
    REAL8 u = 1. / r;
    REAL8 u2 = u * u;
    REAL8 u3 = u2 * u;
    REAL8 u4 = u2 * u2;
    REAL8 u5 = u4 * u;
    // double u6   = u5*u;
    // double u7   = u6*u;
    // double u10  = u5*u5;
    // double u8   = u5*u3;
    // double u9   = u8*u;
    REAL8 logu = log(u);

    REAL8 a5c0 = -4237./60. + 2275./512.*PI2 + 256./5.*LAL_LN2 + 128./5.*LAL_GAMMA;
    REAL8 a5c1 = -221./6. + 41./32.*PI2;
    REAL8 a5 = a5c0 + nu*a5c1;
    REAL8 a6 = 3097.3*nu2 - 1330.6*nu + 81.38;
  
    /* 4PN and 5PN coefficients including all known log terms */
    REAL8 a5tot = a5 + 64./5.*logu;
    REAL8 a6tot = a6 + (-7004./105. - 144./5.*nu)*logu;
    REAL8 a5tot2 = a5tot * a5tot;
  
    /* Coefficients of the Padeed function */
    REAL8 N1 = (-3*(-512 - 32*nu2 + nu*(3520 + 32*a5tot + 8*a6tot - 123*PI2))) / (-768 + nu*(3584 + 24*a5tot - 123*PI2));
    REAL8 D1 = (nu*(-3392 - 48*a5tot - 24*a6tot + 96*nu + 123*PI2)) / (-768 + nu*(3584 + 24*a5tot - 123*PI2));
    REAL8 D2 = (2*nu*(-3392 - 48*a5tot - 24*a6tot + 96*nu + 123*PI2)) / (-768 + nu*(3584 + 24*a5tot - 123*PI2));
    REAL8 D3 = (-2*nu*(6016 + 48*a6tot + 3392*nu + 24*a5tot*(4 + nu) - 246*PI2 - 123*nu*PI2)) / (-768 + nu*(3584 + 24*a5tot - 123*PI2));
    REAL8 D4 = -(nu*(-4608*a6tot*(-4 + nu) + a5tot*(36864 + nu*(72192 - 2952*PI2)) + nu*(2048*(5582 + 9*nu) - 834432*PI2 + 15129*PI4)))/(96.*(-768 + nu*(3584 + 24*a5tot - 123*PI2)));
    REAL8 D5 = (nu*(-24*a6tot*(1536 + nu*(-3776 + 123*PI2)) + nu*(-2304*a5tot2 + 96*a5tot*(-3392 + 123*PI2) - (-3776 + 123*PI2)*(-3008 - 96*nu + 123*PI2))))/(96.*(-768 + nu*(3584 + 24*a5tot - 123*PI2)));
  
    /* First derivatives */
    REAL8 dN1 = (160*nu*(-828672 - 32256*nu2 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*PI2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 174045*PI2)))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*PI2)),2)*u);
    REAL8 dD1 = (160*nu*(-828672 - 32256*nu2 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*PI2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 174045*PI2)))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*PI2)),2)*u);
    REAL8 dD2 = (320*nu*(-828672 - 32256*nu2 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*PI2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 174045*PI2)))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*PI2)),2)*u);
    REAL8 dD3 = (640*nu*(-828672 - 32256*nu2 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*PI2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 174045*PI2)))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*PI2)),2)*u);
    REAL8 dD4 = (-320*(-4 + nu)*nu*(-828672 - 32256*nu2 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*PI2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 174045*PI2)))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*PI2)),2)*u);
    REAL8 dD5 = (nu*(-8400*nu*(-24*(a6 - (4*logu*(1751 + 756*nu))/105.)*(1536 + nu*(-3776 + 123*PI2)) + nu*(-2304*gsl_pow_int(a5 + (64*logu)/5.,2) + 96*(a5 + (64*logu)/5.)*(-3392 + 123*PI2) - (-3776 + 123*PI2)*(-32*(94 + 3*nu) + 123*PI2))) - (1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*PI2)))*(4128768*logu*nu + 5*(-2689536 + nu*(11170624 + 64512*a5 - 380685*PI2) - 756*nu*(1536 + nu*(-3776 + 123*PI2))))))/(2625.*gsl_pow_int(-768 + nu*(3584 + 24*(a5 + (64*logu)/5.) - 123*PI2),2)*u);
    
    /* Numerator and denominator of the Pade */
    REAL8 Num = 1 + N1*u;
    REAL8 Den = 1 + D1*u + D2*u2 + D3*u3 + D4*u4 + D5*u5;
    *Aorb = Num/Den;
    
    /* First derivative */
    REAL8 dNum  = dN1*u + N1;
    REAL8 dDen  = D1 + u*(dD1 + 2*D2) + u2*(dD2 + 3*D3) + u3*(dD3 + 4*D4) + u4*(dD4 + 5*D5) + dD5*u5;
  
    /* Derivative of A function with respect to u */
    REAL8 prefactor = (*Aorb)/(Num*Den);
    REAL8 dABydu = prefactor * (dNum*Den - dDen*Num);

    /* Derivative of A with respect to r */
    /* *dA = -u2*dA_u; */
    *dAorbBydu = dABydu;

    /* Second derivatives of Pade coefficients */
    REAL8 d2N1 = (160*nu*(-3840 + 1536*logu*nu + nu*(20992 + 120*a5 - 615*PI2))*(828672 + nu*(-42024*a5 - 8064*a6 + 3584*(-1397 + 9*nu) + 174045*PI2) + 756*nu*(768 + nu*(-3584 - 24*a5 + 123*PI2))))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*PI2)),3)*u2);
    REAL8 d2D1 = (160*nu*(-3840 + 1536*logu*nu + nu*(20992 + 120*a5 - 615*PI2))*(828672 + nu*(-42024*a5 - 8064*a6 + 3584*(-1397 + 9*nu) + 174045*PI2) + 756*nu*(768 + nu*(-3584 - 24*a5 + 123*PI2))))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*PI2)),3)*u2);
    REAL8 d2D2 = (320*nu*(-3840 + 1536*logu*nu + nu*(20992 + 120*a5 - 615*PI2))*(828672 + nu*(-42024*a5 - 8064*a6 + 3584*(-1397 + 9*nu) + 174045*PI2) + 756*nu*(768 + nu*(-3584 - 24*a5 + 123*PI2))))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*PI2)),3)*u2);
    REAL8 d2D3 = (640*nu*(-3840 + 1536*logu*nu + nu*(20992 + 120*a5 - 615*PI2))*(828672 + nu*(-42024*a5 - 8064*a6 + 3584*(-1397 + 9*nu) + 174045*PI2) + 756*nu*(768 + nu*(-3584 - 24*a5 + 123*PI2))))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*PI2)),3)*u2);
    REAL8 d2D4 = (320*(-4 + nu)*nu*(-828672 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*PI2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 32256*nu - 174045*PI2))*(-3840 + 1536*logu*nu + nu*(20992 + 120*a5 - 615*PI2)))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*PI2)),3)*u2);
    REAL8 d2D5 = (nu*(gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*PI2)),2)*(4128768*logu*nu - 7680*(1751 + 756*nu) + nu*(64*(808193 + 5040*a5 + 223020*nu) - 615*(3095 + 756*nu)*PI2)) + 3072*nu*(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*PI2)))*(4128768*logu*nu - 7680*(1751 + 756*nu) + 5*nu*(64*(174541 + 1008*a5 + 44604*nu) - 123*(3095 + 756*nu)*PI2)) + 25804800*nu2*(-24*(a6 - (4*logu*(1751 + 756*nu))/105.)*(1536 + nu*(-3776 + 123*PI2)) + nu*(-2304*gsl_pow_int(a5 + (64*logu)/5.,2) + 96*(a5 + (64*logu)/5.)*(-3392 + 123*PI2) - (-3776 + 123*PI2)*(-32*(94 + 3*nu) + 123*PI2))) + 42000*nu*(-768 + nu*(3584 + 24*(a5 + (64*logu)/5.) - 123*PI2))*(-24*(a6 - (4*logu*(1751 + 756*nu))/105.)*(1536 + nu*(-3776 + 123*PI2)) + nu*(-2304*gsl_pow_int(a5 + (64*logu)/5.,2) + 96*(a5 + (64*logu)/5.)*(-3392 + 123*PI2) - (-3776 + 123*PI2)*(-32*(94 + 3*nu) + 123*PI2)))))/(13125.*gsl_pow_int(-768 + nu*(3584 + 24*(a5 + (64*logu)/5.) - 123*PI2),3)*u2);
    
    /* Second derivative of numerator and denominator */
    REAL8 d2Num = 2.*dN1 + d2N1*u;
    REAL8 d2Den = 2.*(D2 + dD1) + u*(6.*D3 + 4.*dD2 + d2D1) + u2*(12.*D4 + 6.*dD3 + d2D2) + u3*(20.*D5 + 8.*dD4 + d2D3) + u4*(10.*dD5 + d2D4) + u5*d2D5;

    /* Second derivative with respect of u */
    REAL8 d2A_u = prefactor*(2.*dDen*dDen*(*Aorb) - 2.*dNum*dDen + Den*d2Num - d2Den*Num);

    *d2AorbBydu2 = d2A_u;

	return XLAL_SUCCESS;
}

int
XLALSimInspiralEOBPostAdiabaticsGSDynamics(
	REAL8 *ggm,
	REAL8 r,
	REAL8 rc,
	REAL8 drcBydr,
	REAL8 prstar,
	LALDict *LALParams)
{
	static REAL8 c10;
	static REAL8 c20;
	static REAL8 c30;
	static REAL8 c02;
	static REAL8 c12;
	static REAL8 c04;

    static REAL8 cs10;
    static REAL8 cs20;
    static REAL8 cs30;
    static REAL8 cs40;
    static REAL8 cs02;
    static REAL8 cs12;
    static REAL8 cs04;

    /* Compute the nu-dependent coefficient at first call only */
    static INT2 firstcall = 1;

    if (firstcall)
    {
        firstcall = 0;

        REAL8 nu = XLALDictLookupREAL8Value(LALParams, "nu");
        REAL8 cN3LO = XLALDictLookupREAL8Value(LALParams, "c3");

        REAL8 nu2 = nu*nu;

        /* coefficients of hat{GS} */
        c10 = 5./16. * nu;
        c20 = 51./8.*nu + 41./256.*nu2;
        c30 = nu * cN3LO;
        c02 = 27. / 16. * nu;
        c12 = 12.*nu - 49./128.*nu2;
        c04 = -5./16.*nu + 169./256.*nu2;

        /* coefficients of hat{GS*} */
        cs10 = 3./4. + nu/2.;
        cs20 = 27./16. + 29./4.*nu + 3./8.*nu2;
        cs02 = 5./4. + 3./2.*nu;
        cs12 = 4. + 11.*nu - 7./8.*nu2;
        cs04 = 5./48. + 25./12.*nu + 3./8.*nu2;
        cs30 = nu*cN3LO + 135./32.;
        cs40 = 2835./256.;
    }

    REAL8 u = 1. / r;
    REAL8 u2 = u * u;
  
    REAL8 uc = 1. / rc;
    REAL8 uc2 = uc * uc;
    REAL8 uc3 = uc2 * uc;
    REAL8 uc4 = uc3 * uc;
    REAL8 prstar2 = prstar * prstar;
    REAL8 prstar4 = prstar2 * prstar2;
  
    REAL8 GS0 = 2. * u * uc2;
    REAL8 dGS0_duc = 2.*u2/drcBydr + 4.*u*uc;
  
    REAL8 GSs0 = 3./2.*uc3;
    REAL8 dGSs0_duc = 9./2.*uc2;
    REAL8 dGSs0_dprstar = 0.0;
    REAL8 dGSs0_dpph = 0.0;
  
    REAL8 hGS = 1./(1. + c10*uc + c20*uc2 + c30*uc3 + c02*prstar2 + c12*uc*prstar2 + c04*prstar4);   
    REAL8 hGSs = 1./(1. + cs10*uc + cs20*uc2 + cs30*uc3 + cs40*uc4 + cs02*prstar2 + cs12*uc*prstar2 + cs04*prstar4); 
  
    /* complete gyro-gravitomagnetic functions */
    REAL8 GS = GS0*hGS; 
    REAL8 GSs = GSs0*hGSs; 
  
    /* Get derivatives of gyro-gravitomagnetic functions */
    REAL8 dhGS_dprstar = -2.*prstar*hGS*hGS *( c02 + c12*uc + 2.*c04*prstar2);
    REAL8 dhGSs_dprstar = -2.*prstar*hGSs*hGSs*(cs02 + cs12*uc + 2.*cs04*prstar2);
  
    REAL8 dGS_dprstar = GS0 * dhGS_dprstar; 
    REAL8 dGSs_dprstar = GSs0*dhGSs_dprstar + dGSs0_dprstar*hGSs; 
  
    /* derivatives of hat{G} with respect to uc */
    REAL8 dhGS_duc = -hGS * hGS * (c10 + 2.*c20*uc  + 3.*c30*uc2);
    REAL8 dhGSs_duc = -hGSs * hGSs * (cs10 + 2.*cs20*uc + 3.*cs30*uc2 + 4.*cs40*uc3);
  
    /* derivatives of G with respect to uc */
    REAL8 dGS_duc  =  dGS0_duc*hGS  +  GS0*dhGS_duc;
    REAL8 dGSs_duc = dGSs0_duc*hGSs + GSs0*dhGSs_duc;
  
    /* derivatives of (G,G*) with respect to r */
    REAL8 dGS_dr = -drcBydr * uc2 * dGS_duc; 
    REAL8 dGSs_dr = -drcBydr * uc2 * dGSs_duc; 
  
    /* derivatives of (G,G*) with respect to pph */
    REAL8 dGS_dpph = 0.; 
    REAL8 dGSs_dpph = dGSs0_dpph * hGSs;    
  
    /* For initial data: compute the two ratios of ggm.dG_dprstar/prstar for GS and GSs */
    const REAL8 dGS_dprstarbyprstar = -2. * GS0 * hGS * hGS * (c02 + c12*uc + 2.*c04*prstar2);
    const REAL8 dGSs_dprstarbyprstar = -2. * GSs0 * hGSs * hGSs * (cs02 + cs12*uc + 2.*cs04*prstar2);
  
    /* For NQC: Second derivatives neglecting all pr_star^2 terms */
    const REAL8 d2GS_dprstar20 = GS0 * (-2.*hGS*hGS * (c02 +  c12*uc +  2.*c04*prstar2));
    const REAL8 d2GSs_dprstar20 = GSs0 * (-2.*hGSs*hGSs * (cs02 + cs12*uc + 2.*cs04*prstar2));
  
    ggm[0] = hGS;
    ggm[1] = hGSs;
    ggm[2] = GS;
    ggm[3] = GSs;
    ggm[4] = dGS_dprstar;
    ggm[5] = dGSs_dprstar;
    ggm[6] = dGS_dr;
    ggm[7] = dGSs_dr;
    ggm[8] = dGS_dpph;
    ggm[9] = dGSs_dpph;
    ggm[10] = dGS_dprstarbyprstar;
    ggm[11] = dGSs_dprstarbyprstar;
    ggm[12] = d2GS_dprstar20;
    ggm[13] = d2GSs_dprstar20;

	return XLAL_SUCCESS;
}

// remove S and Sstar from the parameters
int
XLALSimInspiralEOBPostAdiabaticHamiltonianS(
	UNUSED REAL8 *H,
	UNUSED REAL8 *Heff,
	UNUSED REAL8 *Heff_orb,
	UNUSED REAL8 *dHeff_dr,
	UNUSED REAL8 *dHeff_dprstar,
	UNUSED REAL8 *dHeff_dpphi,
	UNUSED REAL8 *d2Heff_dprstar20,
	UNUSED REAL8 r,
	UNUSED REAL8 rc,
	UNUSED REAL8 drc_dr,
	UNUSED REAL8 pphi,
	UNUSED REAL8 prstar,
	UNUSED REAL8 S,
	UNUSED REAL8 Sstar,
	UNUSED REAL8 A,
	UNUSED REAL8 dA,
	UNUSED LALDict *LALParams)
{
	REAL8 nu;
	nu = XLALDictLookupREAL8Value(LALParams, "nu");

	/* Shorthands */
    const REAL8 z3 = XLALSimInspiralEOBPostAdiabaticz3(nu);
    const REAL8 pphi2 = pphi * pphi;
    const REAL8 prstar2 = prstar * prstar;
    const REAL8 prstar4 = prstar2 * prstar2;
    const REAL8 uc = 1. / rc;
    const REAL8 uc2 = uc * uc;
    const REAL8 uc3 = uc2 * uc;
    
    /* Compute spin-related functions */
    REAL8 ggm[14];
    XLALSimInspiralEOBPostAdiabaticsGSDynamics(ggm, r, rc, drc_dr, prstar, LALParams);

    const REAL8 GS = ggm[2];
    const REAL8 GSs = ggm[3];
    const REAL8 dGS_dprstar = ggm[4];
    const REAL8 dGSs_dprstar = ggm[5];
    const REAL8 dGS_dr = ggm[6];
    const REAL8 dGSs_dr = ggm[7];
    const REAL8 dGSs_dpphi = ggm[9];
    const REAL8 d2GS_dprstar20 = ggm[12];
    const REAL8 d2GSs_dprstar20 = ggm[13];
    
    /* Compute Hamiltonian and its derivatives */
    *Heff_orb = sqrt(prstar2 + A*(1.+pphi2*uc2+z3*prstar4*uc2));
    *Heff = *Heff_orb + (GS*S + GSs*Sstar)*pphi;
    *H = sqrt(1.+2.*nu*(*Heff-1.)) / nu;
    *dHeff_dr = pphi*(dGS_dr*S + dGSs_dr*Sstar) + 1./(2.*(*Heff_orb))*( dA*(1. + pphi2*uc2 + z3*prstar4*uc2) - 2.*A*uc3*drc_dr*(pphi2 + z3*prstar4) );
    *dHeff_dprstar = pphi*(dGS_dprstar*S + dGSs_dprstar*Sstar) + (prstar/(*Heff_orb))*(1. + 2.*A*uc2*z3*prstar2);
    *d2Heff_dprstar20 = pphi*(d2GS_dprstar20*S + d2GSs_dprstar20*Sstar) +  (1./(*Heff_orb))*(1. + 2.*A*uc2*z3*prstar2);
    *dHeff_dpphi = GS*S + (GSs + pphi*dGSs_dpphi)*Sstar + pphi*A*uc2/(*Heff_orb);

	return XLAL_SUCCESS;
}

int
XLALSimInspiralEOBPostAdiabaticFluxS(
	UNUSED REAL8 *Flux,
	UNUSED REAL8 x,
	UNUSED REAL8 Omega,
	UNUSED REAL8 r_omega,
	UNUSED REAL8 E,
	UNUSED REAL8 Heff,
	UNUSED REAL8 jhat,
	UNUSED REAL8 r,
	UNUSED REAL8 pr_star,
	UNUSED REAL8 ddotr,
	UNUSED LALDict *LALParams)
{
	const INT4 useSpins = XLALDictLookupINT4Value(LALParams, "useSpins");
	const INT4 useTidal = XLALDictLookupINT4Value(LALParams, "useTidal");
	const REAL8 nu = XLALDictLookupREAL8Value(LALParams, "nu");
 //    const double chi1 = dyn->chi1;
 //    const double chi2 = dyn->chi2;
 //    const double X1 = dyn->X1;
 //    const double X2 = dyn->X2;
 //    const double a1 = dyn->a1;
 //    const double a2 = dyn->a2;
 //    const double C_Q1 = dyn->C_Q1;
 //    const double C_Q2 = dyn->C_Q2;
 //    const double X12  = X1-X2; /* sqrt(1-4nu) */
 //    const double X12sq = SQ(X12); /* (1-4nu) */

 //    const int usetidal = dyn->use_tidal;
 //    const int usespins = dyn->use_spins;

	REAL8 prefact[] = {
			        jhat, Heff,
			        Heff, jhat, Heff,
			        jhat, Heff, jhat, Heff,
			        Heff, jhat, Heff, jhat, Heff,
			        jhat, Heff, jhat, Heff, jhat, Heff,
			        Heff, jhat, Heff, jhat, Heff, jhat, Heff,
			        jhat, Heff, jhat, Heff, jhat, Heff, jhat, Heff};

    REAL8 FNewtlm[KMAX];
    REAL8 MTlm[KMAX];
    REAL8 rholm[KMAX];
	REAL8 flm[KMAX];
	REAL8 FNewt22;
	INT4 k;
	REAL8 hlmNQC[KMAX];
	REAL8 Modhhatlm[KMAX];
	REAL8 sum_k;
	REAL8 hatf;

	/* Newtonian flux */
    XLALSimInspiralEOBPostAdiabaticFlmNewtonian(FNewtlm, x, LALParams);

    /* Correct amplitudes for specific multipoles and cases */
    if (useSpins)
    {
        /* Correct (2,1), (3,1) and (3,3) (sp2 = 1) */
        REAL8 x6 = gsl_pow_int(x, 6);

        FNewtlm[0] = CNlm[0] * x6; /* (2,1) */
        FNewtlm[2] = CNlm[2] * x6; /* (3,1) */
        FNewtlm[4] = CNlm[4] * x6; /* (3,3) */

        /* Correct (4,1), (4,3)  (sp4 = (1-2nu)^2) */
        REAL8 sp4x8 = (1-2*nu) * (1-2*nu) * gsl_pow_int(x, 8);

        FNewtlm[5] = CNlm[5] * sp4x8; /* (4,1) */
        FNewtlm[7] = CNlm[7] * sp4x8; /* (4,3) */
    }
    else
    {
        if (useTidal) {
            /* Correct (2,1), (3,1) and (3,3) (sp2 = 1) */
            double x6 = gsl_pow_int(x, 6);

            FNewtlm[0] = CNlm[0] * x6; /* (2,1) */
            FNewtlm[2] = CNlm[2] * x6; /* (3,1) */
            FNewtlm[4] = CNlm[4] * x6; /* (3,3) */
        }
    }

    /* Tail term */
    XLALSimInspiralEOBPostAdiabaticTlmFlux(MTlm, E*Omega);

    /* Amplitudes */
    if (useSpins)
    {
        XLALSimInspiralEOBPostAdiabaticGetWavFlmS(rholm, flm, x, LALParams);
    }
    else
    {
        XLALSimInspiralEOBPostAdiabaticWavFlm(rholm, flm, x, LALParams);
    }

    FNewt22 = FNewtlm[1];

    /* NQC correction to the modulus of the (l,m) waveform */
    for (k = 0; k < KMAX; k++)
    {
    	/* no NQC */
    	hlmNQC[k] = 1.; 
    }

    /* Compute modulus of hhat_lm (with NQC) */  
    for (k = 0; k < KMAX; k++)
    { 
        Modhhatlm[k] = prefact[k] * MTlm[k] * flm[k] * hlmNQC[k];
    }

    if (useTidal)
    {
        // Include tides code here if needed
        ;
    }

    /* Total multipolar flux */
    sum_k = 0.;

    for (k = KMAX; k--;)
    {
    	sum_k += Modhhatlm[k] * Modhhatlm[k] * FNewtlm[k];
    }

    /* Normalize to the 22 Newtonian multipole */
    hatf = sum_k / FNewt22;

    /* Horizon flux */ 
    if (!(useTidal))
    {
        REAL8 hatFH;
        
        if (useSpins)
        {
            XLALSimInspiralEOBPostAdiabaticHorizonFluxS(&hatFH, x, Heff, jhat, LALParams);
        }
        else
        {
            XLALSimInspiralEOBPostAdiabaticHorizonFlux(&hatFH, x, Heff, jhat, LALParams);
        }
        
        hatf += hatFH;
    }

    *Flux = (-32./5. * nu * gsl_pow_int(Omega*r_omega, 4) * Omega * hatf);

	return XLAL_SUCCESS;
}

/* Newtonian partial fluxes */
int
XLALSimInspiralEOBPostAdiabaticFlmNewtonian(
	UNUSED REAL8 *Nlm,
	UNUSED REAL8 x,
	UNUSED LALDict *LALParams)
{
	const REAL8 nu = XLALDictLookupREAL8Value(LALParams, "nu");

    /* Shorthands */
    const REAL8 nu2 = nu * nu;
    const REAL8 nu3 = nu2 * nu;
    const REAL8 x5  = x * x * x * x * x;
    const REAL8 x6  = x5 * x;
    const REAL8 x7  = x6 * x;
    const REAL8 x8  = x7 * x;
    const REAL8 x9  = x8 * x;
    const REAL8 x10 = x9 * x;
    const REAL8 x11 = x10 * x;
    const REAL8 x12 = x11 * x;
  
    const REAL8 sp2 = 1. - 4.*nu;
    const REAL8 sp4 = (1-4*nu) * (1-2*nu) * (1-2*nu);
    const REAL8 sp3 = (1.-3.*nu) * (1.-3.*nu);
    const REAL8 sp5 = (1.-5.*nu+5.*nu2) * (1.-5.*nu+5.*nu2);
    const REAL8 sp6 = (1-4*nu) * (3*nu2-4*nu +1) * (3*nu2-4*nu +1);
    const REAL8 sp7 = (1-7*nu+14*nu2-7*nu3) * (1-7*nu+14*nu2-7*nu3);
    const REAL8 sp8 = (1-4*nu) * (1-6*nu+10*nu2-4*nu3) * (1-6*nu+10*nu2-4*nu3);

    REAL8 spx[] = {
		        sp2 * x6, x5, 
		        sp2 * x6, sp3 * x7, sp2 * x6, 
		        sp4 * x8, sp3 * x7, sp4 * x8, sp3 * x7, 
		        sp4 * x8, sp5 * x9, sp4 * x8, sp5 * x9, sp4 * x8, 
		        sp6 * x10, sp5 * x9, sp6 * x10, sp5 * x9, sp6 * x10, sp5 * x9, 
		        sp6 * x10, sp7 * x11, sp6 * x10, sp7 * x11, sp6 * x10, sp7 * x11, sp6 * x10,
		        sp8 * x12, sp7 * x11, sp8 * x12, sp7 * x11, sp8 * x12, sp7 * x11, sp8 * x12, (7*nu3-14*nu2+7*nu-1)*(7*nu3-14*nu2+7*nu-1) * x11
    };

    /* Newtonian partial fluxes */
    INT4 k;

    for (k = 0; k < KMAX; k++)
    {
        Nlm[k] = CNlm[k] * spx[k];
    }

    return XLAL_SUCCESS;
}

/** Tail term (modulus) */
int
XLALSimInspiralEOBPostAdiabaticTlmFlux(
	UNUSED REAL8 *MTlm,
	UNUSED REAL8 w)
{
	INT4 k;
    REAL8 hhatk;
    REAL8 x2;
    REAL8 y;
    REAL8 prod;
    INT4 j;

    for (k = 0; k < KMAX; k++)
    {
        hhatk = MINDEX[k] * w;
        x2 = 4. * hhatk * hhatk;
        prod = 1.;
        
        for (j=1; j <= LINDEX[k]; j++)
        {
            prod *= (j*j + x2);
        }

        y = 4. * LAL_PI * hhatk;
        y /= (1. - exp(-y)); 
        
        MTlm[k] = sqrt(1. / (f14[LINDEX[k]]*f14[LINDEX[k]]) * y * prod);
    }

    return XLAL_SUCCESS;
}

int
XLALSimInspiralEOBPostAdiabaticGetWavFlmS(
	UNUSED REAL8 *rholm,
	UNUSED REAL8 *flm,
	UNUSED REAL8 x,
	UNUSED LALDict *LALParams)
{
	const char *useFlm = XLALDictLookupStringValue(LALParams, "useFlm");

	if (strcmp(useFlm, "SSLO") == 0)
	{
		XLALSimInspiralEOBPostAdiabaticGetWavFlmSSLO(rholm, flm, x, LALParams);
	}
	else if (strcmp(useFlm, "SSNLO") == 0)
	{
		XLALPrintError("The option SSNLO is not available yet.\n");
	}
	else
	{
		XLALPrintError("Wrong option for useFlm parameter, use one of \"SSLO\" and \"SSNLO\".\n");
	}

	return XLAL_SUCCESS;
}

int
XLALSimInspiralEOBPostAdiabaticGetWavFlmSSLO(
	UNUSED REAL8 *rholm,
	UNUSED REAL8 *flm,
	UNUSED REAL8 x,
	UNUSED LALDict *LALParams)
{
	const REAL8 a1 = XLALDictLookupREAL8Value(LALParams, "a1");
	const REAL8 a2 = XLALDictLookupREAL8Value(LALParams, "a2");
	const REAL8 X1 = XLALDictLookupREAL8Value(LALParams, "X1");
	const REAL8 X2 = XLALDictLookupREAL8Value(LALParams, "X2");
	const INT4 useTidal = XLALDictLookupINT4Value(LALParams, "useTidal");
	const REAL8 nu = XLALDictLookupREAL8Value(LALParams, "nu");

	/** Orbital part */
    XLALSimInspiralEOBPostAdiabaticWavFlm(rholm, flm, x, LALParams);

    /** Spin corrections */
    REAL8 rho22S;
    REAL8 rho32S;
    REAL8 rho44S;
    REAL8 rho42S;
    REAL8 f21S;
    REAL8 f33S;
    REAL8 f31S;
    REAL8 f43S;
    REAL8 f41S;
      
    const REAL8 a0 = a1 + a2;
    const REAL8 a12 = a1 - a2;
    const REAL8 X12 = X1 - X2;
    const REAL8 a0X12 = a0 * X12;
    const REAL8 a12X12 = a12 * X12;
  
    const REAL8 v = sqrt(x);
    const REAL8 v2 = x;
    const REAL8 v3 = v2 * v;
    const REAL8 v4 = v3 * v;
    const REAL8 v5 = v4 * v;
     
    /* l=m=2 multipole */
    /* spin-orbit */
    const REAL8 cSO_lo = (-0.5*a0 - a12X12/6.);
    const REAL8 cSO_nlo = (-52./63.-19./504.*nu)*a0 - (50./63.+209./504.*nu)*a12X12;
  
    /* SPIN-SPIN contribution */
    REAL8 cSS_lo;

    if (useTidal)
    {
    	// Include tidal code here if needed
    	printf("Wrong");
    	cSS_lo = 0.0;
    } 
    else
    {
        cSS_lo = 0.5 * a0 * a0; 
    }

    /* rho_22^S: Eq. (80) of Damour & Nagar, PRD 90, 044018 (2014) */
    rho22S = cSO_lo*v3 + cSS_lo*v4 + cSO_nlo*v5 ;
    
    /* l>=3, m=even: multipoles rewritten in compact and self-explanatory form */
    rho32S = (a0-a12X12)/(3.*(1.-3.*nu))*v;
    rho44S = (-19./30.*a0 -  (1.-21.*nu)/(30.-90.*nu)*a12X12)*v3;
    rho42S = ( -1./30.*a0 - (19.-39.*nu)/(30.-90.*nu)*a12X12)*v3;
  
    /* l>=2, m=odd: multipoles rewritten in compact and self-explanatory form */
    f21S = -1.5*a12*v + ((110./21.+79./84.*nu)*a12-13./84.*a0X12)*v3;
    f33S = ((-0.25+2.5*nu)*a12-1.75*a0X12) * v3;
    f31S = ((-2.25+6.5*nu)*a12+0.25*a0X12) * v3;
    f43S = ((5.-10.*nu)*a12 - 5.*a0X12) / (-4.+8.*nu)*v;
    f41S = f43S;
    
    /* Amplitudes (correct with spin terms) */
    flm[0] = gsl_pow_int(rholm[0], 2);
    flm[0] = (X12*flm[0] + f21S);
  
    flm[1] = gsl_pow_int(rholm[1]+ rho22S, 2);
  
    flm[2] = gsl_pow_int(rholm[2], 3);
    flm[2] = (X12*flm[2] + f31S);
  
    flm[3] = gsl_pow_int(rholm[3]+ rho32S, 3);
  
    flm[4] = gsl_pow_int(rholm[4], 3);
    flm[4] = (X12*flm[4] + f33S);
  
    flm[5] = gsl_pow_int(rholm[5], 4);
    flm[5] = (X12*flm[5] + f41S);
  
    flm[6] = gsl_pow_int(rholm[6] + rho42S, 4);
  
    flm[7] = gsl_pow_int(rholm[7], 4);
    flm[7] = (X12*flm[7] + f43S);
  
    flm[8] = gsl_pow_int(rholm[8] + rho44S, 4);

	return XLAL_SUCCESS;
}

/** Resummed amplitudes in the general nu-dependent case.
 *  Refs:
 *  . Damour, Iyer & Nagar, PRD 79, 064004 (2009)     [theory]
 *  . Fujita & Iyer, PRD 82, 044051 (2010)            [test-mass 5.5PN]
 *  . Damour, Nagar & Bernuzzi, PRD 87, 084035 (2013) [complete information]
 */
int
XLALSimInspiralEOBPostAdiabaticWavFlm(
	UNUSED REAL8 *rholm,
	UNUSED REAL8 *flm,
	UNUSED REAL8 x,
	UNUSED LALDict *LALParams)
{
	const REAL8 nu = XLALDictLookupREAL8Value(LALParams, "nu");

    /* Coefficients */
    static REAL8 clm[KMAX][6];
  
    const REAL8 nu2 = nu * nu;
    const REAL8 nu3 = nu * nu2;
    const REAL8 nu4 = nu * nu3;
    const REAL8 Pi2 = LAL_PI * LAL_PI;
  
    static REAL8 firstcall = 1;

    if (firstcall)
    {
        firstcall = 0;

        UINT4 k;
    
        for (k=0; k<KMAX; k++)
        {
        	clm[k][0] = 1.;
        }

        for (k=0; k<KMAX; k++)
        {
        	for (int n=1; n<6; n++)
        	{
        		clm[k][n] = 0.;
        	}
        }

        /** (2,1) */
        clm[0][1] = (-1.0535714285714286 + 0.27380952380952384 *nu);
        clm[0][2] = (-0.8327841553287982 - 0.7789824263038548  *nu + 0.13116496598639457*nu2);

        /** (2,2) */
        clm[1][1] = (-1.0238095238095237 + 0.6547619047619048*nu);
        clm[1][2] = (-1.94208238851096   - 1.5601379440665155*nu + 0.4625614134542706*nu2);

        /** (3,1) */
        clm[2][1] = (-0.7222222222222222 - 0.2222222222222222*nu);
        clm[2][2] = (0.014169472502805836 - 0.9455667789001122*nu - 0.46520763187429853*nu2);

        /** (3,2) */
        clm[3][1] = (0.003703703703703704*(328. - 1115.*nu + 320.*nu2))/(-1. + 3.*nu);
        clm[3][2] = (6.235191420376606e-7*(-1.444528e6 + 8.050045e6*nu - 4.725605e6*nu2 - 2.033896e7*nu3 + 3.08564e6*nu4))/((-1. + 3.*nu)*(-1. + 3.*nu));

        /** (3,3) */
        clm[4][1] =  (-1.1666666666666667 + 0.6666666666666666*nu);
        clm[4][2] = (-1.6967171717171716 - 1.8797979797979798*nu + 0.45151515151515154*nu2);

        /** (4,1) */
        clm[5][1] = (0.001893939393939394*(602. - 1385.*nu + 288.*nu2))/(-1. + 2.*nu);
        clm[5][2] = (- 0.36778992787515513);

        /** (4,2) */
        clm[6][1] = (0.0007575757575757576*(1146. - 3530.*nu + 285.*nu2))/(-1. + 3.*nu);
        clm[6][2] = - (3.1534122443213353e-9*(1.14859044e8 - 2.95834536e8*nu - 1.204388696e9*nu2 + 3.04798116e9*nu3 + 3.79526805e8*nu4))/((-1. + 3.*nu)*(-1. + 3.*nu));
        
        /** (4,3) */
        clm[7][1] = (0.005681818181818182*(222. - 547.*nu + 160.*nu2))/(-1. + 2.*nu);
        clm[7][2] = (- 0.9783218202252293);

        /** (4,4) */
        clm[8][1] = (0.0007575757575757576*(1614. - 5870.*nu + 2625.*nu2))/(-1. + 3.*nu);
        clm[8][2] = (3.1534122443213353e-9*(-5.11573572e8 + 2.338945704e9*nu - 3.13857376e8*nu2 - 6.733146e9*nu3 + 1.252563795e9*nu4))/((-1. + 3.*nu)*(-1. + 3.*nu));
       	
        /** (5,1) */
        clm[9][1] = (0.002564102564102564*(319. - 626.*nu + 8.*nu2))/(-1. + 2.*nu);
        clm[9][2] = (- 0.1047896120973044);

        /** (5,2) */
        clm[10][1] = (0.00007326007326007326*(-15828. + 84679.*nu - 104930.*nu2 + 21980.*nu3))/(1. - 5.*nu + 5.*nu2);
        clm[10][2] = (- 0.4629337197600934)*PMTERMS_eps; 

        /** (5,3) */
        clm[11][1] = (0.002564102564102564*(375. - 850.*nu + 176.*nu2))/(-1. + 2.*nu);
        clm[11][2] = (- 0.5788010707241477);

        /** (5,4) */
        clm[12][1] = (0.00007326007326007326*(-17448. + 96019.*nu - 127610.*nu2 + 33320.*nu3))/(1. - 5.*nu + 5.*nu2);
        clm[12][2] = (- 1.0442142414362194)*PMTERMS_eps;

        /** (5,5) */
        clm[13][1] = (0.002564102564102564*(487. - 1298.*nu + 512.*nu2))/(-1. + 2.*nu);
        clm[13][2] = (- 1.5749727622804546);
        
        /** (6,1) */
        clm[14][1] = (0.006944444444444444*(-161. + 694.*nu - 670.*nu2 + 124.*nu3))/(1. - 4.*nu + 3.*nu2);
        clm[14][2] = (- 0.29175486850885135)*PMTERMS_eps;

        /** (6,2) */
        clm[15][1] = (0.011904761904761904*(-74. + 378.*nu - 413.*nu2 + 49.*nu3))/(1. - 5.*nu + 5.*nu2);
        clm[15][2] = ( - 0.24797525070634313)*PMTERMS_eps;

        /** (6,3) */
        clm[16][1] = (0.006944444444444444*(-169. + 742.*nu - 750.*nu2 + 156.*nu3))/(1. - 4.*nu + 3.*nu2);
        clm[16][2] = (- 0.5605554442947213)*PMTERMS_eps;

        /** (6,4) */
        clm[17][1] = (0.011904761904761904*(-86. + 462.*nu - 581.*nu2 + 133.*nu3))/(1. - 5.*nu + 5.*nu2);
        clm[17][2] = (- 0.7228451986855349)*PMTERMS_eps;

        /** (6,5) */
        clm[18][1] = (0.006944444444444444*(-185. + 838.*nu - 910.*nu2 + 220.*nu3))/(1. - 4.*nu + 3.*nu2);
        clm[18][2] = (- 1.0973940686333457)*PMTERMS_eps;

        /** (6,6) */
        clm[19][1] = (0.011904761904761904*(-106. + 602.*nu - 861.*nu2 + 273.*nu3))/(1. - 5.*nu + 5.*nu2); 
        clm[19][2] = (- 1.5543111183867486)*PMTERMS_eps;

        /** (7,1) */
        clm[20][1] = (0.0014005602240896359*(-618. + 2518.*nu - 2083.*nu2 + 228.*nu3))/(1. - 4.*nu + 3.*nu2);
        clm[20][2] = ( - 0.1508235111143767)*PMTERMS_eps;

        /** (7,2) */
        clm[21][1] = (0.00006669334400426837*(16832. - 123489.*nu + 273924.*nu2 - 190239.*nu3 + 32760.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
        clm[21][2] = (- 0.351319484450667)*PMTERMS_eps;
        
        /** (7,3) */
        clm[22][1] = (0.0014005602240896359*(-666. + 2806.*nu - 2563.*nu2 + 420.*nu3))/(1. - 4.*nu + 3.*nu2);
        clm[22][2] = (- 0.37187416047628863)*PMTERMS_eps;

        /** (7,4) */
        clm[23][1] = (0.00006669334400426837*(17756. - 131805.*nu + 298872.*nu2 - 217959.*nu3 + 41076.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
        clm[23][2] = (- 0.6473746896670599)*PMTERMS_eps;
        
        /** (7,5) */
        clm[24][1] = (0.0014005602240896359*(-762. + 3382.*nu - 3523.*nu2 + 804.*nu3))/(1. - 4.*nu + 3.*nu2);
        clm[24][2] = (- 0.8269193364414116)*PMTERMS_eps;

        /** (7,6) */
        clm[25][1] = (0.0006002400960384153*(2144. - 16185.*nu + 37828.*nu2 - 29351.*nu3 + 6104.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
        clm[25][2] = (- 1.1403265020692532)*PMTERMS_eps;
        
        /** (7,7) */
        clm[26][1] = (0.0014005602240896359*(-906. + 4246.*nu - 4963.*nu2 + 1380.*nu3))/(1. - 4.*nu + 3.*nu2);
        clm[26][2] = (- 1.5418467934923434)*PMTERMS_eps;

        /** (8,1) */
        clm[27][1] = (0.00005482456140350877*(20022. - 126451.*nu + 236922.*nu2 - 138430.*nu3 + 21640.*nu4))/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
        clm[27][2] = (- 0.26842133517043704)*PMTERMS_eps;

        /** (8,2) */
        clm[28][1] = (0.0003654970760233918*(2462. - 17598.*nu + 37119.*nu2 - 22845.*nu3 + 3063.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
        clm[28][2] = (- 0.2261796441029474)*PMTERMS_eps;

        /** (8,3) */
        clm[29][1] = (0.00005482456140350877*(20598. - 131059.*nu + 249018.*nu2 - 149950.*nu3 + 24520.*nu4))/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
        clm[29][2] = (- 0.4196774909106648)*PMTERMS_eps;

        /** (8,4) */
        clm[30][1] = (0.0003654970760233918*(2666. - 19434.*nu + 42627.*nu2 - 28965.*nu3 + 4899.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
        clm[30][2] = (- 0.47652059150068155)*PMTERMS_eps;

        /** (8,5) */
        clm[31][1] = (0.00027412280701754384*(4350. - 28055.*nu + 54642.*nu2 - 34598.*nu3 + 6056.*nu4))/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
        clm[31][2] = (- 0.7220789990670207)*PMTERMS_eps;

        /** (8,6) */
        clm[32][1] = (0.0010964912280701754*(1002. - 7498.*nu + 17269.*nu2 - 13055.*nu3 + 2653.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
        clm[32][2] = (- 0.9061610303170207)*PMTERMS_eps;

        /** (8,7) */
        clm[33][1] = (0.00005482456140350877*(23478. - 154099.*nu + 309498.*nu2 - 207550.*nu3 + 38920.*nu4))/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
        clm[33][2] = (- 1.175404252991305)*PMTERMS_eps;

        /** (8,8) */
        clm[34][1] = (0.0003654970760233918*(3482. - 26778.*nu + 64659.*nu2 - 53445.*nu3 + 12243.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
        clm[34][2] = (- 1.5337092502821381)*PMTERMS_eps;
    }

    /** Compute EulerLogs */
    const REAL8 el1 = XLALSimInspiralEOBPostAdiabaticEulerLog(x, 1);
    const REAL8 el2 = XLALSimInspiralEOBPostAdiabaticEulerLog(x, 2);
    const REAL8 el3 = XLALSimInspiralEOBPostAdiabaticEulerLog(x, 3);
    const REAL8 el4 = XLALSimInspiralEOBPostAdiabaticEulerLog(x, 4);
    const REAL8 el5 = XLALSimInspiralEOBPostAdiabaticEulerLog(x, 5);
    const REAL8 el6 = XLALSimInspiralEOBPostAdiabaticEulerLog(x, 6);
    const REAL8 el7 = XLALSimInspiralEOBPostAdiabaticEulerLog(x, 7);

    /** Coefs with Eulerlogs */
    clm[0][3] = (2.9192806270460925  - 1.019047619047619   *el1);
    clm[0][4] = (-1.28235780892213   + 1.073639455782313   *el1);
    clm[0][5] = (-3.8466571723355227 + 0.8486467106683944  *el1)*PMTERMS_eps;
  
    clm[1][3] = (12.736034731834051  - 2.902228713904598 *nu - 1.9301558466099282*nu2 + 0.2715020968103451*nu3 - 4.076190476190476*el2);
    clm[1][4] = (-2.4172313935587004 + 4.173242630385488 *el2);
    clm[1][5] = (-30.14143102836864  + 7.916297736025627 *el2);

    clm[2][3] = (1.9098284139598072 - 0.4126984126984127*el1+ (-4.646868015386534 + (0.21354166666666666)*Pi2)*nu + 2.3020866307903347*nu2 - 0.5813492634480288*nu3);  
    clm[2][4] = (0.5368150316615179 + 0.2980599647266314*el1);
    clm[2][5] = (1.4497991763035063 - 0.0058477188106817735*el1)*PMTERMS_eps;
    
    clm[3][3] = (6.220997955214429 - 1.6507936507936507*el2);
    clm[3][4] = (-3.4527288879001268 + 2.005408583186361*el2)*PMTERMS_eps;
  
    clm[4][3] = (14.10891386831863 - 3.7142857142857144*el3 + (-5.031429681429682 + (0.21354166666666666)*Pi2)*nu - 1.7781727531727531*nu2 + 0.25923767590434255*nu3);
    clm[4][4] = (-6.723375314944128 + 4.333333333333333*el3);
    clm[4][5] = (-29.568699895427518 + 6.302092352092352*el3)*PMTERMS_eps;
  
    clm[5][3] = (0.6981550175535535 - 0.2266955266955267*el1);
    clm[5][4] = (-0.7931524512893319 + 0.2584672482399755*el1)*PMTERMS_eps;
  
    clm[6][3] = 4.550378418934105e-12*(8.48238724511e11 - 1.9927619712e11*el2);
    clm[6][4] = (-0.6621921297263365 + 0.787251738160829*el2)*PMTERMS_eps;
  
    clm[7][3] = (8.519456157072423 - 2.0402597402597404*el3)*PMTERMS_eps;
    clm[7][4] = (-5.353216984886716 + 2.5735094451003544*el3)*PMTERMS_eps;
  
    clm[8][3] = (15.108111214795123 - 3.627128427128427*el4);
    clm[8][4] = (-8.857121657199649 + 4.434988849534304*el4)*PMTERMS_eps;
  
    clm[9][3] = (0.642701885362399 - 0.14414918414918415*el1)*PMTERMS_eps;
    clm[9][4] = (-0.07651588046467575 + 0.11790664036817883*el1)*PMTERMS_eps;
  
    clm[10][3] = (2.354458371550237 - 0.5765967365967366*el2)*PMTERMS_eps;
  
    clm[11][3] = (5.733973288504755 - 1.2973426573426574*el3)*PMTERMS_eps;
    clm[11][4] = (-1.9573287625526001 + 1.2474448628294783*el3)*PMTERMS_eps;
  
    clm[12][3] = (10.252052781721588 - 2.3063869463869464*el4)*PMTERMS_eps;
  
    clm[13][3] = (15.939827047208668 - 3.6037296037296036*el5)*PMTERMS_eps;
    clm[13][4] = (-10.272578060123237 + 4.500041838503377*el5)*PMTERMS_eps;
  
    clm[14][3] = (0.21653486654395454 - 0.10001110001110002*el1)*PMTERMS_eps;
  
    clm[15][3] = (1.7942694138754138 - 0.40004440004440006*el2)*PMTERMS_eps;
  
    clm[16][3] = (4.002558222882566 - 0.9000999000999002*el3)*PMTERMS_eps;
  
    clm[17][3] = (7.359388663371044 - 1.6001776001776002*el4)*PMTERMS_eps;
  
    clm[18][3] = (11.623366217471297 - 2.5002775002775004*el5)*PMTERMS_eps;
  
    clm[19][3] = (16.645950799433503 - 3.6003996003996006*el6)*PMTERMS_eps;
  
    clm[20][3] = (0.2581280702019663 - 0.07355557607658449*el1)*PMTERMS_eps;
  
    clm[22][3] = (3.0835293524055283 - 0.6620001846892604*el3)*PMTERMS_eps;
  
    clm[24][3] = (8.750589067052443 - 1.838889401914612*el5)*PMTERMS_eps;
  
    clm[26][3] = (17.255875091408523 - 3.6042232277526396*el7)*PMTERMS_eps;
    
    /** rho_lm */
    const REAL8 x2  = x*x;
    const REAL8 x3  = x*x2;
    const REAL8 x4  = x*x3;
    const REAL8 x5  = x*x4;
    const REAL8 xn[] = {1.,x,x2,x3,x4,x5};

    for (int k=0; k<KMAX; k++) {
        rholm[k] = clm[k][0];
    
        for (int n=1; n<6; n++) {
            rholm[k] += clm[k][n] * xn[n];
        }
    }

    /* Amplitudes */
    for (int k = 0; k < KMAX; k++) {
        flm[k] = gsl_pow_int(rholm[k], LINDEX[k]);
    }

	return XLAL_SUCCESS;
}

double
XLALSimInspiralEOBPostAdiabaticEulerLog(
	UNUSED REAL8 x,
	UNUSED INT4 m)
{
	REAL8 logm;
	REAL8 EulerLog;

	logm = log((REAL8)m);

	EulerLog = LAL_GAMMA + LAL_LN2 + logm + 0.5*log(x);

	return EulerLog;
}

int
XLALSimInspiralEOBPostAdiabaticHorizonFluxS(
	UNUSED REAL8 *hatFH,
	UNUSED REAL8 x,
	UNUSED REAL8 Heff,
	UNUSED REAL8 jhat,
	UNUSED LALDict *LALParams)
{
	const REAL8 chi1 = XLALDictLookupREAL8Value(LALParams, "chi1");
	const REAL8 chi2 = XLALDictLookupREAL8Value(LALParams, "chi2");
	const REAL8 X1 = XLALDictLookupREAL8Value(LALParams, "X1");
	const REAL8 X2 = XLALDictLookupREAL8Value(LALParams, "X2");

	REAL8 x4;
	REAL8 v5;
	REAL8 FH22_S;
    REAL8 FH22;
    REAL8 FH21;

    REAL8 cv5[2];
    REAL8 cv8[2];

    x4 = x * x * x * x;
    v5 = sqrt(x4 * x);
  
    /* Coefficients of the v^5 term (Alvi leading order) */
    cv5[0] = -1. / 4. * chi1 * (1.+3.*chi1*chi1) * X1 * X1 * X1;
    cv5[1] = -1. / 4. * chi2 * (1.+3.*chi2*chi2) * X2 * X2 * X2;
  
    /** Coefficients of the v^8=x^4 term */
    cv8[0] = 0.5 * (1.+sqrt(1.-chi1*chi1)) * (1.+3.*chi1*chi1) * X1 * X1 * X1 * X1;
    cv8[1] = 0.5 * (1.+sqrt(1.-chi2*chi2)) * (1.+3.*chi2*chi2) * X2 * X2 * X2 * X2;
  
    FH22_S = (cv5[0]+cv5[1]) * v5;
    FH22 = (cv8[0]+cv8[1]) * x4;
    FH21 = 0.0;
  
    /* Newton-normalized horizon flux: use only l=2 fluxes */
    *hatFH = FH22_S + FH22 + FH21;

	return XLAL_SUCCESS;
}

int
XLALSimInspiralEOBPostAdiabaticHorizonFlux(
	UNUSED REAL8 *hatFH,
	UNUSED REAL8 x,
	UNUSED REAL8 Heff,
	UNUSED REAL8 jhat,
	UNUSED LALDict *LALParams)
{
	const REAL8 nu = XLALDictLookupREAL8Value(LALParams, "nu");

    REAL8 rhoHlm[2]; /* only 21,22 multipoles -> k=0,1 */
    REAL8 FlmHLO[2];
    REAL8 FlmH[2];
  
    /* Shorthands */
    REAL8 nu2 = nu * nu;
    REAL8 nu3 = nu2 * nu;

    REAL8 x2 = x * x;
    REAL8 x3 = x2 * x;
    REAL8 x4 = x3 * x;
    REAL8 x5 = x4 * x;
    REAL8 x9 = x5 * x;
    REAL8 x10 = x9 * x;
    
    /* The Newtonian asymptotic contribution */
    const REAL8 FNewt22 = 32./5.*x5;
  
    /* Compute leading-order part (nu-dependent) */
    FlmHLO[1] = 32./5.*(1-4*nu+2*nu2)*x9;
    FlmHLO[0] = 32./5.*(1-4*nu+2*nu2)*x10;
    
    /* Compute rho_lm */
    REAL8 c1[2];
    REAL8 c2[2];
    REAL8 c3[2];
    REAL8 c4[2];
    
    c1[0] = 0.58121;
    c2[0] = 1.01059;
    c3[0] = 7.955729;
    c4[0] = 1.650228;
  
    c1[1] = (4.-21.*nu + 27.*nu2 - 8.*nu3) / (4.*(1.-4.*nu+2.*nu2));
    c2[1] = 4.78752;
    c3[1] = 26.760136;
    c4[1] = 43.861478;
    
    rhoHlm[1] = 1. + c1[1]*x + c2[1]*x2 + c3[1]*x3 + c4[1]*x4;
    rhoHlm[0] = 1. + c1[0]*x + c2[0]*x2 + c3[0]*x3 + c4[0]*x4;
    
    /* Compute horizon multipolar flux (only l=2) */
    const REAL8 Heff2 = Heff * Heff;
    const REAL8 jhat2 = jhat * jhat;
  
    FlmH[1] = FlmHLO[1] * Heff2 * gsl_pow_int(rhoHlm[1], 4);
    FlmH[0] = FlmHLO[0] * jhat2 * gsl_pow_int(rhoHlm[0], 4);
    
    /* Sum over multipoles and normalize to the 22 Newtonian multipole */
    *hatFH = (FlmH[0]+FlmH[1])/FNewt22;

	return XLAL_SUCCESS;
}

double
XLALSimInspiralEOBPostAdiabaticdpphiFunc(
	REAL8 prstar_sol,
	void *params)
{
	struct PostAdiabaticRootSolveParams *prstarParams = (struct PostAdiabaticRootSolveParams *)params;

	const REAL8 nu = XLALDictLookupREAL8Value(prstarParams->LALParams, "nu");
	const REAL8 z3 = XLALDictLookupREAL8Value(prstarParams->LALParams, "z3");
	const REAL8 S = XLALDictLookupREAL8Value(prstarParams->LALParams, "S");
	const REAL8 Sstar = XLALDictLookupREAL8Value(prstarParams->LALParams, "Sstar");

	REAL8 r = prstarParams->r;
	REAL8 rc = prstarParams->rc;
    REAL8 drc_dr = prstarParams->drcBydr;
    REAL8 uc2 = prstarParams->uc2;
    REAL8 duc_dr = prstarParams->ducBydr;
    REAL8 prstar = prstarParams->prstar;
    REAL8 pphi = prstarParams->pphi;
    REAL8 dpphi_dr = prstarParams->dpphiBydr;
    REAL8 A = prstarParams->A;
    REAL8 B = prstarParams->B;
    REAL8 dA = prstarParams->dA;
    LALDict *LALParams = prstarParams->LALParams;

	REAL8 ggm[14];

	UNUSED REAL8 G;
	UNUSED REAL8 dG_dr;
	UNUSED REAL8 dG_dprstar;
	UNUSED REAL8 dG_dprstarbyprstar;

	REAL8 H;
	REAL8 Heff;
	REAL8 Heff_orb;
	REAL8 dHeff_dr;
	REAL8 dHeff_dprstar;
	REAL8 dHeff_dpphi;
	REAL8 d2Heff_dprstar20;

	REAL8 E;
	REAL8 Omg;

	REAL8 Heff_orb_f;
	REAL8 Heff_f;
	REAL8 E_f;

	REAL8 psi;
	REAL8 r_omg2;
	REAL8 v_phi;
	REAL8 x;
	REAL8 jhat;
	REAL8 flux;

	REAL8 prefactor;
	REAL8 main_sol;

	REAL8 total;

	XLALSimInspiralEOBPostAdiabaticsGSDynamics(ggm, r, rc, drc_dr, prstar_sol, LALParams);

	G = ggm[2]*S + ggm[3]*Sstar;
    dG_dr = ggm[6]*S + ggm[7]*Sstar;
    dG_dprstar = ggm[4]*S + ggm[5]*Sstar;
    dG_dprstarbyprstar = ggm[10]*S + ggm[11]*Sstar;

    XLALSimInspiralEOBPostAdiabaticHamiltonianS(&H, &Heff, &Heff_orb, &dHeff_dr, &dHeff_dprstar, &dHeff_dpphi, &d2Heff_dprstar20, r, rc, drc_dr, pphi, prstar, S, Sstar, A, dA, LALParams);

    E = nu * H;
    Omg = dHeff_dpphi / E;

    Heff_orb_f = sqrt(A * (1.+pphi*pphi*uc2));
    Heff_f = G*pphi + Heff_orb_f;
    E_f = sqrt(1 + 2*nu*(Heff_f-1));

    psi = (duc_dr + dG_dr*rc*sqrt(A/(pphi*pphi) + A*uc2)/A)/(-0.5*dA);
    r_omg2 = pow(((1./sqrt(rc*rc*rc*psi))+G)/(E_f), -2./3.);
    v_phi = r_omg2 * Omg;
    x = v_phi * v_phi;
    jhat = pphi / (r_omg2*v_phi);
    XLALSimInspiralEOBPostAdiabaticFluxS(&flux, x, Omg, r_omg2, E, Heff, jhat, r, prstar_sol, 0.0, LALParams);

    prefactor = sqrt(A/B) * 1. / (nu*H*Heff_orb);
    main_sol = prstar_sol*(1+2*z3*A/(rc*rc)*prstar_sol*prstar_sol) + Heff_orb*pphi*dG_dprstar;

    total = prefactor * main_sol;

	return dpphi_dr*total - flux;
}

double
XLALSimInspiralEOBPostAdiabaticdprstarFunc(
	REAL8 pphi_sol,
	void *params)
{
	struct PostAdiabaticRootSolveParams *pphiParams = (struct PostAdiabaticRootSolveParams *)params;

	const REAL8 z3 = XLALDictLookupREAL8Value(pphiParams->LALParams, "z3");

    REAL8 dAuc2Bydr = pphiParams->dAuc2Bydr;
   	REAL8 HeffOrb = pphiParams->HeffOrb;
    REAL8 dGBydr = pphiParams->dGBydr;
    REAL8 dGBydprstar = pphiParams->dGBydprstar;
    REAL8 dprstarBydr = pphiParams->dprstarBydr;
    REAL8 dA = pphiParams->dA;
   	REAL8 prstar = pphiParams->prstar;
   	REAL8 A = pphiParams->A;
    REAL8 uc2 = pphiParams->uc2;

    REAL8 aCoeff;
    REAL8 bCoeff;
    REAL8 cCoeff;

    REAL8 result;

    aCoeff = dAuc2Bydr;
    bCoeff = 2 * HeffOrb * (dGBydr+dGBydprstar*dprstarBydr);
    cCoeff = dA + 2*prstar*dprstarBydr*(1+2*z3*A*uc2*prstar*prstar) + z3*dAuc2Bydr*gsl_pow_int(prstar, 4);

    result = aCoeff*pphi_sol*pphi_sol + bCoeff*pphi_sol + cCoeff;

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
    solver_type = gsl_root_fsolver_falsepos;

	gsl_root_fsolver *solver;
	solver = gsl_root_fsolver_alloc(solver_type);

	REAL8 F_lower = Func(x_lower, params);
	REAL8 F_upper = Func(x_upper, params);

	if (F_lower*F_upper >= 0.0)
	{
		x_lower = -5.e-2;
		x_upper = 5.e-2;
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
XLALFourthOrderFiniteDifferenceDerivative(
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
XLALSecondOrderFiniteDifferenceDerivative(
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
				derivativeVec->data[i] += secondOrderCoeffs[0][j]*YVec->data[j];
			}
		}
		else if (i == vecLength-1)
		{
			for (j = 0; j <= 2; j++)
			{
				derivativeVec->data[i] += secondOrderCoeffs[2][j]*YVec->data[i+j-2];
			}
		}
		else
		{
			for (j = 0; j <= 2; j++)
			{
				derivativeVec->data[i] += secondOrderCoeffs[1][j]*YVec->data[i+j-1];
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

int
XLALSimInspiralEOBPostAdiabatic(
	UNUSED REAL8TimeSeries **dynamics,
	/**<< OUTPUT, real part of the modes */
	// const REAL8 phiC,
	/**<< coalescence orbital phase (rad) */
	UNUSED REAL8 deltaT,
	/**<< sampling time step */
	UNUSED const REAL8 m1SI,
	/**<< mass-1 in SI unit */
	UNUSED const REAL8 m2SI
	/**<< mass-2 in SI unit */)
{
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

	chi1 = -0.99;
	chi2 = -0.99;

	useSpins = 1;
	useTidal = 0;

	fInit = 0.00059105892307;

	rSize = 1000;

	UNUSED const char centrifugalRadius[] = "LO";
	UNUSED const char useFlm[] = "SSLO";

	// dr = 1.9774396742060579e-01;

	PA_order = 8;

	// set by the code
	UNUSED REAL8 m1;
	UNUSED REAL8 m2;
	UNUSED REAL8 nu;
	UNUSED REAL8 X1;
	UNUSED REAL8 X2;

	UNUSED REAL8 a1;
	UNUSED REAL8 a2;
	UNUSED REAL8 aK;
	UNUSED REAL8 aK2;

	UNUSED REAL8 S1;
	UNUSED REAL8 S2;
	UNUSED REAL8 S;
	UNUSED REAL8 Sstar;

	UNUSED REAL8 c3;

	UNUSED REAL8 C_Q1;
	UNUSED REAL8 C_Q2;

	UNUSED REAL8 time_units_factor;

	UNUSED REAL8 f0;
	UNUSED REAL8 r0;
	UNUSED REAL8 rMin;
	UNUSED REAL8 z3;

	// set the above variables

	nu = XLALSimInspiralEOBPostAdiabaticSymmetricMassRatio(q);
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

	time_units_factor = XLALSimInspiralEOBPostAdiabaticTimeUnitsFactor(M);

	f0 = fInit / time_units_factor;
	r0 = XLALSimInspiralEOBPostAdiabaticDynr0Kepler(fInit); // should be (f0)
	rMin = XLALSimInspiralEOBPostAdiabaticFinalRadius(q, a1, a2);

	z3 = XLALSimInspiralEOBPostAdiabaticz3(nu);

	dr = (r0-rMin) / (rSize-1);

	LALDict *LALparams = XLALCreateDict();

	XLALDictInsertREAL8Value(LALparams, "M", M);
	XLALDictInsertREAL8Value(LALparams, "q", q);

	XLALDictInsertINT4Value(LALparams, "useSpins", useSpins);
	XLALDictInsertINT4Value(LALparams, "useTidal", useTidal);

	XLALDictInsertREAL8Value(LALparams, "chi1", chi1);
	XLALDictInsertREAL8Value(LALparams, "chi2", chi2);

	XLALDictInsertUINT4Value(LALparams, "rSize", rSize);

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
    REAL8Vector *fluxesVec = XLALCreateREAL8Vector(rSize);

    memset(rVec->data, 0, rVec->length * sizeof(REAL8));
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
    memset(fluxesVec->data, 0, fluxesVec->length * sizeof(REAL8));

	REAL8 r;

	UINT4 i;

	for (i = 0; i < rSize; i++)
	{
		r = r0 - i*dr;

		rVec->data[i] = r;

		XLALSimInspiralEOBPostAdiabaticMetricS(&AVec->data[i], &BVec->data[i], &dAVec->data[i], &d2AVec->data[i], &dBVec->data[i], r, LALparams);

		XLALSimInspiralEOBPostAdiabaticGetCentrifugalRadiusLO(&rcVec->data[i], &drcBydrVec->data[i], &d2rBydr2Vec->data[i], r, LALparams);

		REAL8 ggm[14];
		XLALSimInspiralEOBPostAdiabaticsGSDynamics(ggm, r, rcVec->data[i], drcBydrVec->data[i], 0.0, LALparams);

	    REAL8 G;
	    G = ggm[2]*S + ggm[3]*Sstar;

        dGBydrVec->data[i] = ggm[6]*S + ggm[7]*Sstar;
        dGBydprstarVec->data[i] = ggm[4]*S + ggm[5]*Sstar;
        dGBydprstarbyprstarVec->data[i] = ggm[10]*S + ggm[11]*Sstar;
        
        /* Circular values for flux */
        G0Vec->data[i] = G;
        dGBydr0Vec->data[i] = dGBydrVec->data[i];

        /* Auxiliary variables */
        REAL8 uc;

        sqrtAByBVec->data[i] = sqrt(AVec->data[i]/BVec->data[i]);
        uc = 1. / rcVec->data[i];
        uc2Vec->data[i] = uc * uc;
        ducBydrVec->data[i] = -uc2Vec->data[i] * drcBydrVec->data[i];
        dAuc2BydrVec->data[i] = uc2Vec->data[i] * (dAVec->data[i]-2*AVec->data[i]*uc*drcBydrVec->data[i]);

        REAL8 aCoeff;
        REAL8 bCoeff;
        REAL8 cCoeff;

        aCoeff = dAuc2BydrVec->data[i]*dAuc2BydrVec->data[i] - 4*AVec->data[i]*uc2Vec->data[i]*dGBydrVec->data[i]*dGBydrVec->data[i];
        bCoeff = 2*dAVec->data[i]*dAuc2BydrVec->data[i] - 4*AVec->data[i]*dGBydrVec->data[i]*dGBydrVec->data[i];
        cCoeff = dAVec->data[i] * dAVec->data[i];

        REAL8 root1;
        REAL8 root2;

        gsl_poly_solve_quadratic(aCoeff, bCoeff, cCoeff, &root1, &root2);

        REAL8 j02;

        // dGdr sign determines choice of solution
        if (((dGBydr0Vec->data[i] > 0.) && (aCoeff > 0.)) || ((dGBydr0Vec->data[i] < 0.) && (aCoeff < 0.)))
        {
            j02 = root2;
        }
        else
        {
            j02 = root1;
        }

        // circular orbital angular momentum
        pphiVec->data[i] = sqrt(j02);
        pphi0Vec->data[i] = sqrt(j02);

        REAL8 H;
        REAL8 dHeff_dr;
        REAL8 dHeff_dprstar;
        REAL8 d2Heff_dprstar20;

        XLALSimInspiralEOBPostAdiabaticHamiltonianS(&H, &HeffVec->data[i], &HeffOrbVec->data[i], &dHeff_dr, &dHeff_dprstar, &dHeffBydpphiVec->data[i], &d2Heff_dprstar20, r, rcVec->data[i], drcBydrVec->data[i], pphiVec->data[i], prstarVec->data[i], S, Sstar, AVec->data[i], dAVec->data[i], LALparams);

        EVec->data[i] = nu * H;

        // Circular orbital frequency 
        OmgVec->data[i] = dHeffBydpphiVec->data[i] / EVec->data[i];

        /* Circular real orbital frequency */
        OmgOrbVec->data[i] = (pphiVec->data[i]*AVec->data[i]*uc2Vec->data[i]) / (EVec->data[i]*HeffOrbVec->data[i]);
    }

    REAL8Vector *rReverseVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *pphiReverseVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *dpphiBydrReverseVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *dpphiBydrVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *dpphiBydr0Vec = XLALCreateREAL8Vector(rSize);

    memset(rReverseVec->data, 0, rReverseVec->length * sizeof(REAL8));
    memset(pphiReverseVec->data, 0, pphiReverseVec->length * sizeof(REAL8));
    memset(dpphiBydrReverseVec->data, 0, dpphiBydrReverseVec->length * sizeof(REAL8));
    memset(dpphiBydrVec->data, 0, dpphiBydrVec->length * sizeof(REAL8));
    memset(dpphiBydr0Vec->data, 0, dpphiBydr0Vec->length * sizeof(REAL8));

    REAL8Vector *dpphiBydrVec1 = XLALCreateREAL8Vector(rSize);
    memset(dpphiBydrVec1->data, 0, dpphiBydrVec1->length * sizeof(REAL8));

    *dpphiBydrVec1 = XLALFourthOrderFiniteDifferenceDerivative(rVec, pphiVec);

    REAL8Vector *dpphiBydrVec2 = XLALCreateREAL8Vector(rSize);
    memset(dpphiBydrVec2->data, 0, dpphiBydrVec2->length * sizeof(REAL8));

    *dpphiBydrVec2 = XLALSecondOrderFiniteDifferenceDerivative(rVec, pphiVec);

    *rReverseVec = XLALReverseREAL8Vector(rVec);
    *pphiReverseVec = XLALReverseREAL8Vector(pphiVec);

    *dpphiBydrReverseVec = XLALPostAdiabaticSplineDerivative(rReverseVec, pphiReverseVec);

    *dpphiBydrVec = XLALReverseREAL8Vector(dpphiBydrReverseVec);
    *dpphiBydr0Vec = XLALReverseREAL8Vector(dpphiBydrReverseVec); // Figure out how to copy REAL8Vecotr here instead of reversing it again

    for (i = 0; i < rSize; i++)
	{
		printf("%.18e\n", fabs(dpphiBydrVec1->data[i] - dpphiBydrVec2->data[i]));
	}

    // test begins here
 //    REAL8Vector *testxVec = XLALCreateREAL8Vector(1000);
 //    REAL8Vector *testx2Vec = XLALCreateREAL8Vector(1000);
 //    REAL8Vector *test2xVec = XLALCreateREAL8Vector(1000);
 //    REAL8Vector *testxSecondVec = XLALCreateREAL8Vector(1000);
 //    REAL8Vector *testxFourthVec = XLALCreateREAL8Vector(1000);

 //    memset(testxVec->data, 0, testxVec->length * sizeof(REAL8));
 //    memset(testx2Vec->data, 0, testx2Vec->length * sizeof(REAL8));
 //    memset(test2xVec->data, 0, test2xVec->length * sizeof(REAL8));
 //    memset(testxSecondVec->data, 0, testxSecondVec->length * sizeof(REAL8));
 //    memset(testxFourthVec->data, 0, testxFourthVec->length * sizeof(REAL8));

 //    for (i = 0; i < 1000; i++)
	// {
	// 	testxVec->data[i] = -20 + i * 0.10000000000000000000;
	// 	testx2Vec->data[i] = gsl_pow_int(testxVec->data[i], 4);
	// 	test2xVec->data[i] = 4 * gsl_pow_int(testxVec->data[i], 3);
	// }

	// *testxSecondVec = XLALSecondOrderFiniteDifferenceDerivative(testxVec, testx2Vec);
	// *testxFourthVec = XLALFourthOrderFiniteDifferenceDerivative(testxVec, testx2Vec);

	// for (i = 0; i < rSize; i++)
	// {
	// 	printf("%.18e %.18e %.18e %.18e %.18e\n", testxVec->data[i], testx2Vec->data[i], test2xVec->data[i], testxSecondVec->data[i], testxFourthVec->data[i]);
	// }

    exit(0);

	REAL8Vector *prstarReverseVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *dprstarBydrReverseVec = XLALCreateREAL8Vector(rSize);

    memset(prstarReverseVec->data, 0, prstarReverseVec->length * sizeof(REAL8));
    memset(dprstarBydrReverseVec->data, 0, dprstarBydrReverseVec->length * sizeof(REAL8));

    REAL8Vector *a1Vec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *a2Vec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *aKVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *SstarVec = XLALCreateREAL8Vector(rSize);

    memset(a1Vec->data, a1, a1Vec->length * sizeof(REAL8));
    memset(a2Vec->data, a2, a2Vec->length * sizeof(REAL8));
    memset(aKVec->data, aK, aKVec->length * sizeof(REAL8));
    memset(SstarVec->data, Sstar, SstarVec->length * sizeof(REAL8));

    UNUSED REAL8 SpinAlignedH;
    SpinEOBHCoeffs Hcoeffs;
    XLALSimIMRCalculateSpinEOBHCoeffs(&Hcoeffs, nu, aK, 4);

    TidalEOBParams tidal1, tidal2;

    tidal1.mByM = X1 / (X1+X2);
	tidal1.lambda2Tidal = 0.0;
	tidal1.omega02Tidal = 0.0;
	tidal1.lambda3Tidal = 0.0;
	tidal1.omega03Tidal = 0.0;
	tidal1.quadparam = 1.0;

	tidal2.mByM = X2 / (X1+X2);
	tidal2.lambda2Tidal = 0.0;
	tidal2.omega02Tidal = 0.0;
	tidal2.lambda3Tidal = 0.0;
	tidal2.omega03Tidal = 0.0;
	tidal2.quadparam = 1.0;

	Hcoeffs.tidal1 = &tidal1;
 	Hcoeffs.tidal2 = &tidal2;

    SpinAlignedH = XLALSimIMRSpinEOBHamiltonian(nu, rVec, prstarVec, a1Vec, a2Vec, aKVec, SstarVec, 0, &Hcoeffs);

 //   	printf("%.18f\n", SpinAlignedH);
	// exit(0);

    UINT4 n;
    UINT4 parity;

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
    		r = rVec->data[i];

    		if (parity)
    		{
    			// Odd PA orders: corrections to prstar only
    			REAL8 Heff_orb_f;
    			REAL8 Heff_f;
    			REAL8 E_f;
    			REAL8 psi;
    			REAL8 r_omg2;
    			REAL8 v_phi;
    			REAL8 x;
    			REAL8 jhat;
    			REAL8 prstar_fake;
    			REAL8 ddotr_fake;
    			REAL8 Fphi;

    			// Variables for which Kepler's law is still valid
    			Heff_orb_f = sqrt(AVec->data[i] * (1.+pphiVec->data[i]*pphiVec->data[i]*uc2Vec->data[i]));
                Heff_f = G0Vec->data[i]*pphiVec->data[i] + Heff_orb_f;
                E_f = sqrt(1 + 2*nu*(Heff_f-1));

                psi = (ducBydrVec->data[i]+dGBydr0Vec->data[i]*rcVec->data[i]*sqrt(AVec->data[i]/(pphiVec->data[i]*pphiVec->data[i]) + AVec->data[i]*uc2Vec->data[i])/AVec->data[i]) / (-0.5*dAVec->data[i]);
                r_omg2 = pow((((1./sqrt(rcVec->data[i]*rcVec->data[i]*rcVec->data[i]*psi))+G0Vec->data[i])/(E_f)), -2./3.);
                v_phi = r_omg2 * OmgVec->data[i];
                x = v_phi * v_phi;
                jhat = pphiVec->data[i] / (r_omg2*v_phi);

                // Add save to REAL8Array for these variables

                ddotr_fake  = 0.0;
                prstar_fake = 0.0;

                XLALSimInspiralEOBPostAdiabaticFluxS(
                	&Fphi,
                	x,
                	OmgVec->data[i],
                	r_omg2,
                	EVec->data[i],
                	HeffVec->data[i],
                	jhat,
                	r,
                	prstar_fake,
                	ddotr_fake,
                	LALparams);

                fluxesVec->data[i] = Fphi;

                struct PostAdiabaticRoot prstarRoot;

                struct PostAdiabaticRootSolveParams prstarParams;
                prstarParams.r = r;
			    prstarParams.rc = rcVec->data[i];
			    prstarParams.drcBydr = drcBydrVec->data[i];
			    prstarParams.uc2 = uc2Vec->data[i];
			    prstarParams.ducBydr = ducBydrVec->data[i];
			    prstarParams.prstar = prstarVec->data[i];
			    prstarParams.pphi = pphiVec->data[i];
			    prstarParams.dpphiBydr = dpphiBydrVec->data[i];
			    prstarParams.A = AVec->data[i];
			    prstarParams.B = BVec->data[i];
			    prstarParams.dA = dAVec->data[i];
			    prstarParams.LALParams = LALparams;

				REAL8 x_lower = 0.9 * prstarVec->data[i];
				REAL8 x_upper = 1.1 * prstarVec->data[i];

			    XLALSimInspiralEOBPostAdiabaticRootFinder(
			    	&prstarRoot,
					XLALSimInspiralEOBPostAdiabaticdpphiFunc,
					&prstarParams,
					x_lower,
					x_upper,
					1.e-17,
					1.e-20
				);

				prstarVec->data[i] = prstarRoot.root;

                // REAL8 dHeff_dprstarbyprstar, dr_dtbyprstar;

				// dHeff_dprstarbyprstar = pphiVec->data[i]*dGBydprstarbyprstarVec->data[i] + (1+2*z3*AVec->data[i]*uc2Vec->data[i]*prstarVec->data[i]*prstarVec->data[i])/HeffOrbVec->data[i];
	   //          dr_dtbyprstar = sqrtAByBVec->data[i]/(EVec->data[i]) * dHeff_dprstarbyprstar;
	   //          prstarVec->data[i] = Fphi / dpphiBydrVec->data[i] / dr_dtbyprstar;

	            REAL8 ggm[14];
				XLALSimInspiralEOBPostAdiabaticsGSDynamics(ggm, r, rcVec->data[i], drcBydrVec->data[i], prstarVec->data[i], LALparams);

                dGBydrVec->data[i] = ggm[6]*S + ggm[7]*Sstar;
                dGBydprstarVec->data[i] = ggm[4]*S + ggm[5]*Sstar;
                dGBydprstarbyprstarVec->data[i] = ggm[10]*S + ggm[11]*Sstar;
    		}
    		else
    		{
    			// Even PA orders: corrections to pphi only
    			struct PostAdiabaticRoot pphiRoot;

    			struct PostAdiabaticRootSolveParams pphiParams;
                pphiParams.dAuc2Bydr = dAuc2BydrVec->data[i];
			    pphiParams.HeffOrb = HeffOrbVec->data[i];
			    pphiParams.dGBydr = dGBydrVec->data[i];
			    pphiParams.dGBydprstar = dGBydprstarVec->data[i];
			    pphiParams.dprstarBydr = dprstarBydrVec->data[i];
			    pphiParams.dA = dAVec->data[i];
			   	pphiParams.prstar = prstarVec->data[i];
			   	pphiParams.A = AVec->data[i];
			    pphiParams.uc2 = uc2Vec->data[i];

			    REAL8 x_lower = 0.9 * pphiVec->data[i];
				REAL8 x_upper = 1.1 * pphiVec->data[i];

			    XLALSimInspiralEOBPostAdiabaticRootFinder(
			    	&pphiRoot,
					XLALSimInspiralEOBPostAdiabaticdprstarFunc,
					&pphiParams,
					x_lower,
					x_upper,
					1.e-17,
					1.e-20
				);

				pphiVec->data[i] = pphiRoot.root;
    		}

    		REAL8 H;
			REAL8 dHeff_dr;
			REAL8 dHeff_dprstar;
			REAL8 d2Heff_dprstar20;

			XLALSimInspiralEOBPostAdiabaticHamiltonianS(&H, &HeffVec->data[i], &HeffOrbVec->data[i], &dHeff_dr, &dHeff_dprstar, &dHeffBydpphiVec->data[i], &d2Heff_dprstar20, r, rcVec->data[i], drcBydrVec->data[i], pphiVec->data[i], prstarVec->data[i], S, Sstar, AVec->data[i], dAVec->data[i], LALparams);

	        EVec->data[i] = nu * H;

	        /* Circular orbital frequency */
	        OmgVec->data[i] = dHeffBydpphiVec->data[i] / EVec->data[i];
	        
	        /* Circular real orbital frequency */
	        OmgOrbVec->data[i] = (pphiVec->data[i]*AVec->data[i]*uc2Vec->data[i]) / (EVec->data[i]*HeffOrbVec->data[i]);

	        dtBydrVec->data[i] = EVec->data[i] / (sqrtAByBVec->data[i]*dHeff_dprstar); // dt_dr = 1/dr_dt 
	        dphiBydrVec->data[i] = OmgVec->data[i] * dtBydrVec->data[i];
    	}

		// if (parity)
		// {
		// 	for (i = 0; i < rSize; i++)
		//     {
		//         prstarReverseVec->data[i] = prstarVec->data[rSize-i-1];
		//     }

	 //        acc = gsl_interp_accel_alloc ();
		//     spline = gsl_spline_alloc (gsl_interp_akima, rVec->length);

		//     gsl_spline_init(spline, rReverseVec->data, prstarReverseVec->data, rReverseVec->length);

		//     for (i = 0; i < rSize; i++)
		//     {
		//         dprstarBydrReverseVec->data[i] = gsl_spline_eval_deriv(spline, rReverseVec->data[i], acc);
		//     }

		//     gsl_spline_free(spline);
		//     gsl_interp_accel_free(acc);

	 //        for (i = 0; i < rSize; i++)
		//     {
		//         dprstarBydrVec->data[i] = dprstarBydrReverseVec->data[rSize-i-1];
		//     }
		// }
		// else
		// {
		// 	for (i = 0; i < rSize; i++)
		//     {
		//         pphiReverseVec->data[i] = pphiVec->data[rSize-i-1];
		//     }

	 //        acc = gsl_interp_accel_alloc ();
		//     spline = gsl_spline_alloc (gsl_interp_akima, rVec->length);

		//     gsl_spline_init(spline, rReverseVec->data, pphiReverseVec->data, rReverseVec->length);

		//     for (i = 0; i < rSize; i++)
		//     {
		//         dpphiBydrReverseVec->data[i] = gsl_spline_eval_deriv(spline, rReverseVec->data[i], acc);
		//     }

		//     gsl_spline_free(spline);
		//     gsl_interp_accel_free(acc);

	 //        for (i = 0; i < rSize; i++)
		//     {
		//         dpphiBydrVec->data[i] = dpphiBydrReverseVec->data[rSize-i-1];
		//     }
		// }

		if (n == 2)
			exit(0);
    }

    REAL8Vector *dtBydrReverseVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *tReverseVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *tVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *dphiBydrReverseVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *phiVec = XLALCreateREAL8Vector(rSize);
    REAL8Vector *phiReverseVec = XLALCreateREAL8Vector(rSize);

    memset(dtBydrReverseVec->data, 0, dtBydrReverseVec->length * sizeof(REAL8));
    memset(tReverseVec->data, 0, tReverseVec->length * sizeof(REAL8));
    memset(tVec->data, 0, tReverseVec->length * sizeof(REAL8));
    memset(dphiBydrReverseVec->data, 0, dphiBydrReverseVec->length * sizeof(REAL8));
    memset(phiVec->data, 0, phiVec->length * sizeof(REAL8));
    memset(phiReverseVec->data, 0, phiReverseVec->length * sizeof(REAL8));

    for (i = 0; i < rSize; i++)
    {
        dtBydrReverseVec->data[i] = dtBydrVec->data[rSize-i-1];
    }

    gsl_interp_accel *acc1 = gsl_interp_accel_alloc ();
    gsl_spline *spline1 = gsl_spline_alloc (gsl_interp_cspline, rVec->length);

    gsl_spline_init(spline1, rReverseVec->data, dtBydrReverseVec->data, rReverseVec->length);

    for (i = 1; i < rSize; i++)
    {
        tReverseVec->data[i] = gsl_spline_eval_integ(spline1, rReverseVec->data[i-1], rReverseVec->data[i], acc1);
    }

    gsl_spline_free(spline1);
    gsl_interp_accel_free(acc1);

    for (i = 0; i < rSize; i++)
    {
        tVec->data[i] = tReverseVec->data[rSize-i-1];
    }

    REAL8 tOffset = fabs(tVec->data[0]);

    for (i = 0; i < rSize; i++)
    {
        tVec->data[i] = tVec->data[i] + tOffset;
    }

    for (i = 0; i < rSize; i++)
    {
        dphiBydrReverseVec->data[i] = dphiBydrVec->data[rSize-i-1];
    }

    gsl_interp_accel *acc2 = gsl_interp_accel_alloc ();
    gsl_spline *spline2 = gsl_spline_alloc (gsl_interp_cspline, rVec->length);

    gsl_spline_init(spline2, rReverseVec->data, dphiBydrReverseVec->data, rReverseVec->length);

    for (i = 1; i < rSize; i++)
    {
        phiReverseVec->data[i] = gsl_spline_eval_integ(spline2, rReverseVec->data[i-1], rReverseVec->data[i], acc2);
    }

    gsl_spline_free(spline2);
    gsl_interp_accel_free(acc2);

    for (i = 0; i < rSize; i++)
    {
        phiVec->data[i] = phiReverseVec->data[rSize-i-1];
    }

    REAL8 phiOffset = fabs(phiVec->data[0]);

    for (i = 0; i < rSize; i++)
    {
        phiVec->data[i] = phiVec->data[i] + phiOffset;
    }

    // Output ready
    // tVec, rVec, phiVec, pphiVec, prstarVec, pphi0Vec, rcVec, AVec, dpphiBydrVec

    return XLAL_SUCCESS;
}