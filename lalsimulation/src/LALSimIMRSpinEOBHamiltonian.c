/*
*  Copyright (C) 2011 Craig Robinson, Enrico Barausse, Yi Pan
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
 * \author Craig Robinson, Yi Pan, Andrea Taracchini
 *
 * Functions for calculating the effective one-body Hamiltonian for
 * spinning binaries, as described in
 * Taracchini et al. ( PRD 86, 024011 (2012), arXiv 1202.0790 ).
 * All equation numbers in this file refer to equations of this paper,
 * unless otherwise specified.
 * This code borrows hugely from a C implementation originally written
 * by Enrico Barausse, following Barausse and Buonanno
 * PRD 81, 084024 (2010) and PRD 84, 104027 (2011), henceforth BB1 and BB2
 */

#ifndef _LALSIMIMRSPINEOBHAMILTONIAN_C
#define _LALSIMIMRSPINEOBHAMILTONIAN_C

#include <stdio.h>
#include <math.h>

#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

#include "LALSimIMRSpinEOB.h"

#include "LALSimIMRSpinEOBHamiltonian.h"

//#include "fresnel.h"

#include <math.h>

#define TwoToOneThird 1.25992104989487316476721060728
#define ThreeToOneThird 1.44224957030740838232163831078

static const double sqrt_pi_2   = 1.2533141373155002512078826424; /* sqrt(pi/2) */
static const double sqrt_2_pi   = 0.7978845608028653558798921199; /* sqrt(2/pi) */
static const double _1_sqrt_2pi = 0.3989422804014326779399460599; /* 1/sqrt(2*pi) */
static const double pi_2        = 1.5707963267948966192313216916; /* pi/2 */

static double f_data_a[18] =
{
		  0.76435138664186000189,
    -0.43135547547660179313,
		  0.43288199979726653054,
    -0.26973310338387111029,
    0.08416045320876935378,
    -0.01546524484461381958,
    0.00187855423439822018,
    -0.00016264977618887547,
    0.00001057397656383260,
    -0.00000053609339889243,
    0.00000002181658454933,
    -0.00000000072901621186,
    0.00000000002037332546,
    -0.00000000000048344033,
    0.00000000000000986533,
    -0.00000000000000017502,
    0.00000000000000000272,
    -0.00000000000000000004
};

static double f_data_b[17] =
{
		  0.63041404314570539241,
    -0.42344511405705333544,
		  0.37617172643343656625,
    -0.16249489154509567415,
		  0.03822255778633008694,
    -0.00564563477132190899,
		  0.00057454951976897367,
    -0.00004287071532102004,
		  0.00000245120749923299,
    -0.00000011098841840868,
		  0.00000000408249731696,
    -0.00000000012449830219,
		  0.00000000000320048425,
    -0.00000000000007032416,
		  0.00000000000000133638,
    -0.00000000000000002219,
		  0.00000000000000000032
};

static double fresnel_cos_0_8(double x)
{
    double x_8 = x/8.0;
    double xx = 2.0*x_8*x_8 - 1.0;

    double t0 = 1.0;
    double t1 = xx;
    double sumC = f_data_a[0] + f_data_a[1]*t1;
    double t2;
    int n;
    for (n=2; n < 18; n++)
    {
        t2 = 2.0*xx*t1 - t0;
        sumC += f_data_a[n]*t2;
        t0 = t1; t1 = t2;
    }
    return _1_sqrt_2pi*sqrt(x)*sumC;
}

static double fresnel_sin_0_8(double x)
{
    double x_8 = x/8.0;
    double xx = 2.0*x_8*x_8 - 1.0;
    double t0 = 1.;
    double t1 = xx;
    double ot1 = x_8;
    double ot2 = 2.0*x_8*t1 - ot1;
    double sumS = f_data_b[0]*ot1 + f_data_b[1]*ot2;
    int n;
    double t2;
    for (n=2; n < 17; n++)
    {
        t2 = 2.0*xx*t1 - t0;
        ot1 = ot2;
        ot2 = 2.0*x_8*t2 - ot1;
        sumS += f_data_b[n]*ot2;
        t0 = t1; t1 = t2;
    }
    return _1_sqrt_2pi*sqrt(x)*sumS;
}

static double f_data_e[41] =
{
    0.97462779093296822410,
    -0.02424701873969321371,
    0.00103400906842977317,
    -0.00008052450246908016,
    0.00000905962481966582,
    -0.00000131016996757743,
    0.00000022770820391497,
    -0.00000004558623552026,
    0.00000001021567537083,
    -0.00000000251114508133,
    0.00000000066704761275,
    -0.00000000018931512852,
    0.00000000005689898935,
    -0.00000000001798219359,
    0.00000000000594162963,
    -0.00000000000204285065,
    0.00000000000072797580,
    -0.00000000000026797428,
    0.00000000000010160694,
    -0.00000000000003958559,
    0.00000000000001581262,
    -0.00000000000000646411,
    0.00000000000000269981,
    -0.00000000000000115038,
    0.00000000000000049942,
    -0.00000000000000022064,
    0.00000000000000009910,
    -0.00000000000000004520,
    0.00000000000000002092,
    -0.00000000000000000982,
    0.00000000000000000467,
    -0.00000000000000000225,
    0.00000000000000000110,
    -0.00000000000000000054,
    0.00000000000000000027,
    -0.00000000000000000014,
    0.00000000000000000007,
    -0.00000000000000000004,
    0.00000000000000000002,
    -0.00000000000000000001,
    0.00000000000000000001
};

static double f_data_f[35] =
{
    0.99461545179407928910,
    -0.00524276766084297210,
    0.00013325864229883909,
    -0.00000770856452642713,
    0.00000070848077032045,
    -0.00000008812517411602,
    0.00000001359784717148,
    -0.00000000246858295747,
    0.00000000050925789921,
    -0.00000000011653400634,
    0.00000000002906578309,
    -0.00000000000779847361,
    0.00000000000222802542,
    -0.00000000000067239338,
    0.00000000000021296411,
    -0.00000000000007041482,
    0.00000000000002419805,
    -0.00000000000000861080,
    0.00000000000000316287,
    -0.00000000000000119596,
    0.00000000000000046444,
    -0.00000000000000018485,
    0.00000000000000007527,
    -0.00000000000000003131,
    0.00000000000000001328,
    -0.00000000000000000574,
    0.00000000000000000252,
    -0.00000000000000000113,
    0.00000000000000000051,
    -0.00000000000000000024,
    0.00000000000000000011,
    -0.00000000000000000005,
    0.00000000000000000002,
    -0.00000000000000000001,
    0.00000000000000000001
};

static double fresnel_cos_8_inf(double x)
{
    double xx = 128.0/(x*x) - 1.0;   /* 2.0*(8/x)^2 - 1 */
    double t0 = 1.0;
    double t1 = xx;
    double sumP = f_data_e[0] + f_data_e[1]*t1;
    double sumQ = f_data_f[0] + f_data_f[1]*t1;
    double t2;
    int n;
    for(n = 2; n < 35; n++)
    {
        t2 = 2.0*xx*t1 - t0;
        sumP += f_data_e[n]*t2; /*  sumP += f_data_e[n]*ChebyshevT(n,xx) */
        sumQ += f_data_f[n]*t2; /*  sumQ += f_data_f[n]*ChebyshevT(n,xx) */
        t0 = t1; t1 = t2;
    }
    for(n = 35; n < 41; n++)
    {
        t2 = 2.0*xx*t1 - t0;
        sumP += f_data_e[n]*t2; /*  sumP += f_data_e[n]*ChebyshevT(n,xx) */
        t0 = t1; t1 = t2;
    }
    return 0.5 - _1_sqrt_2pi*(0.5*sumP*cos(x)/x - sumQ*sin(x))/sqrt(x);
}

static double fresnel_sin_8_inf(double x)
{
    double xx = 128.0/(x*x) - 1.0;   /* 2.0*(8/x)^2 - 1 */
    double t0 = 1.0;
    double t1 = xx;
    double sumP = f_data_e[0] + f_data_e[1]*t1;
    double sumQ = f_data_f[0] + f_data_f[1]*t1;
    double t2;
    int n;
    for(n = 2; n < 35; n++)
    {
        t2 = 2.0*xx*t1 - t0;
        sumP += f_data_e[n]*t2; /*  sumP += f_data_e[n]*ChebyshevT(n,xx) */
        sumQ += f_data_f[n]*t2; /*  sumQ += f_data_f[n]*ChebyshevT(n,xx) */
        t0 = t1; t1 = t2;
    }
    for(n = 35; n < 41; n++)
    {
        t2 = 2.0*xx*t1 - t0;
        sumP += f_data_e[n]*t2; /*  sumQ += f_data_f[n]*ChebyshevT(n,xx) */
        t0 = t1; t1 = t2;
    }
    return 0.5 - _1_sqrt_2pi*(0.5*sumP*sin(x)/x + sumQ*cos(x))/sqrt(x);
}


static double fresnel_c(double x)
{
    double xx = x*x*pi_2;
    double ret_val;
    if(xx<=8.0)
        ret_val = fresnel_cos_0_8(xx);
    else
        ret_val = fresnel_cos_8_inf(xx);
    return (x<0.0) ? -ret_val : ret_val;
}

static double fresnel_s(double x)
{
    double xx = x*x*pi_2;
    double ret_val;
    if(xx<=8.0)
        ret_val = fresnel_sin_0_8(xx);
    else
        ret_val = fresnel_sin_8_inf(xx);
    return (x<0.0) ? -ret_val : ret_val;
}

UNUSED static double fresnel_c1(double x)
{
    return fresnel_c(x*sqrt_2_pi);
}

UNUSED static double fresnel_s1(double x)
{
    return fresnel_s(x*sqrt_2_pi);
}
/*------------------------------------------------------------------------------------------
 *
 *          Prototypes of functions defined in this code.
 *
 *------------------------------------------------------------------------------------------
 */

/**
 * This function calculates the DeltaR potential function in the spin EOB Hamiltonian
 */
static REAL8 XLALSimIMRSpinEOBHamiltonianDeltaR (SpinEOBHCoeffs * coeffs,
				/**<< Pre-computed coefficients which appear in the function */
						 const REAL8 r,
				/**<< Current orbital radius (in units of total mass) */
						 const REAL8 eta,
				/**<< Symmetric mass ratio */
						 const REAL8 a
				/**<< Normalized deformed Kerr spin */
  );

static REAL8 XLALSimIMRSpinEOBHamiltonian (const REAL8 eta,
					   REAL8Vector * restrict x,
					   REAL8Vector * restrict p,
					   REAL8Vector * restrict s1Vec,
					   REAL8Vector * restrict s2Vec,
					   REAL8Vector * restrict sigmaKerr,
					   REAL8Vector * restrict sigmaStar,
					   int tortoise,
					   SpinEOBHCoeffs * coeffs);

static int XLALSimIMRCalculateSpinEOBHCoeffs (SpinEOBHCoeffs * coeffs,
					      const REAL8 eta,
					      const REAL8 a,
					      const UINT4
					      SpinAlignedEOBversion);

static REAL8 XLALSimIMRSpinEOBHamiltonianDeltaT (SpinEOBHCoeffs * coeffs,
						 const REAL8 r,
						 const REAL8 eta,
						 const REAL8 a);

static REAL8 XLALSimIMRSpinAlignedEOBCalcOmega (const REAL8 values[],
						SpinEOBParams * funcParams);

static REAL8 XLALSimIMRSpinAlignedEOBNonKeplerCoeff (const REAL8 values[],
						     SpinEOBParams *
						     funcParams);

static double GSLSpinAlignedHamiltonianWrapper (double x, void *params);


/*------------------------------------------------------------------------------------------
 *
 *          Defintions of functions.
 *
 *------------------------------------------------------------------------------------------
 */

/**
 * Function to compute the enhancement of k2tidal due to the presence of f-mode resonance
 */
static REAL8 XLALSimIMRTEOBk2eff (
                                  REAL8 u, /**<< Inverse of radial separation in units of M */
                                  REAL8 eta, /**<< Symmetric mass ratio */
                                  TidalEOBParams * tidal /**<< Tidal parameters */
)
{
    REAL8 w02 = tidal->omega02Tidal;
    REAL8 eps = 64./5.*TwoToOneThird*pow(w02, 5./3.)*eta;
    REAL8 bigomega = pow(1./u, 1.5)*w02/2.;
    REAL8 factorQ = 4. - TwoToOneThird*pow(1./u, 2.5)*pow(w02, 5./3.);
    REAL8 calR = 1./(bigomega*bigomega - 1.) + 10./3./factorQ;
    REAL8 yval = sqrt(3./LAL_PI)*factorQ/5./sqrt(eps);
    REAL8 k2Tidaleff = 0.25 + 3./4.*bigomega*bigomega*(calR + sqrt(LAL_PI/3.)/sqrt(eps)*((1. + 2.*fresnel_s(yval))*cos(0.5*LAL_PI*yval*yval) - (1. + 2.*fresnel_c(yval))*sin(0.5*LAL_PI*yval*yval)));
    return k2Tidaleff;
}

/**
 * Function to compute the u-derivative of the enhancement of k2tidal due to the presence of f-mode resonance
 */
static REAL8 XLALSimIMRTEOBk2eff_u (
                                  REAL8 u, /**<< Inverse of radial separation in units of M */
                                  REAL8 eta, /**<< Symmetric mass ratio */
                                  TidalEOBParams * tidal /**<< Tidal parameters */
)
{
    REAL8 w02 = tidal->omega02Tidal;
    REAL8 eps = 64./5.*TwoToOneThird*pow(w02, 5./3.)*eta;
    REAL8 bigomega = pow(1./u, 1.5)*w02/2.;
    REAL8 bigomega_u = -3./2.*pow(1./u, 2.5)*w02/2.;
    REAL8 factorQ = 4. - TwoToOneThird*pow(1./u, 2.5)*pow(w02, 5./3.);
    REAL8 factorQ_u = 2.5*TwoToOneThird*pow(1./u, 3.5)*pow(w02, 5./3.);
    REAL8 calR = 1./(bigomega*bigomega - 1.) + 10./3./factorQ;
    REAL8 calR_u = -1./(bigomega*bigomega - 1.)/(bigomega*bigomega - 1.)*2.*bigomega*bigomega_u - 10./3./factorQ/factorQ*factorQ_u;
    REAL8 yval = sqrt(3./LAL_PI)*factorQ/5./sqrt(eps);
    REAL8 yval_u = sqrt(3./LAL_PI)*factorQ_u/5./sqrt(eps);
    REAL8 fresnelS = fresnel_s(yval);
    REAL8 fresnelC = fresnel_c(yval);
    REAL8 costerm = cos(0.5*LAL_PI*yval*yval);
    REAL8 sinterm = sin(0.5*LAL_PI*yval*yval);
    REAL8 k2Tidaleff_u = 3./4.*2.*bigomega*bigomega_u*(calR + sqrt(LAL_PI/3.)/sqrt(eps)*((1. + 2.*fresnelS)*costerm - (1. + 2.*fresnelC)*sinterm))
    + 3./4.*bigomega*bigomega*(calR_u - sqrt(LAL_PI/3.)/sqrt(eps)*LAL_PI*yval*yval_u*(costerm*(1. + 2.*fresnelC) + sinterm*(1. + 2.*fresnelS)));
    return k2Tidaleff_u;
}

/**
 * Function to compute the enhancement of k3tidal due to the presence of f-mode resonance
 */
static REAL8 XLALSimIMRTEOBk3eff (
                                  REAL8 u, /**<< Inverse of radial separation in units of M */
                                  REAL8 eta, /**<< Symmetric mass ratio */
                                  TidalEOBParams * tidal /**<< Tidal parameters */
)
{
    REAL8 w03 = tidal->omega03Tidal;
    REAL8 factorO = 9. - ThreeToOneThird*pow(1./u, 2.5)*pow(w03, 5./3.);
    REAL8 XX = factorO/(4.*ThreeToOneThird*ThreeToOneThird*sqrt(10.)*pow(w03,5./6.)*sqrt(eta));
    REAL8 prefactorO = 5.*sqrt(5.*LAL_PI)/u/u/u*pow(w03,7./6.)/(192.*ThreeToOneThird*ThreeToOneThird*sqrt(eta));
    REAL8 k3Tidaleff = 3./8. + w03*w03/u/u/u*(25./48./factorO + 5./72./(-1. + w03*w03/u/u/u/9.)) + prefactorO*(cos(XX*XX)*(0.5 + fresnel_s(sqrt_2_pi*XX)) - sin(XX*XX)*(0.5 + fresnel_c(sqrt_2_pi*XX)));
//    printf("%.16e %.16e\n",u,k3Tidaleff);
    return k3Tidaleff;
}

/**
 * Function to compute the u-derivative of the enhancement of k3tidal due to the presence of f-mode resonance
 */
static REAL8 XLALSimIMRTEOBk3eff_u (
                                  REAL8 u, /**<< Inverse of radial separation in units of M */
                                  REAL8 eta, /**<< Symmetric mass ratio */
                                  TidalEOBParams * tidal /**<< Tidal parameters */
)
{
    REAL8 u2 = u*u;
    REAL8 u3 = u*u2;
    REAL8 u4 = u*u3;
    REAL8 w03 = tidal->omega03Tidal;
    REAL8 factorO = 9. - ThreeToOneThird*pow(1./u, 2.5)*pow(w03, 5./3.);
    REAL8 factorO_u = -2.5/u*(factorO - 9.);
    REAL8 XX = factorO/(4.*ThreeToOneThird*ThreeToOneThird*sqrt(10.)*pow(w03,5./6.)*sqrt(eta));
    REAL8 XX_u = factorO_u*XX/factorO;
    REAL8 prefactorO = 5.*sqrt(5.*LAL_PI)/u/u/u*pow(w03,7./6.)/(192.*ThreeToOneThird*ThreeToOneThird*sqrt(eta));
    REAL8 prefactorO_u = -3./u*prefactorO;
    REAL8 fresnelArg = sqrt_2_pi*XX;
    REAL8 fresnelS = fresnel_s(fresnelArg);
    REAL8 fresnelC = fresnel_c(fresnelArg);
    REAL8 w03Square = w03*w03;
    REAL8 sinXXSquare = sin(XX*XX);
    REAL8 cosXXSquare = cos(XX*XX);
    REAL8 k3Tidaleff_u = 1./48.*( (810.*u2*w03Square)/(w03Square - 9.*u3)/(w03*w03 - 9.*u3) + 24.* (cosXXSquare*(1. + 2.*fresnelS) - sinXXSquare*(1. + 2.*fresnelC)) * prefactorO_u - 48. * prefactorO * (cosXXSquare*(1. + 2.*fresnelC) + sinXXSquare*(1. + 2.*fresnelS)) * XX * XX_u -  25.*w03Square*(3.*factorO + u*factorO_u)/(u4*factorO*factorO));
    return k3Tidaleff_u;
}

/**
 * Function to compute the quadrupolar tidal contribution to the EOB Delta_u potential for just one NS.
 * This implements the model of Phys.Rev.Lett. 116 (2016) no.18, 181101
 */
static REAL8 XLALSimIMRTEOBdeltaUTidalQuadSingleNS (
                                                    REAL8 u, /**<< Inverse of radial separation in units of M */
                                                    REAL8 u2, /**<< Inverse of radial separation^2 in units of M */
                                                    REAL8 u6, /**<< Inverse of radial separation^6 in units of M */
                                                    REAL8 XNS, /**<< NS mass by M */
                                                    REAL8 XCompanion, /**<< Companion mass by M */
                                                    REAL8 lambda2TidalNS, /**<< NS dimensionless quadrupolar tidal deformability */
                                                    REAL8 k2TidalNSeff /**<< Dynamical enhancement of k2TidalNS */
)
{
    return - 3.*XCompanion/XNS*lambda2TidalNS*k2TidalNSeff*u6* ( 1. + 5./2.*u*XNS + u2*(3. + XNS/8. + 337./28.*XNS*XNS) );
}

/**
 * Function to compute the u-derivative of the quadrupolar tidal contribution to the EOB Delta_u potential for just one NS.
 * This implements the model of Phys.Rev.Lett. 116 (2016) no.18, 181101
 */
static REAL8 XLALSimIMRTEOBdeltaUTidalQuadSingleNS_u (
                                                    REAL8 u, /**<< Inverse of radial separation in units of M */
                                                    REAL8 u2, /**<< Inverse of radial separation^2 in units of M */
                                                    REAL8 u6, /**<< Inverse of radial separation^6 in units of M */
                                                    REAL8 XNS, /**<< NS mass by M */
                                                    REAL8 XCompanion, /**<< Companion mass by M */
                                                    REAL8 lambda2TidalNS, /**<< NS dimensionless quadrupolar tidal deformability */
                                                    REAL8 k2TidalNSeff, /**<< Dynamical enhancement of k2TidalNS */
                                                    REAL8 k2TidalNSeff_u /**<< u-derivative of the dynamical enhancement of k2TidalNS */
)
{
    REAL8 deltaUQSingle = XLALSimIMRTEOBdeltaUTidalQuadSingleNS(u, u2, u6, XNS, XCompanion, lambda2TidalNS, k2TidalNSeff);
    return (6./u +  ( 5./2.*XNS + 2.*u*(3. + XNS/8. + 337./28.*XNS*XNS) )/( 1. + 5./2.*u*XNS + u2*(3. + XNS/8. + 337./28.*XNS*XNS) ) + k2TidalNSeff_u/k2TidalNSeff)*deltaUQSingle;
}


/**
 * Function to compute the quadrupolar tidal contribution to the EOB Delta_u potential.
 * This implements the model of Phys.Rev.Lett. 116 (2016) no.18, 181101
 */
static REAL8 XLALSimIMRTEOBdeltaUTidalQuad (
                                        REAL8 u, /**<< Inverse of radial separation in units of M */
                                        REAL8 u2, /**<< Inverse of radial separation^2 in units of M */
                                        REAL8 u6, /**<< Inverse of radial separation^6 in units of M */
                                        REAL8 eta, /**<< Symmetric mass ratio */
                                        TidalEOBParams *tidal1, /**<< Tidal parameters of body 1 */
                                        TidalEOBParams *tidal2 /**<< Tidal parameters of body 2 */
)
{
    REAL8 k2Tidal1eff = 0.;
    REAL8 k2Tidal2eff = 0.;
    REAL8 deltaUQ = 0.;
    if ( tidal1->lambda2Tidal != 0.) {
        k2Tidal1eff = XLALSimIMRTEOBk2eff(u, eta, tidal1);
        deltaUQ +=  XLALSimIMRTEOBdeltaUTidalQuadSingleNS(u, u2, u6, tidal1->mByM, tidal2->mByM, tidal1->lambda2Tidal, k2Tidal1eff);
    }
    if ( tidal2->lambda2Tidal != 0.) {
        k2Tidal2eff = XLALSimIMRTEOBk2eff(u, eta, tidal2);
        deltaUQ += XLALSimIMRTEOBdeltaUTidalQuadSingleNS(u, u2, u6, tidal2->mByM, tidal1->mByM, tidal2->lambda2Tidal, k2Tidal2eff);
    }
    return deltaUQ;
}

/**
 * Function to compute the u-derivative of the quadrupolar tidal contribution to the EOB Delta_u potential.
 * This implements the model of Phys.Rev.Lett. 116 (2016) no.18, 181101
 */
static REAL8 XLALSimIMRTEOBdeltaUTidalQuad_u (
                                            REAL8 u, /**<< Inverse of radial separation in units of M */
                                            REAL8 u2, /**<< Inverse of radial separation^2 in units of M */
                                            REAL8 u6, /**<< Inverse of radial separation^6 in units of M */
                                            REAL8 eta, /**<< Symmetric mass ratio */
                                            TidalEOBParams *tidal1, /**<< Tidal parameters of body 1 */
                                            TidalEOBParams *tidal2 /**<< Tidal parameters of body 2 */
)
{
    REAL8 k2Tidal1eff = 0.;
    REAL8 k2Tidal2eff = 0.;
    REAL8 k2Tidal1eff_u = 0.;
    REAL8 k2Tidal2eff_u = 0.;
    REAL8 deltaUQ_u = 0.;
    if ( tidal1->lambda2Tidal != 0.) {
        k2Tidal1eff = XLALSimIMRTEOBk2eff(u, eta, tidal1);
        k2Tidal1eff_u = XLALSimIMRTEOBk2eff_u(u, eta, tidal1);
        deltaUQ_u += XLALSimIMRTEOBdeltaUTidalQuadSingleNS_u(u, u2, u6, tidal1->mByM, tidal2->mByM, tidal1->lambda2Tidal, k2Tidal1eff, k2Tidal1eff_u);
    }
    if ( tidal2->lambda2Tidal != 0.) {
        k2Tidal2eff = XLALSimIMRTEOBk2eff(u, eta, tidal2);
        k2Tidal2eff_u = XLALSimIMRTEOBk2eff_u(u, eta, tidal2);
        deltaUQ_u += XLALSimIMRTEOBdeltaUTidalQuadSingleNS_u(u, u2, u6, tidal2->mByM, tidal1->mByM, tidal2->lambda2Tidal, k2Tidal2eff, k2Tidal2eff_u);
    }
    return deltaUQ_u;
}


/**
 * Function to compute the octupolar tidal contribution to the EOB Delta_u potential for just one NS.
 * This implements the model of Phys.Rev.Lett. 116 (2016) no.18, 181101
 */
static REAL8 XLALSimIMRTEOBdeltaUTidalOctuSingleNS (
                                            REAL8 u, /**<< Inverse of radial separation in units of M */
                                            REAL8 u2, /**<< Inverse of radial separation^2 in units of M */
                                            REAL8 u8, /**<< Inverse of radial separation^8 in units of M */
                                            REAL8 XNS, /**<< NS mass by M */
                                            REAL8 XCompanion, /**<< Companion mass by M */
                                            REAL8 lambda3TidalNS, /**<< NS dimensionless octupolar tidal deformability */
                                            REAL8 k3TidalNSeff /**<< Dynamical enhancement of k3TidalNS */
)
{
    return - 15.*XCompanion/XNS*lambda3TidalNS*k3TidalNSeff*u8* (1. + u*(15./2.*XNS - 2.) + u2*(8./3. - 311./24.*XNS + 110./3.*XNS*XNS) );
}

/**
 * Function to compute the u-derivative of the octupolar tidal contribution to the EOB Delta_u potential for just one NS.
 * This implements the model of Phys.Rev.Lett. 116 (2016) no.18, 181101
 */
static REAL8 XLALSimIMRTEOBdeltaUTidalOctuSingleNS_u (
                                                    REAL8 u, /**<< Inverse of radial separation in units of M */
                                                    REAL8 u2, /**<< Inverse of radial separation^2 in units of M */
                                                    REAL8 u8, /**<< Inverse of radial separation^8 in units of M */
                                                    REAL8 XNS, /**<< NS mass by M */
                                                    REAL8 XCompanion, /**<< Companion mass by M */
                                                    REAL8 lambda3TidalNS, /**<< NS dimensionless octupolar tidal deformability */
                                                    REAL8 k3TidalNSeff, /**<< Dynamical enhancement of k3TidalNS */
                                                    REAL8 k3TidalNSeff_u /**<< u-derivative of dynamical enhancement of k3TidalNS */
)
{
    REAL8 deltaUOSingle = XLALSimIMRTEOBdeltaUTidalOctuSingleNS(u, u2, u8, XNS, XCompanion, lambda3TidalNS, k3TidalNSeff);
    return (8./u + ((15./2.*XNS - 2.) + 2.*u*(8./3. - 311./24.*XNS + 110./3.*XNS*XNS) )/(1. + u*(15./2.*XNS - 2.) + u2*(8./3. - 311./24.*XNS + 110./3.*XNS*XNS) ) + k3TidalNSeff_u/k3TidalNSeff)*deltaUOSingle;
}

/**
 * Function to compute the octupolar tidal contribution to the EOB Delta_u potential.
 * This implements the model of Phys.Rev.Lett. 116 (2016) no.18, 181101
 */
static REAL8 XLALSimIMRTEOBdeltaUTidalOctu (
                                        REAL8 u, /**<< Inverse of radial separation in units of M */
                                        REAL8 u2, /**<< Inverse of radial separation^2 in units of M */
                                        REAL8 u8, /**<< Inverse of radial separation^8 in units of M */
                                        REAL8 eta, /**<< Symmetric mass ratio */
                                        TidalEOBParams *tidal1, /**<< Tidal parameters of body 1 */
                                        TidalEOBParams *tidal2 /**<< Tidal parameters of body 2 */
)
{
    REAL8 k3Tidal1eff = 0.;
    REAL8 k3Tidal2eff = 0.;
    REAL8 deltaUO = 0.;
    if ( tidal1->lambda3Tidal != 0.) {
        k3Tidal1eff = XLALSimIMRTEOBk3eff(u, eta, tidal1);
        deltaUO += XLALSimIMRTEOBdeltaUTidalOctuSingleNS(u, u2, u8, tidal1->mByM, tidal2->mByM, tidal1->lambda3Tidal, k3Tidal1eff);
    }
    if ( tidal2->lambda3Tidal != 0.) {
        k3Tidal2eff = XLALSimIMRTEOBk3eff(u, eta, tidal2);
       deltaUO += XLALSimIMRTEOBdeltaUTidalOctuSingleNS(u, u2, u8, tidal2->mByM, tidal1->mByM, tidal2->lambda3Tidal, k3Tidal2eff);
    }
    return deltaUO;
}

/**
 * Function to compute the u-derivative of the octupolar tidal contribution to the EOB Delta_u potential.
 * This implements the model of Phys.Rev.Lett. 116 (2016) no.18, 181101
 */
static REAL8 XLALSimIMRTEOBdeltaUTidalOctu_u (
                                            REAL8 u, /**<< Inverse of radial separation in units of M */
                                            REAL8 u2, /**<< Inverse of radial separation^2 in units of M */
                                            REAL8 u8, /**<< Inverse of radial separation^8 in units of M */
                                            REAL8 eta, /**<< Symmetric mass ratio */
                                            TidalEOBParams *tidal1, /**<< Tidal parameters of body 1 */
                                            TidalEOBParams *tidal2 /**<< Tidal parameters of body 2 */
)
{
    REAL8 k3Tidal1eff = 0.;
    REAL8 k3Tidal2eff = 0.;
    REAL8 k3Tidal1eff_u = 0.;
    REAL8 k3Tidal2eff_u = 0.;
    REAL8 deltaUO_u = 0.;
    if ( tidal1->lambda3Tidal != 0.) {
        k3Tidal1eff = XLALSimIMRTEOBk3eff(u, eta, tidal1);
        k3Tidal1eff_u = XLALSimIMRTEOBk3eff_u(u, eta, tidal1);
        deltaUO_u += XLALSimIMRTEOBdeltaUTidalOctuSingleNS_u(u, u2, u8, tidal1->mByM, tidal2->mByM, tidal1->lambda3Tidal, k3Tidal1eff, k3Tidal1eff_u);
    }
    if ( tidal2->lambda3Tidal != 0.) {
        k3Tidal2eff = XLALSimIMRTEOBk3eff(u, eta, tidal2);
        k3Tidal2eff_u = XLALSimIMRTEOBk3eff_u(u, eta, tidal2);
        deltaUO_u += XLALSimIMRTEOBdeltaUTidalOctuSingleNS_u(u, u2, u8, tidal2->mByM, tidal1->mByM, tidal2->lambda3Tidal, k3Tidal2eff, k3Tidal2eff_u);
    }
    return deltaUO_u;
}

/**
 * Function to compute the tidal contribution to the EOB Delta_u potential.
 * This implements the model of Phys.Rev.Lett. 116 (2016) no.18, 181101
 */
static REAL8 XLALSimIMRTEOBdeltaUTidal (
                                 REAL8 u, /**<< Inverse of radial separation in units of M */
                                 REAL8 eta, /**<< Symmetric mass ratio */
                                 TidalEOBParams *tidal1, /**<< Tidal parameters of body 1 */
                                 TidalEOBParams *tidal2 /**<< Tidal parameters of body 2 */
)
{
    REAL8 u2 = u*u;
    REAL8 u6 = u2*u2*u2;
    REAL8 deltaUQ = XLALSimIMRTEOBdeltaUTidalQuad(u, u2, u6, eta, tidal1, tidal2);
    REAL8 deltaUO = 0.;
    if ( (tidal1->lambda3Tidal != 0. && tidal1->omega03Tidal != 0.) || (tidal2->lambda3Tidal != 0. && tidal2->omega03Tidal != 0.) ) {
        REAL8 u8 = u2*u6;
        deltaUO = XLALSimIMRTEOBdeltaUTidalOctu(u, u2, u8, eta, tidal1, tidal2);
    }
    return deltaUQ + deltaUO;
}

/**
 * Function to compute the u-derivative of the  tidal contribution to the EOB Delta_u potential.
 * This implements the model of Phys.Rev.Lett. 116 (2016) no.18, 181101
 */
static REAL8 XLALSimIMRTEOBdeltaUTidal_u (
                                          REAL8 u, /**<< Inverse of radial separation in units of M */
                                          REAL8 eta, /**<< Symmetric mass ratio */
                                          TidalEOBParams *tidal1, /**<< Tidal parameters of body 1 */
                                          TidalEOBParams *tidal2 /**<< Tidal parameters of body 2 */
                                       )
{
    REAL8 u2 = u*u;
    REAL8 u6 = u2*u2*u2;
    REAL8 deltaUQ_u = XLALSimIMRTEOBdeltaUTidalQuad_u(u, u2, u6, eta, tidal1, tidal2);
    REAL8 deltaUO_u = 0.;
    if ( (tidal1->lambda3Tidal != 0. && tidal1->omega03Tidal != 0.) || (tidal2->lambda3Tidal != 0. && tidal2->omega03Tidal != 0.) ) {
        REAL8 u8 = u2*u6;
        deltaUO_u = XLALSimIMRTEOBdeltaUTidalOctu_u(u, u2, u8, eta, tidal1, tidal2);
    }
    return deltaUQ_u + deltaUO_u;
}

/**
 *
 * Function to calculate the value of the spinning Hamiltonian for given values
 * of the dynamical variables (in a Cartesian co-ordinate system). The inputs are
 * as follows:
 *
 * x - the separation vector r expressed in Cartesian co-ordinates
 * p - the momentum vector (with the radial component tortoise pr*)
 * sigmaKerr - spin of the effective Kerr background (a combination of the individual spin vectors)
 * sigmaStar - spin of the effective particle (a different combination of the individual spins).
 * coeffs - coefficients which crop up in the Hamiltonian. These can be calculated using the
 * XLALCalculateSpinEOBParams() function.
 *
 * The function returns a REAL8, which will be the value of the Hamiltonian if all goes well;
 * otherwise, it will return the XLAL REAL8 failure NaN.
 */
static REAL8
XLALSimIMRSpinEOBHamiltonian (const REAL8 eta,	    /**<< Symmetric mass ratio */
			      REAL8Vector * restrict x,
						    /**<< Position vector */
			      REAL8Vector * restrict p,
						    /**<< Momentum vector (tortoise radial component pr*) */
			      REAL8Vector * restrict s1Vec,
						    /**<< Spin vector 1 */
			      REAL8Vector * restrict s2Vec,
						    /**<< Spin vector 2 */
			      REAL8Vector * restrict sigmaKerr,
						    /**<< Spin vector sigma_kerr */
			      REAL8Vector * restrict sigmaStar,
						    /**<< Spin vector sigma_star */
			      INT4 tortoise,	    /**<< flag to state whether the momentum is the tortoise co-ord */
			      SpinEOBHCoeffs * coeffs
						    /**<< Structure containing various coefficients */
  )
{
  REAL8 r, r2, nx, ny, nz;
  REAL8 sKerr_x, sKerr_y, sKerr_z, a, a2;
  REAL8 sStar_x, sStar_y, sStar_z;
  REAL8 e3_x, e3_y, e3_z;
  REAL8 costheta;		/* Cosine of angle between Skerr and r */
  REAL8 xi2, xi_x, xi_y, xi_z;	/* Cross product of unit vectors in direction of Skerr and r */
  REAL8 vx, vy, vz, pxir, pvr, pn, prT, pr, pf, ptheta2;	/*prT is the tortoise pr */
  REAL8 w2, rho2;
  REAL8 u, u2, u3, u4, u5;
  REAL8 bulk, deltaT, deltaR, Lambda;
  REAL8 D, qq, ww, B, w, MU, nu, BR, wr, nur, mur;
  REAL8 wcos, nucos, mucos, ww_r, Lambda_r;
  REAL8 logTerms, deltaU, deltaU_u, Q, deltaT_r, pn2, pp;
  REAL8 deltaSigmaStar_x, deltaSigmaStar_y, deltaSigmaStar_z;
  REAL8 sx, sy, sz, sxi, sv, sn, s3;
  REAL8 H, Hns, Hs, Hss, Hreal, Hwcos, Hwr, HSOL, HSONL;
  REAL8 m1PlusetaKK;

  /* Terms which come into the 3.5PN mapping of the spins */
  //REAL8 aaa, bbb, a13P5, a23P5, a33P5, b13P5, b23P5, b33P5;
  REAL8 sMultiplier1, sMultiplier2;

  /*Temporary p vector which we will make non-tortoise */
  REAL8 tmpP[3];

  REAL8 csi;

  /* Spin gauge parameters. (YP) simplified, since both are zero. */
  // static const double aa=0., bb=0.;

  //printf( "In Hamiltonian:\n" );
  //printf( "x = %.16e\t%.16e\t%.16e\n", x->data[0], x->data[1], x->data[2] );
  //printf( "p = %.16e\t%.16e\t%.16e\n", p->data[0], p->data[1], p->data[2] );

  r2 =
    x->data[0] * x->data[0] + x->data[1] * x->data[1] +
    x->data[2] * x->data[2];
  r = sqrt (r2);
  nx = x->data[0] / r;
  ny = x->data[1] / r;
  nz = x->data[2] / r;

  sKerr_x = sigmaKerr->data[0];
  sKerr_y = sigmaKerr->data[1];
  sKerr_z = sigmaKerr->data[2];

  sStar_x = sigmaStar->data[0];
  sStar_y = sigmaStar->data[1];
  sStar_z = sigmaStar->data[2];

  a2 = sKerr_x * sKerr_x + sKerr_y * sKerr_y + sKerr_z * sKerr_z;
  a = sqrt (a2);

  if (a != 0.)
    {
      e3_x = sKerr_x / a;
      e3_y = sKerr_y / a;
      e3_z = sKerr_z / a;
    }
  else
    {
      e3_x = 0.;
      e3_y = 0.;
      e3_z = 1.;
    }

  costheta = e3_x * nx + e3_y * ny + e3_z * nz;

  xi2 = 1. - costheta * costheta;

  xi_x = -e3_z * ny + e3_y * nz;
  xi_y = e3_z * nx - e3_x * nz;
  xi_z = -e3_y * nx + e3_x * ny;

  vx = -nz * xi_y + ny * xi_z;
  vy = nz * xi_x - nx * xi_z;
  vz = -ny * xi_x + nx * xi_y;

  w2 = r2 + a2;
  rho2 = r2 + a2 * costheta * costheta;

  u = 1. / r;
  u2 = u * u;
  u3 = u2 * u;
  u4 = u2 * u2;
  u5 = u4 * u;

  //printf( "KK = %.16e\n", coeffs->KK );
  m1PlusetaKK = -1. + eta * coeffs->KK;
  /* Eq. 5.75 of BB1 */
  bulk = 1. / (m1PlusetaKK * m1PlusetaKK) + (2. * u) / m1PlusetaKK + a2 * u2;
  /* Eq. 5.73 of BB1 */
  logTerms =
    1. + eta * coeffs->k0 + eta * log (1. + coeffs->k1 * u + coeffs->k2 * u2 +
				       coeffs->k3 * u3 + coeffs->k4 * u4 +
				       coeffs->k5 * u5 +
				       coeffs->k5l * u5 * log (u));
  //printf( "bulk = %.16e, logTerms = %.16e\n", bulk, logTerms );
  /* Eq. 5.73 of BB1 */
  deltaU = bulk * logTerms;
  if ( (coeffs->tidal1->lambda2Tidal != 0. && coeffs->tidal1->omega02Tidal != 0.) || (coeffs->tidal2->lambda2Tidal != 0. && coeffs->tidal2->omega02Tidal != 0.) ) {
      deltaU += XLALSimIMRTEOBdeltaUTidal(u, eta, coeffs->tidal1, coeffs->tidal2);
  }

  /* Eq. 5.71 of BB1 */
  deltaT = r2 * deltaU;
  /* ddeltaU/du */
  deltaU_u = 2. * (1. / m1PlusetaKK + a2 * u) * logTerms +
    bulk * (eta *
	    (coeffs->k1 +
	     u * (2. * coeffs->k2 +
		  u * (3. * coeffs->k3 +
		       u * (4. * coeffs->k4 +
			    5. * (coeffs->k5 +
				  coeffs->k5l * log (u)) * u))))) / (1. +
								     coeffs->
								     k1 * u +
								     coeffs->
								     k2 * u2 +
								     coeffs->
								     k3 * u3 +
								     coeffs->
								     k4 * u4 +
								     (coeffs->
								      k5 +
								      coeffs->
								      k5l *
								      log (u))
								     * u5);
    if ( (coeffs->tidal1->lambda2Tidal != 0. && coeffs->tidal1->omega02Tidal != 0.) || (coeffs->tidal2->lambda2Tidal != 0. && coeffs->tidal2->omega02Tidal != 0.) ) {
        deltaU_u += XLALSimIMRTEOBdeltaUTidal_u(u, eta, coeffs->tidal1, coeffs->tidal2);
//        FILE *out = fopen("outDImpor.dat","a");
//        fprintf(out, "%.16e %.16e %.16e\n", u, XLALSimIMRTEOBdeltaUTidal(u, eta, coeffs->tidal1, coeffs->tidal2), XLALSimIMRTEOBdeltaUTidal_u(u, eta, coeffs->tidal1, coeffs->tidal2));
//        fclose(out);
    }
  /* ddeltaT/dr */
  deltaT_r = 2. * r * deltaU - deltaU_u;
  /* Eq. 5.39 of BB1 */
  Lambda = w2 * w2 - a2 * deltaT * xi2;
  /* Eq. 5.83 of BB1, inverse */
  D = 1. + log (1. + 6. * eta * u2 + 2. * (26. - 3. * eta) * eta * u3);
  /* Eq. 5.38 of BB1 */
  deltaR = deltaT * D;
  /* See Hns below, Eq. 4.34 of Damour et al. PRD 62, 084011 (2000) */
  qq = 2. * eta * (4. - 3. * eta);
  /* See Hns below. In Sec. II D of BB2 b3 and bb3 coeffs are chosen to be zero. */
  ww = 2. * a * r + coeffs->b3 * eta * a2 * a * u + coeffs->bb3 * eta * a * u;

  /* We need to transform the momentum to get the tortoise co-ord */
  if (tortoise)
    {
      csi = sqrt (deltaT * deltaR) / w2;	/* Eq. 28 of Pan et al. PRD 81, 084041 (2010) */
    }
  else
    {
      csi = 1.0;
    }
  //printf( "csi(miami) = %.16e\n", csi );

  prT = p->data[0] * nx + p->data[1] * ny + p->data[2] * nz;
  /* p->data is BL momentum vector; tmpP is tortoise momentum vector */
  tmpP[0] = p->data[0] - nx * prT * (csi - 1.) / csi;
  tmpP[1] = p->data[1] - ny * prT * (csi - 1.) / csi;
  tmpP[2] = p->data[2] - nz * prT * (csi - 1.) / csi;

  pxir = (tmpP[0] * xi_x + tmpP[1] * xi_y + tmpP[2] * xi_z) * r;
  pvr = (tmpP[0] * vx + tmpP[1] * vy + tmpP[2] * vz) * r;
  pn = tmpP[0] * nx + tmpP[1] * ny + tmpP[2] * nz;

  pr = pn;
  pf = pxir;
  ptheta2 = pvr * pvr / xi2;

  //printf( "pr = %.16e, prT = %.16e\n", pr, prT );

  //printf( " a = %.16e, r = %.16e\n", a, r );
  //printf( "D = %.16e, ww = %.16e, rho = %.16e, Lambda = %.16e, xi = %.16e\npr = %.16e, pf = %.16e, deltaR = %.16e, deltaT = %.16e\n",
  //D, ww, sqrt(rho2), Lambda, sqrt(xi2), pr, pf, deltaR, deltaT );
  /* Eqs. 5.36 - 5.46 of BB1 */
  /* Note that the tortoise prT appears only in the quartic term, explained in Eqs. 14 and 15 of Tarrachini et al. */
  Hns =
    sqrt (1. + prT * prT * prT * prT * qq * u2 + ptheta2 / rho2 +
	  pf * pf * rho2 / (Lambda * xi2) +
	  pr * pr * deltaR / rho2) / sqrt (Lambda / (rho2 * deltaT)) +
    pf * ww / Lambda;

  //printf( "term 1 in Hns: %.16e\n",  prT*prT*prT*prT*qq*u2 );
  //printf( "term 2 in Hns: %.16e\n", ptheta2/rho2 );
  //printf( "term 3 in Hns = %.16e\n", pf*pf*rho2/(Lambda*xi2) );
  //printf( "term 4 in Hns = %.16e\n", pr*pr*deltaR/rho2 );
  //printf( "term 5 in Hns = %.16e\n", Lambda/(rho2*deltaT) );
  //printf( "term 6 in Hns = %.16e\n", pf*ww/Lambda );
  /* Eqs. 5.30 - 5.33 of BB1 */
  B = sqrt (deltaT);
  w = ww / Lambda;
  nu = 0.5 * log (deltaT * rho2 / Lambda);
  MU = 0.5 * log (rho2);
  /* dLambda/dr */
  Lambda_r = 4. * r * w2 - a2 * deltaT_r * xi2;

  ww_r =
    2. * a - (a2 * a * coeffs->b3 * eta) * u2 - coeffs->bb3 * eta * a * u2;
  /* Eqs. 5.47a - 5.47d of BB1 */
  BR =
    (-2. * deltaT + sqrt (deltaR) * deltaT_r) / (2. * sqrt (deltaR * deltaT));
  wr = (-Lambda_r * ww + Lambda * ww_r) / (Lambda * Lambda);
  nur =
    (r / rho2 +
     (w2 * (-4. * r * deltaT + w2 * deltaT_r)) / (2. * deltaT * Lambda));
  mur = (r / rho2 - 1. / sqrt (deltaR));
  /* Eqs. 5.47f - 5.47h of BB1 */
  wcos = -2. * a2 * costheta * deltaT * ww / (Lambda * Lambda);
  nucos = a2 * costheta * w2 * (w2 - deltaT) / (rho2 * Lambda);
  mucos = a2 * costheta / rho2;
  /* Eq. 5.52 of BB1, (YP) simplified */
  //Q = 1. + pvr*pvr/(exp(2.*MU)*xi2) + exp(2.*nu)*pxir*pxir/(B*B*xi2) + pn*pn*deltaR/exp(2.*MU);
  Q =
    1. + pvr * pvr / (rho2 * xi2) +
    deltaT * rho2 / Lambda * pxir * pxir / (B * B * xi2) +
    pn * pn * deltaR / rho2;

  pn2 = pr * pr * deltaR / rho2;
  pp = Q - 1.;

  //printf( "pn2 = %.16e, pp = %.16e\n", pn2, pp );
  //printf( "sigmaKerr = %.16e, sigmaStar = %.16e\n", sKerr_z, sStar_z );
  /* Eq. 5.68 of BB1, (YP) simplified for aa=bb=0. */
  /*
     deltaSigmaStar_x=(- 8.*aa*(1. + 3.*pn2*r - pp*r)*sKerr_x - 8.*bb*(1. + 3.*pn2*r - pp*r)*sStar_x +
     eta*(-8.*sKerr_x - 36.*pn2*r*sKerr_x + 3.*pp*r*sKerr_x + 14.*sStar_x - 30.*pn2*r*sStar_x + 4.*pp*r*sStar_x))/(12.*r);

     deltaSigmaStar_y=(-8.*aa*(1. + 3.*pn2*r - pp*r)*sKerr_y - 8.*bb*(1. + 3.*pn2*r - pp*r)*sStar_y +
     eta*(-8.*sKerr_y - 36.*pn2*r*sKerr_y + 3.*pp*r*sKerr_y + 14.*sStar_y - 30.*pn2*r*sStar_y + 4.*pp*r*sStar_y))/(12.*r);

     deltaSigmaStar_z=(-8.*aa*(1. + 3.*pn2*r - pp*r)*sKerr_z - 8.*bb*(1. + 3.*pn2*r - pp*r)*sStar_z +
     eta*(-8.*sKerr_z - 36.*pn2*r*sKerr_z + 3.*pp*r*sKerr_z + 14.*sStar_z - 30.*pn2*r*sStar_z + 4.*pp*r*sStar_z))/(12.*r);
   */
  deltaSigmaStar_x =
    eta * (-8. * sKerr_x - 36. * pn2 * r * sKerr_x + 3. * pp * r * sKerr_x +
	   14. * sStar_x - 30. * pn2 * r * sStar_x +
	   4. * pp * r * sStar_x) / (12. * r);

  deltaSigmaStar_y =
    eta * (-8. * sKerr_y - 36. * pn2 * r * sKerr_y + 3. * pp * r * sKerr_y +
	   14. * sStar_y - 30. * pn2 * r * sStar_y +
	   4. * pp * r * sStar_y) / (12. * r);

  deltaSigmaStar_z =
    eta * (-8. * sKerr_z - 36. * pn2 * r * sKerr_z + 3. * pp * r * sKerr_z +
	   14. * sStar_z - 30. * pn2 * r * sStar_z +
	   4. * pp * r * sStar_z) / (12. * r);


  /* Now compute the additional 3.5PN terms. */
  /* The following gauge parameters correspond to those given by
   * Eqs. (69) and (70) of BB2 (aaa -> a0, bbb -> b0).
   * In SEOBNRv1 model, we chose to set all of them to zero,
   * described between Eqs. (3) and (4).
   */
  /*
     aaa = -3./2.*eta;
     bbb = -5./4.*eta;
     a1 = eta*eta/2.;
     a2 = -(1./8.)*eta*(-7. + 8.*eta);
     a3 = -((9.*eta*eta)/16.);
     b1 = 1./16.*eta*(9. + 5.*eta);
     b2 = -(1./8.)*eta*(-17. + 5.*eta);
     b3 = -3./8.*eta*eta;
   */
  /*aaa = 0.;
     bbb = 0.;
     a13P5 = 0.;
     a23P5 = 0.;
     a33P5 = 0.;
     b13P5 = 0.;
     b23P5 = 0.;
     b33P5 = 0.;
   */
  /* Eq. 52 of BB2, (YP) simplified for zero gauge parameters */
  /*
     sMultiplier1 =-(2.*(24.*b23P5 + eta*(-353. + 27.*eta) + bbb*(56. + 60.*eta)) +
     2.*(24.*b13P5 - 24.*b23P5 + bbb*(14. - 66.*eta) + 103.*eta - 60.*eta*eta)*pp*
     r + 120.*(2.*b33P5 - 3.*eta*(bbb + eta))*pn2*pn2*r*r +
     (-48.*b13P5 + 4.*bbb*(1. + 3.*eta) + eta*(23. + 3.*eta))*pp*pp*
     r*r + 6.*pn2*r*(16.*b13P5 + 32.*b23P5 + 24.*b33P5 - 47.*eta +
     54.*eta*eta + 24.*bbb*(1. + eta) +
     (24.*b13P5 - 24.*b33P5 - 16.*eta + 21.*eta*eta + bbb*(-2. + 30.*eta))*pp*
     r))/(72.*r*r);
   */
  sMultiplier1 =
    -(2. * eta * (-353. + 27. * eta) +
      2. * (103. * eta - 60. * eta * eta) * pp * r +
      120. * (-3. * eta * eta) * pn2 * pn2 * r * r +
      (eta * (23. + 3. * eta)) * pp * pp * r * r +
      6. * pn2 * r * (-47. * eta + 54. * eta * eta +
		      (-16. * eta +
		       21. * eta * eta) * pp * r)) / (72. * r * r);
  /* Eq. 52 of BB2, (YP) simplified for zero gauge parameters */
  /*
     sMultiplier2 = (-16.*(6.*a23P5 + 7.*eta*(8. + 3.*eta) + aaa*(14. + 15.*eta)) +
     4.*(-24.*a13P5 + 24.*a23P5 - 109.*eta + 51.*eta*eta + 2.*aaa*(-7. + 33.*eta))*
     pp*r + 30.*(-16.*a33P5 + 3.*eta*(8.*aaa + 9.*eta))*pn2*pn2*r*r +
     (96.*a13P5 - 45.*eta - 8.*aaa*(1. + 3.*eta))*pp*pp*r*r -
     6.*pn2*r*(32.*a13P5 + 64.*a23P5 + 48.*a33P5 + 16.*eta + 147.*eta*eta +
     48.*aaa*(1. + eta) + (48.*a13P5 - 48.*a33P5 - 6.*eta + 39.*eta*eta +
     aaa*(-4. + 60.*eta))*pp*r))/(144.*r*r);
   */
  sMultiplier2 =
    (-16. * (7. * eta * (8. + 3. * eta)) +
     4. * (-109. * eta + 51. * eta * eta) * pp * r +
     810. * eta * eta * pn2 * pn2 * r * r - 45. * eta * pp * pp * r * r -
     6. * pn2 * r * (16. * eta + 147. * eta * eta +
		     (-6. * eta +
		      39. * eta * eta) * pp * r)) / (144. * r * r);
  /* Eq. 52 of BB2 */
  deltaSigmaStar_x +=
    sMultiplier1 * sigmaStar->data[0] + sMultiplier2 * sigmaKerr->data[0];
  deltaSigmaStar_y +=
    sMultiplier1 * sigmaStar->data[1] + sMultiplier2 * sigmaKerr->data[1];
  deltaSigmaStar_z +=
    sMultiplier1 * sigmaStar->data[2] + sMultiplier2 * sigmaKerr->data[2];

  /* And now the (calibrated) 4.5PN term */
  deltaSigmaStar_x += coeffs->d1 * eta * sigmaStar->data[0] / (r * r * r);
  deltaSigmaStar_y += coeffs->d1 * eta * sigmaStar->data[1] / (r * r * r);
  deltaSigmaStar_z += coeffs->d1 * eta * sigmaStar->data[2] / (r * r * r);
  deltaSigmaStar_x += coeffs->d1v2 * eta * sigmaKerr->data[0] / (r * r * r);
  deltaSigmaStar_y += coeffs->d1v2 * eta * sigmaKerr->data[1] / (r * r * r);
  deltaSigmaStar_z += coeffs->d1v2 * eta * sigmaKerr->data[2] / (r * r * r);


  //printf( "deltaSigmaStar_x = %.16e, deltaSigmaStar_y = %.16e, deltaSigmaStar_z = %.16e\n",
  //   deltaSigmaStar_x, deltaSigmaStar_y, deltaSigmaStar_z );

  sx = sStar_x + deltaSigmaStar_x;
  sy = sStar_y + deltaSigmaStar_y;
  sz = sStar_z + deltaSigmaStar_z;


  sxi = sx * xi_x + sy * xi_y + sz * xi_z;
  sv = sx * vx + sy * vy + sz * vz;
  sn = sx * nx + sy * ny + sz * nz;

  s3 = sx * e3_x + sy * e3_y + sz * e3_z;
  /* Eq. 3.45 of BB1, second term */
  Hwr =
    (exp (-3. * MU - nu) * sqrt (deltaR) *
     (exp (2. * (MU + nu)) * pxir * pxir * sv -
      B * exp (MU + nu) * pvr * pxir * sxi +
      B * B * xi2 * (exp (2. * MU) * (sqrt (Q) + Q) * sv +
		     pn * pvr * sn * sqrt (deltaR) -
		     pn * pn * sv * deltaR))) / (2. * B * (1. +
							   sqrt (Q)) *
						 sqrt (Q) * xi2);
  /* Eq. 3.45 of BB1, third term */
  Hwcos =
    (exp (-3. * MU - nu) *
     (sn *
      (-(exp (2. * (MU + nu)) * pxir * pxir) +
       B * B * (pvr * pvr - exp (2. * MU) * (sqrt (Q) + Q) * xi2)) -
      B * pn * (B * pvr * sv -
		exp (MU +
		     nu) * pxir * sxi) * sqrt (deltaR))) / (2. * B * (1. +
								      sqrt
								      (Q)) *
							    sqrt (Q));
  /* Eq. 3.44 of BB1, leading term */
  HSOL =
    (exp (-MU + 2. * nu) * (-B + exp (MU + nu)) * pxir * s3) / (B * B *
								sqrt (Q) *
								xi2);
  /* Eq. 3.44 of BB1, next-to-leading term */
  HSONL =
    (exp (-2. * MU + nu) *
     (-(B * exp (MU + nu) * nucos * pxir * (1. + 2. * sqrt (Q)) * sn * xi2) +
      (-(BR * exp (MU + nu) * pxir * (1. + sqrt (Q)) * sv) +
       B * (exp (MU + nu) * nur * pxir * (1. + 2. * sqrt (Q)) * sv +
	    B * mur * pvr * sxi + B * sxi * (-(mucos * pn * xi2) +
					     sqrt (Q) * (mur * pvr -
							 nur * pvr + (-mucos +
								      nucos) *
							 pn * xi2)))) *
      sqrt (deltaR))) / (B * B * (sqrt (Q) + Q) * xi2);
  /* Eq. 3.43 and 3.45 of BB1 */
  Hs = w * s3 + Hwr * wr + Hwcos * wcos + HSOL + HSONL;
  /* Eq. 5.70 of BB1, last term */
  Hss = -0.5 * u3 * (sx * sx + sy * sy + sz * sz - 3. * sn * sn);
  /* Eq. 5.70 of BB1 */
  H = Hns + Hs + Hss;

  /* Add the additional calibrated term */
  H +=
    coeffs->dheffSS * eta * (sKerr_x * sStar_x + sKerr_y * sStar_y +
			     sKerr_z * sStar_z) / (r * r * r * r);
  /* One more calibrated term proportional to S1^2+S2^2. Note that we use symmetric expressions of m1,m2 and S1,S2 */
  /*H += coeffs->dheffSSv2 * eta / (r*r*r*r) / (1.-4.*eta)
   * ( (sKerr_x*sKerr_x + sKerr_y*sKerr_y + sKerr_z*sKerr_z)*(1.-4.*eta+2.*eta*eta)
   +(sKerr_x*sStar_x + sKerr_y*sStar_y + sKerr_z*sStar_z)*(-2.*eta+4.*eta*eta)
   +(sStar_x*sStar_x + sStar_y*sStar_y + sStar_z*sStar_z)*(2.*eta*eta) );*/
  H += coeffs->dheffSSv2 * eta / (r * r * r * r)
    * (s1Vec->data[0] * s1Vec->data[0] + s1Vec->data[1] * s1Vec->data[1] +
       s1Vec->data[2] * s1Vec->data[2] + s2Vec->data[0] * s2Vec->data[0] +
       s2Vec->data[1] * s2Vec->data[1] + s2Vec->data[2] * s2Vec->data[2]);
  //printf( "Hns = %.16e, Hs = %.16e, Hss = %.16e\n", Hns, Hs, Hss );
  //printf( "H = %.16e\n", H );
  /* Real Hamiltonian given by Eq. 2, ignoring the constant -1. */
  Hreal = sqrt (1. + 2. * eta * (H - 1.));

  return Hreal;
}


/**
 *
 * This function is used to calculate some coefficients which will be used in the
 * spinning EOB Hamiltonian. It takes the following inputs:
 *
 * coeffs - a (non-null) pointer to a SpinEOBParams structure. This will be populated
 * with the output.
 * eta - the symmetric mass ratio.
 * sigmaKerr - the spin of the effective Kerr background (a combination of the individual spins).
 *
 * If all goes well, the function will return XLAL_SUCCESS. Otherwise, XLAL_FAILURE is returned.
 */
static int XLALSimIMRCalculateSpinEOBHCoeffs (SpinEOBHCoeffs * coeffs,
				/**<< OUTPUT, EOB parameters including pre-computed coefficients */
					      const REAL8 eta,
				/**<< symmetric mass ratio */
					      const REAL8 a,
				/**<< Normalized deformed Kerr spin */
					      const UINT4 SpinAlignedEOBversion
					      /**<< 1 for SEOBNRv1; 2 for SEOBNRv2 */
  )
{

  REAL8 KK, k0, k1, k2, k3, k4, k5, k5l, k1p2, k1p3;
  REAL8 m1PlusEtaKK;

  coeffs->SpinAlignedEOBversion = SpinAlignedEOBversion;

  /* Constants are fits taken from Eq. 37 */
  static const REAL8 c0 = 1.4467;	/* needed to get the correct self-force results */
  static const REAL8 c1 = -1.7152360250654402;
  static const REAL8 c2 = -3.246255899738242;

  static const REAL8 c20 = 1.712;
  static const REAL8 c21 = -1.803949138004582;
  static const REAL8 c22 = -39.77229225266885;
  static const REAL8 c23 = 103.16588921239249;

  if (!coeffs)
    {
      XLAL_ERROR (XLAL_EINVAL);
    }


  coeffs->b3 = 0.;
  coeffs->bb3 = 0.;
  coeffs->KK = KK = c0 + c1 * eta + c2 * eta * eta;
  if (SpinAlignedEOBversion == 2)
    {
      coeffs->KK = KK =
	c20 + c21 * eta + c22 * eta * eta + c23 * eta * eta * eta;
    }

  REAL8 chi = a / (1. - 2. * eta);
  REAL8 eta2 = eta * eta, eta3 = eta2 * eta;
  REAL8 chi2 = chi * chi, chi3 = chi2 * chi;

  if (SpinAlignedEOBversion == 4)
    {
      coeffs->KK = KK =
        coeff00K + coeff01K * chi + coeff02K * chi2 + coeff03K * chi3 +
        coeff10K * eta + coeff11K * eta * chi + coeff12K * eta * chi2 +
        coeff13K * eta * chi3 + coeff20K * eta2 + coeff21K * eta2 * chi +
        coeff22K * eta2 * chi2 + coeff23K * eta2 * chi3 + coeff30K * eta3 +
        coeff31K * eta3 * chi + coeff32K * eta3 * chi2 + coeff33K * eta3 * chi3;
//      printf("KK %.16e\n", KK);
    }

  m1PlusEtaKK = -1. + eta * KK;
  /* Eqs. 5.77 - 5.81 of BB1 */
  coeffs->k0 = k0 = KK * (m1PlusEtaKK - 1.);
  coeffs->k1 = k1 = -2. * (k0 + KK) * m1PlusEtaKK;
  k1p2 = k1 * k1;
  k1p3 = k1 * k1p2;
  coeffs->k2 = k2 =
    (k1 * (k1 - 4. * m1PlusEtaKK)) / 2. -
    a * a * k0 * m1PlusEtaKK * m1PlusEtaKK;
  coeffs->k3 = k3 =
    -k1 * k1 * k1 / 3. + k1 * k2 + k1 * k1 * m1PlusEtaKK - 2. * (k2 -
								 m1PlusEtaKK)
    * m1PlusEtaKK - a * a * k1 * m1PlusEtaKK * m1PlusEtaKK;
  coeffs->k4 = k4 =
    (24. * k1 * k1 * k1 * k1 - 96. * k1 * k1 * k2 + 48. * k2 * k2 -
     64. * k1 * k1 * k1 * m1PlusEtaKK + 48. * a * a * (k1 * k1 -
						       2. * k2) *
     m1PlusEtaKK * m1PlusEtaKK + 96. * k1 * (k3 + 2. * k2 * m1PlusEtaKK) -
     m1PlusEtaKK * (192. * k3 +
		    m1PlusEtaKK * (-3008. + 123. * LAL_PI * LAL_PI))) / 96.;
  coeffs->k5 = k5 = 0.0;
  coeffs->k5l = k5l = 0.0;
  if (SpinAlignedEOBversion == 2)
    {
      coeffs->k5 = k5 = m1PlusEtaKK * m1PlusEtaKK
	* (-4237. / 60. + 128. / 5. * LAL_GAMMA +
	   2275. * LAL_PI * LAL_PI / 512. - 1. / 3. * a * a * (k1p3 -
							       3. * k1 * k2 +
							       3. * k3) -
	   (k1p3 * k1p2 - 5. * k1p3 * k2 + 5. * k1 * k2 * k2 +
	    5. * k1p2 * k3 - 5. * k2 * k3 -
	    5. * k1 * k4) / 5. / m1PlusEtaKK / m1PlusEtaKK + (k1p2 * k1p2 -
							      4. * k1p2 * k2 +
							      2. * k2 * k2 +
							      4. * k1 * k3 -
							      4. * k4) / 2. /
	   m1PlusEtaKK + 256. / 5. * log (2.));
      coeffs->k5l = k5l = m1PlusEtaKK * m1PlusEtaKK * 64. / 5.;
    }
  if (SpinAlignedEOBversion == 4)
    {
      /* Include eta^2 terms at 4PN from arXiv:1305.4884 */
      coeffs->k5 = k5 = m1PlusEtaKK * m1PlusEtaKK
	* (-4237. / 60. + 128. / 5. * LAL_GAMMA +
	   2275. * LAL_PI * LAL_PI / 512. - 1. / 3. * a * a * (k1p3 -
							       3. * k1 * k2 +
							       3. * k3) -
	   (k1p3 * k1p2 - 5. * k1p3 * k2 + 5. * k1 * k2 * k2 +
	    5. * k1p2 * k3 - 5. * k2 * k3 -
	    5. * k1 * k4) / 5. / m1PlusEtaKK / m1PlusEtaKK + (k1p2 * k1p2 -
							      4. * k1p2 * k2 +
							      2. * k2 * k2 +
							      4. * k1 * k3 -
							      4. * k4) / 2. /
	   m1PlusEtaKK + 256. / 5. * log (2.) + (41. * LAL_PI * LAL_PI / 32. -
						 221. / 6.) * eta);
      coeffs->k5l = k5l = m1PlusEtaKK * m1PlusEtaKK * 64. / 5.;
    }
  /*printf( "a = %.16e, k0 = %.16e, k1 = %.16e, k2 = %.16e, k3 = %.16e, k4 = %.16e, b3 = %.16e, bb3 = %.16e, KK = %.16e\n",
     a, coeffs->k0, coeffs->k1, coeffs->k2, coeffs->k3, coeffs->k4, coeffs->b3, coeffs->bb3, coeffs->KK );
   */

  /* Now calibrated parameters for spin models */
  coeffs->d1 = coeffs->d1v2 = 0.0;
  coeffs->dheffSS = coeffs->dheffSSv2 = 0.0;
  switch (SpinAlignedEOBversion)
    {
    case 1:
      coeffs->d1 = -69.5;
      coeffs->dheffSS = 2.75;
      break;
    case 2:
      coeffs->d1v2 = -74.71 - 156. * eta + 627.5 * eta * eta;
      coeffs->dheffSSv2 = 8.127 - 154.2 * eta + 830.8 * eta * eta;
      break;
    case 4:
      // dSO
      coeffs->d1v2 =
            coeff00dSO + coeff01dSO * chi + coeff02dSO * chi2 + coeff03dSO * chi3 +
            coeff10dSO * eta + coeff11dSO * eta * chi + coeff12dSO * eta * chi2 +
            coeff13dSO * eta * chi3 + coeff20dSO * eta2 + coeff21dSO * eta2 * chi +
            coeff22dSO * eta2 * chi2 + coeff23dSO * eta2 * chi3 + coeff30dSO * eta3 +
            coeff31dSO * eta3 * chi + coeff32dSO * eta3 * chi2 + coeff33dSO * eta3 * chi3;

      // dSS
      coeffs->dheffSSv2 =
            coeff00dSS + coeff01dSS * chi + coeff02dSS * chi2 + coeff03dSS * chi3 +
            coeff10dSS * eta + coeff11dSS * eta * chi + coeff12dSS * eta * chi2 +
            coeff13dSS * eta * chi3 + coeff20dSS * eta2 + coeff21dSS * eta2 * chi +
            coeff22dSS * eta2 * chi2 + coeff23dSS * eta2 * chi3 + coeff30dSS * eta3 +
            coeff31dSS * eta3 * chi + coeff32dSS * eta3 * chi2 + coeff33dSS * eta3 * chi3;
//          printf("dSO %.16e, dSS %.16e\n", coeffs->d1v2,coeffs->dheffSSv2);
      break;
    default:
      XLALPrintError
	("XLAL Error - %s: wrong SpinAlignedEOBversion value, must be 1 or 2!\n",
	 __func__);
      XLAL_ERROR (XLAL_EINVAL);
      break;
    }

  return XLAL_SUCCESS;
}


/**
 * This function calculates the function \f$\Delta_t(r)\f$ which appears in the spinning EOB
 * potential function. Eqs. 7a and 8.
 */
static REAL8
XLALSimIMRSpinEOBHamiltonianDeltaT (SpinEOBHCoeffs * coeffs,
				/**<< Pre-computed coefficients which appear in the function */
				    const REAL8 r,
				/**<< Current orbital radius (in units of total mass) */
				    const REAL8 eta,
				/**<< Symmetric mass ratio */
				    const REAL8 a
				/**<< Normalized deformed Kerr spin */
  )
{

  REAL8 a2;
  REAL8 u, u2, u3, u4, u5;
  REAL8 m1PlusetaKK;

  REAL8 bulk;
  REAL8 logTerms;
  REAL8 deltaU;
  REAL8 deltaT;

  u = 1. / r;
  u2 = u * u;
  u3 = u2 * u;
  u4 = u2 * u2;
  u5 = u4 * u;

  a2 = a * a;

  m1PlusetaKK = -1. + eta * coeffs->KK;

  bulk = 1. / (m1PlusetaKK * m1PlusetaKK) + (2. * u) / m1PlusetaKK + a2 * u2;

  logTerms =
    1. + eta * coeffs->k0 + eta * log (1. + coeffs->k1 * u + coeffs->k2 * u2 +
				       coeffs->k3 * u3 + coeffs->k4 * u4 +
				       coeffs->k5 * u5 +
				       coeffs->k5l * u5 * log (u));
  /*printf(" a = %.16e, u = %.16e\n",a,u);
     printf( "k0 = %.16e, k1 = %.16e, k2 = %.16e, k3 = %.16e , k4 = %.16e, k5 = %.16e, k5l = %.16e\n",coeffs->k0,
     coeffs->k1,coeffs->k2,coeffs->k3,coeffs->k4,coeffs->k5,coeffs->k5l);
     printf( "bulk = %.16e, logTerms = %.16e\n", bulk, logTerms ); */
  deltaU = bulk * logTerms;

    if ( (coeffs->tidal1->lambda2Tidal != 0. && coeffs->tidal1->omega02Tidal != 0.) || (coeffs->tidal2->lambda2Tidal != 0. && coeffs->tidal2->omega02Tidal != 0.) ) {
        deltaU += XLALSimIMRTEOBdeltaUTidal(u, eta, coeffs->tidal1, coeffs->tidal2);
    }

  deltaT = r * r * deltaU;


  return deltaT;
}


/**
 * This function calculates the function \f$\Delta_r(r)\f$ which appears in the spinning EOB
 * potential function. Eqs. 10a and 10b
 */
static REAL8
XLALSimIMRSpinEOBHamiltonianDeltaR (SpinEOBHCoeffs * coeffs,
				/**<< Pre-computed coefficients which appear in the function */
				    const REAL8 r,
				/**<< Current orbital radius (in units of total mass) */
				    const REAL8 eta,
				/**<< Symmetric mass ratio */
				    const REAL8 a
				/**<< Normalized deformed Kerr spin */
  )
{


  REAL8 u2, u3;
  REAL8 D;
  REAL8 deltaT;			/* The potential function, not a time interval... */
  REAL8 deltaR;

  u2 = 1. / (r * r);
  u3 = u2 / r;

  D = 1. + log (1. + 6. * eta * u2 + 2. * (26. - 3. * eta) * eta * u3);

  deltaT = XLALSimIMRSpinEOBHamiltonianDeltaT (coeffs, r, eta, a);

  deltaR = deltaT * D;
  return deltaR;
}

/**
 * Function to calculate the value of omega for the spin-aligned EOB waveform.
 * Can NOT be used in precessing cases. This omega is defined as \f$\dot{y}/r\f$ by setting \f$y=0\f$.
 * The function calculates omega = v/r, by first converting (r,phi,pr,pphi) to Cartesian coordinates
 * in which rVec={r,0,0} and pVec={0,pphi/r,0}, i.e. the effective-test-particle is positioned at x=r,
 * and its velocity along y-axis. Then it computes omega, which is now given by dydt/r = (dH/dp_y)/r.
 */
static REAL8
XLALSimIMRSpinAlignedEOBCalcOmega (const REAL8 values[],/**<< Dynamical variables */
				   SpinEOBParams * funcParams
							/**<< EOB parameters */
  )
{
  static const REAL8 STEP_SIZE = 1.0e-4;

  HcapDerivParams params;

  /* Cartesian values for calculating the Hamiltonian */
  REAL8 cartValues[6];

  gsl_function F;
  INT4 gslStatus;

  REAL8 omega;
  REAL8 r;

  /* The error in a derivative as measured by GSL */
  REAL8 absErr;

  /* Set up pointers for GSL */
  params.values = cartValues;
  params.params = funcParams;

  F.function = &GSLSpinAlignedHamiltonianWrapper;
  F.params = &params;

  /* Populate the Cartesian values vector */
  /* We can assume phi is zero wlog */
  memset (cartValues, 0, sizeof (cartValues));
  cartValues[0] = r = values[0];
  cartValues[3] = values[2];
  cartValues[4] = values[3] / values[0];

  /* Now calculate omega. In the chosen co-ordinate system, */
  /* we need dH/dpy to calculate this, i.e. varyParam = 4   */
  params.varyParam = 4;
  XLAL_CALLGSL (gslStatus = gsl_deriv_central (&F, cartValues[4],
					       STEP_SIZE, &omega, &absErr));

  if (gslStatus != GSL_SUCCESS)
    {
      XLALPrintError ("XLAL Error - %s: Failure in GSL function\n", __func__);
      XLAL_ERROR_REAL8 (XLAL_EFUNC);
    }

  omega = omega / r;

  return omega;
}

/**
 * Function to calculate the non-Keplerian coefficient for the spin-aligned EOB model.
 * radius \f$r\f$ times the cuberoot of the returned number is \f$r_\Omega\f$ defined in Eq. A2.
 * i.e. the function returns \f$(r_{\Omega} / r)^3\f$.
 */
static REAL8
XLALSimIMRSpinAlignedEOBNonKeplerCoeff (const REAL8 values[],
							/**<< Dynamical variables */
					SpinEOBParams * funcParams
							/**<< EOB parameters */
  )
{

  REAL8 omegaCirc;

  REAL8 tmpValues[4];

  REAL8 r3;

  /* We need to find the values of omega assuming pr = 0 */
  memcpy (tmpValues, values, sizeof (tmpValues));
  tmpValues[2] = 0.0;

  omegaCirc = XLALSimIMRSpinAlignedEOBCalcOmega (tmpValues, funcParams);
  if (XLAL_IS_REAL8_FAIL_NAN (omegaCirc))
    {
      XLAL_ERROR_REAL8 (XLAL_EFUNC);
    }

  r3 = values[0] * values[0] * values[0];

  return 1.0 / (omegaCirc * omegaCirc * r3);
}



/**
 * Function to calculate the non-Keplerian coefficient for the spin-aligned EOB model.
 * radius \f$r\f$ times the cuberoot of the returned number is \f$r_\Omega\f$ defined in Eq. A2.
 * i.e. the function returns \f$(r_{\Omega} / r)^3\f$.
 * This is the generic precessing version
 */
static REAL8 UNUSED
XLALSimIMRSpinEOBNonKeplerCoeff (const REAL8 values[],	/**<< Dynamical variables */
				 SpinEOBParams * funcParams
							/**<< EOB parameters */
  )
{

  REAL8 omegaCirc;

  REAL8 tmpValues[4];

  REAL8 r3;

  /* We need to find the values of omega assuming pr = 0 */
  memcpy (tmpValues, values, sizeof (tmpValues));
  tmpValues[0] =
    sqrt (values[0] * values[0] + values[1] * values[1] +
	  values[2] * values[2]);
  tmpValues[1] = 0.0;
  tmpValues[2] = 0.0;
  tmpValues[3] = sqrt ((values[0] * values[4] - values[1] * values[3])
		       * (values[0] * values[4] - values[1] * values[3])
		       + (values[2] * values[3] - values[0] * values[5])
		       * (values[2] * values[3] - values[0] * values[5])
		       + (values[1] * values[5] - values[2] * values[4])
		       * (values[1] * values[5] - values[2] * values[4]));

  omegaCirc = XLALSimIMRSpinAlignedEOBCalcOmega (tmpValues, funcParams);
  if (XLAL_IS_REAL8_FAIL_NAN (omegaCirc))
    {
      XLAL_ERROR_REAL8 (XLAL_EFUNC);
    }

  r3 = tmpValues[0] * tmpValues[0] * tmpValues[0];

  return 1.0 / (omegaCirc * omegaCirc * r3);
}



/* Wrapper for GSL to call the Hamiltonian function */
static double
GSLSpinAlignedHamiltonianWrapper (double x, void *params)
{
  HcapDerivParams *dParams = (HcapDerivParams *) params;

  EOBParams *eobParams = dParams->params->eobParams;

  REAL8 tmpVec[6];

  /* These are the vectors which will be used in the call to the Hamiltonian */
  REAL8Vector r, p;
  REAL8Vector *s1Vec = dParams->params->s1Vec;
  REAL8Vector *s2Vec = dParams->params->s2Vec;
  REAL8Vector *sigmaKerr = dParams->params->sigmaKerr;
  REAL8Vector *sigmaStar = dParams->params->sigmaStar;

  /* Use a temporary vector to avoid corrupting the main function */
  memcpy (tmpVec, dParams->values, sizeof (tmpVec));

  /* Set the relevant entry in the vector to the correct value */
  tmpVec[dParams->varyParam] = x;

  /* Set the LAL-style vectors to point to the appropriate things */
  r.length = p.length = 3;
  r.data = tmpVec;
  p.data = tmpVec + 3;

  return XLALSimIMRSpinEOBHamiltonian (eobParams->eta, &r, &p, s1Vec, s2Vec,
				       sigmaKerr, sigmaStar,
				       dParams->params->tortoise,
				       dParams->params->seobCoeffs) /
    eobParams->eta;
}

#endif /*_LALSIMIMRSPINEOBHAMILTONIAN_C*/
