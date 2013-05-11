/*
*  Copyright (C) 2010 Craig Robinson, Yi Pan
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
 * \author Craig Robinson, Yi Pan
 *
 * \brief The functions contained within this file pre-compute the various
 * coefficients which are required for calculating the factorized waveform
 * in EOBNRv2. Note that for some of the higher modes, the coefficients
 * are changed for the generation of the waveform compared to the generation
 * of the flux. Thus we have a function which adds these additional 
 * contributions to the already computed coefficients.
 */

#include <math.h>
#include <complex.h>
#include "LALSimIMREOBNRv2.h"

/* Include static functions */
#include "LALSimInspiraldEnergyFlux.c"
#include "LALSimIMREOBNewtonianMultipole.c" 
#include "LALSimIMREOBNQCCorrection.c"

#ifndef _LALSIMIMRFACTORIZEDWAVEFORM_C
#define _LALSIMIMRFACTORIZEDWAVEFORM_C

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/**
 * Constant which comes up in some of the EOB models. Its value is
 * (94/3 -41/32*pi*pi)
 */
#define ninty4by3etc 18.687902694437592603


static inline REAL8 XLALCalculateA5( REAL8 eta );

static inline REAL8 XLALCalculateA6( REAL8 eta );

static
REAL8 XLALCalculateEOBD( REAL8    r,
                         REAL8	eta) UNUSED;


/**
 * Calculates the a5 parameter in the A potential function in EOBNRv2
 */
static inline
REAL8 XLALCalculateA5( const REAL8 eta /**<< Symmetric mass ratio */
                     )
{
  return - 5.82827 - 143.486 * eta + 447.045 * eta * eta;
}

/**
 * Calculates the a6 parameter in the A potential function in EOBNRv2
 */
static inline
REAL8 XLALCalculateA6( const REAL8 UNUSED eta /**<< Symmetric mass ratio */
                     )
{
  return 184.0;
}


/**
 * Function to pre-compute the coefficients in the EOB A potential function
 */
UNUSED static
int XLALCalculateEOBACoefficients(
          EOBACoefficients * const coeffs, /**<< A coefficients (populated in function) */
          const REAL8              eta     /**<< Symmetric mass ratio */
          )
{
  REAL8 eta2, eta3;
  REAL8 a4, a5, a6;

  eta2 = eta*eta;
  eta3 = eta2 * eta;

  /* Note that the definitions of a5 and a6 DO NOT correspond to those in the paper */
  /* Therefore we have to multiply the results of our a5 and a6 finctions by eta. */

  a4 = ninty4by3etc * eta;
  a5 = XLALCalculateA5( eta ) * eta;
  a6 = XLALCalculateA6( eta ) * eta;

  coeffs->n4 =  -64. + 12.*a4 + 4.*a5 + a6 + 64.*eta - 4.*eta2;
  coeffs->n5 = 32. -4.*a4 - a5 - 24.*eta;
  coeffs->d0 = 4.*a4*a4 + 4.*a4*a5 + a5*a5 - a4*a6 + 16.*a6
             + (32.*a4 + 16.*a5 - 8.*a6) * eta + 4.*a4*eta2 + 32.*eta3;
  coeffs->d1 = 4.*a4*a4 + a4*a5 + 16.*a5 + 8.*a6 + (32.*a4 - 2.*a6)*eta + 32.*eta2 + 8.*eta3;
  coeffs->d2 = 16.*a4 + 8.*a5 + 4.*a6 + (8.*a4 + 2.*a5)*eta + 32.*eta2;
  coeffs->d3 = 8.*a4 + 4.*a5 + 2.*a6 + 32.*eta - 8.*eta2;
  coeffs->d4 = 4.*a4 + 2.*a5 + a6 + 16.*eta - 4.*eta2;
  coeffs->d5 = 32. - 4.*a4 - a5 - 24. * eta;

  return XLAL_SUCCESS;
}

/**
 * This function calculates the EOB A function which using the pre-computed
 * coefficients which should already have been calculated.
 */
static
REAL8 XLALCalculateEOBA( const REAL8 r,                     /**<< Orbital separation (in units of total mass M) */
                         EOBACoefficients * restrict coeffs /**<< Pre-computed coefficients for the A function */
                       )
{

  REAL8 r2, r3, r4, r5;
  REAL8 NA, DA;

  /* Note that this function uses pre-computed coefficients,
   * and assumes they have been calculated. Since this is a static function,
   * so only used here, I assume it is okay to neglect error checking
   */

  r2 = r*r;
  r3 = r2 * r;
  r4 = r2*r2;
  r5 = r4*r;


  NA = r4 * coeffs->n4
     + r5 * coeffs->n5;

  DA = coeffs->d0
     + r  * coeffs->d1
     + r2 * coeffs->d2
     + r3 * coeffs->d3
     + r4 * coeffs->d4
     + r5 * coeffs->d5;

  return NA/DA;
}

/**
 * Calculated the derivative of the EOB A function with respect to 
 * r, using the pre-computed A coefficients
 */
static
REAL8 XLALCalculateEOBdAdr( const REAL8 r,                     /**<< Orbital separation (in units of total mass M) */
                            EOBACoefficients * restrict coeffs /**<< Pre-computed coefficients for the A function */
                          )
{
  REAL8 r2, r3, r4, r5;

  REAL8 NA, DA, dNA, dDA, dA;

  r2 = r*r;
  r3 = r2 * r;
  r4 = r2*r2;
  r5 = r4*r;

  NA = r4 * coeffs->n4
     + r5 * coeffs->n5;

  DA = coeffs->d0
     + r  * coeffs->d1
     + r2 * coeffs->d2
     + r3 * coeffs->d3
     + r4 * coeffs->d4
     + r5 * coeffs->d5;

  dNA = 4. * coeffs->n4 * r3
      + 5. * coeffs->n5 * r4;

  dDA = coeffs->d1
      + 2. * coeffs->d2 * r
      + 3. * coeffs->d3 * r2
      + 4. * coeffs->d4 * r3
      + 5. * coeffs->d5 * r4;

  dA = dNA * DA - dDA * NA;

  return dA / (DA*DA);
}

/**
 * Calculate the EOB D function.
 */
static REAL8 XLALCalculateEOBD( REAL8   r, /**<< Orbital separation (in units of total mass M) */
                         REAL8 eta  /**<< Symmetric mass ratio */
                       )
{
	REAL8  u, u2, u3;

	u = 1./r;
	u2 = u*u;
	u3 = u2*u;

	return 1./(1.+6.*eta*u2+2.*eta*(26.-3.*eta)*u3);
}


/**
 * Function to calculate the EOB effective Hamiltonian for the
 * given values of the dynamical variables. The coefficients in the
 * A potential function should already have been computed.
 * Note that the pr used here is the tortoise co-ordinate.
 */
static
REAL8 XLALEffectiveHamiltonian( const REAL8 eta,          /**<< Symmetric mass ratio */
                                const REAL8 r,            /**<< Orbital separation */
                                const REAL8 pr,           /**<< Tortoise co-ordinate */
                                const REAL8 pp,           /**<< Momentum pphi */
                                EOBACoefficients *aCoeffs /**<< Pre-computed coefficients in A function */
                              )
{

        /* The pr used in here is the tortoise co-ordinate */
        REAL8 r2, pr2, pp2, z3, eoba;

        r2   = r * r;
        pr2  = pr * pr;
        pp2  = pp * pp;

        eoba = XLALCalculateEOBA( r, aCoeffs );
        z3   = 2. * ( 4. - 3. * eta ) * eta;
        return sqrt( pr2 + eoba * ( 1.  + pp2/r2 + z3*pr2*pr2/r2 ) );
}


/**
 * Function which calculates the various coefficients used in the generation
 * of the factorized waveform. These coefficients depend only on the symmetric
 * mass ratio eta. It should be noted that this function calculates the 
 * coefficients used in calculating the flux. For generating the waveforms
 * themselves, the coefficients have additional terms added which are calculated
 * using XLALModifyFacWaveformCoefficients(). THe non-spinning parts of these
 * coefficients can be found in Pan et al, arXiv:1106.1021v1 [gr-qc].
 */
UNUSED static int XLALSimIMREOBCalcFacWaveformCoefficients(
          FacWaveformCoeffs * const coeffs, /**<< Structure containing coefficients (populated in function) */
          const REAL8               eta     /**<< Symmetric mass ratio */
          )
{

  REAL8 eta2 = eta*eta;
  REAL8 eta3 = eta2 * eta;

  REAL8 dM, dM2; //dM3;

  REAL8 a = 0;
  REAL8 a2 = 0;
  REAL8 a3 = 0;
  REAL8 chiS = 0;
  REAL8 chiA = 0;

  /* Combination which appears a lot */
  REAL8 m1Plus3eta, m1Plus3eta2, m1Plus3eta3;

  dM2 = 1. - 4.*eta;
  
  /* Check that deltaM has a reasonable value */
  if ( dM2 < 0 )
  {
    XLALPrintError( "eta seems to be < 0.25 - this isn't allowed!\n" );
    XLAL_ERROR( XLAL_EINVAL );
  }

  dM  = sqrt( dM2 );
  //dM3 = dM2 * dM;

  m1Plus3eta  = - 1. + 3.*eta;
  m1Plus3eta2 = m1Plus3eta * m1Plus3eta;
  m1Plus3eta3 = m1Plus3eta * m1Plus3eta2;

  /* Initialize all coefficients to zero */
  /* This is important, as we will not set some if dM is zero */
  memset( coeffs, 0, sizeof( FacWaveformCoeffs ) );


  /* l = 2 */

  coeffs->delta22vh3 = 7./3.;
  coeffs->delta22vh6 = (-4.*a)/3. + (428.*LAL_PI)/105.;
  coeffs->delta22v8 = (20.*a)/63.;
  coeffs->delta22vh9 = -2203./81. + (1712.*LAL_PI*LAL_PI)/315.;
  coeffs->delta22v5  = - 24.*eta;

  coeffs->rho22v2   = -43./42. + (55.*eta)/84.;
  coeffs->rho22v3   = (-2.*(chiS + chiA*dM - chiS*eta))/3.;
  coeffs->rho22v4   = -20555./10584. + (chiS*chiS + 2.*chiA*chiS*dM + chiA*chiA*dM2)/2.
       - (33025.*eta)/21168. + (19583.*eta2)/42336.;
  coeffs->rho22v5   = (-34.*a)/21.;
  coeffs->rho22v6   = 1556919113./122245200. + (89.*a2)/252. - (48993925.*eta)/9779616. 
       - (6292061.*eta2)/3259872. + (10620745.*eta3)/39118464.
       + (41.*eta*LAL_PI*LAL_PI)/192.;
  coeffs->rho22v6l  = - 428./105.;
  coeffs->rho22v7   = (18733.*a)/15876. + a*a2/3.;
  coeffs->rho22v8   = -387216563023./160190110080. + (18353.*a2)/21168. - a2*a2/8.;
  coeffs->rho22v8l  =  9202./2205.;
  coeffs->rho22v10  = -16094530514677./533967033600.;
  coeffs->rho22v10l =  439877./55566.;


  if ( dM2 )
  {
    coeffs->delta21vh3 = 2./3.;
    coeffs->delta21vh6 = (-17.*a)/35. + (107.*LAL_PI)/105.;
    coeffs->delta21vh7 = (3.*a2)/140.;
    coeffs->delta21vh9 = -272./81. + (214.*LAL_PI*LAL_PI)/315.;
    coeffs->delta21v5  = - 493. * eta /42.;

    coeffs->rho21v1   = (-3.*(chiS+chiA/dM))/(4.);
    //coeffs->rho21v2   = -59./56 - (9.*chiAPlusChiSdM*chiAPlusChiSdM)/(32.*dM2) + (23.*eta)/84.;
    /*coeffs->rho21v3   = (-567.*chiA*chiA*chiA - 1701.*chiA*chiA*chiS*dM
                        + chiA*(-4708. + 1701.*chiS*chiS - 2648.*eta)*(-1. + 4.*eta)
                        + chiS* dM3 *(4708. - 567.*chiS*chiS
                        + 1816.*eta))/(2688.*dM3);*/
    coeffs->rho21v2   = -59./56. + (23.*eta)/84. - 9./32.*a2;
    coeffs->rho21v3   = 1177./672.*a - 27./128.*a3;
    coeffs->rho21v4   = -47009./56448.- (865.*a2)/1792. - (405.*a2*a2)/2048. - (10993.*eta)/14112.
                        + (617.*eta2)/4704.;
    coeffs->rho21v5   = (-98635.*a)/75264. + (2031.*a*a2)/7168. - (1701.*a2*a3)/8192.;
    coeffs->rho21v6   = 7613184941./2607897600.+ (9032393.*a2)/1806336. + (3897.*a2*a2)/16384.
                        - (15309.*a3*a3)/65536.; 
    coeffs->rho21v6l  = - 107./105.;
    coeffs->rho21v7   = (-3859374457.*a)/1159065600. - (55169.*a3)/16384.
                        + (18603.*a2*a3)/65536. - (72171.*a2*a2*a3)/262144.;
    coeffs->rho21v7l  =  107.*a/140.;
    coeffs->rho21v8   = -1168617463883./911303737344.;
    coeffs->rho21v8l  = 6313./5880.;
    coeffs->rho21v10  = -63735873771463./16569158860800.; 
    coeffs->rho21v10l = 5029963./5927040.;
  }

  /* l = 3 */
  if ( dM2 )
  {
    coeffs->delta33vh3 = 13./10.;
    coeffs->delta33vh6 = (-81.*a)/20. + (39.*LAL_PI)/7.;
    coeffs->delta33vh9 = -227827./3000. + (78.*LAL_PI*LAL_PI)/7.;
    coeffs->delta33v5  = - 80897.*eta / 2430.;

    coeffs->rho33v2 = -7./6. + (2.*eta)/3.;
    coeffs->rho33v3 = (chiS*dM*(-4. + 5.*eta) + chiA*(-4. + 19.*eta))/(6.*dM);
    coeffs->rho33v4 = -6719./3960. + a2/2. - (1861.*eta)/990. + (149.*eta2)/330.;
    coeffs->rho33v5 = (-4.*a)/3.;
    coeffs->rho33v6 = 3203101567./227026800. + (5.*a2)/36.;
    coeffs->rho33v6l = - 26./7.;
    coeffs->rho33v7 = (5297.*a)/2970. + a*a2/3.;
    coeffs->rho33v8 = -57566572157./8562153600.;
    coeffs->rho33v8l = 13./3.;
  }

  coeffs->delta32vh3 = (10. + 33.*eta)/(-15.*m1Plus3eta);
  coeffs->delta32vh4 = 4.*a;
  coeffs->delta32vh6 = (-136.*a)/45. + (52.*LAL_PI)/21.;
  coeffs->delta32vh9 = -9112./405. + (208.*LAL_PI*LAL_PI)/63.;

  coeffs->rho32v   = (4.*chiS*eta)/(-3.*m1Plus3eta);
  coeffs->rho32v2  = (-4.*a2*eta2)/(9.*m1Plus3eta2) + (328. - 1115.*eta
                        + 320.*eta2)/(270.*m1Plus3eta);
  coeffs->rho32v3  = (2.*(45.*a*m1Plus3eta3
                        - a*eta*(328. - 2099.*eta + 5.*(733. + 20.*a2)*eta2
                        - 960.*eta3)))/(405.*m1Plus3eta3);
  coeffs->rho32v4  = a2/3. + (-1444528.
                        + 8050045.*eta - 4725605.*eta2 - 20338960.*eta3
                        + 3085640.*eta2*eta2)/(1603800.*m1Plus3eta2);
  coeffs->rho32v5  = (-2788.*a)/1215.;
  coeffs->rho32v6  = 5849948554./940355325. + (488.*a2)/405.;
  coeffs->rho32v6l =  - 104./63.;
  coeffs->rho32v8  = -10607269449358./3072140846775.;
  coeffs->rho32v8l = 17056./8505.;

  if ( dM2 )
  {
    coeffs->delta31vh3 = 13./30.;
    coeffs->delta31vh6 = (61.*a)/20. + (13.*LAL_PI)/21.;
    coeffs->delta31vh7 = (-24.*a2)/5.;
    coeffs->delta31vh9 = -227827./81000. + (26.*LAL_PI*LAL_PI)/63.;
    coeffs->delta31v5  = - 17.*eta/10.; 
 
    coeffs->rho31v2  = -13./18. - (2.*eta)/9.;
    coeffs->rho31v3  = (chiA*(-4. + 11.*eta) + chiS*dM*(-4. + 13.*eta))/(6.*dM);
    coeffs->rho31v4  = 101./7128.
                        - (5.*a2)/6. - (1685.*eta)/1782. - (829.*eta2)/1782.;
    coeffs->rho31v5  = (4.*a)/9.;
    coeffs->rho31v6  = 11706720301./6129723600. - (49.*a2)/108.;
    coeffs->rho31v6l =  - 26./63.;
    coeffs->rho31v7  = (-2579.*a)/5346. + a*a2/9.;
    coeffs->rho31v8  = 2606097992581./4854741091200.;
    coeffs->rho31v8l = 169./567.;
  }

  /* l = 4 */
  
  coeffs->delta44vh3 = (112. + 219.*eta)/(-120.*m1Plus3eta);
  coeffs->delta44vh6 = (-464.*a)/75. + (25136.*LAL_PI)/3465.;

  coeffs->rho44v2 = (1614. - 5870.*eta + 2625.*eta2)/(1320.*m1Plus3eta);
  coeffs->rho44v3 = (chiA*(10. - 39.*eta)*dM + chiS*(10. - 41.*eta
                        + 42.*eta2))/(15.*m1Plus3eta);
  coeffs->rho44v4 = a2/2. + (-511573572.
                        + 2338945704.*eta - 313857376.*eta2 - 6733146000.*eta3
                        + 1252563795.*eta2*eta2)/(317116800.*m1Plus3eta2);
  coeffs->rho44v5 = (-69.*a)/55.;
  coeffs->rho44v6 = 16600939332793./1098809712000. + (217.*a2)/3960.;
  coeffs->rho44v6l = - 12568./3465.;

  if ( dM2 )
  {
    coeffs->delta43vh3 = (486. + 4961.*eta)/(810.*(1. - 2.*eta));
    coeffs->delta43vh4 = (11.*a)/4.;
    coeffs->delta43vh6 = 1571.*LAL_PI/385.;

    coeffs->rho43v   = (5.*(chiA - chiS*dM)*eta)/(8.*dM*(-1. + 2.*eta));
    coeffs->rho43v2  = (222. - 547.*eta + 160.*eta2)/(176.*(-1. + 2.*eta));
    coeffs->rho43v4  = -6894273./7047040. + (3.*a2)/8.;
    coeffs->rho43v5  = (-12113.*a)/6160.;
    coeffs->rho43v6  = 1664224207351./195343948800.;
    coeffs->rho43v6l = - 1571./770.;
  }

  coeffs->delta42vh3 = (7.*(1. + 6.*eta))/(-15.*m1Plus3eta);
  coeffs->delta42vh6 = (212.*a)/75. + (6284.*LAL_PI)/3465.;

  coeffs->rho42v2  = (1146. - 3530.*eta + 285.*eta2)/(1320.*m1Plus3eta);
  coeffs->rho42v3  = (chiA*(10. - 21.*eta)*dM + chiS*(10. - 59.*eta
                        + 78.*eta2))/(15.*m1Plus3eta);
  coeffs->rho42v4  = a2/2. + (-114859044. + 295834536.*eta + 1204388696.*eta2 - 3047981160.*eta3
                        - 379526805.*eta2*eta2)/(317116800.*m1Plus3eta2);
  coeffs->rho42v5  = (-7.*a)/110.;
  coeffs->rho42v6  = 848238724511./219761942400. + (2323.*a2)/3960.;
  coeffs->rho42v6l = - 3142./3465.;

  if ( dM2 )
  {
    coeffs->delta41vh3 = (2. + 507.*eta)/(10.*(1. - 2.*eta));
    coeffs->delta41vh4 = (11.*a)/12.;
    coeffs->delta41vh6 = 1571.*LAL_PI/3465.;

    coeffs->rho41v   = (5.*(chiA - chiS*dM)*eta)/(8.*dM*(-1. + 2.*eta));
    coeffs->rho41v2  = (602. - 1385.*eta + 288.*eta2)/(528.*(-1. + 2.*eta));
    coeffs->rho41v4  = -7775491./21141120. + (3.*a2)/8.;
    coeffs->rho41v5  = (-20033.*a)/55440. - (5*a*a2)/6.;
    coeffs->rho41v6  = 1227423222031./1758095539200.;
    coeffs->rho41v6l = - 1571./6930.;
  }

  /* l = 5 */
  if ( dM2 )
  {
    coeffs->delta55vh3 = (96875. + 857528.*eta)/(131250.*(1 - 2*eta));

    coeffs->rho55v2 = (487. - 1298.*eta + 512.*eta2)/(390.*(-1. + 2.*eta));
    coeffs->rho55v3 = (-2.*a)/3.;
    coeffs->rho55v4 = -3353747./2129400. + a2/2.;
    coeffs->rho55v5 = - 241. * a / 195.;
  }

  coeffs->delta54vh3 = 8./15.;
  coeffs->delta54vh4 = 12.*a/5.;

  coeffs->rho54v2 = (-17448. + 96019.*eta - 127610.*eta2
                        + 33320.*eta3)/(13650.*(1. - 5.*eta + 5.*eta2));
  coeffs->rho54v3 = (-2.*a)/15.;
  coeffs->rho54v4 = -16213384./15526875. + (2.*a2)/5.;

  if ( dM2 )
  {
    coeffs->delta53vh3 = 31./70.;

    coeffs->rho53v2 = (375. - 850.*eta + 176.*eta2)/(390.*(-1. + 2.*eta));
    coeffs->rho53v3 = (-2.*a)/3.;
    coeffs->rho53v4 = -410833./709800. + a2/2.;
    coeffs->rho53v5 = - 103.*a/325.;
  }

  coeffs->delta52vh3 = 4./15.;
  coeffs->delta52vh4 = 6.*a/5.;

  coeffs->rho52v2 = (-15828. + 84679.*eta - 104930.*eta2
                        + 21980.*eta3)/(13650.*(1. - 5.*eta + 5.*eta2));
  coeffs->rho52v3 = (-2.*a)/15.;
  coeffs->rho52v4 = -7187914./15526875. + (2.*a2)/5.;

  if ( dM2 )
  {
    coeffs->delta51vh3 = 31./210.;

    coeffs->rho51v2 = (319. - 626.*eta + 8.*eta2)/(390.*(-1. + 2.*eta));
    coeffs->rho51v3 = (-2.*a)/3.;
    coeffs->rho51v4 = -31877./304200. + a2/2.;
    coeffs->rho51v5 = 139.*a/975.;
  }

  /* l = 6 */

  coeffs->delta66vh3 = 43./70.;
  
  coeffs->rho66v2 = (-106. + 602.*eta - 861.*eta2
                        + 273.*eta3)/(84.*(1. - 5.*eta + 5.*eta2));
  coeffs->rho66v3 = (-2.*a)/3.;
  coeffs->rho66v4 = -1025435./659736. + a2/2.;

  if ( dM2 )
  {
    coeffs->delta65vh3 = 10./21.;
    
    coeffs->rho65v2 = (-185. + 838.*eta - 910.*eta2
                        + 220.*eta3)/(144.*(dM2 + 3.*eta2));
    coeffs->rho65v3 = - 2.*a/9.;
  }

  coeffs->delta64vh3 = 43./105.;
  
  coeffs->rho64v2 = (-86. + 462.*eta - 581.*eta2
                        + 133.*eta3)/(84.*(1. - 5.*eta + 5.*eta2));
  coeffs->rho64v3 = (-2.*a)/3.;
  coeffs->rho64v4 = -476887./659736. + a2/2.;

  if ( dM2 )
  {
    coeffs->delta63vh3 = 2./7.;

    coeffs->rho63v2 = (-169. + 742.*eta - 750.*eta2
                        + 156.*eta3)/(144.*(dM2 + 3.*eta2));
    coeffs->rho63v3 = - 2.*a/9.;
  }

  coeffs->delta62vh3 = 43./210.;

  coeffs->rho62v2 = (-74. + 378.*eta - 413.*eta2
                        + 49.*eta3)/(84.*(1. - 5.*eta + 5.*eta2));
  coeffs->rho62v3 = (-2.*a)/3.;
  coeffs->rho62v4 = -817991./3298680. + a2/2.;

  if ( dM2 )
  {
    coeffs->delta61vh3 = 2./21.;

    coeffs->rho61v2 = (-161. + 694.*eta - 670.*eta2
                        + 124.*eta3)/(144.*(dM2 + 3.*eta2));
    coeffs->rho61v3 = - 2. * a / 9.;
  }

  /* l = 7 */
  if ( dM2 )
  {
    coeffs->delta77vh3 = 19./36.;

    coeffs->rho77v2 = (-906. + 4246.*eta - 4963.*eta2
                        + 1380.*eta3)/(714.*(dM2 + 3.*eta2));
    coeffs->rho77v3 = - 2.*a/3.;
  }

  coeffs->rho76v2 = (2144. - 16185.*eta + 37828.*eta2 - 29351.*eta3
                        + 6104.*eta2*eta2) / (1666.*(-1 + 7*eta - 14*eta2
                        + 7*eta3));

  if ( dM2 )
  {
    coeffs->delta75vh3 = 95./252.;

    coeffs->rho75v2 = (-762. + 3382.*eta - 3523.*eta2
                        + 804.*eta3)/(714.*(dM2 + 3.*eta2));
    coeffs->rho75v3 = - 2.*a/3.;
  }

  coeffs->rho74v2 = (17756. - 131805.*eta + 298872.*eta2 - 217959.*eta3
                        + 41076.*eta2*eta2) / (14994.*(-1. + 7.*eta - 14.*eta2
                        + 7.*eta3));

  if ( dM2 )
  {
    coeffs->delta73vh3 = 19./84.;

    coeffs->rho73v2 = (-666. + 2806.*eta - 2563.*eta2
                        + 420.*eta3)/(714.*(dM2 + 3.*eta2));
    coeffs->rho73v3 = - 2.*a/3.;
  }

  coeffs->rho72v2 = (16832. - 123489.*eta + 273924.*eta2 - 190239.*eta3
                        + 32760.*eta2*eta2) /(14994.*(-1. + 7.*eta - 14.*eta2
                        + 7.*eta3));

  if ( dM2 )
  {
    coeffs->delta71vh3 = 19./252.;

    coeffs->rho71v2 = (-618. + 2518.*eta - 2083.*eta2
                        + 228.*eta3)/(714.*(dM2 + 3.*eta2));
    coeffs->rho71v3 = - 2.*a/3.;
  }

  /* l = 8 */
  
  coeffs->rho88v2 = (3482. - 26778.*eta + 64659.*eta2 - 53445.*eta3
                        + 12243.*eta2*eta2) / (2736.*(-1. + 7.*eta - 14.*eta2
                        + 7.*eta3));

  if ( dM2 )
  {
    coeffs->rho87v2 = (23478. - 154099.*eta + 309498.*eta2 - 207550.*eta3
                        + 38920*eta2*eta2) / (18240.*(-1 + 6*eta - 10*eta2
                        + 4*eta3));
  }

  coeffs->rho86v2 = (1002. - 7498.*eta + 17269.*eta2 - 13055.*eta3
                        + 2653.*eta2*eta2) / (912.*(-1. + 7.*eta - 14.*eta2
                        + 7.*eta3));

  if ( dM2 )
  {
    coeffs->rho85v2 = (4350. - 28055.*eta + 54642.*eta2 - 34598.*eta3
                        + 6056.*eta2*eta2) / (3648.*(-1. + 6.*eta - 10.*eta2
                        + 4.*eta3));
  }

  coeffs->rho84v2 = (2666. - 19434.*eta + 42627.*eta2 - 28965.*eta3
                        + 4899.*eta2*eta2) / (2736.*(-1. + 7.*eta - 14.*eta2
                        + 7.*eta3));

  if ( dM2 )
  {
    coeffs->rho83v2 = (20598. - 131059.*eta + 249018.*eta2 - 149950.*eta3
                        + 24520.*eta2*eta2) / (18240.*(-1. + 6.*eta - 10.*eta2
                        + 4.*eta3));
  }

  coeffs->rho82v2 = (2462. - 17598.*eta + 37119.*eta2 - 22845.*eta3
                        + 3063.*eta2*eta2) / (2736.*(-1. + 7.*eta - 14.*eta2
                        + 7.*eta3));

  if ( dM2 )
  {
    coeffs->rho81v2 = (20022. - 126451.*eta + 236922.*eta2 - 138430.*eta3
                        + 21640.*eta2*eta2) / (18240.*(-1. + 6.*eta - 10.*eta2
                        + 4.*eta3));
  }

  /* All relevant coefficients should be set, so we return */

  return XLAL_SUCCESS;
}


/**
 * Function which adds the additional terms required for waveform generation
 * to the factorized waveform coefficients. Note that this function only calculates
 * additional terms not present in the flux, so the factorized waveform coefficients
 * SHOULD ALREADY HAVE BEEN CALCULATED using XLALCalcFacWaveformCoefficients() prior
 * to calling this function.
 */
UNUSED static int XLALSimIMREOBModifyFacWaveformCoefficients( 
                                       FacWaveformCoeffs * const coeffs, /**<< Structure containing coefficients */
                                       const REAL8 eta                   /**<< Symmetric mass ratio */
                                     )
{

  if ( !coeffs )
  {
    XLAL_ERROR( XLAL_EINVAL );
  }

  /* Tweak the relevant coefficients for the generation of the waveform */
  coeffs->rho21v6 += -5. * eta;
  coeffs->rho33v6 += -20. * eta;
  coeffs->rho44v6 += -15. * eta;
  coeffs->rho55v6 += 4. * eta;

  coeffs->delta21v7 += 30. * eta;
  coeffs->delta33v7 += -10. * eta;
  coeffs->delta44v5 += -70. * eta;
  coeffs->delta55v5 += 40. * eta;

  return XLAL_SUCCESS;
}

/**
 * Computes the non-Keplerian correction to the velocity as determined from the
 * frequency obtained assuming a circular orbit. In the early stages of the evolution,
 * this should be a number close to 1.
 */
static REAL8
nonKeplerianCoefficient(
                   REAL8Vector * restrict values, /**<< Dynamics r, phi, pr, pphi */
                   const REAL8       eta,         /**<< Symmetric mass ratio */
                   EOBACoefficients *coeffs       /**<< Pre-computed A coefficients */
                   )
{

  REAL8 r    = values->data[0];
  REAL8 pphi = values->data[3];

  REAL8 A  = XLALCalculateEOBA( r, coeffs );
  REAL8 dA = XLALCalculateEOBdAdr( r, coeffs );

  return 2. * (1. + 2. * eta * ( -1. + sqrt( (1. + pphi*pphi/(r*r)) * A ) ) )
          / ( r*r * dA );
}

/**
 * Computes the factorized waveform according to the prescription
 * given in Pan et al, arXiv:1106.1021v1 [gr-qc], for a given
 * mode l,m, for the given values of the dynamics at that point.
 * The function returns XLAL_SUCCESS if everything works out properly,
 * otherwise XLAL_FAILURE will be returned.
 */
UNUSED static int  XLALSimIMREOBGetFactorizedWaveform( 
                                COMPLEX16   * restrict hlm,    /**<< The value of hlm (populated by the function) */
                                REAL8Vector * restrict values, /**<< Vector containing dynamics r, phi, pr, pphi for a given point */
                                const REAL8 v,                 /**<< Velocity (in geometric units) */
                                const INT4  l,                 /**<< Mode l */
                                const INT4  m,                 /**<< Mode m */
                                EOBParams   * restrict params  /**<< Structure containing pre-computed coefficients, etc. */
                                )
{

  /* Status of function calls */
  INT4 status;
  INT4 i;

  REAL8 eta;
  REAL8 r, pr, pp, Omega, v2, vh, vh3, k, hathatk, eulerlogxabs;
  REAL8 Hreal, Heff, Slm, deltalm, rholm, rholmPwrl;
  COMPLEX16 Tlm;
  COMPLEX16 hNewton;
  gsl_sf_result lnr1, arg1, z2;

  /* Non-Keplerian velocity */
  REAL8 vPhi;

  /* Pre-computed coefficients */
  FacWaveformCoeffs *hCoeffs = params->hCoeffs;

  if ( abs(m) > (INT4) l )
  {
    XLAL_ERROR( XLAL_EINVAL );
  }


  eta = params->eta;

  /* Check our eta was sensible */
  if ( eta > 0.25 )
  {
    XLALPrintError("Eta seems to be > 0.25 - this isn't allowed!\n" );
    XLAL_ERROR( XLAL_EINVAL );
  }
  else if ( eta == 0.25 && m % 2 )
  {
    /* If m is odd and dM = 0, hLM will be zero */
    memset( hlm, 0, sizeof( COMPLEX16 ) );
    return XLAL_SUCCESS;
  }

  r  = values->data[0];
  pr = values->data[2];
  pp = values->data[3];

  Heff  = XLALEffectiveHamiltonian( eta, r, pr, pp, params->aCoeffs );
  Hreal = sqrt( 1.0 + 2.0 * eta * ( Heff - 1.0) );
  v2    = v * v;
  Omega = v2 * v;
  vh3   = Hreal * Omega;
  vh    = cbrt(vh3);
  eulerlogxabs = LAL_GAMMA + log( 2.0 * (REAL8)m * v );


  /* Calculate the non-Keplerian velocity */
  /* given by Eq. (18) of Pan et al, PRD84, 124052(2011) */
  /* psi given by Eq. (19) of Pan et al, PRD84, 124052(2011) */
  /* Assign temporarily to vPhi */
  vPhi = nonKeplerianCoefficient( values, eta, params->aCoeffs );
  /* Assign rOmega value temporarily to vPhi */
  vPhi  = r * cbrt(vPhi);
  /* Assign rOmega * Omega to vPhi */
  vPhi *= Omega;

  /* Calculate the newtonian multipole */
  status = XLALSimIMREOBCalculateNewtonianMultipole( &hNewton, vPhi * vPhi, vPhi/Omega,
            values->data[1], (UINT4)l, m, params );
  if ( status == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Calculate the source term */
  if ( ( (l+m)%2 ) == 0)
  {
    Slm = Heff;
  }
  else
  {
    Slm = v * pp;
  }

  /* Calculate the Tail term */
  k  = m * Omega;
  hathatk = Hreal * k;
  XLAL_CALLGSL( status = gsl_sf_lngamma_complex_e( l+1.0, -2.0*hathatk, &lnr1, &arg1 ) );
  if (status != GSL_SUCCESS)
  {
    XLALPrintError("Error in GSL function\n" );
    XLAL_ERROR( XLAL_EFUNC );
  }
  XLAL_CALLGSL( status = gsl_sf_fact_e( l, &z2 ) );
  if ( status != GSL_SUCCESS)
  {
    XLALPrintError("Error in GSL function\n" );
    XLAL_ERROR( XLAL_EFUNC );
  }
  Tlm = cexp( ( lnr1.val + LAL_PI * hathatk ) + I * (
        arg1.val + 2.0 * hathatk * log(4.0*k/sqrt(LAL_E)) ) );
  Tlm /= z2.val;

  /* Calculate the residue phase and amplitude terms */
  switch( l )
  {
    case 2:
      switch( abs(m) )
      {
        case 2:
          deltalm = vh3*(hCoeffs->delta22vh3 + vh3*(hCoeffs->delta22vh6
            + vh*vh*(hCoeffs->delta22vh9*vh)))
            + hCoeffs->delta22v5 *v*v2*v2 + hCoeffs->delta22v8 *v2*v2*v2*v2;
          rholm  = 1. + v2*(hCoeffs->rho22v2 + v*(hCoeffs->rho22v3
            + v*(hCoeffs->rho22v4
            + v*(hCoeffs->rho22v5 + v*(hCoeffs->rho22v6
            + hCoeffs->rho22v6l*eulerlogxabs + v*(hCoeffs->rho22v7
            + v*(hCoeffs->rho22v8 + hCoeffs->rho22v8l*eulerlogxabs
            + (hCoeffs->rho22v10 + hCoeffs->rho22v10l * eulerlogxabs)*v2)))))));
          break;
        case 1:
          deltalm = vh3*(hCoeffs->delta21vh3 + vh3*(hCoeffs->delta21vh6
            + vh*(hCoeffs->delta21vh7 + (hCoeffs->delta21vh9)*vh*vh)))
            + hCoeffs->delta21v5*v*v2*v2 + hCoeffs->delta21v7*v2*v2*v2*v;
          rholm  = 1. + v*(hCoeffs->rho21v1
            + v*( hCoeffs->rho21v2 + v*(hCoeffs->rho21v3 + v*(hCoeffs->rho21v4
            + v*(hCoeffs->rho21v5 + v*(hCoeffs->rho21v6 + hCoeffs->rho21v6l*eulerlogxabs
            + v*(hCoeffs->rho21v7 + hCoeffs->rho21v7l * eulerlogxabs
            + v*(hCoeffs->rho21v8 + hCoeffs->rho21v8l * eulerlogxabs
            + (hCoeffs->rho21v10 + hCoeffs->rho21v10l * eulerlogxabs)*v2))))))));
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    case 3:
      switch (m)
      {
        case 3:
          deltalm = vh3*(hCoeffs->delta33vh3 + vh3*(hCoeffs->delta33vh6 + hCoeffs->delta33vh9*vh3))
            + hCoeffs->delta33v5*v*v2*v2 + hCoeffs->delta33v7*v2*v2*v2*v;
          rholm  = 1. + v2*(hCoeffs->rho33v2 + v*(hCoeffs->rho33v3 + v*(hCoeffs->rho33v4
            + v*(hCoeffs->rho33v5 + v*(hCoeffs->rho33v6 + hCoeffs->rho33v6l*eulerlogxabs
            + v*(hCoeffs->rho33v7 + (hCoeffs->rho33v8 + hCoeffs->rho33v8l*eulerlogxabs)*v))))));
          break;
        case 2:
          deltalm = vh3*(hCoeffs->delta32vh3 + vh*(hCoeffs->delta32vh4 + vh*vh*(hCoeffs->delta32vh6
            + hCoeffs->delta32vh9*vh3)));
          rholm  = 1. + v*(hCoeffs->rho32v
            + v*(hCoeffs->rho32v2 + v*(hCoeffs->rho32v3 + v*(hCoeffs->rho32v4 + v*(hCoeffs->rho32v5
            + v*(hCoeffs->rho32v6 + hCoeffs->rho32v6l*eulerlogxabs
            + (hCoeffs->rho32v8 + hCoeffs->rho32v8l*eulerlogxabs)*v2))))));
          break;
        case 1:
          deltalm = vh3*(hCoeffs->delta31vh3 + vh3*(hCoeffs->delta31vh6
            + vh*(hCoeffs->delta31vh7 + hCoeffs->delta31vh9*vh*vh)))
            + hCoeffs->delta31v5*v*v2*v2;
          rholm  = 1. + v2*(hCoeffs->rho31v2 + v*(hCoeffs->rho31v3 + v*(hCoeffs->rho31v4
            + v*(hCoeffs->rho31v5 + v*(hCoeffs->rho31v6 + hCoeffs->rho31v6l*eulerlogxabs
            + v*(hCoeffs->rho31v7 + (hCoeffs->rho31v8 + hCoeffs->rho31v8l*eulerlogxabs)*v))))));
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    case 4:
      switch (m)
      {
        case 4:
          deltalm = vh3*(hCoeffs->delta44vh3 + hCoeffs->delta44vh6 *vh3)
            + hCoeffs->delta44v5*v2*v2*v;
          rholm  = 1. + v2*(hCoeffs->rho44v2
            + v*( hCoeffs->rho44v3 + v*(hCoeffs->rho44v4
            + v*(hCoeffs->rho44v5 + (hCoeffs->rho44v6
            + hCoeffs->rho44v6l*eulerlogxabs)*v))));
          break;
        case 3:
          deltalm = vh3*(hCoeffs->delta43vh3 + vh*(hCoeffs->delta43vh4
            + hCoeffs->delta43vh6*vh*vh));
          rholm  = 1. + v*(hCoeffs->rho43v
            + v*(hCoeffs->rho43v2
            + v2*(hCoeffs->rho43v4 + v*(hCoeffs->rho43v5
            + (hCoeffs->rho43v6 + hCoeffs->rho43v6l*eulerlogxabs)*v))));
          break;
        case 2:
          deltalm = vh3*(hCoeffs->delta42vh3 + hCoeffs->delta42vh6*vh3);
          rholm  = 1. + v2*(hCoeffs->rho42v2
            + v*(hCoeffs->rho42v3 + v*(hCoeffs->rho42v4 + v*(hCoeffs->rho42v5
            + (hCoeffs->rho42v6 + hCoeffs->rho42v6l*eulerlogxabs)*v))));
          break;
        case 1:
          deltalm = vh3*(hCoeffs->delta41vh3 + vh*(hCoeffs->delta41vh4
            + hCoeffs->delta41vh6*vh*vh));
          rholm  = 1. + v*(hCoeffs->rho41v
            + v*(hCoeffs->rho41v2
            + v2*(hCoeffs->rho41v4 + v*(hCoeffs->rho41v5
            + (hCoeffs->rho41v6 +  hCoeffs->rho41v6l*eulerlogxabs)*v))));
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    case 5:
      switch (m)
      {
        case 5:
          deltalm = hCoeffs->delta55vh3*vh3 + hCoeffs->delta55v5*v2*v2*v;
          rholm  = 1. + v2*( hCoeffs->rho55v2
            + v*(hCoeffs->rho55v3 + v*(hCoeffs->rho55v4
            + v*(hCoeffs->rho55v5 + hCoeffs->rho55v6*v))));
          break;
        case 4:
          deltalm = vh3*(hCoeffs->delta54vh3 + hCoeffs->delta54vh4*vh);
          rholm  = 1. + v2*(hCoeffs->rho54v2 + v*(hCoeffs->rho54v3
            + hCoeffs->rho54v4*v));
          break;
        case 3:
          deltalm = hCoeffs->delta53vh3 * vh3;
          rholm  = 1. + v2*(hCoeffs->rho53v2
            + v*(hCoeffs->rho53v3 + v*(hCoeffs->rho53v4 + hCoeffs->rho53v5*v)));
          break;
        case 2:
          deltalm = vh3*(hCoeffs->delta52vh3 + hCoeffs->delta52vh4*vh);
          rholm  = 1. + v2*(hCoeffs->rho52v2 + v*(hCoeffs->rho52v3
            + hCoeffs->rho52v4*v));
          break;
        case 1:
          deltalm = hCoeffs->delta51vh3 * vh3;
          rholm  = 1. + v2*(hCoeffs->rho51v2
            + v*(hCoeffs->rho51v3 + v*(hCoeffs->rho51v4 + hCoeffs->rho51v5*v)));
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    case 6:
      switch (m)
      {
        case 6:
          deltalm = hCoeffs->delta66vh3*vh3;
          rholm  = 1. + v2*(hCoeffs->rho66v2 + v*(hCoeffs->rho66v3
            + hCoeffs->rho66v4*v));
          break;
        case 5:
          deltalm = hCoeffs->delta65vh3*vh3;
          rholm  = 1. + v2*(hCoeffs->rho65v2 + hCoeffs->rho65v3*v);
          break;
        case 4:
          deltalm = hCoeffs->delta64vh3 * vh3;
          rholm  = 1. + v2*(hCoeffs->rho64v2 + v*(hCoeffs->rho64v3
            + hCoeffs->rho64v4*v));
          break;
        case 3:
          deltalm = hCoeffs->delta63vh3 * vh3;
          rholm  = 1. + v2*(hCoeffs->rho63v2 + hCoeffs->rho63v3*v);
          break;
        case 2:
          deltalm = hCoeffs->delta62vh3 * vh3;
          rholm  = 1. + v2*(hCoeffs->rho62v2 + v*(hCoeffs->rho62v3
            + hCoeffs->rho62v4 * v));
          break;
        case 1:
          deltalm = hCoeffs->delta61vh3 * vh3;
          rholm  = 1. + v2*(hCoeffs->rho61v2 + hCoeffs->rho61v3*v);
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    case 7:
      switch (m)
      {
        case 7:
          deltalm = hCoeffs->delta77vh3 * vh3;
          rholm   = 1. + v2*(hCoeffs->rho77v2 + hCoeffs->rho77v3 * v);
          break;
        case 6:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho76v2 * v2;
          break;
        case 5:
          deltalm = hCoeffs->delta75vh3 * vh3;
          rholm   = 1. + v2*(hCoeffs->rho75v2 + hCoeffs->rho75v3*v);
          break;
        case 4:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho74v2 * v2;
          break;
        case 3:
          deltalm = hCoeffs->delta73vh3 *vh3;
          rholm   = 1. + v2*(hCoeffs->rho73v2 + hCoeffs->rho73v3 * v);
          break;
        case 2:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho72v2 * v2;
          break;
        case 1:
          deltalm = hCoeffs->delta71vh3 * vh3;
          rholm   = 1. + v2*(hCoeffs->rho71v2 +hCoeffs->rho71v3 * v);
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    case 8:
      switch (m)
      {
        case 8:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho88v2 * v2;
          break;
        case 7:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho87v2 * v2;
          break;
        case 6:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho86v2 * v2;
          break;
        case 5:
          deltalm = 0.0;
          rholm   = 1. + hCoeffs->rho85v2 * v2;
          break;
        case 4:
          deltalm = 0.0;
          rholm  = 1. + hCoeffs->rho84v2 * v2;
          break;
        case 3:
          deltalm = 0.0;
          rholm  = 1. + hCoeffs->rho83v2 * v2;
          break;
        case 2:
          deltalm = 0.0;
          rholm  = 1. + hCoeffs->rho82v2 * v2;
          break;
        case 1:
          deltalm = 0.0;
          rholm  = 1. + hCoeffs->rho81v2 * v2;
          break;
        default:
          XLAL_ERROR( XLAL_EINVAL );
          break;
      }
      break;
    default:
      XLAL_ERROR( XLAL_EINVAL );
      break;
  }

  /* Raise rholm to the lth power */
  rholmPwrl = 1.0;
  i = l;
  while ( i-- )
  {
    rholmPwrl *= rholm;
  }

  *hlm = Tlm * cexp(I * deltalm) * Slm * rholmPwrl;
  *hlm *= hNewton;

  return XLAL_SUCCESS;
} 

#endif /*_LALSIMIMRFACTORIZEDWAVEFORM_C*/
