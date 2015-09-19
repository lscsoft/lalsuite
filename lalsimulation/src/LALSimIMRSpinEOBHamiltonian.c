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
 * \author Craig Robinson, Yi Pan
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
#include <gsl/gsl_deriv.h>

#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

#include "LALSimIMRSpinEOB.h"

#include "LALSimIMRSpinEOBAuxFuncs.c"
#include "LALSimIMRSpinEOBFactorizedWaveformCoefficients.c"
/*#include "LALSimIMRSpinEOBFactorizedWaveform.c"*/

/*------------------------------------------------------------------------------------------
 *
 *          Prototypes of functions defined in this code.
 *
 *------------------------------------------------------------------------------------------
 */

/**
 * This function calculates the DeltaR potential function in the spin EOB Hamiltonian
 */
static REAL8 XLALSimIMRSpinEOBHamiltonianDeltaR(
        SpinEOBHCoeffs *coeffs, /**<< Pre-computed coefficients which appear in the function */
        const REAL8    r,       /**<< Current orbital radius (in units of total mass) */
        const REAL8    eta,     /**<< Symmetric mass ratio */
        const REAL8    a        /**<< Normalized deformed Kerr spin */
        );

static REAL8 XLALSimIMRSpinEOBHamiltonian(
               const REAL8    eta,
               REAL8Vector    * restrict x,
               REAL8Vector    * restrict p,
               REAL8Vector    * restrict s1Vec,
               REAL8Vector    * restrict s2Vec,
               REAL8Vector    * restrict sigmaKerr,
               REAL8Vector    * restrict sigmaStar,
               int                       tortoise,
               SpinEOBHCoeffs *coeffs);

static int XLALSimIMRCalculateSpinEOBHCoeffs(
        SpinEOBHCoeffs *coeffs,
        const REAL8    eta,
        const REAL8    a,
        const UINT4    SpinAlignedEOBversion
        );

static REAL8 XLALSimIMRSpinEOBHamiltonianDeltaT( 
        SpinEOBHCoeffs *coeffs,
        const REAL8    r,
        const REAL8    eta,
        const REAL8    a
        );

static REAL8 XLALSimIMRSpinAlignedEOBCalcOmega(
                      const REAL8          values[],
                      SpinEOBParams        *funcParams
                      );

static REAL8 XLALSimIMRSpinAlignedEOBNonKeplerCoeff(
                      const REAL8           values[],
                      SpinEOBParams         *funcParams
                      );

static double GSLSpinAlignedHamiltonianWrapper( double x, void *params );

/* Precessing EOB's function declarations below */
UNUSED static REAL8 inner_product( const REAL8 values1[], 
                             const REAL8 values2[]
                             );

UNUSED static REAL8* cross_product( const REAL8 values1[],
                              const REAL8 values2[],
                              REAL8 result[] );

static REAL8 UNUSED XLALSimIMRSpinEOBNonKeplerCoeff(
                      const REAL8           values[],   /**<< Dynamical variables */
                      SpinEOBParams         *funcParams /**<< EOB parameters */
                      );

static REAL8 XLALSimIMRSpinEOBCalcOmega(
                      const REAL8           values[],   /**<< Dynamical variables */
                      SpinEOBParams         *funcParams /**<< EOB parameters */
                      );

UNUSED static int XLALSpinHcapRvecDerivative(
                 double UNUSED     t,         /**<< UNUSED */
                 const  REAL8      values[],  /**<< Dynamical variables */
                 REAL8             dvalues[], /**<< Time derivatives of variables (returned) */
                 void             *funcParams /**<< EOB parameters */
                               );

static double GSLSpinHamiltonianWrapperForRvecDerivs( double x, void *params );

static double GSLSpinHamiltonianWrapperFordHdpphi( double x, void *params );


/*------------------------------------------------------------------------------------------
 *
 *          Defintions of functions.
 *
 *------------------------------------------------------------------------------------------
 */

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
static REAL8 XLALSimIMRSpinEOBHamiltonian( 
               const REAL8    eta,                  /**<< Symmetric mass ratio */
               REAL8Vector    * restrict x,         /**<< Position vector */
               REAL8Vector    * restrict p,	    /**<< Momentum vector (tortoise radial component pr*) */
               REAL8Vector    * restrict s1Vec,     /**<< Spin vector 1 */
               REAL8Vector    * restrict s2Vec,     /**<< Spin vector 2 */
               REAL8Vector    * restrict sigmaKerr, /**<< Spin vector sigma_kerr */
               REAL8Vector    * restrict sigmaStar, /**<< Spin vector sigma_star */
               INT4                      tortoise,  /**<< flag to state whether the momentum is the tortoise co-ord */
	       SpinEOBHCoeffs *coeffs               /**<< Structure containing various coefficients */
               )
{
//    if(isnan(x->data[0]) || isnan(x->data[1]) || isnan(x->data[2])) {x->data[0] = 0.1;x->data[1] = 0.1;x->data[2] = 0.1;}
//    if(isnan(p->data[0]) || isnan(p->data[1]) || isnan(p->data[2])) {p->data[0] = 0.1;p->data[1] = 0.1;p->data[2] = 0.1;}
  int debugPK = 0;
  /* Update the Hamiltonian coefficients, if spins are evolving */
  int UsePrec = 1;
  if ( UsePrec && coeffs->updateHCoeffs )
  {
    SpinEOBHCoeffs tmpCoeffs; 
    REAL8 tmpa;
    
    tmpa = sqrt(sigmaKerr->data[0]*sigmaKerr->data[0]
                + sigmaKerr->data[1]*sigmaKerr->data[1]
                + sigmaKerr->data[2]*sigmaKerr->data[2]);

    if ( XLALSimIMRCalculateSpinEOBHCoeffs( &tmpCoeffs, eta, 
          tmpa, coeffs->SpinAlignedEOBversion ) == XLAL_FAILURE )
    {
      XLAL_ERROR( XLAL_EFUNC );
    }

    tmpCoeffs.SpinAlignedEOBversion = coeffs->SpinAlignedEOBversion;
    tmpCoeffs.updateHCoeffs = coeffs->updateHCoeffs;

    coeffs = &tmpCoeffs;
  }

  REAL8 r, r2, nx, ny, nz;
  REAL8 sKerr_x, sKerr_y, sKerr_z, a, a2;
  REAL8 sStar_x, sStar_y, sStar_z;
  REAL8 e3_x, e3_y, e3_z;
  REAL8 costheta; /* Cosine of angle between Skerr and r */
  REAL8 xi2, xi_x, xi_y, xi_z; /* Cross product of unit vectors in direction of Skerr and r */
  REAL8 vx, vy, vz, pxir, pvr, pn, prT, pr, pf, ptheta2; /*prT is the tortoise pr */
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
    REAL8 tmpP[3] = {0.};

  REAL8 csi;

  /* Spin gauge parameters. (YP) simplified, since both are zero. */
  // static const double aa=0., bb=0.;
  if(debugPK){
  printf( "In Hamiltonian: tortoise flag = %d\n", (int) tortoise );
  printf( "x = %.16e\t%.16e\t%.16e\n", x->data[0], x->data[1], x->data[2] );
  printf( "p = %.16e\t%.16e\t%.16e\n", p->data[0], p->data[1], p->data[2] );
  printf( "sStar = %.16e\t%.16e\t%.16e\n", sigmaStar->data[0], 
		sigmaStar->data[1], sigmaStar->data[2] );
  printf( "sKerr = %.16e\t%.16e\t%.16e\n", sigmaKerr->data[0], 
		sigmaKerr->data[1], sigmaKerr->data[2] );}
  

  r2 = x->data[0]*x->data[0] + x->data[1]*x->data[1] + x->data[2]*x->data[2];
  r  = sqrt(r2);
  nx = x->data[0] / r;
  ny = x->data[1] / r;
  nz = x->data[2] / r;   
     
  sKerr_x = sigmaKerr->data[0];
  sKerr_y = sigmaKerr->data[1];
  sKerr_z = sigmaKerr->data[2];

  sStar_x = sigmaStar->data[0];
  sStar_y = sigmaStar->data[1];
  sStar_z = sigmaStar->data[2];
     
  a2 = sKerr_x*sKerr_x + sKerr_y*sKerr_y + sKerr_z*sKerr_z;
  a  = sqrt( a2 );

  if(a != 0.) 
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
       
  costheta = e3_x*nx + e3_y*ny + e3_z*nz; 
    
  xi2=1. - costheta*costheta; 

  xi_x = -e3_z*ny + e3_y*nz;
  xi_y =  e3_z*nx - e3_x*nz;
  xi_z = -e3_y*nx + e3_x*ny;

  vx = -nz*xi_y + ny*xi_z;
  vy =  nz*xi_x - nx*xi_z;
  vz = -ny*xi_x + nx*xi_y;

  w2 = r2 + a2;
  rho2 = r2 + a2*costheta*costheta;

  u  = 1./r;
  u2 = u*u;
  u3 = u2*u;
  u4 = u2*u2;
  u5 = u4*u;

  if(debugPK)printf( "KK = %.16e\n", coeffs->KK );
  m1PlusetaKK = -1. + eta * coeffs->KK;
  /* Eq. 5.75 of BB1 */
  bulk = 1./(m1PlusetaKK*m1PlusetaKK) + (2.*u)/m1PlusetaKK + a2*u2;
  /* Eq. 5.73 of BB1 */
  logTerms = 1. + eta*coeffs->k0 + eta*log(fabs(1. + coeffs->k1*u + coeffs->k2*u2 + coeffs->k3*u3 + coeffs->k4*u4
                                              + coeffs->k5*u5 + coeffs->k5l*u5*log(u)));
  if(debugPK)printf( "bulk = %.16e, logTerms = %.16e\n", bulk, logTerms );
  /* Eq. 5.73 of BB1 */
  deltaU = bulk*logTerms;
    deltaU = fabs(deltaU);
  /* Eq. 5.71 of BB1 */
  deltaT = r2*deltaU;
  /* ddeltaU/du */
  deltaU_u = 2.*(1./m1PlusetaKK + a2*u)*logTerms + 
	  bulk * (eta*(coeffs->k1 + u*(2.*coeffs->k2 + u*(3.*coeffs->k3 + u*(4.*coeffs->k4 + 5.*(coeffs->k5+coeffs->k5l*log(u))*u)))))
          / (1. + coeffs->k1*u + coeffs->k2*u2 + coeffs->k3*u3 + coeffs->k4*u4 + (coeffs->k5+coeffs->k5l*log(u))*u5);
  /* ddeltaT/dr */
  deltaT_r = 2.*r*deltaU - deltaU_u;
  /* Eq. 5.39 of BB1 */
  Lambda = fabs(w2*w2 - a2*deltaT*xi2);
  /* Eq. 5.83 of BB1, inverse */
  D = 1. + log(1. + 6.*eta*u2 + 2.*(26. - 3.*eta)*eta*u3);
  /* Eq. 5.38 of BB1 */
  deltaR = deltaT*D;
  /* See Hns below, Eq. 4.34 of Damour et al. PRD 62, 084011 (2000) */
  qq = 2.*eta*(4. - 3.*eta);
  /* See Hns below. In Sec. II D of BB2 b3 and bb3 coeffs are chosen to be zero. */
  ww=2.*a*r + coeffs->b3*eta*a2*a*u + coeffs->bb3*eta*a*u;

  /* We need to transform the momentum to get the tortoise co-ord */
  if ( tortoise )
  {
    csi = sqrt( fabs(deltaT * deltaR) )/ w2; /* Eq. 28 of Pan et al. PRD 81, 084041 (2010) */
  }
  else
  {
    csi = 1.0;
  }
  if(debugPK)printf( "csi(miami) = %.16e\n", csi );


  if ( tortoise != 2 )
  {
	  prT = p->data[0]*nx + p->data[1]*ny + p->data[2]*nz;
      /* p->data is BL momentum vector; tmpP is tortoise momentum vector */ 
      tmpP[0] = p->data[0] - nx * prT * (csi - 1.)/csi;
      tmpP[1] = p->data[1] - ny * prT * (csi - 1.)/csi;
      tmpP[2] = p->data[2] - nz * prT * (csi - 1.)/csi;
  }
  else
  {
	  prT = (p->data[0]*nx + p->data[1]*ny + p->data[2]*nz)*csi;
      /* p->data is BL momentum vector; tmpP is tortoise momentum vector */ 
      tmpP[0] = p->data[0];
      tmpP[1] = p->data[1];
      tmpP[2] = p->data[2];
  }
  
  pxir = (tmpP[0]*xi_x + tmpP[1]*xi_y + tmpP[2]*xi_z) * r;
  pvr  = (tmpP[0]*vx + tmpP[1]*vy + tmpP[2]*vz) * r;
  pn   = tmpP[0]*nx + tmpP[1]*ny + tmpP[2]*nz;
          
  pr = pn;
  pf = pxir;
  ptheta2 = pvr * pvr / xi2;

  if(debugPK)
  {printf( "pr = %.16e, prT = %.16e\n", pr, prT );

  printf( " a = %.16e, r = %.16e\n", a, r );
  printf( "D = %.16e, ww = %.16e, rho = %.16e, Lambda = %.16e, xi = %.16e\npr = %.16e, pf = %.16e, deltaR = %.16e, deltaT = %.16e\n", 
      D, ww, sqrt(rho2), Lambda, sqrt(xi2), pr, pf, deltaR, deltaT );}

  /* Eqs. 5.36 - 5.46 of BB1 */
  /* Note that the tortoise prT appears only in the quartic term, explained in Eqs. 14 and 15 of Tarrachini et al. */
  Hns = sqrt(1. + prT*prT*prT*prT*qq*u2 + ptheta2/rho2 + pf*pf*rho2/(Lambda*xi2) + pr*pr*deltaR/rho2)
      / sqrt(Lambda/(rho2*deltaT)) + pf*ww/Lambda;
  
  if(debugPK){
  printf( "term 1 in Hns: %.16e\n",  prT*prT*prT*prT*qq*u2 );
  printf( "term 2 in Hns: %.16e\n", ptheta2/rho2 );
  printf( "term 3 in Hns = %.16e\n", pf*pf*rho2/(Lambda*xi2) );
  printf( "term 4 in Hns = %.16e\n", pr*pr*deltaR/rho2 );
  printf( "term 5 in Hns = %.16e\n", Lambda/(rho2*deltaT) );
  printf( "term 6 in Hns = %.16e\n", pf*ww/Lambda );}
  /* Eqs. 5.30 - 5.33 of BB1 */
  B = sqrt(deltaT);
  w = ww/Lambda;
  nu = 0.5 * log(deltaT*rho2/Lambda);
  MU = 0.5 * log(rho2);  
  /* dLambda/dr */
  Lambda_r = 4.*r*w2 - a2*deltaT_r*xi2;
     
  ww_r=2.*a - (a2*a*coeffs->b3*eta)*u2 - coeffs->bb3*eta*a*u2;
  /* Eqs. 5.47a - 5.47d of BB1 */
  BR = (-2.*deltaT + sqrt(deltaR)*deltaT_r)/(2.*sqrt(deltaR*deltaT));
  wr = (-Lambda_r*ww + Lambda*ww_r)/(Lambda*Lambda);
  nur = (r/rho2 + (w2 * (-4.*r*deltaT + w2*deltaT_r) ) / (2.*deltaT*Lambda) );
  mur = (r/rho2 - 1./sqrt(deltaR));
  /* Eqs. 5.47f - 5.47h of BB1 */
  wcos  = -2.*a2*costheta*deltaT*ww/(Lambda*Lambda);  
  nucos = a2*costheta*w2*(w2-deltaT)/(rho2*Lambda);  
  mucos = a2*costheta/rho2;
  /* Eq. 5.52 of BB1, (YP) simplified */
  //Q = 1. + pvr*pvr/(exp(2.*MU)*xi2) + exp(2.*nu)*pxir*pxir/(B*B*xi2) + pn*pn*deltaR/exp(2.*MU);
  Q = 1. + pvr*pvr/(rho2*xi2) + deltaT*rho2/Lambda*pxir*pxir/(B*B*xi2) + pn*pn*deltaR/rho2;
   if(debugPK){
       printf( "Q = %.16e, pvr = %.16e, xi2 = %.16e , deltaT = %.16e, rho2 = %.16e, Lambda = %.16e, pxir = %.16e, B = %.16e\n", Q, pvr, xi2, deltaT, rho2, Lambda, pxir, B );
   }
  pn2 = pr * pr * deltaR / rho2;
  pp  = Q - 1.;

  if(debugPK){printf( "pn2 = %.16e, pp = %.16e\n", pn2, pp );
  printf( "sigmaKerr = %.16e, sigmaStar = %.16e\n", sKerr_z, sStar_z );}
  /* Eq. 5.68 of BB1, (YP) simplified for aa=bb=0. */
  /*
  deltaSigmaStar_x=(- 8.*aa*(1. + 3.*pn2*r - pp*r)*sKerr_x - 8.*bb*(1. + 3.*pn2*r - pp*r)*sStar_x + 
        eta*(-8.*sKerr_x - 36.*pn2*r*sKerr_x + 3.*pp*r*sKerr_x + 14.*sStar_x - 30.*pn2*r*sStar_x + 4.*pp*r*sStar_x))/(12.*r);

  deltaSigmaStar_y=(-8.*aa*(1. + 3.*pn2*r - pp*r)*sKerr_y - 8.*bb*(1. + 3.*pn2*r - pp*r)*sStar_y + 
        eta*(-8.*sKerr_y - 36.*pn2*r*sKerr_y + 3.*pp*r*sKerr_y + 14.*sStar_y - 30.*pn2*r*sStar_y + 4.*pp*r*sStar_y))/(12.*r);

  deltaSigmaStar_z=(-8.*aa*(1. + 3.*pn2*r - pp*r)*sKerr_z - 8.*bb*(1. + 3.*pn2*r - pp*r)*sStar_z + 
	eta*(-8.*sKerr_z - 36.*pn2*r*sKerr_z + 3.*pp*r*sKerr_z + 14.*sStar_z - 30.*pn2*r*sStar_z + 4.*pp*r*sStar_z))/(12.*r);
  */
  deltaSigmaStar_x=eta*(-8.*sKerr_x - 36.*pn2*r*sKerr_x + 3.*pp*r*sKerr_x + 14.*sStar_x - 30.*pn2*r*sStar_x + 4.*pp*r*sStar_x)/(12.*r);

  deltaSigmaStar_y=eta*(-8.*sKerr_y - 36.*pn2*r*sKerr_y + 3.*pp*r*sKerr_y + 14.*sStar_y - 30.*pn2*r*sStar_y + 4.*pp*r*sStar_y)/(12.*r);

  deltaSigmaStar_z=eta*(-8.*sKerr_z - 36.*pn2*r*sKerr_z + 3.*pp*r*sKerr_z + 14.*sStar_z - 30.*pn2*r*sStar_z + 4.*pp*r*sStar_z)/(12.*r);


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
  sMultiplier1 = -(2.*eta*(-353. + 27.*eta) + 2.*(103.*eta - 60.*eta*eta)*pp*r 
               + 120.*(-3.*eta*eta)*pn2*pn2*r*r + (eta*(23. + 3.*eta))*pp*pp*r*r 
               + 6.*pn2*r*(- 47.*eta + 54.*eta*eta + (- 16.*eta + 21.*eta*eta)*pp*r))
               / (72.*r*r);                        
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
  sMultiplier2 = (-16.*(7.*eta*(8. + 3.*eta)) + 4.*(- 109.*eta + 51.*eta*eta)*pp*r 
               + 810.*eta*eta*pn2*pn2*r*r - 45.*eta*pp*pp*r*r 
               - 6.*pn2*r*(16.*eta + 147.*eta*eta + (- 6.*eta + 39.*eta*eta)*pp*r))
               / (144.*r*r);
  /* Eq. 52 of BB2 */                     
  deltaSigmaStar_x += sMultiplier1*sigmaStar->data[0] + sMultiplier2*sigmaKerr->data[0];
  deltaSigmaStar_y += sMultiplier1*sigmaStar->data[1] + sMultiplier2*sigmaKerr->data[1];
  deltaSigmaStar_z += sMultiplier1*sigmaStar->data[2] + sMultiplier2*sigmaKerr->data[2];

  /* And now the (calibrated) 4.5PN term */
  deltaSigmaStar_x += coeffs->d1 * eta * sigmaStar->data[0] / (r*r*r);
  deltaSigmaStar_y += coeffs->d1 * eta * sigmaStar->data[1] / (r*r*r);
  deltaSigmaStar_z += coeffs->d1 * eta * sigmaStar->data[2] / (r*r*r);
  deltaSigmaStar_x += coeffs->d1v2 * eta * sigmaKerr->data[0] / (r*r*r);
  deltaSigmaStar_y += coeffs->d1v2 * eta * sigmaKerr->data[1] / (r*r*r);
  deltaSigmaStar_z += coeffs->d1v2 * eta * sigmaKerr->data[2] / (r*r*r);


  if(debugPK)printf( "deltaSigmaStar_x = %.16e, deltaSigmaStar_y = %.16e, deltaSigmaStar_z = %.16e\n", 
     deltaSigmaStar_x, deltaSigmaStar_y, deltaSigmaStar_z );
	
  sx = sStar_x + deltaSigmaStar_x;
  sy = sStar_y + deltaSigmaStar_y;
  sz = sStar_z + deltaSigmaStar_z;     
     
     
  sxi = sx*xi_x + sy*xi_y + sz*xi_z;
  sv  = sx*vx + sy*vy + sz*vz;
  sn  = sx*nx + sy*ny + sz*nz; 
     
  s3 = sx*e3_x + sy*e3_y + sz*e3_z;  
  /* Eq. 3.45 of BB1, second term */        
  Hwr = (exp(-3.*MU - nu)*sqrt(deltaR)*(exp(2.*(MU + nu))*pxir*pxir*sv - B*exp(MU + nu)*pvr*pxir*sxi + 
        B*B*xi2*(exp(2.*MU)*(sqrt(Q) + Q)*sv + pn*pvr*sn*sqrt(deltaR) - pn*pn*sv*deltaR)))/(2.*B*(1. + sqrt(Q))*sqrt(Q)*xi2);
  /* Eq. 3.45 of BB1, third term */     
  Hwcos = (exp(-3.*MU - nu)*(sn*(-(exp(2.*(MU + nu))*pxir*pxir) + B*B*(pvr*pvr - exp(2.*MU)*(sqrt(Q) + Q)*xi2)) - 
        B*pn*(B*pvr*sv - exp(MU + nu)*pxir*sxi)*sqrt(deltaR)))/(2.*B*(1. + sqrt(Q))*sqrt(Q));
  /* Eq. 3.44 of BB1, leading term */     
  HSOL = (exp(-MU + 2.*nu)*(-B + exp(MU + nu))*pxir*s3)/(B*B*sqrt(Q)*xi2);
  /* Eq. 3.44 of BB1, next-to-leading term */
  HSONL = (exp(-2.*MU + nu)*(-(B*exp(MU + nu)*nucos*pxir*(1. + 2.*sqrt(Q))*sn*xi2) + 
        (-(BR*exp(MU + nu)*pxir*(1. + sqrt(Q))*sv) + B*(exp(MU + nu)*nur*pxir*(1. + 2.*sqrt(Q))*sv + B*mur*pvr*sxi + 
        B*sxi*(-(mucos*pn*xi2) + sqrt(Q)*(mur*pvr - nur*pvr + (-mucos + nucos)*pn*xi2))))*sqrt(deltaR)))/(B*B*(sqrt(Q) + Q)*xi2);   
  /* Eq. 3.43 and 3.45 of BB1 */
  Hs = w*s3 + Hwr*wr + Hwcos*wcos + HSOL + HSONL;
  /* Eq. 5.70 of BB1, last term */   
  Hss = -0.5*u3 * (sx*sx + sy*sy + sz*sz - 3.*sn*sn);
  /* Eq. 5.70 of BB1 */
  H = Hns + Hs + Hss;

  /* Add the additional calibrated term */
  H += coeffs->dheffSS * eta * (sKerr_x*sStar_x + sKerr_y*sStar_y + sKerr_z*sStar_z) / (r*r*r*r);
  /* One more calibrated term proportional to S1^2+S2^2. Note that we use symmetric expressions of m1,m2 and S1,S2 */
  /*H += coeffs->dheffSSv2 * eta / (r*r*r*r) / (1.-4.*eta)
                         * ( (sKerr_x*sKerr_x + sKerr_y*sKerr_y + sKerr_z*sKerr_z)*(1.-4.*eta+2.*eta*eta)
                            +(sKerr_x*sStar_x + sKerr_y*sStar_y + sKerr_z*sStar_z)*(-2.*eta+4.*eta*eta)
                            +(sStar_x*sStar_x + sStar_y*sStar_y + sStar_z*sStar_z)*(2.*eta*eta) );*/
  H += coeffs->dheffSSv2 * eta / (r*r*r*r)
                         * (s1Vec->data[0]*s1Vec->data[0] + s1Vec->data[1]*s1Vec->data[1] + s1Vec->data[2]*s1Vec->data[2]
                           +s2Vec->data[0]*s2Vec->data[0] + s2Vec->data[1]*s2Vec->data[1] + s2Vec->data[2]*s2Vec->data[2]);
  if(debugPK){
	  printf( "Hns = %.16e, Hs = %.16e, Hss = %.16e\n", Hns, Hs, Hss );
	  printf( "H = %.16e\n", H );}
  /* Real Hamiltonian given by Eq. 2, ignoring the constant -1. */
  Hreal = sqrt(1. + 2.*eta *(fabs(H) - 1.));
  if(debugPK)printf( "Hreal = %.16e\n", Hreal );
  
  if(isnan(Hreal)) {
    printf(
    "\n\nInside Hamiltonian: Hreal is a NAN. Printing its components below:\n");
      printf( "(deltaU, bulk, logTerms, log arg) = (%.16e, %.16e, %.16e, %.16e)\n", deltaU, bulk, logTerms, 1. + coeffs->k1*u + coeffs->k2*u2 + coeffs->k3*u3 + coeffs->k4*u4
             + coeffs->k5*u5 + coeffs->k5l*u5*log(u));

    printf( "In Hamiltonian: tortoise flag = %d\n", (int) tortoise );
    printf( "x = %.16e\t%.16e\t%.16e\n", x->data[0], x->data[1], x->data[2] );
    printf( "p = %.16e\t%.16e\t%.16e\n", p->data[0], p->data[1], p->data[2] );
    printf( "sStar = %.16e\t%.16e\t%.16e\n", sigmaStar->data[0], 
      sigmaStar->data[1], sigmaStar->data[2] );
    printf( "sKerr = %.16e\t%.16e\t%.16e\n", sigmaKerr->data[0], 
      sigmaKerr->data[1], sigmaKerr->data[2] );
      printf("csi = %.16e, Q = %.16e, pvr = %.16e, xi2 = %.16e , deltaT = %.16e, rho2 = %.16e, Lambda = %.16e, pxir = %.16e, B = %.16e\n", csi,Q, pvr, xi2, deltaT, rho2, Lambda, pxir, B );

    printf( "KK = %.16e\n", coeffs->KK );
    printf( "bulk = %.16e, logTerms = %.16e\n", bulk, logTerms );
    printf("csi(miami) = %.16e\n", csi); 
    printf( " a = %.16e, r = %.16e\n", a, r );
    printf( "D = %.16e, ww = %.16e, rho = %.16e, Lambda = %.16e, xi = %.16e\npr = %.16e, pf = %.16e, deltaR = %.16e, deltaT = %.16e\n", 
        D, ww, sqrt(rho2), Lambda, sqrt(xi2), pr, pf, deltaR, deltaT );
    printf( "pr = %.16e, prT = %.16e\n", pr, prT );
    
    printf( " a = %.16e, r = %.16e\n", a, r );
    printf( "D = %.16e, ww = %.16e, rho = %.16e, Lambda = %.16e, xi = %.16e\npr = %.16e, pf = %.16e, deltaR = %.16e, deltaT = %.16e\n", 
        D, ww, sqrt(rho2), Lambda, sqrt(xi2), pr, pf, deltaR, deltaT );
    printf( "pr = %.16e, prT = %.16e\n", pr, prT );
    printf( "pn2 = %.16e, pp = %.16e\n", pn2, pp );
    printf( "deltaSigmaStar_x = %.16e, deltaSigmaStar_y = %.16e, deltaSigmaStar_z = %.16e\n", 
     deltaSigmaStar_x, deltaSigmaStar_y, deltaSigmaStar_z );

    printf( "term 1 in Hns: %.16e\n",  prT*prT*prT*prT*qq*u2 );
    printf( "term 2 in Hns: %.16e\n", ptheta2/rho2 );
    printf( "term 3 in Hns = %.16e\n", pf*pf*rho2/(Lambda*xi2) );
    printf( "term 4 in Hns = %.16e\n", pr*pr*deltaR/rho2 );
    printf( "term 5 in Hns = %.16e\n", Lambda/(rho2*deltaT) );
    printf( "term 6 in Hns = %.16e\n", pf*ww/Lambda );

    printf( "Hns = %.16e, Hs = %.16e, Hss = %.16e\n", Hns, Hs, Hss );
	  printf( "H = %.16e\n", H );

    printf("Done printing components.\n\n");
    XLALPrintError( "XLAL Error - %s: Hreal = nan in Hamiltonian \n", __func__);
    XLAL_ERROR( XLAL_EINVAL );
  }

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
static int XLALSimIMRCalculateSpinEOBHCoeffs(
        SpinEOBHCoeffs *coeffs, /**<< OUTPUT, EOB parameters including pre-computed coefficients */
        const REAL8    eta,     /**<< symmetric mass ratio */
        const REAL8    a,       /**<< Normalized deformed Kerr spin */
        const UINT4    SpinAlignedEOBversion  /**<< 1 for SEOBNRv1; 2 for SEOBNRv2 */
        )
{

  REAL8 KK, k0, k1, k2, k3, k4, k5, k5l, k1p2, k1p3;
  REAL8 m1PlusEtaKK;

  coeffs->SpinAlignedEOBversion = SpinAlignedEOBversion;
   
  int debugPK = 0;
  if( debugPK )
  {
    printf("In XLALSimIMRCalculateSpinEOBHCoeffs: SpinAlignedEOBversion = %d,%d\n",
        (int) SpinAlignedEOBversion, (int) coeffs->SpinAlignedEOBversion );
    fflush( NULL );
  }
  /* Constants are fits taken from Eq. 37 */
  static const REAL8 c0  = 1.4467; /* needed to get the correct self-force results */
  static const REAL8 c1  = -1.7152360250654402;
  static const REAL8 c2  = -3.246255899738242;

  static const REAL8 c20  = 1.712;
  static const REAL8 c21  = -1.803949138004582;
  static const REAL8 c22  = -39.77229225266885;
  static const REAL8 c23  = 103.16588921239249;

  if ( !coeffs )
  {
    XLAL_ERROR( XLAL_EINVAL );
  }


  coeffs->b3  = 0.;
  coeffs->bb3 = 0.;
  coeffs->KK = KK  = c0 + c1*eta + c2*eta*eta;
  if ( SpinAlignedEOBversion == 2)
  {
     coeffs->KK = KK = c20 + c21*eta + c22*eta*eta + c23*eta*eta*eta;
  }

  m1PlusEtaKK = -1. + eta*KK;
  /* Eqs. 5.77 - 5.81 of BB1 */
  coeffs->k0 = k0 = KK*(m1PlusEtaKK - 1.);
  coeffs->k1 = k1 = - 2.*(k0 + KK)*m1PlusEtaKK;
  k1p2= k1*k1;
  k1p3= k1*k1p2;
  coeffs->k2 = k2 = (k1 * (k1 - 4.*m1PlusEtaKK)) / 2. - a*a*k0*m1PlusEtaKK*m1PlusEtaKK;
  coeffs->k3 = k3 = -k1*k1*k1/3. + k1*k2 + k1*k1*m1PlusEtaKK - 2.*(k2 - m1PlusEtaKK)*m1PlusEtaKK - a*a*k1*m1PlusEtaKK*m1PlusEtaKK;
  coeffs->k4 = k4 = (24.*k1*k1*k1*k1 - 96.*k1*k1*k2 + 48.*k2*k2 - 64.*k1*k1*k1*m1PlusEtaKK
      + 48.*a*a*(k1*k1 - 2.*k2)*m1PlusEtaKK*m1PlusEtaKK +
      96.*k1*(k3 + 2.*k2*m1PlusEtaKK) - m1PlusEtaKK*(192.*k3 + m1PlusEtaKK*(-3008. + 123.*LAL_PI*LAL_PI)))/96.;
  coeffs->k5 = k5 = 0.0;
  coeffs->k5l= k5l= 0.0;
  if ( SpinAlignedEOBversion == 2 )
  {
    coeffs->k5 = k5 = m1PlusEtaKK*m1PlusEtaKK
	       * (-4237./60.+128./5.*LAL_GAMMA+2275.*LAL_PI*LAL_PI/512.
	       - 1./3.*a*a*(k1p3-3.*k1*k2+3.*k3)
	       - (k1p3*k1p2-5.*k1p3*k2+5.*k1*k2*k2+5.*k1p2*k3-5.*k2*k3-5.*k1*k4)/5./m1PlusEtaKK/m1PlusEtaKK
	       + (k1p2*k1p2-4.*k1p2*k2+2.*k2*k2+4.*k1*k3-4.*k4)/2./m1PlusEtaKK+256./5.*log(2.));
    coeffs->k5l = k5l = m1PlusEtaKK*m1PlusEtaKK * 64./5.;
  }

  /*printf( "a = %.16e, k0 = %.16e, k1 = %.16e, k2 = %.16e, k3 = %.16e, k4 = %.16e, b3 = %.16e, bb3 = %.16e, KK = %.16e\n",
            a, coeffs->k0, coeffs->k1, coeffs->k2, coeffs->k3, coeffs->k4, coeffs->b3, coeffs->bb3, coeffs->KK );
  */

  /* Now calibrated parameters for spin models */
  coeffs->d1 = coeffs->d1v2 = 0.0;
  coeffs->dheffSS = coeffs->dheffSSv2 = 0.0;
  switch ( SpinAlignedEOBversion )
  {
     case 1:
       coeffs->d1 = -69.5;
       coeffs->dheffSS = 2.75;
       break;
     case 2:
       coeffs->d1v2 = -74.71 - 156.*eta + 627.5*eta*eta;
       coeffs->dheffSSv2 = 8.127 - 154.2*eta + 830.8*eta*eta;
       break;
     default:
       XLALPrintError( "XLAL Error - %s: wrong SpinAlignedEOBversion value, must be 1 or 2!\n", __func__ );
       XLAL_ERROR( XLAL_EINVAL );
       break;
  }

  return XLAL_SUCCESS;
}


/**
 * This function calculates the function \f$\Delta_t(r)\f$ which appears in the spinning EOB
 * potential function. Eqs. 7a and 8.
 */
static REAL8 XLALSimIMRSpinEOBHamiltonianDeltaT( 
        SpinEOBHCoeffs *coeffs, /**<< Pre-computed coefficients which appear in the function */
        const REAL8    r,       /**<< Current orbital radius (in units of total mass) */
        const REAL8    eta,     /**<< Symmetric mass ratio */
        const REAL8    a        /**<< Normalized deformed Kerr spin */
        )

{

  REAL8 a2;
  REAL8 u, u2, u3, u4, u5;
  REAL8 m1PlusetaKK;

  REAL8 bulk;
  REAL8 logTerms;
  REAL8 deltaU;
  REAL8 deltaT;

  u  = 1./r;
  u2 = u*u;
  u3 = u2*u;
  u4 = u2*u2;
  u5 = u4*u;

  a2 = a*a;

  m1PlusetaKK = -1. + eta * coeffs->KK;

  bulk = 1./(m1PlusetaKK*m1PlusetaKK) + (2.*u)/m1PlusetaKK + a2*u2;

  logTerms = 1. + eta*coeffs->k0 + eta*log(fabs(1. + coeffs->k1*u + coeffs->k2*u2 + coeffs->k3*u3 + coeffs->k4*u4
                                              + coeffs->k5*u5 + coeffs->k5l*u5*log(u)));
  /*printf(" a = %.16e, u = %.16e\n",a,u);
  printf( "k0 = %.16e, k1 = %.16e, k2 = %.16e, k3 = %.16e , k4 = %.16e, k5 = %.16e, k5l = %.16e\n",coeffs->k0,
         coeffs->k1,coeffs->k2,coeffs->k3,coeffs->k4,coeffs->k5,coeffs->k5l);
  printf( "bulk = %.16e, logTerms = %.16e\n", bulk, logTerms );*/
  deltaU = bulk*logTerms;
  deltaU = fabs(deltaU);
  deltaT = r*r*deltaU;

  return deltaT;
}


/**
 * This function calculates the function \f$\Delta_r(r)\f$ which appears in the spinning EOB
 * potential function. Eqs. 10a and 10b
 */
UNUSED static REAL8 XLALSimIMRSpinEOBHamiltonianDeltaR(
        SpinEOBHCoeffs *coeffs, /**<< Pre-computed coefficients which appear in the function */
        const REAL8    r,       /**<< Current orbital radius (in units of total mass) */
        const REAL8    eta,     /**<< Symmetric mass ratio */
        const REAL8    a        /**<< Normalized deformed Kerr spin */
        )
{

  
  REAL8 u2, u3;
  REAL8 D;
  REAL8 deltaT; /* The potential function, not a time interval... */
  REAL8 deltaR;

  u2 = 1./(r*r);
  u3 = u2 / r;

  D = 1. + log(1. + 6.*eta*u2 + 2.*(26. - 3.*eta)*eta*u3);

  deltaT = XLALSimIMRSpinEOBHamiltonianDeltaT( coeffs, r, eta, a );

  deltaR = deltaT*D;
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
XLALSimIMRSpinAlignedEOBCalcOmega(
                      const REAL8           values[],   /**<< Dynamical variables */
                      SpinEOBParams         *funcParams /**<< EOB parameters */
                      )
{
  static const REAL8 STEP_SIZE = 1.0e-4;

  HcapDerivParams params;

  /* Cartesian values for calculating the Hamiltonian */
    REAL8 cartValues[6] = {0.};

  gsl_function F;
  INT4         gslStatus;

  REAL8 omega;
  REAL8 r;

  /* The error in a derivative as measured by GSL */
  REAL8 absErr;

  /* Set up pointers for GSL */
  params.values  = cartValues;
  params.params  = funcParams;

  F.function = &GSLSpinAlignedHamiltonianWrapper;
  F.params   = &params;

  /* Populate the Cartesian values vector */
  /* We can assume phi is zero wlog */
  memset( cartValues, 0, sizeof( cartValues ) );
  cartValues[0] = r = values[0];
  cartValues[3] = values[2];
  cartValues[4] = values[3] / values[0];

  /* Now calculate omega. In the chosen co-ordinate system, */
  /* we need dH/dpy to calculate this, i.e. varyParam = 4   */
  params.varyParam = 4;
  XLAL_CALLGSL( gslStatus = gsl_deriv_central( &F, cartValues[4],
                  STEP_SIZE, &omega, &absErr ) );

  if ( gslStatus != GSL_SUCCESS )
  {
    XLALPrintError( "XLAL Error - %s: Failure in GSL function\n", __func__ );
    XLAL_ERROR_REAL8( XLAL_EFUNC );
  }
  
  omega = omega / r;

  return omega;
}

/**
 * Function to calculate the non-Keplerian coefficient for the spin-aligned EOB model.
 * radius \f$r\f$ times the cuberoot of the returned number is \f$r_\Omega\f$ defined in Eq. A2.
 * i.e. the function returns \f$(r_{\Omega} / r)^3\f$.
 */
UNUSED static REAL8
XLALSimIMRSpinAlignedEOBNonKeplerCoeff(
                      const REAL8           values[],   /**<< Dynamical variables */
                      SpinEOBParams         *funcParams /**<< EOB parameters */
                      )
{

  REAL8 omegaCirc;

  REAL8 tmpValues[4]= {0.};

  REAL8 r3;

  /* We need to find the values of omega assuming pr = 0 */
  memcpy( tmpValues, values, sizeof(tmpValues) );
  tmpValues[2] = 0.0;

  omegaCirc = XLALSimIMRSpinAlignedEOBCalcOmega( tmpValues, funcParams );
  if ( XLAL_IS_REAL8_FAIL_NAN( omegaCirc ) )
  {
    XLAL_ERROR_REAL8( XLAL_EFUNC );
  }

  r3 = values[0]*values[0]*values[0];

  return 1.0/(omegaCirc*omegaCirc*r3);
}
  
/* Wrapper for GSL to call the Hamiltonian function */
static double GSLSpinAlignedHamiltonianWrapper( double x, void *params )
{
  HcapDerivParams *dParams = (HcapDerivParams *)params;

  EOBParams *eobParams = dParams->params->eobParams;

  REAL8 tmpVec[6]= {0.};

  /* These are the vectors which will be used in the call to the Hamiltonian */
  REAL8Vector r, p;
  REAL8Vector *s1Vec = dParams->params->s1Vec;
  REAL8Vector *s2Vec = dParams->params->s2Vec;
  REAL8Vector *sigmaKerr = dParams->params->sigmaKerr;
  REAL8Vector *sigmaStar = dParams->params->sigmaStar;

  /* Use a temporary vector to avoid corrupting the main function */
  memcpy( tmpVec, dParams->values, 
               sizeof(tmpVec) );

  /* Set the relevant entry in the vector to the correct value */
  tmpVec[dParams->varyParam] = x;

  /* Set the LAL-style vectors to point to the appropriate things */
  r.length = p.length = 3;
  r.data     = tmpVec;
  p.data     = tmpVec+3;

  return XLALSimIMRSpinEOBHamiltonian( eobParams->eta, &r, &p, s1Vec, s2Vec, sigmaKerr, sigmaStar, dParams->params->tortoise, dParams->params->seobCoeffs ) / eobParams->eta;
}


/* ************************************************************************* */
/* ************************************************************************* */
/* ***************    PRECESSING EOB FUNCTIONS    ************************** */
/* ************************************************************************* */
/* ************************************************************************* */

/**
 * Functions to compute the inner product and cross products 
 * between vectors
 * */
UNUSED static REAL8 inner_product( const REAL8 values1[], 
                                   const REAL8 values2[] )
{
  REAL8 result = 0;
  for( int i = 0; i < 3 ; i++ )
    result += values1[i] * values2[i];

  return result;
}

UNUSED static REAL8* cross_product( const REAL8 values1[],
                                    const REAL8 values2[],
                                    REAL8 result[] )
{
  result[0] = values1[1]*values2[2] - values1[2]*values2[1];
  result[1] = values1[2]*values2[0] - values1[0]*values2[2];
  result[2] = values1[0]*values2[1] - values1[1]*values2[0];

  return result;
}


/**
 * Function to calculate the value of omega for the PRECESSING EOB waveform.
 * Needs the dynamics in Cartesian coordinates. 
 * 
 * First, the frame is rotated so that L is along the y-axis. 
 * this rotation includes the spins. 
 * 
 * Second, \f$\vec{r}\f$ and \f$\vec{p}\f$ are converted to polar coordinates
 * (and not the spins). As L is along the y-axis, \f$\theta\f$ defined as the
 * angle between L and the y-axis is 0, which is a cyclic coordinate now and
 * that fixes 
 * \f$p_\theta = 0\f$.
 * 
 * Third, \f$p_r\f$ is set to 0. 
 * 
 * Fourth, the polar \f$(r,\phi,p_r=0,p_\phi)\f$ and the Cartesian spin vectors
 * are used to compute the derivative
 * \f$\partial Hreal/\partial p_\phi |p_r=0\f$.
 */
static REAL8 XLALSimIMRSpinEOBCalcOmega(
                      const REAL8           values[],   /**<< Dynamical variables */
                      SpinEOBParams         *funcParams /**<< EOB parameters */
                      )
{
  int debugPK = 1;
  if (debugPK){
    for(int i =0; i < 12; i++)
      if( isnan(values[i]) ) {
        printf("XLALSimIMRSpinEOBCalcOmega::values %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f\n", values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9], values[10], values[11]);
          XLALPrintError( "XLAL Error - %s: nan in input values  \n", __func__);
          XLAL_ERROR( XLAL_EINVAL );
      }
  }
  
  /* ********************************************************************* */
  /* ************ Memory Allocation ************************************** */
  /* ********************************************************************* */
  static const REAL8 STEP_SIZE = 1.0e-4;
  REAL8 tmpvar = 0;

  HcapDerivParams params;

  /* Cartesian values for calculating the Hamiltonian */
  REAL8 cartValues[14] = {0.}, dvalues[14] = {0.};
  REAL8 cartvalues[14] = {0.}, polarvalues[6] = {0.}; /* The rotated cartesian/polar values */
  REAL8 polarRPcartSvalues[14] = {0.};
  memcpy( cartValues, values, 14 * sizeof(REAL8) );  
  
  UNUSED INT4 i, j;
  
  REAL8 rvec[3]  = {0.,0,0}, pvec[3]  = {0.,0,0}; 
  REAL8 s1vec[3] = {0.,0,0}, s2vec[3] = {0.,0,0}; 

  REAL8 rdotvec[3] = {0.,0,0}; 
  REAL8 rvecprime[3] = {0.,0,0}, pvecprime[3] = {0.,0,0}, 
        s1vecprime[3]= {0.,0,0}, s2vecprime[3]= {0.,0,0};
  REAL8 rvectmp[3]   = {0.,0,0}, pvectmp[3] = {0.,0,0},
        s1vectmp[3]  = {0.,0,0}, s2vectmp[3]= {0.,0,0};
  REAL8 LNhatprime[3]= {0.,0,0}, LNhatTmp[3]= {0.,0,0};
  REAL8 rcrossrdot[3] = {0.,0,0};

  REAL8 Rot1[3][3] ={{0.}}; // Rotation matrix for prevention of blowing up
  REAL8 Rot2[3][3] ={{0.}} ;
  REAL8 LNhat[3] = {0.,0,0};

  REAL8        Xhat[3] = {1, 0, 0};
  UNUSED REAL8 Yhat[3] = {0, 1, 0};
  UNUSED REAL8 Zhat[3] = {0, 0, 1};

  REAL8 Xprime[3] = {0.,0,0}, Yprime[3] = {0.,0,0}, Zprime[3] = {0.,0,0};
 
  gsl_function F;
  INT4         gslStatus;

  REAL8 omega;

  /* The error in a derivative as measured by GSL */
  REAL8 absErr;

  /* ********************************************************************* */
  /* ************ Main Logic begins ************************************ */
  /* ********************************************************************* */

  /* Copy over the coordinates and spins */
  memcpy( rvec,  values,   3*sizeof(REAL8) );
  memcpy( pvec,  values+3, 3*sizeof(REAL8) );
  memcpy( s1vec, values+6, 3*sizeof(REAL8) );
  memcpy( s2vec, values+9, 3*sizeof(REAL8) );

  /* Calculate rDot = \f$\partial Hreal / \partial p_r\f$ */
  memset( dvalues, 0, 14 * sizeof(REAL8) );
  if( XLALSpinHcapRvecDerivative( 0, values, dvalues, 
                                  (void*) funcParams) == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }
  memcpy( rdotvec, dvalues, 3*sizeof(REAL8) );

  if (debugPK){
    for(int ii =0; ii < 12; ii++)
      if( isnan(dvalues[ii]) ) {
        printf("XLALSimIMRSpinEOBCalcOmega::dvalues %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f\n", dvalues[0], dvalues[1], dvalues[2], dvalues[3], dvalues[4], dvalues[5], dvalues[6], dvalues[7], dvalues[8], dvalues[9], dvalues[10], dvalues[11]);
          XLALPrintError( "XLAL Error - %s: nan in dvalues \n", __func__);
          XLAL_ERROR( XLAL_EINVAL );
      }
  }
    
  /* Calculate LN = r cross rDot */
  cross_product( rvec, rdotvec, rcrossrdot );
  REAL8 rcrossrdotNorm = sqrt(inner_product( rcrossrdot, rcrossrdot ));
  for( i = 0; i < 3; i++ )
    rcrossrdot[i] /= rcrossrdotNorm;
  memcpy( LNhat, rcrossrdot, 3 * sizeof(REAL8) );


  /* ********************************************************************* */
  /* First, the frame is rotated so that L is along the y-axis. */
  /* this rotation includes the spins. */
  /* ********************************************************************* */

  // For Now , set first rotation matrix to identity
  // Check if LNhat and Xhat are too aligned, in which case rotate LNhat
  if( inner_product(LNhat, Xhat) < 0.9 )
  {
	  Rot1[0][0] = 1.; Rot1[0][1] = 0; Rot1[0][2] = 0;  
	  Rot1[1][0] = 0.; Rot1[1][1] = 1; Rot1[1][2] = 0;  
	  Rot1[2][0] = 0.; Rot1[2][1] = 0; Rot1[2][2] = 1;  
 
 	  memcpy(Xprime, LNhat, 3 * sizeof(REAL8));
	  cross_product( Xprime, Xhat, Yprime );
	  tmpvar = sqrt(inner_product(Yprime, Yprime));
    
	  for( i=0; i<3; i++)
      Yprime[i] /= tmpvar;
	  
    cross_product(Xprime, Yprime, Zprime);
	  tmpvar = sqrt(inner_product(Zprime, Zprime));
	  for( i=0; i<3; i++)
      Zprime[i] /= tmpvar;	  
  }
  else
  {
	  Rot1[0][0] = 1./sqrt(2); Rot1[0][1] = -1/sqrt(2); Rot1[0][2] = 0;  
	  Rot1[1][0] = 1./sqrt(2); Rot1[1][1] = 1./sqrt(2); Rot1[1][2] = 0;  
	  Rot1[2][0] = 0.;         Rot1[2][1] = 0;          Rot1[2][2] = 1;  
	  LNhatTmp[0] = LNhatTmp[1] = LNhatTmp[2] = 0.;
    
	  for(i=0; i<3; i++)
      for(j=0; j<3; j++)
        LNhatTmp[i] += Rot1[i][j]*LNhat[j];
    
	  memcpy(Xprime, LNhatTmp, 3*sizeof(REAL8));
	  cross_product(Xprime, Xhat, Yprime);
	  tmpvar = sqrt(inner_product(Yprime, Yprime));
    
	  for( i=0; i<3; i++)
      Yprime[i] /= tmpvar;
    
    cross_product(Xprime, Yprime, Zprime);
	  tmpvar = sqrt(inner_product(Zprime, Zprime));
	  for( i=0; i<3; i++)
      Zprime[i] /= tmpvar;  
  }
    
  Rot2[0][0] = Xprime[0]; Rot2[0][1] = Xprime[1]; Rot2[0][2] = Xprime[2];
  Rot2[1][0] = Yprime[0]; Rot2[1][1] = Yprime[1]; Rot2[1][2] = Yprime[2];
  Rot2[2][0] = Zprime[0]; Rot2[2][1] = Zprime[1]; Rot2[2][2] = Zprime[2];

  memset( rvectmp,    0, 3 * sizeof(REAL8) );
  memset( pvectmp,    0, 3 * sizeof(REAL8) );
  memset( s1vectmp,   0, 3 * sizeof(REAL8) );
  memset( s2vectmp,   0, 3 * sizeof(REAL8) );
  memset( rvecprime,  0, 3 * sizeof(REAL8) );
  memset( pvecprime,  0, 3 * sizeof(REAL8) );
  memset( s1vecprime, 0, 3 * sizeof(REAL8) );
  memset( s2vecprime, 0, 3 * sizeof(REAL8) );
  memset( LNhatprime, 0, 3 * sizeof(REAL8) );
  memset( LNhatTmp,   0, 3 * sizeof(REAL8) );

  /* Perform the actual rotation */
  for (i=0; i<3; i++)
    for(j=0; j<3; j++)
      {
        rvectmp[i]  += Rot1[i][j]*rvec[j];
        pvectmp[i]  += Rot1[i][j]*pvec[j];
        s1vectmp[i] += Rot1[i][j]*s1vec[j];
        s2vectmp[i] += Rot1[i][j]*s2vec[j];
        LNhatTmp[i] += Rot1[i][j]*LNhat[j];
      }
  for (i=0; i<3; i++)
    for(j=0; j<3; j++)
      {
        rvecprime[i]  += Rot2[i][j]*rvectmp[j];
        pvecprime[i]  += Rot2[i][j]*pvectmp[j];
        s1vecprime[i] += Rot2[i][j]*s1vectmp[j];
        s2vecprime[i] += Rot2[i][j]*s2vectmp[j];
        LNhatprime[i] += Rot2[i][j]*LNhatTmp[j];
      }

  memcpy(cartvalues,   rvecprime,  3*sizeof(REAL8));
  memcpy(cartvalues+3, pvecprime,  3*sizeof(REAL8));
  memcpy(cartvalues+6, s1vecprime, 3*sizeof(REAL8));
  memcpy(cartvalues+9, s2vecprime, 3*sizeof(REAL8));

  /* ********************************************************************* */
  /* Second, \f$\vec{r}\f$ and \f$\vec{p}\f$ are converted to polar 
   * coordinates (and not the spins). 
   * As L is along the y-axis, \f$\theta\f$ defined as the angle between 
   * L and the y-axis is 0, which is a cyclic coordinate now and that fixes 
   * \f$p_\theta = 0\f$. */
  /* ********************************************************************* */
  
  /** the polarvalues, respectively, are 
   * {r, \theta, \phi, p_r, p_\theta, p_\phi} */
  polarvalues[0] = sqrt(inner_product(rvecprime,rvecprime));
  polarvalues[1] = acos(rvecprime[0] / polarvalues[0]);
  polarvalues[2] = atan2(-rvecprime[1], rvecprime[2]);
  //polarvalues[3] = inner_product(rvecprime, pvecprime) / polarvalues[0];
  /* FIX p_r = 0 */
  polarvalues[3] = 0;
  
  REAL8 rvecprime_x_xhat[3] = {0.}, rvecprime_x_xhat_x_rvecprime[3] = {0.};
  cross_product(rvecprime, Xhat, rvecprime_x_xhat);
  cross_product(rvecprime_x_xhat, rvecprime, rvecprime_x_xhat_x_rvecprime);
  
  polarvalues[4] = -inner_product(rvecprime_x_xhat_x_rvecprime, pvecprime)
                              / polarvalues[0] / sin(polarvalues[1]);
  polarvalues[5] = -inner_product(rvecprime_x_xhat, pvecprime);

     
  /* ********************************************************************* */  /* Finally, Differentiate Hamiltonian w.r.t. p_\phi, keeping p_r = 0 */
  /* ********************************************************************* */
  
  /* Populate the vector specifying the dynamical variables in mixed frames */
  memcpy( polarRPcartSvalues, cartvalues, 12*sizeof(REAL8));
  memcpy( polarRPcartSvalues, polarvalues, 6*sizeof(REAL8));
  
  /* Set up pointers for GSL */
  params.values  = polarRPcartSvalues;
  params.params  = funcParams;

  F.function = &GSLSpinHamiltonianWrapperFordHdpphi;
  F.params   = &params;

  /* Now calculate omega. In the chosen co-ordinate system, */
  /* we need dH/dpphi to calculate this, i.e. varyParam = 5 */
  params.varyParam = 5;
  XLAL_CALLGSL( gslStatus = gsl_deriv_central( &F, polarvalues[5],
                  STEP_SIZE, &omega, &absErr ) );

  if ( gslStatus != GSL_SUCCESS )
  {
    XLALPrintError( "XLAL Error - %s: Failure in GSL function\n", __func__ );
    XLAL_ERROR_REAL8( XLAL_EFUNC );
  }
  
  return omega;
}


/**
 * Function to calculate the non-Keplerian coefficient for the PRECESSING EOB 
 * model.
 * 
 * radius \f$r\f$ times the cuberoot of the returned number is \f$r_\Omega\f$
 * defined in Eq. A2, i.e. the function returns 
 * \f$(r_{\Omega} / r)^3\f$ 
 *     = \f$1/(r^3 (\partial Hreal/\partial p_\phi |p_r=0)^2)\f$.
 */
static REAL8 
XLALSimIMRSpinEOBNonKeplerCoeff(
                      const REAL8           values[],   /**<< Dynamical variables */
                      SpinEOBParams         *funcParams /**<< EOB parameters */
                      )
{
  int debugPK = 1;
  if (debugPK){
    for(int i =0; i < 12; i++)
      if( isnan(values[i]) ) {
        printf("XLALSimIMRSpinEOBNonKeplerCoeff::values %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f\n",
        values[0], values[1], values[2], values[3], values[4], values[5],
        values[6], values[7], values[8], values[9], values[10], values[11]);
          XLALPrintError( "XLAL Error - %s: nan in values  \n", __func__);
          XLAL_ERROR( XLAL_EINVAL );
      }
  }

  REAL8 omegaCirc = 0;
  REAL8 tmpValues[14]= {0.}; 
  REAL8 r3;
      
  /* We need to find the values of omega assuming pr = 0 */
  memcpy( tmpValues, values, sizeof(tmpValues) );
  omegaCirc = XLALSimIMRSpinEOBCalcOmega( tmpValues, funcParams );
  
  if ( XLAL_IS_REAL8_FAIL_NAN( omegaCirc ) )
  {
    XLAL_ERROR_REAL8( XLAL_EFUNC );
  }

  r3 = pow(inner_product(values, values), 3./2.);
  return 1.0/(omegaCirc*omegaCirc*r3);
}

/**
 * Function to calculate numerical derivatives of the spin EOB Hamiltonian,
 * to get \f$dr/dt\f$, as decribed in Eqs. 21, 22, 26 and 27 of
 * Pan et al. PRD 81, 084041 (2010)
 * This function is not used by the spin-aligned SEOBNRv1 model.
 */
UNUSED static int XLALSpinHcapRvecDerivative(
                 double UNUSED     t,         /**<< UNUSED */
                 const  REAL8      values[],  /**<< Dynamical variables */
                 REAL8             dvalues[], /**<< Time derivatives of variables (returned) */
                 void             *funcParams /**<< EOB parameters */
                               )
{
  UNUSED int debugPK = 1;
  //if (debugPK){
    for(int i =0; i < 12; i++){
      if( isnan(values[i]) ) {
        printf("XLALSpinHcapRvecDerivative::values %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f\n", values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9], values[10], values[11]);
          XLALPrintError( "XLAL Error - %s: nan in input values \n", __func__);
          XLAL_ERROR( XLAL_EINVAL );
        }
        
      if( isnan(dvalues[i]) ) {
        printf("XLALSpinHcapRvecDerivative::dvalues %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f\n", dvalues[0], dvalues[1], dvalues[2], dvalues[3], dvalues[4], dvalues[5], dvalues[6], dvalues[7], dvalues[8], dvalues[9], dvalues[10], dvalues[11]);
          XLALPrintError( "XLAL Error - %s: nan in the input dvalues \n", __func__);
          XLAL_ERROR( XLAL_EINVAL );
        }
    }
  //}

  static const REAL8 STEP_SIZE = 1.0e-4;

  UNUSED static const INT4 lMax = 8;

  HcapDerivParams params;

  /* Since we take numerical derivatives wrt dynamical variables */
  /* but we want them wrt time, we use this temporary vector in  */
  /* the conversion */
  REAL8           tmpDValues[14] = {0.};

  REAL8           H; //Hamiltonian
  //REAL8           flux;

  gsl_function F;
  INT4         gslStatus;
  UINT4 SpinAlignedEOBversion;

  UINT4 i, j, k;//, l;
  
  REAL8Vector rVec, pVec;
  REAL8 rData[3] = {0.}, pData[3] = {0.};

  /* We need r, phi, pr, pPhi to calculate the flux */
  REAL8       UNUSED r;
  REAL8Vector UNUSED polarDynamics;
  REAL8       polData[4] = {0.};

  REAL8 mass1, mass2, eta;
  REAL8 UNUSED rrTerm2, pDotS1, pDotS2;
  REAL8Vector s1, s2, s1norm, s2norm, sKerr, sStar;
  REAL8       s1Data[3]= {0.}, s2Data[3]= {0.}, s1DataNorm[3]= {0.}, s2DataNorm[3]= {0.};
  REAL8       sKerrData[3]= {0.}, sStarData[3]= {0.};
  REAL8 /*magS1, magS2,*/ chiS, chiA, a, tplspin;
  REAL8	UNUSED s1dotL, s2dotL;
  REAL8	UNUSED	  rcrossrDot[3]= {0.}, rcrossrDotMag, s1dotLN, s2dotLN;


  /* Orbital angular momentum */
  REAL8 Lx, Ly, Lz, magL;
  REAL8 Lhatx, Lhaty, Lhatz;
  //REAL8 dLx, dLy, dLz;
  //REAL8 dLhatx, dLhaty, dMagL;

  //REAL8 alphadotcosi;

  //REAL8 rCrossV_x, rCrossV_y, rCrossV_z, omega;

  /* The error in a derivative as measured by GSL */
  REAL8 absErr;

    REAL8 tmpP[3]= {0.}, rMag, rMag2, prT;
	  REAL8 u, u2, u3, u4, u5, w2, a2;
	  REAL8 D, m1PlusetaKK, bulk, logTerms, deltaU, deltaT, deltaR;
	  REAL8 UNUSED eobD_r, deltaU_u, deltaU_r, deltaT_r;
	  REAL8 dcsi, csi;
	  
	REAL8 tmpValues[12]= {0.};
    REAL8 UNUSED Tmatrix[3][3]= {{0.}}, invTmatrix[3][3]= {{0.}}, dTijdXk[3][3][3]= {{{0.}}};
	//REAL8 tmpPdotT1[3], tmpPdotT2[3], tmpPdotT3[3]; // 3 terms of Eq. A5
	
  /* Set up pointers for GSL */ 
  params.values  = values;
  params.params  = (SpinEOBParams *)funcParams;

  F.function = &GSLSpinHamiltonianWrapperForRvecDerivs;
  F.params   = &params;

  mass1 = params.params->eobParams->m1;
  mass2 = params.params->eobParams->m2;
  eta   = params.params->eobParams->eta;
  SpinAlignedEOBversion = params.params->seobCoeffs->SpinAlignedEOBversion;  
  SpinEOBHCoeffs *coeffs = (SpinEOBHCoeffs*) params.params->seobCoeffs;
  
  /* For precessing binaries, the effective spin of the Kerr 
   * background evolves with time. The coefficients used to compute
   * the Hamiltonian depend on the Kerr spin, and hence need to 
   * be updated for the current spin values */
  if ( 0 )
  {/*{{{*/
    /* Set up structures and calculate necessary (spin-only) PN parameters */
    /* Due to precession, these need to get calculated in every step */
    //memset( params.params->seobCoeffs, 0, sizeof(SpinEOBHCoeffs) );
    
    REAL8 tmps1Data[3]= {0.}, tmps2Data[3]= {0.}; REAL8Vector tmps1Vec, tmps2Vec;
    memcpy( tmps1Data, values+6, 3*sizeof(REAL8) );
    memcpy( tmps2Data, values+9, 3*sizeof(REAL8) );
    tmps1Vec.data   = tmps1Data; tmps2Vec.data   = tmps2Data;
    tmps1Vec.length = tmps2Vec.length = 3;
    
    REAL8Vector *tmpsigmaKerr = NULL;
    REAL8Vector *tmpsigmaStar = NULL;
    if ( !(tmpsigmaKerr = XLALCreateREAL8Vector( 3 )) )
    {
      XLAL_ERROR( XLAL_ENOMEM );
    }
    
    if ( !(tmpsigmaStar = XLALCreateREAL8Vector( 3 )) )
    {
      XLAL_ERROR( XLAL_ENOMEM );
    }

    if ( XLALSimIMRSpinEOBCalculateSigmaKerr( tmpsigmaKerr, mass1, mass2, 
                              &tmps1Vec, &tmps2Vec ) == XLAL_FAILURE )
    {
      XLALDestroyREAL8Vector( tmpsigmaKerr );
      XLAL_ERROR( XLAL_EFUNC );
    }
    
    if ( XLALSimIMRSpinEOBCalculateSigmaStar( tmpsigmaStar, mass1, mass2, 
                              &tmps1Vec, &tmps2Vec ) == XLAL_FAILURE )
    {
      XLALDestroyREAL8Vector( tmpsigmaKerr );
      XLALDestroyREAL8Vector( tmpsigmaStar );
      XLAL_ERROR( XLAL_EFUNC );
    }
    
    /* Update a with the Kerr background spin
     * Pre-compute the Hamiltonian coefficients            */
    //REAL8Vector *delsigmaKerr 	= params.params->sigmaKerr;
    params.params->sigmaKerr 	= tmpsigmaKerr;
    params.params->sigmaStar 	= tmpsigmaStar;
    params.params->a 		= sqrt( tmpsigmaKerr->data[0]*tmpsigmaKerr->data[0]
				+ tmpsigmaKerr->data[1]*tmpsigmaKerr->data[1]
				+ tmpsigmaKerr->data[2]*tmpsigmaKerr->data[2] );
    //tmpsigmaKerr->data[2];
    if ( XLALSimIMRCalculateSpinEOBHCoeffs( params.params->seobCoeffs, eta,
			params.params->a, SpinAlignedEOBversion ) == XLAL_FAILURE )
    {
      XLALDestroyREAL8Vector( params.params->sigmaKerr );
      XLAL_ERROR( XLAL_EFUNC );
    }
   
    params.params->seobCoeffs->SpinAlignedEOBversion = SpinAlignedEOBversion;
    /* Release the old memory */
    //if(0)XLALDestroyREAL8Vector( delsigmaKerr );
  /*}}}*/}

  /* Set the position/momenta vectors to point to the appropriate things */
  rVec.length = pVec.length = 3;
  rVec.data   = rData;
  pVec.data   = pData;
  memcpy( rData, values, sizeof(rData) );
  memcpy( pData, values+3, sizeof(pData) );

  /* We need to re-calculate the parameters at each step as precessing
   * spins will not be constant */
  /* TODO: Modify so that only spin terms get re-calculated */

  /* We cannot point to the values vector directly as it leads to a warning */
  s1.length = s2.length = s1norm.length = s2norm.length = 3;
  s1.data = s1Data;
  s2.data = s2Data;
  s1norm.data = s1DataNorm;
  s2norm.data = s2DataNorm;

  memcpy( s1Data, values+6, 3*sizeof(REAL8) );
  memcpy( s2Data, values+9, 3*sizeof(REAL8) );
  memcpy( s1DataNorm, values+6, 3*sizeof(REAL8) );
  memcpy( s2DataNorm, values+9, 3*sizeof(REAL8) );

  for ( i = 0; i < 3; i++ )
  {
	  s1Data[i] *= (mass1+mass2)*(mass1+mass2);
	  s2Data[i] *= (mass1+mass2)*(mass1+mass2);
  }

  sKerr.length = 3;
  sKerr.data   = sKerrData; 
  XLALSimIMRSpinEOBCalculateSigmaKerr( &sKerr, mass1, mass2, &s1, &s2 );

  sStar.length = 3;
  sStar.data   = sStarData;
  XLALSimIMRSpinEOBCalculateSigmaStar( &sStar, mass1, mass2, &s1, &s2 );

  a = sqrt(sKerr.data[0]*sKerr.data[0] + sKerr.data[1]*sKerr.data[1] 
      + sKerr.data[2]*sKerr.data[2]);

  if (isnan(a)){
      XLALPrintError( "XLAL Error - %s: a = nan   \n", __func__);
      XLAL_ERROR( XLAL_EINVAL );
  } 
  if(debugPK && isnan(a))
    printf("a is nan in XLALSpinHcapRvecDerivative \n");
    
  ///* set the tortoise flag to 2 */
  //INT4 oldTortoise = params.params->tortoise;
  //params.params->tortoise = 2;
	  
  /* Convert momenta to p */
  rMag = sqrt(rData[0]*rData[0] + rData[1]*rData[1] + rData[2]*rData[2]);
  prT = pData[0]*(rData[0]/rMag) + pData[1]*(rData[1]/rMag) 
					+ pData[2]*(rData[2]/rMag);
	
	  rMag2 = rMag * rMag;
	  u  = 1./rMag;
      u2 = u*u;
      u3 = u2*u;
      u4 = u2*u2;
      u5 = u4*u;
      a2 = a*a;
      w2 = rMag2 + a2;
      /* Eq. 5.83 of BB1, inverse */
      D = 1. + log(1. + 6.*eta*u2 + 2.*(26. - 3.*eta)*eta*u3);
      eobD_r =  (u2/(D*D))*(12.*eta*u + 6.*(26. - 3.*eta)*eta*u2)/(1. 
			+ 6.*eta*u2 + 2.*(26. - 3.*eta)*eta*u3);
      m1PlusetaKK = -1. + eta * coeffs->KK;
	  /* Eq. 5.75 of BB1 */
      bulk = 1./(m1PlusetaKK*m1PlusetaKK) + (2.*u)/m1PlusetaKK + a2*u2;
	  /* Eq. 5.73 of BB1 */
	  logTerms = 1. + eta*coeffs->k0 + eta*log(fabs(1. + coeffs->k1*u
		+ coeffs->k2*u2 + coeffs->k3*u3 + coeffs->k4*u4
		+ coeffs->k5*u5 + coeffs->k5l*u5*log(u)));
	  /* Eq. 5.73 of BB1 */
      deltaU = bulk*logTerms;
    deltaU = fabs(deltaU);

      /* Eq. 5.71 of BB1 */
      deltaT = rMag2*deltaU;
      /* ddeltaU/du */
	  deltaU_u = 2.*(1./m1PlusetaKK + a2*u)*logTerms + 
		bulk * (eta*(coeffs->k1 + u*(2.*coeffs->k2 + u*(3.*coeffs->k3 
		+ u*(4.*coeffs->k4 + 5.*(coeffs->k5+coeffs->k5l*log(u))*u)))))
		/ (1. + coeffs->k1*u + coeffs->k2*u2 + coeffs->k3*u3 
		+ coeffs->k4*u4 + (coeffs->k5+coeffs->k5l*log(u))*u5);
	  deltaU_r = -u2 * deltaU_u;
      /* Eq. 5.38 of BB1 */
      deltaR = deltaT*D;	   
	  if ( params.params->tortoise )
		  csi = sqrt( fabs(deltaT * deltaR) )/ w2; /* Eq. 28 of Pan et al. PRD 81, 084041 (2010) */
	  else
	      csi = 1.0;

	  for( i = 0; i < 3; i++ )
	  { 
		  tmpP[i] = pData[i] - (rData[i]/rMag) * prT * (csi-1.)/csi;
	  }
	  
		
  /* Calculate the T-matrix, required to convert P from tortoise to 
   * non-tortoise coordinates, and/or vice-versa. This is given explicitly
   * in Eq. A3 of 0912.3466 */
  for( i = 0; i < 3; i++ )
	for( j = 0; j <= i; j++ )
	{
		Tmatrix[i][j] = Tmatrix[j][i] = (rData[i]*rData[j]/rMag2) 
									* (csi - 1.);
		
		invTmatrix[i][j] = invTmatrix[j][i] = 
				- (csi - 1)/csi * (rData[i]*rData[j]/rMag2);
		
		if( i==j ){ 
			Tmatrix[i][j]++;  
			invTmatrix[i][j]++;  }
		
	}

  dcsi = csi * (2./rMag + deltaU_r/deltaU) + csi*csi*csi 
      / (2.*rMag2*rMag2 * deltaU*deltaU) * ( rMag*(-4.*w2)/D - eobD_r*(w2*w2));
   
  for( i = 0; i < 3; i++ )
	for( j = 0; j < 3; j++ )
		for( k = 0; k < 3; k++ )
		{
			dTijdXk[i][j][k]  = 
		(rData[i]*XLALKronecker(j,k) + XLALKronecker(i,k)*rData[j]) 
		*(csi - 1.)/rMag2 
		+ rData[i]*rData[j]*rData[k]/rMag2/rMag*(-2./rMag*(csi - 1.) + dcsi);
		}
	
  //if (debugPK){
    for(i =0; i < 12; i++){
      if( isnan(values[i]) ) {
        printf("XLALSpinHcapRvecDerivative (just before diff)::values %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f\n", values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9], values[10], values[11]);
          XLALPrintError( "XLAL Error - %s: values = nan   \n", __func__);
          XLAL_ERROR( XLAL_EINVAL );
        }
        
      if( isnan(dvalues[i]) ) {
        printf("XLALSpinHcapRvecDerivative (just before diff)::dvalues %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f\n", dvalues[0], dvalues[1], dvalues[2], dvalues[3], dvalues[4], dvalues[5], dvalues[6], dvalues[7], dvalues[8], dvalues[9], dvalues[10], dvalues[11]);
          XLALPrintError( "XLAL Error - %s: dvalues = nan   \n", __func__);
          XLAL_ERROR( XLAL_EINVAL );
        }
      }
//}

  /* Now calculate derivatives w.r.t. each parameter */
  for ( i = 0; i < 6; i++ )
  {
    params.varyParam = i;
    if ( i >=6 && i < 9 )
    {
      params.params->seobCoeffs->updateHCoeffs = 1;
      XLAL_CALLGSL( gslStatus = gsl_deriv_central( &F, values[i],
                      STEP_SIZE*mass1*mass1, &tmpDValues[i], &absErr ) );
    }
    else if ( i >= 9 )
    {
      params.params->seobCoeffs->updateHCoeffs = 1;
      XLAL_CALLGSL( gslStatus = gsl_deriv_central( &F, values[i],
                      STEP_SIZE*mass2*mass2, &tmpDValues[i], &absErr ) );
    }
    else if ( i < 3 )
    {
	    params.params->seobCoeffs->updateHCoeffs = 1;
      	params.params->tortoise = 2;
		memcpy( tmpValues, params.values, sizeof(tmpValues) );
		tmpValues[3] = tmpP[0]; tmpValues[4] = tmpP[1]; tmpValues[5] = tmpP[2];
		params.values = tmpValues;
		
		XLAL_CALLGSL( gslStatus = gsl_deriv_central( &F, values[i], 
                      STEP_SIZE, &tmpDValues[i], &absErr ) );
                      
        params.values = values;
        params.params->tortoise = 1;
	}
    else
    {
      params.params->seobCoeffs->updateHCoeffs = 1;
      XLAL_CALLGSL( gslStatus = gsl_deriv_central( &F, values[i], 
                      STEP_SIZE, &tmpDValues[i], &absErr ) );
    }
    if ( gslStatus != GSL_SUCCESS )
    {
      XLALPrintError( "XLAL Error %s - Failure in GSL function\n", __func__ );
      XLAL_ERROR( XLAL_EFUNC );
    }
  }
  if (debugPK){
    for( i =0; i < 12; i++)
      if( isnan(tmpDValues[i]) ) {
        printf("XLALSpinHcapRvecDerivative (just after diff)::tmpDValues %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f\n", tmpDValues[0], tmpDValues[1], tmpDValues[2], tmpDValues[3], tmpDValues[4], tmpDValues[5], tmpDValues[6], tmpDValues[7], tmpDValues[8], tmpDValues[9], tmpDValues[10], tmpDValues[11]);
        }
    }

  /* Calculate the orbital angular momentum */
  Lx = values[1]*values[5] - values[2]*values[4];
  Ly = values[2]*values[3] - values[0]*values[5];
  Lz = values[0]*values[4] - values[1]*values[3];

  magL = sqrt( Lx*Lx + Ly*Ly + Lz*Lz );

  Lhatx = Lx/magL;
  Lhaty = Ly/magL;
  Lhatz = Lz/magL;

  /* Calculate the polar data */
  polarDynamics.length = 4;
  polarDynamics.data   = polData;

  r = polData[0] = sqrt( values[0]*values[0] + values[1]*values[1] 
						+ values[2]*values[2] );
  polData[1] = 0;
  polData[2] = (values[0]*values[3] + values[1]*values[4] 
				+ values[2]*values[5]) / polData[0];
  polData[3] = magL;


  /*Compute \vec{S_i} \dot \vec{L}	*/
  s1dotL = (s1Data[0]*Lhatx + s1Data[1]*Lhaty + s1Data[2]*Lhatz)
			/ (mass1*mass1);
  s2dotL = (s2Data[0]*Lhatx + s2Data[1]*Lhaty + s2Data[2]*Lhatz)
			/ (mass2*mass2);

  /*Compute \vec{L_N} = \vec{r} \times \.{\vec{r}}, 
   * \vec{S_i} \dot \vec{L_N} and chiS and chiA		*/
  rcrossrDot[0] = values[1]*tmpDValues[5] - values[2]*tmpDValues[4];
  rcrossrDot[1] = values[2]*tmpDValues[3] - values[0]*tmpDValues[5];
  rcrossrDot[2] = values[0]*tmpDValues[4] - values[1]*tmpDValues[3];
  rcrossrDotMag = sqrt( rcrossrDot[0]*rcrossrDot[0] 
		+ rcrossrDot[1]*rcrossrDot[1]	+ rcrossrDot[2]*rcrossrDot[2] );

  rcrossrDot[0] /= rcrossrDotMag;
  rcrossrDot[1] /= rcrossrDotMag;
  rcrossrDot[2] /= rcrossrDotMag;
  
  s1dotLN = (s1Data[0]*rcrossrDot[0] + s1Data[1]*rcrossrDot[1] 
		        + s1Data[2]*rcrossrDot[2]) / (mass1*mass1);
  s2dotLN = (s2Data[0]*rcrossrDot[0] + s2Data[1]*rcrossrDot[1] 
                	+ s2Data[2]*rcrossrDot[2]) / (mass2*mass2);
 
  chiS = 0.5 * (s1dotLN + s2dotLN);
  chiA = 0.5 * (s1dotLN - s2dotLN);

  /* Compute the test-particle limit spin of the deformed-Kerr background */
  /* TODO: Check this is actually the way it works in LAL */
  switch ( SpinAlignedEOBversion )
  {
     case 1:
       tplspin = 0.0;
       break;
     case 2:
       tplspin = (1.-2.*eta) * chiS + (mass1 - mass2)/(mass1 + mass2) * chiA;
       break;
     default:
       XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
       XLAL_ERROR( XLAL_EINVAL );
       break;
  }

  for( i = 0; i< 3; i++ )
  {
    params.params->s1Vec->data[i]     = s1norm.data[i];
    params.params->s2Vec->data[i]     = s2norm.data[i];
    params.params->sigmaStar->data[i] = sStar.data[i];
    params.params->sigmaKerr->data[i] = sKerr.data[i];
  }
 
  //params.params->s1Vec     = &s1norm;
  //params.params->s2Vec     = &s2norm;
  //params.params->sigmaStar = &sStar;
  //params.params->sigmaKerr = &sKerr;
  params.params->a         = a;
 
  XLALSimIMREOBCalcSpinFacWaveformCoefficients( 
		params.params->eobParams->hCoeffs, mass1, mass2, eta, tplspin, 
		chiS, chiA, SpinAlignedEOBversion );
  XLALSimIMRCalculateSpinEOBHCoeffs( params.params->seobCoeffs, eta, a, 
      SpinAlignedEOBversion );

  H = XLALSimIMRSpinEOBHamiltonian( eta, &rVec, &pVec, &s1norm, &s2norm, 
	&sKerr, &sStar, params.params->tortoise, params.params->seobCoeffs ); 
  H = H * (mass1 + mass2);
  
  /* Now make the conversion */
  /* rVectorDot */
  for( i = 0; i < 3; i++ )
	  for( j = 0, dvalues[i] = 0.; j < 3; j++ )
		  dvalues[i] += tmpDValues[j+3]*Tmatrix[i][j];
    
    //if (debugPK) printf("XLALSpinHcapRvecDerivative::tmpDValues = %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f\n",tmpDValues[0],tmpDValues[1],tmpDValues[2],tmpDValues[3],tmpDValues[4],tmpDValues[5]);


  return XLAL_SUCCESS;
}

/**
 * Wrapper for GSL to call the Hamiltonian function. This is simply the function
 * GSLSpinHamiltonianWrapper copied over. The alternative was to make it non-static
 * which increases runtime as static functions can be better optimized.
 */
static double GSLSpinHamiltonianWrapperForRvecDerivs( double x, void *params )
{
  int debugPK = 1;
  HcapDerivParams *dParams = (HcapDerivParams *)params;

  EOBParams *eobParams = (EOBParams*) dParams->params->eobParams;
  SpinEOBHCoeffs UNUSED *coeffs = (SpinEOBHCoeffs*) dParams->params->seobCoeffs;
  
  REAL8 tmpVec[12]= {0.};
  REAL8 s1normData[3]= {0.}, s2normData[3]= {0.}, sKerrData[3]= {0.}, sStarData[3]= {0.};

  /* These are the vectors which will be used in the call to the Hamiltonian */
  REAL8Vector r, p, spin1, spin2, spin1norm, spin2norm;
  REAL8Vector sigmaKerr, sigmaStar;

  INT4 i;
  REAL8 a;
  REAL8 m1 = eobParams->m1;
  REAL8 m2 = eobParams->m2;
  REAL8 UNUSED mT2 = (m1+m2)*(m1+m2);
  REAL8 UNUSED eta = m1*m2/mT2;

  INT4 oldTortoise = dParams->params->tortoise;
  /* Use a temporary vector to avoid corrupting the main function */
  memcpy( tmpVec, dParams->values, sizeof(tmpVec) );

  if (debugPK){
    for( i =0; i < 12; i++)
      if( isnan(tmpVec[i]) ) {
        printf("GSLSpinHamiltonianWrapperForRvecDerivs (from input)::tmpVec %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f\n", tmpVec[0], tmpVec[1], tmpVec[2], tmpVec[3], tmpVec[4], tmpVec[5], tmpVec[6], tmpVec[7], tmpVec[8], tmpVec[9], tmpVec[10], tmpVec[11]);
        }
    }

  /* Set the relevant entry in the vector to the correct value */
  tmpVec[dParams->varyParam] = x;

  /* Set the LAL-style vectors to point to the appropriate things */
  r.length = p.length = spin1.length = spin2.length = spin1norm.length = spin2norm.length = 3;
  sigmaKerr.length = sigmaStar.length = 3;
  r.data     = tmpVec;
  p.data     = tmpVec+3;
  
  spin1.data = tmpVec+6;
  spin2.data = tmpVec+9;
  spin1norm.data = s1normData;
  spin2norm.data = s2normData;
  sigmaKerr.data = sKerrData;
  sigmaStar.data = sStarData;

  memcpy( s1normData, tmpVec+6, 3*sizeof(REAL8) );
  memcpy( s2normData, tmpVec+9, 3*sizeof(REAL8) );

  /* To compute the SigmaKerr and SigmaStar, we need the non-normalized
   * spin values, i.e. S_i. The spins being evolved are S_i/M^2. */
  for ( i = 0; i < 3; i++ )
  {
	 spin1.data[i]  *= mT2;
	 spin2.data[i]  *= mT2;
	 
     //s1normData[i] /= mT2;
     //s2normData[i] /= mT2;
  }

  /* Calculate various spin parameters */
  XLALSimIMRSpinEOBCalculateSigmaKerr( &sigmaKerr, eobParams->m1, 
				eobParams->m2, &spin1, &spin2 );
  XLALSimIMRSpinEOBCalculateSigmaStar( &sigmaStar, eobParams->m1, 
				eobParams->m2, &spin1, &spin2 );
  a = sqrt( sigmaKerr.data[0]*sigmaKerr.data[0] 
			+ sigmaKerr.data[1]*sigmaKerr.data[1]
            + sigmaKerr.data[2]*sigmaKerr.data[2] );
  //printf( "a = %e\n", a );
  //printf( "aStar = %e\n", sqrt( sigmaStar.data[0]*sigmaStar.data[0] + sigmaStar.data[1]*sigmaStar.data[1] + sigmaStar.data[2]*sigmaStar.data[2]) );
  if ( isnan( a ) )
  {
      printf( "a is nan in GSLSpinHamiltonianWrapperForRvecDerivs!!\n");
      XLALPrintError( "XLAL Error - %s: a = nan   \n", __func__);
      XLAL_ERROR( XLAL_EINVAL );
  }
  //XLALSimIMRCalculateSpinEOBHCoeffs( dParams->params->seobCoeffs, eobParams->eta, a );
  /* If computing the derivative w.r.t. the position vector, we need to 
   * pass p and not p* to the Hamiltonian. This is so because we want to 
   * hold p constant as we compute dH/dx. The way to do this is set the
   * tortoise flag = 2, and convert the momenta being evolved, i.e. p*, 
   * to p and pass that as input */
   // TODO
#if 0
  if ( dParams->varyParam < 3 )
  {
	  /* set the tortoise flag to 2 */
	  dParams->params->tortoise = 2;
	  
	  /* Convert momenta to p */
	  REAL8 tmpP[3]= {0.};
	  REAL8 rMag = sqrt(r.data[0]*r.data[0] + r.data[1]*r.data[1] 
					+ r.data[2]*r.data[2]);
	  REAL8 prT = p.data[0]*(r.data[0]/rMag) + p.data[1]*(r.data[1]/rMag) 
					+ p.data[2]*(r.data[2]/rMag);
					
	  REAL8 u, u2, u3, u4, u5, w2, a2;
	  REAL8 csi; 
	  
	  u  = 1./rMag;
      u2 = u*u;
      u3 = u2*u;
      u4 = u2*u2;
      u5 = u4*u;
      a2 = a*a;
      w2 = rMag*rMag + a2;
      /* Eq. 5.83 of BB1, inverse */
      REAL8 D = 1. + log(1. + 6.*eta*u2 + 2.*(26. - 3.*eta)*eta*u3);
      REAL8 m1PlusetaKK = -1. + eta * coeffs->KK;
	  /* Eq. 5.75 of BB1 */
      REAL8 bulk = 1./(m1PlusetaKK*m1PlusetaKK) + (2.*u)/m1PlusetaKK + a2*u2;
	  /* Eq. 5.73 of BB1 */
	  REAL8 logTerms = 1. + eta*coeffs->k0 + eta*log(fabs(1. + coeffs->k1*u
		+ coeffs->k2*u2 + coeffs->k3*u3 + coeffs->k4*u4
		+ coeffs->k5*u5 + coeffs->k5l*u5*log(u)));
	  /* Eq. 5.73 of BB1 */
      REAL8 deltaU = bulk*logTerms;
      deltaU=fabs(deltaU);
      /* Eq. 5.71 of BB1 */
      REAL8 deltaT = rMag*rMag*deltaU;
      /* Eq. 5.38 of BB1 */
      REAL8 deltaR = deltaT*D;	   
	  if ( oldTortoise )
		  csi = sqrt( fabs(deltaT * deltaR) )/ w2; /* Eq. 28 of Pan et al. PRD 81, 084041 (2010) */
	  else
	      csi = 1.0;

	  for( i = 0; i < 3; i++ )
	  { 
		  tmpP[i] = p.data[i] - (r.data[i]/rMag) * prT * (csi-1.)/csi;
	  }
	  memcpy( p.data, tmpP, sizeof(tmpP) );
  }
#endif
//  printf( "Hamiltonian = %e\n", XLALSimIMRSpinEOBHamiltonian( eobParams->eta, &r, &p, &sigmaKerr, &sigmaStar, dParams->params->seobCoeffs ) );
    double magR = r.data[0]*r.data[0] + r.data[1]*r.data[1] + r.data[2]*r.data[2];

  if(debugPK) {
    if(0 && magR < 1.96 * 1.96) {
      printf("GSLSpinHamiltonianWrapperForRvecDerivs (JUST inputs)::tmpVec %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f\n", tmpVec[0], tmpVec[1], tmpVec[2], tmpVec[3], tmpVec[4], tmpVec[5], tmpVec[6], tmpVec[7], tmpVec[8], tmpVec[9], tmpVec[10], tmpVec[11]);
      
      printf(" R = %3.10f\n\n", sqrt(magR));
    }
  }
  
  REAL8 SpinEOBH = XLALSimIMRSpinEOBHamiltonian( eobParams->eta, &r, &p, &spin1norm, &spin2norm, &sigmaKerr, &sigmaStar, dParams->params->tortoise, dParams->params->seobCoeffs ) / eobParams->eta;
  
  if( isnan(SpinEOBH) )
    {
      printf("H is nan in GSLSpinHamiltonianWrapperForRvecDerivs. \n");
    
      printf("GSLSpinHamiltonianWrapperForRvecDerivs (JUST inputs)::tmpVec %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f\n", tmpVec[0], tmpVec[1], tmpVec[2], tmpVec[3], tmpVec[4], tmpVec[5], tmpVec[6], tmpVec[7], tmpVec[8], tmpVec[9], tmpVec[10], tmpVec[11]);
      
      printf(" R = %3.10f\n\n", sqrt(magR));
      XLALPrintError( "XLAL Error - %s: H = nan   \n", __func__);
      XLAL_ERROR( XLAL_EINVAL );
    }
    
  if ( dParams->varyParam < 3 )dParams->params->tortoise = oldTortoise;
  return SpinEOBH;
}


/**
 * Wrapper for GSL to call the Hamiltonian function. This is simply the function
 * GSLSpinHamiltonianWrapper copied over. The alternative was to make it non-static
 * which increases runtime as static functions can be better optimized.
 */
static double GSLSpinHamiltonianWrapperFordHdpphi( double x, void *params )
{
  HcapDerivParams *dParams = (HcapDerivParams *)params;

  EOBParams *eobParams = (EOBParams*) dParams->params->eobParams;
  SpinEOBHCoeffs UNUSED *coeffs = (SpinEOBHCoeffs*) dParams->params->seobCoeffs;
  
    REAL8 tmpVec[12] = {0.};
  REAL8 rpolar[3] = {0.}, rcart[3] = {0.}, ppolar[3] = {0.}, pcart[3] = {0.};
  REAL8 s1normData[3] = {0.}, s2normData[3] = {0.}, sKerrData[3] = {0.}, sStarData[3] = {0.};

  /* These are the vectors which will be used in the call to the Hamiltonian */
  REAL8Vector r, p, spin1, spin2, spin1norm, spin2norm;
  REAL8Vector sigmaKerr, sigmaStar;

  INT4 i;
  REAL8 a;
  REAL8 m1 = eobParams->m1;
  REAL8 m2 = eobParams->m2;
  REAL8 UNUSED mT2 = (m1+m2)*(m1+m2);
  REAL8 UNUSED eta = m1*m2/mT2;

  /* Use a temporary vector to avoid corrupting the main function */
  memcpy( tmpVec, dParams->values, sizeof(tmpVec) );

  /* Set the relevant entry in the vector to the correct value */
  tmpVec[dParams->varyParam] = x;

  /* Set the LAL-style vectors to point to the appropriate things */
  r.length = p.length = spin1.length = spin2.length = spin1norm.length = spin2norm.length = 3;
  sigmaKerr.length = sigmaStar.length = 3;
  
  /* Now rotate the R and P vectors from polar coordinates to Cartesian */
  memcpy( rpolar, tmpVec, 3*sizeof(REAL8)); 
  memcpy( ppolar, tmpVec+3, 3*sizeof(REAL8));
  
  rcart[0] = rpolar[0] * cos(rpolar[1]);
  rcart[1] =-rpolar[0] * sin(rpolar[1])*sin(rpolar[2]);
  rcart[2] = rpolar[0] * sin(rpolar[1])*cos(rpolar[2]);

  if( rpolar[1]==0. || rpolar[1]==LAL_PI )
  {
    rpolar[1] = LAL_PI/2.;
    
    if( rpolar[1]==0.)
      rpolar[2] = 0.;
    else
      rpolar[2] = LAL_PI;

    pcart[0] = ppolar[0]*sin(rpolar[1])*cos(rpolar[2]) 
					+ ppolar[1]/rpolar[0]*cos(rpolar[1])*cos(rpolar[2]) 
					- ppolar[2]/rpolar[0]/sin(rpolar[1])*sin(rpolar[2]);
    pcart[1] = ppolar[0]*sin(rpolar[1])*sin(rpolar[2]) 
					+ ppolar[1]/rpolar[0]*cos(rpolar[1])*sin(rpolar[2]) 
					+ ppolar[2]/rpolar[0]/sin(rpolar[1])*cos(rpolar[2]);
    pcart[2] = ppolar[0]*cos(rpolar[1])- ppolar[1]/rpolar[0]*sin(rpolar[1]);		
  }
  else
  {
    pcart[0] = ppolar[0]*cos(rpolar[1]) -ppolar[1]/rpolar[0]*sin(rpolar[1]);
    pcart[1] =-ppolar[0]*sin(rpolar[1])*sin(rpolar[2]) 
					-ppolar[1]/rpolar[0]*cos(rpolar[1])*sin(rpolar[2]) 
					-ppolar[2]/rpolar[0]/sin(rpolar[1])*cos(rpolar[2]);
    pcart[2] = ppolar[0]*sin(rpolar[1])*cos(rpolar[2]) 
					+ppolar[1]/rpolar[0]*cos(rpolar[1])*cos(rpolar[2]) 
					-ppolar[2]/rpolar[0]/sin(rpolar[1])*sin(rpolar[2]);
  }
 
  r.data     = rcart;
  p.data     = pcart;
  
  spin1.data = tmpVec+6;
  spin2.data = tmpVec+9;
  spin1norm.data = s1normData;
  spin2norm.data = s2normData;
  sigmaKerr.data = sKerrData;
  sigmaStar.data = sStarData;

  memcpy( s1normData, tmpVec+6, 3*sizeof(REAL8) );
  memcpy( s2normData, tmpVec+9, 3*sizeof(REAL8) );

  /* To compute the SigmaKerr and SigmaStar, we need the non-normalized
   * spin values, i.e. S_i. The spins being evolved are S_i/M^2. */
  for ( i = 0; i < 3; i++ )
  {
	 spin1.data[i]  *= mT2;
	 spin2.data[i]  *= mT2;
	 
     //s1normData[i] /= mT2;
     //s2normData[i] /= mT2;
  }

  /* Calculate various spin parameters */
  XLALSimIMRSpinEOBCalculateSigmaKerr( &sigmaKerr, eobParams->m1, 
				eobParams->m2, &spin1, &spin2 );
  XLALSimIMRSpinEOBCalculateSigmaStar( &sigmaStar, eobParams->m1, 
				eobParams->m2, &spin1, &spin2 );
  a = sqrt( sigmaKerr.data[0]*sigmaKerr.data[0] 
			+ sigmaKerr.data[1]*sigmaKerr.data[1]
            + sigmaKerr.data[2]*sigmaKerr.data[2] );
  //printf( "a = %e\n", a );
  //printf( "aStar = %e\n", sqrt( sigmaStar.data[0]*sigmaStar.data[0] + sigmaStar.data[1]*sigmaStar.data[1] + sigmaStar.data[2]*sigmaStar.data[2]) );
  if ( isnan( a ) )
  {
    printf( "a is nan in GSLSpinHamiltonianWrapperFordHdpphi !!\n");
      printf("rpolar, ppolar = %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f\n", rpolar[0], rpolar[1], rpolar[2], ppolar[0], ppolar[1], ppolar[2]);
      printf("rcart, pcart = %3.10f %3.10f %3.10f %3.10f %3.10f %3.10f\n", rcart[0], rcart[1], rcart[2], pcart[0], pcart[1], pcart[2]);
      XLALPrintError( "XLAL Error - %s: a = nan   \n", __func__);
      XLAL_ERROR( XLAL_EINVAL );
  }
  //XLALSimIMRCalculateSpinEOBHCoeffs( dParams->params->seobCoeffs, eobParams->eta, a );
  /* If computing the derivative w.r.t. the position vector, we need to 
   * pass p and not p* to the Hamiltonian. This is so because we want to 
   * hold p constant as we compute dH/dx. The way to do this is set the
   * tortoise flag = 2, and convert the momenta being evolved, i.e. p*, 
   * to p and pass that as input */
   // TODO
#if 0
  if ( dParams->varyParam < 3 )
  {
	  /* set the tortoise flag to 2 */
	  dParams->params->tortoise = 2;
	  
	  /* Convert momenta to p */
	  REAL8 tmpP[3]= {0.};
	  REAL8 rMag = sqrt(r.data[0]*r.data[0] + r.data[1]*r.data[1] 
					+ r.data[2]*r.data[2]);
	  REAL8 prT = p.data[0]*(r.data[0]/rMag) + p.data[1]*(r.data[1]/rMag) 
					+ p.data[2]*(r.data[2]/rMag);
					
	  REAL8 u, u2, u3, u4, u5, w2, a2;
	  REAL8 csi; 
	  
	  u  = 1./rMag;
      u2 = u*u;
      u3 = u2*u;
      u4 = u2*u2;
      u5 = u4*u;
      a2 = a*a;
      w2 = rMag*rMag + a2;
      /* Eq. 5.83 of BB1, inverse */
      REAL8 D = 1. + log(1. + 6.*eta*u2 + 2.*(26. - 3.*eta)*eta*u3);
      REAL8 m1PlusetaKK = -1. + eta * coeffs->KK;
	  /* Eq. 5.75 of BB1 */
      REAL8 bulk = 1./(m1PlusetaKK*m1PlusetaKK) + (2.*u)/m1PlusetaKK + a2*u2;
	  /* Eq. 5.73 of BB1 */
	  REAL8 logTerms = 1. + eta*coeffs->k0 + eta*log(fabs(1. + coeffs->k1*u
		+ coeffs->k2*u2 + coeffs->k3*u3 + coeffs->k4*u4
		+ coeffs->k5*u5 + coeffs->k5l*u5*log(u)));
	  /* Eq. 5.73 of BB1 */
      REAL8 deltaU = bulk*logTerms;
      deltaU = fabs(deltaU);
      /* Eq. 5.71 of BB1 */
      REAL8 deltaT = rMag*rMag*deltaU;
      /* Eq. 5.38 of BB1 */
      REAL8 deltaR = deltaT*D;	   
	  if ( oldTortoise )
		  csi = sqrt( fabs(deltaT * deltaR) )/ w2; /* Eq. 28 of Pan et al. PRD 81, 084041 (2010) */
	  else
	      csi = 1.0;

	  for( i = 0; i < 3; i++ )
	  { 
		  tmpP[i] = p.data[i] - (r.data[i]/rMag) * prT * (csi-1.)/csi;
	  }
	  memcpy( p.data, tmpP, sizeof(tmpP) );
  }
#endif
  //printf( "Hamiltonian = %e\n", XLALSimIMRSpinEOBHamiltonian( eobParams->eta, &r, &p, &sigmaKerr, &sigmaStar, dParams->params->seobCoeffs ) );
  REAL8 SpinEOBH = XLALSimIMRSpinEOBHamiltonian( eobParams->eta, &r, &p, &spin1norm, &spin2norm, &sigmaKerr, &sigmaStar, dParams->params->tortoise, dParams->params->seobCoeffs ) / eobParams->eta;
  
  return SpinEOBH;
}

#endif /*_LALSIMIMRSPINEOBHAMILTONIAN_C*/
