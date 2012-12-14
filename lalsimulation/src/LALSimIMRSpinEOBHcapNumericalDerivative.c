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
 * \brief In newer versions of the EOBNR approximant, we
 * do not have an analytic expression for the derivative of the waveform.
 * As such, it is necessary to calculate the derivatives numerically. This
 * function provides the means to do just that.
 * 
 */

#ifndef _LALSIMIMRSPINEOBHCAPNUMERICALDERIVATIVE_C
#define _LALSIMIMRSPINEOBHCAPNUMERICALDERIVATIVE_C

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include <unistd.h>
#include <math.h>
#include <gsl/gsl_deriv.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

#include "LALSimIMRSpinEOB.h"

#include "LALSimIMRSpinEOBAuxFuncs.c"
#include "LALSimIMRSpinEOBHamiltonian.c"
#include "LALSimIMRSpinEOBFactorizedFlux.c"
#include "LALSimIMREOBFactorizedWaveform.c"

/*------------------------------------------------------------------------------------------
 *
 *          Prototypes of functions defined in this code.
 *
 *------------------------------------------------------------------------------------------
 */


static double GSLSpinHamiltonianWrapper( double x, void *params );

static int XLALSpinHcapNumericalDerivative(
                          double                t,
                          const REAL8           values[],
                          REAL8                 dvalues[],
                          void                  *funcParams
                               ) UNUSED;

static REAL8 XLALSpinHcapNumDerivWRTParam(
                       const INT4 paramIdx,
                       const REAL8 values[],
                       SpinEOBParams *params
                       );


/*------------------------------------------------------------------------------------------
 *
 *          Defintions of functions.
 *
 *------------------------------------------------------------------------------------------
 */

/**
 * Function to calculate numerical derivatives of the spin EOB Hamiltonian,
 * which correspond to time derivatives of the dynamical variables in conservative dynamcis.
 * All derivatives, including those on two terms of the orbital phase, are returned together.
 * The derivatives are combined with energy flux to give right hand side of the ODEs
 * of a generic spin EOB model, as decribed in Eqs. 21, 22, 26 and 27 of 
 * Pan et al. PRD 81, 084041 (2010)
 * This function is not used by the spin-aligned SEOBNRv1 model.
 */
static int XLALSpinHcapNumericalDerivative(
                 double UNUSED     t,         /**<< UNUSED */
                 const  REAL8      values[],  /**<< Dynamical variables */
                 REAL8             dvalues[], /**<< Time derivatives of variables (returned) */
                 void             *funcParams /**<< EOB parameters */
                               )
{

  static const REAL8 STEP_SIZE = 1.0e-4;

  static const INT4 lMax = 8;

  HcapDerivParams params;

  /* Since we take numerical derivatives wrt dynamical variables */
  /* but we want them wrt time, we use this temporary vector in  */
  /* the conversion */
  REAL8           tmpDValues[14];

  REAL8           H; //Hamiltonian
  REAL8           flux;

  gsl_function F;
  INT4         gslStatus;

  UINT4 i;

  REAL8Vector rVec, pVec;
  REAL8 rData[3], pData[3];

  /* We need r, phi, pr, pPhi to calculate the flux */
  REAL8       r;
  REAL8Vector polarDynamics;
  REAL8       polData[4];

  REAL8 mass1, mass2, eta;
  REAL8 rrTerm2, pDotS1, pDotS2;
  REAL8Vector s1, s2, sKerr, sStar;
  REAL8       s1Data[3], s2Data[3];
  REAL8       sKerrData[3], sStarData[3];
  REAL8 magS1, magS2, chiS, chiA, a;


  /* Orbital angular momentum */
  REAL8 Lx, Ly, Lz, magL;
  REAL8 Lhatx, Lhaty, Lhatz;
  REAL8 dLx, dLy, dLz;
  REAL8 dLhatx, dLhaty, dMagL;

  REAL8 alphadotcosi;

  REAL8 rCrossV_x, rCrossV_y, rCrossV_z, omega;

  /* The error in a derivative as measured by GSL */
  REAL8 absErr;

  /* Set up pointers for GSL */ 
  params.values  = values;
  params.params  = (SpinEOBParams *)funcParams;

  F.function = &GSLSpinHamiltonianWrapper;
  F.params   = &params;

  mass1 = params.params->eobParams->m1;
  mass2 = params.params->eobParams->m2;
  eta   = params.params->eobParams->eta;

  /* Now calculate derivatives w.r.t. each parameter */
  for ( i = 0; i < 12; i++ )
  {
    params.varyParam = i;
    if ( i >=6 && i < 9 )
    {
      XLAL_CALLGSL( gslStatus = gsl_deriv_central( &F, values[i],
                      STEP_SIZE*mass1*mass1, &tmpDValues[i], &absErr ) );
    }
    else if ( i >= 9 )
    {
      XLAL_CALLGSL( gslStatus = gsl_deriv_central( &F, values[i],
                      STEP_SIZE*mass2*mass2, &tmpDValues[i], &absErr ) );
    }
    else
    {
      XLAL_CALLGSL( gslStatus = gsl_deriv_central( &F, values[i], 
                      STEP_SIZE, &tmpDValues[i], &absErr ) );
    }
    if ( gslStatus != GSL_SUCCESS )
    {
      XLALPrintError( "XLAL Error %s - Failure in GSL function\n", __func__ );
      XLAL_ERROR( XLAL_EFUNC );
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

  r = polData[0] = sqrt( values[0]*values[0] + values[1]*values[1] + values[2]*values[2] );
  polData[1] = 0;
  polData[2] = values[0]*values[3] + values[1]*values[4] + values[2]*values[5];
  polData[3] = magL / polData[0];

  /* We need to re-calculate the parameters at each step as spins may not be constant */
  /* TODO: Modify so that only spin terms get re-calculated */

  /* We cannot point to the values vector directly as it leads to a warning */
  s1.length = s2.length = 3;
  s1.data = s1Data;
  s2.data = s2Data;

  memcpy( s1Data, values+6, 3*sizeof(REAL8) );
  memcpy( s2Data, values+9, 3*sizeof(REAL8) );

  magS1 = sqrt(s1.data[0]*s1.data[0] + s1.data[1]*s1.data[1] + s1.data[2]*s1.data[2]);
  magS2 = sqrt(s2.data[0]*s2.data[0] + s2.data[1]*s2.data[1] + s2.data[2]*s2.data[2]);

  chiS = 0.5 * ( magS1 / (mass1*mass1) + magS2 / (mass2*mass2) );
  chiA = 0.5 * ( magS1 / (mass1*mass1) - magS2 / (mass2*mass2) );

  sKerr.length = 3;
  sKerr.data   = sKerrData; 
  XLALSimIMRSpinEOBCalculateSigmaKerr( &sKerr, mass1, mass2, &s1, &s2 );

  sStar.length = 3;
  sStar.data   = sStarData;
  XLALSimIMRSpinEOBCalculateSigmaStar( &sStar, mass1, mass2, &s1, &s2 );

  a = sqrt(sKerr.data[0]*sKerr.data[0] + sKerr.data[1]*sKerr.data[1] + sKerr.data[2]*sKerr.data[2]);

  XLALSimIMREOBCalcSpinFacWaveformCoefficients( params.params->eobParams->hCoeffs, mass1, mass2, eta, a, chiS, chiA );
  XLALSimIMRCalculateSpinEOBHCoeffs( params.params->seobCoeffs, eta, a );
 
  rVec.length = pVec.length = 3;
  rVec.data   = rData;
  pVec.data   = pData;

  memcpy( rData, values, sizeof(rData) );
  memcpy( pData, values+3, sizeof(pData) );

  H =  XLALSimIMRSpinEOBHamiltonian( eta, &rVec, &pVec, &sKerr, &sStar, params.params->tortoise, params.params->seobCoeffs ); 


  /* Now make the conversion */
  /* rDot */
  dvalues[0]  = tmpDValues[3];
  dvalues[1]  = tmpDValues[4];
  dvalues[2]  = tmpDValues[5];

  /* Now calculate omega, and hence the flux */
  rCrossV_x = values[1]*dvalues[2] - values[2]*dvalues[1];
  rCrossV_y = values[2]*dvalues[0] - values[0]*dvalues[2];
  rCrossV_z = values[0]*dvalues[1] - values[1]*dvalues[0];

  omega = sqrt( rCrossV_x*rCrossV_x + rCrossV_y*rCrossV_y + rCrossV_z*rCrossV_z ) / (r*r);
  flux  = XLALInspiralSpinFactorizedFlux( &polarDynamics, omega, params.params, H/(mass1+mass2), lMax );

  /* Looking at the non-spinning model, I think we need to divide the flux by eta */
  flux = flux / eta;

  pDotS1 = pData[0]*s1Data[0] + pData[1]*s1Data[1] + pData[2]*s1Data[2];
  pDotS2 = pData[0]*s2Data[0] + pData[1]*s2Data[1] + pData[2]*s2Data[2];
  rrTerm2 = 8./15. *eta*eta * pow(omega,8./3.)/(magL*magL*r) * ((61.+48.*mass2/mass1)*pDotS1 + (61.+48.*mass1/mass2)*pDotS2);

  //printf( "rrForce = %e %e %e\n", - flux * values[3] / (omega*magL), - flux * values[4] / (omega*magL), - flux * values[5] / (omega*magL)) ;

  /* Now pDot */
  dvalues[3]  = - tmpDValues[0] - flux * values[3] / (omega*magL) + rrTerm2*Lx;
  dvalues[4]  = - tmpDValues[1] - flux * values[4] / (omega*magL) + rrTerm2*Ly;
  dvalues[5]  = - tmpDValues[2] - flux * values[5] / (omega*magL) + rrTerm2*Lz;

  /* spin1 */
  //printf( "Raw spin1 derivatives = %e %e %e\n", tmpDValues[6], tmpDValues[7], tmpDValues[8] );
  //printf( "Raw spin2 derivatives = %e %e %e\n", tmpDValues[9], tmpDValues[10], tmpDValues[11] );
  dvalues[6]  = mass1 * mass2 * (tmpDValues[7]*values[8] - tmpDValues[8]*values[7]);
  dvalues[7]  = mass1 * mass2 * (tmpDValues[8]*values[6] - tmpDValues[6]*values[8]);
  dvalues[8]  = mass1 * mass2 * (tmpDValues[6]*values[7] - tmpDValues[7]*values[6]);

  /* spin2 */
  dvalues[9]  = mass1 * mass2 * (tmpDValues[10]*values[11] - tmpDValues[11]*values[10]);
  dvalues[10] = mass1 * mass2 * (tmpDValues[11]*values[9] - tmpDValues[9]*values[11]);
  dvalues[11] = mass1 * mass2 * (tmpDValues[9]*values[10] - tmpDValues[10]*values[9]);

  /* phase and precessing bit */
  dLx = dvalues[1]*values[5] - dvalues[2]*values[4]
      + values[1]*dvalues[5] - values[2]*dvalues[4];

  dLy = dvalues[2]*values[3] - dvalues[0]*values[5]
      + values[2]*dvalues[3] - values[0]*dvalues[5];

  dLz = dvalues[0]*values[4] - dvalues[1]*values[3]
      + values[0]*dvalues[4] - values[1]*dvalues[3];

  dMagL = (Lx*dLx + Ly*dLy + Lz*dLz)/magL;

  dLhatx = (dLx*magL - Lx*dMagL)/(magL*magL);
  dLhaty = (dLy*magL - Ly*dMagL)/(magL*magL);
  
  /* Finn Chernoff convention is used here. TODO: implement the geometric precessing convention */
  if ( Lhatx == 0.0 && Lhaty == 0.0 )
  {
    alphadotcosi = 0.0;
  }
  else
  {
    alphadotcosi = Lhatz * (Lhatx*dLhaty - Lhaty*dLhatx) / (Lhatx*Lhatx + Lhaty*Lhaty);
  }

  dvalues[12] = omega - alphadotcosi;
  dvalues[13] = alphadotcosi;

  /*printf( " r = %e %e %e (mag = %e)\n", values[0], values[1], values[2], sqrt(values[0]*values[0] + values[1]*values[1] + values[2]*values[2]));
  printf( " p = %e %e %e (mag = %e)\n", values[3], values[4], values[5], sqrt(values[3]*values[3] + values[4]*values[4] + values[5]*values[5]));
  printf( "Derivatives:\n" );
  for ( i = 0; i < 12; i++ )
  {
    printf( "\t%e", dvalues[i] );
  }
  printf( "\n" );
  */
  return XLAL_SUCCESS;
}


/**
 * Calculate the derivative of the Hamiltonian w.r.t. a specific parameter
 * Used by generic spin EOB model, including initial conditions solver.
 */
static REAL8 XLALSpinHcapNumDerivWRTParam(
                 const INT4 paramIdx,      /**<< Index of the parameters */
                 const REAL8 values[],     /**<< Dynamical variables */
                 SpinEOBParams *funcParams /**<< EOB Parameters */
                 )
{
  static const REAL8 STEP_SIZE = 1.0e-3;

  HcapDerivParams params;

  REAL8 result;

  gsl_function F;
  INT4         gslStatus;

  REAL8 mass1, mass2;

  /* The error in a derivative as measured by GSL */
  REAL8 absErr;

  /* Set up pointers for GSL */
  params.values  = values;
  params.params  = funcParams;

  F.function       = &GSLSpinHamiltonianWrapper;
  F.params         = &params;
  params.varyParam = paramIdx;

  mass1 = params.params->eobParams->m1;
  mass2 = params.params->eobParams->m2;

  /* Now calculate derivatives w.r.t. the required parameter */
  if ( paramIdx >=6 && paramIdx < 9 )
  {
    XLAL_CALLGSL( gslStatus = gsl_deriv_central( &F, values[paramIdx],
                    STEP_SIZE*mass1*mass1, &result, &absErr ) );
  }
  else if ( paramIdx >= 9 )
  {
    XLAL_CALLGSL( gslStatus = gsl_deriv_central( &F, values[paramIdx],
                    STEP_SIZE*mass2*mass2, &result, &absErr ) );
  }
  else
  {
    XLAL_CALLGSL( gslStatus = gsl_deriv_central( &F, values[paramIdx],
                    STEP_SIZE, &result, &absErr ) );
  }
  if ( gslStatus != GSL_SUCCESS )
  {
    XLALPrintError( "XLAL Error %s - Failure in GSL function\n", __func__ );
    XLAL_ERROR_REAL8( XLAL_EFUNC );
  }

  //printf( "Abserr = %e\n", absErr );

  return result;
}


/** 
 * Wrapper for GSL to call the Hamiltonian function 
 */
static double GSLSpinHamiltonianWrapper( double x, void *params )
{
  HcapDerivParams *dParams = (HcapDerivParams *)params;

  EOBParams *eobParams = dParams->params->eobParams;

  REAL8 tmpVec[12];
  REAL8 sKerrData[3], sStarData[3];

  /* These are the vectors which will be used in the call to the Hamiltonian */
  REAL8Vector r, p, spin1, spin2;
  REAL8Vector sigmaKerr, sigmaStar;

  REAL8 a;

  /* Use a temporary vector to avoid corrupting the main function */
  memcpy( tmpVec, dParams->values, 
               sizeof(tmpVec) );

  /* Set the relevant entry in the vector to the correct value */
  tmpVec[dParams->varyParam] = x;

  /* Set the LAL-style vectors to point to the appropriate things */
  r.length = p.length = spin1.length = spin2.length = 3;
  sigmaKerr.length = sigmaStar.length = 3;
  r.data     = tmpVec;
  p.data     = tmpVec+3;
  spin1.data = tmpVec+6;
  spin2.data = tmpVec+9;
  sigmaKerr.data = sKerrData;
  sigmaStar.data = sStarData;

  /* Calculate various spin parameters */
  XLALSimIMRSpinEOBCalculateSigmaKerr( &sigmaKerr, eobParams->m1, eobParams->m2, &spin1, &spin2 );
  XLALSimIMRSpinEOBCalculateSigmaStar( &sigmaStar, eobParams->m1, eobParams->m2, &spin1, &spin2 );

  a = sqrt( sigmaKerr.data[0]*sigmaKerr.data[0] + sigmaKerr.data[1]*sigmaKerr.data[1]
              + sigmaKerr.data[2]*sigmaKerr.data[2] );
  //printf( "a = %e\n", a );
  //printf( "aStar = %e\n", sqrt( sigmaStar.data[0]*sigmaStar.data[0] + sigmaStar.data[1]*sigmaStar.data[1] + sigmaStar.data[2]*sigmaStar.data[2]) );
  if ( isnan( a ) )
  {
    printf( "a is nan!!\n");
  }
  XLALSimIMRCalculateSpinEOBHCoeffs( dParams->params->seobCoeffs, eobParams->eta, a );

  //printf( "Hamiltonian = %e\n", XLALSimIMRSpinEOBHamiltonian( eobParams->eta, &r, &p, &sigmaKerr, &sigmaStar, dParams->params->seobCoeffs ) );
  return XLALSimIMRSpinEOBHamiltonian( eobParams->eta, &r, &p, &sigmaKerr, &sigmaStar, dParams->params->tortoise, dParams->params->seobCoeffs ) / eobParams->eta;
}

#endif /* _LALSIMIMRSPINEOBHCAPNUMERICALDERIVATIVE_C */
