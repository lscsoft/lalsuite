/*
*  Copyright (C) 2010 Craig Robinson 
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
 * \author Craig Robinson
 *
 * \brief In newer versions of the EOBNR approximant, we
 * do not have an analytic expression for the derivative of the waveform.
 * As such, it is necessary to calculate the derivatives numerically. This
 * function provides the means to do just that.
 * 
 */
#include <unistd.h>

#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

#include "LALSimIMRSpinEOB.h"
#include "LALSimIMRSpinEOBHamiltonian.c"

#include <gsl/gsl_deriv.h>


static double GSLSpinAlignedHamiltonianWrapper( double x, void *params );

static int XLALSpinAlignedHcapDerivative(
                          double     UNUSED     t,
                          const REAL8           values[],
                          REAL8                 dvalues[],
                          void                  *funcParams
                               )
{

  static const REAL8 STEP_SIZE = 1.0e-4;

  static const INT4 lMax = 8;

  HcapDerivParams params;

  /* Since we take numerical derivatives wrt dynamical variables */
  /* but we want them wrt time, we use this temporary vector in  */
  /* the conversion */
  REAL8           tmpDValues[6];

  /* Cartesian values for calculating the Hamiltonian */
  REAL8           cartValues[6];

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

  /* Spins */
  REAL8Vector *sKerr = NULL;
  REAL8Vector *sStar = NULL;

  REAL8 a;

  REAL8 omega;

  /* EOB potential functions */
  REAL8 DeltaT, DeltaR;
  REAL8 csi;

  /* The error in a derivative as measured by GSL */
  REAL8 absErr;

  /* Set up pointers for GSL */ 
  params.values  = cartValues;
  params.params  = (SpinEOBParams *)funcParams;

  sKerr = params.params->sigmaKerr;
  sStar = params.params->sigmaStar;

  F.function = &GSLSpinAlignedHamiltonianWrapper;
  F.params   = &params;

  mass1 = params.params->eobParams->m1;
  mass2 = params.params->eobParams->m2;
  eta   = params.params->eobParams->eta;

  r = values[0];

  /* Since this is spin aligned, I make the assumption */
  /* that the spin vector is along the z-axis.         */
  a  = sKerr->data[2];

  /* Calculate the potential functions */
  DeltaT = XLALSimIMRSpinEOBHamiltonianDeltaT( params.params->seobCoeffs, r, eta, a );
  DeltaR = XLALSimIMRSpinEOBHamiltonianDeltaR( params.params->seobCoeffs, r, eta, a );
  csi    = sqrt( DeltaT * DeltaR ) / (r*r + a*a);
  //printf( "csi in derivatives function = %.16e\n", csi );

  /* Populate the Cartesian values vector */
  /* We can assume phi is zero wlog */
  memset( cartValues, 0, sizeof( cartValues ) );
  cartValues[0] = values[0];
  cartValues[3] = values[2];
  cartValues[4] = values[3] / values[0];

  /* Now calculate derivatives w.r.t. each parameter */
  for ( i = 0; i < 6; i++ )
  {
    params.varyParam = i;
    XLAL_CALLGSL( gslStatus = gsl_deriv_central( &F, cartValues[i], 
                    STEP_SIZE, &tmpDValues[i], &absErr ) );

    if ( gslStatus != GSL_SUCCESS )
    {
      XLALPrintError( "XLAL Error - %s: Failure in GSL function\n", __func__ );
      XLAL_ERROR( XLAL_EFUNC );
    }
  }

  /* Calculate the polar data */
  polarDynamics.length = 4;
  polarDynamics.data   = polData;

  memcpy( polData, values, sizeof( polData ) );

  rVec.length = pVec.length = 3;
  rVec.data   = rData;
  pVec.data   = pData;

  memset( rData, 0, sizeof(rData) );
  memset( pData, 0, sizeof(pData) );

  rData[0] = values[0];
  pData[0] = values[2];
  pData[1] = values[3] / values[0];

  H =  XLALSimIMRSpinEOBHamiltonian( eta, &rVec, &pVec, sKerr, sStar, params.params->tortoise, params.params->seobCoeffs ); 

  //printf( "csi = %.16e, ham = %.16e ( tortoise = %d)\n", csi, H, params.params->tortoise );
  //exit(1);
  //printf( "Hamiltonian = %e\n", H );
  H = H * (mass1 + mass2);


  //printf( "Cartesian derivatives:\n%.16e %.16e %.16e %.16e %.16e %.16e\n",
  //    tmpDValues[3], tmpDValues[4], tmpDValues[5], -tmpDValues[0], -tmpDValues[1], -tmpDValues[2] );

  /* Now calculate omega, and hence the flux */
  omega = tmpDValues[4] / r;
  flux  = XLALInspiralSpinFactorizedFlux( &polarDynamics, omega, params.params, H/(mass1+mass2), lMax );

  /* Looking at the non-spinning model, I think we need to divide the flux by eta */
  flux = flux / eta;

  //printf( "Flux in derivatives function = %.16e\n", flux );

  /* Now we can calculate the final (spherical) derivatives */
  /* csi is needed because we use the tortoise co-ordinate */
  dvalues[0] = csi * tmpDValues[3];
  dvalues[1] = omega;
  dvalues[2] = - tmpDValues[0] + tmpDValues[4] * values[3] / (r*r);
  dvalues[2] = dvalues[2] * csi - ( values[2] / values[3] ) * flux / omega;
  dvalues[3] = - flux / omega;

  //printf("Values:\n%.16e %.16e %.16e %.16e\n", values[0], values[1], values[2], values[3] );

  //printf("Derivatives:\n%.16e %.16e %.16e %.16e\n", dvalues[0], dvalues[1], dvalues[2], dvalues[3] );

  if ( isnan( dvalues[0] ) || isnan( dvalues[1] ) || isnan( dvalues[2] ) || isnan( dvalues[3] ) )
  {
    printf( "Deriv is nan: %e %e %e %e\n", dvalues[0], dvalues[1], dvalues[2], dvalues[3] );
    return 1;
  }

  return XLAL_SUCCESS;
}


/**
 * Function to calculate the value of omega for the spin-aligned EOB waveform
 */
REAL8
XLALSimIMRSpinAlignedEOBCalcOmega(
                      const REAL8           values[],
                      SpinEOBParams         *funcParams
                      )
{
  static const REAL8 STEP_SIZE = 1.0e-4;

  HcapDerivParams params;

  /* Cartesian values for calculating the Hamiltonian */
  REAL8 cartValues[6];

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

REAL8
XLALSimIMRSpinAlignedEOBNonKeplerCoeff(
                      const REAL8           values[],
                      SpinEOBParams         *funcParams
                      )
{

  REAL8 omegaCirc;

  REAL8 tmpValues[4];

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

  REAL8 tmpVec[6];

  /* These are the vectors which will be used in the call to the Hamiltonian */
  REAL8Vector r, p;
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

  return XLALSimIMRSpinEOBHamiltonian( eobParams->eta, &r, &p, sigmaKerr, sigmaStar, dParams->params->tortoise, dParams->params->seobCoeffs ) / eobParams->eta;
}

