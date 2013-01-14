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
 * function provides the means to do that for the SEOBNRv1 Hamiltonian
 * in the spin-aligned case, i.e. equatorial orbits with four dynamical variables:
 * r, phi, pr and pphi. It then combines energy flux and numerical derivatives of 
 * the Hamiltonian to calculate the right hand side of the ODEs given by
 * Eqs. 10a - 10d of Pan et al. PRD 84, 124052 (2011). 
 * Since SEOBNRv1 is a spin-aligned model, the ODEs are written 
 * in the same format as ODEs of a nonspinning model.
 * 
 */

#ifndef LALSIMIMRSPINALIGNEDEOBHCAPDERIVATIVE_C
#define LALSIMIMRSPINALIGNEDEOBHCAPDERIVATIVE_C

#include <unistd.h>

#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

#include "LALSimIMRSpinEOB.h"
#include "LALSimIMRSpinEOBHamiltonian.c"
#include "LALSimIMRSpinEOBFactorizedFlux.c"

#include <gsl/gsl_deriv.h>

/*------------------------------------------------------------------------------------------
 *
 *          Prototypes of functions defined in this code.
 *
 *------------------------------------------------------------------------------------------
 */

static int XLALSpinAlignedHcapDerivative(
                          double                t,
                          const REAL8           values[],
                          REAL8                 dvalues[],
                          void                  *funcParams
                               );
/*------------------------------------------------------------------------------------------
 *
 *          Defintions of functions.
 *
 *------------------------------------------------------------------------------------------
 */

/**
 * Function to calculate R.H.S. of the ODEs, given dyanmical variables,
 * their derivatives and EOB parameters. Since SEOBNRv1 spin Hamiltonian
 * was implemented for Cartesean coordinates while dynamical evolution was 
 * implemented in polar coordinates, we need to perform a transform. 
 * This is done in a particular transform in which
 * x = r, y = z = 0, px = pr, py = pphi/r, pz = 0, and
 * omega = v/r = (dy/dt)/r = (dH/dpy)/r, dr/dt = dx/dt = dH/dpx, etc.
 */
static int XLALSpinAlignedHcapDerivative(
                  double UNUSED t,          /**< UNUSED */
                  const REAL8   values[],   /**< dynamical varables */
                  REAL8         dvalues[],  /**< time derivative of dynamical variables */
                  void         *funcParams  /**< EOB parameters */
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

  /* Calculate the potential functions and the tortoise coordinate factor csi,
     given by Eq. 28 of Pan et al. PRD 81, 084041 (2010) */
  DeltaT = XLALSimIMRSpinEOBHamiltonianDeltaT( params.params->seobCoeffs, r, eta, a );
  DeltaR = XLALSimIMRSpinEOBHamiltonianDeltaR( params.params->seobCoeffs, r, eta, a );
  csi    = sqrt( DeltaT * DeltaR ) / (r*r + a*a);
  //printf( "csi in derivatives function = %.16e\n", csi );

  /* Populate the Cartesian values vector, using polar coordinate values */
  /* We can assume phi is zero wlog */
  memset( cartValues, 0, sizeof( cartValues ) );
  cartValues[0] = values[0];
  cartValues[3] = values[2];
  cartValues[4] = values[3] / values[0];

  /* Now calculate derivatives w.r.t. each Cartesian variable */
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

  /* Calculate the Cartesian vectors rVec and pVec */
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
  /* Calculate Hamiltonian using Cartesian vectors rVec and pVec */
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
  /* Right hand side of Eqs. 10a - 10d of Pan et al. PRD 84, 124052 (2011) */
  dvalues[0] = csi * tmpDValues[3];
  dvalues[1] = omega;
  /* Note: in this special coordinate setting, namely y = z = 0, dpr/dt = dpx/dt + dy/dt * py/r, where py = pphi/r */ 
  dvalues[2] = - tmpDValues[0] + tmpDValues[4] * values[3] / (r*r);
  dvalues[2] = dvalues[2] * csi - ( values[2] / values[3] ) * flux / omega;
  dvalues[3] = - flux / omega;

  //printf("Values:\n%.16e %.16e %.16e %.16e\n", values[0], values[1], values[2], values[3] );

  //printf("Derivatives:\n%.16e %.16e %.16e %.16e\n", dvalues[0], dvalues[1], dvalues[2], dvalues[3] );

  if ( isnan( dvalues[0] ) || isnan( dvalues[1] ) || isnan( dvalues[2] ) || isnan( dvalues[3] ) )
  {
    //printf( "Deriv is nan: %e %e %e %e\n", dvalues[0], dvalues[1], dvalues[2], dvalues[3] );
    return 1;
  }

  return XLAL_SUCCESS;
}

#endif /* LALSIMIMRSPINALIGNEDEOBHCAPDERIVATIVE_C */
