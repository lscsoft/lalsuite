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

#ifndef _LALSIMIMRSPINEOBHAMILTONIANOPTIMIZED_C
#define _LALSIMIMRSPINEOBHAMILTONIANOPTIMIZED_C

#include <stdio.h>
#include <math.h>

#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

#include "LALSimIMRSpinEOB.h"

#include "LALSimIMRSpinEOBHcapExactDerivative.c"

/*------------------------------------------------------------------------------------------
 *
 *          Prototypes of functions defined in this code.
 *
 *------------------------------------------------------------------------------------------
 */


static REAL8 XLALSimIMRSpinEOBHamiltonianOptimized(
               const REAL8    eta,
               REAL8Vector    * restrict x,
               REAL8Vector    * restrict p,
               REAL8Vector    * restrict s1Vec,
               REAL8Vector    * restrict s2Vec,
               REAL8Vector    * restrict sigmaKerr,
               REAL8Vector    * restrict sigmaStar,
               int                       tortoise,
               SpinEOBHCoeffs *coeffs);

static REAL8 XLALSimIMRSpinAlignedEOBCalcOmegaOptimized(
                      const REAL8          values[],
                      SpinEOBParams        *funcParams
                      );

static REAL8 XLALSimIMRSpinAlignedEOBNonKeplerCoeffOptimized(
                      const REAL8           values[],
                      SpinEOBParams         *funcParams
                      );

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
REAL8 XLALSimIMRSpinEOBHamiltonianOptimized(
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
  /* This is a modified sign function: e3z = 1 if sigmaKerr->data[2]>=0, -1 otherwise */
  REAL8 e3z = (0.0 <= sigmaKerr->data[2]) - (sigmaKerr->data[2] < 0.0);

  if(tortoise==1) {
    {
#include "mathematica_codes/SEOBNRv2_opt_tortoise.h"
      return Hreal;
    }
  } else if(tortoise==0) {
    {
#include "mathematica_codes/SEOBNRv2_opt_.h"
      return Hreal;
    }
  } else {
    XLALPrintError( "XLAL Error - %s: SEOBNRv2 Optimized Hamiltonian called with unknown tortoise value = %d.\n", __func__,tortoise);
    XLAL_ERROR( XLAL_EINVAL );
  }
}

/**
 * Function to calculate the value of omega for the spin-aligned EOB waveform.
 * Can NOT be used in precessing cases. This omega is defined as \f$\dot{y}/r\f$ by setting \f$y=0\f$.
 * The function calculates omega = v/r, by first converting (r,phi,pr,pphi) to Cartesian coordinates
 * in which rVec={r,0,0} and pVec={0,pphi/r,0}, i.e. the effective-test-particle is positioned at x=r,
 * and its velocity along y-axis. Then it computes omega, which is now given by dydt/r = (dH/dp_y)/r.
 */
static REAL8
XLALSimIMRSpinAlignedEOBCalcOmegaOptimized(
                      const REAL8           values[],   /**<< Dynamical variables */
                      SpinEOBParams         *funcParams /**<< EOB parameters */
                      )
{
  HcapDerivParams params;

  /* Cartesian values for calculating the Hamiltonian */
  REAL8 cartValues[6];

  REAL8 omega;
  REAL8 r;

  /* Populate the Cartesian values vector */
  /* We can assume phi is zero wlog */
  memset( cartValues, 0, sizeof( cartValues ) );
  cartValues[0] = r = values[0];
  cartValues[3] = values[2];
  cartValues[4] = values[3] / values[0];

  /* Set up input parameters for exact derivative calculation */
  params.values  = cartValues;
  params.params  = funcParams;

  /* Now calculate omega. In the chosen co-ordinate system, */
  /* we need dH/dpy to calculate this, i.e. varyParam = 4   */
  params.varyParam = 4;
  omega = GSLSpinAlignedHamiltonianWrapper_ExactDeriv( cartValues[4], &params );

  omega = omega / r;

  return omega;
}

/**
 * Function to calculate the non-Keplerian coefficient for the spin-aligned EOB model.
 * radius \f$r\f$ times the cuberoot of the returned number is \f$r_\Omega\f$ defined in Eq. A2.
 * i.e. the function returns \f$(r_{\Omega} / r)^3\f$.
 */
static REAL8 XLALSimIMRSpinAlignedEOBNonKeplerCoeffOptimized(
                      const REAL8           values[],   /**<< Dynamical variables */
                      SpinEOBParams         *funcParams /**<< EOB parameters */
                      )
{

  REAL8 omegaCirc;

  REAL8 tmpValues[4];

  REAL8 r3;

  /* We need to find the values of omega assuming pr = 0 */
  memcpy( tmpValues, values, sizeof(tmpValues) );
  tmpValues[2] = 0.0;

  omegaCirc = XLALSimIMRSpinAlignedEOBCalcOmegaOptimized( tmpValues, funcParams );
  if ( XLAL_IS_REAL8_FAIL_NAN( omegaCirc ) )
  {
    XLAL_ERROR_REAL8( XLAL_EFUNC );
  }

  r3 = values[0]*values[0]*values[0];

  return 1.0/(omegaCirc*omegaCirc*r3);
}

#endif /*_LALSIMIMRSPINEOBHAMILTONIANOPTIMIZED_C*/
