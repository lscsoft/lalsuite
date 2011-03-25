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
 * \brief Function to compute the factorized flux as uses in the new EOBNR_PP
 * model. Flux function given by Phys.Rev.D79:064004,2009.
 */

#include <complex.h>

#include <lal/LALComplex.h>
#include <lal/LALInspiral.h>
#include <lal/LALEOBNRv2Waveform.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_linalg.h>


static inline
REAL8 GetNRPeakAmplitude( REAL8 eta )
{
  return eta * ( 1.43494 + 0.115275 * eta + 1.78719 * eta * eta );
}

static inline
REAL8 GetNRPeakADDot( REAL8 eta )
{
  return -0.01 * eta * ( 0.260034 + 0.0450257 * eta + 2.10027 * eta * eta );
}

static inline 
REAL8 GetNRPeakOmega( REAL8 eta )
{
  return 0.273207 + 0.240261 * eta + 0.425017 * eta * eta;
}

static inline 
REAL8 GetNRPeakOmegaDot( REAL8 eta )
{
  return 0.00586032 + 0.0136607 * eta + 0.0317176 * eta * eta;
}

int  XLALEOBNonQCCorrection(
                      COMPLEX16             * restrict nqc,
                      REAL8Vector           * restrict values,
                      REAL8Vector           * restrict dvalues,
                      EOBNonQCCoeffs        * restrict coeffs
                     )

{

  REAL8 rOmega, rOmegaSq;
  REAL8 r, p;

  REAL8 mag, phase;


  r = values->data[0];
  p = values->data[2];

  rOmega = r * dvalues->data[1];
  rOmegaSq = rOmega*rOmega;

  mag = 1. + (p*p / rOmegaSq) * ( coeffs->a1
     + coeffs->a2 / r + coeffs->a3 / (r*sqrt(r))
     + coeffs->a4 / (r*r) );

  phase = coeffs->b1 * p / rOmega + coeffs->b2 * p*p*p/rOmega;

  nqc->re = mag * cos(phase);
  nqc->im = mag * sin(phase);

  return XLAL_SUCCESS;

}


int XLALCalculateNQCCoefficients(
                 REAL8Vector    * restrict amplitude,
                 REAL8Vector    * restrict phase,
                 REAL8Vector    * restrict q1,
                 REAL8Vector    * restrict q2,
                 REAL8Vector    * restrict q3,
                 REAL8Vector    * restrict p1,
                 REAL8Vector    * restrict p2,
                 UINT4                     peakIdx,
                 REAL8                     deltaT,
                 REAL8                     eta,
                 EOBNonQCCoeffs * restrict coeffs )
{

  /* Code for 2,2 at the moment */
  /* Easily extended later */

  UINT4 i;

  /* For gsl perutation stuff */

  int signum;

  REAL8Vector * restrict time = NULL;

  /* Since the vectors we actually want are q etc * A, we will have to generate them here */
  REAL8Vector *q1LM = NULL;
  REAL8Vector *q2LM = NULL;
  REAL8Vector *q3LM = NULL; 

  REAL8 a, aDot, aDDot;
  REAL8 omega/*, omegaDot*/;

  /* Stuff for finding numerical derivatives */
  gsl_spline    *spline = NULL;
  gsl_interp_accel *acc = NULL;

  /* Matrix stuff for calculating coefficients */
  gsl_matrix *qMatrix = NULL, *pMatrix = NULL;
  gsl_vector *aCoeff  = NULL, *bCoeff  = NULL;

  gsl_vector *amps = NULL, *omegaVec = NULL;

  gsl_permutation *perm1 = NULL, *perm2 = NULL;

  /* Temporary to get rid of unused warning */
  p2->data[0] = p2->data[0];

  /* Populate the time vector */
  /* It is okay to assume initial t = 0 */
  time = XLALCreateREAL8Vector( q1->length );
  q1LM = XLALCreateREAL8Vector( q1->length );
  q2LM = XLALCreateREAL8Vector( q2->length );
  q3LM = XLALCreateREAL8Vector( q3->length );

  /* Populate vectors as necessary */
  for ( i = 0; i < time->length; i++ )
  {
    time->data[i] = i * deltaT;
    q1LM->data[i] = q1->data[i] * amplitude->data[i];
    q2LM->data[i] = q2->data[i] * amplitude->data[i];
    q3LM->data[i] = q3->data[i] * amplitude->data[i];
  }

  /* Allocate all the memory we need */
  XLAL_CALLGSL(
    /* a stuff */
    qMatrix = gsl_matrix_alloc( 3, 3 );
    aCoeff  = gsl_vector_alloc( 3 );
    amps    = gsl_vector_alloc( 3 );
    perm1   = gsl_permutation_alloc( 3 );

    /* b stuff */
    pMatrix  = gsl_matrix_alloc( 2, 2 );
    bCoeff   = gsl_vector_alloc( 2 );
    omegaVec = gsl_vector_alloc( 2 );
    perm2    = gsl_permutation_alloc( 2 );
  );

  if ( !qMatrix || !aCoeff || !amps || !pMatrix || !bCoeff || !omegaVec )
  {
    /* TODO : Free memory */
    XLAL_ERROR( __func__, XLAL_ENOMEM );
  }

  /* We are now in a position to use the interp stuff to calculate the derivatives we need */
  /* We will start with the quantities used in the calculation of the a coefficients */
  spline = gsl_spline_alloc( gsl_interp_cspline, amplitude->length );
  acc    = gsl_interp_accel_alloc();

  /* Q1 */
  gsl_spline_init( spline, time->data, q1LM->data, q1LM->length );
  gsl_matrix_set( qMatrix, 0, 0, gsl_spline_eval( spline, time->data[peakIdx], acc ) );
  gsl_matrix_set( qMatrix, 1, 0, gsl_spline_eval_deriv( spline, time->data[peakIdx], acc ) );
  gsl_matrix_set( qMatrix, 2, 0, gsl_spline_eval_deriv2( spline, time->data[peakIdx], acc ) );

  /* Q2 */
  gsl_spline_init( spline, time->data, q2LM->data, q2LM->length );
  gsl_interp_accel_reset( acc );
  gsl_matrix_set( qMatrix, 0, 1, gsl_spline_eval( spline, time->data[peakIdx], acc ) );
  gsl_matrix_set( qMatrix, 1, 1, gsl_spline_eval_deriv( spline, time->data[peakIdx], acc ) );
  gsl_matrix_set( qMatrix, 2, 1, gsl_spline_eval_deriv2( spline, time->data[peakIdx], acc ) );

  /* Q3 */
  gsl_spline_init( spline, time->data, q3LM->data, q3LM->length );
  gsl_interp_accel_reset( acc );
  gsl_matrix_set( qMatrix, 0, 2, gsl_spline_eval( spline, time->data[peakIdx], acc ) );
  gsl_matrix_set( qMatrix, 1, 2, gsl_spline_eval_deriv( spline, time->data[peakIdx], acc ) );
  gsl_matrix_set( qMatrix, 2, 2, gsl_spline_eval_deriv2( spline, time->data[peakIdx], acc ) );

  /* Amplitude */
  gsl_spline_init( spline, time->data, amplitude->data, amplitude->length );
  gsl_interp_accel_reset( acc );
  a     = gsl_spline_eval( spline, time->data[peakIdx], acc );
  aDot  = gsl_spline_eval_deriv( spline, time->data[peakIdx], acc );
  aDDot = gsl_spline_eval_deriv2( spline, time->data[peakIdx], acc );

  gsl_vector_set( amps, 0, GetNRPeakAmplitude( eta ) - a );
  gsl_vector_set( amps, 1, - aDot );
  gsl_vector_set( amps, 2, GetNRPeakADDot( eta ) - aDDot );

  /* We have now set up all the stuff to calculate the a coefficients */
  /* So let us do it! */
  gsl_linalg_LU_decomp( qMatrix, perm1, &signum );
  gsl_linalg_LU_solve( qMatrix, perm1, amps, aCoeff );

  /* Now we (should) have calculated the a values. Now we can do the b values */

  /* P1 */
/*  gsl_spline_init( spline, time->data, p1->data, p1->length );
  gsl_interp_accel_reset( acc );
  gsl_matrix_set( pMatrix, 0, 0, - gsl_spline_eval_deriv( spline, time->data[peakIdx], acc ) );
  gsl_matrix_set( pMatrix, 1, 0, - gsl_spline_eval_deriv2( spline, time->data[peakIdx], acc ) );
*/
  /* P2 */
/*  gsl_spline_init( spline, time->data, p2->data, p2->length );
  gsl_interp_accel_reset( acc );
  gsl_matrix_set( pMatrix, 0, 1, - gsl_spline_eval_deriv( spline, time->data[peakIdx], acc ) );
  gsl_matrix_set( pMatrix, 1, 1, - gsl_spline_eval_deriv2( spline, time->data[peakIdx], acc ) );
*/
  /* Phase */
  gsl_spline_init( spline, time->data, phase->data, phase->length );
  gsl_interp_accel_reset( acc );
  omega    = gsl_spline_eval_deriv( spline, time->data[peakIdx], acc );
/*  omegaDot = gsl_spline_eval_deriv2( spline, time->data[peakIdx], acc );*/

  /*gsl_vector_set( omegaVec, 0, GetNRPeakOmega( eta ) - omega );
  gsl_vector_set( omegaVec, 1, GetNRPeakOmegaDot( eta ) - omegaDot );
*/
  /* And now solve for the b coefficients */
/*  gsl_linalg_LU_decomp( pMatrix, perm2, &signum );
  gsl_linalg_LU_solve( pMatrix, perm2, omegaVec, bCoeff );
*/
  /* We can now populate the coefficients structure */
  coeffs->a1 = gsl_vector_get( aCoeff, 0 );
  coeffs->a2 = gsl_vector_get( aCoeff, 1 );
  coeffs->a3 = gsl_vector_get( aCoeff, 2 );
/*  coeffs->b1 = gsl_vector_get( bCoeff, 0 );
  coeffs->b2 = gsl_vector_get( bCoeff, 1 );
*/
  /* We think the b-fixing might be dodgy - set b2 to zero, and find b1 */
  coeffs->b2 = 0;

  gsl_spline_init( spline, time->data, p1->data, p1->length );
  gsl_interp_accel_reset( acc );

  coeffs->b1 = ( omega - GetNRPeakOmega( eta ) ) / gsl_spline_eval_deriv( spline, time->data[peakIdx], acc );

  /* Free memory and exit */
  gsl_matrix_free( qMatrix );
  gsl_vector_free( amps );
  gsl_vector_free( aCoeff );
  gsl_permutation_free( perm1 );

  gsl_matrix_free( pMatrix );
  gsl_vector_free( omegaVec );
  gsl_vector_free( bCoeff );
  gsl_permutation_free( perm2 );

  gsl_spline_free( spline );
  gsl_interp_accel_free( acc );

  XLALDestroyREAL8Vector( q1LM );
  XLALDestroyREAL8Vector( q2LM );
  XLALDestroyREAL8Vector( q3LM );
  XLALDestroyREAL8Vector( time );

  return XLAL_SUCCESS;
}
