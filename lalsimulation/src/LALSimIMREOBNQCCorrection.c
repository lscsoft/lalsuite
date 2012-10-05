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
 * \brief More recent versions of the EOB models, such as EOBNRv2 and SEOBNRv1, utilise
 * a non-quasicircular correction (NQC) to bring the peak of the EOB frequency
 * into agreement with that of NR simulations. This file contains the functions
 * used to calculate these NQC corrections, described in DCC document T1100433.
 * The fits to NR peak amplitude, frequency, and their derivatives, are taken 
 * from Pan et al. PRD 84 124052 (2011), for EOBNRv2, and
 * from Taracchini et al. PRD 86, 024011 (2012), for SEOBNRv1.
 */

#include <math.h>

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/LALComplex.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_linalg.h>

#include "LALSimIMREOBNRv2.h"

#ifndef _LALSIMIMRNQCCORRECTION_C 
#define _LALSIMIMRNQCCORRECTION_C 

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/* ------------------------------------------------
 *          Non-spin (EOBNRv2)
 * ------------------------------------------------*/

/**
 * Compute the time offset which should be used in computing the
 * non-quasicircular correction and performing the ringdown attachment.
 * These numbers were tuned to numerical relativity simulations, and
 * are taken from Pan et al, PRD84, 124052(2011), lines 1-5 of Table II. 
 */
 static REAL8 XLALSimIMREOBGetNRPeakDeltaT( 
                         INT4 l,    /**<< Mode l */ 
                         INT4 m,    /**<< Mode m */
                         REAL8 eta  /**<< Symmetric mass ratio */
                         )
{
  switch ( l )
  {
    case 2:
      if ( m == 2 )
      {
        return 0.0;
      }
      else if ( m == 1 )
      {
        return 10.67 - 2.5 + 9.0*eta - 41.41 * eta + 76.1 * eta*eta;
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 3:
      if ( m == 3 )
      {
        return 3.383 + 3.847 * eta + 8.979 * eta * eta;
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 4:
      if ( m == 4 )
      {
        return 5.57 - 49.86 * eta + 154.3 * eta * eta;
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 5:
      if ( m == 5 )
      {
        return 6.693 - 34.47 * eta + 102.7 * eta*eta;
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    default:
      XLAL_ERROR_REAL8( XLAL_EINVAL );
      break;
  }
}

/**
 * Function which returns a value of the expected peak amplitude
 * taken from a fit to numerical relativity simulations. The functions
 * are taken from Pan et al, PRD84, 124052(2011), lines 1-5 of Table II.  
 */
static inline
REAL8 GetNRPeakAmplitude( 
                        INT4 l,   /**<< Mode l */ 
                        INT4 m,   /**<< Mode m */
                        REAL8 eta /**<< Symmetric mass ratio */
                        )
{
  switch ( l )
  {
    case 2:
      if ( m == 2 )
      {
        return eta * ( 1.422 + 0.3013 * eta + 1.246 * eta * eta );
      }
      else if ( m == 1 )
      {
        return eta * sqrt( 1.0 - 4. * eta ) * (0.4832 - 0.01032 * eta);
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 3:
      if ( m == 3 )
      {
        return eta * sqrt(1.-4.*eta) * ( 0.5761 - 0.09638 * eta + 2.715*eta*eta );
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 4:
      if ( m == 4 )
      {
        return eta * (0.354 - 1.779 * eta + 2.834 * eta*eta );
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 5:
      if ( m == 5 )
      {
        return eta * sqrt(1.-4.*eta) * ( 0.1353 - 0.1485 * eta );
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    default:
      XLAL_ERROR_REAL8( XLAL_EINVAL );
      break;
  }
}

/**
 * Function which returns second derivative of the amplitude at the peak
 * taken from a fit to numerical relativity simulations. The functions
 * are taken from Pan et al, PRD84, 124052(2011), lines 1-5 of Table II. 
 */
static inline
REAL8 GetNRPeakADDot( 
                    INT4 l,   /**<< Mode l */ 
                    INT4 m,   /**<< Mode m */
                    REAL8 eta /**<< Symmetric mass ratio */
                    )
{
  switch ( l )
  {
    case 2:
      if ( m == 2 )
      {
        return -0.01 * eta * ( 0.1679 + 1.44 * eta - 2.001 * eta * eta );
      }
      else if ( m == 1 )
      {
        return -0.01 * eta * sqrt(1.-4.*eta) * (0.1867 + 0.6094 * eta );
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 3:
      if ( m == 3 )
      {
        return -0.01 * eta * sqrt(1.-4.*eta) * (0.2518 - 0.8145*eta + 5.731*eta*eta);
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 4:
      if ( m == 4 )
      {
        return -0.01 * eta * (0.1813 - 0.9935 * eta + 1.858 * eta * eta );
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 5:
      if ( m == 5 )
      {
        return -0.01 * eta * sqrt(1.-4.*eta) * ( 0.09051 - 0.1604 * eta );
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    default:
      XLAL_ERROR_REAL8( XLAL_EINVAL );
      break;
  }
}


/**
 * Function which returns a value of the expected peak frequency
 * taken from a fit to numerical relativity simulations. The functions
 * are taken from Pan et al, PRD84, 124052(2011), lines 1-5 of Table II. 
 */
static inline 
REAL8 GetNRPeakOmega( 
                    INT4 l,   /**<< Mode l */
                    INT4 m,   /**<< Mode m */
                    REAL8 eta /**<< Symmetric mass ratio */
                    )
{
  switch ( l )
  {
    case 2:
      if ( m == 2 )
      {
        return 0.2733 + 0.2316 * eta + 0.4463 * eta * eta;
      }
      else if ( m == 1 )
      {
        return 0.2907 - 0.08338 * eta + 0.587 * eta*eta;
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 3:
      if ( m == 3 )
      {
        return 0.4539 + 0.5376 * eta + 1.042 * eta*eta;
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 4:
      if ( m == 4 )
      {
        return 0.6435 - 0.05103 * eta + 2.216 * eta*eta;
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 5:
      if ( m == 5 )
      {
        return 0.8217 + 0.2346 * eta + 2.599 * eta*eta;
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    default:
      XLAL_ERROR_REAL8( XLAL_EINVAL );
      break;
  }
}

/**
 * Function which returns the derivative of the expected peak frequency
 * taken from a fit to numerical relativity simulations. The functions
 * are taken from Pan et al, PRD84, 124052(2011), lines 1-5 of Table II. 
 */
static inline 
REAL8 GetNRPeakOmegaDot( 
                       INT4 l,   /**<< Mode l */ 
                       INT4 m,   /**<< Mode m */
                       REAL8 eta /**<< Symmetric mass ratio */
                       )
{
  switch ( l )
  {
    case 2:
      if ( m == 2 )
      {
        return 0.005862 + 0.01506 * eta + 0.02625 * eta * eta;
      }
      else if ( m == 1 )
      {
        return 0.00149 + 0.09197 * eta - 0.1909 * eta*eta;
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 3:
      if ( m == 3 )
      {
        return 0.01074 + 0.0293 * eta + 0.02066 * eta*eta;
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 4:
      if ( m == 4 )
      {
        return 0.01486 + 0.08529 * eta - 0.2174 * eta * eta;
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 5:
      if ( m == 5 )
      {
        return 0.01775 + 0.09801 * eta - 0.1686 * eta*eta;
      }
      else
      {
        XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    default:
      XLAL_ERROR_REAL8( XLAL_EINVAL );
      break;
  }
}


/**
 * For the 2,2 mode, there are fits available for the NQC coefficients,
 * given in Eqs.(40a)-(40c) of Pan et al, PRD84, 124052(2011).
 * This function provides the values of these coefficients, so the 
 * correction can be used in the dynamics prior to finding the more
 * accurate NQC values later on.
 */
UNUSED static int XLALSimIMREOBGetCalibratedNQCCoeffs( 
                                EOBNonQCCoeffs *coeffs, /**<< OUTPUT, Structure for NQC coeffs */
                                INT4            l,      /**<< Mode l */
                                INT4            m,      /**<< Mode m */
                                REAL8           eta     /**<< Symmetric mass ratio */
                                )
{

#ifndef LAL_NDEBUG
  if ( !coeffs )
  {
    XLAL_ERROR( XLAL_EINVAL );
  }
#endif

  if ( l != 2 || m != 2 )
  {
    XLALPrintError( "Mode %d,%d is not supported by this function.\n", l, m );
    XLAL_ERROR( XLAL_EINVAL );
  }

  /* All NQC coefficients are set to zero here */
  /* including coeffs->a3S, coeffs->a4 and coeffs->a5 that are not used in EOBNRv2 */
  memset( coeffs, 0, sizeof( *coeffs ) );

  coeffs->a1 = -4.55919 + 18.761 * eta - 24.226 * eta*eta;
  coeffs->a2 = 37.683 - 201.468 * eta + 324.591 * eta*eta;
  coeffs->a3 = - 39.6024 + 228.899 * eta - 387.222 * eta * eta;

  return XLAL_SUCCESS;
}

/**
 * This function calculates the non-quasicircular correction to apply to 
 * the waveform. The form of this correction can be found in Pan et al, 
 * PRD84, 124052(2011), Eq.(22), and also in the DCC document T1100433. Note
 * that when calling this function, the NQC coefficients should already 
 * have been pre-computed.
 */
UNUSED static int  XLALSimIMREOBNonQCCorrection(
                      COMPLEX16      * restrict nqc,    /**<< OUTPUT, The NQC correction */
                      REAL8Vector    * restrict values, /**<< Dynamics r, phi, pr, pphi */
                      const REAL8               omega,  /**<< Angular frequency */
                      EOBNonQCCoeffs * restrict coeffs  /**<< NQC coefficients */
                     )

{

  REAL8 rOmega, rOmegaSq;
  REAL8 r, p, sqrtR;

  REAL8 mag, phase;


  r = values->data[0];
  p = values->data[2];

  sqrtR = sqrt(r);

  rOmega = r * omega;
  rOmegaSq = rOmega*rOmega;

  /* In EOBNRv2, coeffs->a3S, coeffs->a4 and coeffs->a5 are set to zero */
  /* through XLALSimIMREOBGetCalibratedNQCCoeffs() */
  /* and XLALSimIMREOBCalculateNQCCoefficients() */
  mag = 1. + (p*p / rOmegaSq) * ( coeffs->a1
     + coeffs->a2 / r + ( coeffs->a3 + coeffs->a3S) / (r*sqrtR)
     + coeffs->a4 / (r*r) + coeffs->a5 / (r*r*sqrtR));

  phase = coeffs->b1 * p / rOmega + p*p*p/rOmega * ( coeffs->b2
     + coeffs->b3 / sqrtR + coeffs->b4 / r );

  nqc->re = mag * cos(phase);
  nqc->im = mag * sin(phase);

  return XLAL_SUCCESS;

}

/**
 * This function computes the coefficients a1, a2, etc. used in the
 * non-quasicircular correction. The details of the calculation of these
 * coefficients are found in the DCC document T1100433. 
 */
UNUSED static int XLALSimIMREOBCalculateNQCCoefficients(
                 EOBNonQCCoeffs * restrict coeffs,    /**<< OUTPUT, NQC coefficients */
                 REAL8Vector    * restrict amplitude, /**<< Waveform amplitude, func of time */
                 REAL8Vector    * restrict phase,     /**<< Waveform phase(rad), func of time */
                 REAL8Vector    * restrict q1,        /**<< Function of dynamics (see DCC doc) */
                 REAL8Vector    * restrict q2,        /**<< Function of dynamics (see DCC doc) */
                 REAL8Vector    * restrict q3,        /**<< Function of dynamics (see DCC doc) */
                 REAL8Vector    * restrict p1,        /**<< Function of dynamics (see DCC doc) */
                 REAL8Vector    * restrict p2,        /**<< Function of dynamics (see DCC doc) */
                 INT4                      l,         /**<< Mode l */
                 INT4                      m,         /**<< Mode m */
                 REAL8                     timePeak,  /**<< Time of peak orbital frequency */
                 REAL8                     deltaT,    /**<< Sampling interval */
                 REAL8                     eta        /**<< Symmetric mass ratio */
                 )
{

  UINT4 i;

  int signum;

  REAL8Vector * restrict timeVec = NULL;

  /* Since the vectors we actually want are q etc * A, we will have to generate them here */
  REAL8Vector *q1LM = NULL;
  REAL8Vector *q2LM = NULL;
  REAL8Vector *q3LM = NULL; 

  REAL8 a, aDot, aDDot;
  REAL8 omega, omegaDot;

  REAL8 nra, nraDDot;
  REAL8 nromega, nromegaDot;

  REAL8 nrDeltaT, nrTimePeak;

  /* Stuff for finding numerical derivatives */
  gsl_spline    *spline = NULL;
  gsl_interp_accel *acc = NULL;

  /* Matrix stuff for calculating coefficients */
  gsl_matrix *qMatrix = NULL, *pMatrix = NULL;
  gsl_vector *aCoeff  = NULL, *bCoeff  = NULL;

  gsl_vector *amps = NULL, *omegaVec = NULL;

  gsl_permutation *perm1 = NULL, *perm2 = NULL;

  /* All NQC coefficients are set to zero here */ 
  /* including coeffs->a4 that is not used in EOBNRv2 */
  memset( coeffs, 0, sizeof( EOBNonQCCoeffs ) );

  /* Populate the time vector */
  /* It is okay to assume initial t = 0 */
  timeVec = XLALCreateREAL8Vector( q1->length );
  q1LM    = XLALCreateREAL8Vector( q1->length );
  q2LM    = XLALCreateREAL8Vector( q2->length );
  q3LM    = XLALCreateREAL8Vector( q3->length );

  /* Populate vectors as necessary */
  for ( i = 0; i < timeVec->length; i++ )
  {
    timeVec->data[i] = i * deltaT;
    q1LM->data[i]    = q1->data[i] * amplitude->data[i];
    q2LM->data[i]    = q2->data[i] * amplitude->data[i];
    q3LM->data[i]    = q3->data[i] * amplitude->data[i];
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
    gsl_matrix_free( qMatrix );
    gsl_vector_free( amps );
    gsl_vector_free( aCoeff );
    gsl_permutation_free( perm1 );
    gsl_matrix_free( pMatrix );
    gsl_vector_free( omegaVec );
    gsl_vector_free( bCoeff );
    gsl_permutation_free( perm2 );
    XLALDestroyREAL8Vector( q1LM );
    XLALDestroyREAL8Vector( q2LM );
    XLALDestroyREAL8Vector( q3LM );
    XLALDestroyREAL8Vector( timeVec );
    XLAL_ERROR( XLAL_ENOMEM );
  }

  /* The time we want to take as the peak time depends on l and m */
  /* Calculate the adjustment we need to make here */
  nrDeltaT = XLALSimIMREOBGetNRPeakDeltaT( l, m, eta );
  if ( XLAL_IS_REAL8_FAIL_NAN( nrDeltaT ) )
  {
    XLALDestroyREAL8Vector( q1LM );
    XLALDestroyREAL8Vector( q2LM );
    XLALDestroyREAL8Vector( q3LM );
    XLALDestroyREAL8Vector( timeVec );
    XLAL_ERROR( XLAL_EFUNC );
  }

  nrTimePeak = timePeak + nrDeltaT;

  /* We are now in a position to use the interp stuff to calculate the derivatives we need */
  /* We will start with the quantities used in the calculation of the a coefficients */
  spline = gsl_spline_alloc( gsl_interp_cspline, amplitude->length );
  acc    = gsl_interp_accel_alloc();

  /* Q1 */
  gsl_spline_init( spline, timeVec->data, q1LM->data, q1LM->length );
  gsl_matrix_set( qMatrix, 0, 0, gsl_spline_eval( spline, nrTimePeak, acc ) );
  gsl_matrix_set( qMatrix, 1, 0, gsl_spline_eval_deriv( spline, nrTimePeak, acc ) );
  gsl_matrix_set( qMatrix, 2, 0, gsl_spline_eval_deriv2( spline, nrTimePeak, acc ) );

  /* Q2 */
  gsl_spline_init( spline, timeVec->data, q2LM->data, q2LM->length );
  gsl_interp_accel_reset( acc );
  gsl_matrix_set( qMatrix, 0, 1, gsl_spline_eval( spline, nrTimePeak, acc ) );
  gsl_matrix_set( qMatrix, 1, 1, gsl_spline_eval_deriv( spline, nrTimePeak, acc ) );
  gsl_matrix_set( qMatrix, 2, 1, gsl_spline_eval_deriv2( spline, nrTimePeak, acc ) );

  /* Q3 */
  gsl_spline_init( spline, timeVec->data, q3LM->data, q3LM->length );
  gsl_interp_accel_reset( acc );
  gsl_matrix_set( qMatrix, 0, 2, gsl_spline_eval( spline, nrTimePeak, acc ) );
  gsl_matrix_set( qMatrix, 1, 2, gsl_spline_eval_deriv( spline, nrTimePeak, acc ) );
  gsl_matrix_set( qMatrix, 2, 2, gsl_spline_eval_deriv2( spline, nrTimePeak, acc ) );

  /* Amplitude */
  gsl_spline_init( spline, timeVec->data, amplitude->data, amplitude->length );
  gsl_interp_accel_reset( acc );
  a     = gsl_spline_eval( spline, nrTimePeak, acc );
  aDot  = gsl_spline_eval_deriv( spline, nrTimePeak, acc );
  aDDot = gsl_spline_eval_deriv2( spline, nrTimePeak, acc );

  nra = GetNRPeakAmplitude( l, m, eta );
  nraDDot = GetNRPeakADDot( l, m, eta );

  if ( XLAL_IS_REAL8_FAIL_NAN( nra ) || XLAL_IS_REAL8_FAIL_NAN( nraDDot ) )
  {
    XLALDestroyREAL8Vector( q1LM );
    XLALDestroyREAL8Vector( q2LM );
    XLALDestroyREAL8Vector( q3LM );
    XLALDestroyREAL8Vector( timeVec );
    XLAL_ERROR( XLAL_EFUNC );
  }

  gsl_vector_set( amps, 0, nra - a );
  gsl_vector_set( amps, 1, - aDot );
  gsl_vector_set( amps, 2, nraDDot - aDDot );

  /* We have now set up all the stuff to calculate the a coefficients */
  /* So let us do it! */
  gsl_linalg_LU_decomp( qMatrix, perm1, &signum );
  gsl_linalg_LU_solve( qMatrix, perm1, amps, aCoeff );

  /* Now we (should) have calculated the a values. Now we can do the b values */

  /* P1 */
  gsl_spline_init( spline, timeVec->data, p1->data, p1->length );
  gsl_interp_accel_reset( acc );
  gsl_matrix_set( pMatrix, 0, 0, - gsl_spline_eval_deriv( spline, nrTimePeak, acc ) );
  gsl_matrix_set( pMatrix, 1, 0, - gsl_spline_eval_deriv2( spline, nrTimePeak, acc ) );

  /* P2 */
  gsl_spline_init( spline, timeVec->data, p2->data, p2->length );
  gsl_interp_accel_reset( acc );
  gsl_matrix_set( pMatrix, 0, 1, - gsl_spline_eval_deriv( spline, nrTimePeak, acc ) );
  gsl_matrix_set( pMatrix, 1, 1, - gsl_spline_eval_deriv2( spline, nrTimePeak, acc ) );

  /* Phase */
  gsl_spline_init( spline, timeVec->data, phase->data, phase->length );
  gsl_interp_accel_reset( acc );
  omega    = gsl_spline_eval_deriv( spline, nrTimePeak, acc );
  omegaDot = gsl_spline_eval_deriv2( spline, nrTimePeak, acc );

  /* Since the phase can be decreasing, we need to take care not to have a -ve frequency */
  if ( omega * omegaDot > 0.0 )
  {
    omega    = fabs( omega );
    omegaDot = fabs( omegaDot );
  }
  else
  {
    omega    = fabs( omega );
    omegaDot = - fabs( omegaDot );
  }

  nromega = GetNRPeakOmega( l, m, eta );
  nromegaDot = GetNRPeakOmegaDot( l, m, eta );

  if ( XLAL_IS_REAL8_FAIL_NAN( nromega ) || XLAL_IS_REAL8_FAIL_NAN( nromegaDot ) )
  {
    XLALDestroyREAL8Vector( q1LM );
    XLALDestroyREAL8Vector( q2LM );
    XLALDestroyREAL8Vector( q3LM );
    XLALDestroyREAL8Vector( timeVec );
    XLAL_ERROR( XLAL_EFUNC );
  }

  gsl_vector_set( omegaVec, 0, nromega - omega );
  gsl_vector_set( omegaVec, 1, nromegaDot - omegaDot );

  /* And now solve for the b coefficients */
  gsl_linalg_LU_decomp( pMatrix, perm2, &signum );
  gsl_linalg_LU_solve( pMatrix, perm2, omegaVec, bCoeff );

  /* We can now populate the coefficients structure */
  coeffs->a1 = gsl_vector_get( aCoeff, 0 );
  coeffs->a2 = gsl_vector_get( aCoeff, 1 );
  coeffs->a3 = gsl_vector_get( aCoeff, 2 );
  coeffs->b1 = gsl_vector_get( bCoeff, 0 );
  coeffs->b2 = gsl_vector_get( bCoeff, 1 );

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
  XLALDestroyREAL8Vector( timeVec );

  return XLAL_SUCCESS;
}

/* ------------------------------------------------
 *          Spin (SEOBNRv1)
 * ------------------------------------------------*/

/**
 * The time difference between the orbital peak and the peak amplitude
 * of the mode in question (currently only 2,2 implemented ).
 * Eq. 33 of Taracchini et al. PRD 86, 024011 (2012).
 */
UNUSED static inline REAL8 XLALSimIMREOBGetNRSpinPeakDeltaT( 
                 INT4 l,           /**<< Mode l */
                 INT4 m,           /**<< Mode m */
                 REAL8 UNUSED eta, /**<< Symmetric mass ratio */
                 REAL8 a           /**<< Dimensionless spin */
                 )
{

  switch ( l )
  {
    case 2:
      switch ( m )
      {
        case 2:
          if ( a <= 0.0 )
          {
            return 2.5;
          }
          else
          {
            return (2.5 + 1.77*a*a*a*a/(0.43655*0.43655*0.43655*0.43655)/(1.0-2.0*eta)/(1.0-2.0*eta)/(1.0-2.0*eta)/(1.0-2.0*eta));
          }
          break;
        default:
          XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    default:
      XLAL_ERROR_REAL8( XLAL_EINVAL );
  }

  /* We should never get here, but I expect a compiler whinge without it... */
  XLALPrintError( "XLAL Error %s - We should never get here!!\n", __func__ );
  XLAL_ERROR_REAL8( XLAL_EINVAL );
}

/* FIXME: Add XLALSimIMREOB to these function names */

/**
 * Peak amplitude predicted by fitting NR results (currently only 2,2 available).
 * Tables IV and V and Eq. 42 of Taracchini et al. PRD 86, 024011 (2012).
 */
UNUSED static inline REAL8 GetNRSpinPeakAmplitude( INT4 UNUSED l, INT4 UNUSED m, REAL8 UNUSED eta, REAL8 UNUSED a )
{
  /* Fit for HOMs missing */
  return 1.3547468629743946*eta + 0.9187885481024214*eta*eta;
}

/**
 * Peak amplitude curvature predicted by fitting NR results (currently only 2,2 available).
 * Tables IV and V and Eq. 42 of Taracchini et al. PRD 86, 024011 (2012).
 */
UNUSED static inline REAL8 GetNRSpinPeakADDot( INT4 UNUSED l, INT4 UNUSED m, REAL8 UNUSED eta, REAL8 UNUSED a )
{
  /* Fit for HOMs missing */
  return eta*(-0.0024971911410897156 + (-0.006128515435641139 + 0.01732656*a/(2.0-4.0*eta))*eta);
}

/**
 * Peak frequency predicted by fitting NR results (currently only 2,2 available).
 * Tables IV and V and Eq. 42 of Taracchini et al. PRD 86, 024011 (2012).
 */
UNUSED static inline REAL8 GetNRSpinPeakOmega( INT4 UNUSED l, INT4 UNUSED m, REAL8 UNUSED eta, REAL8 a )
{
  /* Fit for HOMs missing */
  return 0.27581190323955274 + 0.19347381066059993*eta
       - 0.08898338208573725*log(1.0 - a/(1.0-2.0*eta))
       + eta*eta*(1.78832*(0.2690779744133912 + a/(2.0-4.0*eta))*(1.2056469070395925
       + a/(2.0-4.0*eta)) + 1.423734113371796*log(1.0 - a/(1.0-2.0*eta)));
}

/**
 * Peak frequency slope predicted by fitting NR results (currently only 2,2 available).
 * Tables IV and V and Eq. 42 of Taracchini et al. PRD 86, 024011 (2012).
 */
UNUSED static inline REAL8 GetNRSpinPeakOmegaDot( INT4 UNUSED l, INT4 UNUSED m, REAL8 UNUSED eta, REAL8 UNUSED a )
{
  /* Fit for HOMs missing */
  return 0.006075014646800278 + 0.012040017219351778*eta
       + (0.0007353536801336875 + 0.0015592659912461832*a/(1.0-2.0*eta))*log(1.0-a/(1.0-2.0*eta))
       + eta*eta*(0.03575969677378844 + (-0.011765658882139 - 0.02494825585993893*a/(1.0-2.0*eta))
       * log(1.0 - a/(1.0-2.0*eta)));
}

/**
 * Function to interpolate known amplitude NQC coeffcients of spin terms, 
 * namely a3s, a4 and a5. 
 * The a3s, a4 and a5 values were calculated for 
 * 11 mass ratios q=1,1.5,2,3,4,5,6,10,20,50 and 100, and
 * 19 spin (\f$\chi\f$ defined in Taracchini et al. PRD 86, 024011 (2012)) values
 * chi = -1, -0.9, -0.8, ......, 0.3, 0.4, 0.5, 0.55, 0.6, 0.65.
 * The calculation was done by Andrea Taracchini using a C++ code of the UMaryland group.
 * In principle, these numbers can be automatically calculated iteratively by the LAL code.
 * However, since such calcualtion increase the cost of each waveform generation by
 * about an order of magnitude, we prepare these numbers in advance reduce cost.
 * These number can be verified by confirming that
 * the peak amplitude and frequency agree well with the NR-fits predicted values,
 * and to get exact NR-fits predicted values, corrections on these numbers are ~1%.
 */
UNUSED static int XLALSimIMRGetEOBCalibratedSpinNQC( EOBNonQCCoeffs * restrict coeffs, 
                                    INT4 UNUSED l, 
                                    INT4 UNUSED m, 
                                    REAL8 eta, 
                                    REAL8 a )
{
  const unsigned int qdim = 11;
  const unsigned int adim = 19;
  UINT4 i;
  REAL8 eta2 = eta*eta;
  REAL8 a3slist[qdim], a4list[qdim], a5list[qdim];

  memset( coeffs, 0, sizeof( *coeffs ) );

  const double etalist[11] = {0.009802960494069209, 0.01922337562475971, 0.04535147392290249, 0.08264462809917356, 0.1224489795918367, 0.1388888888888889, 0.16, 0.1875, 0.2222222222222222, 0.24, 0.25};
  const double alist[19]   = {-1., -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.65};
  const double a3stab[11][19] = {
{872.2946684165909, 746.2236084123653, 649.163149941468, 572.6541433830212, 465.8188017225688, 392.3649675720608, 273.0917193587359, 169.5591766990823, 63.68981658441692, -43.88053714217796, 0, -272.4077729342899, -407.3959936134043, -601.9350655151072, -969.4473565113408, -1879.762065510577, -2878.764008622104, -4788.461950454964, -9121.788878332776}, 
{771.54227274491, 698.2568074742545, 623.2537930979673, 523.9160488272244, 451.8201785920314, 353.5866336833333, 262.4023361354035, 168.5588834722693, 75.43076043937653, -19.78515042333346, 0, -222.124184085555, -345.742033226789, -527.5838765272559, -880.1532464805088, -1747.367802632628, -2700.520685884339, -4508.530931365056, -8533.70666018925},
{584.790360938073, 520.3799506107248, 463.3016082992465, 416.8112616076601, 346.4483559089367, 279.6007890923834, 213.9763191701454, 147.2590484055326, 79.39716271746178, 10.54388021850271, 0, -139.3629611914365, -233.9974720900363, -383.2073518275249, -684.7706897288465, -1462.49693792761, -2344.62890234302, -4047.628247528404, -7586.070858003201}, 
{358.393024067137, 315.9910713665202, 283.5448822459095, 250.6617964645584, 215.8716005846251, 178.7758869206099, 141.5921001625318, 102.8844042813548, 62.17683409083224, 19.65380048089165, 0, -75.77079505533445, -138.3818921628671, -237.4789081610928, -441.3903579224361, -1014.797333982442, -1757.366834301527, -3300.837096245112, -6261.106246994135}, 
{172.4502172301744, 158.5094665703111, 144.7737145031524, 129.2974131770823, 111.1148613553597, 94.15668672635928, 77.75309763031652, 58.39881789401274, 37.72991474641616, 15.64789354005334, 0, -37.06966893088249, -71.73757489469652, -122.9636184117563, -223.7175794724878, -537.2335342465955, -1045.546564932042, -2277.1994866329, -4656.964146344593}, 
{138.956425940656, 125.2719770864122, 114.2439547678241, 100.1159745726933, 85.59329075111867, 72.0606997841374, 58.11047183384493, 43.66194754319626, 29.16103741236397, 12.69247383987067, 0, -25.24966987681743, -50.1690186496535, -85.57052809341745, -154.7651285088356, -393.2414438928014, -818.813818093502, -1937.875485658054, -4257.755349339037}, 
{107.3542819894746, 96.45163775625795, 83.05367551805435, 73.66541390670488, 61.00285452993759, 50.47167731568784, 39.00802845343568, 28.85305339143407, 18.63597582551556, 9.482152329860394, 0, -12.29119598908415, -25.51322560066465, -44.77484447801996, -87.41372881073904, -278.1489226068006, -664.8606416592977, -1719.018521321553, -4012.506015355466}, 
{75.29005335033047, 72.32250406953167, 58.53362183893035, 52.06080677002219, 41.60827115985353, 30.80218832505059, 21.77241402931981, 14.26456047427533, 8.234600513300206, 4.565505486274613, 0, 1.922706602892058, 0.3841659186770419, -4.542110317519035, -32.44856494826053, -212.1879218643041, -573.3055734566757, -1490.011721489433, -3363.952598596728}, 
{60.692976823583, 52.85477687110887, 42.82625922078951, 32.95085864555419, 23.43415934606698, 13.5470716494695, 5.265085278666652, -0.8053741207511365, -4.571218113049652, -2.064096170310852, 0, 17.18966429426611, 27.12229185726782, 35.68229478261382, 31.14340126034695, -65.56870048289989, -227.895073238744, -550.8261381951995, -1077.742850153688}, 
{58.31543686820211, 49.08571982189395, 40.67414304823948, 28.80664943836456, 17.5708975225823, 7.138190994477281, -3.011281753628637, -10.76734978377357, -12.38138059751043, -7.157615791667049, 0, 23.04889538706169, 42.99693048409946, 61.39518190726074, 81.95419238411841, 82.37629919425241, 63.58902980799827, 40.88041692069835, 55.12022679785675}, 
{58.68532163446057, 49.55879586941457, 38.96859137125088, 27.34667461300703, 14.90593562710915, 2.574347638902496, -9.100158967414277, -15.15689317221503, -18.23765656057034, -10.43487756797647, 0, 28.64225047483616, 52.85025943662195, 82.87283316068358, 121.6690475515416, 184.8371077453084, 242.7131929263538, 349.5858204855281, 573.964402028233}
};
  const double a4tab[11][19] = {
{-3304.799850690366, -2784.349209575298, -2391.025697134022, -2089.746281865023, -1673.798912163021, -1398.700437123426, -952.9034067988564, -579.6592041964316, -209.2770626411607, 154.4172338379267, 0, 884.5288472809657, 1293.603960974343, 1878.230301261083, 2994.90262286655, 5802.143292527219, 8906.377031852466, 14870.18691371657, 28510.23484335775}, 
{-2889.642454404152, -2589.810076920876, -2288.53773344322, -1895.324851078318, -1619.112769976017, -1246.550298559881, -910.5320847040928, -573.8703373166688, -250.1708225426001, 70.05503090466014, 0, 712.5533871229151, 1085.994414921443, 1631.120158656928, 2700.440647403466, 5367.589187008396, 8325.257440946636, 13971.28687387353, 26666.79702454237},
{-2134.136231872978, -1877.064268185912, -1653.931484260663, -1476.842530996682, -1209.990273735086, -962.7622975860129, -726.1421837252742, -491.8117910776703, -260.3861230995321, -33.29803386526321, 0, 434.165072872907, 714.923885864257, 1156.242658278754, 2058.848992439291, 4425.151920702042, 7141.335186948791, 12439.58451030074, 23595.42126585374}, 
{-1248.948965434483, -1087.019414396883, -966.01256921527, -845.5719725008918, -720.2761986685625, -589.1130064192366, -460.6459037122086, -329.9813959836019, -196.0764833907035, -60.38867455952715, 0, 228.031636471871, 408.2999699714043, 692.0230305275416, 1283.744306953042, 2991.302086945571, 5252.778903506285, 10041.62569460853, 19433.33867280159}, 
{-558.3319803695712, -509.4198489809961, -461.9563390321997, -409.0499446611043, -347.4458047108782, -291.342148681483, -238.3136913325, -176.48330909104, -111.931229130744, -44.7807649896551, 0, 108.0802789243228, 203.6925650072404, 343.2462397232334, 623.4631540005192, 1545.167396199536, 3105.358405620486, 6995.247478862227, 14740.47388540465}, 
{-443.0800026835253, -395.5266756529914, -358.1132299746496, -310.3522657897121, -262.0471407656797, -218.0225941855192, -173.4923690434391, -128.3241666756531, -84.18764217333766, -35.11157686619008, 0, 72.42303122610437, 139.5409415608551, 233.6859598886507, 424.5006079645129, 1137.141268822077, 2470.706048016958, 6064.305250831422, 13711.56479304423}, 
{-339.0235915844027, -301.6768991985205, -256.1954692167641, -225.122747136184, -183.4493757774726, -149.6273380222367, -113.4045680331683, -82.23276412066362, -51.69646716347211, -25.50065274962909, 0, 33.78549696898528, 67.7793532347804, 118.128069496536, 240.6452125471342, 844.8808225671395, 2100.766006843633, 5567.719753915647, 13238.13094099677}, 
{-239.7783765681118, -229.6219954929181, -183.0691615728538, -161.4416215851445, -127.0161762815546, -92.06103046507863, -63.46566969908776, -40.26151156959535, -22.2346091726419, -12.07050074675413, 0, -7.691352010334845, -4.67175092301215, 11.33495560642181, 109.0140231343323, 728.6782625485398, 1951.694339911332, 5051.985500326179, 11446.93133618593}, 
{-207.4649427125497, -180.8819079763306, -146.3668476297406, -113.0010560672213, -80.82081391687568, -48.25429242084964, -21.12530533597964, -1.253140183557529, 11.25171807978124, 4.762279251524356, 0, -49.98063474223146, -74.12605255860477, -85.98932491689175, -37.15343028084588, 355.9530561615545, 963.8327713818729, 2150.203870925892, 4105.560245358544}, 
{-209.3963309944106, -177.4126692607984, -148.650816945404, -108.0577395641345, -69.90482433518206, -34.84500498095979, -1.390551419992443, 24.28317025722451, 31.60107555201323, 18.60037455191576, 0, -64.25912340317869, -114.908145772787, -151.1508355309806, -173.7858106326595, -87.5386870155428, 48.65552876803645, 232.6203762611659, 340.4413234234147}, 
{-216.1423631616688, -184.6508002958967, -148.1966726639886, -108.5725305431534, -66.32430505800025, -25.04165933587504, 13.59473105492496, 34.99230429165542, 47.3941356433856, 27.85684411598417, 0, -78.9492732782613, -140.4593668007498, -210.1824360384644, -288.3810317382483, -406.891227028902, -528.0421625843486, -789.1807718418322, -1415.686288121717}
};
  const double a5tab[11][19] = {
{3150.016247770505, 2611.618198196962, 2211.665504985066, 1913.527882869167, 1506.890259721349, 1247.653350744291, 829.0244538318248, 490.7965411134306, 165.1931841920994, -143.5225807274245, 0, -727.7143449867366, -1036.956889037733, -1475.890989150695, -2325.87892786623, -4496.443259384486, -6914.421651064408, -11581.57297036935, -22336.50959693963}, 
{2723.551650891905, 2415.696720726347, 2111.815454623907, 1721.173632012413, 1455.105394608006, 1100.081866949763, 788.776259155259, 485.1307952677881, 202.3107419090617, -68.21906751310863, 0, -579.4715730486229, -861.19421867575, -1269.829837439378, -2082.678708982266, -4139.764970387218, -6440.727137007285, -10859.60711156812, -20891.14198546364},
{1961.803406006329, 1704.39219882756, 1485.324763249102, 1315.636640578445, 1061.361853965819, 831.4871101463038, 616.8904896253562, 409.8390456628748, 211.3169707595789, 23.01741246018258, 0, -342.7214007224965, -550.9609674162657, -877.7557557023476, -1555.250692260378, -3361.292357694617, -5458.466479276242, -9590.160446008238, -18402.26301790666}, 
{1097.917491712732, 942.7099212993509, 829.2629492433241, 718.3193294518603, 604.8041059999124, 488.1234260327431, 376.412681986617, 265.3811233156874, 154.5158349034548, 45.58889720668928, 0, -173.2728677830992, -303.0696121701943, -506.3608767132163, -936.8574238296534, -2212.691135008344, -3939.913483836887, -7663.621488890864, -15126.13271660938}, 
{455.965892024021, 412.8249887882179, 371.5733488147668, 326.0890985574885, 273.6128958291934, 226.923561951563, 183.7837468152386, 134.0832361070535, 83.36595498944038, 32.02352547177204, 0, -79.11769885068468, -144.8391544665689, -239.5251978170959, -434.4296549859647, -1114.68968302118, -2315.479174848672, -5391.12169966747, -11691.50396157746}, 
{356.3674499528969, 314.8655864112112, 282.9505475735912, 242.3801533575559, 202.0030591498134, 165.995454615505, 130.2563398123429, 94.75627128562655, 60.99029652944273, 24.26183377077411, 0, -51.99920816655854, -96.85454037716565, -158.9042161566986, -290.4078886525589, -825.2312116888984, -1872.041288836954, -4757.018234966477, -11048.42540795819}, 
{269.8489760264307, 237.7356164009191, 199.0022600089494, 173.1643050093598, 138.7416802973238, 111.4580791174302, 82.72004581103909, 58.69499923628382, 35.80186051569432, 17.04900950683294, 0, -23.02315959834502, -44.39224047056896, -76.80780846002243, -165.1100914349566, -644.1021872517031, -1661.762820379876, -4504.413410486936, -10895.25396530118}, 
{191.6269943298458, 182.8905043917497, 143.5657624122904, 125.4573733121321, 97.07279993133585, 68.77121273309662, 46.10443492472263, 28.1639705641721, 14.7127004821589, 7.778381914853679, 0, 7.356961680608462, 6.437401700524503, -6.283584763149717, -89.77008977040909, -616.1545845241945, -1642.665112798123, -4246.527019941841, -9669.468308473082}, 
{174.3283075616382, 152.0021950132451, 122.5167243363872, 94.53237349877479, 67.54064550238134, 40.91702246550422, 18.90695056439407, 2.884324487148479, -7.196918180301969, -2.846247367234787, 0, 36.70406314006223, 51.11063423268788, 50.31837535254078, -12.10728963257653, -385.5303262789694, -936.5327847218327, -2001.10679816206, -3776.641036193068}, 
{181.3374638359205, 154.0901974860778, 129.9340904029422, 95.65646684513354, 63.70048086504236, 34.66980400992438, 7.507589203658138, -13.29466480375627, -20.32109843184418, -12.1907950229181, 0, 45.14242771651794, 76.88295117589163, 90.58383588123166, 78.82154764055255, -52.77855613211774, -217.9374491494236, -450.0841472460368, -658.2583906494992}, 
{189.7236933548747, 163.1987007799873, 132.4489445305476, 99.26806905879204, 63.98267137712187, 29.98763963100855, -1.435655737202815, -19.58160817704562, -30.85589303305039, -18.66717165997227, 0, 54.71289909453592, 93.17391460212183, 130.8975241448977, 161.3660740779782, 196.0770241706157, 245.5768655062068, 392.8648029097516, 820.948713060581}
};
  
  /* Stuff for interpolating the data */
  gsl_spline    *spline = NULL;
  gsl_interp_accel *acc = NULL;

  spline = gsl_spline_alloc( gsl_interp_cspline, adim );
  acc    = gsl_interp_accel_alloc();
  for (i = 0; i < qdim; i++)
  {
    gsl_spline_init( spline, alist, a3stab[i], adim );
    gsl_interp_accel_reset( acc );
    a3slist[i] = gsl_spline_eval( spline, a/(1.0-2.0*eta), acc );
    gsl_spline_init( spline, alist, a4tab[i], adim );
    gsl_interp_accel_reset( acc );
    a4list[i] = gsl_spline_eval( spline, a/(1.0-2.0*eta), acc );
    gsl_spline_init( spline, alist, a5tab[i], adim );
    gsl_interp_accel_reset( acc );
    a5list[i] = gsl_spline_eval( spline, a/(1.0-2.0*eta), acc );
  }
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
 
  coeffs->a1 = -12.67955358602124 + 75.41927959573084 * eta - 106.15933052937714 * eta2;
  coeffs->a2 = 101.45522216901628 - 757.3158549733314 * eta + 1473.314771676588 * eta2;
  coeffs->a3 = -107.6647834845902 + 857.6219519536213 * eta - 1776.2776804623143 * eta2;

  spline = gsl_spline_alloc( gsl_interp_cspline, qdim );
  acc    = gsl_interp_accel_alloc();
  gsl_spline_init( spline, etalist, a3slist, qdim );
  gsl_interp_accel_reset( acc );
  coeffs->a3S = gsl_spline_eval( spline, eta, acc );
  gsl_spline_init( spline, etalist, a4list, qdim );
  gsl_interp_accel_reset( acc );
  coeffs->a4 = gsl_spline_eval( spline, eta, acc );
  gsl_spline_init( spline, etalist, a5list, qdim );
  gsl_interp_accel_reset( acc );
  coeffs->a5 = gsl_spline_eval( spline, eta, acc );
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);

  /* Andrea and I have different sign conventions, so I need to put a minus sign in front */
  coeffs->b1 = - (-1.464129495621165 + 12.81732978488213 * eta - 60.09957767247623 * eta2);
  coeffs->b2 = - ( 7.477426352542122 - 85.26122117590637 * eta + 353.3251639728075 * eta2);

  return XLAL_SUCCESS;

}

/**
 * This function computes the coefficients a3s, a4, etc. used in the
 * non-quasicircular correction. The details of the calculation of these
 * coefficients are found in the DCC document T1100433. 
 * In brief, this function populates and solves the linear equations 
 * Eq. 18 (for amplitude) and Eq. 19 (for phase) of the DCC document T1100433v2.
 */
UNUSED static int XLALSimIMRSpinEOBCalculateNQCCoefficients(
                 REAL8Vector    * restrict amplitude,   /**<< Waveform amplitude, func of time */
                 REAL8Vector    * restrict phase,       /**<< Waveform phase(rad), func of time */
                 REAL8Vector    * restrict rVec,        /**<< Position-vector, function of time */
                 REAL8Vector    * restrict prVec,       /**<< Momentum vector, function of time */
                 REAL8Vector    * restrict orbOmegaVec, /**<< Orbital frequency, func of time */
                 INT4                      l,           /**<< Mode index l */
                 INT4                      m,           /**<< Mode index m */
                 REAL8                     timePeak,    /**<< Time of peak orbital frequency */
                 REAL8                     deltaT,      /**<< Sampling interval */
                 REAL8                     eta,         /**<< Symmetric mass ratio */
                 REAL8                     a,           /**<< Normalized spin of deformed-Kerr */
                 EOBNonQCCoeffs * restrict coeffs       /**<< OUTPUT, NQC coefficients */)
{

  /* For gsl permutation stuff */

  int signum;

  REAL8Vector * restrict timeVec = NULL;

  /* Vectors which are used in the computation of the NQC coefficients */
  REAL8Vector *q3 = NULL, *q4 = NULL, *q5 = NULL;
  REAL8Vector *p3 = NULL, *p4 = NULL;

  REAL8Vector *qNS = NULL, *pNS = NULL;

  /* Since the vectors we actually want are q etc * A, we will have to generate them here */
  REAL8Vector *q3LM  = NULL;
  REAL8Vector *q4LM  = NULL;
  REAL8Vector *q5LM  = NULL; 
  REAL8Vector *qNSLM = NULL;

  REAL8 amp, aDot, aDDot;
  REAL8 omega, omegaDot;

  REAL8 qNSLMPeak, qNSLMDot, qNSLMDDot;
  REAL8 pNSLMDot, pNSLMDDot;

  REAL8 nra, nraDDot;
  REAL8 nromega, nromegaDot;

  REAL8 nrDeltaT, nrTimePeak;

  /* Stuff for finding numerical derivatives */
  gsl_spline    *spline = NULL;
  gsl_interp_accel *acc = NULL;

  /* Matrix stuff for calculating coefficients */
  gsl_matrix *qMatrix = NULL, *pMatrix = NULL;
  gsl_vector *aCoeff  = NULL, *bCoeff  = NULL;

  gsl_vector *amps = NULL, *omegaVec = NULL;

  gsl_permutation *perm1 = NULL, *perm2 = NULL;

  memset( coeffs, 0, sizeof( EOBNonQCCoeffs ) );

  /* Populate the time vector */
  /* It is okay to assume initial t = 0 */
  timeVec = XLALCreateREAL8Vector( rVec->length );
  q3    = XLALCreateREAL8Vector( rVec->length );
  q4    = XLALCreateREAL8Vector( rVec->length );
  q5    = XLALCreateREAL8Vector( rVec->length );
  p3    = XLALCreateREAL8Vector( rVec->length );
  p4    = XLALCreateREAL8Vector( rVec->length );
  qNS   = XLALCreateREAL8Vector( rVec->length );
  pNS   = XLALCreateREAL8Vector( rVec->length );
  q3LM  = XLALCreateREAL8Vector( rVec->length );
  q4LM  = XLALCreateREAL8Vector( rVec->length );
  q5LM  = XLALCreateREAL8Vector( rVec->length );
  qNSLM = XLALCreateREAL8Vector( rVec->length );

  if ( !timeVec || !q3 || !q4 || !q5 || !p3 || !p4 || !qNS || !pNS || !q3LM
          || !q4LM || !q5LM || !qNSLM )
  {
    XLALDestroyREAL8Vector( timeVec );
    XLALDestroyREAL8Vector( q3 );
    XLALDestroyREAL8Vector( q4 );
    XLALDestroyREAL8Vector( q5 );
    XLALDestroyREAL8Vector( p3 );
    XLALDestroyREAL8Vector( p4 );
    XLALDestroyREAL8Vector( qNS );
    XLALDestroyREAL8Vector( pNS );
    XLALDestroyREAL8Vector( q3LM );
    XLALDestroyREAL8Vector( q4LM );
    XLALDestroyREAL8Vector( q5LM );
    XLALDestroyREAL8Vector( qNSLM );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* We need the calibrated non-spinning NQC coefficients */
  if ( XLALSimIMRGetEOBCalibratedSpinNQC( coeffs, l, m, eta, a ) == XLAL_FAILURE )
  {
    XLALDestroyREAL8Vector( timeVec );
    XLALDestroyREAL8Vector( q3 );
    XLALDestroyREAL8Vector( q4 );
    XLALDestroyREAL8Vector( q5 );
    XLALDestroyREAL8Vector( p3 );
    XLALDestroyREAL8Vector( p4 );
    XLALDestroyREAL8Vector( qNS );
    XLALDestroyREAL8Vector( pNS );
    XLALDestroyREAL8Vector( q3LM );
    XLALDestroyREAL8Vector( q4LM );
    XLALDestroyREAL8Vector( q5LM );
    XLALDestroyREAL8Vector( qNSLM );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Populate vectors as necessary. Eqs. 14 - 17 of the LIGO DCC document T1100433v2 */
  for ( unsigned int i = 0; i < timeVec->length; i++ )
  {
    
    REAL8 rootR  = sqrt(rVec->data[i]);
    REAL8 rOmega = rVec->data[i] * orbOmegaVec->data[i];

    /* We don't need these as vectors as their coefficients are calibrated */
    REAL8 q1, q2, p1, p2;

    timeVec->data[i] = i * deltaT;
    q1            = prVec->data[i]*prVec->data[i] / (rOmega*rOmega);
    q2            = q1 / rVec->data[i];
    q3->data[i]   = q2 / rootR;
    q4->data[i]   = q2 / rVec->data[i];
    q5->data[i]   = q3->data[i] / rVec->data[i];

    p1          = prVec->data[i] / rOmega;
    p2          = p1 * prVec->data[i] * prVec->data[i];
    p3->data[i] = p2 / rootR;
    p4->data[i] = p2 / rVec->data[i];

    qNS->data[i]  = coeffs->a1 * q1 + coeffs->a2 * q2 + coeffs->a3 * q3->data[i];
    pNS->data[i]  = coeffs->b1 * p1 + coeffs->b2 * p2;
    q3LM->data[i] = q3->data[i] * amplitude->data[i];
    q4LM->data[i] = q4->data[i] * amplitude->data[i];
    q5LM->data[i] = q5->data[i] * amplitude->data[i];

    qNSLM->data[i] = qNS->data[i] * amplitude->data[i];
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
    XLALDestroyREAL8Vector( timeVec );
    XLALDestroyREAL8Vector( q3 );
    XLALDestroyREAL8Vector( q4 );
    XLALDestroyREAL8Vector( q5 );
    XLALDestroyREAL8Vector( p3 );
    XLALDestroyREAL8Vector( p4 );
    XLALDestroyREAL8Vector( qNS );
    XLALDestroyREAL8Vector( pNS );
    XLALDestroyREAL8Vector( q3LM );
    XLALDestroyREAL8Vector( q4LM );
    XLALDestroyREAL8Vector( q5LM );
    XLALDestroyREAL8Vector( qNSLM );
    XLAL_ERROR( XLAL_ENOMEM );
  }

  /* The time we want to take as the peak time depends on l and m */
  /* Calculate the adjustment we need to make here */
  nrDeltaT   = XLALSimIMREOBGetNRSpinPeakDeltaT( l, m, eta, a );
  if ( XLAL_IS_REAL8_FAIL_NAN( nrDeltaT ) )
  {
    XLALDestroyREAL8Vector( timeVec );
    XLALDestroyREAL8Vector( q3 );
    XLALDestroyREAL8Vector( q4 );
    XLALDestroyREAL8Vector( q5 );
    XLALDestroyREAL8Vector( p3 );
    XLALDestroyREAL8Vector( p4 );
    XLALDestroyREAL8Vector( qNS );
    XLALDestroyREAL8Vector( pNS );
    XLALDestroyREAL8Vector( q3LM );
    XLALDestroyREAL8Vector( q4LM );
    XLALDestroyREAL8Vector( q5LM );
    XLALDestroyREAL8Vector( qNSLM );
    XLAL_ERROR( XLAL_EFUNC );
  }

  nrTimePeak = timePeak - nrDeltaT;

  /* We are now in a position to use the interp stuff to calculate the derivatives we need */
  /* We will start with the quantities used in the calculation of the a coefficients */
  spline = gsl_spline_alloc( gsl_interp_cspline, amplitude->length );
  acc    = gsl_interp_accel_alloc();

  /* Populate the Q matrix in Eq. 18 of the LIGO DCC document T1100433v2 */
  /* Q3 */
  gsl_spline_init( spline, timeVec->data, q3LM->data, q3LM->length );
  gsl_matrix_set( qMatrix, 0, 0, gsl_spline_eval( spline, nrTimePeak, acc ) );
  gsl_matrix_set( qMatrix, 1, 0, gsl_spline_eval_deriv( spline, nrTimePeak, acc ) );
  gsl_matrix_set( qMatrix, 2, 0, gsl_spline_eval_deriv2( spline, nrTimePeak, acc ) );

  /* Q4 */
  gsl_spline_init( spline, timeVec->data, q4LM->data, q4LM->length );
  gsl_interp_accel_reset( acc );
  gsl_matrix_set( qMatrix, 0, 1, gsl_spline_eval( spline, nrTimePeak, acc ) );
  gsl_matrix_set( qMatrix, 1, 1, gsl_spline_eval_deriv( spline, nrTimePeak, acc ) );
  gsl_matrix_set( qMatrix, 2, 1, gsl_spline_eval_deriv2( spline, nrTimePeak, acc ) );

  /* Q5 */
  gsl_spline_init( spline, timeVec->data, q5LM->data, q5LM->length );
  gsl_interp_accel_reset( acc );
  gsl_matrix_set( qMatrix, 0, 2, gsl_spline_eval( spline, nrTimePeak, acc ) );
  gsl_matrix_set( qMatrix, 1, 2, gsl_spline_eval_deriv( spline, nrTimePeak, acc ) );
  gsl_matrix_set( qMatrix, 2, 2, gsl_spline_eval_deriv2( spline, nrTimePeak, acc ) );

  /* Populate the r.h.s vector of Eq. 18 of the LIGO DCC document T1100433v2 */
  /* Amplitude */
  gsl_spline_init( spline, timeVec->data, amplitude->data, amplitude->length );
  gsl_interp_accel_reset( acc );
  amp   = gsl_spline_eval( spline, nrTimePeak, acc );
  aDot  = gsl_spline_eval_deriv( spline, nrTimePeak, acc );
  aDDot = gsl_spline_eval_deriv2( spline, nrTimePeak, acc );

  /* qNSLM */
  gsl_spline_init( spline, timeVec->data, qNSLM->data, qNSLM->length );
  gsl_interp_accel_reset( acc );
  qNSLMPeak = gsl_spline_eval( spline, nrTimePeak, acc );
  qNSLMDot  = gsl_spline_eval_deriv( spline, nrTimePeak, acc );
  qNSLMDDot = gsl_spline_eval_deriv2( spline, nrTimePeak, acc );

  nra = GetNRSpinPeakAmplitude( l, m, eta, a );
  nraDDot = - GetNRSpinPeakADDot( l, m, eta, a );

  if ( XLAL_IS_REAL8_FAIL_NAN( nra ) || XLAL_IS_REAL8_FAIL_NAN( nraDDot ) )
  {
    XLALDestroyREAL8Vector( timeVec );
    XLALDestroyREAL8Vector( q3 );
    XLALDestroyREAL8Vector( q4 );
    XLALDestroyREAL8Vector( q5 );
    XLALDestroyREAL8Vector( p3 );
    XLALDestroyREAL8Vector( p4 );
    XLALDestroyREAL8Vector( qNS );
    XLALDestroyREAL8Vector( pNS );
    XLALDestroyREAL8Vector( q3LM );
    XLALDestroyREAL8Vector( q4LM );
    XLALDestroyREAL8Vector( q5LM );
    XLALDestroyREAL8Vector( qNSLM );
    XLAL_ERROR( XLAL_EFUNC );
  }

  gsl_vector_set( amps, 0, nra - amp - qNSLMPeak );
  gsl_vector_set( amps, 1, - aDot - qNSLMDot );
  gsl_vector_set( amps, 2, nraDDot - aDDot - qNSLMDDot );

  /* We have now set up all the stuff to calculate the a coefficients */
  /* So let us do it! */
  gsl_linalg_LU_decomp( qMatrix, perm1, &signum );
  gsl_linalg_LU_solve( qMatrix, perm1, amps, aCoeff );

  /* Now we (should) have calculated the a values. Now we can do the b values */

  /* Populate the P matrix in Eq. 18 of the LIGO DCC document T1100433v2 */
  /* P3 */
  gsl_spline_init( spline, timeVec->data, p3->data, p3->length );
  gsl_interp_accel_reset( acc );
  gsl_matrix_set( pMatrix, 0, 0, - gsl_spline_eval_deriv( spline, nrTimePeak, acc ) );
  gsl_matrix_set( pMatrix, 1, 0, - gsl_spline_eval_deriv2( spline, nrTimePeak, acc ) );

  /* P4 */
  gsl_spline_init( spline, timeVec->data, p4->data, p4->length );
  gsl_interp_accel_reset( acc );
  gsl_matrix_set( pMatrix, 0, 1, - gsl_spline_eval_deriv( spline, nrTimePeak, acc ) );
  gsl_matrix_set( pMatrix, 1, 1, - gsl_spline_eval_deriv2( spline, nrTimePeak, acc ) );

  /* Populate the r.h.s vector of Eq. 18 of the LIGO DCC document T1100433v2 */
  /* Phase */
  gsl_spline_init( spline, timeVec->data, phase->data, phase->length );
  gsl_interp_accel_reset( acc );
  omega    = gsl_spline_eval_deriv( spline, nrTimePeak, acc );
  omegaDot = gsl_spline_eval_deriv2( spline, nrTimePeak, acc );

  /* pNSLM */
  gsl_spline_init( spline, timeVec->data, pNS->data, pNS->length );
  gsl_interp_accel_reset( acc );
  pNSLMDot  = gsl_spline_eval_deriv( spline, nrTimePeak, acc );
  pNSLMDDot = gsl_spline_eval_deriv2( spline, nrTimePeak, acc );

  /* Since the phase can be decreasing, we need to take care not to have a -ve frequency */
  if ( omega * omegaDot > 0.0 )
  {
    omega    = fabs( omega );
    omegaDot = fabs( omegaDot );
  }
  else
  {
    omega    = fabs( omega );
    omegaDot = - fabs( omegaDot );
  }

  //nromega = GetNRPeakOmega( l, m, eta );
  //nromegaDot = GetNRPeakOmegaDot( l, m, eta );
  nromega = GetNRSpinPeakOmega( l, m, eta, a );
  nromegaDot = GetNRSpinPeakOmegaDot( l, m, eta, a );

  //printf("NR inputs: %.16e, %.16e, %.16e, %.16e\n",nra,nraDDot,nromega,nromegaDot);

  if ( XLAL_IS_REAL8_FAIL_NAN( nromega ) || XLAL_IS_REAL8_FAIL_NAN( nromegaDot ) )
  {
    XLALDestroyREAL8Vector( timeVec );
    XLALDestroyREAL8Vector( q3 );
    XLALDestroyREAL8Vector( q4 );
    XLALDestroyREAL8Vector( q5 );
    XLALDestroyREAL8Vector( p3 );
    XLALDestroyREAL8Vector( p4 );
    XLALDestroyREAL8Vector( qNS );
    XLALDestroyREAL8Vector( pNS );
    XLALDestroyREAL8Vector( q3LM );
    XLALDestroyREAL8Vector( q4LM );
    XLALDestroyREAL8Vector( q5LM );
    XLALDestroyREAL8Vector( qNSLM );
    XLAL_ERROR( XLAL_EFUNC );
  }

  gsl_vector_set( omegaVec, 0, nromega - omega + pNSLMDot );
  gsl_vector_set( omegaVec, 1, nromegaDot - omegaDot + pNSLMDDot );

  /*printf( "P MATRIX\n" );
  for (unsigned int i = 0; i < 2; i++ )
  {
    for (unsigned int j = 0; j < 2; j++ )
    {
      printf( "%.12e\t", gsl_matrix_get( pMatrix, i, j ));
    }
    printf( "= %.12e\n", gsl_vector_get( omegaVec, i ) );
  }*/

  /* And now solve for the b coefficients */
  gsl_linalg_LU_decomp( pMatrix, perm2, &signum );
  gsl_linalg_LU_solve( pMatrix, perm2, omegaVec, bCoeff );

  /* We can now populate the coefficients structure */
/*  coeffs->a3S = gsl_vector_get( aCoeff, 0 );
  coeffs->a4  = gsl_vector_get( aCoeff, 1 );
  coeffs->a5  = gsl_vector_get( aCoeff, 2 );*/
  coeffs->b3  = gsl_vector_get( bCoeff, 0 );
  coeffs->b4  = gsl_vector_get( bCoeff, 1 );

  /*printf( "NQC coefficients:\n" );
  printf( "a1 = %.16e, a2 = %.16e, a3 = %.16e, a3s = %.16e, a4 = %.16e, a5 = %.16e\n",
    coeffs->a1, coeffs->a2, coeffs->a3, coeffs->a3S, coeffs->a4, coeffs->a5 );

  printf( "b1 = %.16e, b2 = %.16e, b3 = %.16e, b4 = %.16e\n",
    coeffs->b1, coeffs->b2, coeffs->b3, coeffs->b4 );*/

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

  XLALDestroyREAL8Vector( timeVec );
  XLALDestroyREAL8Vector( q3 );
  XLALDestroyREAL8Vector( q4 );
  XLALDestroyREAL8Vector( q5 );
  XLALDestroyREAL8Vector( p3 );
  XLALDestroyREAL8Vector( p4 );
  XLALDestroyREAL8Vector( qNS );
  XLALDestroyREAL8Vector( pNS );
  XLALDestroyREAL8Vector( q3LM );
  XLALDestroyREAL8Vector( q4LM );
  XLALDestroyREAL8Vector( q5LM );
  XLALDestroyREAL8Vector( qNSLM );

  return XLAL_SUCCESS;
}

#endif /*_LALSIMIMRNQCCORRECTION_C*/
