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
          /* DeltaT22 defined here is a minus sign different from Eq. (33) of Taracchini et al. */
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
  const unsigned int nsqdim = 100;
  const unsigned int   qdim = 50;
  const unsigned int   adim = 69;
  UINT4 i;
  /* REAL8 eta2 = eta*eta;*/
  REAL8 a3slist[qdim], a4list[qdim], a5list[qdim];

  memset( coeffs, 0, sizeof( *coeffs ) );
  const double nsetalist[100] = {0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175, 0.02, 0.0225, 0.025, 0.0275, 0.03, 0.0325, 0.035, 0.0375, 0.04, 0.0425, 0.045, 0.0475, 0.05, 0.0525, 0.055, 0.0575, 0.06, 0.0625, 0.065, 0.0675, 0.07, 0.0725, 0.075, 0.0775, 0.08, 0.0825, 0.085, 0.0875, 0.09, 0.0925, 0.095, 0.0975, 0.1, 0.1025, 0.105, 0.1075, 0.11, 0.1125, 0.115, 0.1175, 0.12, 0.1225, 0.125, 0.1275, 0.13, 0.1325, 0.135, 0.1375, 0.14, 0.1425, 0.145, 0.1475, 0.15, 0.1525, 0.155, 0.1575, 0.16, 0.1625, 0.165, 0.1675, 0.17, 0.1725, 0.175, 0.1775, 0.18, 0.1825, 0.185, 0.1875, 0.19, 0.1925, 0.195, 0.1975, 0.2, 0.2025, 0.205, 0.2075, 0.21, 0.2125, 0.215, 0.2175, 0.22, 0.2225, 0.225, 0.2275, 0.23, 0.2325, 0.235, 0.2375, 0.24, 0.2425, 0.245, 0.2475, 0.25};
  const double etalist[50] = {0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.105, 0.11, 0.115, 0.12, 0.125, 0.13, 0.135, 0.14, 0.145, 0.15, 0.155, 0.16, 0.165, 0.17, 0.175, 0.18, 0.185, 0.19, 0.195, 0.2, 0.205, 0.21, 0.215, 0.22, 0.225, 0.23, 0.235, 0.24, 0.245, 0.25};
  const double alist[69]   = {-1., -0.975, -0.95, -0.925, -0.9, -0.875, -0.85, -0.825, -0.8, -0.775, -0.75, -0.725, -0.7, -0.675, -0.65, -0.625, -0.6, -0.575, -0.55, -0.525, -0.5, -0.475, -0.45, -0.425, -0.4, -0.375, -0.35, -0.325, -0.3, -0.275, -0.25, -0.225, -0.2, -0.175, -0.15, -0.125, -0.1, -0.075, -0.05, -0.025, 0., 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7};

  const double a1list[100] = {-21.648106754277453, -20.47585498319721, -19.530321799953853, -18.760004570030972, -18.079213572570588, -17.419521042795658, -16.859312290615286, -16.28756533894047, -15.772560753656803, -15.303367029222859, -14.847798623544204, -14.404050939670224, -13.976897526988164, -13.578815297784823, -13.200763216702892, -12.841976980076495, -12.501692284239752, -12.179144825526786, -11.865421415437213, -11.561721630442458, -11.267369305715642, -10.980678895252732, -10.699964853049693, -10.42829543384451, -10.167524207712544, -9.914674799750443, -9.669372681081281, -9.43124332282813, -9.199912196114067, -8.97500477206216, -8.756146521795488, -8.542962916437121, -8.335079427110134, -8.130867290476582, -7.9313984345769715, -7.736518055359128, -7.545988962963793, -7.35957396753171, -7.177035879203623, -6.998068318971551, -6.8224288840577945, -6.649961333309703, -6.480441962265538, -6.313647066463564, -6.149352941442041, -5.9873358827392344, -5.827372185893404, -5.669238146442815, -5.512513869477759, -5.347444298469859, -5.183197664196798, -5.019955297249104, -4.857898528217308, -4.697208687691938, -4.538067106263523, -4.380655114522592, -4.225154043059675, -4.071745222465302, -3.9206099833299994, -3.771929656244297, -3.6258855717987246, -3.4826590605838113, -3.342431453190086, -3.2053840802080775, -3.0798495456033685, -2.957143501994997, -2.837141996613968, -2.7197210766912865, -2.6047567894579564, -2.4921251821449832, -2.381702301983373, -2.273364196204128, -2.1669869120382543, -2.0624464967167566, -1.9596189974706397, -1.8495581080328225, -1.7409749049276515, -1.6340876818580403, -1.5291147325269023, -1.4262743506371505, -1.3257848298916983, -1.2278644639934593, -1.1327315466453465, -1.0406043715502744, -0.9517012324111543, -0.8662404229309005, -0.7844402368124264, -0.7065189677586451, -0.6326949094724702, -0.5631863556568149, -0.4982116000145925, -0.4379889362487165, -0.3827366580621, -0.3326730591576565, -0.28801643323829923, -0.24898507400694198, -0.21579727516649713, -0.18867133041987855, -0.1678255334699995, -0.15347817801977337};

  const double a2list[100] = {178.73204288078207, 167.7345170263427, 158.85457776976878, 151.63702661020528, 145.25886103554546, 139.08071349079492, 133.85300186061994, 128.49833153024582, 123.69581397508206, 119.34304449990263, 115.1225036726121, 111.02138715464824, 107.08437843884623, 103.42908442981195, 99.96941370234742, 96.69762600347904, 93.6059810802332, 90.68673867963636, 87.86042720725196, 85.13659365364732, 82.50827361644501, 79.95961751850369, 77.47477578268207, 75.08044751075573, 72.79299640785891, 70.58542621544011, 68.45392194344085, 66.39466860180258, 64.40385120046682, 62.47765474937501, 60.61226425846867, 58.80386473768926, 57.04864119697828, 55.32723154023535, 53.65375061411054, 52.0271070701423, 50.44518833065014, 48.905881817953585, 47.407074954372156, 45.94482445612376, 44.51687981462722, 43.123270309615144, 41.76224017946817, 40.43203366256695, 39.13089499729211, 37.85706842202429, 36.60879817514412, 35.384328495032264, 34.180608559988464, 32.932965816541305, 31.701640999519693, 30.48755199799916, 29.291616701055233, 28.11475299776345, 26.95787877719934, 25.821911928438432, 24.70777034055626, 23.616371902628366, 22.54863450373026, 21.50547603293748, 20.487814379325563, 19.496567431970035, 18.53265307994643, 17.596989212330282, 16.733055464585682, 15.895477726977264, 15.083579856295431, 14.296685709330593, 13.53411914287315, 12.795204013713512, 12.079264178642092, 11.385623494449275, 10.71360581792548, 10.06253500586111, 9.431734915046569, 8.755551062926262, 8.098379005319565, 7.462063310457031, 6.848448546569211, 6.259379281886656, 5.69670008463992, 5.162255523059555, 4.657890165376111, 4.1854485798201475, 3.7467753346222055, 3.3437149980128416, 2.978112138222609, 2.651811323482059, 2.366657122021744, 2.124494102072216, 1.9271668318640272, 1.77651987962773, 1.674397813593876, 1.6226452019930173, 1.6231066130557064, 1.6776266150124943, 1.7880497760939342, 1.956220664530578, 2.1839838485529777, 2.4731838963916855};

  const double a3list[100] = {-198.6684740510964, -185.71560983427335, -175.26102024642407, -166.7798147654891, -159.28986302859136, -152.03810090459336, -145.91734205946045, -139.6355596843493, -134.01653939524905, -128.93945143613453, -124.02034829337536, -119.24645127356649, -114.67037356959001, -110.43118685036107, -106.4264567068743, -102.6466423946065, -99.08220316903453, -95.72359828563523, -92.4798228915885, -89.36096908139594, -86.35837723133453, -83.45329697039512, -80.62697792756852, -77.90968987612793, -75.32016588735331, -72.82709457524146, -70.42579806564999, -68.11159848443651, -65.87981795745866, -63.72577861057402, -61.64480256964021, -59.63221196051487, -57.683328909055604, -55.7731074201788, -53.92035933426217, -52.12397499502319, -50.38150684315194, -48.69050731933857, -47.04852886427319, -45.45070004394427, -43.89439067703587, -42.37998904664536, -40.90552002179141, -39.46900847149271, -38.068479264767944, -36.70195727063581, -35.36746735811498, -34.063034396224154, -32.78536874375959, -31.468633548362497, -30.173781283368523, -28.901550623295833, -27.652680242662605, -26.427908815987003, -25.227975017787198, -24.053617522581362, -22.905575004887663, -21.78458613922428, -20.691389600109368, -19.6267240620611, -18.59132819959765, -17.585940687237187, -16.61130019949788, -15.668145410897901, -14.794309398041404, -13.950185357101319, -13.135122694543167, -12.348470816832473, -11.589579130434757, -10.857797041815543, -10.152473957440362, -9.47295928377472, -8.818602427284151, -8.188752794434174, -7.582759791690314, -6.932160232642741, -6.304213548518645, -5.700899806714059, -5.124199074625014, -4.576091419647542, -4.0585569091776765, -3.573575610611449, -3.123127591344892, -2.709192918774041, -2.3337516602949204, -1.9987838833035663, -1.706269655196011, -1.458189043368287, -1.2565221152164259, -1.1032489381364605, -1.0003495795244224, -0.9498041067763443, -0.953592587288258, -1.0136950884561957, -1.1320916776761898, -1.31076242234427, -1.5516873898564725, -1.8568466476088277, -2.2282202629973678, -2.667788303418125};

  const double b1list[100] = {-0.5693500504085347, -0.576434151837257, -0.5827940588889807, -0.5875005969333106, -0.5915255507494274, -0.5970658827548452, -0.6057016775611604, -0.6053160995270499, -0.6070988602490128, -0.6110941958474475, -0.6140262971912503, -0.6172989788502661, -0.6206089494421659, -0.6234488672149939, -0.6258813681301192, -0.6277498501118401, -0.628897711084455, -0.6291683489722621, -0.6308371388180571, -0.6331613723912005, -0.6359581432020665, -0.6393457874435003, -0.6434426413083475, -0.6478966252652467, -0.6524003236766366, -0.6571425380351175, -0.6620695856777987, -0.6671277839417898, -0.6722634501642004, -0.6774229016821397, -0.6825524558327174, -0.6875984299530429, -0.6925071413802255, -0.6961049747777461, -0.6996298196832709, -0.7032241904115, -0.7069570372346382, -0.7108973104248908, -0.7151139602544623, -0.720064149835206, -0.7258460785616149, -0.7320745022572427, -0.7387427060855165, -0.7458439752098636, -0.7533715947937112, -0.7613188500004865, -0.7696790259936167, -0.778445407936529, -0.7876035918776259, -0.7967733093076456, -0.8063159350249182, -0.8162406283658405, -0.8265565486668099, -0.8372728552642233, -0.8483987074944775, -0.8599432646939698, -0.871915686199097, -0.884325131346256, -0.8971807594718443, -0.9104917299122585, -0.9242672020038958, -0.938516335083153, -0.9532482884864274, -0.9684722215501159, -0.9837886712042127, -0.9996512322126136, -1.0160843677317775, -1.0331125409181625, -1.0507602149282274, -1.0690518529184303, -1.08801191804523, -1.1076648734650847, -1.1280351823344532, -1.1491473078097938, -1.171025713047565, -1.1939017528688474, -1.2175927015073384, -1.2421149961497486, -1.2674850739827876, -1.2937193721931661, -1.320834327967594, -1.3488463784927816, -1.377771960955439, -1.4076275125422764, -1.4384294704400042, -1.4701942718353325, -1.5029383539149717, -1.5366781538656316, -1.5714301088740226, -1.6072106561268549, -1.6440362328108387, -1.6819232761126843, -1.7208882232191016, -1.7609475113168012, -1.802117577592493, -1.8444148592328868, -1.8878557934246936, -1.9324568173546235, -1.9782343682093864, -2.0252048831756926};

  const double b2list[100] = {1.6514745488753086, 1.6733593678301482, 1.687838328174986, 1.7031979370575185, 1.712831020475929, 1.7266279186283089, 1.7581796869631672, 1.7499867318965905, 1.7518398412177276, 1.7634468469918447, 1.770740685014047, 1.779639998248617, 1.788893228893931, 1.7964389585725973, 1.8024779675983216, 1.8063408969988246, 1.8073583878018264, 1.8048610810350476, 1.8077536017385523, 1.8130592620946404, 1.8200051503317694, 1.8290042591854243, 1.8404695813910905, 1.8531136802761718, 1.8658884601302625, 1.8795206835759333, 1.8938451458175296, 1.908696642059398, 1.923909967505884, 1.9393199173613336, 1.9547612868300928, 1.9700688711165077, 1.9850774654249244, 1.9960816492286722, 2.0069990026439246, 2.0182845185376195, 2.0301606467081585, 2.042849836953944, 2.056574539073377, 2.0726094682881886, 2.0912563552368337, 2.111506538675078, 2.1333773752936596, 2.1568862217833127, 2.1820504348347747, 2.208887371138781, 2.2374143873860683, 2.2676488402673725, 2.299578947658098, 2.3318064402049257, 2.3657435940219895, 2.4014679228653333, 2.439056940491002, 2.478588160655038, 2.520139097113486, 2.56378726362239, 2.609610173937794, 2.657685341815741, 2.708090281012276, 2.760902505283443, 2.8161995283852854, 2.874058864073847, 2.9345580261051722, 2.9977745282353045, 3.0615596376827443, 3.128412230361185, 3.1984931979224718, 3.271963432018449, 3.3489838243009618, 3.429715266421855, 3.5143186500329735, 3.602954866786163, 3.6957848083332676, 3.792969366326133, 3.894669432416603, 4.002044979667201, 4.114256382857298, 4.23142577617022, 4.353675293789295, 4.48112706989785, 4.613903238679211, 4.752125934316707, 4.895917290993665, 5.04539944289341, 5.200694524199272, 5.361924669094577, 5.529212011762653, 5.702678686386827, 5.882446827150425, 6.068638568236777, 6.261376043829206, 6.460781388111044, 6.666976735265615, 6.880084219476247, 7.100225974926268, 7.327524135799001, 7.5621008362777795, 7.804078210545929, 8.053578392786775, 8.310723517183646};

  const double a3stab[50][69] = {
{1298.87913, 1261.29637, 1227.72908, 1197.25947, 1168.96973, 1141.94207, 1115.2587, 1088.00181, 1059.25362, 1022.90618, 985.307903, 947.617046, 910.991876, 881.346609, 853.181173, 825.751447, 798.313312, 767.87765, 736.843339, 705.36426, 673.594293, 642.125047, 610.497581, 578.690684, 546.683144, 514.497622, 482.051486, 449.305976, 416.222332, 382.495195, 348.459044, 314.181759, 279.731222, 245.411071, 210.959123, 176.348956, 141.554145, 106.628243, 71.43286, 35.9095836, 0., -36.1519021, -72.887899, -110.347365, -148.669675, -186.733245, -226.44279, -268.442067, -313.374835, -358.91977, -409.871742, -468.06054, -535.315953, -602.768787, -687.227409, -794.8012, -931.599543, -1073.78457, -1269.39181, -1536.50955, -1893.22609, -2326.65745, -2947.80866, -3741.68954, -4974.29064, -6365.10282, -9538.63496, -15643.1414, -25826.8766}, 
{1196.32002, 1167.06104, 1137.27475, 1107.12301, 1076.76769, 1046.37066, 1016.09379, 986.098937, 956.547978, 928.8469, 901.415801, 873.918898, 846.020409, 814.568078, 783.16919, 752.614552, 723.694973, 702.877684, 683.006503, 662.601671, 640.183429, 607.718815, 572.902552, 536.876163, 500.78117, 469.351999, 438.700105, 408.529845, 378.545579, 347.703945, 316.75611, 285.705522, 254.555628, 223.39562, 192.108901, 160.664621, 129.031928, 97.2453174, 65.1824526, 32.7863434, 0., -33.0547348, -66.6852171, -101.01997, -136.187517, -170.971669, -207.383546, -246.089556, -287.756108, -330.313665, -378.258958, -433.352771, -497.35589, -562.191692, -643.393334, -746.656564, -877.67713, -1013.79966, -1200.41146, -1454.54875, -1793.24772, -2204.40576, -2792.4755, -3542.23306, -4707.36116, -6024.93464, -9029.55323, -14806.4354, -24440.7998}, 
{1111.34089, 1088.04208, 1061.93086, 1033.65502, 1003.86233, 973.200551, 942.317471, 911.86086, 882.478492, 858.393797, 835.248631, 812.260507, 788.646935, 759.48951, 729.796029, 700.438372, 672.288418, 651.387883, 631.370875, 611.041341, 589.203225, 559.502697, 527.964588, 495.455951, 462.843844, 433.627925, 404.989603, 376.74289, 348.701797, 320.308767, 291.898012, 263.432172, 234.873889, 206.130982, 177.242843, 148.194039, 118.969142, 89.6533759, 60.0903922, 30.2244979, 0., -30.4005164, -61.3663335, -93.0484555, -125.597886, -157.893074, -191.866601, -228.178495, -267.488782, -307.871972, -353.607815, -406.390546, -467.914398, -530.2754, -608.605274, -708.437534, -835.305696, -967.651586, -1148.93709, -1395.53239, -1723.80769, -2122.58454, -2690.87904, -3414.09589, -4533.73512, -5812.62405, -8672.18555, -14119.1936, -23160.4223}, 
{1039.53251, 1020.2822, 997.603313, 972.198968, 944.772298, 916.02643, 886.664491, 857.389611, 828.904917, 804.925294, 781.93741, 759.439688, 736.930555, 712.636669, 687.836926, 662.538457, 636.748394, 609.788544, 582.625489, 555.540492, 528.814811, 504.169828, 479.870634, 455.622441, 431.13046, 404.996091, 378.469884, 351.698574, 324.8289, 298.520437, 272.201947, 245.815032, 219.301294, 192.43116, 165.385874, 138.17551, 110.810137, 83.4601202, 55.9111218, 28.1090968, 0., -28.1232207, -56.8004007, -86.2243824, -116.588008, -146.981092, -179.140716, -213.700934, -251.2958, -290.060701, -334.127824, -385.130693, -444.702827, -504.716717, -580.471327, -677.504592, -801.354447, -931.385611, -1109.77852, -1352.54038, -1675.67843, -2068.96254, -2627.11195, -3337.00139, -4427.97475, -5693.23265, -8428.0796, -13552.2941, -21985.6547}, 
{976.485655, 959.824144, 940.197989, 918.098323, 894.016279, 868.442993, 841.869596, 814.787222, 787.687004, 761.819814, 736.613152, 712.254257, 688.930369, 670.53532, 652.067116, 632.23036, 609.729651, 574.459963, 537.459377, 500.956344, 467.179318, 448.413342, 432.807639, 418.568025, 403.900314, 381.127188, 356.690849, 331.150363, 305.064801, 280.549729, 255.985115, 231.307428, 206.453135, 181.110154, 155.562924, 129.847331, 103.999264, 78.2732522, 52.3990841, 26.3251899, 0., -26.1568213, -52.856571, -80.3393122, -108.845108, -137.719354, -168.454647, -201.64892, -237.900101, -275.345863, -318.028498, -367.530041, -425.432527, -483.212445, -556.599593, -651.218222, -772.692584, -901.046982, -1077.7456, -1318.65266, -1639.63241, -2031.30849, -2585.26688, -3290.67293, -4364.64226, -5631.82205, -8258.783, -13076.6149, -20916.4077}, 
{917.791081, 902.710621, 885.620781, 866.696557, 846.112943, 824.044936, 800.66753, 776.155722, 750.684506, 722.455778, 694.406874, 667.502029, 642.705478, 629.711227, 617.261835, 602.829629, 583.886937, 541.782439, 496.561572, 452.146121, 412.457876, 398.926373, 390.962551, 385.485097, 379.412698, 359.691907, 337.202401, 312.851722, 287.547412, 264.607418, 241.564715, 218.362682, 194.944699, 170.981969, 146.79892, 122.447804, 97.9808711, 73.7004736, 49.308722, 24.7578271, 0., -24.435292, -49.4039968, -75.1848062, -102.056412, -129.591492, -159.057156, -191.014502, -226.024627, -262.193469, -303.519346, -351.545419, -407.814847, -463.459392, -534.598171, -626.938905, -746.189313, -872.680952, -1047.64817, -1286.94915, -1606.44207, -1997.39115, -2549.4365, -3254.83387, -4318.29988, -5593.45384, -8125.84343, -12663.0342, -19952.5918}, 
{859.039559, 844.98436, 829.777581, 813.337142, 795.580959, 776.426953, 755.793041, 733.597141, 709.757173, 680.211609, 650.449591, 621.980817, 596.314984, 586.690157, 578.196318, 567.651816, 551.875002, 508.136269, 460.621105, 413.967044, 372.811617, 362.402057, 358.522317, 357.566049, 355.926907, 338.360939, 317.554443, 294.556114, 270.414647, 248.904279, 227.257944, 205.434118, 183.391278, 160.860605, 138.11879, 115.215226, 92.1993096, 69.3494866, 46.3944783, 23.2920582, 0., -22.8926065, -46.3118306, -70.5524255, -95.9091441, -122.081137, -150.197001, -180.789729, -214.392317, -249.069531, -288.809882, -335.133657, -389.561139, -443.154361, -512.075161, -602.027124, -718.713836, -842.332775, -1014.29608, -1250.50977, -1566.87987, -1954.97924, -2503.71346, -3209.20757, -4263.50982, -5543.18962, -7990.80852, -12282.43, -19094.1177}, 
{795.821854, 782.688084, 768.57428, 753.36355, 736.939001, 719.18374, 699.980874, 679.213512, 656.764759, 628.465728, 599.872319, 572.488437, 547.817986, 537.997872, 529.6458, 520.012472, 506.348594, 469.901748, 430.32701, 391.276333, 356.401671, 345.533527, 339.673884, 336.003275, 331.702235, 314.804975, 295.296878, 274.017005, 251.804419, 231.651088, 211.382003, 190.97506, 170.40816, 149.560065, 128.547461, 107.387899, 86.098929, 64.8279929, 43.4107958, 21.8129331, 0., -21.4627385, -43.4492248, -66.2337315, -90.0905311, -114.671921, -141.12294, -169.966649, -201.726112, -234.440059, -272.109617, -316.251583, -368.382751, -419.994158, -486.638661, -573.843361, -687.135356, -806.047703, -972.49916, -1202.41444, -1511.71827, -1891.84153, -2432.19041, -3133.5174, -4174.83432, -5446.09099, -7815.22592, -11905.6805, -18340.896}, 
{723.728735, 711.864516, 697.916767, 682.119252, 664.705738, 645.909991, 625.965777, 605.106863, 583.567013, 560.596558, 537.806075, 515.822704, 495.273586, 480.160139, 466.385516, 453.227147, 439.962463, 423.459176, 406.36832, 388.931212, 371.389169, 355.013918, 338.6042, 321.989168, 304.997976, 286.694707, 267.979608, 248.987859, 229.854639, 211.058618, 192.254089, 173.438833, 154.610635, 135.894352, 117.109862, 98.2041202, 79.1240795, 59.7436944, 40.1121174, 20.2055015, 0., -20.0796619, -40.6853319, -62.0202854, -84.2877981, -106.847476, -131.083733, -157.53731, -186.748953, -216.771064, -251.628063, -292.856026, -341.991034, -391.675587, -455.896772, -539.748099, -648.323078, -759.870989, -917.067264, -1135.7431, -1431.72971, -1795.74674, -2318.96003, -3007.48672, -4026.8356, -5267.21955, -7560.6433, -11503.6636, -17692.8373}, 
{670.702392, 658.608735, 645.231211, 630.627242, 614.85425, 597.969656, 580.030884, 561.095353, 541.220486, 518.602027, 495.903746, 473.927738, 453.476095, 438.689549, 425.696101, 413.962387, 402.955046, 391.824136, 380.479507, 368.51443, 355.522175, 339.17332, 321.752905, 303.623275, 285.146778, 267.626973, 250.10851, 232.577251, 215.019058, 197.328079, 179.618576, 161.913096, 144.234187, 126.856941, 109.450342, 91.9359207, 74.2352057, 56.1190544, 37.7199375, 19.0196532, 0., -18.9413868, -38.4048807, -58.5750179, -79.636335, -100.928376, -123.818666, -148.82974, -176.48413, -204.795257, -237.798413, -277.019773, -323.985518, -371.643004, -433.528758, -514.600485, -619.815891, -726.601121, -878.458066, -1091.35706, -1381.26842, -1738.68566, -2256.0096, -2944.21016, -3956.54775, -5195.43966, -7424.38522, -11191.6133, -17045.3526}, 
{623.046755, 610.336599, 597.375254, 583.988785, 570.003256, 555.244732, 539.539276, 522.712953, 504.591828, 481.903784, 458.812338, 436.382826, 415.680585, 401.790876, 390.151139, 380.218743, 371.451052, 363.730969, 355.920111, 347.30563, 337.174678, 321.458633, 304.142732, 285.856435, 267.229206, 250.700686, 234.366085, 218.130789, 201.900189, 185.228258, 168.512366, 151.798468, 135.13252, 118.913825, 102.69365, 86.3766111, 69.8673234, 52.8685906, 35.5675654, 17.9495883, 0., -17.9060924, -36.3290217, -55.439354, -75.4076554, -95.5838756, -117.287444, -141.017173, -167.271875, -194.031478, -225.321234, -262.64751, -307.516673, -353.072705, -412.529312, -490.737815, -592.549535, -694.603535, -841.248297, -1048.62004, -1332.855, -1684.25015, -2196.45943, -2885.12536, -3891.89364, -5131.56784, -7300.31693, -10898.9685, -16428.3503}, 
{578.777874, 565.343763, 552.617764, 540.222136, 527.779136, 514.911024, 501.240056, 486.388491, 469.978587, 447.525103, 424.400795, 401.870923, 381.200742, 368.162933, 357.712364, 349.311322, 342.422098, 337.556737, 332.707869, 326.917882, 319.229161, 304.187697, 287.130832, 268.899514, 250.334688, 234.780471, 219.573371, 204.553063, 189.559225, 173.868813, 158.109315, 142.345496, 126.64212, 111.492196, 96.360951, 81.1418536, 65.7283731, 49.7810319, 33.5194242, 16.9301977, 0., -16.9275463, -34.3655585, -52.4701786, -71.3975492, -90.5102836, -111.075464, -133.566646, -158.457381, -183.715399, -213.322409, -248.754293, -291.486936, -334.842038, -391.711338, -466.832391, -564.942754, -661.929092, -802.920208, -1004.19402, -1282.02842, -1626.52282, -2132.4907, -2819.95411, -3819.61549, -5057.44347, -7167.66153, -10607.4294, -15833.9069}, 
{537.705661, 523.592259, 510.97538, 499.323946, 488.10688, 476.793102, 464.851536, 451.751104, 436.960726, 415.185168, 392.563172, 370.469324, 350.278208, 338.111834, 328.698393, 321.5135, 316.032768, 313.155404, 310.363995, 306.564723, 300.663765, 286.44217, 269.9813, 252.237386, 234.16666, 219.614044, 205.491601, 191.600088, 177.74026, 163.024266, 148.21691, 133.394391, 118.632905, 104.472791, 90.3404513, 76.1264269, 61.7212598, 46.7790037, 31.5212837, 15.9332368, 0., -15.974157, -32.4508828, -49.5726937, -67.4821057, -85.5612143, -105.017124, -126.29652, -149.846085, -173.651974, -201.605614, -235.1379, -275.67973, -316.701541, -370.778872, -442.526802, -536.56041, -628.112094, -762.932689, -957.39035, -1227.85323, -1564.25632, -2062.26726, -2745.81353, -3735.46447, -4967.31563, -7018.71577, -10306.7003, -15248.3048}, 
{499.640027, 485.044119, 472.464742, 461.290869, 450.911477, 440.715539, 430.092032, 418.42993, 405.118208, 384.603164, 363.193521, 342.255325, 323.154622, 311.94369, 303.427849, 297.098651, 292.447647, 290.380935, 288.409705, 285.45969, 280.456626, 267.30371, 251.958228, 235.354929, 218.428561, 204.949119, 191.88201, 179.027883, 166.187391, 152.469136, 138.642637, 124.785364, 110.97479, 97.7363464, 84.5203574, 71.2251094, 57.7488886, 43.7851316, 29.5189138, 14.9304613, 0., -15.0143326, -30.5213867, -46.652101, -63.5374142, -80.5902821, -98.9468193, -119.025158, -141.243429, -163.646157, -189.974524, -221.596107, -259.878479, -298.401753, -349.435951, -417.463636, -506.967366, -592.686844, -720.744632, -907.520433, -1169.39395, -1496.20329, -1983.95294, -2659.82073, -3635.19172, -4855.43337, -6845.77638, -9986.48562, -14657.826}, 
{464.390881, 449.661375, 437.102489, 426.119557, 416.117917, 406.502905, 396.679857, 386.054108, 374.030995, 355.498276, 336.185896, 317.306223, 300.071624, 289.964612, 282.219351, 276.34015, 271.831319, 269.087297, 266.366212, 262.816322, 257.585886, 245.853978, 232.325711, 217.737017, 202.823827, 190.533415, 178.505829, 166.592464, 154.644712, 141.977941, 129.193979, 116.358629, 103.537694, 91.1635985, 78.7888757, 66.3326796, 53.7141645, 40.722041, 27.4580845, 13.8936269, 0., -14.0164816, -28.5134619, -43.6136023, -59.439564, -75.4511011, -92.6989449, -111.57092, -132.45485, -153.502899, -178.232816, -207.926688, -243.866605, -279.69321, -327.386613, -391.285482, -475.728482, -555.187646, -675.814927, -853.895645, -1105.71512, -1421.11637, -1895.71161, -2559.09281, -3514.54842, -4716.04577, -6641.14011, -9636.48968, -14048.7528}, 
{431.768135, 417.406059, 404.905261, 393.806663, 383.651192, 373.979771, 364.333324, 354.252775, 343.279049, 327.589687, 311.434349, 295.699313, 281.270856, 272.480712, 265.391519, 259.511373, 254.348369, 249.128454, 243.754732, 237.848156, 231.029682, 221.17463, 210.347843, 198.868528, 187.055894, 176.114645, 165.124294, 154.049849, 142.856319, 131.325202, 119.678423, 107.954395, 96.1915321, 84.6352843, 73.0342128, 61.3439161, 49.5199927, 37.5123573, 25.2845654, 12.794489, 0., -12.9490123, -26.3635003, -40.3623994, -55.0646445, -69.9972857, -86.1078972, -103.752168, -123.285788, -143.027154, -166.184163, -193.927421, -227.427532, -260.326451, -304.334894, -363.634929, -442.40862, -515.148803, -627.602467, -795.827371, -1035.88127, -1337.74821, -1795.70711, -2440.7469, -3369.28571, -4543.40188, -6397.10368, -9246.41689, -13407.3673}, 
{397.382476, 384.161063, 372.406617, 361.789178, 351.978785, 342.645477, 333.459295, 324.090277, 314.208462, 301.130105, 287.820545, 274.891334, 262.954027, 255.327872, 248.833648, 242.999831, 237.354894, 230.264221, 222.884616, 215.20979, 207.233455, 198.951905, 190.355236, 181.436127, 172.187257, 162.452347, 152.432615, 142.180322, 131.747728, 121.278447, 110.696847, 100.018645, 89.2595619, 78.5197505, 67.6967231, 56.7724257, 45.7288044, 34.6183092, 23.3241805, 11.8001626, 0., -11.943738, -24.3276619, -37.2695573, -50.8872099, -64.7777681, -79.7879094, -96.2436743, -114.471103, -132.961668, -154.609806, -180.475383, -211.618268, -241.71635, -282.164266, -336.974677, -410.160241, -476.397825, -580.770198, -739.024339, -966.907224, -1254.16741, -1694.54713, -2321.20376, -3222.44556, -4362.79135, -6146.7668, -8868.89588, -12823.7025}, 
{359.974686, 348.352, 338.125236, 328.975588, 320.58425, 312.632416, 304.801281, 296.772038, 288.225882, 276.655174, 264.805473, 253.233508, 242.496004, 235.582587, 229.643928, 224.263593, 219.02515, 212.45994, 205.624648, 198.523735, 191.16166, 183.591405, 175.749496, 167.620982, 159.190914, 150.290229, 141.119734, 131.726121, 122.156087, 112.539351, 102.806371, 92.9706295, 83.0456113, 73.124859, 63.1097732, 52.9818134, 42.7224396, 32.3549175, 21.8021786, 11.0289603, 0., -11.1220613, -22.6495511, -34.6968934, -47.3785121, -60.3411594, -74.354, -89.7185264, -106.736232, -124.027317, -144.247083, -168.369539, -197.368694, -225.064687, -262.446946, -313.351026, -381.612486, -442.034343, -539.097709, -688.251154, -904.943253, -1178.0799, -1602.7377, -2215.93586, -3095.8159, -4202.28901, -5927.06542, -8553.90562, -12366.5701}, 
{324.284553, 314.363487, 305.717204, 298.046437, 291.051921, 284.434388, 277.894571, 271.133205, 263.851022, 253.758508, 243.342743, 233.100562, 223.528797, 217.260142, 211.801225, 206.794534, 201.882558, 195.779693, 189.427757, 182.840476, 176.031575, 169.10102, 161.9418, 154.533147, 146.854291, 138.731414, 130.358013, 121.774536, 113.021431, 104.2128, 95.2859787, 86.2519542, 77.1217145, 67.9818174, 58.7374528, 49.3693806, 39.8583605, 30.1989457, 20.3525856, 10.2945228, 0., -10.342512, -21.0562901, -32.2513829, -44.0378392, -56.1058801, -69.1533129, -83.4581171, -99.2982726, -115.426321, -134.255855, -156.675029, -183.572, -208.924082, -243.294606, -290.336061, -353.700939, -408.371479, -498.138522, -638.12266, -843.444484, -1102.14389, -1510.58355, -2109.2871, -2967.15339, -4038.77567, -5707.59445, -8249.83842, -11941.7363}, 
{290.624535, 282.418333, 275.336854, 269.105352, 263.449083, 258.093303, 252.763265, 247.184226, 241.08144, 232.406317, 223.367496, 214.39977, 205.937932, 200.246338, 195.198394, 190.497067, 185.845325, 180.150535, 174.229506, 168.103447, 161.793567, 155.43567, 148.890531, 142.133519, 135.140006, 127.739783, 120.11203, 112.290346, 104.308333, 96.2624774, 88.09834, 79.8243672, 71.4490055, 63.0513927, 54.5410077, 45.8980203, 37.1026006, 28.1249184, 18.959144, 9.58944784, 0., -9.60482061, -19.549136, -29.9368594, -40.8719042, -52.0789138, -64.1927792, -77.4691216, -92.1635624, -107.160693, -124.635576, -145.392245, -170.234732, -193.315202, -224.750304, -268.004818, -326.543524, -375.558951, -458.097032, -588.931449, -782.835883, -1027.0032, -1418.94953, -2002.20341, -2837.66181, -3874.62585, -5491.41154, -7959.08117, -11548.697}, 
{259.307093, 252.739349, 247.138517, 242.255957, 237.843023, 233.651073, 229.431464, 224.935552, 219.914693, 212.564812, 204.814872, 197.038402, 189.608933, 184.426978, 179.72829, 175.275605, 170.831659, 165.499518, 159.965458, 154.256084, 148.398001, 142.550276, 136.554067, 130.382995, 124.010678, 117.279218, 110.346361, 103.238334, 95.9813658, 88.6520651, 81.2061239, 73.6496167, 65.9886178, 58.2943521, 50.4816833, 42.5306256, 34.4211932, 26.1073599, 17.6055963, 8.90633292, 0., -8.90871776, -18.1293461, -27.7571564, -37.8874201, -48.2672439, -59.4793298, -71.758215, -85.3384368, -99.2324482, -115.385705, -134.521577, -157.363436, -178.258714, -206.857097, -246.432332, -300.258165, -343.746479, -419.177631, -540.970111, -723.542412, -953.301677, -1328.70045, -1895.63073, -2708.54494, -3712.21403, -5281.5743, -7684.02071, -11186.9482}, 
{230.644686, 225.549343, 221.276528, 217.601877, 214.301026, 211.14961, 207.923267, 204.397631, 200.348339, 194.200203, 187.620012, 180.923731, 174.427326, 169.687863, 165.283766, 161.034559, 156.759767, 151.753696, 146.571178, 141.241823, 135.795243, 130.399758, 124.890791, 119.24247, 113.428926, 107.313599, 101.025583, 94.5832827, 88.0051028, 81.3452456, 74.5719997, 67.6894509, 60.7016853, 53.6714628, 46.5207255, 39.2300897, 31.7801716, 24.1207947, 16.2756849, 8.23777562, 0., -8.25393396, -16.7981776, -25.7161074, -35.0911002, -44.6778543, -55.0198957, -66.3320723, -78.829232, -91.6435998, -106.505696, -124.063416, -144.964658, -163.775287, -189.658045, -225.69364, -274.962784, -313.083784, -381.584712, -494.531239, -665.989037, -881.683134, -1240.70113, -1790.51499, -2581.00657, -3553.91473, -5081.14038, -7427.04393, -10855.9858}, 
{204.949773, 201.071126, 197.905217, 195.246737, 192.890376, 190.630825, 188.262774, 185.580914, 182.379935, 177.278699, 171.718057, 165.963028, 160.278636, 155.914796, 151.757677, 147.678343, 143.547856, 138.840123, 133.982228, 129.004099, 123.935661, 118.939037, 113.859081, 108.672842, 103.357369, 97.8068086, 92.1142731, 86.2899726, 80.3441169, 74.3057015, 68.1586367, 61.9056181, 55.5493414, 49.1434918, 42.6193799, 35.9593056, 29.1455691, 22.1397473, 14.9531524, 7.57637353, 0., -7.6401998, -15.5568876, -23.817546, -32.4896574, -41.3177284, -50.8214078, -61.1973687, -72.6422837, -84.3961618, -97.9950056, -114.018154, -133.044944, -149.885588, -173.196203, -205.863779, -250.775306, -283.720583, -345.522668, -449.907425, -610.600723, -812.791405, -1155.8164, -1687.80212, -2456.25049, -3402.10245, -4893.16741, -7190.53768, -10555.3056}, 
{182.534813, 179.527507, 177.178917, 175.294163, 173.678361, 172.136629, 170.474086, 168.49585, 166.007037, 161.766511, 157.044147, 152.063566, 147.048389, 142.99358, 139.042878, 135.111368, 131.114134, 126.685851, 122.134173, 117.486347, 112.76962, 108.123033, 103.417319, 98.6350069, 93.7586259, 88.7227275, 83.5770077, 78.3231858, 72.9629811, 67.4971152, 61.9287041, 56.2598663, 50.4927201, 44.6712063, 38.738892, 32.6811664, 26.4834191, 20.1387421, 13.621741, 6.91472423, 0., -7.06724582, -14.4067335, -22.0653056, -30.0898049, -38.19385, -46.8907973, -56.3607792, -66.783928, -77.4921481, -89.8530911, -104.386181, -121.610841, -136.610285, -157.514631, -187.017785, -227.813654, -255.806599, -311.195892, -407.391261, -557.802432, -747.270316, -1074.91108, -1588.43805, -2335.48049, -3259.15168, -4720.71301, -6976.88882, -10284.4035}, 
{165.290507, 162.286376, 160.078255, 158.450848, 157.188853, 156.076972, 154.899907, 153.442357, 151.489025, 147.77082, 143.54775, 139.026033, 134.411887, 130.567225, 126.78029, 122.995023, 119.155364, 115.021951, 110.795345, 106.492805, 102.131591, 97.7874608, 93.395776, 88.9503958, 84.4451802, 79.8803286, 75.2408251, 70.5179937, 65.7031583, 60.7565616, 55.7130412, 50.5763533, 45.3502544, 40.0687128, 34.693188, 29.2153516, 23.6268752, 17.9944365, 12.2046983, 6.21932983, 0., -6.54855977, -13.3849676, -20.5247794, -27.983551, -35.3990376, -43.315716, -51.9002622, -61.3193524, -70.9677412, -82.0927948, -95.1699579, -110.674675, -124.050902, -142.81816800000001, -169.464514, -206.47798, -229.664784, -278.867519, -367.246955, -507.963863, -685.99481, -999.053177, -1492.11906, -2218.17943, -3129.51904, -4573.22466, -6797.68259, -10051.2791}, 
{152.137356, 148.560618, 146.028799, 144.288529, 143.086442, 142.169167, 141.283337, 140.175584, 138.592538, 135.133503, 131.151371, 126.851704, 122.440065, 118.735667, 115.084965, 111.448061, 107.78506, 103.948328, 100.048801, 96.0896769, 92.0741556, 87.9923766, 83.8658206, 79.7029095, 75.5120655, 71.371298, 67.1916066, 62.9535786, 58.6378009, 54.159387, 49.5905871, 44.9381778, 40.2089354, 35.4261284, 30.5734448, 25.6510642, 20.6591665, 15.7676707, 10.7391215, 5.50580328, 0., -6.0717388, -12.4611006, -19.1455073, -26.1023811, -32.8645712, -40.0319015, -47.7596232, -56.2029876, -64.7948193, -74.7017665, -86.3680509, -100.237894, -112.156793, -129.007185, -153.07278, -186.637291, -205.292071, -248.690132, -329.792132, -461.558725, -629.404881, -928.92831, -1400.89323, -2107.06889, -3013.74267, -4448.12437, -6647.19318, -9847.92831}, 
{140.781753, 136.682234, 133.820464, 131.915018, 130.684473, 129.847403, 129.122385, 128.227993, 126.882804, 123.587991, 119.766493, 115.623843, 111.365579, 107.785616, 104.265758, 100.776187, 97.2870864, 93.7156472, 90.1062419, 86.4502507, 82.7390537, 78.8895359, 74.9973708, 71.0837365, 67.1698115, 63.3987972, 59.6210392, 55.8089063, 51.9347676, 47.8783341, 43.7416956, 39.5342838, 35.2655307, 30.9463936, 26.5841688, 22.1876782, 17.7657435, 13.5945717, 9.30864537, 4.80983222, 0., -5.61379701, -11.5769522, -17.8296726, -24.3121648, -30.4556994, -36.9129944, -43.8278312, -51.3439913, -58.9189753, -67.6573579, -77.9774331, -90.2974945, -100.815114, -115.857596, -137.531524, -167.94348, -182.552124, -220.771131, -295.366252, -419.10324, -577.762732, -865.065819, -1317.76997, -2006.19433, -2910.75903, -4337.91853, -6512.61869, -9659.80539}, 
{130.494647, 126.16847, 123.140568, 121.124242, 119.832792, 118.979516, 118.277715, 117.440689, 116.181737, 112.994264, 109.299423, 105.29847, 101.192663, 97.752079, 94.3816266, 91.0550345, 87.7460315, 84.3994753, 81.0295143, 77.6214257, 74.1604867, 70.521654, 66.8446537, 63.1588911, 59.4937719, 56.0305044, 52.5859699, 49.1288527, 45.6278369, 41.946012, 38.1998949, 34.400408, 30.5584734, 26.6696049, 22.7662965, 18.8656338, 14.9847023, 11.5026577, 7.92968738, 4.13804898, 0., -5.16541025, -10.7114327, -16.5445265, -22.5711508, -28.1368435, -33.9333535, -40.0915086, -46.7421363, -53.3523053, -60.9821059, -70.0278696, -80.8859276, -90.0081075, -103.313046, -122.774876, -150.36773, -161.586559, -195.476351, -264.602911, -381.532045, -532.136331, -809.061255, -1246.06057, -1920.26372, -2824.73561, -4245.42361, -6395.05448, -9486.35496}, 
{120.218374, 116.363283, 113.606148, 111.699338, 110.395219, 109.44616, 108.604528, 107.622691, 106.253017, 103.165789, 99.6282923, 95.8257276, 91.9432956, 88.7062903, 85.5437823, 82.4249355, 79.3189135, 76.1411207, 72.935984, 69.6941707, 66.4063484, 62.9681364, 59.5032696, 56.0404345, 52.6083179, 49.374891, 46.1738422, 42.9781445, 39.7607707, 36.4023776, 33.0051808, 29.5790797, 26.1339737, 22.6437028, 19.1686496, 15.7331373, 12.361489, 9.52475412, 6.62183918, 3.4983769, 0., -4.71442103, -9.83778635, -15.2497584, -20.8299995, -25.8728603, -31.0774398, -36.5575249, -42.4269028, -48.1422156, -54.7372534, -62.5886611, -72.0730838, -79.7205999, -91.2930477, -108.705699, -133.873825, -142.639015, -173.419696, -238.560614, -350.406513, -494.295262, -763.592234, -1191.33966, -1857.42283, -2763.52395, -4177.2945, -6299.89362, -9332.48046}, 
{110.75343, 107.590506, 105.231651, 103.484465, 102.156549, 101.055506, 99.9889362, 98.7644409, 97.1896215, 94.1790465, 90.790563, 87.1889854, 83.539128, 80.5146295, 77.5679503, 74.660375, 71.7531881, 68.7101427, 65.6290678, 62.5102607, 59.3540188, 56.1009217, 52.8348715, 49.5800526, 46.3606494, 43.309762, 40.2990927, 37.3092595, 34.3208804, 31.2458379, 28.1609798, 25.0744183, 21.9942657, 18.872276, 15.7954629, 12.7944818, 9.89998836, 7.66367431, 5.38674428, 2.89143919, 0., -4.26469321, -8.96129429, -13.9478179, -19.0822789, -23.6261418, -28.2725914, -33.1182624, -38.2597896, -43.1364487, -48.7651769, -55.5055526, -63.7171541, -69.9228172, -79.8535597, -95.4036569, -118.467384, -125.376221, -153.812357, -215.895186, -323.744101, -462.041787, -725.217765, -1146.38373, -1806.37838, -2714.14396, -4118.48155, -6208.22748, -9172.21812}, 
{102.109544, 99.7569725, 97.8623984, 96.2978131, 94.9352071, 93.6465714, 92.3038964, 90.7791731, 88.9443922, 85.9985892, 82.7558921, 79.3574736, 75.9445064, 73.1306553, 70.3956043, 67.6915293, 64.9706062, 62.0343038, 59.0457882, 56.0175184, 52.9619535, 49.8802388, 46.8006721, 43.7402377, 40.7159199, 37.8110855, 34.9497834, 32.1224447, 29.3195008, 26.4924822, 23.6862811, 20.9068892, 18.1602981, 15.3777876, 12.6699457, 10.0726484, 7.62177209, 5.93500366, 4.23368424, 2.32096572, 0., -3.81327259, -8.07360397, -12.6229576, -17.3032971, -21.3601287, -25.4704559, -29.714825, -34.1737824, -38.2677164, -43.0013945, -48.7194263, -55.7664212, -60.5911985, -68.9924741, -82.8731736, -104.136223, -109.710419, -136.462466, -196.284943, -301.070425, -434.866515, -693.100714, -1109.41568, -1764.10962, -2672.48073, -4063.17178, -6111.48937, -8992.74014}, 
{94.2964493, 92.7695158, 91.3437153, 89.9575739, 88.5496178, 87.0583731, 85.422366, 83.5801226, 81.4701691, 78.5889691, 75.4939364, 72.300422, 69.1237773, 66.5079261, 63.9682178, 61.4485745, 58.8929184, 56.0413661, 53.1231674, 50.1637667, 47.1886082, 44.2663165, 41.3618836, 38.4834818, 35.6392833, 32.8548298, 30.1139763, 27.4179471, 24.7679666, 22.1583995, 19.6000738, 17.0969577, 14.6530196, 12.1827008, 9.81530731, 7.59061801, 5.54841194, 4.35432756, 3.17194063, 1.79068637, 0., -3.35720491, -7.16636284, -11.2594299, -15.4683622, -19.038262, -22.6226808, -26.2883164, -30.1018666, -33.4687308, -37.3814243, -42.1711641, -48.1691674, -51.7021831, -58.7076835, -71.1186729, -90.8681555, -95.5538487, -121.178158, -179.408204, -281.911104, -412.260057, -666.403945, -1078.65842, -1727.59583, -2634.41935, -4005.5522, -6001.11258, -8781.21869}, 
{87.7713906, 86.4960091, 85.1811683, 83.7968101, 82.3128763, 80.6993088, 78.9260494, 76.9630399, 74.7802222, 72.0511526, 69.1607128, 66.1973987, 63.2497064, 60.7494344, 58.3044555, 55.8659448, 53.3850774, 50.6152646, 47.784551, 44.9232173, 42.0615442, 39.301138, 36.5724235, 33.8771513, 31.2170719, 28.5794753, 25.9863566, 23.4452507, 20.9636922, 18.5573668, 16.2223978, 13.9630598, 11.783627, 9.60000198, 7.54017972, 5.64378337, 3.95043606, 3.11669301, 2.31847245, 1.34862468, 0., -2.86170759, -6.14562894, -9.68805125, -13.3252617, -16.3257928, -19.3207882, -22.3736369, -25.547728, -28.2713039, -31.4969588, -35.5421401, -40.7242953, -43.2310953, -49.1616752, -60.4853932, -79.1716077, -83.2430703, -108.194389, -165.573564, -266.928597, -396.429249, -647.75824, -1054.36803, -1693.5088, -2596.01584, -3934.58632, -5846.45193, -8468.84437}, 
{82.0309902, 80.8879567, 79.6096357, 78.1899457, 76.622805, 74.9021323, 73.0218458, 70.9758641, 68.7581056, 66.1932829, 63.5122026, 60.7764656, 58.0476726, 55.637402, 53.257286, 50.8689344, 48.4339568, 45.7340129, 42.9826428, 40.2134361, 37.4599828, 34.8511951, 32.2872114, 29.7634926, 27.2754996, 24.774658, 22.3180779, 19.918834, 17.5900013, 15.3641527, 13.2270657, 11.1840158, 9.24027873, 7.31904842, 5.54051468, 3.94278571, 2.56396973, 2.03191979, 1.55910133, 0.947724641, 0., -2.3603199, -5.10118416, -8.06899945, -11.1101725, -13.5280581, -15.9293355, -18.3776321, -20.9365754, -23.0678086, -25.6777369, -29.0707811, -33.5513625, -35.1556526, -40.1636214, -50.5869898, -68.4374784, -72.203814, -96.8299089, -153.736681, -254.34505, -383.750348, -632.350248, -1032.43624, -1659.18033, -2552.95931, -3850.45028, -5664.25527, -8106.97631}, 
{76.9707195, 75.8621394, 74.5611871, 73.0790334, 71.4268495, 69.6158063, 67.6570748, 65.5618261, 63.3412311, 60.9450827, 58.4704813, 55.9531493, 53.4288089, 51.087234, 48.7484747, 46.3866325, 43.9758091, 41.3365635, 38.6579569, 35.975508, 33.3247354, 30.8558433, 28.4437906, 26.0782215, 23.7487807, 21.3748538, 19.0444472, 16.775309, 14.585187, 12.5188412, 10.5562024, 8.70421353, 6.96981753, 5.28812718, 3.76664768, 2.44105406, 1.34702134, 1.06977331, 0.875616728, 0.580407118, 0., -1.85821402, -4.04652853, -6.42570206, -8.85649313, -10.6846322, -12.4919171, -14.3451174, -16.311003, -17.8933143, -19.9470619, -22.7642269, -26.6367909, -27.4525786, -31.6693904, -41.3408701, -58.5206616, -62.2399847, -86.7838771, -143.414952, -243.395824, -373.020461, -618.457409, -1010.56686, -1621.66818, -2501.0912, -3748.66476, -5451.59296, -7697.07988}, 
{72.4860498, 71.3353381, 69.967892, 68.4061262, 66.6724552, 64.7892935, 62.7790558, 60.6641567, 58.4670106, 56.2362749, 53.9576244, 51.6429764, 49.3042485, 47.0143357, 44.6997869, 42.3481283, 39.9468865, 37.3618686, 34.7510078, 32.1505183, 29.5966139, 27.2544382, 24.9797041, 22.7570538, 20.5711298, 18.3145381, 16.1007718, 13.9512876, 11.8875423, 9.96151604, 8.15193297, 6.46804054, 4.91908621, 3.45552525, 2.16891421, 1.09201745, 0.257599308, 0.200019021, 0.249808109, 0.239092969, 0., -1.36056209, -2.99516214, -4.78158666, -6.59762218, -7.83508972, -9.05212751, -10.3209083, -11.7136048, -12.7828905, -14.3282371, -16.629617, -19.9670027, -20.0985972, -23.6348501, -32.6644414, -49.2760512, -53.1554874, -77.7554513, -134.125772, -233.316278, -363.036691, -604.357163, -986.463671, -1578.03011, -2436.25295, -3624.75043, -5205.53533, -7240.62045}, 
{68.4724527, 67.2243335, 65.7618201, 64.1132768, 62.3070676, 60.3715567, 58.3351084, 56.2260867, 54.0728559, 51.9965824, 49.8957072, 47.7614734, 45.5851245, 43.3341123, 41.0329879, 38.682511, 36.2834413, 33.7488806, 31.2023098, 28.6795519, 26.21643, 23.9863352, 21.8324951, 19.7357053, 17.6767614, 15.5281868, 13.4223587, 11.3833816, 9.43536035, 7.63226104, 5.95638238, 4.41988452, 3.0349276, 1.76952966, 0.697649771, -0.150895095, -0.746287965, -0.607577646, -0.336535071, -0.0837969507, 0., -0.872536286, -1.96058503, -3.16008079, -4.36695815, -5.01900493, -5.65356114, -6.34982007, -7.18697498, -7.77160669, -8.84456592, -10.6740909, -13.5284199, -13.0704323, -16.0158687, -24.475111, -40.558541, -44.7542269, -69.4437895, -125.386536, -223.341773, -352.596143, -588.326948, -957.830483, -1525.32388, -2354.28602, -3474.22797, -4923.15275, -6739.06338}, 
{68.2547617, 64.4030038, 61.0624053, 58.1628557, 55.6342446, 53.4064616, 51.4093962, 49.572938, 47.8269765, 45.957039, 44.0951222, 42.2288605, 40.3458883, 38.4667234, 36.5329635, 34.5190896, 32.3995828, 29.9385344, 27.4049711, 24.8575296, 22.354847, 20.1021213, 17.9528037, 15.9069066, 13.9644425, 12.1438868, 10.4194037, 8.78362046, 7.22916423, 5.71567078, 4.28195543, 2.93384198, 1.67715426, 0.346381608, -0.812783878, -1.72598459, -2.3188629, -1.86839962, -1.20836336, -0.523861146, 0., -0.627400983, -1.46345147, -2.37105288, -3.21310661, -3.25085515, -3.18952241, -3.13267337, -3.18387301, -2.80547233, -2.9987359, -4.1237143, -6.54045811, -6.83371701, -10.6489628, -19.8563666, -36.3260991, -39.4859101, -62.6253602, -116.591589, -212.231736, -342.977808, -571.922343, -920.991834, -1451.388, -2230.85247, -3261.692, -4562.57208, -6152.15819}, 
{68.4539243, 61.4483189, 55.8628302, 51.4992593, 48.1594069, 45.6450741, 43.7580616, 42.3001706, 41.0732018, 39.3851077, 37.7290771, 36.1044501, 34.5105672, 33.1173068, 31.6852556, 30.1455386, 28.4292807, 26.0790926, 23.5700193, 20.9885914, 18.4213396, 16.0987041, 13.9057423, 11.871421, 10.0247072, 8.59280239, 7.32714494, 6.17740787, 5.09326418, 3.90357687, 2.727153, 1.56198962, 0.406083746, -1.08429602, -2.43273282, -3.50453824, -4.16502385, -3.33725149, -2.20568241, -1.01252808, 0., -0.497549335, -1.22525202, -1.99042368, -2.60037993, -2.04190726, -1.27106206, -0.423371592, 0.365636893, 1.81885925, 2.59897586, 2.22709019, 0.224305732, -1.26648705, -6.39468643, -16.6879037, -33.6737502, -35.6379855, -56.6468132, -107.524585, -199.095654, -329.609444, -547.615087, -868.664624, -1349.75899, -2057.6127, -2978.63703, -4119.31905, -5486.14585}, 
{65.4719479, 57.4208681, 51.1223649, 46.3288876, 42.792886, 40.2668095, 38.5031077, 37.2542302, 36.2726266, 34.6627386, 33.0842266, 31.5487434, 30.0679416, 28.9024083, 27.7152878, 26.4186592, 24.9246012, 22.6910859, 20.2659417, 17.7428902, 15.2156527, 12.9020086, 10.7219986, 8.71972156, 6.93927587, 5.70709501, 4.6720087, 3.76518156, 2.9177782, 1.89110459, 0.854127432, -0.194045212, -1.25430528, -2.72905117, -4.05706578, -5.07863845, -5.63405853, -4.49776766, -3.002242, -1.41410996, 0., -0.232086191, -0.689233061, -1.15585079, -1.41634954, -0.326326355, 1.02947019, 2.49510466, 3.91464163, 6.02483426, 7.41998307, 7.58707717, 6.01310568, 4.25526237, -1.09775015, -11.3870246, -27.9536538, -30.037432, -49.9212697, -97.7867791, -183.815572, -309.038656, -511.089458, -801.715921, -1230.32357, -1857.50504, -2661.44699, -3640.93761, -4794.76507}, 
{60.4371914, 52.8944439, 47.0061111, 42.5344545, 39.2417357, 36.8902165, 35.2421584, 34.0598231, 33.1054723, 31.4888045, 29.8856696, 28.3193545, 26.8131459, 25.6686144, 24.5194497, 23.2776252, 21.8551143, 19.7300162, 17.4217281, 15.0157733, 12.5976748, 10.34573, 8.215578, 6.25563199, 4.51430532, 3.33212144, 2.34853948, 1.49512871, 0.703458391, -0.282898639, -1.26917813, -2.24861224, -3.21443316, -4.5363881, -5.68058815, -6.48965946, -6.80622818, -5.41661253, -3.64226975, -1.74834914, 0., 0.144274795, 0.0900209157, 0.0494304766, 0.23469559, 1.80159795, 3.64130425, 5.58857079, 7.47815384, 9.93910759, 11.6941713, 12.260382, 11.1547768, 9.8495151, 5.12446293, -4.2853914, -19.6450596, -22.6973878, -42.0384195, -86.7420329, -165.882106, -281.689148, -463.767142, -721.985092, -1095.87966, -1636.34431, -2319.17767, -3138.9521, -4090.23996}, 
{54.4780137, 48.4428385, 43.6791707, 39.9986734, 37.21301, 35.1338437, 33.5728378, 32.3416556, 31.2519604, 29.5621783, 27.8585046, 26.1738973, 24.5413145, 23.2625116, 21.9941304, 20.6616102, 19.1903901, 17.1513851, 14.9663681, 12.7025882, 10.427294, 8.26356395, 6.20048559, 4.282976, 2.55595225, 1.31323844, 0.251281779, -0.684563488, -1.54894311, -2.57958539, -3.57482051, -4.51606121, -5.38472022, -6.45481104, -7.29810532, -7.77897552, -7.76179404, -6.16045051, -4.16999329, -2.03498794, 0., 0.606819966, 1.05792541, 1.54218478, 2.24846652, 4.24757601, 6.49367, 8.82284217, 11.0711862, 13.6883341, 15.6514259, 16.5511399, 15.9781543, 15.634255, 12.1545697, 4.28533388, -9.22721724, -13.6309911, -32.5879525, -73.754209, -144.785868, -247.984626, -407.069825, -631.311504, -949.22515, -1399.94534, -1960.88484, -2624.88686, -3384.79461}, 
{48.7227738, 44.6398442, 41.3066454, 38.604258, 36.4137625, 34.6162396, 33.0927698, 31.7244338, 30.3923122, 28.5817329, 26.7278303, 24.8699861, 23.0475818, 21.5306864, 20.0357196, 18.509788, 16.8999986, 14.9106943, 12.8288514, 10.6986821, 8.56439891, 6.48920575, 4.49072666, 2.60557732, 0.870373416, -0.504197236, -1.72521994, -2.82570787, -3.83867419, -4.96010825, -5.99485657, -6.91074186, -7.67558681, -8.43282424, -8.90442271, -8.98796086, -8.58101735, -6.79594598, -4.62964023, -2.2937687, 0., 1.13083567, 2.15989592, 3.23917679, 4.5206743, 6.91731822, 9.51579733, 12.1637342, 14.7087513, 17.3991688, 19.5216327, 20.7634863, 20.8120732, 21.7274658, 19.8751872, 13.9934892, 2.82062372, -2.85138016, -21.1595584, -58.1871699, -120.017474, -208.348795, -342.419193, -531.534524, -793.157935, -1154.12298, -1595.62429, -2110.26622, -2690.65314}, 
{44.2998304, 42.0592533, 40.0536373, 38.233922, 36.5510472, 34.9559528, 33.3995785, 31.832864, 30.2067492, 28.2463412, 26.2187454, 24.1652348, 22.1270823, 20.3197254, 18.5406064, 16.7613322, 14.9535099, 12.9634457, 10.9381676, 8.89940251, 6.86887778, 4.85635082, 2.90030642, 1.02725964, -0.73627446, -2.27482884, -3.6864212, -4.98011725, -6.16498272, -7.38561978, -8.46134313, -9.34700391, -9.99745327, -10.418932, -10.4943457, -10.1579897, -9.34415937, -7.38976337, -5.06543819, -2.54443376, 0., 1.69160825, 3.34134796, 5.05717119, 6.94702998, 9.71653493, 12.6369161, 15.5770622, 18.4058619, 21.1982666, 23.5346769, 25.2015565, 25.9853686, 28.2471315, 28.1689322, 24.5074123, 16.0192136, 9.6283068, -7.34292715, -39.4047781, -91.0675361, -163.20536, -271.236932, -424.493521, -630.475917, -904.692043, -1232.45178, -1606.61452, -2020.03966}, 
{31.3611692, 33.4138344, 34.6123875, 35.0543707, 34.837326, 34.0587958, 32.8163219, 31.2074467, 29.3297123, 27.404283, 25.3556298, 23.2318459, 21.0810247, 19.0533799, 17.0540359, 15.0902378, 13.1692306, 11.3600238, 9.58339213, 7.82187479, 6.05801097, 4.17746025, 2.29839333, 0.442101287, -1.37012481, -3.08825048, -4.73122545, -6.289256, -7.75254846, -9.29802506, -10.6544898, -11.7374626, -12.4624633, -12.5589911, -12.202995, -11.3844033, -10.0931443, -7.98083331, -5.51103668, -2.80900778, 0., 2.33855583, 4.71445469, 7.18331415, 9.80075179, 13.1193625, 16.4989955, 19.7964776, 22.8686353, 25.1976694, 27.1648828, 28.7769525, 30.0405554, 34.3745852, 37.0086155, 36.5844368, 31.7438393, 24.3660043, 8.56037482, -18.326215, -58.9469312, -111.659095, -192.003407, -312.956989, -466.821193, -649.230038, -863.162643, -1109.76185, -1390.17051}, 
{37.107078, 37.9494898, 38.1475584, 37.7748672, 36.9049999, 35.61154, 33.9680713, 32.0481773, 29.9254415, 27.7748251, 25.5279832, 23.2179484, 20.8777534, 18.5825868, 16.3064628, 14.0655515, 11.8760232, 9.83055734, 7.83821112, 5.88455095, 3.9551433, 1.96138975, -0.00731245436, -1.93573091, -3.80863323, -5.6288923, -7.35592834, -8.96726683, -10.4404333, -11.8918334, -13.1045604, -14.0005876, -14.5018883, -14.3327386, -13.6918882, -12.5803894, -10.9992944, -8.70169555, -6.03578926, -3.10181195, 0., 2.87088123, 5.82853593, 8.89213931, 12.0808666, 15.681342, 19.3383122, 22.9639726, 26.4705188, 29.5587418, 32.4368037, 35.1014621, 37.5494744, 42.4526156, 46.0626187, 47.3062344, 45.1102132, 39.9921181, 28.651562, 9.37897069, -19.5352305, -57.2665312, -113.12676, -196.55736, -298.235907, -413.117734, -537.811019, -668.510728, -801.411827}, 
{54.2496337, 50.3266662, 46.8385324, 43.7179324, 40.8975662, 38.3101336, 35.8883348, 33.5648695, 31.2724379, 28.7494751, 26.2006517, 23.6363737, 21.0670471, 18.525831, 15.9912767, 13.4646886, 10.9473713, 8.38916945, 5.86343143, 3.39204573, 0.99690084, -1.26416265, -3.41958917, -5.46187108, -7.3835007, -9.22601578, -10.9132451, -12.4180629, -13.7133433, -14.8183369, -15.6409909, -16.1356292, -16.2565754, -15.8913591, -15.0878157, -13.8269866, -12.089913, -9.575889, -6.66040207, -3.43719234, 0., 3.32808289, 6.77440878, 10.3369783, 14.0137919, 17.7915291, 21.6840401, 25.6938539, 29.8234995, 34.4683192, 39.0809036, 43.5066565, 47.5909816, 52.0892308, 55.5728806, 57.5233557, 57.4220807, 57.0603078, 52.6857029, 42.8557596, 26.1279713, 0.306807571, -33.7911656, -77.5102898, -128.97305, -194.828091, -250.841278, -280.499173, -267.288336}, 
{45.9632263, 43.5948047, 41.2944117, 39.0497927, 36.8486933, 34.6788586, 32.5280343, 30.3839656, 28.2343979, 26.0514214, 23.8446989, 21.6082379, 19.3360459, 16.9792869, 14.5919492, 12.1851779, 9.7701179, 7.3618403, 4.9659933, 2.59215129, 0.249888704, -2.11170034, -4.39836902, -6.57635081, -8.61187918, -10.3777734, -11.9710468, -13.3952985, -14.6541278, -15.978402, -17.0535448, -17.792248, -18.1072036, -17.711688, -16.7975747, -15.3573219, -13.3833878, -10.6397615, -7.43875759, -3.86422188, 0., 3.91966746, 8.02148764, 12.2817728, 16.6768352, 21.1786397, 25.7695849, 30.427722, 35.1311022, 39.788373, 44.4747508, 49.1960481, 53.9580778, 59.5390344, 64.8633958, 69.6280218, 73.5297724, 77.145323, 78.9397916, 78.2581119, 74.4452177, 65.9571768, 54.8055203, 38.1237308, 26.1576508, 25.2296282, 44.2436292, 91.4582154, 175.131949}, 
{51.3870855, 48.6794894, 46.1831565, 43.8403542, 41.59335, 39.3844114, 37.155806, 34.8498012, 32.4086645, 29.4555478, 26.3794804, 23.2503763, 20.1381492, 17.3421612, 14.6110986, 11.9230958, 9.25628726, 6.47079096, 3.70996456, 0.999149184, -1.63631406, -4.15046897, -6.54683554, -8.80831871, -10.9178234, -12.8574803, -14.6112782, -16.1624317, -17.4941555, -18.7345028, -19.6639141, -20.2076685, -20.2910453, -19.6603376, -18.491405, -16.7811208, -14.5263586, -11.5733769, -8.12991022, -4.25307808, 0., 4.55743405, 9.3827832, 14.4248361, 19.6323815, 24.7770149, 30.0555953, 35.4877886, 41.0932604, 46.9053614, 52.9245989, 59.1651649, 65.6412515, 72.5154609, 79.5942109, 86.8323292, 94.184644, 100.942288, 107.989262, 115.545873, 123.832425, 131.321797, 143.476579, 160.747982, 191.051251, 237.852606, 314.950336, 433.559716, 604.896017}, 
{52.5791609, 52.189412, 50.9809292, 49.0755733, 46.5952055, 43.6616866, 40.3968776, 36.9226396, 33.3608334, 30.1281937, 26.9337584, 23.7814391, 20.6751474, 17.664573, 14.689538, 11.7356427, 8.78848729, 5.69032192, 2.6274371, -0.357226835, -3.22072955, -5.77288536, -8.1768974, -10.4487235, -12.6043213, -14.8262476, -16.8972218, -18.766562, -20.3835865, -21.7939441, -22.8120902, -23.348811, -23.3148923, -22.3006067, -20.6654593, -18.4484415, -15.6885449, -12.4879459, -8.79717705, -4.62995592, 0., 5.28008011, 10.9150172, 16.810651, 22.8728212, 28.6611796, 34.566229, 40.6322844, 46.9036604, 53.1753873, 59.8407787, 67.0438633, 74.92867, 83.7004277, 93.417485, 104.199391, 116.165694, 127.53876, 141.094195, 157.710418, 178.265852, 202.213911, 234.708038, 276.259409, 335.688216, 421.713219, 535.1499, 678.789558, 855.423494}
};

  const double a4tab[50][69] = {
{-4795.25483, -4657.11547, -4533.57386, -4420.99082, -4315.72713, -4214.14359, -4112.601, -4007.46016, -3895.08186, -3749.86323, -3598.9142, -3447.38103, -3300.40998, -3183.58631, -3073.44169, -2966.94678, -2861.07226, -2743.59299, -2624.35375, -2504.00351, -2383.19127, -2264.35938, -2145.64608, -2026.98299, -1908.30176, -1789.65691, -1670.80798, -1551.63744, -1432.02772, -1310.85768, -1189.41483, -1067.98309, -946.846341, -827.165104, -707.996037, -589.272406, -470.927472, -353.153181, -235.520643, -117.859652, 0., 117.602012, 235.993705, 355.595891, 476.829381, 595.994952, 719.281464, 848.757742, 986.492612, 1125.29678, 1280.20043, 1456.97565, 1661.39451, 1865.75367, 2122.69078, 2451.36808, 2870.9478, 3306.89371, 3909.54591, 4735.54605, 5841.53576, 7186.8075, 9120.05042, 11591.8611, 15444.2924, 19758.4141, 29767.7712, 49167.7902, 81653.8973}, 
{-4386.2486, -4276.4998, -4166.87176, -4057.27842, -3947.63375, -3837.85169, -3727.8462, -3617.53124, -3506.82077, -3396.24985, -3284.86287, -3172.32536, -3058.30283, -2934.12142, -2811.12175, -2692.30507, -2580.67265, -2501.29049, -2426.26919, -2349.78411, -2266.0106, -2143.26267, -2011.92156, -1876.50719, -1741.53944, -1625.69861, -1513.68005, -1404.33951, -1296.53273, -1186.17111, -1076.23248, -966.750311, -857.758075, -749.581249, -641.844515, -534.46455, -427.358033, -320.645434, -213.958122, -107.131257, 0., 107.056887, 214.965285, 324.107473, 434.865731, 543.216635, 655.710448, 774.491733, 901.705051, 1031.01986, 1176.44587, 1343.51767, 1537.76987, 1733.97589, 1980.73599, 2295.88922, 2697.27467, 3113.77504, 3687.76832, 4472.67611, 5521.92004, 6796.91537, 8627.10278, 10961.4315, 14608.3709, 18698.1916, 28214.4395, 46677.342, 77607.1262}, 
{-4089.70007, -3999.83272, -3900.40845, -3793.52861, -3681.29454, -3565.8076, -3449.16915, -3333.48054, -3220.84311, -3124.77905, -3031.40055, -2938.24064, -2842.83236, -2728.08787, -2612.00939, -2497.97828, -2389.37592, -2309.8317, -2234.37973, -2158.30216, -2076.88115, -1964.99854, -1846.4969, -1724.8185, -1603.40559, -1496.11382, -1391.80673, -1289.76124, -1189.25426, -1088.05749, -987.55515, -887.626246, -788.149784, -688.774549, -589.701849, -490.902776, -392.348421, -294.336138, -196.38025, -98.321342, 0., 97.9874303, 196.858427, 297.074709, 399.097994, 499.259521, 603.803679, 714.844378, 834.49553, 956.8933, 1095.32044, 1255.08195, 1441.48285, 1629.7946, 1867.36916, 2171.52496, 2559.58043, 2963.77126, 3520.53168, 4281.21321, 5297.16736, 6532.52802, 8300.29955, 10550.9045, 14057.1258, 18028.6165, 27105.1878, 44569.1112, 73702.6582}, 
{-3869.29607, -3792.59815, -3700.71666, -3596.83571, -3484.13943, -3365.81194, -3245.03736, -3124.99981, -3008.88343, -2913.32449, -2822.67409, -2734.73548, -2647.31191, -2552.51261, -2456.11249, -2358.19241, -2258.83326, -2155.5432, -2052.00489, -1949.32827, -1848.62329, -1756.61105, -1666.5459, -1577.29332, -1487.71881, -1392.40141, -1296.20764, -1199.71758, -1103.51131, -1010.11081, -917.377477, -825.114634, -733.125587, -640.567837, -548.148831, -455.930204, -363.973591, -272.876208, -181.949879, -91.0420086, 0., 90.1866618, 181.269486, 273.857901, 368.561335, 462.482307, 561.139919, 666.546362, 780.71383, 897.948208, 1031.05052, 1185.11547, 1365.23778, 1545.93007, 1775.10198, 2070.08108, 2448.19491, 2844.61617, 3391.68922, 4139.60356, 5138.54867, 6355.55344, 8090.28923, 10297.6099, 13712.8828, 17642.831, 26328.2418, 42776.3866, 69994.5366}, 
{-3688.72344, -3620.28003, -3534.3291, -3434.29407, -3323.59833, -3205.66528, -3083.91832, -2961.78084, -2842.67626, -2739.75984, -2642.83036, -2551.41848, -2465.05485, -2394.42261, -2323.43893, -2247.17347, -2160.69588, -2024.75159, -1882.46418, -1742.633, -1614.05739, -1545.14423, -1489.24234, -1439.30806, -1388.29774, -1306.06025, -1217.90239, -1126.02348, -1032.62284, -945.925062, -859.694111, -773.719216, -687.789602, -600.783951, -513.766253, -426.889953, -340.308496, -254.916558, -169.829861, -84.9053578, 0., 83.447602, 167.794816, 253.817352, 342.290922, 431.243692, 525.297932, 626.328371, 736.209734, 849.215705, 977.862472, 1127.06518, 1301.73898, 1475.10257, 1696.44611, 1983.36334, 2353.44796, 2744.04363, 3285.09418, 4026.29337, 5017.33496, 6227.8996, 7947.72035, 10138.8771, 13497.9681, 17433.977, 25771.8272, 41232.4571, 66536.8049}, 
{-3511.669, -3448.36228, -3367.77851, -3272.99804, -3167.10117, -3053.16822, -2934.27953, -2813.5154, -2693.95617, -2581.95877, -2476.01625, -2377.89827, -2289.37451, -2240.84482, -2193.9966, -2139.14746, -2066.61499, -1903.78345, -1729.07712, -1557.98694, -1406.00384, -1357.64213, -1331.76003, -1316.23912, -1298.96098, -1228.58923, -1147.91061, -1060.49389, -969.90785, -889.094259, -808.499694, -727.943732, -647.245948, -565.245724, -483.134905, -401.125143, -319.42809, -239.108105, -159.183049, -79.523491, 0., 77.5632713, 156.030772, 236.313367, 319.321921, 403.902372, 493.856483, 590.921089, 696.833024, 805.726906, 929.982676, 1074.37806, 1243.69077, 1410.03238, 1623.91323, 1903.17752, 2265.66943, 2649.78746, 3184.59978, 3919.72886, 4904.79718, 6111.47446, 7823.24142, 10012.0358, 13334.7074, 17295.1967, 25324.1697, 39870.6114, 63383.5064}, 
{-3301.81958, -3242.32883, -3167.59761, -3080.04195, -2982.07786, -2876.12137, -2764.5885, -2649.89526, -2534.45768, -2417.79493, -2306.37861, -2203.78347, -2113.58426, -2078.80618, -2047.79337, -2008.34038, -1948.2418, -1778.96537, -1595.16323, -1415.16067, -1257.28299, -1221.1488, -1211.27277, -1213.4629, -1213.52717, -1151.48722, -1077.25191, -994.943767, -908.685307, -833.212394, -757.788871, -682.291925, -606.598741, -529.775993, -452.835581, -375.978895, -299.407326, -224.101762, -149.172293, -74.5085093, 0., 72.3266902, 145.573709, 220.706247, 298.689499, 378.817046, 464.394334, 557.055201, 658.433482, 762.51293, 881.637499, 1020.50106, 1183.79748, 1343.43977, 1550.015, 1821.32939, 2175.18916, 2549.58149, 3074.05924, 3798.35623, 4772.2063, 5968.186, 7667.50096, 9854.41554, 13145.4268, 17119.632, 24873.4951, 38624.1385, 60588.6845}, 
{-3022.862, -2967.66361, -2900.31911, -2822.52015, -2735.95834, -2642.32532, -2543.31273, -2440.61219, -2335.91534, -2225.14199, -2118.06432, -2018.68268, -1930.99743, -1895.33367, -1864.83711, -1828.97824, -1777.2275, -1636.62394, -1484.04202, -1333.92478, -1200.71524, -1162.70827, -1144.95436, -1136.35581, -1125.81492, -1066.2531, -996.945919, -921.188061, -842.274195, -771.873467, -701.556286, -631.267537, -560.952101, -490.197592, -419.449071, -348.794329, -278.321158, -208.548445, -138.960447, -69.472514, 0., 67.5308789, 136.01998, 206.356297, 279.428824, 354.34641, 434.49025, 521.461393, 616.860889, 714.604893, 827.05331, 958.88115, 1114.76342, 1268.04501, 1467.2631, 1729.62474, 2072.337, 2431.15956, 2937.32579, 3640.6217, 4590.83328, 5759.94219, 7431.14749, 9603.34594, 12852.4522, 16800.4251, 24308.029, 37426.3271, 58206.3824}, 
{-2638.48309, -2589.85055, -2532.47575, -2467.52697, -2396.17251, -2319.58066, -2238.91972, -2155.35797, -2070.06371, -1981.87362, -1895.22025, -1812.20453, -1734.92739, -1677.45422, -1625.13573, -1575.28704, -1525.22331, -1463.08574, -1399.03301, -1334.04983, -1269.12095, -1209.3646, -1149.97859, -1090.29423, -1029.64285, -964.385738, -898.012251, -831.041724, -763.993491, -698.671476, -633.796583, -569.374308, -505.410146, -442.333357, -379.556167, -316.914565, -254.24454, -191.09907, -127.710361, -64.0276063, 0., 62.9688578, 126.965942, 192.62382, 260.575062, 328.849163, 401.722993, 480.870349, 567.965028, 657.033913, 760.456481, 882.965295, 1029.29292, 1176.56839, 1368.1692, 1619.86934, 1947.44278, 2282.2555, 2758.25267, 3424.97148, 4331.9491, 5448.65099, 7064.82952, 9196.15657, 12378.1095, 16230.7179, 23515.9972, 36210.4661, 56290.6436}, 
{-2405.71268, -2357.37371, -2305.09333, -2248.92728, -2188.93131, -2125.16116, -2057.67258, -1986.52131, -1911.7631, -1825.78366, -1739.37678, -1655.66624, -1577.77578, -1521.53728, -1472.28315, -1428.05392, -1386.8901, -1345.84384, -1304.33938, -1260.8126, -1213.69937, -1153.95915, -1090.49479, -1024.73273, -958.099404, -895.75918, -833.905408, -772.469359, -711.382305, -650.167044, -589.326712, -528.955972, -469.149485, -410.875751, -353.00606, -295.28554, -237.459321, -178.716175, -119.580127, -60.018845, 0., 59.1806414, 119.414746, 181.265883, 245.297621, 309.449323, 377.958449, 452.438249, 534.501976, 618.058134, 715.506621, 831.542589, 970.861189, 1111.62937, 1295.68177, 1538.3248, 1854.86492, 2173.99376, 2632.27845, 3279.67134, 4166.12475, 5260.24742, 6856.02255, 8984.53783, 12142.789, 15991.1783, 23070.5741, 35201.7283, 54205.3928}, 
{-2204.0302, -2154.47867, -2105.90828, -2057.23597, -2007.37868, -1955.25335, -1899.77692, -1839.86634, -1774.43854, -1688.83889, -1600.98454, -1515.22105, -1435.89401, -1383.22182, -1339.32807, -1302.20921, -1269.86166, -1242.18385, -1214.50944, -1184.07407, -1148.11338, -1090.76646, -1027.60414, -961.10068, -893.730354, -835.115058, -777.722379, -721.167533, -665.065735, -607.620119, -550.422817, -493.653877, -437.493349, -383.360258, -329.700085, -276.19729, -222.536329, -167.664731, -112.298656, -56.4173358, 0., 55.7448171, 112.55891, 170.95477, 231.444888, 291.989405, 356.674102, 427.03241, 504.59776, 583.170483, 675.11035, 785.044031, 917.598198, 1051.62411, 1227.83402, 1461.16475, 1766.55314, 2070.16998, 2511.22455, 3140.16011, 4007.4199, 5080.84751, 6658.68514, 8786.8724, 11925.8417, 15777.5862, 22663.3632, 34249.6139, 52202.7796}, 
{-2019.11349, -1967.86498, -1922.03119, -1879.54631, -1838.34453, -1796.36004, -1751.52704, -1701.77973, -1645.05228, -1560.9017, -1472.99026, -1386.60303, -1307.02508, -1257.73185, -1218.54189, -1187.46414, -1162.50754, -1146.0165, -1129.93027, -1110.52355, -1084.07109, -1029.27314, -967.008674, -900.58221, -833.298258, -778.328451, -725.163331, -673.160561, -621.677803, -567.845774, -514.139863, -460.808512, -408.100162, -357.77398, -307.963395, -258.312561, -208.46563, -157.218748, -105.403279, -53.0025778, 0., 52.510405, 106.100252, 161.230464, 218.361962, 275.485215, 336.519778, 402.914754, 476.119245, 549.893005, 636.450223, 740.315742, 866.014403, 993.013679, 1160.91873, 1384.27734, 1677.6373, 1964.76951, 2387.3094, 2996.11552, 3842.0464, 4892.15254, 6448.71667, 8571.55319, 11686.9318, 15534.6254, 22233.9388, 33309.2935, 50285.1113}, 
{-1852.906, -1799.92433, -1755.7931, -1717.7438, -1683.00798, -1648.81714, -1612.4028, -1570.99649, -1521.82971, -1440.70191, -1354.84952, -1270.0769, -1192.18842, -1146.35433, -1111.26673, -1084.98361, -1065.56296, -1056.86211, -1048.81999, -1037.17485, -1017.66497, -965.979445, -905.925352, -841.260613, -775.743151, -724.511292, -675.390392, -627.586209, -580.304504, -530.041543, -479.796376, -429.858557, -380.517643, -333.704815, -287.411353, -241.270163, -194.914152, -147.113112, -98.7083065, -49.6778856, 0., 49.3678455, 99.8185964, 151.765844, 205.623179, 259.438588, 316.93751, 379.479777, 448.425223, 517.590152, 598.895338, 696.718027, 815.435464, 935.025994, 1094.02532, 1306.57025, 1586.79759, 1856.40165, 2258.93873, 2845.52263, 3667.26715, 4690.50381, 6220.69328, 8329.94908, 11413.3409, 15245.1024, 21759.8548, 32352.1221, 48416.4288}, 
{-1707.35118, -1653.0484, -1609.52501, -1573.71394, -1542.54814, -1512.96053, -1481.88405, -1446.25164, -1402.99624, -1326.96932, -1246.01786, -1165.90739, -1092.40343, -1050.37625, -1018.84473, -995.932496, -979.763174, -974.240986, -969.396721, -961.041762, -944.987495, -897.385648, -841.571127, -781.219183, -720.005065, -672.775514, -627.56569, -583.582244, -540.031828, -493.404962, -446.710883, -400.242696, -354.293506, -310.740658, -267.659322, -224.708907, -181.548824, -137.082706, -92.0280493, -46.3465739, 0., 46.2075789, 93.4937656, 142.233789, 192.802879, 243.351361, 297.36933, 356.121975, 420.874489, 485.626378, 561.814791, 653.61119, 765.18704, 876.888988, 1026.24324, 1226.95119, 1492.71423, 1743.67573, 2124.5183, 2686.36654, 3480.34503, 4472.24266, 5969.19114, 8053.42895, 11092.3504, 14891.8235, 21218.6652, 31349.4551, 46560.7729}, 
{-1584.39249, -1529.62886, -1485.55795, -1449.34221, -1418.14411, -1389.1261, -1359.45065, -1326.28022, -1286.77728, -1218.43376, -1145.95085, -1074.35921, -1008.68952, -971.084549, -942.618, -921.475684, -905.843417, -897.67345, -889.878586, -879.138066, -862.131131, -819.992009, -771.162956, -718.541216, -665.024034, -622.233049, -580.851352, -540.286431, -499.94577, -457.133564, -414.20191, -371.399611, -328.975471, -288.469407, -248.322665, -208.267601, -168.036574, -126.862414, -85.1768186, -42.9119574, 0., 42.9200457, 86.9055837, 132.30718, 179.4754, 226.72537, 277.257271, 332.235845, 392.825836, 453.366138, 524.577678, 610.355535, 714.594789, 817.830593, 956.661925, 1144.32784, 1394.06738, 1625.20104, 1982.45386, 2516.63232, 3278.54291, 4233.71041, 5688.78638, 7733.36167, 10711.2418, 14457.5951, 20587.924, 30272.6475, 44682.1844}, 
{-1485.97338, -1432.0574, -1386.22294, -1346.5141, -1310.975, -1277.64975, -1244.58246, -1209.81724, -1171.39821, -1113.82504, -1054.10405, -995.697115, -942.066106, -909.766207, -883.928673, -862.778061, -844.53893, -826.679817, -808.483705, -788.47756, -765.188346, -730.298793, -691.917793, -651.310007, -609.740092, -571.99583, -534.409508, -496.836536, -459.132322, -420.424886, -381.587985, -342.767985, -304.111253, -266.47896, -229.016746, -191.585055, -154.044331, -116.18712, -77.9689255, -39.2773507, 0., 39.3956861, 79.8338744, 121.658895, 165.21508, 209.062449, 256.043366, 307.215884, 363.638057, 420.173885, 486.553096, 566.311367, 662.984373, 757.078737, 884.370809, 1057.60788, 1289.53726, 1499.58688, 1831.15114, 2334.30506, 3059.12368, 3971.24838, 5374.05517, 7361.11611, 10257.2967, 13925.2237, 19845.1853, 29093.0545, 42744.7041}, 
{-1376.76382, -1326.68783, -1282.07861, -1241.75639, -1204.5414, -1169.25385, -1134.71398, -1099.742, -1063.15815, -1015.4217, -967.058202, -920.232267, -877.1085, -849.472809, -826.019977, -805.066087, -784.927224, -759.752074, -733.691079, -706.727285, -678.843735, -650.016673, -620.238666, -589.495479, -557.772878, -524.522938, -490.478592, -455.839083, -420.803651, -385.914211, -350.890265, -315.793989, -280.687555, -245.903895, -211.126125, -176.308115, -141.403738, -106.570871, -71.4777782, -35.9967311, 0., 36.0990347, 73.185657, 111.604041, 151.698359, 192.269829, 235.822765, 283.318523, 335.718461, 388.430658, 450.191058, 524.182331, 613.587143, 698.997506, 815.223007, 974.482576, 1188.99514, 1378.73642, 1685.06585, 2157.09963, 2843.95399, 3710.28866, 5058.5893, 6988.85511, 9802.94449, 13373.0307, 19088.1684, 27958.192, 40992.9357}, 
{-1237.58523, -1193.72216, -1155.05249, -1120.43609, -1088.73286, -1058.80268, -1029.50544, -999.701022, -968.249314, -926.239712, -883.41079, -841.730634, -803.167327, -778.318849, -757.071429, -737.941191, -719.444261, -696.334781, -672.39565, -647.647786, -622.112106, -595.968306, -569.015013, -541.209634, -512.509574, -482.320899, -451.372892, -419.843492, -387.910642, -356.060285, -324.03916, -291.902006, -259.703565, -227.751645, -195.74669, -163.642216, -131.391734, -99.0699117, -66.4606474, -33.4689929, 0., 33.4513572, 67.8059178, 103.394599, 140.548317, 178.214416, 218.660816, 262.771863, 311.431906, 360.448557, 417.813589, 486.442042, 569.248956, 647.275805, 754.06062, 901.267868, 1100.56201, 1272.16848, 1555.76639, 1999.59584, 2651.8969, 3474.58878, 4774.87422, 6664.71547, 9415.52636, 12887.0161, 18426.2608, 27008.4947, 39608.9522}, 
{-1105.58134, -1068.28477, -1035.7158, -1006.80631, -980.488224, -955.693435, -931.35385, -906.401371, -879.767902, -843.331603, -805.899619, -769.225353, -735.062207, -712.723629, -693.378956, -675.757571, -658.588856, -637.292406, -615.231308, -592.458861, -569.028363, -545.287037, -520.876688, -495.733044, -469.791833, -442.442738, -414.385952, -385.775624, -356.765902, -327.780257, -298.595783, -269.258899, -239.816024, -210.549023, -181.174687, -151.645255, -121.912965, -91.969331, -61.7116068, -31.0763203, 0., 30.9521313, 62.7237173, 95.6297075, 129.985051, 164.863675, 202.31796, 243.159265, 288.198947, 333.654277, 386.768336, 450.190117, 526.568614, 597.431846, 694.998168, 830.364963, 1014.62961, 1168.42203, 1429.49406, 1845.13006, 2462.6144, 3241.06015, 4492.26561, 6339.07397, 9025.31546, 12395.8862, 17770.4317, 26099.91, 38335.2793}, 
{-981.843359, -951.135893, -924.577433, -901.1915, -880.001617, -860.031308, -840.304093, -819.843497, -797.673041, -766.547528, -734.266689, -702.361533, -672.363072, -652.260502, -634.543371, -618.159416, -602.056372, -582.353316, -561.958105, -540.949937, -519.408011, -497.806608, -475.671811, -452.924788, -429.486703, -404.761057, -379.393749, -353.513012, -327.24708, -300.95003, -274.433913, -247.736624, -220.896058, -194.167337, -167.28424, -140.19777, -112.858933, -85.1877245, -57.178561, -28.7948502, 0., 28.5987857, 57.9394136, 88.3159941, 120.022638, 152.2321, 186.808389, 224.494156, 266.032056, 308.046801, 357.046158, 415.419955, 485.558019, 549.519157, 638.156626, 761.992664, 931.549505, 1067.95242, 1306.87941, 1694.61147, 2277.42964, 3011.70205, 4213.44837, 6014.88302, 8636.06863, 11907.0428, 17130.0951, 25239.4523, 37169.3408}, 
{-867.462487, -843.035749, -822.146291, -803.916104, -787.467177, -771.921499, -756.401059, -740.027848, -721.923855, -695.737645, -668.254002, -640.784284, -614.639853, -596.502818, -580.165486, -564.790917, -549.542169, -531.245875, -512.336089, -492.910439, -473.066553, -453.360756, -433.248503, -412.643943, -391.461228, -369.148461, -346.272259, -322.93319, -299.231825, -275.445509, -251.427325, -227.207134, -202.814794, -178.477895, -153.949475, -129.180301, -104.121142, -78.6436884, -52.8094149, -26.6007194, 0., 26.3887492, 53.4533645, 81.4600863, 110.675155, 140.334184, 172.146291, 206.789968, 244.943705, 283.625111, 328.637912, 382.124952, 446.229074, 503.591266, 583.656969, 696.36977, 851.673255, 971.215005, 1188.55301, 1548.94926, 2097.66575, 2788.51376, 3941.10738, 5695.09497, 8251.54269, 11427.8877, 16514.6653, 24434.1361, 36108.5609}, 
{-763.529923, -744.744569, -728.931269, -715.304576, -703.079038, -691.469208, -679.689636, -666.954874, -652.479471, -630.752109, -607.603558, -584.138715, -561.462481, -545.02393, -529.846115, -515.296264, -500.741608, -483.698449, -466.12531, -448.129792, -429.819492, -411.783222, -393.454883, -374.749589, -355.582454, -335.477555, -314.897458, -293.913693, -272.597788, -251.142599, -229.449798, -207.542384, -185.443357, -163.352002, -141.044517, -118.473386, -95.5910942, -72.2558187, -48.5520731, -24.4700646, 0., 24.3194502, 49.2659279, 75.0686116, 101.95668, 129.18442, 158.345858, 190.060128, 224.946366, 260.38819, 301.534457, 350.298506, 408.59368, 459.701698, 531.62017, 633.715083, 775.352427, 878.665118, 1075.14545, 1409.05263, 1924.64589, 2573.49455, 3677.92749, 5382.66222, 7875.49447, 10965.8228, 15933.5562, 23690.976, 35150.3634}, 
{-671.136869, -657.022582, -645.441263, -635.681364, -627.031337, -618.779636, -610.214713, -600.625022, -589.299015, -571.441079, -552.057359, -532.069936, -512.400888, -497.39719, -483.186069, -469.319648, -455.35005, -439.439403, -423.08582, -406.397423, -389.482331, -372.907744, -356.139071, -339.100802, -321.717425, -303.620942, -285.145324, -266.332054, -247.222617, -227.917203, -208.375107, -188.61433, -168.652875, -148.660967, -128.443493, -107.957565, -87.1602927, -65.9427114, -44.3544406, -22.3790227, 0., 22.3883175, 45.3774618, 69.1481977, 93.8812897, 118.797302, 145.421279, 174.318067, 206.052511, 238.33502, 275.726649, 319.934016, 372.66374, 417.903983, 482.167203, 574.247404, 702.938588, 790.758109, 967.287275, 1275.83075, 1759.69318, 2368.6437, 3426.5936, 5080.53715, 7511.68081, 10528.25, 15396.182, 23016.9864, 34292.1725}, 
{-591.374525, -580.630018, -572.185168, -565.370918, -559.518208, -553.957981, -548.021179, -541.038742, -532.341613, -517.654709, -501.357407, -484.223055, -467.025006, -453.195949, -439.786162, -426.505259, -413.062857, -398.197102, -382.977668, -367.502757, -351.870573, -336.568062, -321.149189, -305.556662, -289.733187, -273.451226, -256.891831, -240.065808, -222.983963, -205.645227, -188.07703, -170.294927, -152.314472, -134.276096, -116.020528, -97.5133744, -78.7202402, -59.6229627, -40.164422, -20.3037303, 0., 20.5927796, 41.7883242, 63.7054721, 86.463062, 109.187321, 133.386744, 159.577213, 188.274612, 217.464584, 251.205346, 291.024877, 338.451155, 378.251645, 435.419044, 518.185534, 634.783301, 707.949324, 865.609071, 1150.19281, 1604.13079, 2175.96048, 3189.79057, 4791.67214, 7163.85853, 10122.571, 14911.9569, 22419.182, 33531.4121}, 
{-531.224895, -520.626258, -512.797608, -506.999762, -502.493538, -498.539755, -494.399231, -489.332783, -482.601229, -469.841541, -455.387921, -439.950727, -424.240315, -411.213282, -398.435247, -385.718072, -372.873616, -359.085644, -345.04535, -330.815832, -316.46019, -302.247325, -287.952211, -273.555623, -259.038339, -244.396642, -229.589599, -214.591784, -199.377772, -183.808107, -168.017004, -152.02465, -135.851231, -119.609428, -103.189935, -86.5759385, -69.7506277, -52.91163, -35.741918, -18.1389039, 0., 18.9692432, 38.6028073, 58.9265352, 79.9662703, 100.621877, 122.495568, 146.063579, 171.802143, 197.878892, 228.002105, 263.571461, 305.986636, 341.038107, 391.978429, 466.450961, 572.099059, 631.266927, 771.016734, 1033.1115, 1459.31424, 1998.35806, 2971.09571, 4515.64123, 6831.24812, 9763.79368, 14511.5205, 21934.5798, 32893.1229}, 
{-486.583679, -473.979814, -465.030685, -458.861989, -454.599418, -451.368669, -448.295435, -444.505412, -439.124294, -427.323975, -413.765471, -399.155998, -384.202771, -371.710721, -359.450261, -347.28952, -335.096629, -322.375158, -309.503618, -296.495961, -283.36614, -270.088482, -256.732415, -243.32774, -229.904261, -216.721713, -203.487992, -190.140926, -176.618344, -162.629655, -148.432476, -134.056005, -119.529439, -104.937331, -90.2313794, -75.4186386, -60.5061629, -45.9950208, -31.2006471, -15.9324906, 0., 17.4807597, 35.7296022, 54.6597253, 74.1843266, 92.8924074, 112.55104, 133.603101, 156.491467, 179.486059, 206.071889, 237.561016, 275.265497, 306.105316, 351.541434, 418.642738, 514.478116, 560.661054, 683.898001, 925.440002, 1326.53811, 1837.0297, 2772.40681, 4258.3764, 6521.60732, 9452.61168, 14185.5666, 21543.733, 32350.3718}, 
{-448.814905, -434.351947, -424.23195, -417.480755, -413.1242, -410.188124, -407.698365, -404.680764, -400.161159, -388.964091, -375.997217, -361.966893, -347.579478, -335.550049, -323.772757, -312.150471, -300.586061, -288.808322, -276.963829, -265.025081, -252.964578, -240.504183, -227.967286, -215.426643, -202.955006, -191.03281, -179.162057, -167.252429, -155.213609, -142.637141, -129.878101, -116.973426, -103.960054, -90.8894741, -77.7782523, -64.6575056, -51.5583513, -39.2929482, -26.798955, -13.7950721, 0., 16.0603954, 32.9974454, 50.6153159, 68.7181726, 85.5958004, 103.172498, 121.858185, 142.062778, 162.114943, 185.338352, 212.975426, 246.268582, 273.110677, 313.433519, 373.819356, 460.850433, 495.64119, 604.4288, 827.982629, 1207.07205, 1692.46558, 2594.93511, 4028.40848, 6246.18361, 9183.8677, 13908.4607, 21204.1145, 31854.9808}, 
{-415.107526, -399.828275, -389.110135, -381.960008, -377.384796, -374.391399, -371.986718, -369.177655, -364.971111, -354.165601, -341.659766, -328.143861, -314.308141, -302.780952, -291.539221, -280.497968, -269.57221, -258.587577, -247.584233, -236.512952, -225.324507, -213.595855, -201.801116, -190.040593, -178.414589, -167.532131, -156.781306, -146.058926, -135.261805, -123.926816, -112.454684, -100.886198, -89.2621454, -77.5906511, -65.9582303, -54.4187358, -43.0260202, -32.8911881, -22.5879395, -11.7472262, 0., 14.6788678, 30.3399124, 46.6897497, 63.4349954, 78.6160763, 94.2722731, 110.776677, 128.50238, 145.787851, 165.854652, 189.889726, 219.080012, 241.990545, 277.478935, 331.780884, 411.132096, 436.62174, 533.690667, 742.633193, 1103.74364, 1567.96956, 2443.64553, 3835.9396, 6019.38928, 8970.19375, 13688.5361, 20918.6033, 31404.5825}, 
{-381.404703, -367.810516, -358.064516, -351.312694, -346.701042, -343.375552, -340.482213, -337.167018, -332.575957, -322.140584, -310.207104, -297.407282, -284.372885, -273.56541, -263.055001, -252.741533, -242.524879, -232.132597, -221.705805, -211.213306, -200.6239, -189.585152, -178.515596, -167.512531, -156.673252, -146.561658, -136.621805, -126.76435, -116.899948, -106.625098, -96.2902787, -85.9318111, -75.5860158, -65.1968036, -54.9298686, -44.858495, -35.0559672, -26.8956431, -18.630704, -9.81440464, 0., 13.2974927, 27.6712859, 42.7528299, 58.1735747, 71.8319533, 85.7856399, 100.359291, 115.877565, 130.624358, 147.781391, 168.489624, 193.890017, 212.684812, 243.429178, 292.239562, 365.232414, 384.336086, 473.530361, 672.606928, 1021.35747, 1469.13545, 2327.04725, 3698.40244, 5866.54264, 8835.96394, 13546.5125, 20704.1969, 31015.0258}, 
{-350.693994, -339.567532, -331.244199, -325.066325, -320.376238, -316.516268, -312.828743, -308.655992, -303.340345, -293.184706, -281.786597, -269.704117, -257.495364, -247.430419, -237.670605, -228.089227, -218.55959, -208.634708, -198.636295, -188.565772, -178.424563, -168.013895, -157.615462, -147.310764, -137.1813, -127.67297, -118.357114, -109.169472, -100.045784, -90.6878242, -81.3588831, -72.0882854, -62.9053565, -53.6877426, -44.6771193, -35.9634834, -27.6368316, -21.2984096, -14.9224656, -7.99449665, 0., 11.9307478, 25.0149408, 38.8249936, 52.933321, 65.1521334, 77.5181313, 90.3078109, 103.797669, 116.194025, 130.671622, 148.335027, 170.288808, 185.113374, 211.447112, 255.404249, 323.099016, 337.713944, 421.467639, 613.647007, 953.538957, 1388.60832, 2233.60823, 3592.34408, 5751.29326, 8739.15611, 13433.0691, 20498.0598, 30599.1557}, 
{-322.965118, -314.722469, -308.053627, -302.529432, -297.72072, -293.198329, -288.533096, -283.29586, -277.057457, -267.132496, -256.250535, -244.884904, -233.50893, -224.173201, -215.142885, -206.260408, -197.368196, -187.807977, -178.123154, -168.356433, -158.550519, -148.713216, -138.936092, -129.275811, -119.789041, -110.754011, -101.917197, -93.246636, -84.7103687, -76.143424, -67.7000552, -59.4015055, -51.2690183, -43.1182337, -35.2582393, -27.7925197, -20.8245593, -16.1400159, -11.4873314, -6.29712101, 0., 10.5704871, 22.3478263, 34.8625739, 47.6452865, 58.4725586, 69.3304936, 80.4512328, 92.0669177, 102.300478, 114.336951, 129.252163, 148.121941, 159.204636, 181.520537, 221.27246, 284.663221, 296.452606, 376.86367, 564.67644, 898.670942, 1324.54377, 2160.32525, 3511.71103, 5663.53469, 8666.0039, 13328.6919, 20271.3306, 30113.6518}, 
{-298.207795, -292.898473, -287.897243, -283.010547, -278.044824, -272.806515, -267.102061, -260.737901, -253.520476, -243.818486, -233.451208, -222.800179, -212.246936, -203.590978, -195.228694, -186.974434, -178.64255, -169.366469, -159.913835, -150.371369, -140.825792, -131.514249, -122.312866, -113.248194, -104.346783, -95.6927269, -87.2320163, -78.9681837, -70.9047623, -63.0203273, -55.353353, -47.9173557, -40.7258518, -33.5430427, -26.7314855, -20.4044225, -14.675096, -11.4609901, -8.34940844, -4.73189658, 0., 9.2085648, 19.6468917, 30.8219035, 42.2405231, 51.6891709, 61.0834733, 70.6185542, 80.4895374, 88.7473442, 98.5889822, 111.067257, 127.234972, 134.887001, 153.637255, 189.841712, 249.85635, 260.249361, 339.079625, 524.618236, 855.136285, 1275.09739, 2104.19508, 3450.4498, 5593.16045, 8602.74098, 13213.8667, 19995.1477, 29515.1939}, 
{-277.692073, -273.266665, -268.661563, -263.791764, -258.572268, -252.918071, -246.744171, -239.965567, -232.497256, -223.318657, -213.654579, -203.794253, -194.026909, -185.767148, -177.728681, -169.75059, -161.671955, -152.674565, -143.517713, -134.303398, -125.133621, -116.351293, -107.721136, -99.2487811, -90.9398631, -82.7556789, -74.7639327, -66.9879926, -59.4512265, -52.2008259, -45.2268059, -38.543005, -32.163262, -25.856573, -19.9795559, -14.6439865, -9.96164035, -7.81861146, -5.84262963, -3.43574307, 0., 7.74471995, 16.6413325, 26.2149221, 35.990573, 43.8290239, 51.5854432, 59.4506533, 67.6154769, 74.1745522, 82.25336, 92.8811966, 107.087358, 112.157665, 128.36228, 162.227891, 220.281186, 230.228588, 308.945153, 494.485671, 824.904937, 1246.65773, 2072.59887, 3408.55568, 5528.94543, 8536.44835, 13052.1028, 19567.3577, 28573.6617}, 
{-259.939683, -255.968472, -251.491441, -246.505072, -241.005847, -234.990248, -228.454756, -221.395854, -213.810024, -205.1907, -196.23863, -187.151517, -178.127059, -170.170813, -162.349484, -154.537632, -146.609818, -137.847271, -128.955213, -120.045537, -111.230135, -102.937193, -94.8357926, -86.9113071, -79.1491114, -71.3967003, -63.8324798, -56.4969763, -49.4307162, -42.7343354, -36.3642068, -30.3368131, -24.6686367, -19.1508851, -14.1154261, -9.66885215, -5.91775591, -4.66091588, -3.63586439, -2.27231968, 0., 6.27642281, 13.5921656, 21.5074915, 29.5826637, 35.7843028, 41.9037718, 48.1387909, 54.6870805, 59.7199925, 66.272163, 75.3518594, 87.9673494, 90.8661607, 105.021597, 137.146221, 193.952597, 204.451318, 284.137707, 470.805114, 802.246893, 1227.94069, 2050.62697, 3373.0398, 5462.06578, 8452.38814, 12843.3511, 19043.1335, 27459.9144}, 
{-244.613452, -240.742451, -236.178876, -230.976389, -225.188652, -218.869326, -212.072073, -204.850553, -197.258429, -189.207498, -180.950031, -172.596435, -164.257118, -156.526127, -148.836772, -141.106004, -133.250774, -124.687935, -116.034573, -107.407678, -98.9242389, -91.0767404, -83.4564783, -76.0302437, -68.7648278, -61.4070679, -54.2316906, -47.2934683, -40.6471739, -34.4312784, -28.5833766, -23.1247619, -18.0767276, -13.2663693, -8.98685713, -5.33716324, -2.41625999, -1.89649508, -1.67411532, -1.2187429, 0., 4.81904605, 10.5395963, 16.7694066, 23.1162328, 27.6727934, 32.1678966, 36.8153135, 41.8288147, 45.4844063, 50.7087301, 58.4906628, 69.8190813, 70.9233303, 83.4556316, 114.308675, 170.37515, 182.25048, 263.643527, 451.965888, 784.629157, 1214.99982, 2032.62481, 3336.30428, 5382.70068, 8336.57764, 12572.4318, 18412.3265, 26178.3252}, 
{-231.376208, -227.327158, -222.515866, -217.031636, -210.963774, -204.401586, -197.434378, -190.151455, -182.642124, -175.141937, -167.535452, -159.853477, -152.126817, -144.557247, -136.936217, -129.226146, -121.389454, -112.999904, -104.564031, -96.1997155, -88.0248379, -80.5747267, -73.3828348, -66.400063, -59.5773124, -52.5780585, -45.755598, -39.1758021, -32.9045422, -27.1020775, -21.7021359, -16.7328335, -12.2222863, -8.04341573, -4.44161024, -1.50706368, 0.670030087, 0.566059177, 0.0976149603, -0.252129215, 0., 3.38796234, 7.52383021, 12.0704624, 16.6907177, 19.6122814, 22.5072553, 25.6125674, 29.1651461, 31.5685348, 35.6264003, 42.3090248, 52.5866904, 52.2400158, 63.5048119, 93.4272263, 149.053406, 162.959, 246.448854, 436.357316, 769.518734, 1203.88868, 2012.93783, 3290.75124, 5281.02929, 8175.03416, 12224.1654, 17664.7881, 24733.2674}, 
{-219.890779, -215.461152, -210.294408, -204.496731, -198.174303, -191.433308, -184.37993, -177.120352, -169.760758, -162.766899, -155.741564, -148.647109, -141.445889, -133.988329, -126.393489, -118.6685, -110.820491, -102.586526, -94.3518247, -86.2315432, -78.3408367, -71.2359444, -64.414504, -57.8152372, -51.3768655, -44.7009486, -38.1982351, -31.9423113, -26.0067637, -20.5571553, -15.5383054, -10.9870099, -6.94006466, -3.32241444, -0.327446578, 1.96330266, 3.46829701, 2.81815514, 1.73432386, 0.650404896, 0., 1.99854432, 4.58507283, 7.48045392, 10.405556, 11.7205528, 13.0512852, 14.6628993, 16.8205414, 18.073119, 21.0885124, 26.8183632, 36.2143128, 34.727059, 45.009565, 74.2138497, 129.491933, 145.909806, 231.539927, 422.368724, 754.382627, 1190.66081, 1985.91147, 3228.78282, 5147.23075, 7953.77503, 11783.3722, 16790.3694, 23129.1142}, 
{-226.724899, -212.764992, -200.613432, -190.033165, -180.787138, -172.638298, -165.34959, -158.683961, -152.404358, -145.81283, -139.317579, -132.865911, -126.40513, -119.900384, -113.273999, -106.466146, -99.4169927, -91.5262098, -83.4906648, -75.4667261, -67.6107623, -60.4573783, -53.6334118, -47.1439369, -40.9940276, -35.2433287, -29.8205152, -24.7088327, -19.8915271, -15.2693682, -10.9410679, -6.92286233, -3.23098758, 0.61081133, 3.89681083, 6.41377843, 7.9484816, 6.40498557, 4.20584103, 1.89089637, 0., 1.36847131, 3.32249934, 5.4837441, 7.4738656, 7.10299269, 6.52892902, 6.09794707, 6.15631933, 4.78894311, 5.50801613, 9.56436095, 18.2088001, 19.6826667, 33.4500685, 65.965624, 123.683951, 138.560516, 221.348751, 408.302935, 735.677349, 1179.71831, 1956.70398, 3138.69873, 4942.56472, 7595.34243, 11138.5612, 15660.1283, 21247.9508}, 
{-236.84171, -211.295632, -190.740811, -174.495373, -161.877439, -152.205134, -144.796581, -138.969903, -134.043225, -127.670992, -121.500477, -115.515274, -109.698978, -104.443078, -99.1601186, -93.6705368, -87.7947709, -80.3881138, -72.6222059, -64.7035429, -56.8386203, -49.4919717, -42.5088397, -35.9925048, -30.0462478, -25.3821555, -21.2511794, -17.5130774, -14.0276074, -10.327812, -6.73085042, -3.2271669, 0.192794552, 4.54191817, 8.42510028, 11.4505655, 13.2265385, 10.5952775, 7.03736001, 3.26739713, 0., 1.10351971, 2.87733143, 4.77455029, 6.24829142, 4.25773493, 1.74750499, -0.831709239, -3.02921862, -7.34619724, -9.19934742, -6.95723472, 1.01157529, 6.82986723, 25.4411853, 62.2804239, 122.782477, 135.652615, 213.747207, 393.192996, 710.116727, 1156.14026, 1900.90499, 2996.52076, 4644.78844, 7069.78957, 10256.7919, 14250.7667, 19096.6851}, 
{-232.280967, -202.805662, -179.493136, -161.488734, -147.937801, -137.985681, -130.77772, -125.459263, -121.175653, -114.86966, -108.770235, -102.903752, -97.2965875, -92.6305548, -88.0144139, -83.2123638, -77.9886036, -70.9785561, -63.5267075, -55.8487674, -48.1604456, -40.8346679, -33.8670415, -27.4103897, -21.6175357, -17.5081853, -14.0215263, -10.963629, -8.14056397, -4.8928611, -1.67834712, 1.51069153, 4.68196843, 9.02829042, 12.8982405, 15.8254948, 17.3437298, 13.841766, 9.25607751, 4.37828269, 0., 0.430831866, 1.53721913, 2.70358662, 3.31435919, -0.0847916974, -4.13518649, -8.31689902, -12.1100031, -18.0667268, -21.3661281, -20.2594189, -12.9978114, -5.65089228, 14.4758508, 52.2585556, 112.57336, 127.359326, 203.604497, 375.359842, 676.676327, 1108.31316, 1804.19659, 2804.49135, 4284.70754, 6449.15256, 9251.59936, 12708.4545, 16836.1246}, 
{-217.723575, -189.98152, -168.069873, -151.165799, -138.446465, -129.089038, -122.270683, -117.168567, -112.959856, -106.589909, -100.360423, -94.3412868, -88.602391, -83.9929754, -79.4918381, -74.8571281, -69.8469946, -63.1281994, -55.9868336, -48.6176009, -41.2152054, -34.0498118, -27.2104788, -20.8617257, -15.168072, -11.1893742, -7.83668007, -4.91637443, -2.23484203, 0.921668114, 4.01858121, 7.04145819, 9.97586002, 13.921929, 17.3048121, 19.6642379, 20.5399346, 16.3430655, 10.9933504, 5.28194367, 0., -0.573489477, -0.530341535, -0.4745363, -1.0100539, -5.63864982, -10.9074187, -16.2612305, -21.1449552, -27.7808742, -31.7254815, -31.3126825, -24.8763827, -18.1852808, 0.835428433, 36.825757, 94.4257167, 113.608134, 189.54708, 352.749442, 633.722106, 1037.49814, 1671.00588, 2568.45021, 3871.51835, 5752.88672, 8153.79254, 11073.3328, 14510.6046}, 
{-197.850438, -175.509646, -157.670492, -143.67912, -132.881674, -124.624299, -118.253141, -113.114343, -108.55405, -102.012813, -95.5046078, -89.1378169, -83.0208222, -78.0605024, -73.2473444, -68.3703316, -63.2184475, -56.6677061, -49.7852481, -42.7252446, -35.6418668, -28.701748, -22.0416126, -15.8106468, -10.1580371, -5.99367853, -2.40176533, 0.772799518, 3.68511306, 7.00195942, 10.1620738, 13.1158783, 15.8137953, 19.0748347, 21.6333959, 23.0924661, 23.0550322, 18.2977904, 12.3805353, 6.03677068, 0., -1.83334156, -3.15785451, -4.50520788, -6.40707065, -12.1179021, -18.3574648, -24.5683118, -30.1929959, -36.8968677, -41.0105633, -41.0875167, -35.6811617, -31.1989674, -15.1987182, 16.907766, 69.7086655, 94.3265236, 170.201412, 323.307766, 579.620021, 944.956329, 1505.75997, 2294.23711, 3414.41715, 5000.44737, 6994.1804, 9385.54253, 12164.4601}, 
{-177.34246, -162.076478, -149.494462, -139.181246, -130.721668, -123.700563, -117.702768, -112.313118, -107.11645, -100.319447, -93.4363589, -86.6032825, -79.9563146, -74.3632984, -68.935886, -63.5174757, -57.9514658, -51.4277386, -44.704615, -37.8868996, -31.0793971, -24.3548212, -17.8629043, -11.7212873, -6.04761157, -1.48905451, 2.57809353, 6.24800591, 9.61485598, 13.2341966, 16.5542697, 19.4846971, 21.9351002, 24.3390085, 25.8725726, 26.2358508, 25.1289017, 19.9045551, 13.5489893, 6.70115434, 0., -3.27262161, -6.17782377, -9.13381751, -12.5588139, -19.2366111, -26.2735981, -33.1417511, -39.3130462, -45.8229357, -49.9545293, -50.5544126, -46.4691713, -45.1176209, -33.3452256, -6.56967931, 39.791324, 69.4419812, 144.193951, 284.980784, 512.73603, 831.948853, 1312.88596, 1987.69175, 2922.60026, 4211.28985, 5803.57192, 7685.22462, 9842.0261}, 
{-160.880546, -152.368455, -144.741251, -137.82473, -131.444689, -125.426926, -119.597237, -113.78142, -107.805271, -100.690886, -93.3892444, -86.0476237, -78.8133015, -72.4315259, -66.2124162, -60.0640621, -53.8945533, -47.2389595, -40.5275983, -33.8177671, -27.1667634, -20.5733762, -14.1768151, -8.05778129, -2.29697582, 2.75654163, 7.39777206, 11.6533579, 15.5499414, 19.5045633, 22.9973084, 25.8986598, 28.0791008, 29.5664513, 30.0109229, 29.2200637, 27.0014222, 21.3619739, 14.6300691, 7.33348529, 0., -4.81522689, -9.42275327, -14.1057546, -19.1474064, -26.7088395, -34.4440918, -41.8851566, -48.5640273, -54.9673069, -59.290535, -60.683861, -58.2974344, -60.3669104, -53.3227302, -32.6808411, 6.04281005, 38.8819919, 110.151155, 235.714466, 431.436092, 699.736845, 1096.81095, 1654.65389, 2405.26399, 3404.86948, 4612.77606, 6012.52008, 7587.63786}, 
{-117.464876, -123.918365, -127.492752, -128.518553, -127.326281, -124.246452, -119.609581, -113.746181, -106.986767, -100.040866, -92.7083752, -85.1682045, -77.5992633, -70.6372629, -63.8215906, -57.1484355, -50.6139865, -44.2722108, -38.038408, -31.8856556, -25.7870314, -19.5202853, -13.3319539, -7.27324596, -1.3953703, 4.18790751, 9.51295815, 14.5535955, 19.2836336, 24.2243411, 28.5830951, 32.1147276, 34.5740706, 35.1779179, 34.4343549, 32.3134286, 28.7851863, 22.7771943, 15.718973, 7.99756168, 0., -6.61155838, -13.2352339, -19.9940327, -27.010961, -35.9903708, -44.8413841, -53.054469, -60.1200932, -64.6781795, -67.4099591, -68.1461179, -66.7173418, -74.5060328, -75.1704745, -63.9206668, -35.9666093, 0.717528466, 68.7915847, 178.151228, 338.692126, 545.47663, 852.90036, 1304.26458, 1880.44809, 2572.04849, 3393.21069, 4352.19278, 5457.25286}, 
{-136.527966, -139.159539, -139.607421, -138.121857, -134.953088, -130.351358, -124.566911, -117.849989, -110.450837, -102.93602, -95.112929, -87.1052792, -79.0367843, -71.261443, -63.5805711, -56.0257688, -48.6286365, -41.5501343, -34.6407587, -27.880366, -21.2488126, -14.5876696, -8.07039264, -1.7321522, 4.3918813, 10.3326244, 15.962385, 21.2195577, 26.0425374, 30.7830152, 34.8007711, 37.8688812, 39.7604217, 39.674461, 38.1876865, 35.3027776, 31.0224139, 24.5660901, 17.0329442, 8.73892973, 0., -8.04750415, -16.2361718, -24.5782045, -33.0858038, -42.6581978, -52.0657511, -60.9658549, -69.0159002, -75.50637, -80.6083267, -84.1259248, -85.8633187, -94.7693226, -97.845567, -91.2383424, -71.0939391, -40.7141918, 13.7723711, 99.0816769, 221.929653, 380.179382, 607.105325, 934.522007, 1336.34321, 1800.19873, 2316.75508, 2875.91963, 3467.59973}, 
{-194.664249, -181.230941, -169.302868, -158.648225, -149.035212, -140.232023, -132.006856, -124.127907, -116.363375, -107.812862, -99.1805959, -90.502209, -81.8133354, -73.2086007, -64.6410495, -56.1227188, -47.6656453, -39.1592333, -30.7872051, -22.6106505, -14.6906592, -7.17655082, -0.00589349249, 6.79551487, 13.2018764, 19.3079585, 24.9191718, 29.9614923, 34.3608958, 38.1894679, 41.1686312, 43.1659177, 44.0488595, 43.4935761, 41.6355773, 38.4189602, 33.7878219, 26.7947351, 18.6319311, 9.6001171, 0., -9.24386883, -18.7041645, -28.3297179, -38.0693597, -47.9114983, -57.7495558, -67.516532, -77.1454268, -88.0404396, -98.0748906, -106.593299, -112.940186, -120.019097, -122.191915, -117.379548, -103.502907, -87.2876613, -54.3280535, 0.976913021, 84.2282345, 202.753802, 356.973927, 552.751261, 786.612743, 1083.52895, 1364.03162, 1574.7622, 1662.36213}, 
{-175.36164, -166.624748, -158.159801, -149.925768, -141.881618, -133.986319, -126.198841, -118.478152, -110.783222, -103.030633, -95.2386948, -87.3833331, -79.440472, -71.2174429, -62.9261997, -54.6101037, -46.3125163, -38.121428, -30.017719, -22.0268988, -14.1744769, -6.30376836, 1.30464509, 8.55237609, 15.3410373, 21.2894899, 26.6951987, 31.5728771, 35.9372383, 40.4691856, 44.2507664, 47.0302182, 48.5557784, 48.0057096, 45.9262139, 42.2935185, 37.0838506, 29.5326914, 20.6533127, 10.7182403, 0., -10.8380537, -22.0705555, -33.5813111, -45.2541263, -57.0379543, -68.7253946, -80.1741938, -91.2420989, -101.801961, -111.69038, -120.759063, -128.859713, -138.865711, -146.398417, -150.100866, -148.616093, -144.359406, -130.692659, -104.749978, -63.6654913, -2.1174713, 75.3923952, 179.345261, 277.973065, 352.592929, 374.43045, 317.234105, 154.752373}, 
{-191.788584, -182.390037, -173.740918, -165.639879, -157.885566, -150.276631, -142.611723, -134.68949, -126.308584, -116.144987, -105.569081, -94.8285798, -84.1712001, -74.6575672, -65.3973217, -56.3130147, -47.3271972, -37.9729249, -28.7180425, -19.6408993, -10.8198446, -2.39999017, 5.63378214, 13.2298281, 20.3365033, 26.8791912, 32.8384087, 38.1717007, 42.8366119, 47.2220462, 50.6816452, 53.00041, 53.9633413, 52.8498506, 50.1527642, 45.859319, 39.9567519, 31.9102235, 22.4378773, 11.7357803, 0., -12.5765892, -25.7934503, -39.453239, -53.3586108, -66.8576235, -80.38937, -93.9383452, -107.489044, -121.262378, -134.911858, -148.327413, -161.398972, -174.938597, -187.545228, -198.739941, -208.043809, -213.510736, -216.715838, -217.767058, -216.77234, -208.338342, -209.076867, -219.416934, -264.909434, -355.793203, -538.433422, -850.163683, -1328.31758}, 
{-198.71236, -197.308699, -193.101836, -186.510063, -177.951671, -167.844951, -156.608195, -144.659693, -132.417739, -121.285996, -110.303233, -99.4935919, -88.8812147, -78.707498, -68.6924271, -58.7732421, -48.8871829, -38.459504, -28.1442253, -18.0833812, -8.41900594, 0.229792022, 8.38888204, 16.1070594, 23.4331194, 30.88596, 37.8562326, 44.2046916, 49.7920911, 54.7738453, 58.5981843, 61.0079983, 61.7461775, 59.6423142, 55.7179154, 50.0811906, 42.840349, 34.2000096, 24.1334081, 12.7101898, 0., -14.5670711, -30.0260009, -46.0513217, -62.3175661, -77.5333094, -92.7254243, -107.954826, -123.28243, -138.151262, -153.487283, -169.598565, -186.793178, -205.872273, -226.45361, -248.648032, -272.566377, -292.582224, -316.83858, -347.741193, -387.695808, -434.404278, -504.384028, -598.974796, -746.35905, -974.171705, -1285.39517, -1689.17496, -2194.65663}
};

  const double a5tab[50][69] = {
{4296.48547, 4145.14906, 4028.39347, 3938.25617, 3866.77463, 3805.98633, 3747.92872, 3684.63929, 3608.15549, 3470.8602, 3320.30733, 3164.39619, 3011.02609, 2895.53131, 2787.40219, 2683.56406, 2580.94224, 2467.08826, 2352.05074, 2236.50453, 2121.12447, 2008.42576, 1896.50673, 1785.30605, 1674.76241, 1564.89358, 1455.52754, 1346.57134, 1237.93205, 1128.57981, 1019.73336, 911.67453, 804.685139, 699.870801, 596.360042, 494.10517, 393.058497, 293.389367, 194.746242, 96.9946212, 0., -95.8821865, -191.472354, -287.100981, -383.098545, -476.415091, -572.113704, -671.877035, -777.387735, -883.098866, -1000.81451, -1135.10914, -1290.55726, -1445.5122, -1641.25806, -1892.85779, -2215.37433, -2550.46591, -3015.96209, -3656.28769, -4515.86756, -5562.5241, -7070.48939, -8999.14401, -12015.737, -15367.8837, -23273.8659, -38709.2987, -64649.7973}, 
{3989.50138, 3890.74144, 3791.28858, 3691.26149, 3590.77885, 3489.95935, 3388.92168, 3287.78452, 3186.66657, 3088.65256, 2989.7087, 2888.76729, 2784.76057, 2664.7319, 2544.25806, 2427.02688, 2316.72622, 2240.60956, 2169.37286, 2097.27769, 2018.58566, 1902.04814, 1777.64101, 1649.82996, 1523.08064, 1415.81637, 1312.96214, 1213.40058, 1116.01432, 1016.79356, 918.670313, 821.68419, 725.874783, 631.536102, 538.351566, 446.259007, 355.196258, 265.266634, 176.17629, 87.7968658, 0., -86.9254916, -173.692018, -260.594811, -347.929102, -432.365887, -519.274328, -610.399354, -707.48589, -805.717825, -916.025541, -1042.77838, -1190.34568, -1339.01104, -1526.86385, -1767.90776, -2076.14641, -2395.68728, -2838.38865, -3446.21264, -4261.12134, -5252.34273, -6680.04137, -8501.55957, -11361.8428, -14541.9955, -22090.8829, -36865.43, -61722.5616}, 
{3757.56237, 3688.98233, 3602.74083, 3502.53272, 3392.05283, 3274.99599, 3155.05703, 3035.9308, 2921.31211, 2833.32287, 2749.86003, 2667.24761, 2581.80961, 2470.83397, 2357.29521, 2245.13178, 2138.28212, 2062.17793, 1990.6671, 1919.09078, 1842.79012, 1736.93836, 1625.1117, 1510.71846, 1397.16693, 1298.16665, 1202.70421, 1110.06741, 1019.54408, 928.906343, 839.563964, 751.411028, 664.341619, 578.018704, 492.659931, 408.25183, 324.780933, 242.505048, 161.030918, 80.2365616, 0., -79.1939341, -158.316759, -237.733181, -317.807904, -395.541068, -476.007772, -560.918545, -651.983919, -744.762196, -849.577029, -970.599841, -1112.00206, -1254.41911, -1434.97281, -1667.24897, -1964.8334, -2274.39487, -2703.20307, -3291.61063, -4079.9702, -5039.50056, -6417.95592, -8173.42849, -10924.8947, -14014.4248, -21230.8298, -35251.2358, -58752.7686}, 
{3573.43954, 3518.13048, 3441.64113, 3347.99242, 3241.20528, 3125.30065, 3004.29946, 2882.22265, 2763.09114, 2665.78538, 2573.523, 2484.38113, 2396.4369, 2303.02573, 2208.86315, 2113.92298, 2018.17905, 1918.81705, 1819.71417, 1721.95946, 1626.642, 1540.32011, 1456.42589, 1373.86068, 1291.52585, 1204.16451, 1116.49954, 1029.09559, 942.517315, 859.174238, 777.048176, 695.965823, 615.753875, 535.637321, 456.285245, 377.765026, 300.144041, 223.94418, 148.596504, 73.9865859, 0., -72.5271837, -145.040088, -218.033339, -292.001563, -364.645482, -440.371185, -520.79086, -607.516694, -696.21608, -796.823916, -913.330305, -1049.72535, -1186.00008, -1359.74331, -1584.54477, -1873.9942, -2177.10138, -2597.86799, -3175.71572, -3950.0663, -4894.3482, -6245.96286, -7966.33323, -10645.6449, -13703.476, -20613.072, -33834.2617, -55826.8739}, 
{3409.90395, 3356.44466, 3286.88038, 3203.56312, 3108.8449, 3005.07774, 2894.61363, 2779.80461, 2663.00268, 2546.95431, 2433.45928, 2324.71182, 2222.90615, 2150.49542, 2081.31138, 2009.44471, 1928.98609, 1797.5506, 1660.29474, 1525.89938, 1403.04545, 1339.4171, 1289.09066, 1245.14575, 1200.66196, 1126.02999, 1046.0939, 963.008853, 878.930002, 801.853306, 725.756809, 650.459352, 575.779776, 500.710664, 426.227618, 352.479982, 279.6171, 208.423598, 138.159426, 68.7198163, 0., -66.7649103, -133.555516, -201.012536, -269.776693, -338.38397, -410.421719, -487.372554, -570.71909, -656.063581, -753.121147, -865.726549, -997.714547, -1128.01761, -1295.3337, -1513.4585, -1796.18769, -2094.31954, -2509.84608, -3081.76197, -3849.06182, -4787.23624, -6125.79199, -7831.85622, -10464.8454, -13527.4537, -20156.9751, -32582.0535, -53031.3329}, 
{3239.72667, 3182.18361, 3117.34947, 3045.16737, 2965.58039, 2878.53165, 2783.96423, 2681.82123, 2572.04576, 2437.74389, 2302.43055, 2172.78365, 2055.48109, 2002.43127, 1956.9894, 1907.74119, 1843.27233, 1685.40227, 1516.18949, 1350.92621, 1204.90465, 1161.45301, 1140.61312, 1130.46279, 1119.07984, 1055.98318, 983.233088, 904.330944, 822.778121, 751.19961, 680.323725, 610.002396, 540.087553, 469.557447, 399.487158, 330.078085, 261.53163, 194.782874, 129.006064, 64.1091301, 0., -61.7467835, -123.556551, -186.188021, -250.399908, -315.46138, -384.216524, -458.019878, -538.225981, -620.2888, -713.823669, -822.545353, -950.168616, -1074.73535, -1235.90234, -1447.6535, -1723.97275, -2016.56203, -2426.60002, -2992.98342, -3754.60894, -4688.5153, -6019.17312, -7721.57993, -10323.2482, -13404.6625, -19781.9045, -31462.1566, -50452.6013}, 
{3035.67878, 2973.60609, 2911.93932, 2848.72771, 2782.02045, 2709.86678, 2630.31591, 2541.41704, 2441.2194, 2299.06835, 2153.19849, 2013.14059, 1888.4254, 1848.02151, 1818.24671, 1784.85662, 1733.60687, 1569.39574, 1391.17913, 1217.05561, 1065.12374, 1033.65155, 1028.50034, 1035.70092, 1041.28404, 986.24413, 919.662878, 845.585609, 768.05765, 701.469213, 635.382783, 569.705732, 504.34543, 438.496382, 373.063973, 308.240723, 244.21915, 181.861577, 120.422796, 59.8274049, 0., -57.3124734, -114.736706, -173.077041, -233.137821, -294.582554, -359.81275, -430.089086, -506.672238, -584.87584, -674.286427, -778.543494, -901.286533, -1020.41696, -1175.60758, -1380.79311, -1649.90829, -1934.34159, -2335.59247, -2892.61414, -3644.35983, -4568.53596, -5887.83607, -7587.08678, -10161.6055, -13253.4069, -19407.2257, -30442.1168, -48177.1347}, 
{2770.53136, 2708.97085, 2649.54083, 2590.16668, 2528.77379, 2463.28755, 2391.63335, 2311.73658, 2221.52263, 2091.84191, 1958.52477, 1830.32661, 1716.0028, 1676.45437, 1647.43281, 1616.83524, 1572.5588, 1436.55469, 1289.04433, 1144.30323, 1016.60689, 983.236401, 970.259435, 966.74926, 961.779143, 909.032918, 847.129054, 779.296591, 708.764568, 646.91818, 585.567846, 524.68014, 464.221637, 403.846181, 343.958173, 284.649282, 226.011177, 168.499276, 111.696, 55.5475182, 0., -53.3016495, -106.789491, -161.196848, -217.257045, -274.452337, -335.267548, -400.93643, -472.692737, -545.808803, -629.864368, -728.477751, -845.267274, -959.326113, -1108.60779, -1306.54069, -1566.55318, -1838.17091, -2224.28609, -2763.88818, -3495.96666, -4397.64884, -5693.51063, -7379.95924, -9920.66925, -12991.9915, -18952.3044, -29489.4796, -46291.3886}, 
{2417.05548, 2366.53666, 2309.0449, 2245.40682, 2176.44909, 2102.99835, 2025.88123, 1945.9244, 1863.95449, 1776.9788, 1691.17107, 1608.88567, 1532.477, 1476.91808, 1426.8972, 1379.72126, 1332.69721, 1273.9028, 1213.56578, 1152.68474, 1092.25826, 1037.43127, 983.397481, 929.496943, 875.06971, 816.569606, 757.3774, 697.987635, 638.894853, 581.802574, 525.512772, 470.036397, 415.384401, 361.925557, 309.169865, 256.985149, 205.239231, 153.535543, 102.112055, 50.9423473, 0., -49.5539819, -99.4084146, -150.064689, -202.024195, -253.775575, -308.638068, -367.918164, -432.922351, -499.071792, -575.912436, -667.104903, -776.309814, -885.726459, -1029.06132, -1218.55955, -1466.46631, -1718.56271, -2080.14355, -2590.0396, -3287.08162, -4146.20453, -5397.92663, -7051.77975, -9541.19154, -12538.7207, -18336.5059, -28571.7906, -44881.8187}, 
{2185.86243, 2129.66864, 2072.54176, 2014.18321, 1954.2944, 1892.57674, 1828.73164, 1762.46052, 1693.46479, 1613.81027, 1533.88821, 1456.45424, 1384.26402, 1331.09789, 1284.27693, 1242.1469, 1203.05357, 1164.97587, 1126.77315, 1086.93792, 1043.96269, 989.074843, 930.938071, 870.950935, 810.511996, 754.725473, 699.802008, 645.657899, 592.209444, 538.934767, 486.363612, 434.587547, 383.69814, 334.548323, 286.163756, 238.331466, 190.838476, 142.966997, 95.210792, 47.5588098, 0., -46.4100171, -93.1756315, -140.734562, -189.524529, -237.936684, -289.27394, -344.79264, -405.74913, -467.474595, -539.520605, -625.513568, -729.079894, -833.298458, -970.562217, -1152.71659, -1391.60701, -1630.84359, -1977.80117, -2471.61931, -3151.43754, -3991.40953, -5225.63246, -6875.89452, -9345.32147, -12339.851, -17972.8354, -27755.7736, -43200.1643}, 
{1981.4486, 1917.74747, 1860.59169, 1808.1705, 1758.6731, 1710.28871, 1661.20656, 1609.61586, 1553.70582, 1477.58662, 1399.15813, 1322.2412, 1250.65666, 1201.62972, 1160.21509, 1124.87184, 1094.05904, 1069.0927, 1044.43221, 1017.39389, 985.294068, 932.67365, 874.734515, 813.90314, 752.606002, 700.319407, 649.600069, 600.054533, 551.289343, 501.502803, 452.272994, 403.769755, 356.162928, 310.714537, 266.065364, 221.948375, 178.096535, 133.578698, 89.0575891, 44.5318191, 0., -43.5656426, -87.5295569, -132.282577, -178.215537, -223.723063, -271.99068, -324.207705, -381.563456, -439.300701, -506.933925, -588.031063, -686.16005, -784.990737, -915.94838, -1090.56015, -1320.35322, -1546.9207, -1879.69744, -2358.18424, -3021.88189, -3844.23415, -5062.91291, -6711.52012, -9164.32682, -12162.1323, -17638.9538, -26980.504, -41572.4958}, 
{1794.4152, 1724.40808, 1667.3188, 1620.06974, 1579.58327, 1542.78174, 1506.58753, 1467.923, 1423.71051, 1351.59344, 1275.48475, 1200.01839, 1129.82833, 1084.63602, 1047.95293, 1018.37801, 994.510209, 980.512175, 967.19371, 950.928297, 928.089424, 877.880279, 820.71276, 759.828472, 698.469016, 649.593536, 602.839079, 557.560233, 513.111583, 466.653192, 420.611983, 375.220355, 330.710706, 288.654201, 247.408964, 206.671883, 166.13985, 124.746758, 83.2576926, 41.6747436, 0., -40.8980775, -82.2303199, -124.341187, -167.575137, -210.340302, -255.692002, -304.749227, -358.630969, -412.54525, -475.886418, -552.13785, -644.782925, -738.021909, -862.334537, -1028.91743, -1248.96722, -1462.16565, -1779.83017, -2241.76334, -2887.76773, -3690.53634, -4891.20049, -6534.3606, -8967.55178, -11963.4481, -17291.233, -26220.9625, -40022.6926}, 
{1628.55148, 1554.65838, 1498.11445, 1454.99224, 1421.36432, 1393.30328, 1366.88168, 1338.17208, 1303.24707, 1235.66334, 1163.01568, 1090.38301, 1022.84422, 981.368108, 948.787774, 923.826191, 905.206329, 898.844712, 893.393335, 884.697746, 868.603494, 821.350866, 766.232776, 706.936878, 647.150825, 601.764906, 558.783084, 517.41196, 476.858131, 433.670357, 390.776215, 348.44544, 306.947768, 268.010996, 229.863574, 192.192013, 154.682824, 116.244139, 77.6521993, 38.9048655, 0., -38.3123138, -77.0873241, -116.628058, -157.237541, -197.364371, -239.907775, -285.912554, -336.423505, -386.690844, -445.871786, -517.328963, -604.425008, -691.796217, -809.022092, -966.955801, -1176.45051, -1375.54392, -1677.03085, -2120.89065, -2747.10266, -3527.644, -4706.50075, -6337.96158, -8745.4887, -11730.9564, -16913.4172, -25457.6291, -38528.3499}, 
{1487.64668, 1413.5063, 1358.36997, 1318.04925, 1288.35567, 1265.10079, 1244.09613, 1221.15325, 1192.0837, 1129.6289, 1061.89856, 993.932245, 930.769552, 893.077305, 864.016948, 842.37717, 826.946657, 823.700716, 821.366769, 815.858856, 803.091017, 759.741547, 708.65453, 653.438306, 597.701215, 556.050561, 516.696131, 478.846677, 441.710953, 401.838723, 362.161324, 322.951103, 284.480408, 248.428604, 213.098217, 178.198789, 143.439861, 107.843807, 72.0822061, 36.1394671, 0., -35.7133438, -71.9099732, -108.860855, -146.836957, -184.371235, -224.167871, -267.193037, -314.412902, -361.220084, -416.38373, -483.099433, -564.562788, -645.717902, -755.312447, -903.842611, -1101.80458, -1286.02099, -1570.131, -1994.10022, -2597.89427, -3352.88507, -4504.81926, -6115.86863, -8488.62993, -11451.815, -16489.2504, -24670.9838, -37067.063}, 
{1375.49005, 1305.95973, 1253.47674, 1214.35207, 1184.89673, 1161.42172, 1140.23803, 1117.65665, 1089.9886, 1033.32271, 972.280995, 911.263313, 854.669523, 821.014935, 794.937775, 775.19172, 760.530448, 754.690595, 749.449696, 741.568244, 727.806733, 689.708458, 645.337991, 597.542705, 549.16997, 511.667545, 475.842264, 441.101348, 406.852016, 370.442713, 334.162944, 298.243438, 262.914923, 229.550707, 196.781911, 164.382231, 132.125366, 99.3187258, 66.3888102, 33.2958309, 0., -33.0061597, -66.5076708, -100.757245, -136.007596, -170.936864, -208.00216, -248.086024, -292.070995, -335.615571, -386.915951, -448.944294, -524.672755, -599.191206, -700.507005, -838.745224, -1024.03094, -1192.56234, -1457.96212, -1859.92612, -2438.15015, -3163.58747, -4282.16159, -5861.62737, -8187.4678, -11113.1819, -16002.4765, -23841.5066, -35616.4271}, 
{1295.87085, 1237.02659, 1188.82609, 1149.01198, 1115.32691, 1085.51353, 1057.31449, 1028.47242, 996.729976, 946.577353, 894.310623, 842.97341, 795.609339, 766.432322, 742.847581, 723.430622, 706.756955, 691.424756, 675.977799, 658.982526, 639.005378, 607.907736, 573.643128, 537.46002, 500.606877, 467.832902, 435.48553, 403.412934, 371.463285, 338.766752, 306.176712, 273.828538, 241.857605, 211.020988, 180.583677, 150.432364, 120.453741, 90.4418588, 60.4131083, 30.291239, 0., -30.0957535, -60.6898206, -92.0348943, -124.383668, -156.637225, -190.940513, -228.086866, -268.869621, -309.359907, -356.962153, -414.358576, -484.231398, -551.62037, -643.907167, -770.830997, -942.131067, -1094.13344, -1339.35572, -1716.90239, -2265.87789, -2957.07913, -4034.53329, -5568.78338, -7832.49467, -10702.2149, -15436.8397, -22949.6777, -34154.0375}, 
{1204.65587, 1157.00006, 1114.47906, 1076.03051, 1040.59206, 1007.10134, 974.496, 941.713669, 907.69199, 863.90708, 819.742711, 777.121132, 737.964591, 712.778866, 691.469265, 672.524625, 654.433783, 631.935693, 608.769028, 584.922579, 560.385137, 535.126853, 509.162613, 482.488661, 455.101244, 426.520612, 397.409404, 367.954264, 338.341835, 309.081426, 279.907949, 250.878982, 222.052103, 193.704483, 165.586269, 137.667202, 109.917023, 82.4541578, 55.0401876, 27.5853794, 0., -27.3906496, -55.2573507, -83.8558505, -113.441896, -143.126009, -174.767252, -209.079461, -246.776473, -284.355627, -328.433857, -381.411599, -445.689291, -506.351102, -590.040242, -706.083651, -863.808272, -999.956678, -1225.47393, -1578.72072, -2098.05773, -2753.3297, -3788.44518, -5278.67022, -7480.52598, -10279.4618, -14863.9593, -22096.742, -32840.5336}, 
{1074.63988, 1033.0459, 996.321555, 963.439721, 933.373293, 905.095166, 877.57823, 849.795378, 820.719503, 782.391091, 743.488403, 705.757293, 670.943616, 648.479622, 629.350211, 612.226682, 595.780332, 575.305119, 554.200613, 532.489048, 510.192655, 487.464298, 464.143325, 440.199715, 415.603449, 389.833046, 363.546529, 336.910461, 310.091403, 283.54333, 257.03043, 230.6043, 204.316538, 178.421378, 152.686729, 127.083133, 101.581135, 76.2399296, 50.9059493, 25.5142782, 0., -25.256491, -50.9435458, -77.3042049, -104.581509, -131.993253, -161.217822, -192.908352, -227.717982, -262.453961, -303.153673, -352.008612, -411.210273, -466.199332, -542.618428, -649.359384, -795.314019, -917.313637, -1125.13479, -1456.4935, -1949.1058, -2570.59931, -3568.95528, -5028.68341, -7183.5431, -9910.74615, -14364.4451, -21379.5571, -31790.9996}, 
{952.045299, 916.813429, 886.000482, 858.645986, 833.789466, 810.47045, 787.728463, 764.603033, 740.133686, 707.0761, 673.26719, 640.260019, 609.607655, 589.582905, 572.331196, 556.717696, 541.607574, 522.902498, 503.616539, 483.800265, 463.504248, 443.03075, 422.07797, 400.595797, 378.534124, 355.357107, 331.694667, 307.69099, 283.49026, 259.485015, 235.471748, 211.495303, 187.600525, 164.018334, 140.533068, 117.115138, 93.7349577, 70.3914253, 47.0150752, 23.564927, 0., -23.2517332, -46.8877878, -71.1367257, -96.2271088, -121.468367, -148.375847, -177.545764, -209.574332, -241.585969, -279.037406, -323.913578, -378.199419, -427.713909, -497.074318, -594.73196, -729.138153, -837.32558, -1027.73164, -1337.3751, -1803.27474, -2390.79442, -3351.91758, -4779.57717, -6886.99293, -9541.19849, -13873.6851, -20699.8296, -30835.0089}, 
{837.812697, 808.936927, 783.919108, 761.884766, 741.959432, 723.268633, 704.937899, 686.092758, 665.858739, 637.79203, 608.815237, 580.281624, 553.544456, 535.682276, 520.03296, 505.659659, 491.625529, 474.469091, 456.787981, 438.655203, 420.143763, 401.668489, 382.823831, 363.546065, 343.771465, 322.977172, 301.742249, 280.186624, 258.430227, 236.798139, 215.123074, 193.4429, 171.795483, 150.388204, 129.021612, 107.665771, 86.2907419, 64.8428209, 43.325344, 21.7178806, 0., -21.3729548, -43.0874978, -65.3543693, -88.3843094, -111.556974, -136.246622, -162.996426, -192.349561, -221.745105, -256.071965, -297.114955, -346.658886, -390.927315, -453.490813, -542.358694, -665.540274, -760.3361, -933.747761, -1222.06808, -1661.58988, -2215.46107, -3139.40919, -4533.64928, -6593.803, -9176.58984, -13398.9314, -20062.6987, -29969.7628}, 
{732.882625, 710.050697, 690.4807, 673.391519, 658.002041, 643.53115, 629.197732, 614.220674, 597.81886, 574.368803, 549.868711, 525.47442, 502.341765, 486.371293, 472.076243, 458.714566, 445.544212, 429.746158, 413.486118, 396.852833, 379.935045, 363.21979, 346.238194, 328.919679, 311.19367, 292.577619, 273.577705, 254.288137, 234.803126, 215.37436, 195.875577, 176.337996, 156.792834, 137.423841, 118.04869, 98.6375877, 79.1607395, 59.5282925, 39.7945341, 19.9536936, 0., -19.6167344, -39.5400971, -59.9580924, -81.0587245, -102.2647, -124.83544, -149.265068, -176.047706, -202.924822, -234.244259, -271.6012, -316.590832, -355.872032, -411.950817, -492.396897, -604.779981, -686.688787, -843.66641, -1111.27496, -1525.07654, -2046.14531, -2933.5072, -4293.19749, -6306.90082, -8822.69125, -12947.4358, -19473.3033, -29192.4625}, 
{638.19564, 620.789034, 606.088527, 593.401702, 582.036144, 571.299435, 560.499161, 548.942905, 535.938251, 516.636341, 496.163778, 475.490721, 455.58733, 441.243514, 428.081786, 415.544409, 403.073644, 388.474958, 373.482125, 358.192127, 342.701942, 327.526931, 312.178341, 296.585799, 280.678933, 264.042825, 247.089466, 229.886301, 212.500777, 195.105336, 177.620429, 160.071501, 142.483999, 125.018098, 107.510627, 89.933143, 72.2572021, 54.3820161, 36.3804242, 18.2529207, 0., -17.9796505, -36.2430068, -54.9488516, -74.2559679, -93.5971696, -114.147596, -136.356417, -160.672804, -185.118576, -213.541195, -247.360771, -287.997414, -322.580542, -372.537235, -445.003884, -547.116876, -616.727231, -757.970853, -1005.69828, -1394.76004, -1884.39319, -2736.28873, -4060.51959, -6029.21391, -8485.27376, -12526.4503, -18936.7826, -28500.3094}, 
{554.692301, 541.786237, 531.145858, 522.150773, 514.180592, 506.614925, 498.833382, 490.215574, 480.14111, 464.424568, 447.436604, 429.98284, 412.8689, 399.892497, 387.670329, 375.811181, 363.923841, 350.396751, 336.547181, 322.472055, 308.268298, 294.432189, 280.501557, 266.413583, 252.105452, 237.257168, 222.165962, 206.871889, 191.415003, 175.882728, 160.248799, 144.53432, 128.760398, 113.063829, 97.3037508, 81.4549918, 65.4923811, 49.3381679, 33.0407928, 16.5961166, 0., -16.4582818, -33.193648, -50.3276035, -67.9816533, -85.560007, -104.188383, -124.275204, -146.228893, -168.319821, -193.949684, -224.382126, -260.880791, -291.085326, -335.332971, -400.336967, -492.810558, -550.795024, -677.144355, -906.040578, -1271.66572, -1731.75074, -2549.83088, -3837.91334, -5763.66979, -8170.1084, -12143.2269, -18458.2756, -27890.5049}, 
{483.313163, 473.676603, 466.055961, 459.874189, 454.554236, 449.519054, 444.191593, 437.994804, 430.351638, 417.563406, 403.423355, 388.60309, 373.77422, 361.911802, 350.462611, 339.176876, 327.804822, 315.252799, 302.452463, 289.491591, 276.457961, 263.777843, 251.065126, 238.272192, 225.351423, 212.105025, 198.695625, 185.135674, 171.437623, 157.598192, 143.651858, 129.617363, 115.513453, 101.453888, 87.3243877, 73.1056889, 58.7785281, 44.3309239, 29.7334185, 14.9638358, 0., -15.0492067, -30.3894418, -46.0953047, -62.2413945, -78.1588373, -94.9630951, -113.026156, -132.72001, -152.52201, -175.456634, -202.653723, -235.243119, -261.418867, -300.420927, -358.55346, -442.120629, -489.235757, -601.67018, -813.004394, -1156.81889, -1589.764, -2376.21074, -3627.6765, -5513.19596, -7882.96623, -11805.0174, -18042.9214, -27360.2502}, 
{430.484418, 421.113456, 414.158086, 408.98146, 404.94673, 401.417046, 397.755561, 393.325426, 387.489794, 376.483753, 364.049742, 350.802139, 337.355319, 326.256844, 315.414632, 304.669786, 293.863408, 282.296992, 270.567093, 258.73066, 246.844638, 235.147337, 223.441797, 211.712417, 199.9436, 188.128389, 176.239085, 164.256634, 152.161977, 139.83341, 127.395587, 114.870511, 102.280187, 89.7188522, 77.1073842, 64.4388941, 51.7064931, 39.0573716, 26.2689294, 13.2726458, 0., -13.7764105, -27.9026719, -42.3837523, -57.2246196, -71.5863824, -86.6554121, -102.77422, -120.285319, -137.798375, -158.081883, -182.171493, -211.102857, -233.791427, -268.24113, -320.335691, -395.958842, -432.859693, -532.310436, -727.448645, -1051.41189, -1460.81344, -2218.36381, -3429.83, -5277.77804, -7636.36413, -11536.5847, -17720.2261, -26929.0746}, 
{392.308466, 381.186439, 373.266059, 367.790599, 364.003329, 361.147518, 358.466437, 355.203357, 350.601548, 340.48549, 328.88476, 316.410145, 303.672432, 293.084004, 282.733413, 272.510809, 262.30634, 251.700434, 241.016849, 230.269621, 219.472788, 208.611231, 197.739805, 186.884208, 176.07014, 165.51361, 154.973883, 144.400534, 133.743138, 122.753431, 111.657963, 100.485445, 89.2645885, 78.0710907, 66.8678814, 55.664877, 44.4719937, 33.6607917, 22.7348862, 11.559536, 0., -12.6123, -25.6647331, -39.0785057, -52.774824, -65.6832532, -79.1132972, -93.3828187, -108.809681, -124.072008, -141.783297, -162.917303, -188.447785, -208.078692, -238.561512, -285.377927, -354.009618, -381.598833, -469.302458, -649.937946, -956.32275, -1345.70917, -2077.61012, -3248.67972, -5062.9382, -7430.12351, -11329.7366, -17473.8377, -26574.487}, 
{360.66773, 347.89344, 338.935119, 332.947164, 329.083975, 326.499951, 324.349491, 321.786994, 317.966858, 308.406182, 297.351585, 285.41239, 273.197914, 263.040569, 253.137348, 243.408333, 233.77361, 224.009152, 214.236795, 204.434269, 194.579301, 184.439194, 174.28627, 164.182425, 154.189557, 144.710428, 135.329724, 125.972995, 116.565793, 106.761765, 96.8671266, 86.916191, 76.9432708, 67.0016617, 57.0990999, 47.2623045, 37.5179949, 28.466852, 19.3320488, 9.91071997, 0., -11.5081985, -23.54943, -35.9644712, -48.5940986, -60.1461305, -72.0474852, -84.5921225, -98.0740026, -111.204635, -126.493411, -144.867268, -167.253148, -184.019583, -210.875282, -252.97055, -315.455689, -335.026177, -412.669073, -580.916609, -872.30102, -1244.6194, -1954.6094, -3090.28309, -4876.46181, -7258.87633, -11162.8104, -17268.6597, -26256.8197}, 
{332.848826, 319.339406, 309.84285, 303.496527, 299.437804, 296.804049, 294.732628, 292.360909, 288.82626, 279.6239, 268.990204, 257.519397, 245.805708, 236.102767, 226.681634, 217.472775, 208.406655, 199.343792, 190.312576, 181.27145, 172.178856, 162.676833, 153.16679, 143.733731, 134.462662, 125.865105, 117.42894, 109.068565, 100.698377, 91.9270899, 83.0970604, 74.2449607, 65.4074631, 56.6069349, 47.9000748, 39.3292765, 30.936934, 23.542257, 16.1000972, 8.34212235, 0., -10.4410908, -21.5043985, -32.9596599, -44.576612, -54.88014, -65.3827734, -76.3521899, -88.0560672, -99.2001149, -112.238766, -128.064485, -147.569738, -161.556655, -185.044169, -222.960877, -280.235378, -293.441953, -363.205243, -521.795574, -801.483269, -1160.10306, -1853.23206, -2962.50825, -4729.38986, -7132.2067, -11041.982, -17106.5658, -25973.8078}, 
{304.952152, 292.952218, 284.33093, 278.349426, 274.268842, 271.350315, 268.854983, 266.043982, 262.17845, 253.319746, 243.208695, 232.386348, 221.393753, 212.329699, 203.554397, 194.985801, 186.541861, 178.001582, 169.477443, 160.942978, 152.371717, 143.465933, 134.578921, 125.792715, 117.189349, 109.241916, 101.484971, 93.8441289, 86.2450012, 78.3468177, 70.4481283, 62.5810991, 54.7778966, 47.0124833, 39.3985105, 31.9914262, 24.8466778, 18.9719631, 13.08958, 6.87407634, 0., -9.38029515, -19.4611898, -29.9592584, -40.5910756, -49.7813155, -59.0552134, -68.6461041, -78.7873225, -88.1276716, -99.1188308, -112.627947, -129.522169, -140.631068, -160.874397, -195.134334, -248.293057, -257.393244, -322.292372, -475.008418, -747.559361, -1096.57811, -1780.23784, -2878.9862, -4641.40046, -7069.0396, -10983.4062, -17001.0171, -25738.3888}, 
{279.78566, 269.98824, 262.640456, 257.177906, 253.036189, 249.650905, 246.457651, 242.892026, 238.389629, 229.788609, 220.160994, 209.981364, 199.724297, 191.311845, 183.192125, 175.260728, 167.413243, 159.281274, 151.129992, 142.960581, 134.774227, 126.404438, 118.087145, 109.890602, 101.883063, 94.4378202, 87.1960756, 80.1040686, 73.1080386, 65.9558819, 58.8715185, 51.8805254, 45.0084794, 38.1787473, 31.56, 25.2186982, 19.2213027, 14.7377201, 10.2895871, 5.50198627, 0., -8.3384809, -17.4423055, -26.9875147, -36.6501492, -44.7971021, -52.9372211, -61.2702062, -69.9957574, -77.6869199, -86.8207102, -98.2474902, -112.817622, -121.191888, -138.48606, -169.626333, -219.538899, -226.031281, -287.995812, -437.206154, -705.435969, -1047.95037, -1726.04867, -2820.76157, -4583.33267, -7035.48163, -10946.811, -16899.9525, -25477.5381}, 
{257.302813, 250.074573, 244.204199, 239.330028, 235.090398, 231.123647, 227.068111, 222.562129, 217.244038, 208.851392, 199.683625, 190.139389, 180.617334, 172.839114, 165.351176, 158.02297, 150.723945, 142.906536, 135.024016, 127.112641, 119.208668, 111.321731, 103.525358, 95.8664567, 88.3919317, 81.3337627, 74.4797527, 67.8027781, 61.2757154, 54.7583607, 48.3819031, 42.1644511, 36.1241133, 30.1364421, 24.4191248, 19.0472927, 14.0960769, 10.8651988, 7.71536327, 4.23186531, 0., -7.30986874, -15.4315422, -24.013753, -32.7052334, -39.8530662, -46.9282936, -54.100308, -61.5385019, -67.7347666, -75.206996, -84.795583, -97.3409206, -103.184816, -117.865682, -146.423345, -193.897633, -199.101319, -259.792105, -407.500641, -673.757574, -1012.59923, -1688.03923, -2782.74657, -4546.80967, -7020.11695, -10915.9374, -16779.1949, -25154.8133}, 
{237.457073, 232.838319, 228.454933, 224.153853, 219.782019, 215.186369, 210.213843, 204.71138, 198.525917, 190.328997, 181.613115, 172.69537, 163.892858, 156.701414, 149.787905, 142.997937, 136.177112, 128.601035, 120.91331, 113.187537, 105.497321, 98.0471923, 90.727455, 83.5593415, 76.5640841, 69.8106881, 63.253503, 56.8946518, 50.7362571, 44.7583326, 38.993954, 33.4540878, 28.1497004, 22.9162824, 18.0104665, 13.5134097, 9.5062687, 7.38006997, 5.38215293, 3.0697266, 0., -6.28867938, -13.4126964, -21.0072975, -28.7077289, -34.8747742, -40.9279278, -47.0122213, -53.2726862, -58.1281187, -64.1402798, -72.1446952, -82.9768906, -86.555554, -98.9997838, -125.511841, -171.293985, -176.348612, -237.157795, -385.00374, -651.168655, -988.904083, -1663.58422, -2759.85339, -4523.45463, -7011.52972, -10874.5267, -16614.567, -24733.7721}, 
{221.112345, 217.273009, 213.245771, 208.971399, 204.390662, 199.444326, 194.07316, 188.21793, 181.819406, 174.07591, 165.967632, 157.732316, 149.607708, 142.759168, 136.12578, 129.574241, 122.971252, 115.636065, 108.201803, 100.754142, 93.3787588, 86.3649399, 79.5133073, 72.8280932, 66.31353, 59.9400399, 53.7591891, 47.7887339, 42.0464305, 36.5675385, 31.3453094, 26.3904981, 21.7138595, 17.1544906, 12.9634673, 9.22020788, 6.00413027, 4.68026363, 3.52817058, 2.11302482, 0., -5.20156086, -11.1883212, -17.612775, -24.1274164, -29.1515342, -34.0641099, -39.0109194, -44.1377391, -47.8775291, -52.7740078, -59.6580777, -69.3606412, -71.3464121, -82.3589571, -107.775654, -152.973881, -158.692157, -220.802263, -370.537118, -639.129644, -981.466477, -1657.81939, -2751.51569, -4503.96735, -6999.36377, -10792.7001, -16321.2699, -24022.367}, 
{207.178203, 203.732178, 199.820678, 195.45373, 190.64136, 185.393593, 179.720456, 173.631975, 167.138176, 159.874527, 152.375435, 144.800749, 137.310317, 130.721435, 124.273526, 117.863461, 111.388109, 104.253972, 97.0444365, 89.8525206, 82.7712422, 76.1562566, 69.7328892, 63.4891028, 57.4128605, 51.3843422, 45.5424068, 39.9181301, 34.5425882, 29.4933277, 24.7363656, 20.2841896, 16.1492874, 12.1875237, 8.63065835, 5.55382834, 3.03217056, 2.36462843, 1.91300975, 1.26292884, 0., -4.1202108, -8.95002157, -14.1717987, -19.4679084, -23.3384285, -27.1209295, -30.9706931, -35.0430011, -37.8103107, -41.7838582, -47.7920555, -56.6633142, -57.3923226, -67.3747053, -92.1723635, -137.347198, -144.531429, -208.788512, -361.25222, -633.056326, -981.531677, -1659.22083, -2747.46148, -4481.02972, -6970.9883, -10669.8149, -15944.6335, -23162.5681}, 
{195.375656, 192.004692, 188.015746, 183.466195, 178.413415, 172.914785, 167.027681, 160.809479, 154.317557, 147.537176, 140.626675, 133.672277, 126.760207, 120.360119, 114.021432, 107.676999, 101.25967, 94.2939214, 87.28433, 80.3270978, 73.5184264, 67.2620822, 61.2236764, 55.3763849, 49.6933835, 43.9759018, 38.4378405, 33.1211539, 28.0677965, 23.3846307, 19.0227396, 14.9981141, 11.3267455, 7.89149902, 4.89474193, 2.40571575, 0.493661944, 0.363790214, 0.494986523, 0.502105057, 0., -3.05610252, -6.72786903, -10.7365941, -14.8035724, -17.5234724, -20.1948652, -22.9896958, -26.0799092, -27.9994607, -31.213481, -36.549111, -44.8334916, -44.6143205, -53.9099591, -78.4593255, -124.001338, -133.312004, -200.278317, -355.824359, -630.8742119999999, -985.890452, -1663.18168, -2741.44858, -4446.50068, -6914.71693, -10493.0898, -15475.9586, -22157.6622}, 
{185.425713, 181.879415, 177.667063, 172.87414, 167.58613, 161.888515, 155.86678, 149.606408, 143.192881, 136.876184, 130.511499, 124.11851, 117.716901, 111.447124, 105.159785, 98.8262595, 92.4179225, 85.5950808, 78.7646057, 72.0222999, 65.4639662, 59.5233565, 53.8231446, 48.3239537, 42.9864068, 37.5470253, 32.2801747, 27.2361189, 22.4651217, 18.090378, 14.0600481, 10.3952234, 7.11699527, 4.1425339, 1.63842043, -0.332685234, -1.7081232, -1.39162518, -0.767583015, -0.186780138, 0., -2.0207093, -4.55193494, -7.35938705, -10.2087758, -11.7946811, -13.3823954, -15.166081, -17.3399, -18.5179766, -21.1065264, -25.9317269, -33.8197554, -32.9334411, -41.8276492, -66.3938965, -112.5237, -124.479458, -194.433453, -352.92885, -630.508812, -991.333571, -1665.09509, -2727.23475, -4392.23916, -6818.86329, -10249.7435, -14906.5457, -21010.9359}, 
{177.049383, 173.145211, 168.610719, 163.542914, 158.038804, 152.195397, 146.109701, 139.878725, 133.599477, 127.703878, 121.820057, 115.911057, 109.93992, 103.754355, 97.4788721, 91.1226459, 84.6948522, 77.9966166, 71.3283849, 64.7825529, 58.4515163, 52.7810191, 47.3687694, 42.1658237, 37.1232385, 31.9300195, 26.9040941, 22.1013388, 17.5776303, 13.4595, 9.70390802, 6.33846907, 3.39079796, 0.816745601, -1.25560381, -2.76992995, -3.66991249, -2.97099193, -1.91638277, -0.821060353, 0., -1.02550445, -2.45229068, -4.09240305, -5.75788596, -6.24006989, -6.77999864, -7.59800213, -8.91441027, -9.4388557, -11.5066445, -15.9423856, -23.5706878, -22.2707196, -30.9907063, -55.7334329, -102.501684, -117.479367, -190.415696, -351.241008, -629.885639, -994.651801, -1660.35421, -2698.57779, -4310.10409, -6671.741, -9926.99437, -14227.6955, -19725.6759}, 
{187.82453, 175.380638, 164.522254, 155.048751, 146.759502, 139.453881, 132.931261, 126.991014, 121.432515, 115.682444, 110.061944, 104.519465, 99.0034562, 93.4094581, 87.7599952, 82.0246817, 76.1731321, 69.8497126, 63.4793845, 57.1618613, 50.9968561, 45.3113449, 39.8868738, 34.732251, 29.8562852, 25.3107344, 21.0442775, 17.0485431, 13.3151599, 9.78749463, 6.52474233, 3.53783623, 0.837709548, -1.92123783, -4.22892593, -5.93180815, -6.87633789, -5.53361613, -3.67558966, -1.69885283, 0., -0.623134421, -1.65833802, -2.8432016, -3.91531598, -3.23888229, -2.47423689, -1.90832646, -1.82809767, -0.571915493, -1.15474102, -4.64295363, -12.1029327, -13.5154792, -25.4667823, -53.4574527, -102.988101, -117.986443, -190.155143, -349.623968, -626.522688, -999.557189, -1653.12889, -2646.10523, -4169.23759, -6411.62763, -9437.09199, -13340.7686, -18217.7952}, 
{202.594033, 179.751196, 161.257808, 146.529047, 134.980092, 126.026122, 119.082316, 113.563853, 108.885911, 103.059218, 97.4651852, 92.0807714, 86.8829361, 82.0677256, 77.3053773, 72.4852159, 67.4965658, 61.6813602, 55.6952715, 49.6465805, 43.6435684, 37.8332563, 32.2696889, 27.0456511, 22.2539279, 18.4601648, 15.0951418, 12.0624993, 9.26587801, 6.39380025, 3.65107215, 1.02738158, -1.48758357, -4.64361367, -7.41575128, -9.51851721, -10.6664323, -8.53514926, -5.69360421, -2.67186512, 0., -0.508669667, -1.50711222, -2.60515855, -3.41263955, -1.63513838, 0.451567238, 2.47594729, 4.06647178, 7.34644229, 8.45156461, 6.01237612, -1.34058582, -6.53591649, -22.7602928, -54.7595243, -107.279421, -121.97668, -191.921867, -347.096438, -617.481845, -993.994698, -1623.81099, -2550.28408, -3950.79094, -6012.26715, -8749.39028, -12219.9024, -16481.5456}, 
{202.282333, 175.871048, 154.829792, 138.423698, 125.9179, 116.577528, 109.667717, 104.453599, 100.200306, 94.3003516, 88.640535, 83.2350361, 78.0980346, 73.6586754, 69.3501876, 65.0207652, 60.5186022, 55.0509575, 49.3633342, 43.5603002, 37.7464236, 31.9534243, 26.3878575, 21.1834305, 16.4738504, 13.064521, 10.1487743, 7.59163891, 5.25814349, 2.701185, 0.222776282, -2.18720142, -4.53886688, -7.71996386, -10.5119361, -12.5738524, -13.5647814, -10.8155658, -7.24479091, -3.44281592, 0., -0.08641868, -0.680828277, -1.34170145, -1.62751084, 1.08043447, 4.17563281, 7.22874608, 9.81043621, 14.0900627, 16.0001108, 14.0717634, 6.83620339, -0.0235436439, -17.9888746, -51.3893436, -104.554505, -121.858854, -191.569027, -341.996602, -601.453156, -968.32838, -1560.69951, -2412.60149, -3678.80089, -5530.33367, -7952.42665, -10976.7589, -14635.0096}, 
{191.434042, 166.516391, 146.685199, 131.231755, 119.447344, 110.623254, 104.050772, 99.021184, 94.8257772, 88.8485881, 83.051054, 77.4873619, 72.211699, 67.8168808, 63.6030143, 59.4088351, 55.0730788, 49.8012396, 44.318591, 38.7171651, 33.088994, 27.4028755, 21.9233699, 16.7918031, 12.1495011, 8.82957003, 6.00484342, 3.539935, 1.29945857, -1.20513783, -3.61480841, -5.9236731, -8.12585181, -11.0434715, -13.5114422, -15.1926812, -15.7501056, -12.5223676, -8.42635531, -4.05469173, 0., 0.584774086, 0.691459186, 0.751561269, 1.1965863, 4.6894944, 8.53775573, 12.2802946, 15.4560354, 19.9833366, 22.0699146, 20.30292, 13.2695034, 6.39295062, -11.3101771, -43.9871831, -95.785371, -117.497179, -187.966722, -332.683249, -577.13601, -923.480797, -1467.20723, -2237.65049, -3360.76265, -4981.95802, -7072.18514, -9646.00898, -12717.9945}, 
{174.593769, 154.463423, 138.271025, 125.452267, 115.442846, 107.678454, 101.594787, 96.6275385, 92.2124029, 86.1466711, 80.1598024, 74.3428524, 68.7868768, 64.1769151, 59.7724458, 55.426931, 50.993833, 45.7749414, 40.3960599, 34.9313197, 29.4548519, 23.9126366, 18.5582165, 13.5169832, 8.91432833, 5.46107933, 2.46301766, -0.188639368, -2.60267448, -5.23995495, -7.71634509, -9.99979375, -12.0582498, -14.5073195, -16.4082314, -17.4698714, -17.4011254, -13.8030567, -9.33550295, -4.55047905, 0., 1.44606418, 2.48069555, 3.47902116, 4.81616805, 8.97369958, 13.3777319, 17.5608176, 21.055509, 25.3506551, 27.2399934, 25.4735579, 18.8013827, 13.0848776, -2.88184972, -33.1933156, -81.9440363, -108.755869, -179.985053, -317.515168, -543.229795, -860.374511, -1346.74691, -2030.02411, -3004.17146, -4383.27103, -6134.6498, -8262.3235, -10770.3079}, 
{156.306122, 142.488343, 131.034262, 121.584286, 113.778824, 107.258283, 101.66307, 96.6335931, 91.81026, 85.6373442, 79.4298407, 73.3066107, 67.3865153, 62.3733517, 57.5670705, 52.8525584, 48.1147021, 42.814798, 37.4307592, 32.0169088, 26.6275694, 21.2137341, 15.9743875, 11.0051847, 6.40178052, 2.66481611, -0.677034553, -3.69011116, -6.44075338, -9.31805307, -11.9364969, -14.2333237, -16.1457723, -18.0046911, -19.1962656, -19.5002906, -18.6965613, -14.805135, -10.0694394, -4.97316436, 0., 2.43860714, 4.5578262, 6.6450698, 8.98775053, 13.7147082, 18.5353573, 23.0005397, 26.6610969, 30.516409, 32.0893643, 30.3513891, 24.2739098, 20.4235486, 7.13845799, -19.6480136, -64.0025181, -95.4991373, -166.494121, -294.851148, -498.433899, -779.932084, -1202.73129, -1794.31537, -2616.52254, -3750.40353, -5165.80466, -6860.37338, -8831.75711}, 
{141.115714, 133.367348, 126.421904, 120.126861, 114.329698, 108.877894, 103.618928, 98.4002787, 93.0694261, 86.7633508, 80.3242294, 73.8837403, 67.5735619, 62.0407637, 56.6954766, 51.4632226, 46.2695237, 40.7635444, 35.2577072, 29.788077, 24.3907188, 19.0371946, 13.8538735, 8.90262193, 4.24530605, 0.146547664, -3.61564479, -7.06050731, -10.2072759, -13.3542189, -16.129927, -18.4420232, -20.1981307, -21.4287693, -21.8695065, -21.3788065, -19.8151337, -15.6761045, -10.7253702, -5.36573418, 0., 3.50355851, 6.79379651, 10.0540987, 13.4678499, 18.6941784, 23.8504278, 28.5296855, 32.3250388, 35.8049895, 37.1970445, 35.7041252, 30.5291532, 28.7802748, 18.5930966, -3.99155006, -42.9328337, -77.5911981, -146.364027, -263.049978, -441.44771, -683.076077, -1038.57315, -1535.11731, -2205.31111, -3099.48633, -4191.63378, -5474.82953, -6942.14962}, 
{103.940274, 109.246159, 112.100836, 112.786499, 111.585343, 108.779563, 104.651353, 99.482908, 93.556423, 87.4462621, 81.0255828, 74.4597118, 67.9139761, 62.0289608, 56.3046312, 50.7162112, 45.2389241, 39.780798, 34.4111305, 29.1320235, 23.945579, 18.7904806, 13.7576159, 8.87445415, 4.16846479, -0.307552597, -4.56159118, -8.57631379, -12.3343833, -16.2187059, -19.6516038, -22.4556421, -24.4533865, -25.0813324, -24.702543, -23.294012, -20.8327326, -16.4901945, -11.3710965, -5.77463332, 0., 4.75478168, 9.44887446, 14.1426142, 18.8963367, 25.0263267, 30.8345917, 35.8790878, 39.7177713, 41.4854709, 41.3325214, 38.9861302, 34.1735046, 36.3321689, 31.5948865, 15.8047379, -15.1951963, -53.3876082, -119.973336, -223.978989, -374.43118, -567.576215, -850.781612, -1259.91447, -1784.82538, -2422.06506, -3187.40028, -4091.7938, -5146.20839}, 
{119.696459, 121.904076, 122.249633, 120.94698, 118.209966, 114.25244, 109.288254, 103.531255, 97.1952948, 90.7420274, 84.0383753, 77.1990662, 70.3388277, 63.8365474, 57.4371289, 51.1496357, 44.9831315, 38.9566103, 33.0652331, 27.3140912, 21.708276, 16.2226188, 10.9045753, 5.77134083, 0.84011102, -3.93325859, -8.44569664, -12.6554717, -16.5208523, -20.3073567, -23.5431039, -26.0634625, -27.7038011, -27.8850365, -27.02277, -25.1181511, -22.1723292, -17.5689495, -12.1736673, -6.23463379, 0., 5.72320321, 11.4692564, 17.2135602, 22.9315156, 29.3276159, 35.3565322, 40.7020283, 45.0478681, 48.0136362, 49.3729471, 48.8352361, 46.1099386, 48.6623724, 45.3437374, 32.7611163, 7.52159136, -26.3360594, -82.6091273, -167.663208, -287.863898, -441.843924, -659.167485, -966.431705, -1346.42536, -1792.66723, -2299.56444, -2861.30205, -3472.06512}, 
{169.150239, 157.676836, 147.513519, 138.460033, 130.316123, 122.881534, 115.956011, 109.3393, 102.831146, 95.6550193, 88.4174499, 81.1486922, 73.8790011, 66.6745961, 59.5153815, 52.4172262, 45.3959992, 38.4031129, 31.5446755, 24.8623387, 18.3977542, 12.2411791, 6.36621736, 0.795078651, -4.45002757, -9.41572308, -13.9834346, -18.10342, -21.7259373, -24.9162069, -27.4635394, -29.2722079, -30.2464854, -30.1545278, -29.0911722, -27.0151386, -23.8851471, -18.9590206, -13.1767348, -6.77736834, 0., 6.50069279, 13.0674696, 19.6274911, 26.1079182, 32.5190318, 38.6716244, 44.4596088, 49.7768977, 55.8502074, 60.7075256, 63.7096435, 64.2173524, 64.9421617, 60.5538573, 49.0729433, 28.5199236, 5.13411772, -34.5703119, -95.8603875, -184.003131, -305.004472, -461.914711, -660.672514, -900.983227, -1205.90322, -1511.31808, -1774.90607, -1954.34547}, 
{157.588269, 149.87012, 142.397846, 135.138514, 128.05919, 121.126941, 114.308835, 107.571936, 100.883312, 94.184065, 87.4776118, 80.7414055, 73.952899, 66.9294522, 59.8726476, 52.8239751, 45.8249243, 38.9794848, 32.2416468, 25.6279002, 19.1547349, 12.7006072, 6.47525399, 0.550378712, -5.00231531, -9.89727766, -14.3621909, -18.4108905, -22.0572119, -25.8017927, -28.9769454, -31.4017847, -32.8954252, -32.8713233, -31.7165151, -29.4123786, -25.9402917, -20.6885535, -14.4688524, -7.49979804, 0., 7.57072875, 15.3314634, 23.1600759, 30.9344379, 38.6372245, 45.9995829, 52.8574638, 59.0468175, 64.5920658, 69.0652995, 72.2270807, 73.8379713, 76.5085846, 76.0094109, 70.9609919, 59.9838691, 45.5410236, 20.8735817, -16.934891, -70.8008287, -145.218014, -238.370836, -359.648158, -484.977861, -600.901629, -684.966921, -716.969752, -676.706137}, 
{170.100056, 162.077358, 154.716959, 147.843675, 141.282325, 134.857726, 128.394696, 121.718051, 114.652611, 106.038997, 97.0798998, 87.993815, 78.9992377, 71.0318317, 63.3060567, 55.7535409, 48.3059124, 40.5721135, 32.9355325, 25.4568718, 18.1968339, 11.2698263, 4.66136436, -1.58933151, -7.44304092, -12.8267762, -17.748591, -22.1827721, -26.1036058, -29.8062716, -32.815806, -34.9781382, -36.1391978, -35.7891867, -34.2720525, -31.5760155, -27.6892961, -22.1624404, -15.5964124, -8.15450233, 0., 8.75252273, 17.8715705, 27.174366, 36.478132, 45.3117409, 53.896106, 62.1637903, 70.0473565, 77.8136669, 84.9272654, 91.1869954, 96.3917, 101.4355, 104.58385, 105.197482, 102.63713, 95.6793504, 84.50272, 68.7016409, 47.8705152, 17.3576345, -10.504268, -36.2030578, -39.5162004, -13.1099957, 81.1941944, 273.863774, 595.366146}, 
{177.272638, 176.171898, 172.67048, 167.127495, 159.902055, 151.353269, 141.840249, 131.722107, 121.357952, 111.932907, 102.649669, 93.5369434, 84.6234388, 76.1690857, 67.8788777, 59.6890322, 51.5357667, 42.9013283, 34.3574928, 26.0220658, 18.0128529, 10.8353733, 4.06463329, -2.336647, -8.40574756, -14.5091726, -20.223288, -25.4536843, -30.1059515, -34.3111544, -37.6592191, -39.965546, -41.0455356, -40.0655795, -37.7496906, -34.1728732, -29.4101313, -23.5468726, -16.6435361, -8.77096437, 0., 10.1187145, 20.7860564, 31.7231031, 42.6509319, 52.6140477, 62.2807294, 71.6426832, 80.6916156, 89.0636486, 97.2483064, 105.379529, 113.591258, 122.654121, 131.810694, 140.94024, 149.922024, 154.37936, 160.149841, 168.815111, 181.956813, 197.339586, 227.99609, 274.768336, 359.778501, 504.715443, 711.721124, 987.824231, 1340.05345}
};
  
  /* Stuff for interpolating the NQC data */
  gsl_spline    *spline = NULL;
  gsl_interp_accel *acc = NULL;
  /* Interpolate the spin NQC data in 2-D parameter space -- spin and mass ratio */
  /* First, interpolate in spin dimension for all mass ratios */
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
//printf("%.15f\n",a3slist[i]);
  }
//printf("%.15f\n",a);
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  /* Second, interpolate in mass ratio dimension */
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
 
  /* Interpolate nonspin NQC data in the mass ratio dimension */
  spline = gsl_spline_alloc( gsl_interp_cspline, nsqdim );
  acc    = gsl_interp_accel_alloc();
  gsl_spline_init( spline, nsetalist, a1list, nsqdim );
  gsl_interp_accel_reset( acc );
  coeffs->a1 = gsl_spline_eval( spline, eta, acc );
  gsl_spline_init( spline, nsetalist, a2list, nsqdim );
  gsl_interp_accel_reset( acc );
  coeffs->a2 = gsl_spline_eval( spline, eta, acc );
  gsl_spline_init( spline, nsetalist, a3list, nsqdim );
  gsl_interp_accel_reset( acc );
  coeffs->a3 = gsl_spline_eval( spline, eta, acc );
  gsl_spline_init( spline, nsetalist, b1list, nsqdim );
  gsl_interp_accel_reset( acc );
  coeffs->b1 = gsl_spline_eval( spline, eta, acc );
  gsl_spline_init( spline, nsetalist, b2list, nsqdim );
  gsl_interp_accel_reset( acc );
  coeffs->b2 = gsl_spline_eval( spline, eta, acc );
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  /* Andrea and I have different sign conventions, so I need to put a minus sign in front */
  coeffs->b1 = - coeffs->b1;
  coeffs->b2 = - coeffs->b2;
  /*coeffs->a1 = -6.292993135543338;
  coeffs->a2 = 40.26795686019975;
  coeffs->a3 = -39.29214347590138;
  coeffs->a3S= 226.3844538692294;
  coeffs->a4 = -747.9690610155075;
  coeffs->a5 = 623.9090070111784;
  coeffs->b1 = 0.7467610798775989;
  coeffs->b2 = -2.159938635119833;
  */
  /* Obsolete polynomial fitting of nonspin NQC coefficients a1, a2, a3, b1 and b2 */
  /*
  coeffs->a1 = -12.67955358602124 + 75.41927959573084 * eta - 106.15933052937714 * eta2;
  coeffs->a2 = 101.45522216901628 - 757.3158549733314 * eta + 1473.314771676588 * eta2;
  coeffs->a3 = -107.6647834845902 + 857.6219519536213 * eta - 1776.2776804623143 * eta2;
  // Andrea and I have different sign conventions, so I need to put a minus sign in front 
  coeffs->b1 = - (-1.464129495621165 + 12.81732978488213 * eta - 60.09957767247623 * eta2);
  coeffs->b2 = - ( 7.477426352542122 - 85.26122117590637 * eta + 353.3251639728075 * eta2);
  */

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

  /* nrDeltaT defined in XLALSimIMREOBGetNRSpinPeakDeltaT is a minus sign different from Eq. (33) of Taracchini et al.
   * Therefore, the plus sign in Eq. (21) of Taracchini et al and Eq. (18) of DCC document T1100433v2 is 
   * changed to a minus sign here.
   */
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
