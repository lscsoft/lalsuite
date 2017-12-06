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

#include <complex.h>
#include <math.h>

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

#include "LALSimIMREOBNQCTables.c"

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
/*printf("a1 = %.16e, a2 = %.16e, a3 = %.16e, a3S = %.16e, a4 = %.16e, a5 = %.16e\n",coeffs->a1,coeffs->a2,coeffs->a3,coeffs->a3S, coeffs->a4,coeffs->a5);
printf("b1 = %.16e, b2 = %.16e, b3 = %.16e, b4 = %.16e\n",coeffs->b1,coeffs->b2,coeffs->b3,coeffs->b4);*/
  /* In EOBNRv2, coeffs->a3S, coeffs->a4 and coeffs->a5 are set to zero */
  /* through XLALSimIMREOBGetCalibratedNQCCoeffs() */
  /* and XLALSimIMREOBCalculateNQCCoefficients() */
  mag = 1. + (p*p / rOmegaSq) * ( coeffs->a1
     + coeffs->a2 / r + ( coeffs->a3 + coeffs->a3S) / (r*sqrtR)
     + coeffs->a4 / (r*r) + coeffs->a5 / (r*r*sqrtR));
//printf("NQC INFO mag = %.16e, r = %.16e, p = %.16e\n",mag,r,p);
  phase = coeffs->b1 * p / rOmega + p*p*p/rOmega * ( coeffs->b2
     + coeffs->b3 / sqrtR + coeffs->b4 / r );

  *nqc = mag * cos(phase);
  *nqc += I * mag * sin(phase);
/*printf("r = %.16e, pr = %.16e, omega = %.16e\n",r,p,omega);
printf("NQC mag = %.16e, arg = %.16e\n",mag,phase);*/
  return XLAL_SUCCESS;

}



/**
 * This function calculates the non-quasicircular correction to apply to
 * the waveform. The form of this correction can be found in Pan et al,
 * PRD84, 124052(2011), Eq.(22), and also in the DCC document T1100433. Note
 * that when calling this function, the NQC coefficients should already
 * have been pre-computed.
 * This version is for generic precesing case where the dynamics variable
 * values are given in Catesean coordinates.
 */
UNUSED static int  XLALSimIMRSpinEOBNonQCCorrection(
                      COMPLEX16      * restrict nqc,    /**<< OUTPUT, The NQC correction */
                      REAL8Vector    * restrict values, /**<< Dynamics r, phi, pr, pphi */
                      const REAL8               omega,  /**<< Angular frequency */
                      EOBNonQCCoeffs * restrict coeffs  /**<< NQC coefficients */
                     )

{

  REAL8 rOmega, rOmegaSq;
  REAL8 r, p, sqrtR;

  REAL8 mag, phase;


  r = sqrt(values->data[0]*values->data[0]+values->data[1]*values->data[1]+values->data[2]*values->data[2]);
  p = (values->data[0]*values->data[3]+values->data[1]*values->data[4]+values->data[2]*values->data[5]) / r;

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

  *nqc = mag * cos(phase);
  *nqc += I * mag * sin(phase);

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
 *          Spin (SEOBNRv1 and SEOBNRv2)
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
 * The time difference between the orbital peak and the peak amplitude
 * of the mode in question (currently only 2,2 implemented ).
 */
UNUSED static inline REAL8 XLALSimIMREOBGetNRSpinPeakDeltaTv2(
                 INT4 UNUSED l,    /**<< Mode l */
                 INT4 UNUSED m,    /**<< Mode m */
                 REAL8 UNUSED m1,  /**<< mass 1 */
                 REAL8 UNUSED m2,  /**<< mass 22 */
                 REAL8 chi1,       /**<< Dimensionless spin1 */
                 REAL8 chi2        /**<< Dimensionless spin2 */
                 )
{
  REAL8 chi, chichi;
  REAL8 eta = m1*m2 / ((m1+m2)*(m1+m2));
  chi    = (chi1+chi2)/2. + (chi1-chi2)/2.*((m1-m2)/(m1+m2))/(1-2.*eta);
    
  chichi = chi*chi;
  if ( chi > 0.8 )
  {
    return (0.75*eta*chi + sqrt(1.-4.*eta)) * (57.1755-48.0564*chi);
  }
  else if ( chi > 0.0 )
  {
    return (0.75*eta*chi + sqrt(1.-4.*eta)) * (2.5+10.*chichi+24.*chichi*chichi);
  }
  else
  {
    return 2.5 + (1.+2.5*chi) * (-2.5+2.5*sqrt(1.-4.*eta));
  }
}

/**
 * Peak frequency predicted by fitting NR results (currently only 2,2 available).
 * Take from unpublished SEOBNRv2 results.
 */
UNUSED static inline REAL8 GetNRSpinPeakOmegav2( INT4 UNUSED l, INT4 UNUSED m, REAL8 UNUSED eta, REAL8 a )
{
  REAL8 chi = a/(1.0 - 2.0*eta);
  REAL8 eta2= eta*eta;
  if ( eta > 50./51./51. )
  {
    return 0.43747541927878864 + (-0.10933208665273314 - 0.007325831113333813*chi)*log(
            4.500844771420863 - 9.681916048928946*eta +
            chi*(-4.254886879579986 + 11.513558950322647*eta));
  }
  else
  {
    return 1.5609526077704716 - 122.25721149839733 * eta +
           3586.2192688666914 * eta2 - 13869.506144441548 * eta*eta2 +
           (eta - 0.25) * (1651.5823693445805 * (-0.01922337588094282 + eta) *
           (-0.01922337536857659 + eta) + 66.87492814925524 * chi *
           (0.0003695381704106058 - 0.03844675124951941 * eta + eta2)) *
           log(5600.67382718678 - 5555.824895398546 * chi) +
           (-1412.8186461833657 + 67.66455403259023 * chi) * (eta - 0.001) *
           (0.0003695381704106056 - 0.038446751249519406 * eta + eta2) *
           log(0.5680439481719505 - 0.36813967358200156 * chi) +
           0.012328326527732041 * log(4.500844771420863 - 9.681916048928946 * eta + 
           chi * (-4.254886879579986 + 11.513558950322647 * eta)) +
           0.0008260634258180991 * chi * log(4.500844771420863 -9.681916048928946 * eta +
           chi * (-4.254886879579986 + 11.513558950322647 * eta)) -
           12.6575493872956 * eta * log(4.500844771420863 -
           9.681916048928946 * eta + chi * (-4.254886879579986 + 11.513558950322647 * eta)) -
           0.8481231078533651 * chi * eta * log(4.500844771420863 - 
           9.681916048928946 * eta + chi * (-4.254886879579986 + 11.513558950322647 * eta)) +
           329.2228595635586 * eta2 * log(4.500844771420863 - 9.681916048928946 * eta +
           chi * (-4.254886879579986 + 11.513558950322647 * eta)) +
           22.05968203526603 * chi * eta2 * log(4.500844771420863 -9.681916048928946 * eta +
           chi * (-4.254886879579986 + 11.513558950322647 * eta));
  }
}

/**
 * Peak frequency slope predicted by fitting NR results (currently only 2,2 available).
 * Take from unpublished SEOBNRv2 results.
 */
UNUSED static inline REAL8 GetNRSpinPeakOmegaDotv2( INT4 UNUSED l, INT4 UNUSED m, REAL8 UNUSED eta, REAL8 UNUSED a )
{
  REAL8 chi = a/(1.0-2.0*eta);
  REAL8 eta2= eta*eta;
  /* Fit for HOMs missing */
  if (chi < 0.8 )
  {
    return -0.07086074186161867 * chi * (-0.26367236731979804 + eta) *
           (-0.0010019969893089581 + eta) + 0.2893863668183948 *
           (-0.16845695144529893 + eta) * (0.23032241797163952 + eta) +
           (0.004086861548547749 - 0.06538978477676398 * eta2 +
           chi * (0.0006334026884930817 - 0.010134443015889307 * eta2)) *
           log(68.47466578101876 - 58.30148755701496 * chi);
  }
  else
  {
    return -0.10069512275335238 * (-0.46107388514323044 + eta) *
           (0.2832795481380979 + eta) + 0.2614619716504706 * chi *
           (-0.24838163750494138 + eta) * (0.320112993649413 + eta) +
           chi * chi * (0.010000160002560042 - 0.16000256004096067 * eta2);
  }
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
                 REAL8                     m1,          /**<< Component mass 1 */
                 REAL8                     m2,          /**<< Component mass 2 */
                 REAL8                     a,           /**<< Normalized spin of deformed-Kerr */
                 REAL8                     chiA,        /**<< Assymmetric dimensionless spin combination */
                 REAL8                     chiS,        /**<< Symmetric dimensionless spin combination */
                 EOBNonQCCoeffs * restrict coeffs,      /**<< OUTPUT, NQC coefficients */
                 UINT4                     SpinAlignedEOBversion  /**<< 1 for SEOBNRv1, 2 for SEOBNRv2 */
)
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

  REAL8 eta = (m1*m2)/((m1+m2)*(m1+m2));
  REAL8 amp, aDot, aDDot;
  REAL8 omega, omegaDot;

  REAL8 qNSLMPeak, qNSLMDot, qNSLMDDot;
  REAL8 pNSLMDot, pNSLMDDot;

  REAL8 nra, nraDDot;
  REAL8 nromega, nromegaDot;

  REAL8 nrDeltaT, nrTimePeak;
  REAL8 chi1 = chiS + chiA;
  REAL8 chi2 = chiS - chiA;

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
  switch ( SpinAlignedEOBversion)
  {
   case 1:
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
     break;
   case 2:
     // if ( XLALSimIMRGetEOBCalibratedSpinNQCv2( coeffs, l, m, eta, a ) == XLAL_FAILURE )
     if ( XLALSimIMRGetEOBCalibratedSpinNQC3D( coeffs, l, m, m1, m2, a, chiA ) == XLAL_FAILURE )
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
     break;
   default:
     XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
     XLAL_ERROR( XLAL_EINVAL );
     break;
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
  switch ( SpinAlignedEOBversion )
  {
   case 1:
     nrDeltaT   = XLALSimIMREOBGetNRSpinPeakDeltaT( l, m, eta, a );
     break;
   case 2: 
     nrDeltaT   = XLALSimIMREOBGetNRSpinPeakDeltaTv2( l, m, m1, m2, chi1, chi2 );
     break;
   default:
     XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
     XLAL_ERROR( XLAL_EINVAL );
     break;
  }


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
  switch ( SpinAlignedEOBversion )
  {
   case 1:
     nromega = GetNRSpinPeakOmega( l, m, eta, a );
     nromegaDot = GetNRSpinPeakOmegaDot( l, m, eta, a );
     break;
   case 2:
     nromega = GetNRSpinPeakOmegav2( l, m, eta, a );
     nromegaDot = GetNRSpinPeakOmegaDotv2( l, m, eta, a );
     break;
   default:
     XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
     XLAL_ERROR( XLAL_EINVAL );
     break;
  }

  /*printf("NR inputs: %.16e, %.16e, %.16e, %.16e\n",nra,nraDDot,nromega,nromegaDot);
  printf("NR inputs: %.16e, %.16e, %.16e, %.16e\n",pNSLMDot, pNSLMDDot,omega,omegaDot);*/

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
  switch ( SpinAlignedEOBversion )
  {
   case 1:
     coeffs->b3  = gsl_vector_get( bCoeff, 0 );
     coeffs->b4  = gsl_vector_get( bCoeff, 1 );
     break;
   case 2:
     //coeffs->b3  = gsl_vector_get( bCoeff, 0 );
     //coeffs->b4  = gsl_vector_get( bCoeff, 1 );
     break;
   default:
     XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
     XLAL_ERROR( XLAL_EINVAL );
     break;
  }
  coeffs->b3  *= 1.0;
  coeffs->b4  *= 1.0;
//  coeffs->b3  = -778.891568845;
//  coeffs->b4  = 1237.46952422;
//  coeffs->b3  = -876.669217307;
//  coeffs->b4  = 1386.13223658;
//  coeffs->b3 = 41583.9402122;
//  coeffs->b4 = 68359.70064;

  /*printf( "NQC coefficients:\n" );
  printf( "{%f,%f,%f,%f,%f,%f}\n",  coeffs->a1, coeffs->a2, coeffs->a3, coeffs->a3S, coeffs->a4, coeffs->a5 );

  printf( "{%f,%f,%f,%f}\n",  coeffs->b1, coeffs->b2, coeffs->b3, coeffs->b4 );*/

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
