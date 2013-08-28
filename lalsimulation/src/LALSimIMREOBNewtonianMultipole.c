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
 * Functions to construct the Newtonian multipolar waveform as given
 * by Damour et al, Phys.Rev.D79:064004,2009.
 * All equation numbers in this file refer to equations of this paper,
 * unless otherwise specified.
 *
 * In addition to the function used to do this,
 * XLALCalculateNewtonianMultipole(), this file also contains a function
 * for calculating the standard scalar spherical harmonics Ylm.
 */

#include <complex.h>

#include <gsl/gsl_sf_gamma.h>

#include "LALSimIMREOBNRv2.h"

#ifndef _LALSIMIMRNEWTONIANMULTIPOLE_C
#define _LALSIMIMRNEWTONIANMULTIPOLE_C

static REAL8
XLALAssociatedLegendreXIsZero( const int l,
                               const int m
                             );

static int
XLALScalarSphHarmThetaPiBy2(COMPLEX16 *y,
                         INT4 l,
                         INT4  m,
                         REAL8 phi);

static int
CalculateThisMultipolePrefix(
               COMPLEX16 *prefix,
               const REAL8 m1,
               const REAL8 m2,
               const INT4 l,
               const INT4 m );


/**
 * Function which computes the various coefficients in the Newtonian
 * multipole. The definition of this can be found in Pan et al,
 * arXiv:1106.1021v1 [gr-qc]. Note that, although this function gets passed
 * masses measured in Solar masses (which is fine since it is a static function),
 * the units of mass won't matter so long as they are consistent. This is because
 * all uses of the mass within the function are normalized by the total mass.
 * The prefixes of all (l,m) modes are pre-computed and stored in a structure
 */
static int XLALSimIMREOBComputeNewtonMultipolePrefixes(
                NewtonMultipolePrefixes *prefix, /**<< OUTPUT Structure containing the coeffs */
                const REAL8             m1,      /**<< Mass of first component */
                const REAL8             m2       /**<< Nass of second component */
                )
{

  INT4 l, m;

  memset( prefix, 0, sizeof( NewtonMultipolePrefixes ) );

  for ( l = 2; l <= LALEOB_MAX_MULTIPOLE; l++ )
  {
    for ( m = 1; m <= l; m++ )
    {
      CalculateThisMultipolePrefix( &(prefix->values[l][m]), m1, m2, l, m );
    }
  }
  return XLAL_SUCCESS;
}

/**
 * This function calculates the Newtonian multipole part of the
 * factorized waveform for the EOBNRv2 model. This is defined in Eq. 4.
 */
UNUSED static int
XLALSimIMREOBCalculateNewtonianMultipole(
                 COMPLEX16 *multipole, /**<< OUTPUT, Newtonian multipole */
                 REAL8 x,              /**<< Dimensionless parameter \f$\equiv v^2\f$ */
                 UNUSED REAL8 r,       /**<< Orbital separation (units of total mass M) */
                 REAL8 phi,            /**<< Orbital phase (in radians) */
                 UINT4  l,             /**<< Mode l */
                 INT4  m,              /**<< Mode m */
                 EOBParams *params     /**<< Pre-computed coefficients, parameters, etc. */
                 )
{

   INT4 xlalStatus;

   COMPLEX16 y;

   INT4 epsilon = (l + m) % 2;

   y = 0.0;

  /* Calculate the necessary Ylm */
  xlalStatus = XLALScalarSphHarmThetaPiBy2( &y, l - epsilon, - m, phi );
  if (xlalStatus != XLAL_SUCCESS )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Special treatment for (2,1) and (4,4) modes, defined in Eq. 17ab of PRD84:124052 2011 */
  if ( (l == 4 && m == 4) || ( l == 2 && m == 1 ) )
  {
    *multipole = params->prefixes->values[l][m] * pow( x, (REAL8)(l+epsilon)/2.0 - 1.0)/r;
  }
  else
  {
    *multipole = params->prefixes->values[l][m] * pow( x, (REAL8)(l+epsilon)/2.0);
  }
  *multipole *= y;

  return XLAL_SUCCESS;
}


/**
 * This function calculates the Newtonian multipole part of the
 * factorized waveform for the SEOBNRv1 model. This is defined in Eq. 4.
 */
UNUSED static int
XLALSimIMRSpinEOBCalculateNewtonianMultipole(
                 COMPLEX16 *multipole, /**<< OUTPUT, Newtonian multipole */
                 REAL8 x,              /**<< Dimensionless parameter \f$\equiv v^2\f$ */
                 UNUSED REAL8 r,       /**<< Orbital separation (units of total mass M */
                 REAL8 phi,            /**<< Orbital phase (in radians) */
                 UINT4  l,             /**<< Mode l */
                 INT4  m,              /**<< Mode m */
                 EOBParams *params     /**<< Pre-computed coefficients, parameters, etc. */
                 )
{
   INT4 xlalStatus;

   COMPLEX16 y;

   INT4 epsilon = (l + m) % 2;

   y = 0.0;

  /* Calculate the necessary Ylm */
  xlalStatus = XLALScalarSphHarmThetaPiBy2( &y, l - epsilon, - m, phi );
  if (xlalStatus != XLAL_SUCCESS )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }

  *multipole = params->prefixes->values[l][m] * pow( x, (REAL8)(l+epsilon)/2.0) ;
  *multipole *= y;

  return XLAL_SUCCESS;
}


/* In the calculation of the Newtonian multipole, we only use
 * the spherical harmonic with theta set to pi/2. Since this
 * is always the case, we can use this information to use a 
 * faster version of the spherical harmonic code
 */

static int
XLALScalarSphHarmThetaPiBy2(
                 COMPLEX16 *y, /**<< OUTPUT, Ylm(0,phi) */
                 INT4 l,       /**<< Mode l */
                 INT4  m,      /**<< Mode m */
                 REAL8 phi     /**<< Orbital phase (in radians) */
                 )
{

  REAL8 legendre;
  INT4 absM = abs( m );

  if ( l < 0 || absM > (INT4) l )
  {
    XLAL_ERROR( XLAL_EINVAL );
  }

  /* For some reason GSL will not take negative m */
  /* We will have to use the relation between sph harmonics of +ve and -ve m */
  legendre = XLALAssociatedLegendreXIsZero( l, absM );
  if ( XLAL_IS_REAL8_FAIL_NAN( legendre ))
  {
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Compute the values for the spherical harmonic */
  *y = legendre * cos(m * phi);
  *y += I * legendre * sin(m * phi);

  /* If m is negative, perform some jiggery-pokery */
  if ( m < 0 && absM % 2  == 1 )
  {
    *y *= -1.0;
  }

  return XLAL_SUCCESS;
}


/**
 * Function to calculate associated Legendre function used by the spherical harmonics function
 */
static REAL8
XLALAssociatedLegendreXIsZero( const int l,
                               const int m )
{

  REAL8 legendre;

  if ( l < 0 )
  {
    XLALPrintError( "l cannot be < 0\n" );
    XLAL_ERROR_REAL8( XLAL_EINVAL );
  }
  
  if ( m < 0 || m > l )
  {
    XLALPrintError( "Invalid value of m!\n" );
    XLAL_ERROR_REAL8( XLAL_EINVAL );
  }

  /* we will switch on the values of m and n */
  switch ( l )
  {
    case 1:
      switch ( m )
      {
        case 1:
          legendre = - 1.;
          break;
        default:
          XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 2:
      switch ( m )
      {
        case 2:
          legendre = 3.;
          break;
        case 1:
          legendre = 0.;
          break;
        default:
          XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 3:
      switch ( m )
      {
        case 3:
          legendre = -15.;
          break;
        case 2:
          legendre = 0.;
          break;
        case 1:
          legendre = 1.5;
          break;
        default:
          XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 4:
      switch ( m )
      {
        case 4:
          legendre = 105.;
          break;
        case 3:
          legendre = 0.;
          break;
        case 2:
          legendre = - 7.5;
          break;
        case 1:
          legendre = 0;
          break;
        default:
          XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 5:
      switch ( m )
      {
        case 5:
          legendre = - 945.;
          break;
        case 4:
          legendre = 0.;
          break;
        case 3:
          legendre = 52.5;
          break;
        case 2:
          legendre = 0;
          break;
        case 1:
          legendre = - 1.875;
          break;
        default:
          XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 6:
      switch ( m )
      {
        case 6:
          legendre = 10395.;
          break;
        case 5:
          legendre = 0.;
          break;
        case 4:
          legendre = - 472.5;
          break;
        case 3:
          legendre = 0;
          break;
        case 2:
          legendre = 13.125;
          break;
        case 1:
          legendre = 0;
          break;
        default:
          XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 7:
      switch ( m )
      {
        case 7:
          legendre = - 135135.;
          break;
        case 6:
          legendre = 0.;
          break;
        case 5:
          legendre = 5197.5;
          break;
        case 4:
          legendre = 0.;
          break;
        case 3:
          legendre = - 118.125;
          break;
        case 2:
          legendre = 0.;
          break;
        case 1:
          legendre = 2.1875;
          break;
        default:
          XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    case 8:
      switch ( m )
      {
        case 8:
          legendre = 2027025.;
          break;
        case 7:
          legendre = 0.;
          break;
        case 6:
          legendre = - 67567.5;
          break;
        case 5:
          legendre = 0.;
          break;
        case 4:
          legendre = 1299.375;
          break;
        case 3:
          legendre = 0.;
          break;
        case 2:
          legendre = - 19.6875;
          break;
        case 1:
          legendre = 0.;
          break;
        default:
          XLAL_ERROR_REAL8( XLAL_EINVAL );
      }
      break;
    default:
      XLALPrintError( "Unsupported (l, m): %d, %d\n", l, m );
      XLAL_ERROR_REAL8( XLAL_EINVAL );
  }

  legendre *= sqrt( (REAL8)(2*l+1)*gsl_sf_fact( l-m ) / (4.*LAL_PI*gsl_sf_fact(l+m)));

  return legendre;
}

/**
 * Function to calculate the numerical prefix in the Newtonian amplitude. Eqs. 5 - 7.
 */
static int
CalculateThisMultipolePrefix(
                 COMPLEX16 *prefix, /**<< OUTPUT, Prefix value */
                 const REAL8 m1,    /**<< mass 1 */
                 const REAL8 m2,    /**<< mass 2 */
                 const INT4 l,      /**<< Mode l */
                 const INT4 m       /**<< Mode m */
                 )
{
   COMPLEX16 n;
   REAL8 c;

   REAL8 x1, x2; /* Scaled versions of component masses */

   REAL8 mult1, mult2;

   REAL8 totalMass;
   REAL8 eta;

   INT4 epsilon;
   INT4 sign; /* To give the sign of some additive terms */


   n = 0.0;

   totalMass = m1 + m2;
 
   epsilon = ( l + m )  % 2;
 
   x1 = m1 / totalMass;
   x2 = m2 / totalMass;

   eta = m1*m2/(totalMass*totalMass);

   if  ( abs( m % 2 ) == 0 )
   {
     sign = 1;
   }
   else
   {
     sign = -1;
   }
   /*
    * Eq. 7 of Damour, Iyer and Nagar 2008. 
    * For odd m, c is proportional to dM = m1-m2. In the equal-mass case, c = dM = 0. 
    * In the equal-mass unequal-spin case, however, when spins are different, the odd m term is generally not zero.
    * In this case, c can be written as c0 * dM, while spins terms in PN expansion may take the form chiA/dM.
    * Although the dM's cancel analytically, we can not implement c and chiA/dM with the possibility of dM -> 0.
    * Therefore, for this case, we give numerical values of c0 for relevant modes, and c0 is calculated as
    * c / dM in the limit of dM -> 0. Consistently, for this case, we implement chiA instead of chiA/dM
    * in LALSimIMRSpinEOBFactorizedWaveform.c. 
    */
   if  ( m1 != m2 || sign == 1 )
   {
     c = pow( x2, l + epsilon - 1 ) + sign * pow(x1, l + epsilon - 1 );
   }
   else
   {
     switch( l )
     {
       case 2:
         c = -1.0;
         break;
       case 3:
         c = -1.0;
         break;
       case 4:
         c = -0.5;
         break;
       default:
         c = 0.0;
         break;
     }
   }

   /* Eqs 5 and 6. Dependent on the value of epsilon (parity), we get different n */
   if ( epsilon == 0 )
   {
  
     n = I * m;
     n = cpow( n, (REAL8)l );
  
     mult1 = 8.0 * LAL_PI / gsl_sf_doublefact(2u*l + 1u);
     mult2 = (REAL8)((l+1) * (l+2)) / (REAL8)(l * ((INT4)l - 1));
     mult2 = sqrt(mult2);

     n *= mult1;
     n *= mult2;
  }
  else if ( epsilon == 1 )
  {
  
     n = I * m;
     n = cpow( n, (REAL8)l );
     n = -n;

     mult1 = 16.*LAL_PI / gsl_sf_doublefact( 2u*l + 1u );

     mult2  = (REAL8)( (2*l + 1) * (l+2) * (l*l - m*m) );
     mult2 /= (REAL8)( (2*l - 1) * (l+1) * l * (l-1) );
     mult2  = sqrt(mult2);

     n *= I * mult1;
     n *= mult2;
  }
  else
  {
    XLALPrintError( "Epsilon must be 0 or 1.\n");
    XLAL_ERROR( XLAL_EINVAL );
  }

  *prefix = n * eta * c;

  return XLAL_SUCCESS;
}

#endif /*_LALSIMIMRNEWTONIANMULTIPOLE_C*/
