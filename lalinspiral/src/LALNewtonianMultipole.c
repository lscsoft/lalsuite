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
 *
 * In addition to the function used to do this, 
 * XLALCalculateNewtonianMultipole(), this file also contains a function
 * for calculating the standard scalar spherical harmonics Ylm.
 */

#include <lal/LALInspiral.h>
#include <lal/LALComplex.h>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>

int
XLALCalculateNewtonianMultipole(
                            COMPLEX16 *multipole,
                            REAL8 x,
                            REAL8 phi,
                            UINT4  l,
                            INT4  m,
                            InspiralDerivativesIn *ak
                            )
{
   static const char func[] = "XLALCalculateNewtonianMultipole";

   INT4 xlalStatus;

   COMPLEX16 n;
   REAL8 c;
   COMPLEX16 y;

   REAL8 x1, x2; /* Scaled versions of component masses */

   REAL8 mult1, mult2;

   INT4 epsilon;
   INT4 sign; /* To give the sign of some additive terms */

   n.re = n.im = 0;
   y.re = y.im = 0;

   epsilon = ( l + m )  % 2;

   x1 = ak->coeffs->m1 / ak->coeffs->totalmass;
   x2 = ak->coeffs->m2 / ak->coeffs->totalmass;

   if  ( abs( m % 2 ) == 0 )
   {
     sign = 1;
   }
   else
   {
     sign = -1;
   }

   c = pow( x2, l + epsilon - 1 ) + sign * pow(x1, l + epsilon - 1 );

   /* Dependent on the value of epsilon, we get different n */
   if ( epsilon == 0 )
   {
     
     n.im = m;
     n = XLALCOMPLEX16PowReal( n, (REAL8)l );
     
     mult1 = 8.0 * LAL_PI / gsl_sf_doublefact(2u*l + 1u);
     mult2 = (REAL8)((l+1) * (l+2)) / (REAL8)(l * ((INT4)l - 1));
     mult2 = sqrt(mult2);

     n = XLALCOMPLEX16MulReal( n, mult1 );
     n = XLALCOMPLEX16MulReal( n, mult2 );
  }
  else if ( epsilon == 1 )
  {
     
     n.im = - m;
     n = XLALCOMPLEX16PowReal( n, (REAL8)l );

     mult1 = 16.*LAL_PI / gsl_sf_doublefact( 2u*l + 1u );

     mult2  = (REAL8)( (2*l + 1) * (l+2) * (l*l + m*m) );
     mult2 /= (REAL8)( (2*l - 1) * (l+1) * l * (l+1) );
     mult2  = sqrt(mult2);

     n = XLALCOMPLEX16MulImag( n, mult1 );
     n = XLALCOMPLEX16MulReal( n, mult2 );
  }
  else
  {
    XLALPrintError( "Epsilon must be 0 or 1.\n");
    XLAL_ERROR( func, XLAL_EINVAL );
  }

  /* Calculate the necessary Ylm */
  xlalStatus = XLALScalarSphericalHarmonic( &y, l - epsilon, - m, LAL_PI_2, phi );
  if (xlalStatus != XLAL_SUCCESS )
  {
    XLAL_ERROR( func, XLAL_EFUNC );
  }

  /* Now we can construct the final answer */
  *multipole = XLALCOMPLEX16MulReal( n, c*pow( x, (REAL8)(l+epsilon)/2.0) );
  *multipole = XLALCOMPLEX16Mul( *multipole, y );

  return XLAL_SUCCESS;
}

int
XLALScalarSphericalHarmonic(
                         COMPLEX16 *y,
                         UINT4 l,
                         INT4  m,
                         REAL8 theta,
                         REAL8 phi)
{

  static const char func[] = "XLALScalarSphericalHarmonic";

  int   gslStatus;
  gsl_sf_result pLm;

  INT4 absM = abs( m );

  if ( absM > (INT4) l )
  {
    XLAL_ERROR( func, XLAL_EINVAL );
  }

  /* For some reason GSL will not take negative m */
  /* We will have to use the relation between sph harmonics of +ve and -ve m */
  XLAL_CALLGSL( gslStatus = gsl_sf_legendre_sphPlm_e((INT4)l, absM, cos(theta), &pLm ) );
  if (gslStatus != GSL_SUCCESS)
  {
    XLALPrintError("Error in GSL function\n" );
    XLAL_ERROR( func, XLAL_EFUNC );
  }

  /* Compute the values for the spherical harmonic */
  y->re = pLm.val * cos(m * phi);
  y->im = pLm.val * sin(m * phi);

  /* If m is negative, perform some jiggery-pokery */
  if ( m < 0 && absM % 2  == 1 )
  {
    y->re = - y->re;
    y->im = - y->im;
  }

  return XLAL_SUCCESS;
}
