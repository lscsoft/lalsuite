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

#include <lal/LALComplex.h>
#include <lal/LALInspiral.h>

REAL8 XLALInspiralFactorizedFlux(
                      REAL8Vector           *values,
                      REAL8Vector           *dvalues,
                      InspiralDerivativesIn *ak,
                      const INT4             lMax
                     )

{

  static const char func[] = "XLALInspiralFactorizedFlux";

  REAL8 flux = 0.0;
  REAL8 omegaSq;
  COMPLEX16 hLM;
  INT4 l, m;

#ifndef LAL_NDEBUG
  if ( !values || !dvalues || !ak )
  {
    XLAL_ERROR_REAL8( func, XLAL_EFAULT );
  }
#endif

  if ( lMax < 2 )
  {
    XLAL_ERROR_REAL8( func, XLAL_EINVAL );
  }

  /* Omegs is the derivative of phi */
  omegaSq = dvalues->data[1];
  omegaSq *= omegaSq;

  for ( l = 2; l <= lMax; l++ )
  {
    for ( m = 1; m <= l; m++ )
    {

      if ( XLALGetFactorizedWaveform( &hLM, values, dvalues, ak, l, m )
             == XLAL_FAILURE )
      {
        XLAL_ERROR_REAL8( func, XLAL_EFUNC );
      }

      flux += (REAL8)(m * m) * omegaSq * XLALCOMPLEX16Abs2( hLM );
    }
  }

  return flux * LAL_1_PI / 8.0;
}
