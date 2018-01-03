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

#include <lal/LALInspiral.h>
#include <lal/LALEOBNRv2Waveform.h>

REAL8 XLALInspiralFactorizedFlux(
                      REAL8Vector           *values,
                      const REAL8           omega,
                      EOBParams             *ak,
                      const INT4             lMax
                     )

{

  REAL8 flux = 0.0;
  REAL8 v;
  REAL8 omegaSq;
  COMPLEX16 hLM;
  INT4 l, m;

  EOBNonQCCoeffs *nqcCoeffs;
  COMPLEX16       hNQC;

#ifndef LAL_NDEBUG
  if ( !values || !ak )
  {
    XLAL_ERROR_REAL8( XLAL_EFAULT );
  }
#endif

  if ( lMax < 2 )
  {
    XLAL_ERROR_REAL8( XLAL_EINVAL );
  }

  nqcCoeffs = ak->nqcCoeffs;

  /* Omegs is the derivative of phi */
  omegaSq = omega*omega;

  v = cbrt( omega );

  /* We need to apply the NQC for the (2,2) mode */
  /* To avoid having an if statement in the loop we will */
  /* deal with (2,2) and (2,1) separately */
  /* (2,2) */
  l = 2;
  m = 2;

  if ( XLALEOBNonQCCorrection( &hNQC, values, omega, nqcCoeffs ) == XLAL_FAILURE )
  {
    XLAL_ERROR_REAL8( XLAL_EFUNC );
  }

  if ( XLALGetFactorizedWaveform( &hLM, values, v, l, m, ak )
           == XLAL_FAILURE )
  {
    XLAL_ERROR_REAL8( XLAL_EFUNC );
  }
  /* For the 2,2 mode, we apply NQC correction to the flux */
  hLM = hNQC * hLM;

  flux = (REAL8)(m * m) * omegaSq * cabs( hLM ) * cabs( hLM );

  /* (2,1) */
  l = 2;
  m = 1;

  if ( XLALGetFactorizedWaveform( &hLM, values, v, l, m, ak )
           == XLAL_FAILURE )
  {
    XLAL_ERROR_REAL8( XLAL_EFUNC );
  }

  flux += (REAL8)(m * m) * omegaSq * cabs( hLM ) * cabs( hLM );

  /* All other modes */
  for ( l = 3; l <= lMax; l++ )
  {
    /*INT4 minM = l-3;
    if ( minM < 1 )
      minM = 1;*/

    for ( m = 1; m <= l; m++ )
    {

      if ( XLALGetFactorizedWaveform( &hLM, values, v, l, m, ak )
             == XLAL_FAILURE )
      {
        XLAL_ERROR_REAL8( XLAL_EFUNC );
      }

      flux += (REAL8)(m * m) * omegaSq * cabs( hLM ) * cabs( hLM );
    }
  }

  return flux * LAL_1_PI / 8.0;
}
