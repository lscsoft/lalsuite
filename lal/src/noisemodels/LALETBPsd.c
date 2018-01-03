/*
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
*
*/

#include <lal/LALNoiseModels.h>

/**
 * \author Craig Robinson
 * \ingroup LALNoiseModels_h
 * \brief Function to calculate the noise power spectral density of the projected ET-B detector.
 *
 * Fit taken from a Matlab script by T. Dent which can be found at:
 * https://workarea.et-gw.eu/et/WG4-Astrophysics/base-sensitivity/
 */
REAL8 XLALETBPsd( REAL8 f )
{

  /* Constants for calculating the fit */
  const REAL8 c1 = 2.39e-27;
  const REAL8 c2 = 0.349;
  const REAL8 c3 = 1.76;
  const REAL8 c4 = 0.409;

  const REAL8 p1 = -15.64;
  const REAL8 p2 = -2.145;
  const REAL8 p3 = -0.12;
  const REAL8 p4 = 1.10;

  REAL8 xt;
  REAL8 psd;

#ifndef LAL_NDEBUG
  if ( f <= 0 )
    XLAL_ERROR_REAL8( XLAL_EINVAL );
#endif

  xt = f / 100.;

  psd = c1 * pow( xt, p1 )
      + c2 * pow( xt, p2 )
      + c3 * pow( xt, p3 )
      + c4 * pow( xt, p4 );

  return 1.0e-50 * psd * psd;

}
