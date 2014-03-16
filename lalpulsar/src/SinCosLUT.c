//
// Copyright (C) 2013 Karl Wette
// Copyright (C) 2009, 2010, 2011, 2012, 2013 Bernd Machenschalk
// Copyright (C) 2005, 2009 Reinhard Prix
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA  02111-1307  USA
//

#include <math.h>

#include <lal/CWFastMath.h>

#define OOTWOPI         (1.0 / LAL_TWOPI)      // 1/2pi

// main definition of lookup table code
#include "SinCosLUT.i"

int
XLALSinCosLUT(
  REAL4 *sinx,
  REAL4 *cosx,
  REAL8 x
  )
{
  return XLALSinCos2PiLUT ( sinx, cosx, x * OOTWOPI );
} // XLALSinCosLUT

int
XLALSinCos2PiLUT(
  REAL4 *sin2pix,
  REAL4 *cos2pix,
  REAL8 x
  )
{

  /* the first time we get called, we set up the lookup-table */
  static BOOLEAN firstCall = 1;
  if ( firstCall )
  {
    local_sin_cos_2PI_LUT_init();
    firstCall = 0;
  }

  /* trim the value x to interval [0..2) */
  REAL8 xt;
  SINCOS_TRIM_X(xt,x);

  /* call lookup table to calculate sin(2*pi*xt) and cos(2*pi*xt) */
  return local_sin_cos_2PI_LUT_trimmed ( sin2pix, cos2pix, xt );

} // XLALSinCos2PiLUT
