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

static BOOLEAN haveLUT = 0;

/* global VARIABLES to be used in (global) macros */
UNUSED REAL4 sincosLUTbase[SINCOS_LUT_RES+SINCOS_LUT_RES/4];
UNUSED REAL4 sincosLUTdiff[SINCOS_LUT_RES+SINCOS_LUT_RES/4];

/*
 * LUT initialization. Normally not required for user, as will
 * be called transparently by sincosLUT functions.
 * Put here for certain specialized low-level usage in hotloops.
*/
void
XLALSinCosLUTInit (void)
{
  if ( haveLUT ) {
    return;
  }

  static const REAL8 step = LAL_TWOPI / (REAL8)SINCOS_LUT_RES;
  static const REAL8 divide  = 1.0 / ( 1 << SINCOS_SHIFT );
  REAL8 start, end, true_mid, linear_mid;
  int i;

  start = 0.0; /* sin(0 * step) */
  for( i = 0; i < SINCOS_LUT_RES + SINCOS_LUT_RES/4; i++ )
    {
      true_mid = sin( ( i + 0.5 ) * step );
      end = sin( ( i + 1 ) * step );
      linear_mid = ( start + end ) * 0.5;
      sincosLUTbase[i] = start + ( ( true_mid - linear_mid ) * 0.5 );
      sincosLUTdiff[i] = ( end - start ) * divide;
      start = end;
    } // for i < LUT_RES

  haveLUT = 1;
  return;

} // XLALSinCosLUTInit()

///
/// Calculate sin(x) and cos(x) to roughly 1e-7 precision using a lookup-table and Taylor-expansion.
///
/// \note This function will fail for arguments larger than |x| > INT4_MAX = 2147483647 ~ 2e9 !!!
///
/// Returns XLAL_SUCCESS or XLAL_FAILURE.
///
int
XLALSinCosLUT ( REAL4 *sinx, REAL4 *cosx, REAL8 x )
{
  return XLALSinCos2PiLUT ( sinx, cosx, x * OOTWOPI );
} // XLALSinCosLUT

///
/// Calculate sin(2*pi*x) and cos(2*pi*x) to roughly 1e-7 precision using a lookup-table and
/// Taylor-expansion.
///
/// \note This function will fail for arguments larger than |x| > INT4_MAX = 2147483647 ~ 2e9 !!!
///
/// Returns XLAL_SUCCESS or XLAL_FAILURE.
///
int
XLALSinCos2PiLUT ( REAL4 *sin2pix, REAL4 *cos2pix, REAL8 x )
{
  /* trim the value x to interval [0..2) */
  REAL8 xt;
  SINCOS_TRIM_X(xt,x);

  /* call lookup table to calculate sin(2*pi*xt) and cos(2*pi*xt) */
  return XLALSinCos2PiLUTtrimmed ( sin2pix, cos2pix, xt );

} // XLALSinCos2PiLUT


/** A function that uses the lookup tables to evaluate sin and cos values
 * of 2*Pi*x, but relies on x being already trimmed to the interval [0..2)
 */
int XLALSinCos2PiLUTtrimmed ( REAL4 *s, REAL4 *c, REAL8 x )
{
  /* check range of input only in DEBUG mode */
#ifndef LAL_NDEBUG
  if(x > SINCOS_ADDS) {
    XLALPrintError("%s: x too large: %22f > %f\n", __func__, x, SINCOS_ADDS);
    return XLAL_FAILURE;
  } else if(x < -SINCOS_ADDS) {
    XLALPrintError("%s: x too small: %22f < %f\n", __func__, x, -SINCOS_ADDS);
    return XLAL_FAILURE;
  }
#endif

  /* the first time we get called, we set up the lookup-table */
  if ( ! haveLUT ) {
    XLALSinCosLUTInit();
  }

  /* use the macros defined above */
  SINCOS_PROLOG
  SINCOS_STEP1(x)
  SINCOS_STEP2
  SINCOS_STEP3
  SINCOS_STEP4
  SINCOS_STEP5(s)
  SINCOS_STEP6(c)
  SINCOS_EPILOG(s,c,x)

  return XLAL_SUCCESS;

} /* XLALSinCos2PiLUTtrimmed() */
