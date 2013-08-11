//
// Copyright (C) 2013 Karl Wette
// Copyright (C) 2005 Reinhard Prix
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
#include <lal/LALConstants.h>

#define OOTWOPI         (1.0 / LAL_TWOPI)	// 1/2pi

#define LUT_RES         64      /* resolution of lookup-table */
#define LUT_RES_F	(1.0 * LUT_RES)
#define OO_LUT_RES	(1.0 / LUT_RES)

#define X_TO_IND	(1.0 * LUT_RES * OOTWOPI )
#define IND_TO_X	(LAL_TWOPI * OO_LUT_RES)

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
  REAL8 xt;
  INT4 i0;
  REAL8 d, d2;
  REAL8 ts, tc;
  REAL8 dummy;

  static BOOLEAN firstCall = 1;
  static REAL4 sinVal[LUT_RES+1], cosVal[LUT_RES+1];

  /* the first time we get called, we set up the lookup-table */
  if ( firstCall )
  {
    UINT4 k;
    for (k=0; k <= LUT_RES; k++)
    {
      sinVal[k] = sin( LAL_TWOPI * k * OO_LUT_RES );
      cosVal[k] = cos( LAL_TWOPI * k * OO_LUT_RES );
    }
    firstCall = 0;
  }

  /* we only need the fractional part of 'x', which is number of cylces,
   * this was previously done using
   *   xt = x - (INT4)x;
   * which is numerically unsafe for x > LAL_INT4_MAX ~ 2e9
   * for saftey we therefore rather use modf(), even if that
   * will be somewhat slower...
   */
  xt = modf(x, &dummy);/* xt in (-1, 1) */

  if ( xt < 0.0 )
    xt += 1.0;			/* xt in [0, 1 ) */
#ifndef LAL_NDEBUG
  if ( xt < 0.0 || xt > 1.0 )
  {
    XLALPrintError("\nFailed numerica in sin_cos_2PI_LUT(): xt = %f not in [0,1)\n\n", xt );
    return XLAL_FAILURE;
  }
#endif

  i0 = (INT4)( xt * LUT_RES_F + 0.5 );	/* i0 in [0, LUT_RES ] */
  d = d2 = LAL_TWOPI * (xt - OO_LUT_RES * i0);
  d2 *= 0.5 * d;

  ts = sinVal[i0];
  tc = cosVal[i0];

  /* use Taylor-expansions for sin/cos around LUT-points */
  (*sin2pix) = ts + d * tc - d2 * ts;
  (*cos2pix) = tc - d * ts - d2 * tc;

  return XLAL_SUCCESS;

} // XLALSinCos2PiLUT
