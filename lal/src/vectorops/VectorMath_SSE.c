//
// Copyright (C) 2015 Reinhard Prix
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

// ---------- INCLUDES ----------
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <config.h>

#include <lal/LALConstants.h>

#define IN_VECTORMATH
#include <lal/VectorMath.h>

#ifdef HAVE_SSE

#ifdef HAVE_SSE2
#define USE_SSE2
#endif

#include "VectorMath_sse_mathfun.h"

// ---------- local prototypes ----------
static int XLALVectorFuncf_SSE ( REAL4 *out, const REAL4 *in, UINT4 length, v4sf (*f)(v4sf x) );


// ========== function definitions ==========
static int
XLALVectorFuncf_SSE ( REAL4 *out, const REAL4 *in, UINT4 length, v4sf (*f)(v4sf x) )
{
  XLAL_CHECK ( isMemAligned(out, 16) && isMemAligned(in, 16), XLAL_EINVAL, "All vectors need to be 16-byte aligned (4xfloats\n");

  // walk through vector in blocks of 4
  UINT4 N = length;
  UINT4 i4Max = N - ( N % 4 );
  for ( UINT4 i4 = 0; i4 < i4Max/4; i4 ++ )
    {
      UINT4 i0 = i4 * 4;

      const v4sf *x4p = (const void *)&in[i0];
      v4sf *ret4p = (void *)&out[i0];

      (*ret4p) = (*f)( (*x4p) );

    } // for i4 in 0 ... i4Max/4-1

  // deal with the remaining (<=3) terms separately
  V4SF x4 = {.f = {0,0,0,0}}, ret4;
  for ( UINT4 i = i4Max, j=0; i < N; i ++, j++ ) {
    x4.f[j] = in[i];
  }
  ret4.v = (*f)( x4.v );
  for ( UINT4 i = i4Max,j=0; i < N; i ++, j++ ) {
    out[i] = ret4.f[j];
  }

  return XLAL_SUCCESS;

} // XLALVectorFuncf_SSE()

static int
XLALVectorFuncf2OUT_SSE ( REAL4 *out, REAL4 *out2, const REAL4 *x, UINT4 length, void (*f)( v4sf, v4sf*, v4sf*) )
{
  XLAL_CHECK ( isMemAligned(x,16) && isMemAligned(out,16) && isMemAligned(out2,16), XLAL_EINVAL, "All vectors need to be 16-byte aligned (4xfloats\n");

  // walk through vector in blocks of 4
  UINT4 N = length;
  UINT4 i4Max = N - ( N % 4 );
  for ( UINT4 i4 = 0; i4 < i4Max/4; i4 ++ )
    {
      UINT4 i0 = i4 * 4;

      const v4sf *x4p = (const void *)&(x[i0]);
      v4sf *out4p   = (void *)&(out[i0]);
      v4sf *out4p_2 = (void *)&(out2[i0]);

      (*f) ( (*x4p), out4p, out4p_2 );

    } // for i4 in 0 ... i4Max/4-1

  // deal with the remaining (<=3) terms separately
  V4SF x4 = {.f={0,0,0,0}}, out4, out4_2;
  for ( UINT4 i = i4Max,j=0; i < N; i ++, j++ ) {
    x4.f[j] = x[i];
  }
  (*f) ( x4.v, &out4.v, &out4_2.v );
  for ( UINT4 i = i4Max,j=0; i < N; i ++, j++ ) {
    out[i]  = out4.f[j];
    out2[i] = out4_2.f[j];
  }

  return XLAL_SUCCESS;

} // XLALVectorFuncf2OUT_SSE()

// convenience wrappers for 1-output specific functions
int XLALVectorSinf_SSE ( REAL4 *out, const REAL4 *in, UINT4 length )
{
  return XLALVectorFuncf_SSE ( out, in, length, sin_ps );
}
int XLALVectorCosf_SSE ( REAL4 *out, const REAL4 *in, UINT4 length )
{
  return XLALVectorFuncf_SSE ( out, in, length, cos_ps );
}
int XLALVectorExpf_SSE ( REAL4 *out, const REAL4 *in, UINT4 length )
{
  return XLALVectorFuncf_SSE ( out, in, length, exp_ps );
}
int XLALVectorLogf_SSE ( REAL4 *out, const REAL4 *in, UINT4 length )
{
  return XLALVectorFuncf_SSE ( out, in, length, log_ps );
}

// convenience wrappers for 2-output specific functions
int XLALVectorSinCosf_SSE ( REAL4 *sinx, REAL4 *cosx, const REAL4 *x, UINT4 length )
{
  return XLALVectorFuncf2OUT_SSE ( sinx, cosx, x, length, sincos_ps );
}
int XLALVectorSinCosf2PI_SSE ( REAL4 *sin2pix, REAL4 *cos2pix, const REAL4 *x, UINT4 length )
{
  return XLALVectorFuncf2OUT_SSE ( sin2pix, cos2pix, x, length, sincos_ps_2pi );
}


#else
int XLALVectorSinf_SSE ( REAL4 *out, const REAL4 *in, UINT4 length ){
  XLAL_CHECK ( out != NULL && in != NULL && length > 0, XLAL_EINVAL );
  XLAL_ERROR ( XLAL_EFAILED, "SSE support not compiled in\n");
}
int XLALVectorCosf_SSE ( REAL4 *out, const REAL4 *in, UINT4 length ){
  XLAL_CHECK ( out != NULL && in != NULL && length > 0, XLAL_EINVAL );
  XLAL_ERROR ( XLAL_EFAILED, "SSE support not compiled in\n");
}
int XLALVectorExpf_SSE ( REAL4 *out, const REAL4 *in, UINT4 length ){
  XLAL_CHECK ( out != NULL && in != NULL && length > 0, XLAL_EINVAL );
  XLAL_ERROR ( XLAL_EFAILED, "SSE support not compiled in\n");
}
int XLALVectorLogf_SSE ( REAL4 *out, const REAL4 *in, UINT4 length ){
  XLAL_CHECK ( out != NULL && in != NULL && length > 0, XLAL_EINVAL );
  XLAL_ERROR ( XLAL_EFAILED, "SSE support not compiled in\n");
}
int XLALVectorSinCosf_SSE ( REAL4 *out, REAL4 *out2, const REAL4 *in, UINT4 length ){
  XLAL_CHECK ( out != NULL && out2 != NULL && in != NULL && length > 0, XLAL_EINVAL );
  XLAL_ERROR ( XLAL_EFAILED, "SSE support not compiled in\n");
}
int XLALVectorSinCosf2PI_SSE ( REAL4 *out, REAL4 *out2, const REAL4 *in, UINT4 length ){
  XLAL_CHECK ( out != NULL && out2 != NULL && in != NULL && length > 0, XLAL_EINVAL );
  XLAL_ERROR ( XLAL_EFAILED, "SSE support not compiled in\n");
}
#endif // HAVE_SSE
