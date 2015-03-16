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
static int XLALVectorFuncf_SSE ( REAL4Vector *out, const REAL4Vector *in, v4sf (*f)(v4sf x) );


// ========== function definitions ==========
static int
XLALVectorFuncf_SSE ( REAL4Vector *out, const REAL4Vector *in, v4sf (*f)(v4sf x) )
{
  XLAL_CHECK ( (out != NULL) && (in != NULL) && (f != NULL), XLAL_EINVAL );
  XLAL_CHECK ( in->length == out->length, XLAL_EINVAL );

  XLAL_CHECK ( isMemAligned(out->data,16) && isMemAligned(in->data,16), XLAL_EINVAL, "All vectors need to be 16-byte aligned (4xfloats\n");

  // walk through vector in blocks of 4
  UINT4 N = in->length;
  UINT4 i4Max = N - ( N % 4 );
  for ( UINT4 i4 = 0; i4 < i4Max/4; i4 ++ )
    {
      UINT4 i0 = i4 * 4;

      __m128 *x4p = (void*)&in->data[i0];
      __m128 *ret4p = (void*)&out->data[i0];

      (*ret4p) = (*f)( (*x4p) );

    } // for i4 in 0 ... i4Max/4-1

  // deal with the remaining (<=3) terms separately
  V4SF x4 = {.f = {0,0,0,0}}, ret4;
  for ( UINT4 i = i4Max, j=0; i < N; i ++, j++ ) {
    x4.f[j] = in->data[i];
  }
  ret4.v = (*f)( x4.v );
  for ( UINT4 i = i4Max,j=0; i < N; i ++, j++ ) {
    out->data[i] = ret4.f[j];
  }

  return XLAL_SUCCESS;

} // XLALVectorFuncf_SSE()

static int
XLALVectorFuncf2OUT_SSE ( REAL4Vector *out, REAL4Vector *out2, const REAL4Vector *x, void (*f)( v4sf, v4sf*, v4sf*) )
{
  XLAL_CHECK ( (out != NULL) && (out2 != NULL) && (x != NULL), XLAL_EINVAL );
  XLAL_CHECK ( (out->length == out2->length) && ( x->length == out->length), XLAL_EINVAL );

  XLAL_CHECK ( isMemAligned(x->data,16) && isMemAligned(out->data,16) && isMemAligned(out2->data,16), XLAL_EINVAL, "All vectors need to be 16-byte aligned (4xfloats\n");

  // walk through vector in blocks of 4
  UINT4 N = x->length;
  UINT4 i4Max = N - ( N % 4 );
  for ( UINT4 i4 = 0; i4 < i4Max/4; i4 ++ )
    {
      UINT4 i0 = i4 * 4;

      __m128 *x4p     = (void*)&(x->data[i0]);
      __m128 *out4p   = (void*)&(out->data[i0]);
      __m128 *out4p_2 = (void*)&(out2->data[i0]);

      (*f) ( (*x4p), out4p, out4p_2 );

    } // for i4 in 0 ... i4Max/4-1

  // deal with the remaining (<=3) terms separately
  V4SF x4 = {.f={0,0,0,0}}, out4, out4_2;
  for ( UINT4 i = i4Max,j=0; i < N; i ++, j++ ) {
    x4.f[j] = x->data[i];
  }
  (*f) ( x4.v, &out4.v, &out4_2.v );
  for ( UINT4 i = i4Max,j=0; i < N; i ++, j++ ) {
    out->data[i]  = out4.f[j];
    out2->data[i] = out4_2.f[j];
  }

  return XLAL_SUCCESS;

} // XLALVectorFuncf2OUT_SSE()

// convenience wrappers for 1-output specific functions
int XLALVectorSinf_SSE ( REAL4Vector *out, const REAL4Vector *in )
{
  return XLALVectorFuncf_SSE ( out, in, sin_ps );
}
int XLALVectorCosf_SSE ( REAL4Vector *out, const REAL4Vector *in )
{
  return XLALVectorFuncf_SSE ( out, in, cos_ps );
}
int XLALVectorExpf_SSE ( REAL4Vector *out, const REAL4Vector *in )
{
  return XLALVectorFuncf_SSE ( out, in, exp_ps );
}
int XLALVectorLogf_SSE ( REAL4Vector *out, const REAL4Vector *in )
{
  return XLALVectorFuncf_SSE ( out, in, log_ps );
}

// convenience wrappers for 2-output specific functions
int XLALVectorSinCosf_SSE ( REAL4Vector *sinx, REAL4Vector *cosx, const REAL4Vector *x )
{
  return XLALVectorFuncf2OUT_SSE ( sinx, cosx, x, sincos_ps );
}
int XLALVectorSinCosf2PI_SSE ( REAL4Vector *sin2pix, REAL4Vector *cos2pix, const REAL4Vector *x )
{
  return XLALVectorFuncf2OUT_SSE ( sin2pix, cos2pix, x, sincos_ps_2pi );
}


#else
int XLALVectorSinf_SSE ( REAL4Vector *out, const REAL4Vector *in ){
  XLAL_CHECK ( (out != NULL) && (in != NULL), XLAL_EINVAL );
  XLAL_ERROR ( XLAL_EFAILED, "SSE support not compiled in\n");
}
int XLALVectorCosf_SSE ( REAL4Vector *out, const REAL4Vector *in ){
  XLAL_CHECK ( (out != NULL) && (in != NULL), XLAL_EINVAL );
  XLAL_ERROR ( XLAL_EFAILED, "SSE support not compiled in\n");
}
int XLALVectorExpf_SSE ( REAL4Vector *out, const REAL4Vector *in ){
  XLAL_CHECK ( (out != NULL) && (in != NULL), XLAL_EINVAL );
  XLAL_ERROR ( XLAL_EFAILED, "SSE support not compiled in\n");
}
int XLALVectorLogf_SSE ( REAL4Vector *out, const REAL4Vector *in ){
  XLAL_CHECK ( (out != NULL) && (in != NULL), XLAL_EINVAL );
  XLAL_ERROR ( XLAL_EFAILED, "SSE support not compiled in\n");
}
int XLALVectorSinCosf_SSE ( REAL4Vector *out, REAL4Vector *out2, const REAL4Vector *in ){
  XLAL_CHECK ( (out != NULL) && (out2 != NULL) && (in != NULL), XLAL_EINVAL );
  XLAL_ERROR ( XLAL_EFAILED, "SSE support not compiled in\n");
}
int XLALVectorSinCosf2PI_SSE ( REAL4Vector *out, REAL4Vector *out2, const REAL4Vector *in ){
  XLAL_CHECK ( (out != NULL) && (out2 != NULL) && (in != NULL), XLAL_EINVAL );
  XLAL_ERROR ( XLAL_EFAILED, "SSE support not compiled in\n");
}
#endif // HAVE_SSE
