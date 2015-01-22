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


int
XLALVectorSinCosf_SSE ( REAL4Vector *sinx, REAL4Vector *cosx, const REAL4Vector *x )
{
  XLAL_CHECK ( (sinx != NULL) && (cosx != NULL) && (x != NULL), XLAL_EINVAL );
  XLAL_CHECK ( (sinx->length == cosx->length) && ( x->length == sinx->length), XLAL_EINVAL );

  XLAL_CHECK ( isMemAligned(x->data,16) && isMemAligned(sinx->data,16) && isMemAligned(cosx->data,16), XLAL_EINVAL, "All vectors need to be 16-byte aligned (4xfloats\n");

  // walk through vector in blocks of 4
  UINT4 N = x->length;
  UINT4 i4Max = N - ( N % 4 );
  for ( UINT4 i4 = 0; i4 < i4Max/4; i4 ++ )
    {
      UINT4 i0 = i4 * 4;

      __m128 *x4p = (void*)&x->data[i0];
      __m128 *sin4p = (void*)&sinx->data[i0];
      __m128 *cos4p = (void*)&cosx->data[i0];

      sincos_ps ( (*x4p), sin4p, cos4p );

    } // for i4 in 0 ... i4Max/4-1

  // deal with the remaining (<=3) terms separately
  V4SF x4 = {.f = {0,0,0,0}}, sin4, cos4;
  for ( UINT4 i = i4Max,j=0; i < N; i ++, j++ ) {
    x4.f[j] = x->data[i];
  }
  sincos_ps ( x4.v, &sin4.v, &cos4.v );
  for ( UINT4 i = i4Max,j=0; i < N; i ++, j++ ) {
    sinx->data[i] = sin4.f[j];
    cosx->data[i] = cos4.f[j];
  }

  return XLAL_SUCCESS;

} // XLALVectorSinCosf_SSE()

int
XLALVectorSinCosf2PI_SSE ( REAL4Vector *sin2pix, REAL4Vector *cos2pix, const REAL4Vector *x )
{
  XLAL_CHECK ( (sin2pix != NULL) && (cos2pix != NULL) && (x != NULL), XLAL_EINVAL );
  XLAL_CHECK ( (sin2pix->length == cos2pix->length) && ( x->length == sin2pix->length), XLAL_EINVAL );

  XLAL_CHECK ( isMemAligned(x->data,16) && isMemAligned(sin2pix->data,16) && isMemAligned(cos2pix->data,16), XLAL_EINVAL, "All vectors need to be 16-byte aligned (4xfloats\n");

  // walk through vector in blocks of 4
  UINT4 N = x->length;
  UINT4 i4Max = N - ( N % 4 );
  for ( UINT4 i4 = 0; i4 < i4Max/4; i4 ++ )
    {
      UINT4 i0 = i4 * 4;

      __m128 *x4p = (void*)&x->data[i0];
      __m128 *sin4p = (void*)&sin2pix->data[i0];
      __m128 *cos4p = (void*)&cos2pix->data[i0];

      sincos_ps_2pi ( (*x4p), sin4p, cos4p );

    } // for i4 in 0 ... i4Max/4-1

  // deal with the remaining (<=3) terms separately
  V4SF x4 = {.f = {0,0,0,0}}, sin4, cos4;
  for ( UINT4 i = i4Max,j=0; i < N; i ++, j++ ) {
    x4.f[j] = x->data[i];
  }
  sincos_ps_2pi ( x4.v, &sin4.v, &cos4.v );
  for ( UINT4 i = i4Max,j=0; i < N; i ++, j++ ) {
    sin2pix->data[i] = sin4.f[j];
    cos2pix->data[i] = cos4.f[j];
  }

  return XLAL_SUCCESS;

} // XLALVectorSinCosf2PI_SSE()


// convenience wrappers for specific functions
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


#endif // HAVE_SSE
