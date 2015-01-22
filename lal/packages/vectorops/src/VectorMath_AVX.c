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

#ifdef HAVE_AVX

#include "VectorMath_avx_mathfun.h"

// ---------- local prototypes ----------
static int XLALVectorFuncf_AVX ( REAL4Vector *out, const REAL4Vector *in, v8sf (*f)(v8sf x) );

// ========== function definitions ==========
static int
XLALVectorFuncf_AVX ( REAL4Vector *out, const REAL4Vector *in, v8sf (*f)(v8sf x) )
{
  XLAL_CHECK ( (out != NULL) && (in != NULL) && (f != NULL), XLAL_EINVAL );
  XLAL_CHECK ( in->length == out->length, XLAL_EINVAL );

  XLAL_CHECK ( isMemAligned(out->data,32) && isMemAligned(in->data,32), XLAL_EINVAL, "All vectors need to be 32-byte aligned (8xfloats\n");

  // walk through vector in blocks of 8
  UINT4 N = in->length;
  UINT4 i8Max = N - ( N % 8 );
  for ( UINT4 i8 = 0; i8 < i8Max/8; i8 ++ )
    {
      UINT4 i0 = i8 * 8;

      __m256 *x8p   = (void*)&(in->data[i0]);
      __m256 *ret8p = (void*)&(out->data[i0]);

      (*ret8p) = (*f)( (*x8p) );

    } // for i8 in 0 ... i8Max/8-1

  // deal with the remaining (<=7) terms separately
  V8SF x8 = {.f={0,0,0,0,0,0,0,0}}, ret8;
  for ( UINT4 i = i8Max,j=0; i < N; i ++, j++ ) {
    x8.f[j] = in->data[i];
  }
  ret8.v = (*f)( x8.v );
  for ( UINT4 i = i8Max,j=0; i < N; i ++, j++ ) {
    out->data[i] = ret8.f[j];
  }

  return XLAL_SUCCESS;

} // XLALVectorFuncf_AVX()

int
XLALVectorSinCosf_AVX ( REAL4Vector *sinx, REAL4Vector *cosx, const REAL4Vector *x )
{
  XLAL_CHECK ( (sinx != NULL) && (cosx != NULL) && (x != NULL), XLAL_EINVAL );
  XLAL_CHECK ( (sinx->length == cosx->length) && ( x->length == sinx->length), XLAL_EINVAL );

  XLAL_CHECK ( isMemAligned(x->data,32) && isMemAligned(sinx->data,32) && isMemAligned(cosx->data,32), XLAL_EINVAL, "All vectors need to be 32-byte aligned (8xfloats\n");

  // walk through vector in blocks of 4
  UINT4 N = x->length;
  UINT4 i8Max = N - ( N % 8 );
  for ( UINT4 i8 = 0; i8 < i8Max/8; i8 ++ )
    {
      UINT4 i0 = i8 * 8;

      __m256 *x8p   = (void*)&(x->data[i0]);
      __m256 *sin8p = (void*)&(sinx->data[i0]);
      __m256 *cos8p = (void*)&(cosx->data[i0]);

      sincos256_ps ( (*x8p), sin8p, cos8p );

    } // for i8 in 0 ... i8Max/8-1

  // deal with the remaining (<=7) terms separately
  V8SF x8 = {.f={0,0,0,0,0,0,0,0}}, sin8, cos8;
  for ( UINT4 i = i8Max,j=0; i < N; i ++, j++ ) {
    x8.f[j] = x->data[i];
  }
  sincos256_ps ( x8.v, &sin8.v, &cos8.v );
  for ( UINT4 i = i8Max,j=0; i < N; i ++, j++ ) {
    sinx->data[i] = sin8.f[j];
    cosx->data[i] = cos8.f[j];
  }

  return XLAL_SUCCESS;

} // XLALVectorSinCosf_AVX()

int
XLALVectorSinCosf2PI_AVX ( REAL4Vector *sin2pix, REAL4Vector *cos2pix, const REAL4Vector *x )
{
  XLAL_CHECK ( (sin2pix != NULL) && (cos2pix != NULL) && (x != NULL), XLAL_EINVAL );
  XLAL_CHECK ( (sin2pix->length == cos2pix->length) && ( x->length == sin2pix->length), XLAL_EINVAL );

  XLAL_CHECK ( isMemAligned(x->data,32) && isMemAligned(sin2pix->data,32) && isMemAligned(cos2pix->data,32), XLAL_EINVAL, "All vectors need to be 32-byte aligned (8xfloats\n");

  // walk through vector in blocks of 4
  UINT4 N = x->length;
  UINT4 i8Max = N - ( N % 8 );
  for ( UINT4 i8 = 0; i8 < i8Max/8; i8 ++ )
    {
      UINT4 i0 = i8 * 8;

      __m256 *x8p   = (void*)&(x->data[i0]);
      __m256 *sin8p = (void*)&(sin2pix->data[i0]);
      __m256 *cos8p = (void*)&(cos2pix->data[i0]);

      sincos256_ps_2pi ( (*x8p), sin8p, cos8p );

    } // for i8 in 0 ... i8Max/8-1

  // deal with the remaining (<=7) terms separately
  V8SF x8 = {.f={0,0,0,0,0,0,0,0}}, sin8, cos8;
  for ( UINT4 i = i8Max,j=0; i < N; i ++, j++ ) {
    x8.f[j] = x->data[i];
  }
  sincos256_ps_2pi ( x8.v, &sin8.v, &cos8.v );
  for ( UINT4 i = i8Max,j=0; i < N; i ++, j++ ) {
    sin2pix->data[i] = sin8.f[j];
    cos2pix->data[i] = cos8.f[j];
  }

  return XLAL_SUCCESS;

} // XLALVectorSinCosf2PI_AVX()

// convenience wrappers for specific functions
int XLALVectorSinf_AVX ( REAL4Vector *out, const REAL4Vector *in )
{
  return XLALVectorFuncf_AVX ( out, in, sin256_ps );
}
int XLALVectorCosf_AVX ( REAL4Vector *out, const REAL4Vector *in )
{
  return XLALVectorFuncf_AVX ( out, in, cos256_ps );
}
int XLALVectorExpf_AVX ( REAL4Vector *out, const REAL4Vector *in )
{
  return XLALVectorFuncf_AVX ( out, in, exp256_ps );
}
int XLALVectorLogf_AVX ( REAL4Vector *out, const REAL4Vector *in )
{
  return XLALVectorFuncf_AVX ( out, in, log256_ps );
}


#endif // HAVE_AVX
