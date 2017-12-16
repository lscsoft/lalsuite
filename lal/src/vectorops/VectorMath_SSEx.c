//
// Copyright (C) 2015 Reinhard Prix, Karl Wette
// Copyright (C) 2015 Evan Goetz
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

#include "VectorMath_internal.h"

#ifdef __SSE2__
#define USE_SSE2
#endif

#include "VectorMath_sse_mathfun.h"

// ---------- local operators and operator-wrappers ----------
UNUSED static inline __m128
local_add_ps ( __m128 in1, __m128 in2 )
{
  return _mm_add_ps ( in1, in2 );
}

UNUSED static inline __m128
local_mul_ps ( __m128 in1, __m128 in2 )
{
  return _mm_mul_ps ( in1, in2 );
}

UNUSED static inline __m128
local_max_ps ( __m128 in1, __m128 in2 )
{
  return _mm_max_ps ( in1, in2 );
}

UNUSED static inline __m128d
local_add_pd ( __m128d in1, __m128d in2 )
{
  return _mm_add_pd ( in1, in2 );
}

UNUSED static inline __m128d
local_mul_pd ( __m128d in1, __m128d in2 )
{
  return _mm_mul_pd ( in1, in2 );
}

// in1: a0,b0,a1,b1, in2: c0,d0,c1,d1
UNUSED static inline __m128
local_cmul_ps ( __m128 in1, __m128 in2 )
{
  __m128 neg = _mm_setr_ps(1.0, -1.0, 1.0, -1.0);

  // Multiply in1 and in2
  // a0c0, b0d0, a1c1, b1d1
  __m128 temp1 = _mm_mul_ps(in1, in2);
  //Negate the second and fourth element of temp1
  temp1 = _mm_mul_ps(temp1, neg);

  //switch elements of in2
  //d0,c0,d1,c1
  in2 = _mm_shuffle_ps(in2,in2,0b10110001 );

  // Multiply in1 and in2
  // a0d0, b0c0, a1d1, c1b1
  __m128 temp2 =  _mm_mul_ps(in1, in2);

  //switch elements of temp1 and temp2
  //a0c0, a1c1, a0d0, a1d1
   __m128 shuf1 = _mm_shuffle_ps(temp1,temp2, 0b10001000);
  //-b0d0, -b1d1, b0c0, b1c1
   __m128 shuf2= _mm_shuffle_ps(temp1,temp2, 0b11011101);

  // Add shuf2 to shuf1
  __m128 result = _mm_add_ps(shuf1, shuf2);

  //reorder result
  return _mm_shuffle_ps(result, result,0b11011000);
}

// ========== internal generic SSEx functions ==========

// ---------- generic SSEx operator with 1 REAL4 vector input to 1 REAL4 vector output (S2S) ----------
static inline int
XLALVectorMath_S2S_SSEx ( REAL4 *out, const REAL4 *in, const UINT4 len, __m128 (*f)(__m128) )
{

  // walk through vector in blocks of 4
  UINT4 i4Max = len - ( len % 4 );
  for ( UINT4 i4 = 0; i4 < i4Max; i4 += 4 )
    {
      __m128 in4p = _mm_loadu_ps(&in[i4]);
      __m128 out4p = (*f)( in4p );
      _mm_storeu_ps(&out[i4], out4p);
    }

  // deal with the remaining (<=3) terms separately
  V4SF in4 = {.f = {0,0,0,0}}, out4;
  for ( UINT4 i = i4Max, j=0; i < len; i ++, j++ ) {
    in4.f[j] = in[i];
  }
  out4.v = (*f)( in4.v );
  for ( UINT4 i = i4Max,j=0; i < len; i ++, j++ ) {
    out[i] = out4.f[j];
  }

  return XLAL_SUCCESS;

} // XLALVectorMath_S2S_SSEx()

// ---------- generic SSEx operator with 1 REAL4 vector input to 2 REAL4 vector outputs (S2SS) ----------
static inline int
XLALVectorMath_S2SS_SSEx ( REAL4 *out1, REAL4 *out2, const REAL4 *in, const UINT4 len, void (*f)(__m128, __m128*, __m128*) )
{

  // walk through vector in blocks of 4
  UINT4 i4Max = len - ( len % 4 );
  for ( UINT4 i4 = 0; i4 < i4Max; i4 += 4 )
    {
      __m128 in4p = _mm_loadu_ps(&in[i4]);
      __m128 out4p_1, out4p_2;
      (*f) ( in4p, &out4p_1, &out4p_2 );
      _mm_storeu_ps(&out1[i4], out4p_1);
      _mm_storeu_ps(&out2[i4], out4p_2);
    }

  // deal with the remaining (<=3) terms separately
  V4SF in4 = {.f={0,0,0,0}}, out4_1, out4_2;
  for ( UINT4 i = i4Max,j=0; i < len; i ++, j++ ) {
    in4.f[j] = in[i];
  }
  (*f) ( in4.v, &out4_1.v, &out4_2.v );
  for ( UINT4 i = i4Max,j=0; i < len; i ++, j++ ) {
    out1[i] = out4_1.f[j];
    out2[i] = out4_2.f[j];
  }

  return XLAL_SUCCESS;

} // XLALVectorMath_S2SS_SSEx()

// ---------- generic SSEx operator with 2 REAL4 vector inputs to 1 REAL4 vector output (SS2S) ----------
static inline int
XLALVectorMath_SS2S_SSEx ( REAL4 *out, const REAL4 *in1, const REAL4 *in2, const UINT4 len, __m128 (*op)(__m128, __m128) )
{

  // walk through vector in blocks of 4
  UINT4 i4Max = len - ( len % 4 );
  for ( UINT4 i4 = 0; i4 < i4Max; i4 += 4 )
    {
      __m128 in4p_1 = _mm_loadu_ps(&in1[i4]);
      __m128 in4p_2 = _mm_loadu_ps(&in2[i4]);
      __m128 out4p = (*op) ( in4p_1, in4p_2 );
      _mm_storeu_ps(&out[i4], out4p);
    }

  // deal with the remaining (<=3) terms separately
  V4SF in4_1 = {.f={0,0,0,0}};
  V4SF in4_2 = {.f={0,0,0,0}};
  V4SF out4;
  for ( UINT4 i = i4Max,j=0; i < len; i ++, j++ ) {
    in4_1.f[j] = in1[i];
    in4_2.f[j] = in2[i];
  }
  out4.v = (*op) ( in4_1.v, in4_2.v );
  for ( UINT4 i = i4Max,j=0; i < len; i ++, j++ ) {
    out[i] = out4.f[j];
  }

  return XLAL_SUCCESS;

} // XLALVectorMath_SS2S_SSEx()

// ---------- generic SSEx operator with 1 REAL4 scalar and 1 REAL4 vector inputs to 1 REAL4 vector output (sS2S) ----------
static inline int
XLALVectorMath_sS2S_SSEx ( REAL4 *out, REAL4 scalar, const REAL4 *in, const UINT4 len, __m128 (*op)(__m128, __m128) )
{
  const V4SF scalar4 = {.f={scalar,scalar,scalar,scalar}};

  // walk through vector in blocks of 4
  UINT4 i4Max = len - ( len % 4 );
  for ( UINT4 i4 = 0; i4 < i4Max; i4 += 4 )
    {
      __m128 in4p = _mm_loadu_ps(&in[i4]);
      __m128 out4p = (*op) ( scalar4.v, in4p );
      _mm_storeu_ps(&out[i4], out4p);
    }

  // deal with the remaining (<=3) terms separately
  V4SF in4 = {.f={0,0,0,0}};
  V4SF out4;
  for ( UINT4 i = i4Max,j=0; i < len; i ++, j++ ) {
    in4.f[j] = in[i];
  }
  out4.v = (*op) ( scalar4.v, in4.v );
  for ( UINT4 i = i4Max,j=0; i < len; i ++, j++ ) {
    out[i] = out4.f[j];
  }

  return XLAL_SUCCESS;

} // XLALVectorMath_sS2S_SSEx()

// ---------- generic SSEx operator with 1 REAL8 scalar and 1 REAL8 vector inputs to 1 REAL8 vector output (dD2D) ----------
static inline int
XLALVectorMath_dD2D_SSEx ( REAL8 *out, REAL8 scalar, const REAL8 *in, const UINT4 len, __m128d (*op)(__m128d, __m128d) )
{
  const V2SF scalar2 = {.f={scalar,scalar}};

  // walk through vector in blocks of 2
  UINT4 i2Max = len - ( len % 2 );
  for ( UINT4 i2 = 0; i2 < i2Max; i2 += 2 )
    {
      __m128d in2p = _mm_loadu_pd(&in[i2]);
      __m128d out2p = (*op) ( scalar2.v, in2p );
      _mm_storeu_pd(&out[i2], out2p);
    }

  // deal with the remaining (<=1) terms separately
  V2SF in2 = {.f={0,0}};
  V2SF out2;
  for ( UINT4 i = i2Max,j=0; i < len; i ++, j++ ) {
    in2.f[j] = in[i];
  }
  out2.v = (*op) ( scalar2.v, in2.v );
  for ( UINT4 i = i2Max,j=0; i < len; i ++, j++ ) {
    out[i] = out2.f[j];
  }

  return XLAL_SUCCESS;

} // XLALVectorMath_dD2D_SSEx()

// ---------- generic SSEx operator with 2 REAL8 vector inputs to 1 REAL8 vector output (DD2D) ----------
static inline int
XLALVectorMath_DD2D_SSEx ( REAL8 *out, const REAL8 *in1, const REAL8 *in2, const UINT4 len, __m128d (*op)(__m128d, __m128d) )
{

  // walk through vector in blocks of 2
  UINT4 i2Max = len - ( len % 2 );
  for ( UINT4 i2 = 0; i2 < i2Max; i2 += 2 )
    {
      __m128d in2p_1 = _mm_loadu_pd(&in1[i2]);
      __m128d in2p_2 = _mm_loadu_pd(&in2[i2]);
      __m128d out2p = (*op) ( in2p_1, in2p_2 );
      _mm_storeu_pd(&out[i2], out2p);
    }

  // deal with the remaining (<=1) terms separately
  V2SF in2_1 = {.f={0,0}};
  V2SF in2_2 = {.f={0,0}};
  V2SF out2;
  for ( UINT4 i = i2Max,j=0; i < len; i ++, j++ ) {
    in2_1.f[j] = in1[i];
    in2_2.f[j] = in2[i];
  }
  out2.v = (*op) ( in2_1.v, in2_2.v );
  for ( UINT4 i = i2Max,j=0; i < len; i ++, j++ ) {
    out[i] = out2.f[j];
  }

  return XLAL_SUCCESS;

} // XLALVectorMath_DD2D_SSEx()

// ---------- generic SSEx operator with 2 COMPLEX8 vector inputs to 1 COMPLEX8 vector output (CC2C) ----------
static inline int
XLALVectorMath_CC2C_SSEx ( COMPLEX8 *out, const COMPLEX8 *in1, const COMPLEX8 *in2, const UINT4 len, __m128 (*op)(__m128, __m128) )
{

  // walk through vector in blocks of 2
  UINT4 i2Max = len - ( len % 2 );
  for ( UINT4 i2 = 0; i2 < i2Max; i2 += 2 )
    {
      __m128 in4p_1 = _mm_loadu_ps( (const REAL4*)&in1[i2] );
      __m128 in4p_2 = _mm_loadu_ps( (const REAL4*)&in2[i2] );
      __m128 out4p = (*op) ( in4p_1, in4p_2 );
      _mm_storeu_ps(( REAL4*)&out[i2], out4p);
    }

  // deal with the remaining (<=1) term separately
  V4SF in4_1 = {.f={0,0,0,0}};
  V4SF in4_2 = {.f={0,0,0,0}};
  V4SF out4;
  for ( UINT4 i = i2Max,j=0; i < len; i ++, j=+2 )
    {
      in4_1.f[j] = crealf ( in1[i] );
      in4_1.f[j+1] = cimagf ( in1[i] );
      in4_2.f[j] = crealf ( in2[i] );
      in4_2.f[j+1] = cimagf ( in2[i] );
    }
  out4.v = (*op) ( in4_1.v, in4_2.v );
  for ( UINT4 i = i2Max,j=0; i < len; i ++, j+=2 )
    {
      out[i] = crect( out4.f[j], out4.f[j+1] );
    }


  return XLAL_SUCCESS;

} // XLALVectorMath_CC2C_SSEx()

// ========== internal SSEx vector math functions ==========

// ---------- define vector math functions with 1 REAL4 vector input to 1 REAL4 vector output (S2S) ----------
#define DEFINE_VECTORMATH_S2S(NAME, SSE_OP)                             \
  DEFINE_VECTORMATH_ANY( XLALVectorMath_S2S_SSEx, NAME ## REAL4, ( REAL4 *out, const REAL4 *in, const UINT4 len ), ( (out != NULL) && (in != NULL) ), ( out, in, len, SSE_OP ) )

DEFINE_VECTORMATH_S2S(Sin, sin_ps)
DEFINE_VECTORMATH_S2S(Cos, cos_ps)
DEFINE_VECTORMATH_S2S(Exp, exp_ps)
DEFINE_VECTORMATH_S2S(Log, log_ps)

// ---------- define vector math functions with 1 REAL4 vector input to 2 REAL4 vector outputs (S2SS) ----------
#define DEFINE_VECTORMATH_S2SS(NAME, SSE_OP)                            \
  DEFINE_VECTORMATH_ANY( XLALVectorMath_S2SS_SSEx, NAME ## REAL4, ( REAL4 *out1, REAL4 *out2, const REAL4 *in, const UINT4 len ), ( (out1 != NULL) && (out2 != NULL) && (in != NULL) ), ( out1, out2, in, len, SSE_OP ) )

DEFINE_VECTORMATH_S2SS(SinCos, sincos_ps)
DEFINE_VECTORMATH_S2SS(SinCos2Pi, sincos_ps_2pi)

// ---------- define vector math functions with 2 REAL4 vector inputs to 1 REAL4 vector output (SS2S) ----------
#define DEFINE_VECTORMATH_SS2S(NAME, SSE_OP)                            \
  DEFINE_VECTORMATH_ANY( XLALVectorMath_SS2S_SSEx, NAME ## REAL4, ( REAL4 *out, const REAL4 *in1, const REAL4 *in2, const UINT4 len ), ( (out != NULL) && (in1 != NULL) && (in2 != NULL) ), ( out, in1, in2, len, SSE_OP ) )

DEFINE_VECTORMATH_SS2S(Add, local_add_ps)
DEFINE_VECTORMATH_SS2S(Multiply, local_mul_ps)
DEFINE_VECTORMATH_SS2S(Max, local_max_ps)

// ---------- define vector math functions with 1 REAL4 scalar and 1 REAL4 vector inputs to 1 REAL4 vector output (sS2S) ----------
#define DEFINE_VECTORMATH_sS2S(NAME, SSE_OP)                            \
  DEFINE_VECTORMATH_ANY( XLALVectorMath_sS2S_SSEx, NAME ## REAL4, ( REAL4 *out, REAL4 scalar, const REAL4 *in, const UINT4 len ), ( (out != NULL) && (in != NULL) ), ( out, scalar, in, len, SSE_OP ) )

DEFINE_VECTORMATH_sS2S(Shift, local_add_ps)
DEFINE_VECTORMATH_sS2S(Scale, local_mul_ps)

// ---------- define vector math functions with 1 REAL8 scalar and 1 REAL8 vector inputs to 1 REAL8 vector output (sS2S) ----------
#define DEFINE_VECTORMATH_dD2D(NAME, SSE_OP)                            \
  DEFINE_VECTORMATH_ANY( XLALVectorMath_dD2D_SSEx, NAME ## REAL8, ( REAL8 *out, REAL8 scalar, const REAL8 *in, const UINT4 len ), ( (out != NULL) && (in != NULL) ), ( out, scalar, in, len, SSE_OP ) )

DEFINE_VECTORMATH_dD2D(Scale, local_mul_pd)
DEFINE_VECTORMATH_dD2D(Shift, local_add_pd)

// ---------- define vector math functions with 2 REAL8 vector inputs to 1 REAL8 vector output (DD2D) ----------
#define DEFINE_VECTORMATH_DD2D(NAME, SSE_OP)                            \
  DEFINE_VECTORMATH_ANY( XLALVectorMath_DD2D_SSEx, NAME ## REAL8, ( REAL8 *out, const REAL8 *in1, const REAL8 *in2, const UINT4 len ), ( (out != NULL) && (in1 != NULL) && (in2 != NULL) ), ( out, in1, in2, len, SSE_OP ) )

DEFINE_VECTORMATH_DD2D(Add, local_add_pd)
DEFINE_VECTORMATH_DD2D(Multiply, local_mul_pd)

// ---------- define vector math functions with 2 COMPLEX8 vector inputs to 1 COMPLEX8 vector output (CC2C) ----------
#define DEFINE_VECTORMATH_CC2C(NAME, SSE_OP)                            \
  DEFINE_VECTORMATH_ANY( XLALVectorMath_CC2C_SSEx, NAME ## COMPLEX8, ( COMPLEX8 *out, const COMPLEX8 *in1, const COMPLEX8 *in2, const UINT4 len ), ( (out != NULL) && (in1 != NULL) && (in2 != NULL) ), ( out, in1, in2, len, SSE_OP ) )

DEFINE_VECTORMATH_CC2C(Multiply, local_cmul_ps)
