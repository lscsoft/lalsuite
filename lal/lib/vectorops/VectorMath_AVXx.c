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

#ifndef __AVX__
#error "VectorMath_AVXx.c requires SIMD instruction set AVX"
#endif

#include "VectorMath_avx_mathfun.h"

// ---------- local operators and operator-wrappers ----------
UNUSED static inline __m256
local_add_ps ( __m256 in1, __m256 in2 )
{
  return _mm256_add_ps ( in1, in2 );
}

UNUSED static inline __m256
local_sub_ps ( __m256 in1, __m256 in2 )
{
  return _mm256_sub_ps ( in1, in2 );
}

UNUSED static inline __m256
local_mul_ps ( __m256 in1, __m256 in2 )
{
  return _mm256_mul_ps ( in1, in2 );
}

UNUSED static inline __m256
local_max_ps ( __m256 in1, __m256 in2 )
{
  return _mm256_max_ps ( in1, in2 );
}

UNUSED static inline __m256d
local_add_pd ( __m256d in1, __m256d in2 )
{
  return _mm256_add_pd ( in1, in2 );
}

UNUSED static inline __m256d
local_sub_pd ( __m256d in1, __m256d in2 )
{
  return _mm256_sub_pd ( in1, in2 );
}

UNUSED static inline __m256d
local_mul_pd ( __m256d in1, __m256d in2 )
{
  return _mm256_mul_pd ( in1, in2 );
}

UNUSED static inline __m256d
local_max_pd ( __m256d in1, __m256d in2 )
{
  return _mm256_max_pd ( in1, in2 );
}

UNUSED static inline __m256
local_round_ps ( __m256 in )
{

  __m256 result = _mm256_round_ps ( in, 0x0 );

  //calculate absolute values from result and in
  __m256 neg = _mm256_set1_ps ( -1.0f );
  __m256 in_neg = _mm256_mul_ps ( in, neg );
  __m256 result_neg = _mm256_mul_ps ( result ,neg );
  __m256 in_abs = _mm256_max_ps ( in, in_neg );
  __m256 result_abs = _mm256_max_ps ( result, result_neg );

  //test if we had round half down to even
  __m256 diff = _mm256_sub_ps ( result_abs, in_abs );
  __m256 test = _mm256_set1_ps ( -0.5f );
  __m256 cmp = _mm256_cmp_ps ( diff, test, _CMP_EQ_OQ );
  int mask = _mm256_movemask_ps ( cmp );

  if ( mask != 0 )
    {
      //add or sub 1.0 to the 'rounded half down to even' values
      __m256 pos = _mm256_set1_ps ( 1.0f );
      __m256 zero = _mm256_setzero_ps ();
      __m256 even_down_rounds = _mm256_and_ps ( in, cmp );
      __m256 cmp_neg = _mm256_cmp_ps ( even_down_rounds, zero, _CMP_LT_OQ );
      __m256 cmp_pos = _mm256_cmp_ps ( even_down_rounds, zero, _CMP_GT_OQ );
      __m256 add = _mm256_and_ps ( pos, cmp_pos );
      __m256 sub = _mm256_and_ps ( neg, cmp_neg );
      __m256 tot = _mm256_add_ps ( add, sub );
      result = _mm256_add_ps ( result, tot );
    }

  return result;
}

UNUSED static inline __m256d
local_round_pd ( __m256d in )
{

  __m256d result = _mm256_round_pd ( in, 0x0 );

  //calculate absolute values from result and in
  __m256d neg = _mm256_set1_pd ( -1.0 );
  __m256d in_neg = _mm256_mul_pd ( in, neg );
  __m256d result_neg = _mm256_mul_pd ( result, neg );
  __m256d in_abs = _mm256_max_pd ( in, in_neg );
  __m256d result_abs = _mm256_max_pd ( result, result_neg );

  //test if we had round half down to even
  __m256d diff = _mm256_sub_pd ( result_abs, in_abs );
  __m256d test = _mm256_set1_pd ( -0.5 );
  __m256d cmp = _mm256_cmp_pd ( diff, test, _CMP_EQ_OQ );
  int mask = _mm256_movemask_pd ( cmp );

  if ( mask != 0 )
    {
      //add or sub 1.0 to the 'rounded half down to even' values
      __m256d pos = _mm256_set1_pd ( 1.0 );
      __m256d zero = _mm256_setzero_pd ();
      __m256d even_down_rounds = _mm256_and_pd ( in, cmp );
      __m256d cmp_neg = _mm256_cmp_pd ( even_down_rounds, zero, _CMP_LT_OQ );
      __m256d cmp_pos = _mm256_cmp_pd ( even_down_rounds, zero, _CMP_GT_OQ );
      __m256d add = _mm256_and_pd ( pos, cmp_pos );
      __m256d sub = _mm256_and_pd ( neg, cmp_neg );
      __m256d tot = _mm256_add_pd ( add, sub );
      result = _mm256_add_pd ( result, tot );
    }

  return result;
}

// in1: a0,b0,a1,b1,a2,b2,a3,b3 in2: c0,d0,c1,d1,c2,d2,c3,d3
UNUSED static inline __m256
local_cmul_ps ( __m256 in1, __m256 in2 )
{
  __m256 neg = _mm256_setr_ps(1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0);

  // Multiply in1 and in2
  // a0c0, b0d0, a1c1, b2d2, a3c3, b3d3
  __m256 temp1 = _mm256_mul_ps(in1, in2);

  // Switch the real and imaginary elements of in2
  in2 = _mm256_permute_ps(in2, 0xb1);

  //Negate the real elements of in2
  in2 = _mm256_mul_ps(in2, neg);

  // Multiply in1 and the modified in2
  // a0d0, -b0c0, a1d1, -b1c1, a2d2, -b2c2, a3d3, -b3c3
  __m256 temp2 = _mm256_mul_ps(in1, in2);

  // Horizontally subtract the elements in temp1 and temp2
  // a0c0-b0d0, a1c1-b1d1, a0d0+b0c0, a1d1+b1c1,a2c2-b2d2, a3c3-b3d3, a2d2+b2c2, a3d3+b3c3
  in2 = _mm256_hsub_ps(temp1, temp2);

  // swap odd numbers real component with even numbers imaginary component
  return _mm256_permute_ps(in2, 0xd8);
}

// ========== internal generic AVXx functions ==========

// ---------- generic AVXx operator with 1 REAL4 vector input to 1 REAL4 vector output (S2S) ----------
static inline int
XLALVectorMath_S2S_AVXx ( REAL4 *out, const REAL4 *in, const UINT4 len, __m256 (*f)(__m256) )
{

  // walk through vector in blocks of 8
  UINT4 i8Max = len - ( len % 8 );
  for ( UINT4 i8 = 0; i8 < i8Max; i8 += 8 )
    {
      __m256 in8p = _mm256_loadu_ps(&in[i8]);
      __m256 out8p = (*f)( in8p );
      _mm256_storeu_ps(&out[i8], out8p);
    }

  // deal with the remaining (<=7) terms separately
  V8SF in8 = {.f={0,0,0,0,0,0,0,0}}, out8;
  for ( UINT4 i = i8Max,j=0; i < len; i ++, j++ ) {
    in8.f[j] = in[i];
  }
  out8.v = (*f)( in8.v );
  for ( UINT4 i = i8Max,j=0; i < len; i ++, j++ ) {
    out[i] = out8.f[j];
  }

  return XLAL_SUCCESS;

} // XLALVectorMath_S2S_AVXx()

// ---------- generic AVXx operator with 1 REAL4 vector input to 2 REAL4 vector outputs (S2SS) ----------
static inline int
XLALVectorMath_S2SS_AVXx ( REAL4 *out1, REAL4 *out2, const REAL4 *in, const UINT4 len, void (*f)(__m256, __m256*, __m256*) )
{

  // walk through vector in blocks of 8
  UINT4 i8Max = len - ( len % 8 );
  for ( UINT4 i8 = 0; i8 < i8Max; i8 += 8)
    {
      __m256 in8p = _mm256_loadu_ps(&in[i8]);
      __m256 out8p_1, out8p_2;
      (*f) ( in8p, &out8p_1, &out8p_2 );
      _mm256_storeu_ps(&out1[i8], out8p_1);
      _mm256_storeu_ps(&out2[i8], out8p_2);
    }

  // deal with the remaining (<=7) terms separately
  V8SF in8 = {.f={0,0,0,0,0,0,0,0}}, out8_1, out8_2;
  for ( UINT4 i = i8Max,j=0; i < len; i ++, j++ ) {
    in8.f[j] = in[i];
  }
  (*f) ( in8.v, &out8_1.v, &out8_2.v );
  for ( UINT4 i = i8Max,j=0; i < len; i ++, j++ ) {
    out1[i]  = out8_1.f[j];
    out2[i] = out8_2.f[j];
  }

  return XLAL_SUCCESS;

} // XLALVectorMath_S2SS_AVXx()

// ---------- generic AVXx operator with 2 REAL4 vector inputs to 1 REAL4 vector output (SS2S) ----------
static inline int
XLALVectorMath_SS2S_AVXx ( REAL4 *out, const REAL4 *in1, const REAL4 *in2, const UINT4 len, __m256 (*op)(__m256, __m256) )
{

  // walk through vector in blocks of 8
  UINT4 i8Max = len - ( len % 8 );
  for ( UINT4 i8 = 0; i8 < i8Max; i8 += 8 )
    {
      __m256 in8p_1 = _mm256_loadu_ps(&in1[i8]);
      __m256 in8p_2 = _mm256_loadu_ps(&in2[i8]);
      __m256 out8p = (*op) ( in8p_1, in8p_2 );
      _mm256_storeu_ps(&out[i8], out8p);
    }

  // deal with the remaining (<=7) terms separately
  V8SF in8_1 = {.f={0,0,0,0,0,0,0,0}};
  V8SF in8_2 = {.f={0,0,0,0,0,0,0,0}};
  V8SF out8;
  for ( UINT4 i = i8Max,j=0; i < len; i ++, j++ ) {
    in8_1.f[j] = in1[i];
    in8_2.f[j] = in2[i];
  }
  out8.v = (*op) ( in8_1.v, in8_2.v );
  for ( UINT4 i = i8Max,j=0; i < len; i ++, j++ ) {
    out[i] = out8.f[j];
  }

  return XLAL_SUCCESS;

} // XLALVectorMath_SS2S_AVXx()

// ---------- generic SSEx operator with 1 REAL4 scalar and 1 REAL4 vector inputs to 1 REAL4 vector output (sS2S) ----------
static inline int
XLALVectorMath_sS2S_AVXx ( REAL4 *out, REAL4 scalar, const REAL4 *in, const UINT4 len, __m256 (*op)(__m256, __m256) )
{
  const V8SF scalar8 = {.f={scalar,scalar,scalar,scalar,scalar,scalar,scalar,scalar}};

  // walk through vector in blocks of 8
  UINT4 i8Max = len - ( len % 8 );
  for ( UINT4 i8 = 0; i8 < i8Max; i8 += 8 )
    {
      __m256 in8p = _mm256_loadu_ps(&in[i8]);
      __m256 out8p = (*op) ( scalar8.v, in8p );
      _mm256_storeu_ps(&out[i8], out8p);
    }

  // deal with the remaining (<=7) terms separately
  V8SF in8 = {.f={0,0,0,0,0,0,0,0}};
  V8SF out8;
  for ( UINT4 i = i8Max,j=0; i < len; i ++, j++ ) {
    in8.f[j] = in[i];
  }
  out8.v = (*op) ( scalar8.v, in8.v );
  for ( UINT4 i = i8Max,j=0; i < len; i ++, j++ ) {
    out[i] = out8.f[j];
  }

  return XLAL_SUCCESS;

} // XLALVectorMath_sS2S_AVXx()

// ---------- generic AVXx operator with 1 REAL8 scalar and 1 REAL8 vector inputs to 1 REAL8 vector output (dD2D) ----------
static inline int
XLALVectorMath_dD2D_AVXx ( REAL8 *out, REAL8 scalar, const REAL8 *in, const UINT4 len, __m256d (*op)(__m256d, __m256d) )
{
  const V4SD scalar4 = {.f={scalar,scalar,scalar,scalar}};

  // walk through vector in blocks of 4
  UINT4 i4Max = len - ( len % 4 );
  for ( UINT4 i4 = 0; i4 < i4Max; i4 += 4 )
    {
      __m256d in4p = _mm256_loadu_pd(&in[i4]);
      __m256d out4p = (*op) ( scalar4.v, in4p );
      _mm256_storeu_pd(&out[i4], out4p);
    }

  // deal with the remaining (<=3) terms separately
  V4SD in4 = {.f={0,0,0,0}};
  V4SD out4;
  for ( UINT4 i = i4Max,j=0; i < len; i ++, j++ ) {
    in4.f[j] = in[i];
  }
  out4.v = (*op) ( scalar4.v, in4.v );
  for ( UINT4 i = i4Max,j=0; i < len; i ++, j++ ) {
    out[i] = out4.f[j];
  }

  return XLAL_SUCCESS;

} // XLALVectorMath_dD2D_AVXx()

// ---------- generic AVXx operator with 2 REAL8 vector inputs to 1 REAL8 vector output (DD2D) ----------
static inline int
XLALVectorMath_DD2D_AVXx ( REAL8 *out, const REAL8 *in1, const REAL8 *in2, const UINT4 len, __m256d (*op)(__m256d, __m256d) )
{

  // walk through vector in blocks of 4
  UINT4 i4Max = len - ( len % 4 );
  for ( UINT4 i4 = 0; i4 < i4Max; i4 += 4 )
    {
      __m256d in4p_1 = _mm256_loadu_pd(&in1[i4]);
      __m256d in4p_2 = _mm256_loadu_pd(&in2[i4]);
      __m256d out4p = (*op) ( in4p_1, in4p_2 );
      _mm256_storeu_pd(&out[i4], out4p);
    }

  // deal with the remaining (<=3) terms separately
  V4SD in4_1 = {.f={0,0,0,0}};
  V4SD in4_2 = {.f={0,0,0,0}};
  V4SD out4;
  for ( UINT4 i = i4Max,j=0; i < len; i ++, j++ ) {
    in4_1.f[j] = in1[i];
    in4_2.f[j] = in2[i];
  }
  out4.v = (*op) ( in4_1.v, in4_2.v );
  for ( UINT4 i = i4Max,j=0; i < len; i ++, j++ ) {
    out[i] = out4.f[j];
  }

  return XLAL_SUCCESS;

} // XLALVectorMath_DD2D_AVXx()

// ---------- generic AVXx operator with 2 COMPLEX8 vector inputs to 1 COMPLEX8 vector output (CC2C) ----------
static inline int
XLALVectorMath_CC2C_AVXx ( COMPLEX8 *out, const COMPLEX8 *in1, const COMPLEX8 *in2, const UINT4 len, __m256 (*op)(__m256, __m256) )
{

  // walk through vector in blocks of 4
  UINT4 i4Max = len - ( len % 4 );
  for ( UINT4 i4 = 0; i4 < i4Max; i4+=4 )
    {
      __m256 in8p_1 = _mm256_loadu_ps( (const REAL4*)&in1[i4] );
      __m256 in8p_2 = _mm256_loadu_ps( (const REAL4*)&in2[i4] );
      __m256 out8p = (*op) ( in8p_1, in8p_2 );
      _mm256_storeu_ps( (REAL4*)&out[i4], out8p );
    }

  // deal with the remaining (<=3) terms separately
  V8SF in8_1 = {.f={0,0,0,0,0,0,0,0}};
  V8SF in8_2 = {.f={0,0,0,0,0,0,0,0}};
  V8SF out8;
  for ( UINT4 i = i4Max,j=0; i < len ; i++,j+=2)
    {
      in8_1.f[j]   = crealf ( in1[i] );
      in8_1.f[j+1] = cimagf ( in1[i] );
      in8_2.f[j]   = crealf ( in2[i] );
      in8_2.f[j+1] = cimagf ( in2[i] );
    }

  out8.v = (*op) ( in8_1.v, in8_2.v );
  for ( UINT4 i = i4Max, j = 0; i < len; i++,j+=2 )
    {
      out[i] = crect( out8.f[j], out8.f[j+1] );
    }

  return XLAL_SUCCESS;

} // XLALVectorMath_CC2C_AVXx()

// ---------- generic AVXx operator with 1 COMPLEX8 scalar and 1 COMPLEX8 vector inputs to 1 COMPLEX8 vector output (cC2C) ----------
static inline int
XLALVectorMath_cC2C_AVXx ( COMPLEX8 *out, COMPLEX8 scalar, const COMPLEX8 *in, const UINT4 len, __m256 (*op)(__m256, __m256) )
{
  const V8SF scalar8 = {.f={crealf(scalar),cimagf(scalar),crealf(scalar),cimagf(scalar),crealf(scalar),cimagf(scalar),crealf(scalar),cimagf(scalar)}};

  // walk through vector in blocks of 4
  UINT4 i4Max = len - ( len % 4 );
  for ( UINT4 i4 = 0; i4 < i4Max; i4+=4 )
    {
      __m256 in8p = _mm256_loadu_ps( (const REAL4*)&in[i4] );
      __m256 out8p = (*op) ( scalar8.v, in8p );
      _mm256_storeu_ps( (REAL4*)&out[i4], out8p );
    }

  // deal with the remaining (<=3) terms separately
  V8SF in8 = {.f={0,0,0,0,0,0,0,0}};
  V8SF out8;
  for ( UINT4 i = i4Max,j=0; i < len ; i++, j+=2 )
    {
      in8.f[j]   = crealf ( in[i] );
      in8.f[j+1] = cimagf ( in[i] );
     }

  out8.v = (*op) ( scalar8.v, in8.v );
  for ( UINT4 i = i4Max,j=0; i < len; i++,j+=2 )
    {
      out[i] = crect( out8.f[j], out8.f[j+1] );
    }

  return XLAL_SUCCESS;

} // XLALVectorMath_cC2C_AVXx()

// ---------- generic AVXx operator with 1 REAL8 vector input to 1 REAL8 vector output (D2D) ----------
static inline int
XLALVectorMath_D2D_AVXx ( REAL8 *out, const REAL8 *in, const UINT4 len, __m256d (*f)(__m256d) )
{

  // walk through vector in blocks of 4
  UINT4 i4Max = len - ( len % 4 );
  for ( UINT4 i4 = 0; i4 < i4Max; i4 += 4 )
    {
      __m256d in4p = _mm256_loadu_pd(&in[i4]);
      __m256d out4p = (*f)( in4p );
      _mm256_storeu_pd(&out[i4], out4p);
    }

  // deal with the remaining (<=3) terms separately
  V4SD in4 = {.f={0,0,0,0}}, out4;
  for ( UINT4 i = i4Max,j=0; i < len; i ++, j++ ) {
    in4.f[j] = in[i];
  }
  out4.v = (*f)( in4.v );
  for ( UINT4 i = i4Max,j=0; i < len; i ++, j++ ) {
    out[i] = out4.f[j];
  }

  return XLAL_SUCCESS;

} // XLALVectorMath_D2D_AVXx()

// ========== internal AVXx vector math functions ==========

// ---------- define vector math functions with 1 REAL4 vector input to 1 REAL4 vector output (S2S) ----------
#define DEFINE_VECTORMATH_S2S(NAME, AVX_OP)                             \
  DEFINE_VECTORMATH_ANY( XLALVectorMath_S2S_AVXx, NAME ## REAL4, ( REAL4 *out, const REAL4 *in, const UINT4 len ), ( (out != NULL) && (in != NULL) ), ( out, in, len, AVX_OP ) )

DEFINE_VECTORMATH_S2S(Sin, sin256_ps)
DEFINE_VECTORMATH_S2S(Cos, cos256_ps)
DEFINE_VECTORMATH_S2S(Exp, exp256_ps)
DEFINE_VECTORMATH_S2S(Log, log256_ps)
DEFINE_VECTORMATH_S2S(Round, local_round_ps)

// ---------- define vector math functions with 1 REAL4 vector input to 2 REAL4 vector outputs (S2SS) ----------
#define DEFINE_VECTORMATH_S2SS(NAME, AVX_OP)                            \
  DEFINE_VECTORMATH_ANY( XLALVectorMath_S2SS_AVXx, NAME ## REAL4, ( REAL4 *out1, REAL4 *out2, const REAL4 *in, const UINT4 len ), ( (out1 != NULL) && (out2 != NULL) && (in != NULL) ), ( out1, out2, in, len, AVX_OP ) )

DEFINE_VECTORMATH_S2SS(SinCos, sincos256_ps)
DEFINE_VECTORMATH_S2SS(SinCos2Pi, sincos256_ps_2pi)

// ---------- define vector math functions with 2 REAL4 vector inputs to 1 REAL4 vector output (SS2S) ----------
#define DEFINE_VECTORMATH_SS2S(NAME, AVX_OP)                            \
  DEFINE_VECTORMATH_ANY( XLALVectorMath_SS2S_AVXx, NAME ## REAL4, ( REAL4 *out, const REAL4 *in1, const REAL4 *in2, const UINT4 len ), ( (out != NULL) && (in1 != NULL) && (in2 != NULL) ), ( out, in1, in2, len, AVX_OP ) )

DEFINE_VECTORMATH_SS2S(Add, local_add_ps)
DEFINE_VECTORMATH_SS2S(Sub, local_sub_ps)
DEFINE_VECTORMATH_SS2S(Multiply, local_mul_ps)
DEFINE_VECTORMATH_SS2S(Max, local_max_ps)

// ---------- define vector math functions with 1 REAL4 scalar and 1 REAL4 vector inputs to 1 REAL4 vector output (sS2S) ----------
#define DEFINE_VECTORMATH_sS2S(NAME, AVX_OP)                            \
  DEFINE_VECTORMATH_ANY( XLALVectorMath_sS2S_AVXx, NAME ## REAL4, ( REAL4 *out, REAL4 scalar, const REAL4 *in, const UINT4 len ), ( (out != NULL) && (in != NULL) ), ( out, scalar, in, len, AVX_OP ) )

DEFINE_VECTORMATH_sS2S(Shift, local_add_ps)
DEFINE_VECTORMATH_sS2S(Scale, local_mul_ps)

// ---------- define vector math functions with 1 REAL8 scalar and 1 REAL8 vector inputs to 1 REAL8 vector output (dD2D) ----------
#define DEFINE_VECTORMATH_dD2D(NAME, GEN_OP)                            \
  DEFINE_VECTORMATH_ANY( XLALVectorMath_dD2D_AVXx, NAME ## REAL8, ( REAL8 *out, REAL8 scalar, const REAL8 *in, const UINT4 len ), ( (out != NULL) && (in != NULL) ), ( out, scalar, in, len, GEN_OP ) )

DEFINE_VECTORMATH_dD2D(Scale, local_mul_pd)
DEFINE_VECTORMATH_dD2D(Shift, local_add_pd)

// ---------- define vector math functions with 2 REAL8 vector inputs to 1 REAL8 vector output (DD2D) ----------
#define DEFINE_VECTORMATH_DD2D(NAME, AVX_OP)                            \
  DEFINE_VECTORMATH_ANY( XLALVectorMath_DD2D_AVXx, NAME ## REAL8, ( REAL8 *out, const REAL8 *in1, const REAL8 *in2, const UINT4 len ), ( (out != NULL) && (in1 != NULL) && (in2 != NULL) ), ( out, in1, in2, len, AVX_OP ) )

DEFINE_VECTORMATH_DD2D(Add, local_add_pd)
DEFINE_VECTORMATH_DD2D(Sub, local_sub_pd)
DEFINE_VECTORMATH_DD2D(Multiply, local_mul_pd)
DEFINE_VECTORMATH_DD2D(Max, local_max_pd)

// ---------- define vector math functions with 2 COMPLEX8 vector inputs to 1 COMPLEX8 vector output (CC2C) ----------
#define DEFINE_VECTORMATH_CC2C(NAME, AVX_OP)                            \
  DEFINE_VECTORMATH_ANY( XLALVectorMath_CC2C_AVXx, NAME ## COMPLEX8, ( COMPLEX8 *out, const COMPLEX8 *in1, const COMPLEX8 *in2, const UINT4 len ), ( (out != NULL) && (in1 != NULL) && (in2 != NULL) ), ( out, in1, in2, len, AVX_OP ) )

DEFINE_VECTORMATH_CC2C(Multiply, local_cmul_ps)
DEFINE_VECTORMATH_CC2C(Add, local_add_ps)

// ---------- define vector math functions with 1 COMPLEX8 scalar and 1 COMPLEX8 vector inputs to 1 COMPLEX8 vector output (cC2C) ----------
#define DEFINE_VECTORMATH_cC2C(NAME, AVX_OP)                            \
  DEFINE_VECTORMATH_ANY( XLALVectorMath_cC2C_AVXx, NAME ## COMPLEX8, ( COMPLEX8 *out, COMPLEX8 scalar, const COMPLEX8 *in, const UINT4 len ), ( (out != NULL) && (in != NULL) ), ( out, scalar, in, len, AVX_OP ) )

DEFINE_VECTORMATH_cC2C(Scale, local_cmul_ps)
DEFINE_VECTORMATH_cC2C(Shift, local_add_ps)

// ---------- define vector math functions with 1 REAL8 vector input to 1 REAL8 vector output (D2D) ----------
#define DEFINE_VECTORMATH_D2D(NAME, AVX_OP)                             \
  DEFINE_VECTORMATH_ANY( XLALVectorMath_D2D_AVXx, NAME ## REAL8, ( REAL8 *out, const REAL8 *in, const UINT4 len ), ( (out != NULL) && (in != NULL) ), ( out, in, len, AVX_OP ) )

DEFINE_VECTORMATH_D2D(Round, local_round_pd)
