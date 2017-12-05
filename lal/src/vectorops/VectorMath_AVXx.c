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

#include "VectorMath_avx_mathfun.h"

// ---------- local operators and operator-wrappers ----------
UNUSED static inline __m256
local_add_ps ( __m256 in1, __m256 in2 )
{
  return _mm256_add_ps ( in1, in2 );
}

UNUSED static inline __m256
local_mul_ps ( __m256 in1, __m256 in2 )
{
  return _mm256_mul_ps ( in1, in2 );
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

// ========== internal AVXx vector math functions ==========

// ---------- define vector math functions with 1 REAL4 vector input to 1 REAL4 vector output (S2S) ----------
#define DEFINE_VECTORMATH_S2S(NAME, AVX_OP)                             \
  DEFINE_VECTORMATH_ANY( XLALVectorMath_S2S_AVXx, NAME ## REAL4, ( REAL4 *out, const REAL4 *in, const UINT4 len ), ( (out != NULL) && (in != NULL) ), ( out, in, len, AVX_OP ) )

DEFINE_VECTORMATH_S2S(Sin, sin256_ps)
DEFINE_VECTORMATH_S2S(Cos, cos256_ps)
DEFINE_VECTORMATH_S2S(Exp, exp256_ps)
DEFINE_VECTORMATH_S2S(Log, log256_ps)

// ---------- define vector math functions with 1 REAL4 vector input to 2 REAL4 vector outputs (S2SS) ----------
#define DEFINE_VECTORMATH_S2SS(NAME, AVX_OP)                            \
  DEFINE_VECTORMATH_ANY( XLALVectorMath_S2SS_AVXx, NAME ## REAL4, ( REAL4 *out1, REAL4 *out2, const REAL4 *in, const UINT4 len ), ( (out1 != NULL) && (out2 != NULL) && (in != NULL) ), ( out1, out2, in, len, AVX_OP ) )

DEFINE_VECTORMATH_S2SS(SinCos, sincos256_ps)
DEFINE_VECTORMATH_S2SS(SinCos2Pi, sincos256_ps_2pi)

// ---------- define vector math functions with 2 REAL4 vector inputs to 1 REAL4 vector output (SS2S) ----------
#define DEFINE_VECTORMATH_SS2S(NAME, AVX_OP)                            \
  DEFINE_VECTORMATH_ANY( XLALVectorMath_SS2S_AVXx, NAME ## REAL4, ( REAL4 *out, const REAL4 *in1, const REAL4 *in2, const UINT4 len ), ( (out != NULL) && (in1 != NULL) && (in2 != NULL) ), ( out, in1, in2, len, AVX_OP ) )

DEFINE_VECTORMATH_SS2S(Add, local_add_ps)
DEFINE_VECTORMATH_SS2S(Multiply, local_mul_ps)

// ---------- define vector math functions with 1 REAL4 scalar and 1 REAL4 vector inputs to 1 REAL4 vector output (sS2S) ----------
#define DEFINE_VECTORMATH_sS2S(NAME, AVX_OP)                            \
  DEFINE_VECTORMATH_ANY( XLALVectorMath_sS2S_AVXx, NAME ## REAL4, ( REAL4 *out, REAL4 scalar, const REAL4 *in, const UINT4 len ), ( (out != NULL) && (in != NULL) ), ( out, scalar, in, len, AVX_OP ) )

DEFINE_VECTORMATH_sS2S(Shift, local_add_ps)
DEFINE_VECTORMATH_sS2S(Scale, local_mul_ps)
