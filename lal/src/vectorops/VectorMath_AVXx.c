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

/* prevent inclusion of conflicting non-C99 symbols */
#define _ISOC99_SOURCE

// ---------- INCLUDES ----------
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <config.h>

#include <lal/LALConstants.h>
#include <lal/VectorMath.h>

#include "VectorMath_internal.h"

#include "VectorMath_avx_mathfun.h"

// ========== internal generic AVXx functions ==========

static inline int
XLALVectorMath_S2S_AVXx ( REAL4 *out, const REAL4 *in, const UINT4 len, __m256 (*f)(__m256) )
{
  XLAL_CHECK ( isMemAligned(out, 32) && isMemAligned(in, 32), XLAL_EINVAL, "All vectors need to be 32-byte aligned (8xfloats\n");

  // walk through vector in blocks of 8
  UINT4 i8Max = len - ( len % 8 );
  const __m256 *in8p = (const void *)&(in[0]);
  __m256 *out8p = (void *)&(out[0]);

  for ( UINT4 i8 = 0; i8 < i8Max; i8 += 8 )
    {
      (*out8p++) = (*f)( (*in8p++) );
    } // for i8 in 0 ... i8Max/8-1

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

static inline int
XLALVectorMath_S2SS_AVXx ( REAL4 *out1, REAL4 *out2, const REAL4 *in, const UINT4 len, void (*f)(__m256, __m256*, __m256*) )
{
  XLAL_CHECK ( isMemAligned(out1, 32) && isMemAligned(out2, 32) && isMemAligned(in, 32), XLAL_EINVAL, "All vectors need to be 32-byte aligned (8xfloats\n");

  // walk through vector in blocks of 8
  UINT4 i8Max = len - ( len % 8 );
  const __m256 *in8p = (const void *)&(in[0]);
  __m256 *out8p_1 = (void *)&(out1[0]);
  __m256 *out8p_2 = (void *)&(out2[0]);

  for ( UINT4 i8 = 0; i8 < i8Max; i8 += 8)
    {
      (*f) ( (*in8p++), out8p_1++, out8p_2++ );
    } // for i8 in 0 ... i8Max/8-1

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
