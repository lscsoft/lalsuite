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

#ifdef __SSE2__
#define USE_SSE2
#endif

#include "VectorMath_sse_mathfun.h"

// ========== internal generic SSEx functions ==========

static inline int
XLALVectorMath_S2S_SSEx ( REAL4 *out, const REAL4 *in, const UINT4 len, __m128 (*f)(__m128) )
{
  XLAL_CHECK ( isMemAligned(out, 16) && isMemAligned(in, 16), XLAL_EINVAL, "All vectors need to be 16-byte aligned (4xfloats\n");

  // walk through vector in blocks of 4
  UINT4 i4Max = len - ( len % 4 );
  const __m128 *in4p = (const void *)&(in[0]);
  __m128 *out4p = (void *)&out[0];

  for ( UINT4 i4 = 0; i4 < i4Max; i4 += 4 )
    {
      (*out4p++) = (*f)( (*in4p++) );
    } // for i4 in 0 ... i4Max/4-1

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

static inline int
XLALVectorMath_S2SS_SSEx ( REAL4 *out1, REAL4 *out2, const REAL4 *in, const UINT4 len, void (*f)(__m128, __m128*, __m128*) )
{
  XLAL_CHECK ( isMemAligned(out1, 16) && isMemAligned(out2, 16) && isMemAligned(in, 16), XLAL_EINVAL, "All vectors need to be 16-byte aligned (4xfloats\n");

  // walk through vector in blocks of 4
  UINT4 i4Max = len - ( len % 4 );
  const __m128 *in4p = (const void *)&(in[0]);
  __m128 *out4p_1 = (void *)&(out1[0]);
  __m128 *out4p_2 = (void *)&(out2[0]);

  for ( UINT4 i4 = 0; i4 < i4Max; i4 += 4 )
    {
      (*f) ( (*in4p++), out4p_1++, out4p_2++ );
    } // for i4 in 0 ... i4Max/4-1

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
