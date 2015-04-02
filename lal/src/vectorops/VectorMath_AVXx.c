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

// ========== internal generic SSEx functions ==========

static inline int
XLALVectorFuncf1T1_AVXx ( REAL4 *out, const REAL4 *in, const UINT4 len, v8sf (*f)(v8sf) )
{
  XLAL_CHECK ( isMemAligned(out, 32) && isMemAligned(in, 32), XLAL_EINVAL, "All vectors need to be 32-byte aligned (8xfloats\n");

  // walk through vector in blocks of 8
  UINT4 i8Max = len - ( len % 8 );
  for ( UINT4 i8 = 0; i8 < i8Max/8; i8 ++ )
    {
      UINT4 i0 = i8 * 8;

      const v8sf *in8p = (const void *)&(in[i0]);
      v8sf *out8p = (void *)&(out[i0]);

      (*out8p) = (*f)( (*in8p) );

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

} // XLALVectorFuncf1T1_AVXx()

static inline int
XLALVectorFuncf1T2_AVXx ( REAL4 *out1, REAL4 *out2, const REAL4 *in, const UINT4 len, void (*f)(v8sf, v8sf*, v8sf*) )
{
  XLAL_CHECK ( isMemAligned(out1, 32) && isMemAligned(out2, 32) && isMemAligned(in, 32), XLAL_EINVAL, "All vectors need to be 32-byte aligned (8xfloats\n");

  // walk through vector in blocks of 4
  UINT4 i8Max = len - ( len % 8 );
  for ( UINT4 i8 = 0; i8 < i8Max/8; i8 ++ )
    {
      UINT4 i0 = i8 * 8;

      const v8sf *in8p = (const void *)&(in[i0]);
      v8sf *out8p_1 = (void *)&(out1[i0]);
      v8sf *out8p_2 = (void *)&(out2[i0]);

      (*f) ( (*in8p), out8p_1, out8p_2 );

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

} // XLALVectorFuncf1T2_AVXx()

// ========== internal AVXx vector math functions ==========

#define DEFINE_VECTORMATH_FUNCF_1T1(NAME, SSE_FUNC) \
  DEFINE_VECTORMATH_ANY( XLALVectorFuncf1T1_AVXx, NAME, ( REAL4 *out, const REAL4 *in, const UINT4 len ), ( (out != NULL) && (in != NULL) ), ( out, in, len, SSE_FUNC ) )

DEFINE_VECTORMATH_FUNCF_1T1(Sinf, sin256_ps)
DEFINE_VECTORMATH_FUNCF_1T1(Cosf, cos256_ps)
DEFINE_VECTORMATH_FUNCF_1T1(Expf, exp256_ps)
DEFINE_VECTORMATH_FUNCF_1T1(Logf, log256_ps)

#define DEFINE_VECTORMATH_FUNCF_1T2(NAME, SSE_FUNC) \
  DEFINE_VECTORMATH_ANY( XLALVectorFuncf1T2_AVXx, NAME, ( REAL4 *out1, REAL4 *out2, const REAL4 *in, const UINT4 len ), ( (out1 != NULL) && (out2 != NULL) && (in != NULL) ), ( out1, out2, in, len, SSE_FUNC ) )

DEFINE_VECTORMATH_FUNCF_1T2(SinCosf, sincos256_ps)
DEFINE_VECTORMATH_FUNCF_1T2(SinCosf2PI, sincos256_ps_2pi)
