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

#include "VectorMath_internal.h"

#define SIMD_INSTRSET GEN

// ---------- local math functions ----------

static inline void local_sincosf(REAL4 in, REAL4 *out1, REAL4 *out2) {
  *out1 = sinf ( in );
  *out2 = cosf ( in );
}

static inline void local_sincosf_2pi(REAL4 in, REAL4 *out1, REAL4 *out2) {
  *out1 = sinf ( (REAL4)LAL_TWOPI * in );
  *out2 = cosf ( (REAL4)LAL_TWOPI * in );
}

static inline REAL4 local_addf ( REAL4 x, REAL4 y ) {
  return x + y;
}

static inline REAL4 local_mulf ( REAL4 x, REAL4 y ) {
  return x * y;
}

// ========== internal generic functions ==========

// ---------- generic operator with 1 REAL4 vector input to 1 REAL4 vector output (S2S) ----------
static inline int
XLALVectorMath_S2S_GEN ( REAL4 *out, const REAL4 *in, const UINT4 len, REAL4 (*op)(REAL4) )
{
  for ( UINT4 i = 0; i < len; i ++ )
    {
      out[i] = (*op) ( in[i] );
    }
  return XLAL_SUCCESS;
}

// ---------- generic operator with 1 REAL4 vector input to 2 REAL4 vector outputs (S2SS) ----------
static inline int
XLALVectorMath_S2SS_GEN ( REAL4 *out1, REAL4 *out2, const REAL4 *in, const UINT4 len, void (*op)(REAL4, REAL4*, REAL4*) )
{
  for ( UINT4 i = 0; i < len; i ++ )
    {
      (*op) ( in[i], &(out1[i]), &(out2[i]) );
    }
  return XLAL_SUCCESS;
}

// ---------- generic operator with 2 REAL4 vector inputs to 1 REAL4 vector output (SS2S) ----------
static inline int
XLALVectorMath_SS2S_GEN ( REAL4 *out, const REAL4 *in1, const REAL4 *in2, const UINT4 len, REAL4 (*op)(REAL4, REAL4) )
{
  for ( UINT4 i = 0; i < len; i ++ )
    {
      out[i] = (*op) ( in1[i], in2[i] );
    }
  return XLAL_SUCCESS;
}

// ---------- generic operator with 1 REAL4 scalar and 1 REAL4 vector inputs to 1 REAL4 vector output (sS2S) ----------
static inline int
XLALVectorMath_sS2S_GEN ( REAL4 *out, REAL4 scalar, const REAL4 *in, const UINT4 len, REAL4 (*op)(REAL4, REAL4) )
{
  for ( UINT4 i = 0; i < len; i ++ )
    {
      out[i] = (*op) ( scalar, in[i] );
    }
  return XLAL_SUCCESS;
}

// ========== internal vector math functions ==========

// ---------- define vector math functions with 1 REAL4 vector input to 1 REAL4 vector output (S2S) ----------
#define DEFINE_VECTORMATH_S2S(NAME, GEN_OP)                             \
  DEFINE_VECTORMATH_ANY( XLALVectorMath_S2S_GEN, NAME ## REAL4, ( REAL4 *out, const REAL4 *in, const UINT4 len ), ( (out != NULL) && (in != NULL) ), ( out, in, len, GEN_OP ) )

DEFINE_VECTORMATH_S2S(Sin, sinf)
DEFINE_VECTORMATH_S2S(Cos, cosf)
DEFINE_VECTORMATH_S2S(Exp, expf)
DEFINE_VECTORMATH_S2S(Log, logf)

// ---------- define vector math functions with 1 REAL4 vector input to 2 REAL4 vector outputs (S2SS) ----------
#define DEFINE_VECTORMATH_S2SS(NAME, GEN_OP)                            \
  DEFINE_VECTORMATH_ANY( XLALVectorMath_S2SS_GEN, NAME ## REAL4, ( REAL4 *out1, REAL4 *out2, const REAL4 *in, const UINT4 len ), ( (out1 != NULL) && (out2 != NULL) && (in != NULL) ), ( out1, out2, in, len, GEN_OP ) )

DEFINE_VECTORMATH_S2SS(SinCos, local_sincosf)
DEFINE_VECTORMATH_S2SS(SinCos2Pi, local_sincosf_2pi)

// ---------- define vector math functions with 2 REAL4 vector inputs to 1 REAL4 vector output (SS2S) ----------
#define DEFINE_VECTORMATH_SS2S(NAME, GEN_OP)                            \
  DEFINE_VECTORMATH_ANY( XLALVectorMath_SS2S_GEN, NAME ## REAL4, ( REAL4 *out, const REAL4 *in1, const REAL4 *in2, const UINT4 len ), ( (out != NULL) && (in1 != NULL) && (in2 != NULL) ), ( out, in1, in2, len, GEN_OP ) )

DEFINE_VECTORMATH_SS2S(Add, local_addf)
DEFINE_VECTORMATH_SS2S(Multiply, local_mulf)
DEFINE_VECTORMATH_SS2S(Max, fmaxf)

// ---------- define vector math functions with 1 REAL4 scalar and 1 REAL4 vector inputs to 1 REAL4 vector output (sS2S) ----------
#define DEFINE_VECTORMATH_sS2S(NAME, GEN_OP)                            \
  DEFINE_VECTORMATH_ANY( XLALVectorMath_sS2S_GEN, NAME ## REAL4, ( REAL4 *out, REAL4 scalar, const REAL4 *in, const UINT4 len ), ( (out != NULL) && (in != NULL) ), ( out, scalar, in, len, GEN_OP ) )

DEFINE_VECTORMATH_sS2S(Shift, local_addf)
DEFINE_VECTORMATH_sS2S(Scale, local_mulf)
