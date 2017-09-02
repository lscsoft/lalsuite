/*
 * Copyright (C) 2015 Reinhard Prix, Karl Wette
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA  02111-1307  USA
 *
 */

/* ---------- internal macros ---------- */

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)

/* define internal SIMD-specific vector math functions, used by VectorMath_xxx.c sources */
#define DEFINE_VECTORMATH_ANY(GENERIC_FUNC, NAME, ARG_DEF, ARG_CHK, ARG_CALL) \
  int CONCAT2(XLALVector##NAME##_, SIMD_INSTRSET) ARG_DEF { \
    \
    XLAL_CHECK ( ( ARG_CHK ) , XLAL_EINVAL ); \
    \
    return GENERIC_FUNC ARG_CALL; \
    \
  }

/* ---------- internal prototypes of SIMD-specific vector math functions ---------- */

#define DECLARE_VECTORMATH_ANY(NAME, ARG_DEF, ISET1, ISET2, ISET3, ISET4) \
  extern const char* XLALVector##NAME##_name; \
  int XLALVector##NAME##_##ISET1 ARG_DEF; \
  int XLALVector##NAME##_##ISET2 ARG_DEF; \
  int XLALVector##NAME##_##ISET3 ARG_DEF; \
  int XLALVector##NAME##_##ISET4 ARG_DEF; \
  int XLALVector##NAME##_GEN     ARG_DEF;

/* declare internal prototypes of SIMD-specific vector math functions with 1 REAL4 vector input to 1 REAL4 vector output (S2S) */
#define DECLARE_VECTORMATH_S2S(NAME, ...)                                    \
  DECLARE_VECTORMATH_ANY( NAME ## REAL4, ( REAL4 *out, const REAL4 *in, const UINT4 len ), __VA_ARGS__ )

DECLARE_VECTORMATH_S2S(Sin, AVX2, AVX, SSE2, SSE)
DECLARE_VECTORMATH_S2S(Cos, AVX2, AVX, SSE2, SSE)
DECLARE_VECTORMATH_S2S(Exp, AVX2, AVX, SSE2, SSE)
DECLARE_VECTORMATH_S2S(Log, AVX2, AVX, SSE2, SSE)

/* declare internal prototypes of SIMD-specific vector math functions with 1 REAL4 vector input to 2 REAL4 vector outputs (S2SS) */
#define DECLARE_VECTORMATH_S2SS(NAME, ...)                                   \
  DECLARE_VECTORMATH_ANY( NAME ## REAL4, ( REAL4 *out1, REAL4 *out2, const REAL4 *in, const UINT4 len ), __VA_ARGS__ )

DECLARE_VECTORMATH_S2SS(SinCos, AVX2, AVX, SSE2, SSE)
DECLARE_VECTORMATH_S2SS(SinCos2Pi, AVX2, AVX, SSE2, SSE)

/* declare internal prototypes of SIMD-specific vector math functions with 2 REAL4 vector inputs to 1 REAL4 vector output (SS2S) */
#define DECLARE_VECTORMATH_SS2S(NAME, ...)                                   \
  DECLARE_VECTORMATH_ANY( NAME ## REAL4, ( REAL4 *out, const REAL4 *in1, const REAL4 *in2, const UINT4 len ), __VA_ARGS__ )

DECLARE_VECTORMATH_SS2S(Add, AVX2, AVX, SSE2, SSE)
DECLARE_VECTORMATH_SS2S(Multiply, AVX2, AVX, SSE2, SSE)
DECLARE_VECTORMATH_SS2S(Max, NONE, NONE, NONE, NONE)

/* declare internal prototypes of SIMD-specific vector math functions with 1 REAL4 scalar and 1 REAL4 vector input to 1 REAL4 vector output (sS2S) */
#define DECLARE_VECTORMATH_sS2S(NAME, ...) \
  DECLARE_VECTORMATH_ANY( NAME ## REAL4, ( REAL4 *out, REAL4 scalar, const REAL4 *in, const UINT4 len ), __VA_ARGS__ )

DECLARE_VECTORMATH_sS2S(Shift, AVX2, AVX, SSE2, SSE)
DECLARE_VECTORMATH_sS2S(Scale, AVX2, AVX, SSE2, SSE)

/* declare internal prototypes of SIMD-specific vector math functions with 2 REAL4 vector inputs to 1 UINT4 scalar and 1 UINT4 vector output (SS2uU) */
#define DECLARE_VECTORMATH_SS2uU(NAME, ...)                            \
  DECLARE_VECTORMATH_ANY( NAME ## REAL4, ( UINT4* count, UINT4 *out, const REAL4 *in1, const REAL4 *in2, const UINT4 len ), __VA_ARGS__ )

DECLARE_VECTORMATH_SS2uU(FindVectorLessEqual, AVX2, SSSE3, NONE, NONE)

/* declare internal prototypes of SIMD-specific vector math functions with 1 REAL4 scalar and 1 REAL4 vector inputs to 1 UINT4 scalar and 1 UINT4 vector output (sS2uU) */
#define DECLARE_VECTORMATH_sS2uU(NAME, ...)                            \
  DECLARE_VECTORMATH_ANY( NAME ## REAL4, ( UINT4* count, UINT4 *out, REAL4 scalar, const REAL4 *in, const UINT4 len ), __VA_ARGS__ )

DECLARE_VECTORMATH_sS2uU(FindScalarLessEqual, AVX2, SSSE3, NONE, NONE)


/* declare internal prototypes of SIMD-specific vector math functions with 1 REAL8 scalar and 1 REAL8 vector input to 1 REAL8 vector output (dD2D) */
#define DECLARE_VECTORMATH_dD2D(NAME, ...) \
  DECLARE_VECTORMATH_ANY( NAME ## REAL8, ( REAL8 *out, REAL8 scalar, const REAL8 *in, const UINT4 len ), __VA_ARGS__ )

DECLARE_VECTORMATH_dD2D(Scale, NONE, NONE, NONE, NONE)
