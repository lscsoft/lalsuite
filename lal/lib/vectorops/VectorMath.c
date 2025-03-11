//
// Copyright (C) 2015 Reinhard Prix, Karl Wette
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
// Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
// MA  02110-1301  USA
//

// ---------- INCLUDES ----------
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include <config.h>
#include <simd_dispatch.h>

#include <lal/LALString.h>
#include <lal/LALConstants.h>
#include <lal/VectorMath.h>

#include "VectorMath_internal.h"

//==================== FUNCTION DEFINITIONS ====================*/

// -------------------- our own failsafe aligned memory handling --------------------

///
/// Create a special \<TYPE\>Vector with n-byte aligned memory \c data array.
///
/// This does not rely on \c posix_memalign() being available, and should compile+run everywhere.
/// Use XLALDestroy\<TYPE\>VectorAligned() to free this.
///
#define DEFINE_ALIGNED_VECT_API(TYPE)                                   \
TYPE##VectorAligned *XLALCreate##TYPE##VectorAligned ( const UINT4 length, const UINT4 align ) \
{                                                                       \
 TYPE##VectorAligned *ret;                                              \
 XLAL_CHECK_NULL ( (ret = XLALCalloc ( 1, sizeof(*ret) )) != NULL, XLAL_ENOMEM ); \
                                                                        \
 XLAL_CHECK_NULL ( (ret = XLALResize##TYPE##VectorAligned ( ret, length, align )) != NULL, XLAL_ENOMEM ); \
                                                                        \
 return ret;                                                            \
} /* XLALCreate\<TYPE\>VectorAligned() */                               \
                                                                        \
TYPE##VectorAligned *XLALResize##TYPE##VectorAligned ( TYPE##VectorAligned *in, const UINT4 length, const UINT4 align ) \
{                                                                       \
 if ( in == NULL ) {                                                    \
  return XLALCreate##TYPE##VectorAligned ( length, align );             \
 }                                                                      \
 if ( length == 0 ) {                                                   \
  XLALDestroy##TYPE##VectorAligned ( in );                              \
  return NULL;                                                          \
 }                                                                      \
                                                                        \
 in->length = length;                                                   \
 UINT4 paddedLength = length + align - 1;                               \
 XLAL_CHECK_NULL ( (in->data0 = XLALRealloc ( in->data0, paddedLength * sizeof(in->data0[0]) )) != NULL, XLAL_ENOMEM ); \
                                                                        \
 size_t remBytes = ((size_t)in->data0) % align;                         \
 size_t offsetBytes = (align - remBytes) % align;                       \
 in->data = (void*)(((char*)in->data0) + offsetBytes);                  \
                                                                        \
 XLAL_CHECK_NULL ( ((size_t)in->data) % align == 0, XLAL_EFAULT, "Failed to allocate %zd-byte aligned memory. Must be a coding error.\n", (size_t)align ); \
                                                                        \
 return in;                                                             \
} /* XLALResize\<TYPE\>VectorAligned() */                               \
                                                                        \
void XLALDestroy##TYPE##VectorAligned ( TYPE##VectorAligned *in )       \
{                                                                       \
  if ( !in ) { return; }                                                \
  if ( in->data0 ) {                                                    \
    XLALFree ( in->data0 );                                             \
  }                                                                     \
  XLALFree ( in );                                                      \
  return;                                                               \
} /* XLALDestroy\<TYPE\>VectorAligned() */

DEFINE_ALIGNED_VECT_API(UINT4);
DEFINE_ALIGNED_VECT_API(REAL4);
DEFINE_ALIGNED_VECT_API(REAL8);
DEFINE_ALIGNED_VECT_API(COMPLEX8);
DEFINE_ALIGNED_VECT_API(COMPLEX16);

// -------------------- export vector-operation functions --------------------

/* Declare the function pointer, define the dispatch function, and export vector math function with supported instruction sets */
#define EXPORT_VECTORMATH_ANY(NAME, ARG_DEF, ARG_CALL, ISET1, ISET2, ISET3, ISET4) \
  \
  static int XLALVector##NAME##_DISPATCH ARG_DEF; \
  \
  static int (*XLALVector##NAME##_ptr) ARG_DEF = XLALVector##NAME##_DISPATCH; \
  const char* XLALVector##NAME##_name = "\0"; \
  \
  int XLALVector##NAME##_DISPATCH ARG_DEF { \
    \
    DISPATCH_SELECT_BEGIN(); \
    CONCAT2(DISPATCH_SELECT_,ISET1)(XLALVector##NAME##_ptr = XLALVector##NAME##_##ISET1, XLALVector##NAME##_name = "XLALVector"#NAME"_"#ISET1); \
    CONCAT2(DISPATCH_SELECT_,ISET2)(XLALVector##NAME##_ptr = XLALVector##NAME##_##ISET2, XLALVector##NAME##_name = "XLALVector"#NAME"_"#ISET2); \
    CONCAT2(DISPATCH_SELECT_,ISET3)(XLALVector##NAME##_ptr = XLALVector##NAME##_##ISET3, XLALVector##NAME##_name = "XLALVector"#NAME"_"#ISET3); \
    CONCAT2(DISPATCH_SELECT_,ISET4)(XLALVector##NAME##_ptr = XLALVector##NAME##_##ISET4, XLALVector##NAME##_name = "XLALVector"#NAME"_"#ISET4); \
    DISPATCH_SELECT_END( XLALVector##NAME##_ptr = XLALVector##NAME##_GEN,     XLALVector##NAME##_name = "XLALVector"#NAME"_GEN"   ); \
    \
    return XLALVector##NAME ARG_CALL; \
    \
  } \
  \
  int XLALVector##NAME ARG_DEF { \
    \
    return (XLALVector##NAME##_ptr) ARG_CALL; \
    \
  }

// ---------- define exported vector math functions with 1 REAL4 vector input to 1 INT4 vector output (S2I) ----------
#define EXPORT_VECTORMATH_S2I(NAME, ...)                                     \
  EXPORT_VECTORMATH_ANY( NAME ## REAL4, (INT4 *out, const REAL4 *in, const UINT4 len), (out, in, len), __VA_ARGS__ )

EXPORT_VECTORMATH_S2I(INT4From, SSE2, NONE, NONE, NONE)

// ---------- define exported vector math functions with 1 REAL4 vector input to 1 REAL4 scalar output (S2s) ----------
#define EXPORT_VECTORMATH_S2s(NAME, ...)                                     \
  EXPORT_VECTORMATH_ANY( NAME ## REAL4, (REAL4 *out, const REAL4 *in, const UINT4 len), (out, in, len), __VA_ARGS__ )

EXPORT_VECTORMATH_S2s(ScalarMax, AVX2, AVX, SSE2, NONE)

// ---------- define exported vector math functions with 1 REAL4 vector input to 1 REAL4 vector output (S2S) ----------
#define EXPORT_VECTORMATH_S2S(NAME, ...)                                     \
  EXPORT_VECTORMATH_ANY( NAME ## REAL4, (REAL4 *out, const REAL4 *in, const UINT4 len), (out, in, len), __VA_ARGS__ )

EXPORT_VECTORMATH_S2S(Sin, AVX2, AVX, SSE2, NONE)
EXPORT_VECTORMATH_S2S(Cos, AVX2, AVX, SSE2, NONE)
EXPORT_VECTORMATH_S2S(Exp, AVX2, AVX, SSE2, NONE)
EXPORT_VECTORMATH_S2S(Log, AVX2, AVX, SSE2, NONE)
EXPORT_VECTORMATH_S2S(Round, AVX2, AVX, NONE, NONE)

// ---------- define exported vector math functions with 1 REAL4 vector input to 2 REAL4 vector outputs (S2SS) ----------
#define EXPORT_VECTORMATH_S2SS(NAME, ...)                                    \
  EXPORT_VECTORMATH_ANY( NAME ## REAL4, (REAL4 *out1, REAL4 *out2, const REAL4 *in, const UINT4 len), (out1, out2, in, len), __VA_ARGS__ )

EXPORT_VECTORMATH_S2SS(SinCos, AVX2, AVX, SSE2, NONE)
EXPORT_VECTORMATH_S2SS(SinCos2Pi, AVX2, AVX, SSE2, NONE)

// ---------- define exported vector math functions with 2 REAL4 vector inputs to 1 REAL4 vector output (SS2S) ----------
#define EXPORT_VECTORMATH_SS2S(NAME, ...)                                    \
  EXPORT_VECTORMATH_ANY( NAME ## REAL4, (REAL4 *out, const REAL4 *in1, const REAL4 *in2, const UINT4 len), (out, in1, in2, len), __VA_ARGS__ )

EXPORT_VECTORMATH_SS2S(Add, AVX2, AVX, SSE2, NONE)
EXPORT_VECTORMATH_SS2S(Sub, AVX2, AVX, SSE2, NONE)
EXPORT_VECTORMATH_SS2S(Multiply, AVX2, AVX, SSE2, NONE)
EXPORT_VECTORMATH_SS2S(Max, AVX2, AVX, SSE2, NONE)

// ---------- define exported vector math functions with 1 REAL4 scalar, 1 REAL4 vector inputs to 1 REAL4 vector output (sS2S) ----------
#define EXPORT_VECTORMATH_sS2S(NAME, ...)                                    \
  EXPORT_VECTORMATH_ANY( NAME ## REAL4, (REAL4 *out, REAL4 scalar, const REAL4 *in, const UINT4 len), (out, scalar, in, len), __VA_ARGS__ )

EXPORT_VECTORMATH_sS2S(Scale, AVX2, AVX, SSE2, NONE)
EXPORT_VECTORMATH_sS2S(Shift, AVX2, AVX, SSE2, NONE)

// ---------- define exported vector math functions with 1 REAL4 scalar, 2 REAL4 vector inputs to 1 REAL4 vector output (sSS2S) ----------
#define EXPORT_VECTORMATH_sSS2S(NAME, ...)                                   \
  EXPORT_VECTORMATH_ANY( NAME ## REAL4, (REAL4 *out, REAL4 scalar, const REAL4 *in1, const REAL4 *in2, const UINT4 len), (out, scalar, in1, in2, len), __VA_ARGS__ )

EXPORT_VECTORMATH_sSS2S(ScaleAdd, AVX2, AVX, SSE2, NONE)

// ---------- define exported vector math functions with 2 REAL4 vector inputs to 1 UINT4 scalar and 1 UINT4 vector output (SS2uU) ----------
#define EXPORT_VECTORMATH_SS2uU(NAME, ...)                            \
  EXPORT_VECTORMATH_ANY( NAME ## REAL4, ( UINT4* count, UINT4 *out, const REAL4 *in1, const REAL4 *in2, const UINT4 len ), (count, out, in1, in2, len), __VA_ARGS__ )

EXPORT_VECTORMATH_SS2uU(FindVectorLessEqual, AVX2, SSSE3, NONE, NONE)

// ---------- define exported vector math functions with 1 REAL4 scalar and 1 REAL4 vector inputs to 1 UINT4 scalar and 1 UINT4 vector output (sS2uU) ----------
#define EXPORT_VECTORMATH_sS2uU(NAME, ...)                            \
  EXPORT_VECTORMATH_ANY( NAME ## REAL4, ( UINT4* count, UINT4 *out, REAL4 scalar, const REAL4 *in, const UINT4 len ), (count, out, scalar, in, len), __VA_ARGS__ )

EXPORT_VECTORMATH_sS2uU(FindScalarLessEqual, AVX2, SSSE3, NONE, NONE)

// ---------- define exported vector math functions with 1 REAL8 scalar, 1 REAL8 vector inputs to 1 REAL8 vector output (dD2D) ----------
#define EXPORT_VECTORMATH_dD2D(NAME, ...)                                    \
  EXPORT_VECTORMATH_ANY( NAME ## REAL8, (REAL8 *out, REAL8 scalar, const REAL8 *in, const UINT4 len), (out, scalar, in, len), __VA_ARGS__ )

EXPORT_VECTORMATH_dD2D(Scale, AVX2, AVX, SSE2, NONE)
EXPORT_VECTORMATH_dD2D(Shift, AVX2, AVX, SSE2, NONE)

// ---------- define exported vector math functions with 1 REAL8 scalar, 2 REAL8 vector inputs to 1 REAL8 vector output (dDD2D) ----------
#define EXPORT_VECTORMATH_dDD2D(NAME, ...)                                   \
  EXPORT_VECTORMATH_ANY( NAME ## REAL8, (REAL8 *out, REAL8 scalar, const REAL8 *in1, const REAL8 *in2, const UINT4 len), (out, scalar, in1, in2, len), __VA_ARGS__ )

EXPORT_VECTORMATH_dDD2D(ScaleAdd, AVX2, AVX, SSE2, NONE)

// ---------- define exported vector math functions with 2 REAL8 vector inputs to 1 REAL8 vector output (DD2D) ----------
#define EXPORT_VECTORMATH_DD2D(NAME, ...)                                    \
  EXPORT_VECTORMATH_ANY( NAME ## REAL8, (REAL8 *out, const REAL8 *in1, const REAL8 *in2, const UINT4 len), (out, in1, in2, len), __VA_ARGS__ )

EXPORT_VECTORMATH_DD2D(Add, AVX2, AVX, SSE2, NONE)
EXPORT_VECTORMATH_DD2D(Sub, AVX2, AVX, SSE2, NONE)
EXPORT_VECTORMATH_DD2D(Multiply, AVX2, AVX, SSE2, NONE)
EXPORT_VECTORMATH_DD2D(Max, AVX2, AVX, SSE2, NONE)

// ---------- define exported vector math functions with 2 COMPLEX8 vector inputs to 1 COMPLEX8 vector output (CC2C) ----------
#define EXPORT_VECTORMATH_CC2C(NAME, ...)                                    \
  EXPORT_VECTORMATH_ANY( NAME ## COMPLEX8, (COMPLEX8 *out, const COMPLEX8 *in1, const COMPLEX8 *in2, const UINT4 len), (out, in1, in2, len), __VA_ARGS__ )

EXPORT_VECTORMATH_CC2C(Multiply, AVX2, AVX, SSE2, NONE)
EXPORT_VECTORMATH_CC2C(Add, AVX2, AVX, SSE2, NONE)

// ---------- define exported vector math functions with 1 COMPLEX8 scalar and 1 COMPLEX8 vector inputs to 1 COMPLEX8 vector output (cC2C) ----------
#define EXPORT_VECTORMATH_cC2C(NAME, ...)                                    \
  EXPORT_VECTORMATH_ANY( NAME ## COMPLEX8, (COMPLEX8 *out, COMPLEX8 scalar, const COMPLEX8 *in, const UINT4 len), (out, scalar, in, len), __VA_ARGS__ )

EXPORT_VECTORMATH_cC2C(Scale, AVX2, AVX, SSE2, NONE)
EXPORT_VECTORMATH_cC2C(Shift, AVX2, AVX, SSE2, NONE)

// ---------- define exported vector math functions with 1 REAL4 scalar, 2 COMPLEX8 vector inputs to 1 COMPLEX8 vector output (sCC2C) ----------
#define EXPORT_VECTORMATH_sCC2C(NAME, ...)                                   \
  EXPORT_VECTORMATH_ANY( NAME ## COMPLEX8, (COMPLEX8 *out, REAL4 scalar, const COMPLEX8 *in1, const COMPLEX8 *in2, const UINT4 len), (out, scalar, in1, in2, len), __VA_ARGS__ )

EXPORT_VECTORMATH_sCC2C(ScaleAdd, AVX2, AVX, SSE2, NONE)

// ---------- define exported vector math functions with 1 REAL8 vector input to 1 REAL8 scalar output (D2d) ----------
#define EXPORT_VECTORMATH_D2d(NAME, ...)                                     \
  EXPORT_VECTORMATH_ANY( NAME ## REAL8, (REAL8 *out, const REAL8 *in, const UINT4 len), (out, in, len), __VA_ARGS__ )

EXPORT_VECTORMATH_D2d(ScalarMax, AVX2, AVX, SSE2, NONE)

// ---------- define exported vector math functions with 1 REAL8 vector input to 1 REAL8 vector output (D2D) ----------
#define EXPORT_VECTORMATH_D2D(NAME, ...)                                     \
  EXPORT_VECTORMATH_ANY( NAME ## REAL8, (REAL8 *out, const REAL8 *in, const UINT4 len), (out, in, len), __VA_ARGS__ )

EXPORT_VECTORMATH_D2D(Round, AVX2, AVX, NONE, NONE)
