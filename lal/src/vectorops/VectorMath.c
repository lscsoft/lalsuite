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
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA  02111-1307  USA
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
/// Create a special REAL4 Vector with n-byte aligned memory \c data array.
///
/// This does not rely on \c posix_memalign() being available, and should compile+run everywhere.
/// Use XLALDestroyREAL4VectorAligned() to free this.
///
REAL4VectorAligned *
XLALCreateREAL4VectorAligned ( const UINT4 length, const UINT4 align )
{
  REAL4VectorAligned *ret;
  XLAL_CHECK_NULL ( (ret = XLALMalloc ( sizeof(*ret))) != NULL, XLAL_ENOMEM );

  ret->length = length;
  UINT4 paddedLength = length + align - 1;
  XLAL_CHECK_NULL ( (ret->data0 = XLALMalloc ( paddedLength * sizeof(REAL4) )) != NULL, XLAL_ENOMEM );

  size_t remBytes = ((size_t)ret->data0) % align;
  size_t offsetBytes = (align - remBytes) % align;
  ret->data = (void*)(((char*)ret->data0) + offsetBytes);

  XLAL_CHECK_NULL ( isMemAligned(ret->data,align), XLAL_EFAULT, "Failed to allocate %zd-byte aligned memory. Must be a coding error.\n", (size_t)align );

  return ret;
} // XLALCreateREAL4VectorAligned()

///
/// Destroy special n-byte aligned  REAL4VectorAligned
///
void
XLALDestroyREAL4VectorAligned ( REAL4VectorAligned *in )
{
  if ( !in ) {
    return;
  }
  if ( in->data0 ) {
    XLALFree ( in->data0 );
  }

  XLALFree ( in );

  return;

} // XLALDestroyREAL4VectorAligned()

// -------------------- export vector-operation functions --------------------

/* Declare the function pointer, define the dispatch function and exported vector math function */
#define EXPORT_VECTORMATH_ANY(NAME, ARG_DEF, ARG_CALL) \
  \
  static int XLALVector##NAME##_DISPATCH ARG_DEF; \
  \
  static int (*XLALVector##NAME##_ptr) ARG_DEF = XLALVector##NAME##_DISPATCH; \
  const char* XLALVector##NAME##_name = "\0"; \
  \
  int XLALVector##NAME##_DISPATCH ARG_DEF { \
    \
    DISPATCH_SELECT_BEGIN(); \
    DISPATCH_SELECT_AVX2(XLALVector##NAME##_ptr = XLALVector##NAME##_AVX2, XLALVector##NAME##_name = "XLALVector"#NAME"_AVX2"); \
    DISPATCH_SELECT_AVX( XLALVector##NAME##_ptr = XLALVector##NAME##_AVX,  XLALVector##NAME##_name = "XLALVector"#NAME"_AVX" ); \
    DISPATCH_SELECT_SSE2(XLALVector##NAME##_ptr = XLALVector##NAME##_SSE2, XLALVector##NAME##_name = "XLALVector"#NAME"_SSE2"); \
    DISPATCH_SELECT_SSE( XLALVector##NAME##_ptr = XLALVector##NAME##_SSE,  XLALVector##NAME##_name = "XLALVector"#NAME"_SSE" ); \
    DISPATCH_SELECT_END( XLALVector##NAME##_ptr = XLALVector##NAME##_GEN,  XLALVector##NAME##_name = "XLALVector"#NAME"_GEN" ); \
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

// ---------- define exported vector math functions with 1 REAL4 vector input to 1 REAL4 vector output (S2S) ----------
#define EXPORT_VECTORMATH_S2S(NAME)                                     \
  EXPORT_VECTORMATH_ANY( NAME ## REAL4, (REAL4 *out, const REAL4 *in, const UINT4 len), (out, in, len) )

EXPORT_VECTORMATH_S2S(Sin)
EXPORT_VECTORMATH_S2S(Cos)
EXPORT_VECTORMATH_S2S(Exp)
EXPORT_VECTORMATH_S2S(Log)

// ---------- define exported vector math functions with 1 REAL4 vector input to 2 REAL4 vector outputs (S2SS) ----------
#define EXPORT_VECTORMATH_S2SS(NAME)                                    \
  EXPORT_VECTORMATH_ANY( NAME ## REAL4, (REAL4 *out1, REAL4 *out2, const REAL4 *in, const UINT4 len), (out1, out2, in, len) )

EXPORT_VECTORMATH_S2SS(SinCos)
EXPORT_VECTORMATH_S2SS(SinCos2Pi)

// ---------- define exported vector math functions with 2 REAL4 vector inputs to 1 REAL4 vector output (SS2S) ----------
#define EXPORT_VECTORMATH_SS2S(NAME)                                    \
  EXPORT_VECTORMATH_ANY( NAME ## REAL4, (REAL4 *out, const REAL4 *in1, const REAL4 *in2, const UINT4 len), (out, in1, in2, len) )

EXPORT_VECTORMATH_SS2S(Add);
EXPORT_VECTORMATH_SS2S(Multiply);

// ---------- define exported vector math functions with 1 REAL4 scalar, 1 REAL4 vector inputs to 1 REAL4 vector output (sS2S) ----------
#define EXPORT_VECTORMATH_sS2S(NAME)                                    \
  EXPORT_VECTORMATH_ANY( NAME ## REAL4, (REAL4 *out, REAL4 scalar, const REAL4 *in, const UINT4 len), (out, scalar, in, len) )

EXPORT_VECTORMATH_sS2S(Scale);
EXPORT_VECTORMATH_sS2S(Shift);

