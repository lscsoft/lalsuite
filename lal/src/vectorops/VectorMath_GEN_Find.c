//
// Copyright (C) 2017 Karl Wette
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
#include <math.h>
#include <config.h>

#include <lal/LALConstants.h>
#include <lal/VectorMath.h>

#include "VectorMath_internal.h"

#define SIMD_INSTRSET GEN

// ---------- local operators and operator-wrappers ----------
static inline int local_cmplef ( REAL4 in1, REAL4 in2 ) {
  return in1 <= in2;
}

// ========== internal generic functions ==========

// ---------- generic operator with 2 REAL4 vector inputs to 1 UINT4 scalar and 1 UINT4 vector output (SS2uU) ----------
static inline int
XLALVectorMath_SS2uU_GEN ( UINT4* count, UINT4 *out, const REAL4 *in1, const REAL4 *in2, const UINT4 len, int (*pred)(REAL4, REAL4) )
{
  *count = 0;
  for (UINT4 i = 0; i < len; i++) {
    if ( (*pred)(in1[i], in2[i]) ) {
      out[*count] = i;
      *count += 1;
    }
  }
  return XLAL_SUCCESS;
}

// ---------- generic operator with 1 REAL4 scalar and 1 REAL4 vector inputs to 1 UINT4 scalar and 1 UINT4 vector output (Ss2uU) ----------
static inline int
XLALVectorMath_sS2uU_GEN ( UINT4* count, UINT4 *out, REAL4 scalar, const REAL4 *in, const UINT4 len, int (*pred)(REAL4, REAL4) )
{
  *count = 0;
  for (UINT4 i = 0; i < len; i++) {
    if ( (*pred)(scalar, in[i]) ) {
      out[*count] = i;
      *count += 1;
    }
  }
  return XLAL_SUCCESS;
}

// ========== internal vector math functions ==========

// ---------- define vector math functions with 2 REAL4 vector inputs to 1 UINT4 scalar and 1 UINT4 vector output (SS2uU) ----------
#define DEFINE_VECTORMATH_SS2uU(NAME, PRED)                            \
  DEFINE_VECTORMATH_ANY( XLALVectorMath_SS2uU_GEN, NAME ## REAL4, ( UINT4* count, UINT4 *out, const REAL4 *in1, const REAL4 *in2, const UINT4 len ), ( (count != NULL) && (out != NULL) && (in1 != NULL) && (in2 != NULL) ), ( count, out, in1, in2, len, PRED ) )

DEFINE_VECTORMATH_SS2uU(FindVectorLessEqual, local_cmplef)

// ---------- define vector math functions with 1 REAL4 scalar and 1 REAL4 vector inputs to 1 UINT4 scalar and 1 UINT4 vector output (sS2uU) ----------
#define DEFINE_VECTORMATH_sS2uU(NAME, PRED)                            \
  DEFINE_VECTORMATH_ANY( XLALVectorMath_sS2uU_GEN, NAME ## REAL4, ( UINT4* count, UINT4 *out, REAL4 scalar, const REAL4 *in, const UINT4 len ), ( (count != NULL) && (out != NULL) && (in != NULL) ), ( count, out, scalar, in, len, PRED ) )

DEFINE_VECTORMATH_sS2uU(FindScalarLessEqual, local_cmplef)
