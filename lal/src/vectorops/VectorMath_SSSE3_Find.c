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

#include <xmmintrin.h>
#include <tmmintrin.h>

/* yes I know, the top of this file is quite ugly */
#ifdef _MSC_VER /* visual c++ */
# define ALIGN16_BEG __declspec(align(16))
# define ALIGN16_END
#else /* gcc or icc */
# define ALIGN16_BEG
# define ALIGN16_END __attribute__((aligned(16)))
#endif

/* __m128 is ugly to write */
typedef __m128 v4sf;  // vector of 4 float (sse1)
typedef __m128i v4si; // vector of 4 int (sse2)

typedef ALIGN16_BEG union {
  float f[4];
  int i[4];
  v4sf  v;
  v4si vi;
} ALIGN16_END V4SF;

// ---------- local constants ----------

// Given a bitmask 'mask', 'perm[mask]' returns the permutation required to move all masked items before all unmasked items
static const UINT4 permutation[16][4] = {
  {0x3020100, 0x7060504, 0xb0a0908, 0xf0e0d0c}, {0x3020100, 0x7060504, 0xb0a0908, 0xf0e0d0c}, {0x7060504, 0x3020100, 0xb0a0908, 0xf0e0d0c}, {0x3020100, 0x7060504, 0xb0a0908, 0xf0e0d0c},
  {0xb0a0908, 0x3020100, 0x7060504, 0xf0e0d0c}, {0x3020100, 0xb0a0908, 0x7060504, 0xf0e0d0c}, {0x7060504, 0xb0a0908, 0x3020100, 0xf0e0d0c}, {0x3020100, 0x7060504, 0xb0a0908, 0xf0e0d0c},
  {0xf0e0d0c, 0x3020100, 0x7060504, 0xb0a0908}, {0x3020100, 0xf0e0d0c, 0x7060504, 0xb0a0908}, {0x7060504, 0xf0e0d0c, 0x3020100, 0xb0a0908}, {0x3020100, 0x7060504, 0xf0e0d0c, 0xb0a0908},
  {0xb0a0908, 0xf0e0d0c, 0x3020100, 0x7060504}, {0x3020100, 0xb0a0908, 0xf0e0d0c, 0x7060504}, {0x7060504, 0xb0a0908, 0xf0e0d0c, 0x3020100}, {0x3020100, 0x7060504, 0xb0a0908, 0xf0e0d0c},
};

// Given a bitmask 'mask', 'popcount[mask]' returns its population count, i.e. the number of bits set
static const UINT4 popcount[16] = {
  0, 0x1, 0x1, 0x2,
  0x1, 0x2, 0x2, 0x3,
  0x1, 0x2, 0x2, 0x3,
  0x2, 0x3, 0x3, 0x4,
};

// ---------- local operators and operator-wrappers ----------
UNUSED static inline __m128
local_cmple_ps ( __m128 in1, __m128 in2 )
{
  return _mm_cmple_ps ( in1, in2 );
}

// ========== internal generic SSSE3 functions ==========

// ---------- generic SSSE3 operator with 2 REAL4 vector inputs to 1 UINT4 scalar and 1 UINT4 vector output (SS2uU) ----------
static inline int
XLALVectorMath_SS2uU_SSSE3 ( UINT4* count, UINT4 *out, const REAL4 *in1, const REAL4 *in2, const UINT4 len, __m128 (*pred)(__m128, __m128) )
{
  *count = 0;

  // walk through vector in blocks of 4
  UINT4 i4Max = len - ( len % 4 );
  for ( UINT4 i4 = 0; i4 < i4Max; i4 += 4 )
    {
      // load vector inputs
      __m128 in4p_1 = _mm_loadu_ps(&in1[i4]);
      __m128 in4p_2 = _mm_loadu_ps(&in2[i4]);
      // 'pred' should set each float to 0xFFFFFFFF if satisfied, 0x0 otherwise
      __m128 res = (*pred)(in4p_1, in4p_2);
      // get bitmask from 'r', bits are set if 'pred' is satisfied
      int mask = _mm_movemask_ps(res);
      // get correct permutation so that masked items are moved before unmasked items
      __m128i perm = _mm_set_epi32(permutation[mask][3], permutation[mask][2], permutation[mask][1], permutation[mask][0]);
      // create packed int of indexes to current vector inputs
      // - i.e. idx = i4 + (3, 2, 1, 0)
      __m128i idx = _mm_add_epi32(_mm_set1_epi32(i4), _mm_set_epi32(3, 2, 1, 0));
      // permute so that indexes of items for which 'pred' is satisfied appear first
      __m128i perm_idx = _mm_shuffle_epi8(idx, perm);
      // store indexes in output
      // - indexes of items for which 'pred' was not satisfied will be erased on the next loop
      _mm_storeu_si128((void *) &out[*count], perm_idx);
      // increment count by number of items for which 'pred' is satisfied
      *count += popcount[mask];
    }

  // deal with the remaining (<=3) terms separately
  V4SF in4_1 = {.f={0,0,0,0}};
  V4SF in4_2 = {.f={0,0,0,0}};
  V4SF out4;
  for (UINT4 i = i4Max, j = 0; i < len; i++, j++) {
    in4_1.f[j] = in1[i];
    in4_2.f[j] = in2[i];
  }
  out4.v = (*pred)(in4_1.v, in4_2.v);
  for (UINT4 i = i4Max, j = 0; i < len; i++, j++) {
    if (out4.i[j]) {
      out[*count] = i;
      *count += 1;
    }
  }

  return XLAL_SUCCESS;

} // XLALVectorMath_SS2uU_SSSE3()

// ---------- generic SSSE3 operator with 1 REAL4 scalar and 1 REAL4 vector inputs to 1 UINT4 scalar and 1 UINT4 vector output (Ss2uU) ----------
static inline int
XLALVectorMath_sS2uU_SSSE3 ( UINT4* count, UINT4 *out, REAL4 scalar, const REAL4 *in, const UINT4 len, __m128 (*pred)(__m128, __m128) )
{
  const V4SF scalar4 = {.f={scalar,scalar,scalar,scalar}};

  *count = 0;

  // walk through vector in blocks of 4
  UINT4 i4Max = len - ( len % 4 );
  for ( UINT4 i4 = 0; i4 < i4Max; i4 += 4 )
    {
      // load vector inputs
      __m128 in4p = _mm_loadu_ps(&in[i4]);
      // 'pred' should set each float to 0xFFFFFFFF if satisfied, 0x0 otherwise
      __m128 res = (*pred)(scalar4.v, in4p);
      // get bitmask from 'r', bits are set if 'pred' is satisfied
      int mask = _mm_movemask_ps(res);
      // get correct permutation so that masked items are moved before unmasked items
      __m128i perm = _mm_set_epi32(permutation[mask][3], permutation[mask][2], permutation[mask][1], permutation[mask][0]);
      // create packed int of indexes to current vector inputs
      // - i.e. idx = i4 + (3, 2, 1, 0)
      __m128i idx = _mm_add_epi32(_mm_set1_epi32(i4), _mm_set_epi32(3, 2, 1, 0));
      // permute so that indexes of items for which 'pred' is satisfied appear first
      __m128i perm_idx = _mm_shuffle_epi8(idx, perm);
      // store indexes in output
      // - indexes of items for which 'pred' was not satisfied will be erased on the next loop
      _mm_storeu_si128((void *) &out[*count], perm_idx);
      // increment count by number of items for which 'pred' is satisfied
      *count += popcount[mask];
    }

  // deal with the remaining (<=3) terms separately
  V4SF in4 = {.f={0,0,0,0}};
  V4SF out4;
  for (UINT4 i = i4Max, j = 0; i < len; i++, j++) {
    in4.f[j] = in[i];
  }
  out4.v = (*pred)(scalar4.v, in4.v);
  for (UINT4 i = i4Max, j = 0; i < len; i++, j++) {
    if (out4.i[j]) {
      out[*count] = i;
      *count += 1;
    }
  }

  return XLAL_SUCCESS;

} // XLALVectorMath_sS2uU_SSSE3()

// ========== internal SSSE3 vector math functions ==========

// ---------- define vector math functions with 2 REAL4 vector inputs to 1 UINT4 scalar and 1 UINT4 vector output (SS2uU) ----------
#define DEFINE_VECTORMATH_SS2uU(NAME, PRED)                            \
  DEFINE_VECTORMATH_ANY( XLALVectorMath_SS2uU_SSSE3, NAME ## REAL4, ( UINT4* count, UINT4 *out, const REAL4 *in1, const REAL4 *in2, const UINT4 len ), ( (count != NULL) && (out != NULL) && (in1 != NULL) && (in2 != NULL) ), ( count, out, in1, in2, len, PRED ) )

DEFINE_VECTORMATH_SS2uU(FindVectorLessEqual, local_cmple_ps)

// ---------- define vector math functions with 1 REAL4 scalar and 1 REAL4 vector inputs to 1 UINT4 scalar and 1 UINT4 vector output (sS2uU) ----------
#define DEFINE_VECTORMATH_sS2uU(NAME, PRED)                            \
  DEFINE_VECTORMATH_ANY( XLALVectorMath_sS2uU_SSSE3, NAME ## REAL4, ( UINT4* count, UINT4 *out, REAL4 scalar, const REAL4 *in, const UINT4 len ), ( (count != NULL) && (out != NULL) && (in != NULL) ), ( count, out, scalar, in, len, PRED ) )

DEFINE_VECTORMATH_sS2uU(FindScalarLessEqual, local_cmple_ps)
