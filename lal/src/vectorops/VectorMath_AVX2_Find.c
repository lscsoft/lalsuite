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

#include <immintrin.h>

/* yes I know, the top of this file is quite ugly */
#ifdef _MSC_VER /* visual c++ */
# define ALIGN16_BEG __declspec(align(32))
# define ALIGN16_END
#else /* gcc or icc */
# define ALIGN32_BEG
# define ALIGN32_END __attribute__((aligned(32)))
#endif

/* __m128 is ugly to write */
typedef __m256  v8sf; // vector of 8 float (avx)
typedef __m256i v8si; // vector of 8 int   (avx)
typedef __m128i v4sii; // vector of 8 int   (avx)

typedef ALIGN32_BEG union {
  float f[8];
  int i[8];
  v8sf  v;
  v8si  vi;
} ALIGN32_END V8SF;

// ---------- local constants ----------

// Given a bitmask 'mask', 'shuffle[mask]' returns the shuffling operation required to move all masked items before all unmasked items
static const UINT4 shuffle[256] = {
  0x76543210, 0x76543210, 0x76543201, 0x76543210, 0x76543102, 0x76543120, 0x76543021, 0x76543210, 0x76542103, 0x76542130, 0x76542031, 0x76542310, 0x76541032, 0x76541320, 0x76540321, 0x76543210,
  0x76532104, 0x76532140, 0x76532041, 0x76532410, 0x76531042, 0x76531420, 0x76530421, 0x76534210, 0x76521043, 0x76521430, 0x76520431, 0x76524310, 0x76510432, 0x76514320, 0x76504321, 0x76543210,
  0x76432105, 0x76432150, 0x76432051, 0x76432510, 0x76431052, 0x76431520, 0x76430521, 0x76435210, 0x76421053, 0x76421530, 0x76420531, 0x76425310, 0x76410532, 0x76415320, 0x76405321, 0x76453210,
  0x76321054, 0x76321540, 0x76320541, 0x76325410, 0x76310542, 0x76315420, 0x76305421, 0x76354210, 0x76210543, 0x76215430, 0x76205431, 0x76254310, 0x76105432, 0x76154320, 0x76054321, 0x76543210,
  0x75432106, 0x75432160, 0x75432061, 0x75432610, 0x75431062, 0x75431620, 0x75430621, 0x75436210, 0x75421063, 0x75421630, 0x75420631, 0x75426310, 0x75410632, 0x75416320, 0x75406321, 0x75463210,
  0x75321064, 0x75321640, 0x75320641, 0x75326410, 0x75310642, 0x75316420, 0x75306421, 0x75364210, 0x75210643, 0x75216430, 0x75206431, 0x75264310, 0x75106432, 0x75164320, 0x75064321, 0x75643210,
  0x74321065, 0x74321650, 0x74320651, 0x74326510, 0x74310652, 0x74316520, 0x74306521, 0x74365210, 0x74210653, 0x74216530, 0x74206531, 0x74265310, 0x74106532, 0x74165320, 0x74065321, 0x74653210,
  0x73210654, 0x73216540, 0x73206541, 0x73265410, 0x73106542, 0x73165420, 0x73065421, 0x73654210, 0x72106543, 0x72165430, 0x72065431, 0x72654310, 0x71065432, 0x71654320, 0x70654321, 0x76543210,
  0x65432107, 0x65432170, 0x65432071, 0x65432710, 0x65431072, 0x65431720, 0x65430721, 0x65437210, 0x65421073, 0x65421730, 0x65420731, 0x65427310, 0x65410732, 0x65417320, 0x65407321, 0x65473210,
  0x65321074, 0x65321740, 0x65320741, 0x65327410, 0x65310742, 0x65317420, 0x65307421, 0x65374210, 0x65210743, 0x65217430, 0x65207431, 0x65274310, 0x65107432, 0x65174320, 0x65074321, 0x65743210,
  0x64321075, 0x64321750, 0x64320751, 0x64327510, 0x64310752, 0x64317520, 0x64307521, 0x64375210, 0x64210753, 0x64217530, 0x64207531, 0x64275310, 0x64107532, 0x64175320, 0x64075321, 0x64753210,
  0x63210754, 0x63217540, 0x63207541, 0x63275410, 0x63107542, 0x63175420, 0x63075421, 0x63754210, 0x62107543, 0x62175430, 0x62075431, 0x62754310, 0x61075432, 0x61754320, 0x60754321, 0x67543210,
  0x54321076, 0x54321760, 0x54320761, 0x54327610, 0x54310762, 0x54317620, 0x54307621, 0x54376210, 0x54210763, 0x54217630, 0x54207631, 0x54276310, 0x54107632, 0x54176320, 0x54076321, 0x54763210,
  0x53210764, 0x53217640, 0x53207641, 0x53276410, 0x53107642, 0x53176420, 0x53076421, 0x53764210, 0x52107643, 0x52176430, 0x52076431, 0x52764310, 0x51076432, 0x51764320, 0x50764321, 0x57643210,
  0x43210765, 0x43217650, 0x43207651, 0x43276510, 0x43107652, 0x43176520, 0x43076521, 0x43765210, 0x42107653, 0x42176530, 0x42076531, 0x42765310, 0x41076532, 0x41765320, 0x40765321, 0x47653210,
  0x32107654, 0x32176540, 0x32076541, 0x32765410, 0x31076542, 0x31765420, 0x30765421, 0x37654210, 0x21076543, 0x21765430, 0x20765431, 0x27654310, 0x10765432, 0x17654320, 0x7654321, 0x76543210,
};

// Given a bitmask 'mask', 'popcount[mask]' returns its population count, i.e. the number of bits set
static const UINT4 popcount[256] = {
  0, 0x1, 0x1, 0x2, 0x1, 0x2, 0x2, 0x3, 0x1, 0x2, 0x2, 0x3, 0x2, 0x3, 0x3, 0x4,
  0x1, 0x2, 0x2, 0x3, 0x2, 0x3, 0x3, 0x4, 0x2, 0x3, 0x3, 0x4, 0x3, 0x4, 0x4, 0x5,
  0x1, 0x2, 0x2, 0x3, 0x2, 0x3, 0x3, 0x4, 0x2, 0x3, 0x3, 0x4, 0x3, 0x4, 0x4, 0x5,
  0x2, 0x3, 0x3, 0x4, 0x3, 0x4, 0x4, 0x5, 0x3, 0x4, 0x4, 0x5, 0x4, 0x5, 0x5, 0x6,
  0x1, 0x2, 0x2, 0x3, 0x2, 0x3, 0x3, 0x4, 0x2, 0x3, 0x3, 0x4, 0x3, 0x4, 0x4, 0x5,
  0x2, 0x3, 0x3, 0x4, 0x3, 0x4, 0x4, 0x5, 0x3, 0x4, 0x4, 0x5, 0x4, 0x5, 0x5, 0x6,
  0x2, 0x3, 0x3, 0x4, 0x3, 0x4, 0x4, 0x5, 0x3, 0x4, 0x4, 0x5, 0x4, 0x5, 0x5, 0x6,
  0x3, 0x4, 0x4, 0x5, 0x4, 0x5, 0x5, 0x6, 0x4, 0x5, 0x5, 0x6, 0x5, 0x6, 0x6, 0x7,
  0x1, 0x2, 0x2, 0x3, 0x2, 0x3, 0x3, 0x4, 0x2, 0x3, 0x3, 0x4, 0x3, 0x4, 0x4, 0x5,
  0x2, 0x3, 0x3, 0x4, 0x3, 0x4, 0x4, 0x5, 0x3, 0x4, 0x4, 0x5, 0x4, 0x5, 0x5, 0x6,
  0x2, 0x3, 0x3, 0x4, 0x3, 0x4, 0x4, 0x5, 0x3, 0x4, 0x4, 0x5, 0x4, 0x5, 0x5, 0x6,
  0x3, 0x4, 0x4, 0x5, 0x4, 0x5, 0x5, 0x6, 0x4, 0x5, 0x5, 0x6, 0x5, 0x6, 0x6, 0x7,
  0x2, 0x3, 0x3, 0x4, 0x3, 0x4, 0x4, 0x5, 0x3, 0x4, 0x4, 0x5, 0x4, 0x5, 0x5, 0x6,
  0x3, 0x4, 0x4, 0x5, 0x4, 0x5, 0x5, 0x6, 0x4, 0x5, 0x5, 0x6, 0x5, 0x6, 0x6, 0x7,
  0x3, 0x4, 0x4, 0x5, 0x4, 0x5, 0x5, 0x6, 0x4, 0x5, 0x5, 0x6, 0x5, 0x6, 0x6, 0x7,
  0x4, 0x5, 0x5, 0x6, 0x5, 0x6, 0x6, 0x7, 0x5, 0x6, 0x6, 0x7, 0x6, 0x7, 0x7, 0x8,
};

// ---------- local operators and operator-wrappers ----------
UNUSED static inline __m256
local_cmple_ps ( __m256 in1, __m256 in2 )
{
  return _mm256_cmp_ps ( in1, in2, _CMP_LE_OQ );
}

// ========== internal generic AVX2 functions ==========

// ---------- generic AVX2 operator with 2 REAL4 vector inputs to 1 UINT4 scalar and 1 UINT4 vector output (SS2uU) ----------
static inline int
XLALVectorMath_SS2uU_AVX2 ( UINT4* count, UINT4 *out, const REAL4 *in1, const REAL4 *in2, const UINT4 len, __m256 (*pred)(__m256, __m256) )
{
  *count = 0;

  // walk through vector in blocks of 8
  UINT4 i8Max = len - ( len % 8 );
  for ( UINT4 i8 = 0; i8 < i8Max; i8 += 8 )
    {
      // load vector inputs
      __m256 in8p_1 = _mm256_loadu_ps(&in1[i8]);
      __m256 in8p_2 = _mm256_loadu_ps(&in2[i8]);
      // 'pred' should set each float to 0xFFFFFFFF if satisfied, 0x0 otherwise
      __m256 res = (*pred)(in8p_1, in8p_2);
      // get bitmask from 'r', bits are set if 'pred' is satisfied
      int mask = _mm256_movemask_ps(res);
      // get correct shuffling operation so that masked items are moved before unmasked items
      int32_t shuf = shuffle[mask];
      // transform 's' into form of permutation required by _mm256_permutevar8x32_epi32()
      // - e.g. shuf = 0x75432061
      //        perm = (0x75432061 >> (28, 24, 20, 16, 12, 8, 4, 0)) & 0xf
      //             = (0x7, 0x75, 0x754, 0x7543, 0x75432, 0x754320, 0x7543206, 0x75432061) & 0xf
      //             = (0x7, 0x5, 0x4, 0x3, 0x2, 0x0, 0x6, 0x1)
      __m256i perm = _mm256_and_si256(_mm256_srlv_epi32(_mm256_set1_epi32(shuf), _mm256_set_epi32(28, 24, 20, 16, 12, 8, 4, 0)), _mm256_set1_epi32(0xf));
      // create packed int of indexes to current vector inputs
      // - i.e. idx = i8 + (7, 6, 5, 4, 3, 2, 1, 0)
      __m256i idx = _mm256_add_epi32(_mm256_set1_epi32(i8), _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0));
      // permute so that indexes of items for which 'pred' is satisfied appear first
      __m256i perm_idx = _mm256_permutevar8x32_epi32(idx, perm);
      // store indexes in output
      // - indexes of items for which 'pred' was not satisfied will be erased on the next loop
      _mm256_storeu_si256((void *) &out[*count], perm_idx);
      // increment count by number of items for which 'pred' is satisfied
      *count += popcount[mask];
    }

  // deal with the remaining (<=7) terms separately
  V8SF in8_1 = {.f={0,0,0,0,0,0,0,0}};
  V8SF in8_2 = {.f={0,0,0,0,0,0,0,0}};
  V8SF out8;
  for (UINT4 i = i8Max, j = 0; i < len; i++, j++) {
    in8_1.f[j] = in1[i];
    in8_2.f[j] = in2[i];
  }
  out8.v = (*pred)(in8_1.v, in8_2.v);
  for (UINT4 i = i8Max, j = 0; i < len; i++, j++) {
    if (out8.i[j]) {
      out[*count] = i;
      *count += 1;
    }
  }

  return XLAL_SUCCESS;

} // XLALVectorMath_SS2uU_AVX2()

// ---------- generic AVX2 operator with 1 REAL4 scalar and 1 REAL4 vector inputs to 1 UINT4 scalar and 1 UINT4 vector output (Ss2uU) ----------
static inline int
XLALVectorMath_sS2uU_AVX2 ( UINT4* count, UINT4 *out, REAL4 scalar, const REAL4 *in, const UINT4 len, __m256 (*pred)(__m256, __m256) )
{
  const V8SF scalar8 = {.f={scalar,scalar,scalar,scalar,scalar,scalar,scalar,scalar}};

  *count = 0;

  // walk through vector in blocks of 8
  UINT4 i8Max = len - ( len % 8 );
  for ( UINT4 i8 = 0; i8 < i8Max; i8 += 8 )
    {
      // load vector inputs
      __m256 in8p = _mm256_loadu_ps(&in[i8]);
      // 'pred' should set each float to 0xFFFFFFFF if satisfied, 0x0 otherwise
      __m256 res = (*pred)(scalar8.v, in8p);
      // get bitmask from 'r', bits are set if 'pred' is satisfied
      int mask = _mm256_movemask_ps(res);
      // get correct shuffling operation so that masked items are moved before unmasked items
      int32_t shuf = shuffle[mask];
      // transform 's' into form of permutation required by _mm256_permutevar8x32_epi32()
      // - e.g. shuf = 0x75432061
      //        perm = (0x75432061 >> (28, 24, 20, 16, 12, 8, 4, 0)) & 0xf
      //             = (0x7, 0x75, 0x754, 0x7543, 0x75432, 0x754320, 0x7543206, 0x75432061) & 0xf
      //             = (0x7, 0x5, 0x4, 0x3, 0x2, 0x0, 0x6, 0x1)
      __m256i perm = _mm256_and_si256(_mm256_srlv_epi32(_mm256_set1_epi32(shuf), _mm256_set_epi32(28, 24, 20, 16, 12, 8, 4, 0)), _mm256_set1_epi32(0xf));
      // create packed int of indexes to current vector inputs
      // - i.e. idx = i8 + (7, 6, 5, 4, 3, 2, 1, 0)
      __m256i idx = _mm256_add_epi32(_mm256_set1_epi32(i8), _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0));
      // permute so that indexes of items for which 'pred' is satisfied appear first
      __m256i perm_idx = _mm256_permutevar8x32_epi32(idx, perm);
      // store indexes in output
      // - indexes of items for which 'pred' was not satisfied will be erased on the next loop
      _mm256_storeu_si256((void *) &out[*count], perm_idx);
      // increment count by number of items for which 'pred' is satisfied
      *count += popcount[mask];
    }

  // deal with the remaining (<=7) terms separately
  V8SF in8 = {.f={0,0,0,0,0,0,0,0}};
  V8SF out8;
  for (UINT4 i = i8Max, j = 0; i < len; i++, j++) {
    in8.f[j] = in[i];
  }
  out8.v = (*pred)(scalar8.v, in8.v);
  for (UINT4 i = i8Max, j = 0; i < len; i++, j++) {
    if (out8.i[j]) {
      out[*count] = i;
      *count += 1;
    }
  }

  return XLAL_SUCCESS;

} // XLALVectorMath_sS2uU_AVX2()

// ========== internal AVX2 vector math functions ==========

// ---------- define vector math functions with 2 REAL4 vector inputs to 1 UINT4 scalar and 1 UINT4 vector output (SS2uU) ----------
#define DEFINE_VECTORMATH_SS2uU(NAME, AVX_PRED)                            \
  DEFINE_VECTORMATH_ANY( XLALVectorMath_SS2uU_AVX2, NAME ## REAL4, ( UINT4* count, UINT4 *out, const REAL4 *in1, const REAL4 *in2, const UINT4 len ), ( (count != NULL) && (out != NULL) && (in1 != NULL) && (in2 != NULL) ), ( count, out, in1, in2, len, AVX_PRED ) )

DEFINE_VECTORMATH_SS2uU(FindVectorLessEqual, local_cmple_ps)

// ---------- define vector math functions with 1 REAL4 scalar and 1 REAL4 vector inputs to 1 UINT4 scalar and 1 UINT4 vector output (sS2uU) ----------
#define DEFINE_VECTORMATH_sS2uU(NAME, AVX_PRED)                            \
  DEFINE_VECTORMATH_ANY( XLALVectorMath_sS2uU_AVX2, NAME ## REAL4, ( UINT4* count, UINT4 *out, REAL4 scalar, const REAL4 *in, const UINT4 len ), ( (count != NULL) && (out != NULL) && (in != NULL) ), ( count, out, scalar, in, len, AVX_PRED ) )

DEFINE_VECTORMATH_sS2uU(FindScalarLessEqual, local_cmple_ps)
