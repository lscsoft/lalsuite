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

#define SIMD_INSTRSET FPU

// ---------- local math functions ----------

static void local_sincosf(float in, float *out1, float *out2) {
  *out1 = sinf ( in );
  *out2 = cosf ( in );
}

static void local_sincosf_2pi(float in, float *out1, float *out2) {
  *out1 = sinf ( (REAL4)LAL_TWOPI * in );
  *out2 = cosf ( (REAL4)LAL_TWOPI * in );
}

// ========== internal generic FPU functions ==========

static inline int
XLALVectorFuncf1T1_FPU ( REAL4 *out, const REAL4 *in, const UINT4 len, float (*f)(float) )
{
  for ( UINT4 i = 0; i < len; i ++ )
    {
      out[i] = (*f) ( in[i] );
    } // for i < len
  return XLAL_SUCCESS;
} // XLALVectorFuncf1T1_FPU()

static inline int
XLALVectorFuncf1T2_FPU ( REAL4 *out1, REAL4 *out2, const REAL4 *in, const UINT4 len, void (*f)(float, float*, float*) )
{
  for ( UINT4 i = 0; i < len; i ++ )
    {
      (*f) ( in[i], &(out1[i]), &(out2[i]) );
    } // for i < len
  return XLAL_SUCCESS;
} // XLALVectorFuncf1T1_FPU()

// ========== internal FPU vector math functions ==========

#define DEFINE_VECTORMATH_FUNCF_1T1(NAME, SSE_FUNC) \
  DEFINE_VECTORMATH_ANY( XLALVectorFuncf1T1_FPU, NAME, ( REAL4 *out, const REAL4 *in, const UINT4 len ), ( (out != NULL) && (in != NULL) ), ( out, in, len, SSE_FUNC ) )

DEFINE_VECTORMATH_FUNCF_1T1(Sinf, sinf)
DEFINE_VECTORMATH_FUNCF_1T1(Cosf, cosf)
DEFINE_VECTORMATH_FUNCF_1T1(Expf, expf)
DEFINE_VECTORMATH_FUNCF_1T1(Logf, logf)

#define DEFINE_VECTORMATH_FUNCF_1T2(NAME, SSE_FUNC) \
  DEFINE_VECTORMATH_ANY( XLALVectorFuncf1T2_FPU, NAME, ( REAL4 *out1, REAL4 *out2, const REAL4 *in, const UINT4 len ), ( (out1 != NULL) && (out2 != NULL) && (in != NULL) ), ( out1, out2, in, len, SSE_FUNC ) )

DEFINE_VECTORMATH_FUNCF_1T2(SinCosf, local_sincosf)
DEFINE_VECTORMATH_FUNCF_1T2(SinCosf2PI, local_sincosf_2pi)
