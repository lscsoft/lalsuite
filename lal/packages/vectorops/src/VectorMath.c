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
#include <stdbool.h>
#include <math.h>

#include <config.h>

#include <lal/LALConstants.h>

#include <lal/VectorMath.h>

// ---------- local DEFINES ----------
// ----- Macros -----
// ---------- internal types ----------
// ---------- Global variables ----------
//---------- internal prototypes ----------

/*==================== FUNCTION DEFINITIONS ====================*/

/// Vector version of sinf(x) to take advantage of parallel-math devices (SSE,AVX,GPU...)
int
XLALVectorSinf ( REAL4Vector *out, const REAL4Vector *in )
{
#ifdef HAVE_AVX
  return XLALVectorSinf_AVX ( out, in );
#elif defined HAVE_SSE
  return XLALVectorSinf_SSE ( out, in );
#else // vanilla LALsuite
  return XLALVectorSinf_FPU ( out, in );
#endif
   return XLAL_SUCCESS;
} // XLALVectorSinf()

/// Vector version of cosf(x) to take advantage of parallel-math devices (SSE,AVX,GPU...)
int
XLALVectorCosf ( REAL4Vector *out, const REAL4Vector *in )
{
#ifdef HAVE_AVX
  return XLALVectorCosf_AVX ( out, in );
#elif defined HAVE_SSE
  return XLALVectorCosf_SSE ( out, in );
#else // vanilla LALsuite
  return XLALVectorCosf_FPU ( out, in );
#endif
   return XLAL_SUCCESS;
} // XLALVectorCosf()

/// Vector version of expf(x) to take advantage of parallel-math devices (SSE,AVX,GPU...)
int
XLALVectorExpf ( REAL4Vector *out, const REAL4Vector *in )
{
#ifdef HAVE_AVX
  return XLALVectorExpf_AVX ( out, in );
#elif defined HAVE_SSE
  return XLALVectorExpf_SSE ( out, in );
#else // vanilla LALsuite
  return XLALVectorExpf_FPU ( out, in );
#endif
   return XLAL_SUCCESS;
} // XLALVectorExpf()

/// Vector version of logf(x) to take advantage of parallel-math devices (SSE,AVX,GPU...)
int
XLALVectorLogf ( REAL4Vector *out, const REAL4Vector *in )
{
#ifdef HAVE_AVX
  return XLALVectorLogf_AVX ( out, in );
#elif defined HAVE_SSE
  return XLALVectorLogf_SSE ( out, in );
#else // vanilla LALsuite
  return XLALVectorLogf_FPU ( out, in );
#endif
   return XLAL_SUCCESS;
} // XLALVectorLogf()

/// Vector version of sincosf(x) to take advantage of parallel-math devices (SSE,AVX,GPU...)
int
XLALVectorSinCosf ( REAL4Vector *sinx, REAL4Vector *cosx, const REAL4Vector *x )
{
#ifdef HAVE_AVX
  return XLALVectorSinCosf_AVX ( sinx, cosx, x );
#elif defined HAVE_SSE
  return XLALVectorSinCosf_SSE ( sinx, cosx, x );
#else // vanilla LALsuite
  return XLALVectorSinCosf_FPU ( sinx, cosx, x );
#endif
   return XLAL_SUCCESS;
} // XLALVectorSinCosf()


// -------------------- our own failsafe aligned memory handling --------------------
///
/// Create a special REAL4 Vector with 32-byte aligned memory 'data' array.
///
/// This does not rely on posix_memalign() being available, and should compile+run everywhere.
/// Use XLALDestroyREAL4VectorAligned32() to free.
///
REAL4VectorAligned32 *
XLALCreateREAL4VectorAligned32 ( UINT4 length )
{
  const size_t align = 32;
  REAL4VectorAligned32 *ret;
  XLAL_CHECK_NULL ( (ret = XLALMalloc ( sizeof(*ret))) != NULL, XLAL_ENOMEM );

  ret->length = length;
  UINT4 paddedLength = length + align - 1;
  XLAL_CHECK_NULL ( (ret->data0 = XLALMalloc ( paddedLength * sizeof(REAL4) )) != NULL, XLAL_ENOMEM );

  size_t remBytes = ((size_t)ret->data0) % align;
  size_t offsetBytes = (align - remBytes) % align;
  ret->data = (void*)(((char*)ret->data0) + offsetBytes);

  XLAL_CHECK_NULL ( isMemAligned(ret->data,align), XLAL_EFAULT, "Failed to allocate %zd-byte aligned memory. Must be a coding error.\n", align );

  return ret;
} // XLALCreateREAL4VectorAligned32()

///
/// Destroy special 32-byte aligned  REAL4VectorAligned32
///
void
XLALDestroyREAL4VectorAligned32 ( REAL4VectorAligned32 *in )
{
  if ( !in ) {
    return;
  }
  if ( in->data0 ) {
    XLALFree ( in->data0 );
  }

  XLALFree ( in );

  return;

} // XLALDestroyREAL4VectorAligned32()
