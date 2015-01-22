/*
 * Copyright (C) 2015 Reinhard Prix
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

#ifndef _VECTORMATH_H
#define _VECTORMATH_H

#include <lal/LALStdlib.h>

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \defgroup VectorMath_h Header VectorMath.h
 * \ingroup lal_vectorops
 * \author Reinhard Prix
 *
 * \brief Various functions for performing fast math on vectors of numbers, using SIMD instructions if available.
 *
 * ### Synopsis ###
 * \code
 * #include <lal/VectorMath.h>
 * \endcode
 */
/** @{ */

/* ---------- Macros ---------- */
#define isMemAligned(x,align)  (((size_t)(x) % (align)) == 0)

/* ---------- exported Types ---------- */
/* we provide our own local aligned-memory handling functions, until this
 * may later be merged upstream into the LALMalloc module
 */
typedef struct tagREAL4VectorAligned32
{
  UINT4 length;		/**< number of 'usable' array entries (fully aligned) */
  REAL4 *data;		/**< start of 32-byte aligned memory block */
  REAL4 *data0;		/**< actual physical start of memory block, possibly not aligned */
} REAL4VectorAligned32;

/* ---------- Prototypes ---------- */
int XLALVectorSinf     ( REAL4Vector *out, const REAL4Vector *in );
int XLALVectorSinf_FPU ( REAL4Vector *out, const REAL4Vector *in );
int XLALVectorSinf_SSE ( REAL4Vector *out, const REAL4Vector *in );
int XLALVectorSinf_AVX ( REAL4Vector *out, const REAL4Vector *in );

int XLALVectorCosf     ( REAL4Vector *out, const REAL4Vector *in );
int XLALVectorCosf_FPU ( REAL4Vector *out, const REAL4Vector *in );
int XLALVectorCosf_SSE ( REAL4Vector *out, const REAL4Vector *in );
int XLALVectorCosf_AVX ( REAL4Vector *out, const REAL4Vector *in );

int XLALVectorExpf     ( REAL4Vector *out, const REAL4Vector *in );
int XLALVectorExpf_FPU ( REAL4Vector *out, const REAL4Vector *in );
int XLALVectorExpf_SSE ( REAL4Vector *out, const REAL4Vector *in );
int XLALVectorExpf_AVX ( REAL4Vector *out, const REAL4Vector *in );

int XLALVectorLogf     ( REAL4Vector *out, const REAL4Vector *in );
int XLALVectorLogf_FPU ( REAL4Vector *out, const REAL4Vector *in );
int XLALVectorLogf_SSE ( REAL4Vector *out, const REAL4Vector *in );
int XLALVectorLogf_AVX ( REAL4Vector *out, const REAL4Vector *in );

int XLALVectorSinCosf     ( REAL4Vector *sinx, REAL4Vector *cosx, const REAL4Vector *x );
int XLALVectorSinCosf_FPU ( REAL4Vector *sinx, REAL4Vector *cosx, const REAL4Vector *x );
int XLALVectorSinCosf_SSE ( REAL4Vector *sinx, REAL4Vector *cosx, const REAL4Vector *x );
int XLALVectorSinCosf_AVX ( REAL4Vector *sinx, REAL4Vector *cosx, const REAL4Vector *x );

int XLALVectorSinCosf2PI     ( REAL4Vector *sin2pix, REAL4Vector *cos2pix, const REAL4Vector *x );
int XLALVectorSinCosf2PI_FPU ( REAL4Vector *sin2pix, REAL4Vector *cos2pix, const REAL4Vector *x );
int XLALVectorSinCosf2PI_SSE ( REAL4Vector *sin2pix, REAL4Vector *cos2pix, const REAL4Vector *x );
int XLALVectorSinCosf2PI_AVX ( REAL4Vector *sin2pix, REAL4Vector *cos2pix, const REAL4Vector *x );

REAL4VectorAligned32 *XLALCreateREAL4VectorAligned32 ( UINT4 length );
void XLALDestroyREAL4VectorAligned32 ( REAL4VectorAligned32 *in );
/* @} */

#ifdef  __cplusplus
}
#endif

#endif /* _VECTORMATH_H */
