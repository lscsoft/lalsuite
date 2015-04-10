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
 * \brief Functions for performing fast math on vectors of numbers, using SIMD (SSE, AVX, ...) instructions if available.
 *
 * ### Synopsis ###
 * \code
 * #include <lal/VectorMath.h>
 * \endcode
 */
/** @{ */

/* -------------------- our own failsafe aligned memory handling -------------------- */

/** A special REAL4 Vector with n-byte aligned memory \c data array */
typedef struct tagREAL4VectorAligned
{
  UINT4 length;		/**< number of 'usable' array entries (fully aligned) */
  REAL4 *data;		/**< start of aligned memory block */
  REAL4 *data0;		/**< actual physical start of memory block, possibly not aligned */
} REAL4VectorAligned;

REAL4VectorAligned *XLALCreateREAL4VectorAligned ( const UINT4 length, const UINT4 align );
void XLALDestroyREAL4VectorAligned ( REAL4VectorAligned *in );

/* -------------------- exported vector math functions -------------------- */

/** Compute \f$y = \sin(in)\f$ over REAL4 vectors \c out, \c in with \c len elements */
int XLALVectorSinREAL4 ( REAL4 *out, const REAL4 *in, const UINT4 len );

/** Compute \f$y = \cos(in)\f$ over REAL4 vectors \c out, \c in with \c len elements */
int XLALVectorCosREAL4 ( REAL4 *out, const REAL4 *in, const UINT4 len );

/** Compute \f$y = \exp(in)\f$ over REAL4 vectors \c out, \c in with \c len elements */
int XLALVectorExpREAL4 ( REAL4 *out, const REAL4 *in, const UINT4 len );

/** Compute \f$y = \log(in)\f$ over REAL4 vectors \c out, \c in with \c len elements */
int XLALVectorLogREAL4 ( REAL4 *out, const REAL4 *in, const UINT4 len );

/** Compute \f$y_1 = \sin(in), out_2 = \cos(in)\f$ over REAL4 vectors \c out1, \c out2, \c in with \c len elements */
int XLALVectorSinCosREAL4 ( REAL4 *out1, REAL4 *out2, const REAL4 *in, const UINT4 len );

/** Compute \f$y_1 = \sin(2\pi in), out_2 = \cos(2\pi in)\f$ over REAL4 vectors \c out1, \c out2, \c in with \c len elements */
int XLALVectorSinCos2PiREAL4 ( REAL4 *out1, REAL4 *out2, const REAL4 *in, const UINT4 len );

/** @} */

#ifdef  __cplusplus
}
#endif

#endif /* _VECTORMATH_H */
