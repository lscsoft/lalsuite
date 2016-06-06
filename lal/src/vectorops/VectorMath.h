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
 * \author Reinhard Prix, Karl Wette
 *
 * \brief Functions for performing fast math on vectors of numbers, using SIMD (SSE, AVX, ...) instructions if available.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/VectorMath.h>
 * \endcode
 *
 * ### Alignment ###
 *
 * Neither input nor output vectors are \b required to have any particular memory alignment. Nevertheless, performance
 * \e may be improved if vectors are 16-byte aligned for SSE, and 32-byte aligned for AVX.
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

/** Compute \f$\text{out} = \sin(\text{in})\f$ over REAL4 vectors \c out, \c in with \c len elements */
int XLALVectorSinREAL4 ( REAL4 *out, const REAL4 *in, const UINT4 len );

/** Compute \f$\text{out} = \cos(\text{in})\f$ over REAL4 vectors \c out, \c in with \c len elements */
int XLALVectorCosREAL4 ( REAL4 *out, const REAL4 *in, const UINT4 len );

/** Compute \f$\text{out} = \exp(\text{in})\f$ over REAL4 vectors \c out, \c in with \c len elements */
int XLALVectorExpREAL4 ( REAL4 *out, const REAL4 *in, const UINT4 len );

/** Compute \f$\text{out} = \log(\text{in})\f$ over REAL4 vectors \c out, \c in with \c len elements */
int XLALVectorLogREAL4 ( REAL4 *out, const REAL4 *in, const UINT4 len );

/** Compute \f$\text{out1} = \sin(\text{in}), \text{out2} = \cos(\text{in})\f$ over REAL4 vectors \c out1, \c out2, \c in with \c len elements */
int XLALVectorSinCosREAL4 ( REAL4 *out1, REAL4 *out2, const REAL4 *in, const UINT4 len );

/** Compute \f$\text{out1} = \sin(2\pi \text{in}), \text{out2} = \cos(2\pi \text{in})\f$ over REAL4 vectors \c out1, \c out2, \c in with \c len elements */
int XLALVectorSinCos2PiREAL4 ( REAL4 *out1, REAL4 *out2, const REAL4 *in, const UINT4 len );

/* -------------------- exported vector by vector operations -------------------- */

/** Compute \f$\text{out} = \text{in1} + \text{in2}\f$ over REAL4 vectors \c in1 and \c in2 with \c len elements */
int XLALVectorAddREAL4 ( REAL4 *out, const REAL4 *in1, const REAL4 *in2, const UINT4 len);

/** Compute \f$\text{out} = \text{in1} \times \text{in2}\f$ over REAL4 vectors \c in1 and \c in2 with \c len elements */
int XLALVectorMultiplyREAL4 ( REAL4 *out, const REAL4 *in1, const REAL4 *in2, const UINT4 len);


/* -------------------- exported vector by scalar operations -------------------- */

/** Compute \f$\text{out} = \text{scalar} + \text{in}\f$ over REAL4 vector \c in with \c len elements */
int XLALVectorShiftREAL4 ( REAL4 *out, REAL4 scalar, const REAL4 *in, const UINT4 len);

/** Compute \f$\text{out} = \text{scalar} \times \text{in}\f$ over REAL4 vector \c in with \c len elements */
int XLALVectorScaleREAL4 ( REAL4 *out, REAL4 scalar, const REAL4 *in, const UINT4 len);


/** @} */

#ifdef  __cplusplus
}
#endif

#endif /* _VECTORMATH_H */
