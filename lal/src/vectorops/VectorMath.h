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
 * \brief Functions for performing fast math on vectors of numbers, using SIMD (SSE, AVX, ...) instructions if available.
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
/** enumerate all supported vector 'devices' (SSE, AVX, ...) that can be chosen from at runtime
 */
typedef enum tagVectorDevice_type
{
  /* \cond DONT_DOXYGEN */
  VECTORDEVICE_START,	/**< start marker for range checking */
  /* \endcond */
  VECTORDEVICE_FPU,	/**< not really a vector device, but added here as failsafe fallback option */
  VECTORDEVICE_SSE,	/**< SSE extension */
  VECTORDEVICE_AVX,	/**< AVX extension */
  /* \cond DONT_DOXYGEN */
  VECTORDEVICE_END	/**< end marker for range checking */
  /* \endcond */
}
VectorDevice_type;

/** We provide our own local aligned-memory handling functions, until this
 * may later be merged upstream into the LALMalloc module
 */
typedef struct tagREAL4VectorAligned
{
  UINT4 length;		/**< number of 'usable' array entries (fully aligned) */
  REAL4 *data;		/**< start of aligned memory block */
  REAL4 *data0;		/**< actual physical start of memory block, possibly not aligned */
} REAL4VectorAligned;

/* ---------- Prototypes ---------- */

/* ----- runtime handling of vector device switching */
int XLALVectorDeviceSet ( VectorDevice_type device );
VectorDevice_type XLALVectorDeviceGet ( void );
int XLALVectorDeviceIsAvailable ( VectorDevice_type device );
const CHAR *XLALVectorDeviceName ( VectorDevice_type device );
CHAR *XLALVectorDeviceHelpString ( void );
int XLALVectorDeviceParseString ( VectorDevice_type *device, const char *s );

/* ----- aligned-memory handling */
REAL4VectorAligned *XLALCreateREAL4VectorAligned ( UINT4 length, UINT4 align );
void XLALDestroyREAL4VectorAligned ( REAL4VectorAligned *in );

/* ----- exported vector math functions */
  int XLALVectorSinf     ( REAL4 *out, const REAL4 *in, UINT4 length );
int XLALVectorCosf     ( REAL4 *out, const REAL4 *in, UINT4 length );
int XLALVectorExpf     ( REAL4 *out, const REAL4 *in, UINT4 length );
int XLALVectorLogf     ( REAL4 *out, const REAL4 *in, UINT4 length );
int XLALVectorSinCosf  ( REAL4 *sinx,    REAL4 *cosx,    const REAL4 *x, UINT4 length );
int XLALVectorSinCosf2PI(REAL4 *sin2pix, REAL4 *cos2pix, const REAL4 *x, UINT4 length );

/* ---------- module internal prototypes ---------- */
/* these should not be used except within this module, that's why the
 * exported API only declared the device-generic XLALVector<Funcf>() functions
 */
#ifdef IN_VECTORMATH
int XLALVectorSinf_FPU ( REAL4 *out, const REAL4 *in, UINT4 length );
int XLALVectorSinf_SSE ( REAL4 *out, const REAL4 *in, UINT4 length );
int XLALVectorSinf_AVX ( REAL4 *out, const REAL4 *in, UINT4 length );

int XLALVectorCosf_FPU ( REAL4 *out, const REAL4 *in, UINT4 length );
int XLALVectorCosf_SSE ( REAL4 *out, const REAL4 *in, UINT4 length );
int XLALVectorCosf_AVX ( REAL4 *out, const REAL4 *in, UINT4 length );

int XLALVectorExpf_FPU ( REAL4 *out, const REAL4 *in, UINT4 length );
int XLALVectorExpf_SSE ( REAL4 *out, const REAL4 *in, UINT4 length );
int XLALVectorExpf_AVX ( REAL4 *out, const REAL4 *in, UINT4 length );

int XLALVectorLogf_FPU ( REAL4 *out, const REAL4 *in, UINT4 length );
int XLALVectorLogf_SSE ( REAL4 *out, const REAL4 *in, UINT4 length );
int XLALVectorLogf_AVX ( REAL4 *out, const REAL4 *in, UINT4 length );

int XLALVectorSinCosf_FPU ( REAL4 *sinx, REAL4 *cosx, const REAL4 *x, UINT4 length );
int XLALVectorSinCosf_SSE ( REAL4 *sinx, REAL4 *cosx, const REAL4 *x, UINT4 length );
int XLALVectorSinCosf_AVX ( REAL4 *sinx, REAL4 *cosx, const REAL4 *x, UINT4 length );

int XLALVectorSinCosf2PI_FPU ( REAL4 *sin2pix, REAL4 *cos2pix, const REAL4 *x, UINT4 length );
int XLALVectorSinCosf2PI_SSE ( REAL4 *sin2pix, REAL4 *cos2pix, const REAL4 *x, UINT4 length );
int XLALVectorSinCosf2PI_AVX ( REAL4 *sin2pix, REAL4 *cos2pix, const REAL4 *x, UINT4 length );
#endif

/** @} */

#ifdef  __cplusplus
}
#endif

#endif /* _VECTORMATH_H */