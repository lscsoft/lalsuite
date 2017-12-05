/*
 * Copyright (C) 2015 Jolien Creighton, Reinhard Prix, Karl Wette
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
 */

#ifndef LALSIMD_H
#define LALSIMD_H

#if defined(__cplusplus)
extern "C" {
#endif

/**
 * \defgroup LALSIMD_h Header LALSIMD.h
 * \ingroup lal_std
 * \author Karl Wette
 * \brief SIMD extension detection and runtime selection for LALSuite
 *
 * ### Synopsis ###
 * \code
 * #include <lal/LALSIMD.h>
 *
 * if (LAL_HAVE_SSE_RUNTIME()) {
 *   perform_sse_magic(out, in, ...);
 * }
 * \endcode
 */
/*@{*/

/**
 * SIMD instruction sets this module can detect
 */
typedef enum tagLAL_SIMD_ISET {

  LAL_SIMD_ISET_GEN,		/**< GENeric floating-point unit */
  LAL_SIMD_ISET_SSE,		/**< SSE (Streaming SIMD Extensions) */
  LAL_SIMD_ISET_SSE2,		/**< SSE version 2 */
  LAL_SIMD_ISET_SSE3,		/**< SSE version 3 */
  LAL_SIMD_ISET_SSSE3,		/**< Supplemental SSE version 3 */
  LAL_SIMD_ISET_SSE4_1,		/**< SSE version 4.1 */
  LAL_SIMD_ISET_SSE4_2,		/**< SSE version 4.2 */
  LAL_SIMD_ISET_AVX,		/**< AVX (Advanced Vector Extensions) */
  LAL_SIMD_ISET_AVX2,		/**< AVX version 2 */

  LAL_SIMD_ISET_MAX
} LAL_SIMD_ISET;

/**
 * Return true if the executing machine supports the given instruction set
 */
int XLALHaveSIMDInstructionSet(LAL_SIMD_ISET iset);

/**
 * Return the name of a given instruction set as a string
 */
const char *XLALSIMDInstructionSetName(LAL_SIMD_ISET iset);

/** \name Convenience macros for SIMD runtime selection */
/*@{*/
#define LAL_HAVE_SSE_RUNTIME()		(XLALHaveSIMDInstructionSet(LAL_SIMD_ISET_SSE))
#define LAL_HAVE_SSE2_RUNTIME()		(XLALHaveSIMDInstructionSet(LAL_SIMD_ISET_SSE2))
#define LAL_HAVE_SSE3_RUNTIME()		(XLALHaveSIMDInstructionSet(LAL_SIMD_ISET_SSE3))
#define LAL_HAVE_SSSE3_RUNTIME()	(XLALHaveSIMDInstructionSet(LAL_SIMD_ISET_SSSE3))
#define LAL_HAVE_SSE4_1_RUNTIME()	(XLALHaveSIMDInstructionSet(LAL_SIMD_ISET_SSE4_1))
#define LAL_HAVE_SSE4_2_RUNTIME()	(XLALHaveSIMDInstructionSet(LAL_SIMD_ISET_SSE4_2))
#define LAL_HAVE_AVX_RUNTIME()		(XLALHaveSIMDInstructionSet(LAL_SIMD_ISET_AVX))
#define LAL_HAVE_AVX2_RUNTIME()		(XLALHaveSIMDInstructionSet(LAL_SIMD_ISET_AVX2))
/*@}*/

/*@}*/

#if defined(__cplusplus)
}
#endif

#endif /* LALSIMD_H */
