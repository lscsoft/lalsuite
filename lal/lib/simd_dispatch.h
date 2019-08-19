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
 */

#if !defined(HAVE_SIMD_COMPILER)
#error "simd_dispatch.h must be included after config.h"
#endif

#include <lal/LALSIMD.h>
#include <lal/XLALError.h>

/*
 * Helper macros for performing runtime dispatch using function pointers
 * Sample usage:
 *
 * int SIMDFunction_AVX(...);
 * int SIMDFunction_SSE(...);
 * int SIMDFunction_GEN(...);
 *
 * int (*SIMDFunctionPtr)(...) = SIMDFunction_DISPATCH;
 *
 * int SIMDFunction_DISPATCH(...) {
 *
 *   DISPATCH_SELECT_BEGIN();
 *   DISPATCH_SELECT_AVX(SIMDFunctionPtr = SIMDFunction_AVX);
 *   DISPATCH_SELECT_SSE(SIMDFunctionPtr = SIMDFunction_SSE);
 *   DISPATCH_SELECT_END(SIMDFunctionPtr = SIMDFunction_GEN);
 *
 *   return SIMDFunction(...);
 *
 * }
 */
#define DISPATCH_SELECT_BEGIN()			do { do { } while(0)
#define DISPATCH_SELECT_END(...)		(__VA_ARGS__); } while (0)
#define DISPATCH_SELECT_NONE(...)		do { } while(0)

#if defined(HAVE_SSE_COMPILER)			/* set by config.h if compiler supports SSE */
#define DISPATCH_SELECT_SSE(...)		if (LAL_HAVE_SSE_RUNTIME()) { (__VA_ARGS__); break; } do { } while(0)
#else
#define DISPATCH_SELECT_SSE(...)		DISPATCH_SELECT_NONE()
#endif

#if defined(HAVE_SSE2_COMPILER)			/* set by config.h if compiler supports SSE2 */
#define DISPATCH_SELECT_SSE2(...)		if (LAL_HAVE_SSE2_RUNTIME()) { (__VA_ARGS__); break; } do { } while(0)
#else
#define DISPATCH_SELECT_SSE2(...)		DISPATCH_SELECT_NONE()
#endif

#if defined(HAVE_SSE3_COMPILER)			/* set by config.h if compiler supports SSE3 */
#define DISPATCH_SELECT_SSE3(...)		if (LAL_HAVE_SSE3_RUNTIME()) { (__VA_ARGS__); break; } do { } while(0)
#else
#define DISPATCH_SELECT_SSE3(...)		DISPATCH_SELECT_NONE()
#endif

#if defined(HAVE_SSSE3_COMPILER)		/* set by config.h if compiler supports SSSE3 */
#define DISPATCH_SELECT_SSSE3(...)		if (LAL_HAVE_SSSE3_RUNTIME()) { (__VA_ARGS__); break; } do { } while(0)
#else
#define DISPATCH_SELECT_SSSE3(...)		DISPATCH_SELECT_NONE()
#endif

#if defined(HAVE_SSE4_1_COMPILER)		/* set by config.h if compiler supports SSE4_1 */
#define DISPATCH_SELECT_SSE4_1(...)		if (LAL_HAVE_SSE4_1_RUNTIME()) { (__VA_ARGS__); break; } do { } while(0)
#else
#define DISPATCH_SELECT_SSE4_1(...)		DISPATCH_SELECT_NONE()
#endif

#if defined(HAVE_SSE4_2_COMPILER)		/* set by config.h if compiler supports SSE4_2 */
#define DISPATCH_SELECT_SSE4_2(...)		if (LAL_HAVE_SSE4_2_RUNTIME()) { (__VA_ARGS__); break; } do { } while(0)
#else
#define DISPATCH_SELECT_SSE4_2(...)		DISPATCH_SELECT_NONE()
#endif

#if defined(HAVE_AVX_COMPILER)			/* set by config.h if compiler supports AVX */
#define DISPATCH_SELECT_AVX(...)		if (LAL_HAVE_AVX_RUNTIME()) { (__VA_ARGS__); break; } do { } while(0)
#else
#define DISPATCH_SELECT_AVX(...)		DISPATCH_SELECT_NONE()
#endif

#if defined(HAVE_AVX2_COMPILER)			/* set by config.h if compiler supports AVX2 */
#define DISPATCH_SELECT_AVX2(...)		if (LAL_HAVE_AVX2_RUNTIME()) { (__VA_ARGS__); break; } do { } while(0)
#else
#define DISPATCH_SELECT_AVX2(...)		DISPATCH_SELECT_NONE()
#endif
