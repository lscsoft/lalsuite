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

/*
 * SIMD extension detection and runtime selection for LALSuite
 *
 * SIMD detection based on the following files from:
 *   http://www.agner.org/optimize/#vectorclass
 ****************************  instrset.h   **********************************
 * Author:        Agner Fog
 * Date created:  2012-05-30
 * Last modified: 2014-10-22
 * Version:       1.16
 * Project:       vector classes
 *
 * (c) Copyright 2012 - 2014 GNU General Public License www.gnu.org/licenses
 **************************  instrset_detect.cpp   ****************************
 * Author:        Agner Fog
 * Date created:  2012-05-30
 * Last modified: 2014-07-23
 * Version:       1.14
 * Project:       vector classes
 *
 * (c) Copyright 2012 - 2014 GNU General Public License http://www.gnu.org/licenses
 ******************************************************************************
 */

#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include <config.h>

#include <lal/LALSIMD.h>
#include <lal/LALConfig.h>
#include <lal/LALError.h>
#include <lal/XLALError.h>
#include <lal/LALString.h>

/* Check that this file is being compiled for x86 */
#if defined(__x86_64__) || defined(_M_X64)
#define HAVE_X86   /* x86 64-bit */
#elif defined(__i386) || defined(_M_IX86)
#define HAVE_X86   /* x86 32-bit */
#endif

#if defined(HAVE_X86) && ( defined(__GNUC__) || defined(__clang__) ) && defined(HAVE_CPUID_H)
#include <cpuid.h>
#define HAVE__GET_CPUID 1
#endif

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/* selected SIMD instruction set */
static LAL_SIMD_ISET selected_iset;

/* array of instruction set names */
static const char *const iset_names[LAL_SIMD_ISET_MAX] = {
  [LAL_SIMD_ISET_GEN]		= "GEN",
  [LAL_SIMD_ISET_SSE]		= "SSE",
  [LAL_SIMD_ISET_SSE2]		= "SSE2",
  [LAL_SIMD_ISET_SSE3]		= "SSE3",
  [LAL_SIMD_ISET_SSSE3]		= "SSSE3",
  [LAL_SIMD_ISET_SSE4_1]	= "SSE4.1",
  [LAL_SIMD_ISET_SSE4_2]	= "SSE4.2",
  [LAL_SIMD_ISET_AVX]		= "AVX",
  [LAL_SIMD_ISET_AVX2]		= "AVX2",
};

/* pthread locking to make SIMD detection thread-safe */
#ifdef LAL_PTHREAD_LOCK
#include <pthread.h>
static pthread_once_t lalOnce = PTHREAD_ONCE_INIT;
#define LAL_ONCE(init) pthread_once(&lalOnce, (init))
#else
static int lalOnce = 1;
#define LAL_ONCE(init) (lalOnce ? (init)(), lalOnce = 0 : 0)
#endif

/*
 * Define interface to 'cpuid' instruction
 * input:  eax = functionnumber, ecx = 0
 * output: eax = output[0], ebx = output[1], ecx = output[2], edx = output[3]
 */
static inline UNUSED void cpuid( uint32_t output[4], UNUSED int functionnumber ) {

#if defined(HAVE_X86)

#if HAVE__GET_CPUID

  __get_cpuid(functionnumber, &output[0], &output[1], &output[2], &output[3]);

#elif defined(__GNUC__) || defined(__clang__)	// weird case: gcc|clang but NO cpuid.h file, can happen on Macs for old gcc's: give up here

  output[0] = output[1] = output[2] = output[3] = 0;

#else

  /* Use MASM/Intel inline assembly */
  __asm__ {
    mov eax, functionnumber
    xor ecx, ecx
    cpuid
    mov esi, output
    mov [esi],    eax
    mov [esi+4],  ebx
    mov [esi+8],  ecx
    mov [esi+12], edx
  }

#endif

#else    /* for non-X86 platforms */

  output[0] = output[1] = output[2] = output[3] = 0;

#endif

  return;

} // cpuid()

/*
 * Define interface to 'xgetbv' instruction
 */
static inline UNUSED int64_t xgetbv( UNUSED int ctr ) {

#if defined(HAVE_X86)

#if defined(__GNUC__) || defined(__clang__)

  /* Use GNU/AT&T inline assembly */
  uint32_t a, d;
  __asm__(".byte 0x0f,0x01,0xd0" \
          : "=a"(a),"=d"(d) \
          : "c"(ctr)
          );
  return a | (((uint64_t) d) << 32);

#else

  /* Use MASM/Intel inline assembly */
  uint32_t a, d;
  __asm__ {
    mov ecx, ctr
    _emit 0x0f
    _emit 0x01
    _emit 0xd0
    mov a, eax
    mov d, edx
  }
  return a | (((uint64_t) d) << 32);

#endif   /* inline assembly */

#else    /* !HAVE_X86 */

  return 0;

#endif   /* HAVE_X86 */

}

/*
 * Detect instruction set
 */
static LAL_SIMD_ISET detect_instruction_set(void) {

  /* cpuid results */
  uint32_t abcd[4] = {0, 0, 0, 0};

  LAL_SIMD_ISET iset = LAL_SIMD_ISET_GEN;

  cpuid(abcd, 0);					/* call cpuid function 0 */
  if (abcd[0] == 0) return iset;			/* no further cpuid function supported */
  cpuid(abcd, 1);					/* call cpuid function 1 for feature flags */
  if ((abcd[3] & (1 <<  0)) == 0) return iset;		/* no floating point */
  if ((abcd[3] & (1 << 23)) == 0) return iset;		/* no MMX */
  if ((abcd[3] & (1 << 15)) == 0) return iset;		/* no conditional move */
  if ((abcd[3] & (1 << 24)) == 0) return iset;		/* no FXSAVE */
  if ((abcd[3] & (1 << 25)) == 0) return iset;		/* no SSE */
  iset = LAL_SIMD_ISET_SSE;				/* SSE detected */

  if ((abcd[3] & (1 << 26)) == 0) return iset;		/* no SSE2 */
  iset = LAL_SIMD_ISET_SSE2;				/* SSE2 detected */

  if ((abcd[2] & (1 <<  0)) == 0) return iset;		/* no SSE3 */
  iset = LAL_SIMD_ISET_SSE3;				/* SSE3 detected */

  if ((abcd[2] & (1 <<  9)) == 0) return iset;		/* no SSSE3 */
  iset = LAL_SIMD_ISET_SSSE3;				/* SSSE3 detected */

  if ((abcd[2] & (1 << 19)) == 0) return iset;		/* no SSE4.1 */
  iset = LAL_SIMD_ISET_SSE4_1;				/* SSE4.1 detected */

  if ((abcd[2] & (1 << 23)) == 0) return iset;		/* no POPCNT */
  if ((abcd[2] & (1 << 20)) == 0) return iset;		/* no SSE4.2 */
  iset = LAL_SIMD_ISET_SSE4_2;				/* SSE4.2 detected */

  if ((abcd[2] & (1 << 27)) == 0) return iset;		/* XSAVE not enabled in O.S. */
  if ((xgetbv(0) & 6) != 6) return iset;		/* AVX not enabled in O.S. */
  if ((abcd[2] & (1 << 28)) == 0) return iset;		/* no AVX */
  iset = LAL_SIMD_ISET_AVX;				/* AVX detected */

  cpuid(abcd, 7);					/* call cpuid function 7 for feature flags */
  if ((abcd[1] & (1 <<  5)) == 0) return iset;		/* no AVX2 */
  iset = LAL_SIMD_ISET_AVX2;				/* AVX2 detected */

  return iset;

}

/*
 * Select instruction set, allowing guru users to down-select
 */
static void select_instruction_set(void) {

  /* Detect instruction set */
  selected_iset = detect_instruction_set();
  if (selected_iset >= LAL_SIMD_ISET_MAX) {
    lalAbortHook("%s: SIMD instruction set detection failed!!\n", __func__);
    return;
  }

  /* Check if user wants to down-select instruction set */
  const char *env = getenv("LAL_SIMD_ISET");
  if (env == NULL || *env == '\0') {
    return;
  }

  /* Try to match LAL_SIMD_ISET to an instruction set name */
  LAL_SIMD_ISET user_iset = LAL_SIMD_ISET_MAX;
  for (LAL_SIMD_ISET i = 0; i < LAL_SIMD_ISET_MAX; ++i) {
    if (XLALStringCaseCompare(env, iset_names[i]) == 0) {
      user_iset = i;
      break;
    }
  }
  if (user_iset == LAL_SIMD_ISET_MAX) {
    lalAbortHook("%s: LAL_SIMD_ISET='%s' does not match a SIMD instruction set\n", __func__, env);
    return;
  }

  /* Check user is not trying to select an unavailable instruction set */
  if (user_iset > selected_iset) {
    lalAbortHook("%s: LAL_SIMD_ISET='%s' is not available on this machine\n", __func__, env);
    return;
  }

  /* select user-requested instruction set */
  selected_iset = user_iset;

  return;

}

int XLALHaveSIMDInstructionSet(LAL_SIMD_ISET iset) {
  LAL_ONCE(select_instruction_set);
  return (iset < LAL_SIMD_ISET_MAX) && (iset <= selected_iset);
}

const char *XLALSIMDInstructionSetName(LAL_SIMD_ISET iset) {
  XLAL_CHECK_NULL(iset < LAL_SIMD_ISET_MAX, XLAL_EINVAL);
  return iset_names[iset];
}
