/*
*  Copyright (C) 2007 Jolien Creighton
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/* ---------- SEE LALDatatypes.dox for doxygen documentation ---------- */

#ifndef _LALATOMICDATATYPES_H
#define _LALATOMICDATATYPES_H

#include <stdint.h>
#include <lal/LALConfig.h>

/* macros for certain keywords */
#if __STDC_VERSION__ >= 199901L
# define _LAL_RESTRICT_ restrict
# define _LAL_INLINE_ inline
#elif defined __GNUC__
# define _LAL_RESTRICT_ __restrict__
# define _LAL_INLINE_ __inline__
#else
# define _LAL_RESTRICT_
# define _LAL_INLINE_
#endif

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/** \addtogroup LALDatatypes */ /*@{*/

typedef char CHAR;		/**< One-byte signed integer, see \ref LALDatatypes for more details */
typedef unsigned char UCHAR;	/**< One-byte unsigned integer, see \ref LALDatatypes for more details */
typedef unsigned char BOOLEAN;	/**< Boolean logical type, see \ref LALDatatypes for more details */

/* If INT8 etc. are already defined, undefine them */
#undef CHAR
#undef UCHAR
#undef INT2
#undef INT4
#undef INT8
#undef UINT2
#undef UINT4
#undef UINT8
#undef REAL4
#undef REAL8
#undef COMPLEX8
#undef COMPLEX16

/* Integer types */
typedef int16_t  INT2;		/**< Two-byte signed integer */
typedef int32_t  INT4;		/**< Four-byte signed integer. */
typedef int64_t  INT8;		/**< Eight-byte signed integer; on some platforms this is equivalent to <tt>long int</tt> instead. */
typedef uint16_t UINT2;		/**< Two-byte unsigned integer. */
typedef uint32_t UINT4;		/**< Four-byte unsigned integer. */
typedef uint64_t UINT8;		/**< Eight-byte unsigned integer; on some platforms this is equivalent to <tt>unsigned long int</tt> instead. */

/* Macros for integer constants */
/**
 * \def LAL_INT8_C(v) (v ## LL)
 * \brief Macro for use in defining \a v as an INT8 constant.
 *
 * This macro affixes the appropriate qualifier to form an INT8 constant.
 * For example:
 * \code
 * const INT8 jan_1_2000_gps_nanosec = LAL_INT8_C(63072001300000000)
 * \endcode
 */
#define LAL_INT8_C INT64_C

/**
 * \def LAL_UINT8_C(v) (v ## ULL)
 * \brief Macro for use in defining \a v as an UINT8 constant.
 *
 * This macro affixes the appropriate qualifier to form an UINT8 constant.
 * For example:
 * \code
 * const UINT8 jan_1_2000_gps_nanosec = LAL_UINT8_C(63072001300000000)
 * \endcode
 */
#define LAL_UINT8_C UINT64_C

/* Real types */
typedef float REAL4;    /**< Single precision real floating-point number (4 bytes). */
typedef double REAL8;   /**< Double precision real floating-point number (8 bytes). */

/* Complex types */
#ifndef SWIG /* exclude from SWIG interface */

/* Use C99 complex numbers where available: C99, gcc with non-ANSI extensions */
#if !defined(__cplusplus)
# if __STDC_VERSION__ >= 199901L || (defined(__GNUC__) && !defined(__STRICT_ANSI__))
#  define _LAL_C99_COMPLEX_NUMBERS_
# endif
#endif

#ifdef _LAL_C99_COMPLEX_NUMBERS_

#include <complex.h>

typedef float  complex COMPLEX8;	/**< Single-precision floating-point complex number (8 bytes total) */
typedef double complex COMPLEX16;	/**< Double-precision floating-point complex number (16 bytes total) */

#define crectf(re, im) (((REAL4)(re)) + _Complex_I*((REAL4)(im)))	/**< Construct a COMPLEX8 from real and imaginary parts */
#define crect(re, im)  (((REAL8)(re)) + _Complex_I*((REAL8)(im)))	/**< Construct a COMPLEX16 from real and imaginary parts */
#define cpolarf(r, th) (((REAL4)(r)) * cexpf(crectf(0, th)))		/**< Construct a COMPLEX8 from polar modulus and argument */
#define cpolar(r, th)  (((REAL8)(r)) * cexp(crect(0, th)))		/**< Construct a COMPLEX16 from polar modulus and argument */

#else /* !_LAL_C99_COMPLEX_NUMBERS_ */
/** \cond DONT_DOXYGEN */

/****************************************************************************/
/* Fall back to GSL complex number types if C99 complex numbers are not     */
/* available.  GSL complex numbers are implemented as a struct containing   */
/* an array of 2 elements of the corresponding real type.  The C99 standard */
/* should guarantee that GSL complex numbers are binary-compatible with C99 */
/* complex numbers; 6.2.5 point 13 of the standard states that C99 complex  */
/* numbers are equivalent to an array of 2 elements of the corresponding    */
/* real type, and 6.7.2.1 point 13 states that padding is never added to    */
/* the beginning of a struct.                                               */
/****************************************************************************/

#include <math.h>
#include <gsl/gsl_complex.h>

typedef gsl_complex_float COMPLEX8;
typedef gsl_complex       COMPLEX16;

static _LAL_INLINE_ COMPLEX8 crectf(const REAL4 re, const REAL4 im);
static _LAL_INLINE_ COMPLEX8 crectf(const REAL4 re, const REAL4 im) {
  COMPLEX8 z; GSL_SET_COMPLEX(&z, re, im); return z;
}
static _LAL_INLINE_ COMPLEX16 crect(const REAL8 re, const REAL8 im);
static _LAL_INLINE_ COMPLEX16 crect(const REAL8 re, const REAL8 im) {
  COMPLEX16 z; GSL_SET_COMPLEX(&z, re, im); return z;
}
static _LAL_INLINE_ COMPLEX8 cpolarf(const REAL4 r, const REAL4 th);
static _LAL_INLINE_ COMPLEX8 cpolarf(const REAL4 r, const REAL4 th) {
  COMPLEX8 z; GSL_SET_COMPLEX(&z, r*cos(th), r*sin(th)); return z;
}
static _LAL_INLINE_ COMPLEX16 cpolar(const REAL8 r, const REAL8 th);
static _LAL_INLINE_ COMPLEX16 cpolar(const REAL8 r, const REAL8 th) {
  COMPLEX16 z; GSL_SET_COMPLEX(&z, r*cos(th), r*sin(th)); return z;
}

#define crealf(z) GSL_REAL(z)
#define cimagf(z) GSL_IMAG(z)
#define creal(z)  GSL_REAL(z)
#define cimag(z)  GSL_IMAG(z)

/** \endcond */
#endif /* _LAL_C99_COMPLEX_NUMBERS_ */

#endif /* SWIG */

/*@}*/

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALATOMICDATATYPES_H */
