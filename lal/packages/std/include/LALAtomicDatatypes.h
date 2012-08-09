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

// ---------- SEE LALDatatypes.dox for doxygen documentation ----------

#ifndef _LALATOMICDATATYPES_H
#define _LALATOMICDATATYPES_H

#include <stdint.h>

#ifndef LAL_USE_OLD_COMPLEX_STRUCTS
#if defined(__cplusplus)
#include <complex>
#else
#include <complex.h>
#endif
#endif /* LAL_USE_OLD_COMPLEX_STRUCTS */

#include <lal/LALConfig.h>

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
/** \def LAL_INT8_C(v) (v ## LL)
 * \brief Macro for use in defining \a v as an INT8 constant.
 *
 * This macro affixes the appropriate qualifier to form an INT8 constant.
 * For example:
 * \code
 * const INT8 jan_1_2000_gps_nanosec = LAL_INT8_C(63072001300000000)
 * \endcode
 */
#define LAL_INT8_C INT64_C

/** \def LAL_UINT8_C(v) (v ## ULL)
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

#ifndef SWIG /* exclude from SWIG interface */

#ifndef LAL_USE_OLD_COMPLEX_STRUCTS

/* Complex types */
#if defined(__cplusplus)
typedef std::complex<float> COMPLEX8;
typedef std::complex<double> COMPLEX16;
#else
typedef float complex COMPLEX8;     /**< Single-precision floating-point complex number (8 bytes total) */
typedef double complex COMPLEX16;   /**< Double-precision floating-point complex number (16 bytes total) */
#endif

/* Complex type constructors */
#if defined(__cplusplus)
#define CX8rect( re, im) COMPLEX8(  re, im )
#define CX16rect(re, im) COMPLEX16( re, im )
#define CX8polar( r, th) ( (r) * std::exp( CX8rect(  0, th ) ) )
#define CX16polar(r, th) ( (r) * std::exp( CX16rect( 0, th ) ) )
#else
#define CX8rect( re, im) ( (re) + _Complex_I * (im) )           /**< Construct a COMPLEX8 from real and imaginary parts */
#define CX16rect(re, im) ( (re) + _Complex_I * (im) )           /**< Construct a COMPLEX16 from real and imaginary parts */
#define CX8polar( r, th) ( (r) * cexpf( CX8rect(  0, th ) ) )   /**< Construct a COMPLEX8 from polar modulus and argument */
#define CX16polar(r, th) ( (r) * cexp(  CX16rect( 0, th ) ) )   /**< Construct a COMPLEX16 from polar modulus and argument */
#endif

/* Real and imaginary part accessors */
#if defined(__cplusplus)
#define CX8re( z) std::real(z)
#define CX16re(z) std::real(z)
#define CX8im( z) std::imag(z)
#define CX16im(z) std::imag(z)
#else
#define CX8re( z) crealf(z)   /**< Get the real part of a COMPLEX8 */
#define CX16re(z) creal( z)   /**< Get the real part of a COMPLEX16 */
#define CX8im( z) cimagf(z)   /**< Get the imaginary part of a COMPLEX8 */
#define CX16im(z) cimag( z)   /**< Get the imaginary part of a COMPLEX16 */
#endif

/* Real and imaginary part assignment */
#if !defined(__cplusplus) && defined(__GNUC__)
#define setCX8re( z, re) ( __real__(z) = (re) )
#define setCX16re(z, re) ( __real__(z) = (re) )
#define setCX8im( z, im) ( __imag__(z) = (im) )
#define setCX16im(z, im) ( __imag__(z) = (im) )
#else
#define setCX8re( z, re) ( (z) = CX8rect(  re, CX8im( z) ) )   /**< Set the real part of a COMPLEX8 */
#define setCX16re(z, re) ( (z) = CX16rect( re, CX16im(z) ) )   /**< Set the real part of a COMPLEX16 */
#define setCX8im( z, im) ( (z) = CX8rect(  CX8re( z), im ) )   /**< Set the imaginary part of a COMPLEX8 */
#define setCX16im(z, im) ( (z) = CX16rect( CX16re(z), im ) )   /**< Set the imaginary part of a COMPLEX16 */
#endif

#else /* LAL_USE_OLD_COMPLEX_STRUCTS */

/** \cond DONT_DOXYGEN */
/* Old LAL complex structs, being phased out ... */
typedef struct tagCOMPLEX8 { REAL4 re; REAL4 im; } COMPLEX8;
typedef struct tagCOMPLEX16 { REAL8 re; REAL8 im; } COMPLEX16;
/** \endcond */

#endif /* LAL_USE_OLD_COMPLEX_STRUCTS */

/*@}*/

#endif /* SWIG */

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALATOMICDATATYPES_H */
