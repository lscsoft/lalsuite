/** \file
 * \ingroup std
 * \author Creighton, J. D. E., and Creighton, T. D.
 * \date $Id$
 * \brief The primative LAL datatypes.
 *
 * This header defines the primative LAL datatypes.  These datatypes 
 * are: CHAR, INT2, INT4, INT8 (signed integer types); UCHAR, UINT2
 * UINT4, UINT8 (unsigned integer types); REAL4, REAL8 (single- and
 * double-precision floating point types); and COMPLEX8 and COMPLEX16
 * (single- and double-precision floating point complex types).
 * Note that the complex types are implemented as structures.  This is
 * because LAL conforms to the C89 standard rather than the C99
 * standard.  The non-complex datatypes are known as <em>atomic</em>
 * datatypes: these may be passed directly to functions or used as
 * the return value of XLAL functions.
 *
 * The number in the name of each of these types (other than
 * CHAR and UCHAR) is the number of 8-bit bytes that the datatype
 * occupies.  E.g., INT4 is a four-byte integer.  In C99 it would be
 * called <tt>int32_t</tt>. While the size of types such as <tt>int</tt> and
 * <tt>long int</tt> are platform dependent, the sizes of INT4 and INT8
 * is always 4- and 8-bytes respectively (platform independent).
 * REAL4 and REAL8 are always type <tt>float</tt> and <tt>double</tt>
 * but LAL is only supported on platforms in which these types conform
 * to the IEEE 754 standard.
 *
 * This documentation shows how the datatypes are defined on today's
 * most common (32-bit) platform.  Be careful in particular about
 * the integer type INT8.  On 64-bit platforms it will usually be
 * equivalent to type <tt>long int</tt> rather than type
 * <tt>long long int</tt>.
 */

/*************************** <lalVerbatim file="LALAtomicDatatypesHV">
Author: J. D. E. Creighton, T. D. Creighton
$Id$
**************************** </lalVerbatim> */

/* <lalLaTeX>

\subsection{Primitive datatypes}
\label{ss:LALAtomicDatatypes.h}

The primitive datatypes are defined in a separate header
\verb@LALAtomicDatatypes.h@, which is included by
\verb@LALDatatypes.h@.  This is done in order to facilitate the
interface between LAL and non-LAL modules.  By including just
\verb@LALAtomicDatatypes.h@, a non-LAL module can ensure that it is
using the same arithmetic standard as LAL, without being burdened by
LAL's more specialized structures.

Primitive datatypes are those that conceptually store a single number
or quantity.  They include both the \emph{atomic} datatypes and the
complex datatypes.

\subsubsection*{Atomic datatypes}

Atomic LAL datatypes are platform-independent datatypes corresponding
to the basic types in the C/C++ language.  However, since the C/C++
types are not necessarily the same across platforms, the actual
mapping between LAL and C/C++ datatypes may be different on different
platforms.  The following table lists the LAL atomic datatypes, their
size and range, and the C/C++ datatype to which they \emph{usually}
correspond.

\begin{center}
\begin{tabular}{|lccc|}
\hline
  Type  & Bytes & Range                        & Usual C/C++ type  \\
\hline
\tt CHAR    & 1 & \verb@'\0'@ to \verb@'\255'@ & \tt char          \\
\tt UCHAR   & 1 & \verb@'\0'@ to \verb@'\255'@ & \tt unsigned char \\
\tt BOOLEAN & 1 & 0 or positive                & \tt unsigned char \\
\hline
\tt INT2    & 2 & $-2^{-15}$ to $2^{15}-1$     & \tt short              \\
\tt INT4    & 4 & $-2^{-31}$ to $2^{31}-1$     & {\tt int} or
						{\tt long}              \\
\tt INT8    & 8 & $-2^{-63}$ to $2^{63}-1$     & \tt long long          \\
\tt UINT2   & 2 &          0 to $2^{16}-1$     & \tt unsigned short     \\
\tt UINT4   & 4 &          0 to $2^{32}-1$     & {\tt unsigned int} or
						{\tt long}              \\
\tt UINT8   & 8 &          0 to $2^{64}-1$     & \tt unsigned long long \\
\hline
\tt REAL4   & 4 &  $-3.4\times10^{38}$ to $3.4\times10^{38}$  & \tt float  \\
\tt REAL8   & 8 & $-1.8\times10^{308}$ to $1.8\times10^{308}$ & \tt double \\
\hline
\end{tabular}
\end{center}

The unsigned character and integer datatypes store their values
according to the usual binary system.  For signed characters and
integers, setting the most-significant bit indicates that the number
formed from the remaining bits should be added to the lower value of
the range.  The \verb@REAL4@ and \verb@REAL8@ datatypes should store
values according to the IEEE Standard 754 for Binary Floating-Point
Arithmetic, which gives them the following precisions and dynamic
ranges:

\begin{center}
\begin{tabular}{|l|cc|}
\hline
                              & \tt     REAL4       & \tt     REAL8        \\
\hline
Minimum positive subnormal    & $1.4\times10^{-45}$ & $4.9\times10^{-324}$ \\
Minimum positive normal       & $1.2\times10^{-38}$ & $2.2\times10^{-308}$ \\
Maximum finite normal         & $3.4\times10^{38}$  & $1.8\times10^{308}$  \\
Minimum fractional difference & $6.0\times10^{-8}$  & $1.1\times10^{-16}$  \\
Significant decimal digits    & 6--9                & 15--17               \\
\hline
\end{tabular}
\end{center}

The minimum positive subnormal is the smallest positive representable
number.  The minimum positive normal is the smallest positive number
that can be represented with full precision; that is, one whose
mantissa lies in the range [0.5,1).  The maximum finite normal is the
largest representable number other than the reserved code for
$+\infty$.  The minimum fractional difference is the smallest
fractional difference between consecutive representable numbers, or
\emph{half} the difference between 1 and the next representable
number.  Significant decimal digits gives the number of decimal digits
used to represent the binary number in decimal notation: the first is
the maximum number of digits that are guaranteed not to change upon
conversion to binary, the second is the number of digits required to
represent a unique binary quantity.

</lalLaTeX> */

/** \var typedef char CHAR
 * \brief One-byte signed integer.
 */
/** \var typedef unsigned char UCHAR
 * \brief One-byte unsigned integer.
 */
/** \var typedef short INT2
 * \brief Two-byte signed integer.
 */
/** \var typedef short UINT2
 * \brief Two-byte unsigned integer.
 */
/** \var typedef int INT4
 * \brief Four-byte signed integer.
 */
/** \var typedef int UINT4
 * \brief Four-byte unsigned integer.
 */
/** \var typedef long long INT8
 * \brief Eight-byte signed integer.
 *
 * On some platforms this is equivalent to <tt>long int</tt> instead.
 */
/** \var typedef long long UINT8
 * \brief Eight-byte unsigned integer.
 *
 * On some platforms this is equivalent to <tt>unsigned long int</tt> instead.
 */
/** \var typedef float REAL4
 * \brief Single precision real floating-point number (4 bytes).
 */
/** \var typedef double REAL8
 * \brief Double precision real floating-point number (8 bytes).
 */
/** \def LAL_INT8_C(v) (v ## LL)
 * \brief Macro for use in defining \a v as an INT8 constant.
 *
 * This macro affixes the appropriate qualifier to form an INT8 constant.
 * For example:
 * \code
 * const INT8 jan_1_2000_gps_nanosec = LAL_INT8_C(63072001300000000)
 * \endcode
 */
/** \def LAL_UINT8_C(v) (v ## ULL)
 * \brief Macro for use in defining \a v as an UINT8 constant.
 *
 * This macro affixes the appropriate qualifier to form an UINT8 constant.
 * For example:
 * \code
 * const UINT8 jan_1_2000_gps_nanosec = LAL_UINT8_C(63072001300000000)
 * \endcode
 */

#ifndef _LALATOMICDATATYPES_H
#define _LALATOMICDATATYPES_H

#ifdef  __cplusplus
extern "C" {
#endif

typedef char CHAR;
typedef unsigned char UCHAR;
typedef unsigned char BOOLEAN;

#ifdef LDAS_BUILD

#include <lal/LDASConfig.h>

/* Integer types */

#if SIZEOF_SHORT == 2
  typedef short INT2;
  typedef unsigned short UINT2;
#elif SIZEOF_INT == 2
  typedef int INT2;
  typedef unsigned int UINT2;
#else
# error "ERROR: NO 2 BYTE INTEGER FOUND"
#endif

#if SIZEOF_INT == 4
  typedef int INT4;
  typedef unsigned int UINT4;
#elif SIZEOF_LONG == 4
  typedef long INT4;
  typedef unsigned long UINT4;
#else
# error "ERROR: NO 4 BYTE INTEGER FOUND"
#endif

#if SIZEOF_LONG == 8
  typedef long INT8;
  typedef unsigned long UINT8;
#elif SIZEOF_LONG_LONG == 8
  typedef long long INT8;
  typedef unsigned long long UINT8;
#else
# error "ERROR: NO 8 BYTE INTEGER FOUND"
#endif

/* Real types */

#if SIZEOF_FLOAT == 4
  typedef float REAL4;
#else
# error "ERROR: NO 4 BYTE REAL FOUND"
#endif

#if SIZEOF_DOUBLE == 8
  typedef double REAL8;
#else
# error "ERROR: NO 8 BYTE REAL FOUND"
#endif

#else /* ! LDAS_BUILD */

#include <lal/LALConfig.h>
#include <lal/LALRCSID.h>
NRCSID( LALATOMICDATATYPESH, "$Id$" );

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

/* could do this... but want this file to be independent of stdint.h */
/*
#ifdef HAVE_STDINT_H
#include <stdint.h>
typedef int16_t  INT2;
typedef int32_t  INT4;
typedef int64_t  INT8;
typedef uint16_t UINT2;
typedef uint32_t UINT4;
typedef uint64_t UINT8;
#else
*/

#if LAL_SIZEOF_SHORT == 2
  typedef short INT2;
  typedef unsigned short UINT2;
#elif LAL_SIZEOF_INT == 2
  typedef int INT2;
  typedef unsigned int UINT2;
#else
  typedef short INT2;
  typedef unsigned short UINT2;
# error "ERROR: NO 2 BYTE INTEGER FOUND"
#endif

#if LAL_SIZEOF_INT == 4
  typedef int INT4;
  typedef unsigned int UINT4;
#elif LAL_SIZEOF_LONG == 4
  typedef long INT4;
  typedef unsigned long UINT4;
#else
  typedef int INT4;
  typedef unsigned int UINT4;
# error "ERROR: NO 4 BYTE INTEGER FOUND"
#endif

#if LAL_SIZEOF_LONG == 8
  typedef long INT8;
  typedef unsigned long UINT8;
#elif LAL_SIZEOF_LONG_LONG == 8
#ifdef __GNUC__
  __extension__ typedef long long INT8;
  __extension__ typedef unsigned long long UINT8;
#else
  typedef long long INT8;
  typedef unsigned long long UINT8;
#endif /* __GNUC__ */
#else
  typedef long long INT8;
  typedef unsigned long long UINT8;
# error "ERROR: NO 8 BYTE INTEGER FOUND"
#endif

/* #endif */ /* HAVE_STDINT_H  -- commented out above */

/* Macros for integer constants */
#if LAL_SIZEOF_LONG == 8
#define LAL_INT8_C(v) (v ## L)
#define LAL_UINT8_C(v) (v ## UL) 
#elif LAL_SIZEOF_LONG_LONG == 8
#ifdef __GNUC__
#define LAL_INT8_C(v) (__extension__ v ## LL)
#define LAL_UINT8_C(v) (__extension__ v ## ULL)
#else
#define LAL_INT8_C(v) (v ## LL)
#define LAL_UINT8_C(v) (v ## ULL)
#endif
#else
#define LAL_INT8_C(v) (v ## LL)
#define LAL_UINT8_C(v) (v ## ULL)
# error "ERROR: NO 8 BYTE INTEGER FOUND"
#endif

/* Real types */

#if LAL_SIZEOF_FLOAT == 4
  typedef float REAL4;
#else
  typedef float REAL4;
# error "ERROR: NO 4 BYTE REAL FOUND"
#endif

#if LAL_SIZEOF_DOUBLE == 8
  typedef double REAL8;
#else
  typedef double REAL8;
# error "ERROR: NO 8 BYTE REAL FOUND"
#endif

#endif /* LDAS_BUILD */

/* <lalLaTeX>

\subsubsection*{Complex datatypes}

LAL represents complex numbers as structures with two floating-point
fields, storing the real and imaginary parts.  These are considered
primitive datatypes (rather than aggregate or structured datatypes)
because they conceptually represent a single number.  Furthermore,
atomic and complex datatypes are treated equivalently by LAL aggregate
and structured datatypes.

\vspace{1ex}
\begin{verbatim}
COMPLEX8
\end{verbatim}
This structure stores a single-precision complex number in 8~bytes of
memory.  The fields are:
\begin{description}
\item[\texttt{REAL4 re}] The real part.
\item[\texttt{REAL4 im}] The imaginary part.
\end{description}

\vspace{1ex}
\begin{verbatim}
COMPLEX16
\end{verbatim}
This structure stores a double-precision complex number in 16~bytes of
memory.  The fields are:
\begin{description}
\item[\texttt{REAL8 re}] The real part.
\item[\texttt{REAL8 im}] The imaginary part.
\end{description}

</lalLaTeX> */

/** Single-precision floating-point complex number (8 bytes total) */
typedef struct
tagCOMPLEX8
{
  REAL4 re; /**< The real part. */
  REAL4 im; /**< The imaginary part. */
}
COMPLEX8;

/** Double-precision floating-point complex number (16 bytes total) */
typedef struct
tagCOMPLEX16
{
  REAL8 re; /**< The real part. */
  REAL8 im; /**< The imaginary part. */
}
COMPLEX16;


/* <lalLaTeX>
\vfill{\footnotesize\input{LALAtomicDatatypesHV}}
</lalLaTeX> */


#ifdef  __cplusplus
}
#endif

#endif /* _LALATOMICDATATYPES_H */
