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
# error "ERROR: NO 2 BYTE INTEGER FOUND"
#endif

#if LAL_SIZEOF_INT == 4
  typedef int INT4;
  typedef unsigned int UINT4;
#elif LAL_SIZEOF_LONG == 4
  typedef long INT4;
  typedef unsigned long UINT4;
#else
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
# error "ERROR: NO 8 BYTE INTEGER FOUND"
#endif

/* Real types */

#if LAL_SIZEOF_FLOAT == 4
  typedef float REAL4;
#else
# error "ERROR: NO 4 BYTE REAL FOUND"
#endif

#if LAL_SIZEOF_DOUBLE == 8
  typedef double REAL8;
#else
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

typedef struct
tagCOMPLEX8
{
  REAL4 re;
  REAL4 im;
}
COMPLEX8;

typedef struct
tagCOMPLEX16
{
  REAL8 re;
  REAL8 im;
}
COMPLEX16;


/* <lalLaTeX>
\vfill{\footnotesize\input{LALAtomicDatatypesHV}}
</lalLaTeX> */


#ifdef  __cplusplus
}
#endif

#endif /* _LALATOMICDATATYPES_H */
