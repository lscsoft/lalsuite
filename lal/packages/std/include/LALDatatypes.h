/********************************* <lalVerbatim file="LALDatatypesHV">
Author: J. D. E. Creighton, T. D. Creighton
$Id$
********************************** </lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{LALDatatypes.h}}
\label{s:LALDatatypes.h}

Provides the basic LAL datatypes.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALDatatypes.h>
\end{verbatim}

\noindent This header defines the standard data types and data
structures that are used throughout LAL.  They fall into three general
categories: \emph{primitive} datatypes, \emph{aggregates} of primitive
datatypes, and \emph{structured} datatypes.  The LAL status structure
is a special case of a structured datatype that is used in every
standard LAL function.

This header file is automatically included by the header
\verb@LALStdlib.h@.  In turn, this header file starts by including the
header \verb@LALAtomicDatatypes.h@, which is discussed in the
following section.

</lalLaTeX> */

#ifndef _LALDATATYPES_H
#define _LALDATATYPES_H

/* <lalLaTeX>
\newpage\input{LALAtomicDatatypesH}
</lalLaTeX> */
#include <lal/LALAtomicDatatypes.h>

/**** <lalLaTeX>
 * \newpage
 * \subsection{Atomic datatypes codes}
 * \label{ss:atomic-datatype-codes}
 * \idx[Constant]{LAL\_1\_BYTE\_TYPE\_SIZE}
 * \idx[Constant]{LAL\_2\_BYTE\_TYPE\_SIZE}
 * \idx[Constant]{LAL\_4\_BYTE\_TYPE\_SIZE}
 * \idx[Constant]{LAL\_8\_BYTE\_TYPE\_SIZE}
 * \idx[Constant]{LAL\_16\_BYTE\_TYPE\_SIZE}
 * \idx[Constant]{LAL\_TYPE\_SIZE\_MASK}
 * \idx[Constant]{LAL\_FLTPT\_TYPE\_FLAG}
 * \idx[Constant]{LAL\_CMPLX\_TYPE\_FLAG}
 * \idx[Constant]{LAL\_UNSGN\_TYPE\_FLAG}
 * \idx[Constant]{LAL\_CHAR\_TYPE\_CODE}
 * \idx[Constant]{LAL\_I2\_TYPE\_CODE}
 * \idx[Constant]{LAL\_I4\_TYPE\_CODE}
 * \idx[Constant]{LAL\_I8\_TYPE\_CODE}
 * \idx[Constant]{LAL\_UCHAR\_TYPE\_CODE}
 * \idx[Constant]{LAL\_U2\_TYPE\_CODE}
 * \idx[Constant]{LAL\_U4\_TYPE\_CODE}
 * \idx[Constant]{LAL\_U8\_TYPE\_CODE}
 * \idx[Constant]{LAL\_S\_TYPE\_CODE}
 * \idx[Constant]{LAL\_D\_TYPE\_CODE}
 * \idx[Constant]{LAL\_C\_TYPE\_CODE}
 * \idx[Constant]{LAL\_Z\_TYPE\_CODE}
 * \idx[Type]{LALTYPECODE}
 *
 * The following constants specify the size, in bytes, of the atomic datatype.
 *
 * \begin{center}
 * \begin{tabular}{|lcl|}
 * \hline
 * Name & Octal Value & Description \\
 * \hline
 * \tt LAL\_1\_BYTE\_TYPE\_SIZE & 000 & 1 byte type \\
 * \tt LAL\_2\_BYTE\_TYPE\_SIZE & 001 & 2 byte type \\
 * \tt LAL\_4\_BYTE\_TYPE\_SIZE & 002 & 4 byte type \\
 * \tt LAL\_8\_BYTE\_TYPE\_SIZE & 003 & 8 byte type \\
 * \tt LAL\_16\_BYTE\_TYPE\_SIZE & 004 & 16 byte type \\
 * \tt LAL\_TYPE\_SIZE\_MASK & 007 & Mask for byte type size fields \\
 * \hline
 * \end{tabular}
 * \end{center}
 *
 * \noindent
 * The constant \verb+LAL_TYPE_SIZE_MASK+ is useful in extracting the size
 * information from other type attributes.  For example, the size, in bytes,
 * of an atomic datatype can be found using something like the following:
 * \begin{verbatim}
 * UINT4 code = LAL_S_TYPE_CODE;
 * UINT4 size = 1U << ( code & LAL_TYPE_SIZE_MASK );
 * \end{verbatim}
 *
 * \vspace{3ex}
 *
 * The following constants are flags describing the type attributes.  A type
 * is either an integer or a floating-point, either purely real or complex,
 * and, if integer, is either signed or unsigned.
 *
 * \begin{center}
 * \begin{tabular}{|lcl|}
 * \hline
 * Name & Octal Value & Description \\
 * \hline
 * \tt LAL\_FLTPT\_TYPE\_FLAG & 010 & Floating-point (not integer) type \\
 * \tt LAL\_CMPLX\_TYPE\_FLAG & 020 & Complex (not purely real) type \\
 * \tt LAL\_UNSGN\_TYPE\_FLAG & 040 & Unsigned (no sign info) type \\
 * \hline
 * \end{tabular}
 * \end{center}
 *
 * To get the actual type, these flags are combined together and with the
 * type size constants using the bitwise-or operator (\verb+|+).  For example,
 * an eight-byte floating point number would be
 * \verb+LAL_8_BYTE_TYPE_SIZE | LAL_FLTPT_TYPE_FLAG+.
 * Conceivably you could have a complex type made from a pair of unsigned
 * one-byte integers that would be specified as
 * \verb+LAL_1_BYTE_TYPE_SIZE | LAL_CMPLX_TYPE_FLAG | LAL_UNSGN_TYPE_FLAG+.
 * Fortunately, there are none of these in LAL.  Attribues of a particular
 * type can be extracted using the bitwise-and operator.  For example:
 * \begin{verbatim}
 * UINT4 code = LAL_S_TYPE_CODE;
 * UINT4 isfloat = ( code & LAL_FLTPT_TYPE_FLAG );
 * UINT4 iscmplx = ( code & LAL_CMPLX_TYPE_FLAG );
 * \end{verbatim}
 *
 * \vspace{3ex}
 *
 * The following constants correspond to the types that actually exist in LAL.
 * Their enumeration is the type \verb+LALTYPECODE+.
 * \begin{center}
 * \begin{tabular}{|lcl|}
 * \hline
 * Name & Octal Value & Corresponding Type \\
 * \hline
 * \tt LAL\_CHAR\_TYPE\_CODE & 000 & \tt CHAR \\
 * \tt LAL\_I2\_TYPE\_CODE & 001 & \tt INT2 \\
 * \tt LAL\_I4\_TYPE\_CODE & 002 & \tt INT4 \\
 * \tt LAL\_I8\_TYPE\_CODE & 003 & \tt INT8 \\
 * \tt LAL\_UCHAR\_TYPE\_CODE & 040 & \tt UCHAR \\
 * \tt LAL\_U2\_TYPE\_CODE & 041 & \tt UINT2 \\
 * \tt LAL\_U4\_TYPE\_CODE & 042 & \tt UINT4 \\
 * \tt LAL\_U8\_TYPE\_CODE & 043 & \tt UINT8 \\
 * \tt LAL\_S\_TYPE\_CODE & 012 & \tt REAL4 \\
 * \tt LAL\_D\_TYPE\_CODE & 013 & \tt REAL8 \\
 * \tt LAL\_C\_TYPE\_CODE & 033 & \tt COMPLEX8 \\
 * \tt LAL\_Z\_TYPE\_CODE & 034 & \tt COMPLEX16 \\
 * \hline
 * \end{tabular}
 * \end{center}
 *
 **** </lalLaTeX> */

/* constants */
/* type size constants */
enum
{
  LAL_1_BYTE_TYPE_SIZE  = 000,   /*      0 = 0 */
  LAL_2_BYTE_TYPE_SIZE  = 001,   /*      1 = 1 */
  LAL_4_BYTE_TYPE_SIZE  = 002,   /*     10 = 2 */
  LAL_8_BYTE_TYPE_SIZE  = 003,   /*     11 = 3 */
  LAL_16_BYTE_TYPE_SIZE = 004,   /*    100 = 4 */
  LAL_TYPE_SIZE_MASK    = 007    /*    111 = 7 */
};

/* type flag constants */
enum
{
  LAL_FLTPT_TYPE_FLAG   = 010,   /*   1000 =  8 */
  LAL_CMPLX_TYPE_FLAG   = 020,   /*  10000 = 16 */
  LAL_UNSGN_TYPE_FLAG   = 040    /* 100000 = 32 */
};

/* type codes */
typedef enum
{
  LAL_CHAR_TYPE_CODE    = LAL_1_BYTE_TYPE_SIZE, /* 0 */
  LAL_I2_TYPE_CODE      = LAL_2_BYTE_TYPE_SIZE, /* 1 */
  LAL_I4_TYPE_CODE      = LAL_4_BYTE_TYPE_SIZE, /* 2 */
  LAL_I8_TYPE_CODE      = LAL_8_BYTE_TYPE_SIZE, /* 3 */
  LAL_UCHAR_TYPE_CODE   = LAL_1_BYTE_TYPE_SIZE | LAL_UNSGN_TYPE_FLAG,   /* 32 */
  LAL_U2_TYPE_CODE      = LAL_2_BYTE_TYPE_SIZE | LAL_UNSGN_TYPE_FLAG,   /* 33 */
  LAL_U4_TYPE_CODE      = LAL_4_BYTE_TYPE_SIZE | LAL_UNSGN_TYPE_FLAG,   /* 34 */
  LAL_U8_TYPE_CODE      = LAL_8_BYTE_TYPE_SIZE | LAL_UNSGN_TYPE_FLAG,   /* 35 */
  LAL_S_TYPE_CODE       = LAL_4_BYTE_TYPE_SIZE | LAL_FLTPT_TYPE_FLAG,   /* 18 */
  LAL_D_TYPE_CODE       = LAL_8_BYTE_TYPE_SIZE | LAL_FLTPT_TYPE_FLAG,   /* 19 */
  LAL_C_TYPE_CODE       = LAL_8_BYTE_TYPE_SIZE | LAL_CMPLX_TYPE_FLAG | LAL_FLTPT_TYPE_FLAG,     /* 27 */
  LAL_Z_TYPE_CODE       = LAL_16_BYTE_TYPE_SIZE | LAL_CMPLX_TYPE_FLAG | LAL_FLTPT_TYPE_FLAG     /* 28 */
}
LALTYPECODE;


#include <lal/LALRCSID.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (LALDATATYPESH, "$Id$");


/* <lalLaTeX>
\newpage
\subsection{Aggregate datatypes}
\label{ss:aggregate-datatypes}

These datatypes store arbitrarily large sets or collections of
primitive datatypes.  At this level there is no physical
interpretation assigned to the objects (such as names or units); the
aggregate datatypes simply collect and arrange the primitive
datatypes.  The following types of aggregate datatypes are defines:
vectors, arrays, sequences, vector sequences, and array sequences.

\vspace{2ex}
\begin{verbatim}
<datatype>Vector
\end{verbatim}
This structure stores an ordered set of $n$ elements of type
\verb@<datatype>@, which can be any primitive datatype.  The data are
to be interpreted as being a point in an $n$-dimensional vector space.
The fields are:
\begin{description}
\item[\texttt{UINT4 length}] The number of data $n$.
\item[\texttt{<datatype> *data}] Pointer to the data array.  The data
are stored sequentially as \verb@data[@$0,\ldots,n-1$\verb@]@.
\end{description}

</lalLaTeX> */

typedef struct
tagCHARVector
{
  UINT4  length;
  CHAR  *data;
}
CHARVector;

typedef struct
tagINT2Vector
{
  UINT4  length;
  INT2  *data;
}
INT2Vector;

typedef struct
tagUINT2Vector
{
  UINT4  length;
  UINT2 *data;
}
UINT2Vector;

typedef struct
tagINT4Vector
{
  UINT4  length;
  INT4  *data;
}
INT4Vector;

typedef struct
tagUINT4Vector
{
  UINT4  length;
  UINT4  *data;
}
UINT4Vector;

typedef struct
tagINT8Vector
{
  UINT4  length;
  INT8  *data;
}
INT8Vector;

typedef struct
tagUINT8Vector
{
  UINT4  length;
  UINT8 *data;
}
UINT8Vector;

typedef struct
tagREAL4Vector
{
  UINT4  length;
  REAL4 *data;
}
REAL4Vector;

typedef struct tagREAL8Vector
{
  UINT4  length;
  REAL8 *data;
}
REAL8Vector;

typedef struct tagCOMPLEX8Vector
{
  UINT4     length;
  COMPLEX8 *data;
}
COMPLEX8Vector;

typedef struct tagCOMPLEX16Vector
{
  UINT4      length;
  COMPLEX16 *data;
}
COMPLEX16Vector;


/* <lalLaTeX>

\vspace{2ex}
\begin{verbatim}
<datatype>Array
\end{verbatim}
This structure stores a set of elements of type \verb@<datatype>@,
which can be any primitive datatype, arranged as an $m$-dimensional
array.  That is, each element can be thought of as having $m$
indecies, $\mathsf{A}_{i_0\cdots i_{m-1}}$, where each index $i_k$
runs over its own range $0,\ldots,n_k-1$.  The total number of
elements is then $N=n_0\times\cdots\times n_{m-1}$.  In memory the
array is ``flattened'' so that the elements are stored sequentially in
a contiguous block.  The fields are:
\begin{description}
\item[\texttt{UINT4Vector *dimLength}] Pointer to a vector of length
$m$, storing the index ranges $(n_0,\ldots,n_{m-1})$.
\item[\texttt{<datatype> *data}] Pointer to the data array.  The data
element $\mathsf{A}_{i_0\cdots i_{m-1}}$ is stored as
\verb@data[@$i_{m-1} + n_{m-2}\times(i_{m-2} +
n_{m-3}\times(\cdots(i_1 + n_0\times i_0)\cdots))$\verb@]@; that is,
the index of \verb@data[]@ runs over the entire range of an index
$i_{k+1}$ before incrementing $i_k$.
\end{description}

</lalLaTeX> */

typedef struct
tagINT2Array
{
  UINT4Vector *dimLength;
  INT2        *data;
}
INT2Array;

typedef struct
tagUINT2Array
{
  UINT4Vector *dimLength;
  UINT2       *data;
}
UINT2Array;

typedef struct
tagINT4Array
{
  UINT4Vector *dimLength;
  INT4        *data;
}
INT4Array;

typedef struct
tagUINT4Array
{
  UINT4Vector *dimLength;
  UINT4       *data;
}
UINT4Array;

typedef struct
tagINT8Array
{
  UINT4Vector *dimLength;
  INT8        *data;
}
INT8Array;

typedef struct
tagUINT8Array
{
  UINT4Vector *dimLength;
  UINT8       *data;
}
UINT8Array;

typedef struct
tagREAL4Array
{
  UINT4Vector *dimLength;
  REAL4       *data;
}
REAL4Array;

typedef struct
tagREAL8Array
{
  UINT4Vector *dimLength;
  REAL8       *data;
}
REAL8Array;

typedef struct
tagCOMPLEX8Array
{
  UINT4Vector *dimLength;
  COMPLEX8    *data;
}
COMPLEX8Array;

typedef struct
tagCOMPLEX16Array
{
  UINT4Vector *dimLength;
  COMPLEX16   *data;
}
COMPLEX16Array;


/* <lalLaTeX>

\vspace{2ex}
\begin{verbatim}
<datatype>Sequence
\end{verbatim}
This structure stores an ordered set of $l$ elements of type
\verb@<datatype>@, which can be any primitive datatype.  It is
identical to \verb@<datatype>Vector@, except that the elements are to
be interpreted as $l$ consecutive elements rather than the components
of an $l$-dimensional vector.  The fields are:
\begin{description}
\item[\texttt{UINT4 length}] The number of data $l$.
\item[\texttt{<datatype> *data}] Pointer to the data array.  The data
are stored sequentially as \verb@data[@$0,\ldots,l-1$\verb@]@.
\end{description}

</lalLaTeX> */

typedef CHARVector      CHARSequence;
typedef INT2Vector      INT2Sequence;
typedef UINT2Vector     UINT2Sequence;
typedef INT4Vector      INT4Sequence;
typedef UINT4Vector     UINT4Sequence;
typedef INT8Vector      INT8Sequence;
typedef UINT8Vector     UINT8Sequence;
typedef REAL4Vector     REAL4Sequence;
typedef REAL8Vector     REAL8Sequence;
typedef COMPLEX8Vector  COMPLEX8Sequence;
typedef COMPLEX16Vector COMPLEX16Sequence;


/* <lalLaTeX>

\vspace{2ex}
\begin{verbatim}
<datatype>VectorSequence
\end{verbatim}
This structure stores an ordered set of $l$ elements of type
\verb@<datatype>Vector@, where \verb@<datatype>@ can be any primitive
datatype.  Mathematically the sequence can be written as
$\{\vec{v}^{(0)},\ldots,\vec{v}^{(l-1)}\}$, where each element
$\vec{v}^{(j)}=(v^{(j)}_0,\ldots,v^{(i)}_{n-1})$ is a vector of length
$n$.  In memory the elements are ``flattened''; that is, they are
stored sequentially in a contiguous block of memory.  The fields are:
\begin{description}
\item[\texttt{UINT4 length}] The number of vectors $l$.
\item[\texttt{UINT4 vectorLength}] The length $n$ of each vector.
\item[\texttt{<datatype> *data}] Pointer to the data array.  The data
element $v^{(j)}_i$ is stored as \verb@data[@$j\times n + i$\verb@]@;
that is, the index of \verb@data[]@ runs over the internal index of
each vector element before incrementing to the next vector element.
\end{description}

</lalLaTeX> */

typedef struct
tagCHARVectorSequence
{
  UINT4  length;
  UINT4  vectorLength;
  CHAR  *data;
}
CHARVectorSequence;

typedef struct
tagINT2VectorSequence
{
  UINT4  length;
  UINT4  vectorLength;
  INT2  *data;
}
INT2VectorSequence;

typedef struct
tagUINT2VectorSequence
{
  UINT4  length;
  UINT4  vectorLength;
  UINT2 *data;
}
UINT2VectorSequence;

typedef struct
tagINT4VectorSequence
{
  UINT4  length;
  UINT4  vectorLength;
  INT4  *data;
}
INT4VectorSequence;

typedef struct
tagUINT4VectorSequence
{
  UINT4  length;
  UINT4  vectorLength;
  UINT4 *data;
}
UINT4VectorSequence;

typedef struct
tagINT8VectorSequence
{
  UINT4  length;
  UINT4  vectorLength;
  INT8  *data;
}
INT8VectorSequence;

typedef struct
tagUINT8VectorSequence
{
  UINT4  length;
  UINT4  vectorLength;
  UINT8 *data;
}
UINT8VectorSequence;

typedef struct
tagREAL4VectorSequence
{
  UINT4  length;
  UINT4  vectorLength;
  REAL4 *data;
}
REAL4VectorSequence;

typedef struct
tagREAL8VectorSequence
{
  UINT4  length;
  UINT4  vectorLength;
  REAL8 *data;
}
REAL8VectorSequence;

typedef struct
tagCOMPLEX8VectorSequence
{
  UINT4     length;
  UINT4     vectorLength;
  COMPLEX8 *data;
}
COMPLEX8VectorSequence;

typedef struct
tagCOMPLEX16VectorSequence
{
  UINT4      length;
  UINT4      vectorLength;
  COMPLEX16 *data;
}
COMPLEX16VectorSequence;


/* <lalLaTeX>

\vspace{2ex}
\begin{verbatim}
<datatype>ArraySequence
\end{verbatim}
This structure stores an ordered set of $l$ elements of type
\verb@<datatype>Array@, where \verb@<datatype>@ can be any primitive
datatype.  The indexing of an array sequence can get quite
complicated; it helps to read first the documentation for data arrays,
above.  Mathematically the data can be written as a set
$\{\mathsf{A}^{(j)}_{i_0\cdots i_{m-1}}$, where the sequence number
$j$ runs from 0 to $l-1$, and each array index $i_k$ runs over its own
range $0,\ldots,n_k-1$.  The total number of data in a given array
element is then $N=n_0\times\cdots\times n_{m-1}$, and the total
number of data in the sequence is $N\times l$.  In memory the array is
``flattened'' so that the elements are stored sequentially in a
contiguous block.  The fields are:
\begin{description}
\item[\texttt{UINT4 length}] The number $l$ of array elements in the
sequence.
\item[\texttt{UINT4 arrayDim}] The number of data $N$ (\emph{not} the
number of indecies $m$) in each array element of the sequence.
\item[\texttt{UINT4Vector *dimLength}] Pointer to a vector of length
$m$, storing the index ranges $(n_0,\ldots,n_{m-1})$.
\item[\texttt{<datatype> *data}] Pointer to the data.  The element
$\mathsf{A}^{(j)}_{i_0\cdots i_{m-1}}$ is stored as
\verb@data[@$j\times N + i_{m-1} + n_{m-2}\times(i_{m-2} +
n_{m-3}\times(\cdots(i_1 + n_0\times i_0)\cdots))$\verb@]@; that is,
the index of \verb@data[]@ runs over the internal indecies of each
array element before incrementing to the next array element.
\end{description}

</lalLaTeX> */

typedef struct
tagINT2ArraySequence
{
  UINT4        length;
  UINT4        arrayDim;
  UINT4Vector *dimLength;
  INT2        *data;
}
INT2ArraySequence;

typedef struct
tagUINT2ArraySequence
{
  UINT4        length;
  UINT4        arrayDim;
  UINT4Vector *dimLength;
  UINT2       *data;
}
UINT2ArraySequence;

typedef struct
tagINT4ArraySequence
{
  UINT4        length;
  UINT4        arrayDim;
  UINT4Vector *dimLength;
  INT4        *data;
}
INT4ArraySequence;

typedef struct
tagUINT4ArraySequence
{
  UINT4        length;
  UINT4        arrayDim;
  UINT4Vector *dimLength;
  UINT4       *data;
}
UINT4ArraySequence;

typedef struct
tagINT8ArraySequence
{
  UINT4        length;
  UINT4        arrayDim;
  UINT4Vector *dimLength;
  INT8        *data;
}
INT8ArraySequence;

typedef struct
tagUINT8ArraySequence
{
  UINT4        length;
  UINT4        arrayDim;
  UINT4Vector *dimLength;
  UINT8       *data;
}
UINT8ArraySequence;

typedef struct
tagREAL4ArraySequence
{
  UINT4        length;
  UINT4        arrayDim;
  UINT4Vector *dimLength;
  REAL4       *data;
}
REAL4ArraySequence;

typedef struct
tagREAL8ArraySequence
{
  UINT4        length;
  UINT4        arrayDim;
  UINT4Vector *dimLength;
  REAL8       *data;
}
REAL8ArraySequence;

typedef struct
tagCOMPLEX8ArraySequence
{
  UINT4        length;
  UINT4        arrayDim;
  UINT4Vector *dimLength;
  COMPLEX8    *data;
}
COMPLEX8ArraySequence;

typedef struct
tagCOMPLEX16ArraySequence
{
  UINT4        length;
  UINT4        arrayDim;
  UINT4Vector *dimLength;
  COMPLEX16   *data;
}
COMPLEX16ArraySequence;


/* <lalLaTeX>
\newpage
\subsection{Structured datatypes}
\label{ss:structured-datatypes}

These datatypes embed primitive and aggregate datatypes inside
structures that define their physical meaning.  Most of these
structures are wrappers for aggregate datatypes that store a physical
quantity as a function of time or frequency.  Other structures store
specific physical information, such as the GPS time, or the factored
response function of a filter.

\vspace{2ex}
\begin{verbatim}
LIGOTimeGPS
\end{verbatim}
This structure stores the time, to nanosecond precision, synchronized
to the Global Positioning System time reference.  The zero time for
the GPS standard is the moment of midnight beginning January~6, 1980,
UTC.  The \verb@LIGOTimeGPS@ structure can represent times up to
68~years on either side of this epoch.  (Note that this is better than
an equivalently-sized \verb@REAL8@ representation of time, which can
maintain nanosecond precision only for times within 104~days of its
reference point.  However, the \verb@REAL8@ representation does allow
one to cover arbitrarily long timescales at correspondingly lower
precision.)  The fields are:
\begin{description}
\item[\texttt{INT4 gpsSeconds}] The number of seconds since the GPS
reference time.
\item[\texttt{INT4 gpsNanoSeconds}] The number of nanoseconds since
the last GPS second.
\end{description}

The macro \verb@LIGOTIMEGPSZERO@ can be used to statically initialize a
\verb@LIGOTimeGPS@ object, for example:
\begin{quote}
\verb@LIGOTimeGPS epoch = LIGOTIMEGPSZERO;@
\end{quote}

</lalLaTeX> */

typedef struct
tagLIGOTimeGPS
{
  INT4 gpsSeconds;
  INT4 gpsNanoSeconds;
}
LIGOTimeGPS;

#define LIGOTIMEGPSZERO { 0, 0 }

/* <lalLaTeX>

\vspace{2ex}
\begin{verbatim}
LALUnit
\end{verbatim}
This structure stores units in the mksA system (plus Kelvin, Strain,
and ADC Count).  It also stores an overall power-of-ten scaling factor.
The fields are:
\begin{description}
\item[\texttt{INT2 powerOfTen}] The power $p$ of ten scaling factor.
\item[\texttt{INT2 unitNumerator[LALNumUnits]}] Array of unit numerators,
  $N_i$, $i=0\ldots\textrm{LALNumUnits}-1$.
\item[\texttt{INT2 unitDenominatorMinusOne[LALNumUnits]}] Array of unit
  denominators-minus-one, $D_i$, $i=0\ldots\textrm{LALNumUnits}-1$.
\end{description}
Thus, the units are given by
\begin{equation}
  10^p\times\textrm{m}^{N_0/(1+D_0)}\times\textrm{kg}^{N_1/(1+D_1)}
  \times\textrm{s}^{N_2/(1+D_2)}\times\textrm{A}^{N_3/(1+D_3)}
  \times\textrm{K}^{N_4/(1+D_4)}\times\textrm{strain}^{N_5/(1+D_5)}
  \times\textrm{count}^{N_6/(1+D_6)}
\end{equation}
The indexes of the units can be specified using the constants
\texttt{LALUnitIndexMeter},
\texttt{LALUnitIndexKiloGram},
\texttt{LALUnitIndexSecond},
\texttt{LALUnitIndexAmpere},
\texttt{LALUnitIndexKelvin},
\texttt{LALUnitIndexStrain},
\texttt{LALUnitIndexADCCount},
while \texttt{LALNumUnits} is the total number of units.

</lalLaTeX> */

enum
{
  LALUnitIndexMeter,
  LALUnitIndexKiloGram,
  LALUnitIndexSecond,
  LALUnitIndexAmpere,
  LALUnitIndexKelvin,
  LALUnitIndexStrain,
  LALUnitIndexADCCount,
  LALNumUnits
};

typedef struct
tagLALUnit
{
  INT2  powerOfTen;
  INT2  unitNumerator[LALNumUnits];
  UINT2 unitDenominatorMinusOne[LALNumUnits];
}
LALUnit;


/* <lalLaTeX>

\vspace{2ex}
\begin{verbatim}
<datatype>TimeSeries
\end{verbatim}
This structure represents a sequence of data of type \verb@<datatype>@
(where \verb@<datatype>@ can be any primitive datatype), sampled over
uniform time intervals $t_0, t_0+\Delta t, \ldots , t_0+l\Delta t$.
Essentially this is a \verb@<datatype>Sequence@ with extra fields
defining the sample times and the type of data being sampled.  The raw
data may also have been \emph{heterodyned}; that is, multiplied by a
sinusoid of some frequency $f_0$, low-pass filtered, and resampled, in
order to extract the behaviour in a small bandwidth about $f_0$.  The
fields are:
\begin{description}
\item[\texttt{CHAR name[LALNameLength]}] The name of the data series (i.e.\
the type of data being sampled).
\item[\texttt{LIGOTimeGPS epoch}] The start time $t_0$ of the data
series.
\item[\texttt{REAL8 deltaT}] The sampling interval $\Delta t$, in
seconds.
\item[\texttt{REAL8 f0}] The heterodyning frequency $f_0$, in hertz.
\item[\texttt{LALUnit sampleUnits}] The physical units of the
quantity being sampled.
\item[\texttt{<datatype>Sequence *data}] The sequence of sampled data.
\end{description}

</lalLaTeX> */

enum { LALNameLength = 64 };

typedef struct
tagINT2TimeSeries
{
  CHAR          name[LALNameLength];
  LIGOTimeGPS   epoch;
  REAL8         deltaT;
  REAL8         f0;
  LALUnit       sampleUnits;
  INT2Sequence *data;
}
INT2TimeSeries;

typedef struct
tagUINT2TimeSeries
{
  CHAR           name[LALNameLength];
  LIGOTimeGPS    epoch;
  REAL8          deltaT;
  REAL8          f0;
  LALUnit        sampleUnits;
  UINT2Sequence *data;
}
UINT2TimeSeries;

typedef struct
tagINT4TimeSeries
{
  CHAR          name[LALNameLength];
  LIGOTimeGPS   epoch;
  REAL8         deltaT;
  REAL8         f0;
  LALUnit       sampleUnits;
  INT4Sequence *data;
}
INT4TimeSeries;

typedef struct
tagUINT4TimeSeries
{
  CHAR           name[LALNameLength];
  LIGOTimeGPS    epoch;
  REAL8          deltaT;
  REAL8          f0;
  LALUnit        sampleUnits;
  UINT4Sequence *data;
}
UINT4TimeSeries;

typedef struct
tagINT8TimeSeries
{
  CHAR          name[LALNameLength];
  LIGOTimeGPS   epoch;
  REAL8         deltaT;
  REAL8         f0;
  LALUnit       sampleUnits;
  INT8Sequence *data;
}
INT8TimeSeries;

typedef struct
tagUINT8TimeSeries
{
  CHAR           name[LALNameLength];
  LIGOTimeGPS    epoch;
  REAL8          deltaT;
  REAL8          f0;
  LALUnit        sampleUnits;
  UINT8Sequence *data;
}
UINT8TimeSeries;

typedef struct
tagREAL4TimeSeries
{
  CHAR           name[LALNameLength];
  LIGOTimeGPS    epoch;
  REAL8          deltaT;
  REAL8          f0;
  LALUnit        sampleUnits;
  REAL4Sequence *data;
}
REAL4TimeSeries;

typedef struct
tagREAL8TimeSeries
{
  CHAR           name[LALNameLength];
  LIGOTimeGPS    epoch;
  REAL8          deltaT;
  REAL8          f0;
  LALUnit        sampleUnits;
  REAL8Sequence *data;
}
REAL8TimeSeries;

typedef struct
tagCOMPLEX8TimeSeries
{
  CHAR              name[LALNameLength];
  LIGOTimeGPS       epoch;
  REAL8             deltaT;
  REAL8             f0;
  LALUnit           sampleUnits;
  COMPLEX8Sequence *data;
}
COMPLEX8TimeSeries;

typedef struct
tagCOMPLEX16TimeSeries
{
  CHAR               name[LALNameLength];
  LIGOTimeGPS        epoch;
  REAL8              deltaT;
  REAL8              f0;
  LALUnit            sampleUnits;
  COMPLEX16Sequence *data;
}
COMPLEX16TimeSeries;


/* <lalLaTeX>

\vspace{2ex}
\begin{verbatim}
<datatype>TimeVectorSeries
\end{verbatim}
Like \verb@<datatype>TimeSeries@, above, except that the sampled data
are of type type \verb@<datatype>Vector@ (where \verb@<datatype>@ can
be any primitive datatype).  The fields are:
\begin{description}
\item[\texttt{CHAR name[LALNameLength]}] The name of the data series (i.e.\
the type of data being sampled).
\item[\texttt{LIGOTimeGPS epoch}] The start time of the data series.
\item[\texttt{REAL8 deltaT}] The sampling interval, in seconds.
\item[\texttt{REAL8 f0}] The heterodyning frequency, in hertz.
\item[\texttt{LALUnit sampleUnits}] The physical units of the
quantity being sampled.
\item[\texttt{<datatype>VectorSequence *data}] The sequence of sampled
data.
\end{description}

</lalLaTeX> */

typedef struct
tagINT2TimeVectorSeries
{
  CHAR                 name[LALNameLength];
  LIGOTimeGPS          epoch;
  REAL8                deltaT;
  REAL8                f0;
  LALUnit              sampleUnits;
  INT2VectorSequence  *data;
}
INT2TimeVectorSeries;

typedef struct
tagUINT2TimeVectorSeries
{
  CHAR                 name[LALNameLength];
  LIGOTimeGPS          epoch;
  REAL8                deltaT;
  REAL8                f0;
  LALUnit              sampleUnits;
  UINT2VectorSequence *data;
}
UINT2TimeVectorSeries;

typedef struct
tagINT4TimeVectorSeries
{
  CHAR                 name[LALNameLength];
  LIGOTimeGPS          epoch;
  REAL8                deltaT;
  REAL8                f0;
  LALUnit              sampleUnits;
  INT4VectorSequence  *data;
}
INT4TimeVectorSeries;

typedef struct
tagUINT4TimeVectorSeries
{
  CHAR                 name[LALNameLength];
  LIGOTimeGPS          epoch;
  REAL8                deltaT;
  REAL8                f0;
  LALUnit              sampleUnits;
  UINT4VectorSequence *data;
}
UINT4TimeVectorSeries;

typedef struct
tagINT8TimeVectorSeries
{
  CHAR                 name[LALNameLength];
  LIGOTimeGPS          epoch;
  REAL8                deltaT;
  REAL8                f0;
  LALUnit              sampleUnits;
  INT8VectorSequence  *data;
}
INT8TimeVectorSeries;

typedef struct
tagUINT8TimeVectorSeries
{
  CHAR                 name[LALNameLength];
  LIGOTimeGPS          epoch;
  REAL8                deltaT;
  REAL8                f0;
  LALUnit              sampleUnits;
  UINT8VectorSequence *data;
}
UINT8TimeVectorSeries;

typedef struct
tagREAL4TimeVectorSeries
{
  CHAR                 name[LALNameLength];
  LIGOTimeGPS          epoch;
  REAL8                deltaT;
  REAL8                f0;
  LALUnit              sampleUnits;
  REAL4VectorSequence *data;
}
REAL4TimeVectorSeries;

typedef struct
tagREAL8TimeVectorSeries
{
  CHAR                 name[LALNameLength];
  LIGOTimeGPS          epoch;
  REAL8                deltaT;
  REAL8                f0;
  LALUnit              sampleUnits;
  REAL8VectorSequence *data;
}
REAL8TimeVectorSeries;

typedef struct
tagCOMPLEX8TimeVectorSeries
{
  CHAR                     name[LALNameLength];
  LIGOTimeGPS              epoch;
  REAL8                    deltaT;
  REAL8                    f0;
  LALUnit                  sampleUnits;
  COMPLEX8VectorSequence  *data;
}
COMPLEX8TimeVectorSeries;

typedef struct
tagCOMPLEX16TimeVectorSeries
{
  CHAR                      name[LALNameLength];
  LIGOTimeGPS               epoch;
  REAL8                     deltaT;
  REAL8                     f0;
  LALUnit                   sampleUnits;
  COMPLEX16VectorSequence  *data;
}
COMPLEX16TimeVectorSeries;


/* <lalLaTeX>

\vspace{2ex}
\begin{verbatim}
<datatype>TimeArraySeries
\end{verbatim}
Like \verb@<datatype>TimeSeries@, above, except that the sampled data
are of type type \verb@<datatype>Array@ (where \verb@<datatype>@ can
be any primitive datatype).  The fields are:
\begin{description}
\item[\texttt{CHAR name[LALNameLength]}] The name of the data series (i.e.\
the type of data being sampled).
\item[\texttt{LIGOTimeGPS epoch}] The start time of the data series.
\item[\texttt{REAL8 deltaT}] The sampling interval, in seconds.
\item[\texttt{REAL8 f0}] The heterodyning frequency, in hertz.
\item[\texttt{LALUnit sampleUnits}] The physical units of the
quantity being sampled.
\item[\texttt{<datatype>ArraySequence *data}] The sequence of sampled
data.
\end{description}

</lalLaTeX> */

typedef struct
tagINT2TimeArraySeries
{
  CHAR                name[LALNameLength];
  LIGOTimeGPS         epoch;
  REAL8               deltaT;
  REAL8               f0;
  LALUnit             sampleUnits;
  INT2ArraySequence  *data;
}
INT2TimeArraySeries;

typedef struct
tagUINT2TimeArraySeries
{
  CHAR                name[LALNameLength];
  LIGOTimeGPS         epoch;
  REAL8               deltaT;
  REAL8               f0;
  LALUnit             sampleUnits;
  UINT2ArraySequence *data;
}
UINT2TimeArraySeries;

typedef struct
tagINT4TimeArraySeries
{
  CHAR                name[LALNameLength];
  LIGOTimeGPS         epoch;
  REAL8               deltaT;
  REAL8               f0;
  LALUnit             sampleUnits;
  INT4ArraySequence  *data;
}
INT4TimeArraySeries;

typedef struct
tagUINT4TimeArraySeries
{
  CHAR                name[LALNameLength];
  LIGOTimeGPS         epoch;
  REAL8               deltaT;
  REAL8               f0;
  LALUnit             sampleUnits;
  UINT4ArraySequence *data;
}
UINT4TimeArraySeries;

typedef struct
tagINT8TimeArraySeries
{
  CHAR                name[LALNameLength];
  LIGOTimeGPS         epoch;
  REAL8               deltaT;
  REAL8               f0;
  LALUnit             sampleUnits;
  INT8ArraySequence  *data;
}
INT8TimeArraySeries;

typedef struct
tagUINT8TimeArraySeries
{
  CHAR                name[LALNameLength];
  LIGOTimeGPS         epoch;
  REAL8               deltaT;
  REAL8               f0;
  LALUnit             sampleUnits;
  UINT8ArraySequence *data;
}
UINT8TimeArraySeries;

typedef struct
tagREAL4TimeArraySeries
{
  CHAR                name[LALNameLength];
  LIGOTimeGPS         epoch;
  REAL8               deltaT;
  REAL8               f0;
  LALUnit             sampleUnits;
  REAL4ArraySequence *data;
}
REAL4TimeArraySeries;

typedef struct
tagREAL8TimeArraySeries
{
  CHAR                name[LALNameLength];
  LIGOTimeGPS         epoch;
  REAL8               deltaT;
  REAL8               f0;
  LALUnit             sampleUnits;
  REAL8ArraySequence *data;
}
REAL8TimeArraySeries;

typedef struct
tagCOMPLEX8TimeArraySeries
{
  CHAR                   name[LALNameLength];
  LIGOTimeGPS            epoch;
  REAL8                  deltaT;
  REAL8                  f0;
  LALUnit                sampleUnits;
  COMPLEX8ArraySequence *data;
}
COMPLEX8TimeArraySeries;

typedef struct
tagCOMPLEX16TimeArraySeries
{
  CHAR                    name[LALNameLength];
  LIGOTimeGPS             epoch;
  REAL8                   deltaT;
  REAL8                   f0;
  LALUnit                 sampleUnits;
  COMPLEX16ArraySequence *data;
}
COMPLEX16TimeArraySeries;


/* <lalLaTeX>

\vspace{2ex}
\begin{verbatim}
<datatype>FrequencySeries
\end{verbatim}
This structure represents a frequency spectrum of data of type
\verb@<datatype>@ (where \verb@<datatype>@ can be any primitive
datatype), sampled over uniform frequency intervals $f_0, f_0+\Delta
f, \ldots , f_0+l\Delta f$.  Essentially this is a
\verb@<datatype>Sequence@ with extra fields defining the sample
frequencies, the timestamp of the spectrum, and the type of data being
sampled.  The fields are:
\begin{description}
\item[\texttt{CHAR name[LALNameLength]}] The name of the data series (i.e.\
the type of data being sampled).
\item[\texttt{LIGOTimeGPS epoch}] The start time of the \emph{time}
series from which the spectrum was calculated.
\item[\texttt{REAL8 f0}] The lowest frequency $f_0$ being sampled, in
hertz.
\item[\texttt{REAL8 deltaF}] The frequency sampling interval $\Delta
f$, in hertz.
\item[\texttt{LALUnit sampleUnits}] The physical units of the
quantity being sampled.
\item[\texttt{<datatype>Sequence *data}] The sequence of sampled data.
\end{description}

</lalLaTeX> */

typedef struct
tagINT2FrequencySeries
{
  CHAR          name[LALNameLength];
  LIGOTimeGPS   epoch;
  REAL8         f0;
  REAL8         deltaF;
  LALUnit       sampleUnits;
  INT2Sequence *data;
}
INT2FrequencySeries;

typedef struct
tagUINT2FrequencySeries
{
  CHAR           name[LALNameLength];
  LIGOTimeGPS    epoch;
  REAL8          f0;
  REAL8          deltaF;
  LALUnit        sampleUnits;
  UINT2Sequence *data;
}
UINT2FrequencySeries;

typedef struct
tagINT4FrequencySeries
{
  CHAR          name[LALNameLength];
  LIGOTimeGPS   epoch;
  REAL8         f0;	
  REAL8         deltaF;
  LALUnit       sampleUnits;
  INT4Sequence *data;
}
INT4FrequencySeries;

typedef struct
tagUINT4FrequencySeries
{
  CHAR           name[LALNameLength];
  LIGOTimeGPS    epoch;
  REAL8          f0;	
  REAL8          deltaF;
  LALUnit        sampleUnits;
  UINT4Sequence *data;
}
UINT4FrequencySeries;

typedef struct
tagINT8FrequencySeries
{
  CHAR          name[LALNameLength];
  LIGOTimeGPS   epoch;
  REAL8         f0;
  REAL8         deltaF;
  LALUnit       sampleUnits;
  INT8Sequence *data;
}
INT8FrequencySeries;

typedef struct
tagUINT8FrequencySeries
{
  CHAR           name[LALNameLength];
  LIGOTimeGPS    epoch;
  REAL8          f0;
  REAL8          deltaF;
  LALUnit        sampleUnits;
  UINT8Sequence *data;
}
UINT8FrequencySeries;

typedef struct
tagREAL4FrequencySeries
{
  CHAR           name[LALNameLength];
  LIGOTimeGPS    epoch;
  REAL8          f0;	
  REAL8          deltaF;
  LALUnit        sampleUnits;
  REAL4Sequence *data;
}
REAL4FrequencySeries;

typedef struct
tagREAL8FrequencySeries
{
  CHAR           name[LALNameLength];
  LIGOTimeGPS    epoch;
  REAL8          f0;
  REAL8          deltaF;
  LALUnit        sampleUnits;
  REAL8Sequence *data;
}
REAL8FrequencySeries;

typedef struct
tagCOMPLEX8FrequencySeries
{
  CHAR              name[LALNameLength];
  LIGOTimeGPS       epoch;
  REAL8             f0;	
  REAL8             deltaF;
  LALUnit           sampleUnits;
  COMPLEX8Sequence *data;
}
COMPLEX8FrequencySeries;

typedef struct
tagCOMPLEX16FrequencySeries
{
  CHAR               name[LALNameLength];
  LIGOTimeGPS        epoch;
  REAL8              f0;	
  REAL8              deltaF;
  LALUnit            sampleUnits;
  COMPLEX16Sequence *data;
}
COMPLEX16FrequencySeries;

/* <lalLaTeX>

\vspace{2ex}
\begin{verbatim}
<datatype>ZPGFilter
\end{verbatim}
This structure stores the complex frequency response of a filter or
transfer function in a factored form, where \verb@<datatype>@ can be
either \verb@COMPLEX8@ or \verb@COMPLEX16@.  One defines a
(dimensionless) complex frequency variable $\zeta(f\Delta t)$, where
$\Delta t$ is the time sampling interval of the data to which the
filter will be applied (in the case of a digital filter), or some
other reference timescale (in the case of an analog filter).  The
complex response function can then be given (or approximated) as
$H(f)=g\times\prod_k(\zeta-z_k)/\prod_l(\zeta-p_l)$, where $z_k$ are
the complex \emph{zeros}, $p_l$ are the complex \emph{poles}, and $g$
is the complex \emph{gain} of the response function.  Some common
complex frequency representations are the $z$-plane representation
$\zeta(f\Delta t)=\exp(2\pi if\Delta t)$, which maps the Nyquist
interval $f\in[0,1/2\Delta t)$ onto the upper-half unit circle in
$\zeta$, and the $w$-plane representation $\zeta(f\Delta t)=\tan(\pi
f\Delta t)$, which maps the Nyquist interval onto the positive real
axis in $\zeta$.  The fields of \verb@<datatype>ZPGFilter@ are:
\begin{description}
\item[\texttt{CHAR name[LALNameLength]}] The name of the filter or transfer
function.  This should also mention its complex frequency
representation.
\item[\texttt{REAL8 deltaT}] The sampling time or reference timescale
$\Delta t$ for the filter, in seconds.  If zero, it will be treated as
being equal to the sampling interval of the data being filtered.
\item[\texttt{<datatype>Vector *zeros}]	Pointer to a vector storing
the zeros $z_k$ of the filter.
\item[\texttt{<datatype>Vector *poles}]	Pointer to a vector storing
the poles $p_k$ of the filter.
\item[\texttt{<datatype> gain}] The gain $g$ of the filter.
\end{description}

</lalLaTeX> */

typedef struct
tagCOMPLEX8ZPGFilter
{
  CHAR            name[LALNameLength];
  REAL8           deltaT;
  COMPLEX8Vector *zeros;
  COMPLEX8Vector *poles;
  COMPLEX8        gain;
}
COMPLEX8ZPGFilter;

typedef struct
tagCOMPLEX16ZPGFilter
{
  CHAR             name[LALNameLength];
  REAL8            deltaT;
  COMPLEX16Vector *zeros;
  COMPLEX16Vector *poles;
  COMPLEX16        gain;
}
COMPLEX16ZPGFilter;


/* <lalLaTeX>
\newpage
\subsection{The LAL universal status structure \texttt{LALStatus}}
\label{ss:status-structure}

This structure is the means by which LAL functions report their
success or failure; it provides a useful mechanism for tracking
progress and errors through nested function calls.  The error
reporting structure is a linked list of \verb@LALStatus@ structures, with
each node corresponding to a given function in the current calling
sequence.  When a function terminates successfully, its node is
dropped from the list.  If a function encounters an error, it must
still return control to the calling routine, reporting the error
through its \verb@LALStatus@.  The calling routine must either deal with
the error (pruning the linked list if it succeeds), or else return an
error itself.  A fatal error will thus return a linked list of
\verb@LALStatus@ structures to the top-level routine, where the tail of
the list identifies the source of the error, and the intermediate
nodes identify the sequence of nested function calls that led to the
error.  The fields of the \verb@LALStatus@ are as follows:
\begin{description}
\item[\texttt{INT4 statusCode}] A numerical code identifying the type
of error, or 0 for nominal status.
\item[\texttt{const CHAR *statusDescription}] A description of the
current status or error.
\item[\texttt{volatile const CHAR *Id}] The RCS ID string of the
source file of the current function.
\item[\texttt{const CHAR *function}] The name of the current function.
\item[\texttt{const CHAR *file}] The name of the source file of the
current function.
\item[\texttt{INT4 line}] The line number in the source file where the
current \verb@statusCode@ was set.
\item[\texttt{LALStatus *statusPtr}] Pointer to the next node in the
list; \verb@NULL@ if this function is not reporting a subroutine
error.
\item[\texttt{INT4 level}] The current level in the nested calling
sequence.
\end{description}

</lalLaTeX> */

typedef struct
tagLALStatus
{
  INT4                 statusCode;
  const CHAR          *statusDescription;
  volatile const CHAR *Id;
  const CHAR          *function;
  const CHAR          *file;
  INT4                 line;
  struct tagLALStatus *statusPtr;
  INT4                 level;
}
LALStatus;


/* <lalLaTeX>
\vfill{\footnotesize\input{LALDatatypesHV}}
</lalLaTeX> */


#ifdef  __cplusplus
}
#endif

#endif /* _LALDATATYPES_H */
