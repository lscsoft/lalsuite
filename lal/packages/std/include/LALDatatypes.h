/********************************* <lalVerbatim file="LALDatatypesHV">
Author: Finn, L. S.
$Id$
********************************** </lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{LALDatatypes.h}}
\label{s:LALDatatypes.h}

Provides the basic LAL datatypes.

\subsection*{Synopsis}
\begin{verbatim}
#include "LALDatatypes.h"
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
#include "LALAtomicDatatypes.h"
#include "LALRCSID.h"

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

</lalLaTeX> */

typedef struct
tagLIGOTimeGPS
{
  INT4 gpsSeconds;
  INT4 gpsNanoSeconds;
}
LIGOTimeGPS;


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
\item[\texttt{CHAR *name}] The name of the data series (i.e.\ the type
of data being sampled).
\item[\texttt{LIGOTimeGPS epoch}] The start time $t_0$ of the data
series.
\item[\texttt{REAL8 deltaT}] The sampling interval $\Delta t$, in
seconds.
\item[\texttt{REAL8 f0}] The heterodyning frequency $f_0$, in hertz.
\item[\texttt{CHARVector *sampleUnits}] The physical units of the
quantity being sampled.
\item[\texttt{<datatype>Sequence *data}] The sequence of sampled data.
\end{description}

</lalLaTeX> */

typedef struct
tagINT2TimeSeries
{
  CHAR         *name;
  LIGOTimeGPS   epoch;
  REAL8         deltaT;
  REAL8         f0;
  CHARVector   *sampleUnits;
  INT2Sequence *data;
}
INT2TimeSeries;

typedef struct
tagUINT2TimeSeries
{
  CHAR          *name;
  LIGOTimeGPS    epoch;
  REAL8          deltaT;
  REAL8          f0;
  CHARVector    *sampleUnits;
  UINT2Sequence *data;
}
UINT2TimeSeries;

typedef struct
tagINT4TimeSeries
{
  CHAR         *name;
  LIGOTimeGPS   epoch;
  REAL8         deltaT;
  REAL8         f0;
  CHARVector   *sampleUnits;
  INT4Sequence *data;
}
INT4TimeSeries;

typedef struct
tagUINT4TimeSeries
{
  CHAR          *name;
  LIGOTimeGPS    epoch;
  REAL8          deltaT;
  REAL8          f0;
  CHARVector    *sampleUnits;
  UINT4Sequence *data;
}
UINT4TimeSeries;

typedef struct
tagINT8TimeSeries
{
  CHAR         *name;
  LIGOTimeGPS   epoch;
  REAL8         deltaT;
  REAL8         f0;
  CHARVector   *sampleUnits;
  INT8Sequence *data;
}
INT8TimeSeries;

typedef struct
tagUINT8TimeSeries
{
  CHAR          *name;
  LIGOTimeGPS    epoch;
  REAL8          deltaT;
  REAL8          f0;
  CHARVector    *sampleUnits;
  UINT8Sequence *data;
}
UINT8TimeSeries;

typedef struct
tagREAL4TimeSeries
{
  CHAR          *name;
  LIGOTimeGPS    epoch;
  REAL8          deltaT;
  REAL8          f0;
  CHARVector    *sampleUnits;
  REAL4Sequence *data;
}
REAL4TimeSeries;

typedef struct
tagREAL8TimeSeries
{
  CHAR          *name;
  LIGOTimeGPS    epoch;
  REAL8          deltaT;
  REAL8          f0;
  CHARVector    *sampleUnits;
  REAL8Sequence *data;
}
REAL8TimeSeries;

typedef struct
tagCOMPLEX8TimeSeries
{
  CHAR             *name;
  LIGOTimeGPS       epoch;
  REAL8             deltaT;
  REAL8             f0;
  CHARVector       *sampleUnits;
  COMPLEX8Sequence *data;
}
COMPLEX8TimeSeries;

typedef struct
tagCOMPLEX16TimeSeries
{
  CHAR              *name;
  LIGOTimeGPS        epoch;
  REAL8              deltaT;
  REAL8              f0;
  CHARVector        *sampleUnits;
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
\item[\texttt{CHAR *name}] The name of the data series (i.e.\ the type
of data being sampled).
\item[\texttt{LIGOTimeGPS epoch}] The start time of the data series.
\item[\texttt{REAL8 deltaT}] The sampling interval, in seconds.
\item[\texttt{REAL8 f0}] The heterodyning frequency, in hertz.
\item[\texttt{CHARVector *sampleUnits}] The physical units of the
quantity being sampled.
\item[\texttt{<datatype>VectorSequence *data}] The sequence of sampled
data.
\end{description}

</lalLaTeX> */

typedef struct
tagINT2TimeVectorSeries
{
  CHAR                *name;
  LIGOTimeGPS          epoch;
  REAL8                deltaT;
  REAL8                f0;
  CHARVector          *sampleUnits;
  INT2VectorSequence  *data;
}
INT2TimeVectorSeries;

typedef struct
tagUINT2TimeVectorSeries
{
  CHAR                *name;
  LIGOTimeGPS          epoch;
  REAL8                deltaT;
  REAL8                f0;
  CHARVector          *sampleUnits;
  UINT2VectorSequence *data;
}
UINT2TimeVectorSeries;

typedef struct
tagINT4TimeVectorSeries
{
  CHAR                *name;
  LIGOTimeGPS          epoch;
  REAL8                deltaT;
  REAL8                f0;
  CHARVector          *sampleUnits;
  INT4VectorSequence  *data;
}
INT4TimeVectorSeries;

typedef struct
tagUINT4TimeVectorSeries
{
  CHAR                *name;
  LIGOTimeGPS          epoch;
  REAL8                deltaT;
  REAL8                f0;
  CHARVector          *sampleUnits;
  UINT4VectorSequence *data;
}
UINT4TimeVectorSeries;

typedef struct
tagINT8TimeVectorSeries
{
  CHAR                *name;
  LIGOTimeGPS          epoch;
  REAL8                deltaT;
  REAL8                f0;
  CHARVector          *sampleUnits;
  INT8VectorSequence  *data;
}
INT8TimeVectorSeries;

typedef struct
tagUINT8TimeVectorSeries
{
  CHAR                *name;
  LIGOTimeGPS          epoch;
  REAL8                deltaT;
  REAL8                f0;
  CHARVector          *sampleUnits;
  UINT8VectorSequence *data;
}
UINT8TimeVectorSeries;

typedef struct
tagREAL4TimeVectorSeries
{
  CHAR                *name;
  LIGOTimeGPS          epoch;
  REAL8                deltaT;
  REAL8                f0;
  CHARVector          *sampleUnits;
  REAL4VectorSequence *data;
}
REAL4TimeVectorSeries;

typedef struct
tagREAL8TimeVectorSeries
{
  CHAR                *name;
  LIGOTimeGPS          epoch;
  REAL8                deltaT;
  REAL8                f0;
  CHARVector          *sampleUnits;
  REAL8VectorSequence *data;
}
REAL8TimeVectorSeries;

typedef struct
tagCOMPLEX8TimeVectorSeries
{
  CHAR                    *name;
  LIGOTimeGPS              epoch;
  REAL8                    deltaT;
  REAL8                    f0;
  CHARVector              *sampleUnits;
  COMPLEX8VectorSequence  *data;
}
COMPLEX8TimeVectorSeries;

typedef struct
tagCOMPLEX16TimeVectorSeries
{
  CHAR                     *name;
  LIGOTimeGPS               epoch;
  REAL8                     deltaT;
  REAL8                     f0;
  CHARVector               *sampleUnits;
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
\item[\texttt{CHAR *name}] The name of the data series (i.e.\ the type
of data being sampled).
\item[\texttt{LIGOTimeGPS epoch}] The start time of the data series.
\item[\texttt{REAL8 deltaT}] The sampling interval, in seconds.
\item[\texttt{REAL8 f0}] The heterodyning frequency, in hertz.
\item[\texttt{CHARVector *sampleUnits}] The physical units of the
quantity being sampled.
\item[\texttt{<datatype>ArraySequence *data}] The sequence of sampled
data.
\end{description}

</lalLaTeX> */

typedef struct
tagINT2TimeArraySeries
{
  CHAR               *name;
  LIGOTimeGPS         epoch;
  REAL8               deltaT;
  REAL8               f0;
  CHARVector         *sampleUnits;
  INT2ArraySequence  *data;
}
INT2TimeArraySeries;

typedef struct
tagUINT2TimeArraySeries
{
  CHAR               *name;
  LIGOTimeGPS         epoch;
  REAL8               deltaT;
  REAL8               f0;
  CHARVector         *sampleUnits;
  UINT2ArraySequence *data;
}
UINT2TimeArraySeries;

typedef struct
tagINT4TimeArraySeries
{
  CHAR               *name;
  LIGOTimeGPS         epoch;
  REAL8               deltaT;
  REAL8               f0;
  CHARVector         *sampleUnits;
  INT4ArraySequence  *data;
}
INT4TimeArraySeries;

typedef struct
tagUINT4TimeArraySeries
{
  CHAR               *name;
  LIGOTimeGPS         epoch;
  REAL8               deltaT;
  REAL8               f0;
  CHARVector         *sampleUnits;
  UINT4ArraySequence *data;
}
UINT4TimeArraySeries;

typedef struct
tagINT8TimeArraySeries
{
  CHAR               *name;
  LIGOTimeGPS         epoch;
  REAL8               deltaT;
  REAL8               f0;
  CHARVector         *sampleUnits;
  INT8ArraySequence  *data;
}
INT8TimeArraySeries;

typedef struct
tagUINT8TimeArraySeries
{
  CHAR               *name;
  LIGOTimeGPS         epoch;
  REAL8               deltaT;
  REAL8               f0;
  CHARVector         *sampleUnits;
  UINT8ArraySequence *data;
}
UINT8TimeArraySeries;

typedef struct
tagREAL4TimeArraySeries
{
  CHAR               *name;
  LIGOTimeGPS         epoch;
  REAL8               deltaT;
  REAL8               f0;
  CHARVector         *sampleUnits;
  REAL4ArraySequence *data;
}
REAL4TimeArraySeries;

typedef struct
tagREAL8TimeArraySeries
{
  CHAR               *name;
  LIGOTimeGPS         epoch;
  REAL8               deltaT;
  REAL8               f0;
  CHARVector         *sampleUnits;
  REAL8ArraySequence *data;
}
REAL8TimeArraySeries;

typedef struct
tagCOMPLEX8TimeArraySeries
{
  CHAR                  *name;
  LIGOTimeGPS            epoch;
  REAL8                  deltaT;
  REAL8                  f0;
  CHARVector            *sampleUnits;
  COMPLEX8ArraySequence *data;
}
COMPLEX8TimeArraySeries;

typedef struct
tagCOMPLEX16TimeArraySeries
{
  CHAR                   *name;
  LIGOTimeGPS             epoch;
  REAL8                   deltaT;
  REAL8                   f0;
  CHARVector             *sampleUnits;
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
\item[\texttt{CHAR *name}] The name of the data series (i.e.\ the type
of data being sampled).
\item[\texttt{LIGOTimeGPS epoch}] The start time of the \emph{time}
series from which the spectrum was calculated.
\item[\texttt{REAL8 f0}] The lowest frequency $f_0$ being sampled, in
hertz.
\item[\texttt{REAL8 deltaF}] The frequency sampling interval $\Delta
f$, in hertz.
\item[\texttt{CHARVector *sampleUnits}] The physical units of the
quantity being sampled.
\item[\texttt{<datatype>Sequence *data}] The sequence of sampled data.
\end{description}

</lalLaTeX> */

typedef struct
tagINT2FrequencySeries
{
  CHAR         *name;
  LIGOTimeGPS   epoch;
  REAL8         f0;
  REAL8         deltaF;
  CHARVector   *sampleUnits;
  INT2Sequence *data;
}
INT2FrequencySeries;

typedef struct
tagUINT2FrequencySeries
{
  CHAR          *name;
  LIGOTimeGPS    epoch;
  REAL8          f0;
  REAL8          deltaF;
  CHARVector    *sampleUnits;
  UINT2Sequence *data;
}
UINT2FrequencySeries;

typedef struct
tagINT4FrequencySeries
{
  CHAR         *name;
  LIGOTimeGPS   epoch;
  REAL8         f0;	
  REAL8         deltaF;
  CHARVector   *sampleUnits;
  INT4Sequence *data;
}
INT4FrequencySeries;

typedef struct
tagUINT4FrequencySeries
{
  CHAR          *name;
  LIGOTimeGPS    epoch;
  REAL8          f0;	
  REAL8          deltaF;
  CHARVector    *sampleUnits;
  UINT4Sequence *data;
}
UINT4FrequencySeries;

typedef struct
tagINT8FrequencySeries
{
  CHAR         *name;
  LIGOTimeGPS   epoch;
  REAL8         f0;
  REAL8         deltaF;
  CHARVector   *sampleUnits;
  INT8Sequence *data;
}
INT8FrequencySeries;

typedef struct
tagUINT8FrequencySeries
{
  CHAR          *name;
  LIGOTimeGPS    epoch;
  REAL8          f0;
  REAL8          deltaF;
  CHARVector    *sampleUnits;
  UINT8Sequence *data;
}
UINT8FrequencySeries;

typedef struct
tagREAL4FrequencySeries
{
  CHAR          *name;
  LIGOTimeGPS    epoch;
  REAL8          f0;	
  REAL8          deltaF;
  CHARVector    *sampleUnits;
  REAL4Sequence *data;
}
REAL4FrequencySeries;

typedef struct
tagREAL8FrequencySeries
{
  CHAR          *name;
  LIGOTimeGPS    epoch;
  REAL8          f0;
  REAL8          deltaF;
  CHARVector    *sampleUnits;
  REAL8Sequence *data;
}
REAL8FrequencySeries;

typedef struct
tagCOMPLEX8FrequencySeries
{
  CHAR             *name;
  LIGOTimeGPS       epoch;
  REAL8             f0;	
  REAL8             deltaF;
  CHARVector       *sampleUnits;
  COMPLEX8Sequence *data;
}
COMPLEX8FrequencySeries;

typedef struct
tagCOMPLEX16FrequencySeries
{
  CHAR              *name;
  LIGOTimeGPS        epoch;
  REAL8              f0;	
  REAL8              deltaF;
  CHARVector        *sampleUnits;
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
\item[\texttt{CHAR *name}] The name of the filter or transfer
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
  CHAR           *name;
  REAL8           deltaT;
  COMPLEX8Vector *zeros;
  COMPLEX8Vector *poles;
  COMPLEX8        gain;
}
COMPLEX8ZPGFilter;

typedef struct
tagCOMPLEX16ZPGFilter
{
  CHAR            *name;
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
