/**************************** <lalVerbatim file="StreamSeriesInputCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{StreamSeriesInput.c}}
\label{ss:StreamSeriesInput.c}

Converts an input stream into a time or frequency series.

\subsubsection*{Prototypes}
\vspace{0.1in}
\begin{verbatim}
void
LAL<typecode>ReadTimeSeries( LALStatus            *stat,
                             <datatype>TimeSeries *series,
                             FILE                 *stream )

void
LAL<typecode>ReadTimeVectorSeries( LALStatus                  *stat,
                                   <datatype>TimeVectorSeries *series,
                                   FILE                       *stream )

void
LAL<typecode>ReadTimeArraySeries( LALStatus                 *stat,
                                  <datatype>TimeArraySeries *series,
                                  FILE                      *stream )

void
LAL<typecode>ReadFrequencySeries( LALStatus                 *stat,
                                  <datatype>FrequencySeries *series,
                                  FILE                      *stream )
\end{verbatim}

\idx{LALI2ReadTimeSeries()}
\idx{LALI4ReadTimeSeries()}
\idx{LALI8ReadTimeSeries()}
\idx{LALU2ReadTimeSeries()}
\idx{LALU4ReadTimeSeries()}
\idx{LALU8ReadTimeSeries()}
\idx{LALSReadTimeSeries()}
\idx{LALDReadTimeSeries()}
\idx{LALCReadTimeSeries()}
\idx{LALZReadTimeSeries()}

\idx{LALI2ReadTimeVectorSeries()}
\idx{LALI4ReadTimeVectorSeries()}
\idx{LALI8ReadTimeVectorSeries()}
\idx{LALU2ReadTimeVectorSeries()}
\idx{LALU4ReadTimeVectorSeries()}
\idx{LALU8ReadTimeVectorSeries()}
\idx{LALSReadTimeVectorSeries()}
\idx{LALDReadTimeVectorSeries()}
\idx{LALCReadTimeVectorSeries()}
\idx{LALZReadTimeVectorSeries()}

\idx{LALI2ReadTimeArraySeries()}
\idx{LALI4ReadTimeArraySeries()}
\idx{LALI8ReadTimeArraySeries()}
\idx{LALU2ReadTimeArraySeries()}
\idx{LALU4ReadTimeArraySeries()}
\idx{LALU8ReadTimeArraySeries()}
\idx{LALSReadTimeArraySeries()}
\idx{LALDReadTimeArraySeries()}
\idx{LALCReadTimeArraySeries()}
\idx{LALZReadTimeArraySeries()}

\idx{LALI2ReadFrequencySeries()}
\idx{LALI4ReadFrequencySeries()}
\idx{LALI8ReadFrequencySeries()}
\idx{LALU2ReadFrequencySeries()}
\idx{LALU4ReadFrequencySeries()}
\idx{LALU8ReadFrequencySeries()}
\idx{LALSReadFrequencySeries()}
\idx{LALDReadFrequencySeries()}
\idx{LALCReadFrequencySeries()}
\idx{LALZReadFrequencySeries()}

\subsubsection*{Description}

These routines parse an input stream \verb@*stream@ to fill in the
data and metadata fields of a time or frequency series \verb@*series@.
The field \verb@series->data@ must be \verb@NULL@, so that it can be
created and filled by the routine.  The other fields may be
initialized or not; they will be overwritten by metadata read from
\verb@*stream@.  If an error occurs, \verb@*series@ will be left
unchanged, but \verb@*stream@ will have been read up to the point
where the error occured.

For each of these prototype templates there are in fact 10 separate
routines corresponding to all the atomic datatypes \verb@<datatype>@
(except \verb@CHAR@) referred to by \verb@<typecode>@:
\begin{center}
\begin{tabular}{|c@{\qquad}c|}
\hline
\tt <typecode> & \tt <datatype> \\
\hline
\tt I2 & \tt     INT2  \\
\tt I4 & \tt     INT4  \\
\tt I8 & \tt     INT8  \\
\tt U2 & \tt    UINT2  \\
\tt U4 & \tt    UINT4  \\
\tt U8 & \tt    UINT8  \\
\tt  S & \tt    REAL4  \\
\tt  D & \tt    REAL8  \\
\tt  C & \tt COMPLEX8  \\
\tt  Z & \tt COMPLEX16 \\
\hline
\end{tabular}
\end{center}

\paragraph{Format for \texttt{*stream}:} The input stream is assumed
to be a text stream (ASCII) consisting of a header containing metadata
followed by numerical data in standard integer or floating-point
format, as recognized by the routines in \verb@StringConvert.c@.  The
header consists of zero or more lines beginning with a \verb@'#'@
character, followed by a metadata field name and value in the format:

\medskip
\begin{tabular}{l}
\verb@# @\textit{fieldname}\verb@=@\textit{value}
\end{tabular}
\medskip

\noindent The \verb@=@ sign in this format is standard but optional;
it may be replaced or surrounded with any amount of any whitespace
except a newline \verb@'\n'@.  If \textit{fieldname} is unrecognized,
it is ignored; if it is recognized, then \textit{value} must be in a
suitable format for the field type, as described below.  Blank lines,
or lines containing just a \verb@#@ character, are skipped.  Once a
line is encountered that contains non-whitespace characters and does
not start with \verb@'#'@, that line is assumed to be the beginning of
the numerical data.  From that point on, all non-whitespace characters
must be part of parseable numbers; no more comments are permitted
(although blank lines will still be skipped).

If a metadata field appears twice in the header, the later one takes
precedence.  At present these routines do not track which fields have
been previously assigned, so no warnings or errors are generated.

How the data is packed into the \verb@series->data@ structure depends
on what metadata has been provided, as described below.

\paragraph{Required, conditional, and optional metadata:} The input
stream need not contain a complete set of metadata, allowing some
metadata to be read from \verb@*stream@ and others to be set
elsewhere.  For each type of series, some metadata will be
\emph{required}, and the routine will abort if the metadata is not
found.  Other metadata are \emph{conditional}, meaning that the
routine will operate differently depending on whether or not these
metadata were found.  The remaining metadata are \emph{optional}; if
they are not found in \verb@*stream@, they will be left unchanged.
The recognized metadata fields are listed below.

\medskip\noindent\verb@<datatype>TimeSeries@:
\begin{description}
\item[Required fields:] none
\item[Conditional fields:] \verb@length@
\item[Optional fields:] \verb@name@, \verb@epoch@, \verb@deltaT@,
\verb@f0@, \verb@sampleUnits@
\end{description}

\medskip\noindent\verb@<datatype>TimeVectorSeries@:
\begin{description}
\item[Required fields:] none
\item[Conditional fields:] \verb@length@, \verb@vectorLength@
\item[Optional fields:] \verb@name@, \verb@epoch@, \verb@deltaT@,
\verb@f0@, \verb@sampleUnits@
\end{description}

\medskip\noindent\verb@<datatype>TimeArraySeries@:
\begin{description}
\item[Required fields:] \verb@dimLength@
\item[Conditional fields:] \verb@length@, \verb@arrayDim@
\item[Optional fields:] \verb@name@, \verb@epoch@, \verb@deltaT@,
\verb@f0@, \verb@sampleUnits@
\end{description}

\medskip\noindent\verb@<datatype>FrequencySeries@:
\begin{description}
\item[Required fields:] none
\item[Conditional fields:] \verb@length@
\item[Optional fields:] \verb@name@, \verb@epoch@, \verb@deltaT@,
\verb@f0@, \verb@deltaF@, \verb@sampleUnits@
\end{description}

Below we describe the required format for the field values, as well as
what occurs if a conditional field is or isn't present.

\subparagraph{Required fields:}
\begin{description}
\item[\texttt{dimLength}] (\verb@TimeArraySeries@ only):
\textit{value} consists of a sequence of \verb@UINT4@s separated by
whitespace (but \emph{not} a newline \verb@'\n'@).  These data are
stored in \verb@series->data->dimLength@: the number of integers gives
the number of array indecies, while the value of each integer gives
the dimension of the corresponding array index.
\end{description}

\subparagraph{Conditional fields:}
\begin{description}
\item[\texttt{arrayDim}] (\verb@TimeArraySeries@ only): \textit{value}
is a single \verb@UINT4@, to be stored in
\verb@series->data->arrayDim@.  This must equal the product of the
index ranges in \verb@dimLength@, above, or an error is returned.  If
not given, the \verb@arrayDim@ field will be set equal to the product
of the index ranges in \verb@dimLength@.  (The \verb@arrayDim@ and
\verb@dimLength@ fields can appear in any order in \verb@*stream@;
checking is done only after all header lines have been read.)

\item[\texttt{vectorLength}] (\verb@TimeVectorSeries@ only):
\textit{value} is a single \verb@UINT4@, to be stored in
\verb@series->data->vectorLength@.  If not specified in the header
portion of \verb@*stream@, it will be taken to be the number of data
on the \emph{first} line of the data portion of \verb@*stream@, or
half the number of real data for a complex-valued
\verb@TimeVectorSeries@; if an odd number of real data are found on
the first line of a complex \verb@TimeVectorSeries@, then an error is
returned.

\item[\texttt{length}:] \textit{value} is a single \verb@UINT4@, to be
stored in \verb@series->data->length@.  If it is specified in the
header portion of \verb@*stream@, data will be read until
\verb@length@ is reached.  Otherwise, \verb@*stream@ will be read to
its end or until an unparseable character is read, and \verb@length@
will then be set accordingly.  (If parsing stops in the middle of
filling a complex, vector, or array valued element, the partly-read
element is discarded.)
\end{description}

\subparagraph{Optional fields:}
\begin{description}
\item[\texttt{name}:] \textit{value} is a string surrounded by quotes
\verb@"@, of length less than \verb@LALNameLength@.  The contents and
a terminating \verb@'\0'@ character are stored in \verb@series->name@.
At present, nested quotes \verb@"@ and newlines are not permitted.

\item[\texttt{epoch}:] \textit{value} is a single \verb@INT8@ number
of GPS nanoseconds, or a pair of \verb@INT4@s representing GPS seconds
and nanoseconds separately, separated by non-newline whitespace.

\item[\texttt{deltaT}] (any time series): \textit{value} is a single
\verb@REAL8@ number.

\item[\texttt{f0}:] \textit{value} is a single \verb@REAL8@ number.

\item[\texttt{deltaF}] (\verb@FrequencySeries@ only): \textit{value}
is a single \verb@REAL8@ number.

\item[\texttt{sampleUnits}:] \textit{value} is string surrounded by
quotes \verb@"@; the quotes are stripped and the string passed to
\verb@LALParseUnitString()@ to determine \verb@series->sampleUnits@.
Since \verb@LALParseUnitString()@ is not very robust, it is
recommended to use only unit strings that have been generated by
\verb@LALUnitAsString()@, or to remove this metadata field and set
\verb@series->sampleUnits@ within the code.
\end{description}

\paragraph{Data format:} The data portion of \verb@*stream@ consists
of whitespace-separated integer or real numbers.  For complex input
routines, the real data are parsed as alternately the real and
imaginary parts of successive complex numbers.  By convention, each
line should correspond to a single base, complex, vector, or array
valued element of the \verb@series->data@ sequence.  However, this is
\emph{required} only in the case of a \verb@TimeVectorSeries@ where
the \verb@vectorLength@ metadata was not set in the header, since in
this case the value of \verb@vectorLength@ will be taken from the
number of elements read on the first data line.  After this, and in
all other cases, newlines are treated as any other whitespace.

If a \verb@length@ value is specified in the header, then data are
read until the required length is acheived; if \verb@fscanf()@ returns
zero or negative before this (representing either the end-of-input or
a character that cannot be interpreted as part of the numerical data),
an error is returned.  If a \verb@length@ value was not specified,
data are read until \verb@fscanf()@ returns zero or negative: at this
point any partially-completed complex, vector, or array valued element
is discarded, and \verb@series->data->length@ set to the number of
elements read.

\subsubsection*{Algorithm}

These routines use \verb@LALCHARReadVector()@ to read the header lines
and the first line of data.  After this, data are parsed directly from
\verb@*stream@ using \verb@fscanf()@.  This is done for efficiency:
repeated calling of the LAL string parsing routines in
\verb@StringConvert.c@ involves far too much computational overhead.

After the first data line has been read, the length of each sequence
element will be known from the atomic type, as well as the specified
\verb@dimLength@ (for arrays), \verb@vectorLength@ (for vectors), or
number of elements on the first data line (for vectors without an
explicitly specified \verb@vectorLength@).  If \verb@length@ is also
specified, a sequence of the appropriate size is allocated, and all
the data is copied or read directly into it.  If \verb@length@ was not
specified, the data read with \verb@fscanf()@ are stored in a linked
list of buffers of size \verb@BUFFSIZE@ (a local \verb@#define@d
constant) until parsing stops.  Then a sequence of the appropriate
size is allocated and the data copied into it.

Although the \verb@CHAR@, \verb@UCHAR@, and \verb@BOOLEAN@ series
datatypes exist in LAL, they are rarely used, and no input routines
are provided for them.  This may change in future if someone expresses
a need for them.

\subsubsection*{Uses}
\begin{verbatim}
LALMalloc()                             LALFree()
LALCHARReadVector()                     LALCHARDestroyVector()
LAL<typecode>CreateVector()             LAL<typecode>DestroyVector()
LAL<typecode>CreateVectorSequence()     LAL<typecode>DestroyVectorSequence()
LAL<typecode>CreateArraySequence()      LAL<typecode>DestroyArraySequence()
LALStringTo<typecode>()                 LALWarning()
\end{verbatim}
where \verb@<typecode>@ is any of \verb@I2@, \verb@I4@, \verb@I8@,
\verb@U2@, \verb@U4@, \verb@U8@, \verb@S@, \verb@D@, \verb@C@, or
\verb@Z@.

\subsubsection*{Notes}

\vfill{\footnotesize\input{StreamSeriesInputCV}}

% This quote will fix the C syntax highlighting: "

******************************************************* </lalLaTeX> */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/Units.h>
#include <lal/StringInput.h>
#include <lal/StreamInput.h>

NRCSID( STREAMSERIESINPUTC, "$Id$" );

/* Define a message string for header parsing errors. */
#define LALREADSERIESC_HEADER "Skipping badly-formatted line for metadata field "

/* Define linked-list of buffers for storing an arbitrary number of
arbitrary datatypes. */
#define BUFFSIZE 24
typedef union tagBuffer {
  INT2 I2[BUFFSIZE/2];
  INT4 I4[BUFFSIZE/4];
  INT8 I8[BUFFSIZE/8];
  UINT2 U2[BUFFSIZE/2];
  UINT4 U4[BUFFSIZE/4];
  UINT8 U8[BUFFSIZE/8];
  REAL4 S[BUFFSIZE/4];
  REAL8 D[BUFFSIZE/8];
} Buffer;
typedef struct tagBufferList {
  Buffer buf;
  struct tagBufferList *next;
} BufferList;

/* Define a macro for freeing the linked list. */
#define FREEBUFFERLIST( headPtr )                                    \
if ( headPtr ) {                                                     \
  BufferList *herePtr = headPtr;                                     \
  while ( herePtr ) {                                                \
    BufferList *nextPtr = herePtr->next;                             \
    LALFree( herePtr );                                              \
    herePtr = nextPtr;                                               \
  }                                                                  \
} else (void)(0)

static const BufferList empty;

m4_define(`TYPECODE',`I2')m4_dnl
m4_include(`LALReadSeries.m4')m4_dnl

m4_define(`TYPECODE',`I4')m4_dnl
m4_include(`LALReadSeries.m4')m4_dnl

m4_define(`TYPECODE',`I8')m4_dnl
m4_include(`LALReadSeries.m4')m4_dnl

m4_define(`TYPECODE',`U2')m4_dnl
m4_include(`LALReadSeries.m4')m4_dnl

m4_define(`TYPECODE',`U4')m4_dnl
m4_include(`LALReadSeries.m4')m4_dnl

m4_define(`TYPECODE',`U8')m4_dnl
m4_include(`LALReadSeries.m4')m4_dnl

m4_define(`TYPECODE',`S')m4_dnl
m4_include(`LALReadSeries.m4')m4_dnl

m4_define(`TYPECODE',`D')m4_dnl
m4_include(`LALReadSeries.m4')m4_dnl

m4_define(`TYPECODE',`C')m4_dnl
m4_include(`LALReadSeries.m4')m4_dnl

m4_define(`TYPECODE',`Z')m4_dnl
m4_include(`LALReadSeries.m4')m4_dnl
