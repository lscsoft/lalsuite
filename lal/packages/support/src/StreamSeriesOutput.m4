changecom(`/*',`*/')dnl
/*************************** <lalVerbatim file="StreamSeriesOutputCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{StreamSeriesOutput.c}}
\label{ss:StreamSeriesOutput.c}

Writes a time or frequency series to an output stream.

\subsubsection*{Prototypes}
\vspace{0.1in}
\begin{verbatim}
void
LAL<typecode>WriteTSeries( LALStatus            *stat,
                           FILE                 *stream,
                           <datatype>TimeSeries *series )

void
LAL<typecode>WriteTVectorSeries( LALStatus                  *stat,
                                 FILE                       *stream,
                                 <datatype>TimeVectorSeries *series )

void
LAL<typecode>WriteTArraySeries( LALStatus                 *stat,
                                FILE                      *stream,
                                <datatype>TimeArraySeries *series )

void
LAL<typecode>WriteFSeries( LALStatus                 *stat,
                           FILE                      *stream,
                           <datatype>FrequencySeries *series )
\end{verbatim}

\idx{LALI2WriteTSeries()}
\idx{LALI4WriteTSeries()}
\idx{LALI8WriteTSeries()}
\idx{LALU2WriteTSeries()}
\idx{LALU4WriteTSeries()}
\idx{LALU8WriteTSeries()}
\idx{LALSWriteTSeries()}
\idx{LALDWriteTSeries()}
\idx{LALCWriteTSeries()}
\idx{LALZWriteTSeries()}

\idx{LALI2WriteTVectorSeries()}
\idx{LALI4WriteTVectorSeries()}
\idx{LALI8WriteTVectorSeries()}
\idx{LALU2WriteTVectorSeries()}
\idx{LALU4WriteTVectorSeries()}
\idx{LALU8WriteTVectorSeries()}
\idx{LALSWriteTVectorSeries()}
\idx{LALDWriteTVectorSeries()}
\idx{LALCWriteTVectorSeries()}
\idx{LALZWriteTVectorSeries()}

\idx{LALI2WriteTArraySeries()}
\idx{LALI4WriteTArraySeries()}
\idx{LALI8WriteTArraySeries()}
\idx{LALU2WriteTArraySeries()}
\idx{LALU4WriteTArraySeries()}
\idx{LALU8WriteTArraySeries()}
\idx{LALSWriteTArraySeries()}
\idx{LALDWriteTArraySeries()}
\idx{LALCWriteTArraySeries()}
\idx{LALZWriteTArraySeries()}

\idx{LALI2WriteFSeries()}
\idx{LALI4WriteFSeries()}
\idx{LALI8WriteFSeries()}
\idx{LALU2WriteFSeries()}
\idx{LALU4WriteFSeries()}
\idx{LALU8WriteFSeries()}
\idx{LALSWriteFSeries()}
\idx{LALDWriteFSeries()}
\idx{LALCWriteFSeries()}
\idx{LALZWriteFSeries()}

\subsubsection*{Description}

These routines write the data and metadata in a time or frequency
series \verb@*series@ to an output stream \verb@*stream@ in a standard
format, described below.  It returns an error if any attempt to write
to the stream failed; \verb@*stream@ may then be left in a
partially-written state.

For each of these prototype templates there are in fact 10 separate
routines corresponding to all the atomic datatypes \verb@<datatype>@
(except \verb@CHAR@) referred to by \verb@<typecode>@:
\begin{center}
\begin{tabular}{|c@{\qquad}c|c@{\qquad}c|}
\hline
\tt <typecode> & \tt <datatype> & \tt <typecode> & \tt <datatype> \\
\hline
\tt I2 & \tt  INT2 & \tt U2 & \tt    UINT2  \\
\tt I4 & \tt  INT4 & \tt U4 & \tt    UINT4  \\
\tt I8 & \tt  INT8 & \tt U8 & \tt    UINT8  \\
\tt  S & \tt REAL4 & \tt  C & \tt COMPLEX8  \\
\tt  D & \tt REAL8 & \tt  Z & \tt COMPLEX16 \\
\hline
\end{tabular}
\end{center}

\paragraph{Format for \texttt{*stream}:} The data written to the
output stream will be formatted in a manner consistent with the input
routines in \verb@StreamSeriesInput.c@.  That is, it will begin with a
metadata header, consisting of multiple lines of the form:

\medskip
\begin{tabular}{l}
\verb@# @\textit{fieldname}\verb@ = @\textit{value}
\end{tabular}
\medskip

\noindent where \textit{fieldname} is the name of a field in
\verb@*series@ and \textit{value} is the value of that metadata field,
in some standard format (below).  The following metadata fields will
be written, one per line, based on the type of \verb@*series@:

\begin{description}
\item[\texttt{<datatype>TimeSeries}:] \verb@datatype@, \verb@name@,
\verb@epoch@, \verb@deltaT@, \verb@f0@, \verb@sampleUnits@,
\verb@length@
\item[\texttt{<datatype>TimeVectorSeries}:] \verb@datatype@,
\verb@name@, \verb@epoch@, \verb@deltaT@, \verb@f0@,
\verb@sampleUnits@, \verb@length@, \verb@vectorLength@
\item[\texttt{<datatype>TimeArraySeries}:] \verb@datatype@,
\verb@name@, \verb@epoch@, \verb@deltaT@, \verb@f0@,
\verb@sampleUnits@, \verb@length@, \verb@dimLength@, \verb@arrayDim@
\item[\texttt{<datatype>FrequencySeries}:] \verb@datatype@,
\verb@name@, \verb@epoch@, \verb@deltaT@, \verb@f0@, \verb@deltaF@,
\verb@sampleUnits@, \verb@length@
\end{description}

\noindent After all metadata have been written, the contents of
\verb@series->data->data@ will be written in standard integer or
floating-point notation, according to \verb@<datatype>@: integers will
be written to full precision, while floating-point numbers will be
written in exponential notation with sufficient digits to ensure that
they represent a unique binary floating-point number under the IEEE
Standard 754 (this means 9 digits for \verb@REAL4@s and 17 digits for
\verb@REAL8@s).  Complex datatypes are represented by pairs of
floating-point numbers representing alternately the real and imaginary
parts.

The body of the file will be formatted with newlines \verb@'\n'@
separating individual base, complex, vector, or array valued elements
of the sequence \verb@series->data@.  Within each element, integer or
floating-point components will be separated by single \verb@' '@
characters.  Thus the value of \verb@series->data->length@ will always
equal the number of lines following the metadata header.

\paragraph{Format for metadata fields:} Here we summarize briefly the
format for the individual field values in the metadata header.

\begin{description}
\item[\texttt{datatype}:] \textit{value} is a string (\emph{not}
surrounded by quotes) corresponding to the type of \verb@*series@;
e.g.\ \verb@COMPLEX8FrequencySeries@.

\item[\texttt{name}:] \textit{value} is a string surrounded by quotes
\verb@"@ representing \verb@series->name@.  Standard C-language string
literal notation is used: printable characters are written directly
except for \verb@"@ and \verb@\@ (rendered as \verb@\"@ and \verb@\\@,
respectively), characters with special C escape sequences are written
as those sequences (e.g.\ \verb@\t@ for tab and \verb@\n@ for
newline), and all other character bytes are written as three-digit
octal codes \verb@\@$ooo$.  Writing stops at the first null byte
\verb@\0@.

\item[\texttt{epoch}:] \textit{value} is a single \verb@INT8@ number
representing \verb@series->epoch@ in GPS nanoseconds.

\item[\texttt{deltaT}] (any time series): \textit{value} is a single
\verb@REAL8@ number representing \verb@series->deltaT@.

\item[\texttt{f0}:] \textit{value} is a single \verb@REAL8@ number
representing \verb@series->f0@.

\item[\texttt{deltaF}] (\verb@FrequencySeries@ only): \textit{value}
is a single \verb@REAL8@ number representing \verb@series->deltaF@.

\item[\texttt{sampleUnits}:] \textit{value} is string surrounded by
quotes \verb@"@; inside the quotes is a unit string corresponding to
\verb@series->sampleUnits@ as converted by the routine
\verb@LALUnitAsString()@.

\item[\texttt{length}:] \textit{value} is a single \verb@UINT4@
representing \verb@series->data->length@.

\item[\texttt{vectorLength}] (\verb@TimeVectorSeries@ only):
\textit{value} is a single \verb@UINT4@ representing
\verb@series->data->vectorLength@.

\item[\texttt{dimLength}] (\verb@TimeArraySeries@ only):
\textit{value} consists of a sequence of \verb@UINT4@s separated by
single \verb@' '@ characters, representing the components of
\verb@series->data->dimLength->data@.  The value of
\verb@series->data->dimLength->length@ must be inferred from the
number of components; it is not given as separate metadata.

\item[\texttt{arrayDim}] (\verb@TimeArraySeries@ only): \textit{value}
is a single \verb@UINT4@ representing \verb@series->data->arrayDim@.
If the array sequence was properly constructed, this will equal the
product of the components of \verb@dimLength@, above.
\end{description}

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel                           LALPrintError()
LALCHARCreateVector()                   LALCHARDestroyVector()
LALUnitAsString()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{StreamSeriesOutputCV}}

******************************************************* </lalLaTeX> */

#include <stdio.h>
#include <ctype.h>
#include <lal/LALStdlib.h>
#include <lal/Units.h>
#include <lal/AVFactories.h>
#include <lal/StringInput.h>
#include <lal/StreamOutput.h>
#include <lal/Date.h>

NRCSID( STREAMSERIESOUTPUTC, "$Id$" );

/* Define a function for printing a string as a literal. */
static int
LALWriteLiteral( FILE *stream, const CHAR *string )
{
  CHAR c; /* Current character being considered. */

  /* Open literal and start parsing string. */
  if ( putc( '"', stream ) == EOF ) return 1;
  while ( ( c = *(string++) ) != '\0' ) {

    /* Characters with named escape sequences. */
    if ( c == '\a' ) {
      if ( fputs( "\\a", stream ) == EOF ) return 1;
    } else if ( c == '\b' ) {
      if ( fputs( "\\b", stream ) == EOF ) return 1;
    } else if ( c == '\f' ) {
      if ( fputs( "\\f", stream ) == EOF ) return 1;
    } else if ( c == '\n' ) {
      if ( fputs( "\\n", stream ) == EOF ) return 1;
    } else if ( c == '\r' ) {
      if ( fputs( "\\r", stream ) == EOF ) return 1;
    } else if ( c == '\t' ) {
      if ( fputs( "\\t", stream ) == EOF ) return 1;
    } else if ( c == '\v' ) {
      if ( fputs( "\\v", stream ) == EOF ) return 1;
    } else if ( c == '"' ) {
      if ( fputs( "\\\"", stream ) == EOF ) return 1;
    } else if ( c == '\\' ) {
      if ( fputs( "\\\\", stream ) == EOF ) return 1;
    }

    /* Printable characters. */
    else if ( isprint( c ) ) {
      if ( putc( c, stream ) == EOF ) return 1;
    }

    /* Other characters are given 3-digit octal escape sequences. */
    else {
      if ( fprintf( stream, "\\%03o", (UCHAR)( c ) ) < 0 )
	return 1;
    }
  }

  /* Close literal and exit. */
  if ( putc( '"', stream ) == EOF ) return 1;
  return 0;
}

/* tell GNU C compiler to ignore warnings about the `ll' length modifier */
#ifdef __GNUC__
#define fprintf __extension__ fprintf
#endif

define(`TYPECODE',`I2')dnl
include(`LALWriteSeries.m4')dnl

define(`TYPECODE',`I4')dnl
include(`LALWriteSeries.m4')dnl

define(`TYPECODE',`I8')dnl
include(`LALWriteSeries.m4')dnl

define(`TYPECODE',`U2')dnl
include(`LALWriteSeries.m4')dnl

define(`TYPECODE',`U4')dnl
include(`LALWriteSeries.m4')dnl

define(`TYPECODE',`U8')dnl
include(`LALWriteSeries.m4')dnl

define(`TYPECODE',`S')dnl
include(`LALWriteSeries.m4')dnl

define(`TYPECODE',`D')dnl
include(`LALWriteSeries.m4')dnl

define(`TYPECODE',`C')dnl
include(`LALWriteSeries.m4')dnl

define(`TYPECODE',`Z')dnl
include(`LALWriteSeries.m4')dnl
