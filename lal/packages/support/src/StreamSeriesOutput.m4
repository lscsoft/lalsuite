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
LAL<typecode>WriteTimeSeries( LALStatus            *stat,
                              FILE                 *stream,
                              <datatype>TimeSeries *series )

void
LAL<typecode>WriteTimeVectorSeries( LALStatus                  *stat,
                                    FILE                       *stream,
                                    <datatype>TimeVectorSeries *series )

void
LAL<typecode>WriteTimeArraySeries( LALStatus                 *stat,
                                   FILE                      *stream,
                                   <datatype>TimeArraySeries *series )

void
LAL<typecode>WriteFrequencySeries( LALStatus                 *stat,
                                   FILE                      *stream,
                                   <datatype>FrequencySeries *series )
\end{verbatim}

\idx{LALI2WriteTimeSeries()}
\idx{LALI4WriteTimeSeries()}
\idx{LALI8WriteTimeSeries()}
\idx{LALU2WriteTimeSeries()}
\idx{LALU4WriteTimeSeries()}
\idx{LALU8WriteTimeSeries()}
\idx{LALSWriteTimeSeries()}
\idx{LALDWriteTimeSeries()}
\idx{LALCWriteTimeSeries()}
\idx{LALZWriteTimeSeries()}

\idx{LALI2WriteTimeVectorSeries()}
\idx{LALI4WriteTimeVectorSeries()}
\idx{LALI8WriteTimeVectorSeries()}
\idx{LALU2WriteTimeVectorSeries()}
\idx{LALU4WriteTimeVectorSeries()}
\idx{LALU8WriteTimeVectorSeries()}
\idx{LALSWriteTimeVectorSeries()}
\idx{LALDWriteTimeVectorSeries()}
\idx{LALCWriteTimeVectorSeries()}
\idx{LALZWriteTimeVectorSeries()}

\idx{LALI2WriteTimeArraySeries()}
\idx{LALI4WriteTimeArraySeries()}
\idx{LALI8WriteTimeArraySeries()}
\idx{LALU2WriteTimeArraySeries()}
\idx{LALU4WriteTimeArraySeries()}
\idx{LALU8WriteTimeArraySeries()}
\idx{LALSWriteTimeArraySeries()}
\idx{LALDWriteTimeArraySeries()}
\idx{LALCWriteTimeArraySeries()}
\idx{LALZWriteTimeArraySeries()}

\idx{LALI2WriteFrequencySeries()}
\idx{LALI4WriteFrequencySeries()}
\idx{LALI8WriteFrequencySeries()}
\idx{LALU2WriteFrequencySeries()}
\idx{LALU4WriteFrequencySeries()}
\idx{LALU8WriteFrequencySeries()}
\idx{LALSWriteFrequencySeries()}
\idx{LALDWriteFrequencySeries()}
\idx{LALCWriteFrequencySeries()}
\idx{LALZWriteFrequencySeries()}

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
\begin{tabular}{|c@{qquad}c|}
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

\begin{itemize}
\item[\texttt{<datatype>TimeSeries}:] \verb@name@, \verb@epoch@,
\verb@deltaT@, \verb@f0@, \verb@sampleUnits@, \verb@length@
\item[\texttt{<datatype>TimeVectorSeries}:] \verb@name@, \verb@epoch@,
\verb@deltaT@, \verb@f0@, \verb@sampleUnits@, \verb@length@,
\verb@vectorLength@
\item[\texttt{<datatype>TimeArraySeries}:] \verb@name@, \verb@epoch@,
\verb@deltaT@, \verb@f0@, \verb@sampleUnits@, \verb@length@,
\verb@dimLength@, \verb@arrayDim@
\item[\texttt{<datatype>FrequencySeries}:] \verb@name@, \verb@epoch@,
\verb@deltaT@, \verb@f0@, \verb@deltaF@, \verb@sampleUnits@,
\verb@length@
\end{itemize}

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

\begin{itemize}
\item[\texttt{name}:] \textit{value} is a string surrounded by quotes
\verb@"@ representing \verb@series->name@.  At present,
\verb@series->name@ will be truncated before any occurence of
\verb@'"'@, \verb@'\n'@, or \verb@'\0'@ characters.

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
\end{itemize}

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
LALCHARReadVector()                     LALCHARDestroyVector()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{StreamSeriesOutputCV}}

% This quote will fix the C syntax highlighting: "

******************************************************* </lalLaTeX> */

#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/Units.h>
#include <lal/AVFactories.h>
#include <lal/StringInput.h>
#include <lal/StreamOutput.h>

NRCSID( STREAMSERIESOUTPUTC, "$Id$" );

m4_define(`TYPECODE',`I2')m4_dnl
m4_include(`LALWriteSeries.m4')m4_dnl

m4_define(`TYPECODE',`I4')m4_dnl
m4_include(`LALWriteSeries.m4')m4_dnl

m4_define(`TYPECODE',`I8')m4_dnl
m4_include(`LALWriteSeries.m4')m4_dnl

m4_define(`TYPECODE',`U2')m4_dnl
m4_include(`LALWriteSeries.m4')m4_dnl

m4_define(`TYPECODE',`U4')m4_dnl
m4_include(`LALWriteSeries.m4')m4_dnl

m4_define(`TYPECODE',`U8')m4_dnl
m4_include(`LALWriteSeries.m4')m4_dnl

m4_define(`TYPECODE',`S')m4_dnl
m4_include(`LALWriteSeries.m4')m4_dnl

m4_define(`TYPECODE',`D')m4_dnl
m4_include(`LALWriteSeries.m4')m4_dnl

m4_define(`TYPECODE',`C')m4_dnl
m4_include(`LALWriteSeries.m4')m4_dnl

m4_define(`TYPECODE',`Z')m4_dnl
m4_include(`LALWriteSeries.m4')m4_dnl
