changecom(`/*',`*/')dnl
/***************************** <lalVerbatim file="StreamGridOutputCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{StreamGridOutput.c}}
\label{ss:StreamGridOutput.c}

Writes a LAL grid structure to an output stream.

\subsubsection*{Prototypes}
\vspace{0.1in}
\begin{verbatim}
void
LAL<typecode>WriteGrid( LALStatus *stat, FILE *stream, <datatype>Grid *grid )
\end{verbatim}

\idx{LALI2WriteGrid()}
\idx{LALI4WriteGrid()}
\idx{LALI8WriteGrid()}
\idx{LALU2WriteGrid()}
\idx{LALU4WriteGrid()}
\idx{LALU8WriteGrid()}
\idx{LALSWriteGrid()}
\idx{LALDWriteGrid()}
\idx{LALCWriteGrid()}
\idx{LALZWriteGrid()}

\subsubsection*{Description}

These routines write the data and metadata in a grid structure
\verb@*grid@ to an output stream \verb@*stream@ in a standard format,
described below.  It returns an error if any attempt to write to the
stream failed; \verb@*grid@ may then be left in a partially-written
state.

For each of these prototype templates there are in fact 10 separate
routines corresponding to all the numeric atomic datatypes
\verb@<datatype>@ referred to by \verb@<typecode>@:
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
routines in \verb@StreamGridInput.c@.  That is, it will begin with a
metadata header, consisting of multiple lines of the form:

\medskip
\begin{tabular}{l}
\verb@# @\textit{fieldname}\verb@ = @\textit{value}
\end{tabular}
\medskip

\noindent where \textit{fieldname} is the name of a field in
\verb@*series@ and \textit{value} is the value of that metadata field,
in some standard format (below).  The following metadata fields will
be written, one per line:

\begin{description}
\item[\texttt{datatype}:] \textit{value} is a string (\emph{not}
surrounded by quotes) corresponding to the type of \verb@*grid@;
e.g.\ \verb@COMPLEX8Grid@.

\item[\texttt{name}:] \textit{value} is a string surrounded by quotes
\verb@"@ representing \verb@grid->name@.  Standard C-language string
literal notation is used: printable characters are written directly
except for \verb@"@ and \verb@\@ (rendered as \verb@\"@ and \verb@\\@,
respectively), characters with special C escape sequences are written
as those sequences (e.g.\ \verb@\t@ for tab and \verb@\n@ for
newline), and all other character bytes are written as three-digit
octal codes \verb@\@$ooo$.  Writing stops at the first null byte
\verb@\0@.

\item[\texttt{sampleUnits}:] \textit{value} is string surrounded by
quotes \verb@"@; inside the quotes is a unit string corresponding to
\verb@grid->sampleUnits@ as converted by the routine
\verb@LALUnitAsString()@.

\item[\texttt{dimUnits}:] \textit{value} is a sequence of $m$ strings,
surrounded by quotes \verb@"@ and separated by a space, where $m$ is
the grid dimension (number of grid axes); inside the quotes is a unit
string corresponding to the elements of the \verb@grid->dimUnits@
array as converted by the routine \verb@LALUnitAsString()@.

\item[\texttt{offset}:] \textit{value} is a sequence of $m$
\verb@REAL8@ numbers separated by single spaces, representing the
elements of the \verb@grid->offset->data@; the number of data $m$ is
the grid dimension and corresponds to the value of
\verb@grid->offset->length@.

\item[\texttt{interval}:] \textit{value} is a sequence of $m$
\verb@REAL8@ numbers separated by single spaces, representing the
elements of the \verb@grid->interval->data@; the number of data $m$ is
the grid dimension and corresponds to the value of
\verb@grid->interval->length@.

\item[\texttt{dimLength}:] \textit{value} is a sequence of $M$
\verb@REAL8@ numbers separated by single spaces, representing the
elements of the \verb@grid->data->dimLength->data@; the number of data
$M$ is the data dimension and corresponds to the value of
\verb@grid->data->dimLength->length@, which must be greater than or
equal to the grid dimension $m$, above.
\end{description}

\noindent After all metadata have been written, the contents of
\verb@grid->data->data@ will be written in standard integer or
floating-point notation, according to \verb@<datatype>@: integers will
be written to full precision, while floating-point numbers will be
written in exponential notation with sufficient digits to ensure that
they represent a unique binary floating-point number under the IEEE
Standard 754 (this means 9 digits for \verb@REAL4@s and 17 digits for
\verb@REAL8@s).

The input format in \verb@StreamGridInput.c@ does not specify how the
numerical data is to be arranged, other than that the numbers be
separated by whitespace, and that complex datatypes be represented by
alternating real and imaginary components.  These routines adopt the
following conventions to improve human-readability: If the data
dimension is equal to the grid dimension, then each line consists of a
single datum (either a single number, or, for complex datatypes, a
pair of numbers separated by whitespace), followed by a newline
\verb@'\n'@.  If the data dimension is greater than the grid
dimension, then each line will consist of a number of data equal to
the length of the last dimension in \verb@grid->data->dimLength@.  If
the data dimension is at least two greater than the grid dimension,
and the dimension lengths are such that a single grid point comprises
multiple lines of data, then an additional blank line \verb@'\n'@ is
inserted to separate subsequent grid points.

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel                           LALPrintError()
LALCHARCreateVector()                   LALCHARDestroyVector()
LALUnitAsString()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{StreamGridOutputCV}}

% a " to fix C prettyprinting

******************************************************* </lalLaTeX> */

#include <stdio.h>
#include <ctype.h>
#include <lal/LALStdlib.h>
#include <lal/Grid.h>
#include <lal/Units.h>
#include <lal/AVFactories.h>
#include <lal/StringInput.h>
#include <lal/StreamOutput.h>

NRCSID( STREAMGRIDOUTPUTC, "$Id$" );

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

define(`TYPECODE',`I2')dnl
include(`LALWriteGrid.m4')dnl

define(`TYPECODE',`I4')dnl
include(`LALWriteGrid.m4')dnl

define(`TYPECODE',`I8')dnl
include(`LALWriteGrid.m4')dnl

define(`TYPECODE',`U2')dnl
include(`LALWriteGrid.m4')dnl

define(`TYPECODE',`U4')dnl
include(`LALWriteGrid.m4')dnl

define(`TYPECODE',`U8')dnl
include(`LALWriteGrid.m4')dnl

define(`TYPECODE',`S')dnl
include(`LALWriteGrid.m4')dnl

define(`TYPECODE',`D')dnl
include(`LALWriteGrid.m4')dnl

define(`TYPECODE',`C')dnl
include(`LALWriteGrid.m4')dnl

define(`TYPECODE',`Z')dnl
include(`LALWriteGrid.m4')dnl
