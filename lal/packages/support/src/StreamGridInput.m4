changecom(`/*',`*/')dnl
/**************************** <lalVerbatim file="StreamGridInputCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{StreamGridInput.c}}
\label{ss:StreamGridInput.c}

Converts an input stream into a LAL grid structure.

\subsubsection*{Prototypes}
\vspace{0.1in}
\begin{verbatim}
void
LAL<typecode>ReadGrid( LALStatus *stat, <datatype>Grid **grid, FILE *stream )
\end{verbatim}

\idx{LALI2ReadGrid()}
\idx{LALI4ReadGrid()}
\idx{LALI8ReadGrid()}
\idx{LALU2ReadGrid()}
\idx{LALU4ReadGrid()}
\idx{LALU8ReadGrid()}
\idx{LALSReadGrid()}
\idx{LALDReadGrid()}
\idx{LALCReadGrid()}
\idx{LALZReadGrid()}

\subsubsection*{Description}

These routines parse an input stream \verb@*stream@ to create a grid
structure \verb@**grid@, and fill in its data and metadata fields.
The output parameter \verb@grid@ must be a non-\verb@NULL@ handle to a
\verb@NULL@-valued pointer \verb@*grid@.

This prototype template in fact refers to 10 separate routines
corresponding to all the numerical atomic datatypes \verb@<datatype>@
referred to by \verb@<typecode>@:
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
Some metadata fields are \emph{required}, and the routine will abort
if it doesn't find them.  Others are \emph{optional}: if they are not
found, the routine will assign some default value.  The various fields
and their required formats are given below:

\bigskip\noindent\textit{Required fields:}
\begin{description}
\item[\texttt{dimLength}:] \textit{value} consists of a sequence of
\verb@UINT4@s separated by whitespace (but \emph{not} a newline
\verb@'\n'@).  These are used to create \verb@(*grid)->data@ with the
appropriate dimensions, and are stored in
\verb@(*grid)->data->dimLength->data@: the number of integers $M$
gives the data dimension number (the number of array indecies), while
the value of each integer gives the length of each dimension (the
range of values of the corresponding array index).

\item[\texttt{offset}:] \textit{value} consists of a sequence of
\verb@REAL8@s separated by whitespace (but \emph{not} a newline
\verb@'\n'@).  These values are stored in
\verb@(*grid)->offset->data@.  The number of data $m$ gives the grid
dimension, which must be less than or equal to the data dimension $M$
of the array defined above, and must be consistent among the
\verb@offset@, \verb@interval@, and \verb@dimUnits@ fields.

\item[\texttt{interval}:] \textit{value} consists of a sequence of
\verb@REAL8@s separated by whitespace (but \emph{not} a newline
\verb@'\n'@).  These values are stored in
\verb@(*grid)->interval->data@.  The number of data $m$ gives the grid
dimension, which must be less than or equal to the data dimension $M$
of the array defined above, and must be consistent among the
\verb@offset@, \verb@interval@, and \verb@dimUnits@ fields.
\end{description}


\medskip\noindent\textit{Optional fields:}
\begin{description}
\item[\texttt{name}:] \textit{value} is a string surrounded by quotes
\verb@"@, which is parsed in the manner of a string literal in C: it
may contain ordinary printable characters (except \verb@"@ and
\verb@\@), escape sequences (such as \verb@\t@ for tab, \verb@\n@ for
newline, or \verb@\\@ and \verb@\"@ for literal backslash and quote
characters), and octal or hexadecimal codes (\verb@\@$ooo$ or
\verb@\x@$hh$ respectively) for arbitrary bytes.  Unlike in C,
literals cannot be split between lines, adjacent literals are not
concatenated, and converted strings longer than
\verb@LALNameLength@$-1$ will be truncated.  The resulting string is
stored in \verb@(*grid)->name@, and will always contain a \verb@\0@
terminator, beyond which the contents are unspecified.  If this field
is not given in \verb@stream@, then the routine will simply assign
\verb@(*grid)->name[0]@=\verb@'\0'@.

\item[\texttt{sampleUnits}:] \textit{value} is string surrounded by
quotes \verb@"@; the quotes are stripped and the string passed to
\verb@LALParseUnitString()@ to determine \verb@(*grid)->sampleUnits@.
Since \verb@LALParseUnitString()@ is not very robust, it is
recommended to use only unit strings that have been generated by
\verb@LALUnitAsString()@, or to remove this metadata field and set
\verb@(*grid)->sampleUnits@ within the code.  If this field is not
given in \verb@stream@, then \verb@lalDimensionlessUnit@ is assumed.

\item[\texttt{dimUnits}:] \textit{value} is a sequence of strings,
each surrounded by quotes \verb@"@, and optionally separated by
whitespace (but \emph{not} a newline \verb@'\n'@); the quotes are
stripped and the string passed to \verb@LALParseUnitString()@ to
determine \verb@(*grid)->dimUnits@.  Since \verb@LALParseUnitString()@
is not very robust, it is recommended to use only unit strings that
have been generated by \verb@LALUnitAsString()@, or to remove this
metadata field and reset \verb@(*grid)->dimUnits@ within the code.  If
this field is not given in \verb@stream@, then
\verb@(*grid)->dimUnits@ will be allocated as an array containing $m$
units assigned the value \verb@lalDimensionlessUnit@, where $m$ is the
length of the \verb@interval@ and \verb@offset@ vectors.  If this
field \emph{is} given, then the number of unit strings must be
consistent with the lengths of the \verb@interval@ and \verb@offset@
vectors, or the routine will abort.

\item[\texttt{datatype}:] \textit{value} is string identifying the
grid type; e.g. \verb@REAL4Grid@ (\emph{not} surrounded by quotes).
This should correspond to the type of \verb@**grid@, not to any field
in \verb@**grid@.  If there is a type mismatch, a warning is generated
(and errors may occur later while parsing the data).
\end{description}

\paragraph{Data format:} The first line that is neither blank nor
beginning with a \verb@'#'@ character is assumed to be the start of
the grid data, and is parsed as a sequence of whitespace-separated
integers or real numbers.  For complex datatypes, the numbers read are
interpreted as alternately the real and imaginary parts of the data.
Since the number of data required is known from the grid metadata, the
routine will read from \verb@stream@ only until the structure is
filled (except that the first line of data will be read in its
entirety, even if it contains more numbers than required).  If the
routine encounters an unparseable number or the end-of-input before
filling the structure, an error is returned.

\subsubsection*{Algorithm}

These routines use \verb@LALCHARReadVector()@ to read the header lines
and the first line of data.  The metadata are stored in temporary
variables or vectors.  After the first data line has been read, the
size of the grid should be known from the metadata: the structure is
created, and the metadata and first line of data are copied in.  After
this, any additional data are parsed directly from \verb@stream@ using
\verb@fscanf()@.  This is done for efficiency: repeated calling of the
LAL string parsing routines in \verb@StringConvert.c@ involves far too
much computational overhead.

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALPrintError()                         LALWarning()
LALMalloc()                             LALFree()
LALCHARReadVector()                     LALCHARDestroyVector()
LALDCreateVector()                      LALDDestroyVector()
LAL<typecode>CreateGrid()               LAL<typecode>DestroyGrid()
LALStringTo<typecode>()                 LALParseUnitString()
\end{verbatim}
where \verb@<typecode>@ is any of \verb@I2@, \verb@I4@, \verb@I8@,
\verb@U2@, \verb@U4@, \verb@U8@, \verb@S@, \verb@D@, \verb@C@, or
\verb@Z@.

\subsubsection*{Notes}

\vfill{\footnotesize\input{StreamGridInputCV}}

% This quote will fix the C syntax highlighting: "

******************************************************* </lalLaTeX> */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/Units.h>
#include <lal/Grid.h>
#include <lal/StringInput.h>
#include <lal/StreamInput.h>

NRCSID( STREAMGRIDINPUTC, "$Id$" );

/* Define a message string for header parsing errors. */
#define LALREADGRIDC_HEADER "Skipping badly-formatted line for metadata field "

/* Define linked-list of buffers for storing an arbitrary number of
UINT4s, REAL8s, or LALUnits. */
#define BUFFSIZE 24
typedef struct tagU4Buffer {
  UINT4 U4[BUFFSIZE];
  struct tagU4Buffer *next;
} U4Buffer;
typedef struct tagDBuffer {
  REAL8 D[BUFFSIZE];
  struct tagDBuffer *next;
} DBuffer;
typedef struct tagUnitBuffer {
  LALUnit Unit[BUFFSIZE];
  struct tagUnitBuffer *next;
} UnitBuffer;

/* Define macros for freeing the linked lists. */
#define FREEU4BUFFER( headPtr )                                      \
if ( headPtr ) {                                                     \
  U4Buffer *herePtr = headPtr;                                       \
  while ( herePtr ) {                                                \
    U4Buffer *nextPtr = herePtr->next;                               \
    LALFree( herePtr );                                              \
    herePtr = nextPtr;                                               \
  }                                                                  \
} else (void)(0)

#define FREEDBUFFER( headPtr )                                       \
if ( headPtr ) {                                                     \
  DBuffer *herePtr = headPtr;                                        \
  while ( herePtr ) {                                                \
    DBuffer *nextPtr = herePtr->next;                                \
    LALFree( herePtr );                                              \
    herePtr = nextPtr;                                               \
  }                                                                  \
} else (void)(0)

#define FREEUNITBUFFER( headPtr )                                    \
if ( headPtr ) {                                                     \
  UnitBuffer *herePtr = headPtr;                                     \
  while ( herePtr ) {                                                \
    UnitBuffer *nextPtr = herePtr->next;                             \
    LALFree( herePtr );                                              \
    herePtr = nextPtr;                                               \
  }                                                                  \
} else (void)(0)

/* Define a macro for freeing temporary metadata. */
#define CLEANUP                                                      \
do {                                                                 \
  TRY( LALCHARDestroyVector( stat->statusPtr, &line ), stat );       \
  if ( dimLength )                                                   \
    LALFree( dimLength );                                            \
  if ( dimUnits )                                                    \
    LALFree( dimUnits );                                             \
  if ( offset )                                                      \
    LALFree( offset );                                               \
  if ( interval )                                                    \
    LALFree( interval );                                             \
} while (0)


/* Define a function for parsing a string literal. */
static void
LALLiteralToString( LALStatus  *stat,
		    CHAR       *string,
		    const CHAR *literal,
		    UINT4      length )
{
  CHAR c;       /* Current character being considered. */
  UINT4 n = 0;  /* Counter of number of characters written. */

  INITSTATUS( stat, "LALLiteralToString", STREAMGRIDINPUTC );

  /* Find open quote. */
  while ( ( c = *literal ) != '"' && c != '\n' && c != '\0' )
    literal++;
  if ( *literal != '"' ) {
    LALWarning( stat, "No open quote found" );
    RETURN( stat );
  }
  literal++;

  /* Start parsing. */
  while ( n < length - 1 ) {

    /* End of literal, either implicit or explicit. */
    if ( ( c = *(literal++) ) == '\0' || c == '\n' ) {
      LALWarning( stat, "No close quote found" );
      string[n] = '\0';
      RETURN( stat );
    } else if ( c == '"' ) {
      string[n] = '\0';
      RETURN( stat );
    }

    /* Escape sequence. */
    else if ( c == '\\' ) {

      /* Do not allow actual end-of-line or end-of-string to be
         escaped. */
      if ( ( c = *(literal++) ) == '\0' || c == '\n' ) {
	LALWarning( stat, "No close quote found" );
	string[n] = '\0';
	RETURN( stat );
      }

      /* Other special escape characters. */
      else if ( c == 'a' || c == 'A' )
	string[n++] = '\a';
      else if ( c == 'b' || c == 'B' )
	string[n++] = '\b';
      else if ( c == 'f' || c == 'F' )
	string[n++] = '\f';
      else if ( c == 'n' || c == 'N' )
	string[n++] = '\n';
      else if ( c == 'r' || c == 'R' )
	string[n++] = '\r';
      else if ( c == 't' || c == 'T' )
	string[n++] = '\t';
      else if ( c == 'v' || c == 'V' )
	string[n++] = '\v';

      /* Hexadecimal character code. */
      else if ( c == 'x' || c == 'X' ) {
	c = *(literal++);   /* first digit */
	if ( isxdigit( c ) ) {
	  UINT2 value;
	  if ( isdigit( c ) )
	    value = c - '0';
	  else
	    value = 10 + tolower( c ) - 'a';
	  c = *(literal++); /* second digit */
	  if ( isxdigit( c ) ) {
	    value *= 16;
	    if ( isdigit( c ) )
	      value += c - '0';
	    else
	      value += 10 + tolower( c ) - 'a';
	  } else            /* no second digit */
	    literal--;
	  string[n++] = (CHAR)( value );
	  if ( value == 0 ) {
	    LALWarning( stat, "Found explicit end-of-string \\0" );
	    RETURN( stat );
	  }
	} else {            /* no first digit */
	  LALWarning( stat, "Treating empty hex cde as explicit"
		      " end-of-string \\0" );
	  string[n] = '\0';
	  RETURN( stat );
	}
      }

      /* Octal character code. */
      else if ( c >= '0' && c < '8' ) {
	UINT2 value = c - '0';
	c = *(literal++);   /* second digit */
	if ( c >= '0' && c < '8' ) {
	  value *= 8;
	  value += c - '0';
	  c = *(literal++); /* third digit */
	  if ( c >= '0' && c < '8' ) {
	    value *= 8;
	    value += c - '0';
	  } else            /* no third digit */
	    literal--;
	} else              /* no second digit */
	  literal--;
	if ( value > 255 )
	  LALWarning( stat, "Ignoring octal character code >= '\\400'" );
	else
	  string[n++] = (CHAR)( value );
	if ( value == 0 ) {
	  LALWarning( stat, "Found explicit end-of-string \\0" );
	  RETURN( stat );
	}
      }

      /* Other escaped character. */
      else {
	if ( c != '\\' && c != '?' && c != '\'' && c != '"' )
	  LALWarning( stat, "Dropping \\ from unrecognized escape"
		      " sequence" );
	string[n++] = c;
      }
    }

    /* Other character. */
    else
      string[n++] = c;
  }

  if ( *literal != '"' )
    LALWarning( stat, "Reached maximum length before reading close"
		" quote" );
  string[n] = '\0';
  RETURN( stat );
}

/* tell the GNU compiler to ignore issues with the `ll' length modifier */
#ifdef __GNUC__
#define fscanf __extension__ fscanf
#endif

define(`TYPECODE',`I2')dnl
include(`LALReadGrid.m4')dnl

define(`TYPECODE',`I4')dnl
include(`LALReadGrid.m4')dnl

define(`TYPECODE',`I8')dnl
include(`LALReadGrid.m4')dnl

define(`TYPECODE',`U2')dnl
include(`LALReadGrid.m4')dnl

define(`TYPECODE',`U4')dnl
include(`LALReadGrid.m4')dnl

define(`TYPECODE',`U8')dnl
include(`LALReadGrid.m4')dnl

define(`TYPECODE',`S')dnl
include(`LALReadGrid.m4')dnl

define(`TYPECODE',`D')dnl
include(`LALReadGrid.m4')dnl

define(`TYPECODE',`C')dnl
include(`LALReadGrid.m4')dnl

define(`TYPECODE',`Z')dnl
include(`LALReadGrid.m4')dnl
