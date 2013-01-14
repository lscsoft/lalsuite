/**
   \defgroup StreamSeriesInput_c Module StreamSeriesInput.c
   \ingroup StreamInput_h
\author Creighton, T. D.

   \brief Converts an input stream into a time or frequency series.

\heading{Prototypes}

\code
void
LAL<typecode>ReadTSeries( LALStatus            *stat,
                          <datatype>TimeSeries *series,
                          FILE                 *stream )

void
LAL<typecode>ReadTVectorSeries( LALStatus                  *stat,
                                <datatype>TimeVectorSeries *series,
                                FILE                       *stream )

void
LAL<typecode>ReadTArraySeries( LALStatus                 *stat,
                               <datatype>TimeArraySeries *series,
                               FILE                      *stream )

void
LAL<typecode>ReadFSeries( LALStatus                 *stat,
                          <datatype>FrequencySeries *series,
                          FILE                      *stream )
\endcode


\heading{Description}

These routines parse an input stream <tt>*stream</tt> to fill in the
data and metadata fields of a time or frequency series <tt>*series</tt>.
The field <tt>series-\>data</tt> must be \c NULL, so that it can be
created and filled by the routine.  The other fields may be
initialized or not; they will be overwritten by metadata read from
<tt>*stream</tt>.  If an error occurs, <tt>*series</tt> will be left
unchanged, but <tt>*stream</tt> will have been read up to the point
where the error occured.

For each of these prototype templates there are in fact 10 separate
routines corresponding to all the atomic datatypes <tt>\<datatype\></tt>
(except \c CHAR) referred to by <tt>\<typecode\></tt>:

<table>
<tr><th>\<typecode\></th><th>\<datatype\></th><th>\<typecode\></th><th>\<datatype\></th></tr>
<tr><td>I2</td><td> INT2</td><td> U2</td><td>   UINT2</td></tr>
<tr><td>I4</td><td> INT4</td><td> U4</td><td>   UINT4</td></tr>
<tr><td>I8</td><td> INT8</td><td> U8</td><td>   UINT8</td></tr>
<tr><td> S</td><td>REAL4</td><td>  C</td><td>COMPLEX8</td></tr>
<tr><td> D</td><td>REAL8</td><td>  Z</td><td>COMPLEX16</td></tr>
</table>

\heading{Format for <tt>*stream</tt>:} The input stream is assumed
to be a text stream (ASCII) consisting of a header containing metadata
followed by numerical data in standard integer or floating-point
format, as recognized by the routines in \ref StringConvert.c.  The
header consists of zero or more lines beginning with a \c \#
character, followed by a metadata field name and value in the format:

\code
# fieldname=value
\endcode

The <tt>=</tt> sign in this format is standard but optional;
it may be replaced or surrounded with any amount of any whitespace
except a newline <tt>\\n</tt>.  If \e fieldname is unrecognized,
it is ignored; if it is recognized, then \e value must be in a
suitable format for the field type, as described below.  Blank lines,
or lines containing just a \c \# character, are skipped.  Once a
line is encountered that contains non-whitespace characters and does
not start with \c \#, that line is assumed to be the beginning of
the numerical data.  From that point on, all non-whitespace characters
must be part of parseable numbers; no more comments are permitted
(although blank lines will still be skipped).

If a metadata field appears twice in the header, the later one takes
precedence.  At present these routines do not track which fields have
been previously assigned, so no warnings or errors are generated.

How the data is packed into the <tt>series-\>data</tt> structure depends
on what metadata has been provided, as described below.

<b>Required, conditional, and optional metadata:</b>

The input stream need not contain a complete set of metadata, allowing some
metadata to be read from <tt>*stream</tt> and others to be set
elsewhere.  For each type of series, some metadata will be
\e required, and the routine will abort if the metadata is not
found.  Other metadata are \e conditional, meaning that the
routine will operate differently depending on whether or not these
metadata were found.  The remaining metadata are \e optional; if
they are not found in <tt>*stream</tt>, they will be left unchanged.
The recognized metadata fields are listed below.

<tt>\<datatype\>TimeSeries</tt>:
<dl>
<dt>Required fields:</dt><dd> none</dd>
<dt>Conditional fields:</dt><dd> \c length</dd>
<dt>Optional fields:</dt><dd> \c name, \c epoch, \c deltaT, \c f0, \c sampleUnits, \c datatype</dd>
</dl>

<tt>\<datatype\>TimeVectorSeries</tt>:
<dl>
<dt>Required fields:</dt><dd> none</dd>
<dt>Conditional fields:</dt><dd> \c length, \c vectorLength</dd>
<dt>Optional fields:</dt><dd> \c name, \c epoch, \c deltaT, \c f0, \c sampleUnits, \c datatype</dd>
</dl>

<tt>\<datatype\>TimeArraySeries</tt>:
<dl>
<dt>Required fields:</dt><dd> \c dimLength</dd>
<dt>Conditional fields:</dt><dd> \c length, \c arrayDim</dd>
<dt>Optional fields:</dt><dd> \c name, \c epoch, \c deltaT, \c f0, \c sampleUnits, \c datatype</dd>
</dl>

<tt>\<datatype\>FrequencySeries</tt>:
<dl>
<dt>Required fields:</dt><dd> none</dd>
<dt>Conditional fields:</dt><dd> \c length</dd>
<dt>Optional fields:</dt><dd> \c name, \c epoch, \c deltaT, \c f0, \c deltaF, \c sampleUnits, \c datatype</dd>
</dl>

Below we describe the required format for the field values, as well as
what occurs if a conditional field is or isn't present.

\heading{Required fields:}
<dl>
<dt>dimLength</dt><dd> (\c TimeArraySeries only):
\e value consists of a sequence of \c UINT4s separated by
whitespace (but \e not a newline <tt>'\\n'</tt>).  These data are
stored in <tt>series-\>data-\>dimLength</tt>: the number of integers gives
the number of array indecies, while the value of each integer gives
the dimension of the corresponding array index.</dd>
</dl>

\heading{Conditional fields:}
<dl>
<dt>arrayDim</dt><dd> (\c TimeArraySeries only): \e value
is a single \c UINT4, to be stored in
<tt>series-\>data-\>arrayDim</tt>.  This must equal the product of the
index ranges in \c dimLength, above, or an error is returned.  If
not given, the \c arrayDim field will be set equal to the product
of the index ranges in \c dimLength.  (The \c arrayDim and
\c dimLength fields can appear in any order in <tt>*stream</tt>;
checking is done only after all header lines have been read.)</dd>

<dt>vectorLength</dt><dd> (\c TimeVectorSeries only):
\e value is a single \c UINT4, to be stored in
<tt>series-\>data-\>vectorLength</tt>.  If not specified in the header
portion of <tt>*stream</tt>, it will be taken to be the number of data
on the \e first line of the data portion of <tt>*stream</tt>, or
half the number of real data for a complex-valued
\c TimeVectorSeries; if an odd number of real data are found on
the first line of a complex \c TimeVectorSeries, then an error is
returned.</dd>

<dt>length:</dt><dd> \e value is a single \c UINT4, to be
stored in <tt>series-\>data-\>length</tt>.  If it is specified in the
header portion of <tt>*stream</tt>, data will be read until
\c length is reached.  Otherwise, <tt>*stream</tt> will be read to
its end or until an unparseable character is read, and \c length
will then be set accordingly.  (If parsing stops in the middle of
filling a complex, vector, or array valued element, the partly-read
element is discarded.)</dd>
</dl>

\heading{Optional fields:}

<dl>
<dt>name:</dt><dd>\c value is a string surrounded by double-quotes,
which is parsed in the manner of a string literal in C: it
may contain ordinary printable characters (except double-quote and \\),
escape sequences (such as \\t for tab, \\n for
newline, or \\ and double-quote literal backslash and quote
characters), and octal or hexadecimal codes (\\\c ooo or
\\x\c hh, respectively) for arbitrary bytes.  Unlike in C,
literals cannot be split between lines, adjacent literals are not
concatenated, and converted strings longer than
\c LALNameLength-1 will be truncated.  The resulting string is
stored in <tt>series-\>name</tt>, and will always contain a \c \\0
terminator, beyond which the contents are unspecified.</dd>

<dt>epoch:</dt><dd> \e value is a single \c INT8 number
of GPS nanoseconds, or a pair of \c INT4s representing GPS seconds
and nanoseconds separately, separated by non-newline whitespace.</dd>

<dt>deltaT</dt><dd> (any time series): \e value is a single
\c REAL8 number.</dd>

<dt>f0:</dt><dd> \e value is a single \c REAL8 number.</dd>

<dt>deltaF</dt><dd> (\c FrequencySeries only): \e value
is a single \c REAL8 number.</dd>

<dt>sampleUnits:</dt><dd> \e value is string surrounded by
double-quotes; the quotes are stripped and the string passed to
<tt>LALParseUnitString()</tt> to determine <tt>series-\>sampleUnits</tt>.
Since <tt>LALParseUnitString()</tt> is not very robust, it is
recommended to use only unit strings that have been generated by
<tt>LALUnitAsString()</tt>, or to remove this metadata field and set
<tt>series-\>sampleUnits</tt> within the code.</dd>

<dt>datatype:</dt><dd> \e value is string identifying the
series type; e.g. \c REAL4TimeSeries (\e not surrounded by
quotes).  This should correspond to the type of <tt>*series</tt>, not to
any field in <tt>*series</tt>.  If there is a type mismatch, a warning
is generated (and errors may occur later while parsing the data).</dd>

</dl>

\heading{Data format:} The data portion of <tt>*stream</tt> consists
of whitespace-separated integer or real numbers.  For complex input
routines, the real data are parsed as alternately the real and
imaginary parts of successive complex numbers.  By convention, each
line should correspond to a single base, complex, vector, or array
valued element of the <tt>series-\>data</tt> sequence.  However, this is
\e required only in the case of a \c TimeVectorSeries where
the \c vectorLength metadata was not set in the header, since in
this case the value of \c vectorLength will be taken from the
number of elements read on the first data line.  After this, and in
all other cases, newlines are treated as any other whitespace.

If a \c length value is specified in the header, then data are
read until the required length is acheived; if <tt>fscanf()</tt> returns
zero or negative before this (representing either the end-of-input or
a character that cannot be interpreted as part of the numerical data),
an error is returned.  If a \c length value was not specified,
data are read until <tt>fscanf()</tt> returns zero or negative: at this
point any partially-completed complex, vector, or array valued element
is discarded, and <tt>series-\>data-\>length</tt> set to the number of
elements read.

\heading{Algorithm}

These routines use <tt>LALCHARReadVector()</tt> to read the header lines
and the first line of data.  After this, data are parsed directly from
<tt>*stream</tt> using <tt>fscanf()</tt>.  This is done for efficiency:
repeated calling of the LAL string parsing routines in
\ref StringConvert.c involves far too much computational overhead.

After the first data line has been read, the length of each sequence
element will be known from the atomic type, as well as the specified
\c dimLength (for arrays), \c vectorLength (for vectors), or
number of elements on the first data line (for vectors without an
explicitly specified \c vectorLength).  If \c length is also
specified, a sequence of the appropriate size is allocated, and all
the data is copied or read directly into it.  If \c length was not
specified, the data read with <tt>fscanf()</tt> are stored in a linked
list of buffers of size \c BUFFSIZE (a local <tt>\# define</tt>d
constant) until parsing stops.  Then a sequence of the appropriate
size is allocated and the data copied into it.

*/

#include <complex.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/Units.h>
#include <lal/StringInput.h>
#include <lal/StreamInput.h>

/* Define a message string for header parsing errors. */
#define LALREADSERIESC_HEADER "Skipping badly-formatted line for metadata field "

/* Define linked-list of buffers for storing an arbitrary number of
   arbitrary datatypes.  BUFFSIZE should be a multiple of 16. */
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

/* Define a function for parsing a string literal. */
static void
LALLiteralToString( LALStatus  *stat,
		    CHAR       *string,
		    const CHAR *literal,
		    UINT4      length )
{
  CHAR c;       /* Current character being considered. */
  UINT4 n = 0;  /* Counter of number of characters written. */

  INITSTATUS(stat);

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

#define TYPECODE I2
#define DATACODE TYPECODE
#define TYPE INT2
#define DATA TYPE
#define SIZE 2
#define COMPLEX 0
#include "StreamSeriesInput_source.c"
#undef TYPECODE
#undef DATACODE
#undef TYPE
#undef DATA
#undef SIZE
#undef COMPLEX

#define TYPECODE I4
#define DATACODE TYPECODE
#define TYPE INT4
#define DATA TYPE
#define SIZE 4
#define COMPLEX 0
#include "StreamSeriesInput_source.c"
#undef TYPECODE
#undef DATACODE
#undef TYPE
#undef DATA
#undef SIZE
#undef COMPLEX

#define TYPECODE I8
#define DATACODE TYPECODE
#define TYPE INT8
#define DATA TYPE
#define SIZE 8
#define COMPLEX 0
#include "StreamSeriesInput_source.c"
#undef TYPECODE
#undef DATACODE
#undef TYPE
#undef DATA
#undef SIZE
#undef COMPLEX

#define TYPECODE U2
#define DATACODE TYPECODE
#define TYPE UINT2
#define DATA TYPE
#define SIZE 2
#define COMPLEX 0
#include "StreamSeriesInput_source.c"
#undef TYPECODE
#undef DATACODE
#undef TYPE
#undef DATA
#undef SIZE
#undef COMPLEX

#define TYPECODE U4
#define DATACODE TYPECODE
#define TYPE UINT4
#define DATA TYPE
#define SIZE 4
#define COMPLEX 0
#include "StreamSeriesInput_source.c"
#undef TYPECODE
#undef DATACODE
#undef TYPE
#undef DATA
#undef SIZE
#undef COMPLEX

#define TYPECODE U8
#define DATACODE TYPECODE
#define TYPE UINT8
#define DATA TYPE
#define SIZE 8
#define COMPLEX 0
#include "StreamSeriesInput_source.c"
#undef TYPECODE
#undef DATACODE
#undef TYPE
#undef DATA
#undef SIZE
#undef COMPLEX

#define TYPECODE S
#define DATACODE TYPECODE
#define TYPE REAL4
#define DATA TYPE
#define SIZE 4
#define COMPLEX 0
#include "StreamSeriesInput_source.c"
#undef TYPECODE
#undef DATACODE
#undef TYPE
#undef DATA
#undef SIZE
#undef COMPLEX

#define TYPECODE D
#define DATACODE TYPECODE
#define TYPE REAL8
#define DATA TYPE
#define SIZE 8
#define COMPLEX 0
#include "StreamSeriesInput_source.c"
#undef TYPECODE
#undef DATACODE
#undef TYPE
#undef DATA
#undef SIZE
#undef COMPLEX

#define TYPECODE Z
#define DATACODE D
#define TYPE COMPLEX16
#define DATA REAL8
#define SIZE 8
#define COMPLEX 1
#include "StreamSeriesInput_source.c"
#undef TYPECODE
#undef DATACODE
#undef TYPE
#undef DATA
#undef SIZE
#undef COMPLEX

#define TYPECODE C
#define DATACODE S
#define TYPE COMPLEX8
#define DATA REAL4
#define SIZE 4
#define COMPLEX 1
#include "StreamSeriesInput_source.c"
#undef TYPECODE
#undef DATACODE
#undef TYPE
#undef DATA
#undef SIZE
#undef COMPLEX
