/**
   \defgroup StreamGridInput_c Module StreamGridInput.c
   \ingroup StreamInput_h
\author Creighton, T. D.

\brief Converts an input stream into a LAL grid structure.

\heading{Prototypes}

\code
void
LAL<typecode>ReadGrid( LALStatus *stat, <datatype>Grid **grid, FILE *stream )
\endcode

\heading{Description}

These routines parse an input stream <tt>*stream</tt> to create a grid
structure <tt>**grid</tt>, and fill in its data and metadata fields.
The output parameter \c grid must be a non-\c NULL handle to a
\c NULL-valued pointer <tt>*grid</tt>.

This prototype template in fact refers to 10 separate routines
corresponding to all the numerical atomic datatypes <tt>\<datatype\></tt>
referred to by <tt>\<typecode\></tt>:

<table>
<tr><th>\<typecode\></th><th>\<datatype\></th><th>\<typecode\></th><th>\<datatype\></th></tr>
<tr><td>I2</td><td> INT2</td><td>U2</td><td>   UINT2</td></tr>
<tr><td>I4</td><td> INT4</td><td>U4</td><td>   UINT4</td></tr>
<tr><td>I8</td><td> INT8</td><td>U8</td><td>   UINT8</td></tr>
<tr><td> S</td><td>REAL4</td><td> C</td><td>COMPLEX8</td></tr>
<tr><td> D</td><td>REAL8</td><td> Z</td><td>COMPLEX16</td></tr>
</table>

\heading{Format for *stream:} The input stream is assumed
to be a text stream (ASCII) consisting of a header containing metadata
followed by numerical data in standard integer or floating-point
format, as recognized by the routines in \ref StringConvert.c.  The
header consists of zero or more lines beginning with a \c \#
character, followed by a metadata field name and value in the format:

\code
# fieldname=value
\endcode

The = sign in this format is standard but optional;
it may be replaced or surrounded with any amount of any whitespace
except a newline \\n.  If \e fieldname is unrecognized,
it is ignored; if it is recognized, then \e value must be in a
suitable format for the field type, as described below.  Blank lines,
or lines containing just a \# character, are skipped.  Once a
line is encountered that contains non-whitespace characters and does
not start with \#, that line is assumed to be the beginning of
the numerical data.  From that point on, all non-whitespace characters
must be part of parseable numbers; no more comments are permitted
(although blank lines will still be skipped).

If a metadata field appears twice in the header, the later one takes
precedence.  At present these routines do not track which fields have
been previously assigned, so no warnings or errors are generated.
Some metadata fields are \e required, and the routine will abort
if it doesn't find them.  Others are \e optional: if they are not
found, the routine will assign some default value.  The various fields
and their required formats are given below:

\heading{Required fields:}
<dl>
<dt>dimLength:</dt><dd> \e value consists of a sequence of
\c UINT4s separated by whitespace (but \e not a newline
<tt>\\n</tt>).  These are used to create <tt>(*grid)-\>data</tt> with the
appropriate dimensions, and are stored in
<tt>(*grid)-\>data-\>dimLength-\>data</tt>: the number of integers \f$M\f$
gives the data dimension number (the number of array indecies), while
the value of each integer gives the length of each dimension (the
range of values of the corresponding array index).</dd>

<dt>offset:</dt><dd> \e value consists of a sequence of
\c REAL8s separated by whitespace (but \e not a newline
<tt>\\n</tt>).  These values are stored in
<tt>(*grid)-\>offset-\>data</tt>.  The number of data \f$m\f$ gives the grid
dimension, which must be less than or equal to the data dimension \f$M\f$
of the array defined above, and must be consistent among the
\c offset, \c interval, and \c dimUnits fields.</dd>

<dt>interval:</dt><dd> \e value consists of a sequence of
\c REAL8s separated by whitespace (but \e not a newline
<tt>\\n</tt>).  These values are stored in
<tt>(*grid)-\>interval-\>data</tt>.  The number of data \f$m\f$ gives the grid
dimension, which must be less than or equal to the data dimension \f$M\f$
of the array defined above, and must be consistent among the
\c offset, \c interval, and \c dimUnits fields.</dd>
</dl>

\heading{Optional fields:}
<dl>

<dt>name:</dt><dd> \e value is a string surrounded by double-quotes,
which is parsed in the manner of a string literal in C: it
may contain ordinary printable characters (except double-quotes and
<tt>\\</tt>), escape sequences (such as <tt>\\t</tt> for tab, <tt>\\n</tt> for
newline, or <tt>\\</tt> and double-quotes for literal backslash and quote
characters), and octal or hexadecimal codes (<tt>\\</tt>ooo or
<tt>\\x</tt>hh respectively) for arbitrary bytes.  Unlike in C,
literals cannot be split between lines, adjacent literals are not
concatenated, and converted strings longer than
\c LALNameLength-1 will be truncated.  The resulting string is
stored in <tt>(*grid)-\>name</tt>, and will always contain a <tt>\\0</tt>
terminator, beyond which the contents are unspecified.  If this field
is not given in \c stream, then the routine will simply assign
<tt>(*grid)-\>name[0]</tt>=<tt>\\0</tt>.
</dd>

<dt>sampleUnits:</dt><dd> \e value is string surrounded by
double-quotes; the quotes are stripped and the string passed to
<tt>LALParseUnitString()</tt> to determine <tt>(*grid)-\>sampleUnits</tt>.
Since <tt>LALParseUnitString()</tt> is not very robust, it is
recommended to use only unit strings that have been generated by
<tt>LALUnitAsString()</tt>, or to remove this metadata field and set
<tt>(*grid)-\>sampleUnits</tt> within the code.  If this field is not
given in \c stream, then \c lalDimensionlessUnit is assumed.
</dd>

<dt>dimUnits:</dt><dd> \e value is a sequence of strings,
each surrounded by double-quotes, and optionally separated by
whitespace (but \e not a newline <tt>\\n</tt>); the quotes are
stripped and the string passed to LALParseUnitString() to
determine <tt>(*grid)-\>dimUnits</tt>.  Since LALParseUnitString()
is not very robust, it is recommended to use only unit strings that
have been generated by LALUnitAsString(), or to remove this
metadata field and reset <tt>(*grid)-\>dimUnits</tt> within the code.  If
this field is not given in \c stream, then
<tt>(*grid)-\>dimUnits</tt> will be allocated as an array containing \f$m\f$
units assigned the value \c lalDimensionlessUnit, where \f$m\f$ is the
length of the \c interval and \c offset vectors.  If this
field \e is given, then the number of unit strings must be
consistent with the lengths of the \c interval and \c offset
vectors, or the routine will abort.
</dd>

<dt>datatype:</dt><dd> \e value is string identifying the
grid type; e.g. REAL4Grid (\e not surrounded by quotes).
This should correspond to the type of <tt>**grid</tt>, not to any field
in <tt>**grid</tt>.  If there is a type mismatch, a warning is generated
(and errors may occur later while parsing the data).
</dd>
</dl>

\heading{Data format:} The first line that is neither blank nor
beginning with a <tt>\#</tt> character is assumed to be the start of
the grid data, and is parsed as a sequence of whitespace-separated
integers or real numbers.  For complex datatypes, the numbers read are
interpreted as alternately the real and imaginary parts of the data.
Since the number of data required is known from the grid metadata, the
routine will read from \c stream only until the structure is
filled (except that the first line of data will be read in its
entirety, even if it contains more numbers than required).  If the
routine encounters an unparseable number or the end-of-input before
filling the structure, an error is returned.


\heading{Algorithm}

These routines use LALCHARReadVector() to read the header lines
and the first line of data.  The metadata are stored in temporary
variables or vectors.  After the first data line has been read, the
size of the grid should be known from the metadata: the structure is
created, and the metadata and first line of data are copied in.  After
this, any additional data are parsed directly from \c stream using
<tt>fscanf()</tt>.  This is done for efficiency: repeated calling of the
LAL string parsing routines in \ref StringConvert.c involves far too
much computational overhead.

*/

#include <complex.h>
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
#include "StreamGridInput_source.c"
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
#include "StreamGridInput_source.c"
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
#include "StreamGridInput_source.c"
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
#include "StreamGridInput_source.c"
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
#include "StreamGridInput_source.c"
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
#include "StreamGridInput_source.c"
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
#include "StreamGridInput_source.c"
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
#include "StreamGridInput_source.c"
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
#include "StreamGridInput_source.c"
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
#include "StreamGridInput_source.c"
#undef TYPECODE
#undef DATACODE
#undef TYPE
#undef DATA
#undef SIZE
#undef COMPLEX
