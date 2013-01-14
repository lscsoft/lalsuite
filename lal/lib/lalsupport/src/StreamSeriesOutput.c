/**
   \defgroup StreamSeriesOutput_c Module StreamSeriesOutput.c
   \ingroup StreamOutput_h
\author Creighton, T. D.

\brief Writes a time or frequency series to an output stream.

\heading{Prototypes}

\code
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
\endcode


\heading{Description}

These routines write the data and metadata in a time or frequency
series <tt>*series</tt> to an output stream <tt>*stream</tt> in a standard
format, described below.  It returns an error if any attempt to write
to the stream failed; <tt>*stream</tt> may then be left in a
partially-written state.

For each of these prototype templates there are in fact 10 separate
routines corresponding to all the atomic datatypes <tt>\<datatype\></tt>
(except \c CHAR) referred to by <tt>\<typecode\></tt>:

<table><tr><th>\<typecode\></th><th>\<datatype\></th><th>\<typecode\></th><th>\<datatype\></th></tr>
<tr><td> I2</td><td>  INT2</td><td> U2</td><td>    UINT2</td></tr>
<tr><td> I4</td><td>  INT4</td><td> U4</td><td>    UINT4</td></tr>
<tr><td> I8</td><td>  INT8</td><td> U8</td><td>    UINT8</td></tr>
<tr><td>  S</td><td> REAL4</td><td>  C</td><td> COMPLEX8</td></tr>
<tr><td>  D</td><td> REAL8</td><td>  Z</td><td> COMPLEX16</td></tr>
</table>

\heading{Format for <tt>*stream</tt>:} The data written to the
output stream will be formatted in a manner consistent with the input
routines in \ref StreamSeriesInput.c.  That is, it will begin with a
metadata header, consisting of multiple lines of the form:

\code
# fieldname=value
\endcode

where \e fieldname is the name of a field in
<tt>*series</tt> and \e value is the value of that metadata field,
in some standard format (below).  The following metadata fields will
be written, one per line, based on the type of <tt>*series</tt>:

<dl>
<dt><tt>\<datatype\>TimeSeries</tt>:</dt><dd> \c datatype, \c name,\c epoch, \c deltaT, \c f0, \c sampleUnits,\c length</dd>
<dt><tt>\<datatype\>TimeVectorSeries</tt>:</dt><dd> \c datatype,\c name, \c epoch, \c deltaT, \c f0,\c sampleUnits, \c length, \c vectorLength</dd>
<dt><tt>\<datatype\>TimeArraySeries</tt>:</dt><dd> \c datatype,\c name, \c epoch, \c deltaT, \c f0,\c sampleUnits, \c length, \c dimLength, \c arrayDim</dd>
<dt><tt>\<datatype\>FrequencySeries</tt>:</dt><dd> \c datatype,\c name, \c epoch, \c deltaT, \c f0, \c deltaF,\c sampleUnits, \c length</dd>
</dl>

After all metadata have been written, the contents of
<tt>series-\>data-\>data</tt> will be written in standard integer or
floating-point notation, according to <tt>\<datatype\></tt>: integers will
be written to full precision, while floating-point numbers will be
written in exponential notation with sufficient digits to ensure that
they represent a unique binary floating-point number under the IEEE
Standard 754 (this means 9 digits for \c REAL4s and 17 digits for
\c REAL8s).  Complex datatypes are represented by pairs of
floating-point numbers representing alternately the real and imaginary
parts.

The body of the file will be formatted with newlines <tt>'\\n'</tt>
separating individual base, complex, vector, or array valued elements
of the sequence <tt>series-\>data</tt>.  Within each element, integer or
floating-point components will be separated by single <tt>' '</tt>
characters.  Thus the value of <tt>series-\>data-\>length</tt> will always
equal the number of lines following the metadata header.

\heading{Format for metadata fields:} Here we summarize briefly the
format for the individual field values in the metadata header.

<dl>
<dt>datatype:</dt><dd> \e value is a string (\e not
surrounded by quotes) corresponding to the type of <tt>*series</tt>;
e.g.\ COMPLEX8FrequencySeries.</dd>

<dt>name:</dt><dd> \e value is a string surrounded by double-quotes
representing <tt>series-\>name</tt>.  Standard C-language string
literal notation is used: printable characters are written directly
except for double-quotes and <tt>\\</tt> (rendered as <tt>\\</tt>),
characters with special C escape sequences are written
as those sequences (e.g.\ <tt>\\t</tt> for tab and <tt>\\n</tt> for
newline), and all other character bytes are written as three-digit
octal codes <tt>\\</tt>ooo.  Writing stops at the first null byte
<tt>\\0</tt>.</dd>

<dt>epoch:</dt><dd> \e value is a single \c INT8 number
representing <tt>series-\>epoch</tt> in GPS nanoseconds.</dd>

<dt>deltaT</dt><dd> (any time series): \e value is a single
\c REAL8 number representing <tt>series-\>deltaT</tt>.</dd>

<dt>f0:</dt><dd> \e value is a single \c REAL8 number
representing <tt>series-\>f0</tt>.</dd>

<dt>deltaF</dt><dd> (\c FrequencySeries only): \e value
is a single \c REAL8 number representing <tt>series-\>deltaF</tt>.</dd>

<dt>sampleUnits:</dt><dd> \e value is string surrounded by
double-quotes; inside the quotes is a unit string corresponding to
<tt>series-\>sampleUnits</tt> as converted by the routine
LALUnitAsString().</dd>

<dt>length:</dt><dd> \e value is a single \c UINT4
representing <tt>series-\>data-\>length</tt>.</dd>

<dt>vectorLength</dt><dd> (\c TimeVectorSeries only):
\e value is a single \c UINT4 representing
<tt>series-\>data-\>vectorLength</tt>.</dd>

<dt>dimLength</dt><dd> (\c TimeArraySeries only):
\e value consists of a sequence of \c UINT4s separated by
single <tt>' '</tt> characters, representing the components of
<tt>series-\>data-\>dimLength-\>data</tt>.  The value of
<tt>series-\>data-\>dimLength-\>length</tt> must be inferred from the
number of components; it is not given as separate metadata.</dd>

<dt>arrayDim</dt><dd> (\c TimeArraySeries only): \e value
is a single \c UINT4 representing <tt>series-\>data-\>arrayDim</tt>.
If the array sequence was properly constructed, this will equal the
product of the components of \c dimLength, above.</dd>
</dl>

*/

#include <complex.h>
#include <stdio.h>
#include <ctype.h>
#include <lal/LALStdlib.h>
#include <lal/Units.h>
#include <lal/AVFactories.h>
#include <lal/StringInput.h>
#include <lal/StreamOutput.h>
#include <lal/Date.h>

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

#define TYPECODE Z
#define TYPE COMPLEX16
#define FMT "%.16e"
#define COMPLEX 1
#include "StreamSeriesOutput_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef COMPLEX

#define TYPECODE C
#define TYPE COMPLEX8
#define FMT "%.8e"
#define COMPLEX 1
#include "StreamSeriesOutput_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef COMPLEX

#define TYPECODE D
#define TYPE REAL8
#define FMT "%.16e"
#define COMPLEX 0
#include "StreamSeriesOutput_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef COMPLEX

#define TYPECODE S
#define TYPE REAL4
#define FMT "%.8e"
#define COMPLEX 0
#include "StreamSeriesOutput_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef COMPLEX

#define TYPECODE I2
#define TYPE INT2
#define FMT "%"LAL_INT2_FORMAT
#define COMPLEX 0
#include "StreamSeriesOutput_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef COMPLEX

#define TYPECODE I4
#define TYPE INT4
#define FMT "%"LAL_INT4_FORMAT
#define COMPLEX 0
#include "StreamSeriesOutput_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef COMPLEX

#define TYPECODE I8
#define TYPE INT8
#define FMT "%"LAL_INT8_FORMAT
#define COMPLEX 0
#include "StreamSeriesOutput_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef COMPLEX

#define TYPECODE U2
#define TYPE UINT2
#define FMT "%"LAL_UINT2_FORMAT
#define COMPLEX 0
#include "StreamSeriesOutput_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef COMPLEX

#define TYPECODE U4
#define TYPE UINT4
#define FMT "%"LAL_UINT4_FORMAT
#define COMPLEX 0
#include "StreamSeriesOutput_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef COMPLEX

#define TYPECODE U8
#define TYPE UINT8
#define FMT "%"LAL_UINT8_FORMAT
#define COMPLEX 0
#include "StreamSeriesOutput_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef COMPLEX
