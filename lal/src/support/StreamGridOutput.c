/**
 * \defgroup StreamGridOutput_c Module StreamGridOutput.c
 * \ingroup StreamOutput_h
 * \author Creighton, T. D.
 *
 * \brief Writes a LAL grid structure to an output stream.
 *
 * ### Prototypes ###
 *
 * \code
 * void
 * LAL<typecode>WriteGrid( LALStatus *stat, FILE *stream, <datatype>Grid *grid )
 * \endcode
 *
 * ### Description ###
 *
 * These routines write the data and metadata in a grid structure
 * <tt>*grid</tt> to an output stream <tt>*stream</tt> in a standard format,
 * described below.  It returns an error if any attempt to write to the
 * stream failed; <tt>*grid</tt> may then be left in a partially-written
 * state.
 *
 * For each of these prototype templates there are in fact 10 separate
 * routines corresponding to all the numeric atomic datatypes
 * <tt>\<datatype\></tt> referred to by <tt>\<typecode\></tt>:
 *
 * <table><tr><th>\<typecode\></th><th>\<datatype\></th><th>\<typecode\></th><th>\<datatype\></th></tr>
 * <tr><td>I2</td><td>  INT2</td><td> U2</td><td>    UINT2</td></tr>
 * <tr><td>I4</td><td>  INT4</td><td> U4</td><td>    UINT4</td></tr>
 * <tr><td>I8</td><td>  INT8</td><td> U8</td><td>    UINT8</td></tr>
 * <tr><td> S</td><td> REAL4</td><td>  C</td><td> COMPLEX8</td></tr>
 * <tr><td> D</td><td> REAL8</td><td>  Z</td><td> COMPLEX16</td></tr>
 * </table>
 *
 * \par Format for <tt>*stream</tt>:
 * The data written to the
 * output stream will be formatted in a manner consistent with the input
 * routines in \ref StreamGridInput_c.  That is, it will begin with a
 * metadata header, consisting of multiple lines of the form:
 *
 * \code
 * # fieldname = value
 * \endcode
 *
 * where \e fieldname is the name of a field in
 * <tt>*series</tt> and \e value is the value of that metadata field,
 * in some standard format (below).  The following metadata fields will
 * be written, one per line:
 *
 * <dl>
 * <dt> datatype:</dt><dd> \e value is a string (\e not
 * surrounded by quotes) corresponding to the type of <tt>*grid</tt>;
 * e.g.\  COMPLEX8Grid.</dd>
 *
 * <dt> name:</dt><dd> \e value is a string surrounded by double-quotes
 * representing <tt>grid-\>name</tt>.  Standard C-language string
 * literal notation is used: printable characters are written directly
 * except for double-quotes and <tt>\\</tt>,
 * characters with special C escape sequences are written
 * as those sequences (e.g.\ <tt>\\t</tt> for tab and <tt>\\n</tt> for
 * newline), and all other character bytes are written as three-digit
 * octal codes <tt>\\</tt>ooo.  Writing stops at the first null byte
 * <tt>\\0</tt>.</dd>
 *
 * <dt> sampleUnits:</dt><dd> \e value is string surrounded by
 * double-quotes; inside the quotes is a unit string corresponding to
 * <tt>grid-\>sampleUnits</tt> as converted by the routine
 * <tt>XLALUnitAsString()</tt>.</dd>
 *
 * <dt> dimUnits:</dt><dd> \e value is a sequence of \f$m\f$ strings,
 * surrounded by double-quotes and separated by a space, where \f$m\f$ is
 * the grid dimension (number of grid axes); inside the quotes is a unit
 * string corresponding to the elements of the <tt>grid-\>dimUnits</tt>
 * array as converted by the routine <tt>XLALUnitAsString()</tt>.</dd>
 *
 * <dt> offset:</dt><dd> \e value is a sequence of \f$m\f$
 * REAL8 numbers separated by single spaces, representing the
 * elements of the <tt>grid-\>offset-\>data</tt>; the number of data \f$m\f$ is
 * the grid dimension and corresponds to the value of
 * <tt>grid-\>offset-\>length</tt>.</dd>
 *
 * <dt> interval:</dt><dd> \e value is a sequence of \f$m\f$
 * REAL8 numbers separated by single spaces, representing the
 * elements of the <tt>grid-\>interval-\>data</tt>; the number of data \f$m\f$ is
 * the grid dimension and corresponds to the value of
 * <tt>grid-\>interval-\>length</tt>.</dd>
 *
 * <dt> dimLength:</dt><dd> \e value is a sequence of \f$M\f$
 * REAL8 numbers separated by single spaces, representing the
 * elements of the <tt>grid-\>data-\>dimLength-\>data</tt>; the number of data
 * \f$M\f$ is the data dimension and corresponds to the value of
 * <tt>grid-\>data-\>dimLength-\>length</tt>, which must be greater than or
 * equal to the grid dimension \f$m\f$, above.</dd>
 * </dl>
 *
 * After all metadata have been written, the contents of
 * <tt>grid-\>data-\>data</tt> will be written in standard integer or
 * floating-point notation, according to <tt>\<datatype\></tt>: integers will
 * be written to full precision, while floating-point numbers will be
 * written in exponential notation with sufficient digits to ensure that
 * they represent a unique binary floating-point number under the IEEE
 * Standard 754 (this means 9 digits for  REAL4s and 17 digits for
 * REAL8s).
 *
 * The input format in \ref StreamGridInput_c does not specify how the
 * numerical data is to be arranged, other than that the numbers be
 * separated by whitespace, and that complex datatypes be represented by
 * alternating real and imaginary components.  These routines adopt the
 * following conventions to improve human-readability: If the data
 * dimension is equal to the grid dimension, then each line consists of a
 * single datum (either a single number, or, for complex datatypes, a
 * pair of numbers separated by whitespace), followed by a newline
 * <tt>'\\n'</tt>.  If the data dimension is greater than the grid
 * dimension, then each line will consist of a number of data equal to
 * the length of the last dimension in <tt>grid-\>data-\>dimLength</tt>.  If
 * the data dimension is at least two greater than the grid dimension,
 * and the dimension lengths are such that a single grid point comprises
 * multiple lines of data, then an additional blank line <tt>'\\n'</tt> is
 * inserted to separate subsequent grid points.
 *
 */

#include <complex.h>
#include <stdio.h>
#include <ctype.h>
#include <lal/LALStdlib.h>
#include <lal/Grid.h>
#include <lal/Units.h>
#include <lal/AVFactories.h>
#include <lal/StringInput.h>
#include <lal/StreamOutput.h>

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
#include "StreamGridOutput_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef COMPLEX

#define TYPECODE C
#define TYPE COMPLEX8
#define FMT "%.8e"
#define COMPLEX 1
#include "StreamGridOutput_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef COMPLEX

#define TYPECODE D
#define TYPE REAL8
#define FMT "%.16e"
#define COMPLEX 0
#include "StreamGridOutput_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef COMPLEX

#define TYPECODE S
#define TYPE REAL4
#define FMT "%.8e"
#define COMPLEX 0
#include "StreamGridOutput_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef COMPLEX

#define TYPECODE I2
#define TYPE INT2
#define FMT "%"LAL_INT2_FORMAT
#define COMPLEX 0
#include "StreamGridOutput_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef COMPLEX

#define TYPECODE I4
#define TYPE INT4
#define FMT "%"LAL_INT4_FORMAT
#define COMPLEX 0
#include "StreamGridOutput_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef COMPLEX

#define TYPECODE I8
#define TYPE INT8
#define FMT "%"LAL_INT8_FORMAT
#define COMPLEX 0
#include "StreamGridOutput_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef COMPLEX

#define TYPECODE U2
#define TYPE UINT2
#define FMT "%"LAL_UINT2_FORMAT
#define COMPLEX 0
#include "StreamGridOutput_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef COMPLEX

#define TYPECODE U4
#define TYPE UINT4
#define FMT "%"LAL_UINT4_FORMAT
#define COMPLEX 0
#include "StreamGridOutput_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef COMPLEX

#define TYPECODE U8
#define TYPE UINT8
#define FMT "%"LAL_UINT8_FORMAT
#define COMPLEX 0
#include "StreamGridOutput_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef COMPLEX
