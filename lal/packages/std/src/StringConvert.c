/******************************** <lalVerbatim file="StringConvertCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{StringConvert.c}}
\label{ss:StringConvert.c}

Converts a string into a numerical value.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{StringConvertCP}
\idx{LALStringToI2()}
\idx{LALStringToI4()}
\idx{LALStringToI8()}
\idx{LALStringToU2()}
\idx{LALStringToU4()}
\idx{LALStringToU8()}
\idx{LALStringToS()}
\idx{LALStringToD()}
\idx{LALStringToC()}
\idx{LALStringToZ()}

\subsubsection*{Description}

These routines parse the string \verb@*string@ and compute a numerical
value \verb@*value@ of the appropriate datatype.  If
\verb@endptr@$\neq$\verb@NULL@, then after conversion \verb@*endptr@
will point to the character after the last character used in the
conversion.  The routine will always return without error if the
arguments are valid, regardless of the contents of \verb@string@;
failure to parse a number is indicated by \verb@*endptr@ being set
equal to \verb@string@ (provided \verb@endptr@$\neq$\verb@NULL@).

For integer or floating-point conversion to occur, \verb@string@ must
consist of zero or more whitespace characters followed by a number in
any standard base-ten integer or floating-point representation
(described in detail below); the conversion will stop as soon as the
routine encounters a character that is not part of the number.  For
instance, parsing the string \verb@"   123bad"@ will return
\verb@*value@=123, and \verb@*endptr@ will point to the substring
\verb@"bad"@; it is up to the calling routine to determine whether
this is an acceptable input.  By contrast, parsing the string
\verb@" bad"@ will leave \verb@*value@ unchanged and \verb@*endptr@
pointing to the start of the original string.  In general, if the
routine returns with \verb@*endptr@$\neq$\verb@string@ and
\verb@**endptr@ is a whitespace or \verb@'\0'@ character, then the
format was unambiguously acceptable.

Complex conversion is essentially equivalent to performing
floating-point conversion on \verb@string@ to get the real part, and
again on \verb@*endptr@ to get the imaginary part.  Normally this
means that an acceptable format is two floating-point representations
separated by (and possibly preceded with) whitespace, although the
intervening whitespace can be omitted if the separation between the
two numbers is unambiguous.  Thus the string \verb@"-1.0+3.5"@ will be
unambiguously read as $-1.0+3.5i$, but \verb@"-1.03.5"@ will be read
as $-1.03+0.5i$ (since the conversion of the real part stops at the
second \verb@'.'@ character), which may or may not be the intended
conversion.

\subsubsection*{Algorithm}

These functions emulate the standard C functions \verb@strtol()@,
\verb@strtoul()@, and \verb@strtod()@, except that they follow LAL
calling conventions and return values of the appropriate LAL
datatypes.  For integer conversion, only base-ten (decimal)
representations are supported.  Otherwise, the valid format is as for
the corresponding C funcions, which we summarize below:

A string to be converted to an \verb@INT@$n$ (where $n$=2, 4, or 8)
consists of zero or more whitespace characters as determined by
\verb@isspace()@, followed optionally by a single \verb@'+'@ or
\verb@'-'@ character, followed by one or more decimal digits;
conversion stops at the first non-digit character after this.  If the
result would overflow or underflow the \verb@INT@$n$ representation,
then the value is set to \verb@LAL_INT@$n$\verb@_MAX@$=2^{8n-1}-1$ or
\verb@LAL_INT@$n$\verb@_MIN@$=-2^{8n-1}$, respectively.

A string to be converted to a \verb@UINT@$n$ follows the same format,
except that a leading negative sign character \verb@'-'@ is
\emph{ignored} (the routine will compute the magnitude of the result),
and the return value is capped at
\verb@LAL_UINT@$n$\verb@_MAX@$=2^{8n}-1$.

A string to be converted to a floating-point number (\verb@REAL4@ or
\verb@REAL8@) consists of zero or more whitespace characters as
determined by \verb@isspace()@, followed optionally by a single
\verb@'+'@ or \verb@'-'@ character, followed by a sequence of one or
more decimal digits optionally containing a decimal point \verb@'.'@,
optionally followed by an exponent.  An exponent consists of a single
\verb@'E'@ or \verb@'e'@ character, followed optionally by a single
\verb@'+'@ or \verb@'-'@ character, followed by a sequence of one or
more decimal digits.  If the converted value would overflow,
$\pm$\verb@LAL_REAL@$n$\verb@_MAX@ is returned, as appropriate.  If
the value would underflow, 0 is returned.

A string to be converted to a complex number (\verb@COMPLEX8@ or
\verb@COMPLEX16@) consists of two floating-point format substrings
concatenated together, where the first character of the second
substring cannot be interpreted as a continuation of the first number.
Usually this means that the second substring will contain at lead one
leading whitespace character, though this is not strictly necessary.
Overflow or underflow is dealt with as above.

Internally, the floating-point conversion routines call
\verb@strtod()@, then cap and cast the result as necessary.  The
complex conversion routines simply call their floating-point
counterpart twice.  The integer routines call an internal function
\verb@LALStringToU8AndSign()@, which does what you would expect, then
cap and cast the result as necessary.  (The C routines \verb@strtol()@
and \verb@strtol()@ are not used as they are not guaranteed to have
8-byte precision.)

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{StringConvertCV}}

******************************************************* </lalLaTeX> */

#include <ctype.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/StringInput.h>

NRCSID( STRINGCONVERTC, "$Id$" );

/* Extremal integer values, all expressed as unsigned long long. */
#define LAL_UINT8_MAX   LAL_UINT8_C(18446744073709551615)
#define LAL_UINT4_MAX   LAL_UINT8_C(4294967295)
#define LAL_UINT2_MAX   LAL_UINT8_C(65535)
#define LAL_INT8_MAX    LAL_UINT8_C(9223372036854775807)
#define LAL_INT4_MAX    LAL_UINT8_C(2147483647)
#define LAL_INT2_MAX    LAL_UINT8_C(32767)
#define LAL_INT8_ABSMIN LAL_UINT8_C(9223372036854775808)
#define LAL_INT4_ABSMIN LAL_UINT8_C(2147483648)
#define LAL_INT2_ABSMIN LAL_UINT8_C(32768)

/* Maximum number of digits we ever need to parse for an integer. */
#define LAL_UINT8_MAXDIGITS (20)

/* Internal function to parse integer strings into UINT8 magnitudes
   plus a sign. */
static UINT8
LALStringToU8AndSign( INT2 *sign, const CHAR *string, CHAR **endptr )
{
  union { char *s; const char *cs; } bad;/* there is a REASON for warnings... */
  const CHAR *here = string;         /* current position in string */
  CHAR c;                            /* current character in string */
  UINT4 n = LAL_UINT8_MAXDIGITS - 1; /* number of worry-free digits */
  UINT8 value;                       /* current converted value */

  /* Skip leading space, and read sign character, if any. */
  *sign = 1;
  while ( isspace( *here ) )
    here++;
  if ( *here == '+' )
    here++;
  else if ( *here == '-' ) {
    *sign = -1;
    here++;
  }

  /* Read first digit.  Abort if it's not a digit. */
  if ( isdigit( (int)( c = *here ) ) ) {
    value = (UINT8)( c - '0' );
    here++;
  } else {
    bad.cs = string; /* ... and this avoids the warnings... BAD! */
    *endptr = bad.s;
    return 0;
  }

  /* Otherwise, start reading number.  Stop if we get close to
     overflowing. */
  while ( isdigit( (int)( c = *here ) ) && --n ) {
    value *= LAL_INT8_C(10);
    value += (UINT8)( c - '0' );
    here++;
  }

  /* Proceed with caution near overflow.  At this point, if n==0, then
     c = *here is the (LAL_UINT8_MAXDIGITS)th digit read, but value
     does not yet incorporate it. */
  if ( !n ) {
    here++;
    if ( isdigit( (int)( *here ) ) ) {
      value = LAL_UINT8_MAX;
      do
	here++;
      while ( isdigit( (int)( *here ) ) );
    } else if ( value > LAL_UINT8_MAX/LAL_INT8_C(10) ) {
      value = LAL_UINT8_MAX;
    } else {
      UINT8 increment = (UINT8)( c - '0' );
      value *= 10;
      if ( value > LAL_UINT8_MAX - increment )
	value = LAL_UINT8_MAX;
      else
	value += increment;
    }
  }

  /* Return appropriate values. */
  bad.cs = here; /* ... and this avoids the warnings... BAD! */
  *endptr = bad.s;
  return value;
}


/* <lalVerbatim file="StringConvertCP"> */
void
LALStringToU2( LALStatus *stat, UINT2 *value, const CHAR *string, CHAR **endptr )
{ /* </lalVerbatim> */
  UINT8 absValue; /* magnitude of parsed number */
  INT2 sign;      /* sign of parsed number */
  CHAR *end;      /* substring following parsed number */

  INITSTATUS( stat, "LALStringToU2", STRINGCONVERTC );

  /* Check for valid input arguments. */
  ASSERT( value, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL );
  ASSERT( string, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL );

  /* Parse string.  Return if nothing was parsed. */
  absValue = LALStringToU8AndSign( &sign, string, &end );
  if ( string == end ) {
    if ( endptr )
      *endptr = end;
    RETURN( stat );
  }

  /* Cap (if necessary), cast, and return. */
  if ( absValue > LAL_UINT2_MAX )
    *value = (UINT2)( LAL_UINT2_MAX );
  else
    *value = (UINT2)( absValue );
  if ( endptr )
    *endptr = end;
  RETURN( stat );
}


/* <lalVerbatim file="StringConvertCP"> */
void
LALStringToU4( LALStatus *stat, UINT4 *value, const CHAR *string, CHAR **endptr )
{ /* </lalVerbatim> */
  UINT8 absValue; /* magnitude of parsed number */
  INT2 sign;      /* sign of parsed number */
  CHAR *end;      /* substring following parsed number */

  INITSTATUS( stat, "LALStringToU4", STRINGCONVERTC );

  /* Check for valid input arguments. */
  ASSERT( value, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL );
  ASSERT( string, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL );

  /* Parse string.  Return if nothing was parsed. */
  absValue = LALStringToU8AndSign( &sign, string, &end );
  if ( string == end ) {
    if ( endptr )
      *endptr = end;
    RETURN( stat );
  }

  /* Cap (if necessary), cast, and return. */
  if ( absValue > LAL_UINT4_MAX )
    *value = (UINT4)( LAL_UINT4_MAX );
  else
    *value = (UINT4)( absValue );
  if ( endptr )
    *endptr = end;
  RETURN( stat );
}


/* <lalVerbatim file="StringConvertCP"> */
void
LALStringToU8( LALStatus *stat, UINT8 *value, const CHAR *string, CHAR **endptr )
{ /* </lalVerbatim> */
  UINT8 absValue; /* magnitude of parsed number */
  INT2 sign;      /* sign of parsed number */
  CHAR *end;      /* substring following parsed number */

  INITSTATUS( stat, "LALStringToU8", STRINGCONVERTC );

  /* Check for valid input arguments. */
  ASSERT( value, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL );
  ASSERT( string, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL );

  /* Parse string.  Return if nothing was parsed. */
  absValue = LALStringToU8AndSign( &sign, string, &end );
  if ( string == end ) {
    if ( endptr )
      *endptr = end;
    RETURN( stat );
  }

  /* Set values and return. */
  *value = absValue;
  if ( endptr )
    *endptr = end;
  RETURN( stat );
}


/* <lalVerbatim file="StringConvertCP"> */
void
LALStringToI2( LALStatus *stat, INT2 *value, const CHAR *string, CHAR **endptr )
{ /* </lalVerbatim> */
  UINT8 absValue; /* magnitude of parsed number */
  INT2 sign;      /* sign of parsed number */
  CHAR *end;      /* substring following parsed number */

  INITSTATUS( stat, "LALStringToI2", STRINGCONVERTC );

  /* Check for valid input arguments. */
  ASSERT( value, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL );
  ASSERT( string, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL );

  /* Parse string.  Return if nothing was parsed. */
  absValue = LALStringToU8AndSign( &sign, string, &end );
  if ( string == end ) {
    if ( endptr )
      *endptr = end;
    RETURN( stat );
  }

  /* Cap (if necessary), cast, and return. */
  if ( sign > 0 ) {
    if ( absValue > LAL_INT2_MAX )
      *value = (INT2)( LAL_INT2_MAX );
    else
      *value = (INT2)( absValue );
  } else {
    if ( absValue > LAL_INT2_ABSMIN )
      *value = (INT2)( -LAL_INT2_ABSMIN );
    else
      *value = (INT2)( -absValue );
  }
  if ( endptr )
    *endptr = end;
  RETURN( stat );
}


/* <lalVerbatim file="StringConvertCP"> */
void
LALStringToI4( LALStatus *stat, INT4 *value, const CHAR *string, CHAR **endptr )
{ /* </lalVerbatim> */
  UINT8 absValue; /* magnitude of parsed number */
  INT2 sign;      /* sign of parsed number */
  CHAR *end;      /* substring following parsed number */

  INITSTATUS( stat, "LALStringToI4", STRINGCONVERTC );

  /* Check for valid input arguments. */
  ASSERT( value, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL );
  ASSERT( string, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL );

  /* Parse string.  Return if nothing was parsed. */
  absValue = LALStringToU8AndSign( &sign, string, &end );
  if ( string == end ) {
    if ( endptr )
      *endptr = end;
    RETURN( stat );
  }

  /* Cap (if necessary), cast, and return. */
  if ( sign > 0 ) {
    if ( absValue > LAL_INT4_MAX )
      *value = (INT4)( LAL_INT4_MAX );
    else
      *value = (INT4)( absValue );
  } else {
    if ( absValue > LAL_INT4_ABSMIN )
      *value = (INT4)( -LAL_INT4_ABSMIN );
    else
      *value = (INT4)( -absValue );
  }
  if ( endptr )
    *endptr = end;
  RETURN( stat );
}


/* <lalVerbatim file="StringConvertCP"> */
void
LALStringToI8( LALStatus *stat, INT8 *value, const CHAR *string, CHAR **endptr )
{ /* </lalVerbatim> */
  UINT8 absValue; /* magnitude of parsed number */
  INT2 sign;      /* sign of parsed number */
  CHAR *end;      /* substring following parsed number */

  INITSTATUS( stat, "LALStringToI8", STRINGCONVERTC );

  /* Check for valid input arguments. */
  ASSERT( value, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL );
  ASSERT( string, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL );

  /* Parse string.  Return if nothing was parsed. */
  absValue = LALStringToU8AndSign( &sign, string, &end );
  if ( string == end ) {
    if ( endptr )
      *endptr = end;
    RETURN( stat );
  }

  /* Cap (if necessary), cast, and return. */
  if ( sign > 0 ) {
    if ( absValue > LAL_INT8_MAX )
      *value = (INT8)( LAL_INT8_MAX );
    else
      *value = (INT8)( absValue );
  } else {
    if ( absValue > LAL_INT8_ABSMIN )
      *value = (INT8)( -LAL_INT8_ABSMIN );
    else
      *value = (INT8)( -absValue );
  }
  if ( endptr )
    *endptr = end;
  RETURN( stat );
}


/* <lalVerbatim file="StringConvertCP"> */
void
LALStringToS( LALStatus *stat, REAL4 *value, const CHAR *string, CHAR **endptr )
{ /* </lalVerbatim> */
  REAL8 myValue; /* internal representation of value */
  CHAR *end;     /* substring following parsed number */

  INITSTATUS( stat, "LALStringToS", STRINGCONVERTC );

  /* Check for valid input arguments. */
  ASSERT( value, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL );
  ASSERT( string, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL );

  /* Parse string.  Return if nothing was parsed. */
  myValue = strtod( string, &end );
  if ( string == end ) {
    if ( endptr )
      *endptr = end;
    RETURN( stat );
  }

  /* Cap (if necessary), cast, and return. */
  if ( myValue > LAL_REAL4_MAX )
    *value = (REAL4)( LAL_REAL4_MAX );
  else if ( myValue < -LAL_REAL4_MAX )
    *value = (REAL4)( -LAL_REAL4_MAX );
  else
    *value = (REAL4)( myValue );
  if ( endptr )
    *endptr = end;
  RETURN( stat );
}


/* <lalVerbatim file="StringConvertCP"> */
void
LALStringToD( LALStatus *stat, REAL8 *value, const CHAR *string, CHAR **endptr )
{ /* </lalVerbatim> */
  REAL8 myValue; /* internal representation of value */
  CHAR *end;     /* substring following parsed number */

  INITSTATUS( stat, "LALStringToD", STRINGCONVERTC );

  /* Check for valid input arguments. */
  ASSERT( value, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL );
  ASSERT( string, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL );

  /* Parse string.  Return if nothing was parsed. */
  myValue = strtod( string, &end );
  if ( string == end ) {
    if ( endptr )
      *endptr = end;
    RETURN( stat );
  }

  /* Set values and return. */
  if ( myValue > LAL_REAL8_MAX )
    *value = LAL_REAL8_MAX;
  else if ( myValue < -LAL_REAL8_MAX )
    *value = -LAL_REAL8_MAX;
  else
    *value = myValue;
  if ( endptr )
    *endptr = end;
  RETURN( stat );
}


/* <lalVerbatim file="StringConvertCP"> */
void
LALStringToC( LALStatus *stat, COMPLEX8 *value, const CHAR *string, CHAR **endptr )
{ /* </lalVerbatim> */
  REAL4 re, im; /* real and imaginary parts */
  CHAR *end;    /* substring following parsed numbers */

  INITSTATUS( stat, "LALStringToC", STRINGCONVERTC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( value, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL );
  ASSERT( string, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL );

  /* Parse string.  Return if nothing was parsed. */
  TRY( LALStringToS( stat->statusPtr, &re, string, &end ), stat );
  TRY( LALStringToS( stat->statusPtr, &im, end, &end ), stat );
  if ( string == end ) {
    if ( endptr )
      *endptr = end;
    DETATCHSTATUSPTR( stat );
    RETURN( stat );
  }

  /* Set values and return. */
  value->re = re;
  value->im = im;
  if ( endptr )
    *endptr = end;
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


/* <lalVerbatim file="StringConvertCP"> */
void
LALStringToZ( LALStatus *stat, COMPLEX16 *value, const CHAR *string, CHAR **endptr )
{ /* </lalVerbatim> */
  REAL8 re, im; /* real and imaginary parts */
  CHAR *end;    /* substring following parsed numbers */

  INITSTATUS( stat, "LALStringToZ", STRINGCONVERTC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( value, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL );
  ASSERT( string, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL );

  /* Parse string.  Return if nothing was parsed. */
  TRY( LALStringToD( stat->statusPtr, &re, string, &end ), stat );
  TRY( LALStringToD( stat->statusPtr, &im, end, &end ), stat );
  if ( string == end ) {
    if ( endptr )
      *endptr = end;
    DETATCHSTATUSPTR( stat );
    RETURN( stat );
  }

  /* Set values and return. */
  value->re = re;
  value->im = im;
  if ( endptr )
    *endptr = end;
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
