/******************************** <lalVerbatim file="StringConvertCV">
Authors: Creighton, T. D.
         Shawhan, P. S.   (LALStringToGPS)
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
\idx{LALStringToGPS()}

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

GPS conversion is similar to floating-point conversion, but the result
is stored in a \verb@LIGOTimeGPS@ structure as two integer values
representing seconds and nanoseconds.  The \verb@LALStringToGPS@
function does {\it not} convert the string to an intermediate
\verb@REAL8@ value, but parses the string specially to retain the
full precision of the string representation to the nearest nanosecond.

\subsubsection*{Algorithm}

These functions (other than \verb@LALStringToGPS@)
emulate the standard C functions \verb@strtol()@,
\verb@strtoul()@, and \verb@strtod()@, except that they follow LAL
calling conventions and return values of the appropriate LAL
datatypes.  For integer conversion, only base-ten (decimal)
representations are supported.  Otherwise, the valid format is as for
the corresponding C functions, which we summarize below:

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
Usually this means that the second substring will contain at least one
leading whitespace character, though this is not strictly necessary.
Overflow or underflow is dealt with as above.

A string to be converted to a GPS time can have the format of an integer
or a floating-point number, as described above.  The optional exponent
in the floating-point form is supported, and both positive and negative
GPS times are permitted.  If the result would overflow the
\verb@LIGOTimeGPS@ representation (too far in the future), then the
\verb@gpsSeconds@ and \verb@gpsNanoSeconds@ fields are set to
\verb@LAL_INT4_MAX@ and $999999999$, respectively.
For an underflow (too far in the past), the fields
are set to \verb@LAL_INT4_MIN@ and $-999999999$~.

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


/* <lalVerbatim file="StringConvertCP"> */
void
LALStringToGPS( LALStatus *stat, LIGOTimeGPS *value, const CHAR *string, CHAR **endptr )
{ /* </lalVerbatim> */

  const CHAR *here = string;   /* current position in string */
  INT4 signval;       /* sign of value (+1 or -1) */
  CHAR mantissa[64];  /* local string to store mantissa digits */
  INT4 mdigits;       /* number of digits in mantissa */
  INT4 dppos;         /* position of decimal point in mantissa, i.e. the
                         number of mantissa digits preceding the decimal point
                         (initially -1 if no decimal point in input.) */
  CHAR intstring[16]; /* local string to store integer part of time */
  const CHAR *ehere;  /* string pointer for parsing exponent */
  INT4 exponent;      /* exponent given in string */
  INT2 esignval;      /* sign of exponent value (+1 or -1) */
  UINT8 absValue;     /* magnitude of parsed number */
  CHAR *eend;         /* substring following parsed number */
  INT4 nanosecSet;    /* flag to indicate if nanoseconds field has been set */
  CHAR *nptr;         /* pointer to where nanoseconds begin in mantissa */
  INT4 idigit;

  INITSTATUS( stat, "LALStringToGPS", STRINGCONVERTC );
  ATTATCHSTATUSPTR( stat );

  /* Check for valid input arguments. */
  ASSERT( value, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL );
  ASSERT( string, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL );

  /* Skip leading space, and read sign character, if any. */
  signval = 1;
  while ( isspace( *here ) )
    here++;
  if ( *here == '+' )
    here++;
  else if ( *here == '-' ) {
    signval = -1;
    here++;
  }

  /* Copy the mantissa into a local string, keeping track of the
     location of the decimal point */
  mdigits = 0;
  dppos = -1;
  while ( (*here >= '0' && *here <= '9') || (*here == '.') ) {

    if ( *here == '.' ) {
      /* This is a decimal point */
      if ( dppos >= 0 ) {
	/* This is second decimal point encountered, so parsing must stop */
	break;
      } else {
	/* Record the position of the decimal point */
	dppos = (INT4) mdigits;
      }

    } else {
      /* This is a digit.  Append it to the local mantissa string, unless
	 the mantissa string is already full, in which case ignore it */
      if ( (UINT4) mdigits < sizeof(mantissa)-1 ) {
	mantissa[mdigits] = *here;
	mdigits++;
      }
    }

    here++;

  }

  /* If there is no mantissa, then return without consuming any characters
     and without modifying 'value' */
  if ( mdigits == 0 ) {
    if ( endptr )
      *endptr = string;
    DETATCHSTATUSPTR( stat );
    RETURN( stat );
  }

  /* Null-terminate the mantissa string */
  mantissa[mdigits] = '\0';

  /* If there was no explicit decimal point, then it is implicitly
     after all of the mantissa digits */
  if ( dppos == -1 ) { dppos = (INT4) mdigits; }

  /* Read the exponent, if present */
  exponent = 0;
  if ( *here == 'E' || *here == 'e' ) {
    /* So far, this looks like an exponent.  Set working pointer. */
    ehere = here + 1;

    /* Parse the exponent value */
    absValue = LALStringToU8AndSign( &esignval, ehere, &eend );
    if ( eend == ehere ) {
      /* Nothing was parsed, so this isn't a valid exponent.  Leave
	 things as they are, with 'here' pointing to the 'E' or 'e'
	 that we thought introduced an exponent. */
    } else {
      /* We successfully parsed the exponent */
      exponent = (INT4) ( esignval * absValue );
      /* Update the 'here' pointer to just after the exponent */
      here = eend;
    }

  }

  /* The exponent simply modifies the decimal point position */
  dppos += exponent;

  /* OK, now we have the sign ('signval'), mantissa string, and
     decimal point location, and the 'here' pointer points to the
     first character after the part of the string we parsed. */

  nanosecSet = 0;

  if ( dppos > 10 ) {
    /* This is an overflow (positive) or underflow (negative) */
    if ( signval == 1 ) {
      value->gpsSeconds = (INT4) LAL_INT4_MAX;
      value->gpsNanoSeconds = 999999999;
    } else {
      value->gpsSeconds = (INT4)( -LAL_INT4_ABSMIN );
      value->gpsNanoSeconds = -999999999;
    }

  } else if ( dppos < -9 ) {
    /* The time is effectively zero */
    value->gpsSeconds = 0;
    value->gpsNanoSeconds = 0;

  } else {

    /* See whether there is an integer part... */
    if ( dppos > 0 ) {

      /* Pick out the integer part */
      memcpy( intstring, mantissa, dppos );
      intstring[dppos] = '\0';
      absValue = LALStringToU8AndSign( &esignval, intstring, &eend );
      /* We ignore the 'esignval' and 'eend' variables */

      /* If the mantissa had too few digits, need to multiply by tens */
      for ( idigit=mdigits; idigit<dppos; idigit++ ) { absValue *= 10; }

      /* Cap (if necessary) and cast */
      if ( signval > 0 ) {
	if ( absValue > LAL_INT4_MAX ) {
	  value->gpsSeconds = (INT4)( LAL_INT4_MAX );
          value->gpsNanoSeconds = 999999999;
	  nanosecSet = 1;
	} else
	  value->gpsSeconds = (INT4)( absValue );
      } else {
	if ( absValue > LAL_INT4_ABSMIN ) {
	  value->gpsSeconds = (INT4)( -LAL_INT4_ABSMIN );
          value->gpsNanoSeconds = -999999999;
	  nanosecSet = 1;
	} else
	  value->gpsSeconds = (INT4)( -absValue );
      }

    } else {
      /* There is no integer part */
      value->gpsSeconds = 0;
    }

    /* Finally, set the nanoseconds field (if not already set) */
    if ( ! nanosecSet ) {

      if ( dppos >= mdigits ) {
	value->gpsNanoSeconds = 0;

      } else {

	if ( dppos >= 0 ) {
	  nptr = mantissa + dppos;
	  /* Revise dppos to refer to the substring pointed to by 'nptr' */
	  dppos = 0;
	} else {
	  nptr = mantissa;
	  /* Leave dppos with its current negative value */
	}
	/* Revise mdigits to be the number of digits in the nanoseconds part */
	mdigits = strlen(nptr);
	if ( mdigits == 0 ) {
	  /* The nanoseconds field is absent, thus equal to zero */
	  absValue = 0;
	} else if ( mdigits >= 9+dppos ) {
	  memcpy( intstring, nptr, 9+dppos );
	  intstring[9+dppos] = '\0';
	  absValue = strtol(intstring,NULL,10);
	  /* If there is another digit, use it to round */
	  if ( mdigits >= 10+dppos ) {
	    if ( *(nptr+9+dppos) >= '5' ) {
	      absValue++;
	    }
	  }
	} else {
	  /* Digits are not given all the way to nanoseconds, so have to
	     multiply by a power of 10 */
	  absValue = strtol(nptr,NULL,10);
	  while ( mdigits < 9+dppos ) {
	    absValue *= 10;
	    mdigits++;
	  }
	}

	value->gpsNanoSeconds = (INT4) ( signval * absValue );

	/* Check for wraparound due to rounding */
	if ( value->gpsNanoSeconds >= 1000000000 ) {
	  value->gpsNanoSeconds -= 1000000000;
	  value->gpsSeconds += 1;
	}
	if ( value->gpsNanoSeconds <= -1000000000 ) {
	  value->gpsNanoSeconds += 1000000000;
	  value->gpsSeconds -= 1;
	}

      }

    }

  }

  /* Set end pointer (if passed) and return. */
  if ( endptr )
    *endptr = here;
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
