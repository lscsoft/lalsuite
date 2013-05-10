/*
*  Copyright (C) 2007 Jolien Creighton, Peter Shawhan
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/**
   \file
   \ingroup StringInput_h
   \authors Creighton, T. D.
   \authors Shawhan, P. S.

   \brief Converts a string into a numerical value.

\heading{Description}

These routines parse the string <tt>*string</tt> and compute a numerical
value <tt>*value</tt> of the appropriate datatype.  If
\c endptr\f$\neq\f$\c NULL, then after conversion <tt>*endptr</tt>
will point to the character after the last character used in the
conversion.  The routine will always return without error if the
arguments are valid, regardless of the contents of \c string;
failure to parse a number is indicated by <tt>*endptr</tt> being set
equal to \c string (provided \c endptr\f$\neq\f$\c NULL).

For integer or floating-point conversion to occur, \c string must
consist of zero or more whitespace characters followed by a number in
any standard base-ten integer or floating-point representation
(described in detail below); the conversion will stop as soon as the
routine encounters a character that is not part of the number.  For
instance, parsing the string <tt>"   123bad"</tt> will return
<tt>*value</tt>=123, and <tt>*endptr</tt> will point to the substring
<tt>"bad"</tt>; it is up to the calling routine to determine whether
this is an acceptable input.  By contrast, parsing the string
<tt>" bad"</tt> will leave <tt>*value</tt> unchanged and <tt>*endptr</tt>
pointing to the start of the original string.  In general, if the
routine returns with <tt>*endptr</tt>\f$\neq\f$\c string and
<tt>**endptr</tt> is a whitespace or <tt>'\0'</tt> character, then the
format was unambiguously acceptable.

Complex conversion is essentially equivalent to performing
floating-point conversion on \c string to get the real part, and
again on <tt>*endptr</tt> to get the imaginary part.  Normally this
means that an acceptable format is two floating-point representations
separated by (and possibly preceded with) whitespace, although the
intervening whitespace can be omitted if the separation between the
two numbers is unambiguous.  Thus the string <tt>"-1.0+3.5"</tt> will be
unambiguously read as \f$-1.0+3.5i\f$, but <tt>"-1.03.5"</tt> will be read
as \f$-1.03+0.5i\f$ (since the conversion of the real part stops at the
second <tt>'.'</tt> character), which may or may not be the intended
conversion.

GPS conversion is similar to floating-point conversion, but the result
is stored in a \c LIGOTimeGPS structure as two integer values
representing seconds and nanoseconds.  The \c LALStringToGPS
function does \e not convert the string to an intermediate
\c REAL8 value, but parses the string specially to retain the
full precision of the string representation to the nearest nanosecond.

\heading{Algorithm}

These functions (other than \c LALStringToGPS())
emulate the standard C functions <tt>strtol()</tt>,
<tt>strtoul()</tt>, and <tt>strtod()</tt>, except that they follow LAL
calling conventions and return values of the appropriate LAL
datatypes.  For integer conversion, only base-ten (decimal)
representations are supported.  Otherwise, the valid format is as for
the corresponding C functions, which we summarize below:

A string to be converted to an \c INT\f$n\f$ (where \f$n\f$=2, 4, or 8)
consists of zero or more whitespace characters as determined by
<tt>isspace()</tt>, followed optionally by a single <tt>'+'</tt> or
<tt>'-'</tt> character, followed by one or more decimal digits;
conversion stops at the first non-digit character after this.  If the
result would overflow or underflow the \c INT\f$n\f$ representation,
then the value is set to \c LAL_INT\f$n\f$\c _MAX\f$=2^{8n-1}-1\f$ or
\c LAL_INT\f$n\f$\c _MIN\f$=-2^{8n-1}\f$, respectively.

A string to be converted to a \c UINT\f$n\f$ follows the same format,
except that a leading negative sign character <tt>'-'</tt> is
\e ignored (the routine will compute the magnitude of the result),
and the return value is capped at
\c LAL_UINT\f$n\f$\c _MAX\f$=2^{8n}-1\f$.

A string to be converted to a floating-point number (\c REAL4 or
\c REAL8) consists of zero or more whitespace characters as
determined by <tt>isspace()</tt>, followed optionally by a single
<tt>'+'</tt> or <tt>'-'</tt> character, followed by a sequence of one or
more decimal digits optionally containing a decimal point <tt>'.'</tt>,
optionally followed by an exponent.  An exponent consists of a single
<tt>'E'</tt> or <tt>'e'</tt> character, followed optionally by a single
<tt>'+'</tt> or <tt>'-'</tt> character, followed by a sequence of one or
more decimal digits.  If the converted value would overflow,
\f$\pm\f$\c LAL_REAL\f$n\f$\c _MAX is returned, as appropriate.  If
the value would underflow, 0 is returned.

A string to be converted to a complex number (\c COMPLEX8 or
\c COMPLEX16) consists of two floating-point format substrings
concatenated together, where the first character of the second
substring cannot be interpreted as a continuation of the first number.
Usually this means that the second substring will contain at least one
leading whitespace character, though this is not strictly necessary.
Overflow or underflow is dealt with as above.

A string to be converted to a GPS time can have the format of an integer
or a floating-point number, as described above.  The optional exponent
in the floating-point form is supported, and both positive and negative
GPS times are permitted.  If the result would overflow the
\c LIGOTimeGPS representation (too far in the future), then the
\c gpsSeconds and \c gpsNanoSeconds fields are set to
\c LAL_INT4_MAX and \f$999999999\f$, respectively.
For an underflow (too far in the past), the fields
are set to \c LAL_INT4_MIN and \f$0\f$.

Internally, the floating-point conversion routines call
<tt>strtod()</tt>, then cap and cast the result as necessary.  The
complex conversion routines simply call their floating-point
counterpart twice.  The integer routines call an internal function
<tt>LALStringToU8AndSign()</tt>, which does what you would expect, then
cap and cast the result as necessary.  (The C routines <tt>strtol()</tt>
and <tt>strtol()</tt> are not used as they are not guaranteed to have
8-byte precision.)

*/

/** \cond DONT_DOXYGEN */
#include <complex.h>
#include <ctype.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/StringInput.h>

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
LALStringToU8AndSign(INT2 * sign, const CHAR * string, CHAR ** endptr)
{
    union {
        char *s;
        const char *cs;
    } bad;      /* there is a REASON for warnings... */
    const CHAR *here = string;  /* current position in string */
    CHAR c;     /* current character in string */
    UINT4 n = LAL_UINT8_MAXDIGITS - 1;  /* number of worry-free digits */
    UINT8 value;        /* current converted value */

    /* Skip leading space, and read sign character, if any. */
    *sign = 1;
    while (isspace(*here))
        here++;
    if (*here == '+')
        here++;
    else if (*here == '-') {
        *sign = -1;
        here++;
    }

    /* Read first digit.  Abort if it's not a digit. */
    if (isdigit((int) (c = *here))) {
        value = (UINT8) (c - '0');
        here++;
    } else {
        bad.cs = string;        /* ... and this avoids the warnings... BAD! */
        *endptr = bad.s;
        return 0;
    }

    /* Otherwise, start reading number.  Stop if we get close to
       overflowing. */
    while (isdigit((int) (c = *here)) && --n) {
        value *= LAL_INT8_C(10);
        value += (UINT8) (c - '0');
        here++;
    }

    /* Proceed with caution near overflow.  At this point, if n==0, then
       c = *here is the (LAL_UINT8_MAXDIGITS)th digit read, but value
       does not yet incorporate it. */
    if (!n) {
        here++;
        if (isdigit((int) (*here))) {
            value = LAL_UINT8_MAX;
            do
                here++;
            while (isdigit((int) (*here)));
        } else if (value > LAL_UINT8_MAX / LAL_INT8_C(10)) {
            value = LAL_UINT8_MAX;
        } else {
            UINT8 increment = (UINT8) (c - '0');
            value *= 10;
            if (value > LAL_UINT8_MAX - increment)
                value = LAL_UINT8_MAX;
            else
                value += increment;
        }
    }

    /* Return appropriate values. */
    bad.cs = here;      /* ... and this avoids the warnings... BAD! */
    *endptr = bad.s;
    return value;
}

/** \endcond */

void
LALStringToU2(LALStatus * stat, UINT2 * value, const CHAR * string,
              CHAR ** endptr)
{
    UINT8 absValue;     /* magnitude of parsed number */
    INT2 sign;  /* sign of parsed number */
    CHAR *end;  /* substring following parsed number */

    INITSTATUS(stat);

    /* Check for valid input arguments. */
    ASSERT(value, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL);
    ASSERT(string, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL);

    /* Parse string.  Return if nothing was parsed. */
    absValue = LALStringToU8AndSign(&sign, string, &end);
    if (string == end) {
        if (endptr)
            *endptr = end;
        RETURN(stat);
    }

    /* Cap (if necessary), cast, and return. */
    if (absValue > LAL_UINT2_MAX)
        *value = (UINT2) (LAL_UINT2_MAX);
    else
        *value = (UINT2) (absValue);
    if (endptr)
        *endptr = end;
    RETURN(stat);
}



void
LALStringToU4(LALStatus * stat, UINT4 * value, const CHAR * string,
              CHAR ** endptr)
{
    UINT8 absValue;     /* magnitude of parsed number */
    INT2 sign;  /* sign of parsed number */
    CHAR *end;  /* substring following parsed number */

    INITSTATUS(stat);

    /* Check for valid input arguments. */
    ASSERT(value, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL);
    ASSERT(string, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL);

    /* Parse string.  Return if nothing was parsed. */
    absValue = LALStringToU8AndSign(&sign, string, &end);
    if (string == end) {
        if (endptr)
            *endptr = end;
        RETURN(stat);
    }

    /* Cap (if necessary), cast, and return. */
    if (absValue > LAL_UINT4_MAX)
        *value = (UINT4) (LAL_UINT4_MAX);
    else
        *value = (UINT4) (absValue);
    if (endptr)
        *endptr = end;
    RETURN(stat);
}



void
LALStringToU8(LALStatus * stat, UINT8 * value, const CHAR * string,
              CHAR ** endptr)
{
    UINT8 absValue;     /* magnitude of parsed number */
    INT2 sign;  /* sign of parsed number */
    CHAR *end;  /* substring following parsed number */

    INITSTATUS(stat);

    /* Check for valid input arguments. */
    ASSERT(value, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL);
    ASSERT(string, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL);

    /* Parse string.  Return if nothing was parsed. */
    absValue = LALStringToU8AndSign(&sign, string, &end);
    if (string == end) {
        if (endptr)
            *endptr = end;
        RETURN(stat);
    }

    /* Set values and return. */
    *value = absValue;
    if (endptr)
        *endptr = end;
    RETURN(stat);
}



void
LALStringToI2(LALStatus * stat, INT2 * value, const CHAR * string,
              CHAR ** endptr)
{
    UINT8 absValue;     /* magnitude of parsed number */
    INT2 sign;  /* sign of parsed number */
    CHAR *end;  /* substring following parsed number */

    INITSTATUS(stat);

    /* Check for valid input arguments. */
    ASSERT(value, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL);
    ASSERT(string, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL);

    /* Parse string.  Return if nothing was parsed. */
    absValue = LALStringToU8AndSign(&sign, string, &end);
    if (string == end) {
        if (endptr)
            *endptr = end;
        RETURN(stat);
    }

    /* Cap (if necessary), cast, and return. */
    if (sign > 0) {
        if (absValue > LAL_INT2_MAX)
            *value = (INT2) (LAL_INT2_MAX);
        else
            *value = (INT2) (absValue);
    } else {
        if (absValue > LAL_INT2_ABSMIN)
            *value = (INT2) (-LAL_INT2_ABSMIN);
        else
            *value = (INT2) (-absValue);
    }
    if (endptr)
        *endptr = end;
    RETURN(stat);
}



void
LALStringToI4(LALStatus * stat, INT4 * value, const CHAR * string,
              CHAR ** endptr)
{
    UINT8 absValue;     /* magnitude of parsed number */
    INT2 sign;  /* sign of parsed number */
    CHAR *end;  /* substring following parsed number */

    INITSTATUS(stat);

    /* Check for valid input arguments. */
    ASSERT(value, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL);
    ASSERT(string, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL);

    /* Parse string.  Return if nothing was parsed. */
    absValue = LALStringToU8AndSign(&sign, string, &end);
    if (string == end) {
        if (endptr)
            *endptr = end;
        RETURN(stat);
    }

    /* Cap (if necessary), cast, and return. */
    if (sign > 0) {
        if (absValue > LAL_INT4_MAX)
            *value = (INT4) (LAL_INT4_MAX);
        else
            *value = (INT4) (absValue);
    } else {
        if (absValue > LAL_INT4_ABSMIN)
            *value = (INT4) (-LAL_INT4_ABSMIN);
        else
            *value = (INT4) (-absValue);
    }
    if (endptr)
        *endptr = end;
    RETURN(stat);
}



void
LALStringToI8(LALStatus * stat, INT8 * value, const CHAR * string,
              CHAR ** endptr)
{
    UINT8 absValue;     /* magnitude of parsed number */
    INT2 sign;  /* sign of parsed number */
    CHAR *end;  /* substring following parsed number */

    INITSTATUS(stat);

    /* Check for valid input arguments. */
    ASSERT(value, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL);
    ASSERT(string, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL);

    /* Parse string.  Return if nothing was parsed. */
    absValue = LALStringToU8AndSign(&sign, string, &end);
    if (string == end) {
        if (endptr)
            *endptr = end;
        RETURN(stat);
    }

    /* Cap (if necessary), cast, and return. */
    if (sign > 0) {
        if (absValue > LAL_INT8_MAX)
            *value = (INT8) (LAL_INT8_MAX);
        else
            *value = (INT8) (absValue);
    } else {
        if (absValue > LAL_INT8_ABSMIN)
            *value = (INT8) (-LAL_INT8_ABSMIN);
        else
            *value = (INT8) (-absValue);
    }
    if (endptr)
        *endptr = end;
    RETURN(stat);
}



void
LALStringToS(LALStatus * stat, REAL4 * value, const CHAR * string,
             CHAR ** endptr)
{
    REAL8 myValue;      /* internal representation of value */
    CHAR *end;  /* substring following parsed number */

    INITSTATUS(stat);

    /* Check for valid input arguments. */
    ASSERT(value, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL);
    ASSERT(string, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL);

    /* Parse string.  Return if nothing was parsed. */
    myValue = strtod(string, &end);
    if (string == end) {
        if (endptr)
            *endptr = end;
        RETURN(stat);
    }

    /* Cap (if necessary), cast, and return. */
    if (myValue > LAL_REAL4_MAX)
        *value = (REAL4) (LAL_REAL4_MAX);
    else if (myValue < -LAL_REAL4_MAX)
        *value = (REAL4) (-LAL_REAL4_MAX);
    else
        *value = (REAL4) (myValue);
    if (endptr)
        *endptr = end;
    RETURN(stat);
}



void
LALStringToD(LALStatus * stat, REAL8 * value, const CHAR * string,
             CHAR ** endptr)
{
    REAL8 myValue;      /* internal representation of value */
    CHAR *end;  /* substring following parsed number */

    INITSTATUS(stat);

    /* Check for valid input arguments. */
    ASSERT(value, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL);
    ASSERT(string, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL);

    /* Parse string.  Return if nothing was parsed. */
    myValue = strtod(string, &end);
    if (string == end) {
        if (endptr)
            *endptr = end;
        RETURN(stat);
    }

    /* Set values and return. */
    if (myValue > LAL_REAL8_MAX)
        *value = LAL_REAL8_MAX;
    else if (myValue < -LAL_REAL8_MAX)
        *value = -LAL_REAL8_MAX;
    else
        *value = myValue;
    if (endptr)
        *endptr = end;
    RETURN(stat);
}



void
LALStringToC(LALStatus * stat, COMPLEX8 * value, const CHAR * string,
             CHAR ** endptr)
{
    REAL4 re, im;       /* real and imaginary parts */
    CHAR *end;  /* substring following parsed numbers */

    INITSTATUS(stat);
    ATTATCHSTATUSPTR(stat);

    /* Check for valid input arguments. */
    ASSERT(value, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL);
    ASSERT(string, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL);

    /* Parse string.  Return if nothing was parsed. */
    TRY(LALStringToS(stat->statusPtr, &re, string, &end), stat);
    TRY(LALStringToS(stat->statusPtr, &im, end, &end), stat);
    if (string == end) {
        if (endptr)
            *endptr = end;
        DETATCHSTATUSPTR(stat);
        RETURN(stat);
    }

    /* Set values and return. */
    *value = re;
    *value += im * I;
    if (endptr)
        *endptr = end;
    DETATCHSTATUSPTR(stat);
    RETURN(stat);
}



void
LALStringToZ(LALStatus * stat, COMPLEX16 * value, const CHAR * string,
             CHAR ** endptr)
{
    REAL8 re, im;       /* real and imaginary parts */
    CHAR *end;  /* substring following parsed numbers */

    INITSTATUS(stat);
    ATTATCHSTATUSPTR(stat);

    /* Check for valid input arguments. */
    ASSERT(value, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL);
    ASSERT(string, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL);

    /* Parse string.  Return if nothing was parsed. */
    TRY(LALStringToD(stat->statusPtr, &re, string, &end), stat);
    TRY(LALStringToD(stat->statusPtr, &im, end, &end), stat);
    if (string == end) {
        if (endptr)
            *endptr = end;
        DETATCHSTATUSPTR(stat);
        RETURN(stat);
    }

    /* Set values and return. */
    *value = re;
    *value += im * I;
    if (endptr)
        *endptr = end;
    DETATCHSTATUSPTR(stat);
    RETURN(stat);
}



void
LALStringToGPS(LALStatus * stat, LIGOTimeGPS * value, const CHAR * string,
               CHAR ** endptr)
{

    /* Trick borrowed from LALStringToU8AndSign() to avoid compiler warnings */
    union {
        char *s;
        const char *cs;
    } bad;      /* there is a REASON for warnings... */

    const CHAR *here = string;  /* current position in string */
    INT4 signval;       /* sign of value (+1 or -1) */
    CHAR mantissa[32];  /* local string to store mantissa digits */
    INT4 mdigits;       /* number of digits in mantissa */
    INT4 dppos; /* position of decimal point in mantissa, i.e. the
                   number of mantissa digits preceding the decimal point
                   (initially -1 if no decimal point in input.) */
    const CHAR *ehere;  /* string pointer for parsing exponent */
    INT4 exponent;      /* exponent given in string */
    INT2 esignval;      /* sign of exponent value (+1 or -1) */
    UINT8 absValue;     /* magnitude of parsed number */
    CHAR *eend; /* substring following parsed number */
    INT4 nanosecSet;    /* flag to indicate if nanoseconds field has been set */
    INT4 idigit;        /* Digit value: 0 for 1s digit, 1 for 10s digit, etc. */

    INITSTATUS(stat);
    ATTATCHSTATUSPTR(stat);

    /* Check for valid input arguments. */
    ASSERT(value, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL);
    ASSERT(string, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL);

    /* Skip leading space, and read sign character, if any. */
    signval = 1;
    while (isspace(*here))
        here++;
    if (*here == '+')
        here++;
    else if (*here == '-') {
        signval = -1;
        here++;
    }

    /* Copy the mantissa into a local string, keeping track of the
       location of the decimal point */
    mdigits = 0;
    dppos = -1;
    while ((*here >= '0' && *here <= '9') || (*here == '.')) {

        if (*here == '.') {
            /* This is a decimal point */
            if (dppos >= 0) {
                /* This is second decimal point encountered, so parsing must stop */
                break;
            } else {
                /* Record the position of the decimal point */
                dppos = mdigits;
            }

        } else if (*here == '0' && mdigits == 1 && mantissa[0] == '0') {
            /* We only want to keep at most one leading zero.  This is an
               additional leading zero, so simply ignore it. */

        } else {
            /* This is a digit.  Append it to the local mantissa string, unless
               the mantissa string is already full, in which case ignore it */
            if ((UINT4) mdigits < sizeof(mantissa)) {
                mantissa[mdigits] = *here;
                mdigits++;
            }
        }

        here++;

    }

    /* If there is no mantissa, then return without consuming any characters
       and without modifying 'value' */
    if (mdigits == 0) {
        if (endptr) {
            bad.cs = string;    /* ... and this avoids the warnings... BAD! */
            *endptr = bad.s;
        }
        DETATCHSTATUSPTR(stat);
        RETURN(stat);
    }

    /* If there was no explicit decimal point, then it is implicitly
       after all of the mantissa digits */
    if (dppos == -1) {
        dppos = mdigits;
    }

    /* Read the exponent, if present */
    exponent = 0;
    if (*here == 'E' || *here == 'e') {
        /* So far, this looks like an exponent.  Set working pointer. */
        ehere = here + 1;

        /* Parse the exponent value */
        absValue = LALStringToU8AndSign(&esignval, ehere, &eend);
        if (eend == ehere) {
            /* Nothing was parsed, so this isn't a valid exponent.  Leave
               things as they are, with 'here' pointing to the 'E' or 'e'
               that we thought introduced an exponent. */
        } else {
            /* We successfully parsed the exponent */
            exponent = (INT4) (esignval * absValue);
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

    if (dppos > 11) {
        /* This is an overflow (positive) or underflow (negative) */
        if (signval == 1) {
            value->gpsSeconds = (INT4) LAL_INT4_MAX;
            value->gpsNanoSeconds = 999999999;
        } else {
            value->gpsSeconds = (INT4) (-LAL_INT4_ABSMIN);
            value->gpsNanoSeconds = 0;
        }

    } else if (dppos < -9) {
        /* The time is effectively zero */
        value->gpsSeconds = 0;
        value->gpsNanoSeconds = 0;

    } else {

        /* Pick out the integer part */
        absValue = 0;
        for (idigit = 0; idigit < dppos && idigit < mdigits; idigit++) {
            absValue = 10 * absValue + (UINT8) (mantissa[idigit] - '0');
        }
        /* Fill in missing powers of ten if not all digits were present */
        for (; idigit < dppos; idigit++) {
            absValue *= 10;
        }

        /* Cap (if necessary) and cast */
        if (signval > 0) {
            if (absValue > LAL_INT4_MAX) {
                value->gpsSeconds = (INT4) (LAL_INT4_MAX);
                value->gpsNanoSeconds = 999999999;
                nanosecSet = 1;
            } else
                value->gpsSeconds = (INT4) (absValue);
        } else {
            if (absValue >= LAL_INT4_ABSMIN) {
                value->gpsSeconds = (INT4) (-LAL_INT4_ABSMIN);
                value->gpsNanoSeconds = 0;
                nanosecSet = 1;
            } else
                value->gpsSeconds = (INT4) (-absValue);
        }

        /* Finally, set nanoseconds field (if not already set by over/underflow) */
        if (!nanosecSet) {

            absValue = 0;
            for (idigit = dppos; idigit < dppos + 9 && idigit < mdigits;
                 idigit++) {
                if (idigit >= 0) {
                    absValue =
                        10 * absValue + (UINT8) (mantissa[idigit] - '0');
                }
            }
            /* If there is another digit, use it to round */
            if (idigit == dppos + 9 && idigit < mdigits) {
                if (mantissa[idigit] >= '5') {
                    absValue++;
                }
            }
            /* Fill in missing powers of ten if not all digits were present */
            for (; idigit < dppos + 9; idigit++) {
                absValue *= 10;
            }

            value->gpsNanoSeconds = (INT4) (signval * absValue);

            /* Check for wraparound due to rounding */
            if (value->gpsNanoSeconds >= 1000000000) {
                value->gpsNanoSeconds -= 1000000000;
                value->gpsSeconds += 1;
            }
            /* Ensure that nanoseconds field is nonnegative */
            while (value->gpsNanoSeconds < 0) {
                value->gpsNanoSeconds += 1000000000;
                value->gpsSeconds -= 1;
            }

        }

    }

    /* Set end pointer (if passed) and return. */
    if (endptr) {
        bad.cs = here;  /* ... and this avoids the warnings... BAD! */
        *endptr = bad.s;
    }
    DETATCHSTATUSPTR(stat);
    RETURN(stat);
}
