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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

/**
 * \file
 * \ingroup StringInput_h
 * \authors Creighton, T. D.
 * \authors Shawhan, P. S.
 *
 * \brief Converts a string into a numerical value.
 *
 * ### Description ###
 *
 * These routines parse the string <tt>*string</tt> and compute a numerical
 * value <tt>*value</tt> of the appropriate datatype.  If
 * \c endptr\f$\neq\f$\c NULL, then after conversion <tt>*endptr</tt>
 * will point to the character after the last character used in the
 * conversion.  The routine will always return without error if the
 * arguments are valid, regardless of the contents of \c string;
 * failure to parse a number is indicated by <tt>*endptr</tt> being set
 * equal to \c string (provided \c endptr\f$\neq\f$\c NULL).
 *
 * For integer or floating-point conversion to occur, \c string must
 * consist of zero or more whitespace characters followed by a number in
 * any standard base-ten integer or floating-point representation
 * (described in detail below); the conversion will stop as soon as the
 * routine encounters a character that is not part of the number.  For
 * instance, parsing the string <tt>"   123bad"</tt> will return
 * <tt>*value</tt>=123, and <tt>*endptr</tt> will point to the substring
 * <tt>"bad"</tt>; it is up to the calling routine to determine whether
 * this is an acceptable input.  By contrast, parsing the string
 * <tt>" bad"</tt> will leave <tt>*value</tt> unchanged and <tt>*endptr</tt>
 * pointing to the start of the original string.  In general, if the
 * routine returns with <tt>*endptr</tt>\f$\neq\f$\c string and
 * <tt>**endptr</tt> is a whitespace or <tt>'\0'</tt> character, then the
 * format was unambiguously acceptable.
 *
 * Complex conversion is essentially equivalent to performing
 * floating-point conversion on \c string to get the real part, and
 * again on <tt>*endptr</tt> to get the imaginary part.  Normally this
 * means that an acceptable format is two floating-point representations
 * separated by (and possibly preceded with) whitespace, although the
 * intervening whitespace can be omitted if the separation between the
 * two numbers is unambiguous.  Thus the string <tt>"-1.0+3.5"</tt> will be
 * unambiguously read as \f$-1.0+3.5i\f$, but <tt>"-1.03.5"</tt> will be read
 * as \f$-1.03+0.5i\f$ (since the conversion of the real part stops at the
 * second <tt>'.'</tt> character), which may or may not be the intended
 * conversion.
 *
 * ### Algorithm ###
 *
 * These functions emulate the standard C functions <tt>strtol()</tt>,
 * <tt>strtoul()</tt>, and <tt>strtod()</tt>, except that they follow LAL
 * calling conventions and return values of the appropriate LAL datatypes.
 * For integer conversion, only base-ten (decimal) representations are
 * supported.  Otherwise, the valid format is as for the corresponding C
 * functions, which we summarize below:
 *
 * A string to be converted to an \c INT\f$n\f$ (where \f$n\f$=2, 4, or 8)
 * consists of zero or more whitespace characters as determined by
 * <tt>isspace()</tt>, followed optionally by a single <tt>'+'</tt> or
 * <tt>'-'</tt> character, followed by one or more decimal digits;
 * conversion stops at the first non-digit character after this.  If the
 * result would overflow or underflow the \c INT\f$n\f$ representation,
 * then the value is set to \c LAL_INT\f$n\f$\c _MAX\f$=2^{8n-1}-1\f$ or \c
 * LAL_INT\f$n\f$\c _MIN\f$=-2^{8n-1}\f$, respectively.
 *
 * A string to be converted to a \c UINT\f$n\f$ follows the same format,
 * except that a leading negative sign character <tt>'-'</tt> is \e ignored
 * (the routine will compute the magnitude of the result), and the return
 * value is capped at \c LAL_UINT\f$n\f$\c _MAX\f$=2^{8n}-1\f$.
 *
 * A string to be converted to a floating-point number (\c REAL4 or \c
 * REAL8) consists of zero or more whitespace characters as determined by
 * <tt>isspace()</tt>, followed optionally by a single <tt>'+'</tt> or
 * <tt>'-'</tt> character, followed by a sequence of one or more decimal
 * digits optionally containing a decimal point <tt>'.'</tt>, optionally
 * followed by an exponent.  An exponent consists of a single <tt>'E'</tt>
 * or <tt>'e'</tt> character, followed optionally by a single <tt>'+'</tt>
 * or <tt>'-'</tt> character, followed by a sequence of one or more decimal
 * digits.  If the converted value would overflow, \f$\pm\f$\c
 * LAL_REAL\f$n\f$\c _MAX is returned, as appropriate.  If the value would
 * underflow, 0 is returned.
 *
 * A string to be converted to a complex number (\c COMPLEX8 or \c
 * COMPLEX16) consists of two floating-point format substrings concatenated
 * together, where the first character of the second substring cannot be
 * interpreted as a continuation of the first number.  Usually this means
 * that the second substring will contain at least one leading whitespace
 * character, though this is not strictly necessary.  Overflow or underflow
 * is dealt with as above.
 */

/** \cond DONT_DOXYGEN */
#include <complex.h>
#include <ctype.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/StringInput.h>

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
