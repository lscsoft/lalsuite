/*
 * Copyright (C) 2007  Kipp Cannon
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */


#include <ctype.h>
#include <errno.h>
#include <locale.h>
#include <stdlib.h>
#include <lal/Date.h>
#include <lal/LALDatatypes.h>


/*
 * Check for a base 10 or base 16 number.
 */


static int isbase10(const char *s, int radix)
{
	if(*s == radix)
		s++;
	if(isdigit(*s))
		return(1);
	return(0);
}


static int isbase16(const char *s, int radix)
{
	if(*s == '0') {
		s++;
		if(*s == 'X' || *s == 'x') {
			s++;
			if(*s == radix)
				s++;
			if(isxdigit(*s))
				return(1);
		}
	}
	return(0);
}


/*
 * Check that a string contains an exponent.
 */


static int isdecimalexp(const char *s)
{
	if(*s == 'E' || *s == 'e') {
		s++;
		if(*s == '+' || *s == '-')
			s++;
		if(isdigit(*s))
			return(1);
	}
	return(0);
}


static int isbinaryexp(const char *s)
{
	if(*s == 'P' || *s == 'p') {
		s++;
		if(*s == '+' || *s == '-')
			s++;
		if(isdigit(*s))
			return(1);
	}
	return(0);
}


/**
 * Parse an ASCII string into a LIGOTimeGPS structure.
 */
int XLALStrToGPS(LIGOTimeGPS *t, const char *nptr, char **endptr)
{
	union { char *s; const char *cs; } pconv; /* this is bad */
	int olderrno;
	int radix;
	char *digits;
	int len=0;
	int sign;
	int base;
	int radixpos;
	int exppart;

	/* save and clear C library errno so we can check for failures */
	olderrno = errno;
	errno = 0;

	/* retrieve the radix character */
	radix = localeconv()->decimal_point[0];

	/* this is bad ... there is a reason for warnings! */
	pconv.cs = nptr;

	/* consume leading white space */
	while(isspace(*(pconv.cs)))
		(pconv.cs)++;
	if(endptr)
		*endptr  = pconv.s;

	/* determine the sign */
	if(*(pconv.cs) == '-') {
		sign = -1;
		(pconv.cs)++;
	} else if(*(pconv.cs) == '+') {
		sign = +1;
		(pconv.cs)++;
	} else
		sign = +1;

	/* determine the base */
	if(isbase16((pconv.cs), radix)) {
		base = 16;
		(pconv.cs) += 2;
	} else if(isbase10((pconv.cs), radix)) {
		base = 10;
	} else {
		/* this isn't a recognized number */
		XLALGPSSet(t, 0, 0);
		return(0);
	}

	/* count the number of digits including the radix but not including
	 * the exponent. */
	radixpos = -1;
	switch(base) {
	case 10:
		for(len = 0; 1; len++) {
			if(isdigit((pconv.cs)[len]))
				continue;
			if((pconv.cs)[len] == radix && radixpos < 0) {
				radixpos = len;
				continue;
			}
			break;
		}
		break;

	case 16:
		for(len = 0; 1; len++) {
			if(isxdigit((pconv.cs)[len]))
				continue;
			if((pconv.cs)[len] == radix && radixpos < 0) {
				radixpos = len;
				continue;
			}
			break;
		}
		break;
	}

	/* copy the digits into a scratch space, removing the radix character
	 * if one was found */
	if(radixpos >= 0) {
		digits = malloc(len + 1);
		memcpy(digits, (pconv.cs), radixpos);
		memcpy(digits + radixpos, (pconv.cs) + radixpos + 1, len - radixpos - 1);
		digits[len - 1] = '\0';
		(pconv.cs) += len;
		len--;
	} else {
		digits = malloc(len + 2);
		memcpy(digits, (pconv.cs), len);
		digits[len] = '\0';
		radixpos = len;
		(pconv.cs) += len;
	}

	/* check for and parse an exponent, performing an adjustment of the
	 * radix position */
	exppart = 1;
	switch(base) {
	case 10:
		/* exponent is the number of powers of 10 */
		if(isdecimalexp((pconv.cs)))
			radixpos += strtol((pconv.cs) + 1, &pconv.s, 10);
		break;

	case 16:
		/* exponent is the number of powers of 2 */
		if(isbinaryexp((pconv.cs))) {
			exppart = strtol((pconv.cs) + 1, &pconv.s, 10);
			radixpos += exppart / 4;
			exppart %= 4;
			if(exppart < 0) {
				radixpos--;
				exppart += 4;
			}
			exppart = 1 << exppart;
		}
		break;
	}

	/* save end of converted characters */
	if(endptr)
		*endptr = pconv.s;

	/* insert the radix character, padding the scratch digits with zeroes
	 * if needed */
	if(radixpos < 2) {
		digits = realloc(digits, len + 2 + (2 - radixpos));
		memmove(digits + (2 - radixpos) + 1, digits, len + 1);
		memset(digits, '0', (2 - radixpos) + 1);
		if(radixpos == 1)
			digits[1] = digits[2];
		radixpos = 2;
	} else if(radixpos > len) {
		digits = realloc(digits, radixpos + 2);
		memset(digits + len, '0', radixpos - len);
		digits[radixpos + 1] = '\0';
	} else {
		memmove(digits + radixpos + 1, digits + radixpos, len - radixpos + 1);
	}
	digits[radixpos] = radix;

	/* parse the integer part */
	XLALINT8NSToGPS(t, sign * strtol(digits, NULL, base) * exppart * XLAL_BILLION_INT8);

	/* parse the fractional part */
	if(errno != ERANGE) {
		switch(base) {
		case 10:
			break;

		case 16:
			digits[radixpos - 2] = '0';
			digits[radixpos - 1] = 'x';
			radixpos -= 2;
			break;
		}
		XLALGPSAdd(t, sign * strtod(digits + radixpos, NULL) * exppart);
	}

	/* free the scratch space */
	free(digits);

	/* check for failures and restore errno if there weren't any */
	if(errno == ERANGE)
		XLAL_ERROR(XLAL_ERANGE);
	errno = olderrno;

	/* success */
	return(0);
}


/**
 * \ingroup Date_h
 * \brief Return a string containing the ASCII base 10 representation of a
 * LIGOTimeGPS.  If s is not NULL, then the string is written to that
 * location which must be large enough to hold the string plus a 0
 * terminator.  If s is NULL then a new buffer is allocated, and the string
 * written to it.  The return value is the address of the string or NULL on
 * failure.
 */
char *XLALGPSToStr(char *s, const LIGOTimeGPS *t)
{
	/* so we can play with it */
	LIGOTimeGPS copy = *t;

	/* make sure we've got a buffer */

	if(!s) {
		/* 21 = 9 digits to the right of the decimal point +
		 * decimal point + upto 10 digits to the left of the
		 * decimal point plus an optional sign + a null */
		s = XLALMalloc(21 * sizeof(*s));
		if(!s)
			XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	/* normalize the fractional part */

	while(labs(copy.gpsNanoSeconds) > XLAL_BILLION_INT4) {
		if(copy.gpsNanoSeconds < 0) {
			copy.gpsSeconds -= 1;
			copy.gpsNanoSeconds += XLAL_BILLION_INT4;
		} else {
			copy.gpsSeconds += 1;
			copy.gpsNanoSeconds -= XLAL_BILLION_INT4;
		}
	}

	/* if both components are non-zero, make sure they have the same
	 * sign */

	if(copy.gpsSeconds > 0 && copy.gpsNanoSeconds < 0) {
		copy.gpsSeconds -= 1;
		copy.gpsNanoSeconds += XLAL_BILLION_INT4;
	} else if(copy.gpsSeconds < 0 && copy.gpsNanoSeconds > 0) {
		copy.gpsSeconds += 1;
		copy.gpsNanoSeconds -= XLAL_BILLION_INT4;
	}

	/* print */

	if(copy.gpsSeconds < 0 || copy.gpsNanoSeconds < 0)
		/* number is negative */
		sprintf(s, "-%ld.%09ld", labs(copy.gpsSeconds), labs(copy.gpsNanoSeconds));
	else
		/* number is non-negative */
		sprintf(s, "%ld.%09ld", (long) copy.gpsSeconds, (long) copy.gpsNanoSeconds);

	return s;
}
