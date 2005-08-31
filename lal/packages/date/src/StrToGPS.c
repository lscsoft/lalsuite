#include <ctype.h>
#include <errno.h>
#include <stdlib.h>
#include <lal/Date.h>
#include <lal/LALDatatypes.h>


/*
 * Return true if x is an octal numeral.
 */

static int isodigit(int x)
{
	return(x >= '0' && x <= '7');
}


/*
 * Check for a base 8 or base 16 number.
 */

static int isbase8(const char *s)
{
	if(*s == '0') {
		while(*s == '0')
			s++;
		if(isdigit(*s))
			return(1);
	}
	return(0);
}

static int isbase16(const char *s)
{
	if(*s == '0') {
		s++;
		if(*s == 'X' || *s == 'x') {
			s++;
			if(isxdigit(*s) || *s == '.')
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


/*
 * Parse an ASCII string into a LIGOTimeGPS structure.
 */

int XLALStrToGPS(LIGOTimeGPS *time, const char *nptr, char **endptr)
{
	const char *func = "XLALStrToGPS";
	int olderrno;
	char *digits;
	int len;
	int sign;
	int base;
	int radixpos;
	int exppart;

	/* save and clear C library errno so we can check for failures */
	olderrno = errno;
	errno = 0;

	/* consume leading white space */
	while(isspace(*nptr))
		nptr++;
	if(endptr)
		*endptr = nptr;

	/* determine the sign */
	if(*nptr == '-') {
		sign = -1;
		nptr++;
	} else if(*nptr == '+') {
		sign = +1;
		nptr++;
	} else
		sign = +1;

	/* determine the base */
	if(isbase8(nptr)) {
		base = 8;
		nptr++;
	} else if(isbase16(nptr)) {
		base = 16;
		nptr += 2;
	} else
		base = 10;
	/*if(*nptr != '0')
		base = 10;
	else {
		*nptr++;
		if(*nptr == 'X' || *nptr == 'x') {
			base = 16;
			nptr++;
		} else
			base = 8;
	}*/

	/* count the number of digits including the radix but not including
	 * the exponent.  FIXME: radix character should be obtained from
	 * the locale */
	radixpos = -1;
	switch(base) {
	case 8:
		for(len = 0; 1; len++) {
			if(isodigit(nptr[len]))
				continue;
			if(nptr[len] == '.' && radixpos < 0) {
				radixpos = len;
				continue;
			}
			break;
		}
		break;
	
	case 10:
		for(len = 0; 1; len++) {
			if(isdigit(nptr[len]))
				continue;
			if(nptr[len] == '.' && radixpos < 0) {
				radixpos = len;
				continue;
			}
			break;
		}
		break;
	
	case 16:
		for(len = 0; 1; len++) {
			if(isxdigit(nptr[len]))
				continue;
			if(nptr[len] == '.' && radixpos < 0) {
				radixpos = len;
				continue;
			}
			break;
		}
		break;
	}

	/* check that we have a number */
	if(!len) {
		XLALGPSSet(time, 0, 0);
		return(0);
	}

	/* copy the digits into a scratch space, removing the radix character
	 * if one was found */
	if(radixpos >= 0) {
		digits = malloc(len + 1);
		memcpy(digits, nptr, radixpos);
		memcpy(digits + radixpos, nptr + radixpos + 1, len - radixpos - 1);
		digits[len - 1] = '\0';
		nptr += len;
		len--;
	} else {
		digits = malloc(len + 2);
		memcpy(digits, nptr, len);
		digits[len] = '\0';
		radixpos = len;
		nptr += len;
	}

	/* check for and parse an exponent, performing an adjustment of the
	 * radix position */
	switch(base) {
	case 8:
		/* C library has no format for octal exponents */
		exppart = 1;
		break;

	case 10:
		/* exponent is the number of powers of 10 */
		if(isdecimalexp(nptr))
			radixpos += strtol(++nptr, &nptr, 10);
		exppart = 1;
		break;

	case 16:
		/* exponent is the number of powers of 2 */
		if(isbinaryexp(nptr)) {
			exppart = strtol(++nptr, &nptr, 10);
			radixpos += exppart / 4;
			exppart = 1 << (exppart % 4);
		} else
			exppart = 1;
		break;
	}

	/* save pointer to first character after converted characters */
	if(endptr)
		*endptr = nptr;

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
	digits[radixpos] = '.';

	/* parse the integer part */
	XLALINT8NSToGPS(time, sign * strtol(digits, NULL, base) * exppart * 1000000000ll);

	/* parse the fractional part */
	if(errno != ERANGE) {
		switch(base) {
		case 8:
			break;

		case 10:
			break;

		case 16:
			digits[radixpos - 2] = '0';
			digits[radixpos - 1] = 'x';
			radixpos -= 2;
			break;
		}
		XLALGPSAdd(time, sign * strtod(digits + radixpos, NULL) * exppart);
	}

	/* free the scratch space */
	free(digits);

	/* check for failures and restore errno if there weren't any */
	if(errno == ERANGE)
		XLAL_ERROR(func, XLAL_ERANGE);
	errno = olderrno;

	/* success */
	return(0);
}
