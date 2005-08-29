#include <ctype.h>
#include <errno.h>
#include <stdlib.h>
#include <lal/Date.h>
#include <lal/LALDatatypes.h>

/*
 * Check if a string is a number in scientific notation.
 */

static const char *skipsign(const char *s)
{
	if(*s == '+' || *s == '-')
		s++;
	return(s);
}

static int isscientific(const char *s)
{
	int radix_count = 0;
	const char *p;

	/* skip leading white space, and an optional sign */
	while(isspace(*s))
		s++;
	s = skipsign(s);

	/* check for a hex prefix */
	if(*s == '0' && (*(s + 1) == 'X' || *(s + 1) == 'x')) {
		/* skip hex prefix */
		s += 2;

		/* mark current location */
		p = s;

		/* allow hex digits and the radix character.  FIXME: radix
		 * character should be obtained from current locale .*/
		while(isxdigit(*s) || *s == '.') {
			radix_count += *s == '.';
			s++;
		}

		/* if no digits were found, or more than one radix was found,
		 * or the next character is not the exponent prefix the answer
		 * is no */
		if(p == s || radix_count > 1 || (*s != 'P' && *s != 'p'))
			return(0);

		/* skip exponent prefix, and an optional sign */
		s = skipsign(++s);

		/* if the next character is not a hex digit then the answer is
		 * no */
		return(isxdigit(*s));
	} else {
		/* mark current location */
		p = s;

		/* allow decimal digits and the radix character.  FIXME: radix
		 * character should be obtained from current locale. */
		while(isdigit(*s) || *s == '.') {
			radix_count += *s == '.';
			s++;
		}

		/* if no digits were found, or more than one radix was found,
		 * or the next character is not the exponent prefix the answer
		 * is no */
		if(p == s || radix_count > 1 || (*s != 'E' && *s != 'e'))
			return(0);

		/* skip exponent prefix, and an optional sign */
		s = skipsign(++s);

		/* if the next character is not a decimal digit then the answer
		 * is no */
		return(isdigit(*s));
	}
}


/*
 * Parse an ASCII string into a LIGOTimeGPS structure.
 */

int XLALStrToGPS(LIGOTimeGPS *time, const char *nptr, char **endptr)
{
	const char *func = "XLALStrToGPS";
	char *p, *q;
	int olderrno, i;

	/* save and clear C library errno so we can check for failures */
	olderrno = errno;
	errno = 0;

	/* parse as a double if in scientific notation */
	if(isscientific(nptr)) {
		XLALFloatToGPS(time, strtod(nptr, &q));
	} else {
		/* parse the integer part */
		time->gpsSeconds = strtol(nptr, &p, 10);
		time->gpsNanoSeconds = 0;

		/* check for and parse a fractional part */
		if(*p == '.') {
			if(isdigit(*++p)) {
				time->gpsNanoSeconds = strtol(p, &q, 10);
				/* move the decimal place to the correct location */
				for(i = q - p - 9; i > 0; i--)
					time->gpsNanoSeconds /= 10;
				for(; i < 0; i++)
					time->gpsNanoSeconds *= 10;
				/* non-integer negative times require a correction */
				if((time->gpsSeconds < 0) && time->gpsNanoSeconds) {
					time->gpsSeconds--;
					time->gpsNanoSeconds = 1000000000 - time->gpsNanoSeconds;
				}
			} else
				q = p;
		} else
			q = p;
	}

	/* save pointer to first character after converted characters */
	if(endptr)
		*endptr = q;

	/* check for failures and restore errno if there weren't any */
	if(errno == ERANGE)
		XLAL_ERROR(func, XLAL_ERANGE);
	errno = olderrno;

	/* success */
	return(0);
}
