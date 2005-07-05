#include <ctype.h>
#include <errno.h>
#include <stdlib.h>
#include <lal/Date.h>
#include <lal/LALDatatypes.h>

int XLALStrToGPS(LIGOTimeGPS *time, const char *nptr, char **endptr)
{
	const char *func = "XLALStrToGPS";
	char *p, *q;
	int olderrno, i;

	/* save and clear C library errno so we can check for failures */
	olderrno = errno;
	errno = 0;

	/* parse the integer part */
	time->gpsSeconds = strtol(nptr, &p, 10);
	time->gpsNanoSeconds = 0;

	/* check for and parse a fractional part */
	if(*p == '.') {
		if(isdigit(*++p)) {
			time->gpsNanoSeconds = strtol(p, &q, 10);
			for(i = q - p - 9; i > 0; i--)
				time->gpsNanoSeconds /= 10;
			for(; i < 0; i++)
				time->gpsNanoSeconds *= 10;
		} else
			q = p;
	} else
		q = p;

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
