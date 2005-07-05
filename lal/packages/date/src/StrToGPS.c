#include <errno.h>
#include <stdlib.h>
#include <lal/Date.h>
#include <lal/LALDatatypes.h>

int XLALStrToGPS(LIGOTimeGPS *time, const char *nptr, char **endptr)
{
	const char *func = "XLALStrToGPS";
	char *p, *q;
	int i;

	/* clear C library errno so we can check for failures */
	errno = 0;

	/* parse the integer part */
	time->gpsSeconds = strtol(nptr, &p, 10);

	/* check for a fractional part */
	if(*p == '.') {
		/* parse the fractional part */
		p++;
		time->gpsNanoSeconds = strtol(p, &q, 10);
		for(i = q - p - 9; i > 0; time->gpsNanoSeconds /= 10, i--);
		for(; i < 0; time->gpsNanoSeconds *= 10, i++);
	} else {
		q = p;
		time->gpsNanoSeconds = 0;
	}

	if(endptr)
		*endptr = q;

	if(errno == ERANGE)
		XLAL_ERROR(func, XLAL_ERANGE);
	
	return(0);
}
