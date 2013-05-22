/*
 * Copyright (C) 2007 David Chin, Jolien Creighton
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with with program; see the file COPYING. If not, write to the Free
 * Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/XLALError.h>


#define SUCCESS              0
#define FAIL_JULIAN_DAY      1
#define FAIL_MOD_JULIAN_DAY  2
#define FAIL_JULIAN_DATE     3
#define FAIL_MOD_JULIAN_DATE 4

/* good to 1 micro-day  */
const REAL8 julian_precision = 1.e-6;

const REAL8 coarse_precision = 0.001;

/*  Output of NASA/AMES code:
 Date/time: 1 1 0 12 0 0
 PDS format: "2000-01-01T12:00:00.000Z"
 SQL format: "Jan 1, 2000 12:00:00:000"

  UTC day  =        0
  UTC secs =    43200.000000
  TAI secs =          32.000000
   ET secs =          64.183927

  JD (UTC) =  2451545.000000
  JD (TAI) =  2451545.000370
  JD (ET)  =  2451545.001115
 MJD (UTC) =    51544.500000
 MJD (TAI) =    51544.500370
 MJD (ET)  =    51544.501115

  Date/time: 94 11 16 0 0 0
 PDS format: "1994-11-16T00:00:00.000Z"
 SQL format: "Nov 16, 1994 00:00:00:000"

  UTC day  =    -1872
  UTC secs =        0.000000
  TAI secs =  -161783971.000000
   ET secs =  -161783938.817245

  JD (UTC) =  2449672.500000
  JD (TAI) =  2449672.500336
  JD (ET)  =  2449672.501081
 MJD (UTC) =    49672.000000
 MJD (TAI) =    49672.000336
 MJD (ET)  =    49672.001081

 Date/time: 01 05 15 02 37 54
 PDS format: "2001-05-15T02:37:54.000Z"
 SQL format: "May 15, 2001 02:37:54:000"

  UTC day  =      500
  UTC secs =     9474.000000
  TAI secs =    43166306.000000
   ET secs =    43166338.185257

  JD (UTC) =  2452044.609653
  JD (TAI) =  2452044.610023
  JD (ET)  =  2452044.610768
 MJD (UTC) =    52044.109653
 MJD (TAI) =    52044.110023
 MJD (ET)  =    52044.110768


 http://tycho.usno.navy.mil/cgi-bin/date
15-May-01
MJD 52044.109653
UTC 02:37:54

*/


/*
 * pass 0 for expected_julian_day or expected_modified_julian_day to disable testing
 */

static int test(const struct tm *utc, double expected_julian_day, double expected_modified_julian_day, int line)
{
	double julian_day;
	double modified_julian_day;
	int result = 0;

	if(lalDebugLevel)
		fprintf(stderr, "Testing %s ...\n", asctime(utc));

	julian_day = XLALJulianDay(utc);
	modified_julian_day = XLALModifiedJulianDay(utc);

	if(expected_julian_day && (julian_day != expected_julian_day)) {
		fprintf(stderr, "XLALJulianDay() failed (line %d):  expected %.17g got %.17g\n", line, expected_julian_day, julian_day);
		result = -1;
	} else if(lalDebugLevel) {
		fprintf(stderr, "XLALJulianDay() returned %.16g\n", julian_day);
	}
	if(expected_modified_julian_day && (modified_julian_day != expected_modified_julian_day)) {
		fprintf(stderr, "XLALModifiedJulianDay() failed (line %d):  expected %.17g got %.17g\n", line, expected_modified_julian_day, modified_julian_day);
		result = -1;
	} else if(lalDebugLevel) {
		fprintf(stderr, "XLALModifiedJulianDay() returned %.17g\n", modified_julian_day);
	}

	return result;
}


int main(void)
{
	time_t now;
	struct tm utc;
#if 0
	REAL8 ref_julian_day;
	REAL8 julian_day;
	REAL8 ref_mod_julian_day;
	REAL8 mod_julian_day;
#endif


	/*
	 * Distinctly not robust
	 */


	/*
	 * Get current local time
	 */

	time(&now);
	if(test(localtime(&now), 0, 0, __LINE__))
		return 1;

	/*
	 * Check Julian Day/Date for special dates and times
	 */

	utc.tm_sec = 0;
	utc.tm_min = 0;
	utc.tm_hour = 12;
	utc.tm_mday = 1;
	utc.tm_mon = 0;
	utc.tm_year = 100;
	utc.tm_wday = 6;
	utc.tm_yday = 0;
	utc.tm_isdst = 0;

	if(test(&utc, 2451545.0, 51544.0, __LINE__))
		return 1;

	utc.tm_sec = 0;
	utc.tm_min = 0;
	utc.tm_hour = 11;
	utc.tm_mday = 1;
	utc.tm_mon = 0;
	utc.tm_year = 100;
	utc.tm_wday = 6;
	utc.tm_yday = 0;
	utc.tm_isdst = 0;

	if(test(&utc, 2451544.9583333333, 51544.0, __LINE__))
		return 1;

	/* */

	if(lalDebugLevel > 1) {
		utc.tm_sec = 0;
		utc.tm_min = 0;
		utc.tm_hour = 11;
		utc.tm_mday = 1;
		utc.tm_mon = 0;
		utc.tm_year = -100;
		utc.tm_wday = 6;
		utc.tm_yday = 0;
		utc.tm_isdst = 0;

		/* here, we EXPECT an error */
		if(test(&utc, XLAL_REAL8_FAIL_NAN, XLAL_REAL8_FAIL_NAN, __LINE__))
			return 1;

	}



	/* */

	utc.tm_sec = 0;
	utc.tm_min = 0;
	utc.tm_hour = 0;
	utc.tm_mday = 16;
	utc.tm_mon = 10;
	utc.tm_year = 94;
	utc.tm_wday = 3;
	utc.tm_yday = 0;
	utc.tm_isdst = 0;

	if(test(&utc, 2449672.5, 49672.0, __LINE__))
		return 1;

	/* */

	utc.tm_sec = 54;
	utc.tm_min = 37;
	utc.tm_hour = 2;
	utc.tm_mday = 15;
	utc.tm_mon = 4;
	utc.tm_year = 101;
	utc.tm_wday = 1;
	utc.tm_yday = 0;
	utc.tm_isdst = 1;

	if(test(&utc, 2452044.6096527777, 52044.0, __LINE__))
		return 1;

	return SUCCESS;
}
