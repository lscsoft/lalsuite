/*
*  Copyright (C) 2007 David Chin, Jolien Creighton
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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include <lal/LALStdlib.h>
#include <lal/Date.h>

INT4 lalDebugLevel = 0;

NRCSID(LALTESTLEAPSECSC, "$Id$");

static int do_test(LIGOTimeGPS gps, int tai_utc_before, int tai_utc_after)
{
	int result = 0;
	int lal_tai_utc;

	if(lalDebugLevel > 2)
		printf("TestLeapSecs: BEFORE LEAP SECOND ADDED\n");
	lal_tai_utc = XLALLeapSeconds(gps.gpsSeconds - 1);
	if(lalDebugLevel > 0)
		printf("\tGPS = %9d;    TAI-UTC = %d\n", gps.gpsSeconds - 1, lal_tai_utc);

	if(XLALGetBaseErrno() && lalDebugLevel > 0) {
		XLAL_PERROR("do_test()");
		result = -1;
	}
	if(lal_tai_utc != tai_utc_before) {
		if(lalDebugLevel > 0)
			fprintf(stderr, "TestLeapSecs: XLALLeapSeconds() returned wrong value: expected %d, got %d\n", tai_utc_before, lal_tai_utc);
		result = -1;
	}

	if(lalDebugLevel > 2)
		printf("TestLeapSecs: AFTER LEAP SECOND ADDED\n");
	lal_tai_utc = XLALLeapSeconds(gps.gpsSeconds);
	if(lalDebugLevel > 0)
		printf("\tGPS = %9d;    TAI-UTC = %d\n\n", gps.gpsSeconds, lal_tai_utc);

	if(XLALGetBaseErrno() && lalDebugLevel > 0) {
		XLAL_PERROR("do_test()");
		result = -1;
	}
	if(lal_tai_utc != tai_utc_after) {
		if(lalDebugLevel > 0)
			fprintf(stderr, "TestLeapSecs: XLALLeapSeconds() returned wrong value: expected %d, got %d\n", tai_utc_after, lal_tai_utc);
		result = -1;
	}
	return result;
}


int main(int argc, char *argv[])
{
	struct {
		LIGOTimeGPS gps;	/* GPS time when leap sec was introduced */
		int tai_utc;	/* GPS-UTC at GPS time gps */
	} leapsec_data[] = {
		{{        0, 0}, 19},	/* 1980-Jan-06 */
		{{ 46828800, 0}, 20},	/* 1981-Jul-01 */
		{{ 78364801, 0}, 21},	/* 1982-Jul-01 */
		{{109900802, 0}, 22},	/* 1983-Jul-01 */
		{{173059203, 0}, 23},	/* 1985-Jul-01 */
		{{252028804, 0}, 24},	/* 1988-Jan-01 */
		{{315187205, 0}, 25},	/* 1990-Jan-01 */
		{{346723206, 0}, 26},	/* 1991-Jan-01 */
		{{393984007, 0}, 27},	/* 1992-Jul-01 */
		{{425520008, 0}, 28},	/* 1993-Jul-01 */
		{{457056009, 0}, 29},	/* 1994-Jul-01 */
		{{504489610, 0}, 30},	/* 1996-Jan-01 */
		{{551750411, 0}, 31},	/* 1997-Jul-01 */
		{{599184012, 0}, 32},	/* 1999-Jan-01 */
	};
	unsigned i;

	if(argc > 1)
		lalDebugLevel = atoi(argv[1]);

	for(i = 1; i < sizeof(leapsec_data) / sizeof(*leapsec_data); i++)
		do_test(leapsec_data[i].gps, leapsec_data[i-1].tai_utc, leapsec_data[i].tai_utc);

	return 0;
}
