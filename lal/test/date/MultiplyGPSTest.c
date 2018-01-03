/*
 * Copyright (C) 2015 Kipp Cannon
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
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <lal/LALStdlib.h>
#include <lal/Date.h>


/* caution:  not thread safe;  each call returns a pointer to one of 4
 * internal buffers that are used cyclically, so the 5th call overwrites
 * the contents of the first buffer;  be aware of this when using for
 * multiple arguments in a printf() call. */
static const char *print_LIGOTimeGPS(const LIGOTimeGPS *gps)
{
	static int i = 0;
	static char buffer[4][512];
	return XLALGPSToStr(buffer[i++ & 3], gps);
}


static LIGOTimeGPS compute_with_bc(const char *expr)
{
	LIGOTimeGPS gps;
	char cmd[512];
	char res[512];
	FILE *stream;

	sprintf(cmd, "echo 'scale=17 ; %s' | bc", expr);
	stream = popen(cmd, "r");
	fgets(res, 512, stream);
	pclose(stream);
	XLALStrToGPS(&gps, res, NULL);

	return gps;
}


static LIGOTimeGPS random_LIGOTimeGPS(void)
{
	LIGOTimeGPS gps;

	gps.gpsSeconds = floor((rand() / (double) RAND_MAX - 0.5) * 200000000);
	gps.gpsNanoSeconds = floor(rand() / ((double) RAND_MAX + 1) * 1000000000);

	return gps;
}


int main(void)
{
	double tolerance = 1e-9;
	int i;

	srand(time(NULL));

	for(i = 0; i < 6000; i++) {
		char expr[512];
		LIGOTimeGPS gps, prod, correct;
		double x;

		/* multiplicands */
		gps = random_LIGOTimeGPS();
		gps.gpsSeconds /= 100;
		x = pow(100., (rand() / (double) RAND_MAX) * 2. - 1.);

		/* compute product with lal */
		prod = gps;
		XLALGPSMultiply(&prod, x);

		/* compute product with bc */
		sprintf(expr, "%s * %.17g", print_LIGOTimeGPS(&gps), x);
		correct = compute_with_bc(expr);

		/* compare */
		if(fabs(XLALGPSDiff(&prod, &correct)) <= tolerance)
			continue;

		/* didn't match */
		printf("LAL said %s * %.17g = %s, but bc says = %s\n", print_LIGOTimeGPS(&gps), x, print_LIGOTimeGPS(&prod), print_LIGOTimeGPS(&correct));
		return 1;
	}

	return 0;
}
