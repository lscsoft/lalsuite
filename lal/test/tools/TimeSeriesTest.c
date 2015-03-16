/*
*  Copyright (C) 2007 Bernd Machenschalk, Jolien Creighton, Kipp Cannon
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

/*
 * Note that low-level sequence manipulation is tested in the Sequence
 * verification code.  The object here is mostly to test meta-data
 * manipulation eg. epoch manipulation.
 */

#include <stdio.h>
#include <stdlib.h>

#include <lal/Date.h>
#include <lal/LALDatatypes.h>
#include <lal/Sequence.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>

static LIGOTimeGPS gps_zero = LIGOTIMEGPSZERO;


static REAL4TimeSeries *random_timeseries(size_t n)
{
	REAL4TimeSeries *s = XLALCreateREAL4TimeSeries("blah", &gps_zero, 0.0, 1.0 / n, &lalDimensionlessUnit, n);

	while(n--)
		s->data->data[n] = rand() / (float) RAND_MAX;
	return(s);
}


static INT4TimeSeries *sequential_timeseries(size_t n)
{
	INT4TimeSeries *s = XLALCreateINT4TimeSeries("blah", &gps_zero, 0.0, 1.0 / n, &lalDimensionlessUnit, n);

	while(n--)
		s->data->data[n] = n;
	return(s);
}


static int cmp_real4(REAL4 *a, REAL4 *b, size_t n)
{
	REAL4 d;

	while(n--) {
		d = *a - *b;
		if(d < 0.0)
			return(-1);
		if(d > 0.0)
			return(+1);
	}
	return(0);
}


int main(void)
{
	REAL4TimeSeries *x, *y;
	INT4TimeSeries *a;
	int i;

	/*
	 * Destroy
	 */

	/* NULL pointer */
	XLALDestroyREAL4TimeSeries(NULL);

	/* Incompletely/Incorrectly initialized structure */
	x = XLALCalloc(1, sizeof(*x));
	XLALDestroyREAL4TimeSeries(x);
	x = XLALCalloc(1, sizeof(*x));
	x->data = XLALCalloc(1, sizeof(*x->data));
	XLALDestroyREAL4TimeSeries(x);

	/*
	 * Create
	 */

	/* try segfaulting on array access */
	x = XLALCreateREAL4TimeSeries("blah", &gps_zero, 0.0, 1.0 / 1, &lalDimensionlessUnit, 1);
	x->data->data[0] = 1.0;
	if((x->f0 != 0.0) || (x->deltaT != 1.0) || (x->data->length != 1) || (x->data->data[0] != 1.0)) {
		fprintf(stderr, "Create test 1 failed\n");
		exit(1);
	}
	XLALDestroyREAL4TimeSeries(x);

	/*
	 * Cut
	 */

	/* check metadata */
	x = random_timeseries(1024);
	y = XLALCutREAL4TimeSeries(x, 256, 512);
	if((y->deltaT != x->deltaT) || (y->f0 != x->f0) || cmp_real4(x->data->data + 256, y->data->data, 512) || (XLALGPSDiff(&y->epoch, &x->epoch) != 256 * x->deltaT)) {
		fprintf(stderr, "Cut test 1 failed\n");
		exit(1);
	}
	XLALDestroyREAL4TimeSeries(x);
	XLALDestroyREAL4TimeSeries(y);

	/*
	 * Resize
	 */

	/* check metadata */
	a = sequential_timeseries(1024);
	XLALResizeINT4TimeSeries(a, 256, 512);
	for(i = 0; i < (int) a->data->length; i++)
		if(a->data->data[i] != i + 256) {
			fprintf(stderr, "Resize test 1a failed\n");
			exit(1);
		}
	if((a->data->length != 512) || (XLALGPSDiff(&a->epoch, &gps_zero) != 256 * a->deltaT)) {
		fprintf(stderr, "Resize test 1b failed\n");
		exit(1);
	}
	XLALDestroyINT4TimeSeries(a);

	/*
	 * Success
	 */

	exit(0);
}
