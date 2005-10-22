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


int main(int argc, char *argv[])
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
	x = malloc(sizeof(*x));
	x->data = NULL;
	XLALDestroyREAL4TimeSeries(x);
	x = malloc(sizeof(*x));
	x->data = malloc(sizeof(*x->data));
	x->data->data = NULL;
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
	if((y->deltaT != x->deltaT) || (y->f0 != x->f0) || cmp_real4(x->data->data + 256, y->data->data, 512) || (XLALGPSDiff(&x->epoch, &y->epoch) != 256 * x->deltaT)) {
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
	if((a->data->length != 512) || (XLALGPSDiff(&gps_zero, &a->epoch) != 256 * a->deltaT)) {
		fprintf(stderr, "Resize test 1b failed\n");
		exit(1);
	}
	XLALDestroyINT4TimeSeries(a);

	/*
	 * Success
	 */

	exit(0);
}
