/*
 * Copyright (C) 2014,2015 Kipp Cannon
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA  02111-1307  USA
 */


#include <math.h>
#include <stdio.h>


#include <lal/Date.h>
#include <lal/LALDatatypes.h>
#include <lal/TimeSeries.h>
#include <lal/TimeSeriesInterp.h>
#include <lal/Units.h>


static LIGOTimeGPS gps_zero = LIGOTIMEGPSZERO;


static REAL8TimeSeries *new_series(double deltaT, unsigned length, double fill)
{
	REAL8TimeSeries *new = XLALCreateREAL8TimeSeries("blah", &gps_zero, 0.0, deltaT, &lalDimensionlessUnit, length);
	unsigned i;
	for(i = 0 ; i < new->data->length; i++)
		new->data->data[i] = fill;
	return new;
}


static REAL8TimeSeries *copy_series(const REAL8TimeSeries *src)
{
	return XLALCutREAL8TimeSeries(src, 0, src->data->length);
}


static LIGOTimeGPS t_i(const REAL8TimeSeries *s, unsigned i)
{
	LIGOTimeGPS t = s->epoch;
	XLALGPSAdd(&t, i * s->deltaT);
	return t;
}


static void add_sine(REAL8TimeSeries *s, LIGOTimeGPS epoch, double ampl, double freq)
{
	double t0 = XLALGPSDiff(&s->epoch, &epoch);
	unsigned i;

	for(i = 0; i < s->data->length; i++) {
		double t = t0 + i * s->deltaT;
		s->data->data[i] = ampl * sin(LAL_TWOPI * freq * t);
	}
}


static void evaluate(REAL8TimeSeries *dst, LALREAL8TimeSeriesInterp *interp, int bounds_check)
{
	unsigned i;

	for(i = 0; i < dst->data->length; i++) {
		LIGOTimeGPS t = t_i(dst, i);
		dst->data->data[i] = XLALREAL8TimeSeriesInterpEval(interp, &t, bounds_check);
	}
}


static REAL8TimeSeries *error(const REAL8TimeSeries *s1, const REAL8TimeSeries *s0)
{
	REAL8TimeSeries *result = copy_series(s1);
	unsigned i;

	/* assumes series are compatible */

	for(i = 0; i < s1->data->length; i++)
		result->data->data[i] -= s0->data->data[i];

	return result;
}


static double RMS(const REAL8TimeSeries *s)
{
	double rms = 0.;
	unsigned i;

	for(i = 0; i < s->data->length; i++)
		rms += s->data->data[i] * s->data->data[i];

	return sqrt(rms / s->data->length);
}


static void minmax(const REAL8TimeSeries *s, double *min, double *max)
{
	unsigned i;

	*min = *max = s->data->data[0];
	for(i = 1; i < s->data->length; i++) {
		if(s->data->data[i] < *min)
			*min = s->data->data[i];
		if(s->data->data[i] > *max)
			*max = s->data->data[i];
	}
}


static void check_result(const REAL8TimeSeries *model, const REAL8TimeSeries *result, double rms_bound, double residual_min, double residual_max)
{
	REAL8TimeSeries *err = error(model, result);
	double rms, min, max;
	rms = RMS(err);
	minmax(err, &min, &max);

	fprintf(stderr, "error vector:  RMS=%g, min=%g, max=%g\n", rms, min, max);
	XLALDestroyREAL8TimeSeries(err);
	if(rms > rms_bound || min < residual_min || max > residual_max) {
		fprintf(stderr, "error vector larger than allowed\n");
		exit(1);
	}
}


int main(void)
{
	REAL8TimeSeries *src, *dst, *mdl;
	LALREAL8TimeSeriesInterp *interp;
	double f;

	/*
	 * single frequency.  include many cycles to approximate source
	 * data with infinite extent.  interpolate 1 cycle at high
	 * resolution and compare to model.
	 */

	f = 4000.;	/* close to 50% of Nyquist */

	src = new_series(1.0 / 16384, 1024 * 1024, 0.0);
	add_sine(src, src->epoch, 1.0, f);

	mdl = new_series(1. / 1e8, round(1. / f * 1e8), 0.0);
	XLALGPSAdd(&mdl->epoch, src->data->length * src->deltaT * .4);
	dst = copy_series(mdl);

	fprintf(stderr, "interpolating unit amplitude %g kHz sine function sampled at %g Hz to %g MHz\n", f / 1000., 1.0 / src->deltaT, 1.0 / dst->deltaT / 1e6);

	add_sine(mdl, src->epoch, 1.0, f);

	/* in the interpolator's current incarnation, a kernel length this
	 * size at a sample rate of 16384 Hz leads to an internal "no-op
	 * threshold" of about 0.2 ns --- if two samples have different GPS
	 * times then they get their own interpolating kernels, i.e. we
	 * defeat the kernel caching mechanism so that we can watch the
	 * behaviour of the kernel sample-by-sample */
	interp = XLALREAL8TimeSeriesInterpCreate(src, 65535);
	evaluate(dst, interp, 1);
	XLALREAL8TimeSeriesInterpDestroy(interp);

#if 0
	{
	unsigned i;
	FILE *output = fopen("input.txt", "w");
	for(i = 0; i < src->data->length; i++) {
		LIGOTimeGPS t = t_i(src, i);
		fprintf(output, "%u.%09u %.16g\n", t.gpsSeconds, t.gpsNanoSeconds, src->data->data[i]);
	}
	fclose(output);
	output = fopen("output.txt", "w");
	for(i = 0; i < mdl->data->length; i++) {
		LIGOTimeGPS t = t_i(mdl, i);
		fprintf(output, "%u.%09u %.16g %.16g\n", t.gpsSeconds, t.gpsNanoSeconds, mdl->data->data[i], dst->data->data[i]);
	}
	fclose(output);
	}
#endif

	check_result(mdl, dst, 1.3e-10, -3.6e-10, +3.6e-10);

	XLALDestroyREAL8TimeSeries(src);
	XLALDestroyREAL8TimeSeries(dst);
	XLALDestroyREAL8TimeSeries(mdl);

	/*
	 * higher single frequency, shorter filter, to check for phase
	 * error.  interpolate three cycles this time.
	 */

	f = 6500.;	/* about 80% of Nyquist */

	src = new_series(1.0 / 16384, 256, 0.0);
	add_sine(src, src->epoch, 1.0, f);

	mdl = new_series(1. / 1e8, round(3. / f * 1e8), 0.0);
	XLALGPSAdd(&mdl->epoch, src->data->length * src->deltaT * .4);
	dst = copy_series(mdl);

	fprintf(stderr, "interpolating unit amplitude %g kHz sine function sampled at %g Hz to %g MHz\n", f / 1000., 1.0 / src->deltaT, 1.0 / dst->deltaT / 1e6);

	add_sine(mdl, src->epoch, 1.0, f);

	interp = XLALREAL8TimeSeriesInterpCreate(src, 9);
	evaluate(dst, interp, 1);
	XLALREAL8TimeSeriesInterpDestroy(interp);

#if 0
	{
	unsigned i;
	FILE *output = fopen("input.txt", "w");
	for(i = 0; i < src->data->length; i++) {
		LIGOTimeGPS t = t_i(src, i);
		fprintf(output, "%u.%09u %.16g\n", t.gpsSeconds, t.gpsNanoSeconds, src->data->data[i]);
	}
	fclose(output);
	output = fopen("output.txt", "w");
	for(i = 0; i < mdl->data->length; i++) {
		LIGOTimeGPS t = t_i(mdl, i);
		fprintf(output, "%u.%09u %.16g %.16g\n", t.gpsSeconds, t.gpsNanoSeconds, mdl->data->data[i], dst->data->data[i]);
	}
	fclose(output);
	}
#endif

	check_result(mdl, dst, 0.03, -0.078, +0.083);

	XLALDestroyREAL8TimeSeries(src);
	XLALDestroyREAL8TimeSeries(dst);
	XLALDestroyREAL8TimeSeries(mdl);

	/*
	 * test behaviour in last sample.  allocate series 1 sample longer
	 * than we need it to be so we can control the value of the data
	 * beyond the end of the array, and make sure the interpolator
	 * doesn't return that value.
	 */

	src = new_series(1.0 / 16384, 257, 2.0);/* all 2s */
	src->data->length--;			/* fake the length */
	src->data->data[src->data->length] = 1;	/* place a 1 beyond the end */

	interp = XLALREAL8TimeSeriesInterpCreate(src, 9);
	{
	LIGOTimeGPS t = src->epoch;
	double result;
	XLALGPSAdd(&t, src->data->length * src->deltaT);
	/* evalute at time of sample beyond end.  this should fail */
	fprintf(stderr, "checking for out-of-bounds failure ...\n");
	result = XLALREAL8TimeSeriesInterpEval(interp, &t, 1);
	if(!XLAL_IS_REAL8_FAIL_NAN(result)) {
		fprintf(stderr, "error:  interpolator failed to report error beyond end of array\n");
		exit(1);
	} else
		fprintf(stderr, "... passed\n");
	/* these should work and return 0., not 1. and not 2 */
	result = XLALREAL8TimeSeriesInterpEval(interp, &t, 0);
	if(result != 0.) {
		fprintf(stderr, "error:  interpolator failed in final sample (expected %.16g got %.16g)\n", 0., result);
		exit(1);
	}
	XLALGPSAdd(&t, -1e-9);
	result = XLALREAL8TimeSeriesInterpEval(interp, &t, 1);
	if(result != 0.) {
		fprintf(stderr, "error:  interpolator failed in final sample (expected %.16g got %.16g)\n", 0., result);
		exit(1);
	}
	}
	XLALREAL8TimeSeriesInterpDestroy(interp);

	XLALDestroyREAL8TimeSeries(src);

	/*
	 * success
	 */

	exit(0);
}
