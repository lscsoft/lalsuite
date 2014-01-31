/*
 * Copyright (C) 2014 Kipp Cannon
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


#include <lal/Date.h>
#include <lal/LALDatatypes.h>
#include <lal/TimeSeriesInterp.h>
#include <lal/XLALError.h>


struct tagLALREAL8TimeSeriesInterp {
	const REAL8TimeSeries *series;
	int kernel_length_over_2;
};


/**
 * Create a new REAL8TimeSeries interpolator associated with the given
 * REAL8TimeSeries object.  The kernel_length parameter sets the length of
 * the interpolating kernel in samples.  Use
 * XLALREAL8TimeSeriesInterpDestroy() to free.
 *
 * Notes:
 *
 * The LALREAL8TimeSeriesInterp object is opaque.  Calling code should not
 * attempt to inspect nor manipulate its contents directly.
 *
 * The REAL8TimeSeries object with which the LALREAL8TimeSeriesInterp is
 * associated must remain valid as long as the interpolator exists.  Do not
 * free the REAL8TimeSeries before freeing the LALREAL8TimeSeriesInterp.
 */


LALREAL8TimeSeriesInterp *XLALREAL8TimeSeriesInterpCreate(const REAL8TimeSeries *series, int kernel_length)
{
	LALREAL8TimeSeriesInterp *interp = XLALMalloc(sizeof(*interp));
	if(!interp)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	interp->series = series;
	interp->kernel_length_over_2 = kernel_length / 2.0;

	return interp;
}


/**
 * Free a LALREAL8TimeSeriesInterp object.  NULL is no-op.
 */


void XLALREAL8TimeSeriesInterpDestroy(LALREAL8TimeSeriesInterp *interp)
{
	XLALFree(interp);
}


/**
 * Evaluate a LALREAL8TimeSeriesInterp at the LIGOTimeGPS t.  Raises a
 * XLAL_EDOM domain error if t is not in [epoch, epoch + length * deltaT)
 * where epoch, length, and deltaT are the start time, sample count, and
 * sample period, respectively, of the time series to which the
 * interpolator is attached.  A Hann-windowed sinc interpolating kernel is
 * used.
 */


REAL8 XLALREAL8TimeSeriesInterpEval(LALREAL8TimeSeriesInterp *interp, const LIGOTimeGPS *t)
{
	const double noop_threshold = 1e-3;	/* samples */
	const REAL8 *data = interp->series->data->data;
	REAL8 val;
	double j = XLALGPSDiff(t, &interp->series->epoch) / interp->series->deltaT;
	int start, stop;
	double hann_rad_per_sample;
	int i;

	if(j < 0 || j >= interp->series->data->length)
		XLAL_ERROR_REAL8(XLAL_EDOM);

	if(fabs(round(j) - j) < noop_threshold)
		return data[(int) round(j)];

	start = round(j - interp->kernel_length_over_2);
	stop = round(j + interp->kernel_length_over_2) + 1;
	hann_rad_per_sample = LAL_TWOPI / (stop - 1 - start);

	val = 0.0;
	if(interp->series->data->length < (unsigned) stop)
		stop = interp->series->data->length;
	for(i = start < 0 ? 0 : start; i < stop; i++) {
		/* kernel is Hann-windowed sinc function.  don't check for
		 * 0/0:  can only occur if j is an integer, which is
		 * trapped by no-op path above */
		double x = LAL_PI * (i - j);
		double kernel = sin(x) / x * (0.5 - cos((i - start) * hann_rad_per_sample) / 2.0);
		val += kernel * data[i];
	}

	return val;
}
