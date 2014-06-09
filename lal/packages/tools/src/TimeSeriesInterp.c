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
	int kernel_length;
	double *cached_kernel;
	double residual;
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
	LALREAL8TimeSeriesInterp *interp;
	double *cached_kernel;

	if(kernel_length < 3)
		XLAL_ERROR_NULL(XLAL_EDOM);
	/* interpolator induces phase shifts unless this is odd */
	kernel_length -= (~kernel_length) & 1;

	interp = XLALMalloc(sizeof(*interp));
	cached_kernel = XLALMalloc(kernel_length * sizeof(*cached_kernel));
	if(!interp || !cached_kernel) {
		XLALFree(interp);
		XLALFree(cached_kernel);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	interp->series = series;
	interp->kernel_length = kernel_length;
	interp->cached_kernel = cached_kernel;
	/* >= 1 --> impossible.  forces kernel init on first eval */
	interp->residual = 2.;

	return interp;
}


/**
 * Free a LALREAL8TimeSeriesInterp object.  NULL is no-op.
 */


void XLALREAL8TimeSeriesInterpDestroy(LALREAL8TimeSeriesInterp *interp)
{
	if(interp)
		XLALFree(interp->cached_kernel);
	XLALFree(interp);
}


/**
 * Evaluate a LALREAL8TimeSeriesInterp at the LIGOTimeGPS t.  Raises a
 * XLAL_EDOM domain error if t is not in [epoch, epoch + length * deltaT)
 * where epoch, length, and deltaT are the start time, sample count, and
 * sample period, respectively, of the time series to which the
 * interpolator is attached.  A top hat-windowed sinc interpolating kernel
 * is used.
 */


REAL8 XLALREAL8TimeSeriesInterpEval(LALREAL8TimeSeriesInterp *interp, const LIGOTimeGPS *t)
{
	const double noop_threshold = 1./1024.;	/* samples */
	const REAL8 *data = interp->series->data->data;
	double *cached_kernel = interp->cached_kernel;
	REAL8 val = 0.0;
	double j = XLALGPSDiff(t, &interp->series->epoch) / interp->series->deltaT;
	double residual = round(j) - j;
	int start, stop;
	int i;

	if(j < 0 || j >= interp->series->data->length)
		XLAL_ERROR_REAL8(XLAL_EDOM);

	if(fabs(residual) < noop_threshold)
		return data[(int) round(j)];

	start = round(j - (interp->kernel_length - 1) / 2.0);
	stop = start + interp->kernel_length;

	if(fabs(residual - interp->residual) >= noop_threshold) {
		/* kernel is top hat-windowed sinc function.
		 *
		 *	x = pi (i - j);
		 *	kern = sin(x) / x
		 *
		 * don't check for 0/0:  can only occur if j is an integer,
		 * which is trapped by no-op path above.  note that
		 * argument of sin(x) increases by pi each iteration, so
		 * just need to compute its value for first iteration then
		 * flip sign for each subsequent iteration. */
		double sinx_over_pi = sin(LAL_PI * (start - j)) / LAL_PI;
		for(i = start; i < stop; i++, sinx_over_pi = -sinx_over_pi)
			*cached_kernel++ = sinx_over_pi / (i - j);
		interp->residual = residual;
		/* reset pointer */
		cached_kernel = interp->cached_kernel;
	}

	if(interp->series->data->length < (unsigned) stop)
		stop = interp->series->data->length;
	if(start < 0) {
		i = 0;
		cached_kernel -= start;
	} else
		i = start;
	for(; i < stop; i++)
		val += *cached_kernel++ * data[i];

	return val;
}
