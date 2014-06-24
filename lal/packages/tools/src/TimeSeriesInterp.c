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
#include <lal/Window.h>
#include <lal/XLALError.h>


#define KAISER_BETA	(1.0 * LAL_PI)


struct tagLALREAL8TimeSeriesInterp {
	const REAL8TimeSeries *series;
	int kernel_length;
	double *kaiser_window;
	double *cached_kernel;
	double residual;
	/* samples.  the length of the kernel sets the bandwidth of the
	 * interpolator:  the longer the kernel, the closer to an ideal
	 * interpolator it becomes.  we tie the interval at which the
	 * kernel is regenerated to this in a heuristic way to hide the
	 * sub-sample residual quantization in the filter's roll-off. */
	double noop_threshold;
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
	REAL8Window *kaiser_window;
	double *cached_kernel;

	if(kernel_length < 3)
		XLAL_ERROR_NULL(XLAL_EDOM);
	/* interpolator induces phase shifts unless this is odd */
	kernel_length -= (~kernel_length) & 1;

	interp = XLALMalloc(sizeof(*interp));
	/* we need the window to be centred on (kernel_length-1)/2.  LAL's
	 * window functions do this. */
	kaiser_window = XLALCreateKaiserREAL8Window(kernel_length, KAISER_BETA);
	cached_kernel = XLALMalloc(kernel_length * sizeof(*cached_kernel));
	if(!interp || !kaiser_window || !cached_kernel) {
		XLALFree(interp);
		XLALDestroyREAL8Window(kaiser_window);
		XLALFree(cached_kernel);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	interp->series = series;
	interp->kernel_length = kernel_length;
	/* grab the data pointer from the REAL8Window object */
	interp->kaiser_window = kaiser_window->data->data;
	kaiser_window->data->data = NULL;
	XLALDestroyREAL8Window(kaiser_window);
	interp->cached_kernel = cached_kernel;
	/* >= 1 --> impossible.  forces kernel init on first eval */
	interp->residual = 2.;
	/* set no-op threshold.  the kernel is recomputed when the residual
	 * changes by this much */
	interp->noop_threshold = 1. / (4 * interp->kernel_length);

	return interp;
}


/**
 * Free a LALREAL8TimeSeriesInterp object.  NULL is no-op.
 */


void XLALREAL8TimeSeriesInterpDestroy(LALREAL8TimeSeriesInterp *interp)
{
	if(interp) {
		XLALFree(interp->kaiser_window);
		XLALFree(interp->cached_kernel);
	}
	XLALFree(interp);
}


/**
 * Evaluate a LALREAL8TimeSeriesInterp at the LIGOTimeGPS t.  Raises a
 * XLAL_EDOM domain error if t is not in [epoch, epoch + length * deltaT)
 * where epoch, length, and deltaT are the start time, sample count, and
 * sample period, respectively, of the time series to which the
 * interpolator is attached.
 *
 * A Kaiser-windowed (beta = KAISER_BETA) sinc interpolating kernel is
 * used.  See
 *
 * Smith, Julius O. Digital Audio Resampling Home Page
 * Center for Computer Research in Music and Acoustics (CCRMA), Stanford
 * University, 2014-01-10.  Web published at
 * http://www-ccrma.stanford.edu/~jos/resample/.
 *
 * for more information.
 */


REAL8 XLALREAL8TimeSeriesInterpEval(LALREAL8TimeSeriesInterp *interp, const LIGOTimeGPS *t)
{
	const REAL8 *data = interp->series->data->data;
	double *cached_kernel = interp->cached_kernel;
	/* the (real-valued) sample index at which we wish to evalute the
	 * source time series */
	double j = XLALGPSDiff(t, &interp->series->epoch) / interp->series->deltaT;
	/* split the real-valued sample index into integer and fractional
	 * parts.  the fractional part (residual) is the offset in samples
	 * from where we want to evaluate the function to where we know its
	 * value.  the interpolating kernel depends only on this quantity.
	 * when we compute a kernel, we record the value of this quantity,
	 * and only recompute the kernel if this quantity differs from the
	 * one for which the kernel was computed by more than the no-op
	 * threshold */
	int start = lround(j);
	int stop;
	double residual = start - j;
	int i;
	REAL8 val;

	if(j < 0 || j >= interp->series->data->length)
		XLAL_ERROR_REAL8(XLAL_EDOM);

	if(fabs(residual) < interp->noop_threshold)
		return data[start];

	start -= (interp->kernel_length - 1) / 2;
	stop = start + interp->kernel_length;

	if(fabs(residual - interp->residual) >= interp->noop_threshold) {
		/* kernel is Kaiser-windowed sinc function.  we don't
		 * bother re-computing the Kaiser window, we consider it to
		 * be approximately independent of the sub-sample shift.
		 * only the sinc component is recomputed, and it takes the
		 * form
		 *
		 *	x = pi (i - j);
		 *	kern = sin(x) / x
		 *
		 * we don't check for 0/0 because that can only occur if j
		 * is an integer, which is trapped by no-op path above.
		 * note that the  argument of sin(x) increases by pi each
		 * iteration, so we just need to compute its value for the
		 * first iteration then flip sign for each subsequent
		 * iteration.  for numerical reasons, it's better to
		 * compute sin(x) from residual rather than from (start -
		 * j), i.e. what it's argument should be for the first
		 * iteration, so we also have to figure out how many
		 * factors of -1 to apply to get its sign right for the
		 * first iteration.
		 */
		const double *kaiser_window = interp->kaiser_window;
		double sinx_over_pi = sin(LAL_PI * residual) / LAL_PI;
		if(interp->kernel_length & 2)
			sinx_over_pi = -sinx_over_pi;
		for(i = start; i < stop; i++, sinx_over_pi = -sinx_over_pi)
			*cached_kernel++ = sinx_over_pi / (i - j) * *kaiser_window++;
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
	for(val = 0.0; i < stop; i++)
		val += *cached_kernel++ * data[i];

	return val;
}
