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


struct tagLALREAL8SequenceInterp {
	const REAL8Sequence *s;
	int kernel_length;
	double welch_factor;
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
 * Create a new REAL8Sequence interpolator associated with the given
 * REAL8Sequence object.  The kernel_length parameter sets the length of
 * the interpolating kernel in samples.  Use
 * XLALREAL8SequenceInterpDestroy() to free.
 *
 * Notes:
 *
 * The LALREAL8SequenceInterp object is opaque.  Calling code should not
 * attempt to inspect nor manipulate its contents directly.
 *
 * The REAL8Sequence object with which the LALREAL8SequenceInterp is
 * associated must remain valid as long as the interpolator exists.  Do not
 * free the REAL8Sequence before freeing the LALREAL8SequenceInterp.
 */


LALREAL8SequenceInterp *XLALREAL8SequenceInterpCreate(const REAL8Sequence *s, int kernel_length)
{
	LALREAL8SequenceInterp *interp;
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

	interp->s = s;
	interp->kernel_length = kernel_length;
	interp->welch_factor = 1.0 / ((kernel_length - 1.) / 2. + 1.);
	interp->cached_kernel = cached_kernel;
	/* >= 1 --> impossible.  forces kernel init on first eval */
	interp->residual = 2.;
	/* set no-op threshold.  the kernel is recomputed when the residual
	 * changes by this much */
	interp->noop_threshold = 1. / (4 * interp->kernel_length);

	return interp;
}


/**
 * Free a LALREAL8SequenceInterp object.  NULL is no-op.
 */


void XLALREAL8SequenceInterpDestroy(LALREAL8SequenceInterp *interp)
{
	if(interp)
		XLALFree(interp->cached_kernel);
	XLALFree(interp);
}


/**
 * Evaluate a LALREAL8SequenceInterp at the real-valued index x.  Raises a
 * XLAL_EDOM domain error if x is not in [0, length) where length is the
 * sample count of the sequence which the interpolator is attached.
 *
 * A Welch-windowed sinc interpolating kernel is used.  See
 *
 * Smith, Julius O. Digital Audio Resampling Home Page
 * Center for Computer Research in Music and Acoustics (CCRMA), Stanford
 * University, 2014-01-10.  Web published at
 * http://www-ccrma.stanford.edu/~jos/resample/.
 *
 * for more information, but note that that reference uses a Kaiser window
 * for the sinc kernel's envelope whereas we use a Welch window here.  The
 * Welch window (inverted parabola) is chosen because it yields results
 * similar in accuracy to the Lanczos window but is much less costly to
 * compute.
 *
 * Be aware that for performance reasons the interpolating kernel is cached
 * and only recomputed if the error estimated to arise from failing to
 * recompute it exceeds the error estimated to arise from using a finite
 * interpolating kernel.  Therefore, if a function is interpolated at very
 * high resolution with a short kernel the result will consist of intervals
 * of constant values in a stair-step pattern.  The stair steps should be a
 * small contribution to the interpolation error but numerical
 * differentiation of the result is likely to be unsatisfactory.  In that
 * case, consider interpolating the derivative or use a longer kernel to
 * force more frequent kernel updates.
 */


REAL8 XLALREAL8SequenceInterpEval(LALREAL8SequenceInterp *interp, double x)
{
	const REAL8 *data = interp->s->data;
	double *cached_kernel = interp->cached_kernel;
	double *stop = cached_kernel + interp->kernel_length;
	/* split the real-valued sample index into integer and fractional
	 * parts.  the fractional part (residual) is the offset in samples
	 * from where we want to evaluate the function to where we know its
	 * value.  the interpolating kernel depends only on this quantity.
	 * when we compute a kernel, we record the value of this quantity,
	 * and only recompute the kernel if this quantity differs from the
	 * one for which the kernel was computed by more than the no-op
	 * threshold */
	int start = lround(x);
	double residual = start - x;
	REAL8 val;

	if(x < 0 || x >= interp->s->length)
		XLAL_ERROR_REAL8(XLAL_EDOM);

	if(fabs(residual) < interp->noop_threshold)
		return data[start];

	start -= (interp->kernel_length - 1) / 2;

	if(fabs(residual - interp->residual) >= interp->noop_threshold) {
		/* kernel is Welch-windowed sinc function.  the sinc
		 * component takes the form
		 *
		 *	x = pi (i - x);
		 *	kern = sin(x) / x
		 *
		 * we don't check for 0/0 because that can only occur if x
		 * is an integer, which is trapped by no-op path above.
		 * note that the  argument of sin(x) increases by pi each
		 * iteration, so we just need to compute its value for the
		 * first iteration then flip sign for each subsequent
		 * iteration.  for numerical reasons, it's better to
		 * compute sin(x) from residual rather than from (start -
		 * x), i.e. what it's argument should be for the first
		 * iteration, so we also have to figure out how many
		 * factors of -1 to apply to get its sign right for the
		 * first iteration.
		 */
		double welch_factor = interp->welch_factor;
		/* put a factor of welch_factor in this.  see below */
		double sinx_over_pi = sin(LAL_PI * residual) / LAL_PI * welch_factor;
		int i;
		if(interp->kernel_length & 2)
			sinx_over_pi = -sinx_over_pi;
		for(i = start; cached_kernel < stop; i++, sinx_over_pi = -sinx_over_pi) {
			double y = welch_factor * (i - x);
			if(fabs(y) < 1.)
				/* the window is
				 *
				 * sinx_over_pi / (i - x) * (1. - y * y)
				 *
				 * but by putting an extra factor of
				 * welch_factor into sinx_over_pi we can
				 * replace (i - x) with y, and then move
				 * the factor of 1/y into the parentheses
				 * to reduce the total number of arithmetic
				 * operations in the loop
				 */
				*cached_kernel++ = sinx_over_pi * (1. / y - y);
			else
				*cached_kernel++ = 0.;
		}
		interp->residual = residual;
		/* reset pointer */
		cached_kernel = interp->cached_kernel;
	}

	if(start + interp->kernel_length > (signed) interp->s->length)
		stop -= start + interp->kernel_length - interp->s->length;
	if(start < 0)
		cached_kernel -= start;
	else
		data += start;
	for(val = 0.0; cached_kernel < stop;)
		val += *cached_kernel++ * *data++;

	return val;
}


struct tagLALREAL8TimeSeriesInterp {
	const REAL8TimeSeries *series;
	LALREAL8SequenceInterp *seqinterp;
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
	LALREAL8SequenceInterp *seqinterp;

	interp = XLALMalloc(sizeof(*interp));
	seqinterp = XLALREAL8SequenceInterpCreate(series->data, kernel_length);
	if(!interp || !seqinterp) {
		XLALFree(interp);
		XLALREAL8SequenceInterpDestroy(seqinterp);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	interp->series = series;
	interp->seqinterp = seqinterp;

	return interp;
}


/**
 * Free a LALREAL8TimeSeriesInterp object.  NULL is no-op.
 */


void XLALREAL8TimeSeriesInterpDestroy(LALREAL8TimeSeriesInterp *interp)
{
	if(interp)
		XLALREAL8SequenceInterpDestroy(interp->seqinterp);
	XLALFree(interp);
}


/**
 * Evaluate a LALREAL8TimeSeriesInterp at the LIGOTimeGPS t.  Raises a
 * XLAL_EDOM domain error if t is not in [epoch, epoch + length * deltaT)
 * where epoch, length, and deltaT are the start time, sample count, and
 * sample period, respectively, of the time series to which the
 * interpolator is attached.
 *
 * A Welch-windowed sinc interpolating kernel is used.  See
 *
 * Smith, Julius O. Digital Audio Resampling Home Page
 * Center for Computer Research in Music and Acoustics (CCRMA), Stanford
 * University, 2014-01-10.  Web published at
 * http://www-ccrma.stanford.edu/~jos/resample/.
 *
 * for more information, but note that that reference uses a Kaiser window
 * for the sinc kernel's envelope whereas we use a Welch window here.  The
 * Welch window (inverted parabola) is chosen because it yields results
 * similar in accuracy to the Lanczos window but is much less costly to
 * compute.
 *
 * Be aware that for performance reasons the interpolating kernel is cached
 * and only recomputed if the error estimated to arise from failing to
 * recompute it exceeds the error estimated to arise from using a finite
 * interpolating kernel.  Therefore, if a function is interpolated at very
 * high resolution with a short kernel the result will consist of intervals
 * of constant values in a stair-step pattern.  The stair steps should be a
 * small contribution to the interpolation error but numerical
 * differentiation of the result is likely to be unsatisfactory.  In that
 * case, consider interpolating the derivative or use a longer kernel to
 * force more frequent kernel updates.
 */


REAL8 XLALREAL8TimeSeriesInterpEval(LALREAL8TimeSeriesInterp *interp, const LIGOTimeGPS *t)
{
	return XLALREAL8SequenceInterpEval(interp->seqinterp, XLALGPSDiff(t, &interp->series->epoch) / interp->series->deltaT);
}
