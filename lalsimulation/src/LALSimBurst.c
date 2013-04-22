/*
 * Copyright (C) 2008 J. Creighton, K. Cannon
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


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */


#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <lal/LALSimBurst.h>
#include <lal/LALConstants.h>
#include <lal/FrequencySeries.h>
#include <lal/Sequence.h>
#include <lal/TimeFreqFFT.h>
#include <lal/TimeSeries.h>
#include <lal/RealFFT.h>
#include <lal/Units.h>
#include <lal/Date.h>
#include "check_series_macros.h"

/*
 * ============================================================================
 *
 *          Fill a time series with stationary white Gaussian noise
 *
 * ============================================================================
 */


static void gaussian_noise(REAL8TimeSeries * series, REAL8 rms, gsl_rng * rng)
{
	unsigned i;

	for(i = 0; i < series->data->length; i++)
		series->data->data[i] = gsl_ran_gaussian(rng, rms);
}


/*
 * ============================================================================
 *
 *                                 Utilities
 *
 * ============================================================================
 */


/**
 * Returns the strain of the sample with the largest magnitude.
 */


REAL8 XLALMeasureHPeak(const REAL8TimeSeries *series)
{
	double hpeak;
	unsigned i;

	if(!series->data->length) {
		XLALPrintError("%s(): length must be > 0\n", __func__);
		XLAL_ERROR_REAL8(XLAL_EBADLEN);
	}

	hpeak = series->data->data[0];
	for(i = 1; i < series->data->length; i++)
		if(fabs(series->data->data[i]) > fabs(hpeak))
			hpeak = series->data->data[i];

	return hpeak;
}


/**
 * From two time series, s1 and s2, computes and returns
 *
 * \f$\int s1(t) s2(t) d t\f$
 */


REAL8 XLALMeasureIntS1S2DT(const REAL8TimeSeries *s1, const REAL8TimeSeries *s2)
{
	double e = 0.0;
	double sum = 0.0;
	unsigned i;

	/* FIXME:  this is overly strict, this function could be smarter */

	LAL_CHECK_CONSISTENT_TIME_SERIES(s1, s2, XLAL_REAL8_FAIL_NAN);

	/* Kahans's compensated summation algorithm */

	for(i = 0; i < s1->data->length; i++) {
		double tmp = sum;
		/* what we want to add = s1 * s2 + "error from last
		 * iteration" */
		double x = s1->data->data[i] * s2->data->data[i] + e;
		/* add */
		sum += x;
		/* negative of what was actually added */
		e = tmp - sum;
		/* what didn't get added, add next time */
		e += x;
	}

	return sum * s1->deltaT;
}


/**
 * Returns what people call the "root-sum-square strain".  Infact, this is
 *
 * \f$\sqrt{\sum (h_{+}^{2} + h_{x}^{2}) \Delta t},\f$
 *
 * which is an approximation of
 *
 * \f$\sqrt{\int (h_{+}^{2} + h_{x}^{2}) d t}.\f$
 */


REAL8 XLALMeasureHrss(
	const REAL8TimeSeries *hplus,
	const REAL8TimeSeries *hcross
)
{
	return sqrt(XLALMeasureIntS1S2DT(hplus, hplus) + XLALMeasureIntS1S2DT(hcross, hcross));
}


/**
 * Given the Fourier transform of a real-valued function h(t), compute and
 * return the integral of the square of its derivative:
 *
 * \f$\int \stackrel{.}{h}^{2} d t\f$.
 *
 * The normalization factors in this function assume that
 * XLALREAL8FreqTimeFFT() will be used to convert the frequency series to
 * the time domain.
 */


REAL8 XLALMeasureIntHDotSquaredDT(const COMPLEX16FrequencySeries *fseries)
{
	unsigned i;
	double e = 0.0;
	double sum = 0.0;

	/* Kahan's compensated summation algorithm. The summation is done
	 * from lowest to highest frequency under the assumption that high
	 * frequency components tend to add more to the magnitude of the
	 * derivative.  Note that because only half the components of the
	 * Fourier transform are stored a factor of 2 is added after the
	 * sum.  The DC component should only count once, but it does not
	 * contribute anything to the sum so no special case is required to
	 * handle it. */

	for(i = 0; i < fseries->data->length; i++) {
		double tmp = sum;
		/* what we want to add = f^{2} |\tilde{s}(f)|^{2} + "error
		 * from last iteration" */
		double x = pow(fseries->f0 + i * fseries->deltaF, 2) * pow(cabs(fseries->data->data[i]), 2) + e;
		/* add */
		sum += x;
		/* negative of what was actually added */
		e = tmp - sum;
		/* what didn't get added, add next time */
		e += x;
	}

	/* because we've only summed the positive frequency components */

	sum *= 2;

	/* 4 \pi^{2} \delta f */

	sum *= LAL_TWOPI * LAL_TWOPI * fseries->deltaF;

	return sum;
}


/**
 * Given h+ and hx in the waveframe, compute and return E/r^2.  The return
 * value is in LAL's native units, computed by evaluating
 *
 * \f$\int [ \stackrel{.}{h}_{+}^{2} + \stackrel{.}{h}_{\times}^{2} ] d t\f$
 *
 * and multiplying by LAL_C_SI\f$^{3}\f$ / (4 LAL_G_SI).
 */


REAL8 XLALMeasureEoverRsquared(REAL8TimeSeries *hplus, REAL8TimeSeries *hcross)
{
	REAL8FFTPlan *plan;
	COMPLEX16FrequencySeries *tilde_hplus, *tilde_hcross;
	double e_over_rsquared;
	unsigned i;

	/* FIXME:  this is overly strict, this function could be smarter */

	LAL_CHECK_CONSISTENT_TIME_SERIES(hplus, hcross, XLAL_REAL8_FAIL_NAN);

	/* transform to the frequency domain */

	plan = XLALCreateForwardREAL8FFTPlan(hplus->data->length, 0);
	tilde_hplus = XLALCreateCOMPLEX16FrequencySeries(NULL, &hplus->epoch, 0.0, 0.0, &lalDimensionlessUnit, hplus->data->length / 2 + 1);
	tilde_hcross = XLALCreateCOMPLEX16FrequencySeries(NULL, &hcross->epoch, 0.0, 0.0, &lalDimensionlessUnit, hcross->data->length / 2 + 1);
	if(!plan || !tilde_hplus || !tilde_hcross) {
		XLALDestroyCOMPLEX16FrequencySeries(tilde_hplus);
		XLALDestroyCOMPLEX16FrequencySeries(tilde_hcross);
		XLALDestroyREAL8FFTPlan(plan);
		XLAL_ERROR(XLAL_EFUNC);
	}
	i = XLALREAL8TimeFreqFFT(tilde_hplus, hplus, plan);
	i |= XLALREAL8TimeFreqFFT(tilde_hcross, hcross, plan);
	XLALDestroyREAL8FFTPlan(plan);
	if(i) {
		XLALDestroyCOMPLEX16FrequencySeries(tilde_hplus);
		XLALDestroyCOMPLEX16FrequencySeries(tilde_hcross);
		XLAL_ERROR(XLAL_EFUNC);
	}

	/* measure E / r^2 */

	e_over_rsquared = (double) LAL_C_SI * LAL_C_SI * LAL_C_SI / (4 * LAL_G_SI) * (XLALMeasureIntHDotSquaredDT(tilde_hplus) + XLALMeasureIntHDotSquaredDT(tilde_hcross));

	/* done */

	XLALDestroyCOMPLEX16FrequencySeries(tilde_hplus);
	XLALDestroyCOMPLEX16FrequencySeries(tilde_hcross);

	return e_over_rsquared;
}


/*
 * ============================================================================
 *
 *                   Construct a \delta Function Injection
 *
 * ============================================================================
 */


/**
 * Places a single non-zero sample into the middle of the time series.
 *
 * Parameters:
 *
 * hpeak:  amplitude of the impulse.
 *
 * Output:
 *
 * The h+ and hx time series both have an odd number of samples all set to
 * 0 except for a single sample with amplitude hpeak in the middle of the
 * h+ time series.
 */


int XLALGenerateImpulseBurst(
	REAL8TimeSeries **hplus,
	REAL8TimeSeries **hcross,
	REAL8 hpeak,
	REAL8 delta_t
)
{
	int length;
	LIGOTimeGPS epoch;

	/* length is 1353 samples, because it's 13:53 right now and it's an
	 * odd integer */

	length = 1353;

	/* the middle sample is t = 0 */

	XLALGPSSetREAL8(&epoch, -(length - 1) / 2 * delta_t);

	/* allocate the time series */

	*hplus = XLALCreateREAL8TimeSeries("Impulse +", &epoch, 0.0, delta_t, &lalStrainUnit, length);
	*hcross = XLALCreateREAL8TimeSeries("Impulse x", &epoch, 0.0, delta_t, &lalStrainUnit, length);
	if(!*hplus || !*hcross) {
		XLALDestroyREAL8TimeSeries(*hplus);
		XLALDestroyREAL8TimeSeries(*hcross);
		*hplus = *hcross = NULL;
		XLAL_ERROR(XLAL_EFUNC);
	}

	/* set to zero */

	memset((*hplus)->data->data, 0, length * sizeof(*(*hplus)->data->data));
	memset((*hcross)->data->data, 0, length * sizeof(*(*hcross)->data->data));

	/* put impulse into middle sample of h+ */

	(*hplus)->data->data[(length - 1) / 2] = hpeak;

	/* done */

	return 0;
}


/*
 * ============================================================================
 *
 *            Construct a Band- and Time-Limited White Noise Burst
 *
 * ============================================================================
 */


/**
 * Parameters:
 *
 * duration
 * 	time domain Gaussian envelope is \f$\propto \exp ( -\frac{1}{2} t^{2} / duration^{2} )\f$
 *      where t and duration are in seconds.
 * frequency
 * bandwidth
 * 	frequency domain Gaussian envelope is \f$\propto \exp ( -\frac{1}{2} (f - f_{0})^{2} / bandwidth^{2} )\f$
 *      where f and bandwidth are in Hertz.
 * int_hdot_squared
 * 	waveform is normalized so that \f$\int (\stackrel{.}{h}_{+}^{2} + \stackrel{.}{h}_{\times}^{2}) d t\f$
 *       equals this
 * delta_t
 * 	the sample rate of the time series to construct
 * rng
 * 	a GSL random number generator to be used to produce Gaussian random
 * 	variables
 *
 * Output:
 *
 * Two time series containing h+(t) and hx(t), with the time-domain
 * Gaussian envelope's peak located at t = 0 (as defined by the epoch and
 * deltaT).  The + and x time series are two independent injections.
 *
 * Note:  because the injection is constructed with a random number
 * generator, any changes to this function that change how random numbers
 * are chosen will indirectly have the effect of altering the relationship
 * between injection waveform and random number seed.  For example,
 * increasing the length of the time series will change the injection
 * waveforms.  There's nothing wrong with this, the waveforms are still
 * correct, but if there is a need to reproduce a waveform exactly then it
 * will be necessary to tag the code before making such changes.
 */


int XLALGenerateBandAndTimeLimitedWhiteNoiseBurst(
	REAL8TimeSeries **hplus,
	REAL8TimeSeries **hcross,
	REAL8 duration,
	REAL8 frequency,
	REAL8 bandwidth,
	REAL8 int_hdot_squared,
	REAL8 delta_t,
	gsl_rng *rng
)
{
	int length;
	LIGOTimeGPS epoch;
	COMPLEX16FrequencySeries *tilde_hplus, *tilde_hcross;
	REAL8Window *window;
	REAL8FFTPlan *plan;
	REAL8 norm_factor;
	/* compensate the width of the time-domain window's envelope for
	 * the broadening will be induced by the subsequent application of
	 * the frequency-domain envelope */
	REAL8 sigma_t_squared = duration * duration / 4.0 - 1.0 / (LAL_PI * LAL_PI * bandwidth * bandwidth);
	unsigned i;

	/* check input.  checking if sigma_t_squared < 0 is equivalent to
	 * checking if duration * bandwidth < LAL_2_PI */

	if(duration < 0 || bandwidth < 0 || sigma_t_squared < 0 || int_hdot_squared < 0 || delta_t <= 0) {
		XLALPrintError("%s(): invalid input parameters\n", __func__);
		*hplus = *hcross = NULL;
		XLAL_ERROR(XLAL_EINVAL);
	}

	/* length of the injection time series is 30 * duration, rounded to
	 * the nearest odd integer */

	length = (int) floor(30.0 * duration / delta_t / 2.0);
	length = 2 * length + 1;

	/* the middle sample is t = 0 */

	XLALGPSSetREAL8(&epoch, -(length - 1) / 2 * delta_t);

	/* allocate the time series */

	*hplus = XLALCreateREAL8TimeSeries("BTLWNB +", &epoch, 0.0, delta_t, &lalStrainUnit, length);
	*hcross = XLALCreateREAL8TimeSeries("BTLWNB x", &epoch, 0.0, delta_t, &lalStrainUnit, length);
	if(!*hplus || !*hcross) {
		XLALDestroyREAL8TimeSeries(*hplus);
		XLALDestroyREAL8TimeSeries(*hcross);
		*hplus = *hcross = NULL;
		XLAL_ERROR(XLAL_EFUNC);
	}

	/* fill with independent zero-mean unit variance Gaussian random
	 * numbers (any non-zero amplitude is OK, it will be adjusted
	 * later) */

	gaussian_noise(*hplus, 1, rng);
	gaussian_noise(*hcross, 1, rng);

	/* apply the time-domain Gaussian window.  the window function's
	 * shape parameter is ((length - 1) * delta_t / 2) / \sigma_{t} where
	 * \sigma_{t} is the compensated time-domain window duration */

	window = XLALCreateGaussREAL8Window((*hplus)->data->length, (((*hplus)->data->length - 1) * delta_t / 2) / sqrt(sigma_t_squared));
	if(!window) {
		XLALDestroyREAL8TimeSeries(*hplus);
		XLALDestroyREAL8TimeSeries(*hcross);
		*hplus = *hcross = NULL;
		XLAL_ERROR(XLAL_EFUNC);
	}
	for(i = 0; i < window->data->length; i++) {
		(*hplus)->data->data[i] *= window->data->data[i];
		(*hcross)->data->data[i] *= window->data->data[i];
	}
	XLALDestroyREAL8Window(window);

	/* transform to the frequency domain */

	plan = XLALCreateForwardREAL8FFTPlan((*hplus)->data->length, 0);
	tilde_hplus = XLALCreateCOMPLEX16FrequencySeries(NULL, &epoch, 0.0, 0.0, &lalDimensionlessUnit, (*hplus)->data->length / 2 + 1);
	tilde_hcross = XLALCreateCOMPLEX16FrequencySeries(NULL, &epoch, 0.0, 0.0, &lalDimensionlessUnit, (*hcross)->data->length / 2 + 1);
	if(!plan || !tilde_hplus || !tilde_hcross) {
		XLALDestroyCOMPLEX16FrequencySeries(tilde_hplus);
		XLALDestroyCOMPLEX16FrequencySeries(tilde_hcross);
		XLALDestroyREAL8FFTPlan(plan);
		XLALDestroyREAL8TimeSeries(*hplus);
		XLALDestroyREAL8TimeSeries(*hcross);
		*hplus = *hcross = NULL;
		XLAL_ERROR(XLAL_EFUNC);
	}
	i = XLALREAL8TimeFreqFFT(tilde_hplus, *hplus, plan);
	i |= XLALREAL8TimeFreqFFT(tilde_hcross, *hcross, plan);
	XLALDestroyREAL8FFTPlan(plan);
	if(i) {
		XLALDestroyCOMPLEX16FrequencySeries(tilde_hplus);
		XLALDestroyCOMPLEX16FrequencySeries(tilde_hcross);
		XLALDestroyREAL8TimeSeries(*hplus);
		XLALDestroyREAL8TimeSeries(*hcross);
		*hplus = *hcross = NULL;
		XLAL_ERROR(XLAL_EFUNC);
	}

	/* apply the frequency-domain Gaussian window.  the window
	 * function's shape parameter is computed similarly to that of the
	 * time-domain window, with \sigma_{f} = \Delta f / 2.  the window
	 * is created with its peak on the middle sample, which we need to
	 * shift to the sample corresponding to the injection's centre
	 * frequency. */

	window = XLALCreateGaussREAL8Window(2 * tilde_hplus->data->length + 1, (tilde_hplus->data->length * tilde_hplus->deltaF) / (bandwidth / 2.0));
	if(!window) {
		XLALDestroyCOMPLEX16FrequencySeries(tilde_hplus);
		XLALDestroyCOMPLEX16FrequencySeries(tilde_hcross);
		XLALDestroyREAL8TimeSeries(*hplus);
		XLALDestroyREAL8TimeSeries(*hcross);
		*hplus = *hcross = NULL;
		XLAL_ERROR(XLAL_EFUNC);
	}
	XLALResizeREAL8Sequence(window->data, tilde_hplus->data->length - (unsigned) floor(frequency / tilde_hplus->deltaF + 0.5), tilde_hplus->data->length);
	for(i = 0; i < window->data->length; i++) {
		tilde_hplus->data->data[i] *= window->data->data[i];
		tilde_hcross->data->data[i] *= window->data->data[i];
	}
	XLALDestroyREAL8Window(window);

	/* normalize the waveform to achieve the desired \int
	 * \f$(\stackrel{.}{h}_{+}^{2} + \stackrel{.}{h}_{\times}^{2}) dt\f$ */

	norm_factor = sqrt((XLALMeasureIntHDotSquaredDT(tilde_hplus) + XLALMeasureIntHDotSquaredDT(tilde_hcross)) / int_hdot_squared);
	if(int_hdot_squared == 0 || norm_factor == 0) {
		XLALDestroyCOMPLEX16FrequencySeries(tilde_hplus);
		XLALDestroyCOMPLEX16FrequencySeries(tilde_hcross);
		XLALDestroyREAL8TimeSeries(*hplus);
		XLALDestroyREAL8TimeSeries(*hcross);
		*hplus = *hcross = NULL;
		XLAL_ERROR(XLAL_EFPDIV0);
	}
	for(i = 0; i < tilde_hplus->data->length; i++) {
		tilde_hplus->data->data[i] /= norm_factor;
		tilde_hcross->data->data[i] /= norm_factor;
	}

	/* transform to the time domain */

	plan = XLALCreateReverseREAL8FFTPlan((*hplus)->data->length, 0);
	if(!plan) {
		XLALDestroyCOMPLEX16FrequencySeries(tilde_hplus);
		XLALDestroyCOMPLEX16FrequencySeries(tilde_hcross);
		XLALDestroyREAL8TimeSeries(*hplus);
		XLALDestroyREAL8TimeSeries(*hcross);
		*hplus = *hcross = NULL;
		XLAL_ERROR(XLAL_EFUNC);
	}
	i = XLALREAL8FreqTimeFFT(*hplus, tilde_hplus, plan);
	i |= XLALREAL8FreqTimeFFT(*hcross, tilde_hcross, plan);
	XLALDestroyREAL8FFTPlan(plan);
	XLALDestroyCOMPLEX16FrequencySeries(tilde_hplus);
	XLALDestroyCOMPLEX16FrequencySeries(tilde_hcross);
	if(i) {
		XLALDestroyREAL8TimeSeries(*hplus);
		XLALDestroyREAL8TimeSeries(*hcross);
		*hplus = *hcross = NULL;
		XLAL_ERROR(XLAL_EFUNC);
	}

	/* force the sample rate incase round-off has shifted it a bit */

	(*hplus)->deltaT = (*hcross)->deltaT = delta_t;

	/* apply a Tukey window for continuity at the start and end of the
	 * injection.  the window's shape parameter sets what fraction of
	 * the window is used by the tapers */

	window = XLALCreateTukeyREAL8Window((*hplus)->data->length, 0.5);
	if(!window) {
		XLALDestroyREAL8TimeSeries(*hplus);
		XLALDestroyREAL8TimeSeries(*hcross);
		*hplus = *hcross = NULL;
		XLAL_ERROR(XLAL_EFUNC);
	}
	for(i = 0; i < window->data->length; i++) {
		(*hplus)->data->data[i] *= window->data->data[i];
		(*hcross)->data->data[i] *= window->data->data[i];
	}
	XLALDestroyREAL8Window(window);

	/* done */

	return 0;
}


/*
 * ============================================================================
 *
 *                         Sine-Gaussian and Friends
 *
 * ============================================================================
 */


/**
 * Input:
 *
 * Q:  the "Q" of the waveform.  The Gaussian envelope is \f$exp(-1/2 t^{2} /
 * \sigma_{t}^{2})\f$ where \f$\sigma_{t} = Q / (2 \pi f)\f$.  High Q --> long
 * duration.
 *
 * centre_frequency:   the frequency of the sinusoidal oscillations that
 * get multiplied by the Gaussian envelope.
 *
 * hrss:  the root-sum-squares strain of the waveform (summed over both
 * polarizations).  See K. Riles, LIGO-T040055-00.pdf.
 *
 * eccentricity:  0 --> circularly polarized, 1 --> linearly polarized.
 *
 * polarization:  the angle from the + axis to the major axis of the
 * waveform ellipsoid.  with the eccentricity set to 1 (output is linearly
 * polarized):  0 --> output contains + polarization only;  pi/2 --> output
 * contains x polarization only.  with the eccentricity set to 0 (output is
 * circularly polarized), the polarization parameter is irrelevant.
 *
 * Output:
 *
 * h+ and hx time series containing a cosine-Gaussian in the + polarization
 * and a sine-Gaussian in the x polarization.  The Gaussian envelope peaks
 * in both at t = 0 as defined by epoch and deltaT.  Note that a Tukey
 * window with tapers covering 50% of the time series is applied to make
 * the waveform go to 0 smoothly at the start and end.
 */


int XLALSimBurstSineGaussian(
	REAL8TimeSeries **hplus,
	REAL8TimeSeries **hcross,
	REAL8 Q,
	REAL8 centre_frequency,
	REAL8 hrss,
	REAL8 eccentricity,
	REAL8 polarization,
	REAL8 delta_t
)
{
	REAL8Window *window;
	/* semimajor and semiminor axes of waveform ellipsoid */
	const double a = 1.0 / sqrt(2.0 - eccentricity * eccentricity);
	const double b = a * sqrt(1.0 - eccentricity * eccentricity);
	/* rss of plus and cross polarizations */
	const double hplusrss  = hrss * (a * cos(polarization) - b * sin(polarization));
	const double hcrossrss = hrss * (b * cos(polarization) + a * sin(polarization));
	/* rss of unit amplitude cosine- and sine-gaussian waveforms.  see
	 * K. Riles, LIGO-T040055-00.pdf */
	const double cgrss = sqrt((Q / (4.0 * centre_frequency * sqrt(LAL_PI))) * (1.0 + exp(-Q * Q)));
	const double sgrss = sqrt((Q / (4.0 * centre_frequency * sqrt(LAL_PI))) * (1.0 - exp(-Q * Q)));
	/* "peak" amplitudes of plus and cross */
	const double h0plus  = hplusrss / cgrss;
	const double h0cross = hcrossrss / sgrss;
	LIGOTimeGPS epoch;
	int length;
	unsigned i;

	/* length of the injection time series is 30 * the width of the
	 * Gaussian envelope (sigma_t in the comments above), rounded to
	 * the nearest odd integer */

	length = (int) floor(30.0 * Q / (LAL_TWOPI * centre_frequency) / delta_t / 2.0);
	length = 2 * length + 1;

	/* the middle sample is t = 0 */

	XLALGPSSetREAL8(&epoch, -(length - 1) / 2 * delta_t);

	/* allocate the time series */

	*hplus = XLALCreateREAL8TimeSeries("sine-Gaussian +", &epoch, 0.0, delta_t, &lalStrainUnit, length);
	*hcross = XLALCreateREAL8TimeSeries("sine-Gaussian x", &epoch, 0.0, delta_t, &lalStrainUnit, length);
	if(!*hplus || !*hcross) {
		XLALDestroyREAL8TimeSeries(*hplus);
		XLALDestroyREAL8TimeSeries(*hcross);
		*hplus = *hcross = NULL;
		XLAL_ERROR(XLAL_EFUNC);
	}

	/* populate */

	for(i = 0; i < (*hplus)->data->length; i++) {
		const double t = ((int) i - (length - 1) / 2) * delta_t;
		const double phi = LAL_TWOPI * centre_frequency * t;
		const double fac = exp(-0.5 * phi * phi / (Q * Q));
		(*hplus)->data->data[i]  = h0plus * fac * cos(phi);
		(*hcross)->data->data[i] = h0cross * fac * sin(phi);
	}

	/* apply a Tukey window for continuity at the start and end of the
	 * injection.  the window's shape parameter sets what fraction of
	 * the window is used by the tapers */

	window = XLALCreateTukeyREAL8Window((*hplus)->data->length, 0.5);
	if(!window) {
		XLALDestroyREAL8TimeSeries(*hplus);
		XLALDestroyREAL8TimeSeries(*hcross);
		*hplus = *hcross = NULL;
		XLAL_ERROR(XLAL_EFUNC);
	}
	for(i = 0; i < window->data->length; i++) {
		(*hplus)->data->data[i] *= window->data->data[i];
		(*hcross)->data->data[i] *= window->data->data[i];
	}
	XLALDestroyREAL8Window(window);

	return 0;
}


/*
 * ============================================================================
 *
 *                                String Cusp
 *
 * ============================================================================
 */


/**
 * Input:
 *	amplitude = waveform's amplitude parameter
 *	f_high = high frequency cutoff
 *	delta_t = sample period of output time series
 *
 * Output:
 * 	h+(t) and hx(t), where the cusp waveform has been placed entirely
 * 	in the + polarization (the x polarization is zeroed), and the
 * 	waveform peaks at t = 0 (as defined by the epoch and deltaT).
 *
 * The low frequency cut-off is fixed at 1 Hz;  there's nothing special
 * about 1 Hz except that it is low compared to the frequency at which we
 * should be high-passing the data
 */


int XLALGenerateStringCusp(
	REAL8TimeSeries **hplus,
	REAL8TimeSeries **hcross,
	REAL8 amplitude,
	REAL8 f_high,
	REAL8 delta_t
)
{
	COMPLEX16FrequencySeries *tilde_h;
	REAL8FFTPlan *plan;
	LIGOTimeGPS epoch;
	int length;
	int i;
	/* low frequency cut-off in Hertz */
	const double f_low = 1.0;

	/* check input */

	if(amplitude < 0 || f_high < f_low || delta_t <= 0) {
		XLALPrintError("%s(): invalid input parameters\n", __func__);
		*hplus = *hcross = NULL;
		XLAL_ERROR(XLAL_EINVAL);
	}

	/* length of the injection time series is 15 / f_low, rounded to
	 * the nearest odd integer */

	length = (int) (15 / f_low / delta_t / 2.0);
	length = 2 * length + 1;

	/* the middle sample is t = 0 */

	XLALGPSSetREAL8(&epoch, -(length - 1) / 2 * delta_t);

	/* allocate time and frequency series and FFT plan */

	*hplus = XLALCreateREAL8TimeSeries("string cusp +", &epoch, 0.0, delta_t, &lalStrainUnit, length);
	*hcross = XLALCreateREAL8TimeSeries("string cusp x", &epoch, 0.0, delta_t, &lalStrainUnit, length);
	tilde_h = XLALCreateCOMPLEX16FrequencySeries("string cusp +", &epoch, 0.0, 1.0 / (length * delta_t), &lalDimensionlessUnit, length / 2 + 1);
	plan = XLALCreateReverseREAL8FFTPlan(length, 0);
	if(!*hplus || !*hcross || !tilde_h || !plan) {
		XLALDestroyREAL8TimeSeries(*hplus);
		XLALDestroyREAL8TimeSeries(*hcross);
		XLALDestroyCOMPLEX16FrequencySeries(tilde_h);
		XLALDestroyREAL8FFTPlan(plan);
		*hplus = *hcross = NULL;
		XLAL_ERROR(XLAL_EFUNC);
	}
	XLALUnitMultiply(&tilde_h->sampleUnits, &(*hplus)->sampleUnits, &lalSecondUnit);

	/* zero the x time series, injection is done in + only */

	memset((*hcross)->data->data, 0, (*hcross)->data->length * sizeof(*(*hcross)->data->data));

	/* construct the waveform in the frequency domain */

	for(i = 0; (unsigned) i < tilde_h->data->length; i++) {
		double f = tilde_h->f0 + i * tilde_h->deltaF;

		/* frequency-domain wave form.  includes taper factor above
		 * h_high, and phase shift to put waveform's peak on the
		 * middle sample of the time series */

		double amp = amplitude * pow((sqrt(1 + f_low * f_low / (f * f))), -8) * pow(f, -4.0 / 3.0) * (f > f_high ? exp(1 - f / f_high) : 1);

		tilde_h->data->data[i] = amp * cexp(-I * LAL_PI * i * (length - 1) / length);
	}

	/* set DC and Nyquist to zero */

	tilde_h->data->data[0] = tilde_h->data->data[tilde_h->data->length - 1] = 0;

	/* transform to time domain */

	i = XLALREAL8FreqTimeFFT(*hplus, tilde_h, plan);
	XLALDestroyCOMPLEX16FrequencySeries(tilde_h);
	XLALDestroyREAL8FFTPlan(plan);
	if(i) {
		XLALDestroyREAL8TimeSeries(*hplus);
		XLALDestroyREAL8TimeSeries(*hcross);
		*hplus = *hcross = NULL;
		XLAL_ERROR(XLAL_EFUNC);
	}

	/* force the sample rate incase round-off has shifted it a bit */

	(*hplus)->deltaT = (*hcross)->deltaT = delta_t;

	/* apodize the time series */
	/* FIXME:  use a Tukey window? */

	for(i = (*hplus)->data->length - 1; i >= 0; i--)
		(*hplus)->data->data[i] -= (*hplus)->data->data[0];

	/* done */

	return 0;
}
