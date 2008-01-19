/*
 * Copyright (C) 2007 Jolien Creighton, Patrick Brady, Saikat Ray-Majumder,
 * Xavier Siemens, Teviet Creighton, Kipp Cannon
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
#include <lal/FrequencySeries.h>
#include <lal/LALSimulation.h>
#include <lal/Sequence.h>
#include <lal/TimeSeries.h>


/* FIXME:  which of these are still needed? */
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Units.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/VectorOps.h>
#include <lal/Inject.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/TimeFreqFFT.h>
#include <lal/GenerateBurst.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/RealFFT.h>


NRCSID(GENERATEBURSTC, "$Id$");


/*
 * ============================================================================
 *
 *           Fill a time series with stationary white Gaussin noise
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
 *                               Normalizations
 *
 * ============================================================================
 */


/*
 * Returns the strain of the sample with the largest magnitude.
 */


static REAL8 XLALMeasureHPeak(REAL8TimeSeries *series)
{
	static const char func[] = "XLALMeasureHPeak";
	double hpeak;
	unsigned i;

	if(!series->data->length)
		XLAL_ERROR_REAL8(func, XLAL_EBADLEN);

	hpeak = series->data->data[0];
	for(i = 1; i < series->data->length; i++)
		if(fabs(series->data->data[i]) > fabs(hpeak))
			hpeak = series->data->data[i];

	return hpeak;
}


/*
 * Returns the root-sum-square strain.
 */


static REAL8 XLALMeasureHrss(REAL8TimeSeries *series)
{
	double e = 0.0;
	double sum = 0.0;
	unsigned i;

	/* Kahans's compensated summation algorithm */

	for(i = 0; i < series->data->length; i++) {
		double tmp = sum;
		/* what we want to add = h^{2} + "error from last
		 * iteration" */
		double x = series->data->data[i] * series->data->data[i] + e;
		/* add */
		sum += x;
		/* negative of what was actually added */
		e = tmp - sum;
		/* what didn't get added, add next time */
		e += x;
	}

	return sqrt(sum);
}


/*
 * Given the Fourier transform of a real-valued function h(t), compute and
 * return the integral of the square of its derivative:
 *
 * \int \dot{h}^{2} \diff t.
 *
 * The normalization factors in this function assume that
 * XLALREAL8FreqTimeFFT() will be used to convert the frequency series to
 * the time domain.
 */


static REAL8 XLALMeasureIntHDotSquaredDT(COMPLEX16FrequencySeries *fseries)
{
	unsigned i;
	double e = 0.0;
	double sum = 0.0;

	/* Kahan's compensated summation algorithm. The summation is done
	 * from lowest to highest frequency under the assumption that high
	 * frequency components tend to add more to the magnitude of the
	 * derivative.  */

	for(i = 0; i < fseries->data->length; i++) {
		double tmp = sum;
		/* what we want to add = f^{2} |\tilde{s}(f)|^{2} + "error
		 * from last iteration" */
		double x = (fseries->f0 + i * fseries->deltaF) * (fseries->f0 + i * fseries->deltaF) * (fseries->data->data[i].re * fseries->data->data[i].re + fseries->data->data[i].im * fseries->data->data[i].im) + e;
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


/*
 * ============================================================================
 *
 *            Construct a Band- and Time-Limited White Noise Burst
 *
 * ============================================================================
 */


/*
 * Parameters:
 *
 * duration
 * 	width of time domain Gaussian envelope in seconds
 * frequency
 * 	centre frequency of waveform in Hertz
 * bandwidth
 * 	width of frequency domain Gaussian envelope in Hertz
 * int_hdot_squared
 * 	waveform is normalized so that \int (\dot{h}_{+}^{2} +
 * 	\dot{h}_{\times}^{2}) \diff t equals this
 * delta_t
 * 	the sample rate of the time series to construct
 * rng
 * 	a GSL random number generator to be used to produce Gaussian random
 * 	variables
 *
 * Output:
 *
 * Two time series containing h+(t) and hx(t), with the injection centred
 * on t = 0 (as defined by the epoch and deltaT).  The + and x time series
 * are two independent injections.
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


int XLALBandAndTimeLimitedWhiteNoiseBurst(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 duration, REAL8 frequency, REAL8 bandwidth, REAL8 int_hdot_squared, REAL8 delta_t, gsl_rng *rng)
{
	static const char func[] = "XLALBandAndTimeLimitedWhiteNoiseBurst";
	int length;
	LIGOTimeGPS epoch;
	COMPLEX16FrequencySeries *tilde_hplus, *tilde_hcross;
	REAL8Window *window;
	REAL8FFTPlan *plan;
	REAL8 norm_factor;
	unsigned i;

	/* check input */

	if(duration < 0 || bandwidth < 0 || duration * bandwidth < LAL_2_PI || int_hdot_squared < 0 || delta_t <= 0) {
		*hplus = *hcross = NULL;
		XLAL_ERROR(func, XLAL_EINVAL);
	}

	/* length of the injection time series is 10 * duration, rounded to
	 * the nearest odd integer */

	length = (int) (10.0 * duration / delta_t / 2.0);
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
		XLAL_ERROR(func, XLAL_EFUNC);
	}

	/* fill with independent zero-mean unit variance Gaussian random
	 * numbers (any non-zero amplitude is OK, it will be adjusted
	 * later) */

	gaussian_noise(*hplus, 1, rng);
	gaussian_noise(*hcross, 1, rng);

	/* apply the time-domain Gaussian window.  the window function's
	 * shape parameter is ((length - 1) * delta_t / 2) / \sigma_{t} where
	 *
	 * \sigma_{t} = \sqrt{duration^{2} / 4 - 1 / (\pi^{2} bandwidth^{2})}
	 *
	 * is the compensated time-domain window duration */

	window = XLALCreateGaussREAL8Window((*hplus)->data->length, (((*hplus)->data->length - 1) * delta_t / 2) / sqrt(duration * duration / 4.0 - 1.0 / (LAL_PI * LAL_PI * bandwidth * bandwidth)));
	if(!window) {
		XLALDestroyREAL8TimeSeries(*hplus);
		XLALDestroyREAL8TimeSeries(*hcross);
		*hplus = *hcross = NULL;
		XLAL_ERROR(func, XLAL_EFUNC);
	}
	for(i = 0; i < window->data->length; i++) {
		(*hplus)->data->data[i] *= window->data->data[i];
		(*hcross)->data->data[i] *= window->data->data[i];
	}
	XLALDestroyREAL8Window(window);

	/* transform to the frequency domain */

	plan = XLALCreateForwardREAL8FFTPlan((*hplus)->data->length, 0);
	tilde_hplus = XLALCreateCOMPLEX16FrequencySeries(NULL, &epoch, 0.0, 0.0, &lalDimensionlessUnit, (*hplus)->data->length / 2 + 1);
	tilde_hcross = XLALCreateCOMPLEX16FrequencySeries(NULL, &epoch, 0.0, 0.0, &lalDimensionlessUnit, (*hplus)->data->length / 2 + 1);
	if(!plan || !tilde_hplus || !tilde_hcross) {
		XLALDestroyCOMPLEX16FrequencySeries(tilde_hplus);
		XLALDestroyCOMPLEX16FrequencySeries(tilde_hcross);
		XLALDestroyREAL8FFTPlan(plan);
		XLALDestroyREAL8TimeSeries(*hplus);
		XLALDestroyREAL8TimeSeries(*hcross);
		*hplus = *hcross = NULL;
		XLAL_ERROR(func, XLAL_EFUNC);
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
		XLAL_ERROR(func, XLAL_EFUNC);
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
		XLAL_ERROR(func, XLAL_EFUNC);
	}
	XLALResizeREAL8Sequence(window->data, tilde_hplus->data->length - (unsigned) floor(frequency / tilde_hplus->deltaF + 0.5), tilde_hplus->data->length);
	for(i = 0; i < window->data->length; i++) {
		tilde_hplus->data->data[i].re *= window->data->data[i];
		tilde_hplus->data->data[i].im *= window->data->data[i];
		tilde_hcross->data->data[i].re *= window->data->data[i];
		tilde_hcross->data->data[i].im *= window->data->data[i];
	}
	XLALDestroyREAL8Window(window);

	/* normalize the waveform to achieve the desired \int
	 * (\dot{h}_{+}^{2} + \dot{h}_{\times}^{2}) dt */

	norm_factor = sqrt(int_hdot_squared / (XLALMeasureIntHDotSquaredDT(tilde_hplus) + XLALMeasureIntHDotSquaredDT(tilde_hcross)));
	for(i = 0; i < tilde_hplus->data->length; i++) {
		tilde_hplus->data->data[i].re *= norm_factor;
		tilde_hplus->data->data[i].im *= norm_factor;
		tilde_hcross->data->data[i].re *= norm_factor;
		tilde_hcross->data->data[i].im *= norm_factor;
	}

	/* transform to the time domain */

	plan = XLALCreateReverseREAL8FFTPlan((*hplus)->data->length, 0);
	if(!plan) {
		XLALDestroyCOMPLEX16FrequencySeries(tilde_hplus);
		XLALDestroyCOMPLEX16FrequencySeries(tilde_hcross);
		XLALDestroyREAL8TimeSeries(*hplus);
		XLALDestroyREAL8TimeSeries(*hcross);
		*hplus = *hcross = NULL;
		XLAL_ERROR(func, XLAL_EFUNC);
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
		XLAL_ERROR(func, XLAL_EFUNC);
	}

	/* apply a Tukey window for continuity at the start and end of the
	 * injection.  the window's shape parameter sets what fraction of
	 * the window is used by the tapers */

	window = XLALCreateTukeyREAL8Window((*hplus)->data->length, 0.5);
	if(!window) {
		XLALDestroyREAL8TimeSeries(*hplus);
		XLALDestroyREAL8TimeSeries(*hcross);
		*hplus = *hcross = NULL;
		XLAL_ERROR(func, XLAL_EFUNC);
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


/*
 * ============================================================================
 *
 *                                String Cusp
 *
 * ============================================================================
 */


/*
 * Input:
 *	amplitude = waveform's amplitude parameter
 *	f_high = high frequency cutoff
 *	delta_t = sample period of output time series
 *
 * Output:
 * 	h(t) with waveform peak at t = 0 (as defined by the epoch and
 * 	deltaT).
 *
 * The low frequency cutoff is fixed at 1 Hz;  there's nothing special
 * about 1 Hz except that it is low compared to the frequency at which we
 * should be high-passing the data
 */


int XLALGenerateStringCusp(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 amplitude, REAL8 f_high, REAL8 delta_t)
{
	static const char func[] = "XLALGenerateStringCusp";
	COMPLEX16FrequencySeries *tilde_h;
	REAL8FFTPlan *plan;
	LIGOTimeGPS epoch;
	int length;
	int i;
	const double f_low = 1.0;

	/* check input */

	if(amplitude < 0 || f_high < f_low || delta_t <= 0) {
		*hplus = *hcross = NULL;
		XLAL_ERROR(func, XLAL_EINVAL);
	}

	/* length of the injection time series is 5 / f_low, rounded to the
	 * nearest odd integer */

	length = (int) (5 / f_low / delta_t / 2.0);
	length = 2 * length + 1;

	/* the middle sample is t = 0 */

	XLALGPSSetREAL8(&epoch, -(length - 1) / 2 * delta_t);

	/* allocate time and frequency series and FFT plan */

	*hplus = XLALCreateREAL8TimeSeries("string cusp +", &epoch, 0.0, delta_t, &lalStrainUnit, length);
	*hcross = XLALCreateREAL8TimeSeries("string cusp x", &epoch, 0.0, delta_t, &lalStrainUnit, length);
	tilde_h = XLALCreateCOMPLEX16FrequencySeries("string cusp", &epoch, 0.0, 1.0 / (length * delta_t), &lalDimensionlessUnit, length / 2 + 1);
	plan = XLALCreateReverseREAL8FFTPlan(length, 0);
	if(!*hplus || !*hcross || !tilde_h || !plan) {
		XLALDestroyREAL8TimeSeries(*hplus);
		XLALDestroyREAL8TimeSeries(*hcross);
		XLALDestroyCOMPLEX16FrequencySeries(tilde_h);
		XLALDestroyREAL8FFTPlan(plan);
		*hplus = *hcross = NULL;
		XLAL_ERROR(func, XLAL_EFUNC);
	}
	XLALUnitMultiply(&tilde_h->sampleUnits, &(*hplus)->sampleUnits, &lalSecondUnit);

	/* zero the cross time series, injection is done in + only */

	memset((*hcross)->data->data, 0, (*hcross)->data->length * sizeof(*(*hcross)->data->data));

	/* construct the waveform in the frequency domain */

	for(i = 0; (unsigned) i < tilde_h->data->length; i++) {
		double f = i * tilde_h->deltaF;

		/* frequency-domain wave form */

		tilde_h->data->data[i].re = amplitude * pow((sqrt(1 + f_low * f_low / (f * f))), -8) * pow(f, -4.0 / 3.0);
		if(f > f_high)
			tilde_h->data->data[i].re *= exp(1 - f / f_high);
		tilde_h->data->data[i].im = tilde_h->data->data[i].re;

		/* phase shift to put waveform's peak on the middle sample
		 * of the time series */

		tilde_h->data->data[i].im *= sin(-LAL_PI * i * (length - 1) / length);
		tilde_h->data->data[i].re *= cos(-LAL_PI * i * (length - 1) / length);
	}

	/* set DC to zero */

	tilde_h->data->data[0].re = 0;
	tilde_h->data->data[0].im = 0;

	/* set Nyquist to zero */

	tilde_h->data->data[tilde_h->data->length - 1].re = 0;
	tilde_h->data->data[tilde_h->data->length - 1].im = 0;

	/* transform to time domain */

	i = XLALREAL8FreqTimeFFT(*hplus, tilde_h, plan);
	XLALDestroyCOMPLEX16FrequencySeries(tilde_h);
	XLALDestroyREAL8FFTPlan(plan);
	if(i) {
		XLALDestroyREAL8TimeSeries(*hplus);
		XLALDestroyREAL8TimeSeries(*hcross);
		*hplus = *hcross = NULL;
		XLAL_ERROR(func, XLAL_EFUNC);
	}

	/* apodize the time series */

	for(i = (*hplus)->data->length - 1; i >= 0; i--)
		(*hplus)->data->data[i] -= (*hplus)->data->data[0];

	/* done */

	return 0;
}


/*
 * ============================================================================
 *
 *                              sim_burst Nexus
 *
 * ============================================================================
 */


/*
 * Convenience wrapper to iterate over the entries in a sim_burst linked
 * list and inject them into a strain time series.
 */


static int XLALBurstInjectSignals(LALDetector *detector, REAL8TimeSeries *h, SimBurstTable *sim_burst)
{
	static const char func[] = "XLALBurstInjectSignals";
	/* + and x time series for injection waveform */
	REAL8TimeSeries *injection_hplus, *injection_hcross;
	/* injection time series as added to detector's */
	REAL8TimeSeries *injection_h;

	for(; sim_burst; sim_burst = sim_burst->next) {
		/* construct the h+ and hx time series for the injection
		 * waveform */

		if(!strcmp(sim_burst->waveform, "BTLWNB")) {
			/* hrss --> int \dot{h}^2 dt, freq --> f_{0},
			 * dtplus --> duration, dtminus --> bandwidth,
			 * zm_number --> seed */
			gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
			if(!rng)
				XLAL_ERROR(func, XLAL_ENOMEM);
			gsl_rng_set(rng, sim_burst->zm_number);
			if(XLALBandAndTimeLimitedWhiteNoiseBurst(&injection_hplus, &injection_hcross, sim_burst->dtplus, sim_burst->freq, sim_burst->dtminus, sim_burst->hrss, h->deltaT, rng)) {
				gsl_rng_free(rng);
				XLAL_ERROR(func, XLAL_EFUNC);
			}
			gsl_rng_free(rng);
		} else if(!strcmp(sim_burst->waveform, "StringCusp")) {
			/* hpeak --> amplitude, freq --> f_{high} */
			if(XLALGenerateStringCusp(&injection_hplus, &injection_hcross, sim_burst->hpeak, sim_burst->freq, h->deltaT))
				XLAL_ERROR(func, XLAL_EFUNC);
		} else
			/* unrecognized waveform */
			XLAL_ERROR(func, XLAL_EINVAL);

		/* project the wave strain onto the detector's response
		 * tensor to produce the injection strain as seen in the
		 * detector.  longitude --> right ascension, latitude -->
		 * declination, polarization --> psi, geocent_peak_time -->
		 * "time" of injection at geocentre */

		injection_h = XLALSimDetectorStrainREAL8TimeSeries(injection_hplus, injection_hcross, sim_burst->longitude, sim_burst->latitude, sim_burst->polarization, detector, &sim_burst->geocent_peak_time);
		XLALDestroyREAL8TimeSeries(injection_hplus);
		XLALDestroyREAL8TimeSeries(injection_hcross);
		if(!injection_h)
			XLAL_ERROR(func, XLAL_EFUNC);

		/* add the injection strain time series to the detector
		 * data */

		if(XLALAddInjectionREAL8TimeSeries(h, injection_h, NULL)) {
			XLALDestroyREAL8TimeSeries(injection_h);
			XLAL_ERROR(func, XLAL_EFUNC);
		}
		XLALDestroyREAL8TimeSeries(injection_h);
	}

	return 0;
}


/*
 * ============================================================================
 *
 *                   Legacy Code --- Please Update to XLAL!
 *
 * ============================================================================
 */


static int XLALGenerateBurst(CoherentGW *output, SimBurstTable *simBurst, BurstParamStruc *params)
{
	static const char func[] = "XLALGenerateBurst";
	UINT4 n, i;		/* number of and index over samples */
	REAL8 t, dt, duration;	/* time, interval */
	REAL8 t0, tau, gtime;	/* central time, decay time, gaussian time */
	REAL8 f0;		/* initial frequency */
	REAL8 twopif0;		/* 2*pi*f0 */
	REAL4 hpeak;		/* peak strain for burst */
	REAL4 *fData;		/* pointer to frequency data */
	REAL8 *phiData;		/* pointer to phase data */
	REAL4 *aData;		/* pointer to frequency data */
	LIGOTimeGPS startTime;	/* start time of injection */

	/* Set up some constants to avoid repeated dereferencing.  notice
	 * the factor of 2 in the definition of n confusingly makes
	 * injections twice as long as the variable duration */
	duration = simBurst->dtplus + simBurst->dtminus;
	dt = params->deltaT;
	n = 2.0 * duration / dt;
	if(!n)
		XLAL_ERROR(func, XLAL_EINVAL);

	/* start time of data is peak time duration */
	startTime = simBurst->geocent_peak_time;
	XLALGPSAdd(&startTime, -duration);

	/* Generic burst parameters */
	hpeak = simBurst->hpeak;
	tau = simBurst->tau;
	f0 = simBurst->freq;
	twopif0 = LAL_TWOPI * f0;

	/* Allocate output structures. */
	output->a = LALCalloc(1, sizeof(*output->a));
	output->f = XLALCreateREAL4TimeSeries("Burst frequency", &startTime, 0.0, params->deltaT, &lalHertzUnit, n);
	output->phi = XLALCreateREAL8TimeSeries("Burst phase", &startTime, 0.0, params->deltaT, &lalDimensionlessUnit, n);
	if(!output->a || !output->f || !output->phi) {
		LALFree(output->a);
		XLALDestroyREAL4TimeSeries(output->f);
		XLALDestroyREAL8TimeSeries(output->phi);
		output->a = NULL;
		output->f = NULL;
		output->phi = NULL;
		XLAL_ERROR(func, XLAL_EFUNC);
	}
	output->a->data = XLALCreateREAL4VectorSequence(n, 2);
	if(!output->a->data) {
		LALFree(output->a);
		XLALDestroyREAL4TimeSeries(output->f);
		XLALDestroyREAL8TimeSeries(output->phi);
		output->a = NULL;
		output->f = NULL;
		output->phi = NULL;
		XLAL_ERROR(func, XLAL_EFUNC);
	}

	/* Set output structure metadata fields. */
	output->position.longitude = simBurst->longitude;
	output->position.latitude = simBurst->latitude;
	output->position.system = params->system;
	output->psi = simBurst->polarization;
	output->a->epoch = startTime;
	output->a->deltaT = params->deltaT;
	output->a->sampleUnits = lalStrainUnit;
	LALSnprintf(output->a->name, LALNameLength, "Burst amplitudes");

	/* Fill frequency and phase arrays. */
	fData = output->f->data->data;
	phiData = output->phi->data->data;
	aData = output->a->data->data;

	/* this depends on the waveform type */
	if(!(strcmp(simBurst->waveform, "SineGaussian"))) {
		/* find the peak time as a REAL8 relative to start of segment */
		t0 = XLALGPSDiff(&(simBurst->geocent_peak_time), &startTime);

		/* construct the signal */
		for(i = 0; i < n; i++) {
			t = i * dt;
			gtime = (t - t0) / tau;
			*(fData++) = f0;
			*(phiData++) = twopif0 * (t - t0);
			*(aData++) = hpeak * exp(-gtime * gtime);
			*(aData++) = 0.0;
		}
	} else if(!(strcmp(simBurst->waveform, "Gaussian"))) {
		/* find the peak time as a REAL8 relative to start of segment */
		t0 = XLALGPSDiff(&(simBurst->geocent_peak_time), &startTime);

		/* construct the signal */
		for(i = 0; i < n; i++) {
			t = i * dt;
			gtime = (t - t0) / tau;
			*(fData++) = 0.0;
			*(phiData++) = 0.0;
			*(aData++) = hpeak * exp(-gtime * gtime);
			*(aData++) = 0.0;
		}
	} else if(!(strcmp(simBurst->waveform, "Ringdown"))) {
		t0 = XLALGPSDiff(&(simBurst->geocent_peak_time), &startTime);
		for(i = 0; i < n; i++) {
			t = i * dt;
			gtime = (t - t0) / tau;
			*fData++ = f0;
			*phiData++ = twopif0 * (t - t0);
			if(gtime > 0)
				*aData++ = hpeak * exp(-gtime);
			else
				*aData++ = 0;
			*aData++ = 0;
		}
	} else if(!(strcmp(simBurst->waveform, "Ringup"))) {
		t0 = XLALGPSDiff(&(simBurst->geocent_peak_time), &startTime);
		for(i = 0; i < n; i++) {
			t = i * dt;
			gtime = (t - t0) / tau;
			*fData++ = f0;
			*phiData++ = twopif0 * (t - t0);
			if(gtime < 0)
				*aData++ = hpeak * exp(gtime);
			else
				*aData++ = 0;
			*aData++ = 0;
		}
	} else if(!(strcmp(simBurst->waveform, "StringCusp"))) {
		/* I use the simburst table as follows: The duration is still
		   dtplus+dtminus; hpeak is the amplitude of the cusp, not the
		   value of the strain at the peak, t0 is the central time of the
		   cusp, which introduces a phase into the waveform. The low
		   frequency cutoff will be fixed at 1Hz there's nothing special
		   about 1Hz except that it low compared to the ferquecny at which
		   we should be high-passing the data; the high frequency cutoff
		   is given by f0 */
		REAL4Vector *vector;
		COMPLEX8Vector *vtilde;
		RealFFTPlan *rplan;
		REAL4 dfreq = 1 / (2 * duration);	/* the factor of two here is becaus the length of injections 
							   is actually twice the value of the variable duration */
		REAL4 flow = 1;

		t0 = XLALGPSDiff(&(simBurst->geocent_peak_time), &startTime);

		/* create vector to store h(t), create vector that will
		 * hold frequency domain template, create fft plan */
		vector = XLALCreateREAL4Sequence(n);
		vtilde = XLALCreateCOMPLEX8Sequence(n / 2 + 1);
		rplan = XLALCreateReverseREAL4FFTPlan(n, 0);
		if(!vector || !vtilde || !rplan) {
			XLALDestroyREAL4Sequence(vector);
			XLALDestroyCOMPLEX8Sequence(vtilde);
			XLALDestroyREAL4FFTPlan(rplan);
			XLAL_ERROR(func, XLAL_EFUNC);
		}

		/* Set the FD template */
		for(i = 0; i < vtilde->length - 1; i++) {
			REAL4 freq = i * dfreq;
			vtilde->data[i].re = hpeak * pow((sqrt(1 + pow(flow, 2) * pow(freq, -2))), -8) * pow(freq, -4.0 / 3.0);

			if(freq >= f0)
				vtilde->data[i].re *= exp(1 - freq / f0);

			vtilde->data[i].im = vtilde->data[i].re * sin(-LAL_TWOPI * freq * duration);
			vtilde->data[i].re = vtilde->data[i].re * cos(-LAL_TWOPI * freq * duration);
		}

		/* set dc to zero */
		vtilde->data[0].re = 0;
		vtilde->data[0].im = 0;
		/* set nyquist to zero */
		vtilde->data[vtilde->length - 1].re = 0;
		vtilde->data[vtilde->length - 1].im = 0;

		/* Reverse FFT */
		if(XLALREAL4ReverseFFT(vector, vtilde, rplan))
			XLAL_ERROR(func, XLAL_EFUNC);

		/* multiply times dfreq to make sure units are correct */
		for(i = 0; i < vector->length; i++)
			vector->data[i] *= dfreq;

		/* make sure injection starts precisely at 0 */
		for(i = 0; i < vector->length; i++)
			vector->data[i] -= vector->data[0];

		for(i = 0; i < n; i++) {
			*fData++ = 0.0;
			*phiData++ = 0.0;
			*aData++ = vector->data[i];
			*aData++ = 0;
		}

		/* free the data */
		XLALDestroyREAL4Sequence(vector);
		XLALDestroyCOMPLEX8Sequence(vtilde);
		XLALDestroyREAL4FFTPlan(rplan);
	} else if(!(strcmp(simBurst->waveform, "warren"))) {
		/* set everything to 0 */
		for(i = 0; i < n; i++) {
			*(fData++) = 0.0;
			*(phiData++) = 0.0;
			*(aData++) = 0.0;
			*(aData++) = 0.0;
		}
	} else {
		/* unknown waveform */
		XLAL_ERROR(func, XLAL_EINVAL);
	}

	return 0;
}


void LALBurstInjectSignals(LALStatus * stat, REAL4TimeSeries * series, SimBurstTable * injections, COMPLEX8FrequencySeries * resp, INT4 calType)
{
	UINT4 k;
	INT4 injStartTime;
	INT4 injStopTime;
	DetectorResponse detector;
	COMPLEX8Vector *unity = NULL;
	CoherentGW waveform;
	BurstParamStruc burstParam;
	REAL4TimeSeries signal;
	SimBurstTable *simBurst = NULL;
	LALDetector *tmpDetector = NULL /*,*nullDetector=NULL */ ;
	COMPLEX8FrequencySeries *transfer = NULL;

	INITSTATUS(stat, "LALBurstInjectSignals", GENERATEBURSTC);
	ATTATCHSTATUSPTR(stat);

	/* set up start and end of injection zone TODO: fix this hardwired 10 */
	injStartTime = series->epoch.gpsSeconds - 10;
	injStopTime = series->epoch.gpsSeconds + 10 + (INT4) (series->data->length * series->deltaT);

	/* 
	 *compute the transfer function 
	 */

	/* allocate memory and copy the parameters describing the freq series */
	memset(&detector, 0, sizeof(DetectorResponse));
	transfer = (COMPLEX8FrequencySeries *)
	    LALCalloc(1, sizeof(COMPLEX8FrequencySeries));
	if(!transfer) {
		ABORT(stat, GENERATEBURSTH_EMEM, GENERATEBURSTH_MSGEMEM);
	}
	memcpy(&(transfer->epoch), &(resp->epoch), sizeof(LIGOTimeGPS));
	transfer->f0 = resp->f0;
	transfer->deltaF = resp->deltaF;

	tmpDetector = detector.site = (LALDetector *) LALMalloc(sizeof(LALDetector));
	/* set the detector site */
	switch (series->name[0]) {
	case 'H':
		*(detector.site) = lalCachedDetectors[LALDetectorIndexLHODIFF];
		LALWarning(stat, "computing waveform for Hanford.");
		break;
	case 'L':
		*(detector.site) = lalCachedDetectors[LALDetectorIndexLLODIFF];
		LALWarning(stat, "computing waveform for Livingston.");
		break;
	default:
		LALFree(detector.site);
		detector.site = NULL;
		tmpDetector = NULL;
		LALWarning(stat, "Unknown detector site, computing plus mode " "waveform with no time delay");
		break;
	}

	/* set up units for the transfer function */
	{
		RAT4 negOne = { -1, 0 };
		LALUnit unit;
		LALUnitPair pair;
		pair.unitOne = &lalADCCountUnit;
		pair.unitTwo = &lalStrainUnit;
		LALUnitRaise(stat->statusPtr, &unit, pair.unitTwo, &negOne);
		CHECKSTATUSPTR(stat);
		pair.unitTwo = &unit;
		LALUnitMultiply(stat->statusPtr, &(transfer->sampleUnits), &pair);
		CHECKSTATUSPTR(stat);
	}

	/* invert the response function to get the transfer function */
	LALCCreateVector(stat->statusPtr, &(transfer->data), resp->data->length);
	CHECKSTATUSPTR(stat);

	LALCCreateVector(stat->statusPtr, &unity, resp->data->length);
	CHECKSTATUSPTR(stat);
	for(k = 0; k < resp->data->length; ++k) {
		unity->data[k].re = 1.0;
		unity->data[k].im = 0.0;
	}

	LALCCVectorDivide(stat->statusPtr, transfer->data, unity, resp->data);
	CHECKSTATUSPTR(stat);

	LALCDestroyVector(stat->statusPtr, &unity);
	CHECKSTATUSPTR(stat);

	/* Set up a time series to hold signal in ADC counts */
	signal.deltaT = series->deltaT;
	if((signal.f0 = series->f0) != 0) {
		ABORT(stat, GENERATEBURSTH_EMEM, GENERATEBURSTH_MSGEMEM);
	}
	signal.sampleUnits = lalADCCountUnit;

	signal.data = NULL;
	LALSCreateVector(stat->statusPtr, &(signal.data), series->data->length);
	CHECKSTATUSPTR(stat);

	/* loop over list of waveforms and inject into data stream */
	for(simBurst = injections; simBurst; simBurst = simBurst->next) {
		/* only do the work if the burst is in injection zone */
		if((injStartTime - simBurst->geocent_peak_time.gpsSeconds) * (injStopTime - simBurst->geocent_peak_time.gpsSeconds) > 0)
			continue;

		/* set the burt params */
		burstParam.deltaT = series->deltaT;
		if(!(strcmp(simBurst->coordinates, "HORIZON"))) {
			burstParam.system = COORDINATESYSTEM_HORIZON;
		} else if(!(strcmp(simBurst->coordinates, "ZENITH"))) {
			/* set coordinate system for completeness */
			burstParam.system = COORDINATESYSTEM_EQUATORIAL;
			detector.site = NULL;
		} else if(!(strcmp(simBurst->coordinates, "GEOGRAPHIC"))) {
			burstParam.system = COORDINATESYSTEM_GEOGRAPHIC;
		} else if(!(strcmp(simBurst->coordinates, "EQUATORIAL"))) {
			burstParam.system = COORDINATESYSTEM_EQUATORIAL;
		} else if(!(strcmp(simBurst->coordinates, "ECLIPTIC"))) {
			burstParam.system = COORDINATESYSTEM_ECLIPTIC;
		} else if(!(strcmp(simBurst->coordinates, "GALACTIC"))) {
			burstParam.system = COORDINATESYSTEM_GALACTIC;
		} else
			burstParam.system = COORDINATESYSTEM_EQUATORIAL;

		/* generate the burst */
		memset(&waveform, 0, sizeof(CoherentGW));
		if(XLALGenerateBurst(&waveform, simBurst, &burstParam)) {
			ABORT(stat, GENERATEBURSTH_EMEM, GENERATEBURSTH_MSGEMEM);
		}

		/* must set the epoch of signal since it's used by coherent GW */
		signal.epoch = waveform.a->epoch;
		memset(signal.data->data, 0, signal.data->length * sizeof(REAL4));

		/* decide which way to calibrate the data; defaul to old way */
		if(calType)
			detector.transfer = NULL;
		else
			detector.transfer = transfer;

		/* convert this into an ADC signal */
		LALSimulateCoherentGW(stat->statusPtr, &signal, &waveform, &detector);
		CHECKSTATUSPTR(stat);

		/* if calibration using RespFilt */
		if(calType == 1)
			XLALRespFilt(&signal, transfer);

		/* inject the signal into the data channel */
		LALSSInjectTimeSeries(stat->statusPtr, series, &signal);
		CHECKSTATUSPTR(stat);

		/* free memory in coherent GW structure.  TODO:  fix this */
		LALSDestroyVectorSequence(stat->statusPtr, &(waveform.a->data));
		CHECKSTATUSPTR(stat);
		LALSDestroyVector(stat->statusPtr, &(waveform.f->data));
		CHECKSTATUSPTR(stat);
		LALDDestroyVector(stat->statusPtr, &(waveform.phi->data));
		CHECKSTATUSPTR(stat);
		LALFree(waveform.a);
		waveform.a = NULL;
		LALFree(waveform.f);
		waveform.f = NULL;
		LALFree(waveform.phi);
		waveform.phi = NULL;

		/* reset the detector site information in case it changed */
		detector.site = tmpDetector;
	}

	/* destroy the signal */
	LALSDestroyVector(stat->statusPtr, &(signal.data));
	CHECKSTATUSPTR(stat);

	LALCDestroyVector(stat->statusPtr, &(transfer->data));
	CHECKSTATUSPTR(stat);

	if(detector.site)
		LALFree(detector.site);
	LALFree(transfer);

	DETATCHSTATUSPTR(stat);
	RETURN(stat);
}
