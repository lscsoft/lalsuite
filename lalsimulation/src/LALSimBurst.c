/*
 * Copyright (C) 2007--2015 J. Creighton, K. Cannon, K. Wette, R. Prix, A. Mercer
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
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALSimBurst.h>
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
 *                              Static Functions
 *
 * ============================================================================
 */


/*
 * Fill a time series with stationary white Gaussian noise
 */


static void gaussian_noise(REAL8TimeSeries * series, double rms, gsl_rng * rng)
{
	unsigned i;

	for(i = 0; i < series->data->length; i++)
		series->data->data[i] = gsl_ran_gaussian(rng, rms);
}


/*
 * compute semimajor and semiminor axes lengths from eccentricity assuming
 * that a^2 + b^2 = 1.  eccentricity is e = \sqrt{1 - (b / a)^2}.  from
 * those two constraints the following expressions are obtained.
 */


static void semi_major_minor_from_e(double e, double *a, double *b)
{
	double e2 = e * e;

	*a = 1.0 / sqrt(2.0 - e2);
	*b = *a * sqrt(1.0 - e2);
}


/*
 * ============================================================================
 *
 *                                 Utilities
 *
 * ============================================================================
 */


/**
 * @brief Return the strain of the sample with the largest magnitude.
 *
 * @details
 * The input must have non-zero length.
 *
 * @param[in] hplus A time series.
 *
 * @param[in] hcross A time series.  Optional.  Pass NULL to disregard.
 *
 * @param[out] index The index of the peak sample will be written to this address.  Optional.  NULL to disregard.
 *
 * @returns The strain of the sample with the largest magnitude.  The \f$+\f$ polarization amplitude is in the real component, the imaginary component contains the \f$\times\f$ polarization's amplitude (or is 0 if hcross is not supplied).
 *
 * @retval XLAL_REAL8_FAIL_NAN Falure.
 */


COMPLEX16 XLALMeasureHPeak(const REAL8TimeSeries *hplus, const REAL8TimeSeries *hcross, unsigned *index)
{
	double complex hpeak;
	double abs_hpeak;
	unsigned i;

	/* check input */

	if(hcross) {
		LAL_CHECK_CONSISTENT_TIME_SERIES(hplus, hcross, XLAL_REAL8_FAIL_NAN);
	}
	if(!hplus->data->length) {
		XLALPrintError("%s(): length must be > 0\n", __func__);
		XLAL_ERROR_REAL8(XLAL_EBADLEN);
	}

	/* find peak */

	hpeak = hplus->data->data[0] + I * (hcross ? hcross->data->data[0] : 0.);
	abs_hpeak = cabs(hpeak);
	for(i = 1; i < hplus->data->length; i++)
		if(cabs(hplus->data->data[i] + I * (hcross ? hcross->data->data[i] : 0.)) > abs_hpeak) {
			hpeak = hplus->data->data[i] + I * (hcross ? hcross->data->data[i] : 0.);
			abs_hpeak = cabs(hpeak);
			if(index)
				*index = i;
		}

	return hpeak;
}


/**
 * @brief Computes the integral of the product of two time series.
 *
 * @details
 * From two time series, \f$s_{1}\f$ and \f$s_{2}\f$, computes and returns
 * \f{equation}{
 *    \int s_{1}(t) s_{2}(t) \diff t.
 * \f}
 *
 * @param[in] s1 A time series.
 *
 * @param[in] s2 A time series.
 *
 * @returns The integral of the product.
 *
 * @retval XLAL_REAL8_FAIL_NAN Failure
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
 * @brief Computes "root-sum-square strain", or \f$h_{\mathrm{rss}}\f$.
 *
 * @details In fact, this is
 * \f{equation}{
 * h_{\mathrm{rss}}
 *    = \sqrt{\sum (h_{+}^{2} + h_{x}^{2}) \Delta t},
 * \f}
 * (includes a factor of \f$\Delta t\f$), which is an approximation of the
 * square root of the square integral,  \f$\sqrt{\int (h_{+}^{2} +
 * h_{x}^{2}) \diff t}\f$.
 *
 * The input time series must start and end at the same times and have the
 * same sample rates.
 *
 * @param[in] hplus \f$h_{+}\f$ time series.
 *
 * @param[in] hcross \f$h_{\times}\f$ time series.
 *
 * @returns The \f$h_{\mathrm{rss}}\f$
 *
 * @retval XLAL_REAL8_FAIL_NAN Failure.
 */


REAL8 XLALMeasureHrss(
	const REAL8TimeSeries *hplus,
	const REAL8TimeSeries *hcross
)
{
	return sqrt(XLALMeasureIntS1S2DT(hplus, hplus) + XLALMeasureIntS1S2DT(hcross, hcross));
}


/**
 * @brief Computes the integral of the square of a real-valued time series'
 * first derivative from its Fourier transform.
 *
 * @details
 * Given the Fourier transform of a real-valued function \f$h(t)\f$,
 * compute and return the integral of the square of its derivative,
 * \f{equation}{
 * \int \stackrel{.}{h}^{2} \diff t.
 * \f}
 * The normalization factors in this function assume that
 * XLALREAL8FreqTimeFFT() will be used to convert the frequency series to
 * the time domain.
 *
 * @param[in] fseries The Fourier transform of a real-valued function of
 * time.  See also XLALREAL8TimeFreqFFT().
 *
 * @returns \f$\int \stackrel{.}{h}^{2} \diff t\f$
 *
 * @retval XLAL_REAL8_FAIL_NAN Failure.
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
 * @brief Computes the areal energy density carried by a gravitational
 * wave.
 *
 * @details
 * The local gravitational wave flux density in the two independent
 * polarizations, \f$h_{+}(t)\f$ and \f$h_{\times}(t)\f$, is [Isaacson
 * 1968]
 * \f{equation}{
 * \frac{\diff E}{\diff A \diff t}
 *    = \frac{1}{16 \pi} \frac{c^{3}}{G} \left( \dot{h}_{+}^{2} +
 *    \dot{h}_{\times}^{2} \right).
 * \f}
 * For a source at non-cosmological distances (distances small enough that
 * for spheres of that radius \f$A = 4 \pi r^{2}\f$), the equivalent
 * isotropic radiated energy in a gravitational wave for a source at a
 * distance \f$r\f$ is
 * \f{equation}{
 * E
 *    = \frac{c^{3}}{4 G} r^{2} \int \left( \dot{h}_{+}^{2}(t) +
 *    \dot{h}_{\times}^{2}(t) \right) \diff t.
 * \f}
 * Given \f$h_{+}(t)\f$ and \f$h_{\times}(t)\f$ in the waveframe, this
 * function returns
 * \f{equation}{
 * \frac{E}{r^{2}}
 *    = \frac{c^{3}}{4 G} \int ( \stackrel{.}{h}_{+}^{2} +
 *    \stackrel{.}{h}_{\times}^{2} ) \diff t.
 * \f}
 *
 * The input time series must start and end at the same times and have the
 * same sample rates.  The square integrals of the derivatives are
 * evaluated in the frequency domain so this function implicitly assumes
 * the input time series are periodic on their intervals.  The calling code
 * must ensure appropriate conditioning (e.g., tapering and padding) has
 * been applied to the time series for this to be a good approximation.
 * Waveforms that have been conditioned for use as software injections are
 * almost certainly suitably conditioned for this function.
 *
 * @param[in] hplus The \f$h_{+}\f$ time series.
 *
 * @param[in] hcross The \f$h_{\times}\f$ time series.
 *
 * @returns
 * Energy per unit area in Joules per square metre.
 *
 * @retval XLAL_REAL8_FAIL_NAN Falure.
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
 * @brief Genereates a single-sample impulse waveform
 *
 * @details
 * Places a single non-zero sample into the middle of the time series.  The
 * \f$h_{+}\f$ and \f$h_{\times}\f$ time series both have an odd number of
 * samples all set to 0 except for a single sample with amplitude hpeak in
 * the middle of the \f$h_{+}\f$ time series.
 *
 * @param[out] hplus Address of a REAL8TimeSeries pointer to be set to the
 * address of the newly allocated \f$h_{+}\f$ time series.  Set to NULL on
 * failure.
 *
 * @param[out] hcross Address of a REAL8TimeSeries pointer to be set to the
 * address of the newly allocated \f$h_{\times}\f$ time series.  Set to NULL
 * on failure.
 *
 * @param[in] hpeak Strain amplitude of the impulse.
 *
 * @param[in] delta_t Sample period of output time series in seconds.
 *
 * @retval 0 Success
 * @retval <0 Failure
 */


int XLALGenerateImpulseBurst(
	REAL8TimeSeries **hplus,
	REAL8TimeSeries **hcross,
	REAL8 hpeak,
	REAL8 delta_t
)
{
	LIGOTimeGPS epoch;
	/* length is 1 sample.  XLALSimDetectorStrainREAL8TimeSeries() will
	 * add sufficient padding to accomodate the interpolation kernel on
	 * either side */
	int length = 1;

	/* the middle sample is t = 0 */

	if(!XLALGPSSetREAL8(&epoch, -(length - 1) / 2 * delta_t))
		XLAL_ERROR(XLAL_EFUNC);

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
 * @brief Generate a band- and time-limited white-noise burst waveform
 * with Gaussian envelopes in the time and frequency domains.
 *
 * @details
 * Generates two time series containing \f$h_{+}(t)\f$ and \f$h_{x}(t)\f$,
 * with the time-domain Gaussian envelope's peak located at \f$t = 0\f$ (as
 * defined by the epoch and deltaT).  The \f$+\f$ and \f$\times\f$ time
 * series are statistically independent.
 *
The construction of a BTLWNB waveform with duration \f$\Delta t\f$ and
bandwidth \f$\Delta f\f$ centred on \f$f_{0}\f$ begins by populating a time
series with independent Gaussian random numbers.  The origin of the time
co-ordinate corresponds to the middle sample in the time series.  We apply
an initial time-limiting window function to the time series by multiplying
the time series with a Gaussian window function
\f{equation}{
w_{1}(t)
   \propto \ee^{-\frac{1}{2} t^{2} / \sigma_{t}^{2}},
\f}
where \f$\sigma_{t}\f$ sets the duration of the window.  The windowed time
series is then Fourier transformed and a second Gaussian window applied in
the frequency domain
\f{equation}{
\tilde{w}_{2}(f)
   \propto \ee^{-\frac{1}{2} (f - f_{0})^{2} / \sigma_{f}^{2}},
\f}
where \f$\sigma_{f} = \frac{1}{2} \Delta f\f$.

Since the inital time series is real-valued, the negative frequency
components of the Fourier transform are the complex conjugates of the
positive frequency components and need not be stored.  The frequency-domain
filter is real-valued (phase preserving), and so when the positive
frequency components are the only ones being stored applying the window
function to them alone achieves the correct result.

The multiplication of the frequency domain data by the window function is
equivalent to convolving the time domain data with the Fourier transform of
the window.  Since the Fourier transform of the frequency window is not a
\f$\delta\f$ function, the application of the band-limiting window has the
effect of spreading the signal in the time domain, i.e.\ increasing its
duration.  We can compensate for this by choosing an appropriate value for
\f$\sigma_{t}\f$ so that the waveform has the correct duration \e after
application of the frequency domain window.  The inverse Fourier transform
of \f$\tilde{w}_{2}(f)\f$ is
\f{equation}{
w_{2}(t)
   \propto \ee^{-2 \pi^{2} \sigma_{f}^{2} t^{2}}.
\f}
The result of convolving two Gaussians with one another is another
Gaussian, so the effective time-domain window is
\f{equation}{
w(t)
   = w_{1}(t) \otimes w_{2}(t)
   \propto \ee^{-\frac{1}{2} t^{2} / \sigma^{2}},
\f}
where
\f{equation}{
\sigma^{2}
   = \sigma_{t}^{2} + \frac{1}{4 \pi^{2} \sigma_{f}^{2}}
   = \sigma_{t}^{2} + \frac{1}{\pi^{2} \Delta f^{2}}
\f}
We wish this Gaussian's width to be \f$\sigma = \frac{1}{2} \Delta t\f$,
therefore
\f{equation}{
\sigma_{t}
   = \sqrt{\frac{1}{4} \Delta t^{2} - \frac{1}{\pi^{2} \Delta f^{2}}}.
\f}
Note that \f$\sigma_{t}\f$ is only real-valued when
\f{equation}{
\Delta t \Delta f
   \geq \frac{2}{\pi}.
\f}

After application of the frequency domain window the data is inverse
transformed to the time domain for injection into the strain data.

### Details

This function produces both \f$h_{+}\f$ and \f$h_{\times}\f$ waveforms.
These are independent waveforms constructed by applying the time series
construction algorithm twice.  The length of the result is \f$21 \Delta
t\f$ rounded to the nearest odd integer,
\f{equation}{
L
   = 2 \left\lfloor \frac{1}{2} \frac{21 \Delta t}{\delta t} \right\rfloor
   + 1
\f}
where \f$\delta t\f$ is the sample period of the time series.  The middle
sample is \f$t = 0\f$, so the first and last samples are at \f$t = \pm \delta
t (L - 1) / 2\f$.  The time-domain Gaussian window is constructed with a
call to <tt>XLALCreateGaussREAL8Window()</tt> with a shape parameter of
\f{equation}{
\beta
   = \frac{(L - 1) \delta t / 2}{\sigma_{t}}.
\f}
The numerator transforms the normalized co-ordinate \f$y \in [-1, +1]\f$ in
the definition of the window function to \f$t\f$. (See the LAL
documentation for more information.  Sample index 0 is \f$y = -1\f$, sample
index \f$L - 1\f$ is \f$y = +1\f$, so there are \f$(L - 1) / 2\f$ sample indexes
per unit of \f$y\f$.)

The time series is transformed to the frequency domain with a call to
<tt>XLALREAL8TimeFreqFFT()</tt>, which populates the metadata of the output
frequency series with the appropriate values.  There are \f$(L + 1) / 2\f$
complex-valued frequency components with a bin spacing of \f$\delta f = (L
\delta t)^{-1}\f$.  The frequency domain Gaussian window is constructed with
a call to <tt>XLALCreateGaussREAL8Window()</tt> requesting a window with a
length of \f$L + 2\f$ (twice the length of the frequency series rounded up to
the next odd integer), and a shape parameter of
\f{equation}{
\beta
   = \frac{(L + 1) \delta f / 2}{\sigma_{f}}.
\f}
The numerator in the shape parameter converts the normalized co-ordinate
\f$y \in [-1, +1]\f$ in the definition of the window function to
frequency. (See the LAL documentation for more information.  The
window has \f$L + 2\f$ samples, sample index 0 is \f$y = -1\f$, sample index
\f$L + 1\f$ is \f$y = +1\f$, so there are \f$(L + 1) / 2\f$ sample indexes per
unit of \f$y\f$.) The window is created with the peak in the middle sample
at index \f$(L + 1) / 2\f$, and we use <tt>XLALResizeREAL8Sequence()</tt> to
extract as many samples as there are in the frequency series with the peak
shifted to the correct bin.  We want the peak to be at sample index \f$f_{0}
/ \delta f\f$, so we extract \f$(L + 1) / 2\f$ samples starting at index \f$(L
+ 1) / 2 - \lfloor f_{0} / \delta f + 0.5 \rfloor\f$.

Following application of the frequency-domain window, the injection is
transformed back to the time domain with a call to
<tt>XLALREAL8FreqTimeFFT()</tt>.  If \f$\tilde{h}_{k}\f$ are the complex
values in the frequency bins, the output time series is
\f{equation}{
h_{j}
   = \delta f \sum_{k = 0}^{L - 1} \tilde{h}_{k} \ee^{2 \pi \aye j k / L}
   = \delta f \sum_{k = 0}^{L - 1} \tilde{h}_{k} \ee^{2 \pi \aye t k / (L
   \delta t)},
\f}
where \f$t = j \delta t\f$.  Differentiating with respect to \f$t\f$,
\f{equation}{
\dot{h}_{j}
   = \delta f \sum_{k = 0}^{L - 1} \left( \frac{2 \pi \aye k}{L \delta t}
   \right) \tilde{h}_{k} \ee^{2 \pi \aye j k / L},
\f}
and so
\f{align}{
\sum_{j = 0}^{L - 1} \dot{h}_{j}^{2} \delta t
   & = \delta f^{2} \delta t \sum_{k = 0}^{L - 1} \sum_{k' = 0}^{L - 1}
   \left( \frac{4 \pi^{2} k k'}{L^{2} \delta t^{2}} \right) \tilde{h}_{k}
   \conj{\tilde{h}_{k'}} \sum_{j = 0}^{L - 1} \ee^{2 \pi \aye j (k - k') /
   L}
   \\
   & = \delta f^{2} L \delta t \sum_{k = 0}^{L - 1} \left( \frac{4 \pi^{2}
   k^{2}}{L^{2} \delta t^{2}} \right) \magnitude{\tilde{h}_{k}}^{2}
   \\
   & = 4 \pi^{2} \delta f \sum_{k = 0}^{L - 1} (k \delta f)^{2}
   \magnitude{\tilde{h}_{k}}^{2}.
\f}
This relationship is used to normalize the injection time series.  The
expression on the left hand side is \f$\int \dot{h}^{2} \diff t\f$.  For both
polarizations the right hand side is computed in the frequency domain
following application of the Gaussian window, and the amplitudes of the
frequency components scaled prior to conversion to the time domain so that
\f$\int (\dot{h}_{+}^{2} + \dot{h}_{\times}^{2}) \diff t\f$ has the desired
value.

To ensure no discontinuities in the strain time series when the injection
is added to it, a final Tukey window is applied to the injection in the
time domain.  The Tukey window is constructed with a call to
<tt>XLALCreateTukeyREAL8Window()</tt> with a shape parameter of \f$\beta =
0.5\f$ so that the tapers span a total of 50\% of the time series.
Because the Tukey window is flat with unit amplitude in the middle, it has
no effect on the injection time series where the bulk of the energy is
concentrated, and the large tapers ensure the Tukey window induces
negligble spread of the injection in the frequency domain.  Because the
injection is normalized in the frequency domain prior to transformation to
the time domain, the application of the Tukey window does bias the
normalization slightly by reducing the total energy in the injection,
however the Tukey window's tapers start several \f$\sigma_{t}\f$ away from
the injection's peak and so this effect is negligble.

In order that the waveforms be reproducable so that an analysis can be
repeated, or the waveforms constructed multiple times for injection into
the strain data from more than one instrument, it is necessary to specify
how the initial time series of independent Gaussian random numbers is to be
constructed.  This is done by specifying the seed to be used with the
random number generator.  The random number generator is not specified, so
the same seed may produce different injections with different versions of
the code, but a seed and revision tag combination should be guaranteed to
produce the same injection.  Note also that changing the length of the
injection time series changes the number of random numbers used to
construct it, so the injection waveform also depends on the time series'
sample rate.  One has to be careful when constructing injection waveforms
for instruments with different sample rates (e.g., LIGO and VIRGO).  The
injection must be constructed at the same sample rate for both instruments
and then up- or down-sampled as needed when injected into each instrument's
time series.

\anchor xlalsimburstbtlwnb_examples
\image html lalsimburst_btlwnbexamples.svg
Example of the \f$+\f$ and \f$\times\f$ polarizations of a band- and
time-limited white-noise burst injection waveforms with different degrees
of freedom.
 *
 * @param[out] hplus Address of a REAL8TimeSeries pointer to be set to the
 * address of the newly allocated \f$h_{+}\f$ time series.  Set to NULL on
 * failure.
 *
 * @param[out] hcross Address of a REAL8TimeSeries pointer to be set to the
 * address of the newly allocated \f$h_{\times}\f$ time series.  Set to NULL
 * on failure.
 *
 * @param[in] duration Width of the Gaussian envelope in the time domain in
 * seconds.  The time domain envelope is \f$\propto \exp ( -\frac{1}{2}
 * t^{2} / \mathrm{duration}^{2} )\f$
 *
 * @param[in] frequency Centre frequency, \f$f_{0}\f$, of the Gaussian
 * envelope in the frequency domain in Hertz.
 *
 * @param[in] bandwidth Width of the Gaussian envelope in the frequency
 * domain in Hertz.  The frequency domain envelope is \f$\propto \exp (
 * -\frac{1}{2} (f - f_{0})^{2} / \mathrm{bandwidth}^{2} )\f$
 *
 * @param[in] eccentricity The eccentricity, \f$\epsilon = \sqrt{1 -
 * ({h_{0}}_{\times} / {h_{0}}_{+})^{2}}\f$, of the polarization ellipse
 * setting the relative amplitudes of the \f$h_{+}\f$ and \f$h_{\times}\f$
 * components' Gaussian envelopes.  With eccentricity = 0 the two
 * components have equal amplitudes (circularly polarized); with
 * eccentricity = 1 the amplitude of the \f$h_{\times}\f$ component is 0
 * (linearly polarized).  Note that this controls the relationship between
 * the expected amplitudes, not the realized amplitudes.
 *
 * @param[in] phase The phase, \f$\phi\f$, of the sinusoidal oscillations
 * that get multiplied by the Gaussian envelope.  With \f$\phi=0\f$,
 * \f$h_{+}\f$ is cosine-like and \f$h_{\times}\f$ is sine-like.  With
 * \f$\phi=\pi/2\f$, \f$h_{+}\f$ is sine-like and \f$h_{\times}\f$ is
 * cosine-like.
 *
 * @param[in] int_hdot_squared The output is normalized so that \f$\int
 * (\stackrel{.}{h}_{+}^{2} + \stackrel{.}{h}_{\times}^{2}) \diff t\f$
 * equals this.  Note that the normalization is not on the expected
 * amplitude of the waveform but on the realized amplitude of the waveform.
 *
 * @param[in] delta_t Sample period of output time series in seconds.
 *
 * @param[in] rng GSL random number generator instance.  Will be used to
 * generate normally distributed random variables to seed the
 * \f$h_{+}(t)\f$ and \f$h_{x}(t)\f$ components.
 *
 * @retval 0 Success
 * @retval <0 Failure
 *
 * @note
 * Because the injection is constructed with a random number generator, any
 * changes to this function that change how random numbers are chosen will
 * indirectly have the effect of altering the relationship between
 * injection waveform and random number seed.  For example, increasing the
 * length of the time series will change the injection waveforms.  There's
 * nothing wrong with this, the waveforms are still correct, but if there
 * is a need to reproduce a waveform exactly then it will be necessary to
 * tag the code before making such changes.
 *
 * @note
 * The algorithm's low degree-of-freedom limit is equivalent to
 * XLALSimBurstSineGaussian() but instead of Q and centre frequency the
 * duration and centre frequency are the degrees of freedom, which allows
 * the algorithm to also be evaluated in the low-frequency limit where
 * (when eccentricity = 1) it yields output equivalent to
 * XLALSimBurstGaussian().  If 2-degree-of-freedom waveforms or Gaussian
 * waveforms are required, the other functions are substantially more
 * efficient ways to generate them, but this function provides an interface
 * that yields sine-Gaussian family waveforms that is valid in all regimes.
 */


int XLALGenerateBandAndTimeLimitedWhiteNoiseBurst(
	REAL8TimeSeries **hplus,
	REAL8TimeSeries **hcross,
	REAL8 duration,
	REAL8 frequency,
	REAL8 bandwidth,
	REAL8 eccentricity,
	REAL8 phase,
	REAL8 int_hdot_squared,
	REAL8 delta_t,
	gsl_rng *rng
)
{
	int length;
	double a, b;
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

	if(duration < 0 || frequency < 0 || bandwidth < 0 || eccentricity < 0 || eccentricity > 1 || sigma_t_squared < 0 || int_hdot_squared < 0 || delta_t <= 0) {
		XLALPrintError("%s(): invalid input parameters\n", __func__);
		*hplus = *hcross = NULL;
		XLAL_ERROR(XLAL_EINVAL);
	}

	/* length of the injection time series is 21 * duration, rounded to
	 * the nearest odd integer.  this length is chosen because it works
	 * well for sine-Gaussians, but I have no metric for testing the
	 * quality of the result here. */

	length = (int) floor(21.0 * duration / delta_t / 2.0);
	length = 2 * length + 1;

	/* the middle sample is t = 0 */

	if(!XLALGPSSetREAL8(&epoch, -(length - 1) / 2 * delta_t))
		XLAL_ERROR(XLAL_EFUNC);

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
	 * shape parameter is ((length - 1) / 2) / (\sigma_{t} / delta_t) where
	 * \sigma_{t} is the compensated time-domain window duration */

	window = XLALCreateGaussREAL8Window((*hplus)->data->length, (((*hplus)->data->length - 1) / 2) / (sqrt(sigma_t_squared) / delta_t));
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
	 * frequency.  we also apply the eccentricity amplitude adjustments
	 * at this stage (last chance before the overall normalization is
	 * computed). */

	semi_major_minor_from_e(eccentricity, &a, &b);
	{
	double beta = -0.5 / (bandwidth * bandwidth / 4.);
	for(i = 0; i < tilde_hplus->data->length; i++) {
		double f = (tilde_hplus->f0 + i * tilde_hplus->deltaF) - frequency;
		double w = f == 0. ? 1. : exp(f * f * beta);
		tilde_hplus->data->data[i] *= a * w;
		tilde_hcross->data->data[i] *= b * w;
		/* rotate phases of non-DC components */
		if(i != 0) {
			tilde_hplus->data->data[i] *= cexp(-I * phase);
			tilde_hcross->data->data[i] *= I * cexp(-I * phase);
		}
	}
	}

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
 *                         SineGaussian and Friends
 *
 * ============================================================================
 */

/**
 * @brief Compute the Q of a sine-Gaussian waveform from the duration and
 * centre frequency.
 *
 * @details The relationship is
 * \f{equation}{
 * Q
 *    = 2 \pi f_{0} \Delta t.
 * \f}
 * The result becomes independent of duration at 0 Hz.
 *
 * @param[in] duration The duration, \f$\Delta t\f$, of the sine-Gaussian in
 * seconds.
 *
 * @param[in] centre_frequency The centre frequency, \f$f_{0}\f$, of the
 * sine-Gaussian in Hertz.
 *
 * @retval Q The \f$Q\f$ of the sine-Gaussian.
 *
 * See also:  XLALSimBurstSineGaussianDuration()
 */


double XLALSimBurstSineGaussianQ(
	double duration,
	double centre_frequency
)
{
	return LAL_TWOPI * duration * centre_frequency;
}


/**
 * @brief Compute the duration of a sine-Gaussian waveform from the Q and centre
 * frequency.
 *
 * @details The relationship is
 * \f{equation}{
 * Q
 *    = 2 \pi f_{0} \Delta t.
 * \f}
 * The relationship is undefined at 0 Hz.
 *
 * @param[in] Q The \f$Q\f$ of the sine-Gaussian.
 *
 * @param[in] centre_frequency The centre frequency, \f$f_{0}\f$, of the
 * sine-Gaussian in Hertz.
 *
 * @retval duration The duration of the sine-Gaussian, \f$\Delta t\f$, in
 * seconds.
 *
 * See also:  XLALSimBurstSineGaussianQ()
 */


double XLALSimBurstSineGaussianDuration(
	double Q,
	double centre_frequency
)
{
	double duration = Q / (LAL_TWOPI * centre_frequency);
	if(!isfinite(duration))
		XLAL_ERROR_REAL8(XLAL_EDOM);
	return duration;
}


/**
 * @brief Generate sine- and cosine-Gaussian waveforms with various
 * polarizations and phases.
 *
 * @details
 * Generates two time series, \f$h_{+}\f$ and \f$h_{\times}\f$, containing
 * add-mixtures of cosine-Gaussian and sine-Gaussian waveforms.  The
 * Gaussian envelope peaks in both at t = 0 as defined by epoch and deltaT.
 * By setting the eccentricity and phase to appropriate values any
 * linearly, elliptically, or cicularly polarized sine- or cosine-Gaussian
 * waveform can be generated.  The dominant polarization is placed in the
 * \f$h_{+}\f$ component.
 *
 * A Tukey window is applied to make the waveform go to 0 smoothly at the
 * start and end.
 *
 * \anchor xlalsimburstsinegaussian_examples
 * \image html lalsimburst_sinegaussianexamples.svg "Sine-Gaussian examples."
 *
 * @param[out] hplus Address of a REAL8TimeSeries pointer to be set to the
 * address of the newly allocated \f$h_{+}\f$ time series.  Set to NULL on
 * failure.
 *
 * @param[out] hcross Address of a REAL8TimeSeries pointer to be set to the
 * address of the newly allocated \f$h_{\times}\f$ time series.  Set to NULL
 * on failure.
 *
 * @param[in] Q The "Q" of the waveform.  The Gaussian envelope is
 * \f$\propto \exp(-\frac{1}{2} t^{2} / \sigma_{t}^{2})\f$ where
 * \f$\sigma_{t} = Q / (2 \pi f)\f$.  See also XLALSimBurstSineGaussianQ()
 * and XLALSimBurstSineGaussianDuration().
 *
 * @param[in] centre_frequency The frequency of the sinusoidal oscillations
 * that get multiplied by the Gaussian envelope.
 *
 * @param[in] hrss The \f$h_{\mathrm{rss}}\f$ of the waveform to be
 * generated.  See K. Riles, LIGO-T040055-00.pdf.  This function normalizes
 * the waveform algebraically assuming it to be an ideal sine-Gaussian
 * (continuous in time, with no time boundaries and no tapering window), so
 * the actual numerical normalization might be slightly different.  See
 * also XLALMeasureHrss().
 *
 * @param[in] eccentricity The eccentricity, \f$\epsilon = \sqrt{1 -
 * ({h_{0}}_{\times} / {h_{0}}_{+})^{2}}\f$, of the polarization ellipse
 * setting the relative amplitudes of the \f$h_{+}\f$ and \f$h_{\times}\f$
 * components' Gaussian envelopes.  With eccentricity = 0 the two
 * components have equal amplitudes (circularly polarized); with
 * eccentricity = 1 the amplitude of the \f$h_{\times}\f$ component is 0
 * (linearly polarized).
 *
 * @param[in] phase The phase, \f$\phi\f$, of the sinusoidal oscillations
 * that get multiplied by the Gaussian envelope.  With \f$\phi=0\f$,
 * \f$h_{+}\f$ is cosine-like and \f$h_{\times}\f$ is sine-like.  With
 * \f$\phi=\pi/2\f$, \f$h_{+}\f$ is sine-like and \f$h_{\times}\f$ is
 * cosine-like.
 *
 * @param[in] delta_t Sample period of output time series in seconds.
 *
 * @retval 0 Success
 * @retval <0 Failure
 */
int XLALSimBurstSineGaussian(
	REAL8TimeSeries **hplus,
	REAL8TimeSeries **hcross,
	REAL8 Q,
	REAL8 centre_frequency,
	REAL8 hrss,
	REAL8 eccentricity,
	REAL8 phase,
	REAL8 delta_t
)
{
	REAL8Window *window;
	/* square integral of unit amplitude cosine- and sine-Gaussian
	 * waveforms.  the sine-Gaussian case is derived in K. Riles,
	 * LIGO-T040055-00.pdf, equation (7).  the cosine-Gaussian case is
	 * obtained by replacing cos^2 with 1-sin^2, using equation (5) and
	 * the result for sine-Gaussians. */
	const double cgsq = Q / (4.0 * centre_frequency * sqrt(LAL_PI)) * (1.0 + exp(-Q * Q));
	const double sgsq = Q / (4.0 * centre_frequency * sqrt(LAL_PI)) * (1.0 - exp(-Q * Q));
	/* semimajor and semiminor axes of waveform ellipsoid. */
	double a, b;
	semi_major_minor_from_e(eccentricity, &a, &b);
	/* peak amplitudes of plus and cross */
	double cosphase = cos(phase);
	double sinphase = sin(phase);
	const double h0plus  = hrss * a / sqrt(cgsq * cosphase * cosphase + sgsq * sinphase * sinphase);
	const double h0cross = hrss * b / sqrt(cgsq * sinphase * sinphase + sgsq * cosphase * cosphase);
	LIGOTimeGPS epoch;
	int length;
	unsigned i;
	/* don't compute these in loops */
	const double negative2Qsquared = -2. * Q * Q;
	const double twopif0 = LAL_TWOPI * centre_frequency;
	/* some pointers */
	double *hp, *hc, *w;

	/* check input. */

	if(Q < 0 || centre_frequency < 0 || hrss < 0 || eccentricity < 0 || eccentricity > 1 || delta_t <= 0) {
		XLALPrintError("%s(): invalid input parameters\n", __func__);
		*hplus = *hcross = NULL;
		XLAL_ERROR(XLAL_EINVAL);
	}

	/* length of the injection time series is 21 * the width of the
	 * Gaussian envelope (sigma_t in the comments above), rounded to
	 * the nearest odd integer.  experiments suggest that that's the
	 * minimum length without the hrss of the output deviating from the
	 * requested hrss by more than numerical noise. */

	length = (int) floor(21.0 * Q / centre_frequency / LAL_TWOPI / delta_t / 2.0);
	length = 2 * length + 1;

	/* the middle sample is t = 0 */

	if(!XLALGPSSetREAL8(&epoch, -(length - 1) / 2 * delta_t))
		XLAL_ERROR(XLAL_EFUNC);

	/* allocate the time series */

	*hplus = XLALCreateREAL8TimeSeries("sine-Gaussian +", &epoch, 0.0, delta_t, &lalStrainUnit, length);
	*hcross = XLALCreateREAL8TimeSeries("sine-Gaussian x", &epoch, 0.0, delta_t, &lalStrainUnit, length);
	window = XLALCreateTukeyREAL8Window((*hplus)->data->length, 0.5);
	if(!*hplus || !*hcross || !window) {
		XLALDestroyREAL8TimeSeries(*hplus);
		XLALDestroyREAL8TimeSeries(*hcross);
		XLALDestroyREAL8Window(window);
		*hplus = *hcross = NULL;
		XLAL_ERROR(XLAL_EFUNC);
	}

	/* populate */

	hp = (*hplus)->data->data;
	hc = (*hcross)->data->data;
	w = window->data->data;
	for(i = 0; i < (*hplus)->data->length; i++) {
		const double t = ((int) i - (length - 1) / 2) * delta_t;
		const double phi = twopif0 * t;
		const complex double fac = cexp(phi * phi / negative2Qsquared + I * (phi - phase)) * w[i];
		hp[i] = h0plus * creal(fac);
		hc[i] = h0cross * cimag(fac);
	}
	XLALDestroyREAL8Window(window);

	/* done */

	return 0;
}


/**
 * @brief Generate Gaussian waveforms.
 *
 * @details
 * The burst working group has traditionally treated these as a distinct
 * class of waveform rather than, say, the low-frequency limit of the
 * sine-Gaussian class of waveform.  Therefore, for convenience, a separate
 * interface is provided to generate these waveforms.
 *
 * Generates two time series, \f$h_{+}\f$ and \f$h_{\times}\f$, containing
 * a Gaussian in \f$h_{+}\f$.  \f$h_{\times}\f$ is set to 0.  The Gaussian
 * peaks at t = 0 as defined by epoch and deltaT.  The degrees of freedom
 * are the duration and the \f$h_{\mathrm{rss}}\f$.  The function is
 * \f{equation}{
 * h_{+}(t)
 *    = \frac{h_{\mathrm{rss}}}{\sqrt{\sqrt{\pi} \Delta t}} \exp -\frac{1}{2} \frac{t^{2}}{\Delta t^{2}}.
 * \f}
 *
 * A Tukey window is applied to make the waveform go to 0 smoothly at the
 * start and end.
 *
 * @param[out] hplus Address of a REAL8TimeSeries pointer to be set to the
 * address of the newly allocated \f$h_{+}\f$ time series.  Set to NULL on
 * failure.
 *
 * @param[out] hcross Address of a REAL8TimeSeries pointer to be set to the
 * address of the newly allocated \f$h_{\times}\f$ time series.  Set to NULL
 * on failure.
 *
 * @param[in] duration The width of the Gaussian, \f$\Delta t\f$.
 *
 * @param[in] hrss The \f$h_{\mathrm{rss}}\f$ of the waveform to be
 * generated.  This function normalizes the waveform algebraically assuming
 * it to be an ideal Gaussian (continuous in time, with no time boundaries
 * and no tapering window), so the actual numerical normalization might be
 * slightly different.  See also XLALMeasureHrss().
 *
 * @param[in] delta_t Sample period of output time series in seconds.
 *
 * @retval 0 Success
 * @retval <0 Failure
 */


int XLALSimBurstGaussian(
	REAL8TimeSeries **hplus,
	REAL8TimeSeries **hcross,
	REAL8 duration,
	REAL8 hrss,
	REAL8 delta_t
)
{
	REAL8Window *window;
	const double h0plus  = hrss / sqrt(sqrt(LAL_PI) * duration);
	LIGOTimeGPS epoch;
	int i, length;

	/* check input. */

	if(duration < 0 || hrss < 0 || !isfinite(h0plus) || delta_t <= 0) {
		XLALPrintError("%s(): invalid input parameters\n", __func__);
		*hplus = *hcross = NULL;
		XLAL_ERROR(XLAL_EINVAL);
	}

	/* length of the injection time series is 21 * the width of the
	 * Gaussian envelope, because that's what works well for
	 * sine-Gaussians */

	length = (int) floor(21.0 * duration / delta_t / 2.0);
	length = 2 * length + 1;

	/* the middle sample is t = 0 */

	if(!XLALGPSSetREAL8(&epoch, -(length - 1) / 2 * delta_t))
		XLAL_ERROR(XLAL_EFUNC);

	/* allocate the time series */

	*hplus = XLALCreateREAL8TimeSeries("Gaussian +", &epoch, 0.0, delta_t, &lalStrainUnit, length);
	*hcross = XLALCreateREAL8TimeSeries("Gaussian x", &epoch, 0.0, delta_t, &lalStrainUnit, length);
	window = XLALCreateTukeyREAL8Window(length, 0.5);
	if(!*hplus || !*hcross || !window) {
		XLALDestroyREAL8TimeSeries(*hplus);
		XLALDestroyREAL8TimeSeries(*hcross);
		XLALDestroyREAL8Window(window);
		*hplus = *hcross = NULL;
		XLAL_ERROR(XLAL_EFUNC);
	}

	/* populate */

	for(i = 0; i < (length - 1) / 2; i++) {
		const double t = ((int) i - (length - 1) / 2) * delta_t;
		(*hplus)->data->data[i] = (*hplus)->data->data[length - 1 - i] = h0plus * exp(-0.5 * t * t / (duration * duration)) * window->data->data[i];
	}
	(*hplus)->data->data[i] = h0plus;
	memset((*hcross)->data->data, 0, (*hcross)->data->length * sizeof(*(*hcross)->data->data));

	XLALDestroyREAL8Window(window);

	/* done */

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
 * @brief Generates cosmic string cusp waveforms.
 *
 * @details
 * Generates the \f$h_{+}\f$ and \f$h_{\times}\f$ components of a cosmic
 * string cusp waveform.  These waveforms are linearly polarized and placed
 * in the \f$h_{+}\f$ compnent.  The \f$h_{\times}\f$ component is set to
 * 0.  The waveform peaks at t = 0 (as defined by the epoch and deltaT).
 *
 * In the frequency domain, the waveform is \f$A f^{-\frac{4}{3}}\f$ with a
 * (non-physical) low-frequency cut-off and a (physical) high-frequency
 * cut-off.
 * \f{equation}{
 * \tilde{h}_{+}(f)
 *    = A f^{-\frac{4}{3}} \left(1 +
 *    \frac{f_{\mathrm{low}}^{2}}{f^{2}}\right)^{-4} \begin{cases} \exp(1 -
 *    f/f_{\mathrm{high}}) & f > f_{\mathrm{high}} \\ 1 & f \leq
 *    f_{\mathrm{high}} \end{cases}
 * \f}
 *
 * The output has a Tukey window applied to force it to go to 0 smoothly at
 * the start and end.  The low frequnecy cut-off is fixed at
 * \f$f_{\mathrm{low}} = 1 \mathrm{Hz}\f$, so these waveforms should be
 * high-pass filtered before being used as injections or search templates.
 *
 * \anchor xlalgeneratestringcusp_examples
 * \image html lalsimburst_stringcuspexamples.svg "String cusp examples."
 *
 * @param[out] hplus Address of a REAL8TimeSeries pointer to be set to the
 * address of the newly allocated \f$h_{+}\f$ time series.  Set to NULL on
 * failure.
 *
 * @param[out] hcross Address of a REAL8TimeSeries pointer to be set to the
 * address of the newly allocated \f$h_{\times}\f$ time series.  Set to NULL
 * on failure.
 *
 * @param[in] amplitude Waveform's amplitude parameter, \f$A\f$, in units
 * of \f$\mathrm{strain}\,\mathrm{s}^{-\frac{1}{3}}\f$.
 *
 * @param[in] f_high High frequency cut-off, \f$f_{\mathrm{high}}\f$, in
 * Hertz.
 *
 * @param[in] delta_t Sample period of output time series in seconds.
 *
 * @retval 0 Success
 * @retval <0 Failure
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
	REAL8Window *window;
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

	/* length of the injection time series is 9 / f_low, rounded to
	 * the nearest odd integer.  at that length the waveform's
	 * amplitude has decayed to the level of numerical noise in the FFT
	 * so there's no advantage in making it longer. */

	length = (int) (9.0 / f_low / delta_t / 2.0);
	length = 2 * length + 1;

	/* the middle sample is t = 0 */

	if(!XLALGPSSetREAL8(&epoch, -(length - 1) / 2 * delta_t))
		XLAL_ERROR(XLAL_EFUNC);

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

		double amp = amplitude * pow(1. + f_low * f_low / (f * f), -4.) * pow(f, -4. / 3.) * (f > f_high ? exp(1. - f / f_high) : 1.);

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
	for(i = 0; i < (int) window->data->length; i++)
		(*hplus)->data->data[i] *= window->data->data[i];
	XLALDestroyREAL8Window(window);

	/* done */

	return 0;
}
