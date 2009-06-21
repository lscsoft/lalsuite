/*
 * $Id$
 *
 * Copyright (C) 2007  Kipp Cannon and Flanagan, E
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */


#include <math.h>


#include <gsl/gsl_matrix.h>


#include <lal/FrequencySeries.h>
#include <lal/LALAtomicDatatypes.h>
#include <lal/LALMalloc.h>
#include <lal/LALStdlib.h>
#include <lal/RealFFT.h>
#include <lal/Sequence.h>
#include <lal/TFTransform.h>
#include <lal/Units.h>
#include <lal/Window.h>
#include <lal/XLALError.h>


#include <lal/LALRCSID.h>
NRCSID(CREATETFPLANEC, "$Id$");


/*
 * ============================================================================
 *
 *                             Timing Arithmetic
 *
 * ============================================================================
 */


/*
 * Round target_length down so that an integer number of intervals of
 * length segment_length, each shifted by segment_shift with respect to the
 * interval preceding it, fits into the result.
 */


INT4 XLALOverlappedSegmentsCommensurate(
	INT4 target_length,
	INT4 segment_length,
	INT4 segment_shift
)
{
	static const char func[] = "XLALOverlappedSegmentsCommensurate";
	UINT4 segments;

	/*
	 * check input
	 */

	if(segment_length < 1) {
		XLALPrintError("segment_length < 1");
		XLAL_ERROR(func, XLAL_EINVAL);
	}
	if(segment_shift < 1) {
		XLALPrintError("segment_shift < 1");
		XLAL_ERROR(func, XLAL_EINVAL);
	}

	/*
	 * trivial case
	 */

	if(target_length < segment_length)
		return 0;

	/*
	 * do the arithmetic
	 */

	segments = (target_length - segment_length) / segment_shift;

	return segments * segment_shift + segment_length;
}


/*
 * Compute and return the timing parameters for an excess power analysis.
 * Pass NULL for any optional pointer to not compute and return that
 * parameter.
 *
 * Input:
 * 	window_length:
 * 		number of samples in a window used for the time-frequency
 * 		plane
 *
 * 	max_tile_length:
 * 		number of samples in the tile of longest duration
 *
 * 	fractional_tile_shift:
 * 		number of samples by which the start of the longest tile is
 * 		shifted from the start of the tile preceding it, as a
 * 		fraction of its length
 *
 * 	psd_length (optional, required for psd_shift):
 * 		user's desired number of samples to use in computing a PSD
 * 		estimate
 *
 * Output:
 * 	psd_length (optional):
 * 		actual number of samples to use in computing a PSD estimate
 * 		(rounded down to be comensurate with the windowing)
 *
 *	psd_shift (optional):
 *		number of samples by which the start of a PSD is to be
 *		shifted from the start of the PSD that preceded it in order
 *		that the tiling pattern continue smoothly across the
 *		boundary
 *
 *	window_shift (optional):
 *		number of samples by which the start of a time-frequency
 *		plane window is shifted from the window preceding it in
 *		order that the tiling pattern continue smoothly across the
 *		boundary
 *
 *	window_pad (optional):
 *		how many samples at the start and end of each window are
 *		treated as padding, and will not be covered by the tiling
 *
 *	tiling_length (optional):
 *		how many samples will be covered by the tiling
 *
 * NOTE:  this function is wrapped in the pyLAL package to teach the
 * Python-based DAG construction scripts how the search code's internal
 * timing works.  If you change this function, you need to update pyLAL.
 */


INT4 XLALEPGetTimingParameters(
	INT4 window_length,
	INT4 max_tile_length,
	REAL8 fractional_tile_shift,
	INT4 *psd_length,
	INT4 *psd_shift,
	INT4 *window_shift,
	INT4 *window_pad,
	INT4 *tiling_length
)
{
	static const char func[] = "XLALEPGetTimingParameters";
	int max_tile_shift = fractional_tile_shift * max_tile_length;
	int wpad;
	int tlength;
	int wshift;

	/*
	 * check input parameters
	 */

	if(window_length % 4 != 0) {
		XLALPrintError("window_length is not a multiple of 4");
		XLAL_ERROR(func, XLAL_EINVAL);
	}
	if(max_tile_length < 1) {
		XLALPrintError("max_tile_length < 1");
		XLAL_ERROR(func, XLAL_EINVAL);
	}
	if(fractional_tile_shift <= 0) {
		XLALPrintError("fractional_tile_shift <= 0");
		XLAL_ERROR(func, XLAL_EINVAL);
	}
	if(fmod(fractional_tile_shift * max_tile_length, 1) != 0) {
		XLALPrintError("fractional_tile_shift * max_tile_length not an integer");
		XLAL_ERROR(func, XLAL_EINVAL);
	}
	if(max_tile_shift < 1) {
		XLALPrintError("fractional_tile_shift * max_tile_length < 1");
		XLAL_ERROR(func, XLAL_EINVAL);
	}

	/*
	 * discard first and last 4096 samples
	 *
	 * FIXME.  this should be tied to the sample frequency and
	 * time-frequency plane's channel spacing.  multiplying the time
	 * series by a window has the effect of convolving the Fourier
	 * transform of the data by the F.T. of the window, and we don't
	 * want this to blur the spectrum by an amount larger than 1
	 * channel --- a delta function in the spectrum should remain
	 * confined to a single bin.  a channel width of 2 Hz means the
	 * notch feature created in the time series by the window must be
	 * at least .5 s long to not result in undesired leakage.  at a
	 * sample frequency of 8192 samples / s, it must be at least 4096
	 * samples long (2048 samples at each end of the time series).  to
	 * be safe, we double that to 4096 samples at each end.
	 */

	wpad = 4096;

	/*
	 * tiling covers the remainder, rounded down to fit an integer
	 * number of tiles
	 */

	tlength = window_length - 2 * wpad;
	tlength = XLALOverlappedSegmentsCommensurate(tlength, max_tile_length, max_tile_shift);
	if(tlength < 1) {
		XLALPrintError("tiling_length < 1");
		XLAL_ERROR(func, XLAL_EINVAL);
	}
	if(tiling_length)
		*tiling_length = tlength;

	/*
	 * now re-compute window_pad from rounded-off tiling_length
	 */

	wpad = (window_length - tlength) / 2;
	if(tlength + 2 * wpad != window_length) {
		XLALPrintError("cannot find window parameters consistent with tiling parameters");
		XLAL_ERROR(func, XLAL_EINVAL);
	}
	if(window_pad)
		*window_pad = wpad;

	/*
	 * adjacent tilings overlap so that their largest tiles overlap the
	 * same as within each tiling
	 */

	wshift = tlength - (max_tile_length - max_tile_shift);
	if(wshift < 1) {
		XLALPrintError("window_shift < 1");
		XLAL_ERROR(func, XLAL_EINVAL);
	}
	if(window_shift)
		*window_shift = wshift;

	/*
	 * compute the adjusted PSD length if desired
	 */

	if(psd_length) {
		*psd_length = XLALOverlappedSegmentsCommensurate(*psd_length, window_length, wshift);
		if(*psd_length < 0)
			XLAL_ERROR(func, XLAL_EFUNC);

		if(psd_shift) {
			*psd_shift = *psd_length - (window_length - wshift);
			if(*psd_shift < 1) {
				XLALPrintError("psd_shift < 1");
				XLAL_ERROR(func, XLAL_EINVAL);
			}
		}
	} else if(psd_shift) {
		/* for safety */
		*psd_shift = -1;
		/* can't compute psd_shift without psd_length input */
		XLAL_ERROR(func, XLAL_EFAULT);
	}

	return 0;
}


/*
 * ============================================================================
 *
 *                   Time-Frequency Plane Create / Destroy
 *
 * ============================================================================
 */


/**
 * Compute the two-point spectral correlation function for a whitened
 * frequency series from the window applied to the original time series.
 *
 * If x_{j} is a stationary process then the components of its Fourier
 * transform, X_{k}, are independent random variables, and let their mean
 * square be <|X_{k}|^{2}> = 1.  If x_{j} is multiplied by the window
 * function w_{j} then it is no longer stationary and the components of its
 * Fourier transform are no longer independent.  Their correlations are
 *
 *	<X_{k} X*_{k'}>
 *
 * and depend only on |k - k'|.
 *
 * Given the window function w_{j}, this function computes and returns a
 * sequence containing <X_{k} X*_{k'}>.  The sequence's indices are |k -
 * k'|.  A straight-forward normalization factor can be applied to convert
 * this for use with a sequence x_{j} whose Fourier transform does not have
 * bins with equal mean square.
 *
 * The FFT plan argument must be a forward plan (time to frequency) whose
 * length equals that of the window.  If the window has length N the return
 * value is the address of a newly-allocated sequence of length floor(N/2 +
 * 1), or NULL on error.  This function assumes the window is symmetric
 * about its midpoint (is an even function of the sample index if the
 * midpoint is index 0), so that its Fourier transform is real-valued.
 */


static REAL8Sequence *XLALREAL8WindowTwoPointSpectralCorrelation(
	const REAL8Window *window,
	const REAL8FFTPlan *plan
)
{
	static const char func[] = "XLALREAL8WindowTwoPointSpectralCorrelation";
	REAL8Sequence *wsquared;
	COMPLEX16Sequence *tilde_wsquared;
	REAL8Sequence *correlation;
	unsigned i;

	if(window->sumofsquares <= 0)
		XLAL_ERROR_NULL(func, XLAL_EDOM);

	/*
	 * Create a sequence to hold the normalized square of the window
	 * and its Fourier transform.
	 */

	wsquared = XLALCopyREAL8Sequence(window->data);
	tilde_wsquared = XLALCreateCOMPLEX16Sequence(window->data->length / 2 + 1);
	if(!wsquared || !tilde_wsquared) {
		XLALDestroyREAL8Sequence(wsquared);
		XLALDestroyCOMPLEX16Sequence(tilde_wsquared);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}

	/*
	 * Compute the normalized square of the window.
	 */

	for(i = 0; i < wsquared->length; i++)
		wsquared->data[i] *= wsquared->data[i] / window->sumofsquares;

	/*
	 * Fourier transform.
	 */

	if(XLALREAL8ForwardFFT(tilde_wsquared, wsquared, plan)) {
		XLALDestroyREAL8Sequence(wsquared);
		XLALDestroyCOMPLEX16Sequence(tilde_wsquared);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}
	XLALDestroyREAL8Sequence(wsquared);

	/*
	 * Create space to hold the two-point correlation function.
	 */

	correlation = XLALCreateREAL8Sequence(tilde_wsquared->length);
	if(!correlation) {
		XLALDestroyCOMPLEX16Sequence(tilde_wsquared);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}

	/*
	 * Extract real components from Fourier transform.
	 */

	for(i = 0; i < correlation->length; i++)
		correlation->data[i] = tilde_wsquared->data[i].re;
	XLALDestroyCOMPLEX16Sequence(tilde_wsquared);

	/*
	 * Done.
	 */

	return correlation;
}


/*
 * Create and initialize a time-frequency plane object.
 */


REAL8TimeFrequencyPlane *XLALCreateTFPlane(
	/* length of time series from which TF plane will be computed */
	UINT4 tseries_length,
	/* sample rate of time series */
	REAL8 tseries_deltaT,
	/* minimum frequency to search for */
	REAL8 flow,
	/* bandwidth of TF plane */
	REAL8 bandwidth,
	/* overlap of adjacent tiles */
	REAL8 tiling_fractional_stride,
	/* largest tile's bandwidth */
	REAL8 max_tile_bandwidth,
	/* largest tile's duration */
	REAL8 max_tile_duration,
	/* forward plan whose length is tseries_length */
	const REAL8FFTPlan *plan
)
{
	static const char func[] = "XLALCreateTFPlane";
	REAL8TimeFrequencyPlane *plane;
	gsl_matrix *channel_data;
	REAL8Sequence *channel_buffer;
	REAL8Sequence *unwhitened_channel_buffer;
	REAL8Window *tukey;
	REAL8Sequence *correlation;

	/*
	 * resolution of FT of input time series
	 */

	const double fseries_deltaF = 1.0 / (tseries_length * tseries_deltaT);

	/*
	 * time-frequency plane's channel spacing
	 */

	const double deltaF = 1 / max_tile_duration * tiling_fractional_stride;

	/*
	 * total number of channels
	 */

	const int channels = floor(bandwidth / deltaF + 0.5);

	/*
	 * stride
	 */

	const unsigned inv_fractional_stride = floor(1.0 / tiling_fractional_stride + 0.5);

	/*
	 * tile size limits
	 */

	const unsigned min_length = floor((1 / max_tile_bandwidth) / tseries_deltaT + 0.5);
	const unsigned max_length = floor(max_tile_duration / tseries_deltaT + 0.5);
	const unsigned min_channels = inv_fractional_stride;
	const unsigned max_channels = floor(max_tile_bandwidth / deltaF + 0.5);

	/*
	 * sample on which tiling starts
	 */

	int tiling_start;

	/*
	 * length of tiling
	 */

	int tiling_length;

	/*
	 * window shift
	 */

	int window_shift;

	/*
	 * Compute window_shift, tiling_start, and tiling_length.
	 */

	if(XLALEPGetTimingParameters(tseries_length, max_tile_duration / tseries_deltaT, tiling_fractional_stride, NULL, NULL, &window_shift, &tiling_start, &tiling_length) < 0)
		XLAL_ERROR_NULL(func, XLAL_EFUNC);

	/*
	 * Make sure that input parameters are reasonable, and that a
	 * complete tiling is possible.
	 *
	 * Note that because all tile durations are integer power of two
	 * multiples of the smallest duration, if the largest duration fits
	 * an integer number of times in the tiling length, then all tile
	 * sizes do so there's no need to test them all.  Likewise for the
	 * bandwidths.
	 *
	 * FIXME:  these tests require an integer number of non-overlapping
	 * tiles to fit, which is stricter than required;  only need an
	 * integer number of overlapping tiles to fit, but then probably
	 * have to test all sizes separately.
	 */

	if((flow < 0) ||
	   (bandwidth <= 0) ||
	   (deltaF <= 0) ||
	   (inv_fractional_stride * tiling_fractional_stride != 1) ||
	   (fmod(max_tile_duration, tseries_deltaT) != 0) ||
	   (fmod(deltaF, fseries_deltaF) != 0) ||
	   (tseries_deltaT <= 0) ||
	   (channels * deltaF != bandwidth) ||
	   (min_length * tseries_deltaT != (1 / max_tile_bandwidth)) ||
	   (min_length % inv_fractional_stride != 0) ||
	   (tiling_length % max_length != 0) ||
	   (channels % max_channels != 0)) {
		XLALPrintError("unable to construct time-frequency tiling from input parameters\n");
		XLAL_ERROR_NULL(func, XLAL_EINVAL);
	}

	/*
	 * Allocate memory.
	 */

	plane = XLALMalloc(sizeof(*plane));
	channel_data = gsl_matrix_alloc(tseries_length, channels);
	channel_buffer = XLALCreateREAL8Sequence(tseries_length);
	unwhitened_channel_buffer = XLALCreateREAL8Sequence(tseries_length);
	tukey = XLALCreateTukeyREAL8Window(tseries_length, (tseries_length - tiling_length) / (double) tseries_length);
	if(tukey)
		correlation = XLALREAL8WindowTwoPointSpectralCorrelation(tukey, plan);
	else
		/* error path */
		correlation = NULL;
	if(!plane || !channel_data || !channel_buffer || !unwhitened_channel_buffer || !tukey || !correlation) {
		XLALFree(plane);
		if(channel_data)
			gsl_matrix_free(channel_data);
		XLALDestroyREAL8Sequence(channel_buffer);
		XLALDestroyREAL8Sequence(unwhitened_channel_buffer);
		XLALDestroyREAL8Window(tukey);
		XLALDestroyREAL8Sequence(correlation);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}

	/*
	 * Initialize the structure
	 */

	plane->name[0] = '\0';
	plane->epoch.gpsSeconds = 0;
	plane->epoch.gpsNanoSeconds = 0;
	plane->deltaT = tseries_deltaT;
	plane->fseries_deltaF = fseries_deltaF;
	plane->deltaF = deltaF;
	plane->flow = flow;
	plane->channel_data = channel_data;
	plane->channel_buffer = channel_buffer;
	plane->unwhitened_channel_buffer = unwhitened_channel_buffer;
	plane->tiles.max_length = max_length;
	plane->tiles.min_channels = min_channels;
	plane->tiles.max_channels = max_channels;
	plane->tiles.tiling_start = tiling_start;
	plane->tiles.tiling_end = tiling_start + tiling_length;
	plane->tiles.inv_fractional_stride = inv_fractional_stride;
	plane->tiles.dof_per_pixel = 2 * tseries_deltaT * deltaF;
	plane->window = tukey;
	plane->window_shift = window_shift;
	plane->two_point_spectral_correlation = correlation;

	/*
	 * Success
	 */

	return plane;
}


/*
 * Free a time-frequency plane object.
 */


void XLALDestroyTFPlane(
	REAL8TimeFrequencyPlane *plane
)
{
	if(plane) {
		if(plane->channel_data)
			gsl_matrix_free(plane->channel_data);
		XLALDestroyREAL8Sequence(plane->channel_buffer);
		XLALDestroyREAL8Sequence(plane->unwhitened_channel_buffer);
		XLALDestroyREAL8Window(plane->window);
		XLALDestroyREAL8Sequence(plane->two_point_spectral_correlation);
	}
	XLALFree(plane);
}


/*
 * ============================================================================
 *
 *                         Channel Filter Management
 *
 * ============================================================================
 */


/*
 * Compute the magnitude of the inner product of two arbitrary channel
 * filters.  Note that the sums are done over only the positive frequency
 * components, so this function multiplies by the required factor of 2.
 * The result is the *full* inner product, not the half inner product.  It
 * is safe to pass the same filter as both arguments.
 *
 * The second version computes the PSD-weighted inner product.  This is
 * used in reconstructing h_rss.
 *
 * NOTE:  checks are not done to ensure that the delta fs are all equal,
 * and that the PSD spans the frequency range of the two filters.   these
 * things are implicit in the code that calls these functions so it would
 * just burn CPU time to do them, but be aware of these requirements if
 * this code is copied and used somewhere else!
 */


static REAL8 filter_inner_product(
	const COMPLEX16FrequencySeries *filter1,
	const COMPLEX16FrequencySeries *filter2,
	const REAL8Sequence *correlation
)
{
	const int k10 = floor(filter1->f0 / filter1->deltaF + 0.5);
	const int k20 = floor(filter2->f0 / filter2->deltaF + 0.5);
	int k1, k2;
	COMPLEX16 sum = {0, 0};

	for(k1 = 0; k1 < (int) filter1->data->length; k1++) {
		const COMPLEX16 *f1data = &filter1->data->data[k1];
		for(k2 = 0; k2 < (int) filter2->data->length; k2++) {
			const COMPLEX16 *f2data = &filter2->data->data[k2];
			const unsigned delta_k = abs(k10 + k1 - k20 - k2);
			const double sksk = (delta_k & 1 ? -1 : +1) * (delta_k < correlation->length ? correlation->data[delta_k] : 0);

			sum.re += sksk * (f1data->re * f2data->re + f1data->im * f2data->im);
			sum.im += sksk * (f1data->im * f2data->re - f1data->re * f2data->im);
		}
	}

	return 2 * sqrt(sum.re * sum.re + sum.im * sum.im);
}


static REAL8 psd_weighted_filter_inner_product(
	const COMPLEX16FrequencySeries *filter1,
	const COMPLEX16FrequencySeries *filter2,
	const REAL8Sequence *correlation,
	const REAL8FrequencySeries *psd
)
{
	const int k10 = floor(filter1->f0 / filter1->deltaF + 0.5);
	const int k20 = floor(filter2->f0 / filter2->deltaF + 0.5);
	const REAL8 *pdata = psd->data->data - (int) (psd->f0 / psd->deltaF);
	int k1, k2;
	COMPLEX16 sum = {0, 0};

	for(k1 = 0; k1 < (int) filter1->data->length; k1++) {
		const COMPLEX16 *f1data = &filter1->data->data[k1];
		for(k2 = 0; k2 < (int) filter2->data->length; k2++) {
			const COMPLEX16 *f2data = &filter2->data->data[k2];
			const unsigned delta_k = abs(k10 + k1 - k20 - k2);
			double sksk = (delta_k & 1 ? -1 : +1) * (delta_k < correlation->length ? correlation->data[delta_k] : 0);

			sksk *= sqrt(pdata[k10 + k1] * pdata[k20 + k2]);

			sum.re += sksk * (f1data->re * f2data->re + f1data->im * f2data->im);
			sum.im += sksk * (f1data->im * f2data->re - f1data->re * f2data->im);
		}
	}

	return 2 * sqrt(sum.re * sum.re + sum.im * sum.im);
}


/*
 * Generate the frequency domain channel filter function.  The filter is
 * nominally a Hann window twice the channel's width, centred on the
 * channel's centre frequency.  The filter is normalized so that its
 * "magnitude" as defined by the inner product function above is N.  Then
 * the filter is divided by the square root of the PSD frequency series
 * prior to normalilization.  This has the effect of de-emphasizing
 * frequency bins with high noise content, and is called "over whitening".
 */


static COMPLEX16FrequencySeries *generate_filter(
	REAL8 channel_flow,
	REAL8 channel_width,
	const REAL8FrequencySeries *psd,
	const REAL8Sequence *correlation
)
{
	static const char func[] = "generate_filter";
	char filter_name[100];
	REAL8Window *hann;
	COMPLEX16FrequencySeries *filter;
	REAL8 *pdata;
	unsigned i;
	REAL8 norm;

	sprintf(filter_name, "channel %g +/- %g Hz", channel_flow + channel_width / 2, channel_width / 2);

	/*
	 * Channel filter is a Hann window twice the channel's width,
	 * centred on the channel's centre frequency.  This makes a sum
	 * across channels equivalent to constructing a Tukey window
	 * spanning the same frequency band.  This trick is one of the
	 * ingredients that allows us to accomplish a multi-resolution
	 * tiling using a single frequency channel projection.  Really,
	 * there's no need for the "effective window" resulting from
	 * summing across channels to be something that has a name, any
	 * channel filter at all would do, but this way the code's
	 * behaviour is more easily understood --- it's easy to say "the
	 * channel filter is a Tukey window of variable central width".
	 *
	 * Note:  the number of samples in the window is odd, being one
	 * more than the number of frequency bins in twice the channel
	 * width.  This gets the Hann windows to super-impose to form a
	 * Tukey window.  (you'll have to draw yourself a picture).
	 */

	filter = XLALCreateCOMPLEX16FrequencySeries(filter_name, &psd->epoch, channel_flow - channel_width / 2, psd->deltaF, &lalDimensionlessUnit, 2 * channel_width / psd->deltaF + 1);
	/* FIXME:  decide what to do about this */
	hann = XLALCreateHannREAL8Window(filter->data->length);
	/*hann = XLALCreateTukeyREAL8Window(filter->data->length, 2 / ((channel_width / 2.0) + 1));*/
	if(!filter || !hann) {
		XLALDestroyCOMPLEX16FrequencySeries(filter);
		XLALDestroyREAL8Window(hann);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}
	for(i = 0; i < filter->data->length; i++) {
		filter->data->data[i].re = hann->data->data[i];
		filter->data->data[i].im = 0.0;
	}
	XLALDestroyREAL8Window(hann);

	/*
	 * divide by square root of PSD to "overwhiten".
	 */

	pdata = psd->data->data + (int) floor((filter->f0 - psd->f0) / psd->deltaF + 0.5);
	for(i = 0; i < filter->data->length; i++) {
		filter->data->data[i].re /= sqrt(pdata[i]);
		filter->data->data[i].im /= sqrt(pdata[i]);
	}

	/*
	 * normalize the filter.  the filter needs to be normalized so that
	 * it's inner product with itself is (width / delta F), the width
	 * of the filter in bins.
	 */

	norm = sqrt((channel_width / filter->deltaF) / filter_inner_product(filter, filter, correlation));
	for(i = 0; i < filter->data->length; i++) {
		filter->data->data[i].re *= norm;
		filter->data->data[i].im *= norm;
	}

	/*
	 * success
	 */

	return filter;
}


/**
 * From the power spectral density function, generate the comb of channel
 * filters for the time-frequency plane --- an excess power filter bank.
 */


LALExcessPowerFilterBank *XLALCreateExcessPowerFilterBank(
	double filter_deltaF,
	double flow,
	double channel_bandwidth,
	int n_channels,
	const REAL8FrequencySeries *psd,
	const REAL8Sequence *two_point_spectral_correlation
)
{
	static const char func[] = "XLALCreateExcessPowerFilterBank";
	LALExcessPowerFilterBank *new;
	struct ExcessPowerFilter *basis_filters;
	REAL8Sequence *twice_channel_overlap;
	REAL8Sequence *unwhitened_cross;
	int i;

	new = malloc(sizeof(*new));
	basis_filters = calloc(n_channels, sizeof(*basis_filters));
	twice_channel_overlap = XLALCreateREAL8Sequence(n_channels - 1);
	unwhitened_cross = XLALCreateREAL8Sequence(n_channels - 1);
	if(!new || !basis_filters || !twice_channel_overlap || !unwhitened_cross) {
		free(new);
		free(basis_filters);
		XLALDestroyREAL8Sequence(twice_channel_overlap);
		XLALDestroyREAL8Sequence(unwhitened_cross);
		XLAL_ERROR_NULL(func, XLAL_ENOMEM);
	}

	new->n_filters = n_channels;
	new->basis_filters = basis_filters;
	new->twice_channel_overlap = twice_channel_overlap;
	new->unwhitened_cross = unwhitened_cross;

	for(i = 0; i < n_channels; i++) {
		basis_filters[i].fseries = generate_filter(flow + i * channel_bandwidth, channel_bandwidth, psd, two_point_spectral_correlation);
		if(!basis_filters[i].fseries) {
			while(i--)
				XLALDestroyCOMPLEX16FrequencySeries(basis_filters[i].fseries);
			free(new);
			XLALDestroyREAL8Sequence(twice_channel_overlap);
			XLALDestroyREAL8Sequence(unwhitened_cross);
			XLAL_ERROR_NULL(func, XLAL_EFUNC);
		}

		/* compute the unwhitened root mean square for this channel */
		basis_filters[i].unwhitened_rms = sqrt(psd_weighted_filter_inner_product(basis_filters[i].fseries, basis_filters[i].fseries, two_point_spectral_correlation, psd) * filter_deltaF / 2);
	}

	/* compute the cross terms for the channel normalizations and
	 * unwhitened mean squares */
	for(i = 0; i < new->n_filters - 1; i++) {
		twice_channel_overlap->data[i] = 2 * filter_inner_product(basis_filters[i].fseries, basis_filters[i + 1].fseries, two_point_spectral_correlation);
		unwhitened_cross->data[i] = psd_weighted_filter_inner_product(basis_filters[i].fseries, basis_filters[i + 1].fseries, two_point_spectral_correlation, psd) * psd->deltaF;
	}

	return new;
}


/**
 * Destroy and excess power filter bank.
 */


void XLALDestroyExcessPowerFilterBank(
	LALExcessPowerFilterBank *bank
)
{
	if(bank) {
		if(bank->basis_filters) {
			int i;
			for(i = 0; i < bank->n_filters; i++)
				XLALDestroyCOMPLEX16FrequencySeries(bank->basis_filters[i].fseries);
			free(bank->basis_filters);
		}
		XLALDestroyREAL8Sequence(bank->twice_channel_overlap);
		XLALDestroyREAL8Sequence(bank->unwhitened_cross);
	}

	free(bank);
}
