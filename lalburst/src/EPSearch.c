/*
 *
 * Copyright (C) 2007  Brady, P. and Kipp Cannon
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


#include <complex.h>
#include <math.h>


#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>


#include <lal/Date.h>
#include <lal/EPSearch.h>
#include <lal/FrequencySeries.h>
#include <lal/LALChisq.h>
#include <lal/LALDatatypes.h>
#include <lal/LALMalloc.h>
#include <lal/LALStdlib.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataBurstUtils.h>
#include <lal/RealFFT.h>
#include <lal/TFTransform.h>
#include <lal/TimeFreqFFT.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/XLALError.h>
#include <lal/Window.h>


static double min(double a, double b)
{
	return a < b ? a : b;
}


static double max(double a, double b)
{
	return a > b ? a : b;
}


/*
 * ============================================================================
 *
 *                                Filter Bank
 *
 * ============================================================================
 */


typedef struct tagExcessPowerFilter {
	COMPLEX16FrequencySeries *fseries;
	REAL8 unwhitened_rms;			/**< root mean square of the unwhitened time series corresponding to this filter */
} ExcessPowerFilter;


typedef struct tagLALExcessPowerFilterBank {
	int n_filters;
	ExcessPowerFilter *basis_filters;
	REAL8Sequence *twice_channel_overlap;	/**< twice the inner product of filters for neighbouring channels;  twice_channel_overlap[0] is twice the inner product of the filters for channels 0 and 1, and so on (for n channels, there are n - 1 channel_overlaps) */
	REAL8Sequence *unwhitened_cross;	/**< the mean square cross terms for wide channels (indices same as for twice_channel_overlap) */
} LALExcessPowerFilterBank;


/**
 * From the power spectral density function, generate the comb of channel
 * filters for the time-frequency plane --- an excess power filter bank.
 */
static LALExcessPowerFilterBank *XLALCreateExcessPowerFilterBank(
	double filter_deltaF,
	double flow,
	double channel_bandwidth,
	int n_channels,
	const REAL8FrequencySeries *psd,
	const REAL8Sequence *two_point_spectral_correlation
)
{
	LALExcessPowerFilterBank *new;
	ExcessPowerFilter *basis_filters;
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
		XLAL_ERROR_NULL(XLAL_ENOMEM);
	}

	new->n_filters = n_channels;
	new->basis_filters = basis_filters;
	new->twice_channel_overlap = twice_channel_overlap;
	new->unwhitened_cross = unwhitened_cross;

	for(i = 0; i < n_channels; i++) {
		basis_filters[i].fseries = XLALCreateExcessPowerFilter(flow + i * channel_bandwidth, channel_bandwidth, psd, two_point_spectral_correlation);
		if(!basis_filters[i].fseries) {
			while(i--)
				XLALDestroyCOMPLEX16FrequencySeries(basis_filters[i].fseries);
			free(new);
			XLALDestroyREAL8Sequence(twice_channel_overlap);
			XLALDestroyREAL8Sequence(unwhitened_cross);
			XLAL_ERROR_NULL(XLAL_EFUNC);
		}

		/* compute the unwhitened root mean square for this channel */
		basis_filters[i].unwhitened_rms = sqrt(XLALExcessPowerFilterInnerProduct(basis_filters[i].fseries, basis_filters[i].fseries, two_point_spectral_correlation, psd) * filter_deltaF / 2);
	}

	/* compute the cross terms for the channel normalizations and
	 * unwhitened mean squares */
	for(i = 0; i < new->n_filters - 1; i++) {
		twice_channel_overlap->data[i] = 2 * XLALExcessPowerFilterInnerProduct(basis_filters[i].fseries, basis_filters[i + 1].fseries, two_point_spectral_correlation, NULL);
		unwhitened_cross->data[i] = XLALExcessPowerFilterInnerProduct(basis_filters[i].fseries, basis_filters[i + 1].fseries, two_point_spectral_correlation, psd) * psd->deltaF;
	}

	return new;
}


/**
 * Destroy and excess power filter bank.
 */
static void XLALDestroyExcessPowerFilterBank(
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


/*
 * ============================================================================
 *
 *                             Timing Arithmetic
 *
 * ============================================================================
 */


/**
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
	UINT4 segments;

	/*
	 * check input
	 */

	if(segment_length < 1) {
		XLALPrintError("segment_length < 1");
		XLAL_ERROR(XLAL_EINVAL);
	}
	if(segment_shift < 1) {
		XLALPrintError("segment_shift < 1");
		XLAL_ERROR(XLAL_EINVAL);
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


/**
 * Compute and return the timing parameters for an excess power analysis.
 * Pass NULL for any optional pointer to not compute and return that
 * parameter.  The return is 0 on success, negative on failure.
 */
int XLALEPGetTimingParameters(
	int window_length,	/**< Number of samples in a window used for the time-frequency plane */
	int max_tile_length,	/**< Number of samples in the tile of longest duration */
	double fractional_tile_shift,	/**< Number of samples by which the start of the longest tile is shifted from the start of the tile preceding it, as a fraction of its length */
	int *psd_length,	/**< (optional) User's desired number of samples to use in computing a PSD estimate.  Will be replaced with actual number of samples to use in computing a PSD estimate (rounded down to be comensurate with the windowing). */
	int *psd_shift,	/**< (optional) Number of samples by which the start of a PSD is to be shifted from the start of the PSD that preceded it in order that the tiling pattern continue smoothly across the boundary. */
	int *window_shift,	/**< Number of samples by which the start of a time-frequency plane window is shifted from the window preceding it in order that the tiling pattern continue smoothly across the boundary. */
	int *window_pad,	/**< How many samples at the start and end of each window are treated as padding, and will not be covered by the tiling. */
	int *tiling_length	/**< How many samples will be covered by the tiling. */
)
{
	int max_tile_shift = fractional_tile_shift * max_tile_length;

	/*
	 * check input parameters
	 */

	if(window_length % 4 != 0) {
		XLALPrintError("window_length is not a multiple of 4");
		XLAL_ERROR(XLAL_EINVAL);
	}
	if(max_tile_length < 1) {
		XLALPrintError("max_tile_length < 1");
		XLAL_ERROR(XLAL_EINVAL);
	}
	if(fractional_tile_shift <= 0) {
		XLALPrintError("fractional_tile_shift <= 0");
		XLAL_ERROR(XLAL_EINVAL);
	}
	if(fmod(fractional_tile_shift * max_tile_length, 1) != 0) {
		XLALPrintError("fractional_tile_shift * max_tile_length not an integer");
		XLAL_ERROR(XLAL_EINVAL);
	}
	if(max_tile_shift < 1) {
		XLALPrintError("fractional_tile_shift * max_tile_length < 1");
		XLAL_ERROR(XLAL_EINVAL);
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

	*window_pad = 4096;

	/*
	 * tiling covers the remainder, rounded down to fit an integer
	 * number of tiles
	 */

	*tiling_length = window_length - 2 * *window_pad;
	*tiling_length = XLALOverlappedSegmentsCommensurate(*tiling_length, max_tile_length, max_tile_shift);
	if(*tiling_length <= 0) {
		XLALPrintError("window_length too small for tiling, must be >= 2 * %d + %d", *window_pad, max_tile_length);
		XLAL_ERROR(XLAL_EINVAL);
	}

	/*
	 * now re-compute window_pad from rounded-off tiling_length
	 */

	*window_pad = (window_length - *tiling_length) / 2;
	if(*tiling_length + 2 * *window_pad != window_length) {
		XLALPrintError("window_length does not permit equal padding before and after tiling");
		XLAL_ERROR(XLAL_EINVAL);
	}

	/*
	 * adjacent tilings overlap so that their largest tiles overlap the
	 * same as within each tiling
	 */

	*window_shift = *tiling_length - (max_tile_length - max_tile_shift);
	if(*window_shift < 1) {
		XLALPrintError("window_shift < 1");
		XLAL_ERROR(XLAL_EINVAL);
	}

	/*
	 * compute the adjusted PSD length if desired
	 */

	if(psd_length) {
		*psd_length = XLALOverlappedSegmentsCommensurate(*psd_length, window_length, *window_shift);
		if(*psd_length < 0)
			XLAL_ERROR(XLAL_EFUNC);

		*psd_shift = *psd_length - (window_length - *window_shift);
		if(*psd_shift < 1) {
			XLALPrintError("psd_shift < 1");
			XLAL_ERROR(XLAL_EINVAL);
		}
	} else if(psd_shift) {
		/* for safety */
		*psd_shift = -1;
		/* can't compute psd_shift without psd_length input */
		XLAL_ERROR(XLAL_EFAULT);
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
 * A time-frequency plane
 */


typedef struct tagREAL8TimeFrequencyPlaneTiles {
	unsigned max_length;
	unsigned min_channels;
	unsigned max_channels;
	unsigned tiling_start;
	unsigned tiling_end;
	unsigned inv_fractional_stride;
	double dof_per_pixel;
} REAL8TimeFrequencyPlaneTiles;


typedef struct tagREAL8TimeFrequencyPlane {
	CHAR name[LALNameLength];	/**< name of data from which this was computed */
	LIGOTimeGPS epoch;		/**< epoch of data from which this was computed */
	REAL8 deltaT;			/**< time resolution of the plane */
	REAL8 fseries_deltaF;		/**< input frequency series' resolution */
	REAL8 deltaF;			/**< TF plane's frequency resolution (channel spacing) */
	REAL8 flow;			/**< low frequency boundary of TF plane */
	gsl_matrix *channel_data;   	/**< channel data.  each channel is placed into its own column.  channel_data[i * channels + j] corresponds to time epoch + i * deltaT and the frequency band [flow + j * deltaF, flow + (j + 1) * deltaF) */
	REAL8Sequence *channel_buffer;	/**< re-usable holding area for the data for a single channel */
	REAL8Sequence *unwhitened_channel_buffer;	/**< UNDOCUMENTED */
	REAL8TimeFrequencyPlaneTiles tiles;	/**< time-frequency plane's tiling information */
	REAL8Window *window;		/**< time-domain window applied to input time series for tapering edges to 0 */
	INT4 window_shift;		/**< by how many samples a window's start should be shifted from the start of the window preceding it */
	REAL8Sequence *two_point_spectral_correlation;	/**< two-point spectral correlation of the whitened frequency series, computed from the time-domain window function */
} REAL8TimeFrequencyPlane;


/**
 * Create and initialize a time-frequency plane object.
 */


static REAL8TimeFrequencyPlane *XLALCreateTFPlane(
	UINT4 tseries_length,		/**< length of time series from which TF plane will be computed */
	REAL8 tseries_deltaT,		/**< sample rate of time series */
	REAL8 flow,			/**< minimum frequency to search for */
	REAL8 bandwidth,		/**< bandwidth of TF plane */
	REAL8 tiling_fractional_stride,	/**< overlap of adjacent tiles */
	REAL8 max_tile_bandwidth,	/**< largest tile's bandwidth */
	REAL8 max_tile_duration,	/**< largest tile's duration */
	const REAL8FFTPlan *plan	/**< forward plan whose length is tseries_length */
)
{
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

	const int channels = round(bandwidth / deltaF);

	/*
	 * stride
	 */

	const unsigned inv_fractional_stride = round(1.0 / tiling_fractional_stride);

	/*
	 * tile size limits
	 */

	const unsigned min_length = round((1 / max_tile_bandwidth) / tseries_deltaT);
	const unsigned max_length = round(max_tile_duration / tseries_deltaT);
	const unsigned min_channels = inv_fractional_stride;
	const unsigned max_channels = round(max_tile_bandwidth / deltaF);

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
		XLAL_ERROR_NULL(XLAL_EFUNC);

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
		XLAL_ERROR_NULL(XLAL_EINVAL);
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
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	/*
	 * Initialize the structure
	 */

	plane->name[0] = '\0';
	XLALGPSSetREAL8(&plane->epoch, 0.0);
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


/**
 * Free a time-frequency plane object.
 */


static void XLALDestroyTFPlane(
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
 *                      Time-Frequency Plane Projection
 *
 * ============================================================================
 */


/*
 * Multiply the data by the filter.  The check that the frequency
 * resolutions and units are compatible is omitted because it is implied by
 * the calling code.
 */


static COMPLEX16Sequence *apply_filter(
	COMPLEX16Sequence *outputseq,
	const COMPLEX16FrequencySeries *inputseries,
	const COMPLEX16FrequencySeries *filterseries
)
{
	/* find bounds of common frequencies */
	const double flo = max(filterseries->f0, inputseries->f0);
	const double fhi = min(filterseries->f0 + filterseries->data->length * filterseries->deltaF, inputseries->f0 + inputseries->data->length * inputseries->deltaF);
	COMPLEX16 *output = outputseq->data + (int) round((flo - inputseries->f0) / inputseries->deltaF);
	COMPLEX16 *last = outputseq->data + (int) round((fhi - inputseries->f0) / inputseries->deltaF);
	const COMPLEX16 *input = inputseries->data->data + (int) round((flo - inputseries->f0) / inputseries->deltaF);
	const COMPLEX16 *filter = filterseries->data->data + (int) round((flo - filterseries->f0) / filterseries->deltaF);

	if(outputseq->length != inputseries->data->length)
		XLAL_ERROR_NULL(XLAL_EBADLEN);

	if(((unsigned) (output - outputseq->data) > outputseq->length) || (last - outputseq->data < 0))
		/* inputseries and filterseries don't intersect */
		memset(outputseq->data, 0, outputseq->length * sizeof(*outputseq->data));
	else {
		/* output = inputseries * conj(filter) */
		memset(outputseq->data, 0, (output - outputseq->data) * sizeof(*outputseq->data));
		for(; output < last; output++, input++, filter++)
			*output = *input * conj(*filter);
		memset(last, 0, (outputseq->length - (last - outputseq->data)) * sizeof(*outputseq->data));
	}

	return outputseq;
}


/*
 * Project a frequency series onto the comb of channel filters
 */


static int XLALFreqSeriesToTFPlane(
	REAL8TimeFrequencyPlane *plane,
	const LALExcessPowerFilterBank *filter_bank,
	const COMPLEX16FrequencySeries *fseries,
	const REAL8FFTPlan *reverseplan
)
{
	COMPLEX16Sequence *fcorr;
	unsigned i;

	/* check input parameters */
	if((fmod(plane->deltaF, fseries->deltaF) != 0.0) ||
	   (fmod(plane->flow - fseries->f0, fseries->deltaF) != 0.0))
		XLAL_ERROR(XLAL_EINVAL);

	/* make sure the frequency series spans an appropriate band */
	if((plane->flow < fseries->f0) ||
	   (plane->flow + plane->channel_data->size2 * plane->deltaF > fseries->f0 + fseries->data->length * fseries->deltaF))
		XLAL_ERROR(XLAL_EDATA);

	/* create temporary vectors */
	fcorr = XLALCreateCOMPLEX16Sequence(fseries->data->length);
	if(!fcorr)
		XLAL_ERROR(XLAL_EFUNC);

#if 0
	/* diagnostic code to dump data for the \hat{s}_{k} histogram */
	{
	unsigned k;
	FILE *f = fopen("sk.dat", "a");
	for(k = plane->flow / fseries->deltaF; k < (plane->flow + plane->channel_data->size2 * plane->deltaF) / fseries->deltaF; k++)
		fprintf(f, "%g\n%g\n", fseries->data->data[k].re, fseries->data->data[k].im);
	fclose(f);
	}
#endif
#if 0
	/* diagnostic code to dump data for the \hat{s}_{k}
	 * \hat{s}^{*}_{k'} histogram */
	{
	unsigned k, dk;
	FILE *f = fopen("sksk.dat", "a");
	for(dk = 0; dk < 100; dk++) {
		double avg_r = 0;
		double avg_i = 0;
	for(k = plane->flow / fseries->deltaF; k + dk < (plane->flow + plane->channel_data->size2 * plane->deltaF) / fseries->deltaF; k++) {
		double dr = fseries->data->data[k].re;
		double di = fseries->data->data[k].im;
		double dkr = fseries->data->data[k + dk].re;
		double dki = fseries->data->data[k + dk].im;
		avg_r += dr * dkr + di * dki;
		avg_i += di * dkr - dr * dki;
	}
		avg_r /= k - plane->flow / fseries->deltaF;
		avg_i /= k - plane->flow / fseries->deltaF;
		fprintf(f, "%d %g %g\n", dk, avg_r, avg_i);
	}
	fclose(f);
	}
#endif

	/* loop over the time-frequency plane's channels */
	for(i = 0; i < plane->channel_data->size2; i++) {
		unsigned j;
		/* cross correlate the input data against the channel
		 * filter by taking their product in the frequency domain
		 * and then inverse transforming to the time domain to
		 * obtain an SNR time series.  Note that
		 * XLALREAL8ReverseFFT() omits the factor of 1 / (N Delta
		 * t) in the inverse transform. */
		apply_filter(fcorr, fseries, filter_bank->basis_filters[i].fseries);
		if(XLALREAL8ReverseFFT(plane->channel_buffer, fcorr, reverseplan)) {
			XLALDestroyCOMPLEX16Sequence(fcorr);
			XLAL_ERROR(XLAL_EFUNC);
		}
		/* interleave the result into the channel_data array */
		for(j = 0; j < plane->channel_buffer->length; j++)
			gsl_matrix_set(plane->channel_data, j, i, plane->channel_buffer->data[j]);
	}

	/* clean up */
	XLALDestroyCOMPLEX16Sequence(fcorr);

	/* set the name and epoch of the TF plane */
	strncpy(plane->name, fseries->name, LALNameLength);
	plane->epoch = fseries->epoch;

	/* success */
	return 0;
}


/*
 * ============================================================================
 *
 *                               Tile Analysis
 *
 * ============================================================================
 */


static double compute_unwhitened_mean_square(
	const LALExcessPowerFilterBank *filter_bank,
	unsigned channel,
	unsigned channels
)
{
	unsigned i;
	double mean_square = 0;

	for(i = channel; i < channel + channels; i++)
		mean_square += pow(filter_bank->basis_filters[i].unwhitened_rms, 2);

	return mean_square;
}


/*
 * Convert time-frequency tile info to a SnglBurst row.
 */


static SnglBurst *XLALTFTileToBurstEvent(
	const REAL8TimeFrequencyPlane *plane,
	double tile_start,	/* in samples, allowed to be non-integer */
	double tile_length,	/* in samples, allowed to be non-integer */
	double f_centre,
	double bandwidth,
	double h_rss,
	double E,
	double d,
	double confidence
)
{
	SnglBurst *event = XLALCreateSnglBurst();

	if(!event)
		XLAL_ERROR_NULL(XLAL_ENOMEM);

	event->next = NULL;
	strncpy(event->ifo, plane->name, 2);
	event->ifo[2] = '\0';
	strncpy(event->search, "excesspower", LIGOMETA_SEARCH_MAX);
	event->search[LIGOMETA_SEARCH_MAX - 1] = '\0';
	strncpy(event->channel, plane->name, LIGOMETA_CHANNEL_MAX);
	event->channel[LIGOMETA_CHANNEL_MAX - 1] = '\0';
	event->start_time = plane->epoch;
	XLALGPSAdd(&event->start_time, tile_start * plane->deltaT);
	event->duration = tile_length * plane->deltaT;
	event->peak_time = event->start_time;
	XLALGPSAdd(&event->peak_time, event->duration / 2);
	event->bandwidth = bandwidth;
	event->central_freq = f_centre;
	/* FIXME: put h_rss into the "hrss" column */
	event->amplitude = h_rss;
	event->snr = E / d - 1;
	/* -ln P(event | stationary Gaussian white noise) */
	event->confidence = confidence;

	return event;
}


static SnglBurst *XLALComputeExcessPower(
	const REAL8TimeFrequencyPlane *plane,
	const LALExcessPowerFilterBank *filter_bank,
	SnglBurst *head,
	double confidence_threshold
)
{
	gsl_vector filter_output = {
		.size = plane->tiles.tiling_end - plane->tiles.tiling_start,
		.stride = plane->channel_data->size2,
		.data = NULL,
		.block = NULL,
		.owner = 0
	};
	gsl_vector_view filter_output_view;
	gsl_vector *channel_buffer;
	gsl_vector *unwhitened_channel_buffer;
	unsigned channel;
	unsigned channels;
	unsigned channel_end;
	double h_rss;
	double confidence;
	/* number of degrees of freedom in tile = number of
	 * "virtual pixels" in tile. */
	double tile_dof;
	unsigned i;

	/*
	 * create work spaces
	 */

	channel_buffer = gsl_vector_alloc(filter_output.size);
	unwhitened_channel_buffer = gsl_vector_alloc(filter_output.size);
	if(!channel_buffer || !unwhitened_channel_buffer) {
		if(channel_buffer)
			gsl_vector_free(channel_buffer);
		if(unwhitened_channel_buffer)
			gsl_vector_free(unwhitened_channel_buffer);
		XLAL_ERROR_NULL(XLAL_ENOMEM);
	}

	/*
	 * enter loops.  note:  this is a quadruply-nested loop that is not
	 * indented to help fit the code on a terminal display.
	 */

	for(channels = plane->tiles.min_channels; channels <= plane->tiles.max_channels; channels *= 2) {
		/* compute distance between "virtual pixels" for this
		 * (wide) channel */
		const unsigned stride = round(1.0 / (channels * plane->tiles.dof_per_pixel));

	for(channel_end = (channel = 0) + channels; channel_end <= plane->channel_data->size2; channel_end = (channel += channels / plane->tiles.inv_fractional_stride) + channels) {
		/* the root mean square of the "virtual channel",
		 * \sqrt{\mu^{2}} in the algorithm description */
		const double sample_rms = sqrt(channels * plane->deltaF / plane->fseries_deltaF + XLALREAL8SequenceSum(filter_bank->twice_channel_overlap, channel, channels - 1));
		/* the root mean square of the "uwapprox" quantity computed
		 * below, which is proportional to an approximation of the
		 * unwhitened time series. */
		double uwsample_rms;
		/* true unwhitened root mean square for this channel.  the
		 * ratio of this squared to uwsample_rms^2 is the
		 * correction factor to be applied to uwapprox^2 to convert
		 * it to an approximation of the square of the unwhitened
		 * channel */
		const double strain_rms = sqrt(compute_unwhitened_mean_square(filter_bank, channel, channels) + XLALREAL8SequenceSum(filter_bank->unwhitened_cross, channel, channels - 1));

		/* compute uwsample_rms */
		uwsample_rms = compute_unwhitened_mean_square(filter_bank, channel, channels);
		for(i = channel; i < channel_end - 1; i++)
			uwsample_rms += filter_bank->twice_channel_overlap->data[i] * filter_bank->basis_filters[i].unwhitened_rms * filter_bank->basis_filters[i + 1].unwhitened_rms * plane->fseries_deltaF / plane->deltaF;
		uwsample_rms = sqrt(uwsample_rms);

		/* reconstruct the time series and unwhitened time series
		 * for this (possibly multi-filter) channel.  both time
		 * series are normalized so that each sample has a mean
		 * square of 1 */
		filter_output.data = plane->channel_data->data + filter_output.stride * plane->tiles.tiling_start + channel;
		filter_output_view = gsl_vector_subvector_with_stride(&filter_output, 0, stride, filter_output.size / stride);
		gsl_vector_set_zero(channel_buffer);
		gsl_vector_set_zero(unwhitened_channel_buffer);
		channel_buffer->size = unwhitened_channel_buffer->size = filter_output_view.vector.size;
		for(i = channel; i < channel_end; filter_output_view.vector.data++, i++) {
			gsl_blas_daxpy(1.0 / sample_rms, &filter_output_view.vector, channel_buffer);
			gsl_blas_daxpy(filter_bank->basis_filters[i].unwhitened_rms * sqrt(plane->fseries_deltaF / plane->deltaF) / uwsample_rms, &filter_output_view.vector, unwhitened_channel_buffer);
		}

#if 0
		/* diagnostic code to dump data for the s_{j} histogram */
		{
		FILE *f = fopen("sj.dat", "a");
		for(i = 0; i < channel_buffer->size; i++)
			fprintf(f, "%g\n", gsl_vector_get(unwhitened_channel_buffer, i));
		fclose(f);
		}
#endif

		/* square the samples in the channel time series because
		 * from now on that's all we'll need */
		for(i = 0; i < channel_buffer->size; i++) {
			gsl_vector_set(channel_buffer, i, pow(gsl_vector_get(channel_buffer, i), 2));
			gsl_vector_set(unwhitened_channel_buffer, i, pow(gsl_vector_get(unwhitened_channel_buffer, i), 2));
		}

	/* start with at least 2 degrees of freedom */
	for(tile_dof = 2; tile_dof <= plane->tiles.max_length / stride; tile_dof *= 2) {
		unsigned start;
	for(start = 0; start + tile_dof <= channel_buffer->size; start += tile_dof / plane->tiles.inv_fractional_stride) {
		/* compute sum of squares, and unwhitened sum of squares
		 * (samples have already been squared) */
		double sumsquares = 0;
		double uwsumsquares = 0;
		for(i = start; i < start + tile_dof; i++) {
			sumsquares += gsl_vector_get(channel_buffer, i);
			uwsumsquares += gsl_vector_get(unwhitened_channel_buffer, i);
		}

		/* compute statistical confidence */
		/* FIXME:  the 0.62 is an empirically determined
		 * degree-of-freedom fudge factor.  figure out what its
		 * origin is, and account for it correctly.  it's most
		 * likely due to the time-frequency plane pixels not being
		 * independent of one another as a consequence of a
		 * non-zero inner product of the time-domain impulse
		 * response of the channel filter for adjacent pixels */
		confidence = -XLALLogChisqCCDF(sumsquares * .62, tile_dof * .62);
		if(XLALIsREAL8FailNaN(confidence)) {
			gsl_vector_free(channel_buffer);
			gsl_vector_free(unwhitened_channel_buffer);
			XLAL_ERROR_NULL(XLAL_EFUNC);
		}

		/* record tiles whose statistical confidence is above
		 * threshold and that have real-valued h_rss */
		if((confidence >= confidence_threshold) && (uwsumsquares >= tile_dof)) {
			SnglBurst *oldhead = head;

			/* compute h_rss */
			h_rss = sqrt((uwsumsquares - tile_dof) * (stride * plane->deltaT)) * strain_rms;

			/* add new event to head of linked list */
			head = XLALTFTileToBurstEvent(plane, plane->tiles.tiling_start + (start - 0.5) * stride, tile_dof * stride, plane->flow + (channel + .5 * channels) * plane->deltaF, channels * plane->deltaF, h_rss, sumsquares, tile_dof, confidence);
			if(!head) {
				gsl_vector_free(channel_buffer);
				gsl_vector_free(unwhitened_channel_buffer);
				XLAL_ERROR_NULL(XLAL_EFUNC);
			}
			head->next = oldhead;
		}
	}
	}
	}
	}

	/* success */
	gsl_vector_free(channel_buffer);
	gsl_vector_free(unwhitened_channel_buffer);
	return head;
}


/*
 * ============================================================================
 *
 *                                Entry Point
 *
 * ============================================================================
 */


/**
 * Generate a linked list of burst events from a time series.
 */
SnglBurst *XLALEPSearch(
	struct XLALEPSearchDiagnostics *diagnostics,
	const REAL8TimeSeries *tseries,
	REAL8Window *window,
	double flow,
	double bandwidth,
	double confidence_threshold,
	double fractional_stride,
	double maxTileBandwidth,
	double maxTileDuration
)
{
	SnglBurst *head = NULL;
	int errorcode = 0;
	int start_sample;
	COMPLEX16FrequencySeries *fseries;
	REAL8FFTPlan *fplan;
	REAL8FFTPlan *rplan;
	REAL8FrequencySeries *psd;
	REAL8TimeSeries *cuttseries = NULL;
	LALExcessPowerFilterBank *filter_bank = NULL;
	REAL8TimeFrequencyPlane *plane = NULL;

	/*
	 * Construct forward and reverse FFT plans, storage for the PSD,
	 * the time-frequency plane, and a tiling.  Note that the flat part
	 * of the Tukey window needs to match the locations of the tiles as
	 * specified by the tiling_start parameter of XLALCreateTFPlane.
	 * The metadata for the two frequency series will be filled in
	 * later, so it doesn't all have to be correct here.
	 */

	fplan = XLALCreateForwardREAL8FFTPlan(window->data->length, 1);
	rplan = XLALCreateReverseREAL8FFTPlan(window->data->length, 1);
	psd = XLALCreateREAL8FrequencySeries("PSD", &tseries->epoch, 0, 0, &lalDimensionlessUnit, window->data->length / 2 + 1);
	fseries = XLALCreateCOMPLEX16FrequencySeries(tseries->name, &tseries->epoch, 0, 0, &lalDimensionlessUnit, window->data->length / 2 + 1);
	if(fplan)
		plane = XLALCreateTFPlane(window->data->length, tseries->deltaT, flow, bandwidth, fractional_stride, maxTileBandwidth, maxTileDuration, fplan);
	if(!fplan || !rplan || !psd || !fseries || !plane) {
		errorcode = XLAL_EFUNC;
		goto error;
	}

#if 0
	/* diagnostic code to replace the input time series with stationary
	 * Gaussian white noise.  the normalization is such that it yields
	 * unit variance frequency components without a call to the
	 * whitening function. */
	{
	unsigned i;
	static RandomParams *rparams = NULL;
	if(!rparams)
		rparams = XLALCreateRandomParams(0);
	XLALNormalDeviates(tseries->data, rparams);
	for(i = 0; i < tseries->data->length; i++)
		tseries->data->data[i] *= sqrt(0.5 / tseries->deltaT);
	}
#endif
#if 0
	/* diagnostic code to disable the tapering window */
	{
	unsigned i = plane->window->data->length;
	XLALDestroyREAL8Window(plane->window);
	plane->window = XLALCreateRectangularREAL8Window(i);
	}
#endif

	/*
	 * Compute the average spectrum.
	 *
	 * FIXME: is using windowShift here correct?  we have to, otherwise
	 * the time series' lengths are inconsistent
	 */

	if(XLALREAL8AverageSpectrumMedian(psd, tseries, plane->window->data->length, plane->window_shift, plane->window, fplan) < 0) {
		errorcode = XLAL_EFUNC;
		goto error;
	}

	if(diagnostics)
		diagnostics->XLALWriteLIGOLwXMLArrayREAL8FrequencySeries(diagnostics->LIGOLwXMLStream, NULL, psd);

	/*
	 * Construct the time-frequency plane's channel filters.
	 */

	XLALPrintInfo("%s(): constructing channel filters\n", __func__);
	filter_bank = XLALCreateExcessPowerFilterBank(psd->deltaF, plane->flow, plane->deltaF, plane->channel_data->size2, psd, plane->two_point_spectral_correlation);
	if(!filter_bank) {
		errorcode = XLAL_EFUNC;
		goto error;
	}

	/*
	 * Loop over data applying excess power method.
	 */

	for(start_sample = 0; start_sample + plane->window->data->length <= tseries->data->length; start_sample += plane->window_shift) {
		/*
		 * Verbosity.
		 */

		XLALPrintInfo("%s(): ", __func__);
		XLALPrintProgressBar(start_sample / (double) (tseries->data->length - plane->window->data->length));
		XLALPrintInfo(" complete\n");

		/*
		 * Extract a window-length of data from the time series.
		 */

		cuttseries = XLALCutREAL8TimeSeries(tseries, start_sample, plane->window->data->length);
		if(!cuttseries) {
			errorcode = XLAL_EFUNC;
			goto error;
		}
		XLALPrintInfo("%s(): analyzing %u samples (%.9lf s) at offset %u (%.9lf s) from epoch %d.%09u s\n", __func__, cuttseries->data->length, cuttseries->data->length * cuttseries->deltaT, start_sample, start_sample * cuttseries->deltaT, tseries->epoch.gpsSeconds, tseries->epoch.gpsNanoSeconds);
		if(diagnostics)
			diagnostics->XLALWriteLIGOLwXMLArrayREAL8TimeSeries(diagnostics->LIGOLwXMLStream, NULL, cuttseries);

		/*
		 * Window and DFT the time series.
		 */

		XLALPrintInfo("%s(): computing the Fourier transform\n", __func__);
		if(!XLALUnitaryWindowREAL8Sequence(cuttseries->data, plane->window)) {
			errorcode = XLAL_EFUNC;
			goto error;
		}
		if(XLALREAL8TimeFreqFFT(fseries, cuttseries, fplan)) {
			errorcode = XLAL_EFUNC;
			goto error;
		}
		XLALDestroyREAL8TimeSeries(cuttseries);
		cuttseries = NULL;

		/*
		 * Normalize the frequency series to the average PSD.
		 */

#if 1
		XLALPrintInfo("%s(): normalizing to the average spectrum\n", __func__);
		if(!XLALWhitenCOMPLEX16FrequencySeries(fseries, psd)) {
			errorcode = XLAL_EFUNC;
			goto error;
		}
		if(diagnostics)
			diagnostics->XLALWriteLIGOLwXMLArrayCOMPLEX16FrequencySeries(diagnostics->LIGOLwXMLStream, "whitened", fseries);
#endif

		/*
		 * Compute the time-frequency plane from the frequency
		 * series and channel filters.
		 */

		XLALPrintInfo("%s(): projecting data onto time-frequency plane\n", __func__);
		if(XLALFreqSeriesToTFPlane(plane, filter_bank, fseries, rplan)) {
			errorcode = XLAL_EFUNC;
			goto error;
		}

		/*
		 * Compute the excess power for each time-frequency tile
		 * using the data in the time-frequency plane, and add
		 * those tiles whose confidence is above threshold to the
		 * trigger list.  Note that because it is possible for
		 * there to be 0 triggers found, we can't check for errors
		 * by testing for head == NULL.
		 */

		XLALPrintInfo("%s(): computing the excess power for each tile\n", __func__);
		XLALClearErrno();
		head = XLALComputeExcessPower(plane, filter_bank, head, confidence_threshold);
		if(xlalErrno) {
			errorcode = XLAL_EFUNC;
			goto error;
		}
	}

	/*
	 * Memory clean-up.
	 */

	XLALPrintInfo("%s(): done\n", __func__);

	error:
	XLALDestroyREAL8FFTPlan(fplan);
	XLALDestroyREAL8FFTPlan(rplan);
	XLALDestroyREAL8FrequencySeries(psd);
	XLALDestroyREAL8TimeSeries(cuttseries);
	XLALDestroyCOMPLEX16FrequencySeries(fseries);
	XLALDestroyExcessPowerFilterBank(filter_bank);
	XLALDestroyTFPlane(plane);
	if(errorcode) {
		XLALDestroySnglBurstTable(head);
		XLAL_ERROR_NULL(errorcode);
	}
	return(head);
}
