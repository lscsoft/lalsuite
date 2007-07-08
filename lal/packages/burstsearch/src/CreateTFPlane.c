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


#include <lal/LALRCSID.h>


NRCSID(CREATETFPLANEC, "$Id$");


#include <math.h>
#include <lal/FrequencySeries.h>
#include <lal/LALAtomicDatatypes.h>
#include <lal/LALMalloc.h>
#include <lal/LALStdlib.h>
#include <lal/RealFFT.h>
#include <lal/Sequence.h>
#include <lal/TFTransform.h>
#include <lal/XLALError.h>


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
 * Compute the two-point spectral correlation function for the whitened
 * frequency series from the window applied to the original time series.
 * The indices of the output sequence are |k - k'|.  The window is, by
 * construction, an even function of the sample index, so the Fourier
 * transform is real-valued (contains only cosine components).
 */


static REAL4Sequence *compute_two_point_spectral_correlation(
	const REAL4Window *window
)
{
	REAL4Sequence *w_squared = XLALCopyREAL4Sequence(window->data);
	COMPLEX8Sequence *tmp = XLALCreateCOMPLEX8Sequence(window->data->length / 2 + 1);
	REAL4Sequence *correlation = XLALCreateREAL4Sequence(window->data->length / 2 + 1);
	RealFFTPlan *plan = XLALCreateForwardREAL4FFTPlan(window->data->length, 0);
	unsigned i;

	if(!w_squared || !tmp || !correlation || !plan) {
		XLALDestroyREAL4Sequence(w_squared);
		XLALDestroyCOMPLEX8Sequence(tmp);
		XLALDestroyREAL4Sequence(correlation);
		XLALDestroyREAL4FFTPlan(plan);
		return NULL;
	}

	/* square and normalize the window */
	for(i = 0; i < w_squared->length; i++)
		w_squared->data[i] *= w_squared->data[i] / window->sumofsquares;

	/* Fourier transform */
	if(XLALREAL4ForwardFFT(tmp, w_squared, plan)) {
		XLALDestroyREAL4Sequence(correlation);
		correlation = NULL;
	} else {
		/* extract real components */
		for(i = 0; i < correlation->length; i++)
			correlation->data[i] = tmp->data[i].re;
	}

	XLALDestroyREAL4Sequence(w_squared);
	XLALDestroyCOMPLEX8Sequence(tmp);
	XLALDestroyREAL4FFTPlan(plan);

	return correlation;
}


/*
 * Create and initialize a time-frequency plane object.
 */


REAL4TimeFrequencyPlane *XLALCreateTFPlane(
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
	REAL8 max_tile_duration
)
{
	static const char func[] = "XLALCreateTFPlane";
	REAL4TimeFrequencyPlane *plane;
	COMPLEX8FrequencySeries **filter;
	REAL8Sequence *twice_channel_overlap;
	REAL8Sequence *unwhitened_rms;
	REAL8Sequence *unwhitened_cross;
	REAL4Sequence **channel;
	REAL4Window *tukey;
	REAL4Sequence *correlation;
	int i;

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

	const int channels = bandwidth / deltaF;

	/*
	 * stride
	 */

	const unsigned inv_fractional_stride = 1 / tiling_fractional_stride;

	/*
	 * tile size limits
	 */

	const unsigned min_length = (1 / max_tile_bandwidth) / tseries_deltaT;
	const unsigned max_length = max_tile_duration / tseries_deltaT;
	const unsigned min_channels = inv_fractional_stride;
	const unsigned max_channels = max_tile_bandwidth / deltaF;

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
	filter = XLALMalloc(channels * sizeof(*filter));
	twice_channel_overlap = XLALCreateREAL8Sequence(channels - 1);
	unwhitened_rms = XLALCreateREAL8Sequence(channels);
	unwhitened_cross = XLALCreateREAL8Sequence(channels - 1);
	channel = XLALMalloc(channels * sizeof(*channel));
	tukey = XLALCreateTukeyREAL4Window(tseries_length, (tseries_length - tiling_length) / (double) tseries_length);
	correlation = tukey ? compute_two_point_spectral_correlation(tukey) : NULL;
	if(!plane || !filter || !twice_channel_overlap || !unwhitened_rms || !unwhitened_cross || !channel || !tukey || !correlation) {
		XLALFree(plane);
		XLALFree(filter);
		XLALDestroyREAL8Sequence(twice_channel_overlap);
		XLALDestroyREAL8Sequence(unwhitened_rms);
		XLALDestroyREAL8Sequence(unwhitened_cross);
		XLALFree(channel);
		XLALDestroyREAL4Window(tukey);
		XLALDestroyREAL4Sequence(correlation);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}
	for(i = 0; i < channels; i++) {
		filter[i] = NULL;
		channel[i] = XLALCreateREAL4Sequence(tseries_length);
		if(!channel[i]) {
			while(--i)
				XLALDestroyREAL4Sequence(channel[i]);
			XLALFree(plane);
			XLALFree(filter);
			XLALDestroyREAL8Sequence(twice_channel_overlap);
			XLALDestroyREAL8Sequence(unwhitened_rms);
			XLALDestroyREAL8Sequence(unwhitened_cross);
			XLALFree(channel);
			XLALDestroyREAL4Window(tukey);
			XLALDestroyREAL4Sequence(correlation);
			XLAL_ERROR_NULL(func, XLAL_ENOMEM);
		}
	}

	/* 
	 * Initialize the structure
	 */

	plane->name[0] = '\0';
	plane->epoch.gpsSeconds = 0;
	plane->epoch.gpsNanoSeconds = 0;
	plane->deltaT = tseries_deltaT;
	plane->fseries_deltaF = fseries_deltaF;
	plane->channels = channels;
	plane->deltaF = deltaF;
	plane->flow = flow;
	plane->filter = filter;
	plane->twice_channel_overlap = twice_channel_overlap;
	plane->unwhitened_rms = unwhitened_rms;
	plane->unwhitened_cross = unwhitened_cross;
	plane->channel = channel;
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
	REAL4TimeFrequencyPlane *plane
)
{
	if(plane) {
		unsigned i;
		XLALDestroyREAL8Sequence(plane->twice_channel_overlap);
		XLALDestroyREAL8Sequence(plane->unwhitened_rms);
		XLALDestroyREAL8Sequence(plane->unwhitened_cross);
		for(i = 0; i < plane->channels; i++) {
			XLALDestroyCOMPLEX8FrequencySeries(plane->filter[i]);
			XLALDestroyREAL4Sequence(plane->channel[i]);
		}
		XLALFree(plane->filter);
		XLALFree(plane->channel);
		XLALDestroyREAL4Window(plane->window);
		XLALDestroyREAL4Sequence(plane->two_point_spectral_correlation);
	}
	XLALFree(plane);
}
