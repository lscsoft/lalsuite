/*
 * $Id$
 *
 * Copyright (C) 2007  Kipp Cannon and Flanagan, E.
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


#include <math.h>


#include <lal/LALRCSID.h>


NRCSID(FREQSERIESTOTFPLANEC, "$Id$");


#include <lal/Date.h>
#include <lal/LALDatatypes.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/RealFFT.h>
#include <lal/Sequence.h>
#include <lal/TFTransform.h>
#include <lal/Thresholds.h>
#include <lal/XLALError.h>


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


static double min(double a, double b)
{
	return a < b ? a : b;
}


static double max(double a, double b)
{
	return a > b ? a : b;
}


static COMPLEX16Sequence *apply_filter(
	COMPLEX16Sequence *outputseq,
	const COMPLEX16FrequencySeries *inputseries,
	const COMPLEX16FrequencySeries *filterseries
)
{
	static const char func[] = "apply_filter";
	/* find bounds of common frequencies */
	const double flo = max(filterseries->f0, inputseries->f0);
	const double fhi = min(filterseries->f0 + filterseries->data->length * filterseries->deltaF, inputseries->f0 + inputseries->data->length * inputseries->deltaF);
	COMPLEX16 *output = outputseq->data + (int) floor((flo - inputseries->f0) / inputseries->deltaF + 0.5);
	COMPLEX16 *last = outputseq->data + (int) floor((fhi - inputseries->f0) / inputseries->deltaF + 0.5);
	const COMPLEX16 *input = inputseries->data->data + (int) floor((flo - inputseries->f0) / inputseries->deltaF + 0.5);
	const COMPLEX16 *filter = filterseries->data->data + (int) floor((flo - filterseries->f0) / filterseries->deltaF + 0.5);

	if(outputseq->length != inputseries->data->length)
		XLAL_ERROR_NULL(func, XLAL_EBADLEN);

	/* zero the product vector */
	memset(outputseq->data, 0, outputseq->length * sizeof(*outputseq->data));

	/* output = inputseries * conj(filter) */
	for(; output < last; output++, input++, filter++) {
		output->re = input->re * filter->re + input->im * filter->im;
		output->im = input->im * filter->re - input->re * filter->im;
	}

	return outputseq;
}


/*
 * Project a frequency series onto the comb of channel filters
 */


/******** <lalVerbatim file="FreqSeriesToTFPlaneCP"> ********/
int XLALFreqSeriesToTFPlane(
	REAL8TimeFrequencyPlane *plane,
	const COMPLEX16FrequencySeries *fseries,
	const REAL8FFTPlan *reverseplan
)
/******** </lalVerbatim> ********/
{
	static const char func[] = "XLALFreqSeriesToTFPlane";
	COMPLEX16Sequence *fcorr;
	unsigned i;

	/* check input parameters */
	if((fmod(plane->deltaF, fseries->deltaF) != 0.0) ||
	   (fmod(plane->flow - fseries->f0, fseries->deltaF) != 0.0))
		XLAL_ERROR(func, XLAL_EINVAL);

	/* make sure the frequency series spans an appropriate band */
	if((plane->flow < fseries->f0) ||
	   (plane->flow + plane->channels * plane->deltaF > fseries->f0 + fseries->data->length * fseries->deltaF))
		XLAL_ERROR(func, XLAL_EDATA);

	/* create temporary vectors */
	fcorr = XLALCreateCOMPLEX16Sequence(fseries->data->length);
	if(!fcorr)
		XLAL_ERROR(func, XLAL_EFUNC);

#if 0
	/* diagnostic code to dump data for the \hat{s}_{k} histogram */
	{
	unsigned k;
	FILE *f = fopen("sk.dat", "a");
	for(k = plane->flow / fseries->deltaF; k < (plane->flow + plane->channels * plane->deltaF) / fseries->deltaF; k++)
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
	for(k = plane->flow / fseries->deltaF; k + dk < (plane->flow + plane->channels * plane->deltaF) / fseries->deltaF; k++) {
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
	for(i = 0; i < plane->channels; i++) {
		/* cross correlate the input data against the channel
		 * filter by taking their product in the frequency domain
		 * and then inverse transforming to the time domain to
		 * obtain an SNR time series.  Note that
		 * XLALREAL4ReverseFFT() omits the factor of 1 / (N Delta
		 * t) in the inverse transform. */
		apply_filter(fcorr, fseries, plane->filter[i]);
		if(XLALREAL8ReverseFFT(plane->channel[i], fcorr, reverseplan)) {
			XLALDestroyCOMPLEX16Sequence(fcorr);
			XLAL_ERROR(func, XLAL_EFUNC);
		}
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
 *                                   Output
 *
 * ============================================================================
 */


/*
 * Convert time-frequency tile info to a SnglBurstTable row.
 */


static SnglBurstTable *XLALTFTileToBurstEvent(
	const REAL8TimeFrequencyPlane *plane,
	unsigned tile_start,
	unsigned tile_length,
	double f_centre,
	double bandwidth,
	double h_rss,
	double E,
	double d,
	double confidence
)
{
	static const char func[] = "XLALTFTileToBurstEvent";
	SnglBurstTable *event = XLALCalloc(1, sizeof(*event));

	if(!event)
		XLAL_ERROR_NULL(func, XLAL_ENOMEM);

	event->next = NULL;
	strncpy(event->ifo, plane->name, 2);
	event->ifo[2] = '\0';
	strncpy(event->search, "excesspower", LIGOMETA_SEARCH_MAX);
	event->search[LIGOMETA_SEARCH_MAX] = '\0';
	strncpy(event->channel, plane->name, LIGOMETA_CHANNEL_MAX);
	event->channel[LIGOMETA_CHANNEL_MAX] = '\0';
	event->start_time = plane->epoch; 
	XLALGPSAdd(&event->start_time, tile_start * plane->deltaT);
	event->duration = tile_length * plane->deltaT;
	event->peak_time = event->start_time;
	XLALGPSAdd(&event->peak_time, 0.5 * event->duration);
	event->bandwidth = bandwidth;
	event->central_freq = f_centre;
	/* FIXME: put h_rss into the "hrss" column */
	event->amplitude = h_rss;
	event->snr = E / d - 1;
	/* -ln P(event | stationary Gaussian white noise) */
	event->confidence = confidence;
	/* for safety */
	event->string_cluster_t = XLAL_REAL4_FAIL_NAN;
	/* will be set correctly later */
	event->event_id = 0;

	return(event);
}


/*
 * ============================================================================
 *
 *                               Tile Analysis
 *
 * ============================================================================
 */


/******** <lalVerbatim file="ComputeExcessPowerCP"> ********/
SnglBurstTable *XLALComputeExcessPower(
	const REAL8TimeFrequencyPlane *plane,
	SnglBurstTable *head,
	double confidence_threshold
)
/******** </lalVerbatim> ********/
{
	static const char func[] = "XLALComputeExcessPower";
	unsigned start;
	unsigned length;
	unsigned end;
	unsigned channel;
	unsigned channels;
	unsigned channel_end;
	double h_rss;
	double confidence;

	for(channels = plane->tiles.min_channels; channels <= plane->tiles.max_channels; channels *= 2) {
	for(channel_end = (channel = 0) + channels; channel_end <= plane->channels; channel_end = (channel += channels / plane->tiles.inv_fractional_stride) + channels) {
		/* the mean square of the "virtual channel", \mu^{2} in the
		 * algorithm description */
		const double sample_mean_square = channels * plane->deltaF / plane->fseries_deltaF + XLALREAL8SequenceSum(plane->twice_channel_overlap, channel, channels - 1);
		/* the mean square of the "uwapprox" quantity computed
		 * below, which is proportional to an approximation of the
		 * unwhitened time series. */
		double uwapprox_mean_square;
		/* true unwhitened mean square for this channel.  the ratio
		 * of this to uwapprox_mean_square is the correction factor to
		 * be applied to uwapprox^{2} to convert it to an
		 * approximation of the square of the unwhitened channel */
		const double strain_mean_square = XLALREAL8SequenceSumSquares(plane->unwhitened_rms, channel, channels) + XLALREAL8SequenceSum(plane->unwhitened_cross, channel, channels - 1);
		unsigned c;

		/* compute uwapprox_mean_square */
		uwapprox_mean_square = XLALREAL8SequenceSumSquares(plane->unwhitened_rms, channel, channels);
		for(c = channel; c < channel_end - 1; c++)
			uwapprox_mean_square += plane->twice_channel_overlap->data[c] * plane->unwhitened_rms->data[c] * plane->unwhitened_rms->data[c + 1] * plane->fseries_deltaF / plane->deltaF;

	/* start with at least 2 degrees of freedom */
	for(length = 2 / (channels * plane->tiles.dof_per_pixel); length <= plane->tiles.max_length; length *= 2) {
		const double tile_dof = (length * channels) * plane->tiles.dof_per_pixel;
		/* number of degrees of freedom in tile = number of
		 * "virtual pixels" in tile.  compute distance between
		 * tile's "virtual pixels" */
		const unsigned stride = length / tile_dof;

	for(end = (start = plane->tiles.tiling_start) + length; end <= plane->tiles.tiling_end; end = (start += length / plane->tiles.inv_fractional_stride) + length) {
		double sumsquares = 0;
		double uwsumsquares = 0;
		unsigned t;

		/* compute sum of squares, and unwhitened sum of squares */
		for(t = start + stride / 2; t < end; t += stride) {
			double sample = 0;
			double uwapprox = 0;

			/* compute the whitened and (approximate)
			 * unwhitened time series samples for this t */
			for(c = channel; c < channel_end; c++) {
				sample += plane->channel[c]->data[t];
				uwapprox += plane->channel[c]->data[t] * plane->unwhitened_rms->data[c] * sqrt(plane->fseries_deltaF / plane->deltaF);
			}

#if 0
			/* diagnostic code to dump data for the s_{j} histogram */
			{
			FILE *f = fopen("sj.dat", "a");
			fprintf(f, "%g\n", uwapprox / sqrt(uwapprox_mean_square));
			fclose(f);
			}
#endif

			sumsquares += sample * sample;
			uwsumsquares += uwapprox * uwapprox;
		}
		/* normalization to give each sample a mean square of 1.0 */
		sumsquares /= sample_mean_square;
		uwsumsquares /= uwapprox_mean_square;

		/* compute statistical confidence */
		/* FIXME:  the 0.62 is an empirically determined
		 * degree-of-freedom fudge factor.  figure out what its
		 * origin is, and account for it correctly.  it's most
		 * likely due to the time-frequency plane pixels not being
		 * independent of one another as a consequence of a
		 * non-zero inner product of the time-domain impulse
		 * response of the channel filter for adjacent pixels */
		confidence = -XLALlnOneMinusChisqCdf(sumsquares * .62, tile_dof * .62);
		if(XLALIsREAL8FailNaN(confidence))
			XLAL_ERROR_NULL(func, XLAL_EFUNC);

		/* record tiles whose statistical confidence is above
		 * threshold and that have real-valued h_rss */
		if((confidence >= confidence_threshold) && (uwsumsquares >= tile_dof)) {
			SnglBurstTable *oldhead = head;

			/* compute h_rss */
			h_rss = sqrt((uwsumsquares - tile_dof) * strain_mean_square * stride * plane->deltaT);

			/* add new event to head of linked list */
			head = XLALTFTileToBurstEvent(plane, start, length, plane->flow + (channel + .5 * channels) * plane->deltaF, channels * plane->deltaF, h_rss, sumsquares, tile_dof, confidence);
			if(!head)
				XLAL_ERROR_NULL(func, XLAL_EFUNC);
			head->next = oldhead;
		}
	}
	}
	}
	}

	/* success */
	return(head);
}


/*
 * ============================================================================
 *
 *                              Experimentation
 *
 * ============================================================================
 */


static void heterodyne_channel(
	COMPLEX16Sequence *fdata,
	unsigned k0,
	unsigned bandwidth
)
{
	unsigned half_width = (bandwidth + 1) / 2;
	unsigned k_centre = k0 + half_width - 1;
	unsigned k;

	for(k = 0; k < half_width; k++) {
		fdata->data[k_centre + k].re += fdata->data[k_centre - k].re;
		fdata->data[k_centre + k].im -= fdata->data[k_centre - k].im;
	}

	memmove(fdata->data, &fdata->data[k_centre], half_width * sizeof(*fdata->data));
	memset(&fdata->data[half_width], 0, (fdata->length - half_width) * sizeof(*fdata->data));
	/* floating point rules allow (x-x) to equal -0.0 rather than +0.0,
	 * but some tests check that the DC component's imaginary part is
	 * +0.0, and will fail if it is -0.0 */
	fdata->data[0].im = 0;
}


SnglBurstTable *XLALExcessPowerProject(
	const COMPLEX16FrequencySeries *fseries,
	REAL8TimeFrequencyPlane *plane,
	const struct ExcessPowerTemplateBank *bank,
	SnglBurstTable *head,
	double confidence_threshold,
	const REAL8FFTPlan *reverseplan
)
{
	static const char func[] = "XLALExcessPowerProject";
	COMPLEX16Sequence *fcorr;
	REAL8Sequence *channel;
	unsigned start;
	unsigned length;
	unsigned end;
	int template;
	double h_rss;
	double confidence;

	/* check input parameters */
	if((fmod(plane->deltaF, fseries->deltaF) != 0.0) ||
	   (fmod(plane->flow - fseries->f0, fseries->deltaF) != 0.0))
		XLAL_ERROR_NULL(func, XLAL_EINVAL);

	/* make sure the frequency series spans an appropriate band */
	if((plane->flow < fseries->f0) ||
	   (plane->flow + plane->channels * plane->deltaF > fseries->f0 + fseries->data->length * fseries->deltaF))
		XLAL_ERROR_NULL(func, XLAL_EDATA);

	/* create temporary vectors */
	fcorr = XLALCreateCOMPLEX16Sequence(fseries->data->length);
	channel = XLALCreateREAL8Sequence(plane->window->data->length);
	if(!fcorr || !channel) {
		XLALDestroyCOMPLEX16Sequence(fcorr);
		XLALDestroyREAL8Sequence(channel);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}

	/* set the name and epoch of the TF plane for tile construction */
	strncpy(plane->name, fseries->name, LALNameLength);
	plane->epoch = fseries->epoch;

	/* loop over time-frequency tiles */
	for(template = 0; template < bank->n_templates; template++) {
		const double sum_mean_square = bank->templates[template].bandwidth / plane->fseries_deltaF;

		/* cross correlate the input data against the channel
		 * filter by taking their product in the frequency domain
		 * and then inverse transforming to the time domain to
		 * obtain an SNR time series.  Note that
		 * XLALREAL4ReverseFFT() omits the factor of 1 / (N Delta
		 * t) in the inverse transform. */
		apply_filter(fcorr, fseries, bank->templates[template].filter);
		heterodyne_channel(fcorr, bank->templates[template].filter->f0 / bank->templates[template].filter->deltaF, bank->templates[template].filter->data->length);
		if(XLALREAL8ReverseFFT(channel, fcorr, reverseplan)) {
			XLALDestroyCOMPLEX16Sequence(fcorr);
			XLALDestroyREAL8Sequence(channel);
			XLAL_ERROR_NULL(func, XLAL_EFUNC);
		}

	/* 2 * template bandwidth == channel samples (degrees of freedom)
	 * per second
	 *
	 * channel samples per second * (seconds per tseries sample) ==
	 * degrees of freedom per tseries sample
	 *
	 * 2 degrees of freedom / (degrees of freedom per tseries sample)
	 * == tseries samples for a 2 d-o-f tile
	 */
	for(length =  2 / (2 * bank->templates[template].bandwidth * plane->deltaT); length <= plane->tiles.max_length; length *= 2) {
		const double tile_dof = length * (2 * bank->templates[template].bandwidth * plane->deltaT);
		/* number of degrees of freedom in tile = number of
		 * "virtual pixels" in tile.  compute distance between
		 * tile's "virtual pixels" */
		const unsigned stride = length / tile_dof;

	for(end = (start = plane->tiles.tiling_start) + length; end <= plane->tiles.tiling_end; end = (start += length / plane->tiles.inv_fractional_stride) + length) {
		double sumsquares = 0;
		unsigned t;

		/* compute sum of squares */
		for(t = start + stride / 2; t < end; t += stride)
			sumsquares += channel->data[t] * channel->data[t];
		/* normalization to give each pixel a mean square of 1.0 */
		sumsquares /= sum_mean_square;

		/* compute statistical confidence, see if it's above
		 * threshold */
		confidence = -XLALlnOneMinusChisqCdf(sumsquares, tile_dof);
		if(XLALIsREAL8FailNaN(confidence)) {
			XLALDestroyCOMPLEX16Sequence(fcorr);
			XLALDestroyREAL8Sequence(channel);
			XLAL_ERROR_NULL(func, XLAL_EFUNC);
		}
		if(confidence >= confidence_threshold) {
			SnglBurstTable *oldhead;

			/* compute h_rss */
			h_rss = sqrt((sumsquares - tile_dof) * bank->templates[template].unwhitened_mean_square * stride * plane->deltaT);

			/* add new event to head of linked list */
			oldhead = head;
			head = XLALTFTileToBurstEvent(plane, start, length, bank->templates[template].f_centre, bank->templates[template].bandwidth, h_rss, sumsquares, tile_dof, confidence);
			if(!head) {
				XLALDestroyCOMPLEX16Sequence(fcorr);
				XLALDestroyREAL8Sequence(channel);
				XLAL_ERROR_NULL(func, XLAL_EFUNC);
			}
			head->next = oldhead;
		}
	}
	}
	}

	/* clean up */
	XLALDestroyCOMPLEX16Sequence(fcorr);
	XLALDestroyREAL8Sequence(channel);

	/* success */
	return(head);
}
