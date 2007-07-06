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


#include <math.h>
#include <lal/LALRCSID.h>


NRCSID (COMPUTEEXCESSPOWERC, "$Id$");


#include <lal/Date.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/Sequence.h>
#include <lal/TFTransform.h>
#include <lal/Thresholds.h>
#include <lal/XLALError.h>


/*
 * Convert an array of tiles to a linked list of burst events.  The
 * threshold cut is applied here.
 */


static SnglBurstTable *XLALTFTileToBurstEvent(
	const REAL4TimeFrequencyPlane *plane,
	unsigned tile_start,
	unsigned tile_length,
	unsigned tile_channel,
	unsigned tile_channels,
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
	event->bandwidth = tile_channels * plane->deltaF;
	event->central_freq = plane->flow + tile_channel * plane->deltaF + (0.5 * event->bandwidth);
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


/******** <lalVerbatim file="ComputeExcessPowerCP"> ********/
SnglBurstTable *XLALComputeExcessPower(
	const REAL4TimeFrequencyPlane *plane,
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
		/* the mean square of the "virtual channel" */
		const double sum_mean_square = channels * plane->deltaF / plane->fseries_deltaF + XLALREAL8SequenceSum(plane->twice_channel_overlap, channel, channels - 1);
		/* the mean square of the "uwsum" quantity computed below,
		 * which is proportional to an approximation of the
		 * unwhitened time series. */
		double uwsum_mean_square;
		/* true unwhitened mean square for this channel.  the ratio
		 * of this to uwsum_mean_square is the correction factor to
		 * be applied to uwsum^{2} to convert it to an
		 * approximation of the square of the unwhitened channel */
		const double unwhitened_mean_square = XLALREAL8SequenceSumSquares(plane->unwhitened_rms, channel, channels) + XLALREAL8SequenceSum(plane->unwhitened_cross, channel, channels - 1);
		unsigned c;

		/* compute uwsum_mean_square */
		uwsum_mean_square = XLALREAL8SequenceSumSquares(plane->unwhitened_rms, channel, channels);
		for(c = channel; c < channel_end - 1; c++)
			uwsum_mean_square += plane->twice_channel_overlap->data[c] * plane->unwhitened_rms->data[c] * plane->unwhitened_rms->data[c + 1] * plane->fseries_deltaF / plane->deltaF;

	/* start with at least 2 degrees of freedom */
	for(length = 2 / (channels * plane->tiles.dof_per_pixel); length <= plane->tiles.max_length; length *= 2) {
		const double tile_dof = (length * channels) * plane->tiles.dof_per_pixel;
		/* number of degrees of freedom in tile = number of
		 * "virtual pixels" in tile.  compute distance between
		 * tile's "virtual pixels" */
		const unsigned tstep = length / tile_dof;
	for(end = (start = plane->tiles.tiling_start) + length; end <= plane->tiles.tiling_end; end = (start += length / plane->tiles.inv_fractional_stride) + length) {
		double sumsquares = 0;
		double uwsumsquares = 0;
		unsigned t;

		/* compute sum of squares, and unwhitened sum of squares */
		for(t = start + tstep / 2; t < end; t += tstep) {
			double sum = 0;
			double uwsum = 0;

			/* compute snr and dewhitened amplitude for this
			 * "virtual pixel". */
			for(c = channel; c < channel_end; c++) {
				sum += plane->channel[c]->data[t];
				uwsum += plane->channel[c]->data[t] * plane->unwhitened_rms->data[c] * sqrt(plane->fseries_deltaF / plane->deltaF);
			}

			sumsquares += sum * sum;
			uwsumsquares += uwsum * uwsum;
		}
		/* normalization to give each pixel a mean square of 1.0 */
		sumsquares /= sum_mean_square;

		/* compute statistical confidence, see if it's above
		 * threshold */
		/* FIXME:  the 0.65 is an empirically determined
		 * degree-of-freedom fudge factor.  figure out what its
		 * origin is, and account for it correctly.  it's most
		 * likely due to the time-frequency plane pixels not being
		 * independent of one another as a consequence of a
		 * non-zero inner product of the time-domain impulse
		 * response of the channel filter for adjacent pixels */
		confidence = -XLALlnOneMinusChisqCdf(sumsquares * .65, tile_dof * .65);
		if(XLALIsREAL8FailNaN(confidence))
			XLAL_ERROR_NULL(func, XLAL_EFUNC);
		if(confidence >= confidence_threshold) {
			SnglBurstTable *oldhead;

			/* normalization to give each unwhitened pixel a
			 * mean square of 1.0 */
			uwsumsquares /= uwsum_mean_square;

			/* compute h_rss */
			h_rss = sqrt((uwsumsquares - tile_dof) * unwhitened_mean_square * tstep * plane->deltaT);

			/* add new event to head of linked list */
			oldhead = head;
			head = XLALTFTileToBurstEvent(plane, start, length, channel, channels, h_rss, sumsquares, tile_dof, confidence);
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
