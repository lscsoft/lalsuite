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
	double excess_power,
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
	event->snr = excess_power;
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
	double excess_power;
	double h_rss;
	double confidence;

	for(channels = plane->tiles.min_channels; channels <= plane->tiles.max_channels; channels *= 2) {
	for(channel_end = (channel = 0) + channels; channel_end <= plane->channels; channel_end = (channel += channels / plane->tiles.inv_fractional_stride) + channels) {
		/* sum of inner products of all pairs of channels
		 * contributing to this tile's "virtual channel" */
		const double twice_channel_overlap = XLALREAL8SequenceSum(plane->twice_channel_overlap, channel, channels - 1);
		/* mean square for this tile's "virtual channel" */
		const double pixel_mean_square = XLALREAL8SequenceSumSquares(plane->channel_rms, channel, channels) / (channels + twice_channel_overlap);
	/* start with at least 2 degrees of freedom */
	for(length = 2 / (channels * plane->tiles.dof_per_pixel); length <= plane->tiles.max_length; length *= 2) {
		const double tile_dof = (length * channels) * plane->tiles.dof_per_pixel;
		/* number of degrees of freedom in tile = number of
		 * "virtual pixels" in tile.  compute distance between
		 * tile's "virtual pixels" */
		const unsigned tstep = length / tile_dof;
	for(end = (start = plane->tiles.tiling_start) + length; end <= plane->tiles.tiling_end; end = (start += length / plane->tiles.inv_fractional_stride) + length) {
		double sumsquares = 0.0;
		double hsumsquares = 0.0;
		unsigned c;
		unsigned t;

		/* compute sum of snr squares (for excess power), and sum
		 * of energies (for h_rss) */
		for(t = start + tstep / 2; t < end; t += tstep) {
			double sum = 0.0;
			double hsum = 0.0;

			/* compute snr and dewhitened amplitude for this
			 * "virtual pixel". */
			for(c = channel; c < channel_end; c++) {
				sum += plane->channel[c]->data[t];
				hsum += plane->channel[c]->data[t] * plane->channel_rms->data[c];
			}

			sumsquares += sum * sum;
			hsumsquares += hsum * hsum;
		}
		/* normalization to give each "virtual pixel" a mean square
		 * snr of 1.0 */
		sumsquares /= channels + twice_channel_overlap;
		hsumsquares /= channels + twice_channel_overlap;

		/* compute statistical confidence, see if it's above
		 * threshold */
		confidence = -XLALlnOneMinusChisqCdf(sumsquares, tile_dof);
		if(XLALIsREAL8FailNaN(confidence))
			XLAL_ERROR_NULL(func, XLAL_EFUNC);
		if(confidence >= confidence_threshold) {
			SnglBurstTable *oldhead;

			/* compute excess power and h_rss */
			excess_power = sumsquares - tile_dof;
			h_rss = sqrt(hsumsquares - tile_dof * pixel_mean_square);

			/* add new event to head of linked list */
			oldhead = head;
			head = XLALTFTileToBurstEvent(plane, start, length, channel, channels, h_rss, excess_power, confidence);
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
