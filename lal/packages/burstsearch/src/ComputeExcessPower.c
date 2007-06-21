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
	unsigned t_start,
	unsigned t_length,
	unsigned t_channel_start,
	unsigned t_channels,
	double h_rss,
	double excess_power,
	double confidence
)
{
	const char func[] = "XLALTFTileToBurstEvent";
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
 
	XLALGPSAdd(&event->start_time, t_start * plane->deltaT);
	event->duration = t_length * plane->deltaT;
	event->peak_time = event->start_time;
	XLALGPSAdd(&event->peak_time, 0.5 * event->duration);
	event->bandwidth = t_channels * plane->deltaF;
	event->central_freq = plane->flow + t_channel_start * plane->deltaF + (0.5 * event->bandwidth);
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
	const char func[] = "XLALComputeExcessPower";
	const TFTiling *tiling = plane->tiling;
	unsigned t_channel_start;
	unsigned t_channels;
	unsigned t_start;
	unsigned t_length;
	double excess_power;
	double h_rss;
	double confidence;

	for(t_length = tiling->min_length; t_length <= tiling->max_length; t_length *= 2)
	for(t_channels = tiling->min_channels; t_channels <= tiling->max_channels; t_channels *= 2)
	for(t_start = tiling->tiling_t_start; t_start + t_length <= tiling->tiling_t_end; t_start += t_length / tiling->inv_fractional_stride)
	for(t_channel_start = 0; t_channel_start + t_channels <= tiling->tiling_n_channels; t_channel_start += t_channels / tiling->inv_fractional_stride) {
		/* number of degrees of freedom in tile */
		double t_dof = (t_length * t_channels) * tiling->dof_per_pixel;
		/* sum of inner products of all pairs of channels
		 * contributing to this tile's "virtual channel" */
		const double channel_overlap = XLALREAL8SequenceSum(plane->channel_overlap, t_channel_start, t_channels - 1);
		/* mean square for this tile's "virtual channel" */
		const double pixel_mean_square = XLALREAL8SequenceSumSquares(plane->channel_rms, t_channel_start, t_channels) / (t_channels + channel_overlap);
		/* distance between tile's time bins */
		const unsigned tstep = t_length / t_dof;
		double sumsquares = 0.0;
		double hsumsquares = 0.0;
		unsigned t;

		/* compute sum of squares, and de-whitened sum of squares
		 * */
		for(t = t_start + tstep / 2; t < t_start + t_length; t += tstep) {
			unsigned channel;
			double sum = 0.0;
			double hsum = 0.0;

			for(channel = t_channel_start; channel < t_channel_start + t_channels; channel++) {
				sum += plane->channel[channel]->data[t];
				hsum += plane->channel[channel]->data[t] * plane->channel_rms->data[channel];
			}

			sumsquares += sum * sum / (t_channels + channel_overlap);
			hsumsquares += hsum * hsum / (t_channels + channel_overlap);
		}

		/* compute excess power, de-whitened root-sum-squares, and
		 * statistical confidence */
		excess_power = sumsquares - t_dof;
		h_rss = sqrt(hsumsquares - t_dof * pixel_mean_square);
		confidence = -XLALlnOneMinusChisqCdf(sumsquares, t_dof);
		if(XLALIsREAL8FailNaN(confidence))
			XLAL_ERROR_NULL(func, XLAL_EFUNC);

		/* test confidence */
		if(confidence >= confidence_threshold) {
			SnglBurstTable *oldhead = head;
			head = XLALTFTileToBurstEvent(plane, t_start, t_length, t_channel_start, t_channels, h_rss, excess_power, confidence);
			if(!head)
				XLAL_ERROR_NULL(func, XLAL_EFUNC);
			head->next = oldhead;
		}
	}

	/* success */
	return(head);
}
