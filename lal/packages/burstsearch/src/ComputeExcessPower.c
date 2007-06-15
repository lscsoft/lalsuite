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


#include <lal/Sequence.h>
#include <lal/TFTransform.h>
#include <lal/Thresholds.h>
#include <lal/XLALError.h>


/******** <lalVerbatim file="ComputeExcessPowerCP"> ********/
int XLALComputeExcessPower(
	const REAL4TimeFrequencyPlane *plane
)
/******** </lalVerbatim> ********/
{
	const char func[] = "XLALComputeExcessPower";
	size_t i;

	for(i = 0; i < plane->tiling->numtiles; i++) {
		TFTile *tile = &plane->tiling->tiles[i];
		const double channel_overlap = XLALREAL8SequenceSum(plane->channel_overlap, tile->channel0, tile->channels - 1);
		const double pixel_mean_square = XLALREAL8SequenceSumSquares(plane->channel_rms, tile->channel0, tile->channels) / (tile->channels + channel_overlap);
		const unsigned tstep = (tile->tend - tile->tstart) / tile->dof;
		double sumsquares = 0.0;
		double hsumsquares = 0.0;
		unsigned t;

		for(t = tile->tstart + tstep / 2; t < tile->tend; t += tstep) {
			unsigned channel;
			double sum = 0.0;
			double hsum = 0.0;

			for(channel = tile->channel0; channel < tile->channel0 + tile->channels; channel++) {
				sum += plane->channel[channel]->data[t];
				hsum += plane->channel[channel]->data[t] * plane->channel_rms->data[channel];
			}

			sumsquares += sum * sum / (tile->channels + channel_overlap);
			hsumsquares += hsum * hsum / (tile->channels + channel_overlap);
		}

		tile->excess_power = sumsquares - tile->dof;
		tile->h_rss = sqrt(hsumsquares - tile->dof * pixel_mean_square);
		tile->confidence = -XLALlnOneMinusChisqCdf(sumsquares, tile->dof);
		if(XLALIsREAL8FailNaN(tile->confidence))
			XLAL_ERROR(func, XLAL_EFUNC);
	}

	/* success */
	return 0;
}
