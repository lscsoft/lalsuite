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


#ifndef _TFTRANSFORM_H
#define _TFTRANSFORM_H


#include <lal/LALDatatypes.h>
#include <lal/Window.h>
#include <lal/RealFFT.h>
#include <lal/LALRCSID.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/Sequence.h>


#ifdef  __cplusplus		/* C++ protection. */
extern "C" {
#endif


NRCSID(TFTRANSFORMH, "$Id$");


/*
 * A time-frequency plane
 */


typedef struct tagREAL8TimeFrequencyPlane {
	/* name of data from which this was computed */
	CHAR name[LALNameLength];
	/* epoch of data from which this was computed */
	LIGOTimeGPS epoch;
	/* time resolution of the plane */
	REAL8 deltaT;
	/* input frequency series' resolution */
	REAL8 fseries_deltaF;
	/* Number of frequency channels in TF plane */
	UINT4 channels;
	/* TF plane's frequency resolution (channel spacing) */
	REAL8 deltaF;
	/* low frequency boundary of TF plane */
	REAL8 flow;
	/* channel filters */
	COMPLEX16FrequencySeries **filter;
	/* twice the inner product of filters for neighbouring channels;
	 * twice_channel_overlap[0] is twice the inner product of the
	 * filters for channels 0 and 1, and so on (for n channels, there
	 * are n - 1 channel_overlaps) */
	REAL8Sequence *twice_channel_overlap;
	/* root mean square for the unwhitened time series corresponding to
	 * each channel, and the mean square cross terms for wide channels
	 * (indices same as for twice_channel_overlap) */
	REAL8Sequence *unwhitened_rms;
	REAL8Sequence *unwhitened_cross;
	/* channel data.  channel[j]->data[i] corresponds to time
	 *
	 * epoch + i * deltaT
	 *
	 * and the frequency band
	 *
	 * [flow + j * deltaF, flow + (j + 1) * deltaF)
	 */
	REAL8Sequence **channel;
	/* time-frequency plane's tiling information */
	struct TFTiling {
		unsigned max_length;
		unsigned min_channels;
		unsigned max_channels;
		unsigned tiling_start;
		unsigned tiling_end;
		unsigned inv_fractional_stride;
		double dof_per_pixel;
	} tiles;
	/* window applied to input time series for tapering edges to 0 */
	REAL8Window *window;
	/* by how many samples a window's start should be shifted from the
	 * start of the window preceding it */
	INT4 window_shift;
	/* two-point spectral correlation of the whitened frequency series,
	 * computed from the time-domain window function */
	REAL8Sequence *two_point_spectral_correlation;
} REAL8TimeFrequencyPlane;


REAL8TimeFrequencyPlane *XLALCreateTFPlane(
	UINT4 tseries_length,
	REAL8 tseries_deltaT,
	REAL8 flow,
	REAL8 bandwidth,
	REAL8 tiling_fractional_stride,
	REAL8 tiling_max_bandwidth,
	REAL8 tiling_max_duration
);


void XLALDestroyTFPlane(
	REAL8TimeFrequencyPlane *plane
);


INT4 XLALTFPlaneMakeChannelFilters(
	REAL8TimeFrequencyPlane *plane,
	const REAL8FrequencySeries *psd
);


int XLALFreqSeriesToTFPlane(
	REAL8TimeFrequencyPlane *tfplane,
	const COMPLEX16FrequencySeries *fseries,
	const REAL8FFTPlan *reverseplan
);


SnglBurstTable *XLALComputeExcessPower(
	const REAL8TimeFrequencyPlane *plane,
	SnglBurstTable *head,
	double confidence_threshold
);


INT4 XLALOverlappedSegmentsCommensurate(
	INT4 target_length,
	INT4 segment_length,
	INT4 segment_shift
);


INT4 XLALEPGetTimingParameters(
	INT4 window_length,
	INT4 max_tile_length,
	REAL8 fractional_tile_stride,
	INT4 *psd_length,
	INT4 *psd_shift,
	INT4 *window_shift,
	INT4 *window_pad,
	INT4 *tiling_length
);


/*
 * An excess power "template bank"
 */


struct ExcessPowerTemplateBank {
	struct ExcessPowerTemplate {
		COMPLEX16FrequencySeries *filter;
		REAL8 f_centre;
		REAL8 bandwidth;
		REAL8 unwhitened_mean_square;
	} *templates;
	int n_templates;
};


struct ExcessPowerTemplateBank *XLALCreateExcessPowerTemplateBank(
	const COMPLEX16FrequencySeries *template,
	const REAL8TimeFrequencyPlane *plane,
	const REAL8FrequencySeries *psd
);


void XLALDestroyExcessPowerTemplateBank(
	struct ExcessPowerTemplateBank *bank
);


SnglBurstTable *XLALExcessPowerProject(
	const COMPLEX16FrequencySeries *fseries,
	REAL8TimeFrequencyPlane *plane,
	const struct ExcessPowerTemplateBank *bank,
	SnglBurstTable *head,
	double confidence_threshold,
	const REAL8FFTPlan *reverseplan
);


#ifdef  __cplusplus
}
#endif				/* C++ protection. */
#endif
