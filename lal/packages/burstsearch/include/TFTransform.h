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


#include <gsl/gsl_matrix.h>


#include <lal/LALDatatypes.h>
#include <lal/Window.h>
#include <lal/RealFFT.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/Sequence.h>


#ifdef  __cplusplus		/* C++ protection. */
extern "C" {
#endif


#include <lal/LALRCSID.h>
NRCSID(TFTRANSFORMH, "$Id$");


/*
 * An excess power filter bank
 */


typedef struct tagLALExcessPowerFilterBank {
	int n_filters;
	struct ExcessPowerFilter {
		COMPLEX16FrequencySeries *fseries;
		/* root mean square of the unwhitened time series
		 * corresponding to this filter */
		REAL8 unwhitened_rms;
	} *basis_filters;
	/* twice the inner product of filters for neighbouring channels;
	 * twice_channel_overlap[0] is twice the inner product of the
	 * filters for channels 0 and 1, and so on (for n channels, there
	 * are n - 1 channel_overlaps) */
	REAL8Sequence *twice_channel_overlap;
	/* the mean square cross terms for wide channels (indices same as
	 * for twice_channel_overlap) */
	REAL8Sequence *unwhitened_cross;
} LALExcessPowerFilterBank;


LALExcessPowerFilterBank *XLALCreateExcessPowerFilterBank(
	double filter_deltaF,
	double flow,
	double channel_bandwidth,
	int n_channels,
	const REAL8FrequencySeries *psd,
	const REAL8Sequence *two_point_spectral_correlation
);


void XLALDestroyExcessPowerFilterBank(
	LALExcessPowerFilterBank *bank
);


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
	/* TF plane's frequency resolution (channel spacing) */
	REAL8 deltaF;
	/* low frequency boundary of TF plane */
	REAL8 flow;
	/* channel data.  each channel is placed into its own column.
	 * channel_data[i * channels + j] corresponds to time
	 *
	 * epoch + i * deltaT
	 *
	 * and the frequency band
	 *
	 * [flow + j * deltaF, flow + (j + 1) * deltaF)
	 */
	gsl_matrix *channel_data;
	/* re-usable holding area for the data for a single channel */
	REAL8Sequence *channel_buffer;
	REAL8Sequence *unwhitened_channel_buffer;
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
	/* time-domain window applied to input time series for tapering
	 * edges to 0 */
	REAL8Window *window;
	/* by how many samples a window's start should be shifted from the
	 * start of the window preceding it */
	INT4 window_shift;
	/* two-point spectral correlation of the whitened frequency series,
	 * computed from the time-domain window function */
	REAL8Sequence *two_point_spectral_correlation;
} REAL8TimeFrequencyPlane;



REAL8Sequence *XLALREAL8WindowTwoPointSpectralCorrelation(
	const REAL8Window *window,
	const REAL8FFTPlan *plan
);


REAL8TimeFrequencyPlane *XLALCreateTFPlane(
	UINT4 tseries_length,
	REAL8 tseries_deltaT,
	REAL8 flow,
	REAL8 bandwidth,
	REAL8 tiling_fractional_stride,
	REAL8 tiling_max_bandwidth,
	REAL8 tiling_max_duration,
	const REAL8FFTPlan *plan
);


void XLALDestroyTFPlane(
	REAL8TimeFrequencyPlane *plane
);


int XLALFreqSeriesToTFPlane(
	REAL8TimeFrequencyPlane *tfplane,
	const LALExcessPowerFilterBank *filter_bank,
	const COMPLEX16FrequencySeries *fseries,
	const REAL8FFTPlan *reverseplan
);


SnglBurst *XLALComputeExcessPower(
	const REAL8TimeFrequencyPlane *plane,
	const LALExcessPowerFilterBank *filter_bank,
	SnglBurst *head,
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


#ifdef  __cplusplus
}
#endif				/* C++ protection. */
#endif
