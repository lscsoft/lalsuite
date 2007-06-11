/*
 * $Id$
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


#include <math.h>
#include <stdio.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/Date.h>
#include <lal/EPSearch.h>
#include <lal/FrequencySeries.h>
#include <lal/LALDatatypes.h>
#include <lal/LALRCSID.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALStdlib.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/RealFFT.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/TFTransform.h>
#include <lal/TimeFreqFFT.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/Window.h>
#include <lal/XLALError.h>


NRCSID(EPSEARCHC, "$Id$");


/*
 * Delete a SnglBurstTable linked list
 */


static void XLALDestroySnglBurstTable(SnglBurstTable *head)
{
	SnglBurstTable *event;

	while(head) {
		event = head;
		head = head->next;
		XLALFree(event);
	}
}


/*
 * Convert an array of tiles to a linked list of burst events.  The
 * threshold cut is applied here.
 */


static SnglBurstTable *XLALTFTileToBurstEvent(
	const TFTile *tile,
	const char *channelName,
	const LIGOTimeGPS *epoch
)
{
	const char func[] = "XLALTFTileToBurstEvent";
	SnglBurstTable *event = XLALCalloc(1, sizeof(*event));
	if(!event)
		XLAL_ERROR_NULL(func, XLAL_ENOMEM);

	event->next = NULL;
	strncpy(event->ifo, channelName, 2);
	event->ifo[2] = '\0';
	strncpy(event->search, "excesspower", LIGOMETA_SEARCH_MAX);
	event->search[LIGOMETA_SEARCH_MAX] = '\0';
	strncpy(event->channel, channelName, LIGOMETA_CHANNEL_MAX);
	event->channel[LIGOMETA_CHANNEL_MAX] = '\0';

	event->start_time = *epoch; 
 
	XLALGPSAdd(&event->start_time, tile->tstart * tile->deltaT);
	event->duration = (tile->tend - tile->tstart) * tile->deltaT;
	event->peak_time = event->start_time;
	XLALGPSAdd(&event->peak_time, 0.5 * event->duration);
	event->bandwidth = tile->channels * tile->deltaF;
	event->central_freq = tile->flow + tile->channel0 * tile->deltaF + (0.5 * event->bandwidth);
	/* FIXME: put hrss into the "hrss" column */
	event->amplitude = tile->hrss;
	event->snr = tile->excessPower;
	event->confidence = tile->confidence;
	event->string_cluster_t = XLAL_REAL4_FAIL_NAN;
	event->event_id = 0;

	return(event);
}


static SnglBurstTable *XLALTFTilesToSnglBurstTable(SnglBurstTable *head, const REAL4TimeFrequencyPlane *plane, const TFTiling *tiling, REAL8 confidence_threshold)
{
	const char func[] = "XLALTFTilesToSnglBurstTable";
	SnglBurstTable *oldhead;
	TFTile *tile;
	size_t i;

	for(i = 0, tile = tiling->tile; i < tiling->numtiles; i++, tile++) {
		if(tile->confidence >= confidence_threshold) {
			oldhead = head;
			head = XLALTFTileToBurstEvent(tile, plane->name, &plane->epoch); 
			if(!head) {
				XLALDestroySnglBurstTable(oldhead);
				XLAL_ERROR_NULL(func, XLAL_EFUNC);
			}
			head->next = oldhead;
		}
	}

	return(head);
}


/*
 * Generate a linked list of burst events from a time series.
 */


/******** <lalVerbatim file="EPSearchCP"> ********/
SnglBurstTable *XLALEPSearch(
	const REAL4TimeSeries *tseries,
	EPSearchParams *params
)
/******** </lalVerbatim> ********/
{ 
	const char func[] = "XLALEPSearch";
	SnglBurstTable *head = NULL;
	int errorcode = 0;
	int start_sample;
	COMPLEX8FrequencySeries *fseries = NULL;
	REAL4Window *tukey;
	RealFFTPlan *fplan;
	RealFFTPlan *rplan;
	REAL4FrequencySeries *psd;
	REAL4TimeSeries *cuttseries;
	REAL4TimeFrequencyPlane *tfplane;

	/*
	 * FreqSeriesToTFPlane() is passed a frequency series, and it needs
	 * to know how many samples the time series from which it was
	 * generated contained.  This can only be done if you assume the
	 * number of samples was even (or odd, but you have to pick one).
	 * So we need to make sure that's true.
	 */

	if(params->window->data->length & 1) {
		/* window length is odd */
		XLAL_ERROR_NULL(func, XLAL_EINVAL);
	}

	/*
	 * Construct forward and reverse FFT plans, storage for the PSD,
	 * the time-frequency plane, and a tiling.  Note that the flat part
	 * of the Tukey window needs to match the locations of the tiles as
	 * specified by the tiling_start parameter of XLALCreateTFPlane.
	 */

	fplan = XLALCreateForwardREAL4FFTPlan(params->window->data->length, 1);
	rplan = XLALCreateReverseREAL4FFTPlan(params->window->data->length, 1);
	psd = XLALCreateREAL4FrequencySeries("PSD", &tseries->epoch, 0, 0, &lalDimensionlessUnit, params->window->data->length / 2 + 1);
	tfplane = XLALCreateTFPlane(params->window->data->length, tseries->deltaT, params->tf_freqBins, params->tf_deltaF, params->tf_flow, params->window->data->length / 4, params->inv_fractional_stride, params->maxTileBandwidth, params->maxTileDuration);
	tukey = XLALCreateTukeyREAL4Window(params->window->data->length, 0.5);
	if(!fplan || !rplan || !psd || !tfplane || !tukey) {
		errorcode = XLAL_EFUNC;
		goto error;
	}

	/*
	 * Adjust the Tukey window's normalization so that it is
	 * RMS-preserving in the flat portion.  This is done by setting its
	 * normalization to be equal to that of a rectangular window
	 * (pretend the tapers aren't present).
	 */

	tukey->sumofsquares = tukey->data->length;

	/*
	 * Compute the average spectrum.
	 */

	switch(params->method) {
	case useMean:
		XLALREAL4AverageSpectrumWelch(psd, tseries, params->window->data->length, params->windowShift, params->window, fplan);
		break;

	case useMedian:
		XLALREAL4AverageSpectrumMedian(psd, tseries, params->window->data->length, params->windowShift, params->window, fplan);
		break;

	default:
		errorcode = XLAL_EINVAL;
		goto error;
	}

	if(params->diagnostics)
		params->diagnostics->XLALWriteLIGOLwXMLArrayREAL4FrequencySeries(params->diagnostics->LIGOLwXMLStream, NULL, psd);

	/*
	 * Loop over data applying excess power method.
	 */

	for(start_sample = 0; start_sample + params->window->data->length <= tseries->data->length; start_sample += params->windowShift) {
		/*
		 * Extract a window-length of data from the time series,
		 * compute its DFT, then free it.
		 */

		cuttseries = XLALCutREAL4TimeSeries(tseries, start_sample, params->window->data->length);
		if(!cuttseries) {
			errorcode = XLAL_EFUNC;
			goto error;
		}
		XLALPrintInfo("XLALEPSearch(): analyzing %u samples (%.9lf s) at offset %u (%.9lf s) from epoch %d.%09u s\n", cuttseries->data->length, cuttseries->data->length * cuttseries->deltaT, start_sample, start_sample * cuttseries->deltaT, tseries->epoch.gpsSeconds, tseries->epoch.gpsNanoSeconds);
		if(params->diagnostics)
			params->diagnostics->XLALWriteLIGOLwXMLArrayREAL4TimeSeries(params->diagnostics->LIGOLwXMLStream, NULL, cuttseries);

		XLALPrintInfo("XLALEPSearch(): computing the Fourier transform\n");
		fseries = XLALWindowedREAL4ForwardFFT(cuttseries, tukey, fplan);
		XLALDestroyREAL4TimeSeries(cuttseries);
		if(!fseries) {
			errorcode = XLAL_EFUNC;
			goto error;
		}

		/*
		 * Normalize the frequency series to the average PSD.
		 * FIXME: the psd has been computed from data that was
		 * windowed in the time domain with a different window than
		 * has been used in optaining fseries.  This means that
		 * whatever spectral leakage is induced by windowing has
		 * not occured identically in the two.  The windows are
		 * pretty close to each other (Hann vs. 50% Tukey) so it's
		 * almost certainly not an issue, but it might be worth one
		 * day reassuring ourselves that the PSD is a good PSD for
		 * the fseries.
		 */

		XLALPrintInfo("XLALEPSearch(): normalizing to the average spectrum\n");
		if(!XLALWhitenCOMPLEX8FrequencySeries(fseries, psd, tfplane->flow, fseries->f0 + fseries->data->length * fseries->deltaF)) {
			errorcode = XLAL_EFUNC;
			goto error;
		}
		if(params->diagnostics)
			params->diagnostics->XLALWriteLIGOLwXMLArrayCOMPLEX8FrequencySeries(params->diagnostics->LIGOLwXMLStream, "whitened", fseries);

		/*
		 * Compute the time-frequency plane from the frequency
		 * series.
		 */

		XLALPrintInfo("XLALEPSearch(): computing the time-frequency decomposition\n");
		if(XLALFreqSeriesToTFPlane(tfplane, fseries, psd, rplan, params->useOverWhitening)) {
			errorcode = XLAL_EFUNC;
			goto error;
		}
		XLALDestroyCOMPLEX8FrequencySeries(fseries);
		fseries = NULL;

		/*
		 * Compute the excess power for each time-frequency tile
		 * using the data in the time-frequency plane.
		 */

		XLALPrintInfo("XLALEPSearch(): computing the excess power for each tile\n");
		if(XLALComputeExcessPower(tfplane)) {
			errorcode = XLAL_EFUNC;
			goto error;
		}

		/*
		 * Convert the TFTiles into sngl_burst events for output.
		 * Note that because it is possible for there to be 0
		 * triggers found, we can't check for errors by testing for
		 * head == NULL
		 */

		XLALPrintInfo("XLALEPSearch(): converting tiles to trigger list\n");
		XLALClearErrno();
		head = XLALTFTilesToSnglBurstTable(head, tfplane, tfplane->tiling, params->confidence_threshold);
		if(xlalErrno) {
			errorcode = XLAL_EFUNC;
			goto error;
		}
	}

	/*
	 * Memory clean-up.
	 */

	XLALPrintInfo("XLALEPSearch(): done\n");

	error:
	XLALDestroyREAL4FFTPlan(fplan);
	XLALDestroyREAL4FFTPlan(rplan);
	XLALDestroyREAL4FrequencySeries(psd);
	XLALDestroyCOMPLEX8FrequencySeries(fseries);
	XLALDestroyTFPlane(tfplane);
	XLALDestroyREAL4Window(tukey);
	if(errorcode) {
		XLALDestroySnglBurstTable(head);
		XLAL_ERROR_NULL(func, errorcode);
	}
	return(head);
}


/*
 * Condition the time series prior to analysis by the power code
 */


/* <lalVerbatim file="EPConditionDataCP"> */
int XLALEPConditionData(
	REAL4TimeSeries  *series,
	REAL8             flow,
	REAL8             resampledeltaT,
	INT4              corruption
)
/* </lalVerbatim> */
{
	const char func[] = "XLALEPConditionData";
	const REAL8         epsilon = 1.0e-3;
	PassBandParamStruc  highpassParam;
	size_t              newlength;

	/*
	 * Resample the time series if necessary
	 */

	if(fabs(resampledeltaT - series->deltaT) / series->deltaT >= epsilon)
		if(XLALResampleREAL4TimeSeries(series, resampledeltaT))
			XLAL_ERROR(func, XLAL_EFUNC);

	/*
	 * High-pass filter the time series.
	 */

	highpassParam.nMax = 8;
	highpassParam.f2 = flow;
	highpassParam.f1 = -1.0;
	highpassParam.a2 = 0.9;
	highpassParam.a1 = -1.0;
	if(XLALButterworthREAL4TimeSeries(series, &highpassParam))
		XLAL_ERROR(func, XLAL_EFUNC);

	/*
	 * The filter corrupts the ends of the time series.  Chop them off.
	 */

	newlength = series->data->length - 2 * corruption;
	if(XLALShrinkREAL4TimeSeries(series, corruption, newlength) != newlength)
		XLAL_ERROR(func, XLAL_EFUNC);

	return(0);
}
