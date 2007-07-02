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
#include <lal/LALAtomicDatatypes.h>
#include <lal/BandPassTimeSeries.h>
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
 * Generate a linked list of burst events from a time series.
 */


SnglBurstTable *XLALEPSearch(
	const REAL4TimeSeries *tseries,
	EPSearchParams *params
)
{ 
	const char func[] = "XLALEPSearch";
	SnglBurstTable *head = NULL;
	int errorcode = 0;
	int start_sample;
	COMPLEX8FrequencySeries *fseries = NULL;
	RealFFTPlan *fplan;
	RealFFTPlan *rplan;
	REAL4FrequencySeries *psd;
	REAL4TimeSeries *cuttseries;
	REAL4TimeFrequencyPlane *plane;

	/*
	 * Construct forward and reverse FFT plans, storage for the PSD,
	 * the time-frequency plane, and a tiling.  Note that the flat part
	 * of the Tukey window needs to match the locations of the tiles as
	 * specified by the tiling_start parameter of XLALCreateTFPlane.
	 */

	fplan = XLALCreateForwardREAL4FFTPlan(params->window->data->length, 1);
	rplan = XLALCreateReverseREAL4FFTPlan(params->window->data->length, 1);
	psd = XLALCreateREAL4FrequencySeries("PSD", &tseries->epoch, 0, 0, &lalDimensionlessUnit, params->window->data->length / 2 + 1);
	plane = XLALCreateTFPlane(params->window->data->length, tseries->deltaT, params->flow, params->bandwidth, params->fractional_stride, params->maxTileBandwidth, params->maxTileDuration);
	if(!fplan || !rplan || !psd || !plane) {
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
	XLALDestroyREAL4Window(plane->window);
	plane->window = XLALCreateRectangularREAL4Window(i);
	}
#endif

	/*
	 * Compute the average spectrum.
	 *
	 * FIXME: is using windowShift here correct?  we have to, otherwise
	 * the time series' lengths are inconsistent
	 */

	switch(params->method) {
	case useMean:
		if(XLALREAL4AverageSpectrumWelch(psd, tseries, plane->window->data->length, plane->window_shift, plane->window, fplan) < 0) {
			errorcode = XLAL_EFUNC;
			goto error;
		}
		break;

	case useMedian:
		if(XLALREAL4AverageSpectrumMedian(psd, tseries, plane->window->data->length, plane->window_shift, plane->window, fplan) < 0) {
			errorcode = XLAL_EFUNC;
			goto error;
		}
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

	for(start_sample = 0; start_sample + plane->window->data->length <= tseries->data->length; start_sample += plane->window_shift) {
		/*
		 * Extract a window-length of data from the time series,
		 * compute its DFT, then free it.
		 */

		cuttseries = XLALCutREAL4TimeSeries(tseries, start_sample, plane->window->data->length);
		if(!cuttseries) {
			errorcode = XLAL_EFUNC;
			goto error;
		}
		XLALPrintInfo("XLALEPSearch(): analyzing %u samples (%.9lf s) at offset %u (%.9lf s) from epoch %d.%09u s\n", cuttseries->data->length, cuttseries->data->length * cuttseries->deltaT, start_sample, start_sample * cuttseries->deltaT, tseries->epoch.gpsSeconds, tseries->epoch.gpsNanoSeconds);
		if(params->diagnostics)
			params->diagnostics->XLALWriteLIGOLwXMLArrayREAL4TimeSeries(params->diagnostics->LIGOLwXMLStream, NULL, cuttseries);

		XLALPrintInfo("XLALEPSearch(): computing the Fourier transform\n");
		fseries = XLALWindowedREAL4ForwardFFT(cuttseries, plane->window, fplan);
		XLALDestroyREAL4TimeSeries(cuttseries);
		if(!fseries) {
			errorcode = XLAL_EFUNC;
			goto error;
		}

		/*
		 * Normalize the frequency series to the average PSD.
		 */

#if 1
		XLALPrintInfo("XLALEPSearch(): normalizing to the average spectrum\n");
		if(!XLALWhitenCOMPLEX8FrequencySeries(fseries, psd, plane->flow, plane->flow + plane->channels * plane->deltaF)) {
			errorcode = XLAL_EFUNC;
			goto error;
		}
		if(params->diagnostics)
			params->diagnostics->XLALWriteLIGOLwXMLArrayCOMPLEX8FrequencySeries(params->diagnostics->LIGOLwXMLStream, "whitened", fseries);
#endif

		/*
		 * Compute the time-frequency plane from the frequency
		 * series.
		 */

		if(XLALFreqSeriesToTFPlane(plane, fseries, psd, rplan, params->useOverWhitening)) {
			errorcode = XLAL_EFUNC;
			goto error;
		}
		XLALDestroyCOMPLEX8FrequencySeries(fseries);
		fseries = NULL;

		/*
		 * Compute the excess power for each time-frequency tile
		 * using the data in the time-frequency plane, and add
		 * those whose confidence is above threshold to the trigger
		 * list.  Note that because it is possible for there to be
		 * 0 triggers found, we can't check for errors by testing
		 * for head == NULL.
		 */

		XLALPrintInfo("XLALEPSearch(): computing the excess power for each tile\n");
		XLALClearErrno();
		head = XLALComputeExcessPower(plane, head, params->confidence_threshold);
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
	XLALDestroyTFPlane(plane);
	if(errorcode) {
		XLALDestroySnglBurstTable(head);
		XLAL_ERROR_NULL(func, errorcode);
	}
	return(head);
}


/*
 * Condition the time series prior to analysis by the power code
 */


int XLALEPConditionData(
	REAL4TimeSeries *series,
	REAL8 flow,
	REAL8 resampledeltaT,
	INT4 corruption
)
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
