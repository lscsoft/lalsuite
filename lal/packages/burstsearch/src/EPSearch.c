/******** <lalVerbatim file="EPSearchCV"> ********
Author: Brady, P. and Cannon, K.
Revision: $Id$
********* </lalVerbatim> ********/

#include <math.h>
#include <stdio.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/Date.h>
#include <lal/EPSearch.h>
#include <lal/ExcessPower.h>
#include <lal/FrequencySeries.h>
#include <lal/LALDatatypes.h>
#include <lal/LALErrno.h>
#include <lal/LALRCSID.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALStdlib.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/RealFFT.h>
#include <lal/ResampleTimeSeries.h>
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
		LALFree(event);
	}
}


/*
 * Convert an array of tiles to a linked list of burst events.  The
 * threshold cut is applied here.
 */

static SnglBurstTable *XLALTFTileToBurstEvent(
	const TFTile *tile,
	const char *channelName,
	const LIGOTimeGPS *epoch,
	const EPSearchParams *params  
)
{
	const char *func = "XLALTFTileToBurstEvent";
	SnglBurstTable *event = LALCalloc(1, sizeof(*event));
	if(!event)
		XLAL_ERROR_NULL(func, XLAL_ENOMEM);

	event->next = NULL;
	strncpy(event->ifo, channelName, 2);
	event->ifo[2] = '\0';
	strncpy(event->search, "power", LIGOMETA_SEARCH_MAX);
	event->search[LIGOMETA_SEARCH_MAX] = '\0';
	strncpy(event->channel, channelName, LIGOMETA_CHANNEL_MAX);
	event->channel[LIGOMETA_CHANNEL_MAX] = '\0';

	event->start_time = *epoch; 
 
	XLALGPSAdd(&event->start_time, tile->tstart * tile->deltaT);
	event->duration = tile->tbins * tile->deltaT;
	event->peak_time = event->start_time;
	XLALGPSAdd(&event->peak_time, 0.5 * event->duration);
	event->bandwidth = tile->fbins * tile->deltaF;
	event->central_freq = params->tfPlaneParams.flow + tile->fstart*tile->deltaF + (0.5 * event->bandwidth);
	event->amplitude = tile->hrss;
	event->snr = tile->excessPower;
	event->confidence =  tile->lnalpha;
	event->tfvolume = XLALTFTileDegreesOfFreedom(tile) / 2.0;
	event->string_cluster_t = XLAL_REAL4_FAIL_NAN;
	event->event_id = 0;

	return(event);
}


static SnglBurstTable *XLALTFTilesToSnglBurstTable(SnglBurstTable *head, const TFTiling *tiling, const char *channelName, const LIGOTimeGPS *epoch, const EPSearchParams *params)
{
	const char *func = "XLALTFTilesToSnglBurstTable";
	SnglBurstTable *oldhead;
	TFTile *tile;
	size_t i;

	for(i = 0, tile = tiling->tile; i < tiling->numtiles; i++, tile++) {
		if(tile->lnalpha <= params->lnalphaThreshold) {
			oldhead = head;
			head = XLALTFTileToBurstEvent(tile, channelName, epoch, params); 
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
 * Print a frequency series.
 */

static void print_real4fseries(const REAL4FrequencySeries *fseries, const char *file)
{
#if 0
	/* FIXME: why can't the linker find this function? */
	LALSPrintFrequencySeries(fseries, file);
#else
	FILE *fp = fopen(file, "w");
	size_t i;

	if(fp) {
		for(i = 0; i < fseries->data->length; i++)
			fprintf(fp, "%f\t%g\n", i * fseries->deltaF, fseries->data->data[i]);
		fclose(fp);
	}
#endif
}

static void print_complex8fseries(const COMPLEX8FrequencySeries *fseries, const char *file)
{
#if 0
	/* FIXME: why can't the linker find this function? */
	LALCPrintFrequencySeries(fseries, file);
#else
	FILE *fp = fopen(file, "w");
	size_t i;

	if(fp) {
		for(i = 0; i < fseries->data->length; i++)
			fprintf(fp, "%f\t%g\n", i * fseries->deltaF, sqrt(fseries->data->data[i].re * fseries->data->data[i].re + fseries->data->data[i].im * fseries->data->data[i].im));
		fclose(fp);
	}
#endif
}


/*
 * Normalize a complex8 fseries to a real4 average psd so that the rms of Re or
 * Im is 1.  (i.e. whiten the data).
 */

static void whiten(COMPLEX8FrequencySeries *fseries, const REAL4FrequencySeries *psd)
{
	REAL4 factor;
	size_t i;

	for(i = 0; i < fseries->data->length; i++) {
		/* FIXME: the computation of the average PSD sometimes
		 * underflows at low frequencies due to the strength of the
		 * high-pass filter(s) used in the data conditioning phase.
		 * This is OK when it happens outside the band of interest,
		 * but it is an error for this to occur in the band of
		 * interest.  We have *never* seen this happen, but it
		 * still might be worth adding some sort of check. */
		if(psd->data->data[i] == 0.0)
			factor = 0.0;
		else
			factor = 2.0 * sqrt(fseries->deltaF / psd->data->data[i]);
		fseries->data->data[i].re *= factor;
		fseries->data->data[i].im *= factor;
	}
}


/*
 * Generate a linked list of burst events from a time series.
 */

/******** <lalVerbatim file="EPSearchCP"> ********/
SnglBurstTable *
XLALEPSearch(
	const COMPLEX8FrequencySeries  *hrssresponse,
	const REAL4TimeSeries  *tseries,
	EPSearchParams   *params
)
/******** </lalVerbatim> ********/
{ 
	static const char *func = "EPSearch";
	SnglBurstTable *head = NULL;
	int errorcode;
	int                      start_sample;
	COMPLEX8FrequencySeries *fseries = NULL;
	REAL4Window             *window = params->window;
	RealFFTPlan             *fplan;
	RealFFTPlan             *rplan;
	REAL4FrequencySeries    *psd;
	REAL4TimeSeries         *cuttseries;
	TFTiling                *Tiling;
	REAL4TimeFrequencyPlane *tfplane;
	REAL4                   *normalisation;
	REAL8                   *hrssfactor;
	const LIGOTimeGPS        gps_zero = LIGOTIMEGPSZERO;

	/*
	 * FreqSeriesToTFPlane() is passed a frequency series, and it needs
	 * to know how many samples the time series from which it was
	 * generated contained.  This can only be done if you assume the
	 * number of samples was even (or odd, but you have to pick one).
	 * So we need to make sure that's true.
	 */

	if(window->data->length & 1) {
		/* window length is odd */
		XLAL_ERROR_NULL(func, XLAL_EINVAL);
	}

	/*
	 * Create an FFT plan, allocate space for the average spectrum,
	 * allocate temporary storage for frequency series data, allocate
	 * storage for the normalisation data, allocate and initialize the
	 * time-frequency plane storage, and construct a time-frequency
	 * tiling of the plane.
	 */

	fplan = XLALCreateForwardREAL4FFTPlan(window->data->length, 1);
	rplan = XLALCreateReverseREAL4FFTPlan(window->data->length, 1);
	psd = XLALCreateREAL4FrequencySeries("anonymous", &gps_zero, 0, 0, &lalDimensionlessUnit, window->data->length / 2 + 1);
	tfplane = XLALCreateTFPlane(&params->tfPlaneParams);
	normalisation = LALMalloc(tfplane->params.freqBins * sizeof(*normalisation));
	hrssfactor = LALMalloc(tfplane->params.freqBins * sizeof(*hrssfactor));
	Tiling = XLALCreateTFTiling(&params->tfTilingInput, &tfplane->params);

	if(!normalisation || !hrssfactor) {
		errorcode = XLAL_ENOMEM;
		goto error;
	}
	if(!fplan || !rplan || !psd || !tfplane || !Tiling) {
		errorcode = XLAL_EFUNC;
		goto error;
	}

	/*
	 * Compute the average spectrum.
	 */

	switch(params->method) {
	case useMean:
		XLALREAL4AverageSpectrumWelch(psd, tseries, window->data->length, params->windowShift, window, fplan);
		break;

	case useMedian:
		XLALREAL4AverageSpectrumMedian(psd, tseries, window->data->length, params->windowShift, window, fplan);
		break;

	default:
		errorcode = XLAL_EINVAL;
		goto error;
	}

	if(params->printSpectrum)
		print_real4fseries(psd, params->printSpectrum);

	/*
	 * If doing hrss estimation, check the response function's
	 * parameters
	 */

	if(hrssresponse) {
		if((hrssresponse->f0 != psd->f0) ||
		   (hrssresponse->deltaF != psd->deltaF) ||
		   (hrssresponse->data->length != psd->data->length)) {
		errorcode = XLAL_EINVAL;
		goto error;
		}
	}

	/*
	 * Loop over data applying excess power method.
	 */

	for(start_sample = 0; start_sample + window->data->length <= tseries->data->length; start_sample += params->windowShift) {
		/*
		 * Extract a window-length of data from the time series,
		 * compute its DFT, then free it.
		 */

		cuttseries = XLALCutREAL4TimeSeries(tseries, start_sample, window->data->length);
		if(!cuttseries) {
			errorcode = XLAL_EFUNC;
			goto error;
		}
		XLALPrintInfo("XLALEPSearch(): analyzing samples %zu -- %zu (%.9lf s -- %.9lf s)\n", start_sample, start_sample + cuttseries->data->length, start_sample * cuttseries->deltaT, (start_sample + cuttseries->data->length) * cuttseries->deltaT);

		XLALPrintInfo("XLALEPSearch(): computing the Fourier transform\n");
		fseries = XLALComputeFrequencySeries(cuttseries, window, fplan);
		XLALDestroyREAL4TimeSeries(cuttseries);
		if(!fseries) {
			errorcode = XLAL_EFUNC;
			goto error;
		}

		/*
		 * Normalize the frequency series to the average PSD.
		 */

		XLALPrintInfo("XLALEPSearch(): normalizing to the average spectrum\n");
		whiten(fseries, psd);
		if(params->printSpectrum)
			print_complex8fseries(fseries, "frequency_series.dat");

		/*
		 * Compute the time-frequency plane from the frequency
		 * series.
		 */

		XLALPrintInfo("XLALEPSearch(): computing the time-frequency decomposition\n");
		if(XLALFreqSeriesToTFPlane(tfplane, normalisation, fseries, params->useOverWhitening ? psd : NULL, rplan)) {
			XLALDestroyCOMPLEX8FrequencySeries(fseries);
			errorcode = XLAL_EFUNC;
			goto error;
		}
		XLALDestroyCOMPLEX8FrequencySeries(fseries);

		/*
		 * Compute the hrss factors
		 */

		if(hrssresponse) {
			XLALPrintInfo("XLALEPSearch(): computing the hrss factors\n");
			XLALTFPlaneEvalHrssFactor(hrssfactor, tfplane, hrssresponse, psd);
		} else
			memset(hrssfactor, 0, tfplane->params.freqBins * sizeof(*hrssfactor));

		/*
		 * Compute the excess power for each time-frequency tile
		 * using the data in the time-frequency plane.
		 */

		XLALPrintInfo("XLALEPSearch(): computing the excess power for each tile\n");
		if(XLALComputeExcessPower(Tiling, tfplane, hrssfactor, normalisation)) {
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
		head = XLALTFTilesToSnglBurstTable(head, Tiling, tseries->name, &tfplane->epoch, params);
		if(xlalErrno) {
			errorcode = XLAL_EFUNC;
			goto error;
		}
	}

	/*
	 * Memory clean-up.
	 */

	XLALPrintInfo("XLALEPSearch(): done\n");
	XLALDestroyREAL4FrequencySeries(psd);
	XLALDestroyREAL4FFTPlan(fplan);
	XLALDestroyREAL4FFTPlan(rplan);
	XLALDestroyTFPlane(tfplane);
	LALFree(normalisation);
	LALFree(hrssfactor);
	XLALDestroyTFTiling(Tiling);
	return(head);

	error:
	XLALDestroyREAL4FFTPlan(fplan);
	XLALDestroyREAL4FFTPlan(rplan);
	XLALDestroyREAL4FrequencySeries(psd);
	XLALDestroyCOMPLEX8FrequencySeries(fseries);
	XLALDestroyTFPlane(tfplane);
	LALFree(normalisation);
	LALFree(hrssfactor);
	XLALDestroyTFTiling(Tiling);
	XLALDestroySnglBurstTable(head);
	XLAL_ERROR_NULL(func, errorcode);
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
	const char *func = "XLALEPConditionData";
	const REAL8         epsilon = 1.0e-8;
	PassBandParamStruc  highpassParam;
	size_t              newlength;

	/*
	 * Resample the time series if necessary
	 */

	if(fabs(resampledeltaT - series->deltaT) >= epsilon)
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
