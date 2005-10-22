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
#include <lal/LALStdlib.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/RealFFT.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/LALStatusMacros.h>
#include <lal/TimeFreqFFT.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/Window.h>
#include <lal/XLALError.h>

NRCSID(EPSEARCHC, "$Id$");


/*
 * Convert an array of tiles to a linked list of burst events.
 */
 
static SnglBurstTable *TFTileToBurstEvent(
	const TFTile *tile,
	const char *channelName,
	const LIGOTimeGPS *epoch,
	const EPSearchParams *params  
)
{
	SnglBurstTable *event = LALMalloc(sizeof(*event));
	if(!event)
		return(NULL);

	event->next = NULL;
	strncpy(event->ifo, channelName, 2);
	event->ifo[2] = '\0';
	strncpy(event->search, "power", LIGOMETA_SEARCH_MAX);
	event->search[LIGOMETA_SEARCH_MAX] = '\0';
	strncpy(event->channel, channelName, LIGOMETA_CHANNEL_MAX);
	event->channel[LIGOMETA_CHANNEL_MAX] = '\0';

	event->start_time = *epoch; 
 
	XLALAddFloatToGPS(&event->start_time, (tile->tstart * tile->deltaT));
	event->duration = tile->tbins * tile->deltaT;
	event->peak_time = event->start_time;
	XLALAddFloatToGPS(&event->peak_time, 0.5 * event->duration);
	event->bandwidth = tile->fbins * tile->deltaF;
	event->central_freq = params->tfPlaneParams.flow + tile->fstart*tile->deltaF + (0.5 * event->bandwidth);
	event->amplitude = tile->hrss;
	event->snr = tile->excessPower;
	event->confidence =  tile->lnalpha;

	event->event_id = NULL;

	return(event);
}


static SnglBurstTable **TFTilesToSnglBurstTable(const TFTiling *tiling, const char *channelName, SnglBurstTable **addpoint, const LIGOTimeGPS *epoch, const EPSearchParams *params)
{
	TFTile *tile;
	size_t i;

	for(i = 0, tile = tiling->tile; (i < tiling->numtiles) && (tile->lnalpha <= params->lnalphaThreshold - tile->lnweight); i++, tile++) {
		*addpoint = TFTileToBurstEvent(tile, channelName, epoch, params); 

		if(*addpoint)
			addpoint = &(*addpoint)->next;
	}

	return(addpoint);
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
int
XLALEPSearch(
	SnglBurstTable  **burstEvent,
	const COMPLEX8FrequencySeries  *hrssresponse,
	const REAL4TimeSeries  *tseries,
	EPSearchParams   *params
)
/******** </lalVerbatim> ********/
{ 
	static const char *func = "EPSearch";
	int errorcode;
	int                      start_sample;
	int                      overwhiten_flag = 0; /* default */
	COMPLEX8FrequencySeries *fseries = NULL;
	const COMPLEX8FrequencySeries *response;
	REAL4Window             *window = params->window;
	RealFFTPlan             *plan;
	REAL4FrequencySeries    *AverageSpec;
	REAL4FrequencySeries    *PsdSpec;
	REAL4TimeSeries         *cutTimeSeries;
	SnglBurstTable         **EventAddPoint = burstEvent;
	TFTiling                *Tiling;
	REAL4TimeFrequencyPlane *tfplane;
	REAL4                   *normalisation;
	REAL8                   *hrssfactor;
	REAL8                   *fachrss;
	const LIGOTimeGPS        gps_zero = LIGOTIMEGPSZERO;

	/*
	 * Create an FFT plan, allocate space for the average spectrum,
	 * allocate temporary storage for frequency series data, allocate
	 * storage for the normalisation data, allocate and initialize the
	 * time-frequency plane storage, and construct a time-frequency
	 * tiling of the plane.
	 */

	plan = XLALCreateForwardREAL4FFTPlan(window->data->length, 0);
	AverageSpec = XLALCreateREAL4FrequencySeries("anonymous", &gps_zero, 0, 0, &lalDimensionlessUnit, window->data->length / 2 + 1);
	tfplane = XLALCreateTFPlane(&params->tfPlaneParams);
	normalisation = LALMalloc(params->tfPlaneParams.freqBins * sizeof(*normalisation));
	hrssfactor = LALMalloc(params->tfPlaneParams.freqBins * sizeof(*hrssfactor));
	Tiling = XLALCreateTFTiling(&params->tfTilingInput, &tfplane->params);

	if(!normalisation || !hrssfactor) {
		errorcode = XLAL_ENOMEM;
		goto error;
	}
	if(!plan || !AverageSpec || !tfplane || !Tiling) {
		errorcode = XLAL_EFUNC;
		goto error;
	}
	if(!Tiling->numtiles) {
		/* couldn't fit any tiles into the TF plane! */
		errorcode = XLAL_EINVAL;
		goto error;
	}

	if(params->useOverWhitening || hrssresponse)
		PsdSpec = AverageSpec;
	else 
		PsdSpec = NULL;
	

	if(params->useOverWhitening)
		overwhiten_flag = 1;

	/*
	 * Compute the average spectrum.
	 */

	switch(params->method) {
		case useMean:
		XLALREAL4AverageSpectrumWelch(AverageSpec, tseries, window->data->length, params->windowShift, window, plan);
		break;

		case useMedian:
		XLALREAL4AverageSpectrumMedian(AverageSpec, tseries, window->data->length, params->windowShift, window, plan);
		break;

		default:
		errorcode = XLAL_EINVAL;
		goto error;
	}

	if(params->printSpectrum)
		print_real4fseries(AverageSpec, "average_spectrum.dat");

	/*
	 * Loop over data applying excess power method.
	 */

	for(start_sample = 0; start_sample + window->data->length <= tseries->data->length; start_sample += params->windowShift) {
		/*
		 * Extract a window-length of data from the time series,
		 * compute its DFT, then free it.
		 */

		cutTimeSeries = XLALCutREAL4TimeSeries(tseries, start_sample, window->data->length);
		if(!cutTimeSeries) {
			errorcode = XLAL_EFUNC;
			goto error;
		}
		XLALPrintInfo("XLALEPSearch(): analyzing samples %zu -- %zu (%.9lf s -- %.9lf s)\n", start_sample, start_sample + cutTimeSeries->data->length, start_sample * cutTimeSeries->deltaT, (start_sample + cutTimeSeries->data->length) * cutTimeSeries->deltaT);

		XLALPrintInfo("XLALEPSearch(): computing the Fourier transform\n");
		fseries = XLALComputeFrequencySeries(cutTimeSeries, window, plan);
		if(!fseries) {
			errorcode = XLAL_EFUNC;
			goto error;
		}
		XLALDestroyREAL4TimeSeries(cutTimeSeries);

		/* tweak the frequency series' epoch */
		XLALAddFloatToGPS(&fseries->epoch, params->windowShift * tseries->deltaT);
		
		/*
		 * Normalize the frequency series to the average PSD.
		 */

		XLALPrintInfo("XLALEPSearch(): normalizing to the average spectrum\n");
		whiten(fseries, AverageSpec);
		if(params->printSpectrum)
			print_complex8fseries(fseries, "frequency_series.dat");

		/*
		 * Apply the phase factor of the reponse here, if estimating h_rss
		 * FIXME: Now we are not multiplying by the phase, will do after
		 * implementing the rest.
		 */

		if (hrssresponse) {
			fachrss = hrssfactor;
			response = hrssresponse;
		}
		else {
			fachrss = NULL;
			response = NULL;
		}

		/*
		 * Compute the time-frequency plane from the frequency
		 * series.
		 */

		XLALPrintInfo("XLALEPSearch(): computing the time-frequency decomposition\n");
		if(XLALFreqSeriesToTFPlane(tfplane, fseries, window->data->length / 2 - params->windowShift, fachrss, normalisation, response, PsdSpec, overwhiten_flag)) {
			XLALDestroyCOMPLEX8FrequencySeries(fseries);
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
		if(XLALComputeExcessPower(Tiling, tfplane, fachrss, normalisation)) {
			errorcode = XLAL_EFUNC;
			goto error;
		}

		/*
		 * Compute the likelihood for slightly better detection
		 * method.
		 */

#if 0
		params->lambda = XLALComputeLikelihood(Tiling);
#endif

		/*
		 * Sort the results.
		 */

		XLALPrintInfo("XLALEPSearch(): sorting the tiles\n");
		XLALSortTFTilingByAlpha(Tiling);

		/*
		 * Convert the TFTiles into sngl_burst events for output.
		 * The threhsold cut determined by alpha is applied here
		 */

		EventAddPoint = TFTilesToSnglBurstTable(Tiling, tseries->name, EventAddPoint, &tfplane->epoch, params);
	}

	/*
	 * Memory clean-up.
	 */

	XLALPrintInfo("XLALEPSearch(): done\n");
	XLALDestroyREAL4FrequencySeries(AverageSpec);
	XLALDestroyREAL4FFTPlan(plan);
	XLALDestroyTFPlane(tfplane);
	LALFree(normalisation);
	LALFree(hrssfactor);
	XLALDestroyTFTiling(Tiling);
	return(0);

	error:
	XLALDestroyREAL4FFTPlan(plan);
	XLALDestroyREAL4FrequencySeries(AverageSpec);
	XLALDestroyCOMPLEX8FrequencySeries(fseries);
	XLALDestroyTFPlane(tfplane);
	LALFree(normalisation);
	LALFree(hrssfactor);
	XLALDestroyTFTiling(Tiling);
	XLAL_ERROR(func, errorcode);
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
