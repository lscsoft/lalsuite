/******** <lalVerbatim file="EPSearchCV"> ********
Author: Brady, P
Revision: $Id$
********* </lalVerbatim> ********/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
#include <lal/PrintFTSeries.h>
#include <lal/RealFFT.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/TimeSeries.h>
#include <lal/Window.h>

NRCSID(EPSEARCHC, "$Id$");

#define FALSE 0

extern INT4 lalDebugLevel;

/*
 * Weight tiles according to the number present of a given time-frequency
 * volume.
 */

static INT4 DegreesOfFreedom(TFTile *tile)
{
	return(2 * (tile->tend - tile->tstart + 1) * (tile->fend - tile->fstart + 1));
}

static void WeighTFTileList(TFTiling *tfTiling, INT4 maxDOF)
{
	TFTile *tile;
	INT4 *weight = LALCalloc(2 * maxDOF, sizeof(*weight));

	for(tile = tfTiling->firstTile; tile; tile = tile->nextTile)
		weight[DegreesOfFreedom(tile)]++;

	for(tile = tfTiling->firstTile; tile; tile = tile->nextTile)
		tile->weight = weight[DegreesOfFreedom(tile)];

	LALFree(weight);
}


static INT4 ModDegreesOfFreedom(TFTile *tile)
{
	return(2 * (tile->tend - tile->tstart)*tile->deltaT * (tile->fend - tile->fstart)*tile->deltaF);
}

static void ModWeighTFTileList(TFTiling *tfTiling, INT4 maxDOF)
{
	TFTile *tile;
	INT4 *weight = LALCalloc(2 * maxDOF, sizeof(*weight));

	for(tile = tfTiling->firstTile; tile; tile = tile->nextTile)
		weight[ModDegreesOfFreedom(tile)]++;

	for(tile = tfTiling->firstTile; tile; tile = tile->nextTile)
		tile->weight = weight[ModDegreesOfFreedom(tile)];

	LALFree(weight);
}


/*
 * Compute the average spectrum for the given time series
 */

static void ComputeAverageSpectrum(
	LALStatus *status,
	REAL4FrequencySeries *spectrum,
	REAL4TimeSeries *tseries,
	EPSearchParams *params
)
{
	LALWindowParams winParams;
	AverageSpectrumParams spec_params;

	INITSTATUS(status, "ComputeAverageSpectrum", EPSEARCHC);
	ATTATCHSTATUSPTR(status);

	memset(&spec_params, 0, sizeof(spec_params));

	winParams.type = params->windowType;
	winParams.length = params->windowLength;
	LALCreateREAL4Window(status->statusPtr, &spec_params.window, &winParams);
	CHECKSTATUSPTR(status);
	LALCreateForwardRealFFTPlan(status->statusPtr, &spec_params.plan, spec_params.window->data->length, 0);
	CHECKSTATUSPTR(status);

	spec_params.overlap = spec_params.window->data->length - params->windowShift;
	spec_params.method = params->method;

	LALREAL4AverageSpectrum(status->statusPtr, spectrum, tseries, &spec_params);
	CHECKSTATUSPTR(status);

	LALDestroyREAL4Window(status->statusPtr, &spec_params.window);
	CHECKSTATUSPTR(status);

	LALDestroyRealFFTPlan(status->statusPtr, &spec_params.plan);
	CHECKSTATUSPTR(status);

	DETATCHSTATUSPTR(status);
	RETURN(status);
}


/*
 * Convert a linked list of tiles to a linked list of burst events.
 */

static SnglBurstTable *TFTileToBurstEvent(
	TFTile *tile,
	LIGOTimeGPS *epoch,
	EPSearchParams *params  
)
{
	SnglBurstTable *event = LALMalloc(sizeof(*event));
	if(!event)
		return(NULL);

	event->next = NULL;
	strncpy(event->ifo, params->channelName, 2);
	event->ifo[2] = '\0';
	strncpy(event->search, "power", LIGOMETA_SEARCH_MAX);
	strncpy(event->channel, params->channelName, LIGOMETA_CHANNEL_MAX);

	event->start_time = *epoch;
	XLALAddFloatToGPS(&event->start_time, tile->tstart * tile->deltaT);
	event->duration = (tile->tend - tile->tstart + 1) * tile->deltaT;
	event->peak_time = event->start_time;
	XLALAddFloatToGPS(&event->peak_time, 0.5 * event->duration);
	event->bandwidth = (tile->fend - tile->fstart + 1) / tile->deltaT;
	event->central_freq = params->tfTilingInput.flow + tile->fstart/tile->deltaT + (0.5 * event->bandwidth);
	event->amplitude = tile->excessPower;
	event->snr = tile->excessPower;

	if(tile->alpha < 0)
		event->confidence =  LAL_REAL4_MAX;
	else if(tile->alpha == 0)
		event->confidence = -LAL_REAL4_MAX;
	else
		event->confidence =  log(tile->alpha);

	event->event_id = NULL;

	return(event);
}

 
static SnglBurstTable *ModTFTileToBurstEvent(
	TFTile *tile,
	LIGOTimeGPS *epoch,
	EPSearchParams *params  
)
{
	SnglBurstTable *event = LALMalloc(sizeof(*event));
	if(!event)
		return(NULL);

	event->next = NULL;
	strncpy(event->ifo, params->channelName, 2);
	event->ifo[2] = '\0';
	strncpy(event->search, "power", LIGOMETA_SEARCH_MAX);
	strncpy(event->channel, params->channelName, LIGOMETA_CHANNEL_MAX);

	event->start_time = *epoch; 

	/* Moving the epoch has to be fixed */ 
	XLALAddFloatToGPS(&event->start_time, (0.5 + (tile->tstart * tile->deltaT)));
	event->duration = (tile->tend - tile->tstart ) * tile->deltaT;
	event->peak_time = event->start_time;
	XLALAddFloatToGPS(&event->peak_time, 0.5 * event->duration);
	event->bandwidth = (tile->fend - tile->fstart ) * tile->deltaF;
	event->central_freq = params->tfTilingInput.flow + tile->fstart*tile->deltaF + (0.5 * event->bandwidth);
	event->amplitude = tile->excessPower;
	event->snr = tile->excessPower;

	if(tile->alpha < 0)
		event->confidence =  LAL_REAL4_MAX;
	else if(tile->alpha == 0)
		event->confidence = -LAL_REAL4_MAX;
	else
		event->confidence =  log(tile->alpha);

	event->event_id = NULL;

	return(event);
}


static SnglBurstTable **TFTilesToSnglBurstTable(TFTile *tile, SnglBurstTable **addpoint, LIGOTimeGPS *epoch, EPSearchParams *params)
{
	INT4 count;

	for(count = 0; tile && (tile->alpha <= params->alphaThreshold/tile->weight) && (count < params->eventLimit); tile = tile->nextTile, count++) {
		*addpoint = TFTileToBurstEvent(tile, epoch, params); 

		if(*addpoint)
			addpoint = &(*addpoint)->next;
	}

	return(addpoint);
}


static SnglBurstTable **ModTFTilesToSnglBurstTable(TFTile *tile, SnglBurstTable **addpoint, LIGOTimeGPS *epoch, EPSearchParams *params)
{
	INT4 count;

	for(count = 0; tile && (tile->alpha <= params->alphaThreshold/tile->weight) && (count < params->eventLimit); tile = tile->nextTile, count++) {
		*addpoint = ModTFTileToBurstEvent(tile, epoch, params); 

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
 * Im is 1.
 */

static void normalize_to_psd(COMPLEX8FrequencySeries *fseries, REAL4FrequencySeries *psd)
{
	REAL4 factor;
	size_t i;

	for(i = 0; i < fseries->data->length; i++) {
		factor = 2.0 * sqrt(fseries->deltaF / psd->data->data[i]);
		fseries->data->data[i].re *= factor;
		fseries->data->data[i].im *= factor;
	}
}


/*
 * Generate a linked list of burst events from a time series.
 */

/******** <lalVerbatim file="EPSearchCP"> ********/
void
EPSearch(
	LALStatus        *status,
	REAL4TimeSeries  *tseries,
	EPSearchParams   *params,
	SnglBurstTable  **burstEvent
)
/******** </lalVerbatim> ********/
{ 
	int                       start_sample, i, j;
	COMPLEX8FrequencySeries  *fseries;
	RealDFTParams            *dftparams = NULL;
	LALWindowParams           winParams;
	REAL4FrequencySeries     *AverageSpec;
	REAL4FrequencySeries     *Psd;
	REAL4TimeSeries          *cutTimeSeries;
	SnglBurstTable          **EventAddPoint = burstEvent;
	TFTiling                 *tfTiling = NULL;
	REAL4                    *normalisation;
	INT4                     tileStartShift;
	FILE                     *fp;

	INITSTATUS(status, "EPSearch", EPSEARCHC);
	ATTATCHSTATUSPTR(status);

	/*
	 * Check arguments.
	 */

	ASSERT(tseries != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	ASSERT(params != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
	ASSERT(burstEvent != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);

	/*
	 * Compute the average spectrum.
	 */

	LALCreateREAL4FrequencySeries(status->statusPtr, &AverageSpec, "anonymous", LIGOTIMEGPSINITIALIZER, 0, 0, LALUNITINITIALIZER, params->windowLength / 2 + 1);
	CHECKSTATUSPTR(status);
	ComputeAverageSpectrum(status->statusPtr, AverageSpec, tseries, params);
	CHECKSTATUSPTR(status);

	LALCreateREAL4FrequencySeries(status->statusPtr, &Psd, "anonymous", LIGOTIMEGPSINITIALIZER, 0, 0, LALUNITINITIALIZER, params->windowLength / 2 + 1);
	CHECKSTATUSPTR(status);
	Psd->deltaF = AverageSpec->deltaF;
	Psd->epoch = AverageSpec->epoch;
	if(params->useOverWhitening){
	  for(i=0; i<Psd->data->length;i++){
	    Psd->data->data[i] = AverageSpec->data->data[i];
	  }
	}
	else{
    	  for(i=0; i<Psd->data->length;i++){
	    Psd->data->data[i] = 1.0;
	  }
	}

	if(params->printSpectrum)
		print_real4fseries(Psd, "psd.dat");

	if(params->printSpectrum)
		print_real4fseries(AverageSpec, "average_spectrum.dat");

	/*
	 * Assign temporary memory for the frequency data.
	 */

	LALCreateCOMPLEX8FrequencySeries(status->statusPtr, &fseries, "anonymous", LIGOTIMEGPSINITIALIZER, 0, 0, LALUNITINITIALIZER, params->windowLength / 2 + 1);
	CHECKSTATUSPTR(status);

	/*
	 * Create the dft params.
	 */

	winParams.type = params->windowType;
	winParams.length = params->windowLength;
	LALCreateRealDFTParams(status->statusPtr , &dftparams, &winParams, 1);
	CHECKSTATUSPTR(status);

	/*
	 * Loop over data applying excess power method.
	 */

	for(start_sample = 0; start_sample + params->windowLength <= tseries->data->length; start_sample += params->windowShift) {
		LALInfo(status->statusPtr, "Analyzing a window...");
		CHECKSTATUSPTR(status);
		/*
		 * Extract a windowLength of data from the time series,
		 * compute its DFT, then free it.
		 */

		LALInfo(status->statusPtr, "Computing the DFT");
		CHECKSTATUSPTR(status);

		LALCutREAL4TimeSeries(status->statusPtr, &cutTimeSeries, tseries,  start_sample, params->windowLength);
		CHECKSTATUSPTR(status);

		LALComputeFrequencySeries(status->statusPtr, fseries, cutTimeSeries, dftparams);
		CHECKSTATUSPTR(status);

		/* calculate the duration for which tiles are to be created
		 * in the Single TFPlane
		 */ 
		if (params->tfPlaneMethod == useSingleTFPlane)
		  params->tfPlaneParams.timeDuration = 2*params->windowShift*cutTimeSeries->deltaT;
		
		LALDestroyREAL4TimeSeries(status->statusPtr, cutTimeSeries);
		CHECKSTATUSPTR(status);
		
		/*
		 * Normalize the frequency series to the average PSD.
		 */

		normalize_to_psd(fseries, AverageSpec);

		if(params->printSpectrum)
			print_complex8fseries(fseries, "frequency_series.dat");

		/*
		 * Create time-frequency tiling of plane.
		 */

		if(!tfTiling) {
			/* the factor of 2 here is to account for the
			 * overlapping */
			params->tfTilingInput.deltaF = 2.0 * fseries->deltaF;
			if (params->tfPlaneMethod == useMultipleTFPlane)
			  LALCreateTFTiling(status->statusPtr, &tfTiling, &params->tfTilingInput);
			else
			  LALModCreateTFTiling(status->statusPtr, &tfTiling, &params->tfTilingInput, &params->tfPlaneParams);
			CHECKSTATUSPTR(status);
		}

		/*
		 * Compute the TFplanes for the data segment.
		 */
		tileStartShift = (INT4)(params->windowLength/2)-params->windowShift;
		LALInfo(status->statusPtr, "Computing the TFPlanes");
		CHECKSTATUSPTR(status);
		if (params->tfPlaneMethod == useMultipleTFPlane)
		  LALComputeTFPlanes(status->statusPtr, tfTiling, fseries, tileStartShift);
		else{
		  normalisation = LALMalloc(params->tfPlaneParams.freqBins * sizeof(REAL4));
		  LALModComputeTFPlanes(status->statusPtr, tfTiling, fseries, tileStartShift, normalisation, Psd);
		}
		CHECKSTATUSPTR(status);
	
                 /*
		 * Search these planes.
		 */
		LALInfo(status->statusPtr, "Computing the excess power");
		CHECKSTATUSPTR(status);
		if (params->tfPlaneMethod == useMultipleTFPlane)
		  LALComputeExcessPower (status->statusPtr, tfTiling, &params->compEPInput);
		else
		  LALModComputeExcessPower (status->statusPtr, tfTiling, &params->compEPInput, normalisation);
		CHECKSTATUSPTR(status);

		/*
		 * Compute the likelihood for slightly better detection
		 * method.
		 */

		/*LALComputeLikelihood(status->statusPtr, &(params->lambda),
		tfTiling); CHECKSTATUSPTR(status);*/

		/*
		 * Sort the results.
		 */

		LALInfo(status->statusPtr, "Sorting TFTiling");
		CHECKSTATUSPTR(status);
		LALSortTFTiling(status->statusPtr, tfTiling);
		CHECKSTATUSPTR(status);

		/*
		 * Determine the weighting for each tile.
		 */
		if (params->tfPlaneMethod == useMultipleTFPlane)
		  WeighTFTileList(tfTiling, 10000);
		else 
		  ModWeighTFTileList(tfTiling, 10000);

		/*
		 * Convert the TFTiles into sngl_burst events for output.
		 * The threhsold cut determined by alpha is applied here
		 */
		if (params->tfPlaneMethod == useMultipleTFPlane)
		  EventAddPoint = TFTilesToSnglBurstTable(tfTiling->firstTile, EventAddPoint, &fseries->epoch, params);
		else
		  EventAddPoint = ModTFTilesToSnglBurstTable(tfTiling->firstTile, EventAddPoint, &fseries->epoch, params);

		/*
		 * Reset the flags on the tftiles.
		 */

		tfTiling->planesComputed=FALSE;
		tfTiling->excessPowerComputed=FALSE;
		tfTiling->tilesSorted=FALSE;
		LALFree(normalisation);
	}

	/*
	 * Memory clean-up.
	 */

	LALDestroyTFTiling(status->statusPtr, &tfTiling);
	CHECKSTATUSPTR(status);
	LALDestroyCOMPLEX8FrequencySeries(status->statusPtr, fseries);
	CHECKSTATUSPTR(status);
	LALDestroyRealDFTParams(status->statusPtr, &dftparams);
	CHECKSTATUSPTR(status);
	LALDestroyREAL4FrequencySeries(status->statusPtr, AverageSpec);
	CHECKSTATUSPTR(status);
	LALDestroyREAL4FrequencySeries(status->statusPtr, Psd);
	CHECKSTATUSPTR(status);
	/*
	 * Normal exit.
	 */

	DETATCHSTATUSPTR(status);
	RETURN(status);
}


/*
 * Condition the time series prior to analysis by the power code
 */

/* <lalVerbatim file="EPConditionDataCP"> */
void EPConditionData(
	LALStatus        *status,
	REAL4TimeSeries  *series,
	REAL4             flow,
	REAL8             resampledeltaT,
	ResampleTSFilter  resampFiltType,
	INT4              corruption
)
/* </lalVerbatim> */
{
	ResampleTSParams    resampleParams;
	const REAL8         epsilon = 1.0e-8;
	PassBandParamStruc  highpassParam;

	INITSTATUS (status, "LALConditionData", EPSEARCHC);
	ATTATCHSTATUSPTR (status);

	ASSERT(series != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);

	/*
	 * Resample the time series if necessary
	 */

	memset(&resampleParams, 0, sizeof(resampleParams));
	if(!(fabs(resampledeltaT - series->deltaT) < epsilon)) {
		resampleParams.deltaT = resampledeltaT;
		resampleParams.filterType = resampFiltType;
		LALResampleREAL4TimeSeries(status->statusPtr, series, &resampleParams);
		CHECKSTATUSPTR (status);
	}

	/*
	 * High-pass filter the time series.
	 */

	highpassParam.nMax = 8;
	highpassParam.f2 = (flow - 10.0) > 150.0 ? 150.0 : (flow - 10.0);
	highpassParam.f1 = -1.0;
	highpassParam.a2 = 0.9;
	highpassParam.a1 = -1.0;
	LALButterworthREAL4TimeSeries(status->statusPtr, series, &highpassParam);
	CHECKSTATUSPTR (status);

	/*
	 * The filter corrupts the ends of the time series.  Chop them off.
	 */

	LALShrinkREAL4TimeSeries(status->statusPtr, series, corruption, series->data->length - 2 * corruption);
	CHECKSTATUSPTR (status);

	DETATCHSTATUSPTR(status);
	RETURN(status);
}
