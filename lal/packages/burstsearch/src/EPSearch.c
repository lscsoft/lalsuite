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
#include <lal/Units.h>
#include <lal/Window.h>
#include <lal/XLALError.h>

NRCSID(EPSEARCHC, "$Id$");

#define FALSE 0


/*
 * Weight tiles according to the number present of a given time-frequency
 * volume.
 */

static INT4 DegreesOfFreedom(TFTile *tile)
{
	return(2 * (tile->tend - tile->tstart)*tile->deltaT * (tile->fend - tile->fstart)*tile->deltaF);
}

static void WeightTFTileList(TFTiling *tfTiling, INT4 maxDOF)
{
	TFTile *tile;
	int *weight = LALCalloc(maxDOF, sizeof(*weight));

	for(tile = tfTiling->firstTile; tile; tile = tile->nextTile)
		weight[DegreesOfFreedom(tile)]++;

	for(tile = tfTiling->firstTile; tile; tile = tile->nextTile)
		tile->weight = weight[DegreesOfFreedom(tile)];

	LALFree(weight);
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

	/* Moving the epoch has to be fixed */ 
	XLALAddFloatToGPS(&event->start_time, (tile->tstart * tile->deltaT));
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
		/* FIXME: it is an error for this to occur in the band of
		 * interest.  The entire PSD window should be thrown out */
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
	REAL4TimeSeries  *tseries,
	EPSearchParams   *params
)
/******** </lalVerbatim> ********/
{ 
	static const char *func = "EPSearch";
	int                       start_sample;
	COMPLEX8FrequencySeries  *fseries;
	REAL4Window              *window;
	RealFFTPlan              *plan;
	REAL4FrequencySeries     *AverageSpec;
	REAL4FrequencySeries     *OverWhiteningSpec;
	REAL4TimeSeries          *cutTimeSeries;
	SnglBurstTable          **EventAddPoint = burstEvent;
	TFTiling                 *tfTiling = NULL;
	REAL4                    *normalisation;
	const LIGOTimeGPS        gps_zero = LIGOTIMEGPSZERO;

	/*
	 * Create a time-domain window, an FFT plan, and allocate space for
	 * the average spectrum, and temporary frequency series data.
	 */

	window = XLALCreateREAL4Window(params->windowLength, params->windowType, 0.0);
	plan = XLALCreateForwardREAL4FFTPlan(window->data->length, 0);
	AverageSpec = XLALCreateREAL4FrequencySeries("anonymous", &gps_zero, 0, 0, &lalDimensionlessUnit, params->windowLength / 2 + 1);
	fseries = XLALCreateCOMPLEX8FrequencySeries("anonymous", &gps_zero, 0, 0, &lalDimensionlessUnit, params->windowLength / 2 + 1);
	if(!window || !plan || !AverageSpec || !fseries) {
		XLALDestroyREAL4Window(window);
		XLALDestroyREAL4FFTPlan(plan);
		XLALDestroyREAL4FrequencySeries(AverageSpec);
		XLALDestroyCOMPLEX8FrequencySeries(fseries);
		XLAL_ERROR(func, XLAL_EFUNC);
	}

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
		XLAL_ERROR(func, XLAL_EINVAL);
		break;
	}

	if(params->printSpectrum)
		print_real4fseries(AverageSpec, "average_spectrum.dat");

	if(params->useOverWhitening)
		OverWhiteningSpec = AverageSpec;
	else
		OverWhiteningSpec = NULL;

	/*
	 * Loop over data applying excess power method.
	 */

	for(start_sample = 0; start_sample + params->windowLength <= tseries->data->length; start_sample += params->windowShift) {
		XLALPrintInfo("Analyzing a window...");

		/*
		 * Calculate the duration for which tiles are to be created
		 * in the Single TFPlane.
		 */ 

		params->tfPlaneParams.timeDuration = 2.0 * params->windowShift * tseries->deltaT;

		/*
		 * Extract a windowLength of data from the time series,
		 * compute its DFT, then free it.
		 */

		XLALPrintInfo("Computing the DFT");
		cutTimeSeries = XLALCutREAL4TimeSeries(tseries,  start_sample, params->windowLength);
		if(!cutTimeSeries)
			XLAL_ERROR(func, XLAL_EFUNC);
		if(XLALComputeFrequencySeries(fseries, cutTimeSeries, window, plan))
			XLAL_ERROR(func, XLAL_EFUNC);
		XLALDestroyREAL4TimeSeries(cutTimeSeries);
		
		/*
		 * Normalize the frequency series to the average PSD.
		 */

		normalize_to_psd(fseries, AverageSpec);
		if(params->printSpectrum)
			print_complex8fseries(fseries, "frequency_series.dat");

		/*
		 * Create time-frequency tiling of plane.  The factor of 2
		 * here is to account for the overlapping
		 */

		params->tfTilingInput.deltaF = 2.0 * fseries->deltaF;
		tfTiling = XLALCreateTFTiling(&params->tfTilingInput, &params->tfPlaneParams);
		if(!tfTiling)
			XLAL_ERROR(func, XLAL_EFUNC);

		/*
		 * Compute the TFplanes for the data segment.
		 */

		XLALPrintInfo("Computing the TFPlanes");
		normalisation = LALMalloc(params->tfPlaneParams.freqBins * sizeof(*normalisation));
		if(XLALComputeTFPlanes(tfTiling, fseries, params->windowLength / 2 - params->windowShift, normalisation, OverWhiteningSpec))
			XLAL_ERROR(func, XLAL_EFUNC);
	
		/*
		 * Search these planes.
		 */

		XLALPrintInfo("Computing the excess power");
		if(XLALComputeExcessPower(tfTiling, &params->compEPInput, normalisation))
			XLAL_ERROR(func, XLAL_EFUNC);
		LALFree(normalisation);
		normalisation = NULL;

		/*
		 * Compute the likelihood for slightly better detection
		 * method.
		 */

#if 0
		params->lambda = XLALComputeLikelihood(tfTiling);
		if(XLALIsREAL8FailNAN(params->lambda))
			XLAL_ERROR(func, XLAL_EFUNC);
#endif

		/*
		 * Sort the results.
		 */

		XLALPrintInfo("Sorting TFTiling");
		if(XLALSortTFTiling(tfTiling))
			XLAL_ERROR(func, XLAL_EFUNC);

		/*
		 * Determine the weighting for each tile.
		 */

		WeightTFTileList(tfTiling, 20000);

		/*
		 * Convert the TFTiles into sngl_burst events for output.
		 * The threhsold cut determined by alpha is applied here
		 */

		EventAddPoint = TFTilesToSnglBurstTable(tfTiling->firstTile, EventAddPoint, &fseries->epoch, params);

		/*
		 * Clean up
		 */

		XLALDestroyTFTiling(tfTiling);
	}

	/*
	 * Memory clean-up.
	 */

	XLALDestroyCOMPLEX8FrequencySeries(fseries);
	XLALDestroyREAL4FrequencySeries(AverageSpec);
	XLALDestroyREAL4Window(window);
	XLALDestroyREAL4FFTPlan(plan);

	return(0);
}


/*
 * Condition the time series prior to analysis by the power code
 */

/* <lalVerbatim file="EPConditionDataCP"> */
void LALEPConditionData(
	LALStatus        *status,
	REAL4TimeSeries  *series,
	REAL8             flow,
	REAL8             resampledeltaT,
	ResampleTSFilter  resampFiltType,
	INT4              corruption
)
/* </lalVerbatim> */
{
	ResampleTSParams    resampleParams;
	const REAL8         epsilon = 1.0e-8;
	PassBandParamStruc  highpassParam;
	size_t              newlength;

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
	highpassParam.f2 = flow;
	highpassParam.f1 = -1.0;
	highpassParam.a2 = 0.9;
	highpassParam.a1 = -1.0;
	LALButterworthREAL4TimeSeries(status->statusPtr, series, &highpassParam);
	CHECKSTATUSPTR (status);

	/*
	 * The filter corrupts the ends of the time series.  Chop them off.
	 */

	newlength = series->data->length - 2 * corruption;
	ASSERT(XLALShrinkREAL4TimeSeries(series, corruption, newlength) == newlength, status, LAL_FAIL_ERR, LAL_FAIL_MSG);

	DETATCHSTATUSPTR(status);
	RETURN(status);
}
