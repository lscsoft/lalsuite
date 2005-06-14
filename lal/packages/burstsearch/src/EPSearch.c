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

extern INT4 lalDebugLevel;

/*
 * Weight tiles according to the number present of a given time-frequency
 * volume.
 */

static INT4 DegreesOfFreedom(TFTile *tile)
{
	return(2 * (tile->tend - tile->tstart)*tile->deltaT * (tile->fend - tile->fstart)*tile->deltaF);
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


/*
 * Compute the average spectrum for the given time series
 */

static int ComputeAverageSpectrum(
	REAL4FrequencySeries *spectrum,
	REAL4TimeSeries *tseries,
	size_t windowLength,
	size_t windowShift,
	WindowType windowType,
	AvgSpecMethod method
)
{
	static const char *func = "ComputeAverageSpectrum";
	REAL4Window *window;
	RealFFTPlan *plan;
	
	window = XLALCreateREAL4Window(windowLength, windowType, 0.0);
	plan = XLALCreateForwardREAL4FFTPlan(windowLength, 0);
	if(!window || !plan)
		XLAL_ERROR(func, XLAL_EFUNC);

	switch(method) {
		case useMean:
		XLALREAL4AverageSpectrumWelch(spectrum, tseries, windowLength, windowShift, window, plan);
		break;

		case useMedian:
		XLALREAL4AverageSpectrumMedian(spectrum, tseries, windowLength, windowShift, window, plan);
		break;

		default:
		XLALDestroyREAL4Window(window);
		XLALDestroyREAL4FFTPlan(plan);
		XLAL_ERROR(func, XLAL_EINVAL);
		break;
	}

	XLALDestroyREAL4Window(window);
	XLALDestroyREAL4FFTPlan(plan);
	return(0);
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
void
LALEPSearch(
	LALStatus        *status,
	REAL4TimeSeries  *tseries,
	EPSearchParams   *params,
	SnglBurstTable  **burstEvent
)
/******** </lalVerbatim> ********/
{ 
	static const char *func = "EPSearch";
	int                       start_sample, i;
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
	const LIGOTimeGPS        gps_zero = LIGOTIMEGPSZERO;

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

	AverageSpec = XLALCreateREAL4FrequencySeries("anonymous", &gps_zero, 0, 0, &lalDimensionlessUnit, params->windowLength / 2 + 1);
	if(!AverageSpec)
		XLAL_ERROR_VOID(func, XLAL_EFUNC);
	if(ComputeAverageSpectrum(AverageSpec, tseries, params->windowLength, params->windowShift, params->windowType, params->method))
		XLAL_ERROR_VOID(func, XLAL_EFUNC);

	Psd = XLALCreateREAL4FrequencySeries("anonymous", &gps_zero, 0, 0, &lalDimensionlessUnit, params->windowLength / 2 + 1);
	if(!Psd)
		XLAL_ERROR_VOID(func, XLAL_EFUNC);
	Psd->deltaF = AverageSpec->deltaF;
	Psd->epoch = AverageSpec->epoch;
	if(params->useOverWhitening)
		for(i = 0; (unsigned) i < Psd->data->length; i++)
			Psd->data->data[i] = AverageSpec->data->data[i];
	else
		for(i = 0; (unsigned) i < Psd->data->length; i++)
			Psd->data->data[i] = 1.0;

	if(params->printSpectrum)
		print_real4fseries(Psd, "psd.dat");

	if(params->printSpectrum)
		print_real4fseries(AverageSpec, "average_spectrum.dat");

	/*
	 * Assign temporary memory for the frequency data.
	 */

	fseries = XLALCreateCOMPLEX8FrequencySeries("anonymous", &gps_zero, 0, 0, &lalDimensionlessUnit, params->windowLength / 2 + 1);
	if(!fseries)
		XLAL_ERROR_VOID(func, XLAL_EFUNC);

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
		/* Initialize the normalisation */
		normalisation = NULL;

		/*
		 * Extract a windowLength of data from the time series,
		 * compute its DFT, then free it.
		 */

		LALInfo(status->statusPtr, "Computing the DFT");
		CHECKSTATUSPTR(status);

		cutTimeSeries = XLALCutREAL4TimeSeries(tseries,  start_sample, params->windowLength);
		if(!cutTimeSeries)
			XLAL_ERROR_VOID(func, XLAL_EFUNC);

		LALComputeFrequencySeries(status->statusPtr, fseries, cutTimeSeries, dftparams);
		CHECKSTATUSPTR(status);

		/* calculate the duration for which tiles are to be created
		 * in the Single TFPlane
		 */ 
		params->tfPlaneParams.timeDuration = 2*params->windowShift*cutTimeSeries->deltaT;
		
		XLALDestroyREAL4TimeSeries(cutTimeSeries);
		
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
			LALCreateTFTiling(status->statusPtr, &tfTiling, &params->tfTilingInput, &params->tfPlaneParams);
			CHECKSTATUSPTR(status);
		}

		/*
		 * Compute the TFplanes for the data segment.
		 */

		tileStartShift = (INT4)(params->windowLength/2)-params->windowShift;
		LALInfo(status->statusPtr, "Computing the TFPlanes");
		CHECKSTATUSPTR(status);
		normalisation = LALMalloc(params->tfPlaneParams.freqBins * sizeof(REAL4));
		if(XLALComputeTFPlanes(tfTiling, fseries, tileStartShift, normalisation, Psd)) {
			XLALClearErrno();
			ABORT(status, LAL_FAIL_ERR, LAL_FAIL_MSG);
		}
	
		/*
		 * Search these planes.
		 */

		LALInfo(status->statusPtr, "Computing the excess power");
		CHECKSTATUSPTR(status);
		if(XLALComputeExcessPower(tfTiling, &params->compEPInput, normalisation))
			XLAL_ERROR_VOID(func, XLAL_EFUNC);

		/*
		 * Compute the likelihood for slightly better detection
		 * method.
		 */

#if 0
		params->lambda = XLALComputeLikelihood(tfTiling);
		if(XLALIsREAL8FailNAN(params->lambda))
			XLAL_ERROR_VOID(func, XLAL_EFUNC);
#endif

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

		WeighTFTileList(tfTiling, 10000);

		/*
		 * Convert the TFTiles into sngl_burst events for output.
		 * The threhsold cut determined by alpha is applied here
		 */

		EventAddPoint = TFTilesToSnglBurstTable(tfTiling->firstTile, EventAddPoint, &fseries->epoch, params);

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

	XLALDestroyTFTiling(tfTiling);
	XLALDestroyCOMPLEX8FrequencySeries(fseries);
	LALDestroyRealDFTParams(status->statusPtr, &dftparams);
	CHECKSTATUSPTR(status);
	XLALDestroyREAL4FrequencySeries(AverageSpec);
	XLALDestroyREAL4FrequencySeries(Psd);

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
