/******** <lalVerbatim file="EPSearchCV"> ********
Author: Brady, P
$Id$
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
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALErrno.h>
#include <lal/LALRCSID.h>
#include <lal/LALStdlib.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/PrintFTSeries.h>
#include <lal/Random.h>
#include <lal/RealFFT.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/SeqFactories.h>
#include <lal/Thresholds.h>
#include <lal/TimeFreqFFT.h>
#include <lal/TimeSeries.h>
#include <lal/Window.h>

NRCSID (EPSEARCHC, "$Id$");

#define TRUE 1
#define FALSE 0

extern INT4 lalDebugLevel;

/****************************************************************
 *
 * Weights tiles according to the number present of a given
 * time-frequency volume.
 *
 ***************************************************************/

static INT4 DegreesOfFreedom(TFTile *tile)
{
	return(2 * (tile->tend - tile->tstart + 1) * (tile->fend - tile->fstart + 1));
}

static void WeighTFTileList (
        TFTiling          *tfTiling,
        INT4               maxDOF
        )
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

static void ComputeAverageSpectrum(
	LALStatus *status,
	REAL4FrequencySeries *spectrum,
	REAL4TimeSeries *tseries,
	EPSearchParams *params
)
{
	int i;
	LALWindowParams winParams;
	AverageSpectrumParams spec_params = {};

	INITSTATUS(status, "ComputeAverageSpectrum", EPSEARCHC);
	ATTATCHSTATUSPTR(status);

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

	/* Adjust the normalization to agree with our old implementation.
	 * Eventually, the rest of the code will be migrated to the LAL
	 * conventions, and this will not be needed. */
	for(i = 0; i < spectrum->data->length; i++)
		spectrum->data->data[i] /= 2 * tseries->deltaT;

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

	event->next = NULL;
	strncpy(event->ifo, params->channelName, 2);
	event->ifo[2] = '\0';
	strncpy(event->search, "power", LIGOMETA_SEARCH_MAX);
	strncpy(event->channel, params->channelName, LIGOMETA_CHANNEL_MAX);

	event->start_time = XLALAddFloatToGPS(*epoch, tile->tstart * tile->deltaT);
	event->duration = (tile->tend - tile->tstart + 1) * tile->deltaT;
	event->peak_time = XLALAddFloatToGPS(event->start_time, 0.5 * event->duration);
	event->bandwidth = (tile->fend - tile->fstart + 1) / tile->deltaT;
	event->central_freq = params->tfTilingInput->flow + tile->fstart/tile->deltaT + (0.5 * event->bandwidth);
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


static SnglBurstTable **TFTilesToSnglBurstTable(LALStatus *status, TFTile *tile, SnglBurstTable **event, LIGOTimeGPS *epoch, EPSearchParams *params)
{
	INT4 numevents;

	LALInfo(status->statusPtr, "Converting times into sngl_burst events");
	CHECKSTATUSPTR (status);

	for(numevents = 0; tile && (tile->alpha <= params->alphaThreshold/tile->weight) && (numevents < params->eventLimit); tile = tile->nextTile, numevents++) {
		*event = TFTileToBurstEvent(tile, epoch, params); 

		event = &(*event)->next;
	}

	return(event);
}


/******** <lalVerbatim file="EPSearchCP"> ********/
void
EPSearch (
        LALStatus               *status,
        REAL4TimeSeries         *tseries,
        EPSearchParams          *params,
        SnglBurstTable         **burstEvent
        )
/******** </lalVerbatim> ********/
{ 
    int                       start_sample;
    INT4                      j,freqcl;
    COMPLEX8FrequencySeries  *fseries;
    RealDFTParams            *dftparams = NULL;
    LALWindowParams           winParams;
    REAL4FrequencySeries     *AverageSpec;
    REAL4TimeSeries          *cutTimeSeries;
    INT4                      nevents, dumevents;
    SnglBurstTable          **EventAddPoint = burstEvent;
    TFTiling                 *tfTiling = NULL;

    INITSTATUS(status, "EPSearch", EPSEARCHC);
    ATTATCHSTATUSPTR(status);

    /* make sure that arguments are not NULL */
    ASSERT(tseries != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
    ASSERT(params != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);
    ASSERT(burstEvent != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);

    /* Compute the average spectrum */
    LALCreateREAL4FrequencySeries(status->statusPtr, &AverageSpec, "anonymous", LIGOTIMEGPSINITIALIZER, 0, 0, LALUNITINITIALIZER, params->windowLength / 2 + 1);
    CHECKSTATUSPTR(status);
    ComputeAverageSpectrum(status->statusPtr, AverageSpec, tseries, params);
    CHECKSTATUSPTR(status);
    if(params->printSpectrum) { 
        FILE *fp = fopen("average_spectrum.dat","w");
        for(j = 0; j < (INT4)AverageSpec->data->length; j++)
            fprintf(fp, "%f\t%g\n", j*AverageSpec->deltaF, AverageSpec->data->data[j]);
        fclose(fp);
    }

    /* assign temporary memory for the frequency data */
    LALCreateCOMPLEX8FrequencySeries(status->statusPtr, &fseries, "anonymous", LIGOTIMEGPSINITIALIZER, 0, 0, LALUNITINITIALIZER, params->windowLength / 2 + 1);
    CHECKSTATUSPTR(status);

    /* create the dft params */
    winParams.type = params->windowType;
    winParams.length = params->windowLength;
    LALCreateRealDFTParams(status->statusPtr , &dftparams, &winParams, 1);
    CHECKSTATUSPTR(status);

    /* loop over data applying excess power method */
    for(start_sample = 0; start_sample + params->windowLength <= tseries->data->length; start_sample += params->windowShift) {
      /* extract a windowLength of data from the time series */
      LALCutREAL4TimeSeries(status->statusPtr, &cutTimeSeries, tseries,  start_sample, params->windowLength);
      CHECKSTATUSPTR(status);

      /* compute its DFT */
      LALInfo(status->statusPtr, "Computing the frequency series");
      CHECKSTATUSPTR(status);
      LALComputeFrequencySeries(status->statusPtr, fseries, cutTimeSeries, dftparams);
      CHECKSTATUSPTR(status);
  
      /* delete it */
      LALDestroyREAL4TimeSeries(status->statusPtr, cutTimeSeries);
      CHECKSTATUSPTR(status);

      /* normalize the spectrum so that rms of Re or Im is 1 */
      for(j=0 ; j<(INT4)fseries->data->length ; j++) {
        REAL4 tmpVar;
        tmpVar = sqrt( 2.0 / AverageSpec->data->data[j] );
        fseries->data->data[j].re *= tmpVar;
        fseries->data->data[j].im *= tmpVar;
      }

      /* write diagnostic info to disk */
      if(params->printSpectrum) { 
        FILE *fp = fopen("frequency_series.dat","w");
        for(j = 0; j < (INT4)fseries->data->length; j++)
          fprintf(fp, "%f\t%g\n", j*fseries->deltaF, sqrt(fseries->data->data[j].re * fseries->data->data[j].re + fseries->data->data[j].im * fseries->data->data[j].im));
        fclose(fp);
      }

      /* create time-frequency tiling of plane.  */
      if(!tfTiling) {
        /* the factor of 2 here is to account for the overlapping */
        params->tfTilingInput->deltaF = 2.0 * fseries->deltaF;
        LALCreateTFTiling(status->statusPtr, &tfTiling, params->tfTilingInput);
        CHECKSTATUSPTR(status);
      }

      /* compute the TFplanes for the data segment */
      LALInfo(status->statusPtr, "Computing the TFPlanes");
      CHECKSTATUSPTR(status);
      LALComputeTFPlanes(status->statusPtr, tfTiling, fseries);
      CHECKSTATUSPTR(status);

      /* search these planes */
      LALInfo(status->statusPtr, "Computing the excess power");
      CHECKSTATUSPTR(status);
      LALComputeExcessPower (status->statusPtr, tfTiling, params->compEPInput);
      CHECKSTATUSPTR(status);

      /* compute the likelihood for slightly better detection method */
      /*
       * LALComputeLikelihood(status->statusPtr, &(params->lambda), tfTiling);
       * CHECKSTATUSPTR(status);
       */

      /* sort the results. */
      LALInfo(status->statusPtr, "Sorting TFTiling");
      CHECKSTATUSPTR(status);
      LALSortTFTiling(status->statusPtr, tfTiling);
      CHECKSTATUSPTR(status);

      /* determine the weighting for each tile */
      WeighTFTileList(tfTiling, 10000);

      /* convert the TFTiles into sngl_burst events for output */
      EventAddPoint = TFTilesToSnglBurstTable(status, tfTiling->firstTile, EventAddPoint, &fseries->epoch, params);

      /* reset the flags on the tftiles */
      tfTiling->planesComputed=FALSE;
      tfTiling->excessPowerComputed=FALSE;
      tfTiling->tilesSorted=FALSE;
    }

    /* cluster the events if requested */
    if(params->cluster && *burstEvent) {
      j = nevents = 0;
      do {
        LALSortSnglBurst(status->statusPtr, burstEvent, LALCompareSnglBurstByTimeAndFreq);
        CHECKSTATUSPTR(status);
        dumevents = nevents;
        LALClusterSnglBurstTable(status->statusPtr, *burstEvent, &nevents);
        CHECKSTATUSPTR(status);
      } while((dumevents != nevents) && (++j < 500));
    }

    /* memory clean-up */
    LALDestroyTFTiling(status->statusPtr, &tfTiling);
    CHECKSTATUSPTR(status);
    LALDestroyCOMPLEX8FrequencySeries(status->statusPtr, fseries);
    CHECKSTATUSPTR(status);
    LALDestroyRealDFTParams(status->statusPtr, &dftparams);
    CHECKSTATUSPTR(status);
    LALDestroyREAL4FrequencySeries(status->statusPtr, AverageSpec);
    CHECKSTATUSPTR(status);

    /* normal exit */
    DETATCHSTATUSPTR(status);
    RETURN(status);
}


/********************************************************************
 *
 * This function conditions the time series prior to analysis
 * by the power code
 *
 ********************************************************************/
/* <lalVerbatim file="EPConditionDataCP"> */
void EPConditionData(
        LALStatus             *status,
        REAL4TimeSeries       *series,
	REAL4                  flow,
	REAL8                  resampledeltaT,
	ResampleTSFilter       resampFiltType,
        EPSearchParams        *params
        )
/* </lalVerbatim> */
{
    UINT4                         i,j;
    INT8                          dataTimeNS  = 0;
    REAL4                        *dummyData    = NULL;  
    REAL4                         fsafety=0;
    PassBandParamStruc            highpassParam;

    /* resample parameters */
    UINT4                  resampleChan = 0;
    ResampleTSParams       resampleParams;
    const REAL8            epsilon = 1.0e-8;

    INITSTATUS (status, "LALConditionData", EPSEARCHC);
    ATTATCHSTATUSPTR (status);

    /****************************************************************
     *
     * identify multiDimData sequences by their names 
     *
     ***************************************************************/

    ASSERT(params != NULL, status, LAL_NULL_ERR, LAL_NULL_MSG);

    /* **************************************************************
     *
     * Resample the data if necessary
     *
     * *************************************************************/

    /* set the time series parameters of the input data and resample params */
    memset( &resampleParams, 0, sizeof(ResampleTSParams) );
    resampleParams.deltaT = resampledeltaT;

    if( ! ( fabs( resampleParams.deltaT - series->deltaT ) < epsilon ) )
      {
	resampleChan = 1;
	
	if ( resampFiltType == 0 )
	  {
	    resampleParams.filterType = LDASfirLP;
	  }
	else if ( resampFiltType == 1 )
	  {
	    resampleParams.filterType = defaultButterworth;
	  }
      }
    
    /*resample if required to */
    if ( resampleChan )
      {
	LALResampleREAL4TimeSeries(status->statusPtr, series, &resampleParams);
	CHECKSTATUSPTR (status);
      }

    /* ***************************************************************
     *
     * translate InPut to DataSegment and allocate memory for response
     * function data
     *
     ****************************************************************/

    /* Set up for a highpass filter */
    highpassParam.nMax = 4;
    fsafety = params->tfTilingInput->flow - 10.0;
    highpassParam.f2 = fsafety > 150.0 ? 150.0 : fsafety;
    highpassParam.f1 = -1.0;
    highpassParam.a2 = 0.1;
    highpassParam.a1 = -1.0;
    LALButterworthREAL4TimeSeries(status->statusPtr, series, &highpassParam);
    CHECKSTATUSPTR (status);
            
    LALShrinkREAL4TimeSeries(status->statusPtr, series, params->windowShift, series->data->length - 2*params->windowShift);
    CHECKSTATUSPTR (status);
    /****************************************************************
     * 
     * clean up and return
     *
     ****************************************************************/
    DETATCHSTATUSPTR (status);
    RETURN (status);

}
