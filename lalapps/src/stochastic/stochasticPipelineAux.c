/* stochasticPipelineAux.c - SGWB Standalone Analysis Pipeline
 *                         - auxiliary functions
 *
 * Adam Mercer <ram@star.sr.bham.ac.uk>
 * Tania Regimbau <Tania.Regimbau@astro.cf.ac.uk>
 *
 * $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>

#include <unistd.h>
#include <getopt.h>

#include <FrameL.h>

#include <lal/AVFactories.h>
#include <lal/Calibration.h>
#include <lal/ComplexFFT.h>
#include <lal/RealFFT.h>
#include <lal/Date.h>
#include <lal/DetectorSite.h>
#include <lal/FrameCache.h>
#include <lal/FrameCalibration.h>
#include <lal/FrameStream.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALRCSID.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/PrintFTSeries.h>
#include <lal/PrintVector.h>
#include <lal/Random.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/StreamInput.h>
#include <lal/StochasticCrossCorrelation.h>
#include <lal/CoarseGrainFrequencySeries.h>
#include <lal/Units.h>
#include <lal/Window.h>
#include <lal/TimeFreqFFT.h>
#include <lal/SimulateSB.h>

#include "stochasticPipeline.h"

NRCSID ( PIPELINEAUXC, "$Id$" );

/* function to read data in frames */
void readDataPair(LALStatus *status,
		StreamPair *streamPair,
		ReadDataPairParams *params)
{
	/* counters */
	INT4 i;

	/* variables */
	FrCache *frCacheOne = NULL;
	FrStream *frStreamOne = NULL;
	FrCache *frCacheTwo = NULL;
	FrStream *frStreamTwo = NULL;
	FrChanIn frChanInOne;
	FrChanIn frChanInTwo;
	REAL4TimeSeries dataStreamOne;
	REAL4TimeSeries dataStreamTwo;
	ResampleTSParams resampleParams;
	LIGOTimeGPS bufferStartTime;
	UINT8 startTime;
	INT4 buffer;
	INT4 resampleRate;
	INT4 sampleRate;

	/* read parameters */
	startTime = params->startTime;
	buffer = params->buffer;
	resampleRate = params->resampleRate;
	sampleRate = params->sampleRate;

	/* initialise status pointer */
	INITSTATUS( status, "readDataPair", PIPELINEAUXC );
	ATTATCHSTATUSPTR( status );

	/* buffer start time */
	bufferStartTime.gpsSeconds = startTime - buffer;
	bufferStartTime.gpsNanoSeconds = 0;

	/* set channels */
	frChanInOne.name = params->channelOne;
	frChanInTwo.name = params->channelTwo;
	frChanInTwo.type = ADCDataChannel;
	frChanInOne.type = ADCDataChannel;

	/* initial data structures */
	dataStreamOne.epoch = bufferStartTime;
	dataStreamTwo.epoch = bufferStartTime;

	/* allocate memory */
	dataStreamOne.data = NULL;
	dataStreamTwo.data = NULL;
	LALSCreateVector(status->statusPtr, &(dataStreamOne.data), \
			sampleRate * (params->duration + (2 * buffer)));
	CHECKSTATUSPTR( status );
	LALSCreateVector(status->statusPtr, &(dataStreamTwo.data), \
			sampleRate * (params->duration + (2 * buffer)));
	CHECKSTATUSPTR( status );
	memset(dataStreamOne.data->data, 0, \
			dataStreamOne.data->length * sizeof(*dataStreamOne.data->data));
	memset(dataStreamTwo.data->data, 0, \
			dataStreamTwo.data->length * sizeof(*dataStreamTwo.data->data));

	/* open first frame cache */
	LALFrCacheImport(status->statusPtr, &frCacheOne, params->dataCacheOne);
	CHECKSTATUSPTR( status );
	LALFrCacheOpen(status->statusPtr, &frStreamOne, frCacheOne);
	CHECKSTATUSPTR( status );

	/* read first channel */
	LALFrSeek(status->statusPtr, &(bufferStartTime), frStreamOne);
	CHECKSTATUSPTR( status );
	LALFrGetREAL4TimeSeries(status->statusPtr, &dataStreamOne, \
			&frChanInOne, frStreamOne);
	CHECKSTATUSPTR( status );

	if (strcmp(params->dataCacheOne, params->dataCacheTwo) == 0)
	{
		/* read in second channel */
		LALFrSeek(status->statusPtr, &(bufferStartTime), frStreamOne);
		CHECKSTATUSPTR( status );
		LALFrGetREAL4TimeSeries(status->statusPtr, &dataStreamTwo, \
				&frChanInTwo, frStreamOne);
		CHECKSTATUSPTR( status );

		/* close frame cache */
		LALFrClose(status->statusPtr, &frStreamOne);
		CHECKSTATUSPTR( status );
	}
	else
	{
		/* close first frame cache */
		LALFrClose(status->statusPtr, &frStreamOne);
		CHECKSTATUSPTR( status );

		/* open second frame cache and read in second channel */
		LALFrCacheImport(status->statusPtr, &frCacheTwo, params->dataCacheTwo);
		CHECKSTATUSPTR( status );
		LALFrCacheOpen(status->statusPtr, &frStreamTwo, frCacheTwo);
		CHECKSTATUSPTR( status );
		LALFrSeek(status->statusPtr, &(bufferStartTime), frStreamTwo);
		CHECKSTATUSPTR( status );
		LALFrGetREAL4TimeSeries(status->statusPtr, &dataStreamTwo, \
				&frChanInTwo, frStreamTwo);
		CHECKSTATUSPTR( status );

		/* close second frame stream */
		LALFrClose(status->statusPtr, &frStreamTwo);
		CHECKSTATUSPTR( status );
	}

	/* resample */
	if (resampleRate != sampleRate)
	{
		/* set resample parameters */
		resampleParams.deltaT = 1.0 / (REAL8)resampleRate;
		resampleParams.filterType = defaultButterworth;

		/* resample */
		LALResampleREAL4TimeSeries(status->statusPtr, &dataStreamOne, \
				&resampleParams);
		CHECKSTATUSPTR( status );
		LALResampleREAL4TimeSeries(status->statusPtr, &dataStreamTwo, \
				&resampleParams);
		CHECKSTATUSPTR( status );
	}

	/* build output */
	strncpy(streamPair->streamOne->name,dataStreamOne.name, LALNameLength);
	strncpy(streamPair->streamTwo->name,dataStreamTwo.name, LALNameLength);
	streamPair->streamOne->epoch.gpsSeconds = startTime;
	streamPair->streamTwo->epoch.gpsSeconds = startTime;
	streamPair->streamOne->epoch.gpsNanoSeconds = 0;
	streamPair->streamTwo->epoch.gpsNanoSeconds = 0;
	streamPair->streamOne->deltaT = 1./(REAL8)resampleRate;
	streamPair->streamTwo->deltaT = 1./(REAL8)resampleRate;
	streamPair->streamOne->f0 = 0;
	streamPair->streamTwo->f0 = 0;
	streamPair->streamOne->sampleUnits = dataStreamOne.sampleUnits;
	streamPair->streamTwo->sampleUnits = dataStreamTwo.sampleUnits;

	/* remove buffer, and hence corruption due to resampling */
	for (i = 0; i < params->duration * resampleRate; i++)
	{
		streamPair->streamOne->data->data[i] = \
			dataStreamOne.data->data[i + (resampleRate * buffer)];
		streamPair->streamTwo->data->data[i] = \
			dataStreamTwo.data->data[i + (resampleRate * buffer)];
	}

	/* clean up */
	LALSDestroyVector(status->statusPtr, &(dataStreamOne.data));
	CHECKSTATUSPTR( status );
	LALSDestroyVector(status->statusPtr, &(dataStreamTwo.data));
	CHECKSTATUSPTR( status );

	/* return status */
	DETATCHSTATUSPTR( status );
	RETURN( status );
}

void monteCarlo(LALStatus *status,
		SSSimStochBGOutput *MCoutput,
		MonteCarloInput *MCinput,
		MonteCarloParams *MCparams)
{
	/* counters */
	INT4 i, loop;

	/* seed for random number generation */
	INT4 seed;

	/* various times, frequencies, sample rates */
	UINT8 startTime;
	UINT8 calTime;
	UINT4 length;
	INT4 lengthSegment;
	INT4 numSegment;
	UINT4 calibOffset;
	UINT4 freqLength;
	REAL8 deltaF;
	REAL8 deltaT;
	REAL8 sampleRate;

	/* IFO parameters */  
	INT4 siteOne;
	INT4 siteTwo;
	CHAR *ifoOne;
	CHAR *ifoTwo;

	/* parameters for calibration */
	CHAR *calCacheOne;
	CHAR *calCacheTwo;
	LIGOTimeGPS duration;
	const LALUnit countPerStrain = {0,{0,0,0,0,0,-1,1},{0,0,0,0,0,0,0}};

	/* omegaGW parameters */
	REAL8 omegaRef;
	REAL8 f0;
	REAL8 fRef;
	REAL8 alpha;
	StochasticOmegaGWParameters parametersOmega;

	/* other structures for SimulateSB */
	SSSimStochBGParams SBParams;
	SSSimStochBGInput SBInput;
	SSSimStochBGOutput SBOutput;
	REAL4TimeSeries whitenedSSimStochBGOne;
	REAL4TimeSeries whitenedSSimStochBGTwo;
	REAL4FrequencySeries omegaGW;
	COMPLEX8FrequencySeries responseOne;
	COMPLEX8FrequencySeries responseTwo;

	/* initialize status pointer */
	INITSTATUS(status, "monteCarlo", PIPELINEAUXC);
	ATTATCHSTATUSPTR(status);

	/* read parameters */
	fRef = MCparams->fRef;
	f0 = MCparams->f0;
	omegaRef = MCparams->omegaRef;
	alpha = MCparams->alpha;
	siteOne = MCparams->siteOne;
	siteTwo = MCparams->siteTwo;
	startTime = MCparams->startTime;
	lengthSegment = MCparams->lengthSegment;
	numSegment = MCparams->numSegment;
	sampleRate = MCparams->sampleRate;
	seed = MCparams->seed;

	/* derive other parameters */
	length = lengthSegment * numSegment;
	freqLength = lengthSegment / 2 + 1;
	calibOffset =  lengthSegment / (2 * sampleRate);
	deltaT = 1.0 / sampleRate;
	deltaF = 1.0 / (deltaT * lengthSegment);
	calTime = startTime;

	/* read input structure for calibration info */
	ifoOne = MCinput->ifoOne;
	ifoTwo = MCinput->ifoTwo;
	calCacheOne = MCinput->calCacheOne;
	calCacheTwo = MCinput->calCacheTwo;

	/* create vectors to store the simulated signal */

	whitenedSSimStochBGOne.data = NULL;
	whitenedSSimStochBGTwo.data = NULL;
	LALSCreateVector(status->statusPtr, &(whitenedSSimStochBGOne.data), \
			lengthSegment);
	CHECKSTATUSPTR( status );
	LALSCreateVector(status->statusPtr, &(whitenedSSimStochBGTwo.data), \
			lengthSegment);
	CHECKSTATUSPTR( status );
	memset(whitenedSSimStochBGOne.data->data, 0, \
			whitenedSSimStochBGOne.data->length * \
			sizeof(*whitenedSSimStochBGOne.data->data));
	memset(whitenedSSimStochBGTwo.data->data, 0, \
			whitenedSSimStochBGTwo.data->length * \
			sizeof(*whitenedSSimStochBGTwo.data->data));

	/* define parameters for SimulateSB */
	SBParams.length = lengthSegment;
	SBParams.deltaT = 1. / sampleRate;
	SBParams.detectorOne = lalCachedDetectors[siteOne];
	SBParams.detectorTwo = lalCachedDetectors[siteTwo];
	SBParams.SSimStochBGTimeSeries1Unit = lalADCCountUnit;
	SBParams.SSimStochBGTimeSeries2Unit = lalADCCountUnit;

	/* omegaGW */
	parametersOmega.length = freqLength;
	parametersOmega.f0 = f0;
	parametersOmega.deltaF = deltaF;
	parametersOmega.alpha = alpha;
	parametersOmega.fRef = fRef;
	parametersOmega.omegaRef = omegaRef;

	/* allocate memory */
	omegaGW.data = NULL;
	LALSCreateVector(status->statusPtr, &(omegaGW.data), freqLength);
	CHECKSTATUSPTR( status );
	memset(omegaGW.data->data, 0, \
			omegaGW.data->length * sizeof(*omegaGW.data->data));

	/* generate omegaGW */
	LALStochasticOmegaGW(status->statusPtr, &omegaGW, &parametersOmega);
	CHECKSTATUSPTR( status );

	/* response functions */
	/* set metadata fields */
	strncpy(responseOne.name,"responseOne", LALNameLength);
	strncpy(responseTwo.name,"responseTwo", LALNameLength);
	responseOne.sampleUnits = countPerStrain;
	responseTwo.sampleUnits = countPerStrain;
	responseOne.epoch.gpsNanoSeconds = 0;
	responseTwo.epoch.gpsNanoSeconds = 0;
	responseOne.deltaF = deltaF;
	responseTwo.deltaF = deltaF;
	responseOne.f0 = 0;
	responseTwo.f0 = 0;

	/* allocate memory */
	responseOne.data = NULL;
	responseTwo.data = NULL;
	LALCCreateVector(status->statusPtr, &(responseOne.data), freqLength);
	CHECKSTATUSPTR( status );
	LALCCreateVector(status->statusPtr, &(responseTwo.data), freqLength);
	CHECKSTATUSPTR( status );
	memset(responseOne.data->data, 0, \
	       responseOne.data->length * sizeof(*responseOne.data->data));
	memset(responseTwo.data->data, 0, \
	       responseTwo.data->length * sizeof(*responseTwo.data->data));

	/* set response function averaging */
	duration.gpsSeconds = 0;
	duration.gpsNanoSeconds = 0;

	/* construct a time series of the desired length */
	for (loop = 0; loop < numSegment; loop++)
	{
		/* increment seed */
		seed = seed + (2 * loop);
		SBParams.seed = seed;

		/* set calibration epoch */
		responseOne.epoch.gpsSeconds = calTime;
		responseTwo.epoch.gpsSeconds = calTime;

		/* compute response function */
		LALExtractFrameResponse(status->statusPtr, &responseOne, calCacheOne, \
				ifoOne, &duration);
		CHECKSTATUSPTR( status );
		LALExtractFrameResponse(status->statusPtr, &responseTwo, calCacheTwo, \
				ifoTwo, &duration);
		CHECKSTATUSPTR( status );

		/* force DC to be 0 and nyquist to be real */
		responseOne.data->data[0].re = 0.;
		responseTwo.data->data[0].re = 0.;
		responseOne.data->data[0].im = 0.;
		responseTwo.data->data[0].im = 0.;
		responseOne.data->data[freqLength-1].im = 0.;
		responseTwo.data->data[freqLength-1].im = 0.;

		/* define input structure for SimulateSB */
		SBInput.omegaGW = &omegaGW;
		SBInput.whiteningFilter1 = &responseOne;
		SBInput.whiteningFilter2 = &responseTwo;

		/* define output structure for SimulateSB */
		SBOutput.SSimStochBG1 = &whitenedSSimStochBGOne;
		SBOutput.SSimStochBG2 = &whitenedSSimStochBGTwo;

		/* generate whitened simulated SB data */
		LALSSSimStochBGTimeSeries(status->statusPtr, &SBOutput, &SBInput, \
				&SBParams);
		CHECKSTATUSPTR( status );

		/* get output */
		for (i = 0; i < lengthSegment; i++)
		{
			MCoutput->SSimStochBG1->data->data[i + (loop * lengthSegment)] = \
				whitenedSSimStochBGOne.data->data[i];
			MCoutput->SSimStochBG2->data->data[i + (loop * lengthSegment)] = \
				whitenedSSimStochBGTwo.data->data[i];
		}

		/* increase calibration time */
		calTime = calTime + calibOffset;
	}

	/* assign parameters and data to output */
	strncpy(MCoutput->SSimStochBG1->name, "Whitened-SimulatedSBOne", \
			LALNameLength);
	strncpy(MCoutput->SSimStochBG2->name, "Whitened-SimulatedSBTwo", \
			LALNameLength );
	MCoutput->SSimStochBG1->f0 = f0;
	MCoutput->SSimStochBG2->f0 = f0;
	MCoutput->SSimStochBG1->deltaT = deltaT;
	MCoutput->SSimStochBG2->deltaT = deltaT;
	MCoutput->SSimStochBG1->epoch.gpsSeconds = startTime;
	MCoutput->SSimStochBG2->epoch.gpsSeconds = startTime;
	MCoutput->SSimStochBG1->epoch.gpsNanoSeconds = 0;
	MCoutput->SSimStochBG2->epoch.gpsNanoSeconds = 0;
	MCoutput->SSimStochBG1->sampleUnits = lalADCCountUnit;
	MCoutput->SSimStochBG2->sampleUnits = lalADCCountUnit;

	/* clean up, and exit */
	LALSDestroyVector(status->statusPtr, &(omegaGW.data));
	CHECKSTATUSPTR( status );
	LALCDestroyVector(status->statusPtr, &(responseOne.data));
	CHECKSTATUSPTR( status );
	LALCDestroyVector(status->statusPtr, &(responseTwo.data));
	CHECKSTATUSPTR( status );
	LALSDestroyVector(status->statusPtr, &(whitenedSSimStochBGOne.data));
	CHECKSTATUSPTR( status );
	LALSDestroyVector(status->statusPtr, &(whitenedSSimStochBGTwo.data));
	CHECKSTATUSPTR( status );

	/* return status */
	DETATCHSTATUSPTR(status);
	RETURN(status);
}

void monteCarloSplice(LALStatus *status,
		SSSimStochBGOutput *MCoutput,
		MonteCarloInput *MCinput,
		MonteCarloParams *MCparams)
{
	/* counters */
	INT4 i, m, k, loop;

	/* seed for random number generation */
	INT4 seed;

	/* various times, frequencies, sample rates */
	UINT8 startTime;
	UINT8 calTime;
	UINT4 length;
	INT4 lengthSegment;
	INT4 numSegment;
	INT4 numSegmentSplice;
	INT4 numSegmentTot;
	UINT4 calibOffset;
	UINT4 spliceOffset;
	UINT4 freqLength;
	REAL8 deltaF;
	REAL8 deltaT;
	REAL8 sampleRate;

	/* IFO parameters */  
	INT4 siteOne;
	INT4 siteTwo;
	CHAR *ifoOne;
	CHAR *ifoTwo;

	/* parameters for calibration */
	CHAR *calCacheOne;
	CHAR *calCacheTwo;
	LIGOTimeGPS duration;
	const LALUnit countPerStrain = {0,{0,0,0,0,0,-1,1},{0,0,0,0,0,0,0}};

	/* omegaGW parameters */
	REAL8 omegaRef;
	REAL8 f0;
	REAL8 fRef;
	REAL8 alpha;
	StochasticOmegaGWParameters parametersOmega;

	/* other structures for SimulateSB */
	SSSimStochBGParams SBParams;
	SSSimStochBGInput SBInput;
	SSSimStochBGOutput SBOutput;
	REAL4TimeSeries whitenedSSimStochBGOne;
	REAL4TimeSeries whitenedSSimStochBGTwo;
	REAL4FrequencySeries omegaGW;
	COMPLEX8FrequencySeries responseOne;
	COMPLEX8FrequencySeries responseTwo;

	/* structures for splice function */
	REAL4Vector **longTrain1 = NULL;
	REAL4Vector **shortTrain1 = NULL;
	REAL4Vector **longTrain2 = NULL;
	REAL4Vector **shortTrain2 = NULL;

	/* initialize status pointer */
	INITSTATUS( status, "monteCarloSplice", PIPELINEAUXC );
	ATTATCHSTATUSPTR( status );

	/* read parameters */
	fRef = MCparams->fRef;
	f0 = MCparams->f0;
	omegaRef = MCparams->omegaRef;
	alpha = MCparams->alpha;
	siteOne = MCparams->siteOne;
	siteTwo = MCparams->siteTwo;
	startTime = MCparams->startTime;
	lengthSegment = MCparams->lengthSegment;
	numSegment = MCparams->numSegment;
	sampleRate = MCparams->sampleRate;
	seed = MCparams->seed;

	/* derive other parameters */
	length = lengthSegment * numSegment;
	numSegmentSplice = numSegment - 1;
	numSegmentTot = numSegment + numSegmentSplice;
	freqLength = (lengthSegment / 2) + 1;
	calibOffset =  lengthSegment / (2 * sampleRate);
	spliceOffset =  calibOffset;
	deltaT = 1.0 / sampleRate;
	deltaF = 1.0 / (deltaT * lengthSegment);
	calTime = startTime;

	/* read input structure for calibration info */
	ifoOne = MCinput->ifoOne;
	ifoTwo = MCinput->ifoTwo;
	calCacheOne = MCinput->calCacheOne;
	calCacheTwo = MCinput->calCacheTwo;

	/* create vectors to store the simulated signal */
	whitenedSSimStochBGOne.data = NULL;
	whitenedSSimStochBGTwo.data = NULL;
	LALSCreateVector(status->statusPtr, &(whitenedSSimStochBGOne.data), \
			lengthSegment);
	CHECKSTATUSPTR( status );
	LALSCreateVector(status->statusPtr, &(whitenedSSimStochBGTwo.data), \
			lengthSegment);
	CHECKSTATUSPTR( status );
	memset(whitenedSSimStochBGOne.data->data, 0, \
			whitenedSSimStochBGOne.data->length * \
			sizeof(*whitenedSSimStochBGOne.data->data));
	memset(whitenedSSimStochBGTwo.data->data, 0, \
			whitenedSSimStochBGTwo.data->length * \
			sizeof(*whitenedSSimStochBGTwo.data->data));

	/* allocate memory for longer data train */
	longTrain1 = (REAL4Vector**)LALMalloc(sizeof(REAL4Vector*) * numSegment);
	for (i = 0; i < numSegment; i++)
	{
		longTrain1[i] = NULL;
		LALSCreateVector(status->statusPtr, &longTrain1[i], lengthSegment);
		CHECKSTATUSPTR( status );
		for (k = 0; k < lengthSegment; k++)
		{
			longTrain1[i]->data[k] = 0;
		}
	}

	longTrain2 = (REAL4Vector**)LALMalloc(sizeof(REAL4Vector*) * numSegment);
	for (i = 0; i < numSegment; i++)
	{
		longTrain2[i] = NULL;
		LALSCreateVector(status->statusPtr, &longTrain2[i], lengthSegment);
		CHECKSTATUSPTR( status );
		for (k = 0; k < lengthSegment; k++)
		{
			longTrain2[i]->data[k] = 0;
		}
		longTrain2[i] = NULL;
	}

	/* allocate memory for shorter data train */
	shortTrain1 = (REAL4Vector**)LALMalloc(sizeof(REAL4Vector*) * \
			(numSegmentSplice));
	for(i = 0; i < numSegmentSplice; i++)
	{
		shortTrain1[i] = NULL;
		LALSCreateVector(status->statusPtr, &shortTrain1[i], lengthSegment);
		CHECKSTATUSPTR( status );
		for (k = 0; k < lengthSegment; k++)
		{
			shortTrain1[i]->data[k] = 0;
		}
	}
	shortTrain2 = (REAL4Vector**)LALMalloc(sizeof(REAL4Vector*) * \
			(numSegmentSplice));
	for(i = 0; i < numSegmentSplice; i++)
	{
		shortTrain2[i] = NULL;
		LALSCreateVector(status->statusPtr, &shortTrain2[i], lengthSegment);
		CHECKSTATUSPTR( status );
		for (k = 0; k < lengthSegment; k++)
		{
			shortTrain2[i]->data[k] = 0;
		}
	}

	/* define parameters for SimulateSB */
	SBParams.length = lengthSegment;
	SBParams.deltaT = 1. / sampleRate;
	SBParams.detectorOne = lalCachedDetectors[siteOne];
	SBParams.detectorTwo = lalCachedDetectors[siteTwo];
	SBParams.SSimStochBGTimeSeries1Unit = lalADCCountUnit;
	SBParams.SSimStochBGTimeSeries2Unit = lalADCCountUnit;

	/* omegaGW */
	parametersOmega.length = freqLength;
	parametersOmega.f0 = f0;
	parametersOmega.deltaF = deltaF;
	parametersOmega.alpha = alpha;
	parametersOmega.fRef = fRef;
	parametersOmega.omegaRef = omegaRef;

	/* allocate memory */
	omegaGW.data = NULL;
	LALSCreateVector(status->statusPtr, &(omegaGW.data), freqLength);
	CHECKSTATUSPTR( status );
	memset(omegaGW.data->data, 0, \
			omegaGW.data->length * sizeof(*omegaGW.data->data));

	/* generate omegaGW */
	LALStochasticOmegaGW(status->statusPtr, &omegaGW, &parametersOmega);
	CHECKSTATUSPTR( status );

	/* response functions */
	/* set metadata fields */
	strncpy(responseOne.name, "responseOne", LALNameLength);
	strncpy(responseTwo.name, "responseTwo", LALNameLength);
	responseOne.sampleUnits = countPerStrain;
	responseTwo.sampleUnits = countPerStrain;
	responseOne.epoch.gpsNanoSeconds = 0;
	responseTwo.epoch.gpsNanoSeconds = 0;
	responseOne.deltaF = deltaF;
	responseTwo.deltaF = deltaF;
	responseOne.f0 = 0;
	responseTwo.f0 = 0;

	/* allocate memory */
	responseOne.data = NULL;
	responseTwo.data = NULL;
	LALCCreateVector(status->statusPtr, &(responseOne.data), freqLength);
	CHECKSTATUSPTR( status );
	LALCCreateVector(status->statusPtr, &(responseTwo.data), freqLength);
	CHECKSTATUSPTR( status );
	memset(responseOne.data->data, 0, \
			responseOne.data->length * sizeof(*responseOne.data->data));
	memset(responseTwo.data->data, 0, \
			responseTwo.data->length * sizeof(*responseTwo.data->data));

	/* set response function averaging */
	duration.gpsSeconds = 0;
	duration.gpsNanoSeconds = 0;

	/* construct a time series of the desired length */
	for (loop = 0; loop < numSegmentTot; loop++)
	{
		/* increment seed */
		seed = seed + (2 * loop);
		SBParams.seed = seed;

		/* set calibration epoch */
		responseOne.epoch.gpsSeconds = calTime;
		responseTwo.epoch.gpsSeconds = calTime;

		/* compute response function */
		LALExtractFrameResponse(status->statusPtr, &responseOne, calCacheOne, \
				ifoOne, &duration);
		CHECKSTATUSPTR( status );
		LALExtractFrameResponse(status->statusPtr, &responseTwo, calCacheTwo, \
				ifoTwo, &duration);
		CHECKSTATUSPTR( status );

		/* force DC to be 0 and nyquist to be real */
		responseOne.data->data[0].re = 0.;
		responseTwo.data->data[0].re = 0.;
		responseOne.data->data[0].im = 0.;
		responseTwo.data->data[0].im = 0.;
		responseOne.data->data[freqLength-1].im = 0.;
		responseTwo.data->data[freqLength-1].im = 0.;

		/* define input structure for SimulateSB */
		SBInput.omegaGW = &omegaGW;
		SBInput.whiteningFilter1 = &responseOne;
		SBInput.whiteningFilter2 = &responseTwo;

		/* define output structure for SimulateSB */
		SBOutput.SSimStochBG1 = &whitenedSSimStochBGOne;
		SBOutput.SSimStochBG2 = &whitenedSSimStochBGTwo;

		/* generate whitened simulated SB data */
		LALSSSimStochBGTimeSeries(status->statusPtr, &SBOutput, &SBInput, \
				&SBParams);
		CHECKSTATUSPTR( status );
		
		/* store in the appropriate vector */
		if (loop % 2 == 0)
		{
			m = (UINT4)(loop / 2);
			longTrain1[m] = whitenedSSimStochBGOne.data;
			longTrain2[m] = whitenedSSimStochBGTwo.data;
		}
		else
		{
			m = (UINT4)((loop-1) / 2);
			shortTrain1[m] = whitenedSSimStochBGOne.data;
			shortTrain2[m] = whitenedSSimStochBGTwo.data;
		}

		/* increase calibration time */
		calTime = calTime + calibOffset;
	}

	/* splice long and short data trains into the output vector */
	SinusoidalSplice(longTrain1, shortTrain1, MCoutput->SSimStochBG1->data, \
			numSegmentSplice, spliceOffset);
	SinusoidalSplice(longTrain2, shortTrain2, MCoutput->SSimStochBG2->data, \
			numSegmentSplice, spliceOffset);

	/* assign parameters and data to output */
	strncpy(MCoutput->SSimStochBG1->name, "Whitened-SimulatedSBOne", \
			LALNameLength);
	strncpy(MCoutput->SSimStochBG2->name, "Whitened-SimulatedSBTwo", \
			LALNameLength);
	MCoutput->SSimStochBG1->f0 = f0;
	MCoutput->SSimStochBG2->f0 = f0;
	MCoutput->SSimStochBG1->deltaT = deltaT;
	MCoutput->SSimStochBG2->deltaT = deltaT;
	MCoutput->SSimStochBG1->epoch.gpsSeconds = startTime;
	MCoutput->SSimStochBG2->epoch.gpsSeconds = startTime;
	MCoutput->SSimStochBG1->epoch.gpsNanoSeconds = 0;
	MCoutput->SSimStochBG2->epoch.gpsNanoSeconds = 0;
	MCoutput->SSimStochBG1->sampleUnits = lalADCCountUnit;
	MCoutput->SSimStochBG2->sampleUnits = lalADCCountUnit;

	/* clean up, and exit */
	LALSDestroyVector(status->statusPtr, &(omegaGW.data));
	CHECKSTATUSPTR( status );
	LALCDestroyVector(status->statusPtr, &(responseOne.data));
	CHECKSTATUSPTR( status );
	LALCDestroyVector(status->statusPtr, &(responseTwo.data));
	CHECKSTATUSPTR( status );
	LALSDestroyVector(status->statusPtr, &(whitenedSSimStochBGOne.data));
	CHECKSTATUSPTR( status );
	LALSDestroyVector(status->statusPtr, &(whitenedSSimStochBGTwo.data));
	CHECKSTATUSPTR( status );

	/* return status */
	DETATCHSTATUSPTR(status);
	RETURN(status);
}

void SinusoidalSplice(REAL4Vector **longData,
		REAL4Vector **shortData,
		REAL4Vector *output,
		UINT4 nSpliceSegs,
		INT4 offset)
{
	/* counters */
	INT4 i, j;
	
	UINT4 segLen;
	INT4 nOutputPts;
	UINT4 nSplicePts;
	INT4 lastSpliceIndexPlusOne;
	UINT4 leftOverlap;
	UINT4 rightOverlap;
	UINT4 iMod;
	UINT4 jMod;
	REAL4 leftScale;
	REAL4 rightScale;
	REAL4 phase;

	segLen = longData[0]->length;
	nOutputPts = output->length;
	nSplicePts = nSpliceSegs * segLen;
	lastSpliceIndexPlusOne = nSplicePts + offset;

	leftOverlap = offset % segLen;
	rightOverlap = segLen - leftOverlap;

	/* leftOverlap and rightOverlap guaranteed to be non-zero by
	 * CheckArguments */
	leftScale = LAL_PI / (2.0 * leftOverlap);
	rightScale = LAL_PI / (2.0 * rightOverlap);

	for (i = 0; i < offset; i++)
	{
		output->data[i] = longData[i / segLen]->data[i % segLen];
	}

	/* check to make sure this is correct later */
	for(i = offset, j = 0; i < lastSpliceIndexPlusOne; i++, j++)
	{
		iMod = i % segLen;
		jMod = j % segLen;

		/* splice from start to middle of shortData segment */
		if (iMod < leftOverlap)
		{
			phase = iMod * leftScale;

			output->data[i] = (longData[i/segLen]->data[iMod] * sin(phase)) + \
				(shortData[j / segLen]->data[jMod] * cos(phase));
		}
		/* splice from middle to end of shortData segment */
		else
		{
			phase = (iMod - leftOverlap) * rightScale;

			output->data[i] = (longData[i/segLen]->data[iMod] * cos(phase)) + \
				(shortData[j / segLen]->data[jMod] * sin(phase));
		}
	}

	for (i = lastSpliceIndexPlusOne; i < nOutputPts; i++)
		output->data[i] = longData[i / segLen]->data[i % segLen];
}
