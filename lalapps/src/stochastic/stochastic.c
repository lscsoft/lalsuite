/*
 * stochasticPipeline.c - SGWB Standalone Analysis Pipeline
 *                      - core pipeline
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

#include <lalapps.h>

#include <lal/AVFactories.h>
#include <lal/Calibration.h>
#include <lal/ComplexFFT.h>
#include <lal/CoarseGrainFrequencySeries.h>
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
#include <lal/RealFFT.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/SimulateSB.h>
#include <lal/StreamInput.h>
#include <lal/StochasticCrossCorrelation.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Units.h>
#include <lal/Window.h>

#include "stochasticPipeline.h"

NRCSID (PIPELINEC, "$Id$");
RCSID ("$Id$");

#define CVS_REVISION "$Revision$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "stochasticPipeline"

/* variables for getopt options parsing */
char *optarg;
int optind;

/* fixed parameters for the stochastic search */

/* sampling parameters */
INT4 sampleRate = 16384;
INT4 resampleRate = 1024;
REAL8 deltaF = 0.25;

/* duration parameters */
INT4 streamDuration = 60;
INT4 segmentDuration = 60;
INT4 calibDuration = 60;

/* frequency band */
INT4 fMax = 300;
INT4 fMin = 40;

/* omegaGW parameters */
REAL4 alpha = 0.0;
REAL4 fRef = 100.0;
REAL4 omegaRef = 1.;

/* monte carlo parameters */
BOOLEAN simulData = 0;
BOOLEAN splice = 0;
REAL4 scaleFactor = 156.25;

/* ifo parameters */
INT4 siteOne;
INT4 siteTwo;

/* misc parameters */
INT4 hannDuration = 1;
INT4 padData = 1;
BOOLEAN applyMask = 0;
BOOLEAN PRINT = 0;

/* default input data parameters */
LIGOTimeGPS gpsStartTime;
UINT8 startTime = 731350400;
UINT4 seed = 173;
CHAR dataCacheOne[LALNameLength] = "cachefiles/H-731350382.cache";
CHAR dataCacheTwo[LALNameLength] = "cachefiles/H-731350382.cache";
CHAR calCacheOne[LALNameLength] = "calibration/H1-CAL-V03-729273600-734367600.cache";
CHAR calCacheTwo[LALNameLength] = "calibration/H2-CAL-V03-729296220-731849040.cache";
CHAR channelOne[LALNameLength] = "H1:LSC-AS_Q";
CHAR channelTwo[LALNameLength] = "H2:LSC-AS_Q";
CHAR ifoOne[LALNameLength] = "H1";
CHAR ifoTwo[LALNameLength] = "H2";

INT4 main(INT4 argc, CHAR *argv[])
{
	/* variable declarations */

	/* status pointer */
	LALStatus status;

	/* output file */
	FILE *out;
	CHAR outputFilename[LALNameLength];

	/* counters */
	INT4 i, segLoop;

	/* results parameters */
	REAL8 y;
	REAL8 yWhitenedSum = 0.;
	REAL8 varTheo;
	REAL8 inVarTheo;
	REAL8 inVarTheoSum = 0.;
	REAL8 pointEstimate;
	REAL8 errorBar;

	/* input data streams */
	INT4 numSegments;
	INT4 streamLength;
	ReadDataPairParams streamParams;
	StreamPair streamPair;
	REAL4TimeSeries streamOne;
	REAL4TimeSeries streamTwo;

	/* simulated signal structures */
	SSSimStochBGOutput MCoutput;
	MonteCarloParams MCparams;
	MonteCarloInput MCinput;
	UINT4 MCLength;

	/* simulated output structures */
	REAL4TimeSeries SimStochBGOne;
	REAL4TimeSeries SimStochBGTwo;

	/* data structures for segments */
	INT4 segmentLength;
	REAL4TimeSeries segmentOne;
	REAL4TimeSeries segmentTwo;

	/* data structures for PSDs */
	INT4 overlapPSDLength;
	INT4 psdTempLength;
	INT4 windowPSDLength;
	INT4 filterLength;
	INT4 numPointInf;
	INT4 numPointSup;
	LALWindowParams winparPSD;
	AverageSpectrumParams specparPSD;
	REAL4FrequencySeries psdTempOne;
	REAL4FrequencySeries psdTempTwo;
	REAL4FrequencySeries psdOne;
	REAL4FrequencySeries psdTwo;
	LALUnit psdUnits = {0,{0,0,1,0,0,0,2},{0,0,0,0,0,0,0}};

	/* response functions */
	COMPLEX8FrequencySeries responseTempOne;
	COMPLEX8FrequencySeries responseTempTwo;
	COMPLEX8FrequencySeries responseOne;
	COMPLEX8FrequencySeries responseTwo;
	INT4 respLength;
	LIGOTimeGPS gpsCalTime;
	LIGOTimeGPS duration;
	LALUnit countPerAttoStrain = {18,{0,0,0,0,0,-1,1},{0,0,0,0,0,0,0}};

	/* calibrated and half calibrated inverse noise data structures */
	REAL4FrequencySeries calInvPSDOne;
	REAL4FrequencySeries calInvPSDTwo;
	COMPLEX8FrequencySeries halfCalPSDOne;
	COMPLEX8FrequencySeries halfCalPSDTwo;

	/* units for inverse noise */
	LALUnit calPSDUnit = {36,{0,0,-1,0,0,-2,0},{0,0,0,0,0,0,0}};
	LALUnit halfCalPSDUnit = {18,{0,0,0,0,0,-1,-1},{0,0,0,0,0,0,0}};

	/* window for segment data streams */
	REAL4TimeSeries dataWindow;

	/* hann window */
	INT4 hannLength;
	LALWindowParams hannParams;
	REAL4Vector *hannWindow;

	/* zeropad and fft structures */
	SZeroPadAndFFTParameters zeroPadParams;
	RealFFTPlan *fftDataPlan = NULL;
	COMPLEX8FrequencySeries hBarTildeOne;
	COMPLEX8FrequencySeries hBarTildeTwo;
	UINT4 zeroPadLength;
	UINT4 fftDataLength;

	/* overlap reduction function */
	LALDetectorPair detectors;
	REAL4FrequencySeries overlap;
	OverlapReductionFunctionParameters ORFparams;
	LALUnit overlapUnit = {0,{0,0,0,0,0,2,0},{0,0,0,0,0,0,0}};

	/* frequency mask structures */
	REAL4FrequencySeries mask;
	REAL4Vector *maskTemp;

	/* spectrum structures */
	StochasticOmegaGWParameters omegaGWParams;
	REAL4FrequencySeries omegaGW;

	/* structures for LALInverseNoise */
	StochasticInverseNoiseInput inverseNoiseInOne;
	StochasticInverseNoiseInput inverseNoiseInTwo;
	StochasticInverseNoiseOutput inverseNoiseOutOne;
	StochasticInverseNoiseOutput inverseNoiseOutTwo;

	/* structures for optimal filter normalisation */
	StochasticOptimalFilterNormalizationInput normInput;
	StochasticOptimalFilterNormalizationOutput normOutput;
	StochasticOptimalFilterNormalizationParameters normParams;
	REAL4WithUnits normLambda;
	REAL4WithUnits normSigma;
	REAL8 lambda;

	/* structures for optimal filter */
	COMPLEX8FrequencySeries optFilter;
	StochasticOptimalFilterInput optFilterIn;

	/* structures for CC spectrum and CC statistics */
	StochasticCrossCorrelationInput ccIn;
	BOOLEAN epochsMatch = 1;
	REAL4WithUnits ccStat;
	COMPLEX8FrequencySeries ccSpectrum;

	/* error handler */
	status.statusPtr = NULL;
	lal_errhandler = LAL_ERR_EXIT;
	set_debug_level( "1" );

	/* parse command line options */
	parseOptions(argc, argv);

	/* open output file */
	LALSnprintf(outputFilename, LALNameLength, "output-%s%s-%d.dat", 
            ifoOne, ifoTwo, (INT4)startTime);
	if ((out = fopen(outputFilename, "w")) == NULL)
	{
		fprintf(stderr, "Can't open file for output...\n");
		exit(1);
	}

	if (PRINT)
	{
		fprintf(stdout, "Calculating number of segments...\n");
	}
	
	/* get number of segments */
	numSegments = streamDuration / segmentDuration;

	/* get stream duration and length */
	streamLength = streamDuration * resampleRate;

	/* set metadata fields for data streams */
	strncpy(streamOne.name, "streamOne", LALNameLength);
	strncpy(streamTwo.name, "streamTwo", LALNameLength);
	streamOne.sampleUnits = lalADCCountUnit;
	streamTwo.sampleUnits = lalADCCountUnit;
	streamOne.epoch = gpsStartTime;
	streamTwo.epoch = gpsStartTime;
	streamOne.deltaT = 1./(REAL8)resampleRate;
	streamTwo.deltaT = 1./(REAL8)resampleRate;
	streamOne.f0 = 0;
	streamTwo.f0 = 0;

	if (PRINT)
	{
		fprintf(stdout, "Allocating memory for data streams...\n");
	}

	/* allocate memory for data streams */
	streamOne.data = NULL;
	streamTwo.data = NULL;
	LAL_CALL( LALSCreateVector(&status, &(streamOne.data), streamLength), \
			&status );
	LAL_CALL( LALSCreateVector(&status, &(streamTwo.data), streamLength), \
			&status );
	memset(streamOne.data->data, 0, \
			streamOne.data->length * sizeof(*streamOne.data->data));
	memset(streamTwo.data->data, 0, \
			streamTwo.data->length * sizeof(*streamTwo.data->data));

	/* set stream input parameters */
	streamParams.duration = streamDuration;
	streamParams.dataCacheOne = dataCacheOne;
	streamParams.dataCacheTwo = dataCacheTwo;
	streamParams.ifoOne = ifoOne;
	streamParams.ifoTwo = ifoTwo;
	streamParams.channelOne = channelOne;
	streamParams.channelTwo = channelTwo;
	streamParams.startTime = startTime;
	streamParams.buffer = padData;
	streamParams.sampleRate = sampleRate;
	streamParams.resampleRate = resampleRate;

	/* set stream data structures */
	streamPair.streamOne = &streamOne;
	streamPair.streamTwo = &streamTwo;

	/* set parameter structure for monte carlo simulations */
	MCparams.lengthSegment = calibDuration * resampleRate;
	if ((streamDuration % calibDuration) == 0)
	{
		MCparams.numSegment = streamDuration / calibDuration;
	}
	else
	{
		MCparams.numSegment = (streamDuration / calibDuration) + 1;
	}
	MCparams.sampleRate = resampleRate;
	MCparams.startTime = startTime;
	MCparams.seed = seed;
	MCparams.fRef = fRef;
	MCparams.f0 = 0.;
	MCparams.omegaRef = omegaRef;
	MCparams.alpha = alpha;
	MCparams.siteOne = siteOne;
	MCparams.siteTwo = siteTwo;

	/* set input structure for monte carlo simulations */
	MCinput.ifoOne = ifoOne;
	MCinput.ifoTwo = ifoTwo;
	MCinput.calCacheOne = calCacheOne;
	MCinput.calCacheTwo = calCacheTwo;
	MCLength = MCparams.numSegment * MCparams.lengthSegment;

	if (PRINT)
	{
		fprintf(stdout, "Allocating memory for MC...\n");
	}

	/* allocate memory for monte carlo */
	SimStochBGOne.data = NULL;
	LAL_CALL( LALSCreateVector(&status, &(SimStochBGOne.data), MCLength), \
			&status );
	SimStochBGTwo.data = NULL;
	LAL_CALL( LALSCreateVector(&status, &(SimStochBGTwo.data), MCLength), \
			&status );
	memset(SimStochBGOne.data->data, 0, \
			SimStochBGOne.data->length * sizeof(*SimStochBGOne.data->data));
	memset(SimStochBGTwo.data->data, 0, \
			SimStochBGTwo.data->length * sizeof(*SimStochBGTwo.data->data));

	/* output structure for LALStochasticMC */
	MCoutput.SSimStochBG1 = &SimStochBGOne;
	MCoutput.SSimStochBG2 = &SimStochBGTwo;

	/* set length for data segments */
	segmentLength = segmentDuration * resampleRate;

	/* set metadata fields for data segments */
	strncpy(segmentOne.name, "segmentOne", LALNameLength);
	strncpy(segmentTwo.name, "segmentTwo", LALNameLength);
	segmentOne.sampleUnits = streamOne.sampleUnits;
	segmentTwo.sampleUnits = streamTwo.sampleUnits;
	segmentOne.epoch = gpsStartTime;
	segmentTwo.epoch = gpsStartTime;
	segmentOne.deltaT = 1./(REAL8)resampleRate;
	segmentTwo.deltaT = 1./(REAL8)resampleRate;
	segmentOne.f0 = 0;
	segmentTwo.f0 = 0;

	if (PRINT)
	{
		fprintf(stdout, "Allocating memory for data segments...\n");
	}

	/* allocate memory for data segments */
	segmentOne.data = NULL;
	segmentTwo.data = NULL;
	LAL_CALL( LALSCreateVector(&status, &(segmentOne.data), segmentLength), \
			&status );
	LAL_CALL( LALSCreateVector(&status, &(segmentTwo.data), segmentLength), \
			&status );
	memset(segmentOne.data->data, 0, \
			segmentOne.data->length * sizeof(*segmentOne.data->data));
	memset(segmentTwo.data->data, 0, \
			segmentTwo.data->length * sizeof(*segmentTwo.data->data));

	/* set PSD window length */
	windowPSDLength = (UINT4)(resampleRate / deltaF);

	/* set parameters for PSD estimation */
	overlapPSDLength = windowPSDLength / 2;
	psdTempLength = (windowPSDLength / 2) + 1;
	numPointInf = (UINT4)(fMin / deltaF);
	numPointSup = (UINT4)(fMax / deltaF);
	filterLength = numPointSup - numPointInf + 1;

	/* set metadata fields for PSDs */
	strncpy(psdTempOne.name, "psdTempOne", LALNameLength);
	strncpy(psdTempTwo.name, "psdTempTwo", LALNameLength);
	psdTempOne.sampleUnits = psdUnits;
	psdTempTwo.sampleUnits = psdUnits;
	psdTempOne.deltaF = deltaF;
	psdTempTwo.deltaF = deltaF;
	psdTempOne.f0 = 0;
	psdTempTwo.f0 = 0;

	if (PRINT)
	{
		fprintf(stdout, "Allocating memory for PSDs...\n");
	}

	/* allocate memory for PSDs */
	psdTempOne.data = NULL;
	psdTempTwo.data = NULL;
	LAL_CALL( LALCreateVector(&status, &(psdTempOne.data), psdTempLength), \
			&status );
	LAL_CALL( LALCreateVector(&status, &(psdTempTwo.data), psdTempLength), \
			&status );
	memset(psdTempOne.data->data, 0, \
			psdTempOne.data->length * sizeof(*psdTempOne.data->data));
	memset(psdTempTwo.data->data, 0, \
			psdTempTwo.data->length * sizeof(*psdTempTwo.data->data));

	/* reduced frequency band PSDs */
	/* set metadata fields for reduced frequency band PSDs */
	strncpy(psdOne.name, "psdOne", LALNameLength);
	strncpy(psdTwo.name, "psdTwo", LALNameLength);
	psdOne.deltaF = deltaF;
	psdTwo.deltaF = deltaF;
	psdOne.f0 = fMin;
	psdTwo.f0 = fMin;
	psdOne.sampleUnits = psdUnits;
	psdTwo.sampleUnits = psdUnits;

	if (PRINT)
	{
		fprintf(stdout, "Allocating memory for reduced frequency band PSDs...\n");
	}

	/* allocate memory for reduced frequency band PSDs */
	psdOne.data = NULL;
	psdTwo.data = NULL;
	LAL_CALL( LALCreateVector(&status, &psdOne.data, filterLength), &status );
	LAL_CALL( LALCreateVector(&status, &psdTwo.data, filterLength), &status );
	memset (psdOne.data->data, 0, \
			psdOne.data->length * sizeof(*psdOne.data->data));
	memset( psdTwo.data->data, 0, \
			psdTwo.data->length * sizeof(*psdTwo.data->data));

	/* set window parameters for PSD estimation */
	winparPSD.length = windowPSDLength;
	winparPSD.type = Hann;

	/* set parameters for PSD estimation */
	specparPSD.method = useMean;
	specparPSD.overlap = overlapPSDLength;
	specparPSD.plan = NULL;
	specparPSD.window = NULL;

	if (PRINT)
	{
		fprintf(stdout, "Creating FFT plan for PSD estimation...\n");
	}

	/* create fft plan */
	LAL_CALL (LALCreateForwardRealFFTPlan(&status, &specparPSD.plan, \
				windowPSDLength, 0), &status );

	if (PRINT)
	{
		fprintf(stdout, "Creating window for PSD estimation...\n");
	}

	/* create window for PSD estimation */
	LAL_CALL( LALCreateREAL4Window(&status, &specparPSD.window, &winparPSD), \
			&status );

	/* set parameters for response functions */
	respLength = (UINT4)(fMax / deltaF) + 1;
	duration.gpsSeconds = 0;
	duration.gpsNanoSeconds = 0;
	gpsCalTime.gpsSeconds = startTime;
	gpsCalTime.gpsNanoSeconds = 0;

	/* set metadata fields for response functions */
	strncpy(responseTempOne.name, "responseTempOne", LALNameLength);
	strncpy(responseTempTwo.name, "responseTempTwo", LALNameLength);
	responseTempOne.sampleUnits = countPerAttoStrain;
	responseTempTwo.sampleUnits = countPerAttoStrain;
	responseTempOne.epoch = gpsCalTime;
	responseTempTwo.epoch = gpsCalTime;
	responseTempOne.deltaF = deltaF;
	responseTempTwo.deltaF = deltaF;
	responseTempOne.f0 = 0;
	responseTempTwo.f0 = 0;

	if (PRINT)
	{
		fprintf(stdout, "Allocating memory for response functions...\n");
	}

	/* allocate memory for response functions */
	responseTempOne.data = NULL;
	responseTempTwo.data = NULL;
	LAL_CALL( LALCCreateVector(&status, &(responseTempOne.data), respLength), \
			&status );
	LAL_CALL( LALCCreateVector(&status, &(responseTempTwo.data), respLength), \
			&status );
	memset(responseTempOne.data->data, 0, \
			responseTempOne.data->length * sizeof(*responseTempOne.data->data));
	memset(responseTempTwo.data->data, 0, \
			responseTempTwo.data->length * sizeof(*responseTempTwo.data->data));

	/* reduced frequency band response functions */
	/* set metadata fields for reduced frequency band response functions */
	strncpy(responseOne.name, "responseOne", LALNameLength);
	strncpy(responseTwo.name, "responseTwo", LALNameLength);
	responseOne.sampleUnits = countPerAttoStrain;
	responseTwo.sampleUnits = countPerAttoStrain;
	responseOne.epoch = gpsCalTime;
	responseTwo.epoch = gpsCalTime;
	responseOne.deltaF = deltaF;
	responseTwo.deltaF = deltaF;
	responseOne.f0 = fMin;
	responseTwo.f0 = fMin;

	if (PRINT)
	{
		fprintf(stdout, "Allocating memory for reduced frequency band response " \
				"functions...\n");
	}

	/* allocate memory for reduced frequency band response functions */
	responseOne.data = NULL;
	responseTwo.data = NULL;
	LAL_CALL( LALCCreateVector(&status, &(responseOne.data), filterLength), \
			&status );
	LAL_CALL( LALCCreateVector(&status, &(responseTwo.data), filterLength), \
			&status );
	memset(responseOne.data->data, 0, \
			responseOne.data->length * sizeof(*responseOne.data->data));
	memset(responseTwo.data->data, 0, \
			responseTwo.data->length * sizeof(*responseTwo.data->data));

	/* set metadata fields for inverse noise structures */
	strncpy(calInvPSDOne.name, "calInvPSDOne", LALNameLength);
	strncpy(calInvPSDTwo.name, "calInvPSDTwo", LALNameLength);
	calInvPSDOne.sampleUnits = calPSDUnit;
	calInvPSDTwo.sampleUnits = calPSDUnit;
	calInvPSDOne.deltaF = deltaF;
	calInvPSDTwo.deltaF = deltaF;
	calInvPSDOne.f0 = fMin;
	calInvPSDTwo.f0 = fMin;
	strncpy(halfCalPSDOne.name, "halfCalPSDOne", LALNameLength);
	strncpy(halfCalPSDTwo.name, "halfCalPSDTwo", LALNameLength);
	halfCalPSDOne.deltaF = deltaF;
	halfCalPSDTwo.deltaF = deltaF;
	halfCalPSDOne.f0 = fMin;
	halfCalPSDTwo.f0 = fMin;
	halfCalPSDOne.sampleUnits = halfCalPSDUnit;
	halfCalPSDTwo.sampleUnits = halfCalPSDUnit;

	if (PRINT)
	{
		fprintf(stdout, "Allocating memory for inverse noise...\n");
	}

	/* allocate memory for inverse noise */
	calInvPSDOne.data = NULL;
	calInvPSDTwo.data = NULL;
	LAL_CALL( LALCreateVector(&status, &(calInvPSDOne.data), filterLength), \
			&status );
	LAL_CALL( LALCreateVector(&status, &(calInvPSDTwo.data), filterLength), \
			&status );
	memset(calInvPSDOne.data->data, 0, \
			calInvPSDOne.data->length * sizeof(*calInvPSDOne.data->data));
	memset(calInvPSDTwo.data->data, 0, \
			calInvPSDTwo.data->length * sizeof(*calInvPSDTwo.data->data));
	halfCalPSDOne.data = NULL;
	halfCalPSDTwo.data = NULL;
	LAL_CALL( LALCCreateVector(&status, &(halfCalPSDOne.data), filterLength), \
			&status );
	LAL_CALL( LALCCreateVector(&status, &(halfCalPSDTwo.data), filterLength), \
			&status );
	memset(halfCalPSDOne.data->data, 0, \
			halfCalPSDOne.data->length * sizeof(*halfCalPSDOne.data->data));
	memset(halfCalPSDTwo.data->data, 0, \
			halfCalPSDTwo.data->length * sizeof(*halfCalPSDTwo.data->data));

	/* set inverse noise inputs */
	inverseNoiseInOne.unCalibratedNoisePSD = &psdOne;
	inverseNoiseInOne.responseFunction = &responseOne;
	inverseNoiseInTwo.unCalibratedNoisePSD = &psdTwo;
	inverseNoiseInTwo.responseFunction = &responseTwo;

	/* set inverse noise outputs */
	inverseNoiseOutOne.calibratedInverseNoisePSD = &calInvPSDOne;
	inverseNoiseOutOne.halfCalibratedInverseNoisePSD = &halfCalPSDOne;
	inverseNoiseOutTwo.calibratedInverseNoisePSD = &calInvPSDTwo;
	inverseNoiseOutTwo.halfCalibratedInverseNoisePSD = &halfCalPSDTwo;

	/* set window parameters for segment data streams */
	strncpy(dataWindow.name, "dataWindow", LALNameLength);
	dataWindow.sampleUnits = lalDimensionlessUnit;
	dataWindow.deltaT = 1./resampleRate;
	dataWindow.f0 = 0;

	if (PRINT)
	{
		fprintf(stdout, "Allocating memory for data segment window...\n");
	}

	/* allocate memory for segment window */
	dataWindow.data = NULL;
	LAL_CALL( LALSCreateVector(&status, &(dataWindow.data), segmentLength), \
			&status );
	memset(dataWindow.data->data, 0, \
			dataWindow.data->length * sizeof(*dataWindow.data->data));

	if (PRINT)
	{
		fprintf(stdout, "Generating data segment window...\n");
	}

	/* generate window */
	for (i = 0; i < segmentLength; i++)
	{
		dataWindow.data->data[i] = 1.;
	}

	if (hannDuration != 0)
	{
		/* generate pure Hann window */
		hannLength = hannDuration * resampleRate;
		hannParams.length = hannLength;
		hannParams.type = Hann;

		/* allocate memory for hann window */
		hannWindow = NULL;
		LAL_CALL( LALSCreateVector(&status, &hannWindow, hannLength), &status );
		memset(hannWindow->data, 0, \
				hannWindow->length * sizeof(*hannWindow->data));

		/* generate hann window */
		LAL_CALL( LALWindow(&status, hannWindow, &hannParams), &status );

		/* construct Tukey window */
		for (i = 0; i < hannLength / 2; i++)
		{
			dataWindow.data->data[i] = hannWindow->data[i];
		}
		for (i = segmentLength - (hannLength / 2); i < segmentLength; i++)
		{
			dataWindow.data->data[i] = \
				hannWindow->data[i - segmentLength + hannLength];
		}
	}

	/* save window */
	if (PRINT)
	{
		LALSPrintTimeSeries(&dataWindow, "dataWindow.dat");
	}

	/* zeropad lengths */
	zeroPadLength = 2 * segmentLength;
	fftDataLength = (zeroPadLength / 2) + 1;

	/* create fft plan */
	LAL_CALL( LALCreateForwardRealFFTPlan(&status, &fftDataPlan, \
				zeroPadLength, 0), &status );

	/* set metadata fields for zeropad ffts */
	strncpy(hBarTildeOne.name, "hBarTildeOne", LALNameLength);
	strncpy(hBarTildeTwo.name, "hBarTildeTwo", LALNameLength);

	if (PRINT)
	{
		fprintf(stdout, "Allocating memory for zeropad...\n");
	}

	/* allocate memory for zeropad */
	hBarTildeOne.data = NULL;
	hBarTildeTwo.data = NULL;
	LAL_CALL( LALCCreateVector(&status, &(hBarTildeOne.data), fftDataLength), \
			&status );
	LAL_CALL( LALCCreateVector(&status, &(hBarTildeTwo.data), fftDataLength), \
			&status );
	memset(hBarTildeOne.data->data, 0, \
			hBarTildeOne.data->length * sizeof(*hBarTildeOne.data->data));
	memset(hBarTildeTwo.data->data, 0, \
			hBarTildeTwo.data->length * sizeof(*hBarTildeTwo.data->data));

	/* set zeropad parameters */
	zeroPadParams.fftPlan = fftDataPlan;
	zeroPadParams.window = dataWindow.data;
	zeroPadParams.length = zeroPadLength;

	/* quantities needed to build the optimal filter */

	/* set metadata fields for overlap reduction function */
	strncpy(overlap.name, "overlap", LALNameLength);
	overlap.sampleUnits = overlapUnit;
	overlap.deltaF = deltaF;
	overlap.f0 = fMin;

	if (PRINT)
	{
		fprintf(stdout, "Allocating memory for the overlap reduction " \
				"function...\n");
	}

	/* allocate memory for overlap reduction function */
	overlap.data = NULL;
	LAL_CALL( LALCreateVector(&status, &(overlap.data), filterLength), &status );
	memset(overlap.data->data, 0, \
			overlap.data->length * sizeof(*overlap.data->data));

	/* set parameters for overlap reduction function */
	ORFparams.length = filterLength;
	ORFparams.f0 = fMin;
	ORFparams.deltaF = deltaF;

	/* set overlap reduction function detector pair */
	detectors.detectorOne = lalCachedDetectors[siteOne];
	detectors.detectorTwo = lalCachedDetectors[siteTwo];

	if (PRINT)
	{
		fprintf(stdout, "Generating the overlap reduction function...\n");
	}

	/* generate overlap reduction function */
	LAL_CALL (LALOverlapReductionFunction(&status, &overlap, &detectors, \
				&ORFparams), &status );

	/* save */
	if (PRINT)
	{
		LALSPrintFrequencySeries(&overlap, "overlap.dat");
	}

	/* set metadata fields for spectrum */
	strncpy(omegaGW.name, "omegaGW", LALNameLength);
	omegaGW.sampleUnits = lalDimensionlessUnit;
	omegaGW.deltaF = deltaF;
	omegaGW.f0 = fMin;

	if (PRINT)
	{
		fprintf(stdout, "Allocating memory for spectrum...\n");
	}

	/* allocate memory for spectrum */
	omegaGW.data = NULL;
	LAL_CALL( LALCreateVector(&status, &(omegaGW.data), filterLength), &status );
	memset(omegaGW.data->data, 0, \
			omegaGW.data->length * sizeof(*omegaGW.data->data));

	/* set omegaGW parameters */
	omegaGWParams.alpha = alpha;
	omegaGWParams.fRef = fRef;
	omegaGWParams.omegaRef = omegaRef;
	omegaGWParams.length = filterLength;
	omegaGWParams.f0 = fMin;
	omegaGWParams.deltaF = deltaF;

	if (PRINT)
	{
		fprintf(stdout, "Generating spectrum for optimal filter...\n");
	}

	/* generage omegaGW */
	LAL_CALL( LALStochasticOmegaGW(&status, &omegaGW, &omegaGWParams), &status );

	/* frequency mask */
	if (applyMask)
	{
		/* set metadata fields for frequency mask */
		strncpy(mask.name, "mask", LALNameLength);
		mask.deltaF = deltaF;
		mask.f0 = fMin;
		mask.sampleUnits = lalDimensionlessUnit;

		if (PRINT)
		{
			fprintf(stdout, "Allocating memory for frequency mask...\n");
		}

		/* allocate memory for frequency mask */
		mask.data = NULL;
		maskTemp = NULL;
		LAL_CALL( LALCreateVector(&status, &(mask.data), filterLength), &status );
		LAL_CALL( LALCreateVector(&status, &maskTemp, respLength), &status );
		memset(mask.data->data, 0, \
				mask.data->length * sizeof(*mask.data->data));
		memset(maskTemp->data, 0, \
				maskTemp->length * sizeof(*maskTemp->data));

		if (PRINT)
		{
			fprintf(stdout, "Generating frequency mask...\n");
		}

		/* set all values to 1 */
		for (i = 0; i < respLength; i++)
		{
			maskTemp->data[i] = 1.;
		}

		if (PRINT)
		{
			fprintf(stdout, "Masking multiples of 16 Hz...\n");
		}

		/* remove multiples of 16 Hz */
		for (i = 0; i < respLength; i += (UINT4)(16 / deltaF))
		{
			maskTemp->data[i] = 0.;
		}

		if (PRINT)
		{
			fprintf(stdout, "Masking multiples of 60 Hz...\n");
		}

		/* remove multiples of 60 Hz */
		for (i = 0; i < respLength; i += (UINT4)(60 / deltaF))
		{
			maskTemp->data[i] = 0.;
		}

		if (PRINT)
		{
			fprintf(stdout, "Getting appropriate frequency band for mask...\n");
		}

		/* get appropriate band */
		for (i = 0; i < filterLength; i++)
		{
			mask.data->data[i] = maskTemp->data[i + numPointInf];
		}

		if (PRINT)
		{
			LALSPrintFrequencySeries(&mask, "mask.dat");
			fprintf(stdout, "Applying frequency mask to spectrum..\n");
		}

		/* apply mask to omegaGW */
		for (i = 0; i < filterLength; i++)
		{
			omegaGW.data->data[i] *= mask.data->data[i];
		}
	}

	/* save */
	if (PRINT)
	{
		LALSPrintFrequencySeries(&omegaGW, "omegaGW.dat");
	}

	/* set normalisation parameters */
	normParams.fRef = fRef;
	normParams.heterodyned = 0;
	normParams.window1 = dataWindow.data;
	normParams.window2 = dataWindow.data;

	/* set normalisation input */
	normInput.overlapReductionFunction = &overlap;
	normInput.omegaGW = &omegaGW;
	normInput.inverseNoisePSD1 = &calInvPSDOne;
	normInput.inverseNoisePSD2 = &calInvPSDTwo;

	/* set normalisation output */
	normOutput.normalization = &normLambda;
	normOutput.variance = &normSigma;

	/* set metadata fields for optimal filter */
	strncpy(optFilter.name, "optFilter", LALNameLength);
	optFilter.epoch = gpsCalTime;
	optFilter.deltaF = deltaF;
	optFilter.f0 = fMin;

	if (PRINT)
	{
		fprintf(stdout, "Allocating memory for optimal filter...\n");
	}

	/* allocate memory for optimal filter */
	optFilter.data = NULL;
	LAL_CALL( LALCCreateVector(&status, &(optFilter.data), filterLength), \
			&status );
	memset(optFilter.data->data, 0, \
			optFilter.data->length * sizeof(*optFilter.data->data));

	/* set optimal filter inputs */
	optFilterIn.overlapReductionFunction = &overlap;
	optFilterIn.omegaGW = &omegaGW;
	optFilterIn.halfCalibratedInverseNoisePSD1 = &halfCalPSDOne;
	optFilterIn.halfCalibratedInverseNoisePSD2 = &halfCalPSDTwo;

	/* set metadata fields for CC spectrum */
	strncpy(ccSpectrum.name, "ccSpectrum", LALNameLength);
	ccSpectrum.epoch = gpsCalTime;
	ccSpectrum.deltaF = deltaF;
	ccSpectrum.f0 = fMin;

	/* allocate memory for CC spectrum*/
	ccSpectrum.data = NULL;
	LAL_CALL( LALCCreateVector(&status, &(ccSpectrum.data), filterLength), \
			&status );
	memset(ccSpectrum.data->data, 0, \
			ccSpectrum.data->length * sizeof(*ccSpectrum.data->data));

	/* set CC inputs */
	ccIn.hBarTildeOne = &hBarTildeOne;
	ccIn.hBarTildeTwo = &hBarTildeTwo;
	ccIn.optimalFilter = &optFilter;

	/*
	 ** read data, downsample and eventually inject simulated signal **
	 */

	if (PRINT)
	{
		fprintf(stdout, "Reading data...\n");
	}

	/* read data */
	readDataPair(&status, &streamPair, &streamParams);

	/* save */
	if (PRINT)
	{
		LALSPrintTimeSeries(&streamOne, "stream1.dat");
		LALSPrintTimeSeries(&streamTwo, "stream2.dat");
	}

	/* simulate signal */
	if (simulData)
	{
		/* perform monte carlo */
		if (splice)
		{
			monteCarloSplice(&status, &MCoutput, &MCinput, &MCparams);
		}
		else
		{
			monteCarlo(&status, &MCoutput, &MCinput, &MCparams);
		}

		/* output the results */
		if (PRINT)
		{
			LALSPrintTimeSeries(&SimStochBGOne, "simStochBG1.dat");
			LALSPrintTimeSeries(&SimStochBGTwo, "simStochBG2.dat");
		}

		/* multiply by scale factor and inject into real data */
		for (i = 0; i < streamLength ; i++)
		{
			streamOne.data->data[i] = streamOne.data->data[i] + \
				(scaleFactor * SimStochBGOne.data->data[i]);
			streamTwo.data->data[i] = streamTwo.data->data[i] + \
				(scaleFactor * SimStochBGTwo.data->data[i]);
		}

		/* clean up */
		LAL_CALL( LALDestroyVector(&status, &SimStochBGOne.data), &status );
		LAL_CALL( LALDestroyVector(&status, &SimStochBGTwo.data), &status );
	}

	/* output the results */
	if (PRINT)
	{
		LALSPrintTimeSeries(&streamOne, "signal1.dat");
		LALSPrintTimeSeries(&streamTwo, "signal2.dat");
	}

	/*
	 ** loop over segments **
	 */

	if (PRINT)
	{
		fprintf(stdout, "Looping over %d segments...\n", numSegments);
	}

	for (segLoop = 0; segLoop < numSegments; segLoop++)
	{
		/* define segment epoch */
		gpsStartTime.gpsSeconds = startTime + (segLoop * segmentDuration);
		gpsCalTime.gpsSeconds = gpsStartTime.gpsSeconds;
		segmentOne.epoch = gpsStartTime;
		segmentTwo.epoch = gpsStartTime;

		if (PRINT)
		{
			fprintf(stdout, "Performing search on segment %d of %d...\n", \
					segLoop + 1, numSegments);
		}

		/* build small segment */
		for (i = 0; i < segmentLength ; i++)
		{
			segmentOne.data->data[i] = streamOne.data->data[i + \
				(segLoop * segmentLength)];
			segmentTwo.data->data[i] = streamTwo.data->data[i + \
				(segLoop * segmentLength)];
		}

		/* output the results */
		if (PRINT)
		{
			LALSPrintTimeSeries(&segmentOne, "segment1.dat");
			LALSPrintTimeSeries(&segmentTwo, "segment2.dat");
		}

		if (PRINT)
		{
			fprintf(stdout, "Zero padding data...\n");
		}

		/* zero pad and fft */
		LAL_CALL( LALSZeroPadAndFFT(&status, &hBarTildeOne, &segmentOne, \
					&zeroPadParams), &status );
		LAL_CALL( LALSZeroPadAndFFT(&status, &hBarTildeTwo, &segmentTwo, \
					&zeroPadParams), &status );

		/* save */
		if (PRINT)
		{
			LALCPrintFrequencySeries(&hBarTildeOne, "hBarTilde1.dat");
			LALCPrintFrequencySeries(&hBarTildeTwo, "hBarTilde2.dat");
		}

		if (PRINT)
		{
			fprintf(stdout, "Estimating PSDs...\n");
		}

		/* compute uncalibrated PSDs */
		LAL_CALL( LALREAL4AverageSpectrum(&status, &psdTempOne, &segmentOne, \
					&specparPSD), &status );
		LAL_CALL( LALREAL4AverageSpectrum(&status, &psdTempTwo, &segmentTwo, \
					&specparPSD), &status );

		if (PRINT)
		{
			fprintf(stdout, "Getting appropriate frequency band for PSDs...\n");
		}

		/* reduce to the optimal filter frequency range */
		for (i = 0; i < filterLength; i++)
		{
			psdOne.data->data[i] = psdTempOne.data->data[i + numPointInf];
			psdTwo.data->data[i] = psdTempTwo.data->data[i + numPointInf];
		}

		/* output the results */
		if (PRINT)
		{
			LALSPrintFrequencySeries(&psdOne, "psd1.dat");
			LALSPrintFrequencySeries(&psdTwo, "psd2.dat");
		}

		if (PRINT)
		{
			fprintf(stdout, "Generating response functions...\n");
		}

		/* compute response function */
		responseTempOne.epoch = gpsCalTime;
		responseTempTwo.epoch = gpsCalTime;
		LAL_CALL( LALExtractFrameResponse(&status, &responseTempOne, calCacheOne, \
				ifoOne, &duration), &status );
		LAL_CALL( LALExtractFrameResponse(&status, &responseTempTwo, calCacheTwo, \
				ifoTwo, &duration), &status );

		if (PRINT)
		{
			fprintf(stdout, "Getting appropriate frequency band for response " \
					"functions...\n");
		}

		/* reduce to the optimal filter frequency range */
		responseOne.epoch = gpsCalTime;
		responseTwo.epoch = gpsCalTime;
		for (i = 0; i < filterLength; i++)
		{
			responseOne.data->data[i] = responseTempOne.data->data[i + numPointInf];
			responseTwo.data->data[i] = responseTempTwo.data->data[i + numPointInf];
		}

		/* output the results */
		if (PRINT)
		{
			LALCPrintFrequencySeries(&responseOne, "response1.dat");
			LALCPrintFrequencySeries(&responseTwo, "response2.dat");
		}

		if (PRINT)
		{
			fprintf(stdout, "Generating inverse noise...\n");
		}

		/* compute inverse calibrate and half calibrated PSDs */
		LAL_CALL( LALStochasticInverseNoise(&status, &inverseNoiseOutOne, \
				&inverseNoiseInOne), &status );
		LAL_CALL( LALStochasticInverseNoise(&status, &inverseNoiseOutTwo, \
				&inverseNoiseInTwo), &status );

		/* save */
		if (PRINT)
		{
			LALSPrintFrequencySeries(&calInvPSDOne, "inPSD1.dat");
			LALSPrintFrequencySeries(&calInvPSDTwo, "inPSD2.dat");
			LALCPrintFrequencySeries(&halfCalPSDOne, "hInPSD1.dat");
			LALCPrintFrequencySeries(&halfCalPSDTwo, "hInPSD2.dat");
		}

		if (PRINT)
		{
			fprintf(stdout, "Normalising optimal filter...\n");
		}

		/* compute variance and normalisation for optimal filter */
		LAL_CALL( LALStochasticOptimalFilterNormalization(&status, &normOutput, \
				&normInput, &normParams), &status );
		lambda = (REAL8)(normLambda.value * \
				pow(10.,normLambda.units.powerOfTen));
		varTheo = (REAL8)(segmentDuration * normSigma.value * \
				pow(10.,normSigma.units.powerOfTen));
		inVarTheo = 1./varTheo;
		inVarTheoSum = inVarTheoSum + inVarTheo;

		/* save */
		if (PRINT)
		{
			fprintf(stdout, "lambda = %e s^-1\n", lambda);
			fprintf(stdout, "varTheo = %e s\n", varTheo);
		}

		if (PRINT)
		{
			fprintf(stdout, "Generating optimal filter...\n");
		}

		/* build optimal filter */
		optFilter.epoch = gpsCalTime;
		LAL_CALL( LALStochasticOptimalFilter(&status, &optFilter, &optFilterIn, \
					&normLambda), &status );

		/* save */
		if (PRINT)
		{
			LALCPrintFrequencySeries(&optFilter, "optFilter.dat");
		}

		if (PRINT)
		{
			fprintf(stdout, "Generating cross correlation spectrum...\n");
		}

		/* cc spectrum */
		LAL_CALL( LALStochasticCrossCorrelationSpectrum(&status, &ccSpectrum, \
					&ccIn, epochsMatch), &status );

		/* save */
		if (PRINT)
		{
			LALCPrintFrequencySeries(&ccSpectrum, "ccSpectrum.dat");
		}

		if (PRINT)
		{
			fprintf(stdout, "Generating cross correlation statistic...\n");
		}

		/* cc statistic */
		LAL_CALL( LALStochasticCrossCorrelationStatistic(&status, &ccStat, &ccIn, \
				epochsMatch), &status );

		y = (REAL8)(ccStat.value * pow(10.,ccStat.units.powerOfTen));
		yWhitenedSum += y * inVarTheo;

		if (PRINT)
		{
			fprintf(stdout, "y = %f\n\n", y);
		}

		/* output to file */
		fprintf(out, "%e %e\n", sqrt(varTheo), y);
	}

	/* compute point estimate, error bar and report results */
	errorBar = sqrt(1./inVarTheoSum) / (REAL8)segmentDuration;
	pointEstimate = yWhitenedSum / (inVarTheoSum * (REAL8)segmentDuration);
	fprintf(stdout, "error bar = %f s\n", errorBar);
	fprintf(stdout, "h_0^2 omega = %f\n\n", pointEstimate);

	/* cleanup */
	LAL_CALL( LALDestroyRealFFTPlan(&status, &(specparPSD.plan)), &status );
	LAL_CALL( LALDestroyRealFFTPlan(&status, &fftDataPlan), &status );
	LAL_CALL( LALDestroyVector(&status, &(streamOne.data)), &status );
	LAL_CALL( LALDestroyVector(&status, &(streamTwo.data)), &status );
	LAL_CALL( LALDestroyVector(&status, &(segmentOne.data)), &status );
	LAL_CALL( LALDestroyVector(&status, &(segmentTwo.data)), &status );
	LAL_CALL( LALDestroyVector(&status, &(psdTempOne.data)), &status );
	LAL_CALL( LALDestroyVector(&status, &(psdTempTwo.data)), &status );
	LAL_CALL( LALDestroyVector(&status, &(psdOne.data)), &status );
	LAL_CALL( LALDestroyVector(&status, &(psdTwo.data)), &status );
	LAL_CALL( LALCDestroyVector(&status, &(responseOne.data)), &status );
	LAL_CALL( LALCDestroyVector(&status, &(responseTwo.data)), &status );
	LAL_CALL( LALCDestroyVector(&status, &(responseTempOne.data)), &status );
	LAL_CALL( LALCDestroyVector(&status, &(responseTempTwo.data)), &status );
	LAL_CALL( LALDestroyVector(&status, &(calInvPSDOne.data)), &status );
	LAL_CALL( LALDestroyVector(&status, &(calInvPSDTwo.data)), &status );
	LAL_CALL( LALCDestroyVector(&status, &(halfCalPSDOne.data)), &status );
	LAL_CALL( LALCDestroyVector(&status, &(halfCalPSDTwo.data)), &status );
	LAL_CALL( LALDestroyVector(&status, &(overlap.data)), &status );
	LAL_CALL( LALDestroyVector(&status, &(omegaGW.data)), &status );
	LAL_CALL( LALDestroyVector(&status, &(dataWindow.data)), &status );
	/*
	LAL_CALL( LALDestroyVector(&status, &(mask.data)), &status );
	LAL_CALL( LALDestroyVector(&status, &maskTemp), &status );
	LAL_CALL( LALDestroyVector(&status, &hannWindow), &status );
	*/
	LAL_CALL( LALCDestroyVector(&status, &(hBarTildeOne.data)), &status );
	LAL_CALL( LALCDestroyVector(&status, &(hBarTildeTwo.data)), &status );

	/* close output file */
	fclose(out);

	return 0;
}

/* parse command line options */
void parseOptions(INT4 argc, CHAR *argv[])
{
	while(1)
	{
		INT4 c = -1;
		c = getopt(argc, argv, "ht:T:l:a:f:F:w:i:I:d:D:r:R:mo:g:kvz:");
		if (c == -1)
		{
			/* end of options, break loop */
			break;
		}

		switch(c)
		{
			case 'h':
				/* HELP!!! */
				displayUsage(0);
				break;

			case 't':
				/* start time */
				startTime = atoi(optarg);
				break;

			case 'T':
				/* stream duration */
				streamDuration = atoi(optarg);
				break;

			case 'l':
				/* duration */
				segmentDuration = atoi(optarg);
				break;

			case 'a':
				/* resampling */
				resampleRate = atoi(optarg);
				break;

			case 'f':
				/* minimal frequency */
				fMin = atoi(optarg);
				break;

			case 'F':
				/* maximal frequency */
				fMax = atoi(optarg);
				break;

			case 'w':
				/* hann window duration */
				hannDuration = atoi(optarg);
				break;

			case 'i':
				/* ifo for first stream */
				strncpy(ifoOne, optarg, LALNameLength);

				/* set site and channel */
				if (strncmp(ifoOne, "H1", 2) == 0)
				{
					siteOne = 0;
					strncpy(channelOne, "H1:LSC-AS_Q", LALNameLength);
				}
				else if (strncmp(ifoOne, "H2", 2) == 0)
				{
					siteOne = 0;
					strncpy(channelOne, "H2:LSC-AS_Q", LALNameLength);
				}
				else if (strncmp(ifoOne, "L1", 2) == 0)
				{
					siteOne = 1;
					strncpy(channelOne, "L1:LSC-AS_Q", LALNameLength);
				}
				else
				{
					fprintf(stderr, "First IFO not recognised...\n");
					exit(1);
				}

				break;

			case 'I':
				/* ifo for second stream */
				strncpy(ifoTwo, optarg, LALNameLength);

				/* set site and channel */
				if (strncmp(ifoTwo, "H1", 2) == 0)
				{
					siteTwo = 0;
					strncpy(channelTwo, "H1:LSC-AS_Q", LALNameLength);
				}
				else if (strncmp(ifoTwo, "H2", 2) == 0)
				{
					siteTwo = 0;
					strncpy(channelTwo, "H2:LSC-AS_Q", LALNameLength);
				}
				else if (strncmp(ifoTwo, "L1", 2) == 0)
				{
					siteTwo = 1;
					strncpy(channelTwo, "L1:LSC-AS_Q", LALNameLength);
				}
				else
				{
					fprintf(stderr, "Second IFO not recognised...\n");
					exit(1);
				}

				break;

			case 'd':
				/* data cache one */
				strncpy(dataCacheOne, optarg, LALNameLength);
				break;

			case 'D':
				/* data cache two */
				strncpy(dataCacheTwo, optarg, LALNameLength);
				break;

			case 'r':
				/* calibration cache one */
				strncpy(calCacheOne, optarg, LALNameLength);
				break;

			case 'R':
				/* calibration cache two */
				strncpy(calCacheTwo, optarg, LALNameLength);
				break;

			case 'm':
				/* injection */
				simulData = 1;
				break;

			case 'o':
				/* scale factor */
				scaleFactor = atof(optarg);
				break;

			case 'g':
				/* seed */
				seed = atoi(optarg);
				break;

			case 'k':
				/* frequency masking */
				applyMask = 1;
				break;

			case 'v':
				/* verbose */
				PRINT = 1;
				break;

			case 'z':
				/* set debug level */
				set_debug_level( optarg );
				break;

			default:
				displayUsage(1);
		}
	}

	if (optind < argc)
	{
		displayUsage(1);
	}

	return;
}

/* display program usage */
void displayUsage(INT4 exitcode)
{
	fprintf(stderr, "Usage: pipeline [options]\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, " -h                  print this message\n");
	fprintf(stderr, " -t startTime        GPS start time\n");
	fprintf(stderr, " -T streamDuration   duration of stream\n");
	fprintf(stderr, " -l segmentDuration  segment duration\n");
	fprintf(stderr, " -a resampleRate     resample rate\n");
	fprintf(stderr, " -f fMin             minimal frequency\n");
	fprintf(stderr, " -F fMax             maximal frequency\n");
	fprintf(stderr, " -w hannDuration     hann duration\n");
	fprintf(stderr, " -i ifoOne           ifo for first stream\n");
	fprintf(stderr, " -I ifoTwo           ifo for second stream\n");
	fprintf(stderr, " -d dataCacheOne     cache file for first stream\n");
	fprintf(stderr, " -D dataCacheTwo     cache file for second stream\n");
	fprintf(stderr, " -r calibCacheOne    first stream calibration cache\n");
	fprintf(stderr, " -R calibCacheTwo    second stream calibration cache\n");
	fprintf(stderr, " -m                  inject a signal into the data\n");
	fprintf(stderr, " -o scaleFactor      scale factor for injection\n");
	fprintf(stderr, " -g seed             seed\n");
	fprintf(stderr, " -k                  apply frequency masking\n");
	fprintf(stderr, " -v                  verbose mode\n");
	fprintf(stderr, " -z debugLevel       set lalDebugLevel\n");

	exit(exitcode);
}
