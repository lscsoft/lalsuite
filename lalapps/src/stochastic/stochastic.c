/*
 * stochastic.c - SGWB Standalone Analysis Pipeline
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
#include <lal/IIRFilter.h>
#include <lal/BandPassTimeSeries.h>
#include "stochastic.h"

NRCSID (STOCHASTICC, "$Id$");
RCSID ("$Id$");

/* cvs info */
#define CVS_ID "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "lalapps_stochastic"

/* variables for getopt options parsing */
char *optarg;
int optind;

/* flags for getopt_long */
static int inject_flag = 0;
static int apply_mask_flag = 0;
static int high_pass_flag = 0;
static int overlap_hann_flag = 0;
static int verbose_flag = 0;
static int condor_flag = 0;

/* parameters for the stochastic search */

/* sampling parameters */
INT4 sampleRate = 16384;
INT4 resampleRate = 1024;
REAL8 deltaF = 0.25;

/* data parameters */
LIGOTimeGPS gpsStartTime;
UINT8 startTime = 729332041;
UINT8 stopTime = 729332101;
INT4 segmentDuration = 60;
INT4 calibDuration = 60;
CHAR frameCacheOne[100] = "cachefiles/H-729332041.cache";
CHAR frameCacheTwo[100] = "cachefiles/L-729332041.cache";
CHAR calCacheOne[100] = "calibration/H1-CAL-V03-729273600-734367600.cache";
CHAR calCacheTwo[100] = "calibration/L1-CAL-V03-729273600-734367600.cache";
CHAR channelOne[LALNameLength]= "H1:LSC-AS_Q";
CHAR channelTwo[LALNameLength]= "L1:LSC-AS_Q";
CHAR ifoOne[LALNameLength] = "H1";
CHAR ifoTwo[LALNameLength] = "L1";
INT4 siteOne = 1;
INT4 siteTwo = 0;

/* frequency band */
INT4 fMin = 64;
INT4 fMax = 265;

/* omegaGW parameters */
REAL4 alpha = 0.0;
REAL4 fRef = 100.0;
REAL4 omegaRef = 1.;

/* monte carlo parameters */
/* at the moment the code cannot do monte carlo with overlapped Hann window */
REAL4 scaleFactor = 1.;
INT4 seed = 173;
INT4 NLoop = 1;

/* window parameters */
/* 60 s for pure Hann, 1 s for Tukey, 0s for rectangular window */
INT4 hannDuration = 1;

/* high pass filtering parameters */
REAL4 highPassFreq = 40.;
REAL4 highPassAt = 0.707;
INT4  highPassOrder = 6;

/* set as 1 if data are downsampled, 0 otherwise */ 
INT4 padData = 1;


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
	REAL8 varTheo;

	/* input data segment */
	INT4 numSegments;
	INT4 segmentLength;
        INT4 segmentShift;
	ReadDataPairParams streamParams;
	StreamPair streamPair;
	REAL4TimeSeries segmentOne, segmentTempOne;
	REAL4TimeSeries segmentTwo, segmentTempTwo;

	/* simulated signal structures */
        StochasticOmegaGWParameters  parametersOmega;
	SSSimStochBGParams                 SBParams;
        SSSimStochBGInput                  SBInput; 
        SSSimStochBGOutput                 SBOutput;
        REAL4TimeSeries                    SimStochBGOne;
        REAL4TimeSeries                    SimStochBGTwo;
  
        REAL4FrequencySeries               MComegaGW;
        COMPLEX8FrequencySeries            MCresponseOne;
        COMPLEX8FrequencySeries            MCresponseTwo;
       
        INT4 MCLoop;
        INT4 MCfreqLength;
        REAL8 MCdeltaF, MCdeltaT;
        LALUnit countPerStrain = {0,{0,0,0,0,0,-1,1},{0,0,0,0,0,0,0}};
    
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

        /* high pass filtering */
        PassBandParamStruc highpassParam;

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

        /* for Condor use, read parameters into input file */
        if (condor_flag == 1)
        { 
                fscanf(stdin,"%d\n%d\n%s\n%s\n",&startTime, &stopTime,&frameCacheOne,&frameCacheTwo);
        }
        
	if (verbose_flag)
	{
		fprintf(stdout, "Calculating number of segments...\n");
	}

	/* get number of segments */
	numSegments = (INT4)((stopTime - startTime) / segmentDuration );
        segmentShift = segmentDuration;

        if (overlap_hann_flag)
	{
	        numSegments = 2 * numSegments - 1;
                segmentShift = segmentDuration / 2;
        }

        /* set length for data segments */
	segmentLength = segmentDuration * resampleRate;

	/* set metadata fields for data segments */
	strncpy(segmentTempOne.name, "segmentTempOne", LALNameLength);
	strncpy(segmentTempTwo.name, "segmentTempTwo", LALNameLength);
	segmentTempOne.sampleUnits = lalADCCountUnit;
	segmentTempTwo.sampleUnits = lalADCCountUnit;
	segmentTempOne.epoch = gpsStartTime;
	segmentTempTwo.epoch = gpsStartTime;
	segmentTempOne.deltaT = 1./(REAL8)resampleRate;
	segmentTempTwo.deltaT = 1./(REAL8)resampleRate;
	segmentTempOne.f0 = 0;
	segmentTempTwo.f0 = 0;
     
        strncpy(segmentOne.name, "segmentOne", LALNameLength);
	strncpy(segmentTwo.name, "segmentTwo", LALNameLength);
	segmentOne.sampleUnits = lalADCCountUnit;
	segmentTwo.sampleUnits = lalADCCountUnit;
	segmentOne.epoch = gpsStartTime;
	segmentTwo.epoch = gpsStartTime;
	segmentOne.deltaT = 1./(REAL8)resampleRate;
	segmentTwo.deltaT = 1./(REAL8)resampleRate;
	segmentOne.f0 = 0;
	segmentTwo.f0 = 0;

	if (verbose_flag)
	{
		fprintf(stdout, "Allocating memory for data segments...\n");
	}

	/* allocate memory for data segments */
	segmentTempOne.data = NULL;
	segmentTempTwo.data = NULL;
	LAL_CALL( LALSCreateVector(&status, &(segmentTempOne.data), segmentLength), \
			&status );
	LAL_CALL( LALSCreateVector(&status, &(segmentTempTwo.data), segmentLength), \
			&status );
	memset(segmentTempOne.data->data, 0, \
			segmentTempOne.data->length * sizeof(*segmentTempOne.data->data));
	memset(segmentTempTwo.data->data, 0, \
			segmentTempTwo.data->length * sizeof(*segmentTempTwo.data->data));

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

	/* set segment input parameters */
	streamParams.duration = segmentDuration;
	streamParams.frameCacheOne = frameCacheOne;
	streamParams.frameCacheTwo = frameCacheTwo;
	streamParams.ifoOne = ifoOne;
	streamParams.ifoTwo = ifoTwo;
	streamParams.channelOne = channelOne;
	streamParams.channelTwo = channelTwo;
	streamParams.startTime = startTime;
	streamParams.buffer = padData;
	streamParams.sampleRate = sampleRate;
	streamParams.resampleRate = resampleRate;

	/* set stream data structures */
	streamPair.streamOne = &segmentTempOne;
	streamPair.streamTwo = &segmentTempTwo;

	if (inject_flag)
	{
		
		if (verbose_flag)
		{
			fprintf(stdout, "Allocating memory for MC...\n");
		}
                MCdeltaT = 1.0 / resampleRate;
                /*MCdeltaF = 1.0 / (MCdeltaT * segmentLength);*/
                MCdeltaF = (REAL8)resampleRate / (REAL8)segmentLength;
                MCfreqLength = segmentLength / 2 + 1;

		/* create vectors to store the simulated signal */
                strncpy(SimStochBGOne.name, "Whitened-SimulatedSBOne",LALNameLength);
		strncpy(SimStochBGTwo.name, "Whitened-SimulatedSBTwo",LALNameLength );
                SimStochBGOne.f0 = 0.;
		SimStochBGTwo.f0 = 0.;
                SimStochBGOne.epoch = gpsStartTime;
	        SimStochBGTwo.epoch = gpsStartTime;
		SimStochBGOne.deltaT = 1./(REAL8)resampleRate;
		SimStochBGTwo.deltaT = 1./(REAL8)resampleRate;
		SimStochBGOne.sampleUnits = lalADCCountUnit;
		SimStochBGTwo.sampleUnits = lalADCCountUnit;

		SimStochBGOne.data = NULL;
		SimStochBGTwo.data = NULL;
		LAL_CALL(LALSCreateVector(&status, &(SimStochBGOne.data),segmentLength), &status);
		LAL_CALL(LALSCreateVector(&status, &(SimStochBGTwo.data),segmentLength), &status);
	
		memset(SimStochBGOne.data->data, 0,SimStochBGOne.data->length *sizeof(*SimStochBGOne.data->data));
		memset(SimStochBGTwo.data->data, 0,SimStochBGTwo.data->length *sizeof(*SimStochBGTwo.data->data));

		/* define parameters for SimulateSB */

		SBParams.length = segmentLength;
		SBParams.deltaT = 1. / resampleRate;
		SBParams.detectorOne = lalCachedDetectors[siteOne];
		SBParams.detectorTwo = lalCachedDetectors[siteTwo];
		SBParams.SSimStochBGTimeSeries1Unit = lalADCCountUnit;
		SBParams.SSimStochBGTimeSeries2Unit = lalADCCountUnit;

		/* omegaGW */

		parametersOmega.length = MCfreqLength;
		parametersOmega.f0 = 0.;
		parametersOmega.deltaF = MCdeltaF;
		parametersOmega.alpha = alpha;
		parametersOmega.fRef = fRef;
		parametersOmega.omegaRef = omegaRef;

		/* allocate memory */
		MComegaGW.data = NULL;
		LAL_CALL(LALSCreateVector(&status, &(MComegaGW.data), MCfreqLength), &status);
       
		memset(MComegaGW.data->data, 0, \
			MComegaGW.data->length * sizeof(*MComegaGW.data->data));

		/* generate omegaGW */
		LAL_CALL(LALStochasticOmegaGW(&status, &MComegaGW, &parametersOmega), &status);
       
		/* response functions */
		
		strncpy(MCresponseOne.name,"MCresponseOne", LALNameLength);
		strncpy(MCresponseTwo.name,"MCresponseTwo", LALNameLength);
		MCresponseOne.sampleUnits = countPerStrain;
		MCresponseTwo.sampleUnits = countPerStrain;
		MCresponseOne.epoch = gpsStartTime;
		MCresponseTwo.epoch = gpsStartTime;
		MCresponseOne.deltaF = MCdeltaF;
		MCresponseTwo.deltaF = MCdeltaF;
		MCresponseOne.f0 = 0;
		MCresponseTwo.f0 = 0;

		/* allocate memory */
		MCresponseOne.data = NULL;
		MCresponseTwo.data = NULL;
		LAL_CALL(LALCCreateVector(&status, &(MCresponseOne.data), MCfreqLength), &status);
	
		LAL_CALL(LALCCreateVector(&status, &(MCresponseTwo.data), MCfreqLength), &status);
   
		memset(MCresponseOne.data->data, 0, \
		       MCresponseOne.data->length * sizeof(*MCresponseOne.data->data));
		memset(MCresponseTwo.data->data, 0, \
		       MCresponseTwo.data->length * sizeof(*MCresponseTwo.data->data));

	}

	

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

	if (verbose_flag)
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

	if (verbose_flag)
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

	if (verbose_flag)
	{
		fprintf(stdout, "Creating FFT plan for PSD estimation...\n");
	}

	/* create fft plan */
	LAL_CALL (LALCreateForwardRealFFTPlan(&status, &specparPSD.plan, \
				windowPSDLength, 0), &status );

	if (verbose_flag)
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

	/* set metadata fields for response functions */
	strncpy(responseTempOne.name, "responseTempOne", LALNameLength);
	strncpy(responseTempTwo.name, "responseTempTwo", LALNameLength);
	responseTempOne.sampleUnits = countPerAttoStrain;
	responseTempTwo.sampleUnits = countPerAttoStrain;
	responseTempOne.epoch = gpsStartTime;
	responseTempTwo.epoch = gpsStartTime;
	responseTempOne.deltaF = deltaF;
	responseTempTwo.deltaF = deltaF;
	responseTempOne.f0 = 0;
	responseTempTwo.f0 = 0;

	if (verbose_flag)
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
	responseOne.epoch = gpsStartTime;
	responseTwo.epoch = gpsStartTime;
	responseOne.deltaF = deltaF;
	responseTwo.deltaF = deltaF;
	responseOne.f0 = fMin;
	responseTwo.f0 = fMin;

	if (verbose_flag)
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

	if (verbose_flag)
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

	if (verbose_flag)
	{
		fprintf(stdout, "Allocating memory for data segment window...\n");
	}

	/* allocate memory for segment window */
	dataWindow.data = NULL;
	LAL_CALL( LALSCreateVector(&status, &(dataWindow.data), segmentLength), \
			&status );
	memset(dataWindow.data->data, 0, \
			dataWindow.data->length * sizeof(*dataWindow.data->data));

	if (verbose_flag)
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
	if (verbose_flag)
	{
		LALSPrintTimeSeries(&dataWindow, "dataWindow.dat");
	}
        
        /* structure for high pass filtering */
        if (high_pass_flag)
	{
	  highpassParam.nMax = highPassOrder;
	  highpassParam.f1 = -1.0;
	  highpassParam.f2 = highPassFreq;
	  highpassParam.a1 = -1.0;
	  highpassParam.a2 = highPassAt;
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

	if (verbose_flag)
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

	if (verbose_flag)
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

	if (verbose_flag)
	{
		fprintf(stdout, "Generating the overlap reduction function...\n");
	}

	/* generate overlap reduction function */
	LAL_CALL (LALOverlapReductionFunction(&status, &overlap, &detectors, \
				&ORFparams), &status );

	/* save */
	if (verbose_flag)
	{
		LALSPrintFrequencySeries(&overlap, "overlap.dat");
	}

	/* set metadata fields for spectrum */
	strncpy(omegaGW.name, "omegaGW", LALNameLength);
	omegaGW.sampleUnits = lalDimensionlessUnit;
	omegaGW.deltaF = deltaF;
	omegaGW.f0 = fMin;

	if (verbose_flag)
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

	if (verbose_flag)
	{
		fprintf(stdout, "Generating spectrum for optimal filter...\n");
	}

	/* generage omegaGW */
	LAL_CALL( LALStochasticOmegaGW(&status, &omegaGW, &omegaGWParams), &status );

	/* frequency mask */
	if (apply_mask_flag)
	{
		/* set metadata fields for frequency mask */
		strncpy(mask.name, "mask", LALNameLength);
		mask.deltaF = deltaF;
		mask.f0 = fMin;
		mask.sampleUnits = lalDimensionlessUnit;

		if (verbose_flag)
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

		if (verbose_flag)
		{
			fprintf(stdout, "Generating frequency mask...\n");
		}

		/* set all values to 1 */
		for (i = 0; i < respLength; i++)
		{
			maskTemp->data[i] = 1.;
		}

		if (verbose_flag)
		{
			fprintf(stdout, "Masking multiples of 16 Hz...\n");
		}

		/* remove multiples of 16 Hz */
		for (i = 0; i < respLength; i += (UINT4)(16 / deltaF))
		{
			maskTemp->data[i] = 0.;
		}

		if (verbose_flag)
		{
			fprintf(stdout, "Masking multiples of 60 Hz...\n");
		}

		/* remove multiples of 60 Hz */
		for (i = 0; i < respLength; i += (UINT4)(60 / deltaF))
		{
			maskTemp->data[i] = 0.;
		}

		if (verbose_flag)
		{
			fprintf(stdout, "Getting appropriate frequency band for mask...\n");
		}

		/* get appropriate band */
		for (i = 0; i < filterLength; i++)
		{
			mask.data->data[i] = maskTemp->data[i + numPointInf];
		}

		if (verbose_flag)
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
	if (verbose_flag)
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
	optFilter.epoch = gpsStartTime;
	optFilter.deltaF = deltaF;
	optFilter.f0 = fMin;

	if (verbose_flag)
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
	ccSpectrum.epoch = gpsStartTime;
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
	 ** loop over segments **
	 */

        if (verbose_flag)
	{
		fprintf(stdout, "Looping over %d segments...\n", numSegments);
	}
        
        lal_errhandler = LAL_ERR_RTRN;

	for (segLoop = 0; segLoop < numSegments; segLoop++)
	{
                                
		/* define segment epoch */
		gpsStartTime.gpsSeconds = startTime + (segLoop * segmentShift);
		segmentOne.epoch = gpsStartTime;
		segmentTwo.epoch = gpsStartTime;

		if (verbose_flag)
		{
			fprintf(stdout, "Performing search on segment %d of %d...\n", \
					segLoop + 1, numSegments);
		}
                
                /* compute response function */

                if (verbose_flag)
		{
			fprintf(stdout, "Getting appropriate frequencresponse " \
					"functions...\n");
		}

		responseTempOne.epoch = gpsStartTime;
		responseTempTwo.epoch = gpsStartTime;
                
		LAL_CALL( LALExtractFrameResponse(&status, &responseTempOne, calCacheOne, \
					ifoOne, &duration), &status );

		if ((status.statusCode !=0)||(responseTempOne.data==NULL))
		  {
		    clear_status(&status);
                    if (segLoop < (numSegments - 1))
		     continue; 
                    else
		     break;   
		  }     
                         
               
		LAL_CALL( LALExtractFrameResponse(&status, &responseTempTwo, calCacheTwo, \
					ifoTwo, &duration), &status );
                
                if ((status.statusCode !=0)||(responseTempTwo.data==NULL))
		  {
		    clear_status(&status);
                    if (segLoop < (numSegments - 1))
		     continue; 
                    else
		     break;   
		  }

                if (verbose_flag)
		{
			fprintf(stdout, "Generating response functions...\n");
		}

		

		/* reduce to the optimal filter frequency range */
		responseOne.epoch = gpsStartTime;
		responseTwo.epoch = gpsStartTime;
		for (i = 0; i < filterLength; i++)
		{
			responseOne.data->data[i] = responseTempOne.data->data[i + numPointInf];
			responseTwo.data->data[i] = responseTempTwo.data->data[i + numPointInf];
		}

		/* output the results */
		if (verbose_flag)
		{
			LALCPrintFrequencySeries(&responseOne, "response1.dat");
			LALCPrintFrequencySeries(&responseTwo, "response2.dat");
		}

		/* compute response function for MC*/
                if (inject_flag)
		{
         	     if (verbose_flag)
        	 	{
	        		fprintf(stdout, "Getting appropriate frequencyresponse " \
					"functions for MC...\n");
        		}

        		MCresponseOne.epoch = gpsStartTime;
        		MCresponseTwo.epoch = gpsStartTime;
                
			LAL_CALL( LALExtractFrameResponse(&status, &MCresponseOne, calCacheOne, \
					ifoOne, &duration), &status );
			
		    
               
			LAL_CALL( LALExtractFrameResponse(&status, &MCresponseTwo, calCacheTwo, \
					ifoTwo, &duration), &status );

                        /* force DC to be 0 and nyquist to be real */
        	       	MCresponseOne.data->data[0].re = 0.;
        		MCresponseTwo.data->data[0].re = 0.;
        		MCresponseOne.data->data[0].im = 0.;
        		MCresponseTwo.data->data[0].im = 0.;
        		MCresponseOne.data->data[MCfreqLength-1].im = 0.;
        		MCresponseTwo.data->data[MCfreqLength-1].im = 0.;
   
                       	/* output the results */
		if (verbose_flag)
		{
			LALCPrintFrequencySeries(&MCresponseOne, "MCresponse1.dat");
			LALCPrintFrequencySeries(&MCresponseTwo, "MCresponse2.dat");
		}
                
		}

        	/* read data and downsample */
        	if (verbose_flag)
        	{
        		fprintf(stdout, "Reading data...\n");
        	}

        	/* read data */
                streamParams.startTime = gpsStartTime.gpsSeconds;
        	LAL_CALL(readDataPair(&status, &streamPair, &streamParams), &status);
                  
                if ((status.statusCode !=0)||(segmentTempOne.data==NULL)||(segmentTempTwo.data==NULL))
		  {
		    clear_status(&status);
                    if (segLoop < (numSegments - 1))
		     continue; 
                    else
		     break;   
		  }
               
        	/* save */
        	if (verbose_flag)
        	{
        		LALSPrintTimeSeries(&segmentTempOne, "segment1.dat");
         		LALSPrintTimeSeries(&segmentTempTwo, "segment2.dat");
        	}

	        for (MCLoop = 0; MCLoop < NLoop; MCLoop ++)
		{	
       		/* open output file */
		LALSnprintf(outputFilename, LALNameLength, "output/stoch-%s%s-%d-%d-%d.dat\
", ifoOne, ifoTwo, (INT4)startTime, (INT4)stopTime, MCLoop);

		/* simulate signal */
		if (inject_flag)
		{
		        /* set parameters for monte carlo */
		        SimStochBGOne.epoch = gpsStartTime;
        		SimStochBGTwo.epoch = gpsStartTime;
      			/* seed = (INT4)(time(NULL)); */
                        SBParams.seed = seed ;
                        
                        /* define input structure for SimulateSB */
			SBInput.omegaGW = &MComegaGW;
			SBInput.whiteningFilter1 = &MCresponseOne;
			SBInput.whiteningFilter2 = &MCresponseTwo;

                        /* define output structure for SimulateSB */
			SBOutput.SSimStochBG1 = &SimStochBGOne;
			SBOutput.SSimStochBG2 = &SimStochBGTwo;
                        
                       /* perform monte carlo */
                        LALSSSimStochBGTimeSeries(&status, &SBOutput, &SBInput, &SBParams);
			/* output the results */
			if (verbose_flag)
			{
				LALSPrintTimeSeries(&SimStochBGOne, "simStochBG1.dat");
				LALSPrintTimeSeries(&SimStochBGTwo, "simStochBG2.dat");
			}

			/* multiply by scale factor and inject into real data */
			for (i = 0; i < segmentLength ; i++)
			{
				segmentOne.data->data[i] = segmentTempOne.data->data[i] + \
					(scaleFactor * SimStochBGOne.data->data[i]);
				segmentTwo.data->data[i] = segmentTempTwo.data->data[i] + \
					(scaleFactor * SimStochBGTwo.data->data[i]);
			}
                
			/* increase seed */
			seed = seed + 2 ;
		}
                else
		{
		  for (i = 0; i < segmentLength; i ++)
                  {
		   segmentOne.data->data[i] = segmentTempOne.data->data[i];
		   segmentTwo.data->data[i] = segmentTempTwo.data->data[i] ;
		   }
		 }
		/* output the results */
		if (verbose_flag)
		{
			LALSPrintTimeSeries(&segmentOne, "segment1.dat");
			LALSPrintTimeSeries(&segmentTwo, "segment2.dat");
		}
                
                if (high_pass_flag)
		{
                         
                         LAL_CALL( LALButterworthREAL4TimeSeries( &status, &segmentOne, &highpassParam ), &status );
                         LAL_CALL( LALButterworthREAL4TimeSeries( &status, &segmentTwo, &highpassParam ), &status );
                }
		/* zero pad and fft */
		LAL_CALL( LALSZeroPadAndFFT(&status, &hBarTildeOne, &segmentOne, \
					&zeroPadParams), &status );
		LAL_CALL( LALSZeroPadAndFFT(&status, &hBarTildeTwo, &segmentTwo, \
					&zeroPadParams), &status );

		/* save */
		if (verbose_flag)
		{
			LALCPrintFrequencySeries(&hBarTildeOne, "hBarTilde1.dat");
			LALCPrintFrequencySeries(&hBarTildeTwo, "hBarTilde2.dat");
		}

		if (verbose_flag)
		{
			fprintf(stdout, "Estimating PSDs...\n");
		}

		/* compute uncalibrated PSDs */
		LAL_CALL( LALREAL4AverageSpectrum(&status, &psdTempOne, &segmentOne, \
					&specparPSD), &status );
		LAL_CALL( LALREAL4AverageSpectrum(&status, &psdTempTwo, &segmentTwo, \
					&specparPSD), &status );

		if (verbose_flag)
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
		if (verbose_flag)
		{
			LALSPrintFrequencySeries(&psdOne, "psd1.dat");
			LALSPrintFrequencySeries(&psdTwo, "psd2.dat");
		}


		if (verbose_flag)
		{
			fprintf(stdout, "Generating inverse noise...\n");
		}

		/* compute inverse calibrate and half calibrated PSDs */
		LAL_CALL( LALStochasticInverseNoise(&status, &inverseNoiseOutOne, \
					&inverseNoiseInOne), &status );
		LAL_CALL( LALStochasticInverseNoise(&status, &inverseNoiseOutTwo, \
					&inverseNoiseInTwo), &status );

		/* save */
		if (verbose_flag)
		{
			LALSPrintFrequencySeries(&calInvPSDOne, "inPSD1.dat");
			LALSPrintFrequencySeries(&calInvPSDTwo, "inPSD2.dat");
			LALCPrintFrequencySeries(&halfCalPSDOne, "hInPSD1.dat");
			LALCPrintFrequencySeries(&halfCalPSDTwo, "hInPSD2.dat");
		}

		if (verbose_flag)
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
		

		/* save */
		if (verbose_flag)
		{
			fprintf(stdout, "lambda = %e s^-1\n", lambda);
			fprintf(stdout, "varTheo = %e s\n", varTheo);
		}

		if (verbose_flag)
		{
			fprintf(stdout, "Generating optimal filter...\n");
		}

		/* build optimal filter */
		optFilter.epoch = gpsStartTime;
		LAL_CALL( LALStochasticOptimalFilter(&status, &optFilter, &optFilterIn, \
					&normLambda), &status );

		/* save */
		if (verbose_flag)
		{
			LALCPrintFrequencySeries(&optFilter, "optFilter.dat");
		}

		if (verbose_flag)
		{
			fprintf(stdout, "Generating cross correlation spectrum...\n");
		}

		/* cc spectrum */
		LAL_CALL( LALStochasticCrossCorrelationSpectrum(&status, &ccSpectrum, \
					&ccIn, epochsMatch), &status );

		/* save */
		if (verbose_flag)
		{
			LALCPrintFrequencySeries(&ccSpectrum, "ccSpectrum.dat");
		}

		if (verbose_flag)
		{
			fprintf(stdout, "Generating cross correlation statistic...\n");
		}

		/* cc statistic */
		LAL_CALL( LALStochasticCrossCorrelationStatistic(&status, &ccStat, &ccIn, \
					epochsMatch), &status );

		y = (REAL8)(ccStat.value * pow(10.,ccStat.units.powerOfTen));

		if (verbose_flag)
		{
			fprintf(stdout, "y = %f\n", y);
		}

                if (status.statusCode !=0)
		  {
		    clear_status(&status);
                    if (MCLoop < (NLoop - 1))
		     continue; 
                    else
		     break;   
		  }

		/* output to file */
                
                out = fopen(outputFilename, "a");
	       	fprintf(out,"%d %e %e\n", gpsStartTime.gpsSeconds, y, varTheo);
                fclose(out);
		}
	}
       
        lal_errhandler = LAL_ERR_EXIT;

	/* cleanup */

	LAL_CALL( LALDestroyRealFFTPlan(&status, &(specparPSD.plan)), &status );
	LAL_CALL( LALDestroyRealFFTPlan(&status, &fftDataPlan), &status );
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
	if (apply_mask_flag)
	{
		LAL_CALL( LALDestroyVector(&status, &(mask.data)), &status );
		LAL_CALL( LALDestroyVector(&status, &maskTemp), &status );
	}
	if (hannDuration != 0)
	{
		LAL_CALL( LALDestroyVector(&status, &hannWindow), &status );
	}
	LAL_CALL( LALCDestroyVector(&status, &(hBarTildeOne.data)), &status );
	LAL_CALL( LALCDestroyVector(&status, &(hBarTildeTwo.data)), &status );
	if (inject_flag)
	{       
		LAL_CALL( LALDestroyVector(&status, &SimStochBGOne.data), &status );
		LAL_CALL( LALDestroyVector(&status, &SimStochBGTwo.data), &status );
		LAL_CALL( LALCDestroyVector(&status, &(MCresponseOne.data)), &status );
		LAL_CALL( LALCDestroyVector(&status, &(MCresponseTwo.data)), &status );
                LAL_CALL( LALDestroyVector(&status, &(MComegaGW.data)), &status ); 
	}
	

	return 0;
}


/* parse command line options */
void parseOptions(INT4 argc, CHAR *argv[])
{
	int c = -1;

	while(1)
	{
		static struct option long_options[] =
		{
			/* options that set a flag */
			{"inject", no_argument, &inject_flag, 1},
			{"apply-mask", no_argument, &apply_mask_flag, 1},
                        {"high-pass-filter", no_argument, &high_pass_flag, 1},
			{"overlap-hann", no_argument, &overlap_hann_flag, 1},
			{"verbose", no_argument, &verbose_flag, 1},
			/* options that don't set a flag */
			{"help", no_argument, 0, 'h'},
			{"gps-start-time", required_argument, 0, 't'},
			{"gps-end-time", required_argument, 0, 'T'},
			{"segment-duration", required_argument, 0, 'l'},
			{"resample-rate", required_argument, 0, 'a'},
			{"f-min", required_argument, 0, 'f'},
			{"f-max", required_argument, 0, 'F'},
			{"hann-duration", required_argument, 0, 'w'},
                        {"hpf-frequency", required_argument, 0, 'k'},
                        {"hpf-attenuation", required_argument, 0, 'p'},
                        {"hpf-order", required_argument, 0, 'P'},
			{"ifo-one", required_argument, 0, 'i'},
			{"ifo-two", required_argument, 0, 'I'},
			{"frame-cache-one", required_argument, 0, 'd'},
			{"frame-cache-two", required_argument, 0, 'D'},
			{"calibration-cache-one", required_argument, 0, 'r'},
			{"calibration-cache-two", required_argument, 0, 'R'},
			{"scale-factor", required_argument, 0, 'o'},
			{"seed", required_argument, 0, 'g'},
                        {"number-of-injection", required_argument, 0, 'n'},
			{"debug-level", required_argument, 0, 'z'},
			{"version", no_argument, 0, 'V'},
			{0, 0, 0, 0}
		};

		/* getopt_long stores the option here */
		int option_index = 0;
		size_t optarg_len;

		c = getopt_long(argc, argv, "ht:T:l:a:f:F:w:k:p:P:i:I:d:D:r:R:o:g:n:z:V", \
				long_options, &option_index);

		if (c == -1)
		{
			/* end of options, break loop */
			break;
		}

		switch(c)
		{
			case 0:
				/* If this option set a flag, do nothing else now. */
				if (long_options[option_index].flag != 0)
					break;
				printf ("option %s", long_options[option_index].name);
				if (optarg)
					printf (" with arg %s", optarg);
				printf ("\n");
				break;

			case 'h':
				/* HELP!!! */
				displayUsage(0);
				break;

			case 't':
				/* start time */
				startTime = atoi(optarg);
				break;

			case 'T':
				/* stop time */
				stopTime = atoi(optarg);
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

                        case 'k':
				/* high pass knee filter frequency  */
				highPassFreq= atoi(optarg);
				break;
                          
                        case 'p':
				/* high pass filter attenuation  */
				highPassAt = atoi(optarg);
				break;

                        case 'P':
				/* high pass filter order  */
				highPassOrder = atoi(optarg);
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
        		  strncpy(frameCacheOne, optarg, LALNameLength);
        		  break;

        		case 'D':
        		  /* data cache two */
        		  strncpy(frameCacheTwo, optarg, LALNameLength);
        		  break;
   
        		case 'r':
         		  /* calibration cache one */
        		  strncpy(calCacheOne, optarg, LALNameLength);
        		  break;

        		case 'R':
        		  /* calibration cache two */
        		  strncpy(calCacheTwo, optarg, LALNameLength);
        		  break;
  
 			case 'o':
				/* scale factor */
				scaleFactor = atof(optarg);
				break;

			case 'g':
				/* seed */
				seed = atoi(optarg);
				break;

                        case 'N':
				/* number of injection */
				NLoop = atoi(optarg);
				break;

			case 'z':
				/* set debug level */
				set_debug_level( optarg );
				break;

			case 'V':
				/* display version info and exit */
				fprintf(stdout, "Standalone SGWB Search Engine\n" CVS_ID "\n");
				exit(0);
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
	fprintf(stderr, " -h                            print this message\n");
	fprintf(stderr, " -V                            display version\n");
	fprintf(stderr, " --verbose                     verbose mode\n");
	fprintf(stderr, " -z                            set lalDebugLevel\n");
	fprintf(stderr, " -t                            GPS start time\n");
	fprintf(stderr, " -T                            GPS stop time\n");
	fprintf(stderr, " -l                            segment duration\n");
	fprintf(stderr, " -a                            resample rate\n");
	fprintf(stderr, " -f                            minimal frequency\n");
	fprintf(stderr, " -F                            maximal frequency\n");
        fprintf(stderr, " -- high-pass-filter           apply high pass filter\n");
	fprintf(stderr, " -k                            high pass filter knee frequency\n");
       	fprintf(stderr, " -p                            high pass filter attenuation\n");
        fprintf(stderr, " -P                            high pass filter order\n");        
        fprintf(stderr, " --overlap-hann                use overlap window\n");                
        fprintf(stderr, " -w                            hann duration\n");
	fprintf(stderr, " -i                            ifo for first stream\n");
	fprintf(stderr, " -I                            ifo for second stream\n");
	fprintf(stderr, " -d                            cache file for first " \
			"stream\n");
	fprintf(stderr, " -D                            cache file for second " \
			"stream\n");
	fprintf(stderr, " -r                            first stream calibration cache\n");
	fprintf(stderr, " -R                            second stream calibration cache\n");
	fprintf(stderr, " --apply-mask                  apply frequency masking\n");
	fprintf(stderr, " --inject                      inject a signal into the data\n");
	fprintf(stderr, " -o                            scale factor for injection\n");
	fprintf(stderr, " -g                            seed\n");
        fprintf(stderr, " -N                            NLoop\n");
        
	exit(exitcode);
}
/*
void displayUsage(INT4 exitcode)
{
	fprintf(stderr, "Usage: pipeline [options]\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, " --help                        print this message\n");
	fprintf(stderr, " --version                     display version\n");
	fprintf(stderr, " --verbose                     verbose mode\n");
	fprintf(stderr, " --debug-level LEVEL           set lalDebugLevel\n");
	fprintf(stderr, " --gps-start-time SEC          GPS start time\n");
	fprintf(stderr, " --gps-end-time SEC            GPS stop time\n");
	fprintf(stderr, " --segment-duration SEC        segment duration\n");
	fprintf(stderr, " --resample-rate F             resample rate\n");
	fprintf(stderr, " --f-min F                     minimal frequency\n");
	fprintf(stderr, " --f-max F                     maximal frequency\n");
	fprintf(stderr, " --hpf-frequency F             high pass filter knee frequency\n");
       	fprintf(stderr, " --hpf-attenuation N           high pass filter attenuation\n");
        fprintf(stderr, " --hpf-order N                 high pass filter order\n");        
                        
        fprintf(stderr, " --hann-duration SEC           hann duration\n");
	fprintf(stderr, " --ifo-one IFO                 ifo for first stream\n");
	fprintf(stderr, " --ifo-two IFO                 ifo for second stream\n");
	fprintf(stderr, " --frame-cache-one FILE        cache file for first " \
			"stream\n");
	fprintf(stderr, " --frame-cache-two FILE        cache file for second " \
			"stream\n");
	fprintf(stderr, " --calibration-cache-one FILE  " \
			"first stream calibration cache\n");
	fprintf(stderr, " --calibration-cache-two FILE  " \
			"second stream calibration cache\n");
	fprintf(stderr, " --apply-mask                  apply frequency masking\n");
	fprintf(stderr, " --inject                      inject a signal into " \
			"the data\n");
	fprintf(stderr, " --scale-factor N              scale factor for "\
			"injection\n");
	fprintf(stderr, " --seed N                      seed\n");
        fprintf(stderr, " --number-of-injection N        NLoop\n");
	exit(exitcode);
}
*/
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
	INITSTATUS( status, "readDataPair", STOCHASTICC );
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

	if (verbose_flag)
	{
		fprintf(stdout, "Allocating memory for raw data streams...\n");
	}

	/* allocate memory */
	dataStreamOne.data = NULL;
	dataStreamTwo.data = NULL;
	LALSCreateVector(status->statusPtr, &(dataStreamOne.data), \
			sampleRate * (params->duration + (2 * buffer)));
	
	LALSCreateVector(status->statusPtr, &(dataStreamTwo.data), \
			sampleRate * (params->duration + (2 * buffer)));
	
	memset(dataStreamOne.data->data, 0, \
			dataStreamOne.data->length * sizeof(*dataStreamOne.data->data));
	memset(dataStreamTwo.data->data, 0, \
			dataStreamTwo.data->length * sizeof(*dataStreamTwo.data->data));

	if (verbose_flag)
	{
		fprintf(stdout, "Opening first frame cache...\n");
	}

	/* open first frame cache */
	LALFrCacheImport(status->statusPtr, &frCacheOne, params->frameCacheOne);
	
	LALFrCacheOpen(status->statusPtr, &frStreamOne, frCacheOne);
	

	if (verbose_flag)
	{
		fprintf(stdout, "Reading in channel \"%s\"...\n", frChanInOne.name);
	}

	/* read first channel */
	LALFrSeek(status->statusPtr, &(bufferStartTime), frStreamOne);
	
	LALFrGetREAL4TimeSeries(status->statusPtr, &dataStreamOne, \
			&frChanInOne, frStreamOne);
       

	if (strcmp(params->frameCacheOne, params->frameCacheTwo) == 0)
	{
		if (verbose_flag)
		{
			fprintf(stdout, "Reading in channel \"%s\" from same cache...\n", \
					frChanInTwo.name);
		}

		/* read in second channel */
		LALFrSeek(status->statusPtr, &(bufferStartTime), frStreamOne);
	    
		LALFrGetREAL4TimeSeries(status->statusPtr, &dataStreamTwo, \
				&frChanInTwo, frStreamOne);
		
		if (verbose_flag)
		{
			fprintf(stdout, "Closing frame cache...\n");
		}

		/* close frame cache */
		LALFrClose(status->statusPtr, &frStreamOne);
		
	}
	else
	{
		if (verbose_flag)
		{
			fprintf(stdout, "Closing first frame cache...\n");
		}

		/* close first frame cache */
		LALFrClose(status->statusPtr, &frStreamOne);
		

		if (verbose_flag)
		{
			fprintf(stdout, "Opening second frame cache...\n");
		}

		/* open second frame cache and read in second channel */
		LALFrCacheImport(status->statusPtr, &frCacheTwo, params->frameCacheTwo);
	    
		LALFrCacheOpen(status->statusPtr, &frStreamTwo, frCacheTwo);
		

		if (verbose_flag)
		{
			fprintf(stdout, "Reading in channel \"%s\"...\n", frChanInTwo.name);
		}

		/* read in second channel */
		LALFrSeek(status->statusPtr, &(bufferStartTime), frStreamTwo);
		
		LALFrGetREAL4TimeSeries(status->statusPtr, &dataStreamTwo, \
				&frChanInTwo, frStreamTwo);
		

		if (verbose_flag)
		{
			fprintf(stdout, "Closing second frame cache...\n");
		}

		/* close second frame stream */
		LALFrClose(status->statusPtr, &frStreamTwo);
		
	}

	/* resample */
	if (resampleRate != sampleRate)
	{
		if (verbose_flag)
		{
			fprintf(stdout, "Resampling to %d Hz...\n", resampleRate);
		}

		/* set resample parameters */
		resampleParams.deltaT = 1.0 / (REAL8)resampleRate;
		resampleParams.filterType = defaultButterworth;

		/* resample */
		LALResampleREAL4TimeSeries(status->statusPtr, &dataStreamOne, \
				&resampleParams);
		
		LALResampleREAL4TimeSeries(status->statusPtr, &dataStreamTwo, \
				&resampleParams);
		
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
       
	LALSDestroyVector(status->statusPtr, &(dataStreamTwo.data));
	
	/* return status */
	DETATCHSTATUSPTR( status );
	RETURN( status );
}

