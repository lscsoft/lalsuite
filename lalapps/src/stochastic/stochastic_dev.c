/*
 * stochastic_dev.c - SGWB Standalone Analysis Pipeline
 *                  - development branch
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

#include <getopt.h>

#include <FrameL.h>

#include <lal/AVFactories.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/Calibration.h>
#include <lal/ComplexFFT.h>
#include <lal/CoarseGrainFrequencySeries.h>
#include <lal/Date.h>
#include <lal/DetectorSite.h>
#include <lal/IIRFilter.h>
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

#include <lalapps.h>

#include "stochastic_dev.h"

NRCSID (STOCHASTICDEVC, "$Id$");
RCSID ("$Id$");

/* cvs info */
#define CVS_ID "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "lalapps_stochastic_dev"

/* variables for getopt option parsing */
char *optarg;
int optind;

/* flags for getopt_long */
static int apply_mask_flag;
static int high_pass_flag;
static int overlap_hann_flag;
static int debug_flag;
extern int vrbflg;

/* sampling parameters */
INT4 sampleRate = -1;
INT4 resampleRate = -1;
REAL8 deltaF = 0.25;

/* data parameters */
UINT8 startTime = 0;
UINT8 endTime = 0;
INT4 streamDuration;
INT4 intervalDuration;
INT4 segmentDuration = -1;
CHAR *frameCacheOne = NULL;
CHAR *frameCacheTwo = NULL;
CHAR *calCacheOne = NULL;
CHAR *calCacheTwo = NULL;
CHAR channelOne[LALNameLength];
CHAR channelTwo[LALNameLength];
CHAR *ifoOne = NULL;
CHAR *ifoTwo = NULL;
INT4 siteOne;
INT4 siteTwo;

/* frequency band */
INT4 fMin = -1;
INT4 fMax = -1;

/* high pass filtering parameters */
REAL4 highPassFreq = -1;
REAL4 highPassAtten = -1;
INT4 highPassOrder = -1;

/* omegaGW parameters */
REAL4 alpha = 0;
REAL4 fRef = 100;
REAL4 omega0 = 1;

/* windowing parameters */
INT4 hannDuration = -1;
INT4 maskBin = -1;

INT4 main(INT4 argc, CHAR *argv[])
{
  /* status pointer */
  LALStatus status;

  /* output file */
  FILE *out;
  CHAR outputFilename[LALNameLength];

  /* counters */
  INT4 i, j, k;

  /* results parameters */
  REAL8 y;
  REAL8 varTheo;

  /* time variables */
  LIGOTimeGPS gpsStartTime;
  LIGOTimeGPS gpsIntervalStart;
  LIGOTimeGPS gpsSegmentAStart;
  LIGOTimeGPS gpsSegmentBStart;
  LIGOTimeGPS gpsSegmentCStart;

  /* input data structures */
  INT4 padData;
  INT4 streamLength;
  ReadDataPairParams streamParams;
  StreamPair streamPair;
  REAL4TimeSeries streamOne;
  REAL4TimeSeries streamTwo;

  /* data structures for intevals */
  INT4 numIntervals;
  INT4 intervalLength;
  REAL4TimeSeries intervalOne;
  REAL4TimeSeries intervalTwo;

  /* data structures for segments */
  INT4 segmentLength;
  INT4 segmentShift;
  REAL4TimeSeries segmentOneA;
  REAL4TimeSeries segmentOneB;
  REAL4TimeSeries segmentOneC;
  REAL4TimeSeries segmentTwoA;
  REAL4TimeSeries segmentTwoB;
  REAL4TimeSeries segmentTwoC;

  /* data structures for PSDs */
  PSDEstimatorInput psdInputOne;
  PSDEstimatorInput psdInputTwo;
  PSDEstimatorParams psdEstParams;
  INT4 overlapPSDLength;
  INT4 psdTempLength;
  INT4 windowPSDLength;
  INT4 filterLength;
  INT4 numFMin;
  INT4 numFMax;
  LALWindowParams psdWindowParams;
  AverageSpectrumParams psdParams;
  REAL4FrequencySeries psdTempOne;
  REAL4FrequencySeries psdTempTwo;
  REAL4FrequencySeries psdOne;
  REAL4FrequencySeries psdTwo;

  /* response functions */
  COMPLEX8FrequencySeries responseTempOneA;
  COMPLEX8FrequencySeries responseTempOneB;
  COMPLEX8FrequencySeries responseTempOneC;
  COMPLEX8FrequencySeries responseTempTwoA;
  COMPLEX8FrequencySeries responseTempTwoB;
  COMPLEX8FrequencySeries responseTempTwoC;
  COMPLEX8FrequencySeries responseOneA;
  COMPLEX8FrequencySeries responseOneB;
  COMPLEX8FrequencySeries responseOneC;
  COMPLEX8FrequencySeries responseTwoA;
  COMPLEX8FrequencySeries responseTwoB;
  COMPLEX8FrequencySeries responseTwoC;
  INT4 respLength;
  CalibrationUpdateParams calfacts;
  LALUnit countPerAttoStrain = {18,{0,0,0,0,0,-1,1},{0,0,0,0,0,0,0}};

  /* inverse noise data structures */
  REAL4FrequencySeries calInvPSDOne;
  REAL4FrequencySeries calInvPSDTwo;
  LALUnit calPSDUnit = {36,{0,0,-1,0,0,-2,0},{0,0,0,0,0,0,0}};

  /* window for segment data streams */
  REAL4TimeSeries dataWindow;

  /* hann window */
  INT4 hannLength;
  LALWindowParams hannParams;
  REAL4Vector *hannWindow;

  /* high pass filtering */
  PassBandParamStruc highPassParam;

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
  REAL4FrequencySeries maskTemp;
  INT4 nBin;

  /* spectrum structures */
  StochasticOmegaGWParameters omegaGWParams;
  REAL4FrequencySeries omegaGW;

  /* structures for optimal filter normalisation */
  StochasticOptimalFilterNormalizationInput normInput;
  StochasticOptimalFilterNormalizationOutput normOutput;
  StochasticOptimalFilterNormalizationParameters normParams;
  REAL4WithUnits normLambda;
  REAL4WithUnits normSigma;
  REAL8 lambda;

  /* structures for optimal filter */
  REAL4FrequencySeries optFilter;
  StochasticOptimalFilterCalInput optFilterIn;

  /* structures for CC spectrum and CC statistics */
  StochasticCrossCorrelationCalInput ccIn;
  BOOLEAN epochsMatch = 1;
  REAL4WithUnits ccStat;
  COMPLEX8FrequencySeries ccSpectrum;

  /* set default error behaviour */
  status.statusPtr = NULL;
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level("33");

  /* parse command line options */
  parseOptions(argc, argv);

  /* should a resample buffer be applied */
  if ((sampleRate == resampleRate) || (high_pass_flag == 0))
  {
    padData = 0;
  }
  else
  {
    padData = 1;
  }

  /* set start time */
  gpsStartTime.gpsSeconds = startTime;
  gpsStartTime.gpsNanoSeconds = 0;

  if (vrbflg)
  {
    fprintf(stdout, "Calculating number of intervals...\n");
  }

  /* get number of intervals */
  streamDuration = endTime - startTime;
  if (overlap_hann_flag)
  {
    numIntervals = (2 * (streamDuration / intervalDuration)) - 1;
    segmentShift = segmentDuration / 2;
  }
  else
  {
    numIntervals = streamDuration / intervalDuration;
    segmentShift = segmentDuration;
  }

  /* get stream length */
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

  if (vrbflg)
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
  streamPair.streamOne = &streamOne;
  streamPair.streamTwo = &streamTwo;

  /* set length for data intervals */
  intervalLength = intervalDuration * resampleRate;

  /* set metadata fields for data intervals */
  strncpy(intervalOne.name, "intervalOne", LALNameLength);
  strncpy(intervalTwo.name, "intervalTwo", LALNameLength);
  intervalOne.sampleUnits = streamOne.sampleUnits;
  intervalTwo.sampleUnits = streamOne.sampleUnits;
  intervalOne.deltaT = 1./(REAL8)resampleRate;
  intervalTwo.deltaT = 1./(REAL8)resampleRate;
  intervalOne.f0 = 0;
  intervalTwo.f0 = 0;

  if (vrbflg)
  {
    fprintf(stdout, "Allocating memory for data intervals...\n");
  }

  /* allocate memory for data intervals */
  intervalOne.data = (REAL4Sequence*)LALCalloc(1, sizeof(REAL4Sequence));
  intervalTwo.data = (REAL4Sequence*)LALCalloc(1, sizeof(REAL4Sequence));
  intervalOne.data->length = intervalLength;
  intervalTwo.data->length = intervalLength;

  /* set length for data segments */
  segmentLength = segmentDuration * resampleRate;

  /* set metadata fields for data segments */
  strncpy(segmentOneA.name, "segmentOneA", LALNameLength);
  strncpy(segmentOneB.name, "segmentOneB", LALNameLength);
  strncpy(segmentOneC.name, "segmentOneC", LALNameLength);
  strncpy(segmentTwoA.name, "segmentTwoA", LALNameLength);
  strncpy(segmentTwoB.name, "segmentTwoB", LALNameLength);
  strncpy(segmentTwoC.name, "segmentTwoC", LALNameLength);
  segmentOneA.sampleUnits = streamOne.sampleUnits;
  segmentOneB.sampleUnits = streamOne.sampleUnits;
  segmentOneC.sampleUnits = streamOne.sampleUnits;
  segmentTwoA.sampleUnits = streamTwo.sampleUnits;
  segmentTwoB.sampleUnits = streamTwo.sampleUnits;
  segmentTwoC.sampleUnits = streamTwo.sampleUnits;
  segmentOneA.deltaT = 1./(REAL8)resampleRate;
  segmentOneB.deltaT = 1./(REAL8)resampleRate;
  segmentOneC.deltaT = 1./(REAL8)resampleRate;
  segmentTwoA.deltaT = 1./(REAL8)resampleRate;
  segmentTwoB.deltaT = 1./(REAL8)resampleRate;
  segmentTwoC.deltaT = 1./(REAL8)resampleRate;
  segmentOneA.f0 = 0;
  segmentOneB.f0 = 0;
  segmentOneC.f0 = 0;
  segmentTwoA.f0 = 0;
  segmentTwoB.f0 = 0;
  segmentTwoC.f0 = 0;

  if (vrbflg)
  {
    fprintf(stdout, "Allocating memory for data segments...\n");
  }

  /* allocate memory for data segments */
  segmentOneA.data = (REAL4Sequence*)LALCalloc(1, sizeof(REAL4Sequence));
  segmentOneB.data = (REAL4Sequence*)LALCalloc(1, sizeof(REAL4Sequence));
  segmentOneC.data = (REAL4Sequence*)LALCalloc(1, sizeof(REAL4Sequence));
  segmentTwoA.data = (REAL4Sequence*)LALCalloc(1, sizeof(REAL4Sequence));
  segmentTwoB.data = (REAL4Sequence*)LALCalloc(1, sizeof(REAL4Sequence));
  segmentTwoC.data = (REAL4Sequence*)LALCalloc(1, sizeof(REAL4Sequence));
  segmentOneA.data->length = segmentLength;
  segmentOneB.data->length = segmentLength;
  segmentOneC.data->length = segmentLength;
  segmentTwoA.data->length = segmentLength;
  segmentTwoB.data->length = segmentLength;
  segmentTwoC.data->length = segmentLength;

  /* set PSD window length */
  windowPSDLength = (UINT4)(resampleRate / deltaF);

  /* set parameters for PSD estimation */
  overlapPSDLength = windowPSDLength / 2;
  psdTempLength = (windowPSDLength / 2) + 1;
  numFMin = (UINT4)(fMin / deltaF);
  numFMax = (UINT4)(fMax / deltaF);
  filterLength = numFMax - numFMin + 1;

  if (vrbflg)
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

  if (vrbflg)
  {
    fprintf(stdout, "Allocating memory for reduced frequency band PSDs...\n");
  }

  /* allocate memory for reduced frequency band PSDs */
  psdOne.data = (REAL4Sequence*)LALCalloc(1, sizeof(REAL4Sequence));
  psdTwo.data = (REAL4Sequence*)LALCalloc(1, sizeof(REAL4Sequence));
  psdOne.data->length = filterLength;
  psdTwo.data->length = filterLength;

  /* set window parameters for PSD estimation */
  psdWindowParams.length = windowPSDLength;
  psdWindowParams.type = Hann;

  /* set parameters for PSD estimation */
  psdParams.method = useMean;
  psdParams.overlap = overlapPSDLength;
  psdParams.plan = NULL;
  psdParams.window = NULL;

  /* set psd inputs */
  psdInputOne.segmentA = &segmentOneA;
  psdInputOne.segmentC = &segmentOneC;
  psdInputOne.responseA = &responseOneA;
  psdInputOne.responseC = &responseOneC;
  psdInputTwo.segmentA = &segmentTwoA;
  psdInputTwo.segmentC = &segmentTwoC;
  psdInputTwo.responseA = &responseTwoA;
  psdInputTwo.responseC = &responseTwoC;
  psdEstParams.psdTempLength = psdTempLength;
  psdEstParams.filterLength = filterLength;
  psdEstParams.psdParams = &psdParams;
  psdEstParams.numFMin = numFMin;

  if (vrbflg)
  {
    fprintf(stdout, "Creating FFT plan for PSD estimation...\n");
  }

  /* create fft plan */
  LAL_CALL (LALCreateForwardRealFFTPlan(&status, &psdParams.plan, \
        windowPSDLength, 0), &status );

  if (vrbflg)
  {
    fprintf(stdout, "Creating window for PSD estimation...\n");
  }

  /* create window for PSD estimation */
  LAL_CALL( LALCreateREAL4Window(&status, &psdParams.window, \
        &psdWindowParams), &status );

  /* set parameters for response functions */
  respLength = (UINT4)(fMax / deltaF) + 1;

  if (vrbflg)
  {
    fprintf(stdout, "Allocating memory for response functions...\n");
  }

  /* allocate memory for response functions */
  responseTempOneA.data = NULL;
  responseTempOneB.data = NULL;
  responseTempOneC.data = NULL;
  responseTempTwoA.data = NULL;
  responseTempTwoB.data = NULL;
  responseTempTwoC.data = NULL;
  LAL_CALL( LALCCreateVector(&status, &(responseTempOneA.data), respLength), \
      &status );
  LAL_CALL( LALCCreateVector(&status, &(responseTempOneB.data), respLength), \
      &status );
  LAL_CALL( LALCCreateVector(&status, &(responseTempOneC.data), respLength), \
      &status );
  LAL_CALL( LALCCreateVector(&status, &(responseTempTwoA.data), respLength), \
      &status );
  LAL_CALL( LALCCreateVector(&status, &(responseTempTwoB.data), respLength), \
      &status );
  LAL_CALL( LALCCreateVector(&status, &(responseTempTwoC.data), respLength), \
      &status );
  memset(responseTempOneA.data->data, 0, \
      responseTempOneA.data->length * sizeof(*responseTempOneA.data->data));
  memset(responseTempOneB.data->data, 0, \
      responseTempOneB.data->length * sizeof(*responseTempOneB.data->data));
  memset(responseTempOneC.data->data, 0, \
      responseTempOneC.data->length * sizeof(*responseTempOneC.data->data));
  memset(responseTempTwoA.data->data, 0, \
      responseTempTwoA.data->length * sizeof(*responseTempTwoA.data->data));
  memset(responseTempTwoB.data->data, 0, \
      responseTempTwoB.data->length * sizeof(*responseTempTwoB.data->data));
  memset(responseTempTwoC.data->data, 0, \
      responseTempTwoC.data->length * sizeof(*responseTempTwoC.data->data));

  /* reduced frequency band response functions */
  /* set metadata fields for reduced frequency band response functions */
  strncpy(responseOneA.name, "responseOneA", LALNameLength);
  strncpy(responseOneB.name, "responseOneB", LALNameLength);
  strncpy(responseOneC.name, "responseOneC", LALNameLength);
  strncpy(responseTwoA.name, "responseTwoA", LALNameLength);
  strncpy(responseTwoB.name, "responseTwoB", LALNameLength);
  strncpy(responseTwoC.name, "responseTwoC", LALNameLength);
  responseOneA.sampleUnits = countPerAttoStrain;
  responseOneB.sampleUnits = countPerAttoStrain;
  responseOneC.sampleUnits = countPerAttoStrain;
  responseTwoA.sampleUnits = countPerAttoStrain;
  responseTwoB.sampleUnits = countPerAttoStrain;
  responseTwoC.sampleUnits = countPerAttoStrain;
  responseOneA.deltaF = deltaF;
  responseOneB.deltaF = deltaF;
  responseOneC.deltaF = deltaF;
  responseTwoA.deltaF = deltaF;
  responseTwoB.deltaF = deltaF;
  responseTwoC.deltaF = deltaF;
  responseOneA.f0 = fMin;
  responseOneB.f0 = fMin;
  responseOneC.f0 = fMin;
  responseTwoA.f0 = fMin;
  responseTwoB.f0 = fMin;
  responseTwoC.f0 = fMin;

  if (vrbflg)
  {
    fprintf(stdout, "Allocating memory for reduced frequency band response " \
        "functions...\n");
  }

  /* allocate memory for reduced frequency band response functions */
  responseOneA.data = (COMPLEX8Sequence*)LALCalloc(1, sizeof(COMPLEX8Sequence));
  responseOneB.data = (COMPLEX8Sequence*)LALCalloc(1, sizeof(COMPLEX8Sequence));
  responseOneC.data = (COMPLEX8Sequence*)LALCalloc(1, sizeof(COMPLEX8Sequence));
  responseTwoA.data = (COMPLEX8Sequence*)LALCalloc(1, sizeof(COMPLEX8Sequence));
  responseTwoB.data = (COMPLEX8Sequence*)LALCalloc(1, sizeof(COMPLEX8Sequence));
  responseTwoC.data = (COMPLEX8Sequence*)LALCalloc(1, sizeof(COMPLEX8Sequence));
  responseOneA.data->length = filterLength;
  responseOneB.data->length = filterLength;
  responseOneC.data->length = filterLength;
  responseTwoA.data->length = filterLength;
  responseTwoB.data->length = filterLength;
  responseTwoC.data->length = filterLength;

  /* set metadata fields for inverse noise structures */
  strncpy(calInvPSDOne.name, "calInvPSDOne", LALNameLength);
  strncpy(calInvPSDTwo.name, "calInvPSDTwo", LALNameLength);
  calInvPSDOne.sampleUnits = calPSDUnit;
  calInvPSDTwo.sampleUnits = calPSDUnit;
  calInvPSDOne.deltaF = deltaF;
  calInvPSDTwo.deltaF = deltaF;
  calInvPSDOne.f0 = fMin;
  calInvPSDTwo.f0 = fMin;

  if (vrbflg)
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

  /* set window parameters for segment data streams */
  strncpy(dataWindow.name, "dataWindow", LALNameLength);
  dataWindow.sampleUnits = lalDimensionlessUnit;
  dataWindow.deltaT = 1./resampleRate;
  dataWindow.f0 = 0;

  if (vrbflg)
  {
    fprintf(stdout, "Allocating memory for data segment window...\n");
  }

  /* allocate memory for segment window */
  dataWindow.data = NULL;
  LAL_CALL( LALSCreateVector(&status, &(dataWindow.data), segmentLength), \
      &status );
  memset(dataWindow.data->data, 0, \
      dataWindow.data->length * sizeof(*dataWindow.data->data));

  if (vrbflg)
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
  if (debug_flag)
  {
    LALSPrintTimeSeries(&dataWindow, "dataWindow.dat");
  }

  /* structure for high pass filtering */
  if (high_pass_flag)
  {
    highPassParam.nMax = highPassOrder;
    highPassParam.f2 = highPassFreq;
    highPassParam.a2 = highPassAtten;
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

  if (vrbflg)
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

  if (vrbflg)
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

  if (vrbflg)
  {
    fprintf(stdout, "Generating the overlap reduction function...\n");
  }

  /* generate overlap reduction function */
  LAL_CALL (LALOverlapReductionFunction(&status, &overlap, &detectors, \
        &ORFparams), &status );

  /* save */
  if (debug_flag)
  {
    LALSPrintFrequencySeries(&overlap, "overlap.dat");
  }

  /* set metadata fields for spectrum */
  strncpy(omegaGW.name, "omegaGW", LALNameLength);
  omegaGW.sampleUnits = lalDimensionlessUnit;
  omegaGW.deltaF = deltaF;
  omegaGW.f0 = fMin;

  if (vrbflg)
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
  omegaGWParams.omegaRef = omega0;
  omegaGWParams.length = filterLength;
  omegaGWParams.f0 = fMin;
  omegaGWParams.deltaF = deltaF;

  if (vrbflg)
  {
    fprintf(stdout, "Generating spectrum for optimal filter...\n");
  }

  /* generage omegaGW */
  LAL_CALL( LALStochasticOmegaGW(&status, &omegaGW, &omegaGWParams), &status );

  /* frequency mask */
  if (apply_mask_flag)
  {
    /* extra bins */
    nBin = (maskBin - 1) / 2;

    /* set metadata fields for frequency mask */
    strncpy(mask.name, "mask", LALNameLength);
    mask.deltaF = deltaF;
    mask.f0 = fMin;
    mask.sampleUnits = lalDimensionlessUnit;

    if (vrbflg)
    {
      fprintf(stdout, "Allocating memory for frequency mask...\n");
    }

    /* allocate memory for frequency mask */
    maskTemp.data = NULL;
    LAL_CALL( LALCreateVector(&status, &(maskTemp.data), respLength), \
        &status );
    memset(maskTemp.data->data, 0, \
        maskTemp.data->length * sizeof(*maskTemp.data->data));

    /* reduced band frequency mask */
    mask.data = (REAL4Sequence*)LALCalloc(1, sizeof(REAL4Sequence));
    mask.data->length = filterLength;

    if (vrbflg)
    {
      fprintf(stdout, "Generating frequency mask...\n");
    }

    /* set all values to 1 */
    for (i = 0; i < respLength; i++)
    {
      maskTemp.data->data[i] = 1.;
    }

    if (vrbflg)
    {
      fprintf(stdout, "Masking multiples of 16 Hz...\n");
    }

    /* remove multiples of 16 Hz */
    for (i = 0; i < respLength; i += (UINT4)(16 / deltaF))
    {
      maskTemp.data->data[i] = 0.;

      for (k = 0; k < nBin; k++)
      {
        if ((i + 1 + k) < respLength)
        {
          maskTemp.data->data[i + 1 + k] = 0.;
        }
        if ((i - 1 - k) > 0)
        {
          maskTemp.data->data[i - 1 - k] = 0.;
        }
      }
    }

    if (vrbflg)
    {
      fprintf(stdout, "Masking multiples of 60 Hz...\n");
    }

    /* remove multiples of 60 Hz */
    for (i = 0; i < respLength; i += (UINT4)(60 / deltaF))
    {
      maskTemp.data->data[i] = 0.;

      for (k = 0; k < nBin; k++)
      {
        if ((i + 1 + k) < respLength)
        {
          maskTemp.data->data[i + 1 + k] = 0.;
        }
        if ((i - 1 - k) > 0)
        {
          maskTemp.data->data[i - 1 - k] = 0.;
        }
      }
    }

    if (vrbflg)
    {
      fprintf(stdout, "Getting appropriate frequency band for mask...\n");
    }

    /* get appropriate band */
    mask.data->data = maskTemp.data->data + numFMin;

    if (debug_flag)
    {
      LALSPrintFrequencySeries(&mask, "mask.dat");
    }

    if (vrbflg)
    {
      fprintf(stdout, "Applying frequency mask to spectrum...\n");
    }

    /* apply mask to omegaGW */
    for (i = 0; i < filterLength; i++)
    {
      omegaGW.data->data[i] *= mask.data->data[i];
    }
  }

  /* save */
  if (debug_flag)
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
  optFilter.deltaF = deltaF;
  optFilter.f0 = fMin;

  if (vrbflg)
  {
    fprintf(stdout, "Allocating memory for optimal filter...\n");
  }

  /* allocate memory for optimal filter */
  optFilter.data = NULL;
  LAL_CALL( LALCreateVector(&status, &(optFilter.data), filterLength), \
      &status );
  memset(optFilter.data->data, 0, \
      optFilter.data->length * sizeof(*optFilter.data->data));

  /* set optimal filter inputs */
  optFilterIn.overlapReductionFunction = &overlap;
  optFilterIn.omegaGW = &omegaGW;
  optFilterIn.calibratedInverseNoisePSD1 = &calInvPSDOne;
  optFilterIn.calibratedInverseNoisePSD2 = &calInvPSDTwo;

  /* set metadata fields for CC spectrum */
  strncpy(ccSpectrum.name, "ccSpectrum", LALNameLength);
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
  ccIn.responseFunctionOne = &responseOneB;
  ccIn.responseFunctionTwo = &responseTwoB;
  ccIn.optimalFilter = &optFilter;

  if (vrbflg)
  {
    fprintf(stdout, "Reading data...\n");
  }

  /* read data */
  readDataPair(&status, &streamPair, &streamParams);

  /* save */
  if (debug_flag)
  {
    LALSPrintTimeSeries(&streamOne, "stream1.dat");
    LALSPrintTimeSeries(&streamTwo, "stream2.dat");
  }

  if (vrbflg)
  {
    fprintf(stdout, "High pass filtering...\n");
  }

  /* high pass filter */
  if (high_pass_flag)
  {
    LAL_CALL( LALButterworthREAL4TimeSeries( &status, &streamOne, \
          &highPassParam ), &status );
    LAL_CALL( LALButterworthREAL4TimeSeries( &status, &streamTwo, \
          &highPassParam ), &status );
  }

  /* open output file */
  LALSnprintf(outputFilename, LALNameLength, "stoch-%s%s-%d-%d.dat", 
      ifoOne, ifoTwo, (INT4)startTime, (INT4)endTime);
  if ((out = fopen(outputFilename, "w")) == NULL)
  {
    fprintf(stderr, "Can't open file for output...\n");
    exit(1);
  }

  if (vrbflg)
  {
    fprintf(stdout, "Looping over %d intervals...\n", numIntervals);
  }

  /* loop over intervals */
  for (j = 0; j < numIntervals; j++)
  {
    /* define interval epoch */
    gpsIntervalStart.gpsSeconds = startTime + (j * segmentShift);
    gpsIntervalStart.gpsNanoSeconds = 0;
    intervalOne.epoch = gpsIntervalStart;
    intervalTwo.epoch = gpsIntervalStart;

    /* define segment epoch */
    gpsSegmentAStart.gpsSeconds = startTime + (j * segmentShift);
    gpsSegmentAStart.gpsNanoSeconds = 0;
    gpsSegmentBStart.gpsSeconds = startTime + (j * segmentShift) + \
      segmentDuration;
    gpsSegmentBStart.gpsNanoSeconds = 0;
    gpsSegmentCStart.gpsSeconds = startTime + (j * segmentShift) + \
      (2 * segmentDuration);
    gpsSegmentCStart.gpsNanoSeconds = 0;
    segmentOneA.epoch = gpsSegmentAStart;
    segmentOneB.epoch = gpsSegmentBStart;
    segmentOneC.epoch = gpsSegmentCStart;
    segmentTwoA.epoch = gpsSegmentAStart;
    segmentTwoB.epoch = gpsSegmentBStart;
    segmentTwoC.epoch = gpsSegmentCStart;

    if (vrbflg)
    {
      fprintf(stdout, "Performing search on interval %d of %d...\n", \
          j + 1, numIntervals);
    }

    /* build interval */
    intervalOne.data->data = streamOne.data->data + (j * \
        segmentShift * resampleRate);
    intervalTwo.data->data = streamTwo.data->data + (j * \
        segmentShift * resampleRate);

    /* build segments */
    segmentOneA.data->data = streamOne.data->data + \
      (resampleRate * j * segmentShift);
    segmentOneB.data->data = streamOne.data->data + \
      (resampleRate * ((j * segmentShift) + segmentDuration));
    segmentOneC.data->data = streamOne.data->data + \
      (resampleRate * ((j * segmentShift) + (2 * segmentDuration)));
    segmentTwoA.data->data = streamTwo.data->data + \
      (resampleRate * j * segmentShift);
    segmentTwoB.data->data = streamTwo.data->data + \
      (resampleRate * ((j * segmentShift) + segmentDuration));
    segmentTwoC.data->data = streamTwo.data->data + \
      (resampleRate * ((j * segmentShift) + (2 * segmentDuration)));

    /* output the results */
    if (debug_flag)
    {
      LALSPrintTimeSeries(&intervalOne, "interval1.dat");
      LALSPrintTimeSeries(&intervalTwo, "interval2.dat");
      LALSPrintTimeSeries(&segmentOneA, "segment1a.dat");
      LALSPrintTimeSeries(&segmentOneB, "segment1b.dat");
      LALSPrintTimeSeries(&segmentOneC, "segment1c.dat");
      LALSPrintTimeSeries(&segmentTwoA, "segment2a.dat");
      LALSPrintTimeSeries(&segmentTwoB, "segment2b.dat");
      LALSPrintTimeSeries(&segmentTwoC, "segment2c.dat");
    }

    if (vrbflg)
    {
      fprintf(stdout, "Zero padding data...\n");
    }

    /* zero pad and fft */
    LAL_CALL( LALSZeroPadAndFFT(&status, &hBarTildeOne, &segmentOneB, \
          &zeroPadParams), &status );
    LAL_CALL( LALSZeroPadAndFFT(&status, &hBarTildeTwo, &segmentTwoB, \
          &zeroPadParams), &status );

    /* save */
    if (debug_flag)
    {
      LALCPrintFrequencySeries(&hBarTildeOne, "hBarTilde1.dat");
      LALCPrintFrequencySeries(&hBarTildeTwo, "hBarTilde2.dat");
    }

    if (vrbflg)
    {
      fprintf(stdout, "Generating response functions...\n");
    }

    /* compute response function */
    responseTempOneA.epoch = gpsSegmentAStart;
    responseTempOneB.epoch = gpsSegmentBStart;
    responseTempOneC.epoch = gpsSegmentCStart;
    responseTempTwoA.epoch = gpsSegmentAStart;
    responseTempTwoB.epoch = gpsSegmentBStart;
    responseTempTwoC.epoch = gpsSegmentCStart;
    memset(&calfacts, 0, sizeof(CalibrationUpdateParams));
    calfacts.ifo = ifoOne;
    LAL_CALL( LALExtractFrameResponse(&status, &responseTempOneA, \
          calCacheOne, &calfacts), &status );
    memset(&calfacts, 0, sizeof(CalibrationUpdateParams));
    calfacts.ifo = ifoOne;
    LAL_CALL( LALExtractFrameResponse(&status, &responseTempOneB, \
          calCacheOne, &calfacts), &status );
    memset(&calfacts, 0, sizeof(CalibrationUpdateParams));
    calfacts.ifo = ifoOne;
    LAL_CALL( LALExtractFrameResponse(&status, &responseTempOneC, \
          calCacheOne, &calfacts), &status );
    memset(&calfacts, 0, sizeof(CalibrationUpdateParams));
    calfacts.ifo = ifoTwo;
    LAL_CALL( LALExtractFrameResponse(&status, &responseTempTwoA, \
          calCacheTwo, &calfacts), &status );
    memset(&calfacts, 0, sizeof(CalibrationUpdateParams));
    calfacts.ifo = ifoTwo;
    LAL_CALL( LALExtractFrameResponse(&status, &responseTempTwoB, \
          calCacheTwo, &calfacts), &status );
    memset(&calfacts, 0, sizeof(CalibrationUpdateParams));
    calfacts.ifo = ifoTwo;
    LAL_CALL( LALExtractFrameResponse(&status, &responseTempTwoC, \
          calCacheTwo, &calfacts), &status );

    if (vrbflg)
    {
      fprintf(stdout, "Getting appropriate frequency band for response " \
          "functions...\n");
    }

    /* reduce to the optimal filter frequency range */
    responseOneA.data->data = responseTempOneA.data->data + numFMin;
    responseOneB.data->data = responseTempOneB.data->data + numFMin;
    responseOneC.data->data = responseTempOneC.data->data + numFMin;
    responseTwoA.data->data = responseTempTwoA.data->data + numFMin;
    responseTwoB.data->data = responseTempTwoB.data->data + numFMin;
    responseTwoC.data->data = responseTempTwoC.data->data + numFMin;

    /* output the results */
    if (debug_flag)
    {
      LALCPrintFrequencySeries(&responseOneA, "response1a.dat");
      LALCPrintFrequencySeries(&responseOneB, "response1b.dat");
      LALCPrintFrequencySeries(&responseOneC, "response1c.dat");
      LALCPrintFrequencySeries(&responseTwoA, "response2a.dat");
      LALCPrintFrequencySeries(&responseTwoB, "response2b.dat");
      LALCPrintFrequencySeries(&responseTwoC, "response2c.dat");
    }

    if (vrbflg)
    {
      fprintf(stdout, "Estimating PSDs...\n");
    }

    /* set epoch */
    calInvPSDOne.epoch = segmentOneB.epoch;
    calInvPSDTwo.epoch = segmentOneB.epoch;

    /* estimate psds */
    LAL_CALL( psdEstimator(&status, calInvPSDOne, psdInputOne, \
          psdEstParams), &status );
    LAL_CALL( psdEstimator(&status, calInvPSDTwo, psdInputTwo, \
          psdEstParams), &status );

    /* output the results */
    if (debug_flag)
    {
      LALSPrintFrequencySeries(&calInvPSDOne, "inPSD1.dat");
      LALSPrintFrequencySeries(&calInvPSDTwo, "inPSD2.dat");
    }

    if (vrbflg)
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
    if (vrbflg)
    {
      fprintf(stdout, "lambda = %e s^-1\n", lambda);
      fprintf(stdout, "varTheo = %e s\n", varTheo);
    }

    if (vrbflg)
    {
      fprintf(stdout, "Generating optimal filter...\n");
    }

    /* build optimal filter */
    optFilter.epoch = gpsSegmentAStart;
    LAL_CALL( LALStochasticOptimalFilterCal(&status, &optFilter, \
          &optFilterIn, &normLambda), &status );

    /* save */
    if (debug_flag)
    {
      LALPrintFrequencySeries(&optFilter, "optFilter.dat");
    }

    /* save */
    if (debug_flag)
    {
      if (vrbflg)
      {
        fprintf(stdout, "Generating cross correlation spectrum...\n");
      }

      /* cc spectrum */
      LAL_CALL( LALStochasticCrossCorrelationSpectrumCal(&status, \
            &ccSpectrum, &ccIn, epochsMatch), &status );

      LALCPrintFrequencySeries(&ccSpectrum, "ccSpectrum.dat");
    }

    if (vrbflg)
    {
      fprintf(stdout, "Generating cross correlation statistic...\n");
    }

    /* cc statistic */
    LAL_CALL( LALStochasticCrossCorrelationStatisticCal(&status, &ccStat, \
          &ccIn, epochsMatch), &status );

    y = (REAL8)(ccStat.value * pow(10.,ccStat.units.powerOfTen));

    if (vrbflg)
    {
      fprintf(stdout, "y = %f\n\n", y);
    }

    /* output to file */
    fprintf(out, "%d %e %e\n", gpsSegmentBStart.gpsSeconds, sqrt(varTheo), y);
  }

  /* close output file */
  fclose(out);

  /* cleanup */
  LAL_CALL( LALDestroyRealFFTPlan(&status, &(psdParams.plan)), &status );
  LAL_CALL( LALDestroyRealFFTPlan(&status, &fftDataPlan), &status );
  LALFree(intervalOne.data);
  LALFree(intervalTwo.data);
  LALFree(segmentOneA.data);
  LALFree(segmentOneB.data);
  LALFree(segmentOneC.data);
  LALFree(segmentTwoA.data);
  LALFree(segmentTwoB.data);
  LALFree(segmentTwoC.data);
  LAL_CALL( LALDestroyVector(&status, &(streamOne.data)), &status );
  LAL_CALL( LALDestroyVector(&status, &(streamTwo.data)), &status );
  LALFree(psdOne.data);
  LALFree(psdTwo.data);
  LAL_CALL( LALDestroyVector(&status, &(psdTempOne.data)), &status );
  LAL_CALL( LALDestroyVector(&status, &(psdTempTwo.data)), &status );
  LALFree(responseOneA.data);
  LALFree(responseOneB.data);
  LALFree(responseOneC.data);
  LALFree(responseTwoA.data);
  LALFree(responseTwoB.data);
  LALFree(responseTwoC.data);
  LAL_CALL( LALCDestroyVector(&status, &(responseTempOneA.data)), &status );
  LAL_CALL( LALCDestroyVector(&status, &(responseTempOneB.data)), &status );
  LAL_CALL( LALCDestroyVector(&status, &(responseTempOneC.data)), &status );
  LAL_CALL( LALCDestroyVector(&status, &(responseTempTwoA.data)), &status );
  LAL_CALL( LALCDestroyVector(&status, &(responseTempTwoB.data)), &status );
  LAL_CALL( LALCDestroyVector(&status, &(responseTempTwoC.data)), &status );
  LAL_CALL( LALDestroyVector(&status, &(calInvPSDOne.data)), &status );
  LAL_CALL( LALDestroyVector(&status, &(calInvPSDTwo.data)), &status );
  LAL_CALL( LALDestroyVector(&status, &(overlap.data)), &status );
  LAL_CALL( LALDestroyVector(&status, &(omegaGW.data)), &status );
  LAL_CALL( LALDestroyVector(&status, &(dataWindow.data)), &status );
  if (apply_mask_flag)
  {
    LALFree(mask.data);
    LAL_CALL( LALDestroyVector(&status, &(maskTemp.data)), &status );
  }
  if (hannDuration != 0)
  {
    LAL_CALL( LALDestroyVector(&status, &hannWindow), &status );
  }
  LAL_CALL( LALCDestroyVector(&status, &(hBarTildeOne.data)), &status );
  LAL_CALL( LALCDestroyVector(&status, &(hBarTildeTwo.data)), &status );

  /* free calloc'd memory */
  free(frameCacheOne);
  free(frameCacheTwo);
  free(calCacheOne);
  free(calCacheTwo);
  free(ifoOne);
  free(ifoTwo);

  return 0;
}

#define USAGE \
  "Usage: pipeline [options]\n"\
  " --help                        print this message\n"\
  " --version                     display version\n"\
  " --verbose                     verbose mode\n"\
  " --debug                       save out intermediate products\n"\
  " --debug-level N               set lalDebugLevel\n"\
  " --gps-start-time N            GPS start time\n"\
  " --gps-end-time N              GPS end time\n"\
  " --segment-duration N          segment duration\n"\
  " --sample-rate N               sample rate\n"\
  " --resample-rate N             resample rate\n"\
  " --f-min N                     minimal frequency\n"\
  " --f-max N                     maximal frequency\n"\
  " --ifo-one IFO                 ifo for first stream\n"\
  " --ifo-two IFO                 ifo for second stream\n"\
  " --frame-cache-one FILE        cache file for first stream\n"\
  " --frame-cache-two FILE        cache file for second stream\n"\
  " --calibration-cache-one FILE  first stream calibration cache\n"\
  " --calibration-cache-two FILE  second stream calibration cache\n"\
  " --apply-mask                  apply frequency masking\n"\
  " --mask-bin N                  number of bins to mask\n"\
  " --overlap-hann                overlaping hann windows\n"\
  " --hann-duration N             hann duration\n"\
  " --high-pass-filter            apply high pass filtering\n"\
  " --hpf-frequency N             high pass filter knee frequency\n"\
  " --hpf-attenuation N           high pass filter attenuation\n"\
  " --hpf-order N                 high pass filter order\n"\
  " --filter-omega-alpha N        omega_gw exponent\n"\
  " --filter-omega-fref N         omega_gw reference frequency\n"\
  " --filter-omega0 N             omega_0\n"

/* parse command line options */
static void parseOptions(INT4 argc, CHAR *argv[])
{
  int c = -1;

  while(1)
  {
    static struct option long_options[] =
    {
      /* options that set a flag */
      {"apply-mask", no_argument, &apply_mask_flag, 1},
      {"high-pass-filter", no_argument, &high_pass_flag, 1},
      {"overlap-hann", no_argument, &overlap_hann_flag, 1},
      {"verbose", no_argument, &vrbflg, 1},
      {"debug", no_argument, &debug_flag, 1},
      /* options that don't set a flag */
      {"mask-bin", required_argument, 0, 'b'},
      {"help", no_argument, 0, 'h'},
      {"gps-start-time", required_argument, 0, 't'},
      {"gps-end-time", required_argument, 0, 'T'},
      {"segment-duration", required_argument, 0, 'l'},
      {"sample-rate", required_argument, 0, 'a'},
      {"resample-rate", required_argument, 0, 'A'},
      {"f-min", required_argument, 0, 'f'},
      {"f-max", required_argument, 0, 'F'},
      {"hann-duration", required_argument, 0, 'w'},
      {"ifo-one", required_argument, 0, 'i'},
      {"ifo-two", required_argument, 0, 'I'},
      {"frame-cache-one", required_argument, 0, 'd'},
      {"frame-cache-two", required_argument, 0, 'D'},
      {"calibration-cache-one", required_argument, 0, 'r'},
      {"calibration-cache-two", required_argument, 0, 'R'},
      {"debug-level", required_argument, 0, 'z'},
      {"hpf-frequency", required_argument, 0, 'k'},
      {"hpf-attenuation", required_argument, 0, 'p'},
      {"hpf-order", required_argument, 0, 'P'},
      {"filter-omega-alpha", required_argument, 0, 'e'},
      {"filter-omega-fref", required_argument, 0, 'E'},
      {"filter-omega0", required_argument, 0, 'O'},
      {"version", no_argument, 0, 'V'},
      {0, 0, 0, 0}
    };

    /* getopt_long stores the option here */
    int option_index = 0;
    size_t optarg_len;

    c = getopt_long(argc, argv, "ht:T:l:a:f:F:w:i:I:d:D:r:R:z:k:p:P:V", \
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
        fprintf(stdout, USAGE);
        exit(0);
        break;

      case 'b':
        /* number of bins to mask */
        maskBin = atoi(optarg);

        /* check */
        if (maskBin <= 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Number of bins to mask must be greater than 0: " \
              "(%d specified)\n", long_options[option_index].name, maskBin);
          exit(1);
        }

        break;

      case 't':
        /* start time */
        startTime = atoi(optarg);

        /* check */
        if (startTime < 441217609)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "GPS start time is prior to 1 January 1994 00:00:00 UTC " \
              "(%lld specified)\n", long_options[option_index].name, \
              startTime);
          exit(1);
        }
        if (startTime > 999999999)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "GPS start time is after 14 September 2011 01:46:26 UTC " \
              "(%lld specified)\n", long_options[option_index].name, \
              startTime);
          exit(1);
        }

        break;

      case 'T':
        /* stop time */
        endTime = atoi(optarg);

        /* check */
        if (endTime < 441217609)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "GPS end time is prior to 1 January 1994 00:00:00 UTC " \
              "(%lld specified)\n", long_options[option_index].name, \
              endTime);
          exit(1);
        }
        if (endTime > 999999999)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "GPS end time is after 14 September 2011 01:46:26 UTC " \
              "(%lld specified)\n", long_options[option_index].name, \
              endTime);
          exit(1);
        }

        break;

      case 'l':
        /* segment duration */
        segmentDuration = atoi(optarg);
        intervalDuration = 3 * segmentDuration;

        /* check */
        if (segmentDuration <= 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Segment duration must be greater than 0: (%d specified)\n", \
              long_options[option_index].name, segmentDuration);
          exit(1);
        }

        break;

      case 'a':
        /* sample rate */
        sampleRate = atoi(optarg);

        /* check */
        if (sampleRate < 2 || sampleRate > 16384 || sampleRate % 2)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Sample rate must be a power of 2 between 2 and 16384: " \
              "inclusive: (%d specified)\n", long_options[option_index].name, \
              sampleRate);
          exit(1);
        }

        break;

      case 'A':
        /* resample rate */
        resampleRate = atoi(optarg);

        /* check */
        if (resampleRate < 2 || resampleRate > 16384 || resampleRate % 2)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Resample rate must be a power of 2 between 2 and 16384: " \
              "inclusive: (%d specified)\n", long_options[option_index].name, \
              resampleRate);
          exit(1);
        }

        break;

      case 'f':
        /* minimal frequency */
        fMin = atoi(optarg);

        /* check */
        if (fMin < 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Minimum frequency is less than 0 Hz (%d specified)\n", \
              long_options[option_index].name, fMin);
          exit(1);
        }

        break;

      case 'F':
        /* maximal frequency */
        fMax = atoi(optarg);

        /* check */
        if (fMax < 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Maximum frequency is less than 0 Hz (%d specified)\n", \
              long_options[option_index].name, fMax);
          exit(1);
        }

        break;

      case 'w':
        /* hann window duration */
        hannDuration = atoi(optarg);

        /* check */
        if (hannDuration < 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Hann duartion is less than 0: (%d specified)\n", \
              long_options[option_index].name, hannDuration);
          exit(1);
        }

        break;

      case 'i':
        /* ifo for first stream */
        optarg_len = strlen(optarg) + 1;
        ifoOne = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        memcpy(ifoOne, optarg, optarg_len);

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
        optarg_len = strlen(optarg) + 1;
        ifoTwo = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        memcpy(ifoTwo, optarg, optarg_len);

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
        optarg_len = strlen(optarg) + 1;
        frameCacheOne = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        memcpy(frameCacheOne, optarg, optarg_len);
        break;

      case 'D':
        /* data cache two */
        optarg_len = strlen(optarg) + 1;
        frameCacheTwo = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        memcpy(frameCacheTwo, optarg, optarg_len);
        break;

      case 'r':
        /* calibration cache one */
        optarg_len = strlen(optarg) + 1;
        calCacheOne = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        memcpy(calCacheOne, optarg, optarg_len);
        break;

      case 'R':
        /* calibration cache two */
        optarg_len = strlen(optarg) + 1;
        calCacheTwo = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        memcpy(calCacheTwo, optarg, optarg_len);
        break;

      case 'z':
        /* set debug level */
        set_debug_level(optarg);

        /* check */
        if (atoi(optarg) < 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Debug level must be greater than 0: (%d specified)\n", \
              long_options[option_index].name, atoi(optarg));
          exit(1);
        }

        break;

      case 'k':
        /* high pass filter knee frequency */
        highPassFreq = atof(optarg);

        /* check */
        if (highPassFreq < 0)
        {
          fprintf(stderr, "Invalid argument tp --%s:\n" \
              "High pass filter knee frequency is less than 0 Hz: "\
              "(%f specified)\n", long_options[option_index].name, \
              highPassFreq);
          exit(1);
        }

        break;

      case 'p':
        /* high pass filter attenuation */
        highPassAtten = atof(optarg);

        /* check */
        if ((highPassAtten < 0.0) || (highPassAtten > 1.0))
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "High pass filter attenuation must be in the range [0:1]: " \
              "(%f specified)\n", long_options[option_index].name, \
              highPassAtten);
          exit(1);
        }

        break;

      case 'P':
        /* high pass filter order */
        highPassOrder = atoi(optarg);

        /* check */
        if (highPassOrder <= 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "High pass filter order must be greater than 0: " \
              "(%d specified)\n", long_options[option_index].name,
              highPassOrder);
          exit(1);
        }

        break;

      case 'e':
        /* filter omega_gw exponent */
        alpha = atof(optarg);
        break;

      case 'E':
        /* filter omega_gw reference frequency */
        fRef = atof(optarg);

        /* check */
        if (fRef < 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Reference frequency is less than 0 Hz: (%f specified)\n", \
              long_options[option_index].name, fRef);
          exit(1);
        }

        break;

      case 'O':
        /* filter omega_0 */
        omega0 = atof(optarg);

        /* check */
        if (omega0 <= 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Omega_0 is less than or equal to 0: (%f specified)\n", \
              long_options[option_index].name, omega0);
          exit(1);
        }

        break;

      case 'V':
        /* display version info and exit */
        fprintf(stdout, "Standalone SGWB Search Engine\n" CVS_ID "\n");
        exit(0);
        break;

      case '?':
        exit(1);
        break;

      default:
        fprintf(stderr, "Unknown error while parsing options\n");
        exit(1);
    }
  }

  if (optind < argc)
  {
    fprintf(stderr, "Extraneous command line arguments:\n");
    while(optind < argc)
    {
      fprintf(stderr, "%s\n", argv[optind++]);
    }
    exit(1);
  }

  /* check for required arguments */

  /* start time */
  if (startTime == 0)
  {
    fprintf(stderr, "--gps-start-time must be specified\n");
    exit(1);
  }

  if (endTime == 0)
  {
    fprintf(stderr, "--gps-end-time must be specified\n");
    exit(1);
  }

  /* segment duration */
  if (segmentDuration == -1)
  {
    fprintf(stderr, "--segment-duration must be specified\n");
    exit(1);
  }

  /* sample rate */
  if (sampleRate == -1)
  {
    fprintf(stderr, "--sample-rate must be specified\n");
    exit(1);
  }

  /* resample rate */
  if (resampleRate == -1)
  {
    fprintf(stderr, "--resample-rate must be specified\n");
    exit(1);
  }

  /* minimum frequency */
  if (fMin == -1)
  {
    fprintf(stderr, "--f-min must be specified\n");
    exit(1);
  }

  /* maximum frequency */
  if (fMax == -1)
  {
    fprintf(stderr, "--f-max must be specified\n");
    exit(1);
  }

  /* hann duration */
  if (hannDuration == -1)
  {
    fprintf(stderr, "--hann-duration must be specified\n");
    exit(1);
  }

  /* ifo */
  if (ifoOne == NULL)
  {
    fprintf(stderr, "--ifo-one must be specified\n");
    exit(1);
  }
  if (ifoTwo == NULL)
  {
    fprintf(stderr, "--ifo-two must be specified\n");
    exit(1);
  }

  /* frame cache */
  if (frameCacheOne == NULL)
  {
    fprintf(stderr, "--frame-cache-one must be specified\n");
    exit(1);
  }
  if (siteOne != siteTwo)
  {
    /* only need second frame cache if ifos differ */
    if (frameCacheOne == NULL)
    {
      fprintf(stderr, "--frame-cache-two must be specified\n");
      exit(1);
    }
  }

  /* calibration cache */
  if (calCacheOne == NULL)
  {
    fprintf(stderr, "--calibration-cache-one must be specified\n");
    exit(1);
  }
  if (calCacheTwo == NULL)
  {
    fprintf(stderr, "--calibration-cache-two must be specified\n");
    exit(1);
  }

  /* mask */
  if (apply_mask_flag)
  {
    if (maskBin == -1)
    {
      fprintf(stderr, "--mask-bin must be specified\n");
      exit(1);
    }
  }

  /* high pass filter */
  if (high_pass_flag)
  {
    if (highPassFreq == -1)
    {
      fprintf(stderr, "--hpf-frequency must be specified\n");
      exit(1);
    }
    if (highPassAtten == -1)
    {
      fprintf(stderr, "--hpf-attenuation must be specified\n");
      exit(1);
    }
    if (highPassOrder == -1)
    {
      fprintf(stderr, "--hpf-order must be specified\n");
      exit(1);
    }
  }

  /* check for sensible arguments */

  /* start time same as stop time */
  if (startTime == endTime)
  {
    fprintf(stderr, "Start time same as stop time; no analysis to perform\n");
    exit(1);
  }

  /* stop time before start time */
  if (startTime > endTime)
  {
    fprintf(stderr, "Invalid start/stop time; stop time (%lld) is before " \
        "start time (%lld)\n", endTime, startTime);
    exit(1);
  }

  /* resample rate greater than sample rate */
  if (resampleRate > sampleRate)
  {
    fprintf(stderr, "Invalid resample rate (%d); must be less than sample " \
        "rate (%d)\n", resampleRate, sampleRate);
    exit(1);
  }

  /* min frequency same as max */
  if (fMin == fMax)
  {
    fprintf(stderr, "Minimum frequency same as maximum; no analysis to " \
        "perform\n");
    exit(1);
  }

  /* max frequency less than min */
  if (fMin > fMax)
  {
    fprintf(stderr, "Invalid frequency band; maximum frequency (%d Hz) is " \
        "before minimum\nfrequency (%d Hz)\n", fMax, fMin);
    exit(1);
  }

  /* filter reference frequency less than min */
  if (fRef < fMin)
  {
    fprintf(stderr, "Reference frequency (%f Hz) is less than minimum " \
        "frequency (%d Hz)\n", fRef, fMin);
    exit(1);
  }

  /* filter reference frequency greater than max */
  if (fRef > fMax)
  {
    fprintf(stderr, "Reference frequency (%f Hz) is less than maximum " \
        "frequency (%d Hz)\n", fRef, fMax);
    exit(1);
  }

  /* hann duration greater than segment duration */
  if (hannDuration > segmentDuration)
  {
    fprintf(stderr, "Invalid hann duration (%d); must be less than, or " \
        "equal to segment\nduration (%d)\n", hannDuration, segmentDuration);
    exit(1);
  }

  return;
}

/* function to read data in frames */
static void readDataPair(LALStatus *status,
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
  INITSTATUS( status, "readDataPair", STOCHASTICDEVC );
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

  /* set resample parameters */
  resampleParams.deltaT = 1.0 / (REAL8)resampleRate;
  resampleParams.filterType = defaultButterworth;

  if (vrbflg)
  {
    fprintf(stdout, "Allocating memory for first raw data stream...\n");
  }

  /* allocate memory for first raw stream */
  dataStreamOne.data = NULL;
  LALSCreateVector(status->statusPtr, &(dataStreamOne.data), \
      sampleRate * (params->duration + (2 * buffer)));
  CHECKSTATUSPTR( status );
  memset(dataStreamOne.data->data, 0, \
      dataStreamOne.data->length * sizeof(*dataStreamOne.data->data));

  if (vrbflg)
  {
    fprintf(stdout, "Opening first frame cache...\n");
  }

  /* open first frame cache */
  LALFrCacheImport(status->statusPtr, &frCacheOne, params->frameCacheOne);
  CHECKSTATUSPTR( status );
  LALFrCacheOpen(status->statusPtr, &frStreamOne, frCacheOne);
  CHECKSTATUSPTR( status );

  /* set frame reading mode for first stream */
  LALFrSetMode(status->statusPtr, LAL_FR_SILENT_MODE, frStreamOne);
  CHECKSTATUSPTR( status );

  if (vrbflg)
  {
    fprintf(stdout, "Reading in channel \"%s\"...\n", frChanInOne.name);
  }

  /* read first channel */
  LALFrSeek(status->statusPtr, &(bufferStartTime), frStreamOne);
  CHECKSTATUSPTR( status );
  LALFrGetREAL4TimeSeries(status->statusPtr, &dataStreamOne, \
      &frChanInOne, frStreamOne);
  CHECKSTATUSPTR( status );

  /* resample */
  if (resampleRate != sampleRate)
  {
    if (vrbflg)
    {
      fprintf(stdout, "Resampling to %d Hz...\n", resampleRate);
    }

    /* resample */
    LALResampleREAL4TimeSeries(status->statusPtr, &dataStreamOne, \
        &resampleParams);
    CHECKSTATUSPTR( status );
  }

  if (vrbflg)
  {
    fprintf(stdout, "Allocating memory for second raw data stream...\n");
  }

  /* allocate memory for second raw stream */
  dataStreamTwo.data = NULL;
  LALSCreateVector(status->statusPtr, &(dataStreamTwo.data), \
      sampleRate * (params->duration + (2 * buffer)));
  CHECKSTATUSPTR( status );
  memset(dataStreamTwo.data->data, 0, \
      dataStreamTwo.data->length * sizeof(*dataStreamTwo.data->data));

  if (strcmp(params->frameCacheOne, params->frameCacheTwo) == 0)
  {
    if (vrbflg)
    {
      fprintf(stdout, "Reading in channel \"%s\" from same cache...\n", \
          frChanInTwo.name);
    }

    /* read in second channel */
    LALFrSeek(status->statusPtr, &(bufferStartTime), frStreamOne);
    CHECKSTATUSPTR( status );
    LALFrGetREAL4TimeSeries(status->statusPtr, &dataStreamTwo, \
        &frChanInTwo, frStreamOne);
    CHECKSTATUSPTR( status );

    if (vrbflg)
    {
      fprintf(stdout, "Closing frame cache...\n");
    }

    /* close frame cache */
    LALFrClose(status->statusPtr, &frStreamOne);
    CHECKSTATUSPTR( status );
  }
  else
  {
    if (vrbflg)
    {
      fprintf(stdout, "Closing first frame cache...\n");
    }

    /* close first frame cache */
    LALFrClose(status->statusPtr, &frStreamOne);
    CHECKSTATUSPTR( status );

    if (vrbflg)
    {
      fprintf(stdout, "Opening second frame cache...\n");
    }

    /* open second frame cache and read in second channel */
    LALFrCacheImport(status->statusPtr, &frCacheTwo, params->frameCacheTwo);
    CHECKSTATUSPTR( status );
    LALFrCacheOpen(status->statusPtr, &frStreamTwo, frCacheTwo);
    CHECKSTATUSPTR( status );

    /* set frame reading mode for second stream */
    LALFrSetMode(status->statusPtr, LAL_FR_SILENT_MODE, frStreamTwo);
    CHECKSTATUSPTR( status );

    if (vrbflg)
    {
      fprintf(stdout, "Reading in channel \"%s\"...\n", frChanInTwo.name);
    }

    /* read in second channel */
    LALFrSeek(status->statusPtr, &(bufferStartTime), frStreamTwo);
    CHECKSTATUSPTR( status );
    LALFrGetREAL4TimeSeries(status->statusPtr, &dataStreamTwo, \
        &frChanInTwo, frStreamTwo);
    CHECKSTATUSPTR( status );

    if (vrbflg)
    {
      fprintf(stdout, "Closing second frame cache...\n");
    }

    /* close second frame stream */
    LALFrClose(status->statusPtr, &frStreamTwo);
    CHECKSTATUSPTR( status );
  }

  /* resample */
  if (resampleRate != sampleRate)
  {
    if (vrbflg)
    {
      fprintf(stdout, "Resampling to %d Hz...\n", resampleRate);
    }

    /* resample */
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

/* function to estimate the psd */
static void psdEstimator(LALStatus *status,
    REAL4FrequencySeries output,
    PSDEstimatorInput input,
    PSDEstimatorParams params)
{
  /* counters */
  INT4 i;

  /* psd data structures */
  REAL4FrequencySeries psdTempA;
  REAL4FrequencySeries psdTempC;
  REAL4FrequencySeries psdA;
  REAL4FrequencySeries psdC;

  /* initialise status pointer */
  INITSTATUS( status, "psdEstimator", STOCHASTICDEVC );
  ATTATCHSTATUSPTR( status );

  /* allocate memory for psd data structures */
  psdTempA.data = NULL;
  psdTempC.data = NULL;
  LALCreateVector(status->statusPtr, &(psdTempA.data), params.psdTempLength);
  CHECKSTATUSPTR( status );
  memset(psdTempA.data->data, 0, \
      psdTempA.data->length * sizeof(*psdTempA.data->data));
  LALCreateVector(status->statusPtr, &(psdTempC.data), params.psdTempLength);
  CHECKSTATUSPTR( status );
  memset(psdTempC.data->data, 0, \
      psdTempC.data->length * sizeof(&psdTempC.data->data));

  /* allocate memory for reduced band psd data structures */
  psdA.data = (REAL4Sequence*)LALCalloc(1, sizeof(REAL4Sequence));
  psdC.data = (REAL4Sequence*)LALCalloc(1, sizeof(REAL4Sequence));
  psdA.data->length = params.filterLength;
  psdC.data->length = params.filterLength;

  /* compute uncalibrated PSDs */
  LALREAL4AverageSpectrum(status->statusPtr, &psdTempA, input.segmentA, \
      params.psdParams);
  CHECKSTATUSPTR( status );
  LALREAL4AverageSpectrum(status->statusPtr, &psdTempC, input.segmentC, \
      params.psdParams);
  CHECKSTATUSPTR( status );

  /* reduce to the optimal filter frequency range */
  psdA.data->data = psdTempA.data->data + params.numFMin;
  psdC.data->data = psdTempC.data->data + params.numFMin;

  for (i = 0; i < params.filterLength; i++)
  {
    /* compute inverse calibrated psds */
    psdA.data->data[i] = (pow(input.responseA->data->data[i].re, 2) + \
        pow(input.responseA->data->data[i].im, 2)) / psdA.data->data[i];
    psdC.data->data[i] = (pow(input.responseC->data->data[i].re, 2) + \
        pow(input.responseC->data->data[i].im, 2)) / psdC.data->data[i];

    /* average calibrated psds */
    output.data->data[i] = (psdA.data->data[i] + psdC.data->data[i]) / 2.0;
  }

  /* clean up */
  LALFree(psdA.data);
  LALFree(psdC.data);
  LALDestroyVector(status->statusPtr, &(psdTempA.data));
  CHECKSTATUSPTR( status );
  LALDestroyVector(status->statusPtr, &(psdTempC.data));
  CHECKSTATUSPTR( status );

  /* return status */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

/*
 * vim: et
 */
