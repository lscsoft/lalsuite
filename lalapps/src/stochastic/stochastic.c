/*
 * stochastic.c - SGWB Standalone Analysis Pipeline
 *
 * Tania Regimbau <Tania.Regimbau@astro.cf.ac.uk>
 * Adam Mercer <ram@star.sr.bham.ac.uk>
 *
 * $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>

#include <unistd.h>
#include <getopt.h>

#include <FrameL.h>

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

#include <lalapps.h>

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
static int middle_segment_flag = 0;
static int inject_flag = 0;
static int apply_mask_flag = 0;
static int high_pass_flag = 0;
static int overlap_hann_flag = 0;
static int verbose_flag = 0;
static int test_flag = 0;
static int post_analysis_flag = 0;
static int recenter_flag = 0;

/* parameters for the stochastic search */

/* sampling parameters */
INT4 sampleRate = 16384;
INT4 resampleRate = 1024;
REAL8 deltaF = 0.25;

/* data parameters */
LIGOTimeGPS gpsStartTime;
UINT8 startTime = 751957183;
UINT8 endTime = 751957783;
INT4 intervalDuration = 180;
INT4 segmentDuration = 60;
INT4 calibDuration = 60;
INT4 calibOffset = 29;
CHAR *frameCacheOne = NULL;
CHAR *frameCacheTwo = NULL;
CHAR *calCacheOne = NULL;
CHAR *calCacheTwo = NULL;
CHAR *channelOne = NULL;
CHAR *channelTwo = NULL;
CHAR *ifoOne = NULL;
CHAR *ifoTwo = NULL;
INT4 siteOne = 0;
INT4 siteTwo = 1;

/* frequency band */
INT4 fMin = 50;
INT4 fMax = 250;

/* omegaGW parameters */
REAL4 alpha = 0;
REAL4 fRef = 100;
REAL4 omegaRef = 1;

/* monte carlo parameters */
/* at the moment the code cannot do monte carlo with overlapped Hann window */
REAL4 scaleFactor = 1;
INT4 seed = 1;
INT4 NLoop = 1;

/* window parameters */
/* 60 s for pure Hann, 1 s for Tukey, 0s for rectangular window */
INT4 hannDuration = 1;

/* high pass filtering parameters */
REAL4 highPassFreq = 40;
REAL4 highPassAtten = 0.25;
INT4  highPassOrder = 6;

REAL4 geoScaleFactor = 1e18;

/* geo high pass filter parameters */
REAL4 geoHighPassFreq = 70;
INT4  geoHighPassOrder = 8;
REAL4 geoHighPassAtten = 0.9;

/* number of bins for frequency masking */
INT4 maskBin = 1;

/* arguments associated with test flag */
INT4 testInter = 0;
INT4 testSeg = 0;
INT4 testTrial = 0;

/* output file */
CHAR outputFilePath[200] = "/usr1/tregimba/";

INT4 main(INT4 argc, CHAR *argv[])
{
  /* variable declarations */

  /* status pointer */
  LALStatus status;

  /* output file */
  FILE *out1, *out2;
  CHAR outputFilename1[200], outputFilename2[200];

  /* counters */
  INT4 i, j, n, segLoop, interLoop, N;
  INT4 firstpass, pathology;

  /* results parameters */
  REAL8 y, yOpt;
  REAL8 varTheo, inVarTheoSum ;
  REAL8 ptEst, error;

  /* input data segment */
  INT4 duration, durationEff, extrasec;
  INT4 numSegments, numIntervals, segMiddle;
  INT4 segmentLength,segmentPadLength,intervalLength;
  INT4 segmentShift;
  INT4 padData;
  LIGOTimeGPS gpsStartPadTime, gpsCalibTime;
  ReadDataPairParams streamParams;
  StreamPair streamPair;
  REAL4TimeSeries segment1, segment2, segmentPad1, segmentPad2;
  REAL4Vector *seg1[100], *segPad1[100], *seg2[100], *segPad2[100];

  /* simulated signal structures */
  StochasticOmegaGWParameters parametersOmega;
  SSSimStochBGParams SBParams;
  SSSimStochBGInput SBInput;
  SSSimStochBGOutput SBOutput;
  REAL4TimeSeries SimStochBG1, SimStochBG2;

  REAL4FrequencySeries MComegaGW;
  COMPLEX8FrequencySeries MCresponse1, MCresponse2;
  COMPLEX8Vector *MCresp1[100], *MCresp2[100];

  INT4 MCLoop;
  INT4 MCfreqLength = 0;
  REAL8 MCdeltaF, MCdeltaT;
  LALUnit countPerStrain = {0,{0,0,0,0,0,-1,1},{0,0,0,0,0,0,0}};

  /* window for segment data streams */
  REAL4TimeSeries dataWindow;

  /* hann window */
  INT4 hannLength;
  LALWindowParams hannParams;
  REAL4Vector *hannWindow;

  /* high pass filtering */
  PassBandParamStruc highpassParam;

  /* response functions */
  COMPLEX8FrequencySeries responseTemp1, responseTemp2, response1, response2;
  COMPLEX8Vector *resp1[100], *resp2[100];
  INT4 respLength;
  CalibrationUpdateParams calfacts;
  LALUnit countPerAttoStrain = {18,{0,0,0,0,0,-1,1},{0,0,0,0,0,0,0}};
  FrCache *calibCache = NULL;

  /* data structures for PSDs */
  INT4 overlapPSDLength;
  INT4 psdTempLength;
  INT4 windowPSDLength;
  INT4 filterLength;
  INT4 numFMin, numFMax;
  LALWindowParams winparPSD;
  AverageSpectrumParams specparPSD;
  REAL4FrequencySeries psdTemp1, psdTemp2, psd1, psd2;
  REAL4Vector *calPsd1, *calPsd2;
  LALUnit psdUnits = {0,{0,0,1,0,0,0,2},{0,0,0,0,0,0,0}};

  /* calibrated inverse noise data structures */
  REAL4FrequencySeries calInvPsd1, calInvPsd2;

  /* units for inverse noise */
  LALUnit calPSDUnit = {36,{0,0,-1,0,0,-2,0},{0,0,0,0,0,0,0}};

  /* structures for LALInverseNoise */
  StochasticInverseNoiseInput inverseNoiseIn1, inverseNoiseIn2;
  StochasticInverseNoiseCalOutput inverseNoiseOut1, inverseNoiseOut2;

  /* zeropad and fft structures */
  SZeroPadAndFFTParameters zeroPadParams;
  RealFFTPlan *fftDataPlan = NULL;
  COMPLEX8FrequencySeries hBarTilde1, hBarTilde2;
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
  INT4 Nbin;

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

  /* spectrum structures */
  StochasticOmegaGWParameters omegaGWParams;
  REAL4FrequencySeries omegaGW;

  /* structures for CC spectrum and CC statistics */
  StochasticCrossCorrelationCalInput ccIn;
  BOOLEAN epochsMatch = 1;
  REAL4WithUnits ccStat;
  COMPLEX8FrequencySeries ccSpectrum;

  /* error handler */
  status.statusPtr = NULL;

  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "1" );

  /* parse command line options */
  parseOptions(argc, argv);

  if ((sampleRate == resampleRate)&&(high_pass_flag == 0))
    padData = 0;
  else
    padData = 1;
  
  /* get number of segments */
  numSegments = (INT4)(intervalDuration / segmentDuration);
  duration = endTime - startTime;
  numIntervals = (INT4)(((duration - (2 * padData)) / segmentDuration) - \
                 numSegments + 1);
  segMiddle = (INT4)((numSegments - 1) / 2);
  segmentShift = segmentDuration;

  /* recenter */
  if (recenter_flag)
  {
    durationEff = (INT4)(((duration - (2 * padData)) / segmentDuration ) * \
                  segmentDuration);
    extrasec = duration - durationEff;
    startTime = startTime + (INT4)(extrasec / 2);
    endTime = startTime + durationEff;
  }

  if (verbose_flag)
    fprintf(stdout, "%d\n", numIntervals);

  if (overlap_hann_flag)
  {
    /* numSegments = 2 * numSegments - 1; */
    segmentShift = segmentDuration / 2;
  }

  /* initialize gps time structure */
  gpsStartTime.gpsSeconds = startTime;
  gpsStartTime.gpsNanoSeconds = 0;
  gpsStartPadTime.gpsSeconds = startTime - padData;
  gpsStartPadTime.gpsNanoSeconds = 0;
  gpsCalibTime.gpsSeconds = startTime + calibOffset;
  gpsCalibTime.gpsNanoSeconds = 0;

  /* set length for data segments */
  intervalLength = (intervalDuration + (2 * padData)) * resampleRate;
  segmentPadLength = (segmentDuration + (2 * padData)) * resampleRate;
  segmentLength = segmentDuration * resampleRate;

  /* set metadata fields for data segments */
  strncpy(segmentPad1.name, "segmentPad1", LALNameLength);
  strncpy(segmentPad2.name, "segmentPad2", LALNameLength);
  segmentPad1.sampleUnits = segmentPad2.sampleUnits = lalADCCountUnit;
  segmentPad1.epoch = segmentPad2.epoch = gpsStartTime;
  segmentPad1.deltaT = segmentPad2.deltaT = 1./(REAL8)resampleRate;
  segmentPad1.f0 = segmentPad2.f0 = 0;
  strncpy(segment1.name, "segment1", LALNameLength);
  strncpy(segment2.name, "segment2", LALNameLength);
  segment1.sampleUnits = segment2.sampleUnits = lalADCCountUnit;
  segment1.epoch = segment2.epoch = gpsStartTime;
  segment1.deltaT = segment2.deltaT = 1./(REAL8)resampleRate;
  segment1.f0 = segment2.f0 = 0;

  if (verbose_flag)
    fprintf(stdout, "Allocating memory for data segments...\n");

  /* allocate memory for data segments */
  segmentPad1.data = segmentPad2.data = NULL;
  LAL_CALL(LALSCreateVector(&status, &(segmentPad1.data), \
        segmentPadLength), &status);
  LAL_CALL(LALSCreateVector(&status, &(segmentPad2.data), \
        segmentPadLength), &status);
  memset(segmentPad1.data->data, 0, \
      segmentPad1.data->length * sizeof(*segmentPad1.data->data));
  memset(segmentPad2.data->data, 0, \
      segmentPad2.data->length * sizeof(*segmentPad2.data->data));
  segment1.data = segment2.data = NULL;
  LAL_CALL(LALSCreateVector(&status, &(segment1.data), \
        segmentLength), &status);
  LAL_CALL(LALSCreateVector(&status, &(segment2.data), \
        segmentLength), &status);
  memset(segment1.data->data, 0, \
      segment1.data->length * sizeof(*segment1.data->data));
  memset(segment2.data->data, 0, \
      segment2.data->length * sizeof(*segment2.data->data));

  for (i = 0; i < numSegments; i++)
  {
    segPad1[i]= segPad2[i] = NULL;
    LAL_CALL(LALCreateVector(&status, &(segPad1[i]), \
          segmentPadLength), &status);
    LAL_CALL(LALCreateVector(&status, &(segPad2[i]), \
          segmentPadLength), &status);
    memset(segPad1[i]->data, 0, \
        segPad1[i]->length * sizeof(*segPad1[i]->data));
    memset(segPad2[i]->data, 0, \
        segPad2[i]->length * sizeof(*segPad2[i]->data));
  }

  for (i = 0; i < numSegments; i++)
  {
    seg1[i]= seg2[i] = NULL;
    LAL_CALL(LALCreateVector(&status, &(seg1[i]), segmentLength), &status);
    LAL_CALL(LALCreateVector(&status, &(seg2[i]), segmentLength), &status);
    memset(seg1[i]->data, 0, seg1[i]->length * sizeof(*seg1[i]->data));
    memset(seg2[i]->data, 0, seg2[i]->length * sizeof(*seg2[i]->data));
  }

  /* set segment input parameters */
  streamParams.frameCacheOne = frameCacheOne;
  streamParams.frameCacheTwo = frameCacheTwo;
  streamParams.ifoOne = ifoOne;
  streamParams.ifoTwo = ifoTwo;
  streamParams.channelOne = channelOne;
  streamParams.channelTwo = channelTwo;
  streamParams.buffer = 0;
  streamParams.sampleRate = sampleRate;
  streamParams.resampleRate = resampleRate;
  streamParams.duration = segmentDuration + 2 * padData;
  streamPair.streamOne = &segmentPad1;
  streamPair.streamTwo = &segmentPad2;

  if (inject_flag)
  {
    if (verbose_flag)
      fprintf(stdout, "Allocating memory for MC...\n");

    MCdeltaT = 1.0 / resampleRate;
    MCdeltaF = (REAL8)resampleRate / (REAL8)segmentPadLength;
    MCfreqLength = (segmentPadLength / 2) + 1;

    /* create vectors to store the simulated signal */
    strncpy(SimStochBG1.name, "Whitened-SimulatedSB1", LALNameLength);
    strncpy(SimStochBG2.name, "Whitened-SimulatedSB2", LALNameLength);
    SimStochBG1.f0 = SimStochBG2.f0 = 0.;
    SimStochBG1.epoch = SimStochBG2.epoch = gpsStartPadTime;
    SimStochBG1.deltaT = SimStochBG2.deltaT = 1./(REAL8)resampleRate;
    SimStochBG1.sampleUnits = SimStochBG2.sampleUnits = lalADCCountUnit;
    SimStochBG1.data = SimStochBG2.data = NULL;
    LAL_CALL(LALSCreateVector(&status, &(SimStochBG1.data), \
          segmentPadLength), &status);
    LAL_CALL(LALSCreateVector(&status, &(SimStochBG2.data), \
          segmentPadLength), &status);
    memset(SimStochBG1.data->data, 0, \
        SimStochBG1.data->length * sizeof(*SimStochBG1.data->data));
    memset(SimStochBG2.data->data, 0, \
        SimStochBG2.data->length * sizeof(*SimStochBG2.data->data));

    /* define parameters for SimulateSB */
    SBParams.length = segmentPadLength;
    SBParams.deltaT = 1. / resampleRate;
    SBParams.detectorOne = lalCachedDetectors[siteOne];
    SBParams.detectorTwo = lalCachedDetectors[siteTwo];
    SBParams.SSimStochBGTimeSeries1Unit = lalADCCountUnit;
    SBParams.SSimStochBGTimeSeries2Unit = lalADCCountUnit;

    /* omegaGW */
    parametersOmega.length = MCfreqLength;
    parametersOmega.f0 = 0;
    parametersOmega.deltaF = MCdeltaF;
    parametersOmega.alpha = alpha;
    parametersOmega.fRef = fRef;
    parametersOmega.omegaRef = omegaRef;

    /* allocate memory */
    MComegaGW.data = NULL;
    LAL_CALL(LALSCreateVector(&status, &(MComegaGW.data), MCfreqLength), \
        &status);
    memset(MComegaGW.data->data, 0, \
        MComegaGW.data->length * sizeof(*MComegaGW.data->data));

    /* generate omegaGW */
    LAL_CALL(LALStochasticOmegaGW(&status, &MComegaGW, &parametersOmega), \
        &status);

    /* response functions */
    memset(&calfacts, 0, sizeof(CalibrationUpdateParams));
    strncpy(MCresponse1.name,"MCresponse1", LALNameLength);
    strncpy(MCresponse2.name,"MCresponse2", LALNameLength);
    MCresponse1.sampleUnits = MCresponse2.sampleUnits = countPerStrain;
    MCresponse1.epoch = MCresponse2.epoch = gpsCalibTime;
    MCresponse1.deltaF = MCresponse2.deltaF = MCdeltaF;
    MCresponse1.f0 = MCresponse2.f0 = 0;
    MCresponse1.data = MCresponse2.data = NULL;
    LAL_CALL(LALCCreateVector(&status, &(MCresponse1.data), \
          MCfreqLength), &status);
    LAL_CALL(LALCCreateVector(&status, &(MCresponse2.data), \
          MCfreqLength), &status);
    memset(MCresponse1.data->data, 0, \
        MCresponse1.data->length * sizeof(*MCresponse1.data->data));
    memset(MCresponse2.data->data, 0, \
        MCresponse2.data->length * sizeof(*MCresponse2.data->data));

    for (i = 0; i < numSegments; i++)
    {
      MCresp1[i]= MCresp2[i] = NULL;
      LAL_CALL(LALCCreateVector(&status, &(MCresp1[i]), MCfreqLength), \
          &status);
      LAL_CALL(LALCCreateVector(&status, &(MCresp2[i]), MCfreqLength), \
          &status);
      memset(MCresp1[i]->data, 0, \
          MCresp1[i]->length * sizeof(*MCresp1[i]->data));
      memset(MCresp2[i]->data, 0, \
          MCresp2[i]->length * sizeof(*MCresp2[i]->data));
    }
  }

  /* set parameters for PSD estimation */
  windowPSDLength = (UINT4)(resampleRate / deltaF);
  overlapPSDLength = windowPSDLength / 2;
  psdTempLength = (windowPSDLength / 2) + 1;
  numFMin = (UINT4)(fMin / deltaF);
  numFMax = (UINT4)(fMax / deltaF);
  filterLength = numFMax - numFMin + 1;

  specparPSD.method = useMean;
  specparPSD.overlap = overlapPSDLength;
  specparPSD.plan = NULL;
  specparPSD.window = NULL;

  /* set window parameters for PSD estimation */
  winparPSD.length = windowPSDLength;
  winparPSD.type = Hann;

  /* set metadata fields for PSDs */
  strncpy(psdTemp1.name, "psdTemp1", LALNameLength);
  strncpy(psdTemp2.name, "psdTemp2", LALNameLength);
  psdTemp1.sampleUnits = psdTemp2.sampleUnits = psdUnits;
  psdTemp1.deltaF = psdTemp2.deltaF = deltaF;
  psdTemp1.f0 = psdTemp2.f0 = 0;

  if (verbose_flag)
    fprintf(stdout, "Allocating memory for PSDs...\n");

  /* allocate memory for PSDs */
  psdTemp1.data = psdTemp2.data = NULL;
  LAL_CALL(LALCreateVector(&status, &(psdTemp1.data), psdTempLength), \
      &status);
  LAL_CALL(LALCreateVector(&status, &(psdTemp2.data), psdTempLength), \
      &status);
  memset(psdTemp1.data->data, 0, \
      psdTemp1.data->length * sizeof(*psdTemp1.data->data));
  memset(psdTemp2.data->data, 0, \
      psdTemp2.data->length * sizeof(*psdTemp2.data->data));

  /* reduced frequency band PSDs */
  /* set metadata fields for reduced frequency band PSDs */
  strncpy(psd1.name, "psd1", LALNameLength);
  strncpy(psd2.name, "psd2", LALNameLength);
  psd1.deltaF = psd2.deltaF = deltaF;
  psd1.f0 = psd2.f0 = fMin;
  psd1.sampleUnits = psd2.sampleUnits = psdUnits;

  if (verbose_flag)
    fprintf(stdout, "Allocating memory for reduced frequency band PSDs...\n");
  
  /* allocate memory for reduced frequency band PSDs */
  psd1.data = psd2.data = NULL;
  LAL_CALL(LALCreateVector(&status, &psd1.data, filterLength), &status);
  LAL_CALL(LALCreateVector(&status, &psd2.data, filterLength), &status);
  memset(psd1.data->data, 0, psd1.data->length * sizeof(*psd1.data->data));
  memset(psd2.data->data, 0, psd2.data->length * sizeof(*psd2.data->data));

  /* allocate memory for calibrated PSDs */
  calPsd1 = calPsd2 = NULL;
  LAL_CALL(LALCreateVector(&status, &calPsd1, filterLength), &status);
  LAL_CALL(LALCreateVector(&status, &calPsd2, filterLength), &status);
  memset(calPsd1->data, 0, calPsd1->length * sizeof(*calPsd1->data));
  memset(calPsd2->data, 0, calPsd2->length * sizeof(*calPsd2->data));

  /* set parameters for response functions */
  respLength = (UINT4)(fMax / deltaF) + 1;

  /* set metadata fields for response functions */
  strncpy(responseTemp1.name, "responseTemp1", LALNameLength);
  strncpy(responseTemp2.name, "responseTemp2", LALNameLength);
  responseTemp1.sampleUnits = responseTemp2.sampleUnits = countPerAttoStrain;
  responseTemp1.epoch = responseTemp2.epoch = gpsCalibTime;
  responseTemp1.deltaF = responseTemp2.deltaF = deltaF;
  responseTemp1.f0 = responseTemp2.f0 = 0;

  if (verbose_flag)
    fprintf(stdout, "Allocating memory for response functions...\n");

  /* allocate memory for response functions */
  responseTemp1.data = responseTemp2.data = NULL;
  LAL_CALL(LALCCreateVector(&status, &(responseTemp1.data), respLength), \
      &status);
  LAL_CALL(LALCCreateVector(&status, &(responseTemp2.data), respLength), \
      &status);
  memset(responseTemp1.data->data, 0, \
      responseTemp1.data->length * sizeof(*responseTemp1.data->data));
  memset(responseTemp2.data->data, 0, \
      responseTemp2.data->length * sizeof(*responseTemp2.data->data));

  /* reduced frequency band response functions */
  /* set metadata fields for reduced frequency band response functions */
  strncpy(response1.name, "response1", LALNameLength);
  strncpy(response2.name, "response2", LALNameLength);
  response1.sampleUnits = response2.sampleUnits = countPerAttoStrain;
  response1.epoch = response2.epoch = gpsCalibTime;
  response1.deltaF = response2.deltaF = deltaF;
  response1.f0 = response2.f0 = fMin;

  if (verbose_flag)
  {
    fprintf(stdout, "Allocating memory for reduced frequency band " \
        "response functions...\n");
  }

  /* allocate memory for reduced frequency band response functions */
  response1.data = response2.data = NULL;
  LAL_CALL(LALCCreateVector(&status, &(response1.data), filterLength), \
      &status);
  LAL_CALL(LALCCreateVector(&status, &(response2.data), filterLength), \
      &status);
  memset(response1.data->data, 0, \
      response1.data->length * sizeof(*response1.data->data));
  memset(response2.data->data, 0, \
      response2.data->length * sizeof(*response2.data->data));

  for (i = 0; i < numSegments; i++)
  {
    resp1[i]= resp2[i] = NULL;
    LAL_CALL(LALCCreateVector(&status, &(resp1[i]), filterLength), &status);
    LAL_CALL(LALCCreateVector(&status, &(resp2[i]), filterLength), &status);
    memset(resp1[i]->data, 0, resp1[i]->length * sizeof(*resp1[i]->data));
    memset(resp2[i]->data, 0, resp2[i]->length * sizeof(*resp2[i]->data));
  }

  if (verbose_flag)
    fprintf(stdout, "Creating FFT plan for PSD estimation...\n");

  /* create fft plan */
  LAL_CALL(LALCreateForwardRealFFTPlan(&status, &specparPSD.plan, \
        windowPSDLength, 0), &status);

  if (verbose_flag)
    fprintf(stdout, "Creating window for PSD estimation...\n");

  /* create window for PSD estimation */
  LAL_CALL(LALCreateREAL4Window(&status, &specparPSD.window, &winparPSD), \
      &status);

  /* set metadata fields for inverse noise structures */
  strncpy(calInvPsd1.name, "calInvPsd1", LALNameLength);
  strncpy(calInvPsd2.name, "calInvPsd2", LALNameLength);
  calInvPsd1.sampleUnits = calInvPsd2.sampleUnits = calPSDUnit;
  calInvPsd1.deltaF = calInvPsd2.deltaF = deltaF;
  calInvPsd1.f0 = calInvPsd2.f0 = fMin;

  if (verbose_flag)
    fprintf(stdout, "Allocating memory for inverse noise...\n");

  /* allocate memory for inverse noise */
  calInvPsd1.data = calInvPsd2.data = NULL;
  LAL_CALL(LALCreateVector(&status, &(calInvPsd1.data), filterLength), \
      &status);
  LAL_CALL(LALCreateVector(&status, &(calInvPsd2.data), filterLength), \
      &status);
  memset(calInvPsd1.data->data, 0, \
      calInvPsd1.data->length * sizeof(*calInvPsd1.data->data));
  memset(calInvPsd2.data->data, 0, \
      calInvPsd2.data->length * sizeof(*calInvPsd2.data->data));

  /* set inverse noise inputs */
  inverseNoiseIn1.unCalibratedNoisePSD = &psd1;
  inverseNoiseIn1.responseFunction = &response1;
  inverseNoiseIn2.unCalibratedNoisePSD = &psd2;
  inverseNoiseIn2.responseFunction = &response2;

  /* set inverse noise outputs */
  inverseNoiseOut1.calibratedInverseNoisePSD = &calInvPsd1;
  inverseNoiseOut2.calibratedInverseNoisePSD = &calInvPsd2;

  /* set window parameters for segment data streams */
  strncpy(dataWindow.name, "dataWindow", LALNameLength);
  dataWindow.sampleUnits = lalDimensionlessUnit;
  dataWindow.deltaT = 1./resampleRate;
  dataWindow.f0 = 0;

  if (verbose_flag)
    fprintf(stdout, "Allocating memory for data segment window...\n");

  /* allocate memory for segment window */
  dataWindow.data = NULL;
  LAL_CALL(LALSCreateVector(&status, &(dataWindow.data), segmentLength), \
      &status);
  memset(dataWindow.data->data, 0, \
      dataWindow.data->length * sizeof(*dataWindow.data->data));

  if (verbose_flag)
    fprintf(stdout, "Generating data segment window...\n");

  /* generate window */
  for (i = 0; i < segmentLength; i++)
    dataWindow.data->data[i] = 1.;
  
  if (overlap_hann_flag)
    hannDuration = segmentDuration;

  if (hannDuration != 0)
  {
    /* generate pure Hann window */
    hannLength = hannDuration * resampleRate;
    hannParams.length = hannLength;
    hannParams.type = Hann;

    /* allocate memory for hann window */
    hannWindow = NULL;
    LAL_CALL(LALSCreateVector(&status, &hannWindow, hannLength), &status);
    memset(hannWindow->data, 0, \
        hannWindow->length * sizeof(*hannWindow->data));

    /* generate hann window */
    LAL_CALL(LALWindow(&status, hannWindow, &hannParams), &status);

    /* construct Tukey window */
    for (i = 0; i < hannLength / 2; i++)
      dataWindow.data->data[i] = hannWindow->data[i];

    for (i = segmentLength - (hannLength / 2); i < segmentLength; i++)
    {
      dataWindow.data->data[i] = \
        hannWindow->data[i - segmentLength + hannLength];
    }
  }

  /* print window */
  if (test_flag)
    LALSPrintTimeSeries(&dataWindow, "dataWindow.dat");

  /* structure for high pass filtering */
  if (high_pass_flag)
  {
    highpassParam.nMax = highPassOrder;
    highpassParam.f1 = -1;
    highpassParam.f2 = highPassFreq;
    highpassParam.a1 = -1;
    highpassParam.a2 = highPassAtten;
  }

  /* zeropad lengths */
  zeroPadLength = 2 * segmentLength;
  fftDataLength = (zeroPadLength / 2) + 1;

  /* create fft plan */
  LAL_CALL(LALCreateForwardRealFFTPlan(&status, &fftDataPlan, \
        zeroPadLength, 0), &status);

  /* set metadata fields for zeropad ffts */
  strncpy(hBarTilde1.name, "hBarTilde1", LALNameLength);
  strncpy(hBarTilde2.name, "hBarTilde2", LALNameLength);

  if (verbose_flag)
    fprintf(stdout, "Allocating memory for zeropad...\n");

  /* allocate memory for zeropad */
  hBarTilde1.data = hBarTilde2.data = NULL;
  LAL_CALL(LALCCreateVector(&status, &(hBarTilde1.data), fftDataLength), \
      &status);
  LAL_CALL(LALCCreateVector(&status, &(hBarTilde2.data), fftDataLength), \
      &status);
  memset(hBarTilde1.data->data, 0, \
      hBarTilde1.data->length * sizeof(*hBarTilde1.data->data));
  memset(hBarTilde2.data->data, 0, \
      hBarTilde2.data->length * sizeof(*hBarTilde2.data->data));

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
  LAL_CALL(LALCreateVector(&status, &(overlap.data), filterLength), &status);
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
    fprintf(stdout, "Generating the overlap reduction function...\n");

  /* generate overlap reduction function */
  LAL_CALL(LALOverlapReductionFunction(&status, &overlap, &detectors, \
        &ORFparams), &status);

  /* print */
  if (test_flag)
    LALSPrintFrequencySeries(&overlap, "overlap.dat");

  /* set metadata fields for spectrum */
  strncpy(omegaGW.name, "omegaGW", LALNameLength);
  omegaGW.sampleUnits = lalDimensionlessUnit;
  omegaGW.deltaF = deltaF;
  omegaGW.f0 = fMin;

  if (verbose_flag)
    fprintf(stdout, "Allocating memory for spectrum...\n");

  /* allocate memory for spectrum */
  omegaGW.data = NULL;
  LAL_CALL(LALCreateVector(&status, &(omegaGW.data), filterLength), &status);
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
    fprintf(stdout, "Generating spectrum for optimal filter...\n");

  /* generage omegaGW */
  LAL_CALL(LALStochasticOmegaGW(&status, &omegaGW, &omegaGWParams), \
      &status);

  /* frequency mask */
  if (apply_mask_flag)
  {
    /* extra bins */
    Nbin = (maskBin - 1) / 2;

    /* set metadata fields for frequency mask */
    strncpy(mask.name, "mask", LALNameLength);
    mask.deltaF = deltaF;
    mask.f0 = fMin;
    mask.sampleUnits = lalDimensionlessUnit;

    if (verbose_flag)
      fprintf(stdout, "Allocating memory for frequency mask...\n");

    /* allocate memory for frequency mask */
    mask.data = maskTemp = NULL;
    LAL_CALL(LALCreateVector(&status, &(mask.data), filterLength), &status);
    LAL_CALL(LALCreateVector(&status, &maskTemp, respLength), &status);
    memset(mask.data->data, 0, mask.data->length * sizeof(*mask.data->data));
    memset(maskTemp->data, 0, maskTemp->length * sizeof(*maskTemp->data));

    if (verbose_flag)
      fprintf(stdout, "Generating frequency mask...\n");

    /* set all values to 1 */
    for (i = 0; i < respLength; i++)
      maskTemp->data[i] = 1.;

    if (verbose_flag)
      fprintf(stdout, "Masking multiples of 16 Hz...\n");

    /* remove multiples of 16 Hz */
    for (i = 0; i < respLength; i += (UINT4)(16 / deltaF))
    {
      maskTemp->data[i]= 0;
      for (j = 0; j < Nbin; j++)
      {
        if ((i + 1 + j) < respLength)
          maskTemp->data[i + 1 + j]= 0;
        if ((i - 1 - j) > 0 )
          maskTemp->data[i - 1 - j]= 0;
      }
    }

    if (verbose_flag)
      fprintf(stdout, "Masking multiples of 60 Hz...\n");

    /* remove multiples of 60 Hz */
    for (i = 0; i < respLength; i += (UINT4)(60 / deltaF))
    {
      maskTemp->data[i] = 0;
      for (j = 0; j < Nbin; j ++)
      {
        if ((i + 1 + j) < respLength)
          maskTemp->data[i + 1 + j]= 0;
        if ((i - 1 - j) > 0 )
          maskTemp->data[i - 1 - j]= 0;
      }
    }

    if (verbose_flag)
      fprintf(stdout, "Getting appropriate frequency band for mask...\n");

    /* get appropriate band */
    for (i = 0; i < filterLength; i++)
      mask.data->data[i] = maskTemp->data[i + numFMin];

    if (test_flag)
      LALSPrintFrequencySeries(&mask, "mask.dat");
    
    if (verbose_flag)
      fprintf(stdout, "Applying frequency mask to spectrum..\n");

    /* apply mask to omegaGW */
    for (i = 0; i < filterLength; i++)
      omegaGW.data->data[i] *= mask.data->data[i];
  }

  /* print */
  if (test_flag)
    LALSPrintFrequencySeries(&omegaGW, "omegaGW.dat");

  /* set normalisation parameters */
  normParams.fRef = fRef;
  normParams.heterodyned = 0;
  normParams.window1 = normParams.window2 = dataWindow.data;

  /* set normalisation input */
  normInput.overlapReductionFunction = &overlap;
  normInput.omegaGW = &omegaGW;
  normInput.inverseNoisePSD1 = &calInvPsd1;
  normInput.inverseNoisePSD2 = &calInvPsd2;

  /* set normalisation output */
  normOutput.normalization = &normLambda;
  normOutput.variance = &normSigma;

  /* set metadata fields for optimal filter */
  strncpy(optFilter.name, "optFilter", LALNameLength);
  optFilter.epoch = gpsStartTime;
  optFilter.deltaF = deltaF;
  optFilter.f0 = fMin;

  if (verbose_flag)
    fprintf(stdout, "Allocating memory for optimal filter...\n");

  /* allocate memory for optimal filter */
  optFilter.data = NULL;
  LAL_CALL(LALCreateVector(&status, &(optFilter.data), filterLength), \
      &status);
  memset(optFilter.data->data, 0, \
      optFilter.data->length * sizeof(*optFilter.data->data));

  /* set optimal filter inputs */
  optFilterIn.overlapReductionFunction = &overlap;
  optFilterIn.omegaGW = &omegaGW;
  optFilterIn.calibratedInverseNoisePSD1 = &calInvPsd1;
  optFilterIn.calibratedInverseNoisePSD2 = &calInvPsd2;

  /* set metadata fields for CC spectrum */
  strncpy(ccSpectrum.name, "ccSpectrum", LALNameLength);
  ccSpectrum.epoch = gpsStartTime;
  ccSpectrum.deltaF = deltaF;
  ccSpectrum.f0 = fMin;

  /* allocate memory for CC spectrum*/
  ccSpectrum.data = NULL;
  LAL_CALL(LALCCreateVector(&status, &(ccSpectrum.data), filterLength), \
      &status);
  memset(ccSpectrum.data->data, 0, \
      ccSpectrum.data->length * sizeof(*ccSpectrum.data->data));

  /* set CC inputs */
  ccIn.hBarTildeOne = &hBarTilde1;
  ccIn.hBarTildeTwo = &hBarTilde2;
  ccIn.responseFunctionOne = &response1;
  ccIn.responseFunctionTwo = &response2;
  ccIn.optimalFilter = &optFilter;

  /*** DONE HERE WITH ALLOCATION ***/

  if (overlap_hann_flag)
    N = 2;
  else
    N = 1;

  for (MCLoop = 0; MCLoop < NLoop; MCLoop++)
  {
    /* initialize parameters for post analysis */
    yOpt = 0;
    inVarTheoSum = 0;

    /* open output file */
    LALSnprintf(outputFilename1, LALNameLength, \
        "%s/stat-%s%s-%lld-%lld-%d.dat", outputFilePath, ifoOne, ifoTwo, \
        startTime, endTime, MCLoop);
    /* open output file */
    LALSnprintf(outputFilename2, LALNameLength, \
        "%s/post-%s%s-%lld-%lld-%d.dat", outputFilePath, ifoOne, ifoTwo, \
        startTime, endTime, MCLoop);

    for (n = 0; n < N; n++)
    {
      startTime = startTime + (n * segmentShift);
      firstpass = 1;
      pathology = 0;

      for (interLoop = 0; interLoop < numIntervals; interLoop++)
      {	
        if (firstpass)
        {
          if (verbose_flag)
          {
            fprintf(stdout, "First Pass\ninterval %d out of %d\n", \
                interLoop, numIntervals);
          }

          /* read first interval and get response functions */
          lal_errhandler = LAL_ERR_RTRN;

          for (segLoop = 0; segLoop < numSegments; segLoop++)
          {
            /* define segment epoch */
            gpsStartTime.gpsSeconds = startTime + (interLoop + segLoop) * \
                                      segmentDuration;
            gpsStartPadTime.gpsSeconds = gpsStartTime.gpsSeconds - padData;
            segmentPad1.epoch = segmentPad2.epoch = gpsStartPadTime;

            if (verbose_flag)
            {
              fprintf(stdout, "request data at GPS time %d\n", \
                  gpsStartTime.gpsSeconds);
            }

            streamParams.startTime = gpsStartPadTime.gpsSeconds;
            LAL_CALL(readDataPair(&status, &streamPair, &streamParams), \
                &status);

            /* skip segment if data not found or corrupted with 0 values */
            if (status.statusCode != 0)
            {
              clear_status(&status);
              firstpass = 1;
              pathology = 1;
              interLoop = interLoop + segLoop;
              break;
            }

            /* store in memory */
            for (i = 0; i < segmentPadLength; i++)
            {
              segPad1[segLoop]->data[i] = segmentPad1.data->data[i];
              segPad2[segLoop]->data[i] = segmentPad2.data->data[i];
            }

            /* compute response functions */
            gpsCalibTime.gpsSeconds = gpsStartTime.gpsSeconds  + calibOffset;
            response1.epoch  = response2.epoch  = gpsCalibTime;

            if (strncmp(ifoOne, "G1", 2) == 0)
            {
              for (i = 0; i < filterLength; i++)
              {
                response1.data->data[i].re = 1;
                response1.data->data[i].im = 0;
              }
            }
            else
            {
              responseTemp1.epoch = gpsCalibTime;
              if (verbose_flag)
              {
                fprintf(stdout, "request GPS time %d for ifo 1\n", \
                    gpsCalibTime.gpsSeconds );
              }
              memset(&calfacts, 0, sizeof(CalibrationUpdateParams));
              calfacts.ifo = ifoOne;
              LAL_CALL(LALFrCacheImport(&status, &calibCache, calCacheOne), \
                  &status);
              LAL_CALL(LALExtractFrameResponse(&status, &responseTemp1, \
                    calibCache, &calfacts), &status);
              LAL_CALL(LALDestroyFrCache(&status, &calibCache), &status);

              /* skip segment if data not found or corrupted with 0 values */
              if ((status.statusCode !=0)||(responseTemp1.data==NULL))
              {
                clear_status(&status);
                firstpass = 1;
                pathology = 1;
                interLoop = interLoop + segLoop;
                break;
              }

              /* reduce to the optimal filter frequency range */
              for (i = 0; i < filterLength; i++)
              {
                response1.data->data[i] = responseTemp1.data->data[i + \
                                          numFMin];
              }
            }

            if (strncmp(ifoTwo, "G1", 2) == 0)
            {
              for (i = 0; i < filterLength; i++)
              {
                response2.data->data[i].re = 1;
                response2.data->data[i].im = 0;
              }
            }
            else
            {
              responseTemp2.epoch = gpsCalibTime;
              if (verbose_flag)
              {
                fprintf(stdout, "request GPS time %d for ifo 2\n", \
                    gpsCalibTime.gpsSeconds );
              }
              memset(&calfacts, 0, sizeof(CalibrationUpdateParams));
              calfacts.ifo = ifoTwo;
              LAL_CALL(LALFrCacheImport(&status, &calibCache, calCacheTwo ), \
                  &status);
              LAL_CALL(LALExtractFrameResponse( &status, &responseTemp2, \
                    calibCache, &calfacts), &status);
              LAL_CALL(LALDestroyFrCache( &status, &calibCache), &status);

              /* skip segment if data not found or corrupted with 0 values */
              if ((status.statusCode !=0)||(responseTemp2.data==NULL))
              {
                clear_status(&status);
                firstpass = 1;
                pathology = 1;
                interLoop = interLoop + segLoop;
                break;
              }

              /* reduce to the optimal filter frequency range */
              for (i = 0; i < filterLength; i++)
              {
                response2.data->data[i] = responseTemp2.data->data[i + \
                                          numFMin];
              }
            }

            /* store in memory */
            for (i = 0; i < filterLength; i++)
            {
              resp1[segLoop]->data[i] = response1.data->data[i];
              resp2[segLoop]->data[i] = response2.data->data[i];
            }

            /* convert response function for use in the MC routine */
            if (inject_flag)
            {
              MCresponse1.epoch = MCresponse2.epoch = gpsCalibTime;
              LAL_CALL(LALResponseConvert(&status, &MCresponse1, \
                    &responseTemp1), &status);
              LAL_CALL(LALResponseConvert(&status, &MCresponse2, \
                    &responseTemp2), &status);

              /* force DC to be 0 and nyquist to be real */
              MCresponse1.data->data[0].re = MCresponse2.data->data[0].re = 0;
              MCresponse1.data->data[0].im = MCresponse2.data->data[0].im = 0;
              MCresponse1.data->data[MCfreqLength-1].im = 0;
              MCresponse2.data->data[MCfreqLength-1].im = 0;

              /* store in memory */
              for (i = 0; i < MCfreqLength; i++)
              {
                MCresp1[segLoop]->data[i] = MCresponse1.data->data[i];
                MCresp2[segLoop]->data[i] = MCresponse2.data->data[i];
              }
            }
          }

          if (pathology == 1)
          {
            if (verbose_flag)
              fprintf(stdout, "\n \n BREAK\n \n");

            pathology = 0;
            firstpass = 1;

            continue;
          }
          else
            firstpass = 0;
        }
        /* if firstpass = 0, shift data and read extra segment */
        else
        {
          lal_errhandler = LAL_ERR_RTRN;
          if (verbose_flag)
          {
            fprintf(stdout, "Shift Segment\ninterval %d out of %d\n", \
                interLoop, numIntervals);
          }
          /* shift segments */
          for (segLoop = 0; segLoop < numSegments - 1; segLoop++)
          {
            for (i = 0; i < segmentPadLength; i++)
            {
              segPad1[segLoop]->data[i] = segPad1[segLoop + 1]->data[i];
              segPad2[segLoop]->data[i] = segPad2[segLoop + 1]->data[i];
            }
            for (i = 0; i < filterLength; i++)
            {
              resp1[segLoop]->data[i] = resp1[segLoop + 1]->data[i];
              resp2[segLoop]->data[i] = resp2[segLoop + 1]->data[i];

              if (inject_flag)
              {
                MCresp1[segLoop]->data[i] = MCresp1[segLoop + 1]->data[i];
                MCresp2[segLoop]->data[i] = MCresp2[segLoop + 1]->data[i];
              }
            }
          }

          /* read extra segment */
          gpsStartTime.gpsSeconds = startTime + (interLoop + \
              numSegments - 1) * segmentDuration;
          gpsStartPadTime.gpsSeconds = gpsStartTime.gpsSeconds - padData;
          streamParams.startTime = gpsStartPadTime.gpsSeconds;

          if (verbose_flag)
          {
            fprintf( stdout, "read end segment at GPS %d... \n", \
                gpsStartPadTime.gpsSeconds);
          }

          LAL_CALL(readDataPair(&status, &streamPair, &streamParams), &status);

          /* skip segment if data not found or corrupted with 0 values */
          if (status.statusCode != 0)
          {
            clear_status(&status);
            firstpass = 1;
            interLoop = interLoop + (numSegments -1);
            continue;
          }

          /* store in memory */
          for (i = 0; i < segmentPadLength ; i++)
          {
            segPad1[numSegments-1]->data[i] = segmentPad1.data->data[i];
            segPad2[numSegments-1]->data[i] = segmentPad2.data->data[i];
          }

          /* compute extra response function */
          gpsCalibTime.gpsSeconds = gpsStartTime.gpsSeconds + calibOffset;
          response1.epoch  = response2.epoch = gpsCalibTime;

          if (strncmp(ifoOne, "G1", 2) == 0)
          {
            for (i = 0; i < filterLength; i++)
            {
              response1.data->data[i].re = 1;
              response1.data->data[i].im = 0;
            }
          }
          else
          {
            responseTemp1.epoch = gpsCalibTime;
            if (verbose_flag)
            {
              fprintf(stdout, "request GPS time %d for ifo 1\n", \
                  gpsCalibTime.gpsSeconds );
            }
            memset(&calfacts, 0, sizeof(CalibrationUpdateParams));
            calfacts.ifo = ifoOne;
            LAL_CALL(LALFrCacheImport(&status, &calibCache, calCacheOne ), \
                &status);
            LAL_CALL(LALExtractFrameResponse(&status, &responseTemp1, \
                  calibCache, &calfacts), &status);
            LAL_CALL(LALDestroyFrCache(&status, &calibCache), &status);

            /* skip segment if data not found or corrupted with 0 values */
            if ((status.statusCode != 0) || (responseTemp1.data == NULL))
            {
              clear_status(&status);
              firstpass = 1;
              pathology = 1;
              interLoop = interLoop + segLoop;
              break;
            }

            /* reduce to the optimal filter frequency range */
            for (i = 0; i < filterLength; i++)
            {
              response1.data->data[i] = responseTemp1.data->data[i + \
                                        numFMin];
            }
          }

          if (strncmp(ifoTwo, "G1", 2) == 0)
          {
            for (i = 0; i < filterLength; i++)
            {
              response1.data->data[i].re = 1;
              response1.data->data[i].im = 0;
            }
          }
          else
          {
            responseTemp2.epoch = gpsCalibTime;
            if (verbose_flag)
            {
              fprintf(stdout, "request GPS time %d for ifo 2\n", \
                  gpsCalibTime.gpsSeconds );
            }
            memset(&calfacts, 0, sizeof(CalibrationUpdateParams));
            calfacts.ifo = ifoTwo;
            LAL_CALL(LALFrCacheImport(&status, &calibCache, calCacheTwo ), \
                &status);
            LAL_CALL(LALExtractFrameResponse(&status, &responseTemp2, \
                  calibCache, &calfacts), &status);
            LAL_CALL(LALDestroyFrCache(&status, &calibCache), &status);

            /* skip segment if data not found or corrupted with 0 values */
            if ((status.statusCode != 0) || (responseTemp2.data == NULL))
            {
              clear_status(&status);
              firstpass = 1;
              pathology = 1;
              interLoop = interLoop + segLoop;
              break;
            }

            /* reduce to the optimal filter frequency range */
            response2.epoch  = gpsCalibTime;
            for (i = 0; i < filterLength; i++)
            {
              response2.data->data[i] = responseTemp2.data->data[i + \
                                        numFMin];
            }
          }

          /* store in memory */
          for (i = 0; i < filterLength; i++)
          {
            resp1[numSegments-1]->data[i] = response1.data->data[i];
            resp2[numSegments-1]->data[i] = response2.data->data[i];
          }

          /* convert response function for use in the MC routine */
          if (inject_flag)
          {
            MCresponse1.epoch = MCresponse2.epoch = gpsCalibTime;
            LAL_CALL(LALResponseConvert(&status, &MCresponse1, \
                  &responseTemp1), &status);
            LAL_CALL(LALResponseConvert(&status, &MCresponse2, \
                  &responseTemp2), &status);

            /* force DC to be 0 and nyquist to be real */
            MCresponse1.data->data[0].re = MCresponse2.data->data[0].re = 0;
            MCresponse1.data->data[0].im = MCresponse2.data->data[0].im = 0;
            MCresponse1.data->data[MCfreqLength-1].im = 0;
            MCresponse2.data->data[MCfreqLength-1].im = 0;

            /* store in memory */
            for (i = 0; i < MCfreqLength; i++)
            {
              MCresp1[numSegments-1]->data[i] = MCresponse1.data->data[i];
              MCresp2[numSegments-1]->data[i] = MCresponse2.data->data[i];
            }
          }
        }

        /* initialize average PSDs */
        for (i = 0; i < filterLength; i++)
        {
          calPsd1->data[i] = 0;
          calPsd2->data[i] = 0;
        }

        for (segLoop = 0; segLoop < numSegments; segLoop++)
        {
          gpsStartTime.gpsSeconds = startTime + (interLoop + segLoop) * \
                                    segmentDuration;
          gpsStartPadTime.gpsSeconds = gpsStartTime.gpsSeconds - padData;
          segmentPad1.epoch = segmentPad2.epoch = gpsStartPadTime;
          segment1.epoch = segment2.epoch = gpsStartTime;
          gpsCalibTime.gpsSeconds = gpsStartTime.gpsSeconds  + calibOffset;
          response1.epoch = response2.epoch = gpsCalibTime;

          for (i = 0; i < filterLength; i++)
          {
            response1.data->data[i] = resp1[segLoop]->data[i];
            response2.data->data[i] = resp2[segLoop]->data[i];
          }

          /* print */
          if ((test_flag) && (interLoop == testInter) && \
              (segLoop == testSeg) && (MCLoop == testTrial))
          {
            LALCPrintFrequencySeries(&response1, "response1.dat");
            LALCPrintFrequencySeries(&response2, "response2.dat");
          }

          /* simulate signal */
          if (inject_flag)
          {
            for (i = 0; i < MCfreqLength; i++)
            {
              MCresponse1.data->data[i] = MCresp1[segLoop]->data[i];
              MCresponse2.data->data[i] = MCresp2[segLoop]->data[i];
            }

            /* set parameters for monte carlo */
            SimStochBG1.epoch = SimStochBG2.epoch = gpsStartPadTime;
            SBParams.seed = seed;

            /* define input structure for SimulateSB */
            SBInput.omegaGW = &MComegaGW;
            SBInput.whiteningFilter1 = &MCresponse1;
            SBInput.whiteningFilter2 = &MCresponse2;

            /* define output structure for SimulateSB */
            SBOutput.SSimStochBG1 = &SimStochBG1;
            SBOutput.SSimStochBG2 = &SimStochBG2;

            /* perform monte carlo */
            LAL_CALL(LALSSSimStochBGTimeSeries(&status, &SBOutput, \
                  &SBInput, &SBParams), &status);

            /* print */
            if ((test_flag) && (interLoop == testInter) && \
                (segLoop == testSeg) && (MCLoop == testTrial))
            {
              LALSPrintTimeSeries(&SimStochBG1, "simStochBG1.dat");
              LALSPrintTimeSeries(&SimStochBG2, "simStochBG2.dat");
            }

            /* multiply by scale factor and inject into real data */
            for (i = 0; i < segmentPadLength ; i++)
            {
              segmentPad1.data->data[i] = segPad1[segLoop]->data[i] + \
                (scaleFactor * SimStochBG1.data->data[i]);
              segmentPad2.data->data[i] = segPad2[segLoop]->data[i] + \
                (scaleFactor * SimStochBG2.data->data[i]);
            }

            /* print */
            if ((test_flag) && (interLoop == testInter) && \
                (segLoop == testSeg) && (MCLoop == testTrial))
            {
              LALSPrintTimeSeries(&segmentPad1, "segmentPadMC1.dat");
              LALSPrintTimeSeries(&segmentPad2, "segmentPadMC2.dat");
            }

            /* increase seed */
            seed = seed + 2;
          }
          else
          {
            for (i = 0; i < segmentPadLength; i++)
            {
              segmentPad1.data->data[i] = segPad1[segLoop]->data[i];
              segmentPad2.data->data[i] = segPad2[segLoop]->data[i];
            }
          }

          /* print */
          if ((test_flag) && (interLoop == testInter) && \
              (segLoop == testSeg) && (MCLoop == testTrial))
          {
            LALSPrintTimeSeries(&segmentPad1, "segmentPad1.dat");
            LALSPrintTimeSeries(&segmentPad2, "segmentPad2.dat");
          }

          /* high pass fitering */
          if (high_pass_flag)
          {
            LAL_CALL(LALButterworthREAL4TimeSeries(&status, &segmentPad1, \
                  &highpassParam), &status);
            LAL_CALL(LALButterworthREAL4TimeSeries(&status, &segmentPad2, \
                  &highpassParam), &status);
          }

          /* throw away pad data on each side of the segment */
          for (i = 0; i < segmentLength; i++)
          {
            segment1.data->data[i] = segmentPad1.data->data[i + padData * \
                                     resampleRate];
            segment2.data->data[i] = segmentPad2.data->data[i + padData * \
                                     resampleRate];
          }

          /* print */
          if ((test_flag) && (interLoop == testInter) && \
              (segLoop == testSeg) && (MCLoop == testTrial))
          {
            LALSPrintTimeSeries(&segment1, "segment1.dat");
            LALSPrintTimeSeries(&segment2, "segment2.dat");
          }

          /* store in memory */

          for (i = 0; i < segmentLength; i++)
          {
            seg1[segLoop]->data[i] = segment1.data->data[i];
            seg2[segLoop]->data[i] = segment2.data->data[i];
          }

          if ((middle_segment_flag == 0) && (segLoop == segMiddle))
          {
            if (verbose_flag)
              fprintf(stdout, "ignoring middle segment..\n");
          }
          else
          {
            if (verbose_flag)
              fprintf(stdout, "Estimating PSDs...\n");

            /* compute uncalibrated PSDs */
            LAL_CALL(LALREAL4AverageSpectrum(&status, &psdTemp1, &segment1, \
                  &specparPSD), &status);
            LAL_CALL(LALREAL4AverageSpectrum(&status, &psdTemp2, &segment2, \
                  &specparPSD), &status);

            if (verbose_flag)
            {
              fprintf(stdout, "Getting appropriate frequency band for " \
                  "PSDs..\n");
            }

            /* reduce to the optimal filter frequency range */
            for (i = 0; i < filterLength; i++)
            {
              psd1.data->data[i] =  psdTemp1.data->data[i + numFMin];
              psd2.data->data[i] =  psdTemp2.data->data[i + numFMin];
            }

            /* print */
            if ((test_flag) && (interLoop == testInter) && \
                (segLoop == testSeg) && (MCLoop == testTrial))
            {
              LALSPrintFrequencySeries(&psd1, "psd1.dat");
              LALSPrintFrequencySeries(&psd2, "psd2.dat");
            }

            if (verbose_flag)
              fprintf(stdout, "Generating inverse noise...\n");

            /* compute inverse calibrate PSDs */
            LAL_CALL(LALStochasticInverseNoiseCal(&status, \
                  &inverseNoiseOut1, &inverseNoiseIn1), &status);
            LAL_CALL(LALStochasticInverseNoiseCal(&status, \
                  &inverseNoiseOut2, &inverseNoiseIn2), &status);

            /* print */
            if ((test_flag) && (interLoop == testInter) && \
                (segLoop == testSeg) && (MCLoop == testTrial))
            {
              LALSPrintFrequencySeries(&calInvPsd1, "calInvPsd1.dat");
              LALSPrintFrequencySeries(&calInvPsd2, "calInvPsd2.dat");
            }

            /* sum over calibrated PSDs for average */
            for (i = 0; i < filterLength; i++)
            {
              calPsd1->data[i] = calPsd1->data[i] + \
                                 1. / calInvPsd1.data->data[i];
              calPsd2->data[i] = calPsd2->data[i] + \
                                 1. / calInvPsd2.data->data[i];
            }
          }
        }

        /* average calibrated PSDs and take inverse */
        for (i = 0; i < filterLength; i++)
        {
          if (middle_segment_flag == 0)
          {
            calPsd1->data[i] = calPsd1->data[i] / (REAL4)(numSegments - 1);
            calPsd2->data[i] = calPsd2->data[i] / (REAL4)(numSegments - 1);
          }
          else
          {
            calPsd1->data[i] = calPsd1->data[i] / (REAL4)numSegments;
            calPsd2->data[i] = calPsd2->data[i] / (REAL4)numSegments;
          }
          calInvPsd1.data->data[i] = 1. / calPsd1->data[i];
          calInvPsd2.data->data[i] = 1. / calPsd2->data[i];
        }

        /* print */
        if ((test_flag) && (interLoop == testInter) && (MCLoop == testTrial))
        {
          LALSPrintFrequencySeries(&calInvPsd1, "calInvPsdAvg1.dat");
          LALSPrintFrequencySeries(&calInvPsd2, "calInvPsdAvg2.dat");
        }

        if (verbose_flag)
          fprintf(stdout, "Normalising optimal filter...\n");

        /* compute variance and normalisation for optimal filter */
        LAL_CALL(LALStochasticOptimalFilterNormalization(&status, \
              &normOutput, &normInput, &normParams), &status);
        lambda = (REAL8)(normLambda.value * \
            pow(10.,normLambda.units.powerOfTen));
        varTheo = (REAL8)(segmentDuration * normSigma.value * \
            pow(10.,normSigma.units.powerOfTen));

        if (verbose_flag)
          fprintf(stdout, "Generating optimal filter...\n");

        /* build optimal filter */
        optFilter.epoch = gpsStartTime;
        LAL_CALL(LALStochasticOptimalFilterCal(&status, &optFilter, \
              &optFilterIn, &normLambda), &status);

        /* print */
        if ((test_flag) && (interLoop == testInter) && (MCLoop == testTrial))
          LALPrintFrequencySeries(&optFilter, "optFilter.dat");

        /* analyse middle segment */
        gpsStartTime.gpsSeconds = startTime + (interLoop + segMiddle) * \
                                  segmentDuration;
        segment1.epoch = segment2.epoch = gpsStartTime;

        if (verbose_flag)
        {
          fprintf(stdout, "analysing segment at GPS %d\n", \
              gpsStartTime.gpsSeconds);
        }

        for (i = 0; i < segmentLength; i++)
        {
          segment1.data->data[i] = seg1[segMiddle]->data[i];
          segment2.data->data[i] = seg2[segMiddle]->data[i];
        }

        /* print */
        if ((test_flag) && (interLoop == testInter) && (MCLoop == testTrial))
        {
          LALSPrintTimeSeries(&segment1, "segmentMiddle1.dat");
          LALSPrintTimeSeries(&segment2, "segmentMiddle2.dat");
        }

        /* zero pad and fft */
        LAL_CALL(LALSZeroPadAndFFT(&status, &hBarTilde1, &segment1, \
              &zeroPadParams), &status);
        LAL_CALL(LALSZeroPadAndFFT(&status, &hBarTilde2, &segment2, \
              &zeroPadParams), &status);

        /* print */
        if ((test_flag) && (interLoop == testInter) && (MCLoop == testTrial))
        {
          LALCPrintFrequencySeries(&hBarTilde1, "hBarTilde1.dat");
          LALCPrintFrequencySeries(&hBarTilde2, "hBarTilde2.dat");
        }

        if (verbose_flag)
          fprintf(stdout, "Generating cross correlation spectrum...\n");

        /* cc spectrum */
        for (i = 0; i < filterLength; i++)
        {
          response1.data->data[i] = resp1[segMiddle]->data[i];
          response2.data->data[i] = resp2[segMiddle]->data[i];
        }

        LAL_CALL(LALStochasticCrossCorrelationSpectrumCal(&status, \
              &ccSpectrum, &ccIn, epochsMatch), &status);

        /* save */
        if (verbose_flag)
          LALCPrintFrequencySeries(&ccSpectrum, "ccSpectrum.dat");

        /* cc statistic */
        LAL_CALL(LALStochasticCrossCorrelationStatisticCal(&status, &ccStat, \
              &ccIn,epochsMatch), &status);

        /* print */
        if ((test_flag) && (interLoop == testInter) && (MCLoop == testTrial))
          LALCPrintFrequencySeries(&ccSpectrum, "ccSpectrum.dat");

        y = (REAL8)(ccStat.value * pow(10., ccStat.units.powerOfTen));

        /* save */
        if (verbose_flag)
        {
          fprintf(stdout, "interval %d:", interLoop);
          fprintf(stdout, "GPS time = %d: y = %e, sigmaTheo = %e s, " \
              "varTheo = %e\n", gpsStartTime.gpsSeconds, y, \
              sqrt(varTheo), varTheo);
        }

        if (post_analysis_flag)
        {
          yOpt = yOpt + (y/varTheo);
          inVarTheoSum = inVarTheoSum + 1. / varTheo;
        }

        /* output to file */
        out1 = fopen(outputFilename1, "a");
        fprintf(out1,"%d %e %e %e\n", gpsStartTime.gpsSeconds, \
            y, sqrt(varTheo), varTheo);
        fclose(out1);
      }
    }

    lal_errhandler = LAL_ERR_EXIT;

    /* need to add post analysis for overlapping hann */
    if (post_analysis_flag)
    {
      ptEst = (yOpt / inVarTheoSum) /(REAL8)segmentDuration;
      error = sqrt (1./inVarTheoSum ) /(REAL8)segmentDuration;
      out2 = fopen(outputFilename2, "a");
      fprintf(out2,"%lld %lld %e %e\n", startTime, endTime, ptEst, error);
      fclose(out2);
      if (verbose_flag)
        fprintf(stdout, "ptEst = %e, error = %e\n", ptEst, error);
    }
  }

  /* cleanup */
  LAL_CALL(LALDestroyRealFFTPlan(&status, &(specparPSD.plan)), &status);
  LAL_CALL(LALDestroyRealFFTPlan(&status, &fftDataPlan), &status);
  LAL_CALL(LALDestroyVector(&status, &(segmentPad1.data)), &status);
  LAL_CALL(LALDestroyVector(&status, &(segmentPad2.data)), &status);
  LAL_CALL(LALDestroyVector(&status, &(segment1.data)), &status);
  LAL_CALL(LALDestroyVector(&status, &(segment2.data)), &status);
  LAL_CALL(LALDestroyVector(&status, &(psdTemp1.data)), &status);
  LAL_CALL(LALDestroyVector(&status, &(psdTemp2.data)), &status);
  LAL_CALL(LALDestroyVector(&status, &(psd1.data)), &status);
  LAL_CALL(LALDestroyVector(&status, &(psd2.data)), &status);
  LAL_CALL(LALDestroyVector(&status, &(calPsd1)), &status);
  LAL_CALL(LALDestroyVector(&status, &(calPsd2)), &status);
  LAL_CALL(LALCDestroyVector(&status, &(responseTemp1.data)), &status);
  LAL_CALL(LALCDestroyVector(&status, &(responseTemp2.data)), &status);
  LAL_CALL(LALCDestroyVector(&status, &(response1.data)), &status);
  LAL_CALL(LALCDestroyVector(&status, &(response2.data)), &status);
  LAL_CALL(LALDestroyVector(&status, &(optFilter.data)), &status);
  LAL_CALL(LALDestroyVector(&status, &(calInvPsd1.data)), &status);
  LAL_CALL(LALDestroyVector(&status, &(calInvPsd2.data)), &status);
  LAL_CALL(LALDestroyVector(&status, &(overlap.data)), &status);
  LAL_CALL(LALDestroyVector(&status, &(omegaGW.data)), &status);
  LAL_CALL(LALDestroyVector(&status, &(dataWindow.data)), &status);
  if (apply_mask_flag)
  {
    LAL_CALL(LALDestroyVector(&status, &(mask.data)), &status);
    LAL_CALL(LALDestroyVector(&status, &maskTemp), &status);
  }
  if (hannDuration != 0)
  {
    LAL_CALL(LALDestroyVector(&status, &hannWindow), &status );
  }
  LAL_CALL(LALCDestroyVector(&status, &(hBarTilde1.data)), &status);
  LAL_CALL(LALCDestroyVector(&status, &(hBarTilde2.data)), &status);
  if (inject_flag)
  {
    LAL_CALL(LALDestroyVector(&status, &SimStochBG1.data), &status);
    LAL_CALL(LALDestroyVector(&status, &SimStochBG2.data), &status);
    LAL_CALL(LALCDestroyVector(&status, &(MCresponse1.data)), &status);
    LAL_CALL(LALCDestroyVector(&status, &(MCresponse2.data)), &status);
    LAL_CALL(LALDestroyVector(&status, &(MComegaGW.data)), &status);
  }
  /*
  for (i = 0; i <numSegments; i++)
  {
    LAL_CALL(LALCDestroyVector(&status, &(resp1[i])), &status);
    LAL_CALL(LALCDestroyVector(&status, &(resp2[i])), &status);
    LAL_CALL(LALDestroyVector(&status, &(segPad1[i])), &status);
    LAL_CALL(LALDestroyVector(&status, &(segPad2[i])), &status);
    LAL_CALL(LALDestroyVector(&status, &(seg1[i])), &status);
    LAL_CALL(LALDestroyVector(&status, &(seg2[i])), &status);
    if (inject_flag)
    {
      LAL_CALL(LALCDestroyVector(&status, &(MCresp1[i])), &status);
      LAL_CALL(LALCDestroyVector(&status, &(MCresp2[i])), &status);
    }
  }
  */

  /* free calloc'd memory */
  free(frameCacheOne);
  free(frameCacheTwo);
  free(calCacheOne);
  free(calCacheTwo);
  free(channelOne);
  free(channelTwo);
  free(ifoOne);
  free(ifoTwo);

  return 0;
}

#define USAGE \
  "Usage: " PROGRAM_NAME " [options]\n"\
  " -h, --help                        print this message\n"\
  " -v, --version                     display version\n"\
  "     --verbose                     verbose mode\n"\
  " -z, --debug-level N               set lalDebugLevel\n"\
  " -t, --gps-start-time N            GPS start time\n"\
  " -T, --gps-end-time N              GPS end time\n"\
  " -L, --interval-duration N         interval duration\n"\
  " -l, --segment-duration N          segment duration\n"\
  " -A, --sample-rate N               sample rate\n"\
  " -a, --resample-rate N             resample rate\n"\
  " -f, --f-min N                     minimal frequency\n"\
  " -F, --f-max N                     maximal frequency\n"\
  " -i, --ifo-one IFO                 ifo for first stream\n"\
  " -I, --ifo-two IFO                 ifo for second stream\n"\
  " -z, --channel-one CHANNEL         channel for first stream\n"\
  " -Z, --channel-two CHANNEL         channel for second stream\n"\
  " -d, --frame-cache-one FILE        cache file for first stream\n"\
  " -D, --frame-cache-two FILE        cache file for second stream\n"\
  " -r, --calibration-cache-one FILE  first stream calibration cache\n"\
  " -R, --calibration-cache-two FILE  second stream calibration cache\n"\
  " -c, --calibration-offset N        calibration offset\n"\
  "     --apply-mask                  apply frequency masking\n"\
  " -b, --mask-bin N                  number of bins to mask\n"\
  "     --overlap-hann                overlaping hann windows\n"\
  " -w, --hann-duration N             hann duration\n"\
  "     --high-pass-filter            apply high pass filtering\n"\
  " -k, --hpf-frequency N             high pass filter knee frequency\n"\
  " -p, --hpf-attenuation N           high pass filter attenuation\n"\
  " -P, --hpf-order N                 high pass filter order\n"\
  "     --recenter                    recenter jobs\n"\
  "     --post-analysis               perform post analysis\n"\
  "     --middle-segment              use middle segment in PSD estimation\n"\
  "     --inject                      inject a signal into the data\n"\
  " -o, --scale-factor N              scale factor for injection\n"\
  " -g, --seed N                      random seed\n"\
  " -N, --trials N                    number of trials for Monte Carlo\n"\
  " -S, --output-dir DIR              directory for output files\n"\
  "     --test                        output intermediate results\n"\
  " -U, --test-interval N             interval to test\n"\
  " -V, --test-segment N              segment to test\n"\
  " -W, --test-trial N                trail to test\n"\
  " -K, --geo-hpf-frequency N         GEO high pass filter knee frequency\n"\
  " -q, --geo-hpf-attenuation N       GEO high pass filter attenuation\n"\
  " -Q, --geo-hpf-order N             GEO high pass filter order\n"

/* parse command line options */
void parseOptions(INT4 argc, CHAR *argv[])
{
  int c = -1;
  struct stat fileStatus;

  /* tempory variables */
  CHAR *channelOneTemp = NULL;
  CHAR *channelTwoTemp = NULL;

  while(1)
  {
    static struct option long_options[] =
    {
      /* options that set a flag */
      {"inject", no_argument, &inject_flag, 1},
      {"middle-segment", no_argument, &middle_segment_flag, 1},
      {"apply-mask", no_argument, &apply_mask_flag, 1},
      {"high-pass-filter", no_argument, &high_pass_flag, 1},
      {"overlap-hann", no_argument, &overlap_hann_flag, 1},
      {"post-analysis", no_argument, &post_analysis_flag,1},
      {"verbose", no_argument, &verbose_flag, 1},
      {"test", no_argument, &test_flag, 1},
      {"recenter", no_argument, &recenter_flag,1},
      /* options that don't set a flag */
      {"help", no_argument, 0, 'h'},
      {"gps-start-time", required_argument, 0, 't'},
      {"gps-end-time", required_argument, 0, 'T'},
      {"segment-big-duration", required_argument, 0, 'L'},
      {"segment-duration", required_argument, 0, 'l'},
      {"sample-rate", required_argument, 0, 'A'},
      {"resample-rate", required_argument, 0, 'a'},
      {"f-min", required_argument, 0, 'f'},
      {"f-max", required_argument, 0, 'F'},
      {"hann-duration", required_argument, 0, 'w'},
      {"hpf-frequency", required_argument, 0, 'k'},
      {"hpf-attenuation", required_argument, 0, 'p'},
      {"hpf-order", required_argument, 0, 'P'},
      {"geo-hpf-frequency", required_argument, 0, 'K'},
      {"geo-hpf-attenuation", required_argument, 0, 'q'},
      {"geo-hpf-order", required_argument, 0, 'Q'},
      {"ifo-one", required_argument, 0, 'i'},
      {"ifo-two", required_argument, 0, 'I'},
      {"channel-one", required_argument, 0, 'y'},
      {"channel-two", required_argument, 0, 'Y'},
      {"frame-cache-one", required_argument, 0, 'd'},
      {"frame-cache-two", required_argument, 0, 'D'},
      {"calibration-cache-one", required_argument, 0, 'r'},
      {"calibration-cache-two", required_argument, 0, 'R'},
      {"calibration-offset", required_argument, 0, 'c'},
      {"mask-bin", required_argument, 0, 'b'},
      {"scale-factor", required_argument, 0, 'o'},
      {"seed", required_argument, 0, 'g'},
      {"trials", required_argument, 0, 'N'},
      {"output-dir", required_argument, 0, 'S'},
      {"test-interval", required_argument, 0, 'U'},
      {"test-segment", required_argument, 0, 'V'},
      {"test-trial", required_argument, 0, 'W'},
      {"debug-level", required_argument, 0, 'z'},
      {"version", no_argument, 0, 'v'},
      {0, 0, 0, 0}
    };

    /* getopt_long stores the option here */
    int option_index = 0;
    size_t optarg_len;

    c = getopt_long(argc, argv, \
        "ht:T:L:l:A:a:f:F:w:k:p:P:K:q:Q:i:I:y:Y:d:D:r:R:c:b:o:g:N:S:U:V:W:z:v", \
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

      case 'L':
        /* interval duration */
        intervalDuration = atoi(optarg);

        /* check */
        if (intervalDuration <= 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Interval duration must be greater than 0: (%d specified)\n", \
              long_options[option_index].name, intervalDuration);
          exit(1);
        }

        break;

      case 'l':
        /* segment duration */
        segmentDuration = atoi(optarg);

        /* check */
        if (segmentDuration <= 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Segment duration must be greater than 0: (%d specified)\n", \
              long_options[option_index].name, segmentDuration);
          exit(1);
        }

        break;

      case 'A':
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

      case 'a':
        /* resampling */
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

      case 'k':
        /* high pass knee filter frequency  */
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
        /* high pass filter attenuation  */
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
        /* high pass filter order  */
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

      case 'K':
        /* GEO high pass knee filter frequency */
        geoHighPassFreq = atof(optarg);

        /* check */
        if (geoHighPassFreq < 0)
        {
          fprintf(stderr, "Invalid argument tp --%s:\n" \
              "GEO high pass filter knee frequency is less than 0 Hz: "\
              "(%f specified)\n", long_options[option_index].name, \
              geoHighPassFreq);
          exit(1);
        }

        break;

      case 'q':
        /*GEO high pass filter attenuation */
        geoHighPassAtten = atof(optarg);

        /* check */
        if ((geoHighPassAtten < 0.0) || (geoHighPassAtten > 1.0))
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "GEO high pass filter attenuation must be in the range [0:1]: " \
              "(%f specified)\n", long_options[option_index].name, \
              geoHighPassAtten);
          exit(1);
        }

        break;

      case 'Q':
        /* GEO high pass filter order */
        geoHighPassOrder = atoi(optarg);

        /* check */
        if (geoHighPassOrder <= 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "GEO high pass filter order must be greater than 0: " \
              "(%d specified)\n", long_options[option_index].name,
              geoHighPassOrder);
          exit(1);
        }

        break;

      case 'i':
        /* ifo for first stream */
        optarg_len = strlen(optarg) + 1;
        ifoOne = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        strncpy(ifoOne, optarg, optarg_len);

        /* set site */
        if (strncmp(ifoOne, "H1", 2) == 0)
        {
          siteOne = 0;
        }
        else if (strncmp(ifoOne, "H2", 2) == 0)
        {
          siteOne = 0;
        }
        else if (strncmp(ifoOne, "L1", 2) == 0)
        {
          siteOne = 1;
        }
        else if (strncmp(ifoOne, "G1", 2) == 0)
        {
          siteOne = 3;
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
        strncpy(ifoTwo, optarg, optarg_len);

        /* set site */
        if (strncmp(ifoTwo, "H1", 2) == 0)
        {
          siteTwo = 0;
        }
        else if (strncmp(ifoTwo, "H2", 2) == 0)
        {
          siteTwo = 0;
        }
        else if (strncmp(ifoTwo, "L1", 2) == 0)
        {
          siteTwo = 1;
        }
        else if (strncmp(ifoTwo, "G1", 2) == 0)
        {
          siteOne = 3;
        }
        else
        {
          fprintf(stderr, "Second IFO not recognised...\n");
          exit(1);
        }

        break;

      case 'y':
        /* channel one */
        optarg_len = strlen(optarg) + 4;
        channelOneTemp = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        channelOne = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        strncpy(channelOneTemp, optarg, optarg_len);
        break;

      case 'Y':
        /* channel two */
        optarg_len = strlen(optarg) + 4;
        channelTwoTemp = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        channelTwo = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        strncpy(channelTwoTemp, optarg, optarg_len);
        break;

      case 'd':
        /* data cache one */
        optarg_len = strlen(optarg) + 1;
        frameCacheOne = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        strncpy(frameCacheOne, optarg, optarg_len);

        /* check that file exists */
        if ( stat(frameCacheOne, &fileStatus) == -1)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "File does not exist: (%s specified)\n", \
              long_options[option_index].name, frameCacheOne);
          exit(1);
        }

        break;

      case 'D':
        /* data cache two */
        optarg_len = strlen(optarg) + 1;
        frameCacheTwo = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        strncpy(frameCacheTwo, optarg, optarg_len);

        /* check that file exists */
        if ( stat(frameCacheTwo, &fileStatus) == -1)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "File does not exist: (%s specified)\n", \
              long_options[option_index].name, frameCacheTwo);
          exit(1);
        }

        break;

      case 'r':
        /* calibration cache one */
        optarg_len = strlen(optarg) + 1;
        calCacheOne = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        strncpy(calCacheOne, optarg, optarg_len);

        /* check that file exists */
        if ( stat(calCacheOne, &fileStatus) == -1)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "File does not exist: (%s specified)\n", \
              long_options[option_index].name, calCacheOne);
          exit(1);
        }

        break;

      case 'R':
        /* calibration cache two */
        optarg_len = strlen(optarg) + 1;
        calCacheTwo = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        strncpy(calCacheTwo, optarg, optarg_len);

        /* check that file exists */
        if ( stat(calCacheTwo, &fileStatus) == -1)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "File does not exist: (%s specified)\n", \
              long_options[option_index].name, calCacheTwo);
          exit(1);
        }

        break;

      case 'c':
        /* calibration time offset */
        calibOffset = atoi(optarg);
        break;

      case 'b':
        /* number of bins to mask for frequency mask */
        maskBin = atoi(optarg);
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

      case 'S':
        /* directory for output files */
        strncpy(outputFilePath, optarg, LALNameLength);
        break;

      case 'U':
        /* interval number for test */
        testInter = atoi(optarg);
        break;

      case 'V':
        /* segment number for test */
        testSeg = atoi(optarg);
        break;

      case 'W':
        /* trial  number for test */
        testTrial = atoi(optarg);
        break;

      case 'z':
        /* set debug level */
        set_debug_level( optarg );
        break;

      case 'v':
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

  /* set channels */
  strcpy(channelOne, ifoOne);
  strcpy(channelTwo, ifoTwo);
  strcat(channelOne, ":");
  strcat(channelTwo, ":");
  strcat(channelOne, channelOneTemp);
  strcat(channelTwo, channelTwoTemp);
  free(channelOneTemp);
  free(channelTwoTemp);

  printf("channelOne = %s\n", channelOne);
  printf("channelTwo = %s\n", channelTwo);
  exit(1);

  return;
}

/* function to read data in frames */
void readDataPair(LALStatus *status,
    StreamPair *streamPair,
    ReadDataPairParams *params)
{
  /* counters */
  UINT4 j;
  INT4 i;

  /* variables */
  FrCache *frCache1 = NULL;
  FrStream *frStream1 = NULL;
  FrCache *frCache2 = NULL;
  FrStream *frStream2 = NULL;
  FrChanIn frChanIn1, frChanIn2;
  REAL8TimeSeries dataStreamGeo;
  REAL4TimeSeries dataStream1, dataStream2;
  ResampleTSParams resampleParams;
  LIGOTimeGPS bufferStartTime;
  UINT8 startTime;
  INT4 buffer;
  INT4 resampleRate, sampleRate;
  PassBandParamStruc geoHighpassParam;

  /* read parameters */
  startTime = params->startTime;
  buffer = params->buffer;
  resampleRate = params->resampleRate;
  sampleRate = params->sampleRate;

  /* initialise status pointer */
  INITSTATUS(status, "readDataPair", STOCHASTICC);
  ATTATCHSTATUSPTR(status);

  /* buffer start time */
  bufferStartTime.gpsSeconds = startTime - buffer;
  bufferStartTime.gpsNanoSeconds = 0;

  /* set channels */
  frChanIn1.name = params->channelOne;
  frChanIn2.name = params->channelTwo;
  frChanIn2.type = ADCDataChannel;
  frChanIn1.type = ADCDataChannel;

  /* initial data structures */
  dataStream1.epoch =  dataStream2.epoch  = bufferStartTime;

  if (verbose_flag)
    fprintf(stdout, "Allocating memory for raw data streams...\n");

  /* allocate memory */
  dataStream1.data = dataStream2.data = NULL;
  LALSCreateVector(status->statusPtr, &(dataStream1.data), \
      sampleRate * (params->duration + (2 * buffer)));
  CHECKSTATUSPTR(status);
  LALSCreateVector(status->statusPtr, &(dataStream2.data), \
      sampleRate * (params->duration + (2 * buffer)));
  CHECKSTATUSPTR(status);
  memset(dataStream1.data->data, 0, \
      dataStream1.data->length * sizeof(*dataStream1.data->data));
  memset(dataStream2.data->data, 0, \
      dataStream2.data->length * sizeof(*dataStream2.data->data));

  if ((strncmp(ifoOne, "G1", 2) == 0) || (strncmp(ifoTwo, "G1", 2) == 0))
  {
    dataStreamGeo.epoch = bufferStartTime;
    dataStreamGeo.data = NULL;
    LALDCreateVector(status->statusPtr, &(dataStreamGeo.data), \
        sampleRate * (params->duration + (2 * buffer)));
    CHECKSTATUSPTR(status);
    memset(dataStreamGeo.data->data, 0, \
        dataStreamGeo.data->length * sizeof(*dataStreamGeo.data->data));
  }

  if (verbose_flag)
    fprintf(stdout, "Opening first frame cache...\n");

  /* open first frame cache */
  LALFrCacheImport(status->statusPtr, &frCache1, params->frameCacheOne);
  CHECKSTATUSPTR(status);
  LALFrCacheOpen(status->statusPtr, &frStream1, frCache1);
  CHECKSTATUSPTR(status);

  if (verbose_flag)
    fprintf(stdout, "Reading in channel \"%s\"...\n", frChanIn1.name);

  /* read first channel */
  LALFrSeek(status->statusPtr, &(bufferStartTime), frStream1);
  CHECKSTATUSPTR(status);

  if (strncmp(ifoOne, "G1", 2) == 0)
  {
    LALFrGetREAL8TimeSeries(status->statusPtr, &dataStreamGeo, \
        &frChanIn1, frStream1);
    CHECKSTATUSPTR(status);

    /* high pass the GEO data */
    geoHighpassParam.nMax = geoHighPassOrder;
    geoHighpassParam.f1 = -1;
    geoHighpassParam.f2 = (REAL8)geoHighPassFreq;
    geoHighpassParam.a1 = -1;
    geoHighpassParam.a2 = geoHighPassAtten;
    LALButterworthREAL8TimeSeries(status->statusPtr, &dataStreamGeo, \
        &geoHighpassParam);
    CHECKSTATUSPTR(status);

    /* cast the GEO data to REAL4 in the channel time series */
    /* which already has the correct amount of memory allocated */
    for (j = 0; j < dataStream1.data->length; j++)
    {
      dataStream1.data->data[j] = geoScaleFactor * \
                                  (REAL4)(dataStreamGeo.data->data[j]);
    }

    /* re-copy the data paramaters from the GEO channel */
    /* to input data channel */
    LALSnprintf(dataStream1.name, LALNameLength * sizeof(CHAR), "%s", \
        dataStreamGeo.name);
    dataStream1.epoch = dataStreamGeo.epoch;
    dataStream1.deltaT = dataStreamGeo.deltaT;
    dataStream1.f0 = dataStreamGeo.f0;
    dataStream1.sampleUnits = dataStreamGeo.sampleUnits;

    /* free the REAL8 GEO input data */
    LALDDestroyVector(status->statusPtr, &dataStreamGeo.data);
    CHECKSTATUSPTR(status);
    dataStreamGeo.data = NULL;
  }
  else
  {
    LALFrGetREAL4TimeSeries(status->statusPtr, &dataStream1, \
        &frChanIn1, frStream1);
    CHECKSTATUSPTR(status);
  }
  if (strcmp(params->frameCacheOne, params->frameCacheTwo) == 0)
  {
    if (verbose_flag)
    {
      fprintf(stdout, "Reading in channel \"%s\" from same cache...\n", \
        frChanIn2.name);
    }

    /* read in second channel */
    LALFrSeek(status->statusPtr, &(bufferStartTime), frStream1);
    CHECKSTATUSPTR(status);
    LALFrGetREAL4TimeSeries(status->statusPtr, &dataStream2, &frChanIn2, \
        frStream1);
    CHECKSTATUSPTR(status);

    if (verbose_flag)
      fprintf(stdout, "Closing frame cache...\n");

    /* close frame cache */
    LALFrClose(status->statusPtr, &frStream1);
    CHECKSTATUSPTR(status);
  }
  else
  {
    if (verbose_flag)
      fprintf(stdout, "Closing first frame cache...\n");

    /* close first frame cache */
    LALFrClose(status->statusPtr, &frStream1);
    CHECKSTATUSPTR(status);

    if (verbose_flag)
      fprintf(stdout, "Opening second frame cache...\n");

    /* open second frame cache and read in second channel */
    LALFrCacheImport(status->statusPtr, &frCache2, params->frameCacheTwo);
    CHECKSTATUSPTR(status);
    LALFrCacheOpen(status->statusPtr, &frStream2, frCache2);
    CHECKSTATUSPTR(status);

    if (verbose_flag)
      fprintf(stdout, "Reading in channel \"%s\"...\n", frChanIn2.name);

    /* read in second channel */
    LALFrSeek(status->statusPtr, &(bufferStartTime), frStream2);
    CHECKSTATUSPTR(status);

    if (strncmp(ifoTwo, "G1", 2) == 0)
    {
      LALFrGetREAL8TimeSeries(status->statusPtr, &dataStreamGeo, \
          &frChanIn2, frStream2);
      CHECKSTATUSPTR(status);
      
      /* high pass the GEO data */
      geoHighpassParam.nMax = geoHighPassOrder;
      geoHighpassParam.f1 = -1;
      geoHighpassParam.f2 = (REAL8) geoHighPassFreq;
      geoHighpassParam.a1 = -1;
      geoHighpassParam.a2 = geoHighPassAtten;
      LALButterworthREAL8TimeSeries(status->statusPtr, &dataStreamGeo, \
          &geoHighpassParam);
      CHECKSTATUSPTR(status);

      /* cast the GEO data to REAL4 in the channel time series */
      /* which already has the correct amount of memory allocated */
      for (j = 0; j < dataStream2.data->length; j++)
      {
        dataStream2.data->data[j] = geoScaleFactor * \
                                    (REAL4)(dataStreamGeo.data->data[j]);
      }

      /* free the REAL8 GEO input data */
      LALDDestroyVector(status->statusPtr, &dataStreamGeo.data);
      CHECKSTATUSPTR(status);
      dataStreamGeo.data = NULL;
    }
    else
    {
      LALFrGetREAL4TimeSeries(status->statusPtr, &dataStream2, \
          &frChanIn2, frStream2);
      CHECKSTATUSPTR(status);
    }

    if (verbose_flag)
      fprintf(stdout, "Closing second frame cache...\n");

    /* close second frame stream */
    LALFrClose(status->statusPtr, &frStream2);
    CHECKSTATUSPTR(status);
  }

  /* resample */
  if (resampleRate != sampleRate)
  {
    if (verbose_flag)
      fprintf(stdout, "Resampling to %d Hz...\n", resampleRate);

    /* set resample parameters */
    resampleParams.deltaT = 1.0 / (REAL8)resampleRate;
    resampleParams.filterType = defaultButterworth;

    /* resample */
    LALResampleREAL4TimeSeries(status->statusPtr, &dataStream1, \
        &resampleParams);
    CHECKSTATUSPTR(status);
    LALResampleREAL4TimeSeries(status->statusPtr, &dataStream2, \
        &resampleParams);
    CHECKSTATUSPTR(status);
  }

  /* build output */
  strncpy(streamPair->streamOne->name,dataStream1.name, LALNameLength);
  strncpy(streamPair->streamTwo->name,dataStream2.name, LALNameLength);
  streamPair->streamOne->epoch.gpsSeconds = startTime;
  streamPair->streamTwo->epoch.gpsSeconds = startTime;
  streamPair->streamOne->epoch.gpsNanoSeconds = 0;
  streamPair->streamTwo->epoch.gpsNanoSeconds = 0;
  streamPair->streamOne->deltaT = 1./(REAL8)resampleRate;
  streamPair->streamTwo->deltaT = 1./(REAL8)resampleRate;
  streamPair->streamOne->f0 = streamPair->streamTwo->f0 = 0;
  streamPair->streamOne->sampleUnits = dataStream1.sampleUnits;
  streamPair->streamTwo->sampleUnits = dataStream2.sampleUnits;

  /* remove buffer, and hence corruption due to resampling */
  for (i = 0; i < params->duration * resampleRate; i++)
  {
    streamPair->streamOne->data->data[i] = \
      dataStream1.data->data[i + (resampleRate * buffer)];
    streamPair->streamTwo->data->data[i] = \
      dataStream2.data->data[i + (resampleRate * buffer)];
  }

  /* clean up */
  LALSDestroyVector(status->statusPtr, &(dataStream1.data));
  CHECKSTATUSPTR(status);
  LALSDestroyVector(status->statusPtr, &(dataStream2.data));
  CHECKSTATUSPTR(status);

  /* return status */
  DETATCHSTATUSPTR(status);
  RETURN(status);
}

/*
 * vim: et
 */
