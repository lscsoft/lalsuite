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
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include <errno.h>

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
#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOMetadataTables.h>

#include <lalapps.h>
#include <processtable.h>

#include "stochastic.h"

NRCSID(STOCHASTICC, "$Id$");
RCSID("$Id$");

/* cvs info */
#define PROGRAM_NAME "lalapps_stochastic"
#define CVS_ID "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_DATE "$Date$"
#define CVS_SOURCE "$Source$"

#define ADD_PROCESS_PARAM(pptype, format, ppvalue) \
this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
  calloc(1, sizeof(ProcessParamsTable)); \
  LALSnprintf(this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
      PROGRAM_NAME); \
  LALSnprintf(this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
      long_options[option_index].name); \
  LALSnprintf(this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype); \
  LALSnprintf(this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue);

/* system error checking */
extern int errno;

/* variables for getopt options parsing */
char *optarg;
int optind;

/* flags for getopt_long */
static int middle_segment_flag;
static int inject_flag;
static int apply_mask_flag;
static int high_pass_flag;
static int overlap_hann_flag;
static int recentre_flag;
extern int vrbflg;

/* xml comment/tags */
CHAR comment[LIGOMETA_COMMENT_MAX];
CHAR *userTag = NULL;

/* xml tables */
MetadataTable proctable;
MetadataTable procparams;
ProcessParamsTable *this_proc_param;
LIGOLwXMLStream xmlStream;

/* parameters for the stochastic search */

/* sampling parameters */
INT4 sampleRate = -1;
INT4 resampleRate = -1;
REAL8 deltaF = 0.25;

/* data parameters */
LIGOTimeGPS gpsStartTime;
INT4 startTime = 0;
INT4 endTime = 0;
INT4 intervalDuration = -1;
INT4 segmentDuration = -1;
INT4 calibOffset = -1;
CHAR *frameCacheOne = NULL;
CHAR *frameCacheTwo = NULL;
CHAR *calCacheOne = NULL;
CHAR *calCacheTwo = NULL;
CHAR *channelOne = NULL;
CHAR *channelTwo = NULL;
CHAR *ifoOne = NULL;
CHAR *ifoTwo = NULL;
INT4 siteOne;
INT4 siteTwo;

/* frequency band */
INT4 fMin = -1;
INT4 fMax = -1;

/* omegaGW parameters */
REAL4 alpha = 0;
REAL4 fRef = 100;
REAL4 omegaRef = 1;

/* monte carlo parameters */
REAL4 scaleFactor = -1;
INT4 seed = -1;
INT4 NLoop = 1;

/* window parameters */
INT4 hannDuration = -1;

/* high pass filtering parameters */
REAL4 highPassFreq = -1;
REAL4 highPassAtten = -1;
INT4  highPassOrder = -1;

/* GEO scale factor */
REAL4 geoScaleFactor = 1e18;

/* GEO high pass filter parameters */
REAL4 geoHighPassFreq = -1; /* 70; */
INT4  geoHighPassOrder = -1; /* 8; */
REAL4 geoHighPassAtten = -1; /* 0.9; */

/* number of bins for frequency masking */
INT4 maskBin = -1;

/* output file */
CHAR *outputFilePath = NULL;

INT4 main(INT4 argc, CHAR *argv[])
{
  /* lal initialisation variables */
  LALStatus status = blank_status;
  LALLeapSecAccuracy accuracy = LALLEAPSEC_LOOSE;

  /* output file */
  FILE *out;
  CHAR outputFilename[200];

  /* xml output file name */
  CHAR xmlFileName[FILENAME_MAX];

  /* counters */
  INT4 i, j, n, segLoop, interLoop, N;
  INT4 firstpass, pathology;

  /* results parameters */
  REAL8 y, yOpt;
  REAL8 varTheo;

  /* input data segment */
  INT4 duration, durationEff, extrasec;
  INT4 numSegments, numIntervals, segMiddle;
  INT4 segmentLength, segmentPadLength, intervalLength;
  INT4 segmentShift;
  INT4 padData;
  LIGOTimeGPS gpsStartPadTime, gpsCalibTime;
  ReadDataPairParams streamParams;
  StreamPair streamPair;
  REAL4TimeSeries *segmentOne;
  REAL4TimeSeries *segmentTwo;
  REAL4TimeSeries *segmentPadOne;
  REAL4TimeSeries *segmentPadTwo;
  REAL4Vector *segOne[100], *segPadOne[100], *segTwo[100], *segPadTwo[100];

  /* simulated signal structures */
  StochasticOmegaGWParameters parametersOmega;
  SSSimStochBGParams SBParams;
  SSSimStochBGInput SBInput;
  SSSimStochBGOutput SBOutput;
  REAL4TimeSeries *SimStochBGOne;
  REAL4TimeSeries *SimStochBGTwo;
  REAL4FrequencySeries MComegaGW;
  COMPLEX8FrequencySeries MCresponseOne, MCresponseTwo;
  COMPLEX8Vector *MCrespOne[100], *MCrespTwo[100];
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
  COMPLEX8FrequencySeries responseTempOne, responseTempTwo;
  COMPLEX8FrequencySeries responseOne, responseTwo;
  COMPLEX8Vector *respOne[100], *respTwo[100];
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
  REAL4FrequencySeries psdTempOne, psdTempTwo, psdOne, psdTwo;
  REAL4Vector *calPsdOne, *calPsdTwo;
  LALUnit psdUnits = {0,{0,0,1,0,0,0,2},{0,0,0,0,0,0,0}};

  /* calibrated inverse noise data structures */
  REAL4FrequencySeries calInvPsdOne, calInvPsdTwo;

  /* units for inverse noise */
  LALUnit calPSDUnit = {36,{0,0,-1,0,0,-2,0},{0,0,0,0,0,0,0}};

  /* structures for LALInverseNoise */
  StochasticInverseNoiseInput inverseNoiseInOne, inverseNoiseInTwo;
  StochasticInverseNoiseCalOutput inverseNoiseOutOne, inverseNoiseOutTwo;

  /* zeropad and fft structures */
  SZeroPadAndFFTParameters zeroPadParams;
  RealFFTPlan *fftDataPlan = NULL;
  COMPLEX8FrequencySeries hBarTildeOne, hBarTildeTwo;
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

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) calloc(1, sizeof(ProcessTable));
  LAL_CALL(LALGPSTimeNow(&status, &(proctable.processTable->start_time), \
        &accuracy), &status);
  LAL_CALL(populate_process_table(&status, proctable.processTable, \
        PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE), &status);
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *) \
                    calloc(1, sizeof(ProcessParamsTable));
  memset(comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR));

  /* parse command line options */
  parseOptions(argc, argv);

  /* get xml file name */
  LALSnprintf(xmlFileName, FILENAME_MAX, "%s%s-stochastic-%d-%d.xml", \
      ifoOne, ifoTwo, startTime, endTime);

  /* only add a buffer if the data is going to be resample and/or high
   * pass filtered */
  if ((sampleRate == resampleRate) && (high_pass_flag == 0))
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

  /* recentre */
  if (recentre_flag)
  {
    durationEff = (INT4)(((duration - (2 * padData)) / segmentDuration ) * \
                  segmentDuration);
    extrasec = duration - durationEff;
    startTime = startTime + (INT4)(extrasec / 2);
    endTime = startTime + durationEff;
  }

  if (overlap_hann_flag)
    segmentShift = segmentDuration / 2;

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

  if (vrbflg)
    fprintf(stdout, "Allocating memory for data segments...\n");

  /* allocate memory for data segments */
  LAL_CALL(LALCreateREAL4TimeSeries(&status, &segmentOne, "segmentOne", \
        gpsStartTime, 0, 1./resampleRate, lalDimensionlessUnit, \
        segmentLength), &status);
  LAL_CALL(LALCreateREAL4TimeSeries(&status, &segmentTwo, "segmentTwo", \
        gpsStartTime, 0, 1./resampleRate, lalDimensionlessUnit, \
        segmentLength), &status);
  LAL_CALL(LALCreateREAL4TimeSeries(&status, &segmentPadOne, "segmentPadOne", \
        gpsStartTime, 0, 1./resampleRate, lalDimensionlessUnit, \
        segmentPadLength), &status);
  LAL_CALL(LALCreateREAL4TimeSeries(&status, &segmentPadTwo, "segmentPadTwo", \
        gpsStartTime, 0, 1./resampleRate, lalDimensionlessUnit, \
        segmentPadLength), &status);
  
  for (i = 0; i < numSegments; i++)
  {
    segPadOne[i]= segPadTwo[i] = NULL;
    LAL_CALL(LALCreateVector(&status, &(segPadOne[i]), \
          segmentPadLength), &status);
    LAL_CALL(LALCreateVector(&status, &(segPadTwo[i]), \
          segmentPadLength), &status);
    memset(segPadOne[i]->data, 0, \
        segPadOne[i]->length * sizeof(*segPadOne[i]->data));
    memset(segPadTwo[i]->data, 0, \
        segPadTwo[i]->length * sizeof(*segPadTwo[i]->data));
  }

  for (i = 0; i < numSegments; i++)
  {
    segOne[i]= segTwo[i] = NULL;
    LAL_CALL(LALCreateVector(&status, &(segOne[i]), segmentLength), &status);
    LAL_CALL(LALCreateVector(&status, &(segTwo[i]), segmentLength), &status);
    memset(segOne[i]->data, 0, segOne[i]->length * sizeof(*segOne[i]->data));
    memset(segTwo[i]->data, 0, segTwo[i]->length * sizeof(*segTwo[i]->data));
  }

  /* set segment input parameters */
  streamParams.buffer = 0;
  streamParams.duration = segmentDuration + 2 * padData;
  streamPair.streamOne = segmentPadOne;
  streamPair.streamTwo = segmentPadTwo;

  if (inject_flag)
  {
    if (vrbflg)
      fprintf(stdout, "Allocating memory for MC...\n");

    MCdeltaT = 1.0 / resampleRate;
    MCdeltaF = (REAL8)resampleRate / (REAL8)segmentPadLength;
    MCfreqLength = (segmentPadLength / 2) + 1;

    /* create vectors to store the simulated signal */
    LAL_CALL(LALCreateREAL4TimeSeries(&status, &SimStochBGOne, \
          "Whitened-SimulatedSB1", gpsStartPadTime, 0, 1./resampleRate, \
          lalDimensionlessUnit, segmentPadLength), &status);
    LAL_CALL(LALCreateREAL4TimeSeries(&status, &SimStochBGTwo, \
          "Whitened-SimulatedSB2", gpsStartPadTime, 0, 1./resampleRate, \
          lalDimensionlessUnit, segmentPadLength), &status);
    
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
    strncpy(MCresponseOne.name,"MCresponseOne", LALNameLength);
    strncpy(MCresponseTwo.name,"MCresponseTwo", LALNameLength);
    MCresponseOne.sampleUnits = MCresponseTwo.sampleUnits = countPerStrain;
    MCresponseOne.epoch = MCresponseTwo.epoch = gpsCalibTime;
    MCresponseOne.deltaF = MCresponseTwo.deltaF = MCdeltaF;
    MCresponseOne.f0 = MCresponseTwo.f0 = 0;
    MCresponseOne.data = MCresponseTwo.data = NULL;
    LAL_CALL(LALCCreateVector(&status, &(MCresponseOne.data), \
          MCfreqLength), &status);
    LAL_CALL(LALCCreateVector(&status, &(MCresponseTwo.data), \
          MCfreqLength), &status);
    memset(MCresponseOne.data->data, 0, \
        MCresponseOne.data->length * sizeof(*MCresponseOne.data->data));
    memset(MCresponseTwo.data->data, 0, \
        MCresponseTwo.data->length * sizeof(*MCresponseTwo.data->data));

    for (i = 0; i < numSegments; i++)
    {
      MCrespOne[i]= MCrespTwo[i] = NULL;
      LAL_CALL(LALCCreateVector(&status, &(MCrespOne[i]), MCfreqLength), \
          &status);
      LAL_CALL(LALCCreateVector(&status, &(MCrespTwo[i]), MCfreqLength), \
          &status);
      memset(MCrespOne[i]->data, 0, \
          MCrespOne[i]->length * sizeof(*MCrespOne[i]->data));
      memset(MCrespTwo[i]->data, 0, \
          MCrespTwo[i]->length * sizeof(*MCrespTwo[i]->data));
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
  strncpy(psdTempOne.name, "psdTempOne", LALNameLength);
  strncpy(psdTempTwo.name, "psdTempTwo", LALNameLength);
  psdTempOne.sampleUnits = psdTempTwo.sampleUnits = psdUnits;
  psdTempOne.deltaF = psdTempTwo.deltaF = deltaF;
  psdTempOne.f0 = psdTempTwo.f0 = 0;

  if (vrbflg)
    fprintf(stdout, "Allocating memory for PSDs...\n");

  /* allocate memory for PSDs */
  psdTempOne.data = psdTempTwo.data = NULL;
  LAL_CALL(LALCreateVector(&status, &(psdTempOne.data), psdTempLength), \
      &status);
  LAL_CALL(LALCreateVector(&status, &(psdTempTwo.data), psdTempLength), \
      &status);
  memset(psdTempOne.data->data, 0, \
      psdTempOne.data->length * sizeof(*psdTempOne.data->data));
  memset(psdTempTwo.data->data, 0, \
      psdTempTwo.data->length * sizeof(*psdTempTwo.data->data));

  /* reduced frequency band PSDs */
  /* set metadata fields for reduced frequency band PSDs */
  strncpy(psdOne.name, "psdOne", LALNameLength);
  strncpy(psdTwo.name, "psdTwo", LALNameLength);
  psdOne.deltaF = psdTwo.deltaF = deltaF;
  psdOne.f0 = psdTwo.f0 = fMin;
  psdOne.sampleUnits = psdTwo.sampleUnits = psdUnits;

  if (vrbflg)
    fprintf(stdout, "Allocating memory for reduced frequency band PSDs...\n");
  
  /* allocate memory for reduced frequency band PSDs */
  psdOne.data = psdTwo.data = NULL;
  LAL_CALL(LALCreateVector(&status, &psdOne.data, filterLength), &status);
  LAL_CALL(LALCreateVector(&status, &psdTwo.data, filterLength), &status);
  memset(psdOne.data->data, 0, psdOne.data->length * \
      sizeof(*psdOne.data->data));
  memset(psdTwo.data->data, 0, psdTwo.data->length * \
      sizeof(*psdTwo.data->data));

  /* allocate memory for calibrated PSDs */
  calPsdOne = calPsdTwo = NULL;
  LAL_CALL(LALCreateVector(&status, &calPsdOne, filterLength), &status);
  LAL_CALL(LALCreateVector(&status, &calPsdTwo, filterLength), &status);
  memset(calPsdOne->data, 0, calPsdOne->length * sizeof(*calPsdOne->data));
  memset(calPsdTwo->data, 0, calPsdTwo->length * sizeof(*calPsdTwo->data));

  /* set parameters for response functions */
  respLength = (UINT4)(fMax / deltaF) + 1;

  /* set metadata fields for response functions */
  strncpy(responseTempOne.name, "responseTempOne", LALNameLength);
  strncpy(responseTempTwo.name, "responseTempTwo", LALNameLength);
  responseTempOne.sampleUnits = countPerAttoStrain;
  responseTempTwo.sampleUnits = countPerAttoStrain;
  responseTempOne.epoch = responseTempTwo.epoch = gpsCalibTime;
  responseTempOne.deltaF = responseTempTwo.deltaF = deltaF;
  responseTempOne.f0 = responseTempTwo.f0 = 0;

  if (vrbflg)
    fprintf(stdout, "Allocating memory for response functions...\n");

  /* allocate memory for response functions */
  responseTempOne.data = responseTempTwo.data = NULL;
  LAL_CALL(LALCCreateVector(&status, &(responseTempOne.data), respLength), \
      &status);
  LAL_CALL(LALCCreateVector(&status, &(responseTempTwo.data), respLength), \
      &status);
  memset(responseTempOne.data->data, 0, \
      responseTempOne.data->length * sizeof(*responseTempOne.data->data));
  memset(responseTempTwo.data->data, 0, \
      responseTempTwo.data->length * sizeof(*responseTempTwo.data->data));

  /* reduced frequency band response functions */
  /* set metadata fields for reduced frequency band response functions */
  strncpy(responseOne.name, "responseOne", LALNameLength);
  strncpy(responseTwo.name, "responseTwo", LALNameLength);
  responseOne.sampleUnits = responseTwo.sampleUnits = countPerAttoStrain;
  responseOne.epoch = responseTwo.epoch = gpsCalibTime;
  responseOne.deltaF = responseTwo.deltaF = deltaF;
  responseOne.f0 = responseTwo.f0 = fMin;

  if (vrbflg)
  {
    fprintf(stdout, "Allocating memory for reduced frequency band " \
        "response functions...\n");
  }

  /* allocate memory for reduced frequency band response functions */
  responseOne.data = responseTwo.data = NULL;
  LAL_CALL(LALCCreateVector(&status, &(responseOne.data), filterLength), \
      &status);
  LAL_CALL(LALCCreateVector(&status, &(responseTwo.data), filterLength), \
      &status);
  memset(responseOne.data->data, 0, \
      responseOne.data->length * sizeof(*responseOne.data->data));
  memset(responseTwo.data->data, 0, \
      responseTwo.data->length * sizeof(*responseTwo.data->data));

  for (i = 0; i < numSegments; i++)
  {
    respOne[i]= respTwo[i] = NULL;
    LAL_CALL(LALCCreateVector(&status, &(respOne[i]), filterLength), &status);
    LAL_CALL(LALCCreateVector(&status, &(respTwo[i]), filterLength), &status);
    memset(respOne[i]->data, 0, respOne[i]->length * \
        sizeof(*respOne[i]->data));
    memset(respTwo[i]->data, 0, respTwo[i]->length * \
        sizeof(*respTwo[i]->data));
  }

  if (vrbflg)
    fprintf(stdout, "Creating FFT plan for PSD estimation...\n");

  /* create fft plan */
  LAL_CALL(LALCreateForwardRealFFTPlan(&status, &specparPSD.plan, \
        windowPSDLength, 0), &status);

  if (vrbflg)
    fprintf(stdout, "Creating window for PSD estimation...\n");

  /* create window for PSD estimation */
  LAL_CALL(LALCreateREAL4Window(&status, &specparPSD.window, &winparPSD), \
      &status);

  /* set metadata fields for inverse noise structures */
  strncpy(calInvPsdOne.name, "calInvPsdOne", LALNameLength);
  strncpy(calInvPsdTwo.name, "calInvPsdTwo", LALNameLength);
  calInvPsdOne.sampleUnits = calInvPsdTwo.sampleUnits = calPSDUnit;
  calInvPsdOne.deltaF = calInvPsdTwo.deltaF = deltaF;
  calInvPsdOne.f0 = calInvPsdTwo.f0 = fMin;

  if (vrbflg)
    fprintf(stdout, "Allocating memory for inverse noise...\n");

  /* allocate memory for inverse noise */
  calInvPsdOne.data = calInvPsdTwo.data = NULL;
  LAL_CALL(LALCreateVector(&status, &(calInvPsdOne.data), filterLength), \
      &status);
  LAL_CALL(LALCreateVector(&status, &(calInvPsdTwo.data), filterLength), \
      &status);
  memset(calInvPsdOne.data->data, 0, \
      calInvPsdOne.data->length * sizeof(*calInvPsdOne.data->data));
  memset(calInvPsdTwo.data->data, 0, \
      calInvPsdTwo.data->length * sizeof(*calInvPsdTwo.data->data));

  /* set inverse noise inputs */
  inverseNoiseInOne.unCalibratedNoisePSD = &psdOne;
  inverseNoiseInOne.responseFunction = &responseOne;
  inverseNoiseInTwo.unCalibratedNoisePSD = &psdTwo;
  inverseNoiseInTwo.responseFunction = &responseTwo;

  /* set inverse noise outputs */
  inverseNoiseOutOne.calibratedInverseNoisePSD = &calInvPsdOne;
  inverseNoiseOutTwo.calibratedInverseNoisePSD = &calInvPsdTwo;

  /* set window parameters for segment data streams */
  strncpy(dataWindow.name, "dataWindow", LALNameLength);
  dataWindow.sampleUnits = lalDimensionlessUnit;
  dataWindow.deltaT = 1./resampleRate;
  dataWindow.f0 = 0;

  if (vrbflg)
    fprintf(stdout, "Allocating memory for data segment window...\n");

  /* allocate memory for segment window */
  dataWindow.data = NULL;
  LAL_CALL(LALSCreateVector(&status, &(dataWindow.data), segmentLength), \
      &status);
  memset(dataWindow.data->data, 0, \
      dataWindow.data->length * sizeof(*dataWindow.data->data));

  if (vrbflg)
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
  strncpy(hBarTildeOne.name, "hBarTildeOne", LALNameLength);
  strncpy(hBarTildeTwo.name, "hBarTildeTwo", LALNameLength);

  if (vrbflg)
    fprintf(stdout, "Allocating memory for zeropad...\n");

  /* allocate memory for zeropad */
  hBarTildeOne.data = hBarTildeTwo.data = NULL;
  LAL_CALL(LALCCreateVector(&status, &(hBarTildeOne.data), fftDataLength), \
      &status);
  LAL_CALL(LALCCreateVector(&status, &(hBarTildeTwo.data), fftDataLength), \
      &status);
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

  if (vrbflg)
    fprintf(stdout, "Generating the overlap reduction function...\n");

  /* generate overlap reduction function */
  LAL_CALL(LALOverlapReductionFunction(&status, &overlap, &detectors, \
        &ORFparams), &status);

  /* set metadata fields for spectrum */
  strncpy(omegaGW.name, "omegaGW", LALNameLength);
  omegaGW.sampleUnits = lalDimensionlessUnit;
  omegaGW.deltaF = deltaF;
  omegaGW.f0 = fMin;

  if (vrbflg)
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

  if (vrbflg)
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

    if (vrbflg)
      fprintf(stdout, "Allocating memory for frequency mask...\n");

    /* allocate memory for frequency mask */
    mask.data = maskTemp = NULL;
    LAL_CALL(LALCreateVector(&status, &(mask.data), filterLength), &status);
    LAL_CALL(LALCreateVector(&status, &maskTemp, respLength), &status);
    memset(mask.data->data, 0, mask.data->length * sizeof(*mask.data->data));
    memset(maskTemp->data, 0, maskTemp->length * sizeof(*maskTemp->data));

    if (vrbflg)
      fprintf(stdout, "Generating frequency mask...\n");

    /* set all values to 1 */
    for (i = 0; i < respLength; i++)
      maskTemp->data[i] = 1.;

    if (vrbflg)
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

    if (vrbflg)
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

    if (vrbflg)
      fprintf(stdout, "Getting appropriate frequency band for mask...\n");

    /* get appropriate band */
    for (i = 0; i < filterLength; i++)
      mask.data->data[i] = maskTemp->data[i + numFMin];

    if (vrbflg)
      fprintf(stdout, "Applying frequency mask to spectrum..\n");

    /* apply mask to omegaGW */
    for (i = 0; i < filterLength; i++)
      omegaGW.data->data[i] *= mask.data->data[i];
  }

  /* set normalisation parameters */
  normParams.fRef = fRef;
  normParams.heterodyned = 0;
  normParams.window1 = normParams.window2 = dataWindow.data;

  /* set normalisation input */
  normInput.overlapReductionFunction = &overlap;
  normInput.omegaGW = &omegaGW;
  normInput.inverseNoisePSD1 = &calInvPsdOne;
  normInput.inverseNoisePSD2 = &calInvPsdTwo;

  /* set normalisation output */
  normOutput.normalization = &normLambda;
  normOutput.variance = &normSigma;

  /* set metadata fields for optimal filter */
  strncpy(optFilter.name, "optFilter", LALNameLength);
  optFilter.epoch = gpsStartTime;
  optFilter.deltaF = deltaF;
  optFilter.f0 = fMin;

  if (vrbflg)
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
  optFilterIn.calibratedInverseNoisePSD1 = &calInvPsdOne;
  optFilterIn.calibratedInverseNoisePSD2 = &calInvPsdTwo;

  /* set metadata fields for CC spectrum */
  strncpy(ccSpectrum.name, "ccSpectrum", LALNameLength);
  ccSpectrum.epoch = gpsStartTime;
  ccSpectrum.deltaF = deltaF;
  ccSpectrum.f0 = fMin;

  if (vrbflg)
    fprintf(stdout, "Allocating memory for CC Spectrum...\n");

  /* allocate memory for CC spectrum*/
  ccSpectrum.data = NULL;
  LAL_CALL(LALCCreateVector(&status, &(ccSpectrum.data), filterLength), \
      &status);
  memset(ccSpectrum.data->data, 0, \
      ccSpectrum.data->length * sizeof(*ccSpectrum.data->data));

  /* set CC inputs */
  ccIn.hBarTildeOne = &hBarTildeOne;
  ccIn.hBarTildeTwo = &hBarTildeTwo;
  ccIn.responseFunctionOne = &responseOne;
  ccIn.responseFunctionTwo = &responseTwo;
  ccIn.optimalFilter = &optFilter;

  if (vrbflg)
    fprintf(stdout, "Done with memory allocation...\n");

  if (overlap_hann_flag)
    N = 2;
  else
    N = 1;

  for (MCLoop = 0; MCLoop < NLoop; MCLoop++)
  {
    /* initialize parameters for post analysis */
    yOpt = 0;

    /* open output file */
    LALSnprintf(outputFilename, LALNameLength, \
        "%s/stat-%s%s-%lld-%lld-%d.dat", outputFilePath, ifoOne, ifoTwo, \
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
          if (vrbflg)
          {
            fprintf(stdout, "First Pass\ninterval %d out of %d\n", \
                interLoop + 1, numIntervals);
          }

          /* read first interval and get response functions */
          lal_errhandler = LAL_ERR_RTRN;

          for (segLoop = 0; segLoop < numSegments; segLoop++)
          {
            /* define segment epoch */
            gpsStartTime.gpsSeconds = startTime + (interLoop + segLoop) * \
                                      segmentDuration;
            gpsStartPadTime.gpsSeconds = gpsStartTime.gpsSeconds - padData;
            segmentPadOne->epoch = gpsStartPadTime;
            segmentPadTwo->epoch = gpsStartPadTime;

            if (vrbflg)
            {
              fprintf(stdout, "request data at GPS time %d\n", \
                  gpsStartTime.gpsSeconds);
            }

            streamParams.start = gpsStartPadTime.gpsSeconds;
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
              segPadOne[segLoop]->data[i] = segmentPadOne->data->data[i];
              segPadTwo[segLoop]->data[i] = segmentPadTwo->data->data[i];
            }

            /* compute response functions */
            gpsCalibTime.gpsSeconds = gpsStartTime.gpsSeconds  + calibOffset;
            responseOne.epoch  = responseTwo.epoch  = gpsCalibTime;

            /* set response function to unity for GEO */
            if (strncmp(ifoOne, "G1", 2) == 0)
            {
              for (i = 0; i < filterLength; i++)
              {
                responseOne.data->data[i].re = 1;
                responseOne.data->data[i].im = 0;
              }
            }
            else
            {
              responseTempOne.epoch = gpsCalibTime;
              if (vrbflg)
              {
                fprintf(stdout, "request GPS time %d for ifo 1\n", \
                    gpsCalibTime.gpsSeconds );
              }
              memset(&calfacts, 0, sizeof(CalibrationUpdateParams));
              calfacts.ifo = ifoOne;
              LAL_CALL(LALFrCacheImport(&status, &calibCache, calCacheOne), \
                  &status);
              LAL_CALL(LALExtractFrameResponse(&status, &responseTempOne, \
                    calibCache, &calfacts), &status);
              LAL_CALL(LALDestroyFrCache(&status, &calibCache), &status);

              /* skip segment if data not found or corrupted with 0 values */
              if ((status.statusCode !=0)||(responseTempOne.data==NULL))
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
                responseOne.data->data[i] = responseTempOne.data->data[i + \
                                          numFMin];
              }
            }

            /* set response function to unity for GEO */
            if (strncmp(ifoTwo, "G1", 2) == 0)
            {
              for (i = 0; i < filterLength; i++)
              {
                responseTwo.data->data[i].re = 1;
                responseTwo.data->data[i].im = 0;
              }
            }
            else
            {
              responseTempTwo.epoch = gpsCalibTime;
              if (vrbflg)
              {
                fprintf(stdout, "request GPS time %d for ifo 2\n", \
                    gpsCalibTime.gpsSeconds );
              }
              memset(&calfacts, 0, sizeof(CalibrationUpdateParams));
              calfacts.ifo = ifoTwo;
              LAL_CALL(LALFrCacheImport(&status, &calibCache, calCacheTwo ), \
                  &status);
              LAL_CALL(LALExtractFrameResponse( &status, &responseTempTwo, \
                    calibCache, &calfacts), &status);
              LAL_CALL(LALDestroyFrCache( &status, &calibCache), &status);

              /* skip segment if data not found or corrupted with 0 values */
              if ((status.statusCode !=0)||(responseTempTwo.data==NULL))
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
                responseTwo.data->data[i] = responseTempTwo.data->data[i + \
                                          numFMin];
              }
            }

            /* store in memory */
            for (i = 0; i < filterLength; i++)
            {
              respOne[segLoop]->data[i] = responseOne.data->data[i];
              respTwo[segLoop]->data[i] = responseTwo.data->data[i];
            }

            /* convert response function for use in the MC routine */
            if (inject_flag)
            {
              MCresponseOne.epoch = MCresponseTwo.epoch = gpsCalibTime;
              LAL_CALL(LALResponseConvert(&status, &MCresponseOne, \
                    &responseTempOne), &status);
              LAL_CALL(LALResponseConvert(&status, &MCresponseTwo, \
                    &responseTempTwo), &status);

              /* force DC to be 0 and nyquist to be real */
              MCresponseOne.data->data[0].re = 0;
              MCresponseTwo.data->data[0].re = 0;
              MCresponseOne.data->data[0].im = 0;
              MCresponseTwo.data->data[0].im = 0;
              MCresponseOne.data->data[MCfreqLength-1].im = 0;
              MCresponseTwo.data->data[MCfreqLength-1].im = 0;

              /* store in memory */
              for (i = 0; i < MCfreqLength; i++)
              {
                MCrespOne[segLoop]->data[i] = MCresponseOne.data->data[i];
                MCrespTwo[segLoop]->data[i] = MCresponseTwo.data->data[i];
              }
            }
          }

          if (pathology == 1)
          {
            if (vrbflg)
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
          if (vrbflg)
          {
            fprintf(stdout, "Shift Segment\ninterval %d out of %d\n", \
                interLoop, numIntervals);
          }
          /* shift segments */
          for (segLoop = 0; segLoop < numSegments - 1; segLoop++)
          {
            for (i = 0; i < segmentPadLength; i++)
            {
              segPadOne[segLoop]->data[i] = segPadOne[segLoop + 1]->data[i];
              segPadTwo[segLoop]->data[i] = segPadTwo[segLoop + 1]->data[i];
            }
            for (i = 0; i < filterLength; i++)
            {
              respOne[segLoop]->data[i] = respOne[segLoop + 1]->data[i];
              respTwo[segLoop]->data[i] = respTwo[segLoop + 1]->data[i];

              if (inject_flag)
              {
                MCrespOne[segLoop]->data[i] = MCrespOne[segLoop + 1]->data[i];
                MCrespTwo[segLoop]->data[i] = MCrespTwo[segLoop + 1]->data[i];
              }
            }
          }

          /* read extra segment */
          gpsStartTime.gpsSeconds = startTime + (interLoop + \
              numSegments - 1) * segmentDuration;
          gpsStartPadTime.gpsSeconds = gpsStartTime.gpsSeconds - padData;
          streamParams.start = gpsStartPadTime.gpsSeconds;

          if (vrbflg)
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
            segPadOne[numSegments-1]->data[i] = segmentPadOne->data->data[i];
            segPadTwo[numSegments-1]->data[i] = segmentPadTwo->data->data[i];
          }

          /* compute extra response function */
          gpsCalibTime.gpsSeconds = gpsStartTime.gpsSeconds + calibOffset;
          responseOne.epoch  = responseTwo.epoch = gpsCalibTime;

          if (strncmp(ifoOne, "G1", 2) == 0)
          {
            for (i = 0; i < filterLength; i++)
            {
              responseOne.data->data[i].re = 1;
              responseOne.data->data[i].im = 0;
            }
          }
          else
          {
            responseTempOne.epoch = gpsCalibTime;
            if (vrbflg)
            {
              fprintf(stdout, "request GPS time %d for ifo 1\n", \
                  gpsCalibTime.gpsSeconds );
            }
            memset(&calfacts, 0, sizeof(CalibrationUpdateParams));
            calfacts.ifo = ifoOne;
            LAL_CALL(LALFrCacheImport(&status, &calibCache, calCacheOne ), \
                &status);
            LAL_CALL(LALExtractFrameResponse(&status, &responseTempOne, \
                  calibCache, &calfacts), &status);
            LAL_CALL(LALDestroyFrCache(&status, &calibCache), &status);

            /* skip segment if data not found or corrupted with 0 values */
            if ((status.statusCode != 0) || (responseTempOne.data == NULL))
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
              responseOne.data->data[i] = responseTempOne.data->data[i + \
                                        numFMin];
            }
          }

          if (strncmp(ifoTwo, "G1", 2) == 0)
          {
            for (i = 0; i < filterLength; i++)
            {
              responseOne.data->data[i].re = 1;
              responseOne.data->data[i].im = 0;
            }
          }
          else
          {
            responseTempTwo.epoch = gpsCalibTime;
            if (vrbflg)
            {
              fprintf(stdout, "request GPS time %d for ifo 2\n", \
                  gpsCalibTime.gpsSeconds );
            }
            memset(&calfacts, 0, sizeof(CalibrationUpdateParams));
            calfacts.ifo = ifoTwo;
            LAL_CALL(LALFrCacheImport(&status, &calibCache, calCacheTwo ), \
                &status);
            LAL_CALL(LALExtractFrameResponse(&status, &responseTempTwo, \
                  calibCache, &calfacts), &status);
            LAL_CALL(LALDestroyFrCache(&status, &calibCache), &status);

            /* skip segment if data not found or corrupted with 0 values */
            if ((status.statusCode != 0) || (responseTempTwo.data == NULL))
            {
              clear_status(&status);
              firstpass = 1;
              pathology = 1;
              interLoop = interLoop + segLoop;
              break;
            }

            /* reduce to the optimal filter frequency range */
            responseTwo.epoch  = gpsCalibTime;
            for (i = 0; i < filterLength; i++)
            {
              responseTwo.data->data[i] = responseTempTwo.data->data[i + \
                                        numFMin];
            }
          }

          /* store in memory */
          for (i = 0; i < filterLength; i++)
          {
            respOne[numSegments-1]->data[i] = responseOne.data->data[i];
            respTwo[numSegments-1]->data[i] = responseTwo.data->data[i];
          }

          /* convert response function for use in the MC routine */
          if (inject_flag)
          {
            MCresponseOne.epoch = MCresponseTwo.epoch = gpsCalibTime;
            LAL_CALL(LALResponseConvert(&status, &MCresponseOne, \
                  &responseTempOne), &status);
            LAL_CALL(LALResponseConvert(&status, &MCresponseTwo, \
                  &responseTempTwo), &status);

            /* force DC to be 0 and nyquist to be real */
            MCresponseOne.data->data[0].re = 0;
            MCresponseTwo.data->data[0].re = 0;
            MCresponseOne.data->data[0].im = 0;
            MCresponseTwo.data->data[0].im = 0;
            MCresponseOne.data->data[MCfreqLength-1].im = 0;
            MCresponseTwo.data->data[MCfreqLength-1].im = 0;

            /* store in memory */
            for (i = 0; i < MCfreqLength; i++)
            {
              MCrespOne[numSegments-1]->data[i] = MCresponseOne.data->data[i];
              MCrespTwo[numSegments-1]->data[i] = MCresponseTwo.data->data[i];
            }
          }
        }

        /* initialize average PSDs */
        for (i = 0; i < filterLength; i++)
        {
          calPsdOne->data[i] = 0;
          calPsdTwo->data[i] = 0;
        }

        for (segLoop = 0; segLoop < numSegments; segLoop++)
        {
          gpsStartTime.gpsSeconds = startTime + (interLoop + segLoop) * \
                                    segmentDuration;
          gpsStartPadTime.gpsSeconds = gpsStartTime.gpsSeconds - padData;
          segmentPadOne->epoch = gpsStartPadTime;
          segmentPadTwo->epoch = gpsStartPadTime;
          segmentOne->epoch = gpsStartTime;
          segmentTwo->epoch = gpsStartTime;
          gpsCalibTime.gpsSeconds = gpsStartTime.gpsSeconds  + calibOffset;
          responseOne.epoch = responseTwo.epoch = gpsCalibTime;

          for (i = 0; i < filterLength; i++)
          {
            responseOne.data->data[i] = respOne[segLoop]->data[i];
            responseTwo.data->data[i] = respTwo[segLoop]->data[i];
          }

          /* simulate signal */
          if (inject_flag)
          {
            for (i = 0; i < MCfreqLength; i++)
            {
              MCresponseOne.data->data[i] = MCrespOne[segLoop]->data[i];
              MCresponseTwo.data->data[i] = MCrespTwo[segLoop]->data[i];
            }

            /* set parameters for monte carlo */
            SimStochBGOne->epoch = gpsStartPadTime;
            SimStochBGTwo->epoch = gpsStartPadTime;
            SBParams.seed = seed;

            /* define input structure for SimulateSB */
            SBInput.omegaGW = &MComegaGW;
            SBInput.whiteningFilter1 = &MCresponseOne;
            SBInput.whiteningFilter2 = &MCresponseTwo;

            /* define output structure for SimulateSB */
            SBOutput.SSimStochBG1 = SimStochBGOne;
            SBOutput.SSimStochBG2 = SimStochBGTwo;

            /* perform monte carlo */
            LAL_CALL(LALSSSimStochBGTimeSeries(&status, &SBOutput, \
                  &SBInput, &SBParams), &status);

            /* multiply by scale factor and inject into real data */
            for (i = 0; i < segmentPadLength ; i++)
            {
              segmentPadOne->data->data[i] = segPadOne[segLoop]->data[i] + \
                (scaleFactor * SimStochBGOne->data->data[i]);
              segmentPadTwo->data->data[i] = segPadTwo[segLoop]->data[i] + \
                (scaleFactor * SimStochBGTwo->data->data[i]);
            }

            /* increase seed */
            seed = seed + 2;
          }
          else
          {
            for (i = 0; i < segmentPadLength; i++)
            {
              segmentPadOne->data->data[i] = segPadOne[segLoop]->data[i];
              segmentPadTwo->data->data[i] = segPadTwo[segLoop]->data[i];
            }
          }

          /* high pass fitering */
          if (high_pass_flag)
          {
            LAL_CALL(LALButterworthREAL4TimeSeries(&status, segmentPadOne, \
                  &highpassParam), &status);
            LAL_CALL(LALButterworthREAL4TimeSeries(&status, segmentPadTwo, \
                  &highpassParam), &status);
          }

          /* throw away pad data on each side of the segment */
          for (i = 0; i < segmentLength; i++)
          {
            segmentOne->data->data[i] = segmentPadOne->data->data[i + \
                                        padData * resampleRate];
            segmentTwo->data->data[i] = segmentPadTwo->data->data[i + \
                                        padData * resampleRate];
          }

          /* store in memory */
          for (i = 0; i < segmentLength; i++)
          {
            segOne[segLoop]->data[i] = segmentOne->data->data[i];
            segTwo[segLoop]->data[i] = segmentTwo->data->data[i];
          }

          if ((middle_segment_flag == 0) && (segLoop == segMiddle))
          {
            if (vrbflg)
              fprintf(stdout, "ignoring middle segment..\n");
          }
          else
          {
            if (vrbflg)
              fprintf(stdout, "Estimating PSDs...\n");

            /* compute uncalibrated PSDs */
            LAL_CALL(LALREAL4AverageSpectrum(&status, &psdTempOne, \
                  segmentOne, &specparPSD), &status);
            LAL_CALL(LALREAL4AverageSpectrum(&status, &psdTempTwo, \
                  segmentTwo, &specparPSD), &status);

            if (vrbflg)
            {
              fprintf(stdout, "Getting appropriate frequency band for " \
                  "PSDs..\n");
            }

            /* reduce to the optimal filter frequency range */
            for (i = 0; i < filterLength; i++)
            {
              psdOne.data->data[i] =  psdTempOne.data->data[i + numFMin];
              psdTwo.data->data[i] =  psdTempTwo.data->data[i + numFMin];
            }

            if (vrbflg)
              fprintf(stdout, "Generating inverse noise...\n");

            /* compute inverse calibrate PSDs */
            LAL_CALL(LALStochasticInverseNoiseCal(&status, \
                  &inverseNoiseOutOne, &inverseNoiseInOne), &status);
            LAL_CALL(LALStochasticInverseNoiseCal(&status, \
                  &inverseNoiseOutTwo, &inverseNoiseInTwo), &status);

            /* sum over calibrated PSDs for average */
            for (i = 0; i < filterLength; i++)
            {
              calPsdOne->data[i] = calPsdOne->data[i] + \
                                 1. / calInvPsdOne.data->data[i];
              calPsdTwo->data[i] = calPsdTwo->data[i] + \
                                 1. / calInvPsdTwo.data->data[i];
            }
          }
        }

        /* average calibrated PSDs and take inverse */
        for (i = 0; i < filterLength; i++)
        {
          if (middle_segment_flag == 0)
          {
            calPsdOne->data[i] = calPsdOne->data[i] / (REAL4)(numSegments - 1);
            calPsdTwo->data[i] = calPsdTwo->data[i] / (REAL4)(numSegments - 1);
          }
          else
          {
            calPsdOne->data[i] = calPsdOne->data[i] / (REAL4)numSegments;
            calPsdTwo->data[i] = calPsdTwo->data[i] / (REAL4)numSegments;
          }
          calInvPsdOne.data->data[i] = 1. / calPsdOne->data[i];
          calInvPsdTwo.data->data[i] = 1. / calPsdTwo->data[i];
        }

        if (vrbflg)
          fprintf(stdout, "Normalising optimal filter...\n");

        /* compute variance and normalisation for optimal filter */
        LAL_CALL(LALStochasticOptimalFilterNormalization(&status, \
              &normOutput, &normInput, &normParams), &status);
        lambda = (REAL8)(normLambda.value * \
            pow(10.,normLambda.units.powerOfTen));
        varTheo = (REAL8)(segmentDuration * normSigma.value * \
            pow(10.,normSigma.units.powerOfTen));

        if (vrbflg)
          fprintf(stdout, "Generating optimal filter...\n");

        /* build optimal filter */
        optFilter.epoch = gpsStartTime;
        LAL_CALL(LALStochasticOptimalFilterCal(&status, &optFilter, \
              &optFilterIn, &normLambda), &status);

        /* analyse middle segment */
        gpsStartTime.gpsSeconds = startTime + (interLoop + segMiddle) * \
                                  segmentDuration;
        segmentOne->epoch = gpsStartTime;
        segmentTwo->epoch = gpsStartTime;

        if (vrbflg)
        {
          fprintf(stdout, "analysing segment at GPS %d\n", \
              gpsStartTime.gpsSeconds);
        }

        for (i = 0; i < segmentLength; i++)
        {
          segmentOne->data->data[i] = segOne[segMiddle]->data[i];
          segmentTwo->data->data[i] = segTwo[segMiddle]->data[i];
        }

        /* zero pad and fft */
        LAL_CALL(LALSZeroPadAndFFT(&status, &hBarTildeOne, segmentOne, \
              &zeroPadParams), &status);
        LAL_CALL(LALSZeroPadAndFFT(&status, &hBarTildeTwo, segmentTwo, \
              &zeroPadParams), &status);

        if (vrbflg)
          fprintf(stdout, "Generating cross correlation spectrum...\n");

        /* cc spectrum */
        for (i = 0; i < filterLength; i++)
        {
          responseOne.data->data[i] = respOne[segMiddle]->data[i];
          responseTwo.data->data[i] = respTwo[segMiddle]->data[i];
        }

        LAL_CALL(LALStochasticCrossCorrelationSpectrumCal(&status, \
              &ccSpectrum, &ccIn, epochsMatch), &status);

        /* save */
        if (vrbflg)
          LALCPrintFrequencySeries(&ccSpectrum, "ccSpectrum.dat");

        /* cc statistic */
        LAL_CALL(LALStochasticCrossCorrelationStatisticCal(&status, &ccStat, \
              &ccIn,epochsMatch), &status);
        y = (REAL8)(ccStat.value * pow(10., ccStat.units.powerOfTen));

        /* save */
        if (vrbflg)
        {
          fprintf(stdout, "interval %d:\n", interLoop + 1);
          fprintf(stdout, "GPS time = %d\n", gpsStartTime.gpsSeconds);
          fprintf(stdout, "y = %e\n", y);
          fprintf(stdout, "sigmaTheo = %e\n", sqrt(varTheo));
          fprintf(stdout, "varTheo = %e\n", varTheo);
        }

        /* output to file */
        out = fopen(outputFilename, "a");
        fprintf(out,"%d %e %e\n", gpsStartTime.gpsSeconds, y, \
            sqrt(varTheo));
        fclose(out);
      }
    }

    lal_errhandler = LAL_ERR_EXIT;
  }

  /* save out any flags to the process params table */
  if (middle_segment_flag)
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
                      calloc(1, sizeof(ProcessParamsTable));
    LALSnprintf(this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
        PROGRAM_NAME);
    LALSnprintf(this_proc_param->param, LIGOMETA_PARAM_MAX, \
        "--middle-segment");
    LALSnprintf(this_proc_param->type, LIGOMETA_TYPE_MAX, "string");
    LALSnprintf(this_proc_param->value, LIGOMETA_VALUE_MAX, " ");
  }
  if (inject_flag)
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
                      calloc(1, sizeof(ProcessParamsTable));
    LALSnprintf(this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
        PROGRAM_NAME);
    LALSnprintf(this_proc_param->param, LIGOMETA_PARAM_MAX, \
        "--inject");
    LALSnprintf(this_proc_param->type, LIGOMETA_TYPE_MAX, "string");
    LALSnprintf(this_proc_param->value, LIGOMETA_VALUE_MAX, " ");
  }
  if (apply_mask_flag)
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
                      calloc(1, sizeof(ProcessParamsTable));
    LALSnprintf(this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
        PROGRAM_NAME);
    LALSnprintf(this_proc_param->param, LIGOMETA_PARAM_MAX, \
        "--apply-mask");
    LALSnprintf(this_proc_param->type, LIGOMETA_TYPE_MAX, "string");
    LALSnprintf(this_proc_param->value, LIGOMETA_VALUE_MAX, " ");
  }
  if (high_pass_flag)
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
                      calloc(1, sizeof(ProcessParamsTable));
    LALSnprintf(this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
        PROGRAM_NAME);
    LALSnprintf(this_proc_param->param, LIGOMETA_PARAM_MAX, \
        "--high-pass-filter");
    LALSnprintf(this_proc_param->type, LIGOMETA_TYPE_MAX, "string");
    LALSnprintf(this_proc_param->value, LIGOMETA_VALUE_MAX, " ");
  }
  if (overlap_hann_flag)
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
                      calloc(1, sizeof(ProcessParamsTable));
    LALSnprintf(this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
        PROGRAM_NAME);
    LALSnprintf(this_proc_param->param, LIGOMETA_PARAM_MAX, \
        "--overlap-hann");
    LALSnprintf(this_proc_param->type, LIGOMETA_TYPE_MAX, "string");
    LALSnprintf(this_proc_param->value, LIGOMETA_VALUE_MAX, " ");
  }
  if (recentre_flag)
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
                      calloc(1, sizeof(ProcessParamsTable));
    LALSnprintf(this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
        PROGRAM_NAME);
    LALSnprintf(this_proc_param->param, LIGOMETA_PARAM_MAX, \
        "--recentre");
    LALSnprintf(this_proc_param->type, LIGOMETA_TYPE_MAX, "string");
    LALSnprintf(this_proc_param->value, LIGOMETA_VALUE_MAX, " ");
  }
  if (vrbflg)
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
                      calloc(1, sizeof(ProcessParamsTable));
    LALSnprintf(this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
        PROGRAM_NAME);
    LALSnprintf(this_proc_param->param, LIGOMETA_PARAM_MAX, \
        "--verbose");
    LALSnprintf(this_proc_param->type, LIGOMETA_TYPE_MAX, "string");
    LALSnprintf(this_proc_param->value, LIGOMETA_VALUE_MAX, " ");
  }

  /* add the xml comment, if specified */
  if (!*comment)
  {
    LALSnprintf(proctable.processTable->comment, LIGOMETA_COMMENT_MAX, " ");
  }
  else
  {
    LALSnprintf(proctable.processTable->comment, LIGOMETA_COMMENT_MAX, "%s", \
        comment);
  }
    
  /* delete empty first entry in process params table */
  this_proc_param = procparams.processParamsTable;
  procparams.processParamsTable = procparams.processParamsTable->next;
  free(this_proc_param);

  /* write out xml */
  if (vrbflg)
    fprintf(stdout, "Writing output XML files...\n");

  /* opening xml file stream */
  memset(&xmlStream, 0, sizeof(LIGOLwXMLStream));
  LAL_CALL(LALOpenLIGOLwXMLFile(&status, &xmlStream, xmlFileName), &status);

  /* write out process and process params tables */
  LAL_CALL(LALGPSTimeNow(&status, &(proctable.processTable->end_time), \
        &accuracy), &status);
  LAL_CALL(LALBeginLIGOLwXMLTable(&status, &xmlStream, process_table), \
      &status);
  LAL_CALL(LALWriteLIGOLwXMLTable(&status, &xmlStream, proctable, \
        process_table), &status);
  LAL_CALL(LALEndLIGOLwXMLTable(&status, &xmlStream), &status);
  free(proctable.processTable);

  /* write the process params table */
  LAL_CALL(LALBeginLIGOLwXMLTable(&status, &xmlStream, \
        process_params_table), &status);
  LAL_CALL(LALWriteLIGOLwXMLTable(&status, &xmlStream, procparams, \
        process_params_table), &status);
  LAL_CALL(LALEndLIGOLwXMLTable(&status, &xmlStream), &status);

  /* close xml file */
  LAL_CALL(LALCloseLIGOLwXMLFile(&status, &xmlStream), &status);

  /* cleanup */
  LAL_CALL(LALDestroyREAL4TimeSeries(&status, segmentOne), &status);
  LAL_CALL(LALDestroyREAL4TimeSeries(&status, segmentTwo), &status);
  LAL_CALL(LALDestroyREAL4TimeSeries(&status, segmentPadOne), &status);
  LAL_CALL(LALDestroyREAL4TimeSeries(&status, segmentPadTwo), &status);

  LAL_CALL(LALDestroyRealFFTPlan(&status, &(specparPSD.plan)), &status);
  LAL_CALL(LALDestroyRealFFTPlan(&status, &fftDataPlan), &status);
  LAL_CALL(LALDestroyVector(&status, &(psdTempOne.data)), &status);
  LAL_CALL(LALDestroyVector(&status, &(psdTempTwo.data)), &status);
  LAL_CALL(LALDestroyVector(&status, &(psdOne.data)), &status);
  LAL_CALL(LALDestroyVector(&status, &(psdTwo.data)), &status);
  LAL_CALL(LALDestroyVector(&status, &(calPsdOne)), &status);
  LAL_CALL(LALDestroyVector(&status, &(calPsdTwo)), &status);
  LAL_CALL(LALCDestroyVector(&status, &(responseTempOne.data)), &status);
  LAL_CALL(LALCDestroyVector(&status, &(responseTempTwo.data)), &status);
  LAL_CALL(LALCDestroyVector(&status, &(responseOne.data)), &status);
  LAL_CALL(LALCDestroyVector(&status, &(responseTwo.data)), &status);
  LAL_CALL(LALDestroyVector(&status, &(optFilter.data)), &status);
  LAL_CALL(LALDestroyVector(&status, &(calInvPsdOne.data)), &status);
  LAL_CALL(LALDestroyVector(&status, &(calInvPsdTwo.data)), &status);
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
  LAL_CALL(LALCDestroyVector(&status, &(hBarTildeOne.data)), &status);
  LAL_CALL(LALCDestroyVector(&status, &(hBarTildeTwo.data)), &status);
  if (inject_flag)
  {
    LAL_CALL(LALDestroyREAL4TimeSeries(&status, SimStochBGOne), &status);
    LAL_CALL(LALDestroyREAL4TimeSeries(&status, SimStochBGTwo), &status);
    LAL_CALL(LALCDestroyVector(&status, &(MCresponseOne.data)), &status);
    LAL_CALL(LALCDestroyVector(&status, &(MCresponseTwo.data)), &status);
    LAL_CALL(LALDestroyVector(&status, &(MComegaGW.data)), &status);
  }
  /*
  for (i = 0; i <numSegments; i++)
  {
    LAL_CALL(LALCDestroyVector(&status, &(respOne[i])), &status);
    LAL_CALL(LALCDestroyVector(&status, &(respTwo[i])), &status);
    LAL_CALL(LALDestroyVector(&status, &(segPadOne[i])), &status);
    LAL_CALL(LALDestroyVector(&status, &(segPadTwo[i])), &status);
    LAL_CALL(LALDestroyVector(&status, &(segOne[i])), &status);
    LAL_CALL(LALDestroyVector(&status, &(segTwo[i])), &status);
    if (inject_flag)
    {
      LAL_CALL(LALCDestroyVector(&status, &(MCrespOne[i])), &status);
      LAL_CALL(LALCDestroyVector(&status, &(MCrespTwo[i])), &status);
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
  " --help                        print this message\n"\
  " --version                     display version\n"\
  " --verbose                     verbose mode\n"\
  " --debug-level N               set lalDebugLevel\n"\
  " --gps-start-time N            GPS start time\n"\
  " --gps-end-time N              GPS end time\n"\
  " --interval-duration N         interval duration\n"\
  " --segment-duration N          segment duration\n"\
  " --sample-rate N               sample rate\n"\
  " --resample-rate N             resample rate\n"\
  " --f-min N                     minimal frequency\n"\
  " --f-max N                     maximal frequency\n"\
  " --ifo-one IFO                 ifo for first stream\n"\
  " --ifo-two IFO                 ifo for second stream\n"\
  " --channel-one CHANNEL         channel for first stream\n"\
  " --channel-two CHANNEL         channel for second stream\n"\
  " --frame-cache-one FILE        cache file for first stream\n"\
  " --frame-cache-two FILE        cache file for second stream\n"\
  " --calibration-cache-one FILE  first stream calibration cache\n"\
  " --calibration-cache-two FILE  second stream calibration cache\n"\
  " --calibration-offset N        calibration offset\n"\
  " --apply-mask                  apply frequency masking\n"\
  " --mask-bin N                  number of bins to mask\n"\
  " --overlap-hann                overlaping hann windows\n"\
  " --hann-duration N             hann duration\n"\
  " --high-pass-filter            apply high pass filtering\n"\
  " --hpf-frequency N             high pass filter knee frequency\n"\
  " --hpf-attenuation N           high pass filter attenuation\n"\
  " --hpf-order N                 high pass filter order\n"\
  " --recentre                    recentre jobs\n"\
  " --middle-segment              use middle segment in PSD estimation\n"\
  " --inject                      inject a signal into the data\n"\
  " --scale-factor N              scale factor for injection\n"\
  " --seed N                      random seed\n"\
  " --trials N                    number of trials for Monte Carlo\n"\
  " --output-dir DIR              directory for output files\n"\
  " --geo-hpf-frequency N         GEO high pass filter knee frequency\n"\
  " --geo-hpf-attenuation N       GEO high pass filter attenuation\n"\
  " --geo-hpf-order N             GEO high pass filter order\n"\
  " --alpha N                     exponent on filter spectrum\n"\
  " --f-ref N                     reference frequency for filter spectrum\n"\
  " --omega0 N                    reference omega_0 for filter spectrum\n"\
  " --comment STRING              set the process table to STRING\n"\
  " --user-tag STRING             set the process_params usertag to STRING\n"

/* parse command line options */
static void parseOptions(INT4 argc, CHAR *argv[])
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
      {"verbose", no_argument, &vrbflg, 1},
      {"recentre", no_argument, &recentre_flag,1},
      /* options that don't set a flag */
      {"help", no_argument, 0, 'a'},
      {"gps-start-time", required_argument, 0, 'b'},
      {"gps-end-time", required_argument, 0, 'c'},
      {"interval-duration", required_argument, 0, 'd'},
      {"segment-duration", required_argument, 0, 'e'},
      {"sample-rate", required_argument, 0, 'f'},
      {"resample-rate", required_argument, 0, 'g'},
      {"f-min", required_argument, 0, 'h'},
      {"f-max", required_argument, 0, 'i'},
      {"hann-duration", required_argument, 0, 'j'},
      {"hpf-frequency", required_argument, 0, 'k'},
      {"hpf-attenuation", required_argument, 0, 'l'},
      {"hpf-order", required_argument, 0, 'm'},
      {"geo-hpf-frequency", required_argument, 0, 'n'},
      {"geo-hpf-attenuation", required_argument, 0, 'o'},
      {"geo-hpf-order", required_argument, 0, 'p'},
      {"ifo-one", required_argument, 0, 'q'},
      {"ifo-two", required_argument, 0, 'r'},
      {"channel-one", required_argument, 0, 's'},
      {"channel-two", required_argument, 0, 't'},
      {"frame-cache-one", required_argument, 0, 'u'},
      {"frame-cache-two", required_argument, 0, 'v'},
      {"calibration-cache-one", required_argument, 0, 'w'},
      {"calibration-cache-two", required_argument, 0, 'x'},
      {"calibration-offset", required_argument, 0, 'y'},
      {"mask-bin", required_argument, 0, 'z'},
      {"scale-factor", required_argument, 0, 'A'},
      {"seed", required_argument, 0, 'B'},
      {"trials", required_argument, 0, 'C'},
      {"output-dir", required_argument, 0, 'D'},
      {"debug-level", required_argument, 0, 'E'},
      {"version", no_argument, 0, 'F'},
      {"alpha", required_argument, 0, 'G'},
      {"f-ref", required_argument, 0, 'H'},
      {"omega0", required_argument, 0, 'I'},
      {"comment", required_argument, 0, 'J'},
      {"user-tag", required_argument, 0, 'K'},
      {0, 0, 0, 0}
    };

    /* getopt_long stores the option here */
    int option_index = 0;
    size_t optarg_len;

    c = getopt_long(argc, argv, \
        "ab:c:d:e:f:g:h:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:" \
        "A:B:C:D:E:FG:H:I:J:", long_options, &option_index);

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

      case 'a':
        /* help */
        fprintf(stdout, USAGE);
        exit(0);
        break;

      case 'b':
        /* start time */
        startTime = atoi(optarg);

        /* check */
        if (startTime < 441217609)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "GPS start time is prior to 1 January 1994 00:00:00 UTC " \
              "(%d specified)\n", long_options[option_index].name, \
              startTime);
          exit(1);
        }
        if (startTime > 999999999)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "GPS start time is after 14 September 2011 01:46:26 UTC " \
              "(%d specified)\n", long_options[option_index].name, \
              startTime);
          exit(1);
        }

        ADD_PROCESS_PARAM("int", "%ld", startTime);

        break;

      case 'c':
        /* end time */
        endTime = atoi(optarg);

        /* check */
        if (endTime < 441217609)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "GPS end time is prior to 1 January 1994 00:00:00 UTC " \
              "(%d specified)\n", long_options[option_index].name, \
              endTime);
          exit(1);
        }
        if (endTime > 999999999)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "GPS end time is after 14 September 2011 01:46:26 UTC " \
              "(%d specified)\n", long_options[option_index].name, \
              endTime);
          exit(1);
        }

        ADD_PROCESS_PARAM("int", "%ld", endTime);

        break;

      case 'd':
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

        ADD_PROCESS_PARAM("int", "%d", intervalDuration);

        break;

      case 'e':
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

        ADD_PROCESS_PARAM("int", "%d", segmentDuration);

        break;

      case 'f':
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

        ADD_PROCESS_PARAM("int", "%d", sampleRate);

        break;

      case 'g':
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

        ADD_PROCESS_PARAM("int", "%d", resampleRate);

        break;

      case 'h':
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

        ADD_PROCESS_PARAM("int", "%d", fMin);

        break;

      case 'i':
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

        ADD_PROCESS_PARAM("int", "%d", fMax);

        break;

      case 'j':
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

        ADD_PROCESS_PARAM("int", "%d", hannDuration);

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

        ADD_PROCESS_PARAM("float", "%e", highPassFreq);

        break;

      case 'l':
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

        ADD_PROCESS_PARAM("float", "%e", highPassAtten);

        break;

      case 'm':
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

        ADD_PROCESS_PARAM("int", "%d", highPassOrder);

        break;

      case 'n':
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

        ADD_PROCESS_PARAM("float", "%e", geoHighPassFreq);

        break;

      case 'o':
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

        ADD_PROCESS_PARAM("float", "%e", geoHighPassAtten);

        break;

      case 'p':
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

        ADD_PROCESS_PARAM("int", "%d", geoHighPassOrder);

        break;

      case 'q':
        /* ifo for first stream */
        optarg_len = strlen(optarg) + 1;
        ifoOne = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        strncpy(ifoOne, optarg, optarg_len);

        /* check and set site */
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

        ADD_PROCESS_PARAM("string", "%s", ifoOne);

        break;

      case 'r':
        /* ifo for second stream */
        optarg_len = strlen(optarg) + 1;
        ifoTwo = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        strncpy(ifoTwo, optarg, optarg_len);

        /* check and set site */
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

        ADD_PROCESS_PARAM("string", "%s", ifoTwo);

        break;

      case 's':
        /* channel one */
        optarg_len = strlen(optarg) + 4;
        channelOneTemp = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        channelOne = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        strncpy(channelOneTemp, optarg, optarg_len);
        ADD_PROCESS_PARAM("string", "%s", channelOneTemp);
        break;

      case 't':
        /* channel two */
        optarg_len = strlen(optarg) + 4;
        channelTwoTemp = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        channelTwo = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        strncpy(channelTwoTemp, optarg, optarg_len);
        ADD_PROCESS_PARAM("string", "%s", channelTwoTemp);
        break;

      case 'u':
        /* data cache one */
        optarg_len = strlen(optarg) + 1;
        frameCacheOne = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        strncpy(frameCacheOne, optarg, optarg_len);

        /* check that file exists */
        if ((stat(frameCacheOne, &fileStatus) == -1) && (errno = ENOENT))
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "File does not exist: (%s specified)\n", \
              long_options[option_index].name, frameCacheOne);
          exit(1);
        }

        ADD_PROCESS_PARAM("string", "%s", frameCacheOne);

        break;

      case 'v':
        /* data cache two */
        optarg_len = strlen(optarg) + 1;
        frameCacheTwo = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        strncpy(frameCacheTwo, optarg, optarg_len);

        /* check that file exists */
        if ((stat(frameCacheTwo, &fileStatus) == -1) && (errno = ENOENT))
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "File does not exist: (%s specified)\n", \
              long_options[option_index].name, frameCacheTwo);
          exit(1);
        }

        ADD_PROCESS_PARAM("string", "%s", frameCacheTwo);

        break;

      case 'w':
        /* calibration cache one */
        optarg_len = strlen(optarg) + 1;
        calCacheOne = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        strncpy(calCacheOne, optarg, optarg_len);

        /* check that file exists */
        if ((stat(calCacheOne, &fileStatus) == -1) && (errno = ENOENT))
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "File does not exist: (%s specified)\n", \
              long_options[option_index].name, calCacheOne);
          exit(1);
        }

        ADD_PROCESS_PARAM("string", "%s", calCacheOne);

        break;

      case 'x':
        /* calibration cache two */
        optarg_len = strlen(optarg) + 1;
        calCacheTwo = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        strncpy(calCacheTwo, optarg, optarg_len);

        /* check that file exists */
        if ((stat(calCacheTwo, &fileStatus) == -1) && (errno = ENOENT))
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "File does not exist: (%s specified)\n", \
              long_options[option_index].name, calCacheTwo);
          exit(1);
        }

        ADD_PROCESS_PARAM("string", "%s", calCacheTwo);

        break;

      case 'y':
        /* calibration time offset */
        calibOffset = atoi(optarg);
        ADD_PROCESS_PARAM("int", "%d", calibOffset);
        break;

      case 'z':
        /* number of bins to mask for frequency mask */
        maskBin = atoi(optarg);

        /* check */
        if (maskBin <= 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Number of bins to mask must be greater than 0: " \
              "(%d specified)\n", long_options[option_index].name, maskBin);
          exit(1);
        }

        ADD_PROCESS_PARAM("int", "%d", maskBin);

        break;

      case 'A':
        /* scale factor */
        scaleFactor = atof(optarg);
        ADD_PROCESS_PARAM("float", "%e", scaleFactor);
        break;

      case 'B':
        /* seed */
        seed = atoi(optarg);
        ADD_PROCESS_PARAM("float", "%e", seed);
        break;

      case 'C':
        /* number of trials */
        NLoop = atoi(optarg);

        /* check */
        if (NLoop <= 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Number of trials to must be greater than 0: " \
              "(%d specified)\n", long_options[option_index].name, NLoop);
          exit(1);
        }

        ADD_PROCESS_PARAM("int", "%d", NLoop);

        break;

      case 'D':
        /* directory for output files */
        optarg_len = strlen(optarg) + 1;
        outputFilePath = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        strncpy(outputFilePath, optarg, optarg_len);

        /* check */
        if ((stat(outputFilePath, &fileStatus) == -1) && (errno = ENOENT))
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Directory does not exist: (%s specified)\n", \
              long_options[option_index].name, outputFilePath);
          exit(1);
        }

        ADD_PROCESS_PARAM("string", "%s", outputFilePath);

        break;

      case 'E':
        /* set debug level */
        set_debug_level( optarg );
        ADD_PROCESS_PARAM("string", "%s", optarg);
        break;

      case 'F':
        /* display version info and exit */
        fprintf(stdout, "Standalone SGWB Search Engine\n" CVS_ID "\n");
        exit(0);
        break;

      case 'G':
        /* filter spectrum exponent */
        alpha = atof(optarg);
        ADD_PROCESS_PARAM("float", "%e", alpha);
        break;

      case 'H':
        /* filter reference frequency */
        fRef = atof(optarg);

        /* check */
        if (fRef < 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Reference frequency must be greater than 0: " \
              "(%f specified)\n", long_options[option_index].name, fRef);
          exit(1);
        }

        ADD_PROCESS_PARAM("float", "%e", fRef);

        break;

      case 'I':
        /* filter reference omega */
        omegaRef = atof(optarg);

        /* check */
        if (omegaRef <= 0)
        {
          fprintf(stderr, "Invalid argument to --%s:\n" \
              "Reference omega_0 must be positive: (%f specified)\n", \
              long_options[option_index].name, omegaRef);
          exit(1);
        }

        ADD_PROCESS_PARAM("float", "%e", omegaRef);

        break;

      case 'J':
        /* xml comment */
        if (strlen(optarg) > LIGOMETA_COMMENT_MAX - 1)
        {
          fprintf(stderr, "invalid argument to --%s:\n" \
              "comment must be less than %d characters\n", \
              long_options[option_index].name, LIGOMETA_COMMENT_MAX);
          exit(1);
        }
        else
        {
          LALSnprintf(comment, LIGOMETA_COMMENT_MAX, "%s", optarg);
        }
        break;

      case 'K':
        /* process_params usertag */
        optarg_len = strlen(optarg) + 1;
        userTag = (CHAR*)calloc(optarg_len, sizeof(CHAR));
        strncpy(userTag, optarg, optarg_len);

        /* add to process_params table */
        this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
                          calloc(1, sizeof(ProcessParamsTable));
        LALSnprintf(this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
            PROGRAM_NAME);
        LALSnprintf(this_proc_param->param, LIGOMETA_PARAM_MAX, "--user-tag");
        LALSnprintf(this_proc_param->type, LIGOMETA_TYPE_MAX, "string");
        LALSnprintf(this_proc_param->value, LIGOMETA_VALUE_MAX, "%s", optarg);
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

  /* start/end time */
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

  /* interval duration */
  if (intervalDuration == -1)
  {
    fprintf(stderr, "--interval-duration must be specified\n");
    exit(1);
  }

  /* segment duration */
  if (segmentDuration == -1)
  {
    fprintf(stderr, "--segment-duration must be specified\n");
    exit(1);
  }

  /* sample rates */
  if (sampleRate == -1)
  {
    fprintf(stderr, "--sample-rate must be specified\n");
    exit(1);
  }
  if (resampleRate == -1)
  {
    fprintf(stderr, "--resample-rate must be specified\n");
    exit(1);
  }

  /* min/max frequency */
  if (fMin == -1)
  {
    fprintf(stderr, "--f-min must be specified\n");
    exit(1);
  }
  if (fMax == -1)
  {
    fprintf(stderr, "--f-max must be specified\n");
    exit(1);
  }

  /* ifos */
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

  /* channels */
  if (channelOne == NULL)
  {
    fprintf(stderr, "--channel-one must be specified\n");
    exit(1);
  }
  if (channelTwo == NULL)
  {
    fprintf(stderr, "--channel-two must be specified\n");
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
    if (frameCacheTwo == NULL)
    {
      fprintf(stderr, "--frame-cache-two must be specified\n");
      exit(1);
    }
  }

  /* calibration cache */
  if (strncmp(ifoOne, "G1", 2) != 0)
  {
    if (calCacheOne == NULL)
    {
      fprintf(stderr, "--calibration-cache-one must be specified\n");
      exit(1);
    }
  }
  if (strncmp(ifoTwo, "G1", 2) != 0)
  {
    if (calCacheTwo == NULL)
    {
      fprintf(stderr, "--calibration-cache-two must be specified\n");
      exit(1);
    }
  }

  /* calibration offset */
  if (calibOffset == -1)
  {
    fprintf(stderr, "--calibration-offset must be specified\n");
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

  /* hann duration */
  if (overlap_hann_flag)
  {
    if (hannDuration != -1)
    {
      fprintf(stderr, "Overlapping Hann windows specified, --hann-duration " \
          "will be ignored\n");
    }
  }
  else
  {
    if (hannDuration == -1)
    {
      fprintf(stderr, "--hann-duration must be specified\n");
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

  /* injections */
  if (inject_flag)
  {
    if (scaleFactor == -1)
    {
      fprintf(stderr, "--scale-factor must be specified\n");
      exit(1);
    }
    if (seed == -1)
    {
      fprintf(stderr, "--seed must be specified\n");
      exit(1);
    }
  }

  /* GEO high pass filter */
  if ((strncmp(ifoOne, "G1", 2) == 0) || (strncmp(ifoTwo, "G1", 2) == 0))
  {
    if (geoHighPassFreq == -1)
    {
      fprintf(stderr, "--geo-hpf-frequency must be specified\n");
      exit(1);
    }
    if (geoHighPassAtten == -1)
    {
      fprintf(stderr, "--geo-hpf-attenuation must be specified\n");
      exit(1);
    }
    if (geoHighPassOrder == -1)
    {
      fprintf(stderr, "--geo-hpf-order must be specified\n");
      exit(1);
    }
  }

  /* check for sensible arguments */

  /* start time same as stop time */
  if (startTime == endTime)
  {
    fprintf(stderr, "Start time same as end time; no analysis to perform\n");
    exit(1);
  }

  /* stop time before start time */
  if (startTime > endTime)
  {
    fprintf(stderr, "Invalid start/end time; end time (%d) is before " \
        "start time (%d)\n", endTime, startTime);
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
    fprintf(stderr, "Reference frequency (%f Hz) is greater than maximum " \
        "frequency (%d Hz)\n", fRef, fMax);
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

  /* output dir */
  if (outputFilePath == NULL)
  {
    fprintf(stderr, "Output file path not set, using current directory\n");
    outputFilePath = "./";
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
  UINT4 j;

  /* variables */
  FrCache *frCacheOne = NULL;
  FrStream *frStreamOne = NULL;
  FrCache *frCacheTwo = NULL;
  FrStream *frStreamTwo = NULL;
  FrChanIn frChanInOne, frChanInTwo;
  REAL8TimeSeries *dataStreamGEO;
  REAL4TimeSeries *dataStreamOne;
  REAL4TimeSeries *dataStreamTwo;
  ResampleTSParams resampleParams;
  LIGOTimeGPS bufferStartTime;
  PassBandParamStruc geoHighpassParam;

  /* initialise status pointer */
  INITSTATUS(status, "readDataPair", STOCHASTICC);
  ATTATCHSTATUSPTR(status);

  /* buffer start time */
  bufferStartTime.gpsSeconds = params->start - params->buffer;
  bufferStartTime.gpsNanoSeconds = 0;

  /* set channels */
  frChanInOne.name = channelOne;
  frChanInTwo.name = channelTwo;
  frChanInTwo.type = ADCDataChannel;
  frChanInOne.type = ADCDataChannel;

  if (vrbflg)
    fprintf(stdout, "Allocating memory for raw data streams...\n");

  /* allocate memory */
  LALCreateREAL4TimeSeries(status->statusPtr, &dataStreamOne, \
      "Data Stream One", bufferStartTime, 0, 1./sampleRate, \
      lalDimensionlessUnit, \
      sampleRate * (params->duration + (2 * params->buffer)));
  CHECKSTATUSPTR(status);
  LALCreateREAL4TimeSeries(status->statusPtr, &dataStreamTwo, \
      "Data Stream Two", bufferStartTime, 0, 1./sampleRate, \
      lalDimensionlessUnit, \
      sampleRate * (params->duration + (2 * params->buffer)));
  CHECKSTATUSPTR(status);

  /* allocate memory for geo double precision h(t) */
  LALCreateREAL8TimeSeries(status->statusPtr, &dataStreamGEO, \
      "GEO Data Stream", bufferStartTime, 0, 1./sampleRate, \
      lalDimensionlessUnit, \
      sampleRate * (params->duration + (2 * params->buffer)));
  CHECKSTATUSPTR(status);

  if (vrbflg)
    fprintf(stdout, "Opening first frame cache...\n");

  /* open first frame cache */
  LALFrCacheImport(status->statusPtr, &frCacheOne, frameCacheOne);
  CHECKSTATUSPTR(status);
  LALFrCacheOpen(status->statusPtr, &frStreamOne, frCacheOne);
  CHECKSTATUSPTR(status);

  if (vrbflg)
    fprintf(stdout, "Reading in channel \"%s\"...\n", frChanInOne.name);

  /* read first channel */
  LALFrSeek(status->statusPtr, &(bufferStartTime), frStreamOne);
  CHECKSTATUSPTR(status);

  if (strncmp(ifoOne, "G1", 2) == 0)
  {
    LALFrGetREAL8TimeSeries(status->statusPtr, dataStreamGEO, \
        &frChanInOne, frStreamOne);
    CHECKSTATUSPTR(status);

    /* high pass the GEO data */
    geoHighpassParam.nMax = geoHighPassOrder;
    geoHighpassParam.f1 = -1;
    geoHighpassParam.f2 = (REAL8)geoHighPassFreq;
    geoHighpassParam.a1 = -1;
    geoHighpassParam.a2 = geoHighPassAtten;
    LALButterworthREAL8TimeSeries(status->statusPtr, dataStreamGEO, \
        &geoHighpassParam);
    CHECKSTATUSPTR(status);

    /* cast the GEO data to REAL4 in the channel time series */
    /* which already has the correct amount of memory allocated */
    for (j = 0; j < dataStreamOne->data->length; j++)
    {
      dataStreamOne->data->data[j] = geoScaleFactor * \
                                  (REAL4)(dataStreamGEO->data->data[j]);
    }

    /* re-copy the data paramaters from the GEO channel */
    /* to input data channel */
    LALSnprintf(dataStreamOne->name, LALNameLength * sizeof(CHAR), "%s", \
        dataStreamGEO->name);
    dataStreamOne->epoch = dataStreamGEO->epoch;
    dataStreamOne->deltaT = dataStreamGEO->deltaT;
    dataStreamOne->f0 = dataStreamGEO->f0;
    dataStreamOne->sampleUnits = dataStreamGEO->sampleUnits;
  }
  else
  {
    LALFrGetREAL4TimeSeries(status->statusPtr, dataStreamOne, \
        &frChanInOne, frStreamOne);
    CHECKSTATUSPTR(status);
  }

  /* if site ids are the same then both ifo channels will be in same
   * file */
  if (siteOne == siteTwo)
  {
    if (vrbflg)
    {
      fprintf(stdout, "Reading in channel \"%s\" from same cache...\n", \
        frChanInTwo.name);
    }

    /* read in second channel */
    LALFrSeek(status->statusPtr, &(bufferStartTime), frStreamOne);
    CHECKSTATUSPTR(status);
    LALFrGetREAL4TimeSeries(status->statusPtr, dataStreamTwo, &frChanInTwo, \
        frStreamOne);
    CHECKSTATUSPTR(status);

    if (vrbflg)
      fprintf(stdout, "Closing frame cache...\n");

    /* close frame cache */
    LALFrClose(status->statusPtr, &frStreamOne);
    CHECKSTATUSPTR(status);
  }
  else
  {
    if (vrbflg)
      fprintf(stdout, "Closing first frame cache...\n");

    /* close first frame cache */
    LALFrClose(status->statusPtr, &frStreamOne);
    CHECKSTATUSPTR(status);

    if (vrbflg)
      fprintf(stdout, "Opening second frame cache...\n");

    /* open second frame cache and read in second channel */
    LALFrCacheImport(status->statusPtr, &frCacheTwo, frameCacheTwo);
    CHECKSTATUSPTR(status);
    LALFrCacheOpen(status->statusPtr, &frStreamTwo, frCacheTwo);
    CHECKSTATUSPTR(status);

    if (vrbflg)
      fprintf(stdout, "Reading in channel \"%s\"...\n", frChanInTwo.name);

    /* read in second channel */
    LALFrSeek(status->statusPtr, &(bufferStartTime), frStreamTwo);
    CHECKSTATUSPTR(status);

    if (strncmp(ifoTwo, "G1", 2) == 0)
    {
      LALFrGetREAL8TimeSeries(status->statusPtr, dataStreamGEO, \
          &frChanInTwo, frStreamTwo);
      CHECKSTATUSPTR(status);
      
      /* high pass the GEO data */
      geoHighpassParam.nMax = geoHighPassOrder;
      geoHighpassParam.f1 = -1;
      geoHighpassParam.f2 = (REAL8) geoHighPassFreq;
      geoHighpassParam.a1 = -1;
      geoHighpassParam.a2 = geoHighPassAtten;
      LALButterworthREAL8TimeSeries(status->statusPtr, dataStreamGEO, \
          &geoHighpassParam);
      CHECKSTATUSPTR(status);

      /* cast the GEO data to REAL4 in the channel time series */
      /* which already has the correct amount of memory allocated */
      for (j = 0; j < dataStreamTwo->data->length; j++)
      {
        dataStreamTwo->data->data[j] = geoScaleFactor * \
                                    (REAL4)(dataStreamGEO->data->data[j]);
      }
    }
    else
    {
      LALFrGetREAL4TimeSeries(status->statusPtr, dataStreamTwo, \
          &frChanInTwo, frStreamTwo);
      CHECKSTATUSPTR(status);
    }

    if (vrbflg)
      fprintf(stdout, "Closing second frame cache...\n");

    /* close second frame stream */
    LALFrClose(status->statusPtr, &frStreamTwo);
    CHECKSTATUSPTR(status);
  }

  /* resample */
  if (resampleRate != sampleRate)
  {
    if (vrbflg)
      fprintf(stdout, "Resampling to %d Hz...\n", resampleRate);

    /* set resample parameters */
    resampleParams.deltaT = 1.0 / resampleRate;
    resampleParams.filterType = defaultButterworth;

    /* resample */
    LALResampleREAL4TimeSeries(status->statusPtr, dataStreamOne, \
        &resampleParams);
    CHECKSTATUSPTR(status);
    LALResampleREAL4TimeSeries(status->statusPtr, dataStreamTwo, \
        &resampleParams);
    CHECKSTATUSPTR(status);
  }

  /* build output */
  strncpy(streamPair->streamOne->name, dataStreamOne->name, LALNameLength);
  strncpy(streamPair->streamTwo->name, dataStreamTwo->name, LALNameLength);
  streamPair->streamOne->epoch.gpsSeconds = params->start;
  streamPair->streamTwo->epoch.gpsSeconds = params->start;
  streamPair->streamOne->epoch.gpsNanoSeconds = 0;
  streamPair->streamTwo->epoch.gpsNanoSeconds = 0;
  streamPair->streamOne->deltaT = 1./(REAL8)resampleRate;
  streamPair->streamTwo->deltaT = 1./(REAL8)resampleRate;
  streamPair->streamOne->f0 = streamPair->streamTwo->f0 = 0;
  streamPair->streamOne->sampleUnits = dataStreamOne->sampleUnits;
  streamPair->streamTwo->sampleUnits = dataStreamTwo->sampleUnits;

  /* remove buffer, and hence corruption due to resampling */
  for (i = 0; i < params->duration * resampleRate; i++)
  {
    streamPair->streamOne->data->data[i] = \
      dataStreamOne->data->data[i + (resampleRate * params->buffer)];
    streamPair->streamTwo->data->data[i] = \
      dataStreamTwo->data->data[i + (resampleRate * params->buffer)];
  }

  /* clean up */
  LALDestroyREAL4TimeSeries(status->statusPtr, dataStreamOne);
  CHECKSTATUSPTR(status);
  LALDestroyREAL4TimeSeries(status->statusPtr, dataStreamTwo);
  CHECKSTATUSPTR(status);
  LALDestroyREAL8TimeSeries(status->statusPtr, dataStreamGEO);
  CHECKSTATUSPTR(status);

  /* return status */
  DETATCHSTATUSPTR(status);
  RETURN(status);
}

/*
 * vim: et
 */
