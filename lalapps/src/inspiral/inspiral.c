/*----------------------------------------------------------------------- 
 * 
 * File Name: inspiral.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <config.h>
#include "inspiral.h"
#include "inspiralfrutils.h"
#include "ligolwbank.h"

RCSID( "$Id$" );

#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"

#define CANDLE_MASS1 1.4
#define CANDLE_MASS2 1.4
#define CANDLE_RHOSQ 8.0

/*
 *
 * variables that control program behaviour
 *
 */


/* debugging */
extern int vrbflg;                      /* verbocity of lal function    */

/* input data parameters */
UINT8  gpsStartTimeNS   = 0;            /* input data GPS start time ns */
LIGOTimeGPS gpsStartTime;               /* input data GPS start time    */
UINT8  gpsEndTimeNS     = 0;            /* input data GPS end time ns   */
LIGOTimeGPS gpsEndTime;                 /* input data GPS end time      */
INT4  padData = 0;                      /* saftety margin on input data */
CHAR  *fqChanName       = NULL;         /* name of data channel         */
CHAR  *frInCacheName    = NULL;         /* cache file containing frames */
INT4  numPoints         = -1;           /* points in a segment          */
INT4  numSegments       = -1;           /* number of segments           */
INT4  ovrlap            = -1;           /* overlap between segments     */
CHAR  site[2];                          /* single character site code   */
CHAR  ifo[3];                           /* two character ifo code       */
CHAR *channelName = NULL;               /* channel string               */
UINT4 inputDataLength = 0;              /* number of points in input    */

/* data conditioning parameters */
INT4   resampFiltType   = -1;           /* low pass filter used for res */
INT4   sampleRate       = -1;           /* sample rate of filter data   */
INT4   highPass         = -1;           /* enable high pass on raw data */
REAL4  highPassFreq     = 0;            /* high pass frequency          */
REAL4  fLow             = -1;           /* low frequency cutoff         */
INT4   specType         = -1;           /* use median or mean psd       */
INT4   badMeanPsd       = 0;            /* use a mean with no overlap   */
INT4   invSpecTrunc     = -1;           /* length of inverse spec (s)   */
REAL4  dynRangeExponent = -1;           /* exponent of dynamic range    */
CHAR  *calCacheName     = NULL;         /* location of calibration data */
CHAR  *injectionFile    = NULL;         /* name of file containing injs */

/* matched filter parameters */
CHAR *bankFileName      = NULL;         /* name of input template bank  */
INT4  startTemplate     = -1;           /* index of first template      */
INT4  stopTemplate      = -1;           /* index of last template       */
INT4  numChisqBins      = -1;           /* number of chisq bins         */
REAL4 snrThresh         = -1;           /* signal to noise thresholds   */
REAL4 chisqThresh       = -1;           /* chisq veto thresholds        */
INT4  eventCluster      = -1;           /* perform chirplen clustering  */

/* output parameters */
int    enableOutput     = -1;           /* write out inspiral events    */
int    writeRawData     = 0;            /* write the raw data to a file */
int    writeFilterData  = 0;            /* write post injection data    */
int    writeResponse    = 0;            /* write response function used */
int    writeSpectrum    = 0;            /* write computed psd to file   */
int    writeRhosq       = 0;            /* write rhosq time series      */
int    writeChisq       = 0;            /* write chisq time series      */

/* other command line args */
CHAR comment[LIGOMETA_COMMENT_MAX];     /* process param comment        */

int main( int argc, char *argv[] )
{
  /* lal function variables */
  LALStatus             status = blank_status;
  LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;

  /* frame input data */
  FrCache      *frInCache = NULL;
  FrStream     *frStream = NULL;
  FrChanIn      frChan;

  /* frame output data */
  struct FrFile *frOutFile  = NULL;
  struct FrameH *outFrame   = NULL;
  UINT4          nRhosqFr = 0;
  UINT4          nChisqFr = 0;

  /* raw input data storage */
  REAL4TimeSeries               chan;
  REAL4FrequencySeries          spec;
  COMPLEX8FrequencySeries       resp;
  DataSegmentVector            *dataSegVec = NULL;

  /* structures for preconditioning */
  COMPLEX8FrequencySeries       injResp;
  COMPLEX8FrequencySeries      *injRespPtr;
  ResampleTSParams              resampleParams;
  LALWindowParams               wpars;
  AverageSpectrumParams         avgSpecParams;

  /* findchirp data structures */
  FindChirpInitParams          *fcInitParams   = NULL;
  FindChirpSegmentVector       *fcSegVec       = NULL;
  FindChirpSPDataParams        *fcDataParams   = NULL;
  FindChirpSPTmpltParams       *fcTmpltParams  = NULL;
  FindChirpFilterParams        *fcFilterParams = NULL;
  FindChirpFilterInput         *fcFilterInput  = NULL;
  FindChirpStandardCandle       candle;

  /* inspiral template structures */
  INT4                          numTmplts    = 0;
  InspiralTemplate             *bankHead     = NULL;
  InspiralTemplate             *bankCurrent  = NULL;
  InspiralTemplateNode         *tmpltHead    = NULL;
  InspiralTemplateNode         *tmpltCurrent = NULL;
  InspiralTemplateNode         *tmpltInsert  = NULL;

  /* inspiral events */
  INT4                          numEvents   = 0;
  SnglInspiralTable            *event       = NULL;
  SnglInspiralTable            *eventList   = NULL;
  MetadataTable                 savedEvents;

  /* output data */
  MetadataTable         proctable;
  MetadataTable         procparams;
  MetadataTable         searchsumm;
  MetadataTable         searchsummvars;
  SearchSummvarsTable  *this_search_summvar;
  MetadataTable         summvalue;
  SummValueTable        candleTable;
  ProcessParamsTable   *this_proc_param;
  LIGOLwXMLStream       results;

  /* counters and other variables */
  const LALUnit strainPerCount = {0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};
  UINT4 i;
  INT4  inserted;
  INT4  currentLevel;
  CHAR  fname[256];
  REAL8 inputLengthNS;
  UINT4 numInputPoints;
  const REAL8 epsilon = 1.0e-8;
  UINT4 resampleChan = 0;
  REAL8 tsLength;


  /*
   *
   * initialization
   *
   */


  /* set up inital debugging values */
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "1" );

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) 
    LALCalloc( 1, sizeof(ProcessTable) );
  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->start_time),
        &accuracy ), &status );
  LAL_CALL( populate_process_table( &status, proctable.processTable, 
        PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE ), &status );
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *) 
    LALCalloc( 1, sizeof(ProcessParamsTable) );
  memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );

  /* create the search summary table */
  searchsumm.searchSummaryTable = (SearchSummaryTable *)
    LALCalloc( 1, sizeof(SearchSummaryTable) );

  /* call the argument parse and check function */
  arg_parse_check( argc, argv, procparams );
  for ( this_proc_param = procparams.processParamsTable; this_proc_param->next;
     this_proc_param = this_proc_param->next );

  /* fill the comment, if a user has specified on, or leave it blank */
  if ( ! *comment )
  {
    snprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX, " " );
    snprintf( searchsumm.searchSummaryTable->comment, LIGOMETA_COMMENT_MAX, 
        " " );
  } 
  else 
  {
    snprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX,
        "%s", comment );
    snprintf( searchsumm.searchSummaryTable->comment, LIGOMETA_COMMENT_MAX,
        "%s", comment );
  }

  /* the number of nodes for a standalone job is always 1 */
  searchsumm.searchSummaryTable->nnodes = 1;

  /* make sure the pointer to the first event is null */
  savedEvents.snglInspiralTable = NULL;

  /* create the standard candle and database table */
  summvalue.summValueTable = &candleTable;
  memset( &candleTable, 0, sizeof(SummValueTable) );
  memset( &candle, 0, sizeof(FindChirpStandardCandle) );
  strncpy( candle.ifo, ifo, 2 * sizeof(CHAR) );
  candle.tmplt.mass1 = CANDLE_MASS1;
  candle.tmplt.mass2 = CANDLE_MASS2;
  candle.rhosq       = CANDLE_RHOSQ;
  candle.tmplt.totalMass = candle.tmplt.mass1 + candle.tmplt.mass2;
  candle.tmplt.mu = candle.tmplt.mass1 * candle.tmplt.mass2 / 
    candle.tmplt.totalMass;
  candle.tmplt.eta = candle.tmplt.mu / candle.tmplt.totalMass;


  /*
   *
   * create a (possibly heirarcical) template bank
   *
   */


  /* read in the template bank from a ligo lw xml file */
  numTmplts = InspiralTmpltBankFromLIGOLw( &bankHead, bankFileName,
      startTemplate, stopTemplate );
  if ( numTmplts < 0 )
  {
    fprintf( stderr, "error: unable to read templates from %s\n", 
        bankFileName );
    exit( 1 );
  }
  else if ( numTmplts == 0 )
  {
    fprintf( stdout, "no templates found in template bank file: %s\n"
        "exiting without searching for events.\n", bankFileName );
    goto cleanexit;
  }

  if ( vrbflg ) fprintf( stdout, "parsed %d templates from %s\n", 
      numTmplts, bankFileName );

  /* save the minimal match of the bank in the process params */
  this_proc_param = this_proc_param->next = (ProcessParamsTable *) 
    LALCalloc( 1, sizeof(ProcessParamsTable) ); 
  snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", 
      PROGRAM_NAME ); \
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--minimal-match" );
  snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "float" ); 
  snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%e", 
      bankHead->minMatch );

  /* create the linked list of template nodes for coarse templates */
  for ( bankCurrent = bankHead; bankCurrent; bankCurrent = bankCurrent->next )
  {
    LAL_CALL( LALFindChirpCreateTmpltNode( &status, 
          bankCurrent, &tmpltCurrent ), &status );
    if ( !tmpltHead ) tmpltHead = tmpltCurrent;
  }


  /*
   *
   * read in the input data channel
   *
   */


  /* set the time series parameters of the input data and resample params */
  memset( &resampleParams, 0, sizeof(ResampleTSParams) );
  resampleParams.deltaT = 1.0 / (REAL8) sampleRate;

  /* set the params of the input data time series */
  memset( &chan, 0, sizeof(REAL4TimeSeries) );
  chan.epoch = gpsStartTime;
  chan.epoch.gpsSeconds -= padData; /* subtract pad seconds from start */

  /* open a frame cache or files, seek to required epoch and set chan name */
  if ( frInCacheName )
  {
    LAL_CALL( LALFrCacheImport( &status, &frInCache, frInCacheName), &status );
    LAL_CALL( LALFrCacheOpen( &status, &frStream, frInCache ), &status );
  }
  else
  {
    LAL_CALL( LALFrOpen( &status, &frStream, NULL, "*.gwf" ), &status );
  }
  LAL_CALL( LALFrSeek( &status, &(chan.epoch), frStream ), &status );
  frChan.name = fqChanName;

  /* determine the sample rate of the raw data and allocate enough memory */
  LAL_CALL( LALFrGetREAL4TimeSeries( &status, &chan, &frChan, frStream ),
      &status );

  /* store the input sample rate */
  this_search_summvar = searchsummvars.searchSummvarsTable = 
    (SearchSummvarsTable *) LALCalloc( 1, sizeof(SearchSummvarsTable) );
  LALSnprintf( this_search_summvar->name, LIGOMETA_NAME_MAX * sizeof(CHAR),
      "raw data sample rate" );
  this_search_summvar->value = chan.deltaT;

  /* determine if we need to resample the channel */
  if ( vrbflg )
  {
    fprintf( stdout, "resampleParams.deltaT = %e\n", resampleParams.deltaT );
    fprintf( stdout, "chan.deltaT = %e\n", chan.deltaT );
  }
  if ( ! ( fabs( resampleParams.deltaT - chan.deltaT ) < epsilon ) )
  {
    resampleChan = 1;
    if ( vrbflg )
      fprintf( stdout, "input channel will be resampled\n" );

    if ( resampFiltType == 0 )
    {
      resampleParams.filterType = LDASfirLP;
    }
    else if ( resampFiltType == 1 )
    {
      resampleParams.filterType = defaultButterworth;
    }
  }

  /* determine the number of points to get and create storage forr the data */
  inputLengthNS = 
    (REAL8) ( gpsEndTimeNS - gpsStartTimeNS + 2000000000LL * padData );
  chan.deltaT *= 1.0e9;
  numInputPoints = (UINT4) floor( inputLengthNS / chan.deltaT + 0.5 );
  LAL_CALL( LALSCreateVector( &status, &(chan.data), numInputPoints ), 
      &status );

  if ( vrbflg ) fprintf( stdout, "input channel %s has sample interval "
      "(deltaT) = %e\nreading %d points from frame stream\n", fqChanName, 
      chan.deltaT / 1.0e9, numInputPoints );

  /* read the data channel time series from frames */
  LAL_CALL( LALFrGetREAL4TimeSeries( &status, &chan, &frChan, frStream ),
      &status );
  memcpy( &(chan.sampleUnits), &lalADCCountUnit, sizeof(LALUnit) );

  /* store the start and end time of the raw channel in the search summary */
  searchsumm.searchSummaryTable->in_start_time = chan.epoch;
  LAL_CALL( LALGPStoFloat( &status, &tsLength, &(chan.epoch) ), 
      &status );
  tsLength += chan.deltaT * (REAL8) chan.data->length;
  LAL_CALL( LALFloatToGPS( &status, 
        &(searchsumm.searchSummaryTable->in_end_time), &tsLength ), &status );

  /* close the frame file stream and destroy the cache */
  LAL_CALL( LALFrClose( &status, &frStream ), &status );
  if ( frInCacheName ) LAL_CALL( LALDestroyFrCache( &status, &frInCache ), 
      &status );

  /* write the raw channel data as read in from the frame files */
  if ( writeRawData ) outFrame = fr_add_proc_REAL4TimeSeries( outFrame, 
      &chan, "ct", "RAW" );

  if ( vrbflg ) fprintf( stdout, "read channel %s from frame stream\n"
      "got %d points with deltaT %e\nstarting at GPS time %d sec %d ns\n", 
      chan.name, chan.data->length, chan.deltaT, 
      chan.epoch.gpsSeconds, chan.epoch.gpsNanoSeconds );


  /*
   *
   * generate the response function for the requested time
   *
   */


  /* create storage for the response function */
  memset( &resp, 0, sizeof(COMPLEX8FrequencySeries) );
  LAL_CALL( LALCCreateVector( &status, &(resp.data), numPoints / 2 + 1 ), 
      &status );

  /* set the parameters of the response to match the data */
  resp.epoch = gpsStartTime;
  resp.deltaF = (REAL8) sampleRate / (REAL8) numPoints;
  resp.sampleUnits = strainPerCount;
  strcpy( resp.name, chan.name );

  /* generate the response function for the current time */
  if ( vrbflg ) fprintf( stdout, "generating response at time %d sec %d ns\n",
      resp.epoch.gpsSeconds, resp.epoch.gpsNanoSeconds );
  LAL_CALL( LALExtractFrameResponse( &status, &resp, calCacheName, ifo ),
      &status );

  if ( writeResponse ) outFrame = fr_add_proc_COMPLEX8FrequencySeries( 
      outFrame, &resp, "strain/ct", "RESPONSE" );


  /*
   *
   * inject signals into the raw, unresampled data
   *
   */


  if ( injectionFile )
  {
    INT4 injSafety = 0; /* get injections within 10 mins of the segment */
    int  numInjections = 0;
    SimInspiralTable    *injections = NULL;
    SimInspiralTable    *thisInj = NULL;
    CHAR tmpChName[LALNameLength];

    /* read in the injection data from XML */
    numInjections = SimInspiralTableFromLIGOLw( &injections, injectionFile,
        gpsStartTime.gpsSeconds - injSafety, 
        gpsEndTime.gpsSeconds + injSafety );
    if ( numInjections < 0 )
    {
      fprintf( stderr, "error: cannot read injection file" );
      exit( 1 );
    }

    /* determine if we need a higher resolution response to do the injections */
    if ( resampleChan )
    {
      /* we need a different resolution of response function for injections */
      REAL8 ratioDeltaT = 
        ( chan.deltaT * 1.0e9 ) / ( resampleParams.deltaT * 1.0e9 );
      UINT4 rateRatio = (UINT4) rint ( ratioDeltaT );
      UINT4 rawNumPoints = rateRatio * numPoints;

      memset( &injResp, 0, sizeof(COMPLEX8FrequencySeries) );
      LAL_CALL( LALCCreateVector( &status, &(injResp.data), 
            rawNumPoints / 2 + 1 ), &status );
      injResp.epoch = gpsStartTime;
      injResp.deltaF = 1.0 / ( rawNumPoints * chan.deltaT );
      injResp.sampleUnits = strainPerCount;
      strcpy( injResp.name, chan.name );

      /* generate the response function for the current time */
      if ( vrbflg ) fprintf( stdout, 
          "generating high resolution response at time %d sec %d ns\n"
          "length = %d points, deltaF = %e Hz\n",
          resp.epoch.gpsSeconds, resp.epoch.gpsNanoSeconds,
          injResp.data->length, injResp.deltaF );
      LAL_CALL( LALExtractFrameResponse( &status, &injResp, calCacheName, ifo ),
          &status );

      injRespPtr = &injResp;

      if ( writeResponse ) outFrame = fr_add_proc_COMPLEX8FrequencySeries( 
          outFrame, &resp, "strain/ct", "RESPONSE_INJ" );
    }
    else
    {
      /* the data is already at the correct sample rate, just do injections */
      injRespPtr = &resp;
      memset( &injResp, 0, sizeof(COMPLEX8FrequencySeries) );
    }

    /* inject the signals, preserving the channel name (Tev mangles it) */
    strcpy( tmpChName, chan.name );
    LAL_CALL( LALFindChirpInjectSignals( &status, &chan, injections, 
          injRespPtr ), &status );
    strcpy( chan.name, tmpChName );

    if ( vrbflg ) fprintf( stdout, "injected %d signals from %s into %s\n", 
        numInjections, injectionFile, chan.name );

    while ( injections )
    {
      thisInj = injections;
      injections = injections->next;
      LALFree( thisInj );
    }

    /* write the raw channel data plus injections to the output frame file */
    if ( writeRawData ) outFrame = fr_add_proc_REAL4TimeSeries( outFrame, 
        &chan, "ct", "RAW_INJ" );

    if ( injResp.data )
      LAL_CALL( LALCDestroyVector( &status, &(injResp.data) ), &status );
  }

  
  /*
   *
   * resample the data to the requested rate
   *
   */


  if ( resampleChan )
  {
    if (vrbflg) fprintf( stdout, "resampling input data from %e to %e\n",
        chan.deltaT, resampleParams.deltaT );

    LAL_CALL( LALResampleREAL4TimeSeries( &status, &chan, &resampleParams ),
        &status );

    if ( vrbflg ) fprintf( stdout, "channel %s resampled:\n"
        "%d points with deltaT %e\nstarting at GPS time %d sec %d ns\n", 
        chan.name, chan.data->length, chan.deltaT, 
        chan.epoch.gpsSeconds, chan.epoch.gpsNanoSeconds );

    /* write the resampled channel data as read in from the frame files */
    if ( writeRawData ) outFrame = fr_add_proc_REAL4TimeSeries( outFrame, 
        &chan, "ct", "RAW_RESAMP" );
  }

  /* store the filter data sample rate */
  this_search_summvar = this_search_summvar->next = 
    (SearchSummvarsTable *) LALCalloc( 1, sizeof(SearchSummvarsTable) );
  LALSnprintf( this_search_summvar->name, LIGOMETA_NAME_MAX * sizeof(CHAR),
      "filter data sample rate" );
  this_search_summvar->value = chan.deltaT;


  /* 
   *
   * high pass the data, removed pad from time series and check length of data
   *
   */


  /* iir filter to remove low frequencies from data channel */
  if ( highPass )
  {
    PassBandParamStruc highpassParam;
    highpassParam.nMax = 4;
    highpassParam.f1 = -1.0;
    highpassParam.f2 = highPassFreq;
    highpassParam.a1 = -1.0;
    highpassParam.a2 = 0.1;

    if ( vrbflg ) fprintf( stdout, "highpassing data above %e\n", 
        highPassFreq );

    LAL_CALL( LALButterworthREAL4TimeSeries( &status, &chan, &highpassParam ),
        &status );
  }

  /* remove pad from requested data from start and end of time series */
  memmove( chan.data->data, chan.data->data + padData * sampleRate, 
      (chan.data->length - 2 * padData * sampleRate) * sizeof(REAL4) );
  LALRealloc( chan.data->data, 
      (chan.data->length - 2 * padData * sampleRate) * sizeof(REAL4) );
  chan.data->length -= 2 * padData * sampleRate;
  chan.epoch.gpsSeconds += padData;

  if ( vrbflg ) fprintf( stdout, "after removal of %d second padding at "
      "start and end:\ndata channel sample interval (deltaT) = %e\n"
      "data channel length = %d\nstarting at %d sec %d ns\n", 
      padData , chan.deltaT , chan.data->length, 
      chan.epoch.gpsSeconds, chan.epoch.gpsNanoSeconds );

  if ( writeFilterData ) outFrame = fr_add_proc_REAL4TimeSeries( outFrame, 
      &chan, "ct", "FILTER" );

  /* check data length */
  if ( chan.data->length != inputDataLength )
  {
    fprintf( stderr, "error: computed channel length and requested\n"
        "input data length do not match:\nchan.data->length = %d\n"
        "inputDataLength = %d\nyou have found a bug in the code.\n"
        "please report this to <duncan@gravity.phys.uwm.edu>\n",
        chan.data->length, inputDataLength );
    exit( 1 );
  }

  /* store the start and end time of the filter channel in the search summ */
  /* noting that we don't look for events in the first and last quarter    */
  /* of each findchirp segment of the input data                           */
  LAL_CALL( LALGPStoFloat( &status, &tsLength, &(chan.epoch) ), 
      &status );
  tsLength += (REAL8) (numPoints / 4) * chan.deltaT;
  LAL_CALL( LALFloatToGPS( &status, 
        &(searchsumm.searchSummaryTable->out_start_time), &tsLength ), 
      &status );
  LAL_CALL( LALGPStoFloat( &status, &tsLength, &(chan.epoch) ), 
      &status );
  tsLength += chan.deltaT * ((REAL8) chan.data->length - (REAL8) (numPoints/4));
  LAL_CALL( LALFloatToGPS( &status, 
        &(searchsumm.searchSummaryTable->out_end_time), &tsLength ), 
      &status );


  /* 
   *
   * create and populate findchip initialization structure 
   *
   */


  if ( ! ( fcInitParams = (FindChirpInitParams *) 
        LALCalloc( 1, sizeof(FindChirpInitParams) ) ) )
  {
    fprintf( stderr, "could not allocate memory for findchirp init params\n" );
    exit( 1 );
  }
  fcInitParams->numPoints      = numPoints;
  fcInitParams->numSegments    = numSegments;
  fcInitParams->numChisqBins   = numChisqBins;
  fcInitParams->createRhosqVec = writeRhosq;
  fcInitParams->ovrlap         = ovrlap;

  /* create the data segment vector */
  memset( &spec, 0, sizeof(REAL4FrequencySeries) );
  LAL_CALL( LALSCreateVector( &status, &(spec.data), numPoints / 2 + 1 ), 
      &status );
  LAL_CALL( LALInitializeDataSegmentVector( &status, &dataSegVec,
        &chan, &spec, &resp, fcInitParams ), &status );


  /*
   *
   * power spectrum estimation and data conditioning
   *
   */


  /* create the findchirp data storage */
  LAL_CALL( LALCreateFindChirpSegmentVector( &status, &fcSegVec, 
        fcInitParams ), &status );

  /* initialize data conditioning routines */
  LAL_CALL( LALFindChirpSPDataInit( &status, &fcDataParams, fcInitParams ), 
      &status );
  fcDataParams->invSpecTrunc = invSpecTrunc * sampleRate;
  fcDataParams->fLow = fLow;

  /* compute the windowed power spectrum for the data channel */
  avgSpecParams.window = NULL;
  avgSpecParams.plan   = fcDataParams->fwdPlan;
  switch ( specType )
  {
    case 0:
      avgSpecParams.method = useMean;
      if ( vrbflg ) fprintf( stdout, "computing mean psd" );
      break;
    case 1:
      avgSpecParams.method = useMedian;
      if ( vrbflg ) fprintf( stdout, "computing median psd" );
      break;
  }

  wpars.type = Hann;
  wpars.length = numPoints;
  if ( badMeanPsd )
  {
    avgSpecParams.overlap = 0;
    if ( vrbflg ) fprintf( stdout, " without overlap\n" );
  }
  else
  {
    avgSpecParams.overlap = numPoints / 2;
    if ( vrbflg ) 
      fprintf( stdout, " with overlap %d\n", avgSpecParams.overlap );
  }

  LAL_CALL( LALCreateREAL4Window( &status, &(avgSpecParams.window),
        &wpars ), &status );
  LAL_CALL( LALREAL4AverageSpectrum( &status, &spec, &chan, &avgSpecParams ),
      &status );
  LAL_CALL( LALDestroyREAL4Window( &status, &(avgSpecParams.window) ), 
      &status );
  strcpy( spec.name, chan.name );

  /* write the spectrum data to a file */
  if ( writeSpectrum ) outFrame = fr_add_proc_REAL4FrequencySeries( outFrame, 
      &spec, "ct/sqrtHz", "PSD" );


  /*
   *
   * create the data structures needed for findchirp
   *
   */


  if ( vrbflg ) fprintf( stdout, "initializing findchirp\n" );

  /* initialize the template functions */
  LAL_CALL( LALFindChirpSPTemplateInit( &status, &fcTmpltParams, 
        fcInitParams ), &status );

  fcDataParams->dynRange = fcTmpltParams->dynRange = 
    pow( 2.0, dynRangeExponent );
  fcDataParams->deltaT = fcTmpltParams->deltaT = 1.0 / (REAL4) sampleRate;
  fcTmpltParams->fLow = fLow;

  /* initialize findchirp filter functions */
  LAL_CALL( LALFindChirpFilterInit( &status, &fcFilterParams, fcInitParams ), 
      &status );
  fcFilterParams->deltaT = 1.0 / (REAL4) sampleRate;
  fcFilterParams->computeNegFreq = 0;

  LAL_CALL( LALCreateFindChirpInput( &status, &fcFilterInput, fcInitParams ), 
      &status );
  LAL_CALL( LALFindChirpChisqVetoInit( &status, fcFilterParams->chisqParams, 
        fcInitParams->numChisqBins, fcInitParams->numPoints ), 
      &status );

  /* parse the thresholds */
  fcFilterParams->rhosqThresh = snrThresh * snrThresh;
  fcFilterParams->chisqThresh = chisqThresh;
  fcFilterParams->maximiseOverChirp = eventCluster;


  /*
   *
   * condition data segments for filtering
   *
   */


  if ( vrbflg ) fprintf( stdout, "findchirp conditioning data\n" );
  LAL_CALL( LALFindChirpSPData (&status, fcSegVec, dataSegVec, fcDataParams),
      &status );

  /* compute the standard candle */
  {
      REAL4 cannonDist = 1.0; /* Mpc */
      REAL4 m  = (REAL4) candle.tmplt.totalMass;
      REAL4 mu = (REAL4) candle.tmplt.mu;
      REAL4 distNorm = 2.0 * LAL_MRSUN_SI / (cannonDist * 1e6 * LAL_PC_SI);
      REAL4 candleTmpltNorm = sqrt( (5.0*mu) / 96.0 ) *
        pow( m / (LAL_PI*LAL_PI) , 1.0/3.0 ) *
        pow( LAL_MTSUN_SI / (REAL4) chan.deltaT, -1.0/6.0 );

      distNorm *= fcTmpltParams->dynRange;
      candleTmpltNorm *= candleTmpltNorm;
      candleTmpltNorm *= distNorm * distNorm;

      candle.sigmasq = 4.0 * ( (REAL4) chan.deltaT / (REAL4) numPoints );
      candle.sigmasq *= candleTmpltNorm * fcSegVec->data->segNorm;

      candle.effDistance = sqrt( candle.sigmasq / candle.rhosq );

      if ( vrbflg ) 
        fprintf( stdout, "candle.effDistance = %e Mpc\n", candle.effDistance );
    }
  

  /*
   *
   * hierarchial search engine
   *
   */


  for ( tmpltCurrent = tmpltHead, inserted = 0; tmpltCurrent; 
      tmpltCurrent = tmpltCurrent->next, inserted = 0 )
  {
    /*  generate template */
    LAL_CALL( LALFindChirpSPTemplate( &status, fcFilterInput->fcTmplt, 
          tmpltCurrent->tmpltPtr, fcTmpltParams ), &status );
    fcFilterInput->tmplt = tmpltCurrent->tmpltPtr;

    /* loop over data segments */
    for ( i = 0; i < fcSegVec->length ; ++i )
    {
      /* filter data segment */ 
      if ( fcSegVec->data[i].level == tmpltCurrent->tmpltPtr->level )
      {
        if ( vrbflg ) fprintf( stdout, 
            "filtering segment %d/%d againt template %d/%d (%e,%e)\n", 
            fcSegVec->data[i].number,  fcSegVec->length,
            tmpltCurrent->tmpltPtr->number, numTmplts,
            fcFilterInput->tmplt->mass1, fcFilterInput->tmplt->mass2 );

        fcFilterInput->segment = fcSegVec->data + i;
        LAL_CALL( LALFindChirpFilterSegment( &status, 
              &eventList, fcFilterInput, fcFilterParams ), &status );

        if ( writeRhosq )
        {
          CHAR snrsqStr[LALNameLength];
          LALSnprintf( snrsqStr, LALNameLength*sizeof(CHAR), 
              "SNRSQ_%d", nRhosqFr++ );
          strcpy( fcFilterParams->rhosqVec->name, chan.name );
          outFrame = fr_add_proc_REAL4TimeSeries( outFrame, 
              fcFilterParams->rhosqVec, "none", snrsqStr );
        }

        if ( writeChisq )
        {
          CHAR chisqStr[LALNameLength];
          REAL4TimeSeries chisqts;
          LALSnprintf( chisqStr, LALNameLength*sizeof(CHAR), 
              "CHISQ_%d", nChisqFr++ );
          chisqts.epoch = fcFilterInput->segment->data->epoch;
          memcpy( &(chisqts.name), fcFilterInput->segment->data->name,
              LALNameLength * sizeof(CHAR) );
          chisqts.deltaT = fcFilterInput->segment->deltaT;
          chisqts.data = fcFilterParams->chisqVec;
          outFrame = fr_add_proc_REAL4TimeSeries( outFrame, 
              &chisqts, "none", chisqStr );
        }

      }

      /*  test if filter returned any events */
      if ( eventList )
      {
        if ( vrbflg ) fprintf( stdout, "segment %d rang template %e,%e\n",
            fcSegVec->data[i].number,
            fcFilterInput->tmplt->mass1, fcFilterInput->tmplt->mass2 );

        if ( tmpltCurrent->tmpltPtr->fine != NULL && inserted == 0 )
        {
          if ( vrbflg ) fprintf( stdout, 
              "inserting fine templates into list\n" );

          tmpltInsert = tmpltCurrent;
          inserted = 1;
          fcSegVec->data[i].level += 1;

          for ( bankCurrent = tmpltCurrent->tmpltPtr->fine ; 
              bankCurrent; bankCurrent = bankCurrent->next )
          {
            LAL_CALL( LALFindChirpCreateTmpltNode( &status, 
                  bankCurrent, &tmpltInsert ), &status );
          }

        }
        else if ( ! tmpltCurrent->tmpltPtr->fine )
        {
          if ( vrbflg ) fprintf( stdout, "***>  dumping events  <***\n" );

          if ( ! savedEvents.snglInspiralTable )
          {
            savedEvents.snglInspiralTable = eventList;
          }
          else
          {
            event->next = eventList;
          }
        }
        else
        {
          if ( vrbflg ) fprintf( stdout, 
              "already inserted fine templates, skipping\n" ); 

          fcSegVec->data[i].level += 1;
        } 

        /* save a pointer to the last event in the list and count the events */
        ++numEvents;
        while ( eventList->next )
        {
          eventList = eventList->next;
          ++numEvents;
        }
        event = eventList;
        eventList = NULL;
      } /* end if ( events ) */

      /* if going up a level, remove inserted nodes, reset segment levels */ 
      if ( tmpltCurrent->next && (tmpltCurrent->next->tmpltPtr->level < 
            tmpltCurrent->tmpltPtr->level) )
      {
        /* record the current number */
        currentLevel = tmpltCurrent->tmpltPtr->level;

        /* decrease segment filter levels if the have been increased */
        for ( i = 0 ; i < fcSegVec->length; i++ )
        {
          if ( fcSegVec->data[i].level == currentLevel )
          {
            fcSegVec->data[i].level -= 1;
          }
        }

        if ( vrbflg ) fprintf( stdout, "removing inserted fine templates\n" ); 

        while ( tmpltCurrent->tmpltPtr->level == currentLevel )
        {
          LAL_CALL( LALFindChirpDestroyTmpltNode( &status, &tmpltCurrent ),
              &status );
        }          
      } /* end if up a level */

    } /* end loop over data segments */

  } /* end loop over linked list */

  /* save the number of events in the search summary table */
  searchsumm.searchSummaryTable->nevents = numEvents;


  /*
   *
   * free memory used by filtering code
   *
   */


  if ( vrbflg ) fprintf( stdout, "freeing memory\n" );

  /* free memory used by findchirp */
  LAL_CALL( LALFindChirpChisqVetoFinalize( &status, 
        fcFilterParams->chisqParams, fcInitParams->numChisqBins ), 
      &status );
  LAL_CALL( LALDestroyFindChirpInput( &status, &fcFilterInput ), 
      &status );
  LAL_CALL( LALFindChirpFilterFinalize( &status, &fcFilterParams ), 
      &status );
  LAL_CALL( LALFindChirpSPTemplateFinalize( &status, &fcTmpltParams ), 
      &status );
  LAL_CALL( LALFindChirpSPDataFinalize( &status, &fcDataParams ),
      &status );
  LAL_CALL( LALDestroyFindChirpSegmentVector( &status, &fcSegVec ),
      &status );
  LALFree( fcInitParams );

  /* free the template bank */
  while ( bankHead )
  {
    bankCurrent = bankHead;
    bankHead = bankHead->next;
    LALFree( bankCurrent );
    bankCurrent = NULL;
  }
  /* destroy linked list of template nodes */
  while ( tmpltHead )
  {
    LAL_CALL( LALFindChirpDestroyTmpltNode( &status, &tmpltHead ), &status );
  }

  /* free the data storage */
  LAL_CALL( LALFinalizeDataSegmentVector( &status, &dataSegVec ), &status );
  LAL_CALL( LALSDestroyVector( &status, &(chan.data) ), &status );
  LAL_CALL( LALSDestroyVector( &status, &(spec.data) ), &status );
  LAL_CALL( LALCDestroyVector( &status, &(resp.data) ), &status );


  /*
   *
   * write the result results to disk
   *
   */


cleanexit:

  if ( vrbflg ) fprintf( stdout, "writing result data to disk\n" );

  /* write the output frame */
  if ( writeRawData || writeFilterData || writeResponse || writeSpectrum ||
      writeRhosq || writeChisq )
  {
    snprintf( fname, sizeof(fname), "%s-INSPIRAL-%d-%d.gwf",
        site, gpsStartTime.gpsSeconds,
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
    frOutFile = FrFileONew( fname, 0 );
    FrameWrite( outFrame, frOutFile );
    FrFileOEnd( frOutFile );
  }

  /* open the output xml file */
  memset( &results, 0, sizeof(LIGOLwXMLStream) );
  snprintf( fname, sizeof(fname), "%s-INSPIRAL-%d-%d.xml",
      site, gpsStartTime.gpsSeconds,
      gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &results, fname), &status );

  /* write the process table */
  snprintf( proctable.processTable->ifos, LIGOMETA_IFOS_MAX, "%s", ifo );
  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->end_time),
        &accuracy ), &status );
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, process_table ), 
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, proctable, 
        process_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  LALFree( proctable.processTable );

  /* write the process params table */
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, process_params_table ), 
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, procparams, 
        process_params_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  while( procparams.processParamsTable )
  {
    this_proc_param = procparams.processParamsTable;
    procparams.processParamsTable = this_proc_param->next;
    LALFree( this_proc_param );
  }
  
  /* write the search summary table */
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, search_summary_table ), 
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, searchsumm, 
        search_summary_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  LALFree( searchsumm.searchSummaryTable );

  /* write the search summvars table */
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, search_summvars_table ), 
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, searchsummvars, 
        search_summvars_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  while( searchsummvars.searchSummvarsTable )
  {
    this_search_summvar = searchsummvars.searchSummvarsTable;
    searchsummvars.searchSummvarsTable = this_search_summvar->next;
    LALFree( this_search_summvar );
  }

  /* write the summvalue table */
  snprintf( summvalue.summValueTable->program, LIGOMETA_PROGRAM_MAX, 
      "%s", PROGRAM_NAME );
  summvalue.summValueTable->version = 0;
  summvalue.summValueTable->start_time = gpsStartTime;
  summvalue.summValueTable->end_time = gpsEndTime;
  snprintf( summvalue.summValueTable->ifo, LIGOMETA_IFO_MAX, "%s", ifo );
  snprintf( summvalue.summValueTable->name, LIGOMETA_SUMMVALUE_NAME_MAX, 
      "%s", "inspiral_effective_distance" );
  snprintf( summvalue.summValueTable->comment, LIGOMETA_SUMMVALUE_COMM_MAX, 
      "%s", "1.4_1.4_8" );
  summvalue.summValueTable->value = candle.effDistance;
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, summ_value_table ), 
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, summvalue, 
        summ_value_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );

  /* write the inspiral events to the file */
  if ( savedEvents.snglInspiralTable )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, sngl_inspiral_table ), 
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, savedEvents, 
          sngl_inspiral_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  }
  while ( savedEvents.snglInspiralTable )
  {
    event = savedEvents.snglInspiralTable;
    savedEvents.snglInspiralTable = savedEvents.snglInspiralTable->next;
    LALFree( event );
  }

  /* close the output xml file */
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &results ), &status );

  /* free the rest of the memory, check for memory leaks and exit */
  if ( injectionFile ) LALFree( injectionFile ); 
  if ( calCacheName ) LALFree( calCacheName );
  if ( frInCacheName ) LALFree( frInCacheName );
  if ( bankFileName ) LALFree( bankFileName );
  if ( channelName ) LALFree( channelName );
  if ( fqChanName ) LALFree( fqChanName );

  if ( vrbflg ) fprintf( stdout, "checking memory leaks and exiting\n" );
  LALCheckMemoryLeaks();
  exit( 0 );
}

/* ------------------------------------------------------------------------- */

#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
  LALCalloc( 1, sizeof(ProcessParamsTable) ); \
  snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
      PROGRAM_NAME ); \
      snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
          long_options[option_index].name ); \
          snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
          snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );


int arg_parse_check( int argc, char *argv[], MetadataTable procparams )
{
  /* getopt arguments */
  struct option long_options[] =
  {
    /* these options set a flag */
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"enable-event-cluster",    no_argument,       &eventCluster,     1 },
    {"disable-event-cluster",   no_argument,       &eventCluster,     0 },
    {"enable-output",           no_argument,       &enableOutput,     1 },
    {"disable-output",          no_argument,       &enableOutput,     0 },
    {"disable-high-pass",       no_argument,       &highPass,         0 },
    /* these options don't set a flag */
    {"gps-start-time",          required_argument, 0,                'a'},
    {"gps-start-time-ns",       required_argument, 0,                'A'},
    {"gps-end-time",            required_argument, 0,                'b'},
    {"gps-end-time-ns",         required_argument, 0,                'B'},
    {"channel-name",            required_argument, 0,                'c'},
    {"segment-length",          required_argument, 0,                'd'},
    {"number-of-segments",      required_argument, 0,                'e'},
    {"segment-overlap",         required_argument, 0,                'f'},
    {"sample-rate",             required_argument, 0,                'g'},
    {"help",                    no_argument,       0,                'h'},
    {"low-frequency-cutoff",    required_argument, 0,                'i'},
    {"spectrum-type",           required_argument, 0,                'j'},
    {"inverse-spec-length",     required_argument, 0,                'k'},
    {"dynamic-range-exponent",  required_argument, 0,                'l'},
    {"start-template",          required_argument, 0,                'm'},
    {"minimal-match",           required_argument, 0,                'M'},
    {"stop-template",           required_argument, 0,                'n'},
    {"chisq-bins",              required_argument, 0,                'o'},
    {"calibration-cache",       required_argument, 0,                'p'},
    {"snr-threshold",           required_argument, 0,                'q'},
    {"chisq-threshold",         required_argument, 0,                'r'},
    {"resample-filter",         required_argument, 0,                'R'},
    {"comment",                 required_argument, 0,                's'},
    {"enable-high-pass",        required_argument, 0,                't'},
    {"frame-cache",             required_argument, 0,                'u'},
    {"bank-file",               required_argument, 0,                'v'},
    {"injection-file",          required_argument, 0,                'w'},
    {"pad-data",                required_argument, 0,                'x'},
    {"debug-level",             required_argument, 0,                'z'},
    {"user-tag",                required_argument, 0,                'Z'},
    {"userTag",                 required_argument, 0,                'Z'},
    /* frame writing options */
    {"write-raw-data",          no_argument,       &writeRawData,     1 },
    {"write-filter-data",       no_argument,       &writeFilterData,  1 },
    {"write-response",          no_argument,       &writeResponse,    1 },
    {"write-spectrum",          no_argument,       &writeSpectrum,    1 },
    {"write-snrsq",             no_argument,       &writeRhosq,       1 },
    {"write-chisq",             no_argument,       &writeChisq,       1 },
    {0, 0, 0, 0}
  };
  int c;
  INT4 haveDynRange = 0;
  ProcessParamsTable *this_proc_param = procparams.processParamsTable;
  LALStatus             status = blank_status;


  /*
   *
   * parse command line arguments
   *
   */


  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;

    c = getopt_long_only( argc, argv, 
        "a:A:b:B:c:d:e:f:g:h:i:j:k:l:m:M:n:o:p:q:r:R:s:t:u:v:w:x:z:Z:", 
        long_options, &option_index );

    /* detect the end of the options */
    if ( c == - 1 )
    {
      break;
    }

    switch ( c )
    {
      case 0:
        /* if this option set a flag, do nothing else now */
        if ( long_options[option_index].flag != 0 )
        {
          break;
        }
        else
        {
          fprintf( stderr, "error parsing option %s with argument %s\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        break;

      case 'a':
        {
          long int gstartt = atol( optarg );
          if ( gstartt < 441417609 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS start time is prior to " 
                "Jan 01, 1994  00:00:00 UTC:\n"
                "(%ld specified)\n",
                long_options[option_index].name, gstartt );
            exit( 1 );
          }
          if ( gstartt > 999999999 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS start time is after " 
                "Sep 14, 2011  01:46:26 UTC:\n"
                "(%ld specified)\n", 
                long_options[option_index].name, gstartt );
            exit( 1 );
          }
          gpsStartTimeNS += (UINT8) gstartt * 1000000000LL;
          ADD_PROCESS_PARAM( "int", "%ld", gstartt );
        }
        break;

      case 'A':
        {
          long int gstarttns = atol( optarg );
          if ( gstarttns < 0 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS start time nanoseconds is negative\n",
                long_options[option_index].name );
            exit( 1 );
          }
          if ( gstarttns > 999999999 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS start time nanoseconds is greater than unity:\n" 
                "Must be <= 999999999 (%ld specified)\n", 
                long_options[option_index].name, gstarttns );
            exit( 1 );
          }
          gpsStartTimeNS += (UINT8) gstarttns;
          ADD_PROCESS_PARAM( "int", "%ld", gstarttns );
        }
        break;

      case 'b':
        {
          long int gendt = atol( optarg );
          if ( gendt > 999999999 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS end time is after " 
                "Sep 14, 2011  01:46:26 UTC:\n"
                "(%ld specified)\n", 
                long_options[option_index].name, gendt );
            exit( 1 );
          }
          else if ( gendt < 441417609 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS end time is prior to " 
                "Jan 01, 1994  00:00:00 UTC:\n"
                "(%ld specified)\n", 
                long_options[option_index].name, gendt );
            exit( 1 );
          }            
          gpsEndTimeNS += (UINT8) gendt * 1000000000LL;
          ADD_PROCESS_PARAM( "int", "%ld", gendt );
        }
        break;

      case 'B':
        {
          long int gendtns = atol( optarg );
          if ( gendtns < 0 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS end time nanoseconds is negative\n",
                long_options[option_index].name );
            exit( 1 );
          }
          else if ( gendtns > 999999999 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS end time nanoseconds is greater than unity:\n" 
                "Must be <= 999999999:\n"
                "(%ld specified)\n", 
                long_options[option_index].name, gendtns );
            exit( 1 );
          }            
          gpsEndTimeNS += (UINT8) gendtns;
          ADD_PROCESS_PARAM( "int", "%ld", gendtns );
        }
        break;

      case 'c':
        {
          /* create storage for the channel name and copy it */
          char *channamptr = NULL;
          size_t chanlen = strlen( optarg ) + 1;
          fqChanName = (CHAR *) LALCalloc( chanlen, sizeof(CHAR) );
          memcpy( fqChanName, optarg, chanlen );
          ADD_PROCESS_PARAM( "string", "%s", optarg );

          /* check that we have a proper channel name */
          if ( ! (channamptr = strstr( fqChanName, ":" ) ) )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "channel name must be a full LIGO channel name "
                "e.g. L1:LSC-AS_Q\n(%s specified)\n",
                long_options[option_index].name, optarg );
            exit( 1 );
          }
          chanlen = strlen( ++channamptr ) + 1;
          channelName = (CHAR *) LALCalloc( chanlen, sizeof(CHAR) );
          memcpy( channelName, channamptr, chanlen );

          /* copy the first character to site and the first two to ifo */
          memset( site, 0, sizeof(site) );
          memset( ifo, 0, sizeof(ifo) );
          memcpy( site, optarg, sizeof(site) - 1 );
          memcpy( ifo, optarg, sizeof(ifo) - 1 );
        }
        break;

      case 'd':
        numPoints = (INT4) atoi( optarg );
        if ( numPoints < 2 || numPoints % 2 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "number of points must be a non-zero power of 2: "
              "(%d specified) \n", 
              long_options[option_index].name, numPoints );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", numPoints );
        break;

      case 'e':
        numSegments = (INT4) atoi( optarg );
        if ( numSegments < 1 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "number of data segment must be greater than 0: "
              "(%d specified)\n", 
              long_options[option_index].name, numSegments );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", numSegments );
        break;

      case 'f':
        ovrlap = (INT4) atoi( optarg );
        if ( ovrlap < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "data segment overlap must be positive: "
              "(%d specified)\n", 
              long_options[option_index].name, ovrlap );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", ovrlap );
        break;

      case 'g':
        sampleRate = (INT4) atoi( optarg );
        if ( sampleRate < 2 || sampleRate > 16384 || sampleRate % 2 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "rate must be power of 2 between 2 and 16384 inclusive: "
              "(%d specified)\n", 
              long_options[option_index].name, sampleRate );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", sampleRate );
        break;

      case 'h':
        fprintf( stdout, USAGE );
        exit( 0 );
        break;

      case 'i':
        fLow = (REAL4) atof( optarg );
        if ( fLow < 40 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "low frequency cutoff is less than 40 Hz: "
              "(%f Hz specified)\n",
              long_options[option_index].name, fLow );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", fLow );
        break;

      case 'j':
        if ( ! strcmp( "mean", optarg ) )
        {
          specType = 0;
        }
        else if ( ! strcmp( "median", optarg ) )
        {
          specType = 1;
        }
        else if ( ! strcmp( "bad-mean", optarg ) )
        {
          specType = 0;
          badMeanPsd = 1;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown power spectrum type: "
              "%s (must be mean or median)\n", 
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'k':
        invSpecTrunc = (INT4) atoi( optarg );
        if ( invSpecTrunc < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "inverse spectrum length must be positive or zero: "
              "(%d specified)\n", 
              long_options[option_index].name, invSpecTrunc );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", invSpecTrunc );
        break;

      case 'l':
        dynRangeExponent = (REAL4) atof( optarg );
        haveDynRange = 1;
        ADD_PROCESS_PARAM( "float", "%e", dynRangeExponent );
        break;

      case 'm':
        startTemplate = (INT4) atoi( optarg );
        if ( startTemplate < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "template bank start index must be positive: "
              "(%d specified)\n", 
              long_options[option_index].name, startTemplate );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", startTemplate );
        break;

      case 'M':
        /* minimal match option is ignored */
        break;

      case 'n':
        stopTemplate = (INT4) atoi( optarg );
        if ( stopTemplate < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "template bank stop index must be positive: "
              "(%d specified)\n", 
              long_options[option_index].name, stopTemplate );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", stopTemplate );
        break;

      case 'o':
        numChisqBins = (INT4) atoi( optarg );
        if ( numChisqBins < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "number of chisq veto bins must be positive: "
              "(%d specified)\n", 
              long_options[option_index].name, numChisqBins );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", numChisqBins );
        break;

      case 'p':
        {
          /* create storage for the calibration frame cache name */
          size_t ccnamelen = strlen(optarg) + 1;
          calCacheName = (CHAR *) LALCalloc( ccnamelen, sizeof(CHAR));
          memcpy( calCacheName, optarg, ccnamelen );
          ADD_PROCESS_PARAM( "string", "%s", optarg );
        }
        break;

      case 'q':
        snrThresh = atof( optarg );
        if ( snrThresh < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "signal to noise threshold must be positive: "
              "(%f specified)\n", 
              long_options[option_index].name, snrThresh );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'r':
        chisqThresh = atof( optarg );
        if ( chisqThresh < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "chi squared threshold must be positive: "
              "(%f specified)\n", 
              long_options[option_index].name, chisqThresh );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'R':
        if ( ! strcmp( "ldas", optarg ) )
        {
          resampFiltType = 0;
        }
        else if ( ! strcmp( "butterworth", optarg ) )
        {
          resampFiltType = 1;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown resampling filter type: "
              "%s (must be ldas or butterworth)\n", 
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;


      case 's':
        if ( strlen( optarg ) > LIGOMETA_COMMENT_MAX - 1 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "comment must be less than %d characters\n",
              long_options[option_index].name, LIGOMETA_COMMENT_MAX );
          exit( 1 );
        }
        else
        {
          snprintf( comment, LIGOMETA_COMMENT_MAX, "%s", optarg);
        }
        break;

      case 't':
        highPass = 1;
        highPassFreq = (REAL4) atof( optarg );
        if ( highPassFreq < 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "low frequency cutoff is less than 0 Hz: "
              "(%f Hz specified)\n",
              long_options[option_index].name, highPassFreq );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", highPassFreq );
        break;

      case 'u':
        {
          /* create storage for the input frame cache name */
          size_t frcnamelen = strlen(optarg) + 1;
          frInCacheName = (CHAR *) LALCalloc( frcnamelen, sizeof(CHAR) );
          memcpy( frInCacheName, optarg, frcnamelen );
          ADD_PROCESS_PARAM( "string", "%s", optarg );
        }
        break;

      case 'v':
        {
          /* create storage for the calibration frame cache name */
          size_t bfnamelen = strlen(optarg) + 1;
          bankFileName = (CHAR *) LALCalloc( bfnamelen, sizeof(CHAR));
          memcpy( bankFileName, optarg, bfnamelen );
          ADD_PROCESS_PARAM( "string", "%s", optarg );
        }
        break;

      case 'w':
        {
          /* create storage for the injection file name */
          size_t ifnamelen = strlen(optarg) + 1;
          injectionFile = (CHAR *) LALCalloc( ifnamelen, sizeof(CHAR));
          memcpy( injectionFile, optarg, ifnamelen );
          ADD_PROCESS_PARAM( "string", "%s", optarg );
        }
        break;

      case 'x':
        padData = (UINT4) atoi( optarg );
        if ( padData < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "number of seconds to pad from input data"
              "must be greater than 0: (%d specified)\n", 
              long_options[option_index].name, padData );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", padData );
        break;

      case 'z':
        set_debug_level( optarg );
        break;

      case 'Z':
        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          LALCalloc( 1, sizeof(ProcessParamsTable) );
        snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", 
            PROGRAM_NAME );
        snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "-userTag" );
        snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%s",
            optarg );
        break;

      case '?':
        exit( 1 );
        break;

      default:
        fprintf( stderr, "unknown error while parsing options\n" );
        exit( 1 );
    }
  }

  if ( optind < argc )
  {
    fprintf( stderr, "extraneous command line arguments:\n" );
    while ( optind < argc )
    {
      fprintf ( stderr, "%s\n", argv[optind++] );
    }
    exit( 1 );
  }

  /* enable output is stored in the first process param row */
  if ( enableOutput == 1 )
  {
    snprintf( procparams.processParamsTable->program, 
        LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME );
    snprintf( procparams.processParamsTable->param,
        LIGOMETA_PARAM_MAX, "--enable-output" );
    snprintf( procparams.processParamsTable->type, 
        LIGOMETA_TYPE_MAX, "int" );
    snprintf( procparams.processParamsTable->value, 
        LIGOMETA_TYPE_MAX, "1" );
  }
  else if ( enableOutput == 0 )
  {
    snprintf( procparams.processParamsTable->program, 
        LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME );
    snprintf( procparams.processParamsTable->param,
        LIGOMETA_PARAM_MAX, "--disable-output" );
    snprintf( procparams.processParamsTable->type, 
        LIGOMETA_TYPE_MAX, "int" );
    snprintf( procparams.processParamsTable->value, 
        LIGOMETA_TYPE_MAX, "1" );
  }
  else
  {
    fprintf( stderr, "--enable-output or --disable-output "
        "argument must be specified\n" );
    exit( 1 );
  }


  /* check event cluster option */
  this_proc_param = this_proc_param->next = (ProcessParamsTable *)
    LALCalloc( 1, sizeof(ProcessParamsTable) );
  if ( eventCluster == 1 )
  {
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--enable-event-cluster" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "int" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, "1" );
  }
  else if ( eventCluster == 0 )
  {
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--disable-event-cluster" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "int" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, "1" );
  }
  else
  {
    fprintf( stderr, "--enable-event-cluster or "
        "--disable-event-cluster argument must be specified\n" );
    exit( 1 );
  }


  /*
   *
   * check validity of arguments
   *
   */


  /* check validity of input data time */
  if ( ! gpsStartTimeNS )
  {
    fprintf( stderr, "--gps-start-time must be specified\n" );
    exit( 1 );
  }
  LAL_CALL( LALINT8toGPS( &status, &gpsStartTime, &gpsStartTimeNS ), 
      &status );
  if ( ! gpsEndTimeNS )
  {
    fprintf( stderr, "--gps-end-time must be specified\n" );
    exit( 1 );
  }
  LAL_CALL( LALINT8toGPS( &status, &gpsEndTime, &gpsEndTimeNS ), 
      &status );
  if ( gpsEndTimeNS <= gpsStartTimeNS )
  {
    fprintf( stderr, "invalid gps time range: "
        "start time: %d, end time %d\n",
        gpsStartTime.gpsSeconds, gpsEndTime.gpsSeconds );
    exit( 1 );
  }

  /* check validity of data length parameters */
  if ( numPoints < 0 )
  {
    fprintf( stderr, "--segment-length must be specified\n" );
    exit( 1 );
  }
  if ( numSegments < 0 )
  {
    fprintf( stderr, "--number-of-segments must be specified\n" );
    exit( 1 );
  }
  if ( ovrlap < 0 )
  {
    fprintf( stderr, "--segment-overlap must be specified\n" );
    exit( 1 );
  }

  /* check sample rate has been given */
  if ( sampleRate < 0 )
  {
    fprintf( stderr, "--sample-rate must be specified\n" );
    exit( 1 );
  }

  /* check high pass option has been given */
  if ( highPass < 0 )
  {
    fprintf( stderr, "--disable-high-pass or --enable-high-pass (freq)"
        " must be specified\n" );
    exit( 1 );
  }

  /* check validity of input data length */
  inputDataLength = numPoints * numSegments - ( numSegments - 1 ) * ovrlap;
  {
    UINT8 gpsChanIntervalNS = gpsEndTimeNS - gpsStartTimeNS;
    UINT8 inputDataLengthNS = (UINT8) inputDataLength * 1000000000LL / 
      (UINT8) sampleRate;

    if ( inputDataLengthNS != gpsChanIntervalNS )
    {
      fprintf( stderr, "length of input data and data chunk do not match\n" );
      fprintf( stderr, "start time: %lld, end time %lld\n",
          gpsStartTimeNS / 1000000000LL, gpsEndTimeNS / 1000000000LL );
      fprintf( stderr, "gps channel time interval: %lld ns\n"
          "computed input data length: %lld ns\n", 
          gpsChanIntervalNS, inputDataLengthNS );
      exit( 1 );
    }
  }

  /* check filter parameters have been specified */
  if ( numChisqBins < 0 )
  {
    fprintf( stderr, "--chisq-bins must be specified\n" );
    exit( 1 );
  }
  if ( fLow < 0 )
  {
    fprintf( stderr, "--low-frequency-cutoff must be specified\n" );
    exit( 1 );
  }
  if ( resampFiltType < 0 )
  {
    fprintf( stderr, "--resample-filter must be specified\n" );
    exit( 1 );
  }
  if ( specType < 0 )
  {
    fprintf( stderr, "--spectrum-type must be specified\n" );
    exit( 1 );
  }
  if ( invSpecTrunc < 0 )
  {
    fprintf( stderr, "--inverse-spec-length must be specified\n" );
    exit( 1 );
  }
  else if ( invSpecTrunc * sampleRate > numPoints )
  {
    fprintf( stderr, "--inverse-spec-length must be less than "
        "--segment-length\n" );
    exit( 1 );
  }

  if ( ! haveDynRange )
  {
    fprintf( stderr, "--dynamic-range-exponent must be specified\n" );
    exit( 1 );
  }

  /* check that a channel has been requested and fill the ifo and site */
  if ( ! fqChanName )
  {
    fprintf( stderr, "--channel-name must be specified\n" );
    exit( 1 );
  }

  /* check that the thresholds have been specified */
  if ( snrThresh < 0 )
  {
    fprintf( stderr, "--snr-threshold must be specified\n" );
    exit( 1 );
  }
  if ( chisqThresh < 0 )
  {
    fprintf( stderr, "--chisq-threshold must be specified\n" );
    exit( 1 );
  }

  /* check that the frame caches have been specified */
  if ( ! frInCacheName )
  {
    fprintf( stderr, "--frame-cache must be specified\n" );
    exit( 1 );
  }
  if ( ! calCacheName )
  {
    fprintf( stderr, "--calibration-cache must be specified\n" );
    exit( 1 );
  }
  if ( ! bankFileName )
  {
    fprintf( stderr, "--bank-file must be specified\n" );
    exit( 1 );
  }

  return 0;
}

#undef ADD_PROCESS_PARAM

/* --- function to graph an array of REAL4 ------------------------------ */
#if 1
void
graphREAL4 (
    REAL4      *array, 
    INT4        n,
    INT4        spacing
    ) 
{
  FILE *fp;
  INT4 i;

  /* open a file for writing */
  if ( !(fp = fopen( "temp.graph", "w" )) )
  {
    printf( "couldn't open file\n" );
  }

  /* print data into the file */
  for ( i = 0; i < n; ++i )
    fprintf( fp, "%d\t%e\n", i, array[i * spacing] );

  /* close the file */
  fclose( fp );

  return;
}
#endif
