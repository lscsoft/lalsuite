/*----------------------------------------------------------------------- 
 * 
 * File Name: coherent_inspiral.c
 *
 * Author: Bose, S., Brown, D. A., and Noel, J. S.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <regex.h>
#include <time.h>
#include <math.h>

#include <FrameL.h>

#include <lalapps.h>
#include <series.h>
#include <processtable.h>
#include <lalappsfrutils.h>

#include <lal/LALConfig.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/LALConstants.h>
#include <lal/FrameStream.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/FrameCalibration.h>
#include <lal/Window.h>
#include <lal/TimeFreqFFT.h>
#include <lal/IIRFilter.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/FindChirpChisq.h>
#include <lal/TwoInterfFindChirp.h>
#include <lal/FindChirpEngine.h>

#include "inspiral.h"


RCSID( "$Id$" );

#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "coherent_inspiral"

/* define the parameters for a 1.4,1.4 solar mass standard candle with snr 8 */
#define CANDLE_MASS1 1.4
#define CANDLE_MASS2 1.4
#define CANDLE_RHOSQ 64.0

#define USAGE \
"lalapps_coherent_inspiral [options]\n\n"\
"See lalapps documentation for more information\n"\
"\n"


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
INT4  padData           = 0;            /* saftety margin on input data */
CHAR  *fqChanName[2]    = {NULL,NULL};  /* name of data channel         */
CHAR  *frInCacheName[2] = {NULL,NULL};  /* cache file containing frames */
INT4  numPoints         = -1;           /* points in a segment          */
INT4  numSegments       = -1;           /* number of segments           */
INT4  ovrlap            = -1;           /* overlap between segments     */
CHAR  ifo[2][3];                        /* ifo code for two detectors   */
UINT4 numDetectors      = 2;            /* hard-wired to TWO            */
CHAR *channelName[2] = {NULL,NULL};     /* channel string               */
UINT4 inputDataLength = 0;              /* number of points in input    */
REAL4 minimalMatch = -1;                /* override bank minimal match  */

/* data conditioning parameters */
INT4   resampFiltType   = -1;           /* low pass filter used for res */
INT4   sampleRate       = -1;           /* sample rate of filter data   */
INT4   highPass         = -1;           /* enable high pass on raw data */
REAL4  highPassFreq[2]  = {0,0};        /* high pass frequency          */
REAL4  fLow[2]          = {-1,-1};      /* low frequency cutoff         */
INT4   specType         = -1;           /* use median or mean psd       */
INT4   badMeanPsd       = 0;            /* use a mean with no overlap   */
INT4   invSpecTrunc[2]  = {-1,-1};      /* length of inverse spec (s)   */
REAL4  dynRangeExponent = -1;           /* exponent of dynamic range    */
CHAR  *calCacheName[2]  = {NULL,NULL};  /* location of calibration data */
CHAR  *injectionFile    = NULL;         /* name of file containing injs */

/* matched filter parameters */
CHAR *bankFileName      = NULL;         /* name of input template bank  */
INT4  startTemplate     = -1;           /* index of first template      */
INT4  stopTemplate      = -1;           /* index of last template       */
INT4  numChisqBins      = -1;           /* number of chisq bins         */
REAL4 snrThresh[2]      = {-1,-1};      /* snr thresholds in each ifo   */
REAL4 coherentSnrThresh = -1;           /* coherent snr threshold       */
REAL4 chisqThresh[2]    = {-1,-1};      /* chisq veto thresholds        */
INT4  eventCluster      = -1;           /* perform chirplen clustering  */
UINT4 site0             = -1;           /* 1st ifo: lalCachedDetector site identifier*/
UINT4 site1             = -1;           /* 2nd ifo: lalCachedDetector site identifier*/

/* output parameters */
CHAR  *userTag          = NULL;         /* string the user can tag with */
int    enableOutput     = -1;           /* write out inspiral events    */
int    writeRawData     = 0;            /* write the raw data to a file */
int    writeFilterData  = 0;            /* write post injection data    */
int    writeResponse    = 0;            /* write response function used */
int    writeSpectrum    = 0;            /* write computed psd to file   */
int    writeRhosq       = 0;            /* write rhosq time series      */
int    writeCoherentRhosq = 0;          /* write network rhosq ts       */
int    writeChisq       = 0;            /* write chisq time series      */

/* other command line args */
CHAR comment[LIGOMETA_COMMENT_MAX];     /* process param comment        */

int main( int argc, char *argv[] )
{
  /* lal function variables */
  LALStatus             status = blank_status;
  LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;

  /* frame input data pair*/
  FrCache      *frInCache[2] = {NULL,NULL};
  FrStream     *frStream[2] = {NULL,NULL};
  FrChanIn      frChan[2];

  /* frame output data */
  struct FrFile *frOutFile  = NULL;
  struct FrameH *outFrame[2]= {NULL,NULL};
  struct FrameH *outFrameNet= NULL;
  UINT4          nRhosqFr = 0;
  UINT4          nChisqFr = 0;

  /* raw input data storage */
  REAL4TimeSeries               chan[2];
  REAL4FrequencySeries          spec[2];
  COMPLEX8FrequencySeries       resp[2];
  DataSegmentVector            *dataSegVec = NULL;
  DataSegmentVector            *dataVecPtr[2] = {NULL, NULL};
  
  /* structures for preconditioning */
#if 0
  COMPLEX8FrequencySeries       injResp;        
  COMPLEX8FrequencySeries      *injRespPtr;     
#endif
  ResampleTSParams              resampleParams; 
  LALWindowParams               wpars;
  AverageSpectrumParams         avgSpecParams[2];
  
  /* findchirp data structures */
  FindChirpInitParams          *fcInitParams[2]   = {NULL,NULL};
  FindChirpSPDataParams        *fcDataParams   = NULL;
  FindChirpSPTmpltParams       *fcTmpltParams[2]  = {NULL,NULL};
  FindChirpFilterParams        *fcFilterParams = NULL;
  FindChirpFilterInput         *fcFilterInput  = NULL;
  FindChirpStandardCandle       candle; 
  
  /* twointerffindchirp data structures */
  TwoInterfFindChirpInitParams          *twoInterfInitParams = NULL;
  TwoInterfFindChirpSegmentVector       *twoInterfFcSegVec = NULL;
  TwoInterfFindChirpFilterInputVector   *twoInterfFilterInputVec = NULL;
  TwoInterfDataSegmentVector            *twoInterfDataSegVec = NULL;
  TwoInterfFindChirpSPDataParamsVector  *twoInterfDataParamsVec = NULL;
  TwoInterfFindChirpFilterParams        *twoInterfFilterParams = NULL;
  /* XXX no event handling at present */
#if 0
  TwoInterfInspiralEvent                *twoInterfInspiralEvent = NULL;
#endif
  LALDetectorPair               detectors;
  
  /* inspiral template structures */
  INT4                          numTmplts    = 0;
  InspiralTemplate             *bankHead     = NULL;
  InspiralTemplate             *bankCurrent  = NULL;
  InspiralTemplateNode         *tmpltHead    = NULL;
  InspiralTemplateNode         *tmpltCurrent = NULL;
  
  /* inspiral events */
  INT4                          numEvents   = 0;
  
  TwoInterfInspiralEvent       *event       = NULL;
  TwoInterfInspiralEvent       *eventList   = NULL;
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
  UINT4 i, n;
  INT4  inserted;
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
    calloc( 1, sizeof(ProcessTable) );
  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->start_time),
        &accuracy ), &status );
  LAL_CALL( populate_process_table( &status, proctable.processTable, 
        PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE ), &status );
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *) 
    calloc( 1, sizeof(ProcessParamsTable) );
  memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );

  /* create the search summary table and summvars table*/
  searchsumm.searchSummaryTable = (SearchSummaryTable *)
    calloc( 1, sizeof(SearchSummaryTable) );
  searchsummvars.searchSummvarsTable = NULL;

  /* call the argument parse and check function */
  arg_parse_check( argc, argv, procparams );
  for ( this_proc_param = procparams.processParamsTable; this_proc_param->next;
      this_proc_param = this_proc_param->next );

  /* OK to use LALCalloc() from here... */

  /* fill the comment, if a user has specified on, or leave it blank */
  if ( ! *comment )
  {
    LALSnprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX, " " );
    LALSnprintf( searchsumm.searchSummaryTable->comment, LIGOMETA_COMMENT_MAX, 
        " " );
  }
  else
  {
    LALSnprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX,
        "%s", comment );
    LALSnprintf( searchsumm.searchSummaryTable->comment, LIGOMETA_COMMENT_MAX,
        "%s", comment );
  }

  /* the number of nodes for a standalone job is always 1 */
  searchsumm.searchSummaryTable->nnodes = 1;

  /* make sure the pointer to the first event is null */
  savedEvents.multiInspiralTable = NULL;


  /* create the standard candle and database table */
  summvalue.summValueTable = &candleTable;
  memset( &candleTable, 0, sizeof(SummValueTable) );
  memset( &candle, 0, sizeof(FindChirpStandardCandle) );
  strncpy( candle.ifo, ifo[0], 2 * sizeof(CHAR) );
  candle.tmplt.mass1 = CANDLE_MASS1;
  candle.tmplt.mass2 = CANDLE_MASS2;
  candle.rhosq       = CANDLE_RHOSQ;
  candle.tmplt.totalMass = candle.tmplt.mass1 + candle.tmplt.mass2;
  candle.tmplt.mu = candle.tmplt.mass1 * candle.tmplt.mass2 / 
    candle.tmplt.totalMass;
  candle.tmplt.eta = candle.tmplt.mu / candle.tmplt.totalMass;


  /*
   *
   * create a template bank
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

  /* override the minimal match of the bank if specified on the command line */
  if ( minimalMatch >= 0 )
  {
    if ( vrbflg )
    {
      fprintf( stdout, "Overriding bank minimal match:\n   value in bank = %e,"
          " new value = %e\n", bankHead->minMatch, minimalMatch );
    }
    for ( bankCurrent = bankHead; bankCurrent; bankCurrent = bankCurrent->next )
    {
      bankCurrent->minMatch = minimalMatch;
    }
  }

  /* save the minimal match of the bank in the process params */
  this_proc_param = this_proc_param->next = (ProcessParamsTable *) 
    calloc( 1, sizeof(ProcessParamsTable) ); 
  LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", 
      PROGRAM_NAME );
  LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--minimal-match" );
  LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "float" ); 
  LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%e", 
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
  for ( n = 0 ; n < numDetectors ; ++n ) 
  { 
    memset( &chan[n], 0, sizeof(REAL4TimeSeries) );
    chan[n].epoch = gpsStartTime;
    chan[n].epoch.gpsSeconds -= padData; /* subtract pad seconds from start */

    /* open a frame cache or files, seek to required epoch and set chan name */
    if ( frInCacheName[n] ) 
    {
      LAL_CALL( LALFrCacheImport( &status, &frInCache[n], frInCacheName[n]), 
          &status );
      LAL_CALL( LALFrCacheOpen( &status, &frStream[n], frInCache[n] ), 
          &status );
    }
    else 
    {
      LAL_CALL( LALFrOpen( &status, &frStream[n], NULL, "*.gwf" ), &status );
    }
    LAL_CALL( LALFrSeek( &status, &(chan[n].epoch), frStream[n] ), &status );
    frChan[n].name = fqChanName[n];

    /* determine the sample rate of the raw data and allocate enough memory */
    LAL_CALL( LALFrGetREAL4TimeSeries( &status, &chan[n], &frChan[n], 
          frStream[n] ), &status );

    /* store the input sample rate */
    this_search_summvar = searchsummvars.searchSummvarsTable = 
      (SearchSummvarsTable *) calloc( 1, sizeof(SearchSummvarsTable) );
    LALSnprintf( this_search_summvar->name, LIGOMETA_NAME_MAX * sizeof(CHAR),
        "raw data sample rate" );
    this_search_summvar->value = chan[n].deltaT;

    /* determine if we need to resample the channel */
    if ( vrbflg )
    {
      fprintf( stdout, "resampleParams.deltaT = %e\n", resampleParams.deltaT );
      fprintf( stdout, "chan.deltaT = %e\n", chan[n].deltaT );
    }
    if ( ! ( fabs( resampleParams.deltaT - chan[n].deltaT ) < epsilon ) )
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

    /* determine the number of points to get and create storage for the data */
    inputLengthNS = 
      (REAL8) ( gpsEndTimeNS - gpsStartTimeNS + 2000000000LL * padData );
    numInputPoints = 
      (UINT4) floor( inputLengthNS / (chan[n].deltaT * 1.0e9) + 0.5 );
    LAL_CALL( LALSCreateVector( &status, &(chan[n].data), numInputPoints ), 
        &status );

    if ( vrbflg ) fprintf( stdout, "input channel %s has sample interval "
        "(deltaT) = %e\nreading %d points from frame stream\n", 
        fqChanName[n], chan[n].deltaT, numInputPoints );

    /* read the data channel time series from frames */
    LAL_CALL( LALFrGetREAL4TimeSeries( &status, &chan[n], &frChan[n], 
          frStream[n] ), &status );
    memcpy( &(chan[n].sampleUnits), &lalADCCountUnit, sizeof(LALUnit) );

    /* store the start and end time of the raw channel in the search summary */
    searchsumm.searchSummaryTable->in_start_time = chan[n].epoch;
    LAL_CALL( LALGPStoFloat( &status, &tsLength, &(chan[n].epoch) ), 
        &status );
    tsLength += chan[n].deltaT * (REAL8) chan[n].data->length;
    LAL_CALL( LALFloatToGPS( &status, 
          &(searchsumm.searchSummaryTable->in_end_time), &tsLength ), 
        &status );

    /* close the frame file stream and destroy the cache */
    LAL_CALL( LALFrClose( &status, &frStream[n] ), &status );
    if ( frInCacheName[n] )
    {
      LAL_CALL( LALDestroyFrCache( &status, &frInCache[n] ), &status );
    }
  }


  /* write the raw channel data as read in from the frame files */
  if ( writeRawData )
  {
    outFrame[0] = fr_add_proc_REAL4TimeSeries( outFrame[0], &chan[0], 
        "ct", "RAW1" );

    outFrame[1] = fr_add_proc_REAL4TimeSeries( outFrame[1], &chan[1], 
        "ct", "RAW2" );
  }

  if ( vrbflg ) 
  {
    fprintf( stdout, "read channel %s from frame stream 1\n"
        "got %d points with deltaT %e\nstarting at GPS time %d sec %d ns\n", 
        chan[0].name, chan[0].data->length, chan[0].deltaT, 
        chan[0].epoch.gpsSeconds, chan[0].epoch.gpsNanoSeconds );
    fprintf( stdout, "read channel %s from frame stream 2\n"
        "got %d points with deltaT %e\nstarting at GPS time %d sec %d ns\n", 
        chan[1].name, chan[1].data->length, chan[1].deltaT, 
        chan[1].epoch.gpsSeconds, chan[1].epoch.gpsNanoSeconds );
  }


  /*
   *
   * generate the response function for the requested time
   *
   */


  /* create storage for the response function */
  for ( n = 0 ; n < numDetectors ; ++n ) 
  { 
    memset( &resp[n], 0, sizeof(COMPLEX8FrequencySeries) );
    LAL_CALL( LALCCreateVector( &status, &(resp[n].data), numPoints / 2 + 1 ), 
        &status );

    /* set the parameters of the response to match the data */
    resp[n].epoch = gpsStartTime;
    resp[n].deltaF = (REAL8) sampleRate / (REAL8) numPoints;
    resp[n].sampleUnits = strainPerCount;
    strcpy( resp[n].name, chan[n].name );

    /* generate the response function for the current time */
    if ( vrbflg ) fprintf( stdout, "generating response at time %d sec %d ns\n",
        resp[n].epoch.gpsSeconds, resp[n].epoch.gpsNanoSeconds );
    LAL_CALL( LALExtractFrameResponse( &status, &resp[n], calCacheName[n], 
          ifo[n] ), &status );
  }

  if ( writeResponse )
  {
    outFrame[0] = fr_add_proc_COMPLEX8FrequencySeries( outFrame[0], 
        &resp[0], "strain/ct", "RESPONSE1" );
    outFrame[1] = fr_add_proc_COMPLEX8FrequencySeries( outFrame[1], 
        &resp[1], "strain/ct", "RESPONSE2" );
  }


  /*
   *
   * resample the data to the requested rate
   *
   */


  if ( resampleChan )
  {
    for ( n = 0 ; n < numDetectors ; ++n )
    { 
      if (vrbflg) fprintf( stdout, "resampling input data from %e to %e\n",
          chan[n].deltaT, resampleParams.deltaT );

      LAL_CALL( LALResampleREAL4TimeSeries( &status, &chan[n], 
            &resampleParams ), &status );

      if ( vrbflg ) fprintf( stdout, "channel %s resampled:\n"
          "%d points with deltaT %e\nstarting at GPS time %d sec %d ns\n", 
          chan[n].name, chan[n].data->length, chan[n].deltaT, 
          chan[n].epoch.gpsSeconds, chan[n].epoch.gpsNanoSeconds );
    }
  }

  /* write the resampled channel data as read in from the frame files */
  if ( writeRawData )
  {
    outFrame[0] = fr_add_proc_REAL4TimeSeries( outFrame[0], &chan[0], 
        "ct", "RAW_RESAMP1" );
    outFrame[1] = fr_add_proc_REAL4TimeSeries( outFrame[1], &chan[1], 
        "ct", "RAW_RESAMP2" );
  }

  /* store the filter data sample rate */
  this_search_summvar = this_search_summvar->next = 
    (SearchSummvarsTable *) calloc( 1, sizeof(SearchSummvarsTable) );
  LALSnprintf( this_search_summvar->name, LIGOMETA_NAME_MAX * sizeof(CHAR),
      "filter data sample rate" );
  this_search_summvar->value = chan[0].deltaT;


  /* 
   *
   * high pass the data, removed pad from time series and check length of data
   *
   */


  /* iir filter to remove low frequencies from data channel */
  if ( highPass )
  {
    if ( vrbflg )
    {
      fprintf( stdout, "highpassing detector 1 data above %e\n", 
          highPassFreq[0] );
      fprintf( stdout, "highpassing detector 2 data above %e\n", 
          highPassFreq[1] );
    }

    for ( n = 0; n < numDetectors ; ++n )
    {
      PassBandParamStruc highpassParam;
      highpassParam.nMax = 4;
      highpassParam.f1 = -1.0;
      highpassParam.f2 = highPassFreq[n];
      highpassParam.a1 = -1.0;
      highpassParam.a2 = 0.1;

      LAL_CALL( LALButterworthREAL4TimeSeries( &status, &chan[n], 
            &highpassParam ), &status );
    }
  }

  /* remove pad from requested data from start and end of time series */
  for ( n = 0 ; n < numDetectors ; ++n ) 
  {
    memmove( chan[n].data->data, chan[n].data->data + padData * sampleRate, 
        (chan[n].data->length - 2 * padData * sampleRate) * sizeof(REAL4) );
    LALRealloc( chan[n].data->data, 
        (chan[n].data->length - 2 * padData * sampleRate) * sizeof(REAL4) );
    chan[n].data->length -= 2 * padData * sampleRate;
    chan[n].epoch.gpsSeconds += padData;

    if ( vrbflg ) fprintf( stdout, "after removal of %d second padding at "
        "start and end:\ndata channel sample interval (deltaT) = %e\n"
        "data channel length = %d\nstarting at %d sec %d ns\n", 
        padData , chan[n].deltaT , chan[n].data->length, 
        chan[n].epoch.gpsSeconds, chan[n].epoch.gpsNanoSeconds );
  }

  if ( writeFilterData ) 
  {
    outFrame[0] = fr_add_proc_REAL4TimeSeries( outFrame[0], &chan[0],
        "ct", "FILTER1" );
    outFrame[1] = fr_add_proc_REAL4TimeSeries( outFrame[1], &chan[1],
        "ct", "FILTER2" );
  }

  /* check data length */
  for ( n = 0 ; n < numDetectors ; ++n )
  {
  
    if ( chan[n].data->length != inputDataLength )
    {
      fprintf( stderr, "error: computed channel length and requested\n"
          "input data length do not match:\nchan.data->length = %d\n"
          "inputDataLength = %d\nyou have found a bug in the code.\n"
          "please report this to <duncan@gravity.phys.uwm.edu>\n",
          chan[n].data->length, inputDataLength );
      exit( 1 );
    }

    /* store the start and end time of the filter channel in the search summ */
    /* noting that we don't look for events in the first and last quarter    */
    /* of each findchirp segment of the input data                           */
    LAL_CALL( LALGPStoFloat( &status, &tsLength, &(chan[n].epoch) ), 
        &status );
    tsLength += (REAL8) (numPoints / 4) * chan[n].deltaT;
    LAL_CALL( LALFloatToGPS( &status, 
          &(searchsumm.searchSummaryTable->out_start_time), 
          &tsLength ), &status );
    LAL_CALL( LALGPStoFloat( &status, &tsLength, &(chan[n].epoch) ), 
        &status );
    tsLength += chan[n].deltaT 
      * ((REAL8) chan[n].data->length - (REAL8) (numPoints/4));
    LAL_CALL( LALFloatToGPS( &status, 
          &(searchsumm.searchSummaryTable->out_end_time), 
          &tsLength ), &status );
  }


  /* 
   *
   * create and populate twointerffindchirp initialization structure 
   *
   */
  
  
  if ( ! (twoInterfInitParams = (TwoInterfFindChirpInitParams *) 
	  LALCalloc (1, sizeof(TwoInterfFindChirpInitParams))))
    {
      fprintf( stderr, 
	       "could not allocate memory for twointerffindchirp init params\n" );
      exit( 1 );
    }
  twoInterfInitParams->numDetectors            = numDetectors;
  twoInterfInitParams->numSegments             = numSegments;
  twoInterfInitParams->numPoints               = numPoints;
  twoInterfInitParams->numChisqBins            = numChisqBins;
  twoInterfInitParams->createRhosqVec          = writeRhosq;
  twoInterfInitParams->createTwoInterfRhosqVec = writeCoherentRhosq;
  twoInterfInitParams->ovrlap                  = ovrlap;

   /* initialize data conditioning routines */
  LAL_CALL( LALTwoInterfFindChirpSPDataInit(&status, &twoInterfDataParamsVec, 
	twoInterfInitParams), &status );
  
  fcDataParams = twoInterfDataParamsVec->data;  

  /* create the twointerf-findchirp data storage */
  LAL_CALL( LALCreateTwoInterfFindChirpSegmentVector(&status, 
        &twoInterfFcSegVec, twoInterfInitParams), &status );

  if ( ! (twoInterfDataSegVec = (TwoInterfDataSegmentVector *) 
	  LALCalloc( 1, sizeof(TwoInterfDataSegmentVector) ) ))
    {    
      fprintf( stderr, 
	       "could not allocate memory for twoInterfDataSegVec\n" );
      exit( 1 );

    }

  /* create the data segment vector */
  twoInterfDataSegVec->length = numDetectors;
  if ( ! (dataSegVec = twoInterfDataSegVec->data = (DataSegmentVector *) 
	  LALCalloc( 1, twoInterfDataSegVec->length * 
		     sizeof(TwoInterfDataSegmentVector) ) ))
    {    
      fprintf( stderr, 
	       "could not allocate memory for twoInterfDataSegVec->data\n" );
      exit( 1 );
    }  
  
  for ( n = 0 ; n < numDetectors ; ++n) {
    if ( ! ( dataSegVec[n].data = (DataSegment *) 
	     LALCalloc( 1, numSegments * sizeof(DataSegment) ) ) )
      {
	fprintf( stderr, 
		 "could not allocate memory for data segments\n" );
	exit( 1 );
      }
  }
  
  /*
   *
   * power spectrum estimation and data conditioning
   *
   */

  for ( n = 0 ; n < numDetectors ; ++n) {
    if ( ! ( fcInitParams[n] = (FindChirpInitParams *) 
	     LALCalloc( 1, sizeof(FindChirpInitParams) ) ) )
      {
	fprintf( stderr, 
		 "could not allocate memory for findchirp init params\n" );
	exit( 1 );
      }
    fcInitParams[n]->numPoints      = numPoints;
    fcInitParams[n]->numSegments    = numSegments;
    fcInitParams[n]->numChisqBins   = numChisqBins;
    fcInitParams[n]->createRhosqVec = writeRhosq;
    fcInitParams[n]->ovrlap         = ovrlap;
    
    /* create the data segment vector */
    memset( &spec[n], 0, sizeof(REAL4FrequencySeries) );
    LAL_CALL( LALSCreateVector( &status, &(spec[n].data), numPoints / 2 + 1 ), 
	      &status );

    LAL_CALL( LALInitializeDataSegmentVector( &status, &dataVecPtr[n],
					      &chan[n], &spec[n], 
					      &resp[n], fcInitParams[n] ), 
	      &status );
    
    
    /*
     *
     * power spectrum estimation and data conditioning
     *
     */
    
    fcDataParams[n].invSpecTrunc = (invSpecTrunc[n] * sampleRate);
    fcDataParams[n].fLow         = fLow[n];
    
    /* compute the windowed power spectrum for the data channel */
    avgSpecParams[n].window = NULL;
    avgSpecParams[n].plan   = fcDataParams[n].fwdPlan;
    switch ( specType )
      {
      case 0:
        avgSpecParams[n].method = useMean;
        if ( vrbflg ) fprintf( stdout, "computing mean psd" );
        break;
      case 1:
        avgSpecParams[n].method = useMedian;
        if ( vrbflg ) fprintf( stdout, "computing median psd" );
        break;
    }

    wpars.type = Hann;
    wpars.length = numPoints;
    if ( badMeanPsd )
    {
      avgSpecParams[n].overlap = 0;
      if ( vrbflg ) fprintf( stdout, " without overlap\n" );
    }
    else
    {
      avgSpecParams[n].overlap = numPoints / 2;
      if ( vrbflg ) 
        fprintf( stdout, " with overlap %d\n", avgSpecParams[n].overlap );
    }

    LAL_CALL( LALCreateREAL4Window( &status, &(avgSpecParams[n].window),
          &wpars ), &status );
    LAL_CALL( LALREAL4AverageSpectrum( &status, &spec[n], 
          &chan[n], &avgSpecParams[n] ),
        &status );
    LAL_CALL( LALDestroyREAL4Window( &status, &(avgSpecParams[n].window) ), 
        &status );
    strcpy( spec[n].name, chan[n].name );
  }
  
  /* write the spectrum data to a file */
  if ( writeSpectrum )
  {
    outFrame[0] = fr_add_proc_REAL4FrequencySeries(outFrame[0], &spec[0],
        "ct/sqrtHz", "PSD1");
    outFrame[1] = fr_add_proc_REAL4FrequencySeries(outFrame[1], &spec[1],
	"ct/sqrtHz", "PSD2");
  }  
  
  for ( n = 0 ; n < numDetectors ; ++n) {
    dataSegVec[n].length = dataVecPtr[n]->length;
    for ( i = 0 ; i < dataSegVec[n].length ; i++ ) {
      dataSegVec[n].data[i].chan = dataVecPtr[n]->data[i].chan;
      dataSegVec[n].data[i].spec = dataVecPtr[n]->data[i].spec;
      dataSegVec[n].data[i].resp = dataVecPtr[n]->data[i].resp;
      
    }
  }
    

  /*
   *
   * create the data structures needed for twointerffindchirp
   *
   */
  
  
  if ( vrbflg ) fprintf( stdout, "initializing twointerffindchirp\n" );
  
  /* initialize the template functions */
  for ( n = 0 ; n < numDetectors ; ++n) {
    LAL_CALL( LALFindChirpSPTemplateInit( &status, &fcTmpltParams[n], 
					  fcInitParams[n] ), &status );
    
    fcDataParams[n].dynRange = fcTmpltParams[n]->dynRange = 
      pow( 2.0, dynRangeExponent );
    fcDataParams[n].deltaT = 
      fcTmpltParams[n]->deltaT = 1.0 / (REAL4) sampleRate;
    fcTmpltParams[n]->fLow = fLow[n];
  }
  
  /* initialize findchirp filter functions */
  LAL_CALL( LALTwoInterfFindChirpFilterInit ( &status, 
        &twoInterfFilterParams, twoInterfInitParams), &status );

  detectors.detectorOne = lalCachedDetectors[site0];
  detectors.detectorTwo = lalCachedDetectors[site1];
  
  twoInterfFilterParams->twoInterfRhosqThresh = coherentSnrThresh;
  twoInterfFilterParams->detectors            = &detectors;
  
  LAL_CALL( LALCreateTwoInterfFindChirpInputVector (&status, 
       &twoInterfFilterInputVec, twoInterfInitParams), &status );
  
  
  fcFilterParams = twoInterfFilterParams->paramsVec->filterParams;
  
  for ( n = 0 ; n < numDetectors ; ++n)
    {
    fcFilterParams[n].deltaT = 1.0 / (REAL4) sampleRate;
    fcFilterParams[n].computeNegFreq = 0;
    LAL_CALL( LALTwoInterfFindChirpChisqVetoInit (&status, 
          fcFilterParams[n].chisqParams, 
          numChisqBins, numPoints),
        &status );

    /* parse the thresholds */
    fcFilterParams[n].rhosqThresh = snrThresh[n] * snrThresh[n];
    fcFilterParams[n].chisqThresh = chisqThresh[n];
    fcFilterParams[n].maximiseOverChirp = eventCluster;
  }


  /*
   *
   * condition data segments for filtering
   *
   */


  if ( vrbflg ) fprintf( stdout, "twointerffindchirp conditioning data\n" );
  LAL_CALL( LALTwoInterfFindChirpSPData (&status, twoInterfFcSegVec, 
					 twoInterfDataSegVec, twoInterfDataParamsVec),
	    &status );

  /* compute the standard candle */
  {
    REAL4 cannonDist = 1.0; /* Mpc */
    REAL4 m  = (REAL4) candle.tmplt.totalMass;
    REAL4 mu = (REAL4) candle.tmplt.mu;
    REAL4 distNorm = 2.0 * LAL_MRSUN_SI / (cannonDist * 1e6 * LAL_PC_SI);
    REAL4 candleTmpltNorm = sqrt( (5.0*mu) / 96.0 ) *
      pow( m / (LAL_PI*LAL_PI) , 1.0/3.0 ) *
      pow( LAL_MTSUN_SI / (REAL4) chan[0].deltaT, -1.0/6.0 );

    /*assumes same value for both detectors*/
    distNorm *= fcTmpltParams[0]->dynRange;
    candleTmpltNorm *= candleTmpltNorm;
    candleTmpltNorm *= distNorm * distNorm;

    candle.sigmasq = 4.0 * ( (REAL4) chan[0].deltaT / (REAL4) numPoints );

    /* candle sigmasq and eff. dist. computed for detector 1 data alone*/
    candle.sigmasq *= 
      candleTmpltNorm * twoInterfFcSegVec->data[0].data->segNorm;
    candle.effDistance = sqrt( candle.sigmasq / candle.rhosq );

    if ( vrbflg ) 
    {
      fprintf( stdout, "candle m = %e\ncandle mu = %e\n"
          "candle.rhosq = %e\nchan.deltaT = %e\n"
          "numPoints = %d\nfcSegVec->data->segNorm = %e\n"
          "candleTmpltNorm = %e\ncandle.effDistance = %e Mpc\n"
          "candle.sigmasq = %e\n",
          m, mu, candle.rhosq, chan[0].deltaT, numPoints, 
          twoInterfFcSegVec->data[0].data->segNorm, candleTmpltNorm, 
          candle.effDistance, candle.sigmasq );
      fflush( stdout );
    }
  }



  /*
   *
   * search engine
   *
   */

  fcFilterInput = twoInterfFilterInputVec->filterInput;
  
  for ( tmpltCurrent = tmpltHead, inserted = 0; tmpltCurrent; 
	tmpltCurrent = tmpltCurrent->next, inserted = 0 )
   {
    /* loop over detectors */
    for ( n = 0; n < numDetectors; ++n )
    {
      /*  generate template */
      LAL_CALL( LALFindChirpSPTemplate( &status, fcFilterInput[n].fcTmplt,
            tmpltCurrent->tmpltPtr, fcTmpltParams[n] ), 
          &status );
      fcFilterInput[n].tmplt = tmpltCurrent->tmpltPtr;

      /*
       *
       * loop over segments
       *
       */

      for ( i = 0; i < twoInterfInitParams->numSegments; ++i )
      {
        twoInterfFilterInputVec->filterInput[n].segment = 
          twoInterfFcSegVec->data[n].data + i;

      } /* end loop over number of segments */

    } /* end loop over detectors */


    LAL_CALL( LALTwoInterfFindChirpFilterSegment (&status, &eventList, 
          twoInterfFilterInputVec, twoInterfFilterParams), 
        &status );

    if ( writeRhosq )
    {
      for ( n = 0; n < numDetectors; ++n ) 
      {
        CHAR snrsqStr[LALNameLength];
        LALSnprintf( snrsqStr, LALNameLength*sizeof(CHAR), 
            "SNRSQ_%d", nRhosqFr++ );
        strcpy( fcFilterParams[n].rhosqVec->name, chan[n].name );
        outFrame[n] = fr_add_proc_REAL4TimeSeries(outFrame[n], 
            fcFilterParams[n].rhosqVec, 
            "none", snrsqStr );
      }
    }

    if ( writeCoherentRhosq )
    {
      CHAR snrsqStr[LALNameLength];
      LALSnprintf( snrsqStr, LALNameLength*sizeof(CHAR), 
          "SNRSQ_%d", nRhosqFr++ );
      outFrameNet = fr_add_proc_REAL4TimeSeries(outFrameNet, 
          twoInterfFilterParams->twoInterfRhosqVec, 
          "none", snrsqStr );
    }

    if ( writeChisq ) 
    {
      for ( n = 0; n < numDetectors; ++n )
      {
        CHAR chisqStr[LALNameLength];
        REAL4TimeSeries chisqts;
        LALSnprintf( chisqStr, LALNameLength*sizeof(CHAR), 
            "CHISQ_%d", nChisqFr++ );
        chisqts.epoch = fcFilterInput[n].segment->data->epoch;
        memcpy( &(chisqts.name), fcFilterInput[n].segment->data->name,
            LALNameLength * sizeof(CHAR) );
        chisqts.deltaT = fcFilterInput[n].segment->deltaT;
        chisqts.data = fcFilterParams[n].chisqVec;
        outFrame[n] = fr_add_proc_REAL4TimeSeries( outFrame[n], 
            &chisqts, "none", chisqStr );
      }
    }


    /*  test if filter returned any events */
    if ( eventList )
    {
      if ( vrbflg )
      {
        fprintf( stdout, "segment %d rang template %e,%e\n",
            twoInterfFcSegVec->data[0].data[i].number,
            fcFilterInput[0].tmplt->mass1, 
            fcFilterInput[0].tmplt->mass2 );
        fprintf( stdout, "***>  dumping events  <***\n" );
      }

        /* XXX event handling code will segfault as it deferences NULL ptr */
#if 0
      if ( eventList && ! savedEvents.multiInspiralTable ) {
	/* Eventually the filter code will be modfified to output */
	/* MultiInspiralTable instead of TwoInterfInspiralEvent   */
        savedEvents.multiInspiralTable->end_time = eventList->time;
        savedEvents.multiInspiralTable->impulse_time = eventList->impulseTime;
        savedEvents.multiInspiralTable->eff_distance = eventList->effDist;
        savedEvents.multiInspiralTable->eta = (REAL4) eventList->tmplt.eta;
        savedEvents.multiInspiralTable->mass1 = (REAL4) eventList->tmplt.mass1;
        savedEvents.multiInspiralTable->mass2 = (REAL4) eventList->tmplt.mass2;
        savedEvents.multiInspiralTable->tau0 = (REAL4) eventList->tmplt.t0;
        savedEvents.multiInspiralTable->tau2 = (REAL4) eventList->tmplt.t2;
        savedEvents.multiInspiralTable->tau3 = (REAL4) eventList->tmplt.t3;
        savedEvents.multiInspiralTable->tau4 = (REAL4) eventList->tmplt.t4;
        savedEvents.multiInspiralTable->tau5 = (REAL4) eventList->tmplt.t5;
        savedEvents.multiInspiralTable->snr = sqrt(eventList->snrsq);

        /*Temporarily set "network" chisq to the average value*/
        savedEvents.multiInspiralTable->chisq = 
          0.5 * (eventList->chisq1+eventList->chisq2);
        
        savedEvents.multiInspiralTable->sigmasq = 
          eventList->sigma * eventList->sigma;
        savedEvents.multiInspiralTable->ligo_axis_ra = 
          eventList->twoInterfAxisRa;
        savedEvents.multiInspiralTable->ligo_axis_dec = 
          eventList->twoInterfAxisDec;

        /*Wave arrival angle with respect to detector1-to-detector2 
          (e.g., Hanford-to-Livingston) ray/baseline...*/           
        savedEvents.multiInspiralTable->ligo_angle 
          = eventList->twoInterfAngle;
      }
      else
      {
        event->next = eventList;
      }
#endif
      /* save a pointer to the last event in the list and count the events */
      ++numEvents;
      while ( eventList->next )
	{
        eventList = eventList->next;
        ++numEvents;
      }
      event = eventList;
      eventList = NULL;
    } /* end if ( eventList ) */

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
  for ( n = 0; n < numDetectors; ++n )
  {
    LAL_CALL( LALTwoInterfFindChirpChisqVetoFinalize (&status, 
          fcFilterParams[n].chisqParams, fcInitParams[n]->numChisqBins),
        &status);
    LALFree( fcInitParams[n] );

    LAL_CALL( LALFindChirpSPTemplateFinalize(&status, &fcTmpltParams[n]), 
	      &status );
  }
  LAL_CALL( LALDestroyTwoInterfFindChirpInputVector (&status, 
	&twoInterfFilterInputVec), &status );
  LAL_CALL( LALTwoInterfFindChirpFilterFinalize (&status, 
						 &twoInterfFilterParams), &status );
  LAL_CALL( LALTwoInterfFindChirpSPDataFinalize (&status, 
						 &twoInterfDataParamsVec), &status );
  LAL_CALL( LALDestroyTwoInterfFindChirpSegmentVector( &status, 
						       &twoInterfFcSegVec ),
	    &status );
  LALFree( twoInterfInitParams );
  
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
  for (n = 0; n < numDetectors; ++n) 
    {
      LAL_CALL( LALFinalizeDataSegmentVector( &status, &dataVecPtr[n] ),&status );
      LAL_CALL( LALSDestroyVector( &status, &(chan[n].data) ), &status );
      LAL_CALL( LALSDestroyVector( &status, &(spec[n].data) ), &status );
      LAL_CALL( LALCDestroyVector( &status, &(resp[n].data) ), &status );
      LALFree( dataSegVec[n].data );
    }
  LALFree( dataSegVec );
  LALFree( twoInterfDataSegVec);
  twoInterfDataSegVec= NULL;

  
  /*
   *
   * write the results to disk
   *
   */


  if ( vrbflg ) fprintf( stdout, "writing frame data to disk\n" );

  /* write the output frame */
  if ( writeRawData || writeFilterData || writeResponse || writeSpectrum ||
      writeRhosq || writeCoherentRhosq || writeChisq )
  {
    for ( n =0 ; n < numDetectors ; ++n )
    {
      LALSnprintf( fname, sizeof(fname), "%s-COHERE_INSP-%d-%d.gwf",
          ifo[n], gpsStartTime.gpsSeconds,
          gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
      frOutFile = FrFileONew( fname, 0 );
      FrameWrite( outFrame[n], frOutFile );
    }
    FrameWrite( outFrameNet, frOutFile );
    FrFileOEnd( frOutFile );
  }

cleanexit:

  if ( vrbflg ) fprintf( stdout, "writing xml data to disk\n" );

  /* open the output xml file */
  memset( &results, 0, sizeof(LIGOLwXMLStream) );
  LALSnprintf( fname, sizeof(fname), "%s%s-COHERE_INSP-%d-%d.xml",
      ifo[0], ifo[1], gpsStartTime.gpsSeconds,
      gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &results, fname), &status );

  /* write the process table */
  LALSnprintf( proctable.processTable->ifos, LIGOMETA_IFOS_MAX, "%s", ifo );
  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->end_time),
        &accuracy ), &status );
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, process_table ), 
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, proctable, 
        process_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  free( proctable.processTable );

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
    free( this_proc_param );
  }

  /* write the search summary table */
  if ( numTmplts )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, 
          search_summary_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, searchsumm, 
          search_summary_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  }
  free( searchsumm.searchSummaryTable );

  /* write the search summvars table */
  if ( numTmplts )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, 
          search_summvars_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, searchsummvars, 
          search_summvars_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  }
  while( searchsummvars.searchSummvarsTable )
  {
    this_search_summvar = searchsummvars.searchSummvarsTable;
    searchsummvars.searchSummvarsTable = this_search_summvar->next;
    free( this_search_summvar );
  }

  /* write the summvalue table */
  if ( numTmplts )
  {
    LALSnprintf( summvalue.summValueTable->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    summvalue.summValueTable->version = 0;
    summvalue.summValueTable->start_time = gpsStartTime;
    summvalue.summValueTable->end_time = gpsEndTime;
    LALSnprintf( summvalue.summValueTable->ifo, LIGOMETA_IFO_MAX, "%s", ifo );
    LALSnprintf( summvalue.summValueTable->name, LIGOMETA_SUMMVALUE_NAME_MAX, 
        "%s", "cohere_insp_effective_distance" );
    LALSnprintf( summvalue.summValueTable->comment, LIGOMETA_SUMMVALUE_COMM_MAX, 
        "%s", "1.4_1.4_8" );
    summvalue.summValueTable->value = candle.effDistance;
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, summ_value_table ), 
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, summvalue, 
          summ_value_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  }

  /* close the output xml file */
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &results ), &status );

  /* free the rest of the memory, check for memory leaks and exit */
  if ( bankFileName ) free( bankFileName );

  for ( n = 0; n < numDetectors; ++n )
  {
    if ( calCacheName[n] ) free( calCacheName[n] );
    if ( frInCacheName[n] ) free( frInCacheName[n] );
    if ( channelName[n] ) free( channelName[n] );
    if ( fqChanName[n] ) free( fqChanName[n] );
  }

  if ( vrbflg ) fprintf( stdout, "checking memory leaks and exiting\n" );
  LALCheckMemoryLeaks();
  exit( 0 );
}


/* ------------------------------------------------------------------------- */

#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
  LALCalloc( 1, sizeof(ProcessParamsTable) ); \
  LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
      PROGRAM_NAME ); \
      LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
          long_options[option_index].name ); \
  LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
  LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );


#define USAGE \
"lalapps_inspiral [options]\n\n"\
"  --help                       display this message\n"\
"  --verbose                    print progress information\n"\
"  --debug-level LEVEL          set the LAL debug level to LEVEL\n"\
"  --user-tag STRING            set the process_params usertag to STRING\n"\
"  --comment STRING             set the process table comment to STRING\n"\
"\n"\
"  --gps-start-time SEC         GPS second of data start time\n"\
"  --gps-start-time-ns NS       GPS nanosecond of data start time\n"\
"  --gps-end-time SEC           GPS second of data end time\n"\
"  --gps-end-time-ns NS         GPS nanosecond of data end time\n"\
"  --pad-data T                 pad the data start and end time by T seconds\n"\
"\n"\
"  --site0                      first detector's site number [0 for LHO]\n"\
"  --site1                      second detector's site number [1 for LLO]\n"\
"  --ifo1-frame-cache                obtain ifo1 frame data from LAL frame cache FILE\n"\
"  --ifo2-frame-cache                obtain ifo2 frame data from LAL frame cache FILE\n"\
"  --ifo1-calibration-cache FILE     obtain ifo1 calibration from LAL frame cache FILE\n"\
"  --ifo2-calibration-cache FILE     obtain ifo2 calibration from LAL frame cache FILE\n"\
"  --channel-name-ifo1 CHAN          read data from the 1st interferometer's channel CHAN\n"\
"  --channel-name-ifo2 CHAN          read data from the 2nd interferometer's channel CHAN\n"\
"\n"\
"\n"\
"  --bank-file FILE             read template bank parameters from FILE\n"\
"  --minimal-match M            override bank minimal match with M (sets delta)\n"\
"  --start-template N           start filtering at template number N in bank\n"\
"  --stop-templateN             stop filtering at template number N in bank\n"\
"\n"\
"  --sample-rate F              filter data at F Hz, downsampling if necessary\n"\
"  --resample-filter TYPE       set resample filter to TYPE (ldas|butterworth)\n"\
"\n"\
"  --disable-high-pass          turn off the IIR highpass filter\n"\
"  --enable-high-pass-ifo1 F    high pass ifo1 data above F Hz using an IIR filter\n"\
"  --enable-high-pass-ifo2 F    high pass ifo2 data above F Hz using an IIR filter\n"\
"  --spectrum-type TYPE         use PSD estimator TYPE (mean|median)\n"\
"\n"\
"  --segment-length N           set data segment length to N points\n"\
"  --number-of-segments N       set number of data segments to N\n"\
"  --segment-overlap N          overlap data segments by N points\n"\
"\n"\
"  --low-frequency-cutoff-ifo1 F     do not filter ifo1 data below F Hz\n"\
"  --low-frequency-cutoff-ifo2 F     do not filter ifo2 below F Hz\n"\
"  --inverse-spec-length-ifo1 T      set length of inverse spectrum of ifo1 to T seconds\n"\
"  --inverse-spec-length-ifo2 T      set length of inverse spectrum of ifo2 to T seconds\n"\
"  --dynamic-range-exponent X   set dynamic range scaling to 2^X\n"\
"\n"\
"  --chisq-bins P               set number of chisq veto bins to P\n"\
"  --ifo1-snr-threshold RHO          set signal-to-noise threshold in ifo1 to RHO\n"\
"  --ifo2-snr-threshold RHO          set signal-to-noise threshold in ifo2 to RHO\n"\
"  --ifo1-chisq-threshold X          threshold on ifo1 chi^2 < X * ( p + rho^2 * delta^2 )\n"\
"  --ifo2-chisq-threshold X          threshold on ifo2 chi^2 < X * ( p + rho^2 * delta^2 )\n"\
"  --enable-event-cluster       turn on maximization over chirp length\n"\
"  --disable-event-cluster      turn off maximization over chirp length\n"\
"\n"\
"  --enable-output              write the results to a LIGO LW XML file\n"\
"  --disable-output             do not write LIGO LW XML output file\n"\
"\n"\
"  --write-raw-data             write raw data to a frame file\n"\
"  --write-filter-data          write data that is passed to filter to a frame\n"\
"  --write-response             write the computed response function to a frame\n"\
"  --write-spectrum             write the uncalibrated psd to a frame\n"\
"  --write-snrsq                write the snr time series for each data segment\n"\
"  --write-chisq                write the r^2 time series for each data segment\n"\
"\n"

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
    {"channel-name-ifo1",       required_argument, 0,                'c'},
    {"channel-name-ifo2",       required_argument, 0,                'C'},
    {"segment-length",          required_argument, 0,                'd'},
    {"number-of-segments",      required_argument, 0,                'e'},
    {"segment-overlap",         required_argument, 0,                'f'},
    {"sample-rate",             required_argument, 0,                'g'},
    {"help",                    no_argument,       0,                'h'},
    {"low-freq-cutoff-ifo1",    required_argument, 0,                'i'},
    {"low-freq-cutoff-ifo2",    required_argument, 0,                'I'},
    {"spectrum-type",           required_argument, 0,                'j'},
    {"inverse-spec-length-ifo1",required_argument, 0,                'k'},
    {"inverse-spec-length-ifo2",required_argument, 0,                'K'},
    {"dynamic-range-exponent",  required_argument, 0,                'l'},
    {"start-template",          required_argument, 0,                'm'},
    {"minimal-match",           required_argument, 0,                'M'},
    {"stop-template",           required_argument, 0,                'n'},
    {"chisq-bins",              required_argument, 0,                'o'},
    {"ifo1-calibration-cache",  required_argument, 0,                'p'},
    {"ifo2-calibration-cache",  required_argument, 0,                'P'},
    {"ifo1-snr-threshold",      required_argument, 0,                'q'},
    {"ifo2-snr-threshold",      required_argument, 0,                'Q'},
    {"ifo1-chisq-threshold",    required_argument, 0,                'r'},
    {"ifo2-chisq-threshold",    required_argument, 0,                'S'},
    {"resample-filter",         required_argument, 0,                'R'},
    {"comment",                 required_argument, 0,                's'},
    {"enable-high-pass-ifo1",   required_argument, 0,                't'},
    {"enable-high-pass-ifo2",   required_argument, 0,                'T'},
    {"ifo1-frame-cache",        required_argument, 0,                'u'},
    {"ifo2-frame-cache",        required_argument, 0,                'U'},
    {"bank-file",               required_argument, 0,                'v'},
    {"injection-file",          required_argument, 0,                'w'},
    {"pad-data",                required_argument, 0,                'x'},
    {"coherent-snr-threshold",  required_argument, 0,                'X'},
    {"site0",                   required_argument, 0,                'y'},
    {"site1",                   required_argument, 0,                'Y'},
    {"debug-level",             required_argument, 0,                'z'},
    {"user-tag",                required_argument, 0,                'Z'},
    {"userTag",                 required_argument, 0,                'Z'},
    /* frame writing options */
    {"write-raw-data",          no_argument,       &writeRawData,     1 },
    {"write-filter-data",       no_argument,       &writeFilterData,  1 },
    {"write-response",          no_argument,       &writeResponse,    1 },
    {"write-spectrum",          no_argument,       &writeSpectrum,    1 },
    {"write-detector-snrsq",    no_argument,       &writeRhosq,       1 },
    {"write-coherent-snrsq",    no_argument,       &writeCoherentRhosq, 1 },
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
    size_t optarg_len;

    c = getopt_long_only( argc, argv, 
        "a:A:b:B:c:C:d:e:f:g:h:i:j:k:K:l:m:M:n:o:p:P:q:Q:r:S:R:s:t:T:u:U:v:w:x:X:y:Y:z:Z:", 
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
          optarg_len = strlen( optarg ) + 1;
          fqChanName[0] = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
          memcpy( fqChanName[0], optarg, optarg_len );
          ADD_PROCESS_PARAM( "string", "%s", optarg );

          /* check that we have a proper channel name */
          if ( ! (channamptr = strstr( fqChanName[0], ":" ) ) )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "channel name must be a full LIGO channel name "
                "e.g. L1:LSC-AS_Q\n(%s specified)\n",
                long_options[option_index].name, optarg );
            exit( 1 );
          }
          optarg_len = strlen( ++channamptr ) + 1;
          channelName[0] = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
          memcpy( channelName[0], channamptr, optarg_len );

          /* copy the first two characters to the ifo name */
          memset( ifo[0], 0, sizeof(ifo[0]) );
          memcpy( ifo[0], optarg, sizeof(ifo[0]) - 1 );
        }
        break;

      case 'C':
        {
          /* create storage for the channel name and copy it */
          char *channamptr = NULL;
          size_t optarg_len = strlen( optarg ) + 1;
          fqChanName[1] = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
          memcpy( fqChanName[1], optarg, optarg_len );
          ADD_PROCESS_PARAM( "string", "%s", optarg );

          /* check that we have a proper channel name */
          if ( ! (channamptr = strstr( fqChanName[1], ":" ) ) )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "channel name must be a full LIGO channel name "
                "e.g. L1:LSC-AS_Q\n(%s specified)\n",
                long_options[option_index].name, optarg );
            exit( 1 );
          }
          optarg_len = strlen( ++channamptr ) + 1;
          channelName[1] = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
          memcpy( channelName[1], channamptr, optarg_len );

          /* copy the first two characters to the ifo name */
          memset( ifo[1], 0, sizeof(ifo[1]) );
          memcpy( ifo[1], optarg, sizeof(ifo[1]) - 1 );
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
	fLow[0] = (REAL4) atof( optarg );
	if ( fLow[0] < 40 )
	{
          fprintf( stdout, "invalid argument to --%s:\n"
              "low frequency cutoff in ifo1 is less than 40 Hz: "
              "(%f Hz specified)\n",
              long_options[option_index].name, fLow[0] );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", fLow[0] );
        break;

      case 'I':
	fLow[1] = (REAL4) atof( optarg );
	if ( fLow[1] < 40 )
	{
          fprintf( stdout, "invalid argument to --%s:\n"
              "low frequency cutoff in ifo2 is less than 40 Hz: "
              "(%f Hz specified)\n",
              long_options[option_index].name, fLow[1] );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", fLow[1] );
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
        invSpecTrunc[0] = (INT4) atoi( optarg );
        if ( invSpecTrunc[0] < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "inverse spectrum length in ifo1 must be positive or zero: "
              "(%d specified)\n", 
              long_options[option_index].name, invSpecTrunc[0] );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", invSpecTrunc[0] );
        break;

      case 'K':
        invSpecTrunc[1] = (INT4) atoi( optarg );
        if ( invSpecTrunc[1] < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "inverse spectrum length in ifo1 must be positive or zero: "
              "(%d specified)\n", 
              long_options[option_index].name, invSpecTrunc[1] );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", invSpecTrunc[1] );
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
        minimalMatch = (REAL4) atof( optarg );
        if ( minimalMatch < 0 || minimalMatch > 1 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "minimal match must be in the range [0,1]: "          
              "(%e specified)\n", 
              long_options[option_index].name, minimalMatch );
        }
        /* process param added after bank is generated so that a */
        /* value in the bank looks like a command line option.   */
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
          /* create storage for detector 1's calibration frame cache name */
          optarg_len = strlen(optarg) + 1;
          calCacheName[0] = (CHAR *) calloc( optarg_len, sizeof(CHAR));
          memcpy( calCacheName[0], optarg, optarg_len );
          ADD_PROCESS_PARAM( "string", "%s", optarg );
        }
        break;

      case 'P':
        {
          /* create storage for detector 2's calibration frame cache name */
          optarg_len = strlen(optarg) + 1;
          calCacheName[1] = (CHAR *) calloc( optarg_len, sizeof(CHAR));
          memcpy( calCacheName[1], optarg, optarg_len );
          ADD_PROCESS_PARAM( "string", "%s", optarg );
        }
        break;

      case 'q':
        snrThresh[0] = atof( optarg );
        if ( snrThresh[0] < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "signal to noise threshold in ifo1 must be positive: "
              "(%f specified)\n", 
              long_options[option_index].name, snrThresh[0] );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'Q':
        snrThresh[1] = atof( optarg );
        if ( snrThresh[1] < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "signal to noise threshold in ifo2 must be positive: "
              "(%f specified)\n", 
              long_options[option_index].name, snrThresh[1] );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'r':
        chisqThresh[0] = atof( optarg );
        if ( chisqThresh[0] < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "chi squared threshold in ifo1 must be positive: "
              "(%f specified)\n", 
              long_options[option_index].name, chisqThresh[0] );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'S':
        chisqThresh[1] = atof( optarg );
        if ( chisqThresh[1] < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "chi squared threshold in ifo2 must be positive: "
              "(%f specified)\n", 
              long_options[option_index].name, chisqThresh[1] );
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
          LALSnprintf( comment, LIGOMETA_COMMENT_MAX, "%s", optarg);
        }
        break;

      case 't':
        highPass = 1;
        highPassFreq[0] = (REAL4) atof( optarg );
        if ( highPassFreq[0] < 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "low frequency cutoff in ifo1 is less than 0 Hz: "
              "(%f Hz specified)\n",
              long_options[option_index].name, highPassFreq[0] );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", highPassFreq[0] );
        break;

      case 'T':
        highPass = 1;
        highPassFreq[1] = (REAL4) atof( optarg );
        if ( highPassFreq[1] < 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "low frequency cutoff in ifo2 is less than 0 Hz: "
              "(%f Hz specified)\n",
              long_options[option_index].name, highPassFreq[1] );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", highPassFreq[1] );
        break;

      case 'u':
        {
          /* create storage for the input frame cache name */
          optarg_len = strlen(optarg) + 1;
          frInCacheName[0] = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
          memcpy( frInCacheName[0], optarg, optarg_len );
          ADD_PROCESS_PARAM( "string", "%s", optarg );
        }
        break;

      case 'U':
        {
          /* create storage for the input frame cache name */
          optarg_len = strlen(optarg) + 1;
          frInCacheName[1] = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
          memcpy( frInCacheName[1], optarg, optarg_len );
          ADD_PROCESS_PARAM( "string", "%s", optarg );
        }
        break;

      case 'v':
        {
          /* create storage for the calibration frame cache name */
          optarg_len = strlen(optarg) + 1;
          bankFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
          memcpy( bankFileName, optarg, optarg_len );
          ADD_PROCESS_PARAM( "string", "%s", optarg );
        }
        break;

      case 'w':
        {
          /* create storage for the injection file name */
          optarg_len = strlen(optarg) + 1;
          injectionFile = (CHAR *) calloc( optarg_len, sizeof(CHAR));
          memcpy( injectionFile, optarg, optarg_len );
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

      case 'X':
        coherentSnrThresh = atof( optarg );
        if ( coherentSnrThresh < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "network signal to noise threshold must be positive: "
              "(%f specified)\n", 
              long_options[option_index].name, coherentSnrThresh  );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'y':
        site0 = (UINT4) atoi( optarg );
        if ( site0 < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "the lalCachedDetector site identifier must be positive: "
              "(%d specified)\n", 
              long_options[option_index].name, site0  );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", optarg );
        break;

      case 'Y':
        site1 = (UINT4) atoi( optarg );
        if ( site1 < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "the lalCachedDetector site identifier must be positive: "
              "(%d specified)\n", 
              long_options[option_index].name, site1  );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", optarg );
        break;

      case 'z':
        set_debug_level( optarg );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'Z':
	/* create storage for the usertag */
        optarg_len = strlen( optarg ) + 1;
        userTag = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
        memcpy( userTag, optarg, optarg_len );

        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          calloc( 1, sizeof(ProcessParamsTable) );
        LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", 
            PROGRAM_NAME );
        LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "-userTag" );
        LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%s",
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
    LALSnprintf( procparams.processParamsTable->program, 
        LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME );
    LALSnprintf( procparams.processParamsTable->param,
        LIGOMETA_PARAM_MAX, "--enable-output" );
    LALSnprintf( procparams.processParamsTable->type, 
        LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( procparams.processParamsTable->value, 
        LIGOMETA_TYPE_MAX, " " );
  }
  else if ( enableOutput == 0 )
  {
    LALSnprintf( procparams.processParamsTable->program, 
        LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME );
    LALSnprintf( procparams.processParamsTable->param,
        LIGOMETA_PARAM_MAX, "--disable-output" );
    LALSnprintf( procparams.processParamsTable->type, 
        LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( procparams.processParamsTable->value, 
        LIGOMETA_TYPE_MAX, " " );
  }
  else
  {
    fprintf( stderr, "--enable-output or --disable-output "
        "argument must be specified\n" );
    exit( 1 );
  }


  /* check event cluster option */
  this_proc_param = this_proc_param->next = (ProcessParamsTable *)
    calloc( 1, sizeof(ProcessParamsTable) );
  if ( eventCluster == 1 )
  {
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--enable-event-cluster" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }
  else if ( eventCluster == 0 )
  {
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--disable-event-cluster" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
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
  else if ( ! highPass )
  {
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--disable-high-pass" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
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
  if ( fLow[0] < 0 )
  {
    fprintf( stderr, "--low-frequency-cutoff in ifo1 must be specified\n" );
    exit( 1 );
  }
  if ( fLow[1] < 0 )
  {
    fprintf( stderr, "--low-frequency-cutoff in ifo2 must be specified\n" );
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
  if ( invSpecTrunc[0] < 0 )
  {
    fprintf( stderr, "--inverse-spec-length in ifo1 must be specified\n" );
    exit( 1 );
  }
  else if ( invSpecTrunc[0] * sampleRate > numPoints )
  {
    fprintf( stderr, "--inverse-spec-length in ifo1 must be less than "
        "--segment-length\n" );
    exit( 1 );
  }
  if ( invSpecTrunc[1] < 0 )
  {
    fprintf( stderr, "--inverse-spec-length in ifo2 must be specified\n" );
    exit( 1 );
  }
  else if ( invSpecTrunc[1] * sampleRate > numPoints )
  {
    fprintf( stderr, "--inverse-spec-length in ifo2 must be less than "
        "--segment-length\n" );
    exit( 1 );
  }

  if ( ! haveDynRange )
  {
    fprintf( stderr, "--dynamic-range-exponent must be specified\n" );
    exit( 1 );
  }

  /* check that a channel has been requested and fill the ifo */
  if ( ! (fqChanName[0] && fqChanName[1]))
  {
    fprintf( stderr, "--channel-names of detectors must be specified\n" );
    exit( 1 );
  }

  /* check that the thresholds have been specified */
  if ( (snrThresh[0] < 0 &&  snrThresh[1] < 0 && coherentSnrThresh < 0) )
  {
    fprintf( stderr, "--snr-thresholds must be specified\n" );
    exit( 1 );
  }
  if ( chisqThresh[0] < 0 || chisqThresh[1] < 0 )
  {
    fprintf( stderr, "--chisq-thresholds must be specified\n" );
    exit( 1 );
  }

  /* check that the frame caches have been specified */
  if ( ! (frInCacheName[0] && frInCacheName[1]) )
  {
    fprintf( stderr, "--frame-caches must be specified\n" );
    exit( 1 );
  }
  if ( ! (calCacheName[0] && calCacheName[1]) )
  {
    fprintf( stderr, "--calibration-caches must be specified\n" );
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
