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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <regex.h>
#include <time.h>

#include <lalapps.h>
#include <processtable.h>
#include <lal/LALConfig.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/AVFactories.h>
#include <lal/LALDatatypes.h>
#include <lal/LALConstants.h>
#include <lal/FrameStream.h>
#include <lal/IIRFilter.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOLwXML.h>
#include <lal/Date.h>

#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/FindChirpChisq.h>

#include "inspiral.h"

#define RESULT_FILE "results.xml"

RCSID( "$Id$" );

#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"

#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
  this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
    LALCalloc( 1, sizeof(ProcessParamsTable) ); \
  LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
      PROGRAM_NAME ); \
  LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
      long_options[option_index].name ); \
  LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
  LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );


long long atoll(const char *nptr);


/*
 *
 * variables that control program behaviour
 *
 */

/* debugging */
extern int vrbflg;                      /* verbocity of lal function    */
int   debugflg = 0;                     /* internal debugging flag      */

/* input data parameters */
UINT8  gpsStartTimeNS   = 0;            /* input data GPS start time    */
UINT8  gpsEndTimeNS     = 0;            /* input data GPS end time      */
CHAR  *channelName      = NULL;         /* name of data channel         */
INT4  numPoints         = -1;           /* points in a segment          */
INT4  numSegments       = -1;           /* number of segments           */
INT4  ovrlap            = -1;           /* overlap between segments     */

/* data conditioning parameters */
INT4   sampleRate       = -1;           /* sample rate of filter data   */
REAL4  fLow             = -1;           /* low frequency cutoff         */
INT4   specType         = -1;           /* given, median or mean psd    */
INT4   invSpecTrunc     = -1;           /* length of inverse spec (s)   */
REAL4  dynRangeExponent = -1;           /* exponent of dynamic range    */

/* matched filter parameters */
INT4  startTemplate     = -1;           /* index of first template      */
INT4  stopTemplate      = -1;           /* index of last template       */
INT4  numChisqBins      = -1;           /* number of chisq bins         */
INT4  hierDepth         = -1;           /* number of levels in search   */
char  *rhosqStr         = NULL;         /* string of rhosq thresholds   */
char  *chisqStr         = NULL;         /* string of chisq thresholds   */
REAL4 *rhosqThresh      = NULL;         /* signal to noise thresholds   */
REAL4 *chisqThresh      = NULL;         /* chisq veto thresholds        */
int    eventCluster     = -1;           /* perform chirplen clustering  */

/* output parameters */
int    enableOutput     = -1;           /* write out inspiral events    */
int    writeRawData     = 0;            /* write the raw data to a file */
int    writeFilterData  = 0;            /* write post injection data    */
int    writeResponse    = 0;            /* write response function used */
int    writeSpec        = 0;            /* write computed psd to file   */
int    writeRhosq       = 0;            /* write rhosq time series      */
int    writeChisq       = 0;            /* write chisq time series      */


int main( int argc, char *argv[] )
{



  /*
   *
   * general variables in execution
   *
   */


  /* current getopt argument */
  int                   c;

  /* lal function variables */
  LALStatus             status = blank_status;
  LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;
  INT4                  inputDataLength = 0;

  /* input data */
  FrStream                     *frStream       = NULL;
  FrChanIn                     *frChan         = NULL;
  REAL4TimeSeries              *dataChannel    = NULL;
  
  /* findchirp data structures */
  DataSegmentVector            *dataSegVec     = NULL;
  FindChirpInitParams          *fcInitParams   = NULL;
  FindChirpSegmentVector       *fcSegVec       = NULL;
  FindChirpSPDataParams        *fcDataParams   = NULL;
  FindChirpSPTmpltParams       *fcTmpltParams  = NULL;
  FindChirpFilterParams        *fcFilterParams = NULL;
  FindChirpFilterInput         *fcFilterInput  = NULL;
  
  /* output data */
  MetadataTable         proctable;
  MetadataTable         procparams;
  ProcessParamsTable   *this_proc_param;
  LIGOLwXMLStream       results;


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
  
  

  /*
   *
   * parse command line arguments
   *
   */


  while ( 1 )
  {
    static struct option long_options[] =
    {
      /* these options set a flag */
      {"verbose",                 no_argument,       &vrbflg,           1 },
      {"debug",                   no_argument,       &debugflg,         1 },
      {"enable-event-cluster",    no_argument,       &eventCluster,     1 },
      {"disable-event-cluster",   no_argument,       &eventCluster,     0 },
      {"enable-output",           no_argument,       &enableOutput,     1 },
      {"disable-output",          no_argument,       &enableOutput,     0 },
      {"write-raw-data",          no_argument,       &writeRawData,     1 },
      {"write-filter-data",       no_argument,       &writeFilterData,  1 },
      {"write-response",          no_argument,       &writeResponse,    1 },
      {"write-spec",              no_argument,       &writeSpec,        1 },
      {"write-rhosq",             no_argument,       &writeRhosq,       1 },
      {"write-chisq",             no_argument,       &writeChisq,       1 },
      /* these options don't set a flag */
      {"gps-start-time",          required_argument, 0,                'a'},
      {"gps-end-time",            required_argument, 0,                'b'},
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
      {"stop-template",           required_argument, 0,                'n'},
      {"chisq-bins",              required_argument, 0,                'o'},
      {"hierarchy-depth",         required_argument, 0,                'p'},
      {"rhosq-thresholds",        required_argument, 0,                'q'},
      {"chisq-thresholds",        required_argument, 0,                'r'},
      {"comment",                 required_argument, 0,                's'},
      {0, 0, 0, 0}
    };

    /* getopt_long stores long option here */
    int option_index = 0;

    c = getopt_long( argc, argv, 
        "a:b:c:d:e:f:g:h:i:j:k:l:m:n:o:p:q:r:s:", long_options, &option_index );

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
          long long int gstartt = atoll( optarg );
          if ( gstartt < 441417609 )
          {
            fprintf( stderr, "GPS start time is prior to " 
                "Jan 01, 1994  00:00:00 UTC: %lld\n", gstartt );
            exit( 1 );
          }
          gpsStartTimeNS = (UINT8) gstartt * 1000000000LL;
          ADD_PROCESS_PARAM( "int", "%lld", gstartt );
        }
        break;

      case 'b':
        {
          long long int gendt = atoll( optarg );
          if ( gendt > 999999999 )
          {
            fprintf( stderr, "GPS end time is after " 
                "Sep 14, 2011  01:46:26 UTC: %lld\n", gendt );
            exit( 1 );
          }
          gpsEndTimeNS = (UINT8) gendt * 1000000000LL;
          ADD_PROCESS_PARAM( "int", "%lld", gendt );
        }
        break;

      case 'c':
        {
          size_t chanlen = strlen( optarg );
          channelName = (CHAR *) LALMalloc( ++chanlen );
          memcpy( channelName, optarg, chanlen );
          ADD_PROCESS_PARAM( "string", "%s", optarg );
        }
        break;

      case 'd':
        numPoints = (INT4) atoi( optarg );
        if ( numPoints < 2 || numPoints % 2 )
        {
          fprintf( stderr, "invalid number of points in data segment: %d\n",
              numPoints );
          exit( 1 );
        }
        break;

      case 'e':
        numSegments = (INT4) atoi( optarg );
        if ( numSegments < 1 )
        {
          fprintf( stderr, "invalid number data segments: %d\n", numSegments );
          exit( 1 );
        }
        break;

      case 'f':
        ovrlap = (INT4) atoi( optarg );
        if ( ovrlap < 0 )
        {
          fprintf( stderr, "invalid data segment overlap: %d\n", ovrlap );
          exit( 1 );
        }
        break;

      case 'g':
        sampleRate = (INT4) atoi( optarg );
        if ( sampleRate < 2 || sampleRate > 16384 || sampleRate % 2 )
        {
          fprintf( stderr, "invalid sample rate: %d\n", sampleRate );
          exit( 1 );
        }
        break;

      case 'h':
        fprintf( stdout, USAGE );
        exit( 0 );
        break;

      case 'i':
        fLow = (REAL4) atof( optarg );
        if ( fLow < 40 )
        {
          fprintf( stdout, "low frequency cutoff is less than 40 Hz: %f Hz\n",
              fLow );
          exit( 1 );
        }
        break;

      case 'j':
        if ( ! strcmp( "file", optarg ) )
        {
          specType = 0;
        }
        else if ( ! strcmp( "mean", optarg ) )
        {
          specType = 1;
        }
        else if ( ! strcmp( "median", optarg ) )
        {
          specType = 2;
        }
        else
        {
          fprintf( stderr, "unknown power spectrum type: %s\n", optarg );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'k':
        invSpecTrunc = (INT4) atoi( optarg );
        if ( invSpecTrunc < 1 )
        {
          fprintf( stderr, "length of inverse spectrum is less than 1 second"
              ": %d\n", invSpecTrunc );
          exit( 1 );
        }
        break;

      case 'l':
        dynRangeExponent = (REAL4) atof( optarg );
        break;
          
      case 'm':
        startTemplate = (INT4) atoi( optarg );
        if ( startTemplate < 0 )
        {
          fprintf( stderr, "template bank start index is invalid: %d\n", 
              startTemplate );
          exit( 1 );
        }
        break;

      case 'n':
        stopTemplate = (INT4) atoi( optarg );
        if ( stopTemplate < 0 )
        {
          fprintf( stderr, "template bank stop index is invalid: %d\n", 
              stopTemplate );
          exit( 1 );
        }
        break;

      case 'o':
        numChisqBins = (INT4) atoi( optarg );
        if ( numChisqBins < 0 )
        {
          fprintf( stderr, "invalid number of chisq veto bins: %d\n", 
              numChisqBins );
          exit( 1 );
        }
        break;

      case 'p':
        hierDepth = (INT4) atoi( optarg );
        if ( hierDepth < 0 || hierDepth > 1 )
        {
          fprintf( stderr, "invalid hierarchical search depth: %d\n", 
              hierDepth );
          exit( 1 );
        }
        break;

      case 'q':
        {
          size_t rhosqlen = strlen( optarg );
          rhosqStr = (char *) LALMalloc( ++rhosqlen );
          memcpy( rhosqStr, optarg, rhosqlen );
          ADD_PROCESS_PARAM( "string", "%s", optarg );
        }
        break;

      case 'r':
        {
          size_t chisqlen = strlen( optarg );
          chisqStr = (char *) LALMalloc( ++chisqlen );
          memcpy( chisqStr, optarg, chisqlen );
          ADD_PROCESS_PARAM( "string", "%s", optarg );
        }
        break;

      case 's':
        if ( strlen( optarg ) > LIGOMETA_COMMENT_MAX - 1 )
        {
          fprintf( stderr, "--comment must be less than %d characters",
              LIGOMETA_COMMENT_MAX );
          exit( 1 );
        }
        else
        {
          LALSnprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX,
              "%s", optarg);
        }
        break;

      case '?':
        fprintf( stderr, USAGE );
        exit( 1 );
        break;

      default:
        fprintf( stderr, "unknown error while parsing options\n\n" );
        fprintf( stderr, USAGE );
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
        LIGOMETA_TYPE_MAX, "int" );
    LALSnprintf( procparams.processParamsTable->value, 
        LIGOMETA_TYPE_MAX, "1" );
  }
  else if ( enableOutput == 0 )
  {
    LALSnprintf( procparams.processParamsTable->program, 
        LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME );
    LALSnprintf( procparams.processParamsTable->param,
        LIGOMETA_PARAM_MAX, "--disable-output" );
    LALSnprintf( procparams.processParamsTable->type, 
        LIGOMETA_TYPE_MAX, "int" );
    LALSnprintf( procparams.processParamsTable->value, 
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
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--enable-event-cluster" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "int" );
    LALSnprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, "1" );
  }
  else if ( eventCluster == 0 )
  {
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--disable-event-cluster" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "int" );
    LALSnprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, "1" );
  }
  else
  {
    fprintf( stderr, "--enable-event-cluster or --disable-event-cluster "
        "argument must be specified\n" );
    exit( 1 );
  }


  /*
   *
   * check validity of arguments
   *
   */


  /* check validity of input data time */
  if ( ! gpsStartTimeNS || ! gpsEndTimeNS )
  {
    fprintf( stderr, "gps start and end time must be specified\n" );
    exit( 1 );
  }

  if ( gpsEndTimeNS <= gpsStartTimeNS )
  {
    fprintf( stderr, "invalid gps time range: "
        "start time: %lld, end time %lld\n",
        gpsStartTimeNS / 1000000000LL, gpsEndTimeNS / 1000000000LL );
    exit( 1 );
  }

  /* check validity of data length parameters */
  if ( numPoints < 0 )
  {
    fprintf( stderr, "segment length must be specified\n" );
    exit( 1 );
  }
  if ( numSegments < 0 )
  {
    fprintf( stderr, "number data segments must be specified\n" );
    exit( 1 );
  }
  if ( ovrlap < 0 )
  {
    fprintf( stderr, "data segment overlap must be specified\n" );
    exit( 1 );
  }

  /* check sample rate has been sepcified */
  if ( sampleRate < 0 )
  {
    fprintf( stderr, "filter sample rate must be specified\n" );
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
    fprintf( stderr, "number of chi squared veto bins must be positive\n" );
    exit( 1 );
  }
  if ( fLow < 0 )
  {
    fprintf( stderr, "low frequency cutoff must be specified\n" );
    exit( 1 );
  }
  if ( specType < 0 )
  {
    fprintf( stderr, "power spectrum estimation type must be specified\n" );
    exit( 1 );
  }
  if ( invSpecTrunc < 0 )
  {
    fprintf( stderr, "duration of inverse power spectrum must be specified\n" );
    exit( 1 );
  }
  if ( dynRangeExponent < 0 )
  {
    fprintf( stderr, "dynamic range exponent must be specified\n" );
    exit( 1 );
  }


  /*
   *
   * create the data structures needed for the filtering
   *
   */


  /* create and populate findchip initialization structure */
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

  /* create data storage */
  if ( ! ( dataChannel = (REAL4TimeSeries *) 
      LALCalloc( 1, sizeof(REAL4TimeSeries) ) ) )
  {
    fprintf( stderr, "could not allocate memory for input data time series\n" );
    exit( 1 );
  }
  LAL_CALL( LALCreateFindChirpSegmentVector( &status, &fcSegVec, 
        fcInitParams ), &status );

  /* initialize findchirp stationary phase routines */
  LAL_CALL( LALFindChirpSPDataInit( &status, &fcDataParams, fcInitParams ), 
      &status );
  LAL_CALL( LALFindChirpSPTemplateInit( &status, &fcTmpltParams, 
        fcInitParams ), &status );

  fcDataParams->invSpecTrunc = invSpecTrunc * sampleRate;
  fcDataParams->fLow = fLow;
  fcDataParams->dynRange = fcTmpltParams->dynRange = 
    pow( 2.0, dynRangeExponent );
  fcTmpltParams->fLow = fLow;

  /* initialize findchirp filter functions */
  LAL_CALL( LALFindChirpFilterInit( &status, &fcFilterParams, fcInitParams ), 
      &status );
  fcFilterParams->computeNegFreq = 0;

  LAL_CALL( LALCreateFindChirpInput( &status, &fcFilterInput, fcInitParams ), 
      &status );
  LAL_CALL( LALFindChirpChisqVetoInit( &status, fcFilterParams->chisqParams, 
        fcInitParams->numChisqBins, fcInitParams->numPoints ), 
      &status );



  /*
   *
   * read in the input data
   *
   */
  

  /* read the data channel time series from frames */

  /* call the magic calibration function to get the calibration */
  
  /* read in the template bank from a text file */


  /*
   *
   * pre-condition the data channel
   *
   */
  

  /* band pass the data to get rid of low frequency crap */

  /* compute the median power spectrum for the data channel */


  /*
   *
   * findchirp engine
   *
   */
  

  /* sleepy bye-byes */
  sleep( 2 );


  /*
   *
   * free the structures used for the data filtering
   *
   */

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
  LALFree( dataChannel );
  

  /*
   *
   * write the result results to disk
   *
   */


  /* store the job end time in the process table */

  /* open the output xml file */
  memset( &results, 0, sizeof(LIGOLwXMLStream) );
  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &results, RESULT_FILE ), &status );

  /* write the process table */
  sprintf( proctable.processTable->ifos, "L1" );
  sprintf( proctable.processTable->comment, "monkey" );
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

  /* write the inspiral events to the file */

  /* close the output xml file */
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &results ), &status );


  /* check for memory leaks and exit */
  LALCheckMemoryLeaks();
  exit( 0 );
}
