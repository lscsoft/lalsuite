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

#include <lalapps.h>
#include <lal/LALConfig.h>
#include <lal/LALStdlib.h>
#include <lal/Date.h>

#include "inspiral.h"
#include "inspiralxml.h"

RCSID( "$Id$" );

#define PROGRAM_NAME "lalapps_inspiral"
#define CVS_VERSION  "$Revision$"
#define CVS_SOURCE   "$Source$"

long long atoll(const char *nptr);


/*
 *
 * command line arguments
 *
 */


/* global debugging options */
extern int vrbflg;                      /* verbocity of lal function    */
int   debugflg = 0;                     /* internal debugging flag      */

/* input data parameters */
UINT8  gpsStartTimeNS   = 0;            /* input data GPS start time    */
UINT8  gpsStopTimeNS    = 0;            /* input data GPS stop time     */
CHAR  *channelName      = NULL;         /* name of data channel         */
INT4  numPoints         = -1;           /* points in a segment          */
INT4  numSegments       = -1;           /* number of segments           */
INT4  ovrlap            = -1;           /* overlap between segments     */

/* data conditioning parameters */
INT4   sampleRate       = -1;           /* sample rate of filter data   */
int    bandPass         = -1;           /* band pass data before psd    */
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

/* monte carlo options */
char  *monteModule      = NULL;         /* optional monte carlo module  */

/* output parameters */
int    enableOutput     = -1;           /* write out inspiral events    */
int    writeRawData     = 0;            /* write the raw data to a file */
int    writeFilterData  = 0;            /* write post injection data    */
int    writeResponse    = 0;            /* write response function used */
int    writeSpec        = 0;            /* write computed psd to file   */
int    writeRhosq       = 0;            /* write rhosq time series      */
int    writeChisq       = 0;            /* write chisq time series      */

/* other parameters */
char   process_comment[PROCESS_COMMENT_LEN];

/*
 *
 * main program
 *
 */


int main( int argc, char *argv[] )
{

  /*
   *
   * variable declarations for main program
   *
   */


  /* miscellaneous declarations */
  LALStatus             status = blank_status;
  int                   c;

  /* process table and associcated variables */
  time_t        ticks;
  LALDate       laldate;
  LALLeapSecAccuracy accuracy = LALLEAPSEC_LOOSE;
  CHAR          process_program[] = "lalapps_inspiral";
  /* CHAR         *process_version_str = NULL; */
  /* CHAR         *process_cvs_repository_str = NULL; */
  /* LIGOTimeGPS   process_cvs_entry_time; */
  INT4          process_isonline = 0;
  CHAR          process_node[] = "medusa";
  CHAR          process_username[] = "duncan";
  INT4          process_unix_procid = 9999;
  LIGOTimeGPS   process_start_time;
  LIGOTimeGPS   process_end_time;
  INT4          process_jobid = 0;
  CHAR          process_domain[] = "grid";
  
  /* output data files */
  char          resultFile[] = "results.xml";
  char          outputFile[] = "output.xml";
  FILE         *resultFp     = NULL;
  FILE         *outputFp     = NULL;

  /* input data parameters */
  UINT4         inputDataLength = 0;


  /*
   *
   * initialization
   *
   */


  /* set up inital debugging values */
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "1" );

  /* get the start gps time */
  ticks = time( NULL );
  gmtime_r( &ticks, &(laldate.unixDate) );
  laldate.residualNanoSeconds = 0;
  LAL_CALL(
      LALUTCtoGPS( &status, &process_start_time, &laldate, &accuracy ),
      &status );
  
  /* open the output file */
  if ( ! (outputFp = fopen( outputFile, "w" )) )
  {
    perror( "could not open output file for writing" );
    exit( 1 );
  }

  /* set the process comment to null */
  memset( process_comment, 0, PROCESS_COMMENT_LEN * sizeof(CHAR) );


  /*
   *
   * parse command line options and write process_params
   *
   */


  fprintf( outputFp, PROCESS_PARAMS_HEADER );

  while ( 1 )
  {
    static struct option long_options[] =
    {
      /* these options set a flag */
      {"verbose",                 no_argument,       &vrbflg,           1 },
      {"debug",                   no_argument,       &debugflg,         1 },
      {"enable-band-pass",        no_argument,       &bandPass,         1 },
      {"disable-band-pass",       no_argument,       &bandPass,         0 },
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
      {"gps-stop-time",           required_argument, 0,                'b'},
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

    c = getopt_long( argc, argv, "a:b:c:d:e:f:g:h:i:j:k:l:m:n:o:p:q:r:s:",
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
          long long int gstartt = atoll( optarg );
          if ( gstartt < 441417609 )
          {
            fprintf( stderr, "GPS start time is prior to " 
                "Jan 01, 1994  00:00:00 UTC: %lld\n", gstartt );
            exit( 1 );
          }
          fprintf( outputFp, PPARAMS "\"--%s\",\"int\",\"%lld\",\n", 
              long_options[option_index].name, gstartt );
          gpsStartTimeNS = (UINT8) gstartt * 1000000000LL;
        }
        break;

      case 'b':
        {
          long long int gstopt = atoll( optarg );
          if ( gstopt > 999999999 )
          {
            fprintf( stderr, "GPS stop time is after " 
                "Sep 14, 2011  01:46:26 UTC: %lld\n", gstopt );
            exit( 1 );
          }
          fprintf( outputFp, PPARAMS "\"--%s\",\"int\",\"%lld\",\n", 
              long_options[option_index].name, gstopt );
          gpsStopTimeNS = (UINT8) gstopt * 1000000000LL;
        }
        break;

      case 'c':
        {
          size_t chanlen = strlen( optarg );
          channelName = (CHAR *) LALMalloc( ++chanlen );
          memcpy( channelName, optarg, chanlen );
        }
        fprintf( outputFp, PPARAMS "\"--%s\",\"string\",\"%s\",\n", 
            long_options[option_index].name, channelName );
        break;

      case 'd':
        numPoints = (INT4) atoi( optarg );
        if ( numPoints < 2 || numPoints % 2 )
        {
          fprintf( stderr, "invalid number of points in data segment: %d\n",
              numPoints );
          exit( 1 );
        }
        fprintf( outputFp, PPARAMS "\"--%s\",\"int\",\"%d\",\n", 
            long_options[option_index].name, numPoints );
        break;

      case 'e':
        numSegments = (INT4) atoi( optarg );
        if ( numSegments < 1 )
        {
          fprintf( stderr, "invalid number data segments: %d\n", numSegments );
          exit( 1 );
        }
        fprintf( outputFp, PPARAMS "\"--%s\",\"int\",\"%d\",\n", 
            long_options[option_index].name, numSegments );
        break;

      case 'f':
        ovrlap = (INT4) atoi( optarg );
        if ( ovrlap < 0 )
        {
          fprintf( stderr, "invalid data segment overlap: %d\n", ovrlap );
          exit( 1 );
        }
        fprintf( outputFp, PPARAMS "\"--%s\",\"int\",\"%d\",\n", 
            long_options[option_index].name, ovrlap );
        break;

      case 'g':
        sampleRate = (INT4) atoi( optarg );
        if ( sampleRate < 2 || sampleRate > 16384 || sampleRate % 2 )
        {
          fprintf( stderr, "invalid sample rate: %d\n", sampleRate );
          exit( 1 );
        }
        fprintf( outputFp, PPARAMS "\"--%s\",\"int\",\"%d\",\n", 
            long_options[option_index].name, sampleRate );
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
        fprintf( outputFp, PPARAMS "\"--%s\",\"real\",\"%f\",\n", 
            long_options[option_index].name, fLow );
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
        fprintf( outputFp, PPARAMS "\"--%s\",\"string\",\"%s\",\n", 
            long_options[option_index].name, optarg );
        break;

      case 'k':
        invSpecTrunc = (INT4) atoi( optarg );
        if ( invSpecTrunc < 1 )
        {
          fprintf( stderr, "length of inverse spectrum is less than 1 second"
              ": %d\n", invSpecTrunc );
          exit( 1 );
        }
        fprintf( outputFp, PPARAMS "\"--%s\",\"int\",\"%d\",\n", 
            long_options[option_index].name, invSpecTrunc );
        break;

      case 'l':
        dynRangeExponent = (REAL4) atof( optarg );
        fprintf( outputFp, PPARAMS "\"--%s\",\"real\",\"%f\",\n", 
            long_options[option_index].name, dynRangeExponent );
        break;
          
      case 'm':
        startTemplate = (INT4) atoi( optarg );
        if ( startTemplate < 0 )
        {
          fprintf( stderr, "template bank start index is invalid: %d\n", 
              startTemplate );
          exit( 1 );
        }
        fprintf( outputFp, PPARAMS "\"--%s\",\"int\",\"%d\",\n", 
            long_options[option_index].name, startTemplate );
        break;

      case 'n':
        stopTemplate = (INT4) atoi( optarg );
        if ( stopTemplate < 0 )
        {
          fprintf( stderr, "template bank stop index is invalid: %d\n", 
              stopTemplate );
          exit( 1 );
        }
        fprintf( outputFp, PPARAMS "\"--%s\",\"int\",\"%d\",\n", 
            long_options[option_index].name, stopTemplate );
        break;

      case 'o':
        numChisqBins = (INT4) atoi( optarg );
        if ( numChisqBins < 0 )
        {
          fprintf( stderr, "invalid number of chisq veto bins: %d\n", 
              numChisqBins );
          exit( 1 );
        }
        fprintf( outputFp, PPARAMS "\"--%s\",\"int\",\"%d\",\n", 
            long_options[option_index].name, numChisqBins);
        break;

      case 'p':
        hierDepth = (INT4) atoi( optarg );
        if ( hierDepth < 0 || hierDepth > 1 )
        {
          fprintf( stderr, "invalid hierarchical search depth: %d\n", 
              hierDepth );
          exit( 1 );
        }
        fprintf( outputFp, PPARAMS "\"--%s\",\"int\",\"%d\",\n", 
            long_options[option_index].name, hierDepth );
        break;

      case 'q':
        {
          size_t rhosqlen = strlen( optarg );
          rhosqStr = (char *) LALMalloc( ++rhosqlen );
          memcpy( rhosqStr, optarg, rhosqlen );
        }
        fprintf( outputFp, PPARAMS "\"--%s\",\"string\",\"%s\",\n", 
            long_options[option_index].name, optarg );
        break;

      case 'r':
        {
          size_t chisqlen = strlen( optarg );
          chisqStr = (char *) LALMalloc( ++chisqlen );
          memcpy( chisqStr, optarg, chisqlen );
        }
        fprintf( outputFp, PPARAMS "\"--%s\",\"string\",\"%s\",\n", 
            long_options[option_index].name, optarg );
        break;

      case 's':
        if ( strlen( optarg ) > PROCESS_COMMENT_LEN - 1 )
        {
          fprintf( stderr, "--comment must be less than %d characters",
              PROCESS_COMMENT_LEN );
          exit( 1 );
        }
        else
        {
          strncpy( process_comment, optarg, 
              PROCESS_COMMENT_LEN * sizeof(CHAR) );
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

  /* check band pass option */
  if ( bandPass == 1 )
  {
    fprintf( outputFp, PPARAMS "\"--enable-band-pass\",\"int\",\"1\",\n" );
  }
  else if ( bandPass == 0 )
  {
    fprintf( outputFp, PPARAMS "\"--disable-band-pass\",\"int\",\"1\",\n" );
  }
  else
  {
    fprintf( stderr, "--enable-band-pass or --disable-band-pass "
        "argument must be specified\n" );
    exit( 1 );
  }

  /* check event cluster option */
  if ( eventCluster == 1 )
  {
    fprintf( outputFp, PPARAMS "\"--enable-event-cluster\",\"int\",\"1\",\n" );
  }
  else if ( eventCluster == 0 )
  {
    fprintf( outputFp, PPARAMS "\"--disable-event-cluster\",\"int\",\"1\",\n" );
  }
  else
  {
    fprintf( stderr, "--enable-event-cluster or --disable-event-cluster "
        "argument must be specified\n" );
    exit( 1 );
  }

  if ( enableOutput == 1 )
  {
    fprintf( outputFp, PPARAMS "\"--enable-output\",\"int\",\"1\"\n" );
  }
  else if ( enableOutput == 0 )
  {
    fprintf( outputFp, PPARAMS "\"--disable-output\",\"int\",\"1\"\n" );
  }
  else
  {
    fprintf( stderr, "--enable-output or --disable-output "
        "argument must be specified\n" );
    exit( 1 );
  }

  /* write the process_params table footer */
  fprintf( outputFp, TABLE_FOOTER );


  /*
   *
   * check validity of arguments
   *
   */


  /* check validity of input data time */
  if ( ! gpsStartTimeNS || ! gpsStopTimeNS )
  {
    fprintf( stderr, "gps start and stop time must be specified\n" );
    exit( 1 );
  }

  if ( gpsStopTimeNS <= gpsStartTimeNS )
  {
    fprintf( stderr, "invalid gps time range: "
        "start time: %lld, stop time %lld\n",
        gpsStartTimeNS / 1000000000LL, gpsStopTimeNS / 1000000000LL );
    exit( 1 );
  }

  /* check validity of data length parameters */
  if ( numPoints < 0 )
  {
    fprintf( stderr, "number of points in segment must be specified\n" );
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
    UINT8 gpsChanIntervalNS = gpsStopTimeNS - gpsStartTimeNS;
    UINT8 inputDataLengthNS = (UINT8) inputDataLength * 1000000000LL / 
      (UINT8) sampleRate;

    if ( inputDataLengthNS != gpsChanIntervalNS )
    {
      fprintf( stderr, "length of input data and data chunk do not match\n" );
      fprintf( stderr, "start time: %lld, stop time %lld\n",
        gpsStartTimeNS / 1000000000LL, gpsStopTimeNS / 1000000000LL );
      fprintf( stderr, "gps channel time interval: %lld ns\n"
          "computed input data length: %lld ns\n", 
          gpsChanIntervalNS, inputDataLengthNS );
      exit( 1 );
    }
  }

  /* check filter parameters have been specified */
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
   * do something
   *
   */


  sleep( 2 );

  
  /*
   *
   * close the output file, write the result file and exit sucessfully
   *
   */


  /* close the output file now everything has been written to it */
  fclose( outputFp );

  /* get the job end gps time */
  ticks = time( NULL );
  gmtime_r( &ticks, &(laldate.unixDate) );
  laldate.residualNanoSeconds = 0;
  LAL_CALL(
      LALUTCtoGPS( &status, &process_end_time, &laldate, &accuracy ),
      &status );
  
  /* write the process table to a file */
  if ( ! (resultFp = fopen( resultFile, "w" )) )
  {
    perror( "could not open result file for writing" );
    exit( 1 );
  }
  fprintf( resultFp, LIGO_LW_HEADER );
  fprintf( resultFp, PROCESS_HEADER );
  if ( process_comment[0] )
  {
    fprintf( resultFp, PROCESS_COMMENT_HEADER );
  }
  fprintf( resultFp, PROCESS_STREAM_HEADER );
  fprintf( resultFp, "         \"%s\",\"%s\",\"%s\",%d,"
      "%d,\"%s\",\"%s\",%d,%d,%d,%d,\"%s\",",
      process_program,
      "process_version_str",
      "process_cvs_repository_str",
      500000000,
      process_isonline,
      process_node,
      process_username,
      process_unix_procid,
      process_start_time.gpsSeconds,
      process_end_time.gpsSeconds,
      process_jobid,
      process_domain );

  if ( process_comment[0] )
  {
    fprintf( resultFp, "\"%s\",", process_comment );
  }
  fprintf( resultFp, "\"process:process_id:0\"\n" );
  fprintf( resultFp, TABLE_FOOTER );
  fclose( resultFp );

  if ( enableOutput )
  {
    /* cat the contents of the output file into the result file */
    int outpp, resup;
    unsigned char *buf = (unsigned char *) malloc( 1024 );
    unsigned char *bufptr;
    ssize_t n_read, len, n_write, written;

    if ( (outpp = open( outputFile, O_RDONLY )) < 1 )
    {
      fprintf( stderr, "error opening output file for concatenation\n" );
      exit( 1 );
    }
    if ( (resup = open( resultFile, O_WRONLY|O_APPEND)) < 1 )
    {
      fprintf( stderr, "error opening output file for concatenation\n" );
      exit( 1 );
    }

    while ( 1 )
    {
      /* read from output file, retrying if interrupted */
      do
      {
        n_read = read( outpp, buf, 1024 );
      }
      while ( n_read < 0 && errno == EINTR );

      if ( n_read < 0 )
      {
        perror( "error reading results from output file" );
        exit( 1 );
      }
      
      /* exit the loop on end of file */
      if ( n_read == 0 )
      {
        free( buf );
        close( outpp );
        close( resup );
        break;
      }

      /* write this block out */
      len = n_read;
      bufptr = buf;
      written = 0;
      while ( len > 0 )
      {
        n_write = write( resup, bufptr, len );
        if ( n_write <= 0 )
        {
          if ( n_write == 0 )
          {
            errno = ENOSPC;
          }
	  if ( errno == EINTR )
          {
	    continue;
          }
          perror( "error writing results from output file" );
          exit( 1 );
        }
        written += n_write;
        bufptr += n_write;
        len -= n_write;
      }

      /* check that wh have written all that we have read */
      if ( written != n_read )
      {
        fprintf( stderr, "error reading/writing results from output file:"
            "read = %d, written = %d\n", written, n_read );
        exit( 1 );
      }
    }
  }
  
  /* write the xml footer */
  if ( ! (resultFp = fopen( resultFile, "a" )) )
  {
    perror( "could not open result file for appending" );
    exit( 1 );
  }
  fprintf( resultFp, LIGO_LW_FOOTER );
  fclose( resultFp );

  /* delete the output file as we don't need it any more */
  if ( unlink( outputFile ) )
  {
    perror( "could not delete output file" );
    exit( 1 );
  }

  /* exit sucessfully */
  exit( 0 );
  }
