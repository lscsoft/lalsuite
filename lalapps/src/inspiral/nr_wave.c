/*----------------------------------------------------------------------- 
 * 
 * File Name: nr_wave.c
 *
 * Author: S.Fairhurst, B. Krishnan, L.Santamaria 
 *
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>

#include <lalapps.h>

#include <lal/LALConfig.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALDatatypes.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/AVFactories.h>
#include <lal/NRWaveIO.h>
#include <lal/NRWaveInject.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/Inject.h>

RCSID( "$Id$" );

#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "nr_wave"

/*
 * 
 * USAGE
 *
 */
static void print_usage(char *program)
{
  fprintf(stderr,
      "Usage:  %s [options] [LIGOLW XML input files]\n" \
      "The following options are recognized.  Options not surrounded in [] are\n" \
      "required.\n" \
      "  [--help]                      display this message\n"\
      "  [--verbose]                   print progress information\n"\
      "  [--version]                   print version information and exit\n"\
      "  --injection-file  inj_file    read injection details from inj_file\n"\
      "  --ifo             ifo         ifo for which to generate injections\n"\
      "  --nr-meta-file    meta_file   file containing details of available\n"\
      "                                numerical relativity waveforms\n"\
      "  --nr-data-file    meta_file   before we have the meta file, specify\n"\
      "                                file with numerical relativity waveform\n"\
      "  --gps-start-time  start       start time of output file\n"\
      "  --gps-end-time    end         end time of output file\n"\
      "  --sample-rate     rate        the sample rate used to generate injections\n"\
      "\n", program);
}



/*
 * 
 * MAIN
 *
 */
int main( int argc, char *argv[] )
{
  extern int vrbflg;

  LALStatus             status = blank_status;

  CHAR  *injectionFile = NULL;         /* name of file containing injs   */
  CHAR  *nrMetaFile    = NULL;         /* name of file with nr meta info */
  CHAR  *nrDataFile    = NULL;         /* name of file with nr waveform  */

  CHAR   ifo[LIGOMETA_IFO_MAX];        /* name of ifo                    */
  int gpsStartSec = -1;                /* start time of data             */
  int gpsEndSec = -1;                  /* end time of data               */
  LIGOTimeGPS gpsStartTime = {0,0};    /* start time GPS                 */
  LIGOTimeGPS gpsEndTime = {0,0};      /* end time GPS                   */

  int sampleRate = -1;                 /* output sample rate             */
  int numPoints = 0;                   /* number of data points          */
  int numInjections = 0;

  SimInspiralTable *injections = NULL; /* list of injections to be done  */
  SimInspiralTable    *thisInj = NULL; /* current injection              */

  REAL4TimeSeries        inj;          /* time series of zeros to which we
                                          add injections                 */
  REAL4TimeVectorSeries *strain = NULL;/* h+,hx time series              */
  REAL4TimeSeries       *htData = NULL;/* h(t) data for given detector   */

  /* getopt arguments */
  struct option long_options[] =
  {
    /* these options set a flag */
    {"verbose",             no_argument,   &vrbflg,                   1 },
    /* these options don't set a flag */
    {"gps-start-time",          required_argument, 0,                'a'},
    {"gps-end-time",            required_argument, 0,                'b'},
    {"injection-file",          required_argument, 0,                'f'},
    {"nr-meta-file",            required_argument, 0,                'm'},
    {"nr-data-file",            required_argument, 0,                'd'},
    {"sample-rate",             required_argument, 0,                'r'},
    {"ifo",                     required_argument, 0,                'i'},
    {"help",                    no_argument,       0,                'h'},
    {"version",                 no_argument,       0,                'V'},
    {0, 0, 0, 0}
  };
  int c;

  /* parse the arguments */
  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    size_t optarg_len;

    c = getopt_long_only( argc, argv, 
        "a:b:d:f:i:m:r:V",
        long_options, &option_index );

    /* detect the end of the options */
    if ( c == -1 )
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
          fprintf( stderr, "Error parsing option %s with argument %s\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        break;

      case 'h':
        /* help message */
        print_usage(argv[0]);
        exit( 0 );
        break;

      case 'V':
        /* print version information and exit */
        fprintf( stdout, "Numerical Relativity Waveform Injection Program\n" 
            "CVS Version: " CVS_ID_STRING "\n"
            "CVS Tag: " CVS_NAME_STRING "\n" );
        exit( 0 );
        break;

      case 'a':
        /* set gps start seconds */
        gpsStartSec = atoi( optarg );
        if ( gpsStartSec < 441417609 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "GPS start time is prior to "
              "Jan 01, 1994  00:00:00 UTC:\n"
              "(%d specified)\n",
              long_options[option_index].name, gpsStartSec );
          exit( 1 );
        }
        if ( gpsStartSec > 999999999 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "GPS start time is after " 
              "Sep 14, 2011  01:46:26 UTC:\n"
              "(%d specified)\n", 
              long_options[option_index].name, gpsStartSec );
          exit( 1 );
        }
        gpsStartTime.gpsSeconds = gpsStartSec;
        break;

      case 'b':
        /* set gps end seconds */
        gpsEndSec = atoi( optarg );
        if ( gpsEndSec > 999999999 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "GPS end time is after " 
              "Sep 14, 2011  01:46:26 UTC:\n"
              "(%d specified)\n", 
              long_options[option_index].name, gpsEndSec );
          exit( 1 );
        }
        else if ( gpsEndSec < 441417609 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "GPS end time is prior to " 
              "Jan 01, 1994  00:00:00 UTC:\n"
              "(%d specified)\n", 
              long_options[option_index].name, gpsEndSec );
          exit( 1 );
        }            
        gpsEndTime.gpsSeconds = gpsEndSec;
        break;

      case 'r':
        /* set the sample rate */
        sampleRate = (INT4) atoi( optarg );
        if ( sampleRate < 1 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "sample rate must be a positive integer: "
              "(%d specified) \n", 
              long_options[option_index].name, sampleRate );
          exit( 1 );
        }
        break;

      case 'i':
        {
          /* create storage for the ifo name and copy it */
          memset( ifo, 0, sizeof(ifo) );
          memcpy( ifo, optarg, sizeof(ifo) - 1 );
        }
        break;

      case 'f':
        /* create storage for the injection file name */
        optarg_len = strlen( optarg ) + 1;
        injectionFile = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( injectionFile, optarg, optarg_len );
        break;

      case 'm':
        /* create storage for the injection file name */
        optarg_len = strlen( optarg ) + 1;
        nrMetaFile = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( nrMetaFile, optarg, optarg_len );
        break;

      case 'd':
        /* create storage for the injection file name */
        optarg_len = strlen( optarg ) + 1;
        nrDataFile = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( nrDataFile, optarg, optarg_len );
        break;

      case '?':
        print_usage(argv[0]);
        exit( 1 );
        break;

      default:
        fprintf( stderr, "unknown error while parsing options\n" );
        print_usage(argv[0]);
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

  /*
   *
   * check validity of arguments
   *
   */


  /* check validity of input data time */

  /* start time specified */
  if ( gpsStartSec < 0 )
  {
    fprintf( stderr, "--gps-start-time must be specified\n" );
    exit( 1 );
  }

  /* end time specified */
  if (  gpsEndSec < 0 )
  {
    fprintf( stderr, "--gps-end-time must be specified\n" );
    exit( 1 );
  }

  /* end after start */
  if ( gpsEndSec <= gpsStartSec )
  {
    fprintf( stderr, "invalid gps time range: "
        "start time: %d, end time %d\n",
        gpsStartSec, gpsEndSec );
    exit( 1 );
  }

  /* check that we have injections */
  if ( !injectionFile )
  {
    fprintf( stderr, 
        "--injection-file must be specified\n");
    exit( 1 );
  }

  if ( !nrMetaFile && !nrDataFile )
  {
    fprintf( stderr, 
        "either --nr-meta-file or --nr-data-file must be specified\n");
    exit( 1 );
  }
 
  if ( nrMetaFile && nrDataFile )
  {
    fprintf( stderr, 
        "only one of --nr-meta-file or --nr-data-file must be specified\n");
    exit( 1 );
  }
 


  /*
   *
   * Main Code
   *
   */
     
  /* set up the injData to be zeros of the correct length, to which we will 
   * add the injections */

  memset( &inj, 0, sizeof(REAL4TimeSeries) );

  inj.epoch = gpsStartTime;
  inj.deltaT = 1./sampleRate;
  LALSnprintf( inj.name, LIGOMETA_CHANNEL_MAX * sizeof(CHAR), "%s:STRAIN", 
          ifo );
  numPoints = sampleRate * (gpsEndSec - gpsStartSec);

  /* create data structure */
  inj.data = XLALCreateREAL4Vector( numPoints );
  if ( ! inj.data )
  {
      fprintf( stderr, "error: could not allocate memory for injections" );
      exit( 1 );
  }
  memset( inj.data->data, 0, numPoints * sizeof(REAL4) );



  /* read the injections */
  numInjections = SimInspiralTableFromLIGOLw( &injections, injectionFile,
        gpsStartSec, gpsEndSec );

  if ( numInjections < 0 )  
  {
    fprintf( stderr, "error: cannot read injection file" );
    exit( 1 );
  }

  for( thisInj = injections; thisInj; thisInj = thisInj ->next )
  {
    if ( nrMetaFile )
    {
      fprintf( stderr, "XLALFindNRFile hasn't been implemented yet.\n"
          "Specify the data file using --nr-data-file" );
      exit( 1 );

      /* nrDataFile = XLALFindNRFile( nrMetaFile, thisInj ); */
    }

    /* compute the h+, hx strain averaged over the sky */
    strain = XLALReadNRWave( (thisInj->mass1 + thisInj->mass2), nrDataFile);

    /* compute the h+, hx strain for the given inclination, coalescence phase*/
    strain = XLALOrientNRWave( *strain, thisInj->inclination,
        thisInj->coa_phase );

    /* compute the strain time series for given ifo, sky location */
    htData = XLALCalculateNRStrain( strain, thisInj, ifo, sampleRate);

    /* XXX inject the htData into our injection time stream XXX */
    /*LAL_CALL( LALSSInjectTimeSeries( &status, &inj, htData ), &status );*/
  }
  exit( 0 );
}

