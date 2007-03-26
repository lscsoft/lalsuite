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
#include <lal/FileIO.h>
#include <lal/Units.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/VectorOps.h>


RCSID( "$Id$" );

#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "nr_wave"

static void output_ht( CHAR *fName, REAL4TimeSeries *injData );

extern int vrbflg;

/*
 * 
 * USAGE
 *
 */
static void print_usage( char *program )
{
  fprintf( stderr,
      "Usage:  %s [options]\n"\
      "The following options are recognized.  Options not surrounded in [] are\n"\
      "required.\n"\
      "  [--help]                      display this message\n"\
      "  [--verbose]                   print progress information\n"\
      "  [--version]                   print version information and exit\n"\
      "  --injection-file  inj_file    read injection details from inj_file\n"\
      "  --ifo             ifo         IFO for which to generate injections\n"\
      "  --all-ifos                    create injections for all IFOs\n"\
      "  --nr-meta-file    meta_file   file containing details of available\n"\
      "                                numerical relativity waveforms\n"\
      "  --nr-data-dir     dir         specify directory containing numerical\n"\
      "                                relativity waveforms\n"\
      "  --gps-start-time  start       start time of output file\n"\
      "  --gps-end-time    end         end time of output file\n"\
      "  --sample-rate     rate        the sample rate used to generate injections\n"\
      "  --write-output                write h(t) to an ascii file in NRwave directory\n"\
      "\n", program );
}

/*
 * 
 * MAIN
 *
 */
int main( int argc, char *argv[] )
{
  LALStatus status = blank_status;

  int i;                                 /* loop counter */
  int num_ifos;                          /* number of ifos */

  CHAR *injectionFile = NULL;            /* name of file containing injs   */
  CHAR *nrMetaFile    = NULL;            /* name of file with nr meta info */
  CHAR *nrDataDir     = NULL;            /* name of dir with nr waveform   */
  NRWaveMetaData thisMetaData;           /* single NR wave metadata struct */

  NRWaveCatalog nrCatalog;               /* NR wave metadata struct        */

  CHAR ifo[LIGOMETA_IFO_MAX];            /* name of ifo                    */
  CHAR fileName[FILENAME_MAX];           /* name of output file            */
  CHAR name[LALNameLength];

  int gpsStartSec          = -1;         /* start time of data             */
  int gpsEndSec            = -1;         /* end time of data               */
  LIGOTimeGPS gpsStartTime = {0, 0};     /* start time GPS                 */
  LIGOTimeGPS gpsEndTime   = {0, 0};     /* end time GPS                   */

  int sampleRate    = -1;                /* output sample rate             */
  int numInjections = 0;                 /* number of injections */

  SimInspiralTable *injections = NULL;   /* list of injections to be done  */
  SimInspiralTable *thisInj    = NULL;   /* current injection              */

  REAL4TimeSeries injData;               /* time series of zeros to which
                                            we add injections              */
  REAL4TimeVectorSeries *strain = NULL;  /* h+, hx time series             */
  REAL4TimeSeries *htData;               /* h(t) data for given detector   */

  int writeFlag = 0;                     /* write h(t) to file?            */
  int ifosFlag  = 0;                     /* injections for all ifos?       */

  /* getopt arguments */
  struct option long_options[] =
  {
    /* these options set a flag */
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"write-output",            no_argument,       &writeFlag,        1 },
    {"all-ifos",                no_argument,       &ifosFlag,         1 },
    /* these options don't set a flag */
    {"gps-start-time",          required_argument, 0,                'a'},
    {"gps-end-time",            required_argument, 0,                'b'},
    {"injection-file",          required_argument, 0,                'f'},
    {"nr-meta-file",            required_argument, 0,                'm'},
    {"nr-data-dir",             required_argument, 0,                'd'},
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
        "a:b:d:f:i:m:r:V:W",
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
        /* create storage for the ifo name and copy it */
        optarg_len = strlen( optarg ) + 1;
        memcpy( ifo, optarg, optarg_len );

        /* check for supported ifo */
        if (XLALIFONumber(ifo) == LAL_UNKNOWN_IFO)
        {
          fprintf(stderr, "IFO not recognised: %s\n", ifo);
          exit(1);
        }
        break;

      case 'f':
        /* create storage for the injection file name */
        optarg_len = strlen( optarg ) + 1;
        injectionFile = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( injectionFile, optarg, optarg_len );
        break;

      case 'm':
        /* create storage for the meta file name */
        optarg_len = strlen( optarg ) + 1;
        nrMetaFile = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( nrMetaFile, optarg, optarg_len );
        break;

      case 'd':
        /* create storage for the nr data directory */
        optarg_len = strlen( optarg ) + 1;
        nrDataDir = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( nrDataDir, optarg, optarg_len );
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

  /* check that sample rate has been specified */
  if ( sampleRate < 0 )
  {
    fprintf( stderr, "--sample-rate must be specified\n" );
    exit( 1 );
  }

  /* check that we have injections */
  if ( injectionFile == NULL )
  {
    fprintf( stderr, "--injection-file must be specified\n" );
    exit( 1 );
  }

  /* ifo specified, or all ifos */
  if (( !ifo ) && ( !ifosFlag ))
  {
    fprintf( stderr, "--ifo, or --all-ifos, must be specifed\n" );
    exit( 1 );
  }

  /* metadata file specified */
  if ( nrMetaFile == NULL )
  {
    fprintf( stderr, "--nr-meta-file must be specified\n" );
    exit( 1 );
  }

  /* data directory specified */
  if ( nrDataDir == NULL )
  {
    fprintf( stderr, "--nr-data-dir must be specified\n" );
    exit( 1 );
  }

  /*
   *
   * Main Code
   *
   */

  /* set channel name */
  LALSnprintf( name, LIGOMETA_CHANNEL_MAX * sizeof(CHAR), "%s:STRAIN", ifo );

  /* set up the injData to be zeros of the correct length, to which we will 
   * add the injections */
  injData = *XLALCreateREAL4TimeSeries( name, &gpsStartTime, 0, 1./sampleRate, 
      &lalADCCountUnit, sampleRate * (gpsEndSec - gpsStartSec) );

  /* read the injections */
  numInjections = SimInspiralTableFromLIGOLw( &injections, injectionFile,
      gpsStartSec, gpsEndSec );

  if ( vrbflg )
  {
    fprintf(stdout,"Read %d injection(s) from the file %s\n",
        numInjections, injectionFile);
  }

  if ( numInjections < 0 )  
  {
    fprintf( stderr, "ERROR: Cannot read injection file\n" );
    exit( 1 );
  }

  /* get catalog of numrel waveforms from metadata file */
  LAL_CALL(LALNRDataFind( &status, &nrCatalog, nrDataDir, nrMetaFile ), 
      &status);

  /* get number of ifos to calculate injections for */
  if (ifosFlag)
  {
    num_ifos = LAL_NUM_IFO;
  }
  else
  {
    num_ifos = 1;
  }

  /* loop over ifos */
  for ( i = 0; i < num_ifos; i++ )
  {
    /* get ifo */
    if (ifosFlag)
    {
      XLALReturnIFO( ifo, i );
    }

    if ( vrbflg )
    {
      fprintf( stdout, "Perfroming injections for IFO: %s\n", ifo);
    }

    /* set output filename */
    LALSnprintf( fileName, FILENAME_MAX, "%s-NR_WAVE-%d-%d.dat", 
        ifo, gpsStartSec, gpsEndSec - gpsStartSec);

    /* loop over injections */
    for ( thisInj = injections; thisInj; thisInj = thisInj->next )
    {
      /* find nearest matching numrel waveform */
      XLALFindNRFile( &thisMetaData, &nrCatalog, thisInj, 2, 2);

      if ( vrbflg )
      {
        fprintf(stdout, "Reading the waveform from the file \"%s\"...",
            thisMetaData.filename );
      }

      /* read numrel waveform */
      LAL_CALL(LALReadNRWave(&status, &strain, thisInj->mass1 + thisInj->mass2, 
            thisMetaData.filename), &status);

      if ( vrbflg )
      {
        fprintf(stdout, "done\n");
      }

      if ( vrbflg )
      {
        fprintf(stdout,
            "Generating waveform for inclination = %f, coa_phase = %f\n",
            thisInj->inclination, thisInj->coa_phase );
      }

      /* compute the h+ and hx for given inclination and coalescence phase*/
      strain = XLALOrientNRWave( strain, thisMetaData.mode[0],
          thisMetaData.mode[1], thisInj->inclination, thisInj->coa_phase);

      if ( vrbflg )
      {
        fprintf(stdout,
            "Generating the strain data for the given sky location\n");
      }

      /* compute strain for given sky location */
      htData = XLALCalculateNRStrain( strain, thisInj, ifo, sampleRate);

      /* inject the htData into injection time stream */
      LAL_CALL( LALSSInjectTimeSeries( &status, &injData, htData ), &status );

      /* clear memory for strain */
      XLALDestroyREAL4VectorSequence ( strain->data );
      LALFree(strain);
      strain = NULL;

    } /* end loop over injections */

    /* output injections */
    if ( writeFlag )
    {
      output_ht( fileName, &injData);
    }

  } /* end loop over ifos */

  /* clear memory */
  LALFree(nrCatalog.data);
  LALCheckMemoryLeaks(); 

  exit( 0 );
}


/* function to output h(t) waveform to file */
static void output_ht(CHAR *fileName,
    REAL4TimeSeries *injData)
{
  FILE *htOut = NULL;
  UINT4 i = 0;
  REAL8 time = 0;

  /* open output file */
  if ((htOut = fopen(fileName, "w")) == NULL)
  {
    fprintf(stderr, "ERROR: Unable to open output file \"%s\"\n", fileName);
    exit( 1 );
  }

  /* get gps time */
  time = XLALGPSGetREAL8( &(injData->epoch) );

  if ( vrbflg )
  {
    fprintf( stdout, "Writing output data to %s\n", fileName );
  }

  /* output h(t) waveform */
  for (i = 0; i < injData->data->length; i++)
  {
    fprintf(htOut, "%.6f\t%e\t\n", time, injData->data->data[i]);

    /* increment time */
    time += injData->deltaT;
  }
  return;
}
