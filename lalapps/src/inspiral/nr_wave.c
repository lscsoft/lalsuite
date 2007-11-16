/*
 * Copyright (C) 2007 Badri Krishnan, Chad Hanna, Lucia Santamaria Lara, Robert Adam Mercer, Stephen Fairhurst
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA  02111-1307  USA
 *
 * Revision: $Id$
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
#include <lal/LIGOMetadataUtils.h>
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
#include <lal/LALDetectors.h>
#include <lal/LALFrameIO.h>


RCSID( "$Id$" );

#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "nr_wave"

/* function prototypes */
static void print_usage( CHAR *program );
static void output_frame( CHAR *ifo, INT4 gpsStart, INT4 gpsEnd,
    REAL4TimeSeries *injData );
static void output_multi_channel_frame( INT4 num_ifos, INT4 gpsStart,
    INT4 gpsEnd, REAL4TimeSeries *injData[LAL_NUM_IFO] );


/* verbose flag */
extern int vrbflg;


/* main program entry */
INT4 main( INT4 argc, CHAR *argv[] )
{
  LALStatus status = blank_status;

  INT4 i;                                 /* loop counter                   */
  INT4 num_ifos;                          /* number of ifos                 */

  INT4 modeLlo = -1;                     /* lowest value of l to inject    */
  INT4 modeLhi = -1;                     /* highest values of l to inject  */

  CHAR *injectionFile = NULL;            /* name of file containing injs   */
  CHAR *nrMetaFile    = NULL;            /* name of file with nr meta info */
  CHAR *nrDataDir     = NULL;            /* name of dir with nr waveform   */

  NumRelInjectParams nrPar;
  NRWaveCatalog nrCatalog;               /* NR wave metadata struct        */

  CHAR ifo[LIGOMETA_IFO_MAX];            /* name of ifo                    */

  INT4 gpsStartSec          = -1;         /* start time of data             */
  INT4 gpsEndSec            = -1;         /* end time of data               */
  LIGOTimeGPS gpsStartTime = {0, 0};     /* start time GPS                 */
  LIGOTimeGPS gpsEndTime   = {0, 0};     /* end time GPS                   */

  INT4 sampleRate    = -1;                /* output sample rate             */
  INT4 numInjections = 0;                 /* number of injections           */

  SimInspiralTable *injections = NULL;   /* list of injections to be done  */

  REAL4TimeSeries *injData[LAL_NUM_IFO]; /* time series of zeros to which
                                            we add injections              */

  INT4 ifosFlag  = 0;                     /* injections for all ifos?       */
  INT4 frameFlag = 0;                     /* write h(t) to a frame?         */
  int c;

  /* default debug level */
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "5" );

  /* getopt arguments */
  struct option long_options[] =
  {
    /* these options set a flag */
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"all-ifos",                no_argument,       &ifosFlag,         1 },
    {"write-frame",             no_argument,       &frameFlag,        1 },
    /* these options don't set a flag */
    {"gps-start-time",          required_argument, 0,                'a'},
    {"gps-end-time",            required_argument, 0,                'b'},
    {"injection-file",          required_argument, 0,                'f'},
    {"nr-meta-file",            required_argument, 0,                'm'},
    {"nr-data-dir",             required_argument, 0,                'd'},
    {"sample-rate",             required_argument, 0,                'r'},
    {"ifo",                     required_argument, 0,                'i'},
    {"modeL-lo",                required_argument, 0,                'L'},
    {"modeL-hi",                required_argument, 0,                'H'},
    {"debug-level",             required_argument, 0,                'D'},
    {"help",                    no_argument,       0,                'h'},
    {"version",                 no_argument,       0,                'V'},
    {0, 0, 0, 0}
  };


  /* parse the arguments */
  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    size_t optarg_len;

    c = getopt_long_only( argc, argv, "a:b:f:m:d:r:i:L:H:hV",
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
        print_usage( argv[0] );
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
        if ( XLALIFONumber( ifo ) == LAL_UNKNOWN_IFO )
        {
          fprintf( stderr, "IFO not recognised: %s\n", ifo );
          exit(1);
        }
        break;

      case 'L':
        /* set lower bound of l */
        modeLlo = (INT4) atoi( optarg );
        if ( modeLlo < 2 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "l value must be a greater than 1: "
              "(%d specified) \n",
              long_options[option_index].name, modeLlo );
          exit( 1 );
        }
        break;

      case 'H':
        /* set lower bound of l */
        modeLhi = (INT4) atoi( optarg );
        if ( modeLhi < 2 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "l value must be a greater than 1: "
              "(%d specified) \n",
              long_options[option_index].name, modeLhi );
          exit( 1 );
        }
        break;

      case 'D':
        /* set debug level */
        set_debug_level( optarg );
        break;

      case '?':
        print_usage( argv[0] );
        exit( 1 );
        break;

      default:
        fprintf( stderr, "unknown error while parsing options\n" );
        print_usage( argv[0] );
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
  if ( gpsEndSec < 0 )
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

  /* lowest value of l */
  if ( modeLlo == -1 )
  {
    fprintf( stderr, "--modeL-lo must be specified\n" );
    exit( 1 );
  }

  /* highest value of l */
  if ( modeLhi == -1 )
  {
    fprintf( stderr, "--modeL-hi must be specified\n" );
    exit( 1 );
  }

  /*
   *
   * Main Code
   *
   */

  /* get number of ifos */
  if ( ifosFlag )
    num_ifos = LAL_NUM_IFO;
  else
    num_ifos = 1;

  /* setup the injection time series to be zeros of the correct length */
  for ( i = 0; i < num_ifos; i++ )
  {
    injData[i] = XLALCreateREAL4TimeSeries( "", &gpsStartTime, 0, 1./sampleRate,
        &lalADCCountUnit, sampleRate * (gpsEndSec - gpsStartSec) );
    memset( injData[i]->data->data, 0.0, injData[i]->data->length * sizeof(REAL4) );
  }

  /* read the injections */
  numInjections = SimInspiralTableFromLIGOLw( &injections, injectionFile,
      gpsStartSec, gpsEndSec );

  if ( vrbflg )
  {
    fprintf( stdout, "Read %d injection(s) from the file %s\n",
        numInjections, injectionFile );
  }

  if ( numInjections < 0 )
  {
    fprintf( stderr, "ERROR: Cannot read injection file\n" );
    exit( 1 );
  }

  /* get catalog of numrel waveforms from metadata file */
  LAL_CALL( LALNRDataFind( &status, &nrCatalog, nrDataDir, nrMetaFile ),
      &status );

  /* set parameters */
  nrPar.modeLlo = modeLlo;
  nrPar.modeLhi = modeLhi;
  nrPar.nrCatalog = &nrCatalog;

  for ( i = 0; i < num_ifos; i++ )
  {
    /* get ifo */
    if ( ifosFlag )
    {
      XLALReturnIFO( ifo, i );
    }

    /* set ifo */
    nrPar.ifo = ifo;

    /* perform injection */
    LAL_CALL( LALDriveNRInject( &status, injData[i], injections, &nrPar), &status );

    /* set strain as unit */
    injData[i]->sampleUnits = lalStrainUnit;
  }

  /* output frame */
  if ( frameFlag )
  {
    if ( ifosFlag )
    {
      output_multi_channel_frame( num_ifos, gpsStartSec, gpsEndSec, injData );
    }
    else
    {
      output_frame( ifo, gpsStartSec, gpsEndSec, injData[0] );
    }
  }

  /* clear memory */
  LALFree( nrCatalog.data );

  if ( injectionFile )
    free( injectionFile );
  if ( nrMetaFile )
    free( nrMetaFile );
  if ( nrDataDir )
    free( nrDataDir );

  for ( i = 0; i < num_ifos; i++ )
  {
    if ( injData[i]->data->data )
      LALFree( injData[i]->data->data );
    if ( injData[i]->data )
      LALFree( injData[i]->data );
    if ( injData[i] )
      LALFree( injData[i] );
  }

  while ( injections )
  {
    SimInspiralTable *thisInj = NULL;
    thisInj = injections;
    injections = injections->next;
    LALFree( thisInj );
  }

  LALCheckMemoryLeaks();

  exit( 0 );
}


/* function to display program usgae */
static void print_usage( CHAR *program )
{
  fprintf( stderr,
      "Usage:  %s [options]\n"\
      "The following options are recognized.  Options not surrounded in [] are\n"\
      "required.\n"\
      "  [--help]                      display this message\n"\
      "  [--verbose]                   print progress information\n"\
      "  [--version]                   print version information and exit\n"\
      "  --debug-level     lvl         set debug level to 'lvl'\n"\
      "  --injection-file  inj_file    read injection details from inj_file\n"\
      "  --ifo             ifo         IFO for which to generate injections\n"\
      "  --all-ifos                    create injections for all IFOs\n"\
      "  --nr-meta-file    meta_file   file containing details of available\n"\
      "                                numerical relativity waveforms\n"\
      "  --nr-data-dir     dir         specify directory containing numerical\n"\
      "                                relativity waveforms\n"\
      "  --gps-start-time  start       start time of output file\n"\
      "  --gps-end-time    end         end time of output file\n"\
      "  --modeL-lo        lo          lowest value of l to inject\n"\
      "  --modeL-hi        hi          highest value of l to inject\n"\
      "  --sample-rate     rate        the sample rate used to generate injections\n"\
      "  --write-frame                 write h(t) waveform to a frame file\n"\
      "\n", program );
}


/* function to output h(t) waveform in a frame file */
static void output_frame(CHAR *ifo,
    INT4 gpsStart,
    INT4 gpsEnd,
    REAL4TimeSeries *injData)
{
  CHAR fname[FILENAME_MAX];
  INT4 duration;
  INT4 detectorFlags;
  FrameH *frame;

  /* get frame filename */
  duration = gpsEnd - gpsStart;
  LALSnprintf( fname, FILENAME_MAX, "%s-NR_WAVE-%d-%d.gwf", ifo, gpsStart,
      duration );

  /* set detector flags */
  if ( strncmp( ifo, "H2", 2 ) == 0 )
    detectorFlags = LAL_LHO_2K_DETECTOR_BIT;
  else if ( strncmp( ifo, "H1", 2 ) == 0 )
    detectorFlags = LAL_LHO_4K_DETECTOR_BIT;
  else if ( strncmp( ifo, "L1", 2 ) == 0 )
    detectorFlags = LAL_LLO_4K_DETECTOR_BIT;
  else if ( strncmp( ifo, "G1", 2 ) == 0 )
    detectorFlags = LAL_GEO_600_DETECTOR_BIT;
  else if ( strncmp( ifo, "V1", 2 ) == 0 )
    detectorFlags = LAL_VIRGO_DETECTOR_BIT;
  else if ( strncmp( ifo, "T1", 2 ) == 0 )
    detectorFlags = LAL_TAMA_300_DETECTOR_BIT;
  else
  {
    fprintf( stderr, "ERROR: Unrecognised IFO: %s\n", ifo );
    exit( 1 );
  }

  /* define frame */
  frame = XLALFrameNew( &injData->epoch, duration, "LIGO", 0, 1,
      detectorFlags );

  /* add channel to frame */
  XLALFrameAddREAL4TimeSeriesSimData( frame, injData );

  if ( vrbflg )
  {
    fprintf( stdout, "Writing injection to frame: %s\n", fname );
  }

  /* write frame */
  if (XLALFrameWrite( frame, fname, 8) != 0)
  {
    fprintf( stderr, "ERROR: Cannot save frame file: %s\n", fname );
    exit( 1 );
  }

  /* clear frame */
  FrameFree( frame );

  return;
}


/* write injections for all ifos into a single frame */
static void output_multi_channel_frame(INT4 num_ifos,
    INT4 gpsStart,
    INT4 gpsEnd,
    REAL4TimeSeries *injData[LAL_NUM_IFO])
{
  CHAR fname[FILENAME_MAX];
  INT4 duration;
  INT4 detectorFlags;
  FrameH *frame;
  INT4 i;

  /* get frame filename */
  duration = gpsEnd - gpsStart;
  LALSnprintf( fname, FILENAME_MAX, "GHLTV-NR_WAVE-%d-%d.gwf", gpsStart,
      duration );

  /* set detector flags */
  detectorFlags = LAL_GEO_600_DETECTOR_BIT | LAL_LHO_4K_DETECTOR_BIT |
    LAL_LHO_2K_DETECTOR_BIT | LAL_LLO_4K_DETECTOR_BIT |
    LAL_TAMA_300_DETECTOR_BIT | LAL_VIRGO_DETECTOR_BIT;

  /* define frame */
  frame = XLALFrameNew( &(injData[0])->epoch, duration, "LIGO", 0, 1,
      detectorFlags );

  /* add channels to frame */
  for( i = 0; i < num_ifos; i++ )
  {
    XLALFrameAddREAL4TimeSeriesSimData( frame, injData[i] );
  }

  if (vrbflg)
  {
    fprintf( stdout, "Writing injections to frame: %s\n", fname );
  }

  /* write frame */
  if ( XLALFrameWrite( frame, fname, 8 ) != 0 )
  {
    fprintf( stderr, "ERROR: Cannot save frame file: %s\n", fname );
    exit( 1 );
  }

  /* clear frame */
  FrameFree( frame );

  return;
}
