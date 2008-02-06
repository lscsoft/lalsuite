/*
 * Copyright (C) 2007, 2008 Badri Krishnan, Chad Hanna, Lucia Santamaria Lara,
 * Robert Adam Mercer, Stephen Fairhurst
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
#include <lal/FrameStream.h>
#include <lal/FindChirp.h>

#include "inspiral.h"

/* cvs information */
RCSID( "$Id$" );
#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "lalapp_mdc_ninja"

/* defines */
#define HISTORY_COMMENT 512

/* function prototypes */
static void print_usage( CHAR *program );
static void output_frame( CHAR *ifo, INT4 gpsStart, INT4 gpsEnd,
    REAL4TimeSeries *injData );
static void output_multi_channel_frame( INT4 num_ifos, INT4 gpsStart,
    INT4 gpsEnd, REAL4TimeSeries *injData[LAL_NUM_IFO] );
static void write_mdc_log_file( CHAR *filename, SimInspiralTable *injections,
    INT4 gps_start, CHAR *set_name );


/* getopt flags */
extern int vrbflg;
INT4 ifosFlag  = 0;
INT4 frameFlag = 0;
INT4 mdcFlag   = 0;
INT4 noNR      = 0;


/* main program entry */
INT4 main( INT4 argc, CHAR *argv[] )
{
  LALStatus status = blank_status;

  /* counters */
  int c;
  INT4 i;
  INT4 num_ifos = 0;

  /* file/directory/set names */
  CHAR *injectionFile = NULL;
  CHAR *setName       = NULL;
  CHAR *mdcFileName   = NULL;

  /* ifo name */
  CHAR *ifo = NULL;

  /* channel */
  CHAR channel[LALNameLength];

  /* start/end times */
  INT4 gpsStartSec          = -1;
  INT4 gpsEndSec            = -1;
  LIGOTimeGPS gpsStartTime  = {0, 0};
  LIGOTimeGPS gpsEndTime    = {0, 0};

  /* injections */
  SimInspiralTable *injections = NULL;
  INT4 numInjections           = 0;

  /* injection waveforms time series */
  INT4 sampleRate = -1;
  REAL4TimeSeries *injData[LAL_NUM_IFO];

  /* response function */
  COMPLEX8FrequencySeries *response = NULL;
  const LALUnit strainPerCount = {0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};

  /* the inspiral pipeline resizes data day 2^dynRange. Set to 1.0 when
   * using as standalone code */
  REAL4 dynRange = 1.0;

  /* getopt arguments */
  struct option long_options[] =
  {
    /* these options set a flag */
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"all-ifos",                no_argument,       &ifosFlag,         1 },
    {"write-frame",             no_argument,       &frameFlag,        1 },
    {"write-mdc-log",           no_argument,       &mdcFlag,          1 },
    {"no-numerical",            no_argument,       &noNR,             1 },
    /* these options don't set a flag */
    {"debug-level",             required_argument, 0,                'D'},
    {"gps-start-time",          required_argument, 0,                'a'},
    {"gps-end-time",            required_argument, 0,                'b'},
    {"injection-file",          required_argument, 0,                'f'},
    {"sample-rate",             required_argument, 0,                'r'},
    {"ifo",                     required_argument, 0,                'i'},
    {"set-name",                required_argument, 0,                'n'},
    {"mdc-log",                 required_argument, 0,                'o'},
    {"help",                    no_argument,       0,                'h'},
    {"version",                 no_argument,       0,                'V'},
    {0, 0, 0, 0}
  };

  /* set default debug level */
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "33" );

  /* parse the arguments */
  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    size_t optarg_len;

    /* parse command line arguments */
    c = getopt_long_only( argc, argv, "D:a:b:f:r:i:n:o:hV",
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
          fprintf( stderr, "Error parsing option '%s' with argument '%s'\n",
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
        fprintf( stdout, "%s - Numerical Relativity MDC Injection Program\n" \
            "CVS Version: %s\nCVS Tag: %s\n", PROGRAM_NAME, CVS_ID_STRING, \
            CVS_NAME_STRING );
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
        ifo = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( ifo, optarg, optarg_len );

        /* check for supported ifo */
        if ( XLALIFONumber( ifo ) == LAL_UNKNOWN_IFO )
        {
          fprintf( stderr, "IFO not recognised: %s\n", ifo );
          exit(1);
        }
        break;

      case 'n':
        /* create storage for the injection set name */
        optarg_len = strlen(optarg) + 1;
        setName = (CHAR *)calloc(optarg_len, sizeof(CHAR));
        memcpy(setName, optarg, optarg_len);
        break;

      case 'o':
        /* create storage for the output mdc log file name */
        optarg_len = strlen(optarg) + 1;
        mdcFileName = (CHAR *)calloc(optarg_len, sizeof(CHAR));
        memcpy(mdcFileName, optarg, optarg_len);
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
        fprintf( stderr, "ERROR: Unknown error while parsing options\n" );
        print_usage( argv[0] );
        exit( 1 );
    }
  }

  if ( optind < argc )
  {
    fprintf( stderr, "ERROR: Extraneous command line arguments:\n" );
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


  if ( frameFlag )
  {
    /* check that sample rate has been specified */
    if ( sampleRate < 0 )
    {
      fprintf( stderr, "ERROR: --sample-rate must be specified\n" );
      exit( 1 );
    }

    /* ifo specified, or all ifos */
    if (( !ifo ) && ( ifosFlag == 0 ))
    {
      fprintf( stderr, "ERROR: --ifo, or --all-ifos, must be specifed\n" );
      exit( 1 );
    }
  }

  if (( frameFlag ) || ( mdcFlag ))
  {
    /* check that we have injections */
    if ( injectionFile == NULL )
    {
      fprintf( stderr, "ERROR: --injection-file must be specified\n" );
      exit( 1 );
    }

    /* start time specified */
    if ( gpsStartSec < 0 )
    {
      fprintf( stderr, "ERROR: --gps-start-time must be specified\n" );
      exit( 1 );
    }

    /* end time specified */
    if ( gpsEndSec < 0 )
    {
      fprintf( stderr, "ERROR: --gps-end-time must be specified\n" );
      exit( 1 );
    }

    /* end after start */
    if ( gpsEndSec <= gpsStartSec )
    {
      fprintf( stderr, "ERROR: Invalid gps time range: "
          "start time: %d, end time %d\n",
          gpsStartSec, gpsEndSec );
      exit( 1 );
    }
  }

  /* mdc log options */
  if ( mdcFlag )
  {
    /* set name */
    if ( setName == NULL )
    {
      fprintf( stderr, "ERROR: --set-name must be specified\n" );
      exit( 1 );
    }

    /* mdc log name */
    if (mdcFileName == NULL )
    {
      fprintf( stderr, "ERROR: --mdc-log must be specified\n" );
      exit( 1 );
    }
  }

  if (( mdcFlag == 0 ) && ( frameFlag == 0 ))
  {
    fprintf( stderr, "Nothing to do, exiting...\n" );
    exit( 1 );
  }

  /*
   *
   * Main Code
   *
   */

  /* read the injections */
  numInjections = SimInspiralTableFromLIGOLw( &injections, injectionFile,
      gpsStartSec, gpsEndSec );

  if ( vrbflg )
  {
    fprintf( stdout, "Read %d injection(s) from the file '%s'\n",
        numInjections, injectionFile );
  }

  if ( numInjections < 0 )
  {
    fprintf( stderr, "ERROR: Cannot read injection file '%s'\n", injectionFile );
    exit( 1 );
  }

  if ( frameFlag )
  {

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

    if (vrbflg)
      fprintf(stdout, "generating unity response...\n");


    /* setup a unity response frequency series */
    response = XLALCreateCOMPLEX8FrequencySeries("", &gpsStartTime, 0,
        1, &strainPerCount, (sampleRate * (gpsEndSec - gpsStartSec))/2 + 1);
    for ( i = 0; i < (INT4)response->data->length; i++)
    {
      response->data->data[i].re = 1.0;
      response->data->data[i].im = 0;
    }

    /* loop over ifos */
    for ( i = 0; i < num_ifos; i++ )
    {
      /* get ifo */
      if ( ifosFlag )
      {
        ifo = (CHAR *) calloc( LIGOMETA_IFO_MAX, sizeof(CHAR));
        XLALReturnIFO( ifo, i );
      }

      if (noNR)
      {
        /* set the channel name */
        LALSnprintf(channel, LALNameLength, "%s:STRAIN", ifo);
        strncpy(injData[i]->name, channel, LALNameLength);

        /* injected specified waveforms */
        LAL_CALL( LALFindChirpInjectSignals( &status, injData[i], injections, response), &status);

        /* reset the channel name to IFO:STRAIN as LALFindChirpInjectSignals()
         * messes with it */
        strncpy(injData[i]->name, channel, LALNameLength);        
      }
      else
      {
        /* now we can finally inject the numerical waveforms */
        LAL_CALL( InjectNumRelWaveforms ( &status, injData[i], injections, ifo, 
  	  		dynRange), &status);
      }

      /* set strain as unit */
      injData[i]->sampleUnits = lalStrainUnit;

    } /* loop over ifos */

    /* output frame */
    if ( ifosFlag )
    {
      output_multi_channel_frame( num_ifos, gpsStartSec, gpsEndSec, injData );
    }
    else
    {
      output_frame( ifo, gpsStartSec, gpsEndSec, injData[0] );
    }
  }

  /* write mdc log */
  if ( mdcFlag )
  {
    write_mdc_log_file(mdcFileName, injections, gpsStartSec, setName);
  }


  if ( injectionFile )
    free( injectionFile );

  if ( ifo )
    free( ifo );

  for ( i = 0; i < num_ifos; i++ )
    XLALDestroyREAL4TimeSeries(injData[i]);

  XLALDestroyCOMPLEX8FrequencySeries(response);

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
      "  --injection-file  inj_file    read injection details from xml sim-inspiral inj_file\n"\
      "  --ifo             ifo         IFO for which to generate injections\n"\
      "  --all-ifos                    create injections for all IFOs\n"\
      "  --gps-start-time  start       start time of output file\n"\
      "  --gps-end-time    end         end time of output file\n"\
      "  --sample-rate     rate        the sample rate used to generate injections\n"\
      "  --write-mdc-log               write an MDC log file\n"\
      "  --set-name        set_name    set the injection set name\n"\
      "  --mdc-log         mdc_log     name of file for MDC log file\n"\
      "  --write-frame                 write h(t) waveform to a frame file\n"\
      "  --no-numerical                the injections are not numerical\n"\
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
  CHAR creator[HISTORY_COMMENT];  

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
    fprintf( stderr, "ERROR: Unrecognised IFO: '%s'\n", ifo );
    exit( 1 );
  }

  /* define frame */
  frame = XLALFrameNew( &injData->epoch, duration, "LIGO", 0, 1,
      detectorFlags );

  /* set creator metadata */
  LALSnprintf(creator, HISTORY_COMMENT, "creator:$Id$");
  XLALFrHistoryAdd(frame, "creator", creator);

  /* add channel to frame */
  XLALFrameAddREAL4TimeSeriesSimData( frame, injData );

  if ( vrbflg )
  {
    fprintf( stdout, "Writing injection to frame: '%s'\n", fname );
  }

  /* write frame */
  if (XLALFrameWrite( frame, fname, 8) != 0)
  {
    fprintf( stderr, "ERROR: Cannot save frame file: '%s'\n", fname );
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
  CHAR creator[HISTORY_COMMENT];

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

  /* set creator metadata */
  LALSnprintf(creator, HISTORY_COMMENT, "creator:$Id$");
  XLALFrHistoryAdd(frame, "creator", creator);

  /* add channels to frame */
  for( i = 0; i < num_ifos; i++ )
  {
    XLALFrameAddREAL4TimeSeriesSimData( frame, injData[i] );
  }

  if (vrbflg)
  {
    fprintf( stdout, "Writing injections to frame: '%s'\n", fname );
  }

  /* write frame */
  if ( XLALFrameWrite( frame, fname, 8 ) != 0 )
  {
    fprintf( stderr, "ERROR: Cannot save frame file: '%s'\n", fname );
    exit( 1 );
  }

  /* clear frame */
  FrameFree( frame );

  return;
}

/* function to write a Burst MDC log file */
static void write_mdc_log_file(CHAR *filename, SimInspiralTable *injections, INT4 gps_start, CHAR *set_name)
{
  /* variables */
  FILE *output;
  SimInspiralTable *thisInj;
  float f_isco;

  /* open output file */
  if ((output = fopen(filename, "w")) == NULL)
  {
    fprintf(stderr, "ERROR: Cannot open MDC log file: '%s'\n", \
        filename);
    exit(1);
  }

  if (vrbflg)
  {
    fprintf(stdout, "Writing MDC log file: '%s'\n", filename);
  }

  /* loop over injections */
  for (thisInj = injections; thisInj; thisInj = thisInj->next)
  {
    /* GravEn_SimID */
    fprintf(output, "file ");
    /* GravEn_Ampl */
    fprintf(output, "1 ");
    /* StartSamp1 */
    fprintf(output, "0 ");
    /* StartSamp2 */
    fprintf(output, "0 ");
    /* Internal_x */
    fprintf(output, "%g ", cos(thisInj->inclination));
    /* Internal_phi */
    fprintf(output, "%g ", thisInj->coa_phase);
    /* External_x */
    fprintf(output, "%g ", cos(thisInj->latitude - LAL_PI_2));
    /* External_phi */
    fprintf(output, "%g ", thisInj->longitude);
    /* External_psi */
    fprintf(output, "%g ", thisInj->polarization);
    /* FrameGPS */
    fprintf(output, "%d ", gps_start);
    /* SimStartGPS */
    fprintf(output, "%d.%09d ", thisInj->h_end_time.gpsSeconds, thisInj->h_end_time.gpsNanoSeconds);
    /* SimName */
    fprintf(output, "%s ", set_name);
    /* SimHpHp */
    fprintf(output, "0 ");
    /* SimHcHc */
    fprintf(output, "0 ");
    /* SimHpHc */
    fprintf(output, "0 ");
    /* GEO GPS F_+ F_x eff_dist */
    fprintf(output, "GEO %d.%09d 1 0 %g ", thisInj->g_end_time.gpsSeconds,
        thisInj->g_end_time.gpsNanoSeconds, thisInj->eff_dist_g);
    /* H1 GPS F_+ F_x eff_dist */
    fprintf(output, "H1 %d.%09d 1 0 %g ", thisInj->h_end_time.gpsSeconds,
        thisInj->h_end_time.gpsNanoSeconds, thisInj->eff_dist_h);
    /* H2 GPS F_+ F_x eff_dist */
    fprintf(output, "H2 %d.%09d 1 0 %g ", thisInj->h_end_time.gpsSeconds,
        thisInj->h_end_time.gpsNanoSeconds, thisInj->eff_dist_h);
    /* L1 GPS F_+ F_x eff_dist */
    fprintf(output, "L1 %d.%09d 1 0 %g ", thisInj->l_end_time.gpsSeconds,
        thisInj->l_end_time.gpsNanoSeconds, thisInj->eff_dist_l);
    /* TAMA GPS F_+ F_x eff_dist */
    fprintf(output, "TAMA %d.%09d 1 0 %g ", thisInj->t_end_time.gpsSeconds,
        thisInj->t_end_time.gpsNanoSeconds, thisInj->eff_dist_t);
    /* VIRGO GPS F_+ F_x eff_dist */
    fprintf(output, "VIRGO %d.%09d 1 0 %g ", thisInj->v_end_time.gpsSeconds,
        thisInj->v_end_time.gpsNanoSeconds, thisInj->eff_dist_v);

    /* calculate frequency of innermost stable circulat orbit */
    /* taken from: gr-qc/0612100 */
    f_isco = 205 * (20 / (thisInj->mass1 + thisInj->mass2));

    /* numerical relativity specific parameters */
    fprintf(output, "insp ");
    fprintf(output, "distance %g ", thisInj->distance);
    fprintf(output, "mass1 %g ", thisInj->mass1);
    fprintf(output, "mass2 %g ", thisInj->mass2);
    fprintf(output, "mchirp %g ", thisInj->mchirp);
    fprintf(output, "spin1 %g %g %g ", thisInj->spin1x, thisInj->spin1y, thisInj->spin1z);
    fprintf(output, "spin2 %g %g %g ", thisInj->spin2x, thisInj->spin2y, thisInj->spin2z);
    fprintf(output, "freq %g %g\n", thisInj->f_lower, f_isco);
  }

  /* close output file */
  fclose(output);
}


