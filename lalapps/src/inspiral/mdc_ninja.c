/*
 * Copyright (C) 2007-2009 Badri Krishnan, Chad Hanna, Lucia Santamaria Lara,
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
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>

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
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/Inject.h>
#include <lal/FileIO.h>
#include <lal/Units.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/VectorOps.h>
#include <lal/LALDetectors.h>
#include <lal/LALFrameIO.h>
#include <lal/LALFrStream.h>
#include <lal/FindChirp.h>
#include <lal/Random.h>
#include <lal/LALNoiseModels.h>
#include <lal/Date.h>

#include <LALAppsVCSInfo.h>

#include "inspiral.h"

/* cvs information */
#define PROGRAM_NAME "lalapp_mdc_ninja"

/* defines */
#define HISTORY_COMMENT 512

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/* function prototypes */
static void print_usage( CHAR *program );
static void output_frame( CHAR *ifo, INT4 gpsStart, INT4 gpsEnd,
    REAL4TimeSeries *injData, CHAR *frameType, CHAR *outDir );
static void output_frame_real8( CHAR *ifo, INT4 gpsStart, INT4 gpsEnd,
    REAL8TimeSeries *injData, CHAR *frameType, CHAR *outDir );
static void output_multi_channel_frame( INT4 num_ifos, INT4 gpsStart,
    INT4 gpsEnd, REAL4TimeSeries *injData[LAL_NUM_IFO], CHAR *frameType, CHAR *outDir );
static void output_multi_channel_frame_real8( INT4 num_ifos, INT4 gpsStart,
    INT4 gpsEnd, REAL8TimeSeries *injData[LAL_NUM_IFO], CHAR *frameType, CHAR *outDir );
static void write_mdc_log_file( CHAR *filename, SimInspiralTable *injections,
    INT4 gps_start, CHAR *set_name );
static int get_spectrum(REAL8Sequence *spectrum, InterferometerNumber ifoNumber,
    REAL8 deltaF, REAL8 strainHighPassFreq, REAL8 dynRange);
static void add_colored_noise(LALStatus *status, REAL4TimeSeries *chan, INT4 ifoNumber,
    RandomParams *randParams, REAL8 dynRange, REAL8 strainHighpassFreq);

/* getopt flags */
extern int vrbflg;
INT4 ifosFlag   = 0;
INT4 frameFlag  = 0;
INT4 mdcFlag    = 0;
INT4 addNoise   = 0;
INT4 doingREAL8 = 0;

/* main program entry */
INT4 main( INT4 argc, CHAR *argv[] )
{
  LALStatus status = blank_status;

  /* counters */
  int c;
  INT4 i;
  INT4 num_ifos = 0;
  int mkdir_result;

  /* file/directory/set names */
  CHAR *injectionFile = NULL;
  CHAR *setName       = NULL;
  CHAR *mdcFileName   = NULL;
  CHAR *frameType     = NULL;
  CHAR *outDir        = NULL;

  /* ifo name */
  CHAR *ifo = NULL;

  /* channel */
  CHAR channel[LALNameLength];

  /* start/end times */
  INT4 gpsStartSec          = -1;
  INT4 gpsEndSec            = -1;
  INT4 injectWindow         = 0;
  LIGOTimeGPS gpsStartTime  = {0, 0};
  LIGOTimeGPS UNUSED gpsEndTime = {0, 0};

  REAL8 freqLowCutoff = -1;
  REAL8 strainLowPassFreq = -1;
  REAL8 snrLow = -1;
  REAL8 snrHigh = -1;

  /* injections */
  SimInspiralTable *injections = NULL;
  INT4 numInjections           = 0;
  CHAR *injectionType = NULL;
  CHAR *fnameOutXML = NULL;

  /* injection waveforms time series */
  INT4 sampleRate = -1;
  REAL4TimeSeries *injDataREAL4[LAL_NUM_IFO];
  REAL8TimeSeries *injDataREAL8[LAL_NUM_IFO];

  /* random seed for producing gaussian noise if required */
  RandomParams  *randParams = NULL;
  INT4  randomSeed = 0;

  /* response function */
  COMPLEX8FrequencySeries *response = NULL;
  const LALUnit strainPerCount = {0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};

  /* the inspiral pipeline resizes data day 2^dynRange. Set to 1.0 when
   * using as standalone code */
  REAL4 dynRange = 1.0;

  /* PSDs and low-frequency cutoffs for SNR calculations */
  CHAR *ligoPsdFile     = NULL;
  CHAR *virgoPsdFile    = NULL;
  REAL8FrequencySeries *ligoPsd  = NULL;
  REAL8FrequencySeries *virgoPsd = NULL;
  REAL8 ligoSnrLowFreq  = 0;
  REAL8 virgoSnrLowFreq = 0;
        
  /* getopt arguments */
  struct option long_options[] =
  {
    /* these options set a flag */
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"all-ifos",                no_argument,       &ifosFlag,         1 },
    {"write-frame",             no_argument,       &frameFlag,        1 },
    {"write-mdc-log",           no_argument,       &mdcFlag,          1 },
    {"simulate-noise",          no_argument,       &addNoise,         1 },
    {"double-precision",        no_argument,       &doingREAL8,       1 },
    /* these options don't set a flag */
    {"injection-type",          required_argument, 0,                'T'},
    {"gps-start-time",          required_argument, 0,                'a'},
    {"gps-end-time",            required_argument, 0,                'b'},
    {"injection-file",          required_argument, 0,                'f'},
    {"sample-rate",             required_argument, 0,                'r'},
    {"ifo",                     required_argument, 0,                'i'},
    {"frame-type",              required_argument, 0,                't'},
    {"set-name",                required_argument, 0,                'n'},
    {"mdc-log",                 required_argument, 0,                'o'},
    {"freq-low-cutoff",         required_argument, 0,                'l'},
    {"strain-lowpass-freq",     required_argument, 0,                'L'},
    {"snr-low",                 required_argument, 0,                's'},
    {"snr-high",                required_argument, 0,                'S'},
    {"ligo-psd-file",           required_argument, 0,                'c'},
    {"ligo-low-freq-cutoff",    required_argument, 0,                'e'},
    {"virgo-psd-file",          required_argument, 0,                'g'},
    {"virgo-low-freq-cutoff",   required_argument, 0,                'j'},
    {"out-xml-file",            required_argument, 0,                'O'},
    {"fr-out-dir",              required_argument, 0,                'd'},
    {"inject-window",           required_argument, 0,                'I'},
    {"help",                    no_argument,       0,                'h'},
    {"version",                 no_argument,       0,                'V'},
    {0, 0, 0, 0}
  };

  lal_errhandler = LAL_ERR_EXIT;

  XLALSetSilentErrorHandler();

  /* parse the arguments */
  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    size_t optarg_len;

    /* parse command line arguments */
    c = getopt_long_only( argc, argv, "T:a:b:f:r:i:I:t:n:o:l:L:s:S:c:e:f:g:O:d:hV",
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
        fprintf( stdout, "%s - Numerical Relativity MDC Injection Program\n", PROGRAM_NAME);
        XLALOutputVersionString(stderr, 0);
        exit( 0 );
        break;

      case 'T':
        /* create storage for the injection type */
        optarg_len = strlen(optarg) + 1;
        injectionType = (CHAR *)calloc(optarg_len, sizeof(CHAR));
        memcpy(injectionType, optarg, optarg_len);
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
        gpsStartTime.gpsSeconds = gpsStartSec;
        break;

      case 'b':
        /* set gps end seconds */
        gpsEndSec = atoi( optarg );
        if ( gpsEndSec < 441417609 )
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

      case 'I':
        /* set gps end seconds */
        injectWindow = atoi( optarg );
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

      case 't':
        /* create storage for the frame type */
        optarg_len = strlen( optarg ) + 1;
        frameType = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( frameType, optarg, optarg_len );

      case 'n':
        /* create storage for the injection set name */
        optarg_len = strlen(optarg) + 1;
        setName = (CHAR *)calloc(optarg_len, sizeof(CHAR));
        memcpy(setName, optarg, optarg_len);
        break;

      case 'o':
        /* create storage for the output mdc log file name */
        optarg_len = strlen(optarg) + 1;
        mdcFileName = (CHAR *)calloc(1, optarg_len*sizeof(CHAR));
        memcpy(mdcFileName, optarg, optarg_len);
        break;

      case 'l':
        /* set lower cutoff frequency */
        freqLowCutoff = atof(optarg);
        if (freqLowCutoff < 0 )
        {
          fprintf(stderr, "invalid argument to --%s:\n"
              "lower cutoff frequency must be positive: "
              "(%f specified) \n",
              long_options[option_index].name, freqLowCutoff);
          exit(1);
        }
        break;

      case 'L':
        /* set low-pass cutoff frequency for producing noise */
        strainLowPassFreq = atof(optarg);
        if (strainLowPassFreq < 0 )
        {
          fprintf(stderr, "invalid argument to --%s:\n"
              "low pass frequency must be positive: "
              "(%f specified) \n",
              long_options[option_index].name, strainLowPassFreq);
          exit(1);
        }
        break;

      case 's':
        /* set low-pass cutoff frequency for producing noise */
        snrLow = atof(optarg);
        break;

      case 'S':
        /* set low-pass cutoff frequency for producing noise */
        snrHigh = atof(optarg);
        break;

      case 'c':
        /* Specify a file to use as the LIGO psd */
        optarg_len = strlen( optarg ) + 1;
        ligoPsdFile = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( ligoPsdFile, optarg, optarg_len );
        break;
      case 'e':
        /* Specify the low-frequency cutoff for LIGO SNR integrals */
        ligoSnrLowFreq = atof(optarg);
        if (ligoSnrLowFreq < 0 )
        {
          fprintf(stderr, "invalid argument to --%s:\n"
              "LIGO low frequency must be positive: "
              "(%f specified) \n",
              long_options[option_index].name, ligoSnrLowFreq);
          exit(1);
        }
        break;
      case 'g':
        /* Specify a file to use as the Virgo psd */
        optarg_len = strlen( optarg ) + 1;
        virgoPsdFile = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( virgoPsdFile, optarg, optarg_len );
        break;
      case 'j':
        /* Specify the low-frequency cutoff for Virgo SNR integrals */
        virgoSnrLowFreq = atof(optarg);
        if (virgoSnrLowFreq < 0 )
        {
          fprintf(stderr, "invalid argument to --%s:\n"
              "Virgo low frequency must be positive: "
              "(%f specified) \n",
              long_options[option_index].name, virgoSnrLowFreq);
          exit(1);
        }
        break;
        break;

      case 'O':
        /* set output xml file name */
        optarg_len = strlen(optarg) + 1;
        fnameOutXML = (CHAR *)calloc(1,optarg_len*sizeof(CHAR));
        memcpy(fnameOutXML, optarg, optarg_len);
        break;

      case 'd':
        /* set frame output directory */
        optarg_len = strlen(optarg) + 1;
        outDir = (CHAR *)calloc(1, optarg_len * sizeof(CHAR));
        memcpy(outDir, optarg, optarg_len);
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

  if (outDir != NULL)
  {
    /* create output directory if needed */
    errno = 0;
    mkdir_result = mkdir(outDir, S_IRWXU | S_IRWXG | S_IRWXO);

    if ( (mkdir_result == -1) && (errno != EEXIST) )
    {
      fprintf(stderr, "unable to create output directory '%s'\n", outDir);
      exit( 1 );
    }
  }
  else
  {
    /* use current directory for output */
    outDir = (CHAR *)calloc(1, 2 * sizeof(CHAR));
    memcpy(outDir, ".", 2);
  }

  if (( mdcFlag == 0 ) && ( frameFlag == 0 ))
  {
    fprintf( stderr, "Nothing to do, exiting...\n" );
    exit( 1 );
  }

  if (addNoise)
  {
    if ( strainLowPassFreq < 0 )
    {
      fprintf( stderr, "ERROR: --strain-lowpass-freq must be specified\n" );
      exit( 1 );
    }
  }

  if ( frameFlag )
  {
    if (injectionType == NULL)
    {
      fprintf(stderr, "ERROR: --injection-type must be specified\n");
      exit(1);
    }
    else
    {
      if (!((strncmp(injectionType, "NR", strlen(injectionType) + 1) == 0) ||
            (strncmp(injectionType, "approximant", strlen(injectionType) + 1) == 0)))
      {
        fprintf( stderr, "ERROR: --injection-type must be 'NR', or 'approximant'\n");
        exit(1);
      }
    }

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

    /* check that frameType has been specified */
    if ( !frameType )
    {
      fprintf( stderr, "ERROR: --frame-type must be specified\n");
      exit( 1 );
    }
  }

  if ((frameFlag) && (strncmp(injectionType, "NR", strlen(injectionType) + 1) == 0))
  {
    if ( freqLowCutoff < 0 )
    {
      fprintf( stderr, "ERROR: --freq-low-cutoff must be specified\n" );
      exit( 1 );
    }

    if (snrLow < 0)
    {
      fprintf(stderr, "ERROR: --snr-low must be be specified\n");
      exit(1);
    }

    if (snrHigh < 0)
    {
      fprintf(stderr, "ERROR: --snr-high must be be specified\n");
      exit(1);
    }

    if (snrLow > snrHigh)
    {
      fprintf( stderr, "ERROR: --snr-low must be lesser than --snr-high\n");
      exit(1);
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

  /* PSD options */
  if ( ligoPsdFile != NULL || ligoSnrLowFreq != 0.0 )
  {
      if ( ligoPsdFile == NULL || ligoSnrLowFreq == 0.0 )
      {
        fprintf( stderr, "ERROR: both --ligo-psd-file and --ligo-low-freq-cutoff must be specified if either is\n" );
        exit( 1 );
      }

      XLALPsdFromFile(&ligoPsd, ligoPsdFile);
  }

  if ( virgoPsdFile != NULL || virgoSnrLowFreq != 0.0 )
  {
      if ( virgoPsdFile == NULL || virgoSnrLowFreq == 0.0 )
      {
        fprintf( stderr, "ERROR: both --virgo-psd-file and --virgo-low-freq-cutoff must be specified if either is\n" );
        exit( 1 );
      }
      XLALPsdFromFile(&virgoPsd, virgoPsdFile);
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

  /*
   *
   * Main Code
   *
   */

  /* read the injections */
  numInjections = SimInspiralTableFromLIGOLw( &injections, injectionFile,
      gpsStartSec-injectWindow, gpsEndSec+injectWindow );

  if ( vrbflg )
  {
    fprintf( stdout, "Read %d injection(s) from the file '%s'\n",
        numInjections, injectionFile );
  }

  if ( numInjections == 0 )
  {
    fprintf( stderr, "ERROR: No injections in specified time\n");
    exit( 1 );
  }

  if ( numInjections < 0 )
  {
    fprintf( stderr, "ERROR: Cannot read injection file '%s'\n", injectionFile );
    exit( 1 );
  }

  /* create random parameters to be used later for adding
     simulated noise if required */
  if ( addNoise )
    LAL_CALL( LALCreateRandomParams( &status, &randParams, randomSeed ), &status );

  if ( frameFlag )
  {

    /* get number of ifos */
    if ( ifosFlag )
      num_ifos = LAL_NUM_IFO;
    else
      num_ifos = 1;

    /* setup the injection time series to be zeros of the correct length */
    if (doingREAL8)
    {
      if (addNoise)  {
        fprintf( stderr, "ERROR: --simulate-noise is only available for single-precision" );
        exit( 1 );
      }

      for ( i = 0; i < num_ifos; i++ )
      {
        injDataREAL8[i] = XLALCreateREAL8TimeSeries( "", &gpsStartTime, 0, 1./sampleRate,
            &lalADCCountUnit, sampleRate * (gpsEndSec - gpsStartSec) );
        memset( injDataREAL8[i]->data->data, 0.0, injDataREAL8[i]->data->length * sizeof(REAL8) );
      }
    }
    else
    {
      for ( i = 0; i < num_ifos; i++ )
      {
        injDataREAL4[i] = XLALCreateREAL4TimeSeries( "", &gpsStartTime, 0, 1./sampleRate,
            &lalADCCountUnit, sampleRate * (gpsEndSec - gpsStartSec) );
        memset( injDataREAL4[i]->data->data, 0.0, injDataREAL4[i]->data->length * sizeof(REAL4) );

        if (addNoise) {
          LAL_CALL( add_colored_noise( &status, injDataREAL4[i], i, randParams, dynRange, strainLowPassFreq), &status);
        }
      }
    }

    /* setup a unity response frequency series */
    if (strncmp(injectionType, "approximant", strlen(injectionType) + 1) == 0)
    {
      if (doingREAL8)
      {
        fprintf( stderr, "ERROR: approximant injections are only available for single-precision" );
        exit( 1 );
      }

      if (vrbflg)
        fprintf(stdout, "generating unity response...\n");

      response = XLALCreateCOMPLEX8FrequencySeries("", &gpsStartTime, 0,
          1, &strainPerCount, (sampleRate * (gpsEndSec - gpsStartSec))/2 + 1);
      for ( i = 0; i < (INT4)response->data->length; i++)
      {
        response->data->data[i] = 1.0;
      }
    }

    if (vrbflg)
      fprintf(stdout, "Generating injection for:");

    /* loop over ifos */
    for ( i = 0; i < num_ifos; i++ )
    {
      /* get ifo */
      if ( ifosFlag )
      {
        ifo = (CHAR *) calloc( LIGOMETA_IFO_MAX, sizeof(CHAR));
        XLALReturnIFO( ifo, i );
      }

      if (vrbflg)
      {
        fprintf(stdout, " %s", ifo);
        fflush(stdout);
      }

      /* set the channel name */
      snprintf(channel, LALNameLength, "%s:STRAIN", ifo);
      if (doingREAL8)
        strncpy(injDataREAL8[i]->name, channel, LALNameLength);
      else
        strncpy(injDataREAL4[i]->name, channel, LALNameLength);

      if (strncmp(injectionType, "approximant", strlen(injectionType) + 1) == 0)
      {
        if (doingREAL8)
        {
          fprintf( stderr, "ERROR: approximant injections are only available for single-precision" );
          exit( 1 );
        }
        
        /* inject the specified waveforms */
        LAL_CALL( LALFindChirpInjectSignals( &status, injDataREAL4[i], injections, response), &status);

        /* reset the channel name to IFO:STRAIN as LALFindChirpInjectSignals()
         * messes with it */
        strncpy(injDataREAL4[i]->name, channel, LALNameLength);
      }
      else
      {
        /* inject the numerical waveforms */
        if (doingREAL8)
        {
          LAL_CALL( InjectNumRelWaveformsUsingPSDREAL8 ( &status, injDataREAL8[i], injections, ifo,
                freqLowCutoff, snrLow, snrHigh,
                ligoPsd, ligoSnrLowFreq,
                virgoPsd, virgoSnrLowFreq,
                fnameOutXML), &status);
        }
        else
        {
          LAL_CALL( InjectNumRelWaveforms ( &status, injDataREAL4[i], injections, ifo,
                dynRange, freqLowCutoff, snrLow, snrHigh,
                fnameOutXML), &status);
        }
      }

      /* set strain as unit */
      if (doingREAL8)
        injDataREAL8[i]->sampleUnits = lalStrainUnit;
      else
        injDataREAL4[i]->sampleUnits = lalStrainUnit;

    } /* loop over ifos */

    if (vrbflg)
      fprintf(stdout, "\n");

    /* output frame */
    if ( ifosFlag )
    {
      if ( doingREAL8 )
         output_multi_channel_frame_real8( num_ifos, gpsStartSec, gpsEndSec, injDataREAL8, frameType, outDir );
      else
         output_multi_channel_frame( num_ifos, gpsStartSec, gpsEndSec, injDataREAL4, frameType, outDir );
    }
    else {
      if ( doingREAL8 )
        output_frame_real8( ifo, gpsStartSec, gpsEndSec, injDataREAL8[0], frameType, outDir );
      else
        output_frame( ifo, gpsStartSec, gpsEndSec, injDataREAL4[0], frameType, outDir );
    }
  }

  /* write mdc log */
  if ( mdcFlag )
    write_mdc_log_file(mdcFileName, injections, gpsStartSec, setName);

  /* free memory */
  if ( injectionFile )
    free( injectionFile );

  if ( randParams )
    LAL_CALL( LALDestroyRandomParams( &status, &randParams ), &status );

  if ( ifo )
    free( ifo );

  for ( i = 0; i < num_ifos; i++ )
    if (doingREAL8)
      XLALDestroyREAL8TimeSeries(injDataREAL8[i]);
    else
      XLALDestroyREAL4TimeSeries(injDataREAL4[i]);

  if (response != NULL)
    XLALDestroyCOMPLEX8FrequencySeries(response);

  while ( injections )
  {
    SimInspiralTable *thisInj = NULL;
    thisInj = injections;
    injections = injections->next;
    LALFree( thisInj );
  }

  if (fnameOutXML) {
    free(fnameOutXML);
  }

  if (ligoPsd) {
    XLALDestroyREAL8FrequencySeries( ligoPsd );
  }

  if (virgoPsd) {
    XLALDestroyREAL8FrequencySeries( virgoPsd );
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
      "  [--help]                          display this message\n"\
      "  [--verbose]                       print progress information\n"\
      "  [--version]                       print version information and exit\n"\
      "  --injection-type      type        set injection type ('approximant' or 'NR')\n"\
      "  --injection-file      inj_file    read inj details from xml sim-insp inj_file\n"\
      "  --inject-window       time        Buffer time in which to generate injections. This covers injections which do not merge within the GPS times given but may still overlap the data segment.\n"\
      "  --ifo                 ifo         IFO for which to generate injections\n"\
      "  --all-ifos                        create injections for all IFOs\n"\
      "  --gps-start-time      start       start time of output file\n"\
      "  --gps-end-time        end         end time of output file\n"\
      "  --sample-rate         rate        the sample rate used to generate injections\n"\
      "  --write-mdc-log                   write an MDC log file\n"\
      "  --frame-type          FR_TYPE     set the name of the output frame\n"\
      "  --set-name            set_name    set the injection set name\n"\
      "  --mdc-log             mdc_log     name of file for MDC log file\n"\
      "  --write-frame                     write h(t) waveform to a frame file\n"\
      "  --simulate-noise                  add simulated colored Gaussian noise\n"\
      "  --freq-low-cutoff     freq        lower cutoff frequency for injections\n"\
      "  --strain-lowpass-freq freq        lowpass frequency when noise is produced\n"\
      "  --snr-low             snr_lo      lower cutoff on snr\n"\
      "  --snr-high            snr_hi      upper cutoff on snr\n"\
      "  --out-xml-file        output xml  output file with list of injections performed\n"\
      "  --fr-out-dir          dir         directory to output frames to\n"\
      "  --double-precision                read REAL8 NR files, produce REAL8 injections\n"\
      "\n", program );
}


/* function to output h(t) waveform in a frame file */
static void output_frame(CHAR *ifo,
    INT4 gpsStart,
    INT4 gpsEnd,
    REAL4TimeSeries *injData,
    CHAR *frameType,
    CHAR *outDir)
{
  CHAR fname[FILENAME_MAX];
  INT4 duration;
  INT4 detectorFlags;
  LALFrameH *frame;
  CHAR creator[HISTORY_COMMENT];
  CHAR channel[LALNameLength];

  /* get frame filename */
  duration = gpsEnd - gpsStart;
  if (outDir)
    snprintf( fname, FILENAME_MAX, "%s/GHLTV-%s-%d-%d.gwf", outDir, frameType, gpsStart, duration );
  else
    snprintf( fname, FILENAME_MAX, "GHLTV-%s-%d-%d.gwf", frameType, gpsStart, duration );

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

  /* set the channel name */
  snprintf(channel, LALNameLength, "%s:STRAIN", ifo);
  strncpy(injData->name, channel, LALNameLength);

  /* define frame */
  frame = XLALFrameNew( &injData->epoch, duration, "LIGO", 0, 1,
      detectorFlags );

  /* set creator metadata */
  /** \deprecated FIXME: the following code uses obsolete CVS ID tags.
   *  It should be modified to use git version information. */
  snprintf(creator, HISTORY_COMMENT, "creator:$Id$");
  XLALFrameAddFrHistory(frame, "creator", creator);

  /* add channel to frame */
  XLALFrameAddREAL4TimeSeriesSimData( frame, injData );

  if ( vrbflg )
    fprintf( stdout, "Writing injection to frame: '%s'\n", fname );

  /* write frame */
  if (XLALFrameWrite( frame, fname) != 0)
  {
    fprintf( stderr, "ERROR: Cannot save frame file: '%s'\n", fname );
    exit( 1 );
  }

  /* clear frame */
  XLALFrameFree( frame );

  return;
}


/* function to output h(t) waveform in a frame file */
static void output_frame_real8(CHAR *ifo,
    INT4 gpsStart,
    INT4 gpsEnd,
    REAL8TimeSeries *injData,
    CHAR *frameType,
    CHAR *outDir)
{
  CHAR fname[FILENAME_MAX];
  INT4 duration;
  INT4 detectorFlags;
  LALFrameH *frame;
  CHAR creator[HISTORY_COMMENT];
  CHAR channel[LALNameLength];

  /* get frame filename */
  duration = gpsEnd - gpsStart;
  if (outDir)
    snprintf( fname, FILENAME_MAX, "%s/GHLTV-%s-%d-%d.gwf", outDir, frameType, gpsStart, duration );
  else
    snprintf( fname, FILENAME_MAX, "GHLTV-%s-%d-%d.gwf", frameType, gpsStart, duration );

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

  /* set the channel name */
  snprintf(channel, LALNameLength, "%s:STRAIN", ifo);
  strncpy(injData->name, channel, LALNameLength);

  /* define frame */
  frame = XLALFrameNew( &injData->epoch, duration, "LIGO", 0, 1,
      detectorFlags );

  /* set creator metadata */
  /** \deprecated FIXME: the following code uses obsolete CVS ID tags.
   *  It should be modified to use git version information. */
  snprintf(creator, HISTORY_COMMENT, "creator:$Id$");
  XLALFrameAddFrHistory(frame, "creator", creator);

  /* add channel to frame */
  XLALFrameAddREAL8TimeSeriesSimData( frame, injData );

  if ( vrbflg )
    fprintf( stdout, "Writing injection to frame: '%s'\n", fname );

  /* write frame */
  if (XLALFrameWrite( frame, fname) != 0)
  {
    fprintf( stderr, "ERROR: Cannot save frame file: '%s'\n", fname );
    exit( 1 );
  }

  /* clear frame */
  XLALFrameFree( frame );

  return;
}

/* write injections for all ifos into a single frame */
static void output_multi_channel_frame(INT4 num_ifos,
    INT4 gpsStart,
    INT4 gpsEnd,
    REAL4TimeSeries *injData[LAL_NUM_IFO],
    CHAR *frameType,
    CHAR *outDir)
{
  CHAR fname[FILENAME_MAX];
  INT4 duration;
  INT4 detectorFlags;
  LALFrameH *frame;
  INT4 i;
  CHAR creator[HISTORY_COMMENT];
  CHAR *ifo = NULL;

  /* get frame filename */
  duration = gpsEnd - gpsStart;
  if (outDir)
    snprintf( fname, FILENAME_MAX, "%s/GHLTV-%s-%d-%d.gwf", outDir, frameType, gpsStart, duration );
  else
    snprintf( fname, FILENAME_MAX, "GHLTV-%s-%d-%d.gwf", frameType, gpsStart, duration );

  /* set detector flags */
  detectorFlags = LAL_GEO_600_DETECTOR_BIT | LAL_LHO_4K_DETECTOR_BIT |
    LAL_LHO_2K_DETECTOR_BIT | LAL_LLO_4K_DETECTOR_BIT |
    LAL_TAMA_300_DETECTOR_BIT | LAL_VIRGO_DETECTOR_BIT;

  /* define frame */
  frame = XLALFrameNew( &(injData[0])->epoch, duration, "LIGO", 0, 1,
      detectorFlags );

  /* set creator metadata */
  /** \deprecated FIXME: the following code uses obsolete CVS ID tags.
   *  It should be modified to use git version information. */
  snprintf(creator, HISTORY_COMMENT, "creator:$Id$");
  XLALFrameAddFrHistory(frame, "creator", creator);

  /* add channels to frame */
  for( i = 0; i < num_ifos; i++ )
  {
    ifo = (CHAR*)calloc(LIGOMETA_IFO_MAX, sizeof(CHAR));
    XLALReturnIFO(ifo, i);
    printf("adding %s channel to frame\n", ifo);
    XLALFrameAddREAL4TimeSeriesSimData( frame, injData[i] );
  }

  if (vrbflg)
    fprintf( stdout, "Writing injections to frame: '%s'\n", fname );

  /* write frame */
  if ( XLALFrameWrite( frame, fname ) != 0 )
  {
    fprintf( stderr, "ERROR: Cannot save frame file: '%s'\n", fname );
    exit( 1 );
  }

  /* clear frame */
  XLALFrameFree( frame );

  return;
}


/* write injections for all ifos into a single frame */
static void output_multi_channel_frame_real8(INT4 num_ifos,
    INT4 gpsStart,
    INT4 gpsEnd,
    REAL8TimeSeries *injData[LAL_NUM_IFO],
    CHAR *frameType,
    CHAR *outDir)
{
  CHAR fname[FILENAME_MAX];
  INT4 duration;
  INT4 detectorFlags;
  LALFrameH *frame;
  INT4 i;
  CHAR creator[HISTORY_COMMENT];
  CHAR *ifo = NULL;

  /* get frame filename */
  duration = gpsEnd - gpsStart;
  if (outDir)
    snprintf( fname, FILENAME_MAX, "%s/GHLTV-%s-%d-%d.gwf", outDir, frameType, gpsStart, duration );
  else
    snprintf( fname, FILENAME_MAX, "GHLTV-%s-%d-%d.gwf", frameType, gpsStart, duration );

  /* set detector flags */
  detectorFlags = LAL_GEO_600_DETECTOR_BIT | LAL_LHO_4K_DETECTOR_BIT |
    LAL_LHO_2K_DETECTOR_BIT | LAL_LLO_4K_DETECTOR_BIT |
    LAL_TAMA_300_DETECTOR_BIT | LAL_VIRGO_DETECTOR_BIT;

  /* define frame */
  frame = XLALFrameNew( &(injData[0])->epoch, duration, "LIGO", 0, 1,
      detectorFlags );

  /* set creator metadata */
  /** \deprecated FIXME: the following code uses obsolete CVS ID tags.
   *  It should be modified to use git version information. */
  snprintf(creator, HISTORY_COMMENT, "creator:$Id$");
  XLALFrameAddFrHistory(frame, "creator", creator);

  /* add channels to frame */
  for( i = 0; i < num_ifos; i++ )
  {
    ifo = (CHAR*)calloc(LIGOMETA_IFO_MAX, sizeof(CHAR));
    XLALReturnIFO(ifo, i);
    printf("adding %s channel to frame\n", ifo);
    XLALFrameAddREAL8TimeSeriesSimData( frame, injData[i] );
  }

  if (vrbflg)
    fprintf( stdout, "Writing injections to frame: '%s'\n", fname );

  /* write frame */
  if ( XLALFrameWrite( frame, fname ) != 0 )
  {
    fprintf( stderr, "ERROR: Cannot save frame file: '%s'\n", fname );
    exit( 1 );
  }

  /* clear frame */
  XLALFrameFree( frame );

  return;
}
/* function to write a Burst MDC log file */
static void write_mdc_log_file(CHAR *filename,
    SimInspiralTable *injections,
    INT4 gps_start,
    CHAR *set_name)
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
    fprintf(stdout, "Writing MDC log file: '%s'\n", filename);

  /* loop over injections */
  for (thisInj = injections; thisInj; thisInj = thisInj->next)
  {
    /* declare variables */
    REAL8 gmst;
    REAL8 longitude;

    /* GravEn_SimID */
    if (strncmp(thisInj->waveform, "NumRel", LIGOMETA_WAVEFORM_MAX) == 0)
      fprintf(output, "%s ", thisInj->numrel_data);
    else
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
    gmst = fmod(XLALGreenwichMeanSiderealTime(&thisInj->geocent_end_time), LAL_TWOPI);
    longitude = fmod(thisInj->longitude - gmst, LAL_TWOPI);
    if (longitude < 0)
      longitude += LAL_TWOPI;
    fprintf(output, "%g ", longitude);
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

static void add_colored_noise(LALStatus *status,
    REAL4TimeSeries *chan,
    INT4 ifoNumber,
    RandomParams *randParams,
    REAL8 dynRange,
    REAL8 strainHighPassFreq)
{

  UINT4 k;
  COMPLEX8Sequence *ntilde         = NULL;
  REAL4Sequence    *ntilde_re      = NULL;
  REAL4Sequence    *ntilde_im      = NULL;
  REAL8Sequence    *spectrum       = NULL;
  REAL4FFTPlan     *invPlan        = NULL;
  INT4              length         = chan->data->length;
  REAL8             deltaT         = chan->deltaT;
  REAL8             deltaF         = 1.0 / (deltaT * (REAL8) length);
  REAL8             tObs           = length * deltaT;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  /* Generate white Gaussian noise with unit variance */
  TRY( LALSCreateVector( status->statusPtr, &ntilde_re, length / 2 + 1 ), status );

  TRY( LALNormalDeviates( status->statusPtr, ntilde_re, randParams ), status );
  TRY( LALSCreateVector( status->statusPtr, &ntilde_im, length / 2 + 1 ), status );
  TRY( LALNormalDeviates( status->statusPtr, ntilde_im, randParams ), status );

  /* create storage for the frequency domain noise and psd*/
  TRY( LALCCreateVector( status->statusPtr, &ntilde, length / 2 + 1 ), status );
  TRY( LALDCreateVector( status->statusPtr, &spectrum, length / 2 + 1 ), status );

  get_spectrum( spectrum, ifoNumber, deltaF, strainHighPassFreq, dynRange);

  /* Color white noise with given psd */
  for ( k=0; k < ntilde->length; k++ )
  {
    ntilde->data[k] = crectf( ntilde_re->data[k] * sqrt( 0.25 * tObs * spectrum->data[k] ),
                              ntilde_im->data[k] * sqrt( 0.25 * tObs * spectrum->data[k] ) );
  }
  /* setting d.c. and Nyquist to zero */
  ntilde->data[0] = crealf(ntilde->data[0]);
  ntilde->data[length / 2] = crealf(ntilde->data[length / 2]);


  /*   fp = fopen("ntilde.dat", "w"); */
  /*   for ( k = 0; k < length/2 + 1; k++) */
  /*     fprintf(fp, "%e   %e   %e   %e    %e   %e\n", k*deltaF,
         ntilde_re->data[k],
         ntilde_im->data[k], ntilde->data[k].re, ntilde->data[k].im,
         spectrum->data[k]); */
  /*   fclose(fp); */


  /* Fourier transform back in the time domain */
  TRY( LALCreateReverseRealFFTPlan( status->statusPtr, &invPlan, length, 0 ), status );
  TRY( LALReverseRealFFT( status->statusPtr, chan->data, ntilde, invPlan ), status);

  /* normalize the noise */
  for ( k = 0; k < (UINT4)length; ++k )
  {
    chan->data->data[k] /= (REAL4) length ;
  }

  if ( vrbflg ) fprintf( stdout, "done\n" );

  TRY( LALSDestroyVector( status->statusPtr, &ntilde_re ), status);
  TRY( LALSDestroyVector( status->statusPtr, &ntilde_im ), status);
  TRY( LALCDestroyVector( status->statusPtr, &ntilde ), status);
  TRY( LALDDestroyVector( status->statusPtr, &spectrum), status);
  TRY( LALDestroyRealFFTPlan( status->statusPtr, &invPlan), status);

  DETATCHSTATUSPTR (status);
  RETURN (status);

}


static int get_spectrum(REAL8Sequence *spectrum,
    InterferometerNumber ifoNumber,
    REAL8 deltaF,
    REAL8 strainHighPassFreq,
    REAL8 dynRange)
{
  UINT4 k;
  INT4 kmin = strainHighPassFreq / deltaF > 1 ? strainHighPassFreq / deltaF : 1 ;

  if ( (ifoNumber == LAL_IFO_H1) || (ifoNumber == LAL_IFO_L1) )
  {
    /* set the spectrum to the Initial LIGO design noise curve */
    REAL8 psd_value;
    LALLIGOIPsd( NULL, &psd_value, strainHighPassFreq );
    for ( k = 0; k < (UINT4)kmin ; ++k )
    {
      spectrum->data[k] = 9.0e-46 * psd_value * dynRange * dynRange;
    }
    for ( k = kmin; k < spectrum->length ; ++k )
    {
      REAL8 psd_freq = (REAL8) k * deltaF;
      LALLIGOIPsd( NULL, &psd_value, psd_freq );
      spectrum->data[k] = 9.0e-46 * psd_value * dynRange * dynRange;
    }
  }
  else if ( ifoNumber == LAL_IFO_V1)
  {
    /* set the spectrum to the Advanced LIGO design noise curve */
    REAL8 psd_value;
    LALVIRGOPsd( NULL, &psd_value, strainHighPassFreq );
    for ( k = 0; k < (UINT4)kmin ; ++k )
    {
      spectrum->data[k] = 9.0e-46 * psd_value * dynRange * dynRange;
    }
    for ( k = kmin; k < spectrum->length ; ++k )
    {
      REAL8 psd_freq = (REAL8) k * deltaF;
      LALAdvLIGOPsd( NULL, &psd_value, psd_freq );
      spectrum->data[k] = 9.0e-46 * psd_value * dynRange * dynRange;
    }
  }
  else if ( ifoNumber == LAL_IFO_T1 )
  {
    /* set the spectrum to the Advanced LIGO design noise curve */
    REAL8 psd_value;
    LALTAMAPsd( NULL, &psd_value, strainHighPassFreq );
    for ( k = 0; k < (UINT4)kmin ; ++k )
    {
      spectrum->data[k] = 9.0e-46 * psd_value * dynRange * dynRange;
    }
    for ( k = kmin; k < spectrum->length ; ++k )
    {
      REAL8 psd_freq = (REAL8) k * deltaF;
      LALAdvLIGOPsd( NULL, &psd_value, psd_freq );
      spectrum->data[k] = 9.0e-46 * psd_value * dynRange * dynRange;
    }
  }
  else if ( ifoNumber == LAL_IFO_G1 )
  {
    /* set the spectrum to the Advanced LIGO design noise curve */
    REAL8 psd_value;
    LALGEOPsd( NULL, &psd_value, strainHighPassFreq );
    for ( k = 0; k < (UINT4)kmin ; ++k )
    {
      spectrum->data[k] = 9.0e-46 * psd_value * dynRange * dynRange;
    }
    for ( k = kmin; k < spectrum->length ; ++k )
    {
      REAL8 psd_freq = (REAL8) k * deltaF;
      LALAdvLIGOPsd( NULL, &psd_value, psd_freq );
      spectrum->data[k] = 9.0e-46 * psd_value * dynRange * dynRange;
    }
  }
  else if ( ifoNumber == LAL_IFO_H2 )
  {
    /* set the spectrum to the Advanced LIGO design noise curve */
    REAL8 psd_value;
    LALLIGOIPsd( NULL, &psd_value, strainHighPassFreq );
    for ( k = 0; k < (UINT4)kmin ; ++k )
    {
      spectrum->data[k] = 4.0 * 9.0e-46 * psd_value * dynRange * dynRange;
    }
    for ( k = kmin; k < spectrum->length ; ++k )
    {
      REAL8 psd_freq = (REAL8) k * deltaF;
      LALAdvLIGOPsd( NULL, &psd_value, psd_freq );
      spectrum->data[k] = 2.0 * 9.0e-46 * psd_value * dynRange * dynRange;
    }
  }
  else
    return -1;

  return 0;

}
