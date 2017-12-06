/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Lisa M. Goggin, Matt Pitkin
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>

#include <lal/LALStdio.h>
#include <lal/LIGOMetadataRingdownUtils.h>
#include "lalapps.h"
#include "errutil.h"
#include "gpstime.h"
#include "ring.h"
#include "injsgnl.h"

extern int vrbflg;
static int ring_usage( const char *program );
static int ring_default_params( struct ring_params *params );


/* parse command line arguments using getopt_long to get ring params */
int ring_parse_options( struct ring_params *params, int argc, char **argv )
{
  static struct ring_params localparams;
  struct option long_options[] =
  {
    { "verbose", no_argument, &vrbflg, 1 },
    { "white-spectrum",          no_argument, &localparams.whiteSpectrum, 1 },
    { "bank-only",               no_argument, &localparams.bankOnly, 1 },
    { "write-raw-data",          no_argument, &localparams.writeRawData, 1 },
    { "write-data",              no_argument, &localparams.writeProcessedData, 1 },
    { "write-response",          no_argument, &localparams.writeResponse, 1 },
    { "write-spectrum",          no_argument, &localparams.writeSpectrum, 1 },
    { "write-inv-spectrum",      no_argument, &localparams.writeInvSpectrum, 1 },
    { "write-segment",           no_argument, &localparams.writeSegment, 1 },
    { "write-template-time-series", no_argument, &localparams.writeTemplateTimeSeries, 1 },
    { "write-template-fft",      no_argument, &localparams.writeTemplateFFT, 1 },
    { "write-filter-output",     no_argument, &localparams.writeFilterOutput, 1 },
    { "write-compress",          no_argument, &localparams.outCompress, 1 },
    { "help",                    no_argument,       0, 'h' },
    { "version",                 no_argument,       0, 'V' },
    { "gps-start-time",          required_argument, 0, 'a' },
    { "gps-start-time-ns",       required_argument, 0, 'A' },
    { "gps-end-time",            required_argument, 0, 'b' },
    { "gps-end-time-ns",         required_argument, 0, 'B' },
    { "channel-name",            required_argument, 0, 'c' },
    { "calibration-cache",       required_argument, 0, 'C' },
    { "frame-cache",             required_argument, 0, 'D' },
    { "cutoff-frequency",        required_argument, 0, 'e' },
    { "highpass-frequency",      required_argument, 0, 'E' },
    { "bank-min-frequency",      required_argument, 0, 'f' },
    { "bank-max-frequency",      required_argument, 0, 'F' },
    { "data-type",               required_argument, 0, 'G' },
    { "spectrum-type",           required_argument, 0, 'L' },
    { "injection-type",          required_argument, 0, 'J' },
    { "injection-file",          required_argument, 0, 'i' },
    { "inject-mdc-frame",        required_argument, 0, 'I' },
    { "user-tag",                required_argument, 0, 'k' },
    { "ifo-tag",                 required_argument, 0, 'K' },
    { "bank-max-mismatch",       required_argument, 0, 'm' },
    { "maximize-duration",       required_argument, 0, 'M' },
    { "only-segment-numbers",    required_argument, 0, 'n' },
    { "only-template-numbers",   required_argument, 0, 'N' },
    { "output-file",             required_argument, 0, 'o' },
    { "bank-file",               required_argument, 0, 'O' },
    { "bank-template-phase",     required_argument, 0, 'p' },
    { "bank-min-quality",        required_argument, 0, 'q' },
    { "bank-max-quality",        required_argument, 0, 'Q' },
    { "random-seed",             required_argument, 0, 'r' },
    { "dynamic-range-factor",    required_argument, 0, 'R' },
    { "sample-rate",             required_argument, 0, 's' },
    { "segment-duration",        required_argument, 0, 'S' },
    { "threshold",               required_argument, 0, 't' },
    { "inverse-spec-length",     required_argument, 0, 'T' },
    { "trig-start-time",         required_argument, 0, 'u' },
    { "trig-end-time",           required_argument, 0, 'U' },
    { "block-duration",          required_argument, 0, 'w' },
    { "pad-data",                required_argument, 0, 'W' },
    { 0, 0, 0, 0 }
  };
  char args[] = "a:A:b:B:c:C:d:D:e:E:f:F:G:hi:I:J:L:m:o:O:p:q:Q:r:R:s:S:t:T:u:U:V:w:W";
  char *program = argv[0];

  /* set default values for parameters before parsing arguments */
  ring_default_params( &localparams );

  while ( 1 )
  {
    int option_index = 0;
    int c;

    c = getopt_long_only( argc, argv, args, long_options, &option_index );
    if ( c == -1 ) /* end of options */
      break;

    switch ( c )
    {
      case 0: /* if option set a flag, nothing else to do */
        if ( long_options[option_index].flag )
          break;
        else
          error( "error parsing option %s with argument %s\n",
              long_options[option_index].name, optarg );
      case 'a': /* gps-start-time */
        localparams.startTime.gpsSeconds = atol( optarg );
        break;
      case 'A': /* gps-start-time-ns */
        localparams.startTime.gpsNanoSeconds = atol( optarg );
        break;
      case 'b': /* gps-end-time */
        localparams.endTime.gpsSeconds = atol( optarg );
        break;
      case 'B': /* gps-end-time-ns */
        localparams.endTime.gpsNanoSeconds = atol( optarg );
        break;
      case 'c': /* channel-name */
        localparams.channel = optarg;
        break;
      case 'C': /* calibration-cache */
        localparams.calibCache = optarg;
        break;
      case 'D': /* frame-cache */
        localparams.dataCache = optarg;
        break;
      case 'e': /* cutoff-frequency */
        localparams.lowCutoffFrequency = atof( optarg );
        break;
      case 'E': /* highpass-frequency */
        localparams.highpassFrequency = atof( optarg );
        break;
      case 'f': /* bank min frequency */
        localparams.bankParams.minFrequency = atof( optarg );
        break;
      case 'F': /* bank max frequency */
        localparams.bankParams.maxFrequency = atof( optarg );
        break;
     case 'h': /* help */
        ring_usage( program );
        exit( 0 );
      case 'i': /* injection-file */
        localparams.injectFile = optarg;
        break;
      case 'G': /* data type */
        if( ! strcmp( "sim", optarg ) )
        {
          localparams.dataType = 0;
        }
        else if( ! strcmp( "zero", optarg ) )
        {
          localparams.dataType = 1;
        }
        else if( ! strcmp( "uncal", optarg ) )
        {
          localparams.dataType = 2;
        }
        else if( ! strcmp( "ht_real4", optarg ) )
        {
          localparams.dataType = 3;
        }
        else if( ! strcmp( "ht_real8", optarg ) )
        {
          localparams.dataType = 4;
        }
        else
        {
          localparams.dataType = -1;
          fprintf( stderr, "invalid --data_type:\n"
              "(must be sim, zero, uncal, ht_real4 or ht_real8)\n" );
          exit( 1 );
        }
        break;
      case 'L': /* spectrum type */
        if( ! strcmp( "median", optarg ) )
        {
          localparams.spectrumType = 0;
        }
        else if( ! strcmp( "median_mean", optarg ) )
        {
          localparams.spectrumType = 1;
        }
        else
        {
          localparams.spectrumType = -1;
          fprintf( stderr, "invalid --spectrum_type:\n"
              "(must be median or median_mean)\n" );
          exit( 1 );
        }
        break;
      case 'J': /* injection type */
        if( ! strcmp( "RINGDOWN", optarg ) )
        { 
          localparams.injectType = 0;
        }
        else if( ! strcmp( "IMR", optarg ) )
        {
          localparams.injectType = 1;
        }
        else if( ! strcmp( "IMR_RINGDOWN", optarg ) )
        {
          localparams.injectType = 2;
        }
        else if( ! strcmp( "EOBNR", optarg ) )
        {
          localparams.injectType = 3;
        }
        else if( ! strcmp( "PHENOM", optarg ) )
        {
          localparams.injectType = 4;
        }
        else
        {
          localparams.injectType = -1;
          fprintf( stderr, "invalid --injection_type:\n"
              "(must be RINGDOWN, IMR, IMR_RINGDOWN, EOBNR or PHENOM)\n" );
          exit( 1 );
        }
        break;
      case 'I': /* inject-mdc-frame */
        error( "currently unsupported option: --inject-mdc-frame\n" );
        break;
      case 'k': /* user-tag */
        strncpy( localparams.userTag, optarg, sizeof( localparams.userTag ) - 1 );
        break;
      case 'K': /* ifo-tag */
        strncpy( localparams.ifoTag, optarg, sizeof( localparams.ifoTag ) - 1 );
        break;
      case 'm': /* bank max mismatch */
        localparams.bankParams.maxMismatch = atof( optarg );
        break;
      case 'M': /* maximize duration */
        localparams.maximizeEventDuration = atof( optarg );
        break;
      case 'n': /* only-segment-numbers */
        localparams.segmentsToDoList = optarg;
        break;
      case 'N': /* only-template-number */
        localparams.templatesToDoList = optarg;
        break;
      case 'o': /* output-file */
        strncpy( localparams.outputFile, optarg, sizeof( localparams.outputFile ) - 1 );
        break;
      case 'O': /* bank-file */
        strncpy( localparams.bankFile, optarg, sizeof( localparams.bankFile ) - 1 );
        break;
      case 'p': /* bank template phase */
        localparams.bankParams.templatePhase = atof( optarg );
        break;
      case 'q': /* bank min quality */
        localparams.bankParams.minQuality = atof( optarg );
        break;
      case 'Q': /* bank max quality */
        localparams.bankParams.maxQuality = atof( optarg );
        break;
      case 'r': /* random seed */
        localparams.randomSeed = atoi( optarg );
        break;
      case 'R': /* dynamic range factor */
        localparams.dynRangeFac = atof( optarg );
        break;
      case 's': /* sample rate */
        localparams.sampleRate = atof( optarg );
        break;
      case 'S': /* segment-duration */
        localparams.segmentDuration = atof( optarg );
        break;
      case 't': /* threshold */
        localparams.threshold = atof( optarg );
        break;
      case 'T': /* inverse-spec-length */
        localparams.invSpecLen = atof( optarg );
        break;
      case 'u': /* trig-start-time */
        localparams.trigStartTimeNS = (INT8) atol( optarg ) * LAL_INT8_C(1000000000);
        break;
      case 'U': /* trig-end-time */
        localparams.trigEndTimeNS = (INT8) atol( optarg ) * LAL_INT8_C(1000000000);
        break;
      case 'w': /* block-duration */
        localparams.duration = atof( optarg );
        break;
      case 'W': /* pad-data */
        localparams.padData = atof( optarg );
        break;
      case 'V': /* version */
        XLALOutputVersionString(stderr, 0);
        exit( 0 );
      case '?':
        error( "unknown error while parsing options\n" );
      default:
        error( "unknown error while parsing options\n" );
    }
  }

  if ( optind < argc )
  {
    fprintf( stderr, "extraneous command line arguments:\n" );
    while ( optind < argc )
      fprintf( stderr, "%s\n", argv[optind++] );
    exit( 1 );
  }

  *params = localparams;

  return 0;
}

/* sets default values for parameters */
static int ring_default_params( struct ring_params *params )
{
  /* overall, default values are zero */
  memset( params, 0, sizeof( *params ) );

  /* dynamic range factor must be greater than zero */
  params->dynRangeFac = 1.0;

  /* generate a template at 1 Mpc with an epsilon of 0.01 */
  params->bankParams.templateDistance = 1.0;
  params->bankParams.templateEpsilon  = 0.01;

  /* negative value means use the "default" values */
  params->highpassFrequency     = -1.0; /* use low-frequency cutoff */
  params->maximizeEventDuration = -1.0; /* use filter duration */

  /* segments and templates to do: all of them */
  params->segmentsToDoList  = "^-$";
  params->templatesToDoList = "^-$";

  /* flags specifying what to do: default is to do everything */
  params->getBank     = 1;
  params->getData     = 1;
  params->getResponse = 1;
  params->getSpectrum = 1;
  params->doFilter    = 1;
  
  params->injectType  = -1;
  params->spectrumType = -1;

  return 0;
}


/* macro for testing validity of a condition that prints an error if invalid */
#define sanity_check( condition ) \
  ( condition ? 0 : ( fputs( #condition " not satisfied\n", stderr ), error( "sanity check failed\n" ) ) ) 

/* check sanity of parameters and sets appropriate values of unset parameters */
int ring_params_sanity_check( struct ring_params *params )
{
  UINT4 recordLength = 0;
  UINT4 segmentLength = 0;
  UINT4 segmentStride = 0;
  UINT4 truncateLength = 0;
  INT8 startTime;
  INT8 endTime;
  int validChannelIFO;

  if ( params->dataType == LALRINGDOWN_DATATYPE_HT_REAL4 || params->dataType == LALRINGDOWN_DATATYPE_HT_REAL8)
      params->getResponse = 0;

  if ( params->bankOnly )
  {
    params->getData     = 0;
    params->getResponse = 0;
    params->getSpectrum = 0;
    params->doFilter    = 0;
  }
  
  if ( params->dataType == LALRINGDOWN_DATATYPE_SIM )
    sanity_check( params->randomSeed );

  if ( params->dataType == LALRINGDOWN_DATATYPE_UNCAL )
    sanity_check( params->calibCache );

  if ( params->doFilter ) /* need data, bank, and spectrum */
    sanity_check( params->getData && params->getBank && params->getSpectrum );

  /* parameters required to get data */
  if ( params->getData )
  {
    /* checks on data duration */
    startTime = epoch_to_ns( &params->startTime );
    endTime   = epoch_to_ns( &params->endTime );
    sanity_check( startTime > 0 );
    sanity_check( endTime > startTime );
    sanity_check( params->duration > 0 );
    sanity_check( 1e9*params->duration == ((endTime - startTime)) );

    /* checks on size of data record */
    sanity_check( params->sampleRate > 0 );
    recordLength = params->duration * params->sampleRate;
    sanity_check( recordLength > 0 );

    sanity_check( params->channel );
    sanity_check( params->dataType == LALRINGDOWN_DATATYPE_SIM || params->dataCache );

    /* record ifo name */
    validChannelIFO = sscanf( params->channel, "%2[A-Z1-9]", params->ifoName );
    sanity_check( validChannelIFO );

    /* check the spectrum type is specified */
    sanity_check( params->spectrumType >= 0.0 );

    /* check that injection type is specified if an injection file is given */

    if ( params->injectFile ) /* geo data parameters */
    {
      sanity_check( params->injectType >= 0.0 );
    }
    
    /* will need response to do injections unless strain data */
    sanity_check( params->injectFile == NULL || (params->dataType == LALRINGDOWN_DATATYPE_HT_REAL4 
      || params->dataType == LALRINGDOWN_DATATYPE_HT_REAL8) || params->getResponse );
  }

  /* parameters required to get spectrum */
  if ( params->getSpectrum )
  {
    /* checks on size of data segments and stride */
    sanity_check( params->segmentDuration > 0 );
    segmentLength = floor(params->segmentDuration * params->sampleRate + 0.5);
    sanity_check( recordLength / segmentLength > 0 );
    params->strideDuration = 0.5 * params->segmentDuration;
    segmentStride = floor(params->strideDuration * params->sampleRate + 0.5);
    sanity_check( segmentStride > 0 );
    params->truncateDuration = 0.25 * params->strideDuration;
    truncateLength = floor(params->truncateDuration * params->sampleRate + 0.5);
    sanity_check( truncateLength > 0 );
    /* record length, segment length and stride need to be commensurate */
    sanity_check( !( (recordLength - segmentLength) % segmentStride ) );
    params->numOverlapSegments = 1 + (recordLength - segmentLength)/segmentStride;

    if ( params->spectrumType > 0 )
      sanity_check( ! (params->numOverlapSegments % 2) ); /* required to be even for median-mean method */

    /* checks on data input information */
    sanity_check( params->channel );
    sanity_check( params->dynRangeFac > 0.0 );
  }

  /* parameters required to get response */
  if ( params->getResponse )
  {
    sanity_check( params->calibCache );
  }

  /* parameters required to do filtering */
  if ( params->doFilter )
  {
    /* checks on low-cutoff and highpass frequencies */
    sanity_check( params->lowCutoffFrequency > 0 );
    sanity_check( params->lowCutoffFrequency < 0.5 * params->sampleRate );
    if ( params->highpassFrequency < 0 )
      params->highpassFrequency = params->lowCutoffFrequency;
    sanity_check( params->lowCutoffFrequency >= params->highpassFrequency );

    /* checks on filter threshold */
    sanity_check( params->threshold > 0.0 );

    /* output file name */
    if ( ! strlen( params->outputFile ) )
    {
      if ( strlen( params->userTag ) && strlen( params->ifoTag ) && params->outCompress )
      {
        snprintf( params->outputFile, sizeof( params->outputFile ),
          "%s-RING_%s_%s-%d-%d.xml.gz", params->ifoName, 
          params->ifoTag, params->userTag, params->startTime.gpsSeconds,
          (int)ceil( params->duration ) );
      }
      else if ( strlen( params->userTag ) && strlen( params->ifoTag )  && !params->outCompress )
      {
        snprintf( params->outputFile, sizeof( params->outputFile ),
          "%s-RING_%s_%s-%d-%d.xml", params->ifoName,
          params->ifoTag, params->userTag, params->startTime.gpsSeconds,
          (int)ceil( params->duration ) );
      }
      else if ( strlen( params->userTag ) && !strlen( params->ifoTag ) && params->outCompress )
      {
        snprintf( params->outputFile, sizeof( params->outputFile ),
          "%s-RING_%s-%d-%d.xml.gz", params->ifoName, 
          params->userTag, params->startTime.gpsSeconds,
          (int)ceil( params->duration ) );
      }
      else if ( strlen( params->userTag ) && !strlen( params->ifoTag ) && !params->outCompress  )
      {
        snprintf( params->outputFile, sizeof( params->outputFile ),
          "%s-RING_%s-%d-%d.xml", params->ifoName,
          params->userTag, params->startTime.gpsSeconds,
          (int)ceil( params->duration ) );
      }
      else if ( !strlen( params->userTag ) && strlen( params->ifoTag ) && params->outCompress )
      {
        snprintf( params->outputFile, sizeof( params->outputFile ),
          "%s-RING_%s-%d-%d.xml.gz", params->ifoName, 
          params->ifoTag, params->startTime.gpsSeconds,
          (int)ceil( params->duration ) );
      }
      else if ( !strlen( params->userTag ) && strlen( params->ifoTag ) && !params->outCompress )
      {
        snprintf( params->outputFile, sizeof( params->outputFile ),
          "%s-RING_%s-%d-%d.xml", params->ifoName,
          params->ifoTag, params->startTime.gpsSeconds,
          (int)ceil( params->duration ) );
      }
      else 
      {
        snprintf( params->outputFile, sizeof( params->outputFile ),
          "%s-RING-%d-%d.xml", params->ifoName, 
          params->startTime.gpsSeconds,
          (int)ceil( params->duration ) );
      }
    }

  }

  /* parameters required to make bank */
  if ( params->getBank )
  {
    /* checks on bank parameters */
    sanity_check( params->bankParams.minFrequency > params->lowCutoffFrequency );
    sanity_check( params->bankParams.maxFrequency > params->bankParams.minFrequency );
    sanity_check( params->bankParams.maxFrequency < 0.5 * params->sampleRate );
    sanity_check( params->bankParams.minQuality >= 2.0 );
    sanity_check( params->bankParams.maxQuality > params->bankParams.minQuality );
    sanity_check( params->bankParams.maxMismatch > 0.0 );
    sanity_check( params->bankParams.maxMismatch < 1.0 );
  }

  return 0;
}


/* prints a help message */
static int ring_usage( const char *program )
{
  fprintf( stderr, "usage: %s options\n", program );
  fprintf( stderr, "\ngeneral options:\n" );
  fprintf( stderr, "--help                     print this message\n" );
  fprintf( stderr, "--version                  print the version of the code\n" );
  fprintf( stderr, "--verbose                  print verbose messages while running\n" );

  fprintf( stderr, "\ndata reading options:\n" );
  fprintf( stderr, "--frame-cache=cachefile    name of the frame cache file\n" );
  fprintf( stderr, "--channel-name             data channel to analyze\n" );
  fprintf( stderr, "--gps-start-time=tstart    GPS start time of data to analyze (sec)\n" );
  fprintf( stderr, "--gps-start-time-ns=tstartns  nanosecond residual of start time\n" );
  fprintf( stderr, "--gps-end-time=tstop       GPS stop time of data to analyze (sec)\n" );
  fprintf( stderr, "--gps-end-time-ns=tstopns  nanosecond residual of stop time\n" );

  fprintf( stderr, "\ndata conditioning options:\n" );
  fprintf( stderr, "--highpass-frequency=fhi   high-pass filter data at frequency fhi (Hz)\n" );
  fprintf( stderr, "--sample-rate=srate        decimate data to be at sample rate srate (Hz)\n" );

  fprintf( stderr, "\nsimulated injection options:\n" );
  fprintf( stderr, "--injection-type       type of injection, must be one of \n \t \t [RINGDOWN, EOBNR, PHENOM, IMR, IMR_RINGDOWN] \n \t and must be accompanied by the appropriate injection file\n" );
  fprintf( stderr, "--injection-file=injfile      XML file with injection parameters\n \t \t should a sim_ringdown table for 'ringdown' injections \n \t \t and a sim_inspiral table for the other types\n" );
  fprintf( stderr, "--inject-mdc-frame=mdcframe  frame file with MDC-frame injections\n" );

  fprintf( stderr, "\n data options \n" );
  fprintf( stderr, "--data-type = type         can be sim (for simulated data), zero (for a data set comprised of zeros), uncal (for uncalibrated data), ht_real4 for single precision strain data, ht_real8 for double precision strain data).\n");
  fprintf( stderr, "--simulated-data           create simulated white Gaussian noise\n" );
  fprintf( stderr, "--random-seed=seed         random number seed for simulated data\n" );
  fprintf( stderr, "--sample-rate=srate        sampling rate of simulated data (Hz)\n" );
  fprintf( stderr, "--calibration-cache=calcache  cache file for calibration frames\n" );
  fprintf( stderr, "--dynamic-range-factor=dynfac  scale calibration by factor dynfac\n" );

  fprintf( stderr, "\ndata segmentation options:\n" );
  fprintf( stderr, "--segment-duration=duration  duration of a data segment (sec) (Subdivisions of analysis block)\n" );
  fprintf( stderr, "--block-duration=duration    duration of an analysis block (sec) (Blocks are subdivided into segments)\n" );
  fprintf( stderr, "--pad-data=duration          input data padding (sec)\n" );

  fprintf( stderr, "\npower spectrum options:\n" );
  fprintf( stderr, "--white-spectrum           use uniform white power spectrum\n" );
  fprintf( stderr, "--cutoff-frequency=fcut    low frequency spectral cutoff (Hz)\n" );
  fprintf( stderr, "--inverse-spec-length=t    set length of inverse spectrum to t seconds\n" );
  fprintf( stderr, "--spectrum-type            specify the algorithm used to calculate the spectrum; must be either median or median_mean\n" );

  fprintf( stderr, "\nbank generation options:\n" );
  fprintf( stderr, "--bank-template-phase=phi  phase of ringdown waveforms (rad, 0=cosine)\n" );
  fprintf( stderr, "--bank-min-quality=qmin    minimum Q of bank\n" );
  fprintf( stderr, "--bank-max-quality=qmax    maximum Q of bank\n" );
  fprintf( stderr, "--bank-min-frequency=fmin  minimum central frequency of bank (Hz)\n" );
  fprintf( stderr, "--bank-max-frequency=fmax  maximum central frequency of bank (Hz)\n" );
  fprintf( stderr, "--bank-max-mismatch=maxmm  maximum template mismatch in bank\n" );
  fprintf( stderr, "--bank-file=name           write template bank to LIGO_LW XML file\n" );
  fprintf( stderr, "--bank-only                generate bank only -- do not read data or filter\n" );

  fprintf( stderr, "\nfiltering options:\n" );
  fprintf( stderr, "--threshold                SNR threshold to identify triggers\n" );
  fprintf( stderr, "--maximize-duration=maxdur  maximize triggers over duration maxdur (sec)\n" );
  fprintf( stderr, "--only-segment-numbers=seglist  list of segment numbers to compute\n" );
  fprintf( stderr, "--only-template-numbers=tmpltlist  list of filter templates to use\n" );

  fprintf( stderr, "\ntrigger output options:\n" );
  fprintf( stderr, "--output-file=outfile      output triggers to file outfile\n" );
  fprintf( stderr, "--trig-start-time=sec      output only triggers after GPS time sec\n" );
  fprintf( stderr, "--trig-end-time=sec        output only triggers before GPS time sec\n" );
  fprintf( stderr, "--ifo-tag=string           set ifotag to string for file naming\n" );
  fprintf( stderr, "--user-tag=string          set the process_params usertag to string\n" );

  fprintf( stderr, "\nintermediate data output options:\n" );
  fprintf( stderr, "--write-raw-data           write raw data before injection or conditioning\n" );
  fprintf( stderr, "--write-data               write data after injection and conditioning\n" );
  fprintf( stderr, "--write-response           write response function used\n" );
  fprintf( stderr, "--write-spectrum           write computed data power spectrum\n" );
  fprintf( stderr, "--write-inv-spectrum       write inverse power spectrum\n" );
  fprintf( stderr, "--write-segment            write overwhitened data segments\n" );
  fprintf( stderr, "--write-template-time-series   write template time series\n");
  fprintf( stderr, "--write-template-fft       write template fft\n");
  fprintf( stderr, "--write-filter-output      write filtered data segments\n" );
  fprintf( stderr, "--write-compress           write a compressed xml file\n");
  return 0;
}
