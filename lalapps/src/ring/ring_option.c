#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>

#include <lal/LALStdio.h>
#include "lalapps.h"
#include "errutil.h"
#include "gpstime.h"
#include "ring.h"

RCSID( "$Id$" );

extern int vrbflg;
static int ring_usage( const char *program );
static int ring_default_params( struct ring_params *params );

int ring_parse_options( struct ring_params *params, int argc, char **argv )
{
  struct option long_options[] =
  {
    { "verbose", no_argument, &vrbflg, 1 },
    { "geo-data", no_argument, &params->geoData, 1 },
    { "strain-data", no_argument, &params->strainData, 1 },
    { "simulated-data", no_argument, &params->simData, 1 },
    { "white-spectrum", no_argument, &params->whiteSpectrum, 1 },
    { "bank-only", no_argument, &params->bankOnly, 1 },
    { "write-raw-data",     no_argument, &params->writeRawData, 1 },
    { "write-data",         no_argument, &params->writeProcessedData, 1 },
    { "write-response",     no_argument, &params->writeResponse, 1 },
    { "write-spectrum",     no_argument, &params->writeSpectrum, 1 },
    { "write-inv-spectrum", no_argument, &params->writeInvSpectrum, 1 },
    { "write-bank",         no_argument, &params->writeBank, 1 },
    { "write-segment",      no_argument, &params->writeSegment, 1 },
    { "write-filter-output",no_argument, &params->writeFilterOutput, 1 },
    { "help",                    no_argument,       0, 'h' },
    { "version",                 no_argument,       0, 'V' },
    { "gps-start-time",          required_argument, 0, 'a' },
    { "gps-start-time-ns",       required_argument, 0, 'A' },
    { "gps-end-time",            required_argument, 0, 'b' },
    { "gps-end-time-ns",         required_argument, 0, 'B' },
    { "channel-name",            required_argument, 0, 'c' },
    { "calibration-cache",       required_argument, 0, 'C' },
    { "debug-level",             required_argument, 0, 'd' },
    { "data-cache",              required_argument, 0, 'D' },
    { "cutoff-frequency",        required_argument, 0, 'e' },
    { "highpass-frequency",      required_argument, 0, 'E' },
    { "bank-min-frequency",      required_argument, 0, 'f' },
    { "bank-max-frequency",      required_argument, 0, 'F' },
    { "geo-highpass-frequency",  required_argument, 0, 'g' },
    { "geo-data-scale",          required_argument, 0, 'G' },
    { "inject-file",             required_argument, 0, 'i' },
    { "inject-mdc-frame",        required_argument, 0, 'I' },
    { "bank-max-mismatch",       required_argument, 0, 'm' },
    { "maximize-duration",       required_argument, 0, 'M' },
    { "only-segment-numbers",    required_argument, 0, 'n' },
    { "only-template-numbers",   required_argument, 0, 'N' },
    { "output-file",             required_argument, 0, 'o' },
    { "output-format",           required_argument, 0, 'O' },
    { "bank-template-phase",     required_argument, 0, 'p' },
    { "bank-min-quality",        required_argument, 0, 'q' },
    { "bank-max-quality",        required_argument, 0, 'Q' },
    { "random-seed",             required_argument, 0, 'r' },
    { "dynamic-range-factor",    required_argument, 0, 'R' },
    { "sample-rate",             required_argument, 0, 's' },
    { "segment-duration",        required_argument, 0, 'S' },
    { "threshold",               required_argument, 0, 't' },
    { 0, 0, 0, 0 }
  };
  char args[] = "a:A:b:B:c:C:d:D:e:E:f:F:g:G:hi:I:m:o:O:p:q:Q:r:R:s:S:t:V";
  char *program = argv[0];

  /* set default values for parameters before parsing arguments */
  ring_default_params( params );

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
        params->startTime.gpsSeconds = atol( optarg );
        break;
      case 'A': /* gps-start-time-ns */
        params->startTime.gpsNanoSeconds = atol( optarg );
        break;
      case 'b': /* gps-end-time */
        params->endTime.gpsSeconds = atol( optarg );
        break;
      case 'B': /* gps-end-time-ns */
        params->endTime.gpsNanoSeconds = atol( optarg );
        break;
      case 'c': /* channel-name */
        params->channel = optarg;
        break;
      case 'C': /* calibration-cache */
        params->calibCache = optarg;
        break;
      case 'd': /* debug-level */
        set_debug_level( optarg );
        break;
      case 'D': /* data-cache */
        params->dataCache = optarg;
        break;
      case 'e': /* cutoff-frequency */
        params->lowCutoffFrequency = atof( optarg );
        break;
      case 'E': /* highpass-frequency */
        params->highpassFrequency = atof( optarg );
        break;
      case 'f': /* bank min frequency */
        params->bankParams.minFrequency = atof( optarg );
        break;
      case 'F': /* bank max frequency */
        params->bankParams.maxFrequency = atof( optarg );
        break;
      case 'g': /* geo-highpass-frequency */
        params->geoHighpassFrequency = atof( optarg );
        break;
      case 'G': /* geo-data-scale */
        params->geoScale = atof( optarg );
        error( "currently geo scale is not implemented correctly\n" );
        break;
      case 'h': /* help */
        ring_usage( program );
        exit( 0 );
      case 'i': /* inject-file */
        params->injectFile = optarg;
        break;
      case 'I': /* inject-mdc-frame */
        error( "currently unsupported option: --inject-mdc-frame\n" );
        break;
      case 'm': /* bank max mismatch */
        params->bankParams.maxMismatch = atof( optarg );
        break;
      case 'M': /* maximize duration */
        params->maximizeEventDuration = atof( optarg );
        break;
      case 'n': /* only-segment-numbers */
        params->segmentsToDoList = optarg;
        break;
      case 'N': /* only-template-number */
        params->templatesToDoList = optarg;
        break;
      case 'o': /* output-file */
        strncpy( params->outputFile, optarg, sizeof( params->outputFile ) - 1 );
        break;
      case 'O': /* output-file */
        if ( strstr( optarg, "xml" ) || strstr( optarg, "XML" ) )
          params->outputFormat = output_ligolw;
        else if ( strstr( optarg, "ligo" ) || strstr( optarg, "LIGO" ) )
          params->outputFormat = output_ligolw;
        else if ( strstr( optarg, "asc" ) || strstr( optarg, "ASC" ) )
          params->outputFormat = output_ascii;
        else
          error( "unknown output format %s: "
              "expect either \"xml\" or \"asc\"\n", optarg );
        break;
      case 'p': /* bank template phase */
        params->bankParams.templatePhase = atof( optarg );
        break;
      case 'q': /* bank min quality */
        params->bankParams.minQuality = atof( optarg );
        break;
      case 'Q': /* bank max quality */
        params->bankParams.maxQuality = atof( optarg );
        break;
      case 'r': /* data-cache */
        params->randomSeed = atoi( optarg );
        break;
      case 'R': /* dynamic range factor */
        params->dynRangeFac = atof( optarg );
        break;
      case 's': /* data-cache */
        params->sampleRate = atof( optarg );
        break;
      case 'S': /* segment-duration */
        params->segmentDuration = atof( optarg );
        break;
      case 't': /* threshold */
        params->threshold = atof( optarg );
        break;
      case 'V': /* version */
        PRINT_VERSION( "ring" );
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

  return 0;
}

static int ring_default_params( struct ring_params *params )
{
  /* overall, default values are zero */
  memset( params, 0, sizeof( *params ) );

  /* dynamic range factor must be greater than zero */
  params->dynRangeFac = 1.0;
  params->geoScale    = 1.0;

  /* negative value means use the "default" values */
  params->highpassFrequency     = -1.0; /* use low-frequency cutoff */
  params->geoHighpassFrequency  = -1.0; /* use low-frequency cutoff */
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

  return 0;
}

#define sanity_check( condition ) \
  ( condition ? 0 : ( fputs( #condition " not satisfied\n", stderr ), error( "sanity check failed\n" ) ) ) 
int ring_params_sanity_check( struct ring_params *params )
{
  UINT4 recordLength = 0;
  UINT4 segmentLength = 0;
  UINT4 segmentStride = 0;
  UINT4 truncateLength = 0;
  INT8 startTime;
  INT8 endTime;
  int validChannelIFO;

  if ( params->geoData )
    params->strainData = 1;

  if ( params->strainData )
      params->getResponse = 0;

  if ( params->bankOnly )
  {
    params->getData     = 0;
    params->getResponse = 0;
    params->getSpectrum = 0;
    params->doFilter    = 0;
  }


  if ( params->getSpectrum ) /* need data and response if not strain data */
    sanity_check( params->getData && (params->strainData || params->getResponse) );
  if ( params->doFilter ) /* need data, bank, and spectrum */
    sanity_check( params->getData && params->getBank && params->getSpectrum );

  if ( params->getData )
  {
    /* checks on data duration */
    startTime = epoch_to_ns( &params->startTime );
    endTime   = epoch_to_ns( &params->endTime );
    sanity_check( startTime > 0 );
    sanity_check( endTime > startTime );
    params->duration = 1e-9*(endTime - startTime);
    sanity_check( params->duration > 0 );

    /* checks on size of data record */
    sanity_check( params->sampleRate > 0 );
    recordLength = params->duration * params->sampleRate;
    sanity_check( recordLength > 0 );

    sanity_check( params->channel );
    sanity_check( params->simData || params->dataCache );

    /* record ifo name */
    validChannelIFO = sscanf( params->channel, "%2[A-Z1-9]", params->ifoName );
    sanity_check( validChannelIFO );
    sanity_check( params->ifoName );

    /* will need response to do injections unless strain data */
    sanity_check( params->injectFile == NULL || params->strainData || params->getResponse );

    if ( params->geoData )
    {
      sanity_check( params->geoScale > 0.0 );
      if ( params->geoHighpassFrequency < 0.0 )
      {
        if ( params->highpassFrequency < 0.0 )
          params->geoHighpassFrequency = params->highpassFrequency;
        else
          params->geoHighpassFrequency = params->lowCutoffFrequency;
      }
      sanity_check( params->geoHighpassFrequency > 0.0 );
    }
  }

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
    sanity_check( ! (params->numOverlapSegments % 2) ); /* required to be even for median-mean method */

    /* checks on data input information */
    sanity_check( params->channel );
    sanity_check( params->dynRangeFac > 0.0 );
  }

  if ( params->getResponse )
  {
    sanity_check( params->calibCache );
  }

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
      LALSnprintf( params->outputFile, sizeof( params->outputFile ),
          "%s-RING-%d-%d.%s", params->ifoName, params->startTime.gpsSeconds,
          (int)ceil( params->duration ),
          params->outputFormat == output_ligolw ? "xml" : "asc" );
  }

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

static int ring_usage( const char *program )
{
  fprintf( stderr, "usage: %s options\n", program );
  fprintf( stderr, "options:\n" );

  return 0;
}
