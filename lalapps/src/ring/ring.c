#include <math.h>
#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <lalapps.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/Units.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/FrameStream.h>
#include <lal/FrameCalibration.h>
#include <lal/RingSearch.h>

#include <lal/PrintFTSeries.h>

#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOLwXML.h>

#include <processtable.h>

RCSID( "$Id$" );

#ifndef HAVE_LIBLALFRAME
int main( void )
{
  lalDebugLevel = 0; /* do something with it so the linker picks it up */
  fputs( "Disabled: LALApps compiled with non-frame-enabled LAL\n", stderr );
  return 77;
}
#else

#define PROGRAM_NAME "ring"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"

extern char *optarg;
extern int optind, opterr, optopt;
extern int vrbflg;

void *my_malloc( size_t size );
void *my_calloc( size_t nobj, size_t size );
char *my_strdup( const char *str );
void reverse_process_params( void );
int add_process_params_prefix( int argc, char **argv, char *prefix );
int add_process_params( int argc, char **argv );
int parse_options( int argc, char **argv );
int import_filter_params( char *fname );
SearchSummaryTable *create_search_summary( LIGOTimeGPS *startepoch,
    REAL8 datadur, REAL8 segdur, UINT4 numseg, UINT4 nevents );
ProcessTable *create_process_table( const char *ifo );
REAL4TimeSeries *get_data( UINT4 segsz );
COMPLEX8FrequencySeries *get_response( UINT4 segsz, double dt, const char *ifo );
int read_response( COMPLEX8FrequencySeries *series, const char *fname );
int vrbmsg( const char *fmt, ... );
int usage( const char *program );

/* global variables: these are the program parameters */
ProcessParamsTable *procpartab;

LIGOTimeGPS tstart;
LIGOTimeGPS tend;
double duration = 64; /* default is to analyze 64 seconds of data */
const char *program;
const char *dbglvl;
const char *frcalib;
const char *frdata;
const char *frpath  = ".";
const char *frptrn  = "*.gwf";
const char *frchan  = "H1:LSC-AS_Q";
const char *rspfile = "response.asc";
const char *outfile;
double srate = -1;
int bmin = -1;
int bmax = -1;
int keep;
int verbose;

const char *fpars;
int         fargc = 1;
char       *fargv[64];


int main( int argc, char *argv[] )
{
  LALStatus status = blank_status;

  static RingSearchData    data;
  static RingSearchInput   input;
  static RingSearchParams *params;
  static SnglBurstTable   *events;

  static MetadataTable     proctable;
  static MetadataTable     procparams;
  static MetadataTable     searchsumm;
  static MetadataTable     ringevents;
  static LIGOLwXMLStream   results;

  char  ofile[64];

  program = argv[0];

  /* add the non-filter arguments to the process param list */
  add_process_params( argc, argv );

  /* parse options */
  parse_options( argc, argv );

  /* reverse process param list so output is in correct order */
  reverse_process_params();
  procparams.processParamsTable = procpartab;

  /* set debug level */
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( dbglvl );

  /* initialize ring search */
  vrbmsg( "call LALRingSearchInit" );
  LAL_CALL( LALRingSearchInit( &status, &params, fargv, fargc ), &status );

  /* set IFO name */
  strncpy( params->ifoName, frchan, 2 );

  /* get data and response function */
  data.channel  = get_data( params->segmentSize );
  data.response = get_response( params->segmentSize, data.channel->deltaT,
      params->ifoName );
  data.spectrum = NULL;

  /* condition data */
  vrbmsg( "call LALRingSearchConditionData" );
  LAL_CALL( LALRingSearchConditionData( &status, params, &data ), &status );

  /* filter entire template bank */
  bmin = bmin < params->templateBank->numTmplt ? bmin :
    params->templateBank->numTmplt;
  bmin = bmin > 0 ? bmin : 0;
  bmax = bmax < params->templateBank->numTmplt ? bmax :
    params->templateBank->numTmplt;
  bmax = bmax >= bmin ? bmax : params->templateBank->numTmplt;
  input.startTemplate = bmin;
  input.templatesToDo = bmax - bmin + 1;
  params->keepResults = keep;
  vrbmsg( "call LALRingSearch" );
  LAL_CALL( LALRingSearch( &status, &events, &input, params ), &status );

  /* output events */

  proctable.processTable = create_process_table( params->ifoName );
  searchsumm.searchSummaryTable = create_search_summary( &tstart,
      duration, params->segmentSize / params->sampleRate,
      params->numSegments, params->numEvents );
  ringevents.snglBurstTable = events;

  /* open results xml file */
  if ( ! outfile )
  {
    LALSnprintf( ofile, sizeof( ofile ), "%s-RING-%d-%d.xml", params->ifoName,
        tstart.gpsSeconds, (int)ceil( duration ) );
    outfile = ofile;
  }
  vrbmsg( "output events to file %s", outfile );
  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &results, outfile ), &status );

  /* output the process table */
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, process_table ),
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, proctable,
        process_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable( &status, &results ), &status );

  /* output process params table */
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, process_params_table ),
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, procparams,
        process_params_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable( &status, &results ), &status );

  /* output search summary table */
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, search_summary_table ),
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, searchsumm,
        search_summary_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable( &status, &results ), &status );

  /* output the events */
  if ( ringevents.snglBurstTable )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, sngl_burst_table ),
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, ringevents,
          sngl_burst_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable( &status, &results ), &status );
  }

  /* close the xml file */
  LAL_CALL( LALCloseLIGOLwXMLFile( &status, &results ), &status );

  /* free event list */
  while ( events )
  {
    SnglBurstTable *next = events->next;
    LALFree( events );
    events = next;
  }

  vrbmsg( "cleaning up" );

  /* deallocate response and channel data vectors */
  LAL_CALL( LALCDestroyVector( &status, &data.response->data ), &status );
  LAL_CALL( LALSDestroyVector( &status, &data.channel->data ), &status );
  free( data.response );
  free( data.channel );

  /* finalize ring search */
  vrbmsg( "call LALRingSearchFini" );
  LAL_CALL( LALRingSearchFini( &status, &params ), &status );

  vrbmsg( "check memory leaks" );
  LALCheckMemoryLeaks();
  return 0;
}

void *my_malloc( size_t size )
{
  void *p;
  p = malloc( size );
  if ( ! p )
  {
    perror( "check_malloc" );
    exit( 1 );
  }
  return p;
}

void *my_calloc( size_t nobj, size_t size )
{
  void *p;
  p = calloc( nobj, size );
  if ( ! p )
  {
    perror( "check_malloc" );
    exit( 1 );
  }
  return p;
}

char *my_strdup( const char *str )
{
  char *s;
  s = my_malloc( strlen( str ) + 1 );
  strcpy( s, str );
  return s;
}

#define is_long_option( s ) \
  ( strlen(s)>2 && (s)[0]=='-' && (s)[1]=='-' && isalpha(s[2]) )
#define is_short_option( s ) \
  ( strlen(s)>1 && (s)[0]=='-' && isalpha(s[1]) )
#define is_option( s ) ( is_long_option( s ) || is_short_option( s ) )
#define is_double_hyphen( s ) ( ! strcmp( s, "--" ) )

void reverse_process_params( void )
{
  ProcessParamsTable *next;
  if ( ! procpartab )
    return;
  next = procpartab->next;
  procpartab->next = NULL;
  while ( next )
  {
    ProcessParamsTable *tmp;
    tmp = next->next;
    next->next = procpartab;
    procpartab = next;
    next = tmp;
  }
  return;
}

int add_process_params_prefix( int argc, char **argv, char *prefix )
{
  int c;
  for ( c = 1; c < argc; ++c )
  {
    char *arg = argv[c];
    char *opt = NULL;
    char *val = NULL;

    /* set option and make sure it is an option */
    opt = my_strdup( arg );
    if ( ! is_option( opt ) )
    {
      fprintf( stderr, "add_process_params: %s is not an option\n", opt );
      exit( 1 );
    }

    /* search for equals sign to identify val */
    if ( ( val = strchr( opt, '=' ) ) )
      *val++ = 0; /* nul terminate opt and set val to value */
    else if ( c < argc - 1 && ! is_option( argv[c+1] ) )/* next arg is value? */
      val = argv[++c]; /* set value and increment counter */
    else /* no value for this option */
      val = NULL;

    /* special cases */
    if ( strstr( opt, "user-tag" ) ) /* change way this is written */
      strncpy( opt, "-userTag", strlen( opt ) );
    if ( strstr( opt, "filter-params" ) ) /* don't write this */
    {
      free( opt );
      opt = NULL;
    }

    /* now write the option and value */
    if ( opt )
    {
      ProcessParamsTable *par = NULL;
      par = my_calloc( 1, sizeof( *par ) );
      par->next = procpartab;
      procpartab = par;
      strncpy( par->program, PROGRAM_NAME, LIGOMETA_PROGRAM_MAX - 1 );
      LALSnprintf( par->param, LIGOMETA_PARAM_MAX, "%s%s",
          prefix ? prefix : "", opt );
      if ( val )
      {
        strncpy( par->type, "string", LIGOMETA_TYPE_MAX - 1 );
        strncpy( par->value, val, LIGOMETA_VALUE_MAX - 1 );
        val = NULL;
      }
      free( opt );
      opt = NULL;
    }
  }
  return c;
}

int add_process_params( int argc, char **argv )
{
  return add_process_params_prefix( argc, argv, NULL );
}

int parse_options( int argc, char **argv )
{
  struct option long_options[] =
  {
    /* these options set a flag */
    { "verbose", no_argument, &vrbflg, 1 },
    /* these options don't set a flag */
    { "help",    no_argument, 0, 'h' },
    { "version", no_argument, 0, 'V' },
    /* these options require an argument */
    { "gps-start-time",          required_argument, 0, 'a' },
    { "gps-start-time-ns",       required_argument, 0, 'A' },
    { "gps-end-time",            required_argument, 0, 'b' },
    { "gps-end-time-ns",         required_argument, 0, 'B' },
    { "channel-name",            required_argument, 0, 'c' },
    { "calibration-frame-cache", required_argument, 0, 'C' },
    { "debug-level",             required_argument, 0, 'd' },
    { "data-frame-cache",        required_argument, 0, 'D' },
    { "frame-files",             required_argument, 0, 'f' },
    { "frame-path",              required_argument, 0, 'F' },
    { "bank-start-template",     required_argument, 0, 'i' },
    { "bank-end-template",       required_argument, 0, 'j' },
    { "output-file",             required_argument, 0, 'o' },
    { "response-file",           required_argument, 0, 'r' },
    { "sample-rate",             required_argument, 0, 's' },
    { "user-tag",                required_argument, 0, 'U' },
    { "userTag",                 required_argument, 0, 'U' },
    /* these are filter parameters */
    { "filter-params",           required_argument, 0, 'z' },
    { "filter-segsz",            required_argument, 0, 'Z' },
    { "filter-speclen",          required_argument, 0, 'Z' },
    { "filter-flow",             required_argument, 0, 'Z' },
    { "filter-fmin",             required_argument, 0, 'Z' },
    { "filter-fmax",             required_argument, 0, 'Z' },
    { "filter-qmin",             required_argument, 0, 'Z' },
    { "filter-qmax",             required_argument, 0, 'Z' },
    { "filter-maxmm",            required_argument, 0, 'Z' },
    { "filter-thresh",           required_argument, 0, 'Z' },
    { "filter-scale",            required_argument, 0, 'Z' }
  };
  char args[] = "hVa:A:b:B:c:C:d:D:f:F:i:j:o:r:s:U:z:Z:";

  program  = argv[0];
  fargv[0] = my_strdup( program );

  while ( 1 )
  {
    int option_index = 0; /* getopt_long stores long option here */
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
        {
          fprintf( stderr, "error parsing option %s with argument %s\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
      case 'h': /* help */
        usage( program );
        exit( 0 );
      case 'V': /* version */
        PRINT_VERSION( "ring" );
        exit( 0 );
      case 'a': /* gps-start-time */
        tstart.gpsSeconds = atol( optarg );
        break;
      case 'A': /* gps-start-time-ns */
        tstart.gpsNanoSeconds = atol( optarg );
        break;
      case 'b': /* gps-end-time */
        tend.gpsSeconds = atol( optarg );
        break;
      case 'B': /* gps-end-time-ns */
        tend.gpsNanoSeconds = atol( optarg );
        break;
      case 'c': /* channel-name */
        frchan = optarg;
        break;
      case 'C': /* calibration-frame-cache */
        frcalib = optarg;
        break;
      case 'd': /* debug-level */
        dbglvl = optarg;
        break;
      case 'D': /* data-frame-cache */
        frdata = optarg;
        break;
      case 'f': /* frame-files */
        frptrn = optarg;
        break;
      case 'F': /* frame-path */
        frpath = optarg;
        break;
      case 'i': /* bank-start-template */
        bmin = atol( optarg );
        break;
      case 'j': /* bank-end-template */
        bmax = atol( optarg );
        break;
      case 'o': /* output-file */
        outfile = optarg;
        break;
      case 'r': /* response-file */
        rspfile = optarg;
        break;
      case 's': /* sample-rate */
        srate = atof( optarg );
        break;
      case 'U': /* user-tag: do nothing with this! */
        break;
      case 'z': /* filter-params: inport a filter params file */
        import_filter_params( optarg );
        add_process_params_prefix( fargc, fargv, "--filter" );
        break;
      case 'Z': /* filter-something: put this in fargc */
        fargv[fargc++] = my_strdup(
            strstr( "filter", long_options[option_index].name ) + 6 );
        fargv[fargc++] = my_strdup( optarg );
        break;
      case '?':
        exit( 1 );
      default:
        fprintf( stderr, "unknown error while parsing options\n" );
        exit( 1 );
    }
  }

  if ( optind < argc )
  {
    fprintf( stderr, "extraneous command line arguments:\n" );
    while ( optind < argc )
      fprintf( stderr, "%s\n", argv[optind++] );
    exit( 1 );
  }

  fargv[fargc] = NULL;
  return 0;
}

int import_filter_params( char *fname )
{
  char  line[1024];
  FILE *fp;
  fp = fname ? fopen( fname, "r" ) : stdin;
  if ( ! fp )
  {
    perror( "import_filter_params" );
    exit( 1 );
  }
  while ( fgets( line, sizeof( line ), fp ) )
  {
    char *opt;
    char *val;
    if ( *line == '#' )
      continue;
    opt = strchr( line, '-' );
    val = opt + strcspn( opt, " =" ); /* span to end of option */
    *val++ = 0; /* nul-terminate opt and point val to option value */
    val += strspn( val, " " ); /* span any spaces */
    val[strcspn( val, " #" )] = 0; /* span to end of value and terminate */
    fargv[fargc++] = my_strdup( opt );
    fargv[fargc++] = my_strdup( val );
  }
  fname ? fclose( fp ) : 0;
  return 0;
}

SearchSummaryTable *create_search_summary( LIGOTimeGPS *startepoch,
    REAL8 datadur, REAL8 segdur, UINT4 numseg, UINT4 nevents )
{
  SearchSummaryTable *summ = NULL;
  INT8 instart;
  INT8 inend;
  INT8 outstart;
  INT8 outend;

  /* setup search summary table */
  summ = my_calloc( 1, sizeof( *summ ) );
  strncpy( summ->comment, program, LIGOMETA_COMMENT_MAX );
  summ->nnodes = 1;

  instart   = (INT8)startepoch->gpsSeconds * (INT8)(1000000000);
  instart  += (INT8)startepoch->gpsNanoSeconds;
  inend     = instart;
  inend    += (INT8)( 1e9 * datadur );
  /* input and output start times differ by a quarter of a segment */
  /* (lost time that is ignored from beginning of first segment) */
  outstart  = instart;
  outstart += (INT8)( 1e9 * 0.25 * segdur );
  /* data analyzed is numseg times half a segment duration (they overlap) */
  outend    = instart;
  outend   += (INT8)( 1e9 * 0.5 * numseg * segdur );

  /* store input start time and end time of raw data in search summary */
  summ->in_start_time = *startepoch;
  summ->in_end_time.gpsSeconds        = inend / 1000000000;
  summ->in_end_time.gpsNanoSeconds    = inend % 1000000000;
  summ->out_start_time.gpsSeconds     = outstart / 1000000000;
  summ->out_start_time.gpsNanoSeconds = outstart % 1000000000;
  summ->out_end_time.gpsSeconds       = outend / 1000000000;
  summ->out_end_time.gpsNanoSeconds   = outend % 1000000000;

  summ->nevents = nevents;

  return summ;
}

ProcessTable *create_process_table( const char *ifo )
{
  LALLeapSecAccuracy accuracy = LALLEAPSEC_LOOSE;
  LALStatus status = blank_status;
  ProcessTable *proc = NULL;
  proc = my_calloc( 1, sizeof( *proc ) );

  LAL_CALL( populate_process_table( &status, proc, PROGRAM_NAME,
        CVS_REVISION, CVS_SOURCE, CVS_DATE ), &status );
  strncpy( proc->comment, " ", LIGOMETA_COMMENT_MAX );
  strncpy( proc->ifos, ifo, LIGOMETA_IFOS_MAX );
  LAL_CALL( LALGPSTimeNow( &status, &proc->end_time, &accuracy ), &status );

  return proc;
}

static FrChanIn blank_chanin;
REAL4TimeSeries *get_data( UINT4 segsz )
{
  LALStatus  status = blank_status;
  FrChanIn   chanin = blank_chanin;
  FrStream  *stream = NULL;
  REAL4TimeSeries *channel;
  UINT4     npts;

  /* allocate memory */
  channel = my_calloc( 1, sizeof( *channel ) );

  /* get frame data */
  chanin.name = frchan;
  chanin.type = ADCDataChannel;
  if ( frdata ) /* open data cache to get frame files */
  {
    FrCache *cache = NULL;
    vrbmsg( "get data from cache file %s", frdata );
    LAL_CALL( LALFrCacheImport( &status, &cache, frdata ), &status );
    LAL_CALL( LALFrCacheOpen( &status, &stream, cache ), &status );
    LAL_CALL( LALDestroyFrCache( &status, &cache ), &status );
  }
  else /* get frame files from specified path and pattern */
  {
    vrbmsg( "get data from frame files %s/%s", frpath, frptrn );
    LAL_CALL( LALFrOpen( &status, &stream, frpath, frptrn ), &status );
  }
  if ( tstart.gpsSeconds || tstart.gpsNanoSeconds )
  {
    vrbmsg( "seek to gps time %d.%09d",
        tstart.gpsSeconds, tstart.gpsNanoSeconds );
    LAL_CALL( LALFrSeek( &status, &tstart, stream ), &status );
    if ( tend.gpsSeconds || tend.gpsNanoSeconds )
    {
      INT8 ns;
      ns  = tstart.gpsSeconds - tend.gpsSeconds;
      ns *= 1000000000;
      ns += tstart.gpsNanoSeconds - tend.gpsNanoSeconds;
      if ( ns > 0 )
        duration = 1e-9 * ns;
    }
  }

  /* get the data */
  vrbmsg( "read %f seconds of data", duration );

  /* call first time to get sample rate */
  LAL_CALL( LALFrGetREAL4TimeSeries( &status, channel, &chanin, stream ),
      &status );
  /* compute how much data we need to get and allocate memory */
  npts = duration / channel->deltaT;
  LAL_CALL( LALSCreateVector( &status, &channel->data, npts ), &status );
  /* get the data */
  LAL_CALL( LALFrGetREAL4TimeSeries( &status, channel, &chanin, stream ),
      &status );
  /* close the frame file */
  LAL_CALL( LALFrClose( &status, &stream ), &status );
  /* set epoch and duration to actual start time and duration */
  tstart   = channel->epoch;
  duration = channel->deltaT * channel->data->length;
  /* resample frame data if required */
  if ( srate > 0 && srate * channel->deltaT < 1 )
  {
    ResampleTSParams resamplepar;
    memset( &resamplepar, 0, sizeof( resamplepar ) );
    resamplepar.deltaT     = 1.0 / srate;
    resamplepar.filterType = defaultButterworth;
    vrbmsg( "resampling data from %f Hz to %f Hz",
        1.0 / channel->deltaT, srate );
    LAL_CALL( LALResampleREAL4TimeSeries( &status, channel, &resamplepar ),
        &status );
  }
  /* compute how many points we will actually analyze */
  npts  = channel->data->length / ( segsz / 2 );
  npts *= segsz / 2;
  if ( channel->data->length > npts ) /* resize data array */
  {
    channel->data->length = npts;
    channel->data->data   = LALRealloc( channel->data->data,
        npts * sizeof( *channel->data->data ) );
    duration = channel->deltaT * channel->data->length;
    vrbmsg( "resizing data duration to %f seconds", duration );
  }
  return channel;
}

COMPLEX8FrequencySeries *get_response( UINT4 segsz, double dt, const char *ifo )
{
  COMPLEX8FrequencySeries *response;
  LALStatus status = blank_status;
  RAT4 minusOne    = { -1, 0 };
  LALUnitPair unitPair;
  LALUnit     unit;

  response = my_calloc( 1, sizeof( *response ) );
  LAL_CALL( LALCCreateVector( &status, &response->data, segsz / 2 + 1 ),
      &status );
  response->f0     = 0.0;
  response->deltaF = 1.0 / ( segsz * dt );
  LAL_CALL( LALUnitRaise( &status, &unit, &lalADCCountUnit, &minusOne ),
      &status );
  unitPair.unitOne = &lalStrainUnit;
  unitPair.unitTwo = &unit;
  LAL_CALL( LALUnitMultiply( &status, &response->sampleUnits, &unitPair ),
      &status );
  strncpy( response->name, "response", sizeof( response->name ) );
  response->epoch = tstart;

  if ( frcalib ) /* use calibration cache file to extract & update response */
  {
    vrbmsg( "get calibration data from frames in cache file %s", frcalib );
    LAL_CALL( LALExtractFrameResponse( &status, response, frcalib, ifo ),
        &status );
  }
  else /* get fixed response function from an ascii file */
  {
    vrbmsg( "get calibration data from ascii file %s", rspfile );
    read_response( response, rspfile );
  }
  return response;
}


int read_response( COMPLEX8FrequencySeries *series, const char *fname )
{
  char line[1024];
  FILE *fp;
  size_t n;
  size_t i, j;
  double *freq;
  double *re;
  double *im;
  fp = fopen( fname, "r" );

  n = 0;
  while ( fgets( line, sizeof( line ), fp ) )
    if ( *line != '#' )
      ++n;
  rewind( fp );

  freq = LALCalloc( n, sizeof( *freq ) );
  re   = LALCalloc( n, sizeof( *re ) );
  im   = LALCalloc( n, sizeof( *im ) );

  for ( i = 0; i < n; ++i )
  {
    fgets( line, sizeof( line ), fp );
    if ( *line == '#' )
      continue;
    sscanf( line, "%le %le %le %*e\n", freq + i, re + i, im + i );
  }

  i = 1;
  for ( j = 0; j < series->data->length; ++j )
  {
    double f, x;
    f = j * series->deltaF;
    while ( i < n )
      if ( freq[i] > f )
        break;
      else
        ++i;
    x = ( f - freq[i - 1] ) / ( freq[i] - freq[i - 1] );
    series->data->data[j].re = re[i - 1] + x * ( re[i] - re[i - 1] );
    series->data->data[j].im = im[i - 1] + x * ( im[i] - im[i - 1] );
  }

  LALFree( freq );
  LALFree( re );
  LALFree( im );
  
  fclose( fp );
  return 0;
}

int vrbmsg( const char *fmt, ... )
{
  int c = 0;
  if ( vrbflg )
  {
    va_list ap;
    va_start( ap, fmt );
    fputs( "verbose: ", stderr );
    c = vfprintf( stderr, fmt, ap );
    fputs( "\n", stderr );
    va_end( ap );
  }
  return c;
}

const char *usgfmt =
"Usage: %s [options] [filterparams]\n"
"Options [default in brackets]:\n\n"
"  --help\tprint this message\n\n"
"  --version\tprint version info\n\n"
"  --verbose\tverbose\n\n"
"  --debug-level dbglvl\n\t\tset debug level to dbglvl [0]\n\n"
"  --gps-start-time startsec\n\t\tstart time in GPS seconds [start of data]\n\n"
"  --gps-start-time-ns startnan\n\t\tstart time in GPS nano-seconds [0]\n\n"
"  --gps-end-time endsec\n\t\tend time in GPS seconds [start + 64]\n\n"
"  --gps-end-time-ns endnan\n\t\tend time in GPS nano-seconds [0]\n\n"
"  --channel-name channel\n\t\tchannel to analyze [H1:LSC-AS_Q]\n\n"
"  --calibration-frame-cache calcache\n\t\tcalib frame cache file [use ascii]\n\n"
"  --response-file respfile\n\t\tascii response file [response.asc]\n\n"
"  --data-frame-cache datcache\n\t\tdata frame cache file [don't use cache]\n\n"
"  --frame-files pattern\n\t\tframe file to use [*.gwf]\n\n"
"  --frame-path path\n\t\tpath to look for frame files [.]\n\n"
"  --bank-start-template bmin\n\t\tfirst template of bank to use [0]\n\n"
"  --bank-end-tempate bmax\n\t\tlast template of bank to use [last]\n\n"
"  --output-file outfile\n\t\tevent output file [IFO-RING-tstart-duration.xml]\n\n"
"  --sample-rate srate\n\t\tdesired sampling rate in Hz [raw rate]\n\n"
"  --user-tag tag\n"
"  --userTag  tag\n\t\tuser tag [none]\n\n"
"Filter Parameters [defaults in brackets; otherwise required]\n\n"
"  --filter-params parfile\n\t\tread filter parameters from file [none]\n\n"
"  --filter-segsz npts\n\t\tset size of segments to analyze to npts\n\n"
"  --filter-speclen len\n\t\tset size of inverse spectrum truncation to len [0]\n\n"
"  --filter-flow flow\n\t\tset low frequency cutoff to flow (Hz)\n\n"
"  --filter-fmin fmin\n\t\tset minimum frequency for bank to fmin (Hz)\n\n"
"  --filter-fmax fmax\n\t\tset maximum frequency for bank to fmax (Hz)\n\n"
"  --filter-qmin qmin\n\t\tset minimum quality for bank to qmin\n\n"
"  --filter-qmax qmax\n\t\tset maximum quality for bank to qmax\n\n"
"  --filter-maxmm maxmm\n\t\tset maximum mismatch of bank to maxmm\n\n"
"  --filter-thresh thrsh\n\t\tset ringdown event snr threshold to thrsh\n\n"
"  --filter-scale scale\n\t\tscale response function by a factor of scale [1]\n\n";

int usage( const char *program )
{
  return fprintf( stderr, usgfmt, program );
}

#endif
