/*
 *
 * Program: lalapps_ring -- search code for ringdown waveforms
 * Authors: Jolien Creighton, Rana Adhikari
 *
 */
#include <math.h>
#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>

/* lal headers */
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/Random.h>
#include <lal/Units.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/GenerateBurst.h>
#include <lal/FrameStream.h>
#include <lal/FrameCalibration.h>
#include <lal/RingSearch.h>

/* lal-support headers */
#include <lal/PrintFTSeries.h>

/* lal-metaio headers */
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>

/* frame library headers */
#include <FrameL.h>

/* lalapps headers */
#include <lalapps.h>
#include <series.h>
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
REAL4TimeSeries *get_data( UINT4 segsz, const char *ifo );
COMPLEX8FrequencySeries *get_response( UINT4 segsz, double dt, const char *ifo );
int read_response( COMPLEX8FrequencySeries *series, const char *fname );
FrameH *fr_add_proc_REAL4TimeSeries( FrameH *frame,
    REAL4TimeSeries *series );
FrameH *fr_add_proc_REAL4FrequencySeries( FrameH *frame,
    REAL4FrequencySeries *series );
FrameH *fr_add_proc_COMPLEX8FrequencySeries( FrameH *frame,
    COMPLEX8FrequencySeries *series );
int vrbmsg( const char *fmt, ... );
int usage( const char *program );

/* possible output formats */
enum { ascii_format, xml_format, frame_format };

/* global variables: these are the program parameters */
ProcessParamsTable *procpartab;

CalibrationUpdateParams calfacts;
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
char       *injfile;
int outfmt = xml_format;
double srate = -1;
int bmin = -1;
int bmax = -1;
int verbose;
int strain_data;
int geo_data;

/* flags for writing output */
int write_raw_data;
int write_data;
int write_response;
int write_inverse_spectrum;
int write_data_segments;
int write_filter_output;
int write_format = frame_format;
FrameH *write_frame;

/* flags for internal tests */
int test_zero_data;
int test_gaussian_data;
int test_white_spectrum;
int test_unit_response;
int test_inject;

/* internal test injection parameters */
int test_inject_sec;
int test_inject_nan;
float test_inject_freq;
float test_inject_qual;
float test_inject_ampl;
float test_inject_phase;

/* geo high pass filter parameters */
float geoHighPassFreq = 70.;
float geoHighPassOrder = 8;
float geoHighPassAtten = 0.1;
float geoDynamicRange = 65.;

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

  char  ofile[64];

  program = argv[0];

  /* add the non-filter arguments to the process param list */
  add_process_params( argc, argv );

  /* reverse process param list so output is in correct order */
  reverse_process_params();

  /* parse options */
  parse_options( argc, argv );
  if (geo_data) {strain_data = 1;}

  /* echo command if verbose */
  if ( vrbflg )
  {
    int arg;
    fprintf( stderr, "command: %s", program );
    for ( arg = 1; arg < argc; ++arg )
      fprintf( stderr, " %s", argv[arg] );
    fprintf( stderr, "\n" );
  }

  /* set debug level */
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( dbglvl );

  /* initialize ring search */
  vrbmsg( "call LALRingSearchInit" );
  LAL_CALL( LALRingSearchInit( &status, &params, fargv, fargc ), &status );

  /* set IFO name and other filter (test) parameters */
  strncpy( params->ifoName, frchan, 2 );
  if ( test_zero_data )
    params->testZeroData = 1;
  if ( test_gaussian_data && test_white_spectrum ) /* use unity spec method */
  {
    if ( ! ( srate > 0 ) )
    {
      fprintf( stderr, "sample rate must be set (use --sample-rate)\n" );
      exit( 1 );
    }
    params->avgSpecMeth = useUnity;
    params->avgSpecNorm = 2.0 / srate;
  }
  if ( test_inject ) /* injection test parameters */
  {
    params->testInject = 1;
    params->testInjectTime.gpsSeconds     = test_inject_sec;
    params->testInjectTime.gpsNanoSeconds = test_inject_nan;
    params->testInjectFreq                = test_inject_freq;
    params->testInjectQual                = test_inject_qual;
    params->testInjectAmpl                = test_inject_ampl;
    params->testInjectPhase               = test_inject_phase;
  }

  /* get data and response function; write to files if requested */
  data.channel  = get_data( params->segmentSize, params->ifoName );
  data.response = get_response( params->segmentSize, data.channel->deltaT,
      params->ifoName );
  data.spectrum = NULL;

  /* condition data */
  vrbmsg( "call LALRingSearchConditionData" );
  LAL_CALL( LALRingSearchConditionData( &status, params, &data ), &status );

  /* write data to file if required */
  if ( write_data )
  {
    if ( write_format == frame_format )
    {
      LALSnprintf( data.channel->name, sizeof( data.channel->name ), "%s_PROC", frchan );
      vrbmsg( "writing data to frame" );
      write_frame = fr_add_proc_REAL4TimeSeries( write_frame, data.channel );
    }
    else
    {
      const char *fname = "ring-proc-data.dat";
      vrbmsg( "writing data to file %s", fname );
      LALSPrintTimeSeries( data.channel, fname );
    }
  }

  /* write inverse spectrum and/or data segments if required */
  if ( write_inverse_spectrum )
  {
    LALSnprintf( params->invSpectrum->name, sizeof( params->invSpectrum->name ),
        "%s_INVSPEC", frchan );
    if ( write_format == frame_format )
    {
      vrbmsg( "writing inverse spectrum to frame" );
      write_frame = fr_add_proc_REAL4FrequencySeries( write_frame,
          params->invSpectrum );
    }
    else
    {
      const char *fname = "ring-invspec.dat";
      vrbmsg( "writing inverse spectrum to file %s", fname );
      LALSPrintFrequencySeries( params->invSpectrum, fname );
    }
  }
  if ( write_data_segments )
  {
    size_t seg;
    for ( seg = 0; seg < params->numSegments; ++seg )
    {
      LALSnprintf( params->dataSegment[seg].name,
          sizeof( params->dataSegment[seg].name ), "%s_SEG_%03d", frchan, seg );
      if ( write_format == frame_format )
      {
        vrbmsg( "writing data segment %d to frame", seg );
        write_frame = fr_add_proc_COMPLEX8FrequencySeries( write_frame,
            params->dataSegment + seg );
      }
      else
      {
        char fname[64];
        LALSnprintf( fname, sizeof( fname ), "ring-segment-%03d.dat", seg );
        vrbmsg( "writing data segment %d to file %s", seg, fname );
        LALCPrintFrequencySeries( params->dataSegment + seg, fname );
      }
    }
  }

  /* filter entire template bank */
  bmin = bmin < params->templateBank->numTmplt ? bmin :
    params->templateBank->numTmplt;
  bmin = bmin > 0 ? bmin : 0;
  bmax = bmax < params->templateBank->numTmplt ? bmax :
    params->templateBank->numTmplt;
  bmax = bmax >= bmin ? bmax : params->templateBank->numTmplt;
  input.startTemplate = bmin;
  input.templatesToDo = bmax - bmin + 1;
  params->keepResults = write_filter_output; /* keep output if to write it */
  vrbmsg( "call LALRingSearch" );
  LAL_CALL( LALRingSearch( &status, &events, &input, params ), &status );

  /* write filter outputs if required */
  if ( write_filter_output )
  {
    size_t rslt;
    for ( rslt = 0; rslt < params->numResults; ++rslt )
    {
      int tmplt;
      int sgmnt;
      sscanf( params->result[rslt].name, "snr-%d.%d", &tmplt, &sgmnt );
      LALSnprintf( params->result[rslt].name,
          sizeof( params->result[rslt].name ), "%s_SNR_%03d_%03d",
          frchan, tmplt, sgmnt );
      if ( write_format == frame_format )
      {
        vrbmsg( "writing filter %d output to frame", rslt );
        write_frame = fr_add_proc_REAL4TimeSeries( write_frame,
            params->result + rslt );
      }
      else
      {
        char fname[64];
        LALSnprintf( fname, sizeof(fname), "snr-%03d-%03d.dat", tmplt, sgmnt );
        vrbmsg( "writing filter %d output to file %s", rslt, fname );
        LALSPrintTimeSeries( params->result + rslt, fname );
      }
    }
  }

  /* actually write frame results if required */
  if ( write_frame )
  {
    FrFile *frfile;
    char fname[64];
    LALSnprintf( fname, sizeof( fname ), "%s-RING-%d-%d.gwf", params->ifoName,
        tstart.gpsSeconds, (int)ceil( duration ) );
    vrbmsg( "writing frame file %s", fname );
    frfile = FrFileONew( fname, 0 );
    FrameWrite( write_frame, frfile );
    FrFileOEnd( frfile );
  }

  /* output events */
  if ( ! outfile )
  {
    LALSnprintf( ofile, sizeof( ofile ), "%s-RING-%d-%d.%s", params->ifoName,
        tstart.gpsSeconds, (int)ceil( duration ),
        outfmt == xml_format ? "xml" : "asc" );
    outfile = ofile;
  }
  vrbmsg( "output events to file %s", outfile );

  if ( outfmt == xml_format )  /* XML output */
  {
    MetadataTable proctable;
    MetadataTable procparams;
    MetadataTable searchsumm;
    MetadataTable ringevents;
    LIGOLwXMLStream  results;

    memset( &proctable, 0, sizeof( proctable ) );
    memset( &procparams, 0, sizeof( procparams ) );
    memset( &searchsumm, 0, sizeof( searchsumm ) );
    memset( &ringevents, 0, sizeof( ringevents ) );
    memset( &results, 0, sizeof( results ) );

    proctable.processTable = create_process_table( params->ifoName );
    procparams.processParamsTable = procpartab;
    searchsumm.searchSummaryTable = create_search_summary( &tstart,
        duration, params->segmentSize / params->sampleRate,
        params->numSegments, params->numEvents );
    ringevents.snglBurstTable = events;

    /* open results xml file */
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
  }
  else /* ASCII output */
  {
    SnglBurstTable *this_event = events;
    FILE *fp;
    fp = fopen( outfile, "w" );
    if ( ! fp )
    {
      perror( "output file" );
      exit( 1 );
    }
    fprintf( fp,
        "# gps start time\tsignal/noise\tamplitude\tfrequency\tbandwidth\n" );
    while ( this_event )
    {
      fprintf( fp, "%9d.%09d\t%e\t%e\t%e\t%e\n",
          (int) this_event->start_time.gpsSeconds,
          (int) this_event->start_time.gpsNanoSeconds,
          this_event->snr,
          this_event->amplitude,
          this_event->central_freq,
          this_event->bandwidth );
      this_event = this_event->next;
    }
    fclose( fp );
  }

  vrbmsg( "cleaning up" );

  /* free event list */
  while ( events )
  {
    SnglBurstTable *next = events->next;
    LALFree( events );
    events = next;
  }

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
    { "strain-data", no_argument, &strain_data, 1 },
    { "geo-data", no_argument, &geo_data, 1 },
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
    { "calibration-cache",       required_argument, 0, 'C' },
    { "debug-level",             required_argument, 0, 'd' },
    { "data-frame-cache",        required_argument, 0, 'D' },
    { "frame-cache",             required_argument, 0, 'D' },
    { "frame-files",             required_argument, 0, 'f' },
    { "frame-path",              required_argument, 0, 'F' },
    { "bank-start-template",     required_argument, 0, 'i' },
    { "bank-end-template",       required_argument, 0, 'j' },
    { "inject-file",             required_argument, 0, 'I' },
    { "output-file",             required_argument, 0, 'o' },
    { "output-format",           required_argument, 0, 'O' },
    { "response-file",           required_argument, 0, 'r' },
    { "sample-rate",             required_argument, 0, 's' },
    { "user-tag",                required_argument, 0, 'U' },
    { "userTag",                 required_argument, 0, 'U' },
    /* these are filter parameters */
    { "filter-params",           required_argument, 0, 'z' },
    { "filter-segsz",            required_argument, 0, 'Z' },
    { "filter-speclen",          required_argument, 0, 'Z' },
    { "filter-flow",             required_argument, 0, 'Z' },
    { "filter-fhighpass",        required_argument, 0, 'Z' },
    { "filter-fmin",             required_argument, 0, 'Z' },
    { "filter-fmax",             required_argument, 0, 'Z' },
    { "filter-qmin",             required_argument, 0, 'Z' },
    { "filter-qmax",             required_argument, 0, 'Z' },
    { "filter-maxmm",            required_argument, 0, 'Z' },
    { "filter-thresh",           required_argument, 0, 'Z' },
    { "filter-scale",            required_argument, 0, 'Z' },
    /* these options are for geo data high pass filter */
    { "geo-high-pass-freq",      required_argument, 0, 'G' },
    { "geo-high-pass-order",     required_argument, 0, 'H' },
    { "geo-high-pass-atten",     required_argument, 0, 'P' },
    { "geo-dyn_range",           required_argument, 0, 'e' },
    /* these options are for writing output */
    { "write-format", required_argument, 0, 'w' },
    { "write-raw-data",         no_argument, &write_raw_data,         1 },
    { "write-data",             no_argument, &write_data,             1 },
    { "write-response",         no_argument, &write_response,         1 },
    { "write-inverse-spectrum", no_argument, &write_inverse_spectrum, 1 },
    { "write-data-segments",    no_argument, &write_data_segments,    1 },
    { "write-filter-output",    no_argument, &write_filter_output,    1 },
    /* these options are for internal tests */
    { "test-zero-data",         no_argument,    &test_zero_data,        1 },
    { "test-gaussian-data",     no_argument,    &test_gaussian_data,    1 },
    { "test-white-spectrum",    no_argument,    &test_white_spectrum,   1 },
    { "test-unit-response",     no_argument,    &test_unit_response,    1 },
    { "test-inject",            required_argument,      0,      'T' },
    /* nul terminate */
    { 0, 0, 0, 0 }
  };
  char args[] = "hVa:A:b:B:c:C:d:D:f:F:i:j:I:o:O:r:s:T:U:w:z:Z:G:H:P:e";

  program  = argv[0];
  fargv[0] = my_strdup( program );

  while ( 1 )
  {
    int option_index = 0; /* getopt_long stores long option here */
    int c;
    int n;

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
      case 'C': /* calibration-frame-cache or calibration-cache */
        frcalib = optarg;
        break;
      case 'd': /* debug-level */
        dbglvl = optarg;
        break;
      case 'D': /* data-frame-cache or frame-cache */
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
      case 'I': /* inject-file */
        injfile = optarg;
        break;
      case 'o': /* output-file */
        outfile = optarg;
        break;
      case 'O' : /* output-format */
        outfmt = strstr( optarg, "asc" ) ? ascii_format : xml_format;
        break;
      case 'r': /* response-file */
        rspfile = optarg;
        break;
      case 's': /* sample-rate */
        srate = atof( optarg );
        break;
      case 'T': /* test-inject */
        n = sscanf( optarg, "%d.%d,%f,%f,%f,%f",
            &test_inject_sec, &test_inject_nan, &test_inject_freq,
            &test_inject_qual, &test_inject_ampl, &test_inject_phase );
        if ( n != 6 )
        {
          fprintf( stderr, "invalid format for --test-inject\n" );
          fprintf( stderr, "required format is sec.nan,freq,qual,ampl,phase\n");
          exit( 1 );
        }
        test_inject = 1;
        break;
      case 'U': /* user-tag: do nothing with this! */
        break;
      case 'w' : /* write-format */
        write_format = strstr( optarg, "asc" ) ? ascii_format : frame_format;
        break;
      case 'z': /* filter-params: inport a filter params file */
        import_filter_params( optarg );
        add_process_params_prefix( fargc, fargv, "--filter" );
        break;
      case 'Z': /* filter-something: put this in fargc */
        fargv[fargc++] = my_strdup(
            strstr( long_options[option_index].name, "filter" ) + 6 );
        fargv[fargc++] = my_strdup( optarg );
        break;
     
      case 'G': /* geo high pass filter knee frequency  */
	geoHighPassFreq= atof(optarg);
	break;
                          
       case 'H': /* geo high pass filter order  */
	 geoHighPassOrder = atoi(optarg);
	 break;

       case 'P': /* geo high pass filter attenuation  */
	  geoHighPassAtten = atof(optarg);
	  break;
      case 'e': /* geo dynamic range  */
        geoDynamicRange = atof(optarg);
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
REAL4TimeSeries *get_data( UINT4 segsz, const char *ifo )
{
  LALStatus  status = blank_status;
  FrChanIn   chanin = blank_chanin;
  FrStream  *stream = NULL;
  REAL4TimeSeries *channel;
  REAL8TimeSeries *geoChannel;
  UINT4     npts, j;
  
  PassBandParamStruc geoHighpassParam;

  /* allocate memory */
  channel = my_calloc( 1, sizeof( *channel ) );
  if ( geo_data )
   {
     geoChannel = my_calloc(1, sizeof( *geoChannel ) );
   }

  if ( test_gaussian_data ) /* generate white gaussian noise */
  {
    RandomParams *rpar = NULL;
    if ( ! ( srate > 0 ) )
    {
      fprintf( stderr, "sample rate must be set (use --sample-rate)\n" );
      exit( 1 );
    }
    npts = duration * srate;
    LAL_CALL( LALSCreateVector( &status, &channel->data, npts ), &status );
    LAL_CALL( LALCreateRandomParams( &status, &rpar, 0 ), &status );
    LAL_CALL( LALNormalDeviates( &status, channel->data, rpar ), &status );
    LAL_CALL( LALDestroyRandomParams( &status, &rpar ), &status );
    strncpy( channel->name, frchan, sizeof( channel->name ) );
    channel->epoch       = tstart;
    channel->deltaT      = 1.0 / srate;
    channel->sampleUnits = lalADCCountUnit;
  }
  else /* read real data */
  {
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
      int mode = LAL_FR_VERBOSE_MODE; /* fails if bad time request/data gaps */
      vrbmsg( "get data from frame files %s/%s", frpath, frptrn );
      LAL_CALL( LALFrOpen( &status, &stream, frpath, frptrn ), &status );
      LAL_CALL( LALFrSetMode( &status, mode, stream ), &status );
    }
    if ( tstart.gpsSeconds || tstart.gpsNanoSeconds )
    {
      vrbmsg( "seek to gps time %d.%09d",
          tstart.gpsSeconds, tstart.gpsNanoSeconds );
      LAL_CALL( LALFrSeek( &status, &tstart, stream ), &status );
      if ( tend.gpsSeconds || tend.gpsNanoSeconds )
      {
        INT8 ns;
        ns  = tend.gpsSeconds - tstart.gpsSeconds;
        ns *= 1000000000;
        ns += tend.gpsNanoSeconds - tstart.gpsNanoSeconds;
        if ( ns > 0 )
          duration = 1e-9 * ns;
      }
    }

    /* get the data */
    vrbmsg( "read %f seconds of data", duration );

    if (geo_data)
     {
      /* call first time to get sample rate */
      LAL_CALL( LALFrGetREAL8TimeSeries( &status, geoChannel, &chanin, stream ),
        &status );
      /* compute how much data we need to get and allocate memory */
      npts = duration / geoChannel->deltaT;
      LAL_CALL( LALDCreateVector( &status, &geoChannel->data, npts ), &status );     LAL_CALL( LALSCreateVector( &status, &channel->data, npts ), &status );
      /* get the data */
      LAL_CALL( LALFrGetREAL8TimeSeries( &status, geoChannel, &chanin, stream ),
        &status );

      /* close the frame file */
      LAL_CALL( LALFrClose( &status, &stream ), &status );

      /* high pass the GEO data using the parameters specified on the cmd line */      
       
       geoHighpassParam.nMax = geoHighPassOrder;
       geoHighpassParam.f1 = -1.0;
       geoHighpassParam.f2 = (REAL8) geoHighPassFreq;
       geoHighpassParam.a1 = -1.0;
       geoHighpassParam.a2 = (REAL8)(1.0 - geoHighPassAtten);

       LAL_CALL( LALButterworthREAL8TimeSeries( &status, geoChannel, 
          &geoHighpassParam ), &status );
      
       /* cast the GEO data to REAL4 in the channel time series       */
       /* which already has the correct amount of memory allocated */
       for ( j = 0 ; j < npts ; ++j )
        {
         channel->data->data[j] = geoDynamicRange * (REAL4) (geoChannel->data->data[j]);
         }
       /* re-copy the data paramaters from the GEO channel to input data channel */
       LALSnprintf( channel->name, LALNameLength * sizeof(CHAR), "%s", geoChannel->name );
       channel->epoch          = geoChannel->epoch;
       channel->deltaT         = geoChannel->deltaT;
       channel->f0             = geoChannel->f0;
       channel->sampleUnits    = geoChannel->sampleUnits;       
       
       /* free the REAL8 GEO input data */
       LAL_CALL( LALDDestroyVector( &status, &geoChannel->data ), &status );
       geoChannel->data = NULL;
     }
 
    else
     {
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
     }
  }

  /* if this is strain data, correct the units */
  if ( strain_data )
    channel->sampleUnits = lalStrainUnit;

  /* write this data to file if required */
  if ( write_raw_data )
  {
    if ( write_format == frame_format )
    {
      LALSnprintf( channel->name, sizeof( channel->name ), "%s_RAW", frchan );
      vrbmsg( "writing raw data to frame" );
      write_frame = fr_add_proc_REAL4TimeSeries( write_frame, channel );
    }
    else
    {
      const char *fname = "ring-raw-data.dat";
      vrbmsg( "writing raw data to file %s", fname );
      LALSPrintTimeSeries( channel, fname );
    }
  }

  /* set epoch and duration to actual start time and duration */
  tstart   = channel->epoch;
  duration = channel->deltaT * channel->data->length;

  /* inject a signal if required */
  if ( injfile )
  {
    INT4 start = tstart.gpsSeconds;
    INT4 stop  = start + (INT4)ceil( 1e-9 * tstart.gpsNanoSeconds + duration );
    COMPLEX8FrequencySeries *response;
    SimBurstTable *inj = NULL;
    response = get_response( channel->data->length, channel->deltaT, ifo );
    vrbmsg( "reading simulated-burst tables from file %s", injfile );
    LAL_CALL( LALSimBurstTableFromLIGOLw( &status, &inj, injfile, start, stop ),
        &status );
    vrbmsg( "injecting signals into time series" );
    LAL_CALL( LALBurstInjectSignals( &status, channel, inj, response ),
        &status );
    while ( inj ) /* free injection list */
    {
      SimBurstTable *next = inj->next;
      LALFree( inj );
      inj = next;
    }
    LAL_CALL( LALCDestroyVector( &status, &response->data ), &status );
    free( response );
  }

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

  if ( strain_data ) /* strain data: response function should be unity */
  {
    UINT4 i;
    for ( i = 0; i < response->data->length; ++i ) /* set response to unity */
    {
      response->data->data[i].re = 1;
      response->data->data[i].im = 0;
    }
    /* response is dimensionless */
    memset( &response->sampleUnits, 0, sizeof( response->sampleUnits ) );
  }
  else if ( frcalib ) /* use calibration cache to extract & update response */
  {
    vrbmsg( "get calibration data from frames in cache file %s", frcalib );
    memset( &calfacts, 0, sizeof(CalibrationUpdateParams) );
    calfacts.ifo = ifo;
    LAL_CALL( LALExtractFrameResponse( &status, response, frcalib, &calfacts ),
        &status );
  }
  else /* get fixed response function from an ascii file */
  {
    vrbmsg( "get calibration data from ascii file %s", rspfile );
    read_response( response, rspfile );
  }

  /* write this data to file if required */
  if ( write_response )
  {
    LALSnprintf( response->name, sizeof( response->name ), "%s_RESPONSE",
        frchan );
    if ( write_format == frame_format )
    {
      vrbmsg( "writing response to frame" );
      write_frame = fr_add_proc_COMPLEX8FrequencySeries( write_frame,
          response );
    }
    else
    {
      const char *fname = "ring-response.dat";
      vrbmsg( "writing response to file %s", fname );
      LALCPrintFrequencySeries( response, fname );
    }
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

struct series blank_series;

FrameH *fr_add_proc_REAL4TimeSeries( FrameH *frame, REAL4TimeSeries *series )
{
  LALStatus status = blank_status;
  struct series data = blank_series;
  CHARVector uvec;
  char unit[256];
  uvec.length = sizeof( unit );
  uvec.data   = unit;
  LAL_CALL( LALUnitAsString( &status, &uvec, &series->sampleUnits ), &status );
  data.name = series->name;
  data.tbeg = series->epoch;
  epoch_add( &data.tend, &series->epoch, series->deltaT*series->data->length );
  data.dom  = Time;
  data.type = FR_VECT_4R;
  data.step = series->deltaT;
  data.unit = unit;
  data.size = series->data->length;
  data.data = (float *)series->data->data;
  return fr_add_proc_data( frame, &data );
}

FrameH *fr_add_proc_REAL4FrequencySeries( FrameH *frame,
    REAL4FrequencySeries *series )
{
  LALStatus status = blank_status;
  struct series data = blank_series;
  CHARVector uvec;
  char unit[256];
  uvec.length = sizeof( unit );
  uvec.data   = unit;
  LAL_CALL( LALUnitAsString( &status, &uvec, &series->sampleUnits ), &status );
  data.name = series->name;
  data.tbeg = series->epoch;
  epoch_add( &data.tend, &series->epoch,
      ( series->data->length - 1 )/( series->deltaF * series->data->length ) );
  data.dom  = Freq;
  data.type = FR_VECT_4R;
  data.step = series->deltaF;
  data.unit = unit;
  data.size = series->data->length;
  data.data = (float *)series->data->data;
  return fr_add_proc_data( frame, &data );
}

FrameH *fr_add_proc_COMPLEX8FrequencySeries( FrameH *frame,
    COMPLEX8FrequencySeries *series )
{
  LALStatus status = blank_status;
  struct series data = blank_series;
  CHARVector uvec;
  char unit[256];
  uvec.length = sizeof( unit );
  uvec.data   = unit;
  LAL_CALL( LALUnitAsString( &status, &uvec, &series->sampleUnits ), &status );
  data.name = series->name;
  data.tbeg = series->epoch;
  epoch_add( &data.tend, &series->epoch,
      ( series->data->length - 1 )/( series->deltaF * series->data->length ) );
  data.dom  = Freq;
  data.type = FR_VECT_8C;
  data.step = series->deltaF;
  data.unit = unit;
  data.size = series->data->length;
  data.data = (float *)series->data->data;
  return fr_add_proc_data( frame, &data );
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
"  --calibration-frame-cache calcache\n"
"  --calibration-cache       calcache\n\t\tcalib frame cache file [use ascii]\n\n"
"  --response-file respfile\n\t\tascii response file [response.asc]\n\n"
"  --data-frame-cache datcache\n"
"  --frame-cache      datcache\n\t\tdata frame cache file [don't use cache]\n\n"
"  --frame-files pattern\n\t\tframe file to use [*.gwf]\n\n"
"  --frame-path path\n\t\tpath to look for frame files [.]\n\n"
"  --geo-data\n\t\tdata is double precision\n\n"
"  --strain-data\n\t\tdata is strain data (use unit response)\n\n"
"  --bank-start-template bmin\n\t\tfirst template of bank to use [0]\n\n"
"  --bank-end-tempate bmax\n\t\tlast template of bank to use [last]\n\n"
"  --inject-file injfile\n\t\tLIGOLw XML file containing injection parameters\n\n"
"  --output-file outfile\n\t\tevent output file [IFO-RING-tstart-duration.xml]\n\n"
"  --output-format format\n\t\tevent output file format (xml or ascii) [xml]\n\n"
"  --sample-rate srate\n\t\tdesired sampling rate in Hz [raw rate]\n\n"
"  --user-tag tag\n"
"  --userTag  tag\n\t\tuser tag [none]\n\n"
"Write Options (used to write intermediate results to files)\n\n"
"  --write-format format\n\t\tformat to write results (frame or ascii) [frame]\n\n"
"  --write-raw-data\n\t\twrite data input before pre-processing\n\n"
"  --write-data\n\t\twrite data after pre-processing (and injections)\n\n"
"  --write-response\n\t\twrite response function used\n\n"
"  --write-inverse-spectrum\n\t\twrite inverse spectrum computed\n\n"
"  --write-data-segments\n\t\twrite conditioned data segments\n\n"
"  --write-filter-output\n\t\twrite matched filter output\n\n"
"Internal Test Options\n\n"
"  --test-zero-data\n\t\tzero the data before performing the search\n\n"
"  --test-gaussian-data\n\t\tgenerate white Gaussian data\n\n"
"  --test-white-spectrum\n\t\tuse an exact white spectrum for Gaussian data\n\n"
"  --test-unit-response\n\t\tuse a constant unit response function\n\n"
"  --test-inject sec.nan,freq,qual,ampl,phase\n\t\tperform internal test injection\n\n"
"Filter Parameters [defaults in brackets; otherwise required]\n\n"
"  --filter-params parfile\n\t\tread filter parameters from file [none]\n\n"
"  --filter-segsz npts\n\t\tset size of segments to analyze to npts\n\n"
"  --filter-speclen len\n\t\tset size of inverse spectrum truncation to len [0]\n\n"
"  --filter-flow flow\n\t\tset low frequency cutoff to flow (Hz)\n\n"
"  --filter-fhighpass fhighpass\n\t\tset highpass frequency to fhighpass (Hz)\n\n"
"  --filter-fmin fmin\n\t\tset minimum frequency for bank to fmin (Hz)\n\n"
"  --filter-fmax fmax\n\t\tset maximum frequency for bank to fmax (Hz)\n\n"
"  --filter-qmin qmin\n\t\tset minimum quality for bank to qmin\n\n"
"  --filter-qmax qmax\n\t\tset maximum quality for bank to qmax\n\n"
"  --filter-maxmm maxmm\n\t\tset maximum mismatch of bank to maxmm\n\n"
"  --filter-thresh thrsh\n\t\tset ringdown event snr threshold to thrsh\n\n"
"  --filter-scale scale\n\t\tscale response function by a factor of scale [1]\n\n"
"  --geo-high-pass-freq scale\n\t\tknee frequency for geo high pass filter[1]\n\n"
"  --geo-high-pass-order scale\n\t\torder of geo high pass filter[1]\n\n"
"  --geo-high-pass-atten scale\n\t\tattenuation of geo high pass filter[1]\n\n";


int usage( const char *program )
{
  return fprintf( stderr, usgfmt, program );
}

#endif
