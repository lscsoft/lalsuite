#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <lalapps.h>
#include <lal/LALStdlib.h>
#include <lal/Units.h>
#include <lal/AVFactories.h>
#include <lal/FrameStream.h>
#include <lal/RingSearch.h>

#include <lal/PrintFTSeries.h>


RCSID( "$Id$" );

#define usgfmt \
  "Usage: %s [options] -- [filterparams]\n" \
  "Options [default in brackets]:\n" \
  "  -h            print this message\n" \
  "  -V            print version info\n" \
  "  -v            verbose\n" \
  "  -d dbglvl     set debug level to dbglvl [0]\n" \
  "  -i infile     use input file infile for filter params [stdin]\n" \
  "  -o outfile    use output file outfile [stdout]\n" \
  "  -k            keep filtering output (use with -s)\n" \
  "  -s            save filtering data as .dat and snr- files\n" \
  "  -f framefile  frame file names [./*.gwf]\n" \
  "  -r respfile   response file name [response.asc]\n" \
  "  -c channel    channel name [H1:LSC-AS_Q]\n" \
  "  -n numpoints  number of points of data to analyze [65536]\n" \
  "  -t starttime  use data beginning at gps time starttime\n" \
  "  -b b0[,b1]    filter only templates numbers b0 to b1 [0,end-of-bank]\n" \
  "Filter Parameters [defaults in brackets; otherwise required]\n" \
  "  -segsz npts   set size of segments to analyze to npts\n" \
  "  -speclen len  set size of inverse spectrum truncation to len [0]\n" \
  "  -flow flow    set low frequency cutoff to flow (Hz)\n" \
  "  -fmin fmin    set minimum frequency for bank to fmin (Hz)\n" \
  "  -fmax fmax    set maximum frequency for bank to fmax (Hz)\n" \
  "  -qmin qmin    set minimum quality for bank to qmin\n" \
  "  -qmax qmax    set maximum quality for bank to qmax\n" \
  "  -maxmm maxmm  set maximum mismatch of bank to maxmm\n" \
  "  -thresh thrsh set ringdown event snr threshold to thrsh\n" \
  "  -scale scale  scale response function by a factor of scale [1]\n"

#define usage( program ) fprintf( stderr, usgfmt, program )

#ifndef LAL_FRAME_ENABLED
int main( void ) { fputs( "LAL not frame-enabled\n", stderr ); return 77; }
#else

extern char *optarg;
extern int optind, opterr, optopt;
extern int vrbflg;

int read_response( COMPLEX8FrequencySeries *series, const char *fname );

int main( int argc, char *argv[] )
{
  enum { maxargs = 256 };

  const char *program = argv[0];
  const char *infile  = NULL;
  const char *outfile = NULL;
  const char *dbglvl  = NULL;
  const char *frpath  = ".";
  const char *frptrn  = "*.gwf";
  const char *frchan  = "H1:LSC-AS_Q";
  const char *rspfile = "response.asc";

  LALStatus status = blank_status;

  static RingSearchData    data;
  static RingSearchInput   input;
  static RingSearchParams *params;
  static RingEventList    *events;

  static FrChanIn          chanin;
  static FrStream         *stream;
  static LIGOTimeGPS       epoch;

  RAT4 minusOne = { -1, 0 };
  LALUnitPair unitPair;
  LALUnit     unit;

  size_t npts = 65536;
  size_t rslt;

  const char *fltrargv[maxargs];
  int fltrargc = 0;

  char *tmpstr;
  int bmin = -1;
  int bmax = -1;
  int keep = 0;
  int save = 0;
  int opt;

  FILE *fp;

  /* parse options */
  while ( 0 < ( opt = getopt( argc, argv, "hVvd:i:o:ksf:r:n:t:b:" ) ) )
  {
    switch ( opt )
    {
      case 'h':
        usage( program );
        return 0;
      case 'V':
        PRINT_VERSION( "ring" );
        return 0;
      case 'v':
        vrbflg = 1;
        break;
      case 'd':
        dbglvl = optarg;
        break;
      case 'i':
        infile = optarg;
        break;
      case 'o':
        outfile = optarg;
        break;
      case 'k':
        keep = 1;
        break;
      case 's':
        save = 1;
        break;
      case 'f':
        if ( ( tmpstr = strrchr( optarg, '/' ) ) )
        {
          *tmpstr++ = 0;
          frptrn = tmpstr;
          frpath = optarg;
        }
        else
          frptrn = optarg;
        break;
      case 'r':
        rspfile = optarg;
        break;
      case 'n':
        npts = atoi( optarg );
        break;
      case 't':
        if ( ( tmpstr = strchr( optarg, '.' ) ) )
        {
          epoch.gpsNanoSeconds = 1e9 * atof( tmpstr );
          *tmpstr = 0;
        }
        epoch.gpsSeconds = atoi( optarg );
        break;
      case 'b':
        if ( ( tmpstr = strchr( optarg, ',' ) ) )
        {
          *tmpstr++ = 0;
          if ( strlen( tmpstr ) )
            bmax = atoi( tmpstr );
        }
        if ( strlen( optarg ) )
          bmin = atoi( optarg );
        break;
      default:
        usage( program );
        return 1;
    }
  }

  /* remaining arguments are filter parameters */
  --optind;
  argc -= optind;
  argv += optind;

  if ( argc > 1 ) /* copy args to fltrargv */
  {
    int arg;
    fltrargc = argc;
    for ( arg = 0; arg < argc; ++arg )
      fltrargv[arg] = argv[arg];
    fltrargv[fltrargc] = NULL;
  }
  else /* no remaining arguments: read from stdin */
  {
    char  line[1024];
    FILE *fpin;
    fpin = infile ? fopen( infile, "r" ) : stdin;
    fltrargc = 0;
    fltrargv[fltrargc++] = infile ? infile : "<stdin>";
    while ( fgets( line, sizeof( line ), fpin ) )
    {
      int n00, n01;
      int n10, n11;
      char *a0;
      char *a1;
      if ( *line == '#' )
        continue;
      if ( EOF == sscanf( line, "%n%*s%n %n%*s%n", &n00, &n01, &n10, &n11 ) )
      {
        fprintf( stderr, "syntax error on input line: %s\n", line );
        continue;
      }
      a0 = malloc( n01 - n00 + 1 ); 
      a1 = malloc( n11 - n10 + 1 ); 
      if ( 0 < sscanf( line, "%s %s", a0, a1 ) )
      {
        fltrargv[fltrargc++] = a0;
        fltrargv[fltrargc++] = a1;
      }
      else
      {
        free( a0 );
        free( a1 );
        fprintf( stderr, "syntax error on input line: %s\n", line );
      }
    }
    fltrargv[fltrargc] = NULL;
    outfile ? fclose( fpin ) : 0;
  }

  /* set debug level */
  set_debug_level( dbglvl );

  /* initialize ring search */
  LAL_CALL( LALRingSearchInit( &status, &params, fltrargv, fltrargc ),
      &status );

  /* allocate memory for data */
  data.channel  = LALCalloc( 1, sizeof( *data.channel ) );
  data.response = LALCalloc( 1, sizeof( *data.response ) );
  data.spectrum = NULL;
  LAL_CALL( LALSCreateVector( &status, &data.channel->data, npts ), &status );

  /* get frame data */
  chanin.name = frchan;
  chanin.type = ADCDataChannel;
  LAL_CALL( LALFrOpen( &status, &stream, frpath, frptrn ), &status );
  if ( epoch.gpsSeconds || epoch.gpsNanoSeconds )
    LAL_CALL( LALFrSeek( &status, &epoch, stream ), &status );
  LAL_CALL( LALFrGetREAL4TimeSeries( &status, data.channel, &chanin, stream ),
        &status );
  if ( save )
    LALSPrintTimeSeries( data.channel, "channel.dat" );
  LAL_CALL( LALFrClose( &status, &stream ), &status );

  /* get response */
  LAL_CALL( LALCCreateVector( &status, &data.response->data,
        params->segmentSize / 2 + 1 ), &status );
  data.response->f0     = 0.0;
  data.response->deltaF = 1.0 / ( params->segmentSize * data.channel->deltaT );
  LAL_CALL( LALUnitRaise( &status, &unit, &lalADCCountUnit, &minusOne ),
      &status );
  unitPair.unitOne = &lalStrainUnit;
  unitPair.unitTwo = &unit;
  LAL_CALL( LALUnitMultiply( &status, &data.response->sampleUnits, &unitPair ),
      &status );
  strncpy( data.response->name, "response", sizeof( data.response->name ) );
  read_response( data.response, rspfile );
  if ( save )
    LALCPrintFrequencySeries( data.response, "response.dat" );

  /* condition data */
  LAL_CALL( LALRingSearchConditionData( &status, params, &data ), &status );
  if ( save )
  {
    UINT4 seg;
    LALSPrintFrequencySeries( params->invSpectrum, "invspec.dat" );
    for ( seg = 0; seg < params->numSegments; ++seg )
    {
      char fname[64];
      sprintf( fname, "segment-%03d.dat", seg );
      LALCPrintFrequencySeries( params->dataSegment + seg, fname );
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
  params->keepResults = keep;
  LAL_CALL( LALRingSearch( &status, &events, &input, params ), &status );

  /* output events */
  fp = outfile ? fopen( outfile, "w" ) : stdout;
  if ( events )
    fprintf( fp, "# gps start time\tsignal/noise\tamplitude\tfrequency\tquality\n" );
  while ( events )
  {
    RingEventList *next = events->next;
    fprintf( fp, "%9d.%09d\t%e\t%e\t%e\t%e\n",
        (int)( events->startTimeNS / 1000000000 ),
        (int)( events->startTimeNS % 1000000000 ),
        events->snr, events->amplitude, events->frequency, events->quality );
    LALFree( events );
    events = next;
  }
  outfile ? fclose( fp ) : 0;

  if ( save )
    for ( rslt = 0; rslt < params->numResults; ++rslt )
    {
      LALSPrintTimeSeries( params->result + rslt, params->result[rslt].name );
    }

  /* deallocate response and channel data vectors */
  LAL_CALL( LALCDestroyVector( &status, &data.response->data ), &status );
  LAL_CALL( LALSDestroyVector( &status, &data.channel->data ), &status );

  /* finalize ring search */
  LAL_CALL( LALRingSearchFini( &status, &params ), &status );

  LALCheckMemoryLeaks();
  return 0;
}


int read_response( COMPLEX8FrequencySeries *series, const char *fname )
{
  char line[1024];
  FILE *fp;
  size_t n;
  size_t i, j;
  double *freq;
  /* double *rad; */
  /* double *phi; */
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
  /*
  rad  = LALCalloc( n, sizeof( *rad ) );
  phi  = LALCalloc( n, sizeof( *phi ) );
  */

  for ( i = 0; i < n; ++i )
  {
    /* double x, y, r; */
    fgets( line, sizeof( line ), fp );
    if ( *line == '#' )
      continue;
    sscanf( line, "%le %le %le %*e\n", freq + i, re + i, im + i );
    /*
    sscanf( line, "%le %le %le %*e\n", freq + i, &x, &y );
    rad[i] = r = sqrt( x * x + y * y );
    phi[i] = atan2( y, x );
    */
  }

  i = 1;
  for ( j = 0; j < series->data->length; ++j )
  {
    double f, x /*, r, q */;
    f = j * series->deltaF;
    while ( i < n )
      if ( freq[i] > f )
        break;
      else
        ++i;
    x = ( f - freq[i - 1] ) / ( freq[i] - freq[i - 1] );
    /*
    r = rad[i - 1] + x * ( rad[i] - rad[i - 1] );
    q = phi[i - 1] + x * ( phi[i] - phi[i - 1] );
    series->data->data[j].re = r * cos( q );
    series->data->data[j].im = r * sin( q );
    */
    series->data->data[j].re = re[i - 1] + x * ( re[i] - re[i - 1] );
    series->data->data[j].im = im[i - 1] + x * ( im[i] - im[i - 1] );
  }

  LALFree( freq );
  /*
  LALFree( rad );
  LALFree( phi );
  */
  LALFree( re );
  LALFree( im );
  
  fclose( fp );
  return 0;
}

#endif
