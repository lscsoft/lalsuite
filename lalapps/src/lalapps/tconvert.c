#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define _GNU_SOURCE
#include "getopt.h"

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include "getdate.h"
#include "lalapps.h"

#define MAX_DATE_STRING_LENGTH 256

const char * skip_space( const char *s )
{
  while ( s && *s && isspace( *s ) )
    ++s;
  return s;
}

static int is_integer( const char *s )
{
  if ( ! s )
    return -1;
  s = skip_space( s );
  if ( ! *s ) /* it was blank, not an integer! */
    return 0;
  while ( *s && isdigit( *s ) )
    ++s;
  s = skip_space( s );
  return *s ? 0 : 1;
}

static int is_blank( const char *s )
{
  s = skip_space( s );
  return *s ? 0 : 1;
}

const char *rcsid = "$Id$";

void output( int gps_sec, int output_type );
void usage( const char *program, int exitcode );
char * parse_options( char *buf, int buflen, int argc, char **argv );
int    unix_to_gps( time_t unix_sec );
time_t gps_to_unix( int gps_sec, int *leap );

/* standard formats */
const char *iso_8601_format = "%Y-%m-%dT%H:%M:%S%z";
const char *rfc_822_format = "%a, %d %b %Y %H:%M:%S %z";
const char *unix_date_format = "%a %b %d %H:%M:%S %Z %Y";
const char *output_date_format = NULL;

enum { GPS_EPOCH, UNIX_EPOCH } epoch = GPS_EPOCH; /* default is GPS epoch */
enum { SIDEREAL_HMS, SIDEREAL_RAD } sidereal_format = SIDEREAL_HMS; /* default isHH:MM:SS */
enum { OUTPUT_DEFAULT, OUTPUT_DATE, OUTPUT_SECS, OUTPUT_JD, OUTPUT_MJD, OUTPUT_LEAPS, OUTPUT_GMST } 
    output_type = OUTPUT_DEFAULT; 
int utc_flag = 0; /* default is local time */

int verbose = 0; /* default is brief output */

FILE *fp = NULL; /* set to stdin to read dates from stdin */

int lalDebugLevel = 7; /* print info, warnings, and errors */

int main( int argc, char *argv[] )
{
  char buf[1024] = "";
  char * arg;  /* argument after options have been parsed */
  time_t unix_sec;
  int gps_sec;
  int leap;

  XLALSetErrorHandler( XLALExitErrorHandler );

  errno = 0;
  arg = parse_options( buf, sizeof( buf ) - 1, argc, argv );
  if ( ! arg )
  {
    if ( errno )
      perror( "parse options" );
    usage( argv[0], 1 );
  }
  
  while ( arg )
  {
    if ( is_integer( arg ) ) /* argument is an integer seconds since epoch */
    {
      if ( output_type == OUTPUT_DEFAULT && ! verbose ) 
        output_type = OUTPUT_DATE;
      switch ( epoch )
      {
        case GPS_EPOCH:
          gps_sec = atoi( arg );
          break;
        case UNIX_EPOCH:
          unix_sec = atoi( arg );
          gps_sec = unix_to_gps( unix_sec );
          break;
        default:
          fprintf( stderr, "internal error: unrecognized epoch\n" );
          exit( 1 );
      }
    }
    else /* assume argument is a date string */
    {
      if ( output_type == OUTPUT_DEFAULT && ! verbose ) 
        output_type = OUTPUT_SECS;
      if ( is_blank( arg ) )
        unix_sec = get_date( "now", 0 );
      else
        unix_sec = get_date( arg, 0 );
      /* need to check to see if there was a leap second */
    leap = 0;
    if ( strstr( arg, ":59:60" ) ) /* aha! there was one! */
    {
      --unix_sec; /* to distinguish this from next day */
      leap = 1; /* add the leap second back in later */
    }

    if ( unix_sec < 0 )
    {
      fprintf( stderr, "could not parse date string: %s\n", arg );
      exit( 1 );
    }
    gps_sec = unix_to_gps( unix_sec );
    gps_sec += leap; /* the leap second */
    }

    output( gps_sec, output_type );

    if ( fp )
      arg = fgets( buf, sizeof( buf ) - 1, fp  );
    else
      break;
  }

  return 0;
}


void output_secs( int gps_sec )
{
  time_t unix_sec;
  int leap;
  switch ( epoch )
  {
    case GPS_EPOCH:
      if ( verbose )
        printf( "GPS Time = " );
      printf( "%d\n", gps_sec );
      break;
    case UNIX_EPOCH:
      if ( verbose )
        printf( "UNIX Time = " );
      unix_sec = gps_to_unix( gps_sec, &leap );
      printf( "%d", (int)unix_sec );
      if ( leap )
      {
        if ( verbose )
          printf( " + 1 leap second" );
        else
          printf( "+1" );
      }
      printf( "\n" );
      break;
    default:
      fprintf( stderr, "internal error: unrecognized epoch\n" );
      exit( 1 );
  }
  return;
}

void output_date( int gps_sec )
{
  char date_string[256];
  struct tm *tp;
  time_t unix_sec;
  const char *fmt;
  int leap;

  fmt = output_date_format ? output_date_format : unix_date_format;
  unix_sec = gps_to_unix( gps_sec, &leap );
  if ( utc_flag )
    tp = gmtime( &unix_sec );
  else
    tp = localtime( &unix_sec );
  if ( leap )
    ++tp->tm_sec;

  strftime( date_string, sizeof( date_string ), fmt, tp );
  puts( date_string );
  return;
}

void output_jd( int gps_sec )
{
  struct tm utc;
  double jd;
  utc = *XLALGPSToUTC( &utc, gps_sec );
  jd = XLALJulianDay( &utc );
  if ( verbose )
    printf( "Julian Day = " );
  printf( "%.6f\n", jd );
  return;
}

void output_mjd( int gps_sec )
{
  struct tm utc;
  double jd;
  double mjd;
  utc = *XLALGPSToUTC( &utc, gps_sec );
  jd = XLALJulianDay( &utc );
  mjd = jd - XLAL_MJD_REF;
  if ( verbose )
    printf( "Modified Julian Day = " );
  printf( "%.6f\n", mjd );
  return;
}

void output_leaps( int gps_sec )
{
  int leaps;
  leaps = XLALLeapSeconds( gps_sec );
  if ( verbose )
    printf( "Leap Seconds (TAI-UTC) = " );
  printf( "%d\n", leaps );
  return;
}

void output_gmst( int gps_sec )
{
  LIGOTimeGPS ligo_time;
  double gmst_rad;
  double gmst_hrs;
  double gmst_min;
  double gmst_sec;
  ligo_time.gpsSeconds = gps_sec;
  ligo_time.gpsNanoSeconds = 0;
  gmst_rad = XLALGreenwichMeanSiderealTime( &ligo_time );
  while ( gmst_rad >= 2.0 * LAL_PI )
    gmst_rad -= 2.0 * LAL_PI;
  while ( gmst_rad < 0 )
    gmst_rad += 2.0 * LAL_PI;
  switch ( sidereal_format )
  {
    case SIDEREAL_HMS:
      gmst_hrs = 12.0 * gmst_rad / LAL_PI;
      gmst_min = 60.0 * modf( gmst_hrs, &gmst_hrs );
      gmst_sec = 60.0 * modf( gmst_min, &gmst_min );
      gmst_sec = floor( gmst_sec + 0.5 ); /* round to nearest second */
      if ( verbose )
        printf( "Greenwich Mean Sidereal Time (HH:MM:SS) = " );
      printf( "%02.0f:%02.0f:%02d\n", gmst_hrs, gmst_min, (int)gmst_sec );
      break;
    case SIDEREAL_RAD:
      if ( verbose )
        printf( "Greenwich Mean Sidereal Time (radians) = " );
      printf( "%.6f\n", gmst_rad );
      break;
    default:
      fprintf( stderr, "internal error: unrecognized sidereal time format\n" );
      exit( 1 );
  }
  return;
}

void output( int gps_sec, int type )
{
  switch ( type )
  {
    case OUTPUT_DATE:
      output_date( gps_sec );
      break;
    case OUTPUT_SECS:
      output_secs( gps_sec );
      break;
    case OUTPUT_JD:
      output_jd( gps_sec );
      break;
    case OUTPUT_MJD:
      output_mjd( gps_sec );
      break;
    case OUTPUT_LEAPS:
      output_leaps( gps_sec );
      break;
    case OUTPUT_GMST:
      output_gmst( gps_sec );
      break;
    case OUTPUT_DEFAULT:
      if ( verbose )
      {
        int save_epoch = epoch;
        output_date( gps_sec );
        epoch = GPS_EPOCH;
        output_secs( gps_sec );
        epoch = UNIX_EPOCH;
        output_secs( gps_sec );
        epoch = save_epoch;
        output_leaps( gps_sec );
        output_jd( gps_sec );
        output_mjd( gps_sec );
        output_gmst( gps_sec );
        break;
      }
      /* fall-through */
    default:
      fprintf( stderr, "internal error: unrecognized output type\n" );
      exit( 1 );
  }
  return;
}

/*
 * If this is a leap second, return the unix time of the previous second
 * and set *leap to 1.  Thus if you add the return value and *leap you'll
 * have the actual unix time.  This is done since the two times, say
 * 16:59:60 of some day and 17:00:00 of the same day have the same unix
 * time.  If this routine gets a gps_sec corresponding to 16:59:60, it will
 * return a unix_sec corresponding to 16:59:59 and will set *leap to 1.
 * However if it gets a gps time corresponding to 17:00:00, it will return
 * a unix_sec corresponding to 17:00:00 and will set *leap to 0.
 */
time_t gps_to_unix( int gps_sec, int *leap )
{
  struct tm utc;
  time_t unix_sec;

  utc = *XLALGPSToUTC( &utc, gps_sec );
  if ( leap ) /* worry about whether this second was a leap second */
  {
    *leap = 0;
    if ( utc.tm_sec == 60 )
    {
      --utc.tm_sec;
      *leap = 1;
    }
  }

  /* posix definition of the unix time */
  unix_sec = utc.tm_sec + utc.tm_min*60 + utc.tm_hour*3600 
    + utc.tm_yday*86400 + (utc.tm_year-70)*31536000
    + ((utc.tm_year-69)/4)*86400 - ((utc.tm_year-1)/100)*86400
    + ((utc.tm_year+299)/400)*86400;

  return unix_sec;
}

int unix_to_gps( time_t unix_sec )
{
  struct tm *tp;
  int gps_sec;

  tp = gmtime( &unix_sec );
  gps_sec = XLALUTCToGPS( tp );

  return gps_sec;
}

void usage( const char *program, int exitcode )
{
  const char *usages = "Usage:\t%s [OPTIONS] DATE_STRING\n"
                       "  or: \t%s [OPTIONS] SEC_SINCE_EPOCH\n";
  const char *synopsis = "Convert between various representations of time\n";
  const char *options = "Options:\n\n\
        -fFORMAT, --format=FORMAT\n\
                output format for date string [see date(1)]\n\
        \n\
        -G, --gps-epoch\n\
                use GPS epoch (1980-01-06 00:00:00 UTC)\n\
        \n\
        -I, --iso-8601\n\
                output date string in ISO 8601 format\n\
        \n\
        -j, --jd, --julian-day\n\
                output Julian day\n\
        \n\
        -l, --leap-seconds\n\
                output number of leap seconds (TAI-UTC)\n\
        \n\
        -m, --mjd, --modified-julian-day\n\
                output modified Julian day\n\
        \n\
        -R, --rfc-2822\n\
                output RFC-2822 compliant date string\n\
        \n\
        -s, --sidereal-time[= HMS|RAD]\n\
                print Greenwich mean sidreal time [default format HH:MM:SS]\n\
        \n\
        -S, --stdin\n\
                read from standard input\n\
        \n\
        -u, --utc\n\
                print Coordinated Universal Time\n\
        \n\
        -U, --unix-epoch\n\
                use Unix epoch (1970-01-01 00:00:00 UTC)\n\
        \n\
        --help  display this help and exit\n\
        \n\
        --version\n\
                output version information and exit\n";
  fprintf( stdout, usages, program, program );
  fprintf( stdout, "%s", synopsis );
  fprintf( stdout, "%s", options );
  exit( exitcode );
}


char * parse_options( char *buf, int buflen, int argc, char **argv )
{
  struct option long_options[] =
  {
    { "help",           no_argument,            0,      'h'     },
    { "utc",            no_argument,            0,      'u'     },
    { "universal",      no_argument,            0,      'u'     },
    { "unix-epoch",     no_argument,            0,      'U'     },
    { "gps-epoch",      no_argument,            0,      'G'     },
    { "jd",             no_argument,            0,      'j'     },
    { "julian-day",     no_argument,            0,      'j'     },
    { "mjd",            no_argument,            0,      'm'     },
    { "modified-julian-day",  no_argument,      0,      'm'     },
    { "leap-seconds",   no_argument,            0,      'l'     },
    { "format",         required_argument,      0,      'f'     },
    { "iso-8601",       no_argument,            0,      'I'     },
    { "rfc-8286",       no_argument,            0,      'R'     },
    { "sidereal-time",  optional_argument,      0,      's'     },
    { "stdin",          no_argument,            0,      'S'     },
    { "verbose",        no_argument,            0,      'v'     },
    { "version",        no_argument,            0,      'V'     },
    { 0, 0, 0, 0 }
  };
  char args[] = "+f:GhIjlmRsSuUvVW:";
  char *program = argv[0];

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
        if ( ! long_options[option_index].flag )
        {
          fprintf( stderr, "error parsing option %s with argument %s\n",
              long_options[option_index].name, optarg );
          usage( program, 1 );
        }
        break;

      case 'f': /* format */
        output_date_format = optarg;
        break;

      case 'G': /* gps-epoch */
        epoch = GPS_EPOCH;
        break;

      case 'I': /* iso-8601 */
        output_date_format = iso_8601_format;
        break;

      case 'j': /* julian-day */
        output_type = OUTPUT_JD;
        break;

      case 'l': /* leap-seconds */
        output_type = OUTPUT_LEAPS;
        break;

      case 'm': /* mjd */
        output_type = OUTPUT_MJD;
        break;

      case 'R': /* rfc-8286 */
        output_date_format = rfc_822_format;
        break;

      case 's': /* sidereal-time */
        output_type = OUTPUT_GMST;
        if ( optarg )
        {
          if ( ! strcmp( optarg, "rad" ) || ! strcmp( optarg, "RAD" ) )
            sidereal_format = SIDEREAL_RAD;
          else if ( ! strcmp( optarg, "hms" ) || ! strcmp( optarg, "HMS" ) )
            sidereal_format = SIDEREAL_HMS;
          else
          {
            fprintf( stderr, "unrecognized sidereal time format %s\n", optarg );
            fprintf( stderr, "expect \"HMS\" (hours:minutes:seconds)" );
            fprintf( stderr, "or \"RAD\" (radians)\n" );
            exit( 1 );
          }
        }
        break;

      case 'S': /* stdin */
        fp = stdin;
        break;

      case 'u': /* utc */
        utc_flag = 1;
        break;

      case 'U': /* unix-epoch */
        epoch = UNIX_EPOCH;
        break;

      case 'v': /* verbose */
        verbose = 1;
        break;

      case 'V': /* version */
        PRINT_VERSION( "tconvert" );
        exit( 0 );
        break;

      case '?':
        /* fall-through */
      default:
        fprintf( stderr, "unknown error while parsing options\n" );
        usage( program, 1 );
    }
  }

  if ( optind == argc && fp )
      fgets( buf, buflen, fp );
  else
  {
    for ( ; optind < argc; ++optind )
    {
      int len;
      len = strlen( buf );
      strncat( buf + len, argv[optind], buflen - len );
      len = strlen( buf );
      if ( len == buflen )
      {
        fprintf( stderr, "error: line too long\n" );
        exit( 1 );
      }
      strncat( buf + len, " ", buflen - len );
      if ( len == buflen )
      {
        fprintf( stderr, "error: line too long\n" );
        exit( 1 );
      }
    }
  }

  return buf;
}
