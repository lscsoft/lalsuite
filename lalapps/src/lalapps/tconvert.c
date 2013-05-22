/*
*  Copyright (C) 2007 Jolien Creighton
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

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define _GNU_SOURCE
#include "getopt.h"

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include "getdate.h"
#include "lalapps.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#define MAX_DATE_STRING_LENGTH 256

/* Longitudes of LHO and LLO taken from LIGO-T980044-08 */
#define LHO_LONGITUDE_RAD_E (LAL_PI*(-119.0+(24.0+27.5657/60.0)/60.0)/180.0)
#define LLO_LONGITUDE_RAD_E (LAL_PI*(-90.0+(46.0+27.2654/60.0)/60.0)/180.0)

static const char * skip_space( const char *s )
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

void output( int gps_sec, int output_type );
void usage( const char *program, int exitcode );
char * parse_options( char *buf, int buflen, int argc, char **argv );
int    unix_to_gps( time_t unix_sec );
time_t gps_to_unix( int gps_sec, int *leap );
void set_zone( const char *zone );
void reset_zone( void );
void set_zone_lho( void );
void set_zone_llo( void );

/* standard formats */
const char *iso_8601_format = "%Y-%m-%dT%H:%M:%S%z";
const char *rfc_822_format = "%a, %d %b %Y %H:%M:%S %z";
const char *unix_date_format = "%a %b %d %H:%M:%S %Z %Y";
const char *output_date_format = NULL;
const char *tz = NULL;

enum { GPS_EPOCH, UNIX_EPOCH } epoch = GPS_EPOCH; /* default is GPS epoch */
enum { SIDEREAL_HMS, SIDEREAL_RAD } sidereal_format = SIDEREAL_HMS; /* default isHH:MM:SS */
enum { OUTPUT_DEFAULT, OUTPUT_DATE, OUTPUT_SECS, OUTPUT_JD, OUTPUT_MJD, OUTPUT_LEAPS, OUTPUT_GMST } output_type = OUTPUT_DEFAULT;
enum { SITE_UNSPECIFIED, SITE_LHO, SITE_LLO } site = SITE_UNSPECIFIED;
int utc_flag = 1; /* default is universal time rather than local time */

int verbose = 0; /* default is brief output */

FILE *fp = NULL; /* set to stdin to read dates from stdin */

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
  if ( tz )
    set_zone( tz );

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

      if ( utc_flag )
        set_zone( "00" );  /* make zone UTC for date parsing unless local */

      if ( is_blank( arg ) )
        unix_sec = get_date( "now", 0 );
      else
        unix_sec = get_date( arg, 0 );

      if ( utc_flag )
        set_zone( tz ); /* if we set zone, set it back now */

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


static void output_secs( int gps_sec )
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

static void output_date( int gps_sec )
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

static void output_date_utc( int gps_sec )
{
  int utc_flag_save = utc_flag;
  utc_flag = 1;
  if ( verbose )
    printf( "Coordinated Universal Time: " );
  output_date( gps_sec );
  utc_flag = utc_flag_save;
  return;
}

static void output_date_local( int gps_sec )
{
  int utc_flag_save = utc_flag;
  utc_flag = 0;
  if ( verbose )
    printf( "Local Time: " );
  output_date( gps_sec );
  utc_flag = utc_flag_save;
  return;
}

static void output_lho_date( int gps_sec )
{
  set_zone_lho();
  if ( verbose )
    printf( "LHO " );
  output_date_local( gps_sec );
  set_zone( tz );
}

static void output_llo_date( int gps_sec )
{
  set_zone_llo();
  if ( verbose )
    printf( "LLO " );
  output_date_local( gps_sec );
  set_zone( tz );
}

static void output_date_site( int gps_sec )
{
  switch ( site )
  {
    case SITE_LHO:
      output_lho_date( gps_sec );
      break;
    case SITE_LLO:
      output_llo_date( gps_sec );
      break;
    case SITE_UNSPECIFIED:
      output_date( gps_sec );
      break;
    default:
      fprintf( stderr, "internal error: unrecognized site\n" );
      exit( 1 );
  }
  return;
}

static void output_jd( int gps_sec )
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

static void output_mjd( int gps_sec )
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

static void output_leaps( int gps_sec )
{
  int leaps;
  leaps = XLALLeapSeconds( gps_sec );
  if ( verbose )
    printf( "Leap Seconds (TAI-UTC) = " );
  printf( "%d\n", leaps );
  return;
}

static void output_sidereal_time( double sidereal_time_rad )
{
  double sidereal_time_hrs;
  double sidereal_time_min;
  double sidereal_time_sec;
  while ( sidereal_time_rad >= 2.0 * LAL_PI )
    sidereal_time_rad -= 2.0 * LAL_PI;
  while ( sidereal_time_rad < 0 )
    sidereal_time_rad += 2.0 * LAL_PI;
  switch ( sidereal_format )
  {
    case SIDEREAL_HMS:
      sidereal_time_hrs = 12.0 * sidereal_time_rad / LAL_PI;
      sidereal_time_min = 60.0 * modf( sidereal_time_hrs, &sidereal_time_hrs );
      sidereal_time_sec = 60.0 * modf( sidereal_time_min, &sidereal_time_min );
      sidereal_time_sec = floor( sidereal_time_sec + 0.5 ); /* round to nearest second */
      if ( verbose )
        printf( "Sidereal Time (HH:MM:SS) = " );
      printf( "%02.0f:%02.0f:%02d\n", sidereal_time_hrs, sidereal_time_min, (int)sidereal_time_sec );
      break;
    case SIDEREAL_RAD:
      if ( verbose )
        printf( "Sidereal Time (radians) = " );
      printf( "%.6f\n", sidereal_time_rad );
      break;
    default:
      fprintf( stderr, "internal error: unrecognized sidereal time format\n" );
      exit( 1 );
  }
  return;
}

static void output_gmst( int gps_sec )
{
  LIGOTimeGPS ligo_time;
  double gmst_rad;
  ligo_time.gpsSeconds = gps_sec;
  ligo_time.gpsNanoSeconds = 0;
  gmst_rad = XLALGreenwichMeanSiderealTime( &ligo_time );
  if ( verbose )
    printf( "Greenwich Mean " );
  output_sidereal_time( gmst_rad );
  return;
}

static void output_local_sidereal_time_lho( int gps_sec )
{
  LIGOTimeGPS ligo_time;
  double sidereal_time_rad;
  ligo_time.gpsSeconds = gps_sec;
  ligo_time.gpsNanoSeconds = 0;
  sidereal_time_rad  = XLALGreenwichMeanSiderealTime( &ligo_time );
  sidereal_time_rad += LHO_LONGITUDE_RAD_E;
  if ( verbose )
    printf( "LHO Local " );
  output_sidereal_time( sidereal_time_rad );
}

static void output_local_sidereal_time_llo( int gps_sec )
{
  LIGOTimeGPS ligo_time;
  double sidereal_time_rad;
  ligo_time.gpsSeconds = gps_sec;
  ligo_time.gpsNanoSeconds = 0;
  sidereal_time_rad  = XLALGreenwichMeanSiderealTime( &ligo_time );
  sidereal_time_rad += LLO_LONGITUDE_RAD_E;
  if ( verbose )
    printf( "LLO Local " );
  output_sidereal_time( sidereal_time_rad );
}

static void output_sidereal_time_site( int gps_sec )
{
  switch ( site )
  {
    case SITE_LHO:
      output_local_sidereal_time_lho( gps_sec );
      break;
    case SITE_LLO:
      output_local_sidereal_time_llo( gps_sec );
      break;
    case SITE_UNSPECIFIED:
      output_gmst( gps_sec );
      break;
    default:
      fprintf( stderr, "internal error: unrecognized site\n" );
      exit( 1 );
  }
  return;
}

void output( int gps_sec, int type )
{
  switch ( type )
  {
    case OUTPUT_DATE:
      output_date_site( gps_sec );
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
      output_sidereal_time_site( gps_sec );
      break;
    case OUTPUT_DEFAULT:
      if ( verbose )
      {
        int save_epoch = epoch;
        output_date_utc( gps_sec );
        output_date_local( gps_sec );
        output_lho_date( gps_sec );
        output_llo_date( gps_sec );
        epoch = GPS_EPOCH;
        output_secs( gps_sec );
        epoch = UNIX_EPOCH;
        output_secs( gps_sec );
        epoch = save_epoch;
        output_leaps( gps_sec );
        output_jd( gps_sec );
        output_mjd( gps_sec );
        output_gmst( gps_sec );
        output_local_sidereal_time_lho( gps_sec );
        output_local_sidereal_time_llo( gps_sec );
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

/* Sets TZ environment to a particular zone string, e.g., "US/Pacific".
 * Lots of static variables here... these need to be static since they
 * become part of the environment.  The argument is copied.
 * Pass NULL to this function to reset TZ to its original value.
 * Note: this function is ugly and a little bit dangerous, perhaps.*/
#define MAX_ENV_SIZE 1024
void set_zone( const char *zone )
{
#ifdef HAVE_SETENV
  static int first = 1;
  static char *tz_env_orig = NULL;
  if ( first )
  {
    first = 0;
    tz_env_orig = getenv( "TZ" );
  }

  if ( zone )
    setenv( "TZ", zone, 1 );
  else /* reset to original */
  {
    if ( tz_env_orig )
      setenv( "TZ", tz_env_orig, 1 );
    else
      unsetenv( "TZ" );
  }

  return;
#else
  static int first = 1;
  static char tz_equals_string[] = "TZ=";
  static char tz_string[] = "TZ";
  static char tz_env[MAX_ENV_SIZE];
  static char *tz_env_orig = NULL;

  if ( first )
  {
    first = 0;
    tz_env_orig = getenv( "TZ" );
  }

  if ( zone )
  {
    strncpy( tz_env, "TZ=", sizeof( tz_env ) - 1 );
    strncat( tz_env, zone, sizeof( tz_env ) - strlen( tz_env ) - 1 );
    putenv( tz_env );
  }
  else /* reset to original */
  {
    if ( tz_env_orig ) /* TZ was originally set */
    {
      strncpy( tz_env, "TZ=", sizeof( tz_env ) - 1 );
      strncat( tz_env, tz_env_orig, sizeof( tz_env ) - strlen( tz_env ) - 1 );
      putenv( tz_env );
    }
    else /* try to reset */
    {
      putenv( tz_equals_string ); /* set TZ to blank in case can't unset */
      putenv( tz_string ); /* this will unset TZ if possible */
      strcpy( tz_equals_string, "" ); /* a real hack at unsetting otherwise */
    }
  }

  return;
#endif
}

void reset_zone( void )
{
  set_zone( NULL );
  return;
}

void set_zone_lho( void )
{
  set_zone( "PST8PDT" );
  return;
}

void set_zone_llo( void )
{
  set_zone( "CST6CDT" );
  return;
}


void usage( const char *program, int exitcode )
{
  const char *name = "Name:\n\t%s\n\n";
  const char *synopsis = "Synopsis:\n\t%s [OPTIONS] DATE_STRING\n\t%s [OPTIONS] TIMER_SECONDS\n\n";
  const char *description = "Description:\n\
        Convert between various representations of time.\n\n";
  const char *options = "Options:\n\n\
        -d, --date\n\
                output a date string [default when input is timer seconds]\n\
        \n\
        -fFORMAT, --format=FORMAT\n\
                output format for date string [see date(1)]\n\
        \n\
        -G, --gps-epoch\n\
                use GPS epoch for timer (1980-01-06 00:00:00 UTC) [default]\n\
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
        -L, --local\n\
                date strings formatted for local time (current TZ setting)\
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
        -t, --timer\n\
                output timer seconds [default when input is a date string]\n\
        \n\
        -u, --utc\n\
                date strings formatted for Coordinated Universal Time [default]\n\
        \n\
        -U, --unix-epoch\n\
                use Unix epoch for timer (1970-01-01 00:00:00 UTC)\n\
        \n\
        -zTZ, --zone=TZ\n\
                set time zone to TZ (use with --local)\n\
        \n\
        -ZSITE, --site=SITE\n\
                set site to SITE which is [valid sites: LHO | LLO]\n\
        \n\
        --help  display this help and exit\n\
        \n\
        --version\n\
                output version information and exit\n\n";
  const char *environment = "Environment Variables:\n\
        TZ      The timezone to use when parsing or displaying dates.\n\n";
  const char *examples = "Examples:\n\
        Get current time in GPS seconds:\n\
        \n\
                %s now\n\
        \n\
        Get local time corresponding to GPS time 800000000 in RFC 2822 format:\n\
        \n\
                %s --rfc-2822 --local 800000000\n\
        \n\
        Find the local sidereal time at LHO three hours ago:\n\
        \n\
                %s --sidereal-time --site=LHO 3 hours ago\n\
        \n\
        Get the Julian day of Y2K (UTC is default unless --local is used):\n\
        \n\
                %s --julian-day Jan 1, 2000\n\
        \n\
        Find out all about now:\n\
        \n\
                %s --verbose\n\
        \n";
  fprintf( stdout, name, program );
  fprintf( stdout, synopsis, program, program );
  fprintf( stdout, "%s", description );
  fprintf( stdout, "%s", options );
  fprintf( stdout, "%s", environment );
  fprintf( stdout, examples, program, program, program, program, program );
  exit( exitcode );
}


char * parse_options( char *buf, int buflen, int argc, char **argv )
{
  struct option long_options[] =
  {
    { "date",           no_argument,            0,      'd'     },
    { "format",         required_argument,      0,      'f'     },
    { "gps-epoch",      no_argument,            0,      'G'     },
    { "help",           no_argument,            0,      'h'     },
    { "iso-8601",       no_argument,            0,      'I'     },
    { "jd",             no_argument,            0,      'j'     },
    { "julian-day",     no_argument,            0,      'j'     },
    { "leap-seconds",   no_argument,            0,      'l'     },
    { "local",          no_argument,            0,      'L'     },
    { "mjd",            no_argument,            0,      'm'     },
    { "modified-julian-day",  no_argument,      0,      'm'     },
    { "rfc-2822",       no_argument,            0,      'R'     },
    { "sidereal-time",  optional_argument,      0,      's'     },
    { "stdin",          no_argument,            0,      'S'     },
    { "utc",            no_argument,            0,      'u'     },
    { "universal",      no_argument,            0,      'u'     },
    { "unix-epoch",     no_argument,            0,      'U'     },
    { "timer",          no_argument,            0,      't'     },
    { "verbose",        no_argument,            0,      'v'     },
    { "version",        no_argument,            0,      'V'     },
    { "zone",           required_argument,      0,      'z'     },
    { "site",           required_argument,      0,      'Z'     },
    { 0, 0, 0, 0 }
  };
  char args[] = "+df:GhIjlLmRsSuUvVW:z:Z:";
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

      case 'd': /* date */
        output_type = OUTPUT_DATE;
        break;

      case 'f': /* format */
        output_date_format = optarg;
        break;

      case 'G': /* gps-epoch */
        epoch = GPS_EPOCH;
        break;

      case 'h': /* help */
        usage( program, 0 );
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

      case 'L': /* local time */
        utc_flag = 0;
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

      case 't': /* timer */
        output_type = OUTPUT_SECS;
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

      case 'z': /* zone */
        tz = optarg;
        break;

      case 'Z': /* site */
        if ( strcmp( optarg, "LHO" ) == 0 || strcmp( optarg, "lho" ) == 0 )
          site = SITE_LHO;
        else if ( strcmp( optarg, "LLO" ) == 0 || strcmp( optarg, "llo" ) == 0 )
          site = SITE_LLO;
        else
        {
          fprintf( stderr, "unrecognized detector site %s\n", optarg );
          fprintf( stderr, "expect either \"LHO\" or \"LLO\"\n" );
          exit( 1 );
        }
        break;

      case '?':
        /* fall-through */
      default:
        fprintf( stderr, "unknown error while parsing options\n" );
        usage( program, 1 );
    }
  }

  if ( optind == argc && fp )
  {
      char UNUSED *c;
      c = fgets( buf, buflen, fp );
  }
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
