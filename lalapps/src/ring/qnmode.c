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

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/BlackHoleMode.h>
#define LAL_USE_COMPLEX_SHORT_MACROS
#include <lal/LALComplex.h>

#define a_invalid -100.0 /* invalid */
#define l_invalid (INT_MIN + 1) /* invalid */
#define m_invalid (INT_MAX) /* invalid */
#define s_invalid (INT_MAX) /* invalid */
REAL8 M = 0.0;
REAL8 a = a_invalid;
int l = l_invalid;
int m = m_invalid;
int s = s_invalid;
int leaver = 0;

int usage( const char *program );
int parseargs( int argc, char **argv );

int main( int argc, char *argv[] )
{
  COMPLEX16 A, omega;
  lalDebugLevel = 7;
  XLALSetErrorHandler( XLALBacktraceErrorHandler );

  parseargs( argc, argv );

  if ( ! leaver ) /* spin was not specified in Leaver conventions */
    a *= 0.5;

  XLALBlackHoleModeEigenvaluesLeaver( &A, &omega, a, l, m, s );
  if ( ! leaver )
    omega = cmulr(omega,0.5);


  /*
   * if mass has been specified, return f and Q
   * otherwise return omega in either leaver or standard conventions
   */

  fprintf( stdout, "using %s conventions\n", leaver ? "Leaver's" : "standard" );
  fprintf( stdout, "mode l = %d, m = %d, s = %d\n", l, m, s );
  fprintf( stdout, "spin a = %g (dimensionless)\n", leaver ? a : 2.0 * a );
  fprintf( stdout, "M * omega = (%+.6f,%+.6f)\n", LAL_REAL(omega), LAL_IMAG(omega) );
  if ( M != 0.0 )
  {
    REAL8 f, Q;
    f = fabs( LAL_REAL( omega ) / ( LAL_TWOPI * M * LAL_MTSUN_SI ) );
    Q = fabs( LAL_REAL( omega ) ) / ( -2.0 * LAL_IMAG( omega ) );

    fprintf( stdout, "mass M = %g solar masses\n", M );
    fprintf( stdout, "frequency f = %g Hz\n", f );
    fprintf( stdout, "quality Q = %g\n", Q );
  }

  LALCheckMemoryLeaks();

  return 0;
}

int parseargs( int argc, char **argv )
{
  struct option long_options[] = 
  {
    { "help", no_argument, 0, 'h' },
    { "leaver", no_argument, 0, 'L' },
    { "mass", required_argument, 0, 'M' },
    { "spin", required_argument, 0, 'a' },
    { "l", required_argument, 0, 'l' },
    { "m", required_argument, 0, 'm' },
    { "s", required_argument, 0, 's' },
    { 0, 0, 0, 0 }
  };
  char args[] = "hLM:a:l:m:s:";
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
        {
          fprintf( stderr, "error parsing option %s with argument %s\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
      case 'h': /* help */
        usage( argv[0] );
        exit( 0 );
      case 'L': /* leaver */
        leaver = 1;
        break;
      case 'M': /* mass */
        M = atof( optarg );
        break;
      case 'a': /* spin */
        a = atof( optarg );
        break;
      case 'l':
        l = atoi( optarg );
        break;
      case 'm':
        m = atoi( optarg );
        break;
      case 's':
        s = atoi( optarg );
        break;
      case '?':
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

  if ( a == a_invalid || l == l_invalid || m == m_invalid || s == s_invalid )
  {
    fprintf( stderr, "must specify a, l, m, and s\n" );
    usage( argv[0] );
    exit( 1 );
  }

  if ( leaver && M != 0.0 )
  {
    fprintf( stderr, "do not use both --leaver and --mass options\n" );
    usage( argv[0] );
    exit( 1 );
  }

  return 0;
}

int usage( const char *program )
{
  fprintf( stderr, "usage: %s [options]\n", program );
  fprintf( stderr, "options:\n" );
  fprintf( stderr, "\t-h, --help     \tprint this message and exit\n" );
  fprintf( stderr, "\t-L, --leaver   \tuse Leaver's conventions\n" );
  fprintf( stderr, "\t-M Msolar      \t(optional) set black hole mass (solar masses)\n" );
  fprintf( stderr, "\t-a a           \t(required) set value of a, -1<a<1\n" );
  fprintf( stderr, "\t-l l           \t(required) set value of l, l>=0\n" );
  fprintf( stderr, "\t-m m           \t(required) set value of m, abs(m)<=l\n" );
  fprintf( stderr, "\t-s s           \t(required) set value of s, s<=0\n" );
  return 0;
}
