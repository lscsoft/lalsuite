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
#include <lal/Units.h>
#include <lal/TimeSeries.h>
#include <lal/BlackHoleMode.h>
#define LAL_USE_COMPLEX_SHORT_MACROS
#include <lal/LALComplex.h>

#define a_invalid -100.0 /* invalid */
#define i_invalid -361.0 /* invalid */
#define l_invalid (INT_MIN + 1) /* invalid */
#define m_invalid (INT_MAX) /* invalid */
REAL8 a = a_invalid;
REAL8 M = 0.0;
REAL8 r = 0.0;
REAL8 e = 0.0;
REAL8 i = i_invalid;
REAL8 q = 0.0;
int l = l_invalid;
int m = m_invalid;

int usage( const char *program );
int parseargs( int argc, char **argv );

int main( int argc, char *argv[] )
{
  const int s = -2;
  BlackHoleMode mode;
  COMPLEX16 amplitudePlus;
  COMPLEX16 amplitudeCross;

  lalDebugLevel = 7;
  XLALSetErrorHandler( XLALBacktraceErrorHandler );

  parseargs( argc, argv );

  XLALSetBlackHoleModeParams( &mode, a, l, m, s );

  XLALBlackHoleRingdownAmplitude( &amplitudePlus, &amplitudeCross, M, e, r, i, q, &mode );
  fprintf( stdout, "A_+ = %e\tphi_+ = %f degrees\n", LAL_CABS( amplitudePlus ), LAL_CARG( amplitudePlus ) / LAL_PI_180 );
  fprintf( stdout, "A_x = %e\tphi_x = %f degrees\n", LAL_CABS( amplitudeCross ), LAL_CARG( amplitudeCross ) / LAL_PI_180 );

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
    { "inclination", required_argument, 0, 'i' },
    { "azimuth", required_argument, 0, 'q' },
    { "energy", required_argument, 0, 'e' },
    { "distance", required_argument, 0, 'r' },
    { "l", required_argument, 0, 'l' },
    { "m", required_argument, 0, 'm' },
    { 0, 0, 0, 0 }
  };
  char args[] = "hM:a:i:q:e:r:l:m:";
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
      case 'M': /* mass */
        M = atof( optarg );
        break;
      case 'a': /* spin */
        a = atof( optarg );
        break;
      case 'i': /* inclination */
        i = LAL_PI_180 * atof( optarg );
        break;
      case 'q': /* azimuth */
        q = LAL_PI_180 * atof( optarg );
        break;
      case 'e': /* inclination */
        e = atof( optarg );
        break;
      case 'r': /* distance */
        r = atof( optarg );
        break;
      case 'l':
        l = atoi( optarg );
        break;
      case 'm':
        m = atoi( optarg );
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

  if ( a == a_invalid || l == l_invalid || m == m_invalid
      || M <= 0.0 || e <= 0.0 || e >= 1.0 || r <= 0.0 || i == i_invalid )
  {
    fprintf( stderr, "must specify mass, spin, distance, frac. energy loss, l, m\n" );
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
  fprintf( stderr, "\t-M Msolar      \t(required) set black hole mass (solar masses)\n" );
  fprintf( stderr, "\t-a a           \t(required) set value of a, -1<a<1\n" );
  fprintf( stderr, "\t-r distanceMpc \t(required) set distance (Mpc)\n" );
  fprintf( stderr, "\t-e fracEnergy  \t(required) set energy radiated (fraction of M)\n" );
  fprintf( stderr, "\t-i inclination \t(required) set inclination angle (degrees)\n" );
  fprintf( stderr, "\t-q azimuth     \t(default=0) set azimuthal angle (degrees)\n" );
  fprintf( stderr, "\t-l l           \t(required) set value of l, l>=0\n" );
  fprintf( stderr, "\t-m m           \t(required) set value of m, abs(m)<=l\n" );
  return 0;
}
