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
REAL8 a = a_invalid;
int l = l_invalid;
int m = m_invalid;
int s = s_invalid;
int spherical = 0;
int plotfmt = 0;

int usage( const char *program );
int parseargs( int argc, char **argv );

int main( int argc, char *argv[] )
{
  BlackHoleMode mode;
  COMPLEX16 norm;
  REAL8 muvec[] = { -0.99, -0.95, -0.75, -0.55, -0.35, -0.15, 0.0, 0.15, 0.35, 0.55, 0.75, 0.95, 0.99};
  const char *fmt;
  int nummu;
  int i;
  lalDebugLevel = 7;
  XLALSetErrorHandler( XLALBacktraceErrorHandler );

  parseargs( argc, argv );

  if ( spherical && plotfmt )
    fmt = "%+e\t%+e\n";
  else if ( plotfmt )
    fmt = "%+e\t%+e\t%+e\n";
  else if ( spherical )
    fmt = "%+0.2f\t%+0.6f\n";
  else
    fmt = "%+0.2f\t(%+0.6f,%+0.6f)\n";

  if ( plotfmt )
  {
    nummu = 101;
    i = -nummu + 1;
  }
  else
  {
    nummu = sizeof(muvec)/sizeof(*muvec);
    i = 0;
  }

  XLALSetBlackHoleModeParams( &mode, a, l, m, s );
  XLALSpheroidalWaveFunctionNorm( &norm, &mode );

  for ( ; i < nummu; ++i )
  {
    REAL8 mu;
    COMPLEX16 swf;
    if ( plotfmt )
      mu = (double)i / (double)(nummu - 1);
    else
      mu = muvec[i];
    XLALSpheroidalWaveFunction1( &swf, mu, &mode );
    swf = cmul( swf, norm );
    if ( spherical )
      fprintf( stdout, fmt, mu, LAL_REAL(swf)/sqrt(2.0*LAL_PI) );
    else
      fprintf( stdout, fmt, mu, LAL_REAL(swf), LAL_IMAG(swf) );
  }

  return 0;
}

int parseargs( int argc, char **argv )
{
  struct option long_options[] = 
  {
    { "help", no_argument, 0, 'h' },
    { "spherical", no_argument, 0, 'Y' },
    { "plotable", no_argument, 0, 'p' },
    { "spin", required_argument, 0, 'a' },
    { "l", required_argument, 0, 'l' },
    { "m", required_argument, 0, 'm' },
    { "s", required_argument, 0, 's' },
    { 0, 0, 0, 0 }
  };
  char args[] = "hYpa:l:m:s:";
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
      case 'Y': /* spherical */
        spherical = 1;
        break;
      case 'p': /* plotfmt */
        plotfmt = 1;
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

  if ( spherical )
  {
    if ( a != a_invalid )
      fprintf( stderr, "a must be zero for spherical harmonics\n" );
    a = 0.0;
  }

  if ( a == a_invalid || l == l_invalid || m == m_invalid || s == s_invalid )
  {
    fprintf( stderr, "must specify a, l, m, and s\n" );
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
  fprintf( stderr, "\t-Y, --spherical\toutput spin weighted spherical harmonics\n" );
  fprintf( stderr, "\t-p, --plotable \toutput in plotable format\n" );
  fprintf( stderr, "\t-a a           \t(required unless spherical) set value of a, -1<a<1\n" );
  fprintf( stderr, "\t-l l           \t(required) set value of l, l>=0\n" );
  fprintf( stderr, "\t-m m           \t(required) set value of m, abs(m)<=l\n" );
  fprintf( stderr, "\t-s s           \t(required) set value of s, s<=0\n" );
  return 0;
}
