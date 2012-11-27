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
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/LALStdlib.h>
#include <lal/BlackHoleMode.h>
#define LAL_USE_COMPLEX_SHORT_MACROS
#include <lal/LALComplex.h>

#define l_invalid (INT_MIN + 1) /* invalid */
#define m_invalid (INT_MAX) /* invalid */
#define s_invalid (INT_MAX) /* invalid */
int l = l_invalid;
int m = m_invalid;
int s = s_invalid;
int leaver = 0;
int plotfmt = 0;

int usage( const char *program );
int parseargs( int argc, char **argv );

int main( int argc, char *argv[] )
{
  REAL8 avec[] = {0.0,0.1,0.2,0.3,0.4,0.45,0.49,0.4999};
  BlackHoleModeTable *modetab;
  size_t numa = sizeof( avec ) / sizeof( *avec );
  size_t i;

  lalDebugLevel = 7;
  XLALSetErrorHandler( XLALBacktraceErrorHandler );

  parseargs( argc, argv );

  fprintf( stdout, "# quasinormal mode table for l=%d m=%d s=%d ", l, m, s );
  if ( leaver )
    fprintf( stdout, "(Leaver's conventions)\n" );
  else
    fprintf( stdout, "(standard conventions)\n" );
  modetab = XLALCreateBlackHoleModeTable( l, m, s );
  if ( plotfmt )
    fprintf( stdout, "# Re(omega)\t  Im(omega)\t     a\n" );
  else
    fprintf( stdout, "#  a   \t        %s        \t        omega\n",
        leaver ? "A" : "E" );

  /* if we are doing a table, do positive spins first
   * otherwise do negative spins first */

  for ( i = 0; i < numa; ++i )
  {
    COMPLEX16 A, omega;
    REAL8 a;
    a = avec[i];
    XLALBlackHoleModeEigenvaluesLeaverT( &A, &omega, a, modetab );
    if ( ! leaver )
    {
      a *= 2.0;
      omega = cmulr(omega,0.5);
      A = caddr(A,s*(s+1));
    }
    if ( plotfmt )
      fprintf( stdout, "%+e\t%+e\t%e\n", LAL_REAL(omega), LAL_IMAG(omega), a );
    else
      fprintf( stdout, "%.4f \t(%.5f,%+.5f)\t(%+.6f,%.6f)\n", a, LAL_REAL(A), LAL_IMAG(A), LAL_REAL(omega), LAL_IMAG(omega) );
  }

  if ( ! plotfmt )
    fprintf( stdout, "\n" );

  for ( i = 0; i < numa; ++i )
  {
    COMPLEX16 A, omega;
    REAL8 a;
    a = avec[i];
    XLALBlackHoleModeEigenvaluesLeaverT( &A, &omega, -a, modetab );
    if ( ! leaver )
    {
      a *= 2.0;
      omega = cmulr(omega,0.5);
      A = caddr(A,s*(s+1));
    }
    if ( plotfmt )
      fprintf( stdout, "%+e\t%+e\t%e\n", -LAL_REAL(omega), LAL_IMAG(omega), a );
    else
      fprintf( stdout, "%.4f \t(%.5f,%+.5f)\t(%+.6f,%.6f)\n", a, LAL_REAL(A), LAL_IMAG(A), LAL_REAL(omega), LAL_IMAG(omega) );
  }

  XLALDestroyBlackHoleModeTable( modetab );
  LALCheckMemoryLeaks();

  return 0;
}

int parseargs( int argc, char **argv )
{
  struct option long_options[] = 
  {
    { "help", no_argument, 0, 'h' },
    { "leaver", no_argument, 0, 'L' },
    { "plotable", no_argument, 0, 'p' },
    { "l", required_argument, 0, 'l' },
    { "m", required_argument, 0, 'm' },
    { "s", required_argument, 0, 's' },
    { 0, 0, 0, 0 }
  };
  char args[] = "hLpl:m:s:";
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
      case 'p': /* plotable */
        plotfmt = 1;
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

  if ( l == l_invalid || m == m_invalid || s == s_invalid )
  {
    fprintf( stderr, "must specify l, m, and s\n" );
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
  fprintf( stderr, "\t-p, --plotable \tuse Leaver's conventions\n" );
  fprintf( stderr, "\t-l l           \t(required) set value of l, l>=0\n" );
  fprintf( stderr, "\t-m m           \t(required) set value of m, abs(m)<=l\n" );
  fprintf( stderr, "\t-s s           \t(required) set value of s, s<=0\n" );
  return 0;
}
