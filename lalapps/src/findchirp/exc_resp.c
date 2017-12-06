/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Patrick Brady
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

#include <config.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lalapps.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

/* Calibration Constants for L1:LSC-ETMX        */
/* Based on mail from Gaby Mon, 24 Jun 2002     */
/* #define L1_DC_RESP 2.63  */      /* nm/count */
/* Based on e-log entry http://www.ligo-la.caltech.edu/ilog/pub/ilog.cgi?group=detector&date_to_view=02/01/2003&anchor_to_scroll_to=2003:02:02:00:28:03-brian */
#define L1_DC_RESP 0.39         /* +/- 0.02 (nm/ct) */
#define L1_PEND_F0 0.76         /* Hz */
#define L1_LENGTH  4000         /* m */

/* Calibration constants for H2:LSC-ETMX        */
/* Based on E7 calibration web page             */
/* #define H2_DC_RESP 1.33 */        /* nm/count */
/* Based on e-mail from Mike Landry Wed, 12 Feb 2003 11:21:12 -0800 */
#define H2_DC_RESP 1.05         /* nm/count */
#define H2_PEND_F0 0.74         /* Hz */
#define H2_LENGTH  2000         /* m */

/* Calibration constants for H1:LSC-ETMX        */
/* Based on e-mail from Mike Landry Wed, 12 Feb 2003 11:21:12 -0800 */
#define H1_DC_RESP 0.93         /* nm/count */
#define H1_PEND_F0 0.74         /* Hz */
#define H1_LENGTH  4000         /* m */

#define usgfmt \
"Usage: %s [options]\n"\
"Options [default in brackets]:\n"\
"  -h            print this message\n"\
"  -V            print version info\n"\
"  -v            verbose\n"\
"  -L flow       low frequency cutoff for response function [25.0]\n"\
"  -H fhigh      high frequency cutoff for response function [3000.0]\n"\
"  -n numpoints  number of points in response function [8192]\n"

#define usage( program ) fprintf( stderr, usgfmt, program )

extern char *optarg;
extern int optind, opterr, optopt;
extern int vrbflg;

int main ( int argc, char *argv[] )
{
  const char *program = argv[0];
  int   k;
  int   opt;
  int   numpts = 8192;
  REAL4 f_min = 25.0;
  REAL4 f_max = 3000.0;
  REAL4 df;
  FILE  *l1_fp = NULL, *h2_fp = NULL, *h1_fp = NULL;
  
  /* parse options */
  while ( 0 < ( opt = getopt( argc, argv, "hVvL:H:n:" ) ) )
  {
    switch ( opt )
    {
      case 'h':
        usage( program );
        return 0;
      case 'V':
        PRINT_VERSION( "hello" );
        return 0;
      case 'v':
        vrbflg = 1;
        break;
      case 'L':
        f_min = (REAL4) atof(optarg);
        break;
      case 'H':
        f_max = (REAL4) atof(optarg);
        break;
      case 'n':
        numpts = atoi(optarg);
        break;
      default:
        usage( program );
        return 1;
    }
  }
  if ( optind < argc )
  {
    usage( program );
    return 1;
  }
  
  df = (f_max - f_min) / ( (REAL4) numpts );

  l1_fp = fopen( "L1:LSC-ETMX_response", "w" );
  h2_fp = fopen( "H2:LSC-ETMX_response", "w" );
  h1_fp = fopen( "H1:LSC-ETMX_response", "w" );
  if ( !l1_fp || !h2_fp || !h1_fp )
  {
    fprintf( stderr, "error opening response function files for writing\n" );
    return 1;
  }

  fprintf( l1_fp, "# epoch = %lli\n", 0LL ); 
  fprintf( l1_fp, "# f0 = %e\n", (REAL8) f_min );
  fprintf( l1_fp, "# deltaF = %e\n", (REAL8) df );
  fprintf( h2_fp, "# epoch = %lli\n", 0LL ); 
  fprintf( h2_fp, "# f0 = %e\n", (REAL8) f_min );
  fprintf( h2_fp, "# deltaF = %e\n", (REAL8) df );
  fprintf( h1_fp, "# epoch = %lli\n", 0LL ); 
  fprintf( h1_fp, "# f0 = %e\n", (REAL8) f_min );
  fprintf( h1_fp, "# deltaF = %e\n", (REAL8) df );

  for ( k = 0; k < numpts; ++k )
  {
    REAL4 f = f_min + (REAL4) k * df;
    REAL4 l1_r = ( L1_PEND_F0 * L1_PEND_F0 * L1_DC_RESP * 1e-9 ) / 
      ( L1_LENGTH * f * f );
    REAL4 h2_r = ( H2_PEND_F0 * H2_PEND_F0 * H2_DC_RESP * 1e-9 ) / 
      ( H2_LENGTH * f * f );
    REAL4 h1_r = ( H1_PEND_F0 * H1_PEND_F0 * H1_DC_RESP * 1e-9 ) / 
      ( H1_LENGTH * f * f );
    fprintf( l1_fp, "%e %e\n", l1_r, 0.0 );
    fprintf( h2_fp, "%e %e\n", h2_r, 0.0 );
    fprintf( h1_fp, "%e %e\n", h1_r, 0.0 );
  }

  fclose( l1_fp );
  fclose( h2_fp );
  fclose( h1_fp );

  return 0;
}
