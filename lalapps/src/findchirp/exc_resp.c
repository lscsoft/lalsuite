#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lalapps.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>

RCSID( "$Id$" );

/* Calibration Constants for L1:LSC-ETMX        */
/* Based on mail from Gaby Mon, 24 Jun 2002     */
#define L1_DC_RESP 2.63         /* nm/count */
#define L1_PEND_F0 0.76         /* Hz */
#define L1_LENGTH  4000         /* m */

/* Calibration constants for H2:LSC-ETMX        */
/* Based on E7 calibration web page             */
#define H2_DC_RESP 1.33         /* nm/count */
#define H2_PEND_F0 0.74         /* Hz */
#define H2_LENGTH  2000         /* m */

#define usgfmt \
"Usage: %s [options]\n"\
"Options [default in brackets]:\n"\
"  -h            print this message\n"\
"  -V            print version info\n"\
"  -v            verbose\n"\
"  -d dbglvl     set debug level to dbglvl [0]\n"\
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
  const char *dbglvl  = NULL;
  int   k;
  int   opt;
  int   numpts = 8192;
  REAL4 f_min = 25.0;
  REAL4 f_max = 3000.0;
  REAL4 df;
  FILE  *l1_fp = NULL, *h2_fp = NULL;
  
  /* parse options */
  while ( 0 < ( opt = getopt( argc, argv, "hVvd:L:H:n:" ) ) )
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
      case 'd':
        dbglvl = optarg;
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
  if ( !l1_fp || !h2_fp )
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

  for ( k = 0; k < numpts; ++k )
  {
    REAL4 f = f_min + (REAL4) k * df;
    REAL4 l1_r = ( L1_PEND_F0 * L1_PEND_F0 * L1_DC_RESP * 1e-9 ) / 
      ( L1_LENGTH * f * f );
    REAL4 h2_r = ( H2_PEND_F0 * H2_PEND_F0 * H2_DC_RESP * 1e-9 ) / 
      ( H2_LENGTH * f * f );
    fprintf( l1_fp, "%e %e\n", l1_r, 0.0 );
    fprintf( h2_fp, "%e %e\n", h2_r, 0.0 );
  }

  fclose( l1_fp );
  fclose( h2_fp );

  return 0;
}
