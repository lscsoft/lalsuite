/*----------------------------------------------------------------------- 
 * 
 * File Name: releff.c
 *
 * Author: Brady, P. R.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <lal/LALError.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALConstants.h>
#include <lal/Thresholds.h>
#include <lalapps.h>

RCSID("$Id$");

#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "releff"

/* Usage format string. */
#define USAGE \
"Usage: %s [options]\n\n"\
"  --help                    display this message\n"\
"  --version                 print version information and exit\n"\
"  --debug-level LEVEL       set the LAL debug level to LEVEL\n"\
"\n"

int verbose = 0;
int computeAmplitude = 0;

int main ( int argc, char *argv[] )
{
  REAL8                 chi2;
  REAL8                 dof;
  REAL8                 faProb=0.0;
  REAL8                 fdProb=0.0;

  /* getopt arguments */
  struct option long_options[] =
  {
    {"verbose",                 no_argument,       &verbose,          1 },
    {"compute-amplitude",       no_argument,       &computeAmplitude, 1 },
    {"chisq",                   required_argument, 0,                'c'},
    {"dof",                     required_argument, 0,                'd'},
    {"false-alarm-prob",        required_argument, 0,                'e'},
    {"false-dismissal-prob",    required_argument, 0,                'f'},
    {"help",                    no_argument,       0,                'h'}, 
    {"debug-level",             required_argument, 0,                'z'},
    {"version",                 no_argument,       0,                'V'},
    {0, 0, 0, 0}
  };
  int c;


  /*
   * 
   * initialize things
   *
   */

  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "1" );
  setvbuf( stdout, NULL, _IONBF, 0 );

  /* default values */
  chi2 = 10.0;
  dof = 2;

  /* parse the arguments */
  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;

    c = getopt_long_only( argc, argv, 
        "hz:V", long_options, 
	&option_index );

    /* detect the end of the options */
    if ( c == -1 )
    {
      break;
    }
    
    switch ( c )
    {
      case 0:
        /* if this option set a flag, do nothing else now */
        if ( long_options[option_index].flag != 0 )
        {
          break;
        }
        else
        {
          fprintf( stderr, "Error parsing option %s with argument %s\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        break;

      case 'c':
        chi2 = atof( optarg );
        break;

      case 'd':
        dof = atof( optarg );
        break;

      case 'e':
        faProb = atof( optarg );
        break;

      case 'f':
        fdProb = atof( optarg );
        break;

      case 'h':
        /* help message */
        fprintf( stderr, USAGE , argv[0]);
        exit( 1 );
        break;

      case 'z':
        set_debug_level( optarg );
        break;

      case 'V':
        /* print version information and exit */
        fprintf( stdout, "Output amplitude comparison of \n" 
            "Patrick Brady and Duncan Brown\n"
            "CVS Version: " CVS_ID_STRING "\n"
            "CVS Tag: " CVS_NAME_STRING "\n" );
        exit( 0 );
        break;

      case '?':
        fprintf( stderr, USAGE , argv[0]);
        exit( 1 );
        break;

      default:
        fprintf( stderr, "Error: Unknown error while parsing options\n" );
        fprintf( stderr, USAGE, argv[0] );
        exit( 1 );
    }
  }

  if ( computeAmplitude )
  {
    REAL8 myChisq;

    myChisq = XLALChi2Threshold(dof, faProb);
    fprintf(stdout,"Threshold %e\n", myChisq);
    fprintf(stdout,"Amplitude %e\n", XLALRhoThreshold(myChisq, dof, fdProb));
  }
  else
  {
    fprintf(stdout,"%e\n", XLALChisqCdf(chi2, dof));
  }

  /***************************************************************************
   * 
   **************************************************************************/
  return 0;
}
