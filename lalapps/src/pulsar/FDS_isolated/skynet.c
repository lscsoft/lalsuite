/*-----------------------------------------------------------------------

 *

 * File Name: skynet.c

 *

 * Author: Essinger-Hileman, Tom

 *

 * Revision: $Id$

 *

 *-----------------------------------------------------------------------

 */



#include <getopt.h>
#include <math.h>
#include <strings.h>
#include <stdio.h>

#include <lal/LALConfig.h>
#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>


char *optarg = NULL; /* option argument for getopt_long() */


int main( int argc, char *argv[] )

{

  /* Define option input variables */
  char metric_code[]   = "ptolemetric";  /* choice of LAL metric code
                           'ptolemetric'   = Ptolemetric(default),
                           'baryptolemaic' = CoherentMetric + DTBarycenter,
                           'ephemeris'     = CoherentMetric + DTEphemeris  */

  char detector[]      = "livingston";   /* choice of detector
                           'hanford'    = LIGO Hanford
                           'livingston' = LIGO Livingston (default)
                           'virgo'      = VIRGO
                           'geo'        = GEO600
                           'TAMA'       = TAMA300          */

  /* Define input variables and set default values */
  int begin            = 731265908;      /* start time of integration */
  int debug_level      = 0;              /* LAL debug level */
  REAL4 duration       = 1e5;            /* duration of integration */
  int min_spindown     = 1;              /* sets minimum spindown age */
  int spindown_order   = 1;              /* sets minimum spindown order */
  REAL4 mismatch       = 0.02;           /* sets maximum mismatch of mesh */
  REAL4 max_frequency  = 1e3;            /* maximum frequency of integration*/



  int option_index = 0; /* getopt_long option index */
  int opt; /* Argument for switch statement with getopt_long */

  /* Set getopt_long option arguments */
  static struct option long_options[] = {
    {"lal-metric",        1, 0, 1},
    {"start-gps-seconds", 1, 0, 2},
    {"detector",          1, 0, 3},
    {"debug-level",       1, 0, 4},
    {"integration-time",  1, 0, 5},
    {"min-spindown-age",  1, 0 ,6},
    {"spindown-order",    1, 0, 7},
    {"mismatch",          1, 0, 8},
    {"max-frequency",     1, 0, 9},
    {0, 0, 0, 0}
  };

  int metric_count;                      /* Counter for lal-metric option */
  char ptolemetric[] = "ptolemetric";    /* input strings, lal-metric option */
  char baryptolemaic[] = "baryptolemaic";
  char ephemeris[] = "ephemeris";

  int detector_count;                /* Counter for detector option */
  char hanford[]     = "hanford";    /* input strings for detector option */
  char livingston[]  = "livingston";
  char virgo[]       = "virgo";
  char geo[]         = "geo";
  char tama[]        = "tama";

  while( (opt = getopt_long( argc, argv, "a:bcdefghjk", long_options, &option_index )) != -1 )
  {
    
    switch ( opt )

    {

    case 1: /* lal-metric option */
      
      for (metric_count = 1; metric_count <= strlen(optarg); metric_count++)
	{
	  if ((optarg[metric_count] != ptolemetric[metric_count]) && 
              (optarg[metric_count] != baryptolemaic[metric_count]) && 
              (optarg[metric_count] != ephemeris[metric_count]))
	    {
	      printf("Invalid option argument for the --lal-metric option\n");
	      break;
	    }
	  else if ( metric_count == strlen(optarg) )
	    {
	      strcpy( metric_code, optarg);
	      break;
	    }
	}

    case 2: /* start-gps-seconds option */
      {
	begin = atoi ( optarg );
	break;
      }

    case 3:
 
      for (detector_count = 1; detector_count <= strlen(optarg); detector_count++)
	{
	  if ((optarg[detector_count] != hanford[detector_count]) && 
              (optarg[detector_count] != livingston[detector_count]) && 
              (optarg[detector_count] != virgo[detector_count]) &&
              (optarg[detector_count] != geo[detector_count]) &&
              (optarg[detector_count] != tama[detector_count]))
	    {
	      printf("Invalid option argument for the --detector option\n");
	      break;
	    }
	  else if ( detector_count == strlen(optarg) )
	    {
              strcpy( detector, optarg);
	      break;
	    }
	}

    case 4:
      {
	debug_level = atoi( optarg );
	break;
      }

    case 5:
      {
	duration = atoi (optarg);
	break;
      }

    case 6:
      {
	min_spindown = atoi ( optarg );
	break;
      }

    case 7:
      {
	spindown_order = atoi ( optarg );
	break;
      }

    case 8:
      {
	mismatch = atof( optarg );
	break;
      }

    case 9:
      {
	max_frequency = atoi( optarg );
	break;
      }

    }/*switch( opt )*/

  }/*while ( getopt... )*/


  /* Return input parameter values */
  printf("\nmetric_code is %s\n", metric_code);
  printf("integration start is %d seconds\n", begin); 
  printf("detector is %s\n", detector); 
  printf("LAL debug level is set to %d\n", debug_level); 
  printf("Integration duration is %f seconds\n", duration); 
  printf("Minimum spindown age is %d seconds\n", min_spindown); 
  printf("spindown order is %d\n", spindown_order); 
  printf("Mismatch of mesh is %f\n", mismatch); 
  printf("Maximum frequency of integration is %f Hz\n\n", max_frequency);



  return 0;

} /* main() */
