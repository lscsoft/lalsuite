/*-----------------------------------------------------------------------
 *
 * File Name: skynet.c
 *
 * Authors: Essinger-Hileman, T. and Owen, B.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */


#include <getopt.h>
#include <math.h>
#include <strings.h>
#include <stdio.h>

#include <lalapps.h>

#include <lal/AVFactories.h>
#include <lal/LALBarycenter.h>
#include <lal/LALConfig.h>
#include <lal/LALDatatypes.h>
#include <lal/LALInitBarycenter.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALXMGRInterface.h>
#include <lal/PtoleMetric.h>
#include <lal/StackMetric.h>
#include <lal/TwoDMesh.h>

#include "DopplerScan.h"


RCSID( "$Id$" );
 
/* type of parameter space metric to use  */
enum {
  undefined,     /* duh */
  ptolemetric,   /* PtoleMetric() */
  baryptolemaic, /* CoherentMetric() plus TBaryPtolemaic()  */
  ephemeris     /* CoherentMetric() plus TEphemeris() */
} metric_type = undefined;


enum {
  undetermined,
  hanford,
  livingston,
  geo,
  tama,
  virgo
} detector = undetermined;


/* Limits of sky search  */
REAL4 ra_min  = 0.0;
REAL4 ra_max  = LAL_TWOPI*0.999; /* BEN: need to fix this */
REAL4 dec_min = -LAL_PI_2;
REAL4 dec_max = LAL_PI_2;
REAL4 MAX_NODES = 8e6; /* limit on number of nodes for TwoDMesh  */

char *optarg = NULL; /* option argument for getopt_long() */
REAL8Vector *tevlambda;


int main( int argc, char *argv[] )
{

  LALStatus stat = blank_status;  /* status structure */
  REAL8Vector outputMetric; /* output argument for PulsarMetric */

  /* Define input variables and set default values */
  int begin            = 731265908;  /* start time of integration */
  REAL4 duration       = 1.8e5;      /* duration of integration (seconds) */
  REAL4 min_spindown   = 1e10;       /* minimum spindown age (seconds) */
  int spindown_order   = 1;          /* minimum spindown order */
  REAL4 mismatch       = 0.05;       /* mismatch threshold of mesh */
  REAL4 max_frequency  = 1e3;        /* maximum frequency of search (Hz) */


  /* define PtoleMetricIn structure and define default values */
  static PtoleMetricIn search;
  search.epoch.gpsSeconds = begin;
  search.duration = duration;
  search.maxFreq = max_frequency;



  /* other useful variables */
  FILE *fp;                       /* where to write the output */
  REAL4 f1;
  int option_index = 0;  /* getopt_long option index */
  int opt;               /* Argument for switch statement with getopt_long */
  int detector_argument; /* setting of detector location */


  /* Set getopt_long option arguments */
  static struct option long_options[] = {
    {"metric-type",       1, 0, 1},
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

  /* Set up. */
  lal_errhandler = LAL_ERR_EXIT;

  /* Parse command-line options. */
  while( (opt = getopt_long( argc, argv, "a:bc:defghjk", long_options, &option_index )) != -1 )
  {
    
    switch ( opt ) {

    case 1: /* Set type of metric for LALMetricWrapper */
      if( !strcmp( optarg, "ptolemetric" ) )
        metric_type = ptolemetric;
      if( !strcmp( optarg, "baryptolemaic" ) )
	metric_type = baryptolemaic;
      if( !strcmp( optarg, "ephemeris" ) )
	metric_type = ephemeris;
      break;

    case 2: /* start-gps-seconds option */
      begin = atoi ( optarg );
      break;

    case 3: /* Set detector site for LALMetricWrapper */
      if( !strcmp( optarg, "hanford" ) )
	detector = hanford;
      if( !strcmp( optarg, "livingston" ) )
	detector = livingston;
      if( !strcmp( optarg, "geo" ) )
	detector = geo;
      if( !strcmp( optarg, "tama" ) )
	detector = tama;
      if( !strcmp( optarg, "virgo" ) )
	detector = virgo;
      break;

    case 4:
      set_debug_level( optarg );
      break;

    case 5:
      duration = atoi (optarg);
      break;

    case 6:
      min_spindown = atoi ( optarg );
      break;

    case 7:
      spindown_order = atoi ( optarg );
      break;

    case 8:
      mismatch = atof( optarg );
      break;

    case 9:
      max_frequency = atoi( optarg );
      break;

    }/*switch( opt )*/

  }/*while ( getopt... )*/
printf( "parsed options...\n" );


 /* Set metric type  */
 switch( metric_type ) {

 case ptolemetric: 
   search.metricType = LAL_PMETRIC_COH_PTOLE_ANALYTIC;
   break;

 case baryptolemaic:
   search.metricType = LAL_PMETRIC_COH_PTOLE_NUMERIC;
   break;

 case ephemeris:
   search.metricType = LAL_PMETRIC_COH_EPHEM;

 default:
   printf( "Invalid metric type\n" );
   exit(1);

 } 
 printf( "Set metric type\n", metric_type );


 /* Set detector location  */
 switch( detector ) {

 case hanford:
   detector_argument = LALDetectorIndexLHODIFF;
   break;

 case livingston:
   detector_argument = LALDetectorIndexLLODIFF;
   break;

 case geo:
   detector_argument = LALDetectorIndexGEO600DIFF;
   break;

 case tama:
   detector_argument = LALDetectorIndexTAMA300DIFF;
   break;

 case virgo:
   detector_argument = LALDetectorIndexVIRGODIFF;
   break;

 default:
   printf( "Invalid detector argument\n" );
   exit(1);

 }
 printf( "Set detector location\n", detector );

 search.site = &lalCachedDetectors[detector_argument];

 LALPulsarMetric( &stat, &outputMetric, &search );

  LAL_CALL( LALCheckMemoryLeaks(), &stat );
  return 0;
} /* main() */
