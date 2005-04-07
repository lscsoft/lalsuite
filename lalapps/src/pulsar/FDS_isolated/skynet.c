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
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALXMGRInterface.h>
#include <lal/PtoleMetric.h>
#include <lal/StackMetric.h>
#include <lal/TwoDMesh.h>


RCSID( "$Id$" );

/* type of parameter space metric to use */
enum {
  undefined,      /* duh */
  ptolemetric,    /* PtoleMetric() */
  baryptolemaic,  /* CoherentMetric() + TBaryPtolemaic() */
  ephemeris       /* CoherentMetric() + TEphemeris() */
} metric_type = undefined; 


/* Limits of sky search  */
REAL4 ra_min  = 0.0;
REAL4 ra_max  = LAL_TWOPI*0.999; /* BEN: need to fix this */
REAL4 dec_min = -LAL_PI_2;
REAL4 dec_max = LAL_PI_2;
REAL4 MAX_NODES = 2e6; /* limit on number of nodes for TwoDMesh  */


void getRange( LALStatus *, REAL4 [2], REAL4, void * );
void getMetric( LALStatus *, REAL4 [3], REAL4 [2], void * );

char *optarg = NULL; /* option argument for getopt_long() */
REAL8Vector *tevlambda;


int main( int argc, char *argv[] )
{

  LALStatus stat = blank_status;  /* status structure */

  /* Define option input variables */
  char detector[]      = "livingston";   /* choice of detector
                           'hanford'    = LIGO Hanford
                           'livingston' = LIGO Livingston (default)
                           'virgo'      = VIRGO
                           'geo'        = GEO600
                           'TAMA'       = TAMA300          */

  /* Define input variables and set default values */
  int begin            = 731265908;  /* start time of integration */
  REAL4 duration       = 1.8e5;      /* duration of integration (seconds) */
  REAL4 min_spindown   = 1e10;       /* minimum spindown age (seconds) */
  int spindown_order   = 1;          /* minimum spindown order */
  REAL4 mismatch       = 0.05;       /* mismatch threshold of mesh */
  REAL4 max_frequency  = 1e3;        /* maximum frequency of search (Hz) */

  /* structures for LAL functions */
  TwoDMeshNode *firstNode;
  static TwoDMeshParamStruc mesh;
  static PtoleMetricIn search;
  static MetricParamStruc tevparam;
  TwoDMeshNode *node;
  EphemerisData *eph;
  PulsarTimesParamStruc tevpulse;

  /* other useful variables */
  FILE *fp;                       /* where to write the output */
  REAL4 f1;
  int option_index = 0; /* getopt_long option index */
  int opt; /* Argument for switch statement with getopt_long */

  /* Set getopt_long option arguments */
  static struct option long_options[] = {
    {"metric-type",        1, 0, 1},
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

  int detector_count;                /* Counter for detector option */
  char hanford[]     = "hanford";    /* input strings for detector option */
  char livingston[]  = "livingston";
  char virgo[]       = "virgo";
  char geo[]         = "geo";
  char tama[]        = "tama";


  /* Set up. */
  lal_errhandler = LAL_ERR_EXIT;

  /* Parse command-line options. */
  while( (opt = getopt_long( argc, argv, "a:bcdefghjk", long_options, &option_index )) != -1 )
  {
    
    switch ( opt ) {

    case 1: /* metric-type option */
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

    case 3:
      for (detector_count = 1; detector_count <= strlen(optarg); detector_count++)
      {
        if ((optarg[detector_count] != hanford[detector_count]) && 
            (optarg[detector_count] != livingston[detector_count]) && 
            (optarg[detector_count] != virgo[detector_count]) &&
            (optarg[detector_count] != geo[detector_count]) &&
            (optarg[detector_count] != tama[detector_count]))
          printf("Invalid option argument for the --detector option\n");
        else if ( detector_count == strlen(optarg) )
          strcpy( detector, optarg);
      }
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


  /* Set TwoDMesh input parameters. */
  mesh.mThresh = mismatch;
  mesh.nIn = MAX_NODES;
  mesh.getRange = getRange;
  mesh.getMetric = getMetric;
  mesh.domain[0] = dec_min;
  mesh.domain[1] = dec_max;


  /* Set metric input parameters. */
  switch( metric_type ) {

  case ptolemetric:
    /* fill PtoleMetric() input structure */
    search.site = &lalCachedDetectors[LALDetectorIndexLLODIFF];
    search.position.system = COORDINATESYSTEM_EQUATORIAL;
    search.spindown = NULL;
    search.epoch.gpsSeconds = begin;
    search.epoch.gpsNanoSeconds = 0;
    search.duration = duration;
    search.maxFreq = max_frequency;
    /* tell TwoDMesh() to use PtoleMetric() */
    mesh.metricParams = (void *) &search;
    mesh.rangeParams = (void *) &search;
    break;

  case baryptolemaic:
    tevlambda = NULL;
    LAL_CALL( LALDCreateVector( &stat, &tevlambda, 3+spindown_order ), &stat );
    tevlambda->data[0] = max_frequency;
    tevparam.constants = &tevpulse;
    tevparam.n = 1;
    tevparam.errors = 0;
    tevparam.start = 0; /* start time relative to epoch */
    LAL_CALL( LALGetEarthTimes( &stat, &tevpulse ), &stat );
    tevpulse.t0 = 0.0; /* relative reference time for spindown defs */
    tevpulse.site = &lalCachedDetectors[LALDetectorIndexLLODIFF];
    tevpulse.epoch.gpsSeconds = begin;
    tevpulse.epoch.gpsNanoSeconds = 0;
    tevparam.deltaT = duration;
    /* set timing function */
    tevparam.dtCanon = LALDTBaryPtolemaic;

  case ephemeris:
    eph = (EphemerisData *) LALMalloc( sizeof(EphemerisData) );
    eph->ephiles.earthEphemeris = "earth00-04.dat";
    eph->ephiles.sunEphemeris = "sun00-04.dat";
    eph->leap = 13; /* for years 2000-2004; shouldn't matter if wrong */
    LAL_CALL( LALInitBarycenter( &stat, eph ), &stat );
    tevpulse.ephemeris = eph;
    tevparam.dtCanon = LALDTEphemeris;
    break;

  default:
    printf( "bad metric type\n" );
    exit(1);

  }
printf( "set input parameters...\n" );

  /* Create 2D mesh. */
  firstNode = NULL;
  LAL_CALL( LALCreateTwoDMesh( &stat, &firstNode, &mesh ), &stat );
  printf( "created %d nodes\n", mesh.nOut );
  if( mesh.nOut == MAX_NODES )
    printf( "This overflowed your limit. Try a smaller search.\n" );


  /* Write what we've got to file mesh.dat, if required */
  fp = fopen( "mesh.dat", "w" );

  for( node = firstNode; node; node = node->next )
    fprintf( fp, "%e %e\n", 
	       (double)((node->y)*180/LAL_PI), (double)((node->x)*180/LAL_PI));
  fclose( fp );


  /* Clean up and leave. */
  LALDestroyTwoDMesh( &stat, &firstNode, &mesh.nOut );
  printf( "destroyed %d nodes\n", mesh.nOut );
 
  LALCheckMemoryLeaks();
  return 0;
} /* main() */


/* This is the parameter range function as required by TwoDMesh. */

void getRange( LALStatus *stat, REAL4 y[2], REAL4 x, void *unused )
{
 
  y[0] = ra_min;
  y[1] = ra_max;

} /* getRange() */


/* This is the wrapped metric function as required by TwoDMesh. */
void getMetric( LALStatus *stat,
                REAL4 g[3],
                REAL4 x[2],
                void *params )
{

  REAL8Vector      *metric = NULL;   /* for output of metric */
  PtoleMetricIn    *Ppatch = NULL;   /* pointer for PtoleMetric params */
  REAL8             determinant;     /* Determinant of projected metric */


  Ppatch = params;

  /* set up shop  */
  INITSTATUS( stat, "getMetric", rcsid );
  ATTATCHSTATUSPTR( stat );
  TRY( LALDCreateVector( stat->statusPtr, &metric, 6 ), stat );
  
 
  Ppatch->position.longitude = x[1];
  Ppatch->position.latitude =  x[0];

  LALPtoleMetric( stat->statusPtr, metric, Ppatch );

  BEGINFAIL( stat )
    TRY( LALDDestroyVector( stat->statusPtr, &metric ), stat );
  ENDFAIL( stat );
  LALProjectMetric( stat->statusPtr, metric, 0 );
  BEGINFAIL( stat )
    TRY( LALDDestroyVector( stat->statusPtr, &metric ), stat );
  ENDFAIL( stat );

  determinant = metric->data[5]*metric->data[2]-pow(metric->data[4],2.0);
 

  /* Translate output. */
  g[1] = metric->data[2];
  g[0] = metric->data[5];
  g[2] = metric->data[4];
 
  
  /* Clean up and leave. */
  TRY( LALDDestroyVector( stat->statusPtr, &metric ), stat );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
  

} /* getMetric() */
