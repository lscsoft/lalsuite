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



NRCSID( SKYNETC, "$$" );


/* Limits of sky search  */
REAL4 ra_min  = 0.0;
REAL4 ra_max  = LAL_PI_2;
REAL4 dec_min = 0.0;
REAL4 dec_max = LAL_PI_2;

REAL4 MAX_NODES = 1e6; /* limit on number of nodes for TwoDMesh  */

void getRange( LALStatus *, REAL4 [2], REAL4, void * );
void getMetric( LALStatus *, REAL4 [3], REAL4 [2], void * );



char *optarg = NULL; /* option argument for getopt_long() */


int main( int argc, char *argv[] )

{

  static LALStatus stat;  /* status structure */
  FILE *fp;               /* where to write the output */

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
  int begin            = 731265908;  /* start time of integration */
  int debug_level      = 0;          /* LAL debug level */
  REAL4 duration       = 1e5;        /* duration of integration (seconds) */
  REAL4 min_spindown   = 1e10;       /* minimum spindown age (seconds) */
  int spindown_order   = 1;          /* minimum spindown order */
  REAL4 mismatch       = 0.05;       /* mismatch threshold of mesh */
  REAL4 max_frequency  = 1e3;        /* maximum frequency of search (Hz) */

  /* Define structures for TwoDMesh and Ptolemetric  */
  TwoDMeshNode *firstNode;
  static TwoDMeshParamStruc mesh;
  static PtoleMetricIn search;


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
  printf("Minimum spindown age is %f seconds\n", min_spindown); 
  printf("spindown order is %d\n", spindown_order); 
  printf("Mismatch of mesh is %f\n", mismatch); 
  printf("Maximum frequency of integration is %f Hz\n\n", max_frequency);


 /* Set the mesh input parameters. */
  mesh.mThresh = mismatch;
  mesh.nIn = MAX_NODES;
  mesh.getRange = getRange;
  mesh.getMetric = getMetric;
  mesh.metricParams = (void *) &search;
  mesh.domain[0] = dec_min;
  mesh.domain[1] = dec_max;
  mesh.rangeParams = (void *) &search;
  

  /* Set detector location */
  search.site = &lalCachedDetectors[LALDetectorIndexLLODIFF];

  /* Ptolemetric constants */
  search.position.system = COORDINATESYSTEM_EQUATORIAL;
  search.spindown = NULL;
  search.epoch.gpsSeconds = begin;
  search.epoch.gpsNanoSeconds = 0;
  search.duration = duration;
  search.maxFreq = max_frequency;

  /* Create mesh */
  firstNode = NULL;
  LALCreateTwoDMesh( &stat, &firstNode, &mesh );
  if( stat.statusCode )
    return stat.statusCode;
  printf( "created %d nodes\n", mesh.nOut );
  if( mesh.nOut == MAX_NODES )
    printf( "This overflowed your limit. Try a smaller search.\n" );


  /* Write what we've got to file mesh.dat, if required */
  TwoDMeshNode *node;
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
  INITSTATUS( stat, "getMetric", SKYNETC );
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
