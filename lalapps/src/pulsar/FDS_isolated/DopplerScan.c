/************************************ <lalVerbatim file="DopplerScanCV">
Author: Prix, Reinhard
$Id$
************************************* </lalVerbatim> */

/* some parts of this file are based on code from PtoleMeshTest, written
   by Ben Owen and Ian Jones */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{DopplerScan}
\label{ss:DopplerScan.c}

Module to generate subsequent sky-positions.

\subsubsection*{Prototypes}
\input{DopplerScanCP}
\idx{NextSkyPosition()}

\subsubsection*{Description}

This module generates subsequent sky-positions corresponding to a 
given search area and resolution 

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{DopplerScanCV}}

******************************************************* </lalLaTeX> */
#include <math.h>

#include <lal/LALStdlib.h>
#include <lal/DetectorSite.h>
#include <lal/StackMetric.h>
#include <lal/AVFactories.h>
#include <lal/LALError.h>
#include <lal/LALXMGRInterface.h>
#include <lal/StringInput.h>

#include <lal/ConfigFile.h>
#include <lal/Velocity.h>


#include "DopplerScan.h"


NRCSID( DOPPLERSCANC, "$Id$" );

/* TwoDMesh() can have either of two preferred directions of meshing: */
enum {
  ORDER_ALPHA_DELTA,
  ORDER_DELTA_ALPHA
};

static int meshOrder = ORDER_DELTA_ALPHA;

#define MIN(x,y) (x < y ? x : y)
#define MAX(x,y) (x > y ? x : y)

#define TRUE (1==1)
#define FALSE (1==0)

extern INT4 lalDebugLevel;

/* some empty structs for initializations */
static TwoDMeshParamStruc empty_meshpar;
static PtoleMetricIn empty_metricpar;
static DopplerScanGrid empty_grid;
static MetricParamStruc empty_MetricParamStruc;
static PulsarTimesParamStruc empty_PulsarTimesParamStruc;

/* internal prototypes */
void getRange( LALStatus *stat, REAL4 y[2], REAL4 x, void *params );
void getMetric( LALStatus *status, REAL4 g[3], REAL4 skypos[2], void *params );
void LALMetricWrapper (LALStatus *stat, REAL8Vector *metric, const PtoleMetricIn *params);
void LALTemplateDistance (LALStatus *stat, REAL8 *dist, const DopplerPosition *pos1, const DopplerPosition *pos2);
	
void ConvertTwoDMesh2Grid ( LALStatus *stat, DopplerScanGrid **grid, const TwoDMeshNode *mesh2d, const SkyRegion *region );

BOOLEAN pointInPolygon ( const SkyPosition *point, const SkyRegion *polygon );

void gridFlipOrder ( TwoDMeshNode *grid );

void buildFlatGrid (LALStatus *stat, DopplerScanGrid **grid, const SkyRegion *region, REAL8 dAlpha, REAL8 dDelta);
void buildIsotropicGrid (LALStatus *stat, DopplerScanGrid **grid, const SkyRegion *skyRegion, REAL8 dAlpha, REAL8 dDelta);
void buildMetricGrid (LALStatus *stat, DopplerScanGrid **grid, SkyRegion *skyRegion,  const DopplerScanInit *init);
void loadSkyGridFile (LALStatus *stat, DopplerScanGrid **grid, const CHAR *fname);

void plotGrid (LALStatus *stat, DopplerScanGrid *grid, const SkyRegion *region, const DopplerScanInit *init);

void freeGrid (DopplerScanGrid *grid);
void printFrequencyShifts ( LALStatus *stat, const DopplerScanState *scan, const DopplerScanInit *init);

/*----------------------------------------------------------------------*/
/* <lalVerbatim file="DopplerScanCP"> */
void
InitDopplerScan( LALStatus *stat, 
		 DopplerScanState *scan, 	/* the un-initialized scan-structure */
		 const DopplerScanInit *init)		/* init-params */
{ /* </lalVerbatim> */
  DopplerScanGrid *node; 

  INITSTATUS( stat, "DopplerScanInit", DOPPLERSCANC );
  ATTATCHSTATUSPTR (stat); 

  /* This traps coding errors in the calling routine. */
  ASSERT ( scan != NULL, stat, DOPPLERSCANH_ENONULL, DOPPLERSCANH_MSGENONULL );  
  ASSERT ( scan->state == STATE_IDLE, stat, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );
  ASSERT ( init->gridType < GRID_LAST, stat, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );
  
  /* trap some abnormal input */
  if ( (init->gridType != GRID_FILE) && (init->skyRegion == NULL) ) 
    {
      LALPrintError ( "\nERROR: No sky-region was specified!\n\n");
      ABORT (stat,  DOPPLERSCANH_ENULL ,  DOPPLERSCANH_MSGENULL );
    } 
  if ( (init->gridType == GRID_FILE) && (init->skyGridFile == NULL) )
    {
      LALPrintError ("\nERROR: no skyGridFile has been specified!\n\n");
      ABORT (stat,  DOPPLERSCANH_ENULL ,  DOPPLERSCANH_MSGENULL );
    }

  /* general initializations */
  scan->grid = NULL;  
  scan->gridNode = NULL;
  
  if (init->gridType != GRID_FILE ) {
    TRY (ParseSkyRegion (stat->statusPtr, &(scan->skyRegion), init->skyRegion ), stat);

    if (scan->skyRegion.numVertices == 2)	{ /* this is an anomaly! Allowed are either 1 or >= 3 */
      ABORT (stat, DOPPLERSCANH_E2DSKY, DOPPLERSCANH_MSGE2DSKY);
    }
  } /* if gridType != GRID_FILE */

 
  switch (init->gridType)
    {
    case GRID_FLAT:		/* flat-grid: constant dAlpha, dDelta */

      TRY ( buildFlatGrid ( stat->statusPtr, &(scan->grid), &(scan->skyRegion), init->dAlpha, init->dDelta), stat);

      break;

    case GRID_ISOTROPIC: 	/* variant of manual stepping: try to produce an isotropic mesh */

      TRY ( buildIsotropicGrid ( stat->statusPtr, &(scan->grid), &(scan->skyRegion), init->dAlpha, init->dDelta), stat);

      break;

    case GRID_METRIC:

      TRY ( buildMetricGrid (stat->statusPtr, &(scan->grid), &(scan->skyRegion), init), stat);

      break;

    case GRID_FILE:

      TRY ( loadSkyGridFile (stat->statusPtr, &(scan->grid), init->skyGridFile), stat);

      break;

    default:
      LALPrintError ("\nUnknown grid-type `%d`\n\n", init->gridType);
      ABORT ( stat, DOPPLERSCANH_EMETRIC, DOPPLERSCANH_MSGEMETRIC);
      break;

    } /* switch (metric) */

  /* NOTE: we want to make sure we return at least one grid-point: 
   * so check if we got one, and if not, we return the
   * first point of the skyRegion-polygon as a grid-point
   */
  if (scan->grid == NULL)
    {
      scan->grid = LALCalloc (1, sizeof(DopplerScanGrid));
      if (scan->grid == NULL) {
	ABORT (stat, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM);
      }
      scan->grid->alpha = scan->skyRegion.vertices[0].longitude;
      scan->grid->delta = scan->skyRegion.vertices[0].latitude;
      
    } /* no points found inside of sky-region */

  /* initialize grid-pointer to first node in list */
  scan->gridNode = scan->grid; 	

  /* count number of nodes in our Dopplerscan-grid */
  scan->numGridPoints = 0;
  node = scan->grid;
  while (node)
    {
      scan->numGridPoints ++;
      node = node->next;
    }


  if (lalDebugLevel >= 1)
    LALPrintError ("\nFinal Scan-grid has %d nodes\n", scan->numGridPoints);
  if (lalDebugLevel >= 3)
    {
      LALPrintError ("\nDEBUG: plotting sky-grid into file 'mesh_debug.agr' ...");
      TRY( plotGrid (stat->statusPtr, scan->grid, &(scan->skyRegion), init), stat);
      LALPrintError (" done. \n");
      LALPrintError ("\nDEBUG: outputing predicted frequency-shifts 'dFreq.pred' \n");
      TRY ( printFrequencyShifts (stat->statusPtr, scan, init), stat);
      LALPrintError (" done. \n");
    }

  scan->state = STATE_READY;

  /* clean up */
  DETATCHSTATUSPTR (stat);

  RETURN( stat );

} /* InitDopplerScan() */


/*----------------------------------------------------------------------*/
/* <lalVerbatim file="DopplerScanCP"> */
void
FreeDopplerScan (LALStatus *stat, DopplerScanState *scan)
{ /* </lalVerbatim> */

  INITSTATUS( stat, "FreeDopplerScan", DOPPLERSCANC);
  ATTATCHSTATUSPTR (stat);

  /* This traps coding errors in the calling routine. */
  ASSERT( scan != NULL, stat, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL );
  
  if ( scan->state == STATE_IDLE )
    {
      LALError (stat, "\nTrying to free an uninitialized DopplerScan.\n");
      ABORT (stat, DOPPLERSCANH_ENOTREADY, DOPPLERSCANH_MSGENOTREADY);
    }
  else if ( scan->state != STATE_FINISHED )
    LALWarning (stat, "freeing unfinished DopplerScan.");

  freeGrid (scan->grid);

  scan->grid = scan->gridNode = NULL;

  if (scan->skyRegion.vertices)
    LALFree (scan->skyRegion.vertices);
  scan->skyRegion.vertices = NULL;
    
  scan->state = STATE_IDLE;

  DETATCHSTATUSPTR (stat);
  RETURN( stat );

} /* FreeDopplerScan() */

/*----------------------------------------------------------------------
 *  helper-function: free the linked list containing the grid 
 *----------------------------------------------------------------------*/
void
freeGrid (DopplerScanGrid *grid)
{
  DopplerScanGrid *node, *next;

  if ( grid == NULL)
    return;

  node = grid;

  while (node)
    {
      next = node->next;
      LALFree (node);
      node = next;
    } /* while node */

  return;
} /* freeGrid() */

/*----------------------------------------------------------------------
 * NextDopplerPos(): step through template-grid 
 *----------------------------------------------------------------------*/
/* <lalVerbatim file="DopplerScanCP"> */
void
NextDopplerPos( LALStatus *stat, DopplerPosition *pos, DopplerScanState *scan)
{ /* </lalVerbatim> */

  INITSTATUS( stat, "NextDopplerPos", DOPPLERSCANC);

  /* This traps coding errors in the calling routine. */
  ASSERT( pos != NULL, stat, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL );

  switch( scan->state )
    {
    case STATE_IDLE:
    case STATE_FINISHED:
      ABORT (stat, DOPPLERSCANH_ENOTREADY, DOPPLERSCANH_MSGENOTREADY);
      break;

    case STATE_READY:  
      if (scan->gridNode == NULL) 	/* we're done */
	scan->state = STATE_FINISHED;
      else
	{
	  pos->skypos.longitude = scan->gridNode->alpha;
	  pos->skypos.latitude =  scan->gridNode->delta;

	  pos->spindowns.length = 0; /*  FIXME: no spindowns for now */

	  scan->gridNode = scan->gridNode->next;
	}
      break;

    } /* switch scan->stat */

  RETURN( stat );

} /* NextSkyPos() */


/* **********************************************************************
   The following 2 helper-functions for TwoDMesh() have been adapted 
   from Ben's PtoleMeshTest.

   NOTE: Currently we are very expicitly restricted to 2D searches!! 
   FIXME: generalize to N-dimensional parameter-searches
********************************************************************** */

/* ----------------------------------------------------------------------
 * This is the parameter range function as required by TwoDMesh. 
 *
 * NOTE: for the moment we only provide a trival range as defined by the  
 * rectangular parameter-area [ a1, a2 ] x [ d1, d2 ]
 * 
 *----------------------------------------------------------------------*/
void getRange( LALStatus *stat, REAL4 y[2], REAL4 x, void *params )
{
  SkyRegion *region = (SkyRegion*)params;
  REAL4 nix;

  /* Set up shop. */
  INITSTATUS( stat, "getRange", DOPPLERSCANC );
  /*   ATTATCHSTATUSPTR( stat ); */

  nix = x;	/* avoid compiler warning */
  
  /* for now: we return the fixed y-range, indendent of x */
  if (meshOrder == ORDER_ALPHA_DELTA)
    {
      y[0] = region->lowerLeft.latitude;
      y[1] = region->upperRight.latitude;
    }
  else
    {
      y[0] = region->lowerLeft.longitude;
      y[1] = region->upperRight.longitude;
    }


  /* Clean up and leave. */
  /*   DETATCHSTATUSPTR( stat ); */

  RETURN( stat );
} /* getRange() */


/* ----------------------------------------------------------------------
 * This is a wrapper for the metric function as required by TwoDMesh. 
 *
 *
 * NOTE: this will be called by TwoDMesh(), therefore
 * skypos is in internalOrder, which is not necessarily ORDER_ALPHA_DELTA!!
 * 
 *----------------------------------------------------------------------*/
void getMetric( LALStatus *stat, REAL4 g[3], REAL4 skypos[2], void *params )
{
  INT2 dim;  	/* dimension of (full) parameter space (incl. freq) */
  REAL8Vector   *metric = NULL;  /* for output of metric */
  PtoleMetricIn *metricpar = (PtoleMetricIn*) params;

  /* Set up shop. */
  INITSTATUS( stat, "getMetric", DOPPLERSCANC );
  ATTATCHSTATUSPTR( stat );

  /* currently we use only f0, alpha, delta: -> 3D metric */
  dim = 3;

  TRY( LALDCreateVector( stat->statusPtr, &metric, dim*(dim+1)/2 ), stat );

  /* Call the metric function. (Ptole or Coherent, which is handled by wrapper) */
  if (meshOrder == ORDER_ALPHA_DELTA)
    {
      metricpar->position.longitude = skypos[0];
      metricpar->position.latitude =  skypos[1];
    }
  else
    {
      metricpar->position.longitude = skypos[1];
      metricpar->position.latitude =  skypos[0];
    }

  /* before we call the metric: make sure the sky-position  is "normalized" */
  TRY ( LALNormalizeSkyPosition (stat->statusPtr, &(metricpar->position), &(metricpar->position)), stat);

  TRY ( LALMetricWrapper( stat->statusPtr, metric, metricpar), stat);

  BEGINFAIL( stat )
    TRY( LALDDestroyVector( stat->statusPtr, &metric ), stat );
  ENDFAIL( stat );

  TRY (LALProjectMetric( stat->statusPtr, metric, 0 ), stat);

  BEGINFAIL( stat )
    TRY( LALDDestroyVector( stat->statusPtr, &metric ), stat );
  ENDFAIL( stat );


  /* the general indexing scheme is g_ab for a>=b: index = b + a*(a+1)/2 */
#define INDEX_AA (1 + 1*(1+1)/2)	/* g_aa */
#define INDEX_DD (2 + 2*(2+1)/2)	/* g_dd */
#define INDEX_AD (1 + 2*(2+1)/2)	/* g_ad */

  /* Translate output. Careful about the coordinate-order here!! */
  if (meshOrder == ORDER_ALPHA_DELTA)
    {
      g[0] = metric->data[INDEX_AA]; /* gxx */
      g[1] = metric->data[INDEX_DD]; /* gyy */
    }
  else
    {
      g[0] = metric->data[INDEX_DD]; /* gxx */
      g[1] = metric->data[INDEX_AA]; /* gyy */
    }

  g[2] = metric->data[INDEX_AD]; /* gxy = g21: 1 + 2*(2+1)/2 = 4; */

  /*
  if (lalDebugLevel >= 3)
    printf ("DEBUG: (ra,dec)=(%f,%f): gxx = %f, gyy = %f, gxy = %f\n", 
	    metricpar->position.longitude, metricpar->position.latitude,  g[0], g[1], g[2]);
  */
 
  /* Clean up and leave. */
  TRY( LALDDestroyVector( stat->statusPtr, &metric ), stat );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
} /* getMetric() */


/*----------------------------------------------------------------------
 * Debug helper for mesh and metric stuff
 *----------------------------------------------------------------------*/
#define SPOKES 60  /* spokes for ellipse-plots */
#define NUM_SPINDOWN 0       /* Number of spindown parameters */

void 
plotGrid (LALStatus *stat, 
	   DopplerScanGrid *grid, 
	   const SkyRegion *region, 
	   const DopplerScanInit *init)
{
  FILE *fp = NULL;
  DopplerScanGrid *node;
  REAL8 alpha, delta;
  UINT4 set, i;
  UINT4 dim;

  const CHAR *xmgrHeader = 
    "@version 50103\n"
    "@title \"Sky-grid\"\n"
    "@world xmin -0.1\n"
    "@world xmax 6.4\n"
    "@world ymin -1.58\n"
    "@world ymax 1.58\n"
    "@xaxis label \"alpha\"\n"
    "@yaxis label \"delta\"\n";

  /* Set up shop. */
  INITSTATUS( stat, "plotGrid", DOPPLERSCANC );
  ATTATCHSTATUSPTR( stat );

  /* currently we use only f0, alpha, delta: -> 3D metric */
  dim = 3;

  fp = fopen ("mesh_debug.agr", "w");

  if( !fp ) {
    ABORT ( stat, DOPPLERSCANH_ESYS, DOPPLERSCANH_MSGESYS );
  }
  
  fprintf (fp, xmgrHeader);

  set = 0;

  /* Plot boundary. (if given) */
  if (region->vertices)
    {
      fprintf( fp, "@target s%d\n@type xy\n", set );

      for( i = 0; i < region->numVertices; i++ )
	{
	  fprintf( fp, "%e %e\n", region->vertices[i].longitude, region->vertices[i].latitude );
	}
      fprintf (fp, "%e %e\n", region->vertices[0].longitude, region->vertices[0].latitude ); /* close contour */
      
      set ++;
    }
      
  /* Plot mesh points. */
  fprintf( fp, "@s%d symbol 9\n@s%d symbol size 0.33\n", set, set );
  fprintf( fp, "@s%d line type 0\n", set );
  fprintf( fp, "@target s%d\n@type xy\n", set );

  for( node = grid; node; node = node->next )
  {
    fprintf( fp, "%e %e\n", node->alpha, node->delta );
  }


  /* if metric is given: plot ellipses */
  if ( ( init->metricType > LAL_PMETRIC_NONE) && (init->metricType < LAL_PMETRIC_LAST) )
    {
      REAL8Vector  *metric = NULL;   /* Parameter-space metric: for plotting ellipses */
      REAL8 gaa, gad, gdd, angle, smaj, smin;
      REAL8 mismatch = init->metricMismatch;
      PtoleMetricIn metricPar;

      set++;

      /* set up the metric parameters common to all sky-points */
      metricPar.metricType = init->metricType;
      metricPar.position.system = COORDINATESYSTEM_EQUATORIAL;
      metricPar.spindown = LALCalloc ( 1, sizeof(REAL4Vector) );
      metricPar.spindown->length=0;
      metricPar.spindown->data=NULL;
      metricPar.epoch = init->obsBegin;
      metricPar.duration = init->obsDuration;
      metricPar.maxFreq = init->fmax;
      metricPar.site = init->Detector;

      TRY( LALDCreateVector( stat->statusPtr, &metric, dim*(dim+1)/2 ), stat);

      node = grid;
      while (node)
	{
	  alpha =  node->alpha;
	  delta =  node->delta;

	  /* Get the metric at this skypos */
	  /* only need the update the position, the other
	   * parameter have been set already ! */
	  metricPar.position.longitude = alpha;
	  metricPar.position.latitude  = delta;

	  /* make sure we "normalize" point before calling metric */
	  TRY( LALNormalizeSkyPosition (stat->statusPtr, &(metricPar.position), &(metricPar.position) ), stat);

	  TRY( LALMetricWrapper( stat->statusPtr, metric, &metricPar), stat);

	  TRY (LALProjectMetric( stat->statusPtr, metric, 0 ), stat);

	  gaa = metric->data[INDEX_AA];
	  gad = metric->data[INDEX_AD];
	  gdd = metric->data[INDEX_DD];
	  
	  /* Semiminor axis from larger eigenvalue of metric. */
	  smin = gaa+gdd + sqrt( pow(gaa-gdd,2) + pow(2*gad,2) );
	  smin = sqrt(2.0* mismatch/smin);
	  /* Semiminor axis from smaller eigenvalue of metric. */
	  smaj = gaa+gdd - sqrt( pow(gaa-gdd,2) + pow(2*gad,2) );
	  smaj = sqrt(2.0* mismatch /smaj);
	  
	  /* Angle of semimajor axis with "horizontal" (equator). */
	  angle = atan2( gad, mismatch /smaj/smaj-gdd );
	  if (angle <= -LAL_PI_2) angle += LAL_PI;
	  if (angle > LAL_PI_2) angle -= LAL_PI;

	  /*
	  printf ("alpha=%f delta=%f\ngaa=%f gdd=%f gad=%f\n", alpha, delta, gaa, gdd, gad);
	  printf ("smaj = %f, smin = %f angle=%f\n", smaj, smin, angle);
	  */
	  set ++;
	  /* Print set header. */
	  fprintf( fp, "@target G0.S%d\n@type xy\n", set);
	  fprintf( fp, "@s%d color (0,0,0)\n", set );
	  
	  /* Loop around patch ellipse. */
	  for (i=0; i<=SPOKES; i++) {
	    float c, r, b, x, y;
	    
	    c = LAL_TWOPI*i/SPOKES;
	    x = smaj*cos(c);
	    y = smin*sin(c);
	    r = sqrt( x*x + y*y );
	    b = atan2 ( y, x );
	    fprintf( fp, "%e %e\n", alpha + r*cos(angle+b), delta + r*sin(angle+b) );
	  }
	  
	  node = node -> next;
	  
	} /* while node */
      
      LALFree (metricPar.spindown);
      TRY( LALDDestroyVector( stat->statusPtr, &metric ), stat );


    } /* if plotEllipses */
      
  fclose(fp);

  DETATCHSTATUSPTR( stat );
  RETURN (stat);

} /* plotGrid */

/*----------------------------------------------------------------------
 * this is a "wrapper" to provide a uniform interface to PtoleMetric() 
 * and CoherentMetric().
 *
 * the parameter structure of PtoleMetric() was used, because it's more compact
 *----------------------------------------------------------------------*/
/* <lalVerbatim file="DopplerScanCP"> */
void
LALMetricWrapper (LALStatus *stat, 
		  REAL8Vector *metric, 		/* output: metric (pre-allocated)*/
		  const PtoleMetricIn *input) 	/* input-params for the metric */
{ /* </lalVerbatim> */
  MetricParamStruc params = empty_MetricParamStruc;
  PulsarTimesParamStruc spinParams = empty_PulsarTimesParamStruc;
  PulsarTimesParamStruc baryParams = empty_PulsarTimesParamStruc;
  PulsarTimesParamStruc compParams = empty_PulsarTimesParamStruc;
  REAL8Vector *lambda = NULL;
  UINT4 i, nSpin;

  INITSTATUS( stat, "LALMetricWrapper", DOPPLERSCANC );
  ATTATCHSTATUSPTR (stat);

  ASSERT ( input, stat, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL );
  ASSERT ( input->spindown, stat, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL );
  if (input->metricType == LAL_PMETRIC_COH_EPHEM) {
    ASSERT ( input->ephemeris != NULL, stat, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL);
  }

  switch (input->metricType)
    {
    case LAL_PMETRIC_COH_PTOLE_ANALYTIC: /* use Ben&Ian's analytic ptolemaic metric */
      TRY ( LALPtoleMetric (stat->statusPtr, metric, input), stat);
      break;

    case LAL_PMETRIC_COH_PTOLE_NUMERIC:   /* use CoherentMetric + Ptolemaic timing */
    case LAL_PMETRIC_COH_EPHEM:   /* use CoherentMetric + ephemeris timing */
      nSpin = input->spindown->length;

      /* Set up constant parameters for barycentre transformation. */
      baryParams.epoch = input->epoch;
      baryParams.t0 = 0;
      baryParams.latitude = input->site->frDetector.vertexLatitudeRadians;	/* FIXME: should be redundant now, with Detector passed */
      baryParams.longitude = input->site->frDetector.vertexLongitudeRadians;
      baryParams.site = input->site;
      TRY( LALGetEarthTimes( stat->statusPtr, &baryParams ), stat );

      /* set timing-function for earth-motion: either ptolemaic or ephemeris */
      if (input->metricType == LAL_PMETRIC_COH_PTOLE_NUMERIC)
	{
	  baryParams.t1 = LALTBaryPtolemaic;
	  baryParams.dt1 = LALDTBaryPtolemaic;
	}
      else	/* use precise ephemeris-timing */
	{
	  baryParams.t1 = NULL;		/* FIXME: no LALTEphemeris() exists currently (not needed here?) */
	  baryParams.dt1 = LALDTEphemeris;
	  baryParams.ephemeris = input->ephemeris;
	}


      /* Set up input structure for CoherentMetric()  */
      if (nSpin)
	{
	  /* Set up constant parameters for spindown transformation. */
	  spinParams.epoch = input->epoch;
	  spinParams.t0 = 0;
	  
	  /* Set up constant parameters for composed transformation. */
	  compParams.epoch = input->epoch;
	  compParams.t1 = baryParams.t1;
	  compParams.dt1 = baryParams.dt1;
	  compParams.t2 = LALTSpin;
	  compParams.dt2 = LALDTSpin;
	  compParams.constants1 = &baryParams;
	  compParams.constants2 = &spinParams;
	  compParams.nArgs = 2;

	  params.dtCanon = LALDTComp;
	  params.constants = &compParams;
	}
      else	/* simple case: just account for earth motion */
	{
	  params.dtCanon = baryParams.dt1;
	  params.constants = &baryParams;
	}

      params.start = 0;
      params.deltaT = (REAL8) input->duration;
      params.n = 1; 	/* only 1 stack */
      params.errors = 0;

      /* Set up the parameter list. */
      TRY ( LALDCreateVector( stat->statusPtr, &lambda, nSpin + 2 + 1 ), stat );

      lambda->data[0] = (REAL8) input->maxFreq;
      lambda->data[1] = (REAL8) input->position.longitude;	/* alpha */
      lambda->data[2] = (REAL8) input->position.latitude;	/* delta */

      if ( nSpin ) 
	{
	  for (i=0; i < nSpin; i++)
	    lambda->data[3 + i] = (REAL8) input->spindown->data[i];
	}

      /* _finally_ we can call the metric */
      TRY ( LALCoherentMetric( stat->statusPtr, metric, lambda, &params ), stat );
      TRY ( LALDDestroyVector( stat->statusPtr, &lambda ), stat );

      break;

    default:
      LALPrintError ("Unknown metric type `%d`\n", input->metricType);
      ABORT (stat, DOPPLERSCANH_EMETRIC,  DOPPLERSCANH_MSGEMETRIC);
      break;
      
    } /* switch type */

  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* LALMetricWrapper() */


/*----------------------------------------------------------------------
 * function for checking if a given point lies inside or outside a given
 * polygon, which is specified by a list of points in a SkyPositionVector
 *
 * NOTE: the list of points must not close on itself, the last point
 * is automatically assumed to be connected to the first 
 * 
 * Alorithm: count the number of intersections of rays emanating to the right
 * from the point with the lines of the polygon: even=> outside, odd=> inside
 *
 * NOTE2: we try to get this algorith to count all boundary-points as 'inside'
 *     we do this by counting intersection to the left _AND_ to the right
 *     and consider the point inside if either of those says its inside...
 *     to this end we even allow points to lie outside by up to eps~1e-14
 *
 * Return : TRUE or FALSE
 *----------------------------------------------------------------------*/
BOOLEAN
pointInPolygon ( const SkyPosition *point, const SkyRegion *polygon )
{
  UINT4 i;
  UINT4 N;
  UINT4 insideLeft, insideRight;
  BOOLEAN inside = 0;
  SkyPosition *vertex;
  REAL8 xinter, v1x, v1y, v2x, v2y, px, py;

  if (!point || !polygon || !polygon->vertices || (polygon->numVertices < 3) )
    return 0;

  vertex = polygon->vertices;
  N = polygon->numVertices; 	/* num of vertices = num of edges */

  insideLeft = insideRight = 0;

  px = point->longitude;
  py = point->latitude;

  for (i=0; i < N; i++)
    {
      v1x = vertex[i].longitude;
      v1y = vertex[i].latitude;
      v2x = vertex[(i+1) % N].longitude;
      v2y = vertex[(i+1) % N].latitude;

      /* pre-select candidate edges */
      if ( (py <  MIN(v1y,  v2y)) || (py >=  MAX(v1y, v2y) ) || (v1y == v2y) )
	continue;

      /* now calculate the actual intersection point of the horizontal ray with the edge in question*/
      xinter = v1x + (py - v1y) * (v2x - v1x) / (v2y - v1y);

      if (xinter > px)	      /* intersection lies to the right of point */
	insideLeft ++;

      if (xinter < px)       /* intersection lies to the left of point */
	insideRight ++;

    } /* for sides of polygon */

  inside = ( ((insideLeft %2) == 1) || (insideRight %2) == 1);
  return inside;
  
} /* pointInPolygon() */


/*----------------------------------------------------------------------
 * parse a string into a SkyRegion structure: the expected string-format is
 *   " (ra1, dec1), (ra2, dec2), (ra3, dec3), ... "
 *----------------------------------------------------------------------*/
/* <lalVerbatim file="DopplerScanCP"> */
void
ParseSkyRegion (LALStatus *stat, SkyRegion *region, const CHAR *input)
{ /* </lalVerbatim> */
  const CHAR *pos;
  UINT4 i;

  INITSTATUS( stat, "ParseSkyRegion", DOPPLERSCANC );

  ASSERT (region != NULL, stat, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL);
  ASSERT (region->vertices == NULL, stat, DOPPLERSCANH_ENONULL,  DOPPLERSCANH_MSGENONULL);
  ASSERT (input != NULL, stat, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL);

  region->numVertices = 0;
  region->lowerLeft.longitude = LAL_TWOPI;
  region->lowerLeft.latitude = LAL_PI;
  region->upperRight.longitude = 0;
  region->upperRight.latitude = 0;

  /* count number of entries (by # of opening parantheses) */
  pos = input;
  while ( (pos = strchr (pos, '(')) != NULL )
    {
      region->numVertices ++;
      pos ++;
    }

  if (region->numVertices == 0) {
    LALPrintError ("Failed to parse sky-region: `%s`\n", input);
    ABORT (stat, DOPPLERSCANH_ESKYREGION, DOPPLERSCANH_MSGESKYREGION);
  }
    
  
  /* allocate list of vertices */
  if ( (region->vertices = LALMalloc (region->numVertices * sizeof (SkyPosition))) == NULL) {
    ABORT (stat, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM);
  }

  region->lowerLeft.longitude = LAL_TWOPI;
  region->lowerLeft.latitude  = LAL_PI/2.0;

  region->upperRight.longitude = 0;
  region->upperRight.latitude  = -LAL_PI/2;


  /* and parse list of vertices from input-string */
  pos = input;
  for (i = 0; i < region->numVertices; i++)
    {
      if ( sscanf (pos, "(%" LAL_REAL8_FORMAT ", %" LAL_REAL8_FORMAT ")", 
		   &(region->vertices[i].longitude), &(region->vertices[i].latitude) ) != 2) 
	{
	  ABORT (stat, DOPPLERSCANH_ESKYREGION, DOPPLERSCANH_MSGESKYREGION);
	}

      /* keep track of min's and max's to get the bounding square */
      region->lowerLeft.longitude = MIN (region->lowerLeft.longitude, region->vertices[i].longitude);
      region->lowerLeft.latitude  = MIN (region->lowerLeft.latitude, region->vertices[i].latitude);

      region->upperRight.longitude = MAX (region->upperRight.longitude, region->vertices[i].longitude);
      region->upperRight.latitude  =  MAX (region->upperRight.latitude, region->vertices[i].latitude);


      pos = strchr (pos + 1, '(');

    } /* for numVertices */

  RETURN (stat);

} /* ParseSkyRegion() */

/*----------------------------------------------------------------------
 * Translate a TwoDMesh into a DopplerScanGrid using a SkyRegion for clipping
 * 
 * NOTE: the returned grid will be NULL if there are no points inside the sky-region
 *----------------------------------------------------------------------*/
void 
ConvertTwoDMesh2Grid ( LALStatus *stat, 
		       DopplerScanGrid **grid, 	/* output: a dopperScan-grid */
		       const TwoDMeshNode *mesh2d, /* input: a 2Dmesh */
		       const SkyRegion *region )   /* a sky-region for clipping */
{
  const TwoDMeshNode *meshpoint;
  DopplerScanGrid head = empty_grid;
  DopplerScanGrid *node;
  SkyPosition point;

  INITSTATUS( stat, "ConvertTwoDMesh2Grid", DOPPLERSCANC );

  ASSERT ( *grid == NULL, stat, DOPPLERSCANH_ENONULL, DOPPLERSCANH_MSGENONULL );
  ASSERT (region, stat, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL);
  ASSERT (mesh2d, stat, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL);

  meshpoint = mesh2d;	/* this is the 2-d mesh from LALTwoDMesh() */

  node = &head;		/* this is our Doppler-grid (empty for now) */
  
  while (meshpoint)
    {
      if (meshOrder == ORDER_ALPHA_DELTA)
	{
	  point.longitude = meshpoint->x;
	  point.latitude = meshpoint->y;
	}
      else
	{
	  point.longitude = meshpoint->y;
	  point.latitude = meshpoint->x;
	}
      
      if (pointInPolygon (&point, region) )
	{
	  /* prepare a new node for this point */
	  if ( (node->next = (DopplerScanGrid*) LALCalloc (1, sizeof(DopplerScanGrid) )) == NULL) {
	    ABORT ( stat, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM );
	  }
	  node = node->next;

	  node->alpha = point.longitude;
	  node->delta = point.latitude;

	} /* if point in polygon */

      meshpoint = meshpoint->next;

    } /* while meshpoint */

  
  *grid = head.next;		/* return the final grid (excluding static head) */

  RETURN (stat);

} /* ConvertTwoDMesh2Grid() */


/*----------------------------------------------------------------------
 *
 * make a "flat" grid, i.e. a grid with fixed mesh-sizes dAlpha, dDelta
 *
 *----------------------------------------------------------------------*/
void
buildFlatGrid (LALStatus *stat, DopplerScanGrid **grid, const SkyRegion *skyRegion, REAL8 dAlpha, REAL8 dDelta)
{
  SkyPosition thisPoint;
  DopplerScanGrid head = empty_grid;  /* empty head to start grid-list */
  DopplerScanGrid *node = NULL;

  INITSTATUS( stat, "buildFlatGrid", DOPPLERSCANC );

  ASSERT ( grid, stat, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL);
  ASSERT ( skyRegion, stat, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL);
  ASSERT ( (dAlpha > 0) && (dDelta > 0), stat, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT);
  ASSERT ( *grid == NULL, stat, DOPPLERSCANH_ENONULL, DOPPLERSCANH_MSGENONULL);

  thisPoint = skyRegion->lowerLeft;	/* start from lower-left corner */

  node = &head;

  while (1)
    {
      if (pointInPolygon (&thisPoint, skyRegion) )
	{
	  /* prepare this node */
	  node->next = LALCalloc (1, sizeof(DopplerScanGrid));
	  if (node->next == NULL) {
	    freeGrid ( head.next );
	    ABORT (stat, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM);
	  }
	  node = node->next;
	  
	  node->alpha = thisPoint.longitude;
	  node->delta = thisPoint.latitude;
	} /* if pointInPolygon() */
	  
      thisPoint.latitude += dDelta;
	  
      if (thisPoint.latitude > skyRegion->upperRight.latitude)
	{
	  thisPoint.latitude = skyRegion->lowerLeft.latitude;
	  thisPoint.longitude += dAlpha;
	} 
      
      /* this it the break-condition: are we done yet? */
      if (thisPoint.longitude >= skyRegion->upperRight.longitude + dAlpha)
	break;
	  
    } /* while(1) */

  *grid = head.next;	/* return final grid-list */

  
  RETURN (stat);

} /* buildFlatGrid */


/*----------------------------------------------------------------------
 *
 * (approx.) isotropic grid with cells of fixed solid-angle dAlpha x dDelta
 *
 *----------------------------------------------------------------------*/
void
buildIsotropicGrid (LALStatus *stat, DopplerScanGrid **grid, const SkyRegion *skyRegion, REAL8 dAlpha, REAL8 dDelta)
{
  SkyPosition thisPoint;
  DopplerScanGrid head = empty_grid;  /* empty head to start grid-list */
  DopplerScanGrid *node = NULL;
  REAL8 step_alpha, step_delta, cos_delta;

  INITSTATUS( stat, "makeIsotropicGrid", DOPPLERSCANC );

  ASSERT ( grid, stat, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL);
  ASSERT ( skyRegion, stat, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL);
  ASSERT ( *grid == NULL, stat, DOPPLERSCANH_ENONULL, DOPPLERSCANH_MSGENONULL);

  thisPoint = skyRegion->lowerLeft;	/* start from lower-left corner */

  step_delta = dDelta;	/* delta step-size is fixed */
  cos_delta = fabs(cos (thisPoint.latitude));

  node = &head;		/* start our grid with an empty head */

  while (1)
    {
      if (pointInPolygon ( &thisPoint, skyRegion ) ) 
	{
	  /* prepare this node */
	  node->next = LALCalloc (1, sizeof(DopplerScanGrid));
	  if (node->next == NULL) {
	    freeGrid ( head.next );
	    ABORT (stat, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM);
	  }
	  node = node->next;
	      
	  node->alpha = thisPoint.longitude;
	  node->delta = thisPoint.latitude;
	} /* if pointInPolygon() */
      
      step_alpha = dAlpha / cos_delta;	/* alpha stepsize depends on delta */

      thisPoint.longitude += step_alpha;
      if (thisPoint.longitude > skyRegion->upperRight.longitude)
	{
	  thisPoint.longitude = skyRegion->lowerLeft.longitude;
	  thisPoint.latitude += step_delta;	
	  cos_delta = fabs (cos (thisPoint.latitude));
	} 

      /* this it the break-condition: are we done yet? */
      if (thisPoint.latitude > skyRegion->upperRight.latitude)
	break;
	  
    } /* while(1) */

  *grid = head.next;	/* set result: could be NULL! */

  RETURN (stat);

} /* buildIsotropicGrid() */

/*----------------------------------------------------------------------
 *
 * make skygrid using a specified metric
 *
 *----------------------------------------------------------------------*/
void
buildMetricGrid (LALStatus *stat, DopplerScanGrid **grid, SkyRegion *skyRegion,  const DopplerScanInit *init)
{
  SkyPosition thisPoint;
  TwoDMeshNode *mesh2d = NULL;
  TwoDMeshParamStruc meshpar = empty_meshpar;
  PtoleMetricIn metricpar = empty_metricpar;

  INITSTATUS( stat, "makeMetricGrid", DOPPLERSCANC );
  ATTATCHSTATUSPTR (stat);

  ASSERT ( grid, stat, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL);
  ASSERT ( skyRegion, stat, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL);
  ASSERT ( *grid == NULL, stat, DOPPLERSCANH_ENONULL, DOPPLERSCANH_MSGENONULL);
  ASSERT ( (init->metricType > LAL_PMETRIC_NONE) && (init->metricType < LAL_PMETRIC_LAST), 
	   stat, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT);

  thisPoint = skyRegion->lowerLeft;	/* start from lower-left corner */
  
  /* some general mesh-settings are needed in metric case */
  meshpar.getRange = getRange;
      
  if (meshOrder == ORDER_ALPHA_DELTA)
    {
      meshpar.domain[0] = skyRegion->lowerLeft.longitude;
      meshpar.domain[1] = skyRegion->upperRight.longitude;
    }
  else
    {
      meshpar.domain[0] = skyRegion->lowerLeft.latitude;
      meshpar.domain[1] = skyRegion->upperRight.latitude;
    }
  
  meshpar.rangeParams = (void*) skyRegion;

  /* Prepare call of TwoDMesh(): the mesh-parameters */
  meshpar.mThresh = init->metricMismatch;
  meshpar.nIn = 1e6;  	/* maximum nodes in mesh */ /*  FIXME: hardcoded */
      
  /* helper-function: range and metric */
  meshpar.getMetric = getMetric;
  /* and its parameters: metric-parms */
  meshpar.metricParams = (void *) &metricpar;

  /* set up the metric parameters proper (using PtoleMetricIn as container-type) */
  metricpar.metricType = init->metricType;
  metricpar.position.system = COORDINATESYSTEM_EQUATORIAL;
  /* currently, CreateVector's are broken as they don't allow length=0 */
  /*      TRY( LALSCreateVector( stat->statusPtr, &(scan->MetricPar.spindown), 0 ), stat ); 	*/
  /* FIXME: replace when fixed in LAL */
  metricpar.spindown = LALCalloc ( 1, sizeof(REAL4Vector) );
  metricpar.spindown->length=0;
  metricpar.spindown->data=NULL;

  metricpar.epoch = init->obsBegin;
  metricpar.duration = init->obsDuration;
  metricpar.maxFreq = init->fmax;
  metricpar.site = init->Detector;

  /* finally: create the mesh! (ONLY 2D for now!) */
  TRY( LALCreateTwoDMesh( stat->statusPtr, &mesh2d, &meshpar ), stat);

  if (metricpar.spindown) {
    /* FIXME: this is currently broken in LAL, as length=0 is not allowed */
    /*    TRY (LALSDestroyVector ( stat->statusPtr, &(scan->MetricPar.spindown) ), stat); */
    LALFree (metricpar.spindown);
    metricpar.spindown = NULL;
  }

  /* convert this 2D-mesh into our grid-structure, including clipping to the skyRegion */
  ConvertTwoDMesh2Grid ( stat->statusPtr, grid, mesh2d, skyRegion );

  /* get rid of 2D-mesh */
  TRY ( LALDestroyTwoDMesh ( stat->statusPtr,  &mesh2d, 0), stat);

  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* buildMetricGrid() */


/*----------------------------------------------------------------------
 *
 * load skygrid from a file
 *
 *----------------------------------------------------------------------*/
void
loadSkyGridFile (LALStatus *stat, DopplerScanGrid **grid, const CHAR *fname)
{
  LALParsedDataFile *data = NULL;
  DopplerScanGrid *node, head = empty_grid;
  UINT4 i;

  INITSTATUS( stat, "loadGridFile", DOPPLERSCANC );
  ATTATCHSTATUSPTR (stat);

  ASSERT ( grid, stat, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL);
  ASSERT ( *grid == NULL, stat, DOPPLERSCANH_ENONULL, DOPPLERSCANH_MSGENONULL);
  ASSERT ( fname, stat, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL);

  TRY (LALParseDataFile (stat->statusPtr, &data, fname), stat);

  /* parse this list of lines into a sky-grid */
  node = &head;
  for (i=0; i < data->lines->nTokens; i++)
    {
      if ( (node->next = LALCalloc (1, sizeof (DopplerScanGrid))) == NULL)
	{
	  freeGrid (head.next);
	  ABORT (stat, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM);
	}
      if ( 2 != sscanf( data->lines->tokens[i], "%" LAL_REAL8_FORMAT "%" LAL_REAL8_FORMAT, &(node->alpha), &(node->delta)) )
	{
	  LALPrintError ("\nERROR: could not parse line %d in skyGrid-file '%s'\n\n", i, fname);
	  freeGrid (head.next);
	  ABORT (stat, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT);
	}
      
      node = node->next;

    } /* for i < nLines */

  TRY ( LALDestroyParsedDataFile (stat->statusPtr, &data), stat);

  *grid = head.next;	/* pass result */

  DETATCHSTATUSPTR (stat);
  RETURN (stat);
  
} /* loadSkyGridFile() */

/*----------------------------------------------------------------------
 * write a sky-grid to a file, including some comments containing the 
 * parameters of the grid (?)
 *----------------------------------------------------------------------*/
void
writeSkyGridFile (LALStatus *stat, const DopplerScanGrid *grid, const CHAR *fname, const DopplerScanInit *init)
{
  FILE *fp;
  const DopplerScanGrid *node;

  INITSTATUS( stat, "writeGridFile", DOPPLERSCANC );
  ATTATCHSTATUSPTR (stat);

  ASSERT ( grid, stat, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL);
  ASSERT ( fname, stat, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL);
  ASSERT ( init, stat, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL);

  if ( (fp = fopen(fname, "w")) == NULL )
    {
      LALPrintError ("\nERROR: could not open %s for writing!\n", fname);
      ABORT (stat, DOPPLERSCANH_ESYS, DOPPLERSCANH_MSGESYS);
    }

  node = grid;
  while (node)
    {
      fprintf (fp, "%g %g \n", node->alpha, node->delta);
      node = node->next;
    }
  
  fclose (fp);

  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* writeSkyGridFile() */

/*----------------------------------------------------------------------
 * 
 * write predicted frequency-shift of Fmax as function of sky-position
 *
 *----------------------------------------------------------------------*/
void
printFrequencyShifts ( LALStatus *stat, const DopplerScanState *scan, const DopplerScanInit *init)
{
  FILE *fp;
  const CHAR *fname = "dFreq.pred";
  const DopplerScanGrid *node = NULL;

  REAL8 v[3], a[3];
  REAL8 np[3], n[3];
  REAL8 fact, kappa0;
  REAL8 alpha, delta, f0;
  REAL8* vel;
  REAL8* acc;
  REAL8 t0e;        /*time since first entry in Earth ephem. table */  
  INT4 ientryE;      /*entry in look-up table closest to current time, tGPS */

  REAL8 tinitE;           /*time (GPS) of first entry in Earth ephem table*/
  REAL8 tdiffE;           /*current time tGPS minus time of nearest entry
                             in Earth ephem look-up table */
  REAL8 tgps[2];
  EphemerisData *edat = init->ephemeris;
  UINT4 j;
  REAL8 corr0, corr1, accN, accDot[3];
  REAL8 Tobs, dT;
  REAL8 V0[3], V1[3], V2[3];

  INITSTATUS( stat, "printFrequencyShifts", DOPPLERSCANC );
  ATTATCHSTATUSPTR (stat);

  if ( (fp = fopen(fname, "w")) == NULL) {
    LALPrintError ("\nERROR: could not open %s for writing!\n", fname);
    ABORT (stat, DOPPLERSCANH_ESYS, DOPPLERSCANH_MSGESYS);
  }

  /* get detector velocity */
  Tobs = init->obsDuration;

  tgps[0] = (REAL8)(init->obsBegin.gpsSeconds);
  tgps[1] = (REAL8)(init->obsBegin.gpsNanoSeconds);

  tinitE = edat->ephemE[0].gps;
  dT = edat->dtEtable;
  t0e = tgps[0] - tinitE;
  ientryE = floor((t0e/dT) +0.5e0);  /*finding Earth table entry closest to arrival time*/

  /* tdiff is arrival time minus closest Earth table entry; tdiff can be pos. or neg. */
  tdiffE = t0e - edat->dtEtable*ientryE + tgps[1]*1.e-9; 


  vel = edat->ephemE[ientryE].vel;
  acc = edat->ephemE[ientryE].acc; 

  /* interpolate v, a to the actual start-time */
  for (j = 0; j < 3; j++)
    {
      accDot[j] = (edat->ephemE[ientryE+1].acc[j] - edat->ephemE[ientryE].acc[j]) / dT;
      v[j] = vel[j] + acc[j]*tdiffE + 0.5 * accDot[j] * tdiffE * tdiffE;
      a[j] = acc[j] + accDot[j] * tdiffE;
    }


  /* calculate successively higher order-expressions for V */
  for ( j=0; j < 3; j++)
    {
      V0[j] = v[j];			/* 0th order */
      V1[j] = v[j] + 0.5 * a[j]*Tobs;	/* 1st order */
      V2[j] = v[j] + 0.5 * a[j]*Tobs + (2.0/5.0)*accDot[j] * Tobs * Tobs; /* 2nd order */
    }

  printf ("dT = %f, tdiffE = %f, Tobs = %f\n", dT, tdiffE, Tobs);
  printf (" vel =  [ %g, %g, %g ]\n", vel[0], vel[1], vel[2]);
  printf (" acc =  [ %g, %g, %g ]\n", acc[0], acc[1], acc[2]);
  printf (" accDot =  [ %g, %g, %g ]\n\n", accDot[0], accDot[1], accDot[2]);

  printf (" v =  [ %g, %g, %g ]\n", v[0], v[1], v[2]);
  printf (" a =  [ %g, %g, %g ]\n", a[0], a[1], a[2]);

  printf ("\nVelocity-expression in circle-equation: \n");
  printf (" V0 = [ %g, %g, %g ]\n", V0[0], V0[1], V0[2] );
  printf (" V1 = [ %g, %g, %g ]\n", V1[0], V1[1], V1[2] );
  printf (" V2 = [ %g, %g, %g ]\n", V2[0], V2[1], V2[2] );

  node = scan->grid;

  /* signal params */
  alpha = 0.8;
  delta = 1.0;
  f0 = 101.0;

  n[0] = cos(delta) * cos(alpha);
  n[1] = cos(delta) * sin(alpha);
  n[2] = sin(delta);
  
  kappa0 = (1.0 + n[0]*v[0] + n[1]*v[1] + n[2]*v[2] );

  accN = (acc[0]*n[0] + acc[1]*n[1] + acc[2]*n[2]);
  corr1 = (1.0/60.0)* (accN * accN) * Tobs * Tobs / kappa0;

  while (node)
    {
      np[0] = cos(node->delta) * cos(node->alpha);
      np[1] = cos(node->delta) * sin(node->alpha);
      np[2] = sin(node->delta);

      fact = 1.0 / (1.0 + np[0]*v[0] + np[1]*v[1] + np[2]*v[2] );
      
      fprintf (fp, "%.7f %.7f %.7f\n", node->alpha, node->delta, fact);

      node = node->next;
    } /* while grid-points */
  


  fclose (fp);
  
  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* printFrequencyShifts() */
