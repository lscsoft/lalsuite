/************************************ <lalVerbatim file="DopplerScanCV">
Author: Prix, Reinhard,  Owen, B. J.,   Jones, D. I.  
$Id$
************************************* </lalVerbatim> */

/* some parts of this file are based on code from PtoleMeshTest, which
   is why I've included Ben & Ian in the Author-list */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{DopplerScan}
\label{ss:DopplerScan.c}

Module to generate subsequent sky-positions.

\subsubsection*{Prototypes}
\input{DopplerScanCP}
\idx{NextSkyPosition()}

\subsubsection*{Description}

This module generates subsequent sky-positions corresponding to a 
given search area and resolution and possibly taking account of the
detector sensitivity pattern (.. in the future).

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

#include "DopplerScan.h"


NRCSID( DOPPLERSCANC, "$Id$" );

/* the SkyScanner can be in one of the following states */
enum {
  STATE_IDLE,   	/* not initialized yet */
  STATE_MANUAL,		/* initialized: 'manual' stepping (i.e. Alpha + i*dAlpha, ...) */
  STATE_GRID,    	/* use pre-calculated template-grid */
  STATE_FINISHED
};


extern INT4 lalDebugLevel;

void getRange( LALStatus *stat, REAL4 y[2], REAL4 x, void *scan );
void getMetric( LALStatus *status, REAL4 g[3], REAL4 skypos[2], void *scan );
void gridStandarizeOrder ( DopplerScanState *scan );

void printGrid (LALStatus *stat, DopplerScanState *scan, BOOLEAN printEllipses);


/*----------------------------------------------------------------------*/
/* <lalVerbatim file="DopplerScanCP"> */
void
InitDopplerScan( LALStatus *stat, DopplerScanState *scan, DopplerScanInit init)
{ /* </lalVerbatim> */

  INITSTATUS( stat, "DopplerScanInit", DOPPLERSCANC );

  ATTATCHSTATUSPTR( stat ); /* prepare for call of LAL-subroutines */

  /* This traps coding errors in the calling routine. */
  ASSERT( scan != NULL, stat, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL );  	/* don't expect output! */

  scan->state = STATE_IDLE;  /* uninitialized */

  scan->Alpha = init.Alpha;
  scan->AlphaBand = init.AlphaBand;
  scan->Delta = init.Delta;
  scan->DeltaBand = init.DeltaBand;

  /* check if there's any scanning of doppler-positions required at all or if it's only one point */
  if ( (scan->AlphaBand == 0.0) && (scan->DeltaBand == 0) )
    { 
      scan->state = STATE_MANUAL;  /* no grid required here */
      
      DETATCHSTATUSPTR (stat);
      RETURN( stat );
    }
  
  /* make sure Bands are always positive (for safty in griding-routines) */
  if ( scan->AlphaBand < 0 )
    {
      scan->Alpha += scan->AlphaBand;
      scan->AlphaBand = - scan->AlphaBand;
    }
  if ( scan->DeltaBand < 0 )
    {
      scan->Delta += scan->DeltaBand;
      scan->DeltaBand = - scan->DeltaBand;
    }

  /* check consistency of bounds for the sky-region */
  /*  FIXME: actually that's nonsense, as the region might stretch across this coordinate "cut" */
  if ( (scan->Alpha < 0) || (scan->Alpha >= LAL_TWOPI) || (fabs(scan->Delta) > LAL_PI_2) 
       ||(scan->Alpha + scan->AlphaBand >= LAL_TWOPI) || (scan->Delta + scan->DeltaBand > LAL_PI_2) )
    {
      ABORT (stat, DOPPLERSCANH_ESKYPARAM, DOPPLERSCANH_MSGESKYPARAM );
    }

  /* init stuff for metric grid stepping */
  scan->grid = NULL;  /* for metric sky-grid: first node */
  scan->gridNode = NULL;

  scan->useMetric = init.useMetric;

  if (init.flipTiling)
    scan->internalOrder = ORDER_ALPHA_DELTA;
  else
    scan->internalOrder = ORDER_DELTA_ALPHA;   /* this is the default for _tiling_ */
  
    
  switch (scan->useMetric)      
    {
    case DOPPLER_MANUAL:
  
      /* stuff for manual stepping */
      scan->dAlpha = fabs(init.dAlpha);
      scan->dDelta = fabs(init.dDelta);
      scan->AlphaCounter = scan->DeltaCounter = 0;

      /* we *need* a 2-dimensional sky-path for now: */
      if (scan->AlphaBand * scan->DeltaBand == 0.0)
	{
	  ABORT( stat, DOPPLERSCANH_E2DSKY, DOPPLERSCANH_MSGE2DSKY );
	}
      
      /* again: we need both or none! */
      if (scan->dAlpha * scan->dDelta == 0.0)
	{
	  ABORT( stat, DOPPLERSCANH_E2DSTEP, DOPPLERSCANH_MSGE2DSTEP);
	}
      
      /* ok: prepare manual sky-stepping */
      scan->AlphaImax = (UINT4)(scan->AlphaBand/scan->dAlpha + 0.5) + 1;
      scan->DeltaImax = (UINT4)(scan->DeltaBand/scan->dDelta + 0.5) + 1;  
      
      scan->state = STATE_MANUAL;
      
      break;
	  
    case DOPPLER_PTOLE_METRIC:
    case DOPPLER_COHERENT_METRIC:

      /* Prepare call of TwoDMesh(): the mesh-parameters */
      scan->meshpar.mThresh = init.metricMismatch;
      scan->meshpar.nIn = 1e6;  	/* maximum nodes in mesh */ /*  FIXME: hardcoded */
      
      /* helper-function: range and metric */
      scan->meshpar.getRange = getRange;
      scan->meshpar.getMetric = getMetric;
      
      /* on the level of DopplerScan we fixed the "standard order"
       * to be ORDER_ALPHA_DELTA, but the meshing+metric functions
       * might need to use the ORDER_DELTA_ALPHA (as is the case for
       * Teviet's TwoDMesh() functions... therefore we have to keep track
       * of the order and put it into normal order after mesh-generation
       * is completed.
       */
      if (scan->internalOrder == ORDER_ALPHA_DELTA)  /* this is the "standard order" */
	{
	  scan->meshpar.domain[0] = scan->Alpha;
	  scan->meshpar.domain[1] = scan->Alpha + scan->AlphaBand;
	}
      else   
	{
	  scan->meshpar.domain[0] = scan->Delta;
	  scan->meshpar.domain[1] = scan->Delta + scan->DeltaBand;
	}
      
      if (scan->useMetric == DOPPLER_PTOLE_METRIC)
	{
	  scan->ptoleMetricPar.position.system = COORDINATESYSTEM_EQUATORIAL;
	  scan->ptoleMetricPar.spindown = NULL;  /*  FIXME: no spindowns for now */
	  scan->ptoleMetricPar.epoch = init.obsBegin;
	  scan->ptoleMetricPar.duration = init.obsDuration;
	  scan->ptoleMetricPar.maxFreq = init.fmax;
	  scan->ptoleMetricPar.site = init.Detector.frDetector;
	}
      else if (scan->useMetric == DOPPLER_COHERENT_METRIC)
	{ /* modelled after StackMetricTest.c */
	  /*  FIXME: currently no spindowns are used */
	  INT2 nSpin = 0;

	  /* Set up start time. */
	  scan->coherentMetricPar.start = 0;
	  scan->coherentMetricPar.deltaT = init.obsDuration;
	  scan->coherentMetricPar.n = 1; /* only 1 stack */
	  scan->coherentMetricPar.errors = 0;

	  /* Set up constant parameters for barycentre transformation. */
	  scan->baryParams.epoch = init.obsBegin;
	  scan->baryParams.latitude = init.Detector.frDetector.vertexLatitudeRadians;
	  scan->baryParams.longitude = init.Detector.frDetector.vertexLongitudeRadians;
	  TRY( LALGetEarthTimes( stat->statusPtr, &(scan->baryParams) ), stat );

	  /* Set up constant parameters for metric calculation. */
	  scan->coherentMetricPar.dtCanon = LALDTBaryPtolemaic;
	  scan->coherentMetricPar.constants = &(scan->baryParams);

	  /* Set up the parameter list. */
	  TRY ( LALDCreateVector( stat->statusPtr, &(scan->lambda), nSpin + 2 + 1 ), stat );

	  scan->lambda->data[0] = init.fmax;
	  scan->lambda->data[1] = 0; /* to be set by grid-function */
	  scan->lambda->data[2] = 0;
	}


      scan->grid = NULL;      
      /* both getMetric and getRange now take the ScanState as parameter */
      scan->meshpar.metricParams = (void *) scan;	 
      scan->meshpar.rangeParams = (void*) scan; 
      
      /* finally: create the mesh! (ONLY 2D for now!) */
      TRY( LALCreateTwoDMesh( stat->statusPtr, &(scan->grid), &(scan->meshpar) ), stat);
      
      /* make sure alpha and delta are in "standard order" in the grid */
      gridStandarizeOrder ( scan );
      
      scan->gridNode = scan->grid; 	/* init to first node in list */
      
      scan->state = STATE_GRID;  /* we are set for a metric-grid scan */
      
      if (lalDebugLevel)
	{
	  TRY( printGrid (stat->statusPtr, scan, 1), stat );
	}
      
      break;

    default:
      ABORT ( stat, DOPPLERSCANH_EMETRIC, DOPPLERSCANH_MSGEMETRIC);
      break;

    } /* switch (useMetric) */


  /* clean up */
  DETATCHSTATUSPTR (stat);

  RETURN( stat );

} /* DopplerScanInit() */


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
    LALWarning (stat, "freeing uninitialized DopplerScan.");
  else if ( scan->state != STATE_FINISHED )
    LALWarning (stat, "freeing unfinished DopplerScan.");

  if ( scan->grid ) {
    TRY (LALDestroyTwoDMesh ( stat->statusPtr,  &(scan->grid), 0), stat);
  }
  scan->grid = scan->gridNode = NULL;

  if (scan->lambda) {
    TRY (LALDDestroyVector ( stat->statusPtr,  &(scan->lambda) ), stat);
  }
  scan->lambda = NULL;
    
  scan->state = STATE_IDLE;

  DETATCHSTATUSPTR (stat);
  RETURN( stat );

} /* FreeDopplerScan() */

/*----------------------------------------------------------------------*/
/* <lalVerbatim file="DopplerScanCP"> */
void
NextDopplerPos( LALStatus *stat, DopplerPosition *pos, DopplerScanState *scan)
{ /* </lalVerbatim> */

  INITSTATUS( stat, "NextSkyPos", DOPPLERSCANC);

  /* This traps coding errors in the calling routine. */
  ASSERT( pos != NULL, stat, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL );


  switch( scan->state )
    {
    case STATE_IDLE:
      ABORT( stat, DOPPLERSCANH_ENOINIT, DOPPLERSCANH_MSGENOINIT );
      break;
      
    case STATE_MANUAL: /* "manual" stepping through templates: Alpha + dAlpha ...*/

      pos->skypos.longitude = scan->Alpha + scan->AlphaCounter * scan->dAlpha;
      pos->skypos.latitude = scan->Delta + scan->DeltaCounter * scan->dDelta;

      pos->spindowns.length = 0; /*  FIXME: not done in here YET */

      /* prepare next step: 'loop' over alpha and delta */ 
      scan->DeltaCounter ++;
      if (scan->DeltaCounter >= scan->DeltaImax)
	{
	  scan->DeltaCounter = 0;
	  scan->AlphaCounter ++;
	  if (scan->AlphaCounter >= scan->AlphaImax)
	    scan->state = STATE_FINISHED;  /* we finish in next call */
	}

      pos->finished = 0; 

      break;

    case STATE_GRID:  /* use pre-determined template-grid */
      /* sanity check: this should NEVER happen! */
      if (scan->gridNode == NULL) {
	ABORT ( stat, DOPPLERSCANH_EGRIDCRPT, DOPPLERSCANH_MSGEGRIDCRPT );
      }

      /* for now we don't use submeshes */
      if (scan->gridNode->subMesh)
	LALWarning ( stat, "The grid contains a sub-mesh: this is currently ignored!" );

      /* NOTE: the grid is assumed to be in "standard order": ORDER_ALPHA_DELTA */
      pos->skypos.longitude = scan->gridNode->x;
      pos->skypos.latitude =  scan->gridNode->y;

      pos->spindowns.length = 0; /*  FIXME */

      /* prepare next step */
      scan->gridNode = scan->gridNode->next;
      if( scan->gridNode == NULL)
	scan->state = STATE_FINISHED; /* finish in next call */

      pos->finished = 0;
      break;

    case STATE_FINISHED:
      pos->finished = 1;  /*  signal to caller that we've finished */
      break;
    }

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
  DopplerScanState *scan = params;
  REAL4 nix;

  /* Set up shop. */
  INITSTATUS( stat, "getRange", DOPPLERSCANC );
  /*   ATTATCHSTATUSPTR( stat ); */

  nix = x;	/* avoid compiler warning */
  
  /* for now: we return the fixed y-range, indendent of x */
  if (scan->internalOrder == ORDER_ALPHA_DELTA)
    {
      y[0] = scan->Delta;
      y[1] = scan->Delta + scan->DeltaBand;
    }
  else
    {
      y[0] = scan->Alpha;
      y[1] = scan->Alpha + scan->AlphaBand;
    }

  /* Clean up and leave. */
  /*   DETATCHSTATUSPTR( stat ); */

  RETURN( stat );
} /* getRange() */


/* ----------------------------------------------------------------------
 * This is a wrapper for the metric function as required by TwoDMesh. 
 *
 * NOTE: at the moment only PtoleMetric() is supported
 *
 * NOTE2: this will be called by TwoDMesh(), therefore
 * skypos is in internalOrder, which is not necessarily ORDER_ALPHA_DELTA!!
 * 
 *----------------------------------------------------------------------*/
void getMetric( LALStatus *stat, REAL4 g[3], REAL4 skypos[2], void *params )
{
  INT2 dim;  	/* dimension of (full) parameter space (incl. freq) */
  DopplerScanState *scan = params;
  REAL8Vector   *metric = NULL;  /* for output of metric */
  REAL4 lon, lat;

  /* Set up shop. */
  INITSTATUS( stat, "getMetric", DOPPLERSCANC );
  ATTATCHSTATUSPTR( stat );


  /* currently we use only f0, alpha, delta: -> 3D metric */
  dim = 3;

  TRY( LALDCreateVector( stat->statusPtr, &metric, dim*(dim+1)/2 ), stat );

  /* Translate input. */
  if (scan->internalOrder == ORDER_DELTA_ALPHA)
    {
      lat = skypos[0];
      lon = skypos[1];
    }
  else
    {
      lon = skypos[0];
      lat = skypos[1];
    }

  /* Call the real metric function. */
  if (scan->useMetric == DOPPLER_PTOLE_METRIC)
    {
      scan->ptoleMetricPar.position.longitude = lon;
      scan->ptoleMetricPar.position.latitude =  lat;
      LALPtoleMetric( stat->statusPtr, metric, &(scan->ptoleMetricPar) );
    }
  else if (scan->useMetric == DOPPLER_COHERENT_METRIC)
    {
      scan->lambda->data[1] = lon;
      scan->lambda->data[2] = lat;
      LALCoherentMetric( stat->statusPtr, metric, scan->lambda, &(scan->coherentMetricPar) );
    }

  BEGINFAIL( stat )
    TRY( LALDDestroyVector( stat->statusPtr, &metric ), stat );
  ENDFAIL( stat );

  LALProjectMetric( stat->statusPtr, metric, 0 );
  BEGINFAIL( stat )
    TRY( LALDDestroyVector( stat->statusPtr, &metric ), stat );
  ENDFAIL( stat );


  /* the general indexing scheme is g_ab for a>=b: index = b + a*(a+1)/2 */
#define INDEX_AA (1 + 1*(1+1)/2)	/* g_aa */
#define INDEX_DD (2 + 2*(2+1)/2)	/* g_dd */
#define INDEX_AD (1 + 2*(2+1)/2)	/* g_ad */

  /* Translate output. Careful about the coordinate-order here!! */
  if (scan->internalOrder == ORDER_ALPHA_DELTA )
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

void printGrid (LALStatus *stat, DopplerScanState *scan, BOOLEAN plotEllipses)
{
  FILE *fp = NULL;
  TwoDMeshNode *node;
  REAL8Vector  *metric = NULL;   /* Parameter-space metric: for plotting ellipses */
  REAL8 gaa, gad, gdd, angle, smaj, smin;
  REAL8 alpha, delta;
  REAL8 mismatch = scan->meshpar.mThresh;
  INT2 set, i;
  INT2 dim;
  CHAR *xmgrHeader = 
    "@version 50103\n"
    "@title \"Sky-grid\"\n"
    "@world xmin -0.1\n"
    "@world xmax 6.4\n"
    "@world ymin -3.2\n"
    "@world ymax 3.2\n"
    "@xaxis label \"alpha\"\n"
    "@yaxis label \"delta\"\n";
  
  /* Set up shop. */
  INITSTATUS( stat, "printGrid", DOPPLERSCANC );
  ATTATCHSTATUSPTR( stat );

  /* currently we use only f0, alpha, delta: -> 3D metric */
  dim = 3;

  if (scan->useMetric == DOPPLER_PTOLE_METRIC)
    fp = fopen ("mesh_ptole.agr", "w");
  else if (scan->useMetric == DOPPLER_COHERENT_METRIC)
    fp = fopen ("mesh_coherent.agr", "w");
  else /* manual stepping: nothing to print here */
    {
      DETATCHSTATUSPTR( stat );
      RETURN (stat);
    }


  if( !fp ) {
    ABORT ( stat, DOPPLERSCANH_ESYS, DOPPLERSCANH_MSGESYS );
  }
  
  fprintf (fp, xmgrHeader);

  /* plot boundary + grid-points using this LAL function */
  TRY (LALXMGRPlotMesh( stat->statusPtr, scan->grid, fp, &scan->meshpar ), stat);

  if (plotEllipses)
    {
      /* make a good guess at which set we're standing at now, say 2 */
      set = 10;
      /* plot ellipses (we need metric for that) */
      TRY( LALDCreateVector( stat->statusPtr, &metric, dim*(dim+1)/2 ), stat);
      node = scan->grid;
      while (node)
	{
	  alpha =  node->x;
	  delta =  node->y;
	  /* Get the metric at this skypos */
	  if (scan->useMetric == DOPPLER_PTOLE_METRIC)
	    {
	      scan->ptoleMetricPar.position.longitude = alpha;
	      scan->ptoleMetricPar.position.latitude  = delta;

	      TRY( LALPtoleMetric( stat->statusPtr, metric, &(scan->ptoleMetricPar) ), stat);

	    }
	  else if (scan->useMetric == DOPPLER_COHERENT_METRIC)
	    {
	      scan->lambda->data[1] = alpha;
	      scan->lambda->data[2] = delta;
	      TRY( LALCoherentMetric( stat->statusPtr, metric, scan->lambda, &(scan->coherentMetricPar) ), stat);
	    }

	  TRY( LALProjectMetric( stat->statusPtr, metric, 0 ), stat);

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

	  printf ("alpha=%f delta=%f\ngaa=%f gdd=%f gad=%f\n", alpha, delta, gaa, gdd, gad);
	  printf ("smaj = %f, smin = %f angle=%f\n", smaj, smin, angle);
	  
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

/* 	/\* Loop around patch ellipse. *\/ */
/* 	  for (i=0; i<=SPOKES; i++) { */
/* 	    float c, r; */
/* 	    c = LAL_TWOPI*i/SPOKES; */
/* 	    r = smaj*smin/sqrt( pow(smaj*sin(c),2) + pow(smin*cos(c),2) ); */
/* 	    fprintf( fp, "%e %e\n", alpha + r*cos(angle-c), delta + r*sin(angle-c) ); */

/*         } /\* for (a...) *\/ */

	  
	  node = node -> next;
	  
	} /* while node */
      
      TRY( LALDDestroyVector( stat->statusPtr, &metric ), stat );

    } /* if plotEllipses */

  /* do some post-processing on the grace-file: if internalOrder is not standardOrder
   * we have to flip the axes of the boundary-drawing ;)
   */
  if (scan->internalOrder != ORDER_ALPHA_DELTA)
    {
      fprintf (fp, "&\n@copy g0.s0 to g0.s1\n");
      fprintf (fp, "@s1.x = s0.y\n");
      fprintf (fp, "@s1.y = s0.x\n");
      fprintf (fp, "@move g0.s1 to g0.s0\n");
    }
      
  fclose(fp);

  DETATCHSTATUSPTR( stat );
  RETURN (stat);

} /* printGrid */

/*-------------------------------------------------------------------------
 * put alpha and delta into "standard order" in the grid: ORDER_ALPHA_DELTA
 *------------------------------------------------------------------------*/
void
gridStandarizeOrder ( DopplerScanState *scan )
{
  REAL4 tmp;
  TwoDMeshNode *node = scan->grid;

  if (scan->internalOrder == ORDER_ALPHA_DELTA)
    return; /* nothing to do here */

  /* otherwise: flip x <--> y */
  while (node)
    {
      tmp = node->x;
      node->x = node->y;
      node->y = tmp;
      
      /* ok, we can't fully preserve the interval-info, as it's incomplete,
	 but I don't see a use for that anyway. Certainly not in ComputeFstatistic
	 or in DopplerScan for that matter...
	 -> we replace dy[0],dy[1] by its mean value
      */
      tmp = (node->dy[1] - node->dy[0])/2.0;
      node->dy[0] = -node->dx;
      node->dy[1] = node->dx;
      node->dx = tmp;

      node = node->next;
    } /* while node */

  return;

}  /* gridStandarizeOrder() */
