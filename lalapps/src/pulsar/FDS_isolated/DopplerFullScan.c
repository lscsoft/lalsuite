/*
 * Copyright (C) 2006 Reinhard Prix
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

/**
 * \author Reinhard Prix
 * \date 2006
 * \file
 * \brief Functions for handling "full" multi-dimensional search-grids for CFS.
 *        as opposed to the "factored" grids implemented in DopplerScan.[ch]
 */

/*---------- INCLUDES ----------*/
#include <math.h>

#include <lal/LALStdlib.h>
#include <lal/DetectorSite.h>
#include <lal/LALError.h>
#include <lal/LatticeCovering.h>
#include <lal/LogPrintf.h>

#include "FlatPulsarMetric.h"

#include "DopplerFullScan.h"


/*---------- DEFINES ----------*/
#define MIN(x,y) (x < y ? x : y)
#define MAX(x,y) (x > y ? x : y)

#define TRUE (1==1)
#define FALSE (1==0)

NRCSID( DOPPLERFULLSCANC, "$Id$" );

/*---------- internal types ----------*/

/** ----- internal [opaque] type to store the state of a FULL multidimensional grid-scan ----- */
struct tagDopplerFullScanState {
  INT2 state;  			/**< idle, ready or finished */
  DopplerGridType gridType;	/**< what type of grid are we dealing with */

  PulsarDopplerParams thisPoint; /**< current doppler-position of the scan */
  REAL8 numTemplates;		/**< total number of templates in the grid */

  /* ----- the following are used to emulate FACTORED grids sky x Freq x f1dot ... */
  DopplerSkyScanState skyScan;	/**< keep track of sky-grid stepping */
  PulsarSpinRange spinRange;	/**< spin-range to search */

  /* ----- extensions for full metric lattice-grids  ----- */
  gsl_matrix *gij;		/**< flat parameter-space metric */
  REAL8VectorList *covering;	/**< multi-dimensional covering lattice */
 
} /* struct DopplerFullScanState */;


/*---------- empty initializers ---------- */
/* some empty structs for initializations */
const DopplerFullScanState empty_DopplerFullScanState;

/*---------- Global variables ----------*/
extern INT4 lalDebugLevel;

/*---------- internal prototypes ----------*/
void initLatticeCovering(LALStatus *, DopplerFullScanState *scan, const MultiDetectorStateSeries *mdetStates, const DopplerFullScanInit *init);
void initFactoredGrid (LALStatus *, DopplerFullScanState *scan,	const MultiDetectorStateSeries *mdetStates, const DopplerFullScanInit *init);
int nextPointInFactoredGrid (PulsarDopplerParams *pos, DopplerFullScanState *scan);
int printGSLmatrix ( FILE *fp, const CHAR *fmt, const gsl_matrix *mat );

/*==================== FUNCTION DEFINITIONS ====================*/

/** Set up a full multi-dimensional grid-scan. 
 * Currently this only emulates a 'factored' grid-scan with 'sky x Freq x f1dot ...' , but
 * keeps all details within the DopplerScan module for future extension to real multidimensional
 * grids. 
 * 
 * NOTE: Use 'NextDopplerPos()' to step through this template grid.
 * 
 */
void 
InitDopplerFullScan(LALStatus *status, 
		    DopplerFullScanState **scan,		/**< [out] initialized Doppler scan state */
		    const MultiDetectorStateSeries *mdetStates,	/**< [in] used for list of integration timestamps and detector-info */
		    const DopplerFullScanInit *init		/**< [in] initialization parameters */
		    )
{
  DopplerFullScanState *thisScan;

  INITSTATUS( status, "InitDopplerFullScan", DOPPLERFULLSCANC );
  ATTATCHSTATUSPTR (status); 

  ASSERT ( scan, status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );
  ASSERT ( *scan == NULL, status, DOPPLERSCANH_ENONULL, DOPPLERSCANH_MSGENONULL );  
  ASSERT ( mdetStates, status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );
  ASSERT ( init, status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );

  if ( (thisScan = LALCalloc (1, sizeof(*thisScan) )) == NULL ) {
    ABORT (status, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM);
  }

  /* which "class" of template grid to generate?: factored, or full-multidim ? */
  switch ( init->gridType )
    {
      /* emulate old 'factored' grids 'sky x f0dot x f1dot x f2dot x f3dot': */
    case GRID_FLAT:
    case GRID_ISOTROPIC:
    case GRID_METRIC:
    case GRID_FILE_SKYGRID:
    case GRID_METRIC_SKYFILE:
      /* backwards-compatibility mode */
      TRY ( initFactoredGrid ( status->statusPtr, thisScan, mdetStates, init ), status );
      break;

      /* ----- multi-dimensional coverings based on flat metric approximation ----- */
    case GRID_METRIC_LATTICE:
      /* NOTE: experimental und under construction */
      TRY ( initLatticeCovering ( status->statusPtr, thisScan, mdetStates, init ), status );
      
      break;

    default:
      LALPrintError("\nInvalid grid type '%d'\n\n", init->gridType );
      ABORT ( status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );
      break;
    } /* switch gridType */

  /* return result */
  (*scan) = thisScan;

  DETATCHSTATUSPTR (status);
  RETURN( status );

} /* InitDopplerFullScan() */

/** Initialize Doppler-scanner to emulate an old-style factored template grid: 'sky x f0dot x f1dot x f2dot x f3dot'.
 *  This is a compatiblity-mode with the previous implementation currently also used in ComputeFStatistic.c.
 */
void 
initFactoredGrid (LALStatus *status, 
		  DopplerFullScanState *scan,			/**< [bout] initialized Doppler scan state */
		  const MultiDetectorStateSeries *mdetStates, 	/**< [in] used for list of integration timestamps and detector-info */
		  const DopplerFullScanInit *init		/**< [in] initialization parameters */
		  )
{
  DopplerSkyScanInit skyScanInit = empty_DopplerSkyScanInit;
  SkyPosition skypos;
  UINT4 i;

  INITSTATUS( status, "initFactoredGrid", DOPPLERFULLSCANC );
  ATTATCHSTATUSPTR (status); 

  ASSERT ( scan, status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );
  ASSERT ( scan->state == STATE_IDLE, status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );
  ASSERT ( mdetStates, status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );
  ASSERT ( init, status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );

  /* prepare initialization of DopplerSkyScanner to step through paramter space */
  skyScanInit.dAlpha = init->stepSizes.Alpha;
  skyScanInit.dDelta = init->stepSizes.Delta;
  skyScanInit.gridType = init->gridType;
  skyScanInit.metricType = init->metricType;
  skyScanInit.metricMismatch = init->metricMismatch;
  skyScanInit.projectMetric = TRUE;
  skyScanInit.obsBegin = mdetStates->startTime;
  skyScanInit.obsDuration = mdetStates->Tspan;

  skyScanInit.Detector = &(mdetStates->data[0]->detector);
  skyScanInit.ephemeris = init->ephemeris;		/* used only by Ephemeris-based metric */
  skyScanInit.skyGridFile = init->skyGridFile;
  skyScanInit.skyRegionString = init->searchRegion.skyRegionString;
  skyScanInit.Freq = init->searchRegion.fkdot[0] + init->searchRegion.fkdotBand[0];

  TRY ( InitDopplerSkyScan ( status->statusPtr, &(scan->skyScan), &skyScanInit), status); 

  scan->spinRange.refTime = init->searchRegion.refTime;
  memcpy ( scan->spinRange.fkdot, init->searchRegion.fkdot, sizeof(PulsarSpins) );
  memcpy ( scan->spinRange.fkdotBand, init->searchRegion.fkdotBand, sizeof(PulsarSpins) );

  /* overload spin step-sizes with user-settings if given */
  for (i=0; i < PULSAR_MAX_SPINS; i ++ )
    if ( init->stepSizes.fkdot[i] )
      scan->skyScan.dfkdot[i] = init->stepSizes.fkdot[i];

  /* ----- set Doppler-scanner to start-point ----- */
  scan->thisPoint.refTime = init->searchRegion.refTime;	/* set proper reference time for spins */

  /* normalize skyposition: correctly map into [0,2pi]x[-pi/2,pi/2] */
  skypos.longitude = scan->skyScan.skyNode->Alpha;
  skypos.latitude  = scan->skyScan.skyNode->Delta;
  skypos.system = COORDINATESYSTEM_EQUATORIAL;
  XLALNormalizeSkyPosition ( &skypos );
  scan->thisPoint.Alpha = skypos.longitude;
  scan->thisPoint.Delta = skypos.latitude;
  /* set spins to start */
  for (i=0; i < PULSAR_MAX_SPINS; i ++ )
    scan->thisPoint.fkdot[i] = scan->spinRange.fkdot[i];

  { /* count total number of templates */
    REAL8 nSky, nTot;
    REAL8 nSpins[PULSAR_MAX_SPINS];
    nSky = scan->skyScan.numSkyGridPoints;
    for ( i=0; i < PULSAR_MAX_SPINS; i ++ )
      nSpins[i] = floor( scan->spinRange.fkdotBand[i] / scan->skyScan.dfkdot[i] ) + 1.0;
    nTot = nSky;
    for ( i=0; i < PULSAR_MAX_SPINS; i ++ )
      nTot *= nSpins[i];
    scan->numTemplates = nTot;
    LogPrintf (LOG_DEBUG, "Template grid: nSky x nFreq x nf1dot = %.0f x %.0f x %.0f = %.0f \n", nSky, nSpins[0], nSpins[1], nTot );
  }
  
  /* we're ready */
  scan->state = STATE_READY;

  DETATCHSTATUSPTR (status);
  RETURN( status );

} /* initFactoredGrid() */


/** Return (and compute if not done so yet) the total number of Doppler templates 
 * of the DopplerScan \a scan
 */
REAL8 
XLALNumDopplerTemplates ( DopplerFullScanState *scan)
{ 
  if ( scan->numTemplates )	/* pre-computed already ? */
    return scan->numTemplates;
  else
    {
      /* FIXME: not implemented yet */
      LogPrintf ( LOG_NORMAL, "template counting not implemented yet!\n"); 
      return -1;
    }

} /* XLALNumDopplerTemplates() */

/** Function to step through the full template grid point by point. 
 * Normal return = 0, 
 * errors return -1, 
 * end of scan is signalled by return = 1
 */
int
XLALNextDopplerPos(PulsarDopplerParams *pos, DopplerFullScanState *scan)
{ 

  /* This traps coding errors in the calling routine. */
  if ( pos == NULL || scan == NULL )
    {
      xlalErrno = XLAL_EINVAL;
      return -1;
    }
  if ( scan->state == STATE_IDLE )
    {
      LALPrintError ("\nCalled XLALNextDopplerPos() on un-initialized DopplerFullScanState !\n\n");
      xlalErrno = XLAL_EINVAL;
      return -1;
    }

  /* is this search finished? then return '1' */
  if (  scan->state == STATE_FINISHED ) 
    return 1;

  /* ----- step foward one template in full grid ----- */
  /* Which "class" of template grid are we dealing with: factored, or full-multidim ? */
  switch ( scan->gridType )
    {
      /* emulate old 'factored' grids 'sky x f0dot x f1dot x f2dot x f3dot': */
    case GRID_FLAT:
    case GRID_ISOTROPIC:
    case GRID_METRIC:
    case GRID_FILE_SKYGRID:
    case GRID_METRIC_SKYFILE:
      /* backwards-compatibility mode */
      nextPointInFactoredGrid ( pos, scan );
      break;

    default:
      LALPrintError("\nInvalid grid type '%d'\n\n", scan->gridType );
      xlalErrno = XLAL_EINVAL;
      return -1;
      break;
    } /* switch gridType */

  return 0;
  
} /* XLALNextDopplerPos() */

/** return current grid-point and step forward one template in 'factored' grids (sky x f0dot x f1dot ... )
 * return 0 = OK, -1 = ERROR
*/
int
nextPointInFactoredGrid (PulsarDopplerParams *pos, DopplerFullScanState *scan)
{
  PulsarDopplerParams nextPos = empty_PulsarDopplerParams;
  PulsarSpinRange *range = &(scan->spinRange);	/* shortcut */
  PulsarSpins fkdotMax;
  SkyPosition skypos;

  if ( pos == NULL || scan == NULL )
    return -1;

  (*pos) = scan->thisPoint;	/* RETURN current Doppler-point (struct-copy) */
  
  nextPos = scan->thisPoint;	/* struct-copy: start from current point to get next one */
  
  /* shortcuts: spin boundaries */
  fkdotMax[3] = range->fkdot[3] + range->fkdotBand[3];
  fkdotMax[2] = range->fkdot[2] + range->fkdotBand[2];
  fkdotMax[1] = range->fkdot[1] + range->fkdotBand[1];
  fkdotMax[0] = range->fkdot[0] + range->fkdotBand[0];
  
  /* try to advance to next template */
  nextPos.fkdot[0] += scan->skyScan.dfkdot[0];		/* f0dot one forward */
  if ( nextPos.fkdot[0] >  fkdotMax[0] )		/* too far? */
    {
      nextPos.fkdot[0] = range->fkdot[0];		/* f0dot return to start */
      nextPos.fkdot[1] += scan->skyScan.dfkdot[1];	/* f1dot one step forward */
      if ( nextPos.fkdot[1] > fkdotMax[1] )
	{
	  nextPos.fkdot[1] = range->fkdot[1];		/* f1dot return to start */
	  nextPos.fkdot[2] += scan->skyScan.dfkdot[2];	/* f2dot one forward */
	  if ( nextPos.fkdot[2] > fkdotMax[2] )
	    {
	      nextPos.fkdot[2] = range->fkdot[2];	/* f2dot return to start */
	      nextPos.fkdot[3] += scan->skyScan.dfkdot[3]; /* f3dot one forward */
	      if ( nextPos.fkdot[3] > fkdotMax[3] )
		{
		  nextPos.fkdot[3] = range->fkdot[3];	/* f3dot return to start */
		  					/* skygrid one forward */
		  if ( (scan->skyScan.skyNode = scan->skyScan.skyNode->next) == NULL ) /* no more sky-points ?*/
		    {
		      scan->skyScan.state = STATE_FINISHED;	/* avoid warning when freeing */
		      scan->state = STATE_FINISHED;	/* we're done */
		    } 
		  else
		    {
		      /* normalize next skyposition: correctly map into [0,2pi]x[-pi/2,pi/2] */
		      skypos.longitude = scan->skyScan.skyNode->Alpha; 
		      skypos.latitude  = scan->skyScan.skyNode->Delta;
		      skypos.system = COORDINATESYSTEM_EQUATORIAL;
		      XLALNormalizeSkyPosition ( &skypos );
		      nextPos.Alpha = skypos.longitude;
		      nextPos.Delta = skypos.latitude;      
		    }
		} /* f3dot */
	    } /* f2dot */
	} /* f1dot */
    } /* f0dot */
  
  /* prepare next step */
  scan->thisPoint = nextPos;

  return 0;

} /* nextPointInFactoredGrid() */

/** Destroy the a full DopplerFullScanState structure  
 */
void
FreeDopplerFullScan (LALStatus *status, DopplerFullScanState **scan)
{ 

  INITSTATUS( status, "FreeDopplerFullScan", DOPPLERFULLSCANC);
  ATTATCHSTATUSPTR (status);

  /* This traps coding errors in the calling routine. */
  ASSERT( scan, status, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL );
  ASSERT( *scan, status, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL );
  
  TRY ( FreeDopplerSkyScan ( status->statusPtr, &((*scan)->skyScan) ), status );

  LALFree ( (*scan) );
  (*scan) = NULL;

  DETATCHSTATUSPTR (status);
  RETURN( status );

} /* FreeDopplerSkyScan() */

/** Initialized and construct an optimal lattice-covering for the given searchRegion.
 */
void
initLatticeCovering ( LALStatus *status, 
		     DopplerFullScanState *scan, 
		     const MultiDetectorStateSeries *mdetStates, 
		     const DopplerFullScanInit *init)
{
  UINT4 numSpins, dim;

  INITSTATUS( status, "initLatticeCovering", DOPPLERFULLSCANC );
  ATTATCHSTATUSPTR ( status );

  ASSERT ( scan, status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );
  ASSERT ( scan->state == STATE_IDLE, status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );
  ASSERT ( mdetStates, status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );
  ASSERT ( init, status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );

  /* determine number of spins to compute metric for (at least 1) */
  numSpins = PULSAR_MAX_SPINS;
  while ( (numSpins > 1) && (init->searchRegion.fkdotBand[numSpins - 1] == 0) )
    numSpins --;

  dim = 2 + numSpins;	/* sky + spins (must be at least 3) */

  /* compute flat metric */
  if ( (scan->gij = gsl_matrix_calloc (dim, dim)) == NULL ) {
    ABORT (status, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM);
  }
  if ( XLALFlatMetricCW ( scan->gij, mdetStates, init->searchRegion.refTime ) != 0 ) {
    LALPrintError ("\nCall to XLALFlatMetricCW() failed!\n\n");
    ABORT ( status, DOPPLERSCANH_EXLAL, DOPPLERSCANH_MSGEXLAL );
  }
  
  printGSLmatrix ( stderr, "%.9f ", scan->gij );


  
  DETATCHSTATUSPTR ( status );
  RETURN ( status );

} /* initLatticeCovering() */

/** output gsl-matrix in octave-compatible format to fp, using format \a fmt 
 * for each element 
 */
int
printGSLmatrix ( FILE *fp, const CHAR *fmt, const gsl_matrix *mat )
{
  UINT4 i, j, rows, cols;

  if ( !fp || !fmt || !mat )
    return -1;

  rows = mat->size1;
  cols = mat->size2;

  fprintf (fp,  "[ ");
  for ( i=0; i < rows; i ++ )
    {
      for ( j=0; j < cols; j ++ )
	fprintf ( fp, fmt, gsl_matrix_get (mat, i, j ) );
      if ( i != rows - 1) fprintf (fp, ";\n");
    }
  fprintf (fp, " ];\n");
  
  return 0;

} /* printGSLmatrix() */
