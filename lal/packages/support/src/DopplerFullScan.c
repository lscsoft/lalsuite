/*
 * Copyright (C) 2006, 2007 Reinhard Prix, Karl Wette
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
#include <lal/FileIO.h>
#include <lal/DetectorSite.h>
#include <lal/LALError.h>
#include <lal/LatticeCovering.h>
#include <lal/ConfigFile.h>
#include <lal/LogPrintf.h>

#include <lal/FlatPulsarMetric.h>
#include <lal/DopplerFullScan.h>
#include <lal/DopplerLatticeCovering.h>
/* #include <lal/FlatLatticeTiling.h> */

/*---------- DEFINES ----------*/
#define MIN(x,y) (x < y ? x : y)
#define MAX(x,y) (x > y ? x : y)

#define TRUE (1==1)
#define FALSE (1==0)

NRCSID( DOPPLERFULLSCANC, "$Id$" );

/*---------- internal types ----------*/
typedef struct {
  PulsarDopplerParams thisPoint; /**< current doppler-position of the scan */
  DopplerSkyScanState skyScan;	/**< keep track of sky-grid stepping */
  PulsarSpinRange spinRange;	/**< spin-range to search */
} factoredGridScan_t;

/** ----- internal [opaque] type to store the state of a FULL multidimensional grid-scan ----- */
struct tagDopplerFullScanState {
  INT2 state;  			/**< idle, ready or finished */
  DopplerGridType gridType;	/**< what type of grid are we dealing with */
  REAL8 numTemplates;		/**< total number of templates in the grid */

  /* ----- full multi-dim parameter-space grid stuff ----- */
  gsl_matrix *gij;			/**< flat parameter-space metric */
  LIGOTimeGPS refTime;			/**< reference time for grid templates */
  REAL8VectorList *covering;		/**< multi-dimensional covering */
  REAL8VectorList *thisGridPoint; 	/**< pointer to current grid-point */
  /* lattice scan state */
  DopplerLatticeScan *latticeScan;	/**< state of lattice Scan */
  /* spindown lattice tiling */
/*   FlatLatticeTiling *spindownTiling;    /\**< state of spindown lattice tiling *\/ */

  /* ----- emulate old-style factored grids */
  factoredGridScan_t *factoredScan;	/**< only used to emulate FACTORED grids sky x Freq x f1dot */

} /* struct DopplerFullScanState */;


/*---------- empty initializers ---------- */
/* some empty structs for initializations */
const DopplerFullScanState empty_DopplerFullScanState;

/*---------- Global variables ----------*/
extern INT4 lalDebugLevel;

/*---------- internal prototypes ----------*/
void initFactoredGrid (LALStatus *, DopplerFullScanState *scan,	const DopplerFullScanInit *init);
int nextPointInFactoredGrid (PulsarDopplerParams *pos, DopplerFullScanState *scan);

/*==================== FUNCTION DEFINITIONS ====================*/

/** Set up a full multi-dimensional grid-scan.
 * Currently this only emulates a 'factored' grid-scan with 'sky x Freq x f1dot ...' , but
 * keeps all details within the DopplerScan module for future extension to real multidimensional
 * grids.
 *
 * NOTE: Use 'XLALNextDopplerPos()' to step through this template grid.
 *
 */
void
InitDopplerFullScan(LALStatus *status,
		    DopplerFullScanState **scan,		/**< [out] initialized Doppler scan state */
		    const DopplerFullScanInit *init		/**< [in] initialization parameters */
		    )
{
  DopplerFullScanState *thisScan;

  INITSTATUS( status, "InitDopplerFullScan", DOPPLERFULLSCANC );
  ATTATCHSTATUSPTR (status);

  ASSERT ( scan, status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );
  ASSERT ( *scan == NULL, status, DOPPLERSCANH_ENONULL, DOPPLERSCANH_MSGENONULL );
  ASSERT ( init, status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );

  if ( (thisScan = LALCalloc (1, sizeof(*thisScan) )) == NULL ) {
    ABORT (status, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM);
  }

  thisScan->gridType = init->gridType;

  /* which "class" of template grid to generate?: factored, or full-multidim ? */
  switch ( thisScan->gridType )
    {
      /* emulate old 'factored' grids 'sky x f0dot x f1dot x f2dot x f3dot': */
    case GRID_FLAT:
    case GRID_ISOTROPIC:
    case GRID_METRIC:
    case GRID_FILE_SKYGRID:
    case GRID_METRIC_SKYFILE:
      /* backwards-compatibility mode */
      TRY ( initFactoredGrid ( status->statusPtr, thisScan, init ), status );
      break;

      /* ----- multi-dimensional covering of full parameter space ----- */
    case GRID_FILE_FULLGRID:
      TRY ( loadFullGridFile ( status->statusPtr, thisScan, init ), status );
      thisScan->refTime = init->startTime;
      break;

    case GRID_METRIC_LATTICE:
      /* NOTE: experimental und under construction */
      {
	DopplerLatticeInit latticeInit;
	latticeInit.searchRegion = init->searchRegion;
	latticeInit.metricMismatch = init->metricMismatch;
	latticeInit.startTime = init->startTime;
	latticeInit.Tspan = init->Tspan;
	latticeInit.ephemeris = init->ephemeris;

	TRY ( InitDopplerLatticeScan ( status->statusPtr, &(thisScan->latticeScan), &latticeInit ), status );
	thisScan->state = STATE_READY;
      }

      break;

    case GRID_SPINDOWN_LATTICE: /* spindown lattice tiling */
      {
      }

      break;

    default:
      LALPrintError("\nInvalid grid type '%d'\n\n", init->gridType );
      ABORT ( status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );
      break;
    } /* switch gridType */

  /* we're ready */
  thisScan->state = STATE_READY;

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
		  const DopplerFullScanInit *init		/**< [in] initialization parameters */
		  )
{
  DopplerSkyScanInit skyScanInit = empty_DopplerSkyScanInit;
  SkyPosition skypos;
  factoredGridScan_t *fscan = NULL;
  UINT4 i;

  INITSTATUS( status, "initFactoredGrid", DOPPLERFULLSCANC );
  ATTATCHSTATUSPTR (status);

  ASSERT ( scan, status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );
  ASSERT ( scan->state == STATE_IDLE, status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );
  ASSERT ( init, status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );

  /* prepare initialization of DopplerSkyScanner to step through paramter space */
  skyScanInit.dAlpha = init->stepSizes.Alpha;
  skyScanInit.dDelta = init->stepSizes.Delta;
  skyScanInit.gridType = init->gridType;
  skyScanInit.metricType = init->metricType;
  skyScanInit.metricMismatch = init->metricMismatch;
  skyScanInit.projectMetric = init->projectMetric;
  skyScanInit.obsBegin = init->startTime;
  skyScanInit.obsDuration = init->Tspan;

  skyScanInit.Detector = init->Detector;
  skyScanInit.ephemeris = init->ephemeris;		/* used only by Ephemeris-based metric */
  skyScanInit.skyGridFile = init->gridFile;
  skyScanInit.skyRegionString = init->searchRegion.skyRegionString;
  skyScanInit.Freq = init->searchRegion.fkdot[0] + init->searchRegion.fkdotBand[0];

  if ( (fscan = LALCalloc ( 1, sizeof( *fscan ))) == NULL ) {
    ABORT ( status, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM );
  }
  scan->factoredScan = fscan;
  TRY ( InitDopplerSkyScan ( status->statusPtr, &(fscan->skyScan), &skyScanInit), status);

  fscan->spinRange.refTime = init->searchRegion.refTime;
  memcpy ( fscan->spinRange.fkdot, init->searchRegion.fkdot, sizeof(PulsarSpins) );
  memcpy ( fscan->spinRange.fkdotBand, init->searchRegion.fkdotBand, sizeof(PulsarSpins) );

  /* overload spin step-sizes with user-settings if given */
  for (i=0; i < PULSAR_MAX_SPINS; i ++ )
    if ( init->stepSizes.fkdot[i] )
      fscan->skyScan.dfkdot[i] = init->stepSizes.fkdot[i];

  /* ----- set Doppler-scanner to start-point ----- */
  fscan->thisPoint.refTime = init->searchRegion.refTime;	/* set proper reference time for spins */

  /* normalize skyposition: correctly map into [0,2pi]x[-pi/2,pi/2] */
  skypos.longitude = fscan->skyScan.skyNode->Alpha;
  skypos.latitude  = fscan->skyScan.skyNode->Delta;
  skypos.system = COORDINATESYSTEM_EQUATORIAL;
  XLALNormalizeSkyPosition ( &skypos );
  fscan->thisPoint.Alpha = skypos.longitude;
  fscan->thisPoint.Delta = skypos.latitude;
  /* set spins to start */
  for (i=0; i < PULSAR_MAX_SPINS; i ++ )
    fscan->thisPoint.fkdot[i] = fscan->spinRange.fkdot[i];

  { /* count total number of templates */
    REAL8 nSky, nTot;
    REAL8 nSpins[PULSAR_MAX_SPINS];
    nSky = fscan->skyScan.numSkyGridPoints;
    for ( i=0; i < PULSAR_MAX_SPINS; i ++ )
      nSpins[i] = floor( fscan->spinRange.fkdotBand[i] / fscan->skyScan.dfkdot[i] ) + 1.0;
    nTot = nSky;
    for ( i=0; i < PULSAR_MAX_SPINS; i ++ )
      nTot *= nSpins[i];
    scan->numTemplates = nTot;
    LogPrintf (LOG_DEBUG, "Template grid: nSky x nFreq x nf1dot = %.0f x %.0f x %.0f = %.0f \n", nSky, nSpins[0], nSpins[1], nTot );
  }

  DETATCHSTATUSPTR (status);
  RETURN( status );

} /* initFactoredGrid() */


/** Return (and compute if not done so yet) the total number of Doppler templates
 * of the DopplerScan \a scan
 */
REAL8
XLALNumDopplerTemplates ( DopplerFullScanState *scan)
{
  if ( ! scan->numTemplates )	/* not pre-computed already ? */
    {
      switch ( scan->gridType )
	{
	case GRID_METRIC_LATTICE:
	  LogPrintf ( LOG_DEBUG, "Now counting number of templates in lattice ... ");
	  scan->numTemplates = XLALCountLatticeTemplates ( scan->latticeScan );
	  LogPrintfVerbatim( LOG_DEBUG, " done. (%.0f)\n", scan->numTemplates );
	  break;

	case GRID_SPINDOWN_LATTICE: /* spindown lattice tiling */
	  break;

	default:
	  /* FIXME: not implemented yet */
	  LogPrintf ( LOG_NORMAL, "template counting not implemented yet for gridType=%d!\n", scan->gridType );
	  return -1;
	  break;
	} /* switch() */
    } /* ! numTemplates */

  return scan->numTemplates;

} /* XLALNumDopplerTemplates() */

/** Function to step through the full template grid point by point.
 * Normal return = 0,
 * errors return -1,
 * end of scan is signalled by return = 1
 */
int
XLALNextDopplerPos(PulsarDopplerParams *pos, DopplerFullScanState *scan)
{
  int ret;
  const CHAR *fn = "XLALNextDopplerPos()";

  /* This traps coding errors in the calling routine. */
  if ( pos == NULL || scan == NULL ) {
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }
  if ( scan->state == STATE_IDLE ) {
    LALPrintError ("\nCalled XLALNextDopplerPos() on un-initialized DopplerFullScanState !\n\n");
    XLAL_ERROR ( fn, XLAL_EINVAL );
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

    case GRID_FILE_FULLGRID:
      pos->refTime = scan->refTime;
      pos->fkdot[0]	= scan->thisGridPoint->entry.data[0];
      pos->Alpha 	= scan->thisGridPoint->entry.data[1];
      pos->Delta 	= scan->thisGridPoint->entry.data[2];
      pos->fkdot[1] 	= scan->thisGridPoint->entry.data[3];
      pos->fkdot[2] 	= scan->thisGridPoint->entry.data[4];
      pos->fkdot[3] 	= scan->thisGridPoint->entry.data[5];
      pos->orbit = NULL;
      /* advance to next grid point */
      if ( ( scan->thisGridPoint = scan->thisGridPoint->next ) == NULL )
	scan->state = STATE_FINISHED;

      break;

    case GRID_METRIC_LATTICE:
      if ( XLALgetCurrentDopplerPos ( pos, scan->latticeScan, COORDINATESYSTEM_EQUATORIAL ) ) {
	XLAL_ERROR ( fn, XLAL_EFUNC );
      }
      /* advance to next point */
      ret = XLALadvanceLatticeIndex ( scan->latticeScan );
      if ( ret < 0 ) {
	XLAL_ERROR ( fn, XLAL_EFUNC );
      }
      else if ( ret == 1 )
	{
	  LALPrintError ( "\n\nXLALadvanceLatticeIndex(): this was the last lattice points!\n\n");
	  scan->state = STATE_FINISHED;
	}
#if 0
      { /* debugging */
	gsl_vector_int *index = NULL;
	XLALgetCurrentLatticeIndex ( &index, scan->latticeScan );
	XLALfprintfGSLvector_int ( stderr, "%d", index );
	gsl_vector_int_free ( index );
      }
#endif

      break;

    case GRID_SPINDOWN_LATTICE: /* spindown lattice tiling */
      {
      }

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
  factoredGridScan_t *fscan;
  PulsarSpinRange *range;
  PulsarSpins fkdotMax;
  SkyPosition skypos;

  if ( pos == NULL || scan == NULL )
    return -1;

  if ( ( fscan = scan->factoredScan ) == NULL )
    return -1;

  range = &(fscan->spinRange);	/* shortcut */

  (*pos) = fscan->thisPoint;	/* RETURN current Doppler-point (struct-copy) */

  nextPos = fscan->thisPoint;	/* struct-copy: start from current point to get next one */

  /* shortcuts: spin boundaries */
  fkdotMax[3] = range->fkdot[3] + range->fkdotBand[3];
  fkdotMax[2] = range->fkdot[2] + range->fkdotBand[2];
  fkdotMax[1] = range->fkdot[1] + range->fkdotBand[1];
  fkdotMax[0] = range->fkdot[0] + range->fkdotBand[0];

  /* try to advance to next template */
  nextPos.fkdot[0] += fscan->skyScan.dfkdot[0];		/* f0dot one forward */
  if ( nextPos.fkdot[0] >  fkdotMax[0] )		/* too far? */
    {
      nextPos.fkdot[0] = range->fkdot[0];		/* f0dot return to start */
      nextPos.fkdot[1] += fscan->skyScan.dfkdot[1];	/* f1dot one step forward */
      if ( nextPos.fkdot[1] > fkdotMax[1] )
	{
	  nextPos.fkdot[1] = range->fkdot[1];		/* f1dot return to start */
	  nextPos.fkdot[2] += fscan->skyScan.dfkdot[2];	/* f2dot one forward */
	  if ( nextPos.fkdot[2] > fkdotMax[2] )
	    {
	      nextPos.fkdot[2] = range->fkdot[2];	/* f2dot return to start */
	      nextPos.fkdot[3] += fscan->skyScan.dfkdot[3]; /* f3dot one forward */
	      if ( nextPos.fkdot[3] > fkdotMax[3] )
		{
		  nextPos.fkdot[3] = range->fkdot[3];	/* f3dot return to start */
		  					/* skygrid one forward */
		  if ( (fscan->skyScan.skyNode = fscan->skyScan.skyNode->next) == NULL ) /* no more sky-points ?*/
		    {
		      fscan->skyScan.state = STATE_FINISHED;	/* avoid warning when freeing */
		      scan->state = STATE_FINISHED;	/* we're done */
		    }
		  else
		    {
		      /* normalize next skyposition: correctly map into [0,2pi]x[-pi/2,pi/2] */
		      skypos.longitude = fscan->skyScan.skyNode->Alpha;
		      skypos.latitude  = fscan->skyScan.skyNode->Delta;
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
  fscan->thisPoint = nextPos;

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

  if ( (*scan)->factoredScan ) {
    TRY ( FreeDopplerSkyScan ( status->statusPtr, &((*scan)->factoredScan->skyScan) ), status );
    LALFree ( (*scan)->factoredScan );
  }

  if ( (*scan)->covering )
    XLALREAL8VectorListDestroy ( (*scan)->covering );

  if ( (*scan)->latticeScan ) {
    XLALFreeDopplerLatticeScan ( &((*scan)->latticeScan) );
  }

/*   if ((*scan)->spindownTiling) { */
/*     XLALDestroyFlatLatticeTiling((*scan)->spindownTiling); */
/*   } */

  LALFree ( (*scan) );
  (*scan) = NULL;

  DETATCHSTATUSPTR (status);
  RETURN( status );

} /* FreeDopplerSkyScan() */



/** load a full multi-dim template grid from the file init->gridFile
 */
void
loadFullGridFile ( LALStatus *status,
		   DopplerFullScanState *scan,
		   const DopplerFullScanInit *init
		   )
{
  REAL8VectorList head = empty_REAL8VectorList;
  REAL8VectorList *tail = NULL;
  REAL8Vector *entry = NULL;
  UINT4 numTemplates, i;
  FILE *fp;

  INITSTATUS( status, "loadFullGridFile", DOPPLERFULLSCANC );
  ATTATCHSTATUSPTR ( status );

  ASSERT ( scan, status, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL);
  ASSERT ( init, status, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL);
  ASSERT ( init->gridFile, status, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL);
  ASSERT ( scan->state == STATE_IDLE, status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT);

  if ( (fp = LALOpenDataFile (init->gridFile)) == NULL) {
    LALPrintError ("Could not open data-file: `%s`\n\n", init->gridFile);
    ABORT (status, CONFIGFILEH_EFILE, CONFIGFILEH_MSGEFILE);
  }

  if ( (entry = XLALCreateREAL8Vector ( 6 ) ) == NULL ) {
    ABORT (status, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM);
  }

  /* parse this list of lines into a full grid */
  numTemplates = 0;
  tail = &head;	/* head will remain empty! */
  while ( ! feof ( fp ) )
    {
      if ( 6 != fscanf( fp, "%lf %lf %lf %lf %lf %lf\n",
			entry->data + 0, entry->data + 1, entry->data + 2, entry->data + 3, entry->data + 4, entry->data + 5 ) )
	{
	  LogPrintf (LOG_CRITICAL,"ERROR: Failed to parse 6 REAL's from line %d in grid-file '%s'\n\n", i, init->gridFile);
	  if ( head.next )
	    XLALREAL8VectorListDestroy (head.next);
	  ABORT (status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT);
	}

      if ( (tail = XLALREAL8VectorListAddEntry (tail, entry)) == NULL )
	{
	  if ( head.next )
	    XLALREAL8VectorListDestroy (head.next);
	  ABORT ( status, DOPPLERSCANH_EXLAL, DOPPLERSCANH_MSGEXLAL );
	}

      numTemplates ++ ;

    } /* while !feof(fp) */

  LALFree ( entry->data );
  LALFree ( entry );

  LogPrintf (LOG_DEBUG, "Template grid: nTot = %.0f\n", 1.0 * numTemplates );
  /* return ready scan-state  */
  scan->numTemplates = numTemplates;
  scan->covering = head.next;	/* pass result (without head!) */
  scan->thisGridPoint = scan->covering;	/* init to start */

  DETATCHSTATUSPTR ( status );
  RETURN ( status );

} /* loadFullGridFile() */
