/*
 *  Copyright (C) 2007, 2008 Karl Wette
 *  Copyright (C) 2006, 2007 Reinhard Prix
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
 * \author Reinhard Prix, Karl Wette
 * \date 2006, 2007, 2008
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
#include <lal/FlatLatticeTiling.h>
#include <lal/FlatLatticeTilingPulsar.h>

/*---------- DEFINES ----------*/
#define MIN(x,y) (x < y ? x : y)
#define MAX(x,y) (x > y ? x : y)

#define TRUE (1==1)
#define FALSE (1==0)

#define INIT_MEM(x) memset(&(x), 0, sizeof((x)))

/*---------- internal types ----------*/
typedef struct {
  PulsarDopplerParams thisPoint; /**< current doppler-position of the scan */
  DopplerSkyScanState skyScan;	/**< keep track of sky-grid stepping */
} factoredGridScan_t;

/** ----- internal [opaque] type to store the state of a FULL multidimensional grid-scan ----- */
struct tagDopplerFullScanState {
  INT2 state;  			/**< idle, ready or finished */
  DopplerGridType gridType;	/**< what type of grid are we dealing with */
  REAL8 numTemplates;		/**< total number of templates in the grid */
  PulsarSpinRange spinRange;	/**< spin-range covered by template bank */
  SkyRegion skyRegion;		/**< sky-range covered by template bank */

  /* ----- full multi-dim parameter-space grid stuff ----- */
  gsl_matrix *gij;			/**< flat parameter-space metric */
  REAL8VectorList *covering;		/**< multi-dimensional covering */
  REAL8VectorList *thisGridPoint; 	/**< pointer to current grid-point */
  /* lattice scan state */
  DopplerLatticeScan *latticeScan;	/**< state of lattice Scan */
  /* spindown lattice tiling */
  FlatLatticeTiling *spindownTiling;    /**< state of spindown lattice tiling */

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
InitDopplerFullScan(LALStatus *status,			/**< pointer to LALStatus structure */
		    DopplerFullScanState **scan,	/**< [out] initialized Doppler scan state */
		    const DopplerFullScanInit *init	/**< [in] initialization parameters */
		    )
{
  DopplerFullScanState *thisScan;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT ( scan, status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );
  ASSERT ( *scan == NULL, status, DOPPLERSCANH_ENONULL, DOPPLERSCANH_MSGENONULL );
  ASSERT ( init, status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );

  if ( (thisScan = LALCalloc (1, sizeof(*thisScan) )) == NULL ) {
    ABORT (status, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM);
  }

  thisScan->gridType = init->gridType;

  /* store the user-input spinRange (includes refTime) in DopplerFullScanState */
  thisScan->spinRange.refTime = init->searchRegion.refTime;
  memcpy ( thisScan->spinRange.fkdot, init->searchRegion.fkdot, sizeof(PulsarSpins) );
  memcpy ( thisScan->spinRange.fkdotBand, init->searchRegion.fkdotBand, sizeof(PulsarSpins) );

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

    case GRID_SPINDOWN_SQUARE: /* square parameter space */
    case GRID_SPINDOWN_AGEBRK: /* age-braking index parameter space */
      {

	int i;
	SkyRegion sky = empty_SkyRegion;

	/* Check that the reference time is the same as the start time */
	if (XLALGPSCmp(&thisScan->spinRange.refTime, &init->startTime) != 0) {
	  XLALPrintError("\nGRID_SPINDOWN_{SQUARE,AGEBRK}: This option currently restricts "
			"the reference time to be the same as the start time.\n");
	  ABORT(status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT);
	}

	/* Create a flat lattice tiling */
	if (NULL == (thisScan->spindownTiling = XLALCreateFlatLatticeTiling(2 + PULSAR_MAX_SPINS))) {
	  XLALPrintError("\nGRID_SPINDOWN_{SQUARE,AGEBRK}: XLALCreateFlatLatticeTiling failed\n");
	  ABORT(status, DOPPLERSCANH_EXLAL, DOPPLERSCANH_MSGEXLAL);
	}

	/* Parse the sky region string and check that it consists of only one point, and set bounds on it */
 	TRY(ParseSkyRegionString(status->statusPtr, &sky, init->searchRegion.skyRegionString), status);
	if (sky.numVertices != 1) {
	  XLALPrintError("\nGRID_SPINDOWN_{SQUARE,AGEBRK}: This option can only handle a single sky position.\n");
	  ABORT(status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT);
	}
	if (sky.vertices[0].system != COORDINATESYSTEM_EQUATORIAL) {
 	  XLALPrintError("\nGRID_SPINDOWN_{SQUARE,AGEBRK}: This option only understands COORDINATESYSTEM_EQUATORIAL\n");
	  ABORT(status, DOPPLERSCANH_ESKYREGION, DOPPLERSCANH_MSGESKYREGION);
	}
	if (XLAL_SUCCESS != XLALSetFlatLatticeConstantBound(thisScan->spindownTiling, 0, sky.vertices[0].longitude, sky.vertices[0].longitude)) {
          XLALPrintError("\nGRID_SPINDOWN_{SQUARE,AGEBRK}: XLALSetFlatLatticeTilingConstantBound failed\n");
	  ABORT(status, DOPPLERSCANH_EXLAL, DOPPLERSCANH_MSGEXLAL);
	}
	if (XLAL_SUCCESS != XLALSetFlatLatticeConstantBound(thisScan->spindownTiling, 1, sky.vertices[0].latitude, sky.vertices[0].latitude)) {
          XLALPrintError("\nGRID_SPINDOWN_{SQUARE,AGEBRK}: XLALSetFlatLatticeTilingConstantBound failed\n");
	  ABORT(status, DOPPLERSCANH_EXLAL, DOPPLERSCANH_MSGEXLAL);
	}
	if (sky.vertices)
	  LALFree (sky.vertices);

	/* Setup parameter space */
	if (thisScan->gridType == GRID_SPINDOWN_SQUARE) { /* square parameter space */

	  /* Set square bounds on the frequency and spindowns */
	  for (i = 0; i < PULSAR_MAX_SPINS; ++i) {
	    if (XLAL_SUCCESS != XLALSetFlatLatticeConstantBound(thisScan->spindownTiling, 2 + i, init->searchRegion.fkdot[i],
                                                                init->searchRegion.fkdot[i] + init->searchRegion.fkdotBand[i])) {
	      XLALPrintError("\nGRID_SPINDOWN_SQUARE: XLALSetFlatLatticeTilingConstantBound failed\n");
	      ABORT(status, DOPPLERSCANH_EXLAL, DOPPLERSCANH_MSGEXLAL);
	    }
	  }

	}
	else if (thisScan->gridType == GRID_SPINDOWN_AGEBRK) { /* age-braking index parameter space */

	  /* Get age and braking index from extra arguments */
	  const REAL8 spindownAge = init->extraArgs[0];
	  const REAL8 minBraking = init->extraArgs[1];
	  const REAL8 maxBraking = init->extraArgs[2];

	  /* Set age-braking index parameter space */
          if (XLAL_SUCCESS != XLALSetFlatLatticeConstantBound(thisScan->spindownTiling, 2, init->searchRegion.fkdot[0],
                                                              init->searchRegion.fkdot[0] + init->searchRegion.fkdotBand[0])) {
            XLALPrintError("\nGRID_SPINDOWN_AGEBRK: XLALSetFlatLatticeTilingConstantBound failed\n");
            ABORT(status, DOPPLERSCANH_EXLAL, DOPPLERSCANH_MSGEXLAL);
          }
	  if (XLAL_SUCCESS != XLALSetFlatLatticeF1DotAgeBrakingBound(thisScan->spindownTiling, 2, 3, spindownAge, minBraking, maxBraking)) {
	    XLALPrintError("\nGRID_SPINDOWN_AGEBRK: XLALSetFlatLatticeF1DotAgeBrakingBound failed\n");
	    ABORT(status, DOPPLERSCANH_EXLAL, DOPPLERSCANH_MSGEXLAL);
	  }
	  if (XLAL_SUCCESS != XLALSetFlatLatticeF2DotBrakingBound(thisScan->spindownTiling, 2, 3, 4, minBraking, maxBraking)) {
	    XLALPrintError("\nGRID_SPINDOWN_AGEBRK: XLALSetFlatLatticeF2DotBrakingBound failed\n");
	    ABORT(status, DOPPLERSCANH_EXLAL, DOPPLERSCANH_MSGEXLAL);
	  }

	  /* This current only goes up to second spindown, so bound higher dimensions */
	  for (i = 3; i < PULSAR_MAX_SPINS; ++i) {
	    if (XLAL_SUCCESS != XLALSetFlatLatticeConstantBound(thisScan->spindownTiling, 2 + i, init->searchRegion.fkdot[i],
                                                                init->searchRegion.fkdot[i] + init->searchRegion.fkdotBand[i])) {
	      XLALPrintError("\nGRID_SPINDOWN_AGEBRK: XLALSetFlatLatticeTilingConstantBound failed\n");
	      ABORT(status, DOPPLERSCANH_EXLAL, DOPPLERSCANH_MSGEXLAL);
	    }
	  }

	}

	/* Set Anstar lattice */
	if (XLAL_SUCCESS != XLALSetFlatLatticeGenerator(thisScan->spindownTiling, XLALAnstarLatticeGenerator)) {
	  XLALPrintError("\nGRID_SPINDOWN_{SQUARE,AGEBRK}: XLALSetFlatTilingAnstarLattice failed\n");
	  ABORT(status, DOPPLERSCANH_EXLAL, DOPPLERSCANH_MSGEXLAL);
	}

	/* Set spindown metric */
        gsl_matrix* metric = gsl_matrix_alloc(2 + PULSAR_MAX_SPINS, 2 + PULSAR_MAX_SPINS);
        if (metric == NULL) {
	  XLALPrintError("\nGRID_SPINDOWN_{SQUARE,AGEBRK}: metric == NULL\n");
	  ABORT(status, DOPPLERSCANH_EXLAL, DOPPLERSCANH_MSGEXLAL);
        }
        gsl_matrix_set_identity(metric);
        gsl_matrix_view spin_metric = gsl_matrix_submatrix(metric, 2, 2, PULSAR_MAX_SPINS, PULSAR_MAX_SPINS);
        if (XLALSpindownMetric(&spin_metric.matrix, init->Tspan) != XLAL_SUCCESS) {
	  XLALPrintError("\nGRID_SPINDOWN_{SQUARE,AGEBRK}: XLALSpindownMetric() failed\n");
	  ABORT(status, DOPPLERSCANH_EXLAL, DOPPLERSCANH_MSGEXLAL);
        }
        if (XLALSetFlatLatticeMetric(thisScan->spindownTiling, metric, init->metricMismatch) != XLAL_SUCCESS) {
	  XLALPrintError("\nGRID_SPINDOWN_{SQUARE,AGEBRK}: XLALSetFlatLatticeMetric() failed\n");
	  ABORT(status, DOPPLERSCANH_EXLAL, DOPPLERSCANH_MSGEXLAL);
        }
        gsl_matrix_free(metric);

      }

      break;

    default:
      XLALPrintError("\nInvalid grid type '%d'\n\n", init->gridType );
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

/** Return the spin-range spanned by the given template bank stored in the
 * *opaque* DopplerFullScanState.
 *
 * \note The user cannot directly access any internal fields of that opaque type,
 * which is why we need this API.
 */
int
XLALGetDopplerSpinRange ( PulsarSpinRange *spinRange, const DopplerFullScanState *scan )
{
  if ( spinRange == NULL )
    XLAL_ERROR ( XLAL_EINVAL, "\nPulsarSpinRange pointer 'spinRange' is NULL\n" );

  if ( scan == NULL )
    XLAL_ERROR ( XLAL_EINVAL, "\nDopplerFullScanState pointer 'scan' is NULL\n" );

  (*spinRange) = scan->spinRange;	// simple struct-copy is all that's needed

  return XLAL_SUCCESS;

} /* XLALGetDopplerSpinRange() */

/** Initialize Doppler-scanner to emulate an old-style factored template grid: 'sky x f0dot x f1dot x f2dot x f3dot'.
 *  This is a compatiblity-mode with the previous implementation currently also used in ComputeFStatistic.c.
 */
void
initFactoredGrid (LALStatus *status,				/**< pointer to LALStatus structure */
		  DopplerFullScanState *scan,			/**< [bout] initialized Doppler scan state */
		  const DopplerFullScanInit *init		/**< [in] initialization parameters */
		  )
{
  DopplerSkyScanInit skyScanInit = empty_DopplerSkyScanInit;
  SkyPosition skypos;
  factoredGridScan_t *fscan = NULL;
  UINT4 i;

  INITSTATUS(status);
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
  XLALNormalizeSkyPosition ( &skypos.longitude, &skypos.latitude );
  fscan->thisPoint.Alpha = skypos.longitude;
  fscan->thisPoint.Delta = skypos.latitude;
  /* set spins to start */
  for (i=0; i < PULSAR_MAX_SPINS; i ++ )
    fscan->thisPoint.fkdot[i] = scan->spinRange.fkdot[i];

  { /* count total number of templates */
    REAL8 nSky, nTot;
    REAL8 nSpins[PULSAR_MAX_SPINS];
    nSky = fscan->skyScan.numSkyGridPoints;
    for ( i=0; i < PULSAR_MAX_SPINS; i ++ )
      nSpins[i] = floor( scan->spinRange.fkdotBand[i] / fscan->skyScan.dfkdot[i] ) + 1.0;
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

	case GRID_SPINDOWN_SQUARE: /* square parameter space */
	case GRID_SPINDOWN_AGEBRK: /* age-braking index parameter space */
	  LogPrintf(LOG_DEBUG, "Counting spindown lattice templates ... ");
	  scan->numTemplates = (REAL8)XLALCountTotalFlatLatticePoints(scan->spindownTiling);
	  LogPrintfVerbatim(LOG_DEBUG, "%0.0f\n", scan->numTemplates);
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

  /* This traps coding errors in the calling routine. */
  if ( pos == NULL || scan == NULL ) {
    XLAL_ERROR ( XLAL_EINVAL );
  }
  if ( scan->state == STATE_IDLE ) {
    XLALPrintError ("\nCalled XLALNextDopplerPos() on un-initialized DopplerFullScanState !\n\n");
    XLAL_ERROR ( XLAL_EINVAL );
  }

  /* is this search finished? then return '1' */
  if (  scan->state == STATE_FINISHED )
    return 1;

  // set refTime in returned template to the refTime of the grid
  pos->refTime  = scan->spinRange.refTime;

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
      INIT_MEM(pos->fkdot);
      pos->fkdot[0] = scan->thisGridPoint->entry.data[0];
      pos->Alpha    = scan->thisGridPoint->entry.data[1];
      pos->Delta    = scan->thisGridPoint->entry.data[2];
      pos->fkdot[1] = scan->thisGridPoint->entry.data[3];
      pos->fkdot[2] = scan->thisGridPoint->entry.data[4];
      pos->fkdot[3] = scan->thisGridPoint->entry.data[5];
      pos->orbit    = NULL;
      /* advance to next grid point */
      if ( ( scan->thisGridPoint = scan->thisGridPoint->next ) == NULL )
	scan->state = STATE_FINISHED;

      break;

    case GRID_METRIC_LATTICE:
      if ( XLALgetCurrentDopplerPos ( pos, scan->latticeScan, COORDINATESYSTEM_EQUATORIAL ) ) {
	XLAL_ERROR ( XLAL_EFUNC );
      }
      /* advance to next point */
      ret = XLALadvanceLatticeIndex ( scan->latticeScan );
      if ( ret < 0 ) {
	XLAL_ERROR ( XLAL_EFUNC );
      }
      else if ( ret == 1 )
	{
	  XLALPrintError ( "\n\nXLALadvanceLatticeIndex(): this was the last lattice points!\n\n");
	  scan->state = STATE_FINISHED;
	}
#if 0
      { /* debugging */
	gsl_vector_int *lal_index = NULL;
	XLALgetCurrentLatticeIndex ( &lal)index, scan->latticeScan );
	XLALfprintfGSLvector_int ( stderr, "%d", lal_index );
	gsl_vector_int_free ( lal_index );
      }
#endif

      break;

    case GRID_SPINDOWN_SQUARE: /* square parameter space */
    case GRID_SPINDOWN_AGEBRK: /* age-braking index parameter space */
      {

	/* Advance to next tile */
        int retn = XLALNextFlatLatticePoint(scan->spindownTiling);
        if (xlalErrno != 0) {
	  XLALPrintError("\nGRID_SPINDOWN_{SQUARE,AGEBRK}: XLALNextFlatLatticeTile failed\n");
	  return -1;
        }

        if (retn >= 0) {

	  /* Found a point */
          const gsl_vector* current = XLALGetFlatLatticePoint(scan->spindownTiling);
          pos->Alpha      = gsl_vector_get(current, 0);
          pos->Delta      = gsl_vector_get(current, 1);
          pos->fkdot[0]   = gsl_vector_get(current, 2);
          for (size_t i = 1; i < PULSAR_MAX_SPINS; ++i)
            pos->fkdot[i] = gsl_vector_get(current, i + 2);

          return 0;

        } else {

	  /* No more points */
	  scan->state = STATE_FINISHED;
	  return 1;

	}

      }

      break;

    default:
      XLALPrintError("\nInvalid grid type '%d'\n\n", scan->gridType );
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

  range = &(scan->spinRange);	/* shortcut */

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
		      XLALNormalizeSkyPosition ( &skypos.longitude, &skypos.latitude );
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

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* This traps coding errors in the calling routine. */
  ASSERT( scan, status, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL );
  ASSERT( *scan, status, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL );

  if ( (*scan)->factoredScan ) {
    TRY ( FreeDopplerSkyScan ( status->statusPtr, &((*scan)->factoredScan->skyScan) ), status );
    XLALFree ( (*scan)->factoredScan );
  }

  if ( (*scan)->covering )
    XLALREAL8VectorListDestroy ( (*scan)->covering );

  if ( (*scan)->latticeScan ) {
    XLALFreeDopplerLatticeScan ( &((*scan)->latticeScan) );
  }

  if ((*scan)->spindownTiling) {
    XLALDestroyFlatLatticeTiling((*scan)->spindownTiling);
  }

  if ( (*scan)->skyRegion.vertices)
    XLALFree ( (*scan)->skyRegion.vertices);


  XLALFree ( (*scan) );
  (*scan) = NULL;

  DETATCHSTATUSPTR (status);
  RETURN( status );

} /* FreeDopplerSkyScan() */



/** load a full multi-dim template grid from the file init->gridFile,
 * the file-format is: lines of 6 columns, which are:
 *
 * Freq   Alpha  Delta  f1dot  f2dot  f3dot
 *
 * \note
 * *) this function returns the effective spinRange covered by the read-in template bank
 *    by storing it in scan->spinRange, potentially overwriting any previous user-input values in there.
 *
 * *) a possible future extension should probably *clip* the template-bank to the user-specified ranges,
 *    then return the effective ranges spanned by the resultant template bank.
 *
 * *) in order to avoid surprises until such a feature is implemented, we currently return an error if
 *    any of the input spinRanges are non-zero
 *
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
  UINT4 numTemplates;
  FILE *fp;

  INITSTATUS(status);
  ATTATCHSTATUSPTR ( status );

  ASSERT ( scan, status, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL);
  ASSERT ( init, status, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL);
  ASSERT ( init->gridFile, status, DOPPLERSCANH_ENULL, DOPPLERSCANH_MSGENULL);
  ASSERT ( scan->state == STATE_IDLE, status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT);

  /* Check that all user-input spin- and sky-ranges are zero, otherwise fail!
   *
   * NOTE: In the future we should allow combining the user-input ranges with
   * those found in the grid-file by forming the intersection, ie *clipping*
   * of the read-in grids to the user-input ranges.
   * Right now we require empty ranges input, and report back the ranges from the grid-file
   *
   */
  if ( init->searchRegion.skyRegionString != NULL ) {
    XLALPrintError ("\n%s: non-NULL skyRegion input currently not supported! skyRegion = '%s'\n\n", __func__, init->searchRegion.skyRegionString );
    ABORT ( status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );
  }
  for ( UINT4 s = 0; s < PULSAR_MAX_SPINS; s ++ )
    {
      if ( (init->searchRegion.fkdot[s] != 0) || (init->searchRegion.fkdotBand[s] != 0 )) {
        XLALPrintError ("\n%s: non-zero input spinRanges currently not supported! fkdot[%d] = %g, fkdotBand[%d] = %g\n\n",
                        __func__, s, init->searchRegion.fkdot[s], s, init->searchRegion.fkdotBand[s] );
        ABORT ( status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );
      }
    } /* for s < max_spins */

  /* open input data file */
  if ( (fp = LALOpenDataFile (init->gridFile)) == NULL) {
    XLALPrintError ("Could not open data-file: `%s`\n\n", init->gridFile);
    ABORT (status, CONFIGFILEH_EFILE, CONFIGFILEH_MSGEFILE);
  }

  /* prepare grid-entry buffer */
  if ( (entry = XLALCreateREAL8Vector ( 6 ) ) == NULL ) {
    ABORT (status, DOPPLERSCANH_EMEM, DOPPLERSCANH_MSGEMEM);
  }

  /* keep track of the sky- and spinRanges spanned by the template bank */
  REAL8 FreqMax  = - LAL_REAL4_MAX, FreqMin  = LAL_REAL4_MAX;	// only using REAL4 ranges to avoid over/under flows, and should be enough
  REAL8 f1dotMax = - LAL_REAL4_MAX, f1dotMin = LAL_REAL4_MAX;
  REAL8 f2dotMax = - LAL_REAL4_MAX, f2dotMin = LAL_REAL4_MAX;
  REAL8 f3dotMax = - LAL_REAL4_MAX, f3dotMin = LAL_REAL4_MAX;

  REAL8 alphaMax = - LAL_REAL4_MAX, alphaMin = LAL_REAL4_MAX;
  REAL8 deltaMax = - LAL_REAL4_MAX, deltaMin = LAL_REAL4_MAX;

  /* parse this list of lines into a full grid */
  numTemplates = 0;
  tail = &head;	/* head will remain empty! */
  while ( ! feof ( fp ) )
    {
      REAL8 Freq, Alpha, Delta, f1dot, f2dot, f3dot;
      // File format expects lines containing 6 columns: Freq   Alpha  Delta  f1dot  f2dot  f3dot
      if ( 6 != fscanf( fp, "%" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT "\n",
                        &Freq, &Alpha, &Delta, &f1dot, &f2dot, &f3dot ) )
        {
          LogPrintf (LOG_CRITICAL,"ERROR: Failed to parse 6 REAL8's from line %d in grid-file '%s'\n\n", numTemplates + 1, init->gridFile);
          if ( head.next )
            XLALREAL8VectorListDestroy (head.next);
          ABORT (status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT);
        }

      /* keep track of maximal spans */
      alphaMin = fmin ( Alpha, alphaMin );
      deltaMin = fmin ( Delta, deltaMin );
      FreqMin  = fmin ( Freq,  FreqMin );
      f1dotMin = fmin ( f1dot, f1dotMin );
      f2dotMin = fmin ( f2dot, f2dotMin );
      f3dotMin = fmin ( f3dot, f3dotMin );

      alphaMax = fmax ( Alpha, alphaMax );
      deltaMax = fmax ( Delta, deltaMax );
      FreqMax  = fmax ( Freq,  FreqMax );
      f1dotMax = fmax ( f1dot, f1dotMax );
      f2dotMax = fmax ( f2dot, f2dotMax );
      f3dotMax = fmax ( f3dot, f3dotMax );

      /* add this entry to template-bank list */
      entry->data[0] = Freq;
      entry->data[1] = Alpha;
      entry->data[2] = Delta;
      entry->data[3] = f1dot;
      entry->data[4] = f2dot;
      entry->data[5] = f3dot;

      if ( (tail = XLALREAL8VectorListAddEntry (tail, entry)) == NULL )
	{
	  if ( head.next )
	    XLALREAL8VectorListDestroy (head.next);
	  ABORT ( status, DOPPLERSCANH_EXLAL, DOPPLERSCANH_MSGEXLAL );
	}

      numTemplates ++ ;

    } /* while !feof(fp) */

  XLALDestroyREAL8Vector ( entry );

  /* ---------- update scan-state  ---------- */

  // ----- report back ranges actually spanned by grid-file

  CHAR *skyRegionString = NULL;
  REAL8 eps = LAL_REAL8_EPS;
  TRY ( SkySquare2String ( status->statusPtr, &skyRegionString, alphaMin, deltaMin, (alphaMax - alphaMin) + eps, (deltaMax - deltaMin) + eps ), status );
  // note: we slight expanded the enclosing sky-square by eps to avoid complaints when a grid-file contains
  // only points in a line, which is perfectly valid here.
  TRY ( ParseSkyRegionString ( status->statusPtr, &scan->skyRegion, skyRegionString ), status );
  XLALFree ( skyRegionString );

  scan->spinRange.fkdot[0]     = FreqMin;
  scan->spinRange.fkdotBand[0] = FreqMax - FreqMin;

  scan->spinRange.fkdot[1]     = f1dotMin;
  scan->spinRange.fkdotBand[1] = f1dotMax - f1dotMin;

  scan->spinRange.fkdot[2]     = f2dotMin;
  scan->spinRange.fkdotBand[2] = f2dotMax - f2dotMin;

  scan->spinRange.fkdot[3]     = f3dotMin;
  scan->spinRange.fkdotBand[3] = f3dotMax - f3dotMin;

  scan->numTemplates = numTemplates;
  scan->covering = head.next;	/* pass result (without head!) */
  scan->thisGridPoint = scan->covering;	/* init to start */

  LogPrintf (LOG_DEBUG, "Template grid: nTot = %.0f\n", 1.0 * numTemplates );
  LogPrintf (LOG_DEBUG, "Spanned ranges: Freq in [%g, %g], f1dot in [%g, %g], f2dot in [%g, %g], f3dot in [%g, %g]\n",
             FreqMin, FreqMax, f1dotMin, f1dotMax, f2dotMin, f2dotMax, f3dotMin, f3dotMax );

  DETATCHSTATUSPTR ( status );
  RETURN ( status );

} /* loadFullGridFile() */
