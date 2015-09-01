/*
 *  Copyright (C) 2007, 2008, 2012 Karl Wette
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
 * as opposed to the "factored" grids implemented in DopplerScan.[ch]
 */

/*---------- INCLUDES ----------*/
#include <math.h>

#include <lal/LALStdlib.h>
#include <lal/FileIO.h>
#include <lal/DetectorSite.h>
#include <lal/LALError.h>
#include <lal/ConfigFile.h>
#include <lal/LogPrintf.h>

#include <lal/FlatPulsarMetric.h>
#include <lal/DopplerFullScan.h>

/*---------- DEFINES ----------*/
#define MIN(x,y) (x < y ? x : y)
#define MAX(x,y) (x > y ? x : y)

#define TRUE (1==1)
#define FALSE (1==0)

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/*---------- internal types ----------*/
typedef struct {
  PulsarDopplerParams thisPoint; /**< current doppler-position of the scan */
  DopplerSkyScanState skyScan;  /**< keep track of sky-grid stepping */
} factoredGridScan_t;

/** doubly linked list of REAL8-vectors (physical vectors) */
typedef struct tagREAL8VectorList
{
  REAL8Vector entry;
  struct tagREAL8VectorList *next;
  struct tagREAL8VectorList *prev;
} REAL8VectorList;

/** ----- internal [opaque] type to store the state of a FULL multidimensional grid-scan ----- */
struct tagDopplerFullScanState {
  INT2 state;                   /**< idle, ready or finished */
  DopplerGridType gridType;     /**< what type of grid are we dealing with */
  REAL8 numTemplates;           /**< total number of templates in the grid */
  PulsarSpinRange spinRange;    /**< spin-range covered by template bank */
  SkyRegion skyRegion;          /**< sky-range covered by template bank */

  /* ----- full multi-dim parameter-space grid stuff ----- */
  REAL8VectorList *covering;            /**< multi-dimensional covering */
  REAL8VectorList *thisGridPoint;       /**< pointer to current grid-point */

  /* ----- spindown lattice tiling ----- */
  LatticeTiling *spindownTiling;              /**< spindown lattice tiling */
  LatticeTilingIterator *spindownTilingItr;   /**< iterator over spindown lattice tiling */
  gsl_vector *spindownTilingPoint;            /**< current point in spindown lattice tiling */

  /* ----- emulate old-style factored grids */
  factoredGridScan_t *factoredScan;     /**< only used to emulate FACTORED grids sky x Freq x f1dot */

} /* struct DopplerFullScanState */;

/*---------- internal prototypes ----------*/
int XLALInitFactoredGrid ( DopplerFullScanState *scan,  const DopplerFullScanInit *init );
int nextPointInFactoredGrid (PulsarDopplerParams *pos, DopplerFullScanState *scan);
int XLALLoadFullGridFile ( DopplerFullScanState *scan, const DopplerFullScanInit *init );

/*==================== FUNCTION DEFINITIONS ====================*/

/**
 * Set up a full multi-dimensional grid-scan.
 * Currently this only emulates a 'factored' grid-scan with 'sky x Freq x f1dot ...' , but
 * keeps all details within the DopplerScan module for future extension to real multidimensional
 * grids.
 *
 * NOTE: Use 'XLALNextDopplerPos()' to step through this template grid.
 *
 */
DopplerFullScanState *
XLALInitDopplerFullScan ( const DopplerFullScanInit *init       /**< [in] initialization parameters */
                          )
{
  XLAL_CHECK_NULL ( init != NULL, XLAL_EINVAL );

  DopplerFullScanState *thisScan;
  XLAL_CHECK_NULL ( (thisScan = LALCalloc (1, sizeof(*thisScan) )) != NULL, XLAL_ENOMEM );

  thisScan->gridType = init->gridType;

  /* store the user-input spinRange (includes refTime) in DopplerFullScanState */
  thisScan->spinRange.refTime = init->searchRegion.refTime;
  memcpy ( thisScan->spinRange.fkdot, init->searchRegion.fkdot, sizeof(PulsarSpins) );
  memcpy ( thisScan->spinRange.fkdotBand, init->searchRegion.fkdotBand, sizeof(PulsarSpins) );

  // check that some old metric-codes aren't used with refTime!=startTime, which they don't handle correctly
  switch ( thisScan->gridType )
    {
    case GRID_METRIC:
    case GRID_METRIC_SKYFILE:
    case GRID_SPINDOWN_SQUARE: /* square parameter space */
    case GRID_SPINDOWN_AGEBRK: /* age-braking index parameter space */

      XLAL_CHECK_NULL ( XLALGPSDiff ( &init->startTime, &init->searchRegion.refTime ) == 0, XLAL_EINVAL,
                        "NOTE: gridType={metric,4,spin-square,spin-age-brk} only work for refTime (%f) == startTime (%f)!\n",
                        XLALGPSGetREAL8(&(init->searchRegion.refTime)), XLALGPSGetREAL8(&(init->startTime)) );;

      break;
    default:
      break;
    }

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
      XLAL_CHECK_NULL ( XLALInitFactoredGrid ( thisScan, init ) == XLAL_SUCCESS, XLAL_EFUNC );
      break;

      /* ----- multi-dimensional covering of full parameter space ----- */
    case GRID_FILE_FULLGRID:
      XLAL_CHECK_NULL ( XLALLoadFullGridFile ( thisScan, init ) == XLAL_SUCCESS, XLAL_EFUNC );
      break;

    case GRID_SPINDOWN_SQUARE: /* square parameter space */
    case GRID_SPINDOWN_AGEBRK: /* age-braking index parameter space */
      {
        const size_t n = 2 + PULSAR_MAX_SPINS;

        /* Check that the reference time is the same as the start time */
        XLAL_CHECK_NULL ( XLALGPSCmp ( &thisScan->spinRange.refTime, &init->startTime) == 0, XLAL_EINVAL,
                          "\nGRID_SPINDOWN_{SQUARE,AGEBRK}: This option currently restricts the reference time to be the same as the start time.\n");

        /* Create a vector to hold lattice tiling parameter-space points */
        XLAL_CHECK_NULL ( (thisScan->spindownTilingPoint = gsl_vector_alloc(n)) != NULL, XLAL_ENOMEM,
                          "\nGRID_SPINDOWN_{SQUARE,AGEBRK}: gsl_vector_alloc failed\n");

        /* Create a lattice tiling */
        XLAL_CHECK_NULL ( (thisScan->spindownTiling = XLALCreateLatticeTiling(n)) != NULL, XLAL_EFUNC );

        /* Parse the sky region string and check that it consists of only one point, and set bounds on it */
        SkyRegion XLAL_INIT_DECL(sky);
        XLAL_CHECK_NULL ( XLALParseSkyRegionString ( &sky, init->searchRegion.skyRegionString ) == XLAL_SUCCESS, XLAL_EFUNC );
        XLAL_CHECK_NULL ( sky.numVertices == 1, XLAL_EINVAL, "\nGRID_SPINDOWN_{SQUARE,AGEBRK}: This option can only handle a single sky position.\n");
        XLAL_CHECK_NULL ( sky.vertices[0].system == COORDINATESYSTEM_EQUATORIAL, XLAL_EINVAL, "\nGRID_SPINDOWN_{SQUARE,AGEBRK}: This option only understands COORDINATESYSTEM_EQUATORIAL\n");

        XLAL_CHECK_NULL ( XLALSetLatticeTilingConstantBound(thisScan->spindownTiling, 0, sky.vertices[0].longitude, sky.vertices[0].longitude) == XLAL_SUCCESS, XLAL_EFUNC );

        XLAL_CHECK_NULL ( XLALSetLatticeTilingConstantBound(thisScan->spindownTiling, 1, sky.vertices[0].latitude, sky.vertices[0].latitude) == XLAL_SUCCESS, XLAL_EFUNC );
        if (sky.vertices) {
          XLALFree (sky.vertices);
        }

        /* Set up parameter space */
        if (thisScan->gridType == GRID_SPINDOWN_SQUARE) { /* square parameter space */

          /* Set square bounds on the frequency and spindowns */
          for (size_t i = 0; i < PULSAR_MAX_SPINS; ++i) {
            XLAL_CHECK_NULL ( XLALSetLatticeTilingConstantBound(thisScan->spindownTiling, 2 + i, init->searchRegion.fkdot[i], init->searchRegion.fkdot[i] + init->searchRegion.fkdotBand[i]) == XLAL_SUCCESS, XLAL_EFUNC );
          }

        } else if (thisScan->gridType == GRID_SPINDOWN_AGEBRK) { /* age-braking index parameter space */

          /* Get age and braking index from extra arguments */
          const REAL8 spindownAge = init->extraArgs[0];
          const REAL8 minBraking = init->extraArgs[1];
          const REAL8 maxBraking = init->extraArgs[2];

          /* Set age-braking index parameter space */
          XLAL_CHECK_NULL ( XLAL_SUCCESS == XLALSetLatticeTilingConstantBound(thisScan->spindownTiling, 2, init->searchRegion.fkdot[0], init->searchRegion.fkdot[0] + init->searchRegion.fkdotBand[0]), XLAL_EFUNC );
          XLAL_CHECK_NULL ( XLAL_SUCCESS == XLALSetLatticeTilingF1DotAgeBrakingBound(thisScan->spindownTiling, 2, 3, spindownAge, minBraking, maxBraking), XLAL_EFUNC );
          XLAL_CHECK_NULL ( XLAL_SUCCESS == XLALSetLatticeTilingF2DotBrakingBound(thisScan->spindownTiling, 2, 3, 4, minBraking, maxBraking), XLAL_EFUNC );

          /* This current only goes up to second spindown, so bound higher dimensions */
          for (size_t i = 3; i < PULSAR_MAX_SPINS; ++i) {
            XLAL_CHECK_NULL ( XLAL_SUCCESS == XLALSetLatticeTilingConstantBound(thisScan->spindownTiling, 2 + i, init->searchRegion.fkdot[i], init->searchRegion.fkdot[i] + init->searchRegion.fkdotBand[i]), XLAL_EFUNC );
          }

        }

        /* Create a lattice tiling with Anstar lattice and spindown metric */
        gsl_matrix* metric;
        XLAL_CHECK_NULL ( (metric = gsl_matrix_alloc(n, n)) != NULL, XLAL_ENOMEM );
        gsl_matrix_set_identity(metric);
        gsl_matrix_view spin_metric = gsl_matrix_submatrix(metric, 2, 2, PULSAR_MAX_SPINS, PULSAR_MAX_SPINS);
        XLAL_CHECK_NULL ( XLALSpindownMetric(&spin_metric.matrix, init->Tspan) == XLAL_SUCCESS, XLAL_EFUNC );
        XLAL_CHECK_NULL ( XLALSetTilingLatticeAndMetric(thisScan->spindownTiling, TILING_LATTICE_ANSTAR, metric, init->metricMismatch) == XLAL_SUCCESS, XLAL_EFUNC );

        /* Create iterator over flat lattice tiling */
        XLAL_CHECK_NULL ( (thisScan->spindownTilingItr = XLALCreateLatticeTilingIterator(thisScan->spindownTiling, n)) != NULL, XLAL_EFUNC );

        /* Cleanup */
        gsl_matrix_free(metric);

      }

      break;

    default:
      XLAL_ERROR_NULL ( XLAL_EINVAL, "\nInvalid grid type '%d'\n\n", init->gridType );
      break;
    } /* switch gridType */

  /* we're ready */
  thisScan->state = STATE_READY;

  /* return result */
  return thisScan;

} // XLALInitDopplerFullScan()


/** \deprecated Use XLALInitDopplerFullScan() instead.
 */
void
InitDopplerFullScan(LALStatus *status,                  /**< pointer to LALStatus structure */
                    DopplerFullScanState **scan,        /**< [out] initialized Doppler scan state */
                    const DopplerFullScanInit *init     /**< [in] initialization parameters */
                    )
{
  INITSTATUS(status);

  XLAL_PRINT_DEPRECATION_WARNING ( "XLALInitDopplerFullScan" );

  ASSERT ( scan, status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );
  ASSERT ( *scan == NULL, status, DOPPLERSCANH_ENONULL, DOPPLERSCANH_MSGENONULL );
  ASSERT ( init, status, DOPPLERSCANH_EINPUT, DOPPLERSCANH_MSGEINPUT );

  if ( ((*scan) = XLALInitDopplerFullScan ( init )) == NULL ) {
    ABORT (status, DOPPLERSCANH_EXLAL, DOPPLERSCANH_MSGEXLAL);
  }

  RETURN( status );

} /* InitDopplerFullScan() */

/**
 * Return the spin-range spanned by the given template bank stored in the
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

  (*spinRange) = scan->spinRange;       // simple struct-copy is all that's needed

  return XLAL_SUCCESS;

} /* XLALGetDopplerSpinRange() */

/**
 * Initialize Doppler-scanner to emulate an old-style factored template grid: 'sky x f0dot x f1dot x f2dot x f3dot'.
 * This is a compatiblity-mode with the previous implementation currently also used in ComputeFstatistic.c.
 */
int
XLALInitFactoredGrid ( DopplerFullScanState *scan,                      /**< [bout] initialized Doppler scan state */
                       const DopplerFullScanInit *init          /**< [in] initialization parameters */
                       )
{
  XLAL_CHECK ( scan, XLAL_EINVAL );
  XLAL_CHECK ( scan->state == STATE_IDLE, XLAL_EINVAL );
  XLAL_CHECK ( init, XLAL_EINVAL );

  DopplerSkyScanInit XLAL_INIT_DECL(skyScanInit);
  SkyPosition skypos;
  factoredGridScan_t *fscan = NULL;

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
  skyScanInit.ephemeris = init->ephemeris;              /* used only by Ephemeris-based metric */
  skyScanInit.skyGridFile = init->gridFile;
  skyScanInit.skyRegionString = init->searchRegion.skyRegionString;
  skyScanInit.Freq = init->searchRegion.fkdot[0] + init->searchRegion.fkdotBand[0];

  XLAL_CHECK ( (fscan = LALCalloc ( 1, sizeof( *fscan ))) != NULL, XLAL_ENOMEM );
  scan->factoredScan = fscan;
  XLAL_CHECK ( XLALInitDopplerSkyScan ( &(fscan->skyScan), &skyScanInit) == XLAL_SUCCESS, XLAL_EFUNC );

  /* overload spin step-sizes with user-settings if given */
  for (UINT4 i=0; i < PULSAR_MAX_SPINS; i ++ ) {
    if ( init->stepSizes.fkdot[i] ) {
      fscan->skyScan.dfkdot[i] = init->stepSizes.fkdot[i];
    }
  }

  /* ----- set Doppler-scanner to start-point ----- */
  fscan->thisPoint.refTime = init->searchRegion.refTime;        /* set proper reference time for spins */

  /* normalize skyposition: correctly map into [0,2pi]x[-pi/2,pi/2] */
  skypos.longitude = fscan->skyScan.skyNode->Alpha;
  skypos.latitude  = fscan->skyScan.skyNode->Delta;
  skypos.system = COORDINATESYSTEM_EQUATORIAL;
  XLALNormalizeSkyPosition ( &skypos.longitude, &skypos.latitude );
  fscan->thisPoint.Alpha = skypos.longitude;
  fscan->thisPoint.Delta = skypos.latitude;
  /* set spins to start */
  for (UINT4 i=0; i < PULSAR_MAX_SPINS; i ++ )
    fscan->thisPoint.fkdot[i] = scan->spinRange.fkdot[i];

  { /* count total number of templates */
    REAL8 nSky, nTot;
    REAL8 nSpins[PULSAR_MAX_SPINS];
    nSky = fscan->skyScan.numSkyGridPoints;
    for ( UINT4 i=0; i < PULSAR_MAX_SPINS; i ++ ) {
      nSpins[i] = floor( scan->spinRange.fkdotBand[i] / fscan->skyScan.dfkdot[i] ) + 1.0;
    }
    nTot = nSky;
    for ( UINT4 i=0; i < PULSAR_MAX_SPINS; i ++ ) {
      nTot *= nSpins[i];
    }
    scan->numTemplates = nTot;
    XLALPrintInfo ("Template grid: nSky x nFreq x nf1dot = %.0f x %.0f x %.0f = %.0f \n", nSky, nSpins[0], nSpins[1], nTot );
  }

  return XLAL_SUCCESS;

} // XLALInitFactoredGrid()


/**
 * Return (and compute if not done so yet) the total number of Doppler templates
 * of the DopplerScan \a scan
 */
REAL8
XLALNumDopplerTemplates ( DopplerFullScanState *scan)
{
  if ( ! scan->numTemplates )   /* not pre-computed already ? */
    {
      switch ( scan->gridType )
        {
        /* case GRID_METRIC_LATTICE: */
        /*   LogPrintf ( LOG_DEBUG, "Now counting number of templates in lattice ... "); */
        /*   scan->numTemplates = XLALCountLatticeTemplates ( scan->latticeScan ); */
        /*   LogPrintfVerbatim( LOG_DEBUG, " done. (%.0f)\n", scan->numTemplates ); */
        /*   break; */

        case GRID_SPINDOWN_SQUARE: /* square parameter space */
        case GRID_SPINDOWN_AGEBRK: /* age-braking index parameter space */
          LogPrintf(LOG_DEBUG, "Counting spindown lattice templates ... ");
          scan->numTemplates = (REAL8)XLALTotalLatticeTilingPoints(scan->spindownTilingItr);
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

/**
 * Function to step through the full template grid point by point.
 * Normal return = 0,
 * errors return -1,
 * end of scan is signalled by return = 1
 */
int
XLALNextDopplerPos(PulsarDopplerParams *pos, DopplerFullScanState *scan)
{

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
      XLAL_INIT_MEM(pos->fkdot);
      pos->fkdot[0] = scan->thisGridPoint->entry.data[0];
      pos->Alpha    = scan->thisGridPoint->entry.data[1];
      pos->Delta    = scan->thisGridPoint->entry.data[2];
      pos->fkdot[1] = scan->thisGridPoint->entry.data[3];
      pos->fkdot[2] = scan->thisGridPoint->entry.data[4];
      pos->fkdot[3] = scan->thisGridPoint->entry.data[5];
      pos->asini = 0;   // isolated pulsar
      /* advance to next grid point */
      if ( ( scan->thisGridPoint = scan->thisGridPoint->next ) == NULL )
        scan->state = STATE_FINISHED;

      break;

/*     case GRID_METRIC_LATTICE: */
/*       if ( XLALgetCurrentDopplerPos ( pos, scan->latticeScan, COORDINATESYSTEM_EQUATORIAL ) ) { */
/*      XLAL_ERROR ( XLAL_EFUNC ); */
/*       } */
/*       /\* advance to next point *\/ */
/*       ret = XLALadvanceLatticeIndex ( scan->latticeScan ); */
/*       if ( ret < 0 ) { */
/*      XLAL_ERROR ( XLAL_EFUNC ); */
/*       } */
/*       else if ( ret == 1 ) */
/*      { */
/*        XLALPrintError ( "\n\nXLALadvanceLatticeIndex(): this was the last lattice points!\n\n"); */
/*        scan->state = STATE_FINISHED; */
/*      } */
/* #if 0 */
/*       { /\* debugging *\/ */
/*      gsl_vector_int *lal_index = NULL; */
/*      XLALgetCurrentLatticeIndex ( &lal)index, scan->latticeScan ); */
/*      XLALfprintfGSLvector_int ( stderr, "%d", lal_index ); */
/*      gsl_vector_int_free ( lal_index ); */
/*       } */
/* #endif */

      break;

    case GRID_SPINDOWN_SQUARE: /* square parameter space */
    case GRID_SPINDOWN_AGEBRK: /* age-braking index parameter space */
      {

        /* Advance to next tile */
        int retn = XLALNextLatticeTilingPoint(scan->spindownTilingItr, scan->spindownTilingPoint);
        if (retn < 0) {
          XLALPrintError("\nGRID_SPINDOWN_{SQUARE,AGEBRK}: XLALNextLatticeTilingPoint() failed\n");
          return -1;
        }

        if (retn > 0) {

          /* Found a point */
          pos->Alpha      = gsl_vector_get(scan->spindownTilingPoint, 0);
          pos->Delta      = gsl_vector_get(scan->spindownTilingPoint, 1);
          pos->fkdot[0]   = gsl_vector_get(scan->spindownTilingPoint, 2);
          for (size_t i = 1; i < PULSAR_MAX_SPINS; ++i) {
            pos->fkdot[i] = gsl_vector_get(scan->spindownTilingPoint, i + 2);
          }

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

/**
 * return current grid-point and step forward one template in 'factored' grids (sky x f0dot x f1dot ... )
 * return 0 = OK, -1 = ERROR
 */
int
nextPointInFactoredGrid (PulsarDopplerParams *pos, DopplerFullScanState *scan)
{
  PulsarDopplerParams XLAL_INIT_DECL(nextPos);
  factoredGridScan_t *fscan;
  PulsarSpinRange *range;
  PulsarSpins fkdotMax;
  SkyPosition skypos;

  if ( pos == NULL || scan == NULL )
    return -1;

  if ( ( fscan = scan->factoredScan ) == NULL )
    return -1;

  range = &(scan->spinRange);   /* shortcut */

  (*pos) = fscan->thisPoint;    /* RETURN current Doppler-point (struct-copy) */

  nextPos = fscan->thisPoint;   /* struct-copy: start from current point to get next one */

  /* shortcuts: spin boundaries */
  fkdotMax[3] = range->fkdot[3] + range->fkdotBand[3];
  fkdotMax[2] = range->fkdot[2] + range->fkdotBand[2];
  fkdotMax[1] = range->fkdot[1] + range->fkdotBand[1];
  fkdotMax[0] = range->fkdot[0] + range->fkdotBand[0];

  /* try to advance to next template */
  nextPos.fkdot[0] += fscan->skyScan.dfkdot[0];                 /* f0dot one forward */
  if ( nextPos.fkdot[0] >  fkdotMax[0] )                /* too far? */
    {
      nextPos.fkdot[0] = range->fkdot[0];               /* f0dot return to start */
      nextPos.fkdot[1] += fscan->skyScan.dfkdot[1];     /* f1dot one step forward */
      if ( nextPos.fkdot[1] > fkdotMax[1] )
        {
          nextPos.fkdot[1] = range->fkdot[1];           /* f1dot return to start */
          nextPos.fkdot[2] += fscan->skyScan.dfkdot[2];         /* f2dot one forward */
          if ( nextPos.fkdot[2] > fkdotMax[2] )
            {
              nextPos.fkdot[2] = range->fkdot[2];       /* f2dot return to start */
              nextPos.fkdot[3] += fscan->skyScan.dfkdot[3]; /* f3dot one forward */
              if ( nextPos.fkdot[3] > fkdotMax[3] )
                {
                  nextPos.fkdot[3] = range->fkdot[3];   /* f3dot return to start */
                                                        /* skygrid one forward */
                  if ( (fscan->skyScan.skyNode = fscan->skyScan.skyNode->next) == NULL ) /* no more sky-points ?*/
                    {
                      fscan->skyScan.state = STATE_FINISHED;    /* avoid warning when freeing */
                      scan->state = STATE_FINISHED;     /* we're done */
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

static void
XLALREAL8VectorListDestroy (REAL8VectorList *head)
{
  REAL8VectorList *ptr, *next;

  if ( !head )
    return;

  next = head;

  do
    {
      /* step to next element */
      ptr = next;
      /* remember pointer to next element */
      next = ptr->next;
      /* free current element */
      LALFree (ptr->entry.data);
      LALFree (ptr);

    } while ( (ptr = next) != NULL );

  return;
} /* XLALREAL8VectorListDestroy() */

/**
 * Destroy the a full DopplerFullScanState structure
 */
void
XLALDestroyDopplerFullScan ( DopplerFullScanState *scan )
{
  if ( scan == NULL ) {
    return;
  }

  if ( scan->factoredScan ) {
    XLALDestroyDopplerSkyScan ( &(scan->factoredScan->skyScan) );
    XLALFree ( scan->factoredScan );
  }

  if ( scan->covering ) {
    XLALREAL8VectorListDestroy ( scan->covering );
  }

  if (scan->spindownTiling) {
    XLALDestroyLatticeTiling(scan->spindownTiling);
  }
  if (scan->spindownTilingItr) {
    XLALDestroyLatticeTilingIterator(scan->spindownTilingItr);
  }
  if (scan->spindownTilingPoint) {
    gsl_vector_free(scan->spindownTilingPoint);
  }

  if ( scan->skyRegion.vertices) {
    XLALFree ( scan->skyRegion.vertices);
  }

  XLALFree ( scan );

  return;

} // XLALDestroyDopplerFullScan()


static REAL8VectorList *
XLALREAL8VectorListAddEntry (REAL8VectorList *head, const REAL8Vector *entry)
{
  UINT4 dim;
  REAL8VectorList *ptr = NULL;  /* running list-pointer */
  REAL8VectorList *newElement = NULL;   /* new list-element */
  /* check illegal input */
  if ( (head == NULL) || (entry == NULL) )
    return NULL;

  /* find tail of list */
  ptr = head;
  while ( ptr->next )
    ptr = ptr->next;

  /* construct new list-element */
  dim = entry->length;
  if ( (newElement = LALCalloc (1, sizeof (*newElement))) == NULL)
    return NULL;
  if ( (newElement->entry.data = LALCalloc (dim, sizeof(entry->data[0]))) == NULL ) {
    LALFree (newElement);
    return NULL;
  }
  newElement->entry.length = dim;
  memcpy (newElement->entry.data, entry->data, dim * sizeof(entry->data[0]) );

  /* link this to the tail of list */
  ptr->next = newElement;
  newElement->prev = ptr;

  return newElement;

} /* XLALREAL8VectorListAddEntry() */

/**
 * load a full multi-dim template grid from the file init->gridFile,
 * the file-format is: lines of 6 columns, which are:
 *
 * Freq   Alpha  Delta  f1dot  f2dot  f3dot
 *
 * \note
 * *) this function returns the effective spinRange covered by the read-in template bank
 * by storing it in scan->spinRange, potentially overwriting any previous user-input values in there.
 *
 * *) a possible future extension should probably *clip* the template-bank to the user-specified ranges,
 * then return the effective ranges spanned by the resultant template bank.
 *
 * *) in order to avoid surprises until such a feature is implemented, we currently return an error if
 * any of the input spinRanges are non-zero
 *
 */
int
XLALLoadFullGridFile ( DopplerFullScanState *scan,
                       const DopplerFullScanInit *init
                       )
{
  XLAL_CHECK ( (scan != NULL) && (init != NULL), XLAL_EINVAL );
  XLAL_CHECK ( init->gridFile != NULL, XLAL_EINVAL );
  XLAL_CHECK ( scan->state == STATE_IDLE, XLAL_EINVAL );

  REAL8VectorList XLAL_INIT_DECL(head);
  REAL8VectorList *tail = NULL;
  REAL8Vector *entry = NULL;
  UINT4 numTemplates;
  FILE *fp;


  /* Check that all user-input spin- and sky-ranges are zero, otherwise fail!
   *
   * NOTE: In the future we should allow combining the user-input ranges with
   * those found in the grid-file by forming the intersection, ie *clipping*
   * of the read-in grids to the user-input ranges.
   * Right now we require empty ranges input, and report back the ranges from the grid-file
   *
   */
  if ( init->searchRegion.skyRegionString != NULL ) {
    XLAL_ERROR ( XLAL_EINVAL, "\nnon-NULL skyRegion input currently not supported! skyRegion = '%s'\n\n", init->searchRegion.skyRegionString );
  }
  for ( UINT4 s = 0; s < PULSAR_MAX_SPINS; s ++ )
    {
      if ( (init->searchRegion.fkdot[s] != 0) || (init->searchRegion.fkdotBand[s] != 0 )) {
        XLAL_ERROR ( XLAL_EINVAL, "\nnon-zero input spinRanges currently not supported! fkdot[%d] = %g, fkdotBand[%d] = %g\n\n",
                     s, init->searchRegion.fkdot[s], s, init->searchRegion.fkdotBand[s] );
      }
    } /* for s < max_spins */

  /* open input data file */
  XLAL_CHECK ( (fp = LALFopen (init->gridFile, "r")) != NULL, XLAL_ESYS, "Could not open data-file: `%s`\n\n", init->gridFile );

  /* prepare grid-entry buffer */
  XLAL_CHECK ( (entry = XLALCreateREAL8Vector ( 6 ) ) != NULL, XLAL_EFUNC );

  /* keep track of the sky- and spinRanges spanned by the template bank */
  REAL8 FreqMax  = - LAL_REAL4_MAX, FreqMin  = LAL_REAL4_MAX;   // only using REAL4 ranges to avoid over/under flows, and should be enough
  REAL8 f1dotMax = - LAL_REAL4_MAX, f1dotMin = LAL_REAL4_MAX;
  REAL8 f2dotMax = - LAL_REAL4_MAX, f2dotMin = LAL_REAL4_MAX;
  REAL8 f3dotMax = - LAL_REAL4_MAX, f3dotMin = LAL_REAL4_MAX;

  REAL8 alphaMax = - LAL_REAL4_MAX, alphaMin = LAL_REAL4_MAX;
  REAL8 deltaMax = - LAL_REAL4_MAX, deltaMin = LAL_REAL4_MAX;

  /* parse this list of lines into a full grid */
  numTemplates = 0;
  tail = &head;         /* head will remain empty! */
  while ( ! feof ( fp ) )
    {
      REAL8 Freq, Alpha, Delta, f1dot, f2dot, f3dot;
      // File format expects lines containing 6 columns: Freq   Alpha  Delta  f1dot  f2dot  f3dot
      if ( 6 != fscanf( fp, "%" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT "\n",
                        &Freq, &Alpha, &Delta, &f1dot, &f2dot, &f3dot ) )
        {
          XLALPrintError ("ERROR: Failed to parse 6 REAL8's from line %d in grid-file '%s'\n\n", numTemplates + 1, init->gridFile);
          if ( head.next ) {
            XLALREAL8VectorListDestroy (head.next);
          }
          XLAL_ERROR ( XLAL_EINVAL );
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
          if ( head.next ) {
            XLALREAL8VectorListDestroy (head.next);
          }
          XLAL_ERROR ( XLAL_EINVAL );
        }

      numTemplates ++ ;

    } /* while !feof(fp) */

  XLALDestroyREAL8Vector ( entry );

  /* ---------- update scan-state  ---------- */

  // ----- report back ranges actually spanned by grid-file

  CHAR *skyRegionString = NULL;
  REAL8 eps = LAL_REAL8_EPS;
  XLAL_CHECK ( (skyRegionString = XLALSkySquare2String ( alphaMin, deltaMin, (alphaMax - alphaMin) + eps, (deltaMax - deltaMin) + eps )) != NULL, XLAL_EFUNC );

  // note: we slight expanded the enclosing sky-square by eps to avoid complaints when a grid-file contains
  // only points in a line, which is perfectly valid here.
  XLAL_CHECK ( XLALParseSkyRegionString ( &scan->skyRegion, skyRegionString ) == XLAL_SUCCESS, XLAL_EFUNC );
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
  scan->covering = head.next;   /* pass result (without head!) */
  scan->thisGridPoint = scan->covering;         /* init to start */

  XLALPrintInfo ( "Template grid: nTot = %.0f\n", 1.0 * numTemplates );
  XLALPrintInfo ( "Spanned ranges: Freq in [%g, %g], f1dot in [%g, %g], f2dot in [%g, %g], f3dot in [%g, %g]\n",
                  FreqMin, FreqMax, f1dotMin, f1dotMax, f2dotMin, f2dotMax, f3dotMin, f3dotMax );

  return XLAL_SUCCESS;

} /* XLALLoadFullGridFile() */

typedef struct {
  size_t freq_dim;
  double scale;
} F1DotAgeBrakingBoundInfo;

static double F1DotAgeBrakingBound(
  const void* data,
  const size_t dim UNUSED,
  const gsl_vector* point
  )
{

  // Get bounds info
  const F1DotAgeBrakingBoundInfo* info = (const F1DotAgeBrakingBoundInfo*)data;

  // Get current value of frequency
  const double freq = gsl_vector_get(point, info->freq_dim);

  // Return first spindown bounds
  return info->scale * freq;

}

int XLALSetLatticeTilingF1DotAgeBrakingBound(
  LatticeTiling* tiling,
  const size_t freq_dimension,
  const size_t f1dot_dimension,
  const double age,
  const double min_braking,
  const double max_braking
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(freq_dimension < f1dot_dimension, XLAL_EINVAL);
  XLAL_CHECK(age > 0.0, XLAL_EINVAL);
  XLAL_CHECK(min_braking > 1.0, XLAL_EINVAL);
  XLAL_CHECK(max_braking > 1.0, XLAL_EINVAL);
  XLAL_CHECK(min_braking <= max_braking, XLAL_EINVAL);

  // Allocate memory
  const size_t info_len = sizeof(F1DotAgeBrakingBoundInfo);
  F1DotAgeBrakingBoundInfo* info_lower = XLALMalloc(info_len);
  XLAL_CHECK(info_lower != NULL, XLAL_ENOMEM);
  F1DotAgeBrakingBoundInfo* info_upper = XLALMalloc(info_len);
  XLAL_CHECK(info_upper != NULL, XLAL_ENOMEM);

  // Set the parameter-space bound
  info_lower->freq_dim = info_upper->freq_dim = freq_dimension;
  info_lower->scale = -1.0 / ((min_braking - 1.0) * age);
  info_upper->scale = -1.0 / ((max_braking - 1.0) * age);
  XLAL_CHECK(XLALSetLatticeTilingBound(tiling, f1dot_dimension, F1DotAgeBrakingBound, info_len, info_lower, info_upper) == XLAL_SUCCESS, XLAL_EFAILED);

  return XLAL_SUCCESS;

}

typedef struct {
  size_t freq_dim;
  size_t f1dot_dim;
  double scale;
} F2DotBrakingBoundInfo;

static double F2DotBrakingBound(
  const void* data,
  const size_t dim UNUSED,
  const gsl_vector* point
  )
{

  // Get bounds info
  const F2DotBrakingBoundInfo* info = (const F2DotBrakingBoundInfo*)data;

  // Get current values of frequency and first spindown
  const double freq = gsl_vector_get(point, info->freq_dim);
  const double f1dot = gsl_vector_get(point, info->f1dot_dim);

  // Return second spindown bounds
  return info->scale * f1dot * f1dot / freq;

}

int XLALSetLatticeTilingF2DotBrakingBound(
  LatticeTiling* tiling,
  const size_t freq_dimension,
  const size_t f1dot_dimension,
  const size_t f2dot_dimension,
  const double min_braking,
  const double max_braking
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(freq_dimension < f1dot_dimension, XLAL_EINVAL);
  XLAL_CHECK(f1dot_dimension < f2dot_dimension, XLAL_EINVAL);
  XLAL_CHECK(min_braking > 0.0, XLAL_EINVAL);
  XLAL_CHECK(max_braking > 0.0, XLAL_EINVAL);
  XLAL_CHECK(min_braking <= max_braking, XLAL_EINVAL);

  // Allocate memory
  const size_t info_len = sizeof(F2DotBrakingBoundInfo);
  F2DotBrakingBoundInfo* info_lower = XLALMalloc(info_len);
  XLAL_CHECK(info_lower != NULL, XLAL_ENOMEM);
  F2DotBrakingBoundInfo* info_upper = XLALMalloc(info_len);
  XLAL_CHECK(info_upper != NULL, XLAL_ENOMEM);

  // Set the parameter-space bound
  info_lower->freq_dim = info_upper->freq_dim = freq_dimension;
  info_lower->f1dot_dim = info_upper->f1dot_dim = f1dot_dimension;
  info_lower->scale = min_braking;
  info_upper->scale = max_braking;
  XLAL_CHECK(XLALSetLatticeTilingBound(tiling, f2dot_dimension, F2DotBrakingBound, info_len, info_lower, info_upper) == XLAL_SUCCESS, XLAL_EFAILED);

  return XLAL_SUCCESS;

}
