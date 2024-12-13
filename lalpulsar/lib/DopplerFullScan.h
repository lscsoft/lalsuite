/*
 * Copyright (C) 2007, 2008, 2012 Karl Wette
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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */

/**
 * \author Reinhard Prix, Karl Wette
 * \date 2006
 * \file
 * \brief Header file defining the API for DopplerFullScan.
 *
 */

#ifndef _DOPPLERFULLSCAN_H  /* Double-include protection. */
#define _DOPPLERFULLSCAN_H

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/*---------- INCLUDES ----------*/
#include <lal/LALDatatypes.h>
#include <lal/SkyCoordinates.h>
#include <lal/PtoleMetric.h>
#include <lal/LALBarycenter.h>
#include <lal/PulsarDataTypes.h>
#include <lal/ComputeFstat.h>
#include <lal/DopplerScan.h>
#include <lal/LatticeTiling.h>

/*---------- DEFINES ----------*/

/*---------- external types ----------*/

/* ==================== FULL MULTIDIMENSIONAL-GRID types ==================== */
/**
 * Structure describing a region in paramter-space (a,d,f,f1dot,..).
 * Currently this is simply a direct product of skyRegion x FreqBand x f1dotBand.
 */

/** initialization struct for full InitDopplerScan() [UNDER CONSTRUCTION] */
#ifdef SWIG /* SWIG interface directives */
SWIGLAL( IMMUTABLE_MEMBERS( tagDopplerFullScanInit, Detector, ephemeris, gridFile ) );
#endif /* SWIG */
typedef struct tagDopplerFullScanInit {
  DopplerRegion searchRegion;           /**< Doppler-space region to be covered + scanned */
  DopplerGridType gridType;             /**< which type of grid to generate */
  LALPulsarMetricType metricType;       /**< which metric to use if GRID_METRIC */
  BOOLEAN projectMetric;                /**< project metric on f=const subspace */
  PulsarDopplerParams stepSizes;        /**< user-settings for stepsizes if GRID_FLAT */
  REAL8 metricMismatch;                 /**< for GRID_METRIC and GRID_ISOTROPIC */
  LIGOTimeGPS startTime;                /**< start-time of the observation */
  REAL8 Tspan;                          /**< total time spanned by observation */
  const LALDetector *Detector;          /**< Current detector */
  const EphemerisData *ephemeris;       /**< ephemeris-data for numerical metrics */
  const CHAR *gridFile;                 /**< filename for sky-grid or full-grid if GRID_FILE_SKYGRID or GRID_FILE_FULLGRID */
  REAL8 extraArgs[3];                   /**< extra grid-specific setup arguments */
} DopplerFullScanInit;

/** opaque type to reflects the current state of a full multi-dimensional DopplerScan */
typedef struct tagDopplerFullScanState DopplerFullScanState;    /* opaque type */


/*---------- Global variables ----------*/
/* some empty structs for initializations */

/*---------- external prototypes [API] ----------*/

/* ----- full-fledged multi-dimensional Doppler scanner functions ----- */
DopplerFullScanState *XLALInitDopplerFullScan( const DopplerFullScanInit *init );

int  XLALNextDopplerPos( PulsarDopplerParams *pos, DopplerFullScanState *scan );
REAL8 XLALNumDopplerTemplates( DopplerFullScanState *scan );
UINT8 XLALNumDopplerPointsAtDimension( DopplerFullScanState *scan, const size_t dim );
int XLALGetDopplerSpinRange( PulsarSpinRange *spinRange, const DopplerFullScanState *scan );
void XLALDestroyDopplerFullScan( DopplerFullScanState *scan );

/* ----- variout utility functions ----- */

///
/// Set a first spindown bound derived from spindown age and braking indices
///
int XLALSetLatticeTilingF1DotAgeBrakingBound(
  LatticeTiling *tiling,                ///< [in] Lattice tiling
  const size_t freq_dimension,          ///< [in] Frequency dimension
  const size_t f1dot_dimension,         ///< [in] First spindown dimension
  const double age,                     ///< [in] Spindown age
  const double min_braking,             ///< [in] Minimum braking index
  const double max_braking              ///< [in] Maximum braking index
);

///
/// Set a second spindown bound derived from braking indices
///
int XLALSetLatticeTilingF2DotBrakingBound(
  LatticeTiling *tiling,                ///< [in] Lattice tiling
  const size_t freq_dimension,          ///< [in] Frequency dimension
  const size_t f1dot_dimension,         ///< [in] First spindown dimension
  const size_t f2dot_dimension,         ///< [in] Second spindown dimension
  const double min_braking,             ///< [in] Minimum braking index
  const double max_braking              ///< [in] Maximum braking index
);

///
/// Return the step size of the spindown lattice tiling in a given dimension, or 0 for non-tiled dimensions.
///
REAL8 XLALGetDopplerLatticeTilingStepSize(
  DopplerFullScanState *scan,           ///< [in] Doppler scan state object
  const size_t dim                      ///< [in] Dimension of which to return step size
);

// ========== deprecated LAL functions ==========
void InitDopplerFullScan( LALStatus *, DopplerFullScanState **scanState, const DopplerFullScanInit *init );

#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
