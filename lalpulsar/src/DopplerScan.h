/*
 * Copyright (C) 2007, 2008 Karl Wette
 * Copyright (C) 2004, 2005, 2006 Reinhard Prix
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
 * \date 2004, 2005, 2006, 2007, 2008
 * \file
 * \brief Header file defining the API for DopplerScan.
 *
 */

#ifndef _DOPPLERSCAN_H  /* Double-include protection. */
#define _DOPPLERSCAN_H

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/*---------- INCLUDES ----------*/
#include <lal/LALDatatypes.h>
#include <lal/SkyCoordinates.h>
#include <lal/PtoleMetric.h>
#include <lal/StackMetric.h>
#include <lal/LALBarycenter.h>
#include <lal/PulsarDataTypes.h>
#include <lal/ComputeFstat.h>

/*---------- DEFINES ----------*/

#define DOPPLERSCANH_ENULL 		1
#define DOPPLERSCANH_ENOTREADY 		2
#define DOPPLERSCANH_ESYS      		3
#define DOPPLERSCANH_E2DSKY		4
#define DOPPLERSCANH_E2DSTEP		5
#define DOPPLERSCANH_EGRIDCRPT		6
#define DOPPLERSCANH_ESKYPARAM		7
#define DOPPLERSCANH_EMETRICTYPE	8
#define DOPPLERSCANH_ENONULL		9
#define DOPPLERSCANH_EMEM		10
#define DOPPLERSCANH_ESKYREGION		11
#define DOPPLERSCANH_EINPUT		12
#define DOPPLERSCANH_ENEGMETRIC		13
#define DOPPLERSCANH_EXLAL		14


#define DOPPLERSCANH_MSGENULL 		"Arguments contained an unexpected null pointer"
#define DOPPLERSCANH_MSGENOTREADY 	"Doppler scan is uninitialized or has finished"
#define DOPPLERSCANH_MSGESYS		"System call failed (probably file IO)"
#define DOPPLERSCANH_MSGE2DSKY		"Either need one sky-point or a polygon. (2 sky-points where given)"
#define DOPPLERSCANH_MSGE2DSTEP		"If not using the metric, you need to specify _both_ dDelta and dAlpha"
#define DOPPLERSCANH_MSGEGRIDCRPT	"Unexpected NULL in grid-list. This points to a bug in the code... "
#define DOPPLERSCANH_MSGESKYPARAM	"Invalid sky region! We need 0<= alpha < 2Pi and -Pi/2 <= delta <= PI/2"
#define DOPPLERSCANH_MSGEMETRICTYPE	"Unknown type of metric specified."
#define DOPPLERSCANH_MSGENONULL		"Output pointer is not NULL"
#define DOPPLERSCANH_MSGEMEM		"Out of memory"
#define DOPPLERSCANH_MSGESKYREGION	"Could not parse sky-region correctly"
#define DOPPLERSCANH_MSGEINPUT		"Invald input parameter"
#define DOPPLERSCANH_MSGENEGMETRIC	"Negative metric encountered"
#define DOPPLERSCANH_MSGEXLAL		"XLAL call failed"

/*---------- external types ----------*/

/** Different 'states' a Doppler-scan can be in */
typedef enum {
  STATE_IDLE = 0,   	/**< not initialized yet */
  STATE_READY,		/**< initialized and ready */
  STATE_FINISHED,	/**< all templates have been read */
  STATE_LAST
} scan_state_t;

/** different types of grids: */
typedef enum
{
  /* ----- factored grid-types: sky x f0dot x f1dot x f2dot x f3dot  */
  GRID_FLAT 		= 0,		/**< "flat" sky-grid: fixed step-size (dAlpha,dDelta) */
  GRID_ISOTROPIC	= 1,		/**< approximately isotropic sky-grid */
  GRID_METRIC		= 2,		/**< generate grid using a 2D sky-metric */
  GRID_FILE_SKYGRID	= 3,		/**< read skygrid from a file */
  GRID_METRIC_SKYFILE	= 4,		/**< 'hybrid' grid-construction: use skygrid from file, metric for others */
  GRID_SKY_LAST,			/**< end-marker for factored grid types */
  /* ----- full multi-dim grid-types ----- */
  GRID_FILE_FULLGRID	= 6,		/**< load the full D-dim grid from a file */
  GRID_METRIC_LATTICE	= 7,		/**< 'optimal' covering using An*-lattice and flat metric */
  GRID_SPINDOWN_SQUARE  = 8,            /**< spindown tiling for a single sky position and square parameter space */
  GRID_SPINDOWN_AGEBRK  = 9,            /**< spindown tiling for a single sky position and non-square parameter space
					   defined by spindown age and braking indices */
  /* ----- end-marker ----- */
  GRID_LAST
} DopplerGridType;

/** structure describing a polygon-region in the sky */
typedef struct tagSkyRegion {
  UINT4 numVertices;		/**< number of polygon-vertices */
  SkyPosition *vertices;	/**< array of vertices */
  SkyPosition lowerLeft;	/**< lower-left point of bounding square */
  SkyPosition upperRight;	/**< upper-right point of bounding square */
} SkyRegion;

typedef struct tagDopplerRegion {
  CHAR *skyRegionString;	/**< sky-region string '(a1,d1), (a2,d2), ..' */
  LIGOTimeGPS refTime;		/** reference time of definition of Doppler parameters */
  PulsarSpins fkdot;		/**< first points of spin-intervals */
  PulsarSpins fkdotBand;	/**< spin-intervals */
} DopplerRegion;

/* ==================== SKYGRID-ONLY types ==================== */
/** sky grid */
typedef struct tagDopplerSkyGrid {
  REAL8 Alpha;
  REAL8 Delta;
  struct tagDopplerSkyGrid *next;
} DopplerSkyGrid;

/** initialization-structure passed to InitDopplerSkyScan() */
#ifdef SWIG /* SWIG interface directives */
SWIGLAL(STRUCT_IMMUTABLE(tagDopplerSkyScanInit, skyGridFile));
#endif /* SWIG */
typedef struct tagDopplerSkyScanInit {
  CHAR *skyRegionString;	/**< sky-region to search: format polygon '(a1,d1), (a2,d2), ..' */
  REAL8 Freq;			/**< Frequency for which to build the skyGrid */
  DopplerGridType gridType;	/**< which type of skygrid to generate */
  LALPulsarMetricType metricType; /**< which metric to use if GRID_METRIC */
  REAL8 dAlpha, dDelta;		/**< sky step-sizes for GRID_FLAT and GRID_ISOTROPIC */
  REAL8 metricMismatch;		/**< for GRID_METRIC */
  LIGOTimeGPS obsBegin; 	/**< GPS start-time of time-series */
  REAL8 obsDuration;		/**< length of time-series in seconds */
  BOOLEAN projectMetric;	/**< project the metric orthogonal to Freq? */
  const LALDetector *Detector; 	/**< Current detector */
  const EphemerisData *ephemeris;/**< ephemeris-data for "exact" metric */
  const CHAR *skyGridFile;		/**< file containing a sky-grid (list of points) if GRID_FILE */
  UINT4 numSkyPartitions;	/**< number of (roughly) equal partitions to split sky-grid into */
  UINT4 partitionIndex;		/**< index of requested sky-grid partition: in [0, numPartitions - 1] */
} DopplerSkyScanInit;

/** this structure reflects the current state of a DopplerSkyScan */
typedef struct tagDopplerSkyScanState {
  scan_state_t state;  			/**< idle, ready or finished */
  SkyRegion skyRegion; 		/**< polygon (and bounding square) defining sky-region  */
  UINT4 numSkyGridPoints;	/**< how many skygrid-points */
  PulsarSpins dfkdot;		/**< fixed-size steps in spins */
  DopplerSkyGrid *skyGrid; 	/**< head of linked list of skygrid nodes */
  DopplerSkyGrid *skyNode;	/**< pointer to current grid-node in skygrid */
} DopplerSkyScanState;

/** a "sky-ellipse", described by the two major axes and it's angle wrt x-axis */
typedef struct tagMetricEllipse {
  REAL8 smajor;
  REAL8 sminor;
  REAL8 angle;
} MetricEllipse;

/*---------- Global variables ----------*/
/* some empty structs for initializations */
extern const DopplerSkyGrid empty_DopplerSkyGrid;
extern const DopplerSkyScanState empty_DopplerSkyScanState;
extern const DopplerSkyScanInit empty_DopplerSkyScanInit;
extern const DopplerRegion empty_DopplerRegion;
extern const SkyRegion empty_SkyRegion;

/*---------- external prototypes [API] ----------*/

/* ------ functions to handle factored grids: 'sky x freq x f1dot...' covering ----- */
void InitDopplerSkyScan(LALStatus *, DopplerSkyScanState *skyScan, const DopplerSkyScanInit *init);
int  XLALNextDopplerSkyPos( PulsarDopplerParams *pos, DopplerSkyScanState *skyScan);
void FreeDopplerSkyScan(LALStatus *, DopplerSkyScanState *skyScan);

void writeSkyGridFile(LALStatus *, const DopplerSkyGrid *grid, const CHAR *fname );

/* ----- various utility functions ----- */
void ParseSkyRegionString (LALStatus *, SkyRegion *region, const CHAR *input);
void SkySquare2String (LALStatus *, CHAR **string, REAL8 Alpha, REAL8 Delta, REAL8 AlphaBand, REAL8 DeltaBand);

void getMCDopplerCube (LALStatus *, DopplerRegion *cube, PulsarDopplerParams lal_signal, UINT4 PointsPerDim, const DopplerSkyScanInit *params);
void getMetricEllipse(LALStatus *, MetricEllipse *ellipse, REAL8 mismatch, const REAL8Vector *metric, UINT4 dim0);

int fprintfDopplerParams ( FILE *fp, const PulsarDopplerParams *params );

DopplerSkyGrid *XLALEquiPartitionSkygrid ( const DopplerSkyGrid *skygrid, UINT4 jPart, UINT4 numPartitions );

#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
