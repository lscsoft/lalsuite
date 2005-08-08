/*
 * Copyright (C) 2004, 2005 Reinhard Prix
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
 * \date 2004, 2005
 * \file 
 * \brief Header file defining the API for DopplerScan.
 *
 * $Id$
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

NRCSID( DOPPLERSCANH, "$Id$" );

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

/*---------- external types ----------*/

/** Different 'states' a Doppler-scan can be in */
enum {
  STATE_IDLE = 0,   	/**< not initialized yet */
  STATE_READY,		/**< initialized and ready */
  STATE_FINISHED,	/**< all templates have been read */
  STATE_LAST
};

/** different types of Sky-grids: */
typedef enum
{
  GRID_FLAT,			/**< "flat" sky-grid: fixed step-size (dAlpha,dDelta) */
  GRID_ISOTROPIC,		/**< approximately isotropic sky-grid */
  GRID_METRIC,			/**< generate grid using a 2D sky-metric */
  GRID_FILE,			/**< read grid from a file */
  GRID_LAST
} DopplerGridType;

/** structure holding one point in (phase-) parameter-space */
typedef struct {
  REAL8 Alpha;			/**< longitude in Radians, EquatorialCoordinates */
  REAL8 Delta;			/**< latitude in Radians, EquatorialCoordinates */
  REAL8 Freq;			/**< frequency */
  REAL8 f1dot;			/**< first frequency-derivative (spindown) */
} DopplerPosition;

/** Structure describing a region in paramter-space (a,d,f,f1dot,..). 
 *  Currently this is simply a direct product of skyRegion x FreqBand x f1dotBand.
 */
typedef struct {
  CHAR *skyRegionString;	/**< sky-region string '(a1,d1), (a2,d2), ..' */
  REAL8 Freq;			/**< first point of Frequency-interval */
  REAL8 FreqBand;		/**< width of frequency-interval */
  REAL8 f1dot;			/**< first spindown-value */
  REAL8 f1dotBand;		/**< width of spindown-interval */
} DopplerRegion;

/** structure describing a polygon-region in the sky */
typedef struct {
  UINT4 numVertices;		/**< number of polygon-vertices */
  SkyPosition *vertices;	/**< array of vertices */
  SkyPosition lowerLeft;	/**< lower-left point of bounding square */
  SkyPosition upperRight;	/**< upper-right point of bounding square */
} SkyRegion;

/** general scan-grid */
typedef struct tagDopplerScanGrid {
  REAL8 Alpha;
  REAL8 Delta;
  REAL8 Freq;
  REAL8 f1dot;
  struct tagDopplerScanGrid *next;
} DopplerScanGrid;

/** a "sky-ellipse", described by the two major axes and it's angle wrt x-axis */
typedef struct {
  REAL8 smajor;
  REAL8 sminor;
  REAL8 angle;
} MetricEllipse;

/** initialization-structure passed to InitDopplerScan() */  
typedef struct {
  DopplerRegion searchRegion;	/**< parameter-space region to search over */
  DopplerGridType gridType;	/**< which type of skygrid to generate */  
  LALPulsarMetricType metricType; /**< which metric to use if GRID_METRIC */
  REAL8 dAlpha;			/**< step-sizes for GRID_FLAT */
  REAL8 dDelta;
  REAL8 metricMismatch;		/**< for GRID_METRIC and GRID_ISOTROPIC */
  LIGOTimeGPS obsBegin; 	/**< GPS start-time of time-series */
  REAL8 obsDuration;		/**< length of time-series in seconds */
  BOOLEAN projectMetric;	/**< project the metric orthogonal to Freq? */
  LALDetector *Detector; 	/**< Current detector */
  EphemerisData *ephemeris;	/**< ephemeris-data for "exact" metric */
  CHAR *skyGridFile;		/**< file containing a sky-grid (list of points) if GRID_FILE */
} DopplerScanInit;

/** this structure reflects the current state of DopplerScan */
typedef struct {
  INT2 state;  			/**< idle, ready or finished */
  SkyRegion skyRegion; 		/**< polygon (and bounding square) defining sky-region  */
  UINT4 numGridPoints;		/**< how many grid-points */
  REAL8 dFreq;			/**< stepsize in frequency */
  REAL8 df1dot;			/**< stepsize in spindown-value f1dot */
  DopplerScanGrid *grid; 	/**< head of linked list of skygrid nodes */  
  DopplerScanGrid *gridNode;	/**< pointer to current grid-node in skygrid */
} DopplerScanState;

/*---------- Global variables ----------*/
/* some empty structs for initializations */
extern const DopplerScanGrid empty_DopplerScanGrid;
extern const DopplerScanState empty_DopplerScanState;
extern const DopplerScanInit empty_DopplerScanInit;
extern const DopplerPosition empty_DopplerPosition;
extern const DopplerRegion empty_DopplerRegion;

/*---------- external prototypes [API] ----------*/

void InitDopplerScan(LALStatus *, DopplerScanState *scan, const DopplerScanInit *init);
void NextDopplerPos(LALStatus *, DopplerPosition *pos, DopplerScanState *scan);
void FreeDopplerScan(LALStatus *, DopplerScanState *scan);

void writeSkyGridFile(LALStatus *, const DopplerScanGrid *grid, const CHAR *fname, 
		      const DopplerScanInit *init);
void ParseSkyRegionString (LALStatus *, SkyRegion *region, const CHAR *input);
void SkySquare2String (LALStatus *, CHAR **string, REAL8 Alpha, REAL8 Delta, 
		       REAL8 AlphaBand, REAL8 DeltaBand);

void getGridSpacings(LALStatus *, DopplerPosition *spacings, DopplerPosition gridpoint, 
		     const DopplerScanInit *params);
void getMCDopplerCube (LALStatus *, DopplerRegion *cube, DopplerPosition signal, UINT4 PointsPerDim, 
		       const DopplerScanInit *params);

void getMetricEllipse(LALStatus *, 
		      MetricEllipse *ellipse, 
		      REAL8 mismatch, 
		      const REAL8Vector *metric, 
		      UINT4 dim0);

#ifdef  __cplusplus
}
#endif  
/* C++ protection. */

#endif  /* Double-include protection. */
