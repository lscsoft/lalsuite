/************************************ <lalVerbatim file="DopplerScanHV">
Author: Prix, Reinhard
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\section{Header \texttt{DopplerScan.h}}
\label{s:DopplerScan.h}

Header file for DopplerScan

\subsection*{Synopsis}
\begin{verbatim}
#include "DopplerScan.h"
\end{verbatim}

******************************************************* </lalLaTeX> */

#ifndef _DOPPLERSCAN_H  /* Double-include protection. */
#define _DOPPLERSCAN_H

#include <lal/LALDatatypes.h>
#include <lal/SkyCoordinates.h>
#include <lal/PtoleMetric.h>
#include <lal/StackMetric.h>
#include <lal/LALBarycenter.h>

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

NRCSID( DOPPLERSCANH, "$Id$" );

/********************************************************** <lalLaTeX>
\subsection*{Error codes}
</lalLaTeX>
***************************************************** <lalErrTable> */
#define DOPPLERSCANH_ENULL 		1
#define DOPPLERSCANH_ENOTREADY 		2
#define DOPPLERSCANH_ESYS      		3
#define DOPPLERSCANH_E2DSKY		4
#define DOPPLERSCANH_E2DSTEP		5
#define DOPPLERSCANH_EGRIDCRPT		6
#define DOPPLERSCANH_ESKYPARAM		7
#define DOPPLERSCANH_EMETRIC		8
#define DOPPLERSCANH_ENONULL		9
#define DOPPLERSCANH_EMEM		10
#define DOPPLERSCANH_ESKYREGION		11
#define DOPPLERSCANH_EINPUT		12

#define DOPPLERSCANH_MSGENULL 		"Arguments contained an unexpected null pointer"
#define DOPPLERSCANH_MSGENOTREADY 	"Doppler scan is uninitialized or has finished"
#define DOPPLERSCANH_MSGESYS		"System call failed (probably file IO)"
#define DOPPLERSCANH_MSGE2DSKY		"Either need one sky-point or a polygon. (2 sky-points where given)"
#define DOPPLERSCANH_MSGE2DSTEP		"If not using the metric, you need to specify _both_ dDelta and dAlpha"
#define DOPPLERSCANH_MSGEGRIDCRPT	"Unexpected NULL in grid-list. This points to a bug in the code... "
#define DOPPLERSCANH_MSGESKYPARAM	"Invalid sky region! We need 0<= alpha < 2Pi and -Pi/2 <= delta <= PI/2"
#define DOPPLERSCANH_MSGEMETRIC		"Unknown type of metric specified."
#define DOPPLERSCANH_MSGENONULL		"Output pointer is not NULL"
#define DOPPLERSCANH_MSGEMEM		"Out of memory"
#define DOPPLERSCANH_MSGESKYREGION	"Could not parse sky-region correctly"
#define DOPPLERSCANH_MSGEINPUT		"Invald input parameter"

/*************************************************** </lalErrTable> */

/** a Doppler-scan can be in one of the following states */
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

/** initialization-structure passed to InitDopplerScan() */  
typedef struct {
  DopplerGridType gridType;	/**< which type of skygrid to generate */  
  LALPulsarMetricType metricType; /**< which metric to use if GRID_METRIC */

  REAL8 dAlpha;			/**< step-sizes for GRID_FLAT */
  REAL8 dDelta;
  REAL8 metricMismatch;		/**< for GRID_METRIC and GRID_ISOTROPIC */
  LIGOTimeGPS obsBegin; 	/**< start-time of time-series */
  REAL8 obsDuration;		/**< length of time-series in seconds */
  REAL8 fmax; 			/**< max frequency of search */
  LALDetector *Detector; 	/**< Our detector*/
  EphemerisData *ephemeris;	/**< ephemeris for "exact" metric */
  CHAR *skyRegion;		/**< list of sky-positions describing a sky-region */
  CHAR *skyGridFile;		/**< file containing a sky-grid (list of points) for GRID_FILE */
} DopplerScanInit;


typedef struct {
  SkyPosition skypos;
  REAL8 freq;
  REAL8Vector spindowns;
} DopplerPosition;


typedef struct {
  UINT4 numVertices;
  SkyPosition *vertices;
  SkyPosition lowerLeft;
  SkyPosition upperRight;
} SkyRegion;

/** general scan-grid */
typedef struct tagDopplerScanGrid {
  REAL8 alpha;
  REAL8 delta;
  REAL8 freq;
  REAL8Vector spindowns;
  struct tagDopplerScanGrid *next;
} DopplerScanGrid;

/** this structure reflects the current state of DopplerScan */
typedef struct {
  INT2 state;  			/**< idle, ready or finished */
  SkyRegion skyRegion; 		/**< polygon (and bounding square) defining sky-region  */
  UINT4 numGridPoints;		/**< how many grid-points */
  REAL8 dFreq;			/**< stepsize for frequency */
  REAL8 df1dot;			/**< stepsize for first spindown-value f1dot */
  DopplerScanGrid *grid; 	/**< head of linked list of skygrid nodes */  
  DopplerScanGrid *gridNode;	/**< pointer to current grid-node in skygrid */
} DopplerScanState;
  
/********************************************************** <lalLaTeX>
\vfill{\footnotesize\input{DopplerScanHV}}
\newpage\input{DopplerScanC}
******************************************************* </lalLaTeX> */

/* Function prototypes */

void InitDopplerScan( LALStatus *stat, DopplerScanState *scan, const DopplerScanInit *init);
void NextDopplerPos ( LALStatus *stat, DopplerPosition *pos, DopplerScanState *scan);
void FreeDopplerScan (LALStatus *stat, DopplerScanState *scan);

void writeSkyGridFile (LALStatus *stat, const DopplerScanGrid *grid, const CHAR *fname, const DopplerScanInit *init);
void ParseSkyRegion (LALStatus *stat, SkyRegion *region, const CHAR *input);
void LALMetricWrapper(LALStatus *stat, REAL8Vector **metric, PtoleMetricIn *params);
void get_dfkdot( LALStatus *lstat, REAL8Vector *dfkdot, SkyPosition skypos, const DopplerScanInit *params);


/********************************************************** <lalLaTeX>
\newpage\input{LALSampleTestC}
******************************************************* </lalLaTeX> */

#ifdef  __cplusplus
}
#endif  
/* C++ protection. */

#endif  /* Double-include protection. */
