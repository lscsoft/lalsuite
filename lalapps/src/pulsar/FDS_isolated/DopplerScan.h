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

\noindent This header provides two trivial functions to divide real
numbers.  It exists primarily to demonstrate documentation and coding
standards for LAL headers.

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
#define DOPPLERSCANH_ENOINIT 		2
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
#define DOPPLERSCANH_MSGENOINIT 	"DopplerScan module has not been initialized"
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

typedef enum
{
  LAL_METRIC_NONE = 0,
  LAL_METRIC_PTOLE_ANALYTIC,	/* analytic ptolemaic approx for the metric */
  LAL_METRIC_COHERENT_PTOLE,	/* numerical metric using ptole-approximation in timing */
  LAL_METRIC_COHERENT_EPHEM,	/* numerical exact metric using ephemeris-timing */
  LAL_METRIC_LAST
} PulsarMetricType;

  /* NOTE: if you change these, you also need to adapt the help-string + cfs_example-file !! */
typedef enum
{
  GRID_FLAT,			/* "flat" sky-grid: fixed step-size (dAlpha,dDelta) */
  GRID_ISOTROPIC,		/* approximately isotropic sky-grid */
  GRID_METRIC,			/* generate grid using a 2D sky-metric */
  GRID_FILE,			/* read grid from a file */
  GRID_LAST
} DopplerGridType;

/* this structure is handed over to InitDopplerScan() */  
typedef struct {
  DopplerGridType gridType;	/* which type of grid to generate */  
  PulsarMetricType metricType; 	/* which metric to use if GRID_METRIC */

  REAL8 dAlpha;			/* step-sizes for GRID_FLAT */
  REAL8 dDelta;

  REAL8 metricMismatch;		/* for GRID_METRIC and GRID_ISOTROPIC */

  LIGOTimeGPS obsBegin; 	/* start-time of time-series */
  REAL8 obsDuration;		/* length of time-series in seconds */
  REAL8 fmax; 			/* max frequency of search */
  LALDetector *Detector; 	/* Our detector*/
  EphemerisData *ephemeris;	/* ephemeris for "exact" metric */

  CHAR *skyRegion;		/* string containing a list of sky-positions describing a sky-region */
  CHAR *skyGridFile;		/* file containing a sky-grid (list of points) for GRID_FILE */
} DopplerScanInit;


typedef struct {
  SkyPosition skypos;
  REAL8Vector spindowns;
  BOOLEAN finished;
} DopplerPosition;


typedef struct {
  UINT4 numVertices;
  SkyPosition *vertices;
  SkyPosition lowerLeft;
  SkyPosition upperRight;
} SkyRegion;

/* general scan-grid */
typedef struct tagDopplerScanGrid {
  REAL8 freq;
  REAL8 alpha;
  REAL8 delta;
  REAL8Vector spindowns;
  struct tagDopplerScanGrid *next;
} DopplerScanGrid;

/* this structure reflects the internal state of DopplerScan */
typedef struct {
  INT2 state;  			/* idle, ready or finished */

  SkyRegion skyRegion; 		/* polygon (and bounding square) defining sky-region  */

  UINT4 numGridPoints;		/* how many grid-points */
  DopplerScanGrid *grid; 	/* head of linked list of nodes */  
  DopplerScanGrid *gridNode;	/* pointer to current grid-node in grid */
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

/********************************************************** <lalLaTeX>
\newpage\input{LALSampleTestC}
******************************************************* </lalLaTeX> */

#ifdef  __cplusplus
}
#endif  
/* C++ protection. */

#endif  /* Double-include protection. */
