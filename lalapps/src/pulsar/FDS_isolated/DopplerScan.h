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

#define DOPPLERSCANH_MSGENULL 		"Arguments contained an unexpected null pointer"
#define DOPPLERSCANH_MSGENOINIT 	"NextDopplerPos() called without prior initialization by InitDopplerScan()"
#define DOPPLERSCANH_MSGESYS		"System call failed (probably file IO)"
#define DOPPLERSCANH_MSGE2DSKY		"Either need one sky-point or a polygon. (2 sky-points where given)"
#define DOPPLERSCANH_MSGE2DSTEP		"If not using the metric, you need to specify _both_ dDelta and dAlpha"
#define DOPPLERSCANH_MSGEGRIDCRPT	"Unexpected NULL in grid-list. This points to a bug in the code... "
#define DOPPLERSCANH_MSGESKYPARAM	"Invalid sky region! We need 0<= alpha < 2Pi and -Pi/2 <= delta <= PI/2"
#define DOPPLERSCANH_MSGEMETRIC		"Unknown type of metric specified. 1=PtoleMetric, 2=CoherentMetric"
#define DOPPLERSCANH_MSGENONULL		"Output pointer is not NULL"
#define DOPPLERSCANH_MSGEMEM		"Out of memory"
#define DOPPLERSCANH_MSGESKYREGION	"Could not parse sky-region properly"

/*************************************************** </lalErrTable> */

typedef enum
{
  LAL_METRIC_NONE = 0,
  LAL_METRIC_PTOLE,		/* analytic ptolemaic approx for the metric */
  LAL_METRIC_COHERENT_PTOLE,	/* numerical 'exact' metric using ptole-timing */
  LAL_METRIC_COHERENT_EXACT,	/* numerical exact metric using ephemeris-timing */
  LAL_METRIC_PSEUDO_ISOTROPIC
} LALMetricType;

/* this structure is handed over to InitDopplerScan() */  
typedef struct {
  REAL8 dAlpha;		/* step-sizes for manual stepping */
  REAL8 dDelta;
  
  LALMetricType metricType;   	/* 0 = manual, 1 = PtoleMetric, 2 = CohMetric_ptole, .. */
  REAL8 metricMismatch;
  LIGOTimeGPS obsBegin; /* start-time of time-series */
  REAL8 obsDuration;	/* length of time-series in seconds */
  REAL8 fmax; 		/* max frequency of search */
  LALDetector *Detector; /* Our detector*/
  EphemerisData *ephemeris;	/* used by Ephemeris-based metric */
  CHAR *skyRegion;	/* string containing a list of sky-positions (ra,dec) which describe a sky-region */
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

  REAL8 dAlpha;			/* step-sizes for manual stepping */
  REAL8 dDelta;

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

void InitDopplerScan( LALStatus *stat, DopplerScanState *scan, DopplerScanInit init);
void NextDopplerPos ( LALStatus *stat, DopplerPosition *pos, DopplerScanState *scan);
void FreeDopplerScan (LALStatus *stat, DopplerScanState *scan);

void ParseSkyRegion (LALStatus *stat, SkyRegion *region, const CHAR *input);

void LALMetricWrapper (LALStatus *stat, REAL8Vector *metric, PtoleMetricIn *input, LALMetricType metricType);

/********************************************************** <lalLaTeX>
\newpage\input{LALSampleTestC}
******************************************************* </lalLaTeX> */

#ifdef  __cplusplus
}
#endif  
/* C++ protection. */

#endif  /* Double-include protection. */
