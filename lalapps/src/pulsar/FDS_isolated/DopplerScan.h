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
#include <lal/TwoDMesh.h>
#include <lal/PtoleMetric.h>
#include <lal/StackMetric.h>

#ifdef  __cplusplus   /* C++ protection. */
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
#define DOPPLERSCANH_MSGE2DSKY		"Cannot deal with 1-D sky-regions"
#define DOPPLERSCANH_MSGE2DSTEP		"If not using the metric, you need to specify _both_ dDelta and dAlpha"
#define DOPPLERSCANH_MSGEGRIDCRPT	"Unexpected NULL in grid-list. This points to a bug in the code... "
#define DOPPLERSCANH_MSGESKYPARAM	"Invalid sky region! We need 0<= alpha < 2Pi and -Pi/2 <= delta <= PI/2"
#define DOPPLERSCANH_MSGEMETRIC		"Unknown type of metric specified. 1=PtoleMetric, 2=CoherentMetric"
#define DOPPLERSCANH_MSGENONULL		"Output pointer is not NULL"
#define DOPPLERSCANH_MSGEMEM		"Out of memory"
#define DOPPLERSCANH_MSGESKYREGION	"Could not parse sky-region properly"

/*************************************************** </lalErrTable> */

/* this structure is handed over to InitDopplerScan() */  
typedef struct {
  REAL8 Alpha;
  REAL8 dAlpha;
  REAL8 AlphaBand;
  REAL8 Delta;
  REAL8 dDelta;
  REAL8 DeltaBand;
  
  INT2 useMetric;   	/* 0 = manual, 1 = PtoleMetric, 2 = CoherentMetric */
  REAL8 metricMismatch;
  LIGOTimeGPS obsBegin; /* start-time of time-series */
  REAL4 obsDuration;	/* length of time-series in seconds */
  REAL4 fmax; 		/* max frequency of search */
  LALDetector Detector; /* Our detector*/
  BOOLEAN flipTiling;	/* use non-standard internal grid order? ORDER_DELTA_ALPHA */

  CHAR *skyRegion;	/* string containing a list of sky-positions (ra,dec) which describe a sky-region */
} DopplerScanInit;

typedef struct {
  SkyPosition skypos;
  REAL8Vector spindowns;
  BOOLEAN finished;
} DopplerPosition;


typedef struct {
  UINT4 length;
  SkyPosition *data;
} SkyPositionVector;


/* this structure reflects the internal state of DopplerScan */
typedef struct {
  INT2 state;  /* idle, manual, grid, or finished */

  REAL8 Alpha;
  REAL8 Delta;
  REAL8 AlphaBand;
  REAL8 DeltaBand;
  REAL8 dAlpha;
  REAL8 dDelta;

  /* these are ONLY used for manual sky-stepping: */
  INT8 AlphaCounter;  /* loop counters for manual stepping */
  INT8 DeltaCounter;

  SkyPositionVector skyRegion; 

  /* these are used for the metric sky-grid: */
  INT2 useMetric;		/* MANUAL, PTOLE_METRIC or COHERENT_METRIC */
  INT2 internalOrder;		/* coord-order used internally in mesh-function: ORDER_ALPHA_DELTA or ORDER_DELTA_ALPHA */

  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   * !! NOTE: independently if the mesh+metric internals
   *  the grid HAS to be in "standard order", which we fix as  
   *  x = alpha, y = delta
   *!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  TwoDMeshNode	*grid; 		/* head of linked list of nodes */  
  TwoDMeshNode  *gridNode;	/* pointer to current grid-node in grid */

  PtoleMetricIn ptoleMetricPar;	/* input parameters needed for PtoleMetric() */
  MetricParamStruc coherentMetricPar; /* input params for CoherentMetric() */
  PulsarTimesParamStruc baryParams; /* used by CoherentMetric() */
  REAL8Vector *lambda;		/* doppler params for CoherentMetric() */

  TwoDMeshParamStruc meshpar;	/* input params for the 2D mesh */

} DopplerScanState;

/* the meshing+metric functions might use a different grid-coordinate order
 * than our "standard order", which we fixed as ORDER_ALPHA_DELTA
 */
enum {
  ORDER_ALPHA_DELTA,
  ORDER_DELTA_ALPHA
};

typedef enum
{
  LAL_METRIC_NONE = 0,
  LAL_METRIC_PTOLE,
  LAL_METRIC_COHERENT,
  LAL_METRIC_LAST
} LALMetricType;
  
/********************************************************** <lalLaTeX>
\vfill{\footnotesize\input{DopplerScanHV}}
\newpage\input{DopplerScanC}
******************************************************* </lalLaTeX> */

/* Function prototypes */

void InitDopplerScan( LALStatus *stat, DopplerScanState *scan, DopplerScanInit init);
void NextDopplerPos ( LALStatus *stat, DopplerPosition *pos, DopplerScanState *scan);
void FreeDopplerScan (LALStatus *stat, DopplerScanState *scan);

void ParseSkyRegion (LALStatus *stat, SkyPositionVector *skylist, const CHAR *input);

void LALMetricWrapper (LALStatus *stat, REAL8Vector *metric, PtoleMetricIn *input, LALMetricType type);

/********************************************************** <lalLaTeX>
\newpage\input{LALSampleTestC}
******************************************************* </lalLaTeX> */

#ifdef  __cplusplus
}
#endif  /* C++ protection. */
#endif  /* Double-include protection. */
