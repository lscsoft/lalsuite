/**
 * \author Jones, D. I.  Owen, B. J.
 * \date 2001 -- 2006
 * \file 
 * \ingroup PulsarMetric
 * \brief Provides routines to compute pulsar parameter-space metrics using the
 * ``Ptolemaic'' approximation.
 *
 * $Id$
 *
 */

#ifndef _PTOLEMETRIC_H
#define _PTOLEMETRIC_H

/** 
    This header covers routines for using a ``Ptolemaic'' (epicyclic)
    approximation to the detector motion to compute the parameter-space metric
    for a pulsar search. (At the moment, the search is assumed to be a single
    coherent integration.) The results should be very similar to those under the
    <tt>StackMetric.h</tt> header, and reading that documention is a good
    background for this documentation.

    Why this extra header? Two words: simplicity and speed. The metric
    components can be expressed analytically in terms of trig functions,
    allowing one to get a feel for what the parameter space will look like
    before using a single CPU cycle. In addition, CPU usage is much reduced
    (compared to the routines in \texttt{StackMetric.h}) in numerical
    explorations such as testing the suitability of various tiling codes. Thus,
    the functions in this header can be very useful in the current stage of
    exploring parameter space and wondering how we can practically take
    advantage of correlations. It's also good at catching bugs and errors in the
    numerical routines under <tt>StackMetric.h</tt>. The effectiveness of the
    tiling at catching signals should be very little reduced by the
    approximation. Jones, Owen, and Whitbeck will write a short paper on this
    and other details.
*/

#include <lal/DetectorSite.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/SkyCoordinates.h>
#include <lal/LALBarycenter.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID( PTOLEMETRICH, "$Id$" );

/** @{ \name Error conditions
 */
#define PTOLEMETRICH_ENULL   1
#define PTOLEMETRICH_EPARM   2
#define PTOLEMETRICH_EDIM    3
#define PTOLEMETRICH_ENONULL 4
#define PTOLEMETRICH_EMETRIC 5

#define PTOLEMETRICH_MSGENULL   "unexpected null pointer"
#define PTOLEMETRICH_MSGEPARM   "bad parameter value"
#define PTOLEMETRICH_MSGEDIM    "bad array length"
#define PTOLEMETRICH_MSGENONULL "unexpected non-null pointer"
/* @} */

/** Constants defining different types of pulsar-metrics. */
typedef enum
{
  LAL_PMETRIC_NONE = 0,
  LAL_PMETRIC_COH_PTOLE_ANALYTIC,
  LAL_PMETRIC_COH_PTOLE_NUMERIC,
  LAL_PMETRIC_COH_EPHEM,
  LAL_PMETRIC_LAST
} LALPulsarMetricType;
/* </lalVerbatim> */

#define PMETRIC_MIN(x,y) ((x) < (y) ? (x) : (y))
#define PMETRIC_MAX(x,y) ((x) > (y) ? (x) : (y))

/** Translate metrix matrix-indices (a,b) into vector-index l */
#define PMETRIC_INDEX(a,b) (PMETRIC_MIN((a),(b))+PMETRIC_MAX((a),(b))*(PMETRIC_MAX((a),(b)) + 1 ) / 2 )


/** This structure will likely be changed to match up better with
    those under the StackMetric.h header. It contains the bare
    necessities, not needing function pointers etc. */
typedef struct
tagPtoleMetricIn
{
  SkyPosition    position;	/**< The equatorial coordinates at which the metric components are evaluated. */
  REAL4Vector   *spindown;	/**< The (dimensionless) spindown parameters for which the metric components are evaluated. */
  LIGOTimeGPS    epoch;		/**< When the coherent integration begins */
  REAL4          duration;	/**< Duration of integration, in seconds. */
  REAL4          maxFreq;	/**< The maximum frequency to be searched, in Hz. */
  const LALDetector    *site;	/**< The detector site, used only for its latitude and longitude. */
  EphemerisData  *ephemeris;	/**< Not used for the Ptolemaic approximation, this is for compatibility with other metrics. */
  LALPulsarMetricType metricType; /**< The type of metric to use: analytic, Ptolemaic or fully ephemeris-based. */
} PtoleMetricIn;


/* ----- prototypes ----- */
void
LALPtoleMetric( LALStatus      *status,
                REAL8Vector    *metric,
                PtoleMetricIn  *input );

void
LALPulsarMetric( LALStatus      *status,
                 REAL8Vector    **metric,
                 PtoleMetricIn  *input );

int XLALFindMetricDim ( const REAL8Vector *metric );


#ifdef  __cplusplus
}
#endif

#endif
