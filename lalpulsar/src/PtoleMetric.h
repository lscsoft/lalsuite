/*
*  Copyright (C) 2007 Thomas Essinger-Hileman, Jolien Creighton, Ian Jones, Benjamin Owen, Reinhard Prix, Karl Wette
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

#ifndef _PTOLEMETRIC_H
#define _PTOLEMETRIC_H

#ifdef  __cplusplus
extern "C" {
#endif

#include <gsl/gsl_matrix.h>
#include <lal/DetectorSite.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/SkyCoordinates.h>
#include <lal/LALBarycenter.h>

/**
 * \defgroup PtoleMetric_h Header PtoleMetric.h
 * \author Jones, D. I.  Owen, B. J.
 * \date 2001 -- 2006
 * \ingroup pkg_pulsarMetric
 * \brief Provides routines to compute pulsar parameter-space metrics using the
 * ``Ptolemaic'' approximation.
 *
 * \heading{Synopsis}
 \code
 #include <lal/PtoleMetric.h>
 \endcode

 This module covers routines for using a ``Ptolemaic'' (epicyclic)
 approximation to the detector motion to compute the parameter-space metric
 for a pulsar search. (At the moment, the search is assumed to be a single
 coherent integration.) The results should be very similar to those under
 \ref StackMetric_h, and reading that documention is a good
 background for this documentation.

 Why this extra module? Two words: simplicity and speed. The metric
 components can be expressed analytically in terms of trig functions,
 allowing one to get a feel for what the parameter space will look like
 before using a single CPU cycle. In addition, CPU usage is much reduced
 (compared to the routines in \ref StackMetric_h) in numerical
 explorations such as testing the suitability of various tiling codes. Thus,
 the functions in this header can be very useful in the current stage of
 exploring parameter space and wondering how we can practically take
 advantage of correlations. It's also good at catching bugs and errors in the
 numerical routines under \ref StackMetric_h. The effectiveness of the
 tiling at catching signals should be very little reduced by the
 approximation. Jones, Owen, and Whitbeck will write a short paper on this
 and other details.

 *
 */
/*@{*/

#define PMETRIC_MIN(x,y) ((x) < (y) ? (x) : (y))
#define PMETRIC_MAX(x,y) ((x) > (y) ? (x) : (y))

/** Translate metrix matrix-indices (a,b) into vector-index l */
#define PMETRIC_INDEX(a,b) (PMETRIC_MIN((a),(b))+PMETRIC_MAX((a),(b))*(PMETRIC_MAX((a),(b)) + 1 ) / 2 )

/** \name Error conditions */
/*@{*/
#define PTOLEMETRICH_ENULL   1
#define PTOLEMETRICH_EPARM   2
#define PTOLEMETRICH_EDIM    3
#define PTOLEMETRICH_ENONULL 4
#define PTOLEMETRICH_EMETRIC 5

#define PTOLEMETRICH_MSGENULL   "unexpected null pointer"
#define PTOLEMETRICH_MSGEPARM   "bad parameter value"
#define PTOLEMETRICH_MSGEDIM    "bad array length"
#define PTOLEMETRICH_MSGENONULL "unexpected non-null pointer"
#define PTOLEMETRICH_MSGEMETRIC "unknown metric type"
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


/** This structure will likely be changed to match up better with
    those under \ref StackMetric_h; it contains the bare
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
  const EphemerisData  *ephemeris;	/**< Not used for the Ptolemaic approximation, this is for compatibility with other metrics. */
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

int XLALSpindownMetric(
  gsl_matrix* metric,
  double Tspan
  );

/*@}*/

#ifdef  __cplusplus
}
#endif

#endif
