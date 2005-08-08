/**************************************** <lalVerbatim file="PtoleMetricHV">
Author: Jones, D. I.    Owen, B. J.
$Id$
********************************************************** </lalVerbatim> */

/**************************************************************** <lalLaTeX>
 
\section{Header \texttt{PtoleMetric.h}}
\label{s:PtoleMetric.h}
 
Provides routines to compute pulsar parameter-space metrics using the
``Ptolemaic'' approximation.
 
\subsection*{Synopsis}
\begin{verbatim}
#include <lal/PtoleMetric.h>
\end{verbatim}
 
\noindent This header covers routines for using a ``Ptolemaic'' (epicyclic)
approximation to the detector motion to compute the parameter-space metric
for a pulsar search. (At the moment, the search is assumed to be a single
coherent integration.) The results should be very similar to those under the
\texttt{StackMetric.h} header, and reading that documention is a good
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
numerical routines under \texttt{StackMetric.h}.  The effectiveness of the
tiling at catching signals should be very little reduced by the
approximation. Owen will write a short paper on this and other details.

**************************************************************</lalLaTeX> */

#ifndef _PTOLEMETRIC_H
#define _PTOLEMETRIC_H

#include <lal/DetectorSite.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/SkyCoordinates.h>
#include <lal/LALBarycenter.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID( PTOLEMETRICH, "$Id$" );

/**************************************************************** <lalLaTeX>
\subsection*{Error conditions}
*********************************************** </lalLaTeX> <lalErrTable> */

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

/************************************************* </lalErrTable> <lalLaTeX>
\subsection*{Structures}

\begin{verbatim}
PtoleMetricIn
\end{verbatim}
\idx[Type]{PtoleMetricIn}

\noindent This structure will likely be changed to match up better with
those under the \texttt{StackMetric.h} header. It contains the bare
necessities, not needing function pointers etc.

\begin{description}

\item[\texttt{SkyPosition position}] The equatorial coordinates at which the
metric components are evaluated.

\item[\texttt{REAL4Vector *spindown}] The spindown parameters for which the
metric components are evaluated. These are dimensionless.

\item[\texttt{LIGOTimeGPS epoch}] When the coherent integration begins.

\item[\texttt{REAL4 duration}] Duration of integration, in seconds.

\item[\texttt{REAL4 maxFreq}] The maximum frequency to be searched, in Hz.

\item[\texttt{LALDetector *site}] The detector site, used only for its
latitude and longitude.

\item[\texttt{EphemerisData *ephemeris}] Not used for the Ptolemaic
approximation, this is for compatibility with other metrics.

\item[\texttt{LALPulsarMetricType metricType}] the type of metric to use: analytic, Ptolemaic or fully ephemeris-based.

\end{description}

\subsection*{Constants}

\begin{verbatim}
enum LALPulsarMetricType
\end{verbatim}
\idx[Type]{LALPulsarMetricType}

\noindent Constants defining different types of pulsar-metrics.
</lalLaTeX> */
/* <lalVerbatim> */
typedef enum
{
  LAL_PMETRIC_NONE = 0,
  LAL_PMETRIC_COH_PTOLE_ANALYTIC,
  LAL_PMETRIC_COH_PTOLE_NUMERIC,
  LAL_PMETRIC_COH_EPHEM,
  LAL_PMETRIC_LAST
} LALPulsarMetricType;
/* </lalVerbatim> */

typedef struct
tagPtoleMetricIn
{
  SkyPosition    position;
  REAL4Vector   *spindown;
  LIGOTimeGPS    epoch;
  REAL4          duration;
  REAL4          maxFreq;
  LALDetector    *site;
  EphemerisData  *ephemeris;
  LALPulsarMetricType metricType;
}
PtoleMetricIn;

/**************************************************************** <lalLaTeX>
\vfill{\footnotesize\input{PtoleMetricHV}}
\newpage\input{PtoleMetricC}
************************************************************* </lalLaTeX> */

void
LALPtoleMetric( LALStatus      *status,
                REAL8Vector    *metric,
                PtoleMetricIn  *input );

void
LALPulsarMetric( LALStatus      *status,
                 REAL8Vector    **metric,
                 PtoleMetricIn  *input );

int XLALFindMetricDim ( const REAL8Vector *metric );

/**************************************************************** <lalLaTeX>
\newpage\input{PtoleMetricTestC}
\newpage\input{PtoleMeshTestC}
\newpage\input{GeneralMetricTestC}
\newpage\input{GeneralMeshTestC}
************************************************************* </lalLaTeX> */

#ifdef  __cplusplus
}
#endif

#endif
