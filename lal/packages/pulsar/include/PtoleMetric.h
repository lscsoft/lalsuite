/**************************************** <lalVerbatim file="PtoleMetricHV">
Author: Owen, B. J.
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

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID( PTOLEMETRICH, "$Id$" );

/**************************************************************** <lalLaTeX>
\subsection*{Error conditions}
*********************************************** </lalLaTeX> <lalErrTable> */

#define PTOLEMETRICH_ENULL 0x01
#define PTOLEMETRICH_EPARM 0x02
#define PTOLEMETRICH_EDIM  0x04

#define PTOLEMETRICH_MSGENULL "unexpected null pointer"
#define PTOLEMETRICH_MSGEPARM "bad parameter value"
#define PTOLEMETRICH_MSGEDIM  "bad array length"

/************************************************* </lalErrTable> <lalLaTeX>
\subsection*{Structures}

\begin{verbatim}
PtoleMetricIn
\end{verbatim}
\index{\texttt{PtoleMetricIn}}

\noindent This structure will likely be changed to match up better with
those under the \texttt{StackMetric.h} header. It contains the bare
necessities, not needing function pointers etc.

\begin{description}

\item[\texttt{SkyPosition position}] The position on the sky at which the
metric components are evaluated.

\item[\texttt{REAL4Vector *spindown}] The spindown parameters for which the
metric components are evaluated. These are dimensionless.

\item[\texttt{LIGOTimeGPS epoch}] When the coherent integration begins.

\item[\texttt{REAL4 duration}] Duration of integration, in seconds.

\item[\texttt{REAL4 maxFreq}] The maximum frequency to be searched, in Hz.

\item[\texttt{LALFrDetector site}] The detector site, used only for its
latitude and longitude.

\end{description}

************************************************************* </lalLaTeX> */

typedef struct
tagPtoleMetricIn
{
  SkyPosition    position;
  REAL4Vector   *spindown;
  LIGOTimeGPS    epoch;
  REAL4          duration;
  REAL4          maxFreq;
  LALFrDetector  site;
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

/**************************************************************** <lalLaTeX>
\newpage\input{PtoleMetricTestC}
************************************************************* </lalLaTeX> */

#ifdef  __cplusplus
}
#endif

#endif
