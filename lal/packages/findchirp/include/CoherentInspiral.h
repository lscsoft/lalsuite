/*----------------------------------------------------------------------- 
 * 
 * File Name: CoherentInspiral.h
 *
 * Author: Bose, S., Seader, S. E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="CoherentInspiralHV">
Author: Bose, S., Seader, S. E.
$Id$
</lalVerbatim> 

<lalLaTeX>
\section{Header \texttt{CoherentInspiral.h}}
\label{s:CoherentInspiral.h}

\noindent Provides core prototypes, structures and functions to filter
data from multiple interferometers coherently for binary inspiral chirps.  

\subsection*{Coherent search statistic for binary neutron stars}

The coherent statistic will be defined here.

\subsubsection*{Synopsis}

\begin{verbatim}
#include <lal/CoherentInspiral.h>
\end{verbatim}

\input{CoherentInspiralHDoc}

</lalLaTeX>
#endif

#ifndef _COHERENTINSPIRALH_H
#define _COHERENTINSPIRALH_H

#include <lal/LALDatatypes.h>
#include <lal/ComplexFFT.h>
#include <lal/DataBuffer.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LALInspiral.h>
#include <lal/LALInspiralBank.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/FindChirpChisq.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif


NRCSID (COHERENTINSPIRALH, "$Id$");

#if 0
<lalLaTeX> 
\subsection*{Error codes} 
</lalLaTeX>
#endif
/* <lalErrTable> */
#define COHERENTINSPIRALH_ENULL 1
#define COHERENTINSPIRALH_ENNUL 2
#define COHERENTINSPIRALH_EALOC 3
#define COHERENTINSPIRALH_ENUMZ 4
#define COHERENTINSPIRALH_ESEGZ 5
#define COHERENTINSPIRALH_ECHIZ 6
#define COHERENTINSPIRALH_EDTZO 7
#define COHERENTINSPIRALH_EFREE 8
#define COHERENTINSPIRALH_ERHOT 9
#define COHERENTINSPIRALH_ECHIT 10
#define COHERENTINSPIRALH_ESMSM 11
#define COHERENTINSPIRALH_MSGENULL "Null pointer"
#define COHERENTINSPIRALH_MSGENNUL "Non-null pointer"
#define COHERENTINSPIRALH_MSGEALOC "Memory allocation error"
#define COHERENTINSPIRALH_MSGENUMZ "Invalid number of points in segment"
#define COHERENTINSPIRALH_MSGESEGZ "Invalid number of segments"
#define COHERENTINSPIRALH_MSGECHIZ "Invalid number of chi squared bins"
#define COHERENTINSPIRALH_MSGEDTZO "deltaT is zero or negative"
#define COHERENTINSPIRALH_MSGEFREE "Error freeing memory"
#define COHERENTINSPIRALH_MSGERHOT "coherentSNR threshold is negative"
#define COHERENTINSPIRALH_MSGECHIT "Chisq threshold is negative"
#define COHERENTINSPIRALH_MSGESMSM "Size mismatch between vectors"
/* </lalErrTable> */


#if 0
<lalLaTeX>
\subsection*{Types}
</lalLaTeX>
#endif

/* --- parameter structure for the coherent inspiral filtering function ---- */
/* <lalVerbatim file="CoherentInspiralHCoherentInspiralFilterParams"> */
typedef struct
tagCoherentInspiralFilterParams
{
  UINT4                         numPoints;
  REAL4                         cohSNRThresh;
  REAL4TimeSeries              *cohSNRVec;
}
CoherentInspiralFilterParams;
/* </lalVerbatim> */
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{CoherentInspiralFilterParams}}
\idx[Type]{CoherentInspiralFilterParams}

\input{CoherentInspiralHCoherentInspiralFilterParams}

\noindent This structure provides the parameters used by the
\texttt{CoherentInspiralFilterSegment()} function.

\begin{description}
\item[\texttt{UINT4 numPoints}] Number of time-points in the $z$ series
from each detector. This determines the number of time-points in the
\texttt{cohSNRVec} time-series.

\item[\texttt{REAL4 cohSNRThresh}] The value to threshold the multi-detector
coherent signal to noise
ratio square, $\rho^2$, on. If the signal to noise exceeds this value, then a
candidate event is generated. Must be $\ge 0$ on entry.

\item[\texttt{REAL4Vector *cohSNRVec}] Pointer to a vector that is set to
$\rho^2(t_j)$ on exit. If NULL $\rho^2(t_j)$ is not stored.

\end{description}
</lalLaTeX>
#endif


/* --- input to the CoherentInspiral filtering functions --------- */
/* <lalVerbatim file="CoherentInspiralHCoherentInspiralZVector"> */
typedef struct
tagCoherentInspiralZVector
{
  UINT4                   length;
  COMPLEX8TimeSeries     *zData;
}
CoherentInspiralZVector;
/* </lalVerbatim> */
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{CoherentInspiralZVector}}
\idx[Type]{CoherentInspiralZVector}

\input{CoherentInspiralHCoherentInspiralZVector}

\noindent This structure groups the $z = x+iy$ outputs of $M$ detectors 
into an ordered set. The FindChirpFilter code, when separately run on the 
data from multiple detectors, outputs a \texttt{COMPLEX8TimeSeries}, $z$, for
each detector. If a coherent search is to be performed on the data from
these $M$ detectors, one of the inputs required is the 
\texttt{CoherentInspiralZVector} structure with a default vector 
\texttt{length} of $M=6$ and with the vector index ordered as 0=H1, 1=L1, 
2=V (Virgo), 3=G (GEO), 4=T (Tama), (just like the lalcached detector siteIDs) 
and 5=H2. If a coherent search is to be performed on, say, the data from 
H1, L1, Virgo, and GEO, then the \texttt{length} 
member above will be set to 6 (by default), but the pointers to the fourth and 
fifth \texttt{COMPLEX8TimeSeries} will be set to NULL; the remainder will 
point to the $z$ outputs from the above 4 detectors, in that order.

\begin{description}
\item[\texttt{UINT4  length}] Length of the vector; set to 6 (by default) 
for the total number of operating (or nearly so) interferometers.

\item[\texttt{COMPLEX8TimeSeries  *zData}] Pointer to the z outputs of 
the 6 interferometers.
\end{description}
</lalLaTeX>
#endif


/* <lalVerbatim file="CoherentInspiralHCoherentInspiralFilterInput"> */
typedef struct
tagCoherentInspiralFilterInput
{
  CoherentInspiralZVector   *multiZData;
}
CoherentInspiralFilterInput;
/* </lalVerbatim> */
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{CoherentInspiralFilterInput}}
\idx[Type]{CoherentInspiralFilterInput}

\input{CoherentInspiralFilterHCoherentInspiralFilterInput}

\noindent This structure provides the essential information for
computing the coherent SNR from the $z$ outputs of multiple detectors.
In addition to this, the code requires the beam-pattern coefficients 
for the different detectors. These coefficients are currently 
computed by a Mathematica code and are read in as ascii files directly
by the coherent code. But there are plans for the future where a new member
will be added to this structure to store these coefficients.

\begin{description}
\item[\texttt{CoherentInspiralZVector   *multiZData}] Pointer to the vector
of COMPLEX8TimeSeries, namely, \texttt{CoherentInspiralZVector}.
\end{description}
</lalLaTeX>
#endif

#if 0
<lalLaTeX>
\vfill{\footnotesize\input{CoherentInspiralHV}}
</lalLaTeX> 
#endif


/*
 *
 * function prototypes for coherent inspiral filter function
 *
 */


#if 0
<lalLaTeX>
\newpage\input{CoherentInspiralFilterC}
</lalLaTeX>
#endif

void
LALCoherentInspiralFilterSegment (
    LALStatus                             *status,
    MultiInspiralTable                   **eventList,
    CoherentInspiralFilterInput           *input,
    CoherentInspiralFilterParams          *params
    );



#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _COHERENTINSPIRALH_H */
