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

#include <lal/LALRCSID.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/LALDatatypes.h>
#include <lal/DataBuffer.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/TwoInterfFindChirp.h>
#include <lal/LALInspiralBank.h>

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
#define COHERENTINSPIRALH_ENDET 12
#define COHERENTINSPIRALH_EGDET 13
#define COHERENTINSPIRALH_EZDET 14
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
#define COHERENTINSPIRALH_MSGENDET "Number of detectors is 1; it should be greater than 1 and less than 4"
#define COHERENTINSPIRALH_MSGEGDET "Number of detectors is > 4; it should be greater than 1 and less than 4"
#define COHERENTINSPIRALH_MSGEZDET "Number of detectors is 0; it should be greater than 1 and less than 4"
/* </lalErrTable> */


#if 0
<lalLaTeX>
\subsection*{Types}
</lalLaTeX>
#endif
/* structure for describing a binary insipral event */
/* <lalVerbatim file="CoherentInspiralHCoherentInspiralEvent"> */
typedef struct
tagInspiralEventVector
{
  UINT4                                  numDetectors;
  InspiralEvent                         *event; /* ordered list of events*/
}
InspiralEventVector;

typedef struct
tagCoherentInspiralEvent
{
  CHAR                                   ifos[LIGOMETA_IFOS_MAX];
  UINT4                                  eventId;
  UINT4                                  timeIndex;
  REAL4                                  mass1;
  REAL4                                  mass2;
  REAL4                                  cohSNR;
  REAL4                                  theta;
  REAL4                                  phi;
  LIGOTimeGPS                            time;
  InspiralEventVector                   *inspEventVec;
  struct tagCoherentInspiralEvent       *next;
}
CoherentInspiralEvent;
/* </lalVerbatim> */
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{CoherentInspiralEvent}}
\idx[Type]{CoherentInspiralEvent}

%\input{CoherentInspiralHCoherentInspiralEvent}

\noindent This structure describes inspiral events in the data of a pair 
of detectors found by \texttt{CoherentInspiralEventFilter}.
The fields are:

\begin{description}
\item[\texttt{UINT4 eventId}] A unique number assigned by the filter 
routine to each network event it finds.

\item[\texttt{UINT4 timeIndex}] The index in the fiducial detector at 
which the event occured in the array containing the filter output.

\item[\texttt{REAL4 cohSNR}] The value of network $\rho^2$ for the event.

\item[\texttt{InspiralEvent *inspEventVec}] A pointer to a structure of type 
\texttt{InspiralEvent} to allow the construction of a linked list of events
in participating detectors.

\item[\texttt{struct tagCoherentInspiralEvent *next}] A pointer to a 
structure of type \texttt{CoherentInspiralEvent} to allow the construction of 
a linked list of network events.
\end{description}
</lalLaTeX>
#endif

/* --- parameter structure for the coherent inspiral filtering function ---- */
/* <lalVerbatim file="CoherentInspiralHCoherentInspiralFilterParams"> */
typedef struct
tagCoherentInspiralInitParams
{
  UINT4                         numDetectors;
  UINT4                         numSegments;
  UINT4                         numPoints;
  UINT4                         numBeamPoints;
  BOOLEAN                       cohSNROut;
}
CoherentInspiralInitParams;
/* </lalVerbatim> */

typedef struct
tagDetectorVector
{
  UINT4                   numDetectors;
  LALDetector            *detector;
}
DetectorVector;
/* </lalVerbatim> */


typedef struct
tagDetectorBeamArray
{
  UINT4                    numBeamPoints;
  REAL4TimeSeries         *thetaPhiVs;/* 4D array: theta,phi,v+,v- */ 
}
DetectorBeamArray;


typedef struct
tagCoherentInspiralBeamVector
{
  UINT4                   numDetectors;
  DetectorBeamArray      *detBeamArray;
}
CoherentInspiralBeamVector;


typedef struct
tagCoherentInspiralFilterParams
{
  INT4                          numTmplts;
  UINT4                         maximiseOverChirp;
  UINT4                         numDetectors;
  UINT4                         numSegments;
  UINT4                         numPoints;
  UINT4                         numBeamPoints;
  REAL4                         fLow;
  REAL4                         deltaT;
  REAL4                         cohSNRThresh;
  BOOLEAN                       cohSNROut;
  UINT2Vector                  *detIDVec; /* Note: H1, H2 are from same site, but are different detectors */
  DetectorVector               *detectorVec; /*stores detectors' site info */
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
\texttt{CoherentInspiralFilter()} function.

\begin{description}
\item[\texttt{UINT4 numPoints}] Number of time-points in the $z$ series
from each detector. This determines the number of time-points in the
\texttt{cohSNRVec} time-series.

\item[\texttt{REAL4 cohSNRThresh}] The value to threshold the multi-detector
coherent signal to noise
ratio square, $\rho^2$, on. If the signal to noise exceeds this value, then a
candidate event is generated. Must be $\ge 0$ on entry.

\item[\texttt{DetectorVector   detectors}] This structure is defined below. It specifies the detectors on which
the coherent search is being performed.

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
  UINT4                   numDetectors;
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
  InspiralTemplate            *tmplt;
  CoherentInspiralBeamVector  *beamVec;
  CoherentInspiralZVector     *multiZData;
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
 * function prototypes for memory management functions
 .*
 */

void
LALCoherentInspiralFilterInputInit (
    LALStatus                       *status,
    CoherentInspiralFilterInput    **input,
    CoherentInspiralInitParams      *params
    );

void
LALCoherentInspiralFilterInputFinalize (
    LALStatus                       *status,
    CoherentInspiralFilterInput    **input
    );

void
LALCoherentInspiralFilterParamsInit (
    LALStatus                       *status,
    CoherentInspiralFilterParams   **output,
    CoherentInspiralInitParams      *params
    );

void
LALCoherentInspiralFilterParamsFinalize (
    LALStatus                       *status,
    CoherentInspiralFilterParams   **output
    );

/*
 *
 * function prototypes for coherent inspiral filter function
 *
 */


#if 0
<lalLaTeX>
\newpage\input{CoherentInspiralHV}
</lalLaTeX>
#endif

void
LALCoherentInspiralFilterSegment (
    LALStatus                             *status,
    CoherentInspiralEvent                **eventList,
    CoherentInspiralFilterInput           *input,
    CoherentInspiralFilterParams          *params
    );



#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _COHERENTINSPIRALH_H */
