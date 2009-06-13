/*
*  Copyright (C) 2007 Sukanta Bose, Jolien Creighton, Sean Seader, Thomas Cokelaer
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
#define COHERENTINSPIRALH_EZDET 12
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
#define COHERENTINSPIRALH_MSGEZDET "Number of detectors should be greater than 1 and less than 5"
/* </lalErrTable> */

#if 0
<lalLaTeX>
\subsection*{Types}
</lalLaTeX>
#endif

/* --- structure for describing a binary insipral event ------------------ */
/* <lalVerbatim file="FindChirpHInspiralEvent"> */


/* </lalVerbatim> */
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{InspiralEvent}}
\idx[Type]{InspiralEvent}

\input{FindChirpHInspiralEvent}

\noindent This structure describes inspiral events found by \texttt{findchirp}.
The fields are:

\begin{description}
\item[\texttt{UINT4 id}] A unique number assigned by the filter routine to
each event it finds.

\item[\texttt{UINT4 segmentNumber}] The id number of the
\texttt{FindChirpDataSegment} in which the event was found.

\item[\texttt{LIGOTimeGPS time}] The GPS time at which the event occoured.

\item[\texttt{UINT4 timeIndex}] The index at which the event occoured in
the array containing the filter output.

\item[\texttt{InspiralTemplate tmplt}] The parameters of the inspiral template
for the event.

\item[\texttt{REAL4 snrsq}] The value of $\rho^2$ for the event.

\item[\texttt{REAL4 chisq}] The value of the $\chi^2$ veto for the event, if
it has been computed.

\item[\texttt{REAL4 sigma}] The value of the normalisation constant $\sigma$
for the event.

\item[\texttt{REAL4 effDist}] The effective distance in megaparsecs to the
event.

 \item[\texttt{REAL4 coaPhase}] The coalescence phase of the chirp.

\item[\texttt{CHAR ifoName[2]}] Array for storing the two character
interferometer name (e.g. L1, H2, etc.)

\item[\texttt{struct tagInspiralEvent *next}] A pointer to a structure of type
\texttt{InspiralEvent} to allow the construction of a linked list of events.
\end{description}
</lalLaTeX>
#endif

#if 0
<lalLaTeX>
\subsection*{Types}
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
  UINT4                         cohSNROut;
  UINT4                         cohH1H2SNROut;
  UINT4                         nullStatH1H2Out;
  UINT4                         nullStatOut;
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
  UINT4                         maximizeOverChirp;
  UINT4                         numDetectors;
  UINT4                         numSegments;
  INT4                          numPoints;
  UINT4                         numBeamPoints;
  REAL4                         fLow;
  REAL8                         deltaT;
  REAL4                         cohSNRThresh;
  REAL8Vector                  *sigmasqVec;
  REAL4                         templateNorm;
  INT4                          segmentLength; /* time points */
  UINT4                         cohSNROut;
  UINT4                         cohH1H2SNROut;
  UINT4                         nullStatH1H2Out;
  UINT4                         nullStatOut;
  UINT2Vector                  *detIDVec; /* Note: H1, H2 are from same site, but are different detectors */
  DetectorVector               *detectorVec; /*stores detectors' site info */
  REAL4TimeSeries              *cohSNRVec;
  REAL4TimeSeries              *cohH1H2SNRVec;
  REAL4TimeSeries              *nullStatH1H2Vec;
  REAL4TimeSeries              *nullStatVec;
  REAL4                         chirpTime;
  double                        decStep;
  double                        raStep;
}
CoherentInspiralFilterParams;
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{CoherentInspiralFilterParams}}
\idx[Type]{CoherentInspiralFilterParams}

\input{CoherentInspiralHCoherentInspiralFilterParams}

\noindent This structure provides the parameters used by the
\texttt{CoherentInspiralFilter()} function.

\begin{description}
\item[\texttt{UINT4 numPoints}] Number of time-points in the $c$ series
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
/* <lalVerbatim file="CoherentInspiralHCoherentInspiralCVector"> */
typedef struct
tagCoherentInspiralCVector
{
  UINT4                   numDetectors;
  COMPLEX8TimeSeries     *cData[4];
}
CoherentInspiralCVector;
/* </lalVerbatim> */

#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{CoherentInspiralCVector}}
\idx[Type]{CoherentInspiralCVector}

\input{CoherentInspiralHCoherentInspiralCVector}

\noindent This structure groups the $c = x+iy$ outputs of $M$ detectors
into an ordered set. The FindChirpFilter code, when separately run on the
data from multiple detectors, outputs a \texttt{COMPLEX8TimeSeries}, $c$, for
each detector. If a coherent search is to be performed on the data from
these $M$ detectors, one of the inputs required is the
\texttt{CoherentInspiralCVector} structure with a default vector
\texttt{length} of $M=6$ and with the vector index ordered as 0=H1, 1=L1,
2=V (Virgo), 3=G (GEO), 4=T (Tama), (just like the lalcached detector siteIDs)
and 5=H2. If a coherent search is to be performed on, say, the data from
H1, L1, Virgo, and GEO, then the \texttt{length}
member above will be set to 6 (by default), but the pointers to the fourth and
fifth \texttt{COMPLEX8TimeSeries} will be set to NULL; the remainder will
point to the $c$ outputs from the above 4 detectors, in that order.

\begin{description}
\item[\texttt{UINT4  length}] Length of the vector; set to 6 (by default)
for the total number of operating (or nearly so) interferometers.

\item[\texttt{COMPLEX8TimeSeries  *cData}] Pointer to the c outputs of
the 6 interferometers.
\end{description}
</lalLaTeX>
#endif


/* <lalVerbatim file="CoherentInspiralHCoherentInspiralFilterInput"> */
typedef struct
tagCoherentInspiralFilterInput
{
  InspiralTemplate            *tmplt;
  CoherentInspiralCVector     *multiCData;
  CoherentInspiralBeamVector  *beamVec;
}
CoherentInspiralFilterInput;
/* </lalVerbatim> */
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{CoherentInspiralFilterInput}}
\idx[Type]{CoherentInspiralFilterInput}

\input{CoherentInspiralFilterHCoherentInspiralFilterInput}

\noindent This structure provides the essential information for
computing the coherent SNR from the $c$ outputs of multiple detectors.
In addition to this, the code requires the beam-pattern coefficients
for the different detectors. These coefficients are currently
computed by a Mathematica code and are read in as ascii files directly
by the coherent code. But there are plans for the future where a new member
will be added to this structure to store these coefficients.

\begin{description}
\item[\texttt{CoherentInspiralCVector   *multiCData}] Pointer to the vector
of COMPLEX8TimeSeries, namely, \texttt{CoherentInspiralCVector}.
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
LALCoherentInspiralEstimatePsiEpsilonCoaPhase (
    LALStatus                             *status,
    INT4                                   caseID[LAL_NUM_IFO],
    REAL8                                 *sigmasq,
    REAL4                                  theta,
    REAL4                                  phi,
    COMPLEX8                               cData[4],
    REAL4                                 *inclination,
    REAL4                                 *polarization,
    REAL4                                 *coaPhase
    );

void
LALCoherentInspiralEstimateDistance (
    LALStatus                             *status,
    REAL8                                 *sigmasq,
    REAL4                                  templateNorm,
    REAL8                                  deltaT,
    INT4                                   segmentLength,  /* time pts */
    REAL4                                  coherentSNR,
    REAL4                                 *distance
    );

double XLALCoherentCBCParamEstim( double *psi_est, double *iota_est, double *coa_phase_est, double a1, double a2, double a3, double a4, double *eff_distance0,double *eff_distance1,double *eff_distance2,double *eff_distanceH1H2,double amplitudeConst, double chirpTime, double C_Real0, double C_Real1, double C_Real2,double C_Im0, double C_Im1, double C_Im2, REAL8 *sigmasq);

void
LALCoherentInspiralFilterSegment (
    LALStatus                             *status,
    MultiInspiralTable                    **eventList,
    CoherentInspiralFilterInput           *input,
    CoherentInspiralFilterParams          *params
    );



#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _COHERENTINSPIRALH_H */
