/*---------------------------------------------------------------------------
 *
 * File Name: TwoInterfFindChirp.h
 *
 * Author: Bose, S., Sorensen, S., and Noel, J. S.
 *
 * Revision: $Id$
 *
 *---------------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="TwoInterfFindChirpHV">
Author: Bose, S. and Sorensen, S.
$Id$
</lalVerbatim> 

<lalLaTeX>
\section{Header \texttt{TwoInterfFindChirp.h}}
\label{s:TwoInterfFindChirp.h}

Provides core prototypes, structures and functions to filter data from
a pair of interferometers for binary inspiral chirps.

\subsection*{Synopsis}

\begin{verbatim}
#include <lal/TwoInterfFindChirp.h>
\end{verbatim}

%\input{TwoInterfFindChirpHDoc}
</lalLaTeX>
#endif

#ifndef _TWOINTERFFINDCHIRPH_H
#define _TWOINTERFFINDCHIRPH_H

#include <lal/LALDatatypes.h>
#include <lal/ComplexFFT.h>
#include <lal/DataBuffer.h>
#ifdef LAL_ENABLE_MPI
#include <lal/Comm.h>
#endif
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/FindChirpChisq.h>
#include <lal/DetectorSite.h>
#include <lal/StochasticCrossCorrelation.h>
#include <lal/SkyCoordinates.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif
  
NRCSID (TWOINTERFFINDCHIRPH, "$Id$");
  
#if 0
<lalLaTeX> 
\subsection*{Error codes} 
</lalLaTeX>
#endif
/***************** <lalErrTable> */
#define TWOINTERFFINDCHIRPH_ENULL 1
#define TWOINTERFFINDCHIRPH_ENNUL 2
#define TWOINTERFFINDCHIRPH_EALOC 3
#define TWOINTERFFINDCHIRPH_ENUMZ 5
#define TWOINTERFFINDCHIRPH_ESEGZ 6
#define TWOINTERFFINDCHIRPH_ECHIZ 7
#define TWOINTERFFINDCHIRPH_EDTZO 8
#define TWOINTERFFINDCHIRPH_EMLZO 9
#define TWOINTERFFINDCHIRPH_ETRNC 10
#define TWOINTERFFINDCHIRPH_EFLOW 11
#define TWOINTERFFINDCHIRPH_EFREE 12
#define TWOINTERFFINDCHIRPH_ENUMF 13
#define TWOINTERFFINDCHIRPH_ERHOT 15
#define TWOINTERFFINDCHIRPH_ECHIT 16
#define TWOINTERFFINDCHIRPH_ECRUP 17
#define TWOINTERFFINDCHIRPH_MSGENULL "Null pointer"
#define TWOINTERFFINDCHIRPH_MSGENNUL "Non-null pointer"
#define TWOINTERFFINDCHIRPH_MSGEALOC "Memory allocation error"
#define TWOINTERFFINDCHIRPH_MSGENUMZ "Invalid number of points in segment"
#define TWOINTERFFINDCHIRPH_MSGESEGZ "Invalid number of segments"
#define TWOINTERFFINDCHIRPH_MSGECHIZ "Invalid number of chi squared bins"
#define TWOINTERFFINDCHIRPH_MSGEDTZO "deltaT is zero or negative or unequal for detector 1 and 2 data segments"
#define TWOINTERFFINDCHIRPH_MSGEMLZO "maxLag is zero or negative"
#define TWOINTERFFINDCHIRPH_MSGETRNC "Duration of inverse spectrum in time domain is negative"
#define TWOINTERFFINDCHIRPH_MSGEFLOW "Inverse spectrum low frequency cutoff is negative"
#define TWOINTERFFINDCHIRPH_MSGEFREE "Memory free error"
#define TWOINTERFFINDCHIRPH_MSGENUMF "Invalid number of points in filter"
#define TWOINTERFFINDCHIRPH_MSGERHOT "Rhosq threshold is zero or negative"
#define TWOINTERFFINDCHIRPH_MSGECHIT "Chisq threshold is zero or negative"
#define TWOINTERFFINDCHIRPH_MSGECRUP "Chirp length or invSpecTrunc too long for length of data segment"
/* </lalErrTable> */

/*
 *
 * typedefs of structures used by the twointerffindchirp functions
 *
 */

#if 0
<lalLaTeX>
\subsection*{Types}
</lalLaTeX>
#endif

/* --- structure for describing a binary insipral event ------------------ */
/* <lalVerbatim file="FindChirpHInspiralEvent"> */
typedef struct
tagInspiralEvent
{
  UINT4                         id;
  UINT4                         segmentNumber;
  LIGOTimeGPS                   time;
  LIGOTimeGPS                   impulseTime;
  REAL8                         templateDuration;
  REAL8                         eventDuration;
  UINT4                         timeIndex;
  InspiralTemplate              tmplt;
  REAL4                         snrsq;
  REAL4                         chisq;
  REAL4                         sigma;
  REAL4                         effDist;
  REAL4                         coaPhase;
  UINT4                         numChisqBins;
  CHAR                          ifoName[2];
  CHAR                          channel[LALNameLength];
  WaveformGenerator             wavGen;
  struct tagInspiralEvent      *next;
}
InspiralEvent;
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

\item[\texttt{WaveformGenerator wavGen}] The particular waveform generator that was used in the matched filter for this chirp.

\item[\texttt{CHAR ifoName[2]}] Array for storing the two character
interferometer name (e.g. L1, H2, etc.)

\item[\texttt{struct tagInspiralEvent *next}] A pointer to a structure of type 
\texttt{InspiralEvent} to allow the construction of a linked list of events.
\end{description}
</lalLaTeX>
#endif

/* structure for describing a binary insipral event */
/* <lalVerbatim file="TwoInterfFindChirpHTwoInterfInspiralEvent"> */
typedef struct
tagTwoInterfInspiralEvent
{
  UINT4                                  twoInterfId;
  UINT4                                  segmentNumber;
  LIGOTimeGPS                            time;
  LIGOTimeGPS                            impulseTime;
  UINT4                                  timeIndex;
  UINT4                                  timeIndex2;
  InspiralTemplate                       tmplt;
  REAL4                                  snrsq;
  REAL4                                  sigma;
  REAL4                                  effDist;
  REAL4                                  chisq1;
  REAL4                                  chisq2;
  REAL4                                  twoInterfAxisRa; /*Direction of detector1-to-detector2 (e.g., Hanford-to-Livingston ray at time of eventthe central axis of the cone on which the source lies)*/
  REAL4                                  twoInterfAxisDec; 
  REAL4                                  twoInterfAngle; /*Wave arrival angle with respect to detector1-to-detector2 (e.g., Hanford-to-Livingston) ray...*/
  REAL4                                  twoInterfAngleSig; /*...and error*/
  InspiralEvent                         *eventIn1;
  InspiralEvent                         *eventIn2;
  struct tagTwoInterfInspiralEvent      *next;
}
TwoInterfInspiralEvent;
/* </lalVerbatim> */
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{TwoInterfInspiralEvent}}
\idx[Type]{TwoInterfInspiralEvent}

%\input{TwoInterfFindChirpHTwoInterfInspiralEvent}

\noindent This structure describes inspiral events in the data of a pair 
of detectors found by \texttt{TwoInterffindchirp}.
The fields are:

\begin{description}
\item[\texttt{UINT4 twoInterfId}] A unique number assigned by the filter 
routine to each network event it finds.

\item[\texttt{UINT4 segmentNumber}] The id number of the 
\texttt{FindChirpDataSegment} in the fiducial detector (chosen
as detector 1 here) in which the event was found.

\item[\texttt{LIGOTimeGPS time}] The GPS time in the fiducial detector at 
which the event occured.

\item[\texttt{UINT4 timeIndex}] The index in the fiducial detector at 
which the event occured in the array containing the filter output.

\item[\texttt{InspiralTemplate tmplt}] The parameters of the inspiral template
for the event.

\item[\texttt{REAL4 snrsq}] The value of network $\rho^2$ for the event.

\item[\texttt{InspiralEvent event1}] A pointer to a structure of type 
\texttt{InspiralEvent} to allow the construction of a linked list of events
in detector 1.

\item[\texttt{InspiralEvent event2}] Similar to the \texttt{event1} structure,
but for constructing a linked list of events
in detector 2.

\item[\texttt{struct tagTwoInterfInspiralEvent *next}] A pointer to a 
structure of type \texttt{InspiralEvent} to allow the construction of 
a linked list of network events.
\end{description}
</lalLaTeX>
#endif

#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{MultiInspiralTable}}
\idx[Type]{MultiInspiralTable}

%\input{TwoInterfFindChirpHMultiInspiralTable}

\noindent This is an alternative structure describing 
inspiral events in the data of a pair of detectors found by \texttt{TwoInterffindchirp}.
This structure is being redesigned to be eventually consistent with the database table
definitions.
</lalLaTeX>
#endif

/* --- vector of DataSegmentVector, as defined in the 
   findchirp/framedata packages --------- */
/* <lalVerbatim file="TwoInterfFindChirpHTwoInterfDataSegmentVector"> */
typedef struct
tagTwoInterfDataSegmentVector
{
  UINT4                         length;
  DataSegmentVector            *data;
}
TwoInterfDataSegmentVector;
/* </lalVerbatim> */
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{TwoInterfDataSegmentVector}}
\idx[Type]{TwoInterfDataSegmentVector}

%\input{TwoInterfFindChirpHTwoInterfDataSegmentVector}

\noindent This structure provides a LAL like vector structure for the
\texttt{DataSegmentVector} structure defined elsewhere
in the package \texttt{findchirp}.

\begin{description}
\item[\texttt{UINT4 length}] Number of \texttt{DataSegmentVector} 
structures in the vector. It is set equal to the number of detectors
in a network (in this case 2).

\item[\texttt{DataSegmentVector *data}] Pointer to the data.
\end{description}
</lalLaTeX>
#endif

/* vector of FindChirpSegmentVector defined in the finchirp package */
/* <lalVerbatim file="TwoInterfFindChirpSegmentVector"> */
typedef struct
tagTwoInterfFindChirpSegmentVector
{
  UINT4                         length;
  FindChirpSegmentVector       *data;
}
TwoInterfFindChirpSegmentVector;
/* </lalVerbatim> */
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{TwoInterfFindChirpSegmentVector}}
\idx[Type]{TwoInterfFindChirpSegmentVector}

%\input{TwoInterfFindChirpHTwoInterfFindChirpSegmentVector}

\noindent This structure provides a LAL like vector structure for the
\texttt{FindChirpSegmentVector} structure defined elsewhere in the
findchirp package.

\begin{description}
\item[\texttt{UINT4 length}] Number of \texttt{FindChirpSegmentVector} 
structures in the vector. It is set equal to the number of detectors
in a network (in this case 2).

\item[\texttt{FindChirpSegmentVector *data}] Pointer to the data.
\end{description}
</lalLaTeX>
#endif

/*
 *
 * typedefs of parameter structures used by functions in twointerffindchirp
 *
 */

/* parameter structure for all init functions */
/* <lalVerbatim file="FindChirpHFindChirpInitParams"> */
typedef struct
tagTwoInterfFindChirpInitParams
{
  UINT4                         numDetectors;
  UINT4                         numSegments;
  UINT4                         numPoints;  
  UINT4                         numChisqBins;
  BOOLEAN                       createRhosqVec;
  UINT4                         ovrlap;
  BOOLEAN                       createTwoInterfRhosqVec;
}
TwoInterfFindChirpInitParams;
/* </lalVerbatim> */
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{TwoInterfFindChirpInitParams}}
\idx[Type]{TwoInterfFindChirpInitParams}

%\input{TwoInterfFindChirpHTwoInterfFindChirpInitParams}

\noindent This structure provides the essential information for the
filter initialization and memory allocation functions used in the
\texttt{TwoInterfFindChirp} codes in the \texttt{FindChirp} package.

\begin{description}
\item[\texttt{UINT4 numSegments}] The number of detectors in the network
(which is set equal to 2 here)

\item[\texttt{UINT4 numSegments}] The number of data segments per
detector to allocate storage for.

\item[\texttt{UINT4 numPoints}] The number of discrete data points in each
data segment. 

\item[\texttt{UINT4 numChisqBins}] The number of bins used to contruct the
$\chi^2$ veto.

\item[\texttt{BOOLEAN createRhosqVec}] Debugging flag that controls whether
or not the function \texttt{TwoInterfFindChirpFilterSegment()} should store 
the output of the filter, $\rho_I^2(t)$, of each detector ($I = 1$, 2)
separately (apart from generating the network events). If the flag is set 
to 1, memory is only allocated for two filter vectors, one for each detector.

\item[\texttt{BOOLEAN createTwoInterfRhosqVec}] Debugging flag that controls 
whether or not the function \texttt{TwoInterfFindChirpFilterSegment()} should 
store the output of the network filter, $\rho^2(t)$. 
Memory is only allocated for this vector if the flag is set to 1.
\end{description}
</lalLaTeX>
#endif

/* parameter structure for the filtering function */
/* <lalVerbatim file="TwoInterfFindChirpHTwoInterfFindChirpFilterParamsVector"> */
typedef struct
tagTwoInterfFindChirpFilterParamsVector
{
  UINT4                         length;
  FindChirpFilterParams        *filterParams;
}
TwoInterfFindChirpFilterParamsVector;
/* </lalVerbatim> */
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{TwoInterfFindChirpFilterParams}}
\idx[Type]{TwoInterfFindChirpFilterParams}

%\input{TwoInterfFindChirpHFindChirpFilterParams}

\noindent This structure provides the parameters 
used by the \texttt{TwoInterfFindChirpFilterParams()} structure defined below.

\begin{description}
\item[\texttt{UINT4 length}] Number of \texttt{FindChirpFilterParams} 
structures in the vector. It is set equal to the number of detectors
in a network (in this case 2).

\item[\texttt{FindChirpFilterParams *data}] Pointer to the data.
\end{description}
</lalLaTeX>
#endif

/* parameter structure for the filtering function */
/* <lalVerbatim file="TwoInterfFindChirpHTwoInterfFindChirpFilterParams"> */
typedef struct
tagTwoInterfFindChirpFilterParams
{
  REAL4                                    twoInterfRhosqThresh;
  const LALDetectorPair                   *detectors;
  TwoInterfFindChirpFilterParamsVector    *paramsVec;
  REAL4TimeSeries                         *twoInterfRhosqVec;
}
TwoInterfFindChirpFilterParams;
/* </lalVerbatim> */
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{TwoInterfFindChirpFilterParams}}
\idx[Type]{TwoInterfFindChirpFilterParams}

%\input{TwoInterfFindChirpHTwoInterfFindChirpFilterParams}

\noindent This structure provides the parameters used by the
\texttt{TwoInterfFindChirpFilterSegment()} function.

\begin{description}
\item[\texttt{REAL4 twoInterfRhosqThresh}] The value to threshold signal to noise
ratio square, $\rho^2$, on. If the network signal to noise exceeds this value, then a
candidate event is generated. A $\chi^2$ veto is performed, if requested,
on each detector event; otherwise a network \texttt{InspiralEvent} is generated. 
Must be $\ge 0$ on entry.

\item[\texttt{const LALDetectorPair *detectors}] Pointer to the detector pair
whose data is filtered.

\item[\texttt{TwoInterfFindChirpFilterParamsVector *paramsVec}] Pointer to the
\texttt{TwoInterfFindChirpFilterParamsVector} structure defined above.

\item[\texttt{REAL4TimeSeries *twoInterfRhosqVec}] Pointer to a network vector that 
is set to $\rho^2(t_j)$ on exit. If NULL, $\rho^2(t_j)$ is not stored.
\end{description}
</lalLaTeX>
#endif

/*
 *
 * typedefs of input structures used by functions in twointerffindchirp
 *
 */


/* input to the filtering functions */
typedef struct
tagTwoInterfFindChirpFilterInputVector
{
  UINT4                         length;
  FindChirpFilterInput         *filterInput;
}
TwoInterfFindChirpFilterInputVector;
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{TwoInterfFindChirpFilterInputVector}}
\idx[Type]{TwoInterfFindChirpFilterInputVector}

%\input{TwoInterfFindChirpHTwoInterFindChirpFilterInputVector}

\noindent This structure groups the input data required for the
 \texttt{TwoInterfFindChirpFilterSegment()} function into a single structure.

\begin{description}
\item[\texttt{UINT4 length}] Number of \texttt{FindChirpFilterInput} 
structures in the vector. It is set equal to the number of detectors
in a network (in this case 2).

\item[\texttt{FindChirpFilterInput *filterInput}] Pointer to the input data.
\end{description}
</lalLaTeX>
#endif

typedef struct
tagTwoInterfFindChirpSPDataParamsVector
{
  UINT4                         length;
  FindChirpSPDataParams        *data;
}
TwoInterfFindChirpSPDataParamsVector;
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{TwoInterfFindChirpSPDataParamsVector}}
\idx[Type]{TwoInterfFindChirpSPDataParamsVector}

%\input{TwoInterfFindChirpHTwoInterfFindChirpSPDataParamsVector}

\noindent This structure contains the parameters needed to call the
\texttt{TwoInterfFindChirpSPData()} function. It should be initialized by
\texttt{TwoInterfFindChirpSPDataInit()} and destroyed by
\texttt{TwoInterfFindChirpSPDataFinalize()}. The fields are:


\begin{description}
\item[\texttt{UINT4 length}] Number of \texttt{FindChirpSPDataParams} 
structures in the vector. It is set equal to the number of detectors
in a network (in this case 2).

\item[\texttt{FindChirpSPDataParams  *data}] Pointer to the params data.
\end{description}
</lalLaTeX>
#endif


#if 0
<lalLaTeX>
\vfill{\footnotesize\input{TwoInterfFindChirpHV}}
</lalLaTeX> 
#endif


/*
 *
 * function prototypes for memory management functions
 .*
 */

void
LALCreateTwoInterfDataSegmentVector (
    LALStatus                       *status,
    TwoInterfDataSegmentVector     **vector,
    TwoInterfFindChirpInitParams    *params
    );

void
LALDestroyTwoInterfDataSegmentVector (
    LALStatus                   *status,
    TwoInterfDataSegmentVector **vector
    );

void
LALCreateTwoInterfFindChirpSegmentVector (
    LALStatus                        *status,
    TwoInterfFindChirpSegmentVector **vector,
    TwoInterfFindChirpInitParams     *params
    );

void
LALDestroyTwoInterfFindChirpSegmentVector (
    LALStatus                        *status,
    TwoInterfFindChirpSegmentVector **vector
    );


/*
 *
 * function prototypes for initialization, finalization and filter functions
 *
 */

void
LALTwoInterfFindChirpFilterInit (
    LALStatus                              *status,
    TwoInterfFindChirpFilterParams        **output,
    TwoInterfFindChirpInitParams           *params
    );

void
LALTwoInterfFindChirpFilterFinalize (
    LALStatus                             *status,
    TwoInterfFindChirpFilterParams       **output
    );

void
LALCreateTwoInterfFindChirpInputVector (
    LALStatus                                 *status,
    TwoInterfFindChirpFilterInputVector      **vector,
    TwoInterfFindChirpInitParams              *params
    );

void
LALDestroyTwoInterfFindChirpInputVector (
    LALStatus                                *status,
    TwoInterfFindChirpFilterInputVector     **vector
    );

void
LALTwoInterfFindChirpFilterSegment (
    LALStatus                            *status,
    TwoInterfInspiralEvent              **eventList,
    TwoInterfFindChirpFilterInputVector  *input,
    TwoInterfFindChirpFilterParams       *params
    );

void
LALTwoInterfFindChirpSPDataInit (
    LALStatus                                 *status,
    TwoInterfFindChirpSPDataParamsVector     **output,
    TwoInterfFindChirpInitParams              *params
    );

void
LALTwoInterfFindChirpSPData (
    LALStatus                                 *status,
    TwoInterfFindChirpSegmentVector           *twoInterfcSegVec,
    TwoInterfDataSegmentVector                *twoInterfDataSegVec,
    TwoInterfFindChirpSPDataParamsVector      *twoInterfDataParamsVec
    );

void
LALTwoInterfFindChirpSPDataFinalize (
    LALStatus                                 *status,
    TwoInterfFindChirpSPDataParamsVector     **output
    );

/* prototypes for chisq vetoing functions */
     
void
LALTwoInterfFindChirpChisqVetoInit (
    LALStatus                  *status,
    FindChirpChisqParams       *params,
    UINT4                       numChisqBins,
    UINT4                       numPoints
    );

void
LALTwoInterfFindChirpChisqVetoFinalize (
    LALStatus                  *status,
    FindChirpChisqParams       *params,
    UINT4                       numChisqBins
    );

void
LALTwoInterfFindChirpChisqVeto (
    LALStatus                  *status,
    REAL4Vector                *chisqVec,
    FindChirpChisqInput        *input,
    FindChirpChisqParams       *params
    );

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _TWOINTERFFINDCHIRPH_H */
