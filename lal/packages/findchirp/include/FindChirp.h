/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirp.h
 *
 * Author: Allen, B., Brown, D. A. and Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpHV">
Author: Allen, B., Brown, D. A. and Creighton, J. D. E.
$Id$
</lalVerbatim> 

<lalLaTeX>
\section{Header \texttt{FindChirp.h}}
\label{s:FindChirp.h}

Provides core protypes, structures and functions to filter interferometer data
for binary inspiral chirps.

\subsection*{Synopsis}

\begin{verbatim}
#include <lal/FindChirp.h>
\end{verbatim}

\input{FindChirpHDoc}
</lalLaTeX>
#endif

#ifndef _FINDCHIRPH_H
#define _FINDCHIRPH_H

#include <lal/LALDatatypes.h>
#include <lal/ComplexFFT.h>
#include <lal/DataBuffer.h>
#ifdef LAL_ENABLE_MPI
#include <lal/Comm.h>
#endif
#include <lal/LALInspiral.h>
#include <lal/FindChirpChisq.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif


NRCSID (FINDCHIRPH, "$Id$");

#if 0
<lalLaTeX> 
\subsection*{Error codes} 
</lalLaTeX>
#endif
/* <lalErrTable> */
#define FINDCHIRPH_ENULL 1
#define FINDCHIRPH_ENNUL 2
#define FINDCHIRPH_EALOC 3
#define FINDCHIRPH_ENUMZ 5
#define FINDCHIRPH_ESEGZ 6
#define FINDCHIRPH_ECHIZ 7
#define FINDCHIRPH_EDTZO 8
#define FINDCHIRPH_ETRNC 10
#define FINDCHIRPH_EFLOW 11
#define FINDCHIRPH_EFREE 12
#define FINDCHIRPH_ERHOT 15
#define FINDCHIRPH_ECHIT 16
#define FINDCHIRPH_ECRUP 17
#define FINDCHIRPH_MSGENULL "Null pointer"
#define FINDCHIRPH_MSGENNUL "Non-null pointer"
#define FINDCHIRPH_MSGEALOC "Memory allocation error"
#define FINDCHIRPH_MSGENUMZ "Invalid number of points in segment"
#define FINDCHIRPH_MSGESEGZ "Invalid number of segments"
#define FINDCHIRPH_MSGECHIZ "Invalid number of chi squared bins"
#define FINDCHIRPH_MSGEDTZO "deltaT is zero or negative"
#define FINDCHIRPH_MSGETRNC "Duration of inverse spectrum in time domain is negative"
#define FINDCHIRPH_MSGEFLOW "Inverse spectrum low frequency cutoff is negative"
#define FINDCHIRPH_MSGEFREE "Error freeing memory"
#define FINDCHIRPH_MSGERHOT "Rhosq threshold is zero or negative"
#define FINDCHIRPH_MSGECHIT "Chisq threshold is zero or negative"
#define FINDCHIRPH_MSGECRUP "Chirp length or invSpecTrunc too long for length of data segment"
/* </lalErrTable> */


#if 0
<lalLaTeX>
\subsection*{Types}
</lalLaTeX>
#endif

/* --- enumeraion type for simulations ----------------------------------- */
#pragma <lalVerbatim file="FindChirpHWaveformGenerator">
typedef enum
{
  findchirpSP,
  injectGenPPN
}
WaveformGenerator;
#pragma </lalVerbatim>
#if 0
<lalLaTeX>
\subsubsection*{Enumeration Type \texttt{WaveformGenerator}}
\idx[Type]{WaveformGenerator}

\input{FindChirpHWaveformGenerator}

\noindent This enumaration type lists the possible methods of generating
inspiral waveforms. The types in the enum should be of the form 
\texttt{packagenameGenType}, where \texttt{packagename} is the name of the
LAL package that the particular waveform generator resides and \texttt{GenType}
is a unique identifier for that waveform generator.  The choices are:

\begin{description} 
\item[\texttt{findchirpSP}] The \texttt{findchirp} built in stationary phase
waveform generator.

\item[\texttt{injectGenPPN}] The \texttt{LALGeneratePPNInspiral()} function
from the \texttt{inject} package.

\item[\texttt{fctFCT}] The \texttt{fct} algorithm was used to find the chirp.
\end{description}
</lalLaTeX>
#endif
  

/*
 *
 * typedefs of structures used by the findchirp functions
 *
 */


#if 0
<lalLaTeX>
\subsubsection*{Input and output data structures}
</lalLaTeX>
#endif

/* --- structure for describing a binary insipral event ------------------ */
#pragma <lalVerbatim file="FindChirpHInspiralEvent">
typedef struct
tagInspiralEvent
{
  UINT4                         id;
  UINT4                         segmentNumber;
  LIGOTimeGPS                   time;
  LIGOTimeGPS                   impulseTime;
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
#pragma </lalVerbatim>
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

/* --- vector of DataSegment, as defined the framedata package ----------- */
#pragma <lalVerbatim file="FindChirpHDataSegmentVector">
typedef struct
tagDataSegmentVector
{
  UINT4                         length;
  DataSegment                  *data;
}
DataSegmentVector;
#pragma </lalVerbatim>
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{DataSegmentVector}}
\idx[Type]{DataSegmentVector}

\input{FindChirpHDataSegmentVector}

\noindent This structure provides a LAL like vector structure for the
\texttt{DataSegment} structure defined in the package \texttt{framedata}

\begin{description}
\item[\texttt{UINT4 length}] Number of \texttt{DataSegment} structres in the
vector

\item[\texttt{DataSegment *data}] Pointer to the data.
\end{description}
</lalLaTeX>
#endif

/* --- processed data segment used by FindChirp filter routine ----------- */
#pragma <lalVerbatim file="FindChirpHFindChirpSegment">
typedef struct
tagFindChirpSegment
{
  COMPLEX8FrequencySeries      *data;
  UINT4Vector                  *chisqBinVec;
  REAL8                         deltaT;
  REAL4                         segNorm;
  REAL4                         fLow;
  UINT4                         invSpecTrunc;
  UINT4                         number;
}
FindChirpSegment;
#pragma </lalVerbatim>
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{FindChirpSegment}}
\idx[Type]{FindChirpSegment}

\input{FindChirpHFindChirpSegment}

\noindent This structure contains the conditioned input data and its
parameters for the \texttt{FindChirpFilter()} function.

\begin{description}
\item[\texttt{COMPLEX8FrequencySeries *data}] The conditioned input data.
The conditioneing performed is as described in the documentation for the
module \texttt{FindChirpSPData.c}

\item[\texttt{UINT4Vector *chisqBinVec}] A vector containing the indices of
the boundaries of the bins of equal power for the $\chi^2$ veto created by 
\texttt{FindChirpSPData()}

\item[\texttt{REAL8 deltaT}] The time step $\Delta$ of the time series 
input data.

\item[\texttt{REAL4 segNorm}] The template independent part of the 
normalisation constant $\sigma$.

\item[\texttt{UINT4 invSpecTrunc}] The number of points to which the inverse 
power spectrum \ospsd is truncated to in the time domain in order to smooth
out high $Q$ features in the power spectrum.

\item[\texttt{UINT4 number}] A unique identification number for the 
\texttt{FindChirpDataSegment}. This will generally correspond to the number in
the \texttt{DataSegment} from which the conditioned data was computed.
\end{description}
</lalLaTeX>
#endif

/* --- vector of FindChirpSegment defined above -------------------------- */
#pragma <lalVerbatim file="FindChirpHFindChirpSegmentVector">
typedef struct
tagFindChirpSegmentVector
{
  UINT4                         length;
  FindChirpSegment             *data;
}
FindChirpSegmentVector;
#pragma </lalVerbatim>
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{FindChirpSegmentVector}}
\idx[Type]{FindChirpSegmentVector}

\input{FindChirpHFindChirpSegmentVector}

\noindent This structure provides a LAL like vector structure for the
\texttt{FindChirpSegment} structure defined above.

\begin{description}
\item[\texttt{UINT4 length}] Number of \texttt{FindChirpSegment} structres in
the vector

\item[\texttt{DataSegment *data}] Pointer to the data.
\end{description}
</lalLaTeX>
#endif

/* --- structure to contain an inspiral template ------------------------- */
#pragma <lalVerbatim file="FindChirpHFindChirpTemplate">
typedef struct
tagFindChirpTemplate
{
  COMPLEX8Vector               *data;
  REAL4                         tmpltNorm;
}
FindChirpTemplate;
#pragma </lalVerbatim>
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{FindChirpTemplate}}
\idx[Type]{FindChirpTemplate}

\input{FindChirpHFindChirpTemplate}

\noindent This structure provides contains the frequency domain representation
of the cosine phase inspiral template $\tilde{h_c}(f)$.

\begin{description}
\item[\texttt{COMPLEX8Vector *data}] A vector containing $\tilde{h_c}(f)$. Note
that in the future, this will be changed to a \texttt{COMPLEX8FrequencySeries}.

\item[\texttt{REAL4 tmpltNorm}] The template dependent part of the 
normalisation constant $\sigma$.
\end{description}
</lalLaTeX>
#endif


/*
 *
 * typedefs of parameter structures used by functions in findchirp
 *
 */


#if 0
<lalLaTeX>
\subsubsection*{Initalisation and parameter structures}
</lalLaTeX>
#endif
  

/* --- parameter structure for all init funtions ------------------------- */
#pragma <lalVerbatim file="FindChirpHFindChirpInitParams">
typedef struct
tagFindChirpInitParams
{
  UINT4                         numSegments;
  UINT4                         numPoints;
  UINT4                         numChisqBins;
  BOOLEAN                       createRhosqVec;
}
FindChirpInitParams;
#pragma </lalVerbatim>
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{FindChirpInitParams}}
\idx[Type]{FindChirpInitParams}

\input{FindChirpHFindChirpInitParams}

\noindent This structure provides the essential information for the
filter initialisation and memory allocation functions used in the
\texttt{FindChirp} package.

\begin{description}
\item[\texttt{UINT4 numSegments}] The number of data segments to allocate
storage for.

\item[\texttt{UINT4 numPoints}] The number of discrete data points in each
data segment. 

\item[\texttt{UINT4 numChisqBins}] The number of bins used to contruct the
$\chi^2$ veto.

\item[\texttt{BOOLEAN createRhosqVec}] Debugging flag that controls whether
or not the function \texttt{FindChirpFilterSegment()} should store the output
of the filter, $\rho^2(t)$, as well as the events. Memory is only allocated
for this vector if the flag is set to 1.
\end{description}
</lalLaTeX>
#endif


/* --- parameter structure for the filtering function -------------------- */
#pragma <lalVerbatim file="FindChirpHFindChirpFilterParams">
typedef struct
tagFindChirpFilterParams
{
  REAL4                         deltaT;
  REAL4                         rhosqThresh;
  REAL4                         chisqThresh;
  REAL4                         norm;
  UINT4                         maximiseOverChirp;
  BOOLEAN                       computeNegFreq;
  COMPLEX8Vector               *qVec;
  COMPLEX8Vector               *qtildeVec;
  ComplexFFTPlan               *invPlan;
  REAL4Vector                  *rhosqVec;
  REAL4Vector                  *chisqVec;
  FindChirpChisqParams         *chisqParams;
  FindChirpChisqInput          *chisqInput;
}
FindChirpFilterParams;
#pragma </lalVerbatim>
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{FindChirpFilterParams}}
\idx[Type]{FindChirpFilterParams}

\input{FindChirpHFindChirpFilterParams}

\noindent This structure provides the parameters used by the
\texttt{FindChirpFilterSegment()} function.

\begin{description}
\item[\texttt{REAL4 deltaT}] The timestep for the sampled data. Must be
set on entry.  {FIXME: This should be a \texttt{REAL8}}

\item[\texttt{REAL4 rhosqThresh}] The value to threshold signal to noise
ratio square, $\rho^2$, on. If the signal to noise exceeds this value, then a
candidate event is generated. A $\chi^2$ veto is performed, if requested,
otherwise an \texttt{InspiralEvent} is generated. Must be $\ge 0$ on entry.

\item[\texttt{REAL4 chisqThresh}] The value to threshold the $\chi^2$ veto on.
If the $chi^2$ veto is below this threshold the candidate event an
\texttt{InspiralEvent} is generated. Must be $\ge 0$ on entry.

\item[\texttt{REAL4 norm}] On exit this contains the normalisation constant
that relates the quantity $q_j$ with the signal to noise squared, 
$\rho^2(t_j)$ by 
\begin{equation}
\rho^2(t_j) = \textrm{norm} \times \left|q_j\right|^2.
\end{equation}

\item[\texttt{UINT4 maximiseOverChirp}] If not zero, use algorithm that
maximised over chirp lengths. Otherwise record all points that pass
the $\rho^2$ threshold as events.

\item[\texttt{BOOLEAN computeNegFreq}] Currently unused. Must be set to
$0$ on entry.

\item[\texttt{COMPLEX8Vector *qVec}] Pointer to vector allocated by 
\texttt{FindChirpFilterInit()} to store the quantity $q_j$. Set to the
value of $q_j$ on exit. Must not be NULL.

\item[\texttt{COMPLEX8Vector *qtildeVec}] Pointer to vector allocated by 
\texttt{FindChirpFilterInit()} to store the quantity $\tilde{q}_k$. Set to the
value of $\tilde{q}_k$ on exit. Must not be NULL

\item[\texttt{ComplexFFTPlan *invPlan}] Pointer to FFTW plan created by 
\texttt{FindChirpFilterInit()} to transform the quantity $\tilde{q}_k$ to
${q}_j$ usimg the inverse DFT. Must not be NULL.

\item[\texttt{REAL4Vector *rhosqVec}] Pointer to a vector that is set to
$\rho^2(t_j)$ on exit. If NULL $\rho^2(t_j)$ is not stored.

\item[\texttt{REAL4Vector *chisqVec}] Workspace vector used to compute and
store $\chi^2(t_j)$. Must not be NULL if \texttt{numChisqBins} is greater than
zero. Contains $\chi^2(t_j)$ on exit.

\item[\texttt{FindChirpChisqParams *chisqParams}] Pointer to parameter
structure for \texttt{FindChirpChisqVeto()} function. Must not be NULL if
\texttt{numChisqBins} is greater than zero.

\item[\texttt{FindChirpChisqInput *chisqInput}] Pointer to input data
structure for \texttt{FindChirpChisqVeto()} function. Must not be NULL if
\texttt{numChisqBins} is greater than zero.
\end{description}
</lalLaTeX>
#endif


/*
 *
 * typedefs of input structures used by functions in findchirp
 *
 */


#if 0
<lalLaTeX>
\subsubsection*{Filter function input structures}
</lalLaTeX>
#endif


/* --- input to the filtering functions --------------------------------- */
#pragma <lalVerbatim file="FindChirpHFindChirpFilterInput">
typedef struct
tagFindChirpFilterInput
{
  InspiralTemplate             *tmplt;
  FindChirpTemplate            *fcTmplt;
  FindChirpSegment             *segment;
}
FindChirpFilterInput;
#pragma </lalVerbatim>
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{FindChirpFilterInput}}
\idx[Type]{FindChirpSegmentVector}

\input{FindChirpHFindChirpFilterInput}

\noindent This structure groups the input data required for the 
\texttt{FindChirpFilterSegment()} function into a single structure.

\begin{description}
\item[\texttt{InspiralTemplate *tmplt}] Pointer the structure that contains
the parameters of the template chirp.

\item[\texttt{FindChirpTemplate *fcTmplt}] Pointer to the input template
in a form that can be used by \texttt{FindChirpFilterSegment()}

\item[\texttt{FindChirpSegment *segment}] Pointer to the input data segment
in a form that can be used by \texttt{FindChirpFilterSegment()}
\end{description}
</lalLaTeX>
#endif

#if 0
<lalLaTeX>
\vfill{\footnotesize\input{FindChirpHV}}
</lalLaTeX> 
#endif

/*
 *
 * function prototypes for memory management functions
 *
 */


#if 0
<lalLaTeX>
\newpage\input{FindChirpMemoryC}
</lalLaTeX>
#endif

void
LALCreateDataSegmentVector (
    LALStatus                  *status,
    DataSegmentVector         **vector,
    FindChirpInitParams        *params
    );

void
LALDestroyDataSegmentVector (
    LALStatus                  *status,
    DataSegmentVector         **vector
    );

void
LALCreateFindChirpSegmentVector (
    LALStatus                  *status,
    FindChirpSegmentVector    **vector,
    FindChirpInitParams        *params
    );

void
LALDestroyFindChirpSegmentVector (
    LALStatus                  *status,
    FindChirpSegmentVector    **vector
    );


/*
 *
 * function prototypes for initialization, finalization and filter functions
 *
 */


#if 0
<lalLaTeX>
\newpage\input{FindChirpFilterC}
</lalLaTeX>
#endif

void
LALFindChirpFilterInit (
    LALStatus                  *status,
    FindChirpFilterParams     **output,
    FindChirpInitParams        *params
    );

void
LALFindChirpFilterFinalize (
    LALStatus                  *status,
    FindChirpFilterParams     **output
    );

void
LALCreateFindChirpInput (
    LALStatus                  *status,
    FindChirpFilterInput      **output,
    FindChirpInitParams        *params
    );

void
LALDestroyFindChirpInput (
    LALStatus                  *status,
    FindChirpFilterInput      **output
    );

void
LALFindChirpFilterSegment (
    LALStatus                  *status,
    InspiralEvent             **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params
    );



#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _FINDCHIRPH_H */
