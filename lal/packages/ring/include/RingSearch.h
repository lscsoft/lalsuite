/**** <lalVerbatim file="RingSearchHV">
 * Author: Jolien Creighton
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 * 
 * \section{Header \texttt{RingSearch.h}}
 *
 * Black hole ringdown search code.
 *
 * \subsection*{Synopsis}
 * \begin{verbatim}
 * #include <lal/RingSearch.h>
 * \end{verbatim}
 * 
 * Routines for searching for black hole ringdown waveforms.
 *
 **** </lalLaTeX> */


#ifndef _RINGSEARCH_H
#define _RINGSEARCH_H

#include <lal/LALDatatypes.h>
#include <lal/RealFFT.h>
#include <lal/Ring.h>

#ifdef __cplusplus
extern "C" {
#pragma }
#endif

NRCSID( RINGSEARCHH, "$Id$" );

/**** <lalLaTeX>
 * \subsection*{Error conditions}
 **** </lalLaTeX> */
/**** <lalErrTable> */
#define RINGSEARCHH_ENULL 00001
#define RINGSEARCHH_ENNUL 00002
#define RINGSEARCHH_EALOC 00004
#define RINGSEARCHH_ESIZE 00010
#define RINGSEARCHH_ESZMM 00020
#define RINGSEARCHH_ENSEG 00040
#define RINGSEARCHH_EIOPT 00100
#define RINGSEARCHH_EFLOW 00200
#define RINGSEARCHH_EFREQ 00400
#define RINGSEARCHH_EQUAL 01000

#define RINGSEARCHH_MSGENULL "Null pointer"
#define RINGSEARCHH_MSGENNUL "Non-null pointer"
#define RINGSEARCHH_MSGEALOC "Memory allocation error"
#define RINGSEARCHH_MSGESIZE "Invalid segment size"
#define RINGSEARCHH_MSGESZMM "Size mismatch"
#define RINGSEARCHH_MSGENSEG "Non integer number of segments in data"
#define RINGSEARCHH_MSGEIOPT "Invalid option"
#define RINGSEARCHH_MSGEFLOW "Invalid low frequency cutoff"
#define RINGSEARCHH_MSGEFREQ "Invalid bank frequency range"
#define RINGSEARCHH_MSGEQUAL "Invalid bank quality range"
/**** </lalErrTable> */

/**** <lalLaTeX>
 * 
 * \subsection*{Structures}
 * 
 * \subsubsection*{Type \texttt{RingSearchParams}}
 * \idx[Type]{RingSearchParams}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct
tagRingSearchParams
{
  UINT4                    segmentSize;
  UINT4                    numSegments;
  COMPLEX8FrequencySeries *dataSegment;
  REAL4FrequencySeries    *invSpectrum;
  RealFFTPlan             *forwardPlan;
  RealFFTPlan             *reversePlan;
  REAL4                    dynRangeFac;
  UINT4                    invSpecTrunc;
  REAL4                    lowFrequency;
  REAL4                    minFrequency;
  REAL4                    maxFrequency;
  REAL4                    minQuality;
  REAL4                    maxQuality;
  REAL4                    maxMismatch;
  REAL4                    sampleRate;
  RingTemplateBank        *templateBank;
  UINT4                    templatesSent;
  UINT4                    templatesDone;
  INT4                     searchMaster;
  INT4                     myProcNumber;
  INT4                     numSlaves;
  REAL4                    threshold;
  CHAR                     ifoName[3];
  INT4                     maximizeEvents;
  INT4                     keepResults;
  UINT4                    numResults;
  REAL4TimeSeries         *result;
}
RingSearchParams;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 *
 * This structure contains the ring search parameters.  The structure is
 * primarily for internal use, but some of the fields can be modified.
 *
 * This structure is created by \verb+LALFindChirpInit()+ and is destroyed
 * by \verb+LALFindChirpFini()+, and modified as needed in the other routines.
 *
 * The fields are:
 * \begin{description}
 * \item[\texttt{segmentSize}]   The number of points in a segment of data.
 * \item[\texttt{numSegments}]   The number of segments of data to filter.
 * \item[\texttt{dataSegment}]   FFTed data segments.
 * \item[\texttt{invSpectrum}]   Inverse noise power spectrum.
 * \item[\texttt{forwardPlan}]   Forward FFT plan.
 * \item[\texttt{reversePlan}]   Reverse FFT plan.
 * \item[\texttt{dynRangeFac}]   Dynamical range factor used to scale strain.
 * \item[\texttt{invSpecTrunc}]  Number of points of the time-domain inverse
 *   spectrum response to keep.
 * \item[\texttt{lowFrequency}]  Low frequency (Hz) cutoff of data and spectrum.
 * \item[\texttt{minFrequency}]  Minimum ring frequency (Hz) for the bank.
 * \item[\texttt{maxFrequency}]  Maximum ring frequency (Hz) for the bank.
 * \item[\texttt{minQuality}]    Minimum ring quality for the bank.
 * \item[\texttt{maxQuality}]    Maximum ring quality for the bank.
 * \item[\texttt{maxMismatch}]   Maximum allowed mismatch for the bank.
 * \item[\texttt{sampleRate}]    Sample rate of the data.
 * \item[\texttt{templateBank}]  The template bank.
 * \item[\texttt{templatesSent}] For MPI code: number of templates sent.
 * \item[\texttt{templatesDone}] For MPI code: number of templates done.
 * \item[\texttt{searchMaster}]  For MPI code: non-zero if search master.
 * \item[\texttt{myProcNumber}]  For MPI code: MPI comm rank.
 * \item[\texttt{numSlaves}]     For MPI code: number of slaves.
 * \item[\texttt{threshold}]     SNR threshold for an event.
 * \item[\texttt{ifoName[3]}]    IFO name (e.g., ``H1,'' ``H2,'' ``L1,'' etc.).
 * \item[\texttt{maximizeEvents}] Non zero if events are maximized over arrival
 *   time.
 * \item[\texttt{keepResults}]   Non zero if SNR output is to be kept.
 * \item[\texttt{numResults}]    Number of SNR outputs kept.
 * \item[\texttt{result}]        The SNR outputs.
 * \end{description}
 *
 *
 * \subsubsection*{Type \texttt{RingSearchData}}
 * \idx[Type]{RingSearchData}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct
tagRingSearchData
{
  REAL4TimeSeries         *channel;
  REAL4FrequencySeries    *spectrum;
  COMPLEX8FrequencySeries *response;
}
RingSearchData;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 *
 * This structure contains the channel data, the noise power spectrum of the
 * channel data, and the response funtion.  It is used as input to the
 * \verb+LALRingSearchConditionData()+ function, which conditions and segments
 * the data, and computes the inverse strain noise spectrum.
 *
 * The fields are:
 *
 * \begin{description}
 * \item[\texttt{channel}]  The raw channel data.
 * \item[\texttt{spectrum}] The power spectrum of the channel data, or
 *   \verb+NULL+ if it is to be estimated from the data.
 * \item[\texttt{response}] The response function.
 * \end{description}
 *
 *
 * \subsubsection*{Type \texttt{RingEventList}}
 * \idx[Type]{RingEventList}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct
tagRingEventList
{
  CHAR  ifoName[3];
  INT8  startTimeNS;
  REAL4 duration;
  REAL4 amplitude;
  REAL4 frequency;
  REAL4 quality;
  REAL4 mass;
  REAL4 snr;
  REAL4 confidence;
  struct tagRingEventList *next;
}
RingEventList;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 *
 * This structure is a node on a linked list of events returned by the
 * \verb+LALRingSearch()+ function.  The fields are:
 *
 * \begin{description}
 * \item[\texttt{ifoName}] The IFO name (e.g., ``H1,'' ``L1,'' etc.).
 * \item[\texttt{startTimeNS}] The start time of the event in GPS nanoseconds.
 * \item[\texttt{duration}] The duration of the event in seconds.
 * \item[\texttt{amplitude}] The amplitude of the event.
 * \item[\texttt{frequency}] The central frequency of the ring filter.
 * \item[\texttt{quality}] The quality factor of the ring filter.
 * \item[\texttt{mass}] The mass of a black hole corresponding to the ring
 *   filter.
 * \item[\texttt{snr}] The amplitude signal-to-noise ratio of the event.
 * \item[\texttt{confidence}] The statistical confidence of this event.
 * \item[\texttt{next}] Pointer to the next element in the linked list.
 * \end{description}
 *
 *
 * \subsubsection*{Type \texttt{RingSearchInput}}
 * \idx[Type]{RingSearchInput}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct
tagRingSearchInput
{
  UINT4 startTemplate;
  UINT4 templatesToDo;
}
RingSearchInput;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 *
 * This structure is the input to the \verb+LALRingSearch()+ function.
 * The fields are:
 *
 * \begin{description}
 * \item[\texttt{startTemplate}] The number of the template in the template
 *   to use as the first filter.
 * \item[\texttt{templatesToDo}] The number of templates from the bank to
 *   to use as filters.
 * \end{description}
 *
 **** </lalLaTeX> */

void LALRingSearchInit(
    LALStatus         *status,
    RingSearchParams **searchParams,
    const CHAR       **argv,
    INT4               argc
    );

void
LALRingSearchFini(
    LALStatus         *status,
    RingSearchParams **searchParams
    );

void
LALRingSearchConditionData(
    LALStatus               *status,
    RingSearchParams        *params,
    RingSearchData          *data
    );

void
LALRingSearch(
    LALStatus         *status,
    RingEventList    **output,
    RingSearchInput   *input,
    RingSearchParams  *params
    );

/**** <lalLaTeX>
 * 
 * \vfill{\footnotesize\input{RingSearchHV}}
 * \newpage\input{RingSearchInitC}
 * \newpage\input{RingSearchConditionDataC}
 * \newpage\input{RingSearchC}
 * \newpage\input{RingSearchTestC}
 * 
 **** </lalLaTeX> */

#ifdef __cplusplus
#pragma {
}
#endif 

#endif /* _RINGSEARCH_H */
