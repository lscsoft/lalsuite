/**** <lalVerbatim file="StringSearchHV">
 * Author: Jocelyn Read
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 * 
 * \section{Header \texttt{StringSearch.h}}
 *
 * Black hole ringdown search code.
 *
 * \subsection*{Synopsis}
 * \begin{verbatim}
 * #include <lal/StringSearch.h>
 * \end{verbatim}
 * 
 * Routines for searching for black hole ringdown waveforms.
 *
 **** </lalLaTeX> */


#ifndef _STRINGSEARCH_H
#define _STRINGSEARCH_H

#include <lal/LALDatatypes.h>
#include <lal/RealFFT.h>
#include <lal/Window.h>
#include <lal/String.h>

#ifdef __cplusplus
extern "C" {
#pragma }
#endif

NRCSID( STRINGSEARCHH, "$Id$" );

/**** <lalLaTeX>
 * \subsection*{Error conditions}
 **** </lalLaTeX> */
/**** <lalErrTable> */
#define STRINGSEARCHH_ENULL 00001
#define STRINGSEARCHH_ENNUL 00002
#define STRINGSEARCHH_EALOC 00004
#define STRINGSEARCHH_ESIZE 00010
#define STRINGSEARCHH_ESZMM 00020
#define STRINGSEARCHH_ENSEG 00040
#define STRINGSEARCHH_EIOPT 00100
#define STRINGSEARCHH_EFLOW 00200
#define STRINGSEARCHH_EFREQ 00400
#define STRINGSEARCHH_ELSEG 02000

#define STRINGSEARCHH_MSGENULL "Null pointer"
#define STRINGSEARCHH_MSGENNUL "Non-null pointer"
#define STRINGSEARCHH_MSGEALOC "Memory allocation error"
#define STRINGSEARCHH_MSGESIZE "Invalid segment size"
#define STRINGSEARCHH_MSGESZMM "Size mismatch"
#define STRINGSEARCHH_MSGENSEG "Non integer number of segments in data"
#define STRINGSEARCHH_MSGEIOPT "Invalid option"
#define STRINGSEARCHH_MSGEFLOW "Invalid low frequency cutoff"
#define STRINGSEARCHH_MSGEFREQ "Invalid bank frequency range"
#define STRINGSEARCHH_MSGELSEG "Less than two segments in data"
/**** </lalErrTable> */

/**** <lalLaTeX>
 * 
 * \subsection*{Structures}
 * 
 * \subsubsection*{Type \texttt{StringSearchParams}}
 * \idx[Type]{StringSearchParams}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct
tagStringSearchParams
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
  REAL4			   freqPower;
  REAL4                    minFrequency;
  REAL4                    maxFrequency;
  REAL4                    maxMismatch;
  REAL4                    sampleRate;
  StringTemplateBank        *templateBank;
  UINT4                    templatesSent;
  UINT4                    templatesDone;
  INT4                     searchMaster;
  INT4                     myProcNumber;
  INT4                     numSlaves;
  UINT4                    numEvents;
  REAL4                    threshold;
  CHAR                     ifoName[3];
  INT4                     maximizeEvents;
  INT4                     keepResults;
  UINT4                    numResults;
  REAL4TimeSeries         *result;
}
StringSearchParams;
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
 * \item[\texttt{numEvents}]     Cumulative number of events found.
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
 * \subsubsection*{Type \texttt{StringSearchData}}
 * \idx[Type]{StringSearchData}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct
tagStringSearchData
{
  REAL4TimeSeries         *channel;
  REAL4FrequencySeries    *spectrum;
  COMPLEX8FrequencySeries *response;
}
StringSearchData;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 *
 * This structure contains the channel data, the noise power spectrum of the
 * channel data, and the response funtion.  It is used as input to the
 * \verb+LALLALStringSearchConditionData()+ function, which conditions and segments
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
 * \subsubsection*{Type \texttt{StringEventList}}
 * \idx[Type]{StringEventList}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct
tagStringEventList
{
  CHAR  ifoName[3];
  INT8  startTimeNS;
  REAL4 duration;
  REAL4 amplitude;
  REAL4 frequency;
  REAL4 mass;
  REAL4 snr;
  REAL4 confidence;
  struct tagStringEventList *next;
}
StringEventList;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 *
 * This structure is a node on a linked list of events returned by the
 * \verb+LALLALStringSearch()+ function.  The fields are:
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
 * \subsubsection*{Type \texttt{StringSearchInput}}
 * \idx[Type]{StringSearchInput}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct
tagStringSearchInput
{
  UINT4 startTemplate;
  UINT4 templatesToDo;
}
StringSearchInput;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 *
 * This structure is the input to the \verb+LALLALStringSearch()+ function.
 * The fields are:
 *
 * \begin{description}
 * \item[\texttt{startTemplate}] The number of the template in the template
 *   to use as the first filter.
 * \item[\texttt{templatesToDo}] The number of templates from the bank to
 *   to use as filters.
 * \end{description}
 *
 * \subsubsection*{Type \texttt{AvgSpecParams}}
 * \idx[Type]{AvgSpecParams}
 *
 **** </lalLaTeX> */
#ifndef AVG_SPEC_DECLARED
#define AVG_SPEC_DECLARED
/**** <lalVerbatim> */
typedef struct
tagAvgSpecParams
{
  UINT4        segsize;
  RealFFTPlan *fwdplan;
  WindowType   wintype;
}
AvgSpecParams;
/**** </lalVerbatim> */
#endif
/**** <lalLaTeX>
 *
 * This structure contains parameters for routines that compute average
 * power spectra.  The fields are:
 *
 * \begin{description}
 * \item[\texttt{segsize}] The size of each segment of data.
 * \item[\texttt{fwdplan}] Forward real FFT plan for that segment size.
 * \item[\texttt{wintype}] type of window to apply to FFT.
 * \end{description}
 *
 **** </lalLaTeX> */

void LALStringSearchInit(
    LALStatus         *status,
    StringSearchParams **searchParams,
    const CHAR       **argv,
    INT4               argc
    );

void
LALStringSearchFini(
    LALStatus         *status,
    StringSearchParams **searchParams
    );

void
LALStringSearchConditionData(
    LALStatus               *status,
    StringSearchParams        *params,
    StringSearchData          *data
    );

void
LALStringSearch(
    LALStatus         *status,
    StringEventList    **output,
    StringSearchInput   *input,
    StringSearchParams  *params
    );

/**** <lalLaTeX>
 * 
 * \vfill{\footnotesize\input{StringSearchHV}}
 * \newpage\input{StringSearchInitC}
 * \newpage\input{StringSearchConditionDataC}
 * \newpage\input{StringSearchC}
 * \newpage\input{StringSearchTestC}
 * 
 **** </lalLaTeX> */

#ifdef __cplusplus
#pragma {
}
#endif 

#endif /* _STRINGSEARCH_H */
