/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpEngine.h
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpEngineHV">
Author: Brown, D. A.
$Id$
</lalVerbatim> 

<lalLaTeX>
\section{Header \texttt{FindChirpEngine.h}}
\label{s:FindChirpEngine.h}

Provides functions to drive the filtering functions in the package
\texttt{findchirp}. Contains the functions necessary to execute a flat
or heirarchical inspiral search using the LDAS \texttt{wrapperAPI} or
the LALapps program \texttt{findchirp}.

\subsection*{Synopsis}

\begin{verbatim}
#include <lal/FindChirpEngine.h>
\end{verbatim}

\input{FindChirpEngineHDoc}

</lalLaTeX>
#endif

#ifndef _FINDCHIRPENGINEH_H
#define _FINDCHIRPENGINEH_H

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <lal/LALStdlib.h>
#include <lal/Random.h>
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/LALInspiralBank.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/Units.h>
#include <lal/Inject.h>
#ifdef LAL_MPI_ENABLED
#include <lal/Comm.h>
#include <lal/FindChirpExch.h>
#endif /* LAL_MPI_ENABLED */


#ifdef  __cplusplus
extern "C" {
#pragma }
#endif


NRCSID (FINDCHIRPENGINEHH, "$Id$");

#if 0
<lalLaTeX> 
\subsection*{Error codes} 
</lalLaTeX>
#endif
/* <lalErrTable> */
#define FINDCHIRPENGINEH_ENULL 1
#define FINDCHIRPENGINEH_ENNUL 2
#define FINDCHIRPENGINEH_ENUMZ 3
#define FINDCHIRPENGINEH_EDELT 4
#define FINDCHIRPENGINEH_ERHOZ 5
#define FINDCHIRPENGINEH_ECHIZ 6
#define FINDCHIRPENGINEH_ETMPL 7
#define FINDCHIRPENGINEH_EALOC 8
#define FINDCHIRPENGINEH_ERANK 9
#define FINDCHIRPENGINEH_EUEXT 10
#define FINDCHIRPENGINEH_ELVEL 11
#define FINDCHIRPENGINEH_ESEGZ 12
#define FINDCHIRPENGINEH_EUSIM 13
#define FINDCHIRPENGINEH_EWAVL 14
#define FINDCHIRPENGINEH_MSGENULL "Null pointer"
#define FINDCHIRPENGINEH_MSGENNUL "Non-null pointer"
#define FINDCHIRPENGINEH_MSGENUMZ "Data segment length is zero"
#define FINDCHIRPENGINEH_MSGEDELT "deltaT is zero or negative"
#define FINDCHIRPENGINEH_MSGERHOZ "snr squared threshold is zero or negative"
#define FINDCHIRPENGINEH_MSGECHIZ "chi squared threshold is zero or negative"
#define FINDCHIRPENGINEH_MSGETMPL "linked list of templates to filter in null"
#define FINDCHIRPENGINEH_MSGEALOC "Memory allocation error"
#define FINDCHIRPENGINEH_MSGERANK "Search node has incorrect rank"
#define FINDCHIRPENGINEH_MSGEUEXT "Unrecognised exchange type"
#define FINDCHIRPENGINEH_MSGELVEL "Invalid heriarchical template bank level"
#define FINDCHIRPENGINEH_MSGESEGZ "Number of data segments is zero"
#define FINDCHIRPENGINEH_MSGEUSIM "Unkown simulation type requested"
#define FINDCHIRPENGINEH_MSGEWAVL "Simulated waveform is longer than dataseg"
/* </lalErrTable> */


#if 0
<lalLaTeX>
\subsection*{Types}
</lalLaTeX>
#endif

/* --- enumeraion type for simulations ----------------------------------- */
#pragma <lalVerbatim file="FindChirpEngineHFindChirpSimulationType">
typedef enum
{
  fcNoSim,
  fcCreateRhosqVec,
  fcGaussianNoise,
  fcGaussianNoiseInject,
  fcRealDataInject,
  fcBankMinimalMatch
}
FindChirpSimulationType;
#pragma </lalVerbatim>
#if 0
<lalLaTeX>
\subsubsection*{Enumeration Type \texttt{FindChirpSimulationType}}
\idx[Type]{FindChirpSimulationType}

\input{FindChirpEngineHFindChirpSimulationType}

\noindent This enumaration type lists the possible types of simulation
that can be performed by the \texttt{FindChirpEngine} functions.  The choices
are:

\begin{description}
\item[\texttt{fcNoSim}] A regular data analysis run; no simulation performed.

\item[\texttt{fcCreateRhosqVec}] Create a vector contain the statistic
$\rho^2(t_j)$ when filtering.

\item[\texttt{fcGaussianNoise}] Replace the input data in the
\texttt{DataSegmentVector} with Gaussian noise and a simulated LIGO I response
function.

\item[\texttt{fcGaussianNoiseInject}] Replace the input data in the
\texttt{DataSegmentVector} with Gaussian noise and a simulated LIGO I response
function. Inject random signals into the data as descibed by the
\texttt{FindChirpSimulationParams}.

\item[\texttt{fcGaussianNoiseInject}] Inject random signals into the input
data as descibed by the \texttt{FindChirpSimulationParams}.

\item[\texttt{fcBankMinimalMatch}] Perform a bank minimal match simulation
using the parameters in the \texttt{FindChirpSimulationParams} structure.
\end{description}
</lalLaTeX>
#endif


/* --- structure for managing a list of inspiral templates --------------- */
#pragma <lalVerbatim file="FindChirpEngineHInspiralTemplateNode">
typedef struct
tagInspiralTemplateNode
{
  struct tagInspiralTemplateNode       *next;
  struct tagInspiralTemplateNode       *prev;
  InspiralTemplate                     *tmpltPtr;
}
InspiralTemplateNode;
#pragma </lalVerbatim>
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{InspiralTemplateNode}}
\idx[Type]{InspiralTemplateNode}

\input{FindChirpEngineHInspiralTemplateNode}

\noindent This structure provides a method of constucting doubly linked
lists of \texttt{InspiralTemplate} structures. The fields are:

\begin{description}
\item[\texttt{struct tagInspiralTemplateNode *next}] The next structure in
the linked list.

\item[\texttt{struct tagInspiralTemplateNode *prev}] The previous structure in
the linked list.

\item[\texttt{InspiralTemplate *tmpltPtr}] A pointer to an \texttt{InspiralTemplate} structure.
\end{description}
</lalLaTeX>
#endif

typedef struct
tagInspiralInjection
{
  INT8                          coaTime;
  REAL4                         totalMass;
  REAL4                         eta;
  REAL4                         dist;
  REAL4                         inclination;
  REAL4                         coaPhase;
  REAL4                         longitude;
  REAL4                         latitude;
  REAL4                         polarization;
  struct tagInspiralInjection  *next;
}
InspiralInjection;

typedef struct
tagFindChirpSimulationParams
{
  FindChirpSimulationType       simType;        /* type of simulation        */
  UINT4                         simCount;       /* number of simulations     */
  UINT4                         injectStartIdx; /* index of the first inj    */
  InspiralInjection            *injectEvent;    /* events to inject          */
  InspiralEvent                *loudestEvent;   /* array of loudest events   */
  REAL4                        *signalNorm;     /* (s|s) for each dataSeg    */
  REAL4                         mMin;           /* minimum mass for binary   */
  REAL4                         mMax;           /* maximum mass for binary   */
  REAL4                         fLowInject;     /* waveform low freq cutoff  */
  REAL4                         gaussianVarsq;  /* variance squared of noise */
  RandomParams                 *randomParams;   /* random seed container     */
  REAL4TimeSeries              *chan;           /* pointer to the data chunk */
}
FindChirpSimulationParams;

typedef struct
tagFindChirpStandardCandle
{
  CHAR                          ifo[2];
  InspiralTemplate              tmplt;
  REAL4                         rhosq;
  REAL4                         sigmasq;
  REAL4                         effDistance;
}
FindChirpStandardCandle;

typedef struct
tagFindChirpSlaveParams
{
  UINT4                         dataConditioned;
  UINT4                        *inspiralDebugFlagPtr;
  REAL4                        *rhosqThreshVec;
  REAL4                        *chisqThreshVec;
  FILE                         *tmpltBankFilePtr;
  FILE                         *eventFilePtr;
  FindChirpStandardCandle      *candlePtr;
  FindChirpSegmentVector       *fcSegVec;
  FindChirpSPDataParams        *dataParams;
  FindChirpSPTmpltParams       *tmpltParams;
  FindChirpFilterParams        *filterParams;
  FindChirpFilterInput         *filterInput;
  FindChirpSimulationParams    *simParams;
  BOOLEAN                      *notFinished;
  BOOLEAN                       useMPI;
  void                         *mpiComm;
}
FindChirpSlaveParams;

typedef struct
tagFindChirpCreateBankParams
{
  INT4                          numCoarse;
  UINT4                         numSegments;
  UINT4                         numLevel;
}
FindChirpCreateBankParams;

/* Only define the master and exchange type if we have MPI enabled */
#ifdef LAL_MPI_ENABLED
enum ExchObjectType
{
  ExchDataSegment,
  ExchFindChirpSegment,
  ExchInspiralTemplate,
  ExchInspiralEvent,
  ExchNumTmpltsFiltered,
  ExchFinished
};

typedef struct
tagFindChirpMasterParams
{
  UINT4                         numTmpltExch;
  UINT4                         numTmpltsTotal;
  UINT4                         numTmpltsToFilter;
  UINT4                         numTmpltsFiltered;
  UINT4                        *inspiralDebugFlagPtr;
  UINT4Vector                  *bankSentVec;
  MPI_Comm                     *mpiComm;
  UINT4                        *numSlaves;
  InspiralTemplate             *tmpltBankHead;
  InspiralTemplateNode         *currentTmpltNode;
  InspiralTemplateNode         *tmpltNodeHead;
}
FindChirpMasterParams;

void
LALFindChirpMaster (
    LALStatus                  *status, 
    InspiralEvent             **eventList,
    FindChirpMasterParams      *params 
    );
#endif /* LAL_MPI_ENABLED */

void
LALFindChirpSlave (
    LALStatus                  *status, 
    InspiralEvent             **outputEventHandle,
    DataSegmentVector          *dataSegVec,
    FindChirpSlaveParams        *params 
    );

void
LALFindChirpCreateInspiralBank (
    LALStatus                  *status,
    InspiralCoarseBankIn       *bankIn,
    InspiralTemplate          **bankHead,
    FindChirpCreateBankParams  *params
    );

void
LALFindChirpDestroyInspiralBank (
    LALStatus                  *status,
    InspiralTemplate          **bankHead
    );

void
LALFindChirpCreateTmpltNode (
    LALStatus                  *status,
    InspiralTemplate           *tmplt,
    InspiralTemplateNode      **tmpltNode
    );

void
LALFindChirpDestroyTmpltNode ( 
    LALStatus                  *status,
    InspiralTemplateNode      **tmpltNode
    );

void
LALFindChirpInjectSignals (
    LALStatus                  *status,
    REAL4TimeSeries            *chan,
    InspiralInjection          *events,
    COMPLEX8FrequencySeries    *resp
    );

#ifdef  __cplusplus
}
#endif

#endif /* _FINDCHIRPENGINEH_H */
