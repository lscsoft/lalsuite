/*********************** <lalVerbatim file="TwoInterfFindChirpHV">
Author: Sukanta Bose
$Id$ 
*********************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>
\section{Header \texttt{TwoInterfFindChirp.h}}
\label{twointerffindchirp:s:TwoInterfFindChirp.h}

Provides prototype and error code information for the modules needed
to search for an inspiral waveform in the outputs of a pair of 
interferometers.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/TwoInterfFindChirp.h> 
\end{verbatim}

\noindent 

\subsection*{Error conditions}
\input{TwoInterfFindChirpHE}

\subsection*{Structures}

*********************************************************** </lalLaTeX> */

#include <lal/LALDatatypes.h>
#include <lal/ComplexFFT.h>
#include <lal/DataBuffer.h>
#include <lal/Comm.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirpChisq.h>
#include <lal/DetectorSite.h>
#include <lal/StochasticCrossCorrelation.h>

#ifdef  __cplusplus
extern "C" {
#endif
  
  
  NRCSID (TWOINTERFFINDCHIRPH, "$Id$");
  
/***************** <lalErrTable file="TwoInterfFindChirpHE"> */
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
#define TWOINTERFFINDCHIRPH_MSGENULL "Null pointer"
#define TWOINTERFFINDCHIRPH_MSGENNUL "Non-null pointer"
#define TWOINTERFFINDCHIRPH_MSGEALOC "Memory allocation error"
#define TWOINTERFFINDCHIRPH_MSGENUMZ "Invalid number of points in segment"
#define TWOINTERFFINDCHIRPH_MSGESEGZ "Invalid number of segments"
#define TWOINTERFFINDCHIRPH_MSGECHIZ "Invalid number of chi squared bins"
#define TWOINTERFFINDCHIRPH_MSGEDTZO "deltaT is zero or negative"
#define TWOINTERFFINDCHIRPH_MSGEMLZO "maxLag is zero or negative"
#define TWOINTERFFINDCHIRPH_MSGETRNC "Duration of inverse spectrum in time domain is negative"
#define TWOINTERFFINDCHIRPH_MSGEFLOW "Inverse spectrum low frequency cutoff is negative"
#define TWOINTERFFINDCHIRPH_MSGEFREE "Memory free error"
#define TWOINTERFFINDCHIRPH_MSGENUMF "Invalid number of points in filter"
#define TWOINTERFFINDCHIRPH_MSGERHOT "Rhosq threshold is zero or negative"
#define TWOINTERFFINDCHIRPH_MSGECHIT "Chisq threshold is zero or negative"

/************************************ </lalErrTable> */

  /*************************************************************
   *                                                           *
   *       Structures and prototypes associated with           *
   *             TwoInterfFindChirp.c                          *
   *                                                           *
   *************************************************************/

/********************************************************** <lalLaTeX>

\subsubsection*{Structures associated with 
  \texttt{TwoInterfFindChirp.c}
  (Sec.~\ref{twointerffindchirp:ss:TwoInterfFindChirp.c})}

\subsubsection*{\texttt{struct TwoInterfInspiralEvent}}
\idx[Type]{TwoInterfInspiralEvent}

\noindent Contains the list and description of events detected by 
\texttt{LALTwoInterfFindChirp()}. 

\begin{description}
\item[\texttt{UINT4  twoInterfId}]
\item[\texttt{REAL4  snrsq}]
\item[\texttt{InspiralTemplate  tmplt}]
\item[\texttt{InspiralEvent  *eventIn1}]
\item[\texttt{InspiralEvent  *eventIn2}]
\end{description}

*********************************************************** </lalLaTeX> */

/*
 *
 * typedefs of structures used by the twointerffindchirp functions
 *
 */


/* structure for describing a binary insipral event */
  typedef struct
  tagTwoInterfInspiralEvent
  {
    UINT4                                  twoInterfId;
    REAL4                                  snrsq;
    InspiralTemplate                       tmplt;
    InspiralEvent                         *eventIn1;
    InspiralEvent                         *eventIn2;
    struct tagTwoInterfInspiralEvent      *next;
  }
  TwoInterfInspiralEvent;
  
/********************************************************** <lalLaTeX>

\subsubsection*{\texttt{struct TwoInterfFindChirpInitParams}}
\idx[Type]{TwoInterfFindChirpInitParams}

\noindent Contains the parameters needed to create the filter-input structure
\texttt{LALTwoInterfFindChirpFilterInput()}. 

\begin{description}
\item[\texttt{UINT4  numSegments}]
\item[\texttt{UINT4  numPoints}]
\item[\texttt{FindChirpInitParams   *initParams1}]
\item[\texttt{FindChirpInitParams   *initParams2}]
\item[\texttt{BOOLEAN        createTwoInterfRhosqVec}]
\end{description}

*********************************************************** </lalLaTeX> */

  /*
   *
   * typedefs of parameter structures used by functions in twointerffindchirp
   *
   */
  
  /* parameter structure for all init functions */
  typedef struct
  tagTwoInterfFindChirpInitParams
  {
    UINT4                         numSegments;
    UINT4                         numPoints;  
    FindChirpInitParams          *initParams1;
    FindChirpInitParams          *initParams2;
    BOOLEAN                       createTwoInterfRhosqVec;
  }
  TwoInterfFindChirpInitParams;
  
/********************************************************** <lalLaTeX>

\subsubsection*{\texttt{struct TwoInterfFindChirpFilterParams}}
\idx[Type]{TwoInterfFindChirpFilterParams}

\noindent Contains the parameters needed for a chirp search by the 
function \texttt{LALTwoInterfFindChirpFilterSegment()}. 

\begin{description}
\item[\texttt{REAL4                         maxLag}]
\item[\texttt{REAL4                         twoInterfRhosqThresh}]
\item[\texttt{REAL4Vector                  *twoInterfRhosqVec}]
\item[\texttt{FindChirpFilterParams        *filterParams1}]
\item[\texttt{FindChirpFilterParams        *filterParams2}]
\end{description}

*********************************************************** </lalLaTeX> */

  /* parameter structure for the filtering function */
  typedef struct
  tagTwoInterfFindChirpFilterParams
  {
    REAL4                         maxLag; 
    REAL4                         twoInterfRhosqThresh;
    REAL4Vector                  *twoInterfRhosqVec;
    FindChirpFilterParams        *filterParams1;
    FindChirpFilterParams        *filterParams2;  
  }
  TwoInterfFindChirpFilterParams;
  
/********************************************************** <lalLaTeX>

\subsubsection*{\texttt{struct TwoInterfFindChirpFilterInput}}
\idx[Type]{TwoInterfFindChirpFilterInput}

\noindent Contains the input structures needed for a chirp 
search by the function \texttt{LALTwoInterfFindChirpFilterSegment()}. 

\begin{description}
\item[\texttt{InspiralTemplate                      *tmplt}]
\item[\texttt{FindChirpTemplate                     *fcTmplt}]
\item[\texttt{LALDetectorPair                       *detectors}]
\item[\texttt{FindChirpSegment                      *segment1}]
\item[\texttt{FindChirpSegment                      *segment2}]
\end{description}

*********************************************************** </lalLaTeX> */
  
  /*
   *
   * typedefs of input structures used by functions in twointerffindchirp
   *
   */
  
  
  /* input to the filtering functions */
  typedef struct
  tagTwoInterfFindChirpFilterInput
  {
    InspiralTemplate                      *tmplt;
    FindChirpTemplate                     *fcTmplt;
    LALDetectorPair                       *detectors;
    FindChirpSegment                      *segment1;
    FindChirpSegment                      *segment2;
  }
  TwoInterfFindChirpFilterInput;

  /*
   *
   * function prototypes for initialization and finalization functions
   *
   */
  
  
  void
  LALTwoInterfFindChirpFilterInit (
				   LALStatus                           *status,
				   TwoInterfFindChirpFilterParams     **output,
				   TwoInterfFindChirpInitParams        *params
				   );
  
  void
  LALTwoInterfFindChirpFilterFinalize (
				       LALStatus                       *status,
				       TwoInterfFindChirpFilterParams **output
				       );
  
  void
  LALCreateTwoInterfFindChirpFilterInput (
				    LALStatus                          *status,
				    TwoInterfFindChirpFilterInput     **output,
				    TwoInterfFindChirpInitParams       *params
				    );
  
  void
  LALDestroyTwoInterfFindChirpInput (
				     LALStatus                         *status,
				     TwoInterfFindChirpFilterInput    **output
				     );
  
  
  /*
   *
   * function prototype for the filtering function
   *
   */
  
  
  void
  LALTwoInterfFindChirpFilterSegment (
				      LALStatus                      *status,
				      TwoInterfInspiralEvent        **eventList,
				      TwoInterfFindChirpFilterInput  *input,
				      TwoInterfFindChirpFilterParams *params
				      );
#ifdef  __cplusplus
}
#endif /* C++ protection */

/********************************************************** <lalLaTeX>
							    
\vfill{\footnotesize\input{TwoInterfFindChirpHV}}
							   
\newpage\input{TwoInterfFindChirpFilterC}
%\newpage\input{TwoInterfFindChirpTestC}
*********************************************************** </lalLaTeX> */

