/***************************** <lalVerbatim file="GenerateBurstHV">
Author: Brady, P R
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\section{Header \texttt{GenerateBurst.h}}
\label{s:GenerateBurst.h}

Provides routines to generate burst waveforms.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/GenerateBurst.h>
\end{verbatim}

This header covers routines to generate a variety of burst waveforms
......


******************************************************* </lalLaTeX> */

#ifndef _GENERATEBURST_H
#define _GENERATEBURST_H

#include <lal/LALStdlib.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/SkyCoordinates.h>
#include <lal/LIGOMetadataTables.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif

NRCSID( GENERATEBURSTH, "$Id$" );

/********************************************************** <lalLaTeX>
\subsection*{Error conditions}
****************************************** </lalLaTeX><lalErrTable> */
#define GENERATEBURSTH_ENUL 1
#define GENERATEBURSTH_EOUT 2
#define GENERATEBURSTH_EMEM 3

#define GENERATEBURSTH_MSGENUL "Unexpected null pointer in arguments"
#define GENERATEBURSTH_MSGEOUT "Output field a, f, phi, or shift already exists"
#define GENERATEBURSTH_MSGEMEM "Out of memory"
/******************************************** </lalErrTable><lalLaTeX>

\subsection*{Types}
\idx[Type]{SimBurstType}
\idx[Type]{BurstParamStruc}

\subsubsection*{Structure \texttt{SimBurstType}}

******************************************************* </lalLaTeX> */

/******************************************** <lalVerbatim> */
typedef enum{
  sineGaussian,
  Gaussian
} SimBurstType;
/******************************************** </lalVerbatim> */


/******************************************** <lalLaTeX>

\subsubsection*{Structure \texttt{BurstParamStruc}}

This structure stores the parameters for constructing a burst gravitational
waveform 

******************************************************* </lalLaTeX> */

/******************************************** <lalVerbatim> */
typedef struct tagBurstParamStruc {
  /* Passed parameters. */
  SkyPosition position; /* location of source on sky */
  REAL4 psi;            /* polarization angle (radians) */
  LIGOTimeGPS epoch;    /* start time of output time series */

  /* Input parameters. */
  REAL8 deltaT;         /* requested sampling interval (s) */
  UINT4 length;         /* length of time series */
  REAL4 hrss;           /* root sum square amplitude */

  /* Waveform dependent parameters */
  SimBurstType burstType;  /* the type of burst to inject */
  REAL8 f0;               /* central frequency for sine-gaussian      */
  REAL8 tau;              /* decay time for (sine-)gaussian envelope  */

} BurstParamStruc;
/******************************************** </lalVerbatim> */


/* <lalLaTeX>
\vfill{\footnotesize\input{GenerateBurstHV}}
</lalLaTeX> */


/* Function prototypes. */

/* <lalLaTeX>
\newpage\input{GenerateBurstC}
</lalLaTeX> */
void
LALGenerateBurst( 
    LALStatus          *stat, 
    CoherentGW         *output, 
    BurstParamStruc    *params 
    );

/* <lalLaTeX>
\newpage\input{SimulateTaylorCWTestC}
</lalLaTeX> */

/* <lalLaTeX>
\newpage\input{BurstInjectSignalsCC}
</lalLaTeX> */
void
LALBurstInjectSignals( 
    LALStatus               *stat, 
    REAL4TimeSeries         *series, 
    SimBurstTable           *injections,
    COMPLEX8FrequencySeries *resp
    );

/* <lalLaTeX>
\newpage\input{SimulateTaylorCWTestC}
</lalLaTeX> */

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _GENERATEBURST_H */
