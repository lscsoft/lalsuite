/***************************** <lalVerbatim file="GenerateRingHV">
Author: Brady, P R
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\section{Header \texttt{GenerateRing.h}}
\label{s:GenerateRing.h}

Provides routines to generate ringdown waveforms.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/GenerateRing.h>
\end{verbatim}

This header covers routines to generate a variety of ringdown waveforms
......


******************************************************* </lalLaTeX> */

#ifndef _GENERATERING_H
#define _GENERATERING_H

#include <lal/LALStdlib.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/SkyCoordinates.h>
#include <lal/LIGOMetadataTables.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif


/********************************************************** <lalLaTeX>
\subsection*{Error conditions}
****************************************** </lalLaTeX><lalErrTable> */
#define GENERATERINGH_ENUL 1
#define GENERATERINGH_EOUT 2
#define GENERATERINGH_EMEM 3
#define GENERATERINGH_ETYP 4
#define GENERATERINGH_ELEN 5

#define GENERATERINGH_MSGENUL "Unexpected null pointer in arguments"
#define GENERATERINGH_MSGEOUT "Output field a, f, phi, or shift already exists"
#define GENERATERINGH_MSGEMEM "Out of memory"
#define GENERATERINGH_MSGETYP "Waveform type not implemented"
#define GENERATERINGH_MSGELEN "Waveform length not correctly specified"
/******************************************** </lalErrTable><lalLaTeX>

\subsection*{Types}
\idx[Type]{SimRingType}
\idx[Type]{RingParamStruc}

\subsubsection*{Structure \texttt{SimRingType}}

******************************************************* </lalLaTeX> */

/******************************************** <lalVerbatim> */
typedef enum{
  Ringdown
} SimRingType;
/******************************************** </lalVerbatim> */


/******************************************** <lalLaTeX>

\subsubsection*{Structure \texttt{RingParamStruc}}

This structure stores the parameters for constructing a burst gravitational
waveform 

******************************************************* </lalLaTeX> */

/******************************************** <lalVerbatim> */
typedef struct tagRingParamStruc {
  REAL8 deltaT;             /* requested sampling interval (s) */
  CoordinateSystem system;  /* coordinate system to assume for simRing */
} RingParamStruc;
/******************************************** </lalVerbatim> */


/* <lalLaTeX>
\vfill{\footnotesize\input{GenerateRingHV}}
</lalLaTeX> */


/* Function prototypes. */

/* <lalLaTeX>
\newpage\input{GenerateRingC}
</lalLaTeX> */
void
LALGenerateRing( 
    LALStatus          *status, 
    CoherentGW         *output, 
    SimRingTable       *simRing,
    RingParamStruc     *params
    );

void
LALRingInjectSignals( 
    LALStatus               *status, 
    REAL4TimeSeries         *series, 
    SimRingTable           *injections,
    COMPLEX8FrequencySeries *resp,
    INT4                     calType
    );

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _GENERATERING_H */
