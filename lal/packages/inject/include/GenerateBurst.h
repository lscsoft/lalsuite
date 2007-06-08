/*
*  Copyright (C) 2007 Jolien Creighton, Patrick Brady, Saikat Ray-Majumder
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
#define GENERATEBURSTH_ETYP 4
#define GENERATEBURSTH_ELEN 5

#define GENERATEBURSTH_MSGENUL "Unexpected null pointer in arguments"
#define GENERATEBURSTH_MSGEOUT "Output field a, f, phi, or shift already exists"
#define GENERATEBURSTH_MSGEMEM "Out of memory"
#define GENERATEBURSTH_MSGETYP "Waveform type not implemented"
#define GENERATEBURSTH_MSGELEN "Waveform length not correctly specified"
/******************************************** </lalErrTable><lalLaTeX>

\subsection*{Types}
\idx[Type]{SimBurstType}
\idx[Type]{BurstParamStruc}

\subsubsection*{Structure \texttt{SimBurstType}}

******************************************************* </lalLaTeX> */

/******************************************** <lalVerbatim> */
typedef enum{
  sineGaussian,
  Gaussian,
  Ringdown,
  Ringup
} SimBurstType;
/******************************************** </lalVerbatim> */


/******************************************** <lalLaTeX>

\subsubsection*{Structure \texttt{BurstParamStruc}}

This structure stores the parameters for constructing a burst gravitational
waveform 

******************************************************* </lalLaTeX> */

/******************************************** <lalVerbatim> */
typedef struct tagBurstParamStruc {
  REAL8 deltaT;             /* requested sampling interval (s) */
  CoordinateSystem system;  /* coordinate system to assume for simBurst */
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
    LALStatus          *status, 
    CoherentGW         *output, 
    SimBurstTable      *simBurst,
    BurstParamStruc    *params 
    );

void
LALBurstInjectSignals( 
    LALStatus               *status, 
    REAL4TimeSeries         *series, 
    SimBurstTable           *injections,
    COMPLEX8FrequencySeries *resp,
    INT4                     calType
    );

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _GENERATEBURST_H */
