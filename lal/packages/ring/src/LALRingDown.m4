/************************************ <lalVerbatim file="LALRingDownCV">
Author: Tibbits, M M
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{LALRingDown.c}}
\label{s:LALRingDown.c}

Routine to fill a frequency series with given parameters.

\subsubsection*{Prototypes}
\input{LALRingDownCP}

\subsubsection*{Description}

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALRingDownCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include "LALRingDown.h"

NRCSID( LALRINGDOWNC, "$Id$");

define(`TYPECODE',`D')
include(`LALRingDownBase.m4')

define(`TYPECODE',`S')
include(`LALRingDownBase.m4')
