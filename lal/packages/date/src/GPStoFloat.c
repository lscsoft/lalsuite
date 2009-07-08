/*
*  Copyright (C) 2007 David Chin, Jolien Creighton, Kipp Cannon, Reinhard Prix
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

/* <lalVerbatim file="GPStoFloatCV">

Author: Berukoff, S.J.  <steveb@aei-potsdam.mpg.de>
$Id$

</lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{GPStoFloat.c}}
\label{ss:GPStoFloat.c}

Converts between \texttt{LALTimeInterval} and \texttt{REAL8} formats.


\subsection*{Prototypes}
\vspace{0.1in}
\input{GPStoFloatCP}
\idx{LALFloatToInterval()}
\idx{LALIntervalToFloat()}

\subsubsection*{Description}

\begin{itemize}
  \item \texttt{LALFloatToInterval()} converts a time interval in \texttt{REAL8} to a \texttt{LALTimeInterval}
  \item \texttt{LALIntervalToFloat()} converts a time interval in \texttt{LALTimeInterval} to \texttt{REAL8}
\end{itemize}

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{GPStoFloatCV}}

</lalLaTeX> */

#include <math.h>

#include "Date.h"

NRCSID(GPSTOFLOATC, "$Id$");

/* D'oh. rint(3) is not in the ANSI std. Thanks Jolien. */
#define rint(x) floor(0.5 + (x))

static const INT4 oneBillion = 1000000000;

/* <lalVerbatim file="GPStoFloatCP"> */
void LALFloatToInterval(LALStatus *status,
                        LALTimeInterval *pInterval,  /* output: deltaT in LALTimeInterval format */
                        const REAL8 *pDeltaT) /* input: time interval in floating point */
{  /* </lalVerbatim> */
  INITSTATUS(status, "LALFloatToInterval", GPSTOFLOATC);
  ATTATCHSTATUSPTR(status);

  ASSERT(pInterval, status, DATEH_ENULLOUTPUT, DATEH_MSGENULLOUTPUT);
  ASSERT(pDeltaT, status, DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  pInterval->seconds     = (INT4)(*pDeltaT);
  pInterval->nanoSeconds = (INT4)rint((*pDeltaT -
                                       (REAL8)(pInterval->seconds)) *
                                      (REAL8)oneBillion);


  DETATCHSTATUSPTR(status);
  RETURN(status);
} /* END: FloatToInterval() */


/* <lalVerbatim file="GPStoFloatCP"> */
void LALIntervalToFloat(LALStatus *status,
                        REAL8 *pDeltaT,    /* output: floating point time interval */
                        const LALTimeInterval *pInterval) /* input: LALTimeInterval */
{ /* </lalVerbatim> */
  INITSTATUS(status, "LALIntervalToFloat", GPSTOFLOATC);
  ATTATCHSTATUSPTR(status);

  ASSERT(pDeltaT, status, DATEH_ENULLOUTPUT, DATEH_MSGENULLOUTPUT);
  ASSERT(pInterval, status, DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  *pDeltaT = (REAL8)(pInterval->seconds) +
    ((REAL8)(pInterval->nanoSeconds) / (REAL8)oneBillion);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}
