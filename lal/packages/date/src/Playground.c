/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Stephen Fairhurst
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

/*-----------------------------------------------------------------------
 *
 * File Name: Playground.c
 *
 * Author: Brown, D. A.
 *
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="PlaygroundCV">
Author: Brown D. A.
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{Playground.c}}
\label{ss:Playground.c}

Determines if a given time (or segment) is playground data.

\subsection*{Prototypes}
\vspace{0.1in}
\input{PlaygroundCP}
\idx{XLALINT8NanoSecIsPlayground()}

\subsubsection*{Description}

This module contains two routines to determine if a given time is in the data
designated as playground or not. The first routines takes input as
\texttt{INT8} nanoseconds and the second as a \texttt{LIGOTimeGPS} structure.
The third routine decides if some or all of a given time interval is
playground or not.

\subsubsection*{Algorithm}

The playground algorithm is given in LIGO techincal document T030020-01.
Briefly, $t$ is playground if
\begin{equation}
t - 729273613 \% 6370 < 600.
\end{equation}

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{PlaygroundCV}}
</lalLaTeX>
#endif

#include <lal/Date.h>
#include <lal/XLALError.h>

/* <lalVerbatim file="PlaygroundCP"> */
int
XLALINT8NanoSecIsPlayground (
    INT8         ns
    )
/* </lalVerbatim> */
{
  const INT8 start = 729273613 * XLAL_BILLION_INT8;
  const INT8 interval = 6370 * XLAL_BILLION_INT8;
  const INT8 length = 600 * XLAL_BILLION_INT8;

  return (ns - start) % interval < length;
}
