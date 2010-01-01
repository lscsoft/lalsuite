/*
*  Copyright (C) 2007 Jolien Creighton
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

/************************* <lalVerbatim file="DestroyResampleRulesCV">
Author: Creighton, T. D.
Revision: $Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{DestroyResampleRules.c}}
\label{ss:DestroyResampleRules.c}

Destroys an object of type \verb@ResampleRules@.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{DestroyResampleRulesCP}
\idx{LALDestroyResampleRules()}

\subsubsection*{Description}

This function destroys an object \verb@**rules@ of type
\texttt{ResampleRules}, and sets \verb@*rules@ to \verb@NULL@.

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
void LALFree()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{DestroyResampleRulesCV}}

******************************************************* </lalLaTeX> */

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/Resample.h>

NRCSID(DESTROYRESAMPLERULESC,"$Id$");


/* <lalVerbatim file="DestroyResampleRulesCP"> */
void
LALDestroyResampleRules( LALStatus     *stat,
			 ResampleRules **rules )
{ /* </lalVerbatim> */
  INITSTATUS(stat,"LALDestroyResampleRules",DESTROYRESAMPLERULESC);

  /* Make sure that the handle is non-null, that it points to a
     non-null pointer, and that the interval and shift fields are
     non-null.) */
  ASSERT(rules,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT(*rules,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT((*rules)->interval,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT((*rules)->shift,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);

  /* Free all the memory. */
  LALFree((*rules)->interval);
  LALFree((*rules)->shift);
  LALFree(*rules);

  /* Point the handle to NULL, then get out of here. */
  *rules=NULL;
  RETURN(stat);
}
