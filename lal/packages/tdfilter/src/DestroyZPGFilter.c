/*
*  Copyright (C) 2007 Jolien Creighton, Teviet Creighton
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

/***************************** <lalVerbatim file="DestroyZPGFilterCV">
Author: Creighton, T. D.
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{DestroyZPGFilter.c}}
\label{ss:DestroyZPGFilter.c}

Destroys ZPG filter objects.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{DestroyZPGFilterCP}
\idx{LALDestroyCOMPLEX8ZPGFilter()}
\idx{LALDestroyCOMPLEX16ZPGFilter()}

\subsubsection*{Description}

These functions destroy an object \verb@**output@ of type
\verb@COMPLEX8ZPGFilter@ or \verb@COMPLEX16ZPGFilter@, and set
\verb@*output@ to \verb@NULL@.

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
LALFree()
LALCDestroyVector()
LALZDestroyVector()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{DestroyZPGFilterCV}}

******************************************************* </lalLaTeX> */

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/ZPGFilter.h>

void XLALDestroyCOMPLEX8ZPGFilter( COMPLEX8ZPGFilter *filter )
{
  if ( filter )
  {
    XLALDestroyCOMPLEX8Vector( filter->zeros );
    XLALDestroyCOMPLEX8Vector( filter->poles );
    LALFree( filter );
  }
  return;
}

void XLALDestroyCOMPLEX16ZPGFilter( COMPLEX16ZPGFilter *filter )
{
  if ( filter )
  {
    XLALDestroyCOMPLEX16Vector( filter->zeros );
    XLALDestroyCOMPLEX16Vector( filter->poles );
    LALFree( filter );
  }
  return;
}

/* <lalVerbatim file="DestroyZPGFilterCP"> */
void
LALDestroyCOMPLEX8ZPGFilter( LALStatus         *stat,
			     COMPLEX8ZPGFilter **input )
{ /* </lalVerbatim> */
  INITSTATUS(stat);

  /* Make sure handle is non-null, and points to a non-null
     pointer. */
  ASSERT(input,stat,ZPGFILTERH_ENUL,ZPGFILTERH_MSGENUL);
  ASSERT(*input,stat,ZPGFILTERH_ENUL,ZPGFILTERH_MSGENUL);

  XLALDestroyCOMPLEX8ZPGFilter( *input );
  *input = NULL;

  /* Normal exit */
  RETURN(stat);
}


/* <lalVerbatim file="DestroyZPGFilterCP"> */
void
LALDestroyCOMPLEX16ZPGFilter( LALStatus          *stat,
			      COMPLEX16ZPGFilter **input )
{ /* </lalVerbatim> */
  INITSTATUS(stat);

  /* Make sure handle is non-null, and points to a non-null
     pointer. */
  ASSERT(input,stat,ZPGFILTERH_ENUL,ZPGFILTERH_MSGENUL);
  ASSERT(*input,stat,ZPGFILTERH_ENUL,ZPGFILTERH_MSGENUL);

  XLALDestroyCOMPLEX16ZPGFilter( *input );
  *input = NULL;

  /* Normal exit */
  RETURN(stat);
}
