/*
*  Copyright (C) 2007 Jolien Creighton, Teviet Creighton, John Whelan
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

/********************************* <lalVerbatim file="InjectVectorCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{InjectVector.c}}
\label{ss:InjectVector.c}

Injects a vector of floating-point numbers into a vector of integers,
with dithering.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{InjectVectorCP}
\idx{LALSI2InjectVector()}
\idx{LALSSInjectVector()}

\subsubsection*{Description}

The function \verb@LALSI2InjectVector()@ (i.e.\ ``Single-precision to
\verb@INT2@'') dithers the contents of \verb@*output@, adds the
contents of \verb@*signalvec@, and rounds to the nearest integer, storing
the result back in \verb@*output@.  If desired, the random parameters
for the dithering can be created outside this routine and passed in as
\verb@*params@ (see \verb@Random.h@); if this pointer is \verb@NULL@,
the parameters will be generated internally.

The function \verb@LALSSInjectVector()@ (i.e.\ ``Single-precision to
single-precision'') simply adds the contents of \verb@*signalvec@ to
\verb@*output@ where they overlap, without performing any dithering.

\subsubsection*{Algorithm}

Dithering is done with a flat random distribution as described in
\verb@Inject.h@.  Injected values outside the dynamic range of the
output force the output to its ``rails'' of $-2^{8N-1}$ or
$2^{8N-1}-1$, where $N$ is the number of bytes in the integer.  The
two vectors need not be of equal length; the injection stops when
either vector reaches its end.

If \verb@params@ is \verb@NULL@, a \verb@RandomParams@ structure will
be generated internally using a seed of zero (i.e.\ the current time
will be used to initialize the pseudorandom sequence).

\subsubsection*{Uses}
\begin{verbatim}
LALCreateRandomParams()
LALDestroyRandomParams()
LALUniformDeviate()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{InjectVectorCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/Random.h>
#include <lal/Inject.h>

NRCSID( INJECTVECTORC, "$Id$" );

/* <lalVerbatim file="InjectVectorCP"> */
void
LALSI2InjectVector( LALStatus    *stat,
		    INT2Vector   *output,
		    REAL4Vector  *signalvec,
		    RandomParams *params )
{ /* </lalVerbatim> */
  UINT4 n;  /* number of samples injected */
  UINT4 i;  /* an index */
  RandomParams *internal = NULL; /* internal random parameters */

  const INT2 max = 32767;  /* largest INT2 */
  const INT2 min = -32768; /* smallest INT2 */

  INITSTATUS( stat, "LALSI2InjectVector", INJECTVECTORC );
  ATTATCHSTATUSPTR( stat );

  /* Make sure parameter structures and their fields exist. */
  ASSERT( signalvec, stat, INJECTH_ENUL, INJECTH_MSGENUL );
  ASSERT( signalvec->data, stat, INJECTH_ENUL, INJECTH_MSGENUL );
  ASSERT( output, stat, INJECTH_ENUL, INJECTH_MSGENUL );
  ASSERT( output->data, stat, INJECTH_ENUL, INJECTH_MSGENUL );

  /* If params = NULL, generate an internal set of parameters. */
  if ( !params )
    TRY( LALCreateRandomParams( stat->statusPtr, &internal, 0 ),
	 stat );
  else
    internal = params;

  /* Find out how many samples will be injected. */
  n = output->length;
  if ( n > signalvec->length )
    n = signalvec->length;

  /* Start injecting... */
  for ( i = 0; i < n; i++ ) {
    REAL4 x = (REAL4)( output->data[i] ); /* current output sample */
    REAL4 d;                              /* current dithering */

    /* Compute the dithering. */
    LALUniformDeviate( stat->statusPtr, &d, internal );
    BEGINFAIL( stat )
      if ( !params ) {
	TRY( LALDestroyRandomParams( stat->statusPtr, &internal ),
	     stat );
      }
    ENDFAIL( stat );

    /* Dither and inject. */
    x += d + signalvec->data[i];
    if ( x > max )
      output->data[i] = max;
    else if ( x < min )
      output->data[i] = min;
    else
      output->data[i] = (INT2)( floor( x ) );
  }

  /* Cleanup and exit. */
  if ( !params ) {
    TRY( LALDestroyRandomParams( stat->statusPtr, &internal ), stat );
  }
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


/* <lalVerbatim file="InjectVectorCP"> */
void
LALSSInjectVector( LALStatus    *stat,
		   REAL4Vector  *output,
		   REAL4Vector  *signalvec )
{ /* </lalVerbatim> */
  UINT4 n;  /* number of samples injected */
  UINT4 i;  /* an index */

  INITSTATUS( stat, "LALSSInjectVector", INJECTVECTORC );

  /* Make sure parameter structures and their fields exist. */
  ASSERT( signalvec, stat, INJECTH_ENUL, INJECTH_MSGENUL );
  ASSERT( signalvec->data, stat, INJECTH_ENUL, INJECTH_MSGENUL );
  ASSERT( output, stat, INJECTH_ENUL, INJECTH_MSGENUL );
  ASSERT( output->data, stat, INJECTH_ENUL, INJECTH_MSGENUL );

  /* Find out how many samples will be injected. */
  n = output->length;
  if ( n > signalvec->length )
    n = signalvec->length;

  /* Inject and exit. */
  for ( i = 0; i < n; i++ )
    output->data[i] += signalvec->data[i];
  RETURN( stat );
}
