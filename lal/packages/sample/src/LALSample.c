/************************************ <lalVerbatim file="LALSampleCV">
Author: Creighton, T. D.
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{LALSample.c}}
\label{ss:LALSample.c}

Example module for LAL.

\subsubsection*{Prototypes}
\input{LALSampleCP}
\index{\texttt{REAL8Invert()}}
\index{\texttt{REAL8Divide()}}

\subsubsection*{Description}

This module exists to demostrate documentation and coding standards
for LAL modules, using two trivial functions.  \verb@REAL8Invert()@
computes \verb@*output@ = 1/\verb@input@; if \verb@input@ = 0, it
leaves \verb@*output@ unchanged and returns an error.
\verb@REAL8Divide()@ computes \verb@*output@ =
\verb@numer@/\verb@denom@, calling \verb@REAL8Invert()@ as a
subroutine.  This allows us to demonstrate LAL error handling through
nested function calls.

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALSampleCV}}

******************************************************* </lalLaTeX> */

#include "LALStdlib.h"
#include "LALSample.h"

NRCSID( LALSAMPLEC, "$Id$" );

/* <lalVerbatim file="LALSampleCP"> */
void
REAL8Invert( Status *stat, REAL8 *output, REAL8 input )
{ /* </lalVerbatim> */
  INITSTATUS( stat, "REAL8Invert", LALSAMPLEC );

  /* This traps coding errors in the calling routine. */
  ASSERT( output != NULL, stat, LALSAMPLEH_ENULL, LALSAMPLEH_MSGENULL );

  /* This traps runtime errors. */
  if ( input == 0.0 )
    ABORT( stat, LALSAMPLEH_EDIV0, LALSAMPLEH_MSGEDIV0 );

  *output = 1.0/input;
  RETURN( stat );
}


/* <lalVerbatim file="LALSampleCP"> */
void
REAL8Divide( Status *stat, REAL8 *output, REAL8 numer, REAL8 denom )
{ /* </lalVerbatim> */
  INITSTATUS( stat, "REAL8Divide", LALSAMPLEC );
  ATTATCHSTATUSPTR( stat );

  TRY( REAL8Invert( stat->statusPtr, output, denom ), stat );
  *output *= numer;

  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
