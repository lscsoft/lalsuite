/************************************ <lalVerbatim file="LALErrorCV">
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{LALError.c}}
\label{ss:LALError.c}

Error handling routines for LAL.  These should \emph{not} be invoked in
production code, except in very specific circumstances.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALErrorCP}
\index{\verb&LALPrintError()&}
\index{\verb&LALAbort()&}

\subsubsection*{Description}

These functions cause LAL to print error messages and abort.  Their
implementation is quite simple but may be altered in the future to provide
reasonable behaviour when integrated with other systems (e.g., LDAS).  Use of
these routines is restricted to a few isolated places within LAL.

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALErrorCV}}

</lalLaTeX> */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "LALError.h"

NRCSID( LALERRORC, "$Id$" );

/* <lalVerbatim file="LALErrorCP"> */
int
LALPrintError( const char *fmt, ... )
{ /* </lalVerbatim> */
  int n;
  va_list ap;
  va_start( ap, fmt );
  n = vfprintf( stderr, fmt, ap );
  va_end( ap );
  return n;
}


/* <lalVerbatim file="LALErrorCP"> */
void
LALAbort( const char *msg )
{ /* </lalVerbatim> */
  LALPrintError( msg );
  abort();
}
