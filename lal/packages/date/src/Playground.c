/*----------------------------------------------------------------------- 
 * 
 * File Name: Playground.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="PlaygroundCV">
Author: Brown D. A.
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{Playground.c}}
\label{ss:Playground.c}

Determines if a given time (or segment) is playground data.

\subsection*{Prototypes}
\vspace{0.1in}
\input{PlaygroundCP}
\idx{LALINT8NanoSecIsPlayground()}
\idx{LALGPSIsPlayground()}
\idx{LALSegmentIsPlayground()}

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

#include <lal/LALStdlib.h>
#include <lal/Date.h>

NRCSID( PLAYGROUNDC, "$Id$" );

/* <lalVerbatim file="PlaygroundCP"> */
void
LALINT8NanoSecIsPlayground (
    LALStatus          *status,
    INT4               *playground,
    INT8               *ns
    )
/* </lalVerbatim> */
{
  const INT8 start = 729273613 * LAL_INT8_C(1000000000);
  const INT8 interval = 6370 * LAL_INT8_C(1000000000);
  const INT8 length = 600 * LAL_INT8_C(1000000000);

  INITSTATUS( status, "LALINT8NanoSecIsPlayground", PLAYGROUNDC );

  if ( (*ns - start) % interval < length )
  {
    *playground = 1;
  }
  else
  {
    *playground = 0;
  }

  RETURN( status );
}

/* <lalVerbatim file="PlaygroundCP"> */
void
LALGPSIsPlayground (
    LALStatus          *status,
    INT4               *playground,
    LIGOTimeGPS        *gpstime
    )
/* </lalVerbatim> */
{
  INT8  ns;

  INITSTATUS( status, "LALINT8NanoSecIsPlayground", PLAYGROUNDC );
  ATTATCHSTATUSPTR( status );

  LALGPStoINT8( status->statusPtr, &ns, gpstime );
  CHECKSTATUSPTR( status );

  LALINT8NanoSecIsPlayground( status->statusPtr, playground, &ns );
  CHECKSTATUSPTR( status );
  
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

/* <lalVerbatim file="PlaygroundCP"> */
void
LALSegmentIsPlayground (
    LALStatus          *status,
    INT4               *playground,
    LIGOTimeGPS        *gpsStart,
    LIGOTimeGPS        *gpsEnd
    )
/* </lalVerbatim> */
{
  INITSTATUS( status, "LALSegmentIsPlayground", PLAYGROUNDC );
  ATTATCHSTATUSPTR( status );

  /* JC: THIS FUNCTION DOESN'T SEEM TO DO MUCH... I GUESS IT SHOULD DO SOMETHING*/
  playground=NULL;
  gpsStart=NULL;
  gpsEnd=NULL;

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
