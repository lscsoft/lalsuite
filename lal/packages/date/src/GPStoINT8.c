/*----------------------------------------------------------------------- 
 * 
 * File Name: GPStoINT8.c
 *
 * Author: Brown, D. A., and Creighton, T. D.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="GPStoINT8CV">
Author: Brown D. A., and Creighton, T. D.
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{GPStoINT8.c}}
\label{ss:GPStoINT8.c}

Converts between \texttt{LIGOTimeGPS} and \texttt{INT8} formats.

\subsection*{Prototypes}
\vspace{0.1in}
\input{GPStoINT8CP}
\idx{LALGPStoINT8()}
\idx{LALINT8toGPS()}

\subsubsection*{Description}

This modules contains two routines, one of which converts from
\texttt{LIGOTimeGPS} to \texttt{INT8} nanoseconds, and the other, from
\texttt{INT8} nanoseconds to \texttt{LIGOTimeGPS}.

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{GPStoINT8CV}}
</lalLaTeX>
#endif

#include <lal/LALStdlib.h>

NRCSID( GPSTOINT8C, "$Id$" );

#pragma <lalVerbatim file="GPStoINT8CP">
void
LALINT8toGPS ( 
    LALStatus          *status,
    LIGOTimeGPS        *output, 
    const INT8         *input 
    )
#pragma </lalVerbatim>
{
  INT8 s = (*input) / 1000000000LL;

  INITSTATUS( status, "LALINT8toGPS", GPSTOINT8C );
  
  output->gpsSeconds = (INT4)( s );
  output->gpsNanoSeconds = (INT4)( (*input) - 1000000000LL*s );

  RETURN( status );
}

/*----------------------------------------------------------------------*/
#pragma <lalVerbatim file="GPStoINT8CP">
void
LALGPStoINT8 ( 
    LALStatus          *status,
    INT8               *output, 
    const LIGOTimeGPS  *input 
    )
#pragma </lalVerbatim>
{
  INITSTATUS( status, "LALGPStoINT8", GPSTOINT8C );
  
  *output = (INT8) input->gpsNanoSeconds 
    + 1000000000LL * (INT8) input->gpsSeconds;

  RETURN( status );
}
