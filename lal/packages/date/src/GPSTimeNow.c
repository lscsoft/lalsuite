/*
 * File Name: GPSTimeNow.c
 *
 * Author: Brown, D. A.
 *
 * Revision: $Id$
 *
 */

#if 0
<lalVerbatim file="GPSTimeNowCV">
Author: Duncan Brown
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{GPSTimeNow.c}}
\label{ss:GPSTimeNow.c}

Routine to convert the current unix system clock time into a
\texttt{LIGOTimeGPS} structure.

\subsection*{Prototypes}
\vspace{0.1in}
\input{GPSTimeNowCP}
\idx{LALGPSTimeNow()}

\subsubsection*{Description}

This module contains a single funtion that converts the current unix system
time as returned by the \texttt{time()} function to GPS seconds. The leap
second accuracy is determined by the \texttt{accuracy} argument.

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}

If the system clock time is incorect, the returned GPS time will obviously be
wrong by the same amount.

\vfill{\footnotesize\input{GPSTimeNowCV}}
</lalLaTeX>
#endif

#include <lal/LALStdlib.h>
#include <lal/Date.h>

NRCSID( GPSTIMENOWC, "$Id$" );

/* <lalVerbatim file="GPSTimeNowCP"> */
void
LALGPSTimeNow (
    LALStatus           *status,
    LIGOTimeGPS         *gpstime,
    const LALLeapSecAccuracy  *accuracy
    )
/* </lalVerbatim> */
{
  time_t ticks;
  LALDate laldate;
  INITSTATUS( status, "LALGPSTimeNow", GPSTIMENOWC );
  ATTATCHSTATUSPTR( status );
  ASSERT( gpstime, status, DATEH_ENULLOUTPUT, DATEH_MSGENULLOUTPUT );
  ASSERT( accuracy, status, DATEH_ENULLINPUT, DATEH_MSGENULLINPUT );
  
  ticks = time( NULL );
  gmtime_r( &ticks, &(laldate.unixDate) );
  laldate.residualNanoSeconds = 0;
  LALUTCtoGPS( status->statusPtr, gpstime, &laldate, accuracy );
  CHECKSTATUSPTR( status );
  
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
