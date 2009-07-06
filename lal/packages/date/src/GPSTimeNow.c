/*
 * $Id$
 *
 * Copyright (C) 2007  Brown, D. A., and Kipp Cannon
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
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

#include <time.h>
#include <lal/Date.h>
#include <lal/XLALError.h>

NRCSID( GPSTIMENOWC, "$Id$" );

/* <lalVerbatim file="GPSTimeNowCP"> */
LIGOTimeGPS *
XLALGPSTimeNow (
    LIGOTimeGPS *gpstime
    )
/* </lalVerbatim> */
{
  static const char func[] = "XLALGPSTimeNow";
  time_t ticks = time(NULL);

  gpstime->gpsSeconds = XLALUTCToGPS(gmtime(&ticks));
  gpstime->gpsNanoSeconds = 0;

  /* XLALUTCToGPS returns < 0 on error, even though of course time did not
   * begin at GPS 0 */
  if(gpstime->gpsSeconds < 0)
    XLAL_ERROR_NULL(func, XLAL_EFUNC);

  return gpstime;
}


/* <lalVerbatim file="GPSTimeNowCP"> */
void
LALGPSTimeNow (
    LALStatus           *status,
    LIGOTimeGPS         *gpstime,
    const LALLeapSecAccuracy  *accuracy
    )
/* </lalVerbatim> */
{
  INITSTATUS( status, "LALGPSTimeNow", GPSTIMENOWC );
  ATTATCHSTATUSPTR( status );
  ASSERT( gpstime, status, DATEH_ENULLOUTPUT, DATEH_MSGENULLOUTPUT );
  ASSERT( accuracy, status, DATEH_ENULLINPUT, DATEH_MSGENULLINPUT );

  ASSERT( XLALGPSTimeNow(gpstime), status, DATEH_EDATETOOEARLY, DATEH_MSGEDATETOOEARLY );

  XLALPrintDeprecationWarning("LALGPSTimeNow", "XLALGPSTimeNow");

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
