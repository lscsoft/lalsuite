/* <lalVerbatim file="UtoGPSCV">
Author: David Chin <dwchin@umich.edu> +1-734-730-1274
$Id$
</lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{UtoGPS.c}}
\label{ss:UtoGPS.c}

Converts time between Unix and GPS epochs.


\subsection*{Prototypes}
\vspace{0.1in}
\input{UtoGPSCP}
\index{\texttt{LALUtoGPS()}}
\index{\texttt{LALGPtoU()}}

\subsubsection*{Description}

\texttt{UtoGPS()} converts time measured in the Unix epoch to GPS
time.  \texttt{GPStoU()} converts GPS time to Unix time. (See
\texttt{ctime(3C)} in the Unix manual.)

\subsubsection*{Algorithms}

The only difference between time measured in the Unix epoch and time
measured in the GPS epoch is the choice of reference, \textit{i.e.}
the time from which an elapsed time is measured.  The Unix reference
is 1970-01-01 00:00:00 UTC, and the GPS reference is 1980-01-06
00:00:00 UTC.  GPS is \emph{not} adjusted for leap seconds.  As of
1999-01-01, 

\begin{tabular}{rl}
  TAI is ahead of UTC by & 32 seconds \\
  TAI is ahead of GPS by & 19 seconds \\
  GPS is ahead of UTC by & 13 seconds
\end{tabular}

For more information, consult the US Naval Observatory
(\url{http://tycho.usno.navy.mil/leapsec.html}), and also the
comments in \texttt{Utime.c}.

\subsubsection*{Uses}

\subsubsection*{Notes}

For more information, see \url{http://tycho.usno.navy.mil/leapsec.html}
and also the comments in \texttt{Utime.c}.

The integer number of seconds between the Unix epoch and the GPS epoch is
315964811: 8 years, 2 leap years, 5 days, and $11 = (19.0 - 8.0)$ leap
seconds between 1970-01-01 00:00:00 and 1980-01-06 00:00:00.  So,
converting from GPS to Unix epoch and vice versa requires adding or
subtracting 315964811 seconds. Strictly speaking, however, the amount of
time between the Unix epoch and the GPS epoch is 315964810.999918 seconds.


</lalLaTeX> */

#include <lal/LALRCSID.h>

NRCSID (UTOGPSC, "$Id$");

#include <lal/Date.h>
#include "date_value.h"


/*
 * Convert UTC seconds to GPS seconds
 *
 * Input:
 *
 * Output:
 */
/* <lalVerbatim file="UtoGPSCP"> */ 
void
LALUtoGPS (LALStatus          *status,
           LIGOTimeGPS        *gpstime,
           const LIGOTimeUnix *unixtime)
{ /* </lalVerbatim> */
  INITSTATUS (status, "LALUtoGPS", UTOGPSC);

  /*
   * Check pointer to input variable
   */
  ASSERT (unixtime != (LIGOTimeUnix *)NULL, status,
          DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  /*
   * Check pointer to output variable:
   * allocate memory if necessary
   */
  ASSERT (gpstime != (LIGOTimeGPS *)NULL, status,
          DATEH_ENULLOUTPUT, DATEH_MSGENULLOUTPUT);


  /* shift the time to GPS time
   * GPS epoch is 1980-01-06, Unix epoch is 1970-01-01, so GPS time must
   * be smaller than UTC time */
  gpstime->gpsSeconds = unixtime->unixSeconds - UNIXGPS;

  RETURN (status);
} /* END LALUtoGPS() */


/*
 * Convert GPS seconds to UTC seconds
 *
 * Input:
 *
 * Output:
 */
/* <lalVerbatim file="UtoGPSCP"> */ 
void
LALGPStoU (LALStatus         *status,
           LIGOTimeUnix      *unixtime,
           const LIGOTimeGPS *gpstime)
{ /* </lalVerbatim> */
  INITSTATUS (status, "LALGPStoU", UTOGPSC);

  /*
   * Check pointer to input variable
   */
  ASSERT (gpstime != (LIGOTimeGPS *)NULL, status,
          DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  /*
   * Check pointer to output variable:
   * allocate memory if necessary
   */
  ASSERT (unixtime != (LIGOTimeUnix *)NULL, status,
          DATEH_ENULLOUTPUT, DATEH_MSGENULLOUTPUT);

  /* shift the time to UTC time
   * GPS epoch is 1980-01-06, Unix epoch is 1970-01-01, so GPS time must
   * be smaller than UTC time */
  unixtime->unixSeconds = gpstime->gpsSeconds + UNIXGPS;

  RETURN (status);
} /* END LALGPStoU() */


