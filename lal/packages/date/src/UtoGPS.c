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

Background information from
\url{http://www.iers.org/iers/earth/rotation/ut1lod/ut1lod.html}, and
\url{http://www.iers.org/iers/earth/rotation/utc/utc.html}.

International Atomic Time (TAI) began at 1958-Jan-01 00:00:00 GMT.  The TAI
clock ticks once every ISO second (with some high precision).  The origin
time was chosen such that the difference between UT1 and TAI at that time
was approximately zero. UT1 is a timescale derived from observations,
corrected for the shift in longitude of the observing station due to polar
motion. Besides being more unstable than the atomic clock, the UT1 scale is
slowing down (\textit{i.e.} the length of its second is increasing) since
the Earth's rotation is slowly decaying.

The Unix clock began at 1970-Jan-01 00:00:00 GMT. In principle, each
computer has a clock that ticks once every ISO second.  In practice,
computer clock chips are not that accurate.

Universal Co-ordinated Time (UTC) began at 1972-Jan-01 00:00:00 GMT.  The
UTC clock ticks in time with the TAI clock. At UTC's start time, TAI was
1972-Jan-01 00:00:10 GMT, i.e. TAI was ahead by 10 seconds. This difference
is due to the \textit{leap seconds}.  UTC is supposed to be the standard
for civil time, so it is meant to correspond closely to the position of the
Sun in the sky.  However, since the Earth's rotation is decaying, UTC and
UT1 are drifting apart: UTC is slightly faster.  To make UTC useful for
civil time keeping, (integer) leap seconds are added to UTC such that the
difference between UTC and UT1 is less than 0.9 seconds, \textit{e.g.}
adding one second to UTC at 2000-Dec-31 23:59:59 makeing 2000-Dec-31
23:59:60 holds the arrival of 2001 off for another second. 

The GPS clock began at 1980-Jan-06 00:00:00 UTC.  It ticks once every ISO
second because its rate is synced to that of UTC to within a few hundred
nanoseconds.  It is \em{not} adjusted for leap seconds. As of 1999-Jan-01, 

\begin{tabular}{rl}
  TAI is ahead of UTC by & 32 seconds \\
  TAI is ahead of GPS by & 19 seconds \\
  GPS is ahead of UTC by & 13 seconds
\end{tabular} 

For more information, consult the US Naval Observatory 
(\url{http://tycho.usno.navy.mil/leapsec.html}, and
\url{http://tycho.usno.navy.mil/gpstt.html}), and 
also the comments in \texttt{Utime.c}.  

The only difference between time measured in the Unix epoch and time 
measured in the GPS epoch is the choice of reference, \textit{i.e.} 
the time from which an elapsed time is measured.  The Unix reference 
is 1970-01-01 00:00:00 UTC, and the GPS reference is
1980-01-06 00:00:00 UTC.  In principle, since both Unix time and GPS time
are measures of elapsed time, the difference between the two should be
an integer that never changes with time. However, there are complications.

Most Unices do not take leap seconds into account, but some do. So, if the
clock of a Unix machine is as accurate as the TAI clock, and if they
were synchronized, and if the Unix machine did \em{not} take leap seconds
into account, then a call to \texttt{gmtime(3)} would produce a time
that was ahead.  However, if the Unix machine \em{did} take leap seconds into
account, \texttt{gmtime(3)} would return an accurate time.  If the Unix
machine were running NTP, then its clock would have leap seconds included,
and \texttt{gmtime(3)} would return an accurate time, as well.

An example will help, I think: say the Unix clock counter at
1999-Dec-31 23:59:59 reads $N$.  Calling \texttt{gmtime()} would result
in the right date.  One second later, since there was a leap second introduced,
the date would be 1999-Dec-31 23:59:60, and the Unix clock counter would read
$N+1$. If the machine's OS did not take leap seconds into account, calling
\texttt{gmtime()} would give 2000-Jan-01 00:00:00.  If the machine's OS took
leap seconds into account, the right date would be produced.  If the machine's
OS did \em{not} take leap seconds into account but ran NTP, the clock
counter would \em{not} increment from $N$ to $N+1$: the leap seconds are
added by holding back the counter.  Hence, calling \texttt{gmtime()} would
produce 1999-Dec-31 23:59:59.  If the machine's OS took leap seconds into
account, and also ran NTP, the counter would increment to $N+1$, and
\texttt{gmtime()} would give 1999-Dec-31 23:59:60.

So, we see that there are two behaviors for the Unix clock counter:
constant increase (if the OS supports leap seconds), or non-constant
increment (if the OS does not support leap seconds and NTP is being used).
In the first case, the difference between GPS and Unix time would be a
constant integer.  In the second, it would change as leap seconds are
added to UTC.  We shall assume that the second case holds, \textit{i.e.} no
OS support for leap seconds (or if there is support, it is not activated)
and NTP is being used.

We care about converting between GPS and UTC because we would like
human-readable time stamps.  Rarely, if ever, would we care about
the actual current time reported by the computer.

GPS to Unix and vice versa can only work for dates > 1980-Jan-06 00:00:00 UTC.
At that time, TAI was ahead of GPS by 19 seconds, and GPS was neither ahead
nor behind UTC.


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


