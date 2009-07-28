/*
*  Copyright (C) 2007 David Chin, Jolien Creighton
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

/* <lalVerbatim file="DateStringCV">
Author: David Chin <dwchin@umich.edu> +1-734-709-9119
$Id$
</lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{DateString.c}}
\label{ss:DateString.c}

Returns a formatted string for the date and time in ISO 8601 format,
given the date in an \texttt{LALDate} structure.

\subsection*{Prototypes}
\vspace{0.1in}
\input{DateStringCP}
\idx{LALDateString()}

\subsubsection*{Description}

Returns a formatted string for the date and time in ISO 8601 format,
given the date in an \texttt{LALDate} structure.  A date and time in
ISO 8601 format looks like 2001-03-04 17:03:52 for March 4, 2001,
5:03:52 pm.  The advantage of this format, besides avoiding Y2K
issues, is that a numerical-order sort of dates will result in a
chronologically ordered sort.  This routine is a replacement for
\texttt{strftime(3c)}.

\subsubsection*{Algorithms}

Trivial.

\subsubsection*{Uses}

Suppose we would like to form a timestamp string for the current time.
The following program would accomplish this:

\begin{verbatim}
#include <lal/Date.h>
int main(int argc, char *argv[])
{
    static LALStatus status;
    LIGOTimeGPS  gpstime;
    LALDate      laldate;
    LALUnixDate  utimestruct;
    CHARVector  *utc = NULL;
    time_t       ticks;

    INITSTATUS (status, "printone", TESTUTOGPSC);

    LALCHARCreateVector(status, &utc, (UINT4)64);

    ticks = time(NULL);
    gmtime_r(&ticks, &(laldate->unixDate));
    DateString(status, utc, &laldate);

    printf("%s\n", utc->data);

    CHARDestroyVector(status, &utc);

    RETURN (status);
}
\end{verbatim}

\subsubsection*{Notes}

See the official ISO document at
\url{http://www.iso.ch/markete/8601.pdf}, and
an overview of the format at
\url{http://www.cl.cam.ac.uk/~mgk25/iso-time.html}.

</lalLaTeX> */

/*#include <time.h>*/	/* for strftime() */
#include <lal/LALRCSID.h>
#include <lal/Date.h>


NRCSID (DATESTRINGC, "$Id$");


/* <lalVerbatim file="DateStringCP"> */
void
LALDateString (LALStatus     *status,
               CHARVector    *timestamp,
               const LALDate *date)
{ /* </lalVerbatim> */
  CHAR tmpmon[3];
  CHAR tmpmday[3];
  CHAR tmphour[3];
  CHAR tmpmin[3];
  CHAR tmpsec[3];
  CHAR tmpwday[4];

  INITSTATUS (status, "LALDateString", DATESTRINGC);

  /*
   * Check pointer to input variable
   */
  ASSERT (date != (LALDate *)NULL, status,
          DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  /*
   * Check that timestamp buffer is large enough
   */
  ASSERT (timestamp->length >= 32, status,
          DATEH_EBUFFTOOSMALL, DATEH_MSGEBUFFTOOSMALL);

  /*
   * Use strftime (3) to form ISO8601-format time stamp, plus day name
   */
  /*
    strftime(timestamp->data, timestamp->length,
           "%Y-%m-%d %H:%M:%S UTC %a", &(date->unixDate));
  */

  /*
   * AVOID STRFTIME() since it causes seg faults on Solaris.  We give
   * the locale-specific day name.
   */
  if (date->unixDate.tm_mon >= LALMONTH_OCT) {
    sprintf(tmpmon, "%2d", date->unixDate.tm_mon + 1);
  } else {
    sprintf(tmpmon, "0%1d", date->unixDate.tm_mon + 1);
  }

  if (date->unixDate.tm_mday >= 10) {
    sprintf(tmpmday, "%2d", date->unixDate.tm_mday);
  } else {
    sprintf(tmpmday, "0%1d", date->unixDate.tm_mday);
  }

  if (date->unixDate.tm_hour >= 10) {
    sprintf(tmphour, "%2d", date->unixDate.tm_hour);
  } else {
    sprintf(tmphour, "0%1d", date->unixDate.tm_hour);
  }

  if (date->unixDate.tm_min >= 10) {
    sprintf(tmpmin, "%2d", date->unixDate.tm_min);
  } else {
    sprintf(tmpmin, "0%1d", date->unixDate.tm_min);
  }

  if (date->unixDate.tm_sec >= 10) {
    sprintf(tmpsec, "%2d", date->unixDate.tm_sec);
  } else {
    sprintf(tmpsec, "0%1d", date->unixDate.tm_sec);
  }

  switch(date->unixDate.tm_wday) {
  case 0:
    sprintf(tmpwday, "Sun");
    break;

  case 1:
    sprintf(tmpwday, "Mon");
    break;

  case 2:
    sprintf(tmpwday, "Tue");
    break;

  case 3:
    sprintf(tmpwday, "Wed");
    break;

  case 4:
    sprintf(tmpwday, "Thu");
    break;

  case 5:
    sprintf(tmpwday, "Fri");
    break;

  case 6:
    sprintf(tmpwday, "Sat");
    break;

  default:
    sprintf(tmpwday, "%c", '\0');
    break;
  }

  sprintf(timestamp->data, "%d-%s-%s %s:%s:%s UTC %s",
          date->unixDate.tm_year + 1900,
          tmpmon, tmpmday, tmphour, tmpmin, tmpsec, tmpwday);

  RETURN (status);
} /* END LALDateString() */


