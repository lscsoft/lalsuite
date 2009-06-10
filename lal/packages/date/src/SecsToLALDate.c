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

/* <lalVerbatim file="SecsToLALDateCV">
Author: David Chin <dwchin@umich.edu> +1-734-709-9119
$Id$
</lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{SecsToLALDate.c}}
\label{ss:SecsToLALDate.c}

Converts time in seconds to time in an \texttt{LALDate} structure.

\subsection*{Prototypes}
\vspace{0.1in}

\input{SecsToLALDateCP}
\idx{LALSecsToLALDate()}

\subsubsection*{Description}

This routine converts a time of day in decimal seconds since 0h (midnight)
to an LALDate structure.  Of course, the date information is not present.


\subsubsection*{Algorithms}

\subsubsection*{Uses}

A simple example:

\begin{verbatim}
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/Date.h>

INT4 debuglevel = 2;

NRCSID (TESTLMSTC, "Id");

int
main(int argc, char *argv[])
{
    LALDate date;
    LALDate mstdate;
    REAL8   gmstsecs;
    CHAR    timestamp[64], tmpstr[64];
    Status  status = {0};


    INITSTATUS(&status, "TestLMST", TESTLMSTC);

    printf("TEST of GMST1 routine\n");
    printf("=====================\n");

    date.unixDate.tm_sec  = 0;
    date.unixDate.tm_min  = 0;
    date.unixDate.tm_hour = 0;
    date.unixDate.tm_mday = 16;
    date.unixDate.tm_mon  = LALMONTH_NOV;
    date.unixDate.tm_year = 94;

    GMST1(&status, &gmstsecs, &date, MST_SEC);
    SecsToLALDate(&status, &mstdate, gmstsecs);
    strftime(timestamp, 64, "%Hh %Mm %S", &(mstdate.unixDate));
    sprintf(tmpstr, "%fs", mstdate.residualNanoSeconds * 1.e-9);
    strcat(timestamp, tmpstr+1);

    printf("gmst = %s\n", timestamp);
    printf("    expect: 03h 39m 20.6222s \n");

    return(0);
}

\end{verbatim}

\subsubsection*{Notes}

</lalLaTeX> */


/*-----------------------------------------------------------------------
 *
 * SYNOPSIS
 *
 * LALSecsToLALDate(): Converts time in seconds to time in an LALDate
 *                  structure.
 *
 * DESCRIPTION
 *
 * LALSecsToLALDate():
 *       Inputs:  REAL8  seconds  -- time in seconds since 0h (midnight)
 *
 *       Outputs: LALDate *date   -- time in LALDate structure.  Of course,
 *                                   only the hour, min, and sec fields are
 *                                   set since there is no information on
 *                                   the date, timezone, etc.
 *
 * DIAGNOSTICS
 * (Abnormal termination conditions, error and warning codes summarized
 * here. More complete descriptions are found in documentation.)
 *
 * CALLS
 *
 * NOTES
 *
 *----------------------------------------------------------------------- */

#include <lal/LALRCSID.h>

NRCSID (SECSTOLALDATEC, "$Id$");

#include <math.h>
#include <lal/Date.h>
#include "date_value.h"

/*
 * Convert time in seconds to LALDate struct
 */

/* <lalVerbatim file="SecsToLALDateCP"> */
void
LALSecsToLALDate(LALStatus  *status,
                 LALDate    *date,    /* output - date */
                 REAL8       seconds) /* input - time in seconds since 0h */
{ /* </lalVerbatim> */
    REAL8 hr, min, sec;
    REAL8 tmpsec;
    REAL8 dum;

    INITSTATUS (status, "SECSTOLALDATE", SECSTOLALDATEC);

    /*
     * Check pointer to output variable
     */
    ASSERT (date != (LALDate *)NULL, status,
            DATEH_ENULLOUTPUT, DATEH_MSGENULLOUTPUT);

    /*
     * Make sure time is positive, i.e. seconds == #secs since 0h
     */
    tmpsec = fmod(seconds, (REAL8)SECS_PER_DAY);
    while (tmpsec < 0)
        tmpsec += (REAL8)SECS_PER_DAY;

    hr  = tmpsec / (REAL8)SECS_PER_HOUR;
    min = fmod(hr * (REAL8)MINS_PER_HOUR, (REAL8)MINS_PER_HOUR);
    sec = fmod(min * (REAL8)SECS_PER_MIN, (REAL8)SECS_PER_MIN);

    /*
     * Set the non-required fields to zero
     */
    date->unixDate.tm_mday = 0;
    date->unixDate.tm_mon  = LALMONTH_JAN;
    date->unixDate.tm_year = 0;
    date->unixDate.tm_wday = 0;
    date->unixDate.tm_yday = 0;
    date->unixDate.tm_isdst = 0;


    /*
     * Set the LALDate structure. Only the time fields matter, of course.
     * Fractional seconds become "residual".
     */
    date->unixDate.tm_hour = (INT4)hr;
    date->unixDate.tm_min  = (INT4)min;
    date->unixDate.tm_sec  = (INT4)sec;
    date->residualNanoSeconds = modf(sec, &dum) * 1.e+09;


    RETURN (status);
} /* END LALSecsToLALDate() */
