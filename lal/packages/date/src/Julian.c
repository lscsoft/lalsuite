/* <lalVerbatim file="JulianCV">
Author: David Chin <dwchin@umich.edu> +1-734-730-1274
$Id$
</lalVerbatim> */

/* <lalLaTeX>
\subsection{Module \texttt{Julian.c}}
\label{ss:Julian.c}

Converts between Gregorian date and Julian Days/Dates.


\subsection*{Prototypes}
\vspace{0.1in}
\input{JulianCP}
\index{\texttt{LALJulianDay()}}
\index{\texttt{LALJulianDate()}}
\index{\texttt{LALModJulianDate()}}


\subsubsection*{Description}

These routines compute Julian Day, Julian Date, and Modified Julian Date
for a given Gregorian date and time UTC.  Julian Day and Modified Julian
Day are integer number of days; Julian Date and Modified Julian Date are
decimal number of days.

\subsubsection*{Algorithms}

See~\cite{esaa:1992} and~\cite{green:1985} for details. First, some
definitions:  

\begin{tabular}{rcl}
  Mean Julian Year & = & 365.25 days \\
  Julian Epoch     & = & 1 Jan 4713BCE, 12:00 GMT (4713 BC Jan 01d.5 GMT) \\
  Fundamental Epoch J2000.0 & = & 2001-01-01.5 TDB
\end{tabular}

\emph{Julian Date} is the amount of time elapsed since the Julian
Epoch, measured in days and fractions of a day.  There are a couple of
complications arising from the length of a year: the \emph{Tropical
Year} is 365.2422 days. First, the Gregorian correction where 10 days
(1582-10-05 through 1582-10-14) were eliminated.  Second, leap years:
years ending with two zeroes (\textit{e.g.} 1700, 1800) are leap only
if divisible by 400; so, 400 civil years contain $(400 \times 365.25)
- 3 = 146097$ days.  So, the Julian Date of J2000.0 is JD 2451545.0,
and thus the Julian Epoch = J2000.0 + (JD - 2451545)/365.25,
\textit{i.e.} number of years elapsed since J2000.0.

One algorithm for computing the Julian Day is from~\cite{vfp:1979}
based on a formula in~\cite{esaa:1992} where the algorithm is due
to~\cite{fvf:1968} and ``compactified'' by P.~M. Muller and
R.~N. Wimberly.  The formula is
%
\begin{displaymath}
  jd = 367 \times y - 7 \times (y + (m + 9)/12)/4 - 
  3 \times ((y + (m - 9)/7)/100 + 1)/4
        + 275 \times m/9 + d + 1721029
\end{displaymath}
%
where $jd$ is the Julian day number, $y$ is the year, $m$ is the month
(1-12), and $d$ is the day (1-31).  This formula is valid only for 
$\mathrm{JD} \ge 0$, \textit{i.e.} after -4713 Nov 23 = 4712 BCE Nov 23.

A shorter formula from the same reference, but which only works for
dates since 1900-March is:
%
\begin{displaymath}
  jd = 367 \times y - 7 \times (y + (m + 9)/12)/4 + 275 \times m/9 + 
  d + 1721014
\end{displaymath}
%
We will use this shorter formula since there is unlikely to be any
analyzable data from before 1900-Mar.

\subsubsection*{Uses}

Suppose we would like to get the Julian Date for
today.  The following program would accomplish this:

\begin{verbatim}
#include <lal/LALStdlib.h>
#include <lal/Date.h>

INT4 debuglevel = 2;

NRCSID (TESTJULIANDAYC, "Id");

int
main(int argc, char *argv[])
{
    time_t        now;
    Status        status = {0};
    LALDate       date;
    REAL8         jDate;

    INITSTATUS (&status, "TestJulianDay", TESTJULIANDAYC);

    now = time(NULL);
    gmtime_r(&now, &(date.unixDate));

    date.unixDate.tm_sec  = ltime->tm_sec;
    date.unixDate.tm_min  = ltime->tm_min;
    date.unixDate.tm_hour = ltime->tm_hour;
    date.unixDate.tm_mday = ltime->tm_mday;
    date.unixDate.tm_mon  = ltime->tm_mon;
    date.unixDate.tm_year = ltime->tm_year;
    date.unixDate.tm_wday = ltime->tm_wday;
    date.unixDate.tm_yday = ltime->tm_yday;
    date.unixDate.tm_isdst = 0; 

    LALJulianDate(&status, &jDate, &date);
    printf("\tJulian Date                = %10.1f\n", jDate);

    return 0;
}


\end{verbatim}

\subsubsection*{Notes}


</lalLaTeX> */



#include <lal/LALRCSID.h>
#include <lal/Date.h>
#include "date_value.h"

NRCSID (JULIANC, "$Id$");

/*
 * Compute Julian Day for given Gregorian date
 */
/* <lalVerbatim file="JulianCP"> */
void
LALJulianDay (LALStatus     *status,
              INT4          *jDay,
              const LALDate *date)
{ /* </lalVerbatim> */
    INT4 y, m, d;
    
    INITSTATUS (status, "LALJulianDay", JULIANC);

    /*
     * Check pointer to input variable
     */
    ASSERT (date != (LALDate *)NULL, status,
            DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

    /*
     * Check pointer to output variable:
     */
    ASSERT (jDay != (INT4 *)NULL, status,
            DATEH_ENULLOUTPUT, DATEH_MSGENULLOUTPUT);

    /*
     * Pull out Year, Month, Day, and convert from the
     * struct tm value ranges
     */
    y = (INT4)((date->unixDate).tm_year) + 1900;
    m = (INT4)((date->unixDate).tm_mon) + 1;

    /* Julian Day begins at noon, so fix day of month if necessary */
    /* if ((date->unixDate).tm_hour < 12)
        d = (date->unixDate).tm_mday - 1;
        else */

    d = (date->unixDate).tm_mday;



    /*
     * Check for legal input: Input date must be after 1900-Mar
     * Recall: !(A && B) == (!A || !B)
     */
    if (y < 1900 || (y == 1900 && m < 3))
    {
        ABORT (status, DATEH_EDATETOOEARLY,
               DATEH_MSGEDATETOOEARLY);
    }

    /*
     * Compute Julian Day
     */
    *jDay = (INT4)(367 * y
                   - 7 * (y + (m + 9)/12) / 4
                   + 275 * m/9
                   + d + 1721014);

    RETURN (status);
} /* END LALJulianDay() */



/*
 * Compute Julian Date for given Gregorian date and UTC time
 */
/* <lalVerbatim file="JulianCP"> */
void
LALJulianDate (LALStatus     *status,
               REAL8         *jDateOut,
               const LALDate *date)
{ /* </lalVerbatim> */
    INT4  hr, min, sec;
    REAL8 rns;          /* residual nanoseconds */
    INT4  jday;
    REAL8 jdate;

    INITSTATUS(status, "LALJulianDate", JULIANC);
    ATTATCHSTATUSPTR(status);

    /*
     * Check pointer to input variable
     */
    ASSERT (date != (LALDate *)NULL, status,
            DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

    /*
     * Check pointer to output variable:
     */
    ASSERT (jDateOut != (REAL8 *)NULL, status, DATEH_ENULLOUTPUT,
            DATEH_MSGENULLOUTPUT);

    /*
     * Extract Hour, Minute, Second, and residual nanoseconds
     */
    hr  = (date->unixDate).tm_hour;
    min = (date->unixDate).tm_min;
    sec = (date->unixDate).tm_sec;
    rns = date->residualNanoSeconds * (REAL8)1.e-09;

    /*
     * Get Julian Day number
     */
    TRY( LALJulianDay( status->statusPtr, &jday, date ), status );

    /*
     * Convert to fractions of a day
     */
    jdate  = (REAL8)jday - (REAL8)0.5;
    jdate += (REAL8)hr / (REAL8)24.;
    jdate += (REAL8)min / (REAL8)1440.;
    jdate += ((REAL8)sec + rns) / (REAL8)86400.;

    *jDateOut = jdate;

    DETATCHSTATUSPTR(status);
    RETURN (status);
} /* END LALJulianDate() */



/*
 * Compute Modified Julian Date for given Gregorian date and UTC time
 */
/* <lalVerbatim file="JulianCP"> */
void
LALModJulianDate (LALStatus     *status,
                  REAL8         *modJDate,
                  const LALDate *date)
{ /* </lalVerbatim> */
  REAL8 mjd;

  INITSTATUS(status, "LALModJulianDate", JULIANC);
  ATTATCHSTATUSPTR(status);

  /*
   * Check pointer to input variable
   */
  ASSERT (date != (LALDate *)NULL, status,
          DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  /*
   * Check pointer to output variable:
   */
  ASSERT (modJDate != (REAL8 *)NULL, status,
          DATEH_ENULLOUTPUT, DATEH_MSGENULLOUTPUT);

  /*
   * Get Julian Date, and modify it
   */
  TRY( LALJulianDate( status->statusPtr, &mjd, date ), status );
  mjd -= MJDREF;
    
  *modJDate = mjd;

  DETATCHSTATUSPTR(status);
  RETURN (status);
} /* END LALModJulianDate() */


