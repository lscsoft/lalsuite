/* <lalVerbatim file="UtimeCV">
Author: David Chin <dwchin@umich.edu> +1-734-730-1274
$Id$
</lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{Utime.c}}
\label{ss:Utime.c}

Converts UTC time in \texttt{LIGOTimeUnix} format to \texttt{LALDate}
format, correcting for leap seconds.

\subsection*{Prototypes}
\vspace{0.1in}
\input{UtimeCP}
\index{\texttt{LALUtime()}}

\subsubsection*{Description}

This routine converts (UTC) time in a \texttt{LIGOTimeUnix} structure
measured in the Unix epoch to the same time in a \texttt{LALDate}
structure, corrected for leap seconds.

\subsubsection*{Algorithms}

\subsubsection*{Uses}
Here is a simple example (adapted from GRASP) to form a timestamp string:

\begin{verbatim}
void printone(Status *status, const LIGOTimeUnix *time1)
{
    LALDate      laldate;
    LIGOTimeUnix tmp;
    LALUnixDate  utimestruct;
    CHARVector *utc = NULL;

    INITSTATUS (status, "printone", TESTUTOGPSC);

    LALCHARCreateVector(status, &utc, (UINT4)64);

    Utime(status, &laldate, time1);
    LALDateString(status, utc, &laldate, LALLEAPSEC_LOOSE);

    printf("%s\n", utc->data);

    CHARDestroyVector(status, &utc);
    
    RETURN (status);
}
\end{verbatim}

\subsubsection*{Notes}
This routine must be updated for every year beyond 1999 that has leap
seconds.  Leaps seconds for any particular year are set by
observation, so there is no formula for predicting them. See
\url{http://tycho.usno.navy.mil/leapsec.html} for more information on
leap seconds.  A table of current leap seconds may be obtained from
\url{ftp://maia.usno.navy.mil/ser7/tai-utc.dat}. The following is from
the USNO:

\begin{quote}
 Coordinated Universal Time (UTC) is defined by the CCIR Recommendation
 460-4 (1986). It differs from TAI by the total number of leap seconds,
 so that UT1-UTC stays smaller than 0.9s in absolute value.  The decision
 to introduce a leap second in UTC is the responsibility of the
 International Earth Rotation Service (IERS). According to the CCIR
 Recommendation, first preference is given to the opportunities at the
 end of December and June, and second preference to those at the end of
 March and September. Since the system was introduced in 1972, only dates
 in June and December have been used.  TAI is expressed in terms of UTC
 by the relation TAI = UTC + dAT, where dAT is the total algebraic sum of
 leap seconds.
 
 The first leap second was introduced on June 30, 1972. Information on
 the most recent leap second can be found here. The historical list of
 leap seconds lists all of them.
 
 \begin{verbatim}
  The Global Positioning System (GPS) epoch is January 6, 1980 and is
  synchronized to UTC. GPS is NOT adjusted for leap seconds. 
  
  As of 1 January 1999,
         TAI is ahead of UTC   by 32 seconds.
         TAI is ahead of GPS   by 19 seconds.
         GPS is ahead of UTC   by 13 seconds.
 \end{verbatim}
 
 Until 1960, Universal Time (UT) was taken as the independent variable of
 astronomical ephemerides.  UT was then replaced by Ephemeris Time (ET),
 based on the motion of the sun.  However, ET did not include
 relativistic effects, such as corrections for the gravitational
 potential and velocity, as required by advances in the accuracy of time
 comparisons.  Thus ET was superseded in 1981 by Terrestrial Dynamical
 Time (TDT) and Barycentric Dynamical Time (TDB), which distinguish
 coordinate systems with origins at the center of the Earth and the
 center of the solar system, respectively, and are consistent with the
 general theory of relativity.  In the language of general relativity,
 TDT is a proper time while TDB is a coordinate time.  In 1991, TDT was
 renamed simply Terrestrial Time (TT) and two additional relativistic
 time scales, Geocentric Coordinate Time (TCG) and Barycentric Coordinate
 Time (TCB) were adopted.  Definitions of these time scales are given in
 Systems of Time.
  
\end{quote}


</lalLaTeX> */


#include <lal/LALRCSID.h>
#include <lal/Date.h>
#include "date_value.h"

NRCSID (UTIMEC, "$Id$");

/*
 * Unix epoch times when leap seconds are in effect.
 * Data from ftp://maia.usno.navy.mil/ser7/tai-utc.dat
 * This static constant must be updated whenever new leap seconds are
 * introduced. 
 */
static const time_t leaps[]={
  ((2440587-2440587)*SECS_PER_DAY),
  ((2440973-2440587)*SECS_PER_DAY),
  ((2441317-2440587)*SECS_PER_DAY),
  ((2441499-2440587)*SECS_PER_DAY),
  ((2441683-2440587)*SECS_PER_DAY),
  ((2442048-2440587)*SECS_PER_DAY),
  ((2442413-2440587)*SECS_PER_DAY),
  ((2442778-2440587)*SECS_PER_DAY),
  ((2443144-2440587)*SECS_PER_DAY),
  ((2443509-2440587)*SECS_PER_DAY),
  ((2443874-2440587)*SECS_PER_DAY),
  ((2444239-2440587)*SECS_PER_DAY),
  ((2444786-2440587)*SECS_PER_DAY),
  ((2445151-2440587)*SECS_PER_DAY),
  ((2445516-2440587)*SECS_PER_DAY),
  ((2446247-2440587)*SECS_PER_DAY),
  ((2447161-2440587)*SECS_PER_DAY),
  ((2447892-2440587)*SECS_PER_DAY),
  ((2448257-2440587)*SECS_PER_DAY),
  ((2448804-2440587)*SECS_PER_DAY),
  ((2449169-2440587)*SECS_PER_DAY),
  ((2449534-2440587)*SECS_PER_DAY),
  ((2450083-2440587)*SECS_PER_DAY),
  ((2450630-2440587)*SECS_PER_DAY),
  ((2451179-2440587)*SECS_PER_DAY)
};

/* Here's the list from ftp://maia.usno.navy.mil/ser7/tai-utc.dat
 1961 JAN  1 =JD 2437300.5  TAI-UTC=   1.4228180 S + (MJD - 37300.) X 0.001296 S
 1961 AUG  1 =JD 2437512.5  TAI-UTC=   1.3728180 S + (MJD - 37300.) X 0.001296 S
 1962 JAN  1 =JD 2437665.5  TAI-UTC=   1.8458580 S + (MJD - 37665.) X 0.0011232S
 1963 NOV  1 =JD 2438334.5  TAI-UTC=   1.9458580 S + (MJD - 37665.) X 0.0011232S
 1964 JAN  1 =JD 2438395.5  TAI-UTC=   3.2401300 S + (MJD - 38761.) X 0.001296 S
 1964 APR  1 =JD 2438486.5  TAI-UTC=   3.3401300 S + (MJD - 38761.) X 0.001296 S
 1964 SEP  1 =JD 2438639.5  TAI-UTC=   3.4401300 S + (MJD - 38761.) X 0.001296 S
 1965 JAN  1 =JD 2438761.5  TAI-UTC=   3.5401300 S + (MJD - 38761.) X 0.001296 S
 1965 MAR  1 =JD 2438820.5  TAI-UTC=   3.6401300 S + (MJD - 38761.) X 0.001296 S
 1965 JUL  1 =JD 2438942.5  TAI-UTC=   3.7401300 S + (MJD - 38761.) X 0.001296 S
 1965 SEP  1 =JD 2439004.5  TAI-UTC=   3.8401300 S + (MJD - 38761.) X 0.001296 S
 1966 JAN  1 =JD 2439126.5  TAI-UTC=   4.3131700 S + (MJD - 39126.) X 0.002592 S
 1968 FEB  1 =JD 2439887.5  TAI-UTC=   4.2131700 S + (MJD - 39126.) X 0.002592 S
 1972 JAN  1 =JD 2441317.5  TAI-UTC=  10.0       S + (MJD - 41317.) X 0.0      S
 1972 JUL  1 =JD 2441499.5  TAI-UTC=  11.0       S + (MJD - 41317.) X 0.0      S
 1973 JAN  1 =JD 2441683.5  TAI-UTC=  12.0       S + (MJD - 41317.) X 0.0      S
 1974 JAN  1 =JD 2442048.5  TAI-UTC=  13.0       S + (MJD - 41317.) X 0.0      S
 1975 JAN  1 =JD 2442413.5  TAI-UTC=  14.0       S + (MJD - 41317.) X 0.0      S
 1976 JAN  1 =JD 2442778.5  TAI-UTC=  15.0       S + (MJD - 41317.) X 0.0      S
 1977 JAN  1 =JD 2443144.5  TAI-UTC=  16.0       S + (MJD - 41317.) X 0.0      S
 1978 JAN  1 =JD 2443509.5  TAI-UTC=  17.0       S + (MJD - 41317.) X 0.0      S
 1979 JAN  1 =JD 2443874.5  TAI-UTC=  18.0       S + (MJD - 41317.) X 0.0      S
 1980 JAN  1 =JD 2444239.5  TAI-UTC=  19.0       S + (MJD - 41317.) X 0.0      S
 1981 JUL  1 =JD 2444786.5  TAI-UTC=  20.0       S + (MJD - 41317.) X 0.0      S
 1982 JUL  1 =JD 2445151.5  TAI-UTC=  21.0       S + (MJD - 41317.) X 0.0      S
 1983 JUL  1 =JD 2445516.5  TAI-UTC=  22.0       S + (MJD - 41317.) X 0.0      S
 1985 JUL  1 =JD 2446247.5  TAI-UTC=  23.0       S + (MJD - 41317.) X 0.0      S
 1988 JAN  1 =JD 2447161.5  TAI-UTC=  24.0       S + (MJD - 41317.) X 0.0      S
 1990 JAN  1 =JD 2447892.5  TAI-UTC=  25.0       S + (MJD - 41317.) X 0.0      S
 1991 JAN  1 =JD 2448257.5  TAI-UTC=  26.0       S + (MJD - 41317.) X 0.0      S
 1992 JUL  1 =JD 2448804.5  TAI-UTC=  27.0       S + (MJD - 41317.) X 0.0      S
 1993 JUL  1 =JD 2449169.5  TAI-UTC=  28.0       S + (MJD - 41317.) X 0.0      S
 1994 JUL  1 =JD 2449534.5  TAI-UTC=  29.0       S + (MJD - 41317.) X 0.0      S
 1996 JAN  1 =JD 2450083.5  TAI-UTC=  30.0       S + (MJD - 41317.) X 0.0      S
 1997 JUL  1 =JD 2450630.5  TAI-UTC=  31.0       S + (MJD - 41317.) X 0.0      S
 1999 JAN  1 =JD 2451179.5  TAI-UTC=  32.0       S + (MJD - 41317.) X 0.0      S

We interpolate a couple of dates:
 1970 JAN  1 =JD 2440587.5  TAI-UTC=   8.0000820 S
 1970 JAN 22 =JD 2440973.5  TAI-UTC=   9.0005940 S
*/
 
 

/*
 * Return UTC time in LALDate structure, taking leap seconds into account.
 * Replaces the C library's gmtime(). Adapted from the GRASP library by
 * Bruce Allen, et al.
 */
/* <lalVerbatim file="UtimeCP">  */
void
LALUtime ( LALStatus                *status,
           LALDate                  *utc,
           const LIGOTimeUnix       *unixtime)
     /*            const LALLeapSecAccuracy *accuracy) */
{ /* </lalVerbatim> */
  /* latest time for which this routine will work: 2002-Mar-31 23:59:00 */
  /* 24 leap seconds because of the two interpolated ones: 1970-Jan-1 and
   *  1970-Jan-22 */
  const INT4   maxtested = (24*365 + 8*366 + 2*31 + 28)*SECS_PER_DAY - 60 + 24;
  /* number of times leap seconds occur */
  const INT4   numleaps = sizeof(leaps)/sizeof(time_t);
  time_t       tmptime;
  LALUnixDate *tmputc;
  INT4         i;

  LALLeapSecAccuracy  acc = LALLEAPSEC_LOOSE;
  LALLeapSecAccuracy *accuracy = &acc;

    
  INITSTATUS (status, "LALUtime", UTIMEC);

  /*
   * Check pointer to input variable
   */
  ASSERT (unixtime != (LIGOTimeUnix *)NULL, status,
          DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  /*
   * Check for valid input
   */
  /*
  tmptime = unixtime->unixSeconds;
  ASSERT (tmptime >= 0 && tmptime <= maxtested, status,
          DATEH_ERANGE, DATEH_MSGERANGE);
  */
  if (lalDebugLevel > 0)
    {
      fprintf(stderr, "LALUtime: info: maxtested = %ld\n", maxtested);
      tmptime = maxtested - 24;  /* kill leap secs */
      fprintf(stderr, "                %s\n", asctime(gmtime(&tmptime)));
    }
    
  /*
   * Check pointer to output variable.
   */
  ASSERT (utc != (LALDate *)NULL, status,
          DATEH_ENULLOUTPUT, DATEH_MSGENULLOUTPUT);


  /*
   * Check for broken gmtime(). If not broken, use gmtime().
   */

  /* a leap second was added at the end of 1998-Dec-31;
   * this is about 29 years from Unix reference time, plus 24 leap seconds;
   * We check by calling gmtime() on 1998-Dec-31 23:59:59, and then adding
   * one second and calling gmtime() again.  If gmtime takes care of leap
   * seconds, the second call will result in a date that is still in
   * 1998. (This code was written on 2001-Aug-8, and the latest leap second
   * was the one at the end of Dec 1998. */

  /* 1998-Dec-31 23:59:59 */
  tmptime = (22*365 + 7*366 + 7*31 + 4*30 + 28)* SECS_PER_DAY - 1 + 24;
  tmputc  = gmtime(&tmptime);

  if (lalDebugLevel > 0)
    {
      fprintf(stderr, "LALUtime: info: tmputc = %s\n", asctime(tmputc));
    }

  /* system gmtime does take leap seconds into account */
  if (tmputc->tm_year == (time_t)98)
    {
      if (lalDebugLevel > 0)
        {
          fprintf(stderr, "LALUtime: info:  gmtime() takes leap seconds ");
          fprintf(stderr, "into account\n");
        }
      
      /* check that date requested is not later than 2002-Mar-31 23:59:59,
       * which is when the next possible leap second will be. IERS has
       * announced that there will be NO leap second at the end of 2001
       * or any time before */

      /* NOTE: this will break if system gmtime() has taken leap seconds
       * into account in the past (i.e. before the test date */

      /*
       * if date is later
       *    check accuracy param
       *        if anal accuracy
       *            die
       *        else
       *            print warning message
       *
       * // date is not later
       * compute date struct
       */

      if (unixtime->unixSeconds > maxtested)
        {
          /* check accuracy param */
          if (*accuracy == LALLEAPSEC_STRICT)  /* strict accuracy */
            {
              if (lalDebugLevel > 0)
                {
                  fprintf(stderr,
                          "LALUtime: error: may be missing leap seconds\n");
                }

              /*
               * DIE
               */
            }
          else if (*accuracy == LALLEAPSEC_LOOSE) /* loose accuracy */
            {
              if (lalDebugLevel > 0)
                {
                  fprintf(stderr,
                          "LALUtime: warning: maybe missing leap seconds\n");
                }
            }
        }

      /* compute date struct */
      tmptime = unixtime->unixSeconds;
      tmputc = gmtime(&tmptime);
      utc->unixDate.tm_sec   = tmputc->tm_sec;
      utc->unixDate.tm_min   = tmputc->tm_min;
      utc->unixDate.tm_hour  = tmputc->tm_hour;
      utc->unixDate.tm_mday  = tmputc->tm_mday;
      utc->unixDate.tm_mon   = tmputc->tm_mon;
      utc->unixDate.tm_year  = tmputc->tm_year;
      utc->unixDate.tm_wday  = tmputc->tm_wday;
      utc->unixDate.tm_yday  = tmputc->tm_yday;
      utc->unixDate.tm_isdst = 0;    /* always ignore tm_isdst field */
    }
  else /* system gmtime() does NOT take leap secs into account */
    {
      if (lalDebugLevel > 0)
        {
          fprintf(stderr, "LALUtime: info: gmtime() does not figure in ");
          fprintf(stderr, "leap seconds\n");
        }
      
      /* fix up leap seconds */
      tmptime = unixtime->unixSeconds;
      i       = 0;
      while (i < numleaps && leaps[i] + i - 1 < tmptime)
        ++i;

      if (tmptime == (leaps[i] + i - 1))
        {
          tmptime -= i;
          tmputc   = gmtime(&tmptime);
          utc->unixDate.tm_sec   = 60;
          utc->unixDate.tm_min   = tmputc->tm_min;
          utc->unixDate.tm_hour  = tmputc->tm_hour;
          utc->unixDate.tm_mday  = tmputc->tm_mday;
          utc->unixDate.tm_mon   = tmputc->tm_mon;
          utc->unixDate.tm_year  = tmputc->tm_year;
          utc->unixDate.tm_wday  = tmputc->tm_wday;
          utc->unixDate.tm_yday  = tmputc->tm_yday;
          utc->unixDate.tm_isdst = 0;
        }
      else
        {
          tmptime -= (i - 1);
          tmputc   = gmtime(&tmptime);
          utc->unixDate.tm_sec   = tmputc->tm_sec;
          utc->unixDate.tm_min   = tmputc->tm_min;
          utc->unixDate.tm_hour  = tmputc->tm_hour;
          utc->unixDate.tm_mday  = tmputc->tm_mday;
          utc->unixDate.tm_mon   = tmputc->tm_mon;
          utc->unixDate.tm_year  = tmputc->tm_year;
          utc->unixDate.tm_wday  = tmputc->tm_wday;
          utc->unixDate.tm_yday  = tmputc->tm_yday;
          utc->unixDate.tm_isdst = 0;
        }
    }
      
  /* set residual nanoseconds */
  utc->residualNanoSeconds = unixtime->unixNanoSeconds;

  RETURN (status);
} /* END LALUtime() */

