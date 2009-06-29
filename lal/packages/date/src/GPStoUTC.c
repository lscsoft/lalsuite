/*
*  Copyright (C) 2007 Bruce Allen, Duncan Brown, David Chin, Jolien Creighton, Patrick Brady, Peter Shawhan, Reinhard Prix
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

/* <lalVerbatim file="GPStoUTCCV">
Author: David Chin <dwchin@umich.edu> +1-734-709-9119
$Id$
</lalVerbatim> */

/* <lalLaTeX>
\subsection{Module \texttt{GPStoUTC.c}}
\label{ss:GPStoUTC.c}

Converts between GPS time (in seconds and nanoseconds) and UTC in a
\texttt{LALDate} structure.

\subsection*{Prototypes}
\vspace{0.1in}
\input{GPStoUTCCP}
\idx{LALGPStoUTC()}
\idx{LALUTCtoGPS()}
\idx{LALLeapSecs()}

\subsubsection*{Description}

\texttt{LALGPStoUTC()} and \texttt{LALUTCtoGPS} convert time in GPS seconds
and nanoseconds (\texttt{LIGOTimeGPS}) and time in UTC (\texttt{LALDate}),
taking into account leap seconds until 2006-Dec-31 23:59:59 UTC. % UPDATEME %

\texttt{LALLeapSecs()} returns the number of leap seconds introduced since
the GPS epoch 1980-Jan-06, abbreviated GPS-UTC.


\subsubsection*{Algorithms}

The conversion from GPS to UTC is copied directly from
GRASP~\cite{grasp:194}.  It does the conversion by counting TAI seconds
starting from the Unix epoch origin, 1970-Jan-01 00:00:00 UTC.  A static
table of leap seconds is compiled in: this \emph{must} be updated whenever
a new leap second is introduced.  The latest leap second included is
2006-Jan-01. % UPDATEME %

The conversion from UTC to GPS is done by counting the amount of elapsed
time since the GPS epoch origin, 1980-Jan-06 00:00:00 UTC.  Again, leap
seconds are accounted for by a static table (different from the one used in
GPS to UTC) which \emph{must} be updated whenever a new leap second is
introduced.  The latest leap second included is 2006-Jan-01. % UPDATEME %

The computation of GPS-UTC is from a static table published by the USNO
at \url{ftp://maia.usno.navy.mil/ser7/tai-utc.dat}.  The latest leap second
included is 2006-Jan-01. % UPDATEME %

\subsubsection*{Uses}

\subsubsection*{Notes}

These routines will not work for times before 1980-01-06 00:00:00 UTC (GPS
0).  The latest leap second that can be accounted for is the one added at
the end of 2005-Dec. % UPDATEME %  These routines have accurate leap second
information until 2006-Dec-31. % UPDATEME %

\textbf{Example:} To convert a GPS time to UTC, and then back to GPS:

\begin{verbatim}

#include <lal/LALStdlib.h>
#include <lal/Date.h>

struct tm *gmtime_r( const time_t *, struct tm * );
char *asctime_r( const struct tm *, char *, int );

int main(int argc, char *argv[])
{
    static LALStatus   status;
    LIGOTimeGPS        gps = {615081613, 123456789};
    LALDate            date;
    LALLeapSecAccuracy accuracy = LALLEAPSEC_STRICT;
    CHARVector        *timestamp = NULL;

    LALCHARCreateVector(&status, &timestamp, (UINT4)128);

    LALGPStoUTC(&status, &date, &gps, &accuracy);

    LALDateString(&status, timestamp, &date);

    printf("GPS (%d, %d) = %s\n", gps.gpsSeconds, gps.gpsNanoSeconds,
           timestamp->data);

    LALUTCtoGPS(&status, &gps, &date, &accuracy);

    printf("%s = GPS (%d, %d)\n", timestamp->data, gps.gpsSeconds,
           gps.gpsNanoSeconds);

    return 0;
}

\end{verbatim}

For an example of how \texttt{LALLeapSecs()} is used, see the test program
\texttt{TestLeapSecs.c} in the \texttt{packages/date/test} directory.


</lalLaTeX> */

#include <lal/LALRCSID.h>

NRCSID (GPSTOUTCC, "$Id$");

#include <lal/LALStdio.h>
#include <lal/Date.h>
#include "date_value.h"
#include <lal/XLALError.h>

#define INFOSTR_LEN 256

struct tm *gmtime_r( const time_t *, struct tm * );
char *asctime_r( const struct tm *, char * );


/* UPDATEME */
/* The international earth rotation service announces leap seconds;
 * their data is posted in bulletin C at
 *    http://www.iers.org/MainDisp.csl?pid=36-9
 * Another useful web site is at
 *    http://maia.usno.navy.mil/
 * which carries these announcements too.  The latest time for which
 * this routine will work: 2006-12-31 23:59:59 UTC
 */

/* GPS for maxtestedGPS computed using lalapps_tconvert (part of ligotools) */
static const INT4 maxtestedGPS = 851644813;

/*
 * Convert GPS seconds to UTC date-time contained in LALDate structure
 */
/* <lalVerbatim file="GPStoUTCCP"> */
void
LALGPStoUTC (LALStatus                *status,
             LALDate                  *p_utcDate,  /* output - date */
             const LIGOTimeGPS        *p_gpsTime,  /* input - GPS seconds */
             const LALLeapSecAccuracy *p_accuracy) /* accuracy of
                                                      leap-second
                                                      accounting:
                                                      LALLEAPSEC_LOOSE, or LALLEAPSEC_STRICT */
{ /* </lalVerbatim> */
  /* UPDATEME -- to update, add an entry at the end listing
   * the Unix epoch time when a leap second is introduced */
  /* this is a table of Unix epoch times when leap
   * seconds were introduced */
  /* What's the funny format? These Unix epoch times are expressed
   * in terms of Julian Day Numbers, i.e. say JD1 = Julian Day number
   * of the day when a leap sec was added, JD0 = Unix epoch origin (i.e. 1970-01-01),
   * then the Unix epoch time is given by
   * (JD1 - JD0) * SECS_PER_DAY
   * FIXME: at some point, we should just dispense with this wacky Julian Day
   * stuff. It's just confusing. */
  /* As of 2005-07-05 08:05 UTC-4, Bulletin C has been released, but the Naval
   * Observatory's tai-utc.dat hasn't been updated. --dwchin */
  static const time_t leaps[] = {
             0,
      33350400,
      63072000,
      78796800,
      94694400,
     126230400,
     157766400,
     189302400,
     220924800,
     252460800,
     283996800,
     315532800,
     362793600,
     394329600,
     425865600,
     489024000,
     567993600,
     631152000,
     662688000,
     709948800,
     741484800,
     773020800,
     820454400,
     867715200,
     915148800,
    1136073600,
  };

  /* number of times leap seconds occur */
  static const INT4   numleaps = sizeof(leaps)/sizeof(time_t);
  time_t       unixTime;
  time_t       tmptime;
  LALUnixDate  tmputc;
  char         tmpstamp[32];
  CHAR         infostr[INFOSTR_LEN];
  INT4         i;

  INITSTATUS (status, "LALGPStoUTC", GPSTOUTCC);

  ASSERT (p_gpsTime != (LIGOTimeGPS *)NULL, status,
          DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  ASSERT (p_gpsTime->gpsSeconds >= 0, status,
          DATEH_ERANGEGPSABS, DATEH_MSGERANGEGPSABS);

  ASSERT (p_accuracy != (LALLeapSecAccuracy *)NULL, status,
          DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  ASSERT ((*p_accuracy == LALLEAPSEC_STRICT ||
           *p_accuracy == LALLEAPSEC_LOOSE), status,
          DATEH_EACCPARAMOUTOFRANGE, DATEH_MSGEACCPARAMOUTOFRANGE);

  ASSERT (p_utcDate != (LALDate *)NULL, status,
          DATEH_ENULLOUTPUT, DATEH_MSGENULLOUTPUT);

  XLALPrintDeprecationWarning("LALGPStoUTC", "XLALGPSToUTC");

  /* we use Unix epoch as our origin */
  unixTime = p_gpsTime->gpsSeconds + UNIXGPS;

  if (lalDebugLevel & LALINFO)
  {
    snprintf(infostr, INFOSTR_LEN, "Max. tested GPS is %d\n", maxtestedGPS);
    LALInfo(status, infostr);
  }

  /* 1998-Dec-31 23:59:60 (if leap seconds are taken into account) */
  tmptime = (22*365 + 7*366)* SECS_PER_DAY + 23;
  gmtime_r(&tmptime, &tmputc);

  if (lalDebugLevel & LALINFO)
  {
    asctime_r(&tmputc, tmpstamp);
    snprintf(infostr, INFOSTR_LEN, "tmputc = %s\n", tmpstamp);
    LALInfo(status, infostr);
  }

  /*
   * if GPS is later than maxtestedGPS
   *    check accuracy param
   *        if anal accuracy
   *            die
   *        else
   *            print warning message
   */
  if (p_gpsTime->gpsSeconds > maxtestedGPS)
    {
      /* check accuracy param */
      if (*p_accuracy == LALLEAPSEC_STRICT)
        {
          ABORT(status, DATEH_ERANGEGPSABS, DATEH_MSGERANGEGPSABS);
        }
      else if (*p_accuracy == LALLEAPSEC_LOOSE)
        {
          LALWarning(status, "may be missing leap seconds");
        }
      else
        {
          LALWarning(status, "may be missing leap seconds");
        }
    }




  /* system gmtime does take leap seconds into account */
  if (tmputc.tm_sec == 60)
    {
      LALInfo(status, "gmtime_r() takes leap seconds into account");

      /* check that date requested is not later than 2002-Mar-31 23:59:59,
       * which is when the next possible leap second will be. IERS has
       * announced that there will be NO leap second at the end of 2001
       * or any time before */

      /* NOTE: this will break if system gmtime() has taken leap seconds
       * into account in the past (i.e. before the test date) */

      /* compute date struct */
      gmtime_r(&unixTime, &tmputc);
      p_utcDate->unixDate.tm_sec   = tmputc.tm_sec;
      p_utcDate->unixDate.tm_min   = tmputc.tm_min;
      p_utcDate->unixDate.tm_hour  = tmputc.tm_hour;
      p_utcDate->unixDate.tm_mday  = tmputc.tm_mday;
      p_utcDate->unixDate.tm_mon   = tmputc.tm_mon;
      p_utcDate->unixDate.tm_year  = tmputc.tm_year;
      p_utcDate->unixDate.tm_wday  = tmputc.tm_wday;
      p_utcDate->unixDate.tm_yday  = tmputc.tm_yday;
      p_utcDate->unixDate.tm_isdst = 0;    /* always ignore tm_isdst field */
    }
  else /* system gmtime() does NOT take leap secs into account */
    {
      LALInfo(status, "gmtime_r() does not figure in leap seconds");

      /* fix up leap seconds */
      i       = 0;
      while (i < numleaps && leaps[i] + i - 1 < unixTime)
        ++i;

      if (lalDebugLevel & LALINFO)
      {
        snprintf(infostr, INFOSTR_LEN, "unixTime = %ld; leaps[%d] = %ld",
            unixTime, i, leaps[i]);
        LALInfo(status, infostr);
      }

      if (unixTime == (leaps[i] + i - 1))
        {
          unixTime -= i;
          gmtime_r(&unixTime, &tmputc);
          p_utcDate->unixDate.tm_sec   = 60;
          p_utcDate->unixDate.tm_min   = tmputc.tm_min;
          p_utcDate->unixDate.tm_hour  = tmputc.tm_hour;
          p_utcDate->unixDate.tm_mday  = tmputc.tm_mday;
          p_utcDate->unixDate.tm_mon   = tmputc.tm_mon;
          p_utcDate->unixDate.tm_year  = tmputc.tm_year;
          p_utcDate->unixDate.tm_wday  = tmputc.tm_wday;
          p_utcDate->unixDate.tm_yday  = tmputc.tm_yday;
          p_utcDate->unixDate.tm_isdst = 0;
        }
      else
        {
          unixTime -= (i - 1);
          gmtime_r(&unixTime, &tmputc);
          p_utcDate->unixDate.tm_sec   = tmputc.tm_sec;
          p_utcDate->unixDate.tm_min   = tmputc.tm_min;
          p_utcDate->unixDate.tm_hour  = tmputc.tm_hour;
          p_utcDate->unixDate.tm_mday  = tmputc.tm_mday;
          p_utcDate->unixDate.tm_mon   = tmputc.tm_mon;
          p_utcDate->unixDate.tm_year  = tmputc.tm_year;
          p_utcDate->unixDate.tm_wday  = tmputc.tm_wday;
          p_utcDate->unixDate.tm_yday  = tmputc.tm_yday;
          p_utcDate->unixDate.tm_isdst = 0;
        }
    }

  /* set residual nanoseconds */
  p_utcDate->residualNanoSeconds = p_gpsTime->gpsNanoSeconds;

  RETURN (status);
} /* END: LALGPStoUTC() */


/*
 * Returns no. of days in year of given date
 */
static int days_in_year(const LALDate *p_utcDate)
{
  int year = p_utcDate->unixDate.tm_year + 1900;

  /* Deal with the years ending with '00: only multiples of 400 are leap */
  if (year % 100  == 0)
    {
      if (year % 400 == 0)
        return 366;
      else
        return 365;
    }

  /* non-'00' years */
  if (year % 4 == 0)
    return 366;

  return 365;
}

/*
 * Returns no. of days in month of given date
 */
static int days_in_month(const LALDate *p_utcDate)
{
  int month = p_utcDate->unixDate.tm_mon;

  switch (month) {
  case LALMONTH_JAN:
  case LALMONTH_MAR:
  case LALMONTH_MAY:
  case LALMONTH_JUL:
  case LALMONTH_AUG:
  case LALMONTH_OCT:
  case LALMONTH_DEC:
    return 31;

  case LALMONTH_APR:
  case LALMONTH_JUN:
  case LALMONTH_SEP:
  case LALMONTH_NOV:
    return 30;

  case 1:
    if (days_in_year(p_utcDate) == 366)
      return 29;
    else
      return 28;
  }

  return -1;
}

/*
 * Struct for leap seconds
 */
typedef struct leap_sec
{
  int    year;       /* year - 1900 */
  int    mon;        /* 0 through 11 */
  INT4   leapsec;
}
leap_sec_t;

/* <lalVerbatim file="GPStoUTCCP"> */
void
LALUTCtoGPS (LALStatus                *status,
             LIGOTimeGPS              *p_gpsTime,  /* output - GPS seconds */
             const LALDate            *p_utcDate,  /* input - date in UTC */
             const LALLeapSecAccuracy *p_accuracy) /* accuracy of
                                                      leap-second
                                                      accounting:
                                                      LALLEAPSEC_LOOSE, or LALLEAPSEC_STRICT */
{ /* </lalVerbatim> */

  /* UPDATEME */
  /*
   * Table of leap seconds
   * Format: Year, Month, and Day of Month in struct tm definition.
   *    (see ctime(3))
   */
  static leap_sec_t leap_sec_data[] =
    {
      {72, LALMONTH_JAN, 1},
      {73, LALMONTH_JAN, 1},
      {74, LALMONTH_JAN, 1},
      {75, LALMONTH_JAN, 1},
      {76, LALMONTH_JAN, 1},
      {77, LALMONTH_JAN, 1},
      {78, LALMONTH_JAN, 1},
      {79, LALMONTH_JAN, 1},
      {80, LALMONTH_JAN, 1},
      {81, LALMONTH_JUL, 1},
      {82, LALMONTH_JUL, 1},
      {83, LALMONTH_JUL, 1},
      {85, LALMONTH_JUL, 1},
      {88, LALMONTH_JAN, 1},
      {90, LALMONTH_JAN, 1},
      {91, LALMONTH_JAN, 1},
      {92, LALMONTH_JUL, 1},
      {93, LALMONTH_JUL, 1},
      {94, LALMONTH_JUL, 1},
      {96, LALMONTH_JAN, 1},
      {97, LALMONTH_JUL, 1},
      {99, LALMONTH_JAN, 1},
      {106, LALMONTH_JAN, 1},
    };


  int ddays = 0;
  int dsecs = 0;
  LALDate tmpdate;
  static LALDate gpsref;
  int i = 0;
  char infostr[256];
  static const int nleaps = sizeof(leap_sec_data)/sizeof(leap_sec_t);

  XLALPrintDeprecationWarning("LALUTCtoGPS", "XLALUTCToGPS");

  /* When GPS began */
  gpsref.unixDate.tm_sec  = 0;
  gpsref.unixDate.tm_min  = 0;
  gpsref.unixDate.tm_hour = 0;
  gpsref.unixDate.tm_mday = 6;
  gpsref.unixDate.tm_mon  = LALMONTH_JAN;
  gpsref.unixDate.tm_year = 80;
  gpsref.unixDate.tm_wday = 0;
  gpsref.unixDate.tm_yday = 0;

  if (lalDebugLevel & LALINFO)
  {
    snprintf(infostr, INFOSTR_LEN, "Date given: %d-%d-%d %d:%d:%d %d\n",
        p_utcDate->unixDate.tm_year+1900, p_utcDate->unixDate.tm_mon+1,
        p_utcDate->unixDate.tm_mday, p_utcDate->unixDate.tm_hour,
        p_utcDate->unixDate.tm_min, p_utcDate->unixDate.tm_sec,
        p_utcDate->residualNanoSeconds);

    LALInfo(status, infostr);
  }

  INITSTATUS(status, "LALUTCtoGPS", GPSTOUTCC);

  ASSERT (p_gpsTime != (LIGOTimeGPS *)NULL, status,
          DATEH_ENULLOUTPUT, DATEH_MSGENULLOUTPUT);

  ASSERT (p_accuracy != (LALLeapSecAccuracy *)NULL, status,
          DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  ASSERT ((*p_accuracy == LALLEAPSEC_STRICT ||
           *p_accuracy == LALLEAPSEC_LOOSE), status,
          DATEH_EACCPARAMOUTOFRANGE, DATEH_MSGEACCPARAMOUTOFRANGE);

  ASSERT (p_utcDate != (LALDate *)NULL, status,
          DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  /* Can't convert dates before 1980-Jan-06 */
  ASSERT (p_utcDate->unixDate.tm_year > 80 ||
          (p_utcDate->unixDate.tm_year == 80 &&
           (p_utcDate->unixDate.tm_mon > LALMONTH_JAN ||
            (p_utcDate->unixDate.tm_mon == LALMONTH_JAN &&
             p_utcDate->unixDate.tm_mday >= 6))), status,
          DATEH_EGPSDATETOOEARLY, DATEH_MSGEGPSDATETOOEARLY);

  /* UPDATEME -- to update, fix the comment and the first if() statement */
  /*
   * Check that time asked for is not after last known leap sec
   * Use by: 2006-Dec-31 23:59:59 UTC
   * Check bulletins such as the following to see if additional ones are needed:
   * http://hpiers.obspm.fr/eoppc/bul/bulc/bulletinc.dat
   * if date is later
   *    check accuracy param
   *        if anal accuracy
   *            die
   *        else
   *            print warning message
   *
   * // date is not later
   * do the conversion
   */
  if (p_utcDate->unixDate.tm_year > 106 ||
      (p_utcDate->unixDate.tm_year == 106 &&
       (p_utcDate->unixDate.tm_mon > LALMONTH_DEC ||
        (p_utcDate->unixDate.tm_mon == LALMONTH_DEC &&
         p_utcDate->unixDate.tm_mday == 31 &&
         p_utcDate->unixDate.tm_hour == 23 &&
         p_utcDate->unixDate.tm_min  == 59 &&
         p_utcDate->unixDate.tm_sec > 59))))
    {
      /* check accuracy param */
      if (*p_accuracy == LALLEAPSEC_STRICT) /* strict accuracy */
        {
          ABORT(status, DATEH_ERANGEGPSABS, DATEH_MSGERANGEGPSABS);
        }
      else if (*p_accuracy == LALLEAPSEC_LOOSE)
        {
          LALWarning(status, "may be missing leap seconds");
        }
      else
        {
          LALWarning(status, "may be missing leap seconds");
        }
    }


  /* start counting from the origin of GPS */
  tmpdate.unixDate.tm_year = 80;
  tmpdate.unixDate.tm_mon  = LALMONTH_JAN;
  tmpdate.unixDate.tm_mday =  6;
  tmpdate.unixDate.tm_hour =  0;
  tmpdate.unixDate.tm_min  =  0;
  tmpdate.unixDate.tm_sec  =  0;
  tmpdate.residualNanoSeconds = 0;

  while (tmpdate.unixDate.tm_year < p_utcDate->unixDate.tm_year)
    {
      ddays += days_in_year(&tmpdate);
      tmpdate.unixDate.tm_year++;
    }
  ddays -= 5; /* 5 days in early Jan 1980 */

  while (tmpdate.unixDate.tm_mon < p_utcDate->unixDate.tm_mon)
    {
      ddays += days_in_month(&tmpdate);
      tmpdate.unixDate.tm_mon++;
    }

  ddays += p_utcDate->unixDate.tm_mday - 1;
  dsecs  = ddays * SECS_PER_DAY;

  dsecs += p_utcDate->unixDate.tm_hour * SECS_PER_HOUR +
    p_utcDate->unixDate.tm_min * SECS_PER_MIN +
    p_utcDate->unixDate.tm_sec;

  /* add in leap seconds */
  i = 9;   /* corresponds to the leap sec data for 1981-Jul-1 */
  while (i < nleaps)
    {
      if (leap_sec_data[i].year < p_utcDate->unixDate.tm_year)
        dsecs++;
      else if (leap_sec_data[i].year == p_utcDate->unixDate.tm_year &&
               leap_sec_data[i].mon <= p_utcDate->unixDate.tm_mon)
        dsecs++;

      ++i;
    }

  p_gpsTime->gpsSeconds = dsecs;
  p_gpsTime->gpsNanoSeconds = p_utcDate->residualNanoSeconds;


  RETURN (status);
} /* END: LALUTCtoGPS() */



/*
 * LALLeapSecs()
 * Compute TAI-UTC for a given GPS time.
 * To reduce complications, we'll only deal with dates since the
 * beginning of GPS time, i.e. 1980-Jan-06  00:00:00 UTC
 * The source code for this function MUST be updated whenever a new leap
 * second is announced.

 ftp://maia.usno.navy.mil/ser7/tai-utc.dat

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
 2006 JAN  1 =JD 2453736.5  TAI-UTC=  33.0       S + (MJD - 41317.) X 0.0      S
 */

typedef struct gps_leap_sec {
  time_t gps;       /* GPS time when leap sec was introduced */
  INT4   tai_utc;   /* TAI-UTC at GPS time gps */
} gps_leap_sec_t;

/* <lalVerbatim file="GPStoUTCCP"> */
void
LALLeapSecs (LALStatus                    *status,
             INT4                         *p_leapSecs,  /* output - GPS-UTC,
                                                           i.e. the number of
                                                           leap seconds introduced
                                                           since the GPS epoch
                                                           1980-Jan-06 */
             const LIGOTimeGPS            *p_gpsTime,   /* input - GPS time */
             const LALLeapSecFormatAndAcc *p_formatAndAcc)  /* format and
                                                               accuracy parameters */
{ /* </lalVerbatim> */

  /* UPDATEME -- to update, add the new (GPS, TAI - UTC) entry at the end
   *             see http://hpiers.obspm.fr/eoppc/bul/bulc/bulletinc.dat */
  /* Table of TAI-UTC */
  static const gps_leap_sec_t gpsLeaps[] =
    {
      {0,         19},  /* 1980-Jan-06 */
      {46828801,  20},  /* 1981-Jul-01 */
      {78364802,  21},  /* 1982-Jul-01 */
      {109900803, 22},  /* 1983-Jul-01 */
      {173059204, 23},  /* 1985-Jul-01 */
      {252028805, 24},  /* 1988-Jan-01 */
      {315187206, 25},  /* 1990-Jan-01 */
      {346723207, 26},  /* 1991-Jan-01 */
      {393984008, 27},  /* 1992-Jul-01 */
      {425520009, 28},  /* 1993-Jul-01 */
      {457056010, 29},  /* 1994-Jul-01 */
      {504489611, 30},  /* 1996-Jan-01 */
      {551750412, 31},  /* 1997-Jul-01 */
      {599184013, 32},  /* 1999-Jan-01 */
      {820108814, 33},  /* 2006-Jan-01 */
    };

  /* number of times leap seconds occur */
  static const INT4   numleaps = sizeof(gpsLeaps)/sizeof(gps_leap_sec_t);
  char   infostr[256];
  INT4   i;


  INITSTATUS (status, "LALLeapSecs", GPSTOUTCC);

  ASSERT (p_leapSecs != (INT4 *)NULL, status,
          DATEH_ENULLOUTPUT, DATEH_MSGENULLOUTPUT);

  ASSERT (p_gpsTime != (LIGOTimeGPS *)NULL, status,
          DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  ASSERT (p_gpsTime->gpsSeconds >= 0, status,
          DATEH_ERANGEGPSABS, DATEH_MSGERANGEGPSABS);

  ASSERT (p_formatAndAcc != (LALLeapSecFormatAndAcc *)NULL, status,
          DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  ASSERT ((p_formatAndAcc->format == LALLEAPSEC_TAIUTC ||
           p_formatAndAcc->format == LALLEAPSEC_GPSUTC),
          status, DATEH_EFORMATPARAMOUTOFRANGE,
          DATEH_MSGEFORMATPARAMOUTOFRANGE);

  ASSERT ((p_formatAndAcc->accuracy == LALLEAPSEC_LOOSE ||
           p_formatAndAcc->accuracy == LALLEAPSEC_STRICT),
          status, DATEH_EACCPARAMOUTOFRANGE,
          DATEH_MSGEACCPARAMOUTOFRANGE);

  XLALPrintDeprecationWarning("LALLeapSecs", "XLALLeapSeconds");

  if (lalDebugLevel & LALINFO)
  {
    snprintf(infostr, INFOSTR_LEN, "Max. tested GPS is %d\n", maxtestedGPS);
    LALInfo(status, infostr);
  }

  /*
   * if GPS is later than maxtestedGPS
   *    check accuracy param
   *        if anal accuracy
   *            die
   *        else
   *            print warning message
   */
  if (p_gpsTime->gpsSeconds > maxtestedGPS)
    {
      /* check accuracy param */
      if (p_formatAndAcc->accuracy == LALLEAPSEC_STRICT)
        {
          ABORT(status, DATEH_ERANGEGPSABS, DATEH_MSGERANGEGPSABS);
        }
      else if (p_formatAndAcc->accuracy == LALLEAPSEC_LOOSE)
        {
          LALWarning(status, "may be missing leap seconds");
        }
      else
        {
          LALWarning(status, "may be missing leap seconds");
        }
    }

  if (p_gpsTime->gpsSeconds == 0)
    {
      *p_leapSecs = gpsLeaps[0].tai_utc;
    }
  else
    {
      i = 1;
      while (p_gpsTime->gpsSeconds > gpsLeaps[i].gps &&
             i < numleaps)
        ++i;

      *p_leapSecs = gpsLeaps[i-1].tai_utc;
    }

  if (lalDebugLevel & LALINFO)
  {
    snprintf(infostr, INFOSTR_LEN, "Format = %d\n", p_formatAndAcc->format);
    LALInfo(status, infostr);
  }

  if (p_formatAndAcc->format == LALLEAPSEC_GPSUTC)
    {
      *p_leapSecs -= 19;
    }

  RETURN(status);
} /* LALLeapSecs */
