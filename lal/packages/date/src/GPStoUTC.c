/* <lalVerbatim file="GPStoUTCCV">
Author: David Chin <dwchin@umich.edu> +1-734-730-1274
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

\subsubsection*{Description}

These routines convert time in GPS seconds and nanoseconds
(\texttt{LIGOTimeGPS}) and time in UTC (\texttt{LALDate}), taking into
account leap seconds until 2002-Mar-31 23:59 UTC.

\subsubsection*{Algorithms}

The conversion from GPS to UTC is copied directly from
GRASP~\cite{grasp:194}.  It does the conversion by counting TAI seconds
starting from the Unix epoch origin, 1970-Jan-01 00:00:00 UTC.  A static
table of leap seconds is compiled in: this \emph{must} be updated whenever
a new leap second is introduced.  The latest leap second included is
1999-Jan-01.

The conversion from UTC to GPS is done by counting the amount of elapsed 
time since the GPS epoch origin, 1980-Jan-06 00:00:00 UTC.  Again, leap
seconds are accounted for by a static table (different from the one used in
GPS to UTC) which \emph{must} be updated whenever a new leap seconds is
introduced.  The latest leap second included is 1999-Jan-01.

\subsubsection*{Uses}

\subsubsection*{Notes}

\textbf{Known deficiencies:} \texttt{LALUTCtoGPS()} is not well-tested for
dates before 1980-Jan-06, although it will produce output (negative values
of \texttt{gpsSeconds},\texttt{gpsNanoSeconds} field set to negative as
well).

\textbf{Example:} To convert a GPS time to UTC, and then back to GPS:

\begin{verbatim}

#include <lal/LALStdlib.h>
#include <lal/Date.h>

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


</lalLaTeX> */

#include <lal/LALRCSID.h>

NRCSID (GPSTOUTCC, "$Id$");

#include <lal/Date.h>
#include "date_value.h"


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
  ((2451179-2440587)*SECS_PER_DAY),
};

/*
 * Convert GPS seconds to UTC date-time contained in LALDate structure
 */
/* <lalVerbatim file="GPStoUTCCP"> */
void
LALGPStoUTC (LALStatus                *status,
             LALDate                  *p_utcDate,
             const LIGOTimeGPS        *p_gpsTime,
             const LALLeapSecAccuracy *p_accuracy)
{ /* </lalVerbatim> */
  time_t unixTime;
  /* latest time for which this routine will work: 2002-Mar-31 23:59:00 */
  /* 24 leap seconds because of the two interpolated ones: 1970-Jan-1 and
   * 1970-Jan-22 */
  /*  const INT4   maxtested = (24*365 + 8*366 + 2*31 + 28)*SECS_PER_DAY -
      60 + 24; */
  const INT4   maxtestedGPS = 701654353;
  /* number of times leap seconds occur */
  const INT4   numleaps = sizeof(leaps)/sizeof(time_t);
  time_t       tmptime;
  LALUnixDate  tmputc;
  char         tmpstamp[32];
  CHAR         infostr[128];
  INT4         i;

  INITSTATUS (status, "LALGPStoUTC", GPSTOUTCC);

  ASSERT (p_gpsTime != (LIGOTimeGPS *)NULL, status,
          DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  ASSERT (p_gpsTime->gpsSeconds >= 0, status,
          DATEH_ERANGEGPSABS, DATEH_MSGERANGEGPSABS);

  ASSERT (p_accuracy != (LALLeapSecAccuracy *)NULL, status,
          DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  ASSERT (p_utcDate != (LALDate *)NULL, status,
          DATEH_ENULLOUTPUT, DATEH_MSGENULLOUTPUT);

  /* we use Unix epoch as our origin */
  unixTime = p_gpsTime->gpsSeconds + UNIXGPS;

  if (lalDebugLevel > 0)
    {
      sprintf(infostr, "Max. tested GPS is %d\n", maxtestedGPS);
      LALInfo(status, infostr);
      
      /* 1998-Dec-31 23:59:60 */
      tmptime = (22*365 + 7*366)* SECS_PER_DAY + 23;
      gmtime_r(&tmptime, &tmputc);

      asctime_r(&tmputc, tmpstamp);
      sprintf(infostr, "tmputc = %s\n", tmpstamp);
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
       * into account in the past (i.e. before the test date */

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

      if (lalDebugLevel > 0)
        sprintf(infostr, "unixTime = %ld; leaps[%d] = %ld", unixTime,
                i, leaps[i]);
      LALInfo(status, infostr);

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
  case 0:
  case 2:
  case 4:
  case 6:
  case 7:
  case 9:
  case 11:
    return 31;
    break;

  case 3:
  case 5:
  case 8:
  case 10:
    return 30;
    break;

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

/*
 * Table of leap seconds
 */
static leap_sec_t leap_sec_data[] =
  {
    {72, 6, 1},
    {73, 0, 1},
    {74, 0, 1},
    {75, 0, 1},
    {76, 0, 1},
    {77, 0, 1},
    {78, 0, 1},
    {79, 0, 1},
    {80, 0, 1},
    {81, 6, 1},
    {82, 6, 1},
    {83, 6, 1},
    {85, 6, 1},
    {88, 0, 1},
    {90, 0, 1},
    {91, 0, 1},
    {92, 6, 1},
    {93, 6, 1},
    {94, 6, 1},
    {96, 0, 1},
    {97, 6, 1},
    {99, 0, 1},
  };

/* <lalVerbatim file="GPStoUTCCP"> */
void
LALUTCtoGPS (LALStatus                *status,
             LIGOTimeGPS              *p_gpsTime,
             const LALDate            *p_utcDate,
             const LALLeapSecAccuracy *p_accuracy)
{ /* </lalVerbatim> */
  int ddays = 0;
  int dsecs = 0;
  LALDate tmpdate;
  static LALDate gpsref;
  int i = 0;
  char infostr[256];
  static const int nleaps = sizeof(leap_sec_data)/sizeof(leap_sec_t);

  /* When GPS began */
  gpsref.unixDate.tm_sec  = 0;
  gpsref.unixDate.tm_min  = 0;
  gpsref.unixDate.tm_hour = 0;
  gpsref.unixDate.tm_mday = 6;
  gpsref.unixDate.tm_mon  = 0;
  gpsref.unixDate.tm_year = 80;
  gpsref.unixDate.tm_wday = 0;
  gpsref.unixDate.tm_yday = 0;

  if (lalDebugLevel >= 8)
    {
      sprintf(infostr, "Date given: %d-%d-%d %d:%d:%d %d\n",
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

  ASSERT (p_utcDate != (LALDate *)NULL, status,
          DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  /* Can't convert dates before 1980-Jan-06 */
  ASSERT (p_utcDate->unixDate.tm_year > 80 ||
          (p_utcDate->unixDate.tm_year == 80 &&
           (p_utcDate->unixDate.tm_mon > 0 ||
            (p_utcDate->unixDate.tm_mon == 0 &&
             p_utcDate->unixDate.tm_mday >= 6))), status,
          DATEH_EGPSDATETOOEARLY, DATEH_MSGEGPSDATETOOEARLY);


  LALInfo(status, ">= 1980-01-06 only");

  /*
   * Check that time asked for is not after last known leap sec
   * Use by: 2002-Mar-31 23:59:00
   *
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
  if (p_utcDate->unixDate.tm_year > 102 ||
      (p_utcDate->unixDate.tm_year == 102 &&
       (p_utcDate->unixDate.tm_mon > 2 ||
        (p_utcDate->unixDate.tm_mon == 2 &&
         p_utcDate->unixDate.tm_mday == 31 &&
         p_utcDate->unixDate.tm_hour == 23 &&
         p_utcDate->unixDate.tm_min  == 59 &&
         p_utcDate->unixDate.tm_sec > 0))))
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
  
  
  /* start counting from the origin */
  tmpdate.unixDate.tm_year = 80;
  tmpdate.unixDate.tm_mon  =  0;
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
