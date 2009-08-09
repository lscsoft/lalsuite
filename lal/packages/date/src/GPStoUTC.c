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
\idx{LALLeapSecs()}

\subsubsection*{Description}

\texttt{LALLeapSecs()} returns the number of leap seconds introduced since
the GPS epoch 1980-Jan-06, abbreviated GPS-UTC.


\subsubsection*{Algorithms}

\subsubsection*{Uses}

\subsubsection*{Notes}

These routines will not work for times before 1980-01-06 00:00:00 UTC (GPS
0).  The latest leap second that can be accounted for is the one added at
the end of 2005-Dec. % UPDATEME %  These routines have accurate leap second
information until 2006-Dec-31. % UPDATEME %

</lalLaTeX> */

#include <lal/LALRCSID.h>

NRCSID (GPSTOUTCC, "$Id$");

#include <lal/LALStdio.h>
#include <lal/Date.h>
#include "date_value.h"
#include <lal/XLALError.h>

#define INFOSTR_LEN 256


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
