/*
*  Copyright (C) 2007 Duncan Brown, David Chin, Jolien Creighton, Kipp Cannon, Peter Shawhan, Reinhard Prix
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

/* <lalVerbatim file="LMST1CV">
Author: David Chin <dwchin@umich.edu> +1-734-709-9119
$Id$
</lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{LMST1.c}}
\label{ss:LMST1.c}

Routines related to Local and Greenwich Mean Sidereal Time (LMST1 and
GMST1) computations.

\subsection*{Prototypes}
\vspace{0.1in}

\input{LMST1CP}
\idx{LALGMST1()}
\idx{LALGPStoGMST1()}
\idx{LALLMST1()}
\idx{LALGPStoLMST1()}


\subsubsection*{Description}

The routines in this module compute Mean Sidereal Time in a choice of
units: seconds, hours, degrees, or radians. LMST1 is offset from GMST1
by the longitude of the observing post.

\begin{itemize}
\item \texttt{LALGMST1()} computes GMST1 given a Gregorian date UTC in an
\texttt{LALDate} structure.
\item  \texttt{LALGPStoGMST1()} computes GMST1 given a GPS time.
\item  \texttt{LALLMST1()} computes LMST1 given an observing location (in a
  \texttt{LALDetector} structure) and a Gregorian date UTC.
\item \texttt{LALGPStoLMST1()} computes LMST1 given an observing location
  and a GPS time.
\end{itemize}


All the routines will output GMST1 or LMST1 in the units and leap second
accuracy specified by the \texttt{pUnitsAndAcc} argument.  The sidereal
time produced is within $\tilde 1$ sidereal second of values published in
the Almanac.

\subsubsection*{Algorithms}

The basic definitions and formulae are from~\cite{esaa:1992}, Ch. 2,
Sec. 24, and also Sec. B of the Astronomical Almanac. The formula
computes GMST1 for 0h UT1.  To compute GMST1 for other times, a simple
linear interpolation is done. The implementation used in
\texttt{LALGMST1()} is from~\cite{novas:1999} .  Since 1984-Jan-01,
GMST has been related to UT1 as follows:
%
\begin{displaymath}
  \mathrm{GMST\>of}\>0^{h}\mathrm{UT1} = 24110^{s}.54841 +
  8640184^{s}.812866\,T_{u} + 0^{s}.093104\,T^{2}_{u} -
  6.2\times10^{-6}\,T^{3}_{u}
\end{displaymath}
%
where $T_{u} = d_{u} / 36525$, $d_{u}$ is the number of days of
Universal Time elapsed since JD 2451545.0 UT1 (January 1, 2000, at 12:00 UT1),
taking on values of $\pm 0.5, \pm 1.5, \pm 2.5,
\pm 3.5, \ldots$



\subsubsection*{Uses}

\subsubsection*{Notes}

Here is a simple example:

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
    REAL8   gmsthours, lmsthours;
    REAL8   gmstsecs;
    REAL8   longitude;
    LALMSTUnitsAndAcc units_and_acc;
    time_t  timer;
    CHAR    timestamp[64], tmpstr[64];
    Status  status = {0};

    if (argc == 1)
    {
        printf("Usage: TestUTCtoGPS debug_level -- debug_level = [0,1,2]\n");
        return 0;
    }

    if (argc == 2)
        debuglevel = atoi(argv[1]);

    INITSTATUS(&status, "TestLMST", TESTLMSTC);

    printf("TEST of GMST1 routine\n");
    printf("=====================\n");


    // Check against the Astronomical Almanac:
    // For 1994-11-16 0h UT - Julian Date 2449672.5, GMST 03h 39m 21.2738s
    date.unixDate.tm_sec  = 0;
    date.unixDate.tm_min  = 0;
    date.unixDate.tm_hour = 0;
    date.unixDate.tm_mday = 16;
    date.unixDate.tm_mon  = LALMONTH_NOV;
    date.unixDate.tm_year = 94;

    longitude = 0.;
    LALGMST1(&status, &gmsthours, &date, MST_HRS);
    LALLMST1(&status, &lmsthours, &date, longitude, MST_HRS);

    LALGMST1(&status, &gmstsecs, &date, MST_SEC);
    LALSecsToLALDate(&status, &mstdate, gmstsecs);
    strftime(timestamp, 64, "%Hh %Mm %S", &(mstdate.unixDate));
    sprintf(tmpstr, "%fs", mstdate.residualNanoSeconds * 1.e-9);
    strcat(timestamp, tmpstr+1);

    printf("gmsthours = %f = %s\n", gmsthours, timestamp);
    printf("    expect: 3.655728 = 03h 39m 20.6222s \n");

    return(0);
}
\end{verbatim}


From~\cite{esaa:1992}:

\begin{quote}
  Universal Time (UT) is the measure of time used as the basis for all
  civil time-keeping.  It conforms closely to the mean diurnal motion
  of the Sun.  The apparent diurnal motion of the Sun involves both
  the nonuniform diurnal rotation of the Earth and the motion of the
  Earth in its orbit around the Sun.  Although it would be possible to
  define a system of time measurement in terms of the hour angle of
  the Sun, such a system could never be related precisely to sidereal
  time and could not, therefore, be determined by observations of star
  transits.  As a result, Universal Time is directly related to
  sidereal time by means of a numerical formula.  It does not refer to
  the motion of the Earth and is not precisely related to the hour
  angle of the Sun.

  Universal Time at any instant can be derived from observations of
  the diurnal motion of the stars or radio sources.  The uncorrected
  observed rotational timescale, which is dependent on the place of
  observation, is designated UT0.  Correcting this timescale for the
  shift in longitude of the observing station caused by polar motion
  produces the UT1 timescale, which is independent of observing
  location, but is influenced by the slightly variable rotation of the
  Earth.
\end{quote}

A useful resource is \url{http://aa.usno.navy.mil/}.

</lalLaTeX> */


/*******************************************************************
 *
 * SYNOPSIS
 *
 * LALGMST1(): Returns LALGMST1 for given date-time
 * LALLMST1(): Returns LALLMST1 given date-time, and longitude
 *
 * DESCRIPTION
 *
 * LALGMST1():
 *      Inputs: LALDate *date      -- date-time to compute GMST1
 *              LALMSTUnits outunits -- units requested:
 *                                      MST_SEC - seconds
 *                                      MST_HRS - hours
 *                                      MST_DEG - degrees
 *                                      MST_RAD - radians
 *
 *      Outputs: REAL8  *gmst      -- LALGMST1 in units requested
 *
 * LALLMST1():
 *      Inputs: LALPlaceAndDate *place_and_date -- location and date-time
 *                                                 to compute GMST1
 *              LALMSTUnits outunits -- units requested:
 *                                      MST_SEC - seconds
 *                                      MST_HRS - hours
 *                                      MST_DEG - degrees
 *                                      MST_RAD - radians
 *
 *      Outputs: REAL8  *lmst      -- LALGMST1 in units requested
 *
 *
 *----------------------------------------------------------------------- */

#include <lal/LALRCSID.h>

NRCSID (LMST1C, "$Id$");

#include <math.h>
#include <lal/LALStdio.h>
#include <lal/Date.h>
#include "date_value.h"
#include <lal/XLALError.h>

#define INFOSTR_LEN 256

/*
 * Compute GMST in requested units given date UTC.
 */

/* <lalVerbatim file="LMST1CP"> */
void
LALGMST1 (LALStatus     *status,
          REAL8         *p_gmst,     /* output - GMST1 */
          const LALDate *p_date,     /* input  - date and time */
          LALMSTUnits    outunits)   /* GMST1 units */
{ /* </lalVerbatim> */
  LIGOTimeGPS gps;
  LALMSTUnitsAndAcc UnitsAndAcc;
  char infostr[INFOSTR_LEN];

  INITSTATUS (status, "LALGMST1", LMST1C);
  ATTATCHSTATUSPTR(status);

  /*
   * Check pointer to input variables
   */
  ASSERT (p_date != (LALDate *)NULL, status,
          DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  /*
   * Check pointer to output variable
   */
  ASSERT (p_gmst != (REAL8 *)NULL, status,
          DATEH_ENULLOUTPUT, DATEH_MSGENULLOUTPUT);

  /*
   * Compute GMST for UTC on given date in seconds
   */

  UnitsAndAcc.accuracy = LALLEAPSEC_STRICT;
  UnitsAndAcc.units = outunits;

  TRY( LALUTCtoGPS( status->statusPtr, &gps, p_date, &UnitsAndAcc.accuracy ), status );
  TRY( LALGPStoGMST1( status->statusPtr, p_gmst, &gps, &UnitsAndAcc ), status );

  if (lalDebugLevel & LALINFO)
  {
    LALSnprintf(infostr, INFOSTR_LEN, "LALGMST1: *p_gmst = %g\n", *p_gmst);
    LALInfo(status, infostr);
  }

  DETATCHSTATUSPTR(status);
  RETURN (status);
} /* END LALGMST1() */



/*
 * Compute GMST1 in requested units given GPS time
 */
/* <lalVerbatim file="LMST1CP"> */
void
LALGPStoGMST1( LALStatus         *status,
               REAL8             *p_gmst,   /* output - GMST1 */
               const LIGOTimeGPS *p_gps,    /* input - GPS time */
               const LALMSTUnitsAndAcc *pUnitsAndAcc) /* GMST1 units and accuracy */
{ /* </lalVerbatim> */
  INITSTATUS (status, "LALGPStoGMST1", LMST1C);
  ATTATCHSTATUSPTR(status);

  /*
   * Check pointer to input variables
   */
  ASSERT (p_gps != (LIGOTimeGPS *)NULL, status,
          DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);
  ASSERT (pUnitsAndAcc != (LALMSTUnitsAndAcc *)NULL, status,
          DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  /*
   * Check pointer to output variable
   */
  ASSERT (p_gmst != (REAL8 *)NULL, status,
          DATEH_ENULLOUTPUT, DATEH_MSGENULLOUTPUT);

  XLALPrintDeprecationWarning("LALGPStoGMST1", "XLALGreenwichMeanSiderealTime");

  /*
   * Compute GMST for GPS on given date in seconds
   */
  *p_gmst = fmod(XLALGreenwichMeanSiderealTime(p_gps) / (2.0 * LAL_PI) * SECS_PER_DAY, SECS_PER_DAY);
  if(*p_gmst < 0.0)
    *p_gmst += SECS_PER_DAY;

  /*
   * Convert GMST to requested units
   */
  switch (pUnitsAndAcc->units)
    {
    case MST_SEC:
      break;
    case MST_HRS:
      *p_gmst /= (REAL8)SECS_PER_HOUR;
      break;
    case MST_DEG:
      *p_gmst /= (REAL8)SECS_PER_HOUR / (REAL8)DEGS_PER_HOUR;
      break;
    case MST_RAD:
      *p_gmst /= (REAL8)SECS_PER_HOUR / (REAL8)DEGS_PER_HOUR * 180. /
        (REAL8)LAL_PI;
      break;
    default:
      break;
    }

  DETATCHSTATUSPTR(status);
  RETURN(status);
} /* END: LALGPStoGMST1() */





/*
 * Compute LMST1 in requested units given date-time
 */
/* <lalVerbatim file="LMST1CP"> */
void
LALLMST1 (LALStatus             *status,
          REAL8                 *p_lmst,            /* output - LMST1 */
          const LALPlaceAndDate *p_place_and_date,  /* input - place and date */
          LALMSTUnits            outunits)          /* LMST1 units */
{ /* </lalVerbatim> */
  REAL8 gmst;
  REAL8 day = 0;
  REAL8 longitude = LAL_180_PI *
    atan2(p_place_and_date->p_detector->location[1],
          p_place_and_date->p_detector->location[0]);

  INITSTATUS (status, "LALLMST1", LMST1C);
  ATTATCHSTATUSPTR(status);

  /*
   * Check pointer to input variables
   */
  ASSERT (p_place_and_date != (LALPlaceAndDate *)NULL, status,
          DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  /*
   * Check pointer to output variable
   */
  ASSERT (p_lmst != (REAL8 *)NULL, status,
          DATEH_ENULLOUTPUT, DATEH_MSGENULLOUTPUT);

  /*
   * Compute LMST1 in seconds since J2000.0
   */

  /* get GMST1 in seconds */
  TRY( LALGMST1(status->statusPtr, &gmst, p_place_and_date->p_date, outunits),
       status);

  /* convert longitude to appropriate units of sidereal time */
  switch (outunits)
        {
    case MST_SEC:
      longitude /= (REAL8)DEGS_PER_HOUR / (REAL8)SECS_PER_HOUR;
      day        = SECS_PER_DAY;
      break;
    case MST_HRS:
      longitude /= (REAL8)DEGS_PER_HOUR;
      day        = 24.;
      break;
    case MST_DEG:
      day        = 360.;
      break;
    case MST_RAD:
      longitude /= 180. / (REAL8)LAL_PI;
      day        = (REAL8)LAL_TWOPI;
      break;
    default:
      break;
    }

  *p_lmst = gmst + longitude;

  while (*p_lmst < 0)
    *p_lmst += day;

  DETATCHSTATUSPTR(status);
  RETURN (status);
} /* END LALLMST1() */



/*
 * Compute LMST1 in requested units given GPS time
 */
/* <lalVerbatim file="LMST1CP"> */
void
LALGPStoLMST1( LALStatus             *status,
               REAL8                 *p_lmst,          /* output - LMST1 */
               const LALPlaceAndGPS  *p_place_and_gps, /* input - place and GPS */
               const LALMSTUnitsAndAcc *pUnitsAndAcc)        /* LMST1 units and accuracy */
{ /* </lalVerbatim> */
  LALDate         date;
  LALPlaceAndDate place_and_date;

  INITSTATUS (status, "LALGPStoLMST1", LMST1C);
  ATTATCHSTATUSPTR(status);


  /*
   * Check pointer to input variables
   */
  ASSERT (p_place_and_gps != (LALPlaceAndGPS *)NULL, status,
          DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  ASSERT (pUnitsAndAcc != (LALMSTUnitsAndAcc *)NULL, status,
          DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  /*
   * Check pointer to output variable
   */
  ASSERT (p_lmst != (REAL8 *)NULL, status,
          DATEH_ENULLOUTPUT, DATEH_MSGENULLOUTPUT);


  /*
   * Convert GPS to date-time structure
   */

  /* first, GPS to Unix */
  TRY( LALGPStoUTC(status->statusPtr, &date, p_place_and_gps->p_gps,
                   &(pUnitsAndAcc->accuracy)),
       status );

  /* stuff it all into a LALPlaceAndDate */
  place_and_date.p_detector = p_place_and_gps->p_detector;
  place_and_date.p_date     = &date;

  TRY( LALLMST1(status->statusPtr, p_lmst, &place_and_date,
                pUnitsAndAcc->units), status );

  DETATCHSTATUSPTR(status);
  RETURN(status);
} /* END: LALGPStoLMST1() */
