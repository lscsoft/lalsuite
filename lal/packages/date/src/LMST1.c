/* <lalVerbatim file="LMST1CV">
Author: David Chin <dwchin@umich.edu> +1-734-730-1274
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
\index{\texttt{LALGMST1()}}
\index{\texttt{LALGPStoGMST1()}}
\index{\texttt{LALLMST1()}}
\index{\texttt{LALGPStoLMST1()}}


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


All the routines will output GMST1 or LMST1 in the units specified by
the \texttt{outunits} argument.

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
Universal Time elapsed since JD 2451545.0 UT1 (J2000.0), taking on
values of $\pm 0.5, \pm 1.5, \pm 2.5,
\pm 3.5, \ldots$



\subsubsection*{Uses}



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
    date.unixDate.tm_mon  = 10;
    date.unixDate.tm_year = 94;

    longitude = 0.; 
    GMST1(&status, &gmsthours, &date, MST_HRS);
    LMST1(&status, &lmsthours, &date, longitude, MST_HRS);

    GMST1(&status, &gmstsecs, &date, MST_SEC);
    SecsToLALDate(&status, &mstdate, gmstsecs);
    strftime(timestamp, 64, "%Hh %Mm %S", &(mstdate.unixDate));
    sprintf(tmpstr, "%fs", mstdate.residualNanoSeconds * 1.e-9);
    strcat(timestamp, tmpstr+1);
    
    printf("gmsthours = %f = %s\n", gmsthours, timestamp);
    printf("    expect: 3.655728 = 03h 39m 20.6222s \n");

    return(0);
}
\end{verbatim}

\subsubsection*{Notes}

From~\cite{esaa:1992}:

\begin{quote}
  Universal Time (UT) is the measure of time used as the basis for all
  civil time-keeping.  It conforms closely to the mean diurnal motion
  of the Sun.  The apparent diurnal motion of the Sun involves both
  the nonuniform diurnal rotation of the Earth and th motion of the
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
#include <lal/Date.h>
#include "date_value.h"


/*
 * Compute GMST1 in requested units given date UTC.
 * Use algorithm as coded in NOVAS-C Version 2.0.1 (10 Dec 99)
 *   Naval Observatory Vector Astrometry Subroutines
 *  C Version
 */

/* <lalVerbatim file="LMST1CP"> */
void
LALGMST1 (LALStatus     *status,
          REAL8         *p_gmst,     /* output - GMST1 */
          const LALDate *p_date,     /* input  - date and time */
          LALMSTUnits    outunits)   /* GMST1 units */
{ /* </lalVerbatim> */
  REAL8 jdate;
  REAL8 jd_hi, jd_lo;
  INT4  tmp;
  REAL8 st;
  REAL8 tu;
  REAL8 tu2;
  const REAL8      a =         -6.2e-6;
  const REAL8      b =          0.093104;
  const REAL8      c =      67310.54841;
  const REAL8      d =    8640184.812866;
  const REAL8      e = 3155760000.0;

  INITSTATUS (status, "LALGMST1", LMST1C);

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
   * Compute GMST1 for 0h UT1 on given date in seconds since J2000.0 
   */
  LALJulianDate(status, &jdate, p_date);
  jdate -= J2000_0;
  tmp    = (INT4)jdate;  /* get the integral part of jdate */
  jd_hi  = (REAL8)tmp;
  jd_lo  = jdate - jd_hi;

  tu     = jdate / JDAYS_PER_CENT;
  tu2    = tu * tu;

  jd_hi /= JDAYS_PER_CENT;
  jd_lo /= JDAYS_PER_CENT;
    
  st = a*tu2*tu + b*tu2 + c + d*jd_lo + e*jd_lo
    + d*jd_hi + e*jd_hi;

  *p_gmst = fmod(st, (REAL8)SECS_PER_DAY);

  if (*p_gmst < 0.)
    *p_gmst += (REAL8)SECS_PER_DAY;

  /*
   * Convert GMST to requested units
   */
  switch (outunits)
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
               LALMSTUnits        outunits) /* GMST1 units */
{ /* </lalVerbatim> */
  LIGOTimeUnix    unixTime;
  LALDate         date;
    
  INITSTATUS (status, "LALGPStoGMST1", LMST1C);

    
  /*
   * Check pointer to input variables
   */
  ASSERT (p_gps != (LIGOTimeGPS *)NULL, status,
          DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  /*
   * Check pointer to output variable
   */
  ASSERT (p_gmst != (REAL8 *)NULL, status,
          DATEH_ENULLOUTPUT, DATEH_MSGENULLOUTPUT);

  /*
   * Convert GPS to date-time structure
   */

  /* first, GPS to Unix */
  LALGPStoU(status, &unixTime, p_gps);

  /* then, Unix to date-time */
  LALUtime(status, &date, &unixTime);

  LALGMST1(status, p_gmst, &date, outunits);

  return;
} /* END: LALGPStoGMST1() */





/*
 * Compute LMST1 in requested units given date-time
 */
/* <lalVerbatim file="LMST1CP"> */
void
LALLMST1 (LALStatus             *status,
          REAL8                 *p_lmst,            /* output - LMST1 */
          const LALPlaceAndDate *p_place_and_date,  /* input - location and date */
          LALMSTUnits            outunits)          /* LMST1 units */
{ /* </lalVerbatim> */
  REAL8 gmst;
  REAL8 day = 0;
  REAL8 longitude =
    p_place_and_date->p_detector->frDetector.vertexLongitudeDegrees;


  INITSTATUS (status, "LALLMST1", LMST1C);

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
  LALGMST1(status, &gmst, p_place_and_date->p_date, outunits);

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

  RETURN (status);
} /* END LALLMST1() */



/*
 * Compute LMST1 in requested units given GPS time
 */
/* <lalVerbatim file="LMST1CP"> */
void
LALGPStoLMST1( LALStatus             *status,
               REAL8                 *p_lmst,          /* output - LMST1 */
               const LALPlaceAndGPS  *p_place_and_gps, /* input - location and GPS */
               LALMSTUnits            outunits)        /* LMST1 units */
{ /* </lalVerbatim> */
  LIGOTimeUnix    unixTime;
  LALDate         date;
  LALPlaceAndDate place_and_date;


  INITSTATUS (status, "LALGPStoLMST1", LMST1C);

    
  /*
   * Check pointer to input variables
   */
  ASSERT (p_place_and_gps != (LALPlaceAndGPS *)NULL, status,
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
  LALGPStoU(status, &unixTime, p_place_and_gps->p_gps);

  /* then, Unix to date-time */
  LALUtime(status, &date, &unixTime);

  /* stuff it all into a LALPlaceAndDate */
  place_and_date.p_detector = p_place_and_gps->p_detector;
  place_and_date.p_date     = &date;

  LALLMST1(status, p_lmst, &place_and_date, outunits);

  return;
} /* END: LALGPStoLMST1() */


              
