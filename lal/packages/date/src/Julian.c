/*----------------------------------------------------------------------- 
 * 
 * File Name: Julian.c 
 * 
 * Author: David Chin <dwchin@umich.edu>
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * JulianDay
 * ModJulianDay
 * JulianDate
 * ModJulianDate
 * 
 * SYNOPSIS 
 *     LALJulianDay(): Returns Julian Day Number for a given Gregorian date UTC
 *     LALModJulianDay(): Returns Modified Julian Day Number for given
 *                     Gregorian date UTC, i.e. Julian Day - 2400000.5
 *     LALJulianDate(): Returns Julian Date in decimal days
 *     LALModJulianDate(): Returns Modified Julian Date in decimal days
 *
 *
 * DESCRIPTION
 * LALJulianDay()
 *       Inputs:  LALDate *date  -- Gregorian date (UTC)
 *
 *       Outputs: INT4 *jDay     -- Julian Day number
 *
 * LALModJulianDay()
 *       Inputs:  LALDate *date   -- Gregorian date (UTC)
 *
 *       Outputs: REAL8 *modJDay  -- modified Julian day
 *
 * LALJulianDate()
 *       Inputs:  LALDate *date -- Gregorian date (UTC)
 *
 *       Outputs: REAL8 *jDate  -- Julian Date in decimal days
 *
 * LALModJulianDate()
 *       Inputs:  LALDate *date   -- Gregorian date (UTC)
 *
 *       Outputs: REAL8 *modJDate -- Modified Julian Date in decimal days 
 * 
 * DIAGNOSTICS 
 * (Abnormal termination conditions, error and warning codes summarized 
 * here. More complete descriptions are found in documentation.)
 *
 * CALLS
 * LALJulianDay(): (none)
 * LALModJulianDay(): LALJulianDay()
 * LALJulianDate(): LALJulianDay()
 * LALModJulianDate(): LALJulianDate()
 * 
 * NOTES
 * (Other notes)
 * See the _Explanatory_Supplement_to_the_Astronomical_Almanac_ for more
 * details.
 *
 * From Green's _Spherical_Astronomy_:
 * 
 *            Mean Julian Year == 365.25 days
 *
 *                Julian Epoch == 1 Jan 4713BCE, 12:00 GMT
 *                                i.e. 4713 BC Jan 01d.5 GMT
 *
 *  Fundamental Epoch: J2000.0 == 2000-01-01.5 TDB
 *
 *  Julian Date is time elapsed since Julian Epoch, measured in days and
 *  fractions of a day
 *  Complications: Tropical Year == 365.2422 days
 *                => Gregorian correction: 1582-10-5 through 1582-10-14
 *                                         were eliminated (10 days)
 *
 *                 Leap years: years ending with two zeroes (e.g. 1700,
 *                 1800) are leap only if divisible by 400. So, 400 civil
 *                 years contains (400 * 365.25) - 3 = 146097 days.
 *
 *  So, the Julian Date of J2000.0 == JD 2451545.0
 *
 *  and, Julian Epoch = J2000.0 + (JD - 2451545)/365.25
 *                      (i.e. no. of years elapsed since J2000.0)
 *
 *
 * One algorithm for computing the Julian day is from
 * Van Flandern, T.C., and Pulkkinen, K.F., Astrophysical Journal
 * Supplement Series, 41, 391-411, 1979 Nov.
 * based on a formula in _Explanatory_Supplement_to_the_Astronomical_Almanac_,
 * 1992, Chapter 12 where the algorithm is due to Fliegen and Van Flandern, "A
 * Machine Algorithm for the Processing of Calendar Dates", Communications
 * of the ACM, 11, 657 (1968), and "compactified" by Muller, P.M., and
 * Wimberly, R.N.
 *
 * jd = 367 * y - 7 * (y + (m + 9)/12)/4 - 3 * ((y + (m - 9)/7)/100 + 1)/4
 *      + 275 * m/9 + d + 1721029
 *
 * This is valid only for JD >= 0, i.e. after -4713 Nov 23 == 4712 BCE Nov 23
 *
 * A shorter formula from the same reference, but which only works for
 * dates since 1900-March
 *
 * jd = 367 * y - 7 * (y + (m + 9)/12)/4 + 275 * m/9 + d + 1721014
 *
 * We will use the one for dates > 1900-Mar since there is no LIGO data
 * before 1900-Mar.
 * 
 *----------------------------------------------------------------------- */

#include <lal/LALRCSID.h>

NRCSID (JULIANC, "$Id$");

#include <lal/Date.h>
#include "date_value.h"

/*
 * Compute Julian Day for given Gregorian date
 */
void
LALJulianDay (LALStatus     *status,
              INT4          *jDay,
              const LALDate *date)
{
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
    y = (date->unixDate).tm_year + 1900;
    m = (date->unixDate).tm_mon + 1;

    /* Julian Day begins at noon, so fix day of month if necessary */
    if ((date->unixDate).tm_hour < 12)
        d = (date->unixDate).tm_mday - 1;
    else
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
 * Compute Modified Julian Day for given Gregorian date
 */
void
LALModJulianDay (LALStatus        *status,
              REAL8         *modJDay,
              const LALDate *date)
{
    INT4 jd;

    INITSTATUS (status, "LALModJulianDay", JULIANC);
    
    /*
     * Check pointer to input variable
     */
    ASSERT (date != (LALDate *)NULL, status,
            DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

    /*
     * Check pointer to output variable:
     */
    ASSERT (modJDay   != (REAL8 *)NULL, status,
DATEH_ENULLOUTPUT, DATEH_MSGENULLOUTPUT);

    LALJulianDay(status, &jd, date);

    *modJDay = (REAL8)jd - MJDREF;

    RETURN (status);
} /* END LALModJulianDay() */



/*
 * Compute Julian Date for given Gregorian date and UTC time
 */
void
LALJulianDate (LALStatus        *status,
            REAL8         *jDateOut,
            const LALDate *date)
{
    INT4  hr, min, sec;
    REAL8 rns;          /* residual nanoseconds */
    INT4  jday;
    REAL8 jdate;

    INITSTATUS(status, "LALJulianDate", JULIANC);

    /*
     * Check pointer to input variable
     */
    ASSERT (date != (LALDate *)NULL, status,
            DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

    /*
     * Check pointer to output variable:
     */
    ASSERT (jDateOut != (REAL8 *)NULL, status,
DATEH_ENULLOUTPUT, DATEH_MSGENULLOUTPUT);

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
    LALJulianDay(status, &jday, date);

    /*
     * Convert to fractions of a day
     */
    jdate  = (REAL8)jday + 0.5;
    jdate += (REAL8)hr/24.;
    jdate += (REAL8)min/1440.;
    jdate += ((REAL8)sec + rns)/86400.;

    *jDateOut = jdate;

    RETURN (status);
} /* END LALJulianDate() */



/*
 * Compute Modified Julian Date for given Gregorian date and UTC time
 */
void
LALModJulianDate (LALStatus     *status,
                  REAL8         *modJDate,
                  const LALDate *date)
{
    INT4  hr, min, sec;
    REAL8 jdate;

    INITSTATUS(status, "LALModJulianDate", JULIANC);

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
     * Extract Hour, Minute, and Second
     */
    hr  = (date->unixDate).tm_hour;
    min = (date->unixDate).tm_min;
    sec = (date->unixDate).tm_sec;

    /*
     * Get Julian Date, and modify it
     */
    LALJulianDate(status, &jdate, date);
    jdate -= MJDREF;
    
    *modJDate = jdate;

    RETURN (status);
} /* END LALModJulianDate() */


