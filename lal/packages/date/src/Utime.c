/*----------------------------------------------------------------------- 
 * 
 * File Name: Utime.c 
 * 
 * Author: David Chin <dwchin@umich.edu>
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * Utime
 * 
 * SYNOPSIS 
 *
 * Utime(): Returns UTC time in LALDate structure, corrected for
 *          leap seconds given UTC time in seconds since Unix epoch
 * 
 * DESCRIPTION
 *
 * Utime():
 *      Inputs:   LIGOTimeUnix *unixtime  -- UTC time in seconds since
 *                                             Unix epoch
 *                                             1970-01-01 00:00:00 UTC
 *
 *      Outputs:  LALDate *utc            -- UTC time
 * 
 * 
 * DIAGNOSTICS 
 * (Abnormal termination conditions, error and warning codes summarized 
 * here. More complete descriptions are found in documentation.)
 *
 * CALLS
 * (list of LLAL, LDAS, other non-system functions/procedures called. 
 * 
 * NOTES
 * From http://tycho.usno.navy.mil/leapsec.html

 * Coordinated Universal Time (UTC) is defined by the CCIR Recommendation
 * 460-4 (1986). It differs from TAI by the total number of leap seconds,
 * so that UT1-UTC stays smaller than 0.9s in absolute value.  The decision
 * to introduce a leap second in UTC is the responsibility of the
 * International Earth Rotation Service (IERS). According to the CCIR
 * Recommendation, first preference is given to the opportunities at the
 * end of December and June, and second preference to those at the end of
 * March and September. Since the system was introduced in 1972, only dates
 * in June and December have been used.  TAI is expressed in terms of UTC
 * by the relation TAI = UTC + dAT, where dAT is the total algebraic sum of
 * leap seconds.
 *
 * The first leap second was introduced on June 30, 1972. Information on
 * the most recent leap second can be found here. The historical list of
 * leap seconds lists all of them.
 *
 * > The Global Positioning System (GPS) epoch is January 6, 1980 and is
 * > synchronized to UTC. GPS is NOT adjusted for leap seconds. 
 * > 
 * > As of 1 January 1999,
 * >        TAI is ahead of UTC   by 32 seconds.
 * >        TAI is ahead of GPS   by 19 seconds.
 * >        GPS is ahead of UTC   by 13 seconds.
 *
 * Until 1960, Universal Time (UT) was taken as the independent variable of
 * astronomical ephemerides.  UT was then replaced by Ephemeris Time (ET),
 * based on the motion of the sun.  However, ET did not include
 * relativistic effects, such as corrections for the gravitational
 * potential and velocity, as required by advances in the accuracy of time
 * comparisons.  Thus ET was superseded in 1981 by Terrestrial Dynamical
 * Time (TDT) and Barycentric Dynamical Time (TDB), which distinguish
 * coordinate systems with origins at the center of the Earth and the
 * center of the solar system, respectively, and are consistent with the
 * general theory of relativity.  In the language of general relativity,
 * TDT is a proper time while TDB is a coordinate time.  In 1991, TDT was
 * renamed simply Terrestrial Time (TT) and two additional relativistic
 * time scales, Geocentric Coordinate Time (TCG) and Barycentric Coordinate
 * Time (TCB) were adopted.  Definitions of these time scales are given in
 * Systems of Time.
 * 
 *----------------------------------------------------------------------- */

#include "LALRCSID.h"

NRCSID (UTIMEC, "$Id$");

#include "Date.h"
#include "date_value.h"


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


/*
 * Return UTC time in LALDate structure, taking leap seconds into account.
 * Replaces the C library's gmtime(). Adapted from the GRASP library by
 * Bruce Allen, et al.
 */
void
Utime (Status              *status,
       LALDate             *utc,
       const LIGOTimeUnix  *unixtime)
{
    /* latest time for which this routine will work */
    const INT4   maxtested = 946684823;
    /* number of times leap seconds occur */
    const INT4   numleaps = sizeof(leaps)/sizeof(time_t);
    time_t       tmptime;
    INT4         i;
    LALUnixDate *gmt;
    
    INITSTATUS (status, "Utime", UTIMEC);

    /*
     * Check pointer to input variable
     */
    ASSERT (unixtime != (LIGOTimeUnix *)NULL, status,
            UTIME_ENULLINPUT, UTIME_MSGENULLINPUT);

    /*
     * Check for valid input
     */
    tmptime = unixtime->unixSeconds;
    ASSERT (tmptime >= 0 && tmptime <= maxtested, status,
            UTIME_ERANGE, UTIME_MSGERANGE);
    
    /*
     * Check pointer to output variable.
     */
    ASSERT (utc != (LALDate *)NULL, status,
            UTIME_ENULLOUTPUT, UTIME_MSGENULLOUTPUT);


    /*
     * Check for broken gmtime(). If not broken, use gmtime().
     */
    tmptime = (2441317 - 2440587) * SECS_PER_DAY; /* pick a special time */
    gmt     = gmtime(&tmptime);

    /* If gmtime() is NOT broken, tm_sec field will be 60 since this is a
     * year with a leap second */
    tmptime = unixtime->unixSeconds;
    if (gmt->tm_sec == (time_t)60)
    {
        gmt = gmtime(&tmptime);
        utc->unixDate.tm_sec   = gmt->tm_sec;
        utc->unixDate.tm_min   = gmt->tm_min;
        utc->unixDate.tm_hour  = gmt->tm_hour;
        utc->unixDate.tm_mday  = gmt->tm_mday;
        utc->unixDate.tm_mon   = gmt->tm_mon;
        utc->unixDate.tm_year  = gmt->tm_year;
        utc->unixDate.tm_wday  = gmt->tm_wday;
        utc->unixDate.tm_yday  = gmt->tm_yday;
        utc->unixDate.tm_isdst = gmt->tm_isdst;
        
        utc->residualNanoSeconds = unixtime->unixNanoSeconds;

        RETURN (status);
    }

    /*
     * Check the given time to see if a leap second is in effect
     */
    tmptime = unixtime->unixSeconds;
    i = 0;
    while (i < numleaps && leaps[i] + i - 1 < tmptime)
        ++i;

    if (tmptime == (leaps[i] + i - 1))
    {
        tmptime -= i;
        gmt      = gmtime(&tmptime);
        utc->unixDate.tm_sec   = gmt->tm_sec;
        utc->unixDate.tm_min   = gmt->tm_min;
        utc->unixDate.tm_hour  = gmt->tm_hour;
        utc->unixDate.tm_mday  = gmt->tm_mday;
        utc->unixDate.tm_mon   = gmt->tm_mon;
        utc->unixDate.tm_year  = gmt->tm_year;
        utc->unixDate.tm_wday  = gmt->tm_wday;
        utc->unixDate.tm_yday  = gmt->tm_yday;
        utc->unixDate.tm_isdst = gmt->tm_isdst;

        /* insert leap second */
        utc->unixDate.tm_sec = 60;
    }
    else
    {
        tmptime -= (i - 1);
        gmt      = gmtime(&tmptime);
        utc->unixDate.tm_sec   = gmt->tm_sec;
        utc->unixDate.tm_min   = gmt->tm_min;
        utc->unixDate.tm_hour  = gmt->tm_hour;
        utc->unixDate.tm_mday  = gmt->tm_mday;
        utc->unixDate.tm_mon   = gmt->tm_mon;
        utc->unixDate.tm_year  = gmt->tm_year;
        utc->unixDate.tm_wday  = gmt->tm_wday;
        utc->unixDate.tm_yday  = gmt->tm_yday;
        utc->unixDate.tm_isdst = gmt->tm_isdst;
    }

    /* set residual nanoseconds */
    utc->residualNanoSeconds = unixtime->unixNanoSeconds;

    RETURN (status);
} /* END Utime() */

