/*----------------------------------------------------------------------- 
 * 
 * File Name: SecsToLALDate.c 
 * 
 * Author: David Chin <dwchin@umich.edu>
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
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
 * To convert a LIGOTimeUTC structure to LALDate, use LALUtime().
 * 
 *----------------------------------------------------------------------- */

#include "LALRCSID.h"

NRCSID (SECSTOLALDATEC, "$Id$");

#include <math.h>
#include "Date.h"
#include "date_value.h"

/*
 * Convert time in seconds to LALDate struct
 */
void
LALSecsToLALDate(LALStatus  *status,
              LALDate *date,
              REAL8    seconds)
{
    REAL8 hr, min, sec;
    REAL8 tmpsec;
    REAL8 dum;

    INITSTATUS (status, "SECSTOLALDATE", SECSTOLALDATEC);

    /*
     * Check pointer to output variable
     */
    ASSERT (date != (LALDate *)NULL, status,
            SECSTOLALDATE_ENULLOUTPUT, SECSTOLALDATE_MSGENULLOUTPUT);

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
     * Set the LALDate structure. Only the time fields matter, of course.
     * Fractional seconds become "residual".
     */
    date->unixDate.tm_hour = (INT4)hr;
    date->unixDate.tm_min  = (INT4)min;
    date->unixDate.tm_sec  = (INT4)sec;
    date->residualNanoSeconds = modf(sec, &dum) * 1.e+09;

    RETURN (status);
} /* END LALSecsToLALDate() */
