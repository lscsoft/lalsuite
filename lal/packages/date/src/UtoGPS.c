/*----------------------------------------------------------------------- 
 * 
 * File Name: GPStoU.c 
 * 
 * Author: David Chin <dwchin@umich.edu>
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * UtoGPS
 * 
 * SYNOPSIS 
 *
 * LALUtoGPS(): Returns UTC in GPS seconds given time in Unix seconds.
 * LALGPStoU(): Returns UTC in Unix seconds given time in GPS seconds.
 * 
 * DESCRIPTION
 *
 * LALUtoGPS():
 *      Inputs:   LIGOTimeUnix *utctime -- UTC time in seconds since Unix epoch
 *
 *      Outputs:  LIGOTimeGPS *gpstime -- UTC time in seconds since GPS epoch
 *
 * LALGPStoU():
 *      Inputs:  LIGOTimeGPS *gpstime -- UTC time in seconds since GPS epoch
 *
 *      Outputs:  LIGOTimeUnix *utctime -- UTC time in seconds since Unix epoch
 *
 * 
 * DIAGNOSTICS 
 * (Abnormal termination conditions, error and warning codes summarized 
 * here. More complete descriptions are found in documentation.)
 *
 * CALLS
 *
 * LALUtoGPS(): none
 * LALGPStoU(): none
 * 
 * NOTES
 * The only difference between LIGOTimeUnix and LIGOTimeGPS is the reference
 * epoch: LIGOTimeUnix uses the Unix epoch 1970-01-01 00:00:00UTC, GPS uses
 * 1980-01-06 00:00:00UTC.
 * US Naval Observatory: GPS is NOT adjusted for leap seconds.  As of
 * 1999-01-01:
 *       TAI is ahead of UTC   by 32 seconds
 *       TAI is ahead of GPS   by 19 seconds
 *       GPS is ahead of UTC   by 13 seconds
 *
 * For more information, see http://tycho.usno.navy.mil/leapsec.html.
 * See also: comments in UTCtime.c
 * 
 *-----------------------------------------------------------------------
 */

#include "LALRCSID.h"

NRCSID (UTOGPSC, "$Id$");

#include "Date.h"
#include "date_value.h"


/*
 * Convert UTC seconds to GPS seconds
 *
 * Input:
 *
 * Output:
 */
void
LALUtoGPS (LALStatus             *status,
        LIGOTimeGPS        *gpstime,
        const LIGOTimeUnix *unixtime)
{
    INITSTATUS (status, "LALUtoGPS", UTOGPSC);

    /*
     * Check pointer to input variable
     */
    ASSERT (unixtime != (LIGOTimeUnix *)NULL, status,
            UTOGPS_ENULLINPUT, UTOGPS_MSGENULLINPUT);

    /*
     * Check pointer to output variable:
     * allocate memory if necessary
     */
    ASSERT (gpstime != (LIGOTimeGPS *)NULL, status,
            UTOGPS_ENULLOUTPUT, UTOGPS_MSGENULLOUTPUT);


    /* shift the time to GPS time
     * GPS epoch is 1980-01-06, Unix epoch is 1970-01-01, so GPS time must
     * be smaller than UTC time */
    gpstime->gpsSeconds = unixtime->unixSeconds - UNIXGPS;

    RETURN (status);
} /* END LALUtoGPS() */


/*
 * Convert GPS seconds to UTC seconds
 *
 * Input:
 *
 * Output:
 */
void
LALGPStoU (LALStatus            *status,
        LIGOTimeUnix      *unixtime,
        const LIGOTimeGPS *gpstime)
{
    INITSTATUS (status, "LALGPStoU", UTOGPSC);

    /*
     * Check pointer to input variable
     */
    ASSERT (gpstime != (LIGOTimeGPS *)NULL, status,
            UTOGPS_ENULLINPUT, UTOGPS_MSGENULLINPUT);

    /*
     * Check pointer to output variable:
     * allocate memory if necessary
     */
    ASSERT (unixtime != (LIGOTimeUnix *)NULL, status,
            UTOGPS_ENULLOUTPUT, UTOGPS_MSGENULLOUTPUT);

    /* shift the time to UTC time
     * GPS epoch is 1980-01-06, Unix epoch is 1970-01-01, so GPS time must
     * be smaller than UTC time */
    unixtime->unixSeconds = gpstime->gpsSeconds + UNIXGPS;

    RETURN (status);
} /* END LALGPStoU() */

