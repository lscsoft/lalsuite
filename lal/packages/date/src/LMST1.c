/*----------------------------------------------------------------------- 
 * 
 * File Name: LMST1.c 
 * 
 * Author: David Chin <dwchin@umich.edu>
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 *
 * SYNOPSIS
 *
 * GMST1(): Returns GMST1 for given date-time
 * LMST1(): Returns LMST1 given date-time, and longitude
 * 
 * DESCRIPTION
 *
 * GMST1():
 *      Inputs: LALDate *date      -- date-time to compute GMST1
 *              INT4    outunits   -- units requested:
 *                                    MST_SEC - seconds
 *                                    MST_HRS - hours
 *                                    MST_DEG - degrees
 *                                    MST_RAD - radians
 *
 *      Outputs: REAL8  *gmst      -- GMST1 in units requested
 *
 * LMST1():
 *      Inputs: LALDate *date      -- date-time to compute GMST1
 *              REAL8   longitude  -- longitude (in decimal degrees East)
 *              INT4    outunits   -- units requested:
 *                                    MST_SEC - seconds
 *                                    MST_HRS - hours
 *                                    MST_DEG - degrees
 *                                    MST_RAD - radians
 *
 *      Outputs: REAL8  *lmst      -- GMST1 in units requested
 *
 * 
 * DIAGNOSTICS 
 * (Abnormal termination conditions, error and warning codes summarized 
 * here. More complete descriptions are found in documentation.)
 *
 * CALLS
 *    GMST1():  JulianDate()
 *    LMST1():  GMST1()
 * 
 * NOTES
 * This is from the Explanatory Supplement to the Astronomical Almanac,
 * 1992, Ch. 2, Sec. 24, and also Section B of the Almanac. The formula
 * computes GMST1 for 0h UT1.  To  compute GMST1 for other times, a simple
 * linear interpolation is done.  
 * 
 *-----------------------------------------------------------------------
 */

#include "LALRCSID.h"

NRCSID (LMST1C, "$Id$");

#include "Date.h"
#include "date_value.h"


/*
 * Compute GMST1 in requested units given date UTC.
 */
void
GMST1 (Status        *status,
       REAL8         *gmst,
       const LALDate *date,
       INT4           outunits)
{
    REAL8 du;
    REAL8 tu;
    REAL8 tu2;
    LALDate tmpdate;
    REAL8   interval;
    REAL8 ut2mst; /* conversion factor: UT1 -> Mean Sidereal Time */
    const REAL8      a =   24110.54841;
    const REAL8      b = 8640185.812866;
    const REAL8      c =       0.093104;
    const REAL8      d =      -6.2e-6;
    
    INITSTATUS (status, "GMST1", LMST1C);

    /*
     * Check pointer to input variables
     */
    ASSERT (date != (LALDate *)NULL, status,
            GMST1_ENULLINPUT, GMST1_MSGENULLINPUT);

    /*
     * Check pointer to output variable
     */
    ASSERT (gmst != (REAL8 *)NULL, status,
            GMST1_ENULLOUTPUT, GMST1_MSGENULLOUTPUT);

    /*
     * Make a date structure for midnight of the input date
     */
    tmpdate.unixDate.tm_sec   = 0;
    tmpdate.unixDate.tm_min   = 0;
    tmpdate.unixDate.tm_hour  = 0;
    tmpdate.unixDate.tm_mday  = date->unixDate.tm_mday;
    tmpdate.unixDate.tm_mon   = date->unixDate.tm_mon;
    tmpdate.unixDate.tm_year  = date->unixDate.tm_year;
    
    /*
     * Compute GMST1 for 0h UT1 on given date in seconds since J2000.0 
     */
    JulianDate(status, &du, &tmpdate);
    du -= J2000_0;
    tu  = du / JDAYS_PER_CENT;
    tu2 = tu * tu;
    
    *gmst = a + b * tu + c * tu2 + d * tu * tu2;

    while (*gmst < 0)
        *gmst += (REAL8)SECS_PER_DAY;

    /*
     * Add the equivalent mean sidereal time interval
     * from 0h to given time UT1
     */
    interval = (REAL8)(date->unixDate.tm_sec) +
        60. * (REAL8)(date->unixDate.tm_min) +
        3600. * (REAL8)(date->unixDate.tm_hour);

    /* Eqn 2.241-3 of Expl. Suppl. to Astr. Almanac (this is probably
     * overkill at REAL8 precision) */
    ut2mst = 1.002737909350795 + 5.9006e-11*tu - 5.9e-15*tu2;
    
    *gmst += ut2mst * interval;

    switch (outunits)
    {
    case MST_SEC:
        break;
    case MST_HRS:
        *gmst /= (REAL8)SECS_PER_HOUR;
        break;
    case MST_DEG:
        *gmst /= (REAL8)SECS_PER_HOUR / (REAL8)DEGS_PER_HOUR;
        break;
    case MST_RAD:
        *gmst /= (REAL8)SECS_PER_HOUR / (REAL8)DEGS_PER_HOUR * 180. /
            LAL_PI;
        break;
    default:
        break;
    }

    RETURN (status);
} /* END GMST1() */



/*
 * Compute LMST1 in requested units given date-time, and longitude in degrees
 */
void
LMST1 (Status        *status,
       REAL8         *lmst,
       const LALDate *date,
       REAL8          longitude,
       INT4           outunits)
{
    REAL8 gmst;
	REAL8 day = 0;

    INITSTATUS (status, "LMST1", LMST1C);

    /*
     * Check pointer to input variables
     */
    ASSERT (date != (LALDate *)NULL, status,
            LMST1_ENULLINPUT, LMST1_MSGENULLINPUT);

    /*
     * Check pointer to output variable
     */
    ASSERT (lmst != (REAL8 *)NULL, status,
            LMST1_ENULLOUTPUT, LMST1_MSGENULLOUTPUT);

    /*
     * Compute LMST1 in seconds since J2000.0
     */

    /* get GMST1 in seconds */
    GMST1(status, &gmst, date, outunits);

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
        longitude /= 180. / LAL_PI;
		day        = LAL_TWOPI;
        break;
    default:
        break;
    }

    *lmst = gmst + longitude;

	while (*lmst < 0)
        *lmst += day;

    RETURN (status);
} /* END LMST1() */
