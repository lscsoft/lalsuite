#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>

/* $Id$ */

extern char *tzname[2];
/* struct tm *gmtime_r( const time_t *, struct tm * ); */

INT4 lalDebugLevel = 2;

NRCSID (TESTUTOGPSC, "$Id$");

/* Difference between UTC and GPS (quoted from GRASP distribution)
 * 
 * difference between UTC and GPS time 315964811 = 
 *     3600 sec/hour x 24 hours/day x (365 days/year x 8 years + 
 *     366 days/year x 2 years + 5 days) + 11 leap seconds 
 */
const time_t UTCGPS = 315964811;

/* Taken from GRASP */
static void printone(LALStatus *status, LIGOTimeUnix *time1)
{
    LIGOTimeGPS  gpstime;
    LALDate      laldate;
    LIGOTimeUnix tmp;
    LALUnixDate  utimestruct;
    CHARVector *utc = NULL;
    CHARVector *gmt = NULL;
    CHARVector *gps = NULL;

    INITSTATUS (status, "printone", TESTUTOGPSC);

    /*
     * Allocate space for CHARVectors
     */
    LALCHARCreateVector(status, &utc, (UINT4)64);
    LALCHARCreateVector(status, &gmt, (UINT4)64);
    LALCHARCreateVector(status, &gps, (UINT4)64);

    /* compute GPS time */
    LALUtoGPS(status, &gpstime, time1);

    /*
     * construct strings with appropriate time stamps
     */
    /* gmtime() */
    gmtime_r((time_t *)&(time1->unixSeconds), &utimestruct);
    strftime(gmt->data, gmt->length,
             "%Y-%m-%d %H:%M:%S UTC %a", &(utimestruct));

    /* LALUtime() */
    LALUtime(status, &laldate, time1);
    /* strcpy(utc, asctime(&(laldate.unixDate))); */
    LALDateString(status, utc, &laldate);

    /* GPS */
    /* tmp.unixSeconds += UTCGPS; */ /* shift back to Unix epoch */
    LALGPStoU(status, &tmp, &gpstime);
    LALUtime(status, &laldate, &tmp);
    LALDateString(status, gps, &laldate);
    /* strcpy(gps, asctime(&(laldate.unixDate))); */

    /* printf("C-time: %10d GPS Time %10d  UTC: %.24s GPS: %.24s gmtime: %s",
       (int)(time1->unixSeconds),(int)(gpstime.gpsSeconds),utc,gps,gmt); */
    printf("%10d | %10d | %.33s | %.33s | %s\n",
           (int)(time1->unixSeconds),(int)(gpstime.gpsSeconds),
           utc->data, gps->data, gmt->data);

    /*
     * House cleaning
     */
    LALCHARDestroyVector(status, &utc);
    LALCHARDestroyVector(status, &gmt);
    LALCHARDestroyVector(status, &gps);
    
    RETURN (status);
}

static void
testunixtime(LALStatus *status, const LIGOTimeUnix *time1)
{
    LIGOTimeUnix tmp;

    INITSTATUS (status, "testunixtime", TESTUTOGPSC);
    
    tmp.unixSeconds     = time1->unixSeconds;
    tmp.unixNanoSeconds = time1->unixNanoSeconds;
    
    tmp.unixSeconds--;
    printone(status, &tmp);
    tmp.unixSeconds++;
    printone(status, &tmp);
    tmp.unixSeconds++;
    printone(status, &tmp);
    printf("           |            |                             |                             |         \n");

    RETURN (status);
}
        
    

int
main(int argc, char *argv[])
{
    LIGOTimeGPS  gpstime;
    LIGOTimeUnix unixtime;
    time_t      tmptime;
    static LALStatus status;


    if (argc == 1)
    {
        /*
         * Print help message and exit
         */
        printf("Usage: TestUtoGPS debug_level -- debug_level = [0,1,2]\n");
        return 0;
    }

    if (argc == 2)
        lalDebugLevel = atoi(argv[1]);

    printf("TEST Unix time and Unix to GPS conversion\n");
    printf("   and LALDateString\n");
    printf("=========================================\n");

    printf("\nLocal timezone: tzname[0] = %s\n                tzname[1] = %s\n", tzname[0], tzname[1]);

    printf("Current time:\t");

    time(&tmptime);
    unixtime.unixSeconds = tmptime;

    printf("unixtime = %13d;\t\t", unixtime.unixSeconds);

    LALUtoGPS(&status, &gpstime, &unixtime);

    printf("gpstime = %13d\n", gpstime.gpsSeconds);
    printf("\n");

    printf("Various times, some which include leap second - look for '60' in the seconds field of the time:\n\n");

    printf(" Unix time | GPS Time   |          UTC                |          GPS                | gmtime_r()\n");
    printf("================================================================================================================\n");

    unixtime.unixSeconds = 1;
    testunixtime(&status, &unixtime);
    unixtime.unixSeconds = 33350400;
    testunixtime(&status, &unixtime);
    unixtime.unixSeconds = 63072000;
    testunixtime(&status, &unixtime);
    unixtime.unixSeconds = 78796801;
    testunixtime(&status, &unixtime);
    unixtime.unixSeconds = 94694402;
    testunixtime(&status, &unixtime);
    unixtime.unixSeconds = 126230403;
    testunixtime(&status, &unixtime);
    unixtime.unixSeconds = 315964811;
    testunixtime(&status, &unixtime);
    unixtime.unixSeconds = 784880277;
    testunixtime(&status, &unixtime);
    unixtime.unixSeconds = 911110677;
    testunixtime(&status, &unixtime);

    return 0;
}
