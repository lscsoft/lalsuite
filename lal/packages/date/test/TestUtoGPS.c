#include "LALStdlib.h"
#include "Date.h"
#include "AVFactories.h"

/* $Id$ */

INT4 debuglevel = 2;

NRCSID (TESTUTOGPSC, "$Id$");

/* Difference between UTC and GPS (quoted from GRASP distribution)
 * 
 * difference between UTC and GPS time 315964811 = 
 *     3600 sec/hour x 24 hours/day x (365 days/year x 8 years + 
 *     366 days/year x 2 years + 5 days) + 11 leap seconds 
 */
const time_t UTCGPS = 315964811;

/* Taken from GRASP */
void printone(Status *status, const LIGOTimeUnix *time1)
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
    CHARCreateVector(status, &utc, (UINT4)64);
    CHARCreateVector(status, &gmt, (UINT4)64);
    CHARCreateVector(status, &gps, (UINT4)64);

    /* compute GPS time */
    UtoGPS(status, &gpstime, time1);

    /*
     * construct strings with appropriate time stamps
     */
    /* gmtime() */
    gmtime_r((time_t *)&(time1->unixSeconds), &utimestruct);
    strftime(gmt->data, gmt->length,
             "%Y-%m-%d %H:%M:%S UTC %a", &(utimestruct));

    /* Utime() */
    Utime(status, &laldate, time1);
    /* strcpy(utc, asctime(&(laldate.unixDate))); */
    DateString(status, utc, &laldate);

    /* GPS */
    /* tmp.unixSeconds += UTCGPS; */ /* shift back to Unix epoch */
    GPStoU(status, &tmp, &gpstime);
    Utime(status, &laldate, &tmp);
    DateString(status, gps, &laldate);
    /* strcpy(gps, asctime(&(laldate.unixDate))); */

    /* printf("C-time: %10d GPS Time %10d  UTC: %.24s GPS: %.24s gmtime: %s",
       (int)(time1->unixSeconds),(int)(gpstime.gpsSeconds),utc,gps,gmt); */
    printf("%10d | %10d | %.33s | %.33s | %s\n",
           (int)(time1->unixSeconds),(int)(gpstime.gpsSeconds),
           utc->data, gps->data, gmt->data);

    /*
     * House cleaning
     */
    CHARDestroyVector(status, &utc);
    CHARDestroyVector(status, &gmt);
    CHARDestroyVector(status, &gps);
    
    RETURN (status);
}

void
testunixtime(Status *status, const LIGOTimeUnix *time1)
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
    Status      status = {0};


    if (argc == 1)
    {
        /*
         * Print help message and exit
         */
        printf("Usage: TestUtoGPS debug_level -- debug_level = [0,1,2]\n");
        return 0;
    }

    if (argc == 2)
        debuglevel = atoi(argv[1]);

    printf("TEST Unix time and Unix to GPS conversion\n");
    printf("   and DateString\n");
    printf("=========================================\n");

    printf("\nLocal timezone: tzname[0] = %s\n                tzname[1] = %s\n", tzname[0], tzname[1]);

    printf("Current time:\t");

    time(&tmptime);
    unixtime.unixSeconds = tmptime;

    printf("unixtime = %13d;\t\t", unixtime.unixSeconds);

    UtoGPS(&status, &gpstime, &unixtime);

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
