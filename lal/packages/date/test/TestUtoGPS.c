#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>


/* $Id$ */

extern char *tzname[2];
/* struct tm *gmtime_r( const time_t *, struct tm * ); */

INT4 lalDebugLevel = LALALLDBG;

NRCSID (TESTUTOGPSC, "$Id$");

/* Difference between UTC and GPS (quoted from GRASP distribution)
 * 
 * difference between UTC and GPS time 315964811 = 
 *     3600 sec/hour x 24 hours/day x (365 days/year x 8 years + 
 *     366 days/year x 2 years + 5 days) + 11 leap seconds 
 */
const time_t UTCGPS = 315964811;

/* Taken from GRASP */
static void printone(LALStatus *status, LIGOTimeUnix *time1, int has_leap_p,
                     int *retval)
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
    if (lalDebugLevel > 2)
      REPORTSTATUS(status);
    LALCHARCreateVector(status, &gmt, (UINT4)64);
    if (lalDebugLevel > 2)
      REPORTSTATUS(status);
    LALCHARCreateVector(status, &gps, (UINT4)64);
    if (lalDebugLevel > 2)
      REPORTSTATUS(status);


    /* compute GPS time */
    LALUtoGPS(status, &gpstime, time1);
    if (lalDebugLevel > 2)
      REPORTSTATUS(status);

    /*
     * construct strings with appropriate time stamps
     */
    /* gmtime() */
    gmtime_r((time_t *)&(time1->unixSeconds), &utimestruct);
    if (strftime(gmt->data, gmt->length,
                 "%Y-%m-%d %H:%M:%S UTC %a", &(utimestruct)) == 0)
      {
        exit(1);
      }

    /* LALUtime() */
    LALUtime(status, &laldate, time1);
    if (lalDebugLevel > 2)
      REPORTSTATUS(status);

    if (has_leap_p > 0)
      if (laldate.unixDate.tm_sec == 60) 
        *retval = 0;
      else
        *retval = 1;
    else
      *retval = 0;

    LALDateString(status, utc, &laldate);
    if (lalDebugLevel > 2)
      REPORTSTATUS(status);

    /* GPS */
    /* tmp.unixSeconds += UTCGPS; */ /* shift back to Unix epoch */
    LALGPStoU(status, &tmp, &gpstime);
    if (lalDebugLevel > 2)
      REPORTSTATUS(status);
    LALUtime(status, &laldate, &tmp);
    if (lalDebugLevel > 2)
      REPORTSTATUS(status);
    LALDateString(status, gps, &laldate);
    if (lalDebugLevel > 2)
      REPORTSTATUS(status);


    /* printf("C-time: %10d GPS Time %10d  UTC: %.24s GPS: %.24s gmtime: %s",
       (int)(time1->unixSeconds),(int)(gpstime.gpsSeconds),utc,gps,gmt); */
    printf("%10d | %10d | %.33s | %.33s | %s\n",
           (int)(time1->unixSeconds),(int)(gpstime.gpsSeconds),
           utc->data, gps->data, gmt->data);

    /*
     * House cleaning
     */
    LALCHARDestroyVector(status, &utc);
    if (lalDebugLevel > 2)
      REPORTSTATUS(status);

    LALCHARDestroyVector(status, &gmt);
    if (lalDebugLevel > 2)
      REPORTSTATUS(status);

    LALCHARDestroyVector(status, &gps);
    if (lalDebugLevel > 2)
      REPORTSTATUS(status);
    
    RETURN (status);
}

static void
testunixtime(LALStatus *status, const LIGOTimeUnix *time1, int has_leap_p,
             int *retval)
{
    LIGOTimeUnix tmp;


    tmp.unixSeconds     = time1->unixSeconds;
    tmp.unixNanoSeconds = time1->unixNanoSeconds;

    tmp.unixSeconds--;
    printone(status, &tmp, has_leap_p, retval);
    tmp.unixSeconds++;
    printone(status, &tmp, has_leap_p, retval);
    tmp.unixSeconds++;
    printone(status, &tmp, has_leap_p, retval);
    printf("           |            |                             |                             |         \n");

    
    RETURN (status);
}
        
    

int
main(int argc, char *argv[])
{
    LIGOTimeGPS  gpstime;
    LIGOTimeUnix unixtime;
    LALDate      date;
    time_t       now;
    CHARVector  *timestamp = NULL;
    time_t       tmptime;
    int          retval = 0;
    static LALStatus status;


    if (argc == 1)
    {
        /*
         * Print help message and exit
         */
        printf("Usage: TestUtoGPS debug_level -- debug_level = [0,1,2]\n");
        lalDebugLevel = 0;
    }

    if (argc == 2)
        lalDebugLevel = atoi(argv[1]);

    LALCHARCreateVector(&status, &timestamp, (UINT4)128);
    if (lalDebugLevel > 2) 
      REPORTSTATUS(&status);

    printf("TEST Unix time and Unix to GPS conversion\n");
    printf("   and LALDateString\n");
    printf("=========================================\n");

    printf("\nLocal timezone: tzname[0] = %s\n                tzname[1] = %s\n", tzname[0], tzname[1]);

    printf("Current time:\t");

    tmptime = time(NULL);
    unixtime.unixSeconds = tmptime;

    printf("unixtime = %13d;\t\t", unixtime.unixSeconds);

    LALUtoGPS(&status, &gpstime, &unixtime);
    if (lalDebugLevel > 2)
      REPORTSTATUS(&status);

    printf("gpstime = %13d\n", gpstime.gpsSeconds);
    printf("\n");

    tmptime = time(NULL);
    unixtime.unixSeconds = tmptime;
    unixtime.unixNanoSeconds = 0;
    LALUtime(&status, &date, &unixtime);
    LALDateString(&status, timestamp, &date);

    printf("Various times, some which include leap second - look for '60' in the seconds field of the time:\n\n");

    printf(" Unix time | GPS Time   |        UTC (LALUtime)       |          GPS                | gmtime_r()\n");
    printf("================================================================================================================\n");

    unixtime.unixSeconds = 1;
    testunixtime(&status, &unixtime, 0, &retval);
    unixtime.unixSeconds = 33350400;
    testunixtime(&status, &unixtime, 1, &retval);
    unixtime.unixSeconds = 63072001;
    testunixtime(&status, &unixtime, 1, &retval);
    unixtime.unixSeconds = 78796802;
    testunixtime(&status, &unixtime, 1, &retval);
    unixtime.unixSeconds = 94694403;
    testunixtime(&status, &unixtime, 1, &retval);
    unixtime.unixSeconds = 126230404;
    testunixtime(&status, &unixtime, 1, &retval);
    unixtime.unixSeconds = 315964811;
    testunixtime(&status, &unixtime, 0, &retval);
    unixtime.unixSeconds = 784880277;
    testunixtime(&status, &unixtime, 0, &retval);
    unixtime.unixSeconds = 911110677;
    testunixtime(&status, &unixtime, 0, &retval);
    unixtime.unixSeconds = 773539211 + 21;
    testunixtime(&status, &unixtime, 0, &retval);
    unixtime.unixSeconds = time(NULL);
    testunixtime(&status, &unixtime, 0, &retval);

    return retval;
}
