#include "LALStdlib.h"
#include "Date.h"

INT4 lalDebugLevel = 2;

NRCSID (TESTJULIANDAYC, "$Id$");

int
main(int argc, char *argv[])
{
    time_t        now;
    LALUnixDate  *ltime;
    LALStatus        status = {0};
    LALDate       date;
    INT4          julian_day;
    REAL8         mod_julian_day;
    REAL8         julian_date;
    /* const REAL8   longitude = -118.133; */ /* CalTech */


    if (argc == 1)
    {
        /*
         * Print help message and exit
         */
        printf("Usage: TestJulianDay debug_level -- debug_level = [0,1,2]\n");
        return 0;
    }

    if (argc == 2)
        lalDebugLevel = atoi(argv[1]);

    /*
     * Get current local time
     */
    time(&now);
    ltime = localtime(&now);
    
    date.unixDate.tm_sec  = ltime->tm_sec;
    date.unixDate.tm_min  = ltime->tm_min;
    date.unixDate.tm_hour = ltime->tm_hour;
    date.unixDate.tm_mday = ltime->tm_mday;
    date.unixDate.tm_mon  = ltime->tm_mon;
    date.unixDate.tm_year = ltime->tm_year;
    date.unixDate.tm_wday = ltime->tm_wday;
    date.unixDate.tm_yday = ltime->tm_yday;
    date.unixDate.tm_isdst = ltime->tm_isdst;

    /*
     * Print out the structure members to check
     */
    printf("TEST OF JULIAN DAY and DATE ROUTINE\n");
    printf("===================================\n");
    printf("\nFirst, test the LALDate structure\n");
    printf("You should see the current date and time:\n");
    printf("\tsec  = %d\t(seconds after the minute - [0,61])\n",
           date.unixDate.tm_sec);
    printf("\tmin  = %d\t(minutes after the hour - [0,59]\n",
           date.unixDate.tm_min);
    printf("\thour = %d\t(hours - [0,23])\n", date.unixDate.tm_hour);
    printf("\tmday = %d\t(day of month - [1,31])\n", date.unixDate.tm_mday);
    printf("\tmon  = %d\t(month of year - [0,11])\n", date.unixDate.tm_mon);
    printf("\tyear = %d\t(yrs since 1900)\n", date.unixDate.tm_year);
    printf("\twday = %d\t(days since Sunday - [0,6])\n", date.unixDate.tm_wday);
    printf("\tyday = %d\t(days since January 1 - [0,365])\n",
           date.unixDate.tm_yday);
    printf("\tisdst = %d\t(daylight savings time flag)\n",
           date.unixDate.tm_isdst);

    /*
     * Convert above date to Julian Day and Date
     */
    printf("\n");
    printf("Find Julian Day/Date for current date and time:\n");
    
    LALJulianDay(&status, &julian_day, &date);
    printf("\n");
    printf("\tJulian Day Number          = %9d\n", julian_day);

    LALModJulianDay(&status, &mod_julian_day, &date);
    printf("\tModified Julian Day Number = %13.3f\n", mod_julian_day);
    
    LALJulianDate(&status, &julian_date, &date);
    printf("\tJulian Date                = %13.3f\n", julian_date);

    LALModJulianDate(&status, &julian_date, &date);
    printf("\tModified Julian Date       = %13.3f\n", julian_date);

    
    /*
     * Check Julian Day/Date for special dates and times
     */
    printf("\n");
    printf("Find Julian Day/Date for 2000-Jan-01 12h 00m 00s:\n");
    printf("(expect J2000.0 = 2451545.0)\n");

    date.unixDate.tm_sec  = 0;
    date.unixDate.tm_min  = 0;
    date.unixDate.tm_hour = 12;
    date.unixDate.tm_mday = 1;
    date.unixDate.tm_mon  = 0;
    date.unixDate.tm_year = 100;
    date.unixDate.tm_wday = 6;
    date.unixDate.tm_yday = 0;
    date.unixDate.tm_isdst = 0;
    
    LALJulianDay(&status, &julian_day, &date);
    printf("\n");
    printf("\tJulian Day Number          = %8d\n", julian_day);

    LALModJulianDay(&status, &mod_julian_day, &date);
    printf("\tModified Julian Day Number = %10.1f\n", mod_julian_day);
    
    LALJulianDate(&status, &julian_date, &date);
    printf("\tJulian Date                = %10.1f\n", julian_date);

    LALModJulianDate(&status, &julian_date, &date);
    printf("\tModified Julian Date       = %10.1f\n", julian_date);

    /* */
    printf("\n");
    printf("Find Julian Day/Date for 2000-Jan-01 11h 00m 00s:\n");

    date.unixDate.tm_sec  = 0;
    date.unixDate.tm_min  = 0;
    date.unixDate.tm_hour = 11;
    date.unixDate.tm_mday = 1;
    date.unixDate.tm_mon  = 0;
    date.unixDate.tm_year = 100;
    date.unixDate.tm_wday = 6;
    date.unixDate.tm_yday = 0;
    date.unixDate.tm_isdst = 0;
    
    LALJulianDay(&status, &julian_day, &date);
    printf("\n");
    printf("\tJulian Day Number          = %8d\n", julian_day);

    LALModJulianDay(&status, &mod_julian_day, &date);
    printf("\tModified Julian Day Number = %10.1f\n", mod_julian_day);
    
    LALJulianDate(&status, &julian_date, &date);
    printf("\tJulian Date                = %10.1f\n", julian_date);

    LALModJulianDate(&status, &julian_date, &date);
    printf("\tModified Julian Date       = %10.1f\n", julian_date);

    /* */
    printf("\n");
    printf("Find Julian Day/Date for 1800-Jan-01 11h 00m 00s:\n");
    printf("(this should produce errors at lalDebugLevel > 0)\n");

    date.unixDate.tm_sec  = 0;
    date.unixDate.tm_min  = 0;
    date.unixDate.tm_hour = 11;
    date.unixDate.tm_mday = 1;
    date.unixDate.tm_mon  = 0;
    date.unixDate.tm_year = -100;
    date.unixDate.tm_wday = 6;
    date.unixDate.tm_yday = 0;
    date.unixDate.tm_isdst = 0;
    
    LALJulianDay(&status, &julian_day, &date);
    printf("\n");
    printf("\tJulian Day Number          = %8d\n", julian_day);

    LALModJulianDay(&status, &mod_julian_day, &date);
    printf("\tModified Julian Day Number = %10.1f\n", mod_julian_day);
    
    LALJulianDate(&status, &julian_date, &date);
    printf("\tJulian Date                = %10.1f\n", julian_date);

    LALModJulianDate(&status, &julian_date, &date);
    printf("\tModified Julian Date       = %10.1f\n", julian_date);

    /* */
    printf("\n");
    printf("Find Julian Day/Date for 1994-Nov-16 00h 00m 00s:\n");
    printf("(Almanac: Julian Date 2449672.5;  ");
    printf("GMST 03h 39m 20.6222)\n");

    date.unixDate.tm_sec  = 0;
    date.unixDate.tm_min  = 0;
    date.unixDate.tm_hour = 0;
    date.unixDate.tm_mday = 16;
    date.unixDate.tm_mon  = 10;
    date.unixDate.tm_year = 94;
    date.unixDate.tm_wday = 3;
    date.unixDate.tm_yday = 0;
    date.unixDate.tm_isdst = 0;
    
    LALJulianDay(&status, &julian_day, &date);
    printf("\n");
    printf("\tJulian Day Number          = %8d\n", julian_day);

    LALModJulianDay(&status, &mod_julian_day, &date);
    printf("\tModified Julian Day Number = %10.1f\n", mod_julian_day);
    
    LALJulianDate(&status, &julian_date, &date);
    printf("\tJulian Date                = %10.1f\n", julian_date);

    LALModJulianDate(&status, &julian_date, &date);
    printf("\tModified Julian Date       = %10.1f\n", julian_date);

    return 0;
}

