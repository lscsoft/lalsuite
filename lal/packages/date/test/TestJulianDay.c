#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/Date.h>

INT4 lalDebugLevel = 0;

NRCSID (LALTESTJULIANDAYC, "$Id$");

#define SUCCESS              0
#define FAIL_JULIAN_DAY      1
#define FAIL_MOD_JULIAN_DAY  2
#define FAIL_JULIAN_DATE     3
#define FAIL_MOD_JULIAN_DATE 4

/* good to 1 micro-day  */
const REAL8 julian_precision = 1.e-6;

const REAL8 coarse_precision = 0.001;

/*  Output of NASA/AMES code:
 Date/time: 1 1 0 12 0 0
 PDS format: "2000-01-01T12:00:00.000Z"
 SQL format: "Jan 1, 2000 12:00:00:000"

  UTC day  =        0
  UTC secs =    43200.000000
  TAI secs =          32.000000
   ET secs =          64.183927

  JD (UTC) =  2451545.000000
  JD (TAI) =  2451545.000370
  JD (ET)  =  2451545.001115
 MJD (UTC) =    51544.500000
 MJD (TAI) =    51544.500370
 MJD (ET)  =    51544.501115

  Date/time: 94 11 16 0 0 0
 PDS format: "1994-11-16T00:00:00.000Z"
 SQL format: "Nov 16, 1994 00:00:00:000"

  UTC day  =    -1872
  UTC secs =        0.000000
  TAI secs =  -161783971.000000
   ET secs =  -161783938.817245

  JD (UTC) =  2449672.500000
  JD (TAI) =  2449672.500336
  JD (ET)  =  2449672.501081
 MJD (UTC) =    49672.000000
 MJD (TAI) =    49672.000336
 MJD (ET)  =    49672.001081

 Date/time: 01 05 15 02 37 54
 PDS format: "2001-05-15T02:37:54.000Z"
 SQL format: "May 15, 2001 02:37:54:000"

  UTC day  =      500
  UTC secs =     9474.000000
  TAI secs =    43166306.000000
   ET secs =    43166338.185257

  JD (UTC) =  2452044.609653
  JD (TAI) =  2452044.610023
  JD (ET)  =  2452044.610768
 MJD (UTC) =    52044.109653
 MJD (TAI) =    52044.110023
 MJD (ET)  =    52044.110768


 http://tycho.usno.navy.mil/cgi-bin/date
15-May-01
MJD 52044.109653
UTC 02:37:54

*/


int
main(int argc, char *argv[])
{
    time_t        now;
    LALUnixDate  *ltime;
    static LALStatus status;
    LALDate       date;
    INT4          ref_julian_day;
    INT4          julian_day;
    INT4          old_debuglvl;
    REAL8         ref_mod_julian_day;
    REAL8         mod_julian_day;
    REAL8         ref_julian_date;
    REAL8         julian_date;
    REAL8         ref_mod_julian_date;
    REAL8         mod_julian_date;


    /*
     * Distinctly not robust
     */
    if (argc == 1)
      fprintf(stderr,
              "Usage: TestJulianDay [debug_level] -- debug_level = [0,1,2]\n");

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
    date.residualNanoSeconds = 0;

    if (lalDebugLevel > 1) {
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
      printf("\twday = %d\t(days since Sunday - [0,6])\n",
             date.unixDate.tm_wday);
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
      if (status.statusCode && lalDebugLevel > 0)
        {
          fprintf(stderr, "TestJulianDay: LALJulianDay() failed, line %i, %s\n",
                  __LINE__, LALTESTJULIANDAYC);
          REPORTSTATUS(&status);
          return status.statusCode;
        }
      
      printf("\n");
      printf("\tJulian Day Number          = %8d\n", julian_day);

      LALJulianDate(&status, &julian_date, &date);
      if (status.statusCode && lalDebugLevel > 0)
        {
          fprintf(stderr,
                  "TestJulianDay: LALJulianDate() failed, line %i, %s\n",
                  __LINE__, LALTESTJULIANDAYC);
          REPORTSTATUS(&status);
          return status.statusCode;
        }
      printf("\tJulian Date                = %19.10f\n", julian_date);

      LALModJulianDate(&status, &julian_date, &date);
      if (status.statusCode && lalDebugLevel > 0)
        {
          fprintf(stderr,
                  "TestJulianDay: LALModJulianDate() failed, line %i, %s\n",
                  __LINE__, LALTESTJULIANDAYC);
          REPORTSTATUS(&status);
          return status.statusCode;
        }
      printf("\tModified Julian Date       = %19.10f\n", julian_date);

    }

    
    /*
     * Check Julian Day/Date for special dates and times
     */

    if (lalDebugLevel > 0) {
      printf("\n");
      printf("Find Julian Day/Date for 2000-Jan-01 12h 00m 00s:\n");
      printf("(expect J2000.0 = 2451545.0)\n");
    }

    date.unixDate.tm_sec  = 0;
    date.unixDate.tm_min  = 0;
    date.unixDate.tm_hour = 12;
    date.unixDate.tm_mday = 1;
    date.unixDate.tm_mon  = 0;
    date.unixDate.tm_year = 100;
    date.unixDate.tm_wday = 6;
    date.unixDate.tm_yday = 0;
    date.unixDate.tm_isdst = 0;

    ref_julian_day = 2451545;
    LALJulianDay(&status, &julian_day, &date);
    if (status.statusCode && lalDebugLevel > 0)
      {
        fprintf(stderr,
                "TestJulianDay: LALJulianDay() failed, line %i, %s\n",
                __LINE__, LALTESTJULIANDAYC);
        REPORTSTATUS(&status);
        return status.statusCode;
      }

      
    if (lalDebugLevel > 0) {
      printf("\n");
      printf("\tJulian Day Number          = %8d\n", julian_day);
    }

    if (julian_day != ref_julian_day) {
      fprintf(stderr, "FAIL: Julian Day: value not as expected; line %i, %s\n",
              __LINE__, LALTESTJULIANDAYC);
      return FAIL_JULIAN_DAY;
    }

    ref_julian_date = 2451545.0;
    LALJulianDate(&status, &julian_date, &date);
    if (status.statusCode && lalDebugLevel > 0)
      {
        fprintf(stderr,
                "TestJulianDay: LALJulianDay() failed, line %i, %s\n",
                __LINE__, LALTESTJULIANDAYC);
        REPORTSTATUS(&status);
        return status.statusCode;
      }

    
    if (lalDebugLevel > 0) 
      printf("\tJulian Date                = %19.10f\n", julian_date);

    if (fabs(julian_date - ref_julian_date) > julian_precision) {

      if (fabs(julian_date - ref_julian_date) > coarse_precision) {  
        fprintf(stderr, "FAIL: Julian Date, line %i, %s\n",
                __LINE__, LALTESTJULIANDAYC);
        return FAIL_JULIAN_DATE;
      } else {
        fprintf(stderr, "WARNING: Inaccuracy in LALJulianDate() > %e\n",
                julian_precision);
      }
    }

    ref_mod_julian_date = 51544.5;
    LALModJulianDate(&status, &mod_julian_date, &date);
    if (status.statusCode && lalDebugLevel > 0)
      {
        fprintf(stderr,
                "TestJulianDay: LALModJulianDate() failed, line %i, %s\n",
                __LINE__, LALTESTJULIANDAYC);
        REPORTSTATUS(&status);
        return status.statusCode;
      }



    if (lalDebugLevel > 0)
      printf("\tModified Julian Date       = %19.10f\n", mod_julian_date);

    if (fabs(mod_julian_date - ref_mod_julian_date) > julian_precision) {
      
      if (fabs(mod_julian_date - ref_mod_julian_date) > coarse_precision) {
        fprintf(stderr, "FAIL: Modified Julian Date, line %i, %s\n",
                __LINE__, LALTESTJULIANDAYC);
        return FAIL_MOD_JULIAN_DATE;
      } else {
        fprintf(stderr, "WARNING: Inaccuracy in LALModJulianDate() > %e\n",
                julian_precision);
      }
    }


        
    /* */
    if (lalDebugLevel > 0) {
      printf("\n");
      printf("Find Julian Day/Date for 2000-Jan-01 11h 00m 00s:\n");
    }

    date.unixDate.tm_sec  = 0;
    date.unixDate.tm_min  = 0;
    date.unixDate.tm_hour = 11;
    date.unixDate.tm_mday = 1;
    date.unixDate.tm_mon  = 0;
    date.unixDate.tm_year = 100;
    date.unixDate.tm_wday = 6;
    date.unixDate.tm_yday = 0;
    date.unixDate.tm_isdst = 0;

    ref_julian_day = 2451545;
    LALJulianDay(&status, &julian_day, &date);
    if (status.statusCode && lalDebugLevel > 0)
      {
        fprintf(stderr,
                "TestJulianDay: LALJulianDay() failed, line %i, %s\n",
                __LINE__, LALTESTJULIANDAYC);
        REPORTSTATUS(&status);
        return status.statusCode;
      }



    if (lalDebugLevel > 0) {
      printf("\n");
      printf("\tJulian Day Number          = %8d\n", julian_day);
    }

    if (julian_day != ref_julian_day) {
      fprintf(stderr, "FAIL: Julian Day, line %i, %s\n", __LINE__,
              LALTESTJULIANDAYC);
      return FAIL_JULIAN_DAY;
    }


    ref_julian_date = 2451544.9583333344;
    LALJulianDate(&status, &julian_date, &date);
    if (status.statusCode && lalDebugLevel > 0)
      {
        fprintf(stderr,
                "TestJulianDay: LALJulianDate() failed, line %i, %s\n",
                __LINE__, LALTESTJULIANDAYC);
        REPORTSTATUS(&status);
        return status.statusCode;
      }


    if (lalDebugLevel > 0)
      printf("\tJulian Date                = %19.10f\n", julian_date);

    if (fabs(julian_date - ref_julian_date) > julian_precision) {

      if (fabs(julian_date - ref_julian_date) > coarse_precision) {
        fprintf(stderr, "FAIL: Julian Date\n");
        return FAIL_JULIAN_DATE;
      } else {
        fprintf(stderr, "WARNING: Inaccuracy in LALJulianDate() > %e\n",
                julian_precision);
      }
    }

    ref_mod_julian_date = 51544.4583333344;
    LALModJulianDate(&status, &mod_julian_date, &date);
    if (status.statusCode && lalDebugLevel > 0)
      {
        fprintf(stderr,
                "TestJulianDay: LALModJulianDate() failed, line %i, %s\n",
                __LINE__, LALTESTJULIANDAYC);
        REPORTSTATUS(&status);
        return status.statusCode;
      }


    if (lalDebugLevel > 0)
      printf("\tModified Julian Date       = %19.10f\n", mod_julian_date);

    if (fabs(mod_julian_date - ref_mod_julian_date) > julian_precision) {

      if (fabs(mod_julian_date - ref_mod_julian_date) > coarse_precision) {  
        fprintf(stderr, "FAIL: Modified Julian Date\n");
        return FAIL_MOD_JULIAN_DATE;
      } else {
        fprintf(stderr, "WARNING: Inaccuracy in LALJulianDate() > %e\n",
                julian_precision);
      }
    }
      

    /* */
    if (lalDebugLevel > 1) {
      old_debuglvl = lalDebugLevel;

      lalDebugLevel = 1;  /* LAL seems to not consider any debug level > 1 */

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
      /* here, we EXPECT an error */
      if (!status.statusCode && lalDebugLevel > 0)
      {
        fprintf(stderr,
                "TestJulianDay: LALJulianDay() failed, line %i, %s\n",
                __LINE__, LALTESTJULIANDAYC);
        REPORTSTATUS(&status);
        return status.statusCode;
      }

      printf("\n");
      printf("\tJulian Day Number          = %8d\n", julian_day);

      LALJulianDate(&status, &julian_date, &date);
      if (status.statusCode && lalDebugLevel > 0)
      {
        fprintf(stderr,
                "TestJulianDay: LALJulianDate() failed, line %i, %s\n",
                __LINE__, LALTESTJULIANDAYC);
        REPORTSTATUS(&status);
        return status.statusCode;
      }

      printf("\tJulian Date                = %19.10f\n", julian_date);

      LALModJulianDate(&status, &julian_date, &date);
      if (status.statusCode && lalDebugLevel > 0)
      {
        fprintf(stderr,
                "TestJulianDay: LALModJulianDate() failed, line %i, %s\n",
                __LINE__, LALTESTJULIANDAYC);
        REPORTSTATUS(&status);
        return status.statusCode;
      }

      printf("\tModified Julian Date       = %19.10f\n", mod_julian_date);

      lalDebugLevel = old_debuglvl;
    }


    
    /* */

    if (lalDebugLevel > 0) {
      printf("\n");
      printf("Find Julian Day/Date for 1994-Nov-16 00h 00m 00s:\n");
      printf("(Almanac: Julian Date 2449672.5;  ");
      printf("GMST 03h 39m 20.6222)\n\n");
    }

    date.unixDate.tm_sec  = 0;
    date.unixDate.tm_min  = 0;
    date.unixDate.tm_hour = 0;
    date.unixDate.tm_mday = 16;
    date.unixDate.tm_mon  = 10;
    date.unixDate.tm_year = 94;
    date.unixDate.tm_wday = 3;
    date.unixDate.tm_yday = 0;
    date.unixDate.tm_isdst = 0;

    ref_julian_day = 2449673;
    LALJulianDay(&status, &julian_day, &date);
    if (status.statusCode && lalDebugLevel > 0)
      {
        fprintf(stderr,
                "TestJulianDay: LALJulianDay() failed, line %i, %s\n",
                __LINE__, LALTESTJULIANDAYC);
        REPORTSTATUS(&status);
        return status.statusCode;
      }


    if (lalDebugLevel > 0)
      printf("\tJulian Day Number          = %8d\n", julian_day);

    if (julian_day != ref_julian_day) {
      fprintf(stderr, "FAIL: Julian Day\n");
      return FAIL_JULIAN_DAY;
    }

    ref_julian_date = 2449672.5;
    LALJulianDate(&status, &julian_date, &date);
    if (status.statusCode && lalDebugLevel > 0)
      {
        fprintf(stderr,
                "TestJulianDay: LALJulianDate() failed, line %i, %s\n",
                __LINE__, LALTESTJULIANDAYC);
        REPORTSTATUS(&status);
        return status.statusCode;
      }


    if (lalDebugLevel > 0)
      printf("\tJulian Date                = %19.10f\n", julian_date);

    if (fabs(julian_date - ref_julian_date) > julian_precision) {

      if (fabs(julian_date - ref_julian_date) > coarse_precision) {
        fprintf(stderr, "FAIL: Julian Date\n");
        return FAIL_JULIAN_DATE;
      } else {
        fprintf(stderr, "WARNING: Inaccuracy in LALJulianDate() > %e\n",
                julian_precision);
      }
    }


    ref_mod_julian_date = 49672.0;
    LALModJulianDate(&status, &mod_julian_date, &date);
    if (status.statusCode && lalDebugLevel > 0)
      {
        fprintf(stderr,
                "TestJulianDay: LALModJulianDate() failed, line %i, %s\n",
                __LINE__, LALTESTJULIANDAYC);
        REPORTSTATUS(&status);
        return status.statusCode;
      }


    if (lalDebugLevel > 0) 
      printf("\tModified Julian Date       = %19.10f\n", mod_julian_date);

    if (fabs(mod_julian_date - ref_mod_julian_date) > julian_precision) {
      
      if (fabs(mod_julian_date - ref_mod_julian_date) > coarse_precision) {
        fprintf(stderr, "FAIL: Modified Julian Date\n");
        return FAIL_MOD_JULIAN_DATE;
      } else {
        fprintf(stderr, "WARNING: Inaccuracy in LALModJulianDate() > %e\n",
                julian_precision);
      }
    }


    /* */
    if (lalDebugLevel > 0) {
      printf("\n");
      printf("Find Julian Day/Date for 2001-May-15 02h 37m 54s:\n");
    }

    date.unixDate.tm_sec  = 54;
    date.unixDate.tm_min  = 37;
    date.unixDate.tm_hour = 2;
    date.unixDate.tm_mday = 15;
    date.unixDate.tm_mon  = 4;
    date.unixDate.tm_year = 101;
    date.unixDate.tm_wday = 1;
    date.unixDate.tm_yday = 0;
    date.unixDate.tm_isdst = 1;

    ref_julian_day = 2452045;
    LALJulianDay(&status, &julian_day, &date);
    if (status.statusCode && lalDebugLevel > 0)
      {
        fprintf(stderr,
                "TestJulianDay: LALJulianDay() failed, line %i, %s\n",
                __LINE__, LALTESTJULIANDAYC);
        REPORTSTATUS(&status);
        return status.statusCode;
      }


    if (lalDebugLevel > 0)
      printf("\tJulian Day Number          = %8d\n", julian_day);

    if (julian_day != ref_julian_day) {
      fprintf(stderr, "FAIL: Julian Day\n");
      return FAIL_JULIAN_DAY;
    }

    ref_julian_date = 2452044.609653;
    LALJulianDate(&status, &julian_date, &date);
    if (status.statusCode && lalDebugLevel > 0)
      {
        fprintf(stderr,
                "TestJulianDay: LALJulianDate() failed, line %i, %s\n",
                __LINE__, LALTESTJULIANDAYC);
        REPORTSTATUS(&status);
        return status.statusCode;
      }


    if (lalDebugLevel > 0)
      printf("\tJulian Date                = %19.10f\n", julian_date);

    if (fabs(julian_date - ref_julian_date) > julian_precision) {

      if (fabs(julian_date - ref_julian_date) > coarse_precision) {  
        fprintf(stderr, "FAIL: Julian Date\n");
        return FAIL_JULIAN_DATE;
      } else {
        fprintf(stderr, "WARNING: Inaccuracy in LALModJulianDate() > %e\n",
                julian_precision);
      }
    }


    ref_mod_julian_date = 52044.109653;
    LALModJulianDate(&status, &mod_julian_date, &date);
    if (status.statusCode && lalDebugLevel > 0)
      {
        fprintf(stderr,
                "TestJulianDay: LALModJulianDate() failed, line %i, %s\n",
                __LINE__, LALTESTJULIANDAYC);
        REPORTSTATUS(&status);
        return status.statusCode;
      }


    if (lalDebugLevel > 0) 
      printf("\tModified Julian Date       = %19.10f\n", mod_julian_date);

    if (fabs(mod_julian_date - ref_mod_julian_date) > julian_precision) {
      
      if (fabs(mod_julian_date - ref_mod_julian_date) > coarse_precision) {
        fprintf(stderr, "FAIL: Modified Julian Date\n");
        return FAIL_MOD_JULIAN_DATE;
      } else {
        fprintf(stderr, "WARNING: Inaccuracy in LALModJulianDate() > %e\n",
                julian_precision);
      }
    }

    return SUCCESS;
}

