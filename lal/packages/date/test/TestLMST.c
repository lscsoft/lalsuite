#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>

INT4 lalDebugLevel = 2;

NRCSID (TESTLMSTC, "$Id$");

#define TESTLMSTC_OK          0
#define TESTLMSTC_DATESTRING  1

int main(int argc, char *argv[])
{
    static LALStatus stat;
    LIGOTimeGPS      gpstime;
    LIGOTimeUnix     unixtime;
    REAL8            gmsthours, lmsthours;
    REAL8            gmstsecs;
    LALDate          date;
    LALDate          mstdate;
    LALDetector      detector;
    LALPlaceAndDate  place_and_date;
    LALPlaceAndGPS   place_and_gps;
    char             refstr[128];
    char             tmpstr[128];
    CHARVector      *timestamp = NULL;

    /*
     * Allocate the CHARVector
     */
    LALCHARCreateVector(&stat, &timestamp, (UINT4)128);

    if (1) {
      unixtime.unixSeconds     = 60858;
      unixtime.unixNanoSeconds =     0;

      LALUtime(&stat, &date, &unixtime);
      LALDateString(&stat, timestamp, &date);

      sprintf(refstr, "1970-01-01 16:54:18 UTC Thu");

      if (lalDebugLevel > 0) {
        printf("refstr  = %s\n", refstr);
        printf("timestamp->data = %s\n", timestamp->data);
      }

      if (strcmp(refstr, timestamp->data) != 0) {
        fprintf(stderr, "LALDateString() not working right. Quitting.\n");
        return TESTLMSTC_DATESTRING;
      }
    }
  
    /*
     * Check against the Astronomical Almanac:
     * For 1994-11-16 0h UT - Julian Date 2449672.5, GMST 03h 39m 21.2738s
     */
    date.unixDate.tm_sec  =  0;
    date.unixDate.tm_min  =  0;
    date.unixDate.tm_hour =  0;
    date.unixDate.tm_mday = 16;
    date.unixDate.tm_mon  = 10;
    date.unixDate.tm_year = 94;

    /* The corresponding GPS time is */
    gpstime.gpsSeconds     = 468979210;
    gpstime.gpsNanoSeconds =         0;

    detector.frDetector.vertexLongitudeDegrees = 0.;  /* Greenwich */
    place_and_date.p_detector = &detector;
    place_and_date.p_date     = &date;
    LALGMST1(&stat, &gmsthours, &date, MST_HRS);
    LALLMST1(&stat, &lmsthours, &place_and_date, MST_HRS);

    LALGMST1(&stat, &gmstsecs, &date, MST_SEC);
    LALSecsToLALDate(&stat, &mstdate, gmstsecs);
    strftime(refstr, 64, "%Hh %Mm %S", &(mstdate.unixDate));
    sprintf(tmpstr, "%fs", mstdate.residualNanoSeconds * 1.e-9);
    strcat(refstr, tmpstr+1); /* remove leading 0 */
    printf("refstr = %s\n", refstr);
    
    LALDateString(&stat, timestamp, &date);

    printf("timestamp->data: %s\n", timestamp->data);

    return TESTLMSTC_OK;
}
