#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>

INT4 lalDebugLevel = 0;

NRCSID (LALTESTLMSTC, "$Id$");

#define SUCCESS               0
#define TESTLMSTC_DATESTRING  1

int main(int argc, char *argv[])
{
    static LALStatus stat;
    LIGOTimeGPS      gpstime;
    /* REAL8            gmsthours, lmsthours; */
    REAL8            gmstsecs;
    LALDate          date;
    LALDate          mstdate;
    LALDetector      detector;
    LALPlaceAndDate  place_and_date;
    /* LALPlaceAndGPS   place_and_gps; */
    char             refstr[256];
    char             tmpstr[128];
    LALLeapSecAccuracy accuracy = LALLEAPSEC_STRICT;


    if (argc == 2)
      lalDebugLevel = atoi(argv[1]);

    /** TEST NO. 1 **/
    
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
    date.residualNanoSeconds = 0;

    /* The corresponding GPS time is */
    gpstime.gpsSeconds     = 468979210;
    gpstime.gpsNanoSeconds =         0;

    LALUTCtoGPS(&stat, &gpstime, &date, &accuracy);
    if (stat.statusCode && lalDebugLevel > 0)
      {
        fprintf(stderr, "TestLMST: LALUTCtoGPS() failed, line %i, %s\n",
                __LINE__, LALTESTLMSTC);
        REPORTSTATUS(&stat);
        return stat.statusCode;
      }

    detector.frDetector.vertexLongitudeRadians = 0.;  /* Greenwich */
    place_and_date.p_detector = &detector;
    place_and_date.p_date     = &date;

    LALGMST1(&stat, &gmstsecs, &date, MST_SEC);
    if (stat.statusCode && lalDebugLevel > 0)
      {
        fprintf(stderr, "TestLMST: LALGMST1() failed, line %i, %s\n",
                __LINE__, LALTESTLMSTC);
        REPORTSTATUS(&stat);
        return stat.statusCode;
      }
    
    LALSecsToLALDate(&stat, &mstdate, gmstsecs);
    if (stat.statusCode && lalDebugLevel > 0)
      {
        fprintf(stderr, "TestLMST: LALSecsToLALDate() failed, line %i, %s\n",
                __LINE__, LALTESTLMSTC);
        REPORTSTATUS(&stat);
        return stat.statusCode;
      }

    /*strftime(refstr, 64, "%Hh %Mm %S", &(mstdate.unixDate)); */
    /*sprintf(tmpstr, "%.4fs", mstdate.residualNanoSeconds * 1.e-9);*/
    /*strcat(refstr, tmpstr+1); / * remove leading 0 */
    /*printf("   get: %s\n", refstr);*/

    sprintf(refstr, "%2dh %2dm %2d", mstdate.unixDate.tm_hour,
            mstdate.unixDate.tm_min, mstdate.unixDate.tm_sec);
    sprintf(tmpstr, "%.4fs", mstdate.residualNanoSeconds * 1.e-9);
    strcat(refstr, tmpstr+1);
    printf("GMST at 0h UTC on 1994-11-16:\n");
    printf("     get: %s\n", refstr);

    sprintf(tmpstr, "03h 39m 21.2738s");
    printf("  expect: %s\n", tmpstr);
    printf("but since we don't have the equation of equinoxes in, can\n");
    printf("expect up to one sidereal second disagreement\n");

    /** TEST NO. 2 **/

    /*
     * Check against the Astronomical Almanac:
     * For 1994-08-08 0h UT - Julian Date 2449572.5, GMST 21h 05m 05.9812s
     */
    date.unixDate.tm_sec  =   0;
    date.unixDate.tm_min  =   0;
    date.unixDate.tm_hour =   0;
    date.unixDate.tm_mday =   8;
    date.unixDate.tm_mon  =   7;
    date.unixDate.tm_year =  94;
    date.residualNanoSeconds = 0;

    /* The corresponding GPS time is */
    gpstime.gpsSeconds     = 460339210;
    gpstime.gpsNanoSeconds =         0;

    LALUTCtoGPS(&stat, &gpstime, &date, &accuracy);
    if (stat.statusCode && lalDebugLevel > 0)
      {
        fprintf(stderr, "TestLMST: LALUTCtoGPS() failed, line %i, %s\n",
                __LINE__, LALTESTLMSTC);
        REPORTSTATUS(&stat);
        return stat.statusCode;
      }

    detector.frDetector.vertexLongitudeRadians = 0.;  /* Greenwich */
    place_and_date.p_detector = &detector;
    place_and_date.p_date     = &date;

    LALGMST1(&stat, &gmstsecs, &date, MST_SEC);
    if (stat.statusCode && lalDebugLevel > 0)
      {
        fprintf(stderr, "TestLMST: LALGMST1() failed, line %i, %s\n",
                __LINE__, LALTESTLMSTC);
        REPORTSTATUS(&stat);
        return stat.statusCode;
      }
    
    LALSecsToLALDate(&stat, &mstdate, gmstsecs);
    if (stat.statusCode && lalDebugLevel > 0)
      {
        fprintf(stderr, "TestLMST: LALSecsToLALDate() failed, line %i, %s\n",
                __LINE__, LALTESTLMSTC);
        REPORTSTATUS(&stat);
        return stat.statusCode;
      }

    /* strftime(refstr, 64, "%Hh %Mm %S", &(mstdate.unixDate)); */
    /* sprintf(tmpstr, "%.4fs", mstdate.residualNanoSeconds * 1.e-9); */
    /* strcat(refstr, tmpstr+1); / * remove leading 0 */
    /* printf("   get: %s\n", refstr); */

    sprintf(refstr, "%2dh %2dm %2d", mstdate.unixDate.tm_hour,
            mstdate.unixDate.tm_min, mstdate.unixDate.tm_sec);
    sprintf(tmpstr, "%.4fs", mstdate.residualNanoSeconds * 1.e-9);
    printf("\n");
    printf("GMST at 0h UTC on 1994-08-08:\n");
    strcat(refstr, tmpstr+1);
    printf("     get: %s\n", refstr);

    sprintf(tmpstr, "21h 05m 05.9812s");
    printf("  expect: %s\n", tmpstr);
    printf("but since we don't have the equation of equinoxes in, can\n");
    printf("expect up to one sidereal second disagreement\n");



    /** TEST NO. 3 **/
    
    /*** check another date/time ***/

    /*
     * Check against the Astronomical Almanac:
     * For 2003-01-11 0h UT - Julian Date 2452650.5, GMST 7h 20m 21.5980s
     */
    date.unixDate.tm_sec  =   0;
    date.unixDate.tm_min  =   0;
    date.unixDate.tm_hour =   0;
    date.unixDate.tm_mday =  11;
    date.unixDate.tm_mon  =   0;
    date.unixDate.tm_year = 103;
    date.residualNanoSeconds = 0;

    /* The corresponding GPS time is */
    gpstime.gpsSeconds     = 726278413;
    gpstime.gpsNanoSeconds =         0;

    LALUTCtoGPS(&stat, &gpstime, &date, &accuracy);
    if (stat.statusCode && lalDebugLevel > 0)
      {
        fprintf(stderr, "TestLMST: LALUTCtoGPS() failed, line %i, %s\n",
                __LINE__, LALTESTLMSTC);
        REPORTSTATUS(&stat);
        return stat.statusCode;
      }

    detector.frDetector.vertexLongitudeRadians = 0.;  /* Greenwich */
    place_and_date.p_detector = &detector;
    place_and_date.p_date     = &date;

    LALGMST1(&stat, &gmstsecs, &date, MST_SEC);
    if (stat.statusCode && lalDebugLevel > 0)
      {
        fprintf(stderr, "TestLMST: LALGMST1() failed, line %i, %s\n",
                __LINE__, LALTESTLMSTC);
        REPORTSTATUS(&stat);
        return stat.statusCode;
      }
    
    LALSecsToLALDate(&stat, &mstdate, gmstsecs);
    if (stat.statusCode && lalDebugLevel > 0)
      {
        fprintf(stderr, "TestLMST: LALSecsToLALDate() failed, line %i, %s\n",
                __LINE__, LALTESTLMSTC);
        REPORTSTATUS(&stat);
        return stat.statusCode;
      }

    /* strftime(refstr, 64, "%Hh %Mm %S", &(mstdate.unixDate)); */
    /* sprintf(tmpstr, "%.4fs", mstdate.residualNanoSeconds * 1.e-9); */
    /* strcat(refstr, tmpstr+1); / * remove leading 0 */
    /* printf("   get: %s\n", refstr); */

    sprintf(refstr, "%2dh %2dm %2d", mstdate.unixDate.tm_hour,
            mstdate.unixDate.tm_min, mstdate.unixDate.tm_sec);
    sprintf(tmpstr, "%.4fs", mstdate.residualNanoSeconds * 1.e-9);
    strcat(refstr, tmpstr+1);
    printf("\n");
    printf("GMST at 0h UTC on 2003-01-11:\n");
    printf("     get: %s\n", refstr);

    sprintf(tmpstr, "07h 20m 21.5980s");
    printf("  expect: %s\n", tmpstr);
    printf("but since we don't have the equation of equinoxes in, can\n");
    printf("expect up to one sidereal second disagreement\n");




    /** TEST NO. 4 **/
    
    /*** check another date/time ***/

    /*
     * Check against the Astronomical Almanac:
     * For 2003-04-01 0h UT - Julian Date 2452730.5, GMST 12h 35m 45.9989s
     */
    date.unixDate.tm_sec  =   0;
    date.unixDate.tm_min  =   0;
    date.unixDate.tm_hour =   0;
    date.unixDate.tm_mday =   1;
    date.unixDate.tm_mon  =   4;
    date.unixDate.tm_year = 103;
    date.residualNanoSeconds = 0;

    /* The corresponding GPS time is */
    gpstime.gpsSeconds     = 733190413;
    gpstime.gpsNanoSeconds =         0;

    LALUTCtoGPS(&stat, &gpstime, &date, &accuracy);
    if (stat.statusCode && lalDebugLevel > 0)
      {
        fprintf(stderr, "TestLMST: LALUTCtoGPS() failed, line %i, %s\n",
                __LINE__, LALTESTLMSTC);
        REPORTSTATUS(&stat);
        return stat.statusCode;
      }

    detector.frDetector.vertexLongitudeRadians = 0.;  /* Greenwich */
    place_and_date.p_detector = &detector;
    place_and_date.p_date     = &date;

    LALGMST1(&stat, &gmstsecs, &date, MST_SEC);
    if (stat.statusCode && lalDebugLevel > 0)
      {
        fprintf(stderr, "TestLMST: LALGMST1() failed, line %i, %s\n",
                __LINE__, LALTESTLMSTC);
        REPORTSTATUS(&stat);
        return stat.statusCode;
      }
    
    LALSecsToLALDate(&stat, &mstdate, gmstsecs);
    if (stat.statusCode && lalDebugLevel > 0)
      {
        fprintf(stderr, "TestLMST: LALSecsToLALDate() failed, line %i, %s\n",
                __LINE__, LALTESTLMSTC);
        REPORTSTATUS(&stat);
        return stat.statusCode;
      }

    /* strftime(refstr, 64, "%Hh %Mm %S", &(mstdate.unixDate)); */
    /* sprintf(tmpstr, "%.4fs", mstdate.residualNanoSeconds * 1.e-9); */
    /* strcat(refstr, tmpstr+1); / * remove leading 0 */
    /* printf("   get: %s\n", refstr); */

    sprintf(refstr, "%2dh %2dm %2d", mstdate.unixDate.tm_hour,
            mstdate.unixDate.tm_min, mstdate.unixDate.tm_sec);
    sprintf(tmpstr, "%.4fs", mstdate.residualNanoSeconds * 1.e-9);
    strcat(refstr, tmpstr+1);
    printf("\n");
    printf("GMST at 0h UTC on 2003-04-01:\n");
    printf("     get: %s\n", refstr);

    sprintf(tmpstr, "14h 34m 02.6093s");
    printf("  expect: %s\n", tmpstr);
    printf("but since we don't have the equation of equinoxes in, can\n");
    printf("expect up to one sidereal second disagreement\n");
    printf("additionally, this should generate an error/warning since\n");
    printf("the requested is after when we have leap-second data for.\n");

    
    return stat.statusCode;
    
}
