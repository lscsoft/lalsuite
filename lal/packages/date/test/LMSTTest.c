/*
*  Copyright (C) 2007 David Chin, Jolien Creighton
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>


#define SUCCESS               0
#define TESTLMSTC_DATESTRING  1

#if 0
static BOOLEAN mstdate_ok_p(const LALDate *p_date,
                            const LALDate *p_expected_date);
#endif

/*int main(int argc, char *argv[])*/
int main(void)
{
/* FIXME:  this needs to be ported to test the XLAL replacements */
#if 0
    static LALStatus stat;
    LIGOTimeGPS      gpstime;
    /* REAL8            gmsthours, lmsthours; */
    REAL8            gmstsecs;
    LALDate          date;
    LALDate          mstdate;
    LALDate          expected_mstdate;
    LALDetector      detector;
    LALPlaceAndDate  place_and_date;
    /* LALPlaceAndGPS   place_and_gps; */
    char             refstr[256];
    char             tmpstr[128];
    LALLeapSecAccuracy accuracy = LALLEAPSEC_STRICT;


    if (argc > 1)

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

    expected_mstdate.residualNanoSeconds = 0.2738 * 1.e9;
    expected_mstdate.unixDate.tm_sec  = 21;
    expected_mstdate.unixDate.tm_min  = 39;
    expected_mstdate.unixDate.tm_hour =  3;

    if  (!mstdate_ok_p(&mstdate, &expected_mstdate))
      {
        fprintf(stderr, "TestLMST: ERROR: wrong sidereal time:\n");

        sprintf(refstr, "  %02dh %02dm %02d", mstdate.unixDate.tm_hour,
                mstdate.unixDate.tm_min, mstdate.unixDate.tm_sec);
        sprintf(tmpstr, "%.4fs", mstdate.residualNanoSeconds * 1.e-9);
        strcat(refstr, tmpstr+1);
        fprintf(stderr, "  GMST at 0h UTC on 1994-11-16:\n");
        fprintf(stderr, "       get: %s\n", refstr);

        sprintf(tmpstr, "  03h 39m 21.2738s");
        fprintf(stderr, "    expect: %s\n", tmpstr);
        fprintf(stderr, "  but since we don't have the equation of equinoxes in, can\n");
        fprintf(stderr, "  expect up to one sidereal second disagreement\n");

        return (1);
      }

    if (lalDebugLevel > 2)
      {
        sprintf(refstr, "%02dh %02dm %02d", mstdate.unixDate.tm_hour,
                mstdate.unixDate.tm_min, mstdate.unixDate.tm_sec);
        sprintf(tmpstr, "%.4fs", mstdate.residualNanoSeconds * 1.e-9);
        strcat(refstr, tmpstr+1);
        printf("GMST at 0h UTC on 1994-11-16:\n");
        printf("     get: %s\n", refstr);

        sprintf(tmpstr, "03h 39m 21.2738s");
        printf("  expect: %s\n", tmpstr);
        printf("but since we don't have the equation of equinoxes in, can\n");
        printf("expect up to one sidereal second disagreement\n");
      }

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

    expected_mstdate.residualNanoSeconds = 0.9812*1.e9;
    expected_mstdate.unixDate.tm_sec  =  5;
    expected_mstdate.unixDate.tm_min  =  5;
    expected_mstdate.unixDate.tm_hour = 21;

    if  (!mstdate_ok_p(&mstdate, &expected_mstdate))
      {
        fprintf(stderr, "TestLMST: ERROR: wrong sidereal time:\n");

        sprintf(refstr, "  %02dh %02dm %02d", mstdate.unixDate.tm_hour,
                mstdate.unixDate.tm_min, mstdate.unixDate.tm_sec);
        sprintf(tmpstr, "%.4fs", mstdate.residualNanoSeconds * 1.e-9);
        strcat(refstr, tmpstr+1);
        fprintf(stderr, "  GMST at 0h UTC on 1994-08-08:\n");
        fprintf(stderr, "       get: %s\n", refstr);

        sprintf(tmpstr, "  21h 05m 05.9812s");
        fprintf(stderr, "    expect: %s\n", tmpstr);
        fprintf(stderr, "  but since we don't have the equation of equinoxes in, can\n");
        fprintf(stderr, "  expect up to one sidereal second disagreement\n");

        return (1);
      }

    if (lalDebugLevel > 2)
      {
        sprintf(refstr, "%02dh %02dm %02d", mstdate.unixDate.tm_hour,
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
      }


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

    expected_mstdate.residualNanoSeconds = 0.5980*1.e9;
    expected_mstdate.unixDate.tm_sec  = 21;
    expected_mstdate.unixDate.tm_min  = 20;
    expected_mstdate.unixDate.tm_hour =  7;

    if  (!mstdate_ok_p(&mstdate, &expected_mstdate))
      {
        fprintf(stderr, "TestLMST: ERROR: wrong sidereal time:\n");

        sprintf(refstr, "  %02dh %02dm %02d", mstdate.unixDate.tm_hour,
                mstdate.unixDate.tm_min, mstdate.unixDate.tm_sec);
        sprintf(tmpstr, "%.4fs", mstdate.residualNanoSeconds * 1.e-9);
        strcat(refstr, tmpstr+1);
        fprintf(stderr, "  GMST at 0h UTC on 2003-01-11:\n");
        fprintf(stderr, "       get: %s\n", refstr);

        sprintf(tmpstr, "  07h 20m 21.5980s");
        fprintf(stderr, "    expect: %s\n", tmpstr);
        fprintf(stderr, "  but since we don't have the equation of equinoxes in, can\n");
        fprintf(stderr, "  expect up to one sidereal second disagreement\n");

        return (1);
      }



    if (lalDebugLevel > 2)
      {
        sprintf(refstr, "%02dh %02dm %02d", mstdate.unixDate.tm_hour,
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
      }



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


    expected_mstdate.residualNanoSeconds = 0.6093*1.e9;
    expected_mstdate.unixDate.tm_sec  =  2;
    expected_mstdate.unixDate.tm_min  = 34;
    expected_mstdate.unixDate.tm_hour = 14;

    if  (!mstdate_ok_p(&mstdate, &expected_mstdate))
      {
        fprintf(stderr, "TestLMST: ERROR: wrong sidereal time:\n");

        sprintf(refstr, "  %02dh %02dm %02d", mstdate.unixDate.tm_hour,
                mstdate.unixDate.tm_min, mstdate.unixDate.tm_sec);
        sprintf(tmpstr, "%.4fs", mstdate.residualNanoSeconds * 1.e-9);
        strcat(refstr, tmpstr+1);
        fprintf(stderr, "  GMST at 0h UTC on 2003-04-01:\n");
        fprintf(stderr, "       get: %s\n", refstr);

        sprintf(tmpstr, "  14h 34m 02.6093s ");
        fprintf(stderr, "    expect: %s\n", tmpstr);
        fprintf(stderr, "  but since we don't have the equation of equinoxes in, can\n");
        fprintf(stderr, "  expect up to one sidereal second disagreement\n");

        return (1);
      }


    if (lalDebugLevel > 2)
      {
        sprintf(refstr, "%02dh %02dm %02d", mstdate.unixDate.tm_hour,
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
      }


    return stat.statusCode;
#else
    return 0;
#endif
}

#if 0
/* allow up to 1 sidereal second difference */
static BOOLEAN mstdate_ok_p(const LALDate *p_date,
                            const LALDate *p_expected_date)
{
  BOOLEAN secs_ok_p;
  double  secs, expected_secs;

  secs = (double)((*p_date).unixDate.tm_sec) +
    (*p_date).residualNanoSeconds / 1.e9;
  expected_secs = (double)((*p_expected_date).unixDate.tm_sec) +
    (*p_expected_date).residualNanoSeconds / 1.e9;

#if 0
  printf("         secs = % 20.14e\n", secs);
  printf("expected_secs = % 20.14e\n", expected_secs);
#endif

  secs_ok_p = (fabs(secs - expected_secs) < 1.);

  return (secs_ok_p &&
          (*p_date).unixDate.tm_min == (*p_expected_date).unixDate.tm_min &&
          (*p_date).unixDate.tm_hour == (*p_expected_date).unixDate.tm_hour);
}
#endif
