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

#include <stdlib.h>
#include <time.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/Date.h>


/*int main(int argc, char *argv[])*/
int main(void)
{
/* FIXME:  this needs to be ported to test the XLAL replacements */
#if 0
    static LALStatus status;
    LALDate          date;
    LALDate          mstdate;
    REAL8            gmsthours, lmsthours;
    REAL8            gmstsecs;
    LALDetector      detector;
    LALPlaceAndDate  place_and_date;
    LALPlaceAndGPS   place_and_gps;
    CHAR             timestr[128];
    CHAR             tmpstr[128];
    CHARVector      *datestamp = NULL;
    INT4             testid;
    INT4             dayofmonth, monthofyear;
    LIGOTimeGPS      gpstime;
    LALLeapSecAccuracy accuracy = LALLEAPSEC_STRICT;
    LALMSTUnitsAndAcc units_and_acc;


    if (argc != 3)
    {
        /*
         * Print help message and exit
         */
        printf("Usage: TestLMST testid debug_level\n");
        printf("              testid      = [1,2]\n");
        printf("              debug_level = [0,1,2]\n");
        return 0;
    }

    testid        = atoi(argv[1]);

    LALCHARCreateVector(&status, &datestamp, (UINT4)64);

    /*
     * Initialize date structures. Prevents strftime() from
     * dumping core under Solaris 2.7/gcc-2.95.2
     */
    date.unixDate.tm_sec  = 0;
    date.unixDate.tm_min  = 0;
    date.unixDate.tm_hour = 0;
    date.unixDate.tm_mday = 0;
    date.unixDate.tm_mon  = 0;
    date.unixDate.tm_year = 0;
    date.unixDate.tm_wday = 0;
    date.unixDate.tm_yday = 0;
    date.unixDate.tm_isdst = 0;

    date.unixDate.tm_sec  = 0;
    date.unixDate.tm_min  = 0;
    date.unixDate.tm_hour = 0;
    date.unixDate.tm_mday = 0;
    date.unixDate.tm_mon  = 0;
    date.unixDate.tm_year = 0;
    date.unixDate.tm_wday = 0;
    date.unixDate.tm_yday = 0;
    date.unixDate.tm_isdst = 0;

    printf("TEST of LALGMST1 routine\n");
    printf("========================\n");

    if (testid == 1)
    {
        /*
         * Check against the Astronomical Almanac:
         * For 1994-11-16 0h UT - Julian Date 2449672.5, GMST 03h 39m 21.2738s
         */
        date.unixDate.tm_sec  = 0;
        date.unixDate.tm_min  = 0;
        date.unixDate.tm_hour = 0;
        date.unixDate.tm_mday = 16;
        date.unixDate.tm_mon  = 10;
        date.unixDate.tm_year = 94;

        /* all other fields to 0 */
        date.unixDate.tm_wday = 0;
        date.unixDate.tm_yday = 0;
        date.unixDate.tm_isdst = 0;

        /* The corresponding GPS time is */
        gpstime.gpsSeconds     = 468979210;
        gpstime.gpsNanoSeconds = 0;

        detector.frDetector.vertexLongitudeRadians = 0.;  /* Greenwich */
        place_and_date.p_detector = &detector;
        place_and_date.p_date     = &date;
        LALGMST1(&status, &gmsthours, &date, MST_HRS);
        LALLMST1(&status, &lmsthours, &place_and_date, MST_HRS);

        LALGMST1(&status, &gmstsecs, &date, MST_SEC);
        LALSecsToLALDate(&status, &mstdate, gmstsecs);
        strftime(timestr, 64, "%Hh %Mm %S", &(mstdate.unixDate));
        sprintf(tmpstr, "%fs", mstdate.residualNanoSeconds * 1.e-9);
        strcat(timestr, tmpstr+1); /* remove leading 0 */
        LALDateString(&status, datestamp, &date);

        printf("     Time = %s\n", datestamp->data);
        printf("gmsthours = %f = %s\n", gmsthours, timestr);
        printf("    expect: 3.655728 = 03h 39m 20.6222s \n");
        printf("lmsthours = %f\n", lmsthours);

        printf("\n");
        printf("Using the GPStoGMST1() and GPStoLMST1() routines instead:\n");
        units_and_acc.units = MST_SEC;
        units_and_acc.accuracy = LALLEAPSEC_STRICT;
        LALGPStoGMST1(&status, &gmstsecs, &gpstime, &units_and_acc);
        units_and_acc.units = MST_HRS;
        LALGPStoGMST1(&status, &gmsthours, &gpstime, &units_and_acc);
        LALSecsToLALDate(&status, &mstdate, gmstsecs);
        strftime(timestr, 64, "%Hh %Mm %S", &(mstdate.unixDate));
        sprintf(tmpstr, "%fs", mstdate.residualNanoSeconds * 1.e-9);
        strcat(timestr, tmpstr+1); /* remove leading 0 */
        /* LALDateString(&status, datestamp, &date); */

        place_and_gps.p_detector = &detector;
        place_and_gps.p_gps      = &gpstime;
        units_and_acc.units      = MST_HRS;
        LALGPStoLMST1(&status, &lmsthours, &place_and_gps, &units_and_acc);
        printf("     Time = %s\n", datestamp->data);
        printf("gmsthours = %f = %s\n", gmsthours, timestr);
        printf("    expect: 3.655728 = 03h 39m 20.6222s \n");
        printf("lmsthours = %f\n", lmsthours);


        /* Doing the GPS to MST1 conversions the long way */
        printf("\n");
        printf("\nConverting GPS to LMST the long way instead:\n");
        LALGPStoUTC(&status, &date, &gpstime, &accuracy);
        LALGMST1(&status, &gmsthours, &date, MST_HRS);
        LALLMST1(&status, &lmsthours, &place_and_date, MST_HRS);

        LALGMST1(&status, &gmstsecs, &date, MST_SEC);
        LALSecsToLALDate(&status, &mstdate, gmstsecs);
        strftime(timestr, 64, "%Hh %Mm %S", &(mstdate.unixDate));
        sprintf(tmpstr, "%fs", mstdate.residualNanoSeconds * 1.e-9);
        strcat(timestr, tmpstr+1); /* remove leading 0 */
        LALDateString(&status, datestamp, &date);

        printf("     Time = %s\n", datestamp->data);
        printf("gmsthours = %f = %s\n", gmsthours, timestr);
        printf("    expect: 3.655728 = 03h 39m 20.6222s \n");
        printf("lmsthours = %f\n", lmsthours);


        /*
         * Another check against the Almanac:
         * For 1994-08-17 2h 19m 03.0736s -> 0h GMST
         */
        printf("\n* * * * * * * * * * * * * * * * * *\n\n");
        date.residualNanoSeconds = 73600000;
        date.unixDate.tm_sec = 3;
        date.unixDate.tm_min = 19;
        date.unixDate.tm_hour = 2;
        date.unixDate.tm_mday = 17;
        date.unixDate.tm_mon  = 7;
        date.unixDate.tm_year = 94;

        gpstime.gpsSeconds = 461125153;
        gpstime.gpsNanoSeconds = date.residualNanoSeconds;

        detector.frDetector.vertexLongitudeRadians = 0.;  /* Greenwich */
        place_and_date.p_detector = &detector;
        place_and_date.p_date     = &date;
        LALGMST1(&status, &gmsthours, &date, MST_HRS);
        LALLMST1(&status, &lmsthours, &place_and_date, MST_HRS);

        LALGMST1(&status, &gmstsecs, &date, MST_SEC);
        LALSecsToLALDate(&status, &mstdate, gmstsecs);
        strftime(timestr, 64, "%Hh %Mm %S", &(mstdate.unixDate));
        sprintf(tmpstr, "%fs", mstdate.residualNanoSeconds * 1.e-9);
        strcat(timestr, tmpstr+1); /* remove leading 0 */
        /* LALDateString(&status, datestamp, &date); */

        printf("     Time = %s\n", datestamp->data);
        printf("gmsthours = %f = %s\n", gmsthours, timestr);
        printf("    expect:        0h = 00h 00m 00s \n");

        printf("\nUsing the GPStoGMST1() and GPStoLMST1() routines instead:\n");
        units_and_acc.units = MST_SEC;
        LALGPStoGMST1(&status, &gmstsecs, &gpstime, &units_and_acc);
        LALSecsToLALDate(&status, &mstdate, gmstsecs);
        strftime(timestr, 64, "%Hh %Mm %S", &(mstdate.unixDate));
        sprintf(tmpstr, "%fs", mstdate.residualNanoSeconds * 1.e-9);
        strcat(timestr, tmpstr+1); /* remove leading 0 */
        /* LALDateString(&status, datestamp, &date); */

        place_and_gps.p_detector = &detector;
        place_and_gps.p_gps      = &gpstime;
        units_and_acc.units      = MST_HRS;
        LALGPStoLMST1(&status, &lmsthours, &place_and_gps, &units_and_acc);
        printf("     Time = %s\n", datestamp->data);
        printf("gmsthours = %f = %s\n", gmsthours, timestr);
        printf("    expect:        0h = 00h 00m 00s \n");
        printf("lmsthours = %f\n", lmsthours);

        /* Doing the GPS to MST1 conversions the long way */
        printf("\n");
        printf("\nConverting GPS to LMST the long way instead:\n");
        LALGPStoUTC(&status, &date, &gpstime, &accuracy);
        LALGMST1(&status, &gmsthours, &date, MST_HRS);
        LALLMST1(&status, &lmsthours, &place_and_date, MST_HRS);

        LALGMST1(&status, &gmstsecs, &date, MST_SEC);
        LALSecsToLALDate(&status, &mstdate, gmstsecs);
        strftime(timestr, 64, "%Hh %Mm %S", &(mstdate.unixDate));
        sprintf(tmpstr, "%fs", mstdate.residualNanoSeconds * 1.e-9);
        strcat(timestr, tmpstr+1); /* remove leading 0 */
        LALDateString(&status, datestamp, &date);

        printf("     Time = %s\n", datestamp->data);
        printf("gmsthours = %f = %s\n", gmsthours, timestr);
        printf("    expect:        0h = 00h 00m 00s \n");
        printf("lmsthours = %f\n", lmsthours);


        /*
         * A third check against the Almanac:
         * For 1994-05-17 0h
         */
        printf("\n* * * * * * * * * * * * * * * * * *\n\n");
        date.residualNanoSeconds = 0;
        date.unixDate.tm_sec = 0;
        date.unixDate.tm_min = 0;
        date.unixDate.tm_hour = 0;
        date.unixDate.tm_mday = 17;
        date.unixDate.tm_mon  = 4;
        date.unixDate.tm_year = 94;

        gpstime.gpsSeconds = 453168010;
        gpstime.gpsNanoSeconds = 0;

        detector.frDetector.vertexLongitudeRadians = 0.;  /* Greenwich */
        place_and_date.p_detector = &detector;
        place_and_date.p_date     = &date;

        LALGMST1(&status, &gmsthours, &date, MST_HRS);
        LALLMST1(&status, &lmsthours, &place_and_date, MST_HRS);

        LALGMST1(&status, &gmstsecs, &date, MST_SEC);
        LALSecsToLALDate(&status, &mstdate, gmstsecs);
        strftime(timestr, 64, "%Hh %Mm %S", &(mstdate.unixDate));
        sprintf(tmpstr, "%fs", mstdate.residualNanoSeconds * 1.e-9);
        strcat(timestr, tmpstr+1); /* remove leading 0 */
        /* LALDateString(&status, datestamp, &date); */

        printf("     Time = %s\n", datestamp->data);
        printf("gmsthours = %f   = %s\n", gmsthours, timestr);
        printf("    expect: 15.63105328 = 15h 37m 51.7918s\n");

        printf("\n");
        printf("Using the GPStoGMST1() and GPStoLMST1() routines instead:\n");
        units_and_acc.units = MST_SEC;
        LALGPStoGMST1(&status, &gmstsecs, &gpstime, &units_and_acc);
        LALSecsToLALDate(&status, &mstdate, gmstsecs);
        strftime(timestr, 64, "%Hh %Mm %S", &(mstdate.unixDate));
        sprintf(tmpstr, "%fs", mstdate.residualNanoSeconds * 1.e-9);
        strcat(timestr, tmpstr+1); /* remove leading 0 */
        /* LALDateString(&status, datestamp, &date); */

        place_and_gps.p_detector = &detector;
        place_and_gps.p_gps      = &gpstime;
        units_and_acc.units      = MST_HRS;
        LALGPStoLMST1(&status, &lmsthours, &place_and_gps, &units_and_acc);
        printf("     Time = %s\n", datestamp->data);
        printf("gmsthours = %f = %s\n", gmsthours, timestr);
        printf("    expect: 15.63105328 = 15h 37m 51.7918s\n");
        printf("lmsthours = %f\n", lmsthours);


        /* And increment the date by 1 hour to see
         * if it makes GMST change */
        printf("\n* * * * * * * * * * * * * * * * * *\n\n");
        date.residualNanoSeconds = 0;
        date.unixDate.tm_sec = 0;
        date.unixDate.tm_min = 0;
        date.unixDate.tm_hour = 1;
        date.unixDate.tm_mday = 17;
        date.unixDate.tm_mon  = 4;
        date.unixDate.tm_year = 94;

        gpstime.gpsSeconds += 3600;

        detector.frDetector.vertexLongitudeRadians = 0.;  /* Greenwich */
        place_and_date.p_detector = &detector;
        place_and_date.p_date     = &date;

        LALGMST1(&status, &gmsthours, &date, MST_HRS);
        LALLMST1(&status, &lmsthours, &place_and_date, MST_HRS);

        LALGMST1(&status, &gmstsecs, &date, MST_SEC);
        LALSecsToLALDate(&status, &mstdate, gmstsecs);
        strftime(timestr, 64, "%Hh %Mm %S", &(mstdate.unixDate));
        sprintf(tmpstr, "%fs", mstdate.residualNanoSeconds * 1.e-9);
        strcat(timestr, tmpstr+1); /* remove leading 0 */
        /* LALDateString(&status, datestamp, &date); */

        printf("     Time = %s\n", datestamp->data);
        printf("gmsthours = %f   = %s\n", gmsthours, timestr);
        printf("    expect: 16.63105328 = 16h 37m 51.7918s\n");

        printf("\nUsing the GPStoGMST1() and GPStoLMST1() routines instead:\n");
        units_and_acc.units = MST_SEC;
        LALGPStoGMST1(&status, &gmstsecs, &gpstime, &units_and_acc);
        LALSecsToLALDate(&status, &mstdate, gmstsecs);
        strftime(timestr, 64, "%Hh %Mm %S", &(mstdate.unixDate));
        sprintf(tmpstr, "%fs", mstdate.residualNanoSeconds * 1.e-9);
        strcat(timestr, tmpstr+1); /* remove leading 0 */
        /* LALDateString(&status, datestamp, &date); */

        place_and_gps.p_detector = &detector;
        place_and_gps.p_gps      = &gpstime;
        units_and_acc.units      = MST_HRS;
        LALGPStoLMST1(&status, &lmsthours, &place_and_gps, &units_and_acc);
        printf("     Time = %s\n", datestamp->data);
        printf("gmsthours = %f = %s\n", gmsthours, timestr);
        printf("    expect: 16.63105328 = 16h 37m 51.7918s\n");
        printf("lmsthours = %f\n", lmsthours);




        /*
         * A fourth check against the Almanac:
         * For 1994-05-17 08:20:46.7448 -> 0h GMST
         */
        printf("\n* * * * * * * * * * * * * * * * * *\n\n");
        date.residualNanoSeconds = 744800000;
        date.unixDate.tm_sec = 46;
        date.unixDate.tm_min = 20;
        date.unixDate.tm_hour = 8;
        date.unixDate.tm_mday = 17;
        date.unixDate.tm_mon  = 4;
        date.unixDate.tm_year = 94;

        gpstime.gpsSeconds = 453168010;
        gpstime.gpsNanoSeconds = date.residualNanoSeconds;

        detector.frDetector.vertexLongitudeRadians = 0.;  /* Greenwich */
        place_and_date.p_detector = &detector;
        place_and_date.p_date     = &date;

        LALGMST1(&status, &gmsthours, &date, MST_HRS);
        LALLMST1(&status, &lmsthours, &place_and_date, MST_HRS);

        LALGMST1(&status, &gmstsecs, &date, MST_SEC);
        LALSecsToLALDate(&status, &mstdate, gmstsecs);
        strftime(timestr, 64, "%Hh %Mm %S", &(mstdate.unixDate));
        sprintf(tmpstr, "%fs", mstdate.residualNanoSeconds * 1.e-9);
        strcat(timestr, tmpstr+1); /* remove leading 0 */
        /* LALDateString(&status, datestamp, &date); */

        printf("     Time = %s\n", datestamp->data);
        printf("gmsthours = %f = %s\n", gmsthours, timestr);
        printf("    expect:        0h = 00h 00m 00s \n");
    }

    if (testid == 2)
    {
        /* Generate column G (Mean) of Sec. B of the Almanac for 1994 */
        date.unixDate.tm_sec = 0;
        date.unixDate.tm_min = 0;
        date.unixDate.tm_hour = 0;
        date.unixDate.tm_year = 94;

        /* all other fields to 0 */
        date.unixDate.tm_wday = 0;
        date.unixDate.tm_yday = 0;
        date.unixDate.tm_isdst = 0;

        printf("\nGMST1 of 0h UT1 for 1994 (check against Almanac):\n");
        for (monthofyear = 0; monthofyear < 12; ++monthofyear)
        {
            date.unixDate.tm_mon = monthofyear;

            for (dayofmonth = 1; dayofmonth < 32; ++dayofmonth)
            {
                date.unixDate.tm_mday = dayofmonth;

                LALGMST1(&status, &gmstsecs, &date, MST_SEC);
                LALSecsToLALDate(&status, &mstdate, gmstsecs);
                strftime(timestr, 64, "%Hh %Mm %S", &(mstdate.unixDate));
                sprintf(tmpstr, "%fs", mstdate.residualNanoSeconds * 1.e-9);
                strcat(timestr, tmpstr+1); /* remove leading 0 */
                LALDateString(&status, datestamp, &date);
                printf("%s: %s\n", datestamp->data, timestr);

                /* February */
                if (monthofyear == 1)
                    if (dayofmonth == 28)
                        break;
                /* the 30-day months */
                if (monthofyear == 3 || monthofyear == 5 ||
                    monthofyear == 8 || monthofyear == 10)
                    if (dayofmonth == 30)
                        break;
            }
        }
    }

    /*
     * Housekeeping
     */
    LALCHARDestroyVector(&status, &datestamp);
#endif
    return(0);
}
