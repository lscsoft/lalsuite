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

INT4 lalDebugLevel = 0;

NRCSID (LALTESTLEAPSECSC, "$Id$");

typedef struct gps_leap_sec {
  time_t gps;       /* GPS time when leap sec was introduced */
  INT4   tai_utc;   /* GPS-UTC at GPS time gps */
} gps_leap_sec_t;

int main(int argc, char *argv[])
{
  static LALStatus    status;
  static const gps_leap_sec_t leapsec_data[] =
    {
      {0,         19},  /* 1980-Jan-06 */
      {46828801,  20},  /* 1981-Jul-01 */
      {78364802,  21},  /* 1982-Jul-01 */
      {109900803, 22},  /* 1983-Jul-01 */
      {173059204, 23},  /* 1985-Jul-01 */
      {252028805, 24},  /* 1988-Jan-01 */
      {315187206, 25},  /* 1990-Jan-01 */
      {346723207, 26},  /* 1991-Jan-01 */
      {393984008, 27},  /* 1992-Jul-01 */
      {425520009, 28},  /* 1993-Jul-01 */
      {457056010, 29},  /* 1994-Jul-01 */
      {504489611, 30},  /* 1996-Jan-01 */
      {551750412, 31},  /* 1997-Jul-01 */
      {599184013, 32},  /* 1999-Jan-01 */
    };
  static const INT4 num_data = sizeof(leapsec_data)/sizeof(gps_leap_sec_t);
  LIGOTimeGPS         gpsTime = {0, 0};
  INT4                gps_utc = 0;
  INT4                tai_utc = 0;
  INT4                i = 0;
  LALLeapSecFormatAndAcc formatAndAcc = {LALLEAPSEC_GPSUTC, LALLEAPSEC_STRICT};


  if (argc > 1)
    lalDebugLevel = atoi(argv[1]);


  /*
   * To run the test, compute GPS-UTC just before and just after the time
   * that the leap second is added
   */

  if (lalDebugLevel > 0)
    printf("TestLeapSecs: TESTING GPS-UTC\n");

  gpsTime.gpsSeconds = 0;
  LALLeapSecs(&status, &gps_utc, &gpsTime, &formatAndAcc);

  if (lalDebugLevel > 0)
    printf("\tGPS = %9d;    GPS-UTC = %d\n", gpsTime.gpsSeconds,
           gps_utc);

  if (status.statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr, "TestLeapSecs: LALLeapSecs() failed, line %i, %s\n",
              __LINE__, LALTESTLEAPSECSC);
      REPORTSTATUS(&status);
      return status.statusCode;
    }

  if (gps_utc != 0)
    {
      if (lalDebugLevel > 0)
        {
          fprintf(stderr,
                  "TestLeapSecs: LALLeapSecs returned wrong value: expected 0, got %d; i = %d\n",
                  gps_utc, i);
          REPORTSTATUS(&status);
        }
      return 1;
    }

  gpsTime.gpsSeconds = 1;
  LALLeapSecs(&status, &gps_utc, &gpsTime, &formatAndAcc);
  if (lalDebugLevel > 0)
    printf("\tGPS = %9d;    GPS-UTC = %d\n", gpsTime.gpsSeconds,
           gps_utc);

  if (status.statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr, "TestLeapSecs: LALLeapSecs() failed, line %i, %s\n",
              __LINE__, LALTESTLEAPSECSC);
      REPORTSTATUS(&status);
      return status.statusCode;
    }

  if (gps_utc != 0)
    {
      if (lalDebugLevel > 0)
        {
          fprintf(stderr,
                  "TestLeapSecs: LALLeapSecs returned wrong value: expected 0, got %d; i = %d\n",
                  gps_utc, i);
          REPORTSTATUS(&status);
        }
      return 1;
    }



  for (i = 1; i < num_data; ++i)
    {
      if (lalDebugLevel > 2)
        {
          printf("TestLeapSecs: BEFORE LEAP SECOND ADDED\n");
        }

      gpsTime.gpsSeconds = leapsec_data[i].gps - 1;
      LALLeapSecs(&status, &gps_utc, &gpsTime, &formatAndAcc);
      if (lalDebugLevel > 0)
        printf("\tGPS = %9d;    GPS-UTC = %d\n", gpsTime.gpsSeconds,
               gps_utc);

      if (status.statusCode && lalDebugLevel > 0)
        {
          fprintf(stderr, "TestLeapSecs: LALLeapSecs() failed, line %i, %s\n",
                  __LINE__, LALTESTLEAPSECSC);
          REPORTSTATUS(&status);
          return status.statusCode;
        }

      if (gps_utc != leapsec_data[i-1].tai_utc - 19)
        {
          if (lalDebugLevel > 0)
            {
              fprintf(stderr,
                      "TestLeapSecs: LALLeapSecs returned wrong value: expected %d, got %d; i = %d\n",
                      leapsec_data[i-1].tai_utc - 19, gps_utc, i);
              REPORTSTATUS(&status);
            }
          return 1;
        }

      if (lalDebugLevel > 2)
        {
          printf("TestLeapSecs: AFTER LEAP SECOND ADDED\n");
        }

      gpsTime.gpsSeconds += 2;
      LALLeapSecs(&status, &gps_utc, &gpsTime, &formatAndAcc);
      if (lalDebugLevel > 0)
        printf("\tGPS = %9d;    GPS-UTC = %d\n\n", gpsTime.gpsSeconds,
               gps_utc);

      if (status.statusCode && lalDebugLevel > 0)
        {
          fprintf(stderr, "TestLeapSecs: LALLeapSecs() failed, line %i, %s\n",
                  __LINE__, LALTESTLEAPSECSC);
          REPORTSTATUS(&status);
          return status.statusCode;
        }

      if (gps_utc != leapsec_data[i].tai_utc - 19)
        {
          if (lalDebugLevel > 0)
            {
              fprintf(stderr,
                      "TestLeapSecs: LALLeapSecs returned wrong value: expected %d, got %d; i = %d\n",
                      leapsec_data[i-1].tai_utc - 19, gps_utc, i);
              REPORTSTATUS(&status);
            }
          return 1;
        }
    } /* for */


  formatAndAcc.accuracy = LALLEAPSEC_LOOSE;
  gpsTime.gpsSeconds = 701654353 + 4;

  LALLeapSecs(&status, &gps_utc, &gpsTime, &formatAndAcc);
  if (lalDebugLevel > 0)
    printf("\tGPS = %9d;    GPS-UTC = %d\n", gpsTime.gpsSeconds,
           gps_utc);

  if (status.statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr, "TestLeapSecs: LALLeapSecs() failed, line %i, %s\n",
              __LINE__, LALTESTLEAPSECSC);
      REPORTSTATUS(&status);
      return status.statusCode;
    }

  if (gps_utc != 13)
    {
      if (lalDebugLevel > 0)
        {
          fprintf(stderr,
                  "TestLeapSecs: LALLeapSecs returned wrong value: expected 13, got %d; i = %d\n",
                  gps_utc, i);
          REPORTSTATUS(&status);
        }
      return 1;
    }

  if (lalDebugLevel > 0)
    printf("\nTestLeapSecs: TESTING TAI-UTC\n\n");

  /*
   * Test TAI-UTC
   */
  formatAndAcc.format = LALLEAPSEC_TAIUTC;

  for (i = 1; i < num_data; ++i)
    {
      if (lalDebugLevel > 2)
        {
          printf("TestLeapSecs: BEFORE LEAP SECOND ADDED\n");
        }

      gpsTime.gpsSeconds = leapsec_data[i].gps - 1;
      LALLeapSecs(&status, &tai_utc, &gpsTime, &formatAndAcc);
      if (lalDebugLevel > 0)
        printf("\tGPS = %9d;    TAI-UTC = %d\n", gpsTime.gpsSeconds,
               tai_utc);

      if (status.statusCode && lalDebugLevel > 0)
        {
          fprintf(stderr, "TestLeapSecs: LALLeapSecs() failed, line %i, %s\n",
                  __LINE__, LALTESTLEAPSECSC);
          REPORTSTATUS(&status);
          return status.statusCode;
        }

      if (tai_utc != leapsec_data[i-1].tai_utc)
        {
          if (lalDebugLevel > 0)
            {
              fprintf(stderr,
                      "TestLeapSecs: LALLeapSecs returned wrong value: expected %d, got %d; i = %d\n",
                      leapsec_data[i-1].tai_utc, tai_utc, i);
              REPORTSTATUS(&status);
            }
          return 1;
        }

      if (lalDebugLevel > 2)
        {
          printf("TestLeapSecs: AFTER LEAP SECOND ADDED\n");
        }

      gpsTime.gpsSeconds += 2;
      LALLeapSecs(&status, &tai_utc, &gpsTime, &formatAndAcc);
      if (lalDebugLevel > 0)
        printf("\tGPS = %9d;    TAI-UTC = %d\n\n", gpsTime.gpsSeconds,
              tai_utc);

      if (status.statusCode && lalDebugLevel > 0)
        {
          fprintf(stderr, "TestLeapSecs: LALLeapSecs() failed, line %i, %s\n",
                  __LINE__, LALTESTLEAPSECSC);
          REPORTSTATUS(&status);
          return status.statusCode;
        }

      if (tai_utc != leapsec_data[i].tai_utc)
        {
          if (lalDebugLevel > 0)
            {
              fprintf(stderr,
                      "TestLeapSecs: LALLeapSecs returned wrong value: expected %d, got %d; i = %d\n",
                      leapsec_data[i-1].tai_utc, tai_utc, i);
              REPORTSTATUS(&status);
            }
          return 1;
        }
    } /* for */

  return 0;
}
