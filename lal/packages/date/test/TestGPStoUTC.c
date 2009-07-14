/*
*  Copyright (C) 2007 David Chin
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

INT4 lalDebugLevel = 0;

NRCSID (LALTESTGPSTOUTCC, "$Id$");

int main(int argc, char *argv[])
{
  static LALStatus    status;
  LIGOTimeGPS         gpsTime = {0, 0};
  LIGOTimeGPS         tmpGps  = {0, 0};
  LALDate             utcDate;
  LALLeapSecAccuracy  accuracy = LALLEAPSEC_LOOSE; /* prevent ABORT */
  CHARVector         *timestamp = NULL;
  char                refstamp[128];
  /* char                infostr[256]; */

  if (argc > 1)
      lalDebugLevel = atoi(argv[1]);

  memset(&utcDate, 0, sizeof(utcDate));
  LALCHARCreateVector(&status, &timestamp, (UINT4)128);
  if (status.statusCode && lalDebugLevel > 0)
    REPORTSTATUS(&status);

  /*
   * GPS 0 == 1980-01-06 00:00:00 UTC Sun
   */
  if(!XLALGPSToUTC(&utcDate.unixDate, gpsTime.gpsSeconds))
    {
      fprintf(stderr, "TestGPStoUTC: XLALGPSToUTC() failed, line %i, %s\n",
              __LINE__, LALTESTGPSTOUTCC);
      return 1;
    }

  LALDateString(&status, timestamp, &utcDate);
  if (status.statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr, "TestGPStoUTC: LALDateString() failed, line %i, %s\n",
              __LINE__, LALTESTGPSTOUTCC);
      REPORTSTATUS(&status);
      return status.statusCode;
    }

  sprintf(refstamp, "1980-01-06 00:00:00 UTC Sun");

  if (lalDebugLevel > 0)
    {
      fprintf(stderr, "refstamp  = %s\n", refstamp);
      fprintf(stderr, "timestamp = %s\n", timestamp->data);
    }

  if (strcmp(refstamp, timestamp->data) != 0)
    {
      LALInfo(&status, "GPStoUTC conversion failed: wrong UTC result");
      fprintf(stderr, "TestGPStoUTC: date strings do not match, line %i, %s\n",
              __LINE__, LALTESTGPSTOUTCC);
      LALCHARDestroyVector(&status, &timestamp);
      if (status.statusCode && lalDebugLevel > 0)
        {
          fprintf(stderr,
                  "TestGPStoUTC: LALCHARDestroyVector() failed, line %i, %s\n",
                  __LINE__, LALTESTGPSTOUTCC);
          REPORTSTATUS(&status);
          return status.statusCode;
        }
      REPORTSTATUS(&status);
      LALCheckMemoryLeaks();
      return 1;
    }

  /*
   * GPS 457574400 == 1994-07-06 23:59:50 UTC Wed
   */
  gpsTime.gpsSeconds = 457574400;
  gpsTime.gpsNanoSeconds = 0;

  if(!XLALGPSToUTC(&utcDate.unixDate, gpsTime.gpsSeconds))
    {
      fprintf(stderr, "TestGPStoUTC: XLALGPSToUTC() failed, line %i, %s\n",
              __LINE__, LALTESTGPSTOUTCC);
      return 1;
    }

  LALDateString(&status, timestamp, &utcDate);
  if (status.statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr, "TestGPStoUTC: LALDateString() failed, line %i, %s\n",
              __LINE__, LALTESTGPSTOUTCC);
      REPORTSTATUS(&status);
      return status.statusCode;
    }

  sprintf(refstamp, "1994-07-06 23:59:50 UTC Wed");

  if (lalDebugLevel > 0)
    {
      fprintf(stderr, "refstamp  = %s\n", refstamp);
      fprintf(stderr, "timestamp = %s\n", timestamp->data);
    }

  if (strcmp(refstamp, timestamp->data) != 0)
    {
      LALInfo(&status, "GPStoUTC conversion failed: wrong UTC result");
      LALCHARDestroyVector(&status, &timestamp);
      if (status.statusCode && lalDebugLevel > 0)
        {
          fprintf(stderr,
                  "TestGPStoUTC: LALCHARDestroyVector() failed, line %i, %s\n",
                  __LINE__, LALTESTGPSTOUTCC);
          REPORTSTATUS(&status);
          return status.statusCode;
        }
      REPORTSTATUS(&status);
      LALCheckMemoryLeaks();
      return 1;
    }

  /*
   * GPS 599184012 == 1998-12-31 23:59:60 UTC Thu (leap second introduced)
   */
  gpsTime.gpsSeconds = 599184012;
  gpsTime.gpsNanoSeconds = 0;

  if(!XLALGPSToUTC(&utcDate.unixDate, gpsTime.gpsSeconds))
    {
      fprintf(stderr, "TestGPStoUTC: XLALGPSToUTC() failed, line %i, %s\n",
              __LINE__, LALTESTGPSTOUTCC);
      return 1;
    }

  LALDateString(&status, timestamp, &utcDate);
  if (status.statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr, "TestGPStoUTC: LALDateString() failed, line %i, %s\n",
              __LINE__, LALTESTGPSTOUTCC);
      REPORTSTATUS(&status);
      return status.statusCode;
    }

  sprintf(refstamp, "1998-12-31 23:59:60 UTC Thu");

  if (lalDebugLevel > 0)
    {
      fprintf(stderr, "refstamp  = %s\n", refstamp);
      fprintf(stderr, "timestamp = %s\n", timestamp->data);
    }

  if (strcmp(refstamp, timestamp->data) != 0)
    {
      LALInfo(&status, "GPStoUTC conversion failed: wrong UTC result");
      LALCHARDestroyVector(&status, &timestamp);
      if (status.statusCode && lalDebugLevel > 0)
        {
          fprintf(stderr,
                  "TestGPStoUTC: LALCHARDestroyVector() failed, line %i, %s\n",
                  __LINE__, LALTESTGPSTOUTCC);
          REPORTSTATUS(&status);
          return status.statusCode;
        }
      REPORTSTATUS(&status);
      LALCheckMemoryLeaks();
      return 1;
    }

  /*
   * UPDATEME
   * GPS 835747214 == 2006-07-01 00:00:00 (one second past expiry)
   * Expect to fail with status code 5
   * (see the static constant maxtestedGPS in src/GPStoUTC.c)
   */
  gpsTime.gpsSeconds     = 835747214; /* use maxtestedGPS + 1 */
  gpsTime.gpsNanoSeconds = 0;

  if(!XLALGPSToUTC(&utcDate.unixDate, gpsTime.gpsSeconds))
    {
      fprintf(stderr, "TestGPStoUTC: XLALGPSToUTC() failed, line %i, %s\n",
              __LINE__, LALTESTGPSTOUTCC);
      return 1;
    }

  LALDateString(&status, timestamp, &utcDate);
  if (status.statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr, "TestGPStoUTC: LALDateString() failed, line %i, %s\n",
              __LINE__, LALTESTGPSTOUTCC);
      REPORTSTATUS(&status);
      LALCHARDestroyVector(&status, &timestamp);
      LALCheckMemoryLeaks();
      return status.statusCode;
    }

  /* UPDATEME */
  /* the date here should be one second after maxtestedGPS */
  sprintf(refstamp, "2006-07-01 00:00:00 UTC Sat");

  if (lalDebugLevel > 0)
    {
      fprintf(stderr, "refstamp  = %s\n", refstamp);
      fprintf(stderr, "timestamp = %s\n", timestamp->data);
    }

  if (strcmp(refstamp, timestamp->data) != 0)
    {
      LALInfo(&status, "GPStoUTC conversion failed: wrong UTC result");
      LALCHARDestroyVector(&status, &timestamp);
      if (status.statusCode && lalDebugLevel > 0)
        {
          fprintf(stderr,
                  "TestGPStoUTC: LALCHARDestroyVector() failed, line %i, %s\n",
                  __LINE__, LALTESTGPSTOUTCC);
          REPORTSTATUS(&status);
          return status.statusCode;
        }
      REPORTSTATUS(&status);
      LALCheckMemoryLeaks();
      return 1;
    }


  /*
   * GPS -44000
   * should fail.  No leap seconds, yet.
   */
  gpsTime.gpsSeconds     = -44000;
  gpsTime.gpsNanoSeconds = 0;

  if(!XLALGPSToUTC(&utcDate.unixDate, gpsTime.gpsSeconds))
    {
    if (XLALGetBaseErrno() != XLAL_EDOM) /* not expected error */
        {
          fprintf(stderr, "TestGPStoUTC: XLALGPSToUTC() failed, line %i, %s\n",
                  __LINE__, LALTESTGPSTOUTCC);
          return 1;
        }
    }
  else /* no error */
    {
      fprintf(stderr, "TestGPStoUTC: XLALGPSToUTC() failed, line %i, %s\n",
              __LINE__, LALTESTGPSTOUTCC);
      return 1;
    }
  XLALClearErrno();

  /*
   * Now, let's try converting GPS to UTC and back
   */
  gpsTime.gpsSeconds     = 701654354;
  gpsTime.gpsNanoSeconds = 0;

  if(!XLALGPSToUTC(&utcDate.unixDate, gpsTime.gpsSeconds))
    {
      fprintf(stderr, "TestGPStoUTC: XLALGPSToUTC() failed, line %i, %s\n",
              __LINE__, LALTESTGPSTOUTCC);
      return 1;
    }

  LALUTCtoGPS(&status, &tmpGps, &utcDate, &accuracy);
  if (status.statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr,
              "TestUTCtoGPS: error in LALUTCtoGPS, line %i, %s\n",
              __LINE__, LALTESTGPSTOUTCC);
      REPORTSTATUS(&status);
      return status.statusCode;
    }

  if (lalDebugLevel > 0)
    {
      fprintf(stderr, "\tgpsTime = {%d, %d}\n", gpsTime.gpsSeconds,
              gpsTime.gpsNanoSeconds);
      fprintf(stderr, "\ttmpGps  = {%d, %d}\n", tmpGps.gpsSeconds,
              tmpGps.gpsNanoSeconds);
    }


  if (tmpGps.gpsSeconds     != gpsTime.gpsSeconds ||
      tmpGps.gpsNanoSeconds != gpsTime.gpsNanoSeconds)
    {
      fprintf(stderr,
              "TestGPStoUTC: conversion from GPS to UTC and back to GPS failed, line %i, %s\n", __LINE__, LALTESTGPSTOUTCC);
      REPORTSTATUS(&status);
      return 1;
    }



  /*
   * Cleanup and exit
   */
  LALCHARDestroyVector(&status, &timestamp);
  if (status.statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr,
              "TestGPStoUTC: LALCHARDestroyVector() failed, line %i, %s\n",
              __LINE__, LALTESTGPSTOUTCC);
      REPORTSTATUS(&status);
      return status.statusCode;
    }
  if (lalDebugLevel > 0)
    REPORTSTATUS(&status);
  LALCheckMemoryLeaks();
  return 0;
}
