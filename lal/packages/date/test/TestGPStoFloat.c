#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include <lal/LALStdlib.h>
#include <lal/Date.h>

INT4 lalDebugLevel = 0;

NRCSID (LALTESTGPSTOFLOATC, "$Id$");

int main(int argc, char *argv[])
{
  static LALStatus   status;
  LIGOTimeGPS        gpsTime = {0, 0};
  REAL8              realTime = 0.;
  LALTimeInterval    interval = {0, 0};
  REAL8              deltaT = 0.;

  if (argc > 1)
    lalDebugLevel = atoi(argv[1]);

  /* 1 */
  gpsTime.gpsSeconds     = 987654321;
  gpsTime.gpsNanoSeconds = 123456789;
  LALGPStoFloat(&status, &realTime, &gpsTime);

  if (status.statusCode && lalDebugLevel)
    {
      fprintf(stderr, "TestGPStoFloat: LALGPStoFloat() failed; line %i, %s\n",
              __LINE__, LALTESTGPSTOFLOATC);
      REPORTSTATUS(&status);
      return status.statusCode;
    }

    if (lalDebugLevel)
    {
      printf("TestGPStoFloat: expected %22.9f\n                got      %22.9f\n",
             987654321.123456789, realTime);
    }

  if (realTime != 987654321.123456789)
    {
      fprintf(stderr, "TestGPStoFloat: LALGPStoFloat() returned wrong value; expected %g, got %g\n",
              987654321.123456789, realTime);
      return 1;
    }
      

  /* 2 */
  realTime = 54321.123456789;
  LALFloatToGPS(&status, &gpsTime, &realTime);

  if (status.statusCode && lalDebugLevel)
    {
      fprintf(stderr, "TestGPStoFloat: LALFloatToGPS() failed; line %i, %s\n",
              __LINE__, LALTESTGPSTOFLOATC);
      REPORTSTATUS(&status);
      return status.statusCode;
    }

  if (lalDebugLevel)
    {
      printf("TestFloatToGPS: expected (%d, %d)\n                got      (%d, %d)\n",
             54321, 123456789, gpsTime.gpsSeconds, gpsTime.gpsNanoSeconds);
    }

  if (gpsTime.gpsSeconds != 54321 || gpsTime.gpsNanoSeconds != 123456789)
    {
      fprintf(stderr, "TestGPStoFloat: LALFloatToGPS() returned wrong value; expected (%d, %d), got (%d, %d)\n",
	      54321, 123456789, gpsTime.gpsSeconds, gpsTime.gpsNanoSeconds);
      return 2;
    }

  /* 3 */
  interval.seconds     = 123;
  interval.nanoSeconds = 654;

  LALIntervalToFloat(&status, &deltaT, &interval);
  
  if (status.statusCode && lalDebugLevel)
    {
      fprintf(stderr, "TestGPStoFloat: LALIntervalToFloat() failed; line %i, %s\n",
              __LINE__, LALTESTGPSTOFLOATC);
      REPORTSTATUS(&status);
      return status.statusCode;
    }

  if (lalDebugLevel)
    {
      printf("TestGPStoFloat: expected (%18.9f)\n", 123.000000654);
      printf("                     got (%18.9f)\n", deltaT);
    }

  if (deltaT != 123.000000654)
    {
      fprintf(stderr, "TestGPStoFloat: LALIntervalToFloat() returned wrong value; ");
      fprintf(stderr, "expected %18.9f, got %18.9f\n", 123.000000654, deltaT);
      return 3;
    }

  /* 4 */
  LALFloatToInterval(&status, &interval, &deltaT);

  if (status.statusCode && lalDebugLevel)
    {
      fprintf(stderr, "TestGPStoFloat: LALIntervalToFloat() failed; line %i, %s\n",
              __LINE__, LALTESTGPSTOFLOATC);
      REPORTSTATUS(&status);
      return status.statusCode;
    }

  if (lalDebugLevel)
    {
      printf("TestGPStoFloat: expected (%9d:%09d)\n", 123, 654);
      printf("                     got (%9d:%09d)\n", interval.seconds, interval.nanoSeconds);
    }

  if (interval.seconds != 123 || interval.nanoSeconds != 654)
    {
      fprintf(stderr, "TestGPStoFloat: LALFloatToInterval() returned wrong value; ");
      fprintf(stderr, "expected (%9d:%09d), got (%9d:%09d)\n", 123, 654, interval.seconds,
              interval.nanoSeconds);
      return 4;
    }

  return 0;
}
