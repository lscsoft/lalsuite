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

  if (argc > 1)
    lalDebugLevel = atoi(argv[1]);

  gpsTime.gpsSeconds     = 987654321;
  gpsTime.gpsNanoSeconds = 123456789;
  LALGPStoFloat(&status, &realTime, &gpsTime);

  if (status.statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr, "TestGPStoFloat: LALGPStoFloat() failed; line %i, %s\n",
              __LINE__, LALTESTGPSTOFLOATC);
      REPORTSTATUS(&status);
      return status.statusCode;
    }

    if (lalDebugLevel > 0)
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
      

  realTime = 987654321.123456;
  LALFloatToGPS(&status, &gpsTime, &realTime);

  if (status.statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr, "TestGPStoFloat: LALFloatToGPS() failed; line %i, %s\n",
              __LINE__, LALTESTGPSTOFLOATC);
      REPORTSTATUS(&status);
      return status.statusCode;
    }

  if (lalDebugLevel > 0)
    {
      printf("TestGPStoFloat: expected (%d, %d)\n                got      (%d, %d)\n",
             987654321, 123456000, gpsTime.gpsSeconds,
             gpsTime.gpsNanoSeconds);
    }

  if (gpsTime.gpsSeconds != 987654321 || gpsTime.gpsNanoSeconds != 123456000)
    {
      fprintf(stderr, "TestGPStoFloat: LALFloatToGPS() returned wrong value; expected (%d, %d), got (%d, %d)\n",
              987654321, 123456000, gpsTime.gpsSeconds,
              gpsTime.gpsNanoSeconds);
      return 2;
    }

  return 0;
}
