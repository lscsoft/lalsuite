#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>

INT4 lalDebugLevel = LALALLDBG;

NRCSID (LALTESTGPSTOUTCC, "$Id$");

int main(void)
{
  static LALStatus    status;
  LIGOTimeGPS         gpsTime = {0, 0};
  LALDate             utcDate;
  LALLeapSecAccuracy  accuracy = LALLEAPSEC_LOOSE;
  CHARVector         *timestamp = NULL;
  char                refstamp[128];

  LALCHARCreateVector(&status, &timestamp, (UINT4)128);
  REPORTSTATUS(&status);

  LALGPStoUTC(&status, &utcDate, &gpsTime, &accuracy);
  if (status.statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr, "TestGPStoUTC: LALGPStoUTC() failed, line %i, %s\n",
              __LINE__, LALTESTGPSTOUTCC);
      REPORTSTATUS(&status);
      LALCHARDestroyVector(&status, &timestamp);
      return status.statusCode;
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


  sprintf(refstamp, "1980-01-06 00:00:00 UTC Sun");
  fprintf(stderr, "refstamp  = %s\n", refstamp);
  fprintf(stderr, "timestamp = %s\n", timestamp->data);

  if (strcmp(refstamp, timestamp->data) != 0)
    {
      fprintf(stderr, "TestGPStoUTC: date strings do not match, line %i, %s\n",
              __LINE__, LALTESTGPSTOUTCC);
      LALCHARDestroyVector(&status, &timestamp);
      REPORTSTATUS(&status);
      LALCheckMemoryLeaks();
      return 1;
    }

  /* check another time */
  gpsTime.gpsSeconds = 457574400;
  gpsTime.gpsNanoSeconds = 0;
  
  LALGPStoUTC(&status, &utcDate, &gpsTime, &accuracy);
  if (status.statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr, "TestGPStoUTC: LALGPStoUTC() failed, line %i, %s\n",
              __LINE__, LALTESTGPSTOUTCC);
      REPORTSTATUS(&status);
      LALCHARDestroyVector(&status, &timestamp);
      return status.statusCode;
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

  sprintf(refstamp, "1994-07-06 23:59:50 UTC Wed");
  fprintf(stderr, "refstamp  = %s\n", refstamp);
  fprintf(stderr, "timestamp = %s\n", timestamp->data);

  if (strcmp(refstamp, timestamp->data) != 0)
    {
      LALCHARDestroyVector(&status, &timestamp);
      REPORTSTATUS(&status);
      LALCheckMemoryLeaks();
      return 1;
    }

  LALCHARDestroyVector(&status, &timestamp);
  REPORTSTATUS(&status);
  LALCheckMemoryLeaks();
  return 0;
}
