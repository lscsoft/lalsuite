#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>

INT4 lalDebugLevel = 2;

NRCSID (TESTDATESTRINGC, "$Id$");

int main(void)
{
    static LALStatus stat;
    LIGOTimeUnix     unixtime;
    LALDate          date;
    CHAR             refstamp[128];
    CHARVector      *timestamp = NULL;
    struct tm        tmdate;

    LALCHARCreateVector(&stat, &timestamp, (UINT4)64);

    /* unixtime.unixSeconds = 60858; */
    /* 24 leap seconds introduced since 1970 */
    unixtime.unixSeconds = (24*365 + 8*366 + 2*31 + 28)*86400 - 1 + 24;
    unixtime.unixNanoSeconds = 0;

    LALUtime(&stat, &date, &unixtime);

    LALDateString(&stat, timestamp, &date);

    sprintf(refstamp, "2002-03-31 23:59:59 UTC Sun");

    if (lalDebugLevel > 0)
    {
      fprintf(stderr, "asctime(gmtime(&(unixtime.unixSeconds))) = %s",
              asctime(gmtime(&(unixtime.unixSeconds))));
      fprintf(stderr, "                                    expect Mon Apr  1 00:00:23 2002\n");
      unixtime.unixSeconds = 0;
      fprintf(stderr, "asctime(gmtime(&(unixtime.unixSeconds))) = %s",
              asctime(gmtime(&(unixtime.unixSeconds))));
      fprintf(stderr, "                                    expect Thu Jan  1 00:00:00 1970\n");
      fprintf(stderr, "refstamp  = %s\n", refstamp);
      fprintf(stderr, "timestamp = %s\n", timestamp->data);
    }

    if (strcmp(refstamp, timestamp->data) == 0)
      {
        return 0;
      }
    else
      {
        return 1;
      }

    LALCheckMemoryLeaks();
    
    return 0;
}
