#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>

INT4 lalDebugLevel = LALALLDBG;

NRCSID (TESTDATESTRINGC, "$Id$");

int main(void)
{
    static LALStatus stat;
    LIGOTimeUnix     unixtime;
    time_t           tmpsecs;
    LALDate          date;
    CHAR             refstamp[128];
    CHAR             tmpstamp[128];
    CHARVector      *timestamp = NULL;
    struct tm        tmdate;

    LALCHARCreateVector(&stat, &timestamp, (UINT4)64);
    REPORTSTATUS (&stat);

    unixtime.unixSeconds = (24*365 + 8*366 + 2*31 + 28)*86400 - 1 + 24;
    unixtime.unixNanoSeconds = 0;

    LALUtime(&stat, &date, &unixtime);
    REPORTSTATUS(&stat);

    LALDateString(&stat, timestamp, &date);
    REPORTSTATUS(&stat);

    sprintf(refstamp, "2002-03-31 23:59:59 UTC Sun");
    fprintf(stderr, "refstamp  = %s\n", refstamp);
    fprintf(stderr, "timestamp = %s\n", timestamp->data);

    if (strcmp(refstamp, timestamp->data) == 0)
      {
        LALCHARDestroyVector(&stat, &timestamp);
        REPORTSTATUS(&stat);
        LALCheckMemoryLeaks();
        return 0;
      }
    else
      {
        LALCHARDestroyVector(&stat, &timestamp);
        REPORTSTATUS(&stat);
        LALCheckMemoryLeaks();
        return 1;
      }
}
