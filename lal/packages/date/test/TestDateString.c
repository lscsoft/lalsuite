#include <stdio.h>
#include <math.h>
#include <stdlib.h>

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
    char             refstamp[128];
    CHARVector      *timestamp = NULL;

    printf("TestDateString: WARNING: This test does nothing, yet\n");

   /*

    LALCHARCreateVector(&stat, &timestamp, (UINT4)64);

    unixtime.unixSeconds = 60858;
    unixtime.unixNanoSeconds = 0;

    LALUtime(&stat, &date, &unixtime);
    LALDateString(&stat, timestamp, &date);

    sprintf(refstamp, "1970-01-01 16:54:18 UTC Thu");

    if (lalDebugLevel > 0) {
      printf("refstamp  = %s\n", refstamp);
      printf("timestamp = %s\n", timestamp->data);
    }

    if (strcmp(refstamp, timestamp->data) == 0) {
      return 0;
    } else {
      return 1;
    }

    */

    return 0;
}
