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
#ifndef _REENTRANT
#define _REENTRANT
#endif
#ifndef __USE_POSIX
#define __USE_POSIX
#endif
#include <time.h>

#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>

INT4 lalDebugLevel = 0;

NRCSID (LALTESTDATESTRINGC, "$Id$");

int main(int argc, char **argv)
{
    static LALStatus stat;
    time_t           tmpsecs;
    LALDate          date;
    CHAR             refstamp[128];
    CHARVector      *timestamp = NULL;

    if (argc > 1)
      lalDebugLevel = atoi(argv[1]);

    LALCHARCreateVector(&stat, &timestamp, (UINT4)64);

    if (stat.statusCode && lalDebugLevel > 0)
      {
        fprintf(stderr,
                "TestDateString: error in LALCHARCreateVector, line %i, %s\n",
                __LINE__, LALTESTDATESTRINGC);
        REPORTSTATUS(&stat);
        return stat.statusCode;
      }
    if (lalDebugLevel > 2)
      REPORTSTATUS(&stat);

    tmpsecs = (24*365 + 8*366 + 2*31 + 28)*86400 - 1;
    gmtime_r(&tmpsecs, &(date.unixDate));

    LALDateString(&stat, timestamp, &date);

    if (stat.statusCode && lalDebugLevel > 0)
      {
        fprintf(stderr,
                "TestDateString: LALDateString() failed, line %i, %s\n",
                __LINE__, LALTESTDATESTRINGC);
        REPORTSTATUS(&stat);
        return stat.statusCode;
      }

    if (!stat.statusCode)
      {
        sprintf(refstamp, "2002-03-31 23:59:59 UTC Sun");
        if (lalDebugLevel > 2)
          {
            fprintf(stderr, "refstamp  = %s\n", refstamp);
            fprintf(stderr, "timestamp = %s\n", timestamp->data);
          }
      }

    if (strcmp(refstamp, timestamp->data) == 0)
      {
        LALCHARDestroyVector(&stat, &timestamp);
        if (lalDebugLevel > 2)
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
