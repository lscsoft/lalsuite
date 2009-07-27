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

NRCSID (LALTESTUTCTOGPSC, "$Id$");

static int test(struct tm *t, int correct_gps, int line)
{
  int gps = XLALUTCToGPS(t);
 
  if (XLALGetBaseErrno())
    {
      fprintf(stderr, "TestUTCtoGPS: error in XLALUTCToGPS(), line %i\n", line);
      return -1;
    }

  if (lalDebugLevel > 0)
    {
      fprintf(stderr, "Input = %s\tOutput =   %d\n\tExpected = %d\n",
      asctime(t), gps, correct_gps);
    }

  if (gps != correct_gps)
    {
      if (lalDebugLevel > 0)
        fprintf(stderr, "TestUTCtoGPS: error, line %i\n", line);
      return -1;
    }

  return 0;
}

#define TEST(t, correct_gps) test(t, correct_gps, __LINE__)


int main(int argc, char *argv[])
{
  struct tm utcDate;
  time_t sec;

  if (argc > 1)
    lalDebugLevel = atoi(argv[1]);

  utcDate.tm_year = 80;
  utcDate.tm_yday =  5;
  utcDate.tm_wday = -1;	/* unused */
  utcDate.tm_mon  =  0;	/* unused */
  utcDate.tm_mday =  6;	/* unused */
  utcDate.tm_hour =  0;
  utcDate.tm_min  =  0;
  utcDate.tm_sec  =  0;
  utcDate.tm_isdst = 0;
  if (TEST(&utcDate, 0))
    return 1;

  utcDate.tm_year =  94;
  utcDate.tm_yday = 186;
  utcDate.tm_wday =  -1;	/* unused */
  utcDate.tm_mon  =   6;	/* unused */
  utcDate.tm_mday =   6;	/* unused */
  utcDate.tm_hour =  23;
  utcDate.tm_min  =  59;
  utcDate.tm_sec  =  50;
  utcDate.tm_isdst =  1;
  if (TEST(&utcDate, 457574400))
    return 1;

  utcDate.tm_year =  94;
  utcDate.tm_yday = 181;
  utcDate.tm_wday =  -1;	/* unused */
  utcDate.tm_mon  =   6;	/* unused */
  utcDate.tm_mday =   1;	/* unused */
  utcDate.tm_hour =   0;
  utcDate.tm_min  =   0;
  utcDate.tm_sec  =   0;
  utcDate.tm_isdst =  1;
  if (TEST(&utcDate, 457056010))
    return 1;

  for (sec = 457056007; sec < 457056012; sec++)
    {
      utcDate.tm_year =  94;
      utcDate.tm_wday =  -1;	/* unused */
      utcDate.tm_sec  =  58 + (sec - 457056007);
      if (utcDate.tm_sec <= 60)
        {
          utcDate.tm_yday = 180;
          utcDate.tm_mon  =   5;	/* unused */
          utcDate.tm_mday =  30;	/* unused */
          utcDate.tm_hour =  23;
          utcDate.tm_min  =  59;
        }
      else
        {
          utcDate.tm_sec -=  61;
          utcDate.tm_yday = 181;
          utcDate.tm_mon  =   6;	/* unused */
          utcDate.tm_mday =   1;	/* unused */
          utcDate.tm_hour =  00;
          utcDate.tm_min  =  00;
        }
      utcDate.tm_isdst =  1;
      if (TEST(&utcDate, sec))
        return 1;
    }

  utcDate.tm_year =  94;
  utcDate.tm_yday = 319;
  utcDate.tm_wday =  -1;	/* unused */
  utcDate.tm_mon  =  10;	/* unused */
  utcDate.tm_mday =  16;	/* unused */
  utcDate.tm_hour =   0;
  utcDate.tm_min  =   0;
  utcDate.tm_sec  =   0;
  utcDate.tm_isdst =  0;
  if (TEST(&utcDate, 468979210))
    return 1;

  LALCheckMemoryLeaks();
  return 0;
}

