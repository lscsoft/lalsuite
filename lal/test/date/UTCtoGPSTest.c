/*
*  Copyright (C) 2016 Karl Wette
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


static int test(struct tm *t, int correct_gps)
{

  {
    struct tm tf = *t;
    tf.tm_wday = tf.tm_yday = -1;
    XLAL_CHECK(XLALFillUTC(&tf) != NULL, XLAL_EFUNC, "UTCtoGPSTest: error in XLALFillUTC()");
    XLAL_CHECK(tf.tm_wday == t->tm_wday, XLAL_EFAILED, "UTCtoGPSTest: incorrect day of week\n  output   = %d\n  expected = %d\n", tf.tm_wday, t->tm_wday);
    XLAL_CHECK(tf.tm_yday == t->tm_yday, XLAL_EFAILED, "UTCtoGPSTest: incorrect day of year\n  output   = %d\n  expected = %d\n", tf.tm_yday, t->tm_yday);
  }

  t->tm_wday = t->tm_yday = -1;
  int gps = XLALUTCToGPS(t);
  XLAL_CHECK(xlalErrno == 0, XLAL_EFUNC, "UTCtoGPSTest: error in XLALUTCToGPS()");

  if (gps != correct_gps)
    {
      char buf[64];
      strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", t);
      XLAL_ERROR(XLAL_EFAILED, "UTCtoGPSTest: incorrect UTC to GPS conversion\n  input    = %s\n  output   = %d\n  expected = %d\n", buf, gps, correct_gps);
    }

  return XLAL_SUCCESS;

}

int main(void)
{

  struct tm utcDate;
  time_t sec;

  utcDate.tm_year = 80;
  utcDate.tm_yday =  5;
  utcDate.tm_wday =  0;
  utcDate.tm_mon  =  0;
  utcDate.tm_mday =  6;
  utcDate.tm_hour =  0;
  utcDate.tm_min  =  0;
  utcDate.tm_sec  =  0;
  utcDate.tm_isdst = 0;
  XLAL_CHECK_MAIN(test(&utcDate, 0) == XLAL_SUCCESS, XLAL_EFAILED);

  utcDate.tm_year =  94;
  utcDate.tm_yday = 186;
  utcDate.tm_wday =   3;
  utcDate.tm_mon  =   6;
  utcDate.tm_mday =   6;
  utcDate.tm_hour =  23;
  utcDate.tm_min  =  59;
  utcDate.tm_sec  =  50;
  utcDate.tm_isdst =  1;
  XLAL_CHECK_MAIN(test(&utcDate, 457574400) == XLAL_SUCCESS, XLAL_EFAILED);

  utcDate.tm_year =  94;
  utcDate.tm_yday = 181;
  utcDate.tm_wday =   5;
  utcDate.tm_mon  =   6;
  utcDate.tm_mday =   1;
  utcDate.tm_hour =   0;
  utcDate.tm_min  =   0;
  utcDate.tm_sec  =   0;
  utcDate.tm_isdst =  1;
  XLAL_CHECK_MAIN(test(&utcDate, 457056010) == XLAL_SUCCESS, XLAL_EFAILED);

  for (sec = 457056007; sec < 457056012; sec++)
    {
      utcDate.tm_year =  94;
      utcDate.tm_sec  =  58 + (sec - 457056007);
      if (utcDate.tm_sec <= 60)
        {
          utcDate.tm_yday = 180;
          utcDate.tm_wday =   4;
          utcDate.tm_mon  =   5;
          utcDate.tm_mday =  30;
          utcDate.tm_hour =  23;
          utcDate.tm_min  =  59;
        }
      else
        {
          utcDate.tm_sec -=  61;
          utcDate.tm_yday = 181;
          utcDate.tm_wday =   5;
          utcDate.tm_mon  =   6;
          utcDate.tm_mday =   1;
          utcDate.tm_hour =  00;
          utcDate.tm_min  =  00;
        }
      utcDate.tm_isdst =  1;
      XLAL_CHECK_MAIN(test(&utcDate, sec) == XLAL_SUCCESS, XLAL_EFAILED);
    }

  utcDate.tm_year =  94;
  utcDate.tm_yday = 319;
  utcDate.tm_wday =   3;
  utcDate.tm_mon  =  10;
  utcDate.tm_mday =  16;
  utcDate.tm_hour =   0;
  utcDate.tm_min  =   0;
  utcDate.tm_sec  =   0;
  utcDate.tm_isdst =  0;
  XLAL_CHECK_MAIN(test(&utcDate, 468979210) == XLAL_SUCCESS, XLAL_EFAILED);

  LALCheckMemoryLeaks();

  return EXIT_SUCCESS;

}
