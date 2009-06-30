/*
*  Copyright (C) 2007 Bernd Machenschalk, Philip Charlton, Jolien Creighton
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


/*Includes for UserTimeVal()*/
#include <limits.h>
#include <sys/times.h>

/*includes for WallTimeVal()*/
#include <time.h>
#include <sys/time.h>

#include <lal/LALRCSID.h>
NRCSID (TIMINGFUNCTIONSC,"$Id$");

double UserTimeVal(void);
double WallTimeVal(void);

double UserTimeVal(void) {
  /*returns the user time difference in seconds between calls.
    It should be used like this:
      UserTimeVal();
      FunctionToBeTimed(ARGS);
      TempTime = UserTime();
      printf("FunctionToBeTimed took %g user seconds\n");
    Please note that this function can only deal with one wrap around.
    Wrap arounds occur about every 35 user minutes.*/

  static clock_t last_usertime = 0;
  struct tms now;
  double val;
  clock_t Check;


  Check = times(&now);

  if (Check < 0) {
    printf("There was a problem. Check = %ld\n", Check);
  }

  if (now.tms_utime >= last_usertime) {

    val = (now.tms_utime - last_usertime) / (double)CLOCKS_PER_SEC;

  } else if (now.tms_utime < last_usertime) {

   printf("Wrapped around. last_usertime = %ld, now.tms_utime = %ld, diff = %ld\n", \
	   last_usertime, now.tms_utime, now.tms_utime + (LONG_MAX - last_usertime));

    val = (now.tms_utime + (LONG_MAX - last_usertime))/(double)CLOCKS_PER_SEC;

  }
  if (val < 0){
    printf("now.tms_utime = %ld, last_usertime = %ld, val = %g\n", \
	   last_usertime, now.tms_utime, val);
  }

  last_usertime = now.tms_utime;

  return(val);
}

double WallTimeVal(void) {
  /*returns the user time difference in seconds between calls.
    It should be used like this:
      WallTimeVal();
      FunctionToBeTimed(ARGS);
      TempTime = WallTime();
      printf("FunctionToBeTimed took %g user seconds\n");
    Please note that this function can only deal with one wrap around.
    Wrap arounds occur about every 35 user minutes.*/

  static clock_t last_walltime = 0;
  static clock_t now_walltime = 0;
  struct tms now;
  double val;
  clock_t Check;


  Check = times(&now);

  if (Check < 0) {
    printf("There was a problem. Check = %ld\n", Check);
  }
  now_walltime = now.tms_utime + now.tms_stime;

  if (now_walltime >= last_walltime) {

    val = (now_walltime - last_walltime) / (double)CLOCKS_PER_SEC;

  } else if (now_walltime < last_walltime) {

    printf("Wrapped around. last_walltime = %ld, now_walltime = %ld, diff = %ld\n", \
           last_walltime, now_walltime, now_walltime + (LONG_MAX - last_walltime));

    val = (now_walltime + (LONG_MAX - last_walltime))/(double)CLOCKS_PER_SEC;

  }
  if (val < 0){
    printf("now_walltime = %ld, last_walltime = %ld, val = %g\n", \
           last_walltime, now_walltime, val);
  }

  last_walltime = now_walltime;

  return(val);
}

