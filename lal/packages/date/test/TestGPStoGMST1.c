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
#include <time.h>

#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>

INT4 lalDebugLevel = 0;

NRCSID (TESTGPSTOGMST1C, "$Id$");

int main(void)
{
  static LALStatus status;
  LIGOTimeGPS      gps = {0., 0.};
  REAL8            gmst;

  gps.gpsSeconds = 61094;

  for (gps.gpsNanoSeconds =99999999; gps.gpsNanoSeconds < 1000000000;
       gps.gpsNanoSeconds+=10000000)
    {
      gmst = XLALGreenwichMeanSiderealTime(&gps);
      printf("nSec = %d\tgmst = %g\n", gps.gpsNanoSeconds, gmst);
    }

  return 0;
}

