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
  LALMSTUnitsAndAcc mstUnitsAndAcc;
  LIGOTimeGPS      gps = {0., 0.};
  REAL8            gmst;

  gps.gpsSeconds = 61094;
  mstUnitsAndAcc.units = MST_RAD;
  mstUnitsAndAcc.accuracy = LALLEAPSEC_LOOSE;

  for (gps.gpsNanoSeconds =99999999; gps.gpsNanoSeconds < 1000000000;
       gps.gpsNanoSeconds+=10000000)
    {
      LALGPStoGMST1(&status, &gmst, &gps, &mstUnitsAndAcc);
      printf("nSec = %d\tgmst = %g\n", gps.gpsNanoSeconds, gmst);
    }

  return 0;
}
  
