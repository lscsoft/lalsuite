/*
*  Copyright (C) 2007 David Chin, Jolien Creighton, Kipp Cannon, Reinhard Prix
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

/*
 * Author: David Chin <dwchin@umich.edu> +1-734-730-1274
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>


#define SUCCESS 0
#define FAILURE 1

#define TRUE  1
#define FALSE 0

#define LONGESTSTR 256

/* Error codes and messages */
#define LALTESTINCREMENTGPSC_ENORM 	0
#define LALTESTINCREMENTGPSC_ESUB  	1

#define LALTESTINCREMENTGPSC_MSGENORM "Normal exit"
#define LALTESTINCREMENTGPSC_MSGESUB  "LAL Subroutine failed"

/* module-scope variable */
BOOLEAN verbose_p = FALSE;
const INT4 oneBillion = 1000000000;



static int print_compare_errmsg_maybe(const char *msg,
                                      int expected_val,
                                      int got_val);

static BOOLEAN compare_gps_ok(const LIGOTimeGPS *p_gps1,
                              const LIGOTimeGPS *p_gps2,
                              int expected_val);

/* static int sprint_time_interval(char *str,
                                const LALTimeInterval *p_time_interval); */

int main(void)
{
  LIGOTimeGPS      gps;


  if (lalDebugLevel >= 4)
    verbose_p = TRUE;

  /*
   * TEST No. 0 -- test LALCompareGPS()
   */
  {
    LIGOTimeGPS gps2;

    gps.gpsSeconds      =   1249389;
    gps.gpsNanoSeconds  = 498512352;

    /* equal */
    gps2.gpsSeconds     =   1249389;
    gps2.gpsNanoSeconds = 498512352;

    if (!compare_gps_ok(&gps, &gps2, 0))
      return FAILURE;

    /* later */
    gps2.gpsNanoSeconds -= 1;

    if (!compare_gps_ok(&gps, &gps2, 1))
      return FAILURE;

    /* earlier */
    gps2.gpsNanoSeconds += 2;

    if (!compare_gps_ok(&gps, &gps2, -1))
      return FAILURE;

    /* later */
    gps2.gpsSeconds -= 1;
    gps2.gpsNanoSeconds -= 1;

    if (!compare_gps_ok(&gps, &gps2, 1))
      return FAILURE;

    /* earlier */
    gps2.gpsSeconds += 2;

    if (!compare_gps_ok(&gps, &gps2, -1))
      return FAILURE;

  } /* END: test of LALCompareGPS() */


  /*----------------------------------------------------------------------*/
  /* test XLALGPSAdd() and XLALGPSDiff() */
  {
    LIGOTimeGPS gps1, gps2, gpsTest, gps3;
    double deltaT1, deltaT0;

    gps1.gpsSeconds = 714153733;
    gps1.gpsNanoSeconds = 23421234;

    gps2.gpsSeconds = 715615153;
    gps2.gpsNanoSeconds = 712343412;
    /* now use XLALGPSDiff() to calculate the difference in seconds */
    deltaT1 = XLALGPSDiff(&gps1, &gps2);
    /* compare to correct result (error must be < 0.5 ns) */
    deltaT0 = -1461420.688922178;
    if (fabs(deltaT1 - deltaT0) >= .5e-9) {
      LALPrintError ("Failure in XLALGPSDiff(): got %.17g instead of %.17g\n", deltaT1, deltaT0);
      return FAILURE;
    }
    /* now see if deltaT1 takes us from gps1 to gps2 and vice-versa */
    gpsTest = gps2;
    XLALGPSAdd(&gpsTest, deltaT1);
    if ( (gpsTest.gpsSeconds != gps1.gpsSeconds) || (gpsTest.gpsNanoSeconds != gps1.gpsNanoSeconds) ) {
      LALPrintError ("Failure in 1.) XLALGPSAdd(): got %d.%09ds instead of %d.%09ds\n",
		     gpsTest.gpsSeconds, gpsTest.gpsNanoSeconds, gps1.gpsSeconds, gps1.gpsNanoSeconds);
      return FAILURE;
    }
    /* no go the other direction..*/
    deltaT1 = -deltaT1;
    gpsTest = gps1;
    XLALGPSAdd(&gpsTest, deltaT1);
    if ( (gpsTest.gpsSeconds != gps2.gpsSeconds) || (gpsTest.gpsNanoSeconds != gps2.gpsNanoSeconds) ) {
      LALPrintError ("Failure in 2.) XLALGPSAdd(): got %d.%09ds instead of %d.%09ds\n",
		     gpsTest.gpsSeconds, gpsTest.gpsNanoSeconds, gps2.gpsSeconds, gps2.gpsNanoSeconds);
      return FAILURE;
    }
    /* test over-run in ns-position is handled properly */
    gpsTest = gps2;
    XLALGPSAdd(&gpsTest, deltaT1);
    gps3.gpsSeconds = 717076574;
    gps3.gpsNanoSeconds = 401265590;
    if ( (gpsTest.gpsSeconds != gps3.gpsSeconds) || (gpsTest.gpsNanoSeconds != gps3.gpsNanoSeconds) ) {
      LALPrintError ("Failure in 3.) XLALGPSAdd(): got %d.%09ds instead of %d.%09ds\n",
		     gpsTest.gpsSeconds, gpsTest.gpsNanoSeconds, gps2.gpsSeconds, gps2.gpsNanoSeconds);
      return FAILURE;
    }
    gpsTest = gps1;
    XLALGPSAdd(&gpsTest, deltaT1);
    gps3.gpsSeconds = 715615153;
    gps3.gpsNanoSeconds = 712343412;
    if ( (gpsTest.gpsSeconds != gps3.gpsSeconds) || (gpsTest.gpsNanoSeconds != gps3.gpsNanoSeconds) ) {
      LALPrintError ("Failure in 4.) XLALGPSAdd(): got %d.%09ds instead of %d.%09ds\n",
		     gpsTest.gpsSeconds, gpsTest.gpsNanoSeconds, gps2.gpsSeconds, gps2.gpsNanoSeconds);
      return FAILURE;
    }


  } /* testing XLALGPSAdd() and XLALGPSDiff() */
  /*----------------------------------------------------------------------*/

  return SUCCESS;
} /* END: main() */



static int print_compare_errmsg_maybe(const char *msg,
                                      int expected_val,
                                      int got_val)
{
  if (verbose_p)
    fprintf(stderr, "%s; expected %d got %d\n\n", msg, expected_val, got_val);
  return 0;
}



static BOOLEAN compare_gps_ok(const LIGOTimeGPS *p_gps1,
                              const LIGOTimeGPS *p_gps2,
                              int  expected_val)
{
  if (XLALGPSCmp(p_gps1, p_gps2) != expected_val)
    {
      print_compare_errmsg_maybe("XLALGPSCmp() failed",
                                 expected_val, XLALGPSCmp(p_gps1, p_gps2));
      return FALSE;
    }
  else
    {
      return TRUE;
    }
}


#if 0 /* NOT USED */
static int sprint_time_interval(char *str,
                                const LALTimeInterval *p_time_interval)
{
  return snprintf(str, LONGESTSTR - 1 - strlen(str),
                  "%9d:%09d", p_time_interval->seconds,
                  p_time_interval->nanoSeconds);
}
#endif
