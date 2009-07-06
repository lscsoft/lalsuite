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
 * Revision: $Id$
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>

INT4 lalDebugLevel = 0;

NRCSID (LALTESTINCREMENTGPSC, "$Id$");


#define SUCCESS 0
#define FAILURE 1

#define TRUE  1
#define FALSE 0

#define LONGESTSTR 256

/* Error codes and messages */

/************** <lalErrTable file="LALTESTINCREMENTGPSCErrorTable"> */
#define LALTESTINCREMENTGPSC_ENORM 	0
#define LALTESTINCREMENTGPSC_ESUB  	1

#define LALTESTINCREMENTGPSC_MSGENORM "Normal exit"
#define LALTESTINCREMENTGPSC_MSGESUB  "LAL Subroutine failed"
/******************************************** </lalErrTable> */


/*********************************************************************/
/* Macros for printing errors & testing subroutines (from Creighton) */
/*********************************************************************/

#define ERROR( code, msg, statement )                                \
do {                                                                 \
  if ( lalDebugLevel & LALERROR )                                    \
    LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n" \
                   "        %s %s\n", (code), *argv, __FILE__,       \
              __LINE__, LALTESTINCREMENTGPSC, statement ? statement :  \
                   "", (msg) );                                      \
} while (0)

#define INFO( statement )                                            \
do {                                                                 \
  if ( lalDebugLevel & LALINFO )                                     \
    LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"     \
                   "        %s\n", *argv, __FILE__, __LINE__,        \
              LALTESTINCREMENTGPSC, (statement) );                     \
} while (0)

#define SUB( func, statusptr )                                       \
do {                                                                 \
  if ( (func), (statusptr)->statusCode ) {                           \
    ERROR( LALTESTINCREMENTGPSC_ESUB, LALTESTINCREMENTGPSC_MSGESUB,      \
           "Function call \"" #func "\" failed:" );                  \
    return LALTESTINCREMENTGPSC_ESUB;                                  \
  }                                                                  \
} while (0)
/******************************************************************/





/* module-scope variable */
BOOLEAN verbose_p = FALSE;
const INT4 oneBillion = 1000000000;



static int print_compare_errmsg_maybe(const char *msg,
                                      LALGPSCompareResult expected_val,
                                      LALGPSCompareResult got_val);

static BOOLEAN compare_gps_ok(LALStatus *status,
                              const LIGOTimeGPS *p_gps1,
                              const LIGOTimeGPS *p_gps2,
                              LALGPSCompareResult expected_val);

static int sprint_gps_time(char *str, const LIGOTimeGPS *p_gps);

/* static int sprint_time_interval(char *str,
                                const LALTimeInterval *p_time_interval); */

static int print_incr_errmsg_maybe(const char *msg,
                                   const LIGOTimeGPS *p_expected_gps,
                                   const LIGOTimeGPS *p_got_gps);

static BOOLEAN increment_gps_ok(LALStatus *status,
                                LALTimeInterval *p_delta,
                                LIGOTimeGPS *p_expected_gps,
                                const LIGOTimeGPS *p_got_gps);

static BOOLEAN decrement_gps_ok(LALStatus *status,
                                LALTimeInterval *p_delta,
                                LIGOTimeGPS *p_expected_gps,
                                const LIGOTimeGPS *p_got_gps);


int main(int argc, char **argv)
{
  static LALStatus status;
  LIGOTimeGPS      gps;
  LIGOTimeGPS      another_gps;
  LIGOTimeGPS      expected_gps;
  LIGOTimeGPS      decremented_gps;
  LALTimeInterval  deltaT;
  LALTimeInterval  expected_deltaT;


  if (argc > 1)
    lalDebugLevel = atoi(argv[1]);

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

    if (!compare_gps_ok(&status, &gps, &gps2, LALGPS_EQUAL))
      return FAILURE;

    /* later */
    gps2.gpsNanoSeconds -= 1;

    if (!compare_gps_ok(&status, &gps, &gps2, LALGPS_LATER))
      return FAILURE;

    /* earlier */
    gps2.gpsNanoSeconds += 2;

    if (!compare_gps_ok(&status, &gps, &gps2, LALGPS_EARLIER))
      return FAILURE;

    /* later */
    gps2.gpsSeconds -= 1;
    gps2.gpsNanoSeconds -= 1;

    if (!compare_gps_ok(&status, &gps, &gps2, LALGPS_LATER))
      return FAILURE;

    /* earlier */
    gps2.gpsSeconds += 2;

    if (!compare_gps_ok(&status, &gps, &gps2, LALGPS_EARLIER))
      return FAILURE;

  } /* END: test of LALCompareGPS() */


  /*
   * TEST No. 1 -- test of LALIncrementGPS()
   */
  gps.gpsSeconds     =   0;
  gps.gpsNanoSeconds =   0;
  deltaT.seconds     =  13;
  deltaT.nanoSeconds = 100;
  expected_gps.gpsSeconds =      13;
  expected_gps.gpsNanoSeconds = 100;

  if (!increment_gps_ok(&status, &deltaT, &gps, &expected_gps))
    return FAILURE;

  gps.gpsSeconds     = 472139;
  gps.gpsNanoSeconds =   1489;
  deltaT.seconds     =  10000;
  deltaT.nanoSeconds = 700000;
  expected_gps.gpsSeconds     = 482139;
  expected_gps.gpsNanoSeconds = 701489;

  if (!increment_gps_ok(&status, &deltaT, &gps, &expected_gps))
    return FAILURE;

  /*----------------------------------------------------------------------*/
  /* test LALAddFloatToGPS() and LALDeltaFloatGPS() */
  {
    LIGOTimeGPS gps1, gps2, gpsTest, gps3;
    REAL8 deltaT1, deltaT0;

    gps1.gpsSeconds = 714153733;
    gps1.gpsNanoSeconds = 23421234;

    gps2.gpsSeconds = 715615153;
    gps2.gpsNanoSeconds = 712343412;
    /* now use DeltaFloatGPS to calculate the difference in seconds */
    SUB ( LALDeltaFloatGPS (&status, &deltaT1, &gps1, &gps2), &status);
    /* compare to correct result */
    deltaT0 = -1461420.688922178;
    if (deltaT1 != deltaT0) {
      LALPrintError ("Failure in LALDeltaFloatGPS(): got %20.9f instead of %20.9f\n", deltaT1, deltaT0);
      return FAILURE;
    }
    /* now see if deltaT1 takes us from gps1 to gps2 and vice-versa */
    SUB (LALAddFloatToGPS (&status, &gpsTest, &gps2, deltaT1), &status);
    if ( (gpsTest.gpsSeconds != gps1.gpsSeconds) || (gpsTest.gpsNanoSeconds != gps1.gpsNanoSeconds) ) {
      LALPrintError ("Failure in 1.) LALAddFloatToGPS(): got %d.%09ds instead of %d.%09ds\n",
		     gpsTest.gpsSeconds, gpsTest.gpsNanoSeconds, gps1.gpsSeconds, gps1.gpsNanoSeconds);
      return FAILURE;
    }
    /* no go the other direction..*/
    deltaT1 = -deltaT1;
    SUB (LALAddFloatToGPS (&status, &gpsTest, &gps1, deltaT1), &status);
    if ( (gpsTest.gpsSeconds != gps2.gpsSeconds) || (gpsTest.gpsNanoSeconds != gps2.gpsNanoSeconds) ) {
      LALPrintError ("Failure in 2.) LALAddFloatToGPS(): got %d.%09ds instead of %d.%09ds\n",
		     gpsTest.gpsSeconds, gpsTest.gpsNanoSeconds, gps2.gpsSeconds, gps2.gpsNanoSeconds);
      return FAILURE;
    }
    /* test over-run in ns-position is handled properly */
    SUB (LALAddFloatToGPS (&status, &gpsTest, &gps2, deltaT1), &status);	/* over-run in positive sense */
    gps3.gpsSeconds = 717076574;
    gps3.gpsNanoSeconds = 401265590;
    if ( (gpsTest.gpsSeconds != gps3.gpsSeconds) || (gpsTest.gpsNanoSeconds != gps3.gpsNanoSeconds) ) {
      LALPrintError ("Failure in 3.) LALAddFloatToGPS(): got %d.%09ds instead of %d.%09ds\n",
		     gpsTest.gpsSeconds, gpsTest.gpsNanoSeconds, gps2.gpsSeconds, gps2.gpsNanoSeconds);
      return FAILURE;
    }
    SUB (LALAddFloatToGPS (&status, &gpsTest, &gps1, deltaT1), &status);	/* over-run in negative sense */
    gps3.gpsSeconds = 715615153;
    gps3.gpsNanoSeconds = 712343412;
    if ( (gpsTest.gpsSeconds != gps3.gpsSeconds) || (gpsTest.gpsNanoSeconds != gps3.gpsNanoSeconds) ) {
      LALPrintError ("Failure in 4.) LALAddFloatToGPS(): got %d.%09ds instead of %d.%09ds\n",
		     gpsTest.gpsSeconds, gpsTest.gpsNanoSeconds, gps2.gpsSeconds, gps2.gpsNanoSeconds);
      return FAILURE;
    }


  } /* testing LALAddFloatToGPS() and LALDeltaFloatGPS() */
  /*----------------------------------------------------------------------*/

  /* try passing the same pointer to struct as input and output */
  gps.gpsSeconds     =    0;
  gps.gpsNanoSeconds =    0;
  deltaT.seconds     =   13;
  deltaT.nanoSeconds =  100;
  expected_gps.gpsSeconds =  13;
  expected_gps.gpsNanoSeconds = 100;


  /* END: test of LALIncrementGPS() */

  /*
   * TEST No. 2 -- test of LAlDecrementGPS()
   */
  gps.gpsSeconds     =     100;
  gps.gpsNanoSeconds =      70;
  deltaT.seconds     =       1;
  deltaT.nanoSeconds =       1;
  expected_gps.gpsSeconds =     99;
  expected_gps.gpsNanoSeconds = 69;

  if (!decrement_gps_ok(&status, &deltaT, &gps, &expected_gps))
    return FAILURE;


  gps.gpsSeconds     =     100;
  gps.gpsNanoSeconds =      70;
  deltaT.seconds     =     100;
  deltaT.nanoSeconds =      70;
  expected_gps.gpsSeconds =     0;
  expected_gps.gpsNanoSeconds = 0;

  if (!decrement_gps_ok(&status, &deltaT, &gps, &expected_gps))
    return FAILURE;

  gps.gpsSeconds     =    100;
  gps.gpsNanoSeconds =     70;
  deltaT.seconds     =     90;
  deltaT.nanoSeconds =     80;
  expected_gps.gpsSeconds     =  9;
  expected_gps.gpsNanoSeconds = oneBillion + 70 - 80;

  if (!decrement_gps_ok(&status, &deltaT, &gps, &expected_gps))
    return FAILURE;

  /* expect an error on this next one */
  gps.gpsSeconds     = 100;
  gps.gpsNanoSeconds =  70;
  deltaT.seconds     = 100;
  deltaT.nanoSeconds =  80;

  LALDecrementGPS(&status, &decremented_gps, &gps, &deltaT);

  if (status.statusCode > 0 && lalDebugLevel > 0)
    {
      if (status.statusCode == DATEH_EDECRTIMETOOLARGE) /* expected error */
        {
          fprintf(stderr, "failed with status code %d as expected",
                  DATEH_EDECRTIMETOOLARGE);
          REPORTSTATUS(&status);
        }
      else /* some other error */
        {
          fprintf(stderr, "TestIncrementGPS: LALDecrementGPS() failed, line %i, %s\n",
                  __LINE__, LALTESTINCREMENTGPSC);
          REPORTSTATUS(&status);
          return status.statusCode;
        }
    }

  /* END: test of LALDecrementGPS() */


  /*
   * TEST 3: test of LALDeltaGPS()
   *  (gee, maybe I should have made LIGOTimeGPS and LALTimeInterval
   *   into a union)
   */
  gps.gpsSeconds     = 1234;
  gps.gpsNanoSeconds = 4321;
  another_gps.gpsSeconds     = 234;
  another_gps.gpsNanoSeconds = 321;
  expected_deltaT.seconds     = 1000;
  expected_deltaT.nanoSeconds = 4000;

  LALDeltaGPS(&status, &deltaT, &gps, &another_gps);

  if (deltaT.seconds != expected_deltaT.seconds ||
      deltaT.nanoSeconds != expected_deltaT.nanoSeconds)
    {
      if (verbose_p)
        {
          fprintf(stderr, "LALDeltaGPS() failed\n");
          fprintf(stderr, "  expected: %9d:%09d\n", expected_deltaT.seconds,
                  expected_deltaT.nanoSeconds);
          fprintf(stderr, "       got: %9d:%09d\n", deltaT.seconds,
                  deltaT.nanoSeconds);
        }

      return FAILURE;
    }


  /* just reverse order of inputs */
  expected_deltaT.seconds     = -1000;
  expected_deltaT.nanoSeconds = -4000;

  LALDeltaGPS(&status, &deltaT, &another_gps, &gps);

  if (deltaT.seconds != expected_deltaT.seconds ||
      deltaT.nanoSeconds != expected_deltaT.nanoSeconds)
    {
      if (verbose_p)
        {
          fprintf(stderr, "LALDeltaGPS() failed\n");
          fprintf(stderr, "  expected %9d:%09d\n", expected_deltaT.seconds,
                  expected_deltaT.nanoSeconds);
          fprintf(stderr, "       got %9d:%09d\n", deltaT.seconds,
                  deltaT.nanoSeconds);
        }

      return FAILURE;
    }

  /* and a trivial case */
  expected_deltaT.seconds     = 0;
  expected_deltaT.nanoSeconds = 0;

  LALDeltaGPS(&status, &deltaT, &gps, &gps);

  if (deltaT.seconds != expected_deltaT.seconds ||
      deltaT.nanoSeconds != expected_deltaT.nanoSeconds)
    {
      if (verbose_p)
        {
          fprintf(stderr, "LALDeltaGPS() failed\n");
          fprintf(stderr, "  expected %9d:%09d\n", expected_deltaT.seconds,
                  expected_deltaT.nanoSeconds);
          fprintf(stderr, "       got %9d:%09d\n", deltaT.seconds,
                  deltaT.nanoSeconds);
        }

      return FAILURE;
    }

  return SUCCESS;
} /* END: main() */



static int print_compare_errmsg_maybe(const char *msg,
                                      LALGPSCompareResult expected_val,
                                      LALGPSCompareResult got_val)
{
  if (verbose_p)
    {
      char long_msg[LONGESTSTR];
      const char *earlier_str = "LALGPS_EARLIER";
      const char *equal_str = "LALGPS_EQUAL";
      const char *later_str = "LALGPS_LATER";

      if (strlen(msg) > LONGESTSTR - 1)
        {
          fprintf(stderr, "AARRGGHHH!! string too long\n");
          exit(69);
        }

      snprintf(long_msg, LONGESTSTR, "%s; expected ", msg);

      switch (expected_val)
        {
        case LALGPS_EARLIER:
          strncat(long_msg, earlier_str, LONGESTSTR - 1 - strlen(long_msg));
          break;

        case LALGPS_EQUAL:
          strncat(long_msg, equal_str, LONGESTSTR - 1 - strlen(long_msg));
          break;

        case LALGPS_LATER:
          strncat(long_msg, later_str, LONGESTSTR - 1 - strlen(long_msg));
          break;
        }

      strncat(long_msg, ", but got ", LONGESTSTR - 1 - strlen(long_msg));

      switch (got_val)
        {
        case LALGPS_EARLIER:
          strncat(long_msg, earlier_str, LONGESTSTR - 1 - strlen(long_msg));
          break;

        case LALGPS_EQUAL:
          strncat(long_msg, equal_str, LONGESTSTR - 1 - strlen(long_msg));
          break;

        case LALGPS_LATER:
          strncat(long_msg, later_str, LONGESTSTR - 1 - strlen(long_msg));
          break;
        }

      return fprintf(stderr, "%s\n\n", long_msg);
    }
  else
    {
      return 0;
    }
}



static BOOLEAN compare_gps_ok(LALStatus *status,
                              const LIGOTimeGPS *p_gps1,
                              const LIGOTimeGPS *p_gps2,
                              LALGPSCompareResult  expected_val)
{
  LALGPSCompareResult rslt;

  LALCompareGPS(status, &rslt, p_gps1, p_gps2);

  if (status->statusCode && lalDebugLevel > 0)
      {
        fprintf(stderr, "TestIncrementGPS: compare_gps_ok: LALCompareGPS() failed, line %i, %s\n",
                __LINE__, LALTESTINCREMENTGPSC);
        REPORTSTATUS(status);
        exit (status->statusCode);
      }

  if (rslt != expected_val)
    {
      print_compare_errmsg_maybe("LALCompareGPS() failed",
                                 expected_val, rslt);
      return FALSE;
    }
  else
    {
      return TRUE;
    }
}


static int sprint_gps_time(char *str, const LIGOTimeGPS *p_gps)
{
  return snprintf(str, LONGESTSTR - 1 - strlen(str),
                  "%9d:%09d", p_gps->gpsSeconds, p_gps->gpsNanoSeconds);
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




static int print_incr_errmsg_maybe(const char *msg,
                                   const LIGOTimeGPS *p_expected_gps,
                                   const LIGOTimeGPS *p_got_gps)
{
  if (verbose_p)
    {
      char long_msg[LONGESTSTR];

      if (strlen(msg) > LONGESTSTR - 1)
        {
          fprintf(stderr, "AARRGGHHH!! string too long\n");
          exit(69);
        }

      snprintf(long_msg, LONGESTSTR, "%s; expected ", msg);
      sprint_gps_time(long_msg, p_expected_gps);
      snprintf(long_msg, LONGESTSTR, ", got ");
      sprint_gps_time(long_msg, p_got_gps);

      return fprintf(stderr, "%s\n\n", long_msg);
    }
  else
    {
      return 0;
    }
} /* END: print_incr_errmsg_maybe() */



static BOOLEAN increment_gps_ok(LALStatus *status,
                                LALTimeInterval *p_delta,
                                LIGOTimeGPS *p_init_gps,
                                const LIGOTimeGPS *p_expected_gps)
{
  LALGPSCompareResult cmprslt;

  if (verbose_p)
    {
      printf("init_gps:     %9d:%09d\n", p_init_gps->gpsSeconds,
             p_init_gps->gpsNanoSeconds);
      printf("delta_t:      %9d:%09d\n", p_delta->seconds,
             p_delta->nanoSeconds);
      printf("expected_gps: %9d:%09d\n", p_expected_gps->gpsSeconds,
             p_expected_gps->gpsNanoSeconds);
    }

  LALIncrementGPS(status, p_init_gps, p_init_gps, p_delta);
  LALCompareGPS(status, &cmprslt, p_expected_gps, p_init_gps);

  if (status->statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr, "TestIncrementGPS: increment_gps_ok: LALIncrementGPS() failed, line %i, %s\n",
              __LINE__, LALTESTINCREMENTGPSC);
      REPORTSTATUS(status);
      exit (status->statusCode);
    }

  if (verbose_p)
    {
      printf("got_gps:      %9d:%09d\n", p_init_gps->gpsSeconds,
             p_init_gps->gpsNanoSeconds);
    }

  if (cmprslt != LALGPS_EQUAL)
    {
      print_incr_errmsg_maybe("LALIncrementGPS() failed",
                              p_expected_gps, p_init_gps);
      return FALSE;
    }
  else
    {
      return TRUE;
    }
} /* END: increment_gps_ok() */



static BOOLEAN decrement_gps_ok(LALStatus *status,
                                LALTimeInterval *p_delta,
                                LIGOTimeGPS *p_init_gps,
                                const LIGOTimeGPS *p_expected_gps)
{
  LALGPSCompareResult cmprslt;

  if (verbose_p)
    {
      printf("HELLO THERE\n");
      printf("lalDebugLevel = %d\n", lalDebugLevel);
      printf("init_gps:     %9d:%09d\n", p_init_gps->gpsSeconds,
             p_init_gps->gpsNanoSeconds);
      printf("delta_t:      %9d:%09d\n", p_delta->seconds,
             p_delta->nanoSeconds);
      printf("expected_gps: %9d:%09d\n", p_expected_gps->gpsSeconds,
             p_expected_gps->gpsNanoSeconds);
    }

  LALDecrementGPS(status, p_init_gps, p_init_gps, p_delta);

  if (verbose_p)
    {
      printf("got_gps:      %9d:%09d\n", p_init_gps->gpsSeconds,
             p_init_gps->gpsNanoSeconds);
    }

  /* um, if this fails.... */
  LALCompareGPS(status, &cmprslt, p_expected_gps, p_init_gps);

  if (cmprslt != LALGPS_EQUAL)
    {
      print_incr_errmsg_maybe("LALIncrementGPS() failed",
                              p_expected_gps, p_init_gps);
      return FALSE;
    }
  else
    {
      return TRUE;
    }
} /* END: decrement_gps_ok() */
