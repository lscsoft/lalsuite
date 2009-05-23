/*
 * Copyright (C) 2009 Reinhard Prix
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

/*********************************************************************************/
/** \author R. Prix
 * \file
 * \brief
 *  unit-test for Fstat_v3 module
 *
 *********************************************************************************/
#include "config.h"

/* System includes */
#include <stdio.h>
#include <time.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

/* GSL includes */
#include <gsl/gsl_rng.h>


/* LAL-includes */
#include <lal/AVFactories.h>
#include <lal/LALStdlib.h>
#include <lal/LogPrintf.h>
#include <lal/TimeSeries.h>
#include <lal/ComplexFFT.h>

#include <lal/ComputeFstat.h>
#include "Fstat_v3.h"


#include <lal/lalGitID.h>
#include <lalappsGitID.h>

#include <lalapps.h>

/* local includes */

/*---------- DEFINES ----------*/
#define TEST_PASSED     0
#define TEST_FAILED     1
#define TEST_ABORTED    2

/* ---------- local types ---------- */

/*---------- Global variables ----------*/

/* ----- User-variables: can be set from config-file or command-line */
/* ---------- local prototypes ---------- */
int main(int argc, char *argv[]);

int test_XLALSFTVectorToCOMPLEX8TimeSeries(void);

/*---------- empty initializers ---------- */
static LALUnit empty_LALUnit;

int lalDebugLevel = 0;

/*----------------------------------------------------------------------*/
/* Main Function starts here */
/*----------------------------------------------------------------------*/
/**
 * MAIN function
 */
int main(int argc, char *argv[])
{
  const CHAR *fn = "testFstat_v3";

  int res1;

   if ( (res1 = test_XLALSFTVectorToCOMPLEX8TimeSeries()) != TEST_PASSED )
    {
      LogPrintf (LOG_CRITICAL, "%s: test of XLALSFTVectorToCOMPLEX8TimeSeries() failed\n\n");
      return res1;
    }
  else
    {
      LogPrintf (LOG_NORMAL, "%s: test of XLALSFTVectorToCOMPLEX8TimeSeries() succeeded\n\n");
    }


  return TEST_PASSED;

} /* main() */



/** Unit test for XLALSFTVectorToCOMPLEX8TimeSeries() function.
 * Returns: TEST_PASSED, TEST_FAILED, or TEST_ABORTED
 */
int
test_XLALSFTVectorToCOMPLEX8TimeSeries ( void )
{
  const CHAR *fn = "test_XLALSFTVectorToCOMPLEX8TimeSeries()";

  REAL4TimeSeries *inputTS;
  LIGOTimeGPS epoch = { 714180733, 0 };
  REAL8 fmax = 10.0;
  REAL8 deltaT = 1.0/(2.0 * fmax);
  REAL8 duration = 10.0 * 3600;	/* 10 hours */
  UINT4 numSamples = (UINT4)(duration / deltaT);	/* round down */
  REAL8 Tsft = 1800.0;
  UINT4 numSFTs = (UINT4)(duration / Tsft) - 2;		/* leave at least 2 Tsft as gaps */
  UINT4 i_gap1, i_gap2, i_gap3;				/* leave 3 gaps at random indices */
  time_t now;
  gsl_rng * r;

  /* prepare random-number generator */
  if ( (r = gsl_rng_alloc (gsl_rng_taus)) == NULL )
    {
      XLALPrintError ("%s: failed to gsl_rng_alloc (gsl_rng_taus)\n");
      return TEST_ABORTED;
    }
  now = time(NULL);
  gsl_rng_set (r, (unsigned long int)now);	/* sets random-number seed */



  /* ----- (non-heterodyned) random noise timeseries with gaps ----- */
  if ( (inputTS = XLALCreateREAL4TimeSeries ("input TS", &epoch, 0, deltaT, &empty_LALUnit, numSamples )) == NULL )
    {
      XLALPrintError ("%s: XLALCreateREAL4TimeSeries() failed for numSamples = %d\n", numSamples );
      return TEST_ABORTED;
    }

  /* initialized to zero (to allow for gaps) */
  memset ( inputTS->data->data, 0, inputTS->data->length * sizeof ( inputTS->data->data[0] ) );



  /* free memory */
  gsl_rng_free (r);
  XLALDestroyREAL4TimeSeries ( inputTS );

  return TEST_PASSED;

} /* test_XLALSFTVectorToCOMPLEX8TimeSeries() */
