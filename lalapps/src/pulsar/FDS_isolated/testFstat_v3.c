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
#include <gsl/gsl_randist.h>


/* LAL-includes */
#include <lal/AVFactories.h>
#include <lal/LALStdlib.h>
#include <lal/LogPrintf.h>
#include <lal/TimeSeries.h>
#include <lal/ComplexFFT.h>

#include <lal/GeneratePulsarSignal.h>
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

/* ----- Macros ----- */
#define MYMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MYMIN(x,y) ( (x) < (y) ? (x) : (y) )


/* ---------- local types ---------- */

/*---------- Global variables ----------*/

/* ----- User-variables: can be set from config-file or command-line */
/* ---------- local prototypes ---------- */
int main(int argc, char *argv[]);

int test_XLALSFTVectorToCOMPLEX8TimeSeries(void);
int XLALgenerateRandomData ( REAL4TimeSeries **ts, SFTVector **sfts );

/*---------- empty initializers ---------- */
static LALUnit empty_LALUnit;

int lalDebugLevel = 1;

/*----------------------------------------------------------------------*/
/* Main Function starts here */
/*----------------------------------------------------------------------*/
/**
 * MAIN function
 */
int main(int argc, char *argv[])
{
  const CHAR *fn = "testFstat_v3";
  REAL4TimeSeries *ts = NULL;
  SFTVector *sfts = NULL;
  int res1;


  LogPrintf (LOG_NORMAL, "%s: Now testing XLALSFTVectorToCOMPLEX8TimeSeries() ... ", fn);
  if ( (res1 = XLALgenerateRandomData(&ts, &sfts)) != TEST_PASSED )
    {
      LogPrintfVerbatim (LOG_CRITICAL, "failed.\n\n");
      return res1;
    }
  else
    {
      LogPrintfVerbatim (LOG_NORMAL, "succeeded.\n\n");
    }

  XLALDestroyREAL4TimeSeries ( ts );
  XLALDestroySFTVector ( sfts );

  LALCheckMemoryLeaks();

  return TEST_PASSED;

} /* main() */



/**
 * function to generate random time-series with gaps, and corresponding SFTs
 */
int
XLALgenerateRandomData ( REAL4TimeSeries **ts, SFTVector **sfts )
{
  const CHAR *fn = "XLALgenerateRandomData()";

  LIGOTimeGPS epoch0 = { 714180733, 0 };
  UINT4 numSFTs = 19;

  REAL8 fmax = 10.0;
  REAL8 Tsft = 1800.0;
  REAL8 dFreq;

  REAL8 deltaT;
  REAL8 duration;
  UINT4 numFreqBins;
  REAL8 BandEff;		/* effective frequency band of SFTs (integer number of bins) */
  UINT4 numSamplesSFT;		/* number of samples in one SFTs (either time-domain or frequency-domain) */
  UINT4 numSamplesTS;		/* number of samples in the full time-series (including gaps) */

  REAL4TimeSeries *outTS;	/* input timeseries, generated randomly, with 2 random gaps */
  REAL4 *TSdata;		/* convenience pointer to outTS->data->data */
  SFTVector *outSFTs;		/* SFT vector created from input timeseries */
  LIGOTimeGPSVector *timestampsSFT;	/* vector of SFT timestamps */

  UINT4 i1, i2, iGap1, iGap2;	/* leave 2 gaps at random indices */
  UINT4 numBinsGap1, numBinsGap2;
  UINT4 thisBin;

  UINT4 iSFT, iBin;

  time_t now;
  gsl_rng * r;

  /* input sanity checks */
  if ( !ts || !sfts ) {
    XLALPrintError ("%s: NULL input pointers\n", fn );
    XLAL_ERROR (fn, XLAL_EINVAL );
  }
  if ( (*ts != NULL) || (*sfts != NULL) ) {
    XLALPrintError ("%s: output pointers not initialized to NULL.\n", fn);
    XLAL_ERROR (fn, XLAL_EINVAL );
  }

  /* prepare sampling constants */
  dFreq         = 1.0 / Tsft;
  numFreqBins   = (UINT4)(fmax / dFreq);
  BandEff       = numFreqBins * dFreq;
  deltaT        = 1.0 / (2.0 * BandEff);
  numSamplesSFT = 2 * numFreqBins;		/* real-valued timeseries */

  /* prepare random-number generator */
  if ( (r = gsl_rng_alloc (gsl_rng_taus2)) == NULL )
    {
      XLALPrintError ("%s: failed to gsl_rng_alloc (gsl_rng_taus2)\n");
      return TEST_ABORTED;
    }
  now = time(NULL);
  gsl_rng_set (r, (unsigned long int)now);	/* sets random-number seed */

  /* ----- set up 2 unique gaps ----- */
  i1 = gsl_rng_uniform_int (r, numSFTs );	/* uniform integer between [0, n-1] inclusive */
  do {
    i2 = gsl_rng_uniform_int (r, numSFTs );
  } while (i2 == i1);
  iGap1 = MYMIN ( i1, i2 );
  iGap2 = MYMAX ( i1, i2 );

  /* 2 gaps of gap-lengths between 0 and 1/3 Tdata */
  numBinsGap1 = gsl_rng_uniform_int (r, numSamplesSFT * numSFTs / 3 );
  numBinsGap2 = gsl_rng_uniform_int (r, numSamplesSFT * numSFTs / 3 );

  numSamplesTS = numSFTs * numSamplesSFT + numBinsGap1 + numBinsGap2;
  duration = numSamplesTS * deltaT;

  /* ----- allocate timeseries ----- */
  if ( (outTS = XLALCreateREAL4TimeSeries ("input TS", &epoch0, 0, deltaT, &empty_LALUnit, numSamplesTS )) == NULL )
    {
      XLALPrintError ("%s: XLALCreateREAL4TimeSeries() failed for numSamples = %d\n", numSamplesTS );
      XLAL_ERROR ( fn, XLAL_EFUNC );
    }

  TSdata = outTS->data->data;
  /* initialize timeseries to zero (to allow for gaps) */
  memset ( TSdata, 0, outTS->data->length * sizeof ( outTS->data->data[0] ) );

  /* also set up corresponding SFT timestamps vector */
  if ( (timestampsSFT = XLALCalloc (1, sizeof( *timestampsSFT )) ) == NULL ) {
    XLALPrintError ("%s: Failed to XLALCalloc %d bytes\n", fn, sizeof( *timestampsSFT ));
    XLAL_ERROR ( fn, XLAL_ENOMEM );
  }

  timestampsSFT->length = numSFTs;
  if ( (timestampsSFT->data = XLALCalloc (numSFTs, sizeof (*timestampsSFT->data) )) == NULL) {
    XLALPrintError ("%s: Failed to allocate %d x %d bytes\n", fn, numSFTs, sizeof (*timestampsSFT->data) );
    XLALFree ( timestampsSFT );
    XLAL_ERROR ( fn, XLAL_ENOMEM );
  }

  /* ----- set up random-noise timeseries with gaps ---------- */
  thisBin = 0;	/* running bin-counter */
  for ( iSFT=0; iSFT < numSFTs; iSFT ++ )
    {
      /* skip gaps */
      if ( iSFT == iGap1 )
	thisBin += numBinsGap1;
      if ( iSFT == iGap2 )
	thisBin += numBinsGap2;

      /* record this SFT's timestamp */
      timestampsSFT->data[iSFT] = epoch0;
      timestampsSFT->data[iSFT].gpsSeconds += thisBin * deltaT;

      /* generate all data-points of this SFT */
      for ( iBin=0; iBin < numSamplesSFT; iBin++ ) {
	TSdata[thisBin++] = (REAL4)gsl_ran_gaussian (r, 1.0);	/* unit-variance Gaussian noise */
	/* TSdata[thisBin++] = (REAL4)( 2.0 * gsl_rng_uniform (r) - 1.0 );	*//* random-number in [-1, 1] */
      }
    } /* for iSFT < numSFTs */

  if ( thisBin != numSamplesTS ) {
    XLALPrintError ("\n%f: sanity check failed, thisBin = %d does not agree with numSamplesTS = %d\n", fn, thisBin, numSamplesTS );
    XLAL_ERROR ( fn, XLAL_EBADLEN );

  }

  /* ----- generate SFTs from this timeseries ---------- */
  {
    SFTParams sftParams = {0};
    LALStatus status = {0};

    sftParams.Tsft = Tsft;
    sftParams.timestamps = timestampsSFT;
    sftParams.noiseSFTs = NULL;
    sftParams.make_v2SFTs = 1;
    sftParams.window = NULL;

    outSFTs = NULL;
    LALSignalToSFTs(&status, &outSFTs, outTS, &sftParams);
    if ( status.statusCode ) {
      XLALPrintError("\n%s: LALSignalToSFTs() failed, with LAL error %d\n", fn, status.statusCode );
      XLAL_ERROR ( fn, XLAL_EFUNC );
    }
  } /* turn timeseries into SFTs */

  /* free memory */
  gsl_rng_free (r);
  XLALFree ( timestampsSFT->data );
  XLALFree ( timestampsSFT );

  /* return timeseries and SFTvector */
  (*ts)   = outTS;
  (*sfts) = outSFTs;

  return XLAL_SUCCESS;

} /* XLALgenerateRandomData() */
