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
/**
 * \author R. Prix
 * \file
 * \brief
 * unit-test for Fstat_v3 module
 *
 */
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
#include <lal/RealFFT.h>
#include <lal/ComplexFFT.h>
#include <lal/SFTfileIO.h>

#include <lal/GeneratePulsarSignal.h>
#include <lal/ComputeFstat.h>
#include "Fstat_v3.h"

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
int main(void);

int test_XLALSFTVectorToCOMPLEX8TimeSeries(void);
int XLALgenerateRandomData ( REAL4TimeSeries **ts, SFTVector **sfts );

void write_timeSeriesR4 (FILE *fp, const REAL4TimeSeries *series);

/*---------- empty initializers ---------- */
static LALUnit empty_LALUnit;
static LALStatus empty_LALStatus;


/*----------------------------------------------------------------------*/
/* Main Function starts here */
/*----------------------------------------------------------------------*/
/**
 * MAIN function
 */
int main(void)
{
  int res1;


  LogPrintf (LOG_DEBUG, "%s: Now testing XLALSFTVectorToCOMPLEX8TimeSeries() ... ", __func__);
  if ( (res1 = test_XLALSFTVectorToCOMPLEX8TimeSeries()) != TEST_PASSED )
    {
      LogPrintfVerbatim (LOG_CRITICAL, "failed.\n\n");
      return res1;
    }
  else
    {
      LogPrintfVerbatim (LOG_DEBUG, "succeeded.\n\n");
    }

  LALCheckMemoryLeaks();

  return TEST_PASSED;

} /* main() */

/**
 * Unit-Test for function XLALSFTVectorToCOMPLEX8TimeSeries().
 * Generates random data (timeseries + SFTs), feeds the SFTs into XLALSFTVectorToCOMPLEX8TimeSeries()
 * and checks correctness of output timeseries.
 *
 * returns TEST_PASSED, TEST_FAILED or TEST_ABORTED
 */
int
test_XLALSFTVectorToCOMPLEX8TimeSeries(void)
{
  SFTVector *sfts = NULL;
  SFTtype *lftTest = NULL;

  REAL4TimeSeries *tsOrig = NULL;
  COMPLEX8TimeSeries *tsTest = NULL;

  COMPLEX8FFTPlan *LFTplanC8;
  LALStatus status = empty_LALStatus;
  UINT4 numSamples;

  if ( XLALgenerateRandomData ( &tsOrig, &sfts ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALgenerateRandomData() failed!\n", __func__);
    return TEST_ABORTED;
  }


  if ( ( tsTest = XLALSFTVectorToCOMPLEX8TimeSeries ( sfts, NULL, NULL ) ) == NULL ) {
    XLALPrintError ("%s: call to XLALSFTVectorToCOMPLEX8TimeSeries() failed. xlalErrrno = %d\n", __func__, xlalErrno );
    return TEST_ABORTED;
  }

  numSamples = tsTest->data->length;	/* number of (complex) time- or frequency samples */

  /* ----- prepare long-timebaseline FFT ---------- */
  status = empty_LALStatus;
  LALCreateSFTtype (&status, &lftTest, numSamples);
  if ( status.statusCode ) {
    XLALPrintError("\n%s: LALCreateSFTtype() failed, with LAL error %d\n", __func__, status.statusCode );
    return TEST_ABORTED;
  }

  strcpy ( lftTest->name, tsTest->name );
  lftTest->f0 = 0;
  lftTest->epoch = tsTest->epoch;
  lftTest->deltaF = 1.0 / ( numSamples * tsTest->deltaT );

  /* ----- compute FFTW on computed C8 timeseries ---------- */
  if ( (LFTplanC8 = XLALCreateForwardCOMPLEX8FFTPlan( numSamples, 0 )) == NULL )
    {
      XLALPrintError ( "%s: XLALCreateForwardCOMPLEX8FFTPlan(%d, ESTIMATE ) failed! errno = %d!\n", __func__, numSamples, xlalErrno );
      return TEST_ABORTED;
    }

  if ( XLALCOMPLEX8VectorFFT( lftTest->data, tsTest->data, LFTplanC8 ) != XLAL_SUCCESS )
    {
      XLALPrintError ( "%s: XLALCOMPLEX8VectorFFT() failed! errno = %d!\n", __func__, xlalErrno );
      return TEST_ABORTED;
    }

  if ( XLALReorderFFTWtoSFT (lftTest->data) != XLAL_SUCCESS )
    {
      XLALPrintError ( "%s: XLALReorderFFTWtoSFT() failed! errno = %d!\n", __func__, xlalErrno );
      return TEST_ABORTED;
    }


  /* debug output */
  {
    const CHAR *fnameTS = "testFstat_v3-timeseries.dat";
    const CHAR *fnameLFT = "testFstat_v3-LFT.sft";
    FILE *fp;

    if ( (fp = fopen(fnameTS, "wb")) == NULL ) {
      LogPrintf (LOG_CRITICAL, "%s: failed to open '%s' for writing.\n", __func__, fnameTS );
      return TEST_ABORTED;
    }
    write_timeSeriesR4 (fp, tsOrig );
    fclose ( fp );

    /*
    LALWriteSFTVector2Dir (&status, sfts, "./", "test SFTs with random data", "testFstat_v3" );
    if ( status.statusCode ) {
      XLALPrintError("\n%s: LALWriteSFTVector2Dir() failed, with LAL error %d\n", __func__, status.statusCode );
      return TEST_ABORTED;
    }
    */

    status = empty_LALStatus;
    LALWriteSFT2file ( &status, lftTest, fnameLFT, "test data LFT from SFTs");
    if ( status.statusCode ) {
      XLALPrintError("\n%s: LALWriteSFT2file() failed, with LAL error %d\n", __func__, status.statusCode );
      return TEST_ABORTED;
    }

  } /* debug output */

  /* free memory */
  XLALDestroyREAL4TimeSeries ( tsOrig );
  XLALDestroyCOMPLEX8TimeSeries ( tsTest );

  XLALDestroySFTVector ( sfts );
  XLALDestroyCOMPLEX8FFTPlan ( LFTplanC8 );

  status = empty_LALStatus;
  LALDestroySFTtype ( &status, &lftTest );

  return TEST_PASSED;

} /* test_XLALSFTVectorToCOMPLEX8TimeSeries() */


/**
 * function to generate random time-series with gaps, and corresponding SFTs
 */
int
XLALgenerateRandomData ( REAL4TimeSeries **ts, SFTVector **sfts )
{
  LIGOTimeGPS epoch0 = { 714180733, 0 };
  UINT4 numSFTs = 19;

  REAL8 f_max = 10.0;
  REAL8 Tsft = 1800.0;
  REAL8 dFreq;

  REAL8 deltaT;
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
    XLALPrintError ("%s: NULL input pointers\n", __func__);
    XLAL_ERROR ( XLAL_EINVAL );
  }
  if ( (*ts != NULL) || (*sfts != NULL) ) {
    XLALPrintError ("%s: output pointers not initialized to NULL.\n", __func__);
    XLAL_ERROR ( XLAL_EINVAL );
  }

  /* prepare sampling constants */
  dFreq         = 1.0 / Tsft;
  numFreqBins   = (UINT4)(f_max / dFreq);
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

  /* ----- allocate timeseries ----- */
  if ( (outTS = XLALCreateREAL4TimeSeries ("H1:test timeseries", &epoch0, 0, deltaT, &empty_LALUnit, numSamplesTS )) == NULL )
    {
      XLALPrintError ("%s: XLALCreateREAL4TimeSeries() failed for numSamples = %d\n", numSamplesTS );
      XLAL_ERROR ( XLAL_EFUNC );
    }

  TSdata = outTS->data->data;
  /* initialize timeseries to zero (to allow for gaps) */
  memset ( TSdata, 0, outTS->data->length * sizeof ( outTS->data->data[0] ) );

  /* also set up corresponding SFT timestamps vector */
  if ( (timestampsSFT = XLALCalloc (1, sizeof( *timestampsSFT )) ) == NULL ) {
    XLALPrintError ("%s: Failed to XLALCalloc %d bytes\n", __func__, sizeof( *timestampsSFT ));
    XLAL_ERROR ( XLAL_ENOMEM );
  }

  timestampsSFT->length = numSFTs;
  if ( (timestampsSFT->data = XLALCalloc (numSFTs, sizeof (*timestampsSFT->data) )) == NULL) {
    XLALPrintError ("%s: Failed to allocate %d x %d bytes\n", __func__, numSFTs, sizeof (*timestampsSFT->data) );
    XLALFree ( timestampsSFT );
    XLAL_ERROR ( XLAL_ENOMEM );
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
    XLALPrintError ("\n%f: sanity check failed, thisBin = %d does not agree with numSamplesTS = %d\n", __func__, thisBin, numSamplesTS );
    XLAL_ERROR ( XLAL_EBADLEN );

  }

  /* ----- generate SFTs from this timeseries ---------- */
  {
    SFTParams sftParams = empty_SFTParams;
    LALStatus status = empty_LALStatus;

    sftParams.Tsft = Tsft;
    sftParams.timestamps = timestampsSFT;
    sftParams.noiseSFTs = NULL;
    sftParams.window = NULL;

    outSFTs = NULL;
    LALSignalToSFTs(&status, &outSFTs, outTS, &sftParams);
    if ( status.statusCode ) {
      XLALPrintError("\n%s: LALSignalToSFTs() failed, with LAL error %d\n", __func__, status.statusCode );
      XLAL_ERROR ( XLAL_EFUNC );
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


/* write a time-series into a text-file */
void
write_timeSeriesR4 (FILE *fp, const REAL4TimeSeries *series)
{
  REAL8 timestamp0, timestamp;
  UINT4 i;

  if (series == NULL)
  {
    printf ("\nempty input!\n");
    return;
  }

  timestamp0 = 1.0*series->epoch.gpsSeconds + series->epoch.gpsNanoSeconds * 1.0e-9;

  for( i = 0; i < series->data->length; i++)
  {
    timestamp = timestamp0 + 1.0 * i * series->deltaT;
    fprintf( fp, "%16.9f %e\n", timestamp, series->data->data[i] );

  }

  return;

} /* write_timeSeriesR4() */
