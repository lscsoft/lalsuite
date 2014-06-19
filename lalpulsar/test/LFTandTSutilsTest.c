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
 * \brief Unit tests for LFTandTSutils module
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
#include <lal/Units.h>

#include <lal/GeneratePulsarSignal.h>
#include <lal/ComputeFstat.h>
#include <lal/LFTandTSutils.h>

/* local includes */

/*---------- DEFINES ----------*/
/* ----- Macros ----- */
#define MYMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MYMIN(x,y) ( (x) < (y) ? (x) : (y) )
#define RELERR(x,y) ( ( (x) - (y) ) / ( 0.5 * (fabs(x) + fabs(y)) + 2) )
#define SQ(x) ( (x) * (x) )
/* ---------- local types ---------- */

/*---------- Global variables ----------*/
static LALUnit emptyLALUnit;

/* ----- User-variables: can be set from config-file or command-line */
/* ---------- local prototypes ---------- */
int main(void);


int test_XLALSFTVectorToLFT(void);
int XLALgenerateRandomData ( REAL4TimeSeries **ts, SFTVector **sfts );
int write_SFTdata (const char *fname, const SFTtype *sft);
SFTVector *XLALDuplicateSFTVector ( const SFTVector *sftsIn );
int XLALCompareSFTs ( const SFTtype *sft1, const SFTtype *sft2, const VectorComparison *tol );

/*----------------------------------------------------------------------*/
/* Main Function starts here */
/*----------------------------------------------------------------------*/
/**
 * MAIN function
 */
int main(void)
{

  XLAL_CHECK ( test_XLALSFTVectorToLFT() == XLAL_SUCCESS, XLAL_EFUNC );

  LALCheckMemoryLeaks();

  return XLAL_SUCCESS;

} /* main() */

/**
 * Unit-Test for function XLALSFTVectorToLFT().
 * Generates random data (timeseries + corresponding SFTs),
 * then feeds the SFTs into XLALSFTVectorToLFT()
 * and checks correctness of output Fourier transform by
 * comparing to FT of original real-valued timeseries.
 *
 * \note This indirectly also checks XLALSFTVectorToCOMPLEX8TimeSeries()
 * which is used by XLALSFTVectorToLFT().
 *
 * returns XLAL_SUCCESS on success, XLAL-error otherwise
 */
int
test_XLALSFTVectorToLFT ( void )
{
  // ----- generate real-valued random timeseries and corresponding SFTs
  REAL4TimeSeries *tsR4 = NULL;
  SFTVector *sfts0 = NULL;
  XLAL_CHECK ( XLALgenerateRandomData ( &tsR4, &sfts0 ) == XLAL_SUCCESS, XLAL_EFUNC );
  UINT4 numSamplesR4 = tsR4->data->length;
  REAL8 dt_R4 = tsR4->deltaT;
  REAL8 TspanR4 = numSamplesR4 * dt_R4;
  // ----- consider only the frequency band [3Hz, 4Hz]
  REAL8 out_fMin = 3;
  REAL8 out_Band = 1;
  SFTVector *sftsBand;
  XLAL_CHECK ( (sftsBand = XLALExtractBandFromSFTVector ( sfts0, out_fMin, out_Band )) != NULL, XLAL_EFUNC );
  XLALDestroySFTVector ( sfts0 );

  // ----- 1) compute FFT on original REAL4 timeseries -----
  REAL4FFTPlan *planR4;
  XLAL_CHECK ( (planR4 = XLALCreateForwardREAL4FFTPlan ( numSamplesR4, 0 )) != NULL, XLAL_EFUNC );

  SFTtype *lft0;
  XLAL_CHECK ( (lft0 = XLALCreateSFT ( numSamplesR4/2+1 )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( XLALREAL4ForwardFFT ( lft0->data, tsR4->data, planR4 ) == XLAL_SUCCESS, XLAL_EFUNC );

  strcpy ( lft0->name, tsR4->name );
  lft0->f0 = 0;
  lft0->epoch = tsR4->epoch;
  REAL8 dfLFT = 1.0 / TspanR4;
  lft0->deltaF = dfLFT;

  // ----- extract frequency band of interest
  SFTVector XLAL_INIT_DECL(lftR4V0);
  lftR4V0.length = 1;
  lftR4V0.data = lft0;
  SFTVector *lftR4V;
  XLAL_CHECK ( (lftR4V = XLALExtractBandFromSFTVector ( &lftR4V0, out_fMin, out_Band )) != NULL, XLAL_EFUNC );
  XLALDestroySFT ( lft0 );
  SFTtype *lftR4 = &(lftR4V->data[0]);

  for ( UINT4 k=0; k < lftR4->data->length; k ++ ) {
    lftR4->data->data[k] *= dt_R4;
  }

  // ----- 2) compute LFT directly from SFTs ----------
  SFTtype *lftSFTs0;
  XLAL_CHECK ( (lftSFTs0 = XLALSFTVectorToLFT ( sftsBand, 1 )) != NULL, XLAL_EFUNC );
  XLALDestroySFTVector ( sftsBand );

  // ----- re-extract frequency band of interest
  SFTVector XLAL_INIT_DECL(lftSFTsV0);
  lftSFTsV0.length = 1;
  lftSFTsV0.data = lftSFTs0;
  SFTVector *lftSFTsV;
  XLAL_CHECK ( (lftSFTsV = XLALExtractBandFromSFTVector ( &lftSFTsV0, out_fMin, out_Band )) != NULL, XLAL_EFUNC );
  XLALDestroySFT ( lftSFTs0 );
  SFTtype *lftSFTs = &(lftSFTsV->data[0]);

  if ( lalDebugLevel & LALINFO )
  {   // ----- write debug output
    XLAL_CHECK ( XLALdumpREAL4TimeSeries ( "TS_R4.dat", tsR4 ) == XLAL_SUCCESS, XLAL_EFUNC );

    XLAL_CHECK ( write_SFTdata ("LFT_R4T.dat", lftR4 ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK ( write_SFTdata ("LFT_SFTs.dat", lftSFTs ) == XLAL_SUCCESS, XLAL_EFUNC );

  } // end: debug output

  // ========== compare resulting LFTs ==========
  VectorComparison XLAL_INIT_DECL(tol);
  tol.relErr_L1 	= 0.5;
  tol.relErr_L2		= 5e-2;
  tol.angleV 		= 4e-2;	// rad
  tol.relErr_atMaxAbsx 	= 2e-3;
  tol.relErr_atMaxAbsy  = 2e-3;

  XLAL_CHECK ( XLALCompareSFTs ( lftR4, lftSFTs, &tol ) == XLAL_SUCCESS, XLAL_EFUNC );

  // ---------- free memory ----------
  XLALDestroyREAL4TimeSeries ( tsR4 );
  XLALDestroyREAL4FFTPlan ( planR4 );

  XLALDestroySFTVector ( lftR4V );
  XLALDestroySFTVector ( lftSFTsV );

  return XLAL_SUCCESS;

} // test_XLALSFTVectorToCOMPLEX8TimeSeries()


/**
 * function to generate random time-series with gaps, and corresponding SFTs
 */
int
XLALgenerateRandomData ( REAL4TimeSeries **ts, SFTVector **sfts )
{
  /* input sanity checks */
  XLAL_CHECK ( (ts != NULL) && ( (*ts) == NULL ), XLAL_EINVAL );
  XLAL_CHECK ( (sfts != NULL) && ( (*sfts) == NULL ), XLAL_EINVAL );

  // test parameters
  LIGOTimeGPS epoch0 = { 714180733, 0 };
  UINT4 numSFTs = 20;
  REAL8 Tsft = 1000.0;

  // noise sigma
  REAL8 sigmaN = 0;

  // frequency and amplitude of sinuoid signal
  REAL8 freqS = LAL_PI;
  REAL8 ampS = 0.01;
  REAL8 phi0S = 1.0;

  /* prepare sampling constants */
  REAL8 numR4SamplesPerSFT = 2 * 5000;
  REAL8 dtR4 = (Tsft / numR4SamplesPerSFT);

  UINT4 numR4SamplesTS = numSFTs * numR4SamplesPerSFT;

  /* ----- allocate timeseries ----- */
  REAL4TimeSeries *outTS;	// input timeseries
  XLAL_CHECK ( (outTS = XLALCreateREAL4TimeSeries ("H1:test timeseries", &epoch0, 0, dtR4, &emptyLALUnit, numR4SamplesTS )) != NULL, XLAL_EFUNC );

  REAL4 *TSdata = outTS->data->data;
  /* initialize timeseries to zero (to allow for gaps) */
  memset ( TSdata, 0, outTS->data->length * sizeof (*TSdata) );

  /* also set up corresponding SFT timestamps vector */
  LIGOTimeGPSVector *timestampsSFT;
  XLAL_CHECK ( (timestampsSFT = XLALCalloc (1, sizeof(*timestampsSFT)) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( (timestampsSFT->data = XLALCalloc (numSFTs, sizeof (*timestampsSFT->data) )) != NULL, XLAL_ENOMEM );
  timestampsSFT->length = numSFTs;

  /* prepare random-number generator */
  gsl_rng * r;
  XLAL_CHECK ( (r = gsl_rng_alloc (gsl_rng_taus2)) != NULL, XLAL_EFAILED );
  //time_t now = time(NULL);
  gsl_rng_set (r, 13);	/* sets random-number seed */

  /* ----- set up random-noise timeseries with gaps ---------- */
  for ( UINT4 alpha=0; alpha < numSFTs; alpha ++ )
    {
      /* record this SFT's timestamp */
      timestampsSFT->data[alpha] = epoch0;
      timestampsSFT->data[alpha].gpsSeconds += lround( alpha * Tsft );

      /* generate all data-points of this SFT */
      for ( UINT4 j=0; j < numR4SamplesPerSFT; j++ )
        {
          UINT4 alpha_j = alpha * numR4SamplesPerSFT + j;
          REAL8 ti = alpha * Tsft + j * dtR4;
          /* unit-variance Gaussian noise + sinusoid */
          TSdata[alpha_j] = gsl_ran_gaussian (r, sigmaN) + ampS * sin ( LAL_TWOPI * freqS * ti + phi0S); // ( alpha_j == 1000 ? 1 : 0 )
        } // for js < numR4SamplesPerSFT

    } /* for alpha < numSFTs */

  /* ----- generate SFTs from this timeseries ---------- */
  SFTParams XLAL_INIT_DECL(sftParams);
  sftParams.Tsft = Tsft;
  sftParams.timestamps = timestampsSFT;
  sftParams.noiseSFTs = NULL;
  sftParams.window = NULL;

  SFTVector *outSFTs;
  XLAL_CHECK ( (outSFTs = XLALSignalToSFTs ( outTS, &sftParams )) != NULL, XLAL_EFUNC );

  /* free memory */
  gsl_rng_free (r);
  XLALFree ( timestampsSFT->data );
  XLALFree ( timestampsSFT );

  /* return timeseries and SFTvector */
  (*ts)   = outTS;
  (*sfts) = outSFTs;

  return XLAL_SUCCESS;

} // XLALgenerateRandomData()


int
write_SFTdata (const char *fname, const SFTtype *sft)
{
  XLAL_CHECK ( fname != NULL, XLAL_EINVAL );
  XLAL_CHECK ( sft != NULL, XLAL_EINVAL );

  REAL8 f0 = sft->f0;
  REAL8 df = sft->deltaF;

  FILE *fp;
  XLAL_CHECK ( (fp = fopen(fname, "wb")) != NULL, XLAL_ESYS );
  for ( UINT4 k = 0; k < sft->data->length; k ++ )
    {
      REAL8 fk = f0 + k * df;
      fprintf ( fp, "%.9f      % 6e  % 6e  \n", fk, crealf(sft->data->data[k]), cimagf(sft->data->data[k]) );
    } // for k < numFreqBins

  fclose(fp);

  return XLAL_SUCCESS;

} // write_SFTdata()

SFTVector *
XLALDuplicateSFTVector ( const SFTVector *sftsIn )
{
  XLAL_CHECK_NULL ( (sftsIn != NULL) && ( sftsIn->length > 0), XLAL_EINVAL );

  UINT4 numSFTs = sftsIn->length;
  UINT4 numBins = sftsIn->data[0].data->length;

  SFTVector *sftsOut;
  XLAL_CHECK_NULL ( (sftsOut = XLALCreateSFTVector ( numSFTs, numBins )) != NULL, XLAL_EFUNC );

  for ( UINT4 alpha=0; alpha < numSFTs; alpha++ )
    {
      SFTtype *thisSFTIn = &sftsIn->data[alpha];
      SFTtype *thisSFTOut = &sftsOut->data[alpha];

      COMPLEX8Vector *tmp = thisSFTOut->data;
      memcpy ( thisSFTOut, thisSFTIn, sizeof(*thisSFTOut) );
      thisSFTOut->data = tmp;
      thisSFTOut->data->length = numBins;
      memcpy ( thisSFTOut->data->data, thisSFTIn->data->data, numBins * sizeof(thisSFTOut->data->data[0]) );

    } // for alpha < numSFTs

  return sftsOut;

} // XLALDuplicateSFTVector()

// compare two SFTs, return XLAL_SUCCESS if OK, error otherwise
int
XLALCompareSFTs ( const SFTtype *sft1, const SFTtype *sft2, const VectorComparison *tol )
{
  // check input sanity
  XLAL_CHECK ( sft1 != NULL, XLAL_EINVAL );
  XLAL_CHECK ( sft2 != NULL, XLAL_EINVAL );
  XLAL_CHECK ( tol  != NULL, XLAL_EINVAL );

  XLAL_CHECK ( XLALGPSCmp ( &(sft1->epoch), &(sft2->epoch) ) == 0, XLAL_ETOL );
  REAL8 tolFreq = 10 * LAL_REAL8_EPS;
  REAL8 err_f0 = sft1->f0 - sft2->f0;
  XLAL_CHECK ( err_f0 < tolFreq, XLAL_ETOL, "f0_1 = %.16g, f0_2 = %.16g, relErr_f0 = %g (tol = %g)\n", sft1->f0, sft2->f0, err_f0, tolFreq );
  REAL8 relErr_df = RELERR ( sft1->deltaF, sft2->deltaF );
  XLAL_CHECK ( relErr_df < tolFreq, XLAL_ETOL, "dFreq1 = %g, dFreq2 = %g, relErr_df = %g (tol = %g)", sft1->deltaF, sft2->deltaF, relErr_df, tolFreq );

  XLAL_CHECK ( XLALUnitCompare( &(sft1->sampleUnits), &(sft2->sampleUnits) ) == 0, XLAL_ETOL );

  XLAL_CHECK ( sft1->data->length == sft2->data->length, XLAL_ETOL, "sft1->length = %d, sft2->length = %d\n", sft1->data->length, sft2->data->length  );

  VectorComparison XLAL_INIT_DECL(cmp);
  XLAL_CHECK ( XLALCompareCOMPLEX8Vectors ( &cmp, sft1->data, sft2->data, tol ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;

} // XLALCompareSFTs()
