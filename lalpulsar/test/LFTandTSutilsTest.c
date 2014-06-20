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

/* ---------- local types ---------- */

/*---------- Global variables ----------*/
static LALUnit emptyLALUnit;

/* ----- User-variables: can be set from config-file or command-line */
/* ---------- local prototypes ---------- */
int main(void);


int test_XLALSFTVectorToLFT(void);
int test_XLALSincInterpolateCOMPLEX8TimeSeries(void);
int test_XLALDirichletInterpolateSFT ( void );

int XLALgenerateRandomData ( REAL4TimeSeries **ts, SFTVector **sfts );
int write_SFTdata (const char *fname, const SFTtype *sft);
SFTVector *XLALDuplicateSFTVector ( const SFTVector *sftsIn );
int XLALCompareSFTs ( const SFTtype *sft1, const SFTtype *sft2, const VectorComparison *tol );
static COMPLEX8 testSignal ( REAL8 t, REAL8 sigmaN );

/*----------------------------------------------------------------------*/
/* Main Function starts here */
/*----------------------------------------------------------------------*/
/**
 * MAIN function
 */
int main(void)
{

  XLAL_CHECK ( test_XLALSFTVectorToLFT() == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK ( test_XLALSincInterpolateCOMPLEX8TimeSeries() == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK ( test_XLALDirichletInterpolateSFT() == XLAL_SUCCESS, XLAL_EFUNC );

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
  SFTtype *lftR4 = NULL;
  XLAL_CHECK ( XLALExtractBandFromSFT ( &lftR4, lft0, out_fMin, out_Band ) == XLAL_SUCCESS, XLAL_EFUNC );
  for ( UINT4 k=0; k < lftR4->data->length; k ++ ) {
    lftR4->data->data[k] *= dt_R4;
  }

  // ----- 2) compute LFT directly from SFTs ----------
  SFTtype *lftSFTs0;
  XLAL_CHECK ( (lftSFTs0 = XLALSFTVectorToLFT ( sftsBand, 1 )) != NULL, XLAL_EFUNC );
  XLALDestroySFTVector ( sftsBand );

  // ----- re-extract frequency band of interest
  SFTtype *lftSFTs = NULL;
  XLAL_CHECK ( XLALExtractBandFromSFT ( &lftSFTs, lftSFTs0, out_fMin, out_Band ) == XLAL_SUCCESS, XLAL_EFUNC );

  if ( lalDebugLevel & LALINFO )
  {   // ----- write debug output
    XLAL_CHECK ( XLALdumpREAL4TimeSeries ( "TS_R4.dat", tsR4 ) == XLAL_SUCCESS, XLAL_EFUNC );

    XLAL_CHECK ( write_SFTdata ("LFT_R4T.dat", lftR4 ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK ( write_SFTdata ("LFT_SFTs.dat", lftSFTs ) == XLAL_SUCCESS, XLAL_EFUNC );

  } // end: debug output

  // ========== compare resulting LFTs ==========
  VectorComparison XLAL_INIT_DECL(tol);
  tol.relErr_L1 	= 4e-2;
  tol.relErr_L2		= 5e-2;
  tol.angleV 		= 5e-2;	// rad
  tol.relErr_atMaxAbsx 	= 1.2e-2;
  tol.relErr_atMaxAbsy  = 1.2e-2;

  XLALPrintInfo ("Comparing LFT from REAL4-timeseries, to LFT from heterodyned COMPLEX8-timeseries:\n");
  XLAL_CHECK ( XLALCompareSFTs ( lftR4, lftSFTs, &tol ) == XLAL_SUCCESS, XLAL_EFUNC );

  // ---------- free memory ----------
  XLALDestroyREAL4TimeSeries ( tsR4 );
  XLALDestroyREAL4FFTPlan ( planR4 );

  XLALDestroySFT ( lft0 );
  XLALDestroySFT ( lftR4 );
  XLALDestroySFT ( lftSFTs );
  XLALDestroySFT ( lftSFTs0 );

  return XLAL_SUCCESS;

} // test_XLALSFTVectorToLFT()


int
test_XLALSincInterpolateCOMPLEX8TimeSeries ( void )
{

  COMPLEX8TimeSeries* tsIn;
  REAL8 f0 = 100;	// heterodyning frequency
  REAL8 dt = 0.1;	// sampling frequency = 10Hz
  LIGOTimeGPS epoch = { 100, 0 };
  REAL8 tStart = XLALGPSGetREAL8 ( &epoch );
  UINT4 numSamples = 1000;
  REAL8 Tspan = numSamples * dt;

  XLAL_CHECK ( (tsIn = XLALCreateCOMPLEX8TimeSeries ( "test TS_in", &epoch, f0, dt, &emptyLALUnit, numSamples )) != NULL, XLAL_EFUNC );
  for ( UINT4 j = 0; j < numSamples; j ++ ) {
    tsIn->data->data[j] = testSignal ( tStart + j * dt, 0 );
  } // for j < numSamples

  // ---------- interpolate this onto new time-samples
  UINT4 Dterms = 16;
  REAL8 safety = (Dterms+1.0) * dt;	// avoid truncated interpolation to minimize errors, set to 0 for seeing boundary-effects [they're not so bad...]
  LIGOTimeGPS epochOut = epoch;
  XLALGPSAdd ( &epochOut, safety );
  REAL8 TspanOut = Tspan - 2 * safety;

  REAL8 dtOut = dt / 10;
  UINT4 numSamplesOut = lround ( TspanOut / dtOut );
  COMPLEX8TimeSeries *tsOut;
  XLAL_CHECK ( (tsOut = XLALCreateCOMPLEX8TimeSeries ( "test TS_out", &epochOut, f0, dtOut, &emptyLALUnit, numSamplesOut )) != NULL, XLAL_EFUNC );


  REAL8 tStartOut = XLALGPSGetREAL8 ( &epochOut );
  REAL8Vector *times_out;
  XLAL_CHECK ( (times_out = XLALCreateREAL8Vector ( numSamplesOut )) != NULL, XLAL_EFUNC );
  for ( UINT4 j = 0; j < numSamplesOut; j ++ )
    {
      REAL8 t_j = tStartOut + j * dtOut;
      times_out->data[j] = t_j;
    } // for j < numSamplesOut

  XLAL_CHECK ( XLALSincInterpolateCOMPLEX8TimeSeries ( &(tsOut->data), times_out, tsIn, Dterms ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLALDestroyREAL8Vector ( times_out );

  // ---------- check accuracy of interpolation
  COMPLEX8TimeSeries *tsFull;
  XLAL_CHECK ( (tsFull = XLALCreateCOMPLEX8TimeSeries ( "test TS_full", &epochOut, f0, dtOut, &emptyLALUnit, numSamplesOut )) != NULL, XLAL_EFUNC );
  for ( UINT4 j = 0; j < numSamplesOut; j ++ ) {
    tsFull->data->data[j] = testSignal ( tStartOut + j * dtOut, 0 );
  } // for j < numSamplesOut

  // ----- out debug info
  if ( lalDebugLevel & LALINFO )
    {
      XLAL_CHECK ( XLALdumpCOMPLEX8TimeSeries ( "TS_in.dat", tsIn ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK ( XLALdumpCOMPLEX8TimeSeries ( "TS_out.dat", tsOut ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK ( XLALdumpCOMPLEX8TimeSeries ( "TS_full.dat", tsFull ) == XLAL_SUCCESS, XLAL_EFUNC );
    } // if LALINFO

  VectorComparison XLAL_INIT_DECL(tol);
  tol.relErr_L1 	= 2e-2;
  tol.relErr_L2		= 2e-2;
  tol.angleV		= 2e-2;
  tol.relErr_atMaxAbsx	= 2e-2;
  tol.relErr_atMaxAbsy	= 2e-2;

  XLALPrintInfo ("Comparing sinc-interpolated timeseries to exact signal timeseries:\n");
  VectorComparison XLAL_INIT_DECL(cmp);
  XLAL_CHECK ( XLALCompareCOMPLEX8Vectors ( &cmp, tsOut->data, tsFull->data, &tol ) == XLAL_SUCCESS, XLAL_EFUNC );

  // ---------- free memory
  XLALDestroyCOMPLEX8TimeSeries ( tsIn );
  XLALDestroyCOMPLEX8TimeSeries ( tsOut );
  XLALDestroyCOMPLEX8TimeSeries ( tsFull );

  return XLAL_SUCCESS;

} // test_XLALSincInterpolateCOMPLEX8TimeSeries()

int
test_XLALDirichletInterpolateSFT ( void )
{
  REAL8 f0 = 0;		// heterodyning frequency
  REAL8 sigmaN = 0.001;
  REAL8 dt = 0.1;	// sampling frequency = 10Hz
  LIGOTimeGPS epoch = { 100, 0 };
  REAL8 tStart = XLALGPSGetREAL8 ( &epoch );
  UINT4 numSamples = 1000;
  REAL8 Tspan = numSamples * dt;
  REAL8 df = 1.0 / Tspan;

  UINT4 numSamples0padded = 3 * numSamples;
  REAL8 Tspan0padded = numSamples0padded * dt;
  REAL8 df0padded = 1.0 / Tspan0padded;

  UINT4 numBins = NhalfPosDC ( numSamples );
  UINT4 numBins0padded = NhalfPosDC ( numSamples0padded );
  // original timeseries
  REAL4TimeSeries* ts;
  XLAL_CHECK ( (ts = XLALCreateREAL4TimeSeries ( "test TS_in", &epoch, f0, dt, &emptyLALUnit, numSamples )) != NULL, XLAL_EFUNC );
  for ( UINT4 j = 0; j < numSamples; j ++ ) {
    ts->data->data[j] = crealf ( testSignal ( tStart + j * dt, sigmaN ) );
  } // for j < numSamples

  // zero-padded to double length
  REAL4TimeSeries* ts0padded;
  XLAL_CHECK ( (ts0padded = XLALCreateREAL4TimeSeries ( "test TS_padded", &epoch, f0, dt, &emptyLALUnit, numSamples0padded )) != NULL, XLAL_EFUNC );
  memcpy ( ts0padded->data->data, ts->data->data, numSamples * sizeof(ts0padded->data->data[0]) );
  memset ( ts0padded->data->data + numSamples, 0, (numSamples0padded - numSamples) * sizeof(ts0padded->data->data[0]) );

  // compute FFT on ts and ts0padded
  REAL4FFTPlan *plan, *plan0padded;
  XLAL_CHECK ( (plan        = XLALCreateForwardREAL4FFTPlan ( numSamples, 0 )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( (plan0padded = XLALCreateForwardREAL4FFTPlan ( numSamples0padded, 0 )) != NULL, XLAL_EFUNC );
  COMPLEX8Vector *fft, *fft0padded;
  XLAL_CHECK ( (fft        = XLALCreateCOMPLEX8Vector ( numBins )) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( (fft0padded = XLALCreateCOMPLEX8Vector ( numBins0padded )) != NULL, XLAL_ENOMEM );

  XLAL_CHECK ( XLALREAL4ForwardFFT ( fft, ts->data, plan ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALREAL4ForwardFFT ( fft0padded, ts0padded->data, plan0padded ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLALDestroyREAL4TimeSeries ( ts );
  XLALDestroyREAL4TimeSeries ( ts0padded );
  XLALDestroyREAL4FFTPlan( plan );
  XLALDestroyREAL4FFTPlan( plan0padded );

  SFTtype XLAL_INIT_DECL(tmp);
  tmp.f0 = f0;
  tmp.deltaF = df;
  tmp.data = fft;

  REAL8 Band = 0.5/dt - f0;
  SFTtype *sft = NULL;
  XLAL_CHECK ( XLALExtractBandFromSFT ( &sft, &tmp, f0, Band ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLALDestroyCOMPLEX8Vector ( fft );


  tmp.f0 = f0;
  tmp.deltaF = df0padded;
  tmp.data = fft0padded;
  SFTtype *sft0padded = NULL;
  XLAL_CHECK ( XLALExtractBandFromSFT ( &sft0padded, &tmp, f0, Band ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLALDestroyCOMPLEX8Vector ( fft0padded );

  // ---------- interpolate input SFT onto frequency bins of padded-ts FFT for comparison
  UINT4 Dterms = 16;
  REAL8 safetyBins = (Dterms + 1.0);	// avoid truncated interpolation to minimize errors, set to 0 for seeing boundary-effects [they're not so bad...]

  REAL8 fMin = f0 + safetyBins * df;
  fMin = round(fMin / df0padded) * df0padded;
  UINT4 numBinsOut = numBins0padded - 3 * safetyBins;
  REAL8 BandOut = (numBinsOut-1) * df0padded;

  SFTtype *sftUpsampled = NULL;
  XLAL_CHECK ( XLALExtractBandFromSFT ( &sftUpsampled, sft0padded, fMin + 0.5*df0padded, BandOut - df0padded ) == XLAL_SUCCESS, XLAL_EFUNC );

  SFTtype *sftInterpolated;
  XLAL_CHECK ( (sftInterpolated = XLALDirichletInterpolateSFT ( sft, fMin, df0padded, numBinsOut, Dterms )) != NULL, XLAL_EFUNC );

  // ----- out debug info
  if ( lalDebugLevel & LALINFO )
    {
      XLAL_CHECK ( write_SFTdata ( "SFT_in.dat", sft ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK ( write_SFTdata ( "SFT_in0padded.dat", sft0padded ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK ( write_SFTdata ( "SFT_upsampled.dat", sftUpsampled ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK ( write_SFTdata ( "SFT_interpolated.dat", sftInterpolated ) == XLAL_SUCCESS, XLAL_EFUNC );
    } // if LALINFO

  // ---------- check accuracy of interpolation
  VectorComparison XLAL_INIT_DECL(tol);
  tol.relErr_L1 	= 8e-2;
  tol.relErr_L2		= 8e-2;
  tol.angleV		= 8e-2;
  tol.relErr_atMaxAbsx	= 4e-3;
  tol.relErr_atMaxAbsy	= 4e-3;

  XLALPrintInfo ("Comparing Dirichlet SFT interpolation with upsampled SFT:\n");
  XLAL_CHECK ( XLALCompareSFTs ( sftUpsampled, sftInterpolated, &tol ) == XLAL_SUCCESS, XLAL_EFUNC );

  // ---------- free memory
  XLALDestroySFT ( sft );
  XLALDestroySFT ( sft0padded );
  XLALDestroySFT ( sftUpsampled );
  XLALDestroySFT ( sftInterpolated );

  return XLAL_SUCCESS;

} // test_XLALDirichletInterpolateSFT()


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
  REAL8 sigmaN = 0.1;

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
          TSdata[alpha_j] = crealf ( testSignal ( ti, sigmaN ) );
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
  XLALFree ( timestampsSFT->data );
  XLALFree ( timestampsSFT );

  /* return timeseries and SFTvector */
  (*ts)   = outTS;
  (*sfts) = outSFTs;

  return XLAL_SUCCESS;

} // XLALgenerateRandomData()

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


static COMPLEX8
testSignal ( REAL8 t, REAL8 sigmaN )
{
  static gsl_rng * r = NULL;
  static BOOLEAN firstCall = 1;
  if ( firstCall )
    {
      XLAL_CHECK ( (r = gsl_rng_alloc (gsl_rng_taus2)) != NULL, XLAL_EFAILED );
      gsl_rng_set (r, 13);	/* sets random-number seed */
      firstCall = 0;
    }

  REAL8 amp1 = 0.001;
  REAL8 amp2 = 1;
  REAL8 freq1 = LAL_PI;
  REAL8 freq2 = 2;

  COMPLEX8 y = amp1 * sin ( LAL_TWOPI * freq1 * t ) + I * amp2 * cos ( LAL_TWOPI * freq2 * t );
  y += gsl_ran_gaussian (r, sigmaN);
  y += I * gsl_ran_gaussian (r, sigmaN);

  return y;
} // testSignal()

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
