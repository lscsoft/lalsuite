/*
 *  Copyright (C) 2009 Chris Messenger, Reinhard Prix, Pinkesh Patel, Xavier Siemens, Holger Pletsch
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

#include "config.h"

/* System includes */
#include <stdio.h>

/* GSL includes */
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

/* LAL-includes */
#include <lal/LFTandTSutils.h>
#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/LogPrintf.h>
#include <lal/TimeSeries.h>
#include <lal/ComplexFFT.h>
#include <lal/ComputeFstat.h>
#include <lal/CWFastMath.h>

/*---------- DEFINES ----------*/
#define MYMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MYMIN(x,y) ( (x) < (y) ? (x) : (y) )
#define SQ(x) ((x)*(x))

#define NhalfPosDC(N) ((UINT4)(ceil ( ((N)/2.0 - 1e-6 ))))	/* round up */
#define NhalfNeg(N) ((UINT4)( (N) - NhalfPosDC(N) ))		/* round down (making sure N+ + N- = (N-1) */

#define GPS2REAL8(gps) (1.0 * (gps).gpsSeconds + 1.e-9 * (gps).gpsNanoSeconds )
#define RELERR(x,y) ( cabs( (x) - (y) ) / ( 0.5 * (cabs(x) + cabs(y)) ) )

/*---------- Global variables ----------*/
static LALUnit emptyLALUnit;

#define NUM_FACT 7
static const REAL8 inv_fact[NUM_FACT] = { 1.0, 1.0, (1.0/2.0), (1.0/6.0), (1.0/24.0), (1.0/120.0), (1.0/720.0) };

/* ---------- local prototypes ---------- */

/*---------- empty initializers ---------- */

/* ---------- function definitions ---------- */

/**
 * Turn the given multi-IFO SFTvectors into one long Fourier transform (LFT) over the total observation time
 */
SFTtype *
XLALSFTVectorToLFT ( const SFTVector *sfts,	/**< input SFT vector */
		     REAL8 upsampling )		/**< upsampling factor >= 1 */
{
  COMPLEX8FFTPlan *SFTplan, *LFTplan;

  COMPLEX8TimeSeries *lTS = NULL;		/* long time-series corresponding to full set of SFTs */
  COMPLEX8TimeSeries *sTS = NULL; 		/* short time-series corresponding to a single SFT */

  REAL8 fHet;				/* heterodyning frequency */
  LIGOTimeGPS epoch = {0,0};

  /* constant quantities for all SFTs */
  SFTtype *firstSFT;
  SFTtype *outputLFT;
  UINT4 numBinsSFT;
  REAL8 dfSFT;
  REAL8 Tsft;

  REAL8 deltaT;

  UINT4 numSFTs;
  UINT4 n;
  REAL8 Tspan;
  UINT4 numTimeSamples;
  UINT4 NnegSFT, NnegLFT;
  REAL8 f0SFT, f0LFT;
  UINT4 numSFTsFit;

  if ( !sfts || (sfts->length == 0) )
    {
      XLALPrintError ("%s: empty SFT input!\n", __func__ );
      XLAL_ERROR_NULL (XLAL_EINVAL);
    }
  if ( upsampling < 1 )
    {
      XLALPrintError ("%s: upsampling factor (%f) must be >= 1 \n", __func__, upsampling );
      XLAL_ERROR_NULL (XLAL_EINVAL);
    }

  /* some useful shorthands */
  numSFTs = sfts->length;
  firstSFT = &(sfts->data[0]);
  numBinsSFT = firstSFT->data->length;
  dfSFT = firstSFT->deltaF;
  Tsft = 1.0 / dfSFT;
  deltaT = 1.0 / ( numBinsSFT * dfSFT );

  f0SFT = firstSFT->f0;

  /* ---------- determine time-span of the final long time-series */
  Tspan = XLALGPSDiff ( &sfts->data[numSFTs-1].epoch, &sfts->data[0].epoch ) + Tsft;
  Tspan *= upsampling;

  /* NOTE: Tspan MUST be an integer multiple of Tsft,
   * in order for the frequency bins of the final FFT
   * to be commensurate with the SFT bins.
   * This is required so that fHet is an exact
   * frequency-bin in both cases
   */
  numSFTsFit = lround ( Tspan / Tsft );
  Tspan = numSFTsFit * Tsft;
  numTimeSamples = lround ( Tspan / deltaT );

  /* determine the heterodyning frequency */
  /* fHet = DC of our internal DFTs */
  NnegSFT = NhalfNeg ( numBinsSFT );
  fHet = f0SFT + 1.0 * NnegSFT * dfSFT;

  /* translate this back into fmin for the LFT (counting down from DC==fHet) */
  NnegLFT = NhalfNeg ( numTimeSamples );
  f0LFT = fHet - NnegLFT / Tspan;

  XLALPrintInfo ("NSFT = %d, NLFT = %d, fminSFT = %.9f, fHet = %.9f, fminLFT = %.9f\n",
	  numBinsSFT, numTimeSamples, f0SFT, fHet, f0LFT );


  /* ----- Prepare invFFT: compute plan for FFTW */
  if ( (SFTplan = XLALCreateReverseCOMPLEX8FFTPlan( numBinsSFT, 0 )) == NULL )
    {
      XLALPrintError ( "%s: XLALCreateReverseCOMPLEX8FFTPlan(%d, ESTIMATE ) failed! errno = %d!\n", __func__, numBinsSFT, xlalErrno );
      XLAL_ERROR_NULL ( XLAL_EFUNC );
    }

  /* ----- prepare long TimeSeries container ---------- */
  if ( (lTS = XLALCreateCOMPLEX8TimeSeries ( firstSFT->name, &firstSFT->epoch, 0, deltaT, &emptyLALUnit, numTimeSamples )) == NULL )
    {
      XLALPrintError ("%s: XLALCreateCOMPLEX8TimeSeries() for %d timesteps failed! errno = %d!\n", __func__, numTimeSamples, xlalErrno );
      XLAL_ERROR_NULL ( XLAL_EFUNC );
    }
  memset ( lTS->data->data, 0, numTimeSamples * sizeof(*lTS->data->data)); /* set all time-samples to zero */

  /* ----- Prepare short time-series holding ONE invFFT of a single SFT */
  if ( (sTS = XLALCreateCOMPLEX8TimeSeries ( "short timeseries", &epoch, 0, deltaT, &emptyLALUnit, numBinsSFT )) == NULL )
    {
      XLALPrintError ( "%s: XLALCreateCOMPLEX8TimeSeries() for %d timesteps failed! errno = %d!\n", __func__, numBinsSFT, xlalErrno );
      XLAL_ERROR_NULL ( XLAL_EFUNC );
    }

  /* ----- prepare output LFT ---------- */
  if ( (outputLFT = LALCalloc ( 1, sizeof(*outputLFT) )) == NULL )
    {
      XLALPrintError ( "%s: LALCalloc ( 1, %d ) failed!\n", __func__, sizeof(*outputLFT) );
      XLAL_ERROR_NULL ( XLAL_ENOMEM );
    }

  /* prepare LFT header */
  strcpy ( outputLFT->name, firstSFT->name );
  strcat ( outputLFT->name, ":long Fourier transform");
  outputLFT->epoch  = firstSFT->epoch;
  outputLFT->f0     = f0LFT;
  outputLFT->deltaF = 1.0 / Tspan;
  outputLFT->sampleUnits = firstSFT->sampleUnits;

  if ( (outputLFT->data = XLALCreateCOMPLEX8Vector ( numTimeSamples )) == NULL )
    {
      XLALPrintError ( "%s: XLALCreateCOMPLEX8Vector(%d) failed! xlalErrno = %d\n", __func__, numTimeSamples, xlalErrno );
      XLAL_ERROR_NULL ( XLAL_ENOMEM );
    }

  /* ---------- loop over all SFTs and inverse-FFT them ---------- */
  for ( n = 0; n < numSFTs; n ++ )
    {
      SFTtype *thisSFT = &(sfts->data[n]);
      UINT4 bin0_n;
      REAL8 offset_n, nudge_n;
      UINT4 copyLen, binsLeft;

      UINT4 k;
      REAL8 offset0, offsetEff, hetCycles;
      REAL4 hetCorr_re, hetCorr_im;
      REAL4 fact_re, fact_im;
      REAL4 norm = 1.0 / numBinsSFT;


      /* find bin in long timeseries corresponding to starttime of *this* SFT */
      offset_n = XLALGPSDiff ( &thisSFT->epoch, &firstSFT->epoch );
      bin0_n = lround ( offset_n / deltaT );	/* round to closest bin */

      nudge_n = bin0_n * deltaT - offset_n;		/* rounding error */
      nudge_n = 1e-9 * round( nudge_n * 1e9 );	/* round to closest nanosecond */
      {
	REAL8 t0 = XLALGPSGetREAL8 ( &firstSFT->epoch );
	XLALPrintInfo ("n = %d: t0_n = %f, sft_tn =(%d,%d), bin-offset = %g s, corresponding to %g timesteps\n",
                       n, t0 + bin0_n * deltaT, sfts->data[n].epoch.gpsSeconds,  sfts->data[n].epoch.gpsNanoSeconds, nudge_n, nudge_n/deltaT );
      }

      if ( (nudge_n != 0) && (XLALTimeShiftSFT ( thisSFT, nudge_n ) != XLAL_SUCCESS) )
	{
	  XLALPrintError ( "%s: XLALTimeShiftSFT(sft-%d, %g) failed! errno = %d!\n", __func__, n, nudge_n, xlalErrno );
	  XLAL_ERROR_NULL ( XLAL_EFUNC );
	}


      /* determine heterodyning phase-correction for this SFT */
      offset0 = XLALGPSDiff ( &thisSFT->epoch, &firstSFT->epoch );
      /* we know that fHet * Tsft = integer, so we only need the remainder of offset_t0 % Tsft */
      offsetEff = fmod ( offset0, Tsft );
      offsetEff = 1e-9 * round ( offsetEff * 1e9 );	/* round to closest integer multiple of nanoseconds */

      hetCycles = fmod ( fHet * offsetEff, 1);	/* required heterodyning phase-correction for this SFT */

      XLAL_CHECK_NULL( XLALSinCos2PiLUT (&hetCorr_im, &hetCorr_re, -hetCycles ) == XLAL_SUCCESS, XLAL_EFUNC );

      /* apply all phase- and normalizing factors */
      fact_re = norm * hetCorr_re;
      fact_im = norm * hetCorr_im;

      XLALPrintInfo ("SFT n = %d: (tn - t0) = %g s EQUIV %g s, hetCycles = %g, ==> fact = %g + i %g\n", n, offset0, offsetEff, hetCycles, fact_re, fact_im );

      for ( k = 0; k < numBinsSFT; k ++ )
	{
	  REAL8 binReal, binImag;

	  binReal = fact_re * crealf(thisSFT->data->data[k]) - fact_im * cimagf(thisSFT->data->data[k]);
	  binImag = fact_re * cimagf(thisSFT->data->data[k]) + fact_im * crealf(thisSFT->data->data[k]);

	  thisSFT->data->data[k] = crectf( binReal, binImag );
	} /* k < numBins */

      if ( XLALReorderSFTtoFFTW (thisSFT->data) != XLAL_SUCCESS )
	{
	  XLALPrintError ( "%s: XLALReorderSFTtoFFTW() failed! errno = %d!\n", __func__, xlalErrno );
	  XLAL_ERROR_NULL ( XLAL_EFUNC );
	}

      if ( XLALCOMPLEX8VectorFFT( sTS->data, thisSFT->data, SFTplan ) != XLAL_SUCCESS )
	{
	  XLALPrintError ( "%s: XLALCOMPLEX8VectorFFT() failed! errno = %d!\n", __func__, xlalErrno );
	  XLAL_ERROR_NULL ( XLAL_EFUNC );
	}

      /* copy short (shifted) timeseries into correct location within long timeseries */
      binsLeft = numTimeSamples - bin0_n - 1;
      copyLen = MYMIN ( numBinsSFT, binsLeft );		/* make sure not to write past the end of the long TS */
      memcpy ( &lTS->data->data[bin0_n], sTS->data->data, copyLen * sizeof(lTS->data->data[0]) );

    } /* for n < numSFTs */

  /* ---------- now FFT the long timeseries ---------- */

  /* ----- compute plan for FFTW */
  if ( (LFTplan = XLALCreateForwardCOMPLEX8FFTPlan( numTimeSamples, 0 )) == NULL )
    {
      XLALPrintError ( "%s: XLALCreateForwardCOMPLEX8FFTPlan(%d, ESTIMATE ) failed! errno = %d!\n", __func__, numTimeSamples, xlalErrno );
      XLAL_ERROR_NULL ( XLAL_EFUNC );
    }

  if ( XLALCOMPLEX8VectorFFT( outputLFT->data, lTS->data, LFTplan ) != XLAL_SUCCESS )
    {
      XLALPrintError ( "%s: XLALCOMPLEX8VectorFFT() failed! errno = %d!\n", __func__, xlalErrno );
      XLAL_ERROR_NULL ( XLAL_EFUNC );
    }

  if ( XLALReorderFFTWtoSFT (outputLFT->data) != XLAL_SUCCESS )
    {
      XLALPrintError ( "%s: XLALReorderFFTWtoSFT() failed! errno = %d!\n", __func__, xlalErrno );
      XLAL_ERROR_NULL ( XLAL_EFUNC );
    }

  /* cleanup memory */
  XLALDestroyCOMPLEX8TimeSeries ( sTS );
  XLALDestroyCOMPLEX8TimeSeries ( lTS );
  XLALDestroyCOMPLEX8FFTPlan ( SFTplan );
  XLALDestroyCOMPLEX8FFTPlan ( LFTplan );

  return outputLFT;

} /* XLALSFTVectorToLFT() */


/**
 * Turn the given SFTvector into one long time-series, properly dealing with gaps.
 *
 * NOTE: this function <b>modifies</b> the input SFTs in the process!
 * If you need to reuse the SFTvector afterwards, you need to copy it before
 * passing it into this function.
 *
 */
COMPLEX8TimeSeries *
XLALSFTVectorToCOMPLEX8TimeSeries ( SFTVector *sfts,                /**< [in/out] SFT vector, gets modified! */
				    const LIGOTimeGPS *start_in,    /**< [in] start time */
				    const LIGOTimeGPS *end_in       /**< [in] input end time */
				    )
{
  COMPLEX8FFTPlan *SFTplan;

  COMPLEX8TimeSeries *lTS = NULL;		/* long time-series corresponding to full set of SFTs */
  COMPLEX8TimeSeries *sTS = NULL; 		/* short time-series corresponding to a single SFT */

  REAL8 fHet;				/* heterodyning frequency */
  LIGOTimeGPS epoch = {0,0};
  LIGOTimeGPS start;
  LIGOTimeGPS end;

  /* constant quantities for all SFTs */
  SFTtype *firstSFT, *lastSFT;
  UINT4 numBinsSFT;
  REAL8 dfSFT;
  REAL8 Tsft;
  REAL8 deltaT;
  UINT4 numSFTs;
  UINT4 n;
  REAL8 Tspan;
  UINT4 numTimeSamples;
  UINT4 NnegSFT;
  REAL8 f0SFT;
  REAL8 SFTFreqBand;

  /* check sanity of input */
  if ( !sfts || (sfts->length == 0) )
    {
      XLALPrintError ("%s: empty SFT input!\n", __func__ );
      XLAL_ERROR_NULL (XLAL_EINVAL);
    }

  /* define some useful shorthands */
  numSFTs = sfts->length;
  firstSFT = &(sfts->data[0]);
  lastSFT = &(sfts->data[numSFTs-1]);
  numBinsSFT = firstSFT->data->length;
  dfSFT = firstSFT->deltaF;
  Tsft = 1.0 / dfSFT;
  SFTFreqBand = numBinsSFT * dfSFT;
  deltaT = 1.0 / SFTFreqBand;	/* we'll put DC into the middle of [f0, f0+Band], so sampling at fSamp=Band is sufficient */
  f0SFT = firstSFT->f0;

  /* if the start and end input pointers are NOT NULL then determine start and time-span of the final long time-series */
  if (start_in && end_in)
    {
      start.gpsSeconds = start_in->gpsSeconds;
      start.gpsNanoSeconds = start_in->gpsNanoSeconds;
      end.gpsSeconds = end_in->gpsSeconds;
      end.gpsNanoSeconds = end_in->gpsNanoSeconds;

      /* do sanity checks */
      if ( (XLALGPSDiff ( &end, &firstSFT->epoch ) ) < 0 )
	{
	  XLALPrintError ("%s: end time before first SFT!\n", __func__ );
	  XLAL_ERROR_NULL (XLAL_EINVAL);
	}
      if ( (XLALGPSDiff ( &start, &sfts->data[numSFTs-1].epoch) ) > Tsft )
	{
	  XLALPrintError ("%s: start time after end of data!\n", __func__ );
	  XLAL_ERROR_NULL (XLAL_EINVAL);
	}
    }
  else {   /* otherwise we use the start and end of the sft vector */
    start.gpsSeconds = firstSFT->epoch.gpsSeconds;
    start.gpsNanoSeconds = firstSFT->epoch.gpsNanoSeconds;
    end.gpsSeconds = lastSFT->epoch.gpsSeconds;
    end.gpsNanoSeconds = lastSFT->epoch.gpsNanoSeconds;
    if ( XLALGPSAdd(&end,Tsft) == NULL )
    {
      XLALPrintError ("%s: NULL pointer returned from XLALGPSAdd()!\n", __func__ );
      XLAL_ERROR_NULL (XLAL_EFAULT);
    }

  }

  /* determine output time span */
  if ( (Tspan = XLALGPSDiff ( &end, &start ) ) < 0 )
    {
      XLALPrintError ("%s: start time after end time!\n", __func__ );
      XLAL_ERROR_NULL (XLAL_EINVAL);
    }

  numTimeSamples = lround ( Tspan / deltaT );

  /* determine the heterodyning frequency */
  /* fHet = DC of our internal DFTs */
  NnegSFT = NhalfNeg ( numBinsSFT );
  fHet = f0SFT + 1.0 * NnegSFT * dfSFT;

  /* ----- Prepare invFFT of SFTs: compute plan for FFTW */
  if ( (SFTplan = XLALCreateReverseCOMPLEX8FFTPlan( numBinsSFT, 0 )) == NULL )
    {
      XLALPrintError ( "%s: XLALCreateReverseCOMPLEX8FFTPlan(%d, ESTIMATE ) failed! errno = %d!\n", __func__, numBinsSFT, xlalErrno );
      goto failed;
    }

  /* ----- Prepare short time-series holding ONE invFFT of a single SFT */
  if ( (sTS = XLALCreateCOMPLEX8TimeSeries ( "short timeseries", &epoch, 0, deltaT, &emptyLALUnit, numBinsSFT )) == NULL )
    {
      XLALPrintError ( "%s: XLALCreateCOMPLEX8TimeSeries() for %d timesteps failed! errno = %d!\n", __func__, numBinsSFT, xlalErrno );
      goto failed;
    }

  /* ----- prepare long TimeSeries container ---------- */
  if ( (lTS = XLALCreateCOMPLEX8TimeSeries ( firstSFT->name, &start, fHet, deltaT, &emptyLALUnit, numTimeSamples )) == NULL )
    {
      XLALPrintError ("%s: XLALCreateCOMPLEX8TimeSeries() for %d timesteps failed! errno = %d!\n", __func__, numTimeSamples, xlalErrno );
      goto failed;
    }
  memset ( lTS->data->data, 0, numTimeSamples * sizeof(*lTS->data->data)); 	/* set all time-samples to zero (in case there are gaps) */


  /* ---------- loop over all SFTs and inverse-FFT them ---------- */
  for ( n = 0; n < numSFTs; n ++ )
    {
      SFTtype *thisSFT = &(sfts->data[n]);
      UINT4 bin0_n;
      REAL8 offset_n, nudge_n;
      UINT4 copyLen, binsLeft;

      REAL8 offset0, offsetEff, hetCycles;
      COMPLEX8 hetCorrection;

      /* find bin in long timeseries corresponding to starttime of *this* SFT */
      offset_n = XLALGPSDiff ( &thisSFT->epoch, &start );
      bin0_n = lround ( offset_n / deltaT );	/* round to closest bin */

      nudge_n = bin0_n * deltaT - offset_n;		/* rounding error */
      nudge_n = 1e-9 * round ( nudge_n * 1e9 );	/* round to closest nanosecond */
      {
	REAL8 t0 = XLALGPSGetREAL8 ( &start );
	XLALPrintInfo ("n = %d: t0_n = %f, sft_tn =(%d,%d), bin-offset = %g s, corresponding to %g timesteps\n",
		n, t0 + bin0_n * deltaT, sfts->data[n].epoch.gpsSeconds,  sfts->data[n].epoch.gpsNanoSeconds, nudge_n, nudge_n/deltaT );
      }

      /* nudge SFT into integer timestep bin if necessary */
      if ( nudge_n != 0 )
	{
	  if  ( XLALTimeShiftSFT ( thisSFT, nudge_n ) != XLAL_SUCCESS )
	    {
	      XLALPrintError ( "%s: XLALTimeShiftSFT(sft-%d, %g) failed! errno = %d!\n", __func__, n, nudge_n, xlalErrno );
	      goto failed;
	    }
	}

      /* determine heterodyning phase-correction for this SFT */
      offset0 = XLALGPSDiff ( &thisSFT->epoch, &start );

      /* fHet * Tsft is an integer, because fHet is a frequency-bin of the input SFTs, so we only need the remainder offset_t0 % Tsft */
      offsetEff = fmod ( offset0, Tsft );
      offsetEff = 1e-9 * round ( offsetEff * 1e9 );	/* round to closest integer multiple of nanoseconds */
      hetCycles = fmod ( fHet * offsetEff, 1);			/* required heterodyning phase-correction for this SFT */

      {
        REAL4 hetCorrection_re, hetCorrection_im;
        XLAL_CHECK_NULL( XLALSinCos2PiLUT( &hetCorrection_im, &hetCorrection_re, -hetCycles ) == XLAL_SUCCESS, XLAL_EFUNC );
        hetCorrection = crectf( hetCorrection_re, hetCorrection_im );
      }

      /* Note: we also bundle the overall normalization of 'df' into the het-correction.
       * This ensures that the resulting timeseries will have the correct normalization, according to
       * x_l = invFT[sft]_l = df * sum_{k=0}^{N-1} xt_k * e^(i 2pi k l / N )
       * where x_l is the l-th timestamp, and xt_k is the k-th frequency bin of the SFT.
       * See the LAL-conventions on FFTs:  http://www.ligo.caltech.edu/docs/T/T010095-00.pdf
       * (the FFTw convention does not contain the factor of 'df', which is why we need to
       * apply it ourselves)
       *
       */
      hetCorrection *= ((REAL4) dfSFT);

      /* FIXME: check how time-critical this step is, using proper profiling! */
      if ( XLALMultiplySFTbyCOMPLEX8 ( thisSFT, hetCorrection ) != XLAL_SUCCESS )
	{
	  XLALPrintError ( "%s: XLALMultiplySFTbyCOMPLEX8(sft-%d) failed! errno = %d!\n", __func__, n, xlalErrno );
	  goto failed;
	}

      XLALPrintInfo ("SFT n = %d: (tn - t0) = %g s EQUIV %g s, hetCycles = %g, ==> fact = %g + i %g\n",
                     n, offset0, offsetEff, hetCycles, crealf(hetCorrection), cimagf(hetCorrection) );

      /* FIXME: check if required */
      if ( XLALReorderSFTtoFFTW (thisSFT->data) != XLAL_SUCCESS )
	{
	  XLALPrintError ( "%s: XLALReorderSFTtoFFTW() failed! errno = %d!\n", __func__, xlalErrno );
	  goto failed;
	}

      if ( XLALCOMPLEX8VectorFFT( sTS->data, thisSFT->data, SFTplan ) != XLAL_SUCCESS )
	{
	  XLALPrintError ( "%s: XLALCOMPLEX8VectorFFT() failed! errno = %d!\n", __func__, xlalErrno );
	  goto failed;
	}

      /* copy short (shifted) timeseries into correct location within long timeseries */
      binsLeft = numTimeSamples - bin0_n - 1;
      copyLen = MYMIN ( numBinsSFT, binsLeft );		/* make sure not to write past the end of the long TS */
      memcpy ( &lTS->data->data[bin0_n], sTS->data->data, copyLen * sizeof(lTS->data->data[0]) );

    } /* for n < numSFTs */

  goto success;

 failed:
  /* cleanup memory */
  XLALDestroyCOMPLEX8TimeSeries ( sTS );	/* short invSFT timeseries */
  XLALDestroyCOMPLEX8FFTPlan ( SFTplan );
  XLALDestroyCOMPLEX8TimeSeries ( lTS );

  XLAL_ERROR_NULL ( XLAL_EFUNC );

 success:
  /* cleanup memory */
  XLALDestroyCOMPLEX8TimeSeries ( sTS );	/* short invSFT timeseries */
  XLALDestroyCOMPLEX8FFTPlan ( SFTplan );

  return lTS;

} /* XLALSFTVectorToCOMPLEX8TimeSeries() */


/**
 * Change frequency-bin order from fftw-convention to a 'SFT'
 * ie. from FFTW: f[0], f[1],...f[N/2], f[-(N-1)/2], ... f[-2], f[-1]
 * to: f[-(N-1)/2], ... f[-1], f[0], f[1], .... f[N/2]
 */
int
XLALReorderFFTWtoSFT (COMPLEX8Vector *X)
{
  UINT4 N, Npos_and_DC, Nneg;

  /* temporary storage for data */
  COMPLEX8 *tmp;


  if ( !X || (X->length==0) )
    {
      XLALPrintError ("%s: empty input vector 'X'!\n", __func__ );
      XLAL_ERROR (XLAL_EINVAL);
    }

  N = X -> length;
  Npos_and_DC = NhalfPosDC ( N );
  Nneg = NhalfNeg ( N );

  /* allocate temporary storage for swap */
  if ( (tmp = XLALMalloc ( N * sizeof(*tmp) )) == NULL )
    {
      XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", __func__, N * sizeof(*tmp));
      XLAL_ERROR (XLAL_ENOMEM);
    }

  memcpy ( tmp, X->data, N * sizeof(*tmp) );

  /* Copy second half from FFTW: 'negative' frequencies */
  memcpy ( X->data, tmp + Npos_and_DC, Nneg * sizeof(*tmp) );

  /* Copy first half from FFTW: 'DC + positive' frequencies */
  memcpy ( X->data + Nneg, tmp, Npos_and_DC * sizeof(*tmp) );

  XLALFree(tmp);

  return XLAL_SUCCESS;

}  /* XLALReorderFFTWtoSFT() */

/**
 * Change frequency-bin order from 'SFT' to fftw-convention
 * ie. from f[-(N-1)/2], ... f[-1], f[0], f[1], .... f[N/2]
 * to FFTW: f[0], f[1],...f[N/2], f[-(N-1)/2], ... f[-2], f[-1]
 */
int
XLALReorderSFTtoFFTW (COMPLEX8Vector *X)
{
  UINT4 N, Npos_and_DC, Nneg;

  /* temporary storage for data */
  COMPLEX8 *tmp;


  if ( !X || (X->length==0) )
    {
      XLALPrintError ("%s: empty input vector 'X'!\n", __func__ );
      XLAL_ERROR (XLAL_EINVAL);
    }

  N = X->length;
  Npos_and_DC = NhalfPosDC ( N );
  Nneg = NhalfNeg ( N );

  /* allocate temporary storage for swap */
  if ( (tmp = XLALMalloc ( N * sizeof(*tmp) )) == NULL )
    {
      XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", __func__, N * sizeof(*tmp));
      XLAL_ERROR (XLAL_ENOMEM);
    }

  memcpy ( tmp, X->data, N * sizeof(*tmp) );

  /* Copy second half of FFTW: 'negative' frequencies */
  memcpy ( X->data + Npos_and_DC, tmp, Nneg * sizeof(*tmp) );

  /* Copy first half of FFTW: 'DC + positive' frequencies */
  memcpy ( X->data, tmp + Nneg, Npos_and_DC * sizeof(*tmp) );

  XLALFree(tmp);

  return XLAL_SUCCESS;

}  /* XLALReorderSFTtoFFTW() */


/**
 * Multiply SFT frequency bins by given complex factor.
 *
 * NOTE: this <b>modifies</b> the given SFT in place
 */
int
XLALMultiplySFTbyCOMPLEX8 ( SFTtype *sft,	/**< [in/out] SFT */
			    COMPLEX8 factor )	/**< [in] complex8 factor to multiply SFT by */
{
  UINT4 k;

  if ( !sft || !sft->data )
    {
      XLALPrintError ("%s: empty input SFT!\n", __func__ );
      XLAL_ERROR (XLAL_EINVAL);
    }

  for ( k=0; k < sft->data->length; k++)
    {
      REAL4 yRe, yIm;

      yRe = crealf(factor) * crealf(sft->data->data[k]) - cimagf(factor) * cimagf(sft->data->data[k]);
      yIm = crealf(factor) * cimagf(sft->data->data[k]) + cimagf(factor) * crealf(sft->data->data[k]);

      sft->data->data[k] = crectf( yRe, yIm );

    } /* for k < numBins */

  return XLAL_SUCCESS;

} /* XLALMultiplySFTbyCOMPLEX8() */


/**
 * Time-shift the given SFT by an amount of 'shift' seconds,
 * using the frequency-domain expression y(f) = x(f) * e^(-i 2pi f tau),
 * which shifts x(t) into y(t) = x(t - tau)
 *
 * NOTE: this <b>modifies</b> the SFT in place
 */
int
XLALTimeShiftSFT ( SFTtype *sft,	/**< [in/out] SFT to time-shift */
		   REAL8 shift )	/**< time-shift in seconds */
{
  UINT4 k;

  if ( !sft || !sft->data )
    {
      XLALPrintError ("%s: empty input SFT!\n", __func__ );
      XLAL_ERROR (XLAL_EINVAL);
    }

  for ( k=0; k < sft->data->length; k++)
    {
      REAL8 fk = sft->f0 + k * sft->deltaF;	/* frequency of k-th bin */
      REAL8 shiftCyles = shift * fk;
      REAL4 fact_re, fact_im;			/* complex phase-shift factor e^(-2pi f tau) */
      REAL4 yRe, yIm;

      XLAL_CHECK( XLALSinCos2PiLUT ( &fact_im, &fact_re, shiftCyles ) == XLAL_SUCCESS, XLAL_EFUNC );

      yRe = fact_re * crealf(sft->data->data[k]) - fact_im * cimagf(sft->data->data[k]);
      yIm = fact_re * cimagf(sft->data->data[k]) + fact_im * crealf(sft->data->data[k]);

      sft->data->data[k] = crectf( yRe, yIm );

    } /* for k < numBins */

  /* adjust SFTs epoch to the shift */
  XLALGPSAdd( &sft->epoch, shift );

  return XLAL_SUCCESS;

} /* XLALTimeShiftSFT() */


/**
 * Turn the given multiSFTvector into multiple long COMPLEX8TimeSeries, properly dealing with gaps.
 * Memory allocation for the output MultiCOMPLEX8TimeSeries is done within this function.
 *
 * NOTE : We enforce that each detectors timeseries has <b>equal</b> start times and time spans.
 * Also, the input MultiSFTs get <b>modified</b> in place.
 */
MultiCOMPLEX8TimeSeries *XLALMultiSFTVectorToCOMPLEX8TimeSeries (
                     MultiSFTVector *multisfts  /**< [in/out] multi SFT vector, gets modified! */					       )
{
  UINT4 i;
  MultiCOMPLEX8TimeSeries *out = NULL;		/* long time-series corresponding to full set of SFTs */
  LIGOTimeGPS start,end;

  /* check sanity of input */
  if ( !multisfts || (multisfts->length == 0) ) {
    XLALPrintError ("%s: empty multi SFT input!\n", __func__ );
    XLAL_ERROR_NULL (XLAL_EINVAL);
  }
  for (i=0;i<multisfts->length;i++)
    {
      if ( !multisfts->data[i] || (multisfts->data[i]->length == 0) ) {
	XLALPrintError ("%s: empty multiSFT->data[%d] input!\n", __func__,i );
	XLAL_ERROR_NULL (XLAL_EINVAL);
      }
    }

  /* determine the start and end times of the multiSFT observation */
  if ( XLALEarliestMultiSFTsample( &start, multisfts) != XLAL_SUCCESS ) {
    XLALPrintError("%s: Failed to run XLALEarliestMultiSFTsample()\n", __func__ );
    XLAL_ERROR_NULL (XLAL_EFAULT );
  }
  if ( XLALLatestMultiSFTsample( &end, multisfts) != XLAL_SUCCESS ) {
    XLALPrintError("%s: Failed to run XLALLatestMultiSFTsample()\n", __func__ );
    XLAL_ERROR_NULL (XLAL_EFAULT );
  }
  /*
  printf("earliest sample at %d %d\n",start.gpsSeconds,start.gpsNanoSeconds);
  printf("latest sample at %d %d\n",end.gpsSeconds,end.gpsNanoSeconds);
  */

  /* check that earliest is before latest */
  if ( (XLALGPSDiff ( &end, &start ) ) < 0 )
    {
      XLALPrintError ("%s: start time after end time!\n", __func__ );
      XLAL_ERROR_NULL (XLAL_EINVAL);
    }

  /* allocate memory for the output structure */
  if ((out = (MultiCOMPLEX8TimeSeries*)XLALMalloc(sizeof(MultiCOMPLEX8TimeSeries))) == NULL) {
    XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", __func__, sizeof(MultiCOMPLEX8TimeSeries));
    XLAL_ERROR_NULL (XLAL_ENOMEM);
  }
  out->length = multisfts->length;
  if ((out->data = XLALMalloc(out->length*sizeof(COMPLEX8TimeSeries*))) == NULL) {
    XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", __func__, out->length*sizeof(COMPLEX8TimeSeries*));
    XLAL_ERROR_NULL (XLAL_ENOMEM);
  }

  /* loop over detectors */
  for (i=0;i<multisfts->length;i++) {

    /* call XLALSFTVectorToCOMPLEX8TimeSeries for each detector */
    if ((out->data[i] = XLALSFTVectorToCOMPLEX8TimeSeries(multisfts->data[i],&start,&end)) == NULL) {
      XLALPrintError ("%s: Failed to run XLALSFTVectorToCOMPLEX8TimeSeries()\n", __func__);
      XLAL_ERROR_NULL (XLAL_EFAULT );
    }

  }

  return out;

} /* XLALMultiSFTVectorToCOMPLEX8TimeSeries() */


/**
 * Computed the weighted timeseries Fa(t) = x(t).a(t) and Fb(t) = x(t).b(t) for a multi-detector timeseries
 */
int XLALAntennaWeightCOMPLEX8TimeSeries (
            COMPLEX8TimeSeries **Faoft,                   /**< [out] the timeseries weighted by a(t) */
					  COMPLEX8TimeSeries **Fboft,                   /**< [out] the timeseries weighted by b(t) */
					  const COMPLEX8TimeSeries *timeseries,         /**< [in] the input timeseries */
					  const AMCoeffs *AMcoef,                       /**< [in] the AM coefficients */
					  const SFTVector *sfts                         /**< [in] the SFT data */
					  )
{
  UINT4 j,k;

  /* do sanity checks */


  /* local copies */
  REAL8 start = GPS2REAL8(timeseries->epoch);
  REAL8 fHet = timeseries->f0;
  REAL8 deltaT = timeseries->deltaT;
  UINT4 numTimeSamples = timeseries->data->length;
  REAL8 dfSFT = sfts->data[0].deltaF;
  REAL8 Tsft = 1.0 / dfSFT;
  UINT4 nbins = (UINT4)floor(0.5 + Tsft/deltaT);

  /* create empty timeseries structures for Fa(t) and Fb(t) */
  if ( ((*Faoft) = XLALCreateCOMPLEX8TimeSeries ( sfts->data[0].name, &(timeseries->epoch), fHet, deltaT, &emptyLALUnit, numTimeSamples )) == NULL ) {
    XLALPrintError ("%s: XLALCreateCOMPLEX8TimeSeries() for %d timesteps failed! errno = %d!\n", __func__, numTimeSamples, xlalErrno );
    XLAL_ERROR (XLAL_ENOMEM);
  }
  if ( ((*Fboft) = XLALCreateCOMPLEX8TimeSeries ( sfts->data[0].name,&(timeseries->epoch) , fHet, deltaT, &emptyLALUnit, numTimeSamples )) == NULL ) {
    XLALPrintError ("%s: XLALCreateCOMPLEX8TimeSeries() for %d timesteps failed! errno = %d!\n", __func__, numTimeSamples, xlalErrno );
    XLAL_ERROR (XLAL_ENOMEM);
  }
  memset ( (*Faoft)->data->data, 0, numTimeSamples * sizeof(*(*Faoft)->data->data)); 	/* set all time-samples to zero (in case there are gaps) */
  memset ( (*Fboft)->data->data, 0, numTimeSamples * sizeof(*(*Fboft)->data->data)); 	/* set all time-samples to zero (in case there are gaps) */


  /* loop over SFTs */
  for (j=0;j<sfts->length;j++) {

    REAL8 t = GPS2REAL8(sfts->data[j].epoch);                         /* the GPS time at the start of the SFT */
    UINT4 start_index = (UINT4)floor(0.5 + (t - start)/deltaT);                     /* index of timesample corresponding to the start of the SFT */
    REAL8 a = (REAL8)AMcoef->a->data[j];                              /* value of the antenna pattern a(t) at the MID-POINT of the SFT */
    REAL8 b = (REAL8)AMcoef->b->data[j];                              /* value of the antenna pattern b(t) at the MID-POINT of the SFT */

    /* loop over samples from this SFT */
    for (k=0;k<nbins;k++) {

      UINT4 time_index = start_index + k;

      /* weight the complex timeseries by the antenna patterns */
      (*Faoft)->data->data[time_index] = (((REAL4) a) * timeseries->data->data[time_index]);
      (*Fboft)->data->data[time_index] = (((REAL4) b) * timeseries->data->data[time_index]);

      }

    }

  /* success */
  return XLAL_SUCCESS;

}


/**
 * Computed the weighted timeseries Fa(t) = x(t).a(t) and Fb(t) = x(t).b(t) for a multi-detector timeseries
 */
int XLALAntennaWeightMultiCOMPLEX8TimeSeries (
                 MultiCOMPLEX8TimeSeries **Faoft,                        /**< [out] the timeseries weighted by a(t) */
					       MultiCOMPLEX8TimeSeries **Fboft,                        /**< [out] the timeseries weighted by b(t) */
					       const MultiCOMPLEX8TimeSeries *multiTimeseries,         /**< [in] the input multi-detector timeseries */
					       const MultiAMCoeffs *multiAMcoef,                       /**< [in] the multi-detector AM coefficients */
					       const MultiSFTVector *multisfts                         /**< [in] the multi-detector SFT data */
					       )
{
  UINT4 i;

  /* do sanity checks */
  if ( !multiTimeseries || (multiTimeseries->length == 0) ) {
    XLALPrintError ("%s: empty multiTimeseries input!\n", __func__ );
    XLAL_ERROR (XLAL_EINVAL);
  }
  if ( multiTimeseries->length != multisfts->length ) {
    XLALPrintError ("%s: incorrect length of multiTimeseries or multisfts input!\n", __func__ );
    XLAL_ERROR (XLAL_EINVAL);
  }

  /* allocate memory for the output structure */
  if ( ((*Faoft) = XLALMalloc(sizeof(MultiCOMPLEX8TimeSeries))) == NULL ) {
    XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", __func__, sizeof(MultiCOMPLEX8TimeSeries));
    XLAL_ERROR (XLAL_ENOMEM);
  }
  if ( ((*Fboft) = XLALMalloc(sizeof(MultiCOMPLEX8TimeSeries))) == NULL ) {
    XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", __func__, sizeof(MultiCOMPLEX8TimeSeries));
    XLAL_ERROR (XLAL_ENOMEM);
  }

  (*Faoft)->length = multisfts->length;
  (*Fboft)->length = multisfts->length;

  if (((*Faoft)->data = XLALMalloc((multisfts->length)*sizeof(COMPLEX8TimeSeries*))) == NULL) {
    XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", __func__, (multisfts->length)*sizeof(COMPLEX8TimeSeries*));
    XLAL_ERROR (XLAL_ENOMEM);
  }
  if (((*Fboft)->data = XLALMalloc((multisfts->length)*sizeof(COMPLEX8TimeSeries*))) == NULL) {
    XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", __func__, (multisfts->length)*sizeof(COMPLEX8TimeSeries*));
    XLAL_ERROR (XLAL_ENOMEM);
  }


  /* loop over detectors */
  for (i=0;i<multisfts->length;i++) {

    /* point to current detector params */
    COMPLEX8TimeSeries *timeseries = multiTimeseries->data[i];
    AMCoeffs *AMcoef = multiAMcoef->data[i];
    SFTVector *SFTs = multisfts->data[i];

    if ( XLALAntennaWeightCOMPLEX8TimeSeries(&((*Faoft)->data[i]),&((*Fboft)->data[i]),timeseries,AMcoef,SFTs) != XLAL_SUCCESS ) {
      XLALPrintError("\nXLALAntennaWeightMultiCOMPLEX8TimeSeries() failed with error = %d\n\n", xlalErrno );
      XLAL_ERROR (XLAL_EFAULT);
    }

  }

  /* success */
  return XLAL_SUCCESS;

}


/**
 * Performs barycentric resampling on a multi-detector timeseries
 */
int XLALBarycentricResampleMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries **Faoft_RS,                         /**< [out] the resampled timeseries Fa(t_SSB) */
						     MultiCOMPLEX8TimeSeries **Fboft_RS,                         /**< [out] the resampled timeseries Fb(t_SSB) */
						     const MultiCOMPLEX8TimeSeries *Faoft,                       /**< [in] the detector frame timeseries Fa(t) */
						     const MultiCOMPLEX8TimeSeries *Fboft,                       /**< [in] the detector frame timeseries Fb(t) */
						     const MultiSSBtimes *multiSSB,                              /**< [in] the multi-detector SSB times data */
						     const MultiSFTVector *multiSFTs,                            /**< [in] the multi-detector SFT data */
						     const REAL8 deltaF                                          /**< [in] the user defined output frequency resolution */
						     )
{
  UINT4 i;
  LIGOTimeGPS earliest,latest;
  UINT4 numTimeSamples;
  UINT4 numDetectors;
  REAL8 Tspan;
  REAL8 deltaT;
  REAL8 Teff;
  REAL8 fHet;
  REAL8 Tsft;

  /* do sanity checks on input */
  if ( !Faoft || (Faoft->length == 0) ) {
    XLALPrintError ("%s: empty Faoft input!\n", __func__ );
    XLAL_ERROR (XLAL_EINVAL);
  }
  numDetectors = Faoft->length;

  if ( !Fboft || (Fboft->length != numDetectors ) ) {
    XLALPrintError ("%s: empty or invalid length of Fboft input!\n", __func__ );
    XLAL_ERROR (XLAL_EINVAL);
  }

  /* check that SFTs exist */
  if ( !multiSFTs || (multiSFTs->length != numDetectors) ) {
    XLALPrintError ("%s: empty or invalid length of SFTs input!\n", __func__ );
    XLAL_ERROR (XLAL_EINVAL);
  }

  /* check if the SSB times vector has the correct number of entries - one for each SFT */
  if ( !multiSSB || (multiSSB->length != numDetectors ) ) {
    XLALPrintError ("%s: empty or incorrect length of multiSSB input!\n", __func__ );
    XLAL_ERROR (XLAL_EINVAL);
  }

  /* define the length of an SFT (assuming 1/T resolution) */
  Tsft = 1.0 / multiSFTs->data[0]->data[0].deltaF;

  /* find earliest SSB time */
  if ( (XLALEarliestMultiSSBtime (&earliest,multiSSB,Tsft)) != XLAL_SUCCESS ) {
    XLALPrintError("\nXLALEarliestMultiSSBtime() failed with error = %d\n\n", xlalErrno );
    XLAL_ERROR (XLAL_EFAULT);
  }

  /* find latest SSB time */
  if ( (XLALLatestMultiSSBtime (&latest,multiSSB,Tsft)) != XLAL_SUCCESS ) {
    XLALPrintError("\nXLALLatestMultiSSBtime() failed with error = %d\n\n", xlalErrno );
    XLAL_ERROR (XLAL_EFAULT);
  }

  /* determine resampled timeseries parameters */
  Tsft = 1.0 / multiSFTs->data[0]->data[0].deltaF;         /* define the length of an SFT (assuming 1/T resolution) */
  Tspan = XLALGPSDiff(&latest,&earliest);                  /* the time span from the earliest sample to the end of the latest sample over *all* detectors */
  deltaT = Faoft->data[0]->deltaT;                         /* the sample rate of the downsampled detector frame timeseries */
  Teff = 1.0/deltaF;                                       /* the effective observation time based on the requested frequency resolution (for zero padding) */
  fHet = Faoft->data[0]->f0;                               /* the input timeseries heterodyne frequency */

  /* if the requested frequency resolution gives an effective observation time less than the actual data span then we exit */
  if (Tspan > Teff) {
    XLALPrintError ("%s: requested frequency resolution too coarse for resampling!\n", __func__ );
    XLAL_ERROR (XLAL_EINVAL);
  }

  /* redefine sample rate and compute number of samples in the new timeseries */
  numTimeSamples = (UINT4)ceil(Teff/deltaT);      /* we use ceil() so that we artificially widen the band rather than reduce it */
  deltaT = Teff/numTimeSamples;

  /* allocate memory for the output resampled timeseries for Fa and Fb */
  if (((*Faoft_RS) = XLALMalloc(sizeof(MultiCOMPLEX8TimeSeries))) == NULL) {
    XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", __func__, sizeof(MultiCOMPLEX8TimeSeries));
    XLAL_ERROR (XLAL_ENOMEM);
  }
  if (((*Fboft_RS) = XLALMalloc(sizeof(MultiCOMPLEX8TimeSeries))) == NULL) {
    XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", __func__, sizeof(MultiCOMPLEX8TimeSeries));
    XLAL_ERROR (XLAL_ENOMEM);
  }
  (*Faoft_RS)->length = multiSSB->length;
  (*Fboft_RS)->length = multiSSB->length;

  if (((*Faoft_RS)->data = XLALMalloc(multiSSB->length*sizeof(COMPLEX8TimeSeries*))) == NULL) {
    XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", __func__, multiSSB->length*sizeof(COMPLEX8TimeSeries*));
    XLAL_ERROR (XLAL_ENOMEM);
  }
  if (((*Fboft_RS)->data = XLALMalloc(multiSSB->length*sizeof(COMPLEX8TimeSeries*))) == NULL) {
    XLALPrintError ("%s: Failed to allocate XLALMalloc(%d)\n", __func__, multiSSB->length*sizeof(COMPLEX8TimeSeries*));
    XLAL_ERROR (XLAL_ENOMEM);
  }

  /* loop over detectors */
  for (i=0;i<multiSSB->length;i++) {

    /* shorthand pointers */
    SSBtimes *SSB = multiSSB->data[i];
    SFTVector *SFTs = multiSFTs->data[i];
    COMPLEX8TimeSeries *Fa = Faoft->data[i];
    COMPLEX8TimeSeries *Fb = Fboft->data[i];

    /* create empty timeseries structures for the resampled Fa(t) and Fb(t) */
    if ( ((*Faoft_RS)->data[i] = XLALCreateCOMPLEX8TimeSeries ( Faoft->data[i]->name, &earliest, fHet, deltaT, &emptyLALUnit, numTimeSamples )) == NULL ) {
      XLALPrintError ("%s: XLALCreateCOMPLEX8TimeSeries() for %d timesteps failed! errno = %d!\n", __func__, numTimeSamples, xlalErrno );
      XLAL_ERROR (XLAL_ENOMEM);
    }
    if ( ((*Fboft_RS)->data[i] = XLALCreateCOMPLEX8TimeSeries ( Fboft->data[i]->name, &earliest , fHet, deltaT, &emptyLALUnit, numTimeSamples )) == NULL ) {
      XLALPrintError ("%s: XLALCreateCOMPLEX8TimeSeries() for %d timesteps failed! errno = %d!\n", __func__, numTimeSamples, xlalErrno );
      XLAL_ERROR (XLAL_ENOMEM);
    }
    memset ( (*Faoft_RS)->data[i]->data->data, 0, numTimeSamples * sizeof(*(*Faoft_RS)->data[i]->data->data)); 	/* set all time-samples to zero (in case there are gaps) */
    memset ( (*Fboft_RS)->data[i]->data->data, 0, numTimeSamples * sizeof(*(*Fboft_RS)->data[i]->data->data)); 	/* set all time-samples to zero (in case there are gaps) */

    /* perform resampling on current detector timeseries */
    if ( XLALBarycentricResampleCOMPLEX8TimeSeries(&((*Faoft_RS)->data[i]),&((*Fboft_RS)->data[i]),Fa,Fb,SSB,SFTs) != XLAL_SUCCESS ) {
      XLALPrintError("\nXLALBarycentricResampleCOMPLEX8TimeSeries() failed with error = %d\n\n", xlalErrno );
      XLAL_ERROR (XLAL_EFAULT);
    }

  }

  /* success */
  return XLAL_SUCCESS;

}


/**
 * Performs barycentric resampling on a COMPLEX8TimeSeries
 * We expect that the output timeseries has already been allocated correctly.
 *
 */
int XLALBarycentricResampleCOMPLEX8TimeSeries ( COMPLEX8TimeSeries **Faoft_RS,                         /**< [out] the resampled timeseries Fa(t_SSB) */
						COMPLEX8TimeSeries **Fboft_RS,                         /**< [out] the resampled timeseries Fb(t_SSB) */
						const COMPLEX8TimeSeries *Faoft,                       /**< [in] the input detector frame timeseries Fa(t) */
						const COMPLEX8TimeSeries *Fboft,                       /**< [in] the input detector frame timeseries Fb(t) */
						const SSBtimes *SSB,                                   /**< [in] the SSB times at the midpoints of each SFT */
						const SFTVector *SFTs                                  /**< [in] the SFT data */
						)
{
  UINT4 i,j,k;
  REAL8Vector *detectortimes = NULL;                /* used to store a vector of *non-uniform* time values in the detector frame (used for interpolation) */
  UINT4 FAFB_LENGTH = 4;                          /* the total number of elements in an instance of Fa and Fb */
  REAL8Vector* FaFb[FAFB_LENGTH];                 /* used to store Fa and Fb real and imaginary parts in 4 seperate vectors */
  REAL8Vector *t_DET = NULL;                        /* used to store a vector of *uniform* time values in the detector frame (used for interpolation) */
  gsl_spline* spline_FaFb[FAFB_LENGTH];           /* used to store spline coefficients for Fa and Fb real and imaginary parts in 4 seperate vectors */
  UINT4 numTimeSamples_DET;
  UINT4 numSFTs;
  REAL8 Tsft;
  REAL8 refTime;
  REAL8 start_DET;
  REAL8 fHet;
  REAL8 deltaT_DET;
  REAL8 start_SSB;
  REAL8 deltaT_SSB;

  /* do sanity checks on input */
  if ( !Faoft || (Faoft->data->length == 0) ) {
    XLALPrintError ("%s: empty Faoft input!\n", __func__ );
    XLAL_ERROR (XLAL_EINVAL);
  }
  numTimeSamples_DET = Faoft->data->length;         /* set the number of time samples in the detector frame as that defined in Fa */

  /* check if the Fa and Fb timeseries have equal lengths */
  if ( !Fboft || (Fboft->data->length != numTimeSamples_DET ) ) {
    XLALPrintError ("%s: empty or incorrect length of Fboft input!\n", __func__ );
    XLAL_ERROR (XLAL_EINVAL);
  }

  /* check that the Fa and Fb epochs are equal */
  if ( (XLALGPSCmp(&Faoft->epoch, &Fboft->epoch) != 0 ) ) {
    XLALPrintError ("%s: Faoft and Fboft epochs do not match!\n", __func__ );
    XLAL_ERROR (XLAL_EINVAL);
  }

  /* check that the output vectors are not NULL */
  if ( !(*Faoft_RS) || ((*Faoft_RS)->data->length == 0) || !(*Fboft_RS) || ((*Fboft_RS)->data->length == 0) ) {
    XLALPrintError ("%s: empty output vectors Faoft_RS and/or Fboft_RS!\n", __func__ );
    XLAL_ERROR (XLAL_EINVAL);
  }

  /* check that SFTs exist */
  if ( !SFTs || (SFTs->length == 0) ) {
    XLALPrintError ("%s: empty SFTs input!\n", __func__ );
    XLAL_ERROR (XLAL_EINVAL);
  }
  numSFTs = SFTs->length;                           /* set the number of SFTs */

  /* check if the SSB times vector has the correct number of entries - one for each SFT */
  if ( !SSB || (SSB->DeltaT->length != numSFTs ) || ( SSB->Tdot->length != numSFTs ) ) {
    XLALPrintError ("%s: empty or incorrect length of SSB input!\n", __func__ );
    XLAL_ERROR (XLAL_EINVAL);
  }

  /* define some useful shorthands */
  Tsft = 1.0 / SFTs->data[0].deltaF;                  /* define the length of an SFT (assuming 1/T resolution) */
  refTime = GPS2REAL8(SSB->refTime);                   /* set the reftime at which the doppler parameters are defined */
  start_DET = GPS2REAL8(Faoft->epoch);                 /* set the start time of the timeseries at the detector */
  fHet = Faoft->f0;                                    /* set the heterodyne frequency of the input timeseries */
  deltaT_DET = Faoft->deltaT;                          /* set the sampling time at the detector */
  start_SSB = GPS2REAL8((*Faoft_RS)->epoch);           /* set the start time of the resampled output SSB timeseries */
  deltaT_SSB = (*Faoft_RS)->deltaT;                    /* set the sampling time of the resampled output SSB timeseries */
  REAL8 end_DET = start_DET + (numTimeSamples_DET-1) * deltaT_DET;	/* time of last sample in detector timeseries */

  /* allocate memory for the uniformly sampled detector time samples (Fa and Fb real and imaginary) */
  for (i=0;i<FAFB_LENGTH;i++) {
    if ( (FaFb[i] = XLALCreateREAL8Vector(numTimeSamples_DET)) == NULL ) {
      XLALPrintError ("%s: XLALCreateREAL8Vector() failed!  errno = %d!\n", __func__, xlalErrno );
      XLAL_ERROR (XLAL_ENOMEM);
    }
  }
  /* allocate memory for the *uniform* detector time vector required for interpolation */
  if ( (t_DET = XLALCreateREAL8Vector(numTimeSamples_DET)) == NULL ) {
    XLALPrintError ("%s: XLALCreateREAL8Vector() failed!  errno = %d!\n", __func__, xlalErrno );
    XLAL_ERROR (XLAL_ENOMEM);
  }

  /* place the timeseries into REAL8Vectors for gsl to be able to interpolate them */
  /* this is annoying because the data is currently stored in a COMPLEX8 format so we can't just point to it */
  for (j=0;j<numTimeSamples_DET;j++) {
    t_DET->data[j] = start_DET + j*deltaT_DET;              /* fill in the uniform detector time vector */
    FaFb[0]->data[j] = crealf(Faoft->data->data[j]);
    FaFb[1]->data[j] = cimagf(Faoft->data->data[j]);
    FaFb[2]->data[j] = crealf(Fboft->data->data[j]);
    FaFb[3]->data[j] = cimagf(Fboft->data->data[j]);
  }

  /* initialise the gsl spline interpolation for each of the 4 timeseries */
  for (i=0;i<FAFB_LENGTH;i++) {
    if ( (XLALGSLInitInterpolateREAL8Vector(&(spline_FaFb[i]),t_DET,FaFb[i])) != XLAL_SUCCESS ) {
      XLALPrintError("\nXLALGSLInitInterpolateREAL8Vector() failed with error = %d\n\n", xlalErrno );
      XLAL_ERROR (XLAL_EFAULT);
    }
  }

  /* loop over SFTs to compute the detector frame time samples corresponding to uniformly sampled SSB time samples */
  for (j=0;j<numSFTs;j++) {

    REAL8Vector *out_FaFb[4];                                                              /* output vectors for the real and imaginary parts of the resampled Fa and Fb */

    /* define some useful shorthands */
    REAL8 Tdot = SSB->Tdot->data[j];                                                       /* the instantaneous time derivitive dt_SSB/dt_DET at the MID-POINT of the SFT */
    REAL8 SFTmid_SSB = refTime + SSB->DeltaT->data[j];                                     /* MID-POINT time of the SFT at the SSB */
    REAL8 SFTstart_SSB = SFTmid_SSB - 0.5*Tsft*Tdot;                                       /* START time of the SFT at the SSB */
    REAL8 SFTend_SSB = SFTmid_SSB + 0.5*Tsft*Tdot;                                         /* END time of the SFT at the SSB */
    REAL8 SFTstart_DET = GPS2REAL8(SFTs->data[j].epoch);                                   /* START time of the SFT at the detector */
    REAL8 SFTmid_DET = SFTstart_DET + 0.5*Tsft;                                            /* MID-POINT time of the SFT at the detector */

    /* define some indices */
    UINT4 idx_start_SSB = floor(0.5 + (SFTstart_SSB - start_SSB)/deltaT_SSB);              /* the index of the resampled timeseries corresponding to the start of the SFT */
    UINT4 idx_end_SSB = floor(0.5 + (SFTend_SSB - start_SSB)/deltaT_SSB);                  /* the index of the resampled timeseries corresponding to the end of the SFT */
    UINT4 numSamples_SSB = idx_end_SSB - idx_start_SSB + 1;                                /* the number of samples in the SSB for this SFT */

    /* allocate memory for the *non-uniform* detector time samples for this SFT */
    /* have to allocate it inside the loop because it may have different lengths for each SFT */
    if ( (detectortimes = XLALCreateREAL8Vector(numSamples_SSB)) == NULL ) {
      XLALPrintError ("%s: XLALCreateREAL8Vector() failed!  errno = %d!\n", __func__, xlalErrno );
      XLAL_ERROR (XLAL_ENOMEM);
    }

    /* for each time sample in the SSB frame for this SFT we estimate the detector time. */
    /* We use a linear approximation expanding around the midpoint of an SFT where  */
    /* t_DET = SFTmid_DET + (t_SSB - SFTmid_SSB)*dt_DET/dt_SSB */
    for (k=0;k<numSamples_SSB;k++) {

      REAL8 t_SSB = start_SSB + (k+idx_start_SSB)*deltaT_SSB;                                  /* the SSB time of the current resampled time sample */
      detectortimes->data[k] = SFTmid_DET + (t_SSB - SFTmid_SSB)/Tdot;                         /* the approximated DET time of the current resampled time sample */

      /*
       * NOTE: we need to be careful that none of the times falls outside
       * of the range of detector timesamples, in order to avoid problems in the interpolation
       * therefore we truncate the detector-times to fully fall within the detector timeseries span
       */
      if ( detectortimes->data[k] > end_DET )
        {
          detectortimes->data[k] = end_DET;
          XLALPrintWarning ("%s: time-sample jSFT=%d, kSample=%d at t=%f to interpolate is *after* detector-timeseries, nudged back to end (end=%f)\n",
                            __func__, j, k, detectortimes->data[k], end_DET );
        }
      if ( detectortimes->data[k] < start_DET )
        {
          detectortimes->data[k] = start_DET;
          XLALPrintWarning ("%s: time-sample jSFT=%d, kSample=%d at t=%f to interpolate is *before* detector-timeseries, nudged to beginning (start=%f)\n",
                            __func__, j, k, detectortimes->data[k], start_DET );
        }

    } /* for k < numSamples_SSB */

    /* interpolate on the non-uniformly sampled detector time vector for this SFT for re and im parts of Fa and Fb */
    /* this function allocates memory for the output vectors */
    for (i=0;i<FAFB_LENGTH;i++) {
      if ( XLALGSLInterpolateREAL8Vector(&(out_FaFb[i]),detectortimes,spline_FaFb[i]) != XLAL_SUCCESS ) {
	XLALPrintError("\nXLALInterpolateMultiCOMPLEX8TimeSeries() failed with error = %d\n\n", xlalErrno );
	XLAL_ERROR (XLAL_EFAULT);
      }
    }

    /* place these interpolated timeseries into the output */
    /* and apply correction due to non-zero heterodyne frequency of input */
    for (k=0;k<numSamples_SSB;k++) {

      UINT4 idx = k + idx_start_SSB;                                                                     /* the full resampled timeseries index */
      UINT4 numSamples_RS = (*Faoft_RS)->data->length;
      if ( idx >= numSamples_RS ) {	// temporary FIX to avoid writing outside of memory bounds (FIXME!)
        break;
      }
      REAL8 tDiff = start_SSB + idx*deltaT_SSB - detectortimes->data[k];                                 /* the difference between t_SSB and t_DET */
      REAL8 cycles = fmod ( fHet*tDiff, 1 );                                                            /* the accumulated heterodyne cycles */
      REAL4 cosphase,sinphase;                                                                           /* the real and imaginary parts of the phase correction */

      /* use a look-up-table for speed to compute real and imaginary phase */
      XLAL_CHECK( XLALSinCos2PiLUT ( &sinphase, &cosphase, -cycles ) == XLAL_SUCCESS, XLAL_EFUNC );

      /* printf("j = %d t = %6.12f tb = %6.12f tDiff = %6.12f\n",j,detectortimes->data[k],start_SSB + idx*deltaT_SSB,tDiff); */
      (*Faoft_RS)->data->data[idx] = crectf( out_FaFb[0]->data[k]*cosphase - out_FaFb[1]->data[k]*sinphase, out_FaFb[1]->data[k]*cosphase + out_FaFb[0]->data[k]*sinphase );
      (*Fboft_RS)->data->data[idx] = crectf( out_FaFb[2]->data[k]*cosphase - out_FaFb[3]->data[k]*sinphase, out_FaFb[3]->data[k]*cosphase + out_FaFb[2]->data[k]*sinphase );

    }

    /* free memory used for this SFT */
    for (i=0;i<FAFB_LENGTH;i++) XLALDestroyREAL8Vector(out_FaFb[i]);
    XLALDestroyREAL8Vector(detectortimes);

  } /* end loop over SFTs */

  /* free memory */
  for (i=0;i<4;i++) {
    XLALDestroyREAL8Vector(FaFb[i]);
    gsl_spline_free(spline_FaFb[i]);
  }
  XLALDestroyREAL8Vector(t_DET);

  /* success */
  return XLAL_SUCCESS;

}


/**
 * Find the latest timestamp in a multi-SSB data structure
 */
int XLALGSLInterpolateREAL8Vector( REAL8Vector **yi,
				   REAL8Vector *xi,
				   gsl_spline *spline
				   )
 {
   /* check input */

   UINT4 i;
   UINT4 numSamples = xi->length;
   gsl_interp_accel *acc = gsl_interp_accel_alloc();

   /* allocate memory for output vector */
   if ( ((*yi) = XLALCreateREAL8Vector(numSamples)) == NULL ) {
     XLALPrintError ("%s: XLALCreateREAL8Vector() failed!  errno = %d!\n", __func__, xlalErrno );
     XLAL_ERROR (XLAL_ENOMEM);
   }

   /* perform inerpolation */
   for (i=0;i<numSamples;i++) {
     (*yi)->data[i] = gsl_spline_eval(spline,xi->data[i],acc);
   }

   /* free memory */
   gsl_interp_accel_free(acc);

   return XLAL_SUCCESS;

}


/**
 * Find the latest timestamp in a multi-SSB data structure
 */
int XLALGSLInitInterpolateREAL8Vector( gsl_spline **spline,
				       REAL8Vector *x,
				       REAL8Vector *y
				       )
 {
   /* check input */

   UINT4 numSamples_in = x->length;
   REAL8 *xtemp = x->data;
   REAL8 *ytemp = y->data;

   /* compute spline interpolation coefficients */
   if ( ((*spline) = gsl_spline_alloc(gsl_interp_cspline,numSamples_in)) == NULL ) {
     XLALPrintError ("%s: XLALInitGSLInterpolateREAL8Vector() failed!  errno = %d!\n", __func__, xlalErrno );
     XLAL_ERROR (XLAL_ENOMEM);
   }
   gsl_spline_init((*spline),xtemp,ytemp,numSamples_in);

   return XLAL_SUCCESS;

}


/**
 * Shifts an FFT output vector such that the Niquest frequency is the central bin
 */
int XLALFFTShiftCOMPLEX8Vector(COMPLEX8Vector **x)
{
  UINT4 N = (*x)->length;
  UINT4 NQ = NhalfPosDC(N);
  UINT4 NminusNQ = N - NQ;
 /*  printf("NQ = %d NminusNQ = %d N = %d\n",NQ,NminusNQ,N); */

  /* allocate temp memory */
  COMPLEX8 *temp = XLALMalloc(NQ*sizeof(COMPLEX8));

  /* copy the bins 0 -> NQ - 1 to temp memory */
  if ( memcpy(temp,(*x)->data,NQ*sizeof(COMPLEX8)) == NULL ) {
    XLALPrintError ("%s: memcpy() failed!  errno = %d!\n", __func__, xlalErrno );
    XLAL_ERROR (XLAL_ENOMEM);
  }

  /* copy the bins NQ -> N - 1 to the start */
  if (memcpy((*x)->data,&((*x)->data[NQ]),NminusNQ*sizeof(COMPLEX8)) == NULL ) {
    XLALPrintError ("%s: memcpy() failed!  errno = %d!\n", __func__, xlalErrno );
    XLAL_ERROR (XLAL_ENOMEM);
  }

  /* copy the temp bins to the end */
  if (memcpy(&((*x)->data[NminusNQ]),temp,NQ*sizeof(COMPLEX8)) == NULL ) {
    XLALPrintError ("%s: memcpy() failed!  errno = %d!\n", __func__, xlalErrno );
    XLAL_ERROR (XLAL_ENOMEM);
  }

  /* free temp memory */
  XLALFree(temp);

  return XLAL_SUCCESS;

}


/**
 * Multi-detector wrapper for XLALFrequencyShiftCOMPLEX8TimeSeries
 * NOTE: this <b>modifies</b> the MultiCOMPLEX8Timeseries in place
 */
int
XLALFrequencyShiftMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries **x,	/**< [in/out] timeseries to time-shift */
					    const REAL8 shift )	                /**< [in] freq-shift in Hz */
{
  UINT4 i;

  if ( !(*x) || (*x)->length == 0 )
    {
      XLALPrintError ("%s: empty input COMPLEX8timeSeries!\n", __func__ );
      XLAL_ERROR (XLAL_EINVAL);
    }

  /* loop over detectors */
  for ( i=0; i < (*x)->length; i++)
    {

      /* shift the frequency of each detector's data */
      if ( XLALFrequencyShiftCOMPLEX8TimeSeries ( &((*x)->data[i]), shift) != XLAL_SUCCESS ) {
	XLALPrintError ("%s: memcpy() failed!  errno = %d!\n", __func__, xlalErrno );
	XLAL_ERROR (XLAL_EFAULT);
      }

    }

  return XLAL_SUCCESS;

} /* XLALFrequencyShiftMultiCOMPLEX8TimeSeries() */


/**
 * Freq-shift the given COMPLEX8Timeseries by an amount of 'shift' Hz,
 * using the time-domain expression y(t) = x(t) * e^(-i 2pi df t),
 * which shifts x(f) into y(f) = x(f - df)
 *
 * NOTE: this <b>modifies</b> the COMPLEX8TimeSeries in place
 */
int
XLALFrequencyShiftCOMPLEX8TimeSeries ( COMPLEX8TimeSeries **x,	        /**< [in/out] timeseries to time-shift */
				       const REAL8 shift )	        /**< [in] freq-shift in Hz */
{
  UINT4 k;
  REAL8 deltat;

  if ( !(*x) || !(*x)->data )
    {
      XLALPrintError ("%s: empty input COMPLEX8TimeSeries!\n", __func__ );
      XLAL_ERROR (XLAL_EINVAL);
    }

  /* get timeseries epoch */
  deltat = (*x)->deltaT;

  /* loop over COMPLEX8TimeSeries elements */
  for ( k=0; k < (*x)->data->length; k++)
    {
      REAL8 tk = k * deltat;	/* time of k-th bin */
      REAL8 shiftCycles = shift * tk;
      REAL4 fact_re, fact_im;			/* complex phase-shift factor e^(-2pi f tau) */
      REAL4 yRe, yIm;

      /* use a sin/cos look-up-table for speed */
      XLAL_CHECK( XLALSinCos2PiLUT ( &fact_im, &fact_re, shiftCycles ) == XLAL_SUCCESS, XLAL_EFUNC );

      /* apply the phase shift */
      yRe = fact_re * crealf((*x)->data->data[k]) - fact_im * cimagf((*x)->data->data[k]);
      yIm = fact_re * cimagf((*x)->data->data[k]) + fact_im * crealf((*x)->data->data[k]);
      (*x)->data->data[k] = crectf( yRe, yIm );

    } /* for k < numBins */

  /* adjust timeseries heterodyne frequency to the shift */
  (*x)->f0 -= shift;

  return XLAL_SUCCESS;

} /* XLALFrequencyShiftCOMPLEX8TimeSeries() */


/**
 * Apply a spin-down correction to the Fa and Fb complex timeseries
 * using the time-domain expression y(t) = x(t) * e^(-i 2pi sum f_k * (t-tref)^(k+1)),
 *
 * NOTE: this <b>modifies</b> the COMPLEX8TimeSeries Fa and Fb in place
 */
int
XLALSpinDownCorrectionMultiFaFb ( MultiCOMPLEX8TimeSeries **Fa,	                /**< [in/out] timeseries to time-shift */
				  MultiCOMPLEX8TimeSeries **Fb,	                /**< [in/out] timeseries to time-shift */
				  const PulsarDopplerParams *doppler		/**< parameter-space point to correct for */
				  )
{
  INT4 nspins = PULSAR_MAX_SPINS - 1;
  LIGOTimeGPS *epoch;
  UINT4 numSamples,numDetectors;
  REAL8 deltaref,deltaT;
  UINT4 i,j,k;

  /* sanity checks */
  if ( !(*Fa) || !(*Fa)->data || !(*Fb) || !(*Fb)->data )
    {
      XLALPrintError ("%s: empty input MultiCOMPLEX8TimeSeries!\n", __func__ );
      XLAL_ERROR (XLAL_EINVAL);
    }

  /* check if Fa and Fb have the same timeseries parameters */
  epoch = &(*Fa)->data[0]->epoch;
  numSamples = (*Fa)->data[0]->data->length;
  numDetectors = (*Fa)->length;
  deltaT = (*Fa)->data[0]->deltaT;
  if ( ( (*Fa)->length != numDetectors ) || ( (*Fb)->length != numDetectors ) )
    {
      XLALPrintError ("%s: Different numbers of detectors within the Fa and Fb MultiCOMPLEX8TimeSeries!\n", __func__ );
      XLAL_ERROR (XLAL_EINVAL);
    }
  for (i=0;i<(*Fa)->length;i++)
    {
      if ( ( XLALGPSCmp(epoch,&((*Fa)->data[i]->epoch)) != 0 ) || ( XLALGPSCmp(epoch,&((*Fb)->data[i]->epoch)) != 0 ) )
	{
	  XLALPrintError ("%s: Different start times for within Fa and Fb MultiCOMPLEX8TimeSeries!\n", __func__ );
	  XLAL_ERROR (XLAL_EINVAL);
	}
      if ( ( (*Fa)->data[i]->data->length != numSamples ) || ( (*Fb)->data[i]->data->length != numSamples ) )
	{
	  XLALPrintError ("%s: Different MultiCOMPLEX8TimeSeries lengths within Fa and Fb!\n", __func__ );
	  XLAL_ERROR (XLAL_EINVAL);
	}
      if ( ( (*Fa)->data[i]->deltaT != deltaT ) || ( (*Fb)->data[i]->deltaT != deltaT ) )
	{
	  XLALPrintError ("%s: Different MultiCOMPLEX8TimeSeries deltaT within Fa and Fb!\n", __func__ );
	  XLAL_ERROR (XLAL_EINVAL);
	}
    }

  /* determine number of spin down's and check if sensible */
  while (doppler->fkdot[nspins]==0.0) nspins--;
  if ( ( nspins < 0 ) || ( nspins > PULSAR_MAX_SPINS - 1 ) ) {
    XLALPrintError ("%s: Invalid number of spin derivitives, nspins = %d!\n", __func__,nspins );
    XLAL_ERROR (XLAL_EINVAL);
  }
 /*  printf("number of spin downs = %d\n",nspins); */

  /* compute the time difference between timeseries epoch and reference time */
  deltaref = XLALGPSGetREAL8(epoch) - XLALGPSGetREAL8(&(doppler->refTime));
 /*  printf("deltaref = %6.12f\n",deltaref); */
/*   printf("f0 = %6.12f\n",doppler->fkdot[0]); */
/*   printf("f1 = %6.12e\n",doppler->fkdot[1]); */

  /* apply spin derivitive correction to resampled timeseries */
  /* loop over spin derivitives (nspins = 1 means first derivitive, = 2 means second derivitive etc.. ) */
  for (j=1;j<=(UINT4)nspins;j++) {

    /* loop over time samples  and compute the spin down phase correction */
    for (k=0;k<numSamples;k++) {

      /* compute fractional number of cycles the spin-derivitive has added since the reftime */
      REAL8 cycles = fmod ( inv_fact[j+1]*doppler->fkdot[j]*pow(deltaref + k*deltaT,(REAL8)(j+1)), 1);
      REAL4 cosphase, sinphase;

      /* use look-up-table for speed to compute real and imaginary phase */
      XLAL_CHECK( XLALSinCos2PiLUT (&sinphase, &cosphase, -cycles ) == XLAL_SUCCESS, XLAL_EFUNC );

      /* loop over detectors */
      for (i=0;i<numDetectors;i++) {

	/* apply phase correction to Fa and Fb */
	REAL8 Fare = crealf((*Fa)->data[i]->data->data[k])*cosphase - cimagf((*Fa)->data[i]->data->data[k])*sinphase;
	REAL8 Faim = cimagf((*Fa)->data[i]->data->data[k])*cosphase + crealf((*Fa)->data[i]->data->data[k])*sinphase;
	REAL8 Fbre = crealf((*Fb)->data[i]->data->data[k])*cosphase - cimagf((*Fb)->data[i]->data->data[k])*sinphase;
	REAL8 Fbim = cimagf((*Fb)->data[i]->data->data[k])*cosphase + crealf((*Fb)->data[i]->data->data[k])*sinphase;

	(*Fa)->data[i]->data->data[k] = crectf( Fare, Faim );
	(*Fb)->data[i]->data->data[k] = crectf( Fbre, Fbim );

      } /* (i<numDetectors) */

    } /* (k<numSamples) */

  } /* (j<nspins) */

  return XLAL_SUCCESS;

 } /* XLALSpinDownCorrectionMultiFaFb */


/* ===== Object creation/destruction functions ===== */

/**
 * Destroy a MultiCOMPLEX8TimeSeries structure.
 * Note, this is "NULL-robust" in the sense that it will not crash
 * on NULL-entries anywhere in this struct, so it can be used
 * for failure-cleanup even on incomplete structs
 */
void
XLALDestroyMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries *multiTimes )
{
  UINT4 X;
  COMPLEX8TimeSeries *tmp;

  if ( ! multiTimes ) {
    return;
  }
  if ( multiTimes->data ) {

    for ( X=0; X < multiTimes->length; X ++ ) {

      if ( (tmp = multiTimes->data[X]) != NULL ) {
	XLALDestroyCOMPLEX8TimeSeries(tmp);
      }

    }
    LALFree ( multiTimes->data );
  }
  LALFree ( multiTimes );

  return;

} /* XLALDestroyMultiCOMPLEX8TimeSeries() */

/**
 * Duplicates a MultiCOMPLEX8TimeSeries structure.
 * Allocates memory and copies contents.
 */
MultiCOMPLEX8TimeSeries *
XLALDuplicateMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries *multiTimes )
{

  if ( ! multiTimes )
    XLAL_ERROR_NULL ( XLAL_EINVAL, "Invalid NULL input for multi-timeseries 'multiTimes'\n" );

  if ( multiTimes->length == 0 || multiTimes->data == NULL )
    XLAL_ERROR_NULL ( XLAL_EINVAL, "Invalid empty input timeseries, 0 detectors or data==NULL\n");


  UINT4 numDet = multiTimes->length;

  // ----- prepare memory for multicomplex8timeseries container
  MultiCOMPLEX8TimeSeries *out;
  if ( ( out = XLALCalloc ( 1, sizeof(*out) ) ) == NULL )
    XLAL_ERROR_NULL ( XLAL_ENOMEM, "Failed to calloc ( 1, %d )\n", sizeof(*out) );
  out->length = numDet;
  if ( ( out->data = XLALCalloc ( numDet, sizeof(*out->data) ) ) == NULL )
    XLAL_ERROR_NULL ( XLAL_ENOMEM, "Failed to calloc ( %d, %d )\n", numDet, sizeof(*out->data) );

  // ----- copy each of the numDet complex8timeseries contents
  for ( UINT4 X = 0; X < numDet; X ++ )
    {
      if ( (out->data[X] = XLALDuplicateCOMPLEX8TimeSeries ( multiTimes->data[X] )) == NULL )
        XLAL_ERROR_NULL ( XLAL_EFUNC );
    }

  return out;

} /* XLALDuplicateMultiCOMPLEX8TimeSeries() */

/**
 * Duplicates a COMPLEX8TimeSeries structure.
 * Allocates memory and copies contents.
 */
COMPLEX8TimeSeries *
XLALDuplicateCOMPLEX8TimeSeries ( COMPLEX8TimeSeries *times )
{
  if ( !times )
    XLAL_ERROR_NULL ( XLAL_EINVAL, "Invalid NULL input for timeseries 'times'\n" );

  if ( times->data == NULL || times->data->length == 0 || times->data->data == NULL )
    XLAL_ERROR_NULL ( XLAL_EINVAL, "Invalid empty input timeseries, 0 bins or data==NULL\n");

  COMPLEX8TimeSeries *out;
  if ( (out = XLALCalloc ( 1, sizeof(*out) ) ) == NULL )
    XLAL_ERROR_NULL ( XLAL_ENOMEM, "Failed to calloc ( 1, %d )\n", sizeof(*out) );

  // copy header info [including data-pointer, will be reset]
  memcpy ( out, times, sizeof(*times) );

  UINT4 numBins = times->data->length;
  if ( ( out->data = XLALCreateCOMPLEX8Vector ( numBins )) == NULL )
    XLAL_ERROR_NULL ( XLAL_EFUNC, "XLALCreateCOMPLEX8Vector( %d ) failed.\n", numBins );

  // copy contents of COMPLEX8 vector
  memcpy ( out->data->data, times->data->data, numBins * sizeof(*times->data->data) );

  return out;

} /* XLALDuplicateCOMPLEX8TimeSeries() */


/** Compare two COMPLEX8 vectors using various different comparison metrics
 */
int
XLALCompareCOMPLEX8Vectors ( VectorComparison *result,		///< [out] return comparison results
                             const COMPLEX8Vector *x,		///< [in] first input vector
                             const COMPLEX8Vector *y,		///< [in] second input vector
                             const VectorComparison *tol	///< [in] accepted tolerances on comparisons, or NULL for no check
                             )
{
  XLAL_CHECK ( result != NULL, XLAL_EINVAL );
  XLAL_CHECK ( x != NULL, XLAL_EINVAL );
  XLAL_CHECK ( y != NULL, XLAL_EINVAL );
  XLAL_CHECK ( x->data != NULL, XLAL_EINVAL );
  XLAL_CHECK ( y->data != NULL, XLAL_EINVAL );
  XLAL_CHECK ( x->length > 0, XLAL_EINVAL );
  XLAL_CHECK ( x->length == y->length, XLAL_EINVAL );

  REAL8 x_L1 = 0, x_L2 = 0;
  REAL8 y_L1 = 0, y_L2 = 0;
  REAL8 diff_L1 = 0, diff_L2 = 0;
  COMPLEX16 scalar = 0;

  REAL8 maxAbsx = 0, maxAbsy = 0;
  COMPLEX8 x_atMaxAbsx = 0, y_atMaxAbsx = 0;
  COMPLEX8 x_atMaxAbsy = 0, y_atMaxAbsy = 0;

  UINT4 numSamples = x->length;
  for ( UINT4 i = 0; i < numSamples; i ++ )
    {
      COMPLEX8 x_i = x->data[i];
      COMPLEX8 y_i = y->data[i];
      REAL8 xAbs_i = cabs ( x_i );
      REAL8 yAbs_i = cabs ( y_i );

      REAL8 absdiff = cabs ( x_i - y_i );
      diff_L1 += absdiff;
      diff_L2 += SQ(absdiff);

      x_L1 += xAbs_i;
      y_L1 += yAbs_i;
      x_L2 += SQ(xAbs_i);
      y_L2 += SQ(yAbs_i);

      scalar += x_i * conj(y_i);

      if ( xAbs_i > maxAbsx ) {
        maxAbsx = xAbs_i;
        x_atMaxAbsx = x_i;
        y_atMaxAbsx = y_i;
      }
      if ( yAbs_i > maxAbsy ) {
        maxAbsy = yAbs_i;
        x_atMaxAbsy = x_i;
        y_atMaxAbsy = y_i;
      }

    } // for i < numSamples

  // complete L2 norms by taking sqrt
  x_L2 = sqrt ( x_L2 );
  y_L2 = sqrt ( y_L2 );
  diff_L2 = sqrt ( diff_L2 );

  // compute and return comparison results
  result->relErr_L1 = diff_L1 / ( 0.5 * (x_L1 + y_L1 ) );
  result->relErr_L2 = diff_L2 / ( 0.5 * (x_L2 + y_L2 ) );
  REAL8 cosTheta = fmin ( 1, creal ( scalar ) / (x_L2 * y_L2) );
  result->angleV = acos ( cosTheta );
  result->relErr_atMaxAbsx = RELERR ( x_atMaxAbsx, y_atMaxAbsx );
  result->relErr_atMaxAbsy = RELERR ( x_atMaxAbsy, y_atMaxAbsy );;

  XLAL_CHECK ( XLALCheckVectorComparisonTolerances ( result, tol ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;

} // XLALCompareCOMPLEX8Vectors()

/** Compare two REAL4 vectors using various different comparison metrics
 */
int
XLALCompareREAL4Vectors ( VectorComparison *result,	///< [out] return comparison results
                          const REAL4Vector *x,		///< [in] first input vector
                          const REAL4Vector *y,		///< [in] second input vector
                          const VectorComparison *tol	///< [in] accepted tolerances on comparisons, or NULL for no check
                          )
{
  XLAL_CHECK ( result != NULL, XLAL_EINVAL );
  XLAL_CHECK ( x != NULL, XLAL_EINVAL );
  XLAL_CHECK ( y != NULL, XLAL_EINVAL );
  XLAL_CHECK ( x->data != NULL, XLAL_EINVAL );
  XLAL_CHECK ( y->data != NULL, XLAL_EINVAL );
  XLAL_CHECK ( x->length > 0, XLAL_EINVAL );
  XLAL_CHECK ( x->length == y->length, XLAL_EINVAL );

  REAL8 x_L1 = 0, x_L2 = 0;
  REAL8 y_L1 = 0, y_L2 = 0;
  REAL8 diff_L1 = 0, diff_L2 = 0;
  REAL8 scalar = 0;

  REAL4 maxAbsx = 0, maxAbsy = 0;
  REAL4 x_atMaxAbsx = 0, y_atMaxAbsx = 0;
  REAL4 x_atMaxAbsy = 0, y_atMaxAbsy = 0;

  UINT4 numSamples = x->length;
  for ( UINT4 i = 0; i < numSamples; i ++ )
    {
      REAL4 x_i = x->data[i];
      REAL4 y_i = y->data[i];
      REAL4 xAbs_i = fabs ( x_i );
      REAL4 yAbs_i = fabs ( y_i );

      REAL8 absdiff = fabs ( x_i - y_i );
      diff_L1 += absdiff;
      diff_L2 += SQ(absdiff);

      x_L1 += xAbs_i;
      y_L1 += yAbs_i;
      x_L2 += SQ(xAbs_i);
      y_L2 += SQ(yAbs_i);

      scalar += x_i * y_i;

      if ( xAbs_i > maxAbsx ) {
        maxAbsx = xAbs_i;
        x_atMaxAbsx = x_i;
        y_atMaxAbsx = y_i;
      }
      if ( yAbs_i > maxAbsy ) {
        maxAbsy = yAbs_i;
        x_atMaxAbsy = x_i;
        y_atMaxAbsy = y_i;
      }

    } // for i < numSamples

  // complete L2 norms by taking sqrt
  x_L2 = sqrt ( x_L2 );
  y_L2 = sqrt ( y_L2 );
  diff_L2 = sqrt ( diff_L2 );

  // compute and return comparison results
  result->relErr_L1 = diff_L1 / ( 0.5 * (x_L1 + y_L1 ) );
  result->relErr_L2 = diff_L2 / ( 0.5 * (x_L2 + y_L2 ) );
  REAL8 cosTheta = fmin ( 1, scalar / (x_L2 * y_L2) );
  result->angleV = acos ( cosTheta );
  result->relErr_atMaxAbsx = RELERR ( x_atMaxAbsx, y_atMaxAbsx );
  result->relErr_atMaxAbsy = RELERR ( x_atMaxAbsy, y_atMaxAbsy );;

  XLAL_CHECK ( XLALCheckVectorComparisonTolerances ( result, tol ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;

} // XLALCompareREAL4Vectors()

/** Check VectorComparison result against specified tolerances,
 * to allow for standardized comparison and reporting.
 */
int
XLALCheckVectorComparisonTolerances ( const VectorComparison *result,	///< [in] comparison result from XLALCompareREAL4Vectors() or XLALCompareCOMPLEX8Vectors()
                                      const VectorComparison *tol	///< [in] accepted tolerances on comparisons: if NULL=> always accept the result
                                      )
{
  XLAL_CHECK ( result != NULL, XLAL_EINVAL );

  VectorComparison XLAL_INIT_DECL(localTol);
  BOOLEAN failed = 0;

  if ( tol != NULL ) {
    localTol = (*tol);
    failed = ( (result->relErr_L1 > tol->relErr_L1) || (result->relErr_L2 > tol->relErr_L2) ||
               (result->angleV > tol->angleV) ||
               (result->relErr_atMaxAbsx > tol->relErr_atMaxAbsx) || (result->relErr_atMaxAbsy > tol->relErr_atMaxAbsy) );
  }

  if ( failed || (lalDebugLevel & LALINFO) )
    {
      XLALPrintInfo ( "relErr_L1        = %.1e (%.1e)\n"
                      "relErr_L2        = %.1e (%.1e)\n"
                      "angleV           = %.1e (%.1e)\n"
                      "relErr_atMaxAbsx = %.1e (%.1e)\n"
                      "relErr_atMaxAbsy = %.1e (%.1e)\n",
                      result->relErr_L1, localTol.relErr_L1, result->relErr_L2, localTol.relErr_L2, result->angleV,
                      localTol.angleV,
                      result->relErr_atMaxAbsx, localTol.relErr_atMaxAbsx,
                      result->relErr_atMaxAbsy, localTol.relErr_atMaxAbsy );
    }

  if ( failed ) {
    XLAL_ERROR ( XLAL_ETOL, "FAILED. Exceeded at least one tolerance level.\n");
  }

  XLALPrintInfo ("OK.\n");

  return XLAL_SUCCESS;

} // XLALCheckVectorComparisonTolerances()
