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

#include "config.h"

/* System includes */
#include <stdio.h>

/* LAL-includes */
#include <lal/LFTandTSutils.h>
#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/LogPrintf.h>
#include <lal/TimeSeries.h>
#include <lal/ComplexFFT.h>
#include <lal/ComputeFstat.h>

/*---------- DEFINES ----------*/
#define MYMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MYMIN(x,y) ( (x) < (y) ? (x) : (y) )
#define SQ(x) ((x)*(x))

#define NhalfPosDC(N) ((UINT4)(ceil ( ((N)/2.0 - 1e-6 ))))	/* round up */
#define NhalfNeg(N) ((UINT4)( (N) - NhalfPosDC(N) ))		/* round down (making sure N+ + N- = (N-1) */


/*---------- Global variables ----------*/

static LALUnit empty_LALUnit;

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
  numSFTsFit = (UINT4)floor(Tspan / Tsft + 0.5);	/* round */
  Tspan = numSFTsFit * Tsft;
  numTimeSamples = (UINT4)floor(Tspan / deltaT + 0.5);	/* round */

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
  if ( (lTS = XLALCreateCOMPLEX8TimeSeries ( firstSFT->name, &firstSFT->epoch, 0, deltaT, &empty_LALUnit, numTimeSamples )) == NULL )
    {
      XLALPrintError ("%s: XLALCreateCOMPLEX8TimeSeries() for %d timesteps failed! errno = %d!\n", __func__, numTimeSamples, xlalErrno );
      XLAL_ERROR_NULL ( XLAL_EFUNC );
    }
  memset ( lTS->data->data, 0, numTimeSamples * sizeof(*lTS->data->data)); /* set all time-samples to zero */

  /* ----- Prepare short time-series holding ONE invFFT of a single SFT */
  if ( (sTS = XLALCreateCOMPLEX8TimeSeries ( "short timeseries", &epoch, 0, deltaT, &empty_LALUnit, numBinsSFT )) == NULL )
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
      bin0_n = (UINT4) ( offset_n / deltaT + 0.5 );	/* round to closest bin */

      nudge_n = bin0_n * deltaT - offset_n;		/* rounding error */
      nudge_n = 1e-9 * (floor)(nudge_n * 1e9 + 0.5);	/* round to closest nanosecond */
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
      offsetEff = 1e-9 * (floor)( offsetEff * 1e9 + 0.5 );	/* round to closest integer multiple of nanoseconds */

      hetCycles = fmod ( fHet * offsetEff, 1);	/* required heterodyning phase-correction for this SFT */

      sin_cos_2PI_LUT (&hetCorr_im, &hetCorr_re, -hetCycles );

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
  
  numTimeSamples = (UINT4)floor(Tspan / deltaT + 0.5);	/* round */

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
  if ( (sTS = XLALCreateCOMPLEX8TimeSeries ( "short timeseries", &epoch, 0, deltaT, &empty_LALUnit, numBinsSFT )) == NULL )
    {
      XLALPrintError ( "%s: XLALCreateCOMPLEX8TimeSeries() for %d timesteps failed! errno = %d!\n", __func__, numBinsSFT, xlalErrno );
      goto failed;
    }

  /* ----- prepare long TimeSeries container ---------- */
  if ( (lTS = XLALCreateCOMPLEX8TimeSeries ( firstSFT->name, &start, fHet, deltaT, &empty_LALUnit, numTimeSamples )) == NULL )
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
      bin0_n = (UINT4) ( offset_n / deltaT + 0.5 );	/* round to closest bin */

      nudge_n = bin0_n * deltaT - offset_n;		/* rounding error */
      nudge_n = 1e-9 * (floor)(nudge_n * 1e9 + 0.5);	/* round to closest nanosecond */
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
      offsetEff = 1e-9 * (floor)( offsetEff * 1e9 + 0.5 );	/* round to closest integer multiple of nanoseconds */
      hetCycles = fmod ( fHet * offsetEff, 1);			/* required heterodyning phase-correction for this SFT */
    
      {
        REAL4 hetCorrection_re, hetCorrection_im;
        sin_cos_2PI_LUT( &hetCorrection_im, &hetCorrection_re, -hetCycles );
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

      sin_cos_2PI_LUT ( &fact_im, &fact_re, shiftCyles );

      yRe = fact_re * crealf(sft->data->data[k]) - fact_im * cimagf(sft->data->data[k]);
      yIm = fact_re * cimagf(sft->data->data[k]) + fact_im * crealf(sft->data->data[k]);

      sft->data->data[k] = crectf( yRe, yIm );

    } /* for k < numBins */

  /* adjust SFTs epoch to the shift */
  XLALGPSAdd( &sft->epoch, shift );

  return XLAL_SUCCESS;

} /* XLALTimeShiftSFT() */


