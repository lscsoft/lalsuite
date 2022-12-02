/*
 *  Copyright (C) 2014 Reinhard Prix
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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */

#include "config.h"

/* System includes */
#include <stdio.h>

/* LAL-includes */
#include <lal/LFTandTSutils.h>
#include <lal/SFTfileIO.h>
#include <lal/LogPrintf.h>
#include <lal/TimeSeries.h>
#include <lal/ComplexFFT.h>
#include <lal/ComputeFstat.h>
#include <lal/SinCosLUT.h>
#include <lal/Factorial.h>
#include <lal/Window.h>

/*---------- DEFINES ----------*/
#define MYMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MYMIN(x,y) ( (x) < (y) ? (x) : (y) )
#define SQ(x) ((x)*(x))

#define FINITE_OR_ONE(x)  (((x) != 0) ? (x) : 1.0)
#define cRELERR(x,y) ( cabsf( (x) - (y) ) / FINITE_OR_ONE( 0.5 * (cabsf(x) + cabsf(y)) ) )
#define fRELERR(x,y) ( fabsf( (x) - (y) ) / FINITE_OR_ONE( 0.5 * (fabsf(x) + fabsf(y)) ) )
#define OOTWOPI         (1.0 / LAL_TWOPI)      // 1/2pi
#define OOPI         (1.0 / LAL_PI)      // 1/pi
#define LD_SMALL4       (2.0e-4)                // "small" number for REAL4: taken from Demod()
/*---------- Global variables ----------*/
static LALUnit emptyLALUnit;

/* ---------- local prototypes ---------- */

/*---------- empty initializers ---------- */

/* ---------- function definitions ---------- */

/**
 * Turn the given multi-IFO SFTvectors into one long Fourier transform (LFT) over the total observation time
 */
SFTtype *
XLALSFTVectorToLFT ( SFTVector *sfts,		/**< input SFT vector (gets modified!) */
		     REAL8 upsampling		/**< upsampling factor >= 1 */
                     )
{
  XLAL_CHECK_NULL ( (sfts != NULL) && (sfts->length > 0), XLAL_EINVAL );
  XLAL_CHECK_NULL ( upsampling >= 1, XLAL_EDOM, "Upsampling factor (%f) must be >= 1 \n", upsampling );

  // ----- some useful SFT infos
  SFTtype *firstSFT = &(sfts->data[0]);
  UINT4 numBinsSFT = firstSFT->data->length;
  REAL8 dfSFT = firstSFT->deltaF;
  REAL8 Tsft = 1.0 / dfSFT;
  REAL8 f0SFT = firstSFT->f0;

  // ----- turn input SFTs into a complex (heterodyned) timeseries
  COMPLEX8TimeSeries *lTS;
  XLAL_CHECK_NULL ( (lTS = XLALSFTVectorToCOMPLEX8TimeSeries ( sfts )) != NULL, XLAL_EFUNC );
  REAL8 dt = lTS->deltaT;
  UINT4 numSamples0 = lTS->data->length;
  REAL8 Tspan0 = numSamples0 * dt;

  // ---------- determine time-span of upsampled time-series
  /* NOTE: Tspan MUST be an integer multiple of Tsft,
   * in order for the frequency bins of the final FFT
   * to be commensurate with the SFT bins.
   * This is required so that fHet is an exact
   * frequency-bin in both cases
   */
  UINT4 numSFTsFit = lround ( (Tspan0 * upsampling) / Tsft );
  REAL8 Tspan = numSFTsFit * Tsft;
  UINT4 numSamples = lround ( Tspan / dt );

  // ----- enlarge TimeSeries container for zero-padding if neccessary
  if ( numSamples > numSamples0 )
    {
      XLAL_CHECK_NULL ( (lTS->data->data = XLALRealloc ( lTS->data->data, numSamples * sizeof(lTS->data->data[0]) )) != NULL, XLAL_ENOMEM );
      lTS->data->length = numSamples;
      memset ( lTS->data->data + numSamples0, 0, (numSamples - numSamples0) * sizeof(lTS->data->data[0])); /* set all new time-samples to zero */
    }

  /* translate this back into fmin for the LFT (counting down from DC==fHet) */
  /* fHet = DC of our internal DFTs */
  UINT4 NnegSFT = NhalfNeg ( numBinsSFT );
  REAL8 fHet = f0SFT + 1.0 * NnegSFT * dfSFT;

  UINT4 NnegLFT = NhalfNeg ( numSamples );
  REAL8 f0LFT = fHet - NnegLFT / Tspan;

  // ----- prepare output LFT ----------
  SFTtype *outputLFT;
  XLAL_CHECK_NULL ( (outputLFT = XLALCreateSFT ( numSamples )) != NULL, XLAL_EFUNC );

  // prepare LFT header
  strcpy ( outputLFT->name, firstSFT->name );
  strncat ( outputLFT->name, ":long Fourier transform", sizeof(outputLFT->name) - 1 - strlen(outputLFT->name));
  outputLFT->epoch  = firstSFT->epoch;
  outputLFT->f0     = f0LFT;
  outputLFT->deltaF = 1.0 / Tspan;
  outputLFT->sampleUnits = firstSFT->sampleUnits;

  // ---------- FFT the long timeseries ----------
  COMPLEX8FFTPlan *LFTplan;
  XLAL_CHECK_NULL ( (LFTplan = XLALCreateForwardCOMPLEX8FFTPlan( numSamples, 0 )) != NULL, XLAL_EFUNC );
  XLAL_CHECK_NULL ( XLALCOMPLEX8VectorFFT( outputLFT->data, lTS->data, LFTplan ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_NULL ( XLALReorderFFTWtoSFT (outputLFT->data) == XLAL_SUCCESS, XLAL_EFUNC );
  // apply proper normalization 'dt'
  for ( UINT4 k = 0; k < outputLFT->data->length; k ++ ) {
    outputLFT->data->data[k] *= dt;
  }

  /* cleanup memory */
  XLALDestroyCOMPLEX8TimeSeries ( lTS );
  XLALDestroyCOMPLEX8FFTPlan ( LFTplan );

  return outputLFT;

} // XLALSFTVectorToLFT()


/**
 * Turn the given SFTvector into one long time-series, properly dealing with gaps.
 */
COMPLEX8TimeSeries *
XLALSFTVectorToCOMPLEX8TimeSeries ( const SFTVector *sftsIn         /**< [in] SFT vector */
                                    )
{
  // check input sanity
  XLAL_CHECK_NULL ( (sftsIn !=NULL) && (sftsIn->length > 0), XLAL_EINVAL );

  // create a local copy of the input SFTs, as they will be locally modified!
  SFTVector *sfts;
  XLAL_CHECK_NULL ( (sfts = XLALDuplicateSFTVector ( sftsIn )) != NULL, XLAL_EFUNC );

  /* define some useful shorthands */
  UINT4 numSFTs = sfts->length;
  SFTtype *firstSFT = &(sfts->data[0]);
  SFTtype *lastSFT = &(sfts->data[numSFTs-1]);
  UINT4 numFreqBinsSFT = firstSFT->data->length;
  REAL8 dfSFT = firstSFT->deltaF;
  REAL8 Tsft = 1.0 / dfSFT;
  REAL8 deltaT = Tsft / numFreqBinsSFT;	// complex FFT: numSamplesSFT = numFreqBinsSFT
  REAL8 f0SFT = firstSFT->f0;

  /* if the start and end input pointers are NOT NULL then determine start and time-span of the final long time-series */
  LIGOTimeGPS start = firstSFT->epoch;
  LIGOTimeGPS end = lastSFT->epoch;
  XLALGPSAdd ( &end, Tsft );

  /* determine output time span */
  REAL8 Tspan;
  XLAL_CHECK_NULL ( (Tspan = XLALGPSDiff ( &end, &start ) ) > 0, XLAL_EINVAL );

  UINT4 numSamples = lround ( Tspan / deltaT );

  /* determine the heterodyning frequency */
  /* fHet = DC of our internal DFTs */
  UINT4 NnegSFT = NhalfNeg ( numFreqBinsSFT );
  REAL8 fHet = f0SFT + 1.0 * NnegSFT * dfSFT;

  /* ----- Prepare invFFT of SFTs: compute plan for FFTW */
  COMPLEX8FFTPlan *SFTplan;
  XLAL_CHECK_NULL ( (SFTplan = XLALCreateReverseCOMPLEX8FFTPlan( numFreqBinsSFT, 0 )) != NULL, XLAL_EFUNC );

  /* ----- Prepare short time-series holding ONE invFFT of a single SFT */
  LIGOTimeGPS XLAL_INIT_DECL(epoch);
  COMPLEX8TimeSeries *sTS;
  XLAL_CHECK_NULL ( (sTS = XLALCreateCOMPLEX8TimeSeries ( "short timeseries", &epoch, 0, deltaT, &emptyLALUnit, numFreqBinsSFT )) != NULL, XLAL_EFUNC );

  /* ----- prepare long TimeSeries container ---------- */
  COMPLEX8TimeSeries *lTS;
  XLAL_CHECK_NULL ( (lTS = XLALCreateCOMPLEX8TimeSeries ( firstSFT->name, &start, fHet, deltaT, &emptyLALUnit, numSamples )) != NULL, XLAL_EFUNC );
  memset ( lTS->data->data, 0, numSamples * sizeof(*lTS->data->data)); 	/* set all time-samples to zero (in case there are gaps) */

  /* ---------- loop over all SFTs and inverse-FFT them ---------- */
  for ( UINT4 n = 0; n < numSFTs; n ++ )
    {
      SFTtype *thisSFT = &(sfts->data[n]);

      /* find bin in long timeseries corresponding to starttime of *this* SFT */
      REAL8 offset_n = XLALGPSDiff ( &(thisSFT->epoch), &start );
      UINT4 bin0_n = lround ( offset_n / deltaT );	/* round to closest bin */

      REAL8 nudge_n = bin0_n * deltaT - offset_n;		/* rounding error */
      nudge_n = 1e-9 * round ( nudge_n * 1e9 );	/* round to closest nanosecond */
      /* nudge SFT into integer timestep bin if necessary */
      XLAL_CHECK_NULL ( XLALTimeShiftSFT ( thisSFT, nudge_n ) == XLAL_SUCCESS, XLAL_EFUNC  );

      /* determine heterodyning phase-correction for this SFT */
      REAL8 offset = XLALGPSDiff ( &thisSFT->epoch, &start );	// updated value after time-shift
      // fHet * Tsft is an integer by construction, because fHet was chosen as a frequency-bin of the input SFTs
      // therefore we only need the remainder (offset % Tsft)
      REAL8 offsetEff = fmod ( offset, Tsft );
      REAL8 hetCycles = fmod ( fHet * offsetEff, 1); // heterodyning phase-correction for this SFT

      if ( nudge_n != 0 ){
        XLALPrintInfo("n = %d, offset_n = %g, nudge_n = %g, offset = %g, offsetEff = %g, hetCycles = %g\n",
                      n, offset_n, nudge_n, offset, offsetEff, hetCycles );
      }

      REAL4 hetCorrection_re, hetCorrection_im;
      XLAL_CHECK_NULL ( XLALSinCos2PiLUT ( &hetCorrection_im, &hetCorrection_re, -hetCycles ) == XLAL_SUCCESS, XLAL_EFUNC );
      COMPLEX8 hetCorrection = crectf( hetCorrection_re, hetCorrection_im );

      /* Note: we also bundle the overall normalization of 'df' into the het-correction.
       * This ensures that the resulting timeseries will have the correct normalization, according to
       * x_l = invFT[sft]_l = df * sum_{k=0}^{N-1} xt_k * e^(i 2pi k l / N )
       * where x_l is the l-th timestamp, and xt_k is the k-th frequency bin of the SFT.
       * See the LAL-conventions on FFTs:  http://www.ligo.caltech.edu/docs/T/T010095-00.pdf
       * (the FFTw convention does not contain the factor of 'df', which is why we need to
       * apply it ourselves)
       *
       */
      hetCorrection *= dfSFT;

      XLAL_CHECK_NULL ( XLALReorderSFTtoFFTW (thisSFT->data) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_NULL ( XLALCOMPLEX8VectorFFT( sTS->data, thisSFT->data, SFTplan ) == XLAL_SUCCESS, XLAL_EFUNC );

      for ( UINT4 j=0; j < sTS->data->length; j++) {
        sTS->data->data[j] *= hetCorrection;
      } // for j < numFreqBinsSFT

      // copy the short (shifted) heterodyned timeseries into correct location within long timeseries
      UINT4 binsLeft = numSamples - bin0_n;
      UINT4 copyLen = MYMIN ( numFreqBinsSFT, binsLeft );		/* make sure not to write past the end of the long TS */
      memcpy ( &lTS->data->data[bin0_n], sTS->data->data, copyLen * sizeof(lTS->data->data[0]) );

    } /* for n < numSFTs */

  // cleanup memory
  XLALDestroySFTVector ( sfts );
  XLALDestroyCOMPLEX8TimeSeries ( sTS );
  XLALDestroyCOMPLEX8FFTPlan ( SFTplan );

  return lTS;

} // XLALSFTVectorToCOMPLEX8TimeSeries()


/**
 * Turn the given multiSFTvector into multiple long COMPLEX8TimeSeries, properly dealing with gaps.
 * Memory allocation for the output MultiCOMPLEX8TimeSeries is done within this function.
 *
 */
MultiCOMPLEX8TimeSeries *
XLALMultiSFTVectorToCOMPLEX8TimeSeries ( const MultiSFTVector *multisfts  /**< [in] multi SFT vector */
                                         )
{
  // check input sanity
  XLAL_CHECK_NULL ( (multisfts != NULL) && (multisfts->length > 0), XLAL_EINVAL );
  UINT4 numDetectors = multisfts->length;

  /* allocate memory for the output structure */
  MultiCOMPLEX8TimeSeries *out;
  XLAL_CHECK_NULL ( (out = XLALMalloc ( sizeof(MultiCOMPLEX8TimeSeries) )) != NULL, XLAL_ENOMEM );
  XLAL_CHECK_NULL ( (out->data = XLALMalloc ( numDetectors * sizeof(COMPLEX8TimeSeries*) )) != NULL, XLAL_ENOMEM );
  out->length = numDetectors;

  /* loop over detectors */
  for ( UINT4 X=0; X < numDetectors; X++ ) {
    XLAL_CHECK_NULL ( (out->data[X] = XLALSFTVectorToCOMPLEX8TimeSeries ( multisfts->data[X] )) != NULL, XLAL_EFUNC );
  }

  return out;

} // XLALMultiSFTVectorToCOMPLEX8TimeSeries()



/**
 * Change frequency-bin order from fftw-convention to a 'SFT'
 * ie. from FFTW: f[0], f[1],...f[N/2], f[-(N-1)/2], ... f[-2], f[-1]
 * to: f[-(N-1)/2], ... f[-1], f[0], f[1], .... f[N/2]
 */
int
XLALReorderFFTWtoSFT ( COMPLEX8Vector *X )
{
  XLAL_CHECK ( (X != NULL) && (X->length > 0), XLAL_EINVAL );

  UINT4 N = X -> length;
  UINT4 Npos_and_DC = NhalfPosDC ( N );
  UINT4 Nneg = NhalfNeg ( N );

  /* allocate temporary storage for swap */
  COMPLEX8 *tmp;
  XLAL_CHECK ( (tmp = XLALMalloc ( N * sizeof(*tmp) )) != NULL, XLAL_ENOMEM );
  memcpy ( tmp, X->data, N * sizeof(*tmp) );

  /* Copy second half from FFTW: 'negative' frequencies */
  memcpy ( X->data, tmp + Npos_and_DC, Nneg * sizeof(*tmp) );

  /* Copy first half from FFTW: 'DC + positive' frequencies */
  memcpy ( X->data + Nneg, tmp, Npos_and_DC * sizeof(*tmp) );

  XLALFree(tmp);

  return XLAL_SUCCESS;

}  // XLALReorderFFTWtoSFT()

/**
 * Change frequency-bin order from 'SFT' to fftw-convention
 * ie. from f[-(N-1)/2], ... f[-1], f[0], f[1], .... f[N/2]
 * to FFTW: f[0], f[1],...f[N/2], f[-(N-1)/2], ... f[-2], f[-1]
 */
int
XLALReorderSFTtoFFTW (COMPLEX8Vector *X)
{
  XLAL_CHECK ( (X != NULL) && (X->length > 0), XLAL_EINVAL );

  UINT4 N = X->length;
  UINT4 Npos_and_DC = NhalfPosDC ( N );
  UINT4 Nneg = NhalfNeg ( N );

  /* allocate temporary storage for swap */
  COMPLEX8 *tmp;
  XLAL_CHECK ( (tmp = XLALMalloc ( N * sizeof(*tmp) )) != NULL, XLAL_ENOMEM );

  memcpy ( tmp, X->data, N * sizeof(*tmp) );

  /* Copy second half of FFTW: 'negative' frequencies */
  memcpy ( X->data + Npos_and_DC, tmp, Nneg * sizeof(*tmp) );

  /* Copy first half of FFTW: 'DC + positive' frequencies */
  memcpy ( X->data, tmp + Nneg, Npos_and_DC * sizeof(*tmp) );

  XLALFree(tmp);

  return XLAL_SUCCESS;

}  // XLALReorderSFTtoFFTW()

/**
 * Time-shift the given SFT by an amount of 'shift' seconds,
 * using the frequency-domain expression
 * \f$\widetilde{y}(f) = \widetilde{x}(f) \, e^{i 2\pi\,f\,\tau}\f$,
 * which shifts \f$x(t)\f$ into \f$y(t) = x(t + \tau)\f$
 *
 * NOTE: this <b>modifies</b> the SFT in place
 */
int
XLALTimeShiftSFT ( SFTtype *sft,	/**< [in/out] SFT to time-shift */
		   REAL8 shift		/**< time-shift \f$\tau\f$ in seconds */
                   )
{
  XLAL_CHECK ( (sft != NULL) && (sft->data != NULL), XLAL_EINVAL );

  if ( shift == 0 ) {
    return XLAL_SUCCESS;
  }

  for ( UINT4 k=0; k < sft->data->length; k++ )
    {
      REAL8 fk = sft->f0 + k * sft->deltaF;	/* frequency of k-th bin */
      REAL8 shiftCyles = shift * fk;
      REAL4 fact_re, fact_im;			/* complex phase-shift factor e^(-2pi f tau) */
      XLAL_CHECK ( XLALSinCos2PiLUT ( &fact_im, &fact_re, shiftCyles ) == XLAL_SUCCESS, XLAL_EFUNC );

      COMPLEX8 fact = crectf(fact_re, fact_im);
      sft->data->data[k] *= fact;

    } /* for k < numBins */

  /* adjust SFTs epoch to the shift */
  XLALGPSAdd ( &sft->epoch, shift );

  return XLAL_SUCCESS;

} // XLALTimeShiftSFT()

/**
 * Multi-detector wrapper for XLALFrequencyShiftCOMPLEX8TimeSeries
 * NOTE: this <b>modifies</b> the MultiCOMPLEX8Timeseries in place
 */
int
XLALFrequencyShiftMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries *x,		/**< [in/out] timeseries to time-shift */
					    const REAL8 shift	                /**< [in] freq-shift in Hz */
                                            )
{
  XLAL_CHECK ( (x != NULL) && ( x->data != NULL ) && (x->length >  0), XLAL_EINVAL );

  if ( shift == 0 ) {
    return XLAL_SUCCESS;
  }

  /* loop over detectors */
  UINT4 numDetectors = x->length;
  for ( UINT4 X=0; X < numDetectors; X++ )
    {
      /* shift the frequency of each detector's data */
      XLAL_CHECK ( XLALFrequencyShiftCOMPLEX8TimeSeries ( x->data[X], shift ) == XLAL_SUCCESS, XLAL_EFUNC );
    }

  return XLAL_SUCCESS;

} /* XLALFrequencyShiftMultiCOMPLEX8TimeSeries() */


/**
 * Freq-shift the given COMPLEX8Timeseries by an amount of 'shift' Hz,
 * using the time-domain expression y(t) = x(t) * e^(-i 2pi df t),
 * which shifts x(f) into y(f) = x(f + df)
 *
 * NOTE: this <b>modifies</b> the COMPLEX8TimeSeries in place
 */
int
XLALFrequencyShiftCOMPLEX8TimeSeries ( COMPLEX8TimeSeries *x,	        /**< [in/out] timeseries to time-shift */
				       const REAL8 shift	        /**< [in] freq-shift in Hz */
                                       )
{
  XLAL_CHECK ( (x != NULL) && ( x->data != NULL ), XLAL_EINVAL );

  if ( shift == 0 ) {
    return XLAL_SUCCESS;
  }

  /* get timeseries time-step */
  REAL8 deltat = x->deltaT;

  /* loop over COMPLEX8TimeSeries elements */
  for ( UINT4 k=0; k < x->data->length; k++)
    {
      REAL8 tk = k * deltat;	/* time of k-th bin */
      REAL8 shiftCycles = shift * tk;
      REAL4 fact_re, fact_im;			/* complex phase-shift factor e^(-2pi f tau) */

      /* use a sin/cos look-up-table for speed */
      XLAL_CHECK( XLALSinCos2PiLUT ( &fact_im, &fact_re, shiftCycles ) == XLAL_SUCCESS, XLAL_EFUNC );
      COMPLEX8 fact = crectf(fact_re, fact_im);

      x->data->data[k] *= fact;

    } /* for k < numBins */

  /* adjust timeseries heterodyne frequency to the shift */
  x->f0 -= shift;

  return XLAL_SUCCESS;

} /* XLALFrequencyShiftCOMPLEX8TimeSeries() */


/**
 * Apply a spin-down correction to the complex8 timeseries
 * using the time-domain expression y(t) = x(t) * e^(-i 2pi sum f_k * (t-tref)^(k+1)),
 *
 * NOTE: this <b>modifies</b> the input COMPLEX8TimeSeries in place
 */
int
XLALSpinDownCorrectionMultiTS ( MultiCOMPLEX8TimeSeries *multiTimeSeries,       /**< [in/out] timeseries to time-shift */
                                const PulsarDopplerParams *doppler		/**< parameter-space point to correct for */
                                )
{
  // check input sanity
  XLAL_CHECK ( (multiTimeSeries != NULL) &&  (multiTimeSeries->data != NULL) && (multiTimeSeries->length > 0), XLAL_EINVAL );
  XLAL_CHECK ( doppler != NULL, XLAL_EINVAL );

  UINT4 numDetectors = multiTimeSeries->length;

  LIGOTimeGPS *epoch = &(multiTimeSeries->data[0]->epoch);
  UINT4 numSamples = multiTimeSeries->data[0]->data->length;

  REAL8 dt = multiTimeSeries->data[0]->deltaT;

  /* determine number of spin down's and check if sensible */
  UINT4 nspins = PULSAR_MAX_SPINS - 1;
  while ( (nspins > 0) && (doppler->fkdot[nspins] == 0) ) {
    nspins--;
  }

  /* apply spin derivitive correction to resampled timeseries */
  REAL8 tk = XLALGPSDiff ( epoch, &(doppler->refTime) );
  for ( UINT4 k=0; k < numSamples; k++ )
    {
      REAL8 tk_pow_jp1 = tk;
      REAL8 cycles_k = 0;
      for ( UINT4 j=1; j <= nspins; j++ )
        {
          tk_pow_jp1 *= tk;
          /* compute fractional number of cycles the spin-derivitive has added since the reftime */
          cycles_k += LAL_FACT_INV[j+1] * doppler->fkdot[j] * tk_pow_jp1;
        } // for j < nspins

      REAL4 cosphase, sinphase;
      XLAL_CHECK( XLALSinCos2PiLUT ( &sinphase, &cosphase, -cycles_k ) == XLAL_SUCCESS, XLAL_EFUNC );
      COMPLEX8 em2piphase = crectf ( cosphase, sinphase );

      /* loop over detectors */
      for ( UINT4 X=0; X < numDetectors; X++ )
        {
          multiTimeSeries->data[X]->data->data[k] *= em2piphase;
        } // for X < numDetectors

      tk += dt;

    } // for k < numSamples

  return XLAL_SUCCESS;

} // XLALSpinDownCorrectionMultiTS()


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
  if ( ! multiTimes ) {
    return;
  }
  if ( multiTimes->data != NULL )
    {
      UINT4 numDetectors = multiTimes->length;
      for ( UINT4 X=0; X < numDetectors; X ++ ) {
        XLALDestroyCOMPLEX8TimeSeries ( multiTimes->data[X] );
      } // for X < numDetectors
      XLALFree ( multiTimes->data );
    } // if multiTimes->data

  XLALFree ( multiTimes );

  return;

} // XLALDestroyMultiCOMPLEX8TimeSeries()


/**
 * Duplicates a MultiCOMPLEX8TimeSeries structure.
 * Allocates memory and copies contents.
 */
MultiCOMPLEX8TimeSeries *
XLALDuplicateMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries *multiTimes )
{
  XLAL_CHECK_NULL ( multiTimes != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( multiTimes->length > 0, XLAL_EINVAL );

  UINT4 numDetectors = multiTimes->length;

  // ----- prepare memory for multicomplex8timeseries container
  MultiCOMPLEX8TimeSeries *out;
  XLAL_CHECK_NULL ( (out = XLALCalloc ( 1, sizeof(*out) )) != NULL, XLAL_ENOMEM );
  out->length = numDetectors;
  XLAL_CHECK_NULL ( (out->data = XLALCalloc ( numDetectors, sizeof(*out->data) )) != NULL, XLAL_ENOMEM );

  // ----- copy each of the numDet complex8timeseries contents
  for ( UINT4 X = 0; X < numDetectors; X ++ ) {
    XLAL_CHECK_NULL ( (out->data[X] = XLALDuplicateCOMPLEX8TimeSeries ( multiTimes->data[X] )) != NULL, XLAL_EFUNC );
  }

  return out;

} // XLALDuplicateMultiCOMPLEX8TimeSeries()

/**
 * Duplicates a COMPLEX8TimeSeries structure.
 * Allocates memory and copies contents.
 */
COMPLEX8TimeSeries *
XLALDuplicateCOMPLEX8TimeSeries ( COMPLEX8TimeSeries *ttimes )
{
  XLAL_CHECK_NULL ( ttimes != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( (ttimes->data != NULL) && (ttimes->data->length > 0) && ( ttimes->data->data != NULL ), XLAL_EINVAL );

  COMPLEX8TimeSeries *out;
  XLAL_CHECK_NULL ( (out = XLALCalloc ( 1, sizeof(*out) )) != NULL, XLAL_ENOMEM );

  // copy header info [including data-pointer, will be reset]
  memcpy ( out, ttimes, sizeof(*ttimes) );

  UINT4 numBins = ttimes->data->length;
  XLAL_CHECK_NULL ( (out->data = XLALCreateCOMPLEX8Vector ( numBins )) != NULL, XLAL_EFUNC );

  // copy contents of COMPLEX8 vector
  memcpy ( out->data->data, ttimes->data->data, numBins * sizeof(ttimes->data->data[0]) );

  return out;

} // XLALDuplicateCOMPLEX8TimeSeries()

/**
 * Copies a MultiCOMPLEX8TimeSeries structure, output must be allocated of identical size as input!
 */
int
XLALCopyMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries *multiTimesOut,
                                  MultiCOMPLEX8TimeSeries *multiTimesIn
                                  )
{
  XLAL_CHECK ( multiTimesIn != NULL, XLAL_EINVAL );
  XLAL_CHECK ( multiTimesOut != NULL, XLAL_EINVAL );

  UINT4 numDetectors = multiTimesIn->length;
  XLAL_CHECK ( multiTimesOut->length == numDetectors, XLAL_EINVAL );

  // ----- copy each of the numDet complex8timeseries contents
  for ( UINT4 X = 0; X < numDetectors; X ++ ) {
    XLAL_CHECK ( XLALCopyCOMPLEX8TimeSeries ( multiTimesOut->data[X], multiTimesIn->data[X] ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  return XLAL_SUCCESS;

} // XLALCopyMultiCOMPLEX8TimeSeries()

/**
 * Copies a COMPLEX8TimeSeries structure, output must be allocated of identical size as input!
 */
int
XLALCopyCOMPLEX8TimeSeries ( COMPLEX8TimeSeries *ts_out,
                             COMPLEX8TimeSeries *ts_in
                             )
{
  XLAL_CHECK ( ts_out != NULL, XLAL_EINVAL );
  XLAL_CHECK ( ts_in != NULL, XLAL_EINVAL );

  UINT4 numSamples = ts_in->data->length;
  XLAL_CHECK ( ts_out->data->length == numSamples, XLAL_EINVAL );

  COMPLEX8Vector *tmp = ts_out->data;

  // copy header info [including data-pointer, will be restored]
  (*ts_out) = (*ts_in);
  ts_out->data = tmp;

  // copy contents of COMPLEX8 vector
  memcpy ( ts_out->data->data, ts_in->data->data, numSamples * sizeof(ts_out->data->data[0]) );

  return XLAL_SUCCESS;

} // XLALCopyCOMPLEX8TimeSeries()


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
  REAL8 x_L2sq = 0;
  REAL8 y_L1 = 0, y_L2 = 0;
  REAL8 y_L2sq = 0;
  REAL8 diff_L1 = 0, diff_L2 = 0;
  COMPLEX16 scalar = 0;

  REAL8 maxAbsx = 0, maxAbsy = 0;
  COMPLEX8 x_atMaxAbsx = 0, y_atMaxAbsx = 0;
  COMPLEX8 x_atMaxAbsy = 0, y_atMaxAbsy = 0;


  UINT4 numSamples = x->length;
  for ( UINT4 i = 0; i < numSamples; i ++ )
    {
      COMPLEX16 x_i = x->data[i];
      COMPLEX16 y_i = y->data[i];
      REAL8 xAbs2_i = x_i * conj( x_i );
      REAL8 yAbs2_i = y_i * conj( y_i );
      REAL8 xAbs_i  = sqrt ( xAbs2_i );
      REAL8 yAbs_i  = sqrt ( yAbs2_i );
      XLAL_CHECK ( isfinite ( xAbs_i ), XLAL_EFPINVAL, "non-finite element: x(%d) = %g + I %g\n", i, creal(x_i), cimag(x_i) );
      XLAL_CHECK ( isfinite ( yAbs_i ), XLAL_EFPINVAL, "non-finite element: y(%d) = %g + I %g\n", i, creal(y_i), cimag(y_i) );

      REAL8 absdiff = cabs ( x_i - y_i );
      diff_L1 += absdiff;
      diff_L2 += SQ(absdiff);

      x_L1 += xAbs_i;
      y_L1 += yAbs_i;
      x_L2sq += xAbs2_i;
      y_L2sq += yAbs2_i;

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
  x_L2 = sqrt ( x_L2sq );
  y_L2 = sqrt ( y_L2sq );
  diff_L2 = sqrt ( diff_L2 );
  REAL8 sinAngle2 = (x_L2sq * y_L2sq - scalar * conj(scalar)) / FINITE_OR_ONE(x_L2sq * y_L2sq);
  REAL8 angle = asin ( sqrt( fabs(sinAngle2) ) );

  // compute and return comparison results
  result->relErr_L1 = diff_L1 / FINITE_OR_ONE( 0.5 * (x_L1 + y_L1 ) );
  result->relErr_L2 = diff_L2 / FINITE_OR_ONE( 0.5 * (x_L2 + y_L2 ) );
  result->angleV = angle;

  result->relErr_atMaxAbsx = cRELERR ( x_atMaxAbsx, y_atMaxAbsx );
  result->relErr_atMaxAbsy = cRELERR ( x_atMaxAbsy, y_atMaxAbsy );;

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
  REAL8 x_L2sq = 0;
  REAL8 y_L1 = 0, y_L2 = 0;
  REAL8 y_L2sq = 0;
  REAL8 diff_L1 = 0, diff_L2 = 0;
  REAL8 scalar = 0;

  REAL4 maxAbsx = 0, maxAbsy = 0;
  REAL4 x_atMaxAbsx = 0, y_atMaxAbsx = 0;
  REAL4 x_atMaxAbsy = 0, y_atMaxAbsy = 0;

  UINT4 numSamples = x->length;
  for ( UINT4 i = 0; i < numSamples; i ++ )
    {
      REAL8 x_i = x->data[i];
      REAL8 y_i = y->data[i];
      XLAL_CHECK ( isfinite ( x_i ), XLAL_EFPINVAL, "non-finite element: x(%d) = %g\n", i, x_i );
      XLAL_CHECK ( isfinite ( y_i ), XLAL_EFPINVAL, "non-finite element: y(%d) = %g\n", i, y_i );
      REAL8 xAbs_i = fabs ( x_i );
      REAL8 yAbs_i = fabs ( y_i );
      REAL8 xAbs2_i = SQ( x_i );
      REAL8 yAbs2_i = SQ( y_i );

      REAL8 absdiff = fabs ( x_i - y_i );
      diff_L1 += absdiff;
      diff_L2 += SQ(absdiff);

      x_L1 += xAbs_i;
      y_L1 += yAbs_i;
      x_L2sq += xAbs2_i;
      y_L2sq += yAbs2_i;

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
  x_L2 = sqrt ( x_L2sq );
  y_L2 = sqrt ( y_L2sq );
  diff_L2 = sqrt ( diff_L2 );

  // compute and return comparison results
  result->relErr_L1 = diff_L1 / FINITE_OR_ONE( 0.5 * (x_L1 + y_L1 ) );
  result->relErr_L2 = diff_L2 / FINITE_OR_ONE( 0.5 * (x_L2 + y_L2 ) );
  REAL8 sinAngle2 = (x_L2sq * y_L2sq - SQ(scalar) ) / FINITE_OR_ONE( x_L2sq * y_L2sq );
  result->angleV = asin ( sqrt ( fabs ( sinAngle2 ) ) );
  result->relErr_atMaxAbsx = fRELERR ( x_atMaxAbsx, y_atMaxAbsx );
  result->relErr_atMaxAbsy = fRELERR ( x_atMaxAbsy, y_atMaxAbsy );;

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
  if ( tol != NULL ) {
    localTol = (*tol);
  }    

  const char* names[5]   = { "relErr_L1",        "relErr_L2",        "angleV",        "relErr_atMaxAbsx",        "relErr_atMaxAbsy"        };
  const REAL4 results[5] = { result->relErr_L1,  result->relErr_L2,  result->angleV,  result->relErr_atMaxAbsx,  result->relErr_atMaxAbsy  };
  const REAL4 tols[5]    = { localTol.relErr_L1, localTol.relErr_L2, localTol.angleV, localTol.relErr_atMaxAbsx, localTol.relErr_atMaxAbsy };

  BOOLEAN failed = 0;

  for (size_t i = 0; i < XLAL_NUM_ELEM(names); ++i) {

    if ( results[i] > tols[i] ) {
      failed = 1;
    }

    if (lalDebugLevel & LALINFO ) {
      XLALPrintError( "%-16s = %.1e (%.1e): %s\n", names[i], results[i], tols[i],
                      ( results[i] > tols[i] ) ? "FAILED. Exceeded tolerance." : "OK." );
    }

  }

  if ( failed ) {
    XLAL_ERROR ( XLAL_ETOL, "FAILED. Exceeded at least one tolerance level.\n");
  }

  XLALPrintInfo ("OK.\n");

  return XLAL_SUCCESS;

} // XLALCheckVectorComparisonTolerances()

/** Interpolate a given regularly-spaced COMPLEX8 timeseries 'ts_in = x_in(j * dt)' onto new samples
 *  'y_out(t_out)' using *windowed* Shannon sinc interpolation, windowed to (2*Dterms+1) terms, namely
 * \f[
 * \newcommand{\Dterms}{\mathrm{Dterms}}
 * \f]
 *
 * \f{equation}{
 * x(t) = \sum_{j = j^* - \Dterms}^{j^* + \Dterms} x_j \, w_j \, \frac{\sin(\pi\delta_j)}{\pi\delta_j}\,,\quad\text{with}\quad
 * \delta_j \equiv \frac{t - t_j}{\Delta t}\,,
 * \f}
 * where \f$j^* \equiv \mathrm{round}(t / \Delta t)\f$, and
 * where \f$w_j\f$ is the window used (here: Hamming)
 *
 * In order to implement this more efficiently, we observe that \f$\sin(\pi\delta_j) = (-1)^{(j-j0)}\sin(\pi\delta_{j0})\f$ for integer \f$j\f$,
 * and therefore
 *
 * \f{equation}{
 * x(t) = \frac{\sin(\pi\,\delta_{j0})}{\pi} \, \sum_{j = j^* - \Dterms}^{j^* + \Dterms} (-1)^{(j-j0)}\frac{x_j \, w_j}{\delta_j}\,,
 * \f}
 *
 * NOTE: Using Dterms=0 corresponds to closest-bin interpolation
 *
 * NOTE2: samples *outside* the original timespan are returned as 0
 *
 * NOTE3: we're using a Hamming window of length L = 2*Dterms + 1, which has a pass-band wiggle of delta_p ~ 0.0022,
 * and a transition bandwidth of (4/L) * Bandwidth.
 * You need to make sure to include sufficient effective sidebands to the input timeseries, so that the transition band can
 * be safely ignored or 'cut out' at the end
 */
int
XLALSincInterpolateCOMPLEX8TimeSeries ( COMPLEX8Vector *y_out,		///< [out] output series of interpolated y-values [must be same size as t_out]
                                        const REAL8Vector *t_out,	///< [in] output time-steps to interpolate input to
                                        const COMPLEX8TimeSeries *ts_in,///< [in] regularly-spaced input timeseries
                                        UINT4 Dterms			///< [in] window sinc kernel sum to +-Dterms around max
                                        )
{
  XLAL_CHECK ( y_out != NULL, XLAL_EINVAL );
  XLAL_CHECK ( t_out != NULL, XLAL_EINVAL );
  XLAL_CHECK ( ts_in != NULL, XLAL_EINVAL );
  XLAL_CHECK ( y_out->length == t_out->length, XLAL_EINVAL );

  UINT4 numSamplesOut = t_out->length;
  UINT4 numSamplesIn = ts_in->data->length;
  REAL8 dt = ts_in->deltaT;
  REAL8 tmin = XLALGPSGetREAL8 ( &(ts_in->epoch) );	// time of first bin in input timeseries

  REAL8Window *win;
  UINT4 winLen = 2 * Dterms + 1;
  XLAL_CHECK ( (win = XLALCreateHammingREAL8Window ( winLen )) != NULL, XLAL_EFUNC );

  const REAL8 oodt = 1.0 / dt;

  for ( UINT4 l = 0; l < numSamplesOut; l ++ )
    {
      REAL8 t = t_out->data[l] - tmin;		// measure time since start of input timeseries

      // samples outside of input timeseries are returned as 0
      if ( (t < 0) || (t > (numSamplesIn-1)*dt) )	// avoid any extrapolations!
        {
          y_out->data[l] = 0;
          continue;
        }

      REAL8 t_by_dt = t  * oodt;
      INT8 jstar = lround ( t_by_dt );		// bin closest to 't', guaranteed to be in [0, numSamples-1]

      if ( fabs ( t_by_dt - jstar ) < LD_SMALL4 )	// avoid numerical problems near peak
        {
          y_out->data[l] = ts_in->data->data[jstar];	// known analytic solution for exact bin
          continue;
        }

      INT4 jStart0 = jstar - Dterms;
      UINT4 jEnd0 = jstar + Dterms;
      UINT4 jStart = MYMAX ( jStart0, 0 );
      UINT4 jEnd   = MYMIN ( jEnd0, numSamplesIn - 1 );

      REAL4 delta_jStart = (t_by_dt - jStart);
      REAL4 sin0, cos0;
      XLALSinCosLUT ( &sin0, &cos0, LAL_PI * delta_jStart );
      REAL4 sin0oopi = sin0 * OOPI;

      COMPLEX8 y_l = 0;
      REAL8 delta_j = delta_jStart;
      for ( UINT8 j = jStart; j <= jEnd; j ++ )
        {
          COMPLEX8 Cj = win->data->data[j - jStart0] * sin0oopi / delta_j;

          y_l += Cj * ts_in->data->data[j];

          sin0oopi = -sin0oopi;		// sin-term flips sign every step
          delta_j --;
        } // for j in [j* - Dterms, ... ,j* + Dterms]

      y_out->data[l] = y_l;

    } // for l < numSamplesOut

  XLALDestroyREAL8Window ( win );

  return XLAL_SUCCESS;

} // XLALSincInterpolateCOMPLEX8TimeSeries()

/** Interpolate a given regularly-spaced COMPLEX8 frequency-series 'fs_in = x_in( k * df)' onto new samples
 *  'y_out(f_out)' using (complex) Sinc interpolation (obtained from Dirichlet kernel in large-N limit), truncated to (2*Dterms+1) terms, namely
 *
 * \f{equation}{
 * \widetilde{x}(f) = \sum_{k = k^* - \Delta k}^{k^* + \Delta k} \widetilde{x}_k \, \frac{ 1 - e^{-i 2\pi\,\delta_k} }{ i 2\pi\,\delta_k}\,,
 * \f}
 * where \f$k^* \equiv \mathrm{round}(f / \Delta f)\f$, and \f$\delta_k \equiv \frac{f - f_k}{\Delta f}\f$.
 *
 * In order to implement this more efficiently, we observe that \f$e^{-i 2\pi\,\delta_k} = e^{-i 2\pi\,\delta_{k*}}\f$ for integer \f$k\f$, and therefore
 * \f{equation}{
 * \widetilde{x}(f) = \frac{\sin(2\pi\delta_{k*}) + i(\cos(2\pi\delta_{k*}) - 1)}{2\pi} \sum_{k = k^* - \Delta k}^{k^* + \Delta k} \frac{\widetilde{x}_k}{\delta_k}\,,
 * \f}
 *
 * NOTE: Using Dterms=0 corresponds to closest-bin interpolation
 *
 * NOTE2: frequencies *outside* the original frequency band are returned as 0
 */
int
XLALSincInterpolateCOMPLEX8FrequencySeries ( COMPLEX8Vector *y_out,			///< [out] output series of interpolated y-values [must be alloc'ed already]
                                                  const REAL8Vector *f_out,		///< [in] output frequency-values to interpolate input to
                                                  const COMPLEX8FrequencySeries *fs_in,	///< [in] regularly-spaced input frequency-series
                                                  UINT4 Dterms				///< [in] truncate interpolation kernel sum to +-Dterms around max
                                                  )
{
  XLAL_CHECK ( y_out != NULL, XLAL_EINVAL );
  XLAL_CHECK ( f_out != NULL, XLAL_EINVAL );
  XLAL_CHECK ( fs_in != NULL, XLAL_EINVAL );

  UINT4 numBinsOut = f_out->length;
  XLAL_CHECK ( y_out->length == numBinsOut, XLAL_EINVAL );

  UINT4 numBinsIn = fs_in->data->length;
  REAL8 df = fs_in->deltaF;
  REAL8 fMin = fs_in->f0;	// frequency of first input bin

  const REAL8 oodf = 1.0 / df;

  for ( UINT4 l = 0; l < numBinsOut; l ++ )
    {
      REAL8 f = f_out->data[l] - fMin;		// measure frequency from lowest input bin

      // bins outside of input frequency-series are returned as 0
      if ( (f < 0) || (f > (numBinsIn-1)*df) )	// avoid any extrapolations!
        {
          y_out->data[l] = 0;
          continue;
        }

      INT8 kstar = lround ( f * oodf );		// bin closest to 't', guaranteed to be in [0, numBins-1]

      if ( fabs ( f - kstar * df ) < 1e-6 )	// avoid numerical problems near peak
        {
          y_out->data[l] = fs_in->data->data[kstar];	// known analytic solution for exact bin
          continue;
        }

      UINT8 k0 = MYMAX ( kstar - Dterms, 0 );
      UINT8 k1 = MYMIN ( kstar + Dterms, numBinsIn - 1 );

      REAL4 delta_kstar = (f - kstar * df) * oodf;	// in [-0.5, 0.5]
      REAL4 sin2pidelta, cos2pidelta;
      XLALSinCos2PiLUTtrimmed ( &sin2pidelta, &cos2pidelta, delta_kstar + 1.0f );

      COMPLEX8 prefact = OOTWOPI * ( sin2pidelta + I * (cos2pidelta - 1.0f) );

      COMPLEX8 y_l = 0;
      REAL4 delta_k = delta_kstar + Dterms;
      for ( UINT8 k = k0; k <= k1; k ++ )
        {
          y_l += fs_in->data->data[k] / delta_k;
          delta_k --;
        } // for k in [k* - Dterms, ... ,k* + Dterms]

      y_out->data[l] = prefact * y_l;

    } // for l < numBinsOut

  return XLAL_SUCCESS;

} // XLALSincInterpolateCOMPLEX8FrequencySeries()


/** (Complex)Sinc-interpolate an input SFT to an output SFT.
 * This is a simple convenience wrapper to XLALSincInterpolateCOMPLEX8FrequencySeries()
 * for the special case of interpolating onto new SFT frequency bins
 */
SFTtype *
XLALSincInterpolateSFT ( const SFTtype *sft_in,		///< [in] input SFT
                              REAL8 f0Out,		///< [in] new start frequency
                              REAL8 dfOut,		///< [in] new frequency step-size
                              UINT4 numBinsOut,		///< [in] new number of bins
                              UINT4 Dterms		///< [in] truncate interpolation kernel sum to +-Dterms around max
                              )
{
  XLAL_CHECK_NULL ( sft_in != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( dfOut > 0, XLAL_EINVAL );
  XLAL_CHECK_NULL ( numBinsOut > 0, XLAL_EINVAL );

  // setup frequency vector
  REAL8Vector *f_out;
  XLAL_CHECK_NULL ( (f_out = XLALCreateREAL8Vector ( numBinsOut )) != NULL, XLAL_EFUNC );
  for ( UINT4 k = 0; k < numBinsOut; k ++ ) {
    f_out->data[k] = f0Out + k * dfOut;
  } // for k < numBinsOut

  SFTtype *out;
  XLAL_CHECK_NULL ( (out = XLALCalloc ( 1, sizeof(*out))) != NULL, XLAL_EFUNC );
  (*out) = (*sft_in);	// copy header
  out->f0 = f0Out;
  out->deltaF = dfOut;
  XLAL_CHECK_NULL ( (out->data = XLALCreateCOMPLEX8Vector ( numBinsOut )) != NULL, XLAL_EFUNC );

  XLAL_CHECK_NULL ( XLALSincInterpolateCOMPLEX8FrequencySeries ( out->data, f_out, sft_in, Dterms ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLALDestroyREAL8Vector ( f_out );

  return out;

} // XLALSincInterpolateSFT()


/**
 * Interpolate frequency-series to newLen frequency-bins.
 * This is using DFT-interpolation (derived from zero-padding).
 */
COMPLEX8Vector *
XLALrefineCOMPLEX8Vector ( const COMPLEX8Vector *in,
                           UINT4 refineby,
                           UINT4 Dterms
                           )
{
  UINT4 newLen, oldLen, l;
  COMPLEX8Vector *ret = NULL;

  if ( !in )
    return NULL;

  oldLen = in->length;
  newLen = oldLen * refineby;

  /* the following are used to speed things up in the innermost loop */
  if ( (ret = XLALCreateCOMPLEX8Vector ( newLen )) == NULL )
    return NULL;

  for (l=0; l < newLen; l++)
    {

      REAL8 kappa_l_k;
      REAL8 remain, kstarREAL;
      UINT4 kstar, kmin, kmax, k;
      REAL8 sink, coskm1;
      REAL8 Yk_re, Yk_im, Xd_re, Xd_im;

      kstarREAL = 1.0 * l  / refineby;
      kstar = lround( kstarREAL );	/* round to closest bin */
      kstar = MYMIN ( kstar, oldLen - 1 );	/* stay within the old SFT index-bounds */
      remain = kstarREAL - kstar;

      /* boundaries for innermost loop */
      kmin = MYMAX( 0, (INT4)kstar - (INT4)Dterms );
      kmax = MYMIN( oldLen, kstar + Dterms );

      Yk_re = Yk_im = 0;
      if ( fabs(remain) > 1e-5 )	/* denominater doens't vanish */
	{
	  /* Optimization: sin(2pi*kappa(l,k)) = sin(2pi*kappa(l,0) and idem for cos */
	  sink = sin ( LAL_TWOPI * remain );
	  coskm1 = cos ( LAL_TWOPI * remain ) - 1.0;

	  /* ---------- innermost loop: k over 2*Dterms around kstar ---------- */
	  for (k = kmin; k < kmax; k++)
	    {
	      REAL8 Plk_re, Plk_im;

	      Xd_re = crealf(in->data[k]);
	      Xd_im = cimagf(in->data[k]);

	      kappa_l_k = kstarREAL - k;

	      Plk_re = sink / kappa_l_k;
	      Plk_im = coskm1 / kappa_l_k;

	      Yk_re += Plk_re * Xd_re - Plk_im * Xd_im;
	      Yk_im += Plk_re * Xd_im + Plk_im * Xd_re;

	    } /* hotloop over Dterms */
	}
      else	/* kappa -> 0: Plk = 2pi delta(k, l) */
	{
	  Yk_re = LAL_TWOPI * crealf(in->data[kstar]);
	  Yk_im = LAL_TWOPI * cimagf(in->data[kstar]);
	}

      ret->data[l] = crectf( OOTWOPI* Yk_re, OOTWOPI * Yk_im );

    }  /* for l < newlen */

  return ret;

} /* XLALrefineCOMPLEX8Vector() */
