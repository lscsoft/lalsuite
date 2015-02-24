//
// Copyright (C) 2014 Reinhard Prix
// Copyright (C) 2012, 2013, 2014 Karl Wette
// Copyright (C) 2009 Chris Messenger, Reinhard Prix, Pinkesh Patel, Xavier Siemens, Holger Pletsch
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA  02111-1307  USA
//

// This file implements the F-statistic resampling algorithm. It is not compiled directly, but
// included from ComputeFstat.c

// ========== Resamp internals ==========

// ----- local macros ----------
#define MYMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MYMIN(x,y) ( (x) < (y) ? (x) : (y) )

// ----- local constants
static LALUnit emptyLALUnit;
#define NUM_FACT 7
static const REAL8 inv_fact[NUM_FACT] = { 1.0, 1.0, (1.0/2.0), (1.0/6.0), (1.0/24.0), (1.0/120.0), (1.0/720.0) };

// ----- local types ----------
typedef struct tagMultiUINT4Vector
{
  UINT4 length;
  UINT4Vector **data;
} MultiUINT4Vector;

struct tagFstatInput_Resamp {
  MultiCOMPLEX8TimeSeries  *multiTimeSeries_DET;	// input SFTs converted into a heterodyned timeseries
  // ----- buffering -----
  UINT4 prev_numSamplesOut;				// keep track of previous SRC-frame samples (ie dFreqOut)

  PulsarDopplerParams prev_doppler;			// buffering: previous phase-evolution ("doppler") parameters
  MultiAMCoeffs *prev_multiAMcoef;			// buffering: previous AM-coeffs, unique to skypos
  MultiSSBtimes *prev_multiTimingSRC;			// buffering: previous sky+binary multiSSB times

  MultiCOMPLEX8TimeSeries *prev_multiTimeSeries_SRC;	// buffering: multi-detector SRC-frame timeseries
  MultiUINT4Vector *prev_multiSFTinds_SRC;		// buffering: SFT timestamps translated into SRC frame

  // ----- workspace ----------
  COMPLEX8TimeSeries *ws_TimeSeries_SRC;		// workspace holding single-detector SRC-frame spindown-corrected timeseries
  COMPLEX8TimeSeries *ws_TimeSeriesPadded_SRC;		// workspace holding zero-padded + antenna-weighted single-detector SRC-frame spindown timeseries

  COMPLEX8Vector *ws_FabX;				// hold raw results of FTT for unnormalized full-band Fa(f) or Fb(b)
  ComplexFFTPlan *ws_fftplan;

  COMPLEX8 *ws_normX_k;
  COMPLEX8 *ws_FaX_k;
  COMPLEX8 *ws_FbX_k;

  COMPLEX8 *ws_Fa_k;
  COMPLEX8 *ws_Fb_k;
};


// ----- local prototypes ----------
static int
XLALApplySpindownAndFreqShift ( COMPLEX8TimeSeries *xOut,
                                const COMPLEX8TimeSeries *xIn,
                                const UINT4Vector *SFTinds,
                                const PulsarDopplerParams *doppler,
                                REAL8 freqShift
                                );
static int
XLALApplyAmplitudeModulation ( COMPLEX8TimeSeries *xOut,
                               const COMPLEX8TimeSeries *xIn,
                               const UINT4Vector *SFTinds,
                               const REAL4Vector *ab
                               );


static int
XLALBarycentricResampleMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries *mTimeSeries_SRC,
                                                 MultiUINT4Vector *mSFTinds_SRC,
                                                 const MultiCOMPLEX8TimeSeries *mTimeSeries_DET,
                                                 const MultiLIGOTimeGPSVector *mTimestamps_DET,
                                                 const MultiSSBtimes *mSRC_timing,
                                                 REAL8 dt_SRC
                                                 );

static int
XLALBarycentricResampleCOMPLEX8TimeSeries ( COMPLEX8TimeSeries *TimeSeries_SRC,
                                            UINT4Vector *SFTinds_SRC,
                                            const COMPLEX8TimeSeries *TimeSeries_DET,
                                            const LIGOTimeGPSVector *Timestamps_DET,
                                            const SSBtimes *SRC_timing,
                                            REAL8 dt_SRC
                                            );

// ==================== function definitions
static void
XLALDestroyMultiUINT4Vector ( MultiUINT4Vector *v)
{
  if ( v == NULL ) {
    return;
  }
  for ( UINT4 X = 0; X < v->length; X ++ ) {
    XLALDestroyUINT4Vector ( v->data[X] );
  }
  XLALFree ( v->data );
  XLALFree ( v );

  return;

} // XLALDestroyMultiUINT4Vector()

static void
DestroyFstatInput_Resamp ( FstatInput_Resamp* resamp )
{
  XLALDestroyMultiCOMPLEX8TimeSeries (resamp->multiTimeSeries_DET );

  // ----- free workspace
  XLALDestroyCOMPLEX8TimeSeries ( resamp->ws_TimeSeries_SRC );
  XLALDestroyCOMPLEX8TimeSeries ( resamp->ws_TimeSeriesPadded_SRC );

  XLALDestroyCOMPLEX8Vector ( resamp->ws_FabX );
  XLALDestroyCOMPLEX8FFTPlan ( resamp->ws_fftplan );

  XLALFree ( resamp->ws_normX_k );
  XLALFree ( resamp->ws_Fa_k );
  XLALFree ( resamp->ws_Fb_k );
  XLALFree ( resamp->ws_FaX_k );
  XLALFree ( resamp->ws_FbX_k );

  // ----- free buffer
  XLALDestroyMultiAMCoeffs ( resamp->prev_multiAMcoef );
  XLALDestroyMultiSSBtimes ( resamp->prev_multiTimingSRC );
  XLALDestroyMultiCOMPLEX8TimeSeries ( resamp->prev_multiTimeSeries_SRC );
  XLALDestroyMultiUINT4Vector ( resamp->prev_multiSFTinds_SRC );

  XLALFree ( resamp );

  return;
} // DestroyFstatInput_Resamp()

static int
SetupFstatInput_Resamp ( FstatInput_Resamp *resamp,
                         const FstatInput_Common *common,
                         MultiSFTVector *multiSFTs
                         )
{
  // Check input
  XLAL_CHECK(common != NULL, XLAL_EFAULT);
  XLAL_CHECK(resamp != NULL, XLAL_EFAULT);
  XLAL_CHECK(multiSFTs != NULL, XLAL_EFAULT);

  // Convert SFTs into heterodyned complex timeseries [in detector frame]
  XLAL_CHECK ( (resamp->multiTimeSeries_DET = XLALMultiSFTVectorToCOMPLEX8TimeSeries ( multiSFTs )) != NULL, XLAL_EFUNC );

  XLALDestroyMultiSFTVector ( multiSFTs );	// don't need them SFTs any more ...

  // ----- prepare memory for the SRC-frame resampled timeseries buffer
  UINT4 numDetectors = resamp->multiTimeSeries_DET->length;

  XLAL_CHECK ( (resamp->prev_multiTimeSeries_SRC = XLALCalloc ( 1, sizeof(MultiCOMPLEX8TimeSeries)) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( (resamp->prev_multiTimeSeries_SRC->data = XLALCalloc ( numDetectors, sizeof(COMPLEX8TimeSeries) )) != NULL, XLAL_ENOMEM );
  resamp->prev_multiTimeSeries_SRC->length = numDetectors;

  XLAL_CHECK ( (resamp->prev_multiSFTinds_SRC = XLALCalloc ( 1, sizeof(MultiUINT4Vector) )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( (resamp->prev_multiSFTinds_SRC->data = XLALCalloc ( numDetectors, sizeof(UINT4Vector) )) != NULL, XLAL_EFUNC );
  resamp->prev_multiSFTinds_SRC->length = numDetectors;

  UINT4 numSamplesInMax = 0;
  LIGOTimeGPS XLAL_INIT_DECL(epoch0);	// will be set to corresponding SRC-frame epoch when barycentering
  const REAL8 fHet0 = 0;
  const REAL8 dt0 = 0;
  for ( UINT4 X = 0; X < numDetectors; X ++ )
    {
      UINT4 numSamplesInX = resamp->multiTimeSeries_DET->data[X]->data->length;
      numSamplesInMax = MYMAX ( numSamplesInMax, numSamplesInX );
      XLAL_CHECK ( (resamp->prev_multiTimeSeries_SRC->data[X] = XLALCreateCOMPLEX8TimeSeries ( "", &epoch0, fHet0, dt0, &emptyLALUnit, numSamplesInX )) != NULL, XLAL_EFUNC );

      UINT4 numTimestamps = common->timestamps->data[X]->length;
      XLAL_CHECK ( (resamp->prev_multiSFTinds_SRC->data[X] = XLALCreateUINT4Vector ( 2 * numTimestamps )) != NULL, XLAL_EFUNC );
    } // for X < numDetectors

  // ---- prepare (fixed-size) workspace timeseries
  XLAL_CHECK ( (resamp->ws_TimeSeries_SRC = XLALCreateCOMPLEX8TimeSeries ( "", &epoch0, fHet0, dt0, &emptyLALUnit, numSamplesInMax )) != NULL, XLAL_EFUNC );

  return XLAL_SUCCESS;

} // SetupFstatInput_Resamp()


static int
GetFstatExtraBins_Resamp ( FstatInput_Resamp* resamp )
{
  XLAL_CHECK(resamp != NULL, XLAL_EFAULT);
  return 8;	// use 8 extra bins to give better agreement with LALDemod(w Dterms=8) near the boundaries
} // GetFstatExtraBins_Resamp()


static int
ComputeFstat_Resamp ( FstatResults* Fstats,
                      const FstatInput_Common *common,
                      FstatInput_Resamp* resamp
                      )
{
  // Check input
  XLAL_CHECK ( Fstats != NULL, XLAL_EFAULT );
  XLAL_CHECK ( common != NULL, XLAL_EFAULT );
  XLAL_CHECK ( resamp != NULL, XLAL_EFAULT );

  const FstatQuantities whatToCompute = Fstats->whatWasComputed;
  XLAL_CHECK ( !(whatToCompute & FSTATQ_ATOMS_PER_DET), XLAL_EINVAL, "Resampling does not currently support atoms per detector" );

  // ----- handy shortcuts ----------
  PulsarDopplerParams thisPoint = Fstats->doppler;
  REAL8 FreqOut0 = thisPoint.fkdot[0];
  const MultiCOMPLEX8TimeSeries *multiTimeSeries_DET = resamp->multiTimeSeries_DET;
  UINT4 numDetectors = multiTimeSeries_DET->length;
  REAL8 dt_DET = multiTimeSeries_DET->data[0]->deltaT;

  MultiAMCoeffs *multiAMcoef;

  // determine resampled timeseries parameters */
  UINT4 numFreqBinsOut = Fstats->numFreqBins;
  REAL8 dFreqOut = ( Fstats->dFreq > 0 ) ? Fstats->dFreq : 1.0 / (multiTimeSeries_DET->data[0]->data->length * dt_DET);
  REAL8 TspanOut = 1.0 / dFreqOut;
  UINT4 numSamplesOut = (UINT4) ceil ( TspanOut / dt_DET );      // we use ceil() so that we artificially widen the band rather than reduce it
  REAL8 dt_SRC = TspanOut / numSamplesOut;			// adjust sampling rate to allow achieving exact requested dFreqOut=1/TspanOut !

  // ============================== BEGIN: handle buffering =============================
  BOOLEAN same_skypos = (resamp->prev_doppler.Alpha == thisPoint.Alpha) && (resamp->prev_doppler.Delta == thisPoint.Delta);
  BOOLEAN same_refTime = ( XLALGPSDiff ( &resamp->prev_doppler.refTime, &thisPoint.refTime ) == 0 );
  BOOLEAN same_binary = (resamp->prev_doppler.asini == thisPoint.asini) &&
    (resamp->prev_doppler.period == thisPoint.period) &&
    (resamp->prev_doppler.ecc == thisPoint.ecc) &&
    (XLALGPSCmp( &resamp->prev_doppler.tp, &thisPoint.tp ) == 0 ) &&
    (resamp->prev_doppler.argp == thisPoint.argp);
  BOOLEAN same_numSamplesOut = ( resamp->prev_numSamplesOut == numSamplesOut );

  SkyPosition skypos;
  skypos.system = COORDINATESYSTEM_EQUATORIAL;
  skypos.longitude = thisPoint.Alpha;
  skypos.latitude  = thisPoint.Delta;

  // ----- same skyposition? --> reuse antenna-patterns
  if ( same_skypos && (resamp->prev_multiAMcoef != NULL) )
    {
      multiAMcoef = resamp->prev_multiAMcoef;
    }
  else
    {
      XLAL_CHECK ( (multiAMcoef = XLALComputeMultiAMCoeffs ( common->detectorStates, common->noiseWeights, skypos )) != NULL, XLAL_EFUNC );
      XLALDestroyMultiAMCoeffs ( resamp->prev_multiAMcoef );
      resamp->prev_multiAMcoef = multiAMcoef;
    }

  // ----- same skypos+binary+refTime? --> reuse SRC timings
  MultiSSBtimes *multiTimingSRC = NULL;
  if ( same_skypos && same_refTime && same_binary )
    {
      multiTimingSRC = resamp->prev_multiTimingSRC;
    }
  else
    {
      XLAL_CHECK ( (multiTimingSRC = XLALGetMultiSSBtimes ( common->detectorStates, skypos, thisPoint.refTime, common->SSBprec )) != NULL, XLAL_EFUNC );
      if ( thisPoint.asini > 0 ) {
        XLAL_CHECK ( XLALAddMultiBinaryTimes ( &multiTimingSRC, multiTimingSRC, &thisPoint ) == XLAL_SUCCESS, XLAL_EFUNC );
      }
      XLALDestroyMultiSSBtimes ( resamp->prev_multiTimingSRC );
      resamp->prev_multiTimingSRC = multiTimingSRC;
    }
  resamp->prev_doppler = thisPoint;

  // ----- if NOT same SRC timing OR same frequency-resolution (ie numSamplesOut)? --> recompute SRC-frame timeseries
  if (  ! ( same_numSamplesOut && same_skypos && same_refTime && same_binary ) )
    {
      XLAL_CHECK ( XLALBarycentricResampleMultiCOMPLEX8TimeSeries ( resamp->prev_multiTimeSeries_SRC, resamp->prev_multiSFTinds_SRC, multiTimeSeries_DET, common->timestamps, multiTimingSRC, dt_SRC )
                   == XLAL_SUCCESS, XLAL_EFUNC );
    }

  MultiCOMPLEX8TimeSeries *multiTimeSeries_SRC = resamp->prev_multiTimeSeries_SRC;
  MultiUINT4Vector *multiSFTinds_SRC = resamp->prev_multiSFTinds_SRC;

  // ============================== check workspace is properly allocated and initialized ===========
  // ----- workspace that depends on SRC-frame time samples 'numSamplesOut' ----------
  if ( !same_numSamplesOut )
    {
      XLALDestroyCOMPLEX8TimeSeries ( resamp->ws_TimeSeriesPadded_SRC );
      LIGOTimeGPS XLAL_INIT_DECL(epoch0); const REAL8 fHet0 = 0; const REAL8 dt0 = 0;
      XLAL_CHECK ( (resamp->ws_TimeSeriesPadded_SRC = XLALCreateCOMPLEX8TimeSeries ( "", &epoch0, fHet0, dt0, &emptyLALUnit, numSamplesOut )) != NULL, XLAL_EFUNC );

      XLALDestroyCOMPLEX8Vector ( resamp->ws_FabX );
      XLAL_CHECK ( (resamp->ws_FabX = XLALCreateCOMPLEX8Vector ( numSamplesOut )) != NULL, XLAL_EFUNC );

      XLALDestroyCOMPLEX8FFTPlan ( resamp->ws_fftplan );
      XLAL_CHECK ( (resamp->ws_fftplan = XLALCreateCOMPLEX8FFTPlan ( numSamplesOut, 1, 0) ) != NULL, XLAL_EFUNC );

      resamp->prev_numSamplesOut = numSamplesOut;
    } // if number of SRC samples has changed

  // ----- workspace that depends on number of output frequency bins 'numFreqBinsOut' ----------
  XLAL_CHECK ( (resamp->ws_normX_k = XLALRealloc ( resamp->ws_normX_k, numFreqBinsOut * sizeof(COMPLEX8))) != NULL, XLAL_ENOMEM );

  // NOTE: we try to use as much existing memory as possible in FstatResults, so we only
  // allocate local 'workspace' storage in case there's not already a vector allocated in FstatResults for it
  // this also avoid having to copy these results in case the user asked for them to be returned
  if ( whatToCompute & FSTATQ_FAFB )
    {
      if ( resamp->ws_Fa_k != NULL ) {
        XLALFree ( resamp->ws_Fa_k ); // avoid memory leak if allocated in previous call
      }
      resamp->ws_Fa_k = Fstats->Fa;
      if ( resamp->ws_Fb_k != NULL ) {
        XLALFree ( resamp->ws_Fb_k ); // avoid memory leak if allocated in previous call
      }
      resamp->ws_Fb_k = Fstats->Fb;
    } // end: if returning FaFb we can use that return-struct as 'workspace'
  else	// otherwise: we (re)allocate it locally
    {
      XLAL_CHECK ( (resamp->ws_Fa_k = XLALRealloc ( resamp->ws_Fa_k, numFreqBinsOut * sizeof(COMPLEX8))) != NULL, XLAL_ENOMEM );
      XLAL_CHECK ( (resamp->ws_Fb_k = XLALRealloc ( resamp->ws_Fb_k, numFreqBinsOut * sizeof(COMPLEX8))) != NULL, XLAL_ENOMEM );
    }

  if ( whatToCompute & FSTATQ_FAFB_PER_DET )
    {
      if ( resamp->ws_FaX_k != NULL ) {
        XLALFree ( resamp->ws_FaX_k ); // avoid memory leak if allocated in previous call
        resamp->ws_FaX_k = NULL;	// will be set in loop over detectors X
      }
      if ( resamp->ws_FbX_k != NULL ) {
        XLALFree ( resamp->ws_FbX_k ); // avoid memory leak if allocated in previous call
        resamp->ws_FbX_k = NULL;	// will be set in loop over detectors X
      }
    } // end: if returning FaFbPerDet we can use that return-struct as 'workspace'
  else	// otherwise: we (re)allocate it locally
    {
      XLAL_CHECK ( (resamp->ws_FaX_k = XLALRealloc ( resamp->ws_FaX_k, numFreqBinsOut * sizeof(COMPLEX8))) != NULL, XLAL_ENOMEM );
      XLAL_CHECK ( (resamp->ws_FbX_k = XLALRealloc ( resamp->ws_FbX_k, numFreqBinsOut * sizeof(COMPLEX8))) != NULL, XLAL_ENOMEM );
    }

  // ====================================================================================================

  // store AM coefficient integrals in local variables
  REAL4 Ad = multiAMcoef->Mmunu.Ad;
  REAL4 Bd = multiAMcoef->Mmunu.Bd;
  REAL4 Cd = multiAMcoef->Mmunu.Cd;
  REAL4 Ed = multiAMcoef->Mmunu.Ed;
  REAL4 Dd = multiAMcoef->Mmunu.Dd;
  REAL4 Dd_inv = 1.0f / Dd;

  // compute frequency shift to align heterodyne frequency with output frequency bins
  REAL8 fHet = multiTimeSeries_SRC->data[0]->f0;
  REAL8 freqShift = remainder ( FreqOut0 - fHet, dFreqOut ); // frequency shift to closest bin

  // lowest FFT frequency bin (after applying freqShift)
  UINT4 NnegBins = NhalfNeg ( numSamplesOut );
  UINT4 NposBinsDC = numSamplesOut - NnegBins;
  REAL8 fMinFFT = fHet + freqShift - dFreqOut * NnegBins;
  UINT4 offset_bins = (UINT4) lround ( ( FreqOut0 - fMinFFT ) / dFreqOut );

  memset ( resamp->ws_Fa_k, 0, numFreqBinsOut * sizeof(resamp->ws_Fa_k[0]) );
  memset ( resamp->ws_Fb_k, 0, numFreqBinsOut * sizeof(resamp->ws_Fb_k[0]) );

  // loop over detectors
  for ( UINT4 X=0; X < numDetectors; X++ )
    {
      const UINT4Vector *SFTinds_SRC = multiSFTinds_SRC->data[X];
      const COMPLEX8TimeSeries *TimeSeries_SRC = multiTimeSeries_SRC->data[X];
      const AMCoeffs *abX = multiAMcoef->data[X];

      // if return-struct contains memory for holding FaFbPerDet: use that directly instead of local memory
      if ( whatToCompute & FSTATQ_FAFB_PER_DET )
        {
          resamp->ws_FaX_k = Fstats->FaPerDet[X];
          resamp->ws_FbX_k = Fstats->FbPerDet[X];
        }

      // apply spindown phase-factors, store result in 'workspace'
      memset ( resamp->ws_TimeSeries_SRC->data->data, 0, resamp->ws_TimeSeries_SRC->data->length * sizeof(resamp->ws_TimeSeries_SRC->data->data[0]));
      XLAL_CHECK ( XLALApplySpindownAndFreqShift ( resamp->ws_TimeSeries_SRC, TimeSeries_SRC, SFTinds_SRC, &thisPoint, freqShift ) == XLAL_SUCCESS, XLAL_EFUNC );

      // ----- compute normalization factors to be applied to Fa and Fb:
      const REAL8 dtauX = XLALGPSDiff ( &TimeSeries_SRC->epoch, &thisPoint.refTime );
      for ( UINT4 k = 0; k < numFreqBinsOut; k++ )
        {
          REAL8 f_k = FreqOut0 + k * dFreqOut;
          REAL8 cycles = - f_k * dtauX;
          REAL4 sinphase, cosphase;
          XLALSinCos2PiLUT ( &sinphase, &cosphase, cycles );
          resamp->ws_normX_k[k] = dt_SRC * crectf ( cosphase, sinphase );
        } // for k < numFreqBinsOut

      // ----- compute FaX_k
      // apply amplitude modulation factors {a,b}, store result in zero-padded timeseries for FFTing
      memset ( resamp->ws_TimeSeriesPadded_SRC->data->data, 0, resamp->ws_TimeSeriesPadded_SRC->data->length * sizeof(resamp->ws_TimeSeriesPadded_SRC->data->data[0]));

      XLAL_CHECK ( XLALApplyAmplitudeModulation ( resamp->ws_TimeSeriesPadded_SRC, resamp->ws_TimeSeries_SRC, SFTinds_SRC, abX->a ) == XLAL_SUCCESS, XLAL_EFUNC );

      // Fourier transform the resampled Fa(t)
      XLAL_CHECK ( XLALCOMPLEX8VectorFFT ( resamp->ws_FabX, resamp->ws_TimeSeriesPadded_SRC->data, resamp->ws_fftplan ) == XLAL_SUCCESS, XLAL_EFUNC );

      for ( UINT4 k = 0; k < numFreqBinsOut; k++ )
        {
          UINT4 idy = k + offset_bins + NposBinsDC;
          UINT4 idyFFT = idy % numSamplesOut;	// physical access pattern in FFT bin ordering!
          COMPLEX8 FaX_k = resamp->ws_normX_k[k] * resamp->ws_FabX->data[idyFFT];

          resamp->ws_FaX_k[k] = FaX_k;
          resamp->ws_Fa_k[k] += FaX_k;

        } // for k < numFreqBinsOut

      // ----- compute FbX_k
      // apply amplitude modulation factors {a,b}, store result in zero-padded timeseries for FFTing
      // NOTE: we don't need to memset-zero the output timeseries here, as we'll write into the exact same bins as for "a * timeseries" above ...
      XLAL_CHECK ( XLALApplyAmplitudeModulation ( resamp->ws_TimeSeriesPadded_SRC, resamp->ws_TimeSeries_SRC, SFTinds_SRC, abX->b ) == XLAL_SUCCESS, XLAL_EFUNC );

      // Fourier transform the resampled Fa(t)
      XLAL_CHECK ( XLALCOMPLEX8VectorFFT ( resamp->ws_FabX, resamp->ws_TimeSeriesPadded_SRC->data, resamp->ws_fftplan ) == XLAL_SUCCESS, XLAL_EFUNC );

      for ( UINT4 k = 0; k < numFreqBinsOut; k++ )
        {
          UINT4 idy = k + offset_bins + NposBinsDC;
          UINT4 idyFFT = idy % numSamplesOut;	// physical access pattern in FFT bin ordering!
          COMPLEX8 FbX_k = resamp->ws_normX_k[k] * resamp->ws_FabX->data[idyFFT];

          resamp->ws_FbX_k[k] = FbX_k;
          resamp->ws_Fb_k[k] += FbX_k;

        } // for k < numFreqBinsOut

      // ----- if requested: compute per-detector Fstat_X_k
      if ( whatToCompute & FSTATQ_2F_PER_DET )
        {
          REAL4 AdX = abX->A;
          REAL4 BdX = abX->B;
          REAL4 CdX = abX->C;
          REAL4 EdX = 0; // FIXME
          REAL4 DdX_inv = 1.0 / abX->D;

          for ( UINT4 k = 0; k < numFreqBinsOut; k ++ )
            {
              Fstats->twoFPerDet[X][k] = XLALComputeFstatFromFaFb ( resamp->ws_FaX_k[k], resamp->ws_FbX_k[k], AdX, BdX, CdX, EdX, DdX_inv );
            }
        } // for k < numFreqBinsOut

    } // for X < numDetectors

  if ( whatToCompute & FSTATQ_2F )
    {
      for ( UINT4 k=0; k < numFreqBinsOut; k++ )
        {
          Fstats->twoF[k] = XLALComputeFstatFromFaFb ( resamp->ws_Fa_k[k], resamp->ws_Fb_k[k], Ad, Bd, Cd, Ed, Dd_inv );
        } // for k < numFreqBinsOut
    } // if FSTATQ_2F

  // Return F-atoms per detector
  if (whatToCompute & FSTATQ_ATOMS_PER_DET) {
    XLAL_ERROR(XLAL_EFAILED, "NOT implemented!");
  }

  Fstats->Mmunu = multiAMcoef->Mmunu;


  // ----- workspace memory management:
  // if we used the return struct directly to store Fa,Fb results,
  // make sure to wipe those pointers to avoid mistakenly considering them as 'local' memory
  // and re-allocing it in another call to this function
  if ( whatToCompute & FSTATQ_FAFB )
    {
      resamp->ws_Fa_k = NULL;
      resamp->ws_Fb_k = NULL;
    }
  if ( whatToCompute & FSTATQ_FAFB_PER_DET )
    {
      resamp->ws_FaX_k = NULL;
      resamp->ws_FbX_k = NULL;
    }

  return XLAL_SUCCESS;

} // ComputeFstat_Resamp()

static int
XLALApplySpindownAndFreqShift ( COMPLEX8TimeSeries *xOut,      		///< [out] the spindown-corrected SRC-frame timeseries
                                const COMPLEX8TimeSeries *xIn,		///< [in] the input SRC-frame timeseries
                                const UINT4Vector *SFTinds,		///< [in] SFT start- and stop indices in the TimeSeries
                                const PulsarDopplerParams *doppler,	///< [in] containing spindown parameters
                                REAL8 freqShift				///< [in] frequency-shift to apply, sign is "new - old"
                                )
{
  // input sanity checks
  XLAL_CHECK ( xOut != NULL, XLAL_EINVAL );
  XLAL_CHECK ( xIn != NULL, XLAL_EINVAL );
  XLAL_CHECK ( SFTinds != NULL, XLAL_EINVAL );
  XLAL_CHECK ( doppler != NULL, XLAL_EINVAL );

  // determine number of spin downs to include
  UINT4 s_max = PULSAR_MAX_SPINS - 1;
  while ( (s_max > 0) && (doppler->fkdot[s_max] == 0) ) {
    s_max --;
  }

  REAL8 dt = xIn->deltaT;
  UINT4 numSamplesIn  = xIn->data->length;
  UINT4 numSamplesOut = xOut->data->length;
  XLAL_CHECK ( numSamplesOut >= numSamplesIn, XLAL_EINVAL );

  const LIGOTimeGPS *epoch = &(xIn->epoch);
  REAL8 Dtau0 = XLALGPSDiff ( epoch, &(doppler->refTime) );

  // copy TS header information
  COMPLEX8Vector *tmp_xOut = xOut->data;
  (*xOut) = (*xIn);
  xOut->data = tmp_xOut;

  UINT4 numSFTs = SFTinds->length / 2;

  // loop over SFTs
  for ( UINT4 alpha=0; alpha < numSFTs; alpha ++ )
    {
      UINT4 start_index = SFTinds->data[2*alpha];
      UINT4 end_index   = SFTinds->data[2*alpha+1];
      XLAL_CHECK ( start_index < numSamplesIn, XLAL_EINVAL );
      XLAL_CHECK ( end_index < numSamplesIn, XLAL_EINVAL );

      // loop over all samples from this SFT
      for ( UINT4 j=start_index; j <= end_index; j ++ )
        {
          REAL8 taup_j = j * dt;
          REAL8 Dtau_alpha_j = Dtau0 + taup_j;

          REAL8 cycles = - freqShift * taup_j;

          REAL8 Dtau_pow_kp1 = Dtau_alpha_j;
          for ( UINT4 k = 1; k <= s_max; k++ )
            {
              Dtau_pow_kp1 *= Dtau_alpha_j;
              cycles += - inv_fact[k+1] * doppler->fkdot[k] * Dtau_pow_kp1;
            } // for k = 1 ... s_max

          REAL4 cosphase, sinphase;
          XLAL_CHECK( XLALSinCos2PiLUT ( &sinphase, &cosphase, cycles ) == XLAL_SUCCESS, XLAL_EFUNC );
          COMPLEX8 em2piphase = crectf ( cosphase, sinphase );

          // weight the complex timeseries by the antenna patterns
          xOut->data->data[j] = em2piphase * xIn->data->data[j];

        } // for j in [start_index, end_index]

    } // for alpha < numSFTs

  // correct output timeseries heterodyne frequency
  xOut->f0 = xIn->f0 + freqShift;

  return XLAL_SUCCESS;

} // XLALApplySpindownAndFreqShift()

static int
XLALApplyAmplitudeModulation ( COMPLEX8TimeSeries *xOut,      		///< [out] the spindown-corrected SRC-frame timeseries
                               const COMPLEX8TimeSeries *xIn,		///< [in] the input SRC-frame timeseries
                               const UINT4Vector *SFTinds,		///< [in] SFT start- and stop indices in the TimeSeries
                               const REAL4Vector *ab			///< [in] amplitude-modulation factors to apply
                               )
{
  // input sanity checks
  XLAL_CHECK ( xOut != NULL, XLAL_EINVAL );
  XLAL_CHECK ( xIn != NULL, XLAL_EINVAL );
  XLAL_CHECK ( SFTinds != NULL, XLAL_EINVAL );
  XLAL_CHECK ( ab != NULL, XLAL_EINVAL );

  UINT4 numSamplesIn  = xIn->data->length;
  UINT4 numSamplesOut = xOut->data->length;
  XLAL_CHECK ( numSamplesOut >= numSamplesIn, XLAL_EINVAL );

  // copy TS header information
  COMPLEX8Vector *tmp_xOut = xOut->data;
  (*xOut) = (*xIn);
  xOut->data = tmp_xOut;

  UINT4 numSFTs = SFTinds->length / 2;
  XLAL_CHECK ( ab->length == numSFTs, XLAL_EINVAL );

  // loop over SFTs
  for ( UINT4 alpha=0; alpha < numSFTs; alpha ++ )
    {
      const COMPLEX8 fact = (COMPLEX8) ab->data[alpha];
      UINT4 start_index = SFTinds->data[2*alpha];
      UINT4 end_index   = SFTinds->data[2*alpha+1];

      // loop over all samples from this SFT
      // and apply amplitude modulation factor to output timeseries
      for ( UINT4 j=start_index; j <= end_index; j ++ )
        {
          xOut->data->data[j] = fact * xIn->data->data[j];
        } // for j in [start_index, end_index]

    } // for alpha < numSFTs

  return XLAL_SUCCESS;

} // XLALApplyAmplitudeModulation()


///
/// Performs barycentric resampling on a multi-detector timeseries
///
int
XLALBarycentricResampleMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries *mTimeSeries_SRC,		///< [out] resampled timeseries in the source (SRC) frame, must be alloced+initialized correctly!
                                                 MultiUINT4Vector *mSFTinds_SRC,			///< [out] start- and end- SFT times in SRC frame, expressed as indices in the SRC timeseries
                                                 const MultiCOMPLEX8TimeSeries *mTimeSeries_DET,	///< [in] detector frame (DET) timeseries
                                                 const MultiLIGOTimeGPSVector *mTimestamps_DET,		///< [in] multi SFT timestamps in the DET frame
                                                 const MultiSSBtimes *mSRC_timing,			///< [in] multi-detector SRC timing data (time offsets+derivatives)
                                                 REAL8 dt_SRC						///< [in] SRC-frame time-step 'dt' to use
                                                 )
{
  // check input sanity
  XLAL_CHECK ( mTimeSeries_SRC != NULL, XLAL_EINVAL );
  XLAL_CHECK ( mSFTinds_SRC != NULL, XLAL_EINVAL );

  XLAL_CHECK ( mTimeSeries_DET != NULL, XLAL_EINVAL );
  XLAL_CHECK ( mTimestamps_DET != NULL, XLAL_EINVAL );
  XLAL_CHECK ( mSRC_timing != NULL, XLAL_EINVAL );

  UINT4 numDetectors = mTimeSeries_DET->length;
  XLAL_CHECK ( numDetectors >0, XLAL_EINVAL );
  XLAL_CHECK ( (mSRC_timing->length == numDetectors) && (mTimestamps_DET->length == numDetectors), XLAL_EINVAL );
  XLAL_CHECK ( (mTimeSeries_SRC->length == numDetectors) && (mSFTinds_SRC->length == numDetectors), XLAL_EINVAL );

  for ( UINT4 X=0; X < numDetectors; X++)
    {
      // shorthand pointers
      SSBtimes *SRCtimingX = mSRC_timing->data[X];
      COMPLEX8TimeSeries *TimeSeries_DETX = mTimeSeries_DET->data[X];
      LIGOTimeGPSVector *Timestamps_DETX = mTimestamps_DET->data[X];

      COMPLEX8TimeSeries *TimeSeries_SRCX = mTimeSeries_SRC->data[X];
      UINT4Vector *SFTinds_SRCX = mSFTinds_SRC->data[X];

      // perform resampling on current detector timeseries */
      XLAL_CHECK ( XLALBarycentricResampleCOMPLEX8TimeSeries ( TimeSeries_SRCX, SFTinds_SRCX, TimeSeries_DETX, Timestamps_DETX, SRCtimingX, dt_SRC ) == XLAL_SUCCESS, XLAL_EFUNC );
    } // for X < numDetectors

  return XLAL_SUCCESS;

} // XLALBarycentricResampleMultiCOMPLEX8TimeSeries()


///
/// Performs barycentric resampling of a timeseries in the detector frame
/// into a uniformly-sampled timeseries at the source frame.
///
/// We expect that the output timeseries has already been allocated correctly,
/// *and* carry the correct start-time epoch for the output! (FIXME!)
///
int
XLALBarycentricResampleCOMPLEX8TimeSeries ( COMPLEX8TimeSeries *TimeSeries_SRC,		///< [out] resampled timeseries in the source (SRC) frame x(t(t_SRC)), must be alloced+initialized correctly!
                                            UINT4Vector *SFTinds_SRC,			///< [out] start- and end- SFT times in SRC frame, expressed as indices in the SRC timeseries
                                            const COMPLEX8TimeSeries *TimeSeries_DET,	///< [in] the input detector-frame timeseries x(t)
                                            const LIGOTimeGPSVector *Timestamps_DET,	///< [in] the SFT timestamps in the detector frame
                                            const SSBtimes *SRC_timing,			///< [in] the source-frame time-shifts and time-derivatives at the SFT midpoints
                                            REAL8 dt_SRC				///< [in] SRC-frame time-step 'dt' to use
                                            )
{
  // check input sanity
  XLAL_CHECK ( (TimeSeries_DET != NULL) && (TimeSeries_DET->data != NULL), XLAL_EINVAL );
  XLAL_CHECK ( (TimeSeries_SRC != NULL) && (TimeSeries_SRC->data != NULL), XLAL_EINVAL );
  XLAL_CHECK ( (Timestamps_DET != NULL) && (Timestamps_DET->data != NULL), XLAL_EINVAL );
  XLAL_CHECK ( (SFTinds_SRC != NULL) && (SFTinds_SRC->data != NULL), XLAL_EINVAL );
  XLAL_CHECK ( (SRC_timing != NULL) && (SRC_timing->DeltaT != NULL) && (SRC_timing->Tdot != NULL), XLAL_EINVAL );
  XLAL_CHECK ( dt_SRC > 0, XLAL_EINVAL );

  UINT4 numSamples_DET = TimeSeries_DET->data->length;
  XLAL_CHECK ( numSamples_DET > 0, XLAL_EINVAL );

  UINT4 numSFTs = Timestamps_DET->length;
  XLAL_CHECK ( (numSFTs > 0) && (SRC_timing->DeltaT->length == numSFTs) && (SRC_timing->Tdot->length == numSFTs), XLAL_EINVAL );
  XLAL_CHECK ( SFTinds_SRC->length == 2 * numSFTs, XLAL_EINVAL );

  // define some useful shorthands
  REAL8 Tsft    = Timestamps_DET->deltaT;
  REAL8 refTime = XLALGPSGetREAL8 ( &SRC_timing->refTime );
  REAL8 fHet    = TimeSeries_DET->f0;

  REAL8 start_DET   = XLALGPSGetREAL8 ( &TimeSeries_DET->epoch );

  UINT4 numSamples_SRC = TimeSeries_SRC->data->length;
  // determine and set time-series start-time in SRC frame
  REAL8 start_SRC   = refTime + SRC_timing->DeltaT->data[0] - (0.5*Tsft) * SRC_timing->Tdot->data[0];
  XLALGPSSetREAL8 ( &(TimeSeries_SRC->epoch), start_SRC );

  TimeSeries_SRC->deltaT = dt_SRC;
  TimeSeries_SRC->f0 = fHet;

  // make sure all output samples are initialized to zero first, in case of gaps
  memset ( TimeSeries_SRC->data->data, 0, numSamples_SRC * sizeof(COMPLEX8) );

  UINT4 maxSFTnumSamples_SRC = ceil ( 1.1 * Tsft / dt_SRC ) + 1;	// guaranteed to be > SFTnumSamples_SRC
  REAL8Vector *detectortimes; // a vector of *non-uniform* time values in the detector frame (used for interpolation)
  XLAL_CHECK ( (detectortimes = XLALCreateREAL8Vector ( maxSFTnumSamples_SRC )) != NULL, XLAL_EFUNC );
  COMPLEX8Vector *ts_SRC;
  XLAL_CHECK ( (ts_SRC = XLALCreateCOMPLEX8Vector ( maxSFTnumSamples_SRC )) != NULL, XLAL_EFUNC );

  // loop over SFT timestamps to compute the detector frame time samples corresponding to uniformly sampled SRC time samples
  for ( UINT4 j=0; j < numSFTs; j++ )
    {
      // define some useful shorthands
      REAL8 Tdot         = SRC_timing->Tdot->data[j];				// the instantaneous time derivitive dt_SRC/dt_DET at the MID-POINT of the SFT
      REAL8 SFTmid_SRC   = refTime + SRC_timing->DeltaT->data[j];		// MID-POINT time of the SFT at the SRC
      REAL8 SFTstart_SRC = SFTmid_SRC - 0.5*Tsft*Tdot;				// START time of the SFT at the SRC
      REAL8 SFTend_SRC   = SFTmid_SRC + 0.5*Tsft*Tdot;				// END time of the SFT at the SRC
      REAL8 SFTstart_DET = XLALGPSGetREAL8 ( &(Timestamps_DET->data[j]) );	// START time of the SFT at the detector
      REAL8 SFTmid_DET   = SFTstart_DET + 0.5*Tsft;				// MID-POINT time of the SFT at the detector

      // indices of first and last SRC-frame sample corresponding to this SFT
      UINT4 SFTidx_start_SRC  = lround ( (SFTstart_SRC - start_SRC) / dt_SRC );	// the index of the resampled timeseries corresponding to the start of the SFT
      UINT4 SFTidx_end_SRC    = lround ( (SFTend_SRC - start_SRC) / dt_SRC );	// the index of the resampled timeseries corresponding to the end of the SFT

      // truncate to actual SRC-frame timeseries
      SFTidx_start_SRC = MYMIN ( SFTidx_start_SRC, numSamples_SRC - 1);
      SFTidx_end_SRC = MYMIN ( SFTidx_end_SRC, numSamples_SRC - 1);
      UINT4 SFTnumSamples_SRC = SFTidx_end_SRC - SFTidx_start_SRC + 1;		// the number of samples in the SRC-frame for this SFT

      XLAL_CHECK ( SFTnumSamples_SRC <= maxSFTnumSamples_SRC, XLAL_EFAILED, "Coding error: maxSFTnumSamples_SRC = %d < %d = SFTnumSamples_SRC in SFT j = %d\n",
                   maxSFTnumSamples_SRC, SFTnumSamples_SRC, j );

      SFTinds_SRC->data[2*j]   = SFTidx_start_SRC;
      SFTinds_SRC->data[2*j+1] = SFTidx_end_SRC;

      // array of *non-uniform* detector time samples for this SFT
      detectortimes->length = SFTnumSamples_SRC;

      // for each time sample in the SRC frame for this SFT we estimate the detector time. */
      // We use a linear approximation expanding around the midpoint of an SFT where
      // t_DET = SFTmid_DET + (t_SRC - SFTmid_SRC)*dt_DET/dt_SRC
      for ( UINT4 k=0; k < SFTnumSamples_SRC; k++ )
        {
          REAL8 t_SRC = start_SRC + ( k + SFTidx_start_SRC ) * dt_SRC;		// the SRC time of the current resampled time sample
          detectortimes->data[k] = SFTmid_DET + ( t_SRC - SFTmid_SRC ) / Tdot;	// the approximated DET time of the current resampled time sample
        } // for k < SFTnumSamples_SRC

      ts_SRC->length = SFTnumSamples_SRC;
      const UINT4 Dterms = 8;
      XLAL_CHECK ( XLALSincInterpolateCOMPLEX8TimeSeries ( ts_SRC, detectortimes, TimeSeries_DET, Dterms ) == XLAL_SUCCESS, XLAL_EFUNC );

      // place these interpolated timeseries into the output
      // and apply correction due to non-zero heterodyne frequency of input
      for ( UINT4 k=0; k < SFTnumSamples_SRC; k++ )
        {
          UINT4 idx = k + SFTidx_start_SRC;                                                                     // the full resampled timeseries index
          if ( idx >= numSamples_SRC ) {	// temporary FIX to avoid writing outside of memory bounds (FIXME!)
            break;
          }
          REAL8 tDiff = idx * dt_SRC - (detectortimes->data[k] - start_DET); 	// tau' - t'(tau')
          REAL8 cycles = fmod ( fHet * tDiff, 1 );                                                          // the accumulated heterodyne cycles

          // use a look-up-table for speed to compute real and imaginary phase
          REAL4 cosphase, sinphase;                                                                         // the real and imaginary parts of the phase correction
          XLAL_CHECK( XLALSinCos2PiLUT ( &sinphase, &cosphase, -cycles ) == XLAL_SUCCESS, XLAL_EFUNC );
          COMPLEX8 ei2piphase = crectf(cosphase,sinphase);
          TimeSeries_SRC->data->data[idx] = ei2piphase * ts_SRC->data[k];
        } // for k < SFTnumSamples_SRC

    } // for j < numSFTs

  // free memory
  XLALDestroyREAL8Vector ( detectortimes );
  XLALDestroyCOMPLEX8Vector ( ts_SRC );

  return XLAL_SUCCESS;

} // XLALBarycentricResampleCOMPLEX8TimeSeries()
