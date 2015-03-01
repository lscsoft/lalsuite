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

#include <complex.h>
#include <fftw3.h>
#include <lal/FFTWMutex.h>
#include <lal/Units.h>
// ========== Resamp internals ==========

// ----- local macros ----------
#define MYMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MYMIN(x,y) ( (x) < (y) ? (x) : (y) )

// ----- local constants
#define NUM_FACT 7
static const REAL8 inv_fact[NUM_FACT] = { 1.0, 1.0, (1.0/2.0), (1.0/6.0), (1.0/24.0), (1.0/120.0), (1.0/720.0) };

// ----- local types ----------
typedef struct tagMultiUINT4Vector
{
  UINT4 length;
  UINT4Vector **data;
} MultiUINT4Vector;

// ----- workspace ----------
typedef struct tagWorkspace_t
{
  COMPLEX8Vector *TimeSeriesSpinCorr_SRC;		// single-detector SRC-frame spindown-corrected timeseries

  // input padded timeseries ts(t) and output Fab(f) of length 'numSamplesFFT' and corresponding fftw plan
  UINT4 numSamplesFFT;					// keep track of previous padded SRC-frame samples (ie dFreq)
  COMPLEX8 *FabX_Raw;					// raw full-band FFT result Fa,Fb and zero-padded, AM weighted + spindown-corr SRC-frame TS
  fftwf_plan fftplan;					// buffer FFT plan for given numSamplesOut length

  // arrays of size numFreqBinsOut over frequency bins f_k:
  UINT4 numFreqBinsOut;					// number of output frequency bins {f_k}
  COMPLEX8 *normX_k;					// normalization factors turning FabX_Raw into final F^X_{a,b} values
  COMPLEX8 *FaX_k;					// properly normalized F_a^X(f_k) over output bins
  COMPLEX8 *FbX_k;					// properly normalized F_b^X(f_k) over output bins
  COMPLEX8 *Fa_k;					// properly normalized F_a(f_k) over output bins
  COMPLEX8 *Fb_k;					// properly normalized F_b(f_k) over output bins
} Workspace_t;

struct tagFstatInput_Resamp
{
  MultiCOMPLEX8TimeSeries  *multiTimeSeries_DET;	// input SFTs converted into a heterodyned timeseries
  // ----- buffering -----
  PulsarDopplerParams prev_doppler;			// buffering: previous phase-evolution ("doppler") parameters
  MultiAMCoeffs *prev_multiAMcoef;			// buffering: previous AM-coeffs, unique to skypos
  MultiSSBtimes *prev_multiTimingSRC;			// buffering: previous sky+binary multiSSB times

  MultiCOMPLEX8TimeSeries *prev_multiTimeSeries_SRC;	// buffering: multi-detector SRC-frame timeseries
  MultiUINT4Vector *prev_multiSFTinds_SRC;		// buffering: SFT timestamps translated into SRC frame

  Workspace_t ws;					// 'workspace': pre-allocated vectors used to store intermediate results
};


// ----- local prototypes ----------
static int
XLALApplySpindownAndFreqShift ( COMPLEX8Vector *xOut,
                                const COMPLEX8TimeSeries *xIn,
                                const UINT4Vector *SFTinds,
                                const PulsarDopplerParams *doppler,
                                REAL8 freqShift
                                );
static int
XLALApplyAmplitudeModulation ( COMPLEX8 *xOut,
                               const COMPLEX8Vector *xIn,
                               const UINT4Vector *SFTinds,
                               const REAL4Vector *ab
                               );


static int
XLALBarycentricResampleMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries *mTimeSeries_SRC,
                                                 MultiUINT4Vector *mSFTinds_SRC,
                                                 const MultiCOMPLEX8TimeSeries *mTimeSeries_DET,
                                                 const MultiLIGOTimeGPSVector *mTimestamps_DET,
                                                 const MultiSSBtimes *mSRC_timing
                                                 );

static int
XLALBarycentricResampleCOMPLEX8TimeSeries ( COMPLEX8TimeSeries *TimeSeries_SRC,
                                            UINT4Vector *SFTinds_SRC,
                                            const COMPLEX8TimeSeries *TimeSeries_DET,
                                            const LIGOTimeGPSVector *Timestamps_DET,
                                            const SSBtimes *SRC_timing
                                            );


static int
XLALComputeFaFb_Resamp ( Workspace_t *ws,
                         const PulsarDopplerParams thisPoint,
                         REAL8 dFreq,
                         const COMPLEX8TimeSeries *TimeSeries_SRC,
                         const UINT4Vector *SFTinds_SRC,
                         const AMCoeffs *ab
                         );


// ==================== function definitions ====================
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

  // ----- free buffer
  XLALDestroyMultiAMCoeffs ( resamp->prev_multiAMcoef );
  XLALDestroyMultiCOMPLEX8TimeSeries ( resamp->prev_multiTimeSeries_SRC );
  XLALDestroyMultiUINT4Vector ( resamp->prev_multiSFTinds_SRC );

  // ----- free workspace
  XLALDestroyCOMPLEX8Vector ( resamp->ws.TimeSeriesSpinCorr_SRC );

  LAL_FFTW_WISDOM_LOCK;
  fftwf_destroy_plan ( resamp->ws.fftplan );
  LAL_FFTW_WISDOM_UNLOCK;

  fftw_free ( resamp->ws.FabX_Raw );

  XLALFree ( resamp->ws.normX_k );
  XLALFree ( resamp->ws.FaX_k );
  XLALFree ( resamp->ws.FbX_k );
  XLALFree ( resamp->ws.Fa_k );
  XLALFree ( resamp->ws.Fb_k );

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

  UINT4 numDetectors = resamp->multiTimeSeries_DET->length;
  REAL8 dt_DET = resamp->multiTimeSeries_DET->data[0]->deltaT;
  REAL8 fHet = resamp->multiTimeSeries_DET->data[0]->f0;

  // determine resampled timeseries parameters */
  REAL8 TspanFFT = 1.0 / common->dFreq;
  UINT4 numSamplesFFT = (UINT4) ceil ( TspanFFT / dt_DET );      // we use ceil() so that we artificially widen the band rather than reduce it
  REAL8 dt_SRC = TspanFFT / numSamplesFFT;			// adjust sampling rate to allow achieving exact requested dFreq=1/TspanFFT !

  // ----- prepare memory for the SRC-frame resampled timeseries buffer
  XLAL_CHECK ( (resamp->prev_multiTimeSeries_SRC = XLALCalloc ( 1, sizeof(MultiCOMPLEX8TimeSeries)) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( (resamp->prev_multiTimeSeries_SRC->data = XLALCalloc ( numDetectors, sizeof(COMPLEX8TimeSeries) )) != NULL, XLAL_ENOMEM );
  resamp->prev_multiTimeSeries_SRC->length = numDetectors;

  XLAL_CHECK ( (resamp->prev_multiSFTinds_SRC = XLALCalloc ( 1, sizeof(MultiUINT4Vector) )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( (resamp->prev_multiSFTinds_SRC->data = XLALCalloc ( numDetectors, sizeof(UINT4Vector) )) != NULL, XLAL_EFUNC );
  resamp->prev_multiSFTinds_SRC->length = numDetectors;

  UINT4 numSamplesInMax = 0;
  for ( UINT4 X = 0; X < numDetectors; X ++ )
    {
      const char *nameX = resamp->multiTimeSeries_DET->data[X]->name;
      UINT4 numSamplesInX = resamp->multiTimeSeries_DET->data[X]->data->length;
      numSamplesInMax = MYMAX ( numSamplesInMax, numSamplesInX );
      XLAL_CHECK ( dt_DET == resamp->multiTimeSeries_DET->data[X]->deltaT, XLAL_EINVAL, "Input timeseries must have identical 'deltaT'\n");

      LIGOTimeGPS epoch0 = {0,0};	// will be set to corresponding SRC-frame epoch when barycentering
      XLAL_CHECK ( (resamp->prev_multiTimeSeries_SRC->data[X] = XLALCreateCOMPLEX8TimeSeries ( nameX, &epoch0, fHet, dt_SRC, &lalDimensionlessUnit, numSamplesInX )) != NULL, XLAL_EFUNC );

      UINT4 numTimestampsX = common->multiTimestamps->data[X]->length;
      XLAL_CHECK ( (resamp->prev_multiSFTinds_SRC->data[X] = XLALCreateUINT4Vector ( 2 * numTimestampsX )) != NULL, XLAL_EFUNC );
    } // for X < numDetectors

  // ---- prepare (fixed-size) workspace timeseries in SRC frame for spin-corrections
  Workspace_t *ws = &resamp->ws;
  XLAL_CHECK ( ( ws->TimeSeriesSpinCorr_SRC = XLALCreateCOMPLEX8Vector ( numSamplesInMax )) != NULL, XLAL_EFUNC );


  XLAL_CHECK ( (ws->FabX_Raw = fftw_malloc ( numSamplesFFT * sizeof(ws->FabX_Raw[0]) )) != NULL, XLAL_EFUNC );

  LAL_FFTW_WISDOM_LOCK;
  XLAL_CHECK ( (ws->fftplan = fftwf_plan_dft_1d ( numSamplesFFT, ws->FabX_Raw, ws->FabX_Raw, FFTW_FORWARD, FFTW_MEASURE )) != NULL, XLAL_EFAILED, "fftwf_plan_dft_1d() failed\n");
  LAL_FFTW_WISDOM_UNLOCK;

  ws->numSamplesFFT = numSamplesFFT;


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
  const MultiCOMPLEX8TimeSeries *multiTimeSeries_DET = resamp->multiTimeSeries_DET;
  UINT4 numDetectors = multiTimeSeries_DET->length;

  // ============================== BEGIN: handle buffering =============================
  BOOLEAN same_skypos = (resamp->prev_doppler.Alpha == thisPoint.Alpha) && (resamp->prev_doppler.Delta == thisPoint.Delta);
  BOOLEAN same_refTime = ( XLALGPSCmp ( &resamp->prev_doppler.refTime, &thisPoint.refTime ) == 0 );
  BOOLEAN same_binary = \
    (resamp->prev_doppler.asini == thisPoint.asini) &&
    (resamp->prev_doppler.period == thisPoint.period) &&
    (resamp->prev_doppler.ecc == thisPoint.ecc) &&
    (XLALGPSCmp( &resamp->prev_doppler.tp, &thisPoint.tp ) == 0 ) &&
    (resamp->prev_doppler.argp == thisPoint.argp);

  SkyPosition skypos;
  skypos.system = COORDINATESYSTEM_EQUATORIAL;
  skypos.longitude = thisPoint.Alpha;
  skypos.latitude  = thisPoint.Delta;

  Workspace_t *ws = &(resamp->ws);
  // ----- same skyposition? --> reuse antenna-patterns
  MultiAMCoeffs *multiAMcoef;
  if ( same_skypos && (resamp->prev_multiAMcoef != NULL) )
    {
      multiAMcoef = resamp->prev_multiAMcoef;
    }
  else
    {
      XLALDestroyMultiAMCoeffs ( resamp->prev_multiAMcoef );
      XLAL_CHECK ( (multiAMcoef = XLALComputeMultiAMCoeffs ( common->multiDetectorStates, common->multiNoiseWeights, skypos )) != NULL, XLAL_EFUNC );
      resamp->prev_multiAMcoef = multiAMcoef;
    }

  // ----- not same skypos+binary+refTime? --> re-compute SRC-frame timeseries
  if ( ! ( same_skypos && same_refTime && same_binary) )
    {
      MultiSSBtimes *multiTimingSRC;
      XLAL_CHECK ( (multiTimingSRC = XLALGetMultiSSBtimes ( common->multiDetectorStates, skypos, thisPoint.refTime, common->SSBprec )) != NULL, XLAL_EFUNC );
      if ( thisPoint.asini > 0 ) {
        XLAL_CHECK ( XLALAddMultiBinaryTimes ( &multiTimingSRC, multiTimingSRC, &thisPoint ) == XLAL_SUCCESS, XLAL_EFUNC );
      }
      XLAL_CHECK ( XLALBarycentricResampleMultiCOMPLEX8TimeSeries ( resamp->prev_multiTimeSeries_SRC, resamp->prev_multiSFTinds_SRC, multiTimeSeries_DET, common->multiTimestamps, multiTimingSRC )
                   == XLAL_SUCCESS, XLAL_EFUNC );
      XLALDestroyMultiSSBtimes ( multiTimingSRC );
    }
  resamp->prev_doppler = thisPoint;

  MultiCOMPLEX8TimeSeries *multiTimeSeries_SRC = resamp->prev_multiTimeSeries_SRC;
  MultiUINT4Vector *multiSFTinds_SRC = resamp->prev_multiSFTinds_SRC;

  // ============================== check workspace is properly allocated and initialized ===========

  // ----- workspace that depends on number of output frequency bins 'numFreqBins' ----------
  UINT4 numFreqBins = Fstats->numFreqBins;
  ws->numFreqBinsOut = numFreqBins;
  XLAL_CHECK ( (ws->normX_k = XLALRealloc ( ws->normX_k, numFreqBins * sizeof(COMPLEX8) )) != NULL, XLAL_ENOMEM );

  // NOTE: we try to use as much existing memory as possible in FstatResults, so we only
  // allocate local 'workspace' storage in case there's not already a vector allocated in FstatResults for it
  // this also avoid having to copy these results in case the user asked for them to be returned
  if ( whatToCompute & FSTATQ_FAFB )
    {
      XLALFree ( ws->Fa_k ); // avoid memory leak if allocated in previous call
      ws->Fa_k = Fstats->Fa;
      XLALFree ( ws->Fb_k ); // avoid memory leak if allocated in previous call
      ws->Fb_k = Fstats->Fb;
    } // end: if returning FaFb we can use that return-struct as 'workspace'
  else	// otherwise: we (re)allocate it locally
    {
      XLAL_CHECK ( (ws->Fa_k = XLALRealloc ( ws->Fa_k, numFreqBins * sizeof(COMPLEX8))) != NULL, XLAL_ENOMEM );
      XLAL_CHECK ( (ws->Fb_k = XLALRealloc ( ws->Fb_k, numFreqBins * sizeof(COMPLEX8))) != NULL, XLAL_ENOMEM );
    }

  if ( whatToCompute & FSTATQ_FAFB_PER_DET )
    {
      XLALFree ( ws->FaX_k ); // avoid memory leak if allocated in previous call
      ws->FaX_k = NULL;	// will be set in loop over detectors X
      XLALFree ( ws->FbX_k ); // avoid memory leak if allocated in previous call
      ws->FbX_k = NULL;	// will be set in loop over detectors X
    } // end: if returning FaFbPerDet we can use that return-struct as 'workspace'
  else	// otherwise: we (re)allocate it locally
    {
      XLAL_CHECK ( (ws->FaX_k = XLALRealloc ( ws->FaX_k, numFreqBins * sizeof(COMPLEX8))) != NULL, XLAL_ENOMEM );
      XLAL_CHECK ( (ws->FbX_k = XLALRealloc ( ws->FbX_k, numFreqBins * sizeof(COMPLEX8))) != NULL, XLAL_ENOMEM );
    }

  // ====================================================================================================

  // store AM coefficient integrals in local variables
  REAL4 Ad = multiAMcoef->Mmunu.Ad;
  REAL4 Bd = multiAMcoef->Mmunu.Bd;
  REAL4 Cd = multiAMcoef->Mmunu.Cd;
  REAL4 Ed = multiAMcoef->Mmunu.Ed;
  REAL4 Dd = multiAMcoef->Mmunu.Dd;
  REAL4 Dd_inv = 1.0f / Dd;

  memset ( ws->Fa_k, 0, numFreqBins * sizeof(COMPLEX8) );
  memset ( ws->Fb_k, 0, numFreqBins * sizeof(COMPLEX8) );

  // loop over detectors
  for ( UINT4 X=0; X < numDetectors; X++ )
    {
      // if return-struct contains memory for holding FaFbPerDet: use that directly instead of local memory
      if ( whatToCompute & FSTATQ_FAFB_PER_DET )
        {
          ws->FaX_k = Fstats->FaPerDet[X];
          ws->FbX_k = Fstats->FbPerDet[X];
        }
      const UINT4Vector *SFTindsX_SRC = multiSFTinds_SRC->data[X];
      const COMPLEX8TimeSeries *TimeSeriesX_SRC = multiTimeSeries_SRC->data[X];
      const AMCoeffs *abX = multiAMcoef->data[X];

      // compute {Fa^X(f_k), Fb^X(f_k)}: results returned via workspace resamp->ws
      XLAL_CHECK ( XLALComputeFaFb_Resamp ( &(resamp->ws), thisPoint, common->dFreq, TimeSeriesX_SRC, SFTindsX_SRC, abX ) == XLAL_SUCCESS, XLAL_EFUNC );

      for ( UINT4 k = 0; k < numFreqBins; k++ )
        {
          ws->Fa_k[k] += ws->FaX_k[k];
          ws->Fb_k[k] += ws->FbX_k[k];
        }

      // ----- if requested: compute per-detector Fstat_X_k
      if ( whatToCompute & FSTATQ_2F_PER_DET )
        {
          REAL4 AdX = abX->A;
          REAL4 BdX = abX->B;
          REAL4 CdX = abX->C;
          REAL4 EdX = 0; // FIXME
          REAL4 DdX_inv = 1.0 / abX->D;

          for ( UINT4 k = 0; k < numFreqBins; k ++ )
            {
              Fstats->twoFPerDet[X][k] = XLALComputeFstatFromFaFb ( ws->FaX_k[k], ws->FbX_k[k], AdX, BdX, CdX, EdX, DdX_inv );
            }
        } // for k < numFreqBins

    } // for X < numDetectors

  if ( whatToCompute & FSTATQ_2F )
    {
      for ( UINT4 k=0; k < numFreqBins; k++ )
        {
          Fstats->twoF[k] = XLALComputeFstatFromFaFb ( ws->Fa_k[k], ws->Fb_k[k], Ad, Bd, Cd, Ed, Dd_inv );
        } // for k < numFreqBins
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
      ws->Fa_k = NULL;
      ws->Fb_k = NULL;
    }
  if ( whatToCompute & FSTATQ_FAFB_PER_DET )
    {
      ws->FaX_k = NULL;
      ws->FbX_k = NULL;
    }

  return XLAL_SUCCESS;

} // ComputeFstat_Resamp()


static int
XLALComputeFaFb_Resamp ( Workspace_t *ws,				//!< [in,out] pre-allocated 'workspace' for temporary and output quantities
                         const PulsarDopplerParams thisPoint,		//!< [in] Doppler point to compute {FaX,FbX} for
                         REAL8 dFreq,					//!< [in] output frequency resolution
                         const COMPLEX8TimeSeries *TimeSeries_SRC,	//!< [in] SRC-frame single-IFO timeseries
                         const UINT4Vector *SFTinds_SRC,		//!< [in] corresponding SFT start/end indices
                         const AMCoeffs *ab				//!< [in] antenna-pattern coefficients aX[k], bX[k]
                         )
{
  XLAL_CHECK ( (ws != NULL) && (TimeSeries_SRC != NULL) && (SFTinds_SRC != NULL) && (ab != NULL), XLAL_EINVAL );
  XLAL_CHECK ( dFreq > 0, XLAL_EINVAL );

  // check workspace properly setup
  XLAL_CHECK ( ws->TimeSeriesSpinCorr_SRC != NULL, XLAL_EINVAL );
  XLAL_CHECK ( ws->FabX_Raw != NULL, XLAL_EINVAL );
  XLAL_CHECK ( ws->normX_k != NULL, XLAL_EINVAL );
  XLAL_CHECK ( ws->FaX_k != NULL, XLAL_EINVAL );
  XLAL_CHECK ( ws->FbX_k != NULL, XLAL_EINVAL );
  XLAL_CHECK ( ws->fftplan != NULL, XLAL_EINVAL );


  UINT4 numSFTs = SFTinds_SRC->length / 2;
  XLAL_CHECK ( (ab->a->length == numSFTs) && (ab->b->length == numSFTs ), XLAL_EINVAL );

  REAL8 FreqOut0 = thisPoint.fkdot[0];

  // compute frequency shift to align heterodyne frequency with output frequency bins
  REAL8 fHet   = TimeSeries_SRC->f0;
  REAL8 dt_SRC = TimeSeries_SRC->deltaT;

  REAL8 freqShift = remainder ( FreqOut0 - fHet, dFreq ); // frequency shift to closest bin

  // lowest FFT frequency bin (after applying freqShift)
  UINT4 NnegBins = NhalfNeg ( ws->numSamplesFFT );
  UINT4 NposBinsDC = ws->numSamplesFFT - NnegBins;
  REAL8 fMinFFT = fHet + freqShift - dFreq * NnegBins;
  UINT4 offset_bins = (UINT4) lround ( ( FreqOut0 - fMinFFT ) / dFreq );

  // apply spindown phase-factors, store result in 'workspace'
  memset ( ws->TimeSeriesSpinCorr_SRC->data, 0, ws->TimeSeriesSpinCorr_SRC->length * sizeof(COMPLEX8));
  XLAL_CHECK ( XLALApplySpindownAndFreqShift ( ws->TimeSeriesSpinCorr_SRC, TimeSeries_SRC, SFTinds_SRC, &thisPoint, freqShift ) == XLAL_SUCCESS, XLAL_EFUNC );

  // ----- compute normalization factors to be applied to Fa and Fb:
  const REAL8 dtauX = XLALGPSDiff ( &TimeSeries_SRC->epoch, &thisPoint.refTime );
  for ( UINT4 k = 0; k < ws->numFreqBinsOut; k++ )
    {
      REAL8 f_k = FreqOut0 + k * dFreq;
      REAL8 cycles = - f_k * dtauX;
      REAL4 sinphase, cosphase;
      XLALSinCos2PiLUT ( &sinphase, &cosphase, cycles );
      ws->normX_k[k] = dt_SRC * crectf ( cosphase, sinphase );
    } // for k < numFreqBinsOut

  // ----- compute FaX_k
  // apply amplitude modulation factors {a,b}, store result in zero-padded timeseries for FFTing
  memset ( ws->FabX_Raw, 0, ws->numSamplesFFT * sizeof(ws->FabX_Raw[0]) );
  XLAL_CHECK ( XLALApplyAmplitudeModulation ( ws->FabX_Raw, ws->TimeSeriesSpinCorr_SRC, SFTinds_SRC, ab->a ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Fourier transform the resampled Fa(t)
  fftwf_execute ( ws->fftplan );

  for ( UINT4 k = 0; k < ws->numFreqBinsOut; k++ )
    {
      UINT4 idy = k + offset_bins + NposBinsDC;
      UINT4 idyFFT = idy % ws->numSamplesFFT;	// physical access pattern in FFT bin ordering!

      ws->FaX_k[k] = ws->normX_k[k] * ws->FabX_Raw[idyFFT];

    } // for k < numFreqBinsOut

  // ----- compute FbX_k
  // apply amplitude modulation factors {a,b}, store result in zero-padded timeseries for FFTing
  memset ( ws->FabX_Raw, 0, ws->numSamplesFFT * sizeof(ws->FabX_Raw[0]) );
  XLAL_CHECK ( XLALApplyAmplitudeModulation ( ws->FabX_Raw, ws->TimeSeriesSpinCorr_SRC, SFTinds_SRC, ab->b ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Fourier transform the resampled Fa(t)
  fftwf_execute ( ws->fftplan );

  for ( UINT4 k = 0; k < ws->numFreqBinsOut; k++ )
    {
      UINT4 idy = k + offset_bins + NposBinsDC;
      UINT4 idyFFT = idy % ws->numSamplesFFT;	// physical access pattern in FFT bin ordering!

      ws->FbX_k[k] = ws->normX_k[k] * ws->FabX_Raw[idyFFT];

    } // for k < numFreqBinsOut

  return XLAL_SUCCESS;

} // XLALComputeFaFb_Resamp()

static int
XLALApplySpindownAndFreqShift ( COMPLEX8Vector *xOut,      		///< [out] the spindown-corrected SRC-frame timeseries
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
  UINT4 numSamplesOut = xOut->length;
  XLAL_CHECK ( numSamplesOut >= numSamplesIn, XLAL_EINVAL );

  const LIGOTimeGPS *epoch = &(xIn->epoch);
  REAL8 Dtau0 = XLALGPSDiff ( epoch, &(doppler->refTime) );

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
          xOut->data[j] = em2piphase * xIn->data->data[j];

        } // for j in [start_index, end_index]

    } // for alpha < numSFTs

  return XLAL_SUCCESS;

} // XLALApplySpindownAndFreqShift()

static int
XLALApplyAmplitudeModulation ( COMPLEX8 *xOut,      		///< [out] the spindown-corrected SRC-frame timeseries of length >= length(xIn)
                               const COMPLEX8Vector *xIn,	///< [in] the input SRC-frame timeseries
                               const UINT4Vector *SFTinds,	///< [in] SFT start- and stop indices in the TimeSeries
                               const REAL4Vector *ab		///< [in] amplitude-modulation factors to apply
                               )
{
  // input sanity checks
  XLAL_CHECK ( xOut != NULL, XLAL_EINVAL );
  XLAL_CHECK ( xIn != NULL, XLAL_EINVAL );
  XLAL_CHECK ( SFTinds != NULL, XLAL_EINVAL );
  XLAL_CHECK ( ab != NULL, XLAL_EINVAL );

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
          xOut[j] = fact * xIn->data[j];
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
                                                 const MultiSSBtimes *mSRC_timing			///< [in] multi-detector SRC timing data (time offsets+derivatives)
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
      XLAL_CHECK ( XLALBarycentricResampleCOMPLEX8TimeSeries ( TimeSeries_SRCX, SFTinds_SRCX, TimeSeries_DETX, Timestamps_DETX, SRCtimingX ) == XLAL_SUCCESS, XLAL_EFUNC );
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
XLALBarycentricResampleCOMPLEX8TimeSeries ( COMPLEX8TimeSeries *TimeSeries_SRC,		///< [in,out] resampled timeseries in the source (SRC) frame x(t(t_SRC)), must be alloced+initialized correctly!
                                            UINT4Vector *SFTinds_SRC,			///< [out] start- and end- SFT times in SRC frame, expressed as indices in the SRC timeseries
                                            const COMPLEX8TimeSeries *TimeSeries_DET,	///< [in] the input detector-frame timeseries x(t)
                                            const LIGOTimeGPSVector *Timestamps_DET,	///< [in] the SFT timestamps in the detector frame
                                            const SSBtimes *SRC_timing			///< [in] the source-frame time-shifts and time-derivatives at the SFT midpoints
                                            )
{
  // check input sanity
  XLAL_CHECK ( (TimeSeries_DET != NULL) && (TimeSeries_DET->data != NULL), XLAL_EINVAL );
  XLAL_CHECK ( (TimeSeries_SRC != NULL) && (TimeSeries_SRC->data != NULL), XLAL_EINVAL );
  XLAL_CHECK ( (Timestamps_DET != NULL) && (Timestamps_DET->data != NULL), XLAL_EINVAL );
  XLAL_CHECK ( (SFTinds_SRC != NULL) && (SFTinds_SRC->data != NULL), XLAL_EINVAL );
  XLAL_CHECK ( (SRC_timing != NULL) && (SRC_timing->DeltaT != NULL) && (SRC_timing->Tdot != NULL), XLAL_EINVAL );

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
  REAL8 dt_SRC         = TimeSeries_SRC->deltaT;
  // determine and set time-series start-time in SRC frame
  REAL8 start_SRC   = refTime + SRC_timing->DeltaT->data[0] - (0.5*Tsft) * SRC_timing->Tdot->data[0];
  XLALGPSSetREAL8 ( &(TimeSeries_SRC->epoch), start_SRC );

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
