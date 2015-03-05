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
#define COLLECT_TIMING 1
#define NUM_FACT 7
static const REAL8 inv_fact[NUM_FACT] = { 1.0, 1.0, (1.0/2.0), (1.0/6.0), (1.0/24.0), (1.0/120.0), (1.0/720.0) };

// ----- local types ----------
typedef struct tagMultiUINT4Vector
{
  UINT4 length;
  UINT4Vector **data;
} MultiUINT4Vector;

// ----- workspace ----------
typedef struct tagResampTimingInfo
{ // NOTE: all times refer to a single-detector timing case
  REAL8 tauTotal;		// total time spent in ComputeFstat_Resamp()
  REAL8 tauBary;		// time spent in barycentric resampling
  REAL8 tauSpin;		// time spent in spindown+frequency correction
  REAL8 tauAM;			// time spent in applying AM coefficients
  REAL8 tauFFT;			// time spent in FFT
  REAL8 tauNorm;		// time spent normalizing the final Fa,Fb
  REAL8 tauFab2F;		// time to compute Fstat from {Fa,Fb}
  REAL8 tauMem;			// time to realloc and memset-0 arrays
  REAL8 tauSumFabX;		// time to sum_X Fab^X
  REAL8 tauF1Buf;		// Resampling timing 'constant': Fstat time per template per detector for a 'buffered' case (same skypos, same numFreqBins)
  REAL8 tauF1NoBuf;		// Resampling timing 'constant': Fstat time per template per detector for an 'unbuffered' usage (different skypos and numFreqBins)
} ResampTimingInfo;

struct tagFstatWorkspace
{
  COMPLEX8Vector *TimeSeriesSpinCorr_SRC;		// single-detector SRC-frame spindown-corrected timeseries

  // input padded timeseries ts(t) and output Fab(f) of length 'numSamplesFFT' and corresponding fftw plan
  UINT4 numSamplesFFT;					// keep track of previous padded SRC-frame samples (ie dFreq)
  COMPLEX8 *FabX_Raw;					// raw full-band FFT result Fa,Fb and zero-padded, AM weighted + spindown-corr SRC-frame TS
  fftwf_plan fftplan;					// buffer FFT plan for given numSamplesOut length

  // arrays of size numFreqBinsOut over frequency bins f_k:
  UINT4 numFreqBinsOut;					// number of output frequency bins {f_k}
  COMPLEX8 *FaX_k;					// properly normalized F_a^X(f_k) over output bins
  COMPLEX8 *FbX_k;					// properly normalized F_b^X(f_k) over output bins
  COMPLEX8 *Fa_k;					// properly normalized F_a(f_k) over output bins
  COMPLEX8 *Fb_k;					// properly normalized F_b(f_k) over output bins
  UINT4 numFreqBinsAlloc;				// internal: keep track of allocated length of frequency-arrays

  ResampTimingInfo timingInfo;	 // temporary storage for collecting timing data
};

struct tagFstatInput_Resamp
{
  MultiCOMPLEX8TimeSeries  *multiTimeSeries_DET;	// input SFTs converted into a heterodyned timeseries
  // ----- buffering -----
  PulsarDopplerParams prev_doppler;			// buffering: previous phase-evolution ("doppler") parameters
  MultiAMCoeffs *prev_multiAMcoef;			// buffering: previous AM-coeffs, unique to skypos

  MultiCOMPLEX8TimeSeries *prev_multiTimeSeries_SRC;	// buffering: multi-detector SRC-frame timeseries
  MultiUINT4Vector *prev_multiSFTinds_SRC;		// buffering: SFT timestamps translated into SRC frame

  BOOLEAN ownThisWorkspace;				// flag whether we 'own' or share this workspace (ie who is responsible for freeing it)
  FstatWorkspace *ws;					// 'workspace': pre-allocated vectors used to store intermediate results
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
XLALComputeFaFb_Resamp ( FstatWorkspace *ws,
                         const PulsarDopplerParams thisPoint,
                         REAL8 dFreq,
                         const COMPLEX8TimeSeries *TimeSeries_SRC,
                         const UINT4Vector *SFTinds_SRC,
                         const AMCoeffs *ab
                         );

int XLALAppendResampInfo2File ( FILE *fp, const FstatInput *input );
static FstatWorkspace *XLALCreateFstatWorkspace ( UINT4 numSamplesSRC, UINT4 numSamplesFFT );

// ==================== function definitions ====================

// ---------- exported API functions ----------
///
/// Create a new workspace with given time samples in SRC frame 'numSamplesSRC' (holds time-series for spindown-correction)
/// and given total number of time-samples for FFTing (includes zero-padding for frequency-resolution)
///
static FstatWorkspace *
XLALCreateFstatWorkspace ( UINT4 numSamplesSRC,
                           UINT4 numSamplesFFT
                           )
{
  FstatWorkspace *ws;
  XLAL_CHECK_NULL ( (ws = XLALCalloc ( 1, sizeof(*ws))) != NULL, XLAL_ENOMEM );
  XLAL_CHECK_NULL ( (ws->TimeSeriesSpinCorr_SRC = XLALCreateCOMPLEX8Vector ( numSamplesSRC )) != NULL, XLAL_EFUNC );
  XLAL_CHECK_NULL ( (ws->FabX_Raw = fftw_malloc ( numSamplesFFT * sizeof(ws->FabX_Raw[0]) )) != NULL, XLAL_EFUNC );

  LAL_FFTW_WISDOM_LOCK;
  XLAL_CHECK_NULL ( (ws->fftplan = fftwf_plan_dft_1d ( numSamplesFFT, ws->FabX_Raw, ws->FabX_Raw, FFTW_FORWARD, FFTW_MEASURE )) != NULL, XLAL_EFAILED, "fftwf_plan_dft_1d() failed\n");
  LAL_FFTW_WISDOM_UNLOCK;
  ws->numSamplesFFT = numSamplesFFT;

  return ws;
} // XLALCreateFstatWorkspace()

///
/// Function to extract a workspace from a resampling setup, which can be passed in FstatOptionalArgs to be shared by various setups
/// in order to save memory. Note, when using this, you need to free this workspace yourself at the end using XLALDestroyFstatWorkspace().
/// Note: Demod methods don't use a workspace, so NULL (without error) is returned in this case.
///
FstatWorkspace *
XLALGetSharedFstatWorkspace ( FstatInput *input		//!< [in,out] Fstat input structure to extract shared workspace from
                              )
{
  XLAL_CHECK_NULL ( input != NULL, XLAL_EINVAL );

  if ( input->resamp == NULL ) {
    return NULL;
  }

  input->resamp->ownThisWorkspace = 0;	// the caller now owns the workspace and has to free it
  return input->resamp->ws;

} // XLALGetSharedFstatWorkspace()


void
XLALDestroyFstatWorkspace ( FstatWorkspace *ws )
{
  if ( ws == NULL ) {
    return;
  }

  XLALDestroyCOMPLEX8Vector ( ws->TimeSeriesSpinCorr_SRC );

  LAL_FFTW_WISDOM_LOCK;
  fftwf_destroy_plan ( ws->fftplan );
  LAL_FFTW_WISDOM_UNLOCK;

  fftw_free ( ws->FabX_Raw );

  XLALFree ( ws->FaX_k );
  XLALFree ( ws->FbX_k );
  XLALFree ( ws->Fa_k );
  XLALFree ( ws->Fb_k );

  XLALFree ( ws );
  return;

} // XLALDestroyFstatWorkspace()

/// debug/optimizer helper function: dump internal info from resampling code into a file
/// if called with input==NULL, output a header-comment-line
int
XLALAppendResampInfo2File ( FILE *fp, const FstatInput *input )
{
  XLAL_CHECK ( fp != NULL, XLAL_EINVAL );

  if ( input == NULL ) {
    fprintf (fp, "%%%%%8s %10s %6s %10s %10s ",
             "Nfreq", "NsFFT", "Nsft0", "Ns_DET0", "Ns_SRC0" );
    fprintf (fp, "%10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n",
             "tauTotal", "tauFFT", "tauBary", "tauSpin", "tauAM", "tauNorm", "tauFab2F", "tauMem", "tauSumFabX", "tauF1NoBuf", "tauF1Buf" );
    return XLAL_SUCCESS;
  }
  XLAL_CHECK ( input->resamp != NULL, XLAL_EINVAL );
  const FstatInput_Resamp *resamp = input->resamp;
  const FstatWorkspace *ws = resamp->ws;

  fprintf (fp, "%10d %10d", ws->numFreqBinsOut, ws->numSamplesFFT );
  UINT4 numSamples_DETX0 = resamp->multiTimeSeries_DET->data[0]->data->length;
  UINT4 numSFTs_X0 = (resamp->prev_multiAMcoef != NULL) ? resamp->prev_multiAMcoef->data[0]->a->length : 0;
  COMPLEX8TimeSeries *ts_SRCX0 = resamp->prev_multiTimeSeries_SRC->data[0];
  UINT4 numSamples_SRCX0 = ts_SRCX0->data->length;
  fprintf (fp, " %6d %10d %10d ", numSFTs_X0, numSamples_DETX0, numSamples_SRCX0 );

  const ResampTimingInfo *ti = &(ws->timingInfo);
  fprintf (fp, "%10.1e %10.1e %10.1e %10.1e %10.1e %10.1e %10.1e %10.1e %10.1e %10.1e %10.1e\n",
           ti->tauTotal, ti->tauFFT, ti->tauBary, ti->tauSpin, ti->tauAM, ti->tauNorm, ti->tauFab2F, ti->tauMem, ti->tauSumFabX, ti->tauF1NoBuf, ti->tauF1Buf );

  return XLAL_SUCCESS;

} // XLALAppendResampInfo2File()

// ---------- internal functions ----------
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
  if ( resamp->ownThisWorkspace ) {
    XLALDestroyFstatWorkspace ( resamp->ws );
  }

  XLALFree ( resamp );

  return;

} // DestroyFstatInput_Resamp()

static int
SetupFstatInput_Resamp ( FstatInput_Resamp *resamp,
                         const FstatInput_Common *common,
                         MultiSFTVector *multiSFTs,
                         FstatWorkspace *sharedWorkspace
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
  // round numSamplesFFT to next power of 2
  numSamplesFFT = (UINT4) pow ( 2, ceil(log2(numSamplesFFT)));
  REAL8 dt_SRC = TspanFFT / numSamplesFFT;			// adjust sampling rate to allow achieving exact requested dFreq=1/TspanFFT !

  // ----- prepare memory for the SRC-frame resampled timeseries buffer
  XLAL_CHECK ( (resamp->prev_multiTimeSeries_SRC = XLALCalloc ( 1, sizeof(MultiCOMPLEX8TimeSeries)) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( (resamp->prev_multiTimeSeries_SRC->data = XLALCalloc ( numDetectors, sizeof(COMPLEX8TimeSeries) )) != NULL, XLAL_ENOMEM );
  resamp->prev_multiTimeSeries_SRC->length = numDetectors;

  XLAL_CHECK ( (resamp->prev_multiSFTinds_SRC = XLALCalloc ( 1, sizeof(MultiUINT4Vector) )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( (resamp->prev_multiSFTinds_SRC->data = XLALCalloc ( numDetectors, sizeof(UINT4Vector) )) != NULL, XLAL_EFUNC );
  resamp->prev_multiSFTinds_SRC->length = numDetectors;

  UINT4 numSamplesSRCMax = 0;
  for ( UINT4 X = 0; X < numDetectors; X ++ )
    {
      REAL8 dt_DETX = resamp->multiTimeSeries_DET->data[X]->deltaT;
      REAL8 fHetX = resamp->multiTimeSeries_DET->data[X]->f0;
      XLAL_CHECK ( dt_DET == dt_DETX, XLAL_EINVAL, "Input timeseries must have identical 'deltaT' (%.3g != %.3g)\n", dt_DET, dt_DETX);
      XLAL_CHECK ( fHet == fHetX, XLAL_EINVAL, "Input timeseries must have identical heterodyning frequency 'f0'\n", fHet, fHetX );
      const char *nameX = resamp->multiTimeSeries_DET->data[X]->name;
      UINT4 numSamplesInX = resamp->multiTimeSeries_DET->data[X]->data->length;
      UINT4 numSamplesSRCX = (UINT4)ceil ( numSamplesInX * dt_DET / dt_SRC );
      numSamplesSRCMax = MYMAX ( numSamplesSRCMax, numSamplesSRCX );
      LIGOTimeGPS epoch0 = {0,0};	// will be set to corresponding SRC-frame epoch when barycentering
      XLAL_CHECK ( (resamp->prev_multiTimeSeries_SRC->data[X] = XLALCreateCOMPLEX8TimeSeries ( nameX, &epoch0, fHet, dt_SRC, &lalDimensionlessUnit, numSamplesSRCX )) != NULL, XLAL_EFUNC );

      UINT4 numTimestampsX = common->multiTimestamps->data[X]->length;
      XLAL_CHECK ( (resamp->prev_multiSFTinds_SRC->data[X] = XLALCreateUINT4Vector ( 2 * numTimestampsX )) != NULL, XLAL_EFUNC );
    } // for X < numDetectors

  // ---- re-use shared workspace, or allocate here ----------

  if ( sharedWorkspace != NULL )
    {
      XLAL_CHECK ( numSamplesFFT == sharedWorkspace->numSamplesFFT, XLAL_EINVAL, "Shared workspace of different frequency resolution: numSamplesFFT = %d != %d\n",
                   sharedWorkspace->numSamplesFFT, numSamplesFFT );

      // adjust maximal SRC-frame timeseries length, if necessary
      if ( numSamplesSRCMax > sharedWorkspace->TimeSeriesSpinCorr_SRC->length ) {
        XLAL_CHECK ( (sharedWorkspace->TimeSeriesSpinCorr_SRC->data = XLALRealloc ( sharedWorkspace->TimeSeriesSpinCorr_SRC->data, numSamplesSRCMax * sizeof(COMPLEX8) )) != NULL, XLAL_ENOMEM );
        sharedWorkspace->TimeSeriesSpinCorr_SRC->length = numSamplesSRCMax;
      }
      resamp->ws = sharedWorkspace;
      resamp->ownThisWorkspace = 0;
    } // end: if shared workspace given
  else
    {
      XLAL_CHECK ( ( resamp->ws = XLALCreateFstatWorkspace ( numSamplesSRCMax, numSamplesFFT )) != NULL, XLAL_EFUNC );
      resamp->ownThisWorkspace = 1;
    } // end: if we create our own workspace

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

#ifdef COLLECT_TIMING
  // collect internal timing info
  XLAL_INIT_MEM ( resamp->ws->timingInfo );
  ResampTimingInfo *ti = &(resamp->ws->timingInfo);
  REAL8 ticStart,tocEnd;
  ticStart = XLALGetCPUTime();
  REAL8 tic,toc;
#endif

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

  FstatWorkspace *ws = resamp->ws;

#ifdef COLLECT_TIMING
  tic = XLALGetCPUTime();
#endif
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
#ifdef COLLECT_TIMING
  toc = XLALGetCPUTime();
  ti->tauBary = (toc-tic);
#endif

  resamp->prev_doppler = thisPoint;

  MultiCOMPLEX8TimeSeries *multiTimeSeries_SRC = resamp->prev_multiTimeSeries_SRC;
  MultiUINT4Vector *multiSFTinds_SRC = resamp->prev_multiSFTinds_SRC;

  // ============================== check workspace is properly allocated and initialized ===========

  // ----- workspace that depends on number of output frequency bins 'numFreqBins' ----------
  UINT4 numFreqBins = Fstats->numFreqBins;

#ifdef COLLECT_TIMING
  tic = XLALGetCPUTime();
#endif

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
      if ( numFreqBins > ws->numFreqBinsAlloc )
        {
          XLAL_CHECK ( (ws->Fa_k = XLALRealloc ( ws->Fa_k, numFreqBins * sizeof(COMPLEX8))) != NULL, XLAL_ENOMEM );
          XLAL_CHECK ( (ws->Fb_k = XLALRealloc ( ws->Fb_k, numFreqBins * sizeof(COMPLEX8))) != NULL, XLAL_ENOMEM );
        } // only increase workspace arrays
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
      if ( numFreqBins > ws->numFreqBinsAlloc )
        {
          XLAL_CHECK ( (ws->FaX_k = XLALRealloc ( ws->FaX_k, numFreqBins * sizeof(COMPLEX8))) != NULL, XLAL_ENOMEM );
          XLAL_CHECK ( (ws->FbX_k = XLALRealloc ( ws->FbX_k, numFreqBins * sizeof(COMPLEX8))) != NULL, XLAL_ENOMEM );
        } // only increase workspace arrays
    }
  if ( numFreqBins > ws->numFreqBinsAlloc ) {
    ws->numFreqBinsAlloc = numFreqBins;	// keep track of allocated array length
  }
  ws->numFreqBinsOut = numFreqBins;
  // ====================================================================================================

  // store AM coefficient integrals in local variables
  REAL4 Ad = multiAMcoef->Mmunu.Ad;
  REAL4 Bd = multiAMcoef->Mmunu.Bd;
  REAL4 Cd = multiAMcoef->Mmunu.Cd;
  REAL4 Ed = multiAMcoef->Mmunu.Ed;
  REAL4 Dd = multiAMcoef->Mmunu.Dd;
  REAL4 Dd_inv = 1.0f / Dd;

#ifdef COLLECT_TIMING
  toc = XLALGetCPUTime();
  ti->tauMem = (toc-tic);	// this one doesn't scale with number of detector!
#endif

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
      XLAL_CHECK ( XLALComputeFaFb_Resamp ( resamp->ws, thisPoint, common->dFreq, TimeSeriesX_SRC, SFTindsX_SRC, abX ) == XLAL_SUCCESS, XLAL_EFUNC );

#ifdef COLLECT_TIMING
      tic = XLALGetCPUTime();
#endif
      if ( X == 0 )
        { // avoid having to memset this array: for the first detector we *copy* results
          for ( UINT4 k = 0; k < numFreqBins; k++ )
            {
              ws->Fa_k[k] = ws->FaX_k[k];
              ws->Fb_k[k] = ws->FbX_k[k];
            }
        } // end: if X==0
      else
        { // for subsequent detectors we *add to* them
          for ( UINT4 k = 0; k < numFreqBins; k++ )
            {
              ws->Fa_k[k] += ws->FaX_k[k];
              ws->Fb_k[k] += ws->FbX_k[k];
            }
        } // end:if X>0
#ifdef COLLECT_TIMING
      toc = XLALGetCPUTime();
      ti->tauSumFabX += (toc-tic);
      tic = toc;
#endif
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
            }  // for k < numFreqBins
        } // end: if compute F_X
#ifdef COLLECT_TIMING
      toc = XLALGetCPUTime();
      ti->tauFab2F += ( toc - tic );
#endif

    } // for X < numDetectors

#ifdef COLLECT_TIMING
  ti->tauSumFabX /= numDetectors;
  ti->tauFab2F /= numDetectors;
  tic = XLALGetCPUTime();
#endif
  if ( whatToCompute & FSTATQ_2F )
    {
      for ( UINT4 k=0; k < numFreqBins; k++ )
        {
          Fstats->twoF[k] = XLALComputeFstatFromFaFb ( ws->Fa_k[k], ws->Fb_k[k], Ad, Bd, Cd, Ed, Dd_inv );
        } // for k < numFreqBins
    } // if FSTATQ_2F
#ifdef COLLECT_TIMING
      toc = XLALGetCPUTime();
      ti->tauFab2F += ( toc - tic );
#endif

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

#ifdef COLLECT_TIMING
  // timings are per-detector
  tocEnd = XLALGetCPUTime();
  ti->tauTotal = (tocEnd - ticStart);
  // rescale all relevant timings to single-IFO case
  ti->tauTotal /= numDetectors;
  ti->tauBary  /= numDetectors;
  ti->tauSpin  /= numDetectors;
  ti->tauAM    /= numDetectors;
  ti->tauFFT   /= numDetectors;
  ti->tauNorm  /= numDetectors;

  // compute 'fundamental' timing numbers per template per detector
  ti->tauF1NoBuf = ti->tauTotal / numFreqBins;
  ti->tauF1Buf   = (ti->tauTotal - ti->tauBary - ti->tauMem) / numFreqBins;
#endif

  return XLAL_SUCCESS;

} // ComputeFstat_Resamp()


static int
XLALComputeFaFb_Resamp ( FstatWorkspace *ws,				//!< [in,out] pre-allocated 'workspace' for temporary and output quantities
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
  REAL8 fMinFFT = fHet + freqShift - dFreq * (ws->numSamplesFFT/2);	// we'll shift DC into the *middle bin* N/2  [N always even!]
  UINT4 offset_bins = (UINT4) lround ( ( FreqOut0 - fMinFFT ) / dFreq );

#ifdef COLLECT_TIMING
  // collect some internal timing info
  ResampTimingInfo *ti = &(ws->timingInfo);
  REAL8 tic,toc;
  tic = XLALGetCPUTime();
#endif

  // apply spindown phase-factors, store result in 'workspace'
  memset ( ws->TimeSeriesSpinCorr_SRC->data, 0, ws->TimeSeriesSpinCorr_SRC->length * sizeof(COMPLEX8));
  XLAL_CHECK ( XLALApplySpindownAndFreqShift ( ws->TimeSeriesSpinCorr_SRC, TimeSeries_SRC, SFTinds_SRC, &thisPoint, freqShift ) == XLAL_SUCCESS, XLAL_EFUNC );

#ifdef COLLECT_TIMING
  toc = XLALGetCPUTime();
  ti->tauSpin += ( toc - tic);
  tic = toc;
#endif
  // ----- compute FaX_k
  // apply amplitude modulation factors {a,b}, store result in zero-padded timeseries for FFTing
  memset ( ws->FabX_Raw, 0, ws->numSamplesFFT * sizeof(ws->FabX_Raw[0]) );
  XLAL_CHECK ( XLALApplyAmplitudeModulation ( ws->FabX_Raw, ws->TimeSeriesSpinCorr_SRC, SFTinds_SRC, ab->a ) == XLAL_SUCCESS, XLAL_EFUNC );

#ifdef COLLECT_TIMING
  toc = XLALGetCPUTime();
  ti->tauAM += ( toc - tic);
  tic = toc;
#endif

  // Fourier transform the resampled Fa(t)
  fftwf_execute ( ws->fftplan );

  for ( UINT4 k = 0; k < ws->numFreqBinsOut; k++ ) {
    ws->FaX_k[k] = ws->FabX_Raw [ offset_bins + k ];
  }

#ifdef COLLECT_TIMING
  toc = XLALGetCPUTime();
  ti->tauFFT += ( toc - tic);
  tic = toc;
#endif

  // ----- compute FbX_k
  // apply amplitude modulation factors {a,b}, store result in zero-padded timeseries for FFTing
  memset ( ws->FabX_Raw, 0, ws->numSamplesFFT * sizeof(ws->FabX_Raw[0]) );
  XLAL_CHECK ( XLALApplyAmplitudeModulation ( ws->FabX_Raw, ws->TimeSeriesSpinCorr_SRC, SFTinds_SRC, ab->b ) == XLAL_SUCCESS, XLAL_EFUNC );

#ifdef COLLECT_TIMING
  toc = XLALGetCPUTime();
  ti->tauAM += ( toc - tic);
  tic = toc;
#endif

  // Fourier transform the resampled Fa(t)
  fftwf_execute ( ws->fftplan );

  for ( UINT4 k = 0; k < ws->numFreqBinsOut; k++ ) {
    ws->FbX_k[k] = ws->FabX_Raw [ offset_bins + k ];
  }

#ifdef COLLECT_TIMING
  toc = XLALGetCPUTime();
  ti->tauFFT += ( toc - tic);
  tic = toc;
#endif

  // ----- normalization factors to be applied to Fa and Fb:
  const REAL8 dtauX = XLALGPSDiff ( &TimeSeries_SRC->epoch, &thisPoint.refTime );
  for ( UINT4 k = 0; k < ws->numFreqBinsOut; k++ )
    {
      REAL8 f_k = FreqOut0 + k * dFreq;
      REAL8 cycles = - f_k * dtauX;
      REAL4 sinphase, cosphase;
      XLALSinCos2PiLUT ( &sinphase, &cosphase, cycles );
      COMPLEX8 normX_k = dt_SRC * crectf ( cosphase, sinphase );
      ws->FaX_k[k] *= normX_k;
      ws->FbX_k[k] *= normX_k;
    } // for k < numFreqBinsOut

#ifdef COLLECT_TIMING
  toc = XLALGetCPUTime();
  ti->tauNorm += ( toc - tic);
  tic = toc;
#endif

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
      COMPLEX8 fact = (COMPLEX8) ab->data[alpha];
      UINT4 start_index = SFTinds->data[2*alpha];
      UINT4 end_index   = SFTinds->data[2*alpha+1];

      // NOTE: this relies on the number of FFT bins being EVEN!
      // multiply input TS by (-1)^j to bring DC into the middle bin (k=N/2)
      fact = ((start_index % 2) == 0) ? fact : -fact;
      // loop over all samples from this SFT
      // and apply amplitude modulation factor to output timeseries
      for ( UINT4 j=start_index; j <= end_index; j ++ )
        {
          xOut[j] = fact * xIn->data[j];
          fact = -fact;
        } // for j in [start_index, end_index]

    } // for alpha < numSFTs

  return XLAL_SUCCESS;

} // XLALApplyAmplitudeModulation()


///
/// Performs barycentric resampling on a multi-detector timeseries
///
static int
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
static int
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
