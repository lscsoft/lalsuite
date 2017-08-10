//
// Copyright (C) 2009, 2014--2015 Reinhard Prix
// Copyright (C) 2012--2015 Karl Wette
// Copyright (C) 2009 Chris Messenger, Pinkesh Patel, Xavier Siemens, Holger Pletsch
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#include "ComputeFstat_internal.h"

#include <lal/FFTWMutex.h>
#include <lal/Factorial.h>
#include <lal/LFTandTSutils.h>
#include <lal/LogPrintf.h>
#include <lal/SinCosLUT.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>

// ========== Resamp internals ==========

// ----- local macros ----------
#define MYMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MYMIN(x,y) ( (x) < (y) ? (x) : (y) )

// local macro versions of library functions to avoid calling external functions in GPU-ready code
#define GPSDIFF(x,y) (1.0*((x).gpsSeconds - (y).gpsSeconds) + ((x).gpsNanoSeconds - (y).gpsNanoSeconds)*1e-9)
#define GPSGETREAL8(x) ( (x)->gpsSeconds + ( (x)->gpsNanoSeconds / XLAL_BILLION_REAL8 ) );
#define GPSSETREAL8(gps,r8) do {                                        \
    (gps).gpsSeconds     = (UINT4)floor(r8);                            \
    (gps).gpsNanoSeconds = (UINT4)round ( ((r8) - (gps).gpsSeconds) * XLAL_BILLION_REAL8 ); \
    if ( (gps).gpsNanoSeconds == XLAL_BILLION_INT4 ) {                  \
      (gps).gpsSeconds += 1;                                            \
      (gps).gpsNanoSeconds = 0;                                         \
    }                                                                   \
  } while(0)


// ----- local constants

// ----- local types ----------

// ---------- BEGIN: Resamp-specific timing model data ----------
typedef struct tagTimings_t
{
  REAL4 Total;		// total time spent in XLALComputeFstatResamp()
  REAL4 Bary;		// time spent (in this call) in barycentric resampling
  REAL4 Spin;		// time spent in spindown+frequency correction
  REAL4 FFT;		// time spent in FFT
  REAL4 Copy;		// time spent copying results from FFT to FabX
  REAL4 Norm;		// time spent normalizing the final Fa,Fb
  REAL4 Fab2F;		// time to compute Fstat from {Fa,Fb}
  REAL4 Mem;		// time to realloc and Memset-0 arrays
  REAL4 SumFabX;	// time to sum_X Fab^X
  BOOLEAN BufferRecomputed; // did we need to recompute the buffer this time?
} Timings_t;

// Resamp-specific timing model data
typedef struct tagFstatTimingResamp
{
  UINT4 NsampFFT0;		// original number of FFT samples (not rounded to power-of-two)
  UINT4 NsampFFT;		// actual number of FFT samples (rounded up to power-of-two if optArgs->resampFFTPowerOf2 == true)
  REAL4 Resolution;	// (internal) frequency resolution 'R' in natural units: df_internal = R / T_FFT\n

  REAL4 tau0_Fbin;      // timing coefficient for all contributions scaling with output frequency-bins
  REAL4 tau0_spin;      // timing coefficient for spindown-correction
  REAL4 tau0_FFT;       // timing coefficient for FFT-time
  REAL4 tau0_bary;      // timing coefficient for barycentering

  Timings_t Tau;

} FstatTimingResamp;

static char FstatTimingResampHelp[] =
  "%%%% ----- Resampling-specific timing model -----\n"
  "%%%% NsampFFT0:      original number of FFT samples (not yet rounded up to power-of-two)\n"
  "%%%% NsampFFT:       actual number of FFT samples (rounded to power-of-two if optArgs->resampFFTPowerOf2 == true)\n"
  "%%%% R:              (internal) frequency resolution in natural units: df_internal = R / T_FFT\n"
  "%%%%\n"
  "%%%% tau0_Fbin:      timing coefficient for all contributions scaling with output frequency-bins\n"
  "%%%% tau0_spin:      timing coefficient for spindown-correction\n"
  "%%%% tau0_FFT:       timing coefficient for FFT-time\n"
  "%%%% tau0_bary:      timing coefficient for barycentering\n"
  "%%%%\n"
  "%%%% Resampling F-statistic timing model:\n"
  "%%%% tauF_core       = tau0_Fbin + (NsampFFT/NFbin) * ( R * tau0_spin + 5 * log2(NsampFFT) * tau0_FFT )\n"
  "%%%% tauF_buffer     = R * NsampFFT * tau0_bary / NFbin\n"
  "%%%%"
  "";
// ---------- END: Resamp-specific timing model data ----------


// ----- workspace ----------
typedef struct tagResampWorkspace
{
  // intermediate quantities to interpolate and operate on SRC-frame timeseries
  COMPLEX8Vector *TStmp1_SRC;	// can hold a single-detector SRC-frame spindown-corrected timeseries [without zero-padding]
  COMPLEX8Vector *TStmp2_SRC;	// can hold a single-detector SRC-frame spindown-corrected timeseries [without zero-padding]
  REAL8Vector *SRCtimes_DET;	// holds uniformly-spaced SRC-frame timesteps translated into detector frame [for interpolation]

  // input padded timeseries ts(t) and output Fab(f) of length 'numSamplesFFT' and corresponding fftw plan
  UINT4 numSamplesFFTAlloc;	// allocated number of zero-padded SRC-frame time samples (related to dFreq)
  COMPLEX8 *TS_FFT;		// zero-padded, spindown-corr SRC-frame TS
  COMPLEX8 *FabX_Raw;		// raw full-band FFT result Fa,Fb

  // arrays of size numFreqBinsOut over frequency bins f_k:
  COMPLEX8 *FaX_k;		// properly normalized F_a^X(f_k) over output bins
  COMPLEX8 *FbX_k;		// properly normalized F_b^X(f_k) over output bins
  COMPLEX8 *Fa_k;		// properly normalized F_a(f_k) over output bins
  COMPLEX8 *Fb_k;		// properly normalized F_b(f_k) over output bins
  UINT4 numFreqBinsAlloc;	// internal: keep track of allocated length of frequency-arrays

} ResampWorkspace;

typedef struct
{
  UINT4 Dterms;						// Number of terms to use (on either side) in Windowed-Sinc interpolation kernel
  MultiCOMPLEX8TimeSeries  *multiTimeSeries_DET;	// input SFTs converted into a heterodyned timeseries
  // ----- buffering -----
  PulsarDopplerParams prev_doppler;			// buffering: previous phase-evolution ("doppler") parameters
  MultiAMCoeffs *multiAMcoef;				// buffered antenna-pattern functions
  MultiSSBtimes *multiSSBtimes;				// buffered SSB times, including *only* sky-position corrections, not binary
  MultiSSBtimes *multiBinaryTimes;			// buffered SRC times, including both sky- and binary corrections [to avoid re-allocating this]

  AntennaPatternMatrix Mmunu;				// combined multi-IFO antenna-pattern coefficients {A,B,C,E}
  AntennaPatternMatrix MmunuX[PULSAR_MAX_DETECTORS];	// per-IFO antenna-pattern coefficients {AX,BX,CX,EX}

  MultiCOMPLEX8TimeSeries *multiTimeSeries_SRC_a;	// multi-detector SRC-frame timeseries, multiplied by AM function a(t)
  MultiCOMPLEX8TimeSeries *multiTimeSeries_SRC_b;	// multi-detector SRC-frame timeseries, multiplied by AM function b(t)

  UINT4 numSamplesFFT;					// length of zero-padded SRC-frame timeseries (related to dFreq)
  UINT4 decimateFFT;					// output every n-th frequency bin, with n>1 iff (dFreq > 1/Tspan), and was internally decreased by n
  fftwf_plan fftplan;					// FFT plan

  // ----- timing -----
  BOOLEAN collectTiming;				// flag whether or not to collect timing information
  FstatTimingGeneric timingGeneric;			// measured (generic) F-statistic timing values
  FstatTimingResamp  timingResamp;			// measured Resamp-specific timing model data

} ResampMethodData;


// ----- local prototypes ----------

int XLALSetupFstatResamp ( void **method_data, FstatCommon *common, FstatMethodFuncs* funcs, MultiSFTVector *multiSFTs, const FstatOptionalArgs *optArgs );

static int
XLALComputeFstatResamp ( FstatResults* Fstats,
                         const FstatCommon *common,
                         void *method_data
                       );

static int
XLALApplySpindownAndFreqShift ( COMPLEX8 *xOut,
                                const COMPLEX8TimeSeries *xIn,
                                const PulsarDopplerParams *doppler,
                                REAL8 freqShift
                                );

static int
XLALBarycentricResampleMultiCOMPLEX8TimeSeries ( ResampMethodData *resamp,
                                                 const PulsarDopplerParams *thisPoint,
                                                 const FstatCommon *common
                                                 );

static int
XLALComputeFaFb_Resamp ( ResampMethodData *resamp,
                         ResampWorkspace *ws,
                         const PulsarDopplerParams thisPoint,
                         REAL8 dFreq,
                         UINT4 numFreqBins,
                         const COMPLEX8TimeSeries *TimeSeries_SRC_a,
                         const COMPLEX8TimeSeries *TimeSeries_SRC_b
                         );

static void
XLALGetFFTPlanHints ( int * planMode,
                      double * planGenTimeoutSeconds
                      );

// ==================== function definitions ====================

static void
XLALDestroyResampWorkspace ( void *workspace )
{
  ResampWorkspace *ws = (ResampWorkspace*) workspace;

  XLALDestroyCOMPLEX8Vector ( ws->TStmp1_SRC );
  XLALDestroyCOMPLEX8Vector ( ws->TStmp2_SRC );
  XLALDestroyREAL8Vector ( ws->SRCtimes_DET );

  fftw_free ( ws->FabX_Raw );
  fftw_free ( ws->TS_FFT );

  XLALFree ( ws->FaX_k );
  XLALFree ( ws->FbX_k );
  XLALFree ( ws->Fa_k );
  XLALFree ( ws->Fb_k );

  XLALFree ( ws );
  return;

} // XLALDestroyResampWorkspace()

// ---------- internal functions ----------
static void
XLALDestroyResampMethodData ( void* method_data )
{

  ResampMethodData *resamp = (ResampMethodData*) method_data;

  XLALDestroyMultiCOMPLEX8TimeSeries (resamp->multiTimeSeries_DET );

  // ----- free buffer
  XLALDestroyMultiCOMPLEX8TimeSeries ( resamp->multiTimeSeries_SRC_a );
  XLALDestroyMultiCOMPLEX8TimeSeries ( resamp->multiTimeSeries_SRC_b );
  XLALDestroyMultiAMCoeffs ( resamp->multiAMcoef );
  XLALDestroyMultiSSBtimes ( resamp->multiSSBtimes );
  XLALDestroyMultiSSBtimes ( resamp->multiBinaryTimes );

  LAL_FFTW_WISDOM_LOCK;
  fftwf_destroy_plan ( resamp->fftplan );
  LAL_FFTW_WISDOM_UNLOCK;

  XLALFree ( resamp );

} // XLALDestroyResampMethodData()

int
XLALSetupFstatResamp ( void **method_data,
                       FstatCommon *common,
                       FstatMethodFuncs* funcs,
                       MultiSFTVector *multiSFTs,
                       const FstatOptionalArgs *optArgs
                     )
{
  // Check input
  XLAL_CHECK ( method_data != NULL, XLAL_EFAULT );
  XLAL_CHECK ( common != NULL, XLAL_EFAULT );
  XLAL_CHECK ( funcs != NULL, XLAL_EFAULT );
  XLAL_CHECK ( multiSFTs != NULL, XLAL_EFAULT );
  XLAL_CHECK ( optArgs != NULL, XLAL_EFAULT );

  // Allocate method data
  ResampMethodData *resamp = *method_data = XLALCalloc( 1, sizeof(*resamp) );
  XLAL_CHECK( resamp != NULL, XLAL_ENOMEM );

  resamp->Dterms = optArgs->Dterms;

  // Set method function pointers
  funcs->compute_func = XLALComputeFstatResamp;
  funcs->method_data_destroy_func = XLALDestroyResampMethodData;
  funcs->workspace_destroy_func = XLALDestroyResampWorkspace;

  // Extra band needed for resampling: Hamming-windowed sinc used for interpolation has a transition bandwith of
  // TB=(4/L)*fSamp, where L=2*Dterms+1 is the window-length, and here fSamp=Band (i.e. the full SFT frequency band)
  // However, we're only interested in the physical band and we'll be throwing away all bins outside of this.
  // This implies that we're only affected by *half* the transition band TB/2 on either side, as the other half of TB is outside of the band of interest
  // (and will actually get aliased, i.e. the region [-fNy - TB/2, -fNy] overlaps with [fNy-TB/2,fNy] and vice-versa: [fNy,fNy+TB/2] overlaps with [-fNy,-fNy+TB/2])
  // ==> therefore we only need to add an extra TB/2 on each side to be able to safely avoid the transition-band effects
  REAL8 f0 = multiSFTs->data[0]->data[0].f0;
  REAL8 dFreq = multiSFTs->data[0]->data[0].deltaF;
  REAL8 Band = multiSFTs->data[0]->data[0].data->length * dFreq;
  REAL8 extraBand = 2.0  / ( 2 * optArgs->Dterms + 1 ) * Band;
  XLAL_CHECK ( XLALMultiSFTVectorResizeBand ( multiSFTs, f0 - extraBand, Band + 2 * extraBand ) == XLAL_SUCCESS, XLAL_EFUNC );
  // Convert SFTs into heterodyned complex timeseries [in detector frame]
  XLAL_CHECK ( (resamp->multiTimeSeries_DET = XLALMultiSFTVectorToCOMPLEX8TimeSeries ( multiSFTs )) != NULL, XLAL_EFUNC );

  XLALDestroyMultiSFTVector ( multiSFTs );	// don't need them SFTs any more ...

  UINT4 numDetectors = resamp->multiTimeSeries_DET->length;
  REAL8 dt_DET       = resamp->multiTimeSeries_DET->data[0]->deltaT;
  REAL8 fHet         = resamp->multiTimeSeries_DET->data[0]->f0;
  REAL8 Tsft         = common->multiTimestamps->data[0]->deltaT;

  // determine resampled timeseries parameters
  REAL8 TspanFFT = 1.0 / common->dFreq;

  // check that TspanFFT >= TspanX for all detectors X, otherwise increase TspanFFT by an integer factor 'decimateFFT' such that this is true
  REAL8 TspanXMax = 0;
  for ( UINT4 X = 0; X < numDetectors; X ++ )
    {
      UINT4 numSamples_DETX = resamp->multiTimeSeries_DET->data[X]->data->length;
      REAL8 TspanX = numSamples_DETX * dt_DET;
      TspanXMax = fmax ( TspanXMax, TspanX );
    }
  UINT4 decimateFFT = (UINT4)ceil ( TspanXMax / TspanFFT );	// larger than 1 means we need to artificially increase dFreqFFT by 'decimateFFT'
  if ( decimateFFT > 1 ) {
    XLALPrintWarning ("WARNING: Frequency spacing larger than 1/Tspan, we'll internally decimate FFT frequency bins by a factor of %" LAL_UINT4_FORMAT "\n", decimateFFT );
  }
  TspanFFT *= decimateFFT;
  resamp->decimateFFT = decimateFFT;

  UINT4 numSamplesFFT0 = (UINT4) ceil ( TspanFFT / dt_DET );      // we use ceil() so that we artificially widen the band rather than reduce it
  UINT4 numSamplesFFT = 0;
  if ( optArgs->resampFFTPowerOf2 ) {
    numSamplesFFT = (UINT4) pow ( 2, ceil ( log2 ( numSamplesFFT0 ) ) );  // round numSamplesFFT up to next power of 2 for most effiecient FFT
  } else {
    numSamplesFFT = (UINT4) 2 * ceil ( numSamplesFFT0 / 2 );	// always ensure numSamplesFFT is even
  }

  REAL8 dt_SRC = TspanFFT / numSamplesFFT;			// adjust sampling rate to allow achieving exact requested dFreq=1/TspanFFT !

  resamp->numSamplesFFT = numSamplesFFT;
  // ----- allocate buffer Memory ----------

  // header for SRC-frame resampled timeseries buffer
  XLAL_CHECK ( (resamp->multiTimeSeries_SRC_a = XLALCalloc ( 1, sizeof(MultiCOMPLEX8TimeSeries)) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( (resamp->multiTimeSeries_SRC_a->data = XLALCalloc ( numDetectors, sizeof(COMPLEX8TimeSeries) )) != NULL, XLAL_ENOMEM );
  resamp->multiTimeSeries_SRC_a->length = numDetectors;

  XLAL_CHECK ( (resamp->multiTimeSeries_SRC_b = XLALCalloc ( 1, sizeof(MultiCOMPLEX8TimeSeries)) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( (resamp->multiTimeSeries_SRC_b->data = XLALCalloc ( numDetectors, sizeof(COMPLEX8TimeSeries) )) != NULL, XLAL_ENOMEM );
  resamp->multiTimeSeries_SRC_b->length = numDetectors;

  LIGOTimeGPS XLAL_INIT_DECL(epoch0);	// will be set to corresponding SRC-frame epoch when barycentering
  UINT4 numSamplesMax_SRC = 0;
  for ( UINT4 X = 0; X < numDetectors; X ++ )
    {
      // ----- check input consistency ----------
      REAL8 dt_DETX = resamp->multiTimeSeries_DET->data[X]->deltaT;
      XLAL_CHECK ( dt_DET == dt_DETX, XLAL_EINVAL, "Input timeseries must have identical 'deltaT(X=%d)' (%.16g != %.16g)\n", X, dt_DET, dt_DETX);

      REAL8 fHetX = resamp->multiTimeSeries_DET->data[X]->f0;
      XLAL_CHECK ( fHet == fHetX, XLAL_EINVAL, "Input timeseries must have identical heterodyning frequency 'f0(X=%d)' (%.16g != %.16g)\n", X, fHet, fHetX );

      REAL8 TsftX = common->multiTimestamps->data[X]->deltaT;
      XLAL_CHECK ( Tsft == TsftX, XLAL_EINVAL, "Input timestamps must have identical stepsize 'Tsft(X=%d)' (%.16g != %.16g)\n", X, Tsft, TsftX );

      // ----- prepare Memory fo SRC-frame timeseries and AM coefficients
      const char *nameX = resamp->multiTimeSeries_DET->data[X]->name;
      UINT4 numSamples_DETX = resamp->multiTimeSeries_DET->data[X]->data->length;
      UINT4 numSamples_SRCX = (UINT4)ceil ( numSamples_DETX * dt_DET / dt_SRC );

      XLAL_CHECK ( (resamp->multiTimeSeries_SRC_a->data[X] = XLALCreateCOMPLEX8TimeSeries ( nameX, &epoch0, fHet, dt_SRC, &lalDimensionlessUnit, numSamples_SRCX )) != NULL, XLAL_EFUNC );
      XLAL_CHECK ( (resamp->multiTimeSeries_SRC_b->data[X] = XLALCreateCOMPLEX8TimeSeries ( nameX, &epoch0, fHet, dt_SRC, &lalDimensionlessUnit, numSamples_SRCX )) != NULL, XLAL_EFUNC );

      numSamplesMax_SRC = MYMAX ( numSamplesMax_SRC, numSamples_SRCX );
    } // for X < numDetectors

  XLAL_CHECK ( numSamplesFFT >= numSamplesMax_SRC, XLAL_EFAILED, "[numSamplesFFT = %d] < [numSamplesMax_SRC = %d]\n", numSamplesFFT, numSamplesMax_SRC );

  // ---- re-use shared workspace, or allocate here ----------
  ResampWorkspace *ws = (ResampWorkspace*) common->workspace;
  if ( ws != NULL )
    {
      if ( numSamplesFFT > ws->numSamplesFFTAlloc )
        {
          fftw_free ( ws->FabX_Raw );
          XLAL_CHECK ( (ws->FabX_Raw = fftw_malloc ( numSamplesFFT * sizeof(COMPLEX8) )) != NULL, XLAL_ENOMEM );
          fftw_free ( ws->TS_FFT );
          XLAL_CHECK ( (ws->TS_FFT   = fftw_malloc ( numSamplesFFT * sizeof(COMPLEX8) )) != NULL, XLAL_ENOMEM );

          ws->numSamplesFFTAlloc = numSamplesFFT;
        }

      // adjust maximal SRC-frame timeseries length, if necessary
      if ( numSamplesMax_SRC > ws->TStmp1_SRC->length ) {
        XLAL_CHECK ( (ws->TStmp1_SRC->data = XLALRealloc ( ws->TStmp1_SRC->data,   numSamplesMax_SRC * sizeof(COMPLEX8) )) != NULL, XLAL_ENOMEM );
        ws->TStmp1_SRC->length = numSamplesMax_SRC;
        XLAL_CHECK ( (ws->TStmp2_SRC->data = XLALRealloc ( ws->TStmp2_SRC->data,   numSamplesMax_SRC * sizeof(COMPLEX8) )) != NULL, XLAL_ENOMEM );
        ws->TStmp2_SRC->length = numSamplesMax_SRC;
        XLAL_CHECK ( (ws->SRCtimes_DET->data = XLALRealloc ( ws->SRCtimes_DET->data, numSamplesMax_SRC * sizeof(REAL8) )) != NULL, XLAL_ENOMEM );
        ws->SRCtimes_DET->length = numSamplesMax_SRC;
      }

    } // end: if shared workspace given
  else
    {
      XLAL_CHECK ( (ws = XLALCalloc ( 1, sizeof(*ws))) != NULL, XLAL_ENOMEM );
      XLAL_CHECK ( (ws->TStmp1_SRC   = XLALCreateCOMPLEX8Vector ( numSamplesMax_SRC )) != NULL, XLAL_EFUNC );
      XLAL_CHECK ( (ws->TStmp2_SRC   = XLALCreateCOMPLEX8Vector ( numSamplesMax_SRC )) != NULL, XLAL_EFUNC );
      XLAL_CHECK ( (ws->SRCtimes_DET = XLALCreateREAL8Vector ( numSamplesMax_SRC )) != NULL, XLAL_EFUNC );

      XLAL_CHECK ( (ws->FabX_Raw = fftw_malloc ( numSamplesFFT * sizeof(COMPLEX8) )) != NULL, XLAL_ENOMEM );
      XLAL_CHECK ( (ws->TS_FFT   = fftw_malloc ( numSamplesFFT * sizeof(COMPLEX8) )) != NULL, XLAL_ENOMEM );
      ws->numSamplesFFTAlloc = numSamplesFFT;

      common->workspace = ws;
    } // end: if we create our own workspace

  // ----- compute and buffer FFT plan ----------
  int fft_plan_flags=FFTW_MEASURE;
  double fft_plan_timeout= FFTW_NO_TIMELIMIT ;
  char *wisdom_filename;
  static int tried_wisdom = 0;

  LAL_FFTW_WISDOM_LOCK;
  // if FFTWF_WISDOM_FILENAME is set, try to import that wisdom
  wisdom_filename = getenv("FFTWF_WISDOM_FILENAME");
  if (wisdom_filename && !tried_wisdom) {
    if (fftwf_import_wisdom_from_filename(wisdom_filename)) {
      XLALPrintInfo("INFO: imported wisdom from file '%s'\n", wisdom_filename);
    } else {
      XLALPrintWarning("WARNING: Couldn't import wisdom from file '%s'\n", wisdom_filename);
    }
    tried_wisdom = -1;
  }
  XLALGetFFTPlanHints (& fft_plan_flags , & fft_plan_timeout);
  fftw_set_timelimit( fft_plan_timeout );
  XLAL_CHECK ( (resamp->fftplan = fftwf_plan_dft_1d ( resamp->numSamplesFFT, ws->TS_FFT, ws->FabX_Raw, FFTW_FORWARD, fft_plan_flags )) != NULL, XLAL_EFAILED, "fftwf_plan_dft_1d() failed\n");
  LAL_FFTW_WISDOM_UNLOCK;

  // turn on timing collection if requested
  resamp->collectTiming = optArgs->collectTiming;

  // initialize struct for collecting timing data, store invariant 'meta' quantities about this setup
  if ( resamp->collectTiming )
    {
      XLAL_INIT_MEM ( resamp->timingGeneric );
      resamp->timingGeneric.Ndet = numDetectors;

      XLAL_INIT_MEM ( resamp->timingResamp );
      resamp->timingResamp.Resolution = TspanXMax / TspanFFT;
      resamp->timingResamp.NsampFFT0  = numSamplesFFT0;
      resamp->timingResamp.NsampFFT   = numSamplesFFT;
    }

  return XLAL_SUCCESS;

} // XLALSetupFstatResamp()


static int
XLALComputeFstatResamp ( FstatResults* Fstats,
                         const FstatCommon *common,
                         void *method_data
                       )
{
  // Check input
  XLAL_CHECK(Fstats != NULL, XLAL_EFAULT);
  XLAL_CHECK(common != NULL, XLAL_EFAULT);
  XLAL_CHECK(method_data != NULL, XLAL_EFAULT);

  ResampMethodData *resamp = (ResampMethodData*) method_data;

  const FstatQuantities whatToCompute = Fstats->whatWasComputed;
  XLAL_CHECK ( !(whatToCompute & FSTATQ_ATOMS_PER_DET), XLAL_EINVAL, "Resampling does not currently support atoms per detector" );

  ResampWorkspace *ws = (ResampWorkspace*) common->workspace;

  // ----- handy shortcuts ----------
  PulsarDopplerParams thisPoint = Fstats->doppler;
  const MultiCOMPLEX8TimeSeries *multiTimeSeries_DET = resamp->multiTimeSeries_DET;
  UINT4 numDetectors = multiTimeSeries_DET->length;

  // collect internal timing info
  BOOLEAN collectTiming = resamp->collectTiming;
  Timings_t *Tau = &(resamp->timingResamp.Tau);
  XLAL_INIT_MEM ( (*Tau) );	// these need to be initialized to 0 for each call

  REAL8 ticStart = 0, tocEnd = 0;
  REAL8 tic = 0, toc = 0;
  if ( collectTiming ) {
    XLAL_INIT_MEM ( (*Tau) );	// re-set all timings to 0 at beginning of each Fstat-call
    ticStart = XLALGetCPUTime();
  }
  // Note: all buffering is done within that function
  XLAL_CHECK ( XLALBarycentricResampleMultiCOMPLEX8TimeSeries ( resamp, &thisPoint, common ) == XLAL_SUCCESS, XLAL_EFUNC );

  if ( whatToCompute == FSTATQ_NONE ) {
    return XLAL_SUCCESS;
  }

  MultiCOMPLEX8TimeSeries *multiTimeSeries_SRC_a = resamp->multiTimeSeries_SRC_a;
  MultiCOMPLEX8TimeSeries *multiTimeSeries_SRC_b = resamp->multiTimeSeries_SRC_b;

  // ============================== check workspace is properly allocated and initialized ===========

  // ----- workspace that depends on maximal number of output frequency bins 'numFreqBins' ----------
  UINT4 numFreqBins = Fstats->numFreqBins;

  if ( collectTiming ) {
    tic = XLALGetCPUTime();
  }

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

  if ( collectTiming ) {
    toc = XLALGetCPUTime();
    Tau->Mem = (toc-tic);	// this one doesn't scale with number of detector!
  }
  // ====================================================================================================

  // loop over detectors
  for ( UINT4 X=0; X < numDetectors; X++ )
    {
      // if return-struct contains memory for holding FaFbPerDet: use that directly instead of local memory
      if ( whatToCompute & FSTATQ_FAFB_PER_DET )
        {
          ws->FaX_k = Fstats->FaPerDet[X];
          ws->FbX_k = Fstats->FbPerDet[X];
        }
      const COMPLEX8TimeSeries *TimeSeriesX_SRC_a = multiTimeSeries_SRC_a->data[X];
      const COMPLEX8TimeSeries *TimeSeriesX_SRC_b = multiTimeSeries_SRC_b->data[X];

      // compute {Fa^X(f_k), Fb^X(f_k)}: results returned via workspace ws
      XLAL_CHECK ( XLALComputeFaFb_Resamp ( resamp, ws, thisPoint, common->dFreq, numFreqBins, TimeSeriesX_SRC_a, TimeSeriesX_SRC_b ) == XLAL_SUCCESS, XLAL_EFUNC );

      if ( collectTiming ) {
        tic = XLALGetCPUTime();
      }
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

      if ( collectTiming ) {
        toc = XLALGetCPUTime();
        Tau->SumFabX += (toc-tic);
        tic = toc;
      }

      // ----- if requested: compute per-detector Fstat_X_k
      if ( whatToCompute & FSTATQ_2F_PER_DET )
        {
          const REAL4 AdX = resamp->MmunuX[X].Ad;
          const REAL4 BdX = resamp->MmunuX[X].Bd;
          const REAL4 CdX = resamp->MmunuX[X].Cd;
          const REAL4 EdX = resamp->MmunuX[X].Ed;
          const REAL4 DdX_inv = 1.0f / resamp->MmunuX[X].Dd;
          for ( UINT4 k = 0; k < numFreqBins; k ++ )
            {
              Fstats->twoFPerDet[X][k] = XLALComputeFstatFromFaFb ( ws->FaX_k[k], ws->FbX_k[k], AdX, BdX, CdX, EdX, DdX_inv );
            }  // for k < numFreqBins
        } // end: if compute F_X

      if ( collectTiming ) {
        toc = XLALGetCPUTime();
        Tau->Fab2F += ( toc - tic );
      }

    } // for X < numDetectors

  if ( collectTiming ) {
    Tau->SumFabX /= numDetectors;
    Tau->Fab2F /= numDetectors;
    tic = XLALGetCPUTime();
  }

  if ( whatToCompute & FSTATQ_2F )
    {
      const REAL4 Ad = resamp->Mmunu.Ad;
      const REAL4 Bd = resamp->Mmunu.Bd;
      const REAL4 Cd = resamp->Mmunu.Cd;
      const REAL4 Ed = resamp->Mmunu.Ed;
      const REAL4 Dd_inv = 1.0f / resamp->Mmunu.Dd;
      for ( UINT4 k=0; k < numFreqBins; k++ )
        {
          Fstats->twoF[k] = XLALComputeFstatFromFaFb ( ws->Fa_k[k], ws->Fb_k[k], Ad, Bd, Cd, Ed, Dd_inv );
        }
    } // if FSTATQ_2F

  if ( collectTiming ) {
      toc = XLALGetCPUTime();
      Tau->Fab2F += ( toc - tic );
  }

  // Return F-atoms per detector
  if (whatToCompute & FSTATQ_ATOMS_PER_DET) {
    XLAL_ERROR(XLAL_EFAILED, "NOT implemented!");
  }

  // Return antenna-pattern matrix
  Fstats->Mmunu = resamp->Mmunu;

  // return per-detector antenna-pattern matrices
  for ( UINT4 X = 0; X < numDetectors; X ++ )
    {
      Fstats->MmunuX[X] = resamp->MmunuX[X];
    }

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

  if ( collectTiming )
    {
      tocEnd = XLALGetCPUTime();

      FstatTimingGeneric *tiGen = &(resamp->timingGeneric);
      FstatTimingResamp  *tiRS  = &(resamp->timingResamp);
      XLAL_CHECK ( numDetectors == tiGen->Ndet, XLAL_EINVAL, "Inconsistent number of detectors between XLALCreateSetup() [%d] and XLALComputeFstat() [%d]\n", tiGen->Ndet, numDetectors );

      Tau->Total = (tocEnd - ticStart);
      // rescale all relevant timings to per-detector
      Tau->Total /= numDetectors;
      Tau->Bary  /= numDetectors;
      Tau->Spin  /= numDetectors;
      Tau->FFT   /= numDetectors;
      Tau->Norm  /= numDetectors;
      Tau->Copy  /= numDetectors;
      REAL8 Tau_buffer = Tau->Bary;
      // compute generic F-stat timing model contributions
      UINT4 NFbin      = Fstats->numFreqBins;
      REAL8 tauF_eff   = Tau->Total / NFbin;
      REAL8 tauF_core  = (Tau->Total - Tau_buffer) / NFbin;

      // compute resampling timing model coefficients
      REAL8 tau0_Fbin  = (Tau->Copy + Tau->Norm + Tau->SumFabX + Tau->Fab2F) / NFbin;
      REAL8 tau0_spin  = Tau->Spin / (tiRS->Resolution * tiRS->NsampFFT );
      REAL8 tau0_FFT   = Tau->FFT / (5.0 * tiRS->NsampFFT * log2(tiRS->NsampFFT));

      // update the averaged timing-model quantities
      tiGen->NCalls ++;	// keep track of number of Fstat-calls for timing
#define updateAvgF(q) tiGen->q = ((tiGen->q *(tiGen->NCalls-1) + q)/(tiGen->NCalls))
      updateAvgF(tauF_eff);
      updateAvgF(tauF_core);
      // we also average NFbin, which can be different between different calls to XLALComputeFstat() (contrary to Ndet)
      updateAvgF(NFbin);

#define updateAvgRS(q) tiRS->q = ((tiRS->q *(tiGen->NCalls-1) + q)/(tiGen->NCalls))
      updateAvgRS(tau0_Fbin);
      updateAvgRS(tau0_spin);
      updateAvgRS(tau0_FFT);

      // buffer-quantities only updated if buffer was actually recomputed
      if ( Tau->BufferRecomputed )
        {
          REAL8 tau0_bary   = Tau_buffer / (tiRS->Resolution * tiRS->NsampFFT);
          REAL8 tauF_buffer = Tau_buffer / NFbin;

          updateAvgF(tauF_buffer);
          updateAvgRS(tau0_bary);
        } // if BufferRecomputed

    } // if collectTiming

  return XLAL_SUCCESS;

} // XLALComputeFstatResamp()


static int
XLALComputeFaFb_Resamp ( ResampMethodData *resamp,				//!< [in,out] buffered resampling data and workspace
                         ResampWorkspace *ws,					//!< [in,out] resampling workspace (memory-sharing across segments)
                         const PulsarDopplerParams thisPoint,			//!< [in] Doppler point to compute {FaX,FbX} for
                         REAL8 dFreq,						//!< [in] output frequency resolution
                         UINT4 numFreqBins,					//!< [in] number of output frequency bins
                         const COMPLEX8TimeSeries * restrict TimeSeries_SRC_a,	//!< [in] SRC-frame single-IFO timeseries * a(t)
                         const COMPLEX8TimeSeries * restrict TimeSeries_SRC_b	//!< [in] SRC-frame single-IFO timeseries * b(t)
                         )
{
  XLAL_CHECK ( (resamp != NULL) && (ws != NULL) && (TimeSeries_SRC_a != NULL) && (TimeSeries_SRC_b != NULL), XLAL_EINVAL );
  XLAL_CHECK ( dFreq > 0, XLAL_EINVAL );
  XLAL_CHECK ( numFreqBins <= ws->numFreqBinsAlloc, XLAL_EINVAL );

  REAL8 FreqOut0 = thisPoint.fkdot[0];

  // compute frequency shift to align heterodyne frequency with output frequency bins
  REAL8 fHet   = TimeSeries_SRC_a->f0;
  REAL8 dt_SRC = TimeSeries_SRC_a->deltaT;

  REAL8 dFreqFFT = dFreq / resamp->decimateFFT;	// internally may be using higher frequency resolution dFreqFFT than requested
  REAL8 freqShift = remainder ( FreqOut0 - fHet, dFreq ); // frequency shift to closest bin
  REAL8 fMinFFT = fHet + freqShift - dFreqFFT * (resamp->numSamplesFFT/2);	// we'll shift DC into the *middle bin* N/2  [N always even!]
  XLAL_CHECK ( FreqOut0 >= fMinFFT, XLAL_EDOM, "Lowest output frequency outside the available frequency band: [FreqOut0 = %.16g] < [fMinFFT = %.16g]\n", FreqOut0, fMinFFT );
  UINT4 offset_bins = (UINT4) lround ( ( FreqOut0 - fMinFFT ) / dFreqFFT );
  UINT4 maxOutputBin = offset_bins + (numFreqBins - 1) * resamp->decimateFFT;
  XLAL_CHECK ( maxOutputBin < resamp->numSamplesFFT, XLAL_EDOM, "Highest output frequency bin outside available band: [maxOutputBin = %d] >= [numSamplesFFT = %d]\n", maxOutputBin, resamp->numSamplesFFT );

  FstatTimingResamp *tiRS = &(resamp->timingResamp);
  BOOLEAN collectTiming = resamp->collectTiming;
  REAL8 tic = 0, toc = 0;

  XLAL_CHECK ( resamp->numSamplesFFT >= TimeSeries_SRC_a->data->length, XLAL_EFAILED, "[numSamplesFFT = %d] < [len(TimeSeries_SRC_a) = %d]\n", resamp->numSamplesFFT, TimeSeries_SRC_a->data->length );
  XLAL_CHECK ( resamp->numSamplesFFT >= TimeSeries_SRC_b->data->length, XLAL_EFAILED, "[numSamplesFFT = %d] < [len(TimeSeries_SRC_b) = %d]\n", resamp->numSamplesFFT, TimeSeries_SRC_b->data->length );

  if ( collectTiming ) {
    tic = XLALGetCPUTime();
  }
  memset ( ws->TS_FFT, 0, resamp->numSamplesFFT * sizeof(ws->TS_FFT[0]) );
  // ----- compute FaX_k
  // apply spindown phase-factors, store result in zero-padded timeseries for 'FFT'ing
  XLAL_CHECK ( XLALApplySpindownAndFreqShift ( ws->TS_FFT, TimeSeries_SRC_a, &thisPoint, freqShift ) == XLAL_SUCCESS, XLAL_EFUNC );

  if ( collectTiming ) {
    toc = XLALGetCPUTime();
    tiRS->Tau.Spin += ( toc - tic);
    tic = toc;
  }

  // Fourier transform the resampled Fa(t)
  fftwf_execute_dft ( resamp->fftplan, ws->TS_FFT, ws->FabX_Raw );

  if ( collectTiming ) {
    toc = XLALGetCPUTime();
    tiRS->Tau.FFT += ( toc - tic);
    tic = toc;
  }

  for ( UINT4 k = 0; k < numFreqBins; k++ ) {
    ws->FaX_k[k] = ws->FabX_Raw [ offset_bins + k * resamp->decimateFFT ];
  }

  if ( collectTiming ) {
    toc = XLALGetCPUTime();
    tiRS->Tau.Copy += ( toc - tic);
    tic = toc;
  }

  // ----- compute FbX_k
  // apply spindown phase-factors, store result in zero-padded timeseries for 'FFT'ing
  XLAL_CHECK ( XLALApplySpindownAndFreqShift ( ws->TS_FFT, TimeSeries_SRC_b, &thisPoint, freqShift ) == XLAL_SUCCESS, XLAL_EFUNC );

  if ( collectTiming ) {
    toc = XLALGetCPUTime();
    tiRS->Tau.Spin += ( toc - tic);
    tic = toc;
  }

  // Fourier transform the resampled Fa(t)
  fftwf_execute_dft ( resamp->fftplan, ws->TS_FFT, ws->FabX_Raw );

  if ( collectTiming ) {
    toc = XLALGetCPUTime();
    tiRS->Tau.FFT += ( toc - tic);
    tic = toc;
  }

  for ( UINT4 k = 0; k < numFreqBins; k++ ) {
    ws->FbX_k[k] = ws->FabX_Raw [ offset_bins + k * resamp->decimateFFT ];
  }

  if ( collectTiming ) {
    toc = XLALGetCPUTime();
    tiRS->Tau.Copy += ( toc - tic);
    tic = toc;
  }

  // ----- normalization factors to be applied to Fa and Fb:
  const REAL8 dtauX = GPSDIFF ( TimeSeries_SRC_a->epoch, thisPoint.refTime );
  for ( UINT4 k = 0; k < numFreqBins; k++ )
    {
      REAL8 f_k = FreqOut0 + k * dFreq;
      REAL8 cycles = - f_k * dtauX;
      REAL4 sinphase, cosphase;
      XLALSinCos2PiLUT ( &sinphase, &cosphase, cycles );
      COMPLEX8 normX_k = dt_SRC * crectf ( cosphase, sinphase );
      ws->FaX_k[k] *= normX_k;
      ws->FbX_k[k] *= normX_k;
    } // for k < numFreqBinsOut

  if ( collectTiming ) {
    toc = XLALGetCPUTime();
    tiRS->Tau.Norm += ( toc - tic);
    tic = toc;
  }

  return XLAL_SUCCESS;

} // XLALComputeFaFb_Resamp()

static int
XLALApplySpindownAndFreqShift ( COMPLEX8 *restrict xOut,      			///< [out] the spindown-corrected SRC-frame timeseries
                                const COMPLEX8TimeSeries *restrict xIn,		///< [in] the input SRC-frame timeseries
                                const PulsarDopplerParams *restrict doppler,	///< [in] containing spindown parameters
                                REAL8 freqShift					///< [in] frequency-shift to apply, sign is "new - old"
                                )
{
  // input sanity checks
  XLAL_CHECK ( xOut != NULL, XLAL_EINVAL );
  XLAL_CHECK ( xIn != NULL, XLAL_EINVAL );
  XLAL_CHECK ( doppler != NULL, XLAL_EINVAL );

  // determine number of spin downs to include
  UINT4 s_max = PULSAR_MAX_SPINS - 1;
  while ( (s_max > 0) && (doppler->fkdot[s_max] == 0) ) {
    s_max --;
  }

  REAL8 dt = xIn->deltaT;
  UINT4 numSamplesIn  = xIn->data->length;

  LIGOTimeGPS epoch = xIn->epoch;
  REAL8 Dtau0 = GPSDIFF ( epoch, doppler->refTime );

  // loop over time samples
  for ( UINT4 j = 0; j < numSamplesIn; j ++ )
    {
      REAL8 taup_j = j * dt;
      REAL8 Dtau_alpha_j = Dtau0 + taup_j;

      REAL8 cycles = - freqShift * taup_j;

      REAL8 Dtau_pow_kp1 = Dtau_alpha_j;
      for ( UINT4 k = 1; k <= s_max; k++ )
        {
          Dtau_pow_kp1 *= Dtau_alpha_j;
          cycles += - LAL_FACT_INV[k+1] * doppler->fkdot[k] * Dtau_pow_kp1;
        } // for k = 1 ... s_max

      REAL4 cosphase, sinphase;
      XLAL_CHECK( XLALSinCos2PiLUT ( &sinphase, &cosphase, cycles ) == XLAL_SUCCESS, XLAL_EFUNC );
      COMPLEX8 em2piphase = crectf ( cosphase, sinphase );

      // weight the complex timeseries by the antenna patterns
      xOut[j] = em2piphase * xIn->data->data[j];

    } // for j < numSamplesIn

  return XLAL_SUCCESS;

} // XLALApplySpindownAndFreqShift()

///
/// Performs barycentric resampling on a multi-detector timeseries, updates resampling buffer with results
///
/// NOTE Buffering: this function does check
/// 1) whether the previously-buffered solution can be completely reused (same sky-position and binary parameters), or
/// 2) if at least sky-dependent quantities can be re-used (antenna-patterns + timings) in case only binary parameters changed
///
static int
XLALBarycentricResampleMultiCOMPLEX8TimeSeries ( ResampMethodData *resamp,		// [in/out] resampling input and buffer (to store resampling TS)
                                                 const PulsarDopplerParams *thisPoint,	// [in] current skypoint and reftime
                                                 const FstatCommon *common		// [in] various input quantities and parameters used here
                                                 )
{
  // check input sanity
  XLAL_CHECK ( thisPoint != NULL, XLAL_EINVAL );
  XLAL_CHECK ( common != NULL, XLAL_EINVAL );
  XLAL_CHECK ( resamp != NULL, XLAL_EINVAL );
  XLAL_CHECK ( resamp->multiTimeSeries_DET != NULL, XLAL_EINVAL );
  XLAL_CHECK ( resamp->multiTimeSeries_SRC_a != NULL, XLAL_EINVAL );
  XLAL_CHECK ( resamp->multiTimeSeries_SRC_b != NULL, XLAL_EINVAL );

  ResampWorkspace *ws = (ResampWorkspace*) common->workspace;

  UINT4 numDetectors = resamp->multiTimeSeries_DET->length;
  XLAL_CHECK ( resamp->multiTimeSeries_SRC_a->length == numDetectors, XLAL_EINVAL, "Inconsistent number of detectors tsDET(%d) != tsSRC(%d)\n", numDetectors, resamp->multiTimeSeries_SRC_a->length );
  XLAL_CHECK ( resamp->multiTimeSeries_SRC_b->length == numDetectors, XLAL_EINVAL, "Inconsistent number of detectors tsDET(%d) != tsSRC(%d)\n", numDetectors, resamp->multiTimeSeries_SRC_b->length );

  // ============================== BEGIN: handle buffering =============================
  BOOLEAN same_skypos = (resamp->prev_doppler.Alpha == thisPoint->Alpha) && (resamp->prev_doppler.Delta == thisPoint->Delta);
  BOOLEAN same_refTime = ( GPSDIFF ( resamp->prev_doppler.refTime, thisPoint->refTime ) == 0 );
  BOOLEAN same_binary = \
    (resamp->prev_doppler.asini == thisPoint->asini) &&
    (resamp->prev_doppler.period == thisPoint->period) &&
    (resamp->prev_doppler.ecc == thisPoint->ecc) &&
    (GPSDIFF( resamp->prev_doppler.tp, thisPoint->tp ) == 0 ) &&
    (resamp->prev_doppler.argp == thisPoint->argp);

  Timings_t *Tau = &(resamp->timingResamp.Tau);
  REAL8 tic, toc;
  BOOLEAN collectTiming = resamp->collectTiming;

  // if same sky-position *and* same binary, we can simply return as there's nothing to be done here
  if ( same_skypos && same_refTime && same_binary ) {
    Tau->BufferRecomputed = 0;
    return XLAL_SUCCESS;
  }
  // else: keep track of 'buffer miss', ie we need to recompute the buffer
  Tau->BufferRecomputed = 1;
  resamp->timingGeneric.NBufferMisses ++;
  if ( collectTiming ) {
    tic = XLALGetCPUTime();
  }

  MultiSSBtimes *multiSRCtimes = NULL;

  // only if different sky-position: re-compute antenna-patterns and SSB timings, re-use from buffer otherwise
  if ( ! ( same_skypos && same_refTime ) )
    {
      SkyPosition skypos;
      skypos.system = COORDINATESYSTEM_EQUATORIAL;
      skypos.longitude = thisPoint->Alpha;
      skypos.latitude  = thisPoint->Delta;

      XLALDestroyMultiAMCoeffs ( resamp->multiAMcoef );
      XLAL_CHECK ( (resamp->multiAMcoef = XLALComputeMultiAMCoeffs ( common->multiDetectorStates, common->multiNoiseWeights, skypos )) != NULL, XLAL_EFUNC );
      resamp->Mmunu = resamp->multiAMcoef->Mmunu;
      for ( UINT4 X = 0; X < numDetectors; X ++ )
        {
          resamp->MmunuX[X].Ad = resamp->multiAMcoef->data[X]->A;
          resamp->MmunuX[X].Bd = resamp->multiAMcoef->data[X]->B;
          resamp->MmunuX[X].Cd = resamp->multiAMcoef->data[X]->C;
          resamp->MmunuX[X].Ed = 0;
          resamp->MmunuX[X].Dd = resamp->multiAMcoef->data[X]->D;
        }

      XLALDestroyMultiSSBtimes ( resamp->multiSSBtimes );
      XLAL_CHECK ( (resamp->multiSSBtimes = XLALGetMultiSSBtimes ( common->multiDetectorStates, skypos, thisPoint->refTime, common->SSBprec )) != NULL, XLAL_EFUNC );

    } // if cannot re-use buffered solution ie if !(same_skypos && same_binary)

  if ( thisPoint->asini > 0 ) { // binary case
    XLAL_CHECK ( XLALAddMultiBinaryTimes ( &resamp->multiBinaryTimes, resamp->multiSSBtimes, thisPoint ) == XLAL_SUCCESS, XLAL_EFUNC );
    multiSRCtimes = resamp->multiBinaryTimes;
  } else { // isolated case
    multiSRCtimes = resamp->multiSSBtimes;
  }

  // record barycenter parameters in order to allow re-usal of this result ('buffering')
  resamp->prev_doppler = (*thisPoint);

  // shorthands
  REAL8 fHet = resamp->multiTimeSeries_DET->data[0]->f0;
  REAL8 Tsft = common->multiTimestamps->data[0]->deltaT;
  REAL8 dt_SRC = resamp->multiTimeSeries_SRC_a->data[0]->deltaT;

  const REAL4 signumLUT[2] = {1, -1};

  // loop over detectors X
  for ( UINT4 X = 0; X < numDetectors; X++)
    {
      // shorthand pointers: input
      const COMPLEX8TimeSeries *TimeSeries_DETX = resamp->multiTimeSeries_DET->data[X];
      const LIGOTimeGPSVector  *Timestamps_DETX = common->multiTimestamps->data[X];
      const SSBtimes *SRCtimesX                 = multiSRCtimes->data[X];
      const AMCoeffs *AMcoefX			= resamp->multiAMcoef->data[X];

      // shorthand pointers: output
      COMPLEX8TimeSeries *TimeSeries_SRCX_a     = resamp->multiTimeSeries_SRC_a->data[X];
      COMPLEX8TimeSeries *TimeSeries_SRCX_b     = resamp->multiTimeSeries_SRC_b->data[X];
      REAL8Vector *ti_DET = ws->SRCtimes_DET;

      // useful shorthands
      REAL8 refTime8        = GPSGETREAL8 ( &SRCtimesX->refTime );
      UINT4 numSFTsX        = Timestamps_DETX->length;
      UINT4 numSamples_DETX = TimeSeries_DETX->data->length;
      UINT4 numSamples_SRCX = TimeSeries_SRCX_a->data->length;

      // sanity checks on input data
      XLAL_CHECK ( numSamples_SRCX == TimeSeries_SRCX_b->data->length, XLAL_EINVAL );
      XLAL_CHECK ( dt_SRC == TimeSeries_SRCX_a->deltaT, XLAL_EINVAL );
      XLAL_CHECK ( dt_SRC == TimeSeries_SRCX_b->deltaT, XLAL_EINVAL );
      XLAL_CHECK ( numSamples_DETX > 0, XLAL_EINVAL, "Input timeseries for detector X=%d has zero samples. Can't handle that!\n", X );
      XLAL_CHECK ( (SRCtimesX->DeltaT->length == numSFTsX) && (SRCtimesX->Tdot->length == numSFTsX), XLAL_EINVAL );
      XLAL_CHECK ( fHet == resamp->multiTimeSeries_DET->data[X]->f0, XLAL_EINVAL );
      XLAL_CHECK ( Tsft == common->multiTimestamps->data[X]->deltaT, XLAL_EINVAL );

      TimeSeries_SRCX_a->f0 = fHet;
      TimeSeries_SRCX_b->f0 = fHet;
      // set SRC-frame time-series start-time
      REAL8 tStart_SRC_0 = refTime8 + SRCtimesX->DeltaT->data[0] - (0.5*Tsft) * SRCtimesX->Tdot->data[0];
      LIGOTimeGPS epoch;
      GPSSETREAL8 ( epoch, tStart_SRC_0 );
      TimeSeries_SRCX_a->epoch = epoch;
      TimeSeries_SRCX_b->epoch = epoch;

      // make sure all output samples are initialized to zero first, in case of gaps
      memset ( TimeSeries_SRCX_a->data->data, 0, TimeSeries_SRCX_a->data->length * sizeof(TimeSeries_SRCX_a->data->data[0]) );
      memset ( TimeSeries_SRCX_b->data->data, 0, TimeSeries_SRCX_b->data->length * sizeof(TimeSeries_SRCX_b->data->data[0]) );
      // make sure detector-frame timesteps to interpolate to are initialized to 0, in case of gaps
      memset ( ws->SRCtimes_DET->data, 0, ws->SRCtimes_DET->length * sizeof(ws->SRCtimes_DET->data[0]) );

      memset ( ws->TStmp1_SRC->data, 0, ws->TStmp1_SRC->length * sizeof(ws->TStmp1_SRC->data[0]) );
      memset ( ws->TStmp2_SRC->data, 0, ws->TStmp2_SRC->length * sizeof(ws->TStmp2_SRC->data[0]) );

      REAL8 tStart_DET_0 = GPSGETREAL8 ( &(Timestamps_DETX->data[0]) );// START time of the SFT at the detector

      // loop over SFT timestamps and compute the detector frame time samples corresponding to uniformly sampled SRC time samples
      for ( UINT4 alpha = 0; alpha < numSFTsX; alpha ++ )
        {
          // define some useful shorthands
          REAL8 Tdot_al       = SRCtimesX->Tdot->data [ alpha ];		// the instantaneous time derivitive dt_SRC/dt_DET at the MID-POINT of the SFT
          REAL8 tMid_SRC_al   = refTime8 + SRCtimesX->DeltaT->data[alpha];	// MID-POINT time of the SFT at the SRC
          REAL8 tStart_SRC_al = tMid_SRC_al - 0.5 * Tsft * Tdot_al;		// approximate START time of the SFT at the SRC
          REAL8 tEnd_SRC_al   = tMid_SRC_al + 0.5 * Tsft * Tdot_al;		// approximate END time of the SFT at the SRC

          REAL8 tStart_DET_al = GPSGETREAL8 ( &(Timestamps_DETX->data[alpha]) );// START time of the SFT at the detector
          REAL8 tMid_DET_al   = tStart_DET_al + 0.5 * Tsft;			// MID-POINT time of the SFT at the detector

          // indices of first and last SRC-frame sample corresponding to this SFT
          UINT4 iStart_SRC_al = lround ( (tStart_SRC_al - tStart_SRC_0) / dt_SRC );	// the index of the resampled timeseries corresponding to the start of the SFT
          UINT4 iEnd_SRC_al   = lround ( (tEnd_SRC_al - tStart_SRC_0) / dt_SRC );	// the index of the resampled timeseries corresponding to the end of the SFT

          // truncate to actual SRC-frame timeseries
          iStart_SRC_al = MYMIN ( iStart_SRC_al, numSamples_SRCX - 1);
          iEnd_SRC_al   = MYMIN ( iEnd_SRC_al, numSamples_SRCX - 1);
          UINT4 numSamplesSFT_SRC_al = iEnd_SRC_al - iStart_SRC_al + 1;		// the number of samples in the SRC-frame for this SFT

          REAL4 a_al = AMcoefX->a->data[alpha];
          REAL4 b_al = AMcoefX->b->data[alpha];
          for ( UINT4 j = 0; j < numSamplesSFT_SRC_al; j++ )
            {
              UINT4 iSRC_al_j  = iStart_SRC_al + j;

              // for each time sample in the SRC frame, we estimate the corresponding detector time,
              // using a linear approximation expanding around the midpoint of each SFT
              REAL8 t_SRC = tStart_SRC_0 + iSRC_al_j * dt_SRC;
              ti_DET->data [ iSRC_al_j ] = tMid_DET_al + ( t_SRC - tMid_SRC_al ) / Tdot_al;

              // pre-compute correction factors due to non-zero heterodyne frequency of input
              REAL8 tDiff = iSRC_al_j * dt_SRC + (tStart_DET_0 - ti_DET->data [ iSRC_al_j ]); 	// tSRC_al_j - tDET(tSRC_al_j)
              REAL8 cycles = fmod ( fHet * tDiff, 1.0 );				// the accumulated heterodyne cycles

              // use a look-up-table for speed to compute real and imaginary phase
              REAL4 cosphase, sinphase;                                   // the real and imaginary parts of the phase correction
              XLAL_CHECK( XLALSinCos2PiLUT ( &sinphase, &cosphase, -cycles ) == XLAL_SUCCESS, XLAL_EFUNC );
              COMPLEX8 ei2piphase = crectf ( cosphase, sinphase );

              // apply AM coefficients a(t), b(t) to SRC frame timeseries [alternate sign to get final FFT return DC in the middle]
              REAL4 signum = signumLUT [ (iSRC_al_j % 2) ];	// alternating sign, avoid branching
              ei2piphase *= signum;
              ws->TStmp1_SRC->data [ iSRC_al_j ] = ei2piphase * a_al;
              ws->TStmp2_SRC->data [ iSRC_al_j ] = ei2piphase * b_al;
            } // for j < numSamples_SRC_al

        } // for  alpha < numSFTsX

      XLAL_CHECK ( ti_DET->length >= TimeSeries_SRCX_a->data->length, XLAL_EINVAL );
      UINT4 bak_length = ti_DET->length;
      ti_DET->length = TimeSeries_SRCX_a->data->length;
      XLAL_CHECK ( XLALSincInterpolateCOMPLEX8TimeSeries ( TimeSeries_SRCX_a->data, ti_DET, TimeSeries_DETX, resamp->Dterms ) == XLAL_SUCCESS, XLAL_EFUNC );
      ti_DET->length = bak_length;

      // apply heterodyne correction and AM-functions a(t) and b(t) to interpolated timeseries
      for ( UINT4 j = 0; j < numSamples_SRCX; j ++ )
        {
          TimeSeries_SRCX_b->data->data[j] = TimeSeries_SRCX_a->data->data[j] * ws->TStmp2_SRC->data[j];
          TimeSeries_SRCX_a->data->data[j] *= ws->TStmp1_SRC->data[j];
        } // for j < numSamples_SRCX

    } // for X < numDetectors

  if ( collectTiming ) {
    toc = XLALGetCPUTime();
    Tau->Bary = (toc-tic);
  }

  return XLAL_SUCCESS;

} // XLALBarycentricResampleMultiCOMPLEX8TimeSeries()

static void
XLALGetFFTPlanHints ( int * planMode,
                      double * planGenTimeoutSeconds
                      )
{
  char * planMode_env = getenv("LAL_FSTAT_FFT_PLAN_MODE");
  char * planGenTimeout_env = getenv("LAL_FSTAT_FFT_PLAN_TIMEOUT");;
  int fft_plan_flags=FFTW_MEASURE;
  double fft_plan_timeout= FFTW_NO_TIMELIMIT ;

  if ( planGenTimeout_env ) {
    char * end;
    fft_plan_timeout=strtod(planGenTimeout_env,& end);
    if(end[0] != '\0') {
      fft_plan_timeout=FFTW_NO_TIMELIMIT;
    }
  }

  if ( planMode_env ) {
    if ( strcmp(planMode_env , "ESTIMATE" ) == 0 ) {
      fft_plan_flags=FFTW_ESTIMATE;
    }

    if ( strcmp(planMode_env , "MEASURE" ) == 0 ) {
      fft_plan_flags=FFTW_MEASURE;
    }

    if ( strcmp(planMode_env , "PATIENT" ) == 0 ) {
      fft_plan_flags=FFTW_PATIENT;
    }
  }
  *planMode=fft_plan_flags;
  *planGenTimeoutSeconds=fft_plan_timeout;
} // XLALGetFFTPlanHints

int
XLALExtractResampledTimeseries_intern ( MultiCOMPLEX8TimeSeries **multiTimeSeries_SRC_a, MultiCOMPLEX8TimeSeries **multiTimeSeries_SRC_b, const void* method_data )
{
  XLAL_CHECK ( method_data != NULL, XLAL_EINVAL );
  XLAL_CHECK ( ( multiTimeSeries_SRC_a != NULL ) && ( multiTimeSeries_SRC_b != NULL ) , XLAL_EINVAL );
  XLAL_CHECK ( method_data != NULL, XLAL_EINVAL );

  const ResampMethodData *resamp = (const ResampMethodData *) method_data;
  *multiTimeSeries_SRC_a = resamp->multiTimeSeries_SRC_a;
  *multiTimeSeries_SRC_b = resamp->multiTimeSeries_SRC_b;

  return XLAL_SUCCESS;

} // XLALExtractResampledTimeseries_intern()

int
XLALGetFstatTiming_Resamp ( const void *method_data, FstatTimingGeneric *timingGeneric, FstatTimingModel *timingModel )
{
  XLAL_CHECK ( method_data != NULL, XLAL_EINVAL );
  XLAL_CHECK ( timingGeneric != NULL, XLAL_EINVAL );
  XLAL_CHECK ( timingModel != NULL, XLAL_EINVAL );

  const ResampMethodData *resamp = (const ResampMethodData*) method_data;
  XLAL_CHECK ( resamp != NULL, XLAL_EINVAL );

  (*timingGeneric) = resamp->timingGeneric; // struct-copy generic timing measurements

  const FstatTimingResamp *tiRS = &(resamp->timingResamp);

  // return method-specific timing model values
  XLAL_INIT_MEM( (*timingModel) );

  UINT4 i = 0;
  timingModel->names[i]  = "NsampFFT0";
  timingModel->values[i] = tiRS->NsampFFT0;

  i++;
  timingModel->names[i]  = "NsampFFT";
  timingModel->values[i] = tiRS->NsampFFT;

  i++;
  timingModel->names[i]  = "Resolution";
  timingModel->values[i] = tiRS->Resolution;

  i++;
  timingModel->names[i]  = "tau0_Fbin";
  timingModel->values[i] = tiRS->tau0_Fbin;

  i++;
  timingModel->names[i]  = "tau0_spin";
  timingModel->values[i] = tiRS->tau0_spin;

  i++;
  timingModel->names[i]  = "tau0_FFT";
  timingModel->values[i] = tiRS->tau0_FFT;

  i++;
  timingModel->names[i]  = "tau0_bary";
  timingModel->values[i] = tiRS->tau0_bary;

  timingModel->numVariables = i+1;
  timingModel->help      = FstatTimingResampHelp;

  return XLAL_SUCCESS;
} // XLALGetFstatTiming_Resamp()
