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
typedef struct tagTimings_t
{
  REAL8 Total;		// total time spent in XLALComputeFstatResamp()
  REAL8 Bary;		// time spent in barycentric resampling
  REAL8 Spin;		// time spent in spindown+frequency correction
  REAL8 FFT;		// time spent in FFT
  REAL8 Copy;		// time spent copying results from FFT to FabX
  REAL8 Norm;		// time spent normalizing the final Fa,Fb
  REAL8 Fab2F;		// time to compute Fstat from {Fa,Fb}
  REAL8 Mem;		// time to realloc and Memset-0 arrays
  REAL8 SumFabX;	// time to sum_X Fab^X
} Timings_t;

typedef struct tagResampTimingInfo
{ // NOTE: all times refer to a single-detector timing case
  BOOLEAN collectTiming;	// turn on/off the collection of F-stat-method-specific timing-data (stored in workspace)

  UINT4 NFbin;		// number of frequency bins to compute F-stat for
  UINT4 Nsft;		// total number of SFTs used
  UINT4 Ndet;		// number of detectors
  REAL8 Resolution;	// frequency resolution in 'natural' units of 1/Tspan: R = T_span/T_FFT in [0, 1]
  UINT4 NsFFT0;		// 'original' length of barycentered timeseries to be FFT'ed
  UINT4 NsFFT;		// actual length of FFT, potentially rounded up to power-of-2 for efficiency

  Timings_t Tau;

  // ------------------------------------------------------------
  // Resampling timing model in terms of 'fundamental' coefficients tau_Fbin, tau_spin, tau_FFT:
  // ==> tau_RS = tau_Fbin + (NsFFT / NFbin) * [ R * tau_spin + tau_FFT ]
  // The total runtime per call (over NFbin for one detector) will generally have an additional contribution from barycentering:
  // ==> Tau.Total = NFbin * tau_RS + b * R * NsFFT * tau_bary,
  // where b = 1/N_{f1dot,f2dot,...} is the buffering weight applied to the barycentering contribution, which generally
  // will make b<<1 if there's many spindown template per sky- and binary-orbital templates.
  REAL8 tau_RS;		// measured resampling F-stat time per output F-stat bin (excluding barycentering): = (Tau.Total - Tau.Bary) / NFbin
  REAL8 tau_Fbin;	// time contribution of operations that scale exactly with numFreqBins, expressible as
                        // tau_Fbin = (Tau.Copy + Tau.Norm + Tau.SumFabX + Tau.Fab2F) / NFbin
  REAL8 tau_FFT;	// FFT time per FFT sample: tau_FFT = Tau.FFT / NsFFT
  REAL8 tau_spin;	// time contribution from applying spindown-correction (and freq-shift), per SRC-frame sample: tau_spin = Tau.Spin /(R* NsFFT)
  REAL8 tau_bary;	// barycentering time per SRC-frame sample, assuming no buffering: tau_bary = Tau.Bary /(R * NsFFT)
  // ------------------------------------------------------------

} ResampTimingInfo;

// ----- workspace ----------
typedef struct tagResampWorkspace
{
  // intermediate quantities to interpolate and operate on SRC-frame timeseries
  COMPLEX8Vector *TStmp1_SRC;	// can hold a single-detector SRC-frame spindown-corrected timeseries [without zero-padding]
  COMPLEX8Vector *TStmp2_SRC;	// can hold a single-detector SRC-frame spindown-corrected timeseries [without zero-padding]
  REAL8Vector *SRCtimes_DET;	// holds uniformly-spaced SRC-frame timesteps translated into detector frame [for interpolation]

  // input padded timeseries ts(t) and output Fab(f) of length 'numSamplesFFT' and corresponding fftw plan
  UINT4 numSamplesFFT;		// allocated number of zero-padded SRC-frame time samples (related to dFreq)
  UINT4 decimateFFT;		// output every n-th frequency bin, with n>1 iff (dFreq > 1/Tspan), and was internally decreased by n
  fftwf_plan fftplan;		// buffer FFT plan for given numSamplesOut length
  COMPLEX8 *TS_FFT;		// zero-padded, spindown-corr SRC-frame TS
  COMPLEX8 *FabX_Raw;		// raw full-band FFT result Fa,Fb

  // arrays of size numFreqBinsOut over frequency bins f_k:
  UINT4 numFreqBinsOut;		// number of output frequency bins {f_k}
  COMPLEX8 *FaX_k;		// properly normalized F_a^X(f_k) over output bins
  COMPLEX8 *FbX_k;		// properly normalized F_b^X(f_k) over output bins
  COMPLEX8 *Fa_k;		// properly normalized F_a(f_k) over output bins
  COMPLEX8 *Fb_k;		// properly normalized F_b(f_k) over output bins
  UINT4 numFreqBinsAlloc;	// internal: keep track of allocated length of frequency-arrays

  ResampTimingInfo *timingInfo;	// pointer to storage for collecting timing data (which lives in ResampMethodData)
} ResampWorkspace;

typedef struct
{
  MultiCOMPLEX8TimeSeries  *multiTimeSeries_DET;	// input SFTs converted into a heterodyned timeseries
  // ----- buffering -----
  PulsarDopplerParams prev_doppler;			// buffering: previous phase-evolution ("doppler") parameters

  AntennaPatternMatrix Mmunu;				// combined multi-IFO antenna-pattern coefficients {A,B,C,E}
  AntennaPatternMatrix MmunuX[PULSAR_MAX_DETECTORS];	// per-IFO antenna-pattern coefficients {AX,BX,CX,EX}

  MultiCOMPLEX8TimeSeries *multiTimeSeries_SRC_a;	// multi-detector SRC-frame timeseries, multiplied by AM function a(t)
  MultiCOMPLEX8TimeSeries *multiTimeSeries_SRC_b;	// multi-detector SRC-frame timeseries, multiplied by AM function b(t)

  ResampTimingInfo timingInfo;				// temporary storage for collecting timing data
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
XLALComputeFaFb_Resamp ( ResampWorkspace *ws,
                         const PulsarDopplerParams thisPoint,
                         REAL8 dFreq,
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

  LAL_FFTW_WISDOM_LOCK;
  fftwf_destroy_plan ( ws->fftplan );
  LAL_FFTW_WISDOM_UNLOCK;

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

  // Set method function pointers
  funcs->compute_func = XLALComputeFstatResamp;
  funcs->method_data_destroy_func = XLALDestroyResampMethodData;
  funcs->workspace_destroy_func = XLALDestroyResampWorkspace;

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

  UINT4 numSamplesFFT0 = (UINT4) ceil ( TspanFFT / dt_DET );      // we use ceil() so that we artificially widen the band rather than reduce it
  UINT4 numSamplesFFT = (UINT4) pow ( 2, ceil ( log2 ( numSamplesFFT0 ) ) );  // round numSamplesFFT up to next power of 2 for most effiecient FFT
  REAL8 dt_SRC = TspanFFT / numSamplesFFT;			// adjust sampling rate to allow achieving exact requested dFreq=1/TspanFFT !

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
      if ( numSamplesFFT > ws->numSamplesFFT )
        {
          int fft_plan_flags=FFTW_MEASURE;
          double fft_plan_timeout= FFTW_NO_TIMELIMIT ;

          fftw_free ( ws->FabX_Raw );
          XLAL_CHECK ( (ws->FabX_Raw = fftw_malloc ( numSamplesFFT * sizeof(COMPLEX8) )) != NULL, XLAL_ENOMEM );
          fftw_free ( ws->TS_FFT );
          XLAL_CHECK ( (ws->TS_FFT   = fftw_malloc ( numSamplesFFT * sizeof(COMPLEX8) )) != NULL, XLAL_ENOMEM );

          LAL_FFTW_WISDOM_LOCK;
          XLALGetFFTPlanHints (& fft_plan_flags , & fft_plan_timeout);

          fftwf_destroy_plan ( ws->fftplan );

          fftw_set_timelimit( fft_plan_timeout );
          XLAL_CHECK ( (ws->fftplan = fftwf_plan_dft_1d ( numSamplesFFT, ws->TS_FFT, ws->FabX_Raw, FFTW_FORWARD, fft_plan_flags )) != NULL, XLAL_EFAILED, "fftwf_plan_dft_1d() failed\n");
          LAL_FFTW_WISDOM_UNLOCK;
          ws->numSamplesFFT = numSamplesFFT;
          ws->decimateFFT = decimateFFT;
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
      int fft_plan_flags=FFTW_MEASURE;
      double fft_plan_timeout= FFTW_NO_TIMELIMIT ;

      XLAL_CHECK ( (ws = XLALCalloc ( 1, sizeof(*ws))) != NULL, XLAL_ENOMEM );
      XLAL_CHECK ( (ws->TStmp1_SRC   = XLALCreateCOMPLEX8Vector ( numSamplesMax_SRC )) != NULL, XLAL_EFUNC );
      XLAL_CHECK ( (ws->TStmp2_SRC   = XLALCreateCOMPLEX8Vector ( numSamplesMax_SRC )) != NULL, XLAL_EFUNC );
      XLAL_CHECK ( (ws->SRCtimes_DET = XLALCreateREAL8Vector ( numSamplesMax_SRC )) != NULL, XLAL_EFUNC );

      XLAL_CHECK ( (ws->FabX_Raw = fftw_malloc ( numSamplesFFT * sizeof(COMPLEX8) )) != NULL, XLAL_ENOMEM );
      XLAL_CHECK ( (ws->TS_FFT   = fftw_malloc ( numSamplesFFT * sizeof(COMPLEX8) )) != NULL, XLAL_ENOMEM );

      LAL_FFTW_WISDOM_LOCK;
      XLALGetFFTPlanHints (& fft_plan_flags , & fft_plan_timeout);
      fftw_set_timelimit( fft_plan_timeout );
      XLAL_CHECK ( (ws->fftplan = fftwf_plan_dft_1d ( numSamplesFFT, ws->TS_FFT, ws->FabX_Raw, FFTW_FORWARD, fft_plan_flags )) != NULL, XLAL_EFAILED, "fftwf_plan_dft_1d() failed\n");
      LAL_FFTW_WISDOM_UNLOCK;
      ws->numSamplesFFT = numSamplesFFT;
      ws->decimateFFT = decimateFFT;

      common->workspace = ws;
    } // end: if we create our own workspace

  // initialize struct for collecting timing data, store invariant 'meta' quantities about this setup
  XLAL_INIT_MEM ( resamp->timingInfo );
  resamp->timingInfo.collectTiming = optArgs->collectTiming;	// whether or not to collect timing info
  resamp->timingInfo.Resolution = TspanXMax / TspanFFT;
  resamp->timingInfo.NsFFT0 = numSamplesFFT0;
  resamp->timingInfo.NsFFT  = numSamplesFFT;

  resamp->timingInfo.Ndet = numDetectors;
  UINT4 numSFTs = 0;
  for ( UINT4 X = 0; X < numDetectors; X ++ ) {
    numSFTs += common->multiDetectorStates->data[X]->length;
  }
  resamp->timingInfo.Nsft = numSFTs;

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
  ResampTimingInfo *ti = &(resamp->timingInfo);
  REAL8 ticStart = 0, tocEnd = 0;
  REAL8 tic = 0, toc = 0;
  if ( ti->collectTiming ) {
    XLAL_INIT_MEM ( ti->Tau );	// re-set all timings to 0 at beginning of each Fstat-call
    ti->NFbin = Fstats->numFreqBins;

    ticStart = XLALGetCPUTime();
  }
  // store pointer to timing-info storage in workspace (for use by XLALComputeFaFb_Resamp() )
  ws->timingInfo = ti;

  // ============================== BEGIN: handle buffering =============================
  BOOLEAN same_skypos = (resamp->prev_doppler.Alpha == thisPoint.Alpha) && (resamp->prev_doppler.Delta == thisPoint.Delta);
  BOOLEAN same_refTime = ( GPSDIFF ( resamp->prev_doppler.refTime, thisPoint.refTime ) == 0 );
  BOOLEAN same_binary = \
    (resamp->prev_doppler.asini == thisPoint.asini) &&
    (resamp->prev_doppler.period == thisPoint.period) &&
    (resamp->prev_doppler.ecc == thisPoint.ecc) &&
    (GPSDIFF( resamp->prev_doppler.tp, thisPoint.tp ) == 0 ) &&
    (resamp->prev_doppler.argp == thisPoint.argp);

  // ----- not same skypos+binary+refTime? --> re-compute SRC-frame timeseries, AM-coeffs and store in buffer
  if ( ti->collectTiming ) {
    tic = XLALGetCPUTime();
  }

  if ( ! ( same_skypos && same_refTime && same_binary) )
    {
      XLAL_CHECK ( XLALBarycentricResampleMultiCOMPLEX8TimeSeries ( resamp, &thisPoint, common ) == XLAL_SUCCESS, XLAL_EFUNC );
      // record barycenter parameters in order to allow re-usal of this result ('buffering')
      resamp->prev_doppler = thisPoint;
    }

  if ( ti->collectTiming ) {
    toc = XLALGetCPUTime();
    ti->Tau.Bary = (toc-tic);
  }
  MultiCOMPLEX8TimeSeries *multiTimeSeries_SRC_a = resamp->multiTimeSeries_SRC_a;
  MultiCOMPLEX8TimeSeries *multiTimeSeries_SRC_b = resamp->multiTimeSeries_SRC_b;

  // ============================== check workspace is properly allocated and initialized ===========

  // ----- workspace that depends on number of output frequency bins 'numFreqBins' ----------
  UINT4 numFreqBins = Fstats->numFreqBins;

  if ( ti->collectTiming ) {
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
  ws->numFreqBinsOut = numFreqBins;

  if ( ti->collectTiming ) {
    toc = XLALGetCPUTime();
    ti->Tau.Mem = (toc-tic);	// this one doesn't scale with number of detector!
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
      XLAL_CHECK ( XLALComputeFaFb_Resamp ( ws, thisPoint, common->dFreq, TimeSeriesX_SRC_a, TimeSeriesX_SRC_b ) == XLAL_SUCCESS, XLAL_EFUNC );

      if ( ti->collectTiming ) {
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

      if ( ti->collectTiming ) {
        toc = XLALGetCPUTime();
        ti->Tau.SumFabX += (toc-tic);
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

      if ( ti->collectTiming ) {
        toc = XLALGetCPUTime();
        ti->Tau.Fab2F += ( toc - tic );
      }

    } // for X < numDetectors

  if ( ti->collectTiming ) {
    ti->Tau.SumFabX /= numDetectors;
    ti->Tau.Fab2F /= numDetectors;
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

  if ( ti->collectTiming ) {
      toc = XLALGetCPUTime();
      ti->Tau.Fab2F += ( toc - tic );
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

  if ( ti->collectTiming ) {
    // timings are per-detector
    tocEnd = XLALGetCPUTime();
    Timings_t *Tau = &(ti->Tau);
    Tau->Total = (tocEnd - ticStart);
    // rescale all relevant timings to single-IFO case
    Tau->Total /= numDetectors;
    Tau->Bary  /= numDetectors;
    Tau->Spin  /= numDetectors;
    Tau->FFT   /= numDetectors;
    Tau->Norm  /= numDetectors;
    Tau->Copy  /= numDetectors;

    // compute 'fundamental' per output bin timing numbers and timing-model coefficients
    ti->tau_RS	= (Tau->Total - Tau->Bary) / ti->NFbin;
    ti->tau_FFT	= Tau->FFT / ti->NsFFT;
    ti->tau_spin= Tau->Spin / (ti->Resolution * ti->NsFFT );
    ti->tau_Fbin= (Tau->Copy + Tau->Norm + Tau->SumFabX + Tau->Fab2F) / ti->NFbin;
    ti->tau_bary= Tau->Bary / (ti->Resolution * ti->NsFFT);
  }

  return XLAL_SUCCESS;

} // XLALComputeFstatResamp()


static int
XLALComputeFaFb_Resamp ( ResampWorkspace *restrict ws,				//!< [in,out] pre-allocated 'workspace' for temporary and output quantities
                         const PulsarDopplerParams thisPoint,			//!< [in] Doppler point to compute {FaX,FbX} for
                         REAL8 dFreq,						//!< [in] output frequency resolution
                         const COMPLEX8TimeSeries * restrict TimeSeries_SRC_a,	//!< [in] SRC-frame single-IFO timeseries * a(t)
                         const COMPLEX8TimeSeries * restrict TimeSeries_SRC_b	//!< [in] SRC-frame single-IFO timeseries * b(t)
                         )
{
  XLAL_CHECK ( (ws != NULL) && (TimeSeries_SRC_a != NULL) && (TimeSeries_SRC_b != NULL), XLAL_EINVAL );
  XLAL_CHECK ( dFreq > 0, XLAL_EINVAL );

  REAL8 FreqOut0 = thisPoint.fkdot[0];

  // compute frequency shift to align heterodyne frequency with output frequency bins
  REAL8 fHet   = TimeSeries_SRC_a->f0;
  REAL8 dt_SRC = TimeSeries_SRC_a->deltaT;

  REAL8 dFreqFFT = dFreq / ws->decimateFFT;	// internally may be using higher frequency resolution dFreqFFT than requested
  REAL8 freqShift = remainder ( FreqOut0 - fHet, dFreq ); // frequency shift to closest bin
  REAL8 fMinFFT = fHet + freqShift - dFreqFFT * (ws->numSamplesFFT/2);	// we'll shift DC into the *middle bin* N/2  [N always even!]
  XLAL_CHECK ( FreqOut0 >= fMinFFT, XLAL_EDOM, "Lowest output frequency outside the available frequency band: [FreqOut0 = %.16g] < [fMinFFT = %.16g]\n", FreqOut0, fMinFFT );
  UINT4 offset_bins = (UINT4) lround ( ( FreqOut0 - fMinFFT ) / dFreqFFT );
  UINT4 maxOutputBin = offset_bins + (ws->numFreqBinsOut-1) * ws->decimateFFT;
  XLAL_CHECK ( maxOutputBin < ws->numSamplesFFT, XLAL_EDOM, "Highest output frequency bin outside available band: [maxOutputBin = %d] >= [numSamplesFFT = %d]\n", maxOutputBin, ws->numSamplesFFT );

  ResampTimingInfo *ti = ws->timingInfo;
  REAL8 tic = 0, toc = 0;

  XLAL_CHECK ( ws->numSamplesFFT >= TimeSeries_SRC_a->data->length, XLAL_EFAILED, "[numSamplesFFT = %d] < [len(TimeSeries_SRC_a) = %d]\n", ws->numSamplesFFT, TimeSeries_SRC_a->data->length );
  XLAL_CHECK ( ws->numSamplesFFT >= TimeSeries_SRC_b->data->length, XLAL_EFAILED, "[numSamplesFFT = %d] < [len(TimeSeries_SRC_b) = %d]\n", ws->numSamplesFFT, TimeSeries_SRC_b->data->length );

  if ( ti->collectTiming ) {
    tic = XLALGetCPUTime();
  }
  memset ( ws->TS_FFT, 0, ws->numSamplesFFT * sizeof(ws->TS_FFT[0]) );
  // ----- compute FaX_k
  // apply spindown phase-factors, store result in zero-padded timeseries for 'FFT'ing
  XLAL_CHECK ( XLALApplySpindownAndFreqShift ( ws->TS_FFT, TimeSeries_SRC_a, &thisPoint, freqShift ) == XLAL_SUCCESS, XLAL_EFUNC );

  if ( ti->collectTiming ) {
    toc = XLALGetCPUTime();
    ti->Tau.Spin += ( toc - tic);
    tic = toc;
  }

  // Fourier transform the resampled Fa(t)
  fftwf_execute ( ws->fftplan );

  if ( ti->collectTiming ) {
    toc = XLALGetCPUTime();
    ti->Tau.FFT += ( toc - tic);
    tic = toc;
  }

  for ( UINT4 k = 0; k < ws->numFreqBinsOut; k++ ) {
    ws->FaX_k[k] = ws->FabX_Raw [ offset_bins + k * ws->decimateFFT ];
  }

  if ( ti->collectTiming ) {
    toc = XLALGetCPUTime();
    ti->Tau.Copy += ( toc - tic);
    tic = toc;
  }

  // ----- compute FbX_k
  // apply spindown phase-factors, store result in zero-padded timeseries for 'FFT'ing
  XLAL_CHECK ( XLALApplySpindownAndFreqShift ( ws->TS_FFT, TimeSeries_SRC_b, &thisPoint, freqShift ) == XLAL_SUCCESS, XLAL_EFUNC );

  if ( ti->collectTiming ) {
    toc = XLALGetCPUTime();
    ti->Tau.Spin += ( toc - tic);
    tic = toc;
  }

  // Fourier transform the resampled Fa(t)
  fftwf_execute ( ws->fftplan );

  if ( ti->collectTiming ) {
    toc = XLALGetCPUTime();
    ti->Tau.FFT += ( toc - tic);
    tic = toc;
  }

  for ( UINT4 k = 0; k < ws->numFreqBinsOut; k++ ) {
    ws->FbX_k[k] = ws->FabX_Raw [ offset_bins + k * ws->decimateFFT ];
  }

  if ( ti->collectTiming ) {
    toc = XLALGetCPUTime();
    ti->Tau.Copy += ( toc - tic);
    tic = toc;
  }

  // ----- normalization factors to be applied to Fa and Fb:
  const REAL8 dtauX = GPSDIFF ( TimeSeries_SRC_a->epoch, thisPoint.refTime );
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

  if ( ti->collectTiming ) {
    toc = XLALGetCPUTime();
    ti->Tau.Norm += ( toc - tic);
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
/// NOTE: this function does NOT check whether the previously-buffered solution can be reused, it assumes the
/// caller has already done so, and simply computes the requested resampled time-series, and AM-coefficients
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

  SkyPosition skypos;
  skypos.system = COORDINATESYSTEM_EQUATORIAL;
  skypos.longitude = thisPoint->Alpha;
  skypos.latitude  = thisPoint->Delta;

  MultiAMCoeffs *multiAMcoef;
  XLAL_CHECK ( (multiAMcoef = XLALComputeMultiAMCoeffs ( common->multiDetectorStates, common->multiNoiseWeights, skypos )) != NULL, XLAL_EFUNC );
  resamp->Mmunu = multiAMcoef->Mmunu;
  for ( UINT4 X = 0; X < numDetectors; X ++ )
    {
      resamp->MmunuX[X].Ad = multiAMcoef->data[X]->A;
      resamp->MmunuX[X].Bd = multiAMcoef->data[X]->B;
      resamp->MmunuX[X].Cd = multiAMcoef->data[X]->C;
      resamp->MmunuX[X].Ed = 0;
      resamp->MmunuX[X].Dd = multiAMcoef->data[X]->D;
    }

  MultiSSBtimes *multiSRCtimes;
  XLAL_CHECK ( (multiSRCtimes = XLALGetMultiSSBtimes ( common->multiDetectorStates, skypos, thisPoint->refTime, common->SSBprec )) != NULL, XLAL_EFUNC );
  if ( thisPoint->asini > 0 ) {
    XLAL_CHECK ( XLALAddMultiBinaryTimes ( &multiSRCtimes, multiSRCtimes, thisPoint ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

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
      const AMCoeffs *AMcoefX			= multiAMcoef->data[X];

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

      const UINT4 Dterms = 8;
      XLAL_CHECK ( ti_DET->length >= TimeSeries_SRCX_a->data->length, XLAL_EINVAL );
      UINT4 bak_length = ti_DET->length;
      ti_DET->length = TimeSeries_SRCX_a->data->length;
      XLAL_CHECK ( XLALSincInterpolateCOMPLEX8TimeSeries ( TimeSeries_SRCX_a->data, ti_DET, TimeSeries_DETX, Dterms ) == XLAL_SUCCESS, XLAL_EFUNC );
      ti_DET->length = bak_length;

      // apply heterodyne correction and AM-functions a(t) and b(t) to interpolated timeseries
      for ( UINT4 j = 0; j < numSamples_SRCX; j ++ )
        {
          TimeSeries_SRCX_b->data->data[j] = TimeSeries_SRCX_a->data->data[j] * ws->TStmp2_SRC->data[j];
          TimeSeries_SRCX_a->data->data[j] *= ws->TStmp1_SRC->data[j];
        } // for j < numSamples_SRCX

    } // for X < numDetectors

  XLALDestroyMultiAMCoeffs ( multiAMcoef );
  XLALDestroyMultiSSBtimes ( multiSRCtimes );

  return XLAL_SUCCESS;

} // XLALBarycentricResampleMultiCOMPLEX8TimeSeries()

// export timing constants from internally-stored 'method data' struct
int
XLALGetFstatTiming_Resamp ( const void* method_data, REAL8 *tauF1Buf, REAL8 *tauF1NoBuf )
{
  XLAL_CHECK ( method_data != NULL, XLAL_EINVAL );
  XLAL_CHECK ( (tauF1Buf != NULL) && (tauF1NoBuf != NULL), XLAL_EINVAL );

  const ResampMethodData *resamp = (const ResampMethodData *)method_data;

  // time per template assuming perfect buffering
  (*tauF1Buf) = resamp->timingInfo.tau_RS;
  // 'raw' time per freq-bin, potentially including barycentering (if no buffering)
  (*tauF1NoBuf) = resamp->timingInfo.Tau.Total / resamp->timingInfo.NFbin;

  return XLAL_SUCCESS;

} // XLALGetFstatTiming_Resamp()

// append detailed (method-specific) timing info to given file
int
AppendFstatTimingInfo2File_Resamp ( const void* method_data, FILE *fp, BOOLEAN printHeader )
{
  XLAL_CHECK ( method_data != NULL, XLAL_EINVAL );
  XLAL_CHECK ( fp != NULL, XLAL_EINVAL );

  const ResampMethodData *resamp = (const ResampMethodData *)method_data;
  // print header if requested
  if ( printHeader ) {
    fprintf (fp, "%%%% ----- Resampling F-stat timing: -----\n");
    fprintf (fp, "%%%% Measured time (in seconds) per F-stat frequency bin per detector (excluding barycentering):\n");
    fprintf (fp, "%%%% tau_RS = (TauTotal - TauBary) / NFbin\n");
    fprintf (fp, "%%%% tau_RS-predicted = tau_Fbin + (NsFFT/NFbin) * ( R * tau_spin + tau_FFT )\n");
    fprintf (fp, "%%%% with the frequency resolution in natural units, R = Tspan / T_FFT = NsSRC / NsFFT,\n");
    fprintf (fp, "%%%% Total time per detector generally contains an additional barycentering contribution:\n");
    fprintf (fp, "%%%% TauTotal = NFbin * tauRS + b * R * NsFFT * tau_bary\n");
    fprintf (fp, "%%%% where the buffering weight b = 1/N_{f1dot,f2dot,..} goes to 0 for many spindowns per sky+binary template\n");

    fprintf (fp, "%%%%%8s %8s %8s %6s %6s", "NFbin", "NsFFT0", "l2NsFFT", "Ndet", "R" );
    fprintf (fp, " %10s %10s %10s %10s %10s %10s", "TauTotal", "tau_RS", "tau_Fbin", "tau_FFT", "tau_spin", "tau_bary" );
    fprintf (fp, "\n");
  }

  const ResampTimingInfo *ti = &(resamp->timingInfo);
  fprintf (fp, "%10d %8d %8.2f %6d %6.3f",
           ti->NFbin, ti->NsFFT0, log2(ti->NsFFT), ti->Ndet, ti->Resolution );

  fprintf (fp, " %10.1e %10.1e %10.1e %10.1e %10.1e %10.1e",
           ti->Tau.Total, ti->tau_RS, ti->tau_Fbin, ti->tau_FFT, ti->tau_spin, ti->tau_bary );

  fprintf (fp, "\n");

  return XLAL_SUCCESS;
} // AppendFstatTimingInfo2File_Resamp()



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

