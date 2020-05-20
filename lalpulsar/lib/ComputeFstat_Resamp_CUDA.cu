//
// Copyright (C) 2019, 2020 Liam Dunn
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
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cuComplex.h>
#include <cufft.h>

#include "ComputeFstat_internal.h"
#include "ComputeFstat_Resamp_internal.h"

#include <lal/Factorial.h>
#include <lal/LFTandTSutils.h>
#include <lal/LogPrintf.h>
#include <lal/SinCosLUT.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/Window.h>

///
/// \defgroup ComputeFstat_Resamp_CUDA_cu Module ComputeFstat_Resamp_CUDA.cu
/// \ingroup ComputeFstat_h
/// \brief Implements a CUDA version of the \a Resamp FFT-based resampling algorithm for
/// computing the \f$\mathcal{F}\f$-statistic \cite JKS98 .
///

// @{

// ========== Resamp internals ==========

// ----- local constants ----------

#define CUDA_BLOCK_SIZE 512

__constant__ REAL8 PI = M_PI;

__device__ __constant__ REAL8 lal_fact_inv[LAL_FACT_MAX];

// ----- local macros ----------

#define CPLX_MULT(x, y) (crectf(crealf(x)*crealf(y) - cimagf(x)*cimagf(y), cimagf(x)*crealf(y) + crealf(x)*cimagf(y)))

// ----- local types ----------

// ----- workspace ----------
typedef struct tagResampWorkspace
{
  // intermediate quantities to interpolate and operate on SRC-frame timeseries
  COMPLEX8Vector *TStmp1_SRC;	// can hold a single-detector SRC-frame spindown-corrected timeseries [without zero-padding]
  COMPLEX8Vector *TStmp2_SRC;	// can hold a single-detector SRC-frame spindown-corrected timeseries [without zero-padding]
  REAL8Vector *SRCtimes_DET;	// holds uniformly-spaced SRC-frame timesteps translated into detector frame [for interpolation]

  // input padded timeseries ts(t) and output Fab(f) of length 'numSamplesFFT'
  UINT4 numSamplesFFTAlloc;	// allocated number of zero-padded SRC-frame time samples (related to dFreq)
  cuComplex *TS_FFT;		// zero-padded, spindown-corr SRC-frame TS
  cuComplex *FabX_Raw;		// raw full-band FFT result Fa,Fb

  // arrays of size numFreqBinsOut over frequency bins f_k:
  cuComplex *FaX_k;		// properly normalized F_a^X(f_k) over output bins
  cuComplex *FbX_k;		// properly normalized F_b^X(f_k) over output bins
  cuComplex *Fa_k;		// properly normalized F_a(f_k) over output bins
  cuComplex *Fb_k;		// properly normalized F_b(f_k) over output bins
  UINT4 numFreqBinsAlloc;	// internal: keep track of allocated length of frequency-arrays

} ResampWorkspace;

typedef struct
{
  UINT4 Dterms;                                         // Number of terms to use (on either side) in Windowed-Sinc interpolation kernel
  MultiCOMPLEX8TimeSeries  *multiTimeSeries_DET;        // input SFTs converted into a heterodyned timeseries
  // ----- buffering -----
  PulsarDopplerParams prev_doppler;                     // buffering: previous phase-evolution ("doppler") parameters
  MultiAMCoeffs *multiAMcoef;                           // buffered antenna-pattern functions
  MultiSSBtimes *multiSSBtimes;                         // buffered SSB times, including *only* sky-position corrections, not binary
  MultiSSBtimes *multiBinaryTimes;                      // buffered SRC times, including both sky- and binary corrections [to avoid re-allocating this]

  AntennaPatternMatrix Mmunu;                           // combined multi-IFO antenna-pattern coefficients {A,B,C,E}
  AntennaPatternMatrix MmunuX[PULSAR_MAX_DETECTORS];    // per-IFO antenna-pattern coefficients {AX,BX,CX,EX}

  MultiCOMPLEX8TimeSeries *multiTimeSeries_SRC_a;       // multi-detector SRC-frame timeseries, multiplied by AM function a(t)
  MultiCOMPLEX8TimeSeries *multiTimeSeries_SRC_b;       // multi-detector SRC-frame timeseries, multiplied by AM function b(t)

  UINT4 numSamplesFFT;                                  // length of zero-padded SRC-frame timeseries (related to dFreq)
  UINT4 decimateFFT;                                    // output every n-th frequency bin, with n>1 iff (dFreq > 1/Tspan), and was internally decreased by n
  void *fftplan;                                        // FFT plan

  // ----- timing -----
  BOOLEAN collectTiming;                                // flag whether or not to collect timing information
  FstatTimingGeneric timingGeneric;                     // measured (generic) F-statistic timing values
  FstatTimingResamp  timingResamp;                      // measured Resamp-specific timing model data

} ResampMethodData;


// ----- local prototypes ----------

extern "C" int XLALSetupFstatResampCUDA ( void **method_data, FstatCommon *common, FstatMethodFuncs* funcs, MultiSFTVector *multiSFTs, const FstatOptionalArgs *optArgs );

static int XLALComputeFstatResampCUDA ( FstatResults* Fstats, const FstatCommon *common, void *method_data );
static int XLALApplySpindownAndFreqShiftCUDA ( cuComplex *xOut, const COMPLEX8TimeSeries *xIn, const PulsarDopplerParams *doppler, REAL8 freqShift );
static int XLALBarycentricResampleMultiCOMPLEX8TimeSeriesCUDA ( ResampMethodData *resamp, const PulsarDopplerParams *thisPoint, const FstatCommon *common );
static int XLALComputeFaFb_ResampCUDA ( ResampMethodData *resamp, ResampWorkspace *ws, const PulsarDopplerParams thisPoint, REAL8 dFreq, UINT4 numFreqBins, const COMPLEX8TimeSeries *TimeSeries_SRC_a, const COMPLEX8TimeSeries *TimeSeries_SRC_b );
static COMPLEX8Vector *CreateCOMPLEX8VectorCUDA(UINT4 length);
static REAL8Vector *CreateREAL8VectorCUDA(UINT4 length);
static void DestroyCOMPLEX8VectorCUDA(COMPLEX8Vector *vec);
static void DestroyREAL8VectorCUDA(REAL8Vector *vec);
static void DestroyCOMPLEX8TimeSeriesCUDA(COMPLEX8TimeSeries *series);
static void MoveCOMPLEX8TimeSeriesHtoD(COMPLEX8TimeSeries *series);
static void MoveMultiCOMPLEX8TimeSeriesHtoD(MultiCOMPLEX8TimeSeries *multi);
static void DestroyMultiCOMPLEX8TimeSeriesCUDA(MultiCOMPLEX8TimeSeries *multi);

// ==================== function definitions ====================

static COMPLEX8Vector *CreateCOMPLEX8VectorCUDA(UINT4 length)
{
  COMPLEX8Vector *vec;
  vec = (COMPLEX8Vector *)XLALMalloc(sizeof(COMPLEX8Vector));
  if(!vec)
    XLAL_ERROR_NULL(XLAL_ENOMEM);
  vec->length = length;
  cudaError_t err = cudaMalloc((void **)&vec->data, sizeof(COMPLEX8)*length);
  if(cudaSuccess != err)
    XLAL_ERROR_NULL(XLAL_ENOMEM);
  return vec;
}

static REAL8Vector *CreateREAL8VectorCUDA(UINT4 length)
{
  REAL8Vector *vec;
  vec = (REAL8Vector *)XLALMalloc(sizeof(REAL8Vector));
  if(!vec)
    XLAL_ERROR_NULL(XLAL_ENOMEM);
  vec->length = length;
  cudaError_t err = cudaMalloc((void **)&vec->data, sizeof(REAL8)*length);
  if(cudaSuccess != err)
    XLAL_ERROR_NULL(XLAL_ENOMEM);
  return vec;
}

static void DestroyCOMPLEX8VectorCUDA(COMPLEX8Vector *vec)
{
  if(!vec)
    return;
  if((!vec->length || !vec->data) && (vec->length || vec->data))
    XLAL_ERROR_VOID(XLAL_EINVAL);

  if(vec->data)
    cudaFree(vec->data);
  vec->data = NULL;
  XLALFree(vec);
}

static void DestroyREAL8VectorCUDA(REAL8Vector *vec)
{
  if(!vec)
    return;
  if((!vec->length || !vec->data) && (vec->length || vec->data))
    XLAL_ERROR_VOID(XLAL_EINVAL);

  if(vec->data)
    cudaFree(vec->data);
  vec->data = NULL;
  XLALFree(vec);
}

static void DestroyCOMPLEX8TimeSeriesCUDA(COMPLEX8TimeSeries *series)
{
    if(!series)
        return;

    DestroyCOMPLEX8VectorCUDA(series->data);
    XLALFree(series);
    return;
}

static void MoveCOMPLEX8TimeSeriesHtoD(COMPLEX8TimeSeries *series)
{
    COMPLEX8 *cpu_data = series->data->data;
    cudaError_t err = cudaMalloc((void **)&series->data->data, sizeof(COMPLEX8)*series->data->length);
    if(cudaSuccess != err)
        XLAL_ERROR_VOID(XLAL_ENOMEM);
    err = cudaMemcpy((void *)series->data->data, cpu_data, sizeof(COMPLEX8)*series->data->length, cudaMemcpyHostToDevice);
    if(cudaSuccess != err)
        XLAL_ERROR_VOID(XLAL_ENOMEM);
    XLALFree(cpu_data);
}

static void MoveMultiCOMPLEX8TimeSeriesHtoD(MultiCOMPLEX8TimeSeries *multi)
{
    for(UINT4 X = 0; X < multi->length; X++)
        MoveCOMPLEX8TimeSeriesHtoD(multi->data[X]);
}

static void DestroyMultiCOMPLEX8TimeSeriesCUDA(MultiCOMPLEX8TimeSeries *multi)
{
    if(!multi)
        return;
    if(multi->data != NULL)
    {
        UINT4 numDetectors = multi->length;
        for(UINT4 X = 0; X < numDetectors; X++)
        {
            DestroyCOMPLEX8TimeSeriesCUDA(multi->data[X]);
        }
        XLALFree(multi->data);
    }
    XLALFree(multi);
}

static void
XLALDestroyResampWorkspace ( void *workspace )
{
  ResampWorkspace *ws = (ResampWorkspace*) workspace;

  DestroyCOMPLEX8VectorCUDA ( ws->TStmp1_SRC );
  DestroyCOMPLEX8VectorCUDA ( ws->TStmp2_SRC );
  DestroyREAL8VectorCUDA ( ws->SRCtimes_DET );

  cudaFree ( ws->FabX_Raw );
  cudaFree ( ws->TS_FFT );

  XLALFree ( ws );
  return;

} // XLALDestroyResampWorkspace()

// ---------- internal functions ----------

static void
XLALDestroyResampMethodData ( void* method_data )
{

  ResampMethodData *resamp = (ResampMethodData*) method_data;

  DestroyMultiCOMPLEX8TimeSeriesCUDA (resamp->multiTimeSeries_DET );

  // ----- free buffer
  DestroyMultiCOMPLEX8TimeSeriesCUDA ( resamp->multiTimeSeries_SRC_a );
  DestroyMultiCOMPLEX8TimeSeriesCUDA ( resamp->multiTimeSeries_SRC_b );
  XLALDestroyMultiAMCoeffs ( resamp->multiAMcoef );
  XLALDestroyMultiSSBtimes ( resamp->multiSSBtimes );
  XLALDestroyMultiSSBtimes ( resamp->multiBinaryTimes );

  XLALFree ( resamp );

} // XLALDestroyResampMethodData()


extern "C" int
XLALSetupFstatResampCUDA ( void **method_data,
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
  *method_data = XLALCalloc(1, sizeof(ResampMethodData));
  ResampMethodData *resamp = (ResampMethodData *)*method_data;
  XLAL_CHECK( resamp != NULL, XLAL_ENOMEM );

  resamp->Dterms = optArgs->Dterms;

  // Set method function pointers
  funcs->compute_func = XLALComputeFstatResampCUDA;
  funcs->method_data_destroy_func = XLALDestroyResampMethodData;
  funcs->workspace_destroy_func = XLALDestroyResampWorkspace;

  // Copy the inverse factorial lookup table to GPU memory
  cudaMemcpyToSymbol(lal_fact_inv, (void*)&LAL_FACT_INV, sizeof(REAL8)*LAL_FACT_MAX, 0, cudaMemcpyHostToDevice);

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

  // Move detector time series over to GPU
  MoveMultiCOMPLEX8TimeSeriesHtoD(resamp->multiTimeSeries_DET);

  // header for SRC-frame resampled timeseries buffer
  XLAL_CHECK ( (resamp->multiTimeSeries_SRC_a = (MultiCOMPLEX8TimeSeries *)XLALCalloc ( 1, sizeof(MultiCOMPLEX8TimeSeries)) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( (resamp->multiTimeSeries_SRC_a->data = (COMPLEX8TimeSeries **)XLALCalloc ( numDetectors, sizeof(COMPLEX8TimeSeries) )) != NULL, XLAL_ENOMEM );
  resamp->multiTimeSeries_SRC_a->length = numDetectors;

  XLAL_CHECK ( (resamp->multiTimeSeries_SRC_b = (MultiCOMPLEX8TimeSeries *)XLALCalloc ( 1, sizeof(MultiCOMPLEX8TimeSeries)) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( (resamp->multiTimeSeries_SRC_b->data = (COMPLEX8TimeSeries **)XLALCalloc ( numDetectors, sizeof(COMPLEX8TimeSeries) )) != NULL, XLAL_ENOMEM );
  resamp->multiTimeSeries_SRC_b->length = numDetectors;

  LIGOTimeGPS XLAL_INIT_DECL(epoch0);	// will be set to corresponding SRC-frame epoch when barycentering
  UINT4 numSamplesMax_SRC = 0;
  for ( UINT4 X = 0; X < numDetectors; X ++ )
    {
      // ----- check input consistency ----------
      REAL8 dt_DETX = resamp->multiTimeSeries_DET->data[X]->deltaT;
      XLAL_CHECK ( dt_DET == dt_DETX, XLAL_EINVAL, "Input timeseries must have identical 'deltaT(X=%d)' (%.16g != %.16g)\n", X, dt_DET, dt_DETX);

      REAL8 fHetX = resamp->multiTimeSeries_DET->data[X]->f0;
      XLAL_CHECK ( fabs( fHet - fHetX ) < LAL_REAL8_EPS * fHet, XLAL_EINVAL, "Input timeseries must have identical heterodyning frequency 'f0(X=%d)' (%.16g != %.16g)\n", X, fHet, fHetX );

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

  MoveMultiCOMPLEX8TimeSeriesHtoD(resamp->multiTimeSeries_SRC_b);
  MoveMultiCOMPLEX8TimeSeriesHtoD(resamp->multiTimeSeries_SRC_a);
  XLAL_CHECK ( numSamplesFFT >= numSamplesMax_SRC, XLAL_EFAILED, "[numSamplesFFT = %d] < [numSamplesMax_SRC = %d]\n", numSamplesFFT, numSamplesMax_SRC );

  // ---- re-use shared workspace, or allocate here ----------
  ResampWorkspace *ws = (ResampWorkspace*) common->workspace;
  if ( ws != NULL )
    {
      if ( numSamplesFFT > ws->numSamplesFFTAlloc )
        {
          cudaFree ( ws->FabX_Raw );
          XLAL_CHECK ( (cudaMalloc((void **)&ws->FabX_Raw, numSamplesFFT*sizeof(COMPLEX8)) == cudaSuccess), XLAL_ENOMEM);
          cudaFree ( ws->TS_FFT );
          XLAL_CHECK ( (cudaMalloc((void **)&ws->TS_FFT, numSamplesFFT*sizeof(COMPLEX8)) == cudaSuccess), XLAL_ENOMEM);

          ws->numSamplesFFTAlloc = numSamplesFFT;
        }

      // adjust maximal SRC-frame timeseries length, if necessary - in the CPU version this was realloc'd, but I don't think we can do that here.
      if ( numSamplesMax_SRC > ws->TStmp1_SRC->length ) {
        DestroyCOMPLEX8VectorCUDA(ws->TStmp1_SRC);
        XLAL_CHECK ( (ws->TStmp1_SRC   = CreateCOMPLEX8VectorCUDA ( numSamplesMax_SRC )) != NULL, XLAL_EFUNC );

        DestroyCOMPLEX8VectorCUDA(ws->TStmp2_SRC);
        XLAL_CHECK ( (ws->TStmp2_SRC   = CreateCOMPLEX8VectorCUDA ( numSamplesMax_SRC )) != NULL, XLAL_EFUNC );

        DestroyREAL8VectorCUDA(ws->SRCtimes_DET);
        XLAL_CHECK ( (ws->SRCtimes_DET   = CreateREAL8VectorCUDA ( numSamplesMax_SRC )) != NULL, XLAL_EFUNC );
      }

    } // end: if shared workspace given
  else
    {
      XLAL_CHECK ( (ws = (ResampWorkspace *)XLALCalloc ( 1, sizeof(*ws))) != NULL, XLAL_ENOMEM );
      XLAL_CHECK ( (ws->TStmp1_SRC   = CreateCOMPLEX8VectorCUDA ( numSamplesMax_SRC )) != NULL, XLAL_EFUNC );
      XLAL_CHECK ( (ws->TStmp2_SRC   = CreateCOMPLEX8VectorCUDA ( numSamplesMax_SRC )) != NULL, XLAL_EFUNC );
      XLAL_CHECK ( (ws->SRCtimes_DET = CreateREAL8VectorCUDA ( numSamplesMax_SRC )) != NULL, XLAL_EFUNC );
      XLAL_CHECK ( (cudaMalloc((void **)&ws->FabX_Raw, numSamplesFFT*sizeof(COMPLEX8)) == cudaSuccess), XLAL_ENOMEM);
      XLAL_CHECK ( (cudaMalloc((void **)&ws->TS_FFT, numSamplesFFT*sizeof(COMPLEX8)) == cudaSuccess), XLAL_ENOMEM);
      ws->numSamplesFFTAlloc = numSamplesFFT;

      common->workspace = ws;
    } // end: if we create our own workspace

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

} // XLALSetupFstatResampCUDA()

__global__ void CUDAAddToFaFb(cuComplex *Fa_k, cuComplex *Fb_k, cuComplex *FaX_k, cuComplex *FbX_k, UINT4 numFreqBins)
{
  int k = threadIdx.x + blockDim.x*blockIdx.x;
  if(k >= numFreqBins)
    return;
  Fa_k[k] = cuCaddf(Fa_k[k], FaX_k[k]);
  Fb_k[k] = cuCaddf(Fb_k[k], FbX_k[k]);
}

__global__ void CUDAComputeTwoF(REAL4 *twoF, cuComplex *Fa_k, cuComplex *Fb_k, REAL4 A, REAL4 B, REAL4 C, REAL4 E, REAL4 Dinv, UINT4 numFreqBins)
{
  int k = threadIdx.x + blockDim.x*blockIdx.x;
  if(k >= numFreqBins)
    return;
  cuComplex Fa = Fa_k[k];
  cuComplex Fb = Fb_k[k];
  REAL4 Fa_re = cuCrealf(Fa);
  REAL4 Fa_im = cuCimagf(Fa);
  REAL4 Fb_re = cuCrealf(Fb);
  REAL4 Fb_im = cuCimagf(Fb);

  REAL4 twoF_k = 4;	// default fallback = E[2F] in noise when Dinv == 0 due to ill-conditionness of M_munu
  if ( Dinv > 0 )
    {
      twoF_k = 2.0f * Dinv * (  B * ( SQ(Fa_re) + SQ(Fa_im) )
                              + A * ( SQ(Fb_re) + SQ(Fb_im) )
                              - 2.0 * C * (   Fa_re * Fb_re + Fa_im * Fb_im )
                              - 2.0 * E * ( - Fa_re * Fb_im + Fa_im * Fb_re )           // nonzero only in RAA case where Ed!=0
                              );
    }

  twoF[k] = twoF_k;
}

static int
XLALComputeFstatResampCUDA ( FstatResults* Fstats,
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
  XLAL_CHECK ( XLALBarycentricResampleMultiCOMPLEX8TimeSeriesCUDA ( resamp, &thisPoint, common ) == XLAL_SUCCESS, XLAL_EFUNC );

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
  //TODO: Reallocating existing memory has not been implemented in the CUDA version.
 /*
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
          XLAL_CHECK ( (ws->Fa_k = (COMPLEX8 *)XLALRealloc ( ws->Fa_k, numFreqBins * sizeof(COMPLEX8))) != NULL, XLAL_ENOMEM );
          XLAL_CHECK ( (ws->Fb_k = (COMPLEX8 *)XLALRealloc ( ws->Fb_k, numFreqBins * sizeof(COMPLEX8))) != NULL, XLAL_ENOMEM );
        } // only increase workspace arrays
    }
  */
  XLAL_CHECK((cudaMalloc((void **)&ws->FaX_k, numFreqBins*sizeof(COMPLEX8))) == cudaSuccess, XLAL_ENOMEM);
  XLAL_CHECK((cudaMalloc((void **)&ws->FbX_k, numFreqBins*sizeof(COMPLEX8))) == cudaSuccess, XLAL_ENOMEM);
  XLAL_CHECK((cudaMalloc((void **)&ws->Fa_k, numFreqBins*sizeof(COMPLEX8))) == cudaSuccess, XLAL_ENOMEM);
  XLAL_CHECK((cudaMalloc((void **)&ws->Fb_k, numFreqBins*sizeof(COMPLEX8))) == cudaSuccess, XLAL_ENOMEM);
  /*if ( whatToCompute & FSTATQ_FAFB_PER_DET )
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
          XLAL_CHECK ( (ws->FaX_k = (COMPLEX8 *)XLALRealloc ( ws->FaX_k, numFreqBins * sizeof(COMPLEX8))) != NULL, XLAL_ENOMEM );
          XLAL_CHECK ( (ws->FbX_k = (COMPLEX8 *)XLALRealloc ( ws->FbX_k, numFreqBins * sizeof(COMPLEX8))) != NULL, XLAL_ENOMEM );
        } // only increase workspace arrays
    }*/
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
      const COMPLEX8TimeSeries *TimeSeriesX_SRC_a = multiTimeSeries_SRC_a->data[X];
      const COMPLEX8TimeSeries *TimeSeriesX_SRC_b = multiTimeSeries_SRC_b->data[X];

      // compute {Fa^X(f_k), Fb^X(f_k)}: results returned via workspace ws
      XLAL_CHECK ( XLALComputeFaFb_ResampCUDA ( resamp, ws, thisPoint, common->dFreq, numFreqBins, TimeSeriesX_SRC_a, TimeSeriesX_SRC_b ) == XLAL_SUCCESS, XLAL_EFUNC );

      if ( collectTiming ) {
        tic = XLALGetCPUTime();
      }

      if ( X == 0 )
        { // avoid having to memset this array: for the first detector we *copy* results
          cudaMemcpy(ws->Fa_k, ws->FaX_k, sizeof(cuComplex)*numFreqBins, cudaMemcpyDeviceToDevice);
          cudaMemcpy(ws->Fb_k, ws->FbX_k, sizeof(cuComplex)*numFreqBins, cudaMemcpyDeviceToDevice);
        } // end: if X==0
      else
        { // for subsequent detectors we *add to* them
          CUDAAddToFaFb<<<(numFreqBins + CUDA_BLOCK_SIZE - 1)/CUDA_BLOCK_SIZE, CUDA_BLOCK_SIZE>>>(ws->Fa_k, ws->Fb_k, ws->FaX_k, ws->FbX_k, numFreqBins);
        } // end:if X>0

      if ( whatToCompute & FSTATQ_FAFB_PER_DET )
        {
          cudaMemcpy(Fstats->FaPerDet[X], ws->FaX_k, sizeof(cuComplex)*numFreqBins, cudaMemcpyDeviceToHost);
          cudaMemcpy(Fstats->FbPerDet[X], ws->FbX_k, sizeof(cuComplex)*numFreqBins, cudaMemcpyDeviceToHost);
        }

      if ( collectTiming ) {
        toc = XLALGetCPUTime();
        Tau->SumFabX += (toc-tic);
        tic = toc;
      }

      // ----- if requested: compute per-detector Fstat_X_k
      if ( whatToCompute & FSTATQ_2F_PER_DET )
        {
            const REAL4 Ad = resamp->MmunuX[X].Ad;
            const REAL4 Bd = resamp->MmunuX[X].Bd;
            const REAL4 Cd = resamp->MmunuX[X].Cd;
            const REAL4 Ed = resamp->MmunuX[X].Ed;
            const REAL4 Dd_inv = 1.0f / resamp->MmunuX[X].Dd;
            REAL4 *twoF_gpu;
            cudaMalloc((void **)&twoF_gpu, sizeof(REAL4)*numFreqBins);
            CUDAComputeTwoF<<<(numFreqBins + CUDA_BLOCK_SIZE - 1)/CUDA_BLOCK_SIZE, CUDA_BLOCK_SIZE>>>(twoF_gpu, ws->FaX_k, ws->FbX_k, Ad, Bd, Cd, Ed, Dd_inv, numFreqBins);
            cudaMemcpy(Fstats->twoFPerDet[X], twoF_gpu, sizeof(REAL4)*numFreqBins, cudaMemcpyDeviceToHost);
            cudaFree(twoF_gpu);

        } // end: if compute F_X

      if ( collectTiming ) {
        toc = XLALGetCPUTime();
        Tau->Fab2F += ( toc - tic );
      }

    } // for X < numDetectors

  if ( whatToCompute & FSTATQ_FAFB )
  {
    cudaMemcpy(Fstats->Fa, ws->Fa_k, sizeof(cuComplex)*numFreqBins, cudaMemcpyDeviceToHost);
    cudaMemcpy(Fstats->Fb, ws->Fb_k, sizeof(cuComplex)*numFreqBins, cudaMemcpyDeviceToHost);
  }

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
      REAL4 *twoF_gpu;
      cudaMalloc((void **)&twoF_gpu, sizeof(REAL4)*numFreqBins);
      CUDAComputeTwoF<<<(numFreqBins + CUDA_BLOCK_SIZE - 1)/CUDA_BLOCK_SIZE, CUDA_BLOCK_SIZE>>>(twoF_gpu, ws->Fa_k, ws->Fb_k, Ad, Bd, Cd, Ed, Dd_inv, numFreqBins);
      //Fstats->twoF = twoF_gpu;
      cudaMemcpy(Fstats->twoF, twoF_gpu, sizeof(REAL4)*numFreqBins, cudaMemcpyDeviceToHost);
      cudaFree(twoF_gpu);
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
  cudaFree(ws->Fa_k);
  cudaFree(ws->Fb_k);
  cudaFree(ws->FaX_k);
  cudaFree(ws->FbX_k);
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

__global__ void CUDANormFaFb(cuComplex *Fa_out,
                            cuComplex *Fb_out,
                            REAL8 FreqOut0,
                            REAL8 dFreq,
                            REAL8 dtauX,
                            REAL8 dt_SRC,
                            UINT4 numFreqBins)
{
  int k = threadIdx.x + blockDim.x*blockIdx.x;
  if(k >= numFreqBins)
    return;
  REAL8 f_k = FreqOut0 + k * dFreq;
  REAL8 cycles = - f_k * dtauX;
  REAL8 sinphase, cosphase;
  sincospi(2*cycles, &sinphase, &cosphase);
  cuComplex normX_k = { (float) (dt_SRC*cosphase), (float) (dt_SRC*sinphase) };

  Fa_out[k] = cuCmulf(Fa_out[k], normX_k);
  Fb_out[k] = cuCmulf(Fb_out[k], normX_k);
}

__global__ void CUDAPopulateFaFbFromRaw(cuComplex *out, cuComplex *in, UINT4 numFreqBins, UINT4 offset_bins, UINT4 decimateFFT)
{
  int k = threadIdx.x + blockIdx.x*blockDim.x;
  if(k >= numFreqBins)
    return;
  out[k] = in[offset_bins + k*decimateFFT];
}

static int
XLALComputeFaFb_ResampCUDA ( ResampMethodData *resamp,					//!< [in,out] buffered resampling data and workspace
                             ResampWorkspace *ws,					//!< [in,out] resampling workspace (memory-sharing across segments)
                             const PulsarDopplerParams thisPoint,			//!< [in] Doppler point to compute {FaX,FbX} for
                             REAL8 dFreq,						//!< [in] output frequency resolution
                             UINT4 numFreqBins,						//!< [in] number of output frequency bins
                             const COMPLEX8TimeSeries * __restrict__ TimeSeries_SRC_a,	//!< [in] SRC-frame single-IFO timeseries * a(t)
                             const COMPLEX8TimeSeries * __restrict__ TimeSeries_SRC_b	//!< [in] SRC-frame single-IFO timeseries * b(t)
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
  cudaMemset ( ws->TS_FFT, 0, resamp->numSamplesFFT * sizeof(ws->TS_FFT[0]) );
  // ----- compute FaX_k
  // apply spindown phase-factors, store result in zero-padded timeseries for 'FFT'ing
  XLAL_CHECK ( XLALApplySpindownAndFreqShiftCUDA ( ws->TS_FFT, TimeSeries_SRC_a, &thisPoint, freqShift ) == XLAL_SUCCESS, XLAL_EFUNC );

  if ( collectTiming ) {
    toc = XLALGetCPUTime();
    tiRS->Tau.Spin += ( toc - tic);
    tic = toc;
  }

  // Fourier transform the resampled Fa(t)
  cufftHandle fft_plan;
  cufftPlan1d(&fft_plan, resamp->numSamplesFFT, CUFFT_C2C, 1);
  cudaDeviceSynchronize();
  cufftExecC2C(fft_plan, ws->TS_FFT, ws->FabX_Raw, CUFFT_FORWARD);
  cudaDeviceSynchronize();
  if ( collectTiming ) {
    toc = XLALGetCPUTime();
    tiRS->Tau.FFT += ( toc - tic);
    tic = toc;
  }

  CUDAPopulateFaFbFromRaw<<<(numFreqBins + CUDA_BLOCK_SIZE - 1)/CUDA_BLOCK_SIZE, CUDA_BLOCK_SIZE>>>(ws->FaX_k, ws->FabX_Raw, numFreqBins, offset_bins, resamp->decimateFFT);
  if ( collectTiming ) {
    toc = XLALGetCPUTime();
    tiRS->Tau.Copy += ( toc - tic);
    tic = toc;
  }

  // ----- compute FbX_k
  // apply spindown phase-factors, store result in zero-padded timeseries for 'FFT'ing
  XLAL_CHECK ( XLALApplySpindownAndFreqShiftCUDA ( ws->TS_FFT, TimeSeries_SRC_b, &thisPoint, freqShift ) == XLAL_SUCCESS, XLAL_EFUNC );

  if ( collectTiming ) {
    toc = XLALGetCPUTime();
    tiRS->Tau.Spin += ( toc - tic);
    tic = toc;
  }

  // Fourier transform the resampled Fb(t)
  cufftExecC2C(fft_plan, ws->TS_FFT, ws->FabX_Raw, CUFFT_FORWARD);
  cufftDestroy(fft_plan);

  if ( collectTiming ) {
    toc = XLALGetCPUTime();
    tiRS->Tau.FFT += ( toc - tic);
    tic = toc;
  }
  CUDAPopulateFaFbFromRaw<<<(numFreqBins + CUDA_BLOCK_SIZE - 1)/CUDA_BLOCK_SIZE, CUDA_BLOCK_SIZE>>>(ws->FbX_k, ws->FabX_Raw, numFreqBins, offset_bins, resamp->decimateFFT);

  if ( collectTiming ) {
    toc = XLALGetCPUTime();
    tiRS->Tau.Copy += ( toc - tic);
    tic = toc;
  }

  // ----- normalization factors to be applied to Fa and Fb:
  const REAL8 dtauX = GPSDIFF ( TimeSeries_SRC_a->epoch, thisPoint.refTime );
  CUDANormFaFb<<<(numFreqBins + CUDA_BLOCK_SIZE - 1)/CUDA_BLOCK_SIZE, CUDA_BLOCK_SIZE>>> \
              (ws->FaX_k, ws->FbX_k, FreqOut0, dFreq, dtauX, dt_SRC, numFreqBins);
  if ( collectTiming ) {
    toc = XLALGetCPUTime();
    tiRS->Tau.Norm += ( toc - tic);
    tic = toc;
  }
  return XLAL_SUCCESS;

} // XLALComputeFaFb_ResampCUDA()


__global__ void CUDAApplySpindownAndFreqShift(cuComplex *out,
                                              cuComplex *in,
                                              PulsarDopplerParams doppler,
                                              REAL8 freqShift,
                                              REAL8 dt,
                                              REAL8 Dtau0,
                                              UINT4 s_max,
                                              UINT4 numSamplesIn)

{
  int j = threadIdx.x + blockIdx.x*blockDim.x;
  if(j >= numSamplesIn)
      return;
  REAL8 taup_j = j * dt;
  REAL8 Dtau_alpha_j = Dtau0 + taup_j;
  REAL8 cycles = -freqShift * taup_j;
  REAL8 Dtau_pow_kp1 = Dtau_alpha_j;
  for(UINT4 k = 1; k <=  s_max; k++)
  {
    Dtau_pow_kp1 *= Dtau_alpha_j;
    cycles += - lal_fact_inv[k+1] * doppler.fkdot[k] * Dtau_pow_kp1;
  }
  REAL8 cosphase, sinphase;
  sincospi(2*cycles, &sinphase, &cosphase);
  cuComplex em2piphase = {(float)cosphase, (float)sinphase};
  out[j] = cuCmulf(in[j], em2piphase);
}

static int
XLALApplySpindownAndFreqShiftCUDA ( cuComplex *__restrict__ xOut,                           ///< [out] the spindown-corrected SRC-frame timeseries
                                    const COMPLEX8TimeSeries *__restrict__ xIn,		///< [in] the input SRC-frame timeseries
                                    const PulsarDopplerParams *__restrict__ doppler,	///< [in] containing spindown parameters
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

  cudaDeviceSynchronize();
  CUDAApplySpindownAndFreqShift<<<(numSamplesIn + CUDA_BLOCK_SIZE - 1)/CUDA_BLOCK_SIZE, CUDA_BLOCK_SIZE>>> \
                              (xOut, reinterpret_cast<cuComplex *>(xIn->data->data), *doppler, freqShift, dt, Dtau0, s_max, numSamplesIn);
  cudaDeviceSynchronize();
  return XLAL_SUCCESS;

} // XLALApplySpindownAndFreqShiftCUDA()

__global__ void CUDAApplyHetAMInterp(cuComplex *TS_a, cuComplex *TS_b, cuComplex *TStmp1, cuComplex *TStmp2, UINT4 numSamples)
{
  int j = threadIdx.x + blockIdx.x*blockDim.x;
  if(j >= numSamples)
    return;
  TS_b[j] = cuCmulf(TS_a[j], TStmp2[j]);
  TS_a[j] = cuCmulf(TS_a[j], TStmp1[j]);
}

__global__ void CUDAApplyHetAMSrc(cuComplex *TStmp1_SRC_data,
                                  cuComplex *TStmp2_SRC_data,
                                  REAL8 *ti_DET_data,
                                  UINT4 iStart_SRC_al,
                                  REAL8 tStart_SRC_0,
                                  REAL8 tStart_DET_0,
                                  REAL8 dt_SRC,
                                  REAL8 tMid_DET_al,
                                  REAL8 tMid_SRC_al,
                                  REAL8 Tdot_al,
                                  REAL8 fHet,
                                  REAL4 a_al,
                                  REAL4 b_al,
                                  UINT4 numSamples)
{
  int j = threadIdx.x + blockIdx.x*blockDim.x;
  if(j >= numSamples)
    return;

  UINT4 iSRC_al_j  = iStart_SRC_al + j;

  // for each time sample in the SRC frame, we estimate the corresponding detector time,
  // using a linear approximation expanding around the midpoint of each SFT
  REAL8 t_SRC = tStart_SRC_0 + iSRC_al_j * dt_SRC;
  REAL8 ti_DET_al_j = tMid_DET_al + ( t_SRC - tMid_SRC_al ) / Tdot_al;
  ti_DET_data [ iSRC_al_j ] = ti_DET_al_j;

  // pre-compute correction factors due to non-zero heterodyne frequency of input
  REAL8 tDiff = iSRC_al_j * dt_SRC + (tStart_DET_0 - ti_DET_al_j);      // tSRC_al_j - tDET(tSRC_al_j)
  REAL8 cycles = fmod ( fHet * tDiff, 1.0 );				// the accumulated heterodyne cycles

  REAL8 cosphase, sinphase;
  sincospi(-2*cycles, &sinphase, &cosphase);                                   // the real and imaginary parts of the phase correction

  cuComplex ei2piphase = {(float) cosphase, (float) sinphase};

  // apply AM coefficients a(t), b(t) to SRC frame timeseries [alternate sign to get final FFT return DC in the middle]
  //REAL4 signum = signumLUT [ (iSRC_al_j % 2) ];
  cuComplex signum = {1.0f - 2*(iSRC_al_j % 2), 0}; // alternating sign, avoid branching
  ei2piphase = cuCmulf(ei2piphase, signum);

  TStmp1_SRC_data[iSRC_al_j] = {ei2piphase.x * a_al, ei2piphase.y * a_al};
  TStmp2_SRC_data[iSRC_al_j] = {ei2piphase.x * b_al, ei2piphase.y * b_al};
}

__global__ void CUDASincInterp(cuComplex *out,
                           REAL8 *t_out,
                           cuComplex *in,
                           REAL8 *win,
                           UINT4 Dterms,
                           UINT4 numSamplesIn,
                           UINT4 numSamplesOut,
                           REAL8 tmin,
                           REAL8 dt,
                           REAL8 oodt)
{
  int l = threadIdx.x + blockDim.x*blockIdx.x;
  if(l >= numSamplesOut)
    return;
  REAL8 t = t_out[l] - tmin;		// measure time since start of input timeseries

  // samples outside of input timeseries are returned as 0
  if ( (t < 0) || (t > (numSamplesIn-1)*dt) )	// avoid any extrapolations!
    {
      out[l] = {0, 0};
      return;
    }

  REAL8 t_by_dt = t  * oodt;
  INT8 jstar = lround ( t_by_dt );		// bin closest to 't', guaranteed to be in [0, numSamples-1]

  if ( fabs ( t_by_dt - jstar ) < 2.0e-4 )	// avoid numerical problems near peak
    {
      out[l] = in[jstar];	// known analytic solution for exact bin
      return;
    }

  INT4 jStart0 = jstar - Dterms;
  UINT4 jEnd0 = jstar + Dterms;
  UINT4 jStart = max ( jStart0, 0 );
  UINT4 jEnd   = min ( jEnd0, numSamplesIn - 1 );

  REAL4 delta_jStart = (t_by_dt - jStart);
  REAL8 sin0 = sinpi(delta_jStart);

  REAL4 sin0oopi = sin0 * 1.0/M_PI;

  cuComplex y_l = {0, 0};
  REAL8 delta_j = delta_jStart;
  for ( UINT8 j = jStart; j <= jEnd; j ++ )
    {
      cuComplex Cj = {(float) (win[j - jStart0] * sin0oopi / delta_j), 0};

      y_l = cuCaddf(y_l, cuCmulf(Cj, in[j]));

      sin0oopi = -sin0oopi;		// sin-term flips sign every step
      delta_j --;
    } // for j in [j* - Dterms, ... ,j* + Dterms]

  out[l] = y_l;
}

// Reimplementation of XLALCreateHammingREAL8Window
__global__ void CUDACreateHammingWindow(REAL8 *out, UINT4 length)
{

  int i = threadIdx.x + blockDim.x*blockIdx.x;
  if(i >= (length + 1) / 2)
      return;

  int length_reduced = length - 1;
  double Y = length_reduced > 0 ? (2 * i - length_reduced) / (double) length_reduced : 0;
  out[i] = 0.08 + 0.92*pow(cos(M_PI/2*Y), 2);
  out[length - i - 1] = 0.08 + 0.92*pow(cos(M_PI/2*Y), 2);
}

int
SincInterp (COMPLEX8Vector *y_out,		///< [out] output series of interpolated y-values [must be same size as t_out]
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

  UINT4 winLen = 2 * Dterms + 1;

  const REAL8 oodt = 1.0 / dt;
  REAL8 *win_gpu;
  cudaMalloc((void **)&win_gpu, sizeof(REAL8)*winLen);
  CUDACreateHammingWindow<<<(winLen + CUDA_BLOCK_SIZE - 1)/CUDA_BLOCK_SIZE, CUDA_BLOCK_SIZE>>>(win_gpu, winLen);
  //cudaMemcpy(win_gpu, win->data->data, sizeof(REAL8)*win->data->length, cudaMemcpyHostToDevice);
  CUDASincInterp<<<(numSamplesOut + CUDA_BLOCK_SIZE - 1)/CUDA_BLOCK_SIZE, CUDA_BLOCK_SIZE>>> \
                (reinterpret_cast<cuComplex *>(y_out->data),
                t_out->data,
                reinterpret_cast<cuComplex *>(ts_in->data->data),
                win_gpu,
                Dterms,
                numSamplesIn,
                numSamplesOut,
                tmin,
                dt,
                oodt);
  cudaFree(win_gpu);
  return XLAL_SUCCESS;
}

///
/// Performs barycentric resampling on a multi-detector timeseries, updates resampling buffer with results
///
/// NOTE Buffering: this function does check
/// 1) whether the previously-buffered solution can be completely reused (same sky-position and binary parameters), or
/// 2) if at least sky-dependent quantities can be re-used (antenna-patterns + timings) in case only binary parameters changed
///
static int
XLALBarycentricResampleMultiCOMPLEX8TimeSeriesCUDA ( ResampMethodData *resamp,		// [in/out] resampling input and buffer (to store resampling TS)
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
  REAL8 tic = 0, toc = 0;
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
      REAL8 fHetX = resamp->multiTimeSeries_DET->data[X]->f0;
      XLAL_CHECK ( fabs( fHet - fHetX ) < LAL_REAL8_EPS * fHet, XLAL_EINVAL, "Input timeseries must have identical heterodyning frequency 'f0(X=%d)' (%.16g != %.16g)\n", X, fHet, fHetX );
      REAL8 TsftX = common->multiTimestamps->data[X]->deltaT;
      XLAL_CHECK ( Tsft == TsftX, XLAL_EINVAL, "Input timestamps must have identical stepsize 'Tsft(X=%d)' (%.16g != %.16g)\n", X, Tsft, TsftX );

      TimeSeries_SRCX_a->f0 = fHet;
      TimeSeries_SRCX_b->f0 = fHet;
      // set SRC-frame time-series start-time
      REAL8 tStart_SRC_0 = refTime8 + SRCtimesX->DeltaT->data[0] - (0.5*Tsft) * SRCtimesX->Tdot->data[0];
      LIGOTimeGPS epoch;
      GPSSETREAL8 ( epoch, tStart_SRC_0 );
      TimeSeries_SRCX_a->epoch = epoch;
      TimeSeries_SRCX_b->epoch = epoch;

      // make sure all output samples are initialized to zero first, in case of gaps
      cudaMemset ( TimeSeries_SRCX_a->data->data, 0, TimeSeries_SRCX_a->data->length * sizeof(TimeSeries_SRCX_a->data->data[0]) );
      cudaMemset ( TimeSeries_SRCX_b->data->data, 0, TimeSeries_SRCX_b->data->length * sizeof(TimeSeries_SRCX_b->data->data[0]) );
      // make sure detector-frame timesteps to interpolate to are initialized to 0, in case of gaps
      cudaMemset ( ws->SRCtimes_DET->data, 0, ws->SRCtimes_DET->length * sizeof(ws->SRCtimes_DET->data[0]) );

      cudaMemset ( ws->TStmp1_SRC->data, 0, ws->TStmp1_SRC->length * sizeof(ws->TStmp1_SRC->data[0]) );
      cudaMemset ( ws->TStmp2_SRC->data, 0, ws->TStmp2_SRC->length * sizeof(ws->TStmp2_SRC->data[0]) );

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
          CUDAApplyHetAMSrc<<<(numSamplesSFT_SRC_al + CUDA_BLOCK_SIZE - 1)/CUDA_BLOCK_SIZE, CUDA_BLOCK_SIZE>>> \
                            (
                              reinterpret_cast<cuComplex *>(ws->TStmp1_SRC->data),
                              reinterpret_cast<cuComplex *>(ws->TStmp2_SRC->data),
                              ti_DET->data,
                              iStart_SRC_al,
                              tStart_SRC_0,
                              tStart_DET_0,
                              dt_SRC,
                              tMid_DET_al,
                              tMid_SRC_al,
                              Tdot_al,
                              fHet,
                              a_al,
                              b_al,
                              numSamplesSFT_SRC_al
                            );
        } // for  alpha < numSFTsX

      XLAL_CHECK ( ti_DET->length >= TimeSeries_SRCX_a->data->length, XLAL_EINVAL );
      UINT4 bak_length = ti_DET->length;
      ti_DET->length = TimeSeries_SRCX_a->data->length;

      //XLAL_CHECK ( XLALSincInterpolateCOMPLEX8TimeSeries ( TimeSeries_SRCX_a->data, ti_DET, TimeSeries_DETX, resamp->Dterms ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK ( SincInterp ( TimeSeries_SRCX_a->data, ti_DET, TimeSeries_DETX, resamp->Dterms ) == XLAL_SUCCESS, XLAL_EFUNC );

      ti_DET->length = bak_length;

      // apply heterodyne correction and AM-functions a(t) and b(t) to interpolated timeseries

      CUDAApplyHetAMInterp<<<(numSamples_SRCX + CUDA_BLOCK_SIZE - 1)/CUDA_BLOCK_SIZE, CUDA_BLOCK_SIZE>>> \
                    (reinterpret_cast<cuComplex *>(TimeSeries_SRCX_a->data->data),
                    reinterpret_cast<cuComplex *>(TimeSeries_SRCX_b->data->data),
                    reinterpret_cast<cuComplex *>(ws->TStmp1_SRC->data),
                    reinterpret_cast<cuComplex *>(ws->TStmp2_SRC->data),
                    numSamples_SRCX);
    } // for X < numDetectors

  if ( collectTiming ) {
    toc = XLALGetCPUTime();
    Tau->Bary = (toc-tic);
  }

  return XLAL_SUCCESS;

} // XLALBarycentricResampleMultiCOMPLEX8TimeSeriesCUDA()

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

// @}

// Local Variables:
// mode: c
// End:
