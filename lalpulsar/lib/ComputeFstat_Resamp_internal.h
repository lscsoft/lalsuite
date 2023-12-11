//
// Copyright (C) 2013--2015, 2020 Karl Wette
// Copyright (C) 2015 Reinhard Prix
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
// Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
// MA  02110-1301  USA
//

#include "config.h"

#include <lal/ComputeFstat.h>

// ============================================================================================================ //
//                                                                                                              //
// This file should **ONLY** contain definitions that **MUST** be shared between F-statistic resampling methods //
//                                                                                                              //
// ============================================================================================================ //

// ---------- Shared constants/defines ---------- //

// ---------- Shared macro definitions ---------- //

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

// ---------- Shared global variables ---------- //

// ---------- Shared struct definitions ---------- //

// ---------- BEGIN: Resamp-specific timing model data ----------
typedef struct tagTimings_t {
  REAL4 Total;          // total time spent in XLALComputeFstatResamp()
  REAL4 Bary;           // time spent (in this call) in barycentric resampling
  REAL4 Spin;           // time spent in spindown+frequency correction
  REAL4 FFT;            // time spent in FFT
  REAL4 Copy;           // time spent copying results from FFT to FabX
  REAL4 Norm;           // time spent normalizing the final Fa,Fb
  REAL4 Fab2F;          // time to compute Fstat from {Fa,Fb}
  REAL4 Mem;            // time to realloc and Memset-0 arrays
  REAL4 SumFabX;        // time to sum_X Fab^X
  BOOLEAN BufferRecomputed; // did we need to recompute the buffer this time?
} Timings_t;

// Resamp-specific timing model data
typedef struct tagFstatTimingResamp {
  UINT4 NsampFFT0;      // original number of FFT samples (not rounded to power-of-two)
  UINT4 NsampFFT;       // actual number of FFT samples (rounded up to power-of-two if optArgs->resampFFTPowerOf2 == true)
  REAL4 Resolution;     // (internal) frequency resolution 'R' in natural units: df_internal = R / T_FFT\n

  REAL4 tau0_Fbin;      // timing coefficient for all contributions scaling with output frequency-bins
  REAL4 tau0_spin;      // timing coefficient for spindown-correction
  REAL4 tau0_FFT;       // timing coefficient for FFT-time
  REAL4 tau0_bary;      // timing coefficient for barycentering

  Timings_t Tau;

} FstatTimingResamp;

static const char FstatTimingResampHelp[] =
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

// ---------- Shared internal functions ---------- //

static int
XLALGetFstatTiming_Resamp_intern( const FstatTimingResamp *tiRS, FstatTimingModel *timingModel )
{
  XLAL_CHECK( tiRS != NULL, XLAL_EINVAL );
  XLAL_CHECK( timingModel != NULL, XLAL_EINVAL );

  // return method-specific timing model values
  XLAL_INIT_MEM( ( *timingModel ) );

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

  timingModel->numVariables = i + 1;
  timingModel->help      = FstatTimingResampHelp;

  return XLAL_SUCCESS;
} // XLALGetFstatTiming_Resamp_intern()
