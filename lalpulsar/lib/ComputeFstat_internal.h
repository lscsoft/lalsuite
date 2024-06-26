//
// Copyright (C) 2013--2015 Karl Wette
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

// ================================================================================================= //
//                                                                                                   //
// This file should **ONLY** contain definitions that **MUST** be shared between F-statistic methods //
//                                                                                                   //
// ================================================================================================= //

// ---------- Shared constants/defines ---------- //

// ---------- Shared macro definitions ---------- //

#define SQ(x) ( (x) * (x) )

#define MYMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MYMIN(x,y) ( (x) < (y) ? (x) : (y) )

#ifdef __GNUC__
#define likely(x)      __builtin_expect(!!(x), 1)
#define unlikely(x)    __builtin_expect(!!(x), 0)
#else
#define likely(x)      (x)
#define unlikely(x)    (x)
#endif

// ---------- Shared global variables ---------- //

// ---------- Shared struct definitions ---------- //

// Common input data for F-statistic methods
typedef struct {
  LIGOTimeGPS midTime;                                  // Mid-time of SFT data
  REAL8 dFreq;                                          // Requested spacing of \f$ \mathcal{F} \f$ -statistic frequency bins.
  MultiLALDetector detectors;                           // List of detectors
  MultiLIGOTimeGPSVector *multiTimestamps;              // Multi-detector list of SFT timestamps
  MultiNoiseWeights *multiNoiseWeights;                 // Multi-detector noise weights: stored unnormalized (divide by Sinv_Tsft to get normalized version)
  MultiDetectorStateSeries *multiDetectorStates;        // Multi-detector state series
  const EphemerisData *ephemerides;                     // Ephemerides for the time-span of the SFTs
  SSBprecision SSBprec;                                 // Barycentric transformation precision
  void *workspace;                                      // F-statistic method workspace
  BOOLEAN isTimeslice;                                  // Flag if this is a timeslice of another FstatInput struct
  REAL8 allowedMismatchFromSFTLength;                   // optional override for XLALFstatCheckSFTLengthMismatch()
} FstatCommon;

// Pointers to function pointers which perform method-specific operations
typedef struct {
  int ( *compute_func )(                                // F-statistic method computation function
    FstatResults *, const FstatCommon *, void *
  );
  void ( *method_data_destroy_func )( void * );         // F-statistic method data destructor function
  void ( *workspace_destroy_func )( void * );           // Workspace destructor function
} FstatMethodFuncs;

// ---------- Shared internal functions ---------- //

static inline REAL4
compute_fstat_from_fa_fb( COMPLEX8 Fa, COMPLEX8 Fb, REAL4 A, REAL4 B, REAL4 C, REAL4 E, REAL4 Dinv )
{
  REAL4 Fa_re = creal( Fa );
  REAL4 Fa_im = cimag( Fa );
  REAL4 Fb_re = creal( Fb );
  REAL4 Fb_im = cimag( Fb );

  REAL4 twoF = 4;       // default fallback = E[2F] in noise when Dinv == 0 due to ill-conditionness of M_munu
  if ( likely( Dinv > 0 ) ) {
    twoF = 2.0f * Dinv * ( B * ( SQ( Fa_re ) + SQ( Fa_im ) )
                           + A * ( SQ( Fb_re ) + SQ( Fb_im ) )
                           - 2.0 * C * ( Fa_re * Fb_re + Fa_im * Fb_im )
                           - 2.0 * E * ( - Fa_re * Fb_im + Fa_im * Fb_re )           // nonzero only in RAA case where Ed!=0
                         );
  }

  return twoF;

} // compute_fstat_from_fa_fb()
