//
// Copyright (C) 2016, 2017 Karl Wette
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA 02111-1307 USA
//

#ifndef _COMPUTE_RESULTS_H
#define _COMPUTE_RESULTS_H

///
/// \file
/// \ingroup lalapps_pulsar_Weave
///

#include "Weave.h"

#include <lal/VectorMath.h>
#include <lal/ComputeFstat.h>

#include "Statistics.h"

#ifdef __cplusplus
extern "C" {
#endif

///
/// \name Functions which compute coherent and semicoherent results
///
/// @{

///
/// Input data required for computing coherent results
///
typedef struct tagWeaveCohInput WeaveCohInput;

///
/// Results of a coherent computation on a single segment
///
typedef struct tagWeaveCohResults WeaveCohResults;

///
/// Final results of a semicoherent computation over many segments
///
typedef struct tagWeaveSemiResults {
  /// Bitflag representing search simulation level
  WeaveSimulationLevel simulation_level;
  /// Struct holding all parameters for which statistics to output and compute, when, and how
  const WeaveStatisticsParams *statistics_params;
  /// Number of detectors
  UINT4 ndetectors;
  /// Number of segments
  UINT4 nsegments;
  /// Frequency spacing for semicoherent results
  double dfreq;
  /// Number of frequencies
  UINT4 nfreqs;
  /// Per-segment coherent template parameters of the first frequency bin (optional)
  PulsarDopplerParams *coh_phys;
  /// Per-segment multi-detector F-statistics per frequency (optional)
  const REAL4 **coh2F;
  /// Per-segment per-detector F-statistics per frequency (optional)
  const REAL4 **coh2F_det[PULSAR_MAX_DETECTORS];
  /// Number of coherent results processed thus far
  UINT4 ncoh_res;
  /// Semicoherent template parameters of the first frequency bin
  PulsarDopplerParams semi_phys;
  /// Maximized-over-segments multi-detector F-statistics per frequency
  REAL4VectorAligned *max2F;
  /// Number of max-comparison of multi-detector F-statistics thus far
  UINT4 nmax2F;
  /// Maximized-over-segments per-detector F-statistics per frequency
  REAL4VectorAligned *max2F_det[PULSAR_MAX_DETECTORS];
  /// Number of max-comparison to per-detector F-statistics thus far
  UINT4 nmax2F_det[PULSAR_MAX_DETECTORS];
  /// Summed multi-detector F-statistics per frequency
  REAL4VectorAligned *sum2F;
  /// Number of additions to multi-detector F-statistics thus far
  UINT4 nsum2F;
  /// Summed per-detector F-statistics per frequency
  REAL4VectorAligned *sum2F_det[PULSAR_MAX_DETECTORS];
  /// Number of additions to per-detector F-statistics thus far
  UINT4 nsum2F_det[PULSAR_MAX_DETECTORS];
  /// Mean multi-detector F-statistics per frequency
  REAL4VectorAligned *mean2F;
  /// Mean per-detector F-statistics per frequency
  REAL4VectorAligned *mean2F_det[PULSAR_MAX_DETECTORS];
  /// Line-robust log10(B_S/GL) statistic
  REAL4VectorAligned *log10BSGL;
  /// Transient Line-robust log10(B_S/GLtL) statistic
  REAL4VectorAligned *log10BSGLtL;
  /// Transient Signal Line-robust log10(B_tS/GLtL) statistic
  REAL4VectorAligned *log10BtSGLtL;

} WeaveSemiResults;

WeaveCohInput *XLALWeaveCohInputCreate(
  const WeaveSimulationLevel simulation_level,
  FstatInput *Fstat_input,
  const WeaveStatisticsParams *statistics_params
  );
void XLALWeaveCohInputDestroy(
  WeaveCohInput *coh_input
  );
int XLALWeaveCohResultsCompute(
  WeaveCohResults **coh_res,
  WeaveCohInput *coh_input,
  const PulsarDopplerParams *coh_phys,
  const UINT4 coh_nfreqs
  );
void XLALWeaveCohResultsDestroy(
  WeaveCohResults *coh_res
  );
int XLALWeaveSemiResultsInit(
  WeaveSemiResults **semi_res,
  const WeaveSimulationLevel simulation_level,
  const UINT4 ndetectors,
  const UINT4 nsegments,
  const PulsarDopplerParams *semi_phys,
  const double dfreq,
  const UINT4 semi_nfreqs,
  const WeaveStatisticsParams *statistics_params
  );
int XLALWeaveSemiResultsAdd(
  WeaveSemiResults *semi_res,
  const WeaveCohResults *coh_res,
  const UINT4 coh_offset
  );
int XLALWeaveSemiResultsComputeMain(
  WeaveSemiResults *semi_res
  );
void XLALWeaveSemiResultsDestroy(
  WeaveSemiResults *semi_res
  );

/// @}

#ifdef __cplusplus
}
#endif

#endif // _COMPUTE_RESULTS_H

// Local Variables:
// c-file-style: "linux"
// c-basic-offset: 2
// End:
