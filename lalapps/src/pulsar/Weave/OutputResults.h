//
// Copyright (C) 2016 Karl Wette
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

#ifndef _OUTPUT_RESULTS_H
#define _OUTPUT_RESULTS_H

///
/// \file
/// \ingroup lalapps_pulsar_Weave
///

#include "Weave.h"
#include "SetupData.h"
#include "ComputeResults.h"

#include <lal/UserInputParse.h>
#include <lal/LFTandTSutils.h>

#ifdef __cplusplus
extern "C" {
#endif

///
/// \name Routines which handle the output results
///
/// @{

///
/// Bitflags representing possible output toplist types
///
typedef enum {
  /// Toplist ranked by mean multi-detector F-statistic
  WEAVE_TOPLIST_RANKED_MEAN2F           = 0001,
  /// Must be exactly one more than maximum possible bitflag value
  WEAVE_TOPLIST_MAX                     = 0002,
} WeaveToplistType;

///
/// Static array of all #WeaveToplistType choices, for use by the UserInput module parsing routines
///
extern const UserChoices WeaveToplistTypeChoices;

///
/// Output results from a search
///
typedef struct tagWeaveOutputResults WeaveOutputResults;

///
/// Various varameters required to output results
///
typedef struct tagWeaveOutputParams {
  /// Reference time at which search is conducted
  LIGOTimeGPS ref_time;
  /// Number of spindown parameters to output
  size_t nspins;
  /// If outputting per-detector quantities, list of detectors
  LALStringVector *per_detectors;
  /// Number of per-segment items being output (may be zero)
  UINT4 per_nsegments;
} WeaveOutputParams;

///
/// Output result item
///
typedef struct tagWeaveOutputResultItem {
  /// Physical right ascension of semicoherent template
  REAL8 semi_alpha;
  /// Physical right ascension of coherent templates (only needed for per-segment output)
  REAL8 *coh_alpha;
  /// Physical declination of semicoherent template
  REAL8 semi_delta;
  /// Physical declination of coherent templates (only needed for per-segment output)
  REAL8 *coh_delta;
  /// Physical frequency and spindowns of semicoherent template
  REAL8 semi_fkdot[PULSAR_MAX_SPINS];
  /// Physical frequency and spindowns of coherent templates (only needed for per-segment output)
  REAL8 *coh_fkdot[PULSAR_MAX_SPINS];
  /// Mean multi-detector F-statistic
  REAL4 mean2F;
  /// Coherent multi-detector F-statistics (only needed for per-segment output)
  REAL4 *coh2F;
  /// Mean per-detector F-statistic (only needed for per-detector output)
  REAL4 mean2F_det[PULSAR_MAX_DETECTORS];
  /// Coherent per-detector F-statistics (only needed for per-detector and per-segment output)
  REAL4 *coh2F_det[PULSAR_MAX_DETECTORS];
} WeaveOutputResultItem;

WeaveOutputResultItem *XLALWeaveOutputResultItemCreate(
  const WeaveOutputParams *par
  );
int XLALWeaveFillOutputResultItem(
  WeaveOutputResultItem **item,
  BOOLEAN *full_init,
  const WeaveOutputParams *par,
  const WeaveSemiResults *semi_res,
  const size_t nspins,
  const size_t freq_idx
  );
void XLALWeaveOutputResultItemDestroy(
  WeaveOutputResultItem *item
  );
WeaveOutputResults *XLALWeaveOutputResultsCreate(
  const LIGOTimeGPS *ref_time,
  const size_t nspins,
  const LALStringVector *per_detectors,
  const UINT4 per_nsegments,
  const WeaveToplistType toplist_types,
  const int toplist_limit
  );
void XLALWeaveOutputResultsDestroy(
  WeaveOutputResults *out
  );
int XLALWeaveOutputResultsAdd(
  WeaveOutputResults *out,
  const WeaveSemiResults *semi_res,
  const UINT4 semi_nfreqs
  );
int XLALWeaveOutputResultsWrite(
  FITSFile *file,
  const WeaveOutputResults *out
  );
int XLALWeaveOutputResultsReadAppend(
  FITSFile *file,
  WeaveOutputResults **out
  );
int XLALWeaveOutputResultsCompare(
  BOOLEAN *equal,
  const WeaveSetupData *setup,
  const REAL8 param_tol_mism,
  const VectorComparison *result_tol,
  const WeaveOutputResults *out_1,
  const WeaveOutputResults *out_2
  );
int XLALWeaveOutputMiscPerSegInfoWrite(
  FITSFile *file,
  const WeaveSetupData *setup,
  const BOOLEAN write_SFT_info,
  const UINT4 nsegments,
  const WeaveOutputMiscPerSegInfo *per_seg_info
  );

/// @}

#ifdef __cplusplus
}
#endif

#endif // _OUTPUT_RESULTS_H
