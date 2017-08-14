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

#ifndef _RESULTS_TOPLIST_H
#define _RESULTS_TOPLIST_H

///
/// \file
/// \ingroup lalapps_pulsar_Weave
///

#include "Weave.h"
#include "SetupData.h"
#include "ComputeResults.h"
#include "OutputResults.h"

#include <lal/LALHeap.h>
#include <lal/LFTandTSutils.h>

#ifdef __cplusplus
extern "C" {
#endif

///
/// \name Functions which handle toplists of results
///
/// @{

///
/// Toplist of output results
///
typedef struct tagWeaveResultsToplist WeaveResultsToplist;

///
/// Toplist item
///
typedef struct tagWeaveResultsToplistItem {
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
  /// Coherent multi-detector F-statistics (only needed for per-segment output)
  REAL4 *coh2F;
  /// Coherent per-detector F-statistics (only needed for per-detector and per-segment output)
  REAL4 *coh2F_det[PULSAR_MAX_DETECTORS];
  /// Maximized-over-segments multi-detector F-statistic
  REAL4 max2F;
  /// Maximized-over-segments per-detector F-statistic
  REAL4 max2F_det[PULSAR_MAX_DETECTORS];
  /// Summed multi-detector F-statistic
  REAL4 sum2F;
  /// Summed per-detector F-statistic (only needed for per-detector output)
  REAL4 sum2F_det[PULSAR_MAX_DETECTORS];
  /// Mean multi-detector F-statistic
  REAL4 mean2F;
  /// Mean per-detector F-statistic (only needed for per-detector output)
  REAL4 mean2F_det[PULSAR_MAX_DETECTORS];
  /// Line-robust log10(B_S/GL) statistic
  REAL4 log10BSGL;
  /// Line- and transient-line robust log10(B_S/GLtL) statistic
  REAL4 log10BSGLtL;
  /// Transient- signal and line robust log10(B_tS/GLtL) statistic
  REAL4 log10BtSGLtL;
} WeaveResultsToplistItem;

///
/// Function which returns pointer to array of statistics by which toplist items are ranked
///
typedef const REAL4 *( *WeaveResultsToplistRankingStats )( const WeaveSemiResults *semi_res );

///
/// Function which returns the value of the statistic by which toplist items are ranked
///
typedef REAL4 ( *WeaveResultsToplistItemGetRankStat )( const WeaveResultsToplistItem *item );

///
/// Function which sets the value of the statistic by which toplist items are ranked
///
typedef void ( *WeaveResultsToplistItemSetRankStat )( WeaveResultsToplistItem *item, const REAL4 value );

WeaveResultsToplist *XLALWeaveResultsToplistCreate(
  const size_t nspins,
  const WeaveStatisticsParams *statistics_params,
  const char *stat_name,
  const char *stat_desc,
  const UINT4 toplist_limit,
  WeaveResultsToplistRankingStats toplist_rank_stats_fcn,
  WeaveResultsToplistItemGetRankStat toplist_item_get_rank_stat_fcn,
  WeaveResultsToplistItemSetRankStat toplist_item_set_rank_stat_fcn
  );
void XLALWeaveResultsToplistDestroy(
  WeaveResultsToplist *toplist
  );
int XLALWeaveResultsToplistAdd(
  WeaveResultsToplist *toplist,
  const WeaveSemiResults *semi_res,
  const UINT4 semi_nfreqs
  );
int XLALWeaveResultsToplistCompletionLoop(
  WeaveResultsToplist *toplist
  );
int XLALWeaveResultsToplistWrite(
  FITSFile *file,
  const WeaveResultsToplist *toplist
  );
int XLALWeaveResultsToplistReadAppend(
  FITSFile *file,
  WeaveResultsToplist *toplist
  );
int XLALWeaveResultsToplistCompare(
  BOOLEAN *equal,
  const WeaveSetupData *setup,
  const REAL8 param_tol_mism,
  const VectorComparison *result_tol,
  const WeaveResultsToplist *toplist_1,
  const WeaveResultsToplist *toplist_2
  );

/// @}

#ifdef __cplusplus
}
#endif

#endif // _RESULTS_TOPLIST_H

// Local Variables:
// c-file-style: "linux"
// c-basic-offset: 2
// End:
