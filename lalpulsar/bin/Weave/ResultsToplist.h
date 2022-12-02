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
// Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
// MA 02110-1301 USA
//

#ifndef _RESULTS_TOPLIST_H
#define _RESULTS_TOPLIST_H

///
/// \file
/// \ingroup lalpulsar_bin_Weave
/// \brief Module which handles toplists of results
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
/// Toplist item
///
struct tagWeaveResultsToplistItem {
  /// Serial number of template
  UINT8 serial;
  /// Index of coherent templates (only needed for per-segment output)
  UINT8 *coh_index;
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
  /// All statistics values computed in this template in 'stage[0]' (first pass) and 'stage[1]' ('recalculation' step without interpolation)
  WeaveStatisticsValues stage[2];
};

///
/// Function which returns pointer to array of statistics by which toplist items are ranked
///
typedef const REAL4 *( *WeaveResultsToplistRankingStats )( const WeaveSemiResults *semi_res );

///
/// Function which returns the value of the statistic by which toplist items are ranked
///
typedef REAL4( *WeaveResultsToplistItemGetRankStat )( const WeaveResultsToplistItem *item );

///
/// Function which sets the value of the statistic by which toplist items are ranked
///
typedef void ( *WeaveResultsToplistItemSetRankStat )( WeaveResultsToplistItem *item, const REAL4 value );

WeaveResultsToplist *XLALWeaveResultsToplistCreate(
  const size_t nspins,
  WeaveStatisticsParams *statistics_params,
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
  const BOOLEAN sort_by_semi_phys,
  const UINT4 round_param_to_dp,
  const UINT4 round_param_to_sf,
  const REAL8 unmatched_item_tol,
  const REAL8 param_tol_mism,
  const VectorComparison *result_tol,
  const UINT4 toplist_compare_limit,
  const WeaveResultsToplist *toplist_1,
  const WeaveResultsToplist *toplist_2
  );

#ifdef __cplusplus
}
#endif

#endif // _RESULTS_TOPLIST_H

// Local Variables:
// c-file-style: "linux"
// c-basic-offset: 2
// End:
