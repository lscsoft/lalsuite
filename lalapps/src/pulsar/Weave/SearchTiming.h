//
// Copyright (C) 2017 Karl Wette
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

#ifndef _SEARCH_TIMING_H
#define _SEARCH_TIMING_H

///
/// \file
/// \ingroup lalapps_pulsar_Weave
/// \brief Module which collects search timings and builds a timing model
///

#include "Weave.h"
#include "Statistics.h"
#include "CacheResults.h"

#ifdef __cplusplus
extern "C" {
#endif

///
/// Search sections which are timed individually
///
enum tagWeaveSearchTimingSection {
  /// Parameter space iteration section
  WEAVE_SEARCH_TIMING_ITER,
  /// Cache queries section
  WEAVE_SEARCH_TIMING_QUERY,
  /// Computation of coherent results section
  WEAVE_SEARCH_TIMING_COH,
  /// Computation of per-segment semicoherent results section
  WEAVE_SEARCH_TIMING_SEMISEG,
  /// Computation of semicoherent results section
  WEAVE_SEARCH_TIMING_SEMI,
  /// Result output section
  WEAVE_SEARCH_TIMING_OUTPUT,
  /// Checkpointing section
  WEAVE_SEARCH_TIMING_CKPT,
  /// Unaccounted section
  WEAVE_SEARCH_TIMING_OTHER,
  WEAVE_SEARCH_TIMING_MAX
};

WeaveSearchTiming *XLALWeaveSearchTimingCreate(
  const BOOLEAN detailed_timing,
  const WeaveStatisticsParams *statistics_params
  );
void XLALWeaveSearchTimingDestroy(
  WeaveSearchTiming *tim
  );
int XLALWeaveSearchTimingStart(
  WeaveSearchTiming *tim
  );
int XLALWeaveSearchTimingElapsed(
  WeaveSearchTiming *tim,
  double *wall_elapsed,
  double *cpu_elapsed
  );
int XLALWeaveSearchTimingStop(
  WeaveSearchTiming *tim,
  double *wall_total,
  double *cpu_total
  );
int XLALWeaveSearchTimingSection(
  WeaveSearchTiming *tim,
  const WeaveSearchTimingSection prev_section,
  const WeaveSearchTimingSection next_section
  );
int XLALWeaveSearchTimingStatistic(
  WeaveSearchTiming *tim,
  const WeaveStatisticType prev_statistic,
  const WeaveStatisticType next_statistic
  );
int XLALWeaveSearchTimingWriteInfo(
  FITSFile *file,
  const WeaveSearchTiming *tim,
  const WeaveCacheQueries *queries
  );

#ifdef __cplusplus
}
#endif

#endif // _SEARCH_TIMING_H

// Local Variables:
// c-file-style: "linux"
// c-basic-offset: 2
// End:
