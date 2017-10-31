//
// Copyright (C) 2017 Reinhard Prix
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

#ifndef _STATISTICS_H
#define _STATISTICS_H

///
/// \file
/// \ingroup lalapps_pulsar_Weave
/// \brief Module which defines statistics that can be computed and their parameters
///

#include <lal/LineRobustStats.h>
#include <lal/StringVector.h>
#include <lal/UserInput.h>

#ifdef __cplusplus
extern "C" {
#endif

///
/// Bitflags representing all possible statistics that can be computed or returned by Weave
///
typedef enum {
  WEAVE_STATISTIC_NONE                                          = 0,
  /// per segment multi-detector F-statistic
  WEAVE_STATISTIC_COH2F                                         = XLAL_IDX2BIT(0),
  /// per segment per-detector F-statistic
  WEAVE_STATISTIC_COH2F_DET                                     = XLAL_IDX2BIT(1),
  // Maximum over segments multi-detector coherent 2F statistic
  WEAVE_STATISTIC_MAX2F                                         = XLAL_IDX2BIT(2),
  // Maximum over segments per-detector coherent 2F statistic
  WEAVE_STATISTIC_MAX2F_DET                                     = XLAL_IDX2BIT(3),
  /// multi-detector sum (over segments) F-statistic
  WEAVE_STATISTIC_SUM2F                                         = XLAL_IDX2BIT(4),
  /// per detector sum F-statistic
  WEAVE_STATISTIC_SUM2F_DET                                     = XLAL_IDX2BIT(5),
  /// multi-detector average (over segments) F-statistic
  WEAVE_STATISTIC_MEAN2F                                        = XLAL_IDX2BIT(6),
  /// per detector average F-statistic
  WEAVE_STATISTIC_MEAN2F_DET                                    = XLAL_IDX2BIT(7),
  /// line-robust log10(B_S/GL) statistic
  WEAVE_STATISTIC_BSGL                                          = XLAL_IDX2BIT(8),
    /// (transient-)line robust log10(B_S/GLtL) statistic
  WEAVE_STATISTIC_BSGLtL                                        = XLAL_IDX2BIT(9),
  /// (transient-)line robust log10(B_tS/GLtL) statistic
  WEAVE_STATISTIC_BtSGLtL                                       = XLAL_IDX2BIT(10),
  /// Hough number count
  WEAVE_STATISTIC_NCOUNT                                        = XLAL_IDX2BIT(11),
  /// Hough number count per detector
  WEAVE_STATISTIC_NCOUNT_DET                                    = XLAL_IDX2BIT(12),
  /// marker +1 of maximal combined valid statistics value
  WEAVE_STATISTIC_MAX                                           = XLAL_IDX2BIT(13)
} WeaveStatisticType;

///
/// Names of all possible statistics
///
#define WEAVE_STATISTIC_NAME(ws) WeaveStatisticNamesByIndex[XLAL_BIT2IDX(ws)]
/// \cond DONT_DOXYGEN
extern const char *const WeaveStatisticNamesByIndex[XLAL_BIT2IDX(WEAVE_STATISTIC_MAX)];
/// \endcond

///
/// User input choices for toplist ranking statistics
///
extern const UserChoices WeaveToplistChoices;

///
/// User input help string for toplist ranking statistics
///
extern const char *const WeaveToplistHelpString;

///
/// User input choices for all supported statistics
///
extern const UserChoices WeaveStatisticChoices;

///
/// User input help string for all supported statistics
///
extern const char *const WeaveStatisticHelpString;

///
/// Struct holding all parameters and status values for computing various statistics
///
typedef struct tagWeaveStatisticsParams {
  /// ---------- elements describing output statistics [read/write from fits files]
  /// list of detector names
  LALStringVector *detectors;
  /// Number of segments
  UINT4 nsegments;

  /// Number of multi-detector 2F summands (should be == number of segments)
  UINT4 nsum2F;
  /// Number of per-detector 2F summands (should be <= number of segments)
  UINT4 nsum2F_det[PULSAR_MAX_DETECTORS];

  /// ---------- statistics dependency map
  /// Bitflag: set of toplist-ranking statistics
  WeaveStatisticType toplist_statistics;
  /// Bitflag: full set of statistics requested for output (toplist + extra-statistics)
  WeaveStatisticType statistics_to_output;

  /// derived from the above: for internal use only [wont read/write these from fits files]
  /// Bitflag: full set of statistics we need compute (toplist + extra + all dependencies)
  WeaveStatisticType statistics_to_compute;
  /// Bitflag: set of "main-loop" statistics that need to be computed on the semi-coherent "fine" grid
  WeaveStatisticType mainloop_statistics;
  /// Bitflag: subset of "main-loop" statistics to keep around after mainloop: either because 1) needed for output, 2) needed for completionloop-stats
  WeaveStatisticType mainloop_statistics_to_keep;
  /// Bitflag: set of "completion-loop" statistics that will be computed only on the final toplist
  WeaveStatisticType completionloop_statistics;

  /// ---------- input parameters for various statistics
  /// setup for line-robust B_*S/GL* family of statistics
  BSGLSetup *BSGL_setup;

  /// per-segment 2F threshold for computing 'Hough' number counts
  REAL4 nc_2Fth;

} WeaveStatisticsParams;


int XLALWeaveStatisticsSetDirectDependencies(
  WeaveStatisticType *deps,
  const WeaveStatisticType stats
  );

int XLALWeaveStatisticsParamsSetDependencyMap(
  WeaveStatisticsParams *statistics_params,
  const WeaveStatisticType toplist_stats,
  const WeaveStatisticType extra_output_stats
  );

void XLALWeaveStatisticsParamsDestroy (
  WeaveStatisticsParams *statistics_params
  );

#ifdef __cplusplus
}
#endif

#endif // _STATISTICS_H

// Local Variables:
// c-file-style: "linux"
// c-basic-offset: 2
// End:
