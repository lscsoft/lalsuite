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

///
/// \file
/// \ingroup lalapps_pulsar_Weave
///

#include "Statistics.h"
#include "ComputeResults.h"

#include <lal/LALString.h>

///
/// Struct defining the global 'statistics map' that contains all the defining properties
/// of the supported statistics
///
typedef struct {
  WeaveStatisticType val;		///< bitflag value for this statistic
  const char *const name;		///< internal name of this statistics
  WeaveStatisticType dependencies;      ///< set of *direct* input dependencies of this statistic
  const char *const help;               ///< help string explaining this statistic
} WeaveStatisticMap;

///
/// Sets of toplists, extra statistics and dependencies handled by this code
///
#define ENTRY_NONE              WEAVE_STATISTIC_NONE,           "none", "         ", 0, \
    "No statistic selected"

#define ENTRY_COH2F             WEAVE_STATISTIC_COH2F,          "coh2F", "        ", 0, \
    "Per-segment multi-detector coherent 2F statistic"

#define ENTRY_COH2F_DET         WEAVE_STATISTIC_COH2F_DET,      "coh2F_det", "    ", 0, \
    "Per-segment per-detector coherent 2F statistic"

#define ENTRY_MAX2F             WEAVE_STATISTIC_MAX2F,          "max2F", "        ", WEAVE_STATISTIC_COH2F, \
    "Maximum over segments multi-detector coherent 2F statistic"

#define ENTRY_MAX2F_DET         WEAVE_STATISTIC_MAX2F_DET,      "max2F_det", "    ", WEAVE_STATISTIC_COH2F_DET, \
    "Maximum over segments per-detector coherent 2F statistic"

#define ENTRY_SUM2F             WEAVE_STATISTIC_SUM2F,          "sum2F", "        ", WEAVE_STATISTIC_COH2F, \
    "Sum over segments of multi-detector coherent 2F statistic"

#define ENTRY_SUM2F_DET         WEAVE_STATISTIC_SUM2F_DET,      "sum2F_det", "    ", WEAVE_STATISTIC_COH2F_DET, \
    "Sum over segments of single-detector coherent 2F statistic"

#define ENTRY_MEAN2F            WEAVE_STATISTIC_MEAN2F,         "mean2F", "       ", WEAVE_STATISTIC_SUM2F, \
    "Average over segments of multi-detector coherent 2F statistic"

#define ENTRY_MEAN2F_DET        WEAVE_STATISTIC_MEAN2F_DET,     "mean2F_det", "   ", WEAVE_STATISTIC_SUM2F_DET, \
    "Average over segments of single-detector coherent 2F statistic"

#define ENTRY_BSGL              WEAVE_STATISTIC_BSGL,           "log10BSGL", "    ", WEAVE_STATISTIC_SUM2F|WEAVE_STATISTIC_SUM2F_DET, \
    "Bayes factor 'Signal' vs 'Gaussian noise' or 'Line'"

#define ENTRY_BSGLtL            WEAVE_STATISTIC_BSGLtL,         "log10BSGLtL", "  ", WEAVE_STATISTIC_SUM2F|WEAVE_STATISTIC_SUM2F_DET|WEAVE_STATISTIC_MAX2F_DET, \
    "Bayes factor 'Signal' vs 'Gaussian noise' or 'Line' or 'transient Line'."

#define ENTRY_BtSGLtL           WEAVE_STATISTIC_BtSGLtL,        "log10BtSGLtL", " ", WEAVE_STATISTIC_MAX2F|WEAVE_STATISTIC_SUM2F_DET|WEAVE_STATISTIC_MAX2F_DET, \
    "Bayes factor 'transient Signal' vs 'Gaussian noise' or 'Line' or 'transient Line'."

#define ENTRY_NCOUNT            WEAVE_STATISTIC_NCOUNT,         "ncount", "       ", WEAVE_STATISTIC_COH2F, \
    "Multi-detector 'Hough' number count of 'threshold crossings' heavyside(2F - 2Fth) over segments"

#define ENTRY_NCOUNT_DET        WEAVE_STATISTIC_NCOUNT_DET,     "ncount_det", "   ", WEAVE_STATISTIC_COH2F_DET, \
    "Per-detector 'Hough' number count of 'threshold crossings' heavyside(2F - 2Fth) over segments"

#define ENTRY_2_NAME(X) ENTRY_2_NAME_X(X)
#define ENTRY_2_NAME_X(v,n,s,d,h)  [XLAL_BIT2IDX(v)] = n

#define ENTRY_2_MAP(X) ENTRY_2_MAP_X(X)
#define ENTRY_2_MAP_X(v,n,s,d,h)  { .val = v, .name = n, .dependencies = d, .help = h }

#define ENTRY_2_CHOICES(X) ENTRY_2_CHOICES_X(X)
#define ENTRY_2_CHOICES_X(v,n,s,d,h) { .val = v, .name = n }

#define ENTRY_2_HELPSTR(X) ENTRY_2_HELPSTR_X(X)
#define ENTRY_2_HELPSTR_X(v,n,s,d,h) " - " n s ": " h ".\n"

const char *const WeaveStatisticNamesByIndex[XLAL_BIT2IDX( WEAVE_STATISTIC_MAX )] = {
  ENTRY_2_NAME( ENTRY_COH2F ),
  ENTRY_2_NAME( ENTRY_COH2F_DET ),
  ENTRY_2_NAME( ENTRY_MAX2F ),
  ENTRY_2_NAME( ENTRY_MAX2F_DET ),
  ENTRY_2_NAME( ENTRY_SUM2F ),
  ENTRY_2_NAME( ENTRY_SUM2F_DET ),
  ENTRY_2_NAME( ENTRY_MEAN2F ),
  ENTRY_2_NAME( ENTRY_MEAN2F_DET ),
  ENTRY_2_NAME( ENTRY_BSGL ),
  ENTRY_2_NAME( ENTRY_BSGLtL ),
  ENTRY_2_NAME( ENTRY_BtSGLtL ),
  ENTRY_2_NAME( ENTRY_NCOUNT ),
  ENTRY_2_NAME( ENTRY_NCOUNT_DET ),
};

///
/// Array of descriptor structs for all statistics supported by Weave
///
const WeaveStatisticMap statistic_map[] = {
  ENTRY_2_MAP( ENTRY_COH2F ),
  ENTRY_2_MAP( ENTRY_COH2F_DET ),
  ENTRY_2_MAP( ENTRY_MAX2F ),
  ENTRY_2_MAP( ENTRY_MAX2F_DET ),
  ENTRY_2_MAP( ENTRY_SUM2F ),
  ENTRY_2_MAP( ENTRY_SUM2F_DET ),
  ENTRY_2_MAP( ENTRY_MEAN2F ),
  ENTRY_2_MAP( ENTRY_MEAN2F_DET ),
  ENTRY_2_MAP( ENTRY_BSGL ),
  ENTRY_2_MAP( ENTRY_BSGLtL ),
  ENTRY_2_MAP( ENTRY_BtSGLtL ),
  ENTRY_2_MAP( ENTRY_NCOUNT ),
  ENTRY_2_MAP( ENTRY_NCOUNT_DET ),
};

// Total set of current supported statistics
#define SUPPORTED_STATISTICS ( \
    0 \
    | WEAVE_STATISTIC_COH2F \
    | WEAVE_STATISTIC_COH2F_DET \
    | WEAVE_STATISTIC_MAX2F \
    | WEAVE_STATISTIC_MAX2F_DET \
    | WEAVE_STATISTIC_SUM2F \
    | WEAVE_STATISTIC_SUM2F_DET \
    | WEAVE_STATISTIC_MEAN2F \
    | WEAVE_STATISTIC_MEAN2F_DET \
    | WEAVE_STATISTIC_BSGL \
    | WEAVE_STATISTIC_BSGLtL \
    | WEAVE_STATISTIC_BtSGLtL \
    | WEAVE_STATISTIC_NCOUNT \
    | WEAVE_STATISTIC_NCOUNT_DET \
    )
const UserChoices WeaveStatisticChoices = {
  ENTRY_2_CHOICES( ENTRY_NONE ),
  ENTRY_2_CHOICES( ENTRY_COH2F ),
  ENTRY_2_CHOICES( ENTRY_COH2F_DET ),
  ENTRY_2_CHOICES( ENTRY_MAX2F ),
  ENTRY_2_CHOICES( ENTRY_MAX2F_DET ),
  ENTRY_2_CHOICES( ENTRY_SUM2F ),
  ENTRY_2_CHOICES( ENTRY_SUM2F_DET ),
  ENTRY_2_CHOICES( ENTRY_MEAN2F ),
  ENTRY_2_CHOICES( ENTRY_MEAN2F_DET ),
  ENTRY_2_CHOICES( ENTRY_BSGL ),
  ENTRY_2_CHOICES( ENTRY_BSGLtL ),
  ENTRY_2_CHOICES( ENTRY_BtSGLtL ),
  ENTRY_2_CHOICES( ENTRY_NCOUNT ),
  ENTRY_2_CHOICES( ENTRY_NCOUNT_DET ),
  { SUPPORTED_STATISTICS, "all" }
};
const char *const WeaveStatisticHelpString = \
  ENTRY_2_HELPSTR( ENTRY_COH2F ) \
  ENTRY_2_HELPSTR( ENTRY_COH2F_DET ) \
  ENTRY_2_HELPSTR( ENTRY_MAX2F ) \
  ENTRY_2_HELPSTR( ENTRY_MAX2F_DET ) \
  ENTRY_2_HELPSTR( ENTRY_SUM2F ) \
  ENTRY_2_HELPSTR( ENTRY_SUM2F_DET ) \
  ENTRY_2_HELPSTR( ENTRY_MEAN2F ) \
  ENTRY_2_HELPSTR( ENTRY_MEAN2F_DET ) \
  ENTRY_2_HELPSTR( ENTRY_BSGL ) \
  ENTRY_2_HELPSTR( ENTRY_BSGLtL ) \
  ENTRY_2_HELPSTR( ENTRY_BtSGLtL ) \
  ENTRY_2_HELPSTR( ENTRY_NCOUNT ) \
  ENTRY_2_HELPSTR( ENTRY_NCOUNT_DET ) \
  ;
  
// Subset of statistics that are supported as toplist ranking statistics
#define SUPPORTED_TOPLISTS ( \
    0 \
    | WEAVE_STATISTIC_MEAN2F \
    | WEAVE_STATISTIC_SUM2F \
    | WEAVE_STATISTIC_BSGL \
    | WEAVE_STATISTIC_BSGLtL \
    | WEAVE_STATISTIC_BtSGLtL \
    )
const UserChoices WeaveToplistChoices = {
  ENTRY_2_CHOICES( ENTRY_MEAN2F ),
  ENTRY_2_CHOICES( ENTRY_SUM2F ),
  ENTRY_2_CHOICES( ENTRY_BSGL ),
  ENTRY_2_CHOICES( ENTRY_BSGLtL ),
  ENTRY_2_CHOICES( ENTRY_BtSGLtL ),
  {SUPPORTED_TOPLISTS, "all" }
};
const char *const WeaveToplistHelpString = \
  ENTRY_2_HELPSTR( ENTRY_MEAN2F ) \
  ENTRY_2_HELPSTR( ENTRY_SUM2F ) \
  ENTRY_2_HELPSTR( ENTRY_BSGL ) \
  ENTRY_2_HELPSTR( ENTRY_BSGLtL ) \
  ENTRY_2_HELPSTR( ENTRY_BtSGLtL ) \
  ;

///
/// Set all bits in 'deps' corresponding to *direct* dependencies of the set of input statistics 'stat'
///
  int XLALWeaveStatisticsSetDirectDependencies(
    WeaveStatisticType *deps,
    const WeaveStatisticType stats
    )
{
  XLAL_CHECK( ( stats & ~SUPPORTED_STATISTICS ) == 0, XLAL_EINVAL );

  WeaveStatisticType tmp = 0;
  for ( size_t i=0; i < XLAL_NUM_ELEM( statistic_map ); ++i ) {
    if ( stats & statistic_map[i].val ) {
      tmp |= statistic_map[i].dependencies;
    }
  }

  ( *deps ) |= tmp;

  return XLAL_SUCCESS;

}


///
/// Fill StatisticsParams logic for given toplist and extra-output stats
///
int XLALWeaveStatisticsParamsSetDependencyMap(
  WeaveStatisticsParams *statistics_params,	///< [out] statstics dependency map
  const WeaveStatisticType toplist_stats,	///< [in] requested toplist statistics
  const WeaveStatisticType extra_output_stats,	///< [in] requested 'extra' (stage0) output statistics
  const WeaveStatisticType recalc_stats		///< [in] requested 'recalc' (stage1) statistics
  )
{
  XLAL_CHECK( statistics_params != NULL, XLAL_EFAULT );
  XLAL_CHECK( ( toplist_stats & ( ~SUPPORTED_TOPLISTS ) ) == 0, XLAL_EINVAL );
  XLAL_CHECK( ( extra_output_stats & ( ~SUPPORTED_STATISTICS ) ) == 0, XLAL_EINVAL );
  XLAL_CHECK( ( recalc_stats & ( ~SUPPORTED_STATISTICS ) ) == 0, XLAL_EINVAL );

  WeaveStatisticType stage0_stats_to_output = ( toplist_stats | extra_output_stats );

  // Work out the total set of all statistics we need to compute by
  // expanding the statistics dependencies until converged [tree fully expanded]
  WeaveStatisticType stats_to_compute = stage0_stats_to_output;    // Start value
  WeaveStatisticType mainloop_stats = toplist_stats;      // Start value
  WeaveStatisticType prev_stats_to_compute, prev_mainloop_stats;
  do {
    prev_stats_to_compute = stats_to_compute;
    prev_mainloop_stats = mainloop_stats;

    XLALWeaveStatisticsSetDirectDependencies( &stats_to_compute, stats_to_compute );
    XLALWeaveStatisticsSetDirectDependencies( &mainloop_stats, mainloop_stats );

  } while ( ( prev_stats_to_compute != stats_to_compute ) && ( prev_mainloop_stats != mainloop_stats ) );

  // Special handling of stage0 'coh2F' and 'coh2F_det': these can *only* be computed as "main-loop" statistics!
  // as they are defined to refer to the 'fine grid with (typically) interpolation', while
  if ( stats_to_compute & WEAVE_STATISTIC_COH2F ) {
    mainloop_stats |= WEAVE_STATISTIC_COH2F;
  }
  if ( stats_to_compute & WEAVE_STATISTIC_COH2F_DET ) {
    mainloop_stats |= WEAVE_STATISTIC_COH2F_DET;
  }

  WeaveStatisticType completionloop_stats_stage0 = stats_to_compute & ( ~mainloop_stats );

  // Figure out which mainloop statistics to keep outside of main loop: either
  // 1) because they have been requested for output, or
  // 2) they are a _direct_ (stage0) completionloop dependency,
  // All other mainloop stats can be thrown away safely after the mainloop.
  WeaveStatisticType mainloop_stats_to_keep = 0;

  // 1) if requested for output:
  mainloop_stats_to_keep |= ( mainloop_stats & stage0_stats_to_output );

  // 2) if *direct* (stage0) completionloop dependencies:
  WeaveStatisticType completionloop_stage0_deps = 0;
  XLALWeaveStatisticsSetDirectDependencies( &completionloop_stage0_deps, completionloop_stats_stage0 );
  mainloop_stats_to_keep |= ( mainloop_stats & completionloop_stage0_deps );

  // 3) determine full recalc-stats dependencies
  WeaveStatisticType recalc_stats_deps = recalc_stats;
  WeaveStatisticType prev_recalc_stats_deps;
  do {
    prev_recalc_stats_deps = recalc_stats_deps;
    XLALWeaveStatisticsSetDirectDependencies( &recalc_stats_deps, recalc_stats_deps );
  } while ( prev_recalc_stats_deps != recalc_stats_deps );

  stats_to_compute |= recalc_stats_deps;

  // Store the resulting statistics logic in the 'statistics_params' struct
  statistics_params->toplist_statistics = toplist_stats;
  statistics_params->statistics_to_output[0] = stage0_stats_to_output;
  statistics_params->statistics_to_output[1] = recalc_stats;

  statistics_params->mainloop_statistics = mainloop_stats;
  statistics_params->mainloop_statistics_to_keep = mainloop_stats_to_keep;
  statistics_params->completionloop_statistics[0] = completionloop_stats_stage0;
  statistics_params->completionloop_statistics[1] = recalc_stats_deps;

  statistics_params->all_statistics_to_compute = stats_to_compute;

  // Count number of toplists
  statistics_params->ntoplists = 0;
  for ( int idx = 0; idx < XLAL_BIT2IDX( WEAVE_STATISTIC_MAX ); ++idx ) {
    if ( statistics_params->toplist_statistics & XLAL_IDX2BIT( idx ) ) {
      ++statistics_params->ntoplists;
    }
  }

  return XLAL_SUCCESS;

}

///
/// Destroy a StatisticsParams struct
///
void XLALWeaveStatisticsParamsDestroy(
  WeaveStatisticsParams *statistics_params
  )
{
  if ( statistics_params == NULL ) {
    return;
  }

  XLALDestroyStringVector( statistics_params->detectors );
  XLALDestroyBSGLSetup( statistics_params->BSGL_setup );
  if ( statistics_params->coh_input != NULL ) {
    for ( size_t i = 0; i < statistics_params->nsegments; ++i ) {
      XLALWeaveCohInputDestroy( statistics_params->coh_input[i] );
    }
    XLALFree( statistics_params->coh_input );
  }
  if ( statistics_params->coh_input_recalc != NULL ) {
    for ( size_t i = 0; i < statistics_params->nsegments; ++i ) {
      XLALWeaveCohInputDestroy( statistics_params->coh_input_recalc[i] );
    }
    XLALFree( statistics_params->coh_input_recalc );
  }

  XLALWeaveCohResultsDestroy( statistics_params->coh_res );

  XLALFree( statistics_params );

  return;

}


// Local Variables:
// c-file-style: "linux"
// c-basic-offset: 2
// End:
