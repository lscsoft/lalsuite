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

///
/// \file
/// \ingroup lalapps_pulsar_Weave
///

#include "SearchTiming.h"

#include <lal/LogPrintf.h>
#include <lal/UserInputParse.h>

///
/// Container for timings and other information for a timing model
///
struct tagWeaveSearchTiming {
  /// Whether to record detailed timing information
  BOOLEAN detailed_timing;
  /// Struct holding all parameters for which statistics to output and compute, when, and how
  const WeaveStatisticsParams *statistics_params;
  /// Number of output results toplists
  size_t ntoplists;
  /// Total wall time
  double wall_total;
  /// Total CPU time
  double cpu_total;
  /// CPU time taken by various sections
  double section_cpu_times[WEAVE_SEARCH_TIMING_MAX];
  /// Current section being timed
  WeaveSearchTimingSection curr_section;
  /// CPU time for current section being timed
  double curr_section_cpu_time;
  /// CPU time taken by various coherent/semicoherent statistics
  double statistic_cpu_times[XLAL_BIT2IDX(WEAVE_STATISTIC_MAX)];
  /// Section in which each statistic is timed
  WeaveSearchTimingSection statistic_section[XLAL_BIT2IDX(WEAVE_STATISTIC_MAX)];
  /// Current statistic being timed
  WeaveStatisticType curr_statistic;
  /// CPU time for current statistic being timed
  double curr_statistic_cpu_time;
};

///
/// Denominator to use for timing constants
///
typedef enum {
  WEAVE_SEARCH_DENOM_NONE,
  /// Per number of computed coherent results
  WEAVE_SEARCH_DENOM_PCOH,
  /// Per number of semicoherent templates
  WEAVE_SEARCH_DENOM_PSEMI,
  /// Per number of semicoherent templates, per number of segments (minus 1)
  WEAVE_SEARCH_DENOM_PSSM1,
  /// Per number of semicoherent templates, per number of segments
  WEAVE_SEARCH_DENOM_PSSEG,
  /// Per semicoherent templates, per number of toplists
  WEAVE_SEARCH_DENOM_PSTOP,
  WEAVE_SEARCH_DENOM_MAX,
} WeaveSearchTimingDenominator;

///
/// Names of denominator to use for timing constants
///
const char *denom_names[WEAVE_SEARCH_DENOM_MAX] = {
  [WEAVE_SEARCH_DENOM_PCOH]  = "pcoh",
  [WEAVE_SEARCH_DENOM_PSEMI] = "psemi",
  [WEAVE_SEARCH_DENOM_PSSM1] = "psemi_psegm1",
  [WEAVE_SEARCH_DENOM_PSSEG] = "psemi_pseg",
  [WEAVE_SEARCH_DENOM_PSTOP] = "psemi_ptopl",
};

///
/// Parameters of search timing sections
///
const struct {
  const char* name;
  const char* comment;
  const WeaveSearchTimingDenominator denom;
} cpu_sections[WEAVE_SEARCH_TIMING_MAX] = {
  [WEAVE_SEARCH_TIMING_ITER]    = {"iter",      "parameter space iteration",                    WEAVE_SEARCH_DENOM_PSEMI},
  [WEAVE_SEARCH_TIMING_QUERY]   = {"query",     "cache queries",                                WEAVE_SEARCH_DENOM_PSSEG},
  [WEAVE_SEARCH_TIMING_COH]     = {"coh",       "computing coherent results",                   WEAVE_SEARCH_DENOM_PCOH},
  [WEAVE_SEARCH_TIMING_SEMISEG] = {"semiseg",   "computing per-segment semicoherent results",   WEAVE_SEARCH_DENOM_PSSM1},
  [WEAVE_SEARCH_TIMING_SEMI]    = {"semi",      "computing semicoherent results",               WEAVE_SEARCH_DENOM_PSEMI},
  [WEAVE_SEARCH_TIMING_OUTPUT]  = {"output",    "result output",                                WEAVE_SEARCH_DENOM_PSTOP},
  [WEAVE_SEARCH_TIMING_CKPT]    = {"ckpt",      "checkpointing",                                WEAVE_SEARCH_DENOM_PSEMI},
  [WEAVE_SEARCH_TIMING_CMPL]    = {"cmpl",      "computing completion-loop results",            WEAVE_SEARCH_DENOM_NONE},
  [WEAVE_SEARCH_TIMING_OTHER]   = {"other",     "unaccounted",                                  WEAVE_SEARCH_DENOM_NONE},
};

///
/// \name Internal functions
///
/// @{

static inline double wall_time(void);
static inline double cpu_time(void);

/// @}

///
/// Return wall time in seconds
///
double wall_time(
  void
  )
{
  return XLALGetTimeOfDay();
}

///
/// Return CPU time in seconds
double cpu_time(
  void
  )
{
  return XLALGetCPUTime();
}

///
/// Create a search timing structure
///
WeaveSearchTiming *XLALWeaveSearchTimingCreate(
  const BOOLEAN detailed_timing,
  const WeaveStatisticsParams *statistics_params
  )
{

  // Allocate memory
  WeaveSearchTiming *tim = XLALCalloc( 1, sizeof( *tim ) );
  XLAL_CHECK_NULL( tim != NULL, XLAL_ENOMEM );

  // Set fields
  tim->detailed_timing = detailed_timing;
  tim->statistics_params = statistics_params;
  tim->curr_section = WEAVE_SEARCH_TIMING_MAX;
  tim->curr_statistic = WEAVE_STATISTIC_NONE;

  return tim;

}

///
/// Destroy a search timing structure
///
void XLALWeaveSearchTimingDestroy(
  WeaveSearchTiming *tim
  )
{
  if ( tim ) {
    XLALFree( tim );
  }
}

///
/// Start timing of search
///
int XLALWeaveSearchTimingStart(
  WeaveSearchTiming *tim
  )
{

  // Check input
  XLAL_CHECK( tim != NULL, XLAL_EFAULT );
  XLAL_CHECK( tim->curr_section == WEAVE_SEARCH_TIMING_MAX, XLAL_EINVAL );

  // Get current wall time
  const double wall_now = wall_time();

  // Get current CPU time
  const double cpu_now = cpu_time();

  // Start timing
  tim->wall_total = wall_now;
  tim->cpu_total = cpu_now;
  for ( size_t i = 0; i < WEAVE_SEARCH_TIMING_MAX; ++i ) {
    tim->section_cpu_times[i] = 0;
  }
  for ( size_t i = 0; i < XLAL_BIT2IDX(WEAVE_STATISTIC_MAX); ++i ) {
    tim->statistic_cpu_times[i] = 0;
    tim->statistic_section[i] = WEAVE_SEARCH_TIMING_MAX;
  }

  // Start timing next section
  tim->curr_section = WEAVE_SEARCH_TIMING_OTHER;
  tim->curr_section_cpu_time = cpu_now;

  return XLAL_SUCCESS;

}

///
/// Return elapsed wall and CPU times since start of search timing
///
int XLALWeaveSearchTimingElapsed(
  WeaveSearchTiming *tim,
  double *wall_elapsed,
  double *cpu_elapsed
  )
{

  // Check input
  XLAL_CHECK( tim != NULL, XLAL_EFAULT );
  XLAL_CHECK( tim->curr_section < WEAVE_SEARCH_TIMING_MAX, XLAL_EINVAL );
  XLAL_CHECK( wall_elapsed != NULL, XLAL_EFAULT );
  XLAL_CHECK( cpu_elapsed != NULL, XLAL_EFAULT );

  // Get current wall time
  const double wall_now = wall_time();

  // Get current CPU time
  const double cpu_now = cpu_time();

  // Return elapsed wall and CPU times
  *wall_elapsed = wall_now - tim->wall_total;
  *cpu_elapsed = cpu_now - tim->cpu_total;

  return XLAL_SUCCESS;

}

///
/// Stop timing of search
///
int XLALWeaveSearchTimingStop(
  WeaveSearchTiming *tim,
  double *wall_total,
  double *cpu_total
  )
{

  // Check input
  XLAL_CHECK( tim != NULL, XLAL_EFAULT );
  XLAL_CHECK( tim->curr_section == WEAVE_SEARCH_TIMING_OTHER, XLAL_EINVAL );
  XLAL_CHECK( wall_total != NULL, XLAL_EFAULT );
  XLAL_CHECK( cpu_total != NULL, XLAL_EFAULT );

  // Get current wall time
  const double wall_now = wall_time();

  // Get current CPU time
  const double cpu_now = cpu_time();

  // Stop timing previous section
  tim->section_cpu_times[tim->curr_section] += cpu_now - tim->curr_section_cpu_time;

  // Stop timing
  tim->wall_total = wall_now - tim->wall_total;
  tim->cpu_total = cpu_now - tim->cpu_total;
  tim->curr_section = WEAVE_SEARCH_TIMING_MAX;

  // Compute remaining CPU time from total CPU time and other timed sections
  tim->section_cpu_times[WEAVE_SEARCH_TIMING_OTHER] = tim->cpu_total;
  for ( int i = 0; i < WEAVE_SEARCH_TIMING_OTHER; ++i ) {
    tim->section_cpu_times[WEAVE_SEARCH_TIMING_OTHER] -= tim->section_cpu_times[i];
  }

  // Return total wall and CPU times
  *wall_total = tim->wall_total;
  *cpu_total = tim->cpu_total;

  return XLAL_SUCCESS;

}

///
/// Change the search section currently being timed
///
int XLALWeaveSearchTimingSection(
  WeaveSearchTiming *tim,
  const WeaveSearchTimingSection prev_section,
  const WeaveSearchTimingSection next_section
  )
{

  // Check input
  XLAL_CHECK( tim != NULL, XLAL_EFAULT );
  XLAL_CHECK( prev_section < WEAVE_SEARCH_TIMING_MAX, XLAL_EINVAL );
  XLAL_CHECK( next_section < WEAVE_SEARCH_TIMING_MAX, XLAL_EINVAL );
  XLAL_CHECK( tim->curr_section == prev_section, XLAL_EINVAL );

  // Get current CPU time
  const double cpu_now = cpu_time();

  // Stop timing previous section
  if ( prev_section != WEAVE_SEARCH_TIMING_OTHER ) {
    tim->section_cpu_times[prev_section] += cpu_now - tim->curr_section_cpu_time;
  }

  // Change section
  tim->curr_section = next_section;

  // Start timing next section
  if ( next_section != WEAVE_SEARCH_TIMING_OTHER ) {
    tim->curr_section_cpu_time = cpu_now;
  }

  return XLAL_SUCCESS;

}

///
/// Change the search statistic currently being timed
///
int XLALWeaveSearchTimingStatistic(
  WeaveSearchTiming *tim,
  const WeaveStatisticType prev_statistic,
  const WeaveStatisticType next_statistic
  )
{

  // Check input
  XLAL_CHECK( tim != NULL, XLAL_EFAULT );
  XLAL_CHECK( prev_statistic < WEAVE_STATISTIC_MAX, XLAL_EINVAL );
  XLAL_CHECK( next_statistic < WEAVE_STATISTIC_MAX, XLAL_EINVAL );
  XLAL_CHECK( tim->curr_statistic == prev_statistic, XLAL_EINVAL );

  // Get current CPU time
  const double cpu_now = cpu_time();

  // Stop timing previous statistic
  if ( prev_statistic & tim->statistics_params->all_statistics_to_compute ) {
    XLAL_CHECK( tim->statistic_section[XLAL_BIT2IDX(prev_statistic)] == tim->curr_section, XLAL_EINVAL );
    tim->statistic_cpu_times[XLAL_BIT2IDX(tim->curr_statistic)] += cpu_now - tim->curr_statistic_cpu_time;
  }

  // Change statistic
  tim->curr_statistic = next_statistic;

  // Start timing next statistic
  if ( next_statistic & tim->statistics_params->all_statistics_to_compute ) {
    tim->statistic_section[XLAL_BIT2IDX(next_statistic)] = tim->curr_section;
    tim->curr_statistic_cpu_time = cpu_now;
  }

  return XLAL_SUCCESS;

}

///
/// Write information from search timing to a FITS file
///
int XLALWeaveSearchTimingWriteInfo(
  FITSFile *file,
  const WeaveSearchTiming *tim,
  const WeaveCacheQueries *queries
  )
{

  // Check input
  XLAL_CHECK( file != NULL, XLAL_EFAULT );
  XLAL_CHECK( tim != NULL, XLAL_EFAULT );
  XLAL_CHECK( queries != NULL, XLAL_EFAULT );

  // Write total wall and CPU time
  XLAL_CHECK( XLALFITSHeaderWriteREAL8( file, "wall total", tim->wall_total, "total wall time" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALFITSHeaderWriteREAL8( file, "cpu total", tim->cpu_total, "total CPU time" ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Return if detailed timing is disabled
  if ( !tim->detailed_timing ) {
    return XLAL_SUCCESS;
  }

  // Get number of computed coherent results, and number of coherent and semicoherent templates
  UINT8 coh_nres = 0, coh_ntmpl = 0, semi_ntmpl = 0;
  XLAL_CHECK( XLALWeaveCacheQueriesGetCounts( queries, &coh_nres, &coh_ntmpl, &semi_ntmpl ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Compute denominators for timing constants
  const UINT8 denoms[WEAVE_SEARCH_DENOM_MAX] = {
    [WEAVE_SEARCH_DENOM_PCOH]   = coh_nres,
    [WEAVE_SEARCH_DENOM_PSEMI]  = semi_ntmpl,
    [WEAVE_SEARCH_DENOM_PSSM1]  = semi_ntmpl * ( tim->statistics_params->nsegments - 1 ),
    [WEAVE_SEARCH_DENOM_PSSEG]  = semi_ntmpl * tim->statistics_params->nsegments,
    [WEAVE_SEARCH_DENOM_PSTOP]  = semi_ntmpl * tim->statistics_params->ntoplists,
  };

  for ( WeaveSearchTimingSection i = 0; i < WEAVE_SEARCH_TIMING_MAX; ++i ) {

    // Get section CPU time and denominator
    const double section_cpu_time = tim->section_cpu_times[i];
    const UINT8 section_denom = denoms[cpu_sections[i].denom];
    const char *section_denom_name = denom_names[cpu_sections[i].denom];

    // Write CPU time taken by timing section
    {
      char keyword[64];
      char comment[1024];
      snprintf( keyword, sizeof( keyword ), "cpu %s", cpu_sections[i].name );
      snprintf( comment, sizeof( comment ), "%s CPU time", cpu_sections[i].comment );
      XLAL_CHECK( XLALFITSHeaderWriteREAL8( file, keyword, section_cpu_time, comment ) == XLAL_SUCCESS, XLAL_EFUNC );
    }

    // Find any statistics computed within timing section
    BOOLEAN statistics_computed = 0;
    double section_cpu_time_remain = section_cpu_time;
    for ( int j = 0; j < XLAL_BIT2IDX(WEAVE_STATISTIC_MAX); ++j ) {
      const double statistic_cpu_time = tim->statistic_cpu_times[j];
      if ( tim->statistic_section[j] == i && statistic_cpu_time > 0 ) {
        statistics_computed = 1;

        // Get statistics name and determine statistic-specific denominator
        UINT8 statistic_denom = 1;
        const char *statistic_name = WEAVE_STATISTIC_NAME(XLAL_IDX2BIT(j));
        const size_t statistic_name_len = strlen( statistic_name );
        if ( statistic_name_len > 4 && strcmp( statistic_name + statistic_name_len - 4, "_det") == 0 ) {

          // Normalise per-detector statistics by number of detectors
          statistic_denom = tim->statistics_params->detectors->length;

        }

        // Write CPU time taken to compute statistic
        {
          char keyword[64];
          char comment[1024];
          snprintf( keyword, sizeof( keyword ), "cpu %s %s", cpu_sections[i].name, statistic_name );
          snprintf( comment, sizeof( comment ), "%s CPU time", statistic_name );
          XLAL_CHECK( XLALFITSHeaderWriteREAL8( file, keyword, statistic_cpu_time, comment ) == XLAL_SUCCESS, XLAL_EFUNC );
        }

        // Write timing constant for computation of statistic
        if ( section_denom > 0 ) {
          XLAL_CHECK( section_denom_name != NULL, XLAL_EFAILED );
          char keyword[64];
          char comment[1024];
          snprintf( keyword, sizeof( keyword ), "tau %s %s_%s", cpu_sections[i].name, statistic_name, section_denom_name );
          snprintf( comment, sizeof( comment ), "%s timing constant", statistic_name );
          XLAL_CHECK( XLALFITSHeaderWriteREAL8( file, keyword, statistic_cpu_time / section_denom / statistic_denom, comment ) == XLAL_SUCCESS, XLAL_EFUNC );
        }

        // Compute remaining CPU time in timing section
        section_cpu_time_remain -= statistic_cpu_time;

      }
    }

    if ( statistics_computed ) {

      // Write remaining CPU time in timing section
      char keyword[64];
      char comment[1024];
      snprintf( keyword, sizeof( keyword ), "cpu %s other", cpu_sections[i].name );
      snprintf( comment, sizeof( comment ), "other %s CPU time", cpu_sections[i].comment );
      XLAL_CHECK( XLALFITSHeaderWriteREAL8( file, keyword, section_cpu_time_remain, comment ) == XLAL_SUCCESS, XLAL_EFUNC );

    } else if ( section_denom > 0 ) {

      // Write timing constant for overall section
      char keyword[64];
      char comment[1024];
      snprintf( keyword, sizeof( keyword ), "tau %s_%s", cpu_sections[i].name, section_denom_name );
      snprintf( comment, sizeof( comment ), "%s timing constant", cpu_sections[i].comment );
      XLAL_CHECK( XLALFITSHeaderWriteREAL8( file, keyword, section_cpu_time / section_denom, comment ) == XLAL_SUCCESS, XLAL_EFUNC );

    }

  }

  return XLAL_SUCCESS;

}

// Local Variables:
// c-file-style: "linux"
// c-basic-offset: 2
// End:
