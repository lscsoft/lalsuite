//
// Copyright (C) 2020 Karl Wette
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
// MA  02110-1301  USA
//

///
/// \file
/// \ingroup lalpulsar_bin_Weave
///

#include "Weave.h"
#include "SetupData.h"
#include "OutputResults.h"

#include <lal/LogPrintf.h>
#include <lal/UserInput.h>

// User input choices for F-statistic consistency veto
extern const UserChoices WeaveFstatConsistencyVetoChoices;

// User help string for F-statistic consistency veto
extern const char *const WeaveFstatConsistencyVetoHelpString;

typedef struct {
  // Statistics to apply F-stat consistency veto to
  WeaveStatisticType apply_Fstat_consistency_veto;
} Fstat_consistency_veto_params;

static BOOLEAN Fstat_consistency_veto(
  const REAL4 twoF,
  const REAL4 twoF_det[PULSAR_MAX_DETECTORS]
)
{

  // Compute maximum twoF over detectors
  REAL4 twoF_det_max = 0.0;
  size_t ndet = 0;
  for ( size_t X = 0; X < PULSAR_MAX_DETECTORS; ++X ) {
    if ( isfinite( twoF_det[X] ) && twoF_det[X] > 0 ) {
      ++ndet;
      if ( twoF_det[X] > twoF_det_max ) {
        twoF_det_max = twoF_det[X];
      }
    }
  }

  // Always pass if there is only one detector, or no per-detector F-statistics
  if ( ndet < 2 ) {
    return 1;
  }

  // F-statistics are consitent if multi-detector twoF is
  // greater than maximum of single-detector twoFs
  return twoF > twoF_det_max;

}

static BOOLEAN results_toplist_item_select_Fstat_consistency_veto(
  const WeaveResultsToplistItem *item,
  const WeaveStatisticsParams *params,
  const void *extra_params
)
{

  // Unpack extra params
  const Fstat_consistency_veto_params *par = ( const Fstat_consistency_veto_params * ) extra_params;
  const WeaveStatisticType apply = par->apply_Fstat_consistency_veto;

  // Iterate over stages
  for ( size_t s = 0; s < 2; ++s ) {
    const WeaveStatisticType stats = params->statistics_to_output[s];
    const WeaveStatisticsValues *const stage = &item->stage[s];

    // Apply F-stat consistency veto to coherent F-statistics
    if (
      ( stats & WEAVE_STATISTIC_COH2F ) &&
      ( stats & WEAVE_STATISTIC_COH2F_DET ) &&
      ( apply & WEAVE_STATISTIC_COH2F )
    ) {
      for ( UINT4 i = 0; i < params->nsegments; ++i ) {
        REAL4 twoF = stage->coh2F[i];
        REAL4 twoF_det[PULSAR_MAX_DETECTORS] = {0.0};
        for ( size_t X = 0; X < PULSAR_MAX_DETECTORS; ++X ) {
          if ( stage->coh2F_det[X] != NULL ) {
            twoF_det[X] = stage->coh2F_det[X][i];
          }
        }
        if ( !Fstat_consistency_veto( twoF, twoF_det ) ) {
          return 0;
        }
      }
    }

    // Apply F-stat consistency veto to summed F-statistic
    if (
      ( stats & WEAVE_STATISTIC_SUM2F ) &&
      ( stats & WEAVE_STATISTIC_SUM2F_DET ) &&
      ( apply & WEAVE_STATISTIC_SUM2F )
    ) {
      if ( !Fstat_consistency_veto( stage->sum2F, stage->sum2F_det ) ) {
        return 0;
      }
    }

    // Apply F-stat consistency veto to mean F-statistic
    if (
      ( stats & WEAVE_STATISTIC_MEAN2F ) &&
      ( stats & WEAVE_STATISTIC_MEAN2F_DET ) &&
      ( apply & WEAVE_STATISTIC_MEAN2F )
    ) {
      if ( !Fstat_consistency_veto( stage->mean2F, stage->mean2F_det ) ) {
        return 0;
      }
    }

  }

  return 1;

}

int main( int argc, char *argv[] )
{

  // Set help information
  lalUserVarHelpBrief = "concatenate result files produced by lalpulsar_Weave";

  ////////// Parse user input //////////

  // Initialise user input variables
  struct uvar_type {
    CHAR *output_result_file;
    LALStringVector *input_result_files;
    int Fstat_consistency_veto;
    UINT4 toplist_limit;
  } uvar_struct = {
    .Fstat_consistency_veto = WEAVE_STATISTIC_NONE,
    .toplist_limit = 0,
  };
  struct uvar_type *const uvar = &uvar_struct;

  // Register user input variables:
  //
  // - General
  //
  XLALRegisterUvarMember(
    input_result_files, STRINGVector, 'i', REQUIRED,
    "Input result files produced by lalpulsar_Weave for concatenation. "
  );
  XLALRegisterUvarMember(
    output_result_file, STRING, 'o', REQUIRED,
    "Output concatenated result file. "
  );

  XLALRegisterUvarAuxDataMember(
    Fstat_consistency_veto, UserFlag, &WeaveFstatConsistencyVetoChoices, 'F', OPTIONAL,
    "Apply the multi-detector F-statistic consistency veto to the given statistics:\n"
    "%s", WeaveFstatConsistencyVetoHelpString
  );
  XLALRegisterUvarMember(
    toplist_limit, UINT4, 'n', OPTIONAL,
    "Maximum number of candidates to return in an output toplist; if 0, all candidates are returned. "
  );

  // Parse user input
  XLAL_CHECK_MAIN( xlalErrno == 0, XLAL_EFUNC, "A call to XLALRegisterUvarMember() failed" );
  BOOLEAN should_exit = 0;
  XLAL_CHECK_MAIN( XLALUserVarReadAllInput( &should_exit, argc, argv, lalPulsarVCSInfoList ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Check user input:
  //
  // - General
  //

  // Exit if required
  if ( should_exit ) {
    return EXIT_FAILURE;
  }
  LogPrintf( LOG_NORMAL, "Parsed user input successfully\n" );

  ////////// Set up F-statistic consistency veto //////////

  WeaveResultsToplistItemSelect toplist_select_fcn = NULL;
  Fstat_consistency_veto_params toplist_select_extra_params = {
    .apply_Fstat_consistency_veto = uvar->Fstat_consistency_veto,
  };
  if ( uvar->Fstat_consistency_veto != 0 ) {
    toplist_select_fcn = results_toplist_item_select_Fstat_consistency_veto;
  }

  ////////// Concatenate output results //////////

  // Number of computed coherent results, and number of coherent and semicoherent templates
  UINT8 coh_nres = 0;
  UINT8 coh_ntmpl = 0;
  UINT8 semi_ntmpl = 0;

  // Total wall and CPU time
  REAL8 wall_total = 0;
  REAL8 cpu_total = 0;

  // Output results
  WeaveOutputResults *out = NULL;

  // Concatenate input result files
  for ( size_t i = 0; i < uvar->input_result_files->length; ++i ) {
    LogPrintf( LOG_NORMAL, "Opening input result file '%s' for reading ...\n", uvar->input_result_files->data[i] );
    FITSFile *file = XLALFITSFileOpenRead( uvar->input_result_files->data[i] );
    XLAL_CHECK_MAIN( file != NULL, XLAL_EFUNC, "Could not open input result file '%s'", uvar->input_result_files->data[i] );
    {
      UINT8 coh_nres_i = 0;
      XLAL_CHECK_MAIN( XLALFITSHeaderReadUINT8( file, "ncohres", &coh_nres_i ) == XLAL_SUCCESS, XLAL_EFUNC );
      coh_nres += coh_nres_i;
    }
    {
      UINT8 coh_ntmpl_i = 0;
      XLAL_CHECK_MAIN( XLALFITSHeaderReadUINT8( file, "ncohtpl", &coh_ntmpl_i ) == XLAL_SUCCESS, XLAL_EFUNC );
      coh_ntmpl += coh_ntmpl_i;
    }
    {
      UINT8 semi_ntmpl_i = 0;
      XLAL_CHECK_MAIN( XLALFITSHeaderReadUINT8( file, "nsemitpl", &semi_ntmpl_i ) == XLAL_SUCCESS, XLAL_EFUNC );
      semi_ntmpl += semi_ntmpl_i;
    }
    {
      REAL8 wall_total_i = 0;
      XLAL_CHECK( XLALFITSHeaderReadREAL8( file, "wall total", &wall_total_i ) == XLAL_SUCCESS, XLAL_EFUNC );
      wall_total += wall_total_i;
    }
    {
      REAL8 cpu_total_i = 0;
      XLAL_CHECK( XLALFITSHeaderReadREAL8( file, "cpu total", &cpu_total_i ) == XLAL_SUCCESS, XLAL_EFUNC );
      cpu_total += cpu_total_i;
    }
    XLAL_CHECK_MAIN( XLALWeaveOutputResultsReadAppend( file, &out, uvar->toplist_limit, toplist_select_fcn, &toplist_select_extra_params ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLALFITSFileClose( file );
    LogPrintf( LOG_NORMAL, "Closed input result file '%s'\n", uvar->input_result_files->data[i] );
  }

  // Write output concatenated result file
  {
    LogPrintf( LOG_NORMAL, "Opening output result file '%s' for writing ...\n", uvar->output_result_file );
    FITSFile *file = XLALFITSFileOpenWrite( uvar->output_result_file );
    XLAL_CHECK_MAIN( file != NULL, XLAL_EFUNC, "Could not open output result file '%s'", uvar->output_result_file );
    XLAL_CHECK_MAIN( XLALFITSHeaderWriteUINT8( file, "ncohres", coh_nres, "number of computed coherent results" ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK_MAIN( XLALFITSHeaderWriteUINT8( file, "ncohtpl", coh_ntmpl, "number of coherent templates" ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK_MAIN( XLALFITSHeaderWriteUINT8( file, "nsemitpl", semi_ntmpl, "number of semicoherent templates" ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( XLALFITSHeaderWriteREAL8( file, "wall total", wall_total, "total wall time" ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( XLALFITSHeaderWriteREAL8( file, "cpu total", cpu_total, "total CPU time" ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK_MAIN( XLALWeaveOutputResultsWrite( file, out ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLALFITSFileClose( file );
    LogPrintf( LOG_NORMAL, "Closed output result file '%s'\n", uvar->output_result_file );
  }

  ////////// Cleanup memory and exit //////////

  // Cleanup memory from output results
  XLALWeaveOutputResultsDestroy( out );

  // Cleanup memory from user input
  XLALDestroyUserVars();

  // Check for memory leaks
  LALCheckMemoryLeaks();

  return EXIT_SUCCESS;

}

// Local Variables:
// c-file-style: "linux"
// c-basic-offset: 2
// End:
