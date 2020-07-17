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
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA 02111-1307 USA
//

///
/// \file
/// \ingroup lalapps_pulsar_Weave
///

#include "Weave.h"
#include "SetupData.h"
#include "OutputResults.h"

#include <lal/LogPrintf.h>
#include <lal/UserInput.h>

int main( int argc, char *argv[] )
{

  // Set help information
  lalUserVarHelpBrief = "concatenate result files produced by lalapps_Weave";

  ////////// Parse user input //////////

  // Initialise user input variables
  struct uvar_type {
    CHAR *output_result_file;
    LALStringVector *input_result_files;
    UINT4 toplist_limit;
  } uvar_struct = {
    .output_result_file = NULL,
  };
  struct uvar_type *const uvar = &uvar_struct;

  // Register user input variables:
  //
  // - General
  //
  XLALRegisterUvarMember(
    input_result_files, STRINGVector, 'i', REQUIRED,
    "Input result files produced by lalapps_Weave for concatenation. "
    );
  XLALRegisterUvarMember(
    output_result_file, STRING, 'o', REQUIRED,
    "Output concatenated result file. "
    );
  XLALRegisterUvarMember(
    toplist_limit, UINT4, 'n', OPTIONAL,
    "Maximum number of candidates to return in an output toplist; if 0, all candidates are returned. "
    );

  // Parse user input
  XLAL_CHECK_MAIN( xlalErrno == 0, XLAL_EFUNC, "A call to XLALRegisterUvarMember() failed" );
  BOOLEAN should_exit = 0;
  XLAL_CHECK_MAIN( XLALUserVarReadAllInput( &should_exit, argc, argv, lalAppsVCSInfoList ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Check user input:
  //
  // - General
  //

  // Exit if required
  if ( should_exit ) {
    return EXIT_FAILURE;
  }
  LogPrintf( LOG_NORMAL, "Parsed user input successfully\n" );

  ////////// Concatenate output results //////////

  // Initialise output results
  WeaveOutputResults *out = NULL;

  // Concatenate input result files
  for ( size_t i = 0; i < uvar->input_result_files->length; ++i ) {
    LogPrintf( LOG_NORMAL, "Opening input result file '%s' for reading ...\n", uvar->input_result_files->data[i] );
    FITSFile *file = XLALFITSFileOpenRead( uvar->input_result_files->data[i] );
    XLAL_CHECK_MAIN( file != NULL, XLAL_EFUNC, "Could not open input result file '%s'", uvar->input_result_files->data[i] );
    XLAL_CHECK_MAIN( XLALWeaveOutputResultsReadAppend( file, &out, uvar->toplist_limit ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLALFITSFileClose( file );
    LogPrintf( LOG_NORMAL, "Closed input result file '%s'\n", uvar->input_result_files->data[i] );
  }

  // Write output concatenated result file
  {
    LogPrintf( LOG_NORMAL, "Opening output result file '%s' for writing ...\n", uvar->output_result_file );
    FITSFile *file = XLALFITSFileOpenWrite( uvar->output_result_file );
    XLAL_CHECK_MAIN( file != NULL, XLAL_EFUNC, "Could not open output result file '%s'", uvar->output_result_file );
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
