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

///
/// \file
/// \ingroup lalapps_pulsar_Weave
///

#include "Weave.h"
#include "LALAppsVCSInfo.h"

#include <lal/LALString.h>
#include <lal/LogPrintf.h>
#include <lal/UserInput.h>
#include <lal/GSLHelpers.h>
#include <lal/LALInitBarycenter.h>

int main( int argc, char *argv[] )
{

  // Check VCS information
  XLAL_CHECK_MAIN( XLALAppsVCSInfoCheck() == XLAL_SUCCESS, XLAL_EFUNC );

  // Set help information
  lalUserVarHelpBrief = "create setup file for use with lalapps_Weave";

  ////////// Parse user input //////////

  // Initialise user input variables
  struct uvar_type {
    CHAR *segment_list, *ephem_earth, *ephem_sun, *detector_motion, *output_file;
    INT4 segment_count, spindowns;
    LALStringVector *detectors;
    LIGOTimeGPS ref_time;
    LIGOTimeGPSRange first_segment;
    REAL8 segment_gap;
  } uvar_struct = {
    .detector_motion = XLALStringDuplicate( "spin+orbit" ),
    .ephem_earth = XLALStringDuplicate( "earth00-19-DE405.dat.gz" ),
    .ephem_sun = XLALStringDuplicate( "sun00-19-DE405.dat.gz" ),
    .segment_count = 1,
    .spindowns = 1,
  };
  struct uvar_type *const uvar = &uvar_struct;

  // Register user input variables:
  //
  // - Segment list input/generation
  //
  XLALRegisterUvarMember(
    segment_list, STRING, 'L', OPTIONAL,
    "Loads the start and end times of each segment from this file. "
    "Format is:\n"
    "  # comment\n  <segment-start-time-GPS> <segment-end-time-GPS>\n  ..."
    );
  XLALRegisterUvarMember(
    first_segment, EPOCHRange, 't', OPTIONAL,
    "Generate segments; the range of the first segment is specified by this option. "
    );
  XLALRegisterUvarMember(
    segment_count, INT4, 'n', OPTIONAL,
    "Generate this many segments by translating the first segment in time by its span, i.e. so that segments are contiguous and non-overlapping. "
    );
  XLALRegisterUvarMember(
    segment_gap, REAL8, 'g', DEVELOPER,
    "When generating segments, increase the translation of the first segment by this amount (in seconds). "
    "A positive value gives non-contiguous segments; a negative value gives overlapping segments. "
    );
  //
  // - Parameter-space metric computation
  //
  XLALRegisterUvarMember(
    ref_time, EPOCH, 'r', NODEFAULT,
    "Reference time for the search, including the parameter-space metrics computed here, and the parameter space and output of lalapps_Weave. "
    "If omitted, the mid-point between the start of the first segment and the end of the last segment is used. "
    );
  XLALRegisterUvarMember(
    detectors, STRINGVector, 'I', REQUIRED,
    "Comma-separated list of 2-character detector names (e.g. H1,L1,...) for which the parameter-space metrics are computed. "
    );
  XLALRegisterUvarMember(
    detector_motion, STRING, 'm', DEVELOPER,
    "Specify what detector motion to assume when computing the parameter-space metrics. "
    "The only interesting options are:\n"
    " - 'spin+orbit'       use the full ephemeris of the Earth's orbit;\n"
    " - 'spin+ptoleorbit': use a Ptolemaic approximation of the Earth's orbit."
    );
  XLALRegisterUvarMember(
    ephem_earth, STRING, 'E', DEVELOPER,
    "Earth ephemeris file, used to compute the parameter-space metrics and by lalapps_Weave. "
    );
  XLALRegisterUvarMember(
    ephem_sun, STRING, 'S', DEVELOPER,
    "Sun ephemeris file, used to compute the parameter-space metrics and by lalapps_Weave. "
    );
  XLALRegisterUvarMember(
    spindowns, INT4, 's', OPTIONAL,
    "Maximum number of spindowns for which the parameter-space metrics are computed; must be at least 1. "
    "This option limits the size of the spindown parameter space give to lalapps_Weave. "
    );
  //
  // - Output
  //
  XLALRegisterUvarMember(
    output_file, STRING, 'o', REQUIRED,
    "Output file while stores the segment list, parameter-space metrics, and other data required by lalapps_Weave, e.g. ephemerides. "
    );

  // Parse user input
  XLAL_CHECK_MAIN( xlalErrno == 0, XLAL_EFUNC, "A call to XLALRegisterUvarMember() failed" );
  BOOLEAN should_exit = 0;
  XLAL_CHECK_MAIN( XLALUserVarReadAllInput( &should_exit, argc, argv ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Check user input:
  //
  // - Segment list input/generation
  //
  XLALUserVarCheck( &should_exit,
                    UVAR_SET2( segment_list, first_segment ) == 1,
                    "Exactly one of " UVAR_STR2OR( segment_list, first_segment ) " must be specified" );
  XLALUserVarCheck( &should_exit,
                    !UVAR_ALLSET2( segment_list, segment_gap ),
                    "At most one of " UVAR_STR2AND( segment_list, segment_gap ) " may be specified" );
  XLALUserVarCheck( &should_exit,
                    !UVAR_ALLSET2( segment_list, segment_count ),
                    "At most one of " UVAR_STR2AND( segment_list, segment_count ) " may be specified" );
  XLALUserVarCheck( &should_exit,
                    uvar->segment_count > 0,
                    UVAR_STR( segment_count ) " must be strictly positive" );
  //
  // - Parameter-space metric computation
  //
  XLALUserVarCheck( &should_exit,
                    uvar->spindowns >= 0,
                    UVAR_STR( spindowns ) " must be positive" );

  // Exit if required
  if ( should_exit ) {
    return EXIT_FAILURE;
  }
  LogPrintf( LOG_NORMAL, "Parsed user input successfully\n" );

  ////////// Create setup data //////////

  // Initialise setup data
  WeaveSetupData XLAL_INIT_DECL( setup );

  // Copy and sort list of detector names
  setup.detectors = XLALCopyStringVector( uvar->detectors );
  XLAL_CHECK_MAIN( setup.detectors != NULL, XLAL_EFUNC );
  XLAL_CHECK_MAIN( XLALSortStringVector( setup.detectors ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Load ephemerides
  LogPrintf( LOG_NORMAL, "Loading ephemerides from '%s' and '%s' ...\n", uvar->ephem_earth, uvar->ephem_sun );
  setup.ephemerides = XLALInitBarycenter( uvar->ephem_earth, uvar->ephem_sun );
  XLAL_CHECK_MAIN( setup.ephemerides != NULL, XLAL_EFUNC );

  // Create segment list
  if ( UVAR_SET( segment_list ) ) {

    // Load segment list from file
    LogPrintf( LOG_NORMAL, "Loading segment list from '%s' ...\n", uvar->segment_list );
    setup.segments = XLALReadSegmentsFromFile( uvar->segment_list );
    XLAL_CHECK_MAIN( setup.segments != NULL, XLAL_EFUNC );

  } else {

    // Generate segment list
    LogPrintf( LOG_NORMAL, "Generating segment list; first segment = [%" LAL_GPS_FORMAT ", %" LAL_GPS_FORMAT "] GPS, segment gap = %g sec ...\n", LAL_GPS_PRINT( uvar->first_segment[0] ), LAL_GPS_PRINT( uvar->first_segment[1] ), uvar->segment_gap );
    setup.segments = XLALSegListCreate();
    XLAL_CHECK_MAIN( setup.segments != NULL, XLAL_EFUNC );

    // Set first segment
    LALSeg seg;
    XLALSegSet( &seg, &uvar->first_segment[0], &uvar->first_segment[1], 0 );

    // Increment to add to each successive segment
    const REAL8 dt = XLALGPSDiff( &seg.end, &seg.start ) + uvar->segment_gap;

    // Create segments
    for ( INT4 n = 0; n < uvar->segment_count; ++n ) {
      XLAL_CHECK_MAIN( XLALSegListAppend( setup.segments, &seg ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLALGPSAdd( &seg.start, dt );
      XLALGPSAdd( &seg.end, dt );
    }

  }

  // Compute segment list range
  LIGOTimeGPS segments_start, segments_end;
  XLAL_CHECK_MAIN( XLALSegListRange( setup.segments, &segments_start, &segments_end ) == XLAL_SUCCESS, XLAL_EFUNC );
  LogPrintf( LOG_NORMAL, "Segment list range = [%" LAL_GPS_FORMAT ", %" LAL_GPS_FORMAT "] GPS, segment count = %i\n", LAL_GPS_PRINT( segments_start ), LAL_GPS_PRINT( segments_end ), setup.segments->length );

  // Set reference time
  if ( UVAR_SET( ref_time ) ) {
    setup.ref_time = uvar->ref_time;
  } else {
    setup.ref_time = segments_start;
    XLALGPSAdd( &setup.ref_time, 0.5 * XLALGPSDiff( &segments_end, &segments_start ) );
  }
  LogPrintf( LOG_NORMAL, "Reference time = %" LAL_GPS_FORMAT " GPS\n", LAL_GPS_PRINT( setup.ref_time ) );

  // Restrict ephemerides to range of segment list +/- 1 day
  {
    LIGOTimeGPS ephem_start = segments_start, ephem_end = segments_end;
    XLALGPSAdd( &ephem_start, -LAL_DAYSID_SI );
    XLALGPSAdd( &ephem_end, +LAL_DAYSID_SI );
    XLAL_CHECK( XLALRestrictEphemerisData( setup.ephemerides, &ephem_start, &ephem_end ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Parse list of detectors
  MultiLALDetector detector_info;
  {
    char *detectors_string = XLALConcatStringVector( uvar->detectors, "," );
    XLAL_CHECK_MAIN( detectors_string != NULL, XLAL_EFUNC );
    XLAL_CHECK_MAIN( XLALParseMultiLALDetector( &detector_info, uvar->detectors ) == XLAL_SUCCESS, XLAL_EINVAL, "Invalid value '%s' for " UVAR_STR( detectors ), detectors_string );
    XLALFree( detectors_string );
  }

  // Parse detector motion string
  const DetectorMotionType detector_motion = XLALParseDetectorMotionString( uvar->detector_motion );
  XLAL_CHECK_MAIN( xlalErrno == 0, XLAL_EINVAL, "Invalid value '%s' for " UVAR_STR( detector_motion ), uvar->detector_motion );

  // Compute reduced supersky metrics at fiducial frequency of 100Hz
  // - Fiducial frequency is stored in coordinate transform data, so
  //   metrics can later be rescaled by search code
  LogPrintf( LOG_NORMAL, "Computing reduced supersky metrics ...\n" );
  const double fiducial_freq = 100.0;
  setup.metrics = XLALComputeSuperskyMetrics( uvar->spindowns, &setup.ref_time, setup.segments, fiducial_freq, &detector_info, NULL, detector_motion, setup.ephemerides );
  XLAL_CHECK_MAIN( setup.metrics != NULL, XLAL_EFUNC );
  LogPrintf( LOG_NORMAL, "Finished computing reduced supersky metrics\n" );

  ////////// Output setup data //////////

  // Open output file
  LogPrintf( LOG_NORMAL, "Opening output file '%s' for writing ...\n", uvar->output_file );
  FITSFile *file = XLALFITSFileOpenWrite( uvar->output_file );
  XLAL_CHECK_MAIN( file != NULL, XLAL_EFUNC );
  XLAL_CHECK_MAIN( XLALFITSFileWriteVCSInfo( file, lalAppsVCSInfoList ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( XLALFITSFileWriteUVarCmdLine( file ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Write setup data
  XLAL_CHECK_MAIN( XLALWeaveSetupDataWrite( file, &setup ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Close output file
  XLALFITSFileClose( file );
  LogPrintf( LOG_NORMAL, "Closed output file '%s'\n", uvar->output_file );

  ////////// Cleanup memory and exit //////////

  // Cleanup memory from setup data
  XLALWeaveSetupDataClear( &setup );

  // Cleanup memory from user input
  XLALDestroyUserVars();

  // Check for memory leaks
  LALCheckMemoryLeaks();

  LogPrintf( LOG_NORMAL, "Finished successfully!\n" );

  return EXIT_SUCCESS;

}
