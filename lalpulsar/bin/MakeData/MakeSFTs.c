/*
 * Copyright (C) 2011, 2015, 2022 Karl Wette
 * Copyright (C) 2007, 2009, 2011, 2015, 2016 Bernd Machenschalk
 * Copyright (C) 2005--2007, 2012 Gregory Mendell
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA
 */

/**
 * \file
 * \ingroup lalpulsar_bin_SFTTools
 * \author Karl Wette, Evan Goetz, Gregory Mendell, Bernd Machenschalk, Xavier Siemens, Bruce Allen
 * \brief Generate SFTs
 */

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>

#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/UserInput.h>
#include <lal/LogPrintf.h>
#include <lal/Units.h>
#include <lal/TimeSeries.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/Window.h>
#include <lal/RealFFT.h>
#include <lal/SFTfileIO.h>
#include <lal/LALFrStream.h>
#include <lal/LALVCSInfo.h>
#include <lal/LALPulsarVCSInfo.h>

int main( int argc, char *argv[] )
{

  // Set help information
  lalUserVarHelpBrief = "generate SFTs";

  ////////// Parse user input //////////

  // Initialise user input variables
  struct uvar_type {
    char *frame_cache;
    char *channel_name;
    INT4 gps_start_time;
    INT4 gps_end_time;
    INT4 sft_duration;
    REAL8 overlap_fraction;
    REAL8 high_pass_freq;
    char *window_type;
    REAL8 window_radius;
    REAL8 start_freq;
    REAL8 band;
    char *sft_write_path;
    UINT4 observing_run;
    char *observing_kind;
    UINT4 observing_version;
    char *comment_field;
  } uvar_struct = {
    .sft_duration = 1800,
    .overlap_fraction = 0,
    .high_pass_freq = 0,
    .window_type = XLALStringDuplicate("tukey"),
    .window_radius = 0.001,
    .start_freq = 48,
    .band = 2000,
    .sft_write_path = XLALStringDuplicate("."),
  };
  struct uvar_type *const uvar = &uvar_struct;

  // Register user input variables
  //
  // - Input data
  //
  lalUserVarHelpOptionSubsection = "Input data";
  XLALRegisterUvarMember(
    frame_cache, STRING, 'C', REQUIRED,
    "Path to frame cache file to read frames from. "
    );
  XLALRegisterUvarMember(
    channel_name, STRING, 'N', REQUIRED,
    "Name of channel to read within a frame. "
    );
  XLALRegisterUvarMember(
    gps_start_time, INT4, 's', REQUIRED,
    "GPS time to start generating SFTs. "
    );
  XLALRegisterUvarMember(
    gps_end_time, INT4, 'e', REQUIRED,
    "GPS time to end generating SFTs. "
    );
  //
  // - SFT generation
  //
  lalUserVarHelpOptionSubsection = "SFT generation";
  XLALRegisterUvarMember(
    sft_duration, INT4, 't', OPTIONAL,
    "SFT duration in seconds. "
    );
  XLALRegisterUvarMember(
    overlap_fraction, REAL8, 'P', OPTIONAL,
    "Fraction of SFT duration to overlap SFTs by. "
    "Example: use " UVAR_STR( overlap_fraction ) "=0.5 with " UVAR_STR( window_type ) "=hann windows. " );
  XLALRegisterUvarMember(
    high_pass_freq, REAL8, 'f', REQUIRED,
    "High pass filtering frequency in Hertz. " );
  XLALRegisterUvarMember(
    window_type, STRING, 'w', OPTIONAL,
    "Window to apply to SFTs. "
    );
  XLALRegisterUvarMember(
    window_radius, REAL8, 'r', OPTIONAL,
    "Parameter (if required) of window to apply to SFTs. "
    );
  XLALRegisterUvarMember(
    start_freq, REAL8, 'F', OPTIONAL,
    "Start frequency of the SFTs, in Hertz. "
    );
  XLALRegisterUvarMember(
    band, REAL8, 'B', OPTIONAL,
    "Frequency band of the SFTs, in Hertz. "
    );
  //
  // - SFT output
  //
  lalUserVarHelpOptionSubsection = "SFT output";
  XLALRegisterUvarMember(
    sft_write_path, STRING, 'p', OPTIONAL,
    "Path to write SFTs to. "
    );
  XLALRegisterUvarMember(
    observing_run, UINT4, 'O', REQUIRED,
    "Observing run SFTs are generated from. "
    );
  XLALRegisterUvarMember(
    observing_kind, STRING, 'K', REQUIRED,
    "Kind of SFTs being generated: 'RUN', 'AUX', 'SIM', or 'DEV'. "
    );
  XLALRegisterUvarMember(
    observing_version, UINT4, 'V', REQUIRED,
    "Version number of the SFT production. "
    );
  XLALRegisterUvarMember(
    comment_field, STRING, 'c', OPTIONAL,
    "Comment for SFT header. "
    );
  //
  // - Defunct options
  //
  XLALRegisterNamedUvar(
    NULL, "frame_struct_type", STRING, 'u', DEFUNCT,
    "No longer required; the frame channel type is determined automatically. "
    );
  XLALRegisterNamedUvar(
    NULL, "ht_data", BOOLEAN, 'H', DEFUNCT,
    "No longer required; the frame channel type is determined automatically. "
    );
  XLALRegisterNamedUvar(
    NULL, "ifo", STRING, 'i', DEFUNCT,
    "No longer required; the detector prefix is deduced from the channel name. "
    );
  XLALRegisterNamedUvar(
    NULL, "make_gps_dirs", BOOLEAN, 'D', DEFUNCT,
    "No longer supported. "
    );
  XLALRegisterNamedUvar(
    NULL, "make_tmp_file", BOOLEAN, 'Z', DEFUNCT,
    "Default behaviour. "
    );
  XLALRegisterNamedUvar(
    NULL, "misc_desc", STRING, 'X', DEFUNCT,
    "No longer supported in the public SFT filename specification. "
    );
  XLALRegisterNamedUvar(
    NULL, "sft_version", INT4, 0, DEFUNCT,
    "No longer supported. "
    );
  XLALRegisterNamedUvar(
    NULL, "use_single", BOOLEAN, 'S', DEFUNCT,
    "No longer supported. "
    );

  // Parse user input
  XLAL_CHECK_MAIN( xlalErrno == 0, XLAL_EFUNC, "A call to XLALRegisterUvarMember() failed" );
  BOOLEAN should_exit = 0;
  XLAL_CHECK_MAIN( XLALUserVarReadAllInput( &should_exit, argc, argv, lalPulsarVCSInfoList ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( should_exit ) {
    return EXIT_FAILURE;
  }

  // Check user input:
  //
  // - Input data
  //
  XLALUserVarCheck( &should_exit,
                    uvar->gps_start_time > 0,
                    UVAR_STR( gps_start_time ) " must be strictly positive" );
  XLALUserVarCheck( &should_exit,
                    uvar->gps_end_time > 0,
                    UVAR_STR( gps_end_time ) " must be strictly positive" );
  //
  // - SFT generation
  //
  XLALUserVarCheck( &should_exit,
                    uvar->sft_duration > 0,
                    UVAR_STR( sft_duration ) " must be strictly positive" );
  XLALUserVarCheck( &should_exit,
                    uvar->gps_end_time - uvar->gps_start_time >= uvar->sft_duration,
                    "Cannot fit an SFT of duration %d between GPS times %d and %d\n",
                    uvar->sft_duration, uvar->gps_start_time, uvar->gps_end_time );
  XLALUserVarCheck( &should_exit,
                    0 <= uvar->overlap_fraction && uvar->overlap_fraction < 1,
                    UVAR_STR( overlap_fraction ) " must be within range [0.0, 1.0)" );
  XLALUserVarCheck( &should_exit,
                    uvar->high_pass_freq >= 0,
                    UVAR_STR( high_pass_freq ) " must be positive" );
  if ( strcmp( uvar->window_type, "0" ) == 0 ) {
    XLALUserVarCheck( &should_exit, 0,
                      UVAR_STR( window_type ) "=0 is deprecated; use " UVAR_STR( window_type ) "=rectangular" );
  } else if ( strcmp( uvar->window_type, "1" ) == 0 ) {
    XLALUserVarCheck( &should_exit, 0,
                      UVAR_STR( window_type ) "=1 is the old 'Matlab style Tukey window'; this EXACT window is "
                      "no longer supported BUT " UVAR_STR( window_type ) "=tukey will give you almost exactly "
                      "the same window; see https://dcc.ligo.org/LIGO-T040164/public for details" );
  } else if ( strcmp( uvar->window_type, "2" ) == 0 ) {
    XLALUserVarCheck( &should_exit, 0,
                      UVAR_STR( window_type ) "=2 is the old 'make_sfts.c Tukey window'; this window is "
                      "no longer supported. You probably want to use " UVAR_STR( window_type ) "=tukey, "
                      "although this will NOT give you exactly the same window" );
  } else if ( strcmp( uvar->window_type, "3" ) == 0 ) {
    XLALUserVarCheck( &should_exit, 0,
                      UVAR_STR( window_type ) "=3 is deprecated; use " UVAR_STR( window_type ) "=hann" );
  } else if ( XLALUserVarWasSet( &uvar->window_type ) ) {
    XLALUserVarCheck( &should_exit,
                      XLALCheckNamedWindow( uvar->window_type, XLALUserVarWasSet( &uvar->window_radius ) ) == XLAL_SUCCESS,
                      "Invalid/inconsistent " UVAR_STR( window_type ) " and/or " UVAR_STR( window_radius ) );
  }
  XLALUserVarCheck( &should_exit,
                    uvar->window_radius >= 0,
                    UVAR_STR( window_radius ) " must be positive" );
  XLALUserVarCheck( &should_exit,
                    uvar->start_freq >= 0,
                    UVAR_STR( start_freq ) " must be positive" );
  XLALUserVarCheck( &should_exit,
                    uvar->band > 0,
                    UVAR_STR( band ) " must be strictly positive" );
  //
  // - SFT output
  //
  XLALUserVarCheck( &should_exit,
                    uvar->observing_run > 0,
                    UVAR_STR( observing_run ) " must be strictly positive" );
  XLALUserVarCheck( &should_exit,
                    strcmp( uvar->observing_kind, "RUN" ) == 0
                    || strcmp( uvar->observing_kind, "AUX" ) == 0
                    || strcmp( uvar->observing_kind, "SIM" ) == 0
                    || strcmp( uvar->observing_kind, "DEV" ) == 0,
                    UVAR_STR( observing_kind ) " must be one of 'RUN', 'AUX', 'SIM', or 'DEV'" );
  XLALUserVarCheck( &should_exit,
                    uvar->observing_version > 0,
                    UVAR_STR( observing_version ) " must be strictly positive" );

  // Exit if required
  if ( should_exit ) {
    return EXIT_FAILURE;
  }
  LogPrintf( LOG_NORMAL, "Parsed user input successfully\n" );

  ////////// Set up SFT generation //////////

  // Build comment string from program name, VCS information, and full command line option log
  // (including default values); this will include anything passed to --comment-field
  char *SFT_comment = NULL;
  {
    char *vcs_info_str = XLALVCSInfoString( lalPulsarVCSInfoList, 0, NULL );
    XLAL_CHECK_MAIN( vcs_info_str != NULL, XLAL_EFUNC );
    char *cmd_line_str = XLALUserVarGetLogEx( UVAR_LOGFMT_CFGFILE, 0 );
    XLAL_CHECK_MAIN( cmd_line_str != NULL, XLAL_EFUNC );
    SFT_comment = XLALStringAppendFmt( SFT_comment, "Generated by %s\nVersion:\n%s\nOptions:\n%s\n", argv[0], vcs_info_str, cmd_line_str );
    XLAL_CHECK_MAIN( SFT_comment != NULL, XLAL_EFUNC );
    XLALFree( vcs_info_str );
    XLALFree( cmd_line_str );
  }

  // Create frame cache
  LALCache *framecache = XLALCacheImport( uvar->frame_cache );
  XLAL_CHECK_MAIN( framecache != NULL, XLAL_EFUNC,
                   "Failed to open frame cache '%s'", uvar->frame_cache );

  // Open frame stream
  LALFrStream *framestream = XLALFrStreamCacheOpen( framecache );
  XLAL_CHECK_MAIN( framestream != NULL, XLAL_EFUNC,
                   "Failed to open frame stream from cache '%s'", uvar->frame_cache );

  // Set frame stream mode
  {
    const LALFrStreamMode mode = 0
      | LAL_FR_STREAM_VERBOSE_MODE      /* Display warnings and info */
      | LAL_FR_STREAM_IGNOREGAP_MODE    /* Ignore gaps in data */
      | LAL_FR_STREAM_IGNORETIME_MODE   /* Ignore invalid times requested */
      | LAL_FR_STREAM_CHECKSUM_MODE     /* Ensure that file checksums are OK */
      ;
    XLAL_CHECK_MAIN( XLALFrStreamSetMode( framestream, mode ) == XLAL_SUCCESS, XLAL_EFUNC,
                     "Failed to set frame stream mode to %i", mode );
  }

  // Seek frame stream to starting GPS time
  const LIGOTimeGPS gps_start = { .gpsSeconds = uvar->gps_start_time, .gpsNanoSeconds = 0 };
  XLAL_CHECK_MAIN( XLALFrStreamSeek( framestream, &gps_start ) == XLAL_SUCCESS, XLAL_EFUNC,
                   "Failed to seek frame stream to GPS start time %d", gps_start.gpsSeconds );
  LogPrintf( LOG_NORMAL, "Starting SFT generation at GPS time %" LAL_GPS_FORMAT "\n", LAL_GPS_PRINT( gps_start ) );

  // Allocate (empty) time series for SFT data
  REAL8TimeSeries *SFT_time_series = XLALCreateREAL8TimeSeries( uvar->channel_name, &gps_start, 0.0, 0.0, &lalDimensionlessUnit, 0 );
  XLAL_CHECK_MAIN( SFT_time_series != NULL, XLAL_EFUNC,
                   "Failed to allocate SFT time series" );

  // Update time series metadata from frames
  XLAL_CHECK_MAIN( XLALFrStreamGetREAL8TimeSeriesMetadata( SFT_time_series, framestream ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Resize SFT time series to required size
  SFT_time_series = XLALResizeREAL8TimeSeries( SFT_time_series, 0, uvar->sft_duration / SFT_time_series->deltaT );
  XLAL_CHECK_MAIN( SFT_time_series != NULL, XLAL_EFUNC,
                   "Failed to allocate SFT time series to %g elements", uvar->sft_duration / SFT_time_series->deltaT );

  // Create SFT time series window
  const REAL8 window_beta = XLALUserVarWasSet( &uvar->window_radius ) ? uvar->window_radius : 0;
  REAL8Window *SFT_window = XLALCreateNamedREAL8Window( uvar->window_type, window_beta, SFT_time_series->data->length );
  XLAL_CHECK_MAIN( SFT_window != NULL, XLAL_EFUNC,
                   "Failed to allocate SFT time series window of %u elements", SFT_time_series->data->length );

  // Create SFT FFT data vector and plan
  COMPLEX16Vector *SFT_fft_data = XLALCreateCOMPLEX16Vector( SFT_time_series->data->length / 2 + 1 );
  XLAL_CHECK_MAIN( SFT_fft_data != NULL, XLAL_EFUNC,
                   "Failed to allocate SFT FFT data vector of %u elements", SFT_time_series->data->length / 2 + 1 );
  REAL8FFTPlan *SFT_fft_plan = XLALCreateForwardREAL8FFTPlan( SFT_time_series->data->length, 0 );
  XLAL_CHECK_MAIN( SFT_fft_plan != NULL, XLAL_EFUNC,
                   "Failed to allocate SFT FFT plan" );

  // Create SFT type for writing output
  const UINT4 SFT_first_bin = lrint( uvar->start_freq * uvar->sft_duration );
  const UINT4 SFT_bins = lrint( uvar->band * uvar->sft_duration );
  SFTtype *SFT = XLALCreateSFT( SFT_bins );
  XLAL_CHECK_MAIN( SFT != NULL, XLAL_EFUNC,
                   "Failed to allocate SFT of %u bins", SFT_bins );

  ////////// Generate SFTs //////////

  // Read data from frame stream
  while ( 1 ) {

    // Get current GPS time
    LIGOTimeGPS gps_tell;
    XLALFrStreamTell( &gps_tell, framestream );

    // Break if not enough data remaining
    if ( XLALFrStreamEnd( framestream ) ) {
      LogPrintf( LOG_NORMAL, "Reached end of frame stream at GPS time %" LAL_GPS_FORMAT "\n", LAL_GPS_PRINT( gps_tell ) );
      break;
    }

    // Try to read in time series data for the next SFT
    XLAL_CHECK_MAIN( XLALFrStreamGetREAL8TimeSeries( SFT_time_series, framestream ) == XLAL_SUCCESS, XLAL_EFUNC,
                     "XLALFrStreamGetREAL8TimeSeries() failed at GPS time %" LAL_GPS_FORMAT, LAL_GPS_PRINT( gps_tell ) );

    // Break if time series data for next SFT is after GPS end time
    {
      const LIGOTimeGPS gps_end = { .gpsSeconds = SFT_time_series->epoch.gpsSeconds + uvar->sft_duration, .gpsNanoSeconds = 0 };
      if ( gps_end.gpsSeconds > uvar->gps_end_time ) {
        LogPrintf( LOG_NORMAL, "Ending SFT generation at GPS time %" LAL_GPS_FORMAT "\n", LAL_GPS_PRINT( gps_end ) );
        break;
      }
    }

    // Align SFTs on integer GPS seconds
    if ( SFT_time_series->epoch.gpsNanoSeconds > 0 ) {
      SFT_time_series->epoch.gpsNanoSeconds = 0;
      SFT_time_series->epoch.gpsSeconds += 1;
      XLAL_CHECK_MAIN( XLALFrStreamSeek( framestream, &SFT_time_series->epoch ) == XLAL_SUCCESS, XLAL_EFUNC,
                       "Failed to seek frame stream to aligned GPS time %" LAL_GPS_FORMAT, LAL_GPS_PRINT( SFT_time_series->epoch ) );
      continue;
    }

    // High-pass SFT time series data with Butterworth High Pass filter
    if ( uvar->high_pass_freq > 0.0 ) {
      char filtername[] = "Butterworth High Pass";
      PassBandParamStruc filterpar = {
        .name = filtername,
        .nMax = 10,
        .f2   = uvar->high_pass_freq,
        .a2   = 0.5,
        .f1   = -1.0,
        .a1   = -1.0,
      };
      XLAL_CHECK_MAIN( XLALButterworthREAL8TimeSeries( SFT_time_series, &filterpar ) == XLAL_SUCCESS, XLAL_EFUNC,
                       "Failed to apply %s filter to SFT time series data at GPS time %" LAL_GPS_FORMAT, filterpar.name, LAL_GPS_PRINT( gps_tell ) );
    }

    // Window SFT time series data
    XLAL_CHECK_MAIN( XLALUnitaryWindowREAL8Sequence( SFT_time_series->data, SFT_window ) != NULL, XLAL_EFUNC,
                     "Failed to apply window to SFT time series data at GPS time %" LAL_GPS_FORMAT, LAL_GPS_PRINT( gps_tell ) );

    // Fourier transform SFT time series data
    XLAL_CHECK_MAIN( XLALREAL8ForwardFFT( SFT_fft_data, SFT_time_series->data, SFT_fft_plan ) == XLAL_SUCCESS, XLAL_EFUNC,
                     "Failed to Fourier transform SFT time series data at GPS time %" LAL_GPS_FORMAT, LAL_GPS_PRINT( gps_tell ) );

    // Initialise SFT type
    SFT->name[0] = uvar->channel_name[0];
    SFT->name[1] = uvar->channel_name[1] == '0' ? '1' : uvar->channel_name[1];
    SFT->name[2] = 0;
    SFT->epoch = SFT_time_series->epoch;
    SFT->f0 = uvar->start_freq;
    SFT->deltaF = 1.0 / ((REAL8) uvar->sft_duration);

    // Copy and normalise SFT frequency series data into SFT
    const REAL8 SFT_normalisation = SFT_time_series->deltaT;
    for ( UINT4 k = 0; k < SFT_bins; ++k ) {
      const REAL8 SFT_fft_re = creal( SFT_fft_data->data[k + SFT_first_bin] );
      const REAL8 SFT_fft_im = cimag( SFT_fft_data->data[k + SFT_first_bin] );
      SFT->data->data[k] = crectf( SFT_normalisation * SFT_fft_re, SFT_normalisation * SFT_fft_im );
    }

    // Build SFT filename spec
    SFTFilenameSpec XLAL_INIT_DECL(spec);
    XLAL_CHECK_MAIN( XLALFillSFTFilenameSpecStrings( &spec, uvar->sft_write_path, "sft_TO_BE_VALIDATED", NULL, uvar->window_type, NULL, uvar->observing_kind, uvar->channel_name ) == XLAL_SUCCESS, XLAL_EFUNC );
    spec.window_beta = window_beta;
    spec.pubObsRun = uvar->observing_run;
    spec.pubVersion = uvar->observing_version;

    // Write out SFT and retrieve filename
    char *temp_SFT_filename = NULL;
    {
      int retn = XLALWriteSFT2StandardFile( SFT, &spec, SFT_comment );
      temp_SFT_filename = XLALBuildSFTFilenameFromSpec( &spec );
      XLAL_CHECK_MAIN( retn == XLAL_SUCCESS, XLAL_EFUNC,
                       "Failed to write SFT '%s' at GPS time %" LAL_GPS_FORMAT, temp_SFT_filename ? temp_SFT_filename : "<unknown>", LAL_GPS_PRINT( gps_tell ) );
      XLAL_CHECK_MAIN( temp_SFT_filename != NULL, XLAL_EFUNC );
    }

    // Validate SFT
    XLAL_CHECK_MAIN( XLALCheckSFTFileIsValid( temp_SFT_filename ) == XLAL_SUCCESS, XLAL_EFUNC,
                     "Failed to validate SFT '%s' at GPS time %" LAL_GPS_FORMAT, temp_SFT_filename, LAL_GPS_PRINT( gps_tell ) );

    // Move SFT to final filename
    char *final_SFT_filename = XLALStringDuplicate( temp_SFT_filename );
    char *const extn = strrchr( final_SFT_filename, '.' );
    XLAL_CHECK_MAIN( extn != NULL, XLAL_ESYS );
    XLAL_CHECK_MAIN( strlen( extn ) >= 4, XLAL_EFAILED );
    strcpy( extn, ".sft" );
    XLAL_CHECK_MAIN( rename( temp_SFT_filename, final_SFT_filename ) == 0, XLAL_ESYS,
                     "Failed to rename '%s' to '%s': %s", temp_SFT_filename, final_SFT_filename, strerror(errno) );
    LogPrintf( LOG_DEBUG, "Wrote SFT '%s'\n", final_SFT_filename );

    // Overlap SFTs if requested
    if ( uvar->overlap_fraction > 0 ) {
      const REAL8 dt = uvar->overlap_fraction * uvar->sft_duration;
      XLAL_CHECK_MAIN( XLALFrStreamSeekO( framestream, -dt, SEEK_CUR ) == XLAL_SUCCESS, XLAL_EFUNC,
                       "Failed to rewind frame stream by %g seconds at GPS time %" LAL_GPS_FORMAT, dt, LAL_GPS_PRINT( gps_tell ) );
    }

    // Free memory
    XLALFree( temp_SFT_filename );
    XLALFree( final_SFT_filename );

  }


  ////////// Free memory //////////

  XLALFree( SFT_comment );

  XLALDestroyCache( framecache );
  XLALFrStreamClose( framestream );

  XLALDestroyREAL8TimeSeries( SFT_time_series );
  XLALDestroyREAL8Window( SFT_window );
  XLALDestroyCOMPLEX16Vector( SFT_fft_data );
  XLALDestroyREAL8FFTPlan( SFT_fft_plan );
  XLALDestroySFT( SFT );

  XLALDestroyUserVars();

  LALCheckMemoryLeaks();

  return 0;
}
