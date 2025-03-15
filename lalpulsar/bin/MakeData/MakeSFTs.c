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
#include <stdbool.h>

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

/* Prototypes */

int main( int argc, char *argv[] )
{

  // Set help information
  lalUserVarHelpBrief = "generate SFTs";

  ////////// Parse user input //////////

  // default output
  LALStringVector *default_sft_write_path = NULL;
  XLAL_CHECK_MAIN( ( default_sft_write_path = XLALCreateStringVector( ".", NULL ) ) != NULL, XLAL_EFUNC );

  // Initialise user input variables
  struct uvar_type {
    char *frame_cache;
    BOOLEAN frame_checksums;
    LALStringVector *channel_name;
    BOOLEAN allow_skipping;
    INT4 gps_start_time;
    INT4 gps_end_time;
    INT4 sft_duration;
    REAL8 overlap_fraction;
    REAL8 high_pass_freq;
    char *window_type;
    REAL8 window_param;
    REAL8 start_freq;
    REAL8 band;
    LALStringVector *sft_write_path;
    UINT4 observing_run;
    char *observing_kind;
    UINT4 observing_revision;
    char *misc_desc;
    char *comment_field;
  } uvar_struct = {
    .sft_duration = 1800,
    .overlap_fraction = 0,
    .high_pass_freq = 0,
    .window_param = 0,
    .start_freq = 48,
    .band = 2000,
    .sft_write_path = default_sft_write_path,
    .allow_skipping = false,
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
    frame_checksums, BOOLEAN, 'k', OPTIONAL,
    "Validate frame checksums. Default is to only validate checksums for " UVAR_STR( observing_kind ) "=RUN SFTs. "
  );
  XLALRegisterUvarMember(
    channel_name, STRINGVector, 'N', REQUIRED,
    "Name(s) of channel(s) to read within a frame. "
  );
  XLALRegisterUvarMember(
    gps_start_time, INT4, 's', REQUIRED,
    "GPS time to start generating SFTs. "
  );
  XLALRegisterUvarMember(
    gps_end_time, INT4, 'e', REQUIRED,
    "GPS time to end generating SFTs. "
  );
  XLALRegisterUvarMember(
    allow_skipping, BOOLEAN, 'x', OPTIONAL,
    "Channel is allowed to be skipped if it is not in frames or has too low sampling frequency. "
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
    "Example: use " UVAR_STR( overlap_fraction ) "=0.5 with " UVAR_STR( window_type ) "=hann windows. "
  );
  XLALRegisterUvarMember(
    high_pass_freq, REAL8, 'f', REQUIRED,
    "High pass filtering frequency in Hertz. "
  );
  XLALRegisterUvarMember(
    window_type, STRING, 'w', REQUIRED,
    "Window to apply to SFTs. See https://dcc.ligo.org/T040164/public for supported options.\n\n"
    "The standard window applied to LVK production SFTs is documented in `lalpulsar_MakeSFTDAG`. "
  );
  XLALRegisterUvarMember(
    window_param, REAL8, 'r', OPTIONAL,
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
    sft_write_path, STRINGVector, 'p', OPTIONAL,
    "Path(s) to write SFTs to, same order as channel name(s). "
  );
  XLALRegisterUvarMember(
    observing_run, UINT4, 'O', REQUIRED,
    "For >=1, SFTs will be named following the 'public' scheme from T040164-v2, and this number should correspond to the observing run the data comes from (without the leading 'O').\n\n"
    "If set to 0, 'private' SFTs will be generated. "
  );
  XLALRegisterUvarMember(
    observing_kind, STRING, 'K', OPTIONAL,
    "For public SFTs (" UVAR_STR( observing_run ) ">=1) only, the kind of SFTs being generated: 'RUN', 'AUX', 'SIM', or 'DEV'. "
  );
  XLALRegisterUvarMember(
    observing_revision, UINT4, 'R', OPTIONAL,
    "For public SFTs (" UVAR_STR( observing_run ) ">=1) only, the revision number of the SFT production. "
  );
  XLALRegisterUvarMember(
    misc_desc, STRING, 'X', OPTIONAL,
    "For private SFTs only, a miscellaneous description field that will be included in the filenames. "
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
    NULL, "sft_version", INT4, 'V', DEFUNCT,
    "No longer supported. "
  );
  XLALRegisterNamedUvar(
    NULL, "use_single", BOOLEAN, 'S', DEFUNCT,
    "No longer supported. "
  );
  XLALRegisterNamedUvar(
    NULL, "window_radius", REAL8, 0, DEFUNCT,
    "Use " UVAR_STR( window_param ) " instead. "
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
  } else {
    XLALUserVarCheck( &should_exit,
                      XLALCheckNamedWindow( uvar->window_type, XLALUserVarWasSet( &uvar->window_param ) ) == XLAL_SUCCESS,
                      "Invalid/inconsistent " UVAR_STR( window_type ) " and/or " UVAR_STR( window_param ) );
  }
  XLALUserVarCheck( &should_exit,
                    uvar->window_param >= 0,
                    UVAR_STR( window_param ) " must be positive" );
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
                    uvar->observing_run == 0
                    || ( XLALUserVarWasSet( &uvar->observing_kind )
                         && XLALUserVarWasSet( &uvar->observing_revision ) ),
                    "Must set " UVAR_STR( observing_kind ) " and " UVAR_STR( observing_revision ) " when using " UVAR_STR( observing_run ) ">0" );
  XLALUserVarCheck( &should_exit,
                    uvar->observing_run == 0
                    || ( XLALUserVarWasSet( &uvar->observing_kind )
                         && ( strcmp( uvar->observing_kind, "RUN" ) == 0
                              || strcmp( uvar->observing_kind, "AUX" ) == 0
                              || strcmp( uvar->observing_kind, "SIM" ) == 0
                              || strcmp( uvar->observing_kind, "DEV" ) == 0 ) ),
                    UVAR_STR( observing_kind ) " must be one of 'RUN', 'AUX', 'SIM', or 'DEV'" );
  XLALUserVarCheck( &should_exit,
                    uvar->observing_run == 0 || uvar->observing_revision > 0,
                    UVAR_STR( observing_revision ) " must be strictly positive" );
  XLALUserVarCheck( &should_exit,
                    uvar->observing_run == 0 || !XLALUserVarWasSet( &uvar->misc_desc ),
                    UVAR_STR( observing_run ) "=0 is mutually exclusive with " UVAR_STR( misc_desc ) );
  XLALUserVarCheck( &should_exit,
                    uvar->observing_run > 0
                    || !( XLALUserVarWasSet( &uvar->observing_kind )
                          || XLALUserVarWasSet( &uvar->observing_revision ) ),
                    "Setting " UVAR_STR( observing_kind ) " or " UVAR_STR( observing_revision ) " is not allowed with " UVAR_STR( observing_run ) "=0" );
  XLALUserVarCheck( &should_exit,
                    uvar->channel_name->length == uvar->sft_write_path->length
                    || uvar->sft_write_path->length == 1,
                    "Number of channels in " UVAR_STR( channel_name ) " must be the same as the number of output paths in " UVAR_STR( sft_write_path ) "or a single path" );
  XLALUserVarCheck( &should_exit,
                    uvar->observing_run > 0
                    || uvar->channel_name->length == 1,
                    "For private SFTs (" UVAR_STR( observing_run ) "=0), can only provide one channel at a time" );

  // Exit if required
  if ( should_exit ) {
    return EXIT_FAILURE;
  }
  LogPrintf( LOG_NORMAL, "Parsed user input successfully\n" );

  // Decide whether to validate frame checksums
  bool frame_checksums = false;
  if ( XLALUserVarWasSet( &uvar->frame_checksums ) ) {
    if ( uvar->frame_checksums ) {
      frame_checksums = true;
      LogPrintf( LOG_NORMAL, "Validating frame checksums as " UVAR_STR( frame_checksums ) "=true\n" );
    } else {
      frame_checksums = false;
      LogPrintf( LOG_NORMAL, "Not validating frame checksums as " UVAR_STR( frame_checksums ) "=false\n" );
    }
  } else if ( XLALUserVarWasSet( &uvar->observing_kind ) ) {
    if ( strcmp( uvar->observing_kind, "RUN" ) == 0 ) {
      frame_checksums = true;
      LogPrintf( LOG_NORMAL, "Validating frame checksums as " UVAR_STR( observing_kind ) "=%s\n", uvar->observing_kind );
    } else {
      frame_checksums = false;
      LogPrintf( LOG_NORMAL, "Not validating frame checksums as " UVAR_STR( observing_kind ) "=%s!=RUN\n", uvar->observing_kind );
    }
  } else {
    frame_checksums = false;
    LogPrintf( LOG_NORMAL, "Not validating frame checksums as neither " UVAR_STR( frame_checksums ) " nor " UVAR_STR( observing_kind ) " were given\n" );
  }

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
  LogPrintf( LOG_NORMAL, "Opened frame stream from cache '%s'\n", uvar->frame_cache );

  // Set frame stream mode
  {
    LALFrStreamMode mode = 0;
    mode |= LAL_FR_STREAM_VERBOSE_MODE;           /* Display warnings and info */
    if ( frame_checksums ) {
      mode |= LAL_FR_STREAM_CHECKSUM_MODE;        /* Ensure that file checksums are OK */
    }
    XLAL_CHECK_MAIN( XLALFrStreamSetMode( framestream, mode ) == XLAL_SUCCESS, XLAL_EFUNC,
                     "Failed to set frame stream mode to %i", mode );
  }

  // Create SFT type for writing output
  const UINT4 SFT_first_bin = lrint( uvar->start_freq * uvar->sft_duration );
  const UINT4 SFT_bins = lrint( uvar->band * uvar->sft_duration );
  SFTtype *SFT = XLALCreateSFT( SFT_bins );
  XLAL_CHECK_MAIN( SFT != NULL, XLAL_EFUNC,
                   "Failed to allocate SFT of %u bins", SFT_bins );

  ////////// Generate SFTs //////////

  // SFT time series window, allocated once first time series data is read
  REAL8Window *SFT_window = NULL;

  // SFT FFT data vector and plan, allocated once first time series data is read
  COMPLEX16Vector *SFT_fft_data = NULL;
  REAL8FFTPlan *SFT_fft_plan = NULL;

  // Read data from frame stream until SFT end time is reached
  INT4 SFT_epoch_sec = uvar->gps_start_time;
  INT4 num_SFTs_made = 0;
  while ( 1 ) {

    const LIGOTimeGPS SFT_epoch = { .gpsSeconds = SFT_epoch_sec, .gpsNanoSeconds = 0 };

    // Loop over channels for this SFT interval
    for ( UINT4 n = 0; n < uvar->channel_name->length; n++ ) {
      // Check the channel in the framestream. It should return a type if it exists.
      // Need to check both the beginning and end of the time interval.
      // By default (uvar->require_channel = true), if the channel doesn't return a type, then this program fails
      // Optionally, a user could specify that the channel is not required to be in the frame in which case no SFT is made
      {
        int errnum;
        LALTYPECODE XLAL_INIT_DECL( laltype );

        // Set the frame stream to the start of the SFT
        XLAL_CHECK_MAIN( XLALFrStreamSeek( framestream, &SFT_epoch ) == 0, XLAL_EFUNC );

        // First check the beginning
        XLAL_TRY( laltype = XLALFrStreamGetTimeSeriesType( uvar->channel_name->data[n], framestream ), errnum );
        if ( errnum != 0 && uvar->allow_skipping ) {
          LogPrintf( LOG_CRITICAL, "--require-channel=FALSE and %s is not in frames\n", uvar->channel_name->data[n] );
          LogPrintf( LOG_CRITICAL, "==> No SFTs will be made for channel %s\n", uvar->channel_name->data[n] );
          continue;
        } else if ( errnum != 0 ) {
          XLAL_ERROR_MAIN( errnum );
        }

        // Shift the frame stream by SFT length minus one microsecond so that we don't hit the end of the framestream
        XLAL_CHECK_MAIN( XLALFrStreamSeekO( framestream, uvar->sft_duration - 1e-6, SEEK_CUR ) == 0, XLAL_EFUNC );

        // Now check the end
        XLAL_TRY( laltype = XLALFrStreamGetTimeSeriesType( uvar->channel_name->data[n], framestream ), errnum );
        if ( errnum != 0 && uvar->allow_skipping ) {
          LogPrintf( LOG_CRITICAL, "--require-channel=FALSE and %s is not in frames\n", uvar->channel_name->data[n] );
          LogPrintf( LOG_CRITICAL, "==> No SFTs will be made for channel %s\n", uvar->channel_name->data[n] );
          continue;
        } else if ( errnum != 0 ) {
          XLAL_ERROR_MAIN( errnum );
        }

        // Return to the original position
        XLAL_CHECK_MAIN( XLALFrStreamSeek( framestream, &SFT_epoch ) == 0, XLAL_EFUNC );
      }

      // Try to read in time series data for the next SFT
      REAL8TimeSeries *SFT_time_series = NULL;
      {
        SFT_time_series = XLALFrStreamInputREAL8TimeSeries( framestream, uvar->channel_name->data[n], &SFT_epoch, uvar->sft_duration, 0 );
        XLAL_CHECK_MAIN( SFT_time_series != NULL, XLAL_EFUNC,
                         "XLALFrStreamInputREAL8TimeSeries() failed at GPS time %" LAL_INT4_FORMAT, SFT_epoch_sec );
      }

      // Create SFT time series window
      if ( SFT_window != NULL && SFT_window->data->length != SFT_time_series->data->length ) {
        XLALDestroyREAL8Window( SFT_window );
        SFT_window = NULL;
      }
      if ( SFT_window == NULL ) {
        SFT_window = XLALCreateNamedREAL8Window( uvar->window_type, uvar->window_param, SFT_time_series->data->length );
        XLAL_CHECK_MAIN( SFT_window != NULL, XLAL_EFUNC,
                         "Failed to allocate SFT time series window of %u elements at GPS time %" LAL_INT4_FORMAT, SFT_time_series->data->length, SFT_epoch_sec );
      }

      // Create SFT FFT data vector and plan
      if ( SFT_fft_data != NULL && SFT_fft_data->length != SFT_time_series->data->length / 2 + 1 ) {
        XLALDestroyCOMPLEX16Vector( SFT_fft_data );
        XLALDestroyREAL8FFTPlan( SFT_fft_plan );
        SFT_fft_data = NULL;
        SFT_fft_plan = NULL;
      }
      if ( SFT_fft_data == NULL ) {
        SFT_fft_data = XLALCreateCOMPLEX16Vector( SFT_time_series->data->length / 2 + 1 );
        XLAL_CHECK_MAIN( SFT_fft_data != NULL, XLAL_EFUNC,
                         "Failed to allocate SFT FFT data vector of %u elements at GPS time %" LAL_INT4_FORMAT, SFT_time_series->data->length / 2 + 1, SFT_epoch_sec );
        SFT_fft_plan = XLALCreateForwardREAL8FFTPlan( SFT_time_series->data->length, 0 );
        XLAL_CHECK_MAIN( SFT_fft_plan != NULL, XLAL_EFUNC,
                         "Failed to allocate SFT FFT plan at GPS time %" LAL_INT4_FORMAT, SFT_epoch_sec );

        // If the sampling rate is too low for the requested band, skip this channel if allowed
        if ( SFT_fft_data->length < SFT_bins && uvar->allow_skipping ) {
          LogPrintf( LOG_CRITICAL, "Sampling rate is too low for band requested, skipping %s\n", uvar->channel_name->data[n] );
          XLALDestroyREAL8TimeSeries( SFT_time_series );
          XLALDestroyREAL8Window( SFT_window );
          XLALDestroyCOMPLEX16Vector( SFT_fft_data );
          XLALDestroyREAL8FFTPlan( SFT_fft_plan );
          SFT_window = NULL;
          SFT_fft_data = NULL;
          SFT_fft_plan = NULL;
          continue;
        }

        // Exit with error if the data length is too short
        XLAL_CHECK_MAIN( SFT_fft_data->length >= SFT_bins, XLAL_ESIZE );
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
                         "Failed to apply %s filter to SFT time series data at GPS time %" LAL_INT4_FORMAT, filterpar.name, SFT_epoch_sec );
      }

      // Window SFT time series data
      XLAL_CHECK_MAIN( XLALUnitaryWindowREAL8Sequence( SFT_time_series->data, SFT_window ) != NULL, XLAL_EFUNC,
                       "Failed to apply window to SFT time series data at GPS time %" LAL_INT4_FORMAT, SFT_epoch_sec );

      // Fourier transform SFT time series data
      XLAL_CHECK_MAIN( XLALREAL8ForwardFFT( SFT_fft_data, SFT_time_series->data, SFT_fft_plan ) == XLAL_SUCCESS, XLAL_EFUNC,
                       "Failed to Fourier transform SFT time series data at GPS time %" LAL_INT4_FORMAT, SFT_epoch_sec );

      // Initialise SFT type
      SFT->name[0] = uvar->channel_name->data[n][0];
      SFT->name[1] = uvar->channel_name->data[n][1] == '0' ? '1' : uvar->channel_name->data[n][1];
      SFT->name[2] = 0;
      SFT->epoch = SFT_time_series->epoch;
      SFT->f0 = uvar->start_freq;
      SFT->deltaF = 1.0 / ( ( REAL8 ) uvar->sft_duration );

      // Copy and normalise SFT frequency series data into SFT
      const REAL8 SFT_normalisation = SFT_time_series->deltaT;
      for ( UINT4 k = 0; k < SFT_bins; ++k ) {
        const REAL8 SFT_fft_re = creal( SFT_fft_data->data[k + SFT_first_bin] );
        const REAL8 SFT_fft_im = cimag( SFT_fft_data->data[k + SFT_first_bin] );
        SFT->data->data[k] = crectf( SFT_normalisation * SFT_fft_re, SFT_normalisation * SFT_fft_im );
      }

      // Build SFT filename spec
      SFTFilenameSpec XLAL_INIT_DECL( spec );
      if ( uvar->sft_write_path->length > 1 ) {
        XLAL_CHECK_MAIN( XLALFillSFTFilenameSpecStrings( &spec, uvar->sft_write_path->data[n], "sft_TO_BE_VALIDATED", NULL, uvar->window_type, uvar->misc_desc, uvar->observing_kind, uvar->channel_name->data[n] ) == XLAL_SUCCESS, XLAL_EFUNC );
      } else {
        XLAL_CHECK_MAIN( XLALFillSFTFilenameSpecStrings( &spec, uvar->sft_write_path->data[0], "sft_TO_BE_VALIDATED", NULL, uvar->window_type, uvar->misc_desc, uvar->observing_kind, uvar->channel_name->data[n] ) == XLAL_SUCCESS, XLAL_EFUNC );
      }
      spec.window_param = uvar->window_param;
      spec.pubObsRun = uvar->observing_run;
      spec.pubRevision = uvar->observing_revision;

      // Write out SFT and retrieve filename
      char *temp_SFT_filename = NULL;
      {
        int retn = XLALWriteSFT2StandardFile( SFT, &spec, SFT_comment );
        temp_SFT_filename = XLALBuildSFTFilenameFromSpec( &spec );
        XLAL_CHECK_MAIN( retn == XLAL_SUCCESS, XLAL_EFUNC,
                         "Failed to write SFT '%s' at GPS time %" LAL_INT4_FORMAT, temp_SFT_filename ? temp_SFT_filename : "<unknown>", SFT_epoch_sec );
        XLAL_CHECK_MAIN( temp_SFT_filename != NULL, XLAL_EFUNC );
      }

      // Validate SFT
      XLAL_CHECK_MAIN( XLALCheckSFTFileIsValid( temp_SFT_filename ) == XLAL_SUCCESS, XLAL_EFUNC,
                       "Failed to validate SFT '%s' at GPS time %" LAL_INT4_FORMAT, temp_SFT_filename, SFT_epoch_sec );

      // Move SFT to final filename
      char *final_SFT_filename = XLALStringDuplicate( temp_SFT_filename );
      char *const extn = strrchr( final_SFT_filename, '.' );
      XLAL_CHECK_MAIN( extn != NULL, XLAL_ESYS );
      XLAL_CHECK_MAIN( strlen( extn ) >= 4, XLAL_EFAILED );
      strcpy( extn, ".sft" );
      XLAL_CHECK_MAIN( rename( temp_SFT_filename, final_SFT_filename ) == 0, XLAL_ESYS,
                       "Failed to rename '%s' to '%s': %s", temp_SFT_filename, final_SFT_filename, strerror( errno ) );
      LogPrintf( LOG_DEBUG, "Wrote SFT '%s'\n", final_SFT_filename );

      // Increment number of SFTs made
      if ( ++num_SFTs_made == 1 ) {
        LogPrintf( LOG_NORMAL, "Generated first SFT at GPS time %" LAL_INT4_FORMAT "\n", SFT_epoch_sec );
      }

      // Free memory
      XLALDestroyREAL8TimeSeries( SFT_time_series );
      XLALFree( temp_SFT_filename );
      XLALFree( final_SFT_filename );

      // Set LALFrStream to correct location, if needed
      if ( n < uvar->channel_name->length - 1 ) {
        XLAL_CHECK_MAIN( XLALFrStreamSeek( framestream, &SFT_epoch ) == 0, XLAL_EFUNC );
      }

    }

    // Advance to next SFT, overlapping if requested
    {
      const INT4 last_SFT_epoch_sec = SFT_epoch_sec;
      SFT_epoch_sec += ( 1.0 - uvar->overlap_fraction ) * ( ( REAL8 ) uvar->sft_duration );
      if ( SFT_epoch_sec + uvar->sft_duration > uvar->gps_end_time ) {
        LogPrintf( LOG_NORMAL, "Generated last SFT at GPS time %" LAL_INT4_FORMAT "\n", last_SFT_epoch_sec );
        LogPrintf( LOG_NORMAL, "Generated %" LAL_INT4_FORMAT " SFTs\n", num_SFTs_made );
        break;
      }
    }
  }

  ////////// Free memory //////////

  XLALFree( SFT_comment );

  XLALDestroyCache( framecache );
  XLALFrStreamClose( framestream );

  XLALDestroyREAL8Window( SFT_window );
  XLALDestroyCOMPLEX16Vector( SFT_fft_data );
  XLALDestroyREAL8FFTPlan( SFT_fft_plan );
  XLALDestroySFT( SFT );

  XLALDestroyUserVars();

  LALCheckMemoryLeaks();

  return 0;
}
