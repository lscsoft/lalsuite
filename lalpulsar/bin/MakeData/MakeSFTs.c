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

#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/UserInput.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/Window.h>
#include <lal/RealFFT.h>
#include <lal/SFTfileIO.h>
#include <lal/LALFrStream.h>
#include <lal/LALVCSInfo.h>
#include <lal/LALPulsarVCSInfo.h>

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { REPORTSTATUS(pstat); return 100; } else ((void)0)

/***************************************************************************/

/* GLOBAL VARIABLES */
extern REAL8 winFncRMS;
extern REAL8TimeSeries dataDouble;
extern REAL4TimeSeries dataSingle;

static LALStatus status;
LALCache *framecache;         /* frame reading variables */
LALFrStream *framestream = NULL;

INT4 SegmentDuration;
LIGOTimeGPS gpsepoch;

REAL8FFTPlan *fftPlanDouble;           /* fft plan and data container, double precision case */
COMPLEX16Vector *fftDataDouble = NULL;

/***************************************************************************/

/* FUNCTION PROTOTYPES */
/* Windows data */
int WindowData( REAL8 r );
int WindowDataTukey2( void );
int WindowDataHann( void );

/* prototypes */
void getSFTDescField( CHAR *sftDescField, CHAR *numSFTs, CHAR *ifo, CHAR *stringT, CHAR *typeMisc );
void mkSFTFilename( CHAR *sftFilename, CHAR *site, CHAR *numSFTs, CHAR *ifo, CHAR *stringT, CHAR *typeMisc, CHAR *gpstime );
void mvFilenames( CHAR *filename1, CHAR *filename2 );


/* -------------------- function definitions -------------------- */


/* 12/27/05 gam; function to create sftDescField for output directory or filename. */
void getSFTDescField( CHAR *sftDescField, CHAR *numSFTs, CHAR *ifo, CHAR *stringT, CHAR *typeMisc )
{
  strcpy( sftDescField, numSFTs );
  strcat( sftDescField, "_" );
  strcat( sftDescField, ifo );
  strcat( sftDescField, "_" );
  strcat( sftDescField, stringT );
  strcat( sftDescField, "SFT" );
  if ( typeMisc != NULL ) {
    strcat( sftDescField, "_" );
    strcat( sftDescField, typeMisc );
  }
}

/* 12/27/05 gam; make SFT file name according to LIGO T040164-01 specification */
void mkSFTFilename( CHAR *sftFilename, CHAR *site, CHAR *numSFTs, CHAR *ifo, CHAR *stringT, CHAR *typeMisc, CHAR *gpstime )
{
  CHAR sftDescField[256];
  strcpy( sftFilename, site );
  strcat( sftFilename, "-" );
  getSFTDescField( sftDescField, numSFTs, ifo, stringT, typeMisc );
  strcat( sftFilename, sftDescField );
  strcat( sftFilename, "-" );
  strcat( sftFilename, gpstime );
  strcat( sftFilename, "-" );
  strcat( sftFilename, stringT );
  strcat( sftFilename, ".sft" );
}

/* 01/09/06 gam; move filename1 to filename2 */
void mvFilenames( CHAR *filename1, CHAR *filename2 )
{
  CHAR mvFilenamesCommand[512 + 4];
  sprintf( mvFilenamesCommand, "mv %s %s", filename1, filename2 );
  if ( system( mvFilenamesCommand ) ) {
    XLALPrintError( "system() returned non-zero status\n" );
  }
}



/************************************* MAIN PROGRAM *************************************/

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
    INT4 window_type;
    REAL8 window_radius;
    REAL8 start_freq;
    REAL8 band;
    char *sft_write_path;
    char *comment_field;
  } uvar_struct = {
    .sft_duration = 1800,
    .overlap_fraction = 0,
    .high_pass_freq = 0,
    .window_type = 1,
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
    "Example: use " UVAR_STR( overlap_fraction ) "=0.5 with " UVAR_STR( window_type ) "=3 Hann windows. " );
  XLALRegisterUvarMember(
    high_pass_freq, REAL8, 'f', REQUIRED,
    "High pass filtering frequency in Hertz. " );
  XLALRegisterUvarMember(
    window_type, INT4, 'w', OPTIONAL,
    "Window to apply to SFTs. "
    "0 = no window, 1 = Matlab style Tukey window, 2 = make_sfts.c Tukey window, 3 = Hann window. "
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
                    0 <= uvar->overlap_fraction && uvar->overlap_fraction < 1,
                    UVAR_STR( overlap_fraction ) " must be within range [0.0, 1.0)" );
  XLALUserVarCheck( &should_exit,
                    uvar->high_pass_freq >= 0,
                    UVAR_STR( high_pass_freq ) " must be positive" );
  XLALUserVarCheck( &should_exit,
                    0 <= uvar->window_type && uvar->window_type <= 3,
                    UVAR_STR( >window_type ) " must be 0, 1, 2, or 3" );
  XLALUserVarCheck( &should_exit,
                    uvar->start_freq >= 0,
                    UVAR_STR( start_freq ) " must be positive" );
  XLALUserVarCheck( &should_exit,
                    uvar->band > 0,
                    UVAR_STR( band ) " must be strictly positive" );

  // Exit if required
  if ( should_exit ) {
    return EXIT_FAILURE;
  }

  ////////// Set up //////////

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

  SegmentDuration = uvar->gps_end_time - uvar->gps_start_time ;

  /* create Frame cache, open frame stream and delete frame cache */
  framecache = XLALCacheImport( uvar->frame_cache );
  LALFrCacheOpen( &status, &framestream, framecache );
  TESTSTATUS( &status );
  XLALDestroyCache( framecache );

  if ( SegmentDuration < uvar->sft_duration ) {
    fprintf( stderr, "Cannot fit an SFT of duration %d between %d and %d\n",
             uvar->sft_duration, uvar->gps_start_time, uvar->gps_end_time );
    return 0;;
  }

  gpsepoch.gpsSeconds = uvar->gps_start_time;
  gpsepoch.gpsNanoSeconds = 0;

  /* Allocates space for data */
  {
    static FrChanIn chanin;

    chanin.name  = uvar->channel_name;

    /* These calls just return deltaT for the channel */
    chanin.type  = ProcDataChannel;
    /* Get channel time step size by calling LALFrGetREAL8TimeSeries */
    LALFrSeek( &status, &gpsepoch, framestream );
    TESTSTATUS( &status );
    LALFrGetREAL8TimeSeries( &status, &dataDouble, &chanin, framestream );
    TESTSTATUS( &status );
    dataSingle.deltaT = dataDouble.deltaT;

    LALDCreateVector( &status, &dataDouble.data, ( UINT4 )( uvar->sft_duration / dataDouble.deltaT + 0.5 ) );
    TESTSTATUS( &status );

    fftPlanDouble = XLALCreateForwardREAL8FFTPlan( dataDouble.data->length, 0 );
    XLAL_CHECK( fftPlanDouble != NULL, XLAL_EFUNC );

  }

  while ( gpsepoch.gpsSeconds + uvar->sft_duration <= uvar->gps_end_time ) {

    /* Reads T seconds of data */
    {
      static FrChanIn chanin;
      chanin.name  = uvar->channel_name;
      LALFrSeek( &status, &gpsepoch, framestream );
      TESTSTATUS( &status );

      chanin.type  = ProcDataChannel;
      LALFrGetREAL8TimeSeries( &status, &dataDouble, &chanin, framestream );
      TESTSTATUS( &status );
    }

    /* High-pass data with Butterworth filter */
    {
      PassBandParamStruc filterpar;
      char tmpname[] = "Butterworth High Pass";

      filterpar.name  = tmpname;
      filterpar.nMax  = 10;
      filterpar.f2    = uvar->high_pass_freq;
      filterpar.a2    = 0.5;
      filterpar.f1    = -1.0;
      filterpar.a1    = -1.0;

      if ( uvar->high_pass_freq > 0.0 ) {

        /* High pass the time series */
        LALButterworthREAL8TimeSeries( &status, &dataDouble, &filterpar );
        TESTSTATUS( &status );

      }

    }

    /* Window data; 12/28/05 gam; add options */
    if ( uvar->window_type == 1 ) {
      if ( WindowData( uvar->window_radius ) ) {
        return 5;  /* uvar->window_type==1 is the default */
      }
    } else if ( uvar->window_type == 2 ) {
      if ( WindowDataTukey2( ) ) {
        return 5;
      }
    } else if ( uvar->window_type == 3 ) {
      if ( WindowDataHann( ) ) {
        return 5;
      }
    } else {
      /* Continue with no windowing; parsing of command line args makes sure options are one of the above or 0 for now windowing. */
    }

    /* create an SFT */
    /* 11/02/05 gam; allocate container for SFT data */
    LALZCreateVector( &status, &fftDataDouble, dataDouble.data->length / 2 + 1 );
    TESTSTATUS( &status );

    /* compute sft */
    XLAL_CHECK( XLALREAL8ForwardFFT( fftDataDouble, dataDouble.data, fftPlanDouble ) == XLAL_SUCCESS, XLAL_EFUNC );

    /* write out sft */
    {
      char sftname[256];
      char sftFilename[256];
      char sftnameFinal[256]; /* 01/09/06 gam */
      char numSFTs[2]; /* 12/27/05 gam */
      char site[2];    /* 12/27/05 gam */
      char ifo[3];     /* 12/27/05 gam; allow 3rd charactor for null termination */
      int firstbin = ( INT4 )( uvar->start_freq * uvar->sft_duration + 0.5 ), k;
      char gpstime[11]; /* 12/27/05 gam; allow for 10 digit GPS times and null termination */
      SFTtype *oneSFT = NULL;
      INT4 nBins = ( INT4 )( uvar->band * uvar->sft_duration + 0.5 );
      REAL8 doubleDeltaT = 0.0; /* 01/05/06 gam */

      /* 12/27/05 gam; set up the number of SFTs, site, and ifo as null terminated strings */
      numSFTs[0] = '1';
      numSFTs[1] = '\0'; /* null terminate */
      strncpy( site, uvar->channel_name, 1 );
      site[1] = '\0'; /* null terminate */
      strncpy( ifo, uvar->channel_name, 2 );
      ifo[2] = '\0'; /* null terminate */
      sprintf( gpstime, "%09d", gpsepoch.gpsSeconds );

      strcpy( sftname, uvar->sft_write_path );

      strcat( sftname, "/" );
      char stringT[256];
      snprintf( stringT, sizeof(stringT), "%d", uvar->sft_duration );
      mkSFTFilename( sftFilename, site, numSFTs, ifo, stringT, NULL, gpstime );
      /* 01/09/06 gam; sftname will be temporary; will move to sftnameFinal. */
      /* set up sftnameFinal with usual SFT name */
      strcpy( sftnameFinal, sftname );
      strcat( sftnameFinal, sftFilename );
      /* sftname begins with . and ends in .tmp */
      strcat( sftname, "." );
      strcat( sftname, sftFilename );
      strcat( sftname, ".tmp" );

      /* make container to store the SFT data */
      XLAL_CHECK( ( oneSFT = XLALCreateSFT( ( ( UINT4 )nBins ) ) ) != NULL, XLAL_EFUNC );

      /* copy the data to oneSFT */
      strcpy( oneSFT->name, ifo );
      oneSFT->epoch.gpsSeconds = gpsepoch.gpsSeconds;
      oneSFT->epoch.gpsNanoSeconds = gpsepoch.gpsNanoSeconds;
      oneSFT->f0 = uvar->start_freq;
      oneSFT->deltaF = 1.0 / ( ( REAL8 )uvar->sft_duration );
      oneSFT->data->length = nBins;

      doubleDeltaT = ( REAL8 )( dataDouble.deltaT / winFncRMS ); /* include 1 over window function RMS */
      for ( k = 0; k < nBins; k++ ) {
        oneSFT->data->data[k] = crectf( doubleDeltaT * creal( fftDataDouble->data[k + firstbin] ), doubleDeltaT * cimag( fftDataDouble->data[k + firstbin] ) );
        /* 06/26/07 gam; use finite to check that data does not contains a non-FINITE (+/- Inf, NaN) values */
        if ( !isfinite( crealf( oneSFT->data->data[k] ) ) || !isfinite( cimagf( oneSFT->data->data[k] ) ) ) {
          fprintf( stderr, "Infinite or NaN data at freq bin %d.\n", k );
          return 7;
        }
      }
      LALZDestroyVector( &status, &fftDataDouble );
      TESTSTATUS( &status );

      /* write the SFT */
      XLAL_CHECK( XLALWriteSFT2NamedFile(oneSFT, sftname, "unknown" /* FIXME */, 0, SFT_comment) == XLAL_SUCCESS, XLAL_EFUNC );

      /* 01/09/06 gam; sftname is temporary; move to sftnameFinal. */
      mvFilenames( sftname, sftnameFinal );

      XLALDestroySFT( oneSFT );
    }

    gpsepoch.gpsSeconds = gpsepoch.gpsSeconds + ( INT4 )( ( 1.0 - uvar->overlap_fraction ) * ( ( REAL8 )uvar->sft_duration ) );
    gpsepoch.gpsNanoSeconds = 0;
  }

  LALFrClose( &status, &framestream );
  TESTSTATUS( &status );

  LALDDestroyVector( &status, &dataDouble.data );
  TESTSTATUS( &status );
  XLALDestroyREAL8FFTPlan( fftPlanDouble );

  ////////// Free memory //////////

  XLALFree( SFT_comment );

  XLALDestroyUserVars();

  LALCheckMemoryLeaks();

  return 0;
}
