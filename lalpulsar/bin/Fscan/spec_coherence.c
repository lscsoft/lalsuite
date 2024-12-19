/*
*  Copyright (C) 2021 Evan Goetz
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

/**
* \file
* \ingroup lalpulsar_bin_fscan
*/

#include <libgen.h>
#include <unistd.h>
#include <lal/SFTfileIO.h>
#include <lal/LALStdio.h>
#include <lal/LogPrintf.h>
#include <lal/UserInput.h>
#include <lal/LALPulsarVCSInfo.h>

#include "fscanutils.h"

int main( int argc, char **argv )
{
  FILE *COHOUT = NULL;
  int fopenerr = 0;

  SFTCatalog *catalog_a = NULL, *catalog_b = NULL;
  SFTVector *sft_vect_a = NULL, *sft_vect_b = NULL;
  SFTConstraints XLAL_INIT_DECL( constraints );
  LIGOTimeGPS startTime, endTime;
  REAL8Vector *psd_a = NULL, *psd_b = NULL;
  COMPLEX16Vector *coh = NULL;
  REAL8 f0 = 0, deltaF = 0;
  CHAR outbase[256], outfile0[512];

  CHAR *SFTpattA = NULL, *SFTpattB = NULL, *outputDir = NULL, *outputBname = NULL;
  INT4 startGPS = 0, endGPS = 0;
  REAL8 f_min = 0.0, f_max = 0.0, timebaseline = 0;

  /* Default for output directory */
  XLAL_CHECK_MAIN( ( outputDir = XLALStringDuplicate( "." ) ) != NULL, XLAL_EFUNC );

  BOOLEAN allow_skipping = 0;

  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &SFTpattA,      "ChASFTs",         STRING, 'p', REQUIRED, "SFT location/pattern. Possibilities are:\n"
                                          " - '<SFT file>;<SFT file>;...', where <SFT file> may contain wildcards\n - 'list:<file containing list of SFT files>'" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &SFTpattB,      "ChBSFTs",         STRING, 'q', REQUIRED, "SFT location/pattern. Possibilities are:\n"
                                          " - '<SFT file>;<SFT file>;...', where <SFT file> may contain wildcards\n - 'list:<file containing list of SFT files>'" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &startGPS,     "startGPS",     INT4,   's', REQUIRED, "Starting GPS time (SFT timestamps must be >= this)" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &endGPS,       "endGPS",       INT4,   'e', REQUIRED, "Ending GPS time (SFT timestamps must be < this)" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &f_min,        "fMin",         REAL8,  'f', REQUIRED, "Minimum frequency" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &f_max,        "fMax",         REAL8,  'F', REQUIRED, "Maximum frequency" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &outputDir,    "outputDir",    STRING, 'd', OPTIONAL, "Output directory for data files" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &outputBname,  "outputBname",  STRING, 'o', OPTIONAL, "Base name of output files" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &timebaseline, "timeBaseline", REAL8,  't', REQUIRED, "The time baseline of sfts" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &allow_skipping, "allow_skipping", BOOLEAN, 'x', OPTIONAL, "Allow to exit without an error if no SFTs are found" ) == XLAL_SUCCESS, XLAL_EFUNC );

  BOOLEAN should_exit = 0;
  XLAL_CHECK_MAIN( XLALUserVarReadAllInput( &should_exit, argc, argv, lalPulsarVCSInfoList ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( should_exit ) {
    return ( 1 );
  }

  printf( "Starting spec_coherence...\n" );

  /* Provide the constraints to the catalog */
  startTime.gpsSeconds = startGPS;
  startTime.gpsNanoSeconds = 0;
  constraints.minStartTime = &startTime;
  endTime.gpsSeconds = endGPS;
  endTime.gpsNanoSeconds = 0;
  constraints.maxStartTime = &endTime;

  int errnumA, errnumB;
  printf( "Calling XLALSFTdataFind with SFTpattA=%s\n", SFTpattA );
  XLAL_TRY( catalog_a = XLALSFTdataFind( SFTpattA, &constraints ), errnumA );
  printf( "Calling XLALSFTdataFind with SFTpattB=%s\n", SFTpattB );
  XLAL_TRY( catalog_b = XLALSFTdataFind( SFTpattB, &constraints ), errnumB );

  /* Ensure that some SFTs were found given the start and end time and IFO constraints unless --allow_skipping is true */
  if ( errnumA != 0 || errnumB != 0 ) {
    if ( errnumA != 0 && allow_skipping ) {
      LogPrintf( LOG_CRITICAL, "Channel A %s found no SFTs, exiting with code %d due to --allow_skipping=true\n", SFTpattA, XLAL_EUSR0 );
      if ( catalog_a != NULL ) {
        XLALDestroySFTCatalog( catalog_a );
      }
      if ( catalog_b != NULL ) {
        XLALDestroySFTCatalog( catalog_b );
      }
      XLALDestroyUserVars();
      exit( XLAL_EUSR0 );
    } else if ( errnumA != 0 ) {
      XLAL_ERROR_MAIN( errnumA );
    }
    if ( errnumB != 0 && allow_skipping ) {
      LogPrintf( LOG_CRITICAL, "Channel B %s found no SFTs, exiting with code %d due to --allow_skipping=true\n", SFTpattB, XLAL_EUSR0 );
      if ( catalog_a != NULL ) {
        XLALDestroySFTCatalog( catalog_a );
      }
      if ( catalog_b != NULL ) {
        XLALDestroySFTCatalog( catalog_b );
      }
      XLALDestroyUserVars();
      exit( XLAL_EUSR0 );
    } else if ( errnumB != 0 ) {
      XLAL_ERROR_MAIN( errnumB );
    }
  }

  XLAL_CHECK_MAIN( catalog_a->length > 0, XLAL_EFAILED, "No SFTs found for Ch A, please examine start time, end time, frequency range, etc." );
  XLAL_CHECK_MAIN( catalog_b->length > 0, XLAL_EFAILED, "No SFTs found for Ch B, please examine start time, end time, frequency range, etc." );

  LogPrintf( LOG_NORMAL, "Channel A %s has length of %u SFT files\n", SFTpattA, catalog_a->length );
  LogPrintf( LOG_NORMAL, "Channel B %s has length of %u SFT files\n", SFTpattB, catalog_b->length );

  if ( XLALUserVarWasSet( &outputBname ) ) {
    snprintf( outbase, sizeof( outbase ), "%s/%s", outputDir, outputBname );
  } else {
    snprintf( outbase, sizeof( outbase ), "%s/spec_%.2f_%.2f_%d_%d_coh", outputDir, f_min, f_max, startTime.gpsSeconds, endTime.gpsSeconds );
  }

  snprintf( outfile0, sizeof( outfile0 ), "%s.txt", outbase );

  COHOUT = fopen( outfile0, "w" );
  fopenerr = errno;
  XLAL_CHECK_MAIN( COHOUT != NULL, XLAL_EIO, "Failed to open '%s' for writing: %s", outfile0, strerror( fopenerr ) );

  UINT4 nAve = 0;

  printf( "Looping over SFTs to compute coherence\n" );
  for ( UINT4 j = 0; j < catalog_a->length; j++ ) {
    /* Extract one SFT at a time from the catalog */
    fprintf( stderr, "Extracting SFT %d...\n", j );
    XLAL_CHECK_MAIN( ( sft_vect_a = extract_one_sft( catalog_a, catalog_a->data[j].header.epoch, f_min, f_max ) ) != NULL, XLAL_EFUNC );

    /* If no SFT from the B list was found, then just continue with the next SFT in the A list */
    XLAL_TRY( sft_vect_b = extract_one_sft( catalog_b, catalog_a->data[j].header.epoch, f_min, f_max ), errnumB );
    if ( errnumB != 0 ) {
      LogPrintf( LOG_CRITICAL, "Failed to find B channel SFT at time %d, [%.9f, %.9f) Hz\n", catalog_a->data[j].header.epoch.gpsSeconds, f_min, f_max );
      continue;
    }

    /* Check time baseline of the SFTs to confirm they match */
    XLAL_CHECK_MAIN( sft_vect_a->data[0].deltaF * timebaseline == 1.0, XLAL_EINVAL, "Time baseline of SFTs and the request do not match" );
    XLAL_CHECK_MAIN( sft_vect_b->data[0].deltaF * timebaseline == 1.0, XLAL_EINVAL, "Time baseline of SFTs and the request do not match" );

    /* For the first time through the loop, we allocate some vectors */
    if ( nAve == 0 ) {
      UINT4 numBins = sft_vect_a->data->data->length;
      f0 = sft_vect_a->data->f0;
      deltaF = sft_vect_a->data->deltaF;

      XLAL_CHECK_MAIN( ( coh = XLALCreateCOMPLEX16Vector( numBins ) ) != NULL, XLAL_EFUNC );
      XLAL_CHECK_MAIN( ( psd_a = XLALCreateREAL8Vector( numBins ) ) != NULL, XLAL_EFUNC );
      XLAL_CHECK_MAIN( ( psd_b = XLALCreateREAL8Vector( numBins ) ) != NULL, XLAL_EFUNC );
    }

    /* Loop over the SFT bins computing cross spectrum (AB) and auto spectrum (AA and BB) */
    if ( nAve == 0 ) {
      for ( UINT4 i = 0; i < sft_vect_a->data->data->length; i++ ) {
        coh->data[i] = sft_vect_a->data[0].data->data[i] * conj( sft_vect_b->data[0].data->data[i] );
        psd_a->data[i] = sft_vect_a->data[0].data->data[i] * conj( sft_vect_a->data[0].data->data[i] );
        psd_b->data[i] = sft_vect_b->data[0].data->data[i] * conj( sft_vect_b->data[0].data->data[i] );
      }
    } else {
      for ( UINT4 i = 0; i < sft_vect_a->data->data->length; i++ ) {
        coh->data[i] += sft_vect_a->data[0].data->data[i] * conj( sft_vect_b->data[0].data->data[i] );
        psd_a->data[i] += sft_vect_a->data[0].data->data[i] * conj( sft_vect_a->data[0].data->data[i] );
        psd_b->data[i] += sft_vect_b->data[0].data->data[i] * conj( sft_vect_b->data[0].data->data[i] );
      }
    }
    /* Destroys current SFT Vectors */
    XLALDestroySFTVector( sft_vect_a );
    XLALDestroySFTVector( sft_vect_b );
    sft_vect_a = sft_vect_b = NULL;

    nAve++;
  }

  XLAL_CHECK_MAIN( nAve > 0, XLAL_EFAILED, "No SFTs were found to be matching" );

  /* compute the final coherence (|AB|**2 / (AA * BB)) and print to file */
  for ( UINT4 i = 0; i < coh->length; i++ ) {
    REAL8 f = f0 + ( ( REAL4 )i ) * deltaF;
    REAL8 COH = coh->data[i] * conj( coh->data[i] ) / ( psd_a->data[i] * psd_b->data[i] );
    fprintf( COHOUT, "%16.8f %g\n", f, COH );
  }

  fprintf( stderr, "Destroying Variables\n" );
  XLALDestroySFTCatalog( catalog_a );
  XLALDestroySFTCatalog( catalog_b );

  XLALDestroyCOMPLEX16Vector( coh );
  XLALDestroyREAL8Vector( psd_a );
  XLALDestroyREAL8Vector( psd_b );

  fprintf( stderr, "Closing Files\n" );
  fclose( COHOUT );

  XLALDestroyUserVars();
  fprintf( stderr, "Done Destroying Variables\n" );
  fprintf( stderr, "end of spec_coherence\n" );
  fprintf( stderr, "Spec_coherence_done!\n" );

  return ( 0 );

}
/* END main */
