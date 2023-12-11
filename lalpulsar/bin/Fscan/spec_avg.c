/*
*  Copyright (C) 2007 Gregory Mendell
*                2021, 2023 Evan Goetz
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
 * \ingroup lalpulsar_bin_SFTTools
 */

#include <lal/UserInput.h>
#include <lal/FrequencySeries.h>
#include <lal/SFTfileIO.h>
#include <lal/NormalizeSFTRngMed.h>
#include <lal/LALPulsarVCSInfo.h>

#include "fscanutils.h"


int main( int argc, char **argv )
{
  FILE *fp  = NULL, *fp2 = NULL, *fp3 = NULL, *fp4 = NULL;
  int fopenerr = 0;

  SFTCatalog *catalog = NULL;
  SFTVector *sft_vect = NULL;
  UINT4 NumBinsAvg = 1;
  SFTConstraints XLAL_INIT_DECL( constraints );
  LIGOTimeGPS XLAL_INIT_DECL( startTime );
  LIGOTimeGPS XLAL_INIT_DECL( endTime );
  REAL4Vector *timeavg = NULL;
  CHAR outbase[256], outfile[512], outfile2[512], outfile3[512], outfile4[512];

  CHAR *SFTpatt = NULL, *IFO = NULL, *outputDir = NULL, *outputBname = NULL;
  INT4 startGPS = 0, endGPS = 0;
  REAL8 f_min = 0.0, f_max = 0.0, freqres = 0.1, subband = 100.0, timebaseline = 0;
  INT4 blocksRngMed = 101, cur_epoch = 0;

  /* these varibales are for converting GPS seconds into UTC time and date*/
  struct tm date;
  INT4Vector *timestamps = NULL;

  /* Default for output directory */
  XLAL_CHECK_MAIN( ( outputDir = XLALStringDuplicate( "." ) ) != NULL, XLAL_EFUNC );

  /*========================================================================================================================*/


  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &SFTpatt,      "SFTs",         STRING, 'p', REQUIRED, "SFT location/pattern. Possibilities are:\n"
                                          " - '<SFT file>;<SFT file>;...', where <SFT file> may contain wildcards\n - 'list:<file containing list of SFT files>'" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &IFO,          "IFO",          STRING, 'I', REQUIRED, "Detector" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &startGPS,     "startGPS",     INT4,   's', REQUIRED, "Starting GPS time (SFT timestamps must be >= this)" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &endGPS,       "endGPS",       INT4,   'e', REQUIRED, "Ending GPS time (SFT timestamps must be < this)" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &f_min,        "fMin",         REAL8,  'f', REQUIRED, "Minimum frequency in Hz" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &f_max,        "fMax",         REAL8,  'F', REQUIRED, "Maximum frequency in Hz" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &blocksRngMed, "blocksRngMed", INT4,   'w', OPTIONAL, "Running Median window size" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &outputDir,    "outputDir",    STRING, 'd', OPTIONAL, "Output directory for data files" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &outputBname,  "outputBname",  STRING, 'o', OPTIONAL, "Base name of output files" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &freqres,      "freqRes",      REAL8,  'r', OPTIONAL, "Spectrogram freq resolution in Hz" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &subband,      "subband",      REAL8,  'b', OPTIONAL, "Subdivide the output normalized average spectra txt files into these subbands" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &timebaseline, "timeBaseline", REAL8,  't', REQUIRED, "The time baseline of sfts in seconds" ) == XLAL_SUCCESS, XLAL_EFUNC );

  BOOLEAN should_exit = 0;
  XLAL_CHECK_MAIN( XLALUserVarReadAllInput( &should_exit, argc, argv, lalPulsarVCSInfoList ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( should_exit ) {
    return ( 1 );
  }

  // Populate the startTime and endTime LIGOTimeGPS variables
  startTime.gpsSeconds = startGPS;
  endTime.gpsSeconds = endGPS;
  startTime.gpsNanoSeconds = endTime.gpsNanoSeconds = 0;

  // Populate the SFT catalog constraints
  constraints.minStartTime = &startTime;
  constraints.maxStartTime = &endTime;
  constraints.detector = IFO;

  // Load SFT catalog
  XLAL_CHECK_MAIN( ( catalog = XLALSFTdataFind( SFTpatt, &constraints ) ) != NULL, XLAL_EFUNC );

  // Ensure that some SFTs were found given the start and end time and IFO constraints
  XLAL_CHECK_MAIN( catalog->length > 0, XLAL_EFAILED, "No SFTs found, please examine start time, end time, frequency range, etc." );

  // If output base name was set by user, use that, otherwise use a default pattern:
  // spec_<f_min>_<f_max>_<detector>_<GPS-start>_<GPS-end> as the basename
  if ( XLALUserVarWasSet( &outputBname ) ) {
    snprintf( outbase, sizeof( outbase ), "%s/%s", outputDir, outputBname );
  } else {
    snprintf( outbase, sizeof( outbase ), "%s/spec_%.2f_%.2f_%s_%d_%d", outputDir, f_min, f_max, constraints.detector, startTime.gpsSeconds, endTime.gpsSeconds );
  }

  // Create filenames from the outbase value
  snprintf( outfile, sizeof( outfile ),  "%s_spectrogram.txt", outbase );
  snprintf( outfile2, sizeof( outfile2 ), "%s_timestamps.txt", outbase );
  snprintf( outfile4, sizeof( outfile4 ), "%s_date.txt", outbase );

  // Open the files for writing using the "w" mode, which overwrites files if they exist
  fp = fopen( outfile, "w" );
  fopenerr = errno;
  XLAL_CHECK_MAIN( fp != NULL, XLAL_EIO, "Failed to open '%s' for writing: %s", outfile, strerror( fopenerr ) );
  fp2 = fopen( outfile2, "w" );
  fopenerr = errno;
  XLAL_CHECK_MAIN( fp2 != NULL, XLAL_EIO, "Failed to open '%s' for writing: %s", outfile2, strerror( fopenerr ) );
  fp4 = fopen( outfile4, "w" );
  fopenerr = errno;
  XLAL_CHECK_MAIN( fp4 != NULL, XLAL_EIO, "Failed to open '%s' for writing: %s", outfile4, strerror( fopenerr ) );

  // Record timestamps for each SFT or gap
  // This is not known a priori so we end up resizing this vector as we go
  XLAL_CHECK_MAIN( ( timestamps = XLALCreateINT4Vector( 0 ) ) != NULL, XLAL_EFUNC );

  // Need to save the f0 and deltaF from the first SFT
  REAL8 f0 = 0, deltaF = 0;

  // Initialize a temporary frequency series for XLALNormalizeSFT
  REAL8FrequencySeries *rngmed = NULL;

  printf( "Looping over SFTs to compute average spectra and spectrograms\n" );
  for ( UINT4 j = 0; j < catalog->length; j++ ) {
    fprintf( stderr, "Extracting SFT %d...\n", j );

    //Extract one SFT at a time from the catalog
    //we do this by using a catalog timeslice to get just the current SFT
    XLAL_CHECK_MAIN( ( sft_vect = extract_one_sft( catalog, catalog->data[j].header.epoch, f_min, f_max ) ) != NULL, XLAL_EFUNC );
    XLAL_CHECK_MAIN( sft_vect->length == 1, XLAL_EINVAL, "Extracted zero SFTs but should have extracted one" );

    //Make sure the SFTs are the same length as what we're expecting from user input
    XLAL_CHECK_MAIN( fabs( timebaseline * sft_vect->data->deltaF - 1.0 ) <= 10.*LAL_REAL8_EPS, XLAL_EINVAL, "Expected SFTs with length %f but got %f", timebaseline, 1 / sft_vect->data->deltaF );

    // first time through the loop need to allocate timeavg, rngmed
    // save some scalar values
    if ( j == 0 ) {
      // Allocate the vector for the normalized data to be averaged together
      XLAL_CHECK_MAIN( ( timeavg = XLALCreateREAL4Vector( sft_vect->data->data->length ) ) != NULL, XLAL_EFUNC );
      memset( timeavg->data, 0, sizeof( REAL4 )*timeavg->length );

      // Compute the number of bins in an average; at minimum this should be 1 bin
      // It produces the same frequency resolution as specified in the arguments passed to
      // fscanDriver.py
      REAL8 spectrogram_blocksize = round( freqres * sft_vect->data->data->length / ( f_max - f_min ) );
      if ( spectrogram_blocksize < 1.0 ) {
        NumBinsAvg = 1;
      } else {
        NumBinsAvg = ( UINT4 )spectrogram_blocksize;
      }

      f0 = sft_vect->data->f0;
      deltaF = sft_vect->data->deltaF;

      // Allocate frequency series for XLALNormalizeSFT; meta data doesn't matter
      XLAL_CHECK_MAIN( ( rngmed = XLALCreateREAL8FrequencySeries( sft_vect->data->name, &catalog->data[j].header.epoch, f0, deltaF, &sft_vect->data->sampleUnits, sft_vect->data->data->length ) ) != NULL, XLAL_EFUNC );

      fprintf( stderr, "SFTs = %d\tSFT bins = %d\tf0 = %f\n", catalog->length, sft_vect->data->data->length, sft_vect->data->f0 );
    }

    // GPS time of the current SFT and print to file
    cur_epoch = sft_vect->data->epoch.gpsSeconds;
    fprintf( fp2, "%d\t%d\n", timestamps->length, cur_epoch );

    // Get the current UTC time from the GPS seconds of the current SFT and print to file
    XLAL_CHECK_MAIN( XLALGPSToUTC( &date, cur_epoch ) != NULL, XLAL_EFUNC );
    fprintf( fp4, "%d\t %i\t %i\t %i\t %i\t %i\t %i\n", timestamps->length, ( date.tm_year + 1900 ), date.tm_mon + 1, date.tm_mday, date.tm_hour, date.tm_min, date.tm_sec );

    // Here is where the timestamps vector is resized and this epoch is recorded
    XLAL_CHECK_MAIN( XLALResizeINT4Vector( timestamps, timestamps->length + 1 ) != NULL, XLAL_EFUNC );
    timestamps->data[timestamps->length - 1] = cur_epoch;

    /* Spectrogram (course grained) */
    // First line will be the frequencies of the spectrogram indicating what is in each column
    // The first value is 0 simply because column 0 are the GPS times. It keeps the entire
    // output numeric, in case that matters.
    if ( j == 0 ) {
      UINT4 step = 0;
      for ( UINT4 i = NumBinsAvg - 1; i < sft_vect->data->data->length; i += NumBinsAvg ) {
        if ( step == 0 ) {
          fprintf( fp, "0\t" );
        }
        fprintf( fp, "%.6f\t", sft_vect->data->f0 + NumBinsAvg * sft_vect->data->deltaF * step );
        step++;
      }
      fprintf( fp, "\n" );
    }

    // Loop over the number of bins in each SFT, jumping by the average number
    // Start at the highest bin index and average down. This avoids running off the end of the SFT
    for ( UINT4 i = NumBinsAvg - 1; i < sft_vect->data->data->length; i += NumBinsAvg ) {
      // First value in each row is the GPS time
      if ( i == NumBinsAvg - 1 ) {
        fprintf( fp, "%i\t", sft_vect->data->epoch.gpsSeconds );
      }
      REAL8 avg = 0.0;
      for ( UINT4 k = 0; k < NumBinsAvg; k++ ) {
        const REAL8 re = ( REAL8 )crealf( sft_vect->data->data->data[i - k] );
        const REAL8 im = ( REAL8 )cimagf( sft_vect->data->data->data[i - k] );
        avg += 2.0 * ( re * re + im * im ) / ( REAL8 )timebaseline;
      }
      fprintf( fp, "%e\t", sqrt( avg / ( REAL8 )NumBinsAvg ) ); /* 06/15/2017 gam; then take sqrt here. */
    }
    fprintf( fp, "\n" );

    // Fill in gaps where there is no SFT data with zeros
    if ( j < ( catalog->length - 1 ) ) { /*in all cases except when we are examining the last sft, check that there is no gap to the next sft*/
      /*test to see if the next SFT immediately follows, if not entries in the matrix until there is one*/
      while ( cur_epoch + timebaseline < catalog->data[j + 1].header.epoch.gpsSeconds ) {
        cur_epoch += timebaseline;
        fprintf( fp2, "%d.\t%d\n", timestamps->length, cur_epoch );

        XLAL_CHECK_MAIN( XLALGPSToUTC( &date, cur_epoch ) != NULL, XLAL_EFUNC );
        fprintf( fp4, "%d\t %i\t %i\t %i\t %i\t %i\t %i\n", timestamps->length, ( date.tm_year + 1900 ), date.tm_mon + 1, date.tm_mday, date.tm_hour, date.tm_min, date.tm_sec );
        XLAL_CHECK_MAIN( XLALResizeINT4Vector( timestamps, timestamps->length + 1 ) != NULL, XLAL_EFUNC );
        timestamps->data[timestamps->length - 1] = cur_epoch;

        for ( UINT4 i = NumBinsAvg - 1; i < sft_vect->data->data->length; i += NumBinsAvg ) {
          // First value in each row is the GPS time
          if ( i == NumBinsAvg - 1 ) {
            fprintf( fp, "%i\t", cur_epoch );
          }
          REAL8 avg = 0.0;
          fprintf( fp, "%e\t", avg );
        }
        fprintf( fp, "\n" );
      }

    }

    /* Normalized, averaged spectra */
    // Normalize the SFT
    XLAL_CHECK_MAIN( XLALNormalizeSFT( rngmed, sft_vect->data, blocksRngMed, 0.0 ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Loop over the frequency bins in the SFT
    for ( UINT4 i = 0; i < sft_vect->data->data->length; i++ ) {
      const REAL8 re = ( REAL8 )crealf( sft_vect->data->data->data[i] );
      const REAL8 im = ( REAL8 )cimagf( sft_vect->data->data->data[i] );
      timeavg->data[i] += re * re + im * im;
    }

    // Destroys current SFT Vector
    XLALDestroySFTVector( sft_vect );
    sft_vect = NULL;
  }
  fprintf( stderr, "finished checking for missing sfts, l=%d\n", timestamps->length );


  /*----------------------------------------------------------------------------------------------------------------*/

  // Write the averaged data to files broken up by the user option subband (default: 100 Hz)
  for ( UINT4 i = 0; i < timeavg->length; i++ ) {
    REAL8 f_min_subband = f0 + ( ( REAL4 )i ) * deltaF;
    REAL8 f_max_subband = f_min + subband;
    if ( f_max_subband > f_max ) {
      f_max_subband = f_max;
    }

    if ( fp3 == NULL ) {
      if ( XLALUserVarWasSet( &outputBname ) ) {
        snprintf( outbase, sizeof( outbase ), "%s/%s_%.2f_%.2f", outputDir, outputBname, f_min_subband, f_max_subband );
      } else {
        snprintf( outbase, sizeof( outbase ), "%s/spec_%.2f_%.2f_%s_%d_%d", outputDir, f_min_subband, f_max_subband, constraints.detector, startTime.gpsSeconds, endTime.gpsSeconds );
      }
      snprintf( outfile3, sizeof( outfile3 ), "%s_timeaverage.txt", outbase );
      fp3 = fopen( outfile3, "w" );
      fopenerr = errno;
      XLAL_CHECK_MAIN( fp3 != NULL, XLAL_EIO, "Failed to open '%s' for writing: %s", outfile3, strerror( fopenerr ) );
    }

    REAL8 f = f0 + ( ( REAL4 )i ) * deltaF;
    REAL4 PWR = timeavg->data[i] / ( ( REAL4 )catalog->length );
    fprintf( fp3, "%16.6f %16.3f \n", f, PWR );

    if ( f + deltaF >= f_max_subband ) {
      fclose( fp3 );
      fp3 = NULL;
    }
  }

  /*------------------------------------------------------------------------------------------------------------------------*/
  /*End of normal spec_avg code.*/

  XLALDestroySFTCatalog( catalog );
  XLALDestroyREAL4Vector( timeavg );
  XLALDestroyINT4Vector( timestamps );
  XLALDestroyREAL8FrequencySeries( rngmed );

  XLALDestroyUserVars();

  /*close all the files, spec_avg.c is done, all info written to the files.*/
  fclose( fp );
  fclose( fp2 );
  fclose( fp4 );

  fprintf( stderr, "end of spec_avg\n" );

  return ( 0 );


}
/* END main */
