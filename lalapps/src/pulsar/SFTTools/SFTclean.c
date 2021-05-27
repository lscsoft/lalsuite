/*
*  Copyright (C) 2021 David Keitel
*  Copyright (C) 2007 Badri Krishnan, Alicia Sintes Olives
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
 * \ingroup lalapps_pulsar_SFTTools
 * \author Badri Krishnan, Alicia Sintes Olives, David Keitel
 *
 * NOTE: The modern Advanced LIGO linefiles format is not supported.
 * To run this code, any linefiles should first be converted into the legacy
 * format as described in the lalpulsar SFTClean module.
 */

#include "config.h"

#include <lal/SFTClean.h>
#include <lal/SFTutils.h>
#include <lal/SFTReferenceLibrary.h>
#include <glob.h> 
#include <LALAppsVCSInfo.h>

/* Error codes and messages */

/**\name Error Codes */ /** @{ */
#define SFTCLEANC_ENORM 0
#define SFTCLEANC_ESUB  1
#define SFTCLEANC_EARG  2
#define SFTCLEANC_EBAD  3
#define SFTCLEANC_EFILE 4

#define SFTCLEANC_MSGENORM "Normal exit"
#define SFTCLEANC_MSGESUB  "Subroutine failed"
#define SFTCLEANC_MSGEARG  "Error parsing arguments"
#define SFTCLEANC_MSGEBAD  "Bad argument values"
#define SFTCLEANC_MSGEFILE "Could not create output file"
/** @} */

/*---------- local defines ---------- */
#define MAXFILENAMELENGTH 256
#define TRUE (1==1)
#define FALSE (1==0)

#define CMT_NONE 0
#define CMT_OLD  1
#define CMT_FULL 2

int main(int argc, char *argv[]){ 

  static LALStatus       status;  /* LALStatus pointer */ 

  static MultiSFTVector  *inputSFTs = NULL;
  static SFTCatalog *catalog = NULL;
  static SFTCatalog thisCatalog;
  static SFTConstraints constraints;
  const char *misc = "cleaned"; /* used for setting output filenames */

  /* 09/09/05 gam; randPar now a parameter for LALCleanCOMPLEX8SFT */
  FILE *fp=NULL;   
  INT4 seed, ranCount;
  RandomParams *randPar=NULL; 

  /* user input variables */
  LALStringVector *uvar_linefiles=NULL; /* files with harmonics info */
  CHAR *uvar_sftDir;    /* directory for unclean sfts */
  CHAR *uvar_outDir;   /* directory for cleaned sfts */
  REAL8 uvar_fMin, uvar_fMax;
  INT4  uvar_window, uvar_maxBins, uvar_addComment;
  BOOLEAN uvar_outSingleSFT;

  /* set defaults */
  uvar_sftDir = NULL;
  uvar_outDir = NULL;
  uvar_fMin = -1;
  uvar_fMax = -1;
  uvar_window = 100;
  uvar_maxBins = 20;
  uvar_addComment = CMT_FULL; /* add VCS ID and full command-line to every SFT file */
  uvar_outSingleSFT = FALSE;

  /* register user input variables */
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_sftDir,    "sftDir",    STRING,       'i', REQUIRED,  "Input SFT file pattern") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_outDir,    "outDir",    STRING,       'o', REQUIRED,  "Output SFT Directory") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_fMin,      "fMin",      REAL8,        0,   OPTIONAL, "start Frequency (default: full input SFTs width)") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_fMax,      "fMax",      REAL8,        0,   OPTIONAL, "Max Frequency  (default: full input SFTs width)") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_window,    "window",    INT4,         'w',OPTIONAL,  "Window size for noise floor estimation in vicinity of a line. WARNING: estimation will be compromised if SFTs are narrower than this, or if it is not much wider than any line wings.") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_maxBins,   "maxBins",   INT4,         'm',OPTIONAL,  "Max. bins to clean. WARNING: If your linefiles contain very wide lines, make sure this is large enough to accommodate them.") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_linefiles, "linefiles", STRINGVector, 0,   OPTIONAL, "List of per-detector files with list of lines (each full path must start with a canonical IFO name). NOTE: The modern Advanced LIGO linefiles format is not supported. Files should first be converted into the legacy format as described in the lalpulsar SFTClean module.") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar( &uvar_addComment, "addComment", INT4,       'c', OPTIONAL, "How to deal with comments - 0 means no comment is written at all, 1 means that the comment is taken unmodified from the input SFTs, 2 (default) means that the program appends its RCS id and command-line to the comment.") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN( XLALRegisterNamedUvar ( &uvar_outSingleSFT, "outSingleSFT", BOOLEAN, 's', OPTIONAL, "Write a single concatenated SFT file per IFO, instead of individual files") == XLAL_SUCCESS, XLAL_EFUNC );

  /* read all command line variables */
  BOOLEAN should_exit = 0;
  XLAL_CHECK_MAIN( XLALUserVarReadAllInput(&should_exit, argc, argv, lalAppsVCSInfoList) == XLAL_SUCCESS, XLAL_EFUNC);
  if (should_exit)
    exit(1);

  /* record VCS ID and command-line for the comment */
  char *cmdline = NULL;
  if ( uvar_addComment > CMT_OLD ) {
    XLAL_CHECK_MAIN( ( cmdline = ( char * )XLALMalloc( strlen( lalAppsVCSIdentInfo.vcsId ) + strlen( lalAppsVCSIdentInfo.vcsStatus ) + 2 ) ) != NULL, XLAL_ENOMEM, "out of memory allocating cmdline" );
    strcpy( cmdline, lalAppsVCSIdentInfo.vcsId );
    strcat( cmdline, lalAppsVCSIdentInfo.vcsStatus );
    strcat( cmdline, "\n" );
    for ( int arg = 0; arg < argc; arg++ ) {
        XLAL_CHECK_MAIN( ( cmdline = ( char * )XLALRealloc( ( void * )cmdline, strlen( cmdline ) + strlen( argv[arg] ) + 2 ) ) != NULL, XLAL_ENOMEM, "out of memory allocating cmdline" );
        strcat( cmdline, argv[arg] );
        if ( arg == argc - 1 ) {
        strcat( cmdline, "\n" );
        } else {
        strcat( cmdline, " " );
        }
    }
  }

  /* set detector constraint */
  constraints.detector = NULL;

  /* get SFT catalog - this can be for multiple detectors */
  XLAL_CHECK_MAIN( ( catalog = XLALSFTdataFind( uvar_sftDir, &constraints) ) != NULL, XLAL_EFUNC);
  if ( (catalog == NULL) || (catalog->length == 0) ) {
    fprintf (stderr,"Unable to match any SFTs with pattern '%s'\n", uvar_sftDir );
    exit(1);
  }
  fprintf(stdout, "Created catalog with %d SFTs.\n", catalog->length);

  /* get a new seed value, and use it to create a new random parameter structure */
  fp=fopen("/dev/urandom", "r");
  if (!fp) 
    { 
      fprintf(stderr,"Error in opening /dev/urandom \n"); 
      exit(1); 
    } 
  ranCount = fread(&seed, sizeof(seed), 1, fp);
  fclose(fp);
  if (!(ranCount==1)) 
    { 
      fprintf(stderr,"Error in reading random seed \n"); 
      exit(1); 
    } 
  LAL_CALL ( LALCreateRandomParams (&status, &randPar, seed), &status );

  /* loop over all SFTs in the full catalog, but sorted by detector */
  MultiSFTCatalogView* multiCatalogView = NULL;
  XLAL_CHECK ( (multiCatalogView = XLALGetMultiSFTCatalogView ( catalog )) != NULL, XLAL_EFUNC );
  for (UINT4 X = 0; X < multiCatalogView->length; X++) {
  
    /* optionally prepare merged-file output:
     * if uvar_outSingleSFT, a file pointer will be opened and all SFTs for X
     * will be appended to it
     */
    char *outpath = NULL;
    FILE *fpout = NULL;
    UINT4 numSFTs = multiCatalogView->data[X].length;
    if ( uvar_outSingleSFT ) {
      /* grab the first and last entry from the single-IFO catalog,
       * relying here on XLALSFTdataFind returning a catalogue with SFTs sorted by increasing GPS-epochs
       */
      SFTtype *sftStart = &(multiCatalogView->data[X].data[0].header);
      SFTtype *sftEnd = &(multiCatalogView->data[X].data[numSFTs-1].header);
      LIGOTimeGPS *epochStart = &(sftStart->epoch);
      LIGOTimeGPS *epochEnd   = &(sftEnd->epoch);
      const char *name = sftStart->name;
      UINT4 Tsft = (UINT4) round ( 1.0 / sftStart->deltaF );
      /* calculate time interval covered -- may be different from timebase if nanosecond of sft-epochs are non-zero */
      UINT4 Tspan = epochEnd->gpsSeconds - epochStart->gpsSeconds + Tsft;
      if ( epochStart->gpsNanoSeconds > 0) {
        Tspan += 1;
      }
      if ( epochEnd->gpsNanoSeconds > 0) {
        Tspan += 1;
      }
      /* get the official-format filename for the merged output file for this detector */
      char *filename;
      XLAL_CHECK_MAIN ( (filename = XLALOfficialSFTFilename ( name[0], name[1], numSFTs, Tsft, epochStart->gpsSeconds, Tspan, misc )) != NULL, XLAL_EFUNC );
      /* append to output directory name */
      int len = strlen ( uvar_outDir ) + 1 + strlen ( filename ) + 1;
      XLAL_CHECK_MAIN ( (outpath = XLALCalloc ( 1, len )) != NULL, XLAL_ENOMEM );
      sprintf ( outpath, "%s/%s", uvar_outDir, filename );
      XLALFree ( filename );
      printf ( "Writing %d cleaned SFTs for %s to merged output file %s ...\n", numSFTs, name, outpath );
      XLAL_CHECK_MAIN ( (fpout = fopen ( outpath, "wb" )) != NULL, XLAL_EIO, "Failed to open '%s' for writing: %s\n\n", outpath, strerror(errno));
    } else {
      printf ( "Writing %d cleaned SFTs for %s to individual output files in %s ...\n", numSFTs, multiCatalogView->data[X].data[0].header.name, uvar_outDir );
    }

    /* loop over SFTs for this IFO and clean them -- load one SFT at a time */
    thisCatalog.length = 1;
    for (UINT4 j=0; j<numSFTs; j++) {

      thisCatalog.data = multiCatalogView->data[X].data + j;

      /* read a multi-SFT vector, but it will only cover a single IFO and a single SFT */
      XLAL_CHECK_MAIN( ( inputSFTs = XLALLoadMultiSFTs ( &thisCatalog, uvar_fMin, uvar_fMax) ) != NULL, XLAL_EFUNC);
      XLAL_CHECK_MAIN( inputSFTs->length == 1, XLAL_EIO, "We should have data from a single IFO here, but have inputSFTs->length==%d", inputSFTs->length );

      /* construct comment for output SFTs */
      char *comment = NULL;
      char *oldcomment = thisCatalog.data[0].comment;
      int comment_length = strlen(oldcomment) + 1;
      if ( uvar_addComment > CMT_OLD ) {
        /* allocate space for new comment */
        XLAL_CHECK_MAIN( ( comment = ( char * )XLALMalloc( comment_length + strlen( cmdline ) + 1 ) ) != NULL, XLAL_ENOMEM, "out of memory allocating comment" );
        /* append the commandline of this program to the old comment */
        if ( oldcomment ) {
            strcpy( comment, oldcomment );
        } else {
            *comment = '\0';
        }
        strcat( comment, cmdline );
      } else if ( uvar_addComment == CMT_OLD ) {
        /* only copied existing comment, no additional space needed */
        comment = oldcomment;
      } /* else (uvar_addComment == CMT_NONE) and (comment == NULL) i.e. no comment at all */

      /* clean lines
       * Here we still pass the full, possible multi-IFO linesfiles,
       * but LALRemoveKnownLinesInMultiSFTVector() has logic to properly apply
       * only those files with canonic IFO names matching the single-IFO inputSFTs.
       */
      if ( XLALUserVarWasSet( &uvar_linefiles ) ) {
        LAL_CALL( LALRemoveKnownLinesInMultiSFTVector ( &status, inputSFTs, uvar_maxBins,
                      uvar_window, uvar_linefiles, randPar), &status);
      }
  
      /* write output - technically looping over all SFTs for a single detector here,
       * but actually the current implementation has it always be a single SFT too
       */
      if ( uvar_outSingleSFT ) {
        for ( UINT4 k = 0; k < inputSFTs->data[0]->length; k++ ) {
          SFTtype *this_sft = &( inputSFTs->data[0]->data[k] );
          XLAL_CHECK ( XLALWriteSFT2fp ( this_sft, fpout, comment ) == XLAL_SUCCESS, XLAL_EFUNC );
        }
      } else {
        XLAL_CHECK_MAIN( XLALWriteSFTVector2Dir ( inputSFTs->data[0], uvar_outDir, comment, misc ) == XLAL_SUCCESS, XLAL_EFUNC);
      }

      XLALDestroyMultiSFTVector( inputSFTs);
      if ( uvar_addComment > CMT_OLD ) {
        XLALFree( comment );
      }

    } /* end loop over sfts */

    if (uvar_outSingleSFT ) {
      fclose( fpout );
      /* validate the merged output file now that we're finished with it:
       * in principle we have made sure all SFT timestamps are ascending etc,
       * but best to make sure
       */
      XLAL_CHECK_MAIN ( ValidateSFTFile(outpath) == 0, XLAL_EFUNC, "Output file %s failed validation.", outpath );
      XLALFree( outpath );
    }

  } /* for (UINT4 X = 0; X < multiCatalogView->length; X++) */

  /* Free memory */
  XLALDestroySFTCatalog(catalog );
  XLALDestroyMultiSFTCatalogView ( multiCatalogView );
  XLALDestroyUserVars();
  LAL_CALL( LALDestroyRandomParams (&status, &randPar), &status);
  if ( uvar_addComment > CMT_OLD ) {
    XLALFree( cmdline );
  }

  LALCheckMemoryLeaks(); 

  return 0;
}
