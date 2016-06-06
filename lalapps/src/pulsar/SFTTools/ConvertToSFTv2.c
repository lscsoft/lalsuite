/*
 * Copyright (C) 2006 Reinhard Prix
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
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */


/**
 * \author Reinhard Prix
 * \date 2006
 * \file
 * \ingroup lalapps_pulsar_SFTTools
 * \brief Code to convert given input-SFTs (v1 or v2) to v2-SFTs with given extra-comment,
 * and write them out following the SFTv2 naming conventions (see LIGO-T040164-01-Z)
 */

/* ---------- includes ---------- */
#include <sys/stat.h>
#include <sys/types.h>

#include <lalapps.h>

#include <lal/UserInput.h>
#include <lal/PulsarDataTypes.h>
#include <lal/SFTfileIO.h>
#include <lal/SFTutils.h>
#include <lal/LogPrintf.h>

/*---------- DEFINES ----------*/
/*----- Macros ----- */
/*---------- internal types ----------*/

/*---------- empty initializers ---------- */

/*---------- Global variables ----------*/

/* User variables */
typedef struct
{
  CHAR *inputSFTs;
  CHAR *outputDir;
  CHAR *outputSingleSFT;
  CHAR *extraComment;
  CHAR *descriptionMisc;
  CHAR *IFO;
  REAL8 fmin;
  REAL8 fmax;
  REAL8 mysteryFactor;
  INT4 minStartTime;
  INT4 maxStartTime;
  CHAR *timestampsFile;
} UserInput_t;


/*---------- internal prototypes ----------*/
int initUserVars ( UserInput_t *uvar );
int applyFactor2SFTs ( SFTVector *SFTs, REAL8 factor );

/*==================== FUNCTION DEFINITIONS ====================*/

/*----------------------------------------------------------------------
 * main function
 *----------------------------------------------------------------------*/
int
main(int argc, char *argv[])
{
  SFTConstraints XLAL_INIT_DECL(constraints);
  LIGOTimeGPS XLAL_INIT_DECL(minStartTimeGPS);
  LIGOTimeGPS XLAL_INIT_DECL(maxStartTimeGPS);
  SFTCatalog *FullCatalog = NULL;
  CHAR *add_comment = NULL;
  UINT4 i;
  REAL8 fMin, fMax;

  UserInput_t XLAL_INIT_DECL(uvar);

  /* register all user-variables */
  XLAL_CHECK_MAIN ( initUserVars ( &uvar ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* read cmdline & cfgfile  */
  BOOLEAN should_exit = 0;
  XLAL_CHECK_MAIN( XLALUserVarReadAllInput( &should_exit, argc, argv ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( should_exit ) {
    exit (1);
  }

  /* ----- make sure output directory exists ---------- */
  if ( uvar.outputDir )
  {
    int ret;
    ret = mkdir ( uvar.outputDir, 0777);
    if ( (ret == -1) && ( errno != EEXIST ) )
      {
	int errsv = errno;
	LogPrintf (LOG_CRITICAL, "Failed to create directory '%s': %s\n", uvar.outputDir, strerror(errsv) );
	return -1;
      }
  }

  LIGOTimeGPSVector *timestamps = NULL;
  if ( uvar.timestampsFile )
    {
      if ( (timestamps = XLALReadTimestampsFile ( uvar.timestampsFile )) == NULL ) {
        XLALPrintError ("XLALReadTimestampsFile() failed to load timestamps from file '%s'\n", uvar.timestampsFile );
        return -1;
      }
    }

  /* use IFO-contraint if one given by the user */
  if ( LALUserVarWasSet ( &uvar.IFO ) ) {
    XLAL_CHECK_MAIN ( (constraints.detector = XLALGetChannelPrefix ( uvar.IFO )) != NULL, XLAL_EINVAL );
  }

  minStartTimeGPS.gpsSeconds = uvar.minStartTime;
  maxStartTimeGPS.gpsSeconds = uvar.maxStartTime;
  constraints.minStartTime = &minStartTimeGPS;
  constraints.maxStartTime = &maxStartTimeGPS;
  constraints.timestamps = timestamps;

  /* get full SFT-catalog of all matching (multi-IFO) SFTs */
  XLAL_CHECK_MAIN ( (FullCatalog = XLALSFTdataFind ( uvar.inputSFTs, &constraints )) != NULL, XLAL_EFUNC );

  if ( constraints.detector ) {
    XLALFree ( constraints.detector );
  }

  XLAL_CHECK_MAIN ( (FullCatalog != NULL) && (FullCatalog->length > 0), XLAL_EINVAL, "\nSorry, didn't find any matching SFTs with pattern '%s'!\n\n", uvar.inputSFTs );

  /* build up full comment-string to be added to SFTs: 1) converted by ConvertToSFTv2, VCS ID 2) user extraComment */
  {
    UINT4 len = 128;
    len += strlen ( uvar.inputSFTs );
    if ( uvar.extraComment )
      len += strlen ( uvar.extraComment );

    XLAL_CHECK_MAIN ( ( add_comment = LALCalloc ( 1, len )) != NULL, XLAL_ENOMEM );

    /** \deprecated FIXME: the following code uses obsolete CVS ID tags.
     *  It should be modified to use git version information. */
    sprintf ( add_comment, "Converted by $Id$, inputSFTs = '%s';", uvar.inputSFTs );
    if ( uvar.extraComment )
      {
	strcat ( add_comment, "\n");
	strcat ( add_comment, uvar.extraComment );
      }
  } /* construct comment-string */

  /* which frequency-band to extract? */
  fMin = -1;	/* default: all */
  fMax = -1;
  if ( LALUserVarWasSet ( &uvar.fmin ) )
    fMin = uvar.fmin;
  if ( LALUserVarWasSet ( &uvar.fmax ) )
    fMax = uvar.fmax;

  FILE *fpSingleSFT = NULL;
  if ( uvar.outputSingleSFT )
    XLAL_CHECK ( ( fpSingleSFT = fopen ( uvar.outputSingleSFT, "wb" )) != NULL,
                 XLAL_EIO, "Failed to open singleSFT file '%s' for writing\n", uvar.outputSingleSFT );

  /* loop over all SFTs in SFTCatalog */
  for ( i=0; i < FullCatalog->length; i ++ )
    {
      SFTCatalog oneSFTCatalog;
      SFTVector *thisSFT = NULL;
      const CHAR *sft_comment;
      CHAR *new_comment;
      UINT4 comment_len = 0;

      /* set catalog containing only one SFT */
      oneSFTCatalog.length = 1;
      oneSFTCatalog.data = &(FullCatalog->data[i]);

      comment_len = strlen ( add_comment ) + 10;
      sft_comment = oneSFTCatalog.data->comment;
      if ( sft_comment )
	comment_len += strlen ( sft_comment );

      XLAL_CHECK_MAIN ( ( new_comment  = LALCalloc (1, comment_len )) != NULL, XLAL_ENOMEM );

      if ( sft_comment ) {
	strcpy ( new_comment, sft_comment );
	strcat ( new_comment, ";\n");
      }
      strcat ( new_comment, add_comment );

      XLAL_CHECK_MAIN ( (thisSFT = XLALLoadSFTs ( &oneSFTCatalog, fMin, fMax )) != NULL, XLAL_EFUNC );

      if ( uvar.mysteryFactor != 1.0 ) {
	XLAL_CHECK_MAIN ( applyFactor2SFTs ( thisSFT, uvar.mysteryFactor ) == XLAL_SUCCESS, XLAL_EFUNC );
      }

      // if user asked for single-SFT output, add this SFT to the open file
      if ( uvar.outputSingleSFT )
        XLAL_CHECK ( XLAL_SUCCESS == XLALWriteSFT2fp( &(thisSFT->data[0]), fpSingleSFT, new_comment ),
                     XLAL_EFUNC,  "XLALWriteSFT2fp() failed to write SFT to '%s'!\n", uvar.outputSingleSFT );

      // if user asked for directory output, write this SFT into that directory
      if ( uvar.outputDir ) {
        XLAL_CHECK_MAIN ( XLALWriteSFTVector2Dir ( thisSFT, uvar.outputDir, new_comment, uvar.descriptionMisc ) == XLAL_SUCCESS, XLAL_EFUNC );
      }

      XLALDestroySFTVector ( thisSFT );

      XLALFree ( new_comment );

    } /* for i < numSFTs */

  if ( fpSingleSFT ) {
    fclose ( fpSingleSFT );
  }

  /* free memory */
  XLALFree ( add_comment );
  XLALDestroySFTCatalog ( FullCatalog );
  XLALDestroyUserVars();
  XLALDestroyTimestampVector ( timestamps );

  LALCheckMemoryLeaks();

  return 0;
} /* main */


/*----------------------------------------------------------------------*/
/* register all our "user-variables" */
int
initUserVars ( UserInput_t *uvar )
{
  XLAL_CHECK ( uvar != NULL, XLAL_EINVAL );

  /* set defaults */
  uvar->outputDir = NULL;
  uvar->outputSingleSFT = NULL;

  uvar->extraComment = NULL;
  uvar->descriptionMisc = NULL;
  uvar->IFO = NULL;

  uvar->minStartTime = 0;
  uvar->maxStartTime = LAL_INT4_MAX;

  uvar->mysteryFactor = 1.0;

  uvar->timestampsFile = NULL;

  /* now register all our user-variable */
  XLALRegisterUvarMember( inputSFTs,	STRING, 'i', REQUIRED, "File-pattern for input SFTs");
  XLALRegisterUvarMember( IFO,		STRING, 'I', OPTIONAL, "IFO of input SFTs: 'G1', 'H1', 'H2', ...(required for v1-SFTs)");

  XLALRegisterUvarMember( outputSingleSFT,	STRING, 'O', OPTIONAL, "Output all SFTs into a single concatenated SFT-file with this name");
  XLALRegisterUvarMember( outputDir,	STRING, 'o', OPTIONAL, "Output directory for SFTs");

  XLALRegisterUvarMember( extraComment,	STRING, 'C', OPTIONAL, "Additional comment to be added to output-SFTs");

  XLALRegisterUvarMember( descriptionMisc,	STRING, 'D', OPTIONAL, "'Misc' entry in the SFT filename description-field (see SFTv2 naming convention)");
  XLALRegisterUvarMember(   fmin,		REAL8, 'f', OPTIONAL, "Lowest frequency to extract from SFTs. [Default: lowest in inputSFTs]");
  XLALRegisterUvarMember(   fmax,		REAL8, 'F', OPTIONAL, "Highest frequency to extract from SFTs. [Default: highest in inputSFTs]");


  XLALRegisterUvarMember(  	minStartTime, 	 INT4, 0,  OPTIONAL, "Only use SFTs with timestamps starting from (including) this GPS time");
  XLALRegisterUvarMember(  	maxStartTime, 	 INT4, 0,  OPTIONAL, "Only use SFTs with timestamps up to (excluding) this GPS time");

  XLALRegisterUvarMember( timestampsFile,	 STRING, 0, OPTIONAL, "Timestamps file to use as a constraint for SFT loading");

  /* developer-options */
  XLALRegisterUvarMember(   mysteryFactor,	 REAL8, 0, DEVELOPER, "Change data-normalization by applying this factor (for E@H)");

  return XLAL_SUCCESS;

} /* initUserVars() */

int
applyFactor2SFTs ( SFTVector *SFTs, REAL8 factor )
{

  XLAL_CHECK ( SFTs != NULL, XLAL_EINVAL );

  UINT4 numSFTs = SFTs->length;

  for ( UINT4 i=0; i < numSFTs; i ++ )
    {
      SFTtype *thisSFT = &(SFTs->data[i]);
      UINT4 k, numBins = thisSFT->data->length;

      for ( k=0; k < numBins; k ++ )
	{
	  thisSFT->data->data[k] *= ((REAL4) factor);
	} /* for k < numBins */

    } /* for i < numSFTs */

  return XLAL_SUCCESS;

} /* applyFactor2SFTs() */
