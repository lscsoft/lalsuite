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
 * \ingroup pulsarApps
 * \brief Code to convert given input-SFTs (v1 or v2) to v2-SFTs with given extra-comment,
 *        and write them out following the SFTv2 naming conventions (see LIGO-T040164-01-Z)
 */

/* ---------- includes ---------- */
#include <sys/stat.h>
#include <sys/types.h>

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lalapps.h>

#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/LogPrintf.h>

/** \name Error codes */
/*@{*/
#define CONVERTSFT_ENORM 	0
#define CONVERTSFT_EINPUT  	1
#define CONVERTSFT_EMEM		2

#define CONVERTSFT_MSGENORM 	"Normal exit"
#define CONVERTSFT_MSGEINPUT  	"Bad argument values"
#define CONVERTSFT_MSGEMEM	"Out of memory"
/*@}*/

/*---------- DEFINES ----------*/

#define TRUE    (1==1)
#define FALSE   (1==0)

#define LAL_INT4_MAX 2147483647
/*----- Macros ----- */
/*---------- internal types ----------*/

/*---------- empty initializers ---------- */

/*---------- Global variables ----------*/

/* User variables */
BOOLEAN uvar_help;
CHAR *uvar_inputSFTs;
CHAR *uvar_outputDir;
CHAR *uvar_outputSingleSFT;
CHAR *uvar_extraComment;
CHAR *uvar_descriptionMisc;
CHAR *uvar_IFO;
REAL8 uvar_fmin;
REAL8 uvar_fmax;
REAL8 uvar_mysteryFactor;

INT4 uvar_minStartTime;
INT4 uvar_maxEndTime;


/*---------- internal prototypes ----------*/
void initUserVars (LALStatus *status);
void applyFactor2SFTs ( LALStatus *status, SFTVector *SFTs, REAL8 factor );

/*==================== FUNCTION DEFINITIONS ====================*/

/*----------------------------------------------------------------------
 * main function
 *----------------------------------------------------------------------*/
int
main(int argc, char *argv[])
{
  LALStatus status = blank_status;	/* initialize status */
  SFTConstraints constraints = empty_SFTConstraints;
  LIGOTimeGPS minStartTimeGPS = empty_LIGOTimeGPS;
  LIGOTimeGPS maxEndTimeGPS = empty_LIGOTimeGPS;
  SFTCatalog *FullCatalog = NULL;
  CHAR *add_comment = NULL;
  UINT4 i;
  REAL8 fMin, fMax;


  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;	/* exit with returned status-code on error */

  /* register all user-variables */
  LAL_CALL (initUserVars (&status), &status);

  /* read cmdline & cfgfile  */
  LAL_CALL (LALUserVarReadAllInput (&status, argc,argv), &status);

  if (uvar_help) 	/* help requested: we're done */
    exit (0);

  /* ----- make sure output directory exists ---------- */
  if ( uvar_outputDir )
  {
    int ret;
    ret = mkdir ( uvar_outputDir, 0777);
    if ( (ret == -1) && ( errno != EEXIST ) )
      {
	int errsv = errno;
	LogPrintf (LOG_CRITICAL, "Failed to create directory '%s': %s\n", uvar_outputDir, strerror(errsv) );
	return -1;
      }
  }

  /* use IFO-contraint if one given by the user */
  if ( LALUserVarWasSet ( &uvar_IFO ) ) {
    if ( (constraints.detector = XLALGetChannelPrefix ( uvar_IFO )) == NULL ) {
      return CONVERTSFT_EINPUT;
    }
  }

  minStartTimeGPS.gpsSeconds = uvar_minStartTime;
  maxEndTimeGPS.gpsSeconds = uvar_maxEndTime;
  constraints.startTime = &minStartTimeGPS;
  constraints.endTime = &maxEndTimeGPS;

  /* get full SFT-catalog of all matching (multi-IFO) SFTs */
  LAL_CALL ( LALSFTdataFind ( &status, &FullCatalog, uvar_inputSFTs, &constraints ), &status);
  if ( constraints.detector )
    LALFree ( constraints.detector );

  if ( !FullCatalog || (FullCatalog->length == 0)  )
    {
      XLALPrintError ("\nSorry, didn't find any matching SFTs with pattern '%s'!\n\n", uvar_inputSFTs );
      return CONVERTSFT_EINPUT;
    }

  /* build up full comment-string to be added to SFTs: 1) converted by ConvertToSFTv2, VCS ID 2) user extraComment */
  {
    UINT4 len = 128;
    len += strlen ( uvar_inputSFTs );
    if ( uvar_extraComment )
      len += strlen ( uvar_extraComment );

    if ( ( add_comment = LALCalloc ( 1, len )) == NULL ) {
      XLALPrintError ( "\nOut of memory!\n");
      return CONVERTSFT_EMEM;
    }

    /** \deprecated FIXME: the following code uses obsolete CVS ID tags.
     *  It should be modified to use git version information. */
    sprintf ( add_comment, "Converted by $Id$, inputSFTs = '%s';", uvar_inputSFTs );
    if ( uvar_extraComment )
      {
	strcat ( add_comment, "\n");
	strcat ( add_comment, uvar_extraComment );
      }
  } /* construct comment-string */

  /* which frequency-band to extract? */
  fMin = -1;	/* default: all */
  fMax = -1;
  if ( LALUserVarWasSet ( &uvar_fmin ) )
    fMin = uvar_fmin;
  if ( LALUserVarWasSet ( &uvar_fmax ) )
    fMax = uvar_fmax;

  FILE *fpSingleSFT = NULL;
  if ( uvar_outputSingleSFT )
    XLAL_CHECK ( ( fpSingleSFT = fopen ( uvar_outputSingleSFT, "wb" )) != NULL,
                 XLAL_EIO, "Failed to open singleSFT file '%s' for writing\n", uvar_outputSingleSFT );

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

      if ( ( new_comment  = LALCalloc (1, comment_len )) == NULL ) {
	XLALPrintError ( CONVERTSFT_MSGEMEM );
	return CONVERTSFT_EMEM;
      }
      if ( sft_comment ) {
	strcpy ( new_comment, sft_comment );
	strcat ( new_comment, ";\n");
      }
      strcat ( new_comment, add_comment );

      LAL_CALL ( LALLoadSFTs ( &status, &thisSFT, &oneSFTCatalog, fMin, fMax ), &status );

      if ( uvar_mysteryFactor != 1.0 ) {
	LAL_CALL ( applyFactor2SFTs ( &status, thisSFT, uvar_mysteryFactor ), &status );
      }

      // if user asked for single-SFT output, add this SFT to the open file
      if ( uvar_outputSingleSFT )
        XLAL_CHECK ( XLAL_SUCCESS == XLALWriteSFT2fp( &(thisSFT->data[0]), fpSingleSFT, new_comment ),
                     XLAL_EFUNC,  "XLALWriteSFT2fp() failed to write SFT to '%s'!\n", uvar_outputSingleSFT );

      // if user asked for directory output, write this SFT into that directory
      if ( uvar_outputDir )
        LAL_CALL ( LALWriteSFTVector2Dir (&status, thisSFT, uvar_outputDir, new_comment, uvar_descriptionMisc ), &status );

      LAL_CALL ( LALDestroySFTVector ( &status, &thisSFT ), &status );

      LALFree ( new_comment );

    } /* for i < numSFTs */

  if ( fpSingleSFT ) fclose ( fpSingleSFT );

  /* free memory */
  LALFree ( add_comment );
  LAL_CALL (LALDestroySFTCatalog (&status, &FullCatalog), &status );
  LAL_CALL (LALDestroyUserVars (&status), &status);

  LALCheckMemoryLeaks();

  return 0;
} /* main */


/*----------------------------------------------------------------------*/
/* register all our "user-variables" */
void
initUserVars (LALStatus *status)
{
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* set defaults */
  uvar_outputDir = NULL;
  uvar_outputSingleSFT = NULL;

  uvar_extraComment = NULL;
  uvar_descriptionMisc = NULL;
  uvar_IFO = NULL;

  uvar_minStartTime = 0;
  uvar_maxEndTime = LAL_INT4_MAX;

  uvar_mysteryFactor = 1.0;

  /* now register all our user-variable */
  LALregBOOLUserVar(status,   help,		'h', UVAR_HELP,     "Print this help/usage message");
  LALregSTRINGUserVar(status, inputSFTs,	'i', UVAR_REQUIRED, "File-pattern for input SFTs");
  LALregSTRINGUserVar(status, IFO,		'I', UVAR_OPTIONAL, "IFO of input SFTs: 'G1', 'H1', 'H2', ...(required for v1-SFTs)");

  LALregSTRINGUserVar(status, outputSingleSFT,	'O', UVAR_OPTIONAL, "Output all SFTs into a single concatenated SFT-file with this name");
  LALregSTRINGUserVar(status, outputDir,	'o', UVAR_OPTIONAL, "Output directory for SFTs");

  LALregSTRINGUserVar(status, extraComment,	'C', UVAR_OPTIONAL, "Additional comment to be added to output-SFTs");

  LALregSTRINGUserVar(status, descriptionMisc,	'D', UVAR_OPTIONAL, "'Misc' entry in the SFT filename description-field (see SFTv2 naming convention)");
  LALregREALUserVar(status,   fmin,		'f', UVAR_OPTIONAL, "Lowest frequency to extract from SFTs. [Default: lowest in inputSFTs]");
  LALregREALUserVar(status,   fmax,		'F', UVAR_OPTIONAL, "Highest frequency to extract from SFTs. [Default: highest in inputSFTs]");


  LALregINTUserVar ( status, 	minStartTime, 	 0,  UVAR_OPTIONAL, "Earliest GPS start-time to include");
  LALregINTUserVar ( status, 	maxEndTime, 	 0,  UVAR_OPTIONAL, "Latest GPS end-time to include");


  /* developer-options */
  LALregREALUserVar(status,   mysteryFactor,	 0, UVAR_DEVELOPER, "Change data-normalization by applying this factor (for E@H)");




  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* initUserVars() */

void
applyFactor2SFTs ( LALStatus *status, SFTVector *SFTs, REAL8 factor )
{
  UINT4 i, numSFTs;

  INITSTATUS(status);

  ASSERT ( SFTs, status, CONVERTSFT_EINPUT, CONVERTSFT_MSGEINPUT );

  numSFTs = SFTs->length;

  for ( i=0; i < numSFTs; i ++ )
    {
      SFTtype *thisSFT = &(SFTs->data[i]);
      UINT4 k, numBins = thisSFT->data->length;

      for ( k=0; k < numBins; k ++ )
	{
	  thisSFT->data->data[k].realf_FIXME *= factor;
	  thisSFT->data->data[k].imagf_FIXME *= factor;
	} /* for k < numBins */

    } /* for i < numSFTs */

  RETURN (status);

} /* applyFactor2SFTs() */
