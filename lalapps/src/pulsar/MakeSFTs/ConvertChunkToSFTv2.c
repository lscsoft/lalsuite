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
 * \author Reinhard Prix, Bernd Machenschalk
 * \date 2006
 * \file 
 * \brief Code to convert given input-SFTs (v1 or v2) to v2-SFTs with given extra-comment, 
 *        and write them out following the SFTv2 naming conventions (see LIGO-T040164-01-Z)
 *        This was copied & modified from ConvertToSFTv2.c to not read one SFT at a time,
 *        but now all of themm at once to allow for re-combinig (frequency-) segmented SFTs. 
 * $Id$
 *
 */

/* ---------- includes ---------- */
#include <lalapps.h>

#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>

RCSID ("$Id$");

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
CHAR *uvar_extraComment;
CHAR *uvar_descriptionMisc;
CHAR *uvar_IFO;
REAL8 uvar_fmin;
REAL8 uvar_fmax;

INT4 uvar_minStartTime;
INT4 uvar_maxEndTime;


/*---------- internal prototypes ----------*/
void initUserVars (LALStatus *status);

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
  CHAR *new_comment;
  SFTVector *thisSFT = NULL;
  const CHAR *sft_comment;
  UINT4 comment_len = 0;


  lalDebugLevel = 0;

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;	/* exit with returned status-code on error */
  
  /* set debug level */
  LAL_CALL (LALGetDebugLevel (&status, argc, argv, 'v'), &status);

  /* register all user-variables */
  LAL_CALL (initUserVars (&status), &status);	  

  /* read cmdline & cfgfile  */	
  LAL_CALL (LALUserVarReadAllInput (&status, argc,argv), &status);  

  if (uvar_help) 	/* help requested: we're done */
    exit (0);

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
      LALPrintError ("\nSorry, didn't find any matching SFTs with pattern '%s'!\n\n", uvar_inputSFTs );
      return CONVERTSFT_EINPUT;
    }

  /* build up full comment-string to be added to SFTs: 1) converted by ConvertToSFTv2, RCSID 2) user extraComment */
  {
    UINT4 len = 128;
    len += strlen ( uvar_inputSFTs );
    if ( uvar_extraComment )
      len += strlen ( uvar_extraComment );

    if ( ( add_comment = LALCalloc ( 1, len )) == NULL ) {
      LALPrintError ( "\nOut of memory!\n");
      return CONVERTSFT_EMEM;
    }

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

  comment_len = strlen ( add_comment ) + 10;
  sft_comment = FullCatalog->data->comment;
  if ( sft_comment )
    comment_len += strlen ( sft_comment );
  
  if ( ( new_comment  = LALCalloc (1, comment_len )) == NULL ) {
    LALPrintError ( CONVERTSFT_MSGEMEM );
    return CONVERTSFT_EMEM;
  }
  if ( sft_comment ) {
    strcpy ( new_comment, sft_comment );
    strcat ( new_comment, ";\n");
  }
  strcat ( new_comment, add_comment );

  /* switch this back to LALLoadSFTs once LALLoadSegmntedSFTs has replaced it in lalsupport */
  /* LAL_CALL ( LALLoadSegmentedSFTs ( &status, &thisSFT, FullCatalog, fMin, fMax ), &status ); /**/
  LAL_CALL ( LALLoadSFTs ( &status, &thisSFT, FullCatalog, fMin, fMax ), &status ); /**/

  LAL_CALL ( LALWriteSFTVector2Dir (&status, thisSFT, uvar_outputDir,
				    new_comment, uvar_descriptionMisc ), &status );

  LAL_CALL ( LALDestroySFTVector ( &status, &thisSFT ), &status );

  LALFree ( new_comment );

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
  INITSTATUS( status, "initUserVars", rcsid );
  ATTATCHSTATUSPTR (status);

  /* set defaults */
  uvar_outputDir = LALMalloc (2);
  strcpy ( uvar_outputDir, "." );

  uvar_extraComment = NULL;
  uvar_descriptionMisc = NULL;
  uvar_IFO = NULL;

  uvar_minStartTime = 0;
  uvar_maxEndTime = LAL_INT4_MAX;

  /* now register all our user-variable */
  LALregSTRINGUserVar(status, inputSFTs,	'i', UVAR_REQUIRED, "File-pattern for input SFTs");
  LALregSTRINGUserVar(status, IFO,		'I', UVAR_OPTIONAL, "IFO of input SFTs: 'G1', 'H1', 'H2', ...(required for v1-SFTs)");

  LALregSTRINGUserVar(status, outputDir,	'o', UVAR_OPTIONAL, "Output directory for SFTs");

  LALregSTRINGUserVar(status, extraComment,	'C', UVAR_OPTIONAL, "Additional comment to be added to output-SFTs");

  LALregSTRINGUserVar(status, descriptionMisc,	'D', UVAR_OPTIONAL, "'Misc' entry in the SFT filename description-field (see SFTv2 naming convention)");
  LALregREALUserVar(status,   fmin,		'f', UVAR_OPTIONAL, "Lowest frequency to extract from SFTs. [Default: lowest in inputSFTs]");
  LALregREALUserVar(status,   fmax,		'F', UVAR_OPTIONAL, "Highest frequency to extract from SFTs. [Default: highest in inputSFTs]");


  LALregINTUserVar ( status, 	minStartTime, 	 0,  UVAR_OPTIONAL, "Earliest GPS start-time to include");
  LALregINTUserVar ( status, 	maxEndTime, 	 0,  UVAR_OPTIONAL, "Latest GPS end-time to include");


  LALregBOOLUserVar(status,   help,		'h', UVAR_HELP,     "Print this help/usage message");
  
  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* initUserVars() */
