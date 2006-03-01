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
 * \brief Code to convert given input-SFTs (v1 or v2) to v2-SFTs with given extra-comment, 
 *        and write them out following the SFTv2 naming conventions (see LIGO-T040164-01-Z)
 *
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

/*----- Macros ----- */
/*---------- internal types ----------*/

/*---------- empty initializers ---------- */
static const LALStatus empty_status;
static const SFTConstraints empty_constraints;
static const SFTCatalog empty_catalog;

/*---------- Global variables ----------*/

/* User variables */
BOOLEAN uvar_help;
CHAR *uvar_inputSFTs;
CHAR *uvar_outputDir;
CHAR *uvar_extraComment;
CHAR *uvar_descriptionMisc;
CHAR *uvar_IFO;

/*---------- internal prototypes ----------*/
void initUserVars (LALStatus *stat);

/*==================== FUNCTION DEFINITIONS ====================*/

/*----------------------------------------------------------------------
 * main function 
 *----------------------------------------------------------------------*/
int
main(int argc, char *argv[]) 
{
  LALStatus status = empty_status;	/* initialize status */
  SFTConstraints constraints = empty_constraints;
  SFTCatalog *FullCatalog = NULL;
  CHAR *add_comment = NULL;
  UINT4 i;

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
  if ( LALUserVarWasSet ( &uvar_IFO ) )
    if ( (constraints.detector = XLALGetChannelPrefix ( uvar_IFO )) == NULL ) {
      return CONVERTSFT_EINPUT;
    }

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
  }

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
	LALPrintError ( CONVERTSFT_MSGEMEM );
	return CONVERTSFT_EMEM;
      }
      if ( sft_comment ) {
	strcpy ( new_comment, sft_comment );
	strcat ( new_comment, ";\n");
      }
      strcat ( new_comment, add_comment );

      LAL_CALL ( LALLoadSFTs ( &status, &thisSFT, &oneSFTCatalog, -1, -1 ), &status );

      LAL_CALL ( LALWriteSFTVector2Dir (&status, thisSFT, uvar_outputDir, new_comment, uvar_descriptionMisc ), &status );

      LAL_CALL ( LALDestroySFTVector ( &status, &thisSFT ), &status );
      
    } /* for i < numSFTs */

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
initUserVars (LALStatus *stat)
{
  INITSTATUS( stat, "initUserVars", rcsid );
  ATTATCHSTATUSPTR (stat);

  /* set defaults */
  uvar_outputDir = LALMalloc (2);
  strcpy ( uvar_outputDir, "." );

  uvar_extraComment = NULL;
  uvar_descriptionMisc = NULL;
  uvar_IFO = NULL;

  /* now register all our user-variable */
  LALregSTRINGUserVar(stat, inputSFTs,		'i', UVAR_REQUIRED, "File-pattern for input SFTs");
  LALregSTRINGUserVar(stat, IFO,		'I', UVAR_REQUIRED, "IFO of input SFTs: 'G1', 'H1', 'H2', ...");

  LALregSTRINGUserVar(stat, outputDir,		'o', UVAR_OPTIONAL, "Output directory for SFTs");

  LALregSTRINGUserVar(stat, extraComment,	'C', UVAR_OPTIONAL, "Additional comment to be added to output-SFTs");

  LALregSTRINGUserVar(stat, descriptionMisc,	'D', UVAR_OPTIONAL, "'Misc' entry in the SFT filename description-field (see SFTv2 naming convention)");
  LALregBOOLUserVar(stat,   help,		'h', UVAR_HELP,     "Print this help/usage message");
  
  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* initUserVars() */

