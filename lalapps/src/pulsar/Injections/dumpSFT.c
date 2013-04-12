/*
 * Copyright (C) 2004, 2005 Reinhard Prix
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
 * \date 2005
 * \file 
 * \ingroup pulsarApps
 * \brief Code to dump various SFT-info in human-readable form to stdout.
 */

/* ---------- includes ---------- */
#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lalapps.h>

#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include "sft_extra.h"

/** \name Error codes */
/*@{*/
#define MAKEFAKEDATAC_ENORM 	0
#define MAKEFAKEDATAC_ESUB  	1
#define MAKEFAKEDATAC_EARG  	2
#define MAKEFAKEDATAC_EBAD  	3
#define MAKEFAKEDATAC_EFILE 	4
#define MAKEFAKEDATAC_ENOARG 	5


#define MAKEFAKEDATAC_MSGENORM "Normal exit"
#define MAKEFAKEDATAC_MSGESUB  "Subroutine failed"
#define MAKEFAKEDATAC_MSGEARG  "Error parsing arguments"
#define MAKEFAKEDATAC_MSGEBAD  "Bad argument values"
#define MAKEFAKEDATAC_MSGEFILE "File IO error"
#define MAKEFAKEDATAC_MSGENOARG "Missing argument"
/*@}*/


/*---------- DEFINES ----------*/

#define TRUE    (1==1)
#define FALSE   (1==0)

/*----- Macros ----- */
/*---------- internal types ----------*/

/*---------- empty initializers ---------- */
static const LALStatus empty_status;
static const SFTConstraints empty_constraints;
/*---------- Global variables ----------*/

/* User variables */
BOOLEAN uvar_help;
CHAR *uvar_SFTfiles;
BOOLEAN uvar_headerOnly;
BOOLEAN uvar_noHeader;

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
  CHAR detector[2] = "??";	/* allow reading v1-SFTs without detector-info */
  SFTCatalog *catalog = NULL;
  SFTVector *sfts = NULL;
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


  constraints.detector = detector;  /* set '??' detector-constraint, as we don't care about detector-info */
  LAL_CALL ( LALSFTdataFind (&status, &catalog, uvar_SFTfiles, &constraints ), &status );

  if ( !catalog )
    {
      XLALPrintError ("\nNo SFTs seemed to have matched your query!\n\n");
      return 1;
    }

  if ( !uvar_headerOnly ) {
    LAL_CALL ( LALLoadSFTs ( &status, &sfts, catalog, -1, -1 ), &status );
  }

  for ( i=0; i < catalog->length; i ++ )
    {
      SFTDescriptor *ptr = &(catalog->data[i]);

      if ( ! uvar_noHeader )
	{
	  printf ("\n");
	  printf ( "Locator:     '%s'\n", XLALshowSFTLocator ( ptr->locator ) );
	  printf ( "Name:        '%s'\n", ptr->header.name );
	  printf ( "epoch:       [%d, %d]\n", ptr->header.epoch.gpsSeconds, ptr->header.epoch.gpsNanoSeconds ); 
	  printf ( "f0:          %.9f\n", ptr->header.f0 );
	  printf ( "deltaF:      %.9g\n", ptr->header.deltaF );
	  printf ( "comment:     %s\n", (ptr->comment)?(ptr->comment) : "<none>" );
	  printf ( "numBins:     %d\n", ptr->numBins );
	  printf ("\n");
	}

      if ( !uvar_headerOnly )
	{
	  UINT4 k;
	  SFTtype *sft = &(sfts->data[i]);;
	  if ( ! uvar_noHeader ) printf (" Frequency_Hz     Real           Imaginary \n");
	  for ( k=0; k < sft->data->length; k ++ )
	    printf ( "%.9f      % 6e  % 6e  \n", sft->f0 + k * sft->deltaF, crealf(sft->data->data[k]), cimagf(sft->data->data[k]) );
	  
	  printf ("\n");
	  
	} /* if !headerOnly */

    } /* for i < numSFTs */

  if ( sfts ) {
    LAL_CALL ( LALDestroySFTVector (&status, &sfts ), &status );
  }

  /* free memory */
  LAL_CALL (LALDestroySFTCatalog (&status, &catalog), &status );
  LAL_CALL (LALDestroyUserVars (&status), &status);

  LALCheckMemoryLeaks(); 

  return 0;
} /* main */


/*----------------------------------------------------------------------*/
/* register all our "user-variables" */
void
initUserVars (LALStatus *stat)
{
  INITSTATUS(stat);
  ATTATCHSTATUSPTR (stat);

  /* now register all our user-variable */
  LALregBOOLUserVar(stat,   help,	'h', UVAR_HELP,     "Print this help/usage message");
  LALregSTRINGUserVar(stat, SFTfiles,	'i', UVAR_REQUIRED, "File-pattern for input SFTs");
  LALregBOOLUserVar(stat,   headerOnly,	'H', UVAR_OPTIONAL, "Output only header-info");
  LALregBOOLUserVar(stat,   noHeader,	'n', UVAR_OPTIONAL, "Output only data, no header");
  
  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* initUserVars() */

