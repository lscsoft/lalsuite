 /*-----------------------------------------------------------------------
 *
 * File Name: dumpSFT.c
 *
 * Authors: Prix, R.
 *
 * Revision: $Id$
 *           
 *-----------------------------------------------------------------------
 */

/* ---------- includes ---------- */
#include <lalapps.h>

#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include "sft_extra.h"

RCSID ("$Id");

/* Error codes and messages */
/* <lalErrTable file="MAKEFAKEDATACErrorTable"> */
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


/* </lalErrTable> */

/***************************************************/
#define TRUE (1==1)
#define FALSE (1==0)

/* local prototypes */
/* Prototypes for the functions defined in this file */
void initUserVars (LALStatus *stat);

/*----------------------------------------------------------------------*/
static const LALStatus empty_status;
/*----------------------------------------------------------------------*/
/* User variables */
CHAR *uvar_SFTfname;
CHAR *uvar_SFToutname;
BOOLEAN uvar_help;
INT4 uvar_format;
/*----------------------------------------------------------------------
 * main function 
 *----------------------------------------------------------------------*/
int
main(int argc, char *argv[]) 
{
  LALStatus status = empty_status;	/* initialize status */
  SFTtype *sft = NULL;
  SFTHeader header;
  REAL8 deltaF, f_min, f_max;
  FILE *fp;

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

  /* figure out what frequency-band these SFTs contain */
  LAL_CALL (LALReadSFTheader (&status, &header, uvar_SFTfname), &status);

  deltaF = 1.0 / header.timeBase;
  f_min = header.fminBinIndex * deltaF;
  f_max = f_min + (header.length - 1) * deltaF;

  /* read in the binary sft-file */
  LAL_CALL (LALReadSFTfile (&status, &sft, f_min, f_max, uvar_SFTfname), &status);

  /* now dump it into a text-file or to stdout */
  if ( uvar_SFToutname && LALUserVarWasSet (&uvar_SFToutname) )
    {
      if ( (fp = fopen (uvar_SFToutname, "w")) == NULL) 
	{
	  LALPrintError ("Could not open file `%s` for writing.\n", uvar_SFToutname);
	  exit (-1);
	}
    }
  else
    fp = stdout;
  
  dump_SFT (fp, sft, uvar_format);

  if (fp != stdout)
    fclose (fp);


  /* free memory */
  LAL_CALL (LALDestroySFTtype (&status, &sft), &status);
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

#define DEFAULT_OUTFILE "stdout"
  uvar_SFToutname = LALCalloc (strlen(DEFAULT_OUTFILE)+1, 1);
  strcpy (uvar_SFToutname, DEFAULT_OUTFILE);

  uvar_format = 0;	/* defaults to openDX */

  /* now register all our user-variable */
  LALregSTRINGUserVar(stat, SFTfname,	'i', UVAR_REQUIRED, "Path and filename for binary SFTs-file");
  LALregSTRINGUserVar(stat, SFToutname,	'o', UVAR_OPTIONAL, "Path and filename for output SFT text-file");
  LALregINTUserVar(stat,    format,	'F', UVAR_OPTIONAL, "Output format: 0 = openDX, 1 = xmgrace");
  LALregBOOLUserVar(stat,   help,	'h', UVAR_HELP,     "Print this help/usage message");
  
  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* initUserVars() */

