 /*-----------------------------------------------------------------------
 *
 * File Name: compareFstats.c
 *
 * Authors: Prix, R.
 *
 * Revision: $Id$
 *           
 *-----------------------------------------------------------------------
 */

/* ---------- includes ---------- */
#include <math.h>

#include <lalapps.h>

#include <lal/UserInput.h>
#include <lal/LALConstants.h>

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
CHAR *uvar_Fname1;
CHAR *uvar_Fname2;
BOOLEAN uvar_help;

#define relError(x,y) (fabs( (x)-(y) )/(x))

/*----------------------------------------------------------------------
 * main function 
 *----------------------------------------------------------------------*/
int
main(int argc, char *argv[]) 
{
  LALStatus status = empty_status;	/* initialize status */
  LALParsedDataFile *Fstats1 =NULL, *Fstats2 = NULL;
  const CHAR *line1, *line2;
  REAL8 freq1, freq2, a1, a2, d1, d2;
  REAL4 mean1, mean2, std1, std2, max1, max2;
  INT4 N1, N2;
  UINT4 i;
  REAL8 eps8 = 10.0 * LAL_REAL8_EPS;
  REAL4 eps4 = 100.0 * LAL_REAL4_EPS;
  REAL8 relErr;

  lalDebugLevel = 0;

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;	/* exit with returned status-code on error */

  /* register all user-variables */
  LAL_CALL (LALGetDebugLevel (&status, argc, argv, 'v'), &status);
  LAL_CALL (initUserVars (&status), &status);	  

  /* read cmdline & cfgfile  */	
  LAL_CALL (LALUserVarReadAllInput (&status, argc,argv), &status);  

  if (uvar_help) 	/* help requested: we're done */
    exit (0);

  /* read in the two Fstats-files (we use LALParseDataFile() for that purpose) */
  LAL_CALL (LALParseDataFile (&status, &Fstats1, uvar_Fname1), &status);
  LAL_CALL (LALParseDataFile (&status, &Fstats2, uvar_Fname2), &status);

  if ( Fstats2->lines->nTokens != Fstats2->lines->nTokens )
    {
      LALPrintError ("\nFstats files '%s' and '%s' have different length\n", uvar_Fname1, uvar_Fname2);
      exit(1);
    }

  /* step through the two files and compare (trying to avoid stumbling on roundoff-errors ) */
  for (i=0; i < Fstats1->lines->nTokens; i++)
    {
      line1 = Fstats1->lines->tokens[i];
      line2 = Fstats2->lines->tokens[i];
      if ( 7 != sscanf (line1, "%" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT 
			" %" LAL_INT4_FORMAT " %" LAL_REAL4_FORMAT " %" LAL_REAL4_FORMAT " %" LAL_REAL4_FORMAT, 
			&freq1, &a1, &d1, &N1, &mean1, &std1, &max1) ) 
	{
	  LALPrintError ("\nFailed to parse line %d in file '%s' \n", i+1, uvar_Fname1);
	  exit(-1);
	}
      if ( 7 != sscanf (line2, "%" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT 
			" %" LAL_INT4_FORMAT " %" LAL_REAL4_FORMAT " %" LAL_REAL4_FORMAT " %" LAL_REAL4_FORMAT, 
			&freq2, &a2, &d2, &N2, &mean2, &std2, &max2) ) 
	{
	  LALPrintError ("\nFailed to parse line %d in file '%s' \n", i+1, uvar_Fname2);
	  exit(-1);
	}



      if ( (relErr=relError( freq1, freq2)) > eps8 )
	{
	  LALPrintError ("\nRelative frequency-error %g ecceeds %g in line %d\n", relErr, eps8, i+1);
	  exit (1);
	} 
      if ( (relErr=relError( a1, a2)) > eps8 )
	{
	  LALPrintError ("\nRelative error %g in alpha ecceeds %g in line %d\n", relErr, eps8, i+1);
	  exit (1);
	} 
      if ( (relErr=relError( d1, d2)) > eps8 )
	{
	  LALPrintError ("\nRelative error %g in delta ecceeds %g in line %d\n", relErr, eps8, i+1);
	  exit (1);
	} 
      if ( N1 != N2)
	{
	  LALPrintError ("\nDifferent cluster-sizes in line %d\n", i+1);
	  exit (1);
	} 
      if ( (relErr=relError( mean1, mean2)) > eps4 )
	{
	  LALPrintError ("\nRelative error %g in mean ecceeds %g in line %d\n", relErr, eps4, i+1);
	  exit (1);
	}
      if ( (relErr=relError( std1, std2)) > eps4 )
	{
	  LALPrintError ("\nRelative error %g in std-deviation ecceeds %g in line %d\n", relErr, eps4, i+1);
	  exit (1);
	}
      if ( (relErr=relError( max1, max2)) > eps4 )
	{
	  LALPrintError ("\nRelative error %g in std-deviation ecceeds %g in line %d\n", relErr, eps4, i+1);
	  exit (1);
	}
      
 
    } /* for i < nLines1 */

  LAL_CALL ( LALDestroyParsedDataFile ( &status, &Fstats1), &status);
  LAL_CALL ( LALDestroyParsedDataFile ( &status, &Fstats2), &status);
  LAL_CALL ( LALDestroyUserVars (&status), &status);

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

  /* now register all our user-variable */

  LALregSTRINGUserVar(stat, Fname1,	'1', UVAR_REQUIRED, "Path and basefilename for first Fstats file");
  LALregSTRINGUserVar(stat, Fname2,	'2', UVAR_REQUIRED, "Path and basefilename for second Fstats file");
  LALregBOOLUserVar(stat,   help,	'h', UVAR_HELP,     "Print this help/usage message");

  DETATCHSTATUSPTR (stat);
  RETURN (stat);

} /* initUserVars() */
