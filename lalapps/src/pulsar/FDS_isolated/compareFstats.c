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
#define max(x,y) ( (x) > (y) ? (x) : (y) )
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
  UINT4 nlines1, nlines2, minlines;
  BOOLEAN diff_found = 0;
  
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

  nlines1 = Fstats1->lines->nTokens;
  nlines2 = Fstats2->lines->nTokens;

  /* last line now HAS to contain 'DONE'-marker (which is now 7 zeroes..) */
  line1 = Fstats1->lines->tokens[nlines1-1];
  line2 = Fstats2->lines->tokens[nlines2-1];

#define DONE_MARKER "%DONE"
  if ( strcmp(line1, DONE_MARKER ) ) 
    {
      LALPrintError ("ERROR: File '%s' is not properly terminated by '%s' marker!\n\n", uvar_Fname1, DONE_MARKER);
      diff_found=1;
    }
  else
    nlines1 --; /* avoid stepping on DONE-marker in comparison */

  if ( strcmp(line2, DONE_MARKER ) )  
    {
      LALPrintError ("ERROR: File '%s' is not properly terminated by '%s' marker!\n\n", uvar_Fname2, DONE_MARKER);
      diff_found=1;
    }
  else
    nlines2 --;	/* avoid stepping on DONE-marker in comparison */
    
  
  if ( nlines1 != nlines2 )
    {
      LALPrintError ("Fstats files '%s' and '%s' have different length.\n", uvar_Fname1, uvar_Fname2);
      LALPrintError (" len1 = %d, len2 = %d\n\n", nlines1, nlines2);
      diff_found=1;
    }

  /* step through the two files and compare (trying to avoid stumbling on roundoff-errors ) */
  minlines = (nlines1 < nlines2) ? nlines1 : nlines2;
  for (i=0; i < minlines ; i++)
    {
      line1 = Fstats1->lines->tokens[i];
      line2 = Fstats2->lines->tokens[i];

      /* scan ordinary Fstats-lines */
      if ( 7 != sscanf (line1, "%" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT 
			" %" LAL_INT4_FORMAT " %" LAL_REAL4_FORMAT " %" LAL_REAL4_FORMAT " %" LAL_REAL4_FORMAT, 
			&freq1, &a1, &d1, &N1, &mean1, &std1, &max1) ) 
	{
	  LALPrintError ("Failed to parse line %d in file '%s' \n", i+1, uvar_Fname1);
	  exit(-1);
	}
      if ( 7 != sscanf (line2, "%" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT 
			" %" LAL_INT4_FORMAT " %" LAL_REAL4_FORMAT " %" LAL_REAL4_FORMAT " %" LAL_REAL4_FORMAT, 
			&freq2, &a2, &d2, &N2, &mean2, &std2, &max2) ) 
	{
	  LALPrintError ("Failed to parse line %d in file '%s' \n", i+1, uvar_Fname2);
	  exit(-1);
	}



      if ( (relErr=relError( freq1, freq2)) > eps8 )
	{
	  LALPrintError ("Relative frequency-error %g ecceeds %g in line %d\n", relErr, eps8, i+1);
	  diff_found=1;
	} 
      if ( (relErr=relError( a1, a2)) > eps8 )
	{
	  LALPrintError ("Relative error %g in alpha ecceeds %g in line %d\n", relErr, eps8, i+1);
	  diff_found=1;
	} 
      if ( (relErr=relError( d1, d2)) > eps8 )
	{
	  LALPrintError ("Relative error %g in delta ecceeds %g in line %d\n", relErr, eps8, i+1);
	  diff_found=1;
	} 
      if ( N1 != N2)
	{
	  LALPrintError ("Different cluster-sizes in line %d\n", i+1);
	  diff_found=1;
	} 
      if ( (relErr=relError( mean1, mean2)) > eps4 )
	{
	  LALPrintError ("Relative error %g in mean ecceeds %g in line %d\n", relErr, eps4, i+1);
	  diff_found=1;
	}
      if ( (relErr=relError( std1, std2)) > eps4 )
	{
	  LALPrintError ("Relative error %g in std-deviation ecceeds %g in line %d\n", relErr, eps4, i+1);
	  diff_found=1;
	}
      if ( (relErr=relError( max1, max2)) > eps4 )
	{
	  LALPrintError ("Relative error %g in F ecceeds %g in line %d\n", relErr, eps4, i+1);
	  diff_found=1;
	}
      
 
    } /* for i < nLines1 */

  LAL_CALL ( LALDestroyParsedDataFile ( &status, &Fstats1), &status);
  LAL_CALL ( LALDestroyParsedDataFile ( &status, &Fstats2), &status);
  LAL_CALL ( LALDestroyUserVars (&status), &status);

  LALCheckMemoryLeaks(); 

  return diff_found;
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
