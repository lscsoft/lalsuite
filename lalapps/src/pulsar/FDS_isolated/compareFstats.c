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

RCSID ("$Id$");

/* Error codes and messages */
/* <lalErrTable file="MAKEFAKEDATACErrorTable"> */
#define MAKEFAKEDATAC_ENORM 	0
#define MAKEFAKEDATAC_ESUB  	1
#define MAKEFAKEDATAC_EARG  	2
#define MAKEFAKEDATAC_EBAD  	3
#define MAKEFAKEDATAC_EFILE 	4
#define MAKEFAKEDATAC_ENOARG 	5
#define MAKEFAKEDATAC_EFORMAT	6

#define MAKEFAKEDATAC_MSGENORM "Normal exit"
#define MAKEFAKEDATAC_MSGESUB  "Subroutine failed"
#define MAKEFAKEDATAC_MSGEARG  "Error parsing arguments"
#define MAKEFAKEDATAC_MSGEBAD  "Bad argument values"
#define MAKEFAKEDATAC_MSGEFILE "File IO error"
#define MAKEFAKEDATAC_MSGENOARG "Missing argument"
#define MAKEFAKEDATAC_MSGEFORMAT "Bad file format"

/* </lalErrTable> */

/***************************************************/
#define TRUE (1==1)
#define FALSE (1==0)

/* local prototypes */
/* Prototypes for the functions defined in this file */
void initUserVars (LALStatus *);
void compareClusterFiles (LALStatus *status, UINT4 *diff, LALParsedDataFile *f1, LALParsedDataFile *f2 );
void compareFstatFiles (LALStatus *status, UINT4 *diff, LALParsedDataFile *f1, LALParsedDataFile *f2,
			REAL8 Ftol);

/*----------------------------------------------------------------------*/
static const LALStatus empty_status;
/*----------------------------------------------------------------------*/
/* User variables */
CHAR *uvar_Fname1;
CHAR *uvar_Fname2;
BOOLEAN uvar_help;
BOOLEAN uvar_clusterFiles;
REAL8 uvar_Ftolerance;

#define max(x,y) ( (x) > (y) ? (x) : (y) )
#define relError(x,y) (fabs( (x)-(y) )/ (max(fabs(x), fabs(y))))

/*----------------------------------------------------------------------
 * main function 
 *----------------------------------------------------------------------*/
int
main(int argc, char *argv[]) 
{
  LALStatus status = empty_status;	/* initialize status */
  LALParsedDataFile *Fstats1 =NULL, *Fstats2 = NULL;
  UINT4 diffs = 0;
  UINT4 nlines1, nlines2;

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
  if ( nlines1 != nlines2 )
    {
      LALPrintError ("Fstats files have different length.\n");
      LALPrintError (" len1 = %d, len2 = %d\n\n", nlines1, nlines2);
      diffs = 1;
    }

  if ( uvar_clusterFiles ) {
    LAL_CALL ( compareClusterFiles ( &status, &diffs, Fstats1, Fstats2 ), &status );
  } else {
    LAL_CALL ( compareFstatFiles ( &status, &diffs, Fstats1, Fstats2, uvar_Ftolerance ), 
	       &status );
  }

  if ( diffs)
    fprintf(stderr,"\nFStat files differ! (found %d differences) \n\n", diffs);
  
  LAL_CALL ( LALDestroyParsedDataFile ( &status, &Fstats1), &status);
  LAL_CALL ( LALDestroyParsedDataFile ( &status, &Fstats2), &status);
  LAL_CALL ( LALDestroyUserVars (&status), &status);
  
  LALCheckMemoryLeaks(); 
  
  return diffs;
} /* main */


/*----------------------------------------------------------------------*/
/* register all our "user-variables" */
void
initUserVars (LALStatus *status)
{
  INITSTATUS( status, "initUserVars", rcsid );
  ATTATCHSTATUSPTR (status);

  uvar_clusterFiles = TRUE;	/* default: compare output-files from "cluster" */
  uvar_Ftolerance = 100.0 * LAL_REAL4_EPS;

  /* now register all our user-variable */
  LALregSTRINGUserVar(status, Fname1,	'1', UVAR_REQUIRED, "Path and basefilename for first Fstats file");
  LALregSTRINGUserVar(status, Fname2,	'2', UVAR_REQUIRED, "Path and basefilename for second Fstats file");
  LALregBOOLUserVar(status,   help,	'h', UVAR_HELP,     "Print this help/usage message");

  LALregBOOLUserVar(status,   clusterFiles, 0,UVAR_OPTIONAL,"Comparing cluster results-files or pure Fstat-files"); 
  LALregREALUserVar(status,   Ftolerance, 0, UVAR_OPTIONAL, "tolerance of relative-error in F" );

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* initUserVars() */


/** comparison specific to cluster-output files (7 entries ) 
 */    
void
compareClusterFiles (LALStatus *status, UINT4 *diff, LALParsedDataFile *f1, LALParsedDataFile *f2 )
{
  const CHAR *line1, *line2;
  UINT4 nlines1, nlines2, minlines;
  UINT4 i;
  REAL8 freq1, freq2, a1, a2, d1, d2;
  REAL4 mean1, mean2, std1, std2;
  REAL4 Fstat1, Fstat2;
  INT4 N1, N2;
  REAL8 relErr;
  REAL8 eps8 = 10.0 * LAL_REAL8_EPS;	/* tolerances */
  REAL4 eps4 = 100.0 * LAL_REAL4_EPS;

  INITSTATUS (status, "compareClusterFiles", rcsid );

  nlines1 = f1->lines->nTokens;
  nlines2 = f2->lines->nTokens;

  if ( !nlines1 || !nlines2 )	/* one file is empty.. don't continue */
    return;

  /* cluster-files: last line now HAS to contain 'DONE'-marker */
  line1 = f1->lines->tokens[nlines1-1];
  line2 = f2->lines->tokens[nlines2-1];

#define DONE_MARKER "%DONE"
  if ( strcmp(line1, DONE_MARKER ) ) 
    {
      LALPrintError ("ERROR: File 1 is not properly terminated by '%s' marker!\n\n", 
		     DONE_MARKER);
      (*diff) ++;	/* increase difference-counter */
    }
  else
    nlines1 --; /* avoid stepping on DONE-marker in comparison */

  if ( strcmp(line2, DONE_MARKER ) )  
    {
      LALPrintError ("ERROR: File 2 is not properly terminated by '%s' marker!\n\n", 
		     DONE_MARKER);
      (*diff) ++;
    }
  else
    nlines2 --;	/* avoid stepping on DONE-marker in comparison */


  /* step through the two files and compare (trying to avoid stumbling on roundoff-errors ) */
  minlines = (nlines1 < nlines2) ? nlines1 : nlines2;

  for (i=0; i < minlines ; i++)
    {
      line1 = f1->lines->tokens[i];
      line2 = f2->lines->tokens[i];

      /* scan Fstats-lines of cluster-output */
      if ( 7 != sscanf (line1, "%" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT 
			" %" LAL_INT4_FORMAT " %" LAL_REAL4_FORMAT 
			" %" LAL_REAL4_FORMAT " %" LAL_REAL4_FORMAT, 
			&freq1, &a1, &d1, &N1, &mean1, &std1, &Fstat1) ) 
	{
	  LALPrintError ("Failed to parse line %d in file 1\n", i+1);
	  ABORT (status, MAKEFAKEDATAC_EFORMAT, MAKEFAKEDATAC_MSGEFORMAT);
	}
      if ( 7 != sscanf (line2, "%" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT 
			" %" LAL_INT4_FORMAT " %" LAL_REAL4_FORMAT 
			" %" LAL_REAL4_FORMAT " %" LAL_REAL4_FORMAT, 
			&freq2, &a2, &d2, &N2, &mean2, &std2, &Fstat2) ) 
	{
	  LALPrintError ("Failed to parse line %d in file 2\n", i+1);
	  ABORT (status, MAKEFAKEDATAC_EFORMAT, MAKEFAKEDATAC_MSGEFORMAT);
	}
      
      /* now compare all 7 entries */
      if ( (relErr = relError( freq1, freq2)) > eps8 )
	{
	  LALPrintError ("Relative frequency-error %g ecceeds %g in line %d\n", relErr, eps8, i+1);
	  (*diff) ++;
	} 
      if ( (relErr = relError( a1, a2)) > eps8 )
	{
	  LALPrintError ("Relative error %g in alpha ecceeds %g in line %d\n", relErr, eps8, i+1);
	  (*diff) ++;
	} 
      if ( (relErr = relError( d1, d2)) > eps8 )
	{
	  LALPrintError ("Relative error %g in delta ecceeds %g in line %d\n", relErr, eps8, i+1);
	  (*diff) ++;
	} 
      if ( (relErr = relError( Fstat1, Fstat2)) > eps4 )
	{
	  LALPrintError ("Relative error %g in F ecceeds %g in line %d\n", relErr, eps4, i+1);
	  (*diff) ++;
	}
      if ( N1 != N2)
	{
	  LALPrintError ("Different cluster-sizes in line %d\n", i+1);
	  (*diff) ++;
	} 
      if ( (relErr = relError( mean1, mean2)) > eps4 )
	{
	  LALPrintError ("Relative error %g in mean ecceeds %g in line %d\n", relErr, eps4, i+1);
	  (*diff) ++;
	}
      if ( (relErr = relError( std1, std2)) > eps4 )
	{
	  LALPrintError ("Relative error %g in std-deviation ecceeds %g in line %d\n", 
			 relErr, eps4, i+1);
	  (*diff)++;
	}
    } /* for i < nlines */

  RETURN (status);

} /* compareClusterFiles() */


/** comparison specific to pure Fstat-output files (5 entries ) 
 */    
void
compareFstatFiles (LALStatus *status, UINT4 *diff, LALParsedDataFile *f1, LALParsedDataFile *f2, 
		   REAL8 Ftol)
{
  const CHAR *line1, *line2;
  UINT4 nlines1, nlines2, minlines;
  UINT4 i;
  REAL8 freq1, freq2, a1, a2, d1, d2;
  REAL4 Fstat1, Fstat2;
  REAL8 relErr;
  REAL4 f1dot1, f1dot2;
  REAL4 eps4 = 100.0 * LAL_REAL4_EPS;

  INITSTATUS (status, "compareFstatFiles", rcsid );

  nlines1 = f1->lines->nTokens;
  nlines2 = f2->lines->nTokens;

  /* step through the two files and compare (trying to avoid stumbling on roundoff-errors ) */
  minlines = (nlines1 < nlines2) ? nlines1 : nlines2;

  for (i=0; i < minlines ; i++)
    {
      line1 = f1->lines->tokens[i];
      line2 = f2->lines->tokens[i];

      /* read pure Fstats files */
      if ( 5 != sscanf (line1, "%" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT 
			" %" LAL_REAL4_FORMAT " %" LAL_REAL4_FORMAT, 
			&freq1, &a1, &d1, &f1dot1, &Fstat1) ) 
	{
	  LALPrintError ("Failed to parse line %d in file 1\n", i+1);
	  ABORT (status, MAKEFAKEDATAC_EFORMAT, MAKEFAKEDATAC_MSGEFORMAT);
	}

      if ( 5 != sscanf (line2, "%" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT 
			" %" LAL_REAL4_FORMAT " %" LAL_REAL4_FORMAT, 
			&freq2, &a2, &d2, &f1dot2, &Fstat2) ) 
	{
	  LALPrintError ("Failed to parse line %d in file 2\n", i+1);
	  ABORT (status, MAKEFAKEDATAC_EFORMAT, MAKEFAKEDATAC_MSGEFORMAT);
	}

      /* now compare all 5 entries */
      if ( (relErr = relError( freq1, freq2)) > eps4 )
	{
	  LALPrintError ("Relative frequency-error %g ecceeds %g in line %d\n", relErr, eps4, i+1);
	  (*diff) ++;
	} 
      if ( (relErr = relError( a1, a2)) > eps4 )
	{
	  LALPrintError ("Relative error %g in alpha ecceeds %g in line %d\n", relErr, eps4, i+1);
	  (*diff) ++;
	} 
      if ( (relErr = relError( d1, d2)) > eps4 )
	{
	  LALPrintError ("Relative error %g in delta ecceeds %g in line %d\n", relErr, eps4, i+1);
	  (*diff) ++;
	}
      if ( (relErr = relError( f1dot1, f1dot2)) > eps4 )
	{
	  LALPrintError ("Relative error %g in f1dot ecceeds %g in line %d\n", relErr, eps4, i+1);
	  (*diff) ++;
	} 
      if ( (relErr = relError( Fstat1, Fstat2)) > Ftol )
	{
	  LALPrintError ("Relative error %g in F ecceeds %g in line %d\n", relErr, Ftol, i+1);
	  (*diff) ++;
	}

    } /* for i < minlines */

  RETURN (status);

} /* compareFstatFiles() */


