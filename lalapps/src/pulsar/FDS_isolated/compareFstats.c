/*
*  Copyright (C) 2007 Bernd Machenschalk, Reinhard Prix
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
 * \file
 * \ingroup pulsarApps
 * \author R. Prix
 */

/* ---------- includes ---------- */
#include <math.h>

#include <lalapps.h>

#include <lal/UserInput.h>
#include <lal/LALConstants.h>

/* Error codes and messages */
/**\name Error Codes */ /*@{*/
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

/*@}*/

/***************************************************/
#define TRUE (1==1)
#define FALSE (1==0)

/* (possible) fields of the output Fstat-file */
typedef struct {
  REAL8 Freq;
  REAL8 Alpha;
  REAL8 Delta;
  REAL8 f1dot;
  REAL8 f2dot;
  REAL8 f3dot;
  REAL8 TwoF;
} FstatLine_t;

static const FstatLine_t empty_FstatLine;
/*----------------------------------------------------------------------*/
/* User variables */
CHAR *uvar_Fname1;
CHAR *uvar_Fname2;
BOOLEAN uvar_help;
BOOLEAN uvar_clusterFiles;
REAL8 uvar_Ftolerance;
BOOLEAN uvar_sigFtolerance;
INT4 uvar_Nseg;

#define max(x,y) ( (x) > (y) ? (x) : (y) )

/* local prototypes */
/* Prototypes for the functions defined in this file */
void initUserVars (LALStatus *);
void compareClusterFiles (LALStatus *status, UINT4 *diff, LALParsedDataFile *f1, LALParsedDataFile *f2 );
void compareFstatFiles (LALStatus *status, UINT4 *diff, LALParsedDataFile *f1, LALParsedDataFile *f2,
			REAL8 Ftol);
int parse_Fstat_line ( const CHAR *line, FstatLine_t *FstatLine );
REAL8 relError(REAL8 x, REAL8 y);

/*----------------------------------------------------------------------
 * main function
 *----------------------------------------------------------------------*/
int
main(int argc, char *argv[])
{
  LALStatus status = blank_status;	/* initialize status */
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
      printf ("Fstats files have different length.\n");
      printf (" len1 = %d, len2 = %d\n\n", nlines1, nlines2);
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
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  uvar_clusterFiles = TRUE;	/* default: compare output-files from "cluster" */
  uvar_Ftolerance = 100.0 * LAL_REAL4_EPS;
  uvar_sigFtolerance = FALSE;
  uvar_Nseg = 1;

  /* now register all our user-variable */
  LALregSTRINGUserVar(status, Fname1,	'1', UVAR_REQUIRED, "Path and basefilename for first Fstats file");
  LALregSTRINGUserVar(status, Fname2,	'2', UVAR_REQUIRED, "Path and basefilename for second Fstats file");
  LALregBOOLUserVar(status,   help,	'h', UVAR_HELP,     "Print this help/usage message");

  LALregBOOLUserVar(status,   clusterFiles, 0,UVAR_OPTIONAL,"Comparing cluster results-files or pure Fstat-files");
  LALregREALUserVar(status,   Ftolerance, 0, UVAR_OPTIONAL, "tolerance of error in 2F (relative or sigmas, depending on --sigFtolerance)" );
  LALregBOOLUserVar(status,   sigFtolerance, 0,UVAR_OPTIONAL, "Use error in 2F relative to chi^2 std-deviation 'sigma' instead of relative error");
  LALregINTUserVar(status,    Nseg,           0,UVAR_OPTIONAL, "Number of segments Fstat '2F' is averaged over");

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
  REAL8 eps8 = 100.0 * LAL_REAL8_EPS;	/* tolerances */
  REAL4 eps4 = 1000.0 * LAL_REAL4_EPS;

  INITSTATUS(status);

  nlines1 = f1->lines->nTokens;
  nlines2 = f2->lines->nTokens;

  if ( !nlines1 || !nlines2 )	/* one file is empty.. don't continue */
    return;

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
	  printf ("Failed to parse line %d in file 1\n", i+1);
	  ABORT (status, MAKEFAKEDATAC_EFORMAT, MAKEFAKEDATAC_MSGEFORMAT);
	}
      if ( 7 != sscanf (line2, "%" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT
			" %" LAL_INT4_FORMAT " %" LAL_REAL4_FORMAT
			" %" LAL_REAL4_FORMAT " %" LAL_REAL4_FORMAT,
			&freq2, &a2, &d2, &N2, &mean2, &std2, &Fstat2) )
	{
	  printf ("Failed to parse line %d in file 2\n", i+1);
	  ABORT (status, MAKEFAKEDATAC_EFORMAT, MAKEFAKEDATAC_MSGEFORMAT);
	}

      /* now compare all 7 entries */
      if ( fabs(relErr = relError( freq1, freq2)) > eps8 )
	{
	  printf ("Relative frequency-error %g ecceeds %g in line %d\n", relErr, eps8, i+1);
	  (*diff) ++;
	}
      if ( fabs(relErr = relError( a1, a2)) > eps8 )
	{
	  printf ("Relative error %g in alpha ecceeds %g in line %d\n", relErr, eps8, i+1);
	  (*diff) ++;
	}
      if ( fabs(relErr = relError( d1, d2)) > eps8 )
	{
	  printf ("Relative error %g in delta ecceeds %g in line %d\n", relErr, eps8, i+1);
	  (*diff) ++;
	}
      if ( fabs(relErr = relError( Fstat1, Fstat2)) > eps4 )
	{
	  printf ("Relative error %g in F ecceeds %g in line %d\n", relErr, eps4, i+1);
	  (*diff) ++;
	}
      if ( N1 != N2)
	{
	  printf ("Different cluster-sizes in line %d\n", i+1);
	  (*diff) ++;
	}
      if ( fabs(relErr = relError( mean1, mean2)) > eps4 )
	{
	  printf ("Relative error %g in mean ecceeds %g in line %d\n", relErr, eps4, i+1);
	  (*diff) ++;
	}
      if ( fabs(relErr = relError( std1, std2)) > eps4 )
	{
	  printf ("Relative error %g in std-deviation ecceeds %g in line %d\n",
			 relErr, eps4, i+1);
	  (*diff)++;
	}
    } /* for i < nlines */

  RETURN (status);

} /* compareClusterFiles() */


/** comparison specific to pure Fstat-output files (5 entries )
 */
void
compareFstatFiles (LALStatus *status, UINT4 *diff, LALParsedDataFile *f1, LALParsedDataFile *f2, REAL8 Ftol)
{
  const CHAR *line1, *line2;
  UINT4 nlines1, nlines2, minlines;
  UINT4 i;
  REAL8 relErr;
  REAL4 eps4 = 100.0 * LAL_REAL4_EPS;
  FstatLine_t parsed1 = empty_FstatLine, parsed2 = empty_FstatLine;

  INITSTATUS(status);

  nlines1 = f1->lines->nTokens;
  nlines2 = f2->lines->nTokens;

  /* last line HAS to contain 'DONE'-marker */
  line1 = f1->lines->tokens[nlines1-1];
  line2 = f2->lines->tokens[nlines2-1];

  /* step through the two files and compare (trying to avoid stumbling on roundoff-errors ) */
  minlines = (nlines1 < nlines2) ? nlines1 : nlines2;

  for (i=0; i < minlines ; i++)
    {
      line1 = f1->lines->tokens[i];
      line2 = f2->lines->tokens[i];

      /* read pure Fstats files */
      if ( parse_Fstat_line ( line1, &parsed1 ) ) {
	printf ("Failed to parse line 1 in file 1\n");
	ABORT (status, MAKEFAKEDATAC_EFORMAT, MAKEFAKEDATAC_MSGEFORMAT);
      }
      if ( parse_Fstat_line ( line2, &parsed2 ) ) {
	printf ("Failed to parse line 1 in file 2\n");
	ABORT (status, MAKEFAKEDATAC_EFORMAT, MAKEFAKEDATAC_MSGEFORMAT);
      }

      /* now compare all entries */
      if ( fabs(relErr = relError( parsed1.Freq, parsed2.Freq)) > eps4 )
	{
	  printf ("Relative frequency-error %g ecceeds %g in line %d\n", relErr, eps4, i+1);
	  (*diff) ++;
	}
      if ( fabs(relErr = relError( parsed1.Alpha, parsed2.Alpha)) > eps4 )
	{
	  printf ("Relative error %g in alpha ecceeds %g in line %d\n", relErr, eps4, i+1);
	  (*diff) ++;
	}
      if ( fabs(relErr = relError( parsed1.Delta, parsed2.Delta)) > eps4 )
	{
	  printf ("Relative error %g in delta ecceeds %g in line %d\n", relErr, eps4, i+1);
	  (*diff) ++;
	}
      if ( fabs(relErr = relError( parsed1.f1dot, parsed2.f1dot)) > eps4 )
	{
	  printf ("Relative error %g in f1dot ecceeds %g in line %d\n", relErr, eps4, i+1);
	  (*diff) ++;
	}
      if ( fabs(relErr = relError( parsed1.f2dot, parsed2.f2dot)) > eps4 )
	{
	  printf ("Relative error %g in f2dot ecceeds %g in line %d\n", relErr, eps4, i+1);
	  (*diff) ++;
	}
      if ( fabs(relErr = relError( parsed1.f3dot, parsed2.f3dot)) > eps4 )
	{
	  printf ("Relative error %g in f3dot ecceeds %g in line %d\n", relErr, eps4, i+1);
	  (*diff) ++;
	}
      REAL8 err2F;
      if ( uvar_sigFtolerance )	// measure error in Nseg*2F compared to sigmas of chi^2_(4*Nseg) distribution
        {
          REAL8 mean2F = 0.5 * ( parsed1.TwoF + parsed2.TwoF );
          REAL8 noncent = fmax ( 0, uvar_Nseg * ( mean2F - 4 ) );
          REAL8 sigma = sqrt ( 2 * ( 4*uvar_Nseg + 2 * noncent ) );	// std-dev for noncentral chi^2 distribution with dof degrees of freedom
          err2F = uvar_Nseg * fabs ( parsed1.TwoF - parsed2.TwoF ) / sigma;
        }
      else	// relative error between F1 and F2
        {
          err2F = relError( parsed1.TwoF, parsed2.TwoF );
        }

      if ( err2F > Ftol )
	{
	  printf ("%s Error %g in 2F ecceeds threshold %g in line %d\n", ( uvar_sigFtolerance ? "Sigma" : "Relative" ) , err2F, Ftol, i+1);
	  (*diff) ++;
	}

    } /* for i < minlines */

  RETURN (status);

} /* compareFstatFiles() */

/* parse one Fstat-line into the FstatLine_t struct
 *  This function is flexible concerning the number of spindown-entries found
 *  as CFS_v2 now returns second and third spindown also, while CFS_v1 only
 * has one spindown.
 *
 * return: -1 error, 0 = OK
 */
#define MAX_ENTRIES 7
int
parse_Fstat_line ( const CHAR *line, FstatLine_t *FstatLine )
{
  int ret;
  REAL8 e[MAX_ENTRIES];

  ret = sscanf ( line, "%lf %lf %lf %lf %lf %lf %lf",
		 &e[0], &e[1], &e[2], &e[3], &e[4], &e[5], &e[6] );

  if ( ret < 5 )
    {
      printf("\nFailed to parse Fstat-line (less than 5 entries):\n'%s'\n\n", line );
      return -1;
    }

  if ( ret > 7 )
    {
      printf("\nFailed to parse Fstat-line (more than 7 entries):\n'%s'\n\n", line );
      return -1;
    }

  FstatLine->Freq = e[0];
  FstatLine->Alpha = e[1];
  FstatLine->Delta = e[2];
  FstatLine->f1dot = e[3];

  FstatLine->f2dot = 0;
  FstatLine->f3dot = 0;

  switch ( ret )
    {
    case 5:
      FstatLine->TwoF = e[4];
      break;
    case 6:
      FstatLine->f2dot = e[4];
      FstatLine->TwoF = e[5];
      break;
    case 7:
      FstatLine->f2dot = e[4];
      FstatLine->f3dot = e[5];
      FstatLine->TwoF = e[6];
      break;
    } /* switch(ret) */

  return 0;

} /* parse_Fstat_line() */

REAL8
relError(REAL8 x, REAL8 y)
{
  if ( x == y )
    return 0;

  return ( fabs(x - y )/ (0.5 * ( fabs(x) + fabs(y) ) ) );
} /* relError() */
