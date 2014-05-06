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
 * \ingroup pkg_pulsarApps
 * \author R. Prix
 */

/* ---------- includes ---------- */
#include <math.h>

#include <lalapps.h>

#include <lal/UserInput.h>
#include <lal/LALConstants.h>
#include <lal/PulsarDataTypes.h>

// ----- defines & macros ----------
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

/* User variables */
typedef struct
{
  BOOLEAN help;
  CHAR *Fname1;
  CHAR *Fname2;
  BOOLEAN clusterFiles;
  REAL8 Ftolerance;
  BOOLEAN sigFtolerance;
  INT4 Nseg;
} UserVariables_t;

/* ---------- local prototypes ---------- */
int XLALinitUserVars ( UserVariables_t *uvar );
int XLALcompareClusterFiles ( UINT4 *diff, const LALParsedDataFile *f1, const LALParsedDataFile *f2 );
int XLALcompareFstatFiles ( UINT4 *diff, const LALParsedDataFile *f1, const LALParsedDataFile *f2, REAL8 Ftol, const UserVariables_t *uvar );
int XLALParseFstatLine ( FstatLine_t *FstatLine, const CHAR *line );
REAL8 relError ( REAL8 x, REAL8 y );

/*----------------------------------------------------------------------
 * main function
 *----------------------------------------------------------------------*/
int
main (int argc, char *argv[] )
{
  /* register all user-variables */
  UserVariables_t XLAL_INIT_DECL(uvar);
  XLAL_CHECK ( XLALinitUserVars ( &uvar ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* read cmdline & cfgfile  */
  XLAL_CHECK ( XLALUserVarReadAllInput ( argc, argv ) == XLAL_SUCCESS, XLAL_EFUNC );

  if (uvar.help) { 	/* help requested: we're done */
    exit (0);
  }

  /* read in the two Fstats-files (we use XLALParseDataFile() for that purpose) */
  LALParsedDataFile *Fstats1 = NULL, *Fstats2 = NULL;
  XLAL_CHECK ( XLALParseDataFile ( &Fstats1, uvar.Fname1 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALParseDataFile ( &Fstats2, uvar.Fname2 ) == XLAL_SUCCESS, XLAL_EFUNC );

  UINT4 diffs = 0;
  UINT4 nlines1 = Fstats1->lines->nTokens;
  UINT4 nlines2 = Fstats2->lines->nTokens;
  if ( nlines1 != nlines2 )
    {
      XLALPrintError ("Fstats files have different length.\n");
      XLALPrintError (" len1 = %d, len2 = %d\n\n", nlines1, nlines2);
      diffs = 1;
    }

  if ( uvar.clusterFiles )
    {
      XLAL_CHECK ( XLALcompareClusterFiles ( &diffs, Fstats1, Fstats2 ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
  else
    {
      XLAL_CHECK ( XLALcompareFstatFiles ( &diffs, Fstats1, Fstats2, uvar.Ftolerance, &uvar ) == XLAL_SUCCESS, XLAL_EFUNC );
    }

  if ( diffs) {
    XLALPrintError ("FStat files differ! (found %d differences) \n\n", diffs );
  }

  XLALDestroyParsedDataFile ( Fstats1 );
  XLALDestroyParsedDataFile ( Fstats2 );
  XLALDestroyUserVars ();

  LALCheckMemoryLeaks();

  return diffs;

} /* main */


/* register all our "user-variables" */
int
XLALinitUserVars ( UserVariables_t *uvar )
{
  XLAL_CHECK ( uvar != NULL, XLAL_EINVAL );

  uvar->clusterFiles = TRUE;	/* default: compare output-files from "cluster" */
  uvar->Ftolerance = 100.0 * LAL_REAL4_EPS;
  uvar->sigFtolerance = FALSE;
  uvar->Nseg = 1;

  /* now register all user-variables */
  XLALregSTRINGUserStruct ( Fname1,	'1', UVAR_REQUIRED, "Path and basefilename for first Fstats file");
  XLALregSTRINGUserStruct ( Fname2,	'2', UVAR_REQUIRED, "Path and basefilename for second Fstats file");
  XLALregBOOLUserStruct (   help,	'h', UVAR_HELP,     "Print this help/usage message");

  XLALregBOOLUserStruct (   clusterFiles, 0, UVAR_OPTIONAL, "Comparing cluster results-files or pure Fstat-files");
  XLALregREALUserStruct (   Ftolerance,   0, UVAR_OPTIONAL, "tolerance of error in 2F (relative or sigmas, depending on --sigFtolerance)" );
  XLALregBOOLUserStruct (   sigFtolerance, 0,UVAR_OPTIONAL, "Use error in 2F relative to chi^2 std-deviation 'sigma' instead of relative error");
  XLALregINTUserStruct  (   Nseg,          0,UVAR_OPTIONAL, "Number of segments Fstat '2F' is averaged over");

  return XLAL_SUCCESS;

} /* XLALinitUserVars() */


/**
 * comparison specific to cluster-output files (7 entries )
 */
int
XLALcompareClusterFiles ( UINT4 *diff, const LALParsedDataFile *f1, const LALParsedDataFile *f2 )
{
  XLAL_CHECK ( (diff != NULL) && ( f1 != NULL ) && ( f2 != NULL ), XLAL_EINVAL );

  REAL8 eps8 = 100.0 * LAL_REAL8_EPS;	/* tolerances */
  REAL4 eps4 = 1000.0 * LAL_REAL4_EPS;

  UINT4 nlines1 = f1->lines->nTokens;
  UINT4 nlines2 = f2->lines->nTokens;

  /* step through the two files and compare (trying to avoid stumbling on roundoff-errors ) */
  UINT4 minlines = (nlines1 < nlines2) ? nlines1 : nlines2;

  for ( UINT4 i=0; i < minlines ; i++ )
    {
      const char *line1 = f1->lines->tokens[i];
      const char *line2 = f2->lines->tokens[i];

      /* scan Fstats-lines of cluster-output */
      REAL8 freq1, freq2, a1, a2, d1, d2;
      REAL4 mean1, mean2, std1, std2;
      INT4 N1, N2;
      REAL4 Fstat1, Fstat2;

      int numRead = sscanf (line1, "%" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT
                            " %" LAL_INT4_FORMAT " %" LAL_REAL4_FORMAT
                            " %" LAL_REAL4_FORMAT " %" LAL_REAL4_FORMAT,
                            &freq1, &a1, &d1, &N1, &mean1, &std1, &Fstat1 );
      XLAL_CHECK ( numRead == 7, XLAL_EDATA, "Failed to parse line %d in file 1\n", i+1 );

      numRead = sscanf (line2, "%" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT " %" LAL_REAL8_FORMAT
			" %" LAL_INT4_FORMAT " %" LAL_REAL4_FORMAT
			" %" LAL_REAL4_FORMAT " %" LAL_REAL4_FORMAT,
			&freq2, &a2, &d2, &N2, &mean2, &std2, &Fstat2 );
      XLAL_CHECK ( numRead == 7, XLAL_EDATA, "Failed to parse line %d in file 2\n", i+1 );

      /* now compare all 7 entries */
      REAL8 relErr;
      if ( (relErr = relError( freq1, freq2)) > eps8 )
	{
	  XLALPrintError ( "Relative frequency-error %g ecceeds %g in line %d\n", relErr, eps8, i+1);
	  (*diff) ++;
	}
      if ( (relErr = relError( a1, a2)) > eps8 )
	{
	  XLALPrintError ("Relative error %g in alpha ecceeds %g in line %d\n", relErr, eps8, i+1);
	  (*diff) ++;
	}
      if ( (relErr = relError( d1, d2)) > eps8 )
	{
	  XLALPrintError ("Relative error %g in delta ecceeds %g in line %d\n", relErr, eps8, i+1);
	  (*diff) ++;
	}
      if ( (relErr = relError( Fstat1, Fstat2)) > eps4 )
	{
	  XLALPrintError ("Relative error %g in F ecceeds %g in line %d\n", relErr, eps4, i+1);
	  (*diff) ++;
	}
      if ( N1 != N2)
	{
	  XLALPrintError ("Different cluster-sizes in line %d\n", i+1 );
	  (*diff) ++;
	}
      if ( (relErr = relError( mean1, mean2)) > eps4 )
	{
	  XLALPrintError ("Relative error %g in mean ecceeds %g in line %d\n", relErr, eps4, i+1);
	  (*diff) ++;
	}
      if ( (relErr = relError( std1, std2)) > eps4 )
	{
	  XLALPrintError ("Relative error %g in std-deviation ecceeds %g in line %d\n", relErr, eps4, i+1);
	  (*diff)++;
	}
    } /* for i < nlines */

  return XLAL_SUCCESS;

} /* XLALcompareClusterFiles() */


/**
 * comparison specific to pure Fstat-output files (5 entries )
 */
int
XLALcompareFstatFiles ( UINT4 *diff, const LALParsedDataFile *f1, const LALParsedDataFile *f2, REAL8 Ftol, const UserVariables_t *uvar )
{
  XLAL_CHECK ( (diff != NULL) && (f1 != NULL) && ( f2 != NULL ) && (Ftol > 0), XLAL_EINVAL );

  REAL4 eps4 = 100.0 * LAL_REAL4_EPS;
  FstatLine_t XLAL_INIT_DECL(parsed1);
  FstatLine_t XLAL_INIT_DECL(parsed2);

  UINT4 nlines1 = f1->lines->nTokens;
  UINT4 nlines2 = f2->lines->nTokens;

  /* step through the two files and compare (trying to avoid stumbling on roundoff-errors ) */
  UINT4 minlines = (nlines1 < nlines2) ? nlines1 : nlines2;

  for (UINT4 i=0; i < minlines ; i++)
    {
      const char *line1 = f1->lines->tokens[i];
      const char *line2 = f2->lines->tokens[i];

      /* read pure Fstats files */
      XLAL_CHECK ( XLALParseFstatLine ( &parsed1, line1 ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK ( XLALParseFstatLine ( &parsed2, line2 ) == XLAL_SUCCESS, XLAL_EFUNC );

      /* now compare all entries */
      REAL8 relErr;
      if ( (relErr = relError( parsed1.Freq, parsed2.Freq)) > eps4 )
	{
	  XLALPrintError ("Relative frequency-error %g ecceeds %g in line %d\n", relErr, eps4, i+1);
	  (*diff) ++;
	}
      if ( (relErr = relError( parsed1.Alpha, parsed2.Alpha)) > eps4 )
	{
	  XLALPrintError ("Relative error %g in alpha ecceeds %g in line %d\n", relErr, eps4, i+1);
	  (*diff) ++;
	}
      if ( (relErr = relError( parsed1.Delta, parsed2.Delta)) > eps4 )
	{
	  XLALPrintError ("Relative error %g in delta ecceeds %g in line %d\n", relErr, eps4, i+1);
	  (*diff) ++;
	}
      if ( (relErr = relError( parsed1.f1dot, parsed2.f1dot)) > eps4 )
	{
	  XLALPrintError ("Relative error %g in f1dot ecceeds %g in line %d\n", relErr, eps4, i+1);
	  (*diff) ++;
	}
      if ( (relErr = relError( parsed1.f2dot, parsed2.f2dot)) > eps4 )
	{
	  XLALPrintError ("Relative error %g in f2dot ecceeds %g in line %d\n", relErr, eps4, i+1);
	  (*diff) ++;
	}
      if ( (relErr = relError( parsed1.f3dot, parsed2.f3dot)) > eps4 )
	{
	  XLALPrintError ("Relative error %g in f3dot ecceeds %g in line %d\n", relErr, eps4, i+1);
	  (*diff) ++;
	}

      REAL8 err2F;
      if ( uvar->sigFtolerance )	// measure error in Nseg*2F compared to sigmas of chi^2_(4*Nseg) distribution
        {
          REAL8 mean2F = fmax ( 4, 0.5 * ( parsed1.TwoF + parsed2.TwoF ) );
          REAL8 noncent = uvar->Nseg * ( mean2F - 4 );
          REAL8 sigma = sqrt ( 2 * ( 4*uvar->Nseg + 2 * noncent ) );	// std-dev for noncentral chi^2 distribution with dof degrees of freedom
          err2F = uvar->Nseg * fabs ( parsed1.TwoF - parsed2.TwoF ) / sigma;
        }
      else	// relative error between F1 and F2
        {
          err2F = relError( parsed1.TwoF, parsed2.TwoF );
        }

      if ( err2F > Ftol )
	{
	  XLALPrintError ("%s Error %g in 2F ecceeds threshold %g in line %d\n", ( uvar->sigFtolerance ? "Sigma" : "Relative" ) , err2F, Ftol, i+1);
	  (*diff) ++;
	}

    } /* for i < minlines */

  return XLAL_SUCCESS;

} /* XLALcompareFstatFiles() */

/* parse one Fstat-line into the FstatLine_t struct
 * This function is flexible concerning the number of spindown-entries found
 * as CFS_v2 now returns second and third spindown also, while CFS_v1 only
 * has one spindown.
 *
 */
int
XLALParseFstatLine ( FstatLine_t *FstatLine, const CHAR *line )
{
  XLAL_CHECK ( (line != NULL) && (FstatLine != NULL), XLAL_EINVAL );


  REAL8 e[7];
  int ret = sscanf ( line, "%lf %lf %lf %lf %lf %lf %lf",
                     &e[0], &e[1], &e[2], &e[3], &e[4], &e[5], &e[6] );

  XLAL_CHECK ( ret >= 5, XLAL_EDATA, "Failed to parse Fstat-line (less than 5 entries):\n'%s'\n", line );

  FstatLine->Freq  = e[0];
  FstatLine->Alpha = e[1];
  FstatLine->Delta = e[2];
  FstatLine->f1dot = e[3];


  switch ( ret )
    {
    case 5:
      FstatLine->TwoF = e[4];
      FstatLine->f2dot = 0;
      FstatLine->f3dot = 0;
      break;
    case 6:
      FstatLine->f2dot = e[4];
      FstatLine->TwoF  = e[5];
      FstatLine->f3dot = 0;
      break;
    case 7:
      FstatLine->f2dot = e[4];
      FstatLine->f3dot = e[5];
      FstatLine->TwoF  = e[6];
      break;

    } /* switch(ret) */

  return XLAL_SUCCESS;

} /* XLALParseFstatLine() */

REAL8
relError(REAL8 x, REAL8 y)
{
  if ( x == y ) {
    return 0;
  }

  REAL8 denom = fmax ( 4, 0.5*(x+y) );
  return fabs ( (x - y ) / denom );

} /* relError() */
