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
#include <lal/LFTandTSutils.h>

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

  REAL8 tol_L1;		// tolerance on relative error between vectors using L1 norm
  REAL8 tol_L2;		// tolerance on relative error between vectors using L2 norm
  REAL8 tol_angle;	// tolerance on angle between the two vectors, in radians
  REAL8 tol_atMax;	// tolerance on single-sample relative error *at* respective maximum

} UserVariables_t;

/* ---------- local prototypes ---------- */
int XLALinitUserVars ( UserVariables_t *uvar );
int XLALcompareFstatFiles ( const LALParsedDataFile *f1, const LALParsedDataFile *f2, VectorComparison tol );
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

  VectorComparison XLAL_INIT_DECL(tol);
  tol.relErr_L1 	= uvar.tol_L1;
  tol.relErr_L2		= uvar.tol_L2;
  tol.angleV 		= uvar.tol_angle;
  tol.relErr_atMaxAbsx 	= uvar.tol_atMax;
  tol.relErr_atMaxAbsy  = uvar.tol_atMax;

  XLAL_CHECK ( XLALcompareFstatFiles ( Fstats1, Fstats2, tol ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLALDestroyParsedDataFile ( Fstats1 );
  XLALDestroyParsedDataFile ( Fstats2 );
  XLALDestroyUserVars ();

  LALCheckMemoryLeaks();

  return XLAL_SUCCESS;

} /* main */


/* register all our "user-variables" */
int
XLALinitUserVars ( UserVariables_t *uvar )
{
  XLAL_CHECK ( uvar != NULL, XLAL_EINVAL );

  uvar->tol_L1 		= 5.5e-2;
  uvar->tol_L2 		= 4.5e-2;
  uvar->tol_angle	= 0.04;  // rad
  uvar->tol_atMax 	= 5e-2;

  /* now register all user-variables */
  XLALregSTRINGUserStruct ( Fname1,	'1', UVAR_REQUIRED, "Path and basefilename for first Fstats file");
  XLALregSTRINGUserStruct ( Fname2,	'2', UVAR_REQUIRED, "Path and basefilename for second Fstats file");
  XLALregBOOLUserStruct (   help,	'h', UVAR_HELP,     "Print this help/usage message");

  XLALregREALUserStruct (   tol_L1,   	0, UVAR_OPTIONAL, "tolerance on relative error between vectors using L1 norm, between [0,2]");
  XLALregREALUserStruct (   tol_L2,   	0, UVAR_OPTIONAL, "tolerance on relative error between vectors using L2 norm, between [0,2]");
  XLALregREALUserStruct (   tol_angle, 	0, UVAR_OPTIONAL, "tolerance on angle between the two vectors in radians, between [0,pi]");
  XLALregREALUserStruct (   tol_atMax, 	0, UVAR_OPTIONAL, "tolerance on single-sample relative error *at* respective maximum, between [0,2]");

  return XLAL_SUCCESS;

} /* XLALinitUserVars() */

/**
 * comparison specific to pure Fstat-output files (5 entries )
 */
int
XLALcompareFstatFiles ( const LALParsedDataFile *f1, const LALParsedDataFile *f2, VectorComparison tol )
{
  XLAL_CHECK ( (f1 != NULL) && ( f2 != NULL ), XLAL_EINVAL );

  REAL4 eps4 = 100.0 * LAL_REAL4_EPS;
  FstatLine_t XLAL_INIT_DECL(parsed1);
  FstatLine_t XLAL_INIT_DECL(parsed2);

  XLAL_CHECK ( f1->lines->nTokens == f2->lines->nTokens, XLAL_ETOL, "Different number of lines: %d != %d\n", f1->lines->nTokens, f2->lines->nTokens );
  UINT4 nlines = f1->lines->nTokens;

  REAL4Vector *twoF1, *twoF2;
  XLAL_CHECK ( (twoF1 = XLALCreateREAL4Vector ( nlines )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( (twoF2 = XLALCreateREAL4Vector ( nlines )) != NULL, XLAL_EFUNC );

  for (UINT4 i=0; i < nlines ; i++)
    {
      const char *line1 = f1->lines->tokens[i];
      const char *line2 = f2->lines->tokens[i];

      /* read pure Fstats files */
      XLAL_CHECK ( XLALParseFstatLine ( &parsed1, line1 ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK ( XLALParseFstatLine ( &parsed2, line2 ) == XLAL_SUCCESS, XLAL_EFUNC );

      /* compare all template parameters */
      REAL8 relErr;
      if ( (relErr = relError( parsed1.Freq, parsed2.Freq)) > eps4 ) {
	  XLAL_ERROR (XLAL_ETOL, "Relative frequency-error %g ecceeds %g in line %d\n", relErr, eps4, i+1);
      }
      if ( (relErr = relError( parsed1.Alpha, parsed2.Alpha)) > eps4 ) {
        XLAL_ERROR (XLAL_ETOL, "Relative error %g in alpha ecceeds %g in line %d\n", relErr, eps4, i+1);
      }
      if ( (relErr = relError( parsed1.Delta, parsed2.Delta)) > eps4 ) {
        XLAL_ERROR (XLAL_ETOL, "Relative error %g in delta ecceeds %g in line %d\n", relErr, eps4, i+1);
      }
      if ( (relErr = relError( parsed1.f1dot, parsed2.f1dot)) > eps4 ) {
        XLAL_ERROR (XLAL_ETOL, "Relative error %g in f1dot ecceeds %g in line %d\n", relErr, eps4, i+1);
      }
      if ( (relErr = relError( parsed1.f2dot, parsed2.f2dot)) > eps4 ) {
        XLAL_ERROR (XLAL_ETOL, "Relative error %g in f2dot ecceeds %g in line %d\n", relErr, eps4, i+1);
      }
      if ( (relErr = relError( parsed1.f3dot, parsed2.f3dot)) > eps4 ) {
        XLAL_ERROR (XLAL_ETOL, "Relative error %g in f3dot ecceeds %g in line %d\n", relErr, eps4, i+1);
      }

      // and store respective 2F values in vectors for comparison
      twoF1->data[i] = (REAL4)parsed1.TwoF;
      twoF2->data[i] = (REAL4)parsed2.TwoF;

    } /* for i < nlines */

  // ----- finally vector-compare 2F values against given tolerances ----------
  VectorComparison XLAL_INIT_DECL(cmp);
  XLAL_CHECK ( XLALCompareREAL4Vectors ( &cmp, twoF1, twoF2, &tol ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLALDestroyREAL4Vector ( twoF1 );
  XLALDestroyREAL4Vector ( twoF2 );

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
