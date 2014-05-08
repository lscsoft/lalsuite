/*
 * Copyright (C) 2011, 2014 David Keitel
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

/* ---------- Includes -------------------- */
#include <lal/AVFactories.h>
#include <lal/LineRobustStats.h>
#include <lal/LogPrintf.h>
#include <lal/UserInput.h>

/******************************************************
 *  Error codes and messages.
 */

/* ---------- Defines -------------------- */
#define TRUE (1==1)
#define FALSE (1==0)

/*---------- internal prototypes ----------*/
int
XLALCompareLRSComputations ( const REAL4 TwoF,
			     const UINT4 numDetectors,
			     const REAL4Vector *TwoFX,
			     const REAL8 LVrho,
			     const REAL8Vector *oLGX,
			     const REAL4 tolerance
			   );

int
XLALCheckLRSDifferences ( const REAL4 diff,
			  const REAL4 tolerance,
			  const CHAR *casestring
			);

/* ###################################  MAIN  ################################### */

int main( int argc, char *argv[]) {

  /* sanity check for input arguments */
  XLAL_CHECK ( argc == 1, XLAL_EINVAL, "The executable '%s' doesn't support any input arguments right now.\n", argv[0] );

  printf ("Starting test...\n");

  /* set up single- and multi-IFO F-stat input */
  REAL4 TwoF = 7.0;
  UINT4 numDetectors = 2;
  REAL4Vector *TwoFX = NULL;
  XLAL_CHECK ( (TwoFX = XLALCreateREAL4Vector ( numDetectors )) != NULL, XLAL_EFUNC );
  TwoFX->data[0] = 4.0;
  TwoFX->data[1] = 12.0;
  REAL8 LVrho = 0.0; /* prior from LR-stat derivation, 0 means pure line veto, +inf means pure multi-Fstat */
  REAL8Vector *oLGX = NULL; /* per-IFO prior odds ratio for line vs. Gaussian noise, NULL is interpreted as oLG[X]=1 for all X */

  /* maximum allowed difference between recalculated and XLAL result */
  REAL4 tolerance = 1e-03;

  /* compute and compare the results for one set of rho, oLGX values */
  printf ("Computing LR-stat for TwoF_multi=%f, TwoFX=(%f,%f), priors rho=%f and oLGX=NULL...\n", TwoF, TwoFX->data[0], TwoFX->data[1], LVrho );
  XLAL_CHECK ( XLALCompareLRSComputations( TwoF, numDetectors, TwoFX, LVrho, NULL, tolerance ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* change the priors to catch more possible problems */
  LVrho = 5.0;
  XLAL_CHECK ( (oLGX = XLALCreateREAL8Vector ( numDetectors )) != NULL, XLAL_EFUNC );
  oLGX->data[0] = 0.5;
  oLGX->data[1] = 0.8;

  /* compute and compare the results for second set of LVrho, oLGX values */
  printf ("Computing LV-stat for TwoF_multi=%f, TwoFX=(%f,%f), priors rho=%f and oLGX=(%f,%f)...\n", TwoF, TwoFX->data[0], TwoFX->data[1], LVrho, oLGX->data[0], oLGX->data[1] );
  XLAL_CHECK ( XLALCompareLRSComputations( TwoF, numDetectors, TwoFX, LVrho, oLGX, tolerance ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* free memory */
  XLALDestroyREAL4Vector(TwoFX);
  XLALDestroyREAL8Vector(oLGX);

  LALCheckMemoryLeaks();

  return XLAL_SUCCESS;

} /* main */


/**
 * Test function to compute LV-stat values both from scratch and by XLALComputeLineVeto,
 * compare the results and exit if tolerance is violated.
 */
int
XLALCompareLRSComputations ( const REAL4 TwoF,			/**< multi-detector  Fstat */
			     const UINT4 numDetectors,		/**< number of detectors */
			     const REAL4Vector *TwoFX,		/**< vector of single-detector Fstats */
			     const REAL8 LVrho,			/**< amplitude prior normalization for lines */
			     const REAL8Vector *oLGX,		/**< vector of single-detector prior line odds ratio, default to oLGX=1 for all X if NULL */
			     const REAL4 tolerance		/**< tolerance for comparisons */
                          )
{

  /* compute LR-stat "the pedestrian way", explicit formula for numDet=2:
   * LRS = F - log( LVrho^4/70 + oLG(1)*exp(F1) + oLG(2)*exp(F2) )
   */
  REAL4 LRS_extcomp_terms[3];
  LRS_extcomp_terms[0] = pow(LVrho,4)/70.0;
  REAL4 LRS_extcomp_maxterm = LRS_extcomp_terms[0];
  for (UINT4 X = 0; X < numDetectors; X++) {
    LRS_extcomp_terms[1+X] = exp(0.5*TwoFX->data[X]);
    if ( oLGX ) {
      LRS_extcomp_terms[1+X] *= oLGX->data[X];
    }
    if ( LRS_extcomp_terms[1+X] > LRS_extcomp_maxterm ) {
      LRS_extcomp_maxterm = LRS_extcomp_terms[1+X];
    }
  }
  REAL4 LRS_extcomp_notallterms = 0.5*TwoF - log(LRS_extcomp_maxterm);
  REAL4 LRS_extcomp_denom = 0.0;
  for (UINT4 X = 0; X < 1+numDetectors; X++) {
    LRS_extcomp_denom += LRS_extcomp_terms[X];
  }
  REAL4 LRS_extcomp_allterms = 0.5*TwoF - log( LRS_extcomp_denom );

  /* parameter pre-conversion for XLALComputeLineVetoArray() */
  REAL8 logRhoTerm = 4.0 * log(LVrho) - log(70.0);
  REAL8 *logoLGX = NULL;
  REAL8 logoLGXtemp[numDetectors];
  if ( oLGX ) {
    for (UINT4 X = 0; X < numDetectors; X++) {
      logoLGXtemp[X] = log(oLGX->data[X]);
    }
    logoLGX = logoLGXtemp;
  }

  /* faster version: use only the leading term of the LV denominator sum */
  xlalErrno = 0;
  REAL4 LRS_XLAL_notallterms = XLALComputeLineVeto ( TwoF, TwoFX, LVrho, oLGX, FALSE );
  XLAL_CHECK ( xlalErrno == 0, XLAL_EFUNC, "XLALComputeLineVeto() failed with xlalErrno = %d\n", xlalErrno );
  xlalErrno = 0;
  REAL4 LRS_XLAL_array_notallterms = XLALComputeLineVetoArray ( TwoF, numDetectors, TwoFX->data, logRhoTerm, logoLGX, FALSE );
  XLAL_CHECK ( xlalErrno == 0, XLAL_EFUNC, "XLALComputeLineVetoArray() failed with xlalErrno = %d\n", xlalErrno );

  /* more precise version: use all terms of the LRS denominator sum */
  xlalErrno = 0;
  REAL4 LRS_XLAL_allterms = XLALComputeLineVeto ( TwoF, TwoFX, LVrho, oLGX, TRUE );
  XLAL_CHECK ( xlalErrno == 0, XLAL_EFUNC, "XLALComputeLineVeto() failed with xlalErrno = %d\n", xlalErrno );
  xlalErrno = 0;
  REAL4 LRS_XLAL_array_allterms = XLALComputeLineVetoArray ( TwoF, numDetectors, TwoFX->data, logRhoTerm, logoLGX, TRUE );
  XLAL_CHECK ( xlalErrno == 0, XLAL_EFUNC, "XLALComputeLineVetoArray() failed with xlalErrno = %d\n", xlalErrno );

  /* compute relative deviations */
  REAL4 diff_allterms          = fabs( LRS_XLAL_allterms          - LRS_extcomp_allterms )    / ( 0.5 * ( LRS_XLAL_allterms          + LRS_extcomp_allterms ));
  REAL4 diff_array_allterms    = fabs( LRS_XLAL_array_allterms    - LRS_extcomp_allterms )    / ( 0.5 * ( LRS_XLAL_array_allterms    + LRS_extcomp_allterms ));
  REAL4 diff_notallterms       = fabs( LRS_XLAL_notallterms       - LRS_extcomp_notallterms ) / ( 0.5 * ( LRS_XLAL_notallterms       + LRS_extcomp_notallterms ));
  REAL4 diff_array_notallterms = fabs( LRS_XLAL_array_notallterms - LRS_extcomp_notallterms ) / ( 0.5 * ( LRS_XLAL_array_notallterms + LRS_extcomp_notallterms ));

  /* output results and deviations and return with error when tolerances are violated */
  printf ("Externally recomputed      with  allterms: LV=%f\n", LRS_extcomp_allterms);
  printf ("Externally recomputed      with !allterms: LV=%f\n", LRS_extcomp_notallterms);
  printf ("XLALComputeLineVeto()      with  allterms: LV=%f (rel. dev.: %f)", LRS_XLAL_allterms,          diff_allterms);
  XLAL_CHECK ( XLALCheckLRSDifferences ( diff_allterms,          tolerance, "XLALComputeLineVeto() with useAllTerms=TRUE" )       == XLAL_SUCCESS, XLAL_EFUNC );
  printf ("XLALComputeLineVetoArray() with  allterms: LV=%f (rel. dev.: %f)", LRS_XLAL_array_allterms,    diff_array_allterms);
  XLAL_CHECK ( XLALCheckLRSDifferences ( diff_array_allterms,    tolerance, "XLALComputeLineVetoArray() with useAllTerms=TRUE" )  == XLAL_SUCCESS, XLAL_EFUNC );
  printf ("XLALComputeLineVeto()      with !allterms: LV=%f (rel. dev.: %f)", LRS_XLAL_notallterms,       diff_notallterms );
  XLAL_CHECK ( XLALCheckLRSDifferences ( diff_notallterms,       tolerance, "XLALComputeLineVeto() with useAllTerms=FALSE" )      == XLAL_SUCCESS, XLAL_EFUNC );
  printf ("XLALComputeLineVetoArray() with !allterms: LV=%f (rel. dev.: %f)", LRS_XLAL_array_notallterms, diff_array_notallterms );
  XLAL_CHECK ( XLALCheckLRSDifferences ( diff_array_notallterms, tolerance, "XLALComputeLineVetoArray() with useAllTerms=FALSE" ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;

} /* XLALCompareLRSComputations() */


int
XLALCheckLRSDifferences ( const REAL4 diff,
			  const REAL4 tolerance,
			  const CHAR *casestring
			)
{

  if ( fabs(diff) <= tolerance ) {
    printf (" ==> OK!\n");
  }
  else {
    printf (" ==> BAD!\n");
    XLAL_ERROR ( XLAL_EFAILED, "\nTolerance %f exceeded for %s!\n", tolerance, casestring );
    return XLAL_FAILURE;
  }

  return XLAL_SUCCESS;

} /* XLALCheckLRSDifferences() */
