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
			     const REAL4 Fstar0,
			     const REAL4 *oLGX,
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
  /* maximum allowed difference between recalculated and XLAL result */
  REAL4 tolerance = 1e-06;

  /* compute and compare the results for one set of Fstar0, oLGX values */
  REAL4 Fstar0 = -LAL_REAL4_MAX; /* prior from LR-stat derivation, -Inf means pure line veto, +Inf means pure multi-Fstat */
  REAL4 *oLGX = NULL; /* per-IFO prior odds ratio for line vs. Gaussian noise, NULL is interpreted as oLG[X]=1 for all X */
  printf ("Computing LR-stat for TwoF_multi=%f, TwoFX=(%f,%f), priors F*0=%f and oLGX=NULL...\n", TwoF, TwoFX->data[0], TwoFX->data[1], Fstar0 );
  XLAL_CHECK ( XLALCompareLRSComputations( TwoF, numDetectors, TwoFX, Fstar0, oLGX, tolerance ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* change the priors to catch more possible problems */
  Fstar0 = 10.0;
  REAL4 oLGXarray[numDetectors];
  oLGXarray[0] = 0.5;
  oLGXarray[1] = 0.8;
  oLGX = oLGXarray;

  printf ("Computing LR-stat for TwoF_multi=%f, TwoFX=(%f,%f), priors F*0=%f and oLGX=(%f,%f)...\n", TwoF, TwoFX->data[0], TwoFX->data[1], Fstar0, oLGX[0], oLGX[1] );
  XLAL_CHECK ( XLALCompareLRSComputations( TwoF, numDetectors, TwoFX, Fstar0, oLGX, tolerance ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* free memory */
  XLALDestroyREAL4Vector(TwoFX);

  LALCheckMemoryLeaks();

  return XLAL_SUCCESS;

} /* main */


/**
 * Test function to compute LR-stat values from:
 * XLALComputeLineRobustStat, from scratch
 * and from deprecated XLALComputeLineVeto and XLALComputeLineVetoArray
 * compare the results and exit if tolerance is violated.
 */
int
XLALCompareLRSComputations ( const REAL4 TwoF,			/**< multi-detector  Fstat */
			     const UINT4 numDetectors,		/**< number of detectors */
			     const REAL4Vector *TwoFX,		/**< vector of single-detector Fstats */
			     const REAL4 Fstar0,		/**< amplitude prior normalization for lines */
			     const REAL4 *oLGX,			/**< array of single-detector prior line odds ratio, can be NULL */
			     const REAL4 tolerance		/**< tolerance for comparisons */
                          )
{

  /* conversions between old (rho) and new (F*) notation, REAL4 and REAL8 */
  REAL8 LVrho = exp( 0.25 * ( Fstar0 + log(70.0) ) );
  REAL4 oLG = 0.0;
  REAL8Vector *oLGXREAL8 = NULL;
  XLAL_CHECK ( (oLGXREAL8 = XLALCreateREAL8Vector ( numDetectors )) != NULL, XLAL_EFUNC );

  if ( oLGX ) {
    for ( UINT4 X = 0; X < numDetectors; X++ ) {
      oLGXREAL8->data[X] = (REAL8)oLGX[X];
      oLG += oLGX[X];
    }
  }
  else { /* if oLGX == NULL, assume oLGX=1/numDetectors for all X  ==> oLG = sumX oLGX = 1*/
    oLG = 1.0;
    for (UINT4 X = 0; X < numDetectors; X++) {
      oLGXREAL8->data[X] = 1.0/numDetectors; /* need to set this manually, as old functions still assume oLGX=1 instead */
    }
  } // if ( oLGX == NULL )

  /* further parameter pre-conversions for XLALComputeLineVetoArray() */
  REAL8 logRhoTerm = Fstar0;
  REAL8 logoLGX[numDetectors];
  for (UINT4 X = 0; X < numDetectors; X++) {
    logoLGX[X] = log(oLGXREAL8->data[X]);
  }

  /* compute LR-stat "the pedestrian way", from Eq. (40) of Keitel, Prix, Papa, Leaci, Siddiqi, PR D 89, 064023 (2014),
   * explicit formula for numDet=2:
   * LRS = F - log ( e^F* + e^{F1}*oLG1/oLG + e^{F2}*oLG2/oLG )
   */
  REAL4 LRS_extcomp_terms[3];
  LRS_extcomp_terms[0] = exp(Fstar0)/oLG;
  REAL4 LRS_extcomp_maxterm = LRS_extcomp_terms[0];
  for (UINT4 X = 0; X < numDetectors; X++) {
    LRS_extcomp_terms[1+X] = exp(0.5*TwoFX->data[X]);
    if ( oLGX ) {
      LRS_extcomp_terms[1+X] *= oLGX[X]/oLG;
    }
    else {  /* oLGX=NULL is interpreted as oLGX[X]=1/numDetectors=0.5 for all X ==> oLG=1 */
      LRS_extcomp_terms[1+X] *= 0.5/oLG;
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

  /* these are not the log-Bayes-factor, as computed by XLALComputeLineRobustStat(), so need to correct by log(1+1/oLG) */
  LRS_extcomp_allterms    += log(1+1/oLG);
  LRS_extcomp_notallterms += log(1+1/oLG);

  /* faster version: use only the leading term of the LRS denominator sum */
  LRstatSetup *setup;
  XLAL_CHECK ( (setup = XLALCreateLRstatSetup ( numDetectors, Fstar0, oLGX )) != NULL, XLAL_EFUNC );
  REAL4 LRS_XLAL_notallterms = XLALComputeLRstat ( TwoF, TwoFX->data, setup, FALSE ) / LAL_LOG10E;
  XLAL_CHECK ( xlalErrno == 0, XLAL_EFUNC, "XLALComputeLRstat() failed with xlalErrno = %d\n", xlalErrno );

  /* more precise version: use all terms of the LRS denominator sum */
  REAL4 LRS_XLAL_allterms = XLALComputeLRstat ( TwoF, TwoFX->data, setup, TRUE ) / LAL_LOG10E;
  XLAL_CHECK ( xlalErrno == 0, XLAL_EFUNC, "XLALComputeLRstat() failed with xlalErrno = %d\n", xlalErrno );

  XLALFree ( setup ); setup = NULL;

  /* compute relative deviations */
  REAL4 diff_allterms              = fabs( LRS_XLAL_allterms             - LRS_extcomp_allterms    ) / ( 0.5 * ( LRS_XLAL_allterms             + LRS_extcomp_allterms    ));
  REAL4 diff_notallterms           = fabs( LRS_XLAL_notallterms          - LRS_extcomp_notallterms ) / ( 0.5 * ( LRS_XLAL_notallterms          + LRS_extcomp_notallterms ));

  /* output results and deviations and return with error when tolerances are violated */
  printf ( "Externally recomputed       with  allterms: LRS=%f\n", LRS_extcomp_allterms );
  printf ( "Externally recomputed       with !allterms: LRS=%f\n", LRS_extcomp_notallterms );
  printf ( "XLALComputeLineRobustStat() with  allterms: LRS=%f (rel. dev.: %f)", LRS_XLAL_allterms,             diff_allterms );
  XLAL_CHECK ( XLALCheckLRSDifferences ( diff_allterms,             tolerance, "XLALComputeLineRobustStat() with useAllTerms=TRUE" )  == XLAL_SUCCESS, XLAL_EFUNC );
  printf ( "XLALComputeLineRobustStat() with !allterms: LRS=%f (rel. dev.: %f)", LRS_XLAL_notallterms,          diff_notallterms );
  XLAL_CHECK ( XLALCheckLRSDifferences ( diff_notallterms,          tolerance, "XLALComputeLineRobustStat() with useAllTerms=TRUE" )  == XLAL_SUCCESS, XLAL_EFUNC );

  /* also test against deprecated rho-notation functions for consistency
   * need to correct for different prior parametrization,
   * LRS_new=LRS_old+numDetectors*log(1+1/oLG)
   */
  REAL4 old_LV_corr = log(1+oLG);

  xlalErrno = 0;
  REAL4 LRS_XLAL_rho_notallterms      = XLALComputeLineVeto       ( TwoF, TwoFX, LVrho, oLGXREAL8, FALSE ) + old_LV_corr;
  XLAL_CHECK ( xlalErrno == 0, XLAL_EFUNC, "XLALComputeLineVeto() failed with xlalErrno = %d\n", xlalErrno );

  xlalErrno = 0;
  REAL4 LRS_XLAL_rhoarray_notallterms = XLALComputeLineVetoArray  ( TwoF, numDetectors, TwoFX->data, logRhoTerm, logoLGX, FALSE ) + old_LV_corr;
  XLAL_CHECK ( xlalErrno == 0, XLAL_EFUNC, "XLALComputeLineVetoArray() failed with xlalErrno = %d\n", xlalErrno );

  xlalErrno = 0;
  REAL4 LRS_XLAL_rho_allterms         = XLALComputeLineVeto       ( TwoF, TwoFX, LVrho, oLGXREAL8, TRUE ) + old_LV_corr;
  XLAL_CHECK ( xlalErrno == 0, XLAL_EFUNC, "XLALComputeLineVeto() failed with xlalErrno = %d\n", xlalErrno );

  xlalErrno = 0;
  REAL4 LRS_XLAL_rhoarray_allterms    = XLALComputeLineVetoArray  ( TwoF, numDetectors, TwoFX->data, logRhoTerm, logoLGX, TRUE ) + old_LV_corr;
  XLAL_CHECK ( xlalErrno == 0, XLAL_EFUNC, "XLALComputeLineVetoArray() failed with xlalErrno = %d\n", xlalErrno );

  REAL4 diff_rho_allterms          = fabs( LRS_XLAL_rho_allterms         - LRS_extcomp_allterms    ) / ( 0.5 * ( LRS_XLAL_rho_allterms         + LRS_extcomp_allterms ));
  REAL4 diff_rhoarray_allterms     = fabs( LRS_XLAL_rhoarray_allterms    - LRS_extcomp_allterms    ) / ( 0.5 * ( LRS_XLAL_rhoarray_allterms    + LRS_extcomp_allterms ));
  REAL4 diff_rho_notallterms       = fabs( LRS_XLAL_rho_notallterms      - LRS_extcomp_notallterms ) / ( 0.5 * ( LRS_XLAL_rho_notallterms      + LRS_extcomp_notallterms ));
  REAL4 diff_rhoarray_notallterms  = fabs( LRS_XLAL_rhoarray_notallterms - LRS_extcomp_notallterms ) / ( 0.5 * ( LRS_XLAL_rhoarray_notallterms + LRS_extcomp_notallterms ));

  printf( "Legacy functions with rho-notation, corrected as LRS_new=LRS_old+log(1+oLG)=LRS_old+%f:\n", old_LV_corr );
  printf ( "XLALComputeLineVeto()       with  allterms: LRS=%f (rel. dev.: %f)", LRS_XLAL_rho_allterms,         diff_rho_allterms );
  XLAL_CHECK ( XLALCheckLRSDifferences ( diff_rho_allterms,         tolerance, "XLALComputeLineVeto()       with useAllTerms=TRUE" )  == XLAL_SUCCESS, XLAL_EFUNC );
  printf ( "XLALComputeLineVetoArray()  with  allterms: LRS=%f (rel. dev.: %f)", LRS_XLAL_rhoarray_allterms,    diff_rhoarray_allterms );
  XLAL_CHECK ( XLALCheckLRSDifferences ( diff_rhoarray_allterms,    tolerance, "XLALComputeLineVetoArray()  with useAllTerms=TRUE" )  == XLAL_SUCCESS, XLAL_EFUNC );
  printf ( "XLALComputeLineVeto()       with !allterms: LRS=%f (rel. dev.: %f)", LRS_XLAL_rho_notallterms,      diff_rho_notallterms );
  XLAL_CHECK ( XLALCheckLRSDifferences ( diff_rho_notallterms,      tolerance, "XLALComputeLineVeto()       with useAllTerms=FALSE" ) == XLAL_SUCCESS, XLAL_EFUNC );
  printf ( "XLALComputeLineVetoArray()  with !allterms: LRS=%f (rel. dev.: %f)", LRS_XLAL_rhoarray_notallterms, diff_rhoarray_notallterms );
  XLAL_CHECK ( XLALCheckLRSDifferences ( diff_rhoarray_notallterms, tolerance, "XLALComputeLineVetoArray()  with useAllTerms=FALSE" ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLALDestroyREAL8Vector(oLGXREAL8);

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
