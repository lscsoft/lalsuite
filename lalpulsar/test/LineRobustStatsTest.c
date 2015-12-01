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
XLALCompareBSGLComputations ( const REAL4 TwoF,
			     const UINT4 numDetectors,
			     const REAL4Vector *TwoFX,
			     const REAL4 Fstar0,
			     const REAL4 *oLGX,
			     const REAL4 tolerance
			   );

REAL4
XLALComputeBSGLPedestrian ( const REAL4 TwoF,
			    const UINT4 numDetectors,
			    const REAL4Vector *TwoFX,
			    const REAL4 Fstar0,
			    const REAL4 *oLGX,
			    const BOOLEAN useLogCorrection
			   );

int
XLALCheckBSGLDifferences ( const REAL4 diff,
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
  REAL4 Fstar0 = -LAL_REAL4_MAX; /* prior from BSGL derivation, -Inf means pure line veto, +Inf means pure multi-Fstat */
  REAL4 *oLGX = NULL; /* per-IFO prior odds ratio for line vs. Gaussian noise, NULL is interpreted as oLG[X]=1 for all X */
  printf ("Computing BSGL for TwoF_multi=%f, TwoFX=(%f,%f), priors F*0=%f and oLGX=NULL...\n", TwoF, TwoFX->data[0], TwoFX->data[1], Fstar0 );
  XLAL_CHECK ( XLALCompareBSGLComputations( TwoF, numDetectors, TwoFX, Fstar0, oLGX, tolerance ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* change the priors to catch more possible problems */
  Fstar0 = 10.0;
  REAL4 oLGXarray[numDetectors];
  oLGXarray[0] = 0.5;
  oLGXarray[1] = 0.8;
  oLGX = oLGXarray;

  printf ("Computing BSGL for TwoF_multi=%f, TwoFX=(%f,%f), priors F*0=%f and oLGX=(%f,%f)...\n", TwoF, TwoFX->data[0], TwoFX->data[1], Fstar0, oLGX[0], oLGX[1] );
  XLAL_CHECK ( XLALCompareBSGLComputations( TwoF, numDetectors, TwoFX, Fstar0, oLGX, tolerance ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* free memory */
  XLALDestroyREAL4Vector(TwoFX);

  LALCheckMemoryLeaks();

  return XLAL_SUCCESS;

} /* main */


/**
 * Test function to compute BSGL values from XLALComputeBSGL
 * against those recomputed from scratch ("pedestrian" formula);
 * compare the results and exit if tolerance is violated.
 */
int
XLALCompareBSGLComputations ( const REAL4 TwoF,			/**< multi-detector  Fstat */
			     const UINT4 numDetectors,		/**< number of detectors */
			     const REAL4Vector *TwoFX,		/**< vector of single-detector Fstats */
			     const REAL4 Fstar0,		/**< amplitude prior normalization for lines */
			     const REAL4 *oLGX,			/**< array of single-detector prior line odds ratio, can be NULL */
			     const REAL4 tolerance		/**< tolerance for comparisons */
                          )
{

  /* pedestrian version, with and without log corrections
  REAL4 log10BSGL_extcomp_notallterms = XLALComputeBSGLPedestrian ( TwoF, numDetectors, TwoFX, Fstar0, oLGX, FALSE );
  REAL4 log10BSGL_extcomp_allterms    = XLALComputeBSGLPedestrian ( TwoF, numDetectors, TwoFX, Fstar0, oLGX, TRUE );

  /* faster version: use only the leading term of the BSGL denominator sum */
  UINT4 numSegs = 1;
  BSGLSetup *setup_noLogCorrection;
  XLAL_CHECK ( (setup_noLogCorrection = XLALCreateBSGLSetup ( numDetectors, Fstar0, oLGX, FALSE, numSegs )) != NULL, XLAL_EFUNC );
  REAL4 log10BSGL_XLAL_notallterms = XLALComputeBSGL ( TwoF, TwoFX->data, setup_noLogCorrection );
  XLAL_CHECK ( xlalErrno == 0, XLAL_EFUNC, "XLALComputeBSGL() failed with xlalErrno = %d\n", xlalErrno );
  XLALFree ( setup_noLogCorrection ); setup_noLogCorrection = NULL;

  /* more precise version: use all terms of the BSGL denominator sum */
  BSGLSetup *setup_withLogCorrection;
  XLAL_CHECK ( (setup_withLogCorrection = XLALCreateBSGLSetup ( numDetectors, Fstar0, oLGX, TRUE, numSegs )) != NULL, XLAL_EFUNC );
  REAL4 log10BSGL_XLAL_allterms = XLALComputeBSGL ( TwoF, TwoFX->data, setup_withLogCorrection );
  XLAL_CHECK ( xlalErrno == 0, XLAL_EFUNC, "XLALComputeBSGL() failed with xlalErrno = %d\n", xlalErrno );
  XLALFree ( setup_withLogCorrection ); setup_withLogCorrection = NULL;

  /* compute relative deviations */
  REAL4 diff_allterms    = fabs( log10BSGL_XLAL_allterms    - log10BSGL_extcomp_allterms    ) / ( 0.5 * ( log10BSGL_XLAL_allterms    + log10BSGL_extcomp_allterms    ));
  REAL4 diff_notallterms = fabs( log10BSGL_XLAL_notallterms - log10BSGL_extcomp_notallterms ) / ( 0.5 * ( log10BSGL_XLAL_notallterms + log10BSGL_extcomp_notallterms ));

  /* output results and deviations and return with error when tolerances are violated */
  printf ( "Externally recomputed     with  allterms: log10BSGL=%f\n",               log10BSGL_extcomp_allterms );
  printf ( "Externally recomputed     with !allterms: log10BSGL=%f\n",               log10BSGL_extcomp_notallterms );
  printf ( "XLALComputeBSGL()         with  allterms: log10BSGL=%f (rel. dev.: %f)", log10BSGL_XLAL_allterms,    diff_allterms );
  XLAL_CHECK ( XLALCheckBSGLDifferences ( diff_allterms,    tolerance, "XLALComputeBSGL() with useAllTerms=TRUE" ) == XLAL_SUCCESS, XLAL_EFUNC );
  printf ( "XLALComputeBSGL()         with !allterms: log10BSGL=%f (rel. dev.: %f)", log10BSGL_XLAL_notallterms, diff_notallterms );
  XLAL_CHECK ( XLALCheckBSGLDifferences ( diff_notallterms, tolerance, "XLALComputeBSGL() with useAllTerms=TRUE" ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;

} /* XLALCompareBSGLComputations() */

/**
 * compute BSGL "the pedestrian way", from Eq. (40) of Keitel, Prix, Papa, Leaci, Siddiqi, PR D 89, 064023 (2014),
 * explicit formula for numDet=2:
 * log10 BSGL = F - log ( e^F* + e^{F1}*oLG1/oLG + e^{F2}*oLG2/oLG )
 */
REAL4
XLALComputeBSGLPedestrian ( const REAL4 TwoF,			/**< multi-detector  Fstat */
			    const UINT4 numDetectors,		/**< number of detectors */
			    const REAL4Vector *TwoFX,		/**< vector of single-detector Fstats */
			    const REAL4 Fstar0,			/**< amplitude prior normalization for lines */
			    const REAL4 *oLGX,			/**< array of single-detector prior line odds ratio, can be NULL */
			    const BOOLEAN useLogCorrection	//!< include log-term correction or not
                          )
{

  REAL4 log10BSGL_extcomp = 0;

  /* sum up overall prior odds */
  REAL4 oLG = 0.0;
  if ( oLGX ) {
    for ( UINT4 X = 0; X < numDetectors; X++ ) {
      oLG += oLGX[X];
    }
  }
  else { /* if oLGX == NULL, assume oLGX=1/numDetectors for all X  ==> oLG = sumX oLGX = 1*/
    oLG = 1.0;
  } // if ( oLGX == NULL )

  REAL4 denomterms[3];
  denomterms[0] = exp(Fstar0)/oLG;
  REAL4 maxterm = denomterms[0];
  for (UINT4 X = 0; X < numDetectors; X++) {
    denomterms[1+X] = exp(0.5*TwoFX->data[X]);
    if ( oLGX ) {
      denomterms[1+X] *= oLGX[X]/oLG;
    }
    else {  /* oLGX=NULL is interpreted as oLGX[X]=1/numDetectors=0.5 for all X ==> oLG=1 */
      denomterms[1+X] *= 0.5/oLG;
    }
    if ( denomterms[1+X] > maxterm ) {
      maxterm = denomterms[1+X];
    }
  }

  if ( !useLogCorrection ) {
    log10BSGL_extcomp = 0.5*TwoF - log(maxterm);
  }
  else {
    REAL4 BSGL_extcomp_denom = 0.0;
    for (UINT4 X = 0; X < 1+numDetectors; X++) {
      BSGL_extcomp_denom += denomterms[X];
    }
    log10BSGL_extcomp = 0.5*TwoF - log( BSGL_extcomp_denom );
  }

  /* this is not yet the log-Bayes-factor, as computed by XLALComputeBSGL(), so need to correct by log(1+1/oLG) */
  log10BSGL_extcomp += log(1+1/oLG);
  /* and actually switch to log10 */
  log10BSGL_extcomp *= LAL_LOG10E;

  return log10BSGL_extcomp;

} /* XLALComputeBSGLPedestrian() */

int
XLALCheckBSGLDifferences ( const REAL4 diff,
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

} /* XLALCheckBSGLDifferences() */
