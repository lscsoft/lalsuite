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

int
XLALCheckBSGLDifferences ( const REAL4 diff,
			  const REAL4 tolerance,
			  const CHAR *casestring
			);

// ----- deprecated API, kept here locally for comparisob -----

REAL4
XLALComputeLineVeto ( const REAL4 TwoF,
		      const REAL4Vector *TwoFX,
		      const REAL8 rhomaxline,
		      const REAL8Vector *lX,
		      const BOOLEAN useAllTerms
);

REAL4
XLALComputeLineVetoArray ( const REAL4 TwoF,
			   const UINT4 numDetectors,
			   const REAL4 *TwoFX,
			   const REAL8 logRhoTerm,
			   const REAL8 *loglX,
			   const BOOLEAN useAllTerms
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
 * Test function to compute BSGL values from:
 * XLALComputeBSGL, from scratch
 * and from deprecated XLALComputeLineVeto and XLALComputeLineVetoArray
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

  /* compute BSGL "the pedestrian way", from Eq. (40) of Keitel, Prix, Papa, Leaci, Siddiqi, PR D 89, 064023 (2014),
   * explicit formula for numDet=2:
   * log10 BSGL = F - log ( e^F* + e^{F1}*oLG1/oLG + e^{F2}*oLG2/oLG )
   */
  REAL4 BSGL_extcomp_terms[3];
  BSGL_extcomp_terms[0] = exp(Fstar0)/oLG;
  REAL4 BSGL_extcomp_maxterm = BSGL_extcomp_terms[0];
  for (UINT4 X = 0; X < numDetectors; X++) {
    BSGL_extcomp_terms[1+X] = exp(0.5*TwoFX->data[X]);
    if ( oLGX ) {
      BSGL_extcomp_terms[1+X] *= oLGX[X]/oLG;
    }
    else {  /* oLGX=NULL is interpreted as oLGX[X]=1/numDetectors=0.5 for all X ==> oLG=1 */
      BSGL_extcomp_terms[1+X] *= 0.5/oLG;
    }
    if ( BSGL_extcomp_terms[1+X] > BSGL_extcomp_maxterm ) {
      BSGL_extcomp_maxterm = BSGL_extcomp_terms[1+X];
    }
  }
  REAL4 BSGL_extcomp_notallterms = 0.5*TwoF - log(BSGL_extcomp_maxterm);
  REAL4 BSGL_extcomp_denom = 0.0;
  for (UINT4 X = 0; X < 1+numDetectors; X++) {
    BSGL_extcomp_denom += BSGL_extcomp_terms[X];
  }
  REAL4 BSGL_extcomp_allterms = 0.5*TwoF - log( BSGL_extcomp_denom );

  /* these are not the log-Bayes-factor, as computed by XLALComputeBSGL(), so need to correct by log(1+1/oLG) */
  BSGL_extcomp_allterms    += log(1+1/oLG);
  BSGL_extcomp_notallterms += log(1+1/oLG);

  /* faster version: use only the leading term of the BSGL denominator sum */
  BSGLSetup *setup_noLogCorrection;
  XLAL_CHECK ( (setup_noLogCorrection = XLALCreateBSGLSetup ( numDetectors, Fstar0, oLGX, FALSE )) != NULL, XLAL_EFUNC );
  REAL4 BSGL_XLAL_notallterms = XLALComputeBSGL ( TwoF, TwoFX->data, setup_noLogCorrection ) / LAL_LOG10E;
  XLAL_CHECK ( xlalErrno == 0, XLAL_EFUNC, "XLALComputeBSGL() failed with xlalErrno = %d\n", xlalErrno );
  XLALFree ( setup_noLogCorrection ); setup_noLogCorrection = NULL;

  /* more precise version: use all terms of the BSGL denominator sum */
  BSGLSetup *setup_withLogCorrection;
  XLAL_CHECK ( (setup_withLogCorrection = XLALCreateBSGLSetup ( numDetectors, Fstar0, oLGX, TRUE )) != NULL, XLAL_EFUNC );
  REAL4 BSGL_XLAL_allterms = XLALComputeBSGL ( TwoF, TwoFX->data, setup_withLogCorrection ) / LAL_LOG10E;
  XLAL_CHECK ( xlalErrno == 0, XLAL_EFUNC, "XLALComputeBSGL() failed with xlalErrno = %d\n", xlalErrno );
  XLALFree ( setup_withLogCorrection ); setup_withLogCorrection = NULL;

  /* compute relative deviations */
  REAL4 diff_allterms              = fabs( BSGL_XLAL_allterms             - BSGL_extcomp_allterms    ) / ( 0.5 * ( BSGL_XLAL_allterms             + BSGL_extcomp_allterms    ));
  REAL4 diff_notallterms           = fabs( BSGL_XLAL_notallterms          - BSGL_extcomp_notallterms ) / ( 0.5 * ( BSGL_XLAL_notallterms          + BSGL_extcomp_notallterms ));

  /* output results and deviations and return with error when tolerances are violated */
  printf ( "Externally recomputed       with  allterms: BSGL=%f\n", BSGL_extcomp_allterms );
  printf ( "Externally recomputed       with !allterms: BSGL=%f\n", BSGL_extcomp_notallterms );
  printf ( "XLALComputeBSGL()         with  allterms: BSGL=%f (rel. dev.: %f)", BSGL_XLAL_allterms,             diff_allterms );
  XLAL_CHECK ( XLALCheckBSGLDifferences ( diff_allterms,             tolerance, "XLALComputeBSGL() with useAllTerms=TRUE" )  == XLAL_SUCCESS, XLAL_EFUNC );
  printf ( "XLALComputeBSGL()         with !allterms: BSGL=%f (rel. dev.: %f)", BSGL_XLAL_notallterms,          diff_notallterms );
  XLAL_CHECK ( XLALCheckBSGLDifferences ( diff_notallterms,          tolerance, "XLALComputeBSGL() with useAllTerms=TRUE" )  == XLAL_SUCCESS, XLAL_EFUNC );

  /* also test against deprecated rho-notation functions for consistency
   * need to correct for different prior parametrization,
   * BSGL_new=LV_old+numDetectors*log(1+oLG)
   */
  REAL4 old_LV_corr = log(1+oLG);

  xlalErrno = 0;
  REAL4 BSGL_XLAL_rho_notallterms      = XLALComputeLineVeto       ( TwoF, TwoFX, LVrho, oLGXREAL8, FALSE ) + old_LV_corr;
  XLAL_CHECK ( xlalErrno == 0, XLAL_EFUNC, "XLALComputeLineVeto() failed with xlalErrno = %d\n", xlalErrno );

  xlalErrno = 0;
  REAL4 BSGL_XLAL_rhoarray_notallterms = XLALComputeLineVetoArray  ( TwoF, numDetectors, TwoFX->data, logRhoTerm, logoLGX, FALSE ) + old_LV_corr;
  XLAL_CHECK ( xlalErrno == 0, XLAL_EFUNC, "XLALComputeLineVetoArray() failed with xlalErrno = %d\n", xlalErrno );

  xlalErrno = 0;
  REAL4 BSGL_XLAL_rho_allterms         = XLALComputeLineVeto       ( TwoF, TwoFX, LVrho, oLGXREAL8, TRUE ) + old_LV_corr;
  XLAL_CHECK ( xlalErrno == 0, XLAL_EFUNC, "XLALComputeLineVeto() failed with xlalErrno = %d\n", xlalErrno );

  xlalErrno = 0;
  REAL4 BSGL_XLAL_rhoarray_allterms    = XLALComputeLineVetoArray  ( TwoF, numDetectors, TwoFX->data, logRhoTerm, logoLGX, TRUE ) + old_LV_corr;
  XLAL_CHECK ( xlalErrno == 0, XLAL_EFUNC, "XLALComputeLineVetoArray() failed with xlalErrno = %d\n", xlalErrno );

  REAL4 diff_rho_allterms          = fabs( BSGL_XLAL_rho_allterms         - BSGL_extcomp_allterms    ) / ( 0.5 * ( BSGL_XLAL_rho_allterms         + BSGL_extcomp_allterms ));
  REAL4 diff_rhoarray_allterms     = fabs( BSGL_XLAL_rhoarray_allterms    - BSGL_extcomp_allterms    ) / ( 0.5 * ( BSGL_XLAL_rhoarray_allterms    + BSGL_extcomp_allterms ));
  REAL4 diff_rho_notallterms       = fabs( BSGL_XLAL_rho_notallterms      - BSGL_extcomp_notallterms ) / ( 0.5 * ( BSGL_XLAL_rho_notallterms      + BSGL_extcomp_notallterms ));
  REAL4 diff_rhoarray_notallterms  = fabs( BSGL_XLAL_rhoarray_notallterms - BSGL_extcomp_notallterms ) / ( 0.5 * ( BSGL_XLAL_rhoarray_notallterms + BSGL_extcomp_notallterms ));

  printf( "Legacy functions with rho-notation, corrected as BSGL_new=LV_old+log(1+oLG)=LV_old+%f:\n", old_LV_corr );
  printf ( "XLALComputeLineVeto()       with  allterms: BSGL=%f (rel. dev.: %f)", BSGL_XLAL_rho_allterms,         diff_rho_allterms );
  XLAL_CHECK ( XLALCheckBSGLDifferences ( diff_rho_allterms,         tolerance, "XLALComputeLineVeto()       with useAllTerms=TRUE" )  == XLAL_SUCCESS, XLAL_EFUNC );
  printf ( "XLALComputeLineVetoArray()  with  allterms: BSGL=%f (rel. dev.: %f)", BSGL_XLAL_rhoarray_allterms,    diff_rhoarray_allterms );
  XLAL_CHECK ( XLALCheckBSGLDifferences ( diff_rhoarray_allterms,    tolerance, "XLALComputeLineVetoArray()  with useAllTerms=TRUE" )  == XLAL_SUCCESS, XLAL_EFUNC );
  printf ( "XLALComputeLineVeto()       with !allterms: BSGL=%f (rel. dev.: %f)", BSGL_XLAL_rho_notallterms,      diff_rho_notallterms );
  XLAL_CHECK ( XLALCheckBSGLDifferences ( diff_rho_notallterms,      tolerance, "XLALComputeLineVeto()       with useAllTerms=FALSE" ) == XLAL_SUCCESS, XLAL_EFUNC );
  printf ( "XLALComputeLineVetoArray()  with !allterms: BSGL=%f (rel. dev.: %f)", BSGL_XLAL_rhoarray_notallterms, diff_rhoarray_notallterms );
  XLAL_CHECK ( XLALCheckBSGLDifferences ( diff_rhoarray_notallterms, tolerance, "XLALComputeLineVetoArray()  with useAllTerms=FALSE" ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLALDestroyREAL8Vector(oLGXREAL8);

  return XLAL_SUCCESS;

} /* XLALCompareBSGLComputations() */


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

// ---------- OLD API functions - deprecated, but kept here locally for comparison ------------------------------

/**
 * Deprecated function to compute Line Veto statistics from multi- and single-detector \f$ \mathcal{F} \f$-stats,
 * using outdated \f$ \rho \f$ notation and prior normalization.
 * This is not the log-Bayes-factor!
 * \f$ \mathrm{LV} = \mathcal{F} - \log \left( \frac{\rho_{\mathrm{max,line}}^4}{70} + \sum_X l^X e^{\mathcal{F}^X} \right) \f$
 *
 * \deprecated use XLALComputeBSGL() instead.
 *
 * Also this is just a wrapper for XLALComputeLineVetoArray, which is faster for many LV values at identical priors.
 * This function here just translates REAL4Vectors to fixed REAL4 arrays,
 * and linear priors rhomaxline and lX to logarithmic values.
 */
REAL4 XLALComputeLineVeto ( const REAL4 TwoF,			/**< multi-detector  \f$ \mathcal{F} \f$-stat */
                            const REAL4Vector *TwoFXvec,	/**< vector of single-detector \f$ \mathcal{F} \f$-stats */
                            const REAL8 rhomaxline,		/**< amplitude prior normalization for lines */
                            const REAL8Vector *lXvec,		/**< vector of single-detector prior line odds ratio, default to \f$ l^X=1 \f$ for all \f$ X \f$ if NULL */
                            const BOOLEAN useAllTerms		/**< only use leading term (FALSE) or all terms (TRUE) in log sum exp formula? */
                          )
{

  /* check input parameters and report errors */
  XLAL_CHECK_REAL4 ( TwoF && TwoFXvec && TwoFXvec->data, XLAL_EFAULT, "Empty TwoF or TwoFX pointer as input parameter." );
  UINT4 numDetectors = TwoFXvec->length;
  XLAL_CHECK_REAL4 ( !lXvec || ( lXvec->length == numDetectors ), XLAL_EBADLEN, "Input lX (%d) and TwoFX (%d) vectors have different lengths.", lXvec->length, numDetectors );

  XLAL_CHECK_REAL4 ( rhomaxline >= 0, XLAL_EDOM, "Negative prior range 'rhomaxline' = %g. Must be >= 0!", rhomaxline );
  REAL8 logRhoTerm = 0.0;
  if ( rhomaxline > 0.0 ) {
   logRhoTerm = 4.0 * log(rhomaxline) - log(70.0);
  }
  else { /* if rhomaxline == 0.0, logRhoTerm should become irrelevant in summation */
    logRhoTerm = - LAL_REAL8_MAX;
  }

  REAL8 *loglX = NULL;
  REAL8 loglXtemp[numDetectors];
  if ( lXvec ) {
    for (UINT4 X = 0; X < numDetectors; X++) {
      if ( lXvec->data[X] > 0 ) {
        loglXtemp[X] = log(lXvec->data[X]);
      }
      else if ( lXvec->data[X] == 0 ) { /* if zero prior ratio, approximate log(0)=-inf by -LAL_REAL4_MAX to avoid raising underflow exceptions */
        loglXtemp[X] = - LAL_REAL8_MAX;
      }
      else { /* negative prior ratio is a mistake! */
       XLAL_ERROR_REAL4 ( XLAL_EDOM, "Negative input prior-ratio for detector X=%d: lX[X]=%g\n", X, lXvec->data[X] );
      }
    }  /* for X < numDetectors */
    loglX = loglXtemp;
  } /* if lXvec */

  REAL4 LV = XLALComputeLineVetoArray ( TwoF, numDetectors, TwoFXvec->data, logRhoTerm, loglX, useAllTerms );

  return LV;

} /* XLALComputeLineVeto() */


/**
 * Deprecated function to compute Line Veto statistics from multi- and single-detector \f$ \mathcal{F} \f$-stats,
 * using outdated \f$ \rho \f$ notation and prior normalization.
 * This is not the log-Bayes-factor!
 * \f$ \mathrm{LV} = \mathcal{F} - \log \left( \frac{\rho_{\mathrm{max,line}}^4}{70} + \sum_X l^X e^{\mathcal{F}^X} \right) \f$
 *
 * \deprecated use XLALComputeBSGL() instead.
 *
 * Implemented by log sum exp formula:
 * \f$ \mathrm{LV} = \mathcal{F} - \max(\mathrm{denom. terms}) - \log \left( \sum e^{\mathrm{denom. terms}-\max} \right) \f$.
 *
 * From the analytical derivation, there should be an extra term \f$ + o_{SN} + 4 \log \left( \frac{\rho_{\mathrm{max,line}}}{\rho_{\mathrm{max,sig}}} \right) \f$,
 * but this is irrelevant for toplist sorting, only a normalization which can be replaced arbitrarily.
 *
 * NOTE: priors logRhoTerm, loglX have to be logarithmized already.
 */
REAL4
XLALComputeLineVetoArray ( const REAL4 TwoF,		/**< multi-detector \f$ \mathcal{F} \f$-stat */
                           const UINT4 numDetectors,	/**< number of detectors */
                           const REAL4 *TwoFX,		/**< array of single-detector \f$ \mathcal{F} \f$-stats */
                           const REAL8 logRhoTerm,	/**< extra term coming from prior normalization: \f$ \log \left( \frac{\rho_{\mathrm{max,line}}^4}{70} \right) \f$ */
                           const REAL8 *loglX,		/**< array of logs of single-detector prior line odds ratios, default to \f$ \log(l^X)=\log(1)=0 \f$ for all \f$ X \f$ if NULL */
                           const BOOLEAN useAllTerms	/**< only use leading term (FALSE) or all terms (TRUE) in log sum exp formula? */
                         )
{

  /* check input parameters and report errors */
  XLAL_CHECK_REAL4 ( TwoF && TwoFX, XLAL_EFAULT, "Empty TwoF or TwoFX pointer as input parameter." );

  /* set up temporary variables and structs */
  REAL4 log0  = - LAL_REAL4_MAX;	/* approximates -inf */

  REAL4 maxInSum = log0;           /* keep track of largest summand in denominator, for logsumexp formula below */
  REAL4 FXprior[numDetectors];     /* FXprior equiv log(lX * e^(FX)) = FX + loglX */

  for (UINT4 X = 0; X < numDetectors; X++) {
    FXprior[X] = 0.5 * TwoFX[X];
    if (loglX) { /* if no priors given, just use lX=1 => loglX=0 for all X => do not add anything */
      FXprior[X] += loglX[X];
    }
    /* keep track of maximum value in denominator sum  */
    if ( FXprior[X] > maxInSum ) {
      maxInSum = FXprior[X];
    }
  } /* for X < numDetectors */

  /* special treatment for additional denominator term 'rho^4/70' */
  if ( logRhoTerm > maxInSum ) {
    maxInSum = logRhoTerm;
  }

  REAL4 LV = 0.0;	/* output variable for Line Veto statistics */

  LV = 0.5 * TwoF - maxInSum;	/* dominant term to LV-statistic */

  if ( useAllTerms ) { /* optionally add logsumexp term (possibly negligible in many cases) */

    REAL4 extraSum=0;	/* will be:  e^[-(maxInSum - logRhoTerm)] + sum_X e^[ -(maxInSum - FXprior) ] >= 1 */

    /* need to treat (rho^4/70) term separately */
    extraSum += exp ( logRhoTerm - maxInSum );

    /* now add all FX-contributions */
    for (UINT4 X = 0; X < numDetectors; X++) {
      extraSum += exp ( FXprior[X] - maxInSum );
    }

    LV -= log ( extraSum );

  } /* if useAllTerms */

  return LV;

} /* XLALComputeLineVetoArray() */
