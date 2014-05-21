/*
 * Copyright (C) 2011 David Keitel
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
XLALCompareLVComputations ( const REAL4 TwoF,
                            const REAL4Vector *TwoFX,
                            const REAL8 rhomaxline,
                            const REAL8Vector *lX,
                            const REAL4 tolerance_allterms,
                            const REAL4 tolerance_leadterm );

/* ###################################  MAIN  ################################### */

int main( int argc, char *argv[]) {


  /* sanity check for input arguments */
  if ( argc != 1 )
    XLAL_ERROR ( XLAL_EINVAL, "The executable '%s' doesn't support any input arguments right now.\n", argv[0] );

  printf ("Starting test...\n");

  /* set up single- and multi-IFO F-stat input */
  REAL4 TwoF = 7.0;
  UINT4 numDetectors = 2;
  REAL4Vector *TwoFX = NULL;
  if ( (TwoFX = XLALCreateREAL4Vector ( numDetectors )) == NULL ) {
    XLAL_ERROR ( XLAL_EFUNC, "failed to XLALCreateREAL4Vector( %d )\n", numDetectors );
    return XLAL_EFAILED;
  }
  TwoFX->data[0] = 4.0;
  TwoFX->data[1] = 12.0;
  REAL8 rhomaxline = 0.0; /* prior from LV-stat derivation, 0 means pure line veto, +inf means pure multi-Fstat */
  REAL8Vector *lX = NULL; /* per-IFO prior odds ratio for line vs. Gaussian noise, NULL is interpreted as l[X]=1 for all X */

  /* maximum allowed difference between recalculated and XLAL result */
  REAL4 tolerance_allterms = 2e-04;
  REAL4 tolerance_leadterm = 2e-02;

  /* compute and compare the results for one set of rhomaxline, lX values */
  printf ("Computing LV-stat for TwoF_multi=%f, TwoFX[0]=%f, TwoFX[1]=%f, rhomaxline=%f, priors lX=NULL...\n", TwoF, TwoFX->data[0], TwoFX->data[1], rhomaxline );
  if ( XLALCompareLVComputations( TwoF, TwoFX, rhomaxline, lX, tolerance_allterms, tolerance_leadterm ) != XLAL_SUCCESS ) {
    XLAL_ERROR ( XLAL_EFUNC, "Test failed.\n" );
    return XLAL_EFAILED;
  }

  /* change the priors to catch more possible problems */
  rhomaxline = 5.0;
  if ( (lX = XLALCreateREAL8Vector ( numDetectors )) == NULL ) {
    XLAL_ERROR ( XLAL_EFUNC, "failed to XLALCreateREAL8Vector( %d )\n", numDetectors );
    return XLAL_EFAILED;
  }
  lX->data[0] = 0.5;
  lX->data[1] = 0.8;

  /* compute and compare the results for second set of rhomaxline, lX values */
  printf ("Computing LV-stat for TwoF_multi=%f, TwoFX[0]=%f, TwoFX[1]=%f, rhomaxline=%f, priors lX=(%f,%f)...\n", TwoF, TwoFX->data[0], TwoFX->data[1], rhomaxline, lX->data[0], lX->data[1] );
  if ( XLALCompareLVComputations( TwoF, TwoFX, rhomaxline, lX, tolerance_allterms, tolerance_leadterm ) != XLAL_SUCCESS ) {
    XLAL_ERROR ( XLAL_EFUNC, "Test failed.\n" );
    return XLAL_EFAILED;
  }

  /* free memory */
  XLALDestroyREAL4Vector(TwoFX);
  XLALDestroyREAL8Vector(lX);

  LALCheckMemoryLeaks();

  return XLAL_SUCCESS;
} /* main */


/**
 * Test function to compute LV-stat values both from scratch and by XLALComputeLineVeto,
 * compare the results and exit if tolerance is violated.
 */
int
XLALCompareLVComputations ( const REAL4 TwoF,          /**< multi-detector  Fstat */
                            const REAL4Vector *TwoFX,  /**< vector of single-detector Fstats */
                            const REAL8 rhomaxline,    /**< amplitude prior normalization for lines */
                            const REAL8Vector *lX,     /**< vector of single-detector prior line odds ratio, default to lX=1 for all X if NULL */
                            const REAL4 tolerance_allterms, /**< tolerance for useAllTerms=TRUE */
                            const REAL4 tolerance_leadterm  /**< tolerance for useAllTerms=FALSE (usually higher) */
                          )
{

  /* compute LV-stat "the pedestrian way", explicit formula for numDet=2 */
  REAL4 LV_extcomp = 0.5*TwoF;
  if ( lX )
    LV_extcomp -= log( pow(rhomaxline,4)/70.0 + lX->data[0]*exp(0.5*TwoFX->data[0]) + lX->data[1]*exp(0.5*TwoFX->data[1]) );
  else /* lX=NULL is interpreted as l[X]=1 for all X */
    LV_extcomp -= log( pow(rhomaxline,4)/70.0 + exp(0.5*TwoFX->data[0]) + exp(0.5*TwoFX->data[1]) );

  /* faster version: use only the leading term of the LV denominator sum */
  BOOLEAN useAllTerms = FALSE;
  xlalErrno = 0;
  REAL4 LV_XLAL_leadterm = XLALComputeLineVeto ( TwoF, TwoFX, rhomaxline, lX, useAllTerms );
  if ( xlalErrno != 0 ) {
    XLAL_ERROR ( XLAL_EFUNC, "XLALComputeLineVeto() failed with xlalErrno = %d\n", xlalErrno );
    return XLAL_FAILURE;
  }

  /* more precise version: use all terms of the LV denominator sum */
  useAllTerms = TRUE;
  xlalErrno = 0;
  REAL4 LV_XLAL_allterms = XLALComputeLineVeto ( TwoF, TwoFX, rhomaxline, lX, useAllTerms );
  if ( xlalErrno != 0 ) {
    XLAL_ERROR ( XLAL_EFUNC, "XLALComputeLineVeto() failed with xlalErrno = %d\n", xlalErrno );
    return XLAL_FAILURE;
  }

  /* compute relative deviations */
  REAL4 diff_allterms = fabs( LV_XLAL_allterms - LV_extcomp ) / ( 0.5 * ( LV_XLAL_allterms + LV_extcomp ));
  REAL4 diff_leadterm = fabs( LV_XLAL_leadterm - LV_extcomp ) / ( 0.5 * ( LV_XLAL_leadterm + LV_extcomp ));

  /* output results and deviations and return with error when tolerances are violated */
  printf ("Externally recomputed             : LV=%f\n", LV_extcomp);
  printf ("XLALComputeLineVeto with allterms : LV=%f (rel. dev.: %f)", LV_XLAL_allterms, diff_allterms);
  if ( fabs(diff_allterms) <= tolerance_allterms )
    printf (" ==> OK!\n");
  else {
    printf (" ==> BAD!\n");
    XLAL_ERROR ( XLAL_EFAILED, "\nTolerance %f exceeded for useAllTerms=TRUE!\n", tolerance_allterms );
    return XLAL_FAILURE;
  }
  printf ("XLALComputeLineVeto with !allterms: LV=%f (rel. dev.: %f)", LV_XLAL_leadterm, diff_leadterm);
  if ( fabs(diff_leadterm) <= tolerance_leadterm )
    printf (" ==> OK!\n");
  else {
    printf (" ==> BAD!\n");
    XLAL_ERROR ( XLAL_EFAILED, "\nTolerance %f exceeded for useAllTerms=FALSE!\n", tolerance_leadterm );
    return XLAL_FAILURE;
  }

  return XLAL_SUCCESS;

} /* XLALStringVector_TEST() */
