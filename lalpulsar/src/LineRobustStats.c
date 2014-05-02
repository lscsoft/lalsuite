/*
 *  Copyright (C) 2011-2014 David Keitel
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

/*---------- INCLUDES ----------*/
#define __USE_ISOC99 1
#include "LineRobustStats.h"

/*---------- local DEFINES ----------*/
#define TRUE (1==1)
#define FALSE (1==0)

/*----- Macros ----- */

/*----- SWITCHES -----*/

/*---------- internal types ----------*/

/*---------- Global variables ----------*/

/*---------- internal prototypes ----------*/

/*==================== FUNCTION DEFINITIONS ====================*/


/**
 * XLAL function to compute Line Veto statistics from multi- and single-detector Fstats:
 * this is now a wrapper for XLALComputeLineVetoArray which just translates REAL4Vectors to fixed REAL4 arrays
 * and linear to logarithmic priors rhomaxline and lX
 * NOTE: if many LV values at identical priors are required, directly call XLALComputeLineVetoArray for better performance
 */
REAL4 XLALComputeLineVeto ( const REAL4 TwoF,          /**< multi-detector  Fstat */
                            const REAL4Vector *TwoFXvec,  /**< vector of single-detector Fstats */
                            const REAL8 rhomaxline,    /**< amplitude prior normalization for lines */
                            const REAL8Vector *lXvec, /**< vector of single-detector prior line odds ratio, default to lX=1 for all X if NULL */
                            const BOOLEAN useAllTerms  /**< only use leading term (FALSE) or all terms (TRUE) in log sum exp formula? */
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
      else if ( lXvec->data[X] == 0 ) { /* if zero prior ratio, approximate log(0)=-inf by -LAL_REA4_MAX to avoid raising underflow exceptions */
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
 * XLAL function to compute Line Veto statistics from multi- and single-detector Fstats:
 * LV = F - log ( rhomaxline^4/70 + sum(e^FX) )
 * implemented by log sum exp formula:
 * LV = F - max(denom_terms) - log( sum(e^(denom_term-max)) )
 * from the analytical derivation, there should be a term LV += O_SN^0 + 4.0*log(rhomaxline/rhomaxsig)
 * but this is irrelevant for toplist sorting, only a normalization which can be replaced arbitrarily
 * NOTE: priors logRhoTerm, loglX have to be logarithmized already
 */
REAL4
XLALComputeLineVetoArray ( const REAL4 TwoF,   /**< multi-detector Fstat */
                           const UINT4 numDetectors, /**< number of detectors */
                           const REAL4 *TwoFX,       /**< array of single-detector Fstats */
                           const REAL8 logRhoTerm,   /**< extra term coming from prior normalization: log(rho_max_line^4/70) */
                           const REAL8 *loglX,       /**< array of logs of single-detector prior line odds ratios, default to loglX=log(1)=0 for all X if NULL */
                           const BOOLEAN useAllTerms /**< only use leading term (FALSE) or all terms (TRUE) in log sum exp formula? */
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
