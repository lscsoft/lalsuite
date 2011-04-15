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

/**
 * \author David Keitel
 * \date 2011
 * \ingroup pulsarCoherent
 * \file
 * \brief Header-file defining functions related to GCT Line Veto followups
 *
 * This code is partly based on work done by
 * Reinhard Prix, Maria Alessandra Papa, M. Siddiqi
 *
 * $Id$
 *
 */

/*---------- INCLUDES ----------*/
#define __USE_ISOC99 1
#include "LineVeto.h"

NRCSID( LINEVETOC, "$Id$");

/*---------- local DEFINES ----------*/
#define TRUE (1==1)
#define FALSE (1==0)

/*----- Macros ----- */
#define INIT_MEM(x) memset(&(x), 0, sizeof((x)))
#define SQUARE(x) ( (x) * (x) )

/*----- SWITCHES -----*/

/*---------- internal types ----------*/

/*---------- Global variables ----------*/

/* empty initializers  */
const LVcomponents empty_LVcomponents;

/*---------- internal prototypes ----------*/

/*==================== FUNCTION DEFINITIONS ====================*/


/** XLAL function to go through toplist and compute Line Veto statistics for each candidate */
int XLALComputeExtraStatsForToplist ( toplist_t *list,                                        /**< list of cancidates with f, sky position etc. - no output so far */
				      const MultiSFTVectorSequence *multiSFTsV,               /**< data files (SFTs) for all detectors and segments */
				      const MultiNoiseWeightsSequence *multiNoiseWeightsV,    /**< noise weights for all detectors and segments */
				      const MultiDetectorStateSeriesSequence *multiDetStatesV,/**< some state info for all detectors */
				      const ComputeFParams *CFparams,                         /**< additional parameters needed for ComputeFStat */
				      const LIGOTimeGPS refTimeGPS,                           /**< reference time, needed for ExtrapolatePulsarSpins */
				      const LIGOTimeGPS tMidGPS,                              /**< reference time, needed for ExtrapolatePulsarSpins */
				      const BOOLEAN SignalOnly                                /**< flag for case with no noise, makes some extra checks necessary */
				    )
{
  const char *fn = __func__;

  /* check input parameters and report errors */
  if ( !list || !multiSFTsV || !multiNoiseWeightsV || !multiDetStatesV || !CFparams ) {
    XLALPrintError ("\nError in function %s, line %d : Empty pointer as input parameter!\n\n", fn, __LINE__);
    XLAL_ERROR ( fn, XLAL_EFAULT);
  }

  if ( !list->data || !list->heap ) {
    XLALPrintError ("\nError in function %s, line %d : Input toplist has no elements!\n\n", fn, __LINE__);
    XLAL_ERROR ( fn, XLAL_EFAULT);
  }

  if ( list->elems == 0 ) {
    XLALPrintError ("\nError in function %s, line %d : Input toplist has zero length!\n\n", fn, __LINE__);
    XLAL_ERROR ( fn, XLAL_EBADLEN );
  }

  if ( !multiSFTsV->data[0] || (multiSFTsV->data[0]->length == 0) ) {
    XLALPrintError ("\nError in function %s, line %d : Input multiSFT vector has no elements!\n\n", fn, __LINE__);
    XLAL_ERROR ( fn, XLAL_EFAULT);
  }

  /* set up temporary variables and structs */
  PulsarDopplerParams candidateDopplerParams = empty_PulsarDopplerParams; /* struct containing sky position, frequency and fdot for the current candidate */
  PulsarSpins fkdotTMP; /* temporary spin parameters for XLALExtrapolatePulsarSpins */
  REAL8 deltaTau;       /* temporary variable to convert LIGOTimeGPS into real number difference for XLALExtrapolatePulsarSpins */

  /* initialize doppler parameters */
  candidateDopplerParams.refTime = tMidGPS;  /* spin parameters will be given to ComputeFStat at this refTime */
  INIT_MEM( fkdotTMP );
  deltaTau = XLALGPSDiff( &candidateDopplerParams.refTime, &refTimeGPS );

  /* initialise LVcomponents structure and allocate memory */
  UINT4          numDetectors = multiSFTsV->data[0]->length;  /* number of different single detectors */
  LVcomponents   lineVeto;      /* struct containing multi-detector Fstat, single-detector Fstats, Line Veto stat */
  if ( (lineVeto.TwoFX = XLALCreateREAL8Vector ( numDetectors )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCreateREAL8Vector( %d )\n", fn, numDetectors );
    XLAL_ERROR ( fn, XLAL_EFUNC );
  }

  UINT4 j;
  UINT4 numElements = list->elems;
  /* loop over toplist: re-compute sumTwoF and sumTwoFX for all candidates */
  for (j = 0; j < numElements; j++ )
    {
      GCTtopOutputEntry *elem = toplist_elem ( list, j );

      if ( (elem->sumTwoFX = XLALCreateREAL4Vector ( numDetectors )) == NULL ) {
        XLALPrintError ("%s: failed to XLALCreateREAL4Vector( %d )\n", fn, numDetectors );
        XLAL_ERROR ( fn, XLAL_EFUNC );
      }

      /* get frequency, sky position, doppler parameters from toplist candidate and save to dopplerParams */
      candidateDopplerParams.Alpha = elem->Alpha;
      candidateDopplerParams.Delta = elem->Delta;
      fkdotTMP[0] = elem->Freq;
      fkdotTMP[1] = elem->F1dot;

      /* extrapolate pulsar spins to correct time (more stable against large deltaTau than directly resetting refTime) */
      if ( XLALExtrapolatePulsarSpins( candidateDopplerParams.fkdot, fkdotTMP, deltaTau ) != XLAL_SUCCESS ) {
        XLALPrintError ("\n%s, line %d : XLALExtrapolatePulsarSpins() failed.\n\n", fn, __LINE__);
        XLAL_ERROR ( fn, XLAL_EFUNC );
      }

      /*  recalculate multi- and single-IFO Fstats for all segments for this candidate */
      XLALComputeExtraStatsSemiCoherent( &lineVeto, &candidateDopplerParams, multiSFTsV, multiNoiseWeightsV, multiDetStatesV, CFparams, SignalOnly );
      if ( xlalErrno != 0 ) {
        XLALPrintError ("\nError in function %s, line %d : Failed call to XLALComputeLineVetoSemiCoherent().\n\n", fn, __LINE__);
        XLAL_ERROR ( fn, XLAL_EFUNC );
      }

      /* save values in toplist */
      elem->sumTwoFnew         = lineVeto.TwoF;
      UINT4 X;
      for ( X = 0; X < numDetectors; X ++ )
        elem->sumTwoFX->data[X]  = lineVeto.TwoFX->data[X];

    } /* for j < numElements */

  /* free temporary structures */
  XLALDestroyREAL8Vector ( lineVeto.TwoFX );

  return (XLAL_SUCCESS);

} /* XLALComputeExtraStatsForToplist() */






/** XLAL Function to recalculate single-IFO Fstats for all semicoherent search segments, and use them to compute Line Veto statistics
*/
int XLALComputeExtraStatsSemiCoherent ( LVcomponents *lineVeto,                                 /**< [out] structure containing multi TwoF, single TwoF, LV stat */
					const PulsarDopplerParams *dopplerParams,               /**< sky position, frequency and fdot for a given candidate */
					const MultiSFTVectorSequence *multiSFTsV,               /**< data files (SFTs) for all detectors and segments */
					const MultiNoiseWeightsSequence *multiNoiseWeightsV,    /**< noise weights for all detectors and segments */
					const MultiDetectorStateSeriesSequence *multiDetStatesV,/**< some state info for all detectors */
					const ComputeFParams *CFparams,                         /**< additional parameters needed for ComputeFStat */
					const BOOLEAN SignalOnly                                /**< flag for case with no noise, makes some extra checks necessary */
				      )
{
  const char *fn = __func__;

  /* check input parameters and report errors */
  if ( !lineVeto || !lineVeto->TwoFX || !lineVeto->TwoFX->data || !dopplerParams || !multiSFTsV || !multiNoiseWeightsV || !multiDetStatesV || !CFparams ) {
    XLALPrintError ("\nError in function %s, line %d : Empty pointer as input parameter!\n\n", fn, __LINE__);
    XLAL_ERROR ( fn, XLAL_EFAULT);
  }
  if ( !multiSFTsV->data[0] ) {
    XLALPrintError ("\nError in function %s, line %d : Input multiSFT vector has no elements!\n\n", fn, __LINE__);
    XLAL_ERROR ( fn, XLAL_EFAULT);
  }
  if ( multiSFTsV->length == 0 ) {
    XLALPrintError ("\nError in function %s, line %d : Input multiSFT vector over segments has zero length!\n\n", fn, __LINE__);
    XLAL_ERROR ( fn, XLAL_EBADLEN);
  }

  UINT4 numSegments  = multiSFTsV->length;
  UINT4 numDetectors = multiSFTsV->data[0]->length;

  if ( lineVeto->TwoFX->length != numDetectors ) {
    XLALPrintError ("\%s, line %d : Inconsistent number of detector in TwoFX: %d, while multiSFTsV: %d!\n\n", fn, __LINE__, lineVeto->TwoFX->length, numDetectors );
    XLAL_ERROR ( fn, XLAL_EBADLEN );
  }

  /* fake LAL status structure, needed as long as ComputeFStat is LAL function and not XLAL */
  LALStatus fakeStatus = blank_status;
  Fcomponents    Fstat;                           /* temporary struct for ComputeFStat */

  /* temporary copy of Fstatistic parameters structure, needed to change returnAtoms for function scope only */
  ComputeFParams CFparams_internal = (*CFparams);
  CFparams_internal.returnAtoms   = TRUE;

  /* initialiase LVcomponents structure */
  lineVeto->TwoF = 0.0;
  lineVeto->LV   = 0.0;
  UINT4 X;
  for (X = 0; X < numDetectors; X++) {
    lineVeto->TwoFX->data[X] = 0.0;
  }

  REAL8 Tsft = 1.0 / multiSFTsV->data[0]->data[0]->data[0].deltaF;	/* get SFT duration */

  /* compute single- and multi-detector Fstats for each data segment and sum up */
  UINT4 k;
  for (k = 0; k < numSegments; k++)
    {
      /* recompute multi-detector Fstat and atoms */
      fakeStatus = blank_status;
      ComputeFStat ( &fakeStatus, &Fstat, dopplerParams, multiSFTsV->data[k], multiNoiseWeightsV->data[k], multiDetStatesV->data[k], &CFparams_internal, NULL );
      if ( fakeStatus.statusCode ) {
        XLALPrintError ("\%s, line %d : Failed call to LAL function ComputeFStat(). statusCode=%d\n\n", fn, __LINE__, fakeStatus.statusCode);
        XLAL_ERROR ( fn, XLAL_EFUNC );
      }

      if ( SignalOnly ) {      /* normalization factor correction */
        Fstat.F *= 2.0 / Tsft;
        Fstat.F += 2;
      }

      lineVeto->TwoF  += 2.0 * Fstat.F; /* sum up multi-detector Fstat for this segment*/

      /* recompute single-detector Fstats from atoms */
      for (X = 0; X < numDetectors; X++)
        {
          REAL8 twoFX = 2.0 * XLALComputeFstatFromAtoms ( &Fstat, X );
          if ( xlalErrno != 0 ) {
	    XLALPrintError ("\nError in function %s, line %d : Failed call to XLALComputeFstatFromAtoms().\n\n", fn, __LINE__);
	    XLAL_ERROR ( fn, XLAL_EFUNC );
	  }

          if ( SignalOnly ) {                      /* normalization factor correction (TwoF=2.0*F has been done before, this time!) */
            twoFX *= 4.0 / Tsft;
            twoFX += 4;
          }

          lineVeto->TwoFX->data[X]  += twoFX; /* sum up single-detector Fstat for this segment*/

        } /* for X < numDetectors X */

      /* free memory for atoms that was allocated within ComputeFStat  */
      XLALDestroyMultiFstatAtomVector ( Fstat.multiFstatAtoms );

    } /* for k < numSegments */

  /* get average stats over all segments */
  lineVeto->TwoF /= numSegments;
  for (X = 0; X < numDetectors; X++) {
    lineVeto->TwoFX->data[X] /= numSegments;
  }

  return(XLAL_SUCCESS);

} /* XLALComputeExtraStatsSemiCoherent() */


/** XLAL function to compute single-IFO Fstat from multi-IFO Atoms: */
REAL8 XLALComputeFstatFromAtoms ( const Fcomponents *Fstat,   /**< multi-detector Fstat */
				  const UINT4       X        /**< detector number */
				  )
{
  const char *fn = __func__;

  /* check input parameters and report errors */
  if ( !Fstat || !Fstat->multiFstatAtoms || !Fstat->multiFstatAtoms->data || !Fstat->multiFstatAtoms->data[0]->data ) {
    XLALPrintError ("\nError in function %s, line %d : Empty pointer as input parameter!\n\n", fn, __LINE__);
    XLAL_ERROR ( fn, XLAL_EFAULT);
  }

  if ( Fstat->multiFstatAtoms->length == 0 ) {
    XLALPrintError ("\nError in function %s, line %d : Input MultiFstatAtomVector has zero length!\n\n", fn, __LINE__);
    XLAL_ERROR ( fn, XLAL_EBADLEN );
  }

  if ( X > Fstat->multiFstatAtoms->length-1 ) {
    XLALPrintError ("\nError in function %s, line %d : Invalid detector number!\n\n", fn, __LINE__);
    XLAL_ERROR ( fn, XLAL_EDOM );
  }

  if ( Fstat->multiFstatAtoms->data[X]->length == 0 ) {
    XLALPrintError ("\nError in function %s, line %d : Input FstatAtomVector has zero length!\n\n", fn, __LINE__);
    XLAL_ERROR ( fn, XLAL_EDOM );
  }

  /* set up temporary variables and structs */
  REAL8 mmatrixA = 0.0, mmatrixB = 0.0, mmatrixC = 0.0;
  REAL8 FX = 0.0;
  COMPLEX8 Fa, Fb;
  UINT4 alpha;
  UINT4 numSFTs = Fstat->multiFstatAtoms->data[X]->length;

  /* sum up matrix elements and Fa, Fb */
  Fa.re = 0.0;
  Fa.im = 0.0;
  Fb.re = 0.0;
  Fb.im = 0.0;
  for ( alpha = 0; alpha < numSFTs; alpha++)
    {
      FstatAtom *thisAtom = &Fstat->multiFstatAtoms->data[X]->data[alpha];

      mmatrixA += thisAtom->a2_alpha;
      mmatrixB += thisAtom->b2_alpha;
      mmatrixC += thisAtom->ab_alpha;
      Fa.re    += thisAtom->Fa_alpha.re;
      Fa.im    += thisAtom->Fa_alpha.im;
      Fb.re    += thisAtom->Fb_alpha.re;
      Fb.im    += thisAtom->Fb_alpha.im;

    } /* for alpha < numSFTs */

  /* compute determinant and final Fstat (not twoF!) */
  REAL8 Dinv = 1.0 / ( mmatrixA * mmatrixB - SQUARE(mmatrixC) );
  FX = Dinv * ( mmatrixB * ( SQUARE(Fa.re) + SQUARE(Fa.im) ) + mmatrixA * ( SQUARE(Fb.re) + SQUARE(Fb.im) ) - 2.0 * mmatrixC * (Fa.re*Fb.re + Fa.im*Fb.im) );

  return FX;

} /* XLALComputeFstatFromAtoms() */




/** XLAL function to compute Line Veto statistics from multi- and single-detector Fstats:
 *  LV = log( e^2F / sum(e^2FX) )
 *  implemented by log sum exp formula:
 *  LV = 2F - max(2FX) - log( sum(e^(2FX-max(2FX))) )
*/
REAL8 XLALComputeLineVeto ( const REAL8 TwoF,          /**< multi-detector  Fstat */
                            const REAL8Vector *TwoFX,  /**< vector of single-detector Fstats */
                            const REAL8 rhomax,        /**< amplitude prior normalization, necessary for signal-noise veto, set to 0 for pure signal-line veto */
                            const REAL8Vector *priorX  /**< vector of single-detector prior line odds ratio, set all to 1/numDetectors for neutral analysis */
		          )
{
  const char *fn = __func__;

  /* check input parameters and report errors */
  if ( !TwoF || !TwoFX || !TwoFX->data || !priorX || !priorX->data ) {
    XLALPrintError ("\nError in function %s, line %d : Empty pointer as input parameter!\n\n", fn, __LINE__);
    XLAL_ERROR_REAL8 ( fn, XLAL_EFAULT);
  }

  if ( TwoFX->length == 0 ) {
    XLALPrintError ("\nError in function %s, line %d :Input TwoFX vector has zero length!\n\n", fn, __LINE__);
    XLAL_ERROR_REAL8 ( fn, XLAL_EBADLEN);
  }

  if ( priorX->length == 0 ) {
    XLALPrintError ("\nError in function %s, line %d :Input priorX vector has zero length!\n\n", fn, __LINE__);
    XLAL_ERROR_REAL8 ( fn, XLAL_EBADLEN);
  }

  /* set up temporary variables and structs */
  UINT4 X;                            /* loop summation variable */
  UINT4 numDetectors = TwoFX->length; /* loop summation upper limit for detectors */
  REAL8 maxSum = -1000.0;             /* maximum of terms in denominator, for logsumexp formula */
  REAL8 LV = 0.0;                     /* output variable for Line Veto statistics */

  if ( rhomax < 0.0 ) {
    XLALPrintError ("\nError in function %s, line %d : nonpositive input rhomax!\n\n", fn, __LINE__);
    XLAL_ERROR_REAL8 ( fn, XLAL_EFPINVAL);
  }
  /* for rhomax = 0.0, just ignore in summation */
  if ( rhomax > 0.0 ) {
    maxSum = log(rhomax);
  }

  for (X = 0; X < numDetectors; X++) {
    if ( TwoFX->data[X] + log(priorX->data[X]) > maxSum ) {
      maxSum = TwoFX->data[X] + log(priorX->data[X]);
    }
  }

  /* logsumexp formula */
  if ( rhomax > 0.0 )   LV =  exp( log(rhomax) - maxSum );
  for (X = 0; X < numDetectors; X++) {
    LV += exp( TwoFX->data[X]  + log(priorX->data[X]) - maxSum );
  }
  if ( LV <= 0 ) { /* return error code for log (0) */
    XLALPrintError ("\nError in function %s, line %d : log(nonpositive) in LV denominator. \n\n", fn, __LINE__);
    XLAL_ERROR_REAL8 ( fn, XLAL_EFPINVAL );
  }
  else {
    LV = TwoF - maxSum - log( LV );
    return(LV);
  }


} /* XLALComputeLineVeto() */
