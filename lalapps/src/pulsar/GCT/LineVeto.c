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

/*----- SWITCHES -----*/

/*---------- internal types ----------*/

/*---------- Global variables ----------*/

/* empty initializers  */
const LVcomponents empty_LVcomponents;

/*---------- internal prototypes ----------*/

/*==================== FUNCTION DEFINITIONS ====================*/


/** XLAL function to go through toplist and compute Line Veto statistics for each candidate */
int XLALComputeExtraStatsForToplist ( toplist_t *list,                                        /**< list of cancidates with f, sky position etc. - no output so far */
				      const MultiSFTVectorSequence *multiSFTs,                /**< data files (SFTs) for all detectors and segments */
				      const MultiNoiseWeightsSequence *multiNoiseWeights,     /**< noise weights for all detectors and segments */
				      const MultiDetectorStateSeriesSequence *multiDetStates, /**< some state info for all detectors */
				      const ComputeFParams *CFparams,                         /**< additional parameters needed for ComputeFStat */
				      const LIGOTimeGPS refTimeGPS,                           /**< reference time, needed for ExtrapolatePulsarSpins */
				      const LIGOTimeGPS tMidGPS,                              /**< reference time, needed for ExtrapolatePulsarSpins */
				      const BOOLEAN SignalOnly                                /**< flag for case with no noise, makes some extra checks necessary */
				    )
{
  const char *fn = __func__;

  /* set up temporary variables and structs */
  UINT4          j, X;          /* loop counting variables */
  UINT4          numDetectors;  /* number of different single detectors */
  LVcomponents   lineVeto;      /* struct containing multi-detector Fstat, single-detector Fstats, Line Veto stat */
  PulsarDopplerParams candidateDopplerParams; /* struct containing sky position, frequency and fdot for the current candidate */
  PulsarSpins fkdotTMP;
  REAL8 deltaTau;  /* temporary variable to convert LIGOTimeGPS into real number difference for XLALExtrapolatePulsarSpins */

  /* check input parameters and report errors */
  if ( !list || !multiSFTs || !multiNoiseWeights || !multiDetStates || !CFparams ) {
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

  /* initialize doppler parameters */
  candidateDopplerParams.orbit = NULL;
  INIT_MEM ( candidateDopplerParams.fkdot );
  candidateDopplerParams.refTime = tMidGPS;
  INIT_MEM( fkdotTMP );
  deltaTau = XLALGPSDiff( &candidateDopplerParams.refTime, &refTimeGPS );

  /* initialise LVcomponents structure and allocate memory */
  if ( !multiSFTs->data[0] ) {
    XLALPrintError ("\nError in function %s, line %d : Input multiSFT vector has no elements!\n\n", fn, __LINE__);
    XLAL_ERROR ( fn, XLAL_EFAULT);
  }
  if ( multiSFTs->data[0]->length == 0 ) {
    XLALPrintError ("\nError in function %s, line %d : Input multiSFT vector over detectors has zero length!\n\n", fn, __LINE__);
    XLAL_ERROR ( fn, XLAL_EBADLEN);
  }

  numDetectors = multiSFTs->data[0]->length;
  lineVeto.TwoFX = XLALCreateREAL8Vector ( numDetectors );

  /* allocate FX vectors in toplist */
  for (j = 0; j < list->elems; j++ ) { /* loop over toplist elements */
    (*(GCTtopOutputEntry*)list->heap[j]).sumTwoFX = XLALCreateREAL4Vector ( numDetectors );
    (*(GCTtopOutputEntry*)list->heap[j]).sumTwoFnew = 0.0;
    for (X = 0; X < numDetectors; X++ ) {
      (*(GCTtopOutputEntry*)list->heap[j]).sumTwoFX->data[X] = 0.0;
    }
  }

  for (j = 0; j < list->elems; j++ ) { /* loop over toplist elements */

    /* get frequency, sky position, doppler parameters from toplist candidate and save to dopplerParams */
    candidateDopplerParams.Alpha = (*(GCTtopOutputEntry*)list->heap[j]).Alpha;
    candidateDopplerParams.Delta = (*(GCTtopOutputEntry*)list->heap[j]).Delta;
    fkdotTMP[0] = (*(GCTtopOutputEntry*)list->heap[j]).Freq;
    fkdotTMP[1] = (*(GCTtopOutputEntry*)list->heap[j]).F1dot;

    /* extrapolate pulsar spins to correct time (really constant for all candidates and segments?) */
    xlalErrno = 0;
    XLALExtrapolatePulsarSpins( candidateDopplerParams.fkdot, fkdotTMP, deltaTau );
    if ( xlalErrno != 0 ) {
      XLALPrintError ("\nError in function %s, line %d : Failed call to XLALExtrapolatePulsarSpins.\n\n", fn, __LINE__);
      XLAL_ERROR ( fn, XLAL_EFUNC );
    }

    /* compute Line Veto statistic for this candidate by recalculating single-IFO Fstats for all segments */
    xlalErrno = 0;
    XLALComputeExtraStatsSemiCoherent( &lineVeto, &candidateDopplerParams, multiSFTs, multiNoiseWeights, multiDetStates, CFparams, SignalOnly );
    if ( xlalErrno != 0 ) {
      XLALPrintError ("\nError in function %s, line %d : Failed call to XLALComputeLineVetoSemiCoherent_v2.\n\n", fn, __LINE__);
      XLAL_ERROR ( fn, XLAL_EFUNC );
    }

    /* save values in toplist */
    (*(GCTtopOutputEntry*)list->heap[j]).sumTwoFnew         = lineVeto.TwoF;
    (*(GCTtopOutputEntry*)list->heap[j]).sumTwoFX->data[0]  = lineVeto.TwoFX->data[0];
    (*(GCTtopOutputEntry*)list->heap[j]).sumTwoFX->data[1]  = lineVeto.TwoFX->data[1];

  }

  /* free temporary structures */
  XLALDestroyREAL8Vector ( lineVeto.TwoFX );

  return (XLAL_SUCCESS);

} /* XLALComputeExtraStatsForToplist() */






/** XLAL Function to recalculate single-IFO Fstats for all semicoherent search segments, and use them to compute Line Veto statistics
*/
int XLALComputeExtraStatsSemiCoherent ( LVcomponents *lineVeto,                                 /**< [out] structure containing multi TwoF, single TwoF, LV stat */
					const PulsarDopplerParams *dopplerParams,               /**< sky position, frequency and fdot for a given candidate */
					const MultiSFTVectorSequence *multiSFTs,                /**< data files (SFTs) for all detectors and segments */
					const MultiNoiseWeightsSequence *multiNoiseWeights,     /**< noise weights for all detectors and segments */
					const MultiDetectorStateSeriesSequence *multiDetStates, /**< some state info for all detectors */
					const ComputeFParams *CFparams,                         /**< additional parameters needed for ComputeFStat */
					const BOOLEAN SignalOnly                                /**< flag for case with no noise, makes some extra checks necessary */
				      )
{
  const char *fn = __func__;

  /* check input parameters and report errors */
  if ( !lineVeto || !lineVeto->TwoFX->data || !dopplerParams || !multiSFTs || !multiNoiseWeights || !multiDetStates || !CFparams ) {
    XLALPrintError ("\nError in function %s, line %d : Empty pointer as input parameter!\n\n", fn, __LINE__);
    XLAL_ERROR ( fn, XLAL_EFAULT);
  }

  if ( lineVeto->TwoFX->length == 0 ) {
    XLALPrintError ("\nError in function %s, line %d : Input TwoFX vector has zero length!\n\n", fn, __LINE__);
    XLAL_ERROR ( fn, XLAL_EBADLEN );
  }

  /* set up temporary variables and structs */
  UINT4 k, X, numSegments, numDetectors;          /* loop counting variables and upper limits for segments and detectors */
  Fcomponents    Fstat;                           /* temporary struct for ComputeFStat */
  LVcomponents   TempStats;                       /* LV struct for temporary stats inside loops */
  REAL8          Tsft;                            /* SFT duration */

  /* fake LAL status structure, needed as long as ComputeFStat is LAL function and not XLAL */
  LALStatus fakeStatus = blank_status;
  ComputeFBuffer CFbuffer = empty_ComputeFBuffer;

  /* temporary copy of Fstatistic parameters structure, needed to change returnAtoms for function scope only */
  ComputeFParams CFparams_internal;
  CFparams_internal.Dterms        = CFparams->Dterms;
  CFparams_internal.upsampling    = CFparams->upsampling;
  CFparams_internal.SSBprec       = CFparams->SSBprec;
  CFparams_internal.useRAA        = CFparams->useRAA;
  CFparams_internal.bufferedRAA   = CFparams->bufferedRAA;
  CFparams_internal.edat          = CFparams->edat;
  CFparams_internal.returnAtoms   = TRUE;


  /* get number of segments and detectors for loop scopes */
  if ( multiSFTs->length == 0 ) {
    XLALPrintError ("\nError in function %s, line %d : Input multiSFT vector over segments has zero length!\n\n", fn, __LINE__);
    XLAL_ERROR ( fn, XLAL_EBADLEN);
  }
  numSegments = multiSFTs->length;

  if ( !multiSFTs->data[0] ) {
    XLALPrintError ("\nError in function %s, line %d : Input multiSFT vector has no elements!\n\n", fn, __LINE__);
    XLAL_ERROR ( fn, XLAL_EFAULT);
  }
  if ( multiSFTs->data[0]->length == 0 ) {
    XLALPrintError ("\nError in function %s, line %d : Input multiSFT vector over detectors has zero length!\n\n", fn, __LINE__);
    XLAL_ERROR ( fn, XLAL_EBADLEN);
  }
  numDetectors = multiSFTs->data[0]->length;

  /* initialiase temporary LVcomponents structures */
  TempStats.TwoFX = XLALCreateREAL8Vector ( numDetectors );
  lineVeto->TwoF = 0.0;
  lineVeto->LV   = 0.0;
  for (X = 0; X < numDetectors; X++) {
   lineVeto->TwoFX->data[X] = 0.0;
  }

  /* get SFT duration */
  if ( multiSFTs->data[0]->data[0]->data[0].deltaF == 0 ) {
    XLALPrintError ("\nError in function %s, line %d : deltaF=0, cannot compute Tsft!\n\n", fn, __LINE__);
    XLAL_ERROR ( fn, XLAL_EFAILED);
  }
  Tsft = 1.0 / multiSFTs->data[0]->data[0]->data[0].deltaF; /* define the length of an SFT (assuming 1/Tsft resolution) */
  if ( Tsft == 0 ) {
    XLALPrintError ("\nError in function %s, line %d : Got Tsft=0!\n\n", fn, __LINE__);
    XLAL_ERROR ( fn, XLAL_EFAILED);
  }

  /* compute single-detector Fstats for each data segment, sum up and get line veto stat */
  for (k = 0; k < numSegments; k++)
    { /* loop over segments */

      /* recompute multi-detector Fstat and atoms */
      fakeStatus = blank_status;
      xlalErrno = 0;
      ComputeFStat ( &fakeStatus, &Fstat, dopplerParams, multiSFTs->data[k], multiNoiseWeights->data[k], multiDetStates->data[k], &CFparams_internal, &CFbuffer );
      if ( fakeStatus.statusCode ) {
        XLALPrintError ("\nError in function %s, line %d : Failed call to LAL function ComputeFStat. statusCode=%d\n\n", fn, __LINE__, fakeStatus.statusCode);
        XLAL_ERROR ( fn, XLAL_EFUNC );
      }

      if ( SignalOnly ) {      /* normalization factor correction */
        Fstat.F *= 2.0 / Tsft;
        Fstat.F += 2;
      }

      TempStats.TwoF  = 2.0*Fstat.F;     /* ComputeFStat gives F, we need 2F*/
      lineVeto->TwoF  += TempStats.TwoF; /* sum up multi-detector Fstat for this segment*/

      /* recompute single-detector Fstats from atoms */
      for (X = 0; X < numDetectors; X++)
        { /* loop over detectors */

          xlalErrno = 0;
          TempStats.TwoFX->data[X] = 2.0*XLALComputeFstatFromAtoms ( &Fstat, X );
          if ( xlalErrno != 0 ) {
	    XLALPrintError ("\nError in function %s, line %d : Failed call to XLALComputeFstatFromAtoms.\n\n", fn, __LINE__);
	    XLAL_ERROR ( fn, XLAL_EFUNC );
	  }

          if ( SignalOnly ) {                      /* normalization factor correction */
            TempStats.TwoFX->data[X] *= 2.0 / Tsft;
            TempStats.TwoFX->data[X] += 2;
          }

          lineVeto->TwoFX->data[X]  += TempStats.TwoFX->data[X]; /* sum up single-detector Fstat for this segment*/

        } /* end loop over detectors X */

      /* free memory for atoms that was allocated within ComputeFStat  */
      if (Fstat.multiFstatAtoms) {
	XLALDestroyMultiFstatAtomVector ( Fstat.multiFstatAtoms );
      }

    } /* end loop over segments k */

    /* get average stats over all segments */
    lineVeto->TwoF = lineVeto->TwoF/numSegments;
    for (X = 0; X < numDetectors; X++) {
      lineVeto->TwoFX->data[X] = lineVeto->TwoFX->data[X]/numSegments;
    }

  /* free memory of temporary structures  */
  XLALDestroyREAL8Vector ( TempStats.TwoFX );
  XLALEmptyComputeFBuffer ( &CFbuffer );

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

  Fa.re = 0.0;
  Fa.im = 0.0;
  Fb.re = 0.0;
  Fb.im = 0.0;
  for ( alpha = 0; alpha < numSFTs; alpha++) {
    mmatrixA += Fstat->multiFstatAtoms->data[X]->data[alpha].a2_alpha;
    mmatrixB += Fstat->multiFstatAtoms->data[X]->data[alpha].b2_alpha;
    mmatrixC += Fstat->multiFstatAtoms->data[X]->data[alpha].ab_alpha;
    Fa.re    += Fstat->multiFstatAtoms->data[X]->data[alpha].Fa_alpha.re;
    Fa.im    += Fstat->multiFstatAtoms->data[X]->data[alpha].Fa_alpha.im;
    Fb.re    += Fstat->multiFstatAtoms->data[X]->data[alpha].Fb_alpha.re;
    Fb.im    += Fstat->multiFstatAtoms->data[X]->data[alpha].Fb_alpha.im;
  }
  FX = (mmatrixB*(pow(Fa.re,2)+pow(Fa.im,2)) + mmatrixA*(pow(Fb.re,2)+pow(Fb.im,2)) - 2.0*mmatrixC*(Fa.re*Fb.re + Fa.im*Fb.im) ) / (mmatrixA*mmatrixB-pow(mmatrixC,2));

  return( FX );

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
