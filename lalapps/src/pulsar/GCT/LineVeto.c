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
 */

/*---------- INCLUDES ----------*/
#define __USE_ISOC99 1
#include <lal/StringVector.h>
#include "LineVeto.h"

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
				      const char *listEntryTypeName,                          /**< type of toplist entries, give as name string */
				      const MultiSFTVectorSequence *multiSFTsV,               /**< data files (SFTs) for all detectors and segments */
				      const MultiNoiseWeightsSequence *multiNoiseWeightsV,    /**< noise weights for all detectors and segments */
				      const MultiDetectorStateSeriesSequence *multiDetStatesV,/**< some state info for all detectors */
				      const ComputeFParams *CFparams,                         /**< additional parameters needed for ComputeFStat */
				      const LIGOTimeGPS refTimeGPS,                           /**< reference time for fkdot values in toplist */
				      const BOOLEAN SignalOnly,                               /**< flag for case with no noise, makes some extra checks necessary */
				      const char* outputSingleSegStats                        /**< base filename to output Fstats for each segment individually */
				    )
{
  /* check input parameters and report errors */
  if ( !list || !multiSFTsV || !listEntryTypeName || !multiDetStatesV || !CFparams ) {
    XLALPrintError ("\nError in function %s, line %d : Empty pointer as input parameter!\n\n", __func__, __LINE__);
    XLAL_ERROR ( XLAL_EFAULT);
  }

  if ( !list->data || !list->heap ) {
    XLALPrintError ("\nError in function %s, line %d : Input toplist has no elements!\n\n", __func__, __LINE__);
    XLAL_ERROR ( XLAL_EFAULT);
  }

  if ( list->elems == 0 ) {
    XLALPrintError ("\nError in function %s, line %d : Input toplist has zero length!\n\n", __func__, __LINE__);
    XLAL_ERROR ( XLAL_EBADLEN );
  }

  /* check listEntryTypeName only once by strcmp, afterwards by int, to be faster */
  UINT4 listEntryType = 0;
  if (strcmp(listEntryTypeName, "GCTtop") == 0 )
    listEntryType = 1;
  if (strcmp(listEntryTypeName, "HoughFStat") == 0 )
    listEntryType = 2;
  if ( listEntryType == 0 ) {
    XLALPrintError ("\nError in function %s, line %d : Unsupported entry type for input toplist! Supported types currently are: GCTtop, HoughFStat.\n\n", __func__, __LINE__);
    XLAL_ERROR ( XLAL_EBADLEN );
  }

  if ( !multiSFTsV->data[0] || (multiSFTsV->data[0]->length == 0) ) {
    XLALPrintError ("\nError in function %s, line %d : Input multiSFT vector has no elements!\n\n", __func__, __LINE__);
    XLAL_ERROR ( XLAL_EFAULT);
  }

  /* set up temporary variables and structs */
  PulsarDopplerParams candidateDopplerParams = empty_PulsarDopplerParams; /* struct containing sky position, frequency and fdot for the current candidate */
  UINT4 X;

  /* initialize doppler parameters */
  candidateDopplerParams.refTime = refTimeGPS;  /* spin parameters in toplist refer to this refTime */

  /* initialise detector name vector for later identification */
  LALStringVector *detectorIDs;
  if ( ( detectorIDs = XLALGetDetectorIDs ( multiSFTsV )) == NULL ) /* fill detector name vector with all detectors present in any data sements */
    XLAL_ERROR ( XLAL_EFUNC );

  UINT4 numDetectors = detectorIDs->length;

  /* initialise LVcomponents structure and allocate memory */
  LVcomponents   lineVeto;      /* struct containing multi-detector Fstat, single-detector Fstats, Line Veto stat */
  if ( (lineVeto.TwoFX = XLALCreateREAL8Vector ( numDetectors )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCreateREAL8Vector( %d )\n", __func__, numDetectors );
    XLAL_ERROR ( XLAL_EFUNC );
  }

  UINT4 j;
  UINT4 numElements = list->elems;
  /* loop over toplist: re-compute sumTwoF and sumTwoFX for all candidates */
  for (j = 0; j < numElements; j++ )
    {
      /* set up file for individual-segment Fstat output */
      FILE *singleSegStatsFile = NULL;
      if ( outputSingleSegStats ) {
        char jstring[16];
        snprintf ( jstring, sizeof(jstring), "%d", j );
        char *singleSegStatsFileName = NULL;
        singleSegStatsFileName = LALCalloc( strlen(outputSingleSegStats) + 6 +  strlen(jstring) + 4 + 1, sizeof(CHAR) );
        strcpy(singleSegStatsFileName, outputSingleSegStats);
        strcat(singleSegStatsFileName, "_cand_");
        strcat(singleSegStatsFileName, jstring);
        strcat(singleSegStatsFileName, ".dat");
        /* open the file for writing */
        if ((singleSegStatsFile = fopen(singleSegStatsFileName, "wb")) == NULL) {
          fprintf(stderr, "Unable to open file %s for writing\n", singleSegStatsFileName);
          LALFree(singleSegStatsFile);
          /*exit*/
          XLAL_ERROR ( XLAL_EIO );
        }
        LALFree(singleSegStatsFileName);
      } /* if outputSingleSegStats */

      REAL4Vector *sumTwoFX;
      if ( (sumTwoFX = XLALCreateREAL4Vector ( numDetectors )) == NULL ) {
        XLALPrintError ("%s: failed to XLALCreateREAL4Vector( %d )\n", __func__, numDetectors );
        XLAL_ERROR ( XLAL_EFUNC );
      }

      void *elemV;
      if ( listEntryType == 1 ) {
        GCTtopOutputEntry *elem = toplist_elem ( list, j );
        elemV = elem;

        elem->sumTwoFX = sumTwoFX;
        /* get frequency, sky position, doppler parameters from toplist candidate and save to dopplerParams */
        candidateDopplerParams.Alpha = elem->Alpha;
        candidateDopplerParams.Delta = elem->Delta;
        candidateDopplerParams.fkdot[0] = elem->Freq;
        candidateDopplerParams.fkdot[1] = elem->F1dot;
      } else if ( listEntryType == 2 ) {
        HoughFStatOutputEntry *elem = toplist_elem ( list, j );
        elemV = elem;

        elem->sumTwoFX = sumTwoFX;
        /* get frequency, sky position, doppler parameters from toplist candidate and save to dopplerParams */
        candidateDopplerParams.Alpha = elem->AlphaBest;
        candidateDopplerParams.Delta = elem->DeltaBest;
        candidateDopplerParams.fkdot[0] = elem->Freq;
        candidateDopplerParams.fkdot[1] = elem->f1dot;
      } /* if listEntryType 2 */

      /* write header information into segment-Fstats file */
      if ( singleSegStatsFile )
        fprintf ( singleSegStatsFile, "%%%% Freq: %.16g\n%%%% RA: %.13g\n%%%% Dec: %.13g\n%%%% f1dot: %.13g\n%%%% reftime: %d\n",
                  candidateDopplerParams.fkdot[0], candidateDopplerParams.Alpha, candidateDopplerParams.Delta, candidateDopplerParams.fkdot[1], refTimeGPS.gpsSeconds );

      /*  recalculate multi- and single-IFO Fstats for all segments for this candidate */
      XLALComputeExtraStatsSemiCoherent( &lineVeto, &candidateDopplerParams, multiSFTsV, multiNoiseWeightsV, multiDetStatesV, detectorIDs, CFparams, SignalOnly, singleSegStatsFile );
      if ( xlalErrno != 0 ) {
        XLALPrintError ("\nError in function %s, line %d : Failed call to XLALComputeLineVetoSemiCoherent().\n\n", __func__, __LINE__);
        XLAL_ERROR ( XLAL_EFUNC );
      }

      /* save values in toplist */
      if ( listEntryType == 1 )
        {
          GCTtopOutputEntry *elem = elemV;

          elem->sumTwoFnew         = lineVeto.TwoF;
          for ( X = 0; X < numDetectors; X ++ )
            elem->sumTwoFX->data[X]  = lineVeto.TwoFX->data[X];
        }
      else if ( listEntryType == 2 )
        {
          HoughFStatOutputEntry *elem = elemV;

          elem->sumTwoF         = lineVeto.TwoF;
          for ( X = 0; X < numDetectors; X ++ )
            elem->sumTwoFX->data[X]  = lineVeto.TwoFX->data[X];
        }

      /* close single-segment Fstat file */
      if ( singleSegStatsFile )
        fclose (singleSegStatsFile);

    } /* for j < numElements */

  /* free temporary structures */
  XLALDestroyREAL8Vector ( lineVeto.TwoFX );
  XLALDestroyStringVector ( detectorIDs );

  return (XLAL_SUCCESS);

} /* XLALComputeExtraStatsForToplist() */






/** XLAL Function to recalculate single-IFO Fstats for all semicoherent search segments, and use them to compute Line Veto statistics
*/
int XLALComputeExtraStatsSemiCoherent ( LVcomponents *lineVeto,                                 /**< [out] structure containing multi TwoF, single TwoF, LV stat */
					const PulsarDopplerParams *dopplerParams,               /**< sky position, frequency and fdot for a given candidate */
					const MultiSFTVectorSequence *multiSFTsV,               /**< data files (SFTs) for all detectors and segments */
					const MultiNoiseWeightsSequence *multiNoiseWeightsV,    /**< noise weights for all detectors and segments */
					const MultiDetectorStateSeriesSequence *multiDetStatesV,/**< some state info for all detectors */
					const LALStringVector *detectorIDs,                     /**< name strings of all detectors present in multiSFTsV */
					const ComputeFParams *CFparams,                         /**< additional parameters needed for ComputeFStat */
					const BOOLEAN SignalOnly,                               /**< flag for case with no noise, makes some extra checks necessary */
					FILE *singleSegStatsFile                                /**< pointer to file to output Fstats for each segment individually */
				      )
{
  /* check input parameters and report errors */
  if ( !lineVeto || !lineVeto->TwoFX || !lineVeto->TwoFX->data || !dopplerParams || !multiSFTsV || !multiDetStatesV || !CFparams ) {
    XLALPrintError ("\nError in function %s, line %d : Empty pointer as input parameter!\n\n", __func__, __LINE__);
    XLAL_ERROR ( XLAL_EFAULT);
  }
  if ( multiSFTsV->length == 0 ) {
    XLALPrintError ("\nError in function %s, line %d : Input multiSFT vector over segments has zero length!\n\n", __func__, __LINE__);
    XLAL_ERROR ( XLAL_EBADLEN);
  }

  /* fake LAL status structure, needed as long as ComputeFStat is LAL function and not XLAL */
  LALStatus fakeStatus = blank_status;
  Fcomponents    Fstat;                           /* temporary struct for ComputeFStat */
  UINT4 X, Y;

  UINT4 numSegments  = multiSFTsV->length;
  UINT4 numDetectors = detectorIDs->length;

  if ( lineVeto->TwoFX->length != numDetectors ) {
    XLALPrintError ("\%s, line %d : Inconsistent number of detectors: TwoFX vector has length %d, while detectorID list contains %d elements!\n\n", __func__, __LINE__, lineVeto->TwoFX->length, numDetectors );
    XLAL_ERROR ( XLAL_EBADLEN );
  }

  /* temporary copy of Fstatistic parameters structure, needed to change returnAtoms for function scope only */
  ComputeFParams CFparams_internal = (*CFparams);
  CFparams_internal.returnAtoms   = TRUE;

  /* initialiase LVcomponents structure */
  lineVeto->TwoF = 0.0;
  lineVeto->LV   = 0.0;
  for (X = 0; X < numDetectors; X++) {
    lineVeto->TwoFX->data[X] = 0.0;
  }

  REAL8 Tsft = 1.0 / multiSFTsV->data[0]->data[0]->data[0].deltaF;	/* get SFT duration */

  /* variables necessary to catch segments where not all detectors have data */
  INT4 detid = -1; /* count through detector IDs for matching with name strings */
  UINT4 numSegmentsX[numDetectors];  /* number of segments with data might be different for each detector */
  for (X = 0; X < numDetectors; X++)
    numSegmentsX[X] = 0;


  REAL8Vector *twoFXseg = NULL;
  if ( (twoFXseg = XLALCreateREAL8Vector ( numDetectors )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCreateREAL8Vector( %d )\n", __func__, numDetectors );
    XLAL_ERROR ( XLAL_EFUNC );
  }

  /* internal dopplerParams structure, for extrapolating to correct reftimes for each segment */
  PulsarDopplerParams dopplerParams_temp = empty_PulsarDopplerParams; /* struct containing sky position, frequency and fdot for the current candidate */
  dopplerParams_temp.Alpha = dopplerParams->Alpha;
  dopplerParams_temp.Delta = dopplerParams->Delta;
  INIT_MEM( dopplerParams_temp.fkdot );

  /* compute single- and multi-detector Fstats for each data segment and sum up */
  UINT4 k;
  for (k = 0; k < numSegments; k++)
    {
      UINT4 numDetectorsSeg = multiSFTsV->data[k]->length; /* for each segment, number of detectors with data might be smaller than overall number */

      /* initialize temporary single-IFO Fstat vector */
      for (X = 0; X < numDetectors; X++)
        twoFXseg->data[X] = 0.0;

      /* starttime of segment: could be different for each detector, take minimum */
      dopplerParams_temp.refTime.gpsSeconds = multiSFTsV->data[k]->data[0]->data[0].epoch.gpsSeconds;
      for (X = 0; X < numDetectorsSeg; X++) {
        if ( multiSFTsV->data[k]->data[X]->data[0].epoch.gpsSeconds < dopplerParams_temp.refTime.gpsSeconds )
          dopplerParams_temp.refTime = multiSFTsV->data[k]->data[X]->data[0].epoch;
      }
      REAL8 deltaTau = XLALGPSDiff( &dopplerParams_temp.refTime, &dopplerParams->refTime ); /* convert LIGOTimeGPS into real number difference for XLALExtrapolatePulsarSpins */
      /* extrapolate pulsar spins to correct time for this segment */
      if ( XLALExtrapolatePulsarSpins( dopplerParams_temp.fkdot, dopplerParams->fkdot, deltaTau ) != XLAL_SUCCESS ) {
        XLALPrintError ("\n%s, line %d : XLALExtrapolatePulsarSpins() failed.\n\n", __func__, __LINE__);
        XLAL_ERROR ( XLAL_EFUNC );
      }

      MultiNoiseWeights *multiNoiseWeightsThisSeg;
      if ( multiNoiseWeightsV )
        multiNoiseWeightsThisSeg = multiNoiseWeightsV->data[k];
      else
        multiNoiseWeightsThisSeg = NULL;

      /* recompute multi-detector Fstat and atoms */
      if ( singleSegStatsFile )
        fprintf ( singleSegStatsFile, "%%%% Reftime: %d %%%% Freq: %.16g %%%% RA: %.13g %%%% Dec: %.13g %%%% f1dot: %.13g\n", dopplerParams_temp.refTime.gpsSeconds, dopplerParams_temp.fkdot[0], dopplerParams_temp.Alpha, dopplerParams_temp.Delta, dopplerParams_temp.fkdot[1] );
      fakeStatus = blank_status;
      ComputeFStat ( &fakeStatus, &Fstat, &dopplerParams_temp, multiSFTsV->data[k], multiNoiseWeightsThisSeg, multiDetStatesV->data[k], &CFparams_internal, NULL );
      if ( fakeStatus.statusCode ) {
        XLALPrintError ("\%s, line %d : Failed call to LAL function ComputeFStat(). statusCode=%d\n\n", __func__, __LINE__, fakeStatus.statusCode);
        XLAL_ERROR ( XLAL_EFUNC );
      }

      if ( SignalOnly ) {      /* normalization factor correction */
        Fstat.F *= 2.0 / Tsft;
        Fstat.F += 2;
      }

      lineVeto->TwoF  += 2.0 * Fstat.F; /* sum up multi-detector Fstat for this segment*/

      if ( singleSegStatsFile )
        fprintf ( singleSegStatsFile, "%.6f", 2.0*Fstat.F );

      /* recompute single-detector Fstats from atoms */
      for (X = 0; X < numDetectorsSeg; X++)
        {

          /* match detector ID in this segment to one from detectorIDs list, sum up the corresponding twoFX */
          detid = -1;
          for (Y = 0; Y < numDetectors; Y++)
          {
            if ( strcmp( multiSFTsV->data[k]->data[X]->data[0].name, detectorIDs->data[Y] ) == 0 )
              detid = Y;
          }
          if ( detid == -1 )
          {
            XLALPrintError ("\nError in function %s, line %d : For segment k=%d, detector X=%d, could not match detector ID %s.\n\n", __func__, __LINE__, k, X, multiSFTsV->data[k]->data[X]->data[0].name);
            XLAL_ERROR ( XLAL_EFAILED );
          }
          numSegmentsX[detid] += 1; /* have to keep this for correct averaging */

          twoFXseg->data[detid] = 2.0 * XLALComputeFstatFromAtoms ( Fstat.multiFstatAtoms, X );
          if ( xlalErrno != 0 ) {
            XLALPrintError ("\nError in function %s, line %d : Failed call to XLALComputeFstatFromAtoms().\n\n", __func__, __LINE__);
            XLAL_ERROR ( XLAL_EFUNC );
          }

          if ( SignalOnly ) {                      /* normalization factor correction (TwoF=2.0*F has been done before, this time!) */
            twoFXseg->data[detid] *= 4.0 / Tsft;
            twoFXseg->data[detid] += 4;
          }

          lineVeto->TwoFX->data[detid]  += twoFXseg->data[detid]; /* sum up single-detector Fstat for this segment*/

        } /* for X < numDetectorsSeg */

      if ( singleSegStatsFile ) {
        for (X = 0; X < numDetectors; X++)
          fprintf ( singleSegStatsFile, " %.6f",twoFXseg->data[X] );
        fprintf ( singleSegStatsFile, "\n" );
      }

      /* free memory for atoms that was allocated within ComputeFStat  */
      XLALDestroyMultiFstatAtomVector ( Fstat.multiFstatAtoms );

    } /* for k < numSegments */

  /* get average stats over all segments */
  lineVeto->TwoF /= numSegments;
  for (X = 0; X < numDetectors; X++) {
    lineVeto->TwoFX->data[X] /= numSegmentsX[X];
  }

  XLALDestroyREAL8Vector(twoFXseg);

  return(XLAL_SUCCESS);

} /* XLALComputeExtraStatsSemiCoherent() */


/** XLAL function to compute single-or multi-IFO Fstat from multi-IFO Atoms: */
REAL8 XLALComputeFstatFromAtoms ( const MultiFstatAtomVector *multiFstatAtoms,   /**< multi-detector atoms */
				  const INT4                 X                   /**< detector number, give -1 for multi-Fstat */
				  )
{
  /* check input parameters and report errors */
  if ( !multiFstatAtoms || !multiFstatAtoms->data || !multiFstatAtoms->data[0]->data ) {
    XLALPrintError ("\nError in function %s, line %d : Empty pointer as input parameter!\n\n", __func__, __LINE__);
    XLAL_ERROR ( XLAL_EFAULT);
  }

  if ( multiFstatAtoms->length == 0 ) {
    XLALPrintError ("\nError in function %s, line %d : Input MultiFstatAtomVector has zero length! (no detectors)\n\n", __func__, __LINE__);
    XLAL_ERROR ( XLAL_EBADLEN );
  }

  if ( X < -1 ) {
    XLALPrintError ("\nError in function %s, line %d : Invalid detector number X=%d, only nonnegative numbers or -1 for multi-F are allowed!\n\n", __func__, __LINE__, X);
    XLAL_ERROR ( XLAL_EDOM );
  }

  if ( ( X >= 0 ) && ( (UINT4)(X) > multiFstatAtoms->length-1 ) ) {
    XLALPrintError ("\nError in function %s, line %d : Invalid detector number!\nRequested X=%d, but FstatAtoms only have length %d.\n\n", __func__, __LINE__, X, multiFstatAtoms->length);
    XLAL_ERROR ( XLAL_EDOM );
  }

  /* internal detector index Y to do both single- and multi-F case */
  UINT4 Y, Ystart, Yend;
  if ( X == -1 ) { /* loop through all detectors to get multi-Fstat */
    Ystart = 0;
    Yend   = multiFstatAtoms->length-1;
  }
  else { /* just compute single-Fstat for 1 IFO */
    Ystart = X;
    Yend   = X;
  }

  /* set up temporary Fatoms and matrix elements for summations */
  REAL8 mmatrixA = 0.0, mmatrixB = 0.0, mmatrixC = 0.0;
  REAL8 F = 0.0;
  COMPLEX8 Fa, Fb;
  Fa.re = 0.0;
  Fa.im = 0.0;
  Fb.re = 0.0;
  Fb.im = 0.0;

  for (Y = Ystart; Y <= Yend; Y++) {  /* loop through detectors */

    UINT4 alpha, numSFTs;
    numSFTs = multiFstatAtoms->data[Y]->length;
    if ( numSFTs == 0 ) {
      XLALPrintError ("\nError in function %s, line %d : Input FstatAtomVector has zero length! (no timestamps for detector X=%d)\n\n", __func__, __LINE__, Y);
      XLAL_ERROR ( XLAL_EDOM );
    }

    for ( alpha = 0; alpha < numSFTs; alpha++) { /* loop through SFTs */
      FstatAtom *thisAtom = &multiFstatAtoms->data[Y]->data[alpha];
      /* sum up matrix elements and Fa, Fb */
      mmatrixA += thisAtom->a2_alpha;
      mmatrixB += thisAtom->b2_alpha;
      mmatrixC += thisAtom->ab_alpha;
      Fa.re    += thisAtom->Fa_alpha.re;
      Fa.im    += thisAtom->Fa_alpha.im;
      Fb.re    += thisAtom->Fb_alpha.re;
      Fb.im    += thisAtom->Fb_alpha.im;
    } /* loop through SFTs */

  } /* loop through detectors */

  /* compute determinant and final Fstat (not twoF!) */
  REAL8 Dinv = 1.0 / ( mmatrixA * mmatrixB - SQUARE(mmatrixC) );
  F = Dinv * ( mmatrixB * ( SQUARE(Fa.re) + SQUARE(Fa.im) ) + mmatrixA * ( SQUARE(Fb.re) + SQUARE(Fb.im) ) - 2.0 * mmatrixC * (Fa.re*Fb.re + Fa.im*Fb.im) );

  return(F);

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
  /* check input parameters and report errors */
  if ( !TwoF || !TwoFX || !TwoFX->data || !priorX || !priorX->data ) {
    XLALPrintError ("\nError in function %s, line %d : Empty pointer as input parameter!\n\n", __func__, __LINE__);
    XLAL_ERROR_REAL8 ( XLAL_EFAULT);
  }

  if ( TwoFX->length == 0 ) {
    XLALPrintError ("\nError in function %s, line %d :Input TwoFX vector has zero length!\n\n", __func__, __LINE__);
    XLAL_ERROR_REAL8 ( XLAL_EBADLEN);
  }

  if ( priorX->length == 0 ) {
    XLALPrintError ("\nError in function %s, line %d :Input priorX vector has zero length!\n\n", __func__, __LINE__);
    XLAL_ERROR_REAL8 ( XLAL_EBADLEN);
  }

  /* set up temporary variables and structs */
  UINT4 X;                            /* loop summation variable */
  UINT4 numDetectors = TwoFX->length; /* loop summation upper limit for detectors */
  REAL8 maxSum = -1000.0;             /* maximum of terms in denominator, for logsumexp formula */
  REAL8 LV = 0.0;                     /* output variable for Line Veto statistics */

  if ( rhomax < 0.0 ) {
    XLALPrintError ("\nError in function %s, line %d : nonpositive input rhomax!\n\n", __func__, __LINE__);
    XLAL_ERROR_REAL8 ( XLAL_EFPINVAL);
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
    XLALPrintError ("\nError in function %s, line %d : log(nonpositive) in LV denominator. \n\n", __func__, __LINE__);
    XLAL_ERROR_REAL8 ( XLAL_EFPINVAL );
  }
  else {
    LV = TwoF - maxSum - log( LV );
    return(LV);
  }


} /* XLALComputeLineVeto() */




/** XLAL function to get a list of detector IDs from multi-segment multiSFT vectors
* returns all unique detector IDs for cases with some detectors switching on and off
*/
LALStringVector *
XLALGetDetectorIDs ( const MultiSFTVectorSequence *multiSFTsV /**< data files (SFTs) for all detectors and segments */
                     )
{
  LALStringVector *IFOList = NULL;	// IFO string vector for returning

  /* check input parameters and report errors */
  if ( !multiSFTsV || !multiSFTsV->data ) {
    XLALPrintError ("\nError in function %s, line %d : Empty pointer as input parameter!\n\n", __func__, __LINE__);
    XLAL_ERROR_NULL ( XLAL_EFAULT);
  }

  if ( multiSFTsV->length == 0 ) {
    XLALPrintError ("\nError in function %s, line %d : Input multiSFT vector contains no segments (zero length)!\n\n", __func__, __LINE__);
    XLAL_ERROR_NULL ( XLAL_EBADLEN);
  }

  /* set up variables */
  UINT4 k, X;
  UINT4 numSegments  = multiSFTsV->length; /* number of segments with SFTs from any detector */
  UINT4 numDetectorsSeg = 0;               /* number of detectors with data might be different for each segment */

  for (k = 0; k < numSegments; k++)
  {
    if ( !multiSFTsV->data[k] || (multiSFTsV->data[k]->length == 0) ) {
    XLALPrintError ("\nError in function %s, line %d : Input multiSFT vector, segment k=%d has no elements (0 detectors)!\n\n", __func__, __LINE__, k);
    XLAL_ERROR_NULL ( XLAL_EFAULT);
    }

    numDetectorsSeg = multiSFTsV->data[k]->length; /* for each segment, could be smaller than overall number */

    for (X = 0; X < numDetectorsSeg; X++)
    {
      /* check if this segment k, IFO X contains a new detector */
      const char *thisIFO = multiSFTsV->data[k]->data[X]->data[0].name;

      if ( XLALFindStringInVector ( thisIFO, IFOList ) >= 0 )
        { // if already in list, do nothing
          continue;
        }
      else
        {       // otherwise, append to IFOList
          if ( (IFOList = XLALAppendString2Vector ( IFOList, thisIFO )) == NULL )
            XLAL_ERROR_NULL ( XLAL_EFUNC );
        }

    } /* for X < numDetectorsSeg */

  } /* for k < numSegments */

  return IFOList;

} /* XLALGetDetectorIDs() */
