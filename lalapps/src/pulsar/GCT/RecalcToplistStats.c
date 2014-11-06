/*
 *  Copyright (C) 2013 Karl Wette
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
#include "RecalcToplistStats.h"

/*---------- local DEFINES ----------*/
#define TRUE (1==1)
#define FALSE (1==0)

/*----- Macros ----- */

/*----- SWITCHES -----*/

/*---------- internal types ----------*/

/*---------- Global variables ----------*/

/*---------- internal prototypes ----------*/

/*==================== FUNCTION DEFINITIONS ====================*/


/** XLAL function to go through a (Hough or GCT) toplist and compute line-robust statistics for each candidate */
int XLALComputeExtraStatsForToplist ( toplist_t *list,				/**< list of cancidates with f, sky position etc. - no output so far */
				      const RecalcStatsParams *recalcParams	/**< additional input values and parameters */
				    )
{

  /* check input parameters and report errors */
  XLAL_CHECK ( list && recalcParams->listEntryTypeName && recalcParams->Fstat_in_vec && recalcParams->detectorIDs, XLAL_EFAULT, "Empty pointer as input parameter." );
  XLAL_CHECK ( list->data && list->heap, XLAL_EFAULT, "Input toplist has no elements." );
  XLAL_CHECK ( list->elems > 0, XLAL_EBADLEN, "Input toplist has zero length." );

  /* check listEntryTypeName only once by strcmp, afterwards by int, to be faster */
  UINT4 listEntryType = 0;
  if (strcmp(recalcParams->listEntryTypeName, "GCTtop") == 0 ) {
    listEntryType = 1;
  }
  if (strcmp(recalcParams->listEntryTypeName, "HoughFStat") == 0 ) {
    listEntryType = 2;
  }
  XLAL_CHECK ( listEntryType != 0, XLAL_EBADLEN, "Unsupported entry type for input toplist! Supported types currently are: GCTtop, HoughFStat." );

  /* set up temporary variables and structs */
  PulsarDopplerParams XLAL_INIT_DECL(candidateDopplerParams); /* struct containing sky position, frequency and fdot for the current candidate */
  UINT4 X;

  /* initialize doppler parameters */
  candidateDopplerParams.refTime = recalcParams->refTimeGPS;  /* spin parameters in toplist refer to this refTime */

  UINT4 numDetectors = recalcParams->detectorIDs->length;

  UINT4 j;
  UINT4 numElements = list->elems;
  /* loop over toplist: re-compute TwoF and TwoFX for all candidates (average over segments) */
  for (j = 0; j < numElements; j++ )
    {

      void *elemV = NULL;
      if ( listEntryType == 1 ) {
        GCTtopOutputEntry *elem = toplist_elem ( list, j );
        elemV = elem;

        /* get frequency, sky position, doppler parameters from toplist candidate and save to dopplerParams */
        candidateDopplerParams.Alpha = elem->Alpha;
        candidateDopplerParams.Delta = elem->Delta;
        candidateDopplerParams.fkdot[0] = elem->Freq;
        candidateDopplerParams.fkdot[1] = elem->F1dot;
        candidateDopplerParams.fkdot[2] = elem->F2dot;
      }
      else if ( listEntryType == 2 ) {
        HoughFStatOutputEntry *elem = toplist_elem ( list, j );
        elemV = elem;

        XLAL_CHECK ( (elem->sumTwoFX = XLALCreateREAL4Vector ( numDetectors )) != NULL, XLAL_EFUNC, "Failed call to XLALCreateREAL4Vector( %d ).", numDetectors );

        /* get frequency, sky position, doppler parameters from toplist candidate and save to dopplerParams */
        candidateDopplerParams.Alpha = elem->AlphaBest;
        candidateDopplerParams.Delta = elem->DeltaBest;
        candidateDopplerParams.fkdot[0] = elem->Freq;
        candidateDopplerParams.fkdot[1] = elem->f1dot;
        /* no 2nd spindown in HoughFStatOutputEntry */
      } /* if listEntryType 2 */

      /*  recalculate multi- and single-IFO Fstats for all segments for this candidate */
      RecalcStatsComponents XLAL_INIT_DECL(recalcStats); /* struct containing multi-detector F-stat, single-detector F-stats, BSGL */
      recalcStats.log10BSGL = -LAL_REAL4_MAX; /* proper initialization here is not 0 */
      XLAL_CHECK ( XLALComputeExtraStatsSemiCoherent( &recalcStats, &candidateDopplerParams, recalcParams ) == XLAL_SUCCESS, XLAL_EFUNC, "Failed call to XLALComputeExtraStatsSemiCoherent()." );

      /* save values in toplist */
      if ( listEntryType == 1 ) {
          GCTtopOutputEntry *elem = elemV;
          elem->numDetectors = numDetectors;
          elem->avTwoFrecalc = recalcStats.avTwoF; /* average over segments */
          for ( X = 0; X < numDetectors; X ++ ) {
            elem->avTwoFXrecalc[X] = recalcStats.avTwoFX[X];
          }
          elem->log10BSGLrecalc = recalcStats.log10BSGL;
          if ( recalcParams->loudestSegOutput ) {
            elem->loudestSeg      = recalcStats.loudestSeg;
            elem->twoFloudestSeg  = recalcStats.twoFloudestSeg;
            for ( X = 0; X < numDetectors; X ++ ) {
              elem->twoFXloudestSeg[X] = recalcStats.twoFXloudestSeg[X];
            }
          }
      } /* if ( listEntryType == 1 ) */
      else if ( listEntryType == 2 ) {
          HoughFStatOutputEntry *elem = elemV;

          elem->sumTwoF = recalcStats.avTwoF; /* this is also the average over segments, the field is only called "sumTwoF" due to Hough legacy */
          for ( X = 0; X < numDetectors; X ++ ) {
            elem->sumTwoFX->data[X]  = recalcStats.avTwoFX[X];
          }
      } /* if ( listEntryType == 2 ) */

    } /* for j < numElements */

  return (XLAL_SUCCESS);

} /* XLALComputeExtraStatsForToplist() */


/**
 * XLAL Function to recalculate multi-IFO F-stat 2F and single-IFO 2FX for all semicoherent search segments
 * This returns AVERAGE F-stats over segments, not sums.
 */
int XLALComputeExtraStatsSemiCoherent ( RecalcStatsComponents *recalcStats,		/**< [out] structure containing multi TwoF, single TwoF, BSGL */
					const PulsarDopplerParams *dopplerParams,	/**< sky position, frequency and fdot for a given candidate */
					const RecalcStatsParams *recalcParams		/**< additional input values and parameters */
				      )
{

  /* check input parameters and report errors */
  XLAL_CHECK ( recalcStats && recalcStats->avTwoFX && dopplerParams && recalcParams->Fstat_in_vec && recalcParams->detectorIDs && recalcParams->startTstack, XLAL_EFAULT, "Empty pointer as input parameter." );

  UINT4 numSegments  = recalcParams->Fstat_in_vec->length;
  UINT4 numDetectors = recalcParams->detectorIDs->length;

  recalcStats->numDetectors = numDetectors;

  /* variables necessary to catch segments where not all detectors have data */
  INT4 detid = -1; /* count through detector IDs for matching with name strings */
  UINT4 numSegmentsX[numDetectors];  /* number of segments with data might be different for each detector */
  for (UINT4 X = 0; X < numDetectors; X++) {
    numSegmentsX[X] = 0;
  }

  /* internal dopplerParams structure, for extrapolating to correct reftimes for each segment */
  PulsarDopplerParams XLAL_INIT_DECL(dopplerParams_temp); /* struct containing sky position, frequency and fdot for the current candidate */
  dopplerParams_temp.Alpha = dopplerParams->Alpha;
  dopplerParams_temp.Delta = dopplerParams->Delta;
  XLAL_INIT_MEM( dopplerParams_temp.fkdot );

  /* just in case the caller hasn't properly initialized recalcStats, make sure everything is 0 before the loop */
  REAL4 sumTwoF = 0.0;
  REAL4 sumTwoFX[numDetectors];
  for (UINT4 X = 0; X < numDetectors; X++) {
    sumTwoFX[X] = 0.0;
    recalcStats->twoFXloudestSeg[X] = 0.0;
  }
  recalcStats->twoFloudestSeg = 0.0;

  /* compute single- and multi-detector Fstats for each data segment and sum up */
  UINT4 k;
  FstatResults* Fstat_res = NULL;
  for (k = 0; k < numSegments; k++)
    {

      /* starttime of segment */
      dopplerParams_temp.refTime = recalcParams->startTstack->data[k];
      /* convert LIGOTimeGPS into real number difference for XLALExtrapolatePulsarSpins */
      REAL8 deltaTau = XLALGPSDiff( &dopplerParams_temp.refTime, &dopplerParams->refTime );
      /* extrapolate pulsar spins to correct time for this segment */
      XLAL_CHECK ( XLALExtrapolatePulsarSpins( dopplerParams_temp.fkdot, dopplerParams->fkdot, deltaTau ) == XLAL_SUCCESS, XLAL_EFUNC, "XLALExtrapolatePulsarSpins() failed." );

      /* recompute multi-detector Fstat and atoms */
      XLAL_CHECK ( XLALComputeFstat(&Fstat_res, recalcParams->Fstat_in_vec->data[k], &dopplerParams_temp, 0.0, 1, FSTATQ_2F | FSTATQ_2F_PER_DET) == XLAL_SUCCESS, XLAL_EFUNC, "XLALComputeFstat() failed with errno=%d", xlalErrno );

      sumTwoF  += Fstat_res->twoF[0]; /* sum up multi-detector Fstat for this segment*/

      BOOLEAN update_loudest = FALSE;
      if ( recalcParams->loudestSegOutput && ( Fstat_res->twoF[0] > recalcStats->twoFloudestSeg ) )
        {
          update_loudest = TRUE;
          recalcStats->loudestSeg = k;
          recalcStats->twoFloudestSeg = Fstat_res->twoF[0];
        }

      /* for each segment, number of detectors with data might be smaller than overall number */
      const UINT4 numDetectorsSeg = Fstat_res->numDetectors;

      /* recompute single-detector Fstats from atoms */
      for (UINT4 X = 0; X < numDetectorsSeg; X++)
        {

          /* match detector ID in this segment to one from detectorIDs list, sum up the corresponding twoFX */
          detid = -1;
          for (UINT4 Y = 0; Y < numDetectors; Y++) {
            if ( strcmp( Fstat_res->detectorNames[X], recalcParams->detectorIDs->data[Y] ) == 0 ) {
              detid = Y;
            }
          }
          if ( detid == -1 ) {
            XLALPrintError ("\nError in function %s, line %d : For segment k=%d, detector X=%d, could not match detector ID %s.\n\n",
                            __func__, __LINE__, k, X, Fstat_res->detectorNames[X] );
            XLAL_ERROR ( XLAL_EFAILED );
          }
          numSegmentsX[detid] += 1; /* have to keep this for correct averaging */

          sumTwoFX[detid] += Fstat_res->twoFPerDet[X][0]; /* sum up single-detector Fstat for this segment*/

          if ( update_loudest ) {
            recalcStats->twoFXloudestSeg[X] = Fstat_res->twoFPerDet[X][0];
          }

        } /* for X < numDetectorsSeg */

    } /* for k < numSegments */

  if ( recalcParams->BSGLsetup )
    {
      recalcStats->log10BSGL = XLALComputeBSGL ( sumTwoF, sumTwoFX, recalcParams->BSGLsetup );
      XLAL_CHECK ( xlalErrno == 0, XLAL_EFUNC, "XLALComputeBSGL() failed with xlalErrno = %d\n", xlalErrno );
    }

  /* get average stats over all segments */
  recalcStats->avTwoF = sumTwoF/numSegments;
  for (UINT4 X = 0; X < numDetectors; X++) {
    recalcStats->avTwoFX[X] = sumTwoFX[X]/numSegmentsX[X];
  }

  XLALDestroyFstatResults(Fstat_res);

  return(XLAL_SUCCESS);

} /* XLALComputeExtraStatsSemiCoherent() */
