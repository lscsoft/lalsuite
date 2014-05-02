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
#define INIT_MEM(x) memset(&(x), 0, sizeof((x)))

/*----- SWITCHES -----*/

/*---------- internal types ----------*/

/*---------- Global variables ----------*/

/*---------- internal prototypes ----------*/

/*==================== FUNCTION DEFINITIONS ====================*/


/** XLAL function to go through a (Hough or GCT) toplist and compute line-robust statistics for each candidate */
int XLALComputeExtraStatsForToplist ( toplist_t *list,                                        /**< list of cancidates with f, sky position etc. - no output so far */
				      const char *listEntryTypeName,                          /**< type of toplist entries, give as name string */
				      const FstatInputVector *Fstat_in_vec,               /**< vector of input data for XLALComputeFstat() */
				      const LALStringVector *detectorIDs,                     /**< detector name vector with all detectors present in any data sements */
				      const LIGOTimeGPSVector *startTstack,                   /**< starting GPS time of each stack */
				      const LIGOTimeGPS refTimeGPS,                           /**< reference time for fkdot values in toplist */
				      const char* outputSingleSegStats                        /**< base filename to output Fstats for each segment individually */
				    )
{

  /* check input parameters and report errors */
  XLAL_CHECK ( list && listEntryTypeName && Fstat_in_vec && detectorIDs, XLAL_EFAULT, "Empty pointer as input parameter." );
  XLAL_CHECK ( list->data && list->heap, XLAL_EFAULT, "Input toplist has no elements." );
  XLAL_CHECK ( list->elems > 0, XLAL_EBADLEN, "Input toplist has zero length." );

  /* check listEntryTypeName only once by strcmp, afterwards by int, to be faster */
  UINT4 listEntryType = 0;
  if (strcmp(listEntryTypeName, "GCTtop") == 0 ) {
    listEntryType = 1;
  }
  if (strcmp(listEntryTypeName, "HoughFStat") == 0 ) {
    listEntryType = 2;
  }
  XLAL_CHECK ( listEntryType != 0, XLAL_EBADLEN, "Unsupported entry type for input toplist! Supported types currently are: GCTtop, HoughFStat." );

  /* set up temporary variables and structs */
  PulsarDopplerParams candidateDopplerParams = empty_PulsarDopplerParams; /* struct containing sky position, frequency and fdot for the current candidate */
  UINT4 X;

  /* initialize doppler parameters */
  candidateDopplerParams.refTime = refTimeGPS;  /* spin parameters in toplist refer to this refTime */

  UINT4 numDetectors = detectorIDs->length;

  /* initialise LVcomponents structure and allocate memory */
  LVcomponents   lineVeto = empty_LVcomponents; /* struct containing multi-detector Fstat, single-detector Fstats, Line Veto stat */
  XLAL_CHECK ( (lineVeto.TwoFX = XLALCreateREAL4Vector ( numDetectors )) != NULL, XLAL_EFUNC, "Failed call to XLALCreateREAL4Vector( %d ).", numDetectors );

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

      /* write header information into segment-Fstats file */
      if ( singleSegStatsFile ) {
        fprintf ( singleSegStatsFile, "%%%% Freq: %.16g\n%%%% RA: %.13g\n%%%% Dec: %.13g\n%%%% f1dot: %.13g\n%%%% f2dot: %.13g\n%%%% reftime: %d\n",
                  candidateDopplerParams.fkdot[0], candidateDopplerParams.Alpha, candidateDopplerParams.Delta, candidateDopplerParams.fkdot[1], candidateDopplerParams.fkdot[2], refTimeGPS.gpsSeconds );
      }

      /*  recalculate multi- and single-IFO Fstats for all segments for this candidate */
      XLAL_CHECK ( XLALComputeExtraStatsSemiCoherent( &lineVeto, &candidateDopplerParams, Fstat_in_vec, detectorIDs, startTstack, singleSegStatsFile ) == XLAL_SUCCESS, XLAL_EFUNC, "Failed call to XLALComputeLineVetoSemiCoherent()." );

      /* save values in toplist */
      if ( listEntryType == 1 ) {
          GCTtopOutputEntry *elem = elemV;
          elem->numDetectors = numDetectors;
          elem->sumTwoFrecalc  = lineVeto.TwoF;
          for ( X = 0; X < numDetectors; X ++ ) {
            elem->sumTwoFXrecalc[X]  = lineVeto.TwoFX->data[X];
          }
      }
      else if ( listEntryType == 2 ) {
          HoughFStatOutputEntry *elem = elemV;

          elem->sumTwoF         = lineVeto.TwoF;
          for ( X = 0; X < numDetectors; X ++ ) {
            elem->sumTwoFX->data[X]  = lineVeto.TwoFX->data[X];
          }
      }

      /* close single-segment Fstat file */
      if ( singleSegStatsFile ) {
        fclose (singleSegStatsFile);
      }

    } /* for j < numElements */

  /* free temporary structures */
  XLALDestroyREAL4Vector ( lineVeto.TwoFX );

  return (XLAL_SUCCESS);

} /* XLALComputeExtraStatsForToplist() */


/**
 * XLAL Function to recalculate single-IFO Fstats for all semicoherent search segments, and use them to compute line-robust statistics
 */
int XLALComputeExtraStatsSemiCoherent ( LVcomponents *lineVeto,                                 /**< [out] structure containing multi TwoF, single TwoF, LV stat */
					const PulsarDopplerParams *dopplerParams,               /**< sky position, frequency and fdot for a given candidate */
					const FstatInputVector *Fstat_in_vec,               /**< vector of input data for XLALComputeFstat() */
					const LALStringVector *detectorIDs,                     /**< detector name vector with all detectors present in any data sements */
					const LIGOTimeGPSVector *startTstack,                   /**< starting GPS time of each stack */
					FILE *singleSegStatsFile                                /**< pointer to file to output Fstats for each segment individually */
				      )
{

  /* check input parameters and report errors */
  XLAL_CHECK ( lineVeto && lineVeto->TwoFX && lineVeto->TwoFX->data && dopplerParams && Fstat_in_vec && detectorIDs && startTstack, XLAL_EFAULT, "Empty pointer as input parameter." );

  UINT4 numSegments  = Fstat_in_vec->length;
  UINT4 numDetectors = detectorIDs->length;

  XLAL_CHECK ( lineVeto->TwoFX->length == numDetectors, XLAL_EBADLEN, "Inconsistent number of detectors: TwoFX vector has length %d, while detectorID list contains %d elements.", lineVeto->TwoFX->length, numDetectors );

  /* initialiase LVcomponents structure */
  lineVeto->TwoF = 0.0;
  lineVeto->LV   = 0.0;
  for (UINT4 X = 0; X < numDetectors; X++) {
    lineVeto->TwoFX->data[X] = 0.0;
  }

  /* variables necessary to catch segments where not all detectors have data */
  INT4 detid = -1; /* count through detector IDs for matching with name strings */
  UINT4 numSegmentsX[numDetectors];  /* number of segments with data might be different for each detector */
  for (UINT4 X = 0; X < numDetectors; X++) {
    numSegmentsX[X] = 0;
  }

  REAL4Vector *twoFXseg = NULL;
  XLAL_CHECK ( ( twoFXseg = XLALCreateREAL4Vector ( numDetectors )) != NULL, XLAL_EFUNC, "Failed call to XLALCreateREAL4Vector( %d ).", numDetectors );

  /* internal dopplerParams structure, for extrapolating to correct reftimes for each segment */
  PulsarDopplerParams dopplerParams_temp = empty_PulsarDopplerParams; /* struct containing sky position, frequency and fdot for the current candidate */
  dopplerParams_temp.Alpha = dopplerParams->Alpha;
  dopplerParams_temp.Delta = dopplerParams->Delta;
  INIT_MEM( dopplerParams_temp.fkdot );

  /* compute single- and multi-detector Fstats for each data segment and sum up */
  UINT4 k;
  FstatResults* Fstat_res = NULL;
  for (k = 0; k < numSegments; k++)
    {

      /* initialize temporary single-IFO Fstat vector */
      for (UINT4 X = 0; X < numDetectors; X++) {
        twoFXseg->data[X] = 0.0;
      }

      /* starttime of segment */
      dopplerParams_temp.refTime = startTstack->data[k];
      /* convert LIGOTimeGPS into real number difference for XLALExtrapolatePulsarSpins */
      REAL8 deltaTau = XLALGPSDiff( &dopplerParams_temp.refTime, &dopplerParams->refTime );
      /* extrapolate pulsar spins to correct time for this segment */
      XLAL_CHECK ( XLALExtrapolatePulsarSpins( dopplerParams_temp.fkdot, dopplerParams->fkdot, deltaTau ) == XLAL_SUCCESS, XLAL_EFUNC, "XLALExtrapolatePulsarSpins() failed." );

      /* recompute multi-detector Fstat and atoms */
      if ( singleSegStatsFile ) {
        fprintf ( singleSegStatsFile, "%%%% Reftime: %d %%%% Freq: %.16g %%%% RA: %.13g %%%% Dec: %.13g %%%% f1dot: %.13g %%%% f2dot: %.13g\n", dopplerParams_temp.refTime.gpsSeconds, dopplerParams_temp.fkdot[0], dopplerParams_temp.Alpha, dopplerParams_temp.Delta, dopplerParams_temp.fkdot[1], dopplerParams_temp.fkdot[2] );
      }
      XLAL_CHECK ( XLALComputeFstat(&Fstat_res, Fstat_in_vec->data[k], &dopplerParams_temp, 0.0, 1, FSTATQ_2F | FSTATQ_2F_PER_DET) == XLAL_SUCCESS, XLAL_EFUNC, "XLALComputeFstat() failed with errno=%d", xlalErrno );

      lineVeto->TwoF  += Fstat_res->twoF[0]; /* sum up multi-detector Fstat for this segment*/

      if ( singleSegStatsFile ) {
        fprintf ( singleSegStatsFile, "%.6f", Fstat_res->twoF[0] );
      }

      /* for each segment, number of detectors with data might be smaller than overall number */
      const UINT4 numDetectorsSeg = Fstat_res->numDetectors;

      /* recompute single-detector Fstats from atoms */
      for (UINT4 X = 0; X < numDetectorsSeg; X++)
        {

          /* match detector ID in this segment to one from detectorIDs list, sum up the corresponding twoFX */
          detid = -1;
          for (UINT4 Y = 0; Y < numDetectors; Y++) {
            if ( strcmp( Fstat_res->detectorNames[X], detectorIDs->data[Y] ) == 0 ) {
              detid = Y;
            }
          }
          if ( detid == -1 ) {
            XLALPrintError ("\nError in function %s, line %d : For segment k=%d, detector X=%d, could not match detector ID %s.\n\n",
                            __func__, __LINE__, k, X, Fstat_res->detectorNames[X] );
            XLAL_ERROR ( XLAL_EFAILED );
          }
          numSegmentsX[detid] += 1; /* have to keep this for correct averaging */

          twoFXseg->data[detid] = Fstat_res->twoFPerDet[X][0];

          lineVeto->TwoFX->data[detid]  += twoFXseg->data[detid]; /* sum up single-detector Fstat for this segment*/

        } /* for X < numDetectorsSeg */

      if ( singleSegStatsFile ) {
        for (UINT4 X = 0; X < numDetectors; X++) {
          fprintf ( singleSegStatsFile, " %.6f",twoFXseg->data[X] );
        }
        fprintf ( singleSegStatsFile, "\n" );
      }

    } /* for k < numSegments */

  /* get average stats over all segments */
  lineVeto->TwoF /= numSegments;
  for (UINT4 X = 0; X < numDetectors; X++) {
    lineVeto->TwoFX->data[X] /= numSegmentsX[X];
  }

  XLALDestroyREAL4Vector(twoFXseg);
  XLALDestroyFstatResults(Fstat_res);

  return(XLAL_SUCCESS);

} /* XLALComputeExtraStatsSemiCoherent() */
