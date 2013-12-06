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
#include "LineVeto.h"
#include <lal/TransientCW_utils.h> /* for XLALFastNegExp */

/*---------- local DEFINES ----------*/
#define TRUE (1==1)
#define FALSE (1==0)

/* Hooks for Einstein@Home / BOINC
   These are defined to do nothing special in the standalone case
   and will be set in boinc_extras.h if EAH_BOINC is set
*/
#ifdef HS_OPTIMIZATION
extern void
LocalComputeFStat ( LALStatus *status,
		    Fcomponents *Fstat,
		    const PulsarDopplerParams *doppler,
		    const MultiSFTVector *multiSFTs,
		    const MultiNoiseWeights *multiWeights,
		    const MultiDetectorStateSeries *multiDetStates,
		    const ComputeFParams *params,
		    ComputeFBuffer *cfBuffer);
#define COMPUTEFSTAT LocalComputeFStat
#else
#define COMPUTEFSTAT ComputeFStat
#endif

/*----- Macros ----- */
#define INIT_MEM(x) memset(&(x), 0, sizeof((x)))
#define SQUARE(x) ( (x) * (x) )

/*----- SWITCHES -----*/

/*---------- internal types ----------*/

/*---------- Global variables ----------*/

/* empty initializers  */
const LVcomponents empty_LVcomponents;

/*---------- internal prototypes ----------*/

/* ----- module-local fast lookup-table handling of negative exponentials ----- */

/**
 * Lookup-table for logarithms log(x)
 * Holds an array 'data' of 'length' for values log(x) for x in the range (0, xmax]
 */
#define LOGLUT_XMAX 	3.0	// LUT for range (0,numDetectors+1), currently numDetectors = 2 FIXME: get this dynamically
#define LOGLUT_LENGTH 	2000	// number of LUT values to pre-compute
static gsl_vector *logLUT = NULL; 	/**< module-global lookup-table for logarithms log(x) */
#define LOGLUT_DXINV  ((LOGLUT_LENGTH)/(LOGLUT_XMAX))	// 1/dx with dx = xmax/length

static int XLALCreateLogLUT ( void );	/* only ever used internally, destructor is in exported API */

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

  /* temporary copy of Fstatistic parameters structure, needed to change returnSingleF for function scope only */
  ComputeFParams CFparams_internal = (*CFparams);
  CFparams_internal.returnSingleF  = TRUE;

  /* initialize doppler parameters */
  candidateDopplerParams.refTime = refTimeGPS;  /* spin parameters in toplist refer to this refTime */

  /* initialise detector name vector for later identification */
  LALStringVector *detectorIDs;
  if ( ( detectorIDs = XLALGetDetectorIDs ( multiSFTsV )) == NULL ) /* fill detector name vector with all detectors present in any data sements */
    XLAL_ERROR ( XLAL_EFUNC );

  UINT4 numDetectors = detectorIDs->length;

  /* initialise LVcomponents structure and allocate memory */
  LVcomponents   lineVeto = empty_LVcomponents; /* struct containing multi-detector Fstat, single-detector Fstats, Line Veto stat */
  if ( (lineVeto.TwoFX = XLALCreateREAL4Vector ( numDetectors )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCreateREAL4Vector( %d )\n", __func__, numDetectors );
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
      } else if ( listEntryType == 2 ) {
        HoughFStatOutputEntry *elem = toplist_elem ( list, j );
        elemV = elem;

        if ( (elem->sumTwoFX = XLALCreateREAL4Vector ( numDetectors )) == NULL ) {
          XLALPrintError ("%s: failed to XLALCreateREAL4Vector( %d )\n", __func__, numDetectors );
          XLAL_ERROR ( XLAL_EFUNC );
        }

        /* get frequency, sky position, doppler parameters from toplist candidate and save to dopplerParams */
        candidateDopplerParams.Alpha = elem->AlphaBest;
        candidateDopplerParams.Delta = elem->DeltaBest;
        candidateDopplerParams.fkdot[0] = elem->Freq;
        candidateDopplerParams.fkdot[1] = elem->f1dot;
        /* no 2nd spindown in HoughFStatOutputEntry */
      } /* if listEntryType 2 */

      /* write header information into segment-Fstats file */
      if ( singleSegStatsFile )
        fprintf ( singleSegStatsFile, "%%%% Freq: %.16g\n%%%% RA: %.13g\n%%%% Dec: %.13g\n%%%% f1dot: %.13g\n%%%% f2dot: %.13g\n%%%% reftime: %d\n",
                  candidateDopplerParams.fkdot[0], candidateDopplerParams.Alpha, candidateDopplerParams.Delta, candidateDopplerParams.fkdot[1], candidateDopplerParams.fkdot[2], refTimeGPS.gpsSeconds );

      /*  recalculate multi- and single-IFO Fstats for all segments for this candidate */
      XLALComputeExtraStatsSemiCoherent( &lineVeto, &candidateDopplerParams, multiSFTsV, multiNoiseWeightsV, multiDetStatesV, detectorIDs, &CFparams_internal, SignalOnly, singleSegStatsFile );
      if ( xlalErrno != 0 ) {
        XLALPrintError ("\nError in function %s, line %d : Failed call to XLALComputeLineVetoSemiCoherent().\n\n", __func__, __LINE__);
        XLAL_ERROR ( XLAL_EFUNC );
      }

      /* save values in toplist */
      if ( listEntryType == 1 )
        {
          GCTtopOutputEntry *elem = elemV;
          elem->numDetectors = numDetectors;
          elem->sumTwoFrecalc  = lineVeto.TwoF;
          for ( X = 0; X < numDetectors; X ++ )
            elem->sumTwoFXrecalc[X]  = lineVeto.TwoFX->data[X];
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
  XLALDestroyREAL4Vector ( lineVeto.TwoFX );
  XLALDestroyStringVector ( detectorIDs );

  return (XLAL_SUCCESS);

} /* XLALComputeExtraStatsForToplist() */






/**
 * XLAL Function to recalculate single-IFO Fstats for all semicoherent search segments, and use them to compute Line Veto statistics
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


  REAL4Vector *twoFXseg = NULL;
  if ( (twoFXseg = XLALCreateREAL4Vector ( numDetectors )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCreateREAL4Vector( %d )\n", __func__, numDetectors );
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
        fprintf ( singleSegStatsFile, "%%%% Reftime: %d %%%% Freq: %.16g %%%% RA: %.13g %%%% Dec: %.13g %%%% f1dot: %.13g %%%% f2dot: %.13g\n", dopplerParams_temp.refTime.gpsSeconds, dopplerParams_temp.fkdot[0], dopplerParams_temp.Alpha, dopplerParams_temp.Delta, dopplerParams_temp.fkdot[1], dopplerParams_temp.fkdot[2] );
      fakeStatus = blank_status;
      COMPUTEFSTAT ( &fakeStatus, &Fstat, &dopplerParams_temp, multiSFTsV->data[k], multiNoiseWeightsThisSeg, multiDetStatesV->data[k], CFparams, NULL );
      if ( fakeStatus.statusCode ) {
        XLALPrintError ("\%s, line %d : Failed call to LAL function ComputeFStat(). statusCode=%d\n\n", __func__, __LINE__, fakeStatus.statusCode);
        XLAL_ERROR ( XLAL_EFUNC );
      }

      if ( SignalOnly ) {      /* normalization factor correction */
        Fstat.F *= 2.0 / Tsft;
        Fstat.F += 2;          /* E[F] = SNR^2/2 + 2 */
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

          twoFXseg->data[detid] = 2.0 * Fstat.FX[X];

          if ( SignalOnly ) {                      /* normalization factor correction */
            twoFXseg->data[detid] *= 2.0 / Tsft;
            twoFXseg->data[detid] += 4;            /* TwoF=2.0*F has been done before, so we have to add 2*2: E[2F] = SNR^2 + 4 */
          }

          lineVeto->TwoFX->data[detid]  += twoFXseg->data[detid]; /* sum up single-detector Fstat for this segment*/

        } /* for X < numDetectorsSeg */

      if ( singleSegStatsFile ) {
        for (X = 0; X < numDetectors; X++)
          fprintf ( singleSegStatsFile, " %.6f",twoFXseg->data[X] );
        fprintf ( singleSegStatsFile, "\n" );
      }

    } /* for k < numSegments */

  /* get average stats over all segments */
  lineVeto->TwoF /= numSegments;
  for (X = 0; X < numDetectors; X++) {
    lineVeto->TwoFX->data[X] /= numSegmentsX[X];
  }

  XLALDestroyREAL4Vector(twoFXseg);

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
  Fa = 0.0;
  Fb = 0.0;

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
      Fa += thisAtom->Fa_alpha;
      Fb += thisAtom->Fb_alpha;
    } /* loop through SFTs */

  } /* loop through detectors */

  /* compute determinant and final Fstat (not twoF!) */
  REAL8 Dinv = 1.0 / ( mmatrixA * mmatrixB - SQUARE(mmatrixC) );
  F = Dinv * ( mmatrixB * ( SQUARE(crealf(Fa)) + SQUARE(cimagf(Fa)) ) + mmatrixA * ( SQUARE(crealf(Fb)) + SQUARE(cimagf(Fb)) ) - 2.0 * mmatrixC * (crealf(Fa)*crealf(Fb) + cimagf(Fa)*cimagf(Fb)) );

  return(F);

} /* XLALComputeFstatFromAtoms() */




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
  if ( !TwoF || !TwoFXvec || !TwoFXvec->data )
    XLAL_ERROR_REAL4 ( XLAL_EFAULT, "Empty TwoF or TwoFX pointer as input parameter!\n\n");

  UINT4 numDetectors = TwoFXvec->length;

  if ( lXvec && ( lXvec->length != numDetectors ) )
    XLAL_ERROR_REAL4 ( XLAL_EBADLEN, "Input lX and TwoFX vectors have different length!\n\n" );

  if ( rhomaxline < 0 )
    XLAL_ERROR_REAL4 ( XLAL_EDOM, "Negative prior range 'rhomaxline' = %g! Must be >= 0!\n", rhomaxline );
  REAL8 logRhoTerm = 0.0;
  if ( rhomaxline > 0.0 )
   logRhoTerm = 4.0 * log(rhomaxline) - log(70.0);
  else /* if rhomaxline == 0.0, logRhoTerm should become irrelevant in summation */
    logRhoTerm = - LAL_REAL8_MAX;

  REAL8 *loglX = NULL;
  REAL8 loglXtemp[numDetectors];
  if ( lXvec ) {
    for (UINT4 X = 0; X < numDetectors; X++) {
      if ( lXvec->data[X] > 0 )
        loglXtemp[X] = log(lXvec->data[X]);
      else if ( lXvec->data[X] == 0 ) /* if zero prior ratio, approximate log(0)=-inf by -LAL_REA4_MAX to avoid raising underflow exceptions */
        loglXtemp[X] = - LAL_REAL8_MAX;
      else /* negative prior ratio is a mistake! */
       XLAL_ERROR_REAL4 ( XLAL_EDOM, "Negative input prior-ratio for detector X=%d: lX[X]=%g\n", X, lXvec->data[X] );
    }
    loglX = loglXtemp;
  }

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
  if ( !TwoF || !TwoFX )
    XLAL_ERROR_REAL4 ( XLAL_EFAULT, "Empty TwoF or TwoFX pointer as input parameter!\n\n");

  /* set up temporary variables and structs */
  REAL4 log0  = - LAL_REAL4_MAX;	/* approximates -inf */

  REAL4 maxInSum = log0;           /* keep track of largest summand in denominator, for logsumexp formula below */
  REAL4 FXprior[numDetectors];     /* FXprior equiv log(lX * e^(FX)) = FX + loglX */

  for (UINT4 X = 0; X < numDetectors; X++)
    {
      FXprior[X] = 0.5 * TwoFX[X];
      if (loglX) /* if no priors given, just use lX=1 => loglX=0 for all X => do not add anything */
        FXprior[X] += loglX[X];
      /* keep track of maximum value in denominator sum  */
      if ( FXprior[X] > maxInSum )
        maxInSum = FXprior[X];
    } /* for X < numDetectors */

  /* special treatment for additional denominator term 'rho^4/70' */
  if ( logRhoTerm > maxInSum )
    maxInSum = logRhoTerm;

  REAL4 LV = 0.0;	/* output variable for Line Veto statistics */

  LV = 0.5 * TwoF - maxInSum;	/* dominant term to LV-statistic */

  if ( useAllTerms )	/* optionally add logsumexp term (possibly negligible in many cases) */
    {
      REAL4 extraSum=0;	/* will be:  e^[-(maxInSum - logRhoTerm)] + sum_X e^[ -(maxInSum - FXprior) ] >= 1 */

      /* need to treat (rho^4/70) term separately */
      extraSum += exp ( logRhoTerm - maxInSum );

      /* now add all FX-contributions */
      for (UINT4 X = 0; X < numDetectors; X++)
        extraSum += exp ( FXprior[X] - maxInSum );

      LV -= log ( extraSum );

    } /* if useAllTerms */

  return LV;

} /* XLALComputeLineVetoArray() */


/**
 * XLAL function to get a list of detector IDs from multi-segment multiSFT vectors
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

  /* sort final list by detector-name */
  XLALSortStringVector ( IFOList );
  if ( xlalErrno != 0 ) {
    XLALPrintError ("\nError in function %s, line %d : Failed call to XLALSortStringVector().\n\n", __func__, __LINE__);
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }

  return IFOList;

} /* XLALGetDetectorIDs() */


/** Generate a lookup-table logLUT for log(x) over the interval x in (0, xmax], using 'length' points. */
int
XLALCreateLogLUT ( void )
{
  /* create empty output LUT */
  gsl_vector *ret;
  if ( ( ret = gsl_vector_alloc ( LOGLUT_LENGTH + 1)) == NULL ) {
    XLALPrintError ("%s: failed to gsl_vector_alloc (%s)\n", __func__, LOGLUT_LENGTH +1 );
    XLAL_ERROR ( XLAL_ENOMEM );
  }

  /* fill output LUT */
  REAL8 dx = LOGLUT_XMAX / LOGLUT_LENGTH;
  UINT4 i;
  for ( i=0; i <= LOGLUT_LENGTH; i ++ )
    {
      REAL8 xi = i * dx;

      gsl_vector_set ( ret, i, log( xi ) );

    } /* for i < length() */

  /* 'return' this by setting the global vector */
  logLUT = ret;

  return XLAL_SUCCESS;

} /* XLALCreateLogLUT() */

/**
 * Destructor function for logLUT_t lookup table
 */
void
XLALDestroyLogLUT ( void )
{
  if ( !logLUT )
    return;

  gsl_vector_free ( logLUT );

  logLUT = NULL;

  return;

} /* XLALDestroyLogLUT() */

/**
 * Fast logarithmic function log(x) using lookup-table (LUT).
 * We need to compute log(x) for x in (0,xmax], typically in a B-stat
 * integral of the form int e^-x dx: this means that small values e^(-x)
 * will not contribute much to the integral and are less important than
 * values close to 1. Therefore we pre-compute a LUT of e^(-x) for x in [0, xmax],
 * in Npoints points, and set e^(-x) = 0 for x < xmax.
 *
 * NOTE: if module-global logLUT=NULL, we create it here
 * NOTE: if argument is outside of (0,xmax], we use math-lib log(x) instead of LUT
 */
REAL8
XLALFastLog ( REAL8 x )
{
  if ( x > LOGLUT_XMAX )	/* for values bigger than xmax, use normal log function */
    return log(x);

  if ( x < 0 )
    XLAL_ERROR_REAL8 ( XLAL_EDOM, "Negative argument: x=%f\n", x );

  /* if lookup table doesn't exist yet: generate it now */
  if ( !logLUT && ( XLALCreateLogLUT() != XLAL_SUCCESS) ) {
    XLAL_ERROR_REAL8 ( XLAL_EFUNC );
  }

  /* find index of closest point xp in LUT to xm */
  UINT4 i0 = (UINT4) ( x * LOGLUT_DXINV + 0.5 );

  return gsl_vector_get ( logLUT, i0 );

} /* XLALFastLog() */
