/*
 *  Copyright (C) 2012, 2013 John Whelan, Shane Larson and Badri Krishnan
 *  Copyright (C) 2013, 2014 Badri Krishnan, John Whelan, Yuanhao Zhang
 *  Copyright (C) 2016, 2017 Grant David Meadors
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

#include <lal/PulsarCrossCorr_v2.h>

#define SQUARE(x) ((x)*(x))
#define QUAD(x) ((x)*(x)*(x)*(x))
#define GPSDIFF(x,y) (1.0*((x).gpsSeconds - (y).gpsSeconds) + ((x).gpsNanoSeconds - (y).gpsNanoSeconds)*1e-9)
#define MYMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MYMIN(x,y) ( (x) < (y) ? (x) : (y) )
#define USE_ALIGNED_MEMORY_ROUTINES
#define TRUE (1==1)
#define FALSE (1==0)


// ----- local prototypes ----------
static int
XLALApplyCrossCorrFreqShiftResamp
    ( 
     COMPLEX8                   *restrict xOut,
     const COMPLEX8TimeSeries   *restrict xIn,
     const PulsarDopplerParams  *restrict doppler,
     const REAL8                          freqShift,
     const UINT4                          indexStartResamp,
     const UINT4                          indexEndResamp,
     const UINT4                          numSamplesIn,
     const UINT4                          insertPoint
    );

static int
XLALComputeFaFb_CrossCorrResamp(
    ResampCrossCorrWorkspace                *restrict ws2L,
    COMPLEX8                                *restrict wsFaX_k,
    COMPLEX8                                *restrict wsFbX_k,
    const MultiResampSFTPairMultiIndexList  *restrict resampMultiPairs, 
    const MultiCOMPLEX8TimeSeries           *restrict multiTimeSeries_SRC_a,
    const MultiCOMPLEX8TimeSeries           *restrict multiTimeSeries_SRC_b, 
    const PulsarDopplerParams               *restrict dopplerpos,  
    const PulsarDopplerParams               *restrict binaryTemplateSpacings,
    const REAL8                                       SRCsampPerTcoh,
    const UINT4                                       detX,
    const UINT4                                       sftK,
    const UINT4                                       detY,
    const BOOLEAN                                     isL
);
// ==================== function definitions ====================



/** Calculate the Doppler-shifted frequency associated with each SFT in a list */
/* This is according to Eqns 2.11 and 2.12 of Dhurandhar et al 2008 */
/* Also returns the signal phase according to eqn 2.4 */
int XLALGetDopplerShiftedFrequencyInfo
  (
   REAL8Vector             *shiftedFreqs, /**< Output list of shifted frequencies */
   UINT4Vector               *lowestBins, /**< Output list of bin indices */
   COMPLEX8Vector       *expSignalPhases, /**< Output list of signal phases */
   REAL8VectorSequence         *sincList, /**< Output list of sinc factors */
   UINT4                         numBins, /**< Number of frequency bins to use */
   PulsarDopplerParams             *dopp, /**< Doppler parameters for signal */
   SFTIndexList              *sftIndices, /**< List of indices for SFTs */
   MultiSFTVector             *inputSFTs, /**< SFT data (needed for f0) */
   MultiSSBtimes             *multiTimes, /**< SSB or Binary times */
   MultiUINT4Vector             *badBins, /**< List of bin indices contaminated by known lines */
   REAL8                            Tsft  /**< SFT duration */
  )
{
  UINT4 numSFTs = sftIndices->length;
  if ( expSignalPhases->length !=numSFTs
       || shiftedFreqs->length !=numSFTs
       || lowestBins->length !=numSFTs
       || sincList->length !=numSFTs) {
    XLALPrintError("Lengths of SFT-indexed lists don't match!");
    XLAL_ERROR(XLAL_EBADLEN );
  }

  UINT4 numDets = inputSFTs->length;
  if ( multiTimes->length !=numDets
       || ( badBins && badBins->length !=numDets )
       ) {
    XLALPrintError("Lengths of detector-indexed lists don't match!");
    XLAL_ERROR(XLAL_EBADLEN );
  }

  if ( numBins < 1 ) {
    XLALPrintError("Must specify a positive number of bins to use!");
    XLAL_ERROR(XLAL_EBADLEN );
  }

  /* now calculate the intrinsic signal frequency in the SFT */
  /* fhat = f_0 + f_1(t-t0) + f_2(t-t0)^2/2 + ... */

  /* this is the sft reference time  - the pulsar reference time */
  for (UINT4 sftNum=0; sftNum < numSFTs; sftNum++) {
    UINT4 detInd = sftIndices->data[sftNum].detInd;
    XLAL_CHECK ( ( detInd < inputSFTs->length ),
		 XLAL_EINVAL,
		 "SFT asked for detector index off end of list:\n sftNum=%"LAL_UINT4_FORMAT", detInd=%"LAL_UINT4_FORMAT", inputSFTs->length=%d\n",
		 sftNum, detInd, inputSFTs->length );

    UINT4 numSFTsDet = inputSFTs->data[detInd]->length;
    SSBtimes *times;
    times = multiTimes->data[detInd];
    XLAL_CHECK ( ( times->DeltaT->length == numSFTsDet )
		 && ( times->Tdot->length == numSFTsDet ),
		 XLAL_EBADLEN,
		 "Lengths of multilists don't match!" );

    UINT4 sftInd = sftIndices->data[sftNum].sftInd;
    XLAL_CHECK ( ( sftInd < numSFTsDet ),
		 XLAL_EINVAL,
		 "SFT asked for SFT index off end of list:\n sftNum=%"LAL_UINT4_FORMAT", detInd=%"LAL_UINT4_FORMAT", sftInd=%"LAL_UINT4_FORMAT", numSFTsDet=%"LAL_UINT4_FORMAT"\n",
		 sftNum, detInd, sftInd, numSFTsDet );
    REAL8 timeDiff = times->DeltaT->data[sftInd]
      + XLALGPSDiff( &(times->refTime), &(dopp->refTime));
    REAL8 fhat = dopp->fkdot[0]; /* initialization */
    REAL8 phiByTwoPi = fmod ( fhat * timeDiff , 1.0 );
    REAL8 factor = timeDiff;

    for (UINT4 k = 1;  k < PULSAR_MAX_SPINS; k++) {
      fhat += dopp->fkdot[k] * factor;
      factor *= timeDiff / (k+1);
      phiByTwoPi += dopp->fkdot[k] * factor;
    }
    REAL4 sinPhi, cosPhi; /*Phi -> Phase of each SFT*/
    if(XLALSinCos2PiLUT(&sinPhi, &cosPhi, phiByTwoPi)!= XLAL_SUCCESS){
      LogPrintf ( LOG_CRITICAL, "%s: XLALSinCos2PiLUT() failed with errno=%d in XLALGetDopplerShiftedFrequencyInfo\n", __func__, xlalErrno );
      XLAL_ERROR( XLAL_EFUNC );
    }
    expSignalPhases->data[sftNum] = cosPhi + I * sinPhi;
    shiftedFreqs->data[sftNum] = fhat * times->Tdot->data[sftInd];
    REAL8 fminusf0 = shiftedFreqs->data[sftNum] - inputSFTs->data[detInd]->data[sftInd].f0;
    lowestBins->data[sftNum] = ceil( fminusf0 * Tsft - 0.5*numBins );
#define SINC_SAFETY 1e-5
    for (UINT4 l = 0; l < numBins; l++) {
      sincList->data[sftNum*numBins + l] = 1.0;
      if ( badBins && badBins->data[detInd] ) {
	for (UINT4 j = 0;
	     sincList->data[sftNum*numBins + l] != 0.0
	       && j < badBins->data[detInd]->length;
	     j++) {
	  if ( lowestBins->data[sftNum] + l
	       == badBins->data[detInd]->data[j] ) {
	    sincList->data[sftNum*numBins + l] = 0.0;
	  }
	}
      }

      if ( !badBins || !(badBins->data[detInd])
	   || sincList->data[sftNum*numBins + l] != 0.0 ) {
	/* Calculate normalized sinc, i.e., sin(pi*x)/(pi*x) */
	REAL4 sinPiX, cosPiX;
	REAL8 X = lowestBins->data[sftNum] - fminusf0 * Tsft + l;
	if(X > SINC_SAFETY || (X < - SINC_SAFETY)){
	  XLAL_CHECK( XLALSinCos2PiLUT( &sinPiX, &cosPiX, 0.5 * X ) == XLAL_SUCCESS, XLAL_EFUNC ); /*sin(2*pi*0.5*x)=sin(pi*x)*/
	  sincList->data[sftNum*numBins + l] = LAL_1_PI * sinPiX / X;/*1/(pi*x) =1/pi*1/x*/
	}
      }
    }

    /* printf("f=%.7f, f0=%.7f, Tsft=%g, numbins=%d, lowestbin=%d, kappa=%g\n",
	   shiftedFreqs->data[sftNum],
	   inputSFTs->data[detInd]->data[sftInd].f0,
	   Tsft, numBins, lowestBins->data[sftNum],
	   kappaValues->data[sftNum]); */

  }

  return XLAL_SUCCESS;

}

/** Construct flat SFTIndexList out of a MultiSFTVector */
/* Allocates memory as well */
int XLALCreateSFTIndexListFromMultiSFTVect
  (
   SFTIndexList        **indexList,   /* Output: flat list of indices to locate SFTs */
   MultiSFTVector            *sfts    /* Input: set of per-detector SFT vectors */
  )
{
  SFTIndexList *ret = NULL;

  UINT4 numDets = numDets = sfts->length;
  UINT4 numSFTs = 0;
  UINT4 j = 0;

  for (UINT4 k=0; k < numDets; k++) {
    numSFTs += sfts->data[k]->length;
  }

  if ( ( ret = XLALCalloc( 1, sizeof( *ret ) )) == NULL ) {
    XLAL_ERROR ( XLAL_ENOMEM );
  }
  ret->length = numSFTs;
  if ( ( ret->data = XLALCalloc ( numSFTs, sizeof ( *ret->data ) )) == NULL ) {
    XLALFree ( ret );
    XLAL_ERROR ( XLAL_ENOMEM );
  }

  for (UINT4 k=0; k < numDets; k++) {
    UINT4 numForDet = sfts->data[k]->length;
    for (UINT4 l=0; l < numForDet; l++) {
      ret->data[j].detInd = k;
      ret->data[j].sftInd = l;
      ++j;
    }
  }
  /* should sort list by GPS time if possible */
  /* qsort(ret->data, ret->length, sizeof(ret->data[0]), CompareGPSTime ) */

  (*indexList) = ret;

  return XLAL_SUCCESS;
}

/** Construct list of SFT pairs for inclusion in statistic */
/* Allocates memory as well */
int XLALCreateSFTPairIndexList
  (
   SFTPairIndexList  **pairIndexList, /* Output: list of SFT pairs */
   SFTIndexList           *indexList, /* Input: list of indices to locate SFTs */
   MultiSFTVector              *sfts, /* Input: set of per-detector SFT vectors */
   REAL8                      maxLag, /* Maximum allowed lag time */
   BOOLEAN              inclAutoCorr  /* Flag indicating whether a "pair" of an SFT with itself is allowed */
  )
{
  SFTPairIndexList *ret = NULL;

  UINT4 numSFTs = indexList->length;
  UINT4 j = 0;

  if ( ( ret = XLALCalloc( 1, sizeof( *ret ) )) == NULL ) {
    XLAL_ERROR ( XLAL_ENOMEM );
  }

  /* do two passes, one to count the number of pairs so the list can be allocated, and one to actually populate the list. */

  /* maximum possible number of pairs */

  for (UINT4 k=0; k < numSFTs; k++) {
    UINT4 lMin;
    if ( inclAutoCorr ) {
      lMin = k;
    } else {
      lMin = k+1;
    }
    LIGOTimeGPS gps1 = sfts->data[indexList->data[k].detInd]->data[indexList->data[k].sftInd].epoch;
    for (UINT4 l=lMin; l < numSFTs; l++) {
      LIGOTimeGPS gps2 = sfts->data[indexList->data[l].detInd]->data[indexList->data[l].sftInd].epoch;
      REAL8 timeDiff = XLALGPSDiff(&gps1,&gps2);
      if (fabs(timeDiff) <= maxLag) {
	++j;
      }
    }
  }
  ret->length = j;
  if ( ( ret->data = XLALCalloc ( ret->length , sizeof ( *ret->data ) )) == NULL ) {
    XLALFree ( ret );
    XLAL_ERROR ( XLAL_ENOMEM );
  }
  j = 0;
  for (UINT4 k=0; k < numSFTs; k++) {
    UINT4 lMin;
    if ( inclAutoCorr ) {
      lMin = k;
    } else {
      lMin = k+1;
    }
    LIGOTimeGPS gps1 = sfts->data[indexList->data[k].detInd]->data[indexList->data[k].sftInd].epoch;
    for (UINT4 l=lMin; l < numSFTs; l++) {
      LIGOTimeGPS gps2 = sfts->data[indexList->data[l].detInd]->data[indexList->data[l].sftInd].epoch;
      REAL8 timeDiff = XLALGPSDiff(&gps1,&gps2);
      if (fabs(timeDiff) <= maxLag) {
	ret->data[j].sftNum[0] = k;
	ret->data[j].sftNum[1] = l;
	++j;
      }
    }
  }


  (*pairIndexList) = ret;

  return XLAL_SUCCESS;
}


/** Resampling-modified: construct list of SFT pairs for inclusion in statistic */
/* Allocates memory as well */
int XLALCreateSFTPairIndexListResamp
  (
   MultiResampSFTPairMultiIndexList  **resampMultiPairIndexList, /**< [out] resamp list of SFT pairs */
   SFTPairIndexList  **pairIndexList, /**< [out] list of SFT pairs */
   SFTIndexList           *indexList, /**< [in] list of indices to locate SFTs */
   MultiSFTVector              *sfts, /**< [in] set of per-detector SFT vectors */
   REAL8                      maxLag, /**< [in] Maximum allowed lag time */
   BOOLEAN              inclAutoCorr, /**< [in] Flag indicating whether a "pair" of an SFT with itself is allowed */
   BOOLEAN          inclSameDetector, /**< [in] Flag indicating whether a "pair" of a detector with itself is allowed */
   REAL8                        Tsft, /**< [in] duration of a single SFT */
   REAL8                      Tshort  /**< [in] resampling Tshort */
  )
{
  /* Note, inclSameDetect not same as inclAutoCorr; is effectively "Yes" in demod search;
   * should make user variable*/

  SFTPairIndexList *ret = NULL;
  MultiResampSFTPairMultiIndexList *resampMultiRet = NULL;

  UINT4 numSFTs = indexList->length;
  UINT4 numDets = sfts->length;
  UINT4 j = 0;

  fprintf(stdout, "Number of detectors in SFT list: %u\n", numDets);
  if ( ( ret = XLALCalloc( 1, sizeof( *ret ) )) == NULL ) {
    XLAL_ERROR ( XLAL_ENOMEM );
  }
  if ( ( resampMultiRet = XLALCalloc( 1, sizeof( *resampMultiRet ) )) == NULL ) {
    XLAL_ERROR ( XLAL_ENOMEM );
  }

  /* Resamp block: */
  /* first, make an array to store the number of matching L */
  MultiResampSFTMultiCountList *MultiListOfLmatchingGivenMultiK;
 if ( ( MultiListOfLmatchingGivenMultiK = XLALCalloc( 1, sizeof( *MultiListOfLmatchingGivenMultiK ) )) == NULL ) {
    XLAL_ERROR ( XLAL_ENOMEM );
  }
  MultiListOfLmatchingGivenMultiK->length = numDets; /* Because it is a multi-multi-list */
  /* next, Allocate detectors X*/
  if ( ( MultiListOfLmatchingGivenMultiK->data = XLALCalloc ( MultiListOfLmatchingGivenMultiK->length, sizeof ( *MultiListOfLmatchingGivenMultiK->data ) )) == NULL ) {
    XLALFree ( MultiListOfLmatchingGivenMultiK);
    XLAL_ERROR ( XLAL_ENOMEM );
  }
  /* furthermore, scroll through each detector X to read how many SFTs K_X it has */
  for (UINT4 detX=0; detX < numDets; detX++){
    MultiListOfLmatchingGivenMultiK->data[detX].length = sfts->data[detX]->length;
    if ( ( MultiListOfLmatchingGivenMultiK->data[detX].data = XLALCalloc ( MultiListOfLmatchingGivenMultiK->data[detX].length , sizeof ( *MultiListOfLmatchingGivenMultiK->data[detX].data ) )) == NULL ) {
        XLALFree ( MultiListOfLmatchingGivenMultiK->data[detX].data);
        XLAL_ERROR ( XLAL_ENOMEM );
    }
    /* Now we can count how many matching detectors Y are available to K_X*/  
    /* Go through all possible Y given K_X */
    for (UINT4 k=0; k < MultiListOfLmatchingGivenMultiK->data[detX].length; k++){
      /* Note that we may have to adjust this to handle gaps in the data */
      MultiListOfLmatchingGivenMultiK->data[detX].data[k].length = numDets;
      if ( ( MultiListOfLmatchingGivenMultiK->data[detX].data[k].data = XLALCalloc ( MultiListOfLmatchingGivenMultiK->data[detX].data[k].length , sizeof ( *MultiListOfLmatchingGivenMultiK->data[detX].data[k].data ) )) == NULL ) {
      XLALFree ( MultiListOfLmatchingGivenMultiK->data[detX].data[k].data);
      XLAL_ERROR ( XLAL_ENOMEM );
      } 
    }
  }
  /* end resamp block */

  /* do two passes, one to count the number of pairs so the list can be allocated, and one to actually populate the list. */
  /* maximum possible number of pairs */
  for (UINT4 k=0; k < numSFTs; k++) {
    UINT4 lMin;
    if ( inclAutoCorr ) {
      lMin = k;
    } else {
      lMin = k+1;
    }
    LIGOTimeGPS gps1 = sfts->data[indexList->data[k].detInd]->data[indexList->data[k].sftInd].epoch;
    for (UINT4 l=lMin; l < numSFTs; l++) {
      LIGOTimeGPS gps2 = sfts->data[indexList->data[l].detInd]->data[indexList->data[l].sftInd].epoch;
      REAL8 timeDiff = XLALGPSDiff(&gps1,&gps2);
      if (fabs(timeDiff) <= maxLag) {
	++j;
      }
    }
    /* Assign the number of matching L to the list */ 
  }
  /* MULTI-RESAMPLING maximum possible number of pairs */
  UINT4 lMinMulti;
  for (UINT4 detX=0; detX < numDets; detX++){
    MultiListOfLmatchingGivenMultiK->data[detX].detInd = detX;
    for (UINT4 k=0; k <  MultiListOfLmatchingGivenMultiK->data[detX].length; k++){
      MultiListOfLmatchingGivenMultiK->data[detX].data[k].sftInd = k;
      for (UINT4 detY=0; detY < numDets; detY++){
        UINT4 LmatchingGivenKMulti = 0;
        if ( ((detY == detX ) && inclSameDetector) || (detY > detX)  ){
          if (detY == detX){
            if ( inclAutoCorr) {
              lMinMulti = k;
            } else {
              lMinMulti = k+1;
            }
          } else {
            lMinMulti = 0; 
          }
          LIGOTimeGPS gps1 = sfts->data[detX]->data[k].epoch;
          for (UINT4 l=lMinMulti; l < sfts->data[detY]->length; l++) {
            LIGOTimeGPS gps2 = sfts->data[detY]->data[l].epoch;
            REAL8 timeDiff = XLALGPSDiff(&gps1, &gps2);
            if (fabs(timeDiff) <= maxLag) {
              ++LmatchingGivenKMulti;
            }
          }
          MultiListOfLmatchingGivenMultiK->data[detX].data[k].data[detY].detInd = detY;
          MultiListOfLmatchingGivenMultiK->data[detX].data[k].data[detY].sftCount = LmatchingGivenKMulti;
        }
      }
    }
  }

  ret->length = j;
  if ( ( ret->data = XLALCalloc ( ret->length , sizeof ( *ret->data ) )) == NULL ) {
    XLALFree ( ret );
    XLAL_ERROR ( XLAL_ENOMEM );
  }
  j = 0;
  /* Now assign the lists of lists */

  /* Assign memory to the return structure based on what we found about how
   * many matching indices exist for the multi-multi struct */
  /* first, tell the resampling-return how many detectors there are */
  resampMultiRet->length = numDets;
  resampMultiRet->maxLag = maxLag;
  resampMultiRet->Tsft = Tsft;
  resampMultiRet->Tshort = Tshort;
  resampMultiRet->inclAutoCorr = inclAutoCorr;
  resampMultiRet->inclSameDetector = inclSameDetector;
  if ( ( resampMultiRet->data = XLALCalloc ( resampMultiRet->length , sizeof ( *resampMultiRet->data ) )) == NULL ) {
    XLALFree ( resampMultiRet );
    XLAL_ERROR ( XLAL_ENOMEM );
  }
  /* Now assign the multi-list of multi-lists */
  for (UINT4 detX=0; detX < resampMultiRet->length; detX++){
    resampMultiRet->data[detX].length = MultiListOfLmatchingGivenMultiK->data[detX].length; 
    /* Explicitly save detX index. Likely safe regardless because detX loops over numDets; could use detX */
    resampMultiRet->data[detX].detInd = MultiListOfLmatchingGivenMultiK->data[detX].detInd;
    if ( ( resampMultiRet->data[detX].data = XLALCalloc ( resampMultiRet->data[detX].length , sizeof ( *resampMultiRet->data[detX].data ) )) == NULL ) {
      XLALFree ( resampMultiRet->data[detX].data );
      XLAL_ERROR ( XLAL_ENOMEM );
    }    
    for (UINT4 k=0; k <  resampMultiRet->data[detX].length; k++){
      resampMultiRet->data[detX].data[k].length = MultiListOfLmatchingGivenMultiK->data[detX].data[k].length;
      /* Explicitly save K_X index. Likely safe regardless because K_X loops over all SFTs in X; could use k */
      resampMultiRet->data[detX].data[k].sftInd = MultiListOfLmatchingGivenMultiK->data[detX].data[k].sftInd; 
      if ( ( resampMultiRet->data[detX].data[k].data = XLALCalloc ( resampMultiRet->data[detX].data[k].length , sizeof ( *resampMultiRet->data[detX].data[k].data ) )) == NULL ) {
        XLALFree ( resampMultiRet->data[detX].data[k].data );
        XLAL_ERROR ( XLAL_ENOMEM );
      }
      for (UINT4 detY=0; detY < resampMultiRet->data[detX].data[k].length; detY++){
        resampMultiRet->data[detX].data[k].data[detY].length = MultiListOfLmatchingGivenMultiK->data[detX].data[k].data[detY].sftCount;
        /* Explicitly save L_K_X index. Needs explicit safety because not all detectors may be online */
        resampMultiRet->data[detX].data[k].data[detY].detInd = MultiListOfLmatchingGivenMultiK->data[detX].data[k].data[detY].detInd;
        if ( ( resampMultiRet->data[detX].data[k].data[detY].data = XLALCalloc ( resampMultiRet->data[detX].data[k].data[detY].length , sizeof ( *resampMultiRet->data[detX].data[k].data[detY].data ) )) == NULL ) {
          XLALFree ( resampMultiRet->data[detX].data[k].data[detY].data );
          XLAL_ERROR ( XLAL_ENOMEM );
        }
      }
    }
  }

  XLALDestroyMultiMatchList( MultiListOfLmatchingGivenMultiK );
 
  /* Ascertain and transcribe the matches into the return structure*/
  for (UINT4 k=0; k < numSFTs; k++) {
    UINT4 lMin;
    if ( inclAutoCorr ) {
      lMin = k;
    } else {
      lMin = k+1;
    }
    LIGOTimeGPS gps1 = sfts->data[indexList->data[k].detInd]->data[indexList->data[k].sftInd].epoch;
    for (UINT4 l=lMin; l < numSFTs; l++) {
      LIGOTimeGPS gps2 = sfts->data[indexList->data[l].detInd]->data[indexList->data[l].sftInd].epoch;
      REAL8 timeDiff = XLALGPSDiff(&gps1,&gps2);
      if (fabs(timeDiff) <= maxLag) {
	ret->data[j].sftNum[0] = k;
	ret->data[j].sftNum[1] = l;
	++j;
      }
    }
  }
  
  /* Multi-Ascertain and transcribe the matches into the return structure for resampling*/
  UINT4 allPairCounter = 0;
  UINT4 kFlatCounter = 0;
  UINT4 lFlatCounter = 0;
  for (UINT4 detX=0; detX < resampMultiRet->length; detX++){
    for (UINT4 k=0; k <  resampMultiRet->data[detX].length; k++){
      /* Reinitialize the l-flat counter because now we are
       * searching against the OTHER SFTs, with their own count */
      if ( inclAutoCorr ) {
        lFlatCounter = kFlatCounter;
      } else {
        lFlatCounter = kFlatCounter+1;
      }
      for (UINT4 detY=0; detY < resampMultiRet->data[detX].data[k].length; detY++){
        /* note that naively, resampMultiRet->data[detX].data[k].length = numDets,
         * so we search all detectors. This choice may need to be revisited when
         * we have gaps in the data */
        UINT4 outLmatchingGivenKMulti = 0;
        /* Repeating the following check is necessary, because in general
         * our index "l" is not the internal SFT vector index for detector Y.
         * This is in constrast to k, which is (because it is the first index).
         */
        if ( ((detY == detX ) && inclSameDetector) || (detY > detX)  ){
          if (detY == detX){
            if ( inclAutoCorr) {
              lMinMulti = k;
            } else {
              lMinMulti = k+1;
            }
          } else {
            lMinMulti = 0;
          }
          LIGOTimeGPS gps1 = sfts->data[detX]->data[k].epoch;
          for (UINT4 l=lMinMulti; l < sfts->data[detY]->length; l++) {
            LIGOTimeGPS gps2 = sfts->data[detY]->data[l].epoch;
            REAL8 timeDiff = XLALGPSDiff(&gps1, &gps2);
            if (fabs(timeDiff) <= maxLag) {
              /* Must refer to sft vectors, not the flat index*/
              resampMultiRet->data[detX].data[k].data[detY].data[outLmatchingGivenKMulti].sftInd = l;
              resampMultiRet->data[detX].data[k].data[detY].data[outLmatchingGivenKMulti].detInd = detY;
              resampMultiRet->data[detX].data[k].data[detY].data[outLmatchingGivenKMulti].flatInd = lFlatCounter;
              resampMultiRet->data[detX].data[k].data[detY].data[outLmatchingGivenKMulti].pairInd = allPairCounter;
              resampMultiRet->data[detX].data[k].data[detY].data[outLmatchingGivenKMulti].sciFlag = 1.0;
              ++outLmatchingGivenKMulti;
              ++allPairCounter;
            }
            ++lFlatCounter;
          } // end l
        } // similarity check: end
      } // end detY
      resampMultiRet->data[detX].data[k].flatInd = kFlatCounter;
      resampMultiRet->data[detX].data[k].sciFlag = 1.0;
      ++kFlatCounter;
    } // end k
  } // end detX
  resampMultiRet->oldPairCount = ret->length;
  resampMultiRet->allPairCount = allPairCounter;

  (*pairIndexList) = ret;
  (*resampMultiPairIndexList) = resampMultiRet;

  return XLAL_SUCCESS;
} // end XLALCreateSFTPairIndexListResamp


/** With Tshort, and Resampling-modified: construct list of SFT pairs for inclusion in statistic */
/* Allocates memory as well, with Tshort */
int XLALCreateSFTPairIndexListShortResamp
  (
   MultiResampSFTPairMultiIndexList  **        resampMultiPairIndexList, /**< [out] resamp list of SFT pairs */
   const REAL8                                 maxLag,                   /**< [in] Maximum allowed lag time */
   const BOOLEAN                               inclAutoCorr,             /**< [in] Flag indicating whether a "pair" of an SFT with itself is allowed */
   const BOOLEAN                               inclSameDetector,         /**< [in] Flag indicating whether a "pair" of a detector with itself is allowed */
   const REAL8                                 Tsft,                     /**< [in] duration of a single SFT */
   const MultiLIGOTimeGPSVector      *restrict multiTimes                /**< [in] timestamps containing Tshort times for each detector */
  )
{
  /* Note, inclSameDetector not same as inclAutoCorr; is effectively "Yes" in demod search;
   * have made user variable*/

  XLAL_CHECK ( multiTimes != NULL, XLAL_EINVAL );
  MultiResampSFTPairMultiIndexList *restrict resampMultiRet = NULL;

  const UINT4 numDets = multiTimes->length;
  const UINT4 numShortPerDet = multiTimes->data[0]->length;
  const UINT4 numSFTs = numShortPerDet*numDets;
  const REAL8 Tshort = multiTimes->data[0]->deltaT;

  if ( ( resampMultiRet = XLALCalloc( 1, sizeof( *resampMultiRet ) )) == NULL ) {
    XLAL_ERROR ( XLAL_ENOMEM );
  }

  /* First, let us provide a consistent zero-padded length to all
   *  detectors for their time series of Tshorts. We will then add
   *  supplementary information to indicate was is and is not in science */

  /* Resamp block: */
  /* first, make an array to store the number of matching L */
  MultiResampSFTMultiCountList *restrict MultiListOfLmatchingGivenMultiK;
 if ( ( MultiListOfLmatchingGivenMultiK = XLALCalloc( 1, sizeof( *MultiListOfLmatchingGivenMultiK ) )) == NULL ) {
    XLAL_ERROR ( XLAL_ENOMEM );
  }
  MultiListOfLmatchingGivenMultiK->length = numDets; /* Because it is a multi-multi-list */
  /* next, Allocate detectors X*/
  if ( ( MultiListOfLmatchingGivenMultiK->data = XLALCalloc ( MultiListOfLmatchingGivenMultiK->length, sizeof ( *MultiListOfLmatchingGivenMultiK->data ) )) == NULL ) {
    XLALFree ( MultiListOfLmatchingGivenMultiK);
    XLAL_ERROR ( XLAL_ENOMEM );
  }
  /* furthermore, scroll through each detector X to read how many SFTs K_X it has */
  for (UINT4 detX=0; detX < numDets; detX++){
    MultiListOfLmatchingGivenMultiK->data[detX].length = numShortPerDet;
    if ( ( MultiListOfLmatchingGivenMultiK->data[detX].data = XLALCalloc ( MultiListOfLmatchingGivenMultiK->data[detX].length , sizeof ( *MultiListOfLmatchingGivenMultiK->data[detX].data ) )) == NULL ) {
        XLALFree ( MultiListOfLmatchingGivenMultiK->data[detX].data);
        XLAL_ERROR ( XLAL_ENOMEM );
    }
    /* Now we can count how many matching detectors Y are available to K_X*/  
    /* Go through all possible Y given K_X */
    for (UINT4 k=0; k < MultiListOfLmatchingGivenMultiK->data[detX].length; k++){
      /* Note that we may have to adjust this to handle gaps in the data */
      MultiListOfLmatchingGivenMultiK->data[detX].data[k].length = numDets;
      if ( ( MultiListOfLmatchingGivenMultiK->data[detX].data[k].data = XLALCalloc ( MultiListOfLmatchingGivenMultiK->data[detX].data[k].length , sizeof ( *MultiListOfLmatchingGivenMultiK->data[detX].data[k].data ) )) == NULL ) {
      XLALFree ( MultiListOfLmatchingGivenMultiK->data[detX].data[k].data);
      XLAL_ERROR ( XLAL_ENOMEM );
      } 
    }
  }
  /* end resamp block */

  /* do two passes, one to count the number of pairs so the list can be allocated, and one to actually populate the list. */
  /* maximum possible number of pairs */
  UINT4 numMatchMultiResampLists = 0;
  /* MULTI-RESAMPLING maximum possible number of pairs */
  UINT4 lMinMulti;
  UINT4 lMaxMulti =0;
  for (UINT4 detX=0; detX < numDets; detX++){
    MultiListOfLmatchingGivenMultiK->data[detX].detInd = detX;
    for (UINT4 k=0; k <  MultiListOfLmatchingGivenMultiK->data[detX].length; k++){
      MultiListOfLmatchingGivenMultiK->data[detX].data[k].sftInd = k;
      for (UINT4 detY=0; detY < numDets; detY++){
        UINT4 LmatchingGivenKMulti = 0;
        if ( ((detY == detX ) && inclSameDetector) || (detY > detX)  ){
          if (detY == detX){
            if ( inclAutoCorr) {
              lMinMulti = k;
            } else {
              lMinMulti = k+1;
            }
          } else {
            lMinMulti = (UINT4)MYMAX(0, k - round(2*maxLag/Tshort) - 2); 
          }
          /* new style, for Tshort */
          LIGOTimeGPS gps1 = multiTimes->data[detX]->data[k];
          /* For speed */
          lMaxMulti = (UINT4)MYMIN(numShortPerDet, lMinMulti + round(4*maxLag/Tshort) + 4);
          for (UINT4 l=lMinMulti; l < lMaxMulti; l++) {
            LIGOTimeGPS gps2 = multiTimes->data[detY]->data[l];
            REAL8 timeDiff = XLALGPSDiff(&gps1, &gps2);
            if (fabs(timeDiff) <= maxLag) {
              ++LmatchingGivenKMulti;
              ++numMatchMultiResampLists;
            }
          }
          MultiListOfLmatchingGivenMultiK->data[detX].data[k].data[detY].detInd = detY;
          MultiListOfLmatchingGivenMultiK->data[detX].data[k].data[detY].sftCount = LmatchingGivenKMulti;
        }
      }
    }
  }

  /* Now assign the lists of lists */

  /* Assign memory to the return structure based on what we found about how
   * many matching indices exist for the multi-multi struct */
  /* first, tell the resampling-return how many detectors there are */
  resampMultiRet->length = numDets;
  resampMultiRet->maxLag = maxLag;
  resampMultiRet->Tsft = Tsft;
  resampMultiRet->Tshort = Tshort;
  resampMultiRet->inclAutoCorr = inclAutoCorr;
  resampMultiRet->inclSameDetector = inclSameDetector;
  if ( ( resampMultiRet->data = XLALCalloc ( resampMultiRet->length , sizeof ( *resampMultiRet->data ) )) == NULL ) {
    XLALFree ( resampMultiRet );
    XLAL_ERROR ( XLAL_ENOMEM );
  }
  /* Now assign the multi-list of multi-lists */
  for (UINT4 detX=0; detX < resampMultiRet->length; detX++){
    resampMultiRet->data[detX].length = MultiListOfLmatchingGivenMultiK->data[detX].length; 
    /* Explicitly save detX index. Likely safe regardless because detX loops over numDets; could use detX */
    resampMultiRet->data[detX].detInd = MultiListOfLmatchingGivenMultiK->data[detX].detInd;
    if ( ( resampMultiRet->data[detX].data = XLALCalloc ( resampMultiRet->data[detX].length , sizeof ( *resampMultiRet->data[detX].data ) )) == NULL ) {
      XLALFree ( resampMultiRet->data[detX].data );
      XLAL_ERROR ( XLAL_ENOMEM );
    }    
    for (UINT4 k=0; k <  resampMultiRet->data[detX].length; k++){
      resampMultiRet->data[detX].data[k].length = MultiListOfLmatchingGivenMultiK->data[detX].data[k].length;
      /* Explicitly save K_X index. Likely safe regardless because K_X loops over all SFTs in X; could use k */
      resampMultiRet->data[detX].data[k].sftInd = MultiListOfLmatchingGivenMultiK->data[detX].data[k].sftInd; 
      if ( ( resampMultiRet->data[detX].data[k].data = XLALCalloc ( resampMultiRet->data[detX].data[k].length , sizeof ( *resampMultiRet->data[detX].data[k].data ) )) == NULL ) {
        XLALFree ( resampMultiRet->data[detX].data[k].data );
        XLAL_ERROR ( XLAL_ENOMEM );
      }
      for (UINT4 detY=0; detY < resampMultiRet->data[detX].data[k].length; detY++){
        resampMultiRet->data[detX].data[k].data[detY].length = MultiListOfLmatchingGivenMultiK->data[detX].data[k].data[detY].sftCount;
        /* Explicitly save L_K_X index. Needs explicit safety because not all detectors may be online */
        resampMultiRet->data[detX].data[k].data[detY].detInd = MultiListOfLmatchingGivenMultiK->data[detX].data[k].data[detY].detInd;
        if ( ( resampMultiRet->data[detX].data[k].data[detY].data = XLALCalloc ( resampMultiRet->data[detX].data[k].data[detY].length , sizeof ( *resampMultiRet->data[detX].data[k].data[detY].data ) )) == NULL ) {
          XLALFree ( resampMultiRet->data[detX].data[k].data[detY].data );
          XLAL_ERROR ( XLAL_ENOMEM );
        }
      }
    }
  }

  XLALDestroyMultiMatchList( MultiListOfLmatchingGivenMultiK );
 
  /* Prepare the outer, flat list as well */
  SFTIndexList *restrict flatIndexList = NULL;
  if ( ( flatIndexList = XLALCalloc( 1, sizeof( *flatIndexList ) )) == NULL ) {
    XLAL_ERROR ( XLAL_ENOMEM );
  }
  flatIndexList->length = numSFTs;
  if ( ( flatIndexList->data = XLALCalloc ( numSFTs, sizeof ( *flatIndexList->data ) )) == NULL ) {
    XLALFree ( flatIndexList );
    XLAL_ERROR ( XLAL_ENOMEM );
  }

  /* Multi-Ascertain and transcribe the matches into the return structure for resampling*/
  UINT4 allPairCounter = 0;
  UINT4 kFlatCounter = 0;
  UINT4 lFlatCounter = 0;
  for (UINT4 detX=0; detX < resampMultiRet->length; detX++){
    for (UINT4 k=0; k <  resampMultiRet->data[detX].length; k++){
      kFlatCounter = numShortPerDet * detX + k;
      /* Reinitialize the l-flat counter because now we are
         searching against the OTHER SFTs, with their own count */
      for (UINT4 detY=0; detY < resampMultiRet->data[detX].data[k].length; detY++){
        /* note that naively, resampMultiRet->data[detX].data[k].length = numDets,
         * so we search all detectors. This choice may need to be revisited when
         * we have gaps in the data */
        UINT4 outLmatchingGivenKMulti = 0;
        /* Repeating the following check is necessary, because in general
         * our index "l" is not the internal SFT vector index for detector Y.
         * This is in contrast to k, which is (because it is the first index).
         */
        if ( ((detY == detX ) && inclSameDetector) || (detY > detX)  ){
          if (detY == detX){
            if ( inclAutoCorr) {
              lMinMulti = k;
            } else {
              lMinMulti = k+1;
            }
          } else {
            lMinMulti = (UINT4)MYMAX(0, k - round(2*maxLag/Tshort) - 2);
          }
          /* new style, for Tshort */
          LIGOTimeGPS gps1 = multiTimes->data[detX]->data[k];
          lMaxMulti = (UINT4)MYMIN(numShortPerDet, lMinMulti + round(4*maxLag/Tshort) + 4);
          for (UINT4 l=lMinMulti; l < lMaxMulti; l++) {
            LIGOTimeGPS gps2 = multiTimes->data[detY]->data[l];
            REAL8 timeDiff = XLALGPSDiff(&gps1, &gps2);
            if (fabs(timeDiff) <= maxLag) {
              resampMultiRet->data[detX].data[k].data[detY].data[outLmatchingGivenKMulti].sftInd = l;
              resampMultiRet->data[detX].data[k].data[detY].data[outLmatchingGivenKMulti].detInd = detY;
              lFlatCounter = numShortPerDet * detY + l;
              resampMultiRet->data[detX].data[k].data[detY].data[outLmatchingGivenKMulti].flatInd = lFlatCounter;
              /* note that pairInd = alpha in the Whelan et al 2015 paper's notation */
              resampMultiRet->data[detX].data[k].data[detY].data[outLmatchingGivenKMulti].pairInd = allPairCounter;
              resampMultiRet->data[detX].data[k].data[detY].data[outLmatchingGivenKMulti].sciFlag = 1.0;
              ++outLmatchingGivenKMulti;
              ++allPairCounter;
            }
          } // end l
        } // similarity check: end
      } // end detY
      resampMultiRet->data[detX].data[k].flatInd = kFlatCounter;
      resampMultiRet->data[detX].data[k].sciFlag = 1.0;
      flatIndexList->data[kFlatCounter].detInd = detX;
      /* Note, the SFT index is probably just k, but best to be safe */
      flatIndexList->data[kFlatCounter].sftInd = resampMultiRet->data[detX].data[k].sftInd;
      /* kFlatCounter should serve as a flat index of all SFTs, or Shorts, 
       * because k passes through all, for all detectors, even if a given
       * one does not match anything (perhaps b/c it was already assigned
       * to all appropriate pairs */
    } // end k
  } // end detX
  resampMultiRet->allPairCount = allPairCounter;
  /* Add one to counts because this is C, with 0-indexing;
   * remember we assigned kFlatCounter not incrementally but directly  */
  resampMultiRet->sftTotalCount = 1+ kFlatCounter;
  resampMultiRet->indexList = flatIndexList;

  /* Prepare the standard, old-school pair list too */
  SFTPairIndexList *restrict standardPairIndexList = NULL;
  if ( ( standardPairIndexList = XLALCalloc( 1, sizeof( *standardPairIndexList ) )) == NULL ) {
    XLAL_ERROR ( XLAL_ENOMEM );
  }
  standardPairIndexList->length = resampMultiRet->allPairCount;
  if ( ( standardPairIndexList->data = XLALCalloc ( standardPairIndexList->length , sizeof ( *standardPairIndexList->data ) )) == NULL ) {
    XLALFree ( standardPairIndexList );
    XLAL_ERROR ( XLAL_ENOMEM );
  }

  UINT4 sftNum1;
  UINT4 sftNum2;

  UINT4 standardPairCount = 0;
  for (UINT4 detX=0; detX < resampMultiRet->length; detX++){
    for (UINT4 k=0; k <  resampMultiRet->data[detX].length; k++){
      for (UINT4 detY=0; detY < resampMultiRet->data[detX].data[k].length; detY++){
        for (UINT4 l=0; l < resampMultiRet->data[detX].data[k].data[detY].length; l++){
          /* Here we must us the flat SFT number */
          sftNum1 = numShortPerDet * detX + resampMultiRet->data[detX].data[k].sftInd;
          sftNum2 = numShortPerDet * detY + resampMultiRet->data[detX].data[k].data[detY].data[l].sftInd;
          standardPairIndexList->data[standardPairCount].sftNum[0] = sftNum1;
          standardPairIndexList->data[standardPairCount].sftNum[1] = sftNum2;
          ++standardPairCount;
        }
      }
    }
  }
  resampMultiRet->pairIndexList = standardPairIndexList;
  resampMultiRet->oldPairCount = standardPairCount;
  (*resampMultiPairIndexList) = resampMultiRet;

  return XLAL_SUCCESS;
} // end XLALCreateSFTPairIndexListShortResamp

/** Check that the contents of a resampling multi-pair index list are sensible by inspection */
int XLALTestResampPairIndexList
  (
    MultiResampSFTPairMultiIndexList  * resampMultiPairIndexList /**< [in] resamp list of SFT pairs */
  )
{
  MultiResampSFTPairMultiIndexList *resampMultiRet = resampMultiPairIndexList;
  /* SANITY CHECK */
  UINT4 grandTotalCounter = 0;
  for (UINT4 detX=0; detX < resampMultiRet->length; detX++){
    for (UINT4 k=0; k <  resampMultiRet->data[detX].length; k++){
      for (UINT4 detY=0; detY < resampMultiRet->data[detX].data[k].length; detY++){
        for (UINT4 l=0; l < resampMultiRet->data[detX].data[k].data[detY].length; l++){
          ++grandTotalCounter;
          fprintf(stdout, "SFT INDEX: %u\n", resampMultiRet->data[detX].data[k].data[detY].data[l].sftInd);
          fprintf(stdout, "det INDEX: %u\n", resampMultiRet->data[detX].data[k].data[detY].data[l].detInd);
          fprintf(stdout, "pair INDEX: %u\n", resampMultiRet->data[detX].data[k].data[detY].data[l].pairInd);
        }
      }
    }
  }
  fprintf(stdout, "Sanity check: GRAND TOTAL of passes through loop: %u\n", grandTotalCounter);
  fprintf(stdout, "Pairing completed\n");
  /* SANITY CHECK 2 */
  UINT4 newKFlatCounter = 0;
  for (UINT4 detX=0; detX < resampMultiRet->length; detX++){
    for (UINT4 k=0; k <  resampMultiRet->data[detX].length; k++){
      printf("detInd: %u\n", resampMultiRet->indexList->data[newKFlatCounter].detInd);
      printf("sftInd: %u\n", resampMultiRet->indexList->data[newKFlatCounter].sftInd);
      newKFlatCounter++;
    }
  }
  printf("newKFlatCounter: %u\n", newKFlatCounter);
  return XLAL_SUCCESS;
} // end XLALTestResampPairIndexList


/** Construct vector of G_alpha amplitudes for each SFT pair */
/* This is averaged over unknown cosi and psi */
/* Allocates memory as well */
int XLALCalculateCrossCorrGammas
  (
   REAL8Vector          **Gamma_ave, /* Output: vector of aa+bb values */
   REAL8Vector         **Gamma_circ, /* Output: vector of ab-ba values */
   SFTPairIndexList  *pairIndexList, /* Input: list of SFT pairs */
   SFTIndexList          *indexList, /* Input: list of SFTs */
   MultiAMCoeffs       *multiCoeffs  /* Input: AM coefficients */
  )
{

  UINT4 numPairs = pairIndexList->length;

  REAL8Vector *ret1 = NULL;
  XLAL_CHECK ( ( ret1 = XLALCreateREAL8Vector ( numPairs ) ) != NULL, XLAL_EFUNC, "XLALCreateREAL8Vector ( %"LAL_UINT4_FORMAT" ) failed.", numPairs );
  REAL8Vector *ret2 = NULL;
  XLAL_CHECK ( ( ret2 = XLALCreateREAL8Vector ( numPairs ) ) != NULL, XLAL_EFUNC, "XLALCreateREAL8Vector ( %"LAL_UINT4_FORMAT" ) failed.", numPairs );

  for (UINT4 j=0; j < numPairs; j++) {
    UINT4 sftNum1 = pairIndexList->data[j].sftNum[0];
    UINT4 sftNum2 = pairIndexList->data[j].sftNum[1];
    UINT4 detInd1 = indexList->data[sftNum1].detInd;
    UINT4 detInd2 = indexList->data[sftNum2].detInd;
    UINT4 sftInd1 = indexList->data[sftNum1].sftInd;
    UINT4 sftInd2 = indexList->data[sftNum2].sftInd;
    ret1->data[j] = 0.1 * ( multiCoeffs->data[detInd1]->a->data[sftInd1]
			   * multiCoeffs->data[detInd2]->a->data[sftInd2]
			   + multiCoeffs->data[detInd1]->b->data[sftInd1]
			   * multiCoeffs->data[detInd2]->b->data[sftInd2] );
    ret2->data[j] = 0.1 * ( multiCoeffs->data[detInd1]->a->data[sftInd1]
			   * multiCoeffs->data[detInd2]->b->data[sftInd2]
			   - multiCoeffs->data[detInd1]->b->data[sftInd1]
			   * multiCoeffs->data[detInd2]->a->data[sftInd2] );
  }

  (*Gamma_ave) = ret1;
  (*Gamma_circ) = ret2;
  return XLAL_SUCCESS;
}

/** test function for RESAMPLING */
/** Construct multi-dimensional array of G_L_Y_K_X amplitudes for each SFT pair */
/* This is averaged over unknwon cosi and psi */
/* Allocates memory as well */
int XLALCalculateCrossCorrGammasResamp
  (
   REAL8Vector          **Gamma_ave, /**< Output: vector of aa+bb values */
   REAL8Vector         **Gamma_circ, /**< Output: vector of ab-ba values */
   MultiResampSFTPairMultiIndexList  *resampMultiPairIndexList, /**< Input: resamp list of SFT pairs */
   MultiAMCoeffs       *multiCoeffs  /**< Input: AM coefficients */
  )
{
  UINT4 numPairs = resampMultiPairIndexList->allPairCount;

  REAL8Vector *ret1 = NULL;
  XLAL_CHECK ( ( ret1 = XLALCreateREAL8Vector ( numPairs ) ) != NULL, XLAL_EFUNC, "XLALCreateREAL8Vector ( %"LAL_UINT4_FORMAT" ) failed.", numPairs );
  REAL8Vector *ret2 = NULL;
  XLAL_CHECK ( ( ret2 = XLALCreateREAL8Vector ( numPairs ) ) != NULL, XLAL_EFUNC, "XLALCreateREAL8Vector ( %"LAL_UINT4_FORMAT" ) failed.", numPairs );


  /* beware segfault from stack overflow here for large numPairs */
  UINT4 detInd1[numPairs];
  UINT4 detInd2[numPairs];
  UINT4 sftInd1[numPairs];
  UINT4 sftInd2[numPairs];

  UINT4 pairCount = 0;
  for (UINT4 detX=0; detX < resampMultiPairIndexList->length; detX++){
    for (UINT4 k=0; k <  resampMultiPairIndexList->data[detX].length; k++){
      for (UINT4 detY=0; detY < resampMultiPairIndexList->data[detX].data[k].length; detY++){
        for (UINT4 l=0; l < resampMultiPairIndexList->data[detX].data[k].data[detY].length; l++){
          /* This way works if multiCoeffs is allocated for SFTs wrt original
           * vector (as is the case) */
          detInd1[pairCount] = resampMultiPairIndexList->data[detX].detInd;
          detInd2[pairCount] = resampMultiPairIndexList->data[detX].data[k].data[detY].detInd;
          sftInd1[pairCount] = resampMultiPairIndexList->data[detX].data[k].sftInd;
          sftInd2[pairCount] = resampMultiPairIndexList->data[detX].data[k].data[detY].data[l].sftInd;
          ++pairCount;
        }
      }
    }
  }
  for (UINT4 j=0; j < numPairs; j++) {
    ret1->data[j] = 0.1 * (  multiCoeffs->data[detInd1[j]]->a->data[sftInd1[j]]
                           * multiCoeffs->data[detInd2[j]]->a->data[sftInd2[j]]
                           + multiCoeffs->data[detInd1[j]]->b->data[sftInd1[j]]
                           * multiCoeffs->data[detInd2[j]]->b->data[sftInd2[j]] );
    ret2->data[j] = 0.1 * (  multiCoeffs->data[detInd1[j]]->a->data[sftInd1[j]]
                           * multiCoeffs->data[detInd2[j]]->b->data[sftInd2[j]]
                           - multiCoeffs->data[detInd1[j]]->b->data[sftInd1[j]]
                           * multiCoeffs->data[detInd2[j]]->a->data[sftInd2[j]] );
  }

  (*Gamma_ave) = ret1;
  (*Gamma_circ) = ret2;
  return XLAL_SUCCESS;
} // end XLALCalculateCrossCorrGammasResamp

/** test function for RESAMPLING with tShort*/
/** Construct multi-dimensional array of G_L_Y_K_X amplitudes for each SFT pair */
/* This is averaged over unknwon cosi and psi */
/* Allocates memory as well */
int XLALCalculateCrossCorrGammasResampShort
  (
   REAL8Vector          **Gamma_ave, /**< Output: vector of aa+bb values */
   REAL8Vector         **Gamma_circ, /**< Output: vector of ab-ba values */
   MultiResampSFTPairMultiIndexList  *resampMultiPairIndexList, /**< Input: resamp list of SFT pairs */
   MultiAMCoeffs       *multiCoeffs  /**< Input: AM coefficients */
  )
{
  UINT4 numPairs = resampMultiPairIndexList->allPairCount;

  REAL8Vector *ret1 = NULL;
  XLAL_CHECK ( ( ret1 = XLALCreateREAL8Vector ( numPairs ) ) != NULL, XLAL_EFUNC, "XLALCreateREAL8Vector ( %"LAL_UINT4_FORMAT" ) failed.", numPairs );
  REAL8Vector *ret2 = NULL;
  XLAL_CHECK ( ( ret2 = XLALCreateREAL8Vector ( numPairs ) ) != NULL, XLAL_EFUNC, "XLALCreateREAL8Vector ( %"LAL_UINT4_FORMAT" ) failed.", numPairs );

  REAL4 eleMultiCoeffs1a = 0;
  REAL4 eleMultiCoeffs2a = 0;
  REAL4 eleMultiCoeffs1b = 0;
  REAL4 eleMultiCoeffs2b = 0;

  REAL4 sciMultiCoeffs1a = 0;
  REAL4 sciMultiCoeffs2a = 0;
  REAL4 sciMultiCoeffs1b = 0;
  REAL4 sciMultiCoeffs2b = 0;

  /* beware segfault from stack overflow here for large numPairs */
  REAL4 castSciFlag1[numPairs];
  REAL4 castSciFlag2[numPairs];

  UINT4 detInd1[numPairs];
  UINT4 detInd2[numPairs];
  UINT4 sftInd1[numPairs];
  UINT4 sftInd2[numPairs];

  UINT4 pairCount = 0;
  for (UINT4 detX=0; detX < resampMultiPairIndexList->length; detX++){
    for (UINT4 k=0; k <  resampMultiPairIndexList->data[detX].length; k++){
      for (UINT4 detY=0; detY < resampMultiPairIndexList->data[detX].data[k].length; detY++){
        for (UINT4 l=0; l < resampMultiPairIndexList->data[detX].data[k].data[detY].length; l++){
          /* This way works if multiCoeffs is allocated for SFTs wrt original
           * vector (as is the case) */
          detInd1[pairCount] = resampMultiPairIndexList->data[detX].detInd;
          detInd2[pairCount] = resampMultiPairIndexList->data[detX].data[k].data[detY].detInd;
          sftInd1[pairCount] = resampMultiPairIndexList->data[detX].data[k].sftInd;
          sftInd2[pairCount] = resampMultiPairIndexList->data[detX].data[k].data[detY].data[l].sftInd;
          /* Can use the science flags; will produce slightly different 
           * results due to sharper edges than A and B time series;
           * however, by themselves, A and B are zero-padded and so
           * can handle gaps. */
          //castSciFlag1[pairCount] = resampMultiPairIndexList->data[detX].data[k].sciFlag;
          //castSciFlag2[pairCount] = resampMultiPairIndexList->data[detX].data[k].data[detY].data[l].sciFlag;
          castSciFlag1[pairCount] = 1.0;
          castSciFlag2[pairCount] = 1.0;
          ++pairCount;
        }
      }
    }
  }
  for (UINT4 j=0; j < numPairs; j++) {
    eleMultiCoeffs1a = multiCoeffs->data[detInd1[j]]->a->data[sftInd1[j]];
    eleMultiCoeffs2a = multiCoeffs->data[detInd2[j]]->a->data[sftInd2[j]];
    eleMultiCoeffs1b = multiCoeffs->data[detInd1[j]]->b->data[sftInd1[j]];
    eleMultiCoeffs2b = multiCoeffs->data[detInd2[j]]->b->data[sftInd2[j]];

    sciMultiCoeffs1a = eleMultiCoeffs1a * castSciFlag1[j];
    sciMultiCoeffs2a = eleMultiCoeffs2a * castSciFlag2[j];
    sciMultiCoeffs1b = eleMultiCoeffs1b * castSciFlag1[j];
    sciMultiCoeffs2b = eleMultiCoeffs2b * castSciFlag2[j];
    ret1->data[j] = 0.1 * (  sciMultiCoeffs1a
			   * sciMultiCoeffs2a
			   + sciMultiCoeffs1b
			   * sciMultiCoeffs2b );
    ret2->data[j] = 0.1 * (  sciMultiCoeffs1a
			   * sciMultiCoeffs2b
			   - sciMultiCoeffs1b
			   * sciMultiCoeffs2a );
  }

  (*Gamma_ave) = ret1;
  (*Gamma_circ) = ret2;
  return XLAL_SUCCESS;
} // end XLALCalculateCrossCorrGammasResampShort


/** Calculate multi-bin cross-correlation statistic */
/* This assumes rectangular or nearly-rectangular windowing */
int XLALCalculatePulsarCrossCorrStatistic
(
 REAL8                         *ccStat, /* Output: cross-correlation statistic rho */
 REAL8                      *evSquared, /* Output: (E[rho]/h0^2)^2 */
 REAL8Vector                *curlyGAmp, /* Input: Amplitude of curly G for each pair */
 COMPLEX8Vector       *expSignalPhases, /* Input: Phase of signal for each SFT */
 UINT4Vector               *lowestBins, /* Input: Bin index to start with for each SFT */
 REAL8VectorSequence         *sincList, /* Input: input the sinc factors*/
 SFTPairIndexList            *sftPairs, /* Input: flat list of SFT pairs */
 SFTIndexList              *sftIndices, /* Input: flat list of SFTs */
 MultiSFTVector             *inputSFTs, /* Input: SFT data */
 MultiNoiseWeights       *multiWeights, /* Input: nomalizeation factor S^-1 & weights for each SFT */
 UINT4                         numBins  /* Input Number of frequency bins to be taken into calc */
 )
{

  UINT4 numSFTs = sftIndices->length;
  if ( expSignalPhases->length !=numSFTs
       || lowestBins->length !=numSFTs
       || sincList->length !=numSFTs) {
    XLALPrintError("Lengths of SFT-indexed lists don't match!");
    XLAL_ERROR(XLAL_EBADLEN );
  }

  UINT4 numPairs = sftPairs->length;
  if ( curlyGAmp->length !=numPairs ) {
    XLALPrintError("Lengths of pair-indexed lists don't match!");
    XLAL_ERROR(XLAL_EBADLEN );
  }
  REAL8 nume = 0;
  REAL8 curlyGSqr = 0;
  *ccStat = 0.0;
  *evSquared = 0.0;
  for (UINT4 alpha = 0; alpha < numPairs; alpha++) {
    UINT4 sftNum1 = sftPairs->data[alpha].sftNum[0];
    UINT4 sftNum2 = sftPairs->data[alpha].sftNum[1];

    XLAL_CHECK ( ( sftNum1 < numSFTs ) && ( sftNum2 < numSFTs ),
		 XLAL_EINVAL,
		 "SFT pair asked for SFT index off end of list:\n alpha=%"LAL_UINT4_FORMAT", sftNum1=%"LAL_UINT4_FORMAT", sftNum2=%"LAL_UINT4_FORMAT", numSFTs=%"LAL_UINT4_FORMAT"\n",
		 alpha,  sftNum1, sftNum2, numSFTs );

    UINT4 detInd1 = sftIndices->data[sftNum1].detInd;
    UINT4 detInd2 = sftIndices->data[sftNum2].detInd;

    XLAL_CHECK ( ( detInd1 < inputSFTs->length )
		 && ( detInd2 < inputSFTs->length ),
		 XLAL_EINVAL,
		 "SFT asked for detector index off end of list:\n sftNum1=%"LAL_UINT4_FORMAT", sftNum2=%"LAL_UINT4_FORMAT", detInd1=%"LAL_UINT4_FORMAT", detInd2=%"LAL_UINT4_FORMAT", inputSFTs->length=%d\n",
		 sftNum1, sftNum2, detInd1, detInd2, inputSFTs->length );

    UINT4 sftInd1 = sftIndices->data[sftNum1].sftInd;
    UINT4 sftInd2 = sftIndices->data[sftNum2].sftInd;

    XLAL_CHECK ( ( sftInd1 < inputSFTs->data[detInd1]->length )
		 && ( sftInd2 < inputSFTs->data[detInd2]->length ),
		 XLAL_EINVAL,
		 "SFT asked for SFT index off end of list:\n sftNum1=%"LAL_UINT4_FORMAT", sftNum2=%"LAL_UINT4_FORMAT", detInd1=%"LAL_UINT4_FORMAT", detInd2=%"LAL_UINT4_FORMAT", sftInd1=%"LAL_UINT4_FORMAT", sftInd2=%"LAL_UINT4_FORMAT", inputSFTs->data[detInd1]->length=%d, inputSFTs->data[detInd2]->length=%d\n",
		 sftNum1, sftNum2, detInd1, detInd2, sftInd1, sftInd2,
		 inputSFTs->data[detInd1]->length,
		 inputSFTs->data[detInd2]->length );

    COMPLEX8 *dataArray1 = inputSFTs->data[detInd1]->data[sftInd1].data->data;
    COMPLEX8 *dataArray2 = inputSFTs->data[detInd2]->data[sftInd2].data->data;
    UINT4 lenDataArray1 = inputSFTs->data[detInd1]->data[sftInd1].data->length;
    UINT4 lenDataArray2 = inputSFTs->data[detInd2]->data[sftInd2].data->length;
    COMPLEX8 GalphaCC = curlyGAmp->data[alpha]
      * expSignalPhases->data[sftNum1]
      * conj( expSignalPhases->data[sftNum2] );
    INT4 baseCCSign = 1; /* Alternating sign is (-1)**(k1-k2) */
    if ( ( (lowestBins->data[sftNum1]-lowestBins->data[sftNum2]) % 2) != 0 ) {
      baseCCSign = -1;
    }

    UINT4 lowestBin1 = lowestBins->data[sftNum1];
    XLAL_CHECK ( ((lowestBin1 + numBins - 1) < lenDataArray1),
		 XLAL_EINVAL,
		 "Loop would run off end of array:\n lowestBin1=%d, numBins=%d, len(dataArray1)=%d\n",
		 lowestBin1, numBins, lenDataArray1 );
    for (UINT4 j = 0; j < numBins; j++) {
      COMPLEX8 data1 = dataArray1[lowestBin1+j];

      INT4 ccSign = baseCCSign;
      UINT4 lowestBin2 = lowestBins->data[sftNum2];
      XLAL_CHECK ( ((lowestBin2 + numBins - 1) < lenDataArray2),
		   XLAL_EINVAL,
		   "Loop would run off end of array:\n lowestBin2=%d, numBins=%d, len(dataArray2)=%d\n",
		   lowestBin2, numBins, lenDataArray2 );
      for (UINT4 k = 0; k < numBins; k++) {
	COMPLEX8 data2 = dataArray2[lowestBins->data[sftNum2]+k];
	REAL8 sincFactor =1;

	sincFactor = sincList->data[sftNum1 * numBins + j] * sincList->data[sftNum2 * numBins + k];
	nume +=  ccSign * sincFactor * creal(GalphaCC * conj(data1) * data2);
	/*nume += creal ( GalphaCC * ccSign * sincFactor * conj(data1) * data2 );=> abs(GalphaCC)*abs(data1)*abs(data2) * cos(arg(GalphaCC)-arg(data1)+arg(data2))*/
	/*multiWeights->data[detInd1]->data[sftNum1] *  multiWeights->data[detInd2]->data[sftNum2] **/
	REAL8 GalphaAmp = curlyGAmp->data[alpha] * sincFactor ;
	/** multiWeights->data[detInd1]->data[sftNum1] *  multiWeights->data[detInd2]->data[sftNum2]*/
	curlyGSqr += SQUARE( GalphaAmp );
	ccSign *= -1;
      }
      baseCCSign *= -1;
    }
  }
  if (curlyGSqr == 0.0)
    {
      *evSquared = 0.0;
      *ccStat = 0.0;
    }
  else
    {
      *evSquared = 8 * SQUARE(multiWeights->Sinv_Tsft) * curlyGSqr;
      *ccStat = 4 * multiWeights->Sinv_Tsft * nume / sqrt(*evSquared);
    }
  return XLAL_SUCCESS;
}

/** Calculate multi-bin cross-correlation statistic using resampling */
/* This assumes rectangular or nearly-rectangular windowing */
int XLALCalculatePulsarCrossCorrStatisticResamp
(
 REAL8Vector                             *restrict ccStatVector,           /**< [out] vector cross-correlation statistic rho */
 REAL8Vector                             *restrict evSquaredVector,        /**< [out] vector (E[rho]/h0^2) */
 REAL8Vector                             *restrict numeEquivAve,           /**< [out] vector for intermediate average statistic */
 REAL8Vector                             *restrict numeEquivCirc,          /**< [out] vector for intermediate circular statistic */
 const REAL8Vector                       *restrict resampCurlyGAmp,        /**< [in] amplitude of curly G for each L_Y_K_X (alpha-indexed) */
 const MultiResampSFTPairMultiIndexList  *restrict resampMultiPairs,       /**< [in] resamp multi list of SFT pairs */
 const MultiNoiseWeights                 *restrict multiWeights,           /**< [in] normalization factor S^-1 & weights for each SFT */
 const PulsarDopplerParams               *restrict binaryTemplateSpacings, /**< [in] Set of spacings for search */
 const PulsarDopplerParams               *restrict dopplerpos,             /**< [in] Doppler point to search */
 const MultiCOMPLEX8TimeSeries           *restrict multiTimeSeries_SRC_a,  /**< [in] resampled time series A */
 const MultiCOMPLEX8TimeSeries           *restrict multiTimeSeries_SRC_b,  /**< [in] resampled time series B */
 ResampCrossCorrWorkspace                *restrict ws,                     /**< [in/out] first workspace */
 COMPLEX8                                *restrict ws1KFaX_k,              /**< [in/out] holder for detector 1 Fa */
 COMPLEX8                                *restrict ws1KFbX_k,              /**< [in/out] holder for detector 1 Fb */
 COMPLEX8                                *restrict ws2LFaX_k,              /**< [in/out] holder for detector 2 Fa */
 COMPLEX8                                *restrict ws2LFbX_k               /**< [in/out] holder for detector 2 Fb */
 )
{
  XLAL_CHECK ( ws != NULL, XLAL_EINVAL );
  XLAL_CHECK ( ws1KFaX_k != NULL, XLAL_EINVAL );
  XLAL_CHECK ( ws1KFbX_k != NULL, XLAL_EINVAL );
  XLAL_CHECK ( ws2LFaX_k != NULL, XLAL_EINVAL );
  XLAL_CHECK ( ws2LFbX_k != NULL, XLAL_EINVAL );

  /* Initialize nume holders for resamp */
  REAL8 curlyEquivGSqrSum = 0;
  for (UINT4 j = 0; j < ws->numFreqBinsOut; j++){
    numeEquivAve->data[j] = 0;
    numeEquivCirc->data[j]= 0;
    ccStatVector->data[j] = 0;
    evSquaredVector->data[j] = 0;
  }

  const REAL8 dt_SRC = multiTimeSeries_SRC_b->data[0]->deltaT;
  const REAL8 SRCsampPerTcoh = resampMultiPairs->Tshort/dt_SRC;

  /* MAIN LOOP (RESAMPLING) */
  for (UINT4 detX=0; detX < resampMultiPairs->length; detX++){
    for (UINT4 sftK=0; sftK <  resampMultiPairs->data[detX].length; sftK++){
        if ( ( XLALComputeFaFb_CrossCorrResamp(ws, ws1KFaX_k, ws1KFbX_k, resampMultiPairs, multiTimeSeries_SRC_a, multiTimeSeries_SRC_b, dopplerpos, binaryTemplateSpacings, SRCsampPerTcoh, detX, sftK, 0, FALSE) ) != XLAL_SUCCESS) {
          LogPrintf ( LOG_CRITICAL, "%s: XLALComputeFaFb_CrossCorrResamp() failed with errno=%d\n", __func__, xlalErrno );
          XLAL_ERROR( XLAL_EFUNC );
        }
        for (UINT4 detY=0; detY < resampMultiPairs->data[detX].data[sftK].length; detY++){
          if ( ( XLALComputeFaFb_CrossCorrResamp(ws, ws2LFaX_k, ws2LFbX_k, resampMultiPairs, multiTimeSeries_SRC_a, multiTimeSeries_SRC_b, dopplerpos, binaryTemplateSpacings, SRCsampPerTcoh, detX, sftK, detY, TRUE) ) != XLAL_SUCCESS) {
            LogPrintf ( LOG_CRITICAL, "%s: XLALComputeFaFb_CrossCorrResamp() failed with errno=%d\n", __func__, xlalErrno );
            XLAL_ERROR( XLAL_EFUNC );
          }
          /* -...-...-...- NEW CORE -...-...-...- */
          for (UINT8 j = 0; j < ws->numFreqBinsOut; j++){
            numeEquivAve->data[j]  += creal(0.1 * (  conj(ws1KFaX_k[j]) * ws2LFaX_k[j] + conj(ws1KFbX_k[j]) * ws2LFbX_k[j]  ) );
            /* Can use for circular polarization:
             * numeEquivCirc->data[j] += creal(0.1 * (  conj(ws1KFaX_k[j]) * ws2LFbX_k[j] - conj(ws1KFbX_k[j]) * ws2LFaX_k[j]  ) ); */
          }
          /* -...-...-...- END NEW CORE -...-...-...- */
        } /* detY */
    } /* sftK */
  } /* detX */ /* end main loop (NOW WITH RESAMPLING) */

  /* ENDING RESAMPLING SECTION */

  /* Normalization factors moved outside of core: */
  for (UINT4 alpha=0; alpha < resampMultiPairs->allPairCount; alpha++) {
      REAL8 GalphaResampAmp = resampCurlyGAmp->data[alpha];
      curlyEquivGSqrSum += SQUARE(GalphaResampAmp);
  }
  /* Because FstatInput in the outer for-loop read in and normalized
   * sfts according to their original SFT length, and that determines
   * the scale of Fa, Fb, and therefore numeEquivAve (according to
   * to Wiener-Kinchine, E[|data|^2] = Tsft * Sn / 2 for spectral
   * noise Sn, and compare XLALNormalizeMultiSFTVect and
   * XLALNormalizeSFT we need to access the original
   * multi-weights S inverse Tsft for normalize numeEquivAve. However,
   * we use the modified multi-weights S inverse Tsft for curlyEquivGSqrSum,
   * because that was created with Tshort and the modified S inverse Tsft.
   *   Note that 
   *     multiWeights->Sinv_Tsft =...
   *       sum_dets(sum_sfts( length_sft/(running median_sft ) ) ),
   * ergo 1/Sinv_Tsft is the harmonic mean of running median noise.
   * Therefore this normalization should be valid even when used with
   * detectors with different noise levels.
   */
  const REAL8 originalMultiWeightsSinvTsft = multiWeights->Sinv_Tsft * (resampMultiPairs->Tsft / resampMultiPairs->Tshort);
  for (UINT8 j = 0; j < ws->numFreqBinsOut; j++){ 
    evSquaredVector->data[j] = 8 * SQUARE(multiWeights->Sinv_Tsft) * curlyEquivGSqrSum;
    ccStatVector->data[j] = 4 * originalMultiWeightsSinvTsft * numeEquivAve->data[j] / sqrt(evSquaredVector->data[j]);
  }
  return XLAL_SUCCESS;
} // end XLALCalculatePulsarCrossCorrStatisticResamp

/** calculate signal phase derivatives wrt Doppler coords, for each SFT */
/* allocates memory as well */
int XLALCalculateCrossCorrPhaseDerivatives
  (
   REAL8VectorSequence        **phaseDerivs, /**< Output: dPhi_K/dlambda_i; i is the "sequence" index, K is the "vector" index */
   const PulsarDopplerParams  *dopplerPoint, /**< Input: pulsar/binary orbit paramaters */
   const EphemerisData                *edat, /**< Input: Earth/Sun ephemeris */
   SFTIndexList                  *indexList, /**< Input: list of SFT indices */
   MultiSSBtimes                *multiTimes, /**< Input: barycentered times of SFTs */
   const DopplerCoordinateSystem  *coordSys  /**< Input: coordinates with which to differentiate */
   )
{
  XLAL_CHECK ( dopplerPoint != NULL, XLAL_EINVAL );
  XLAL_CHECK ( edat != NULL, XLAL_EINVAL );
  XLAL_CHECK ( indexList != NULL, XLAL_EINVAL );
  XLAL_CHECK ( multiTimes != NULL, XLAL_EINVAL );
  XLAL_CHECK ( coordSys != NULL, XLAL_EINVAL );

  const UINT4 numCoords = coordSys->dim;
  const UINT4 numSFTs = indexList->length;

  REAL8VectorSequence *ret = NULL;

  XLAL_CHECK ( ( ret = XLALCreateREAL8VectorSequence ( numCoords, numSFTs ) ) != NULL, XLAL_EFUNC, "XLALCreateREAL8VectorSequence ( %"LAL_UINT4_FORMAT", %"LAL_UINT4_FORMAT" ) failed.", numCoords,numSFTs );

  for ( UINT4 coordNum=0; coordNum < numCoords; coordNum++ ) {
    for ( UINT4 sftNum=0; sftNum < numSFTs; sftNum++ ) {
      UINT4 detInd = indexList->data[sftNum].detInd;
      UINT4 sftInd = indexList->data[sftNum].sftInd;
      SSBtimes *times;
      times = multiTimes->data[detInd];
      REAL8 refTime8 = XLALGPSGetREAL8 ( &(times->refTime) );
      UINT4 numSFTsDet = times->DeltaT->length;
      XLAL_CHECK ( ( sftInd < numSFTsDet ), XLAL_EINVAL, "SFT asked for SFT index off end of list:\n sftNum=%"LAL_UINT4_FORMAT", detInd=%"LAL_UINT4_FORMAT", sftInd=%"LAL_UINT4_FORMAT", numSFTsDet=%"LAL_UINT4_FORMAT"\n", sftNum, detInd, sftInd, numSFTsDet );
      REAL8 tSSB = refTime8 + times->DeltaT->data[sftInd];
      ret->data[coordNum*numSFTs+sftNum] = XLALComputePhaseDerivative ( tSSB, dopplerPoint, (coordSys->coordIDs[coordNum]), edat, NULL, 0 );
      XLAL_CHECK ( xlalErrno == 0, XLAL_EFUNC, "XLALComputePhaseDerivative() failed with xlalErrno = %d\n", xlalErrno );
    }
  }

  (*phaseDerivs) = ret;

  return XLAL_SUCCESS;

}

/** (test function) MODIFIED for Tshort */
/** calculate signal phase derivatives wrt Doppler coords, for each SFT */
/* allocates memory as well */
int XLALCalculateCrossCorrPhaseDerivativesShort
  (
   REAL8VectorSequence        **resampPhaseDerivs, /**< [out] dPhi_K/dlambda_i; i is the "sequence" index, K is the "vector" index */
   const PulsarDopplerParams  *dopplerPoint, /**< [in] pulsar/binary orbit paramaters */
   const EphemerisData                *edat, /**< [in Earth/Sun ephemeris */
   SFTIndexList                  *indexList, /**< [in] list of SFT indices */
   MultiResampSFTPairMultiIndexList  *resampMultiPairs, /**< [in] resamp list of SFT pairs */
   MultiSSBtimes                *multiTimes, /**< [in] barycentered times of SFTs */
   const DopplerCoordinateSystem  *coordSys  /**< [in] coordinates with which to differentiate */
   )
{
  /* Note: the use of tSSBApprox rather than the true tSSB contributes 
   * a very small but observable difference between this function and
   * a correctly-configured call to the normal phase derivatives function. */
  XLAL_CHECK ( dopplerPoint != NULL, XLAL_EINVAL );
  XLAL_CHECK ( edat != NULL, XLAL_EINVAL );
  XLAL_CHECK ( indexList != NULL, XLAL_EINVAL );
  XLAL_CHECK ( resampMultiPairs != NULL, XLAL_EINVAL );
  XLAL_CHECK ( multiTimes != NULL, XLAL_EINVAL );
  XLAL_CHECK ( coordSys != NULL, XLAL_EINVAL );

  /* The actual SFT indices from the original catalog are no longer
   * relevant; in resampling, we use new SFTs. These are not evenly
   * spaced in the detector frame, but for the purpose of taking the
   * values for gamma ave that is needed for the metric, and hence
   * for the returned resampPhaseDerivs, this should average out. */
  const UINT4 numCoords = coordSys->dim;
  const UINT4 numShorts = resampMultiPairs->sftTotalCount;
  const REAL8 Tshort = resampMultiPairs->Tshort; 
  printf("numShorts: %u\n", numShorts);

  REAL8VectorSequence *retNew = NULL;
  XLAL_CHECK ( ( retNew = XLALCreateREAL8VectorSequence ( numCoords, numShorts ) ) != NULL, XLAL_EFUNC, "XLALCreateREAL8VectorSequence ( %"LAL_UINT4_FORMAT", %"LAL_UINT4_FORMAT" ) failed.", numCoords,numShorts );

  for ( UINT4 coordNum=0; coordNum < numCoords; coordNum++ ) {
    for (UINT4 detX=0; detX < resampMultiPairs->length; detX++){
      for (UINT4 k=0; k <  resampMultiPairs->data[detX].length; k++){
        /* Both detectors' resampled time series have the same start time now, 
         * or should be aligned to be so. We just need the nominal time in SSB,
         * so we could probably calculate the start time of the first SFT and
         * add the number of Tshort segments that follow */
        REAL8 startTime8 = XLALGPSGetREAL8( &(multiTimes->data[0]->refTime)) + multiTimes->data[0]->DeltaT->data[0];
        /* Cannot use Tshort * shortNum -- it needs to be the number
         * within THAT detector */
        UINT4 shortSFTIndOwn = resampMultiPairs->data[detX].data[k].sftInd;
        UINT4 shortSFTInd = resampMultiPairs->data[detX].data[k].flatInd;
        REAL8 tSSBApprox = startTime8 + Tshort*shortSFTIndOwn;
        retNew->data[coordNum*numShorts+shortSFTInd] = XLALComputePhaseDerivative ( tSSBApprox, dopplerPoint, (coordSys->coordIDs[coordNum]), edat, NULL, 0 );
        XLAL_CHECK ( xlalErrno == 0, XLAL_EFUNC, "XLALComputePhaseDerivative() failed with xlalErrno = %d\n", xlalErrno );
      }
    }
  }
  (*resampPhaseDerivs) = retNew;

  return XLAL_SUCCESS;
} /* end XLALCalculateCrossCorrPhaseDerivativesShort */

/** calculate phase metric for CW cross-correlation search, as well as vector used for parameter offsets */
/** This calculates the metric defined in (4.7) of Whelan et al 2015 and the parameter offset epsilon_i */
/** (not including the cosi-dependent prefactor) in defined in (4.8) */
/* allocates memory as well */
int XLALCalculateCrossCorrPhaseMetric
  (
   gsl_matrix                        **g_ij, /**< Output: parameter space metric */
   gsl_vector                       **eps_i, /**< Output: parameter offset vector from (4.8) of WSZP15 */
   REAL8                        *sumGammaSq, /**< Output: sum of (Gamma_ave)^2 for normalization and sensitivity */
   const REAL8VectorSequence   *phaseDerivs, /**< Input: dPhi_K/dlambda_i; i is the "sequence" index, K is the "vector" */
   const SFTPairIndexList    *pairIndexList, /**< Input: list of SFT pairs */
   const REAL8Vector             *Gamma_ave, /**< Input: vector of aa+bb values */
   const REAL8Vector            *Gamma_circ, /**< Input: vector of ab-ba values */
   const DopplerCoordinateSystem  *coordSys  /**< Input: coordinate directions for metric */
   )
{
  XLAL_CHECK ( sumGammaSq != NULL, XLAL_EINVAL );
  XLAL_CHECK ( phaseDerivs != NULL, XLAL_EINVAL );
  XLAL_CHECK ( pairIndexList != NULL, XLAL_EINVAL );
  XLAL_CHECK ( Gamma_ave != NULL, XLAL_EINVAL );
  XLAL_CHECK ( Gamma_circ != NULL, XLAL_EINVAL );
  XLAL_CHECK ( coordSys != NULL, XLAL_EINVAL );

  const UINT4 numCoords = coordSys->dim;
  XLAL_CHECK ( ( phaseDerivs->length == numCoords ), XLAL_EINVAL, "Length mismatch: phaseDerivs->length=%"LAL_UINT4_FORMAT", numCoords=%"LAL_UINT4_FORMAT"\n", phaseDerivs->length, numCoords );
  const UINT4 numSFTs = phaseDerivs->vectorLength;
  const UINT4 numPairs = pairIndexList->length;
  XLAL_CHECK ( ( Gamma_ave->length == numPairs ), XLAL_EINVAL, "Length mismatch: Gamma_ave->length=%"LAL_UINT4_FORMAT", numPairs=%"LAL_UINT4_FORMAT"\n", Gamma_ave->length, numPairs );
  XLAL_CHECK ( ( Gamma_circ->length == numPairs ), XLAL_EINVAL, "Length mismatch: Gamma_circ->length=%"LAL_UINT4_FORMAT", numPairs=%"LAL_UINT4_FORMAT"\n", Gamma_circ->length, numPairs );

  /* ---------- prepare output metric ---------- */
  gsl_matrix *ret_g;
  if ( (ret_g = gsl_matrix_calloc ( numCoords, numCoords )) == NULL ) {
    XLALPrintError ("%s: gsl_matrix_calloc(%d, %d) failed.\n\n", __func__, numCoords, numCoords );
    XLAL_ERROR ( XLAL_ENOMEM );

  }
  gsl_vector *ret_e;
  if ( (ret_e = gsl_vector_calloc ( numCoords )) == NULL ) {
    XLALPrintError ("%s: gsl_vector_calloc(%d) failed.\n\n", __func__, numCoords );
    XLAL_ERROR ( XLAL_ENOMEM );
  }

  REAL8Vector *dDeltaPhi_i = NULL;
  REAL8 denom = 0;
  XLAL_CHECK ( ( dDeltaPhi_i = XLALCreateREAL8Vector ( numCoords ) ) != NULL, XLAL_EFUNC, "XLALCreateREAL8Vector ( %"LAL_UINT4_FORMAT" ) failed.", numCoords );
  for ( UINT4 pairNum=0; pairNum < numPairs; pairNum++ ) {
    UINT4 sftNum1 = pairIndexList->data[pairNum].sftNum[0];
    UINT4 sftNum2 = pairIndexList->data[pairNum].sftNum[1];
    REAL8 aveWeight = SQUARE(Gamma_ave->data[pairNum]);
    denom += aveWeight;
    REAL8 circWeight = Gamma_ave->data[pairNum] * Gamma_circ->data[pairNum];
    for ( UINT4 i=0; i < numCoords; i++ ) {
      dDeltaPhi_i->data[i] = phaseDerivs->data[i*numSFTs+sftNum1] - phaseDerivs->data[i*numSFTs+sftNum2];
      REAL8 epsi = gsl_vector_get( ret_e, i );
      epsi += circWeight * dDeltaPhi_i->data[i];
      gsl_vector_set( ret_e, i, epsi );
      for ( UINT4 j=0; j<=i; j++ ) { /* Doing the loop this way ensures dDeltaPhi_i[j] has been set */
	REAL8 gij = gsl_matrix_get( ret_g, i, j );
	gij += aveWeight * dDeltaPhi_i->data[i] * dDeltaPhi_i->data[j];
	gsl_matrix_set ( ret_g, i, j, gij );
      }
    }
  }

  XLALDestroyREAL8Vector ( dDeltaPhi_i );

  for ( UINT4 i=0; i < numCoords; i++ ) {
    REAL8 epsi = gsl_vector_get( ret_e, i );
    epsi /= denom;
    gsl_vector_set( ret_e, i, epsi );
    for ( UINT4 j=0; j<=i; j++ ) { /* Doing the loop the same way as above */
      REAL8 gij = gsl_matrix_get( ret_g, i, j );
      gij /= (2.*denom);
      gsl_matrix_set ( ret_g, i, j, gij );
      if ( i != j ) gsl_matrix_set ( ret_g, j, i, gij );
    }
  }

  (*g_ij) = ret_g;
  (*eps_i) = ret_e;
  (*sumGammaSq) = denom;

  return XLAL_SUCCESS;

}

/** (test function) Redesigning to use Tshort instead */
/** calculate phase metric for CW cross-correlation search, as well as vector used for parameter offsets */
/** This calculates the metric defined in (4.7) of Whelan et al 2015 and the parameter offset epsilon_i */
/** (not including the cosi-dependent prefactor) in defined in (4.8) */
/* allocates memory as well */
int XLALCalculateCrossCorrPhaseMetricShort
  (
   gsl_matrix                        **g_ij, /**< [out] parameter space metric */
   gsl_vector                       **eps_i, /**< [out] parameter offset vector from (4.8) of WSZP15 */
   REAL8                        *sumGammaSq, /**< [out] sum of (Gamma_ave)^2 for normalization and sensitivity */
   const REAL8VectorSequence   *phaseDerivs, /**< [in] dPhi_K/dlambda_i; i is the "sequence" index, K is the "vector" */
   const MultiResampSFTPairMultiIndexList *resampMultiPairs, /**< [in] resamp multi list of SFT pairs */
   const REAL8Vector             *Gamma_ave, /**< [in]: vector of aa+bb values */
   const REAL8Vector            *Gamma_circ, /**< [in] vector of ab-ba values */
   const DopplerCoordinateSystem  *coordSys  /**< [in] coordinate directions for metric */
   )
{
  XLAL_CHECK ( sumGammaSq != NULL, XLAL_EINVAL );
  XLAL_CHECK ( phaseDerivs != NULL, XLAL_EINVAL );
  XLAL_CHECK ( Gamma_ave != NULL, XLAL_EINVAL );
  XLAL_CHECK ( Gamma_circ != NULL, XLAL_EINVAL );
  XLAL_CHECK ( coordSys != NULL, XLAL_EINVAL );

  const UINT4 numCoords = coordSys->dim;
  XLAL_CHECK ( ( phaseDerivs->length == numCoords ), XLAL_EINVAL, "Length mismatch: phaseDerivs->length=%"LAL_UINT4_FORMAT", numCoords=%"LAL_UINT4_FORMAT"\n", phaseDerivs->length, numCoords );
  const UINT4 numSFTs = phaseDerivs->vectorLength;
  const UINT4 numShorts = resampMultiPairs->sftTotalCount;
  XLAL_CHECK ( ( numSFTs == numShorts ), XLAL_EINVAL, "Length mismatch: phaseDerivs->vectorLength=%"LAL_UINT4_FORMAT", resampMultiPairs->sftTotalCount=%"LAL_UINT4_FORMAT"\n", numSFTs, numShorts );
  const UINT4 numShortPairs = resampMultiPairs->allPairCount;
  XLAL_CHECK ( ( Gamma_ave->length == numShortPairs ), XLAL_EINVAL, "Length mismatch: Gamma_ave->length=%"LAL_UINT4_FORMAT", numPairs=%"LAL_UINT4_FORMAT"\n", Gamma_ave->length, numShortPairs );
  XLAL_CHECK ( ( Gamma_circ->length == numShortPairs ), XLAL_EINVAL, "Length mismatch: Gamma_circ->length=%"LAL_UINT4_FORMAT", numPairs=%"LAL_UINT4_FORMAT"\n", Gamma_circ->length, numShortPairs );

  /* ---------- prepare output metric ---------- */
  gsl_matrix *ret_gNew;
  if ( (ret_gNew = gsl_matrix_calloc ( numCoords, numCoords )) == NULL ) {
    XLALPrintError ("%s: gsl_matrix_calloc(%d, %d) failed.\n\n", __func__, numCoords, numCoords );
    XLAL_ERROR ( XLAL_ENOMEM );

  }
  gsl_vector *ret_eNew;
  if ( (ret_eNew = gsl_vector_calloc ( numCoords )) == NULL ) {
    XLALPrintError ("%s: gsl_vector_calloc(%d) failed.\n\n", __func__, numCoords );
    XLAL_ERROR ( XLAL_ENOMEM );
  }

  REAL8Vector *dDeltaPhiNew_i = NULL;
  XLAL_CHECK ( ( dDeltaPhiNew_i = XLALCreateREAL8Vector ( numCoords ) ) != NULL, XLAL_EFUNC, "XLALCreateREAL8Vector ( %"LAL_UINT4_FORMAT" ) failed.", numCoords );

  /* new style, with proper sftIndices Tshort */
  UINT4 pairMetric = 0;
  UINT4 sftInd1 = 0;
  UINT4 sftInd2 = 0;
  REAL8 aveNewWeight = 0.0;
  REAL8 denomNew = 0.0;
  REAL8 circNewWeight = 0.0;
  for (UINT4 detX=0; detX < resampMultiPairs->length; detX++){
    for (UINT4 k=0; k <  resampMultiPairs->data[detX].length; k++){
      for (UINT4 detY=0; detY < resampMultiPairs->data[detX].data[k].length; detY++){
        for (UINT4 l=0; l < resampMultiPairs->data[detX].data[k].data[detY].length; l++){
          sftInd1 = resampMultiPairs->data[detX].data[k].flatInd;
          sftInd2 = resampMultiPairs->data[detX].data[k].data[detY].data[l].flatInd;
          pairMetric = resampMultiPairs->data[detX].data[k].data[detY].data[l].pairInd;
          aveNewWeight = SQUARE(Gamma_ave->data[pairMetric]);
          denomNew += aveNewWeight;
          circNewWeight = Gamma_ave->data[pairMetric] * Gamma_circ->data[pairMetric];
          for ( UINT4 i=0; i < numCoords; i++ ) {
            dDeltaPhiNew_i->data[i] = phaseDerivs->data[i*numShorts+sftInd1] - phaseDerivs->data[i*numShorts+sftInd2];      
            REAL8 epsNewi = gsl_vector_get( ret_eNew, i);
            epsNewi += circNewWeight * dDeltaPhiNew_i->data[i];
            gsl_vector_set( ret_eNew, i, epsNewi);
            for ( UINT4 j=0; j<=i; j++ ) { /* Doing the loop this way ensures dDeltaPhi_i[j] has been set */
              REAL8 gijNew = gsl_matrix_get ( ret_gNew, i, j);
              gijNew += aveNewWeight * dDeltaPhiNew_i->data[i] * dDeltaPhiNew_i->data[j];
              gsl_matrix_set ( ret_gNew, i, j, gijNew);
            } // end j loop
          } // end i loop over coords
        }
      }
    }
  }

  XLALDestroyREAL8Vector ( dDeltaPhiNew_i );

  for ( UINT4 i=0; i < numCoords; i++ ) {
    REAL8 epsNewi = gsl_vector_get( ret_eNew, i );
    epsNewi /= denomNew;
    gsl_vector_set( ret_eNew, i, epsNewi );
    for ( UINT4 j=0; j<=i; j++ ) { /* Doing the loop the same way as above */
      REAL8 gijNew = gsl_matrix_get( ret_gNew, i, j );
      gijNew /= (2.*denomNew);
      gsl_matrix_set ( ret_gNew, i, j, gijNew );
      if ( i != j ) gsl_matrix_set ( ret_gNew, j, i, gijNew );
    }
  }

  (*g_ij) = ret_gNew;
  (*eps_i) = ret_eNew;
  (*sumGammaSq) = denomNew;

  return XLAL_SUCCESS;

} // end XLALCalculateCrossCorrPhaseMetricShort

/*calculate metric diagonal components, also include the estimation of sensitivity E[rho]/(h_0)^2*/
int XLALCalculateLMXBCrossCorrDiagMetric
  (
   REAL8                      *hSens, /* Output: sensitivity*/
   REAL8                       *g_ff, /* Output: Diagonal frequency metric element */
   REAL8                       *g_aa, /* Output: Diagonal binary projected semimajor axis metric element*/
   REAL8                       *g_TT, /* Output: Diagonal reference time metric element*/
   REAL8                       *g_pp, /* Output: Diagonal orbital period metric element */
   REAL8             *weightedMuTAve, /* output: weighred T mean*/
   PulsarDopplerParams DopplerParams, /*  Input: pulsar/binary orbit paramaters*/
   REAL8Vector              *G_alpha, /*  Input: vector of curlyGunshifted values */
   SFTPairIndexList   *pairIndexList, /*  Input: list of SFT pairs */
   SFTIndexList           *indexList, /*  Input: list of SFTs */
   MultiSFTVector              *sfts, /*  Input: set of per-detector SFT vectors */
   MultiNoiseWeights   *multiWeights  /*  Input: Input: nomalizeation factor S^-1 & weights for each SFT*/
   )

{
  UINT4 sftNum1 = 0;
  UINT4 sftNum2 = 0;
  REAL8 TDiff = 0;
  REAL8 TMean = 0;
  REAL8 denom = 0;
  REAL8 TSquaWeightedAve = 0;
  REAL8 SinSquaWeightedAve = 0;
  REAL8 muT = 0;
  REAL8 muTSqr = 0;
  REAL8 muTAve = 0;
  REAL8 muTAveSqr = 0;
  REAL8 muTSqrAve = 0;
  REAL8 sinSquare = 0;
  REAL8 tSquare = 0;
  REAL8 rhosum = 0;
  LIGOTimeGPS *T1 = NULL;
  LIGOTimeGPS *T2 = NULL;
  UINT4 numalpha = G_alpha->length;


  for (UINT4 j = 0; j < numalpha; j++) {
    sftNum1 = pairIndexList->data[j].sftNum[0];
    sftNum2 = pairIndexList->data[j].sftNum[1];
    UINT4 detInd1 = indexList->data[sftNum1].detInd;
    UINT4 detInd2 = indexList->data[sftNum2].detInd;
    UINT4 sftInd1 = indexList->data[sftNum1].sftInd;
    UINT4 sftInd2 = indexList->data[sftNum2].sftInd;
    T1 = &(sfts->data[detInd1]->data[sftInd1].epoch);
    T2 = &(sfts->data[detInd2]->data[sftInd2].epoch);
    TDiff = XLALGPSDiff(T1, T2);
    TMean = 0.5 * (XLALGPSGetREAL8(T1) + XLALGPSGetREAL8(T2));
    REAL8 sqrG_alpha = SQUARE(G_alpha->data[j]); /*(curlyG_{\alpha})^2*/
    muT +=  sqrG_alpha * TMean; /*(curlyG_\alpha)^2 * \bar{t}_{\alpha}*/
    muTSqr += sqrG_alpha * SQUARE(TMean); /*(curlyG_\alpha)^2 * (\bar{t}_{\alpha})^2*/
    sinSquare += sqrG_alpha * SQUARE(sin(LAL_PI * TDiff/(DopplerParams.period))); /*(G_{\alpha})^2*(sin(\pi*T/T_orbit))^2*/
    tSquare += sqrG_alpha * SQUARE(TDiff); /*(\curlyg_{\alpha}*)^2*T^2*/
    denom += sqrG_alpha; /*calculate the denominator*/
    rhosum += 2*sqrG_alpha;
  }

  muTAve = muT / denom;
  muTAveSqr = SQUARE(muTAve);
  muTSqrAve = muTSqr / denom;
  REAL8 sigmaTSqr = muTSqrAve - muTAveSqr;
  TSquaWeightedAve = tSquare / denom;
  SinSquaWeightedAve = sinSquare / denom;
  *hSens = 4 * SQUARE(multiWeights->Sinv_Tsft) * rhosum;
  *g_ff = TSquaWeightedAve * 2 * SQUARE(LAL_PI);
  *g_aa = SinSquaWeightedAve * SQUARE(2. * LAL_PI * DopplerParams.fkdot[0]);
  *g_TT = SinSquaWeightedAve * SQUARE(SQUARE(2. * LAL_PI) * (DopplerParams.fkdot[0]) * (DopplerParams.asini) / (DopplerParams.period));
  *g_pp = SinSquaWeightedAve * sigmaTSqr * 16 * QUAD(LAL_PI) * SQUARE(DopplerParams.fkdot[0]) * SQUARE(DopplerParams.asini) / (QUAD(DopplerParams.period));
  *weightedMuTAve = muTAve;

  return XLAL_SUCCESS;

}

/** MODIFIED for Tshort:
 * calculate metric diagonal components, also include the estimation of sensitivity E[rho]/(h_0)^2*/
int XLALCalculateLMXBCrossCorrDiagMetricShort
  (
   REAL8                                           * hSens,                    /**< [out] sensitivity */
   REAL8                                           * g_ff,                     /**< [out] Diagonal frequency metric element */
   REAL8                                           * g_aa,                     /**< [out] Diagonal binary projected semimajor axis metric element*/
   REAL8                                           * g_TT,                     /**< [out] Diagonal reference time metric element */
   REAL8                                           * g_pp,                     /**< [out] Diagonal orbital period metric element */
   const PulsarDopplerParams                         DopplerParams,            /**< [in] pulsar/binary orbit paramaters */
   const REAL8Vector                       *restrict G_alpha,                  /**< [in] vector of curlyGunshifted values */
   const MultiResampSFTPairMultiIndexList  *restrict resampMultiPairIndexList, /**< [in] resamp list of SFT pairs */
   const MultiLIGOTimeGPSVector            *restrict timestamps,               /**< [in] timestamps for resampling */
   const MultiNoiseWeights                 *restrict multiWeights              /**< [in] normalization factor S^-1 & weights for each SFT */
   )

{
  UINT4 sftNum1 = 0;
  UINT4 sftNum2 = 0;
  REAL8 TDiff = 0;
  REAL8 TMean = 0;
  REAL8 denom = 0;
  REAL8 TSquaWeightedAve = 0;
  REAL8 SinSquaWeightedAve = 0;
  REAL8 muT = 0;
  REAL8 muTSqr = 0;
  REAL8 muTAve = 0;
  REAL8 muTAveSqr = 0;
  REAL8 muTSqrAve = 0;
  REAL8 sinSquare = 0;
  REAL8 tSquare = 0;
  REAL8 rhosum = 0;
  LIGOTimeGPS *T1 = NULL;
  LIGOTimeGPS *T2 = NULL;

  UINT4 numPairs = resampMultiPairIndexList->allPairCount;
  SFTPairIndexList *pairIndexList = resampMultiPairIndexList->pairIndexList;
  SFTIndexList *indexList = resampMultiPairIndexList->indexList;

  for (UINT4 j = 0; j < numPairs; j++) {
    sftNum1 = pairIndexList->data[j].sftNum[0];
    sftNum2 = pairIndexList->data[j].sftNum[1];
    UINT4 detInd1 = indexList->data[sftNum1].detInd;
    UINT4 detInd2 = indexList->data[sftNum2].detInd;
    UINT4 sftInd1 = indexList->data[sftNum1].sftInd;
    UINT4 sftInd2 = indexList->data[sftNum2].sftInd;
    T1 = &(timestamps->data[detInd1]->data[sftInd1]) ; 
    T2 = &(timestamps->data[detInd2]->data[sftInd2]) ; 
    TDiff = XLALGPSDiff(T1, T2);
    TMean = 0.5 * (XLALGPSGetREAL8(T1) + XLALGPSGetREAL8(T2));
    REAL8 sqrG_alpha = SQUARE(G_alpha->data[j]); /*(curlyG_{\alpha})^2*/
    muT +=  sqrG_alpha * TMean; /*(curlyG_\alpha)^2 * \bar{t}_{\alpha}*/
    muTSqr += sqrG_alpha * SQUARE(TMean); /*(curlyG_\alpha)^2 * (\bar{t}_{\alpha})^2*/
    sinSquare += sqrG_alpha * SQUARE(sin(LAL_PI * TDiff/(DopplerParams.period))); /*(G_{\alpha})^2*(sin(\pi*T/T_orbit))^2*/
    tSquare += sqrG_alpha * SQUARE(TDiff); /*(\curlyg_{\alpha}*)^2*T^2*/
    denom += sqrG_alpha; /*calculate the denominator*/
    rhosum += 2*sqrG_alpha;
  }
  muTAve = muT / denom;
  muTAveSqr = SQUARE(muTAve);
  muTSqrAve = muTSqr / denom;
  REAL8 sigmaTSqr = muTSqrAve - muTAveSqr;
  TSquaWeightedAve = tSquare / denom;
  SinSquaWeightedAve = sinSquare / denom;
  *hSens = 4 * SQUARE(multiWeights->Sinv_Tsft) * rhosum;
  *g_ff = TSquaWeightedAve * 2 * SQUARE(LAL_PI);
  *g_aa = SinSquaWeightedAve * SQUARE(2. * LAL_PI * DopplerParams.fkdot[0]);
  *g_TT = SinSquaWeightedAve * SQUARE(SQUARE(2. * LAL_PI) * (DopplerParams.fkdot[0]) * (DopplerParams.asini) / (DopplerParams.period));
  *g_pp = SinSquaWeightedAve * sigmaTSqr * 16 * QUAD(LAL_PI) * SQUARE(DopplerParams.fkdot[0]) * SQUARE(DopplerParams.asini) / (QUAD(DopplerParams.period));

  return XLAL_SUCCESS;

} // end XLALCalculateLMXBCrossCorrDiagMetricShort

///** (possible future function) Wrapper for Bessel Orbital Space stepper in the manner of V. Dergachev's loosely-coherent search */
//int 
//XLALBesselCrossCorrOrbitalSpaceStep( COMPLEX8             * xTildeOut,              /**< output Fa or Fb (FFT of resampled time series) */
//                                     COMPLEX8             * xTildeIn,               /**< input Fa or Fb (FFT of resampled time series) */
//                                     PulsarDopplerParams  * dopplerPoint,           /**< Input: Doppler point to search */
//                                     PulsarDopplerParams  * binaryTemplateSpacings, /**< Input: Set of spacings for search */
//                                     INT4                   aStep,                  /**< how many steps in a sin i */
//                                     INT4                   pStep,                  /**< how many steps in period */
//                                     INT4                   tStep,                  /**< how many steps in time of ascension */
//                                     INT4                   nFFT,                   /**< number of samples in the FFT */
//                                     UINT4                  nTermsMax               /**< Maximum number of terms in convolution */
//)
//{
//  /* Experimental shift */
//  /* Example call to this function from Compute Fa Fb: */
//  //XLALBesselCrossCorrOrbitalSpaceStep(ws->FaX_k, ws->FaX_k, dopplerpos, binaryTemplateSpacings, 0, 0, 0, ws->numFreqBinsOut, 10);
//  //XLALBesselCrossCorrOrbitalSpaceStep(ws->FbX_k, ws->FbX_k, dopplerpos, binaryTemplateSpacings, 0, 0, 0, ws->numFreqBinsOut, 10);
//
//  XLAL_CHECK ( xTildeOut != NULL, XLAL_EINVAL );
//  XLAL_CHECK ( xTildeIn != NULL, XLAL_EINVAL );
//  XLAL_CHECK ( dopplerPoint != NULL, XLAL_EINVAL );
//  if (pStep > (INT4)nTermsMax){
//    fprintf(stdout, "pStep seems suspiciously large\n");
//  }
//  if (tStep > (INT4)nTermsMax){
//    fprintf(stdout, "tStep seems suspiciously large\n");
//  }
//  //fprintf(stdout, "aStep, pStep, tStep, nFFT, nTermsMax: %i, %i, %i, %i, %u\n", aStep, pStep, tStep, nFFT, nTermsMax);
//
//  COMPLEX8 *besselFX_k;
//  XLAL_CHECK ( (besselFX_k = fftw_malloc ( nFFT * sizeof(COMPLEX8) )) != NULL, XLAL_ENOMEM );
//  memset ( besselFX_k, 0, nFFT * sizeof(besselFX_k[0]) );
//  //memcpy ( besselFX_k, xTildeIn, nFFT);
//
//  /* Act upon the input */
//  /* Simple test */
//  //for ( INT4 k = 0; k < nFFT; k++ )
//  //{
//  //  COMPLEX8 normX_k = 10.0 ;
//  //  besselFX_k[k] = normX_k * xTildeIn[k];
//  //  //if (k == nFFT-1){
//  //  //    fprintf(stdout, "Made it!\n");
//  //  //}
//  //} // for k < nFFT
//
//  //REAL8 asiniPos = dopplerPoint->asini;
//  REAL8 minFreq = dopplerPoint->fkdot[0];
//  REAL8 nominalPeriod = dopplerPoint->period;
//  REAL8 nominalTp = XLALGPSGetREAL8(&dopplerPoint->tp);
//  REAL8 asiniSpace = binaryTemplateSpacings->asini;
//  //REAL8 tpSpace = XLALGPSGetREAL8(&binaryTemplateSpacings->tp);
//  //REAL8 periodSpace = binaryTemplateSpacings->period;
//  REAL8 asiniShift = aStep * asiniSpace;
//  //REAL8 tpShift = tStep * tpSpace;
//  //REAL8 periodShift = pStep * periodSpace;  
//  REAL8 argBessel = 2 * LAL_PI * minFreq * asiniShift;
//  /* Number of sidebands is 1 + 2*(2*pi*f0*asini)*/
//  INT4 oneSidedNsidebands = (INT4)floor(argBessel);
//  if ( (UINT4)oneSidedNsidebands > nTermsMax){
//    fprintf(stdout, "Warning! large number of sidebands! Is this the right part of space?\n");
//  }
//  INT4 nSidebands = 1 + 2 * oneSidedNsidebands;
//  /* Orbital phase shift */
//  REAL8 orbitalPhaseShift = fmod( (2 * LAL_PI)/(nominalPeriod) * (nominalTp), 2 * LAL_PI );
//
//  /* For reference */
//          //REAL8 f_k = (offset_bins2L + (UINT4)floor(k*DedecimateFFT2L)) * dFreqFFT2L+0.0;
//          //REAL8 cycles = - (f_k) * (dtauX2L - dtauX1K); /* Back to usual FFT trick */
//          //REAL4 sinphase, cosphase;
//          //XLALSinCos2PiLUT ( &sinphase, &cosphase, cycles );
//          //COMPLEX8 normX_k = dt_SRC * crectf ( cosphase, sinphase );
//
//  INT4 sidebandIndex;
//  REAL8 cycles;
//  REAL8Vector *besselCoeffs = NULL;
//  XLAL_CHECK ( (besselCoeffs = XLALCreateREAL8Vector ( nSidebands )) != NULL, XLAL_EFUNC );
//  for ( INT4 j = 0; j < nSidebands; j++){
//    /* Wrap argBessel in J_s to get true filter coefficients */
//    /* Note that J_s has s = (j - oneSidedNsidebands), it this notation */
//    /* TECHNICALLY, this should depend on the frequency too */
//    besselCoeffs->data[j] = argBessel;
//  }
//  COMPLEX8Vector *shiftsInPhase = NULL;
//  REAL4 sinphase, cosphase;
//  XLAL_CHECK ( (shiftsInPhase = XLALCreateCOMPLEX8Vector ( nSidebands )) != NULL, XLAL_EFUNC );
//  for ( INT4 j = 0; j < nSidebands; j++){
//    sidebandIndex = j - oneSidedNsidebands;
//    cycles = sidebandIndex * orbitalPhaseShift;
//    XLALSinCos2PiLUT ( &sinphase, &cosphase, cycles );
//    shiftsInPhase->data[j] = crectf( cosphase, sinphase );
//  }
//
//   
//  //fprintf(stdout, "nSidebands, tpShift, periodShift: %u, %f, %f\n", nSidebands, tpShift, periodShift);
//  //fprintf(stdout, "nSidebands, argBessel, orbitalPhaseShift: %u, %f, %f\n", nSidebands, argBessel, orbitalPhaseShift);
//  //fprintf(stdout, "nSidebands, besselCoeffs, cycles: %u, %f, %f\n", nSidebands, besselCoeffs->data[0], cycles);
//  fprintf(stdout, "nSidebands, besselCoeffs, shiftsInPhase: %u, %f, %f\n", nSidebands, besselCoeffs->data[0], cimag(shiftsInPhase->data[0]));
//  //fprintf(stdout, "asini, spacing: %f, %f\n", asiniPos, asiniSpace);
//  /* More complex test */
//  //COMPLEX8 partialSumK = 0;
//  INT4 edgeOfSidebands = 0;
//  //INT4 localNsidebands = 0;
//  INT4 diffInEdge = 0;
//  for ( INT4 k = 0; k < nFFT; k++ )
//  {
//    //COMPLEX8 normX_k = 1.0 ;
//    //besselFX_k[k] = normX_k * xTildeIn[k];
//    /* start bessel application */
//    edgeOfSidebands = MYMIN(k, oneSidedNsidebands);
//    diffInEdge = oneSidedNsidebands - edgeOfSidebands;
//    //localNsidebands = MYMIN(2*edgeOfSidebands+1, nSidebands);
//    for ( INT4 j = 0; j < nSidebands; j++){
//        sidebandIndex = j - edgeOfSidebands;
//        /* May need to take frequency spacing into account very soon*/
//        //fprintf(stdout, "diffInEdge!: %i\n", diffInEdge);
//        besselFX_k[k] += besselCoeffs->data[j+diffInEdge] * xTildeIn[k+sidebandIndex] * shiftsInPhase->data[j+diffInEdge];
//    }
//    /* end bessel application */
//    //if (k == nFFT-1){
//    //    fprintf(stdout, "Made it!\n");
//    //}
//  } // for k < nFFT
//  
//  
//  /* for testing */
//  //xTildeOut = xTildeIn;
//  memcpy (xTildeOut, besselFX_k, nFFT * sizeof(COMPLEX8));
//  //for ( INT4 l = 0; l < nFFT; l++ )
//  //{
//  //    xTildeOut[l] = besselFX_k[l];
//  //}
//  free (besselFX_k);
//  XLALDestroyREAL8Vector(besselCoeffs);
//  XLALDestroyCOMPLEX8Vector(shiftsInPhase);
//  return XLAL_SUCCESS;
//} // XLALBesselCrossCorrOrbitalSpaceStep

/** Function to extract timestamps with tShort */
LIGOTimeGPSVector *
XLALExtractTimestampsFromSFTsShort (
                        REAL8TimeSeries  ** sciFlag,       /**< [out] science or not flag */
                        const SFTVector   * sfts ,         /**< [in] input SFT-vector  */
                        REAL8               tShort,        /**< [in] new time baseline, Tshort */
                        UINT4               numShortPerDet /**< [in] number of Tshort per detector */
)
{

  REAL8TimeSeries *retFlag = NULL;
  /* check input consistency */
  if ( !sfts ) {
    XLALPrintError ("%s: invalid NULL input 'sfts'\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  UINT4 numSFTsNew = numShortPerDet;
  /* create output vector */
  LIGOTimeGPSVector *ret = NULL;
  if ( ( ret = XLALCreateTimestampVector ( numSFTsNew )) == NULL ) {
    XLALPrintError ("%s: XLALCreateTimestampVector(%d) failed.\n", __func__, numSFTsNew );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }
  ret->deltaT = tShort;
  ret->length = numSFTsNew;

  UINT4 i;
  LIGOTimeGPS epochHolder = {0,0};
  LIGOTimeGPS naiveEpoch = sfts->data[0].epoch;
  REAL8 naiveEpochReal = XLALGPSGetREAL8(&( naiveEpoch ));
  /* Wanting a consistent baseline, calculate this externally */
  for ( i=0; i < numShortPerDet; i ++ ){
    epochHolder.gpsSeconds = 0;
    epochHolder.gpsNanoSeconds = 0;
    XLALGPSAdd(&(epochHolder), naiveEpochReal);
    XLALGPSAdd(&(epochHolder), i * tShort );
    ret->data[i] = epochHolder;
  }

  retFlag = XLALCrossCorrGapFinderResampAlt(ret, sfts);
  /* done: return Ts-vector */
  (*sciFlag) = retFlag;
  return ret;

} /* XLALExtractTimestampsFromSFTsShort() */

/** Modify timestamps from one detector with tShort */
LIGOTimeGPSVector *
XLALModifyTimestampsFromSFTsShort ( 
                        REAL8TimeSeries         **         sciFlag,       /**< [out] science or not flag */
                        const LIGOTimeGPSVector *restrict  Times,         /**< [in] input SFT-vector */
                        const REAL8                        tShort,        /**< [in] new time baseline, Tshort */
                        const UINT4                        numShortPerDet /**< [in] number of Tshort per detector */
)
{

  REAL8TimeSeries *retFlag = NULL;
  /* check input consistency */
  if ( !Times ) {
    XLALPrintError ("%s: invalid NULL input 'Times'\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  UINT4 numSFTsNew = numShortPerDet;
  /* create output vector */
  LIGOTimeGPSVector *ret = NULL;
  if ( ( ret = XLALCreateTimestampVector ( numSFTsNew )) == NULL ) {
    XLALPrintError ("%s: XLALCreateTimestampVector(%d) failed.\n", __func__, numSFTsNew );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }
  ret->deltaT = tShort;
  ret->length = numSFTsNew;

  UINT4 i;
  LIGOTimeGPS epochHolder = {0,0};
  LIGOTimeGPS naiveEpoch = Times->data[0];
  REAL8 naiveEpochReal = XLALGPSGetREAL8(&( naiveEpoch ));
  /* Wanting a consistent baseline, calculate this externally */
  for ( i=0; i < numShortPerDet; i ++ ){
    epochHolder.gpsSeconds = 0;
    epochHolder.gpsNanoSeconds = 0;
    XLALGPSAdd(&(epochHolder), naiveEpochReal);
    XLALGPSAdd(&(epochHolder), i * tShort ); 
    ret->data[i] = epochHolder;
  }
  retFlag = XLALCrossCorrGapFinderResamp(ret, Times);
  /* done: return Ts-vector */
  (*sciFlag) = retFlag;
  return ret;

} /* XLALModifyTimestampsFromSFTsShort() */

/**
 * Given a multi-SFT vector, return a MultiLIGOTimeGPSVector holding the
 * modified SFT timestamps (production function using tShort)
 */
MultiLIGOTimeGPSVector *
XLALModifyMultiTimestampsFromSFTs ( 
                        MultiREAL8TimeSeries         **         scienceFlagVect, /**< [out] science flag vector */
                        const MultiLIGOTimeGPSVector  *restrict multiTimes,      /**< [in] multiTimes, standard */
                        const REAL8                             tShort,          /**< [in] new time baseline, Tshort */
                        const UINT4                             numShortPerDet   /**< [in] number of Tshort per detector */
)
{

  MultiREAL8TimeSeries *retFlagVect = NULL;
  /* check input consistency */
  if ( !multiTimes || multiTimes->length == 0 ) {
    XLALPrintError ("%s: illegal NULL or empty input 'multiTimes'.\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }
  UINT4 numIFOs = multiTimes->length;

  /* create output vector */
  MultiLIGOTimeGPSVector *ret = NULL;
  if ( (ret = XLALCalloc ( 1, sizeof(*ret) )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCalloc ( 1, %zu ).\n", __func__, sizeof(*ret));
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  if ( (ret->data = XLALCalloc ( numIFOs, sizeof(*ret->data) )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCalloc ( %d, %zu ).\n", __func__, numIFOs, sizeof(ret->data[0]) );
    XLALFree (ret);
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }
  ret->length = numIFOs;

  /* new: duplicate above for the science flag vector */
  if ( (retFlagVect = XLALCalloc ( 1, sizeof(*retFlagVect) )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCalloc ( 1, %zu ).\n", __func__, sizeof(*retFlagVect));
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  if ( (retFlagVect->data = XLALCalloc ( numIFOs, sizeof(*retFlagVect->data) )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCalloc ( %d, %zu ).\n", __func__, numIFOs, sizeof(retFlagVect->data[0]) );
    XLALFree (retFlagVect);
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }
  retFlagVect->length = numIFOs;

  /* now extract timestamps vector from each SFT-vector */
  UINT4 X;
  for ( X=0; X < numIFOs; X ++ )
    {
      if ( (ret->data[X] = XLALModifyTimestampsFromSFTsShort ( &retFlagVect->data[X], multiTimes->data[X], tShort, numShortPerDet )) == NULL ) {
        XLALPrintError ("%s: XLALExtractTimestampsFromSFTs() failed for X=%d\n", __func__, X );
        XLALDestroyMultiTimestamps ( ret );
        XLAL_ERROR_NULL ( XLAL_EFUNC );
      }

    } /* for X < numIFOs */

  (*scienceFlagVect) = retFlagVect;
  return ret;

} /* XLALModifyMultiTimestampsFromSFTsShort() */

/**
 * Given a multi-SFT vector, return a MultiLIGOTimeGPSVector holding the
 * SFT timestamps (test function using tShort)
 */
MultiLIGOTimeGPSVector *
XLALExtractMultiTimestampsFromSFTsShort (
                        MultiREAL8TimeSeries  ** scienceFlagVect, /**< [out] science flag vector */
                        const MultiSFTVector  *  multiSFTs,       /**< [in] multiSFT vector */
                        REAL8                    tShort,          /**< [in] new time baseline, Tshort */
                        UINT4                    numShortPerDet   /**< [in] number of Tshort per detector */
)
{

  MultiREAL8TimeSeries *retFlagVect = NULL;
  /* check input consistency */
  if ( !multiSFTs || multiSFTs->length == 0 ) {
    XLALPrintError ("%s: illegal NULL or empty input 'multiSFTs'.\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }
  UINT4 numIFOs = multiSFTs->length;

  /* create output vector */
  MultiLIGOTimeGPSVector *ret = NULL;
  if ( (ret = XLALCalloc ( 1, sizeof(*ret) )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCalloc ( 1, %zu ).\n", __func__, sizeof(*ret));
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  if ( (ret->data = XLALCalloc ( numIFOs, sizeof(*ret->data) )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCalloc ( %d, %zu ).\n", __func__, numIFOs, sizeof(ret->data[0]) );
    XLALFree (ret);
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }
  ret->length = numIFOs;

  /* new: duplicate above for the science flag vector */
  if ( (retFlagVect = XLALCalloc ( 1, sizeof(*retFlagVect) )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCalloc ( 1, %zu ).\n", __func__, sizeof(*retFlagVect));
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  if ( (retFlagVect->data = XLALCalloc ( numIFOs, sizeof(*retFlagVect->data) )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCalloc ( %d, %zu ).\n", __func__, numIFOs, sizeof(retFlagVect->data[0]) );
    XLALFree (retFlagVect);
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }
  retFlagVect->length = numIFOs;

  /* now extract timestamps vector from each SFT-vector */
  UINT4 X;
  for ( X=0; X < numIFOs; X ++ )
    {

      if ( (ret->data[X] = XLALExtractTimestampsFromSFTsShort ( &retFlagVect->data[X], multiSFTs->data[X], tShort, numShortPerDet )) == NULL ) {
        XLALPrintError ("%s: XLALExtractTimestampsFromSFTsShort() failed for X=%d\n", __func__, X );
        XLALDestroyMultiTimestamps ( ret );
        XLAL_ERROR_NULL ( XLAL_EFUNC );
      }

    } /* for X < numIFOs */

  (*scienceFlagVect) = retFlagVect;
  return ret;

} /* XLALExtractMultiTimestampsFromSFTsShort() */

/** (test function) fill detector state with tShort, importing various slightly-modified LALPulsar functions for testing */
int
XLALFillDetectorTensorShort (DetectorState *detState,   /**< [out,in]: detector state: fill in detector-tensor */
                        const LALDetector *detector     /**< [in]: which detector */
                        )
{
  const CHAR *prefix;

  if ( !detState || !detector ) {
    xlalErrno = XLAL_EINVAL;
    return -1;
  }

  prefix = detector->frDetector.prefix;

  /* we need to distinguish two cases: space-borne (i.e. LISA) and Earth-based detectors */
  if ( prefix[0] == 'Z' )       /* LISA */
    {
      //if ( XLALprecomputeLISAarms ( detState ) != 0 ) {
      //  XLALPrintError ("\nXLALprecomputeLISAarms() failed !\n\n");
      //  xlalErrno = XLAL_EINVAL;
      //  return -1;
      //}
      //
      //if ( XLALgetLISADetectorTensorLWL ( &(detState->detT), detState->detArms, prefix[1] ) != 0 ) {
      //  XLALPrintError ("\nXLALgetLISADetectorTensorLWL() failed !\n\n");
      //  xlalErrno = XLAL_EINVAL;
      //  return -1;
      //}

    } /* if LISA */
  else
    {
      REAL4 sinG, cosG, sinGcosG, sinGsinG, cosGcosG;
      SymmTensor3 *detT = &(detState->detT);

      XLAL_CHECK( XLALSinCosLUT ( &sinG, &cosG, detState->earthState.gmstRad ) == XLAL_SUCCESS, XLAL_EFUNC );
      sinGsinG = sinG * sinG;
      sinGcosG = sinG * cosG;
      cosGcosG = cosG * cosG;

      /*
      printf("GMST = %fdeg; cosG = %f, sinG= %f\n",
             LAL_180_PI * atan2(sinG,cosG), cosG, sinG);
      */

      detT->d11 = detector->response[0][0] * cosGcosG
            - 2 * detector->response[0][1] * sinGcosG
                + detector->response[1][1] * sinGsinG;
      detT->d22 = detector->response[0][0] * sinGsinG
            + 2 * detector->response[0][1] * sinGcosG
                + detector->response[1][1] * cosGcosG;
      detT->d12 = (detector->response[0][0] - detector->response[1][1])
                                           * sinGcosG
                + detector->response[0][1] * (cosGcosG - sinGsinG);
      detT->d13 = detector->response[0][2] * cosG
                - detector->response[1][2] * sinG;
      detT->d23 = detector->response[0][2] * sinG
                + detector->response[1][2] * cosG;
      detT->d33 = detector->response[2][2];

      /*
      printf("d = (%f %f %f\n",detT->d11,detT->d12,detT->d13);
      printf("     %f %f %f\n",detT->d12,detT->d22,detT->d23);
      printf("     %f %f %f)\n",detT->d13,detT->d23,detT->d33);

      printf("d*= (%f %f %f\n",detector->response[0][0],
             detector->response[0][1],detector->response[0][2]);
      printf("     %f %f %f\n",detector->response[1][0],
             detector->response[1][1],detector->response[1][2]);
      printf("     %f %f %f)\n",detector->response[2][0],
             detector->response[2][1],detector->response[2][2]);
      */

    } /* if Earth-based */

  return 0;

} /* XLALFillDetectorTensorShort() */

/** (test function) get detector states for tShort */
DetectorStateSeries *
XLALGetDetectorStatesShort ( const LIGOTimeGPSVector *timestamps, /**< array of GPS timestamps t_i */
                        const LALDetector *detector,              /**< detector info */
                        const EphemerisData *edat,                /**< ephemeris file data */
                        REAL8 tOffset,                            /**< compute detector states at timestamps SHIFTED by tOffset */
                        REAL8 tShort,                             /**< [in] new time baseline, Tshort */         
                        UINT4 numShortPerDet                      /**< [in] number of Tshort per detector */
                        )
{
  /* check input consistency */
  if ( !timestamps || !detector || !edat ) {
    XLALPrintError ("%s: invalid NULL input, timestamps=%p, detector=%p, edat=%p\n", __func__, timestamps, detector, edat );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  /* prepare return vector */
  UINT4 numSteps = numShortPerDet;
  printf("numSteps in GetDetectorStatesShort: %u\n", numSteps);
  DetectorStateSeries *ret = NULL;
  if ( ( ret = XLALCreateDetectorStateSeries ( numSteps )) == NULL ) {
    XLALPrintError ("%s: XLALCreateDetectorStateSeries(%d) failed.\n", __func__, numSteps );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }

  /* enter detector-info into the head of the state-vector */
  ret->detector = (*detector);

  /* set 'time-span' associated with each timestamp */
  //ret->deltaT = timestamps->deltaT;
  ret->deltaT = tShort;

  /* set SSB coordinate system used: EQUATORIAL for Earth-based, ECLIPTIC for LISA */
  if ( detector->frDetector.prefix[0] == 'Z' )  /* LISA */
    ret->system = COORDINATESYSTEM_ECLIPTIC;
  else  /* Earth-based */
    ret->system = COORDINATESYSTEM_EQUATORIAL;

  /* now fill all the vector-entries corresponding to different timestamps */
  UINT4 i;
  for ( i=0; i < numSteps; i++ )
    {
      BarycenterInput baryinput;
      EmissionTime emit;
      DetectorState *state = &(ret->data[i]);
      EarthState *earth = &(state->earthState);
      LIGOTimeGPS tgps;

      /* shift timestamp by tOffset */
      tgps = timestamps->data[i];
      XLALGPSAdd(&tgps, tOffset);
      /*----- first get earth-state */
      if ( XLALBarycenterEarth ( earth, &tgps, edat ) != XLAL_SUCCESS ) {
        XLALDestroyDetectorStateSeries ( ret );
        XLALPrintError("%s: XLALBarycenterEarth() failed with xlalErrno=%d\n", __func__, xlalErrno );
        XLAL_ERROR_NULL ( XLAL_EFAILED );
      }

      /*----- then get detector-specific info */
      baryinput.tgps = tgps;
      baryinput.site = (*detector);
      baryinput.site.location[0] /= LAL_C_SI;
      baryinput.site.location[1] /= LAL_C_SI;
      baryinput.site.location[2] /= LAL_C_SI;
      baryinput.alpha = baryinput.delta = 0;    /* irrelevant */
      baryinput.dInv = 0;

      if ( XLALBarycenter ( &emit, &baryinput, earth) != XLAL_SUCCESS ) {
        XLALDestroyDetectorStateSeries( ret );
        XLALPrintError("%s: XLALBarycenterEarth() failed with xlalErrno=%d\n", __func__, xlalErrno );
        XLAL_ERROR_NULL ( XLAL_EFAILED );
      }

      /*----- extract the output-data from this */
      UINT4 j;
      for (j=0; j < 3; j++)     /* copy detector's position and velocity */
        {
          state->rDetector[j] = emit.rDetector[j];
          state->vDetector[j] = emit.vDetector[j];
        } /* for j < 3 */

      /* local mean sidereal time = GMST + longitude */
      state->LMST = earth->gmstRad + detector->frDetector.vertexLongitudeRadians;
      state->LMST = fmod (state->LMST, LAL_TWOPI );     /* normalize */

      /* insert timestamp */
      state->tGPS = tgps;

      /* compute the detector-tensor at this time-stamp in SSB-fixed Cartesian coordinates
       * [EQUATORIAL for Earth-based, ECLIPTIC for LISA]
       */
      /* XLALFillDetectorTensorShort should be identical to non-short
         version, simply a copy because non-short was inaccessible
         without editing header files (i.e., it is internal to other file) */
      if ( XLALFillDetectorTensorShort ( state, detector ) != 0 ) {
        XLALDestroyDetectorStateSeries(ret);
        XLALPrintError ( "%s: XLALFillDetectorTensorShort() failed ... errno = %d\n\n", __func__, xlalErrno );
        XLAL_ERROR_NULL ( XLAL_EFUNC );
      }

    } /* for i < numSteps */

  /* return result */
  return ret;

} /* XLALGetDetectorStatesShort() */

/** (test function) get multi detector states for tShort */
MultiDetectorStateSeries *
XLALGetMultiDetectorStatesShort( const MultiLIGOTimeGPSVector *multiTS, /**< [in] multi-IFO timestamps */
                            const MultiLALDetector *multiIFO,      /**< [in] multi-IFO array holding detector info */
                            const EphemerisData *edat,             /**< [in] ephemeris data */
                            REAL8 tOffset,                         /**< [in] shift all timestamps by this amount */
                            REAL8 tShort,                          /**< [in] new time baseline, Tshort */
                            UINT4 numShortPerDet                   /**< [in] number of Tshort segments per detector */
                            )
{
  /* check input consistency */
  if ( !multiIFO || !multiTS || !edat ) {
    XLALPrintError ("%s: invalid NULL input (multiIFO=%p, multiTS=%p or edat=%p)\n", __func__, multiIFO, multiTS, edat );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  UINT4 numDetectors;
  numDetectors = multiIFO->length;
  if ( numDetectors != multiTS->length ) {
    XLALPrintError ("%s: inconsistent number of IFOs in 'multiIFO' (%d) and 'multiTS' (%d)\n", __func__, multiIFO->length, multiTS->length );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  /* prepare return-structure */
  MultiDetectorStateSeries *ret = NULL;
  if ( ( ret = LALCalloc ( 1, sizeof( *ret ) )) == NULL ) {
    XLALPrintError ("%s: LALCalloc ( 1, %zu ) failed\n", __func__, sizeof(*ret) );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }
  if ( ( ret->data = LALCalloc ( numDetectors, sizeof( *(ret->data) ) )) == NULL ) {
    XLALFree ( ret );
    XLALPrintError ("%s: LALCalloc ( %d, %zu ) failed\n", __func__, numDetectors, sizeof(*(ret->data)) );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }
  ret->length = numDetectors;

  REAL8 t0=LAL_REAL4_MAX;
  REAL8 t1=0;
  //REAL8 deltaT = tShort;
  //LIGOTimeGPS startTime = {0, 0};
  /* loop over detectors */
  UINT4 X;
  for ( X=0; X < numDetectors; X ++ )
    {
      LIGOTimeGPSVector *tsX = multiTS->data[X];
      const LALDetector *detX = &(multiIFO->sites[X]);

      if ( !tsX || !detX ) {
        XLALPrintError ("%s: invalid NULL data-vector tsX[%d] = %p, detX[%d] = %p\n", __func__, X, tsX, X, detX );
        XLAL_ERROR_NULL ( XLAL_EINVAL );
      }

      /* fill in the detector-state series for this detector */
      if ( ( ret->data[X] = XLALGetDetectorStatesShort ( tsX, detX, edat, tOffset, tShort, numShortPerDet )) == NULL ) {
        XLALPrintError ("%s: XLALGetDetectorStates() failed.\n", __func__ );
        XLAL_ERROR_NULL ( XLAL_EFUNC );
      }

      /* keep track of earliest/latest timestamp in order to determine total Tspan */
      /* This remains very close to true even with Tshort */
      UINT4 numTS = tsX->length;
      REAL8 t0_X = XLALGPSGetREAL8( &tsX->data[0] );
      REAL8 t1_X = XLALGPSGetREAL8( &tsX->data[numTS-1] );

      if ( t0_X < t0 ) {
        t0 = t0_X;
        //startTime = tsX->data[0];
      }
      if ( t1_X > t1 ) t1 = t1_X;

    } /* for X < numDetectors */

  //ret->Tspan = t1 - t0 + deltaT;        /* total time spanned by all SFTs */
  //ret->startTime = startTime;           /* earliest start-time of observation */

  return ret;

} /* XLALGetMultiDetectorStatesShort() */


/** Modify amplitude weight coefficients for tShort */
int
XLALModifyAMCoeffsWeights (
                            REAL8Vector                   **        resampMultiWeightsX,     /**< [out] new weights */
                            const MultiNoiseWeights       *restrict multiWeights,            /**< [in] old weights */
                            const REAL8                             tShort,                  /**< [in] new time baseline, Tshort */
                            const REAL8                             tSFTOld,                 /**< [in] old time baseline, tSFTOld */
                            const UINT4                             numShortPerDet,          /**< [in] number of tShort segments per detector */
                            const MultiLIGOTimeGPSVector  *restrict multiTimes,              /**< [in] multi-times vector to tell us when the SFTs were */
                            const UINT4                             maxNumStepsOldIfGapless, /**< [in] how many SFTs there would be without gaps */
                            const UINT4                             X                        /**< [in] detector number */
)
{
      REAL8Vector * ret = NULL;
      UINT4 numStepsXNew = numShortPerDet; 
      COMPLEX8Vector *newWeightsX = NULL;
      REAL8Vector *newWeightStamps = NULL;
      REAL8Vector *newWeightsXReal = NULL;

      XLAL_CHECK ( (newWeightsXReal = XLALCreateREAL8Vector ( numStepsXNew )) != NULL, XLAL_EFUNC );
      XLAL_CHECK ( (newWeightsX = XLALCreateCOMPLEX8Vector ( numStepsXNew )) != NULL, XLAL_EFUNC );
      XLAL_CHECK ( (newWeightStamps = XLALCreateREAL8Vector ( numStepsXNew )) != NULL, XLAL_EFUNC );
      COMPLEX8TimeSeries *oldTSWeightsX = NULL;

      for ( UINT4 alphaFuture = 0; alphaFuture < numStepsXNew; alphaFuture++ ){
          newWeightStamps->data[alphaFuture] = alphaFuture*tShort;
      }
      REAL8Vector *weightsX = multiWeights->data[X];
      LIGOTimeGPS weightEpoch = {0, 0};
      /* The idea: we are going to sinc interpolate to find what the
       * weights should be for the coefficients. However, we need to
       * first make the old weights, in weightsX, into an evenly-sampled
       * series, with zeroes filling any gaps. That is oldTSWeightX. We
       * then can interpolate oldTSWeightX into newWeightsX, with
       * spacing according to tShort instead of tSFTOld. */
      oldTSWeightsX = XLALCreateCOMPLEX8TimeSeries("old TS weights", &weightEpoch, 0, tSFTOld, &lalDimensionlessUnit, maxNumStepsOldIfGapless);
      UINT4 jj = 0;
      REAL8 tAlpha = 0;
      REAL8 tAlphaOldWithGaps = 0;
      /* The for-loop is over alphaLoad, an index at a cadence of the
       * old SFTs, as if they did not have gaps */
      for ( UINT4 alphaLoad = 0; alphaLoad < maxNumStepsOldIfGapless; alphaLoad++){
          /* tAlpha is the time in this gapless sequence*/
          tAlpha = alphaLoad * tSFTOld;
          /* tAlphaOldWithGaps is the actual SFT time from timestamps */
          tAlphaOldWithGaps = XLALGPSDiff( &(multiTimes->data[X]->data[jj]), &(multiTimes->data[X]->data[0]));
          for (UINT4 kk = jj; kk < multiTimes->data[X]->length ; kk++){
              /* Will overshoot gap to next SFT when reached */
              tAlphaOldWithGaps = XLALGPSDiff( &(multiTimes->data[X]->data[kk]), &(multiTimes->data[X]->data[0]));            
              if ( tAlphaOldWithGaps <= tAlpha){
                  /* In this branch if we are at least the time of first old SFT */
                  jj = kk;
              }else{
                  /* in this branch, we do not match any old sfts */
                  break;
              }
              if ( (tAlpha - tAlphaOldWithGaps) < tSFTOld){
                  /* Standard progression without gaps, assuming constant cadence */
                  oldTSWeightsX->data->data[alphaLoad] = weightsX->data[jj] + 0.0*I;
              }else{
                  /* in this branch, we do not match match the next sft:
                   * probably there was a gap */
                  oldTSWeightsX->data->data[alphaLoad] = 0.0 + 0.0*I;
              }
          } // for kk
    } // for alphaLoad
    XLAL_CHECK ( XLALSincInterpolateCOMPLEX8TimeSeries( newWeightsX, newWeightStamps, oldTSWeightsX, 8) == XLAL_SUCCESS, XLAL_EFUNC ); 
    /* Note that this above sequence differs from that in the 
     * timestamps, because here we simply transcribed the old weights
     * into a zero-padded sequence at the same cadence and used the 
     * sinc interpolation to do the heavy-lifting into the new
     * timeseries at new cadence. In the timestamps, we were creating, 
     * in effect, a Boolean sequence in the next time series from 
     * branch checks alone, with no sinc interpolation. */
    for ( UINT4 alphaReal = 0; alphaReal < numStepsXNew; alphaReal++){
          /* Weights must never be negative or else next loop yields nans*/
          newWeightsXReal->data[alphaReal] = MYMAX(0, creal(newWeightsX->data[alphaReal]));
    }
    ret = newWeightsXReal;
    XLALDestroyCOMPLEX8TimeSeries( oldTSWeightsX );
    XLALDestroyCOMPLEX8Vector( newWeightsX );
    XLALDestroyREAL8Vector( newWeightStamps );
    (*resampMultiWeightsX) = ret;
    return XLAL_SUCCESS;
} // XLALModifyAMCoeffsWeights


/** Modify multiple detectors' amplitude weight coefficients for tShort */
int
XLALModifyMultiAMCoeffsWeights (  
                            MultiNoiseWeights            **         multiWeights,   /**< [in/out] old and new weights */
                            const REAL8                             tShort,         /**< [in] new time baseline, Tshort */
                            const REAL8                             tSFTOld,        /**< [in] old time baseline, tSFTOld */
                            const UINT4                             numShortPerDet, /**< [in] number of tShort segments per detector */
                            const MultiLIGOTimeGPSVector  *restrict multiTimes      /**< [in] multi-times vector to tell us when the SFTs were */
    )
{

  /* ----- input sanity checks ----- */
  MultiNoiseWeights *multiWeightsLocal = (*multiWeights);
  UINT4 numDetectors = multiWeightsLocal->length;
  UINT4 numIFOs = numDetectors;

  REAL8 ratioSFTs = tShort / tSFTOld;

  /* Prepare output values */
  /* create multi noise weights for output */

  MultiNoiseWeights *ret = NULL;
  XLAL_CHECK ( (ret = XLALCalloc(1, sizeof(*ret))) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( (ret->data = XLALCalloc ( numIFOs, sizeof(*ret->data))) != NULL, XLAL_EFUNC );
  ret->length = numIFOs;
  /* As documented for MultiNoiseWeights, the S inv Tsft field is a
   * normalization factor equal to S^(-1) Tsft. Since we changed the
   * effective Tsft from Tsft to Tshort = ratioSFTs*Tsft, we also
   * need to multiply this normalization factor by ratioSFTs. */
  ret->Sinv_Tsft = multiWeightsLocal->Sinv_Tsft * ratioSFTs;

  REAL8 Tobs = numShortPerDet * tShort;
  UINT4 maxNumStepsOldIfGapless = lround(Tobs/tSFTOld);

  /* ---------- main loop over detectors X ---------- */
  UINT4 X;

  for ( X=0; X < numDetectors; X ++)
    {
      XLALModifyAMCoeffsWeights(&ret->data[X], multiWeightsLocal, tShort, tSFTOld, numShortPerDet, multiTimes, maxNumStepsOldIfGapless, X);
    }
   
  XLALDestroyMultiNoiseWeights ( (*multiWeights) );
  (*multiWeights) = ret;
  return XLAL_SUCCESS;

} /* XLALModifyMultiAMCoeffsWeights() */


/** (test function) used for weighting multi amplitude modulation coefficients */
int
XLALWeightMultiAMCoeffsShort (  
                            MultiAMCoeffs *multiAMcoef,            /**< [in/out] amplitude coefficients */
                            const MultiNoiseWeights *multiWeights, /**< [in/out] weights */
                            REAL8 tShort,                          /**< [in] new time baseline, Tshort */
                            REAL8 tSFTOld,                         /**< [in] old time baseline, tSFTOld */
                            UINT4 numShortPerDet,                  /**< [in] number of tShort segments per detector */
                            MultiLIGOTimeGPSVector *multiTimes /**< [in] multi-times vector to tell us when the SFTs were */
    )
{

  /* ----- input sanity checks ----- */
  if ( !multiAMcoef ) {
    XLALPrintError ("%s: illegal NULL input received in 'multiAMcoefs'.\n", __func__ );
    XLAL_ERROR( XLAL_EINVAL );
  }
  UINT4 numDetectors = multiAMcoef->length;
  /* make sure identical number of detectors in amCoefs and weights */
  if ( multiWeights && (multiWeights->length != numDetectors) ) {
    XLALPrintError("%s: multiWeights must be NULL or have the same number of detectors (numDet=%d) as mulitAMcoef (numDet=%d)!\n", __func__, multiWeights->length, numDetectors );
    XLAL_ERROR( XLAL_EINVAL );
  }

  REAL8 ratioSFTs = tShort / tSFTOld;

  REAL4 Ad = 0, Bd = 0, Cd = 0; // multi-IFO values

  REAL8 Tobs = numShortPerDet * tShort;
  UINT4 maxNumStepsOldIfGapless = lround(Tobs/tSFTOld);

  /* ---------- main loop over detectors X ---------- */
  UINT4 X;
  for ( X=0; X < numDetectors; X ++)
    {
      AMCoeffs *amcoeX = multiAMcoef->data[X];
      UINT4 numStepsXNew = numShortPerDet; 
      if ( multiWeights )
        {
          COMPLEX8Vector *newWeightsX = NULL;
          REAL8Vector *newWeightStamps = NULL;
          REAL8Vector *newWeightsXReal = NULL;
          XLAL_CHECK ( (newWeightsXReal = XLALCreateREAL8Vector ( numStepsXNew )) != NULL, XLAL_EFUNC );
          XLAL_CHECK ( (newWeightsX = XLALCreateCOMPLEX8Vector ( numStepsXNew )) != NULL, XLAL_EFUNC );
          XLAL_CHECK ( (newWeightStamps = XLALCreateREAL8Vector ( numStepsXNew )) != NULL, XLAL_EFUNC );
          for ( UINT4 alphaFuture = 0; alphaFuture < numStepsXNew; alphaFuture++ ){
              newWeightStamps->data[alphaFuture] = alphaFuture*tShort;
          }
          REAL8Vector *weightsX = multiWeights->data[X];
          COMPLEX8TimeSeries *oldTSWeightsX = NULL;
          LIGOTimeGPS weightEpoch = {0, 0};
          /* The idea: we are going to sinc interpolate to find what the
           * weights should be for the coefficients. However, we need to
           * first make the old weights, in weightsX, into an evenly-sampled
           * series, with zeroes filling any gaps. That is oldTSWeightX. We
           * then can interpolate oldTSWeightX into newWeightsX, with
           * spacing according to tShort instead of tSFTOld. */
          oldTSWeightsX = XLALCreateCOMPLEX8TimeSeries("old TS weights", &weightEpoch, 0, tSFTOld, &lalDimensionlessUnit, maxNumStepsOldIfGapless);
          UINT4 jj = 0;
          REAL8 tAlpha = 0;
          REAL8 tAlphaOldWithGaps = 0;
          /* The for-loop is over alphaLoad, an index at a cadence of the
           * old SFTs, as if they did not have gaps */
          for ( UINT4 alphaLoad = 0; alphaLoad < maxNumStepsOldIfGapless; alphaLoad++){
              /* tAlpha is the time in this gapless sequence*/
              tAlpha = alphaLoad * tSFTOld;
              /* tAlphaOldWithGaps is the actual SFT time from timestamps */
              tAlphaOldWithGaps = XLALGPSDiff( &(multiTimes->data[X]->data[jj]), &(multiTimes->data[X]->data[0]));
              for (UINT4 kk = jj; kk < multiTimes->data[X]->length ; kk++){
                  /* Will overshoot gap to next SFT when reached */
                  tAlphaOldWithGaps = XLALGPSDiff( &(multiTimes->data[X]->data[kk]), &(multiTimes->data[X]->data[0])); 
                  
                  if ( tAlphaOldWithGaps <= tAlpha){
                      /* In this branch if we are at least the time of first old SFT */
                      jj = kk;
                  }else{
                      /* in this branch, we do not match any old sfts */
                      //oldTSWeightsX->data->data[alphaLoad] = 0.0 + 0.0*I;
                      break;
                  }
                  if ( (tAlpha - tAlphaOldWithGaps) < tSFTOld){
                      /* Standard progression without gaps, assuming constant cadence */
                      oldTSWeightsX->data->data[alphaLoad] = weightsX->data[jj] + 0.0*I;
                  }else{
                      /* in this branch, we do not match match the next sft:
                       * probably there was a gap */
                      oldTSWeightsX->data->data[alphaLoad] = 0.0 + 0.0*I;
                  }
              } // for kk

          } // for alphaLoad
          XLAL_CHECK ( XLALSincInterpolateCOMPLEX8TimeSeries( newWeightsX, newWeightStamps, oldTSWeightsX, 8) == XLAL_SUCCESS, XLAL_EFUNC ); 
          /* Note that this above sequence differs from that in the 
           * timestamps, because here we simply transcribed the old weights
           * into a zero-padded sequence at the same cadence and used the 
           * sinc interpolation to do the heavy-lifting into the new
           * timeseries at new cadence. In the timestamps, we were creating, 
           * in effect, a Boolean sequence in the next time series from 
           * branch checks alone, with no sinc interpolation. */

          for ( UINT4 alphaReal = 0; alphaReal < numStepsXNew; alphaReal++){
              /* Weights must never be negative or else next loop yields nans*/
              newWeightsXReal->data[alphaReal] = MYMAX(0, creal(newWeightsX->data[alphaReal]));
          }
          UINT4 alpha;
          for(alpha = 0; alpha < numStepsXNew; alpha++)
            {
              REAL8 Sqwi = sqrt ( newWeightsXReal->data[alpha] );
              amcoeX->a->data[alpha] *= Sqwi;
              amcoeX->b->data[alpha] *= Sqwi;
            } // for alphaOld < numSteps
          XLALDestroyCOMPLEX8TimeSeries( oldTSWeightsX );
          XLALDestroyCOMPLEX8Vector( newWeightsX );
          XLALDestroyREAL8Vector( newWeightsXReal );
          XLALDestroyREAL8Vector( newWeightStamps );
        } // if weights


      UINT4 alpha;      // SFT-index
      REAL4 AdX = 0, BdX = 0, CdX = 0;  // single-IFO values
      /* compute single-IFO antenna-pattern coefficients AX,BX,CX, by summing over time-steps 'alpha' */
      for(alpha = 0; alpha < numStepsXNew; alpha++)
        {
          REAL4 ahat = amcoeX->a->data[alpha];
          REAL4 bhat = amcoeX->b->data[alpha];

          AdX += ahat * ahat;
          BdX += bhat * bhat;
          CdX += ahat * bhat;
          //printf("Final alpha! %u\n", alpha);
        } /* for alpha < numStepsXNew */

      /* store those */
      amcoeX->A = AdX;
      amcoeX->B = BdX;
      amcoeX->C = CdX;
      amcoeX->D = AdX * BdX - CdX * CdX;

      // in the unlikely event of a degenerate M-matrix with D = det[A, C; C, B] <= 0,
      // we set D->inf, in order to set the corresponding F-value to zero rather than >>1
      // By setting 'D=inf', we also allow upstream catching/filtering on such singular cases
      if ( amcoeX->D <= 0 ) {
        amcoeX->D = INFINITY;
      }
      /* compute multi-IFO antenna-pattern coefficients A,B,C by summing over IFOs X */
      Ad += AdX;
      Bd += BdX;
      Cd += CdX;

    } /* for X < numDetectors */
    //printf("Final A, B, C, D: %f %f %f %f\n", Ad, Bd, Cd, Ad * Bd - Cd * Cd);

  multiAMcoef->Mmunu.Ad = Ad;
  multiAMcoef->Mmunu.Bd = Bd;
  multiAMcoef->Mmunu.Cd = Cd;
  multiAMcoef->Mmunu.Dd = Ad * Bd - Cd * Cd;

  // in the unlikely event of a degenerate M-matrix with D = det[A, C; C, B] <= 0,
  // we set D->inf, in order to set the corresponding F-value to zero rather than >>1
  // By setting 'D=inf', we also allow upstream catching/filtering on such singular cases
  if ( multiAMcoef->Mmunu.Dd <= 0 ) {
    multiAMcoef->Mmunu.Dd = INFINITY;
  }


  if ( multiWeights ) {
    //multiAMcoef->Mmunu.Sinv_Tsft = multiWeights->Sinv_Tsft;
    multiAMcoef->Mmunu.Sinv_Tsft = multiWeights->Sinv_Tsft * ratioSFTs;
  }

  return XLAL_SUCCESS;

} /* XLALWeightMultiAMCoeffsShort() */

/** (test function) used for computing amplitude modulation weights */
AMCoeffs *
XLALComputeAMCoeffsShort ( const DetectorStateSeries *DetectorStates,        /**< timeseries of detector states */
                      SkyPosition skypos                                /**< {alpha,delta} of the source */
                      )
{
  /* ---------- check input consistency ---------- */
  if ( !DetectorStates ) {
    XLALPrintError ("%s: invalid NULL input 'DetectorStates'\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  /* currently requires sky-pos to be in equatorial coordinates (FIXME) */
  if ( skypos.system != COORDINATESYSTEM_EQUATORIAL ) {
    XLALPrintError ("%s: only equatorial coordinates currently supported in 'skypos'\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  /*---------- We write components of xi and eta vectors in SSB-fixed coords */
  REAL4 alpha = skypos.longitude;
  REAL4 delta = skypos.latitude;

  REAL4 sin1delta, cos1delta;
  REAL4 sin1alpha, cos1alpha;
  XLAL_CHECK_NULL( XLALSinCosLUT (&sin1delta, &cos1delta, delta ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_NULL( XLALSinCosLUT (&sin1alpha, &cos1alpha, alpha ) == XLAL_SUCCESS, XLAL_EFUNC );

  REAL4 xi1 = - sin1alpha;
  REAL4 xi2 =  cos1alpha;
  REAL4 eta1 = sin1delta * cos1alpha;
  REAL4 eta2 = sin1delta * sin1alpha;
  REAL4 eta3 = - cos1delta;

  /* prepare output vector */
  UINT4 numSteps = DetectorStates->length;
  //printf("numSteps in ComputeAMCoeffsShort: %u\n", numSteps);
  AMCoeffs *coeffs;
  if ( ( coeffs = XLALCreateAMCoeffs ( numSteps ) ) == NULL ) {
    XLALPrintError ("%s: XLALCreateAMCoeffs(%d) failed\n", __func__, numSteps );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }

  /*---------- Compute the a(t_i) and b(t_i) ---------- */
  UINT4 i;
  for ( i=0; i < numSteps; i++ )
    {
      REAL4 ai, bi;

      SymmTensor3 *d = &(DetectorStates->data[i].detT);

      ai =    d->d11 * ( xi1 * xi1 - eta1 * eta1 )
        + 2 * d->d12 * ( xi1*xi2 - eta1*eta2 )
        - 2 * d->d13 *             eta1 * eta3
        +     d->d22 * ( xi2*xi2 - eta2*eta2 )
        - 2 * d->d23 *             eta2 * eta3
        -     d->d33 *             eta3*eta3;

      bi =    d->d11 * 2 * xi1 * eta1
        + 2 * d->d12 *   ( xi1 * eta2 + xi2 * eta1 )
        + 2 * d->d13 *     xi1 * eta3
        +     d->d22 * 2 * xi2 * eta2
        + 2 * d->d23 *     xi2 * eta3;

      coeffs->a->data[i] = ai;
      coeffs->b->data[i] = bi;
      //printf("ai, bi: %f %f\n", ai, bi);

    } /* for i < numSteps */

  /* return the result */
  return coeffs;

} /* XLALComputeAMCoeffsShort() */


/** (test function) used for computing multi amplitude modulation weights */
MultiAMCoeffs *
XLALComputeMultiAMCoeffsShort ( 
                           const MultiDetectorStateSeries *multiDetStates,      /**< [in] detector-states at timestamps t_i */
                           const MultiNoiseWeights *multiWeights,               /**< [in] noise-weights at timestamps t_i (can be NULL) */
                           SkyPosition skypos,                                   /**< source sky-position [in equatorial coords!] */
                           REAL8 tShort,                          /**< [in] new time baseline, Tshort */
                           REAL8 tSFTOld,                         /**< [in] old time baseline, tSFTOld */
                           UINT4 numShortPerDet,                  /**< [in] number of tShort segments per detector */
                           MultiLIGOTimeGPSVector *multiTimes /**< [in] multi-times vector to tell us when the SFTs were */
                           )
{
  /* check input consistency */
  if ( !multiDetStates ) {
    XLALPrintError ("%s: invalid NULL input argument 'multiDetStates'\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  UINT4 numDetectors = multiDetStates->length;

  /* prepare output vector */
  MultiAMCoeffs *ret;
  if ( ( ret = XLALCalloc( 1, sizeof( *ret ) )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCalloc( 1, %zu)\n", __func__, sizeof( *ret ) );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  ret->length = numDetectors;
  if ( ( ret->data = XLALCalloc ( numDetectors, sizeof ( *ret->data ) )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCalloc(%d, %zu)\n", __func__, numDetectors, sizeof ( *ret->data ) );
    XLALFree ( ret );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  /* loop over detectors and generate AMCoeffs for each one */
  UINT4 X;
  for ( X=0; X < numDetectors; X ++ )
    {
      if ( (ret->data[X] = XLALComputeAMCoeffsShort ( multiDetStates->data[X], skypos )) == NULL ) {
        XLALPrintError ("%s: call to XLALComputeAMCoeffs() failed with xlalErrno = %d\n", __func__, xlalErrno );
        XLALDestroyMultiAMCoeffs ( ret );
        XLAL_ERROR_NULL ( XLAL_EFUNC );
      }

    } /* for X < numDetectors */

  /* Note that the weighting is the only place where multiTimes are used,
   * and for that, the old vector is the right one, because we need to
   * know when the old SFTs were */
  /* apply noise-weights and compute antenna-pattern matrix {A,B,C} */
  if ( XLALWeightMultiAMCoeffsShort (  ret, multiWeights, tShort, tSFTOld, numShortPerDet, multiTimes ) != XLAL_SUCCESS ) {
    /* Turns out tShort and tSFTOld are, surprisingly, irrelevant (I think) */
    XLALPrintError ("%s: call to XLALWeightMultiAMCoeffs() failed with xlalErrno = %d\n", __func__, xlalErrno );
    XLALDestroyMultiAMCoeffs ( ret );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }

  /* return result */
  return ret;

} /* XLALComputeMultiAMCoeffsShort() */


/** Compute the number of tShort segments per detector */
UINT4
XLALCrossCorrNumShortPerDetector(
    const REAL8 resampTshort, /**< [in] resampling Tshort */
    const INT4  startTime,    /**< [in] start time of observation */
    const INT4  endTime       /**< [in] end time of observation */
)
{
  UINT4 numShortPerDet = 0;
  REAL8 Tobs = (REAL8)(endTime - startTime);
  numShortPerDet = lround( Tobs /resampTshort);
  return numShortPerDet;
} /* XLALCrossCorrNumShortPerDetector */

/** Find gaps in the data given the SFTs*/
REAL8TimeSeries *
XLALCrossCorrGapFinderResamp(
                        LIGOTimeGPSVector        *restrict timestamps, /**< [in] timestamps vector */
                        const LIGOTimeGPSVector  *restrict Times       /**< [in] input SFT-vector  */
)
{
  UINT4 numSFTsNew = timestamps->length;
  REAL8TimeSeries *retFlag = NULL;
  /* creation and memory allocation for the science-or-not vector */
  if ( ( retFlag = XLALCreateREAL8TimeSeries("science-or-not Boolean-as-real vector", &(timestamps->data[0]), 0, (timestamps->deltaT), &lalDimensionlessUnit, numSFTsNew )  ) == NULL ) {
    XLALPrintError ("%s: XLALCreateREAL8TimeSeries(%d) failed.\n", __func__, numSFTsNew );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }
  REAL8 Tsft = Times->deltaT;
  UINT4 jj = 0;
  REAL8 tAlpha = 0;
  REAL8 tAlphaOldWithGaps = 0;
  /* Now make the science-or-not vector */
  for ( UINT4 ii = 0; ii < numSFTsNew; ii++){
    /* tAlpha is the time in this gapless sequence,
     * equivalent to corresponding epoch holder */
    /* for tAlpha, use tShort*/
    tAlpha = XLALGPSDiff( &(timestamps->data[ii]), &(timestamps->data[0]) );
    /* tAlphaOldWithGaps is the actual SFT time from timestamps */
    tAlphaOldWithGaps = XLALGPSDiff( &(Times->data[jj]), &(Times->data[0]));
    /* A science segment begins with an old stamp, continues, and ends
     * when the next stamp is more than 1 usual SFT time away -- the
     * ending time is that one SFT time */
    /* First goal: find out what is the last tAlphaOldWithGaps <= tAlpha;
     * then can come the second goal: is tAlpha - tAlphaOldWithGaps < Tsft?
     * If yes, then is science (flagHolder = 1), else not (flagHolder = 0) */
    for (UINT4 kk = jj; kk < Times->length ; kk++){
        tAlphaOldWithGaps = XLALGPSDiff( &(Times->data[kk]), &(Times->data[0]));
        if ( tAlphaOldWithGaps <= tAlpha){
            jj = kk;
            /* At this point, tAlphaOldWithGaps is right for this kk,
             * and jj has been set */
        } else {
            /* Then it has gone too far */
            break;
        } /* Goal 1: done */
        if ( tAlpha - tAlphaOldWithGaps < Tsft ) {
            /* time is within an SFT and thus science */
            retFlag->data->data[ii] = 1.0;
        } else {
            /* time is outside an SFT and thus not science */
            retFlag->data->data[ii] = 0.0;
        } /* Goal 2: done */
    }
  }
  return retFlag;
} /* XLALCrossCorrGapFinderResamp */

/** (test function) find gaps in the data given the SFTs */
//LIGOTimeGPSVector *
REAL8TimeSeries *
XLALCrossCorrGapFinderResampAlt(
                        LIGOTimeGPSVector  *restrict timestamps, /**< [in] timestamps vector */
                        const SFTVector    *restrict sfts        /**< [in] input SFT-vector  */
)
{
  UINT4 numSFTsNew = timestamps->length;
  REAL8TimeSeries *retFlag = NULL;
  /* creation and memory allocation for the science-or-not vector */
  if ( ( retFlag = XLALCreateREAL8TimeSeries("science-or-not Boolean-as-real vector", &(timestamps->data[0]), 0, (timestamps->deltaT), &lalDimensionlessUnit, numSFTsNew )  ) == NULL ) {
    XLALPrintError ("%s: XLALCreateREAL8TimeSeries(%d) failed.\n", __func__, numSFTsNew );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }
  REAL8 Tsft = 1.0 / sfts->data[0].deltaF;
  UINT4 jj = 0;
  REAL8 tAlpha = 0;
  REAL8 tAlphaOldWithGaps = 0;
  /* Now make the science-or-not vector */
  for ( UINT4 ii = 0; ii < numSFTsNew; ii++){
    /* tAlpha is the time in this gapless sequence,
     * equivalent to corresponding epoch holder */
    /* for tAlpha, use tShort*/
    tAlpha = XLALGPSDiff( &(timestamps->data[ii]), &(timestamps->data[0]) );
    /* tAlphaOldWithGaps is the actual SFT time from timestamps */
    tAlphaOldWithGaps = XLALGPSDiff( &(sfts->data[jj].epoch), &(sfts->data[0]).epoch);
    /* A science segment begins with an old stamp, continues, and ends
     * when the next stamp is more than 1 usual SFT time away -- the
     * ending time is that one SFT time */
    /* First goal: find out what is the last tAlphaOldWithGaps <= tAlpha;
     * then can come the second goal: is tAlpha - tAlphaOldWithGaps < Tsft?
     * If yes, then is science (flagHolder = 1), else not (flagHolder = 0) */
    for (UINT4 kk = jj; kk < sfts->length ; kk++){
        tAlphaOldWithGaps = XLALGPSDiff( &(sfts->data[kk].epoch), &(sfts->data[0]).epoch);
        if ( tAlphaOldWithGaps <= tAlpha){
            jj = kk;
            /* At this point, tAlphaOldWithGaps is right for this kk,
             * and jj has been set */
        } else {
            /* Then it has gone too far */
            break;
        } /* Goal 1: done */
        if ( tAlpha - tAlphaOldWithGaps < Tsft ) {
            /* time is within an SFT and thus science */
            retFlag->data->data[ii] = 1.0; 
        } else {
            /* time is outside an SFT and thus not science */
            retFlag->data->data[ii] = 0.0; 
        } /* Goal 2: done */
    }    
  }
  return retFlag;
} /* XLALCrossCorrGapFinderResampAlt */


/** Demarcate pairs with flags about whether data exists in zero-padded timeseries */
int
XLALEquipCrossCorrPairsWithScienceFlags( 
            MultiResampSFTPairMultiIndexList  * resampMultiPairs,  /**< [in] resampling pairs */ 
            MultiREAL8TimeSeries              * scienceFlagVect    /**< [in] science flags */
)
{
  /* Note that some functions see the science segment gaps implicitly, 
   * as with the gammas function being interpolated with zeroes at gaps,
   * and the gamma function thus affecting, via the coefficient of pairs
   * of phase derivatives, the metric. For these purposes, the science
   * flags serve simply as a sanity check. */
  REAL4 castSciFlag1 = 0;
  REAL4 castSciFlag2 = 0;
  UINT4 detInd1;
  UINT4 detInd2;
  UINT4 sftInd1;
  UINT4 sftInd2;
  for (UINT4 detX=0; detX < resampMultiPairs->length; detX++){
      detInd1 = resampMultiPairs->data[detX].detInd;
    for (UINT4 k=0; k <  resampMultiPairs->data[detX].length; k++){
      sftInd1 = resampMultiPairs->data[detX].data[k].sftInd;
      castSciFlag1 = (REAL4)scienceFlagVect->data[detInd1]->data->data[sftInd1];
      resampMultiPairs->data[detX].data[k].sciFlag = castSciFlag1;
      for (UINT4 detY=0; detY < resampMultiPairs->data[detX].data[k].length; detY++){
        detInd2 = resampMultiPairs->data[detX].data[k].data[detY].detInd;
        for (UINT4 l=0; l < resampMultiPairs->data[detX].data[k].data[detY].length; l++){
          sftInd2 = resampMultiPairs->data[detX].data[k].data[detY].data[l].sftInd;
          castSciFlag2 = (REAL4)scienceFlagVect->data[detInd2]->data->data[sftInd2];
          resampMultiPairs->data[detX].data[k].data[detY].data[l].sciFlag = castSciFlag2;
        }
      }
    }
  }
  return XLAL_SUCCESS;
} /* XLALEquipCrossCorrPairsWithScienceFlags */

///** (possible future function) does not work -- would adjust timestamps of an SFT vector */
//MultiSFTVector
//*XLALModifyCrossCorrTimestampsIntoSFTVector(
//    const MultiLIGOTimeGPSVector *multiTimes /**< [in] timestamps */
//)
//{
//   //XLAL_CHECK ( multiTimes != NULL, XLAL_EINVAL );
//   //XLAL_CHECK ( multiTimes->length != 0, XLAL_EINVAL );
//  
//    UINT4 numIFOs = multiTimes->length;
//    UINT4 nSFTs = 0;
//  
//    /* create multi sft vector */
//    MultiSFTVector *multiSFTs = NULL;
//    //XLAL_CHECK_NULL ( (multiSFTs = XLALCalloc(1, sizeof(*multiSFTs))) != NULL, XLAL_ENOMEM );
//    //XLAL_CHECK_NULL ( (multiSFTs->data = XLALCalloc ( numIFOs, sizeof(*multiSFTs->data))) != NULL, XLAL_ENOMEM );
//    multiSFTs = XLALCalloc(1, sizeof(*multiSFTs));
//    multiSFTs->data = XLALCalloc ( numIFOs, sizeof(*multiSFTs->data));
//    multiSFTs->length = numIFOs;
//  
//    for ( UINT4 X = 0; X < numIFOs; X++ )
//    {
//      nSFTs = multiTimes->data[X]->length;
//      multiSFTs->data[X]->length = nSFTs;
//      //SFTVector* sftVector = NULL;
//      //sftVector = XLALCreateSFTVector (nSFTs, 1);
//      //if (!(sftVector = XLALCreateSFTVector (nSFTs, 1))) {
//      //  XLALPrintError("ERROR: Couldn't create short sftVector\n");
//      //  XLALLOADSFTSERROR(XLAL_EINVAL);
//      //}
//      for ( UINT4 jj = 0; jj < nSFTs; jj++)
//      {
//        multiSFTs->data[X]->data[jj].epoch = multiTimes->data[X]->data[jj];
//      }
//    } // for X < numIFOs
//
//  return multiSFTs;
//} /* XLALModifyCrossCorrTimestampsIntoSFTVector */


/** Generates a resampling workspace for CrossCorr */
int
XLALCreateCrossCorrWorkspace( 
    ResampCrossCorrWorkspace  **        wsOut,                    /**< [out] workspace for one cross-correlation */
    COMPLEX8                  **        ws1KFaX_kOut,             /**< [out] holder for detector 1 Fa */
    COMPLEX8                  **        ws1KFbX_kOut,             /**< [out] holder for detector 1 Fb */
    COMPLEX8                  **        ws2LFaX_kOut,             /**< [out] holder for detector 2 Fa */
    COMPLEX8                  **        ws2LFbX_kOut,             /**< [out] holder for detector 2 Fb */
    MultiCOMPLEX8TimeSeries   **        multiTimeSeries_SRC_aOut, /**< [out] resampling A time series */
    MultiCOMPLEX8TimeSeries   **        multiTimeSeries_SRC_bOut, /**< [out] resampling B time series */
    const PulsarDopplerParams           binaryTemplateSpacings,   /**< [in ]binary template spacings */
    const FstatInput          *restrict resampFstatInput,         /**< [in] resampling f-statistic input */
    const UINT4                         numFreqBins,              /**< [in] number of frequency bins */
    const REAL8                         tCoh,                     /**< [in] Tcoh = 2 * max lag + resampling tShort */
    const BOOLEAN                       treatWarningsAsErrors     /**< [in] abort program if any warnings encountered */
)
{

    MultiCOMPLEX8TimeSeries  * multiTimeSeries_SRC_a = (*multiTimeSeries_SRC_aOut); 
    MultiCOMPLEX8TimeSeries  * multiTimeSeries_SRC_b = (*multiTimeSeries_SRC_bOut); 

    COMPLEX8 *restrict ws1KFaX_k = (*ws1KFaX_kOut );
    COMPLEX8 *restrict ws1KFbX_k = (*ws1KFbX_kOut );
    COMPLEX8 *restrict ws2LFaX_k = (*ws2LFaX_kOut );
    COMPLEX8 *restrict ws2LFbX_k = (*ws2LFbX_kOut );
    /* Extract base info from the resampled time series.
     * Use both a and b time series structs to make vars used */
    XLALExtractResampledTimeseries ( &multiTimeSeries_SRC_a, &multiTimeSeries_SRC_b, resampFstatInput );
    const REAL8 dt_SRC = multiTimeSeries_SRC_b->data[0]->deltaT;
    (*multiTimeSeries_SRC_aOut) = multiTimeSeries_SRC_a; 
    (*multiTimeSeries_SRC_bOut) = multiTimeSeries_SRC_b; 
    /* extract more by using the workspace */

    /* Compare with the timing model document, T1600531, for the F-stat.
     * Tcoh here matches Tcoh there, but note that the metric is defined in
     * terms of maxLag instead. Tcoh = 2*maxLag + tShort. In this stage, we
     * no longer access dt_DET, only dt_SRC.
     *   In source frame, SRCsampPerTcoh = tCoh / dt_SRC. But we can also
     * define the metric timescale,
     *   metricSRCsampPerTcoh = 1/(dFreqMetric*dt_SRC),
     * then round each up to the nearest power of 2 to define,
     *   metricNumSamplesFFT = (UINT4)pow(2, ceil(log2(metricSRCsampPerTcoh))),
     *   shortNumSamplesFFT = (UINT4)pow(2, ceil(log2(metricSRCsampPeTcoh))),
     * and choose the maximum between the two to be numSampledFFT. This is 
     * equivalent to how Fstat calculates it.
     *   To see why, note: our SRCsampPerTcoh is Fstat's N^SRC_samp. We can
     * access it immediately because the resampling is already done for us.
     * If decimateFFT = (UINT4)ceil( tCoh / TspanFFT), 
     *                = (UINT4)ceil( tCoh * dFreqMetric ),
     * then numSamplesFFT0 = (UINT4)ceil(decimateFFT * TspanFFT / dt_SRC)
     * for us -- note, dt_SRC, because we have no access to the detector.
     * and, numSamplesFFT0 = (UINT4)ceil(decimateFFT / (dFreqMetric*dt_SRC)),
     * then rounded up to the nearest power of 2, numSamplesFFT0.
     * Ergo, comingling our two notations,
     *   numSamplesFFT0 = (UINT4)ceil(decimateFFT * metricSRCsampPerTcoh),
     * and the exact same error should be thrown if decimateFFT > 1 as if
     * shortNumSamplesFFT <= metricNumSamplesFFT.
     *   Beware: we divide the resampled time-series up into sections, and
     * in the statistic this plays a crucial role.
     *   The equivalent of decimate is the quantity,
     * RedecimateFFT
     *  = binaryTemplateSpacings->fkdot[0] * (dt_SRC) * (ws->numSamplesFFT)
     *  = dFreqMetric * dt_SRC * (UINT4)pow(2, ceil(log2(...
     *      (UINT4)ceil(decimateFFT / (dFreqMetric * dt_SRC))    )))
     * which although proportional to decimateFFT is not equal to it, and
     * is a real quantity rather than integer.
     * (GDM)
     */
    /* Metric spacing: */
    const REAL8 dFreqMetric = binaryTemplateSpacings.fkdot[0];
    // determine resampled timeseries parameters
    const REAL8 TspanFFT = 1.0 / dFreqMetric;
    const UINT4 decimateFFT = (UINT4)ceil ( tCoh / TspanFFT );     // larger than 1 means we need to artificially increase dFreqFFT by 'decimateFFT'
    if ( 1 <= (dFreqMetric*dt_SRC) ) {
       printf("Warning! Metric spacing less than one FFT sample,...\n ...is maxLag < tShort? \n ...proceeding in radiometer mode...\n ...intended for non-production comparison tests only!\n");
      if ( treatWarningsAsErrors ) {
        XLALPrintError("Error! (treating warnings as errors) metric spacing less than one FFT sample.\n");
        XLAL_CHECK( 1 >= (dFreqMetric*dt_SRC), XLAL_EFUNC );
      }
    }
    if ( decimateFFT > 1 ) {
      printf("Warning! Frequency spacing larger than 1/Tcoh, decimation being applied of %" LAL_UINT4_FORMAT "\n", decimateFFT );
      if ( treatWarningsAsErrors == TRUE ){
        XLALPrintError("Error! (treating warnings as errors) FFT samples limited by tCoh (not metric), unexpected unless high mismatchF.\n");
        XLAL_CHECK( (decimateFFT < 1 ), XLAL_EFUNC);
      }
    }
    const REAL8 TspanFFT1 = decimateFFT * TspanFFT;
    /* Note difference from Fstat in using dt_SRC rather than dt_DET here */
    const UINT4 numSamplesFFT0 = (UINT4) ceil ( TspanFFT1 / dt_SRC );      // we use ceil() so that we artificially widen the band rather than reduce it
    const UINT4 numSamplesFFT = (UINT4) pow ( 2, ceil ( log2 ( numSamplesFFT0 ) ) );  // round numSamplesFFT up to next power of 2 for most effiecient FFT
    ResampCrossCorrWorkspace *restrict ws = NULL;
    /* CP'd from ComputeFstat_resamp:538 */
    const int fft_plan_flags=FFTW_MEASURE;
    //int fft_plan_flags=FFTW_ESTIMATE;
    const double fft_plan_timeout= FFTW_NO_TIMELIMIT ;
    /* Memory pre-allocate:  */
    XLAL_CHECK ( (ws = XLALCalloc ( 1, sizeof(*ws))) != NULL, XLAL_ENOMEM );
    XLAL_CHECK ( (ws->TStmp1_SRC   = XLALCreateCOMPLEX8Vector ( numSamplesFFT )) != NULL, XLAL_EFUNC );
    XLAL_CHECK ( (ws->TStmp2_SRC   = XLALCreateCOMPLEX8Vector ( numSamplesFFT )) != NULL, XLAL_EFUNC );
    XLAL_CHECK ( (ws->SRCtimes_DET = XLALCreateREAL8Vector ( numSamplesFFT )) != NULL, XLAL_EFUNC );
    XLAL_CHECK ( (ws->FaX_k = fftw_malloc ( numFreqBins * sizeof(COMPLEX8) )) != NULL, XLAL_ENOMEM );
    XLAL_CHECK ( (ws->FbX_k = fftw_malloc ( numFreqBins * sizeof(COMPLEX8) )) != NULL, XLAL_ENOMEM );
    XLAL_CHECK ( (ws->FabX_Raw = fftw_malloc ( numSamplesFFT * sizeof(COMPLEX8) )) != NULL, XLAL_ENOMEM );
    XLAL_CHECK ( (ws->TS_FFT   = fftw_malloc ( numSamplesFFT * sizeof(COMPLEX8) )) != NULL, XLAL_ENOMEM );
    ws->decimateFFT = decimateFFT;
    ws->numSamplesFFT = numSamplesFFT;
    ws->numFreqBinsOut = numFreqBins;
    /* -- create FFT plan with FFTW */
    LAL_FFTW_WISDOM_LOCK;
    //XLALGetFFTPlanHints (& fft_plan_flags , & fft_plan_timeout );
    fftw_set_timelimit( fft_plan_timeout );
    XLAL_CHECK ( (ws->fftplan = fftwf_plan_dft_1d ( numSamplesFFT, ws->TS_FFT, ws->FabX_Raw, FFTW_FORWARD, fft_plan_flags )) != NULL, XLAL_EFAILED, "fftwf_plan_dft_1d() failed\n");
    LAL_FFTW_WISDOM_UNLOCK;
    /* -- finish creating FFT plan with FFTW */

    XLAL_CHECK ( (ws1KFaX_k = fftw_malloc ( numFreqBins * sizeof(COMPLEX8) )) != NULL, XLAL_ENOMEM );
    XLAL_CHECK ( (ws1KFbX_k = fftw_malloc ( numFreqBins * sizeof(COMPLEX8) )) != NULL, XLAL_ENOMEM );
    XLAL_CHECK ( (ws2LFaX_k = fftw_malloc ( numFreqBins * sizeof(COMPLEX8) )) != NULL, XLAL_ENOMEM );
    XLAL_CHECK ( (ws2LFbX_k = fftw_malloc ( numFreqBins * sizeof(COMPLEX8) )) != NULL, XLAL_ENOMEM );
    (*ws1KFaX_kOut ) = ws1KFaX_k;
    (*ws1KFbX_kOut ) = ws1KFbX_k;
    (*ws2LFaX_kOut ) = ws2LFaX_k;
    (*ws2LFbX_kOut ) = ws2LFbX_k;
    (*wsOut) = ws;
    return XLAL_SUCCESS;
} /* XLALCreateCrossCorrWorkspace */

/* ===== Object destruction functions ===== */

/**
 * Destroy a SFTIndexList structure.
 * Note, this is "NULL-robust" in the sense that it will not crash
 * on NULL-entries anywhere in this struct, so it can be used
 * for failure-cleanup even on incomplete structs
 */
void
XLALDestroySFTIndexList ( SFTIndexList *sftIndices )
{
  if ( ! sftIndices )
    return;

  if ( sftIndices->data )
    XLALFree(sftIndices->data);

  XLALFree ( sftIndices );

  return;

} /* XLALDestroySFTIndexList() */

/**
 * Destroy a SFTPairIndexList structure.
 * Note, this is "NULL-robust" in the sense that it will not crash
 * on NULL-entries anywhere in this struct, so it can be used
 * for failure-cleanup even on incomplete structs
 */
void
XLALDestroySFTPairIndexList ( SFTPairIndexList *sftPairs )
{
  if ( ! sftPairs )
    return;

  if ( sftPairs->data )
    XLALFree(sftPairs->data);

  XLALFree ( sftPairs );

  return;

} /* XLALDestroySFTPairIndexList() */

void XLALDestroyResampSFTIndexList( ResampSFTIndexList *sftResampList )
{
  if ( !sftResampList )
    return;
  
  if ( sftResampList->data )
    XLALFree(sftResampList->data);

  //XLALFree( sftResampList );
 
  return;
} /* XLALDestroyResampSFTIndexList */

void XLALDestroyResampSFTMultiIndexList( ResampSFTMultiIndexList *sftResampMultiList )
{
  if ( !sftResampMultiList )
    return;

      for (UINT4 detY=0; detY < sftResampMultiList->length; detY++){
          XLALDestroyResampSFTIndexList( & sftResampMultiList->data[detY] );
      }

  if ( sftResampMultiList->data )
    XLALFree(sftResampMultiList->data);


  return;
} /* XLALDestroyResampSFTMultiIndexList */

void XLALDestroyResampSFTPairMultiIndexList(  ResampSFTPairMultiIndexList *sftResampPairMultiList  )
{
  if ( !sftResampPairMultiList )
    return;

    for (UINT4 k=0; k <  sftResampPairMultiList->length; k++){
      XLALDestroyResampSFTMultiIndexList(  & sftResampPairMultiList->data[k]  );
    }

  if ( sftResampPairMultiList->data )
    XLALFree(sftResampPairMultiList->data);


  return;
} /* XLALDestroyResampSFTPairMultiIndexList */

void
XLALDestroyMultiResampSFTPairMultiIndexList ( MultiResampSFTPairMultiIndexList *sftMultiPairsResamp)
{
  if ( ! sftMultiPairsResamp )
    return;

  /* Iterate through the structure and destroy completely */
  for (UINT4 detX=0; detX < sftMultiPairsResamp->length; detX++){
    XLALDestroyResampSFTPairMultiIndexList(  & sftMultiPairsResamp->data[detX]  );
  }

  /* Rely on external functions to clear pairs and flat index
     for now, to make robust */

  if ( sftMultiPairsResamp->data ){
    XLALFree( sftMultiPairsResamp->data);
 }

  XLALFree( sftMultiPairsResamp );

  return;

} /* XLALDestroyMultiResampSFTPairMultiIndexList() */

/* Destroy the struct that matched pairs */
void
XLALDestroyMultiMatchList( MultiResampSFTMultiCountList *localMultiListOfLmatchingGivenMultiK ){

  if ( ! localMultiListOfLmatchingGivenMultiK )
    return;

  UINT4 numDets = localMultiListOfLmatchingGivenMultiK->length; /* Because it is a multi-multi-list */
  /* furthermore, scroll through each detector X to read how many SFTs K_X it has */
  for (UINT4 detX=0; detX < numDets; detX++){
    for (UINT4 k=0; k < localMultiListOfLmatchingGivenMultiK->data[detX].length; k++){
      if ( localMultiListOfLmatchingGivenMultiK->data[detX].data[k].data ) {
        XLALFree ( localMultiListOfLmatchingGivenMultiK->data[detX].data[k].data);
      }
    }
    if (  localMultiListOfLmatchingGivenMultiK->data[detX].data ) {
      XLALFree ( localMultiListOfLmatchingGivenMultiK->data[detX].data);
    }
  }

  if( localMultiListOfLmatchingGivenMultiK->data ){
    XLALFree ( localMultiListOfLmatchingGivenMultiK->data );
  }

  XLALFree ( localMultiListOfLmatchingGivenMultiK );

}

void
XLALDestroyResampCrossCorrWorkspace ( void *workspace )
{
  ResampCrossCorrWorkspace *ws = (ResampCrossCorrWorkspace*) workspace;

  XLALDestroyCOMPLEX8Vector ( ws->TStmp1_SRC );
  XLALDestroyCOMPLEX8Vector ( ws->TStmp2_SRC );
  XLALDestroyREAL8Vector ( ws->SRCtimes_DET );

  LAL_FFTW_WISDOM_LOCK;
  fftwf_destroy_plan ( ws->fftplan );
  LAL_FFTW_WISDOM_UNLOCK;

  fftw_free ( ws->FabX_Raw );
  fftw_free ( ws->TS_FFT );

  XLALFree ( ws->FaX_k );
  XLALFree ( ws->FbX_k );
  XLALFree ( ws->Fa_k );
  XLALFree ( ws->Fb_k );

  XLALFree ( ws );
  return;

} // XLALDestroyResampWorkspace()

// ---------- internal functions ----------

/** Imported and modified from ComputeFstat_Resamp.c */
static int
XLALApplyCrossCorrFreqShiftResamp
  ( 
    COMPLEX8                   *restrict xOut,             ///< [out] the spindown-corrected SRC-frame timeseries
    const COMPLEX8TimeSeries   *restrict xIn,              ///< [in] the input SRC-frame timeseries
    const PulsarDopplerParams  *restrict doppler,          ///< [in] containing spindown parameters
    const REAL8                          freqShift,        ///< [in] frequency-shift to apply, sign is "new - old"
    const UINT4                          indexStartResamp, ///< [in] index in resampling time series to start
    const UINT4                          indexEndResamp,   ///< [in] index in resampling time series to end
    const UINT4                          numSamplesIn,     ///< [in] number of samples in the output
    const UINT4                          insertPoint       ///< [in] number of sample indices offset to begin output
  )
{
  // input sanity checks
  XLAL_CHECK ( xOut != NULL, XLAL_EINVAL );
  XLAL_CHECK ( xIn != NULL, XLAL_EINVAL );
  XLAL_CHECK ( doppler != NULL, XLAL_EINVAL );

  //// determine number of spin downs to include
  //UINT4 s_max = PULSAR_MAX_SPINS - 1;
  //while ( (s_max > 0) && (doppler->fkdot[s_max] == 0) ) {
  //  s_max --;
  //}

  const REAL8 dt = xIn->deltaT;
  if ( insertPoint > numSamplesIn ) {
      XLALPrintError("ERROR! Would insert off end of matrix, numSamplesIn < insertPoint: %i, %i\n", numSamplesIn, insertPoint);
      XLAL_ERROR( XLAL_EINVAL );
  }
  /* original */ //UINT4 numSamplesIn  = xIn->data->length;
  /* resampling: use user input */

  const LIGOTimeGPS epoch = xIn->epoch;
  const REAL8 Dtau0 = GPSDIFF ( epoch, doppler->refTime );

  // loop over time samples
  if ( (indexEndResamp - indexStartResamp) > numSamplesIn){
      XLALPrintError("indexStartResamp, indexEndResamp: %u %u\n", indexStartResamp, indexEndResamp);
      XLALPrintError("ERROR! numSamplesIn only: %u\n", numSamplesIn);
      XLALPrintError("diff in indices, and num samples: %i %u\n", (indexEndResamp - indexStartResamp), numSamplesIn);
      XLAL_ERROR( XLAL_EINVAL );
  }

  REAL8 taup_j = 0.0;
  REAL8 cycles = 0.0;
  REAL4 cosphase, sinphase;
  COMPLEX8 em2piphase;
  for ( UINT4 j = 0; j < (indexEndResamp - indexStartResamp); j ++ )
    {
      /* j does appear to need indexStartResamp: count from overall start */
      taup_j = (j + indexStartResamp) * dt + Dtau0;
      //REAL8 Dtau_alpha_j = Dtau0 + taup_j;

      cycles = - freqShift * taup_j;

      //REAL8 Dtau_pow_kp1 = Dtau_alpha_j;
      /* Unused because we ignore spindown in CrossCorr */
      //for ( UINT4 k = 1; k <= s_max; k++ )
      //  {
      //    Dtau_pow_kp1 *= Dtau_alpha_j;
      //    cycles += - LAL_FACT_INV[k+1] * doppler->fkdot[k] * Dtau_pow_kp1;
      //  } // for k = 1 ... s_max

      XLAL_CHECK( XLALSinCos2PiLUT ( &sinphase, &cosphase, cycles ) == XLAL_SUCCESS, XLAL_EFUNC );
      em2piphase = crectf ( cosphase, sinphase );

      // weight the complex timeseries by the antenna patterns
      /* old comment above? I think it means by any freqshift */
      xOut[j+insertPoint] = em2piphase * xIn->data->data[(j+indexStartResamp)];

    } // for j < numSamplesIn

  return XLAL_SUCCESS;

} // XLALApplyCrossCorrFreqShiftResamp()

/** Compute the equivalent of Fa and Fb for CrossCorr's rho statistic */
static int
XLALComputeFaFb_CrossCorrResamp
  ( 
    ResampCrossCorrWorkspace                *restrict ws,                     /**< [out] contains modified Fa and Fb for cross-correlation */
    COMPLEX8                                *restrict wsFaX_k,                /**< [out] contains modified and normalized Fa for cross-correlation statistic */
    COMPLEX8                                *restrict wsFbX_k,                /**< [out] contains modified and normalized Fb for cross-correlation statistic */
    const MultiResampSFTPairMultiIndexList  *restrict resampMultiPairs,       /**< [in] resamp multi list of SFT pairs */
    const MultiCOMPLEX8TimeSeries           *restrict multiTimeSeries_SRC_a,  /**< [in] Resampled, heterodyned A(t)x(t) */
    const MultiCOMPLEX8TimeSeries           *restrict multiTimeSeries_SRC_b,  /**< [in] Resampled, heterodyned B(t)x(t) */
    const PulsarDopplerParams               *restrict dopplerpos,             /**< [in] Doppler point to search */
    const PulsarDopplerParams               *restrict binaryTemplateSpacings, /**< [in] spacings for the search */
    const REAL8                                       SRCsampPerTcoh,          /**< [in] floating point number of samples per tShort equivalent */
    const UINT4                                       detX,                   /**< [in] detector X index */
    const UINT4                                       sftK,                   /**< [in] tShort K index*/
    const UINT4                                       detY,                   /**< [in] detector Y index */
    const BOOLEAN                                     isL                     /**< [in] Indicates that is not K but L, the second half of the cross-corr pair */
  )
{
        XLAL_CHECK ( ws != NULL, XLAL_EINVAL );
        XLAL_CHECK ( wsFaX_k != NULL, XLAL_EINVAL );
        XLAL_CHECK ( wsFbX_k != NULL, XLAL_EINVAL );
        const REAL8 FreqOut0 = dopplerpos->fkdot[0];
        const REAL8 fHet = multiTimeSeries_SRC_a->data[0]->f0;
        const REAL8 dt_SRC = multiTimeSeries_SRC_b->data[0]->deltaT;
        const REAL8 dFreqMetric = binaryTemplateSpacings->fkdot[0];

        /* Compare to the F-statistic timing model, T1600531. Here, we cannot
         * change dt_SRC -- it is fixed and given to us. We cannot change it
         * because dt_SRC, after ceiling and power-of-two rounding, determined
         * by TspanFFT in SetupFstatResamp, and that TspanFFT must be = Tobs
         * for us to get the resampling data.
         *   Where the F-statistic handles power-of-two rounding by reducing 
         * dt_SRC, we must instead increase T_FFT by zero-padding.
         *   Yet the number of bins, numSamplesFFT, is fixed.
         *   What does this imply? Namely, we have a finer frequency resolution,
         * which we must use to evaluate cycles for times inside the FFT (as 
         * well as for bin shifts):
         *   dFreqMetric -> dFreqFFT
         * How do we find dFreqFFT? It is not dFreqMetric/decimateFFT now,
         * because T_FFT was lengthened. Whereas F-stat has a higher maximum
         * frequency of the band, we have a higher frequency resolution.
         * Nevertheless dFreqFFT remains equal to
         * dFreqFFT = 1/T_FFT = 1/(dt_SRC * numSamplesFFT), definitionally.
         * So if we defined a variable RedecimateFFT, by
         *   dFreqFFT := dFreqMetric/RedecimateFFT 
         * ->RedecimateFFT = dFreqMetric/(dFreqFFT)
         *                 = dFreqMetric * dt_SRC * numSamplesFFT
         * The finer frequency resolution needed by our FFT means that
         * RedecimateFFT is bigger than decimateFFT, by between 1x and 2x.
         * Thus, we must also change multiplications,
         *   x * decimateFFT -> (UINT4)floor( x * RedecimateFFT ) 
         * - - - Comparison to the timing model document - - - 
         * In the workspace function, we note that
         * our SRCsampPerTcoh corresponds to N^SRC_samp from the document.
         * However, our N^SRC_samp for a given Tcoh will be smaller than
         * The F-statistic, because for fixed number of FFT bins, dt_SRC
         * will be bigger than in the F-statistic (because it was not
         * reducable to fit the next power of 2). Correspondingly, our R
         * is smaller. Specifically, 
         *   R_crosscorr / R_fstat = D_fstat / D_crosscorr,
         * Where D_fstat is decimateFFT and D_crosscorr is RedecimateFFT.
         * This ratio is simply the power-of-2 rounding, i.e.,
         *   R_crosscorr / R_fstat = numSamplesFFT0 / numSamplesFFT,
         * - - - CrossCorr timing model by analogy - - -
         * Whereas the F-stat timing model had,
         *   tau^eff_F = tau^(0)_Fbin +
         *    N^FFT_samp / N_fbin * (
         *       R(tau^(0)_spin +b tau^(0)_bary) +5log_2(N^FFT_samp) tau^(0)_FFT
         *     ),
         * we have no search over spin-down and our b = 1 always.
         *   Our equivalent of tau^(0)_Fbin, however, is different, but the
         * more critical difference is the way that semicoherent segments are
         * combined. Namely, instead of multiplying 
         *   tau_semi_F = tau^eff_F by (Tobs/Tcoh) Ndet, [not in timing document]
         * we must calculate the above and multiply by the overlap factor,
         * (tCoh/tShort), after adding in any paired detectors. So we will have
         *   tau_semi_CC ~= tau^eff_CC (Tobs/Tshort) * (Ndet + (Ndet+1)Ndet/2 ) ,
         * where
         *   tau^eff_CC ~= tau^(0)_CCbin +
         *     N^FFT_samp / N_CCbin * (
         *         R_CC b tau^(0)_bary + 5log_2(N^FFT_samp) tau^(0)_FFT
         *       ).
         * More accurately, the Ndet term must go inside, because barycentering
         * occurs once per detector, the first FFT of a pair happens once per
         * detector, and the second FFT of a pair happens Triangular(Ndet),
         * where Triangular is the triangular number sequence,
         *    Triangular(Ndet) = (Ndet)(Ndet+1)/2.
         * Thus the full semicoherent number should be
         *   tau_semi_CC = (Tobs/Tshort) tau^(0)_CCbin + N^FFT_samp/N_CCbin *(
         *       (R_CC b tau^(0)_bary) * Ndet +
         *       (5log_2(N^FFT_samp) tau^(0)_FFT) * (Ndet + Ndet(Ndet+1)/2)
         *     ).
         * Accurately calculating N^FFT_samp = numSamplesFFT is therefore key:
         * following from the timing model document for the F-stat,
         *    numSamplesFFT0 = (Delta f_SFT) decimateFFT / dFreqMetric,
         * where,
         *    delta f_SFT = (1+4/(2 Dterms +1))(Delta f+Delta f_drift+16/Tsft)
         * with Delta_f = f_max - f_min,
         * for Sco X-1,
         *    Delta f_drift = v_earth / c + 2 pi f_max / (c Period)
         * Putting it altogether, where 
         *   N^FFT_samp = numSamplesFFT,
         *              = (UINT4)pow(2, ceil(log2(numSamplesFFT0))),
         *              = (UINT4)pow(2, ceil(log2(
         *                  (decimateFFT/dFreqMetric) *
         *                  (
         *                    (1+4/(2 Dterms + 1)) *
         *                    (  f_max - f_min + 16/Tsft +
         *                       v_earth/c +2 pi f_max /(c Period)
         *                    )
         *                  )
         *                ))).
         * That forms the basis of the CrossCorr timing model.
         * (GDM)
         */
        const REAL4 RedecimateFFT = dFreqMetric * dt_SRC * (ws->numSamplesFFT); //Equals metric freq res divided by freq res of our FFT
        const REAL8 dFreqFFT = dFreqMetric / RedecimateFFT;
        const REAL8 freqShiftInFFT  = remainder(FreqOut0-fHet,dFreqFFT); // frequency shift to closest bin in final resolution
        const REAL8 fMinFFT = fHet + freqShiftInFFT - dFreqFFT * (ws->numSamplesFFT/2);  // we'll shift DC into the *middle bin* N/2  [N always even!]
        XLAL_CHECK ( FreqOut0 >= fMinFFT, XLAL_EDOM, "Lowest output frequency outside the available frequency band: [FreqOut0 = %.16g] < [fMinFFT = %.16g]\n", FreqOut0, fMinFFT );
        /* -> THIS AFFECTS indices below, in FaX_k loop */
        const UINT4 offset_bins = (UINT4) lround ( ( FreqOut0 - fMinFFT ) / dFreqFFT );
        const UINT4 maxOutputBin = offset_bins + (UINT4)floor( (ws->numFreqBinsOut-1) * RedecimateFFT );
        XLAL_CHECK ( maxOutputBin < ws->numSamplesFFT, XLAL_EDOM, "Highest output frequency bin outside available band: [maxOutputBin = %d] >= [ws->numSamplesFFT = %d]\n", maxOutputBin, ws->numSamplesFFT );

        /* start and stop variables */
        UINT4 startFirstInd = 0;
        UINT4 startInd = 0;
        UINT4 endInd = 0;
        UINT4 sftIndFirst = 0;
        UINT4 detInd = 0;
        UINT4 sftInd = 0;

        /* set sftInd and dftInd values using the like-named fields encoded in
         * resampMultiPairs, which are with respect to the sft vector,
         * inputSFTs. This should handle gaps correctly.
         * For any loops over derived values, determined only for matching sfts,
         * instead use detX, sftK, detY, sftL, otherwise index may be out-of-bound 
         */
        const ResampSFTPairMultiIndexList  *restrict resampMultiPairsDetX = &(resampMultiPairs->data[detX]);
        const ResampSFTMultiIndexList  *restrict resampMultiPairsDetXsftK = &(resampMultiPairsDetX->data[sftK]);

        if (isL == TRUE){
          detInd = resampMultiPairsDetXsftK->data[detY].detInd;
        } else {
          detInd = resampMultiPairsDetX->detInd;
        }

        /* ------- begin FFT replacement using RESAMP ------- */
        /* resamp data */
        const COMPLEX8TimeSeries  *restrict resampDataArrayA = multiTimeSeries_SRC_a->data[detInd];
        const COMPLEX8TimeSeries  *restrict resampDataArrayB = multiTimeSeries_SRC_b->data[detInd];

        if (isL == TRUE){
          /* Need to consider handling if no sftInd exists because none match */
          if (resampMultiPairsDetXsftK->data[detY].length > 0){
            sftIndFirst = resampMultiPairsDetXsftK->data[detY].data[0].sftInd;
            startFirstInd = (UINT4)round(sftIndFirst * SRCsampPerTcoh);
          }
        } else {
            sftIndFirst = resampMultiPairsDetXsftK->sftInd;
            startFirstInd = (UINT4)round(sftIndFirst * SRCsampPerTcoh);
        }
        if (startFirstInd > resampDataArrayA->data->length){
            startFirstInd = resampDataArrayA->data->length;
        }

        memset ( ws->TS_FFT, 0, ws->numSamplesFFT * sizeof(ws->TS_FFT[0]) );
        UINT4 sftLength = 0;
        if (isL == TRUE){
            sftLength = resampMultiPairsDetXsftK->data[detY].length;
        } else {
            sftLength = 1;
        }   
        // Load and FFT A time series
        for (UINT4 sft=0; sft < sftLength; sft++){
          //if (resampMultiPairsDetXsftK->data[detY].data[sft].sciFlag > 0){
            /* Needs to be defined wrt resamp time series, so:*/
            if (isL == TRUE){
              sftInd = resampMultiPairsDetXsftK->data[detY].data[sft].sftInd;
            } else {
              sftInd = resampMultiPairsDetXsftK->sftInd;
            }
            startInd = (UINT4)round(sftInd * SRCsampPerTcoh);
            endInd = (UINT4)round((sftInd+1) * SRCsampPerTcoh);
            if (startInd > resampDataArrayA->data->length){
                startInd = resampDataArrayA->data->length;
            }
            if (endInd > resampDataArrayA->data->length){
                endInd = resampDataArrayA->data->length;
            }
            if (startInd > endInd){
                startInd = endInd;
            }
            UINT4 headOfSliceIndex = startInd - startFirstInd;
            XLAL_CHECK ( XLALApplyCrossCorrFreqShiftResamp ( ws->TS_FFT, resampDataArrayA, dopplerpos, freqShiftInFFT, startInd, endInd, ws->numSamplesFFT, headOfSliceIndex ) == XLAL_SUCCESS, XLAL_EFUNC );
          //}
        }
        fftwf_execute ( ws->fftplan );
        for ( UINT4 k = 0; k < ws->numFreqBinsOut; k++ ) {
          ws->FaX_k[k] = ws->FabX_Raw [ offset_bins + (UINT4)floor(k * RedecimateFFT)  ];
        }
        // END load and FFT A time series

        // Load and FFT B time series
        memset ( ws->TS_FFT, 0, ws->numSamplesFFT * sizeof(ws->TS_FFT[0]) );
        for (UINT4 sft=0; sft < sftLength; sft++){
          //if (resampMultiPairsDetXsftK->data[detY].data[sft].sciFlag > 0){
            if (isL == TRUE){
              sftInd = resampMultiPairsDetXsftK->data[detY].data[sft].sftInd;
            } else {
              sftInd = resampMultiPairsDetXsftK->sftInd;
            }
            startInd = (UINT4)round(sftInd * SRCsampPerTcoh);
            endInd = (UINT4)round( (sftInd+1) * SRCsampPerTcoh);
            if (startInd > resampDataArrayB->data->length){
                startInd = resampDataArrayB->data->length;
            }
            if (endInd > resampDataArrayB->data->length){
                endInd = resampDataArrayB->data->length;
            }
            if (startInd > endInd){
                startInd = endInd;
            }
            UINT4 headOfSliceIndex = startInd - startFirstInd;
            XLAL_CHECK ( XLALApplyCrossCorrFreqShiftResamp ( ws->TS_FFT, resampDataArrayB, dopplerpos, freqShiftInFFT, startInd, endInd, ws->numSamplesFFT, headOfSliceIndex ) == XLAL_SUCCESS, XLAL_EFUNC );
          //}
        }
        fftwf_execute ( ws->fftplan );
        for ( UINT4 k = 0; k < ws->numFreqBinsOut; k++ ) {
          ws->FbX_k[k] = ws->FabX_Raw [ offset_bins + (UINT4)floor(k * RedecimateFFT)  ];
        }
        // End load and FFT B time series 

        const REAL8 dtauX = GPSDIFF ( resampDataArrayB->epoch, dopplerpos->refTime );
        const REAL8 timeFromResampStartToSFTFirst = startFirstInd * dt_SRC;
        REAL8 f_kOut;
        REAL8 f_kIn;
        REAL8 cyclesOut;
        REAL8 cyclesIn;
        REAL8 cycles;
        REAL4 sinphase, cosphase;
        COMPLEX8 normX_k;
        for ( UINT4 k = 0; k < ws->numFreqBinsOut; k++ )
        {
          f_kOut = FreqOut0 + k * dFreqMetric;
          f_kIn = (offset_bins + (UINT4)floor(k*RedecimateFFT)) * dFreqFFT;
          cyclesOut = - f_kOut * dtauX;
          cyclesIn = -f_kIn * timeFromResampStartToSFTFirst;
          cycles = cyclesOut + cyclesIn;
          XLALSinCos2PiLUT ( &sinphase, &cosphase, cycles );
          normX_k = dt_SRC * crectf ( cosphase, sinphase );
          wsFaX_k[k] = normX_k * ws->FaX_k[k];
          wsFbX_k[k] = normX_k * ws->FbX_k[k];
        } // for k < numFreqBinsOut

        return XLAL_SUCCESS;
} // XLALComputeFaFb_CrossCorrResamp()
