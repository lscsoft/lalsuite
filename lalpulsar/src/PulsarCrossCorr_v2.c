/*
 *  Copyright (C) 2012, 2013 John Whelan, Shane Larson and Badri Krishnan
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

/** Calculate the Doppler-shifted frequency associate with each SFT in a list */
/* This is according to Eqns 2.11 and 2.12 of Dhurandhar et al 2008 */
int XLALGetDopplerShiftedFrequencyInfo
  (
   REAL8Vector         *shiftedFreqs, /**< Output list of shifted frequencies */
   UINT4Vector         *lowestBins,   /**< Output list of bin indices */
   REAL8Vector         *kappaValues,  /**< Output list of bin offsets */
   UINT4               numBins,       /**< Number of frequency bins to use */
   PulsarDopplerParams *dopp,         /**< Doppler parameters for signal */
   SFTIndexList        *sfts,         /**< List of indices for SFTs */
   MultiSSBtimes       *multiTimes,   /**< SSB or Binary times */
   REAL8               Tsft           /**< SFT duration */
  )
{
  UINT8 numSFTs;
  UINT8 indI;
  UINT4 k;
  REAL8 timeDiff, factor, fhat;
  SFTIndex sftInd;
  SSBtimes *times;

  numSFTs = sfts->length;
  if ( shiftedFreqs->length !=numSFTs
       || lowestBins->length !=numSFTs
       || kappaValues->length !=numSFTs ) {
    XLALPrintError("Lengths of SFT-indexed lists don't match!");
    XLAL_ERROR(XLAL_EBADLEN );
  }

  if ( numBins < 1 ) {
    XLALPrintError("Must specify a positive number of bins to use!");
    XLAL_ERROR(XLAL_EBADLEN );
  }

  /* now calculate the intrinsic signal frequency in the SFT */
  /* fhat = f_0 + f_1(t-t0) + f_2(t-t0)^2/2 + ... */

  /* this is the sft reference time  - the pulsar reference time */
  for (indI=0; indI < numSFTs; indI++) {
    sftInd = sfts->data[indI];
    times = multiTimes->data[sftInd.detInd];
    timeDiff = times->DeltaT->data[sftInd.sftInd]
      + XLALGPSDiff( &(times->refTime), &(dopp->refTime));
    fhat = dopp->fkdot[0]; /* initialization */
    factor = 1.0;
    for (k = 1;  k < PULSAR_MAX_SPINS; k++) {
      factor *= timeDiff / k;
      fhat += dopp->fkdot[k] * factor;
    }
    shiftedFreqs->data[indI] = fhat * times->Tdot->data[sftInd.sftInd];
    lowestBins->data[indI]
      = ceil(shiftedFreqs->data[indI] * Tsft - 0.5*numBins);
    kappaValues->data[indI] = lowestBins->data[indI]
      - shiftedFreqs->data[indI] * Tsft;
  }

  return XLAL_SUCCESS;

}

/** Construct flat SFTIndexList out of a MultiSFTVector */
/* Allocates memory as well */
int XLALCreateSFTIndexListFromMultiSFTVect
  (
   SFTIndexList        **indexList,   /* Output: flat list of indices to locate SFTs */
   MultiSFTVector      *sfts         /* Input: set of per-detector SFT vectors */
  )
{
  SFTIndexList *ret = NULL;
  UINT8 numSFTs;
  UINT8 j, k, l, numForDet;
  UINT4 numDets;

  numDets = sfts->length;

  numSFTs = 0;
  for (k=0; k < numDets; k++) {
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

  j = 0;
  for (k=0; k < numDets; k++) {
    numForDet = sfts->data[k]->length;
    for (l=0; l < numForDet; l++) {
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
   SFTPairIndexList  **pairIndexList,  /* Output: list of SFT pairs */
   SFTIndexList       *indexList,      /* Input: list of indices to locate SFTs */
   MultiSFTVector     *sfts,           /* Input: set of per-detector SFT vectors */
   REAL8               maxLag,         /* Maximum allowed lag time */
   BOOLEAN             inclAutoCorr    /* Flag indicating whether a "pair" of an SFT with itself is allowed */
  )
{
  SFTPairIndexList *ret = NULL;
  UINT8 numSFTs;
  UINT8 numPairs;
  UINT8 j, k, l, lMin;
  REAL8 timeDiff;
  LIGOTimeGPS gps1, gps2;

  numSFTs = indexList->length;

  if ( ( ret = XLALCalloc( 1, sizeof( *ret ) )) == NULL ) {
    XLAL_ERROR ( XLAL_ENOMEM );
  }

  /* maximum possible number of pairs */

  if ( inclAutoCorr ) {
    numPairs = numSFTs*(numSFTs+1)/2;
  } else {
    numPairs = numSFTs*(numSFTs-1)/2;
  }
  ret->length = numPairs;
  if ( ( ret->data = XLALCalloc ( numPairs, sizeof ( *ret->data ) )) == NULL ) {
    XLALFree ( ret );
    XLAL_ERROR ( XLAL_ENOMEM );
  }

  j = 0;
  for (k=0; k < numSFTs; k++) {
    if ( inclAutoCorr ) {
      lMin = k;
    } else {
      lMin = k+1;
    }
    gps1 = sfts->data[indexList->data[k].detInd]->data[indexList->data[k].sftInd].epoch;
    for (l=lMin; l < numSFTs; l++) {
      gps2 = sfts->data[indexList->data[l].detInd]->data[indexList->data[l].sftInd].epoch;
      timeDiff = XLALGPSDiff(&gps1,&gps2);
      if (abs(timeDiff) <= maxLag) {
	ret->data[j].sftNum[0] = k;
	ret->data[j].sftNum[1] = l;
	++j;
      }
    }
  }
  ret->length = j;
  if ( ( ret->data = XLALRealloc ( ret->data, j * sizeof ( *ret->data ) )) == NULL ) {
    XLALFree ( ret );
    XLAL_ERROR ( XLAL_ENOMEM );
  }
  

  (*pairIndexList) = ret;
  
  return XLAL_SUCCESS;
}

/** Construct vector of sigma_alpha values for each SFT pair */
/* This version uses a single frequency rather than the doppler-shifted ones */
/* Allocates memory as well */
int XLALCalculateCrossCorrSigmaUnshifted
  (
   REAL8Vector      **sigma_alpha,    /* Output: vector of sigma_alpha values */
   SFTPairIndexList  *pairIndexList,  /* Input: list of SFT pairs */
   SFTIndexList      *indexList,      /* Input: list of SFTs */
   MultiPSDVector    *psds,           /* Input: PSD estimate (Sn*Tsft/2) for each SFT */
   REAL8              freq,           /* Frequency to extract from PSD */
   REAL8              Tsft            /**< SFT duration */
  )
{

  UINT8 j, numPairs, numSFTs, index;
  REAL8Vector *psdData;
  REAL8FrequencySeries *psd;
  SFTIndex sftIndex;
  REAL8Vector *ret = NULL;
  REAL8 Tsft4;

  Tsft4 = SQUARE(SQUARE(Tsft));

  numPairs = pairIndexList->length;
  numSFTs = indexList->length;

  XLAL_CHECK ( ( psdData = XLALCreateREAL8Vector ( numSFTs ) ) != NULL, XLAL_EFUNC, "XLALCreateREAL8Vector ( %d ) failed.", numSFTs );

  for (j=0; j < numSFTs; j++) {
    sftIndex = indexList->data[j];
    psd = &(psds->data[sftIndex.detInd]->data[sftIndex.sftInd]);
    index = (UINT8) floor( (freq - psd->f0)/psd->deltaF + 0.5 );
    psdData->data[j] = psd->data->data[index];
  }

  XLAL_CHECK ( ( ret = XLALCreateREAL8Vector ( numPairs ) ) != NULL, XLAL_EFUNC, "XLALCreateREAL8Vector ( %d ) failed.", numPairs );

  for (j=0; j < numPairs; j++) {
    ret->data[j] = ( psdData->data[pairIndexList->data[j].sftNum[0]]
		     * psdData->data[pairIndexList->data[j].sftNum[1]]
		     ) / Tsft4 ;
  }

  (*sigma_alpha) = ret;
  XLALDestroyREAL8Vector ( psdData );
  return XLAL_SUCCESS;
}
