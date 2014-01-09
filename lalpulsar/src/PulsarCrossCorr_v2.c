/*
 *  Copyright (C) 2012, 2013 John Whelan, Shane Larson and Badri Krishnan
 *  Copyright (C) 2013, 2014 Badri Krishnan, John Whelan, Yuanhao Zhang
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

/** Calculate the Doppler-shifted frequency associated with each SFT in a list */
/* This is according to Eqns 2.11 and 2.12 of Dhurandhar et al 2008 */
/* Also returns the signal phase according to eqn 2.4 */
int XLALGetDopplerShiftedFrequencyInfo
  (
   REAL8Vector         *shiftedFreqs, /**< Output list of shifted frequencies */
   UINT4Vector         *lowestBins,   /**< Output list of bin indices */
   REAL8Vector         *kappaValues,  /**< Output list of bin offsets */
   REAL8Vector         *signalPhases, /**< Output list of signal phases */
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
  REAL8 timeDiff, factor, fhat, phiByTwoPi;
  SFTIndex sftInd;
  SSBtimes *times;

  numSFTs = sfts->length;
  if ( signalPhases->length !=numSFTs
       || shiftedFreqs->length !=numSFTs
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
    phiByTwoPi = fmod ( fhat * timeDiff , 1.0 );
    factor = timeDiff;
    for (k = 1;  k < PULSAR_MAX_SPINS; k++) {
      fhat += dopp->fkdot[k] * factor;
      factor *= timeDiff / (k+1);
      phiByTwoPi += dopp->fkdot[k] * factor;
    }
    signalPhases->data[indI] = LAL_TWOPI * fmod ( phiByTwoPi , 1.0 );
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
/* Note this is probably obsolete because the noise-weighting takes
   care of these factors */
int XLALCalculateCrossCorrSigmaUnshifted
  (
   REAL8Vector      **sigma_alpha,    /**< Output: vector of sigma_alpha values */
   SFTPairIndexList  *pairIndexList,  /**< Input: list of SFT pairs */
   SFTIndexList      *indexList,      /**< Input: list of SFTs */
   MultiPSDVector    *psds,           /**< Input: PSD estimate (Sn*Tsft/2) for each SFT */
   REAL8              freq,           /**< Frequency to extract from PSD */
   REAL8              Tsft            /**< SFT duration */
  )
{

  UINT8 j, numPairs, numSFTs, freqInd;
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
    freqInd = (UINT8) lrintf( (freq - psd->f0)/psd->deltaF );
    psdData->data[j] = psd->data->data[freqInd];
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

/** Construct vector of G_alpha amplitudes for each SFT pair */
/* This is averaged over unknown cosi and psi */
/* Allocates memory as well */
int XLALCalculateAveCurlyGAmpUnshifted
  (
   REAL8Vector      **G_alpha,       /* Output: vector of sigma_alpha values */
   SFTPairIndexList  *pairIndexList, /* Input: list of SFT pairs */
   SFTIndexList      *indexList,     /* Input: list of SFTs */
   MultiAMCoeffs     *multiCoeffs    /* Input: AM coefficients */
  )
{

  UINT8 j, numPairs;
  UINT8 detInd1, detInd2;
  UINT8 sftInd1, sftInd2;
  UINT8 sftNum1, sftNum2;
  REAL8Vector *ret = NULL;

  numPairs = pairIndexList->length;

  XLAL_CHECK ( ( ret = XLALCreateREAL8Vector ( numPairs ) ) != NULL, XLAL_EFUNC, "XLALCreateREAL8Vector ( %d ) failed.", numPairs );

  for (j=0; j < numPairs; j++) {
    sftNum1 = pairIndexList->data[j].sftNum[0];
    sftNum2 = pairIndexList->data[j].sftNum[1];
    detInd1 = indexList->data[sftNum1].detInd;
    detInd2 = indexList->data[sftNum2].detInd;
    sftInd1 = indexList->data[sftNum1].sftInd;
    sftInd2 = indexList->data[sftNum2].sftInd;
    ret->data[j] = 0.1 * ( multiCoeffs->data[detInd1]->a->data[sftInd1]
			   * multiCoeffs->data[detInd2]->a->data[sftInd2]
			   + multiCoeffs->data[detInd1]->b->data[sftInd1]
			   * multiCoeffs->data[detInd2]->b->data[sftInd2] );
  }

  (*G_alpha) = ret;
  return XLAL_SUCCESS;
}

/** Calculate multi-bin cross-correlation statistic */
/* This assumes rectangular or nearly-rectangular windowing */
int XLALCalculatePulsarCrossCorrStatistic
(
 REAL8              *ccStat,   /* Output: cross-correlation statistic rho */
 REAL8           *evSquared,   /* Output: (E[rho]/h0^2)^2 */
 REAL8Vector     *curlyGAmp,   /* Input: Amplitude of curly G for each pair */
 REAL8Vector  *signalPhases,   /* Input: Phase of signal for each SFT */
 UINT4Vector    *lowestBins,   /* Input: Bin index to start with for each SFT */
 REAL8Vector   *kappaValues,   /* Input: Fractional offset of signal freq from best bin center */
 UINT4              numBins,   /* Input: Number of bins to include in calc */
 SFTPairIndexList *sftPairs,   /* Input: flat list of SFT pairs */
 SFTIndexList   *sftIndices,   /* Input: flat list of SFTs */
 MultiSFTVector  *inputSFTs    /* Input: SFT data */
 )
{

  UINT8 numSFTs = sftIndices->length;
  if ( signalPhases->length !=numSFTs
       || lowestBins->length !=numSFTs
       || kappaValues->length !=numSFTs ) {
    XLALPrintError("Lengths of SFT-indexed lists don't match!");
    XLAL_ERROR(XLAL_EBADLEN );
  }

  UINT8 numPairs = sftPairs->length;
  if ( curlyGAmp->length !=numPairs ) {
    XLALPrintError("Lengths of pair-indexed lists don't match!");
    XLAL_ERROR(XLAL_EBADLEN );
  }

  *ccStat = 0.0;
  *evSquared = 0.0;
  for (UINT8 alpha=0; alpha < numPairs; alpha++) {
    UINT8 sftNum1 = sftPairs->data[alpha].sftNum[0];
    UINT8 sftNum2 = sftPairs->data[alpha].sftNum[1];
    UINT8 detInd1 = sftIndices->data[sftNum1].detInd;
    UINT8 detInd2 = sftIndices->data[sftNum2].detInd;
    UINT8 sftInd1 = sftIndices->data[sftNum1].sftInd;
    UINT8 sftInd2 = sftIndices->data[sftNum2].sftInd;
    COMPLEX8 *dataArray1 = inputSFTs->data[detInd1]->data[sftInd1].data->data;
    COMPLEX8 *dataArray2 = inputSFTs->data[detInd2]->data[sftInd2].data->data;
    COMPLEX16 GalphaCC = curlyGAmp->data[alpha]
      * cexp( I * ( signalPhases->data[sftNum1]
		   - signalPhases->data[sftNum2] )
	      );
    UINT4 baseCCSign = 1; /* Alternating sign is (-1)**(k1-k2) */
    if ( ( (lowestBins->data[sftNum1]-lowestBins->data[sftNum2]) % 2) != 0 ) {
      baseCCSign = -1;
    }

    for (UINT8 j=0; j < numBins; j++) {
      COMPLEX16 data1 = dataArray1[lowestBins->data[sftNum1]+j];
      REAL8 sincFactor = gsl_sf_sinc(kappaValues->data[lowestBins->data[sftNum1]+j]);
      /* Normalized sinc, i.e., sin(pi*x)/(pi*x) */
      UINT8 ccSign = baseCCSign;
      for (UINT8 k=0; k < numBins; k++) {
	COMPLEX16 data2 = dataArray2[lowestBins->data[sftNum2]+k];
	sincFactor *= gsl_sf_sinc(kappaValues->data[lowestBins->data[sftNum2]+k]);
	*ccStat += crealf ( GalphaCC * ccSign * sincFactor
			    * conj(data1) * data2 );
	REAL8 GalphaAmp = curlyGAmp->data[alpha] * sincFactor;
	*evSquared += SQUARE( GalphaAmp );
	ccSign *= -1;
      }
      baseCCSign *= -1;
    }
  }
  return XLAL_SUCCESS;
}
/*calculate the weighted factors, and metric diagnol components, also include the estimation of sensitivity E[rho]/(h_0)^2*/
int XLALCalculateMetricElements
  (
   REAL8             *TSquaWeightedAve, /*Output: weighted factors*/
   REAL8             *SinSquaWeightedAve,  
   REAL8             *hSens,            /*Output:sensitivity*/
   REAL8             *g_ff,             /*Output:metric elements*/
   REAL8             *g_aa, 
   REAL8             *g_TT, 
   REAL8Vector       *G_alpha,       /* Input: vector of sigma_alpha values */ 
   SFTPairIndexList  *pairIndexList, /* Input: list of SFT pairs */
   SFTIndexList      *indexList,     /* Input: list of SFTs */   
   MultiSFTVector    *sfts,          /* Input: set of per-detector SFT vectors */
   REAL8             pOrb,           /* Input: orbit period in second*/
   REAL8             aPro,           /* Input: projected semimajor*/
   REAL8             f              /* Input: frequency*/ 
   /*REAL8             *devTsq,     */      /*Output: mean time deviation^2*/
   /*REAL8             *g_pp,*/
   )
{
  UINT8 sftNum1=0;
  UINT8 sftNum2=0;
  UINT8 j=0;
  REAL8 T=0;
  REAL8 denom=0;
  REAL8 sinSquare=0;
  REAL8 tSquare=0;
  REAL8 rhosum=0;
  LIGOTimeGPS *T1=NULL;
  LIGOTimeGPS *T2=NULL;
   /*  REAL8 hfT=0;
      REAL8 Tmean=0;
      REAL8 muT=0;
      REAL8 sumDev=0;
      UINT8 k=0*/

   UINT8 numalpha = G_alpha->length;

  for (j=0; j < numalpha; j++) {
    REAL8 sqrG_alpha=SQUARE(G_alpha->data[j]);
    sftNum1 = pairIndexList->data[j].sftNum[0];
    sftNum2 = pairIndexList->data[j].sftNum[1];
    UINT8 detInd1 = indexList->data[sftNum1].detInd;
    UINT8 detInd2 = indexList->data[sftNum2].detInd;
    UINT8 sftInd1 = indexList->data[sftNum1].sftInd;
    UINT8 sftInd2 = indexList->data[sftNum2].sftInd;
    T1 = &(sfts->data[detInd1]->data[sftInd1].epoch);
    T2 = &(sfts->data[detInd2]->data[sftInd2].epoch);
    T = XLALGPSDiff( T1, T2 );
    sinSquare += sqrG_alpha*SQUARE(sin(LAL_PI*T/pOrb));/*(G_alpha)^2*(sin(\pi*T/T_orbit))^2*/
    tSquare += sqrG_alpha*SQUARE( G_alpha->data[j]*T); /*(\curlyg_alpha*)^2*T^2*/
    denom += sqrG_alpha;                               /*calculate the denominator*/
    rhosum += 2*sqrG_alpha;
    /*hfT=0.5*T;
      Tmean=XLALGPSAdd(&T2, hfT);*/
    /*muT +=Tmean/numalpha;*/                            /*calculate the average of Tmean*/
      }

  *TSquaWeightedAve =(tSquare/denom);
  *SinSquaWeightedAve =(sinSquare/denom);
  *hSens = sqrt(rhosum);
  *g_ff= *TSquaWeightedAve* 2 * SQUARE(LAL_PI);
  *g_aa= *SinSquaWeightedAve* SQUARE(LAL_PI*f);
  *g_TT= *SinSquaWeightedAve* SQUARE(2*SQUARE(LAL_PI)*f*aPro/pOrb);
  return XLAL_SUCCESS;

  /* *g_pp=SQUARE(2*SQUARE(LAL_PI)*f*aPro/SQUARE(pOrb))*devTsq*SinSquaWeightedAve;*/
  /*for(k=0;k < numalpha;k++){
    sftNum1 = pairIndexList->data[k].sftNum[0];
    sftNum2 = pairIndexList->data[k].sftNum[1];
    T1 = sfts->data[indexList->data[sftNum1].detInd]->data[indexList->data[sftNum1].sftInd].epoch;
    T2 = sfts->data[indexList->data[sftNum2].detInd]->data[indexList->data[sftNum2].sftInd].epoch;
    Tmean=XLALGPSAdd(&T2, hfT);*/
  /*sumDev +=SQUARE (G_alpha->data[k]*(Tmean-muT));  */      /*calculate the mean time deviation squared*/
    /* }  */
  /**devTsq=sumDev/denom;*/
}

/*int XLALCalculateMetricElements
( 
   REAL8             *g_ff, Output:metric elements
   REAL8             *g_aa, 
   REAL8             *g_TT, 
   REAL8             *g_pp,
   REAL8             aPro,            Input: variables
   REAL8             f,
   REAL8             pOrb,
   REAL8             devTsq,           Input: T deviation^2 function(5.25)
   REAL8             TSquaWeightedAve, Input: weighted factors
   REAL8             SinSquaWeightedAve 
    )
{
  *g_ff=2*SQUARE(LAL_PI)*TSquaWeightedAve;
  *g_aa=SQUARE(LAL_PI*f)*SinSquaWeightedAve;
  *g_TT=SQUARE(2*SQUARE(LAL_PI)*f*aPro/pOrb)*SinSquaWeightedAve;
  *g_pp=SQUARE(2*SQUARE(LAL_PI)*f*aPro/SQUARE(pOrb))*devTsq*SinSquaWeightedAve;
  return XLAL_SUCCESS;
}
*/
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
