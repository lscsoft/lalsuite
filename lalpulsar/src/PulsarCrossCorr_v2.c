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
   REAL8Vector        *shiftedFreqs, /**< Output list of shifted frequencies */
   UINT4Vector          *lowestBins, /**< Output list of bin indices */
   REAL8Vector         *kappaValues, /**< Output list of bin offsets */
   COMPLEX8Vector  *expSignalPhases, /**< Output list of signal phases */
   REAL8VectorSequence    *sincList, /**< Output list of sinc factors */
   UINT4                    numBins, /**< Number of frequency bins to use */
   PulsarDopplerParams        *dopp, /**< Doppler parameters for signal */
   SFTIndexList         *sftIndices, /**< List of indices for SFTs */
   MultiSFTVector        *inputSFTs, /**< SFT data (needed for f0) */
   MultiSSBtimes        *multiTimes, /**< SSB or Binary times */
   REAL8                       Tsft  /**< SFT duration */
  )
{
  UINT8 numSFTs;
  UINT4 k;
  REAL8 timeDiff, factor, fhat, phiByTwoPi;

  numSFTs = sftIndices->length;
  if ( expSignalPhases->length !=numSFTs
       || shiftedFreqs->length !=numSFTs
       || lowestBins->length !=numSFTs
       || kappaValues->length !=numSFTs
       || sincList->length !=numSFTs ) {
    XLALPrintError("Lengths of SFT-indexed lists don't match!");
    XLAL_ERROR(XLAL_EBADLEN );
  }

  XLAL_CHECK ( ( inputSFTs->length == multiTimes->length ),
	       XLAL_EBADLEN,
	       "Lengths of detector-indexed lists don't match!" );

  if ( numBins < 1 ) {
    XLALPrintError("Must specify a positive number of bins to use!");
    XLAL_ERROR(XLAL_EBADLEN );
  }

  /* now calculate the intrinsic signal frequency in the SFT */
  /* fhat = f_0 + f_1(t-t0) + f_2(t-t0)^2/2 + ... */

  /* this is the sft reference time  - the pulsar reference time */
  for (UINT8 sftNum=0; sftNum < numSFTs; sftNum++) {
    UINT8 detInd = sftIndices->data[sftNum].detInd;
    XLAL_CHECK ( ( detInd < inputSFTs->length ),
		 XLAL_EINVAL,
		 "SFT asked for detector index off end of list:\n sftNum=%d, detInd=%d, inputSFTs->length=%d\n",
		 sftNum, detInd, inputSFTs->length );

    UINT8 numSFTsDet = inputSFTs->data[detInd]->length;
    SSBtimes *times;
    times = multiTimes->data[detInd];
    XLAL_CHECK ( ( times->DeltaT->length == numSFTsDet )
		 && ( times->Tdot->length == numSFTsDet ),
		 XLAL_EBADLEN,
		 "Lengths of multilists don't match!" );

    UINT8 sftInd = sftIndices->data[sftNum].sftInd;
    XLAL_CHECK ( ( sftInd < numSFTsDet ),
		 XLAL_EINVAL,
		 "SFT asked for SFT index off end of list:\n sftNum=%d, detInd=%d, sftInd=%d, numSFTsDet=%d\n",
		 sftNum, detInd, sftInd, numSFTsDet );
    timeDiff = times->DeltaT->data[sftInd]
      + XLALGPSDiff( &(times->refTime), &(dopp->refTime));
    fhat = dopp->fkdot[0]; /* initialization */
    phiByTwoPi = fmod ( fhat * timeDiff , 1.0 );
    factor = timeDiff;
    for (k = 1;  k < PULSAR_MAX_SPINS; k++) {
      fhat += dopp->fkdot[k] * factor;
      factor *= timeDiff / (k+1);
      phiByTwoPi += dopp->fkdot[k] * factor;
    }
    expSignalPhases->data[sftNum] = cexp(I * LAL_TWOPI * fmod ( phiByTwoPi , 1.0 ));
    shiftedFreqs->data[sftNum] = fhat * times->Tdot->data[sftInd];
    REAL8 fminusf0 = shiftedFreqs->data[sftNum] - inputSFTs->data[detInd]->data[sftInd].f0;
    lowestBins->data[sftNum]
      = ceil( fminusf0 * Tsft - 0.5*numBins );
    kappaValues->data[sftNum] = lowestBins->data[sftNum] - fminusf0 * Tsft;

    for (UINT8 l=0; l < numBins; l++) {
      sincList->data[sftNum*numBins + l] = gsl_sf_sinc(kappaValues->data[sftNum]+l);
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
  UINT8 j, k, l, lMin;
  REAL8 timeDiff;
  LIGOTimeGPS gps1, gps2;

  numSFTs = indexList->length;

  if ( ( ret = XLALCalloc( 1, sizeof( *ret ) )) == NULL ) {
    XLAL_ERROR ( XLAL_ENOMEM );
  }

  /* do two passes, one to count the number of pairs so the list can be allocated, and one to actually populate the list. */

  /* maximum possible number of pairs */

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
 COMPLEX8Vector  *expSignalPhases,   /* Input: Phase of signal for each SFT */
 UINT4Vector    *lowestBins,   /* Input: Bin index to start with for each SFT */
 REAL8Vector   *kappaValues,   /* Input: Fractional offset of signal freq from best bin center */
 REAL8VectorSequence  *sincList,   /* Input: input the sinc factors*/
 SFTPairIndexList *sftPairs,   /* Input: flat list of SFT pairs */
 SFTIndexList   *sftIndices,   /* Input: flat list of SFTs */
 MultiSFTVector  *inputSFTs,   /* Input: SFT data */
 MultiNoiseWeights *multiWeights,  /* Input: nomalizeation factor S^-1 & weights for each SFT */
 UINT4              numBins   /**Input Number of frequency bins to be taken into calc */
 )
{

  UINT8 numSFTs = sftIndices->length;
  if ( expSignalPhases->length !=numSFTs
       || lowestBins->length !=numSFTs
       || kappaValues->length !=numSFTs
       || sincList->length !=numSFTs ) {
    XLALPrintError("Lengths of SFT-indexed lists don't match!");
    XLAL_ERROR(XLAL_EBADLEN );
  }

  UINT8 numPairs = sftPairs->length;
  if ( curlyGAmp->length !=numPairs ) {
    XLALPrintError("Lengths of pair-indexed lists don't match!");
    XLAL_ERROR(XLAL_EBADLEN );
  }
  REAL8 nume=0;
  REAL8 curlyGSqr=0;
  *ccStat = 0.0;
  *evSquared = 0.0;
  for (UINT8 alpha=0; alpha < numPairs; alpha++) {
    UINT8 sftNum1 = sftPairs->data[alpha].sftNum[0];
    UINT8 sftNum2 = sftPairs->data[alpha].sftNum[1];

    XLAL_CHECK ( ( sftNum1 < numSFTs ) && ( sftNum2 < numSFTs ),
		 XLAL_EINVAL,
		 "SFT pair asked for SFT index off end of list:\n alpha=%d, sftNum1=%d, sftNum2=%d, numSFTs=%d\n",
		 alpha,  sftNum1, sftNum2, numSFTs );

    UINT8 detInd1 = sftIndices->data[sftNum1].detInd;
    UINT8 detInd2 = sftIndices->data[sftNum2].detInd;

    XLAL_CHECK ( ( detInd1 < inputSFTs->length )
		 && ( detInd2 < inputSFTs->length ),
		 XLAL_EINVAL,
		 "SFT asked for detector index off end of list:\n sftNum1=%d, sftNum2=%d, detInd1=%d, detInd2=%d, inputSFTs->length=%d\n",
		 sftNum1, sftNum2, detInd1, detInd2, inputSFTs->length );

    UINT8 sftInd1 = sftIndices->data[sftNum1].sftInd;
    UINT8 sftInd2 = sftIndices->data[sftNum2].sftInd;

    XLAL_CHECK ( ( sftInd1 < inputSFTs->data[detInd1]->length )
		 && ( sftInd2 < inputSFTs->data[detInd2]->length ),
		 XLAL_EINVAL,
		 "SFT asked for SFT index off end of list:\n sftNum1=%d, sftNum2=%d, detInd1=%d, detInd2=%d, sftInd1=%d, sftInd2=%d, inputSFTs->data[detInd1]->length=%d, inputSFTs->data[detInd2]->length=%d\n",
		 sftNum1, sftNum2, detInd1, detInd2, sftInd1, sftInd2,
		 inputSFTs->data[detInd1]->length,
		 inputSFTs->data[detInd2]->length );

    COMPLEX8 *dataArray1 = inputSFTs->data[detInd1]->data[sftInd1].data->data;
    COMPLEX8 *dataArray2 = inputSFTs->data[detInd2]->data[sftInd2].data->data;
    UINT4 lenDataArray1 = inputSFTs->data[detInd1]->data[sftInd1].data->length;
    UINT4 lenDataArray2 = inputSFTs->data[detInd2]->data[sftInd2].data->length;
    COMPLEX8 GalphaCC = curlyGAmp->data[alpha] * (expSignalPhases->data[sftNum1]/expSignalPhases->data[sftNum2]);
    INT4 baseCCSign = 1; /* Alternating sign is (-1)**(k1-k2) */
    if ( ( (lowestBins->data[sftNum1]-lowestBins->data[sftNum2]) % 2) != 0 ) {
      baseCCSign = -1;
    }

    UINT4 lowestBin1 = lowestBins->data[sftNum1];
    XLAL_CHECK ( ((lowestBin1 + numBins - 1) < lenDataArray1),
		 XLAL_EINVAL,
		 "Loop would run off end of array:\n lowestBin1=%d, numBins=%d, len(dataArray1)=%d\n",
		 lowestBin1, numBins, lenDataArray1 );
    for (UINT8 j=0; j < numBins; j++) {
      COMPLEX8 data1 = dataArray1[lowestBin1+j];
      /* Normalized sinc, i.e., sin(pi*x)/(pi*x) */
      INT4 ccSign = baseCCSign;
      UINT4 lowestBin2 = lowestBins->data[sftNum2];
      XLAL_CHECK ( ((lowestBin2 + numBins - 1) < lenDataArray2),
		   XLAL_EINVAL,
		   "Loop would run off end of array:\n lowestBin2=%d, numBins=%d, len(dataArray2)=%d\n",
		   lowestBin2, numBins, lenDataArray2 );
      for (UINT8 k=0; k < numBins; k++) {
	COMPLEX8 data2 = dataArray2[lowestBins->data[sftNum2]+k];
	REAL8 sincFactor =1;
	sincFactor = sincList->data[sftNum1 * numBins + k] * sincList->data[sftNum2 * numBins + k];
	nume +=  ccSign * sincFactor * creal ( GalphaCC * conj(data1) * data2 ); /*multiWeights->data[detInd1]->data[sftNum1] *  multiWeights->data[detInd2]->data[sftNum2] **/
	REAL8 GalphaAmp = curlyGAmp->data[alpha] * sincFactor ; /** multiWeights->data[detInd1]->data[sftNum1] *  multiWeights->data[detInd2]->data[sftNum2]*/
	curlyGSqr += SQUARE( GalphaAmp );
	ccSign *= -1;
      }
      baseCCSign *= -1;
    }
  }
  *evSquared = 8 * SQUARE(multiWeights->Sinv_Tsft) * curlyGSqr;
  *ccStat = 4 * multiWeights->Sinv_Tsft * nume / sqrt(*evSquared);
  return XLAL_SUCCESS;
}
/*calculate metric diagnol components, also include the estimation of sensitivity E[rho]/(h_0)^2*/
int XLALFindLMXBCrossCorrDiagMetric
  (
   REAL8             *hSens,         /*Output:sensitivity*/
   REAL8             *g_ff,          /*Output:metric elements */
   REAL8             *g_aa,          /*binary projected semimajor axis*/
   REAL8             *g_TT,          /*reference time*/
   PulsarDopplerParams DopplerParams, /* Input: pulsar/binary orbit paramaters*/
   REAL8Vector       *G_alpha,        /* Input: vector of curlyGunshifted values */
   SFTPairIndexList  *pairIndexList,  /* Input: list of SFT pairs */
   SFTIndexList      *indexList,      /* Input: list of SFTs */
   MultiSFTVector    *sfts            /* Input: set of per-detector SFT vectors */
   /* REAL8Vector       *kappaValues */    /* Input: Fractional offset of signal freq from best bin center */

   /*REAL8             *devTsq,  */   /*Output: mean time deviation^2*/
   /*REAL8             *g_pp,*/
   )
{
  UINT8 sftNum1=0;
  UINT8 sftNum2=0;
  UINT8 j=0;
  REAL8 T=0;
  REAL8 denom=0;
  REAL8 TSquaWeightedAve=0;
  REAL8 SinSquaWeightedAve=0;
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
    sftNum1 = pairIndexList->data[j].sftNum[0];
    sftNum2 = pairIndexList->data[j].sftNum[1];
    UINT8 detInd1 = indexList->data[sftNum1].detInd;
    UINT8 detInd2 = indexList->data[sftNum2].detInd;
    UINT8 sftInd1 = indexList->data[sftNum1].sftInd;
    UINT8 sftInd2 = indexList->data[sftNum2].sftInd;
    T1 = &(sfts->data[detInd1]->data[sftInd1].epoch);
    T2 = &(sfts->data[detInd2]->data[sftInd2].epoch);
    T = XLALGPSDiff( T1, T2 );
    REAL8 sqrG_alpha = SQUARE(G_alpha->data[j]);/*(curlyG_\alpha)^2*/
    sinSquare += sqrG_alpha*SQUARE(sin(LAL_PI*T/(DopplerParams.period)));/*(G_\alpha)^2*(sin(\pi*T/T_orbit))^2*/
    tSquare += sqrG_alpha*SQUARE(T); /*(\curlyg_alpha*)^2*T^2*/
    denom += sqrG_alpha;                               /*calculate the denominator*/
    rhosum += 2*sqrG_alpha;
    /*hfT=0.5*T;
      Tmean=XLALGPSAdd(&T2, hfT);*/
    /*muT +=Tmean/numalpha;*/                            /*calculate the average of Tmean*/
      }
  TSquaWeightedAve =(tSquare/denom);
  SinSquaWeightedAve =(sinSquare/denom);
  *hSens = sqrt(rhosum);
  *g_ff= TSquaWeightedAve* 2 * SQUARE(LAL_PI);
  *g_aa= SinSquaWeightedAve* SQUARE(LAL_PI*DopplerParams.fkdot[0]);
  *g_TT= SinSquaWeightedAve* SQUARE(2*SQUARE(LAL_PI)*(DopplerParams.fkdot[0])*(DopplerParams.asini)/(DopplerParams.period));


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
