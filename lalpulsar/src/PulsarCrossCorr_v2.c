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
#define QUAD(x) (x*x*x*x)
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
   REAL8                            Tsft  /**< SFT duration */
  )
{
  REAL8 timeDiff, factor, fhat, phiByTwoPi;

  UINT4 numSFTs = sftIndices->length;
  if ( expSignalPhases->length !=numSFTs
       || shiftedFreqs->length !=numSFTs
       || lowestBins->length !=numSFTs
       || sincList->length !=numSFTs) {
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
    timeDiff = times->DeltaT->data[sftInd]
      + XLALGPSDiff( &(times->refTime), &(dopp->refTime));
    fhat = dopp->fkdot[0]; /* initialization */
    phiByTwoPi = fmod ( fhat * timeDiff , 1.0 );
    factor = timeDiff;

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
      REAL4 sinPiX, cosPiX;
      REAL8 X;  /* Normalized sinc, i.e., sin(pi*x)/(pi*x) */
      X =  lowestBins->data[sftNum] - fminusf0 * Tsft + l;
      if(X > SINC_SAFETY || (X < - SINC_SAFETY)){
	     XLAL_CHECK( XLALSinCos2PiLUT( &sinPiX, &cosPiX, 0.5 * X ) == XLAL_SUCCESS, XLAL_EFUNC ); /*sin(2*pi*0.5*x)=sin(pi*x)*/
	     sincList->data[sftNum*numBins + l] = LAL_1_PI * sinPiX / X;/*1/(pi*x) =1/pi*1/x*/
	   }
	   else{
	     sincList->data[sftNum*numBins + l] = 1;
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
  *evSquared = 8 * SQUARE(multiWeights->Sinv_Tsft) * curlyGSqr;
  *ccStat = 4 * multiWeights->Sinv_Tsft * nume / sqrt(*evSquared);
  return XLAL_SUCCESS;
}

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


/*calculate metric diagonal components, also include the estimation of sensitivity E[rho]/(h_0)^2*/
int XLALCalculateLMXBCrossCorrDiagMetric
  (
   REAL8                      *hSens, /* Output: sensitivity*/
   REAL8                       *g_ff, /* Output: Diagonal frequency metric element */
   REAL8                       *g_aa, /* Output: Diagonal binary projected semimajor axis metric element*/
   REAL8                       *g_TT, /* Output: Diagonal reference time metric element*/
   REAL8                       *g_pp, /* Output: Diagonal orbital period metric element */
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

  return XLAL_SUCCESS;

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
