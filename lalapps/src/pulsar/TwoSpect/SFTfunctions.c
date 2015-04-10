/*
*  Copyright (C) 2014 Evan Goetz
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

#include <lal/UserInput.h>
#include <lal/CWMakeFakeData.h>
#include <lal/SinCosLUT.h>
#include "SFTfunctions.h"
#include "TwoSpectSpecFunc.h"
#include "statistics.h"

/**
 * Find the SFT data specified by user input
 * \param [in] params Pointer to the UserInput_t
 * \return Pointer to SFTCatalog
 */
SFTCatalog * findSFTdata(UserInput_t *params)
{

   XLAL_CHECK_NULL( params != NULL, XLAL_EINVAL );

   fprintf(LOG, "Finding SFTs... ");
   fprintf(stderr, "Finding SFTs... ");

   //Set the start and end times in the LIGO GPS format
   LIGOTimeGPS start = LIGOTIMEGPSZERO, end = LIGOTIMEGPSZERO;
   XLALGPSSetREAL8(&start, params->t0);
   XLAL_CHECK_NULL( xlalErrno == 0, XLAL_EFUNC );
   XLALGPSSetREAL8(&end, params->t0+params->Tobs-params->SFToverlap);
   XLAL_CHECK_NULL( xlalErrno == 0, XLAL_EFUNC );

   //Setup the constraints
   SFTConstraints XLAL_INIT_DECL(constraints);
   if (params->IFO->length == 1) constraints.detector = params->IFO->data[0];
   constraints.minStartTime = &start;
   constraints.maxStartTime = &end;

   //Find SFT files
   SFTCatalog *catalog = NULL;
   XLAL_CHECK_NULL( (catalog = XLALSFTdataFind(params->inputSFTs, &constraints)) != NULL, XLAL_EFUNC );

   fprintf(LOG, "done\n");
   fprintf(stderr, "done\n");

   return catalog;

}


/**
 * Extract the SFT coefficients from the band of interest
 * \param [in] catalog Pointer to the SFTCatalog
 * \param [in] minfbin Frequency value of the minimum frequency bin
 * \param [in] maxfbin Frequency value of the maximum frequency bin
 * \return Pointer to MultiSFTVector
 */
MultiSFTVector * extractSFTband(SFTCatalog *catalog, REAL8 minfbin, REAL8 maxfbin)
{

   XLAL_CHECK_NULL( catalog != NULL, XLAL_EINVAL );

   fprintf(LOG, "Extracting band from SFTs... ");
   fprintf(stderr, "Extracting band from SFTs... ");

   //Extract the data
   MultiSFTVector *sftvector = NULL;
   XLAL_CHECK_NULL( (sftvector = XLALLoadMultiSFTs(catalog, minfbin, maxfbin)) != NULL, XLAL_EFUNC );

   fprintf(LOG, "done\n");
   fprintf(stderr, "done\n");

   return sftvector;

}


/**
 * Get a MultiSFTVector from the user-input values
 * \param [in] params Pointer to UserInput_t
 * \param [in] minfbin Frequency value of the minimum frequency bin
 * \param [in] maxfbin Frequency value of the maximum frequency bin
 * \return Pointer to MultiSFTVector
 */
MultiSFTVector * getMultiSFTVector(UserInput_t *params, REAL8 minfbin, REAL8 maxfbin)
{
   //Get catalog of SFTs
   SFTCatalog *catalog = NULL;
   XLAL_CHECK_NULL( (catalog = findSFTdata(params)) != NULL, XLAL_EFUNC );

   MultiSFTVector *sftvector = NULL;
   if (params->markBadSFTs && !params->signalOnly) {
      //Extract band to get a MultiSFTVector
      MultiSFTVector *tmpsftvector = NULL;
      XLAL_CHECK_NULL( (tmpsftvector = extractSFTband(catalog, minfbin, maxfbin)) != NULL, XLAL_EFUNC );

      //Get the timestamps of the SFTs applying the KS/Kuipers tests if desired
      MultiLIGOTimeGPSVector *multiTimestamps = NULL;
      XLAL_CHECK_NULL( (multiTimestamps = getMultiTimeStampsFromSFTs(tmpsftvector, params)) != NULL, XLAL_EFUNC );

      //Get the SFT subset
      XLAL_CHECK_NULL( (sftvector = XLALExtractMultiSFTVectorWithMultiTimestamps(tmpsftvector, multiTimestamps)) != NULL, XLAL_EFUNC );

      XLALDestroyMultiSFTVector(tmpsftvector);
      XLALDestroyMultiTimestamps(multiTimestamps);
   } else {
      XLAL_CHECK_NULL( (sftvector = extractSFTband(catalog, minfbin, maxfbin)) != NULL, XLAL_EFUNC );
   }

   XLALDestroySFTCatalog(catalog);

   return sftvector;
}


/**
 * Add SFTs together from a MultiSFTVector
 * \param [in] multiSFTvector      Pointer to a MultiSFTVector containing the SFT data
 * \param [in] multissb            Pointer to a MultiSSBtimes structure
 * \param [in] multiAMcoefficients Pointer to a MultiAMCoeffs structure
 * \param [in] cosiSign            Value of 1, 0, -1 where 1 means average cosi over [0,1], 0 means average over [-1,1], and -1 means average over [-1,0]
 * \param [in] assumeNScosi        Pointer to the assumed cosi value, or NULL if no value assumed
 * \param [in] assumeNSpsi         Pointer to the assumed psi value, or NULL if no value assumed
 * \param [in] params              Pointer to UserInput_t
 * \return REAL4Vector of powers of the coherently combined SFTs
 */
REAL4Vector * coherentlyAddSFTs(MultiSFTVector *multiSFTvector, MultiSSBtimes *multissb, MultiAMCoeffs *multiAMcoefficients, INT4 cosiSign, REAL8 *assumeNScosi, REAL8 *assumeNSpsi, UserInput_t *params) {

   //Sanity check the values
   XLAL_CHECK_NULL( multiSFTvector != NULL, XLAL_EINVAL );
   XLAL_CHECK_NULL( multissb != NULL, XLAL_EINVAL );
   XLAL_CHECK_NULL( multiAMcoefficients != NULL, XLAL_EINVAL );
   XLAL_CHECK_NULL( params != NULL, XLAL_EINVAL );

   fprintf(stderr, "Coherently adding SFT data... ");

   //Determine length of the largest SFTvector and create an SFTvector large enough to hold this
   UINT4 maxSFTs = multiSFTvector->data[0]->length;
   for (UINT4 ii=1; ii<multiSFTvector->length; ii++) {
      if ( multiSFTvector->data[ii]->length > maxSFTs ) {
         maxSFTs = multiSFTvector->data[ii]->length;
      }
   }
   SFTVector *combinedSFTs = NULL;
   XLAL_CHECK_NULL( (combinedSFTs = XLALCreateSFTVector(maxSFTs, 0)) != NULL, XLAL_EFUNC );

   //Create an INT4Vector to determine at which vector we are using in the multiSFTvector
   INT4Vector *whichSFTinMultiSFTvector = NULL, *whichIFOsToBeUsed = NULL;
   XLAL_CHECK_NULL( (whichSFTinMultiSFTvector = XLALCreateINT4Vector(multiSFTvector->length)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (whichIFOsToBeUsed = XLALCreateINT4Vector(multiSFTvector->length)) != NULL, XLAL_EFUNC );
   memset(whichSFTinMultiSFTvector->data, 0, whichSFTinMultiSFTvector->length * sizeof(INT4));  //All index values are set to zero at the beginning

   //Pre-compute the different sin(2*psi) and cos(2*psi) values
   REAL4VectorAligned *twoPsiVector = NULL, *sin2psiVector = NULL, *cos2psiVector = NULL, *Fplus0s = NULL, *Fcross0s = NULL, *FplusXs = NULL, *FcrossXs = NULL, *aVals = NULL, *bVals = NULL;
   XLAL_CHECK_NULL( (twoPsiVector = XLALCreateREAL4VectorAligned(16, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (sin2psiVector = XLALCreateREAL4VectorAligned(twoPsiVector->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (cos2psiVector = XLALCreateREAL4VectorAligned(twoPsiVector->length, 32)) != NULL, XLAL_EFUNC );
   for (UINT4 ii=0; ii<twoPsiVector->length; ii++) twoPsiVector->data[ii] = 0.0625*ii;
   XLAL_CHECK_NULL( XLALVectorSinCos2PiREAL4(sin2psiVector->data, cos2psiVector->data, twoPsiVector->data, twoPsiVector->length) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK_NULL( (Fplus0s = XLALCreateREAL4VectorAligned(twoPsiVector->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (Fcross0s = XLALCreateREAL4VectorAligned(twoPsiVector->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (FplusXs = XLALCreateREAL4VectorAligned(twoPsiVector->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (FcrossXs = XLALCreateREAL4VectorAligned(twoPsiVector->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (aVals = XLALCreateREAL4VectorAligned(twoPsiVector->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (bVals = XLALCreateREAL4VectorAligned(twoPsiVector->length, 32)) != NULL, XLAL_EFUNC );

   //Pre-compute the frequencies and the bin number of the SFT bins
   alignedREAL8Vector *frequencies = NULL, *frequencyBins = NULL;
   XLAL_CHECK_NULL( (frequencies = createAlignedREAL8Vector(multiSFTvector->data[0]->data[0].data->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (frequencyBins = createAlignedREAL8Vector(multiSFTvector->data[0]->data[0].data->length, 32)) != NULL, XLAL_EFUNC );
   for (UINT4 ii=0; ii<frequencyBins->length; ii++) frequencies->data[ii] = multiSFTvector->data[0]->data[0].f0 + multiSFTvector->data[0]->data[0].deltaF*ii;
   if (params->vectorMath==1) XLAL_CHECK_NULL( sseScaleREAL8Vector(frequencyBins, frequencies, params->Tsft) == XLAL_SUCCESS, XLAL_EFUNC );
   else if (params->vectorMath==2) XLAL_CHECK_NULL( avxScaleREAL8Vector(frequencyBins, frequencies, params->Tsft) == XLAL_SUCCESS, XLAL_EFUNC );
   else for (UINT4 ii=0; ii<frequencyBins->length; ii++) frequencyBins->data[ii] = frequencies->data[ii]*params->Tsft;

   //Pre-allocate delta values for IFO_0 and IFO_X, also 2*pi*f_k*tau
   alignedREAL8Vector *delta0vals = NULL, *deltaXvals = NULL, *TwoPiFrequenciesTau = NULL;
   XLAL_CHECK_NULL( (delta0vals = createAlignedREAL8Vector(frequencyBins->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (deltaXvals = createAlignedREAL8Vector(frequencyBins->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (TwoPiFrequenciesTau = createAlignedREAL8Vector(frequencyBins->length, 32)) != NULL, XLAL_EFUNC );

   //Mid-point of the SFT
   REAL8 Tmid = 0.5*params->Tsft;

   //Pre-compute 2*pi*{[DeltaTs(jj) - Tmid*Tdots(jj)]-[DeltaTs(0) - Tmid*Tdots(0)]} and [Tdots(jj) - 1] for jj=0,...,N-1
   //Note that twoPiTauVals->data[0] is unused later in the code; it is only necessary for the intermediate calculation
   alignedREAL8VectorArray *twoPiTauVals = NULL, *TdotMinus1s = NULL;
   alignedREAL8Vector *Tdots = NULL, *DeltaTs = NULL;
   XLAL_CHECK_NULL( (twoPiTauVals = createAlignedREAL8VectorArray(multissb->length, multissb->data[0]->DeltaT->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (TdotMinus1s = createAlignedREAL8VectorArray(multissb->length, multissb->data[0]->Tdot->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (Tdots = createAlignedREAL8Vector(multissb->data[0]->Tdot->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (DeltaTs = createAlignedREAL8Vector(multissb->data[0]->DeltaT->length, 32)) != NULL, XLAL_EFUNC );
   for (UINT4 ii=0; ii<twoPiTauVals->length; ii++) {
      if (params->vectorMath==1) {
         memcpy(Tdots->data, multissb->data[ii]->Tdot->data, sizeof(REAL8)*multissb->data[ii]->Tdot->length);
         memcpy(DeltaTs->data, multissb->data[ii]->DeltaT->data, sizeof(REAL8)*multissb->data[ii]->DeltaT->length);
         XLAL_CHECK_NULL( sseScaleREAL8Vector(twoPiTauVals->data[ii], Tdots, -Tmid) == XLAL_SUCCESS, XLAL_EFUNC );
         XLAL_CHECK_NULL( sseDDVectorSum(twoPiTauVals->data[ii], twoPiTauVals->data[ii], DeltaTs) == XLAL_SUCCESS, XLAL_EFUNC );
         XLAL_CHECK_NULL( sseAddScalarToREAL8Vector(TdotMinus1s->data[ii], Tdots, -1.0) == XLAL_SUCCESS, XLAL_EFUNC );
      } else if (params->vectorMath==2) {
         memcpy(Tdots->data, multissb->data[ii]->Tdot->data, sizeof(REAL8)*multissb->data[ii]->Tdot->length);
         memcpy(DeltaTs->data, multissb->data[ii]->DeltaT->data, sizeof(REAL8)*multissb->data[ii]->DeltaT->length);
         XLAL_CHECK_NULL( avxScaleREAL8Vector(twoPiTauVals->data[ii], Tdots, -Tmid) == XLAL_SUCCESS, XLAL_EFUNC );
         XLAL_CHECK_NULL( avxDDVectorSum(twoPiTauVals->data[ii], twoPiTauVals->data[ii], DeltaTs) == XLAL_SUCCESS, XLAL_EFUNC );
         XLAL_CHECK_NULL( avxAddScalarToREAL8Vector(TdotMinus1s->data[ii], Tdots, -1.0) == XLAL_SUCCESS, XLAL_EFUNC );
      } else {
         for (UINT4 jj=0; jj<Tdots->length; jj++) {
            twoPiTauVals->data[ii]->data[jj] = multissb->data[ii]->DeltaT->data[jj] - Tmid*multissb->data[ii]->Tdot->data[jj];
            TdotMinus1s->data[ii]->data[jj] = multissb->data[ii]->Tdot->data[jj] - 1.0;
         }
      }
   }
   for (UINT4 ii=twoPiTauVals->length; ii>0; ii--) {
      if (params->vectorMath==1) XLAL_CHECK_NULL( sseDDVectorSubtract(twoPiTauVals->data[ii-1], twoPiTauVals->data[ii-1], twoPiTauVals->data[0]) == XLAL_SUCCESS, XLAL_EFUNC );
      else if (params->vectorMath==2) XLAL_CHECK_NULL( avxDDVectorSubtract(twoPiTauVals->data[ii-1], twoPiTauVals->data[ii-1], twoPiTauVals->data[0]) == XLAL_SUCCESS, XLAL_EFUNC );
      else for (UINT4 jj=0; jj<Tdots->length; jj++) twoPiTauVals->data[ii-1]->data[jj] -= twoPiTauVals->data[0]->data[jj];
   }
   for (UINT4 ii=1; ii<twoPiTauVals->length; ii--) {
      if (params->vectorMath==1) XLAL_CHECK_NULL( sseScaleREAL8Vector(twoPiTauVals->data[ii], twoPiTauVals->data[ii], LAL_TWOPI) == XLAL_SUCCESS, XLAL_EFUNC );
      else if (params->vectorMath==2) XLAL_CHECK_NULL( avxScaleREAL8Vector(twoPiTauVals->data[ii], twoPiTauVals->data[ii], LAL_TWOPI) == XLAL_SUCCESS, XLAL_EFUNC );
      else for (UINT4 jj=0; jj<Tdots->length; jj++) twoPiTauVals->data[ii]->data[jj] *= LAL_TWOPI;
   }
   destroyAlignedREAL8Vector(Tdots);
   destroyAlignedREAL8Vector(DeltaTs);

   //Precompute cosi and (cosi*cosi+1)/2 values
   REAL4VectorAligned *cosiVals = NULL, *onePlusCosiSqOver2Vals = NULL;
   XLAL_CHECK_NULL( (cosiVals = XLALCreateREAL4VectorAligned(21, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (onePlusCosiSqOver2Vals = XLALCreateREAL4VectorAligned(21, 32)) != NULL, XLAL_EFUNC );
   if (cosiSign==1) for (UINT4 ii=0; ii<cosiVals->length; ii++) cosiVals->data[ii] = 1.0 - 0.05*ii;
   else if (cosiSign==-1) for (UINT4 ii=0; ii<cosiVals->length; ii++) cosiVals->data[ii] = -0.05*ii;
   else for (UINT4 ii=0; ii<cosiVals->length; ii++) cosiVals->data[ii] = 1.0 - 0.05*ii*2.0;
   if (params->vectorMath==1) {
      XLAL_CHECK_NULL( sseSSVectorMultiply(onePlusCosiSqOver2Vals, cosiVals, cosiVals) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_NULL( sseAddScalarToREAL4Vector(onePlusCosiSqOver2Vals, onePlusCosiSqOver2Vals, -1.0) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_NULL( sseScaleREAL4Vector(onePlusCosiSqOver2Vals, onePlusCosiSqOver2Vals, 0.5) == XLAL_SUCCESS, XLAL_EFUNC );
   } else if (params->vectorMath==2) {
      XLAL_CHECK_NULL( avxSSVectorMultiply(onePlusCosiSqOver2Vals, cosiVals, cosiVals) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_NULL( avxAddScalarToREAL4Vector(onePlusCosiSqOver2Vals, onePlusCosiSqOver2Vals, 1.0) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_NULL( avxScaleREAL4Vector(onePlusCosiSqOver2Vals, onePlusCosiSqOver2Vals, 0.5) == XLAL_SUCCESS, XLAL_EFUNC );
   } else {
      for (UINT4 ii=0; ii<cosiVals->length; ii++) onePlusCosiSqOver2Vals->data[ii] = 0.5*(cosiVals->data[ii]*cosiVals->data[ii] + 1.0);
   }

   //Pre-allocate Aplus and Across
   REAL4VectorAligned *Aplus0s = NULL, *AplusXs = NULL, *Across0s = NULL, *AcrossXs = NULL;
   XLAL_CHECK_NULL( (Aplus0s = XLALCreateREAL4VectorAligned(twoPsiVector->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (Across0s = XLALCreateREAL4VectorAligned(twoPsiVector->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (AplusXs = XLALCreateREAL4VectorAligned(twoPsiVector->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (AcrossXs = XLALCreateREAL4VectorAligned(twoPsiVector->length, 32)) != NULL, XLAL_EFUNC );

   //Pre-allocate DirichletScaling0 and DirichletScaling1
   alignedREAL8Vector *DirichletScaling0 = NULL, *DirichletScalingX = NULL;
   XLAL_CHECK_NULL( (DirichletScaling0 = createAlignedREAL8Vector(delta0vals->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (DirichletScalingX = createAlignedREAL8Vector(delta0vals->length, 32)) != NULL, XLAL_EFUNC );

   //Loop over the combinedSFTs vector to fill it with single or coherently combined SFTs
   for (UINT4 ii=0; ii<combinedSFTs->length; ii++) {
      //Loop over the interferometers, determining which ones to use
      memset(whichIFOsToBeUsed->data, 0, whichIFOsToBeUsed->length * sizeof(INT4));  //All index values are set to zero each time
      LIGOTimeGPS smallestGPS = multiSFTvector->data[0]->data[whichSFTinMultiSFTvector->data[0]].epoch;
      INT4 ifoWithSmallestGPS = 0;
      for (UINT4 jj=1; jj<multiSFTvector->length; jj++) {
         INT4 compareT = XLALGPSCmp(&(multiSFTvector->data[ifoWithSmallestGPS]->data[whichSFTinMultiSFTvector->data[ifoWithSmallestGPS]].epoch), &(multiSFTvector->data[jj]->data[whichSFTinMultiSFTvector->data[jj]].epoch));
         if (compareT>0) {
            smallestGPS = multiSFTvector->data[jj]->data[whichSFTinMultiSFTvector->data[jj]].epoch;
            ifoWithSmallestGPS = jj;
         }
      }
      for (UINT4 jj=0; jj<multiSFTvector->length; jj++) {
         INT4 compareT = XLALGPSCmp(&(smallestGPS), &(multiSFTvector->data[jj]->data[whichSFTinMultiSFTvector->data[jj]].epoch));
         if (compareT==0) whichIFOsToBeUsed->data[jj] = 1;
      }

      BOOLEAN createSFT = 1, computeAntenna0 = 1, computeDelta0vals = 1;
      for (UINT4 jj=0; jj<multiSFTvector->length; jj++) {
         if (jj==0 && whichIFOsToBeUsed->data[jj]==1) {
            //Copy the data from the multiSFTvector into the combinedSFTs vector
            XLAL_CHECK_NULL( XLALCopySFT(&(combinedSFTs->data[ii]), &(multiSFTvector->data[jj]->data[whichSFTinMultiSFTvector->data[jj]])) == XLAL_SUCCESS, XLAL_EFUNC );
            createSFT = 0;
            whichSFTinMultiSFTvector->data[jj]++;
         } else if (jj>0 && whichIFOsToBeUsed->data[jj]==1) {
            //Create a copy of the SFT to be shifted since we will manipulate the SFT coefficients
            SFTtype *sftcopy = NULL;
            XLAL_CHECK_NULL( (sftcopy = XLALCreateSFT(0)) != NULL, XLAL_EFUNC );
            XLAL_CHECK_NULL( XLALCopySFT(sftcopy, &(multiSFTvector->data[jj]->data[whichSFTinMultiSFTvector->data[jj]])) == XLAL_SUCCESS, XLAL_EFUNC );
            REAL8 sftstart = XLALGPSGetREAL8(&(sftcopy->epoch));
            XLAL_CHECK_NULL( xlalErrno == 0, XLAL_EFUNC );
            INT4 fftnum = (INT4)round((sftstart - params->t0)/params->SFToverlap);

            //First do time of arrival: 2*pi*fk*tau
            if (params->vectorMath==1) XLAL_CHECK_NULL( sseScaleREAL8Vector(TwoPiFrequenciesTau, frequencies, twoPiTauVals->data[jj]->data[fftnum]) == XLAL_SUCCESS, XLAL_EFUNC );
            else if (params->vectorMath==2) XLAL_CHECK_NULL( avxScaleREAL8Vector(TwoPiFrequenciesTau, frequencies, twoPiTauVals->data[jj]->data[fftnum]) == XLAL_SUCCESS, XLAL_EFUNC );
            else for (UINT4 kk=0; kk<frequencies->length; kk++) TwoPiFrequenciesTau->data[kk] = frequencies->data[kk]*twoPiTauVals->data[jj]->data[fftnum];

            //If assumNSpsi is not given compute Fplus and Fcross for detector ratio
            if (assumeNSpsi==NULL) {
               if (params->vectorMath==1) {
                  //IFO 0
                  if (computeAntenna0) {
                     XLAL_CHECK_NULL( sseScaleREAL4Vector(aVals, cos2psiVector, multiAMcoefficients->data[0]->a->data[fftnum]) == XLAL_SUCCESS, XLAL_EFUNC );
                     XLAL_CHECK_NULL( sseScaleREAL4Vector(bVals, sin2psiVector, multiAMcoefficients->data[0]->b->data[fftnum]) == XLAL_SUCCESS, XLAL_EFUNC );
                     XLAL_CHECK_NULL( sseSSVectorSum(Fplus0s, aVals, bVals) == XLAL_SUCCESS, XLAL_EFUNC );
                     XLAL_CHECK_NULL( sseScaleREAL4Vector(bVals, cos2psiVector, multiAMcoefficients->data[0]->b->data[fftnum]) == XLAL_SUCCESS, XLAL_EFUNC );
                     XLAL_CHECK_NULL( sseScaleREAL4Vector(aVals, sin2psiVector, -multiAMcoefficients->data[0]->a->data[fftnum]) == XLAL_SUCCESS, XLAL_EFUNC );
                     XLAL_CHECK_NULL( sseSSVectorSum(Fcross0s, bVals, aVals) == XLAL_SUCCESS, XLAL_EFUNC );
                     computeAntenna0 = 0;
                  }
                  //IFO jj
                  XLAL_CHECK_NULL( sseScaleREAL4Vector(aVals, cos2psiVector, multiAMcoefficients->data[jj]->a->data[fftnum]) == XLAL_SUCCESS, XLAL_EFUNC );
                  XLAL_CHECK_NULL( sseScaleREAL4Vector(bVals, sin2psiVector, multiAMcoefficients->data[jj]->b->data[fftnum]) == XLAL_SUCCESS, XLAL_EFUNC );
                  XLAL_CHECK_NULL( sseSSVectorSum(FplusXs, aVals, bVals) == XLAL_SUCCESS, XLAL_EFUNC );
                  XLAL_CHECK_NULL( sseScaleREAL4Vector(bVals, cos2psiVector, multiAMcoefficients->data[jj]->b->data[fftnum]) == XLAL_SUCCESS, XLAL_EFUNC );
                  XLAL_CHECK_NULL( sseScaleREAL4Vector(aVals, sin2psiVector, -multiAMcoefficients->data[jj]->a->data[fftnum]) == XLAL_SUCCESS, XLAL_EFUNC );
                  XLAL_CHECK_NULL( sseSSVectorSum(FcrossXs, bVals, aVals) == XLAL_SUCCESS, XLAL_EFUNC );
               } else if (params->vectorMath==2) {
                  //IFO 0
                  if (computeAntenna0) {
                     XLAL_CHECK_NULL( avxScaleREAL4Vector(aVals, cos2psiVector, multiAMcoefficients->data[0]->a->data[fftnum]) == XLAL_SUCCESS, XLAL_EFUNC );
                     XLAL_CHECK_NULL( avxScaleREAL4Vector(bVals, sin2psiVector, multiAMcoefficients->data[0]->b->data[fftnum]) == XLAL_SUCCESS, XLAL_EFUNC );
                     XLAL_CHECK_NULL( avxSSVectorSum(Fplus0s, aVals, bVals) == XLAL_SUCCESS, XLAL_EFUNC );
                     XLAL_CHECK_NULL( avxScaleREAL4Vector(bVals, cos2psiVector, multiAMcoefficients->data[0]->b->data[fftnum]) == XLAL_SUCCESS, XLAL_EFUNC );
                     XLAL_CHECK_NULL( avxScaleREAL4Vector(aVals, sin2psiVector, -multiAMcoefficients->data[0]->a->data[fftnum]) == XLAL_SUCCESS, XLAL_EFUNC );
                     XLAL_CHECK_NULL( avxSSVectorSum(Fcross0s, bVals, aVals) == XLAL_SUCCESS, XLAL_EFUNC );
                     computeAntenna0 = 0;
                  }
                  //IFO jj
                  XLAL_CHECK_NULL( avxScaleREAL4Vector(aVals, cos2psiVector, multiAMcoefficients->data[jj]->a->data[fftnum]) == XLAL_SUCCESS, XLAL_EFUNC );
                  XLAL_CHECK_NULL( avxScaleREAL4Vector(bVals, sin2psiVector, multiAMcoefficients->data[jj]->b->data[fftnum]) == XLAL_SUCCESS, XLAL_EFUNC );
                  XLAL_CHECK_NULL( avxSSVectorSum(FplusXs, aVals, bVals) == XLAL_SUCCESS, XLAL_EFUNC );
                  XLAL_CHECK_NULL( avxScaleREAL4Vector(bVals, cos2psiVector, multiAMcoefficients->data[jj]->b->data[fftnum]) == XLAL_SUCCESS, XLAL_EFUNC );
                  XLAL_CHECK_NULL( avxScaleREAL4Vector(aVals, sin2psiVector, -multiAMcoefficients->data[jj]->a->data[fftnum]) == XLAL_SUCCESS, XLAL_EFUNC );
                  XLAL_CHECK_NULL( avxSSVectorSum(FcrossXs, bVals, aVals) == XLAL_SUCCESS, XLAL_EFUNC );
               } else {
                  for (UINT4 kk=0; kk<Fplus0s->length; kk++) {
                     if (computeAntenna0) {
                        Fplus0s->data[kk] = multiAMcoefficients->data[0]->a->data[fftnum]*cos2psiVector->data[kk] + multiAMcoefficients->data[0]->b->data[fftnum]*sin2psiVector->data[kk];
                        Fcross0s->data[kk] = multiAMcoefficients->data[0]->b->data[fftnum]*cos2psiVector->data[kk] - multiAMcoefficients->data[0]->a->data[fftnum]*sin2psiVector->data[kk];
                     }
                     FplusXs->data[kk] = multiAMcoefficients->data[jj]->a->data[fftnum]*cos2psiVector->data[kk] + multiAMcoefficients->data[jj]->b->data[fftnum]*sin2psiVector->data[kk];
                     FcrossXs->data[kk] = multiAMcoefficients->data[jj]->b->data[fftnum]*cos2psiVector->data[kk] - multiAMcoefficients->data[jj]->a->data[fftnum]*sin2psiVector->data[kk];
                  }
                  if (computeAntenna0) computeAntenna0 = 0;
               }
            }

            //Average of detector-signal phase relation, unless assumed NS orientation and GW polarization angle
            //Need to compute (Aplus^X+i*Across^X)/(Aplus^0+i*Across^0) to determine an average value of magnitude and phase
            //where exponent indicates different IFO <X>.
            REAL4 detPhaseArg = 0.0, detPhaseMag = 0.0;
            if (assumeNScosi==NULL && assumeNSpsi==NULL) {
               BOOLEAN loopbroken = 0;
               for (UINT4 kk=0; kk<cosiVals->length && !loopbroken; kk++) {
                  if (params->vectorMath==1) {
                     XLAL_CHECK_NULL( sseScaleREAL4Vector(AplusXs, FplusXs, onePlusCosiSqOver2Vals->data[kk]) == XLAL_SUCCESS, XLAL_EFUNC );
                     XLAL_CHECK_NULL( sseScaleREAL4Vector(AcrossXs, FcrossXs, cosiVals->data[kk]) == XLAL_SUCCESS, XLAL_EFUNC );
                     XLAL_CHECK_NULL( sseScaleREAL4Vector(Aplus0s, Fplus0s, onePlusCosiSqOver2Vals->data[kk]) == XLAL_SUCCESS, XLAL_EFUNC );
                     XLAL_CHECK_NULL( sseScaleREAL4Vector(Across0s, Fcross0s, cosiVals->data[kk]) == XLAL_SUCCESS, XLAL_EFUNC );
                  } else if (params->vectorMath==2) {
                     XLAL_CHECK_NULL( avxScaleREAL4Vector(AplusXs, FplusXs, onePlusCosiSqOver2Vals->data[kk]) == XLAL_SUCCESS, XLAL_EFUNC );
                     XLAL_CHECK_NULL( avxScaleREAL4Vector(AcrossXs, FcrossXs, cosiVals->data[kk]) == XLAL_SUCCESS, XLAL_EFUNC );
                     XLAL_CHECK_NULL( avxScaleREAL4Vector(Aplus0s, Fplus0s, onePlusCosiSqOver2Vals->data[kk]) == XLAL_SUCCESS, XLAL_EFUNC );
                     XLAL_CHECK_NULL( avxScaleREAL4Vector(Across0s, Fcross0s, cosiVals->data[kk]) == XLAL_SUCCESS, XLAL_EFUNC );
                  } else {
                     for (UINT4 ll=0; ll<cos2psiVector->length; ll++) {
                        AplusXs->data[ll] = onePlusCosiSqOver2Vals->data[kk]*FplusXs->data[ll];
                        AcrossXs->data[ll] = cosiVals->data[kk]*FcrossXs->data[ll];
                        Aplus0s->data[ll] = onePlusCosiSqOver2Vals->data[kk]*Fplus0s->data[ll];
                        Across0s->data[ll] = cosiVals->data[kk]*Fcross0s->data[ll];
                     }
                  }
                  for (UINT4 ll=0; ll<cos2psiVector->length && !loopbroken; ll++) {
                     COMPLEX8 complexnumerator = crectf(AplusXs->data[ll], AcrossXs->data[ll]);
                     COMPLEX8 complexdenominator = crectf(Aplus0s->data[ll], Across0s->data[ll]);
                     if (cabsf(complexdenominator)>1.0e-6) {
                        COMPLEX8 complexval = complexnumerator/complexdenominator;
                        detPhaseMag += fminf(cabsf(complexval), 10.0);  //fmin here because sometimes the magnitude value is being divided by a small value and causing biases
                        detPhaseArg += (REAL4)gsl_sf_angle_restrict_pos((REAL8)cargf(complexval));
                     } else {
                        loopbroken = 1;
                        detPhaseMag = 0.0;
                        detPhaseArg = 0.0;
                     }
                  }
               }
               detPhaseMag /= (REAL4)(cosiVals->length*cos2psiVector->length);
               detPhaseArg /= (REAL4)(cosiVals->length*cos2psiVector->length);
            } else if (assumeNScosi!=NULL && assumeNSpsi==NULL) {
               BOOLEAN loopbroken = 0;
               REAL4 cosi = *assumeNScosi;
               REAL4 onePlusCosiSqOver2 = 0.5*(1.0 + cosi*cosi);
               if (params->vectorMath==1) {
                  XLAL_CHECK_NULL( sseScaleREAL4Vector(AplusXs, FplusXs, onePlusCosiSqOver2) == XLAL_SUCCESS, XLAL_EFUNC );
                  XLAL_CHECK_NULL( sseScaleREAL4Vector(AcrossXs, FcrossXs, cosi) == XLAL_SUCCESS, XLAL_EFUNC );
                  XLAL_CHECK_NULL( sseScaleREAL4Vector(Aplus0s, Fplus0s, onePlusCosiSqOver2) == XLAL_SUCCESS, XLAL_EFUNC );
                  XLAL_CHECK_NULL( sseScaleREAL4Vector(Across0s, Fcross0s, cosi) == XLAL_SUCCESS, XLAL_EFUNC );
               } else if (params->vectorMath==2) {
                  XLAL_CHECK_NULL( avxScaleREAL4Vector(AplusXs, FplusXs, onePlusCosiSqOver2) == XLAL_SUCCESS, XLAL_EFUNC );
                  XLAL_CHECK_NULL( avxScaleREAL4Vector(AcrossXs, FcrossXs, cosi) == XLAL_SUCCESS, XLAL_EFUNC );
                  XLAL_CHECK_NULL( avxScaleREAL4Vector(Aplus0s, Fplus0s, onePlusCosiSqOver2) == XLAL_SUCCESS, XLAL_EFUNC );
                  XLAL_CHECK_NULL( avxScaleREAL4Vector(Across0s, Fcross0s, cosi) == XLAL_SUCCESS, XLAL_EFUNC );
               } else {
                  for (UINT4 kk=0; kk<cos2psiVector->length; kk++) {
                     AplusXs->data[kk] = onePlusCosiSqOver2*FplusXs->data[kk];
                     AcrossXs->data[kk] = cosi*FcrossXs->data[kk];
                     Aplus0s->data[kk] = onePlusCosiSqOver2*Fplus0s->data[kk];
                     Across0s->data[kk] = cosi*Fcross0s->data[kk];
                  }
               }
               for (UINT4 kk=0; kk<cos2psiVector->length && !loopbroken; kk++) {
                  COMPLEX16 complexnumerator = crect(AplusXs->data[kk], AcrossXs->data[kk]);
                  COMPLEX16 complexdenominator = crect(Aplus0s->data[kk], Across0s->data[kk]);
                  if (cabs(complexdenominator)>1.0e-7) {
                     COMPLEX16 complexval = complexnumerator/complexdenominator;
                     detPhaseMag += fmin(cabs(complexval), 10.0);  //fmin here because sometimes the magnitude value is being divided by a small value and causing biases
                     detPhaseArg += gsl_sf_angle_restrict_pos(carg(complexval));
                  } else {
                     loopbroken = 1;
                     detPhaseMag = 0.0;
                     detPhaseArg = 0.0;
                  }
               }
               detPhaseMag /= (REAL4)(cos2psiVector->length);
               detPhaseArg /= (REAL4)(cos2psiVector->length);
            } else if (assumeNScosi==NULL && assumeNSpsi!=NULL) {
               REAL4 psi = *assumeNSpsi;
               REAL4 sin2psi = 0.0, cos2psi = 0.0;
               XLAL_CHECK_NULL( XLALSinCosLUT(&sin2psi, &cos2psi, 2.0*psi) == XLAL_SUCCESS, XLAL_EFUNC );
               if (sin2psi>1.0) sin2psi = 1.0;
               if (sin2psi<-1.0) sin2psi = -1.0;
               if (cos2psi>1.0) cos2psi = 1.0;
               if (cos2psi<-1.0) cos2psi = -1.0;
               REAL4 Fplus0 = multiAMcoefficients->data[0]->a->data[fftnum]*cos2psi + multiAMcoefficients->data[0]->b->data[fftnum]*sin2psi;
               REAL4 Fcross0 = multiAMcoefficients->data[0]->b->data[fftnum]*cos2psi - multiAMcoefficients->data[0]->a->data[fftnum]*sin2psi;
               REAL4 FplusX = multiAMcoefficients->data[jj]->a->data[fftnum]*cos2psi + multiAMcoefficients->data[jj]->b->data[fftnum]*sin2psi;
               REAL4 FcrossX = multiAMcoefficients->data[jj]->b->data[fftnum]*cos2psi - multiAMcoefficients->data[jj]->a->data[fftnum]*sin2psi;
               BOOLEAN loopbroken = 0;
               for (UINT4 kk=0; kk<cosiVals->length && !loopbroken; kk++) {
                  COMPLEX16 complexnumerator = crect(FplusX*onePlusCosiSqOver2Vals->data[kk], FcrossX*cosiVals->data[kk]);
                  COMPLEX16 complexdenominator = crect(Fplus0*onePlusCosiSqOver2Vals->data[kk], Fcross0*cosiVals->data[kk]);
                  if (cabs(complexdenominator)>1.0e-7) {
                     COMPLEX16 complexval = complexnumerator/complexdenominator;
                     detPhaseMag += fmin(cabs(complexval), 10.0);  //fmin here because sometimes the magnitude value is being divided by a small value and causing biases
                     detPhaseArg += gsl_sf_angle_restrict_pos(carg(complexval));
                  } else {
                     loopbroken = 1;
                     detPhaseMag = 0.0;
                     detPhaseArg = 0.0;
                  }
               }
               detPhaseMag /= (REAL4)(cosiVals->length);
               detPhaseArg /= (REAL4)(cosiVals->length);
            } else {
               REAL4 psi = *assumeNSpsi;
               REAL4 sin2psi = 0.0, cos2psi = 0.0;
               XLAL_CHECK_NULL( XLALSinCosLUT(&sin2psi, &cos2psi, 2.0*psi) == XLAL_SUCCESS, XLAL_EFUNC );
               if (sin2psi>1.0) sin2psi = 1.0;
               else if (sin2psi<-1.0) sin2psi = -1.0;
               if (cos2psi>1.0) cos2psi = 1.0;
               else if (cos2psi<-1.0) cos2psi = -1.0;
               REAL4 Fplus0 = multiAMcoefficients->data[0]->a->data[fftnum]*cos2psi + multiAMcoefficients->data[0]->b->data[fftnum]*sin2psi;
               REAL4 Fcross0 = multiAMcoefficients->data[0]->b->data[fftnum]*cos2psi - multiAMcoefficients->data[0]->a->data[fftnum]*sin2psi;
               REAL4 FplusX = multiAMcoefficients->data[jj]->a->data[fftnum]*cos2psi + multiAMcoefficients->data[jj]->b->data[fftnum]*sin2psi;
               REAL4 FcrossX = multiAMcoefficients->data[jj]->b->data[fftnum]*cos2psi - multiAMcoefficients->data[jj]->a->data[fftnum]*sin2psi;
               REAL4 cosi = *assumeNScosi;
               REAL4 onePlusCosiSqOver2 = 0.5*(1.0 + cosi*cosi);
               COMPLEX16 complexnumerator = crect(FplusX*onePlusCosiSqOver2, FcrossX*cosi);
               COMPLEX16 complexdenominator = crect(Fplus0*onePlusCosiSqOver2, Fcross0*cosi);
               if (cabs(complexdenominator)>1.0e-7) {
                  COMPLEX16 complexval = complexnumerator/complexdenominator;
                  detPhaseMag = fmin(cabs(complexval), 10.0);  //fmin here because sometimes the magnitude value is being divided by a small value and causing biases
                  detPhaseArg = gsl_sf_angle_restrict_pos(carg(complexval));
               } else {
                  detPhaseMag = 0.0;
                  detPhaseArg = 0.0;
               }
            }

            //When reference detector is off, don't add data that is less sensitive than the reference detector
            if (detPhaseMag!=0.0 && whichIFOsToBeUsed->data[0]==0 && detPhaseMag<1.0) {
               detPhaseMag = 0.0;
               detPhaseArg = 0.0;
            }

            //Compute delta values for the frequency difference for the Dirichlet kernel: fk*Tsft*[Tdot(jj) - 1] where jj=0,...,N-1
            if (params->vectorMath==1) {
               if (computeDelta0vals) {
                  XLAL_CHECK_NULL( sseScaleREAL8Vector(delta0vals, frequencyBins, TdotMinus1s->data[0]->data[fftnum]) == XLAL_SUCCESS, XLAL_EFUNC );
                  XLAL_CHECK_NULL( sseDDVectorMultiply(DirichletScaling0, delta0vals, delta0vals) == XLAL_SUCCESS, XLAL_EFUNC );
                  XLAL_CHECK_NULL( sseAddScalarToREAL8Vector(DirichletScaling0, DirichletScaling0, -1.0) == XLAL_SUCCESS, XLAL_EFUNC );
                  XLAL_CHECK_NULL( sseDDVectorMultiply(DirichletScaling0, DirichletScaling0, delta0vals) == XLAL_SUCCESS, XLAL_EFUNC );
               }
               XLAL_CHECK_NULL( sseScaleREAL8Vector(deltaXvals, frequencyBins, TdotMinus1s->data[jj]->data[fftnum]) == XLAL_SUCCESS, XLAL_EFUNC );
               XLAL_CHECK_NULL( sseDDVectorMultiply(DirichletScalingX, deltaXvals, deltaXvals) == XLAL_SUCCESS, XLAL_EFUNC );
               XLAL_CHECK_NULL( sseAddScalarToREAL8Vector(DirichletScalingX, DirichletScalingX, -1.0) == XLAL_SUCCESS, XLAL_EFUNC );
               XLAL_CHECK_NULL( sseDDVectorMultiply(DirichletScalingX, DirichletScalingX, deltaXvals) == XLAL_SUCCESS, XLAL_EFUNC );
            } else if (params->vectorMath==2) {
               if (computeDelta0vals) {
                  XLAL_CHECK_NULL( avxScaleREAL8Vector(delta0vals, frequencyBins, TdotMinus1s->data[0]->data[fftnum]) == XLAL_SUCCESS, XLAL_EFUNC );
                  XLAL_CHECK_NULL( avxDDVectorMultiply(DirichletScaling0, delta0vals, delta0vals) == XLAL_SUCCESS, XLAL_EFUNC );
                  XLAL_CHECK_NULL( avxAddScalarToREAL8Vector(DirichletScaling0, DirichletScaling0, -1.0) == XLAL_SUCCESS, XLAL_EFUNC );
                  XLAL_CHECK_NULL( avxDDVectorMultiply(DirichletScaling0, DirichletScaling0, delta0vals) == XLAL_SUCCESS, XLAL_EFUNC );
               }
               XLAL_CHECK_NULL( avxScaleREAL8Vector(deltaXvals, frequencyBins, TdotMinus1s->data[jj]->data[fftnum]) == XLAL_SUCCESS, XLAL_EFUNC );
               XLAL_CHECK_NULL( avxDDVectorMultiply(DirichletScalingX, deltaXvals, deltaXvals) == XLAL_SUCCESS, XLAL_EFUNC );
               XLAL_CHECK_NULL( avxAddScalarToREAL8Vector(DirichletScalingX, DirichletScalingX, -1.0) == XLAL_SUCCESS, XLAL_EFUNC );
               XLAL_CHECK_NULL( avxDDVectorMultiply(DirichletScalingX, DirichletScalingX, deltaXvals) == XLAL_SUCCESS, XLAL_EFUNC );
            } else {
               for (UINT4 kk=0; kk<frequencyBins->length; kk++) {
                  if (computeDelta0vals) {
                     delta0vals->data[kk] = frequencyBins->data[kk]*TdotMinus1s->data[0]->data[fftnum];
                     DirichletScaling0->data[kk] = delta0vals->data[kk]*(delta0vals->data[kk]*delta0vals->data[kk] - 1.0);
                  }
                  deltaXvals->data[kk] = frequencyBins->data[kk]*TdotMinus1s->data[jj]->data[fftnum];
                  DirichletScalingX->data[kk] = deltaXvals->data[kk]*(deltaXvals->data[kk]*deltaXvals->data[kk] - 1.0);
               }
            }
            if (computeDelta0vals) computeDelta0vals = 0;

            //Now finish the computation with the Dirichlet kernel ratio and final correction
            for (UINT4 kk=0; kk<sftcopy->data->length; kk++) {
               //COMPLEX16 DirichletRatio = conj(DirichletKernelLargeNHann(deltaXvals->data[kk])/DirichletKernelLargeNHann(delta0vals->data[kk]));
               //REAL8 detArgVal = (carg(DirichletRatio));
               REAL4 detArgVal = 0.0;
               COMPLEX8 DirichletRatio;
               if (DirichletKernalLargeNHannRatio(&DirichletRatio, delta0vals->data[kk], deltaXvals->data[kk], DirichletScaling0->data[kk]/DirichletScalingX->data[kk]) == 0) {
                  detArgVal = (REAL4)gsl_sf_angle_restrict_pos((REAL8)cargf(conjf(DirichletRatio)));
               } else {
                  detPhaseMag = 0.0;
                  detArgVal = 0.0;
               }

               //When signals are offset by a bin
               if (llabs((INT8)delta0vals->data[kk]-(INT8)deltaXvals->data[kk])>=1) detArgVal += LAL_PI;

               COMPLEX8 complexfactor = cpolarf(detPhaseMag, detArgVal+detPhaseArg-TwoPiFrequenciesTau->data[kk]);
               sftcopy->data->data[kk] *= complexfactor;
            }

            if (createSFT) {
               XLAL_CHECK_NULL( XLALCopySFT(&(combinedSFTs->data[ii]), sftcopy) == XLAL_SUCCESS, XLAL_EFUNC );
               createSFT = 0;
            }
            else XLAL_CHECK_NULL( XLALSFTAdd(&(combinedSFTs->data[ii]), sftcopy) == XLAL_SUCCESS, XLAL_EFUNC );
            XLALDestroySFT(sftcopy);
            whichSFTinMultiSFTvector->data[jj]++;
         }
      }
   }

   XLALDestroyINT4Vector(whichSFTinMultiSFTvector);
   XLALDestroyINT4Vector(whichIFOsToBeUsed);
   XLALDestroyREAL4VectorAligned(cosiVals);
   XLALDestroyREAL4VectorAligned(onePlusCosiSqOver2Vals);
   XLALDestroyREAL4VectorAligned(twoPsiVector);
   XLALDestroyREAL4VectorAligned(sin2psiVector);
   XLALDestroyREAL4VectorAligned(cos2psiVector);
   XLALDestroyREAL4VectorAligned(aVals);
   XLALDestroyREAL4VectorAligned(bVals);
   XLALDestroyREAL4VectorAligned(Fplus0s);
   XLALDestroyREAL4VectorAligned(Fcross0s);
   XLALDestroyREAL4VectorAligned(FplusXs);
   XLALDestroyREAL4VectorAligned(FcrossXs);
   destroyAlignedREAL8Vector(frequencies);
   destroyAlignedREAL8Vector(frequencyBins);
   destroyAlignedREAL8VectorArray(twoPiTauVals);
   destroyAlignedREAL8VectorArray(TdotMinus1s);
   destroyAlignedREAL8Vector(delta0vals);
   destroyAlignedREAL8Vector(deltaXvals);
   destroyAlignedREAL8Vector(TwoPiFrequenciesTau);
   XLALDestroyREAL4VectorAligned(Aplus0s);
   XLALDestroyREAL4VectorAligned(Across0s);
   XLALDestroyREAL4VectorAligned(AplusXs);
   XLALDestroyREAL4VectorAligned(AcrossXs);
   destroyAlignedREAL8Vector(DirichletScaling0);
   destroyAlignedREAL8Vector(DirichletScalingX);

   fprintf(stderr, "done\n");

   REAL4Vector *tfdata = NULL;
   XLAL_CHECK_NULL( (tfdata = convertSFTdataToPowers(combinedSFTs, params, 2.0/params->Tsft/(params->avesqrtSh*params->avesqrtSh))) != NULL, XLAL_EFUNC );
   XLALDestroySFTVector(combinedSFTs);

   return tfdata;
}


/**
 * Convert a SFTVector of sfts into powers
 * \param [in] sfts          Pointer to the SFTVector
 * \param [in] params        Pointer to the UserInput_t
 * \param [in] normalization Normalization value to prevent underflow
 * \return Pointer to REAL4Vector containing powers
 */
REAL4Vector * convertSFTdataToPowers(SFTVector *sfts, UserInput_t *params, REAL8 normalization)
{

   XLAL_CHECK_NULL( sfts != NULL && params != NULL, XLAL_EINVAL );

   fprintf(LOG, "Converting band to powers... ");
   fprintf(stderr, "Converting band to powers... ");

   UINT4 numffts = (UINT4)floor(params->Tobs/(params->Tsft-params->SFToverlap)-1);
   UINT4 sftlength = sfts->data[0].data->length;
   XLAL_CHECK_NULL( sftlength != 0, XLAL_EINVAL, "SFT has length of 0!\n" );

   //Allocate output
   REAL4Vector *tfdata = NULL;
   XLAL_CHECK_NULL( (tfdata = XLALCreateREAL4Vector(numffts*sftlength)) != NULL, XLAL_EFUNC );

   //Non-existant SFT counter
   UINT4 nonexistantsft = 0;

   REAL8 starttime = params->t0;

   //Load the data into the output vector, roughly normalizing as we go along from the input value
   REAL8 sqrtnorm = sqrt(normalization);
   for (UINT4 ii=0; ii<numffts; ii++) {
      if (ii-nonexistantsft < sfts->length) {
         SFTtype *sft = &(sfts->data[ii - nonexistantsft]);
         if (sft->epoch.gpsSeconds == (INT4)round(ii*(params->Tsft-params->SFToverlap)+starttime)) {
            for (UINT4 jj=0; jj<sftlength; jj++) {
               COMPLEX8 sftcoeff = sft->data->data[jj];
               tfdata->data[ii*sftlength + jj] = (REAL4)((sqrtnorm*crealf(sftcoeff))*(sqrtnorm*crealf(sftcoeff)) + (sqrtnorm*cimagf(sftcoeff))*(sqrtnorm*cimagf(sftcoeff)));  //power, normalized
            } /* for jj < sftLength */
         } else {
            memset(&(tfdata->data[ii*sftlength]), 0, sizeof(REAL4)*sftlength);
            nonexistantsft++;    //increment the nonexistantsft counter
         }
      } else {
         memset(&(tfdata->data[ii*sftlength]), 0, sizeof(REAL4)*sftlength);
         nonexistantsft++;    //increment the nonexistantsft counter
      }
   } /* for ii < numffts */

   fprintf(LOG, "done\n");
   fprintf(stderr, "done\n");

   fprintf(LOG, "Duty factor = %f\n", 1.0-(REAL4)nonexistantsft/(REAL4)numffts);
   fprintf(stderr, "Duty factor = %f\n", 1.0-(REAL4)nonexistantsft/(REAL4)numffts);

   REAL4 meanTFdata = calcMean(tfdata);
   REAL4 stddev = 0.0;
   XLAL_CHECK_NULL( calcStddev(&stddev, tfdata) == XLAL_SUCCESS, XLAL_EFUNC );
   fprintf(LOG, "TF before weighting, mean subtraction: mean = %g, std. dev. = %g\n", meanTFdata, stddev);
   fprintf(stderr, "TF before weighting, mean subtraction: mean = %g, std. dev. = %g\n", meanTFdata, stddev);

   return tfdata;

}


/**
 * Read in the data SFTs in one function
 * \param [in] params        Pointer to UserInput_t
 * \param [in] normalization Normalization value determined from expected noise background
 * \param [in] minfbin       Frequency value of the minimum frequency bin
 * \param [in] maxfbin       Frequency value of the maximum frequency bin
 * \return REAL4Vector of SFT powers
 */
REAL4Vector * readInSFTs(UserInput_t *params, REAL8 normalization, REAL8 minfbin, REAL8 maxfbin)
{

   XLAL_CHECK_NULL( params != NULL, XLAL_EINVAL );

   MultiSFTVector *sftvector = NULL;
   XLAL_CHECK_NULL( (sftvector = getMultiSFTVector(params, minfbin, maxfbin)) != NULL, XLAL_EFUNC );

   REAL4Vector *tfdata = NULL;
   XLAL_CHECK_NULL( (tfdata = convertSFTdataToPowers(sftvector->data[0], params, normalization)) != NULL, XLAL_EFUNC );

   //Destroy stuff
   XLALDestroyMultiSFTVector(sftvector);

   return tfdata;

} /* readInSFTs() */


/**
 * Create a list of timestamps from an SFTCatalog
 * \param [in] catalog Pointer to an SFTCatalog
 * \return Pointer to a list of GPS timestamps in a MultiLIGOTimeGPSVector
 */
MultiLIGOTimeGPSVector * getMultiTimeStampsFromSFTCatalog(SFTCatalog *catalog)
{

   XLAL_CHECK_NULL( catalog != NULL, XLAL_EINVAL );

   //Get the MultiSFTCatalogView
   MultiSFTCatalogView *catalogView = NULL;
   XLAL_CHECK_NULL( (catalogView = XLALGetMultiSFTCatalogView(catalog)) != NULL, XLAL_EFUNC );

   MultiLIGOTimeGPSVector *multiTimestamps = NULL;
   XLAL_CHECK_NULL( (multiTimestamps = XLALTimestampsFromMultiSFTCatalogView(catalogView)) != NULL, XLAL_EFUNC );

   XLALDestroyMultiSFTCatalogView(catalogView);

   return multiTimestamps;

}


/**
 * Create a list of timestamps from SFTs that might be a subset from those in an SFTCatalog, applying KS/Kuipers test if desired
 * \param [in] multiSFTvector Pointer to a MultiSFTVector
 * \param [in] params         Pointer to UserInput_t
 * \return Pointer to a list of GPS timestamps in a MultiLIGOTimeGPSVector
 */
MultiLIGOTimeGPSVector * getMultiTimeStampsFromSFTs(MultiSFTVector *multiSFTvector, UserInput_t *params)
{

   XLAL_CHECK_NULL( params != NULL, XLAL_EINVAL );

   MultiLIGOTimeGPSVector *multiTimestamps = NULL;
   XLAL_CHECK_NULL( (multiTimestamps = XLALCalloc(1, sizeof(*multiTimestamps))) != NULL, XLAL_ENOMEM );
   XLAL_CHECK_NULL( (multiTimestamps->data = XLALCalloc(multiSFTvector->length, sizeof(*multiTimestamps->data))) != NULL, XLAL_ENOMEM );
   multiTimestamps->length = multiSFTvector->length;

   if (params->markBadSFTs && !params->signalOnly) {
      REAL8 tfnormval = 2.0/(params->Tsft*(params->avesqrtSh*params->avesqrtSh));

      for (UINT4 ii=0; ii<multiSFTvector->length; ii++) {
         REAL4Vector *tfdata = NULL;
         XLAL_CHECK_NULL( (tfdata = convertSFTdataToPowers(multiSFTvector->data[ii], params, tfnormval)) != NULL, XLAL_EFUNC );

         INT4Vector *removeTheseSFTs = NULL;
         XLAL_CHECK_NULL( (removeTheseSFTs = markBadSFTs(tfdata, params)) != NULL, XLAL_EFUNC );

         INT4 numberofsfts = 0, sftlength = (INT4)multiSFTvector->data[ii]->data[0].data->length;
         for (UINT4 jj=0; jj<removeTheseSFTs->length; jj++) if (removeTheseSFTs->data[jj]==0 && tfdata->data[jj*sftlength]!=0.0) numberofsfts++;

         XLAL_CHECK_NULL( (multiTimestamps->data[ii] = XLALCreateTimestampVector(numberofsfts)) != NULL, XLAL_EFUNC );

         INT4 kk = 0;
         for (UINT4 jj=0; jj<removeTheseSFTs->length; jj++) {
            if (removeTheseSFTs->data[jj]==0 && tfdata->data[jj*sftlength]!=0.0) {
               XLALGPSSetREAL8(&(multiTimestamps->data[ii]->data[kk]), params->t0+jj*(params->Tsft-params->SFToverlap));
               XLAL_CHECK_NULL( xlalErrno == 0, XLAL_EFUNC );
               kk++;
            }
         }

         multiTimestamps->data[ii]->deltaT = params->Tsft;

         XLALDestroyINT4Vector(removeTheseSFTs);
         XLALDestroyREAL4Vector(tfdata);
      }
      params->markBadSFTs = 0;
   } else {
      XLAL_CHECK_NULL( (multiTimestamps = XLALExtractMultiTimestampsFromSFTs(multiSFTvector)) != NULL, XLAL_EFUNC );
   }

   return multiTimestamps;

}


/**
 * Create a list of timestamps from a segment list
 * \param [in] filenames  LALStringVector of filenames
 * \param [in] t0         Start time of the search
 * \param [in] Tsft       SFT timespan in seconds
 * \param [in] SFToverlap Overlap of the SFTs in seconds
 * \param [in] dur        Duration of the search in seconds
 * \return Pointer to a list of GPS timestamps in a MultiLIGOTimeGPSVector
 */
MultiLIGOTimeGPSVector * getMultiTimeStampsFromSegmentsFile(LALStringVector *filenames, REAL8 t0, REAL8 Tsft, REAL8 SFToverlap, REAL8 dur)
{

   XLAL_CHECK_NULL( filenames != NULL, XLAL_EINVAL );

   MultiLIGOTimeGPSVector *multiTimestamps = NULL;
   XLAL_CHECK_NULL( (multiTimestamps = XLALCalloc(1, sizeof(*multiTimestamps))) != NULL, XLAL_ENOMEM );
   multiTimestamps->length = filenames->length;
   XLAL_CHECK_NULL( (multiTimestamps->data = XLALCalloc(multiTimestamps->length, sizeof(*(multiTimestamps->data)))) != NULL, XLAL_ENOMEM );

   for (UINT4 ii=0; ii<filenames->length; ii++) {
      LIGOTimeGPSVector *timestamps = NULL;
      XLAL_CHECK_NULL( (timestamps = XLALTimestampsFromSegmentFile(filenames->data[ii], Tsft, SFToverlap, 0, 1)) != NULL, XLAL_EFUNC );

      //Shift to find the correct range of start to end SFTs
      //startSFTindex will be at the correct first SFT
      //endSFTindex will be one past the last SFT
      //Thus, we don't need to be inclusive when taking the difference between these two index values
      UINT4 startSFTindex = 0, endSFTindex = 0;
      for (UINT4 jj=0; jj<timestamps->length; jj++) {
         REAL8 timestamp = XLALGPSGetREAL8(&(timestamps->data[jj]));
         if (timestamp<t0) startSFTindex++;
         else if (timestamp>=t0 && timestamp+Tsft<=t0+dur) endSFTindex++;
         else if (timestamp+Tsft>t0+dur) break;
      }

      XLAL_CHECK_NULL( (multiTimestamps->data[ii] = XLALCreateTimestampVector(endSFTindex - startSFTindex)) != NULL, XLAL_EFUNC );
      multiTimestamps->data[ii]->deltaT = Tsft;
      for (UINT4 jj=0; jj<multiTimestamps->data[ii]->length; jj++) multiTimestamps->data[ii]->data[jj] = timestamps->data[startSFTindex + jj];

      XLALDestroyTimestampVector(timestamps);
   }

   return multiTimestamps;

}


/* Critical values of KS test (from Bickel and Doksum). Does not apply directly (mean determined from distribution)
alpha=0.01
n       10      20      30      40      50      60      80      n>80
        .489    .352    .290    .252    .226    .207    .179    1.628/(sqrt(n)+0.12+0.11/sqrt(n))

alpha=0.05
n       10      20      30      40      50      60      80      n>80
        .409    .294    .242    .210    .188    .172    .150    1.358/(sqrt(n)+0.12+0.11/sqrt(n))

alpha=0.1 (E.G derived using root finding)
n       10      20      30      40      50      60      80      n>80
        .369    .265    .218    .189    .170    .155    .135    1.224/(sqrt(n)+0.12+0.11/sqrt(n))

alpha=0.2 (E.G derived using root finding)
n       10      20      30      40      50      60      80      n>80
                                                                1.073/(sqrt(n)+0.12+0.11/sqrt(n))

alpha=0.5 (E.G derived using root finding)
n       10      20      30      40      50      60      80      n>80
        .249    .179    .147    .128    .115    .105    .091    0.828/(sqrt(n)+0.12+0.11/sqrt(n))

alpha=0.9 (E.G derived using root finding)
n                                                               n>80
                                                                0.571/(sqrt(n)+0.12+0.11/sqrt(n))

Critical values of Kuiper's test using root finding by E.G.
alpha=0.05
n                                                               n>80
                                                                1.747/(sqrt(n)+0.155+0.24/sqrt(n))

alpha=0.1
n                                                               n>80
                                                                1.620/(sqrt(n)+0.155+0.24/sqrt(n))

alpha=0.2
n                                                               n>80
                                                                1.473/(sqrt(n)+0.155+0.24/sqrt(n))
*/
/**
 * Mark the non-Gaussian SFTs using K-S and Kuiper's tests
 * \param [in] tfdata Pointer to REAL4Vector of SFT powers
 * \param [in] params Pointer to UserInput_t
 * \return Pointer to an INT4Vector with marked SFTs to be removed with 1 and keep SFT with 0
 */
INT4Vector * markBadSFTs(REAL4Vector *tfdata, UserInput_t *params)
{

   XLAL_CHECK_NULL( tfdata != NULL && params != NULL, XLAL_EINVAL );

   fprintf(stderr, "Marking bad SFTs... ");

   INT4 numffts = (INT4)floor(params->Tobs/(params->Tsft-params->SFToverlap)-1);    //Number of FFTs
   //INT4 numfbins = (INT4)(round(params->fspan*params->Tsft+2.0*params->dfmax*params->Tsft)+12+1)+2*params->maxbinshift+params->blksize-1;     //Number of frequency bins
   INT4 numfbins = (INT4)tfdata->length/numffts;

   //Allocate output data vector and a single SFT data vector
   INT4Vector *output = NULL;
   XLAL_CHECK_NULL( (output = XLALCreateINT4Vector(numffts)) != NULL, XLAL_EFUNC );
   memset(output->data, 0, sizeof(INT4)*output->length);
   REAL4Vector *tempvect = NULL;
   XLAL_CHECK_NULL( (tempvect = XLALCreateREAL4Vector(numfbins)) != NULL, XLAL_EFUNC );

   //Do the KS and Kuiper test on each SFT
   REAL8 ksthreshold = 1.358/(sqrt(numfbins)+0.12+0.11/sqrt(numfbins));
   //REAL8 ksthreshold = 1.224/(sqrt(numfbins)+0.12+0.11/sqrt(numfbins));  //This is a tighter restriction
   //REAL8 ksthreshold = 1.073/(sqrt(numfbins)+0.12+0.11/sqrt(numfbins));  //This is an even tighter restriction
   REAL8 kuiperthreshold = 1.747/(sqrt(numfbins)+0.155+0.24/sqrt(numfbins));
   //REAL8 kuiperthreshold = 1.620/(sqrt(numfbins)+0.155+0.24/sqrt(numfbins));  //This is a tighter restriction
   //REAL8 kuiperthreshold = 1.473/(sqrt(numfbins)+0.155+0.24/sqrt(numfbins));  //This is an even tighter restriction
   INT4 badsfts = 0, totalsfts = 0;
   REAL8 kstest = 0.0, kuipertest = 0.0;
   for (INT4 ii=0; ii<numffts; ii++) {
      if (tfdata->data[ii*numfbins]!=0.0) {
         totalsfts++;
         memcpy(tempvect->data, &(tfdata->data[ii*numfbins]), sizeof(REAL4)*tempvect->length);
         XLAL_CHECK_NULL( ks_test_exp(&kstest, tempvect) == XLAL_SUCCESS, XLAL_EFUNC );
         XLAL_CHECK_NULL( kuipers_test_exp(&kuipertest, tempvect) == XLAL_SUCCESS, XLAL_EFUNC );
         if (kstest>ksthreshold || kuipertest>kuiperthreshold) {
            output->data[ii] = 1;
            badsfts++;
         }
      }
   }

   //Destroy stuff
   XLALDestroyREAL4Vector(tempvect);

   fprintf(stderr, "done\n");
   fprintf(stderr, "Fraction excluded in K-S and Kuiper's tests = %f\n", (REAL4)badsfts/(REAL4)totalsfts);

   return output;

}


/**
 * Remove the marked SFTs as bad by setting values to 0
 * \param [in,out] tfdata  Pointer to REAL4Vector of SFT powers
 * \param [in]     badsfts Poienter to INT4Vector of bad SFTs
 */
void removeBadSFTs(REAL4Vector *tfdata, INT4Vector *badsfts)
{

   XLAL_CHECK_VOID( tfdata != NULL && badsfts != NULL, XLAL_EINVAL );

   fprintf(stderr, "Removing bad SFTs... ");

   INT4 numfbins_tfdata = tfdata->length/badsfts->length;

   for (UINT4 ii=0; ii<badsfts->length; ii++) if (badsfts->data[ii]==1) memset(&(tfdata->data[ii*numfbins_tfdata]), 0, sizeof(REAL4)*numfbins_tfdata);

   fprintf(stderr, "done.\n");

}

MultiSFTVector * generateSFTdata(UserInput_t *uvar, MultiLALDetector *detectors, EphemerisData *edat, INT4 maxbinshift, gsl_rng *rng)
{
   MultiSFTVector *multiSFTvector = NULL;
   
   //Determine band size to get the SFT data (remember to get extra bins because of the running median and the bin shifts due to detector velocity) with nudge of 0.1/Tsft for rounding issues
   REAL8 minfbin = round(uvar->fmin*uvar->Tsft - uvar->dfmax*uvar->Tsft - 0.5*(uvar->blksize-1) - (REAL8)(maxbinshift) - 6.0)/uvar->Tsft + 0.1/uvar->Tsft;
   REAL8 maxfbin = round((uvar->fmin + uvar->fspan)*uvar->Tsft + uvar->dfmax*uvar->Tsft + 0.5*(uvar->blksize-1) + (REAL8)(maxbinshift) + 6.0)/uvar->Tsft - 0.1/uvar->Tsft;
   
   //Start by getting timestamps or creating them
   MultiLIGOTimeGPSVector *multiTimestamps = NULL;
   if (XLALUserVarWasSet(&uvar->inputSFTs)) {
      SFTCatalog *catalog = NULL;
      XLAL_CHECK_NULL( (catalog = findSFTdata(uvar)) != NULL, XLAL_EFUNC );
      MultiSFTVector *tmpsftvector = NULL;
      XLAL_CHECK_NULL( (tmpsftvector = extractSFTband(catalog, minfbin, maxfbin)) != NULL, XLAL_EFUNC );
      XLAL_CHECK_NULL( (multiTimestamps = getMultiTimeStampsFromSFTs(tmpsftvector, uvar)) != NULL, XLAL_EFUNC );
      XLALDestroyMultiSFTVector(tmpsftvector);
      XLALDestroySFTCatalog(catalog);
   } else if (XLALUserVarWasSet(&uvar->timestampsFile)) {
      XLAL_CHECK_NULL( (multiTimestamps = XLALReadMultiTimestampsFiles(uvar->timestampsFile)) != NULL, XLAL_EFUNC );
      for (UINT4 ii=0; ii<multiTimestamps->length; ii++) multiTimestamps->data[ii]->deltaT = uvar->Tsft;
   }
   else if (XLALUserVarWasSet(&uvar->segmentFile)) XLAL_CHECK_NULL( (multiTimestamps = getMultiTimeStampsFromSegmentsFile(uvar->segmentFile, uvar->t0, uvar->Tsft, uvar->SFToverlap, uvar->Tobs)) != NULL, XLAL_EFUNC );
   else {
      LIGOTimeGPS tStart;
      XLALGPSSetREAL8 ( &tStart, uvar->t0 );
      XLAL_CHECK_NULL( xlalErrno == 0, XLAL_EFUNC, "XLALGPSSetREAL8 failed\n" );
      XLAL_CHECK_NULL( (multiTimestamps = XLALMakeMultiTimestamps(tStart, uvar->Tobs, uvar->Tsft, uvar->SFToverlap, detectors->length)) != NULL, XLAL_EFUNC );
   }

   //TwoSpect to analyze:
   REAL8 TwoSpectFmin = round(uvar->fmin*uvar->Tsft - uvar->dfmax*uvar->Tsft - 0.5*(uvar->blksize-1) - (REAL8)maxbinshift - 6.0)/uvar->Tsft;
   REAL8 TwoSpectBand = round(uvar->fspan*uvar->Tsft + 2.0*uvar->dfmax*uvar->Tsft + (uvar->blksize-1) + (REAL8)(2.0*maxbinshift) + 12.0)/uvar->Tsft;

   //Setup the MFD data parameters
   CWMFDataParams DataParams;
   if (XLALUserVarWasSet(&uvar->injFmin) && XLALUserVarWasSet(&uvar->injBand) && uvar->injFmin<=TwoSpectFmin && uvar->injFmin+uvar->injBand>=TwoSpectFmin+TwoSpectBand) {
      DataParams.fMin = uvar->injFmin;
      DataParams.Band = uvar->injBand;
   } else {
      DataParams.fMin = TwoSpectFmin;
      DataParams.Band = TwoSpectBand;
   }
   DataParams.multiIFO.length = detectors->length;
   for (UINT4 ii=0; ii<detectors->length; ii++) DataParams.multiIFO.sites[ii] = detectors->sites[ii];
   DataParams.multiNoiseFloor.length = detectors->length;
   for (UINT4 ii=0; ii<detectors->length; ii++) DataParams.multiNoiseFloor.sqrtSn[ii] = 0.0;
   DataParams.multiTimestamps = *multiTimestamps;
   DataParams.randSeed = uvar->injRandSeed;
   DataParams.SFTWindowType = "Hann";
   DataParams.SFTWindowBeta = 0;

   MultiSFTVector *signalSFTs = NULL;
   PulsarParamsVector *injectionSources = NULL;
   //If injection sources then read them and make signal sfts
   if (XLALUserVarWasSet(&uvar->injectionSources)) {
      XLAL_CHECK_NULL( (injectionSources =  XLALPulsarParamsFromUserInput(uvar->injectionSources)) != NULL, XLAL_EFUNC );

      fprintf(stderr, "Generating signal SFTs with %d signals... ", (INT4)injectionSources->length);

      if (!uvar->signalOnly) XLAL_CHECK_NULL( XLALCWMakeFakeMultiData(&signalSFTs, NULL, injectionSources, &DataParams, edat) == XLAL_SUCCESS, XLAL_EFUNC );
      else XLAL_CHECK_NULL( XLALCWMakeFakeMultiData(&multiSFTvector, NULL, injectionSources, &DataParams, edat) == XLAL_SUCCESS, XLAL_EFUNC );

      if (DataParams.fMin != TwoSpectFmin) XLAL_CHECK_NULL( XLALMultiSFTVectorResizeBand(signalSFTs, TwoSpectFmin, TwoSpectBand) == XLAL_SUCCESS, XLAL_EFUNC );

      fprintf(stderr, "done\n");
   } // if there are injections

   //If not signal only, create sfts that include noise or extract a band from real data
   if (!uvar->signalOnly) {
      if (uvar->gaussNoiseWithSFTgaps || XLALUserVarWasSet(&uvar->timestampsFile) || XLALUserVarWasSet(&uvar->segmentFile) || !XLALUserVarWasSet(&uvar->inputSFTs)) {
         for (UINT4 ii=0; ii<detectors->length; ii++) DataParams.multiNoiseFloor.sqrtSn[ii] = uvar->avesqrtSh;

         fprintf(stderr, "Generating noise SFTs... ");

         XLAL_CHECK_NULL( XLALCWMakeFakeMultiData(&multiSFTvector, NULL, NULL, &DataParams, edat) == XLAL_SUCCESS, XLAL_EFUNC );

         if (DataParams.fMin != TwoSpectFmin) XLAL_CHECK_NULL( XLALMultiSFTVectorResizeBand(multiSFTvector, TwoSpectFmin, TwoSpectBand) == XLAL_SUCCESS, XLAL_EFUNC );

         fprintf(stderr, "done\n");
      } else {
         SFTCatalog *catalog = NULL;
         XLAL_CHECK_NULL( (catalog = findSFTdata(uvar)) != NULL, XLAL_EFUNC );
         XLAL_CHECK_NULL( (multiSFTvector = extractSFTband(catalog, minfbin, maxfbin)) != NULL, XLAL_EFUNC );
         XLALDestroySFTCatalog(catalog);
      }
   } // if not signal only SFTs

   //Add the SFT vectors together
   if (XLALUserVarWasSet(&uvar->injectionSources) && !uvar->signalOnly) {
      XLAL_CHECK_NULL( XLALMultiSFTVectorAdd(multiSFTvector, signalSFTs) == XLAL_SUCCESS, XLAL_EFUNC );
      XLALDestroyMultiSFTVector(signalSFTs);
   }

   //If printing the data outputs, then do that here
   if ((XLALUserVarWasSet(&uvar->printSignalData) || XLALUserVarWasSet(&uvar->printMarginalizedSignalData)) && XLALUserVarWasSet(&uvar->injectionSources)) {
      for (UINT4 ii=0; ii<detectors->length; ii++) DataParams.multiNoiseFloor.sqrtSn[ii] = 0.0;
      PulsarParamsVector *oneSignal = NULL;
      XLAL_CHECK_NULL( (oneSignal = XLALCreatePulsarParamsVector(1)) != NULL, XLAL_EFUNC );

      FILE *SIGNALOUT = NULL, *MARGINALIZEDSIGNALOUT = NULL;
      if (XLALUserVarWasSet(&uvar->printSignalData)) XLAL_CHECK_NULL( (SIGNALOUT = fopen(uvar->printSignalData, "w")) != NULL, XLAL_EIO, "Failed to open %s for writing\n", uvar->printSignalData );
      if (XLALUserVarWasSet(&uvar->printMarginalizedSignalData)) XLAL_CHECK_NULL( (MARGINALIZEDSIGNALOUT = fopen(uvar->printMarginalizedSignalData, "w")) != NULL, XLAL_EIO, "Failed to open %s for writing\n", uvar->printMarginalizedSignalData );

      for (UINT4 ii=0; ii<injectionSources->length; ii++) {
         memcpy(oneSignal->data, &(injectionSources->data[ii]), sizeof(injectionSources->data[0]));
         if (XLALUserVarWasSet(&uvar->printSignalData)) {
            MultiSFTVector *oneSignalSFTs = NULL;
            XLAL_CHECK_NULL( XLALCWMakeFakeMultiData(&oneSignalSFTs, NULL, oneSignal, &DataParams, edat) == XLAL_SUCCESS, XLAL_EFUNC );
            if (DataParams.fMin != TwoSpectFmin) XLAL_CHECK_NULL( XLALMultiSFTVectorResizeBand(oneSignalSFTs, TwoSpectFmin, TwoSpectBand) == XLAL_SUCCESS, XLAL_EFUNC );

            REAL8Vector *aveSFTsPower = NULL;
            XLAL_CHECK_NULL( (aveSFTsPower = XLALCreateREAL8Vector(multiSFTvector->data[0]->data->data->length)) != NULL, XLAL_EFUNC );
            memset(aveSFTsPower->data, 0, sizeof(REAL8)*aveSFTsPower->length);

            for (UINT4 jj=0; jj<oneSignalSFTs->data[0]->length; jj++) {
               SFTtype *sft = &(oneSignalSFTs->data[0]->data[jj]);
               for (UINT4 kk=0; kk<aveSFTsPower->length; kk++) {
                  REAL8 powerval = 2.0*(creal(sft->data->data[kk])*creal(sft->data->data[kk]) + cimag(sft->data->data[kk])*cimag(sft->data->data[kk]))/uvar->Tsft;
                  aveSFTsPower->data[kk] += powerval;
               }
            }
            for (UINT4 jj=0; jj<aveSFTsPower->length; jj++) fprintf(SIGNALOUT,"%.9g %.9g\n", DataParams.fMin+jj/uvar->Tsft, aveSFTsPower->data[jj]/multiSFTvector->data[0]->length);
            XLALDestroyMultiSFTVector(oneSignalSFTs);
            XLALDestroyREAL8Vector(aveSFTsPower);
         }
         if (XLALUserVarWasSet(&uvar->printMarginalizedSignalData)) {
            REAL8Vector *marginalizedSignalData = NULL;
            XLAL_CHECK_NULL ( (marginalizedSignalData = XLALCreateREAL8Vector(multiSFTvector->data[0]->data->data->length)) != NULL, XLAL_EFUNC );
            memset(marginalizedSignalData->data, 0, sizeof(REAL8)*marginalizedSignalData->length);
            for (UINT4 jj=0; jj<300; jj++) {
               oneSignal->data[0].Amp.cosi = 2.0*gsl_rng_uniform(rng) - 1.0;
               oneSignal->data[0].Amp.psi = LAL_TWOPI*gsl_rng_uniform(rng);
               oneSignal->data[0].Amp.phi0 = LAL_TWOPI*gsl_rng_uniform(rng);
               oneSignal->data[0].Doppler.argp = LAL_TWOPI*gsl_rng_uniform(rng);
               MultiSFTVector *oneSignalSFTs = NULL;
               XLAL_CHECK_NULL( XLALCWMakeFakeMultiData(&oneSignalSFTs, NULL, oneSignal, &DataParams, edat) == XLAL_SUCCESS, XLAL_EFUNC );
               if (DataParams.fMin != TwoSpectFmin) XLAL_CHECK_NULL( XLALMultiSFTVectorResizeBand(oneSignalSFTs, TwoSpectFmin, TwoSpectBand) == XLAL_SUCCESS, XLAL_EFUNC );

               for (UINT4 kk=0; kk<oneSignalSFTs->data[0]->length; kk++) {
                  SFTtype *sft = &(oneSignalSFTs->data[0]->data[kk]);
                  for (UINT4 ll=0; ll<marginalizedSignalData->length; ll++) marginalizedSignalData->data[ll] += (2.0*(creal(sft->data->data[ll])*creal(sft->data->data[ll]) + cimag(sft->data->data[ll])*cimag(sft->data->data[ll]))/uvar->Tsft);
               }
               XLALDestroyMultiSFTVector(oneSignalSFTs);
            } //Loop over trials
            for (UINT4 jj=0; jj<marginalizedSignalData->length; jj++) {
               marginalizedSignalData->data[jj] /= 300.0*multiSFTvector->data[0]->length;
               fprintf(MARGINALIZEDSIGNALOUT,"%.9g %.9g\n", DataParams.fMin+jj/uvar->Tsft, marginalizedSignalData->data[jj]);
            }
            XLALDestroyREAL8Vector(marginalizedSignalData);
         } //If printing marginalized data
      } //loop over the number of injected sources
      memset(oneSignal->data, 0, sizeof(injectionSources->data[0]));
      XLALDestroyPulsarParamsVector(oneSignal);
      if (XLALUserVarWasSet(&uvar->printSignalData)) fclose(SIGNALOUT);
      if (XLALUserVarWasSet(&uvar->printMarginalizedSignalData)) fclose(MARGINALIZEDSIGNALOUT);
   } //end printing data

   if (XLALUserVarWasSet(&uvar->injectionSources)) XLALDestroyPulsarParamsVector(injectionSources);
   XLALDestroyMultiTimestamps(multiTimestamps);

   return multiSFTvector;
}


/**
 * Slide the time-frequency data to account for detector motion
 * \param [out] output    Pointer to REAL4Vector of SFT powers that have been corrected
 * \param [in]  params    Pointer to UserInput_t
 * \param [in]  tfdata    Pointer to REAL4Vector of SFT powers
 * \param [in]  binshifts Pointer to INT4Vector of bin shift values
 * \return Status value
 */
INT4 slideTFdata(REAL4Vector *output, UserInput_t *params, REAL4Vector *tfdata, INT4Vector *binshifts)
{

   XLAL_CHECK( output != NULL && params != NULL && tfdata != NULL && binshifts != NULL, XLAL_EINVAL );

   fprintf(stderr, "Sliding TF data... ");

   UINT4 numffts = (UINT4)floor(params->Tobs/(params->Tsft-params->SFToverlap)-1);
   UINT4 numfbins = output->length/numffts;
   UINT4 maxbinshift = (tfdata->length/numffts - numfbins)/2;

   for (UINT4 ii=0; ii<numffts; ii++) {
      XLAL_CHECK( binshifts->data[ii]<=(INT4)maxbinshift, XLAL_EFAILED, "SFT slide value %d is greater than maximum value predicted (%d)", binshifts->data[ii], maxbinshift );
      memcpy(&(output->data[ii*numfbins]), &(tfdata->data[ii*(numfbins+2*maxbinshift) + maxbinshift + binshifts->data[ii]]), sizeof(REAL4)*numfbins);
   }

   fprintf(stderr, "done\n");

   return XLAL_SUCCESS;

} /* slideTFdata() */


/**
 * Determine the running mean of each SFT
 * \param [out] output   Pointer to REAL4Vector of running mean values of each SFT
 * \param [in]  tfdata   Pointer to REAL4Vector of SFT powers
 * \param [in]  numffts  Number of SFTs in the observation time
 * \param [in]  numfbins Number of frequency bins
 * \param [in]  blksize  Number of bins in the running median
 * \return Status value
 */
INT4 tfRngMeans(REAL4Vector *output, REAL4Vector *tfdata, INT4 numffts, INT4 numfbins, INT4 blksize)
{

   XLAL_CHECK( output != NULL && tfdata != NULL && numffts > 0  && numfbins > 0 && blksize > 0, XLAL_EINVAL );

   fprintf(LOG, "Assessing SFT background... ");
   fprintf(stderr, "Assessing SFT background... ");

   LALStatus XLAL_INIT_DECL(status);
   REAL8 bias;
   INT4 totalfbins = numfbins + blksize - 1;

   //Blocksize of running median
   LALRunningMedianPar block = {blksize};

   //Running median bias calculation
   if (blksize<1000) {
      bias = XLALRngMedBias(blksize);
      XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
   } else  bias = LAL_LN2;
   //REAL8 invbias = 1.0/(bias*1.0099993480677538);  //StackSlide normalization for 101 bins
   REAL8 invbias = 1.0/bias;

   //Allocate for a single SFT data and the medians out of each SFT
   REAL4Vector *inpsd = NULL, *mediansout = NULL;
   XLAL_CHECK( (inpsd = XLALCreateREAL4Vector(totalfbins)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (mediansout = XLALCreateREAL4Vector(numfbins)) != NULL, XLAL_EFUNC );

   //Now do the running median
   for (INT4 ii=0; ii<numffts; ii++) {
      //If the SFT values were not zero, then compute the running median
      if (tfdata->data[ii*totalfbins]!=0.0) {
         //Determine running median value, convert to mean value
         memcpy(inpsd->data, &(tfdata->data[ii*inpsd->length]), sizeof(REAL4)*inpsd->length);

         //calculate running median
         LALSRunningMedian2(&status, mediansout, inpsd, block);
         XLAL_CHECK( status.statusCode == 0, XLAL_EFUNC );

         //Now make the output medians into means by multiplying by 1/bias
         for (UINT4 jj=0; jj<mediansout->length; jj++) output->data[ii*numfbins + jj] = (REAL4)(mediansout->data[jj]*invbias);
      } else {
         //Otherwise, set means to zero
         for (UINT4 jj=0; jj<mediansout->length; jj++) output->data[ii*numfbins + jj] = 0.0;
      }
   } /* for ii < numffts */

   //Destroy stuff
   XLALDestroyREAL4Vector(inpsd);
   XLALDestroyREAL4Vector(mediansout);

   fprintf(LOG, "done\n");
   fprintf(stderr, "done\n");
   fprintf(stderr,"Mean of running means = %g\n", calcMean(output));

   return XLAL_SUCCESS;

} /* tfRngMeans() */


/**
 * Subtract running mean values from the SFT data, modifying input time-frequency data
 * \param [in,out] tfdata   Pointer to REAL4Vector time-frequency data (modified by this function!)
 * \param [in]     rngMeans Pointer to REAL4Vector of running mean values
 * \param [in]     numffts  Number of SFTs from observation time
 * \param [in]     numfbins Number of SFT frequency bins
 * \return Status value
 */
INT4 tfMeanSubtract(REAL4Vector *tfdata, REAL4Vector *rngMeans, INT4 numffts, INT4 numfbins)
{

   XLAL_CHECK( tfdata != NULL && rngMeans != NULL && numffts > 0 && numfbins > 0, XLAL_EINVAL );

   fprintf(stderr, "Subtracting expected background... ");

   for (INT4 ii=0; ii<numffts; ii++) if (rngMeans->data[ii*numfbins]!=0.0) for (INT4 jj=0; jj<numfbins; jj++) tfdata->data[ii*numfbins+jj] -= rngMeans->data[ii*numfbins+jj];

   fprintf(stderr, "done\n");

   return XLAL_SUCCESS;

} /* tfMeanSubtract() */


/**
 * Weight the SFTs based on antenna pattern and noise variance (Equation 11, assuming the input time-frequency data is already mean subtracted)
 * \param [out] output                    Pointer to REAL4Vector of mean subtracted, noise and antenna pattern weighted SFTs
 * \param [in]  tfdata                    Pointer to REAL4Vector of mean subtracted SFTs
 * \param [in]  rngMeans                  Pointer to REAL4Vector of running mean values
 * \param [in]  antPatternWeights         Pointer to REAL4VectorAligned of antenna pattern weights
 * \param [in]  indexValuesOfExistingSFTs Pointer to INT4Vector of the index values of the existing SFTs
 * \param [in]  params                    Pointer to UserInput_t
 * \return Status value
 */
INT4 tfWeight(REAL4Vector *output, REAL4Vector *tfdata, REAL4Vector *rngMeans, REAL4VectorAligned *antPatternWeights, INT4Vector *indexValuesOfExistingSFTs, UserInput_t *params)
{

   XLAL_CHECK( output != NULL && tfdata != NULL && rngMeans != NULL && antPatternWeights != NULL && indexValuesOfExistingSFTs != NULL && params != NULL, XLAL_EINVAL );

   fprintf(stderr, "Applying weighting to SFT data... ");

   UINT4 numffts = antPatternWeights->length;
   UINT4 numfbins = tfdata->length/numffts;

   //Initially set output to zero
   memset(output->data, 0, sizeof(REAL4)*output->length);

   REAL4VectorAligned *antweightssq = NULL;
   REAL4Vector *rngMeanssq = NULL;
   XLAL_CHECK( (antweightssq = XLALCreateREAL4VectorAligned(numffts, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (rngMeanssq = XLALCreateREAL4Vector(numffts)) != NULL, XLAL_EFUNC );

   //If user specified that SFTs contain signal only, then the values of the backgrnd vector will be zeros.
   //We must set them equal to 1.0 here, then return them to zeros at the end
   if (params->signalOnly) {
      for (UINT4 ii=0; ii<indexValuesOfExistingSFTs->length; ii++) {
         for (UINT4 jj=0; jj<numfbins; jj++) rngMeans->data[numfbins*indexValuesOfExistingSFTs->data[ii] + jj] = params->avesqrtSh;
      }
   }

   //User specifies whether to use SSE to do the multiplication or not
   if (params->vectorMath==1) XLAL_CHECK( sseSSVectorMultiply(antweightssq, antPatternWeights, antPatternWeights) == XLAL_SUCCESS, XLAL_EFUNC );
   else if (params->vectorMath==2) XLAL_CHECK( avxSSVectorMultiply(antweightssq, antPatternWeights, antPatternWeights) == XLAL_SUCCESS, XLAL_EFUNC );
   else {
      XLAL_CHECK( XLALSSVectorMultiply((REAL4Vector*)antweightssq, (REAL4Vector*)antPatternWeights, (REAL4Vector*)antPatternWeights) != NULL, XLAL_EFUNC );
      XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
   }

   //Loop through the SFT frequency bins and weight the data
   for (UINT4 ii=0; ii<numfbins; ii++) {

      fastSSVectorMultiply_with_stride_and_offset(rngMeanssq, rngMeans, rngMeans, numfbins, numfbins, ii, ii);

      //If noiseWeightOff is given, then set all the noise weights to be 1.0
      if (params->noiseWeightOff) for (UINT4 jj=0; jj<rngMeanssq->length; jj++) if (rngMeanssq->data[jj]!=0.0) rngMeanssq->data[jj] = 1.0;

      //Get sum of antenna pattern weight/variances for each frequency bin as a function of time (only for existant SFTs)
      REAL8 sumofweights = determineSumOfWeights((REAL4Vector*)antweightssq, rngMeanssq);
      REAL8 invsumofweights = 1.0/sumofweights;

      //Now do noise weighting, antenna pattern weighting
      for (UINT4 jj=0; jj<indexValuesOfExistingSFTs->length; jj++) {
         output->data[indexValuesOfExistingSFTs->data[jj]*numfbins+ii] = (REAL4)(invsumofweights*antPatternWeights->data[indexValuesOfExistingSFTs->data[jj]]*tfdata->data[indexValuesOfExistingSFTs->data[jj]*numfbins+ii]/rngMeanssq->data[indexValuesOfExistingSFTs->data[jj]]);
      } /* for jj < indexValuesOfExisitingSFTs->length */
   } /* for ii < numfbins */

   //Remember to reset the backgrnd vector to zero
   if (params->signalOnly) memset(rngMeans->data, 0, sizeof(REAL4)*rngMeans->length);

   //Destroy stuff
   XLALDestroyREAL4VectorAligned(antweightssq);
   XLALDestroyREAL4Vector(rngMeanssq);

   fprintf(stderr,"done\n");

   return XLAL_SUCCESS;

} /* tfWeight() */


/**
 * Determine the sum of the weights
 * \param [in] antweightssq Antenna pattern weights squared
 * \param [in] rngMeanssq   Running mean values squared
 * \return Sum of the weights
 */
REAL8 determineSumOfWeights(REAL4Vector *antweightssq, REAL4Vector *rngMeanssq)
{

   XLAL_CHECK_REAL8( antweightssq != NULL && rngMeanssq != NULL, XLAL_EINVAL );

   REAL8 sumofweights = 0.0;
   for (UINT4 ii=0; ii<antweightssq->length; ii++) if (rngMeanssq->data[ii] != 0.0) sumofweights += antweightssq->data[ii]/rngMeanssq->data[ii];

   return sumofweights;

} /* determineSumOfWeights */


/**
 * Determine if the SFTs are existing based on whether the first data point of power in each SFT is 0
 * \param [in] tfdata   Pointer to REAL4Vector of time-frequency data
 * \param [in] numffts  Number of SFTs from the observation time
 * \return Pointer to INT4Vector containing 0 for non-present SFT or 1 for present SFT
 */
INT4Vector * existingSFTs(REAL4Vector *tfdata, UINT4 numffts)
{

   XLAL_CHECK_NULL( tfdata != NULL && numffts > 0, XLAL_EINVAL );

   fprintf(stderr, "Determining existing SFTs... ");

   UINT4 numfbins = tfdata->length/numffts;     //Number of frequency bins

   INT4Vector *sftexist = NULL;
   XLAL_CHECK_NULL( (sftexist = XLALCreateINT4Vector(numffts)) != NULL, XLAL_EFUNC );

   for (UINT4 ii=0; ii<numffts; ii++) {
      if (tfdata->data[ii*numfbins] == 0.0) sftexist->data[ii] = 0;
      else sftexist->data[ii] = 1;
   }

   fprintf(stderr, "done\n");

   return sftexist;

} /* existingSFTs() */
