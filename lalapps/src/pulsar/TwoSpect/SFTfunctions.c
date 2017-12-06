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
#include <lal/GenerateSpinOrbitCW.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include "SFTfunctions.h"
#include "TwoSpectSpecFunc.h"
#include "statistics.h"

/**
 * Find the SFT data specified by user input
 * \param [in] params Pointer to the UserInput_t
 * \return Pointer to SFTCatalog
 */
SFTCatalog * findSFTdata(const UserInput_t *params)
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

} // findSFTdata()

/**
 * Extract the SFT coefficients from the band of interest
 * \param [in] catalog Pointer to the SFTCatalog
 * \param [in] minfbin Frequency value of the minimum frequency bin
 * \param [in] maxfbin Frequency value of the maximum frequency bin
 * \return Pointer to MultiSFTVector
 */
MultiSFTVector * extractSFTband(const SFTCatalog *catalog, const REAL8 minfbin, const REAL8 maxfbin)
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

} // extractSFTband()

/**
 * Get a MultiSFTVector from the user-input values
 * \param [in] params Pointer to UserInput_t
 * \param [in] minfbin Frequency value of the minimum frequency bin
 * \param [in] maxfbin Frequency value of the maximum frequency bin
 * \return Pointer to MultiSFTVector
 */
MultiSFTVector * getMultiSFTVector(UserInput_t *params, const REAL8 minfbin, const REAL8 maxfbin)
{
   //Get catalog of SFTs
   SFTCatalog *catalog = NULL;
   XLAL_CHECK_NULL( (catalog = findSFTdata(params)) != NULL, XLAL_EFUNC );

   MultiSFTVector *sftvector = NULL;
   if (params->markBadSFTs && !params->signalOnly) {
      //Extract band to get a MultiSFTVector
      MultiSFTVector *tmpsftvector = NULL;
      XLAL_CHECK_NULL( (tmpsftvector = extractSFTband(catalog, minfbin, maxfbin)) != NULL, XLAL_EFUNC );
      XLAL_CHECK_NULL( XLALReorderMultiSFTVector(tmpsftvector, params->IFO) == XLAL_SUCCESS, XLAL_EFUNC );

      //Get the timestamps of the SFTs applying the KS/Kuipers tests if desired
      MultiLIGOTimeGPSVector *multiTimestamps = NULL;
      XLAL_CHECK_NULL( (multiTimestamps = getMultiTimeStampsFromSFTs(tmpsftvector, params)) != NULL, XLAL_EFUNC );

      //Get the SFT subset
      XLAL_CHECK_NULL( (sftvector = XLALExtractMultiSFTVectorWithMultiTimestamps(tmpsftvector, multiTimestamps)) != NULL, XLAL_EFUNC );

      XLALDestroyMultiSFTVector(tmpsftvector);
      XLALDestroyMultiTimestamps(multiTimestamps);
   } else {
      XLAL_CHECK_NULL( (sftvector = extractSFTband(catalog, minfbin, maxfbin)) != NULL, XLAL_EFUNC );
      XLAL_CHECK_NULL( XLALReorderMultiSFTVector(sftvector, params->IFO) == XLAL_SUCCESS, XLAL_EFUNC );
   }

   XLALDestroySFTCatalog(catalog);

   return sftvector;
} // getMultiSFTVector()

/**
 * Add SFTs together from a MultiSFTVector
 * \param [in]  multiSFTvector    Pointer to a MultiSFTVector containing the SFT data
 * \param [in]  multissb          Pointer to a MultiSSBtimes structure
 * \param [in]  multiAMcoeffs     Pointer to a MultiAMCoeffs structure
 * \param [in]  jointTimestamps   Pointer to a LIGOTimeGPSVector of joint SFT times from all detectors
 * \param [in]  backgroundRatio   Pointer to a REAL4VectorAlignedArray of the running means of each IFO divided by the running mean of IFO_0
 * \param [in]  cosiSign          Value of 1, 0, -1 where 1 means average cosi over [0,1], 0 means average over [-1,1], and -1 means average over [-1,0]
 * \param [in]  NSparams          Values for a set of assumed NS spin/binary motion parameters
 * \param [in]  params            Pointer to UserInput_t
 * \param [out] backgroundScaling Pointer to REAL4VectorAligned of background scaling values
 * \return REAL4VectorAligned of powers of the coherently combined SFTs
 */
REAL4VectorAligned * coherentlyAddSFTs(const MultiSFTVector *multiSFTvector, const MultiSSBtimes *multissb, const MultiAMCoeffs *multiAMcoeffs, const LIGOTimeGPSVector *jointTimestamps, const REAL4VectorAlignedArray *backgroundRatio, const INT4 cosiSign, const assumeNSparams *NSparams, const UserInput_t *params, REAL4VectorAligned *backgroundScaling) {

   //Sanity check the values
   XLAL_CHECK_NULL( multiSFTvector != NULL && multissb != NULL && multiAMcoeffs != NULL && jointTimestamps != NULL && params != NULL, XLAL_EINVAL );

   fprintf(stderr, "Coherently adding SFT data... ");

   //parse noise list
   MultiNoiseFloor multiNoiseFloor;
   XLAL_CHECK_NULL( XLALParseMultiNoiseFloor(&multiNoiseFloor, params->avesqrtSh, params->IFO->length) == XLAL_SUCCESS, XLAL_EFUNC );

   //initialize backgroundScaling and normalizationScaling to values of 0
   memset(backgroundScaling->data, 0, sizeof(REAL4)*backgroundScaling->length);

   //get total number of SFTs possible, number of frequency bins in each FFT of the background, and number of original SFT bins to skip because of background calculation
   UINT4 numFbinsInBackground = multiSFTvector->data[0]->data[0].data->length - 20 - (params->blksize - 1);
   UINT4 numSFTbins2skip = (multiSFTvector->data[0]->data[0].data->length - numFbinsInBackground)/2;

   //Allocate the combined SFTs SFTVector
   SFTVector *combinedSFTs = NULL;
   XLAL_CHECK_NULL( (combinedSFTs = XLALCreateSFTVector(jointTimestamps->length, 0)) != NULL, XLAL_EFUNC );

   //Create an INT4Vector to determine at which vector we are using in the multiSFTvector
   INT4Vector *whichSFTinMultiSFTvector = NULL, *whichIFOsToBeUsed = NULL;
   XLAL_CHECK_NULL( (whichSFTinMultiSFTvector = XLALCreateINT4Vector(multiSFTvector->length)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (whichIFOsToBeUsed = XLALCreateINT4Vector(multiSFTvector->length)) != NULL, XLAL_EFUNC );
   memset(whichSFTinMultiSFTvector->data, 0, whichSFTinMultiSFTvector->length * sizeof(INT4));  //All index values are set to zero at the beginning

   //Pre-compute the different sin(2*psi) and cos(2*psi) values
   REAL4VectorAligned *twoPsiVec = NULL, *sin2psiVec = NULL, *cos2psiVec = NULL, *Fplus0s = NULL, *Fcross0s = NULL, *FplusXs = NULL, *FcrossXs = NULL, *aVals = NULL, *bVals = NULL;
   XLAL_CHECK_NULL( (twoPsiVec = XLALCreateREAL4VectorAligned(16, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (sin2psiVec = XLALCreateREAL4VectorAligned(twoPsiVec->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (cos2psiVec = XLALCreateREAL4VectorAligned(twoPsiVec->length, 32)) != NULL, XLAL_EFUNC );
   for (UINT4 ii=0; ii<twoPsiVec->length; ii++) twoPsiVec->data[ii] = 0.03125*ii;
   XLAL_CHECK_NULL( XLALVectorSinCos2PiREAL4(sin2psiVec->data, cos2psiVec->data, twoPsiVec->data, twoPsiVec->length) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK_NULL( (Fplus0s = XLALCreateREAL4VectorAligned(twoPsiVec->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (Fcross0s = XLALCreateREAL4VectorAligned(twoPsiVec->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (FplusXs = XLALCreateREAL4VectorAligned(twoPsiVec->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (FcrossXs = XLALCreateREAL4VectorAligned(twoPsiVec->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (aVals = XLALCreateREAL4VectorAligned(twoPsiVec->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (bVals = XLALCreateREAL4VectorAligned(twoPsiVec->length, 32)) != NULL, XLAL_EFUNC );

   //Pre-compute the GW frequencies and the bin number of the GW frequencies and frequencies of the SFT bins
   alignedREAL8Vector *GWfrequencies = NULL, *GWfrequencyBins = NULL, *SFTfrequencies = NULL, *scaledSFTfrequencies = NULL;
   XLAL_CHECK_NULL( (GWfrequencies = createAlignedREAL8Vector(multiSFTvector->data[0]->data[0].data->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (GWfrequencyBins = createAlignedREAL8Vector(GWfrequencies->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (SFTfrequencies = createAlignedREAL8Vector(GWfrequencies->length, 32)) != NULL, XLAL_EFUNC );
   SFTfrequencies->data[0] = multiSFTvector->data[0]->data[0].f0;
   for (UINT4 ii=1; ii<SFTfrequencies->length; ii++) SFTfrequencies->data[ii] = SFTfrequencies->data[ii-1] + multiSFTvector->data[0]->data[0].deltaF;
   if (NSparams->assumeNSGWfreq == NULL) memcpy(GWfrequencies->data, SFTfrequencies->data, sizeof(REAL8)*SFTfrequencies->length);
   else {
      for (UINT4 ii=0; ii<GWfrequencies->length; ii++) GWfrequencies->data[ii] = *(NSparams->assumeNSGWfreq);
      XLAL_CHECK_NULL( (scaledSFTfrequencies = createAlignedREAL8Vector(SFTfrequencies->length, 32)) != NULL, XLAL_EFUNC );
   }
   XLAL_CHECK_NULL( VectorScaleREAL8(GWfrequencyBins, GWfrequencies, params->Tsft, params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );

   //Pre-allocate delta values for IFO_0 and IFO_X, also 2*pi*f_k*tau
   alignedREAL8Vector *delta0vals = NULL, *deltaXvals = NULL, *delta0valsSubset = NULL, *deltaXvalsSubset = NULL, *floorDelta0valsSubset = NULL, *floorDeltaXvalsSubset = NULL, *diffFloorDeltaValsSubset = NULL, *roundDelta0valsSubset = NULL, *roundDeltaXvalsSubset = NULL, *diffRoundDeltaValsSubset = NULL, *TwoPiGWfrequenciesTau = NULL, *absDiffFloorDeltaValsSubset = NULL;
   REAL4VectorAligned *TwoPiGWfrequenciesTau_f = NULL, *sinTwoPiGWfrequenciesTau = NULL, *cosTwoPiGWfrequenciesTau = NULL;
   XLAL_CHECK_NULL( (delta0vals = createAlignedREAL8Vector(GWfrequencyBins->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (deltaXvals = createAlignedREAL8Vector(GWfrequencyBins->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (delta0valsSubset = createAlignedREAL8Vector(GWfrequencyBins->length - 20, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (deltaXvalsSubset = createAlignedREAL8Vector(delta0valsSubset->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (floorDelta0valsSubset = createAlignedREAL8Vector(delta0valsSubset->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (floorDeltaXvalsSubset = createAlignedREAL8Vector(delta0valsSubset->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (roundDelta0valsSubset = createAlignedREAL8Vector(delta0valsSubset->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (roundDeltaXvalsSubset = createAlignedREAL8Vector(delta0valsSubset->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (diffFloorDeltaValsSubset = createAlignedREAL8Vector(delta0valsSubset->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (absDiffFloorDeltaValsSubset = createAlignedREAL8Vector(delta0valsSubset->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (diffRoundDeltaValsSubset = createAlignedREAL8Vector(delta0valsSubset->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (TwoPiGWfrequenciesTau = createAlignedREAL8Vector(GWfrequencyBins->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (TwoPiGWfrequenciesTau_f = XLALCreateREAL4VectorAligned(GWfrequencyBins->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (sinTwoPiGWfrequenciesTau = XLALCreateREAL4VectorAligned(GWfrequencyBins->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (cosTwoPiGWfrequenciesTau = XLALCreateREAL4VectorAligned(GWfrequencyBins->length, 32)) != NULL, XLAL_EFUNC );

   INT4Vector *shiftVector = NULL;
   XLAL_CHECK_NULL( (shiftVector = XLALCreateINT4Vector(delta0valsSubset->length)) != NULL, XLAL_EFUNC );

   //Mid-point of the SFT
   REAL8 Tmid = 0.5*params->Tsft;

   //Pre-compute 2*pi*{[DeltaTs(jj) - Tmid*Tdots(jj)]-[DeltaTs(0) - Tmid*Tdots(0)]}=twoPiTauVals and [Tdots(jj) - 1]=TdotMinus1s for jj=0,...,N-1
   //Note that twoPiTauVals->data[0] is unused later in the code; it is only necessary for the intermediate calculation
   alignedREAL8VectorArray *twoPiTauVals = NULL, *TdotMinus1s = NULL;
   alignedREAL8Vector *Tdots = NULL, *DeltaTs = NULL;
   XLAL_CHECK_NULL( (twoPiTauVals = createAlignedREAL8VectorArray(multissb->length, multissb->data[0]->DeltaT->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (TdotMinus1s = createAlignedREAL8VectorArray(multissb->length, multissb->data[0]->Tdot->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (Tdots = createAlignedREAL8Vector(multissb->data[0]->Tdot->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (DeltaTs = createAlignedREAL8Vector(multissb->data[0]->DeltaT->length, 32)) != NULL, XLAL_EFUNC );
   for (UINT4 ii=0; ii<twoPiTauVals->length; ii++) {
      memcpy(Tdots->data, multissb->data[ii]->Tdot->data, sizeof(REAL8)*multissb->data[ii]->Tdot->length);
      memcpy(DeltaTs->data, multissb->data[ii]->DeltaT->data, sizeof(REAL8)*multissb->data[ii]->DeltaT->length);
      XLAL_CHECK_NULL( VectorScaleREAL8(twoPiTauVals->data[ii], Tdots, -Tmid, params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_NULL( VectorAddREAL8(twoPiTauVals->data[ii], twoPiTauVals->data[ii], DeltaTs, params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_NULL( VectorShiftREAL8(TdotMinus1s->data[ii], Tdots, -1.0, params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
   }
   for (UINT4 ii=twoPiTauVals->length; ii>0; ii--) XLAL_CHECK_NULL( VectorSubtractREAL8(twoPiTauVals->data[ii-1], twoPiTauVals->data[ii-1], twoPiTauVals->data[0], params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
   for (UINT4 ii=1; ii<twoPiTauVals->length; ii++) XLAL_CHECK_NULL( VectorScaleREAL8(twoPiTauVals->data[ii], twoPiTauVals->data[ii], LAL_TWOPI, params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
   destroyAlignedREAL8Vector(Tdots);
   destroyAlignedREAL8Vector(DeltaTs);

   //Pre-compute cosi and (cosi*cosi+1)/2 values
   REAL4VectorAligned *cosiVals = NULL, *onePlusCosiSqOver2Vals = NULL;
   XLAL_CHECK_NULL( (cosiVals = XLALCreateREAL4VectorAligned(21, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (onePlusCosiSqOver2Vals = XLALCreateREAL4VectorAligned(21, 32)) != NULL, XLAL_EFUNC );
   if (cosiSign==1) for (UINT4 ii=0; ii<cosiVals->length; ii++) cosiVals->data[ii] = 1.0 - 0.05*ii;
   else if (cosiSign==-1) for (UINT4 ii=0; ii<cosiVals->length; ii++) cosiVals->data[ii] = -0.05*ii;
   else for (UINT4 ii=0; ii<cosiVals->length; ii++) cosiVals->data[ii] = 1.0 - 0.05*ii*2.0;
   XLAL_CHECK_NULL( XLALVectorMultiplyREAL4(onePlusCosiSqOver2Vals->data, cosiVals->data, cosiVals->data, cosiVals->length) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK_NULL( XLALVectorShiftREAL4(onePlusCosiSqOver2Vals->data, (REAL4)1.0, onePlusCosiSqOver2Vals->data, onePlusCosiSqOver2Vals->length) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK_NULL( XLALVectorScaleREAL4(onePlusCosiSqOver2Vals->data, (REAL4)0.5, onePlusCosiSqOver2Vals->data, onePlusCosiSqOver2Vals->length) == XLAL_SUCCESS, XLAL_EFUNC );

   //Pre-allocate Aplus and Across
   REAL4VectorAligned *Aplus0s = NULL, *AplusXs = NULL, *Across0s = NULL, *AcrossXs = NULL;
   XLAL_CHECK_NULL( (Aplus0s = XLALCreateREAL4VectorAligned(twoPsiVec->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (Across0s = XLALCreateREAL4VectorAligned(twoPsiVec->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (AplusXs = XLALCreateREAL4VectorAligned(twoPsiVec->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (AcrossXs = XLALCreateREAL4VectorAligned(twoPsiVec->length, 32)) != NULL, XLAL_EFUNC );

   //Pre-allocate DirichletScaling0 and DirichletScaling1
   alignedREAL8Vector *DirichletScaling0 = NULL, *DirichletScalingX = NULL, *scaling = NULL;
   XLAL_CHECK_NULL( (DirichletScaling0 = createAlignedREAL8Vector(delta0valsSubset->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (DirichletScalingX = createAlignedREAL8Vector(delta0valsSubset->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (scaling = createAlignedREAL8Vector(delta0valsSubset->length, 32)) != NULL, XLAL_EFUNC );

   //Pre-allocate the Dirichlet ratio vector and its magnitude
   COMPLEX8Vector *Dratio = NULL;
   REAL4VectorAligned *cabsDratio = NULL;
   XLAL_CHECK_NULL( (Dratio = XLALCreateCOMPLEX8Vector(delta0valsSubset->length)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (cabsDratio = XLALCreateREAL4VectorAligned(Dratio->length, 32)) != NULL, XLAL_EFUNC );

   //Determine new band size
   REAL8 newFmin = SFTfrequencies->data[10];
   REAL8 newBand = round((SFTfrequencies->data[SFTfrequencies->length-11]-SFTfrequencies->data[10])*params->Tsft)/params->Tsft;

   //If using knowledge of a source, then we simulate it using pulsar simulation routines in order to determine the SSB GW frequency
   REAL4TimeSeries *GWfreqInSSB = NULL;
   if (NSparams->assumeNSGWfreq != NULL) {
      PulsarParams XLAL_INIT_DECL(pulsarParams);
      LIGOTimeGPS SSBepoch = multissb->data[0]->refTime, spinEpoch;
      XLALGPSAdd(&SSBepoch, multissb->data[0]->DeltaT->data[0]-params->SFToverlap);
      pulsarParams.Amp.h0 = 0.0;   //arbitrary, we only care about frequency
      pulsarParams.Amp.cosi = 0.0; //arbitrary, we only care about frequency
      pulsarParams.Amp.psi = 0.0;  //arbitrary, we only care about frequency
      pulsarParams.Amp.phi0 = 0.0; //arbitrary, we only care about frequency
      if (NSparams->assumeNSrefTime == NULL) {
	 XLALGPSSetREAL8( &spinEpoch, params->t0 );
	 pulsarParams.Doppler.refTime = spinEpoch;
      } else pulsarParams.Doppler.refTime = *(NSparams->assumeNSrefTime);
      pulsarParams.Doppler.Alpha = NSparams->assumeNSpos.longitude;
      pulsarParams.Doppler.Delta = NSparams->assumeNSpos.latitude;
      pulsarParams.Doppler.fkdot[0] = *(NSparams->assumeNSGWfreq);
      pulsarParams.Doppler.asini = *(NSparams->assumeNSasini);
      pulsarParams.Doppler.period = *(NSparams->assumeNSorbitP);
      pulsarParams.Doppler.ecc = 0.0;
      pulsarParams.Doppler.tp = *(NSparams->assumeNSorbitTp);
      pulsarParams.Doppler.argp = 0.0;
      XLAL_CHECK_NULL( (GWfreqInSSB = computeNSfreqTS(&pulsarParams, SSBepoch, params->Tsft, params->SFToverlap, params->Tobs)) != NULL, XLAL_EFUNC );
   }

   //Loop over the combinedSFTs vector to fill it with single or coherently combined SFTs
   for (UINT4 ii=0; ii<combinedSFTs->length; ii++) {
      //Loop over the interferometers, determining which ones to use and don't go off the end of a list
      memset(whichIFOsToBeUsed->data, 0, whichIFOsToBeUsed->length * sizeof(INT4));  //All index values are set to zero each time
      BOOLEAN foundAnSFT = 0;
      for (UINT4 jj=0; jj<multiSFTvector->length; jj++) {
         if (whichSFTinMultiSFTvector->data[jj]<(INT4)multiSFTvector->data[jj]->length && XLALGPSCmp(&(jointTimestamps->data[ii]), &(multiSFTvector->data[jj]->data[whichSFTinMultiSFTvector->data[jj]].epoch)) == 0) {
            whichIFOsToBeUsed->data[jj] = 1;
            foundAnSFT = 1;
         }
      }
      XLAL_CHECK_NULL( foundAnSFT == 1, XLAL_EFAILED );
      REAL8 sftstart2 = XLALGPSGetREAL8(&(jointTimestamps->data[ii]));
      XLAL_CHECK_NULL( xlalErrno == 0, XLAL_EFUNC );
      INT4 fftnum2 = (INT4)round((sftstart2 - params->t0)/params->SFToverlap);

      //Set backgroundScaling and normalizationScaling to 1 for this SFT, and then we will add to it as more SFTs are coherently summed
      for (UINT4 jj=0; jj<numFbinsInBackground; jj++) backgroundScaling->data[fftnum2*numFbinsInBackground + jj] = 1.0;

      //If there is an assumed NS GW frequency, then determine the frequency as a function of time and scale the SFT frequency bins by -1/assumeNSGWfreq
      if (NSparams->assumeNSGWfreq != NULL) {
	 for (UINT4 jj=0; jj<GWfrequencies->length; jj++) GWfrequencies->data[jj] = GWfreqInSSB->data->data[fftnum2];
	 XLAL_CHECK_NULL( VectorScaleREAL8(GWfrequencyBins, GWfrequencies, params->Tsft, params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
	 XLAL_CHECK_NULL( VectorScaleREAL8(scaledSFTfrequencies, SFTfrequencies, -1.0/GWfrequencies->data[0], params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
      }

      BOOLEAN createSFT = 1, computeAntenna0 = 1, computeDelta0vals = 1;
      for (UINT4 jj=0; jj<multiSFTvector->length; jj++) {
         if (jj==0 && whichIFOsToBeUsed->data[jj]==1) {
            //Copy the data from the multiSFTvector into the combinedSFTs vector, ignoring 10 bins on both ends of the SFT (newFmin and newBand)
            SFTtype *thisSFT = &(combinedSFTs->data[ii]);
            XLAL_CHECK_NULL( XLALExtractBandFromSFT(&thisSFT, &(multiSFTvector->data[0]->data[whichSFTinMultiSFTvector->data[0]]), newFmin, newBand) == XLAL_SUCCESS, XLAL_EFUNC );
            createSFT = 0;
            whichSFTinMultiSFTvector->data[0]++;
         } else if (jj>0 && whichIFOsToBeUsed->data[jj]==1) {
            //Create a copy of the SFT to be shifted since we will manipulate the SFT coefficients
            SFTtype *sftcopySubset = NULL;
            XLAL_CHECK_NULL( XLALExtractBandFromSFT(&sftcopySubset, &(multiSFTvector->data[jj]->data[whichSFTinMultiSFTvector->data[jj]]), newFmin, newBand) == XLAL_SUCCESS, XLAL_EFUNC );
            REAL8 sftstart = XLALGPSGetREAL8(&(sftcopySubset->epoch));
            XLAL_CHECK_NULL( xlalErrno == 0, XLAL_EFUNC );
            INT4 fftnum = (INT4)round((sftstart - params->t0)/params->SFToverlap);
            XLAL_CHECK_NULL( fftnum == fftnum2, XLAL_EFAILED );

            //First do time of arrival: 2*pi*f*tau where f is the GW frequency either assumed (assumeNSGWfreq) or unassumed (f = fk)
            XLAL_CHECK_NULL( VectorScaleREAL8(TwoPiGWfrequenciesTau, GWfrequencies, twoPiTauVals->data[jj]->data[fftnum], params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );

            //If assumeNSpsi is not given compute Fplus and Fcross for detector ratio
            if (NSparams->assumeNSpsi==NULL) {
               if (computeAntenna0) {
                  XLAL_CHECK_NULL( XLALVectorScaleREAL4(aVals->data, multiAMcoeffs->data[0]->a->data[fftnum], cos2psiVec->data, cos2psiVec->length) == XLAL_SUCCESS, XLAL_EFUNC );
                  XLAL_CHECK_NULL( XLALVectorScaleREAL4(bVals->data, multiAMcoeffs->data[0]->b->data[fftnum], sin2psiVec->data, sin2psiVec->length) == XLAL_SUCCESS, XLAL_EFUNC );
                  XLAL_CHECK_NULL( XLALVectorAddREAL4(Fplus0s->data, aVals->data, bVals->data, aVals->length) == XLAL_SUCCESS, XLAL_EFUNC );
                  XLAL_CHECK_NULL( XLALVectorScaleREAL4(bVals->data, multiAMcoeffs->data[0]->b->data[fftnum], cos2psiVec->data, cos2psiVec->length) == XLAL_SUCCESS, XLAL_EFUNC );
                  XLAL_CHECK_NULL( XLALVectorScaleREAL4(aVals->data, -multiAMcoeffs->data[0]->a->data[fftnum], sin2psiVec->data, sin2psiVec->length) == XLAL_SUCCESS, XLAL_EFUNC );
                  XLAL_CHECK_NULL( XLALVectorAddREAL4(Fcross0s->data, bVals->data, aVals->data, bVals->length) == XLAL_SUCCESS, XLAL_EFUNC );
                  computeAntenna0 = 0;
               }
               XLAL_CHECK_NULL( XLALVectorScaleREAL4(aVals->data, multiAMcoeffs->data[jj]->a->data[fftnum], cos2psiVec->data, cos2psiVec->length) == XLAL_SUCCESS, XLAL_EFUNC );
               XLAL_CHECK_NULL( XLALVectorScaleREAL4(bVals->data, multiAMcoeffs->data[jj]->b->data[fftnum], sin2psiVec->data, sin2psiVec->length) == XLAL_SUCCESS, XLAL_EFUNC );
               XLAL_CHECK_NULL( XLALVectorAddREAL4(FplusXs->data, aVals->data, bVals->data, aVals->length) == XLAL_SUCCESS, XLAL_EFUNC );
               XLAL_CHECK_NULL( XLALVectorScaleREAL4(bVals->data, multiAMcoeffs->data[jj]->b->data[fftnum], cos2psiVec->data, cos2psiVec->length) == XLAL_SUCCESS, XLAL_EFUNC );
               XLAL_CHECK_NULL( XLALVectorScaleREAL4(aVals->data, -multiAMcoeffs->data[jj]->a->data[fftnum], sin2psiVec->data, sin2psiVec->length) == XLAL_SUCCESS, XLAL_EFUNC );
               XLAL_CHECK_NULL( XLALVectorAddREAL4(FcrossXs->data, bVals->data, aVals->data, bVals->length) == XLAL_SUCCESS, XLAL_EFUNC );
            }

            //Average of detector-signal phase relation, unless assumed NS orientation and GW polarization angle
            //Need to compute (Aplus^X+i*Across^X)/(Aplus^0+i*Across^0) to determine an average value of magnitude and phase
            //where exponent indicates different IFO <X>.
            REAL4 detPhaseArg = 0.0, detPhaseMag = 0.0;
            if (NSparams->assumeNScosi==NULL && NSparams->assumeNSpsi==NULL) {
               BOOLEAN loopbroken = 0;
               for (UINT4 kk=0; kk<cosiVals->length && !loopbroken; kk++) {
                  XLAL_CHECK_NULL( XLALVectorScaleREAL4(AplusXs->data, onePlusCosiSqOver2Vals->data[kk], FplusXs->data, FplusXs->length) == XLAL_SUCCESS, XLAL_EFUNC );
                  XLAL_CHECK_NULL( XLALVectorScaleREAL4(AcrossXs->data, cosiVals->data[kk], FcrossXs->data, FcrossXs->length) == XLAL_SUCCESS, XLAL_EFUNC );
                  XLAL_CHECK_NULL( XLALVectorScaleREAL4(Aplus0s->data, onePlusCosiSqOver2Vals->data[kk], Fplus0s->data, Fplus0s->length) == XLAL_SUCCESS, XLAL_EFUNC );
                  XLAL_CHECK_NULL( XLALVectorScaleREAL4(Across0s->data, cosiVals->data[kk], Fcross0s->data, Fcross0s->length) == XLAL_SUCCESS, XLAL_EFUNC );
                  for (UINT4 ll=0; ll<cos2psiVec->length && !loopbroken; ll++) {
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
               detPhaseMag /= (REAL4)(cosiVals->length*cos2psiVec->length);
               detPhaseArg /= (REAL4)(cosiVals->length*cos2psiVec->length);
            } else if (NSparams->assumeNScosi!=NULL && NSparams->assumeNSpsi==NULL) {
               BOOLEAN loopbroken = 0;
               REAL4 cosi = *(NSparams->assumeNScosi);
               REAL4 onePlusCosiSqOver2 = 0.5*(1.0 + cosi*cosi);
               XLAL_CHECK_NULL( XLALVectorScaleREAL4(AplusXs->data, onePlusCosiSqOver2, FplusXs->data, FplusXs->length) == XLAL_SUCCESS, XLAL_EFUNC );
               XLAL_CHECK_NULL( XLALVectorScaleREAL4(AcrossXs->data, cosi, FcrossXs->data, FcrossXs->length) == XLAL_SUCCESS, XLAL_EFUNC );
               XLAL_CHECK_NULL( XLALVectorScaleREAL4(Aplus0s->data, onePlusCosiSqOver2, Fplus0s->data, Fplus0s->length) == XLAL_SUCCESS, XLAL_EFUNC );
               XLAL_CHECK_NULL( XLALVectorScaleREAL4(Across0s->data, cosi, Fcross0s->data, Fcross0s->length) == XLAL_SUCCESS, XLAL_EFUNC );
               for (UINT4 kk=0; kk<cos2psiVec->length && !loopbroken; kk++) {
                  COMPLEX16 complexnumerator = crect(AplusXs->data[kk], AcrossXs->data[kk]);
                  COMPLEX16 complexdenominator = crect(Aplus0s->data[kk], Across0s->data[kk]);
                  if (cabs(complexdenominator)>1.0e-6) {
                     COMPLEX16 complexval = complexnumerator/complexdenominator;
                     detPhaseMag += fmin(cabs(complexval), 10.0);  //fmin here because sometimes the magnitude value is being divided by a small value and causing biases
                     detPhaseArg += gsl_sf_angle_restrict_pos(carg(complexval));
                  } else {
                     loopbroken = 1;
                     detPhaseMag = 0.0;
                     detPhaseArg = 0.0;
                  }
               }
               detPhaseMag /= (REAL4)(cos2psiVec->length);
               detPhaseArg /= (REAL4)(cos2psiVec->length);
            } else if (NSparams->assumeNScosi==NULL && NSparams->assumeNSpsi!=NULL) {
               REAL4 psi = *(NSparams->assumeNSpsi);
               REAL4 sin2psi = 0.0, cos2psi = 0.0;
               XLAL_CHECK_NULL( XLALSinCosLUT(&sin2psi, &cos2psi, 2.0*psi) == XLAL_SUCCESS, XLAL_EFUNC );
               if (sin2psi>1.0) sin2psi = 1.0;
               if (sin2psi<-1.0) sin2psi = -1.0;
               if (cos2psi>1.0) cos2psi = 1.0;
               if (cos2psi<-1.0) cos2psi = -1.0;
               REAL4 Fplus0 = multiAMcoeffs->data[0]->a->data[fftnum]*cos2psi + multiAMcoeffs->data[0]->b->data[fftnum]*sin2psi;
               REAL4 Fcross0 = multiAMcoeffs->data[0]->b->data[fftnum]*cos2psi - multiAMcoeffs->data[0]->a->data[fftnum]*sin2psi;
               REAL4 FplusX = multiAMcoeffs->data[jj]->a->data[fftnum]*cos2psi + multiAMcoeffs->data[jj]->b->data[fftnum]*sin2psi;
               REAL4 FcrossX = multiAMcoeffs->data[jj]->b->data[fftnum]*cos2psi - multiAMcoeffs->data[jj]->a->data[fftnum]*sin2psi;
               BOOLEAN loopbroken = 0;
               for (UINT4 kk=0; kk<cosiVals->length && !loopbroken; kk++) {
                  COMPLEX16 complexnumerator = crect(FplusX*onePlusCosiSqOver2Vals->data[kk], FcrossX*cosiVals->data[kk]);
                  COMPLEX16 complexdenominator = crect(Fplus0*onePlusCosiSqOver2Vals->data[kk], Fcross0*cosiVals->data[kk]);
                  if (cabs(complexdenominator)>1.0e-6) {
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
               REAL4 psi = *(NSparams->assumeNSpsi);
               REAL8 sin2psi = sin(2.0*psi), cos2psi = cos(2.0*psi);
               REAL8 Fplus0 = multiAMcoeffs->data[0]->a->data[fftnum]*cos2psi + multiAMcoeffs->data[0]->b->data[fftnum]*sin2psi;
               REAL8 Fcross0 = multiAMcoeffs->data[0]->b->data[fftnum]*cos2psi - multiAMcoeffs->data[0]->a->data[fftnum]*sin2psi;
               REAL8 FplusX = multiAMcoeffs->data[jj]->a->data[fftnum]*cos2psi + multiAMcoeffs->data[jj]->b->data[fftnum]*sin2psi;
               REAL8 FcrossX = multiAMcoeffs->data[jj]->b->data[fftnum]*cos2psi - multiAMcoeffs->data[jj]->a->data[fftnum]*sin2psi;
               REAL4 cosi = *(NSparams->assumeNScosi);
               REAL8 onePlusCosiSqOver2 = 0.5*(1.0 + cosi*cosi);
               COMPLEX16 complexnumerator = crect(FplusX*onePlusCosiSqOver2, FcrossX*cosi);
               COMPLEX16 complexdenominator = crect(Fplus0*onePlusCosiSqOver2, Fcross0*cosi);
               if (cabs(complexdenominator)>1.0e-6) {
                  COMPLEX16 complexval = complexnumerator/complexdenominator;
                  detPhaseMag = (REAL4)fmin(cabs(complexval), 10.0);  //fmin here because sometimes the magnitude value is being divided by a small value and causing biases
                  detPhaseArg = (REAL4)gsl_sf_angle_restrict_pos(carg(complexval));
               } else {
                  detPhaseMag = 0.0;
                  detPhaseArg = 0.0;
               }
            }

            //When reference detector is off, don't add data that is less sensitive than if the reference detector data would be when it is on
            if (detPhaseMag!=0.0 && whichIFOsToBeUsed->data[0]==0 && detPhaseMag<1.0) {
               detPhaseMag = 0.0;
               detPhaseArg = 0.0;
            }

            //Compute delta values for the frequency difference for the Dirichlet kernel.
	    //This is done either with no assumption on GW frequency: fk*Tsft*[Tdot(jj) - 1] where jj=0,...,N-1
	    //Or it is done with assumption of GW frequency: f*Tsft*[Tdot(jj) - fk/f] where jj=0,...,N-1
            //Also do shifting and computing floor(delta) and round(delta)
	    //Compute scaling values delta*[delta^2 - 1], and Dirichlet scaling = scaling0/scalingX
            if (computeDelta0vals) {
	       if (NSparams->assumeNSGWfreq == NULL) XLAL_CHECK_NULL( VectorScaleREAL8(delta0vals, GWfrequencyBins, TdotMinus1s->data[0]->data[fftnum], params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
	       else {
		  XLAL_CHECK_NULL( VectorShiftREAL8(delta0vals, scaledSFTfrequencies, multissb->data[0]->Tdot->data[fftnum], params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
		  XLAL_CHECK_NULL( VectorScaleREAL8(delta0vals, delta0vals, GWfrequencyBins->data[0], params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
	       }
               memcpy(delta0valsSubset->data, &(delta0vals->data[10]), sizeof(REAL8)*delta0valsSubset->length);
               XLAL_CHECK_NULL( VectorFloorREAL8(floorDelta0valsSubset, delta0valsSubset, params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
               XLAL_CHECK_NULL( VectorRoundREAL8(roundDelta0valsSubset, delta0valsSubset, params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
               XLAL_CHECK_NULL( VectorMultiplyREAL8(DirichletScaling0, delta0valsSubset, delta0valsSubset, params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
               XLAL_CHECK_NULL( VectorShiftREAL8(DirichletScaling0, DirichletScaling0, -1.0, params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
               XLAL_CHECK_NULL( VectorMultiplyREAL8(DirichletScaling0, DirichletScaling0, delta0valsSubset, params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
            }
            if (computeDelta0vals) computeDelta0vals = 0;
            if (NSparams->assumeNSGWfreq == NULL) XLAL_CHECK_NULL( VectorScaleREAL8(deltaXvals, GWfrequencyBins, TdotMinus1s->data[jj]->data[fftnum], params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
	    else {
	       XLAL_CHECK_NULL( VectorShiftREAL8(deltaXvals, scaledSFTfrequencies, multissb->data[jj]->Tdot->data[fftnum], params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
	       XLAL_CHECK_NULL( VectorScaleREAL8(deltaXvals, deltaXvals, GWfrequencyBins->data[0], params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
	    }
            memcpy(deltaXvalsSubset->data, &(deltaXvals->data[10]), sizeof(REAL8)*deltaXvalsSubset->length);
            XLAL_CHECK_NULL( VectorRoundREAL8(roundDeltaXvalsSubset, deltaXvalsSubset, params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
            XLAL_CHECK_NULL( VectorSubtractREAL8(diffRoundDeltaValsSubset, roundDelta0valsSubset, roundDeltaXvalsSubset, params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
            XLAL_CHECK_NULL( VectorAddREAL8(deltaXvalsSubset, deltaXvalsSubset, diffRoundDeltaValsSubset, params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
            for (UINT4 kk=0; kk<diffRoundDeltaValsSubset->length; kk++) {
               INT4 shiftVal = (INT4)diffRoundDeltaValsSubset->data[kk];
               shiftVector->data[kk] = shiftVal;
            }
            XLAL_CHECK_NULL( VectorFloorREAL8(floorDeltaXvalsSubset, deltaXvalsSubset, params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
            XLAL_CHECK_NULL( VectorSubtractREAL8(diffFloorDeltaValsSubset, floorDelta0valsSubset, floorDeltaXvalsSubset, params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
            XLAL_CHECK_NULL( VectorMultiplyREAL8(DirichletScalingX, deltaXvalsSubset, deltaXvalsSubset, params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
            XLAL_CHECK_NULL( VectorShiftREAL8(DirichletScalingX, DirichletScalingX, -1.0, params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
            XLAL_CHECK_NULL( VectorMultiplyREAL8(DirichletScalingX, DirichletScalingX, deltaXvalsSubset, params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );

            for (UINT4 kk=0; kk<scaling->length; kk++) scaling->data[kk] = DirichletScaling0->data[kk]/DirichletScalingX->data[kk];
            XLAL_CHECK_NULL( DirichletRatioVector(Dratio, delta0valsSubset, deltaXvalsSubset, scaling, params) == XLAL_SUCCESS, XLAL_EFUNC );

	    //Compute values outside inner loop
	    //If no assumed GW frequency, then normalize Dirichlet kernel ratio values
	    XLAL_CHECK_NULL( VectorCabsfCOMPLEX8(cabsDratio, Dratio, params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
	    if (NSparams->assumeNSGWfreq == NULL) {
	       for (UINT4 kk=0; kk<Dratio->length; kk++) if (cabsDratio->data[kk]!=0.0) Dratio->data[kk] /= cabsDratio->data[kk];
	    } else {
	       for (UINT4 kk=0; kk<Dratio->length; kk++) {
		  if (cabsDratio->data[kk]!=0.0) Dratio->data[kk] /= cabsDratio->data[kk];
	       }
	    }
	    //compute the complex detector-phase relationship
	    COMPLEX8 detPhaseVal = cpolarf(detPhaseMag, detPhaseArg);
	    //If no assumed GW frequency, compute absDiffFloorDeltaValsSubset and phase shift by pi for when signals are on different "sides" of a bin, to reduce error and improve detection efficiency
	    if (NSparams->assumeNSGWfreq == NULL) {
	       XLAL_CHECK_NULL( VectorAbsREAL8(absDiffFloorDeltaValsSubset, diffFloorDeltaValsSubset, params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
	       for (UINT4 kk=0; kk<absDiffFloorDeltaValsSubset->length; kk++) if (absDiffFloorDeltaValsSubset->data[kk]>=1.0) TwoPiGWfrequenciesTau->data[kk+10] -= LAL_PI;
	    }
	    //Cast TwoPiGWfrequenciesTau into a REAL4VectorAligned and take sin and cos
	    //We want to compute exp[-i*(2*pi*f*tau)] = cos[-2*pi*f*tau]+i*sin[-2*pi*f*tau] = cos[2*pi*f*tau]-i*sin[2*pi*f*tau]
	    //Sometimes this is exp[-i*(2*pi*f*tau)+i*pi] = exp[-i*(2*pi*f*tau-pi)] = cos[2*pi*f*tau-pi]-i*sin[2*pi*f*tau-pi]
	    for (UINT4 kk=0; kk<TwoPiGWfrequenciesTau->length; kk++) TwoPiGWfrequenciesTau_f->data[kk] = (REAL4)TwoPiGWfrequenciesTau->data[kk];
	    XLAL_CHECK_NULL( XLALVectorSinCosREAL4(sinTwoPiGWfrequenciesTau->data, cosTwoPiGWfrequenciesTau->data, TwoPiGWfrequenciesTau_f->data, TwoPiGWfrequenciesTau_f->length) == XLAL_SUCCESS, XLAL_EFUNC );

            //Now finish the computation with the Dirichlet kernel ratio and final correction
            for (UINT4 kk=0; kk<sftcopySubset->data->length; kk++) {
               //The complex coefficient to scale SFT bins
	       COMPLEX8 complexfactor = detPhaseVal*crectf(cosTwoPiGWfrequenciesTau->data[kk+10], -sinTwoPiGWfrequenciesTau->data[kk+10])*conj(Dratio->data[kk]);

               REAL4 noiseWeighting = 1.0;
               if (kk>=numSFTbins2skip && kk<numFbinsInBackground+numSFTbins2skip && createSFT) {
                  if (NSparams->assumeNSGWfreq==NULL) backgroundScaling->data[fftnum2*numFbinsInBackground + (kk-numSFTbins2skip)] = detPhaseMag*detPhaseMag;
		  else backgroundScaling->data[fftnum2*numFbinsInBackground + (kk-numSFTbins2skip)] = detPhaseMag*detPhaseMag;
               } else if (kk>=numSFTbins2skip && kk<numFbinsInBackground+numSFTbins2skip && !createSFT) {
                  noiseWeighting = backgroundRatio->data[jj]->data[fftnum];
                  if (NSparams->assumeNSGWfreq==NULL) backgroundScaling->data[fftnum2*numFbinsInBackground + (kk-numSFTbins2skip)] += noiseWeighting*detPhaseMag*detPhaseMag;
		  else backgroundScaling->data[fftnum2*numFbinsInBackground + (kk-numSFTbins2skip)] += noiseWeighting*detPhaseMag*detPhaseMag;
               }

               //Apply the corrections and shift to correct bin with the shiftVector
               sftcopySubset->data->data[kk] = multiSFTvector->data[jj]->data[whichSFTinMultiSFTvector->data[jj]].data->data[kk+10-shiftVector->data[kk]]*noiseWeighting*complexfactor;
            }

            if (createSFT) {
               XLAL_CHECK_NULL( XLALCopySFT(&(combinedSFTs->data[ii]), sftcopySubset) == XLAL_SUCCESS, XLAL_EFUNC );
               createSFT = 0;
            }
            else XLAL_CHECK_NULL( XLALSFTAdd(&(combinedSFTs->data[ii]), sftcopySubset) == XLAL_SUCCESS, XLAL_EFUNC );
            XLALDestroySFT(sftcopySubset);
            whichSFTinMultiSFTvector->data[jj]++;
         }
      } //loop over detectors
   } //loop over SFT times

   XLALDestroyINT4Vector(whichSFTinMultiSFTvector);
   XLALDestroyINT4Vector(whichIFOsToBeUsed);
   XLALDestroyREAL4VectorAligned(cosiVals);
   XLALDestroyREAL4VectorAligned(onePlusCosiSqOver2Vals);
   XLALDestroyREAL4VectorAligned(twoPsiVec);
   XLALDestroyREAL4VectorAligned(sin2psiVec);
   XLALDestroyREAL4VectorAligned(cos2psiVec);
   XLALDestroyREAL4VectorAligned(aVals);
   XLALDestroyREAL4VectorAligned(bVals);
   XLALDestroyREAL4VectorAligned(Fplus0s);
   XLALDestroyREAL4VectorAligned(Fcross0s);
   XLALDestroyREAL4VectorAligned(FplusXs);
   XLALDestroyREAL4VectorAligned(FcrossXs);
   XLALDestroyREAL4VectorAligned(TwoPiGWfrequenciesTau_f);
   XLALDestroyREAL4VectorAligned(sinTwoPiGWfrequenciesTau);
   XLALDestroyREAL4VectorAligned(cosTwoPiGWfrequenciesTau);
   destroyAlignedREAL8Vector(GWfrequencies);
   destroyAlignedREAL8Vector(GWfrequencyBins);
   destroyAlignedREAL8Vector(SFTfrequencies);
   if (NSparams->assumeNSGWfreq != NULL) {
      destroyAlignedREAL8Vector(scaledSFTfrequencies);
      XLALDestroyREAL4TimeSeries(GWfreqInSSB);
   }
   destroyAlignedREAL8VectorArray(twoPiTauVals);
   destroyAlignedREAL8VectorArray(TdotMinus1s);
   destroyAlignedREAL8Vector(delta0vals);
   destroyAlignedREAL8Vector(deltaXvals);
   destroyAlignedREAL8Vector(delta0valsSubset);
   destroyAlignedREAL8Vector(deltaXvalsSubset);
   destroyAlignedREAL8Vector(floorDelta0valsSubset);
   destroyAlignedREAL8Vector(floorDeltaXvalsSubset);
   destroyAlignedREAL8Vector(roundDelta0valsSubset);
   destroyAlignedREAL8Vector(roundDeltaXvalsSubset);
   destroyAlignedREAL8Vector(diffFloorDeltaValsSubset);
   destroyAlignedREAL8Vector(absDiffFloorDeltaValsSubset);
   destroyAlignedREAL8Vector(diffRoundDeltaValsSubset);
   XLALDestroyINT4Vector(shiftVector);
   destroyAlignedREAL8Vector(TwoPiGWfrequenciesTau);
   XLALDestroyREAL4VectorAligned(Aplus0s);
   XLALDestroyREAL4VectorAligned(Across0s);
   XLALDestroyREAL4VectorAligned(AplusXs);
   XLALDestroyREAL4VectorAligned(AcrossXs);
   destroyAlignedREAL8Vector(DirichletScaling0);
   destroyAlignedREAL8Vector(DirichletScalingX);
   destroyAlignedREAL8Vector(scaling);
   XLALDestroyCOMPLEX8Vector(Dratio);
   XLALDestroyREAL4VectorAligned(cabsDratio);

   fprintf(stderr, "done\n");

   REAL4VectorAligned *tfdata = NULL;
   XLAL_CHECK_NULL( (tfdata = convertSFTdataToPowers(combinedSFTs, params, 2.0/params->Tsft/(multiNoiseFloor.sqrtSn[0]*multiNoiseFloor.sqrtSn[0]))) != NULL, XLAL_EFUNC );
   XLALDestroySFTVector(combinedSFTs);

   return tfdata;

} // coherentlyAddSFTs()

/**
 * Compute the GW frequency time series (in the SSB) for a set of assumed parameters. This is supposed to be the same as what is done for Makefakedata
 * \param [in] pulsarParams Pointer to PulsarParams struct
 * \param [in] tepoch       Start time in SSB
 * \param [in] Tsft         SFT length (s)
 * \param [in] SFToverlap   Overlap of the SFTs (s)
 * \param [in] duration     Total duration of the search
 * \return Pointer to REAL4TimeSeries containing the resulting frequency as a function of time
 */
REAL4TimeSeries * computeNSfreqTS(const PulsarParams *pulsarParams, LIGOTimeGPS tepoch, REAL8 Tsft, REAL8 SFToverlap, REAL8 duration)
{
   XLAL_CHECK_NULL( pulsarParams != NULL && duration > 0, XLAL_EINVAL );

   REAL8 h0 = pulsarParams->Amp.h0;
   REAL8 cosi = pulsarParams->Amp.cosi;

   SpinOrbitCWParamStruc XLAL_INIT_DECL(sourceParams);
   sourceParams.psi = pulsarParams->Amp.psi;
   sourceParams.aPlus = 0.5*h0*(1.0+cosi*cosi);
   sourceParams.aCross = h0*cosi;
   sourceParams.phi0 = pulsarParams->Amp.phi0;
   sourceParams.f0 = pulsarParams->Doppler.fkdot[0];
   sourceParams.position.latitude = pulsarParams->Doppler.Delta;
   sourceParams.position.longitude = pulsarParams->Doppler.Alpha;
   sourceParams.position.system = COORDINATESYSTEM_EQUATORIAL;

   sourceParams.orbitEpoch = pulsarParams->Doppler.tp;
   sourceParams.omega = pulsarParams->Doppler.argp;
   sourceParams.oneMinusEcc = 1.0 - pulsarParams->Doppler.ecc;
   sourceParams.rPeriNorm = pulsarParams->Doppler.asini*sourceParams.oneMinusEcc;
   sourceParams.angularSpeed = (LAL_TWOPI/pulsarParams->Doppler.period)*sqrt((1.0 + pulsarParams->Doppler.ecc)/pow(sourceParams.oneMinusEcc, 3.0));

   sourceParams.spinEpoch = pulsarParams->Doppler.refTime;

   sourceParams.deltaT = 5;

   sourceParams.epoch = tepoch;
   sourceParams.length = (UINT4)round(duration/sourceParams.deltaT);

   PulsarCoherentGW XLAL_INIT_DECL(sourceSignal);

   XLAL_CHECK_NULL( XLALGenerateSpinOrbitCW(&sourceSignal, &sourceParams) == XLAL_SUCCESS, XLAL_EFUNC );

   UINT4 downsampledSeriesLength = (UINT4)floor(duration/SFToverlap-1);

   REAL4TimeSeries *out = NULL;
   XLAL_CHECK_NULL( (out = XLALCreateREAL4TimeSeries("Time series", &(sourceParams.spinEpoch), 0, SFToverlap, &lalHertzUnit, downsampledSeriesLength)) != NULL, XLAL_EFUNC );

   REAL4VectorAligned *tempVector = NULL;
   XLAL_CHECK_NULL( (tempVector = XLALCreateREAL4VectorAligned((UINT4)round(Tsft/5), 32)) != NULL, XLAL_EFUNC );
   UINT4 stepsize = tempVector->length/2;
   XLAL_CHECK_NULL( (downsampledSeriesLength-1)*stepsize+tempVector->length <= sourceSignal.f->data->length, XLAL_EBADLEN ); //Check to make sure we don't go past the end of the array

   //compute the mean value of the frequency during an SFT
   for (UINT4 ii=0; ii<downsampledSeriesLength; ii++) {
      memcpy(tempVector->data, &(sourceSignal.f->data->data[ii*stepsize]), sizeof(REAL4)*tempVector->length);
      out->data->data[ii] = calcMean(tempVector);
   }

   //Destroy things
   XLALDestroyREAL4VectorAligned(tempVector);
   XLALDestroyREAL4VectorSequence(sourceSignal.a->data);
   XLALFree( sourceSignal.a );
   XLALDestroyREAL4TimeSeries(sourceSignal.f);
   XLALDestroyREAL8TimeSeries(sourceSignal.phi);

   return out;
}

/**
 * Convert a SFTVector of sfts into powers
 * \param [in] sfts          Pointer to the SFTVector
 * \param [in] params        Pointer to the UserInput_t
 * \param [in] normalization Normalization value to prevent underflow
 * \return Pointer to REAL4VectorAligned containing powers
 */
REAL4VectorAligned * convertSFTdataToPowers(const SFTVector *sfts, const UserInput_t *params, const REAL8 normalization)
{

   XLAL_CHECK_NULL( sfts != NULL && params != NULL, XLAL_EINVAL );

   fprintf(LOG, "Converting band to powers... ");
   fprintf(stderr, "Converting band to powers... ");

   UINT4 numffts = (UINT4)floor(params->Tobs/(params->Tsft-params->SFToverlap)-1);
   UINT4 sftlength = sfts->data[0].data->length;
   XLAL_CHECK_NULL( sftlength != 0, XLAL_EINVAL, "SFT has length of 0!\n" );

   //Allocate output
   REAL4VectorAligned *tfdata = NULL;
   XLAL_CHECK_NULL( (tfdata = XLALCreateREAL4VectorAligned(numffts*sftlength, 32)) != NULL, XLAL_EFUNC );

   //Non-existant SFT counter
   UINT4 nonexistantsft = 0;

   REAL8 starttime = params->t0;

   //Allocate temporary normalized abs[sftcoeff]^2 vector
   alignedREAL8Vector *normAbsSFTcoeffSq = NULL;
   XLAL_CHECK_NULL( (normAbsSFTcoeffSq = createAlignedREAL8Vector(sftlength, 32)) != NULL, XLAL_EFUNC );

   //Load the data into the output vector, roughly normalizing as we go along from the input value
   //REAL8 sqrtnorm = sqrt(normalization);
   for (UINT4 ii=0; ii<numffts; ii++) {
      if (ii-nonexistantsft < sfts->length) {
         SFTtype *sft = &(sfts->data[ii - nonexistantsft]);
         if (sft->epoch.gpsSeconds == (INT4)round(ii*(params->Tsft-params->SFToverlap)+starttime)) {
	    XLAL_CHECK_NULL( VectorCabsCOMPLEX8(normAbsSFTcoeffSq, sft->data, params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
	    XLAL_CHECK_NULL( VectorMultiplyREAL8(normAbsSFTcoeffSq, normAbsSFTcoeffSq, normAbsSFTcoeffSq, params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
	    XLAL_CHECK_NULL( VectorScaleREAL8(normAbsSFTcoeffSq, normAbsSFTcoeffSq, normalization, params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
            for (UINT4 jj=0; jj<sftlength; jj++) {
	       tfdata->data[ii*sftlength + jj] = (REAL4)normAbsSFTcoeffSq->data[jj];
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

   destroyAlignedREAL8Vector(normAbsSFTcoeffSq);

   fprintf(LOG, "done\n");
   fprintf(stderr, "done\n");

   fprintf(LOG, "Duty factor = %f\n", 1.0-(REAL4)nonexistantsft/(REAL4)numffts);
   fprintf(stderr, "Duty factor = %f\n", 1.0-(REAL4)nonexistantsft/(REAL4)numffts);

   //REAL4 meanTFdata = calcMean(tfdata);
   //REAL4 stddev = 0.0;
   //XLAL_CHECK_NULL( calcStddev(&stddev, tfdata) == XLAL_SUCCESS, XLAL_EFUNC );
   //fprintf(LOG, "TF before weighting, mean subtraction: mean = %g, std. dev. = %g\n", meanTFdata, stddev);
   //fprintf(stderr, "TF before weighting, mean subtraction: mean = %g, std. dev. = %g\n", meanTFdata, stddev);

   return tfdata;

} // convertSFTdataToPowers()

/**
 * Read in the data SFTs in one function
 * \param [in] params        Pointer to UserInput_t
 * \param [in] normalization Normalization value determined from expected noise background
 * \param [in] minfbin       Frequency value of the minimum frequency bin
 * \param [in] maxfbin       Frequency value of the maximum frequency bin
 * \return REAL4VectorAligned of SFT powers
 */
REAL4VectorAligned * readInSFTs(UserInput_t *params, const REAL8 normalization, const REAL8 minfbin, const REAL8 maxfbin)
{

   XLAL_CHECK_NULL( params != NULL, XLAL_EINVAL );

   MultiSFTVector *sftvector = NULL;
   XLAL_CHECK_NULL( (sftvector = getMultiSFTVector(params, minfbin, maxfbin)) != NULL, XLAL_EFUNC );

   REAL4VectorAligned *tfdata = NULL;
   XLAL_CHECK_NULL( (tfdata = convertSFTdataToPowers(sftvector->data[0], params, normalization)) != NULL, XLAL_EFUNC );

   //Destroy stuff
   XLALDestroyMultiSFTVector(sftvector);

   return tfdata;

} // readInSFTs()

/**
 * Create a list of timestamps from an SFTCatalog
 * \param [in] catalog Pointer to an SFTCatalog
 * \return Pointer to a list of GPS timestamps in a MultiLIGOTimeGPSVector
 */
MultiLIGOTimeGPSVector * getMultiTimeStampsFromSFTCatalog(const SFTCatalog *catalog)
{

   XLAL_CHECK_NULL( catalog != NULL, XLAL_EINVAL );

   //Get the MultiSFTCatalogView
   MultiSFTCatalogView *catalogView = NULL;
   XLAL_CHECK_NULL( (catalogView = XLALGetMultiSFTCatalogView(catalog)) != NULL, XLAL_EFUNC );

   MultiLIGOTimeGPSVector *multiTimestamps = NULL;
   XLAL_CHECK_NULL( (multiTimestamps = XLALTimestampsFromMultiSFTCatalogView(catalogView)) != NULL, XLAL_EFUNC );

   XLALDestroyMultiSFTCatalogView(catalogView);

   return multiTimestamps;

} // getMultiTimeStampsFromSFTCatalog()

/**
 * Create a list of timestamps from SFTs that might be a subset from those in an SFTCatalog, applying KS/Kuipers test if desired
 * \param [in] multiSFTvector Pointer to a MultiSFTVector
 * \param [in] params         Pointer to UserInput_t
 * \return Pointer to a list of GPS timestamps in a MultiLIGOTimeGPSVector
 */
MultiLIGOTimeGPSVector * getMultiTimeStampsFromSFTs(const MultiSFTVector *multiSFTvector, const UserInput_t *params)
{

   XLAL_CHECK_NULL( params != NULL, XLAL_EINVAL );

   MultiLIGOTimeGPSVector *multiTimestamps = NULL;

   //parse noise list
   MultiNoiseFloor multiNoiseFloor;
   XLAL_CHECK_NULL( XLALParseMultiNoiseFloor(&multiNoiseFloor, params->avesqrtSh, params->IFO->length) == XLAL_SUCCESS, XLAL_EFUNC );

   if (params->markBadSFTs && !params->signalOnly) {
      //REAL8 tfnormval = 2.0/(params->Tsft*(params->avesqrtSh*params->avesqrtSh));

      XLAL_CHECK_NULL( (multiTimestamps = XLALCalloc(1, sizeof(*multiTimestamps))) != NULL, XLAL_ENOMEM );
      XLAL_CHECK_NULL( (multiTimestamps->data = XLALCalloc(multiSFTvector->length, sizeof(*multiTimestamps->data))) != NULL, XLAL_ENOMEM );
      multiTimestamps->length = multiSFTvector->length;

      for (UINT4 ii=0; ii<multiSFTvector->length; ii++) {
         REAL8 tfnormval = 2.0/(params->Tsft*(multiNoiseFloor.sqrtSn[ii]*multiNoiseFloor.sqrtSn[ii]));

         REAL4VectorAligned *tfdata = NULL;
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
         XLALDestroyREAL4VectorAligned(tfdata);
      }
   } else {
      XLAL_CHECK_NULL( (multiTimestamps = XLALExtractMultiTimestampsFromSFTs(multiSFTvector)) != NULL, XLAL_EFUNC );
   }

   return multiTimestamps;

} // getMultiTimeStampsFromSFTs()

/**
 * Create a list of timestamps from a segment list
 * \param [in] filenames  LALStringVector of filenames
 * \param [in] t0         Start time of the search
 * \param [in] Tsft       SFT timespan in seconds
 * \param [in] SFToverlap Overlap of the SFTs in seconds
 * \param [in] dur        Duration of the search in seconds
 * \return Pointer to a list of GPS timestamps in a MultiLIGOTimeGPSVector
 */
MultiLIGOTimeGPSVector * getMultiTimeStampsFromSegmentsFile(const LALStringVector *filenames, const REAL8 t0, const REAL8 Tsft, const REAL8 SFToverlap, const REAL8 dur)
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

} // getMultiTimeStampsFromSegmentsFile

LIGOTimeGPSVector * jointTimestampsFromMultiTimestamps(const MultiLIGOTimeGPSVector *multiTimestamps)
{
   XLAL_CHECK_NULL( multiTimestamps!=NULL, XLAL_EINVAL );

   MultiLIGOTimeGPSVector *squeezedMultiTS = NULL;
   XLAL_CHECK_NULL( (squeezedMultiTS = squeezeMultiTimestamps(multiTimestamps)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( squeezedMultiTS->length>0, XLAL_EFAILED, "Squeezed multiTimestamps to zero length" );

   LIGOTimeGPSVector *output = NULL;
   if (squeezedMultiTS->length<2) {
      XLAL_CHECK_NULL( (output = XLALCreateTimestampVector(squeezedMultiTS->data[0]->length)) != NULL, XLAL_EFUNC );
      memcpy(output->data, squeezedMultiTS->data[0]->data, sizeof(LIGOTimeGPS)*output->length);
      output->deltaT = squeezedMultiTS->data[0]->deltaT;
      XLALDestroyMultiTimestamps(squeezedMultiTS);
      return output;
   }

   //First determine the earliest of the timestamps and the last of the timestamps, adding one timestride to the last of the timestamps
   LIGOTimeGPS earliestTimestamp, lastTimesample, currentTimestamp;
   LIGOTimeGPS smallestGPS = LIGOTIMEGPSZERO, largestGPS = LIGOTIMEGPSZERO;
   INT4 ifoWithSmallestGPS = -1, ifoWithLargestGPS = -1;
   for (UINT4 ii=0; ii<squeezedMultiTS->length; ii++) {
      if (ifoWithSmallestGPS < 0) {
         smallestGPS = squeezedMultiTS->data[ii]->data[0];
         ifoWithSmallestGPS = ii;
      } else if (XLALGPSCmp(&(squeezedMultiTS->data[ifoWithSmallestGPS]->data[0]), &(squeezedMultiTS->data[ii]->data[0]))>0) {
         smallestGPS = squeezedMultiTS->data[ii]->data[0];
         ifoWithSmallestGPS = ii;
      }
      if (ifoWithLargestGPS < 0) {
         largestGPS = squeezedMultiTS->data[ii]->data[squeezedMultiTS->data[ii]->length-1];
         ifoWithLargestGPS = ii;
      } else if (XLALGPSCmp(&(squeezedMultiTS->data[ifoWithLargestGPS]->data[squeezedMultiTS->data[ifoWithLargestGPS]->length-1]), &(squeezedMultiTS->data[ii]->data[squeezedMultiTS->data[ii]->length-1]))<0) {
         largestGPS = squeezedMultiTS->data[ii]->data[squeezedMultiTS->data[ii]->length-1];
         ifoWithLargestGPS = ii;
      }
   }
   earliestTimestamp = smallestGPS;
   lastTimesample = largestGPS;
   XLAL_CHECK_NULL( XLALGPSAdd(&lastTimesample, squeezedMultiTS->data[0]->deltaT) != NULL, XLAL_EFUNC );
   currentTimestamp = earliestTimestamp;

   //Now determine the number of unique timestamps for the length of the output vector
   UINT4Vector *indexInTimestampVector = NULL;
   XLAL_CHECK_NULL( (indexInTimestampVector = XLALCreateUINT4Vector(squeezedMultiTS->length)) != NULL, XLAL_EFUNC );
   memset(indexInTimestampVector->data, 0, sizeof(UINT4)*indexInTimestampVector->length);
   UINT4 totalnumTimestamps = 0;
   while ( XLALGPSCmp(&currentTimestamp, &lastTimesample) < 0 ) {
      BOOLEAN foundTimestamp = 0;
      for (UINT4 ii=0; ii<squeezedMultiTS->length; ii++) {
         if (indexInTimestampVector->data[ii]<squeezedMultiTS->data[ii]->length && XLALGPSCmp(&currentTimestamp, &(squeezedMultiTS->data[ii]->data[indexInTimestampVector->data[ii]]))==0 && !foundTimestamp) {
            foundTimestamp = 1;
            totalnumTimestamps++;
            (indexInTimestampVector->data[ii])++;
         } else if (indexInTimestampVector->data[ii]<squeezedMultiTS->data[ii]->length && XLALGPSCmp(&currentTimestamp, &(squeezedMultiTS->data[ii]->data[indexInTimestampVector->data[ii]]))==0 && foundTimestamp) {
            (indexInTimestampVector->data[ii])++;
         }
      }
      XLAL_CHECK_NULL( XLALGPSAdd(&currentTimestamp, 0.5*squeezedMultiTS->data[0]->deltaT) != NULL, XLAL_EFUNC );
   }

   //Reset index vector
   memset(indexInTimestampVector->data, 0, sizeof(UINT4)*indexInTimestampVector->length);

   //Populate the output vector
   XLAL_CHECK_NULL( (output = XLALCreateTimestampVector(totalnumTimestamps)) != NULL, XLAL_EFUNC );
   output->deltaT = squeezedMultiTS->data[0]->deltaT;
   for (UINT4 ii=0; ii<output->length; ii++) {
      smallestGPS.gpsSeconds = 0;
      smallestGPS.gpsNanoSeconds = 0;
      ifoWithSmallestGPS = -1;
      for (UINT4 jj=0; jj<squeezedMultiTS->length; jj++) {
         if (indexInTimestampVector->data[jj]<squeezedMultiTS->data[jj]->length && ifoWithSmallestGPS < 0) {
            smallestGPS = squeezedMultiTS->data[jj]->data[indexInTimestampVector->data[jj]];
            ifoWithSmallestGPS = jj;
         } else if (indexInTimestampVector->data[jj]<squeezedMultiTS->data[jj]->length && XLALGPSCmp(&(squeezedMultiTS->data[ifoWithSmallestGPS]->data[indexInTimestampVector->data[ifoWithSmallestGPS]]), &(squeezedMultiTS->data[jj]->data[indexInTimestampVector->data[jj]]))>0) {
            smallestGPS = squeezedMultiTS->data[jj]->data[indexInTimestampVector->data[jj]];
            ifoWithSmallestGPS = jj;
         }
      }
      output->data[ii] = smallestGPS;
      for (UINT4 jj=0; jj<squeezedMultiTS->length; jj++) {
         if (indexInTimestampVector->data[jj]<squeezedMultiTS->data[jj]->length && XLALGPSCmp(&(squeezedMultiTS->data[jj]->data[indexInTimestampVector->data[jj]]), &smallestGPS)==0) (indexInTimestampVector->data[jj])++;
      }
   }

   XLALDestroyUINT4Vector(indexInTimestampVector);
   XLALDestroyMultiTimestamps(squeezedMultiTS);

   return output;
}

MultiLIGOTimeGPSVector * squeezeMultiTimestamps(const MultiLIGOTimeGPSVector *multiTS)
{
   XLAL_CHECK_NULL( multiTS != NULL, XLAL_EINVAL );
   UINT4 nonZeroLength = 0;
   for (UINT4 ii=0; ii<multiTS->length; ii++) if (multiTS->data[ii]->length>0) nonZeroLength++;
   XLAL_CHECK_NULL( nonZeroLength != 0, XLAL_EFAILED );
   MultiLIGOTimeGPSVector *ret = NULL;
   XLAL_CHECK_NULL( (ret = XLALCalloc(1, sizeof(*ret))) != NULL, XLAL_ENOMEM );
   XLAL_CHECK_NULL( (ret->data = XLALCalloc(nonZeroLength, sizeof(*ret->data))) != NULL, XLAL_ENOMEM );
   ret->length = nonZeroLength;
   UINT4 whichTS = 0;
   for (UINT4 ii=0; ii<nonZeroLength; ii++) {
      while (multiTS->data[whichTS]->length==0) whichTS++;
      XLAL_CHECK_NULL( (ret->data[ii] = XLALCreateTimestampVector(multiTS->data[whichTS]->length)) != NULL, XLAL_EFUNC );
      memcpy(ret->data[ii]->data, multiTS->data[whichTS]->data, sizeof(LIGOTimeGPS)*ret->data[ii]->length);
      ret->data[ii]->deltaT = multiTS->data[whichTS]->deltaT;
      whichTS++;
   }
   return ret;
}

/**
 * \param [in] multiSFTvector Pointer to MultiSFTVector of the SFTs
 * \param [in] params         Pointer to the user input
 * \return Status value
 */
INT4 printSFTtimestamps2File(const MultiSFTVector *multiSFTvector, const UserInput_t *params)
{
   XLAL_CHECK( multiSFTvector!=NULL && params!=NULL, XLAL_EINVAL );

   MultiLIGOTimeGPSVector *GPStimes = NULL;
   XLAL_CHECK( (GPStimes = getMultiTimeStampsFromSFTs(multiSFTvector, params)) != NULL, XLAL_EFUNC );

   CHARVector *outputfile = NULL;
   XLAL_CHECK( (outputfile = XLALCreateCHARVector(strlen(params->outdirectory)+25)) != NULL, XLAL_EFUNC );

   for (UINT4 ii=0; ii<GPStimes->length; ii++) {
      memset(outputfile->data, 0, sizeof(CHAR)*outputfile->length);
      sprintf(outputfile->data, "%s/%s-%s", params->outdirectory, params->IFO->data[ii], "timestamps.dat");

      FILE *INSFTTIMES = NULL;
      XLAL_CHECK( (INSFTTIMES = fopen(outputfile->data, "w")) != NULL, XLAL_EIO, "Couldn't fopen %s", outputfile->data );

      for (UINT4 jj=0; jj<GPStimes->data[ii]->length; jj++) fprintf(INSFTTIMES, "%d 0\n", GPStimes->data[ii]->data[jj].gpsSeconds);

      fclose(INSFTTIMES);
   }

   XLALDestroyMultiTimestamps(GPStimes);
   XLALDestroyCHARVector(outputfile);

   return XLAL_SUCCESS;
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
 * \param [in] tfdata Pointer to REAL4VectorAligned of SFT powers
 * \param [in] params Pointer to UserInput_t
 * \return Pointer to an INT4Vector with marked SFTs to be removed with 1 and keep SFT with 0
 */
INT4Vector * markBadSFTs(const REAL4VectorAligned *tfdata, const UserInput_t *params)
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
   REAL4VectorAligned *tempvect = NULL;
   XLAL_CHECK_NULL( (tempvect = XLALCreateREAL4VectorAligned(numfbins, 32)) != NULL, XLAL_EFUNC );

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
   XLALDestroyREAL4VectorAligned(tempvect);

   fprintf(stderr, "done\n");
   fprintf(stderr, "Fraction excluded in K-S and Kuiper's tests = %f\n", (REAL4)badsfts/(REAL4)totalsfts);

   return output;

} // markBadSFTs()

/**
 * Remove the marked SFTs as bad by setting values to 0
 * \param [in,out] tfdata  Pointer to REAL4VectorAligned of SFT powers
 * \param [in]     badsfts Poienter to INT4Vector of bad SFTs
 */
void removeBadSFTs(REAL4VectorAligned *tfdata, const INT4Vector *badsfts)
{

   XLAL_CHECK_VOID( tfdata != NULL && badsfts != NULL, XLAL_EINVAL );
   fprintf(stderr, "Removing bad SFTs... ");
   UINT4 numfbins_tfdata = tfdata->length/badsfts->length;
   for (UINT4 ii=0; ii<badsfts->length; ii++) if (badsfts->data[ii]==1) memset(&(tfdata->data[ii*numfbins_tfdata]), 0, sizeof(REAL4)*numfbins_tfdata);
   fprintf(stderr, "done.\n");

} // removeBadSFTs()

MultiSFTVector * generateSFTdata(UserInput_t *uvar, const MultiLALDetector *detectors, const EphemerisData *edat, const INT4 maxbinshift, const gsl_rng *rng)
{

   XLAL_CHECK_NULL( uvar!=NULL && detectors!=NULL && edat!=NULL && rng!=NULL, XLAL_EINVAL );

   MultiSFTVector *multiSFTvector = NULL;

   //parse noise list
   MultiNoiseFloor multiNoiseFloor;
   XLAL_CHECK_NULL( XLALParseMultiNoiseFloor(&multiNoiseFloor, uvar->avesqrtSh, uvar->IFO->length) == XLAL_SUCCESS, XLAL_EFUNC );

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
   CWMFDataParams XLAL_INIT_DECL(DataParams);
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
         for (UINT4 ii=0; ii<detectors->length; ii++) DataParams.multiNoiseFloor.sqrtSn[ii] = multiNoiseFloor.sqrtSn[ii];

         fprintf(stderr, "Generating noise SFTs... ");

         XLAL_CHECK_NULL( XLALCWMakeFakeMultiData(&multiSFTvector, NULL, NULL, &DataParams, edat) == XLAL_SUCCESS, XLAL_EFUNC );

         if (DataParams.fMin != TwoSpectFmin) XLAL_CHECK_NULL( XLALMultiSFTVectorResizeBand(multiSFTvector, TwoSpectFmin, TwoSpectBand) == XLAL_SUCCESS, XLAL_EFUNC );

         fprintf(stderr, "done\n");
      } else {
         SFTCatalog *catalog = NULL;
         XLAL_CHECK_NULL( (catalog = findSFTdata(uvar)) != NULL, XLAL_EFUNC );
         MultiSFTVector *tmpmultiSFTvector = NULL;
         XLAL_CHECK_NULL( (tmpmultiSFTvector = extractSFTband(catalog, minfbin, maxfbin)) != NULL, XLAL_EFUNC );
         XLAL_CHECK_NULL( (multiSFTvector = XLALExtractMultiSFTVectorWithMultiTimestamps(tmpmultiSFTvector, multiTimestamps)) != NULL, XLAL_EFUNC );
         XLALDestroyMultiSFTVector(tmpmultiSFTvector);
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

            alignedREAL8Vector *aveSFTsPower = NULL;
            XLAL_CHECK_NULL( (aveSFTsPower = createAlignedREAL8Vector(multiSFTvector->data[0]->data->data->length, 32)) != NULL, XLAL_EFUNC );
            memset(aveSFTsPower->data, 0, sizeof(REAL8)*aveSFTsPower->length);

	    alignedREAL8Vector *SFTpower = NULL;
	    XLAL_CHECK_NULL( (SFTpower = createAlignedREAL8Vector(aveSFTsPower->length, 32)) != NULL, XLAL_EFUNC );

            for (UINT4 jj=0; jj<oneSignalSFTs->data[0]->length; jj++) {
               SFTtype *sft = &(oneSignalSFTs->data[0]->data[jj]);
	       XLAL_CHECK_NULL( VectorCabsCOMPLEX8(SFTpower, sft->data, uvar->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
	       XLAL_CHECK_NULL( VectorMultiplyREAL8(SFTpower, SFTpower, SFTpower, uvar->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
	       XLAL_CHECK_NULL( VectorScaleREAL8(SFTpower, SFTpower, 2.0/uvar->Tsft, uvar->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
	       //for (UINT4 kk=0; kk<aveSFTsPower->length; kk++) {
               //   REAL8 powerval = 2.0*(creal(sft->data->data[kk])*creal(sft->data->data[kk]) + cimag(sft->data->data[kk])*cimag(sft->data->data[kk]))/uvar->Tsft;
	       //  SFTpower->data[kk] = powerval;
	       //  aveSFTsPower->data[kk] += powerval;
               //}
	       XLAL_CHECK_NULL( VectorAddREAL8(aveSFTsPower, aveSFTsPower, SFTpower, uvar->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
	       //UINT4 max_index = max_index_double((REAL8Vector*)SFTpower);
	       //fprintf(SIGNALOUT, "%d %.6g\n", max_index, SFTpower->data[max_index]);
            }
            for (UINT4 jj=0; jj<aveSFTsPower->length; jj++) fprintf(SIGNALOUT,"%.9g %.9g\n", DataParams.fMin+jj/uvar->Tsft, aveSFTsPower->data[jj]/multiSFTvector->data[0]->length);
            XLALDestroyMultiSFTVector(oneSignalSFTs);
            destroyAlignedREAL8Vector(aveSFTsPower);
	    destroyAlignedREAL8Vector(SFTpower);
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

} // generateSFTdata()

/**
 * Slide the time-frequency data to account for detector motion
 * \param [out] output    Pointer to REAL4VectorAligned of SFT powers that have been corrected
 * \param [in]  params    Pointer to UserInput_t
 * \param [in]  tfdata    Pointer to REAL4VectorAligned of SFT powers
 * \param [in]  binshifts Pointer to INT4Vector of bin shift values
 * \return Status value
 */
INT4 slideTFdata(REAL4VectorAligned *output, const UserInput_t *params, const REAL4VectorAligned *tfdata, const INT4Vector *binshifts)
{

   XLAL_CHECK( output != NULL && params != NULL && tfdata != NULL && binshifts != NULL, XLAL_EINVAL );

   fprintf(stderr, "Sliding TF data... ");

   UINT4 numffts = (UINT4)floor(params->Tobs/(params->Tsft-params->SFToverlap)-1);
   UINT4 numfbins = output->length/numffts;
   UINT4 maxbinshift = (tfdata->length/numffts - numfbins)/2;

   for (UINT4 ii=0; ii<numffts; ii++) {
      XLAL_CHECK( abs(binshifts->data[ii])<(INT4)maxbinshift, XLAL_EFAILED, "SFT slide value %d is greater than maximum value predicted (%d)", binshifts->data[ii], maxbinshift );
      memcpy(&(output->data[ii*numfbins]), &(tfdata->data[ii*(numfbins+2*maxbinshift) + maxbinshift + binshifts->data[ii]]), sizeof(REAL4)*numfbins);
   }

   //fprintf(stderr, "Mean = %g ", calcMean(output));

   fprintf(stderr, "done\n");

   return XLAL_SUCCESS;

} // slideTFdata()

/**
 * Determine the running mean of each SFT
 * \param [out] output   Pointer to REAL4VectorAligned of running mean values of each SFT
 * \param [in]  tfdata   Pointer to REAL4VectorAligned of SFT powers
 * \param [in]  numffts  Number of SFTs in the observation time
 * \param [in]  numfbins Number of frequency bins
 * \param [in]  blksize  Number of bins in the running median
 * \return Status value
 */
INT4 tfRngMeans(REAL4VectorAligned *output, const REAL4VectorAligned *tfdata, const UINT4 numffts, const UINT4 numfbins, const UINT4 blksize)
{

   XLAL_CHECK( output != NULL && tfdata != NULL && numffts > 0  && numfbins > 0 && blksize > 0, XLAL_EINVAL );

   fprintf(LOG, "Assessing SFT background... ");
   fprintf(stderr, "Assessing SFT background... ");

   LALStatus XLAL_INIT_DECL(status);
   REAL8 bias;
   UINT4 totalfbins = numfbins + blksize - 1;

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
   REAL4VectorAligned *inpsd = NULL, *mediansout = NULL;
   XLAL_CHECK( (inpsd = XLALCreateREAL4VectorAligned(totalfbins, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (mediansout = XLALCreateREAL4VectorAligned(numfbins, 32)) != NULL, XLAL_EFUNC );

   //Now do the running median
   for (UINT4 ii=0; ii<numffts; ii++) {
      //If the SFT values were not zero, then compute the running median
      if (tfdata->data[ii*totalfbins]!=0.0) {
         //Determine running median value, convert to mean value
         memcpy(inpsd->data, &(tfdata->data[ii*inpsd->length]), sizeof(REAL4)*inpsd->length);

         //calculate running median
         LALSRunningMedian2(&status, (REAL4Vector*)mediansout, (REAL4Vector*)inpsd, block);
         XLAL_CHECK( status.statusCode == 0, XLAL_EFUNC );

         //Now make the output medians into means by multiplying by 1/bias
         for (UINT4 jj=0; jj<mediansout->length; jj++) output->data[ii*numfbins + jj] = (REAL4)(mediansout->data[jj]*invbias);
      } else {
         //Otherwise, set means to zero
         //for (UINT4 jj=0; jj<mediansout->length; jj++) output->data[ii*numfbins + jj] = 0.0;
         memset(&(output->data[ii*numfbins]), 0, sizeof(REAL4)*mediansout->length);
      }
   } /* for ii < numffts */

   //Destroy stuff
   XLALDestroyREAL4VectorAligned(inpsd);
   XLALDestroyREAL4VectorAligned(mediansout);

   fprintf(LOG, "done\n");
   fprintf(stderr, "done\n");
   fprintf(stderr,"Mean of running means = %g\n", calcMean(output));

   return XLAL_SUCCESS;

} // tfRngMeans()

/**
 * \param [in,out] tfdataarray Pointer to a REAL4VectorAlignedArray that has the time-frequency data of powers
 * \param [in]     numffts     Number of FFTs in the total observation time
 * \return Status value
 */
INT4 replaceTFdataWithSubsequentTFdata(REAL4VectorAlignedArray *tfdataarray, const UINT4 numffts)
{
   XLAL_CHECK( tfdataarray!=NULL && numffts>0, XLAL_EINVAL );
   if (tfdataarray->length==1) return XLAL_SUCCESS;
   UINT4 sftlength = tfdataarray->data[0]->length/numffts;
   for (UINT4 ii=0; ii<numffts; ii++) {
      if (tfdataarray->data[0]->data[ii*sftlength] == 0.0) {
         for (UINT4 jj=1; jj<tfdataarray->length; jj++) {
            if (tfdataarray->data[jj]->data[ii*sftlength] != 0.0) {
               memcpy(&(tfdataarray->data[0]->data[ii*sftlength]), &(tfdataarray->data[jj]->data[ii*sftlength]), sizeof(REAL4)*sftlength);
               break;
            }
         }
      }
   }
   return XLAL_SUCCESS;
}

/**
 * Subtract running mean values from the SFT data, modifying input time-frequency data
 * \param [in,out] tfdata            Pointer to REAL4VectorAligned time-frequency data (modified by this function!)
 * \param [in]     rngMeans          Pointer to REAL4VectorAligned of running mean values
 * \param [in]     backgroundScaling Pointer to REAL4VectorAligned of background scaling values
 * \param [in]     numffts           Number of SFTs from observation time
 * \param [in]     numfbins          Number of SFT frequency bins
 * \return Status value
 */
INT4 tfMeanSubtract(REAL4VectorAligned *tfdata, const REAL4VectorAligned *rngMeans, const REAL4VectorAligned *backgroundScaling, const UINT4 numffts, const UINT4 numfbins)
{

   XLAL_CHECK( tfdata != NULL && rngMeans != NULL && backgroundScaling != NULL && numffts > 0 && numfbins > 0, XLAL_EINVAL );
   fprintf(stderr, "Subtracting expected background... ");
   for (UINT4 ii=0; ii<numffts; ii++) if (rngMeans->data[ii*numfbins]!=0.0) for (UINT4 jj=0; jj<numfbins; jj++) tfdata->data[ii*numfbins+jj] -= rngMeans->data[ii*numfbins+jj]*backgroundScaling->data[ii*numfbins+jj];
   //fprintf(stderr, "TF mean after subtraction = %g ", calcMean(tfdata));
   fprintf(stderr, "done\n");
   return XLAL_SUCCESS;

} // tfMeanSubtract()

/**
 * Weight the SFTs based on antenna pattern and noise variance (Equation 11, assuming the input time-frequency data is already mean subtracted)
 * \param [out] output                    Pointer to REAL4VectorAligned of mean subtracted, noise and antenna pattern weighted SFTs
 * \param [in]  tfdata                    Pointer to REAL4VectorAligned of mean subtracted SFTs
 * \param [in]  rngMeans                  Pointer to REAL4VectorAligned of running mean values
 * \param [in]  antPatternWeights         Pointer to REAL4VectorAligned of antenna pattern weights
 * \param [in]  backgroundScaling         Pointer to REAL4VectorAligned of background scaling values
 * \param [in]  indexValuesOfExistingSFTs Pointer to INT4Vector of the index values of the existing SFTs
 * \param [in]  params                    Pointer to UserInput_t
 * \return Status value
 */
INT4 tfWeight(REAL4VectorAligned *output, const REAL4VectorAligned *tfdata, REAL4VectorAligned *rngMeans, REAL4VectorAligned *antPatternWeights, const REAL4VectorAligned *backgroundScaling, const INT4Vector *indexValuesOfExistingSFTs, const UserInput_t *params)
{

   XLAL_CHECK( output!=NULL && tfdata!=NULL && rngMeans!=NULL && antPatternWeights!=NULL && backgroundScaling!=NULL && indexValuesOfExistingSFTs != NULL && params != NULL, XLAL_EINVAL );

   fprintf(stderr, "Applying weighting to SFT data... ");

   UINT4 numffts = antPatternWeights->length;
   UINT4 numfbins = tfdata->length/numffts;

   //Initially set output to zero
   memset(output->data, 0, sizeof(REAL4)*output->length);

   REAL4VectorAligned *antweightssq0 = NULL, *antweightssq = NULL, *rngMeanssq = NULL, *backgroundScalingSq = NULL;
   XLAL_CHECK( (antweightssq0 = XLALCreateREAL4VectorAligned(numffts, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (antweightssq = XLALCreateREAL4VectorAligned(numffts, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (rngMeanssq = XLALCreateREAL4VectorAligned(numffts, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (backgroundScalingSq = XLALCreateREAL4VectorAligned(numffts, 32)) != NULL, XLAL_EFUNC );

   //If user specified that SFTs contain signal only, then the values of the backgrnd vector will be zeros.
   //We must set them equal to 1.0 here, then return them to zeros at the end
   if (params->signalOnly) {
      for (UINT4 ii=0; ii<indexValuesOfExistingSFTs->length; ii++) {
         for (UINT4 jj=0; jj<numfbins; jj++) rngMeans->data[numfbins*indexValuesOfExistingSFTs->data[ii] + jj] = 1.0;
      }
   }

   //Determine antenna pattern weights squared
   XLAL_CHECK( XLALVectorMultiplyREAL4(antweightssq0->data, antPatternWeights->data, antPatternWeights->data, antPatternWeights->length) == XLAL_SUCCESS, XLAL_EFUNC );


   //Loop through the SFT frequency bins and weight the data and normalize
   for (UINT4 ii=0; ii<numfbins; ii++) {

      //scale antenna pattern weights by the right SFT frequency bin value in backgroundScaling**2
      fastSSVectorMultiply_with_stride_and_offset(backgroundScalingSq, backgroundScaling, backgroundScaling, numfbins, numfbins, ii, ii);
      memcpy(antweightssq->data, antweightssq0->data, sizeof(REAL4)*antweightssq0->length);
      XLAL_CHECK( XLALVectorMultiplyREAL4(antweightssq->data, antweightssq->data, backgroundScalingSq->data, antweightssq->length) == XLAL_SUCCESS, XLAL_EFUNC );

      //get the background squared
      fastSSVectorMultiply_with_stride_and_offset(rngMeanssq, rngMeans, rngMeans, numfbins, numfbins, ii, ii);

      //If noiseWeightOff is given, then set all the noise weights to be 1.0
      if (params->noiseWeightOff) for (UINT4 jj=0; jj<rngMeanssq->length; jj++) if (rngMeanssq->data[jj]!=0.0) rngMeanssq->data[jj] = 1.0;

      //Get sum of antenna pattern weight/variances for each frequency bin as a function of time (only for existant SFTs)
      REAL8 sumofweights = determineSumOfWeights(antweightssq, rngMeanssq);
      REAL8 invsumofweights = 1.0/sumofweights;

      //Now do noise weighting, antenna pattern weighting
      for (UINT4 jj=0; jj<indexValuesOfExistingSFTs->length; jj++) {
         output->data[indexValuesOfExistingSFTs->data[jj]*numfbins+ii] = (REAL4)(invsumofweights*antPatternWeights->data[indexValuesOfExistingSFTs->data[jj]]*tfdata->data[indexValuesOfExistingSFTs->data[jj]*numfbins+ii]/rngMeanssq->data[indexValuesOfExistingSFTs->data[jj]]);
      } /* for jj < indexValuesOfExisitingSFTs->length */
   } /* for ii < numfbins */

   //Remember to reset the backgrnd vector to zero
   if (params->signalOnly) memset(rngMeans->data, 0, sizeof(REAL4)*rngMeans->length);

   //Destroy stuff
   XLALDestroyREAL4VectorAligned(antweightssq0);
   XLALDestroyREAL4VectorAligned(antweightssq);
   XLALDestroyREAL4VectorAligned(rngMeanssq);
   XLALDestroyREAL4VectorAligned(backgroundScalingSq);

   fprintf(stderr,"done\n");

   return XLAL_SUCCESS;

} // tfWeight()

/**
 * Determine the sum of the weights
 * \param [in] antweightssq Antenna pattern weights squared
 * \param [in] rngMeanssq   Running mean values squared
 * \return Sum of the weights
 */
REAL8 determineSumOfWeights(const REAL4VectorAligned *antweightssq, const REAL4VectorAligned *rngMeanssq)
{

   XLAL_CHECK_REAL8( antweightssq != NULL && rngMeanssq != NULL, XLAL_EINVAL );

   REAL8 sumofweights = 0.0;
   for (UINT4 ii=0; ii<antweightssq->length; ii++) if (rngMeanssq->data[ii] != 0.0) sumofweights += antweightssq->data[ii]/rngMeanssq->data[ii];

   return sumofweights;

} /* determineSumOfWeights */


/**
 * Determine if the SFTs are existing based on whether the first data point of power in each SFT is 0
 * \param [in] tfdata   Pointer to REAL4VectorAligned of time-frequency data
 * \param [in] numffts  Number of SFTs from the observation time
 * \return Pointer to INT4Vector containing 0 for non-present SFT or 1 for present SFT
 */
INT4Vector * existingSFTs(const REAL4VectorAligned *tfdata, const UINT4 numffts)
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

/**
 * Go through the backgroundScaling vector and zero out if the SFTexistVector has a 0 in that SFT location instead of 1
 * \param [in,out] background        Pointer to REAL4VectorAligned of background data
 * \param [in,out] backgroundScaling Pointer to REAL4VectorAligned of background scaling data
 * \param [in]     SFTexistVector    Pointer to INT4Vector of 0 or 1 values representing if an SFT is present (1) or not (0)
 * \return Status value
 */
INT4 checkBackgroundScaling(const REAL4VectorAligned *background, const REAL4VectorAligned *backgroundScaling, const INT4Vector *SFTexistVector)
{
   XLAL_CHECK( backgroundScaling != NULL && SFTexistVector != NULL, XLAL_EINVAL );

   UINT4 numfbins = backgroundScaling->length/SFTexistVector->length;
   for (UINT4 ii=0; ii<SFTexistVector->length; ii++) {
      if (SFTexistVector->data[ii]==0 && backgroundScaling->data[ii*numfbins]!=0.0) {
         memset(&(backgroundScaling->data[ii*numfbins]), 0, sizeof(REAL4)*numfbins);
      }
      if (SFTexistVector->data[ii]==0 && background->data[ii*numfbins]!=0.0) {
         memset(&(background->data[ii*numfbins]), 0, sizeof(REAL4)*numfbins);
      }
   }
   return XLAL_SUCCESS;
}
