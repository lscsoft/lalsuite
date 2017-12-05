/*
*  Copyright (C) 2010--2012, 2014, 2015 Evan Goetz
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

#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "IHS.h"
#include "statistics.h"
#include "candidates.h"
#include "cdfdist.h"
#include "vectormath.h"
#include "TwoSpect.h"


/**
 * Create vectors for IHS maxima struct
 * \param [in] fbins Number of frequency bins
 * \param [in] rows  Number of neighboring rows to be summed
 * \return Pointer to a newly created ihsMaximaStruct
 */
ihsMaximaStruct * createihsMaxima(const UINT4 fbins, const UINT4 rows)
{

   XLAL_CHECK_NULL( fbins > 0 && rows > 0, XLAL_EINVAL );

   ihsMaximaStruct *ihsmaxima = NULL;
   XLAL_CHECK_NULL( (ihsmaxima = XLALMalloc(sizeof(*ihsmaxima))) != NULL, XLAL_ENOMEM );

   //The number of ihs maxima = Sum[fbins - (i-1), {i, 2, rows}]
   // = fbins*rows - fbins - (rows**2 - rows)/2
   UINT4 numberofmaxima = fbins*rows - fbins - (UINT4)((rows*rows-rows)/2);

   XLAL_CHECK_NULL( (ihsmaxima->maxima = XLALCreateREAL4VectorAligned(numberofmaxima, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (ihsmaxima->locations = XLALCreateINT4Vector(numberofmaxima)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (ihsmaxima->foms = XLALCreateREAL4VectorAligned(numberofmaxima, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (ihsmaxima->maximaForEachFbin = XLALCreateREAL4VectorAligned(fbins, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (ihsmaxima->locationsForEachFbin = XLALCreateINT4Vector(fbins)) != NULL, XLAL_EFUNC );
   ihsmaxima->rows = rows;

   return ihsmaxima;

} /* createihsMaxima() */


/**
 * Destroy vectors and the IHS maxima struct
 * \param [in] data Pointer to an ihsMaximaStruct to be freed
 */
void destroyihsMaxima(ihsMaximaStruct *data)
{
   if (data) {
      XLALDestroyREAL4VectorAligned(data->maxima);
      XLALDestroyINT4Vector(data->locations);
      XLALDestroyREAL4VectorAligned(data->foms);
      XLALDestroyREAL4VectorAligned(data->maximaForEachFbin);
      XLALDestroyINT4Vector(data->locationsForEachFbin);
      XLALFree((ihsMaximaStruct*)data);
   }
} /* destroyihsMaxima() */


/**
 * Run the IHS algorithm
 * \param [out] output      Pointer to the ihsMaximaStruct
 * \param [in]  input       Pointer to the ffdataStruct
 * \param [in]  ihsfarinput Pointer to the ihsfarStruct
 * \param [in]  params      Pointer to UserInput_t
 * \param [in]  rows        Number of neighboring rows to be summed
 * \param [in]  aveNoise    Pointer to a REAL4VectorAligned of 2nd FFT background powers
 * \param [in]  FbinMean    Pointer to a REAL4VectorAligned of normalized SFT background powers
 * \return Status value
 */
INT4 runIHS(ihsMaximaStruct *output, const ffdataStruct *input, const ihsfarStruct *ihsfarinput, const UserInput_t *params, const UINT4 rows, const REAL4VectorAligned *aveNoise, const REAL4VectorAligned *FbinMean)
{

   XLAL_CHECK( output != NULL && input != NULL && ihsfarinput != NULL && params != NULL && rows > 0 && aveNoise != NULL && FbinMean != NULL, XLAL_EINVAL );

   INT4 numfbins = input->numfbins;
   INT4 numfprbins = input->numfprbins;

   //Allocate memory for the necessary vectors
   REAL4VectorAligned *row = NULL, *ihss = NULL, *ihsvector = NULL;
   INT4Vector *locs = NULL;
   ihsVals *ihsvals = NULL;
   REAL4VectorAlignedArray *ihsvectorarray = NULL;
   XLAL_CHECK( (row = XLALCreateREAL4VectorAligned(numfprbins, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (ihss = XLALCreateREAL4VectorAligned(numfbins, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (ihsvector = XLALCreateREAL4VectorAligned((INT4)floor((1.0/(REAL8)params->ihsfactor)*numfprbins)-5, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (locs = XLALCreateINT4Vector(numfbins)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (ihsvals = createihsVals()) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (ihsvectorarray = createREAL4VectorAlignedArray(numfbins, ihsvector->length, 32)) != NULL, XLAL_EFUNC );

   //We want to ignore daily and sidereal harmonics, so mark the values
   REAL8 dailyharmonic = params->Tobs/(24.0*3600.0);
   REAL8 siderealharmonic = params->Tobs/86164.0905;
   REAL8 dailyharmonic2 = dailyharmonic*2.0, dailyharmonic3 = dailyharmonic*3.0, dailyharmonic4 = dailyharmonic*4.0;
   REAL8 siderealharmonic2 = siderealharmonic*2.0, siderealharmonic3 = siderealharmonic*3.0, siderealharmonic4 = siderealharmonic*4.0;
   INT4Vector *markedharmonics = NULL;
   XLAL_CHECK( (markedharmonics = XLALCreateINT4Vector(row->length)) != NULL, XLAL_EFUNC );
   memset(markedharmonics->data, 0, sizeof(INT4)*markedharmonics->length);
   //If the user has specified not to notch the harmonics, then we skip the step to mark the notched values
   if (!params->noNotchHarmonics) {
      for (UINT4 ii=0; ii<markedharmonics->length; ii++) {
         if (fabs(dailyharmonic-(REAL8)ii)<=1.0 || fabs(dailyharmonic2-(REAL8)ii)<=1.0 || fabs(dailyharmonic3-(REAL8)ii)<=1.0 || fabs(dailyharmonic4-(REAL8)ii)<=1.0 || fabs(siderealharmonic-(REAL8)ii)<=1.0 || fabs(siderealharmonic2-(REAL8)ii)<=1.0 || fabs(siderealharmonic3-(REAL8)ii)<=1.0 || fabs(siderealharmonic4-(REAL8)ii)<=1.0) markedharmonics->data[ii] = 1;
      }
   }

   //Loop through the rows, 1 frequency at a time
   for (UINT4 ii=0; ii<ihss->length; ii++) {

      //For each row, populate it with the data for that frequency bin, excluding harmonics of antenna pattern modulation
      memcpy(row->data, &(input->ffdata->data[ii*numfprbins]), sizeof(REAL4)*numfprbins);
      if (!params->noNotchHarmonics) for (UINT4 jj=0; jj<row->length; jj++) if (markedharmonics->data[jj]==1) row->data[jj] = 0.0;

      //Run the IHS algorithm on the row
      if (!params->weightedIHS) XLAL_CHECK( incHarmSumVector(ihsvector, row, params->ihsfactor) == XLAL_SUCCESS, XLAL_EFUNC );
      else XLAL_CHECK( incHarmSumVectorWeighted(ihsvector, row, aveNoise, params->ihsfactor) == XLAL_SUCCESS, XLAL_EFUNC );

      //Copy the result into the ihsvector array
      memcpy(ihsvectorarray->data[ii]->data, ihsvector->data, sizeof(REAL4)*ihsvector->length);

   } /* for ii < ihss->length */

   //Now do the summing of the IHS values
   XLAL_CHECK( sumIHSarray(output, ihsfarinput, ihsvectorarray, rows, FbinMean, params) == XLAL_SUCCESS, XLAL_EFUNC );

   //Destroy stuff
   destroyREAL4VectorAlignedArray(ihsvectorarray);
   XLALDestroyREAL4VectorAligned(row);
   XLALDestroyREAL4VectorAligned(ihss);
   XLALDestroyREAL4VectorAligned(ihsvector);
   XLALDestroyINT4Vector(locs);
   XLALDestroyINT4Vector(markedharmonics);
   destroyihsVals(ihsvals);

   return XLAL_SUCCESS;

} /* runIHS() */


/**
 * Allocate memory for ihsVals struct
 * \return Pointer to a newly created ihsVals structure
 */
ihsVals * createihsVals(void)
{
   ihsVals *ihsvals = NULL;
   XLAL_CHECK_NULL( (ihsvals = XLALMalloc(sizeof(*ihsvals))) != NULL, XLAL_ENOMEM );
   return ihsvals;
} /* createihsVals() */


/**
 * Destroy ihsVals struct
 * \param [in] ihsvals Pointer to an ihsVals structure to be freed
 */
void destroyihsVals(ihsVals *ihsvals)
{
   if (ihsvals) XLALFree((ihsVals*)ihsvals);
} /* destroyihsVals() */


/**
 * Compute the IHS sum maximum
 * \param [out] output    Pointer to the ihsVals structure
 * \param [in]  input     Pointer to a REAL4VectorAligned
 * \param [in]  ihsfactor Number of folds of the 2nd FFT
 * \return Status value
 */
INT4 incHarmSum(ihsVals *output, const REAL4VectorAligned *input, const UINT4 ihsfactor)
{

   XLAL_CHECK( output != NULL && input != NULL && ihsfactor > 0, XLAL_EINVAL );

   output->ihs = 0.0;

   REAL4VectorAligned *tempvect = NULL;
   XLAL_CHECK( (tempvect = XLALCreateREAL4VectorAligned((INT4)floor((1.0/(REAL8)ihsfactor)*input->length)-5, 32)) != NULL, XLAL_EFUNC );

   XLAL_CHECK( incHarmSumVector(tempvect, input, ihsfactor) == XLAL_SUCCESS, XLAL_EFUNC );

   for (UINT4 ii=0; ii<tempvect->length; ii++) {
      if (tempvect->data[ii]>output->ihs) {
         output->ihs = tempvect->data[ii];
         output->loc = ii+5;
      }
   } /* for ii < tempvect->length */

   XLALDestroyREAL4VectorAligned(tempvect);

   return XLAL_SUCCESS;

} /* incHarmSum() */


/**
 * Compute the IHS vector -- does not compute the maximum value
 * \param [out] output    Pointer to a REAL4VectorAligned to contain the folded values
 * \param [in]  input     Pointer to a REAL4VectorAligned for the input to the IHS
 * \param [in]  ihsfactor Number of folds of the 2nd FFT
 * \return Status value
 */
INT4 incHarmSumVector(REAL4VectorAligned *output, const REAL4VectorAligned *input, const UINT4 ihsfactor)
{

   XLAL_CHECK( output != NULL && input != NULL && ihsfactor > 0, XLAL_EINVAL );

   INT4 highval = (INT4)floor((1.0/(REAL8)ihsfactor)*input->length);

   //From 5 up to the highval
   for (INT4 ii=5; ii<highval; ii++) {
      output->data[ii-5] = 0.0;  //first set the value to zero
      for (UINT4 jj=1; jj<=ihsfactor; jj++) output->data[ii-5] += input->data[ii*jj];

      //check that we're not going outside of the allowed limits
      XLAL_CHECK( ihsfactor*(ii-1)<=input->length-1, XLAL_EBADLEN );
   } /* for ii=5 --> highval */

   return XLAL_SUCCESS;

} /* incHarmSumVector() */


/**
 * Compute the noise weighted IHS vector -- does not compute the maximum value
 * \param [out] output    Pointer to a REAL4VectorAligned to contain the folded values
 * \param [in]  input     Pointer to a REAL4VectorAligned for the input to the IHS
 * \param [in]  aveNoise  Pointer to a REAL4VectorAligned of 2nd FFT background powers
 * \param [in]  ihsfactor Number of folds of the 2nd FFT
 * \return Status value
 */
INT4 incHarmSumVectorWeighted(REAL4VectorAligned *output, const REAL4VectorAligned *input, const REAL4VectorAligned *aveNoise, const UINT4 ihsfactor)
{

   XLAL_CHECK( output != NULL && input != NULL && aveNoise != NULL && ihsfactor > 0, XLAL_EINVAL );

   INT4 highval = (INT4)floor((1.0/(REAL8)ihsfactor)*input->length);

   //From 5 up to the highval
   for (INT4 ii=5; ii<highval; ii++) {
      output->data[ii-5] = 0.0;  //first set the value to zero
      for (UINT4 jj=1; jj<=ihsfactor; jj++) {
         REAL4 weight = 1.0/(aveNoise->data[ii*jj]*aveNoise->data[ii*jj]);
         output->data[ii-5] += input->data[ii*jj]*weight;
      }

      //check that we're not going outside of the allowed limits
      XLAL_CHECK( ihsfactor*(ii-1)<=input->length-1, XLAL_EBADLEN );
   } /* for ii=5 --> highval */

   return XLAL_SUCCESS;

} // incHarmSumVectorWeighted()


/**
 * Allocate memory for ihsfarStruct struct
 * \param [in] rows   Number of neighbors to sum
 * \param [in] params Pointer to UserInput_t
 * \return Pointer to a newly allocated ihsfarStruct
 */
ihsfarStruct * createihsfarStruct(const UINT4 rows, const UserInput_t *params)
{

   XLAL_CHECK_NULL( rows > 0 && params != NULL, XLAL_EINVAL );

   ihsfarStruct *ihsfarstruct = NULL;
   XLAL_CHECK_NULL( (ihsfarstruct = XLALMalloc(sizeof(*ihsfarstruct))) != NULL, XLAL_ENOMEM );

   XLAL_CHECK_NULL( (ihsfarstruct->ihsfar = XLALCreateREAL4VectorAligned(rows-1, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (ihsfarstruct->ihsdistMean = XLALCreateREAL4VectorAligned(rows-1, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (ihsfarstruct->ihsdistSigma = XLALCreateREAL4VectorAligned(rows-1, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (ihsfarstruct->fomfarthresh = XLALCreateREAL4VectorAligned(rows-1, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (ihsfarstruct->ihsfomdistMean = XLALCreateREAL4VectorAligned(rows-1, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (ihsfarstruct->ihsfomdistSigma = XLALCreateREAL4VectorAligned(rows-1, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (ihsfarstruct->expectedIHSVector = XLALCreateREAL4VectorAligned((INT4)floor((1.0/(REAL8)params->ihsfactor)*((INT4)floor(floor(params->Tobs/(params->Tsft-params->SFToverlap)-1)*0.5)+1))-5, 32)) != NULL, XLAL_EFUNC );

   memset(ihsfarstruct->expectedIHSVector->data, 0, sizeof(REAL4)*ihsfarstruct->expectedIHSVector->length);

   return ihsfarstruct;

} // createihsfarStruct()


/**
 * Destroy ihsfarStruct struct
 * \param [in] ihsfarstruct Pointer to the ihsfarStruct to be destroyed
 */
void destroyihsfarStruct(ihsfarStruct *ihsfarstruct)
{
   if (ihsfarstruct) {
      XLALDestroyREAL4VectorAligned(ihsfarstruct->ihsfar);
      XLALDestroyREAL4VectorAligned(ihsfarstruct->ihsdistMean);
      XLALDestroyREAL4VectorAligned(ihsfarstruct->ihsdistSigma);
      XLALDestroyREAL4VectorAligned(ihsfarstruct->fomfarthresh);
      XLALDestroyREAL4VectorAligned(ihsfarstruct->ihsfomdistMean);
      XLALDestroyREAL4VectorAligned(ihsfarstruct->ihsfomdistSigma);
      XLALDestroyREAL4VectorAligned(ihsfarstruct->expectedIHSVector);
      XLALFree((ihsfarStruct*)ihsfarstruct);
   }
} /* destroyihsfarStruct() */


/**
 * Compute the IHS FAR for a sum of a number of rows
 * \param [out] output   Pointer to the output ihsfarStruct
 * \param [in]  params   Pointer to UserInput_t
 * \param [in]  rows     Number of neighbors to sum
 * \param [in]  aveNoise Pointer to REAL4VectorAligned of 2nd FFT background powers
 * \param [in]  rng      Pointer to GSL random number generator
 * \return Status value
 */
INT4 genIhsFar(ihsfarStruct *output, const UserInput_t *params, const UINT4 rows, const REAL4VectorAligned *aveNoise, const gsl_rng *rng)
{

   XLAL_CHECK( output != NULL && params != NULL && rows > 0 && aveNoise != NULL, XLAL_EINVAL );

   fprintf(stderr, "Determining IHS FAR values... ");
   fprintf(LOG, "Determining IHS FAR values... ");

   REAL8 Tobs = params->Tobs;

   UINT4 trials = 5*rows;
   if (trials<1000) {
      trials = 1000;
   }
   if (trials>5000) {
      fprintf(stderr, "Warning: number of trials may be insufficient given the number of rows to sum\n");
      trials = 5000;
   }

   //Allocations for IHS values for the number of trials
   REAL4VectorAligned *ihsvector = NULL, *ihss = NULL, *noise = NULL;
   REAL4VectorAlignedArray *ihsvectorarray = NULL;
   INT4Vector *locs = NULL;
   XLAL_CHECK( (noise = XLALCreateREAL4VectorAligned(aveNoise->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (ihsvector = XLALCreateREAL4VectorAligned((INT4)floor((1.0/(REAL8)params->ihsfactor)*aveNoise->length)-5, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (ihss = XLALCreateREAL4VectorAligned(trials, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (ihsvectorarray = createREAL4VectorAlignedArray(trials, ihsvector->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (locs = XLALCreateINT4Vector(trials)) != NULL, XLAL_EFUNC );

   //Determine the locations of the harmonics of the earth's rotation in the IHS vector
   //Amplitude modulations caused by the varying antenna pattern can sometimes cause excess power, so we ignore these harmonics
   REAL8 dailyharmonic = Tobs/(24.0*3600.0);
   REAL8 siderealharmonic = Tobs/86164.0905;
   REAL8 dailyharmonic2 = dailyharmonic*2.0, dailyharmonic3 = dailyharmonic*3.0, dailyharmonic4 = dailyharmonic*4.0;
   REAL8 siderealharmonic2 = siderealharmonic*2.0, siderealharmonic3 = siderealharmonic*3.0, siderealharmonic4 = siderealharmonic*4.0;
   INT4Vector *markedharmonics = NULL;
   XLAL_CHECK( (markedharmonics = XLALCreateINT4Vector(aveNoise->length)) != NULL, XLAL_EFUNC );
   memset(markedharmonics->data, 0, sizeof(INT4)*markedharmonics->length);
   if (!params->noNotchHarmonics) {
      for (UINT4 ii=0; ii<markedharmonics->length; ii++) {
         if (fabs(dailyharmonic-(REAL8)ii)<=1.0 || fabs(dailyharmonic2-(REAL8)ii)<=1.0 || fabs(dailyharmonic3-(REAL8)ii)<=1.0 || fabs(dailyharmonic4-(REAL8)ii)<=1.0 || fabs(siderealharmonic-(REAL8)ii)<=1.0 || fabs(siderealharmonic2-(REAL8)ii)<=1.0 || fabs(siderealharmonic3-(REAL8)ii)<=1.0 || fabs(siderealharmonic4-(REAL8)ii)<=1.0) markedharmonics->data[ii] = 1;
      }
   }

   //Do the expected IHS values from the expected background here
   memcpy(noise->data, aveNoise->data, sizeof(REAL4)*aveNoise->length);
   for (UINT4 ii=0; ii<aveNoise->length; ii++) if (markedharmonics->data[ii]==1) noise->data[ii] = 0.0;
   if (!params->weightedIHS) XLAL_CHECK( incHarmSumVector(output->expectedIHSVector, noise, params->ihsfactor) == XLAL_SUCCESS, XLAL_EFUNC );
   else XLAL_CHECK( incHarmSumVectorWeighted(output->expectedIHSVector, noise, aveNoise, params->ihsfactor) == XLAL_SUCCESS, XLAL_EFUNC );

   //Now do a number of trials
   for (UINT4 ii=0; ii<trials; ii++) {
      //Make exponential noise, removing harmonics of 24 hours to match with the same method as real analysis
      for (UINT4 jj=0; jj<aveNoise->length; jj++) {
         if (markedharmonics->data[jj]==0) noise->data[jj] = (REAL4)(gsl_ran_exponential(rng, aveNoise->data[jj]));
         else noise->data[jj] = 0.0;
      } /* for jj < aveNoise->length */

      //Make a random number 1 +/- 0.2 to create the variations in the nosie that we typically observe
      //This number needs to be positive
      REAL8 randval = 1.0;
      randval = 1.0 + gsl_ran_gaussian(rng, 0.2);
      while (randval<=0.0 || randval>=2.0) randval = 1.0 + gsl_ran_gaussian(rng, 0.2);     //limit range of variation

      //scale the exponential noise
      XLAL_CHECK( XLALVectorScaleREAL4(noise->data, randval, noise->data, noise->length) == XLAL_SUCCESS, XLAL_EFUNC );

      //Compute IHS value on exponential noise
      if (!params->weightedIHS) XLAL_CHECK( incHarmSumVector(ihsvector, noise, params->ihsfactor) == XLAL_SUCCESS, XLAL_EFUNC );
      else XLAL_CHECK( incHarmSumVectorWeighted(ihsvector, noise, aveNoise, params->ihsfactor) == XLAL_SUCCESS, XLAL_EFUNC );

      //Copy the result into the IHS vector array
      memcpy(ihsvectorarray->data[ii]->data, ihsvector->data, sizeof(REAL4)*ihsvector->length);
   } /* for ii < trials */

   //Destroy stuff
   XLALDestroyREAL4VectorAligned(noise);
   XLALDestroyREAL4VectorAligned(ihsvector);
   XLALDestroyINT4Vector(markedharmonics);

   //Create a fake vector with the same average value in each bin = 1.0
   REAL4VectorAligned *FbinMean = NULL;
   XLAL_CHECK( (FbinMean = XLALCreateREAL4VectorAligned(trials, 32)) != NULL, XLAL_EFUNC );
   for (UINT4 ii=0; ii<trials; ii++) FbinMean->data[ii] = 1.0;

   //Calculate the IHS sum values for the IHS trials
   XLAL_CHECK( sumIHSarrayFAR(output, ihsvectorarray, rows, FbinMean, params, rng) == XLAL_SUCCESS, XLAL_EFUNC );

   //Destroy stuff
   destroyREAL4VectorAlignedArray(ihsvectorarray);
   XLALDestroyREAL4VectorAligned(ihss);
   XLALDestroyREAL4VectorAligned(FbinMean);
   XLALDestroyINT4Vector(locs);

   fprintf(stderr, "done\n");
   fprintf(LOG, "done.\n");

   return XLAL_SUCCESS;

} /* genIhsFar() */


/**
 * Compute the IHS sums for a number of rows used for the FAR calculation
 * \param [out] outputfar      Pointer to the output ihsfarStruct
 * \param [in]  ihsvectorarray Pointer to REAL4VectorAlignedArray to be summed
 * \param [in]  rows           Number of neighbors to sum
 * \param [in]  FbinMean       Pointer to REAL4VectorAligned of normalized SFT background powers
 * \param [in]  params         Pointer to UserInput_t
 * \param [in]  rng            Pointer to random number generator
 * \return Status value
 */
INT4 sumIHSarrayFAR(ihsfarStruct *outputfar, REAL4VectorAlignedArray *ihsvectorarray, const UINT4 rows, const REAL4VectorAligned *FbinMean, const UserInput_t *params, const gsl_rng *rng)
{

   XLAL_CHECK( outputfar != NULL && ihsvectorarray != NULL && rows > 0 && FbinMean != NULL && params != NULL, XLAL_EINVAL );

   //The minimum and maximum index to search in the IHS vector
   INT4 maxIndexForIHS = (INT4)ceil(fmin(params->Tobs/params->Pmin, fmin(params->Tobs/minPeriod(params->dfmin, params->Tsft), params->Tobs/(4.0*params->Tsft)))) - 5;
   INT4 minIndexForIHS = (INT4)floor(fmax(5.0, params->Tobs/params->Pmax)) - 5;

   //Allocate a vector array that holds the summed values of at least two nearest neighbor rows
   //On the first iteration this holds the nearest neighbor sums, but on subsequent iterations, this holds nearest 3 neighbors
   //sums, then nearest 4 neighbor sums, and so on.
   REAL4VectorAlignedArray *tworows = NULL;
   XLAL_CHECK( (tworows = createREAL4VectorAlignedArray(ihsvectorarray->length-1, ihsvectorarray->data[0]->length, 32)) != NULL, XLAL_EFUNC );
   for (UINT4 ii=0; ii<tworows->length; ii++) memset(tworows->data[ii]->data, 0, sizeof(REAL4)*tworows->data[ii]->length);  //Set everything to 0 at the start

   //Allocate vectors of the ihs values and locations of the maximums
   //Vectors for values above the noise and scaling the noise
   REAL4VectorAligned *ihsvalues = NULL;
   INT4Vector *ihslocations = NULL;
   XLAL_CHECK( (ihsvalues = XLALCreateREAL4VectorAligned(ihsvectorarray->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (ihslocations = XLALCreateINT4Vector(ihsvectorarray->length)) != NULL, XLAL_EFUNC );
   REAL4VectorAligned *excessabovenoise = NULL, *scaledExpectedIHSVectorValues = NULL;
   XLAL_CHECK( (excessabovenoise = XLALCreateREAL4VectorAligned(ihsvectorarray->data[0]->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (scaledExpectedIHSVectorValues = XLALCreateREAL4VectorAligned(ihsvectorarray->data[0]->length, 32)) != NULL, XLAL_EFUNC );

   //Finding the maximum for each IHS vector and the location
   for (UINT4 ii=0; ii<ihsvalues->length; ii++) {
      //Scale the expected IHS vector by the value in FbinMean (in this function, it is 1.0)
      XLAL_CHECK( XLALVectorScaleREAL4(scaledExpectedIHSVectorValues->data, FbinMean->data[ii], outputfar->expectedIHSVector->data, outputfar->expectedIHSVector->length) == XLAL_SUCCESS, XLAL_EFUNC );
      //subtract the noise from the data
      XLAL_CHECK( VectorSubtractREAL4(excessabovenoise, ihsvectorarray->data[ii], scaledExpectedIHSVectorValues, params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
      //Find the location of the maximum IHS value
      ihslocations->data[ii] = max_index_in_range(excessabovenoise, minIndexForIHS, maxIndexForIHS) + 5;
      //And get the value
      ihsvalues->data[ii] = ihsvectorarray->data[ii]->data[ihslocations->data[ii]-5];
      //fprintf(stderr, "%d %f\n", ihslocations->data[ii], ihsvalues->data[ii]);
   } /* for ii=0 --> ihsvalues->length */

   //Some useful variables
   INT4Vector *rowarraylocs = NULL;
   REAL4VectorAligned *foms = NULL;

   //Starting from a minimum of 2 rows, start determining the FAR for each nearest neighbor sum, up to the maximum number of
   //rows to be summed
   for (UINT4 ii=2; ii<=rows; ii++) {
      //First allocate the necessary vectors
      XLAL_CHECK( (rowarraylocs = XLALCreateINT4Vector(ii)) != NULL, XLAL_EFUNC );
      XLAL_CHECK( (foms = XLALCreateREAL4VectorAligned(ihsvectorarray->length-(ii-1), 32)) != NULL, XLAL_EFUNC );

      //If the user has specified that we should use SSE operations, then do the nearest neighbor summing.
      //The result is in the tworows variable
      if (params->vectorMath!=0) {
         if (ii>2) XLAL_CHECK( VectorArraySum(tworows, tworows, ihsvectorarray, 0, ii-1, 0, (INT4)(ihsvectorarray->length-(ii-1))) == XLAL_SUCCESS, XLAL_EFUNC );
         else XLAL_CHECK( VectorArraySum(tworows, ihsvectorarray, ihsvectorarray, 0, ii-1, 0, (INT4)(ihsvectorarray->length-(ii-1))) == XLAL_SUCCESS, XLAL_EFUNC );
      }

      //Now we are going to loop through the input ihsvectorarray up to the number of rows-1
      for (UINT4 jj=0; jj<ihsvectorarray->length-(ii-1); jj++) {
         //A bit tricky here: sum IHS values across SFT frequency bins into the tworows variable only if we didn't do it with SSE above
         if (params->vectorMath==0) {
            if (ii>2) for (UINT4 kk=0; kk<tworows->data[jj]->length; kk++) tworows->data[jj]->data[kk] += ihsvectorarray->data[ii-1+jj]->data[kk];
            else for (UINT4 kk=0; kk<tworows->data[jj]->length; kk++) tworows->data[jj]->data[kk] = ihsvectorarray->data[jj]->data[kk] + ihsvectorarray->data[jj+1]->data[kk];
         }

         //Compute IHS FOM value
         memcpy(rowarraylocs->data, &(ihslocations->data[jj]), sizeof(INT4)*ii);  //First copy the necessary values
         foms->data[jj] = ihsFOM(rowarraylocs, (INT4)outputfar->expectedIHSVector->length);  //then compute the FOM
         XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
      } /* for jj< ihsvectorarray->length - (ii-1) */

      //Sample the IHS values that have been summed to compute mean, standard deviation, and FAR threshold values.
      //We have an if-else statement for when there are fewer than 10000 entries that will be in the tworows varaible
      //(only considering the number of rows we have summed together).
      REAL4VectorAligned *sampledtempihsvals = NULL;
      if ((ihsvectorarray->length-(ii-1))*ihsvectorarray->data[0]->length>10000) {
         //We sample the tworows array (up to the number of rows-1) without accepting any zeros
         //because zeros would come from the notched harmonics, which we won't want anyway
         sampledtempihsvals = sampleREAL4VectorAlignedArray_nozerosaccepted(tworows, ihsvectorarray->length-(ii-1), 10000, rng);
         outputfar->ihsdistMean->data[ii-2] = calcMean(sampledtempihsvals);  //And then calculate the mean value
         XLAL_CHECK( calcStddev(&(outputfar->ihsdistSigma->data[ii-2]), sampledtempihsvals) == XLAL_SUCCESS, XLAL_EFUNC );  //We also calculate the standard deviation
      } else {
         //If there were fewer than 10000 entries, then we will keep all those that are part of the nearest neighbor
         //sum up to this point
         XLAL_CHECK( (sampledtempihsvals = XLALCreateREAL4VectorAligned((ihsvectorarray->length-(ii-1))*ihsvectorarray->data[0]->length, 32)) != NULL, XLAL_EFUNC );
         for (UINT4 jj=0; jj<(ihsvectorarray->length-(ii-1)); jj++) memcpy(sampledtempihsvals->data, tworows->data[jj]->data, sizeof(REAL4)*tworows->data[0]->length);
         outputfar->ihsdistMean->data[ii-2] = calcMean_ignoreZeros(sampledtempihsvals);  //Calculate the mean value (remember, don't accept zeros)
         XLAL_CHECK( calcStddev_ignoreZeros(&(outputfar->ihsdistSigma->data[ii-2]), sampledtempihsvals) == XLAL_SUCCESS, XLAL_EFUNC );  //We also calculate the standard deviation
      }

      //If the user has specified the IHS FAR == 1.0, then we don't need to compute the threshold (it is = 0.0)
      INT4 numvals = 0;
      REAL8 farave = 0.0;
      if (params->ihsfar != 1.0) {
         //Looping through the sampled values, we are going to compute the average FAR
         for (UINT4 jj=0; jj<sampledtempihsvals->length; jj++) {
            //When the user has not specified using faster chisq inversion, use the GSL function
            if (!params->fastchisqinv && sampledtempihsvals->data[jj]!=0.0) {
               numvals++;
               farave += 0.5*gsl_cdf_chisq_Qinv(params->ihsfar, 2.0*sampledtempihsvals->data[jj]);
            } else if (params->fastchisqinv && sampledtempihsvals->data[jj]!=0.0) {
               numvals++;
               farave += 0.5*cdf_chisq_Qinv(params->ihsfar, 2.0*sampledtempihsvals->data[jj]);
               XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
            } //fastchisqinv?
         } //for jj=0 --> sampledtempihsval->length
         outputfar->ihsfar->data[ii-2] = farave/(REAL8)numvals;
      } //if params->ihsfar != 1.0
      else {
         outputfar->ihsfar->data[ii-2] = 0.0;
      }

      XLALDestroyREAL4VectorAligned(sampledtempihsvals);

      //FOM part
      outputfar->ihsfomdistMean->data[ii-2] = calcMean(foms);
      XLAL_CHECK( calcStddev(&(outputfar->ihsfomdistSigma->data[ii-2]), foms) == XLAL_SUCCESS, XLAL_EFUNC );
      //We need to find the smallest values
      if (params->ihsfomfar!=1.0 && params->ihsfom==0.0) {
         REAL4VectorAligned *smallestfomvals = NULL;
         XLAL_CHECK( (smallestfomvals = XLALCreateREAL4VectorAligned((INT4)round((ihsvalues->length-ii+1)*params->ihsfomfar)+1, 32)) != NULL, XLAL_EFUNC );
         XLAL_CHECK( sort_float_smallest(smallestfomvals, foms) == XLAL_SUCCESS, XLAL_EFUNC );  //Sort, finding the smallest values
         outputfar->fomfarthresh->data[ii-2] = smallestfomvals->data[smallestfomvals->length-1];  //Pick off the last element
         XLALDestroyREAL4VectorAligned(smallestfomvals);
      } else if (params->ihsfom!=0.0) outputfar->fomfarthresh->data[ii-2] = params->ihsfom;
      else outputfar->fomfarthresh->data[ii-2] = -1.0;

      XLALDestroyINT4Vector(rowarraylocs);
      rowarraylocs = NULL;
      XLALDestroyREAL4VectorAligned(foms);
      foms = NULL;
   } /* for ii <= rows */

   destroyREAL4VectorAlignedArray(tworows);
   XLALDestroyREAL4VectorAligned(ihsvalues);
   XLALDestroyREAL4VectorAligned(excessabovenoise);
   XLALDestroyREAL4VectorAligned(scaledExpectedIHSVectorValues);
   XLALDestroyINT4Vector(ihslocations);

   return XLAL_SUCCESS;

} /*sumIHSarrayFAR() */


/**
 * \brief Compute the IHS sums for a number of rows
 *
 * In the function we will select the the location which is the maximum above the noise
 *
 * \param [out] output         Pointer to the output ihsMaximaStruct
 * \param [in]  inputfar       Pointer to ihsfarStruct
 * \param [in]  ihsvectorarray Pointer to REAL4VectorAlignedArray to be summed
 * \param [in]  rows           Number of neighbors to sum
 * \param [in]  FbinMean       Pointer to REAL4VectorAligned of normalized SFT background powers
 * \param [in]  params         Pointer to UserInput_t
 * \return Status value
 */
INT4 sumIHSarray(ihsMaximaStruct *output, const ihsfarStruct *inputfar, REAL4VectorAlignedArray *ihsvectorarray, const UINT4 rows, const REAL4VectorAligned *FbinMean, const UserInput_t *params)
{

   XLAL_CHECK( output != NULL && inputfar != NULL && ihsvectorarray != NULL && rows > 0 && FbinMean != NULL && params != NULL, XLAL_EINVAL );

   //Again, we start off by allocating a "tworows" vector array of IHS nearest neighbor sums
   REAL4VectorAlignedArray *tworows = NULL;
   XLAL_CHECK( (tworows = createREAL4VectorAlignedArray(ihsvectorarray->length-1, ihsvectorarray->data[0]->length, 32)) != NULL, XLAL_EFUNC );
   for (UINT4 ii=0; ii<tworows->length; ii++) memset(tworows->data[ii]->data, 0, sizeof(REAL4)*tworows->data[0]->length);      //Set everything to 0.0

   //Allocation of ihs values and locations
   //Vectors for values above the noise and scaling the noise
   REAL4VectorAligned *ihsvalues = NULL, *excessabovenoise = NULL, *scaledExpectedIHSVectorValues = NULL;
   INT4Vector *ihslocations = NULL;
   XLAL_CHECK( (ihsvalues = XLALCreateREAL4VectorAligned(ihsvectorarray->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (excessabovenoise = XLALCreateREAL4VectorAligned(ihsvectorarray->data[0]->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (scaledExpectedIHSVectorValues = XLALCreateREAL4VectorAligned(ihsvectorarray->data[0]->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (ihslocations = XLALCreateINT4Vector(ihsvectorarray->length)) != NULL, XLAL_EFUNC );

   //The minimum and maximum index to search in the IHS vector
   INT4 maxIndexForIHS = (INT4)ceil(fmin(params->Tobs/params->Pmin, fmin(params->Tobs/minPeriod(params->dfmin, params->Tsft), params->Tobs/(4.0*params->Tsft)))) - 5;
   INT4 minIndexForIHS = (INT4)floor(fmax(5.0, params->Tobs/params->Pmax)) - 5;

   //Finding the maximum for each IHS vector and the location using SSE or not
   for (UINT4 ii=0; ii<ihsvalues->length; ii++) {
      //Scale the expected IHS vector by the FbinMean data value
      XLAL_CHECK( XLALVectorScaleREAL4(scaledExpectedIHSVectorValues->data, FbinMean->data[ii], inputfar->expectedIHSVector->data, inputfar->expectedIHSVector->length) == XLAL_SUCCESS, XLAL_EFUNC );
      //subtract the noise from the data
      XLAL_CHECK( VectorSubtractREAL4(excessabovenoise, ihsvectorarray->data[ii], scaledExpectedIHSVectorValues, params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );

      //search over the range of Pmin-->Pmax and higher harmonics the user has specified
      for (INT4 jj=0; jj<params->harmonicNumToSearch; jj++) {
         if (jj==0) {
            ihslocations->data[ii] = max_index_in_range(excessabovenoise, minIndexForIHS, maxIndexForIHS) + 5;
            ihsvalues->data[ii] = ihsvectorarray->data[ii]->data[ihslocations->data[ii]-5];
         } else {
            INT4 newIHSlocation = max_index_in_range(excessabovenoise, (jj+1)*minIndexForIHS, (jj+1)*maxIndexForIHS) + 5;
            REAL4 newIHSvalue = ihsvectorarray->data[ii]->data[newIHSlocation-5];
            if (newIHSvalue > ihsvalues->data[ii]) {
               ihslocations->data[ii] = newIHSlocation;
               ihsvalues->data[ii] = newIHSvalue;
            } /* if the new value is better than the previous value */
         }
      } /* for jj=0 --> jj<harmonicNumToSearch */
   }

   //Useful variables
   INT4Vector *rowarraylocs = NULL;

   //Start with the single IHS vector and march up with nearest neighbor sums up to the total number of row sums
   for (UINT4 ii=1; ii<=rows; ii++) {
      if (ii==1) {
         //Copy the data into the output
         memcpy(output->maximaForEachFbin->data, ihsvalues->data, sizeof(REAL4)*ihsvalues->length);
         memcpy(output->locationsForEachFbin->data, ihslocations->data, sizeof(INT4)*ihslocations->length);
      } else {
         //For everything 2 nearest neighbors and higher summed
         XLAL_CHECK( (rowarraylocs = XLALCreateINT4Vector(ii)) != NULL, XLAL_EFUNC );

         //The maximum index to search in the IHS vector
         maxIndexForIHS = (INT4)ceil(fmin(params->Tobs/params->Pmin, fmin(params->Tobs/minPeriod(0.5*(ii-1)/params->Tsft, params->Tsft), params->Tobs/(4.0*params->Tsft)))) - 5;

         REAL4 sumofnoise = 0.0;    //To scale the expected IHS background
         INT4 endloc = ((ii-1)*(ii-1)-(ii-1))/2;

         //Sum up the IHS vectors using SSE functions
         if (params->vectorMath!=0) {
            if (ii>2) XLAL_CHECK( VectorArraySum(tworows, tworows, ihsvectorarray, 0, ii-1, 0, (INT4)(ihsvectorarray->length-(ii-1))) == XLAL_SUCCESS, XLAL_EFUNC );
            else XLAL_CHECK( VectorArraySum(tworows, ihsvectorarray, ihsvectorarray, 0, ii-1, 0, (INT4)(ihsvectorarray->length-(ii-1))) == XLAL_SUCCESS, XLAL_EFUNC );
         }

         //Loop through the IHS vector neighbor sums
         for (UINT4 jj=0; jj<ihsvectorarray->length-(ii-1); jj++) {
            //If we didn't use SSE to sum the vector array (see lines above)
            if (params->vectorMath==0) {
               if (ii>2) for (UINT4 kk=0; kk<tworows->data[jj]->length; kk++) tworows->data[jj]->data[kk] += ihsvectorarray->data[ii-1+jj]->data[kk];
               else for (UINT4 kk=0; kk<tworows->data[jj]->length; kk++) tworows->data[jj]->data[kk] = ihsvectorarray->data[jj]->data[kk] + ihsvectorarray->data[jj+1]->data[kk];
            }

            //To scale the background efficiently
            if (jj==0) for (UINT4 kk=0; kk<ii; kk++) sumofnoise += FbinMean->data[kk];
            else {
               sumofnoise -= FbinMean->data[jj-1];
               sumofnoise += FbinMean->data[jj+(ii-1)];
            }

            //If using SSE, scale the expected IHS vector, subtract the noise from the data
            XLAL_CHECK( XLALVectorScaleREAL4(scaledExpectedIHSVectorValues->data, sumofnoise, inputfar->expectedIHSVector->data, inputfar->expectedIHSVector->length) == XLAL_SUCCESS, XLAL_EFUNC );
            XLAL_CHECK( VectorSubtractREAL4(excessabovenoise, tworows->data[jj], scaledExpectedIHSVectorValues, params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );

            //Compute the maximum IHS value in the second FFT frequency direction
            //search over the range of Pmin-->Pmax and higher harmonics the user has specified
            for (INT4 kk=0; kk<params->harmonicNumToSearch; kk++) {
               if (kk==0) {
                  output->locations->data[(ii-2)*ihsvalues->length-endloc+jj] = max_index_in_range(excessabovenoise, minIndexForIHS, maxIndexForIHS) + 5;
                  output->maxima->data[(ii-2)*ihsvalues->length-endloc+jj] = tworows->data[jj]->data[(output->locations->data[(ii-2)*ihsvalues->length-endloc+jj]-5)];
               } else {
                  INT4 newIHSlocation = max_index_in_range(excessabovenoise, (kk+1)*minIndexForIHS, (kk+1)*maxIndexForIHS) + 5;
                  REAL4 newIHSvalue = tworows->data[ii]->data[newIHSlocation-5];
                  if (newIHSvalue > output->maxima->data[(ii-2)*ihsvalues->length-endloc+jj]) {
                     output->locations->data[(ii-2)*ihsvalues->length-endloc+jj] = newIHSlocation;
                     output->maxima->data[(ii-2)*ihsvalues->length-endloc+jj] = newIHSvalue;
                  } /* if the new value is better than the previous value */
               }
            } /* for kk=0 --> kk<harmonicNumToSearch */

            memcpy(rowarraylocs->data, &(ihslocations->data[jj]), sizeof(INT4)*ii);
            output->foms->data[(ii-2)*ihsvalues->length-endloc+jj] = ihsFOM(rowarraylocs, (INT4)inputfar->expectedIHSVector->length);
            XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
         } /* for jj< ihsvectorarray->length - (ii-1) */

         XLALDestroyINT4Vector(rowarraylocs);
         rowarraylocs = NULL;
      }

   } /* for ii <= rows */

   destroyREAL4VectorAlignedArray(tworows);
   XLALDestroyREAL4VectorAligned(scaledExpectedIHSVectorValues);
   XLALDestroyREAL4VectorAligned(excessabovenoise);
   XLALDestroyREAL4VectorAligned(ihsvalues);
   XLALDestroyINT4Vector(ihslocations);

   return XLAL_SUCCESS;

} /*sumIHSarray() */


/**
 * Calculate the IHS FOM for a number of rows
 * \param [in] locs    Pointer to INT4Vector of location values
 * \param [in] fomnorm Normalization value and is the number of XXX
 * \return IHS FOM value
 */
REAL4 ihsFOM(const INT4Vector *locs, const INT4 fomnorm)
{

   XLAL_CHECK_REAL4( locs != NULL && fomnorm > 0.0, XLAL_EINVAL );

   REAL4 fom = 0.0;

   for (INT4 ii=0; ii<(INT4)(locs->length*0.5); ii++) fom += (REAL4)((locs->data[ii]-locs->data[locs->length-ii-1])*(locs->data[ii]-locs->data[locs->length-ii-1]))/(fomnorm*fomnorm);

   return fom;

} /* ihsFOM() */


/**
 * Find IHS candidates above thresholds
 * \param [out] candlist     Pointer to a pointer containing the candidate list
 * \param [in]  ihsfarstruct Pointer to ihsfarStruct
 * \param [in]  params       Pointer to UserInput_t
 * \param [in]  ffdata       Pointer to ffdataStruct
 * \param [in]  ihsmaxima    Pointer to ihsMaximaStruct containing the data to be tested above thresholds
 * \param [in]  fbinavgs     Pointer to REAL4VectorAligned of normalized SFT background powers
 * \param [in]  trackedlines Pointer to REAL4VectorSequence of lines (allowed to be NULL if no lines)
 * \return Status value
 */
INT4 findIHScandidates(candidateVector **candlist, const ihsfarStruct *ihsfarstruct, const UserInput_t *params, const ffdataStruct *ffdata, const ihsMaximaStruct *ihsmaxima, const REAL4VectorAligned *fbinavgs, const REAL4VectorSequence *trackedlines)
{

   XLAL_CHECK( ihsfarstruct != NULL && params != NULL && ffdata != NULL && ihsmaxima != NULL && fbinavgs != NULL, XLAL_EINVAL );

   REAL8 fsig, per0, B;
   INT4 numberofIHSvalsChecked = 0, numberofIHSvalsExceededThresh = 0, numberPassingBoth = 0, linesinterferewithnum = 0, skipped = 0, notskipped = 0;

   INT4 numfbins = ffdata->numfbins;  //number of fbins
   INT4 minrows = (INT4)round(2.0*params->dfmin*params->Tsft)+1;  //minimum number of sequential rows to search

   REAL4VectorAligned *ihss = NULL, *avgsinrange = NULL;
   INT4Vector *locs = NULL;

   //Check the IHS values against the FAR, checking between IHS width values
   //FILE *IHSVALSOUTPUT = fopen("./output/allihsvalspassthresh.dat","w");
   for (UINT4 ii=minrows; ii<=ihsfarstruct->ihsfar->length+1; ii++) {
      XLAL_CHECK( (ihss = XLALCreateREAL4VectorAligned(ii, 32)) != NULL, XLAL_EFUNC );
      XLAL_CHECK( (locs = XLALCreateINT4Vector(ii)) != NULL, XLAL_EFUNC );
      XLAL_CHECK( (avgsinrange = XLALCreateREAL4VectorAligned(ii, 32)) != NULL, XLAL_EFUNC );

      REAL8 highestval = 0.0, highestsignificance = 0.0; //highestvalnoise = 0.0
      INT4 highestvalloc = -1, jjloc = 0;
      for (UINT4 jj=0; jj<numfbins-(ii-1); jj++) {

         //Noise in the range of the rows, mean for IHS
         memcpy(avgsinrange->data, &(fbinavgs->data[jj]), sizeof(REAL4)*ii);
         REAL4 meanNoise = calcMean(avgsinrange);

         numberofIHSvalsChecked++;  //count the number of ihs values checked

         INT4 locationinmaximastruct = (ii-2)*numfbins-((ii-1)*(ii-1)-(ii-1))/2+jj;  //the location in the ihsmaximastruct

         //Check the IHS sum against the FAR (scaling FAR with mean of the noise in the range of rows)
         if (ihsmaxima->maxima->data[locationinmaximastruct] > ihsfarstruct->ihsfar->data[ii-2]*meanNoise) {

            numberofIHSvalsExceededThresh++;  //count the number of values exceeding the IHS value threshold

            if (ihsfarstruct->fomfarthresh->data[ii-2]==-1.0 || ihsmaxima->foms->data[locationinmaximastruct]<=ihsfarstruct->fomfarthresh->data[ii-2]) {

              numberPassingBoth++;  //count number passing both IHS value and FOM value thresholds

               INT4 loc = ihsmaxima->locations->data[locationinmaximastruct];
               per0 = params->Tobs/loc;                                                         //Candidate period
               fsig = params->fmin - params->dfmax + ((0.5*(ii-1) + jj) - 6.0)/params->Tsft;    //Candidate frequency
               B = 0.5*(ii-1)/params->Tsft;                                                     //Candidate modulation depth

               //Test to see if any tracked lines are overlapping the candidate signal
               INT4 nolinesinterfering = 1;
               if (trackedlines!=NULL) {
                  UINT4 kk = 0;
                  while (kk<trackedlines->length && nolinesinterfering==1) {
                     if (2.0*B>=(trackedlines->data[kk*3+2]-trackedlines->data[kk*3+1])) {
                        if ((trackedlines->data[kk*3+2]>=(REAL4)(fsig-B) && trackedlines->data[kk*3+2]<=(REAL4)(fsig+B)) ||
                            (trackedlines->data[kk*3+1]>=(REAL4)(fsig-B) && trackedlines->data[kk*3+1]<=(REAL4)(fsig+B))) {
                           nolinesinterfering = 0;
                        }
                     } // if the band spanned by the line is smaller than the band spanned by the signal
                     else {
                        if (((REAL4)(fsig+B)>=trackedlines->data[kk*3+1] && (REAL4)(fsig+B)<=trackedlines->data[kk*3+2]) ||
                            ((REAL4)(fsig-B)>=trackedlines->data[kk*3+1] && (REAL4)(fsig-B)<=trackedlines->data[kk*3+2])) {
                           nolinesinterfering = 0;
                        }
                     } // instead if the band spanned by the line is larger than the band spanned by the signal
                     kk++;
                  } // while kk < trackedlines->length && nolinesinterfering==1
               } // if trackedlines != NULL

               if (!nolinesinterfering) {
                  linesinterferewithnum++;
               } else {
                  REAL8 noise = ihsfarstruct->ihsdistMean->data[ii-2];
                  REAL8 totalnoise = meanNoise*noise;

                  REAL8 significance = gsl_cdf_chisq_Q(2.0*ihsmaxima->maxima->data[locationinmaximastruct], 2.0*totalnoise);
                  if (significance==0.0) {
                     significance = log10(LAL_E)*ihsmaxima->maxima->data[locationinmaximastruct] - (totalnoise - 1.0)*log10(ihsmaxima->maxima->data[locationinmaximastruct]) + lgamma(totalnoise)/log(10.0);
                  } else significance = -log10(significance);

                  if ( significance > highestsignificance && (params->followUpOutsideULrange || (!params->followUpOutsideULrange && fsig>=params->ULfmin && fsig<(params->ULfmin+params->ULfspan) && B>=params->ULminimumDeltaf && B<=params->ULmaximumDeltaf)) ) {
                     highestval = ihsmaxima->maxima->data[locationinmaximastruct]-totalnoise;
                     highestvalloc = locationinmaximastruct;
                     highestsignificance = significance;
                     jjloc = jj;
                     notskipped++;
                  } else {
                     skipped++;
                  }

               } // if no lines are interfering
            } // if fom is below or equal to threshold fom
         } // if val exceeds threshold
      } // for jj < numfbins-(ii-1)

      if (highestvalloc != -1) {
         INT4 loc = ihsmaxima->locations->data[highestvalloc];
         fsig = params->fmin - params->dfmax + (0.5*(ii-1) + jjloc - 6.0)/params->Tsft;  //Candidate frequency
         B = 0.5*(ii-1)/params->Tsft;                                                    //Candidate modulation depth
         per0 = params->Tobs/loc;                                                        //Candidate period
         REAL8 h0 = ihs2h0(2.0*highestval, params);  //Candidate h0, need factor of 2 for the degrees of freedom counting

         if ((*candlist)->numofcandidates == (*candlist)->length-1) XLAL_CHECK( (*candlist = resizecandidateVector(*candlist, 2*(*candlist)->length)) != NULL, XLAL_EFUNC );
         loadCandidateData(&((*candlist)->data[(*candlist)->numofcandidates]), fsig, per0, B, 0.0, 0.0, ihsmaxima->maxima->data[highestvalloc], h0, highestsignificance, 0, ffdata->tfnormalization, -1, 0);
         (*candlist)->numofcandidates++;
      }

      //Destroy
      XLALDestroyREAL4VectorAligned(ihss);
      XLALDestroyREAL4VectorAligned(avgsinrange);
      XLALDestroyINT4Vector(locs);
      ihss = NULL;
      locs = NULL;
      avgsinrange = NULL;

   } /* for ii < ihsfarstruct->ihsfar->length */

   //The outer loop is over "frequency" (the order in the ihs maxima vector, first 2 row sums, then three, and so on)
   //for (ii=0; ii<(INT4)numfbins; ii++) {
   /* for (ii=6+(INT4)round(params->dfmax*params->Tsft); ii<(numfbins-6-(INT4)round(params->dfmax*params->Tsft)); ii++) {
      REAL8 highestval = 0.0, highestsignificance = 0.0; //highestvalnoise = 0.0
      INT4 highestvalloc = -1, jjloc = 0;

      //This controls the maximum modulation depth to search which depends on the position in the "frequency" loop
      //INT4 maxrows = numfbins-ii;
      INT4 jjmax = (INT4)ihsfarstruct->ihsfar->length+1;
      //if (maxrows<(INT4)ihsfarstruct->ihsfar->length+1) jjmax = maxrows;
      //else jjmax = (INT4)ihsfarstruct->ihsfar->length+1;

      //Inner loop over modulation depth
      for (jj=minrows; jj<=jjmax; jj++) {
         numberofIHSvalsChecked++;

         ihss = XLALCreateREAL4Vector(jj);
         locs = XLALCreateINT4Vector(jj);
         avgsinrange = XLALCreateREAL4Vector(jj);
         if (ihss==NULL) {
            fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, jj);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         } else if (locs==NULL) {
            fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", __func__, jj);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         } else if (avgsinrange==NULL) {
            fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, jj);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }

         //Noise in the range of the rows, mean for IHS
         memcpy(avgsinrange->data, &(fbinavgs->data[ii]), sizeof(REAL4)*jj);
         REAL4 meanNoise = calcMean(avgsinrange);
         if (XLAL_IS_REAL4_FAIL_NAN(meanNoise)) {
            fprintf(stderr,"%s: calcMean() failed.\n", __func__);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }

         INT4 locationinmaximastruct = (jj-2)*numfbins-((jj-1)*(jj-1)-(jj-1))/2 + ii;

         //Check the IHS sum against the FAR (scaling FAR with mean of the noise in the range of rows)
         if (ihsmaxima->maxima->data[locationinmaximastruct] > ihsfarstruct->ihsfar->data[jj-2]*meanNoise) {
            numberofIHSvalsExceededThresh++;
            if (ihsfarstruct->fomfarthresh->data[jj-2]==-1.0 || ihsmaxima->foms->data[locationinmaximastruct]<=ihsfarstruct->fomfarthresh->data[jj-2]) {

               numberPassingBoth++;

               INT4 loc = ihsmaxima->locations->data[locationinmaximastruct];
               per0 = params->Tobs/loc;                                          //Candidate period
               fsig = params->fmin - params->dfmax + (0.5*(jj-1) + ii - 6)/params->Tsft;             //Candidate frequency
               B = 0.5*(jj-1)/params->Tsft;                                      //Candidate modulation depth

               //Test to see if any tracked lines are overlapping the candidate signal
               INT4 nolinesinterfering = 1;
               if (trackedlines!=NULL) {
                  kk = 0;
                  while (kk<(INT4)trackedlines->length && nolinesinterfering==1) {
                     if (2.0*B>=(trackedlines->data[kk*3+2]-trackedlines->data[kk*3+1])) {
                        if ((trackedlines->data[kk*3+2]>=fsig-B && trackedlines->data[kk*3+2]<=fsig+B) ||
                            (trackedlines->data[kk*3+1]>=fsig-B && trackedlines->data[kk*3+1]<=fsig+B)) {
                           nolinesinterfering = 0;
                        }
                     } // if the band spanned by the line is smaller than the band spanned by the signal
                     else {
                        if ((fsig+B>=trackedlines->data[kk*3+1] && fsig+B<=trackedlines->data[kk*3+2]) ||
                            (fsig-B>=trackedlines->data[kk*3+1] && fsig-B<=trackedlines->data[kk*3+2])) {
                           nolinesinterfering = 0;
                        }
                     } // instead if the band spanned by the line is larger than the band spanned by the signal
                     kk++;
                  } // while kk < trackedlines->length && nolinesinterfering==1
               } // if trackedlines != NULL

               if (!nolinesinterfering) {
                  linesinterferewithnum++;
               } else {
                  REAL8 noise = ihsfarstruct->ihsdistMean->data[jj-2];
                  REAL8 totalnoise = meanNoise*noise;
                  REAL8 sigma = calcRms(avgsinrange)*ihsfarstruct->ihsdistSigma->data[jj-2];

                  REAL8 significance = (2.0*ihsmaxima->maxima->data[locationinmaximastruct] - 2.0*totalnoise)/sqrt(2.0*2.0*sigma);

                  //if (ihsmaxima->maxima->data[locationinmaximastruct]-totalnoise > highestval) {
                  //if ( ihsmaxima->maxima->data[locationinmaximastruct]-totalnoise > highestval &&
                      //(params->followUpOutsideULrange || (!params->followUpOutsideULrange &&
                        //fsig>=params->ULfmin && fsig<=(params->ULfmin+params->ULfspan) &&
                        //B>=params->ULmindf && B<=params->ULmaxdf)) ) {
                  if ( significance > highestsignificance &&
                      (params->followUpOutsideULrange || (!params->followUpOutsideULrange &&
                      fsig>=params->ULfmin && fsig<=(params->ULfmin+params->ULfspan) &&
                      B>=params->ULmindf && B<=params->ULmaxdf)) ) {
                     //highestval = ihsmaxima->maxima->data[locationinmaximastruct];
                     highestval = ihsmaxima->maxima->data[locationinmaximastruct]-totalnoise;
                     //highestvalnoise = totalnoise;
                     highestvalloc = locationinmaximastruct;
                     highestsignificance = significance;
                     jjloc = jj;
                     notskipped++;
                  } else {
                     skipped++;
                  }
               } // if no lines are interfering
            } //If exceeding the FOM threshold
         } //If exceeding the IHS FAR threshold

         XLALDestroyREAL4Vector(ihss);
         XLALDestroyREAL4Vector(avgsinrange);
         XLALDestroyINT4Vector(locs);
      } //loop over modulation depths

      if (highestvalloc != -1) {
         INT4 loc = ihsmaxima->locations->data[highestvalloc];
         //Candidate frequency
         fsig = params->fmin - params->dfmax + (0.5*(jjloc-1) + ii - 6.0)/params->Tsft;
         //Candidate modulation depth
         B = 0.5*(jjloc-1)/params->Tsft;
         //Candidate period
         per0 = params->Tobs/loc;
         //Candidate h0
         REAL8 h0 = ihs2h0(2.0*highestval, params);  //Need factor of 2 for the degrees of freedom counting
         REAL8 significance = highestsignificance;
         //fprintf(stderr, "%d %d %f\n", ii, jjloc, significance);     //remove this

         if (candlist->numofcandidates == candlist->length-1) {
            candlist = resizecandidateVector(candlist, 2*(candlist->length));
            if (candlist->data==NULL) {
               fprintf(stderr,"%s: resizecandidateVector() failed.\n", __func__);
               XLAL_ERROR_VOID(XLAL_EFUNC);
            }
         }
         //loadCandidateData(&candlist->data[candlist->numofcandidates], fsig, per0, B, 0.0, 0.0, ihsmaxima->maxima->data[locationinmaximastruct], h0, 0.0, 0, sqrt(ffdata->tfnormalization/2.0*params->Tsft));
         loadCandidateData(&candlist->data[candlist->numofcandidates], fsig, per0, B, 0.0, 0.0, ihsmaxima->maxima->data[highestvalloc], h0, significance, 0, ffdata->tfnormalization);
         (candlist->numofcandidates)++;
      }
   } //loop over "frequency" */

   //fclose(IHSVALSOUTPUT);

   fprintf(stderr,"Number of IHS vals checked = %d, number exceeding IHS threshold = %d, number passing both = %d, but lines interfere with %d, number not skipped = %d and number of skipped candidates is %d\n", numberofIHSvalsChecked, numberofIHSvalsExceededThresh, numberPassingBoth, linesinterferewithnum, notskipped, skipped);
   fprintf(stderr, "Candidates found in IHS step = %d\n", (*candlist)->numofcandidates);
   fprintf(LOG,"Number of IHS vals checked = %d, number exceeding IHS threshold = %d, number passing both = %d, but lines interfere with %d, number not skipped = %d and number of skipped candidates is %d\n", numberofIHSvalsChecked, numberofIHSvalsExceededThresh, numberPassingBoth, linesinterferewithnum, notskipped, skipped);
   fprintf(LOG, "Candidates found in IHS step = %d\n", (*candlist)->numofcandidates);

   return XLAL_SUCCESS;

} /* findIHScandidates() */


/**
 * Convert the IHS statistic to an estimated h0, based on injections
 * \param [in] ihsval The IHS statistic value
 * \param [in] params Pointer to UserInput_t
 * \return Estimated h0
 */
REAL8 ihs2h0(const REAL8 ihsval, const UserInput_t *params)
{

   if (ihsval<=0.0) return 0.0;
   REAL8 prefact = 1.0;
   //prefact = 7.2;  //old value for when IHS --> h0 when in the upper limits, the value of the noise subtracted was incorrectly calculated and for >99.999% CL
   prefact = 44.7; //new value for IHS --> h0
   return prefact*pow(ihsval/(params->Tsft*params->Tobs),0.25);

}
