/*
*  Copyright (C) 2010, 2011, 2012, 2014 Evan Goetz
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

#include <math.h>
#include <time.h>

#include <lal/LALMalloc.h>
#include <lal/SeqFactories.h>
#include <lal/LALConstants.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "IHS.h"
#include "statistics.h"
#include "candidates.h"
#include "fastchisqinv.h"
#include "vectormath.h"
#include "TwoSpect.h"


/**
 * Create vectors for IHS maxima struct
 * \param [in] fbins Number of frequency bins
 * \param [in] rows  Number of neighboring rows to be summed
 * \return Pointer to a newly created ihsMaximaStruct
 */
ihsMaximaStruct * new_ihsMaxima(INT4 fbins, INT4 rows)
{

   XLAL_CHECK_NULL( fbins > 0 && rows > 0, XLAL_EINVAL );

   ihsMaximaStruct *ihsmaxima = NULL;
   XLAL_CHECK_NULL( (ihsmaxima = XLALMalloc(sizeof(*ihsmaxima))) != NULL, XLAL_ENOMEM );

   //The number of ihs maxima = Sum[fbins - (i-1), {i, 2, rows}]
   // = fbins*rows - fbins - (rows**2 - rows)/2
   INT4 numberofmaxima = fbins*rows - fbins - (INT4)((rows*rows-rows)/2);

   XLAL_CHECK_NULL( (ihsmaxima->maxima = XLALCreateREAL4Vector(numberofmaxima)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (ihsmaxima->locations = XLALCreateINT4Vector(numberofmaxima)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (ihsmaxima->foms = XLALCreateREAL4Vector(numberofmaxima)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (ihsmaxima->maximaForEachFbin = XLALCreateREAL4Vector(fbins)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (ihsmaxima->locationsForEachFbin = XLALCreateINT4Vector(fbins)) != NULL, XLAL_EFUNC );
   ihsmaxima->rows = rows;

   return ihsmaxima;

} /* new_ihsMaxima() */


/**
 * Destroy vectors and the IHS maxima struct
 * \param [in] data Pointer to an ihsMaximaStruct to be freed
 */
void free_ihsMaxima(ihsMaximaStruct *data)
{

   XLALDestroyREAL4Vector(data->maxima);
   XLALDestroyINT4Vector(data->locations);
   XLALDestroyREAL4Vector(data->foms);
   XLALDestroyREAL4Vector(data->maximaForEachFbin);
   XLALDestroyINT4Vector(data->locationsForEachFbin);
   XLALFree((ihsMaximaStruct*)data);

} /* free_ihsMaxima() */


/**
 * Run the IHS algorithm
 * \param [out] output      Pointer to the ihsMaximaStruct
 * \param [in]  input       Pointer to the ffdataStruct
 * \param [in]  ihsfarinput Pointer to the ihsfarStruct
 * \param [in]  params      Pointer to inputParamsStruct
 * \param [in]  rows        Number of neighboring rows to be summed
 * \param [in]  aveNoise    Pointer to a REAL4Vector of 2nd FFT background powers
 * \param [in]  FbinMean    Pointer to a REAL4Vector of normalized SFT background powers
 * \return Status value
 */
INT4 runIHS(ihsMaximaStruct *output, ffdataStruct *input, ihsfarStruct *ihsfarinput, inputParamsStruct *params, INT4 rows, REAL4Vector *aveNoise, REAL4Vector *FbinMean)
{

   XLAL_CHECK( output != NULL && input != NULL && ihsfarinput != NULL && params != NULL && rows > 0 && aveNoise != NULL && FbinMean != NULL, XLAL_EINVAL );

   INT4 ii, jj;

   INT4 numfbins = input->numfbins;
   INT4 numfprbins = input->numfprbins;

   //Allocate memory for the necessary vectors
   REAL4Vector *row = NULL, *ihss = NULL, *ihsvector = NULL;
   INT4Vector *locs = NULL;
   ihsVals *ihsvals = NULL;
   REAL4VectorSequence *ihsvectorsequence = NULL;
   XLAL_CHECK( (row = XLALCreateREAL4Vector(numfprbins)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (ihss = XLALCreateREAL4Vector(numfbins)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (ihsvector = XLALCreateREAL4Vector((INT4)floor((1.0/(REAL8)params->ihsfactor)*numfprbins)-5)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (locs = XLALCreateINT4Vector(numfbins)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (ihsvals = new_ihsVals()) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (ihsvectorsequence = XLALCreateREAL4VectorSequence(numfbins, ihsvector->length)) != NULL, XLAL_EFUNC );

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
      for (ii=0; ii<(INT4)markedharmonics->length; ii++) {
         if (fabs(dailyharmonic-(REAL8)ii)<=1.0 || fabs(dailyharmonic2-(REAL8)ii)<=1.0 || fabs(dailyharmonic3-(REAL8)ii)<=1.0 || fabs(dailyharmonic4-(REAL8)ii)<=1.0 || fabs(siderealharmonic-(REAL8)ii)<=1.0 || fabs(siderealharmonic2-(REAL8)ii)<=1.0 || fabs(siderealharmonic3-(REAL8)ii)<=1.0 || fabs(siderealharmonic4-(REAL8)ii)<=1.0) markedharmonics->data[ii] = 1;
      }
   }

   //Loop through the rows, 1 frequency at a time
   for (ii=0; ii<(INT4)ihss->length; ii++) {

      //For each row, populate it with the data for that frequency bin, excluding harmonics of antenna pattern modulation
      memcpy(row->data, &(input->ffdata->data[ii*numfprbins]), sizeof(REAL4)*numfprbins);
      if (!params->noNotchHarmonics) for (jj=0; jj<(INT4)row->length; jj++) if (markedharmonics->data[jj]==1) row->data[jj] = 0.0;

      //Run the IHS algorithm on the row
      if (!params->weightedIHS) XLAL_CHECK( incHarmSumVector(ihsvector, row, params->ihsfactor) == XLAL_SUCCESS, XLAL_EFUNC );
      else XLAL_CHECK( incHarmSumVectorWeighted(ihsvector, row, aveNoise, params->ihsfactor) == XLAL_SUCCESS, XLAL_EFUNC );

      //Copy the result into the ihsvector sequence
      memcpy(&(ihsvectorsequence->data[ii*ihsvector->length]), ihsvector->data, sizeof(REAL4)*ihsvector->length);

   } /* for ii < ihss->length */

   //Now do the summing of the IHS values
   XLAL_CHECK( sumIHSSequence(output, ihsfarinput, ihsvectorsequence, rows, FbinMean, params) == XLAL_SUCCESS, XLAL_EFUNC );

   //Destroy stuff
   XLALDestroyREAL4VectorSequence(ihsvectorsequence);
   XLALDestroyREAL4Vector(row);
   XLALDestroyREAL4Vector(ihss);
   XLALDestroyREAL4Vector(ihsvector);
   XLALDestroyINT4Vector(locs);
   XLALDestroyINT4Vector(markedharmonics);
   free_ihsVals(ihsvals);

   return XLAL_SUCCESS;

} /* runIHS() */


/**
 * Allocate memory for ihsVals struct
 * \return Pointer to a newly created ihsVals structure
 */
ihsVals * new_ihsVals(void)
{
   ihsVals *ihsvals = NULL;
   XLAL_CHECK_NULL( (ihsvals = XLALMalloc(sizeof(*ihsvals))) != NULL, XLAL_ENOMEM );
   return ihsvals;
} /* new_ihsVals() */


/**
 * Destroy ihsVals struct
 * \param [in] Pointer to an ihsVals structure to be freed
 */
void free_ihsVals(ihsVals *ihsvals)
{
   XLALFree((ihsVals*)ihsvals);
} /* free_ihsVals() */


/**
 * Compute the IHS sum maximum
 * \param [out] output    Pointer to the ihsVals structure
 * \param [in]  input     Pointer to a REAL4Vector
 * \param [in]  ihsfactor Number of folds of the 2nd FFT
 * \return Status value
 */
INT4 incHarmSum(ihsVals *output, REAL4Vector *input, INT4 ihsfactor)
{

   XLAL_CHECK( output != NULL && input != NULL && ihsfactor > 0, XLAL_EINVAL );

   INT4 ii;

   output->ihs = 0.0;

   REAL4Vector *tempvect = NULL;
   XLAL_CHECK( (tempvect = XLALCreateREAL4Vector((INT4)floor((1.0/(REAL8)ihsfactor)*input->length)-5)) != NULL, XLAL_EFUNC );

   XLAL_CHECK( incHarmSumVector(tempvect, input, ihsfactor) == XLAL_SUCCESS, XLAL_EFUNC );

   for (ii=0; ii<(INT4)tempvect->length; ii++) {
      if (tempvect->data[ii]>output->ihs) {
         output->ihs = tempvect->data[ii];
         output->loc = ii+5;
      }
   } /* for ii < tempvect->length */

   XLALDestroyREAL4Vector(tempvect);

   return XLAL_SUCCESS;

} /* incHarmSum() */


/**
 * Compute the IHS vector -- does not compute the maximum value
 * \param [out] output    Pointer to a REAL4Vector to contain the folded values
 * \param [in]  input     Pointer to a REAL4Vector for the input to the IHS
 * \param [in]  ihsfactor Number of folds of the 2nd FFT
 * \return Status value
 */
INT4 incHarmSumVector(REAL4Vector *output, REAL4Vector *input, INT4 ihsfactor)
{

   XLAL_CHECK( output != NULL && input != NULL && ihsfactor > 0, XLAL_EINVAL );

   INT4 ii, jj, highval = (INT4)floor((1.0/(REAL8)ihsfactor)*input->length);

   //From 5 up to the highval
   for (ii=5; ii<highval; ii++) {
      output->data[ii-5] = 0.0;  //first set the value to zero
      for (jj=1; jj<=ihsfactor; jj++) output->data[ii-5] += input->data[ii*jj];

      //check that we're not going outside of the allowed limits
      XLAL_CHECK( ihsfactor*(ii-1)<=(INT4)input->length-1, XLAL_EBADLEN );
   } /* for ii=5 --> highval */

   return XLAL_SUCCESS;

} /* incHarmSumVector() */


/**
 * Compute the noise weighted IHS vector -- does not compute the maximum value
 * \param [out] output    Pointer to a REAL4Vector to contain the folded values
 * \param [in]  input     Pointer to a REAL4Vector for the input to the IHS
 * \param [in]  aveNoise  Pointer to a REAL4Vector of 2nd FFT background powers
 * \param [in]  ihsfactor Number of folds of the 2nd FFT
 * \return Status value
 */
INT4 incHarmSumVectorWeighted(REAL4Vector *output, REAL4Vector *input, REAL4Vector *aveNoise, INT4 ihsfactor)
{

   XLAL_CHECK( output != NULL && input != NULL && aveNoise != NULL && ihsfactor > 0, XLAL_EINVAL );

   INT4 ii, jj, highval = (INT4)floor((1.0/(REAL8)ihsfactor)*input->length);

   //From 5 up to the highval
   for (ii=5; ii<highval; ii++) {
      output->data[ii-5] = 0.0;  //first set the value to zero
      for (jj=1; jj<=ihsfactor; jj++) {
         REAL4 weight = 1.0/(aveNoise->data[ii*jj]*aveNoise->data[ii*jj]);
         output->data[ii-5] += input->data[ii*jj]*weight;
      }

      //check that we're not going outside of the allowed limits
      XLAL_CHECK( ihsfactor*(ii-1)<=(INT4)input->length-1, XLAL_EBADLEN );
   } /* for ii=5 --> highval */

   return XLAL_SUCCESS;

} /* incHarmSumVectorWeighted() */


/**
 * Allocate memory for ihsfarStruct struct
 * \param [in] rows   Number of neighbors to sum
 * \param [in] params Pointer to inputParamsStruct
 * \return Pointer to a newly allocated ihsfarStruct
 */
ihsfarStruct * new_ihsfarStruct(INT4 rows, inputParamsStruct *params)
{

   XLAL_CHECK_NULL( rows > 0 && params != NULL, XLAL_EINVAL );

   ihsfarStruct *ihsfarstruct = NULL;
   XLAL_CHECK_NULL( (ihsfarstruct = XLALMalloc(sizeof(*ihsfarstruct))) != NULL, XLAL_ENOMEM );

   XLAL_CHECK_NULL( (ihsfarstruct->ihsfar = XLALCreateREAL4Vector(rows-1)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (ihsfarstruct->ihsdistMean = XLALCreateREAL4Vector(rows-1)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (ihsfarstruct->ihsdistSigma = XLALCreateREAL4Vector(rows-1)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (ihsfarstruct->fomfarthresh = XLALCreateREAL4Vector(rows-1)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (ihsfarstruct->ihsfomdistMean = XLALCreateREAL4Vector(rows-1)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (ihsfarstruct->ihsfomdistSigma = XLALCreateREAL4Vector(rows-1)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_NULL( (ihsfarstruct->expectedIHSVector = XLALCreateREAL4Vector((INT4)floor((1.0/(REAL8)params->ihsfactor)*((INT4)floor(floor(params->Tobs/(params->Tcoh-params->SFToverlap)-1)*0.5)+1))-5)) != NULL, XLAL_EFUNC );

   memset(ihsfarstruct->expectedIHSVector->data, 0, sizeof(REAL4)*ihsfarstruct->expectedIHSVector->length);

   return ihsfarstruct;

} /* new_ihsfarStruct() */


/**
 * Destroy ihsfarStruct struct
 * \param [in] ihsfarstruct Pointer to the ihsfarStruct to be destroyed
 */
void free_ihsfarStruct(ihsfarStruct *ihsfarstruct)
{
   XLALDestroyREAL4Vector(ihsfarstruct->ihsfar);
   XLALDestroyREAL4Vector(ihsfarstruct->ihsdistMean);
   XLALDestroyREAL4Vector(ihsfarstruct->ihsdistSigma);
   XLALDestroyREAL4Vector(ihsfarstruct->fomfarthresh);
   XLALDestroyREAL4Vector(ihsfarstruct->ihsfomdistMean);
   XLALDestroyREAL4Vector(ihsfarstruct->ihsfomdistSigma);
   XLALDestroyREAL4Vector(ihsfarstruct->expectedIHSVector);
   XLALFree((ihsfarStruct*)ihsfarstruct);
} /* free_ihsfarStruct() */


/**
 * Compute the IHS FAR for a sum of a number of rows
 * \param [out] output Pointer to the output ihsfarStruct
 * \param [in]  params Pointer to inputParamsStruct
 * \param [in]  rows   Number of neighbors to sum
 * \param [in]  aveNoise Pointer to REAL4Vector of 2nd FFT background powers
 * \return Status value
 */
INT4 genIhsFar(ihsfarStruct *output, inputParamsStruct *params, INT4 rows, REAL4Vector *aveNoise)
{

   XLAL_CHECK( output != NULL && params != NULL && rows > 0 && aveNoise != NULL, XLAL_EINVAL );

   fprintf(stderr, "Determining IHS FAR values... ");
   fprintf(LOG, "Determining IHS FAR values... ");

   INT4 ii, jj;
   REAL8 Tobs = params->Tobs;

   INT4 trials = 5*rows;
   if (trials<1000) {
      trials = 1000;
   }
   if (trials>5000) {
      fprintf(stderr, "Warning: number of trials may be insufficient given the number of rows to sum\n");
      trials = 5000;
   }

   //Allocations for IHS values for the number of trials
   REAL4Vector *noise = NULL, *ihsvector = NULL, *ihss = NULL;
   REAL4VectorSequence *ihsvectorsequence = NULL;
   INT4Vector *locs = NULL;
   XLAL_CHECK( (noise = XLALCreateREAL4Vector(aveNoise->length)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (ihsvector = XLALCreateREAL4Vector((INT4)floor((1.0/(REAL8)params->ihsfactor)*aveNoise->length)-5)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (ihss = XLALCreateREAL4Vector(trials)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (ihsvectorsequence = XLALCreateREAL4VectorSequence(trials, ihsvector->length)) != NULL, XLAL_EFUNC );
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
      for (ii=0; ii<(INT4)markedharmonics->length; ii++) {
         if (fabs(dailyharmonic-(REAL8)ii)<=1.0 || fabs(dailyharmonic2-(REAL8)ii)<=1.0 || fabs(dailyharmonic3-(REAL8)ii)<=1.0 || fabs(dailyharmonic4-(REAL8)ii)<=1.0 || fabs(siderealharmonic-(REAL8)ii)<=1.0 || fabs(siderealharmonic2-(REAL8)ii)<=1.0 || fabs(siderealharmonic3-(REAL8)ii)<=1.0 || fabs(siderealharmonic4-(REAL8)ii)<=1.0) markedharmonics->data[ii] = 1;
      }
   }

   //Do the expected IHS values from the expected background here
   memcpy(noise->data, aveNoise->data, sizeof(REAL4)*aveNoise->length);
   for (ii=0; ii<(INT4)aveNoise->length; ii++) if (markedharmonics->data[ii]==1) noise->data[ii] = 0.0;
   if (!params->weightedIHS) XLAL_CHECK( incHarmSumVector(output->expectedIHSVector, noise, params->ihsfactor) == XLAL_SUCCESS, XLAL_EFUNC );
   else XLAL_CHECK( incHarmSumVectorWeighted(output->expectedIHSVector, noise, aveNoise, params->ihsfactor) == XLAL_SUCCESS, XLAL_EFUNC );

   //Now do a number of trials
   for (ii=0; ii<trials; ii++) {
      //Make exponential noise, removing harmonics of 24 hours to match with the same method as real analysis
      for (jj=0; jj<(INT4)aveNoise->length; jj++) {
         if (markedharmonics->data[jj]==0) noise->data[jj] = (REAL4)(gsl_ran_exponential(params->rng, aveNoise->data[jj]));
         else noise->data[jj] = 0.0;
      } /* for jj < aveNoise->length */

      //Make a random number 1 +/- 0.2 to create the variations in the nosie that we typically observe
      //This number needs to be positive
      REAL8 randval = 1.0;
      randval = 1.0 + gsl_ran_gaussian(params->rng, 0.2);
      while (randval<=0.0 || randval>=2.0) randval = 1.0 + gsl_ran_gaussian(params->rng, 0.2);     //limit range of variation

      //scale the exponential noise
      if (params->useSSE)  XLAL_CHECK( sseScaleREAL4Vector(noise, noise, randval) == XLAL_SUCCESS, XLAL_EFUNC );
      else if (params->useAVX) XLAL_CHECK( avxScaleREAL4Vector(noise, noise, randval) == XLAL_SUCCESS, XLAL_EFUNC );
      else for (jj=0; jj<(INT4)noise->length; jj++) noise->data[jj] *= randval;

      //Compute IHS value on exponential noise
      if (!params->weightedIHS) XLAL_CHECK( incHarmSumVector(ihsvector, noise, params->ihsfactor) == XLAL_SUCCESS, XLAL_EFUNC );
      else XLAL_CHECK( incHarmSumVectorWeighted(ihsvector, noise, aveNoise, params->ihsfactor) == XLAL_SUCCESS, XLAL_EFUNC );

      //Copy the result into the IHS vector sequence
      memcpy(&(ihsvectorsequence->data[ii*ihsvector->length]), ihsvector->data, sizeof(REAL4)*ihsvector->length);
   } /* for ii < trials */

   //Destroy stuff
   XLALDestroyREAL4Vector(noise);
   XLALDestroyREAL4Vector(ihsvector);
   XLALDestroyINT4Vector(markedharmonics);

   //Create a fake vector with the same average value in each bin = 1.0
   REAL4Vector *FbinMean = NULL;
   XLAL_CHECK( (FbinMean = XLALCreateREAL4Vector(trials)) != NULL, XLAL_EFUNC );
   for (ii=0; ii<trials; ii++) FbinMean->data[ii] = 1.0;

   //Calculate the IHS sum values for the IHS trials
   XLAL_CHECK( sumIHSSequenceFAR(output, ihsvectorsequence, rows, FbinMean, params) == XLAL_SUCCESS, XLAL_EFUNC );

   //Destroy stuff
   XLALDestroyREAL4VectorSequence(ihsvectorsequence);
   XLALDestroyREAL4Vector(ihss);
   XLALDestroyREAL4Vector(FbinMean);
   XLALDestroyINT4Vector(locs);

   fprintf(stderr, "done\n");
   fprintf(LOG, "done.\n");

   return XLAL_SUCCESS;

} /* genIhsFar() */


/**
 * Compute the IHS sums for a number of rows used for the FAR calculation
 * \param [out] outputfar         Pointer to the output ihsfarStruct
 * \param [in]  ihsvectorsequence Pointer to REAL4VectorSequence to be summed
 * \param [in]  rows              Number of neighbors to sum
 * \param [in]  FbinMean          Pointer to REAL4Vector of normalized SFT background powers
 * \param [in]  params            Pointer to inputParamsStruct
 * \return Status value
 */
INT4 sumIHSSequenceFAR(ihsfarStruct *outputfar, REAL4VectorSequence *ihsvectorsequence, INT4 rows, REAL4Vector *FbinMean, inputParamsStruct *params)
{

   XLAL_CHECK( outputfar != NULL && ihsvectorsequence != NULL && rows > 0 && FbinMean != NULL && params != NULL, XLAL_EINVAL );

   INT4 ii = 0, jj;

   //The minimum and maximum index to search in the IHS vector
   INT4 maxIndexForIHS = (INT4)ceil(fmin(params->Tobs/params->Pmin, fmin(params->Tobs/minPeriod(params->dfmin, params->Tcoh), params->Tobs/(4.0*params->Tcoh)))) - 5;
   INT4 minIndexForIHS = (INT4)floor(fmax(5.0, params->Tobs/params->Pmax)) - 5;

   //Allocate a vector sequence that holds the summed values of at least two nearest neighbor rows
   //On the first iteration this holds the nearest neighbor sums, but on subsequent iterations, this holds nearest 3 neighbors
   //sums, then nearest 4 neighbor sums, and so on.
   REAL4VectorSequence *tworows = NULL;
   XLAL_CHECK( (tworows = XLALCreateREAL4VectorSequence(ihsvectorsequence->length-1, ihsvectorsequence->vectorLength)) != NULL, XLAL_EFUNC );
   memset(tworows->data, 0, sizeof(REAL4)*tworows->length*tworows->vectorLength);  //Set everything to 0 at the start

   //Allocate vectors of the ihs values and locations of the maximums
   //Vectors for values above the noise and scaling the noise
   REAL4Vector *ihsvalues = NULL, *excessabovenoise = NULL, *scaledExpectedIHSVectorValues = NULL;
   INT4Vector *ihslocations = NULL;
   XLAL_CHECK( (ihsvalues = XLALCreateREAL4Vector(ihsvectorsequence->length)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (excessabovenoise = XLALCreateREAL4Vector(ihsvectorsequence->vectorLength)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (scaledExpectedIHSVectorValues = XLALCreateREAL4Vector(ihsvectorsequence->vectorLength)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (ihslocations = XLALCreateINT4Vector(ihsvectorsequence->length)) != NULL, XLAL_EFUNC );

   //Finding the maximum for each IHS vector and the location
   for (ii=0; ii<(INT4)ihsvalues->length; ii++) {
      //Use SSE or not
      if (params->useSSE) {
         XLAL_CHECK( sseScaleREAL4Vector(scaledExpectedIHSVectorValues, outputfar->expectedIHSVector, FbinMean->data[ii]) == XLAL_SUCCESS, XLAL_EFUNC ); //Scale the expected IHS vector by the value in FbinMean (in this function, it is 1.0)
         XLAL_CHECK( sseSSVectorSequenceSubtract(excessabovenoise, ihsvectorsequence, scaledExpectedIHSVectorValues, ii) == XLAL_SUCCESS, XLAL_EFUNC );  //subtract the noise from the data
      } else if (params->useAVX) {
         XLAL_CHECK( avxScaleREAL4Vector(scaledExpectedIHSVectorValues, outputfar->expectedIHSVector, FbinMean->data[ii]) == XLAL_SUCCESS, XLAL_EFUNC ); //Scale the expected IHS vector by the value in FbinMean (in this function, it is 1.0)
         XLAL_CHECK( avxSSVectorSequenceSubtract(excessabovenoise, ihsvectorsequence, scaledExpectedIHSVectorValues, ii) == XLAL_SUCCESS, XLAL_EFUNC );  //subtract the noise from the data
      } else {
	 //If not using SSE, do the scaling and subtraction
         for (jj=0; jj<(INT4)scaledExpectedIHSVectorValues->length; jj++) {
            scaledExpectedIHSVectorValues->data[jj] = FbinMean->data[ii]*outputfar->expectedIHSVector->data[jj];
            excessabovenoise->data[jj] = ihsvectorsequence->data[ii*ihsvectorsequence->vectorLength + jj] - scaledExpectedIHSVectorValues->data[jj];
         }
      }
      ihslocations->data[ii] = max_index_in_range(excessabovenoise, minIndexForIHS, maxIndexForIHS) + 5;  //Find the location of the maximum IHS value
      ihsvalues->data[ii] = ihsvectorsequence->data[ii*ihsvectorsequence->vectorLength + ihslocations->data[ii]-5];  //And get the value
      //fprintf(stderr, "%d %f\n", ihslocations->data[ii], ihsvalues->data[ii]);
   } /* for ii=0 --> ihsvalues->length */

   //Some useful variables
   INT4Vector *rowsequencelocs = NULL;
   REAL4Vector *foms = NULL;

   //Starting from a minimum of 2 rows, start determining the FAR for each nearest neighbor sum, up to the maximum number of
   //rows to be summed
   for (ii=2; ii<=rows; ii++) {
      //First allocate the necessary vectors
      XLAL_CHECK( (rowsequencelocs = XLALCreateINT4Vector(ii)) != NULL, XLAL_EFUNC );
      XLAL_CHECK( (foms = XLALCreateREAL4Vector(ihsvectorsequence->length-(ii-1))) != NULL, XLAL_EFUNC );

      //If the user has specified that we should use SSE operations, then do the nearest neighbor summing.
      //The result is in the tworows variable
      if (params->useSSE) {
         if (ii>2) XLAL_CHECK( sseSSVectorSequenceSum(tworows, tworows, ihsvectorsequence, 0, ii-1, 0, (INT4)ihsvectorsequence->length-(ii-1)) == XLAL_SUCCESS, XLAL_EFUNC );
         else XLAL_CHECK( sseSSVectorSequenceSum(tworows, ihsvectorsequence, ihsvectorsequence, 0, ii-1, 0, (INT4)ihsvectorsequence->length-(ii-1)) == XLAL_SUCCESS, XLAL_EFUNC );
      } else if (params->useAVX) {
         if (ii>2) XLAL_CHECK( avxSSVectorSequenceSum(tworows, tworows, ihsvectorsequence, 0, ii-1, 0, (INT4)ihsvectorsequence->length-(ii-1)) == XLAL_SUCCESS, XLAL_EFUNC );
         else XLAL_CHECK( avxSSVectorSequenceSum(tworows, ihsvectorsequence, ihsvectorsequence, 0, ii-1, 0, (INT4)ihsvectorsequence->length-(ii-1)) == XLAL_SUCCESS, XLAL_EFUNC );
      }

      //Now we are going to loop through the input ihsvectorsequence up to the number of rows-1
      for (jj=0; jj<(INT4)ihsvectorsequence->length-(ii-1); jj++) {
         //A bit tricky here: sum IHS values across SFT frequency bins into the tworows variable only if we didn't do it with SSE above
         //We do this with a fast addition function fastSSVectorSequenceSum() and with no error checking, so be careful!
         if (!(params->useSSE || params->useAVX)) {
            if (ii>2) fastSSVectorSequenceSum(tworows, tworows, ihsvectorsequence, jj, ii-1+jj, jj);
            else fastSSVectorSequenceSum(tworows, ihsvectorsequence, ihsvectorsequence, jj, ii-1+jj, jj);
         }

         //Compute IHS FOM value
         memcpy(rowsequencelocs->data, &(ihslocations->data[jj]), sizeof(INT4)*ii);  //First copy the necessary values
         foms->data[jj] = ihsFOM(rowsequencelocs, (INT4)outputfar->expectedIHSVector->length);  //then compute the FOM
         XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
      } /* for jj< ihsvectorsequence->length - (ii-1) */

      //Sample the IHS values that have been summed to compute mean, standard deviation, and FAR threshold values.
      //We have an if-else statement for when there are fewer than 10000 entries that will be in the tworows varaible
      //(only considering the number of rows we have summed together).
      REAL4Vector *sampledtempihsvals = NULL;
      if ((ihsvectorsequence->length-(ii-1))*ihsvectorsequence->vectorLength>10000) {
         //We sample the tworows sequence (up to the number of rows-1) without accepting any zeros
         //because zeros would come from the notched harmonics, which we won't want anyway
         sampledtempihsvals = sampleREAL4VectorSequence_nozerosaccepted(tworows, ihsvectorsequence->length-(ii-1), 10000, params->rng);
         outputfar->ihsdistMean->data[ii-2] = calcMean(sampledtempihsvals);  //And then calculate the mean value
         XLAL_CHECK( calcStddev(&(outputfar->ihsdistSigma->data[ii-2]), sampledtempihsvals) == XLAL_SUCCESS, XLAL_EFUNC );  //We also calculate the standard deviation
      } else {
         //If there were fewer than 10000 entries, then we will keep all those that are part of the nearest neighbor
         //sum up to this point
         XLAL_CHECK( (sampledtempihsvals = XLALCreateREAL4Vector((ihsvectorsequence->length-(ii-1))*ihsvectorsequence->vectorLength)) != NULL, XLAL_EFUNC );
         memcpy(sampledtempihsvals->data, tworows->data, sizeof(REAL4)*sampledtempihsvals->length);
         outputfar->ihsdistMean->data[ii-2] = calcMean_ignoreZeros(sampledtempihsvals);  //Calculate the mean value (remember, don't accept zeros)
         XLAL_CHECK( calcStddev_ignoreZeros(&(outputfar->ihsdistSigma->data[ii-2]), sampledtempihsvals) == XLAL_SUCCESS, XLAL_EFUNC );  //We also calculate the standard deviation
      }

      //If the user has specified the IHS FAR == 1.0, then we don't need to compute the threshold (it is = 0.0)
      INT4 numvals = 0;
      REAL8 farave = 0.0;
      if (params->ihsfar != 1.0) {
         //Looping through the sampled values, we are going to compute the average FAR
         for (jj=0; jj<(INT4)sampledtempihsvals->length; jj++) {
            //When the user has not specified using faster chisq inversion, use the GSL function
            if (!params->fastchisqinv && sampledtempihsvals->data[jj]!=0.0) {
               numvals++;
               //farave += gsl_cdf_chisq_Qinv(params->ihsfar, 0.5*sampledtempihsvals->data[jj]) + 0.5*sampledtempihsvals->data[jj];  //Old
               farave += 0.5*gsl_cdf_chisq_Qinv(params->ihsfar, 2.0*sampledtempihsvals->data[jj]);
            } else if (params->fastchisqinv && sampledtempihsvals->data[jj]!=0.0) {
               numvals++;
               //farave += cdf_chisq_Qinv(params->ihsfar, 0.5*sampledtempihsvals->data[jj]) + 0.5*sampledtempihsvals->data[jj];  //Old
               farave += 0.5*cdf_chisq_Qinv(params->ihsfar, 2.0*sampledtempihsvals->data[jj]);
               XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
            } //fastchisqinv?
         } //for jj=0 --> sampledtempihsval->length
         outputfar->ihsfar->data[ii-2] = farave/(REAL8)numvals;
      } //if params->ihsfar != 1.0
      else {
         outputfar->ihsfar->data[ii-2] = 0.0;
      }

      XLALDestroyREAL4Vector(sampledtempihsvals);

      //FOM part
      outputfar->ihsfomdistMean->data[ii-2] = calcMean(foms);
      XLAL_CHECK( calcStddev(&(outputfar->ihsfomdistSigma->data[ii-2]), foms) == XLAL_SUCCESS, XLAL_EFUNC );
      //We need to find the smallest values
      if (params->ihsfomfar!=1.0 && params->ihsfom==0.0) {
         REAL4Vector *smallestfomvals = NULL;
         XLAL_CHECK( (smallestfomvals = XLALCreateREAL4Vector((INT4)round((ihsvalues->length-ii+1)*params->ihsfomfar)+1)) != NULL, XLAL_EFUNC );
         XLAL_CHECK( sort_float_smallest(smallestfomvals, foms) == XLAL_SUCCESS, XLAL_EFUNC );  //Sort, finding the smallest values
         outputfar->fomfarthresh->data[ii-2] = smallestfomvals->data[smallestfomvals->length-1];  //Pick off the last element
         XLALDestroyREAL4Vector(smallestfomvals);
      } else if (params->ihsfom!=0.0) outputfar->fomfarthresh->data[ii-2] = params->ihsfom;
      else outputfar->fomfarthresh->data[ii-2] = -1.0;

      XLALDestroyINT4Vector(rowsequencelocs);
      rowsequencelocs = NULL;
      XLALDestroyREAL4Vector(foms);
      foms = NULL;
   } /* for ii <= rows */

   XLALDestroyREAL4VectorSequence(tworows);
   XLALDestroyREAL4Vector(ihsvalues);
   XLALDestroyREAL4Vector(excessabovenoise);
   XLALDestroyREAL4Vector(scaledExpectedIHSVectorValues);
   XLALDestroyINT4Vector(ihslocations);

   return XLAL_SUCCESS;

} /*sumIHSSequenceFAR() */


/**
 * \brief Compute the IHS sums for a number of rows
 *
 * In the function we will select the the location which is the maximum above the noise
 *
 * \param [out] output            Pointer to the output ihsMaximaStruct
 * \param [in]  inputfar          Pointer to ihsfarStruct
 * \param [in]  ihsvectorsequence Pointer to REAL4VectorSequence to be summed
 * \param [in]  rows              Number of neighbors to sum
 * \param [in]  FbinMean          Pointer to REAL4Vector of normalized SFT background powers
 * \param [in]  params            Pointer to inputParamsStruct
 * \return Status value
 */
INT4 sumIHSSequence(ihsMaximaStruct *output, ihsfarStruct *inputfar, REAL4VectorSequence *ihsvectorsequence, INT4 rows, REAL4Vector *FbinMean, inputParamsStruct *params)
{

   XLAL_CHECK( output != NULL && inputfar != NULL && ihsvectorsequence != NULL && rows > 0 && FbinMean != NULL && params != NULL, XLAL_EINVAL );

   INT4 ii = 0, jj, kk;

   //Again, we start off by allocating a "tworows" vector sequence of IHS nearest neighbor sums
   REAL4VectorSequence *tworows = NULL;
   XLAL_CHECK( (tworows = XLALCreateREAL4VectorSequence(ihsvectorsequence->length-1, ihsvectorsequence->vectorLength)) != NULL, XLAL_EFUNC );
   memset(tworows->data, 0, sizeof(REAL4)*tworows->length*tworows->vectorLength);      //Set everything to 0.0

   //Allocation of ihs values and locations
   //Vectors for values above the noise and scaling the noise
   REAL4Vector *ihsvalues = NULL, *excessabovenoise = NULL, *scaledExpectedIHSVectorValues = NULL;
   INT4Vector *ihslocations = NULL;
   XLAL_CHECK( (ihsvalues = XLALCreateREAL4Vector(ihsvectorsequence->length)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (excessabovenoise = XLALCreateREAL4Vector(ihsvectorsequence->vectorLength)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (scaledExpectedIHSVectorValues = XLALCreateREAL4Vector(ihsvectorsequence->vectorLength)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (ihslocations = XLALCreateINT4Vector(ihsvectorsequence->length)) != NULL, XLAL_EFUNC );

   //The minimum and maximum index to search in the IHS vector
   INT4 maxIndexForIHS = (INT4)ceil(fmin(params->Tobs/params->Pmin, fmin(params->Tobs/minPeriod(params->dfmin, params->Tcoh), params->Tobs/(4.0*params->Tcoh)))) - 5;
   INT4 minIndexForIHS = (INT4)floor(fmax(5.0, params->Tobs/params->Pmax)) - 5;

   //Finding the maximum for each IHS vector and the location using SSE or not
   for (ii=0; ii<(INT4)ihsvalues->length; ii++) {
      if (params->useSSE) {
         XLAL_CHECK( sseScaleREAL4Vector(scaledExpectedIHSVectorValues, inputfar->expectedIHSVector, FbinMean->data[ii]) == XLAL_SUCCESS, XLAL_EFUNC );  //Scale the expected IHS vector by the FbinMean data value
         XLAL_CHECK( sseSSVectorSequenceSubtract(excessabovenoise, ihsvectorsequence, scaledExpectedIHSVectorValues, ii) == XLAL_SUCCESS, XLAL_EFUNC );  //subtract the noise from the data
      } else if (params->useAVX) {
         XLAL_CHECK( avxScaleREAL4Vector(scaledExpectedIHSVectorValues, inputfar->expectedIHSVector, FbinMean->data[ii]) == XLAL_SUCCESS, XLAL_EFUNC );  //Scale the expected IHS vector by the FbinMean data value
         XLAL_CHECK( avxSSVectorSequenceSubtract(excessabovenoise, ihsvectorsequence, scaledExpectedIHSVectorValues, ii) == XLAL_SUCCESS, XLAL_EFUNC );  //subtract the noise from the data
      } else {
         for (jj=0; jj<(INT4)scaledExpectedIHSVectorValues->length; jj++) {
            scaledExpectedIHSVectorValues->data[jj] = FbinMean->data[ii]*inputfar->expectedIHSVector->data[jj];
            excessabovenoise->data[jj] = ihsvectorsequence->data[ii*ihsvectorsequence->vectorLength + jj] - scaledExpectedIHSVectorValues->data[jj];
         }
      }

      //search over the range of Pmin-->Pmax and higher harmonics the user has specified
      for (jj=0; jj<params->harmonicNumToSearch; jj++) {
         if (jj==0) {
            ihslocations->data[ii] = max_index_in_range(excessabovenoise, minIndexForIHS, maxIndexForIHS) + 5;
            ihsvalues->data[ii] = ihsvectorsequence->data[ii*ihsvectorsequence->vectorLength + ihslocations->data[ii]-5];
         } else {
            INT4 newIHSlocation = max_index_in_range(excessabovenoise, (jj+1)*minIndexForIHS, (jj+1)*maxIndexForIHS) + 5;
            REAL4 newIHSvalue = ihsvectorsequence->data[ii*ihsvectorsequence->vectorLength + newIHSlocation-5];
            if (newIHSvalue > ihsvalues->data[ii]) {
               ihslocations->data[ii] = newIHSlocation;
               ihsvalues->data[ii] = newIHSvalue;
            } /* if the new value is better than the previous value */
         }
      } /* for jj=0 --> jj<harmonicNumToSearch */
   }

   //Useful variables
   INT4Vector *rowsequencelocs = NULL;

   //Start with the single IHS vector and march up with nearest neighbor sums up to the total number of row sums
   for (ii=1; ii<=rows; ii++) {
      if (ii==1) {
         //Copy the data into the output
         memcpy(output->maximaForEachFbin->data, ihsvalues->data, sizeof(REAL4)*ihsvalues->length);
         memcpy(output->locationsForEachFbin->data, ihslocations->data, sizeof(INT4)*ihslocations->length);
      } else {
         //For everything 2 nearest neighbors and higher summed
         XLAL_CHECK( (rowsequencelocs = XLALCreateINT4Vector(ii)) != NULL, XLAL_EFUNC );

         //The maximum index to search in the IHS vector
         maxIndexForIHS = (INT4)ceil(fmin(params->Tobs/params->Pmin, fmin(params->Tobs/minPeriod(0.5*(ii-1)/params->Tcoh, params->Tcoh), params->Tobs/(4.0*params->Tcoh)))) - 5;

         REAL4 sumofnoise = 0.0;    //To scale the expected IHS background
         INT4 endloc = ((ii-1)*(ii-1)-(ii-1))/2;

         //Sum up the IHS vectors using SSE functions
         if (params->useSSE) {
            if (ii>2) XLAL_CHECK( sseSSVectorSequenceSum(tworows, tworows, ihsvectorsequence, 0, ii-1, 0, (INT4)ihsvectorsequence->length-(ii-1)) == XLAL_SUCCESS, XLAL_EFUNC );
            else XLAL_CHECK( sseSSVectorSequenceSum(tworows, ihsvectorsequence, ihsvectorsequence, 0, ii-1, 0, (INT4)ihsvectorsequence->length-(ii-1)) == XLAL_SUCCESS, XLAL_EFUNC );
         } else if (params->useAVX) {
            if (ii>2) XLAL_CHECK( avxSSVectorSequenceSum(tworows, tworows, ihsvectorsequence, 0, ii-1, 0, (INT4)ihsvectorsequence->length-(ii-1)) == XLAL_SUCCESS, XLAL_EFUNC );
            else XLAL_CHECK( avxSSVectorSequenceSum(tworows, ihsvectorsequence, ihsvectorsequence, 0, ii-1, 0, (INT4)ihsvectorsequence->length-(ii-1)) == XLAL_SUCCESS, XLAL_EFUNC );
         } /* use SSE or AVX code */

         //Loop through the IHS vector neighbor sums
         for (jj=0; jj<(INT4)ihsvectorsequence->length-(ii-1); jj++) {
            //If we didn't use SSE to sum the vector sequence (see lines above)
            if (!(params->useSSE || params->useAVX)) {
               if (ii>2) fastSSVectorSequenceSum(tworows, tworows, ihsvectorsequence, jj, ii-1+jj, jj);
               else fastSSVectorSequenceSum(tworows, ihsvectorsequence, ihsvectorsequence, jj, jj+1, jj);
            }

            //To scale the background efficiently
            if (jj==0) for (kk=0; kk<ii; kk++) sumofnoise += FbinMean->data[kk];
            else {
               sumofnoise -= FbinMean->data[jj-1];
               sumofnoise += FbinMean->data[jj+(ii-1)];
            }

            //If using SSE, scale the expected IHS vector, subtract the noise from the data
            if (params->useSSE) {
               XLAL_CHECK( sseScaleREAL4Vector(scaledExpectedIHSVectorValues, inputfar->expectedIHSVector, sumofnoise) == XLAL_SUCCESS, XLAL_EFUNC );
               XLAL_CHECK( sseSSVectorSequenceSubtract(excessabovenoise, tworows, scaledExpectedIHSVectorValues, jj) == XLAL_SUCCESS, XLAL_EFUNC );
            } else {
               REAL4 scaleval = sumofnoise;
               for (kk=0; kk<(INT4)inputfar->expectedIHSVector->length; kk++) scaledExpectedIHSVectorValues->data[kk] = scaleval*inputfar->expectedIHSVector->data[kk];
               fastSSVectorSequenceSubtract(excessabovenoise, tworows, scaledExpectedIHSVectorValues, jj);
            }

            //Compute the maximum IHS value in the second FFT frequency direction
            //search over the range of Pmin-->Pmax and higher harmonics the user has specified
            for (kk=0; kk<params->harmonicNumToSearch; kk++) {
               if (kk==0) {
                  output->locations->data[(ii-2)*ihsvalues->length-endloc+jj] = max_index_in_range(excessabovenoise, minIndexForIHS, maxIndexForIHS) + 5;
                  output->maxima->data[(ii-2)*ihsvalues->length-endloc+jj] = tworows->data[jj*tworows->vectorLength + (output->locations->data[(ii-2)*ihsvalues->length-endloc+jj]-5)];
               } else {
                  INT4 newIHSlocation = max_index_in_range(excessabovenoise, (kk+1)*minIndexForIHS, (kk+1)*maxIndexForIHS) + 5;
                  REAL4 newIHSvalue = tworows->data[ii*tworows->vectorLength + newIHSlocation-5];
                  if (newIHSvalue > output->maxima->data[(ii-2)*ihsvalues->length-endloc+jj]) {
                     output->locations->data[(ii-2)*ihsvalues->length-endloc+jj] = newIHSlocation;
                     output->maxima->data[(ii-2)*ihsvalues->length-endloc+jj] = newIHSvalue;
                  } /* if the new value is better than the previous value */
               }
            } /* for kk=0 --> kk<harmonicNumToSearch */

            memcpy(rowsequencelocs->data, &(ihslocations->data[jj]), sizeof(INT4)*ii);
            output->foms->data[(ii-2)*ihsvalues->length-endloc+jj] = ihsFOM(rowsequencelocs, (INT4)inputfar->expectedIHSVector->length);
            XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
         } /* for jj< ihsvectorsequence->length - (ii-1) */

         XLALDestroyINT4Vector(rowsequencelocs);
         rowsequencelocs = NULL;
      }

   } /* for ii <= rows */

   XLALDestroyREAL4VectorSequence(tworows);
   XLALDestroyREAL4Vector(scaledExpectedIHSVectorValues);
   XLALDestroyREAL4Vector(excessabovenoise);
   XLALDestroyREAL4Vector(ihsvalues);
   XLALDestroyINT4Vector(ihslocations);

   return XLAL_SUCCESS;

} /*sumIHSSequence() */


/**
 * Calculate the IHS FOM for a number of rows
 * \param [in] locs    Pointer to INT4Vector of location values
 * \param [in] fomnorm Normalization value and is the number of XXX
 * \return IHS FOM value
 */
REAL4 ihsFOM(INT4Vector *locs, INT4 fomnorm)
{

   XLAL_CHECK_REAL4( locs != NULL && fomnorm > 0.0, XLAL_EINVAL );

   INT4 ii;
   REAL4 fom = 0.0;

   for (ii=0; ii<(INT4)(locs->length*0.5); ii++) fom += (REAL4)((locs->data[ii]-locs->data[locs->length-ii-1])*(locs->data[ii]-locs->data[locs->length-ii-1]))/(fomnorm*fomnorm);

   return fom;

} /* ihsFOM() */


/**
 * Find IHS candidates above thresholds
 * \param [out] candlist     Pointer to a pointer containing the candidate list
 * \param [in]  ihsfarstruct Pointer to ihsfarStruct
 * \param [in]  params       Pointer to inputParamsStruct
 * \param [in]  ffdata       Pointer to ffdataStruct
 * \param [in]  ihsmaxima    Pointer to ihsMaximaStruct containing the data to be tested above thresholds
 * \param [in]  fbinavgs     Pointer to REAL4Vector of normalized SFT background powers
 * \param [in]  trackedlines Pointer to REAL4VectorSequence of lines (allowed to be NULL)
 * \return Status value
 */
INT4 findIHScandidates(candidateVector **candlist, ihsfarStruct *ihsfarstruct, inputParamsStruct *params, ffdataStruct *ffdata, ihsMaximaStruct *ihsmaxima, REAL4Vector *fbinavgs, REAL4VectorSequence *trackedlines)
{

   XLAL_CHECK( ihsfarstruct != NULL && params != NULL && ffdata != NULL && ihsmaxima != NULL && fbinavgs != NULL, XLAL_EINVAL );

   INT4 ii, jj, kk;
   REAL8 fsig, per0, B;

   INT4 numberofIHSvalsChecked = 0, numberofIHSvalsExceededThresh = 0, numberPassingBoth = 0, linesinterferewithnum = 0, skipped = 0, notskipped = 0;

   INT4 numfbins = ffdata->numfbins;  //number of fbins
   INT4 minrows = (INT4)round(2.0*params->dfmin*params->Tcoh)+1;  //minimum number of sequential rows to search

   REAL4Vector *ihss = NULL, *avgsinrange = NULL;
   INT4Vector *locs = NULL;

   //Check the IHS values against the FAR, checking between IHS width values
   //FILE *IHSVALSOUTPUT = fopen("./output/allihsvalspassthresh.dat","w");
   for (ii=minrows; ii<=(INT4)ihsfarstruct->ihsfar->length+1; ii++) {
      XLAL_CHECK( (ihss = XLALCreateREAL4Vector(ii)) != NULL, XLAL_EFUNC );
      XLAL_CHECK( (locs = XLALCreateINT4Vector(ii)) != NULL, XLAL_EFUNC );
      XLAL_CHECK( (avgsinrange = XLALCreateREAL4Vector(ii)) != NULL, XLAL_EFUNC );

      REAL8 highestval = 0.0, highestsignificance = 0.0; //highestvalnoise = 0.0
      INT4 highestvalloc = -1, jjloc = 0;
      for (jj=0; jj<(INT4)numfbins-(ii-1); jj++) {

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
               fsig = params->fmin - params->dfmax + ((0.5*(ii-1) + jj) - 6.0)/params->Tcoh;    //Candidate frequency
               B = 0.5*(ii-1)/params->Tcoh;                                                     //Candidate modulation depth

               //Test to see if any tracked lines are overlapping the candidate signal
               INT4 nolinesinterfering = 1;
               if (trackedlines!=NULL) {
                  kk = 0;
                  while (kk<(INT4)trackedlines->length && nolinesinterfering==1) {
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

                  if ( significance > highestsignificance && (params->followUpOutsideULrange || (!params->followUpOutsideULrange && fsig>=params->ULfmin && fsig<(params->ULfmin+params->ULfspan) && B>=params->ULmindf && B<=params->ULmaxdf)) ) {
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
         fsig = params->fmin - params->dfmax + (0.5*(ii-1) + jjloc - 6.0)/params->Tcoh;  //Candidate frequency
         B = 0.5*(ii-1)/params->Tcoh;                                                    //Candidate modulation depth
         per0 = params->Tobs/loc;                                                        //Candidate period
         REAL8 h0 = ihs2h0(2.0*highestval, params);  //Candidate h0, need factor of 2 for the degrees of freedom counting

         if ((*candlist)->numofcandidates == (*candlist)->length-1) XLAL_CHECK( (*candlist = resize_candidateVector(*candlist, 2*(*candlist)->length)) != NULL, XLAL_EFUNC );
         loadCandidateData(&((*candlist)->data[(*candlist)->numofcandidates]), fsig, per0, B, 0.0, 0.0, ihsmaxima->maxima->data[highestvalloc], h0, highestsignificance, 0, ffdata->tfnormalization);
         (*candlist)->numofcandidates++;
      }

      //Destroy
      XLALDestroyREAL4Vector(ihss);
      XLALDestroyREAL4Vector(avgsinrange);
      XLALDestroyINT4Vector(locs);
      ihss = NULL;
      locs = NULL;
      avgsinrange = NULL;

   } /* for ii < ihsfarstruct->ihsfar->length */

   //The outer loop is over "frequency" (the order in the ihs maxima vector, first 2 row sums, then three, and so on)
   //for (ii=0; ii<(INT4)numfbins; ii++) {
   /* for (ii=6+(INT4)round(params->dfmax*params->Tcoh); ii<(numfbins-6-(INT4)round(params->dfmax*params->Tcoh)); ii++) {
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
               fsig = params->fmin - params->dfmax + (0.5*(jj-1) + ii - 6)/params->Tcoh;             //Candidate frequency
               B = 0.5*(jj-1)/params->Tcoh;                                      //Candidate modulation depth

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
         fsig = params->fmin - params->dfmax + (0.5*(jjloc-1) + ii - 6.0)/params->Tcoh;
         //Candidate modulation depth
         B = 0.5*(jjloc-1)/params->Tcoh;
         //Candidate period
         per0 = params->Tobs/loc;
         //Candidate h0
         REAL8 h0 = ihs2h0(2.0*highestval, params);  //Need factor of 2 for the degrees of freedom counting
         REAL8 significance = highestsignificance;
         //fprintf(stderr, "%d %d %f\n", ii, jjloc, significance);     //remove this

         if (candlist->numofcandidates == candlist->length-1) {
            candlist = resize_candidateVector(candlist, 2*(candlist->length));
            if (candlist->data==NULL) {
               fprintf(stderr,"%s: resize_candidateVector() failed.\n", __func__);
               XLAL_ERROR_VOID(XLAL_EFUNC);
            }
         }
         //loadCandidateData(&candlist->data[candlist->numofcandidates], fsig, per0, B, 0.0, 0.0, ihsmaxima->maxima->data[locationinmaximastruct], h0, 0.0, 0, sqrt(ffdata->tfnormalization/2.0*params->Tcoh));
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
 * \param [in] params Pointer to inputParamsStruct
 * \return Estimated h0
 */
REAL8 ihs2h0(REAL8 ihsval, inputParamsStruct *params)
{

   if (ihsval<=0.0) return 0.0;
   REAL8 prefact = 1.0;
   //prefact = 7.2;  //old value for when IHS --> h0 when in the upper limits, the value of the noise subtracted was incorrectly calculated and for >99.999% CL
   prefact = 44.7; //new value for IHS --> h0
   return prefact*pow(ihsval/(params->Tcoh*params->Tobs),0.25);

}
