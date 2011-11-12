/*
*  Copyright (C) 2010, 2011 Evan Goetz
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

#ifdef __SSE__
#include <xmmintrin.h>
#endif

#include <lal/LALMalloc.h>
#include <lal/SeqFactories.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "IHS.h"
#include "statistics.h"
#include "candidates.h"
#include "fastchisqinv.h"
#include "TwoSpect.h"


//////////////////////////////////////////////////////////////
// Create vectors for IHS maxima struct
ihsMaximaStruct * new_ihsMaxima(INT4 fbins, INT4 rows)
{
   
   ihsMaximaStruct *ihsmaxima = XLALMalloc(sizeof(*ihsmaxima));
   if (ihsmaxima==NULL) {
      fprintf(stderr,"%s: XLALMalloc(%zu) failed.\n", __func__, sizeof(*ihsmaxima));
      XLAL_ERROR_NULL(XLAL_ENOMEM);
   }
   
   INT4 numberofmaxima = fbins*rows - fbins - (INT4)((rows*rows-rows)/2);
   
   ihsmaxima->maxima = XLALCreateREAL4Vector(numberofmaxima);
   ihsmaxima->locations = XLALCreateINT4Vector(numberofmaxima);
   ihsmaxima->foms = XLALCreateREAL4Vector(numberofmaxima);
   ihsmaxima->maximaForEachFbin = XLALCreateREAL4Vector(fbins);
   ihsmaxima->locationsForEachFbin = XLALCreateINT4Vector(fbins);
   ihsmaxima->rows = rows;
   
   //Fail if any of the allocations return NULL pointers
   if (ihsmaxima->maxima==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, numberofmaxima);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   } else if (ihsmaxima->locations==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", __func__, numberofmaxima);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   } else if (ihsmaxima->foms==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, numberofmaxima);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   } else if (ihsmaxima->maximaForEachFbin==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, fbins);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   } else if (ihsmaxima->locationsForEachFbin==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", __func__, fbins);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
      
   return ihsmaxima;

} /* new_ihsMaxima() */


//////////////////////////////////////////////////////////////
// Destroy vectors for IHS maxima struct
void free_ihsMaxima(ihsMaximaStruct *data)
{

   XLALDestroyREAL4Vector(data->maxima);
   XLALDestroyINT4Vector(data->locations);
   XLALDestroyREAL4Vector(data->foms);
   XLALDestroyREAL4Vector(data->maximaForEachFbin);
   XLALDestroyINT4Vector(data->locationsForEachFbin);
   XLALFree((ihsMaximaStruct*)data);

} /* free_ihsMaxima() */


//////////////////////////////////////////////////////////////
// Run the IHS algorithm
void runIHS(ihsMaximaStruct *output, ffdataStruct *input, ihsfarStruct *ihsfarinput, inputParamsStruct *params, INT4 rows, REAL4Vector *FbinMean)
{
   
   INT4 ii, jj;
   
   INT4 numfbins = input->numfbins;
   INT4 numfprbins = input->numfprbins;
   
   REAL4Vector *row = XLALCreateREAL4Vector(numfprbins);
   REAL4Vector *ihss = XLALCreateREAL4Vector(numfbins);
   REAL4Vector *ihsvector = XLALCreateREAL4Vector((INT4)floor((1.0/(REAL8)params->ihsfactor)*numfprbins)-5);
   INT4Vector *locs = XLALCreateINT4Vector(numfbins);
   ihsVals *ihsvals = new_ihsVals();
   if (row==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, numfprbins);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (ihss==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, numfbins);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (ihsvector==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, (INT4)floor((1.0/(REAL8)params->ihsfactor)*numfprbins)-5);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (locs==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", __func__, numfbins);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (ihsvals==NULL) {
      fprintf(stderr,"%s: new_ihsVals() failed.\n", __func__);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   
   REAL4VectorSequence *ihsvectorsequence = XLALCreateREAL4VectorSequence(numfbins, ihsvector->length);
   if (ihsvectorsequence==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4VectorSequence(%d,%d) failed.\n", __func__, numfbins, ihsvector->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   
   //Loop through the rows, 1 frequency at a time
   REAL8 dailyharmonic = params->Tobs/(24.0*3600.0);
   REAL8 dailyharmonic2 = dailyharmonic*2.0, dailyharmonic3 = dailyharmonic*3.0, dailyharmonic4 = dailyharmonic*4.0;
   for (ii=0; ii<(INT4)ihss->length; ii++) {
   
      //For each row, populate it with the data for that frequency bin, excluding harmonics of antenna pattern modulation
      memcpy(row->data, &(input->ffdata->data[ii*numfprbins]), sizeof(REAL4)*numfprbins);
      for (jj=0; jj<(INT4)row->length; jj++) if (fabs(dailyharmonic-(REAL8)jj)<=1.0 || fabs(dailyharmonic2-(REAL8)jj)<=1.0 || fabs(dailyharmonic3-(REAL8)jj)<=1.0 || fabs(dailyharmonic4-(REAL8)jj)<=1.0) row->data[jj] = 0.0;
      
      //Run the IHS algorithm on the row
      incHarmSumVector(ihsvector, row, params->ihsfactor);
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: incHarmSumVector() failed.\n", __func__);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      }
      //locs->data[ii] = max_index(ihsvector) + 5;
      //ihss->data[ii] = ihsvector->data[locs->data[ii]-5];
      memcpy(&(ihsvectorsequence->data[ii*ihsvector->length]), ihsvector->data, sizeof(REAL4)*ihsvector->length);
      
   } /* for ii < ihss->length */
   
   //ihsSums2_withFAR(output, NULL, ihsvectorsequence, ihss, locs, aveNoise, rows, FbinMean, input->numfprbins, params, 0);
   //ihsSums2_withFAR_withnoise(output, NULL, ihsvectorsequence, ihss, locs, aveNoise, rows, FbinMean, input->numfprbins, params, 0);
   sumIHSSequence(output, ihsfarinput, ihsvectorsequence, rows, FbinMean, input->numfprbins, params);
   XLALDestroyREAL4VectorSequence(ihsvectorsequence);
   
   
   //Destroy variables
   XLALDestroyREAL4Vector(row);
   XLALDestroyREAL4Vector(ihss);
   XLALDestroyREAL4Vector(ihsvector);
   XLALDestroyINT4Vector(locs);
   free_ihsVals(ihsvals);

} /* runIHS() */


//////////////////////////////////////////////////////////////
// Allocate memory for ihsVals struct
ihsVals * new_ihsVals(void)
{
   
   ihsVals *ihsvals = XLALMalloc(sizeof(*ihsvals));
   if (ihsvals==NULL) {
      fprintf(stderr,"%s: XLALMalloc(%zu) failed.\n", __func__, sizeof(*ihsvals));
      XLAL_ERROR_NULL(XLAL_ENOMEM);
   }

   return ihsvals;

} /* new_ihsVals() */

//////////////////////////////////////////////////////////////
// Destroy ihsVals struct
void free_ihsVals(ihsVals *ihsvals)
{

   XLALFree((ihsVals*)ihsvals);

} /* free_ihsVals() */


//////////////////////////////////////////////////////////////
// Compute the IHS sum
void incHarmSum(ihsVals *output, REAL4Vector *input, INT4 ihsfactor)
{
   
   INT4 ii;
   
   output->ihs = 0.0;
   
   //Start ii >= 15
   /*for (ii=15; ii<(INT4)input->length; ii++) {
      REAL4 sum = input->data[ii] + input->data[(INT4)(ii*0.5)] + input->data[(INT4)(ii/3.0)] + input->data[(INT4)(ii*0.25)] + input->data[(INT4)(ii*0.2)];

      if (sum > output->ihs) {
         output->ihs = sum;
         output->loc = (INT4)round(ii/3.0);
      }
   } */ /* for ii < input->length */
   /* for (ii=5; ii<(INT4)(0.2*input->length); ii++) {
      REAL4 sum = input->data[ii] + input->data[(INT4)(ii*2.0)] + input->data[(INT4)(ii*3.0)] + input->data[(INT4)(ii*4.0)] + input->data[(INT4)(ii*5.0)];
      if (sum > output->ihs) {
         output->ihs = sum;
         output->loc = ii;
      }
   } */ /* for ii = 5 --> 0.2*input->length */
   
   //FILE *ALLIHSVALS = fopen("./output/ihsvect.dat","a");
   REAL4Vector *tempvect = XLALCreateREAL4Vector((INT4)floor((1.0/(REAL8)ihsfactor)*input->length)-5);
   if (tempvect==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, (INT4)floor((1.0/(REAL8)ihsfactor)*input->length)-5);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   incHarmSumVector(tempvect, input, ihsfactor);
   if (xlalErrno!=0) {
      fprintf(stderr, "%s: incHarmSumVector() failed.\n", __func__);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   for (ii=0; ii<(INT4)tempvect->length; ii++) {
      //fprintf(ALLIHSVALS,"%.8g\n",tempvect->data[ii]);
      if (tempvect->data[ii]>output->ihs) {
         output->ihs = tempvect->data[ii];
         output->loc = ii+5;
      }
   } /* for ii < tempvect->length */
   XLALDestroyREAL4Vector(tempvect);
   //fclose(ALLIHSVALS);

} /* incHarmSum() */


//////////////////////////////////////////////////////////////
// Compute the IHS vector
void incHarmSumVector(REAL4Vector *output, REAL4Vector *input, INT4 ihsfactor)
{
   
   INT4 ii, jj, highval = (INT4)floor((1.0/(REAL8)ihsfactor)*input->length);
   
   for (ii=5; ii<highval; ii++) {
      output->data[ii-5] = 0.0;
      for (jj=1; jj<=ihsfactor; jj++) output->data[ii-5] += input->data[ii*jj];
      if (ihsfactor*(ii-1)>(INT4)input->length-1) {
         fprintf(stderr, "%s: final point exceeds the allowed limits!\n", __func__);
         XLAL_ERROR_VOID(XLAL_EBADLEN);
      }
   }
   
} /* incHarmSumVector() */


//////////////////////////////////////////////////////////////
// Allocate memory for ihsfarStruct struct
ihsfarStruct * new_ihsfarStruct(INT4 rows, inputParamsStruct *params)
{
   
   ihsfarStruct *ihsfarstruct = XLALMalloc(sizeof(*ihsfarstruct));
   if (ihsfarstruct == NULL) {
      fprintf(stderr,"%s: XLALMalloc(%zu) failed.\n", __func__, sizeof(*ihsfarstruct));
      XLAL_ERROR_NULL(XLAL_ENOMEM);
   }
   
   ihsfarstruct->ihsfar = XLALCreateREAL4Vector(rows-1);
   ihsfarstruct->ihsdistMean = XLALCreateREAL4Vector(rows-1);
   ihsfarstruct->ihsdistSigma = XLALCreateREAL4Vector(rows-1);
   ihsfarstruct->fomfarthresh = XLALCreateREAL4Vector(rows-1);
   ihsfarstruct->ihsfomdistMean = XLALCreateREAL4Vector(rows-1);
   ihsfarstruct->ihsfomdistSigma = XLALCreateREAL4Vector(rows-1);
   ihsfarstruct->expectedIHSVector = XLALCreateREAL4Vector((INT4)floor((1.0/(REAL8)params->ihsfactor)*((INT4)floor(floor(params->Tobs/(params->Tcoh-params->SFToverlap)-1)*0.5)+1))-5);
   if (ihsfarstruct->ihsfar==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, rows-1);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   } else if (ihsfarstruct->ihsdistMean==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, rows-1);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   } else if( ihsfarstruct->ihsdistSigma==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, rows-1);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   } else if( ihsfarstruct->fomfarthresh==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, rows-1);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   } else if( ihsfarstruct->ihsfomdistMean==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, rows-1);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   } else if( ihsfarstruct->ihsfomdistSigma==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, rows-1);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   } else if (ihsfarstruct->expectedIHSVector==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, (INT4)floor((1.0/(REAL8)params->ihsfactor)*((INT4)floor(floor(params->Tobs/(params->Tcoh-params->SFToverlap)-1)*0.5)+1))-5);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
   memset(ihsfarstruct->expectedIHSVector->data, 0, sizeof(REAL4)*ihsfarstruct->expectedIHSVector->length);
   
   return ihsfarstruct;

} /* new_ihsfarStruct() */


//////////////////////////////////////////////////////////////
// Destroy ihsfarStruct struct
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


//////////////////////////////////////////////////////////////
// Compute the IHS FAR for a sum of a number of rows
void genIhsFar(ihsfarStruct *output, inputParamsStruct *params, INT4 rows, REAL4Vector *aveNoise)
{
   
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
   
   /* INT4 trials = (INT4)round(1.0e-11/params->ihsfar);    //Number of trials to determine FAR value
   if (params->ihsfomfar!=0.0 && trials<(INT4)round(1.0e-11/params->ihsfomfar)) {
      trials = (INT4)round(1.0e-11/params->ihsfomfar);
   }
   trials += rows; */
   
   //Initialize random number generator
   gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
   if (rng==NULL) {
      fprintf(stderr,"%s: gsl_rng_alloc() failed.\n", __func__);
      XLAL_ERROR_VOID(XLAL_ENOMEM);
   }
   srand(time(NULL));
   UINT8 randseed = rand();
   gsl_rng_set(rng, randseed);
   //gsl_rng_set(rng, 0);
   
   //Allocations for IHS values for the number of trials
   REAL4Vector *noise = XLALCreateREAL4Vector(aveNoise->length);
   REAL4Vector *ihsvector = XLALCreateREAL4Vector((INT4)floor((1.0/(REAL8)params->ihsfactor)*aveNoise->length)-5);
   REAL4VectorSequence *ihsvectorsequence = XLALCreateREAL4VectorSequence(trials, ihsvector->length);
   REAL4Vector *ihss = XLALCreateREAL4Vector(trials);
   INT4Vector *locs = XLALCreateINT4Vector(trials);
   if (noise==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, aveNoise->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (ihsvector==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, (INT4)floor((1.0/(REAL8)params->ihsfactor)*aveNoise->length)-5);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (ihsvectorsequence==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4VectorSequence(%d,%d) failed.\n", __func__, trials, ihsvector->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (ihss==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, trials);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (locs==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", __func__, trials);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   REAL8 singleIHSsigma = sqrt(5.0*0.05*0.05);
   //REAL8 sigmaval = 0.05;
   REAL8 dailyharmonic = Tobs/(24.0*3600.0);
   REAL8 dailyharmonic2 = dailyharmonic*2.0, dailyharmonic3 = dailyharmonic*3.0, dailyharmonic4 = dailyharmonic*4.0;
   INT4Vector *markedharmonics = XLALCreateINT4Vector(aveNoise->length);
   if (markedharmonics==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", __func__, aveNoise->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   memset(markedharmonics->data, 0, sizeof(INT4)*markedharmonics->length);
   for (ii=0; ii<(INT4)markedharmonics->length; ii++) {
      if (fabs(dailyharmonic-(REAL8)ii)<=1.0 || fabs(dailyharmonic2-(REAL8)ii)<=1.0 || fabs(dailyharmonic3-(REAL8)ii)<=1.0 || fabs(dailyharmonic4-(REAL8)ii)<=1.0) markedharmonics->data[ii] = 1;
   }
   for (ii=0; ii<trials; ii++) {
      REAL8 randval = 1.0 + 2.0*gsl_ran_gaussian(rng, singleIHSsigma);
      while (randval<0.0) randval = 1.0 + 2.0*gsl_ran_gaussian(rng, singleIHSsigma);
      
      //Make exponential noise removing harmonics of 24 hours to match with the same method as real analysis
      for (jj=0; jj<(INT4)aveNoise->length; jj++) {
         if (markedharmonics->data[jj]==0) noise->data[jj] = (REAL4)(gsl_ran_exponential(rng, aveNoise->data[jj])*randval);
         
         //REAL8 randval = 1.0 + 6.0*gsl_ran_gaussian(rng, sigmaval);
         //while (randval<0.0) randval = 1.0 + 6.0*gsl_ran_gaussian(rng, sigmaval);
         //if (markedharmonics->data[jj]==0) noise->data[jj] = (REAL4)(gsl_ran_exponential(rng, aveNoise->data[jj]*randval));
         
         else noise->data[jj] = 0.0;
      } /* for jj < aveNoise->length */
      
      //Compute IHS value on exponential noise
      incHarmSumVector(ihsvector, noise, params->ihsfactor);
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: incHarmSumVector() failed.\n", __func__);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      }
      
      memcpy(&(ihsvectorsequence->data[ii*ihsvector->length]), ihsvector->data, sizeof(REAL4)*ihsvector->length);
   } /* for ii < trials */
   XLALDestroyREAL4Vector(noise);
   XLALDestroyREAL4Vector(ihsvector);
   XLALDestroyINT4Vector(markedharmonics);
   gsl_rng_free(rng);
   
   
   //Calculate the IHS sum values for the IHS trials
   ihsMaximaStruct *ihsmaxima = new_ihsMaxima(trials, rows);
   if (ihsmaxima==NULL) {
      fprintf(stderr,"%s: new_ihsMaxima(%d, %d) failed.\n", __func__, trials, rows);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   //Create a fake vector with the same average value in each bin
   REAL4Vector *FbinMean = XLALCreateREAL4Vector(trials);
   if (FbinMean==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, trials);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   for (ii=0; ii<trials; ii++) FbinMean->data[ii] = 1.0;
   
   sumIHSSequenceFAR(output, ihsvectorsequence, rows, FbinMean, (INT4)aveNoise->length, params);
   XLALDestroyREAL4VectorSequence(ihsvectorsequence);
   
   //Destroy variables
   XLALDestroyREAL4Vector(ihss);
   XLALDestroyREAL4Vector(FbinMean);
   XLALDestroyINT4Vector(locs);
   free_ihsMaxima(ihsmaxima);
   

} /* genIhsFar() */



//////////////////////////////////////////////////////////////
// Compute the IHS sums for a number of rows
void sumIHSSequenceFAR(ihsfarStruct *outputfar, REAL4VectorSequence *ihsvectorsequence, INT4 rows, REAL4Vector *FbinMean, INT4 locationnormfactor, inputParamsStruct *params)
{
   
   INT4 ii, jj;
   
   REAL4VectorSequence *tworows = XLALCreateREAL4VectorSequence(ihsvectorsequence->length-1, ihsvectorsequence->vectorLength);
   if (tworows==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4VectorSequence(%d,%d) failed.\n", __func__, ihsvectorsequence->length-1, ihsvectorsequence->vectorLength);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   memset(tworows->data, 0, sizeof(REAL4)*tworows->length*tworows->vectorLength);
   
   REAL4Vector *ihsvalues = XLALCreateREAL4Vector(ihsvectorsequence->length);
   INT4Vector *ihslocations = XLALCreateINT4Vector(ihsvalues->length);
   if (ihsvalues==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ihsvectorsequence->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (ihslocations==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ihsvalues->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   for (ii=0; ii<(INT4)ihsvalues->length; ii++) {
      ihslocations->data[ii] = max_index_from_vector_in_REAL4VectorSequence(ihsvectorsequence, ii) + 5;
      ihsvalues->data[ii] = ihsvectorsequence->data[ii*ihsvectorsequence->vectorLength + ihslocations->data[ii]-5];
   }
   
   //get average
   for (ii=0; ii<(INT4)ihsvectorsequence->vectorLength; ii++) {
      for (jj=0; jj<(INT4)ihsvectorsequence->length; jj++) {
         outputfar->expectedIHSVector->data[ii] += ihsvectorsequence->data[jj*tworows->vectorLength + ii];
      }
      outputfar->expectedIHSVector->data[ii] /= (REAL4)jj;
   }   
   
   REAL4Vector *rowsequencemaxima = NULL;
   REAL4Vector *fbinmeanvals = NULL;
   INT4Vector *rowsequencelocs = NULL;
   REAL4Vector *foms = NULL;
   for (ii=2; ii<=rows; ii++) {
      if (ii==2) {
         rowsequencemaxima = XLALCreateREAL4Vector(ii);
         fbinmeanvals = XLALCreateREAL4Vector(ii);
         rowsequencelocs = XLALCreateINT4Vector(ii);
         foms = XLALCreateREAL4Vector(ihsvectorsequence->length-(ii-1));
         if (rowsequencemaxima==NULL) {
            fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ii);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         } else if (fbinmeanvals==NULL) {
            fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ii);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         } else if (rowsequencelocs==NULL) {
            fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", __func__, ii);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         } else if (foms==NULL) {
            fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ihsvectorsequence->length-(ii-1));
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
         
         
         if (params->useSSE) {
            sseSSVectorSequenceSum(tworows, ihsvectorsequence, ihsvectorsequence, 0, 1, 0, (INT4)ihsvectorsequence->length-(ii-1));
            if (xlalErrno!=0) {
               fprintf(stderr, "%s: sseSSVectorSequenceSum() failed.\n", __func__);
               XLAL_ERROR_VOID(XLAL_EFUNC);
            }
         }
         for (jj=0; jj<(INT4)ihsvectorsequence->length-(ii-1); jj++) {
            //Sum IHS values across SFT frequency bins
            if (!params->useSSE) fastSSVectorSequenceSum(tworows, ihsvectorsequence, ihsvectorsequence, jj, jj+1, jj);
            
            //Compute IHS FOM value
            memcpy(rowsequencemaxima->data, &(ihsvalues->data[jj]), sizeof(REAL4)*ii);
            memcpy(fbinmeanvals->data, &(FbinMean->data[jj]), sizeof(REAL4)*ii);
            memcpy(rowsequencelocs->data, &(ihslocations->data[jj]), sizeof(INT4)*ii);
            foms->data[jj] = ihsFOM(rowsequencemaxima, rowsequencelocs, fbinmeanvals, locationnormfactor);
         } /* for jj < ihsvectorsequence->length-(ii-1) */
         
         //sample the IHS values to compute mean, standard deviation, and FAR threshold values
         REAL4Vector *sampledtempihsvals = NULL;
         REAL8 averageval = 0.0, farave = 0.0;
         if ((ihsvectorsequence->length-(ii-1))*ihsvectorsequence->vectorLength>10000) {
            //FILE *tworowvals = fopen("./output/tworowexpectedsample.dat","w");
            sampledtempihsvals = sampleREAL4VectorSequence_nozerosaccepted(tworows, ihsvectorsequence->length-(ii-1), 10000);
            outputfar->ihsdistMean->data[ii-2] = calcMean(sampledtempihsvals);
            for (jj=0; jj<(INT4)sampledtempihsvals->length; jj++) {
               averageval += 1.0;
               if (params->ihsfar != 1.0 && !params->fastchisqinv) farave += gsl_cdf_chisq_Qinv(params->ihsfar, sampledtempihsvals->data[jj]);
               if (params->ihsfar != 1.0 && params->fastchisqinv) farave += cdf_chisq_Qinv(params->ihsfar, sampledtempihsvals->data[jj]);
               if (xlalErrno!=0) {
                  fprintf(stderr, "%s: XXX_chisq_inv() failed.\n", __func__);
                  XLAL_ERROR_VOID(XLAL_EFUNC);
               }
               //fprintf(tworowvals, "%f\n", sampledtempihsvals->data[jj]);
            }
            //fclose(tworowvals);
         } else {
            sampledtempihsvals = XLALCreateREAL4Vector((ihsvectorsequence->length-(ii-1))*ihsvectorsequence->vectorLength);
            memcpy(sampledtempihsvals->data, tworows->data, sizeof(REAL4)*sampledtempihsvals->length);
            outputfar->ihsdistMean->data[ii-2] = calcMean(sampledtempihsvals);
            for (jj=0; jj<(INT4)sampledtempihsvals->length; jj++) {
               if (sampledtempihsvals->data[jj]!=0.0) {
                  averageval += 1.0;
                  if (params->ihsfar != 1.0 && !params->fastchisqinv) farave += gsl_cdf_chisq_Qinv(params->ihsfar, sampledtempihsvals->data[jj]);
                  if (params->ihsfar != 1.0 && params->fastchisqinv) farave += cdf_chisq_Qinv(params->ihsfar, sampledtempihsvals->data[jj]);
                  if (xlalErrno!=0) {
                     fprintf(stderr, "%s: XXX_chisq_inv() failed.\n", __func__);
                     XLAL_ERROR_VOID(XLAL_EFUNC);
                  }
               } /* if sampledtempihsvals->data[jj] != 0.0 */
            } /* for jj < sampledtempihsvals->length */
         }
         outputfar->ihsdistSigma->data[ii-2] = calcStddev(sampledtempihsvals);
         outputfar->ihsfar->data[ii-2] = farave/averageval;
         XLALDestroyREAL4Vector(sampledtempihsvals);
         
         //FOM part
         outputfar->ihsfomdistMean->data[ii-2] = calcMean(foms);
         outputfar->ihsfomdistSigma->data[ii-2] = calcStddev(foms);
         if (params->ihsfomfar!=1.0 && params->ihsfom==0.0) {
            REAL4Vector *smallestfomvals = XLALCreateREAL4Vector((UINT4)roundf((ihsvalues->length-ii+1)*params->ihsfomfar)+1);
            if (smallestfomvals==NULL) {
               fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, (INT4)roundf((ihsvalues->length-ii)*params->ihsfomfar)+1);
               XLAL_ERROR_VOID(XLAL_EFUNC);
            }
            sort_float_smallest(smallestfomvals, foms);
            if (xlalErrno!=0) {
               fprintf(stderr, "%s: sort_float_smallest() failed.\n", __func__);
               XLAL_ERROR_VOID(XLAL_EFUNC);
            }
            outputfar->fomfarthresh->data[ii-2] = smallestfomvals->data[smallestfomvals->length-1];
            XLALDestroyREAL4Vector(smallestfomvals);
         } else if (params->ihsfom!=0.0) {
            outputfar->fomfarthresh->data[ii-2] = params->ihsfom;
         } else {
            outputfar->fomfarthresh->data[ii-2] = -1.0;
         }
         
         XLALDestroyREAL4Vector(rowsequencemaxima);
         rowsequencemaxima = NULL;
         XLALDestroyREAL4Vector(fbinmeanvals);
         fbinmeanvals = NULL;
         XLALDestroyINT4Vector(rowsequencelocs);
         rowsequencelocs = NULL;
         XLALDestroyREAL4Vector(foms);
         foms = NULL;
      } else {
         rowsequencemaxima = XLALCreateREAL4Vector(ii);
         fbinmeanvals = XLALCreateREAL4Vector(ii);
         rowsequencelocs = XLALCreateINT4Vector(ii);
         foms = XLALCreateREAL4Vector(ihsvectorsequence->length-(ii-1));
         if (rowsequencemaxima==NULL) {
            fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ii);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         } else if (fbinmeanvals==NULL) {
            fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ii);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         } else if (rowsequencelocs==NULL) {
            fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", __func__, ii);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         } else if (foms==NULL) {
            fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ihsvectorsequence->length-(ii-1));
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
         
         
         if (params->useSSE) {
            sseSSVectorSequenceSum(tworows, tworows, ihsvectorsequence, 0, ii-1, 0, (INT4)ihsvectorsequence->length-(ii-1));
            if (xlalErrno!=0) {
               fprintf(stderr, "%s: sseSSVectorSequenceSum() failed.\n", __func__);
               XLAL_ERROR_VOID(XLAL_EFUNC);
            }
         }
         for (jj=0; jj<(INT4)ihsvectorsequence->length-(ii-1); jj++) {
            if (!params->useSSE) fastSSVectorSequenceSum(tworows, tworows, ihsvectorsequence, jj, ii-1+jj, jj); //If we didn't use SSE to sum the vector sequence (see lines above)
            
            memcpy(rowsequencemaxima->data, &(ihsvalues->data[jj]), sizeof(REAL4)*ii);
            memcpy(fbinmeanvals->data, &(FbinMean->data[jj]), sizeof(REAL4)*ii);
            memcpy(rowsequencelocs->data, &(ihslocations->data[jj]), sizeof(INT4)*ii);
            foms->data[jj] = ihsFOM(rowsequencemaxima, rowsequencelocs, fbinmeanvals, locationnormfactor);
         } /* for jj< ihsvectorsequence->length - (ii-1) */
         
         
         REAL4Vector *sampledtempihsvals = NULL;
         REAL8 averageval = 0.0, farave = 0.0;
         if ((ihsvectorsequence->length-(ii-1))*ihsvectorsequence->vectorLength>10000) {
            //FILE *row360expect = NULL;
            //if (ii==360) row360expect = fopen("./output/row360expect.dat","w");
            sampledtempihsvals = sampleREAL4VectorSequence_nozerosaccepted(tworows, ihsvectorsequence->length-(ii-1), 10000);
            outputfar->ihsdistMean->data[ii-2] = calcMean(sampledtempihsvals);
            for (jj=0; jj<(INT4)sampledtempihsvals->length; jj++) {
               averageval += 1.0;
               if (params->ihsfar != 1.0 && !params->fastchisqinv) farave += gsl_cdf_chisq_Qinv(params->ihsfar, sampledtempihsvals->data[jj]);
               if (params->ihsfar != 1.0 && params->fastchisqinv) farave += cdf_chisq_Qinv(params->ihsfar, sampledtempihsvals->data[jj]);
               if (xlalErrno!=0) {
                  fprintf(stderr, "%s: XXX_chisq_inv() failed.\n", __func__);
                  XLAL_ERROR_VOID(XLAL_EFUNC);
               }
               //if (ii==360) fprintf(row360expect, "%f\n", sampledtempihsvals->data[jj]);
            }
            //if (ii==360) fclose(row360expect);
         } else {
            sampledtempihsvals = XLALCreateREAL4Vector((ihsvectorsequence->length-(ii-1))*ihsvectorsequence->vectorLength);
            memcpy(sampledtempihsvals->data, tworows->data, sizeof(REAL4)*sampledtempihsvals->length);
            outputfar->ihsdistMean->data[ii-2] = calcMean(sampledtempihsvals);
            for (jj=0; jj<(INT4)sampledtempihsvals->length; jj++) {
               if (sampledtempihsvals->data[jj]!=0.0) {
                  averageval += 1.0;
                  if (params->ihsfar != 1.0 && !params->fastchisqinv) farave += gsl_cdf_chisq_Qinv(params->ihsfar, sampledtempihsvals->data[jj]);
                  if (params->ihsfar != 1.0 && params->fastchisqinv) farave += cdf_chisq_Qinv(params->ihsfar, sampledtempihsvals->data[jj]);
                  if (xlalErrno!=0) {
                     fprintf(stderr, "%s: XXX_chisq_inv() failed.\n", __func__);
                     XLAL_ERROR_VOID(XLAL_EFUNC);
                  }
               }
            }
         }
         outputfar->ihsdistSigma->data[ii-2] = calcStddev(sampledtempihsvals);
         outputfar->ihsfar->data[ii-2] = farave/averageval;
         XLALDestroyREAL4Vector(sampledtempihsvals);
         
         //FOM part
         outputfar->ihsfomdistMean->data[ii-2] = calcMean(foms);
         outputfar->ihsfomdistSigma->data[ii-2] = calcStddev(foms);
         if (params->ihsfomfar!=1.0 && params->ihsfom==0.0) {
            REAL4Vector *smallestfomvals = XLALCreateREAL4Vector((UINT4)roundf((ihsvalues->length-ii+1)*params->ihsfomfar)+1);
            if (smallestfomvals==NULL) {
               fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, (INT4)roundf((ihsvalues->length-ii)*params->ihsfomfar)+1);
               XLAL_ERROR_VOID(XLAL_EFUNC);
            }
            sort_float_smallest(smallestfomvals, foms);
            if (xlalErrno!=0) {
               fprintf(stderr, "%s: sort_float_smallest() failed.\n", __func__);
               XLAL_ERROR_VOID(XLAL_EFUNC);
            }
            outputfar->fomfarthresh->data[ii-2] = smallestfomvals->data[smallestfomvals->length-1];
            XLALDestroyREAL4Vector(smallestfomvals);
         } else if (params->ihsfom!=0.0) {
            outputfar->fomfarthresh->data[ii-2] = params->ihsfom;
         } else {
            outputfar->fomfarthresh->data[ii-2] = -1.0;
         }
         
         XLALDestroyREAL4Vector(rowsequencemaxima);
         rowsequencemaxima = NULL;
         XLALDestroyREAL4Vector(fbinmeanvals);
         fbinmeanvals = NULL;
         XLALDestroyINT4Vector(rowsequencelocs);
         rowsequencelocs = NULL;
         XLALDestroyREAL4Vector(foms);
         foms = NULL;
      }
      
   } /* for ii <= rows */
   
   XLALDestroyREAL4VectorSequence(tworows);
   XLALDestroyREAL4Vector(ihsvalues);
   XLALDestroyINT4Vector(ihslocations);
   
} /*sumIHSSequenceFAR() */
void sumIHSSequence(ihsMaximaStruct *output, ihsfarStruct *inputfar, REAL4VectorSequence *ihsvectorsequence, INT4 rows, REAL4Vector *FbinMean, INT4 locationnormfactor, inputParamsStruct *params)
{
   
   INT4 ii, jj;
   
   REAL4VectorSequence *tworows = XLALCreateREAL4VectorSequence(ihsvectorsequence->length-1, ihsvectorsequence->vectorLength);
   if (tworows==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4VectorSequence(%d,%d) failed.\n", __func__, ihsvectorsequence->length-1, ihsvectorsequence->vectorLength);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   memset(tworows->data, 0, sizeof(REAL4)*tworows->length*tworows->vectorLength);
   
   REAL4Vector *ihsvalues = XLALCreateREAL4Vector(ihsvectorsequence->length);
   INT4Vector *ihslocations = XLALCreateINT4Vector(ihsvalues->length);
   if (ihsvalues==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ihsvectorsequence->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (ihslocations==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ihsvalues->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   for (ii=0; ii<(INT4)ihsvalues->length; ii++) {
      ihslocations->data[ii] = max_index_from_vector_in_REAL4VectorSequence(ihsvectorsequence, ii) + 5;
      ihsvalues->data[ii] = ihsvectorsequence->data[ii*ihsvectorsequence->vectorLength + ihslocations->data[ii]-5];
   }
   
   REAL4Vector *excessabovenoise = XLALCreateREAL4Vector(ihsvectorsequence->vectorLength);
   REAL4Vector *scaledExpectedIHSVectorValues = XLALCreateREAL4Vector(ihsvectorsequence->vectorLength);
   if (excessabovenoise==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ihsvectorsequence->vectorLength);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (scaledExpectedIHSVectorValues==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ihsvectorsequence->vectorLength);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   
   REAL4Vector *rowsequencemaxima = NULL;
   REAL4Vector *fbinmeanvals = NULL;
   INT4Vector *rowsequencelocs = NULL;
   for (ii=1; ii<=rows; ii++) {
      if (ii==1) {
         memcpy(output->maximaForEachFbin->data, ihsvalues->data, sizeof(REAL4)*ihsvalues->length);
         memcpy(output->locationsForEachFbin->data, ihslocations->data, sizeof(INT4)*ihslocations->length);
      } else if (ii==2) {
         rowsequencemaxima = XLALCreateREAL4Vector(ii);
         fbinmeanvals = XLALCreateREAL4Vector(ii);
         rowsequencelocs = XLALCreateINT4Vector(ii);
         if (rowsequencemaxima==NULL) {
            fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ii);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         } else if (fbinmeanvals==NULL) {
            fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ii);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         } else if (rowsequencelocs==NULL) {
            fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", __func__, ii);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
         
         REAL4 sumofnoise = 0.0;    //To scale the expected IHS background
         
         
         if (params->useSSE || params->validateSSE) {
            //Sum oup the IHS vectors using SSE function
            sseSSVectorSequenceSum(tworows, ihsvectorsequence, ihsvectorsequence, 0, 1, 0, (INT4)ihsvectorsequence->length-(ii-1));
            if (xlalErrno!=0) {
               fprintf(stderr, "%s: sseSSVectorSequenceSum() failed.\n", __func__);
               XLAL_ERROR_VOID(XLAL_EFUNC);
            }
            /* FILE *tworowreal = fopen("./output/tworowsumreal.dat","w");
            for (jj=0; jj<(INT4)(tworows->length*tworows->vectorLength); jj++) fprintf(tworowreal, "%f\n", tworows->data[jj]);
            fclose(tworowreal); */
            
            /* validate SSE code */
            if (params->validateSSE) {
               REAL4VectorSequence *tworows_valid = XLALCreateREAL4VectorSequence(tworows->length, tworows->vectorLength);
               if (tworows_valid==NULL) {
                  fprintf(stderr,"%s: XLALCreateREAL4VectorSequence(%d,%d) failed.\n", __func__, tworows->length, tworows->vectorLength);
                  XLAL_ERROR_VOID(XLAL_EFUNC);
               }
               for (jj=0; jj<(INT4)ihsvectorsequence->length-(ii-1); jj++) fastSSVectorSequenceSum(tworows_valid, ihsvectorsequence, ihsvectorsequence, jj, jj+1, jj);
               for (jj=0; jj<(INT4)(tworows->length*tworows->vectorLength); jj++) {
                  if (tworows->data[jj] != tworows_valid->data[jj]) {
                     fprintf(stderr,"%s: sseSSVectorSequenceSum() failed to produce valid results.\n", __func__);
                     XLAL_ERROR_VOID(XLAL_EFUNC);
                  }
               }
               XLALDestroyREAL4VectorSequence(tworows_valid);
            } /* validate SSE code */
         }
         for (jj=0; jj<(INT4)ihsvectorsequence->length-(ii-1); jj++) {
            //Sum IHS values across SFT frequency bins if the SSE function wasn't used
            if (!params->useSSE) fastSSVectorSequenceSum(tworows, ihsvectorsequence, ihsvectorsequence, jj, jj+1, jj);
            
            memcpy(fbinmeanvals->data, &(FbinMean->data[jj]), sizeof(REAL4)*ii);
            
            //To scale the background efficiently
            if (jj==0) {
               INT4 kk;
               for (kk=0; kk<ii; kk++) sumofnoise += fbinmeanvals->data[kk];
            } else {
               sumofnoise -= FbinMean->data[jj-1];
               sumofnoise += FbinMean->data[jj+(ii-1)];
            }
            
            if (params->useSSE || (params->validateSSE && jj==0)) {
               //Scale the expected IHS vector
               sseScaleREAL4Vector(scaledExpectedIHSVectorValues, inputfar->expectedIHSVector, sumofnoise);
               if (xlalErrno!=0) {
                  fprintf(stderr, "%s: sseScaleREAL4Vector() failed.\n", __func__);
                  XLAL_ERROR_VOID(XLAL_EFUNC);
               }
               
               /* validate SSE code */
               if (params->validateSSE && jj==0) {
                  REAL4Vector *scaledExpectedIHSVectorValues_valid = XLALCreateREAL4Vector(scaledExpectedIHSVectorValues->length);
                  if (scaledExpectedIHSVectorValues_valid==NULL) {
                     fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, scaledExpectedIHSVectorValues->length);
                     XLAL_ERROR_VOID(XLAL_EFUNC);
                  }
                  INT4 kk;
                  REAL4 scaleval = sumofnoise;
                  for (kk=0; kk<(INT4)inputfar->expectedIHSVector->length; kk++) {
                     scaledExpectedIHSVectorValues_valid->data[kk] = scaleval*inputfar->expectedIHSVector->data[kk];
                     if (scaledExpectedIHSVectorValues_valid->data[kk] != scaledExpectedIHSVectorValues->data[kk]) {
                        fprintf(stderr,"%s: sseScaleREAL4Vector() failed to produce valid results.\n", __func__);
                        XLAL_ERROR_VOID(XLAL_EFUNC);
                     }
                  }
                  XLALDestroyREAL4Vector(scaledExpectedIHSVectorValues_valid);
               } /* validate SSE code */
               
               //subtract the noise from the data
               sseSSVectorSequenceSubtract(excessabovenoise, tworows, scaledExpectedIHSVectorValues, jj);
               if (xlalErrno!=0) {
                  fprintf(stderr, "%s: sseSSVectorSequenceSubtract() failed.\n", __func__);
                  XLAL_ERROR_VOID(XLAL_EFUNC);
               }
               
               /* validate SSE code */
               if (params->validateSSE && jj==0) {
                  REAL4Vector *excessabovenoise_valid = XLALCreateREAL4Vector(excessabovenoise->length);
                  if (excessabovenoise_valid==NULL) {
                     fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, excessabovenoise->length);
                     XLAL_ERROR_VOID(XLAL_EFUNC);
                  }
                  fastSSVectorSequenceSubtract(excessabovenoise, tworows, scaledExpectedIHSVectorValues, jj);
                  INT4 kk;
                  for (kk=0; kk<(INT4)excessabovenoise->length; kk++) {
                     if (excessabovenoise_valid->data[kk] != excessabovenoise->data[kk]) {
                        fprintf(stderr,"%s: sseSSVectorSequenceSubtract() failed to produce valid results.\n", __func__);
                        XLAL_ERROR_VOID(XLAL_EFUNC);
                     }
                  }
                  XLALDestroyREAL4Vector(excessabovenoise_valid);
               } /* validate SSE code */
               
            } else {
               INT4 kk;
               REAL4 scaleval = sumofnoise;
               for (kk=0; kk<(INT4)inputfar->expectedIHSVector->length; kk++) scaledExpectedIHSVectorValues->data[kk] = scaleval*inputfar->expectedIHSVector->data[kk];
               fastSSVectorSequenceSubtract(excessabovenoise, tworows, scaledExpectedIHSVectorValues, jj);
            }
            
            //Compute the maximum IHS value in the second FFT frequency direction
            output->locations->data[jj] = max_index(excessabovenoise) + 5;
            output->maxima->data[jj] = tworows->data[jj*tworows->vectorLength + (output->locations->data[jj]-5)];
            
            //Compute IHS FOM value
            memcpy(rowsequencemaxima->data, &(ihsvalues->data[jj]), sizeof(REAL4)*ii);
            memcpy(rowsequencelocs->data, &(ihslocations->data[jj]), sizeof(INT4)*ii);
            output->foms->data[jj] = ihsFOM(rowsequencemaxima, rowsequencelocs, fbinmeanvals, locationnormfactor);
         } /* for jj < ihsvectorsequence->length-(ii-1) */
         
         XLALDestroyREAL4Vector(rowsequencemaxima);
         rowsequencemaxima = NULL;
         XLALDestroyREAL4Vector(fbinmeanvals);
         fbinmeanvals = NULL;
         XLALDestroyINT4Vector(rowsequencelocs);
         rowsequencelocs = NULL;
      } else {
         rowsequencemaxima = XLALCreateREAL4Vector(ii);
         fbinmeanvals = XLALCreateREAL4Vector(ii);
         rowsequencelocs = XLALCreateINT4Vector(ii);
         if (rowsequencemaxima==NULL) {
            fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ii);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         } else if (fbinmeanvals==NULL) {
            fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ii);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         } else if (rowsequencelocs==NULL) {
            fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", __func__, ii);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
         
         REAL4 sumofnoise = 0.0;    //To scale the expected IHS background
         INT4 endloc = ((ii-1)*(ii-1)-(ii-1))/2;
         
         if (params->useSSE) {
            sseSSVectorSequenceSum(tworows, tworows, ihsvectorsequence, 0, ii-1, 0, (INT4)ihsvectorsequence->length-(ii-1));
            if (xlalErrno!=0) {
               fprintf(stderr, "%s: sseSSVectorSequenceSum() failed.\n", __func__);
               XLAL_ERROR_VOID(XLAL_EFUNC);
            }
            /* if (ii==360) {
               FILE *row360real = fopen("./output/row360sumreal.dat","w");
               for (jj=0; jj<(INT4)((tworows->length-(ii-2))*tworows->vectorLength); jj++) fprintf(row360real, "%f\n", tworows->data[jj]);
               fclose(row360real);
            } */
         }
         for (jj=0; jj<(INT4)ihsvectorsequence->length-(ii-1); jj++) {
            if (!params->useSSE) fastSSVectorSequenceSum(tworows, tworows, ihsvectorsequence, jj, ii-1+jj, jj); //If we didn't use SSE to sum the vector sequence (see lines above)
            
            memcpy(fbinmeanvals->data, &(FbinMean->data[jj]), sizeof(REAL4)*ii);
            
            //To scale the background efficiently
            if (jj==0) {
               INT4 kk;
               for (kk=0; kk<ii; kk++) sumofnoise += fbinmeanvals->data[kk];
            } else {
               sumofnoise -= FbinMean->data[jj-1];
               sumofnoise += FbinMean->data[jj+(ii-1)];
            }
            
            if (params->useSSE) {
               sseScaleREAL4Vector(scaledExpectedIHSVectorValues, inputfar->expectedIHSVector, sumofnoise);
               if (xlalErrno!=0) {
                  fprintf(stderr, "%s: sseScaleREAL4Vector() failed.\n", __func__);
                  XLAL_ERROR_VOID(XLAL_EFUNC);
               }
               sseSSVectorSequenceSubtract(excessabovenoise, tworows, scaledExpectedIHSVectorValues, jj);
               if (xlalErrno!=0) {
                  fprintf(stderr, "%s: sseSSVectorSequenceSubtract() failed.\n", __func__);
                  XLAL_ERROR_VOID(XLAL_EFUNC);
               }
            } else {
               INT4 kk;
               REAL4 scaleval = sumofnoise;
               for (kk=0; kk<(INT4)inputfar->expectedIHSVector->length; kk++) scaledExpectedIHSVectorValues->data[kk] = scaleval*inputfar->expectedIHSVector->data[kk];
               fastSSVectorSequenceSubtract(excessabovenoise, tworows, scaledExpectedIHSVectorValues, jj);
            }
            
            output->locations->data[(ii-2)*ihsvalues->length-endloc+jj] = max_index(excessabovenoise) + 5;
            output->maxima->data[(ii-2)*ihsvalues->length-endloc+jj] = tworows->data[jj*tworows->vectorLength + (output->locations->data[(ii-2)*ihsvalues->length-endloc+jj]-5)];
            
            memcpy(rowsequencemaxima->data, &(ihsvalues->data[jj]), sizeof(REAL4)*ii);
            memcpy(rowsequencelocs->data, &(ihslocations->data[jj]), sizeof(INT4)*ii);
            output->foms->data[(ii-2)*ihsvalues->length-endloc+jj] = ihsFOM(rowsequencemaxima, rowsequencelocs, fbinmeanvals, locationnormfactor);
         } /* for jj< ihsvectorsequence->length - (ii-1) */
         
         
         XLALDestroyREAL4Vector(rowsequencemaxima);
         rowsequencemaxima = NULL;
         XLALDestroyREAL4Vector(fbinmeanvals);
         fbinmeanvals = NULL;
         XLALDestroyINT4Vector(rowsequencelocs);
         rowsequencelocs = NULL;
      }
      
   } /* for ii <= rows */
   
   XLALDestroyREAL4VectorSequence(tworows);
   XLALDestroyREAL4Vector(scaledExpectedIHSVectorValues);
   XLALDestroyREAL4Vector(excessabovenoise);
   XLALDestroyREAL4Vector(ihsvalues);
   XLALDestroyINT4Vector(ihslocations);
   
} /*sumIHSSequence() */

void SSVectorSequenceSum(REAL4VectorSequence *output, REAL4VectorSequence *input1, REAL4VectorSequence *input2, INT4 vectorpos1, INT4 vectorpos2, INT4 outputvectorpos)
{
   
   INT4 ii, vec1 = vectorpos1*input1->vectorLength, vec2 = vectorpos2*input2->vectorLength, outvec = outputvectorpos*output->vectorLength;
   for (ii=0; ii<(INT4)input1->vectorLength; ii++) output->data[outvec + ii] = input1->data[vec1 + ii] + input2->data[vec2 + ii];
   
}
void fastSSVectorSequenceSum(REAL4VectorSequence *output, REAL4VectorSequence *input1, REAL4VectorSequence *input2, INT4 vectorpos1, INT4 vectorpos2, INT4 outputvectorpos)
{
   
   REAL4 *a, *b, *c;
   INT4 n;
   
   a = &(input1->data[vectorpos1*input1->vectorLength]);
   b = &(input2->data[vectorpos2*input2->vectorLength]);
   c = &(output->data[outputvectorpos*output->vectorLength]);
   n = output->vectorLength;
   
   while (n-- > 0) {
      *c = (*a)+(*b);
      a++;
      b++;
      c++;
   }
   
}
void sseSSVectorSequenceSum(REAL4VectorSequence *output, REAL4VectorSequence *input1, REAL4VectorSequence *input2, INT4 vectorpos1, INT4 vectorpos2, INT4 outputvectorpos, INT4 numvectors)
{
   
#ifdef __SSE__
   INT4 roundedvectorlength = (INT4)input1->vectorLength / 4;
   
   REAL4* allocinput1 = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15);
   REAL4* allocinput2 = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15);
   REAL4* allocoutput = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15);
   if (allocinput1==NULL) {
      fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 4*roundedvectorlength*sizeof(REAL4) + 15);
      XLAL_ERROR_VOID(XLAL_ENOMEM);
   } else if (allocinput2==NULL) {
      fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 4*roundedvectorlength*sizeof(REAL4) + 15);
      XLAL_ERROR_VOID(XLAL_ENOMEM);
   } else if (allocoutput==NULL) {
      fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 4*roundedvectorlength*sizeof(REAL4) + 15);
      XLAL_ERROR_VOID(XLAL_ENOMEM);
   }
   REAL4* alignedinput1 = (void*)(((UINT8)allocinput1+15) & ~15);
   REAL4* alignedinput2 = (void*)(((UINT8)allocinput2+15) & ~15);
   REAL4* alignedoutput = (void*)(((UINT8)allocoutput+15) & ~15);
   
   INT4 ii, jj;
   for (ii=0; ii<numvectors; ii++) {
      INT4 vec1 = (vectorpos1+ii)*input1->vectorLength, vec2 = (vectorpos2+ii)*input2->vectorLength, outvec = (outputvectorpos+ii)*output->vectorLength;
      
      INT4 input1vecaligned = 0, input2vecaligned = 0, outputvecaligned = 0;
      __m128 *arr1, *arr2, *result;
      if ( &(input1->data[vec1])==(void*)(((UINT8)&(input1->data[vec1])+15) & ~15) ) {
         input1vecaligned = 1;
         arr1 = (__m128*)&(input1->data[vec1]);
      } else {
         memcpy(alignedinput1, &(input1->data[vec1]), sizeof(REAL4)*4*roundedvectorlength);
         arr1 = (__m128*)alignedinput1;
      }
      if ( &(input2->data[vec2])==(void*)(((UINT8)&(input2->data[vec2])+15) & ~15) ) {
         input2vecaligned = 1;
         arr2 = (__m128*)&(input2->data[vec2]);
      } else {
         memcpy(alignedinput2, &(input2->data[vec2]), sizeof(REAL4)*4*roundedvectorlength);
         arr2 = (__m128*)alignedinput2;
      }
      if ( &(output->data[outvec])==(void*)(((UINT8)&(output->data[outvec])+15) & ~15) ) {
         outputvecaligned = 1;
         result = (__m128*)&(output->data[outvec]);
      } else result = (__m128*)alignedoutput;
      
      for (jj=0; jj<roundedvectorlength; jj++) {
         *result = _mm_add_ps(*arr1, *arr2);
         arr1++;
         arr2++;
         result++;
      }
      
      if (!outputvecaligned) memcpy(&(output->data[outvec]), alignedoutput, sizeof(REAL4)*4*roundedvectorlength);
      
      REAL4 *a = &(input1->data[vec1+4*roundedvectorlength]);
      REAL4 *b = &(input2->data[vec2+4*roundedvectorlength]);
      REAL4 *c = &(output->data[outvec+4*roundedvectorlength]);
      INT4 n = output->vectorLength-4*roundedvectorlength;
      while (n-- > 0) {
         *c = (*a)+(*b);
         a++;
         b++;
         c++;
      }
   }
   
   XLALFree(allocinput1);
   XLALFree(allocinput2);
   XLALFree(allocoutput);
#else
   fprintf(stderr, "%s: Failed because SSE is not supported, possibly because -msse flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR_VOID(XLAL_EFAILED);
#endif
   
}
void fastSSVectorSequenceSubtract(REAL4Vector *output, REAL4VectorSequence *input1, REAL4Vector *input2, INT4 vectorpos1)
{
   
   REAL4 *a, *b, *c;
   INT4 n;
   
   a = &(input1->data[vectorpos1*input1->vectorLength]);
   b = input2->data;
   c = output->data;
   n = output->length;
   
   while (n-- > 0) {
      *c = (*a)-(*b);
      a++;
      b++;
      c++;
   }
   
}
void sseSSVectorSequenceSubtract(REAL4Vector *output, REAL4VectorSequence *input1, REAL4Vector *input2, INT4 vectorpos1)
{
   
#ifdef __SSE__
   INT4 roundedvectorlength = (INT4)input1->vectorLength / 4;
   INT4 vec1 = vectorpos1*input1->vectorLength, ii = 0;
   INT4 vec1aligned = 0, vec2aligned = 0, outputaligned = 0;
   REAL4 *allocinput1 = NULL, *allocinput2 = NULL, *allocoutput = NULL, *alignedinput1 = NULL, *alignedinput2 = NULL, *alignedoutput = NULL;
   __m128 *arr1, *arr2, *result;
   
   //Allocate memory for aligning input vector 1 if necessary
   if ( &(input1->data[vec1])==(void*)(((UINT8)&(input1->data[vec1])+15) & ~15) ) {
      vec1aligned = 1;
      arr1 = (__m128*)&(input1->data[vec1]);
   } else {
      allocinput1 = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15);
      if (allocinput1==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 4*roundedvectorlength*sizeof(REAL4) + 15);
         XLAL_ERROR_VOID(XLAL_ENOMEM);
      }
      alignedinput1 = (void*)(((UINT8)allocinput1+15) & ~15);
      memcpy(alignedinput1, &(input1->data[vec1]), sizeof(REAL4)*4*roundedvectorlength);
      arr1 = (__m128*)alignedinput1;
   }
   
   //Allocate memory for aligning input vector 2 if necessary
   if ( input2->data==(void*)(((UINT8)input2->data+15) & ~15) ) {
      vec2aligned = 1;
      arr2 = (__m128*)input2->data;
   } else {
      allocinput2 = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15);
      if (allocinput2==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 4*roundedvectorlength*sizeof(REAL4) + 15);
         XLAL_ERROR_VOID(XLAL_ENOMEM);
      }
      alignedinput2 = (void*)(((UINT8)allocinput2+15) & ~15);
      memcpy(alignedinput2, input2->data, sizeof(REAL4)*4*roundedvectorlength);
      arr2 = (__m128*)alignedinput2;
   }
   
   //Allocate memory for aligning output vector if necessary
   if ( output->data==(void*)(((UINT8)output->data+15) & ~15) ) {
      outputaligned = 1;
      result = (__m128*)output->data;
   } else {
      allocoutput = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15);
      if (allocoutput==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 4*roundedvectorlength*sizeof(REAL4) + 15);
         XLAL_ERROR_VOID(XLAL_ENOMEM);
      }
      alignedoutput = (void*)(((UINT8)allocoutput+15) & ~15);
      result = (__m128*)alignedoutput;
   }
   
   for (ii=0; ii<roundedvectorlength; ii++) {
      *result = _mm_sub_ps(*arr1, *arr2);
      arr1++;
      arr2++;
      result++;
   }
   
   if (!outputaligned) memcpy(output->data, alignedoutput, sizeof(REAL4)*4*roundedvectorlength);
   
   REAL4 *a = &(input1->data[vec1+4*roundedvectorlength]);
   REAL4 *b = &(input2->data[4*roundedvectorlength]);
   REAL4 *c = &(output->data[4*roundedvectorlength]);
   INT4 n = output->length-4*roundedvectorlength;
   while (n-- > 0) {
      *c = (*a)-(*b);
      a++;
      b++;
      c++;
   }
   
   if (!vec1aligned) XLALFree(allocinput1);
   if (!vec2aligned) XLALFree(allocinput2);
   if (!outputaligned) XLALFree(allocoutput);
#else
   fprintf(stderr, "%s: Failed because SSE is not supported, possibly because -msse flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR_VOID(XLAL_EFAILED);
#endif
   
}
/* void addScaledVectorIntoVectorSequence(REAL4VectorSequence *output, REAL4Vector *input, REAL4 scale, INT4 outputvectorpos)
{
   
   INT4 ii;
   for (ii=0; ii<(INT4)input->length; ii++) output->data[outputvectorpos*output->vectorLength + ii] += (REAL4)(scale*input->data[ii]);
   
}
void sseAddScaledVectorIntoVectorSequence(REAL4VectorSequence *output, REAL4Vector *input, REAL4 scale, INT4 outputvectorpos)
{
   
#ifdef __SSE__
   INT4 roundedvectorlength = (INT4)input->length / 4;
   INT4 vectoraligned = 0, outputaligned = 0, ii = 0, outvec = outputvectorpos*output->vectorLength;
   
   REAL4 *allocinput = NULL, *allocoutput = NULL, *alignedinput = NULL, *alignedoutput = NULL;
   __m128 *arr1, scaledarr1, *result;
   
   __m128 scalefactor = _mm_set1_ps(scale);
   
   //Allocate memory for aligning input vector if necessary
   if ( input->data==(void*)(((UINT8)input->data+15) & ~15) ) {
      vectoraligned = 1;
      arr1 = (__m128*)input->data;
   } else {
      allocinput = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15);
      if (allocinput==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 4*roundedvectorlength*sizeof(REAL4) + 15);
         XLAL_ERROR_VOID(XLAL_ENOMEM);
      }
      alignedinput = (void*)(((UINT8)allocinput+15) & ~15);
      memcpy(alignedinput, input->data, sizeof(REAL4)*4*roundedvectorlength);
      arr1 = (__m128*)alignedinput;
   }
   
   //Allocate memory for aligning output vector if necessary
   if ( &(output->data[outvec])==(void*)(((UINT8)&(output->data[outvec])+15) & ~15) ) {
      outputaligned = 1;
      result = (__m128*)&(output->data[outvec]);
   } else {
      allocoutput = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15);
      if (allocoutput==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 4*roundedvectorlength*sizeof(REAL4) + 15);
         XLAL_ERROR_VOID(XLAL_ENOMEM);
      }
      alignedoutput = (void*)(((UINT8)allocoutput+15) & ~15);
      memcpy(alignedoutput, output->data, sizeof(REAL4)*4*roundedvectorlength);
      result = (__m128*)alignedoutput;
   }
   
   //Scale the input array, and add the scaled array into the output
   for (ii=0; ii<roundedvectorlength; ii++) {
      scaledarr1 = _mm_mul_ps(*arr1, scalefactor);
      *result = _mm_add_ps(*result, scaledarr1);
      arr1++;
      result++;
   }
   
   //Copy output aligned memory to non-aligned memory if necessary
   if (!outputaligned) memcpy(&(output->data[outvec]), alignedoutput, 4*roundedvectorlength*sizeof(REAL4));
   
   //Finish up the remaining part
   for (ii=4*roundedvectorlength; ii<(INT4)input->length; ii++) output->data[outvec + ii] += (REAL4)(scale*input->data[ii]);
   
   //Free memory if necessary
   if (!vectoraligned) XLALFree(allocinput);
   if (!outputaligned) XLALFree(allocoutput);
   
#else
   fprintf(stderr, "%s: Failed because SSE is not supported, possibly because -msse flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR_VOID(XLAL_EFAILED);
#endif
   
} */


REAL4VectorSequence * ihsVectorSums(REAL4VectorSequence *input, INT4 rows)
{
   
   INT4 ii, jj, kk;
   
   REAL4VectorSequence *output = XLALCreateREAL4VectorSequence(input->length-(rows-1), input->vectorLength);
   REAL8VectorSequence *tworows = XLALCreateREAL8VectorSequence(input->length-1, input->vectorLength);
   
   //Start with 2 rows
   for (jj=0; jj<(INT4)input->length-(2-1); jj++) {
      for (kk=0; kk<(INT4)input->vectorLength; kk++) tworows->data[jj*input->vectorLength + kk] = input->data[jj*input->vectorLength + kk] + input->data[(jj+1)*input->vectorLength + kk];
   }
   
   //contintue with more rows
   for (ii=3; ii<=rows; ii++) {
      for (jj=0; jj<(INT4)input->length-(ii-1); jj++) {
         for (kk=0; kk<(INT4)input->vectorLength; kk++) tworows->data[jj*input->vectorLength + kk] = tworows->data[jj*input->vectorLength + kk] + input->data[(ii-1+jj)*input->vectorLength + kk];
      }
   }
   
   //copy data summed into the output vector
   //memcpy(output->data, tworows->data, sizeof(REAL8)*input->vectorLength*(input->length-(rows-1)));
   for (ii=0; ii<(INT4)(input->vectorLength*(input->length-(rows-1))); ii++) output->data[ii] = (REAL4)tworows->data[ii];
   XLALDestroyREAL8VectorSequence(tworows);
   
   return output;
   
}


//////////////////////////////////////////////////////////////
// Calculate the IHS FOM for a number of rows  -- 
REAL4 ihsFOM(REAL4Vector *ihss, INT4Vector *locs, REAL4Vector *sigma, INT4 locationnormfactor)
//REAL4 ihsFOM(REAL4Vector *ihss, INT4Vector *locs, REAL4 sigma)
{
   
   
   //INT4 ii, maxsnrloc;
   //REAL4 maxsnr;
   INT4 ii;
   
   /* //Find which pair has the best combined SNR (RMS) and the location
   maxsnr = sqrtf(ihss->data[0]*ihss->data[0]/(sigma->data[0]*sigma->data[0]) + ihss->data[ihss->length-1]*ihss->data[ihss->length-1]/(sigma->data[sigma->length-1]*sigma->data[sigma->length-1]));
   maxsnrloc = 0;
   for (ii=1; ii<(INT4)(ihss->length*0.5); ii++) {
      REAL4 snr = sqrtf(ihss->data[ii]*ihss->data[ii]/(sigma->data[ii]*sigma->data[ii]) + ihss->data[ihss->length-ii-1]*ihss->data[ihss->length-ii-1]/(sigma->data[sigma->length-ii-1]*sigma->data[sigma->length-ii-1]));
      if (snr>maxsnr) {
         maxsnr = snr;
         maxsnrloc = ii;
      }
   }
   //For the highest SNR pair, compute the FOM
   REAL4 fom = 12.0*(ihss->data[maxsnrloc]*ihss->data[ihss->length-maxsnrloc-1]/(ihss->data[maxsnrloc]*sigma->data[sigma->length-maxsnrloc-1] + ihss->data[ihss->length-maxsnrloc-1]*sigma->data[maxsnrloc]))*(REAL4)((locs->data[maxsnrloc] - locs->data[locs->length-maxsnrloc-1]) * (locs->data[maxsnrloc] - locs->data[locs->length-maxsnrloc-1])); */
   
   REAL4 fom = 0.0;
   for (ii=0; ii<(INT4)(ihss->length*0.5); ii++) {
      if (locs->data[ii] == locs->data[locs->length-ii-1]) fom += ihss->data[ii]*ihss->data[ihss->length-ii-1]/(ihss->data[ii]*sigma->data[ihss->length-ii-1]+ihss->data[ihss->length-ii-1]*sigma->data[ii])*((0.5/locationnormfactor) * (0.5/locationnormfactor));
      else fom += ihss->data[ii]*ihss->data[ihss->length-ii-1]/(ihss->data[ii]*sigma->data[ihss->length-ii-1]+ihss->data[ihss->length-ii-1]*sigma->data[ii])*(((REAL4)locs->data[ii]/locationnormfactor - (REAL4)locs->data[locs->length-ii-1]/locationnormfactor) * ((REAL4)locs->data[ii]/locationnormfactor - (REAL4)locs->data[locs->length-ii-1]/locationnormfactor));
   }
   fom *= 12.0;
   
   return fom;

} /* ihsFOM() */


//////////////////////////////////////////////////////////////
// Calculate a guess for the location of the brightest pixels
REAL4 ihsLoc(REAL4Vector *ihss, INT4Vector *locs, REAL4Vector *sigma)
//REAL4 ihsLoc(REAL4Vector *ihss, INT4Vector *locs, REAL4 sigma)
{
   
   INT4 ii, maxsnrloc;
   REAL4 maxsnr;
   
   //Find which pair has the best combined SNR (RMS) and the location
   maxsnr = sqrtf(ihss->data[0]*ihss->data[0]/(sigma->data[0]*sigma->data[0]) + ihss->data[ihss->length-1]*ihss->data[ihss->length-1]/(sigma->data[sigma->length-1]*sigma->data[sigma->length-1]));
   maxsnrloc = 0;
   for (ii=1; ii<(INT4)(ihss->length*0.5); ii++) {
      REAL4 snr = sqrtf(ihss->data[ii]*ihss->data[ii]/(sigma->data[ii]*sigma->data[ii]) + ihss->data[ihss->length-ii-1]*ihss->data[ihss->length-ii-1]/(sigma->data[sigma->length-ii-1]*sigma->data[sigma->length-ii-1]));
      if (snr>maxsnr) {
         maxsnr = snr;
         maxsnrloc = ii;
      }
   }
   //For the highest SNR pair, compute the location
   REAL4 location = 0.5*(locs->data[maxsnrloc] + locs->data[locs->length-maxsnrloc-1]);
   
   /* INT4 ii, mindiff = -1, minloc = 0;
   for (ii=0; ii<(INT4)(ihss->length*0.5); ii++) {
      INT4 diff = abs(locs->data[ii]-locs->data[locs->length-1-ii]);
      if (mindiff==-1 || (mindiff>-1 && diff<mindiff)) {
         mindiff = diff;
         minloc = ii;
      }
   }
   REAL4 location = 0.5*(locs->data[minloc] + locs->data[locs->length-minloc-1]); */
   
   return location;

} /* ihsLoc() */



void findIHScandidates(candidateVector *candlist, ihsfarStruct *ihsfarstruct, inputParamsStruct *params, ffdataStruct *ffdata, ihsMaximaStruct *ihsmaxima, REAL4Vector *fbinavgs, REAL4VectorSequence *trackedlines)
{
   
   INT4 ii, jj, kk;
   REAL8 fsig, per0, B;
   
   INT4 numberofIHSvalsChecked = 0, numberofIHSvalsExceededThresh = 0, numberPassingBoth = 0, linesinterferewithnum = 0;
   
   INT4 numfbins = ffdata->numfbins;
   
   INT4 minrows = (INT4)round(2.0*params->dfmin*params->Tcoh)+1;
   
   REAL4Vector *ihss, *avgsinrange;
   INT4Vector *locs;
   
   //Check the IHS values against the FAR, checking between IHS width values
   //FILE *IHSVALSOUTPUT = fopen("./output/allihsvalspassthresh.dat","w");
   for (ii=minrows; ii<=(INT4)ihsfarstruct->ihsfar->length+1; ii++) {
      ihss = XLALCreateREAL4Vector(ii);
      locs = XLALCreateINT4Vector(ii);
      avgsinrange = XLALCreateREAL4Vector(ii);
      if (ihss==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ii);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      } else if (locs==NULL) {
         fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", __func__, ii);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      } else if (avgsinrange==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ii);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      }
      
      //REAL8 highestval = 0.0;
      //INT4 highestvalloc = -1, jjloc = 0;
      for (jj=0; jj<(INT4)numfbins-(ii-1); jj++) {
      
         //Noise in the range of the rows, mean for IHS
         memcpy(avgsinrange->data, &(fbinavgs->data[jj]), sizeof(REAL4)*ii);
         REAL4 meanNoise = calcMean(avgsinrange);
         if (XLAL_IS_REAL4_FAIL_NAN(meanNoise)) {
            fprintf(stderr,"%s: calcMean() failed.\n", __func__);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
         
         numberofIHSvalsChecked++;
         
         INT4 locationinmaximastruct = (ii-2)*numfbins-((ii-1)*(ii-1)-(ii-1))/2+jj;
         
         //Check the IHS sum against the FAR (scaling FAR with mean of the noise in the range of rows)
         if (ihsmaxima->maxima->data[locationinmaximastruct] > ihsfarstruct->ihsfar->data[ii-2]*meanNoise) {
            
            numberofIHSvalsExceededThresh++;
            
            if (ihsfarstruct->fomfarthresh->data[ii-2]==-1.0 || ihsmaxima->foms->data[locationinmaximastruct]<=ihsfarstruct->fomfarthresh->data[ii-2]) {
               
               numberPassingBoth++;
               //fprintf(IHSVALSOUTPUT, "%.6f %.6f\n", 0.5*(ii-1)/params->Tcoh, ihsmaxima->maxima->data[locationinmaximastruct]);
               
               INT4 loc = ihsmaxima->locations->data[locationinmaximastruct];
               per0 = params->Tobs/loc;                                          //Candidate period
               fsig = params->fmin + (0.5*(ii-1) + jj)/params->Tcoh;             //Candidate frequency
               B = 0.5*(ii-1)/params->Tcoh;                                      //Candidate modulation depth
               
               //Test to see if any tracked lines are overlapping the candidate signal
               INT4 nolinesinterfering = 1;
               if (trackedlines!=NULL) {
                  kk = 0;
                  while (kk<(INT4)trackedlines->length && nolinesinterfering==1) {
                     INT4 ll = 0;
                     while (ll<(INT4)trackedlines->vectorLength && nolinesinterfering==1) {
                        if (trackedlines->data[kk*trackedlines->vectorLength + ll]<=fsig+B && trackedlines->data[kk*trackedlines->vectorLength + ll]>=fsig-B) {
                           nolinesinterfering = 0;
                        }
                        ll++;
                     } /* while trackedlines->vectorLength && nolinesinterfering==1 */
                     kk++;
                  } /* while kk < trackedlines->length && nolinesinterfering==1 */
               } /* if trackedlines != NULL */
               
               if (!nolinesinterfering) {
                  linesinterferewithnum++;
               } else {
                  //Test this here
                  REAL8 noise = ihsfarstruct->expectedIHSVector->data[loc-5];
                  REAL8 totalnoise = 0.0;
                  for (kk=0; kk<ii; kk++) totalnoise += noise*fbinavgs->data[jj+kk];
                  
                  /* if (ihsmaxima->maxima->data[locationinmaximastruct]-totalnoise > highestval) {
                     //highestval = ihsmaxima->maxima->data[locationinmaximastruct];
                     highestval = ihsmaxima->maxima->data[locationinmaximastruct]-totalnoise;
                     highestvalloc = locationinmaximastruct;
                     jjloc = jj;
                  } */

                  //Candidate h0
                  //REAL8 h0 = ihs2h0_withNoiseSubtraction(ihsmaxima->maxima->data[locationinmaximastruct], loc, jj, ii, params, aveNoise, fbinavgs);
                  REAL8 h0 = ihs2h0(2.0*(ihsmaxima->maxima->data[locationinmaximastruct]-totalnoise), params);
                  if (candlist->numofcandidates == candlist->length-1) {
                     candlist = resize_candidateVector(candlist, 2*(candlist->length));
                     if (candlist->data==NULL) {
                        fprintf(stderr,"%s: resize_candidateVector() failed.\n", __func__);
                        XLAL_ERROR_VOID(XLAL_EFUNC);
                     }
                  }
                  loadCandidateData(&candlist->data[candlist->numofcandidates], fsig, per0, B, 0.0, 0.0, ihsmaxima->maxima->data[locationinmaximastruct], h0, 0.0, 0, ffdata->tfnormalization);
                  (candlist->numofcandidates)++;
               } /* if no lines are interfering */
            } /* if fom is below or equal to threshold fom */
         } /* if val exceeds threshold */
      } /* for jj < numfbins-(ii-1) */
      
      /* if (highestvalloc != -1) {
         INT4 loc = ihsmaxima->locations->data[highestvalloc];
         //Candidate frequency
         fsig = params->fmin + (0.5*(ii-1) + jjloc)/params->Tcoh;
         //Candidate modulation depth
         B = 0.5*(ii-1)/params->Tcoh;
         //Candidate period
         per0 = params->Tobs/loc;
         //Candidate h0
         //REAL8 h0 = ihs2h0_withNoiseSubtraction(ihsmaxima->maxima->data[highestvalloc], loc, jjloc, ii, params, aveNoise, fbinavgs);
         REAL8 h0 = ihs2h0(2.0*highestval, params);  //Need factor of 2 for the degrees of freedom counting
         //REAL8 significance = 0.0;
         //significance = log10(significance_of_IHSval(ihsmaxima->maxima->data[highestvalloc], loc, jjloc, ii, params, aveNoise, fbinavgs));
         
         if (candlist->numofcandidates == candlist->length-1) {
            candlist = resize_candidateVector(candlist, 2*(candlist->length));
            if (candlist->data==NULL) {
               fprintf(stderr,"%s: resize_candidateVector() failed.\n", __func__);
               XLAL_ERROR_VOID(XLAL_EFUNC);
            }
         }
         //loadCandidateData(&candlist->data[candlist->numofcandidates], fsig, per0, B, 0.0, 0.0, ihsmaxima->maxima->data[locationinmaximastruct], h0, 0.0, 0, sqrt(ffdata->tfnormalization/2.0*params->Tcoh));
         loadCandidateData(&candlist->data[candlist->numofcandidates], fsig, per0, B, 0.0, 0.0, ihsmaxima->maxima->data[highestvalloc], h0, 0.0, 0, ffdata->tfnormalization);
         (candlist->numofcandidates)++;
      } */
      
      //Destroy
      XLALDestroyREAL4Vector(ihss);
      XLALDestroyREAL4Vector(avgsinrange);
      XLALDestroyINT4Vector(locs);
      ihss = NULL;
      locs = NULL;
      avgsinrange = NULL;
      
   } /* for ii < ihsfarstruct->ihsfar->length */
   
   //fclose(IHSVALSOUTPUT);
   
   fprintf(stderr,"Number of IHS vals checked = %d, number exceeding IHS threshold = %d, number passing both = %d, but lines interfere with %d\n", numberofIHSvalsChecked, numberofIHSvalsExceededThresh, numberPassingBoth, linesinterferewithnum);
   
} /* findIHScandidates() */



REAL8 ihs2h0_withNoiseSubtraction(REAL8 ihsval, INT4 location, INT4 lowestfrequencybin, INT4 rows, inputParamsStruct *params, REAL4Vector *aveNoise, REAL4Vector *fbinavgs)
{
   
   INT4 ii;
   
   REAL8 dailyharmonic = params->Tobs/(24.0*3600.0);
   REAL8 dailyharmonic2 = dailyharmonic*2.0, dailyharmonic3 = dailyharmonic*3.0, dailyharmonic4 = dailyharmonic*4.0;
   REAL8 noise = 0.0;
   for (ii=1; ii<=params->ihsfactor; ii++) {
      if (!(fabs(dailyharmonic-(REAL8)(ii*location))<=1.0 || fabs(dailyharmonic2-(REAL8)(ii*location))<=1.0 || fabs(dailyharmonic3-(REAL8)(ii*location))<=1.0 || fabs(dailyharmonic4-(REAL8)(ii*location))<=1.0)) {
         noise += aveNoise->data[ii*location];
      }
   }
   
   REAL8 totalnoise = 0.0;
   for (ii=0; ii<rows; ii++) totalnoise += noise*fbinavgs->data[lowestfrequencybin+ii];
   
   if (ihsval-totalnoise<=0.0) {
      fprintf(stderr, "%s: IHS value is less than expected noise\n", __func__);
      return 0.0;
   }
   
   REAL8 h0 = ihs2h0(2.0*ihsval-2.0*totalnoise, params);  //With 2.0 for chi-square with 2 d.o.f.
   return h0;
   
}
REAL8 ihs2h0(REAL8 ihsval, inputParamsStruct *params)
{
   
   if (ihsval<=0.0) return 0.0;
   //return 4.6*pow(ihsval/(params->Tcoh*params->Tobs),0.25);
   REAL8 prefact = 1.0;
   prefact = 7.2;
   return prefact*pow(ihsval/(params->Tcoh*params->Tobs),0.25);
   
}


REAL8 significance_of_IHSval(REAL8 ihsval, INT4 location, INT4 lowestfrequencybin, INT4 rows, inputParamsStruct *params, REAL4Vector *aveNoise, REAL4Vector *fbinavgs)
{
   
   INT4 ii;
   
   REAL8 dailyharmonic = params->Tobs/(24.0*3600.0);
   REAL8 dailyharmonic2 = dailyharmonic*2.0, dailyharmonic3 = dailyharmonic*3.0, dailyharmonic4 = dailyharmonic*4.0;
   REAL8 noise = 0.0;
   for (ii=1; ii<=params->ihsfactor; ii++) {
      if (!(fabs(dailyharmonic-(REAL8)(ii*location))<=1.0 || fabs(dailyharmonic2-(REAL8)(ii*location))<=1.0 || fabs(dailyharmonic3-(REAL8)(ii*location))<=1.0 || fabs(dailyharmonic4-(REAL8)(ii*location))<=1.0)) {
         noise += aveNoise->data[ii*location];
      }
   }
   
   REAL8 totalnoise = 0.0;
   for (ii=0; ii<rows; ii++) totalnoise += noise*fbinavgs->data[lowestfrequencybin+ii];
   
   REAL8 significance = gsl_cdf_chisq_Q(2.0*ihsval, 2.0*totalnoise);
   
   return significance;
   
}



