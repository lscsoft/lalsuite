/*
*  Copyright (C) 2010, 2011, 2012 Evan Goetz
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


//////////////////////////////////////////////////////////////
// Create vectors for IHS maxima struct
ihsMaximaStruct * new_ihsMaxima(INT4 fbins, INT4 rows)
{
   
   ihsMaximaStruct *ihsmaxima = XLALMalloc(sizeof(*ihsmaxima));
   if (ihsmaxima==NULL) {
      fprintf(stderr,"%s: XLALMalloc(%zu) failed.\n", __func__, sizeof(*ihsmaxima));
      XLAL_ERROR_NULL(XLAL_ENOMEM);
   }
   
   //The number of ihs maxima = Sum[fbins - (i-1), {i, 2, rows}]
   // = fbins*rows - fbins - (rows**2 - rows)/2
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
void runIHS(ihsMaximaStruct *output, ffdataStruct *input, ihsfarStruct *ihsfarinput, inputParamsStruct *params, INT4 rows, REAL4Vector *aveNoise, REAL4Vector *FbinMean)
{
   
   INT4 ii, jj;
   
   INT4 numfbins = input->numfbins;
   INT4 numfprbins = input->numfprbins;
   
   //Allocate memory for the necessary vectors
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

   //We want to ignore daily and sidereal harmonics, so mark the values
   REAL8 dailyharmonic = params->Tobs/(24.0*3600.0);
   REAL8 siderealharmonic = params->Tobs/86164.0905;
   REAL8 dailyharmonic2 = dailyharmonic*2.0, dailyharmonic3 = dailyharmonic*3.0, dailyharmonic4 = dailyharmonic*4.0;
   REAL8 siderealharmonic2 = siderealharmonic*2.0, siderealharmonic3 = siderealharmonic*3.0, siderealharmonic4 = siderealharmonic*4.0;
   INT4Vector *markedharmonics = XLALCreateINT4Vector(row->length);
   if (markedharmonics==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", __func__, row->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
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
      if (!params->weightedIHS) incHarmSumVector(ihsvector, row, params->ihsfactor);
      else incHarmSumVectorWeighted(ihsvector, row, aveNoise, params->ihsfactor);
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: incHarmSumVector() failed.\n", __func__);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      }
      
      //Copy the result into the ihsvector sequence
      memcpy(&(ihsvectorsequence->data[ii*ihsvector->length]), ihsvector->data, sizeof(REAL4)*ihsvector->length);
      
   } /* for ii < ihss->length */

   /* FILE *IHSSEQUENCE = fopen("./outputtemp/ihssequence.dat","w");
   for (ii=0; ii<(INT4)(ihsvectorsequence->length*ihsvectorsequence->vectorLength); ii++) {
      fprintf(IHSSEQUENCE, "%.6f\n", ihsvectorsequence->data[ii]);
   }
   fclose(IHSSEQUENCE); */

   //Now do the summing of the IHS values
   sumIHSSequence(output, ihsfarinput, ihsvectorsequence, rows, FbinMean, params);

   //Destroy stuff
   XLALDestroyREAL4VectorSequence(ihsvectorsequence);
   XLALDestroyREAL4Vector(row);
   XLALDestroyREAL4Vector(ihss);
   XLALDestroyREAL4Vector(ihsvector);
   XLALDestroyINT4Vector(locs);
   XLALDestroyINT4Vector(markedharmonics);
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
// Compute the IHS sum maximum
void incHarmSum(ihsVals *output, REAL4Vector *input, INT4 ihsfactor)
{
   
   INT4 ii;
   
   output->ihs = 0.0;
   
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
// Compute the IHS vector -- does not compute the maximum value
void incHarmSumVector(REAL4Vector *output, REAL4Vector *input, INT4 ihsfactor)
{

   INT4 ii, jj, highval = (INT4)floor((1.0/(REAL8)ihsfactor)*input->length);
   
   //From 5 up to the highval
   for (ii=5; ii<highval; ii++) {
      output->data[ii-5] = 0.0;  //first set the value to zero
      for (jj=1; jj<=ihsfactor; jj++) output->data[ii-5] += input->data[ii*jj];

      //check that we're not going outside of the allowed limits
      if (ihsfactor*(ii-1)>(INT4)input->length-1) {
         fprintf(stderr, "%s: final point exceeds the allowed limits!\n", __func__);
         XLAL_ERROR_VOID(XLAL_EBADLEN);
      } /* check within limits */
   } /* for ii=5 --> highval */

} /* incHarmSumVector() */


void incHarmSumVectorWeighted(REAL4Vector *output, REAL4Vector *input, REAL4Vector *aveNoise, INT4 ihsfactor)
{

   INT4 ii, jj, highval = (INT4)floor((1.0/(REAL8)ihsfactor)*input->length);

   //From 5 up to the highval
   for (ii=5; ii<highval; ii++) {
      output->data[ii-5] = 0.0;  //first set the value to zero
      for (jj=1; jj<=ihsfactor; jj++) {
         REAL4 weight = 1.0/(aveNoise->data[ii*jj]*aveNoise->data[ii*jj]);
         output->data[ii-5] += input->data[ii*jj]*weight;
      }

      //check that we're not going outside of the allowed limits
      if (ihsfactor*(ii-1)>(INT4)input->length-1) {
         fprintf(stderr, "%s: final point exceeds the allowed limits!\n", __func__);
         XLAL_ERROR_VOID(XLAL_EBADLEN);
      } /* check within limits */
   } /* for ii=5 --> highval */

} /* incHarmSumVectorWeighted() */


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
   //comment this
   //trials = 20000;
   
   //Allocations for IHS values for the number of trials
   REAL4Vector *noise = XLALCreateREAL4Vector(aveNoise->length);
   REAL4Vector *ihsvector = XLALCreateREAL4Vector((INT4)floor((1.0/(REAL8)params->ihsfactor)*aveNoise->length)-5);
   if (noise==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, aveNoise->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (ihsvector==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, (INT4)floor((1.0/(REAL8)params->ihsfactor)*aveNoise->length)-5);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   REAL4VectorSequence *ihsvectorsequence = XLALCreateREAL4VectorSequence(trials, ihsvector->length);
   REAL4Vector *ihss = XLALCreateREAL4Vector(trials);
   INT4Vector *locs = XLALCreateINT4Vector(trials);
   if (ihsvectorsequence==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4VectorSequence(%d,%d) failed.\n", __func__, trials, ihsvector->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (ihss==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, trials);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (locs==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", __func__, trials);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   
   //Determine the locations of the harmonics of the earth's rotation in the IHS vector
   //Amplitude modulations caused by the varying antenna pattern can sometimes cause excess power, so we ignore these harmonics
   REAL8 dailyharmonic = Tobs/(24.0*3600.0);
   REAL8 siderealharmonic = Tobs/86164.0905;
   REAL8 dailyharmonic2 = dailyharmonic*2.0, dailyharmonic3 = dailyharmonic*3.0, dailyharmonic4 = dailyharmonic*4.0;
   REAL8 siderealharmonic2 = siderealharmonic*2.0, siderealharmonic3 = siderealharmonic*3.0, siderealharmonic4 = siderealharmonic*4.0;
   INT4Vector *markedharmonics = XLALCreateINT4Vector(aveNoise->length);
   if (markedharmonics==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", __func__, aveNoise->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   memset(markedharmonics->data, 0, sizeof(INT4)*markedharmonics->length);
   if (!params->noNotchHarmonics) {
      for (ii=0; ii<(INT4)markedharmonics->length; ii++) {
         if (fabs(dailyharmonic-(REAL8)ii)<=1.0 || fabs(dailyharmonic2-(REAL8)ii)<=1.0 || fabs(dailyharmonic3-(REAL8)ii)<=1.0 || fabs(dailyharmonic4-(REAL8)ii)<=1.0 || fabs(siderealharmonic-(REAL8)ii)<=1.0 || fabs(siderealharmonic2-(REAL8)ii)<=1.0 || fabs(siderealharmonic3-(REAL8)ii)<=1.0 || fabs(siderealharmonic4-(REAL8)ii)<=1.0) markedharmonics->data[ii] = 1;
      }
   }

   //Do the expected IHS values from the expected background here
   memcpy(noise->data, aveNoise->data, sizeof(REAL4)*aveNoise->length);
   for (ii=0; ii<(INT4)aveNoise->length; ii++) if (markedharmonics->data[ii]==1) noise->data[ii] = 0.0;
   if (!params->weightedIHS) incHarmSumVector(output->expectedIHSVector, noise, params->ihsfactor);
   else incHarmSumVectorWeighted(output->expectedIHSVector, noise, aveNoise, params->ihsfactor);
   //for (ii=0; ii<(INT4)output->expectedIHSVector->length; ii++) fprintf(stderr, "%g\n", output->expectedIHSVector->data[ii]);
   
   //Now do a number of trials
   for (ii=0; ii<trials; ii++) {
      //Make exponential noise, removing harmonics of 24 hours to match with the same method as real analysis
      for (jj=0; jj<(INT4)aveNoise->length; jj++) {
         if (markedharmonics->data[jj]==0) {
            noise->data[jj] = (REAL4)(gsl_ran_exponential(params->rng, aveNoise->data[jj]));
         } else noise->data[jj] = 0.0;
      } /* for jj < aveNoise->length */

      //Make a random number 1 +/- 0.2 to create the variations in the nosie that we typically observe
      //This number needs to be positive
      REAL8 randval = 1.0;
      randval = 1.0 + gsl_ran_gaussian(params->rng, 0.2);
      while (randval<=0.0 || randval>=2.0) randval = 1.0 + gsl_ran_gaussian(params->rng, 0.2);     //limit range of variation

      //scale the exponential noise
      if (params->useSSE) {
         sseScaleREAL4Vector(noise, noise, randval);
         if (xlalErrno!=0) {
           fprintf(stderr, "%s: sseScaleREAL4Vector() failed.\n", __func__);
           XLAL_ERROR_VOID(XLAL_EFUNC);
         }
      } else for (jj=0; jj<(INT4)noise->length; jj++) noise->data[jj] *= randval;

      //Compute IHS value on exponential noise
      if (!params->weightedIHS) incHarmSumVector(ihsvector, noise, params->ihsfactor);
      else incHarmSumVectorWeighted(ihsvector, noise, aveNoise, params->ihsfactor);
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: incHarmSumVector() failed.\n", __func__);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      }

      //Copy the result into the IHS vector sequence
      memcpy(&(ihsvectorsequence->data[ii*ihsvector->length]), ihsvector->data, sizeof(REAL4)*ihsvector->length);
   } /* for ii < trials */
   
   //Destroy stuff
   XLALDestroyREAL4Vector(noise);
   XLALDestroyREAL4Vector(ihsvector);
   XLALDestroyINT4Vector(markedharmonics);
   
   //Create a fake vector with the same average value in each bin = 1.0
   REAL4Vector *FbinMean = XLALCreateREAL4Vector(trials);
   if (FbinMean==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, trials);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   for (ii=0; ii<trials; ii++) FbinMean->data[ii] = 1.0;

   //Calculate the IHS sum values for the IHS trials
   sumIHSSequenceFAR(output, ihsvectorsequence, rows, FbinMean, params);

   //Destroy stuff
   XLALDestroyREAL4VectorSequence(ihsvectorsequence);
   XLALDestroyREAL4Vector(ihss);
   XLALDestroyREAL4Vector(FbinMean);
   XLALDestroyINT4Vector(locs);

} /* genIhsFar() */



//////////////////////////////////////////////////////////////
// Compute the IHS sums for a number of rows used for the FAR calculation
// This is a bit complicated, so read the comments through the function
void sumIHSSequenceFAR(ihsfarStruct *outputfar, REAL4VectorSequence *ihsvectorsequence, INT4 rows, REAL4Vector *FbinMean, inputParamsStruct *params)
{

   INT4 ii, jj;

   //The minimum and maximum index to search in the IHS vector
   INT4 maxIndexForIHS = (INT4)ceil(fmin(params->Tobs/7200.0, params->Tobs/params->Pmin)) - 5;
   INT4 minIndexForIHS = (INT4)floor(fmax(5.0, params->Tobs/params->Pmax)) - 5;

   //Allocate a vector sequence that holds the summed values of at least two nearest neighbor rows
   //On the first iteration this holds the nearest neighbor sums, but on subsequent iterations, this holds nearest 3 neighbors
   //sums, then nearest 4 neighbor sums, and so on.
   REAL4VectorSequence *tworows = XLALCreateREAL4VectorSequence(ihsvectorsequence->length-1, ihsvectorsequence->vectorLength);
   if (tworows==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4VectorSequence(%d,%d) failed.\n", __func__, ihsvectorsequence->length-1, ihsvectorsequence->vectorLength);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   memset(tworows->data, 0, sizeof(REAL4)*tworows->length*tworows->vectorLength);  //Set everything to 0 at the start

   //Allocate vectors of the ihs values and locations of the maximums
   REAL4Vector *ihsvalues = XLALCreateREAL4Vector(ihsvectorsequence->length);
   INT4Vector *ihslocations = XLALCreateINT4Vector(ihsvectorsequence->length);
   if (ihsvalues==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ihsvectorsequence->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (ihslocations==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", __func__, ihsvectorsequence->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }

   //Vectors for values above the noise and scaling the noise
   REAL4Vector *excessabovenoise = XLALCreateREAL4Vector(ihsvectorsequence->vectorLength);
   REAL4Vector *scaledExpectedIHSVectorValues = XLALCreateREAL4Vector(ihsvectorsequence->vectorLength);
   if (excessabovenoise==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ihsvectorsequence->vectorLength);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (scaledExpectedIHSVectorValues==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ihsvectorsequence->vectorLength);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }

   //Finding the maximum for each IHS vector and the location
   for (ii=0; ii<(INT4)ihsvalues->length; ii++) {
      //Use SSE or not
      if (params->useSSE) {
         sseScaleREAL4Vector(scaledExpectedIHSVectorValues, outputfar->expectedIHSVector, FbinMean->data[ii]); //Scale the expected IHS vector by the value in FbinMean (in this function, it is 1.0)
         if (xlalErrno!=0) {
            fprintf(stderr, "%s: sseScaleREAL4Vector() failed.\n", __func__);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
         sseSSVectorSequenceSubtract(excessabovenoise, ihsvectorsequence, scaledExpectedIHSVectorValues, ii);  //subtract the noise from the data
         if (xlalErrno!=0) {
            fprintf(stderr, "%s: sseSSVectorSequenceSubtract() failed.\n", __func__);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
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
      rowsequencelocs = XLALCreateINT4Vector(ii);
      foms = XLALCreateREAL4Vector(ihsvectorsequence->length-(ii-1));
      if (rowsequencelocs==NULL) {
         fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", __func__, ii);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      } else if (foms==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ihsvectorsequence->length-(ii-1));
         XLAL_ERROR_VOID(XLAL_EFUNC);
      }
         
      //If the user has specified that we should use SSE operations, then do the nearest neighbor summing.
      //The result is in the tworows variable
      if (params->useSSE) {
         if (ii>2) sseSSVectorSequenceSum(tworows, tworows, ihsvectorsequence, 0, ii-1, 0, (INT4)ihsvectorsequence->length-(ii-1));
         else sseSSVectorSequenceSum(tworows, ihsvectorsequence, ihsvectorsequence, 0, ii-1, 0, (INT4)ihsvectorsequence->length-(ii-1));
         if (xlalErrno!=0) {
            fprintf(stderr, "%s: sseSSVectorSequenceSum() failed.\n", __func__);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
      }
         
      //Now we are going to loop through the input ihsvectorsequence up to the number of rows-1
      for (jj=0; jj<(INT4)ihsvectorsequence->length-(ii-1); jj++) {
         //A bit tricky here: sum IHS values across SFT frequency bins into the tworows variable only if we didn't do it with SSE above
         //We do this with a fast addition function fastSSVectorSequenceSum() and with no error checking, so be careful!
         if (!params->useSSE && ii>2) fastSSVectorSequenceSum(tworows, tworows, ihsvectorsequence, jj, ii-1+jj, jj); 
         else if (!params->useSSE) fastSSVectorSequenceSum(tworows, ihsvectorsequence, ihsvectorsequence, jj, ii-1+jj, jj);

         //Compute IHS FOM value
         memcpy(rowsequencelocs->data, &(ihslocations->data[jj]), sizeof(INT4)*ii);  //First copy the necessary values
         foms->data[jj] = ihsFOM(rowsequencelocs, (INT4)outputfar->expectedIHSVector->length);  //then compute the FOM
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
         outputfar->ihsdistSigma->data[ii-2] = calcStddev(sampledtempihsvals);  //We also calculate the standard deviation
         if (xlalErrno!=0) {
            fprintf(stderr, "%s: calcStddev() failed.\n", __func__);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
      } else {
         //If there were fewer than 10000 entries, then we will keep all those that are part of the nearest neighbor
         //sum up to this point
         sampledtempihsvals = XLALCreateREAL4Vector((ihsvectorsequence->length-(ii-1))*ihsvectorsequence->vectorLength);
         if (sampledtempihsvals==NULL) {
            fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, (ihsvectorsequence->length-(ii-1))*ihsvectorsequence->vectorLength);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
         memcpy(sampledtempihsvals->data, tworows->data, sizeof(REAL4)*sampledtempihsvals->length);
         outputfar->ihsdistMean->data[ii-2] = calcMean_ignoreZeros(sampledtempihsvals);  //Calculate the mean value (remember, don't accept zeros)
         outputfar->ihsdistSigma->data[ii-2] = calcStddev_ignoreZeros(sampledtempihsvals);  //We also calculate the standard deviation
         if (xlalErrno!=0) {
            fprintf(stderr, "%s: calcStddev_ignoreZeros() failed.\n", __func__);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
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
               if (xlalErrno!=0) {
                  fprintf(stderr, "%s: gsl_cdf_chisq_Qinv() failed.\n", __func__);
                  XLAL_ERROR_VOID(XLAL_EFUNC);
               }
            } else if (params->fastchisqinv && sampledtempihsvals->data[jj]!=0.0) {
               numvals++;
               //farave += cdf_chisq_Qinv(params->ihsfar, 0.5*sampledtempihsvals->data[jj]) + 0.5*sampledtempihsvals->data[jj];  //Old
               farave += 0.5*cdf_chisq_Qinv(params->ihsfar, 2.0*sampledtempihsvals->data[jj]);
               if (xlalErrno!=0) {
                  fprintf(stderr, "%s: cdf_chisq_Qinv() failed.\n", __func__);
                  XLAL_ERROR_VOID(XLAL_EFUNC);
               }
            } //fastchisqinv?
         } //for jj=0 --> sampledtempihsval->length
         outputfar->ihsfar->data[ii-2] = farave/(REAL8)numvals;
      } //if params->ihsfar != 1.0

      //Comment this out
      /* if (ii==360) {
         FILE *IHSOUTPUT = fopen("./output/ihs2rows.dat","w");
         for (jj=0; jj<(INT4)sampledtempihsvals->length; jj++) fprintf(IHSOUTPUT, "%g\n", sampledtempihsvals->data[jj]);
         fclose(IHSOUTPUT);
         } */

      XLALDestroyREAL4Vector(sampledtempihsvals);
         
      //FOM part
      outputfar->ihsfomdistMean->data[ii-2] = calcMean(foms);
      outputfar->ihsfomdistSigma->data[ii-2] = calcStddev(foms);
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: calcStddev() failed.\n", __func__);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      }
      //We need to find the smallest values
      if (params->ihsfomfar!=1.0 && params->ihsfom==0.0) {
         REAL4Vector *smallestfomvals = XLALCreateREAL4Vector((INT4)round((ihsvalues->length-ii+1)*params->ihsfomfar)+1);
         if (smallestfomvals==NULL) {
            fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, (INT4)round((ihsvalues->length-ii)*params->ihsfomfar)+1);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
         sort_float_smallest(smallestfomvals, foms);  //Sort, finding the smallest values
         if (xlalErrno!=0) {
            fprintf(stderr, "%s: sort_float_smallest() failed.\n", __func__);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
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
   
} /*sumIHSSequenceFAR() */


// We are going to find the nearest neighbor sums from some minimum up to a maximum given by rows
// In the function we will select the the location which is the **maximum above the noise**
// The sumIHSSequenceFAR() function is based on this one, and is complicated, so read comments carefully
void sumIHSSequence(ihsMaximaStruct *output, ihsfarStruct *inputfar, REAL4VectorSequence *ihsvectorsequence, INT4 rows, REAL4Vector *FbinMean, inputParamsStruct *params)
{

   INT4 ii, jj, kk;

   //Again, we start off by allocating a "tworows" vector sequence of IHS nearest neighbor sums
   REAL4VectorSequence *tworows = XLALCreateREAL4VectorSequence(ihsvectorsequence->length-1, ihsvectorsequence->vectorLength);
   if (tworows==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4VectorSequence(%d,%d) failed.\n", __func__, ihsvectorsequence->length-1, ihsvectorsequence->vectorLength);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   memset(tworows->data, 0, sizeof(REAL4)*tworows->length*tworows->vectorLength);      //Set everything to 0.0
   
   //Allocation of ihs values and locations
   REAL4Vector *ihsvalues = XLALCreateREAL4Vector(ihsvectorsequence->length);
   INT4Vector *ihslocations = XLALCreateINT4Vector(ihsvectorsequence->length);
   if (ihsvalues==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ihsvectorsequence->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (ihslocations==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", __func__, ihsvectorsequence->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   
   //Vectors for values above the noise and scaling the noise
   REAL4Vector *excessabovenoise = XLALCreateREAL4Vector(ihsvectorsequence->vectorLength);
   REAL4Vector *scaledExpectedIHSVectorValues = XLALCreateREAL4Vector(ihsvectorsequence->vectorLength);
   if (excessabovenoise==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ihsvectorsequence->vectorLength);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (scaledExpectedIHSVectorValues==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, ihsvectorsequence->vectorLength);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   
   //The minimum and maximum index to search in the IHS vector
   INT4 maxIndexForIHS = (INT4)ceil(fmin(params->Tobs/7200.0, params->Tobs/params->Pmin)) - 5;
   INT4 minIndexForIHS = (INT4)floor(fmax(5.0, params->Tobs/params->Pmax)) - 5;
   
   //Finding the maximum for each IHS vector and the location using SSE or not
   for (ii=0; ii<(INT4)ihsvalues->length; ii++) {
      if (params->useSSE) {
         sseScaleREAL4Vector(scaledExpectedIHSVectorValues, inputfar->expectedIHSVector, FbinMean->data[ii]);  //Scale the expected IHS vector by the FbinMean data value
         if (xlalErrno!=0) {
            fprintf(stderr, "%s: sseScaleREAL4Vector() failed.\n", __func__);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
         sseSSVectorSequenceSubtract(excessabovenoise, ihsvectorsequence, scaledExpectedIHSVectorValues, ii);  //subtract the noise from the data
         if (xlalErrno!=0) {
            fprintf(stderr, "%s: sseSSVectorSequenceSubtract() failed.\n", __func__);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
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
         rowsequencelocs = XLALCreateINT4Vector(ii);
         if (rowsequencelocs==NULL) {
            fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", __func__, ii);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
         
         //The maximum index to search in the IHS vector
         //maxIndexForIHS = (INT4)ceil(fmin(2.0*params->Tobs/minPeriod(0.5*(ii-1)/params->Tcoh, params->Tcoh), 2.0*params->Tobs/7200.0)) - 5;
         maxIndexForIHS = (INT4)ceil(fmin(params->Tobs/params->Pmin, fmin(params->Tobs/minPeriod(0.5*(ii-1)/params->Tcoh, params->Tcoh), params->Tobs/7200.0))) - 5;
         
         REAL4 sumofnoise = 0.0;    //To scale the expected IHS background
         INT4 endloc = ((ii-1)*(ii-1)-(ii-1))/2;
         
         //Sum up the IHS vectors using SSE functions
         if (params->useSSE) {
            if (ii>2) sseSSVectorSequenceSum(tworows, tworows, ihsvectorsequence, 0, ii-1, 0, (INT4)ihsvectorsequence->length-(ii-1));
            else sseSSVectorSequenceSum(tworows, ihsvectorsequence, ihsvectorsequence, 0, ii-1, 0, (INT4)ihsvectorsequence->length-(ii-1));
            if (xlalErrno!=0) {
               fprintf(stderr, "%s: sseSSVectorSequenceSum() failed.\n", __func__);
               XLAL_ERROR_VOID(XLAL_EFUNC);
            }
         } /* use SSE code */

         //Loop through the IHS vector neighbor sums
         for (jj=0; jj<(INT4)ihsvectorsequence->length-(ii-1); jj++) {
            //If we didn't use SSE to sum the vector sequence (see lines above)
            if (!params->useSSE) {
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
         } /* for jj< ihsvectorsequence->length - (ii-1) */

         /* if (ii==360) {
            FILE *IHSOUTPUT = fopen("./output/ihs2rows_real.dat","w");
            for (jj=0; jj<(INT4)tworows->vectorLength*(ihsvectorsequence->length-(ii-1)); jj++) fprintf(IHSOUTPUT, "%g\n", tworows->data[jj]);
            fclose(IHSOUTPUT);
            } */

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



//Sum a specific vector to another specific vector in two vector sequences
//not as quick as the fastSSVectorSequenceSum() function or the sseSSVectorSequenceSum function()
void SSVectorSequenceSum(REAL4VectorSequence *output, REAL4VectorSequence *input1, REAL4VectorSequence *input2, INT4 vectorpos1, INT4 vectorpos2, INT4 outputvectorpos)
{
   INT4 ii, vec1 = vectorpos1*input1->vectorLength, vec2 = vectorpos2*input2->vectorLength, outvec = outputvectorpos*output->vectorLength;
   for (ii=0; ii<(INT4)input1->vectorLength; ii++) output->data[outvec + ii] = input1->data[vec1 + ii] + input2->data[vec2 + ii];
}


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
   for (ii=0; ii<(INT4)(input->vectorLength*(input->length-(rows-1))); ii++) output->data[ii] = (REAL4)tworows->data[ii];
   XLALDestroyREAL8VectorSequence(tworows);

   return output;

}


//////////////////////////////////////////////////////////////
// Calculate the IHS FOM for a number of rows
REAL4 ihsFOM(INT4Vector *locs, INT4 fomnorm)
{
   
   INT4 ii;
   REAL4 fom = 0.0;

   for (ii=0; ii<(INT4)(locs->length*0.5); ii++) {
      fom += (REAL4)((locs->data[ii]-locs->data[locs->length-ii-1])*(locs->data[ii]-locs->data[locs->length-ii-1]))/(fomnorm*fomnorm);
   }

   return fom;

} /* ihsFOM() */




//Finds IHS candidates above thresholds
void findIHScandidates(candidateVector *candlist, ihsfarStruct *ihsfarstruct, inputParamsStruct *params, ffdataStruct *ffdata, ihsMaximaStruct *ihsmaxima, REAL4Vector *fbinavgs, REAL4VectorSequence *trackedlines)
{

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
      
      REAL8 highestval = 0.0, highestsignificance = 0.0; //highestvalnoise = 0.0
      INT4 highestvalloc = -1, jjloc = 0;
      for (jj=0; jj<(INT4)numfbins-(ii-1); jj++) {
      
         //Noise in the range of the rows, mean for IHS
         memcpy(avgsinrange->data, &(fbinavgs->data[jj]), sizeof(REAL4)*ii);
         REAL4 meanNoise = calcMean(avgsinrange);
         if (XLAL_IS_REAL4_FAIL_NAN(meanNoise)) {
            fprintf(stderr,"%s: calcMean() failed.\n", __func__);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
         
         numberofIHSvalsChecked++;  //count the number of ihs values checked
         
         INT4 locationinmaximastruct = (ii-2)*numfbins-((ii-1)*(ii-1)-(ii-1))/2+jj;  //the location in the ihsmaximastruct
         
         //Check the IHS sum against the FAR (scaling FAR with mean of the noise in the range of rows)
         if (ihsmaxima->maxima->data[locationinmaximastruct] > ihsfarstruct->ihsfar->data[ii-2]*meanNoise) {
            
	    numberofIHSvalsExceededThresh++;  //count the number of values exceeding the IHS value threshold
            
            if (ihsfarstruct->fomfarthresh->data[ii-2]==-1.0 || ihsmaxima->foms->data[locationinmaximastruct]<=ihsfarstruct->fomfarthresh->data[ii-2]) {
               
	      numberPassingBoth++;  //count number passing both IHS value and FOM value thresholds
               //fprintf(IHSVALSOUTPUT, "%.6f %.6f\n", 0.5*(ii-1)/params->Tcoh, ihsmaxima->maxima->data[locationinmaximastruct]);
               
               INT4 loc = ihsmaxima->locations->data[locationinmaximastruct];
               per0 = params->Tobs/loc;                                          //Candidate period
               fsig = params->fmin - params->dfmax + ((0.5*(ii-1) + jj) - 6.0)/params->Tcoh;             //Candidate frequency
               B = 0.5*(ii-1)/params->Tcoh;                                      //Candidate modulation depth
               
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
                  //if (ii==2) fprintf(stderr, "%g %g\n", meanNoise, calcRms(avgsinrange));     //remove this
                  
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
         //Candidate frequency
         fsig = params->fmin - params->dfmax + (0.5*(ii-1) + jjloc - 6.0)/params->Tcoh;
         //Candidate modulation depth
         B = 0.5*(ii-1)/params->Tcoh;
         //Candidate period
         per0 = params->Tobs/loc;
         //Candidate h0
         //REAL8 h0 = ihs2h0_withNoiseSubtraction(ihsmaxima->maxima->data[highestvalloc], loc, jjloc, ii, params, aveNoise, fbinavgs);
         REAL8 h0 = ihs2h0(2.0*highestval, params);  //Need factor of 2 for the degrees of freedom counting
         
         if (candlist->numofcandidates == candlist->length-1) {
            candlist = resize_candidateVector(candlist, 2*(candlist->length));
            if (candlist->data==NULL) {
               fprintf(stderr,"%s: resize_candidateVector() failed.\n", __func__);
               XLAL_ERROR_VOID(XLAL_EFUNC);
            }
         }
         loadCandidateData(&candlist->data[candlist->numofcandidates], fsig, per0, B, 0.0, 0.0, ihsmaxima->maxima->data[highestvalloc], h0, highestsignificance, 0, ffdata->tfnormalization);
         (candlist->numofcandidates)++;
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
         //REAL8 h0 = ihs2h0_withNoiseSubtraction(ihsmaxima->maxima->data[highestvalloc], loc, jjloc, ii, params, aveNoise, fbinavgs);
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
   
} /* findIHScandidates() */



REAL8 ihs2h0_withNoiseSubtraction(REAL8 ihsval, INT4 location, INT4 lowestfrequencybin, INT4 rows, inputParamsStruct *params, REAL4Vector *aveNoise, REAL4Vector *fbinavgs)
{
   
   INT4 ii;
   
   REAL8 dailyharmonic = params->Tobs/(24.0*3600.0);
   REAL8 siderealharmonic = params->Tobs/86164.0905;
   REAL8 dailyharmonic2 = dailyharmonic*2.0, dailyharmonic3 = dailyharmonic*3.0, dailyharmonic4 = dailyharmonic*4.0;
   REAL8 siderealharmonic2 = siderealharmonic*2.0, siderealharmonic3 = siderealharmonic*3.0, siderealharmonic4 = siderealharmonic*4.0;
   REAL8 noise = 0.0;
   for (ii=1; ii<=params->ihsfactor; ii++) {
      if (!params->noNotchHarmonics || !(fabs(dailyharmonic-(REAL8)(ii*location))<=1.0 || fabs(dailyharmonic2-(REAL8)(ii*location))<=1.0 || fabs(dailyharmonic3-(REAL8)(ii*location))<=1.0 || fabs(dailyharmonic4-(REAL8)(ii*location))<=1.0 || fabs(siderealharmonic-(REAL8)(ii*location))<=1.0 || fabs(siderealharmonic2-(REAL8)(ii*location))<=1.0 || fabs(siderealharmonic3-(REAL8)(ii*location))<=1.0 || fabs(siderealharmonic4-(REAL8)(ii*location))<=1.0)) {
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
   REAL8 prefact = 1.0;
   prefact = 7.2;
   return prefact*pow(ihsval/(params->Tcoh*params->Tobs),0.25);
   
}




