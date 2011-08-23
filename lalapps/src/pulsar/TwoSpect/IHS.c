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

#include <lal/LALMalloc.h>
#include <lal/SeqFactories.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "IHS.h"
#include "statistics.h"
#include "candidates.h"

//////////////////////////////////////////////////////////////
// Create vectors for IHS maxima struct
ihsMaximaStruct * new_ihsMaxima(INT4 fbins, INT4 rows)
{
   
   const char *fn = __func__;
   
   ihsMaximaStruct *ihsmaxima = XLALMalloc(sizeof(*ihsmaxima));
   if (ihsmaxima==NULL) {
      fprintf(stderr,"%s: XLALMalloc(%zu) failed.\n", fn, sizeof(*ihsmaxima));
      XLAL_ERROR_NULL(fn, XLAL_ENOMEM);
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
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numberofmaxima);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   } else if (ihsmaxima->locations==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", fn, numberofmaxima);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   } else if (ihsmaxima->foms==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numberofmaxima);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   } else if (ihsmaxima->maximaForEachFbin==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, fbins);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   } else if (ihsmaxima->locationsForEachFbin==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", fn, fbins);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
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
void runIHS(ihsMaximaStruct *output, ffdataStruct *input, inputParamsStruct *params, INT4 rows, REAL4Vector *FbinMean)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii, jj;
   
   INT4 numfbins = input->numfbins;
   INT4 numfprbins = input->numfprbins;
   
   REAL4Vector *row = XLALCreateREAL4Vector(numfprbins);
   REAL4Vector *ihss = XLALCreateREAL4Vector(numfbins);
   REAL4Vector *ihsvector = XLALCreateREAL4Vector((INT4)floor((1.0/(REAL8)params->ihsfactor)*numfprbins)-5);
   INT4Vector *locs = XLALCreateINT4Vector(numfbins);
   ihsVals *ihsvals = new_ihsVals();
   if (row==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numfprbins);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (ihss==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numfbins);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (ihsvector==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, (INT4)floor((1.0/(REAL8)params->ihsfactor)*numfprbins)-5);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (locs==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", fn, numfbins);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (ihsvals==NULL) {
      fprintf(stderr,"%s: new_ihsVals() failed.\n", fn);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   
   REAL4VectorSequence *ihsvectorsequence = XLALCreateREAL4VectorSequence(numfbins, ihsvector->length);
   if (ihsvectorsequence==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4VectorSequence(%d,%d) failed.\n", fn, numfbins, ihsvector->length);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   
   //Loop through the rows, 1 frequency at a time
   REAL8 dailyharmonic = params->Tobs/(24.0*3600.0);
   REAL8 dailyharmonic2 = dailyharmonic*2.0, dailyharmonic3 = dailyharmonic*3.0, dailyharmonic4 = dailyharmonic*4.0;
   for (ii=0; ii<(INT4)ihss->length; ii++) {
   
      //For each row, populate it with the data for that frequency bin, excluding harmonics of antenna pattern modulation
      memcpy(row->data, &(input->ffdata->data[ii*numfprbins]), sizeof(REAL4)*row->length);
      for (jj=0; jj<(INT4)row->length; jj++) {
         if (fabs(dailyharmonic-(REAL8)jj)<=1.0 || fabs(dailyharmonic2-(REAL8)jj)<=1.0 || fabs(dailyharmonic3-(REAL8)jj)<=1.0 || fabs(dailyharmonic4-(REAL8)jj)<=1.0) row->data[jj] = 0.0;
      }
      
      //Run the IHS algorithm on the row
      incHarmSumVector(ihsvector, row, params->ihsfactor);
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: incHarmSumVector() failed.\n", fn);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
      locs->data[ii] = max_index(ihsvector) + 5;
      ihss->data[ii] = ihsvector->data[locs->data[ii]-5];
      memcpy(&(ihsvectorsequence->data[ii*ihsvector->length]), ihsvector->data, sizeof(REAL4)*ihsvector->length);
      
   } /* for ii < ihss->length */
   
   //REAL4VectorSequence *ihssummedvals = ihsVectorSums(ihsvectorsequence, 360);
   //FILE *IHSSUMVALS = fopen("./output/ihssumvals_mfd.dat","w");
   //for (ii=0; ii<(INT4)ihssummedvals->length; ii++) for (jj=0; jj<(INT4)ihssummedvals->vectorLength; jj++) fprintf(IHSSUMVALS,"%.6f\n", ihssummedvals->data[ii*ihssummedvals->vectorLength + jj]);
   //for (ii=0; ii<(INT4)ihsvectorsequence->length; ii++) for (jj=0; jj<(INT4)ihsvectorsequence->vectorLength; jj++) fprintf(IHSSUMVALS,"%.6f\n", ihsvectorsequence->data[ii*ihsvectorsequence->vectorLength + jj]);
   //fclose(IHSSUMVALS);
   //XLALDestroyREAL4VectorSequence(ihssummedvals);
   ihsSums2_withFAR(output, NULL, ihsvectorsequence, ihss, locs, NULL, rows, FbinMean, input->numfprbins, params, 0);
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
   
   const char *fn = __func__;
   
   ihsVals *ihsvals = XLALMalloc(sizeof(*ihsvals));
   if (ihsvals==NULL) {
      fprintf(stderr,"%s: XLALMalloc(%zu) failed.\n", fn, sizeof(*ihsvals));
      XLAL_ERROR_NULL(fn, XLAL_ENOMEM);
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
   
   const CHAR *fn = __func__;
   
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
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", fn, (INT4)floor((1.0/(REAL8)ihsfactor)*input->length)-5);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   incHarmSumVector(tempvect, input, ihsfactor);
   if (xlalErrno!=0) {
      fprintf(stderr, "%s: incHarmSumVector() failed.\n", fn);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
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
   
   const CHAR *fn = __func__;
   
   INT4 ii, jj, highval = (INT4)floor((1.0/(REAL8)ihsfactor)*input->length);
   
   for (ii=5; ii<highval; ii++) {
      output->data[ii-5] = 0.0;
      for (jj=1; jj<=ihsfactor; jj++) output->data[ii-5] += input->data[ii*jj];
      if (ihsfactor*(ii-1)>(INT4)input->length-1) {
         fprintf(stderr, "%s: final point exceeds the allowed limits!\n", fn);
         XLAL_ERROR_VOID(fn, XLAL_EBADLEN);
      }
   }
   
} /* incHarmSumVector() */


//////////////////////////////////////////////////////////////
// Allocate memory for ihsfarStruct struct
ihsfarStruct * new_ihsfarStruct(INT4 rows)
{
   
   const char *fn = __func__;
   
   ihsfarStruct *ihsfarstruct = XLALMalloc(sizeof(*ihsfarstruct));
   if (ihsfarstruct == NULL) {
      fprintf(stderr,"%s: XLALMalloc(%zu) failed.\n", fn, sizeof(*ihsfarstruct));
      XLAL_ERROR_NULL(fn, XLAL_ENOMEM);
   }
   
   ihsfarstruct->ihsfar = XLALCreateREAL4Vector(rows-1);
   ihsfarstruct->ihsdistMean = XLALCreateREAL4Vector(rows-1);
   ihsfarstruct->ihsdistSigma = XLALCreateREAL4Vector(rows-1);
   ihsfarstruct->fomfarthresh = XLALCreateREAL4Vector(rows-1);
   ihsfarstruct->ihsfomdistMean = XLALCreateREAL4Vector(rows-1);
   ihsfarstruct->ihsfomdistSigma = XLALCreateREAL4Vector(rows-1);
   if (ihsfarstruct->ihsfar==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, rows-1);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   } else if (ihsfarstruct->ihsdistMean==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, rows-1);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   } else if( ihsfarstruct->ihsdistSigma==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, rows-1);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   } else if( ihsfarstruct->fomfarthresh==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, rows-1);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   } else if( ihsfarstruct->ihsfomdistMean==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, rows-1);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   } else if( ihsfarstruct->ihsfomdistSigma==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, rows-1);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   }
   
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
   XLALFree((ihsfarStruct*)ihsfarstruct);

} /* free_ihsfarStruct() */


//////////////////////////////////////////////////////////////
// Compute the IHS FAR for a sum of a number of rows
void genIhsFar(ihsfarStruct *output, inputParamsStruct *params, INT4 rows, REAL4Vector *aveNoise)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii, jj;
   REAL8 Tobs = params->Tobs;
   
   INT4 trials = (INT4)round(1.0*0.000001/params->ihsfar);    //Number of trials to determine FAR value
   if (params->ihsfomfar!=0.0 && trials<(INT4)round(1.0*0.000001/params->ihsfomfar)) {
      trials = (INT4)round(1.0*0.000001/params->ihsfomfar);
   }
   trials += rows;
   
   //Initialize random number generator
   gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
   if (rng==NULL) {
      fprintf(stderr,"%s: gsl_rng_alloc() failed.\n", fn);
      XLAL_ERROR_VOID(fn, XLAL_ENOMEM);
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
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, aveNoise->length);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (ihsvector==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", fn, (INT4)floor((1.0/(REAL8)params->ihsfactor)*aveNoise->length)-5);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (ihsvectorsequence==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4VectorSequence(%d,%d) failed.\n", fn, trials, ihsvector->length);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (ihss==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, trials);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (locs==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", fn, trials);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   REAL8 singleIHSsigma = sqrt(5.0*0.05*0.05);
   REAL8 dailyharmonic = Tobs/(24.0*3600.0);
   REAL8 dailyharmonic2 = dailyharmonic*2.0, dailyharmonic3 = dailyharmonic*3.0, dailyharmonic4 = dailyharmonic*4.0;
   INT4Vector *markedharmonics = XLALCreateINT4Vector(aveNoise->length);
   for (ii=0; ii<(INT4)markedharmonics->length; ii++) {
      if (fabs(dailyharmonic-(REAL8)ii)<=1.0 || fabs(dailyharmonic2-(REAL8)ii)<=1.0 || fabs(dailyharmonic3-(REAL8)ii)<=1.0 || fabs(dailyharmonic4-(REAL8)ii)<=1.0) markedharmonics->data[ii] = 1;
      else markedharmonics->data[ii] = 0;
   }
   for (ii=0; ii<trials; ii++) {
      REAL8 randval = 1.0 + 2.0*gsl_ran_gaussian(rng, singleIHSsigma);
      while (randval<0.0) randval = 1.0 + 2.0*gsl_ran_gaussian(rng, singleIHSsigma);
      //Make exponential noise removing harmonics of 24 hours to match with the same method as real analysis
      for (jj=0; jj<(INT4)aveNoise->length; jj++) {
         if (markedharmonics->data[jj]==0) noise->data[jj] = (REAL4)(gsl_ran_exponential(rng, aveNoise->data[jj])*randval);
         else noise->data[jj] = 0.0;
      } /* for jj < aveNoise->length */
      
      //Compute IHS value on exponential noise
      incHarmSumVector(ihsvector, noise, params->ihsfactor);
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: incHarmSumVector() failed.\n", fn);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
      locs->data[ii] = max_index(ihsvector) + 5;
      ihss->data[ii] = (REAL4)ihsvector->data[locs->data[ii]-5];
      memcpy(&(ihsvectorsequence->data[ii*ihsvector->length]), ihsvector->data, sizeof(REAL4)*ihsvector->length);
   } /* for ii < trials */
   XLALDestroyREAL4Vector(noise);
   XLALDestroyREAL4Vector(ihsvector);
   XLALDestroyINT4Vector(markedharmonics);
   gsl_rng_free(rng);
   
   
   //Calculate the IHS sum values for the IHS trials
   ihsMaximaStruct *ihsmaxima = new_ihsMaxima(trials, rows);
   if (ihsmaxima==NULL) {
      fprintf(stderr,"%s: new_ihsMaxima(%d, %d) failed.\n", fn, trials, rows);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   //Create a fake vector with the same average value in each bin
   REAL4Vector *FbinMean = XLALCreateREAL4Vector(trials);
   if (FbinMean==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", fn, trials);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   for (ii=0; ii<trials; ii++) FbinMean->data[ii] = 1.0;
   
   ihsSums2_withFAR(ihsmaxima, output, ihsvectorsequence, ihss, locs, aveNoise, rows, FbinMean, (INT4)aveNoise->length, params, 1);
   XLALDestroyREAL4VectorSequence(ihsvectorsequence);
   
   //Destroy variables
   XLALDestroyREAL4Vector(ihss);
   XLALDestroyREAL4Vector(FbinMean);
   XLALDestroyINT4Vector(locs);
   free_ihsMaxima(ihsmaxima);
   

} /* genIhsFar() */



//////////////////////////////////////////////////////////////
// Compute the IHS sums for a number of rows
void ihsSums(ihsMaximaStruct *output, REAL4Vector *ihss, INT4Vector *locs, INT4 rows, REAL4Vector *FbinMean, INT4 locationnormfactor)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii, jj, locInMaximaVector;
   INT4 startPosition = 0;
   
   //Start with the vector of single row IHS values
   memcpy(output->maximaForEachFbin->data, ihss->data, sizeof(REAL4)*ihss->length);
   memcpy(output->locations->data, locs->data, sizeof(INT4)*locs->length);
   
   output->rows = rows;
   
   //Efficiently sum the IHS values
   //Sum the pairs
   for (ii=0; ii<(INT4)ihss->length-1; ii++) output->maxima->data[ii] = ihss->data[ii] + ihss->data[ii+1];
   locInMaximaVector = ii; //save the position in the maxima vector
   //Now continue summing only needing to add the next row IHS value to the recent summed value
   for (ii=3; ii<=rows; ii++) {
      for (jj=ii-1; jj<(INT4)ihss->length; jj++) {
         output->maxima->data[locInMaximaVector] = output->maxima->data[startPosition+jj-(ii-1)] + ihss->data[jj];
         locInMaximaVector++; //Increment position in maxima vector
      }
      startPosition = locInMaximaVector - (ihss->length-(ii-1)); //Position to start is the location in the maxima vector minus slots moved
   } /* for ii <= rows */
   
   INT4Vector *templocs = NULL;
   REAL4Vector *tempihss = NULL, *tempfbinmean = NULL;
   locInMaximaVector = 0;
   for (ii=2; ii<=rows; ii++) {
      templocs = XLALCreateINT4Vector(ii);
      tempihss = XLALCreateREAL4Vector(ii);
      tempfbinmean = XLALCreateREAL4Vector(ii);
      if (templocs==NULL) {
         fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", fn, ii);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      } else if (tempihss==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ii);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      } else if (tempfbinmean==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ii);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
      for (jj=0; jj<(INT4)locs->length-(ii-1); jj++) {
         memcpy(templocs->data, &(locs->data[jj]), sizeof(INT4)*templocs->length);
         memcpy(tempihss->data, &(ihss->data[jj]), sizeof(REAL4)*tempihss->length);
         memcpy(tempfbinmean->data, &(FbinMean->data[jj]), sizeof(REAL4)*tempfbinmean->length);
         output->foms->data[locInMaximaVector] = ihsFOM(tempihss, templocs, tempfbinmean, locationnormfactor);
         locInMaximaVector++;
      }
      XLALDestroyINT4Vector(templocs);
      XLALDestroyREAL4Vector(tempihss);
      XLALDestroyREAL4Vector(tempfbinmean);
      templocs = NULL;
      tempihss = NULL;
      tempfbinmean = NULL;
   }
   
} /* ihsSums() */
void ihsSums2(ihsMaximaStruct *output, REAL4VectorSequence *ihsvectorsequence, REAL4Vector *ihss, INT4Vector *locs, INT4 rows, REAL4Vector *FbinMean, INT4 locationnormfactor)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii, jj, kk;
   
   REAL4VectorSequence *tworows = XLALCreateREAL4VectorSequence(ihsvectorsequence->length-1, ihsvectorsequence->vectorLength);
   REAL4Vector *singlerow = XLALCreateREAL4Vector(ihsvectorsequence->vectorLength);
   if (tworows==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4VectorSequence(%d,%d) failed.\n", fn, ihsvectorsequence->length-1, ihsvectorsequence->vectorLength);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (singlerow==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ihsvectorsequence->vectorLength);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   REAL4Vector *rowsequencemaxima = NULL;
   REAL4Vector *fbinmeanvals = NULL;
   INT4Vector *rowsequencelocs = NULL;
   INT4 maximapos = 0;
   for (ii=1; ii<=rows; ii++) {
      if (ii==1) {
         memcpy(output->maximaForEachFbin->data, ihss->data, sizeof(REAL4)*ihss->length);
         memcpy(output->locations->data, locs->data, sizeof(INT4)*locs->length);
      } else if (ii==2) {
         rowsequencemaxima = XLALCreateREAL4Vector(ii);
         fbinmeanvals = XLALCreateREAL4Vector(ii);
         rowsequencelocs = XLALCreateINT4Vector(ii);
         if (rowsequencemaxima==NULL) {
            fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ii);
            XLAL_ERROR_VOID(fn, XLAL_EFUNC);
         } else if (fbinmeanvals==NULL) {
            fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ii);
            XLAL_ERROR_VOID(fn, XLAL_EFUNC);
         } else if (rowsequencelocs==NULL) {
            fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", fn, ii);
            XLAL_ERROR_VOID(fn, XLAL_EFUNC);
         }
         
         for (jj=0; jj<(INT4)ihsvectorsequence->length-(2-1); jj++) {
            for (kk=0; kk<(INT4)ihsvectorsequence->vectorLength; kk++) {
               tworows->data[jj*ihsvectorsequence->vectorLength + kk] = ihsvectorsequence->data[jj*ihsvectorsequence->vectorLength + kk] + ihsvectorsequence->data[(jj+1)*ihsvectorsequence->vectorLength + kk];
            }
            memcpy(singlerow->data, &(tworows->data[jj*ihsvectorsequence->vectorLength]), sizeof(REAL4)*ihsvectorsequence->vectorLength);
            //output->maxima->data[maximapos] = singlerow->data[max_index(singlerow)];
            output->maxima->data[jj] = singlerow->data[max_index(singlerow)];
            
            memcpy(rowsequencemaxima->data, &(ihss->data[jj]), sizeof(REAL4)*ii);
            memcpy(fbinmeanvals->data, &(FbinMean->data[jj]), sizeof(REAL4)*ii);
            memcpy(rowsequencelocs->data, &(locs->data[jj]), sizeof(INT4)*ii);
            output->foms->data[maximapos] = ihsFOM(rowsequencemaxima, rowsequencelocs, fbinmeanvals, locationnormfactor);
            
            maximapos++;
         }
         
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
            fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ii);
            XLAL_ERROR_VOID(fn, XLAL_EFUNC);
         } else if (fbinmeanvals==NULL) {
            fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ii);
            XLAL_ERROR_VOID(fn, XLAL_EFUNC);
         } else if (rowsequencelocs==NULL) {
            fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", fn, ii);
            XLAL_ERROR_VOID(fn, XLAL_EFUNC);
         }
         
         for (jj=0; jj<(INT4)ihsvectorsequence->length-(ii-1); jj++) {
            for (kk=0; kk<(INT4)ihsvectorsequence->vectorLength; kk++) {
               tworows->data[jj*ihsvectorsequence->vectorLength + kk] = tworows->data[jj*ihsvectorsequence->vectorLength + kk] + ihsvectorsequence->data[(ii-1+jj)*ihsvectorsequence->vectorLength + kk];
            }
            memcpy(singlerow->data, &(tworows->data[jj*ihsvectorsequence->vectorLength]), sizeof(REAL4)*ihsvectorsequence->vectorLength);
            //output->maxima->data[maximapos] = singlerow->data[max_index(singlerow)];
            output->maxima->data[(ii-2)*ihss->length-(INT4)(((ii-1)*(ii-1)-(ii-1))/2)+jj] = singlerow->data[max_index(singlerow)];
            
            memcpy(rowsequencemaxima->data, &(ihss->data[jj]), sizeof(REAL4)*ii);
            memcpy(fbinmeanvals->data, &(FbinMean->data[jj]), sizeof(REAL4)*ii);
            memcpy(rowsequencelocs->data, &(locs->data[jj]), sizeof(INT4)*ii);
            output->foms->data[maximapos] = ihsFOM(rowsequencemaxima, rowsequencelocs, fbinmeanvals, locationnormfactor);
            
            maximapos++;
         }
         
         XLALDestroyREAL4Vector(rowsequencemaxima);
         rowsequencemaxima = NULL;
         XLALDestroyREAL4Vector(fbinmeanvals);
         fbinmeanvals = NULL;
         XLALDestroyINT4Vector(rowsequencelocs);
         rowsequencelocs = NULL;
      }
   }
   XLALDestroyREAL4Vector(singlerow);
   XLALDestroyREAL4VectorSequence(tworows);
   
} /*ihsSums2() */
void ihsSums2_withFAR(ihsMaximaStruct *output, ihsfarStruct *outputfar, REAL4VectorSequence *ihsvectorsequence, REAL4Vector *ihss, INT4Vector *locs, REAL4Vector *aveNoise, INT4 rows, REAL4Vector *FbinMean, INT4 locationnormfactor, inputParamsStruct *params, INT4 calcPInvVals)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii, jj, kk;
   
   gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
   if (rng==NULL) {
      fprintf(stderr,"%s: gsl_rng_alloc() failed.\n", fn);
      XLAL_ERROR_VOID(fn, XLAL_ENOMEM);
   }
   gsl_rng_set(rng, 0);
   
   REAL4VectorSequence *tworows = XLALCreateREAL4VectorSequence(ihsvectorsequence->length-1, ihsvectorsequence->vectorLength);
   REAL4Vector *singlerow = XLALCreateREAL4Vector(ihsvectorsequence->vectorLength);
   if (tworows==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4VectorSequence(%d,%d) failed.\n", fn, ihsvectorsequence->length-1, ihsvectorsequence->vectorLength);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (singlerow==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ihsvectorsequence->vectorLength);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   REAL4VectorSequence *tworows2 = NULL;
   REAL4Vector *ihsvalsfromaveNoise = NULL, *randvals = NULL;
   if (calcPInvVals!=0) {
      tworows2 = XLALCreateREAL4VectorSequence(ihsvectorsequence->length-1, ihsvectorsequence->vectorLength);
      ihsvalsfromaveNoise = XLALCreateREAL4Vector(ihsvectorsequence->vectorLength);
      REAL4Vector *tempaveNoise = XLALCreateREAL4Vector(aveNoise->length);
      if (tworows2==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4VectorSequence(%d,%d) failed.\n", fn, ihsvectorsequence->length-1, ihsvectorsequence->vectorLength);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      } else if (ihsvalsfromaveNoise==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ihsvectorsequence->vectorLength);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      } else if (tempaveNoise==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, aveNoise->length);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
      
      REAL8 dailyharmonic = params->Tobs/(24.0*3600.0);
      REAL8 dailyharmonic2 = dailyharmonic*2.0, dailyharmonic3 = dailyharmonic*3.0, dailyharmonic4 = dailyharmonic*4.0;
      for (ii=0; ii<(INT4)aveNoise->length; ii++) {
         if (fabs(dailyharmonic-(REAL8)ii)<=1.0 || fabs(dailyharmonic2-(REAL8)ii)<=1.0 || fabs(dailyharmonic3-(REAL8)ii)<=1.0 || fabs(dailyharmonic4-(REAL8)ii)<=1.0) tempaveNoise->data[ii] = 0.0;
         else tempaveNoise->data[ii] = aveNoise->data[ii];
      }
      incHarmSumVector(ihsvalsfromaveNoise, tempaveNoise, params->ihsfactor);
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: incHarmSumVector() failed.\n", fn);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
      
      randvals = XLALCreateREAL4Vector(ihss->length);
      if (randvals==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ihss->length);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
      REAL8 sigma = sqrt(5.0*0.05*0.05);
      for (ii=0; ii<(INT4)ihss->length; ii++) {
         randvals->data[ii] = (REAL4)(1.0 + 2.0*gsl_ran_gaussian(rng, sigma));
         while (randvals->data[ii]<0.0) randvals->data[ii] = (REAL4)(1.0 + 2.0*gsl_ran_gaussian(rng, sigma));
      }
      
      XLALDestroyREAL4Vector(tempaveNoise);
   } /* if calcPInvVals != 0 */
   
   REAL4Vector *rowsequencemaxima = NULL;
   REAL4Vector *fbinmeanvals = NULL;
   INT4Vector *rowsequencelocs = NULL;
   for (ii=1; ii<=rows; ii++) {
      if (ii==1) {
         memcpy(output->maximaForEachFbin->data, ihss->data, sizeof(REAL4)*ihss->length);
         memcpy(output->locationsForEachFbin->data, locs->data, sizeof(INT4)*locs->length);
         //fprintf(stderr, "Finished row %d\n", ii);
      } else if (ii==2) {
         rowsequencemaxima = XLALCreateREAL4Vector(ii);
         fbinmeanvals = XLALCreateREAL4Vector(ii);
         rowsequencelocs = XLALCreateINT4Vector(ii);
         if (rowsequencemaxima==NULL) {
            fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ii);
            XLAL_ERROR_VOID(fn, XLAL_EFUNC);
         } else if (fbinmeanvals==NULL) {
            fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ii);
            XLAL_ERROR_VOID(fn, XLAL_EFUNC);
         } else if (rowsequencelocs==NULL) {
            fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", fn, ii);
            XLAL_ERROR_VOID(fn, XLAL_EFUNC);
         }
         
         for (jj=0; jj<(INT4)ihsvectorsequence->length-(ii-1); jj++) {
            for (kk=0; kk<(INT4)ihsvectorsequence->vectorLength; kk++) tworows->data[jj*ihsvectorsequence->vectorLength + kk] = ihsvectorsequence->data[jj*ihsvectorsequence->vectorLength + kk] + ihsvectorsequence->data[(jj+1)*ihsvectorsequence->vectorLength + kk];
            memcpy(singlerow->data, &(tworows->data[jj*ihsvectorsequence->vectorLength]), sizeof(REAL4)*ihsvectorsequence->vectorLength);
            output->locations->data[jj] = max_index(singlerow) + 5;
            output->maxima->data[jj] = singlerow->data[output->locations->data[jj]-5];
            
            if (calcPInvVals!=0) for (kk=0; kk<(INT4)ihsvalsfromaveNoise->length; kk++) tworows2->data[jj*ihsvalsfromaveNoise->length + kk] = ihsvalsfromaveNoise->data[kk]*(randvals->data[jj] + randvals->data[jj+1]);
            
            memcpy(rowsequencemaxima->data, &(ihss->data[jj]), sizeof(REAL4)*ii);
            memcpy(fbinmeanvals->data, &(FbinMean->data[jj]), sizeof(REAL4)*ii);
            memcpy(rowsequencelocs->data, &(locs->data[jj]), sizeof(INT4)*ii);
            output->foms->data[jj] = ihsFOM(rowsequencemaxima, rowsequencelocs, fbinmeanvals, locationnormfactor);
         } /* for jj < ihsvectorsequence->length-(ii-1) */
         
         if (calcPInvVals!=0) {
            REAL4Vector *tempihsvals = XLALCreateREAL4Vector((ihsvectorsequence->length-(ii-1))*ihsvectorsequence->vectorLength);
            if (tempihsvals==NULL) {
               fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, (ihsvectorsequence->length-(ii-1))*ihsvectorsequence->vectorLength);
               XLAL_ERROR_VOID(fn, XLAL_EFUNC);
            }
            memcpy(tempihsvals->data, tworows->data, sizeof(REAL4)*tempihsvals->length);
            REAL4Vector *sampledtempihsvals = NULL;
            if (tempihsvals->length>10000) sampledtempihsvals = sampleREAL4Vector(tempihsvals, 10000);
            else {
               sampledtempihsvals = XLALCreateREAL4Vector(tempihsvals->length);
               memcpy(sampledtempihsvals->data, tempihsvals->data, sizeof(REAL4)*tempihsvals->length);
            }
            outputfar->ihsdistMean->data[ii-2] = calcMean(sampledtempihsvals);
            outputfar->ihsdistSigma->data[ii-2] = calcStddev(sampledtempihsvals);
            sort_float_ascend(sampledtempihsvals);
            //REAL4 firstval = tempihsvals->data[0];
            //REAL4 lastval = tempihsvals->data[tempihsvals->length-1];
            /* if (ii==2) {
               FILE *IHSSUMVALS = fopen("./output/ihssumvals.dat","w");
               for (jj=0; jj<(INT4)tempihsvals->length; jj++) fprintf(IHSSUMVALS,"%.6f\n",tempihsvals->data[jj]);
               fclose(IHSSUMVALS);
            } */ //To here
            //outputfar->amtToAddToReach95CL->data[ii-2] = tempihsvals->data[(INT4)round(0.95*tempihsvals->length)]-outputfar->ihsdistMean->data[ii-2];
            //if (ihsfar != 1.0) outputfar->ihsfar->data[ii-2] = tempihsvals->data[(INT4)round((1.0-ihsfar)*tempihsvals->length)];
            //else outputfar->ihsfar->data[ii-2] = 0.0;
            XLALDestroyREAL4Vector(tempihsvals);
            XLALDestroyREAL4Vector(sampledtempihsvals);
            
            //REAL8Vector *orderedihsvals = XLALCreateREAL8Vector(101);
            //REAL8Vector *cdftot = XLALCreateREAL8Vector(orderedihsvals->length);
            REAL8 averageval = 0.0, farave = 0.0;
            /* for (jj=0; jj<(INT4)orderedihsvals->length; jj++) {
               orderedihsvals->data[jj] = firstval + (lastval-firstval)*jj*(1.0/orderedihsvals->length);
               cdftot->data[jj] = 0.0;
            } */
            /* REAL4Vector *tempihsvals = XLALCreateREAL4Vector((ihsvectorsequence->length-(ii-1))*ihsvectorsequence->vectorLength);
            if (tempihsvals==NULL) {
               fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, (ihsvectorsequence->length-(ii-1))*ihsvectorsequence->vectorLength);
               XLAL_ERROR_VOID(fn, XLAL_EFUNC);
            } */
            tempihsvals = XLALCreateREAL4Vector((ihsvectorsequence->length-(ii-1))*ihsvectorsequence->vectorLength);
            if (tempihsvals==NULL) {
               fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, (ihsvectorsequence->length-(ii-1))*ihsvectorsequence->vectorLength);
               XLAL_ERROR_VOID(fn, XLAL_EFUNC);
            }
            memcpy(tempihsvals->data, tworows2->data, sizeof(REAL4)*tempihsvals->length);
            if (tempihsvals->length>10000) sampledtempihsvals = sampleREAL4Vector(tempihsvals, 10000);
            else {
               sampledtempihsvals = XLALCreateREAL4Vector(tempihsvals->length);
               memcpy(sampledtempihsvals->data, tempihsvals->data, sizeof(REAL4)*tempihsvals->length);
            }
            for (jj=0; jj<(INT4)sampledtempihsvals->length; jj++) {
               if (sampledtempihsvals->data[jj]!=0.0) {
                  averageval += 1.0;
                  if (params->ihsfar != 1.0) farave += 0.5*gsl_cdf_chisq_Pinv(1.0-params->ihsfar, 2.0*sampledtempihsvals->data[jj]);
                  //for (kk=0; kk<(INT4)orderedihsvals->length; kk++) cdftot->data[kk] += gsl_cdf_chisq_P(2.0*orderedihsvals->data[kk], 2.0*tempihsvals->data[jj]);
               } /* if sampledtempihsvals->data[jj] != 0.0 */
            } /* for jj < sampledtempihsvals->length */
            outputfar->ihsfar->data[ii-2] = farave/averageval;
            /* FILE *CHISQSUMVALS = fopen("./output/chisqsumvals.dat","w");
            for (jj=0; jj<(INT4)orderedihsvals->length; jj++) fprintf(CHISQSUMVALS,"%.6f %.6f\n",orderedihsvals->data[jj],cdftot->data[jj]/averageval);
            fclose(CHISQSUMVALS); */
            //XLALDestroyREAL8Vector(orderedihsvals);
            //XLALDestroyREAL8Vector(cdftot);
            XLALDestroyREAL4Vector(tempihsvals);
            XLALDestroyREAL4Vector(sampledtempihsvals);
            
            REAL4Vector *tempfomvals = XLALCreateREAL4Vector(ihss->length-ii+1);
            if (tempfomvals==NULL) {
               fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ihss->length-ii+1);
               XLAL_ERROR_VOID(fn, XLAL_EFUNC);
            }
            memcpy(tempfomvals->data, &(output->foms->data[0]), sizeof(REAL4)*tempfomvals->length);
            outputfar->ihsfomdistMean->data[ii-2] = calcMean(tempfomvals);
            outputfar->ihsfomdistSigma->data[ii-2] = calcStddev(tempfomvals);
            if (params->ihsfomfar!=1.0 && params->ihsfom==0.0) {
               REAL4Vector *smallestfomvals = XLALCreateREAL4Vector((UINT4)roundf((ihss->length-ii+1)*params->ihsfomfar)+1);
               if (smallestfomvals==NULL) {
                  fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, (INT4)roundf((ihss->length-ii)*params->ihsfomfar)+1);
                  XLAL_ERROR_VOID(fn, XLAL_EFUNC);
               }
               sort_float_smallest(smallestfomvals, tempfomvals);
               if (xlalErrno!=0) {
                  fprintf(stderr, "%s: sort_float_smallest() failed.\n", fn);
                  XLAL_ERROR_VOID(fn, XLAL_EFUNC);
               }
               outputfar->fomfarthresh->data[ii-2] = smallestfomvals->data[smallestfomvals->length-1];
               XLALDestroyREAL4Vector(smallestfomvals);
            } else if (params->ihsfom!=0.0) {
               outputfar->fomfarthresh->data[ii-2] = params->ihsfom;
            } else {
               outputfar->fomfarthresh->data[ii-2] = -1.0;
            }
            XLALDestroyREAL4Vector(tempfomvals);
         } /* if calcPInvVals != 0 */
         
         XLALDestroyREAL4Vector(rowsequencemaxima);
         rowsequencemaxima = NULL;
         XLALDestroyREAL4Vector(fbinmeanvals);
         fbinmeanvals = NULL;
         XLALDestroyINT4Vector(rowsequencelocs);
         rowsequencelocs = NULL;
         
         //fprintf(stderr, "Finished row %d\n", ii);
      } else {
         rowsequencemaxima = XLALCreateREAL4Vector(ii);
         fbinmeanvals = XLALCreateREAL4Vector(ii);
         rowsequencelocs = XLALCreateINT4Vector(ii);
         if (rowsequencemaxima==NULL) {
            fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ii);
            XLAL_ERROR_VOID(fn, XLAL_EFUNC);
         } else if (fbinmeanvals==NULL) {
            fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ii);
            XLAL_ERROR_VOID(fn, XLAL_EFUNC);
         } else if (rowsequencelocs==NULL) {
            fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", fn, ii);
            XLAL_ERROR_VOID(fn, XLAL_EFUNC);
         }
         
         INT4 endloc = ((ii-1)*(ii-1)-(ii-1))/2;
         for (jj=0; jj<(INT4)ihsvectorsequence->length-(ii-1); jj++) {
            for (kk=0; kk<(INT4)ihsvectorsequence->vectorLength; kk++) tworows->data[jj*ihsvectorsequence->vectorLength + kk] += ihsvectorsequence->data[(ii-1+jj)*ihsvectorsequence->vectorLength + kk];
            memcpy(singlerow->data, &(tworows->data[jj*ihsvectorsequence->vectorLength]), sizeof(REAL4)*ihsvectorsequence->vectorLength);
            output->locations->data[(ii-2)*ihss->length-endloc+jj] = max_index(singlerow) + 5;
            output->maxima->data[(ii-2)*ihss->length-endloc+jj] = singlerow->data[output->locations->data[(ii-2)*ihss->length-endloc+jj]-5];
            
            if (calcPInvVals!=0) for (kk=0; kk<(INT4)ihsvalsfromaveNoise->length; kk++) tworows2->data[jj*ihsvalsfromaveNoise->length + kk] += ihsvalsfromaveNoise->data[kk]*randvals->data[ii-1+jj];
            
            memcpy(rowsequencemaxima->data, &(ihss->data[jj]), sizeof(REAL4)*ii);
            memcpy(fbinmeanvals->data, &(FbinMean->data[jj]), sizeof(REAL4)*ii);
            memcpy(rowsequencelocs->data, &(locs->data[jj]), sizeof(INT4)*ii);
            output->foms->data[(ii-2)*ihss->length-endloc+jj] = ihsFOM(rowsequencemaxima, rowsequencelocs, fbinmeanvals, locationnormfactor);
         } /* for jj< ihsvectorsequence->length - (ii-1) */
         
         if (calcPInvVals!=0) {
            REAL4Vector *tempihsvals = XLALCreateREAL4Vector((ihsvectorsequence->length-(ii-1))*ihsvectorsequence->vectorLength);
            if (tempihsvals==NULL) {
               fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, (ihsvectorsequence->length-(ii-1))*ihsvectorsequence->vectorLength);
               XLAL_ERROR_VOID(fn, XLAL_EFUNC);
            }
            memcpy(tempihsvals->data, tworows->data, sizeof(REAL4)*tempihsvals->length);
            REAL4Vector *sampledtempihsvals = NULL;
            if (tempihsvals->length>10000) sampledtempihsvals = sampleREAL4Vector(tempihsvals, 10000);
            else {
               sampledtempihsvals = XLALCreateREAL4Vector(tempihsvals->length);
               memcpy(sampledtempihsvals->data, tempihsvals->data, sizeof(REAL4)*tempihsvals->length);
            }
            outputfar->ihsdistMean->data[ii-2] = calcMean(sampledtempihsvals);
            outputfar->ihsdistSigma->data[ii-2] = calcStddev(sampledtempihsvals);
            sort_float_ascend(sampledtempihsvals);
            //REAL4 firstval = sampledtempihsvals->data[0];
            //REAL4 lastval = sampledtempihsvals->data[sampledtempihsvals->length-1];
            /* if (ii==360) {
               FILE *IHSSUMVALS = fopen("./output/ihssumvals.dat","w");
               for (jj=0; jj<(INT4)sampledtempihsvals->length; jj++) fprintf(IHSSUMVALS,"%.6f\n",sampledtempihsvals->data[jj]);
               fclose(IHSSUMVALS);
            } */
            REAL8 averageval = 0.0, farave = 0.0;
            for (jj=0; jj<(INT4)sampledtempihsvals->length; jj++) {
               if (sampledtempihsvals->data[jj]!=0.0) {
                  averageval += 1.0;
                  if (params->ihsfar != 1.0) farave += 0.5*gsl_cdf_chisq_Pinv(1.0-params->ihsfar, 2.0*sampledtempihsvals->data[jj]);
               }
            }
            outputfar->ihsfar->data[ii-2] = farave/averageval;
            //outputfar->amtToAddToReach95CL->data[ii-2] = sampledtempihsvals->data[(INT4)round(0.95*sampledtempihsvals->length)]-outputfar->ihsdistMean->data[ii-2];
            //if (ihsfar != 1.0) outputfar->ihsfar->data[ii-2] = sampledtempihsvals->data[(INT4)round((1.0-ihsfar)*sampledtempihsvals->length)];
            //else outputfar->ihsfar->data[ii-2] = 0.0;
            XLALDestroyREAL4Vector(tempihsvals);
            XLALDestroyREAL4Vector(sampledtempihsvals);
            
            /* if (ii==360) {
               REAL8Vector *orderedihsvals = XLALCreateREAL8Vector(101);
               REAL8Vector *cdftot = XLALCreateREAL8Vector(orderedihsvals->length);
               REAL8 averageval = 0.0, amtToAddToReach95CL = 0.0, farave = 0.0;
               for (jj=0; jj<(INT4)orderedihsvals->length; jj++) {
                  orderedihsvals->data[jj] = firstval + (lastval-firstval)*jj*(1.0/orderedihsvals->length);
                  cdftot->data[jj] = 0.0;
               }
               tempihsvals = XLALCreateREAL4Vector((ihsvectorsequence->length-(ii-1))*ihsvectorsequence->vectorLength);
               if (tempihsvals==NULL) {
                  fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, (ihsvectorsequence->length-(ii-1))*ihsvectorsequence->vectorLength);
                  XLAL_ERROR_VOID(fn, XLAL_EFUNC);
               }
               memcpy(tempihsvals->data, tworows2->data, sizeof(REAL4)*tempihsvals->length);
               sampledtempihsvals = sampleREAL4Vector(tempihsvals, 10000);
               for (jj=0; jj<(INT4)sampledtempihsvals->length; jj++) {
                  if (sampledtempihsvals->data[jj]!=0.0) {
                     averageval += 1.0;
                     amtToAddToReach95CL += 0.5*gsl_cdf_chisq_Pinv(0.95, 2.0*sampledtempihsvals->data[jj])-sampledtempihsvals->data[jj];
                     if (ihsfar != 1.0) farave += 0.5*gsl_cdf_chisq_Pinv(1.0-ihsfar, 2.0*sampledtempihsvals->data[jj]);
                     for (kk=0; kk<(INT4)orderedihsvals->length; kk++) cdftot->data[kk] += gsl_cdf_chisq_P(2.0*orderedihsvals->data[kk], 2.0*sampledtempihsvals->data[jj]);
                  }
               }
               outputfar->amtToAddToReach95CL->data[ii-2] = amtToAddToReach95CL/averageval;
               outputfar->ihsfar->data[ii-2] = farave/averageval;
               FILE *CHISQSUMVALS = fopen("./output/chisqsumvals.dat","w");
               for (jj=0; jj<(INT4)orderedihsvals->length; jj++) fprintf(CHISQSUMVALS,"%.6f %.6f\n",orderedihsvals->data[jj],cdftot->data[jj]/averageval);
               fclose(CHISQSUMVALS);
               XLALDestroyREAL8Vector(orderedihsvals);
               XLALDestroyREAL8Vector(cdftot);
               XLALDestroyREAL4Vector(tempihsvals);
               XLALDestroyREAL4Vector(sampledtempihsvals);
            } */ /* if ii == 360 */
            
            REAL4Vector *tempfomvals = XLALCreateREAL4Vector(ihss->length-ii+1);
            if (tempfomvals==NULL) {
               fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ihss->length-ii+1);
               XLAL_ERROR_VOID(fn, XLAL_EFUNC);
            }
            memcpy(tempfomvals->data, &(output->foms->data[(ii-2)*ihss->length-endloc]), sizeof(REAL4)*tempfomvals->length);
            outputfar->ihsfomdistMean->data[ii-2] = calcMean(tempfomvals);
            outputfar->ihsfomdistSigma->data[ii-2] = calcStddev(tempfomvals);
            if (params->ihsfomfar!=1.0 && params->ihsfom==0.0) {
               REAL4Vector *smallestfomvals = XLALCreateREAL4Vector((UINT4)roundf((ihss->length-ii+1)*params->ihsfomfar)+1);
               if (smallestfomvals==NULL) {
                  fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, (INT4)roundf((ihss->length-ii)*params->ihsfomfar)+1);
                  XLAL_ERROR_VOID(fn, XLAL_EFUNC);
               }
               sort_float_smallest(smallestfomvals, tempfomvals);
               if (xlalErrno!=0) {
                  fprintf(stderr, "%s: sort_float_smallest() failed.\n", fn);
                  XLAL_ERROR_VOID(fn, XLAL_EFUNC);
               }
               outputfar->fomfarthresh->data[ii-2] = smallestfomvals->data[smallestfomvals->length-1];
               XLALDestroyREAL4Vector(smallestfomvals);
            } else if (params->ihsfom!=0.0) {
               outputfar->fomfarthresh->data[ii-2] = params->ihsfom;
            } else {
               outputfar->fomfarthresh->data[ii-2] = -1.0;
            }
            XLALDestroyREAL4Vector(tempfomvals);
         } /* if calcPInvVals != 0 */
         
         XLALDestroyREAL4Vector(rowsequencemaxima);
         rowsequencemaxima = NULL;
         XLALDestroyREAL4Vector(fbinmeanvals);
         fbinmeanvals = NULL;
         XLALDestroyINT4Vector(rowsequencelocs);
         rowsequencelocs = NULL;
         
         //fprintf(stderr, "Finished row %d\n", ii);
      }
   } /* for ii <= rows */
   XLALDestroyREAL4Vector(singlerow);
   XLALDestroyREAL4VectorSequence(tworows);
   gsl_rng_free(rng);
   
   if (calcPInvVals!=0) {
      XLALDestroyREAL4VectorSequence(tworows2);
      XLALDestroyREAL4Vector(ihsvalsfromaveNoise);
      XLALDestroyREAL4Vector(randvals);
   }
   
} /*ihsSums2_withFAR() */


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



void findIHScandidates(candidateVector *candlist, ihsfarStruct *ihsfarstruct, inputParamsStruct *params, ffdataStruct *ffdata, ihsMaximaStruct *ihsmaxima, REAL4Vector *aveNoise, REAL4Vector *fbinavgs)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii, jj;
   REAL8 fsig, per0, B;
   
   INT4 numberofIHSvalsChecked = 0, numberofIHSvalsExceededThresh = 0, numberPassingBoth = 0;
   
   INT4 numfbins = ffdata->numfbins;
   
   INT4 minrows = (INT4)round(2.0*params->dfmin*params->Tcoh)+1;
   
   REAL4Vector *ihss, *avgsinrange;
   INT4Vector *locs;
   
   //Check the IHS values against the FAR, checking between IHS width values
   for (ii=minrows; ii<=(INT4)ihsfarstruct->ihsfar->length+1; ii++) {
      ihss = XLALCreateREAL4Vector(ii);
      locs = XLALCreateINT4Vector(ii);
      avgsinrange = XLALCreateREAL4Vector(ii);
      if (ihss==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ii);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      } else if (locs==NULL) {
         fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", fn, ii);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      } else if (avgsinrange==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ii);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
      
      REAL8 highestval = 0.0;
      INT4 highestvalloc = -1, jjloc = 0;
      for (jj=0; jj<(INT4)numfbins-(ii-1); jj++) {
      
         //Noise in the range of the rows, mean and rms values for IHS
         memcpy(avgsinrange->data, &(fbinavgs->data[jj]), sizeof(REAL4)*ii);
         REAL4 meanNoise = calcMean(avgsinrange);
         if (XLAL_IS_REAL4_FAIL_NAN(meanNoise)) {
            fprintf(stderr,"%s: calcMean() failed.\n", fn);
            XLAL_ERROR_VOID(fn, XLAL_EFUNC);
         }
         /* REAL4 rmsNoise = calcRms(avgsinrange);
         if (XLAL_IS_REAL4_FAIL_NAN(rmsNoise)) {
            fprintf(stderr,"%s: calcRms() failed.\n", fn);
            XLAL_ERROR_VOID(fn, XLAL_EFUNC);
         } */
         
         numberofIHSvalsChecked++;
         
         INT4 locationinmaximastruct = (ii-2)*numfbins-((ii-1)*(ii-1)-(ii-1))/2+jj;
         
         //Check the IHS sum against the FAR (scaling FAR with mean of the noise in the range of rows)
         if (ihsmaxima->maxima->data[locationinmaximastruct] > ihsfarstruct->ihsfar->data[ii-2]*meanNoise) {
            
            numberofIHSvalsExceededThresh++;
            
            if (ihsfarstruct->fomfarthresh->data[ii-2]==-1.0 || ihsmaxima->foms->data[locationinmaximastruct]<=ihsfarstruct->fomfarthresh->data[ii-2]) {
               
               numberPassingBoth++;
               
               if (ihsmaxima->maxima->data[locationinmaximastruct] > highestval) {
                  highestval = ihsmaxima->maxima->data[locationinmaximastruct];
                  highestvalloc = locationinmaximastruct;
                  jjloc = jj;
               }
               /* INT4 loc = ihsmaxima->locations->data[locationinmaximastruct];
               //Candidate frequency
               fsig = params->fmin + (0.5*(ii-1) + jj)/params->Tcoh;
               //Candidate modulation depth
               B = 0.5*(ii-1)/params->Tcoh;
               //Candidate period
               per0 = params->Tobs/loc;
               //Candidate h0
               REAL8 h0 = ihs2h0_withNoiseSubtraction(ihsmaxima->maxima->data[locationinmaximastruct], loc, jj, ii, params, aveNoise, fbinavgs);
               if (candlist->numofcandidates == candlist->length-1) {
                  candlist = resize_candidateVector(candlist, 2*(candlist->length));
                  if (candlist->data==NULL) {
                     fprintf(stderr,"%s: resize_candidateVector() failed.\n", fn);
                     XLAL_ERROR_VOID(fn, XLAL_EFUNC);
                  }
               }
               loadCandidateData(&candlist->data[candlist->numofcandidates], fsig, per0, B, 0.0, 0.0, ihsmaxima->maxima->data[locationinmaximastruct], h0, 0.0, 0, ffdata->tfnormalization);
               (candlist->numofcandidates)++; */
               
            } /* if fom is below or equal to threshold fom */
         } /* if val exceeds threshold */
      } /* for jj < numfbins-(ii-1) */
      
      if (highestvalloc != -1) {
         INT4 loc = ihsmaxima->locations->data[highestvalloc];
         //Candidate frequency
         fsig = params->fmin + (0.5*(ii-1) + jjloc)/params->Tcoh;
         //Candidate modulation depth
         B = 0.5*(ii-1)/params->Tcoh;
         //Candidate period
         per0 = params->Tobs/loc;
         //Candidate h0
         REAL8 h0 = ihs2h0_withNoiseSubtraction(ihsmaxima->maxima->data[highestvalloc], loc, jjloc, ii, params, aveNoise, fbinavgs);
         
         if (candlist->numofcandidates == candlist->length-1) {
            candlist = resize_candidateVector(candlist, 2*(candlist->length));
            if (candlist->data==NULL) {
               fprintf(stderr,"%s: resize_candidateVector() failed.\n", fn);
               XLAL_ERROR_VOID(fn, XLAL_EFUNC);
            }
         }
         //loadCandidateData(&candlist->data[candlist->numofcandidates], fsig, per0, B, 0.0, 0.0, ihsmaxima->maxima->data[locationinmaximastruct], h0, 0.0, 0, sqrt(ffdata->tfnormalization/2.0*params->Tcoh));
         loadCandidateData(&candlist->data[candlist->numofcandidates], fsig, per0, B, 0.0, 0.0, ihsmaxima->maxima->data[highestvalloc], h0, 0.0, 0, ffdata->tfnormalization);
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
   
   fprintf(stderr,"Number of IHS vals checked = %d, number exceeding IHS threshold = %d, number passing both = %d\n", numberofIHSvalsChecked, numberofIHSvalsExceededThresh, numberPassingBoth);
   
} /* findIHScandidates() */



REAL8 ihs2h0_withNoiseSubtraction(REAL8 ihsval, INT4 location, INT4 lowestfrequencybin, INT4 rows, inputParamsStruct *params, REAL4Vector *aveNoise, REAL4Vector *fbinavgs)
{
   
   const CHAR *fn = __func__;
   
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
      fprintf(stderr, "%s: IHS value is less than expected noise\n", fn);
      return 0.0;
   }
   
   //REAL8 h0 = 1.0*pow((ihsval-totalnoise)/(params->Tcoh*params->Tobs),0.25);
   //REAL8 h0 = 4.6*pow((ihsval-totalnoise)/(params->Tcoh*params->Tobs),0.25);
   //return h0;
   
   REAL8 h0 = ihs2h0(2.0*ihsval-2.0*totalnoise, params, rows);  //With 2.0 for chi-square with 2 d.o.f.
   return h0;
   
}
REAL8 ihs2h0(REAL8 ihsval, inputParamsStruct *params, INT4 rows)
{
   
   //return 4.6*pow(ihsval/(params->Tcoh*params->Tobs),0.25);
   //REAL8 prefact = 4.4*pow(0.5*((REAL8)rows-1.0)/params->Tcoh/7.334e-3,.1);
   REAL8 dfinmHz = 0.5*(rows-1.0)/params->Tcoh*1000.0;
   REAL8 prefact = 1.0;
   prefact = -2.4e-4*dfinmHz*dfinmHz + 4.96e-2*dfinmHz + 3.49;
   return prefact*pow(ihsval/(params->Tcoh*params->Tobs),0.25);
   
}






