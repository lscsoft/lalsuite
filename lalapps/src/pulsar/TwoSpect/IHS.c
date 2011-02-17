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

#include <lal/LALMalloc.h>

#include "IHS.h"
#include "statistics.h"
#include "candidates.h"

//////////////////////////////////////////////////////////////
// Create vectors for IHS maxima struct  -- done
ihsMaximaStruct * new_ihsMaxima(INT4 fbins, INT4 columns)
{
   
   const char *fn = __func__;
   
   ihsMaximaStruct *ihsmaxima = XLALMalloc(sizeof(*ihsmaxima));
   if (ihsmaxima==NULL) {
      fprintf(stderr,"%s: XLALMalloc(%zu) failed.\n", fn, sizeof(*ihsmaxima));
      XLAL_ERROR_NULL(fn, XLAL_ENOMEM);
   }
   
   INT4 ii;
   UINT4 numToRemove = 0;
   for (ii=1; ii<=columns; ii++) numToRemove += (UINT4)(ii-1);
   
   ihsmaxima->maxima = XLALCreateREAL4Vector((UINT4)(fbins*columns) - numToRemove);
   ihsmaxima->foms = XLALCreateREAL4Vector(ihsmaxima->maxima->length);
   ihsmaxima->ihsForEachFbin = XLALCreateREAL4Vector((UINT4)fbins);
   ihsmaxima->locations = XLALCreateINT4Vector((UINT4)fbins);
   ihsmaxima->columns = columns;
   
   if (ihsmaxima->maxima==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, fbins*columns-(INT4)numToRemove);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   } else if (ihsmaxima->foms==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, fbins*columns-(INT4)numToRemove);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   } else if (ihsmaxima->ihsForEachFbin==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, fbins);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   } else if (ihsmaxima->locations==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", fn, fbins);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   }
      
   return ihsmaxima;

} /* new_ihsMaxima() */


//////////////////////////////////////////////////////////////
// Destroy vectors for IHS maxima struct  -- done
void free_ihsMaxima(ihsMaximaStruct *data)
{

   XLALDestroyREAL4Vector(data->maxima);
   XLALDestroyREAL4Vector(data->foms);
   XLALDestroyREAL4Vector(data->ihsForEachFbin);
   XLALDestroyINT4Vector(data->locations);
   XLALFree((ihsMaximaStruct*)data);

} /* free_ihsMaxima() */


//////////////////////////////////////////////////////////////
// Run the IHS algorithm  -- done
void runIHS(ihsMaximaStruct *output, ffdataStruct *input, inputParamsStruct *params, INT4 columns, REAL4Vector *FbinMean)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii, jj;
   
   INT4 numffts = (INT4)floor(params->Tobs/(params->Tcoh-params->SFToverlap)-1);
   INT4 numfbins = (INT4)(round(params->fspan*params->Tcoh)+1);
   INT4 numfprbins = (INT4)floorf(numffts*0.5) + 1;
   
   REAL4Vector *column = XLALCreateREAL4Vector((UINT4)numfprbins);
   REAL4Vector *ihss = XLALCreateREAL4Vector((UINT4)numfbins);
   INT4Vector *locs = XLALCreateINT4Vector((UINT4)numfbins);
   ihsVals *ihsvals = new_ihsVals();
   if (column==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numfprbins);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (ihss==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numfbins);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (locs==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", fn, numfbins);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (ihsvals==NULL) {
      fprintf(stderr,"%s: new_ihsVals() failed.\n", fn);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   
   //Loop through the columns, 1 frequency at a time
   for (ii=0; ii<(INT4)ihss->length; ii++) {
   
      //For each column, populate it with the data for that frequency bin, excluding harmonics of antenna pattern modulation
      memcpy(column->data, &(input->ffdata->data[ii*numfprbins]), sizeof(REAL4)*column->length);
      for (jj=0; jj<(INT4)column->length; jj++) {
         if (fabs(params->Tobs/(24.0*3600.0)-jj)<=1.0 || fabs(params->Tobs/(12.0*3600.0)-jj)<=1.0 || fabs(params->Tobs/(8.0*3600.0)-jj)<=1.0 || fabs(params->Tobs/(6.0*3600.0)-jj)<=1.0) column->data[jj] = 0.0;
      }
      
      //Run the IHS algorithm on the column
      incHarmSum(ihsvals, column);
      
      //Save the maximum IHS value
      ihss->data[ii] = ihsvals->ihs;
      
      //Save the IHS maximum location value for each column
      locs->data[ii] = ihsvals->loc;
      
   } /* for ii < ihss->length */
   
   //Save the maxima for all the column sums
   ihsSums(output, ihss, locs, columns, FbinMean);
   /* FILE *IHSDATA = fopen("./ihsdata.dat","w");
   for (ii=0; ii<(INT4)out->maxima->length; ii++) fprintf(IHSDATA,"%g %d\n",out->maxima->data[ii],out->locations->data[ii]);
   fclose(IHSDATA); */
   
   //Save the column widths
   //output->columns = columns;
   
   //Destroy variables
   XLALDestroyREAL4Vector(column);
   XLALDestroyREAL4Vector(ihss);
   XLALDestroyINT4Vector(locs);
   free_ihsVals(ihsvals);

} /* runIHS() */


//////////////////////////////////////////////////////////////
// Allocate memory for ihsVals struct  -- done
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
// Destroy ihsVals struct  -- done
void free_ihsVals(ihsVals *ihsvals)
{

   XLALFree((ihsVals*)ihsvals);

} /* free_ihsVals() */


//////////////////////////////////////////////////////////////
// Compute the IHS sum  -- Done
void incHarmSum(ihsVals *output, REAL4Vector *input)
{
   
   INT4 ii, loc;
   REAL4 ihs;
   
   ihs = 0.0;
   loc = 0;
   //Start ii >= 15
   for (ii=15; ii<(INT4)input->length; ii++) {
      REAL4 sum = input->data[ii] + input->data[(INT4)(ii*0.5)] + input->data[(INT4)(ii/3.0)] + input->data[(INT4)(ii*0.25)] + input->data[(INT4)(ii*0.2)];

      if (sum > ihs) {
         ihs = sum;
         loc = (INT4)round(ii/3.0);
      }
   } /* for ii < input->length */
   
   //Load the outputs into the structure
   output->ihs = ihs;
   output->loc = loc;

} /* incHarmSum() */


//////////////////////////////////////////////////////////////
// Allocate memory for ihsfarStruct struct  -- done
ihsfarStruct * new_ihsfarStruct(INT4 columns)
{
   
   const char *fn = __func__;
   
   ihsfarStruct *ihsfarstruct = XLALMalloc(sizeof(*ihsfarstruct));
   if (ihsfarstruct == NULL) {
      fprintf(stderr,"%s: XLALMalloc(%zu) failed.\n", fn, sizeof(*ihsfarstruct));
      XLAL_ERROR_NULL(fn, XLAL_ENOMEM);
   }
   
   ihsfarstruct->ihsfar = XLALCreateREAL4Vector((UINT4)columns-1);
   ihsfarstruct->ihsdistMean = XLALCreateREAL4Vector((UINT4)columns-1);
   ihsfarstruct->ihsdistSigma = XLALCreateREAL4Vector((UINT4)columns-1);
   ihsfarstruct->fomfarthresh = XLALCreateREAL4Vector((UINT4)columns-1);
   if (ihsfarstruct->ihsfar==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, columns-1);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   } else if (ihsfarstruct->ihsdistMean==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, columns-1);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   } else if( ihsfarstruct->ihsdistSigma==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, columns-1);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   } else if( ihsfarstruct->fomfarthresh==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, columns-1);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   }
   
   return ihsfarstruct;

} /* new_ihsfarStruct() */


//////////////////////////////////////////////////////////////
// Destroy ihsfarStruct struct  -- done
void free_ihsfarStruct(ihsfarStruct *ihsfarstruct)
{

   XLALDestroyREAL4Vector(ihsfarstruct->ihsfar);
   XLALDestroyREAL4Vector(ihsfarstruct->ihsdistMean);
   XLALDestroyREAL4Vector(ihsfarstruct->ihsdistSigma);
   XLALDestroyREAL4Vector(ihsfarstruct->fomfarthresh);
   XLALFree((ihsfarStruct*)ihsfarstruct);

} /* free_ihsfarStruct() */


//////////////////////////////////////////////////////////////
// Compute the IHS FAR for a sum of a number of columns  --
void genIhsFar(ihsfarStruct *output, inputParamsStruct *params, INT4 columns, REAL4Vector *aveNoise)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii, jj, length;
   REAL4Vector *noise = NULL;
   REAL8 Tobs = params->Tobs;
   
   length = (INT4)aveNoise->length;
   
   INT4 trials = (INT4)round(10000*0.01/params->ihsfar);    //Number of trials to determine FAR value
   if (params->ihsfomfar!=0.0 && trials<(INT4)round(10000*0.01/params->ihsfomfar)) {
      trials = (INT4)round(10000*0.01/params->ihsfomfar);
   }
   trials += columns;
   
   REAL4Vector *ihss = XLALCreateREAL4Vector((UINT4)trials);
   INT4Vector *locs = XLALCreateINT4Vector((UINT4)trials);
   if (ihss==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, trials);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (locs==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", fn, trials);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   
   //Initialize random number generator
   gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
   if (rng==NULL) {
      fprintf(stderr,"%s: gsl_rng_alloc() failed.\n", fn);
      XLAL_ERROR_VOID(fn, XLAL_ENOMEM);
   }
   /* srand(time(NULL));
   UINT8 randseed = rand();
   gsl_rng_set(rng, randseed); */
   gsl_rng_set(rng, 0);
   
   ihsVals *ihsvals = new_ihsVals();
   if (ihsvals==NULL) {
      fprintf(stderr,"%s: new_ihsVals() failed.\n", fn);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   
   //Determine IHS values for the number of trials
   noise = XLALCreateREAL4Vector((UINT4)length);
   if (noise==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, length);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   for (ii=0; ii<trials; ii++) {
      //Make exponential noise removing harmonics of 24 hours to match with the same method as real analysis
      for (jj=0; jj<length; jj++) {
         if (fabs(Tobs/(24.0*3600.0)-jj)<=1.0 || fabs(Tobs/(12.0*3600.0)-jj)<=1.0 || fabs(Tobs/(8.0*3600.0)-jj)<=1.0 || fabs(Tobs/(6.0*3600.0)-jj)<=1.0) noise->data[jj] = 0.0;
         else noise->data[jj] = (REAL4)expRandNum(aveNoise->data[jj], rng);
      } /* for jj < length */
      
      //Compute IHS value on exponential noise
      incHarmSum(ihsvals, noise);
      
      ihss->data[ii] = ihsvals->ihs;
      locs->data[ii] = ihsvals->loc;
      
   } /* for ii < trials */
   XLALDestroyREAL4Vector(noise);
   
   //Calculate the IHS sum values for the IHS trials
   ihsMaximaStruct *ihsmaxima = new_ihsMaxima(trials, columns);
   if (ihsmaxima==NULL) {
      fprintf(stderr,"%s: new_ihsMaxima(%d, %d) failed.\n", fn, trials, columns);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   
   //Create a fake vector with the same average value in each bin
   REAL4Vector *FbinMean = XLALCreateREAL4Vector(trials);
   if (FbinMean==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", fn, trials);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   for (ii=0; ii<trials; ii++) FbinMean->data[ii] = 1.0;
   
   //Compute IHS sums and FOMs
   ihsSums(ihsmaxima, ihss, locs, columns, FbinMean);
   
   //Now determine distribution values and FAR for the different IHS sum values for each set of columns
   //Also determine distribution of FOM values
   REAL4Vector *tempihsvals = NULL, *topihsvals = NULL;
   REAL4Vector *tempfomvals = NULL, *smallestfomvals = NULL;
   INT4 numToRemove = 0;
   for (ii=2; ii<=columns; ii++) {
      
      if (ii>2) numToRemove += ii-2;
      
      //Temporary vector to hold the trial values of IHS column sums
      tempihsvals = XLALCreateREAL4Vector((UINT4)(trials-(ii-1)));
      if (tempihsvals==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, trials-(ii-1));
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
      //for (jj=0; jj<(INT4)tempihsvals->length; jj++) tempihsvals->data[jj] = ihsmaxima->maxima->data[(ii-1)*trials + jj - numToRemove];
      memcpy(tempihsvals->data, &(ihsmaxima->maxima->data[(ii-2)*trials - numToRemove]), sizeof(REAL4)*tempihsvals->length);
      //Mean and sigma of the various trials
      output->ihsdistMean->data[ii-2] = calcMean(tempihsvals);
      output->ihsdistSigma->data[ii-2] = calcStddev(tempihsvals);
      //Launch qsort method to find the threshold value if the threshold value is not equal to 1.0
      //If equal to 1.0, set the FAR value equal to 0
      if (params->ihsfar!=1.0) {
         topihsvals = XLALCreateREAL4Vector((UINT4)roundf((trials-ii+1)*params->ihsfar)+1);
         if (topihsvals==NULL) {
            fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, (INT4)roundf((trials-ii)*params->ihsfar)+1);
            XLAL_ERROR_VOID(fn, XLAL_EFUNC);
         }
         sort_float_largest(topihsvals, tempihsvals);
         if (xlalErrno!=0) {
            fprintf(stderr, "%s: sort_float_largest() failed.\n", fn);
            XLAL_ERROR_VOID(fn, XLAL_EFUNC);
         }
         output->ihsfar->data[ii-2] = topihsvals->data[topihsvals->length-1];
         XLALDestroyREAL4Vector(topihsvals);
         topihsvals = NULL;
      } else {
         output->ihsfar->data[ii-2] = 0.0;
      }
      //Reset temporary vector
      XLALDestroyREAL4Vector(tempihsvals);
      tempihsvals = NULL;
      
      if (params->ihsfom!=0.0) {
         output->fomfarthresh->data[ii-2] = params->ihsfom;
      } else {
         //Temporary vector to hold the trial values of IHS FOMs
         tempfomvals = XLALCreateREAL4Vector((UINT4)(trials-(ii-1)));
         if (tempfomvals==NULL) {
            fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, trials-(ii-1));
            XLAL_ERROR_VOID(fn, XLAL_EFUNC);
         }
         memcpy(tempfomvals->data, &(ihsmaxima->foms->data[(ii-2)*trials - numToRemove]), sizeof(REAL4)*tempfomvals->length);
         //Launch qsort method to find the threshold value if the threshold value is not equal to 1.0
         //If equal to 1.0, set the FAR value equal to -1
         if (params->ihsfomfar!=1.0) {
            smallestfomvals = XLALCreateREAL4Vector((UINT4)roundf((trials-ii+1)*params->ihsfomfar)+1);
            if (smallestfomvals==NULL) {
               fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, (INT4)roundf((trials-ii)*params->ihsfomfar)+1);
               XLAL_ERROR_VOID(fn, XLAL_EFUNC);
            }
            sort_float_smallest(smallestfomvals, tempfomvals);
            if (xlalErrno!=0) {
               fprintf(stderr, "%s: sort_float_smallest() failed.\n", fn);
               XLAL_ERROR_VOID(fn, XLAL_EFUNC);
            }
            output->fomfarthresh->data[ii-2] = smallestfomvals->data[smallestfomvals->length-1];
            XLALDestroyREAL4Vector(smallestfomvals);
            smallestfomvals = NULL;
         } else {
            output->fomfarthresh->data[ii-2] = -1.0;
         }
         //Reset temporary vector
         XLALDestroyREAL4Vector(tempfomvals);
         tempfomvals = NULL;
      }
      
   } /* for ii <= columns */
   
   //Destroy variables
   XLALDestroyREAL4Vector(ihss);
   free_ihsVals(ihsvals);
   free_ihsMaxima(ihsmaxima);
   gsl_rng_free(rng);
   

} /* genIhsFar() */



//////////////////////////////////////////////////////////////
// Compute the IHS sums for a number of columns  -- 
void ihsSums(ihsMaximaStruct *output, REAL4Vector *ihss, INT4Vector *locs, INT4 cols, REAL4Vector *FbinMean)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii, jj, locInMaximaVector;
   INT4 startPosition = 0;
   INT4 slotsmoved = 0;
   
   //Start with the vector of single column IHS values
   memcpy(output->ihsForEachFbin->data, ihss->data, sizeof(REAL4)*ihss->length);
   memcpy(output->locations->data, locs->data, sizeof(INT4)*locs->length);
   
   output->columns = cols;
   
   //Efficiently sum the IHS values
   //Sum the pairs
   for (ii=0; ii<(INT4)ihss->length-1; ii++) output->maxima->data[ii] = ihss->data[ii] + ihss->data[ii+1];
   locInMaximaVector = ii; //save the position in the maxima vector
   //Now continue summing only needing to add the next column IHS value to the recent summed value
   for (ii=3; ii<=cols; ii++) {
      slotsmoved = 0; //Reset the number of slots moved
      for (jj=ii-1; jj<(INT4)ihss->length; jj++) {
         output->maxima->data[locInMaximaVector] = output->maxima->data[startPosition+jj-(ii-1)] + ihss->data[jj];
         locInMaximaVector++; //Increment position in maxima vector
         slotsmoved++; //Also increment the number of slots moved
      }
      startPosition = locInMaximaVector - slotsmoved; //Position to start is the location in the maxima vector minus slots moved
   } /* for ii<= cols */
   
   INT4Vector *templocs = NULL;
   REAL4Vector *tempihss = NULL, *tempfbinmean = NULL;
   locInMaximaVector = 0;
   for (ii=2; ii<=cols; ii++) {
      templocs = XLALCreateINT4Vector(ii);
      tempihss = XLALCreateREAL4Vector(ii);
      tempfbinmean = XLALCreateREAL4Vector(ii);
      if (templocs==NULL) {
         fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", fn, ii);
         XLAL_ERROR_VOID(fn, XLAL_ENOMEM);
      } else if (tempihss==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ii);
         XLAL_ERROR_VOID(fn, XLAL_ENOMEM);
      } else if (tempfbinmean==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ii);
         XLAL_ERROR_VOID(fn, XLAL_ENOMEM);
      }
      for (jj=0; jj<(INT4)locs->length-(ii-1); jj++) {
         memcpy(templocs->data, &(locs->data[jj]), sizeof(INT4)*templocs->length);
         memcpy(tempihss->data, &(ihss->data[jj]), sizeof(REAL4)*tempihss->length);
         memcpy(tempfbinmean->data, &(FbinMean->data[jj]), sizeof(REAL4)*tempfbinmean->length);
         output->foms->data[locInMaximaVector] = ihsFOM(tempihss, templocs, tempfbinmean);
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


//////////////////////////////////////////////////////////////
// Calculate the IHS FOM for a number of columns  -- 
REAL4 ihsFOM(REAL4Vector *ihss, INT4Vector *locs, REAL4Vector *sigma)
//REAL4 ihsFOM(REAL4Vector *ihss, INT4Vector *locs, REAL4 sigma)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii, maxsnrloc;
   REAL4 maxsnr, fom;
   
   //Create normalized SNR of IHS values
   REAL4Vector *snrs = XLALCreateREAL4Vector(ihss->length);
   if (snrs==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ihss->length);
      XLAL_ERROR_REAL4(fn, XLAL_ENOMEM);
   }
   for (ii=0; ii<(INT4)snrs->length; ii++) snrs->data[ii] = ihss->data[ii]/sigma->data[ii];
   //for (ii=0; ii<(INT4)snrs->length; ii++) snrs->data[ii] = ihss->data[ii]/sigma;
   
   //Find which pair has the best combined SNR (RMS) and the location
   maxsnr = sqrtf(snrs->data[0]*snrs->data[0] + snrs->data[snrs->length-1]*snrs->data[snrs->length-1]);
   maxsnrloc = 0;
   for (ii=1; ii<(INT4)(snrs->length*0.5); ii++) {
      if (sqrtf(snrs->data[ii]*snrs->data[ii] + snrs->data[snrs->length-ii-1]*snrs->data[snrs->length-ii-1])>maxsnr) {
         maxsnr = sqrtf(snrs->data[ii]*snrs->data[ii] + snrs->data[snrs->length-ii-1]*snrs->data[snrs->length-ii-1]);
         maxsnrloc = ii;
      }
   }
   
   //For the highest SNR pair, compute the FOM
   fom = 12.0*(locs->data[maxsnrloc] - locs->data[locs->length-maxsnrloc-1]) * (locs->data[maxsnrloc] - locs->data[locs->length-maxsnrloc-1]);
   
   //Destroy used variables
   XLALDestroyREAL4Vector(snrs);
   
   return fom;

} /* ihsFOM() */


//////////////////////////////////////////////////////////////
// Calculate a guess for the location of the brightest pixels
REAL4 ihsLoc(REAL4Vector *ihss, INT4Vector *locs, REAL4Vector *sigma)
//REAL4 ihsLoc(REAL4Vector *ihss, INT4Vector *locs, REAL4 sigma)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii, maxsnrloc;
   REAL4 maxsnr;
   
   //Create normalized SNR of IHS values
   REAL4Vector *snrs = XLALCreateREAL4Vector(ihss->length);
   if (snrs==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ihss->length);
      XLAL_ERROR_REAL4(fn, XLAL_EFUNC);
   }
   for (ii=0; ii<(INT4)snrs->length; ii++) snrs->data[ii] = ihss->data[ii]/sigma->data[ii];
   //for (ii=0; ii<(INT4)snrs->length; ii++) snrs->data[ii] = ihss->data[ii]/sigma;
   
   //Find which pair has the best combined SNR (RMS) and the location
   maxsnr = sqrtf(snrs->data[0]*snrs->data[0] + snrs->data[snrs->length-1]*snrs->data[snrs->length-1]);
   maxsnrloc = 0;
   for (ii=1; ii<(INT4)(snrs->length*0.5); ii++) {
      if (sqrtf(snrs->data[ii]*snrs->data[ii] + snrs->data[snrs->length-ii-1]*snrs->data[snrs->length-ii-1])>maxsnr) {
         maxsnr = sqrtf(snrs->data[ii]*snrs->data[ii] + snrs->data[snrs->length-ii-1]*snrs->data[snrs->length-ii-1]);
         maxsnrloc = ii;
      }
   }
   
   //For the highest SNR pair, compute the location
   REAL4 location = 0.5*(locs->data[maxsnrloc] + locs->data[locs->length-maxsnrloc-1]);
   
   //Destroy used variables
   XLALDestroyREAL4Vector(snrs);
   
   return location;

} /* ihsLoc() */



void findIHScandidates(candidateVector *candlist, ihsfarStruct *ihsfarstruct, inputParamsStruct *params, ffdataStruct *ffdata, ihsMaximaStruct *ihsmaxima, REAL4Vector *fbinavgs)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii, jj, kk, checkbin;
   REAL8 fsig, per0, B;
   
   //INT4 numberofIHSvalsChecked = 0, numberofIHSvalsExceededThresh = 0;
   
   INT4 numfbins = (INT4)round(params->fspan*params->Tcoh)+1;
   
   //Need to shift the start bin location by the number of removed bins in the maxima struct
   INT4 mincols = (INT4)(2.0*params->dfmin*params->Tcoh)+1;
   INT4 removedbins = 0;
   for (ii=2; ii<=mincols-1; ii++) removedbins += ii-1;
   
   REAL4Vector *ihss, *avgsinrange, *ihsexpect;//, *ihsstddev;
   //REAL4Vector *ihss, *avgsinrange;
   INT4Vector *locs;
   checkbin = numfbins*(mincols-2) - removedbins;  //Starting position in the ihsmaxima vector
   //Check the IHS values against the FAR, checking between IHS width values
   for (ii=mincols; ii<=(INT4)ihsfarstruct->ihsfar->length+1; ii++) {
      ihss = XLALCreateREAL4Vector(ii);
      locs = XLALCreateINT4Vector(ii);
      ihsexpect = XLALCreateREAL4Vector(ii);
      avgsinrange = XLALCreateREAL4Vector(ii);
      if (ihss==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ii);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      } else if (locs==NULL) {
         fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", fn, ii);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      } else if (ihsexpect==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ii);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      } else if (avgsinrange==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ii);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
      
      for (jj=0; jj<(INT4)numfbins-(ii-1); jj++) {
      
         //Noise in the range of the columns, mean and rms values for IHS
         for (kk=0; kk<ii; kk++) {
            ihsexpect->data[kk] = 0.5*ihsfarstruct->ihsdistMean->data[0]*fbinavgs->data[jj+kk];
            //ihsexpect->data[kk] = ihsfarstruct->ihsdistMean->data[0];
         } /* for kk <= ii */
         memcpy(avgsinrange->data, &(fbinavgs->data[jj]), sizeof(REAL4)*ii);
         
         REAL4 meanNoise = calcMean(avgsinrange);
         if (XLAL_IS_REAL4_FAIL_NAN(meanNoise)) {
            fprintf(stderr,"%s: calcMean() failed.\n", fn);
            XLAL_ERROR_VOID(fn, XLAL_EFUNC);
         }
         
         //numberofIHSvalsChecked++;
         
         //Check the IHS sum against the FAR (scaling FAR with mean of the noise in the range of columns)
         if (ihsmaxima->maxima->data[checkbin] > ihsfarstruct->ihsfar->data[ii-2]*meanNoise) {
            
            //numberofIHSvalsExceededThresh++;
            
            if (ihsfarstruct->fomfarthresh->data[ii-2]==-1.0 || ihsmaxima->foms->data[checkbin]<=ihsfarstruct->fomfarthresh->data[ii-2]) {
               for (kk=0; kk<ii; kk++) {
                  ihss->data[kk] = ihsmaxima->ihsForEachFbin->data[jj + kk];
                  locs->data[kk] = ihsmaxima->locations->data[jj + kk];
               }
               REAL4 loc = ihsLoc(ihss, locs, ihsexpect);
               if (XLAL_IS_REAL4_FAIL_NAN(loc)) {
                  fprintf(stderr,"%s: ihsLoc() failed.\n", fn);
                  XLAL_ERROR_VOID(fn, XLAL_EFUNC);
               }               
               if (loc>=5.0 && params->Tobs/loc>=2.0*3600.0) {
                  //Candidate frequency
                  fsig = params->fmin + (0.5*ii + jj)/params->Tcoh;
                  //Candidate modulation depth
                  B = 0.5*ii/params->Tcoh;
                  //Candidate period
                  per0 = params->Tobs/loc;
                  if (candlist->numofcandidates == candlist->length-1) {
                     candlist = resize_candidateVector(candlist, 2*(candlist->length));
                     if (candlist->data==NULL) {
                        fprintf(stderr,"%s: resize_candidateVector() failed.\n", fn);
                        XLAL_ERROR_VOID(fn, XLAL_EFUNC);
                     }
                  }
                  loadCandidateData(&candlist->data[candlist->numofcandidates], fsig, per0, B, 0.0, 0.0, ihsmaxima->maxima->data[checkbin], ihsmaxima->foms->data[checkbin], 0.0, 0, sqrt(ffdata->tfnormalization/2.0*params->Tcoh));
                  (candlist->numofcandidates)++;
               }
            } /* if fom is below or equal to threshold fom */
         } /* if val exceeds threshold */
         
         checkbin++;
      } /* for jj < numfbins-(ii-1) */
      
      //Destroy
      XLALDestroyREAL4Vector(ihss);
      XLALDestroyREAL4Vector(avgsinrange);
      XLALDestroyREAL4Vector(ihsexpect);
      XLALDestroyINT4Vector(locs);
      ihss = NULL;
      locs = NULL;
      avgsinrange = NULL;
      ihsexpect = NULL;
      
   } /* for ii < ihsfarstruct->ihsfar->length */
   
   //fprintf(stderr,"Number of IHS vals checked = %d, number exceeding threshold = %d\n", numberofIHSvalsChecked, numberofIHSvalsExceededThresh);
   
} /* findIHScandidates() */




