/*
*  Copyright (C) 2010 Evan Goetz
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

#include <gsl/gsl_sort.h>

#include "IHS.h"
#include "TwoSpect.h"
#include "candidates.h"

//////////////////////////////////////////////////////////////
// Create vectors for IHS maxima struct  -- done
ihsMaximaStruct * new_ihsMaxima(INT4 fbins, INT4 columns)
{
   
   const char *fn = __func__;
   
   ihsMaximaStruct *ihsmaxima = XLALMalloc(sizeof(*ihsmaxima));
   if (ihsmaxima==NULL) {
      fprintf(stderr,"%s: XLALMalloc(%lu) failed.\n", fn, sizeof(*ihsmaxima));
      XLAL_ERROR_NULL(fn, XLAL_ENOMEM);
   }
   
   INT4 ii;
   UINT4 numToRemove = 0;
   for (ii=2; ii<=columns; ii++) numToRemove += (UINT4)(ii-1);
   
   ihsmaxima->maxima = XLALCreateREAL4Vector((UINT4)(fbins*columns) - numToRemove);
   ihsmaxima->locations = XLALCreateINT4Vector((UINT4)fbins);
   ihsmaxima->columns = columns;
   
   if (ihsmaxima->maxima==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, fbins*columns-(INT4)numToRemove);
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
   XLALDestroyINT4Vector(data->locations);
   XLALFree((ihsMaximaStruct*)data);

} /* free_ihsMaxima() */


//////////////////////////////////////////////////////////////
// Run the IHS algorithm  -- done
void runIHS(ihsMaximaStruct *output, ffdataStruct *input, inputParamsStruct *params, INT4 columns)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii, jj;
   
   INT4 numffts = (INT4)floor(params->Tobs/(params->Tcoh-params->SFToverlap)-1);
   INT4 numfbins = (INT4)(round(params->fspan*params->Tcoh)+1);
   INT4 numfprbins = (INT4)floorf(numffts*0.5) + 1;
   
   REAL4Vector *column = XLALCreateREAL4Vector((UINT4)numfprbins);
   REAL4Vector *ihss = XLALCreateREAL4Vector((UINT4)numfbins);
   ihsVals *ihsvals = new_ihsVals();
   if (column==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numfprbins);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (ihss==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numfbins);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (ihsvals==NULL) {
      fprintf(stderr,"%s: new_ihsVals() failed.\n", fn);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   
   //Loop through the columns, 1 frequency at a time
   for (ii=0; ii<(INT4)ihss->length; ii++) {
   
      //For each column, populate it with the data for that frequency bin, excluding harmonics of antenna pattern modulation
      for (jj=0; jj<(INT4)column->length; jj++) {
         
         if (fabs(params->Tobs/(24.0*3600.0)-jj)<=1.0 || fabs(params->Tobs/(12.0*3600.0)-jj)<=1.0 || fabs(params->Tobs/(8.0*3600.0)-jj)<=1.0 || fabs(params->Tobs/(6.0*3600.0)-jj)<=1.0) column->data[jj] = 0.0;
         else column->data[jj] = input->ffdata->data[ii*numfprbins + jj];
         
      }
      
      //Run the IHS algorithm on the column
      incHarmSum(ihsvals, column);
      
      //Temporarily save the IHS value
      ihss->data[ii] = ihsvals->ihs;
      
      //Save the IHS maximum location value for each column
      output->locations->data[ii] = ihsvals->loc;
      
   } /* for ii < ihss->length */
   
   //Save the maxima for all the column sums
   ihsSums(output->maxima, ihss, columns);
   /* FILE *IHSDATA = fopen("./ihsdata.dat","w");
   for (ii=0; ii<(INT4)out->maxima->length; ii++) fprintf(IHSDATA,"%g %d\n",out->maxima->data[ii],out->locations->data[ii]);
   fclose(IHSDATA); */
   
   //Save the column widths
   output->columns = columns;
   
   //Destroy variables
   XLALDestroyREAL4Vector(column);
   XLALDestroyREAL4Vector(ihss);
   free_ihsVals(ihsvals);

} /* runIHS() */


//////////////////////////////////////////////////////////////
// Allocate memory for ihsVals struct  -- done
ihsVals * new_ihsVals(void)
{
   
   const char *fn = __func__;
   
   ihsVals *ihsvals = XLALMalloc(sizeof(*ihsvals));
   if (ihsvals==NULL) {
      fprintf(stderr,"%s: XLALMalloc(%lu) failed.\n", fn, sizeof(*ihsvals));
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
         loc = (INT4)(ii/3.0);
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
      fprintf(stderr,"%s: XLALMalloc(%lu) failed.\n", fn, sizeof(*ihsfarstruct));
      XLAL_ERROR_NULL(fn, XLAL_ENOMEM);
   }
   
   ihsfarstruct->ihsfar = XLALCreateREAL4Vector((UINT4)columns);
   ihsfarstruct->ihsdistMean = XLALCreateREAL4Vector((UINT4)columns);
   ihsfarstruct->ihsdistSigma = XLALCreateREAL4Vector((UINT4)columns);
   if (ihsfarstruct->ihsfar==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, columns);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   } else if (ihsfarstruct->ihsdistMean==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, columns);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   } else if( ihsfarstruct->ihsdistSigma==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, columns);
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
   XLALFree((ihsfarStruct*)ihsfarstruct);

} /* free_ihsfarStruct() */


//////////////////////////////////////////////////////////////
// Compute the IHS FAR for a sum of a number of columns  --
void genIhsFar(ihsfarStruct *output, INT4 columns, REAL4 threshold, REAL4Vector *aveNoise, REAL8 Tobs)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii, jj, length;
   REAL4Vector *noise = NULL;
   
   length = (INT4)aveNoise->length;
   
   INT4 trials = (INT4)roundf(10000*0.01/threshold);    //Number of trials to determine FAR value
   trials += columns;
   
   REAL4Vector *ihss = XLALCreateREAL4Vector((UINT4)trials);
   if (ihss==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, trials);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   
   //Initialize random number generator
   gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
   if (rng==NULL) {
      fprintf(stderr,"%s: gsl_rng_alloc() failed.\n", fn);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   //srand(time(NULL));
   //UINT8 randseed = rand();
   //gsl_rng_set(rng, randseed);
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
      
   } /* for ii < trials */
   XLALDestroyREAL4Vector(noise);
   
   //Calculate the IHS sum values for the IHS trials
   ihsMaximaStruct *ihsmaxima = new_ihsMaxima(trials, columns);
   if (ihsmaxima==NULL) {
      fprintf(stderr,"%s: new_ihsMaxima(%d, %d) failed.\n", fn, trials, columns);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   
   //Compute IHS sums
   ihsSums(ihsmaxima->maxima, ihss, columns);
   
   //Now determine distribution values and FAR for the different IHS sum values for each set of columns
   REAL4Vector *tempihsvals = NULL, *topihsvals = NULL;
   INT4 numToRemove = 0;
   for (ii=1; ii<=columns; ii++) {
      
      if (ii>2) numToRemove += ii-2;
      
      //Temporary vector to hold the trial values of IHS column sums
      tempihsvals = XLALCreateREAL4Vector((UINT4)(trials-(ii-1)));
      if (tempihsvals==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, trials-(ii-1));
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
      
      for (jj=0; jj<(INT4)tempihsvals->length; jj++) tempihsvals->data[jj] = ihsmaxima->maxima->data[(ii-1)*trials + jj - numToRemove];
      
      //Mean and sigma of the various trials
      output->ihsdistMean->data[ii-1] = calcMean(tempihsvals);
      output->ihsdistSigma->data[ii-1] = calcStddev(tempihsvals);
      
      //Launch insertion sort method to find the threshold value
      topihsvals = XLALCreateREAL4Vector((UINT4)roundf((trials-ii)*threshold)+1);
      if (topihsvals==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, (INT4)roundf((trials-ii)*threshold)+1);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
      if( (gsl_sort_float_largest((float*)topihsvals->data, topihsvals->length, (float*)tempihsvals->data, 1, tempihsvals->length)) != 0) {
         fprintf(stderr,"%s: gsl_sort_float_largest() failed.\n", fn);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
      
      output->ihsfar->data[ii-1] = topihsvals->data[topihsvals->length-1];
      XLALDestroyREAL4Vector(topihsvals);
      topihsvals = NULL;
      
      //Reset temporary vector
      XLALDestroyREAL4Vector(tempihsvals);
      tempihsvals = NULL;
   } /* for ii <= columns */
   
   /* FILE *IHSFAR = fopen("./ihsfar.dat","w");
   for (ii=0; ii<(INT4)out->ihsfar->length; ii++) fprintf(IHSFAR,"%g\n",out->ihsfar->data[ii]);
   fclose(IHSFAR); */
   
   //Destroy variables
   XLALDestroyREAL4Vector(ihss);
   free_ihsVals(ihsvals);
   free_ihsMaxima(ihsmaxima);
   gsl_rng_free(rng);
   

} /* genIhsFar() */



//////////////////////////////////////////////////////////////
// Compute the IHS sums for a number of columns  -- 
void ihsSums(REAL4Vector *output, REAL4Vector *ihss, INT4 cols)
{
   
   INT4 ii, jj, locInMaximaVector;
   INT4 startPosition = 0;
   
   //Start with the vector of single column IHS values
   for (ii=0; ii<(INT4)ihss->length; ii++) output->data[ii] = ihss->data[ii];
   
   //Now make the sums. This is more efficient than summing each value separately.
   //We can just use the previous sum and the next value of the single column
   locInMaximaVector = ihss->length;
   for (ii=1; ii<cols; ii++) {
      //startPosition is the start number of the previous width to be summed with the individual IHS value
      startPosition = locInMaximaVector - (INT4)ihss->length + (ii-1); 
      for (jj=0; jj<(INT4)ihss->length-ii; jj++) {
         //out->data[ii+jj] is the single column IHS values needed to be added to the total sum
         output->data[locInMaximaVector] = output->data[startPosition+jj] + output->data[ii+jj];
         locInMaximaVector++;
      } /* for jj < ihss->length */
   } /* for ii < cols */
   
} /* ihsSums() */


//////////////////////////////////////////////////////////////
// Calculate the IHS FOM for a number of columns  -- 
REAL4 ihsFOM(REAL4Vector *ihss, INT4Vector *locs, REAL4Vector *sigma)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii, maxsnrloc;
   REAL4 maxsnr, fom;
   
   //Create normalized SNR of IHS values
   REAL4Vector *snrs = XLALCreateREAL4Vector(ihss->length);
   if (snrs==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ihss->length);
      XLAL_ERROR_REAL4(fn, XLAL_EFUNC);
   }
   for (ii=0; ii<(INT4)snrs->length; ii++) snrs->data[ii] = ihss->data[ii]/sigma->data[ii];
   
   //Find which pair has the best combined SNR (RMS) and the location
   maxsnr = sqrtf(snrs->data[0]*snrs->data[0] + snrs->data[snrs->length-1]*snrs->data[snrs->length-1]);
   maxsnrloc = 0;
   //for (ii=1; ii<(INT4)floorf(snrs->length*0.5); ii++) {
   for (ii=1; ii<(INT4)(snrs->length*0.5); ii++) {
      if (sqrtf(snrs->data[ii]*snrs->data[ii] + snrs->data[snrs->length-ii-1]*snrs->data[snrs->length-ii-1])>maxsnr) {
         maxsnr = sqrtf(snrs->data[ii]*snrs->data[ii] + snrs->data[snrs->length-ii-1]*snrs->data[snrs->length-ii-1]);
         maxsnrloc = ii;
      }
   }
   
   //For the highest SNR pair, compute the FOM
   fom = 6.0*(locs->data[maxsnrloc] - locs->data[locs->length-maxsnrloc-1]) * (locs->data[maxsnrloc] - locs->data[locs->length-maxsnrloc-1]);
   
   //Destroy used variables
   XLALDestroyREAL4Vector(snrs);
   
   return fom;

} /* ihsFOM() */


//////////////////////////////////////////////////////////////
// Calculate a guess for the location of the brightest pixels
REAL4 ihsLoc(REAL4Vector *ihss, INT4Vector *locs, REAL4Vector *sigma)
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
   
   //Find which pair has the best combined SNR (RMS) and the location
   maxsnr = sqrtf(snrs->data[0]*snrs->data[0] + snrs->data[snrs->length-1]*snrs->data[snrs->length-1]);
   maxsnrloc = 0;
   //for (ii=1; ii<(INT4)floorf(snrs->length*0.5); ii++) {
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



void findIHScandidates(candidateVector *candlist, ihsfarStruct *ihsfarstruct, inputParamsStruct *params, ffdataStruct *ffdata, ihsMaximaStruct *ihsmaxima, REAL4Vector *fbinavgratios, REAL4Vector *fbinrmsratios)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii, jj, kk, checkbin;
   REAL8 fsig, per0, B;
   REAL4 ihsfomfar = 6.0;
   
   //INT4 numberofIHSvalsChecked = 0;
   
   INT4 numfbins = (INT4)round(params->fspan*params->Tcoh)+1;
   
   /* FILE *IHSVALS = fopen("./realihsvals.dat","w");
   for (ii=0; ii<(INT4)ffdata->f->length; ii++) fprintf(IHSVALS,"%g\n",ihsmaxima->maxima->data[ii]);
   fclose(IHSVALS); */
   
   //Need to shift the start bin location by the number of removed bins in the maxima struct
   INT4 mincols = (INT4)floorf(2.0*params->dfmin*params->Tcoh)+1;
   INT4 removedbins = 0;
   for (ii=2; ii<=mincols-1; ii++) removedbins += ii-1;
   
   REAL4Vector *ihss, *avgsinrange, *rmssinrange, *ihsexpect, *ihsstddev;
   INT4Vector *locs;
   checkbin = numfbins*(mincols-1) - removedbins;  //Starting position in the ihsmaxima vector
   //Check the IHS values against the FAR, checking between IHS width values
   for (ii=mincols-1; ii<(INT4)ihsfarstruct->ihsfar->length; ii++) {
      ihss = XLALCreateREAL4Vector((UINT4)(ii+1));
      locs = XLALCreateINT4Vector((UINT4)(ii+1));
      avgsinrange = XLALCreateREAL4Vector((UINT4)(ii+1));
      rmssinrange = XLALCreateREAL4Vector((UINT4)(ii+1));
      ihsexpect = XLALCreateREAL4Vector((UINT4)(ii+1));
      ihsstddev = XLALCreateREAL4Vector(ihsexpect->length);
      if (ihss==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ii+1);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      } else if (locs==NULL) {
         fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", fn, ii+1);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      } else if (avgsinrange==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ii+1);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      } else if (rmssinrange==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ii+1);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      } else if (ihsexpect==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ii+1);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      } else if (ihsstddev==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, ihsexpect->length);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
      
      for (jj=0; jj<(INT4)numfbins-ii; jj++) {
      
         //Noise in the range of the columns, mean and rms values for IHS
         for (kk=0; kk<=ii; kk++) {
            ihsexpect->data[kk] = ihsfarstruct->ihsdistMean->data[0]*fbinavgratios->data[jj+kk];
            ihsstddev->data[kk] = ihsfarstruct->ihsdistSigma->data[0]*fbinrmsratios->data[jj+kk];
            ihsexpect->data[kk] = ihsfarstruct->ihsdistMean->data[0];
            ihsstddev->data[kk] = ihsfarstruct->ihsdistSigma->data[0];
            avgsinrange->data[kk] = fbinavgratios->data[jj+kk];
            rmssinrange->data[kk] = fbinrmsratios->data[jj+kk];
         } /* for kk <= ii */
         
         REAL4 meanNoise = calcMean(avgsinrange);
         if (XLAL_IS_REAL4_FAIL_NAN(meanNoise)) {
            fprintf(stderr,"%s: calcMean() failed.\n", fn);
            XLAL_ERROR_VOID(fn, XLAL_EFUNC);
         }
         
         //Check the IHS sum against the FAR (scaling FAR with mean of the noise in the range of columns)
         if (ihsmaxima->maxima->data[checkbin] > ihsfarstruct->ihsfar->data[ii]*meanNoise) {
         
            //Load temporary vectors for determining the FOM
            for (kk=0; kk<=ii; kk++) {
               ihss->data[kk] = ihsmaxima->maxima->data[jj + kk];
               locs->data[kk] = ihsmaxima->locations->data[jj + kk];
            }
            
            //Compute the IHS FOM
            REAL4 fom = ihsFOM(ihss, locs, ihsstddev);
            if (XLAL_IS_REAL4_FAIL_NAN(fom)) {
               fprintf(stderr,"%s: ihsFOM() failed.\n", fn);
               XLAL_ERROR_VOID(fn, XLAL_EFUNC);
            }
            
            //Compute the best location
            REAL4 loc = ihsLoc(ihss, locs, ihsstddev);
            if (XLAL_IS_REAL4_FAIL_NAN(loc)) {
               fprintf(stderr,"%s: ihsLoc() failed.\n", fn);
               XLAL_ERROR_VOID(fn, XLAL_EFUNC);
            }
         
            //Check the IHS FOM against the FAR, if smaller, and the location is non-zero,
            //and the location is within range then we have a candidate
            if  (fom<=ihsfomfar && loc>=5.0 && params->Tobs/loc>=2.0*3600.0) {
               
               //Candidate frequency
               fsig = params->fmin + (0.5*ii + jj)/params->Tcoh;
               
               //Candidate modulation depth
               B = 0.5*ii/params->Tcoh;
               
               //Candidate period
               per0 = params->Tobs/loc;
               
               //fprintf(stderr,"IHS candidate %d: f0 = %g, P = %g, df = %g\n",(*numofcandidates),fsig,per0,B);
               
               if (candlist->numofcandidates == candlist->length-1) {
                  candlist = resize_candidateVector(candlist, 2*(candlist->length));
                  if (candlist->numofcandidates == candlist->length-1) XLAL_ERROR_VOID(fn, XLAL_EFUNC);
               }
               loadCandidateData(&candlist->data[candlist->numofcandidates], fsig, per0, B, 0.0, 0.0, 0.0, 0.0, 0.0, 0, sqrt(ffdata->tfnormalization/2.0*params->Tcoh));
               (candlist->numofcandidates)++;
            } /* if fom is below threshold and within period limits */
         } /* if val exceeds threshold */
         
         checkbin++;
      } /* for jj < numfbins-ii */
      
      //Destroy
      XLALDestroyREAL4Vector(ihss);
      XLALDestroyREAL4Vector(ihsexpect);
      XLALDestroyREAL4Vector(ihsstddev);
      XLALDestroyREAL4Vector(avgsinrange);
      XLALDestroyREAL4Vector(rmssinrange);
      XLALDestroyINT4Vector(locs);
      ihss = NULL;
      locs = NULL;
      ihsexpect = NULL;
      ihsstddev = NULL;
      avgsinrange = NULL;
      rmssinrange = NULL;
      
   } /* for ii < ihsfarstruct->ihsfar->length */
   
} /* findIHScandidates() */




