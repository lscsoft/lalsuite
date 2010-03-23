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

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <lal/LALMalloc.h>

#include "IHS.h"
#include "candidates.h"

//////////////////////////////////////////////////////////////
// Create vectors for IHS maxima struct  -- done
ihsMaximaStruct * new_ihsMaxima(ffdataStruct *ffdata, INT4 columns)
{

   ihsMaximaStruct *ihsmaxima = (ihsMaximaStruct*)XLALMalloc(sizeof(ihsMaximaStruct));
   
   ihsmaxima->maxima = XLALCreateREAL4Vector(ffdata->f->length*(UINT4)columns);
   ihsmaxima->locations = XLALCreateINT4Vector(ffdata->f->length);
   ihsmaxima->columns = columns;
   
   return ihsmaxima;

}

//////////////////////////////////////////////////////////////
// Destroy vectors for IHS maxima struct  -- done
void free_ihsMaxima(ihsMaximaStruct *data)
{

   XLALDestroyREAL4Vector(data->maxima);
   XLALDestroyINT4Vector(data->locations);
   XLALFree((ihsMaximaStruct*)data);

}


//////////////////////////////////////////////////////////////
// Run the IHS algorithm  -- done
void runIHS(ihsMaximaStruct *out, ffdataStruct *in, INT4 columns)
{

   INT4 ii, jj;
   
   REAL4Vector *column = XLALCreateREAL4Vector(in->fpr->length);
   REAL4Vector *ihss = XLALCreateREAL4Vector(in->f->length);
   ihsVals *ihsvals = new_ihsVals();
   
   //Loop through the columns, 1 frequency at a time
   for (ii=0; ii<(INT4)in->f->length; ii++) {
   
      //For each column, populate it with the data for that frequency bin
      for (jj=0; jj<(INT4)column->length; jj++) column->data[jj] = in->ffdata->data[ii*in->fpr->length + jj];
      
      //Run the IHS algorithm on the column
      incHarmSum(ihsvals, column);
      
      //Temporarily save the IHS value
      ihss->data[ii] = ihsvals->ihs;
      
      //Save the IHS maximum location value for each column
      out->locations->data[ii] = ihsvals->loc;
   }
   
   //Save the maxima for all the column sums
   out->maxima = ihsSums(ihss, columns);
   
   //Save the column widths
   out->columns = columns;
   
   //Destroy variables
   XLALDestroyREAL4Vector(column);
   XLALDestroyREAL4Vector(ihss);
   free_ihsVals(ihsvals);

}


//////////////////////////////////////////////////////////////
// Allocate memory for ihsVals struct  -- done
ihsVals * new_ihsVals(void)
{

   ihsVals *ihsvals = (ihsVals*)XLALMalloc(sizeof(ihsVals));

   return ihsvals;

}

//////////////////////////////////////////////////////////////
// Destroy ihsVals struct  -- done
void free_ihsVals(ihsVals *ihsvals)
{

   XLALFree((ihsVals*)ihsvals);

}


//////////////////////////////////////////////////////////////
// Compute the IHS sum  -- Done
void incHarmSum(ihsVals *out, REAL4Vector *vector)
{
   
   INT4 ii, loc;
   REAL4 ihs;
   
   //Load the stretched spectra
   ihs = 0.0;
   loc = 0;
   for (ii=0; ii<(INT4)vector->length; ii++) {
      //REAL4 sum = vector->data[ii] + 0.5*vector->data[(INT4)floorf(ii*0.5)] + vector->data[(INT4)floorf(ii/3)]/3.0 + 0.25*vector->data[(INT4)floorf(ii*0.25)] + 0.2*vector->data[(INT4)floorf(ii*0.2)];
      REAL4 sum = vector->data[ii] + vector->data[(INT4)floorf(ii*0.5)] + vector->data[(INT4)floorf(ii/3)] + vector->data[(INT4)floorf(ii*0.25)] + vector->data[(INT4)floorf(ii*0.2)];
      if (ii > 5 && sum > ihs) {
         ihs = sum;
         loc = ii;
      }
   }
   
   //Load the outputs into the structure
   out->ihs = ihs;
   out->loc = loc;

}



//////////////////////////////////////////////////////////////
// Allocate memory for ihsfarStruct struct  -- done
ihsfarStruct * new_ihsfarStruct(INT4 columns)
{

   ihsfarStruct *ihsfarstruct = (ihsfarStruct*)XLALMalloc(sizeof(ihsfarStruct));
   
   ihsfarstruct->ihsfar = XLALCreateREAL4Vector((UINT4)columns);
   ihsfarstruct->ihsdistMean = XLALCreateREAL4Vector((UINT4)columns);
   ihsfarstruct->ihsdistSigma = XLALCreateREAL4Vector((UINT4)columns);

   return ihsfarstruct;

}

//////////////////////////////////////////////////////////////
// Destroy ihsfarStruct struct  -- done
void free_ihsfarStruct(ihsfarStruct *ihsfarstruct)
{

   XLALDestroyREAL4Vector(ihsfarstruct->ihsfar);
   XLALDestroyREAL4Vector(ihsfarstruct->ihsdistMean);
   XLALDestroyREAL4Vector(ihsfarstruct->ihsdistSigma);
   XLALFree((ihsfarStruct*)ihsfarstruct);

}


//////////////////////////////////////////////////////////////
// Compute the IHS FAR for a sum of a number of columns  -- done
void genIhsFar(ihsfarStruct *out, ffdataStruct *ffdata, INT4 columns, REAL4 threshold)
{
   
   INT4 ii, jj, kk, length;
   REAL4Vector *noise = NULL;
   
   length = ffdata->fpr->length;
   
   INT4 trials = (INT4)roundf(10000*0.01/threshold);    //Number of trials to determine FAR value
   trials += columns;
   REAL4Vector *ihss = XLALCreateREAL4Vector((UINT4)trials);
   
   //Initialize random number generator
   gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
   srand(time(NULL));
   UINT8 randseed = rand();
   gsl_rng_set(rng, randseed);
   
   ihsVals *ihsvals = new_ihsVals();
   
   //Determine IHS values for the number of trials
   noise = XLALCreateREAL4Vector((UINT4)length);
   for (ii=0; ii<trials; ii++) {
      //Make exponential noise
      for (jj=0; jj<length; jj++) noise->data[jj] = expRandNum(1.0, rng);
      
      //Compute IHS value on exponential noise
      incHarmSum(ihsvals, noise);
      ihss->data[ii] = ihsvals->ihs;
      
   }
   XLALDestroyREAL4Vector(noise);
   
   //Calculate the IHS sum values for the IHS trials
   REAL4Vector *ihssumvals = ihsSums(ihss, columns);
   
   //Now determine distribution values and FAR for the different IHS sum values for each set of columns
   REAL4Vector *tempihsvals = NULL, *topihsvals = NULL;
   INT4 numToRemove = 0;
   for (ii=1; ii<=columns; ii++) {
      
      if (ii>2) numToRemove += ii-1;
      
      //Temporary vector to hold the trial values of IHS column sums
      tempihsvals = XLALCreateREAL4Vector((UINT4)(trials-(ii-1)));
      for (jj=0; jj<(INT4)tempihsvals->length; jj++) {
         tempihsvals->data[jj] = ihssumvals->data[(ii-1)*trials + jj - numToRemove];
         //if (ii==0) tempihsvals->data[jj] = ihssumvals->data[jj];
         //else tempihsvals->data[jj] = ihssumvals->data[ii*trials-(ii-1)+jj];
      }
      
      //Mean and sigma of the various trials
      out->ihsdistMean->data[ii-1] = calcMean(tempihsvals);
      out->ihsdistSigma->data[ii-1] = calcStddev(tempihsvals);
      
      //Launch insertion sort method to find the threshold value
      topihsvals = XLALCreateREAL4Vector((UINT4)roundf((trials-ii)*threshold)+1);
      for (jj=0; jj<(INT4)topihsvals->length; jj++) topihsvals->data[jj] = 0.0;
      topihsvals->data[0] = tempihsvals->data[0];
      for (jj=1; jj<(INT4)topihsvals->length; jj++) {
         INT4 insertionpoint = jj;
         while (insertionpoint > 0 && tempihsvals->data[jj] > topihsvals->data[insertionpoint - 1]) insertionpoint--;
         
         for (kk=topihsvals->length-1; kk>insertionpoint; kk--) topihsvals->data[kk] = topihsvals->data[kk-1];
         topihsvals->data[insertionpoint] = tempihsvals->data[jj];
      }
      for (jj=topihsvals->length; jj<(INT4)tempihsvals->length; jj++) {
         if (tempihsvals->data[jj] > topihsvals->data[topihsvals->length - 1]) {
            INT4 insertionpoint = topihsvals->length - 1;
            while (insertionpoint > 0 && tempihsvals->data[jj] > topihsvals->data[insertionpoint - 1]) insertionpoint--;
            
            for (kk=topihsvals->length-1; kk>insertionpoint; kk--) topihsvals->data[kk] = topihsvals->data[kk-1];
            topihsvals->data[insertionpoint] = tempihsvals->data[jj];
         }
      }
      out->ihsfar->data[ii-1] = topihsvals->data[topihsvals->length-1];
      XLALDestroyREAL4Vector(topihsvals);
      topihsvals = NULL;
      
      //Reset temporary vector
      XLALDestroyREAL4Vector(tempihsvals);
      tempihsvals = NULL;
   }
   
   //Destroy variables
   XLALDestroyREAL4Vector(ihssumvals);
   XLALDestroyREAL4Vector(ihss);
   free_ihsVals(ihsvals);
   gsl_rng_free(rng);
   

}



//////////////////////////////////////////////////////////////
// Compute the IHS sums for a number of columns  -- 
REAL4Vector * ihsSums(REAL4Vector *ihss, INT4 cols)
{
   
   INT4 ii, jj, locInMaximaVector;
   INT4 startPosition = 0;
   
   UINT4 numToRemove = 0;
   for (ii=2; ii<=cols; ii++) numToRemove += (UINT4)(ii-1);
   
   
   //Initialize maxima vector
   REAL4Vector *maxima = XLALCreateREAL4Vector((UINT4)(ihss->length * cols)-numToRemove);
   
   //Start with the vector of single column IHS values
   for (ii=0; ii<(INT4)ihss->length; ii++) maxima->data[ii] = ihss->data[ii];
   
   //Now make the sums. This is more efficient than summing each value separately.
   //We can just use the previous sum and the next value of the single column
   locInMaximaVector = ihss->length;
   for (ii=1; ii<cols; ii++) {
      //startPosition is the start number of the previous width to be summed with the individual IHS value
      startPosition = locInMaximaVector - (INT4)ihss->length + (ii-1); 
      for (jj=0; jj<(INT4)ihss->length-ii; jj++) {
         //maxima->data[ii+jj] is the single column IHS values needed to be added to the total sum
         maxima->data[locInMaximaVector] = maxima->data[startPosition+jj] + maxima->data[ii+jj];
         locInMaximaVector++;
      }
   }
   
   return maxima;

}


//////////////////////////////////////////////////////////////
// Calculate the IHS FOM for a number of columns  -- 
REAL4 ihsFOM(REAL4Vector *ihss, INT4Vector *locs, REAL4Vector *expect, REAL4Vector *sigma)
{

   INT4 ii, maxsnrloc;
   REAL4 maxsnr, fom;
   
   //Create normalized SNR of IHS values
   REAL4Vector *snrs = XLALCreateREAL4Vector(ihss->length);
   for (ii=0; ii<(INT4)snrs->length; ii++) snrs->data[ii] = (ihss->data[ii]-expect->data[ii])/sigma->data[ii];
   
   //Find which pair has the best combined SNR (RMS) and the location
   maxsnr = sqrt(snrs->data[0]*snrs->data[0] + snrs->data[snrs->length-1]*snrs->data[snrs->length-1]);
   maxsnrloc = 0;
   for (ii=1; ii<(INT4)floor(snrs->length*0.5); ii++) {
      if (sqrt(snrs->data[ii]*snrs->data[ii] + snrs->data[snrs->length-ii-1]*snrs->data[snrs->length-ii-1])>maxsnr) {
         maxsnr = sqrt(snrs->data[ii]*snrs->data[ii] + snrs->data[snrs->length-ii-1]*snrs->data[snrs->length-ii-1]);
         maxsnrloc = ii;
      }
   }
   
   //For the highest SNR pair, compute the FOM
   fom = 6.0*(locs->data[maxsnrloc] - locs->data[locs->length-maxsnrloc-1]) * (locs->data[maxsnrloc] - locs->data[locs->length-maxsnrloc-1]);
   
   //Destroy used variables
   XLALDestroyREAL4Vector(snrs);
   
   return fom;

}


//////////////////////////////////////////////////////////////
// Calculate a guess for the location of the brightest pixels
REAL4 ihsLoc(REAL4Vector *ihss, INT4Vector *locs, REAL4Vector *expect, REAL4Vector *sigma)
{

   INT4 ii, maxsnrloc;
   REAL4 maxsnr;
   
   //Create normalized SNR of IHS values
   REAL4Vector *snrs = XLALCreateREAL4Vector(ihss->length);
   for (ii=0; ii<(INT4)snrs->length; ii++) snrs->data[ii] = (ihss->data[ii]-expect->data[ii])/sigma->data[ii];
   
   //Find which pair has the best combined SNR (RMS) and the location
   maxsnr = sqrt(snrs->data[0]*snrs->data[0] + snrs->data[snrs->length-1]*snrs->data[snrs->length-1]);
   maxsnrloc = 0;
   for (ii=1; ii<(INT4)floor(snrs->length*0.5); ii++) {
      if (sqrt(snrs->data[ii]*snrs->data[ii] + snrs->data[snrs->length-ii-1]*snrs->data[snrs->length-ii-1])>maxsnr) {
         maxsnr = sqrt(snrs->data[ii]*snrs->data[ii] + snrs->data[snrs->length-ii-1]*snrs->data[snrs->length-ii-1]);
         maxsnrloc = ii;
      }
   }
   
   //For the highest SNR pair, compute the location
   REAL4 location = 0.5*(locs->data[maxsnrloc] + locs->data[locs->length-maxsnrloc-1]);
   
   //Destroy used variables
   XLALDestroyREAL4Vector(snrs);
   
   return location;

}



void findIHScandidates(candidate *candlist[], INT4 *numofcandidates, ihsfarStruct *ihsfarstruct, REAL4Vector *aveFFnoise, inputParamsStruct *inputParams, ffdataStruct *ffdata, ihsMaximaStruct *ihsmaxima, REAL4 ra, REAL4 dec)
{
   
   INT4 ii, jj, kk, checkbin;
   REAL4 fsig, per0, B;
   REAL4 ihsfomfar = 6.0;
   
   INT4 mincols = (INT4)floorf(2.0*inputParams->dfmin*inputParams->Tcoh)+1;
   INT4 removedbins = 0;
   for (ii=2; ii<=mincols-1; ii++) removedbins += ii-1;
   
   REAL4Vector *ihss, *noiseinrange, *ihsexpect, *ihsstddev;
   INT4Vector *locs;
   checkbin = (INT4)ffdata->f->length*(mincols-1) - removedbins;  //Starting position in the ihsmaxima vector
   //Check the IHS values against the FAR, checking for >1 column widths
   for (ii=mincols-1; ii<(INT4)ihsfarstruct->ihsfar->length; ii++) {
      ihss = XLALCreateREAL4Vector((UINT4)(ii+1));
      locs = XLALCreateINT4Vector((UINT4)(ii+1));
      noiseinrange = XLALCreateREAL4Vector((UINT4)(ii+1));
      ihsexpect = XLALCreateREAL4Vector((UINT4)(ii+1));
      ihsstddev = XLALCreateREAL4Vector(ihsexpect->length);
      for (jj=0; jj<(INT4)ffdata->f->length-ii; jj++) {
      
         //Noise in the range of the columns, mean and rms values for IHS
         for (kk=0; kk<=ii; kk++) {
            noiseinrange->data[kk] = aveFFnoise->data[jj + kk];
            ihsexpect->data[kk] = aveFFnoise->data[jj + kk]*ihsfarstruct->ihsdistMean->data[ii];
            ihsstddev->data[kk] = aveFFnoise->data[jj + kk]*ihsfarstruct->ihsdistSigma->data[ii];
         }
         REAL4 meanNoise = calcMean(noiseinrange);
         REAL4 rmsNoise = calcRms(noiseinrange);
         
         //Check the IHS sum against the FAR (scaling FAR with mean of the noise in the range of columns)
         if (ihsmaxima->maxima->data[checkbin] > ihsfarstruct->ihsfar->data[ii]*meanNoise) {
         
            //Load temporary vectors for determining the FOM
            for (kk=0; kk<=ii; kk++) {
               ihss->data[kk] = ihsmaxima->maxima->data[jj + kk];
               locs->data[kk] = ihsmaxima->locations->data[jj + kk];
            }
            
            //Compute the IHS FOM
            REAL4 fom = ihsFOM(ihss, locs, ihsexpect, ihsstddev);
            
            //Compute the best location
            REAL4 loc = ihsLoc(ihss, locs, ihsexpect, ihsstddev);
         
            //Check the IHS FOM against the FAR, if smaller, and the location is non-zero,
            //and the location is within range then we have a candidate
            if  (fom<=ihsfomfar && loc>=5.0 && inputParams->Tobs/loc>=2.0*3600.0) {
               fprintf(stderr,"IHS statistic = %g (FAR = %g), FOM statistic = %g\n",ihsmaxima->maxima->data[checkbin], ihsfarstruct->ihsfar->data[ii]*meanNoise, fom);
               
               //Candidate frequency
               fsig = inputParams->fmin + (0.5*ii + jj)/inputParams->Tcoh;
               
               //Candidate modulation depth
               B = 0.5*ii/inputParams->Tcoh;
               
               //Candidate period
               per0 = inputParams->Tobs/loc;
               
               REAL4 ihs_sum = ihsmaxima->maxima->data[checkbin];
               REAL4 ihsSnr = (ihs_sum - meanNoise*ihsfarstruct->ihsdistMean->data[ii])/(rmsNoise*ihsfarstruct->ihsdistSigma->data[ii]);
               candlist[(*numofcandidates)] = new_candidate();
               loadCandidateData(candlist[(*numofcandidates)], fsig, per0, B, ra, dec, ihs_sum, ihsSnr, 0.0);
               (*numofcandidates)++;
            }
         }
         
         checkbin++;
      }
      XLALDestroyREAL4Vector(ihss);
      XLALDestroyREAL4Vector(ihsexpect);
      XLALDestroyREAL4Vector(ihsstddev);
      XLALDestroyREAL4Vector(noiseinrange);
      XLALDestroyINT4Vector(locs);
      ihss = NULL;
      locs = NULL;
      ihsexpect = NULL;
      ihsstddev = NULL;
      noiseinrange = NULL;
   }
   
}




