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

#include "IHS.h"
#include "TwoSpect.h"
#include "candidates.h"

//////////////////////////////////////////////////////////////
// Create vectors for IHS maxima struct  -- done
ihsMaximaStruct * new_ihsMaxima(INT4 fbins, INT4 columns)
{

   ihsMaximaStruct *ihsmaxima = (ihsMaximaStruct*)XLALMalloc(sizeof(ihsMaximaStruct));
   
   INT4 ii;
   UINT4 numToRemove = 0;
   for (ii=2; ii<=columns; ii++) numToRemove += (UINT4)(ii-1);
   
   ihsmaxima->maxima = XLALCreateREAL8Vector((UINT4)(fbins*columns) - numToRemove);
   ihsmaxima->locations = XLALCreateINT4Vector((UINT4)fbins);
   ihsmaxima->columns = columns;
   
   return ihsmaxima;

}

//////////////////////////////////////////////////////////////
// Destroy vectors for IHS maxima struct  -- done
void free_ihsMaxima(ihsMaximaStruct *data)
{

   XLALDestroyREAL8Vector(data->maxima);
   XLALDestroyINT4Vector(data->locations);
   XLALFree((ihsMaximaStruct*)data);

}


//////////////////////////////////////////////////////////////
// Run the IHS algorithm  -- done
void runIHS(ihsMaximaStruct *out, ffdataStruct *in, REAL8 Tobs, INT4 columns)
{

   INT4 ii, jj;
   
   REAL8Vector *column = XLALCreateREAL8Vector(in->fpr->length);
   REAL8Vector *ihss = XLALCreateREAL8Vector(in->f->length);
   ihsVals *ihsvals = new_ihsVals();
   
   //Loop through the columns, 1 frequency at a time
   for (ii=0; ii<(INT4)in->f->length; ii++) {
   
      //For each column, populate it with the data for that frequency bin (TRY dividing the expected noise)
      for (jj=0; jj<(INT4)column->length; jj++) {
         
         if (fabs(Tobs/(24.0*3600.0)-jj)<=1.0 || fabs(Tobs/(12.0*3600.0)-jj)<=1.0 || fabs(Tobs/(8.0*3600.0)-jj)<=1.0 || fabs(Tobs/(6.0*3600.0)-jj)<=1.0) column->data[jj] = 0.0;
         else column->data[jj] = in->ffdata->data[ii*in->fpr->length + jj];
         //else column->data[jj] = in->ffdata->data[ii*in->fpr->length + jj]/(tfnoiseratio->data[ii]*ffnoise->data[jj]);
         //column->data[jj] = in->ffdata->data[ii*in->fpr->length + jj]/(tfnoiseratio->data[ii]*ffnoise->data[jj]);
         
      }
      
      //Run the IHS algorithm on the column
      incHarmSum(ihsvals, column);
      
      //Temporarily save the IHS value
      ihss->data[ii] = ihsvals->ihs;
      
      //Save the IHS maximum location value for each column
      out->locations->data[ii] = ihsvals->loc;
   }
   
   //Save the maxima for all the column sums
   //out->maxima = ihsSums(ihss, columns);
   ihsSums(out->maxima, ihss, columns);
   /* FILE *IHSDATA = fopen("./ihsdata.dat","w");
   for (ii=0; ii<(INT4)out->maxima->length; ii++) fprintf(IHSDATA,"%g %d\n",out->maxima->data[ii],out->locations->data[ii]);
   fclose(IHSDATA); */
   
   //Save the column widths
   out->columns = columns;
   
   //Destroy variables
   XLALDestroyREAL8Vector(column);
   XLALDestroyREAL8Vector(ihss);
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
void incHarmSum(ihsVals *out, REAL8Vector *in)
{
   
   INT4 ii, loc;
   REAL8 ihs;
   
   ihs = 0.0;
   loc = 0;
   //Start ii >= 15
   for (ii=15; ii<(INT4)in->length; ii++) {
      //REAL4 sum = in->data[ii] + 0.5*in->data[(INT4)floorf(ii*0.5)] + in->data[(INT4)floorf(ii/3)]/3.0 + 0.25*in->data[(INT4)floorf(ii*0.25)] + 0.2*in->data[(INT4)floorf(ii*0.2)];
      REAL8 sum = in->data[ii] + in->data[(INT4)floor(ii*0.5)] + in->data[(INT4)floor(ii/3.0)] + in->data[(INT4)floor(ii*0.25)] + in->data[(INT4)floor(ii*0.2)];
      if (sum > ihs) {
         ihs = sum;
         loc = (INT4)floor(ii/3.0);
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
   
   ihsfarstruct->ihsfar = XLALCreateREAL8Vector((UINT4)columns);
   ihsfarstruct->ihsdistMean = XLALCreateREAL8Vector((UINT4)columns);
   ihsfarstruct->ihsdistSigma = XLALCreateREAL8Vector((UINT4)columns);

   return ihsfarstruct;

}

//////////////////////////////////////////////////////////////
// Destroy ihsfarStruct struct  -- done
void free_ihsfarStruct(ihsfarStruct *ihsfarstruct)
{

   XLALDestroyREAL8Vector(ihsfarstruct->ihsfar);
   XLALDestroyREAL8Vector(ihsfarstruct->ihsdistMean);
   XLALDestroyREAL8Vector(ihsfarstruct->ihsdistSigma);
   XLALFree((ihsfarStruct*)ihsfarstruct);

}


//////////////////////////////////////////////////////////////
// Compute the IHS FAR for a sum of a number of columns  --
void genIhsFar(ihsfarStruct *out, INT4 columns, REAL8 threshold, REAL8Vector *aveNoise, REAL8 Tobs)
{
   
   INT4 ii, jj, kk, length;
   REAL8Vector *noise = NULL;
   
   length = aveNoise->length;
   
   INT4 trials = (INT4)roundf(10000*0.01/threshold);    //Number of trials to determine FAR value
   trials += columns;
   
   REAL8Vector *ihss = XLALCreateREAL8Vector((UINT4)trials);
   
   //Initialize random number generator
   gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
   srand(time(NULL));
   UINT8 randseed = rand();
   gsl_rng_set(rng, randseed);
   
   ihsVals *ihsvals = new_ihsVals();
   
   //Determine IHS values for the number of trials
   noise = XLALCreateREAL8Vector((UINT4)length);
   for (ii=0; ii<trials; ii++) {
      //Make exponential noise (TRY dividing the expected noise)
      //for (jj=0; jj<(INT4)aveNoise->length; jj++) noise->data[jj] = expRandNum(aveNoise->data[jj], rng);
      for (jj=0; jj<(INT4)aveNoise->length; jj++) {
         if (fabs(Tobs/(24.0*3600.0)-jj)<=1.0 || fabs(Tobs/(12.0*3600.0)-jj)<=1.0 || fabs(Tobs/(8.0*3600.0)-jj)<=1.0 || fabs(Tobs/(6.0*3600.0)-jj)<=1.0) noise->data[jj] = 0.0;
         else noise->data[jj] = expRandNum(aveNoise->data[jj], rng);
      }
      
      //Compute IHS value on exponential noise
      incHarmSum(ihsvals, noise);
      
      ihss->data[ii] = ihsvals->ihs;
      
   }
   XLALDestroyREAL8Vector(noise);
   
   //Calculate the IHS sum values for the IHS trials
   ihsMaximaStruct *ihsmaxima = new_ihsMaxima(trials, columns);
   ihsSums(ihsmaxima->maxima, ihss, columns);
   
   //Now determine distribution values and FAR for the different IHS sum values for each set of columns
   REAL8Vector *tempihsvals = NULL, *topihsvals = NULL;
   INT4 numToRemove = 0;
   for (ii=1; ii<=columns; ii++) {
      
      if (ii>2) numToRemove += ii-2;
      
      //Temporary vector to hold the trial values of IHS column sums
      tempihsvals = XLALCreateREAL8Vector((UINT4)(trials-(ii-1)));
      for (jj=0; jj<(INT4)tempihsvals->length; jj++) tempihsvals->data[jj] = ihsmaxima->maxima->data[(ii-1)*trials + jj - numToRemove];
      
      //Mean and sigma of the various trials
      out->ihsdistMean->data[ii-1] = calcMean(tempihsvals);
      out->ihsdistSigma->data[ii-1] = calcStddev(tempihsvals);
      
      //Launch insertion sort method to find the threshold value
      topihsvals = XLALCreateREAL8Vector((UINT4)roundf((trials-ii)*threshold)+1);
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
      XLALDestroyREAL8Vector(topihsvals);
      topihsvals = NULL;
      
      //Reset temporary vector
      XLALDestroyREAL8Vector(tempihsvals);
      tempihsvals = NULL;
   }
   
   /* FILE *IHSFAR = fopen("./ihsfar.dat","w");
   for (ii=0; ii<(INT4)out->ihsfar->length; ii++) fprintf(IHSFAR,"%g\n",out->ihsfar->data[ii]);
   fclose(IHSFAR); */
   
   //Destroy variables
   XLALDestroyREAL8Vector(ihss);
   free_ihsVals(ihsvals);
   free_ihsMaxima(ihsmaxima);
   gsl_rng_free(rng);
   

}



//////////////////////////////////////////////////////////////
// Compute the IHS sums for a number of columns  -- 
void ihsSums(REAL8Vector *out, REAL8Vector *ihss, INT4 cols)
{
   
   INT4 ii, jj, locInMaximaVector;
   INT4 startPosition = 0;
   
   //Start with the vector of single column IHS values
   for (ii=0; ii<(INT4)ihss->length; ii++) out->data[ii] = ihss->data[ii];
   
   //Now make the sums. This is more efficient than summing each value separately.
   //We can just use the previous sum and the next value of the single column
   locInMaximaVector = ihss->length;
   for (ii=1; ii<cols; ii++) {
      //startPosition is the start number of the previous width to be summed with the individual IHS value
      startPosition = locInMaximaVector - (INT4)ihss->length + (ii-1); 
      for (jj=0; jj<(INT4)ihss->length-ii; jj++) {
         //out->data[ii+jj] is the single column IHS values needed to be added to the total sum
         out->data[locInMaximaVector] = out->data[startPosition+jj] + out->data[ii+jj];
         locInMaximaVector++;
      }
   }
   
}


//////////////////////////////////////////////////////////////
// Calculate the IHS FOM for a number of columns  -- 
REAL8 ihsFOM(REAL8Vector *ihss, INT4Vector *locs, REAL8Vector *sigma)
{

   INT4 ii, maxsnrloc;
   REAL8 maxsnr, fom;
   
   //Create normalized SNR of IHS values
   REAL8Vector *snrs = XLALCreateREAL8Vector(ihss->length);
   for (ii=0; ii<(INT4)snrs->length; ii++) snrs->data[ii] = ihss->data[ii]/sigma->data[ii];
   
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
   XLALDestroyREAL8Vector(snrs);
   
   return fom;

}


//////////////////////////////////////////////////////////////
// Calculate a guess for the location of the brightest pixels
REAL8 ihsLoc(REAL8Vector *ihss, INT4Vector *locs, REAL8Vector *sigma)
{

   INT4 ii, maxsnrloc;
   REAL8 maxsnr;
   
   //Create normalized SNR of IHS values
   REAL8Vector *snrs = XLALCreateREAL8Vector(ihss->length);
   for (ii=0; ii<(INT4)snrs->length; ii++) snrs->data[ii] = ihss->data[ii]/sigma->data[ii];
   
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
   REAL8 location = 0.5*(locs->data[maxsnrloc] + locs->data[locs->length-maxsnrloc-1]);
   
   //Destroy used variables
   XLALDestroyREAL8Vector(snrs);
   
   return location;

}



void findIHScandidates(candidate *candlist[], INT4 *numofcandidates, ihsfarStruct *ihsfarstruct, inputParamsStruct *inputParams, ffdataStruct *ffdata, ihsMaximaStruct *ihsmaxima, REAL8Vector *fbinavgratios, REAL8Vector *fbinrmsratios)
{
   
   INT4 ii, jj, kk, checkbin;
   REAL8 fsig, per0, B;
   REAL8 ihsfomfar = 6.0;
   
   /* FILE *IHSVALS = fopen("./realihsvals.dat","w");
   for (ii=0; ii<(INT4)ffdata->f->length; ii++) fprintf(IHSVALS,"%g\n",ihsmaxima->maxima->data[ii]);
   fclose(IHSVALS); */
   
   //Need to shift the start bin location by the number of removed bins in the maxima struct
   INT4 mincols = (INT4)floorf(2.0*inputParams->dfmin*inputParams->Tcoh)+1;
   INT4 removedbins = 0;
   for (ii=2; ii<=mincols-1; ii++) removedbins += ii-1;
   
   REAL8Vector *ihss, *avgsinrange, *rmssinrange, *ihsexpect, *ihsstddev;
   INT4Vector *locs;
   checkbin = (INT4)ffdata->f->length*(mincols-1) - removedbins;  //Starting position in the ihsmaxima vector
   //Check the IHS values against the FAR, checking between IHS width values
   for (ii=mincols-1; ii<(INT4)ihsfarstruct->ihsfar->length; ii++) {
      ihss = XLALCreateREAL8Vector((UINT4)(ii+1));
      locs = XLALCreateINT4Vector((UINT4)(ii+1));
      avgsinrange = XLALCreateREAL8Vector((UINT4)(ii+1));
      rmssinrange = XLALCreateREAL8Vector((UINT4)(ii+1));
      ihsexpect = XLALCreateREAL8Vector((UINT4)(ii+1));
      ihsstddev = XLALCreateREAL8Vector(ihsexpect->length);
      for (jj=0; jj<(INT4)ffdata->f->length-ii; jj++) {
      
         //Noise in the range of the columns, mean and rms values for IHS
         for (kk=0; kk<=ii; kk++) {
            ihsexpect->data[kk] = ihsfarstruct->ihsdistMean->data[0]*fbinavgratios->data[jj+kk];
            ihsstddev->data[kk] = ihsfarstruct->ihsdistSigma->data[0]*fbinrmsratios->data[jj+kk];
            ihsexpect->data[kk] = ihsfarstruct->ihsdistMean->data[0];
            ihsstddev->data[kk] = ihsfarstruct->ihsdistSigma->data[0];
            avgsinrange->data[kk] = fbinavgratios->data[jj+kk];
            rmssinrange->data[kk] = fbinrmsratios->data[jj+kk];
         }
         
         REAL8 meanNoise = calcMean(avgsinrange);
         //REAL8 rmsNoise = calcRms(rmssinrange);
         
         //Check the IHS sum against the FAR (scaling FAR with mean of the noise in the range of columns)
         if (ihsmaxima->maxima->data[checkbin] > ihsfarstruct->ihsfar->data[ii]*meanNoise) {
         //if (ihsmaxima->maxima->data[checkbin] > ihsfarstruct->ihsfar->data[ii]) {
         
            //Load temporary vectors for determining the FOM
            for (kk=0; kk<=ii; kk++) {
               ihss->data[kk] = ihsmaxima->maxima->data[jj + kk];
               locs->data[kk] = ihsmaxima->locations->data[jj + kk];
            }
            
            //Compute the IHS FOM
            REAL8 fom = ihsFOM(ihss, locs, ihsstddev);
            
            //Compute the best location
            REAL8 loc = ihsLoc(ihss, locs, ihsstddev);
         
            //Check the IHS FOM against the FAR, if smaller, and the location is non-zero,
            //and the location is within range then we have a candidate
            if  (fom<=ihsfomfar && loc>=5.0 && inputParams->Tobs/loc>=2.0*3600.0) {
               
               //Candidate frequency
               fsig = inputParams->fmin + (0.5*ii + jj)/inputParams->Tcoh;
               
               //Candidate modulation depth
               B = 0.5*ii/inputParams->Tcoh;
               
               //Candidate period
               per0 = inputParams->Tobs/loc;
               
               //fprintf(stderr,"IHS candidate %d: f0 = %g, P = %g, df = %g\n",(*numofcandidates),fsig,per0,B);
               
               //REAL4 ihs_sum = ihsmaxima->maxima->data[checkbin];
               //REAL4 ihsSnr = (ihs_sum - meanNoise*ihsfarstruct->ihsdistMean->data[ii])/(rmsNoise*ihsfarstruct->ihsdistSigma->data[ii]);
               candlist[(*numofcandidates)] = new_candidate();
               loadCandidateData(candlist[(*numofcandidates)], fsig, per0, B, 0.0, 0.0, 0.0, 0.0, 0.0, 0, sqrt(ffdata->ffnormalization)*sqrt(ffdata->tfnormalization));
               (*numofcandidates)++;
            }
         }
         
         checkbin++;
      }
      XLALDestroyREAL8Vector(ihss);
      XLALDestroyREAL8Vector(ihsexpect);
      XLALDestroyREAL8Vector(ihsstddev);
      XLALDestroyREAL8Vector(avgsinrange);
      XLALDestroyREAL8Vector(rmssinrange);
      XLALDestroyINT4Vector(locs);
      ihss = NULL;
      locs = NULL;
      ihsexpect = NULL;
      ihsstddev = NULL;
      avgsinrange = NULL;
      rmssinrange = NULL;
   }
   
}




