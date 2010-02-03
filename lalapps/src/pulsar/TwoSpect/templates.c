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
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <lal/LALConstants.h>
#include <lal/LALMalloc.h>
#include <lal/Window.h>
#include "templates.h"

//////////////////////////////////////////////////////////////
// Allocate memory for farStruct struct  -- done
farStruct * new_farStruct(void)
{

   farStruct *farstruct = (farStruct*)XLALMalloc(sizeof(farStruct));

   return farstruct;

}

//////////////////////////////////////////////////////////////
// Destroy farStruct struct  -- done
void free_farStruct(farStruct *farstruct)
{

   XLALFree((farStruct*)farstruct);

}


//////////////////////////////////////////////////////////////
// Estimate the FAR of the R statistic from the weights
void estimateFAR(farStruct *out, REAL4Vector *weights, topbinsStruct *topbinsstruct, REAL4 thresh, REAL4Vector *ffplanenoise)
{
   
   INT4 ii, jj;
   
   REAL4 sumofsqweights = 0;
   for (ii=0; ii<topbinsstruct->topbins->length; ii++) {
      sumofsqweights += (weights->data[ topbinsstruct->topbins->data[ii] ]*
         weights->data[ topbinsstruct->topbins->data[ii] ]);
   }
   REAL4 sumofsqweightsinv = 1.0/sumofsqweights;
   
   INT4 trials = 10000;    //Number of trials to determine FAR value
   REAL4Vector *Rs = XLALCreateREAL4Vector((UINT4)trials);
   
   //RandomParams *param = XLALCreateRandomParams(0);
   gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
   
   for (ii=0; ii<trials; ii++) {
      //Create noise value and R value
      REAL4 R = 0.0;
      for (jj=0; jj<topbinsstruct->topbins->length; jj++) {
         REAL4 noise = expRandNum(ffplanenoise->data[ topbinsstruct->freqoftopbins->data[jj] ], rng);
         R += ((noise-ffplanenoise->data[ topbinsstruct->freqoftopbins->data[jj] ])*
            weights->data[ topbinsstruct->topbins->data[jj] ]);
      }
      Rs->data[ii] = R*sumofsqweightsinv;
   }
   REAL4 mean = calcMean(Rs);
   REAL4 sigma = calcStddev(Rs);
   
   //Find the maximums and get to the FAR
   REAL4 max;
   INT4 maxloc;
   for (ii=0; ii<(INT4)roundf(thresh*trials)+1; ii++) {
      max = 0.0;
      maxloc = 0;
      for (jj=0; jj<trials; jj++) {
         if (Rs->data[jj] > max) {
            max = Rs->data[jj];
            maxloc = jj;
         }
      }
      Rs->data[maxloc] = 0.0;
   }
   
   //Destroy
   XLALDestroyREAL4Vector(Rs);
   gsl_rng_free(rng);
   
   
   out->far = max;
   out->distMean = mean;
   out->distSigma = sigma;

}




//////////////////////////////////////////////////////////////
// Make an estimated template based on FFT of train of Gaussians  -- done
void makeTemplateGaussians(ffdataStruct *out, candidate *in)
{

   INT4 ii, jj, numfbins, numffts, N;
   
   numfbins = (INT4)(roundf(in->fspan*in->Tcoh)+1);   //Number of frequency bins
   numffts = (INT4)floor(2*(in->Tobs/in->Tcoh)-1);     //Number of FFTs
   N = (INT4)floor(in->Tobs/in->period);     //Number of Gaussians
   
   REAL4 periodf = 1.0/in->period;
   
   //Set up frequencies and determine separation in time of peaks for each frequency
   REAL4Vector *phi_actual = XLALCreateREAL4Vector((UINT4)numfbins);
   for (ii=0; ii<phi_actual->length; ii++) {
      out->f->data[ii] = in->fmin + ii/in->Tcoh;
      if ( fabs(out->f->data[ii]-in->fsig)/in->moddepth <= 1 ) {
         phi_actual->data[ii] = 0.5*in->period - asin(fabs(out->f->data[ii]-in->fsig)/
            in->moddepth)*LAL_1_PI*in->period;
      } else {
         phi_actual->data[ii] = 0.0;
      }
   }
   
   //Create second FFT frequencies
   for (ii=0; ii<out->fpr->length; ii++) out->fpr->data[ii] = ii/in->Tobs;
   //Need to add a small number to first frequency or else get a divide by zero problem
   //out->fpr->data[0] = 1e-9;
   
   //Scale used for "spillover" into bins outside of phi_actual
   REAL4 k = in->moddepth*in->Tcoh;    //amplitude of modulation in units of bins
   REAL4Vector *scale = XLALCreateREAL4Vector(out->f->length);      //the scaling factor
   INT4 m0 = (INT4)roundf(in->fsig*in->Tcoh) - (INT4)roundf(in->fmin*in->Tcoh) + 1;   //central frequency bin
   //printf("Round sig bin = %d, rounded min bin = %d\n",(INT4)floor(in->fsig*in->Tcoh+0.5),
   //   (INT4)floor(in->fmin*in->Tcoh+0.5));
   //printf("Warning: adding 1 bin to move center bin of template!\n");
   //TODO: figure out why an extra bin needs to be added to move the signal
   INT4 mextent = (INT4)fabs(floor(in->moddepth*in->Tcoh));   //Bins filled by modulation
   REAL4 overage = (k-mextent)-1;
   for (ii=0; ii<scale->length; ii++) {
      if (mextent != 0) {
         if (ii < m0-mextent-2 || ii > m0+mextent+2) {
            scale->data[ii] = 0;
         } else if (ii == m0-mextent-2 || ii == m0+mextent+2) {
            scale->data[ii] = sincxoverxsqminusone(overage-1)*sincxoverxsqminusone(overage-1);
         } else if (ii == m0-mextent-1 || ii == m0+mextent+1) {
            scale->data[ii] = sincxoverxsqminusone(overage)*sincxoverxsqminusone(overage);
         } else {
            scale->data[ii] = 1;
         }
      } else {
         if (ii < m0-mextent-3 || ii > m0+mextent+2) {
            scale->data[ii] = 0;
         } else if (ii == m0-mextent-3 || ii == m0+mextent+2) {
            scale->data[ii] = sincxoverxsqminusone(overage-1)*sincxoverxsqminusone(overage-1);
         } else if (ii == m0-mextent-2 || ii == m0+mextent+1) {
            scale->data[ii] = sincxoverxsqminusone(overage)*sincxoverxsqminusone(overage);
         } else {
            scale->data[ii] = 1;
         }
      }
   }
   
   //Make sigmas for each frequency
   REAL4Vector *sigmas = XLALCreateREAL4Vector((UINT4)numfbins);
   REAL4Vector *wvals = XLALCreateREAL4Vector((UINT4)floor(2*in->period/in->Tcoh));
   REAL4Vector *allsigmas = XLALCreateREAL4Vector((UINT4)(wvals->length * numfbins));
   for (ii=0; ii<wvals->length; ii++) {         //t = (ii+1)*in->Tcoh*0.5
      REAL4 sigbin = (in->moddepth*cos(LAL_TWOPI*periodf*((ii+1)*in->Tcoh*0.5))+in->fsig)*in->Tcoh;
      REAL4 sigbinvelocity = fabs(-in->moddepth*sin(LAL_TWOPI*periodf*((ii+1)*in->Tcoh*0.5))*
         in->Tcoh*0.5*in->Tcoh*LAL_TWOPI*periodf);
      REAL4 sigma = 0.5*in->Tcoh*((383.85*LAL_1_PI)*(0.5*6.1e-3)/((sigbinvelocity+0.1769)*
         (sigbinvelocity+0.1769)+(0.5*6.1e-3)*(0.5*6.1e-3))+0.3736);   //Derived fit from simulation
      for (jj=0; jj<numfbins; jj++) {
         allsigmas->data[ii*numfbins + jj] = sincxoverxsqminusone(sigbin-
            (INT4)roundf(in->fmin*in->Tcoh+jj))*sincxoverxsqminusone(sigbin-
            (INT4)roundf(in->fmin*in->Tcoh+jj))*sigma;
      }
   }
   for (ii=0; ii<numfbins; ii++) {
      for (jj=0; jj<wvals->length; jj++) {
         wvals->data[jj] = allsigmas->data[ii + jj*numfbins]*allsigmas->data[ii + jj*numfbins];
      }
      sigmas->data[ii] = sqrt(calcMean(wvals));
   }
   
   //Create template
   REAL4 sum = 0;
   for (ii=0; ii<out->f->length; ii++) {
      REAL4 s = sigmas->data[ii];
      REAL4 scale1 = 1.0/(1.0+exp(-phi_actual->data[ii]*phi_actual->data[ii]*0.5/(s*s)));
      for (jj=0; jj<out->fpr->length; jj++) {
      
         if (scale->data[ii] == 0) {
            out->ffdata->data[ii*out->fpr->length + jj] = 0;
         } else {
            if (jj==0) {
               out->ffdata->data[ii*out->fpr->length + jj] = scale->data[ii]*scale1*
                  4.0*LAL_TWOPI*s*s*N*N;
            } else if (fabs(cos(in->period*LAL_TWOPI*out->fpr->data[jj])-1.0)<1e-6) {
               out->ffdata->data[ii*out->fpr->length + jj] = scale->data[ii]*scale1*
                  2.0*LAL_TWOPI*s*s*exp(-s*s*LAL_TWOPI*LAL_TWOPI*out->fpr->data[jj]*
                  out->fpr->data[jj])*(cos(phi_actual->data[ii]*LAL_TWOPI*out->fpr->data[jj])+1.0)*N*N;
            } else {
               out->ffdata->data[ii*out->fpr->length + jj] = scale->data[ii]*scale1*
                  2.0*LAL_TWOPI*s*s*exp(-s*s*LAL_TWOPI*LAL_TWOPI*out->fpr->data[jj]*
                  out->fpr->data[jj])*(cos(N*in->period*LAL_TWOPI*out->fpr->data[jj])-1.0)*
                  (cos(phi_actual->data[ii]*LAL_TWOPI*out->fpr->data[jj])+1.0)/
                  (cos(in->period*LAL_TWOPI*out->fpr->data[jj])-1.0);
            }
         
            /* out->ffdata->data[ii*out->fpr->length + jj] = scale->data[ii]*scale1*LAL_TWOPI*s*s;
            
            gsl_complex val1 = gsl_complex_polar(exp(-s*s*LAL_TWOPI*LAL_PI*out->fpr->data[jj]*
               out->fpr->data[jj]),-phi_actual->data[ii]*LAL_TWOPI*out->fpr->data[jj]);
            
            gsl_complex val2 = gsl_complex_polar(1,phi_actual->data[ii]*LAL_TWOPI*out->fpr->data[jj]);
            val2 = gsl_complex_add_real(val2,1);
            
            gsl_complex val3 = gsl_complex_polar(1,-in->period*N*LAL_TWOPI*out->fpr->data[jj]);
            val3 = gsl_complex_sub_real(val3,1);
            
            gsl_complex val4 = gsl_complex_polar(1,in->period*LAL_TWOPI*out->fpr->data[jj]);
            val4 = gsl_complex_sub_real(val4,1);
            
            gsl_complex val = gsl_complex_mul(val1,val2);
            val = gsl_complex_mul(val,val3);
            val = gsl_complex_div(val,val4);
            
            REAL4 absval2 = gsl_complex_abs2(val);
            out->ffdata->data[ii*out->fpr->length + jj] *= absval2; */
            
            if (out->ffdata->data[ii*out->fpr->length + jj] <= 1e-6) out->ffdata->data[ii*out->fpr->length + jj] = 0.0;
         }
         
         sum += out->ffdata->data[ii*out->fpr->length + jj];
      }
   }
   
   //Normalize
   for (ii=0; ii<out->f->length; ii++) {
      for (jj=0; jj<out->fpr->length; jj++) out->ffdata->data[ii*out->fpr->length + jj] /= sum;
   }
   
   
   //Destroy variables
   XLALDestroyREAL4Vector(phi_actual);
   XLALDestroyREAL4Vector(scale);
   XLALDestroyREAL4Vector(sigmas);
   XLALDestroyREAL4Vector(allsigmas);
   XLALDestroyREAL4Vector(wvals);

}


//////////////////////////////////////////////////////////////
// Make an template based on FFT of sinc squared functions  -- done
void makeTemplate(ffdataStruct *out, candidate *in, REAL4FFTPlan *plan)
{
   
   INT4 ii, jj, numfbins, numffts;
   
   numfbins = (INT4)(roundf(in->fspan*in->Tcoh)+1);   //Number of frequency bins
   numffts = (INT4)floor(2*(in->Tobs/in->Tcoh)-1);     //Number of FFTs
   
   REAL4Vector *psd1 = XLALCreateREAL4Vector((UINT4)(numfbins*numffts));
   INT4Vector *freqbins = XLALCreateINT4Vector((UINT4)numfbins);
   
   REAL4 periodf = 1.0/in->period;
   REAL4 B = in->moddepth*in->Tcoh;
   
   //Create first FFT frequencies
   for (ii=0; ii<out->f->length; ii++) out->f->data[ii] = in->fmin + ii/in->Tcoh;
   
   //Create second FFT frequencies
   for (ii=0; ii<out->fpr->length; ii++) out->fpr->data[ii] = ii/in->Tobs;
   
   //Bin numbers of the frequencies
   for (ii=0; ii<numfbins; ii++) freqbins->data[ii] = (INT4)roundf(in->fmin*in->Tcoh) + ii;
   
   //Determine the signal modulation in bins with time at center of coherence time and create
   //Hann windowed PSDs
   for (ii=0; ii<numffts; ii++) {
      REAL4 t = 0.5*in->Tcoh*ii;
      REAL4 n0 = B*sin(LAL_TWOPI*periodf*t) + in->fsig*in->Tcoh;
      for (jj=0; jj<numfbins; jj++) {
         //Create windowed PSD values
         if ( fabsf(n0-freqbins->data[jj]) <= 5.0 ) psd1->data[ii*numfbins + jj] = 2.0/3.0*in->Tcoh*
            sincxoverxsqminusone(n0-freqbins->data[jj])*sincxoverxsqminusone(n0-freqbins->data[jj]);
         else psd1->data[ii*numfbins + jj] = 0.0;
      }
   }
   
   //Do the second FFT
   REAL4Vector *x = XLALCreateREAL4Vector((UINT4)numffts);
   REAL4Window *win = XLALCreateHannREAL4Window(x->length);
   REAL4 winFactor = 8.0/3.0;
   REAL4Vector *psd = XLALCreateREAL4Vector((UINT4)floor(x->length*0.5)+1);
   REAL4 sum = 0.0;
   INT4 doSecondFFT;
   //First loop over frequencies
   for (ii=0; ii<numfbins; ii++) {
      //Set doSecondFFT check flag to 0. Value becomes 1 if at least one element in frequency row is non-zero
      doSecondFFT = 0;
   
      //Next, loop over times
      for (jj=0; jj<x->length; jj++) {
         //Pick the right frequency bin for each FFT
         x->data[jj] = psd1->data[ii+jj*numfbins];
         
         //Check, do we need to do the second FFT...?
         if (doSecondFFT==0 && x->data[jj]>0.0) doSecondFFT = 1;
         
         //window
         x->data[jj] *= win->data->data[jj];
      }
      
      
      
      //Make the FFT
      INT4 check = 0;
      if (doSecondFFT==1) check = XLALREAL4PowerSpectrum(psd,x,plan);
      if (check != 0) {
         printf("Something wrong with second PSD...\n");
         //fprintf(LOG,"Something wrong with second PSD...\n");
      }
      
      //Scale the data points by 1/N and window factor and (1/fs)
      //Order of vector is by second frequency then first frequency
      if (doSecondFFT==1) {
         for (jj=0; jj<psd->length; jj++) {
            out->ffdata->data[psd->length*ii + jj] = psd->data[jj]*winFactor/x->length*0.5*in->Tcoh;
            sum += out->ffdata->data[psd->length*ii + jj];
         }
      } else {
         for (jj=0; jj<psd->length; jj++) {
            out->ffdata->data[psd->length*ii + jj] = 0.0;
         }
      }
      
   }
   
   //Normalize
   for (ii=0; ii<out->f->length; ii++) {
      for (jj=0; jj<out->fpr->length; jj++) out->ffdata->data[ii*out->fpr->length + jj] /= sum;
   }
   
   
   XLALDestroyREAL4Vector(psd1);
   XLALDestroyINT4Vector(freqbins);
   XLALDestroyREAL4Vector(x);
   XLALDestroyREAL4Window(win);
   XLALDestroyREAL4Vector(psd);
   
}








