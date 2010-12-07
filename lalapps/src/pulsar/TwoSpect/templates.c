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

#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sort.h>

#include <lal/LALConstants.h>
#include <lal/LALMalloc.h>
#include <lal/Window.h>

#include "templates.h"
#include "cdfwchisq.h"

//////////////////////////////////////////////////////////////
// Allocate memory for farStruct struct  -- done
farStruct * new_farStruct(void)
{
   
   const CHAR *fn = __func__;
   
   farStruct *farstruct = XLALMalloc(sizeof(*farstruct));
   if (farstruct==NULL) {
      fprintf(stderr,"%s: XLALMalloc(%lu) failed.\n", fn, sizeof(*farstruct));
      XLAL_ERROR_NULL(fn, XLAL_ENOMEM);
   }
   
   farstruct->topRvalues = NULL;

   return farstruct;

} /* new_farStruct() */


//////////////////////////////////////////////////////////////
// Destroy farStruct struct  -- done
void free_farStruct(farStruct *farstruct)
{
   
   XLALDestroyREAL4Vector(farstruct->topRvalues);
   farstruct->topRvalues = NULL;
   
   XLALFree((farStruct*)farstruct);

} /* free_farStruct() */


//////////////////////////////////////////////////////////////
// Estimate the FAR of the R statistic from the weights
void estimateFAR(farStruct *output, templateStruct *templatestruct, INT4 trials, REAL4 thresh, REAL4Vector *ffplanenoise, REAL4Vector *fbinaveratios)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii, jj;
   INT4 numofweights = 0;
   for (ii=0; ii<(INT4)templatestruct->templatedata->length; ii++) if (templatestruct->templatedata->data[ii]!=0.0) numofweights++;
   
   REAL8 sumofsqweights = 0.0;
   for (ii=0; ii<numofweights; ii++) sumofsqweights += (templatestruct->templatedata->data[ii]*templatestruct->templatedata->data[ii]);
   REAL8 sumofsqweightsinv = 1.0/sumofsqweights;
   
   REAL4Vector *Rs = XLALCreateREAL4Vector((UINT4)trials);
   if (Rs==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, trials);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   
   gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
   if (rng==NULL) {
      fprintf(stderr,"%s: gsl_rng_alloc() failed.\n", fn);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   //srand(time(NULL));
   //UINT8 randseed = rand();
   //gsl_rng_set(rng, randseed);
   gsl_rng_set(rng, 0);
   
   for (ii=0; ii<trials; ii++) {
      //Create noise value and R value
      REAL8 R = 0.0;
      for (jj=0; jj<numofweights; jj++) {
         REAL8 noise = expRandNum(ffplanenoise->data[ templatestruct->secondfftfrequencies->data[jj] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[jj] ], rng);
         R += (noise - ffplanenoise->data[ templatestruct->secondfftfrequencies->data[jj] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[jj] ])*templatestruct->templatedata->data[jj];
      }
      Rs->data[ii] = (REAL4)(R*sumofsqweightsinv);
   } /* for ii < trials */
   REAL4 mean = calcMean(Rs);
   REAL4 sigma = calcStddev(Rs);
   if (XLAL_IS_REAL4_FAIL_NAN(mean)) {
      fprintf(stderr,"%s: calcMean() failed.\n", fn);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (XLAL_IS_REAL4_FAIL_NAN(sigma)) {
      fprintf(stderr,"%s: calcStddev() failed.\n", fn);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   
   //Do an insertion sort. At best this is O(thresh*trials), at worst this is O(thresh*trials*trials).
   if (output->topRvalues == NULL) {
      output->topRvalues = XLALCreateREAL4Vector((UINT4)roundf(thresh*trials)+1);
      if (output->topRvalues==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, (INT4)roundf(thresh*trials)+1);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
   }
   if ((gsl_sort_float_largest((float*)output->topRvalues->data, output->topRvalues->length, (float*)Rs->data, 1, Rs->length)) != 0) {
      fprintf(stderr,"%s: gsl_sort_float_largest() failed.\n", fn);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   
   output->far = output->topRvalues->data[output->topRvalues->length - 1];
   output->distMean = mean;
   output->distSigma = sigma;
   
   //Destroy
   XLALDestroyREAL4Vector(Rs);
   gsl_rng_free(rng);

} /* estimateFAR() */


//////////////////////////////////////////////////////////////
// Numerically solve for the FAR of the R statistic from the weights
void numericFAR(farStruct *output, templateStruct *templatestruct, REAL4 thresh, REAL4Vector *ffplanenoise, REAL4Vector *fbinaveratios)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii;
   
   INT4 numweights = 0;
   for (ii=0; ii<(INT4)templatestruct->templatedata->length; ii++) if (templatestruct->templatedata->data[ii]!=0) numweights++;
   
   REAL8 sumwsq = 0.0;
   for (ii=0; ii<numweights; ii++) sumwsq += templatestruct->templatedata->data[ii]*templatestruct->templatedata->data[ii];
   
   INT4 errcode = 0;
   
   //Set up solver
   const gsl_root_fdfsolver_type *T = gsl_root_fdfsolver_newton;
   gsl_root_fdfsolver *s = gsl_root_fdfsolver_alloc(T);
   gsl_function_fdf FDF;
   
   //Include the various parameters in the struct required by GSL
   struct gsl_probR_pars params = {templatestruct, ffplanenoise, fbinaveratios, thresh, errcode};
   
   //Assign GSL function the necessary parts
   FDF.f = &gsl_probR;
   FDF.df = &gsl_dprobRdR;
   FDF.fdf = &gsl_probRandDprobRdR;
   FDF.params = &params;
   
   //Start off with an initial guess
   REAL8 root = 10.0;
   REAL8 initialroot = root;
   
   //Set the solver at the beginning
   if ((gsl_root_fdfsolver_set(s, &FDF, initialroot)) != 0) {
      fprintf(stderr,"%s: Unable to initialize root solver to first guess.\n", fn);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   
   //And now find the root
   ii = 0;
   INT4 max_iter = 100;
   INT4 status = GSL_CONTINUE;
   REAL8 prevroot;
   while (status==GSL_CONTINUE && ii<max_iter) {
      ii++;
      status = gsl_root_fdfsolver_iterate(s);
      if (status!=GSL_CONTINUE && status!=GSL_SUCCESS) {
         fprintf(stderr,"%s: gsl_root_fdfsolver_iterate() failed with code %d.\n", fn, status);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
      prevroot = root;
      root = gsl_root_fdfsolver_root(s);
      status = gsl_root_test_delta(prevroot, root, 0.0, 0.002);
      if (status!=GSL_CONTINUE && status!=GSL_SUCCESS) {
         fprintf(stderr,"%s: gsl_root_test_delta() failed with code %d.\n", fn, status);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
   } /* while status==GSL_CONTINUE && ii < max_iter */
   
   if (status != GSL_SUCCESS) {
      fprintf(stderr,"%s: Root finding iteration (%d/%d) failed with failure code %d. Previous root = %f, current root = %f\n", fn, ii, max_iter, status, prevroot, root);
      XLAL_ERROR_VOID(fn, XLAL_FAILURE);
   } else if (ii==max_iter) {
      fprintf(stderr,"%s: Root finding iteration (%d/%d) failed with failure code %d. Previous root = %f, current root = %f\n", fn, ii, max_iter, status, prevroot, root);
      XLAL_ERROR_VOID(fn, XLAL_EMAXITER);
   } else if (root<=0.0) {
      fprintf(stderr,"%s: Threshold value found (%f) is less than 0.0!\n", fn, root);
      XLAL_ERROR_VOID(fn, XLAL_ERANGE);
   }
   
   output->far = root;
   output->distMean = 0.0;
   output->distSigma = 0.0; //Fake the value of sigma
   output->farerrcode = errcode;
   
   //Cleanup
   gsl_root_fdfsolver_free(s);
   
} /* numericFAR() */
REAL8 gsl_probR(REAL8 R, void *param)
{
   
   struct gsl_probR_pars *pars = (struct gsl_probR_pars*)param;
   
   REAL8 prob = probR(pars->templatestruct, pars->ffplanenoise, pars->fbinaveratios, R, &pars->errcode);
   
   REAL8 returnval = prob - log10(pars->threshold);
   
   return returnval;
   
} /* gsl_probR() */
REAL8 gsl_dprobRdR(REAL8 R, void *param)
{
   
   struct gsl_probR_pars *pars = (struct gsl_probR_pars*)param;
   
   REAL8 dR = 0.001;
   REAL8 slope;
   
   //Explicit computation of slope
   REAL8 R1 = (1.0+dR)*R;
   REAL8 R2 = (1.0-dR)*R;
   REAL8 prob1 = probR(pars->templatestruct, pars->ffplanenoise, pars->fbinaveratios, R1, &pars->errcode);
   REAL8 prob2 = probR(pars->templatestruct, pars->ffplanenoise, pars->fbinaveratios, R2, &pars->errcode);
   while (fabs(prob1-prob2)<10.0*LAL_REAL4_MIN) {
      dR *= 2.0;
      R1 = (1.0+dR)*R;
      R2 = (1.0-dR)*R;
      prob1 = probR(pars->templatestruct, pars->ffplanenoise, pars->fbinaveratios, R1, &pars->errcode);
      prob2 = probR(pars->templatestruct, pars->ffplanenoise, pars->fbinaveratios, R2, &pars->errcode);
   }
   slope = (prob1-prob2)/(R1-R2);
   
   //One more point
   dR *= 2.0;
   REAL8 R3 = (1.0-dR)*R;
   REAL8 prob3 = probR(pars->templatestruct, pars->ffplanenoise, pars->fbinaveratios, R3, &pars->errcode);
   
   //Compute average
   slope = (slope + (prob3-prob1)/(R3-R1) + (prob2-prob3)/(R2-R3))/3.0;
   
   return slope;
   
} /* gsl_dprobRdR() */
void gsl_probRandDprobRdR(REAL8 R, void *param, REAL8 *probabilityR, REAL8 *dprobRdR)
{
   
   struct gsl_probR_pars *pars = (struct gsl_probR_pars*)param;
   
   *probabilityR = gsl_probR(R, pars);
   
   *dprobRdR = gsl_dprobRdR(R, pars);
   
} /* gsl_probRandDprobRdR() */


//////////////////////////////////////////////////////////////
// Analytically calculate the probability of a true signal output is log10(prob)
REAL8 probR(templateStruct *templatestruct, REAL4Vector *ffplanenoise, REAL4Vector *fbinaveratios, REAL8 R, INT4 *errcode)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii;
   REAL8 prob = 0.0;
   REAL8 sumwsq = 0.0;
   INT4 numweights = 0;
   for (ii=0; ii<(INT4)templatestruct->templatedata->length; ii++) {
      if (templatestruct->templatedata->data[ii]!=0.0) numweights++;
      sumwsq += templatestruct->templatedata->data[ii]*templatestruct->templatedata->data[ii];
   }
   
   REAL8Vector *newweights = XLALCreateREAL8Vector((UINT4)numweights);
   REAL8Vector *noncentrality = XLALCreateREAL8Vector((UINT4)numweights);
   INT4Vector *dofs = XLALCreateINT4Vector((UINT4)numweights);
   INT4Vector *sorting = XLALCreateINT4Vector((UINT4)numweights);
   if (newweights==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL8Vector(%d) failed.\n", fn, numweights);
      XLAL_ERROR_REAL8(fn, XLAL_EFUNC);
   } else if (noncentrality==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL8Vector(%d) failed.\n", fn, numweights);
      XLAL_ERROR_REAL8(fn, XLAL_EFUNC);
   } else if (dofs==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", fn, numweights);
      XLAL_ERROR_REAL8(fn, XLAL_EFUNC);
   } else if (sorting==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", fn, numweights);
      XLAL_ERROR_REAL8(fn, XLAL_EFUNC);
   }
   
   REAL8 Rpr = R;
   for (ii=0; ii<(INT4)newweights->length; ii++) {
      newweights->data[ii] = 0.5*templatestruct->templatedata->data[ii]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[ii] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[ii] ]/sumwsq;
      noncentrality->data[ii] = 0.0;
      dofs->data[ii] = 2;
      Rpr += templatestruct->templatedata->data[ii]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[ii] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[ii] ]/sumwsq;
   }
   
   //INT4 errcode;
   qfvars vars;
   vars.weights = newweights;
   vars.noncentrality = noncentrality;
   vars.dofs = dofs;
   vars.sorting = sorting;
   vars.lim = 10000;
   vars.c = Rpr;
   REAL8 sigma = 0.0;
   REAL8 accuracy = 1.0e-5;   //don't change this value
   
   //cdfwchisq(algorithm variables, sigma, accuracy, error code)
   prob = 1.0 - cdfwchisq(&vars, sigma, accuracy, errcode);
   
   //Large R values can cause a problem when computing the probability. We run out of accuracy quickly even using double precision
   //Potential fix: compute log10(prob) for smaller values of R, for when slope is linear between log10 probabilities
   //Use slope to extend the computation and then compute the exponential of the found log10 probability.
   REAL8 c1, c2, c3, logprob1, logprob2, logprob3, probslope, logprobest;
   INT4 estimatedTheProb = 0;
   if (prob<=1.0e-4) {
      estimatedTheProb = 1;
      
      INT4 errcode1 = 0, errcode2 = 0, errcode3 = 0;
      
      c1 = 0.95*vars.c;
      vars.c = c1;
      REAL8 tempprob = 1.0-cdfwchisq(&vars, sigma, accuracy, &errcode1);
      while (tempprob<=1.0e-4) {
         c1 *= 0.95;
         vars.c = c1;
         tempprob = 1.0-cdfwchisq(&vars, sigma, accuracy, &errcode1);
      }
      logprob1 = log10(tempprob);
      
      c2 = 0.95*c1;
      vars.c = c2;
      REAL8 tempprob2 = 1.0-cdfwchisq(&vars, sigma, accuracy, &errcode2);
      while (tempprob2-tempprob <= 2.0e-4) {
         c2 *= 0.95;
         vars.c = c2;
         tempprob2 = 1.0-cdfwchisq(&vars, sigma, accuracy, &errcode2);
      }
      logprob2 = log10(tempprob2);
      
      c3 = 0.95*c2;
      vars.c = c3;
      REAL8 tempprob3 = 1.0-cdfwchisq(&vars, sigma, accuracy, &errcode3);
      while (tempprob3-tempprob2 <= 2.0e-4) {
         c3 *= 0.95;
         vars.c = c3;
         tempprob3 = 1.0-cdfwchisq(&vars, sigma, accuracy, &errcode3);
      }
      logprob3 = log10(tempprob3);
      
      //If either point along the slope had a problem at the end, then better fail.
      //Otherwise, set errcode = 0;
      if (errcode1!=0 || errcode2!=0 || errcode3!=0) {
         fprintf(stderr,"%s: cdfwchisq() failed.\n", fn);
         XLAL_ERROR_REAL8(fn, XLAL_EFUNC);
      } else {
         *errcode = errcode1;
      }
      
      //Calculating slope
      probslope = ((logprob1-logprob2)/(c1-c2) + (logprob1-logprob3)/(c1-c3) + (logprob2-logprob3)/(c2-c3))/3.0;
      
      //Find the log10(prob) of the original Rpr value
      logprobest = probslope*(Rpr-c2) + logprob2;
      
   } /* if prob <= 1e-4 */
   
   //If errcode is still 0, better fail
   if (*errcode!=0) {
      fprintf(stderr,"%s: cdfwchisq() failed.\n", fn);
      XLAL_ERROR_REAL8(fn, XLAL_EFUNC);
   }
   
   //Cleanup
   XLALDestroyREAL8Vector(newweights);
   XLALDestroyREAL8Vector(noncentrality);
   XLALDestroyINT4Vector(dofs);
   XLALDestroyINT4Vector(sorting);
   
   //return prob;
   if (estimatedTheProb==1) return logprobest;
   else return log10(prob);
   
} /* probR() */


templateStruct * new_templateStruct(INT4 length)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii;
   
   templateStruct *templatestruct = XLALMalloc(sizeof(*templatestruct));
   if (templatestruct==NULL) {
      fprintf(stderr,"%s: XLALMalloc(%lu) failed.\n", fn, sizeof(*templatestruct));
      XLAL_ERROR_NULL(fn, XLAL_ENOMEM);
   }
   
   templatestruct->templatedata = XLALCreateREAL4Vector((UINT4)length);
   templatestruct->pixellocations = XLALCreateINT4Vector((UINT4)length);
   templatestruct->firstfftfrequenciesofpixels = XLALCreateINT4Vector((UINT4)length);
   templatestruct->secondfftfrequencies = XLALCreateINT4Vector((UINT4)length);
   if (templatestruct->templatedata==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, length);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   } else if (templatestruct->pixellocations==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", fn, length);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   } else if (templatestruct->firstfftfrequenciesofpixels==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", fn, length);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   } else if (templatestruct->secondfftfrequencies==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", fn, length);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   } 
   
   for (ii=0; ii<length; ii++) {
      templatestruct->templatedata->data[ii] = 0.0;
      templatestruct->pixellocations->data[ii] = 0;
      templatestruct->firstfftfrequenciesofpixels->data[ii] = 0;
      templatestruct->secondfftfrequencies->data[ii] = 0;
   }
   
   return templatestruct;
   
} /* new_templateStruct() */


void free_templateStruct(templateStruct *nameoftemplate)
{
   
   XLALDestroyREAL4Vector(nameoftemplate->templatedata);
   XLALDestroyINT4Vector(nameoftemplate->pixellocations);
   XLALDestroyINT4Vector(nameoftemplate->firstfftfrequenciesofpixels);
   XLALDestroyINT4Vector(nameoftemplate->secondfftfrequencies);
   
   XLALFree((templateStruct*)nameoftemplate);
   
} /* free_templateStruct() */


//////////////////////////////////////////////////////////////
// Make an estimated template based on FFT of train of Gaussians
//void makeTemplateGaussians(ffdataStruct *out, candidate *in)
void makeTemplateGaussians(templateStruct *output, candidate input, inputParamsStruct *params)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii, jj, numfbins, numffts, N;
   
   numfbins = (INT4)(round(params->fspan*params->Tcoh)+1);   //Number of frequency bins
   numffts = (INT4)floor(params->Tobs/(params->Tcoh-params->SFToverlap)-1);   //Number of FFTs
   N = (INT4)floor(params->Tobs/input.period);     //Number of Gaussians
   
   REAL8 periodf = 1.0/input.period;
   
   //Determine separation in time of peaks for each frequency
   REAL8Vector *phi_actual = XLALCreateREAL8Vector((UINT4)numfbins);
   if (phi_actual==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL8Vector(%d) failed.\n", fn, numfbins);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   for (ii=0; ii<(INT4)phi_actual->length; ii++) {
      if ( fabs(params->fmin + ii/params->Tcoh - input.fsig)/input.moddepth <= 1.0 ) {
         phi_actual->data[ii] = 0.5*input.period - asin(fabs(params->fmin + ii/params->Tcoh - input.fsig)/
            input.moddepth)*LAL_1_PI*input.period;
      } else {
         phi_actual->data[ii] = 0.0;
      }
   } /* for ii < phi_actual->length */
   
   //Create second FFT frequencies
   REAL8Vector *fpr = XLALCreateREAL8Vector((UINT4)floor(numffts*0.5)+1);
   if (fpr==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL8Vector(%d) failed.\n", fn, (INT4)floor(numffts*0.5)+1);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   for (ii=0; ii<(INT4)fpr->length; ii++) fpr->data[ii] = (REAL8)ii/params->Tobs;
   
   //Scale used for "spillover" into bins outside of phi_actual
   REAL8 k = input.moddepth*params->Tcoh;    //amplitude of modulation in units of bins
   REAL8Vector *scale = XLALCreateREAL8Vector((UINT4)numfbins);      //the scaling factor
   if (scale==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL8Vector(%d) failed.\n", fn, numfbins);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   INT4 m0 = (INT4)round(input.fsig*params->Tcoh) - (INT4)round(params->fmin*params->Tcoh);   //central frequency bin
   INT4 mextent = (INT4)floor(input.moddepth*params->Tcoh);   //Bins filled by modulation
   REAL8 overage = (k-(REAL8)mextent)-1.0;
   INT4 fnumstart = -1;
   INT4 fnumend = -1;
   for (ii=0; ii<(INT4)scale->length; ii++) {
      if (mextent != 0) {
         if (ii < m0-mextent-2 || ii > m0+mextent+2) {
            scale->data[ii] = 0.0;
         } else if (ii == m0-mextent-2 || ii == m0+mextent+2) {
            scale->data[ii] = sincxoverxsqminusone(overage-1)*sincxoverxsqminusone(overage-1);
         } else if (ii == m0-mextent-1 || ii == m0+mextent+1) {
            scale->data[ii] = sincxoverxsqminusone(overage)*sincxoverxsqminusone(overage);
         } else {
            scale->data[ii] = 1.0;
         }
      } else {
         if (ii < m0-2 || ii > m0+2) {
            scale->data[ii] = 0.0;
         } else if (ii == m0-2 || ii == m0+2) {
            scale->data[ii] = sincxoverxsqminusone(overage-1)*sincxoverxsqminusone(overage-1);
         } else if (ii == m0-1 || ii == m0+1) {
            scale->data[ii] = sincxoverxsqminusone(overage)*sincxoverxsqminusone(overage);
         } else {
            scale->data[ii] = 1.0;
         }
      } /* if mextent!=0 else ... */
   } /* for ii < scale->length */
   for (ii=0; ii<(INT4)scale->length; ii++) {
      if (scale->data[ii] != 0.0 && fnumstart == -1) fnumstart = ii;
      if (scale->data[ii] == 0.0 && fnumstart != -1 && fnumend==-1) fnumend = ii-1;
   }
   
   //Exit with failure if there is a problem
   if (fnumend==-1) {
      fprintf(stderr, "%s: Failed because fnumend was -1.\n", fn);
      XLAL_ERROR_VOID(fn, XLAL_EFAILED);
   }
   
   //Make sigmas for each frequency
   REAL8Vector *sigmas = XLALCreateREAL8Vector((UINT4)(fnumend-fnumstart+1));
   REAL8Vector *wvals = XLALCreateREAL8Vector((UINT4)floor(2.0*input.period/params->Tcoh));
   REAL8Vector *allsigmas = XLALCreateREAL8Vector(wvals->length * sigmas->length);
   if (sigmas==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL8Vector(%d) failed.\n", fn, (UINT4)(fnumend-fnumstart+1));
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (wvals==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL8Vector(%d) failed.\n", fn, (UINT4)floor(2.0*input.period/params->Tcoh));
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (allsigmas==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL8Vector(%d) failed.\n", fn, wvals->length * sigmas->length);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   for (ii=0; ii<(INT4)wvals->length; ii++) {         //t = (ii+1)*in->Tcoh*0.5
      REAL8 sigbin = (input.moddepth*cos(LAL_TWOPI*periodf*((ii+1)*params->Tcoh*0.5))+input.fsig)*params->Tcoh;
      REAL8 sigbinvelocity = fabs(-input.moddepth*sin(LAL_TWOPI*periodf*((ii+1)*params->Tcoh*0.5))*params->Tcoh*0.5*params->Tcoh*LAL_TWOPI*periodf);
      REAL8 sigma = 0.5 * params->Tcoh * ((383.85*LAL_1_PI)*(0.5*6.1e-3) / ((sigbinvelocity+0.1769)*(sigbinvelocity+0.1769)+(0.5*6.1e-3)*(0.5*6.1e-3)) + 0.3736);   //Derived fit from simulation
      for (jj=0; jj<(INT4)sigmas->length; jj++) {
         allsigmas->data[ii*sigmas->length + jj] = sincxoverxsqminusone(sigbin-round(params->fmin*params->Tcoh+jj+fnumstart))*sincxoverxsqminusone(sigbin-round(params->fmin*params->Tcoh+jj+fnumstart))*sigma;
      }
   } /* for ii < wvals->length */
   for (ii=0; ii<(INT4)sigmas->length; ii++) {
      for (jj=0; jj<(INT4)wvals->length; jj++) wvals->data[jj] = allsigmas->data[ii + jj*sigmas->length]*allsigmas->data[ii + jj*sigmas->length];
      sigmas->data[ii] = sqrt(calcMeanD(wvals));
   }
   
   
   //Create template
   REAL8 sum = 0.0;
   REAL8 dataval;
   for (ii=0; ii<(INT4)sigmas->length; ii++) {
      REAL8 s = sigmas->data[ii];
      REAL8 scale1 = 1.0/(1.0+exp(-phi_actual->data[ii+fnumstart]*phi_actual->data[ii+fnumstart]*0.5/(s*s)));
      for (jj=0; jj<(INT4)fpr->length; jj++) {
         
         if (jj==0 || jj==1) {
            dataval = 0.0;
         } else if (fabs(cos(input.period*LAL_TWOPI*fpr->data[jj])-1.0)<1e-5) {
            dataval = scale->data[ii+fnumstart] * scale1 * 2.0 * LAL_TWOPI * s * s * exp(-s * s * LAL_TWOPI * LAL_TWOPI * fpr->data[jj] * fpr->data[jj]) * (cos(phi_actual->data[ii+fnumstart] * LAL_TWOPI * fpr->data[jj]) + 1.0) * N * N;
         } else {
            dataval = scale->data[ii+fnumstart] * scale1 * 2.0 * LAL_TWOPI * s * s * exp(-s * s * LAL_TWOPI * LAL_TWOPI * fpr->data[jj] * fpr->data[jj]) * (cos(N * input.period * LAL_TWOPI * fpr->data[jj]) - 1.0) * (cos(phi_actual->data[ii+fnumstart] * LAL_TWOPI * fpr->data[jj]) + 1.0) / (cos(input.period * LAL_TWOPI * fpr->data[jj]) - 1.0);
         }
         
         //Set any bin below 1e-7 to 0.0
         if (dataval <= 1.0e-7) dataval = 0.0;
         
         //Sum up the weights in total
         sum += dataval;
         
         //Compare with weakest top bins and if larger, launch a search to find insertion spot (insertion sort)
         if (jj>1 && dataval > output->templatedata->data[output->templatedata->length-1]) {
            insertionSort_template(output, (REAL4)dataval, (ii+fnumstart)*fpr->length+jj, ii+fnumstart, jj);
         }
      } /* for jj < fpr->legnth */
   } /* for ii < sigmas->length */
   
   //Normalize
   for (ii=0; ii<(INT4)output->templatedata->length; ii++) if (output->templatedata->data[ii]!=0.0) output->templatedata->data[ii] /= (REAL4)sum;
   
   //Destroy variables
   XLALDestroyREAL8Vector(phi_actual);
   XLALDestroyREAL8Vector(scale);
   XLALDestroyREAL8Vector(sigmas);
   XLALDestroyREAL8Vector(allsigmas);
   XLALDestroyREAL8Vector(wvals);
   XLALDestroyREAL8Vector(fpr);

} /* mateTemplateGaussians() */


//////////////////////////////////////////////////////////////
// Make an template based on FFT of sinc squared functions  -- done
void makeTemplate(templateStruct *output, candidate input, inputParamsStruct *params, REAL4FFTPlan *plan)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii, jj, numfbins, numffts;
   
   numfbins = (INT4)(round(params->fspan*params->Tcoh)+1);   //Number of frequency bins
   numffts = (INT4)floor(params->Tobs/(params->Tcoh-params->SFToverlap)-1);   //Number of FFTs
   
   REAL4Vector *psd1 = XLALCreateREAL4Vector((UINT4)(numfbins*numffts));
   INT4Vector *freqbins = XLALCreateINT4Vector((UINT4)numfbins);
   if (psd1==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numfbins*numffts);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (freqbins==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", fn, numfbins);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   
   REAL4 periodf = 1.0/input.period;
   REAL4 B = input.moddepth*params->Tcoh;
   
   //Bin numbers of the frequencies
   for (ii=0; ii<numfbins; ii++) freqbins->data[ii] = (INT4)roundf(params->fmin*params->Tcoh) + ii;
   
   //Determine the signal modulation in bins with time at center of coherence time and create
   //Hann windowed PSDs
   for (ii=0; ii<numffts; ii++) {
      REAL4 t = 0.5*params->Tcoh*ii;
      REAL4 n0 = B*sin(LAL_TWOPI*periodf*t) + input.fsig*params->Tcoh;
      for (jj=0; jj<numfbins; jj++) {
         //Create windowed PSD values
         if ( fabs(n0-freqbins->data[jj]) <= 5.0 ) psd1->data[ii*numfbins + jj] = 2.0/3.0*params->Tcoh*sincxoverxsqminusone(n0-freqbins->data[jj])*sincxoverxsqminusone(n0-freqbins->data[jj]);
         else psd1->data[ii*numfbins + jj] = 0.0;
      } /* for jj < numfbins */
   } /* for ii < numffts */
   
   //Do the second FFT
   REAL4Vector *x = XLALCreateREAL4Vector((UINT4)numffts);
   REAL4Window *win = XLALCreateHannREAL4Window(x->length);
   REAL4Vector *psd = XLALCreateREAL4Vector((UINT4)floor(x->length*0.5)+1);
   if (x==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numffts);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (win==NULL) {
      fprintf(stderr,"%s: XLALCreateHannREAL4Window(%d) failed.\n", fn, x->length);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   } else if (psd==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", fn, (UINT4)floor(x->length*0.5)+1);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   REAL4 winFactor = 8.0/3.0;
   REAL8 sum = 0.0;
   INT4 doSecondFFT;
   //First loop over frequencies
   for (ii=0; ii<numfbins; ii++) {
      //Set doSecondFFT check flag to 0. Value becomes 1 if at least one element in frequency row is non-zero
      doSecondFFT = 0;
   
      //Next, loop over times
      for (jj=0; jj<(INT4)x->length; jj++) {
         //Check, do we need to do the second FFT...?
         if (doSecondFFT==0 && psd1->data[ii+jj*numfbins]>0.0) {
            doSecondFFT = 1;
            jj = (INT4)x->length;  //End the search for bins with power
         }
      }
      
      //If there was power in the frequency bin of the template, then do the FFT
      if (doSecondFFT==1) {
         //Obtain and window the time series
         x = SSVectorMultiply_with_stride_and_offset(x, psd1, win->data, numfbins, 1, ii, 0);
         if (xlalErrno!=0) {
            fprintf(stderr,"%s, SSVectorMultiply_with_stride_and_offset() failed.\n", fn);
            XLAL_ERROR_VOID(fn, XLAL_EFUNC);
         }
         
         //Do the FFT
         if (XLALREAL4PowerSpectrum(psd,x,plan) != 0) {
            fprintf(stderr,"%s: XLALREAL4PowerSpectrum() failed.\n", fn);
            XLAL_ERROR_VOID(fn, XLAL_EFUNC);
         }
      } /* if doSecondFFT */
      
      //Scale the data points by 1/N and window factor and (1/fs)
      //Order of vector is by second frequency then first frequency
      //Ignore the DC and 1st frequency bins
      if (doSecondFFT==1) {
         for (jj=0; jj<(INT4)psd->length; jj++) {
            REAL4 correctedValue = (REAL4)(psd->data[jj]*winFactor/x->length*0.5*params->Tcoh);
            if (correctedValue<=1.0e-7) correctedValue = 0.0;
            
            //Sum the total weights
            sum += correctedValue;
            
            //Sort the weights, insertion sort technique
            if (jj>1 && correctedValue > output->templatedata->data[output->templatedata->length-1]) {
               insertionSort_template(output, correctedValue, ii*psd->length+jj, ii, jj);
            }
         } /* for jj < psd->length */
      } /* if(doSecondFFT) */
      
   } /* if ii < numfbins */
   
   //Normalize
   for (ii=0; ii<(INT4)output->templatedata->length; ii++) if (output->templatedata->data[ii]!=0.0) output->templatedata->data[ii] /= (REAL4)sum;
   
   //Destroy
   XLALDestroyREAL4Vector(psd1);
   XLALDestroyINT4Vector(freqbins);
   XLALDestroyREAL4Vector(x);
   XLALDestroyREAL4Window(win);
   XLALDestroyREAL4Vector(psd);
   
}


//////////////////////////////////////////////////////////////
// Does the insertion sort for the template weights
void insertionSort_template(templateStruct *output, REAL4 weight, INT4 pixelloc, INT4 firstfftfreq, INT4 secfftfreq)
{
   
   INT4 ii;
   
   INT4 insertionpoint = (INT4)output->templatedata->length-1;
   while (insertionpoint > 0 && weight > output->templatedata->data[insertionpoint-1]) insertionpoint--;
   
   for (ii=output->templatedata->length-1; ii>insertionpoint; ii--) {
      output->templatedata->data[ii] = output->templatedata->data[ii-1];
      output->pixellocations->data[ii] = output->pixellocations->data[ii-1];
      output->firstfftfrequenciesofpixels->data[ii] = output->firstfftfrequenciesofpixels->data[ii-1];
      output->secondfftfrequencies->data[ii] = output->secondfftfrequencies->data[ii-1];
   }
   
   output->templatedata->data[insertionpoint] = weight;
   output->pixellocations->data[insertionpoint] = pixelloc;
   output->firstfftfrequenciesofpixels->data[insertionpoint] = firstfftfreq;
   output->secondfftfrequencies->data[insertionpoint] = secfftfreq;
   
} /* insertionSort_template() */



//////////////////////////////////////////////////////////////
// Calculates y = sin(pi*x)/(pi*x)/(x^2-1)
REAL8 sincxoverxsqminusone(REAL8 x)
{
   
   REAL8 val;
   
   if (fabs(x*x-1.0)<1.0e-5) val = -0.5;
   else val = gsl_sf_sinc(x)/(x*x-1.0);
   
   return val;
   
} /* sincxoverxsqminusone() */



