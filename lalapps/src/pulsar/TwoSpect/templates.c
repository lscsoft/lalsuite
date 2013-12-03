/*
*  Copyright (C) 2010, 2011, 2012, 2013 Evan Goetz
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

#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_sort.h>

#include <lal/LALConstants.h>
#include <lal/LALMalloc.h>
#include <lal/Window.h>
#include <lal/VectorOps.h>

#include "templates.h"
#include "cdfwchisq.h"
#include "statistics.h"
#include "candidates.h"
#include "vectormath.h"


//////////////////////////////////////////////////////////////
// Allocate memory for farStruct struct
farStruct * new_farStruct(void)
{
   farStruct *farstruct = XLALMalloc(sizeof(*farstruct));
   if (farstruct==NULL) {
      fprintf(stderr,"%s: XLALMalloc(%zu) failed.\n", __func__, sizeof(*farstruct));
      XLAL_ERROR_NULL(XLAL_ENOMEM);
   }
   farstruct->far = 1.0;
   farstruct->topRvalues = NULL;
   return farstruct;
} /* new_farStruct() */


//////////////////////////////////////////////////////////////
// Destroy farStruct struct
void free_farStruct(farStruct *farstruct)
{
   XLALDestroyREAL4Vector(farstruct->topRvalues);
   farstruct->topRvalues = NULL;
   XLALFree((farStruct*)farstruct);
} /* free_farStruct() */


//////////////////////////////////////////////////////////////
// Estimate the FAR of the R statistic from the weights
// We do this by a number of trials
void estimateFAR(farStruct *output, templateStruct *templatestruct, INT4 trials, REAL8 thresh, REAL4Vector *ffplanenoise, REAL4Vector *fbinaveratios)
{
   
   INT4 ii, jj;
   INT4 numofweights = 0;
   for (ii=0; ii<(INT4)templatestruct->templatedata->length; ii++) if (templatestruct->templatedata->data[ii]!=0.0) numofweights++;
   
   REAL8 sumofsqweights = 0.0;
   for (ii=0; ii<numofweights; ii++) sumofsqweights += (templatestruct->templatedata->data[ii]*templatestruct->templatedata->data[ii]);
   REAL8 sumofsqweightsinv = 1.0/sumofsqweights;
   
   REAL4Vector *Rs = XLALCreateREAL4Vector(trials);
   if (Rs==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, trials);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   
   gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
   if (rng==NULL) {
      fprintf(stderr,"%s: gsl_rng_alloc() failed.\n", __func__);
      XLAL_ERROR_VOID(XLAL_EFUNC);
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
      fprintf(stderr,"%s: calcMean() failed.\n", __func__);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (XLAL_IS_REAL4_FAIL_NAN(sigma)) {
      fprintf(stderr,"%s: calcStddev() failed.\n", __func__);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   
   //Do an insertion sort. At best this is O(thresh*trials), at worst this is O(thresh*trials*trials).
   if (output->topRvalues == NULL) {
      output->topRvalues = XLALCreateREAL4Vector((INT4)round(thresh*trials)+1);
      if (output->topRvalues==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, (INT4)roundf(thresh*trials)+1);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      }
   }
   if ((gsl_sort_float_largest((float*)output->topRvalues->data, output->topRvalues->length, (float*)Rs->data, 1, Rs->length)) != 0) {
      fprintf(stderr,"%s: gsl_sort_float_largest() failed.\n", __func__);
      XLAL_ERROR_VOID(XLAL_EFUNC);
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
// This is done using the Davies algorithm and a root finding algorithm
// method = 0: Brent's method
// method = 1: Newton's method
void numericFAR(farStruct *output, templateStruct *templatestruct, REAL8 thresh, REAL4Vector *ffplanenoise, REAL4Vector *fbinaveratios, inputParamsStruct *inputParams, INT4 method)
{

   INT4 ii;
   INT4 errcode = 0;
   
   //Set up solver: method 0 is Brent's method, method 1 is Newton's method
   const gsl_root_fsolver_type *T1 = gsl_root_fsolver_brent;
   gsl_root_fsolver *s1 = gsl_root_fsolver_alloc(T1);
   if (s1==NULL) {
      fprintf(stderr,"%s: gsl_root_fsolver_alloc() failed.\n", __func__);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   gsl_function F;
   const gsl_root_fdfsolver_type *T0 = gsl_root_fdfsolver_newton;
   gsl_root_fdfsolver *s0 = gsl_root_fdfsolver_alloc(T0);
   if (s0==NULL) {
      fprintf(stderr,"%s: gsl_root_fdfsolver_alloc() failed.\n", __func__);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   gsl_function_fdf FDF;
   
   
   //Include the various parameters in the struct required by GSL
   struct gsl_probR_pars params = {templatestruct, ffplanenoise, fbinaveratios, thresh, inputParams, errcode};
   
   //Assign GSL function the necessary parts
   if (method != 0) {
      F.function = &gsl_probR;
      F.params = &params;
   } else {
      FDF.f = &gsl_probR;
      FDF.df = &gsl_dprobRdR;
      FDF.fdf = &gsl_probRandDprobRdR;
      FDF.params = &params;
   }
   
   //Start off with an initial guess and set the solver at the beginning
   REAL8 Rlow = 0.0, Rhigh = 10000.0, root = 400.0;
   if (method != 0) {
      if ( (gsl_root_fsolver_set(s1, &F, Rlow, Rhigh)) != 0 ) {
         fprintf(stderr,"%s: Unable to initialize root solver to bracketed positions.\n", __func__);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      }
   } else {
      if ( (gsl_root_fdfsolver_set(s0, &FDF, root)) != 0 ) {
         fprintf(stderr,"%s: Unable to initialize root solver to first guess.\n", __func__);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      } 
   }
   
   //And now find the root
   ii = 0;
   INT4 max_iter = 100, jj = 0, max_retries = 10;
   INT4 status = GSL_CONTINUE;
   REAL8 prevroot = 0.0;
   while (status==GSL_CONTINUE && ii<max_iter) {
      
      ii++;
      
      if (method != 0) {
         status = gsl_root_fsolver_iterate(s1);
         if (status!=GSL_CONTINUE && status!=GSL_SUCCESS) {
            fprintf(stderr,"%s: gsl_root_fsolver_iterate() failed with code %d.\n", __func__, status);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
         if (ii>0) prevroot = root;
         root = gsl_root_fsolver_root(s1);
         Rlow = gsl_root_fsolver_x_lower(s1);
         Rhigh = gsl_root_fsolver_x_upper(s1);
         status = gsl_root_test_interval(Rlow, Rhigh, 0.0, 0.001);
         if (status!=GSL_CONTINUE && status!=GSL_SUCCESS) {
            fprintf(stderr,"%s: gsl_root_test_interval() failed with code %d.\n", __func__, status);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
      } else {
         status = gsl_root_fdfsolver_iterate(s0);
         if (status!=GSL_CONTINUE && status!=GSL_SUCCESS) {
            fprintf(stderr,"%s: gsl_root_fdfsolver_iterate() failed with code %d.\n", __func__, status);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
         prevroot = root;
         root = gsl_root_fdfsolver_root(s0);
         status = gsl_root_test_delta(prevroot, root, 0.0, 0.001);
         if (status!=GSL_CONTINUE && status!=GSL_SUCCESS) {
            fprintf(stderr,"%s: gsl_root_test_delta() failed with code %d.\n", __func__, status);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
         
         //If there is an issue that the root is negative, try a new initial guess
         if (root<0.0 && jj<max_retries) {
            ii = 0;
            jj++;
            status = GSL_CONTINUE;
            if ( (gsl_root_fdfsolver_set(s0, &FDF, gsl_rng_uniform_pos(inputParams->rng)*Rhigh)) != 0 ) {
               fprintf(stderr,"%s: Unable to initialize root solver to first guess.\n", __func__);
               XLAL_ERROR_VOID(XLAL_EFUNC);
            }
         } else if (root<0.0 && jj==max_retries) {
            status = GSL_FAILURE;
         } //Up to here
         
      }
      
   } /* while status==GSL_CONTINUE && ii < max_iter */
   
   //Failure modes
   if (method != 0) {
      if (status != GSL_SUCCESS) {
         fprintf(stderr,"%s: Root finding iteration (%d/%d) failed with failure code %d. Previous root = %f, current root = %f\n", __func__, ii, max_iter, status, prevroot, root);
         XLAL_ERROR_VOID(XLAL_FAILURE);
      } else if (ii==max_iter) {
         fprintf(stderr,"%s: Root finding iteration (%d/%d) failed with failure code %d. Previous root = %f, current root = %f\n", __func__, ii, max_iter, status, prevroot, root);
         XLAL_ERROR_VOID(XLAL_EMAXITER);
      } else if (root == 0.0) {
         fprintf(stderr,"%s: Root finding iteration (%d/%d) converged to 0.0.\n", __func__, ii, max_iter);
         XLAL_ERROR_VOID(XLAL_ERANGE);
      } else if (root == 1000.0) {
         fprintf(stderr,"%s: Root finding iteration (%d/%d) converged to 1000.0.\n", __func__, ii, max_iter);
         XLAL_ERROR_VOID(XLAL_ERANGE);
      }
   } else {
      if (status != GSL_SUCCESS) {
         fprintf(stderr,"%s: Root finding iteration (%d/%d) failed with failure code %d. Previous root = %f, current root = %f\n", __func__, ii, max_iter, status, prevroot, root);
         XLAL_ERROR_VOID(XLAL_FAILURE);
      } else if (ii==max_iter) {
         fprintf(stderr,"%s: Root finding iteration (%d/%d) failed with failure code %d. Previous root = %f, current root = %f\n", __func__, ii, max_iter, status, prevroot, root);
         XLAL_ERROR_VOID(XLAL_EMAXITER);
      } else if (root<=0.0) {
         fprintf(stderr,"%s: Threshold value found (%f) is less than 0.0!\n", __func__, root);
         XLAL_ERROR_VOID(XLAL_ERANGE);
      }
   }
   
   
   output->far = root;
   output->distMean = 0.0;
   output->distSigma = 0.0; //Fake the value of sigma
   output->farerrcode = errcode;
   
   //Cleanup
   gsl_root_fsolver_free(s1);
   gsl_root_fdfsolver_free(s0);
   
   
} /* numericFAR() */

//For the root finding, calculating the false alarm probability of R
//Takes an average of 3 values of close by R values for stability
REAL8 gsl_probR(REAL8 R, void *param)
{
   
   struct gsl_probR_pars *pars = (struct gsl_probR_pars*)param;
   
   REAL8 dR = 0.005;
   REAL8 R1 = (1.0+dR)*R;
   REAL8 R2 = (1.0-dR)*R;
   INT4 errcode1 = 0, errcode2 = 0, errcode3 = 0;
   
   REAL8 prob = (probR(pars->templatestruct, pars->ffplanenoise, pars->fbinaveratios, R, pars->inputParams, &errcode1) + probR(pars->templatestruct, pars->ffplanenoise, pars->fbinaveratios, R1, pars->inputParams, &errcode2) + probR(pars->templatestruct, pars->ffplanenoise, pars->fbinaveratios, R2, pars->inputParams, &errcode3))/3.0;
   
   if (errcode1!=0) {
      pars->errcode = errcode1;
   } else if (errcode2!=0) {
      pars->errcode = errcode2;
   } else if (errcode3!=0) {
      pars->errcode = errcode3;
   }
   
   REAL8 returnval = prob - log10(pars->threshold);
   
   return returnval;
   
} /* gsl_probR() */

//When doing Newton's method, we need the slope
REAL8 gsl_dprobRdR(REAL8 R, void *param)
{
   
   struct gsl_probR_pars *pars = (struct gsl_probR_pars*)param;
   
   REAL8 dR = 0.005;
   
   INT4 errcode1 = 0, errcode2 = 0;
   
   //Explicit computation of slope
   REAL8 R1 = (1.0+dR)*R;
   REAL8 R2 = (1.0-dR)*R;
   REAL8 prob1 = gsl_probR(R1, pars);
   REAL8 prob2 = gsl_probR(R2, pars);
   while (fabs(prob1-prob2)<100.0*LAL_REAL8_EPS) {
      dR *= 2.0;
      R1 = (1.0+dR)*R;
      R2 = (1.0-dR)*R;
      prob1 = gsl_probR(R1, pars);
      prob2 = gsl_probR(R2, pars);
   }
   REAL8 diffR = R1 - R2;
   REAL8 slope = (prob1-prob2)/diffR;

   if (errcode1!=0)  pars->errcode = errcode1;
   else if (errcode2!=0) pars->errcode = errcode2;
   
   return slope;
   
} /* gsl_dprobRdR() */

//For Newton's method, we need the slope
void gsl_probRandDprobRdR(REAL8 R, void *param, REAL8 *probabilityR, REAL8 *dprobRdR)
{
   
   struct gsl_probR_pars *pars = (struct gsl_probR_pars*)param;
   
   *probabilityR = gsl_probR(R, pars);
   
   *dprobRdR = gsl_dprobRdR(R, pars);
   
} /* gsl_probRandDprobRdR() */


//////////////////////////////////////////////////////////////
// Analytically calculate the probability of a true signal using the Davies' method
// output is log10(prob)
REAL8 probR(templateStruct *templatestruct, REAL4Vector *ffplanenoise, REAL4Vector *fbinaveratios, REAL8 R, inputParamsStruct *params, INT4 *errcode)
{
   
   INT4 ii = 0;
   REAL8 prob = 0.0;
   REAL8 sumwsq = 0.0;
   INT4 numweights = 0;
   for (ii=0; ii<(INT4)templatestruct->templatedata->length; ii++) {
      if (templatestruct->templatedata->data[ii]!=0.0) {
         numweights++;
         sumwsq += templatestruct->templatedata->data[ii]*templatestruct->templatedata->data[ii];
      } else {
         break;
      }
   }
   
   REAL8Vector *newweights = XLALCreateREAL8Vector(numweights);
   INT4Vector *sorting = XLALCreateINT4Vector(numweights);
   if (newweights==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL8Vector(%d) failed.\n", __func__, numweights);
      XLAL_ERROR_REAL8(XLAL_EFUNC);
   } else if (sorting==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", __func__, numweights);
      XLAL_ERROR_REAL8(XLAL_EFUNC);
   }
   
   REAL8 Rpr = R;
   for (ii=0; ii<(INT4)newweights->length; ii++) {
      newweights->data[ii] = 0.5*templatestruct->templatedata->data[ii]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[ii] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[ii] ]/sumwsq;
      Rpr += templatestruct->templatedata->data[ii]*ffplanenoise->data[ templatestruct->secondfftfrequencies->data[ii] ]*fbinaveratios->data[ templatestruct->firstfftfrequenciesofpixels->data[ii] ]/sumwsq;
      sorting->data[ii] = ii;  //This is for the fact that a few steps later (before using Davies' algorithm, we sort the weights)
   }
   
   qfvars vars;
   vars.weights = newweights;
   vars.sorting = sorting;
   vars.dofs = NULL;
   vars.noncentrality = NULL;
   vars.ndtsrt = 0;           //Set because we do the sorting outside of Davies' algorithm with qsort
   vars.lim = 50000000;
   vars.c = Rpr;
   vars.useSSE = params->useSSE;
   REAL8 sigma = 0.0;
   REAL8 accuracy = 1.0e-11;   //(1e-5) old value
   
   //sort the weights here so we don't have to do it later (qsort)
   sort_double_ascend(newweights);
   
   //cdfwchisq(algorithm variables, sigma, accuracy, error code)
   prob = 1.0 - cdfwchisq_twospect(&vars, sigma, accuracy, errcode);
   
   //Large R values can cause a problem when computing the probability. We run out of accuracy quickly even using double precision
   //Potential fix: compute log10(prob) for smaller values of R, for when slope is linear between log10 probabilities
   //Use slope to extend the computation and then compute the exponential of the found log10 probability.
   REAL8 logprobest = 0.0;
   INT4 estimatedTheProb = 0;
   if (prob<=1.0e-8 || *errcode!=0) {
      estimatedTheProb = 1;
      
      INT4 errcode1 = 0;
      REAL8 probslope=0.0, tempprob, c1;
      REAL8 lowerend = 0.0;
      REAL8 upperend = Rpr;
      REAL8Vector *probvals = XLALCreateREAL8Vector(20);
      REAL8Vector *cvals = XLALCreateREAL8Vector(20);
      if (probvals==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL8Vector(%d) failed.\n", __func__, 20);
         XLAL_ERROR_REAL8(XLAL_EFUNC);
      } else if (cvals==NULL) {
         fprintf(stderr,"%s: XLALCreateREAL8Vector(%d) failed.\n", __func__, 20);
         XLAL_ERROR_REAL8(XLAL_EFUNC);
      }
      
      for (ii=0; ii<(INT4)probvals->length; ii++) {
         c1 = gsl_rng_uniform_pos(params->rng)*(upperend-lowerend)+lowerend;
         vars.c = c1;
         tempprob = 1.0-cdfwchisq_twospect(&vars, sigma, accuracy, &errcode1);
         while (tempprob<=1.0e-8 || tempprob>=1.0e-6) {
            if (tempprob<=1.0e-8) upperend = c1;
            else if (tempprob>=1.0e-6) lowerend = c1;
            c1 = gsl_rng_uniform_pos(params->rng)*(upperend-lowerend)+lowerend;
            vars.c = c1;
            tempprob = 1.0-cdfwchisq_twospect(&vars, sigma, accuracy, &errcode1);
            
            INT4 tries = 1;
            while (tries<10 && errcode1 != 0) {
               tries++;
               c1 = gsl_rng_uniform_pos(params->rng)*(upperend-lowerend)+lowerend;
               vars.c = c1;
               tempprob = 1.0-cdfwchisq_twospect(&vars, sigma, accuracy, &errcode1);
            }
            if (tries>=10 && errcode1!=0) {
               fprintf(stderr,"%s: cdfwchisq_twospect() failed with code %d after making %d tries.\n", __func__, errcode1, tries);
               XLAL_ERROR_REAL8(XLAL_EFUNC);
            }
            
         }
         if (errcode1!=0) {
            fprintf(stderr,"%s: cdfwchisq_twospect() failed with code %d.\n", __func__, errcode1);
            XLAL_ERROR_REAL8(XLAL_EFUNC);
         }
         probvals->data[ii] = log10(tempprob);
         cvals->data[ii] = c1;
      }
      
      REAL8 yintercept, cov00, cov01, cov11, sumsq;
      if (gsl_fit_linear(cvals->data, 1, probvals->data, 1, cvals->length, &yintercept, &probslope, &cov00, &cov01, &cov11, &sumsq)!=GSL_SUCCESS) {
         fprintf(stderr,"%s: gsl_fit_linear() failed.\n", __func__);
         XLAL_ERROR_REAL8(XLAL_EFUNC);
      }
      logprobest = probslope*Rpr + yintercept;
      if (logprobest>-0.5) {
         fprintf(stderr, "%s: Failure calculating accurate interpolated value.\n", __func__);
         XLAL_ERROR_REAL8(XLAL_ERANGE);
      }
      
      XLALDestroyREAL8Vector(probvals);
      XLALDestroyREAL8Vector(cvals);
      
      *errcode = errcode1;
      
   }
   
   //If errcode is still != 0, better fail
   if (*errcode!=0) {
      fprintf(stderr,"%s: cdfwchisq_twospect() failed at the end with code %d.\n", __func__, *errcode);
      XLAL_ERROR_REAL8(XLAL_EFUNC);
   }
   
   //Cleanup
   XLALDestroyREAL8Vector(newweights);
   XLALDestroyINT4Vector(sorting);
   
   if (estimatedTheProb==1) return logprobest;
   else return log10(prob);

} /* probR() */


//Create a new template structure
templateStruct * new_templateStruct(INT4 length)
{

   templateStruct *templatestruct = XLALMalloc(sizeof(*templatestruct));
   if (templatestruct==NULL) {
      fprintf(stderr,"%s: XLALMalloc(%zu) failed.\n", __func__, sizeof(*templatestruct));
      XLAL_ERROR_NULL(XLAL_ENOMEM);
   }
   
   templatestruct->templatedata = XLALCreateREAL4Vector(length);
   templatestruct->pixellocations = XLALCreateINT4Vector(length);
   templatestruct->firstfftfrequenciesofpixels = XLALCreateINT4Vector(length);
   templatestruct->secondfftfrequencies = XLALCreateINT4Vector(length);
   if (templatestruct->templatedata==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, length);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   } else if (templatestruct->pixellocations==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", __func__, length);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   } else if (templatestruct->firstfftfrequenciesofpixels==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", __func__, length);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   } else if (templatestruct->secondfftfrequencies==NULL) {
      fprintf(stderr,"%s: XLALCreateINT4Vector(%d) failed.\n", __func__, length);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
   
   memset(templatestruct->templatedata->data, 0, sizeof(REAL4)*length);
   memset(templatestruct->pixellocations->data, 0, sizeof(INT4)*length);
   memset(templatestruct->firstfftfrequenciesofpixels->data, 0, sizeof(INT4)*length);
   memset(templatestruct->secondfftfrequencies->data, 0, sizeof(INT4)*length);
   
   templatestruct->f0 = 0.0;
   templatestruct->period = 0.0;
   templatestruct->moddepth = 0.0;
   
   return templatestruct;
   
} /* new_templateStruct() */


//Reset the values in the template structure
void resetTemplateStruct(templateStruct *templatestruct)
{
   
   INT4 length = (INT4)templatestruct->templatedata->length;
   memset(templatestruct->templatedata->data, 0, sizeof(REAL4)*length);
   memset(templatestruct->pixellocations->data, 0, sizeof(INT4)*length);
   memset(templatestruct->firstfftfrequenciesofpixels->data, 0, sizeof(INT4)*length);
   memset(templatestruct->secondfftfrequencies->data, 0, sizeof(INT4)*length);
   
   templatestruct->f0 = 0.0;
   templatestruct->period = 0.0;
   templatestruct->moddepth = 0.0;
   
} /* resetTemplateStruct() */


//Free the memory of a template structure
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
// It is, essentially, eq. 18 of E. Goetz and K. Riles (2011)
// I've made it so that it handles spillage of power into neighboring bins a little more gracefully
// Numerical stability issues means that we need to compute exp(log(eq. 18)) = eq. 18
// exp(log(eq. 18)) = exp(log(4*pi*sigma^2)-sigma^2*omegapr^2)*(1+cos(delta*omegapr))*exp(log(1-cos(N*P*omegapr))-log(P*omegapr))
void makeTemplateGaussians(templateStruct *output, candidate input, inputParamsStruct *params, INT4 numfbins, INT4 numfprbins)
{
   
   //Set data for output template
   output->f0 = input.fsig;
   output->period = input.period;
   output->moddepth = input.moddepth;
   if (input.fsig==0.0 || input.period==0.0 || input.moddepth==0.0) {
      fprintf(stderr, "%s: Invalid input (%f, %f, %f).\n", __func__, input.fsig, input.period, input.moddepth);
      XLAL_ERROR_VOID(XLAL_EINVAL);
   }
   
   INT4 ii, jj, N;
   
   //Reset the data values to zero, just in case
   memset(output->templatedata->data, 0, sizeof(REAL4)*output->templatedata->length);
   
   N = (INT4)floor(params->Tobs/input.period);     //Number of Gaussians = observation time / period
   
   REAL8 periodf = 1.0/input.period;
   
   //Determine separation in time of peaks for each frequency
   //phi = P/2 - P/pi * asin[(f-f0)/modulation depth]
   //When abs(f-f0)>modulation depth, phi := 0
   REAL4Vector *phi_actual = XLALCreateREAL4Vector(numfbins);
   if (phi_actual==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, numfbins);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   for (ii=0; ii<(INT4)phi_actual->length; ii++) {
      if ( fabs(params->fmin*params->Tcoh - params->dfmax*params->Tcoh - 6.0 + ii - input.fsig*params->Tcoh)/(input.moddepth*params->Tcoh)-1.0e-14 <= 1.0 ) phi_actual->data[ii] = 0.5*input.period - asin(fabs(params->fmin*params->Tcoh - params->dfmax*params->Tcoh - 6.0 + ii - input.fsig*params->Tcoh)/(input.moddepth*params->Tcoh)-1.0e-14)*LAL_1_PI*input.period;
      else phi_actual->data[ii] = 0.0;
   } /* for ii < phi_actual->length */
   
   //Create second FFT frequencies and other useful values
   REAL4Vector *fpr = XLALCreateREAL4Vector(numfprbins);
   if (fpr==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, numfprbins);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   for (ii=0; ii<(INT4)fpr->length; ii++) fpr->data[ii] = (REAL4)ii*(1.0/params->Tobs);
   
   //For speed, we will precompute a number of useful vectors described by their names
   //This part is the allocation
   REAL4Vector *omegapr = XLALCreateREAL4Vector(fpr->length);
   REAL4Vector *omegapr_squared = XLALCreateREAL4Vector(fpr->length);
   REAL4Vector *cos_ratio = XLALCreateREAL4Vector(fpr->length);
   if (omegapr==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, fpr->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (omegapr_squared==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, fpr->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (cos_ratio==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, fpr->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   
   //Doing the precomputation of the useful values
   if (params->useSSE) {
      sseScaleREAL4Vector(omegapr, fpr, (REAL4)LAL_TWOPI);
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: sseScaleREAL4Vector() failed.\n", __func__);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      }
      
      sseSSVectorMultiply(omegapr_squared, omegapr, omegapr);
      if (xlalErrno!=0) {
         fprintf(stderr,"%s: sseSSVectorMultiply() failed.\n", __func__);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      }
      
      for (ii=0; ii<(INT4)fpr->length; ii++) {
         REAL4 cos_omegapr_times_period = cosf((REAL4)(input.period*omegapr->data[ii]));
         REAL4 cos_N_times_omegapr_times_period = cosf((REAL4)(N*input.period*omegapr->data[ii]));
         if (cos_omegapr_times_period<=(1.0-100.0*LAL_REAL4_EPS)) cos_ratio->data[ii] = (1.0 - cos_N_times_omegapr_times_period)/(1.0 - cos_omegapr_times_period);
         else cos_ratio->data[ii] = (REAL4)(N*N);
      }
   } else {
      for (ii=0; ii<(INT4)fpr->length; ii++) {
         omegapr->data[ii] = (REAL4)LAL_TWOPI*fpr->data[ii];
         omegapr_squared->data[ii] = omegapr->data[ii]*omegapr->data[ii];
         REAL4 cos_omegapr_times_period = cosf((REAL4)(input.period*omegapr->data[ii]));
         REAL4 cos_N_times_omegapr_times_period = cosf((REAL4)(N*input.period*omegapr->data[ii]));
         if (cos_omegapr_times_period<=(1.0-100.0*LAL_REAL4_EPS)) cos_ratio->data[ii] = (1.0 - cos_N_times_omegapr_times_period)/(1.0 - cos_omegapr_times_period);
         else cos_ratio->data[ii] = (REAL4)(N*N);
      }
   }
   
   //Scale used for "spillover" into bins outside of phi_actual
   //REAL4 k = input.moddepth*params->Tcoh;    //amplitude of modulation in units of bins
   REAL4Vector *scale = XLALCreateREAL4Vector(numfbins);      //the scaling factor
   if (scale==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, numfbins);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   INT4 bin0 = (INT4)round(params->fmin*params->Tcoh - params->dfmax*params->Tcoh - 6.0);      //bin number of fmin
   REAL4 m0 = input.fsig*params->Tcoh - bin0;  //central frequency bin
   REAL4 mextent = input.moddepth*params->Tcoh;  //Bins filled by modulation                
   REAL4 overage = 1.0;  //spillage
   INT4 fnumstart = (INT4)round(m0-mextent-overage), fnumend = (INT4)round(m0+mextent+overage);  //start and end bins
   memset(scale->data, 0, numfbins*sizeof(*scale->data));
   for (ii=fnumstart; ii<=fnumend; ii++) {
      if ((REAL4)ii>=(m0-mextent) && (REAL4)ii<=(m0+mextent)) scale->data[ii] = 1.0;
      else scale->data[ii] = sqsincxoverxsqminusone(fmin(fabs((REAL4)ii-(m0-mextent)), fabs((REAL4)ii-(m0+mextent))));
   }
   
   //Make sigmas for each frequency
   //First, allocate vectors
   REAL4Vector *sigmas = XLALCreateREAL4Vector((UINT4)(fnumend-fnumstart+1));
   REAL4Vector *wvals = XLALCreateREAL4Vector((UINT4)floor(4.0*input.period/params->Tcoh));
   if (sigmas==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, (UINT4)(fnumend-fnumstart+1));
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (wvals==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, (UINT4)floor(4.0*input.period/params->Tcoh));
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   REAL4Vector *allsigmas = XLALCreateREAL4Vector(wvals->length * sigmas->length);
   REAL4Vector *weightvals = XLALCreateREAL4Vector(wvals->length * sigmas->length);
   if (allsigmas==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, wvals->length * sigmas->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (weightvals==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, wvals->length * sigmas->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   
   //Here is where the sigmas are computed. It is a weighted average. t = (ii+1)*in->Tcoh*0.5
   REAL8 sin2pix = 0.0, cos2pix = 0.0;
   for (ii=0; ii<(INT4)wvals->length; ii++) {
      //calculate sin and cos of 2*pi*t/P and then the bin the signal is in and the signal velocity
      twospect_sin_cos_2PI_LUT(&sin2pix, &cos2pix, periodf*((ii+1)*params->Tcoh*0.5));
      REAL4 sigbin = (input.moddepth*cos2pix+input.fsig)*params->Tcoh;
      REAL4 sigbinvelocity = fabs(-input.moddepth*sin2pix*params->Tcoh*params->Tcoh*LAL_PI*periodf);
      
      //if the velocity approaches zero, the sigma calculation will diverge (which it should do) but this is bad numerically, so we cap it
      if (sigbinvelocity<1.0e-4) sigbinvelocity = 1.0e-4;

      REAL4 sigma = 0.5*params->Tcoh * (0.5346 * powf(sigbinvelocity, -1.0213f));   //Derived fit from simulation
      //REAL4 sigma = 0.5*params->Tcoh * (0.5979 / (sigbinvelocity - 3.2895e-5));  //Could think about using this fit in the future
      
      for (jj=0; jj<(INT4)sigmas->length; jj++) {
         weightvals->data[ii*sigmas->length + jj] = sqsincxoverxsqminusone(sigbin-(bin0+jj+fnumstart));
         allsigmas->data[ii*sigmas->length + jj] = weightvals->data[ii*sigmas->length + jj]*sigma;
      }
      
   } /* for ii < wvals->length */
   for (ii=0; ii<(INT4)sigmas->length; ii++) {
      REAL8 wavesigma = 0.0;
      REAL8 totalw = 0.0;
      for (jj=0; jj<(INT4)wvals->length; jj++) {
         wavesigma += allsigmas->data[ii + jj*sigmas->length];
         totalw += weightvals->data[ii + jj*sigmas->length];
      }
      sigmas->data[ii] = (REAL4)(wavesigma/totalw);
   } /* for ii < sigmas->length */

   //Allocate more useful data vectors. These get computed for each different first FFT frequency bin in the F-F plane
   REAL4Vector *exp_neg_sigma_sq_times_omega_pr_sq = XLALCreateREAL4Vector(omegapr_squared->length);
   if (exp_neg_sigma_sq_times_omega_pr_sq==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, omegapr_squared->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   REAL4Vector *sin_phi_times_omega_pr = XLALCreateREAL4Vector(omegapr->length);
   if (sin_phi_times_omega_pr==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, omegapr->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   REAL4Vector *cos_phi_times_omega_pr = XLALCreateREAL4Vector(omegapr->length);
   if (cos_phi_times_omega_pr==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, omegapr->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   REAL4Vector *phi_times_fpr = XLALCreateREAL4Vector(fpr->length);
   if (phi_times_fpr==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, fpr->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   REAL4Vector *datavector = XLALCreateREAL4Vector(fpr->length);
   if (datavector==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, fpr->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   
   //Create template. We are going to do exp(log(Eq. 18))
   REAL8 sum = 0.0;
   REAL4 log4pi = 2.53102424697f;
   for (ii=0; ii<(INT4)sigmas->length; ii++) {
      REAL4 s = sigmas->data[ii];      //sigma
      
      //Scaling factor for leakage
      REAL4 scale1 = sqrtf((REAL4)(1.0/(1.0+expf((REAL4)(-phi_actual->data[ii+fnumstart]*phi_actual->data[ii+fnumstart]*0.5/(s*s))))));
      
      //pre-factor
      //REAL8 prefact0 = scale1 * 2.0 * LAL_TWOPI * s * s;
      //REAL4 prefact0 = log(scale1 * 2.0 * LAL_TWOPI * s * s);     //We are going to do exp(log(Eq. 18))
      REAL4 prefact0 = log4pi + 2.0*logf((REAL4)(scale1*s));

      if (params->useSSE) {
         //Compute exp(log(4*pi*s*s*exp(-s*s*omegapr_squared))) = exp(log(4*pi*s*s)-s*s*omegapr_squared)
         sseScaleREAL4Vector(exp_neg_sigma_sq_times_omega_pr_sq, omegapr_squared, -s*s);
         if (xlalErrno!=0) {
            fprintf(stderr, "%s: sseScaleREAL4Vector() failed.\n", __func__);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
         sseAddScalarToREAL4Vector(exp_neg_sigma_sq_times_omega_pr_sq, exp_neg_sigma_sq_times_omega_pr_sq, prefact0);
         if (xlalErrno!=0) {
            fprintf(stderr, "%s: sseAddScalarToREAL4Vector() failed.\n", __func__);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
         sse_exp_REAL4Vector(exp_neg_sigma_sq_times_omega_pr_sq, exp_neg_sigma_sq_times_omega_pr_sq);
         if (xlalErrno!=0) {
            fprintf(stderr, "%s: sse_exp_REAL4Vector() failed.\n", __func__);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
         
         //Compute phi_actual*fpr
         sseScaleREAL4Vector(phi_times_fpr, fpr, phi_actual->data[ii+fnumstart]);
         if (xlalErrno!=0) {
            fprintf(stderr, "%s: sseScaleREAL4Vector() failed.\n", __func__);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
         
         //Start computing the datavector values
         INT4 maxindex = max_index(phi_times_fpr);
         if (phi_times_fpr->data[maxindex]<=2.147483647e9) {
            //Compute cos(2*pi*phi_actual*fpr) using LUT and SSE
            sse_sin_cos_2PI_LUT_REAL4Vector(sin_phi_times_omega_pr, cos_phi_times_omega_pr, phi_times_fpr);
            if (xlalErrno!=0) {
               fprintf(stderr, "%s: sse_sin_cos_2PI_LUT_REAL4Vector() failed.\n", __func__);
               XLAL_ERROR_VOID(XLAL_EFUNC);
            }
         } else {
            //Compute cos(2*pi*phi_actual*fpr) using LUT
            for (jj=0; jj<(INT4)omegapr_squared->length; jj++) {
               twospect_sin_cos_2PI_LUT(&sin2pix, &cos2pix, phi_times_fpr->data[jj]);
               cos_phi_times_omega_pr->data[jj] = (REAL4)cos2pix;
            }
         }
         //datavector = cos(phi_actual*omega_pr) + 1.0
         sseAddScalarToREAL4Vector(datavector, cos_phi_times_omega_pr, 1.0);
         if (xlalErrno!=0) {
            fprintf(stderr, "%s: sseAddScalarToREAL4Vector() failed.\n", __func__);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
         //datavector = prefact0 * exp(-s*s*omega_pr*omega_pr) * [cos(phi_actual*omega_pr) + 1.0]
         sseSSVectorMultiply(datavector, datavector, exp_neg_sigma_sq_times_omega_pr_sq);
         if (xlalErrno!=0) {
            fprintf(stderr, "%s: sseSSVectorMultiply() failed.\n", __func__);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
         //datavector = scale * exp(-s*s*omega_pr*omega_pr) * [cos(phi_actual*omega_pr) + 1.0]
         sseScaleREAL4Vector(datavector, datavector, scale->data[ii+fnumstart]);
         if (xlalErrno!=0) {
            fprintf(stderr, "%s: sseScaleREAL4Vector() failed.\n", __func__);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
         //datavector *= cos_ratio
         sseSSVectorMultiply(datavector, datavector, cos_ratio);
         if (xlalErrno!=0) {
            fprintf(stderr, "%s: sseSSVectorMultiply() failed.\n", __func__);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
      } else {
         for (jj=0; jj<(INT4)omegapr_squared->length; jj++) {
            //Do all or nothing if the exponential is too negative
            if ((prefact0-s*s*omegapr_squared->data[jj])>-88.0) {
               exp_neg_sigma_sq_times_omega_pr_sq->data[jj] = expf((REAL4)(prefact0-s*s*omegapr_squared->data[jj]));
               twospect_sin_cos_2PI_LUT(&sin2pix, &cos2pix, phi_actual->data[ii+fnumstart]*fpr->data[jj]);
               cos_phi_times_omega_pr->data[jj] = (REAL4)cos2pix;
               datavector->data[jj] = scale->data[ii+fnumstart]*exp_neg_sigma_sq_times_omega_pr_sq->data[jj]*(cos_phi_times_omega_pr->data[jj]+1.0)*cos_ratio->data[jj];
            } else {
               datavector->data[jj] = 0.0;
            }
         } /* for jj = 0 --> omegapr_squared->length */
      } /* use SSE or not */
      
      //Now loop through the second FFT frequencies, starting with index 4
      for (jj=4; jj<(INT4)omegapr->length; jj++) {
         //Sum up the weights in total
         sum += (REAL8)(datavector->data[jj]);
         
         //Compare with weakest top bins and if larger, launch a search to find insertion spot (insertion sort)
         if (datavector->data[jj] > output->templatedata->data[output->templatedata->length-1]) {
            insertionSort_template(output, datavector->data[jj], (ii+fnumstart)*fpr->length+jj, ii+fnumstart, jj);
         }
      } /* for jj < omegapr->length */
   } /* for ii < sigmas->length */
   
   //Normalize
   REAL4 invsum = (REAL4)(1.0/sum);
   if (!params->useSSE) {
      for (ii=0; ii<(INT4)output->templatedata->length; ii++) if (output->templatedata->data[ii]!=0.0) output->templatedata->data[ii] *= invsum;
   } else {
      output->templatedata = sseScaleREAL4Vector(output->templatedata, output->templatedata, invsum);
      if (xlalErrno!=0) {
         fprintf(stderr,"%s, sseScaleREAL4Vector() failed.\n", __func__);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      }
   }
   
   //Truncate weights when they don't add much to the total sum of weights
   sum = 0.0;
   for (ii=0; ii<params->mintemplatelength; ii++) sum += (REAL8)output->templatedata->data[ii];
   ii = params->mintemplatelength;
   while (ii<(INT4)output->templatedata->length && output->templatedata->data[ii]>=epsval_float((REAL4)sum)) {
      sum += (REAL8)output->templatedata->data[ii];
      ii++;
   }
   for (/* last ii val */; ii<(INT4)output->templatedata->length; ii++) output->templatedata->data[ii] = 0.0;
   
   //Destroy variables
   XLALDestroyREAL4Vector(phi_actual);
   XLALDestroyREAL4Vector(scale);
   XLALDestroyREAL4Vector(sigmas);
   XLALDestroyREAL4Vector(allsigmas);
   XLALDestroyREAL4Vector(weightvals);
   XLALDestroyREAL4Vector(wvals);
   XLALDestroyREAL4Vector(fpr);
   XLALDestroyREAL4Vector(omegapr);
   XLALDestroyREAL4Vector(omegapr_squared);
   XLALDestroyREAL4Vector(cos_ratio);
   XLALDestroyREAL4Vector(exp_neg_sigma_sq_times_omega_pr_sq);
   XLALDestroyREAL4Vector(phi_times_fpr);
   XLALDestroyREAL4Vector(sin_phi_times_omega_pr);
   XLALDestroyREAL4Vector(cos_phi_times_omega_pr);
   XLALDestroyREAL4Vector(datavector);

} /* mateTemplateGaussians() */



//////////////////////////////////////////////////////////////
// Make an template based on FFT of sinc squared functions
// This is eq. 20 of E. Goetz and K. Riles (2011)
void makeTemplate(templateStruct *output, candidate input, inputParamsStruct *params, INT4Vector *sftexist, REAL4FFTPlan *plan)
{
   
   //Set data for output template
   output->f0 = input.fsig;
   output->period = input.period;
   output->moddepth = input.moddepth;
   
   INT4 ii, jj, numfbins, numffts;
   
   //Reset to zero, just in case
   memset(output->templatedata->data, 0, sizeof(REAL4)*output->templatedata->length);
   
   numfbins = (INT4)(round(params->fspan*params->Tcoh+2.0*params->dfmax*params->Tcoh)+12+1);   //Number of frequency bins
   numffts = (INT4)sftexist->length;   //Number of FFTs
   
   REAL4Vector *psd1 = XLALCreateREAL4Vector(numfbins*numffts);
   REAL8Vector *freqbins = XLALCreateREAL8Vector(numfbins);
   REAL8Vector *bindiffs = XLALCreateREAL8Vector(numfbins);
   if (psd1==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, numfbins*numffts);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (freqbins==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL8Vector(%d) failed.\n", __func__, numfbins);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (bindiffs==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL8Vector(%d) failed.\n", __func__, numfbins);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   memset(psd1->data, 0, sizeof(REAL4)*psd1->length);
   
   REAL8 periodf = 1.0/input.period;
   REAL8 B = input.moddepth*params->Tcoh;
   
   //Bin numbers of the frequencies
   for (ii=0; ii<numfbins; ii++) freqbins->data[ii] = round(params->fmin*params->Tcoh - params->dfmax*params->Tcoh - 6.0) + ii;
   
   //Determine the signal modulation in bins with time at center of coherence time and create
   //Hann windowed PSDs
   REAL8 sin2pix = 0.0, cos2pix = 0.0;
   REAL8 PSDprefact = 2.0/3.0;
   for (ii=0; ii<numffts; ii++) {
      REAL8 t = 0.5*params->Tcoh*ii;  //Assumed 50% overlapping SFTs
      twospect_sin_cos_2PI_LUT(&sin2pix, &cos2pix, periodf*t);
      REAL8 n0 = B*sin2pix + input.fsig*params->Tcoh;
      if (params->useSSE) {
         sseAddScalarToREAL8Vector(bindiffs, freqbins, -n0);
         if (xlalErrno!=0) {
            fprintf(stderr,"%s: sseAddScalarToREAL8Vector() failed.\n", __func__);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
      } else for (jj=0; jj<numfbins; jj++) bindiffs->data[jj] = freqbins->data[jj] - n0;
      for (jj=0; jj<numfbins; jj++) {
         //Create PSD values organized by f0 => psd1->data[0...numffts-1], sft1 => psd1->data[numffts...2*numffts-1]
         //Restricting to +/- 1.75 bins means >99.9% of the total power is included in the template calculation
         if ( fabs(bindiffs->data[jj]) <= 1.75 ) psd1->data[ii + jj*numffts] = sqsincxoverxsqminusone(bindiffs->data[jj])*PSDprefact;
      } /* for jj < numfbins */
   } /* for ii < numffts */
   
   //Do the second FFT
   REAL4Vector *x = XLALCreateREAL4Vector(numffts);
   if (x==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, numffts);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   REAL4Window *win = XLALCreateHannREAL4Window(x->length);
   REAL4Vector *psd = XLALCreateREAL4Vector((INT4)floor(x->length*0.5)+1);
   if (win==NULL) {
      fprintf(stderr,"%s: XLALCreateHannREAL4Window(%d) failed.\n", __func__, x->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   } else if (psd==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, (INT4)floor(x->length*0.5)+1);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   REAL8 winFactor = 8.0/3.0;
   REAL8 secPSDfactor = winFactor/x->length*0.5*params->Tcoh;
   REAL8 sum = 0.0;
   INT4 doSecondFFT;
   //First loop over frequencies
   for (ii=0; ii<numfbins; ii++) {
      //Set doSecondFFT check flag to 0. Value becomes 1 if we are to do the second FFT
      doSecondFFT = 0;
   
      //Next, loop over times and check to see if we need to do second FFT
      //Sum up the power in the row and see if it exceeds 5.0*(sinc(3.0)/(3.0^2-1))^2
      REAL4 rowpowersum = 0.0;
      for (jj=0; jj<(INT4)x->length; jj++) rowpowersum += psd1->data[ii*numffts+jj];
      if (rowpowersum > 1.187167e-34) doSecondFFT = 1;
      
      //If we are to do the second FFT then do it!
      if (doSecondFFT) {
         //Obtain and window the time series
         memcpy(x->data, &(psd1->data[ii*numffts]), sizeof(REAL4)*x->length);
         if (!params->useSSE) {
            //x = fastSSVectorMultiply_with_stride_and_offset(x, x, win->data, 1, 1, 0, 0);
            x = XLALSSVectorMultiply(x, x, win->data);
            if (xlalErrno!=0) {
               fprintf(stderr,"%s, XLALSSVectorMultiply() failed.\n", __func__);
               XLAL_ERROR_VOID(XLAL_EFUNC);
            }
         } else {
            x = sseSSVectorMultiply(x, x, win->data);
            if (xlalErrno!=0) {
               fprintf(stderr,"%s, sseSSVectorMultiply() failed.\n", __func__);
               XLAL_ERROR_VOID(XLAL_EFUNC);
            }
         }
         
         //Do the FFT
         if ( XLALREAL4PowerSpectrum(psd, x, plan) != 0 ) {
            fprintf(stderr,"%s: XLALREAL4PowerSpectrum() failed.\n", __func__);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
         
         //Scale the data points by 1/N and window factor and (1/fs)
         //Order of vector is by second frequency then first frequency
         if (params->useSSE) {
            psd = sseScaleREAL4Vector(psd, psd, secPSDfactor);
            if (xlalErrno!=0) {
               fprintf(stderr,"%s, sseScaleREAL4Vector() failed.\n", __func__);
               XLAL_ERROR_VOID(XLAL_EFUNC);
            }
         } else {
            for (jj=0; jj<(INT4)psd->length; jj++) psd->data[jj] *= secPSDfactor;
         }
         
         //Ignore the DC to 3rd frequency bins in sum
         for (jj=4; jj<(INT4)psd->length; jj++) {
            sum += (REAL8)psd->data[jj];     //sum up the total weight
            
            //Sort the weights, insertion sort technique
            //if (correctedValue > output->templatedata->data[output->templatedata->length-1]) insertionSort_template(output, correctedValue, ii*psd->length+jj, ii, jj);
            if (psd->data[jj] > output->templatedata->data[output->templatedata->length-1]) insertionSort_template(output, psd->data[jj], ii*psd->length+jj, ii, jj);
         } /* for jj < psd->length */
      } /* if doSecondFFT */
   } /* if ii < numfbins */
   
   //Normalize
   REAL4 invsum = (REAL4)(1.0/sum);
   if (!params->useSSE) {
      for (ii=0; ii<(INT4)output->templatedata->length; ii++) if (output->templatedata->data[ii]!=0.0) output->templatedata->data[ii] *= invsum;
   } else {
      output->templatedata = sseScaleREAL4Vector(output->templatedata, output->templatedata, invsum);
      if (xlalErrno!=0) {
         fprintf(stderr,"%s, sseScaleREAL4Vector() failed.\n", __func__);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      }
   }
   
   //Truncate weights if they don't contribute much to the sum
   sum = 0.0;
   for (ii=0; ii<params->mintemplatelength; ii++) sum += (REAL8)output->templatedata->data[ii];
   ii = params->mintemplatelength;
   while (ii<(INT4)output->templatedata->length && output->templatedata->data[ii]>=epsval_float((REAL4)sum)) {
      sum += (REAL8)output->templatedata->data[ii];
      ii++;
   }
   for (/* last ii val */; ii<(INT4)output->templatedata->length; ii++) output->templatedata->data[ii] = 0.0;
   
   //Destroy stuff
   XLALDestroyREAL4Vector(psd1);
   XLALDestroyREAL8Vector(freqbins);
   XLALDestroyREAL8Vector(bindiffs);
   XLALDestroyREAL4Vector(x);
   XLALDestroyREAL4Window(win);
   XLALDestroyREAL4Vector(psd);
   
}


void analyzeOneTemplate(candidate *output, candidate *input, ffdataStruct *ffdata, REAL4Vector *aveNoise, REAL4Vector *aveTFnoisePerFbinRatio, inputParamsStruct *params, INT4Vector *sftexist, REAL4FFTPlan *plan)
{

   INT4 proberrcode = 0;

   //Allocate and make the template
   templateStruct *template = new_templateStruct(params->maxtemplatelength);
   if (template==NULL) {
      fprintf(stderr,"%s: new_templateStruct(%d) failed.\n", __func__, params->maxtemplatelength);
      XLAL_ERROR_VOID(XLAL_EFUNC); 
   }
   resetTemplateStruct(template);
   makeTemplate(template, *input, params, sftexist, plan);
   if (xlalErrno!=0) {
      fprintf(stderr,"%s: makeTemplate() failed.\n", __func__);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }

   //Print out data product if requested
   /* if (args_info.printData_given) {
      char w[1000];
      snprintf(w, 1000, "%s/%s", args_info.outdirectory_arg, "templatedata.dat");
      FILE *TEMPLATEDATA = fopen(w, "w");
      if (TEMPLATEDATA==NULL) {
        fprintf(stderr, "%s: fopen %s failed.\n", __func__, w);
        XLAL_ERROR(XLAL_EFUNC);
     }
     for (jj=0; jj<(INT4)template->templatedata->length; jj++) fprintf(TEMPLATEDATA, "%g %d %d %d\n", template->templatedata->data[jj], template->pixellocations->data[jj], template->firstfftfrequenciesofpixels->data[jj], template->secondfftfrequencies->data[jj]);
     fclose(TEMPLATEDATA);
     } */

   //Calculate R from the template and the data
   REAL8 R = calculateR(ffdata->ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
   if (XLAL_IS_REAL8_FAIL_NAN(R)) {
      fprintf(stderr,"%s: calculateR() failed.\n", __func__);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }

   //Calculate FAP
   REAL8 prob = probR(template, aveNoise, aveTFnoisePerFbinRatio, R, params, &proberrcode);
   if (XLAL_IS_REAL8_FAIL_NAN(prob)) {
      fprintf(stderr,"%s: probR() failed.\n", __func__);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }

   //Estimate the h0 if R>0.0
   REAL8 h0 = 0.0;
   if ( R > 0.0 ) h0 = 2.7426*pow(R/(params->Tcoh*params->Tobs),0.25);

   loadCandidateData(output, input->fsig, input->period, input->moddepth, input->ra, input->dec, R, h0, prob, proberrcode, 1.0);

}


//A brute force template search to find the most significant template around a candidate
void bruteForceTemplateSearch(candidate *output, candidate input, REAL8 fminimum, REAL8 fmaximum, INT4 numfsteps, INT4 numperiodslonger, INT4 numperiodsshorter, REAL8 dfmin, REAL8 dfmax, INT4 numdfsteps, inputParamsStruct *params, REAL4Vector *ffdata, INT4Vector *sftexist, REAL4Vector *aveNoise, REAL4Vector *aveTFnoisePerFbinRatio, REAL4FFTPlan *secondFFTplan, INT4 useExactTemplates)
{
   
   INT4 ii, jj, kk;
   REAL8Vector *trialf, *trialb, *trialp;
   REAL8 fstepsize, dfstepsize;
   REAL4 tcohfactor = 1.49e-3*params->Tcoh + 1.76;    //From in-text equation after Eq. 23 of E.G. and K.R. 2011
   REAL8 log10templatefar = params->log10templatefar;
   
   //Set up parameters of modulation depth search
   if (dfmin<(0.5/params->Tcoh-1.0e-9)) dfmin = 0.5/params->Tcoh;
   trialb = XLALCreateREAL8Vector(numdfsteps);
   if (trialb==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL8Vector(%d) failed.\n", __func__, numdfsteps);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   dfstepsize = (dfmax-dfmin)/(REAL8)(numdfsteps-1);
   for (ii=0; ii<numdfsteps; ii++) trialb->data[ii] = dfmin + dfstepsize*ii;
   
   //Set up parameters of signal frequency search
   if (fminimum<params->fmin) fminimum = params->fmin;
   if (fmaximum>params->fmin+params->fspan) fmaximum = params->fmin+params->fspan;
   trialf = XLALCreateREAL8Vector(numfsteps);
   if (trialf==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL8Vector(%d) failed.\n", __func__, numfsteps);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   fstepsize = (fmaximum-fminimum)/(REAL8)(numfsteps-1);
   for (ii=0; ii<numfsteps; ii++) trialf->data[ii] = fminimum + fstepsize*ii;
   
   //Search over numperiods different periods
   trialp = XLALCreateREAL8Vector(numperiodslonger+numperiodsshorter+1);
   if (trialp==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL8Vector(%d) failed.\n", __func__, numperiodslonger+numperiodsshorter+1);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   
   //Now search over the parameter space. Frequency, then modulation depth, then period
   //Initialze best values as the initial point we are searching around
   INT4 bestproberrcode = 0;
   REAL8 bestf = 0.0, bestp = 0.0, bestdf = 0.0, bestR = 0.0, besth0 = 0.0, bestProb = 0.0;
   candidate cand;
   templateStruct *template = new_templateStruct(params->maxtemplatelength);
   if (template==NULL) {
      fprintf(stderr,"%s: new_templateStruct(%d) failed.\n", __func__, params->maxtemplatelength);
      XLAL_ERROR_VOID(XLAL_EFUNC); 
   }
   farStruct *farval = NULL;
   if (params->calcRthreshold) {
      farval = new_farStruct();
      if (farval==NULL) {
         fprintf(stderr,"%s: new_farStruct() failed.\n", __func__);
         XLAL_ERROR_VOID(XLAL_EFUNC); 
      }
   }
   
   INT4 startposition = numperiodsshorter, proberrcode = 0;
   //Search over frequency
   for (ii=0; ii<(INT4)trialf->length; ii++) {
      //Search over modulation depth
      for (jj=0; jj<(INT4)trialb->length; jj++) {
         //Start with period of the first guess, then determine nearest neighbor from the
         //modulation depth amplitude to find the other period guesses. These parameters 
         //are determined from simulation to scale the N.N. distance w.r.t. mod. depth with
         //20% mismatch parameter
         trialp->data[startposition] = input.period;
         for (kk=0; kk<numperiodsshorter; kk++) {
            REAL8 nnp = trialp->data[startposition-kk]*trialp->data[startposition-kk]*(1+trialp->data[startposition-kk]/tcohfactor/params->Tobs)/tcohfactor/params->Tobs*sqrt(3.6e-3/trialb->data[jj]);
            trialp->data[startposition-(kk+1)] = trialp->data[startposition-kk] - nnp;
         }
         for (kk=0; kk<numperiodslonger; kk++) {
            REAL8 nnp = trialp->data[startposition+kk]*trialp->data[startposition+kk]*(1+trialp->data[startposition+kk]/tcohfactor/params->Tobs)/tcohfactor/params->Tobs*sqrt(3.6e-3/trialb->data[jj]);
            trialp->data[startposition+(kk+1)] = trialp->data[startposition+kk] + nnp;
         }
         
         //Search over period
         for (kk=0; kk<(INT4)trialp->length; kk++) {
            //Within boundaries?
            if ( trialf->data[ii]>=params->fmin && 
                trialf->data[ii]<(params->fmin+params->fspan) && 
                trialb->data[jj]<maxModDepth(trialp->data[kk], params->Tcoh) && 
                trialp->data[kk]>minPeriod(trialb->data[jj], params->Tcoh) && 
                trialp->data[kk]<=(0.2*params->Tobs) && 
                trialp->data[kk]>=(2.0*3600.0) && 
                trialb->data[jj]>=params->dfmin && 
                trialb->data[jj]<=params->dfmax && 
                trialp->data[kk]<=params->Pmax && 
                trialp->data[kk]>=params->Pmin ) {
               
               loadCandidateData(&cand, trialf->data[ii], trialp->data[kk], trialb->data[jj], input.ra, input.dec, 0, 0, 0.0, 0, 0.0);
               
               resetTemplateStruct(template);
               
               if (useExactTemplates!=0) {
                  makeTemplate(template, cand, params, sftexist, secondFFTplan);
                  if (xlalErrno!=0) {
                     fprintf(stderr,"%s: makeTemplate() failed.\n", __func__);
                     XLAL_ERROR_VOID(XLAL_EFUNC);
                  }
               } else {
                  makeTemplateGaussians(template, cand, params, (INT4)aveTFnoisePerFbinRatio->length, (INT4)aveNoise->length);
                  if (xlalErrno!=0) {
                     fprintf(stderr,"%s: makeTemplateGaussians() failed.\n", __func__);
                     XLAL_ERROR_VOID(XLAL_EFUNC);
                  }
               }
               
               if (params->calcRthreshold && bestProb==0.0) {
                  numericFAR(farval, template, params->templatefar, aveNoise, aveTFnoisePerFbinRatio, params, params->rootFindingMethod);
                  if (xlalErrno!=0) {
                     fprintf(stderr,"%s: numericFAR() failed.\n", __func__);
                     XLAL_ERROR_VOID(XLAL_EFUNC);
                  }
               }
               
               REAL8 R = calculateR(ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
               if (XLAL_IS_REAL8_FAIL_NAN(R)) {
                  fprintf(stderr,"%s: calculateR() failed.\n", __func__);
                  XLAL_ERROR_VOID(XLAL_EFUNC);
               }
               REAL8 prob = probR(template, aveNoise, aveTFnoisePerFbinRatio, R, params, &proberrcode);
               if (XLAL_IS_REAL8_FAIL_NAN(prob)) {
                  fprintf(stderr,"%s: probR() failed.\n", __func__);
                  XLAL_ERROR_VOID(XLAL_EFUNC);
               }
               REAL8 h0 = 0.0;
               if ( R > 0.0 ) h0 = 2.7426*pow(R/(params->Tcoh*params->Tobs),0.25);
               
               if ( (bestProb!=0.0 && prob < bestProb) || (bestProb==0.0 && !params->calcRthreshold && prob<log10templatefar) || (bestProb==0.0 && params->calcRthreshold && R > farval->far) ) {
                  bestf = trialf->data[ii];
                  bestp = trialp->data[kk];
                  bestdf = trialb->data[jj];
                  bestR = R;
                  besth0 = h0;
                  bestProb = prob;
                  bestproberrcode = proberrcode;
               }
               
            } /* if within boundaries */
         } /* for kk < trialp */
      } /* for jj < trialb */
   } /* for ii < trialf */
   free_templateStruct(template);
   template = NULL;
   if (params->calcRthreshold) {
      free_farStruct(farval);
      farval = NULL;
   }
   XLALDestroyREAL8Vector(trialf);
   XLALDestroyREAL8Vector(trialb);
   XLALDestroyREAL8Vector(trialp);
   trialf = NULL;
   trialb = NULL;
   trialp = NULL;
   
   if (bestProb==0.0) {
      loadCandidateData(output, input.fsig, input.period, input.moddepth, input.ra, input.dec, input.stat, input.h0, input.prob, input.proberrcode, input.normalization);
   } else {
      loadCandidateData(output, bestf, bestp, bestdf, input.ra, input.dec, bestR, besth0, bestProb, bestproberrcode, input.normalization);
   }
   
}


//A brute force template search in a region of parameter space
/// Testing in progress
void templateSearch_scox1Style(candidateVector **output, REAL8 fminimum, REAL8 fspan, REAL8 period, REAL8 asini, inputParamsStruct *params, REAL4Vector *ffdata, INT4Vector *sftexist, REAL4Vector *aveNoise, REAL4Vector *aveTFnoisePerFbinRatio, REAL4FFTPlan *secondFFTplan, INT4 useExactTemplates)
{
   
   INT4 ii;
   REAL8Vector *trialf;
   REAL8 fstepsize;
   
   //Set up parameters of signal frequency search
   INT4 numfsteps = (INT4)round(2.0*fspan*params->Tcoh)+1;
   trialf = XLALCreateREAL8Vector(numfsteps);
   if (trialf==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL8Vector(%d) failed.\n", __func__, numfsteps);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   fstepsize = fspan/(REAL8)(numfsteps-1);
   for (ii=0; ii<numfsteps; ii++) trialf->data[ii] = fminimum + fstepsize*ii;
   
   //Now search over the frequencies
   INT4 proberrcode = 0;
   candidate cand;
   templateStruct *template = new_templateStruct(params->maxtemplatelength);
   if (template==NULL) {
      fprintf(stderr,"%s: new_templateStruct(%d) failed.\n", __func__, params->maxtemplatelength);
      XLAL_ERROR_VOID(XLAL_EFUNC); 
   }
   
   //Search over frequency
   for (ii=0; ii<(INT4)trialf->length; ii++) {
      //Determine modulation depth
      REAL8 moddepth = 0.8727*(trialf->data[ii]/1000.0)*(7200.0/period)*asini;

      //load candidate
      loadCandidateData(&cand, trialf->data[ii], period, moddepth, 0.0, 0.0, 0, 0, 0.0, 0, 0.0);

      //Make the template
      resetTemplateStruct(template);
      if (useExactTemplates!=0) {
         makeTemplate(template, cand, params, sftexist, secondFFTplan);
         if (xlalErrno!=0) {
            fprintf(stderr,"%s: makeTemplate() failed.\n", __func__);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
      } else {
         makeTemplateGaussians(template, cand, params, (INT4)aveTFnoisePerFbinRatio->length, (INT4)aveNoise->length);
         if (xlalErrno!=0) {
            fprintf(stderr,"%s: makeTemplateGaussians() failed.\n", __func__);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
      }

      REAL8 R = calculateR(ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
      if (XLAL_IS_REAL8_FAIL_NAN(R)) {
        fprintf(stderr,"%s: calculateR() failed.\n", __func__);
        XLAL_ERROR_VOID(XLAL_EFUNC);
      }
      REAL8 prob = probR(template, aveNoise, aveTFnoisePerFbinRatio, R, params, &proberrcode);
      if (XLAL_IS_REAL8_FAIL_NAN(prob)) {
         fprintf(stderr,"%s: probR() failed.\n", __func__);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      }
      REAL8 h0 = 0.0;
      if ( R > 0.0 ) h0 = 2.7426*pow(R/(params->Tcoh*params->Tobs),0.25);

      //Resize the output candidate vector if necessary
      if ((*output)->numofcandidates == (*output)->length-1) {
         *output = resize_candidateVector(*output, 2*((*output)->length));
         if (*output==NULL) {
            fprintf(stderr,"%s: resize_candidateVector(%d) failed.\n", __func__, 2*((*output)->length));
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
      }

      loadCandidateData(&((*output)->data[(*output)->numofcandidates]), trialf->data[ii], period, moddepth, 0.0, 0.0, R, h0, prob, proberrcode, 0.0);
      (*output)->numofcandidates++;
      
   } /* for ii < trialf */
   free_templateStruct(template);
   template = NULL;
   XLALDestroyREAL8Vector(trialf);
   trialf = NULL;
   
}


//Untested "efficient" template search. Not ready for prime-time
void efficientTemplateSearch(candidate *output, candidate input, REAL8 fminimum, REAL8 fmaximum, REAL8 minfstep, INT4 numperiods, REAL8 dfmin, REAL8 dfmax, REAL8 minDfstep, inputParamsStruct *params, REAL4Vector *ffdata, INT4Vector *sftexist, REAL4Vector *aveNoise, REAL4Vector *aveTFnoisePerFbinRatio, REAL4FFTPlan *secondFFTplan, INT4 useExactTemplates)
{
   
   INT4 bestproberrcode = 0, ii, jj, kk;
   REAL8 bestf = input.fsig, bestp = input.period, bestdf = input.moddepth, bestR = input.stat, besth0 = input.h0, bestProb = input.prob, fstepsize = 0.25*(fmaximum-fminimum), dfstepsize = 0.25*(dfmax-dfmin);
   REAL4 tcohfactor = 1.49e-3*params->Tcoh + 1.76;
   candidate cand;
   
   templateStruct *template = new_templateStruct(params->maxtemplatelength);
   if (template==NULL) {
      fprintf(stderr,"%s: new_templateStruct(%d) failed.\n", __func__, params->maxtemplatelength);
      XLAL_ERROR_VOID(XLAL_EFUNC); 
   }
   farStruct *farval = NULL;
   if (params->calcRthreshold) {
      farval = new_farStruct();
      if (farval==NULL) {
         fprintf(stderr,"%s: new_farStruct() failed.\n", __func__);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      }
   }
   
   //Search over numperiods different periods
   REAL8Vector *trialp = XLALCreateREAL8Vector(numperiods);
   if (trialp==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL8Vector(%d) failed.\n", __func__, numperiods);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   INT4 midposition = (INT4)roundf((numperiods-1)*0.5), proberrcode = 0;
   trialp->data[midposition] = input.period;
   
   REAL8Vector *trialf = XLALCreateREAL8Vector(2);
   if (trialf==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL8Vector(%d) failed.\n", __func__, 2);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   REAL8Vector *trialb = XLALCreateREAL8Vector(2);
   if (trialb==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL8Vector(%d) failed.\n", __func__, 2);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   
   while (fstepsize-minfstep>1.0e-8 && dfstepsize-minDfstep>1.0e-8) {
      //Initial point with template
      loadCandidateData(&cand, bestf, bestp, bestdf, input.ra, input.dec, 0, 0, 0.0, 0, 0.0);
      if (useExactTemplates!=0) {
         makeTemplate(template, cand, params, sftexist, secondFFTplan);
         if (xlalErrno!=0) {
            fprintf(stderr,"%s: makeTemplate() failed.\n", __func__);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
      } else {
         makeTemplateGaussians(template, cand, params, (INT4)aveTFnoisePerFbinRatio->length, (INT4)aveNoise->length);
         if (xlalErrno!=0) {
            fprintf(stderr,"%s: makeTemplateGaussians() failed.\n", __func__);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
      }
      if (params->calcRthreshold && bestProb==0.0) {
         numericFAR(farval, template, params->templatefar, aveNoise, aveTFnoisePerFbinRatio, params, params->rootFindingMethod);
         if (xlalErrno!=0) {
            fprintf(stderr,"%s: numericFAR() failed.\n", __func__);
            XLAL_ERROR_VOID(XLAL_EFUNC);
         }
      }
      bestR = calculateR(ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
      if (XLAL_IS_REAL8_FAIL_NAN(bestR)) {
         fprintf(stderr,"%s: calculateR() failed.\n", __func__);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      }
      bestProb = probR(template, aveNoise, aveTFnoisePerFbinRatio, bestR, params, &bestproberrcode);
      if (XLAL_IS_REAL8_FAIL_NAN(bestProb)) {
         fprintf(stderr,"%s: probR() failed.\n", __func__);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      }
      besth0 = 2.7426*pow(bestR/(params->Tcoh*params->Tobs),0.25);
      
      for (ii=0; ii<(INT4)trialf->length; ii++) {
         if (ii==0) trialf->data[ii] = bestf - fstepsize;
         else trialf->data[ii] = bestf + fstepsize;
      }
      for (ii=0; ii<(INT4)trialb->length; ii++) {
         if (ii==0) trialb->data[ii] = bestdf - dfstepsize;
         else trialb->data[ii] = bestdf + dfstepsize;
      }
      
      INT4 movedtoabetterpoint = 0;
      for (ii=0; ii<(INT4)trialf->length; ii++) {
         for (jj=0; jj<(INT4)trialb->length; jj++) {
            for (kk=0; kk<midposition; kk++) {
               REAL8 nnp = trialp->data[midposition+kk]*trialp->data[midposition+kk]*(1+trialp->data[midposition+kk]/tcohfactor/params->Tobs)/tcohfactor/params->Tobs*sqrt(3.6e-3/trialb->data[jj]);
               trialp->data[midposition+(kk+1)] = trialp->data[midposition+kk] + nnp;
               nnp = trialp->data[midposition-kk]*trialp->data[midposition-kk]*(1+trialp->data[midposition-kk]/tcohfactor/params->Tobs)/tcohfactor/params->Tobs*sqrt(3.6e-3/trialb->data[jj]);
               trialp->data[midposition-(kk+1)] = trialp->data[midposition-kk] - nnp;
            }
            for (kk=0; kk<(INT4)trialp->length; kk++) {
               if ( (trialf->data[ii]-trialb->data[jj]-6.0/params->Tcoh)>params->fmin && 
                   (trialf->data[ii]+trialb->data[jj]+6.0/params->Tcoh)<(params->fmin+params->fspan) && 
                   trialb->data[jj]<maxModDepth(trialp->data[kk], params->Tcoh) && 
                   trialp->data[kk]>minPeriod(trialb->data[jj], params->Tcoh) && 
                   trialp->data[kk]<=(0.2*params->Tobs) && 
                   trialp->data[kk]>=(2.0*3600.0) && 
                   trialb->data[jj]>=params->dfmin && 
                   trialb->data[jj]<=params->dfmax && 
                   trialp->data[kk]<=params->Pmax && 
                   trialp->data[kk]>=params->Pmin && 
                   trialf->data[ii]>=fminimum && 
                   trialf->data[ii]<=fmaximum && 
                   trialb->data[jj]>=dfmin && 
                   trialb->data[jj]<=dfmax ) {
                  
                  loadCandidateData(&cand, trialf->data[ii], trialp->data[kk], trialb->data[jj], input.ra, input.dec, 0, 0, 0.0, 0, 0.0);
                  
                  resetTemplateStruct(template);
                  
                  if (useExactTemplates!=0) {
                     makeTemplate(template, cand, params, sftexist, secondFFTplan);
                     if (xlalErrno!=0) {
                        fprintf(stderr,"%s: makeTemplate() failed.\n", __func__);
                        XLAL_ERROR_VOID(XLAL_EFUNC);
                     }
                  } else {
                     makeTemplateGaussians(template, cand, params, (INT4)aveTFnoisePerFbinRatio->length, (INT4)aveNoise->length);
                     if (xlalErrno!=0) {
                        fprintf(stderr,"%s: makeTemplateGaussians() failed.\n", __func__);
                        XLAL_ERROR_VOID(XLAL_EFUNC);
                     }
                  }
                  
                  if (params->calcRthreshold && bestProb==0.0) {
                     numericFAR(farval, template, params->templatefar, aveNoise, aveTFnoisePerFbinRatio, params, params->rootFindingMethod);
                     if (xlalErrno!=0) {
                        fprintf(stderr,"%s: numericFAR() failed.\n", __func__);
                        XLAL_ERROR_VOID(XLAL_EFUNC);
                     }
                  }
                  
                  REAL8 R = calculateR(ffdata, template, aveNoise, aveTFnoisePerFbinRatio);
                  if (XLAL_IS_REAL8_FAIL_NAN(R)) {
                     fprintf(stderr,"%s: calculateR() failed.\n", __func__);
                     XLAL_ERROR_VOID(XLAL_EFUNC);
                  }
                  REAL8 prob = probR(template, aveNoise, aveTFnoisePerFbinRatio, R, params, &proberrcode);
                  if (XLAL_IS_REAL8_FAIL_NAN(prob)) {
                     fprintf(stderr,"%s: probR() failed.\n", __func__);
                     XLAL_ERROR_VOID(XLAL_EFUNC);
                  }
                  
                  REAL8 h0 = 2.7426*pow(R/(params->Tcoh*params->Tobs),0.25);
                  
                  if ( prob < bestProb ) {
                     bestf = trialf->data[ii];
                     bestp = trialp->data[kk];
                     bestdf = trialb->data[jj];
                     bestR = R;
                     besth0 = h0;
                     bestProb = prob;
                     bestproberrcode = proberrcode;
                     movedtoabetterpoint = 1;
                  }
                  
               } // if within boundaries
            } // for kk < trialp->length
         } // for jj < trialb->length
      } // for ii < trialf->length
      
      if (movedtoabetterpoint==0) {
         fstepsize *= 0.5;
         dfstepsize *= 0.5;
      }
      
   } // while
   
   
   loadCandidateData(output, bestf, bestp, bestdf, input.ra, input.dec, bestR, besth0, bestProb, bestproberrcode, input.normalization);
   
   
   free_templateStruct(template);
   template = NULL;
   if (params->calcRthreshold) {
      free_farStruct(farval);
      farval = NULL;
   }
   XLALDestroyREAL8Vector(trialp);
   XLALDestroyREAL8Vector(trialf);
   XLALDestroyREAL8Vector(trialb);
   
}


//////////////////////////////////////////////////////////////
// Does the insertion sort for the template weights
void insertionSort_template(templateStruct *output, REAL4 weight, INT4 pixelloc, INT4 firstfftfreq, INT4 secfftfreq)
{
   
   //INT4 ii;
   INT4 insertionpoint = (INT4)output->templatedata->length-1;
   INT4 numbertomove = 0;
   while (insertionpoint > 0 && weight > output->templatedata->data[insertionpoint-1]) {
      insertionpoint--;
      numbertomove++;
   }
   
   /* for (ii=output->templatedata->length-1; ii>insertionpoint; ii--) {
      output->templatedata->data[ii] = output->templatedata->data[ii-1];
      output->pixellocations->data[ii] = output->pixellocations->data[ii-1];
      output->firstfftfrequenciesofpixels->data[ii] = output->firstfftfrequenciesofpixels->data[ii-1];
      output->secondfftfrequencies->data[ii] = output->secondfftfrequencies->data[ii-1];
   } */
   if (insertionpoint<(INT4)output->templatedata->length-1) {
      memmove(&(output->templatedata->data[insertionpoint+1]), &(output->templatedata->data[insertionpoint]), sizeof(REAL4)*numbertomove);
      memmove(&(output->pixellocations->data[insertionpoint+1]), &(output->pixellocations->data[insertionpoint]), sizeof(INT4)*numbertomove);
      memmove(&(output->firstfftfrequenciesofpixels->data[insertionpoint+1]), &(output->firstfftfrequenciesofpixels->data[insertionpoint]), sizeof(INT4)*numbertomove);
      memmove(&(output->secondfftfrequencies->data[insertionpoint+1]), &(output->secondfftfrequencies->data[insertionpoint]), sizeof(INT4)*numbertomove);
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
   
   if (fabs(x*x-1.0)<1.0e-8) return -0.5;
   if (fabs(x)<1.0e-8) return -1.0;
   
   REAL8 pix = LAL_PI*x;
   return sin(pix)/(pix*(x*x-1.0));
   
} /* sincxoverxsqminusone() */
REAL8 sqsincxoverxsqminusone(REAL8 x)
{
   
   REAL8 val = sincxoverxsqminusone(x);
   return val*val;
   
} /* sqsincxoverxsqminusone() */


//Stolen from computeFstat.c with higher resolution
#define OOTWOPI         (1.0 / LAL_TWOPI)
int twospect_sin_cos_LUT(REAL8 *sinx, REAL8 *cosx, REAL8 x)
{
   return twospect_sin_cos_2PI_LUT( sinx, cosx, x * OOTWOPI );
} /* twospect_sin_cos_LUT() */
#define LUT_RES         1024      /* resolution of lookup-table */
#define LUT_RES_F       (1.0 * LUT_RES)
#define OO_LUT_RES      (1.0 / LUT_RES)
#define X_TO_IND        (1.0 * LUT_RES * OOTWOPI )
#define IND_TO_X        (LAL_TWOPI * OO_LUT_RES)
#define TRUE (1==1)
#define FALSE (1==0)
int twospect_sin_cos_2PI_LUT(REAL8 *sin2pix, REAL8 *cos2pix, REAL8 x)
{
   REAL8 xt;
   INT4 i0;
   REAL8 d, d2;
   REAL8 ts, tc;
   REAL8 dummy;
   
   static BOOLEAN firstCall = TRUE;
   static REAL8 sinVal[LUT_RES+1], cosVal[LUT_RES+1];
   
   /* the first time we get called, we set up the lookup-table */
   if ( firstCall ) {
      UINT4 k;
      for (k=0; k <= LUT_RES; k++) {
         sinVal[k] = sin( LAL_TWOPI * k * OO_LUT_RES );
         cosVal[k] = cos( LAL_TWOPI * k * OO_LUT_RES );
      }
      firstCall = FALSE;
   }

   /* we only need the fractional part of 'x', which is number of cylces,
   * this was previously done using
   *   xt = x - (INT4)x;
   * which is numerically unsafe for x > LAL_INT4_MAX ~ 2e9
   * for saftey we therefore rather use modf(), even if that
   * will be somewhat slower...
   */
   //xt = modf(x, &dummy);/* xt in (-1, 1) */
   if (x<2.147483647e9 && x>-2.147483647e9) xt = x - (INT4)x;  // if x < LAL_INT4_MAX
   else xt = modf(x, &dummy);
   
   if ( xt < 0.0 ) xt += 1.0;                  /* xt in [0, 1 ) */
   #ifndef LAL_NDEBUG
      if ( xt < 0.0 || xt > 1.0 ) {
         XLALPrintError("\nFailed numerica in twospect_sin_cos_2PI_LUT(): xt = %f not in [0,1)\n\n", xt );
         return XLAL_FAILURE;
      }
   #endif
   
   i0 = (INT4)( xt * LUT_RES_F + 0.5 );  /* i0 in [0, LUT_RES ] */
   d = d2 = LAL_TWOPI * (xt - OO_LUT_RES * i0);
   d2 *= 0.5 * d;
   
   ts = sinVal[i0];
   tc = cosVal[i0];
   
   /* use Taylor-expansions for sin/cos around LUT-points */
   (*sin2pix) = ts + d * tc - d2 * ts;
   (*cos2pix) = tc - d * ts - d2 * tc;
   
   return XLAL_SUCCESS;
} /* twospect_sin_cos_2PI_LUT() */




