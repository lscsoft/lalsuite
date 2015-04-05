/*
*  Copyright (C) 2010, 2011, 2014 Evan Goetz
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

#include <gsl/gsl_sort.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_fit.h>
#include "falsealarm.h"
#include "statistics.h"
#include "cdfwchisq.h"
#include "vectormath.h"

/**
 * Allocate memory for farStruct
 * \return Pointer to a farStruct
 */
farStruct * new_farStruct(void)
{
   farStruct *farstruct = NULL;
   XLAL_CHECK_NULL( (farstruct = XLALMalloc(sizeof(*farstruct))) != NULL, XLAL_ENOMEM );
   farstruct->far = 1.0;
   farstruct->topRvalues = NULL;
   return farstruct;
} /* new_farStruct() */


/**
 * Destroy an farStruct
 * \param [in] farstruct Pointer to a farStruct
 */
void free_farStruct(farStruct *farstruct)
{
   XLALDestroyREAL4Vector(farstruct->topRvalues);
   farstruct->topRvalues = NULL;
   XLALFree((farStruct*)farstruct);
} /* free_farStruct() */


/**
 * Estimate the FAR of the R statistic from the weights by a number of trials
 * \param [out] output         Pointer to a farStruct
 * \param [in]  template Pointer to a TwoSpectTemplate containing a template
 * \param [in]  trials         Number of trials to estimate the FAR
 * \param [in]  thresh         Threshold value
 * \param [in]  ffplanenoise   Pointer to REAL4Vector containing the expected 2nd FFT background estimates
 * \param [in]  fbinaveratios  Pointer to REAL4Vector containing the noise floor variations
 * \return Status value
 */
INT4 estimateFAR(farStruct *output, TwoSpectTemplate *template, INT4 trials, REAL8 thresh, REAL4Vector *ffplanenoise, REAL4Vector *fbinaveratios)
{

   UINT4 numofweights = 0;
   for (UINT4 ii=0; ii<template->templatedata->length; ii++) if (template->templatedata->data[ii]!=0.0) numofweights++;

   REAL8 sumofsqweights = 0.0;
   for (UINT4 ii=0; ii<numofweights; ii++) sumofsqweights += (template->templatedata->data[ii]*template->templatedata->data[ii]);
   REAL8 sumofsqweightsinv = 1.0/sumofsqweights;

   REAL4Vector *Rs = NULL;
   XLAL_CHECK( (Rs = XLALCreateREAL4Vector(trials)) != NULL, XLAL_EFUNC );

   gsl_rng *rng = NULL;
   XLAL_CHECK( (rng = gsl_rng_alloc(gsl_rng_mt19937)) != NULL, XLAL_EFUNC );
   //srand(time(NULL));
   //UINT8 randseed = rand();
   //gsl_rng_set(rng, randseed);
   gsl_rng_set(rng, 0);

   UINT4 numfprbins = ffplanenoise->length;

   for (INT4 ii=0; ii<trials; ii++) {
      //Create noise value and R value
      REAL8 R = 0.0;
      for (UINT4 jj=0; jj<numofweights; jj++) {
         UINT4 firstfreqbin = template->pixellocations->data[jj]/numfprbins;
         UINT4 secfreqbin = template->pixellocations->data[jj] - firstfreqbin*numfprbins;
         REAL8 noise = expRandNum(ffplanenoise->data[secfreqbin]*fbinaveratios->data[firstfreqbin], rng);
         XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
         R += (noise - ffplanenoise->data[secfreqbin]*fbinaveratios->data[firstfreqbin])*template->templatedata->data[jj];
      }
      Rs->data[ii] = (REAL4)(R*sumofsqweightsinv);
   } /* for ii < trials */
   REAL4 mean = calcMean(Rs);
   REAL4 sigma = 0.0;
   XLAL_CHECK( calcStddev(&sigma, Rs) == XLAL_SUCCESS, XLAL_EFUNC );

   //Do an insertion sort. At best this is O(thresh*trials), at worst this is O(thresh*trials*trials).
   if (output->topRvalues == NULL) XLAL_CHECK( (output->topRvalues = XLALCreateREAL4Vector((INT4)round(thresh*trials)+1)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( gsl_sort_float_largest((float*)output->topRvalues->data, output->topRvalues->length, (float*)Rs->data, 1, Rs->length) == 0, XLAL_EFUNC );

   output->far = output->topRvalues->data[output->topRvalues->length - 1];
   output->distMean = mean;
   output->distSigma = sigma;

   //Destroy
   XLALDestroyREAL4Vector(Rs);
   gsl_rng_free(rng);

   return XLAL_SUCCESS;

} /* estimateFAR() */


/**
 * Numerically solve for the FAR of the R statistic from the weights using the Davies algorithm and a root finding algorithm
 * \param [out] output         Pointer to a farStruct
 * \param [in]  template       Pointer to a TwoSpectTemplate containing a template
 * \param [in]  thresh         Threshold value
 * \param [in]  ffplanenoise   Pointer to REAL4Vector containing the expected 2nd FFT background estimates
 * \param [in]  fbinaveratios  Pointer to REAL4Vector containing the noise floor variations
 * \param [in]  inputParams    Pointer to UserInput_t
 * \param [in]  rng            Pointer to gsl_rng
 * \param [in]  method         Integer value of 0 (Brent's method) or 1 (Newton's method)
 * \return Status value
 */
INT4 numericFAR(farStruct *output, TwoSpectTemplate *template, REAL8 thresh, REAL4Vector *ffplanenoise, REAL4Vector *fbinaveratios, UserInput_t *inputParams, gsl_rng *rng, INT4 method)
{

   XLAL_CHECK( output != NULL && template != NULL && ffplanenoise != NULL && fbinaveratios != NULL && inputParams != NULL, XLAL_EINVAL );

   INT4 ii;
   INT4 errcode = 0;

   //Set up solver: method 0 is Brent's method, method 1 is Newton's method
   const gsl_root_fsolver_type *T1 = gsl_root_fsolver_brent;
   gsl_root_fsolver *s1 = NULL;
   XLAL_CHECK( (s1 = gsl_root_fsolver_alloc(T1)) != NULL, XLAL_EFUNC );
   gsl_function F;
   const gsl_root_fdfsolver_type *T0 = gsl_root_fdfsolver_newton;
   gsl_root_fdfsolver *s0 = NULL;
   XLAL_CHECK( (s0 = gsl_root_fdfsolver_alloc(T0)) != NULL, XLAL_EFUNC );
   gsl_function_fdf FDF;


   //Include the various parameters in the struct required by GSL
   struct gsl_probR_pars params = {template, ffplanenoise, fbinaveratios, thresh, inputParams, rng, errcode};

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
   if (method != 0) XLAL_CHECK( gsl_root_fsolver_set(s1, &F, Rlow, Rhigh) == GSL_SUCCESS, XLAL_EFUNC );
   else XLAL_CHECK( gsl_root_fdfsolver_set(s0, &FDF, root) == GSL_SUCCESS, XLAL_EFUNC );

   //And now find the root
   ii = 0;
   INT4 max_iter = 100, jj = 0, max_retries = 10;
   INT4 status = GSL_CONTINUE;
   REAL8 prevroot = 0.0;
   while (status==GSL_CONTINUE && ii<max_iter) {

      ii++;

      if (method != 0) {
         status = gsl_root_fsolver_iterate(s1);
         XLAL_CHECK( status == GSL_CONTINUE || status == GSL_SUCCESS, XLAL_EFUNC, "gsl_root_fsolver_iterate() failed with code %d\n", status );
         if (ii>0) prevroot = root;
         root = gsl_root_fsolver_root(s1);
         Rlow = gsl_root_fsolver_x_lower(s1);
         Rhigh = gsl_root_fsolver_x_upper(s1);
         status = gsl_root_test_interval(Rlow, Rhigh, 0.0, 0.001);
         XLAL_CHECK( status == GSL_CONTINUE || status == GSL_SUCCESS, XLAL_EFUNC, "gsl_root_test_interval() failed with code %d\n", status );
      } else {
         status = gsl_root_fdfsolver_iterate(s0);
         XLAL_CHECK( status == GSL_CONTINUE || status == GSL_SUCCESS, XLAL_EFUNC, "gsl_root_fdfsolver_iterate() failed with code %d\n", status );
         prevroot = root;
         root = gsl_root_fdfsolver_root(s0);
         status = gsl_root_test_delta(prevroot, root, 0.0, 0.001);
         XLAL_CHECK( status == GSL_CONTINUE || status == GSL_SUCCESS, XLAL_EFUNC, "gsl_root_test_delta() failed with code %d\n", status );

         //If there is an issue that the root is negative, try a new initial guess
         if (root<0.0 && jj<max_retries) {
            ii = 0;
            jj++;
            status = GSL_CONTINUE;
            XLAL_CHECK( gsl_root_fdfsolver_set(s0, &FDF, gsl_rng_uniform_pos(rng)*Rhigh) == GSL_SUCCESS, XLAL_EFUNC );
         } else if (root<0.0 && jj==max_retries) {
            status = GSL_FAILURE;
         } //Up to here
      }
   } /* while status==GSL_CONTINUE && ii < max_iter */

   //Failure modes
   if (method != 0) {
      XLAL_CHECK( status == GSL_SUCCESS, XLAL_FAILURE, "Root finding iteration (%d/%d) failed with failure code %d. Previous root = %f, current root = %f\n", ii, max_iter, status, prevroot, root );
      XLAL_CHECK( ii<max_iter, XLAL_EMAXITER, "Root finding iteration (%d/%d) failed with failure code %d. Previous root = %f, current root = %f\n", ii, max_iter, status, prevroot, root );
      XLAL_CHECK( root != 0.0, XLAL_ERANGE, "Root finding iteration (%d/%d) converged to 0.0\n", ii, max_iter );
   } else {
      XLAL_CHECK( status == GSL_SUCCESS, XLAL_FAILURE, "Root finding iteration (%d/%d) failed with failure code %d. Previous root = %f, current root = %f\n", ii, max_iter, status, prevroot, root );
      XLAL_CHECK( ii<max_iter, XLAL_EMAXITER, "Root finding iteration (%d/%d) failed with failure code %d. Previous root = %f, current root = %f\n", ii, max_iter, status, prevroot, root );
      XLAL_CHECK( root > 0.0, XLAL_ERANGE, "Threshold value found (%f) is less than 0.0\n", root );
   }

   output->far = root;
   output->distMean = 0.0;
   output->distSigma = 0.0; //Fake the value of sigma
   output->farerrcode = errcode;

   //Cleanup
   gsl_root_fsolver_free(s1);
   gsl_root_fdfsolver_free(s0);

   return XLAL_SUCCESS;

} /* numericFAR() */


/**
 * For the root finding, calculating the false alarm probability of R. This method takes the average of 3 values of close by R values for stability
 * \param [in] R     R value of interest
 * \param [in] param A gsl_probR_pars struct
 * \return Difference between the false alarm probability and the log10 threshold value
 */
REAL8 gsl_probR(REAL8 R, void *param)
{

   struct gsl_probR_pars *pars = (struct gsl_probR_pars*)param;

   REAL8 dR = 0.005;
   REAL8 R1 = (1.0+dR)*R;
   REAL8 R2 = (1.0-dR)*R;
   INT4 errcode1 = 0, errcode2 = 0, errcode3 = 0;

   REAL8 prob = (probR(pars->template, pars->ffplanenoise, pars->fbinaveratios, R, pars->inputParams, pars->rng, &errcode1) + probR(pars->template, pars->ffplanenoise, pars->fbinaveratios, R1, pars->inputParams, pars->rng, &errcode2) + probR(pars->template, pars->ffplanenoise, pars->fbinaveratios, R2, pars->inputParams, pars->rng, &errcode3))/3.0;

   if (errcode1!=0) pars->errcode = errcode1;
   else if (errcode2!=0) pars->errcode = errcode2;
   else if (errcode3!=0) pars->errcode = errcode3;

   REAL8 returnval = prob - log10(pars->threshold);

   return returnval;

} /* gsl_probR() */


/**
 * Determine the slope of the inverse cumulative distribution function
 * \param [in] R     R value of interest
 * \param [in] param A gsl_probR_pars struct
 * \return The slope of the inverse distribution function
 */
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


/**
 * Determine the difference between the probability and log10 of the threshold as well as the slope of the inverse cumulative distribution function
 * \param [in]  R            R value of interest
 * \param [in]  param        A gsl_probR_pars struct
 * \param [out] probabilityR The difference between the threshold and the probability of the R value in question
 * \param [out] dprobRdR     The slope of the inverse cumulative distribution function
 * \return The slope of the inverse distribution function
 */
void gsl_probRandDprobRdR(REAL8 R, void *param, REAL8 *probabilityR, REAL8 *dprobRdR)
{

   struct gsl_probR_pars *pars = (struct gsl_probR_pars*)param;

   *probabilityR = gsl_probR(R, pars);

   *dprobRdR = gsl_dprobRdR(R, pars);

} /* gsl_probRandDprobRdR() */


/**
 * Analytically calculate the probability of a true signal using the Davies' method
 * \param [in]  template      Pointer to a TwoSpectTemplate with a template
 * \param [in]  ffplanenoise  Pointer to a REAL4Vector with an estimate of the background of 2nd FFT powers
 * \param [in]  fbinaveratios Pointer to a REAL4Vector of frequency bin average ratios
 * \param [in]  R             The value of R for a given template
 * \param [in]  params        Pointer to UserInput_t
 * \param [in]  rng           Pointer to gsl_rng
 * \param [out] errcode       Pointer to the error code value from the Davies algorithm
 * \return log10 false alarm probability value
 */
REAL8 probR(TwoSpectTemplate *template, REAL4Vector *ffplanenoise, REAL4Vector *fbinaveratios, REAL8 R, UserInput_t *params, gsl_rng *rng, INT4 *errcode)
{

   XLAL_CHECK_REAL8( template != NULL && ffplanenoise != NULL && fbinaveratios != NULL && params != NULL, XLAL_EINVAL );

   REAL8 prob = 0.0;
   REAL8 sumwsq = 0.0;
   INT4 numweights = 0;
   for (UINT4 ii=0; ii<template->templatedata->length; ii++) {
      if (template->templatedata->data[ii]!=0.0) {
         numweights++;
         sumwsq += template->templatedata->data[ii]*template->templatedata->data[ii];
      } else {
         break;
      }
   }

   UINT4 numfprbins = ffplanenoise->length;

   alignedREAL8Vector *newweights = NULL;
   INT4Vector *sorting = NULL;
   XLAL_CHECK_REAL8( (newweights = createAlignedREAL8Vector(numweights, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK_REAL8( (sorting = XLALCreateINT4Vector(numweights)) != NULL, XLAL_EFUNC );

   REAL8 Rpr = R;
   for (UINT4 ii=0; ii<newweights->length; ii++) {
      UINT4 firstfreqbin = template->pixellocations->data[ii]/numfprbins;
      UINT4 secfreqbin = template->pixellocations->data[ii] - firstfreqbin*numfprbins;
      newweights->data[ii] = 0.5*template->templatedata->data[ii]*ffplanenoise->data[secfreqbin]*fbinaveratios->data[firstfreqbin]/sumwsq;
      Rpr += template->templatedata->data[ii]*ffplanenoise->data[secfreqbin]*fbinaveratios->data[firstfreqbin]/sumwsq;
      sorting->data[ii] = ii;  //This is for the fact that a few steps later (before using Davies' algorithm, we sort the weights)
   }

   qfvars vars;
   vars.weights = newweights;
   vars.sorting = sorting;
   vars.dofs = NULL;
   vars.noncentrality = NULL;
   vars.arrayNotSorted = 0;           //Set because we do the sorting outside of Davies' algorithm with qsort
   vars.lim = 50000000;
   vars.c = Rpr;
   vars.vectorMath = params->vectorMath;
   REAL8 sigma = 0.0;
   REAL8 accuracy = 1.0e-11;   //(1e-5) old value

   //sort the weights here so we don't have to do it later (qsort)
   sort_double_ascend((REAL8Vector*)newweights);

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
      REAL8Vector *probvals = NULL, *cvals = NULL;
      XLAL_CHECK_REAL8( (probvals = XLALCreateREAL8Vector(20)) != NULL, XLAL_EFUNC );
      XLAL_CHECK_REAL8( (cvals = XLALCreateREAL8Vector(20)) != NULL, XLAL_EFUNC );

      for (UINT4 ii=0; ii<probvals->length; ii++) {
         c1 = gsl_rng_uniform_pos(rng)*(upperend-lowerend)+lowerend;
         vars.c = c1;
         tempprob = 1.0-cdfwchisq_twospect(&vars, sigma, accuracy, &errcode1);
         while (tempprob<=1.0e-8 || tempprob>=1.0e-6) {
            if (tempprob<=1.0e-8) upperend = c1;
            else if (tempprob>=1.0e-6) lowerend = c1;
            c1 = gsl_rng_uniform_pos(rng)*(upperend-lowerend)+lowerend;
            vars.c = c1;
            tempprob = 1.0-cdfwchisq_twospect(&vars, sigma, accuracy, &errcode1);

            INT4 tries = 1;
            while (tries<10 && errcode1 != 0) {
               tries++;
               c1 = gsl_rng_uniform_pos(rng)*(upperend-lowerend)+lowerend;
               vars.c = c1;
               tempprob = 1.0-cdfwchisq_twospect(&vars, sigma, accuracy, &errcode1);
            }
            if (tries>=10 && errcode1!=0) {
               fprintf(stderr,"%s: cdfwchisq_twospect() failed with code %d after making %d tries.\n", __func__, errcode1, tries);
               XLAL_ERROR_REAL8(XLAL_EFUNC);
            }

         }
         XLAL_CHECK_REAL8( errcode1 == 0, XLAL_EFUNC, "cdfwchisq_twospect() failed with code %d\n", errcode1 );
         probvals->data[ii] = log10(tempprob);
         cvals->data[ii] = c1;
      }

      REAL8 yintercept, cov00, cov01, cov11, sumsq;
      XLAL_CHECK_REAL8( gsl_fit_linear(cvals->data, 1, probvals->data, 1, cvals->length, &yintercept, &probslope, &cov00, &cov01, &cov11, &sumsq) == GSL_SUCCESS, XLAL_EFUNC );
      logprobest = probslope*Rpr + yintercept;
      XLAL_CHECK_REAL8( logprobest<=-0.5, XLAL_ERANGE, "Failure calculating accurate interpolated value\n" );

      XLALDestroyREAL8Vector(probvals);
      XLALDestroyREAL8Vector(cvals);

      *errcode = errcode1;

   }

   //If errcode is still != 0, better fail
   XLAL_CHECK_REAL8( *errcode == 0, XLAL_EFUNC, "cdfwchisq_twospect() failed at the end with code %d\n", *errcode );

   //Cleanup
   destroyAlignedREAL8Vector(newweights);
   XLALDestroyINT4Vector(sorting);

   if (estimatedTheProb==1) return logprobest;
   else return log10(prob);

} /* probR() */
