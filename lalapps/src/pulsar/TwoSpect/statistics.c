/*
 *  Copyright (C) 2011 Evan Goetz
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

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_statistics_double.h>

#include <lal/LALConstants.h>
#include <lal/VectorOps.h>

#include "statistics.h"
#include "fastchisqinv.h"


//////////////////////////////////////////////////////////////
// Create a exponentially distributed noise value  -- done
REAL8 expRandNum(REAL8 mu, gsl_rng *ptrToGenerator)
{
   
   if (mu<=0.0) {
      fprintf(stderr,"%s: expRandNum(%f, %p) failed.\n", __func__, mu, ptrToGenerator);
      XLAL_ERROR_REAL8(XLAL_EINVAL);
   } else if (ptrToGenerator==NULL) {
      fprintf(stderr,"%s: expRandNum(%f, %p) failed.\n", __func__, mu, ptrToGenerator);
      XLAL_ERROR_REAL8(XLAL_EFAULT);
   }
   
   return gsl_ran_exponential(ptrToGenerator, mu);
   
} /* expRandNum() */


//Compute the CDF P value at value x of a chi squared PDF with nu degrees of freedom
//Rougly REAL4 precision
REAL8 twospect_cdf_chisq_P(REAL8 x, REAL8 nu)
{
   if (XLAL_IS_REAL8_FAIL_NAN(x) || XLAL_IS_REAL8_FAIL_NAN(nu)) {
      fprintf(stderr,"%s: Invalid arguments x=%f, nu = %f.\n", __func__, x, nu);
      XLAL_ERROR_REAL8(XLAL_EINVAL);
   }
   REAL8 val = cdf_gamma_P(x, 0.5*nu, 2.0);
   if (XLAL_IS_REAL8_FAIL_NAN(val)) {
      fprintf(stderr,"%s: cdf_gamma_P(%f, %f, 2.0) failed.\n", __func__, x, 0.5*nu);
      XLAL_ERROR_REAL8(XLAL_EFUNC);
   }
   return val;
}

//Compute the CDF P value at value x of a chi squared PDF with nu degrees of freedom using the Matlab-based function
REAL8 matlab_cdf_chisq_P(REAL8 x, REAL8 nu)
{
   if (XLAL_IS_REAL8_FAIL_NAN(x) || XLAL_IS_REAL8_FAIL_NAN(nu)) {
      fprintf(stderr,"%s: Invalid arguments x=%f, nu = %f.\n", __func__, x, nu);
      XLAL_ERROR_REAL8(XLAL_EINVAL);
   }
   REAL8 val = cdf_gamma_P_usingmatlab(x, 0.5*nu, 2.0);
   if (XLAL_IS_REAL8_FAIL_NAN(val)) {
      fprintf(stderr,"%s: cdf_gamma_P_usingmatlab(%f, %f, 2.0) failed.\n", __func__, x, 0.5*nu);
      XLAL_ERROR_REAL8(XLAL_EFUNC);
   }
   return val;
}


//Matlab's version of the non-central chih squared CDF with nu degrees of freedom and non-centrality delta at value x.
REAL8 ncx2cdf(REAL8 x, REAL8 dof, REAL8 delta)
{

   REAL8 prob = 0.0;

   //Fail for bad inputs or return 0 if x<=0
   if (dof<0.0 || delta<0.0) {
      fprintf(stderr,"%s: Invalid arguments x=%f, dof=%f, delta=%f.\n", __func__, x, dof, delta);
      XLAL_ERROR_REAL8(XLAL_EINVAL);
   } else if (x<=0.0) {
      return prob;
   }

   REAL8 err = LAL_REAL8_EPS;
   REAL8 halfdelta = 0.5*delta;
   INT8 counter = (INT8)floor(halfdelta);
   REAL8 P = gsl_ran_poisson_pdf(counter, halfdelta);
   REAL8 C = gsl_cdf_chisq_P(x, dof+2.0*counter);
   if (XLAL_IS_REAL8_FAIL_NAN(C)) {
      fprintf(stderr, "%s: gsl_cdf_chisq_P(%f, %f) failed.\n", __func__, x, dof+2.0*counter);
      XLAL_ERROR_REAL4(XLAL_EFUNC);
   }
   REAL8 E = exp((dof*0.5+counter-1.0)*log(x*0.5) - x*0.5 - lgamma(dof*0.5+counter));

   sumseries_eg(&prob, P, C, E, counter, x, dof, halfdelta, err, 0);
   if (xlalErrno!=0) {
      fprintf(stderr,"%s: sumseries_eg(%f, %f, %f, %f, %f, %f, %f, %f, %f, 0) failed.\n", __func__, prob, P, C, E, floor(halfdelta), x, dof, halfdelta, err);
      XLAL_ERROR_REAL8(XLAL_EFUNC);
   }
   counter--;
   if (counter<0) return fmin(prob, 1.0);

   sumseries_eg(&prob, P, C, E, counter, x, dof, halfdelta, err, 1);
   if (xlalErrno!=0) {
      fprintf(stderr,"%s: sumseries_eg(%f, %f, %f, %f, %f, %f, %f, %f, %f, 1) failed.\n", __func__, prob, P, C, E, floor(halfdelta), x, dof, halfdelta, err);
      XLAL_ERROR_REAL8(XLAL_EFUNC);
   }

   //This part computes small probabilities
   INT4 fromzero = 0;
   if (prob==0.0) fromzero = 1;
   if (fromzero==1) {
      counter = 0;
      REAL8 pk = gsl_ran_poisson_pdf(0, halfdelta)*gsl_cdf_chisq_P(x, dof);
      REAL8 dp = 0.0;
      INT4 ok = 0;
      if ((REAL8)counter<halfdelta) ok = 1;
      while (ok==1) {
         counter++;
         P = gsl_ran_poisson_pdf(counter, halfdelta);
         C = gsl_cdf_chisq_P(x, dof+2.0*counter);
         if (XLAL_IS_REAL8_FAIL_NAN(C)) {
            fprintf(stderr, "%s: gsl_cdf_chisq_P(%f, %f) failed.\n", __func__, x, dof+2.0*counter);
            XLAL_ERROR_REAL4(XLAL_EFUNC);
         }
         dp = P*C;
         pk += dp;
         if (!(ok==1 && (REAL8)counter<halfdelta && dp>=err*pk)) ok = 0;
      }
      prob = pk;
   }

   return fmin(prob, 1.0);

}

//Matlab's sumseries function
void sumseries(REAL8 *computedprob, REAL8 P, REAL8 C, REAL8 E, INT8 counter, REAL8 x, REAL8 dof, REAL8 halfdelta, REAL8 err, INT4 countdown)
{

   //Exit with error if halfdelta = 0.0
   if (halfdelta==0.0) {
      fprintf(stderr, "%s: halfdelta=0.0\n", __func__);
      XLAL_ERROR_VOID(XLAL_EINVAL);
   }

   REAL8 Pint = P, Cint = C, Eint = E;
   INT8 counterint = counter;
   INT8 j = 0;
   if (countdown!=0) {
      if (counterint>=0) j = 1;
      if (j==1) {
         Pint *= (counterint+1.0)/halfdelta;
         Cint += E;
      } else counterint = -1;
   }

   while (counterint!=-1) {
      REAL8 pplus = Pint*Cint;
      if (XLAL_IS_REAL8_FAIL_NAN(pplus)) {
         fprintf(stderr, "%s: pplus is NaN.\n", __func__);
         XLAL_ERROR_VOID(XLAL_EFPOVRFLW);
      }
      *(computedprob) += pplus;
      
      if (pplus > *(computedprob)*err) j = 1;
      else j = 0;
      if (countdown!=0 && counterint<0) j = 0;
      if (j==0) return;
      
      if (countdown!=0) {
         counterint--;
         Pint *= (counterint+1.0)/halfdelta;
         Eint *= (0.5*dof + counterint+1.0)/(x*0.5);
         Cint += Eint;
      } else {
         counterint++;
         Pint *= halfdelta/counterint;
         Eint *= (0.5*x)/(0.5*dof+counterint-1.0);
         Cint -= Eint;
      }
   }

}


//Evan's sumseries function based on matlab's sumseries() version above, but faster
void sumseries_eg(REAL8 *computedprob, REAL8 P, REAL8 C, REAL8 E, INT8 counter, REAL8 x, REAL8 dof, REAL8 halfdelta, REAL8 err, INT4 countdown)
{

   //If halfdelta = 0.0, then exit with error
   if (halfdelta==0.0) {
      fprintf(stderr, "%s: halfdelta=0.0\n", __func__);
      XLAL_ERROR_VOID(XLAL_EINVAL);
   }

   REAL8 Pint = P, Cint = C, Eint = E;
   INT8 counterint = counter;
   REAL8 oneoverhalfdelta = 1.0/halfdelta;   //pre-compute
   REAL8 halfdof = 0.5*dof;                  //pre-compute
   REAL8 halfx = 0.5*x;                      //pre-compute
   
   if (countdown!=0) {
      if (counterint>=0) {
         Pint *= (counterint+1.0)*oneoverhalfdelta;
         Cint += E;
      } else counterint = -1;
   }
   
   if (counterint==-1) return;
   else if (countdown!=0) {
      REAL8 oneoverhalfx = 1.0/halfx;
      while (counterint!=-1) {
         REAL8 pplus = Pint*Cint;
         *(computedprob) += pplus;

         if (pplus<=*(computedprob)*err || counterint<0) return;

         counterint--;
         Pint *= (counterint+1)*oneoverhalfdelta;
         Eint *= (halfdof + counterint+1)*oneoverhalfx;
         Cint += Eint;

         //This part below is wrong!
         /* Pint *= (counterint)*oneoverhalfdelta;
         Eint *= (halfdof + counterint)*oneoverhalfx;
         Cint += Eint;
         counterint--; */
      }
   } else {
      while (counterint!=-1) {
         REAL8 pplus = Pint*Cint;
         *(computedprob) += pplus;

         if (pplus<=*(computedprob)*err) return;

         counterint++;
         Pint *= halfdelta/counterint;
         Eint *= halfx/(halfdof+counterint-1);
         Cint -= Eint;

         //This part below is wrong!
         /* Pint *= halfdelta/counterint;
         Eint *= halfx/(halfdof+counterint);
         Cint -= Eint;
         counterint++; */
      }
   }

}

//Matlab's non-central chi square CDF up to REAL4 precision
REAL4 ncx2cdf_float(REAL4 x, REAL4 dof, REAL4 delta)
{

   REAL8 prob = 0.0;

   //Fail for bad inputs or return 0 if x<=0
   if (dof<0.0 || delta<0.0) {
      fprintf(stderr,"%s: Invalid arguments x=%f, dof=%f, delta=%f.\n", __func__, x, dof, delta);
      XLAL_ERROR_REAL4(XLAL_EINVAL);
   } else if (x<=0.0) {
      return (REAL4)prob;
   }

   REAL8 err = (REAL8)LAL_REAL4_EPS;
   REAL8 halfdelta = 0.5*delta;
   INT8 counter = (INT8)floor(halfdelta);
   REAL8 P = gsl_ran_poisson_pdf(counter, halfdelta);
   REAL8 C = twospect_cdf_chisq_P((REAL8)x, (REAL8)(dof+2.0*counter));
   if (XLAL_IS_REAL8_FAIL_NAN(C)) {
      fprintf(stderr, "%s: twospect_cdf_chisq_P(%f, %f) failed.\n", __func__, x, dof+2.0*counter);
      XLAL_ERROR_REAL4(XLAL_EFUNC);
   }
   REAL8 E = exp((dof*0.5+counter-1.0)*log(x*0.5) - x*0.5 - lgamma(dof*0.5+counter));

   sumseries_eg(&prob, P, C, E, counter, x, dof, halfdelta, err, 0);
   if (xlalErrno!=0) {
      fprintf(stderr,"%s: sumseries_eg(%f, %f, %f, %f, %f, %f, %f, %f, %f, 0) failed.\n", __func__, prob, P, C, E, floor(halfdelta), x, dof, halfdelta, err);
      XLAL_ERROR_REAL4(XLAL_EFUNC);
   }
   counter--;
   if (counter<0) return (REAL4)fmin(prob, 1.0);

   sumseries_eg(&prob, P, C, E, counter, x, dof, halfdelta, err, 1);
   if (xlalErrno!=0) {
      fprintf(stderr,"%s: sumseries_eg(%f, %f, %f, %f, %f, %f, %f, %f, %f, 1) failed.\n", __func__, prob, P, C, E, floor(halfdelta), x, dof, halfdelta, err);
      XLAL_ERROR_REAL4(XLAL_EFUNC);
   }

   //This part computes small probabilities
   INT4 fromzero = 0;
   if (prob==0.0) fromzero = 1;
   if (fromzero==1) {
      counter = 0;
      REAL8 pk = gsl_ran_poisson_pdf(0, halfdelta)*twospect_cdf_chisq_P(x, dof);
      if (xlalErrno!=0) {
         fprintf(stderr, "%s: twospect_cdf_chisq_P(%f, %f) failed.\n", __func__, x, dof);
         XLAL_ERROR_REAL4(XLAL_EFUNC);
      }
      REAL8 dp = 0.0;
      INT4 ok = 0;
      if ((REAL8)counter<halfdelta) ok = 1;
      while (ok==1) {
         counter++;
         P = gsl_ran_poisson_pdf(counter, halfdelta);
         C = twospect_cdf_chisq_P(x, dof+2.0*counter);
         if (XLAL_IS_REAL8_FAIL_NAN(C)) {
            fprintf(stderr, "%s: twospect_cdf_chisq_P(%f, %f) failed.\n", __func__, x, dof+2.0*counter);
            XLAL_ERROR_REAL4(XLAL_EFUNC);
         }
         dp = P*C;
         pk += dp;
         if (!(ok==1 && (REAL8)counter<halfdelta && dp>=err*pk)) ok = 0;
      }
      prob = pk;
   }

   return (REAL4)fmin(prob, 1.0);

}

//Matlab's non-central chi-square tries to compute very small probabilities. We don't normally need this,
//so this function leaves out the last part to compute small probabilities.
REAL8 ncx2cdf_withouttinyprob(REAL8 x, REAL8 dof, REAL8 delta)
{

   REAL8 prob = 0.0;

   //Fail for bad inputs or return 0 if x<=0
   if (dof<0.0 || delta<0.0) {
      fprintf(stderr,"%s: Invalid arguments x=%f, dof=%f, delta=%f.\n", __func__, x, dof, delta);
      XLAL_ERROR_REAL8(XLAL_EINVAL);
   } else if (x<=0.0) return prob;

   REAL8 err = LAL_REAL8_EPS;
   REAL8 halfdelta = 0.5*delta;
   INT8 counter = (INT8)floor(halfdelta);
   REAL8 P = gsl_ran_poisson_pdf(counter, halfdelta);
   REAL8 C = gsl_cdf_chisq_P(x, dof+2.0*counter);
   if (XLAL_IS_REAL8_FAIL_NAN(C)) {
      fprintf(stderr, "%s: gsl_cdf_chisq_P(%f, %f) failed.\n", __func__, x, dof+2.0*counter);
      XLAL_ERROR_REAL8(XLAL_EFUNC);
   }
   REAL8 E = exp((dof*0.5+counter-1.0)*log(x*0.5) - x*0.5 - lgamma(dof*0.5+counter));

   sumseries_eg(&prob, P, C, E, counter, x, dof, halfdelta, err, 0);
   if (xlalErrno!=0) {
      fprintf(stderr,"%s: sumseries_eg(%f, %f, %f, %f, %f, %f, %f, %f, %f, 0) failed.\n", __func__, prob, P, C, E, floor(halfdelta), x, dof, halfdelta, err);
      XLAL_ERROR_REAL8(XLAL_EFUNC);
   }
   counter--;
   if (counter<0) return fmin(prob, 1.0);

   sumseries_eg(&prob, P, C, E, counter, x, dof, halfdelta, err, 1);
   if (xlalErrno!=0) {
      fprintf(stderr,"%s: sumseries_eg(%f, %f, %f, %f, %f, %f, %f, %f, %f, 1) failed.\n", __func__, prob, P, C, E, floor(halfdelta), x, dof, halfdelta, err);
      XLAL_ERROR_REAL8(XLAL_EFUNC);
   }

   return fmin(prob, 1.0);

}

//Without small probabilities up to REAL4 precision
REAL4 ncx2cdf_float_withouttinyprob(REAL4 x, REAL4 dof, REAL4 delta)
{

   REAL8 prob = 0.0;

   //Fail for bad inputs or return 0 if x<=0
   if (dof<0.0 || delta<0.0) {
      fprintf(stderr,"%s: Invalid arguments x=%f, dof=%f, delta=%f.\n", __func__, x, dof, delta);
      XLAL_ERROR_REAL4(XLAL_EINVAL);
   } else if (x<=0.0) return (REAL4)prob;

   REAL8 err = (REAL8)LAL_REAL4_EPS;
   REAL8 halfdelta = 0.5*delta;
   INT8 counter = (INT8)floor(halfdelta);
   REAL8 P = gsl_ran_poisson_pdf(counter, halfdelta);
   REAL8 C = twospect_cdf_chisq_P((REAL8)x, (REAL8)(dof+2.0*counter));
   if (XLAL_IS_REAL8_FAIL_NAN(C)) {
      fprintf(stderr, "%s: twospect_cdf_chisq_P(%f, %f) failed.\n", __func__, x, dof+2.0*counter);
      XLAL_ERROR_REAL4(XLAL_EFUNC);
   }
   REAL8 E = exp((dof*0.5+counter-1.0)*log(x*0.5) - x*0.5 - lgamma(dof*0.5+counter));

   sumseries_eg(&prob, P, C, E, counter, x, dof, halfdelta, err, 0);
   if (xlalErrno!=0) {
      fprintf(stderr,"%s: sumseries_eg(%f, %f, %f, %f, %f, %f, %f, %f, %f, 0) failed.\n", __func__, prob, P, C, E, floor(halfdelta), x, dof, halfdelta, err);
      XLAL_ERROR_REAL4(XLAL_EFUNC);
   }
   counter--;
   if (counter<0) return (REAL4)fmin(prob, 1.0);

   sumseries_eg(&prob, P, C, E, counter, x, dof, halfdelta, err, 1);
   if (xlalErrno!=0) {
      fprintf(stderr,"%s: sumseries_eg(%f, %f, %f, %f, %f, %f, %f, %f, %f, 1) failed.\n", __func__, prob, P, C, E, floor(halfdelta), x, dof, halfdelta, err);
      XLAL_ERROR_REAL4(XLAL_EFUNC);
   }

   return (REAL4)fmin(prob, 1.0);

}


//This is ncx2cdf function like in Matlab, but using the Matlab version of the central chi square calculation instead of the GSL version
REAL8 ncx2cdf_withouttinyprob_withmatlabchi2cdf(REAL8 x, REAL8 dof, REAL8 delta)
{
   
   REAL8 prob = 0.0;
   
   //Fail for bad inputs or return 0 if x<=0
   if (dof<0.0 || delta<0.0) {
      fprintf(stderr,"%s: Invalid arguments x=%f, dof=%f, delta=%f.\n", __func__, x, dof, delta);
      XLAL_ERROR_REAL8(XLAL_EINVAL);
   } else if (x<=0.0) return prob;

   REAL8 err = LAL_REAL8_EPS;
   REAL8 halfdelta = 0.5*delta;
   INT8 counter = (INT8)floor(halfdelta);
   REAL8 P = gsl_ran_poisson_pdf(counter, halfdelta);
   REAL8 C = matlab_cdf_chisq_P(x, dof+2.0*counter);  //Matlab chi square cdf calculation
   if (XLAL_IS_REAL8_FAIL_NAN(C)) {
      fprintf(stderr, "%s: matlab_cdf_chisq_P(%f, %f) failed.\n", __func__, x, dof+2.0*counter);
      XLAL_ERROR_REAL8(XLAL_EFUNC);
   }
   REAL8 E = exp((dof*0.5+counter-1.0)*log(x*0.5) - x*0.5 - lgamma(dof*0.5+counter));

   sumseries_eg(&prob, P, C, E, counter, x, dof, halfdelta, err, 0);
   if (xlalErrno!=0) {
      fprintf(stderr,"%s: sumseries_eg(%f, %f, %f, %f, %f, %f, %f, %f, %f, 0) failed.\n", __func__, prob, P, C, E, floor(halfdelta), x, dof, halfdelta, err);
      XLAL_ERROR_REAL8(XLAL_EFUNC);
   }
   counter--;
   if (counter<0) return fmin(prob, 1.0);

   sumseries_eg(&prob, P, C, E, counter, x, dof, halfdelta, err, 1);
   if (xlalErrno!=0) {
      fprintf(stderr,"%s: sumseries_eg(%f, %f, %f, %f, %f, %f, %f, %f, %f, 1) failed.\n", __func__, prob, P, C, E, floor(halfdelta), x, dof, halfdelta, err);
      XLAL_ERROR_REAL8(XLAL_EFUNC);
   }

   return fmin(prob, 1.0);

}

//This is ncx2cdf function like in Matlab, but using the Matlab version of the central chi square calculation instead of the GSL version; up to REAL4 precision
REAL4 ncx2cdf_float_withouttinyprob_withmatlabchi2cdf(REAL4 x, REAL4 dof, REAL4 delta)
{

   REAL8 prob = 0.0;

   //Fail for bad inputs or return 0 if x<=0
   if (dof<0.0 || delta<0.0) {
      fprintf(stderr,"%s: Invalid arguments x=%f, dof=%f, delta=%f.\n", __func__, x, dof, delta);
      XLAL_ERROR_REAL8(XLAL_EINVAL);
   } else if (x<=0.0) return (REAL4)prob;
   
   REAL8 err = (REAL8)LAL_REAL4_EPS;
   REAL8 halfdelta = 0.5*delta;
   INT8 counter = (INT8)floor(halfdelta);
   REAL8 P = gsl_ran_poisson_pdf(counter, halfdelta);
   REAL8 C = matlab_cdf_chisq_P((REAL8)x, (REAL8)(dof+2.0*counter));  //Matlab chi2cdf
   if (XLAL_IS_REAL8_FAIL_NAN(C)) {
      fprintf(stderr, "%s: matlab_cdf_chisq_P(%f, %f) failed.\n", __func__, x, dof+2.0*counter);
      XLAL_ERROR_REAL4(XLAL_EFUNC);
   }
   REAL8 E = exp((dof*0.5+counter-1.0)*log(x*0.5) - x*0.5 - lgamma(dof*0.5+counter));

   sumseries_eg(&prob, P, C, E, counter, x, dof, halfdelta, err, 0);
   if (xlalErrno!=0) {
      fprintf(stderr,"%s: sumseries_eg(%f, %f, %f, %f, %f, %f, %f, %f, %f, 0) failed.\n", __func__, prob, P, C, E, floor(halfdelta), x, dof, halfdelta, err);
      XLAL_ERROR_REAL4(XLAL_EFUNC);
   }
   counter--;
   if (counter<0) return (REAL4)fmin(prob, 1.0);

   sumseries_eg(&prob, P, C, E, counter, x, dof, halfdelta, err, 1);
   if (xlalErrno!=0) {
      fprintf(stderr,"%s: sumseries_eg(%f, %f, %f, %f, %f, %f, %f, %f, %f, 1) failed.\n", __func__, prob, P, C, E, floor(halfdelta), x, dof, halfdelta, err);
      XLAL_ERROR_REAL4(XLAL_EFUNC);
   }

   return (REAL4)fmin(prob, 1.0);

}


//Like Matlabs ncx2pdf
REAL8 ncx2pdf(REAL8 x, REAL8 dof, REAL8 delta)
{

   //Fail if bad inputs
   if (XLAL_IS_REAL8_FAIL_NAN(x) || XLAL_IS_REAL8_FAIL_NAN(dof) || XLAL_IS_REAL8_FAIL_NAN(delta)) {
      fprintf(stderr,"%s: Invalid arguments x=%f, dof=%f, delta=%f.\n", __func__, x, dof, delta);
      XLAL_ERROR_REAL8(XLAL_EINVAL);
   }

   REAL8 dofint = 0.5*dof-1.0;
   REAL8 x1 = sqrt(x);
   REAL8 delta1 = sqrt(delta);
   REAL8 logreal8min = -708.3964185322641;

   REAL8 ul = 0.0;
   if (dofint<=-0.5) ul = -0.5*(delta+x) + 0.5*x1*delta1/(dofint+1.0) + dofint*(log(x)-LAL_LN2) - LAL_LN2 - lgamma(dofint+1.0);
   else ul = -0.5*(delta1-x1)*(delta1-x1) + dofint*(log(x)-LAL_LN2) - LAL_LN2 - lgamma(dofint+1.0) + (dofint+0.5)*log((dofint+0.5)/(x1*delta1+dofint+0.5));

   if (ul<logreal8min) return 0.0;

   //Scaled Bessel function?
   REAL8 sbes = gsl_sf_bessel_Inu_scaled(dofint, delta1*x1);
   if (!XLAL_IS_REAL8_FAIL_NAN(sbes) && sbes>0) {
      return exp(-LAL_LN2 - 0.5*(x1-delta1)*(x1-delta1) + dofint*log(x1/delta1))*sbes;
   }

   //Bessel function without scaling?
   REAL8 bes = gsl_sf_bessel_Inu(dofint, delta1*x1);
   if (!XLAL_IS_REAL8_FAIL_NAN(bes) && bes>0) {
      return exp(-LAL_LN2 - 0.5*(x+delta) + dofint*log(x1/delta1))*bes;
   }

   //Okay, now recursion
   REAL8 lnsr2pi = log(sqrt(LAL_TWOPI));
   REAL8 dx = delta*x*0.25;
   INT8 K = GSL_MAX_INT(0, (INT8)floor(0.5*(sqrt(dofint*dofint+4.0*dx) - dofint)));
   REAL8 lntK = 0.0;
   if (K==0) {
      lntK = -lnsr2pi - 0.5*(delta+log(dofint)) - (lgamma(dofint+1)-0.5*log(LAL_TWOPI*dofint)+dofint*log(dofint)-dofint) - binodeviance(dofint, 0.5*x);
   } else {
      lntK = -2.0*lnsr2pi - 0.5*(log(K) + log(dofint+K)) - (lgamma(K+1)-0.5*log(LAL_TWOPI*K)+K*log(K)-K) - (lgamma(dofint+K+1)-0.5*log(LAL_TWOPI*(dofint+K))+(dofint+K)*log(dofint+K)-(dofint+K)) - binodeviance(K, 0.5*delta) - binodeviance(dofint+K, 0.5*x);
   }
   REAL8 sumK = 1.0;
   INT4 keep = 0;
   if (K>0) keep = 1;
   REAL8 term = 1.0;
   REAL8 k = K;
   while (keep==1) {
      term *= (dofint+k)*k/dx;
      sumK += term;
      if (k<=0 || term<=epsval(sumK) || keep!=1) keep = 0;
      k--;
   }
   keep = 1;
   term = 1.0;
   k = K+1;
   while (keep==1) {
      term /= (dofint+k)*k/dx;
      sumK += term;
      if (term<=epsval(sumK) || keep!=1) keep = 0;
      k++;
   }
   return 0.5*exp(lntK + log(sumK));

}

//Matlab's binodeviance, a "hidden" function
REAL8 binodeviance(REAL8 x, REAL8 np)
{

   if (XLAL_IS_REAL8_FAIL_NAN(x) || XLAL_IS_REAL8_FAIL_NAN(np)) {
      fprintf(stderr,"%s: Invalid arguments x=%f, np=%f.\n", __func__, x, np);
      XLAL_ERROR_REAL8(XLAL_EINVAL);
   }

   if (fabs(x-np)<0.1*(x+np)) {
      REAL8 s = (x-np)*(x-np)/(x+np);
      REAL8 v = (x-np)/(x+np);
      REAL8 ej = 2.0*x*v;
      REAL8 s1 = 0.0;
      INT4 jj = 0;
      INT4 ok = 1;
      while (ok==1) {
         ej *= v*v;
         jj++;
         s1 = s + ej/(2.0*jj+1.0);
         if (s1!=s) {
            s = s1;
         } else {
            ok = 0;
         }
      }
      return s;
   } else {
      return x*log(x/np)+np-x;
   }

}

//Matlab's eps function for REAL8, but written in C
REAL8 epsval(REAL8 val)
{
   if (XLAL_IS_REAL8_FAIL_NAN(val)) {
      fprintf(stderr,"%s: Invalid arguments val=%f\n", __func__, val);
      XLAL_ERROR_REAL8(XLAL_EINVAL);
   }
   //Same as matlab
   REAL8 absval = fabs(val);
   int exponentval = 0;
   frexp(absval, &exponentval);
   exponentval -= LAL_REAL8_MANT;
   return ldexp(1.0, exponentval);
}

//Matlab's eps function for REAL4, but written in C
REAL4 epsval_float(REAL4 val)
{
   if (XLAL_IS_REAL4_FAIL_NAN(val)) {
      fprintf(stderr,"%s: Invalid arguments val=%f\n", __func__, val);
      XLAL_ERROR_REAL4(XLAL_EINVAL);
   }
   //Same as matlab
   REAL4 absval = fabsf(val);
   int exponentval = 0;
   frexpf(absval, &exponentval);
   exponentval -= LAL_REAL4_MANT;
   return ldexpf(1.0, exponentval);
}

//Matlab's ncx2inv() function
REAL8 ncx2inv(REAL8 p, REAL8 dof, REAL8 delta)
{

   //Fail if bad input
   if (XLAL_IS_REAL8_FAIL_NAN(p) || XLAL_IS_REAL8_FAIL_NAN(dof) || XLAL_IS_REAL8_FAIL_NAN(delta) || delta<0.0) {
      fprintf(stderr,"%s: Invalid arguments p=%f, dof=%f, delta=%f.\n", __func__, p, dof, delta);
      XLAL_ERROR_REAL8(XLAL_EINVAL);
   }

   if (delta==0.0) return gsl_cdf_chisq_Pinv(p, dof);

   REAL8 pk = p;
   INT4 count_limit = 100;
   INT4 count = 0;
   REAL8 crit = sqrt(LAL_REAL8_EPS);
   REAL8 mn = dof + delta;
   REAL8 variance = 2.0*(dof + 2.0*delta);
   REAL8 temp = log(variance + mn*mn);
   REAL8 mu = 2.0*log(mn) - 0.5*temp;
   REAL8 sigma = -2.0*log(mn) + temp;
   REAL8 xk = exp(norminv(pk, mu, sigma));
   REAL8 h = 0.0;
   REAL8 F = ncx2cdf(xk, dof, delta);
   if (XLAL_IS_REAL8_FAIL_NAN(F)) {
      fprintf(stderr, "%s: ncx2cdf(%f, %f, %f) failed.\n", __func__, xk, dof, delta);
      XLAL_ERROR_REAL8(XLAL_EFUNC);
   }
   while (count < count_limit) {
      count++;
      REAL8 f = ncx2pdf(xk, dof, delta);
      if (XLAL_IS_REAL8_FAIL_NAN(f)) {
         fprintf(stderr, "%s: ncx2pdf(%f, %f, %f) failed.\n", __func__, xk, dof, delta);
         XLAL_ERROR_REAL8(XLAL_EFUNC);
      }
      h = (F-pk)/f;
      REAL8 xnew = fmax(0.2*xk, fmin(5.0*xk, xk-h));
      REAL8 newF = ncx2cdf(xnew, dof, delta);
      if (XLAL_IS_REAL8_FAIL_NAN(newF)) {
         fprintf(stderr, "%s: ncx2cdf(%f, %f, %f) failed.\n", __func__, xnew, dof, delta);
         XLAL_ERROR_REAL8(XLAL_EFUNC);
      }
      INT4 worse = 0;
      while (worse==0) {
         if (!(fabs(newF-pk)>fabs(F-pk)*(1.0+crit) && fabs(xk-xnew)>crit*xk)) worse = 1;
         else {
            xnew = 0.5*(xnew + xk);
            newF = ncx2cdf(xnew, dof, delta);
            if (XLAL_IS_REAL8_FAIL_NAN(newF)) {
               fprintf(stderr, "%s: ncx2cdf(%f, %f, %f) failed.\n", __func__, xnew, dof, delta);
               XLAL_ERROR_REAL8(XLAL_EFUNC);
            }
         }
      }
      h = xk-xnew;
      if (!(fabs(h)>crit*fabs(xk) && fabs(h)>crit)) return xk;
      xk = xnew;
      F = newF;
   }

   fprintf(stderr, "%s: Warning! ncx2inv(%g, %g, %g) failed to converge!\n", __func__, p, dof, delta);
   return xk;

}


//Matlab's ncx2inv() function to REAL4 precision
REAL4 ncx2inv_float(REAL8 p, REAL8 dof, REAL8 delta)
{

   //Fail if bad input
   if (XLAL_IS_REAL8_FAIL_NAN(p) || XLAL_IS_REAL8_FAIL_NAN(dof) || XLAL_IS_REAL8_FAIL_NAN(delta) || delta<0.0) {
      fprintf(stderr,"%s: Invalid arguments p=%f, dof=%f, delta=%f.\n", __func__, p, dof, delta);
      XLAL_ERROR_REAL8(XLAL_EINVAL);
   }

   if (delta==0.0) return (REAL4)gsl_cdf_chisq_Pinv(p, dof);
   
   REAL8 pk = p;
   INT4 count_limit = 100;
   INT4 count = 0;
   REAL8 crit = sqrt(LAL_REAL4_EPS);
   REAL8 mn = dof + delta;
   REAL8 variance = 2.0*(dof + 2.0*delta);
   REAL8 temp = log(variance + mn*mn);
   REAL8 mu = 2.0*log(mn) - 0.5*temp;
   REAL8 sigma = -2.0*log(mn) + temp;
   REAL8 xk = exp(norminv(pk, mu, sigma));
   REAL8 h = 0.0;
   REAL8 F = ncx2cdf_float_withouttinyprob_withmatlabchi2cdf(xk, dof, delta);
   if (XLAL_IS_REAL8_FAIL_NAN(F)) {
      fprintf(stderr, "%s: ncx2cdf_float_withouttinyprob_withmatlabchi2cdf(%f, %f, %f) failed.\n", __func__, xk, dof, delta);
      XLAL_ERROR_REAL4(XLAL_EFUNC);
   }
   while (count < count_limit) {
      count++;
      REAL8 f = ncx2pdf(xk, dof, delta);
      if (XLAL_IS_REAL8_FAIL_NAN(f)) {
         fprintf(stderr, "%s: ncx2pdf(%f, %f, %f) failed.\n", __func__, xk, dof, delta);
         XLAL_ERROR_REAL4(XLAL_EFUNC);
      }
      h = (F-pk)/f;
      REAL8 xnew = fmax(0.2*xk, fmin(5.0*xk, xk-h));
      REAL8 newF = ncx2cdf_float_withouttinyprob_withmatlabchi2cdf(xnew, dof, delta);
      if (XLAL_IS_REAL8_FAIL_NAN(newF)) {
         fprintf(stderr, "%s: ncx2cdf_float_withouttinyprob_withmatlabchi2cdf(%f, %f, %f) failed.\n", __func__, xnew, dof, delta);
         XLAL_ERROR_REAL4(XLAL_EFUNC);
      }
      INT4 worse = 0;
      while (worse==0) {
         if (!(fabs(newF-pk)>fabs(F-pk)*(1.0+crit) && fabs(xk-xnew)>crit*xk)) worse = 1;
         else {
            xnew = 0.5*(xnew + xk);
            newF = ncx2cdf_float_withouttinyprob_withmatlabchi2cdf(xnew, dof, delta);
            if (XLAL_IS_REAL8_FAIL_NAN(newF)) {
               fprintf(stderr, "%s: ncx2cdf_float_withouttinyprob_withmatlabchi2cdf(%f, %f, %f) failed.\n", __func__, xnew, dof, delta);
               XLAL_ERROR_REAL4(XLAL_EFUNC);
            }
         }
      }
      h = xk-xnew;
      if (!(fabs(h)>crit*fabs(xk) && fabs(h)>crit)) return xk;
      xk = xnew;
      F = newF;
   }

   fprintf(stderr, "%s: Warning! ncx2inv_float() failed to converge!\n", __func__);
   return xk;

}


//Matlab's norminv function
REAL8 norminv(REAL8 p, REAL8 mu, REAL8 sigma)
{
   if (XLAL_IS_REAL8_FAIL_NAN(p) || XLAL_IS_REAL8_FAIL_NAN(mu) || XLAL_IS_REAL8_FAIL_NAN(sigma)) {
      fprintf(stderr,"%s: Invalid arguments p=%f, mu=%f, delta=%f.\n", __func__, p, mu, sigma);
      XLAL_ERROR_REAL8(XLAL_EINVAL);
   }
   return mu - sigma*gsl_cdf_ugaussian_Qinv(p);
}


//For the normal distribution, what is the SNR of a given value
REAL8 unitGaussianSNR(REAL8 value, REAL8 dof)
{
   REAL8 snr = (value - dof) / sqrt(2.0*dof);
   return snr;
}




/* Critical values of KS test (from Bickel and Doksum). Does not apply directly (mean determined from distribution)
 alpha=0.01
 n       10      20      30      40      50      60      80      n>80
 .489    .352    .290    .252    .226    .207    .179    1.628/(sqrt(n)+0.12+0.11/sqrt(n))
 
 alpha=0.05
 n       10      20      30      40      50      60      80      n>80
 .409    .294    .242    .210    .188    .172    .150    1.358/(sqrt(n)+0.12+0.11/sqrt(n))
 */
REAL8 ks_test_exp(REAL4Vector *vector)
{
   
   INT4 ii;
   
   //First the mean value needs to be calculated from the median value
   REAL4Vector *tempvect = XLALCreateREAL4Vector(vector->length);
   if (tempvect==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, vector->length);
      XLAL_ERROR_REAL8(XLAL_EFUNC);
   }
   memcpy(tempvect->data, vector->data, sizeof(REAL4)*vector->length);
   qsort(tempvect->data, tempvect->length, sizeof(REAL4), qsort_REAL4_compar);   //tempvect becomes sorted
   REAL4 vector_median = 0.0;
   if (tempvect->length % 2 != 1) vector_median = 0.5*(tempvect->data[(INT4)(0.5*tempvect->length)-1] + tempvect->data[(INT4)(0.5*tempvect->length)]);
   else vector_median = tempvect->data[(INT4)(0.5*tempvect->length)];
   REAL4 vector_mean = (REAL4)(vector_median/LAL_LN2);
   
   //Now start doing the K-S testing
   REAL8 ksvalue = 0.0, testval1, testval2, testval;
   REAL8 oneoverlength = 1.0/tempvect->length;
   for (ii=0; ii<(INT4)tempvect->length; ii++) {
      REAL8 pval = gsl_cdf_exponential_P(tempvect->data[ii], vector_mean);
      testval1 = fabs((1.0+ii)*oneoverlength - pval);
      testval2 = fabs(ii*oneoverlength - pval);
      testval = fmax(testval1, testval2);
      if (testval>ksvalue) ksvalue = testval;
   }
   
   //Destroy stuff
   XLALDestroyREAL4Vector(tempvect);
   
   return ksvalue;
   
}


/* Critical values of Kuiper's test using root finding by E.G.
alpha=0.05
n                                                               n>80
                                                                1.747/(sqrt(n)+0.155+0.24/sqrt(n))

alpha=0.1
n                                                               n>80
                                                                1.620/(sqrt(n)+0.155+0.24/sqrt(n)) */
REAL8 kuipers_test_exp(REAL4Vector *vector)
{
   
   INT4 ii;
   
   REAL4Vector *tempvect = XLALCreateREAL4Vector(vector->length);
   if (tempvect==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, vector->length);
      XLAL_ERROR_REAL8(XLAL_EFUNC);
   }
   
   memcpy(tempvect->data, vector->data, sizeof(REAL4)*vector->length);
   
   qsort(tempvect->data, tempvect->length, sizeof(REAL4), qsort_REAL4_compar);
   
   REAL4 vector_median = 0.0;
   if (tempvect->length % 2 != 1) vector_median = 0.5*(tempvect->data[(INT4)(0.5*tempvect->length)-1] + tempvect->data[(INT4)(0.5*tempvect->length)]);
   else vector_median = tempvect->data[(INT4)(0.5*tempvect->length)];
   
   REAL4 vector_mean = (REAL4)(vector_median/LAL_LN2);
   
   //Now the Kuiper's test calculation is made
   REAL8 loval = 0.0, hival = 0.0;
   REAL8 oneoverlength = 1.0/tempvect->length;
   /* for (ii=0; ii<(INT4)tempvect->length; ii++) {
      REAL8 pval = gsl_cdf_exponential_P(tempvect->data[ii], vector_mean);
      REAL8 testval1 = (1.0+ii)*oneoverlength - pval;
      REAL8 testval2 = ii*oneoverlength - pval;
      if (hival<testval1) hival = testval1;
      if (loval<testval2) loval = testval2;
   }
   REAL8 kuiperval1 = hival + loval; */

   loval = -1.0, hival = -1.0;
   for (ii=0; ii<(INT4)tempvect->length; ii++) {
      REAL8 pval = gsl_cdf_exponential_P(tempvect->data[ii], vector_mean);
      REAL8 testval1 = (1.0+ii)*oneoverlength - pval;
      REAL8 testval2 = ii*oneoverlength - pval;
      if (hival<testval1) hival = testval1;
      if (hival<testval2) hival = testval2;
      if (loval<-testval1) loval = -testval1;
      if (loval<-testval2) loval = -testval2;
   }
   REAL8 kuiperval2 = hival + loval;
   
   XLALDestroyREAL4Vector(tempvect);
   
   return kuiperval2;
   
}

//Sort a REAL4Vector, keeping the largest of the values in the output vector
void sort_float_largest(REAL4Vector *output, REAL4Vector *input)
{
   //Copy of the input vector
   REAL4Vector *tempvect = XLALCreateREAL4Vector(input->length);
   if (tempvect==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, input->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   memcpy(tempvect->data, input->data, sizeof(REAL4)*input->length);
   
   //qsort rearranges original vector, so sort the copy of the input vector
   qsort(tempvect->data, tempvect->length, sizeof(REAL4), qsort_REAL4_compar);
   
   INT4 ii;
   for (ii=0; ii<(INT4)output->length; ii++) output->data[ii] = tempvect->data[tempvect->length-1-ii];
   
   XLALDestroyREAL4Vector(tempvect);
   
}

//Sort a REAL4Vector, keeping the smallest of the values in the output vector
void sort_float_smallest(REAL4Vector *output, REAL4Vector *input)
{
   //Copy of the input vector
   REAL4Vector *tempvect = XLALCreateREAL4Vector(input->length);
   if (tempvect==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, input->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   memcpy(tempvect->data, input->data, sizeof(REAL4)*input->length);
   
   //qsort rearranges original vector, so sort the copy of the input vector
   qsort(tempvect->data, tempvect->length, sizeof(REAL4), qsort_REAL4_compar);
   
   memcpy(output->data, tempvect->data, sizeof(REAL4)*output->length);
   
   XLALDestroyREAL4Vector(tempvect);
   
}

//Sort a REAL8Vector from highest to lowest
/* !!!!This modifies the input vector!!!! */
void sort_double_descend(REAL8Vector *vector)
{

   INT4 ii;

   //Since the function will modify the input vector, we can do qsort() on the input vector
   qsort(vector->data, vector->length, sizeof(REAL8), qsort_REAL8_compar);

   //Make a copy of the sorted vector because we will take the reverse of the sorted vector (descending order)
   REAL8Vector *tempvect = XLALCreateREAL8Vector(vector->length);
   if (vector==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL8Vector(%d) failed.\n", __func__, vector->length);
      XLAL_ERROR_VOID(XLAL_EFUNC);
   }
   memcpy(tempvect->data, vector->data, sizeof(REAL8)*vector->length);

   //This does the descending part
   for (ii=0; ii<(INT4)vector->length; ii++) vector->data[ii] = tempvect->data[tempvect->length-1-ii];
   
   XLALDestroyREAL8Vector(tempvect);
   
}

//Sort a REAL8Vector from lowest to highest
/* !!!!This modifies the input vector!!!! */
void sort_double_ascend(REAL8Vector *vector)
{
   qsort(vector->data, vector->length, sizeof(REAL8), qsort_REAL8_compar);
}

//Sort a REAL4Vector from lowest to highest
/* !!!!This modifies the input vector!!!! */
void sort_float_ascend(REAL4Vector *vector)
{
   qsort(vector->data, vector->length, sizeof(REAL4), qsort_REAL4_compar);
}


//Sample a number (sampleSize) of values from a REAL4Vector (input) randomly
REAL4Vector * sampleREAL4Vector(REAL4Vector *input, INT4 sampleSize, gsl_rng *rng)
{

   REAL4Vector *output = XLALCreateREAL4Vector(sampleSize);
   if (output==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, sampleSize);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }
   
   INT4 ii;
   for (ii=0; ii<sampleSize; ii++) output->data[ii] = input->data[(INT4)floor(gsl_rng_uniform(rng)*input->length)];
   
   return output;
   
}

//Sample a number (sampleSize) of values from a REAL4VectorSequence (input) randomly from vector 0 up to numberofvectors
//Needs this numberofvectors limit because of the IHS algorithm
REAL4Vector * sampleREAL4VectorSequence(REAL4VectorSequence *input, INT4 numberofvectors, INT4 sampleSize, gsl_rng *rng)
{

   REAL4Vector *output = XLALCreateREAL4Vector(sampleSize);
   if (output==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, sampleSize);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }

   INT4 ii;
   for (ii=0; ii<sampleSize; ii++) output->data[ii] = input->data[(INT4)floor(gsl_rng_uniform(rng)*numberofvectors*input->vectorLength)];

   return output;

}

//Sample a number (sampleSize) of values from a REAL4VectorSequence (input) randomly from vector 0 up to numberofvectors
//Needs this numberofvectors limit because of the IHS algorithm
//This function doesn't accept zeros in the samples
REAL4Vector * sampleREAL4VectorSequence_nozerosaccepted(REAL4VectorSequence *input, INT4 numberofvectors, INT4 sampleSize, gsl_rng *rng)
{

   REAL4Vector *output = XLALCreateREAL4Vector(sampleSize);
   if (output==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, sampleSize);
      XLAL_ERROR_NULL(XLAL_EFUNC);
   }

   INT4 ii;
   for (ii=0; ii<sampleSize; ii++) {
      output->data[ii] = input->data[(INT4)floor(gsl_rng_uniform(rng)*numberofvectors*input->vectorLength)];
      while (output->data[ii]==0.0) output->data[ii] = input->data[(INT4)floor(gsl_rng_uniform(rng)*numberofvectors*input->vectorLength)];
   }

   return output;

}


//////////////////////////////////////////////////////////////
// Compute the mean value of a vector of values
REAL4 calcMean(REAL4Vector *vector)
{

   //Calculate mean from recurrance relation. Same as GSL
   INT4 ii;
   REAL8 meanval = 0.0;
   for (ii=0; ii<(INT4)vector->length; ii++) meanval += (vector->data[ii] - meanval)/(ii+1);
   
   return (REAL4)meanval;
   
} /* calcMean() */


REAL4 calcMean_ignoreZeros(REAL4Vector *vector)
{
   
   INT4 ii, values = 0;
   REAL8 meanval = 0.0;
   for (ii=0; ii<(INT4)vector->length; ii++) {
      if (vector->data[ii]!=0.0) {
         meanval += vector->data[ii];
         values++;
      }
   }

   return (REAL4)(meanval/values);

} /* calcMean_ignoreZeros() */


//////////////////////////////////////////////////////////////
// Compute the standard deviation of a vector of values
REAL4 calcStddev(REAL4Vector *vector)
{

   INT4 ii;

   double *gslarray = XLALMalloc(sizeof(double)*vector->length);
   if (gslarray==NULL) {
      fprintf(stderr,"%s: XLALMalloc(%d) failed.\n", __func__, vector->length);
      XLAL_ERROR_REAL4(XLAL_ENOMEM);
   }
   for (ii=0; ii<(INT4)vector->length; ii++) gslarray[ii] = (double)vector->data[ii];
   REAL4 stddev = (REAL4)gsl_stats_sd(gslarray, 1, vector->length);

   XLALFree((double*)gslarray);

   return stddev;

} /* calcStddev() */


REAL4 calcStddev_ignoreZeros(REAL4Vector *vector)
{

   REAL4 meanval = calcMean_ignoreZeros(vector);
   INT4 ii, values = 0;
   REAL8 sumtotal = 0.0;
   for (ii=0; ii<(INT4)vector->length; ii++) {
      if (vector->data[ii]!=0.0) {
         sumtotal += (vector->data[ii] - meanval)*(vector->data[ii] - meanval);
         values++;
      }
   }
   REAL4 stddev = sqrtf((REAL4)(sumtotal/(values-1)));

   return stddev;

} /* calcStddev() */


//////////////////////////////////////////////////////////////
// Compute the RMS of a vector of values
REAL4 calcRms(REAL4Vector *vector)
{

   REAL4Vector *sqvector = XLALCreateREAL4Vector(vector->length);
   if (sqvector==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, vector->length);
      XLAL_ERROR_REAL4(XLAL_EFUNC);
   }
   sqvector = XLALSSVectorMultiply(sqvector, vector, vector);
   if (xlalErrno!=0) {
      fprintf(stderr,"%s: XLALSSVectorMultiply() failed.\n", __func__);
      XLAL_ERROR_REAL4(XLAL_EFUNC);
   }
   REAL4 rms = sqrtf(calcMean(sqvector));

   XLALDestroyREAL4Vector(sqvector);

   return rms;

} /* calcRms() */



//////////////////////////////////////////////////////////////
// Compute the mean value of a vector of REAL8 values
REAL8 calcMeanD(REAL8Vector *vector)
{
   REAL8 meanval = gsl_stats_mean((double*)vector->data, 1, vector->length);
   return meanval;
} /* calcMeanD() */


//////////////////////////////////////////////////////////////
// Compute the standard deviation of a vector of REAL8 values
REAL8 calcStddevD(REAL8Vector *vector)
{
   REAL8 stddev = gsl_stats_sd((double*)vector->data, 1, vector->length);
   return stddev;
} /* calcStddevD() */


//Return the index value of the maximum value in a REAL4Vector
INT4 max_index(REAL4Vector *vector)
{
   
   INT4 ii = 0, indexval = 0;
   REAL4 maxval = vector->data[0];

   for (ii=1; ii<(INT4)vector->length; ii++) {
      if (vector->data[ii]>maxval) {
         maxval = vector->data[ii];
         indexval = ii;
      }
   }

   return indexval;

} /* max_index() */

//Return the index value of the maximum value in a REAL8Vector
INT4 max_index_double(REAL8Vector *vector)
{
   INT4 indexval = gsl_stats_max_index(vector->data, 1, vector->length);
   return indexval;
} /* max_index_double() */


//Return the index value of the maximum value in a REAL4Vector between startlocation and lastlocation (inclusive)
INT4 max_index_in_range(REAL4Vector *vector, INT4 startlocation, INT4 lastlocation)
{

   if (startlocation<0) {
      startlocation = 0;
   }

   INT4 ii = startlocation, indexval = ii;
   REAL4 maxval = vector->data[ii];

   if (lastlocation>=(INT4)vector->length) {
      lastlocation = (INT4)vector->length-1;
   }

   for (ii=startlocation+1; ii<=lastlocation; ii++) {
      if (vector->data[ii]>maxval) {
         maxval = vector->data[ii];
         indexval = ii;
      }
   }

   return indexval;

} /* max_index_in_range() */


//Return the index value of the maximum value from a vector (vectornum) in a REAL4VectorSequence
INT4 max_index_from_vector_in_REAL4VectorSequence(REAL4VectorSequence *vectorsequence, INT4 vectornum)
{

   INT4 ii = 0, indexval = 0;
   REAL4 maxval = vectorsequence->data[vectornum*vectorsequence->vectorLength];

   for (ii=1; ii<(INT4)vectorsequence->vectorLength; ii++) {
      if (vectorsequence->data[vectornum*vectorsequence->vectorLength+ii]>maxval) {
         maxval = vectorsequence->data[vectornum*vectorsequence->vectorLength+ii];
         indexval = ii;
      }
   }

   return indexval;

} /* max_index_from_vector_in_REAL4VectorSequence() */


//Return the index value of the maximum value and the minimum value from an INT4Vector
void min_max_index_INT4Vector(INT4Vector *inputvector, INT4 *min_index_out, INT4 *max_index_out)
{

   INT4 ii = 0;
   *min_index_out = 0, *max_index_out = 0;
   INT4 minval = inputvector->data[0];
   INT4 maxval = inputvector->data[0];

   for (ii=1; ii<(INT4)inputvector->length; ii++) {
      if (inputvector->data[ii]<minval) {
         minval = inputvector->data[ii];
         *min_index_out = ii;
      }
      if (inputvector->data[ii]>maxval) {
         maxval = inputvector->data[ii];
         *max_index_out = ii;
      }
   }

} /* min_max_index_INT4Vector() */


//Calculate the median of a REAL4Vector
REAL4 calcMedian(REAL4Vector *vector)
{
   //Make a copy of the original vector
   REAL4Vector *tempvect = XLALCreateREAL4Vector(vector->length);
   if (tempvect==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, vector->length);
      XLAL_ERROR_REAL4(XLAL_EFUNC);
   }
   memcpy(tempvect->data, vector->data, sizeof(REAL4)*vector->length);

   //qsort() on the copied data
   qsort(tempvect->data, tempvect->length, sizeof(REAL4), qsort_REAL4_compar);

   REAL4 ffdata_median = 0.0;
   if (tempvect->length % 2 != 1) ffdata_median = 0.5*(tempvect->data[(INT4)(0.5*tempvect->length)-1] + tempvect->data[(INT4)(0.5*tempvect->length)]);
   else ffdata_median = tempvect->data[(INT4)(0.5*tempvect->length)];

   XLALDestroyREAL4Vector(tempvect);

   return ffdata_median;

} /* calcMedian() */


//Calculate the median value of a REAL4Vector ignoring zero values
REAL4 calcMedian_ignoreZeros(REAL4Vector *vector)
{
   //Make a copy of the original vector
   REAL4Vector *tempvect = XLALCreateREAL4Vector(vector->length);
   if (tempvect==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", __func__, vector->length);
      XLAL_ERROR_REAL4(XLAL_EFUNC);
   }
   memcpy(tempvect->data, vector->data, sizeof(REAL4)*vector->length);
   
   //qsort() on copied data
   qsort(tempvect->data, tempvect->length, sizeof(REAL4), qsort_REAL4_compar);

   //Find forst non-zero element
   INT4 firstnonzeroelement = -1, ii = 0;
   while (firstnonzeroelement<0) {
      if (tempvect->data[ii]!=0.0) firstnonzeroelement = ii;
      else ii++;
   }

   //From the non-zero elements, find the median value
   REAL4 ffdata_median = 0.0;
   if ((tempvect->length-firstnonzeroelement) % 2 != 1) ffdata_median = 0.5*(tempvect->data[(INT4)(0.5*(tempvect->length-firstnonzeroelement))-1+firstnonzeroelement] + tempvect->data[(INT4)(0.5*(tempvect->length-firstnonzeroelement))+firstnonzeroelement]);
   else ffdata_median = tempvect->data[(INT4)(0.5*(tempvect->length-firstnonzeroelement))+firstnonzeroelement];

   XLALDestroyREAL4Vector(tempvect);

   return ffdata_median;

} /* calcMedian_ignoreZeros() */

//Comparison functions for qsort
INT4 qsort_REAL4_compar(const void *a, const void *b)
{
   const REAL4 *y = a;
   const REAL4 *z = b;

   if ( *y < *z ) return -1;
   if ( *y > *z ) return 1;
   return 0;
} /* qsort_REAL4_compar() */
INT4 qsort_REAL8_compar(const void *a, const void *b)
{
   const REAL8 *y = a;
   const REAL8 *z = b;

   if ( *y < *z ) return -1;
   if ( *y > *z ) return 1;
   return 0;
} /* qsort_REAL8_compar() */

