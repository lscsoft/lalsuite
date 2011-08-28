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
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_statistics_double.h>

#include <lal/LALConstants.h>

#include "statistics.h"



//////////////////////////////////////////////////////////////
// Create a exponentially distributed noise value  -- done
REAL8 expRandNum(REAL8 mu, gsl_rng *ptrToGenerator)
{
   
   const CHAR *fn = __func__;
   
   if (mu<=0.0) {
      fprintf(stderr,"%s: expRandNum(%f, %p) failed.\n", fn, mu, ptrToGenerator);
      XLAL_ERROR_REAL8(fn, XLAL_EINVAL);
   } else if (ptrToGenerator==NULL) {
      fprintf(stderr,"%s: expRandNum(%f, %p) failed.\n", fn, mu, ptrToGenerator);
      XLAL_ERROR_REAL8(fn, XLAL_EFAULT);
   }
   
   return gsl_ran_exponential(ptrToGenerator, mu);
   
} /* expRandNum() */



//Matlab's version
REAL8 ncx2cdf(REAL8 x, REAL8 dof, REAL8 delta)
{
   
   const CHAR *fn = __func__;
   
   REAL8 prob = 0.0;
   REAL8 err = LAL_REAL8_EPS;
   REAL8 halfdelta = 0.5*delta;
   INT4 counter = (INT4)floor(halfdelta);
   REAL8 P = gsl_ran_poisson_pdf(counter, halfdelta);
   REAL8 C = gsl_cdf_chisq_P(x, dof+2.0*counter);
   REAL8 E = exp((dof*0.5+counter-1.0)*log(x*0.5) - x*0.5 - lgamma(dof*0.5+counter));
   
   sumseries(&prob, P, C, E, counter, x, dof, halfdelta, err, 0);
   if (xlalErrno!=0) {
      fprintf(stderr,"%s: sumseries() failed.\n", fn);
      XLAL_ERROR_REAL8(fn, XLAL_EFUNC);
   }
   counter--;
   if (counter<0) return fmin(prob, 1.0);
   
   sumseries(&prob, P, C, E, counter, x, dof, halfdelta, err, 1);
   if (xlalErrno!=0) {
      fprintf(stderr,"%s: sumseries() failed.\n", fn);
      XLAL_ERROR_REAL8(fn, XLAL_EFUNC);
   }
   
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
         dp = P*C;
         pk += dp;
         if (!(ok==1 && (REAL8)counter<halfdelta && dp>=err*pk)) ok = 0;
         //if ((REAL8)counter>=halfdelta || dp<err*pk || ok!=1) ok = 0;
      }
      prob = pk;
   }
   
   return fmin(prob, 1.0);
   
}
void sumseries(REAL8 *computedprob, REAL8 P, REAL8 C, REAL8 E, INT4 counter, REAL8 x, REAL8 dof, REAL8 halfdelta, REAL8 err, INT4 countdown)
{
   
   const CHAR *fn = __func__;
   
   REAL8 Pint = P, Cint = C, Eint = E;
   INT4 counterint = counter;
   INT4 j = 0;
   if (countdown!=0) {
      if (counterint>=0) j = 1;
      if (j==1) {
         Pint *= (counterint+1.0)/halfdelta;
         Cint += E;
      } else {
         counterint = -1;
      }
   }
   
   while (counterint!=-1) {
      REAL8 pplus = Pint*Cint;
      if (XLAL_IS_REAL8_FAIL_NAN(pplus)) {
         fprintf(stderr, "%s: pplus is NaN.\n", fn);
         XLAL_ERROR_VOID(fn, XLAL_EFPOVRFLW);
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
REAL4 ncx2cdf_float(REAL4 x, REAL4 dof, REAL4 delta)
{
   
   const CHAR *fn = __func__;
   
   REAL8 prob = 0.0;
   REAL8 err = (REAL8)LAL_REAL4_EPS;
   REAL8 halfdelta = 0.5*delta;
   INT4 counter = (INT4)floor(halfdelta);
   REAL8 P = gsl_ran_poisson_pdf(counter, halfdelta);
   REAL8 C = gsl_cdf_chisq_P(x, dof+2.0*counter);
   REAL8 E = exp((dof*0.5+counter-1.0)*log(x*0.5) - x*0.5 - lgamma(dof*0.5+counter));
   
   sumseries(&prob, P, C, E, counter, x, dof, halfdelta, err, 0);
   if (xlalErrno!=0) {
      fprintf(stderr,"%s: sumseries() failed.\n", fn);
      XLAL_ERROR_REAL8(fn, XLAL_EFUNC);
   }
   counter--;
   if (counter<0) return fminf(prob, 1.0);
   
   sumseries(&prob, P, C, E, counter, x, dof, halfdelta, err, 1);
   if (xlalErrno!=0) {
      fprintf(stderr,"%s: sumseries() failed.\n", fn);
      XLAL_ERROR_REAL8(fn, XLAL_EFUNC);
   }
   
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
         dp = P*C;
         pk += dp;
         if (!(ok==1 && (REAL8)counter<halfdelta && dp>=err*pk)) ok = 0;
         //if ((REAL8)counter>=halfdelta || dp<err*pk || ok!=1) ok = 0;
      }
      prob = pk;
   }
   
   return fminf(prob, 1.0);
   
}


//Like Matlabs ncx2pdf
REAL8 ncx2pdf(REAL8 x, REAL8 dof, REAL8 delta)
{
   
   REAL8 dofint = 0.5*dof-1.0;
   REAL8 x1 = sqrt(x);
   REAL8 delta1 = sqrt(delta);
   
   REAL8 logreal8min = -708.3964185322641;
   
   REAL8 ul = 0.0;
   if (dofint<=-0.5) {
      ul = -0.5*(delta+x) + 0.5*x1*delta1/(dofint+1.0) + dofint*(log(x)-LAL_LN2) - LAL_LN2 - lgamma(dofint+1.0);
   } else {
      ul = -0.5*(delta1-x1)*(delta1-x1) + dofint*(log(x)-LAL_LN2) - LAL_LN2 - lgamma(dofint+1.0) + (dofint+0.5)*log((dofint+0.5)/(x1*delta1+dofint+0.5));
   }
   if (ul<logreal8min) {
      return 0.0;
   }
   
   //Scaled Bessel function?
   REAL8 sbes = gsl_sf_bessel_Inu_scaled(dofint, delta1*x1);
   if (!XLAL_IS_REAL8_FAIL_NAN(sbes) && sbes>0) {
      return exp(-LAL_LN2 - 0.5*(x1-delta1)*(x1-delta1) + dofint*log(x1/delta1))*sbes;
   }
   
   //Bessel function without scaling?
   REAL8 bes = gsl_sf_bessel_Inu(dofint, delta1*x1);
   if (XLAL_IS_REAL8_FAIL_NAN(bes) && bes>0) {
      return exp(-LAL_LN2 - 0.5*(x+delta) + dofint*log(x1/delta1))*bes;
   }
   
   //Okay, now recursion
   REAL8 lnsr2pi = log(sqrt(LAL_TWOPI));
   REAL8 dx = delta*x*0.25;
   INT4 K = GSL_MAX_INT(0, (INT4)floor(0.5*(sqrt(dofint*dofint+4.0*dx) - dofint)));
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
REAL8 binodeviance(REAL8 x, REAL8 np)
{
   
   //From matlab's "hidden" function binodeviance
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
REAL8 epsval(REAL8 val)
{
   
   //Same as matlab
   REAL8 absval = fabs(val);
   int exponentval = 0;
   frexp(absval, &exponentval);
   exponentval -= LAL_REAL8_MANT;
   return ldexp(1.0, exponentval);
   
}

//Matlab's ncx2inv() function
REAL8 ncx2inv(REAL8 p, REAL8 dof, REAL8 delta)
{
   
   const CHAR *fn = __func__;
   
   REAL8 x = 0.0;
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
   while (count < count_limit) {
      count++;
      REAL8 f = ncx2pdf(xk, dof, delta);
      h = (F-pk)/f;
      REAL8 xnew = fmax(0.2*xk, fmin(5.0*xk, xk-h));
      REAL8 newF = ncx2cdf(xnew, dof, delta);
      INT4 worse = 0;
      while (worse==0) {
         if (!(fabs(newF-pk)>fabs(F-pk)*(1.0+crit) && fabs(xk-xnew)>crit*xk)) worse = 1;
         else {
            xnew = 0.5*(xnew + xk);
            newF = ncx2cdf(xnew, dof, delta);
         }
      }
      h = xk-xnew;
      x = xnew;
      if (!(fabs(h)>crit*fabs(xk) && fabs(h)>crit)) return xk;
      xk = xnew;
      F = newF;
   }
   
   fprintf(stderr, "%s: Warning! ncx2inv() failed to converge!\n", fn);
   return xk;
   
}
REAL8 norminv(REAL8 p, REAL8 mu, REAL8 sigma)
{
   
   return mu - sigma*gsl_cdf_ugaussian_Qinv(p);
   
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
   
   const CHAR *fn = __func__;
   
   INT4 ii;
   
   REAL4Vector *tempvect = XLALCreateREAL4Vector(vector->length);
   if (tempvect==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", fn, vector->length);
      XLAL_ERROR_REAL8(fn, XLAL_EFUNC);
   }
   
   memcpy(tempvect->data, vector->data, sizeof(REAL4)*vector->length);
   
   qsort(tempvect->data, tempvect->length, sizeof(REAL4), qsort_REAL4_compar);
   
   REAL4 vector_median = 0.0;
   if (tempvect->length % 2 != 1) vector_median = 0.5*(tempvect->data[(INT4)(0.5*tempvect->length)-1] + tempvect->data[(INT4)(0.5*tempvect->length)]);
   else vector_median = tempvect->data[(INT4)(0.5*tempvect->length)];
   
   REAL4 vector_mean = (REAL4)(vector_median/LAL_LN2);
   
   REAL8 ksvalue = 0.0, testval1, testval2, testval;
   REAL8 oneoverlength = 1.0/tempvect->length;
   for (ii=0; ii<(INT4)tempvect->length; ii++) {
      testval1 = fabs((1.0+ii)*oneoverlength - gsl_cdf_exponential_P(tempvect->data[ii], vector_mean));
      testval2 = fabs(ii*oneoverlength - gsl_cdf_exponential_P(tempvect->data[ii], vector_mean));
      testval = fmax(testval1, testval2);
      if (testval>ksvalue) ksvalue = testval;
   }
   
   XLALDestroyREAL4Vector(tempvect);
   
   return ksvalue;
   
}


void sort_float_largest(REAL4Vector *output, REAL4Vector *input)
{
   
   const CHAR *fn = __func__;
   
   REAL4Vector *tempvect = XLALCreateREAL4Vector(input->length);
   if (tempvect==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", fn, input->length);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   
   memcpy(tempvect->data, input->data, sizeof(REAL4)*input->length);
   
   qsort(tempvect->data, tempvect->length, sizeof(REAL4), qsort_REAL4_compar);
   
   INT4 ii;
   for (ii=0; ii<(INT4)output->length; ii++) output->data[ii] = tempvect->data[tempvect->length-1-ii];
   
   XLALDestroyREAL4Vector(tempvect);
   
}
void sort_float_smallest(REAL4Vector *output, REAL4Vector *input)
{
   
   const CHAR *fn = __func__;
   
   REAL4Vector *tempvect = XLALCreateREAL4Vector(input->length);
   if (tempvect==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", fn, input->length);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   
   memcpy(tempvect->data, input->data, sizeof(REAL4)*input->length);
   
   qsort(tempvect->data, tempvect->length, sizeof(REAL4), qsort_REAL4_compar);
   
   memcpy(output->data, tempvect->data, sizeof(REAL4)*output->length);
   
   XLALDestroyREAL4Vector(tempvect);
   
}

/* !!!!This modifies the input vector!!!! */
void sort_double_descend(REAL8Vector *vector)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii;
   
   qsort(vector->data, vector->length, sizeof(REAL8), qsort_REAL8_compar);
   
   REAL8Vector *tempvect = XLALCreateREAL8Vector(vector->length);
   if (vector==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL8Vector(%d) failed.\n", fn, vector->length);
      XLAL_ERROR_VOID(fn, XLAL_EFUNC);
   }
   
   memcpy(tempvect->data, vector->data, sizeof(REAL8)*vector->length);
   
   for (ii=0; ii<(INT4)vector->length; ii++) {
      vector->data[ii] = tempvect->data[tempvect->length-1-ii];
   }
   
   XLALDestroyREAL8Vector(tempvect);
   
}
/* !!!!This modifies the input vector!!!! */
void sort_double_ascend(REAL8Vector *vector)
{
   
   qsort(vector->data, vector->length, sizeof(REAL8), qsort_REAL8_compar);
   
}
/* !!!!This modifies the input vector!!!! */
void sort_float_ascend(REAL4Vector *vector)
{
   
   qsort(vector->data, vector->length, sizeof(REAL4), qsort_REAL4_compar);
   
}


REAL4Vector * sampleREAL4Vector(REAL4Vector *input, INT4 sampleSize)
{
   
   const CHAR *fn = __func__;
   
   gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
   if (rng==NULL) {
      fprintf(stderr,"%s: gsl_rng_alloc() failed.\n", fn);
      XLAL_ERROR_NULL(fn, XLAL_ENOMEM);
   }
   srand(time(NULL));
   UINT8 randseed = rand();
   gsl_rng_set(rng, randseed);
   //gsl_rng_set(rng, 0);
   
   REAL4Vector *output = XLALCreateREAL4Vector(sampleSize);
   if (output==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", fn, sampleSize);
      XLAL_ERROR_NULL(fn, XLAL_EFUNC);
   }
   
   INT4 ii;
   for (ii=0; ii<sampleSize; ii++) output->data[ii] = input->data[(INT4)floor(gsl_rng_uniform(rng)*input->length)];
   
   gsl_rng_free(rng);
   
   return output;
   
}


//////////////////////////////////////////////////////////////
// Compute the mean value of a vector of values
REAL4 calcMean(REAL4Vector *vector)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii;
   
   double *gslarray = XLALMalloc(sizeof(double)*vector->length);
   if (gslarray==NULL) {
      fprintf(stderr,"%s: XLALMalloc(%d) failed.\n", fn, vector->length);
      XLAL_ERROR_REAL4(fn, XLAL_ENOMEM);
   }
   for (ii=0; ii<(INT4)vector->length; ii++) gslarray[ii] = (double)vector->data[ii];
   REAL4 meanval = (REAL4)gsl_stats_mean(gslarray, 1, vector->length);
   
   XLALFree((double*)gslarray);
   
   return meanval;
   
} /* calcMean() */


//////////////////////////////////////////////////////////////
// Compute the standard deviation of a vector of values
REAL4 calcStddev(REAL4Vector *vector)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii;
   
   double *gslarray = XLALMalloc(sizeof(double)*vector->length);
   if (gslarray==NULL) {
      fprintf(stderr,"%s: XLALMalloc(%d) failed.\n", fn, vector->length);
      XLAL_ERROR_REAL4(fn, XLAL_ENOMEM);
   }
   for (ii=0; ii<(INT4)vector->length; ii++) gslarray[ii] = (double)vector->data[ii];
   REAL4 stddev = (REAL4)gsl_stats_sd(gslarray, 1, vector->length);
   
   XLALFree((double*)gslarray);
   
   return stddev;
   
} /* calcStddev() */



//////////////////////////////////////////////////////////////
// Compute the RMS of a vector of values
REAL4 calcRms(REAL4Vector *vector)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii;
   REAL8Vector *sqvector = XLALCreateREAL8Vector(vector->length);
   if (sqvector==NULL) {
      fprintf(stderr,"%s: XLALCreateREAL8Vector(%d) failed.\n", fn, vector->length);
      XLAL_ERROR_REAL4(fn, XLAL_EFUNC);
   }
   for (ii=0; ii<(INT4)vector->length; ii++) sqvector->data[ii] = (REAL8)(vector->data[ii]*vector->data[ii]);
   REAL4 rms = (REAL4)sqrt(calcMeanD(sqvector));
   
   /* double *gslarray = (double*)XLALMalloc(sizeof(double)*vector->length);
    for (ii=0; ii<(INT4)vector->length; ii++) gslarray[ii] = (double)vector->data[ii];
    REAL4 rms = (REAL4)sqrt(gsl_stats_tss_m(gslarray, 1, vector->length, 0.0)/vector->length); */
   
   XLALDestroyREAL8Vector(sqvector);
   //XLALFree((double*)gslarray);
   
   return rms;
   
} /* calcRms() */



//////////////////////////////////////////////////////////////
// Compute the mean value of a vector of REAL8 values
REAL8 calcMeanD(REAL8Vector *vector)
{
   
   REAL8 meanval = gsl_stats_mean((double*)vector->data, 1, vector->length);
   
   return meanval;
   
} /* calcMeanD */


//////////////////////////////////////////////////////////////
// Compute the standard deviation of a vector of REAL8 values
REAL8 calcStddevD(REAL8Vector *vector)
{
   
   REAL8 stddev = gsl_stats_sd((double*)vector->data, 1, vector->length);
   
   return stddev;
   
} /* calcStddevD */


INT4 max_index(REAL4Vector *vector)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii;
   
   double *gslarray = XLALMalloc(sizeof(double)*vector->length);
   if (gslarray==NULL) {
      fprintf(stderr,"%s: XLALMalloc(%d) failed.\n", fn, vector->length);
      XLAL_ERROR_REAL4(fn, XLAL_ENOMEM);
   }
   for (ii=0; ii<(INT4)vector->length; ii++) gslarray[ii] = (double)vector->data[ii];
   
   INT4 indexval = gsl_stats_max_index(gslarray, 1, vector->length);
   
   XLALFree((double*)gslarray);
   
   return indexval;
   
}
INT4 max_index_double(REAL8Vector *vector)
{
   
   INT4 indexval = gsl_stats_max_index(vector->data, 1, vector->length);
   
   return indexval;
   
}



REAL4 calcMedian(REAL4Vector *vector)
{
   
   const CHAR *fn = __func__;
   
   REAL4Vector *tempvect = XLALCreateREAL4Vector(vector->length);
   if (tempvect==NULL) {
      fprintf(stderr, "%s: XLALCreateREAL4Vector(%d) failed.\n", fn, vector->length);
      XLAL_ERROR(fn, XLAL_EFUNC);
   }
   
   memcpy(tempvect->data, vector->data, sizeof(REAL4)*vector->length);
   
   qsort(tempvect->data, tempvect->length, sizeof(REAL4), qsort_REAL4_compar);
   
   REAL4 ffdata_median = 0.0;
   if (tempvect->length % 2 != 1) ffdata_median = 0.5*(tempvect->data[(INT4)(0.5*tempvect->length)-1] + tempvect->data[(INT4)(0.5*tempvect->length)]);
   else ffdata_median = tempvect->data[(INT4)(0.5*tempvect->length)];
   
   XLALDestroyREAL4Vector(tempvect);
   
   return ffdata_median;
   
}


INT4 qsort_REAL4_compar(const void *a, const void *b)
{
   const REAL4 *y = a;
   const REAL4 *z = b;
   
   if ( *y < *z ) return -1;
   if ( *y > *z ) return 1;
   return 0;
   
}
INT4 qsort_REAL8_compar(const void *a, const void *b)
{
   const REAL8 *y = a;
   const REAL8 *z = b;
   
   if ( *y < *z ) return -1;
   if ( *y > *z ) return 1;
   return 0;
   
}

