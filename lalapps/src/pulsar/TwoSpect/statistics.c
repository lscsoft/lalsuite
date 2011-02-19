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

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
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
   
   REAL8 noise = gsl_ran_exponential(ptrToGenerator, mu);
   
   return noise;
   
} /* expRandNum() */



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
      testval = GSL_MAX(testval1, testval2);
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
