/*
 *  Copyright (C) 2011, 2014, 2015 Evan Goetz
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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */

//Some functions based from Matab 2012a functions, but optimized for TwoSpect analysis

#include <math.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics_double.h>

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/VectorOps.h>

#include "statistics.h"

/**
 * Create a exponentially distributed noise value
 * \param [in] mu             Mean value of the distribution
 * \param [in] ptrToGenerator Pointer to a gsl_rng generator
 * \return A random value drawn from the exponential distribution
 */
REAL8 expRandNum(const REAL8 mu, const gsl_rng *ptrToGenerator)
{
   XLAL_CHECK_REAL8( mu > 0.0 && ptrToGenerator != NULL, XLAL_EINVAL );
   return gsl_ran_exponential(ptrToGenerator, mu);
} /* expRandNum() */

REAL4VectorAligned * expRandNumVector(const UINT4 length, const REAL8 mu, const gsl_rng *ptrToGenerator)
{
   XLAL_CHECK_NULL( mu>0.0 && ptrToGenerator!=NULL, XLAL_EINVAL );
   REAL4VectorAligned *output = NULL;
   XLAL_CHECK_NULL( (output = XLALCreateREAL4VectorAligned(length, 32)) != NULL, XLAL_EFUNC );
   for (UINT4 ii=0; ii<length; ii++) output->data[ii] = (REAL4)expRandNum(mu, ptrToGenerator);
   return output;
}

/* Critical values of KS test (from Bickel and Doksum). Does not apply directly (mean determined from distribution)
 alpha=0.01
 n       10      20      30      40      50      60      80      n>80
 .489    .352    .290    .252    .226    .207    .179    1.628/(sqrt(n)+0.12+0.11/sqrt(n))

 alpha=0.05
 n       10      20      30      40      50      60      80      n>80
 .409    .294    .242    .210    .188    .172    .150    1.358/(sqrt(n)+0.12+0.11/sqrt(n))
 */

/**
 * KS test of data against an expected exponential distribution
 * \param [out] ksvalue Pointer to the KS value
 * \param [in]  vector  Pointer to the REAL4VectorAligned to compare against an exponential distribution
 * \return Status value
 */
INT4 ks_test_exp(REAL8 *ksvalue, const REAL4VectorAligned *vector)
{

   INT4 ii;

   //First the mean value needs to be calculated from the median value
   REAL4VectorAligned *tempvect = NULL;
   XLAL_CHECK( (tempvect = XLALCreateREAL4VectorAligned(vector->length, 32)) != NULL, XLAL_EFUNC );
   memcpy(tempvect->data, vector->data, sizeof(REAL4)*vector->length);
   sort_float_ascend(tempvect);  //tempvect becomes sorted
   REAL4 vector_median = 0.0;
   if (tempvect->length % 2 != 1) vector_median = 0.5*(tempvect->data[(INT4)(0.5*tempvect->length)-1] + tempvect->data[(INT4)(0.5*tempvect->length)]);
   else vector_median = tempvect->data[(INT4)(0.5*tempvect->length)];
   REAL4 vector_mean = (REAL4)(vector_median/LAL_LN2);

   //Now start doing the K-S testing
   *ksvalue = 0.0;
   REAL8 testval1, testval2, testval;
   REAL8 oneoverlength = 1.0/tempvect->length;
   for (ii=0; ii<(INT4)tempvect->length; ii++) {
      REAL8 pval = gsl_cdf_exponential_P(tempvect->data[ii], vector_mean);
      testval1 = fabs((1.0+ii)*oneoverlength - pval);
      testval2 = fabs(ii*oneoverlength - pval);
      testval = fmax(testval1, testval2);
      if (testval>(*ksvalue)) *ksvalue = testval;
   }

   //Destroy stuff
   XLALDestroyREAL4VectorAligned(tempvect);

   return XLAL_SUCCESS;

} /* ks_test_exp() */


/* Critical values of Kuiper's test using root finding by E.G.
alpha=0.05
n                                                               n>80
                                                                1.747/(sqrt(n)+0.155+0.24/sqrt(n))

alpha=0.1
n                                                               n>80
                                                                1.620/(sqrt(n)+0.155+0.24/sqrt(n)) */

/**
 * Kuiper's test of data against an expected exponential distribution
 * \param [out] kuipervalue Pointer to the Kuiper's test value
 * \param [in]  vector      Pointer to the REAL4VectorAligned to compare against an exponential distribution
 * \return Status value
 */
INT4 kuipers_test_exp(REAL8 *kuipervalue, const REAL4VectorAligned *vector)
{

   INT4 ii;

   REAL4VectorAligned *tempvect = NULL;
   XLAL_CHECK( (tempvect = XLALCreateREAL4VectorAligned(vector->length, 32)) != NULL, XLAL_EFUNC );

   memcpy(tempvect->data, vector->data, sizeof(REAL4)*vector->length);

   sort_float_ascend(tempvect);

   REAL4 vector_median = 0.0;
   if (tempvect->length % 2 != 1) vector_median = 0.5*(tempvect->data[(INT4)(0.5*tempvect->length)-1] + tempvect->data[(INT4)(0.5*tempvect->length)]);
   else vector_median = tempvect->data[(INT4)(0.5*tempvect->length)];

   REAL4 vector_mean = (REAL4)(vector_median/LAL_LN2);

   //Now the Kuiper's test calculation is made
   REAL8 loval = 0.0, hival = 0.0;
   REAL8 oneoverlength = 1.0/tempvect->length;
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
   *kuipervalue = hival + loval;

   XLALDestroyREAL4VectorAligned(tempvect);

   return XLAL_SUCCESS;

} /* kuipers_test_exp() */


/**
 * Sort a REAL4VectorAligned, keeping the smallest of the values in the output vector
 * \param [out] output Pointer to a REAL4VectorAligned containing the output vector
 * \param [in]  input  Pointer to the REAL4VectorAligned from which to find the smallest values
 * \return Status value
 */
INT4 sort_float_smallest(REAL4VectorAligned *output, const REAL4VectorAligned *input)
{
   //Copy of the input vector
   REAL4VectorAligned *tempvect = NULL;
   XLAL_CHECK( (tempvect = XLALCreateREAL4VectorAligned(input->length, 32)) != NULL, XLAL_EFUNC );
   memcpy(tempvect->data, input->data, sizeof(REAL4)*input->length);

   //qsort rearranges original vector, so sort the copy of the input vector
   qsort(tempvect->data, tempvect->length, sizeof(REAL4), qsort_REAL4_compar);

   memcpy(output->data, tempvect->data, sizeof(REAL4)*output->length);

   XLALDestroyREAL4VectorAligned(tempvect);

   return XLAL_SUCCESS;

} /* sort_float_smallest() */


/**
 * Sort a REAL8Vector in ascending order, modifying the input vector
 * \param [in,out] vector Pointer to a REAL8Vector to be sorted
 */
void sort_double_ascend(REAL8Vector *vector)
{
   qsort(vector->data, vector->length, sizeof(REAL8), qsort_REAL8_compar);
} /* sort_double_ascend() */


/**
 * Sort a REAL4VectorAligned in ascending order, modifying the input vector
 * \param [in,out] vector Pointer to a REAL4VectorAligned to be sorted
 */
void sort_float_ascend(REAL4VectorAligned *vector)
{
   qsort(vector->data, vector->length, sizeof(REAL4), qsort_REAL4_compar);
} /* sort_float_ascend() */


/**
 * Sample a number (sampleSize) of values from a REAL4VectorAligned (input) randomly
 * \param [out] output     Pointer to output REAL4VectorAligned with length less than input
 * \param [in]  input      Pointer to a REAL4VectorAligned to be sampled from
 * \param [in]  rng        Pointer to a gsl_rng generator
 * \return Newly allocated REAL4VectorAligned of sampled values from the input vector
 */
INT4 sampleREAL4VectorAligned(REAL4VectorAligned *output, const REAL4VectorAligned *input, const gsl_rng *rng)
{
   XLAL_CHECK( output!=NULL && input!=NULL && output->length<input->length, XLAL_EINVAL );
   for (UINT4 ii=0; ii<output->length; ii++) output->data[ii] = input->data[gsl_rng_uniform_int(rng, input->length)];
   return XLAL_SUCCESS;
} /* sampleREAL4VectorAligned() */

/**
 * Sample a number (sampleSize) of values from an REAL4VectorAlignedArray (input) randomly from vector 0 up to numberofvectors without accepting any values of zero
 * Needs this numberofvectors limit because of the IHS algorithm
 * \param [in] input           Pointer to a REAL4VectorAlignedArray to be sampled from
 * \param [in] numberofvectors Number of vectors from the start from which to sample from
 * \param [in] sampleSize      Integer value for the length of the output vector
 * \param [in] rng             Pointer to a gsl_rng generator
 * \return Newly allocated REAL4VectorAligned of sampled values from the input vector
 */
REAL4VectorAligned * sampleREAL4VectorAlignedArray_nozerosaccepted(const REAL4VectorAlignedArray *input, const UINT4 numberofvectors, const UINT4 sampleSize, const gsl_rng *rng)
{
   XLAL_CHECK_NULL( input != NULL, XLAL_EINVAL );

   REAL4VectorAligned *output = NULL;
   XLAL_CHECK_NULL( (output = XLALCreateREAL4VectorAligned(sampleSize, 32)) != NULL, XLAL_EFUNC );

   for (UINT4 ii=0; ii<sampleSize; ii++) {
      output->data[ii] = input->data[(UINT4)floor(gsl_rng_uniform(rng)*numberofvectors)]->data[(UINT4)floor(gsl_rng_uniform(rng)*input->data[0]->length)];
      while (output->data[ii]==0.0) output->data[ii] = input->data[(UINT4)floor(gsl_rng_uniform(rng)*numberofvectors)]->data[(UINT4)floor(gsl_rng_uniform(rng)*input->data[0]->length)];
   }

   return output;

} /* sampleREAL4VectorAlignedArray_nozerosaccepted() */

/**
 * Compute the mean value of a REAL4VectorAligned, computed via recursion like in GSL
 * \param [in] vector Pointer to a REAL4VectorAligned of values
 * \return The mean value
 */
REAL4 calcMean(const REAL4VectorAligned *vector)
{

   //Calculate mean from recurrance relation. Same as GSL
   REAL8 meanval = 0.0;
   for (INT4 ii=0; ii<(INT4)vector->length; ii++) meanval += (vector->data[ii] - meanval)/(ii+1);

   return (REAL4)meanval;

} /* calcMean() */


/**
 * Compute the mean value of a REAL4VectorAligned without accepting values of zero
 * \param [in] vector Pointer to a REAL4VectorAligned of values
 * \return The mean value
 */
REAL4 calcMean_ignoreZeros(const REAL4VectorAligned *vector)
{

   INT4 values = 0;
   REAL8 meanval = 0.0;
   for (INT4 ii=0; ii<(INT4)vector->length; ii++) {
      if (vector->data[ii]!=0.0) {
         meanval += vector->data[ii];
         values++;
      }
   }

   if (values>0) return (REAL4)(meanval/values);
   else return 0.0;

} /* calcMean_ignoreZeros() */


/**
 * \brief Compute the harmonic mean value of a REAL4VectorAligned of SFT values
 *
 * The harmonic mean is computed from the mean values of the SFTs
 * \param [out] harmonicMean Pointer to the output REAL4 harmonic mean value
 * \param [in]  vector       Pointer to a REAL4Value from which to compute the harmonic mean
 * \param [in]  numfbins     Number of frequency bins in the SFTs
 * \param [in]  numffts      Number of SFTs during the observation time
 * \return Status value
 */
INT4 calcHarmonicMean(REAL4 *harmonicMean, const REAL4VectorAligned *vector, const UINT4 numfbins, const UINT4 numffts)
{

   UINT4 values = 0;
   REAL4VectorAligned *tempvect = NULL;
   XLAL_CHECK( (tempvect = XLALCreateREAL4VectorAligned(numfbins, 32)) != NULL, XLAL_EFUNC );

   for (UINT4 ii=0; ii<numffts; ii++) {
      if (vector->data[ii*numfbins]!=0.0) {
         memcpy(tempvect->data, &(vector->data[ii*numfbins]), sizeof(REAL4)*numfbins);
         *harmonicMean += 1.0/calcMean(tempvect);
         values++;
      }
   }
   if (values>0) *harmonicMean = (REAL4)values/(*harmonicMean);
   else *harmonicMean = 0.0;

   XLALDestroyREAL4VectorAligned(tempvect);

   return XLAL_SUCCESS;

} /* calcHarmonicMean() */


/**
 * Compute the standard deviation of a REAL4VectorAligned
 * \param [out] sigma  Pointer to the output standard deviation value
 * \param [in]  vector Pointer to a REAL4VectorAligned of values
 * \return Status value
 */
INT4 calcStddev(REAL4 *sigma, const REAL4VectorAligned *vector)
{

   double *gslarray = NULL;
   XLAL_CHECK( (gslarray = XLALMalloc(sizeof(double)*vector->length)) != NULL, XLAL_ENOMEM );
   for (INT4 ii=0; ii<(INT4)vector->length; ii++) gslarray[ii] = (double)vector->data[ii];
   *sigma = (REAL4)gsl_stats_sd(gslarray, 1, vector->length);

   XLALFree((double*)gslarray);

   return XLAL_SUCCESS;

} /* calcStddev() */


/**
 * Compute the standard deviation of a REAL4VectorAligned ignoring zero values
 * \param [out] sigma  Pointer to the output standard deviation value
 * \param [in]  vector Pointer to a REAL4VectorAligned of values
 * \return Status value
 */
INT4 calcStddev_ignoreZeros(REAL4 *sigma, const REAL4VectorAligned *vector)
{

   REAL4 meanval = calcMean_ignoreZeros(vector);
   if (meanval==0.0) {
      *sigma = 0.0;
      return XLAL_SUCCESS;
   }

   INT4 values = 0;
   REAL8 sumtotal = 0.0;
   for (INT4 ii=0; ii<(INT4)vector->length; ii++) {
      if (vector->data[ii]!=0.0) {
         sumtotal += (vector->data[ii] - meanval)*(vector->data[ii] - meanval);
         values++;
      }
   }

   if (values>1) {
      *sigma = sqrtf((REAL4)(sumtotal/(values-1)));
      return XLAL_SUCCESS;
   }
   else if (values==1) XLAL_ERROR( XLAL_EFPDIV0 );
   else XLAL_ERROR( XLAL_EFPINVAL );

} /* calcStddev_ignoreZeros() */


/**
 * Compute the RMS value of a REAL4VectorAligned
 * \param [out] rms    Pointer to the output RMS value
 * \param [in]  vector Pointer to a REAL4VectorAligned of values
 * \return Status value
 */
INT4 calcRms(REAL4 *rms, const REAL4VectorAligned *vector)
{

   REAL4VectorAligned *sqvector = NULL;
   XLAL_CHECK( (sqvector = XLALCreateREAL4VectorAligned(vector->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( XLALVectorMultiplyREAL4(sqvector->data, vector->data, vector->data, vector->length) == XLAL_SUCCESS, XLAL_EFUNC );
   *rms = sqrtf(calcMean(sqvector));

   XLALDestroyREAL4VectorAligned(sqvector);

   return XLAL_SUCCESS;

} /* calcRms() */


/**
 * Compute the mean value of a REAL8Vector
 * \param [in] vector Pointer to a REAL8Vector of values
 * \return Mean value
 */
REAL8 calcMeanD(const REAL8Vector *vector)
{
   REAL8 meanval = gsl_stats_mean((double*)vector->data, 1, vector->length);
   return meanval;
} /* calcMeanD() */


/**
 * Compute the standard deviation of a REAL8Vector
 * \param [in] vector Pointer to a REAL8Vector of values
 * \return Standard deviation value
 */
REAL8 calcStddevD(const REAL8Vector *vector)
{
   REAL8 stddev = gsl_stats_sd((double*)vector->data, 1, vector->length);
   return stddev;
} /* calcStddevD() */


/**
 * Determine the index value of the maximum value in a REAL4VectorAligned
 * \param [in] vector Pointer to REAL4VectorAligned of values
 * \return Index value of the largest element
 */
UINT4 max_index(const REAL4VectorAligned *vector)
{

   UINT4 indexval = 0;
   REAL4 maxval = vector->data[0];

   for (INT4 ii=1; ii<(INT4)vector->length; ii++) {
      if (vector->data[ii]>maxval) {
         maxval = vector->data[ii];
         indexval = ii;
      }
   }

   return indexval;

} /* max_index() */

/**
 * Determine the index value of the maximum value in a REAL8Vector
 * \param [in] vector Pointer to REAL8Vector of values
 * \return Index value of the largest element
 */
UINT4 max_index_double(const REAL8Vector *vector)
{
   UINT4 indexval = (UINT4)gsl_stats_max_index(vector->data, 1, vector->length);
   return indexval;
} /* max_index_double() */


/**
 * Determine the index value of the maximum value between elements of a REAL4VectorAligned (inclusive)
 * \param [in] vector        Pointer to REAL4VectorAligned of values
 * \param [in] startlocation Index value to start from in the REAL4VectorAligned
 * \param [in] lastlocation  Index value to end at in the REAL4VectorAligned
 * \return Index value of the largest element
 */
UINT4 max_index_in_range(const REAL4VectorAligned *vector, const UINT4 startlocation, const UINT4 lastlocation)
{
   UINT4 last = lastlocation;
   if (last>=vector->length) last = vector->length-1;

   UINT4 indexval = startlocation;
   REAL4 maxval = vector->data[startlocation];

   for (UINT4 ii=startlocation+1; ii<=last; ii++) {
      if (vector->data[ii]>maxval) {
         maxval = vector->data[ii];
         indexval = ii;
      }
   }

   return indexval;

} /* max_index_in_range() */

/**
 * Determine the index value of the maximum and minimum values in an INT4Vector
 * \param [in]  inputvector   Pointer to INT4Vector
 * \param [out] min_index_out Pointer to index value of smallest element
 * \param [out] max_index_out Pointer to index value of largest element
 * \return Status value
 */
INT4 min_max_index_INT4Vector(const INT4Vector *inputvector, UINT4 *min_index_out, UINT4 *max_index_out)
{

   *min_index_out = 0, *max_index_out = 0;
   INT4 minval = inputvector->data[0];
   INT4 maxval = inputvector->data[0];

   for (INT4 ii=1; ii<(INT4)inputvector->length; ii++) {
      if (inputvector->data[ii]<minval) {
         minval = inputvector->data[ii];
         *min_index_out = ii;
      }
      if (inputvector->data[ii]>maxval) {
         maxval = inputvector->data[ii];
         *max_index_out = ii;
      }
   }

   return XLAL_SUCCESS;

} /* min_max_index_INT4Vector() */


/**
 * Calculate the median value from a REAL4VectorAligned
 * \param [out] median Pointer to the median value
 * \param [in]  vector Pointer to a REAL4VectorAligned
 * \return Status value
 */
INT4 calcMedian(REAL4 *median, const REAL4VectorAligned *vector)
{
   //Make a copy of the original vector
   REAL4VectorAligned *tempvect = NULL;
   XLAL_CHECK( (tempvect = XLALCreateREAL4VectorAligned(vector->length, 32)) != NULL, XLAL_EFUNC );
   memcpy(tempvect->data, vector->data, sizeof(REAL4)*vector->length);

   //qsort() on the copied data
   qsort(tempvect->data, tempvect->length, sizeof(REAL4), qsort_REAL4_compar);

   if (tempvect->length % 2 != 1) *median = 0.5*(tempvect->data[(INT4)(0.5*tempvect->length)-1] + tempvect->data[(INT4)(0.5*tempvect->length)]);
   else *median = tempvect->data[(INT4)(0.5*tempvect->length)];

   XLALDestroyREAL4VectorAligned(tempvect);

   return XLAL_SUCCESS;

} /* calcMedian() */


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
