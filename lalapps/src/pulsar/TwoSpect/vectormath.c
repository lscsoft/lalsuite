/*
 *  Copyright (C) 2011, 2012 Evan Goetz
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



#ifdef __SSE__
#include <xmmintrin.h>
#endif
#ifdef __SSE2__
#include <emmintrin.h>
#endif

#include <math.h>

#include <lal/LALConstants.h>

#include "vectormath.h"
#include "templates.h"


#define OOTWOPI         (1.0 / LAL_TWOPI)
#define LUT_RES         1024                 /* resolution of lookup-table */
#define LUT_RES_F       (1.0 * LUT_RES)
#define OO_LUT_RES      (1.0 / LUT_RES)
#define X_TO_IND        (1.0 * LUT_RES * OOTWOPI )
#define IND_TO_X        (LAL_TWOPI * OO_LUT_RES)
#define TRUE            (1==1)
#define FALSE           (1==0)


//Computes a multiplication of two vectors with a stride and initial offset
//Be sure you know what you are doing or this could go wrong (no error checking for speed!)
REAL4Vector * fastSSVectorMultiply_with_stride_and_offset(REAL4Vector *output, REAL4Vector *input1, REAL4Vector *input2, INT4 stride1, INT4 stride2, INT4 offset1, INT4 offset2)
{
   
   REAL4 *a, *b, *c;
   INT4   n;
   
   a = input1->data + offset1;
   b = input2->data + offset2;
   c = output->data;
   n = output->length;
   
   while (n-- > 0) {
      *c = (*a)*(*b);
      a = a + stride1;
      b = b + stride2;
      c++;
   }
   
   return output;
   
} /* SSVectorMultiply_with_stride_and_offset() */


//Sums a sequence of vector values within two vector sequences
//vectorpos1 = vector number of first vector
//vectorpos2 = vector number of second vector
//outputvectorpos = vector number to load the sum into
void fastSSVectorSequenceSum(REAL4VectorSequence *output, REAL4VectorSequence *input1, REAL4VectorSequence *input2, INT4 vectorpos1, INT4 vectorpos2, INT4 outputvectorpos)
{
   
   REAL4 *a, *b, *c;
   INT4 n;
   
   a = &(input1->data[vectorpos1*input1->vectorLength]);
   b = &(input2->data[vectorpos2*input2->vectorLength]);
   c = &(output->data[outputvectorpos*output->vectorLength]);
   n = output->vectorLength;
   
   while (n-- > 0) {
      *c = (*a)+(*b);
      a++;
      b++;
      c++;
   }
   
}

//Does a fast subtraction of 1 vector for a specific vector in a vector sequence (labeled by vectorpos1)
//Output is a single (REAL4) vector
void fastSSVectorSequenceSubtract(REAL4Vector *output, REAL4VectorSequence *input1, REAL4Vector *input2, INT4 vectorpos1)
{
   
   REAL4 *a, *b, *c;
   INT4 n;
   
   a = &(input1->data[vectorpos1*input1->vectorLength]);
   b = input2->data;
   c = output->data;
   n = output->length;
   
   while (n-- > 0) {
      *c = (*a)-(*b);
      a++;
      b++;
      c++;
   }
   
}



//Sum two REAL4Vectors using SSE
REAL4Vector * sseSSVectorSum(REAL4Vector *output, REAL4Vector *input1, REAL4Vector *input2)
{
   
#ifdef __SSE__
   INT4 roundedvectorlength = (INT4)input1->length / 4;
   INT4 vec1aligned = 0, vec2aligned = 0, outputaligned = 0, ii = 0;
   
   REAL4 *allocinput1 = NULL, *allocinput2 = NULL, *allocoutput = NULL, *alignedinput1 = NULL, *alignedinput2 = NULL, *alignedoutput = NULL;
   __m128 *arr1, *arr2, *result;
   
   //Allocate memory for aligning input vector 1 if necessary
   if ( input1->data==(void*)(((UINT8)input1->data+15) & ~15) ) {
      vec1aligned = 1;
      arr1 = (__m128*)(void*)input1->data;
   } else {
      allocinput1 = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15);
      if (allocinput1==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 4*roundedvectorlength*sizeof(REAL4) + 15);
         XLAL_ERROR_NULL(XLAL_ENOMEM);
      }
      alignedinput1 = (void*)(((UINT8)allocinput1+15) & ~15);
      memcpy(alignedinput1, input1->data, sizeof(REAL4)*4*roundedvectorlength);
      arr1 = (__m128*)(void*)alignedinput1;
   }
   
   //Allocate memory for aligning input vector 2 if necessary
   if ( input2->data==(void*)(((UINT8)input2->data+15) & ~15) ) {
      vec2aligned = 1;
      arr2 = (__m128*)(void*)input2->data;
   } else {
      allocinput2 = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15);
      if (allocinput2==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 4*roundedvectorlength*sizeof(REAL4) + 15);
         XLAL_ERROR_NULL(XLAL_ENOMEM);
      }
      alignedinput2 = (void*)(((UINT8)allocinput2+15) & ~15);
      memcpy(alignedinput2, input2->data, sizeof(REAL4)*4*roundedvectorlength);
      arr2 = (__m128*)(void*)alignedinput2;
   }
   
   //Allocate memory for aligning output vector if necessary
   if ( output->data==(void*)(((UINT8)output->data+15) & ~15) ) {
      outputaligned = 1;
      result = (__m128*)(void*)output->data;
   } else {
      allocoutput = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15);
      if (allocoutput==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 4*roundedvectorlength*sizeof(REAL4) + 15);
         XLAL_ERROR_NULL(XLAL_ENOMEM);
      }
      alignedoutput = (void*)(((UINT8)allocoutput+15) & ~15);
      result = (__m128*)(void*)alignedoutput;
   }
   
   //multiply the two vectors into the output
   for (ii=0; ii<roundedvectorlength; ii++) {
      *result = _mm_add_ps(*arr1, *arr2);
      arr1++;
      arr2++;
      result++;
   }
   
   //Copy output aligned memory to non-aligned memory if necessary
   if (!outputaligned) memcpy(output->data, alignedoutput, 4*roundedvectorlength*sizeof(REAL4));
   
   //Finish up the remaining part
   for (ii=4*roundedvectorlength; ii<(INT4)input1->length; ii++) output->data[ii] = input1->data[ii] + input2->data[ii];
   
   //Free memory if necessary
   if (!vec1aligned) XLALFree(allocinput1);
   if (!vec2aligned) XLALFree(allocinput2);
   if (!outputaligned) XLALFree(allocoutput);
   
   return output;
#else
   fprintf(stderr, "%s: Failed because SSE is not supported, possibly because -msse flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR_NULL(XLAL_EFAILED);
#endif
   
}



//Multiply two REAL4Vectors using SSE
REAL4Vector * sseSSVectorMultiply(REAL4Vector *output, REAL4Vector *input1, REAL4Vector *input2)
{
   
#ifdef __SSE__
   INT4 roundedvectorlength = (INT4)input1->length / 4;
   INT4 vec1aligned = 0, vec2aligned = 0, outputaligned = 0, ii = 0;
   
   REAL4 *allocinput1 = NULL, *allocinput2 = NULL, *allocoutput = NULL, *alignedinput1 = NULL, *alignedinput2 = NULL, *alignedoutput = NULL;
   __m128 *arr1, *arr2, *result;
   
   //Allocate memory for aligning input vector 1 if necessary
   if ( input1->data==(void*)(((UINT8)input1->data+15) & ~15) ) {
      vec1aligned = 1;
      arr1 = (__m128*)(void*)input1->data;
   } else {
      allocinput1 = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15);
      if (allocinput1==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 4*roundedvectorlength*sizeof(REAL4) + 15);
         XLAL_ERROR_NULL(XLAL_ENOMEM);
      }
      alignedinput1 = (void*)(((UINT8)allocinput1+15) & ~15);
      memcpy(alignedinput1, input1->data, sizeof(REAL4)*4*roundedvectorlength);
      arr1 = (__m128*)(void*)alignedinput1;
   }
   
   //Allocate memory for aligning input vector 2 if necessary
   if ( input2->data==(void*)(((UINT8)input2->data+15) & ~15) ) {
      vec2aligned = 1;
      arr2 = (__m128*)(void*)input2->data;
   } else {
      allocinput2 = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15);
      if (allocinput2==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 4*roundedvectorlength*sizeof(REAL4) + 15);
         XLAL_ERROR_NULL(XLAL_ENOMEM);
      }
      alignedinput2 = (void*)(((UINT8)allocinput2+15) & ~15);
      memcpy(alignedinput2, input2->data, sizeof(REAL4)*4*roundedvectorlength);
      arr2 = (__m128*)(void*)alignedinput2;
   }
   
   //Allocate memory for aligning output vector if necessary
   if ( output->data==(void*)(((UINT8)output->data+15) & ~15) ) {
      outputaligned = 1;
      result = (__m128*)(void*)output->data;
   } else {
      allocoutput = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15);
      if (allocoutput==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 4*roundedvectorlength*sizeof(REAL4) + 15);
         XLAL_ERROR_NULL(XLAL_ENOMEM);
      }
      alignedoutput = (void*)(((UINT8)allocoutput+15) & ~15);
      result = (__m128*)(void*)alignedoutput;
   }
   
   //multiply the two vectors into the output
   for (ii=0; ii<roundedvectorlength; ii++) {
      *result = _mm_mul_ps(*arr1, *arr2);
      arr1++;
      arr2++;
      result++;
   }
   
   //Copy output aligned memory to non-aligned memory if necessary
   if (!outputaligned) memcpy(output->data, alignedoutput, 4*roundedvectorlength*sizeof(REAL4));
   
   //Finish up the remaining part
   for (ii=4*roundedvectorlength; ii<(INT4)input1->length; ii++) output->data[ii] = input1->data[ii] * input2->data[ii];
   
   //Free memory if necessary
   if (!vec1aligned) XLALFree(allocinput1);
   if (!vec2aligned) XLALFree(allocinput2);
   if (!outputaligned) XLALFree(allocoutput);
   
   return output;
#else
   fprintf(stderr, "%s: Failed because SSE is not supported, possibly because -msse flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR_NULL(XLAL_EFAILED);
#endif
   
}


//Add a REAL4 to all REAL4Vector elements using SSE
REAL4Vector * sseAddScalarToREAL4Vector(REAL4Vector *output, REAL4Vector *input, REAL4 scalar)
{
   
#ifdef __SSE__
   INT4 roundedvectorlength = (INT4)input->length / 4;
   INT4 vecaligned = 0, outputaligned = 0, ii = 0;
   
   REAL4 *allocinput = NULL, *allocoutput = NULL, *alignedinput = NULL, *alignedoutput = NULL;
   __m128 *arr1, *result;
   
   __m128 scalefactor = _mm_set1_ps(scalar);
   
   //Allocate memory for aligning input vector 1 if necessary
   if ( input->data==(void*)(((UINT8)input->data+15) & ~15) ) {
      vecaligned = 1;
      arr1 = (__m128*)(void*)input->data;
   } else {
      allocinput = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15);
      if (allocinput==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 4*roundedvectorlength*sizeof(REAL4) + 15);
         XLAL_ERROR_NULL(XLAL_ENOMEM);
      }
      alignedinput = (void*)(((UINT8)allocinput+15) & ~15);
      memcpy(alignedinput, input->data, sizeof(REAL4)*4*roundedvectorlength);
      arr1 = (__m128*)(void*)alignedinput;
   }
   
   //Allocate memory for aligning output vector if necessary
   if ( output->data==(void*)(((UINT8)output->data+15) & ~15) ) {
      outputaligned = 1;
      result = (__m128*)(void*)output->data;
   } else {
      allocoutput = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15);
      if (allocoutput==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 4*roundedvectorlength*sizeof(REAL4) + 15);
         XLAL_ERROR_NULL(XLAL_ENOMEM);
      }
      alignedoutput = (void*)(((UINT8)allocoutput+15) & ~15);
      result = (__m128*)(void*)alignedoutput;
   }
   
   //Add the value to the vector and put in the output
   for (ii=0; ii<roundedvectorlength; ii++) {
      *result = _mm_add_ps(*arr1, scalefactor);
      arr1++;
      result++;
   }
   
   //Copy output aligned memory to non-aligned memory if necessary
   if (!outputaligned) memcpy(output->data, alignedoutput, 4*roundedvectorlength*sizeof(REAL4));
   
   //Finish up the remaining part
   for (ii=4*roundedvectorlength; ii<(INT4)input->length; ii++) output->data[ii] = input->data[ii] + scalar;
   
   //Free memory if necessary
   if (!vecaligned) XLALFree(allocinput);
   if (!outputaligned) XLALFree(allocoutput);
   
   return output;
#else
   fprintf(stderr, "%s: Failed because SSE is not supported, possibly because -msse flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR_NULL(XLAL_EFAILED);
#endif
   
}


//Add a REAL8 to all REAL8Vector elements using SSE
REAL8Vector * sseAddScalarToREAL8Vector(REAL8Vector *output, REAL8Vector *input, REAL8 scalar)
{
   
#ifdef __SSE2__
   INT4 roundedvectorlength = (INT4)input->length / 2;
   INT4 vecaligned = 0, outputaligned = 0, ii = 0;
   
   REAL8 *allocinput = NULL, *allocoutput = NULL, *alignedinput = NULL, *alignedoutput = NULL;
   __m128d *arr1, *result;
   
   __m128d scalefactor = _mm_set1_pd(scalar);
   
   //Allocate memory for aligning input vector 1 if necessary
   if ( input->data==(void*)(((UINT8)input->data+15) & ~15) ) {
      vecaligned = 1;
      arr1 = (__m128d*)(void*)input->data;
   } else {
      allocinput = (REAL8*)XLALMalloc(2*roundedvectorlength*sizeof(REAL8) + 15);
      if (allocinput==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 2*roundedvectorlength*sizeof(REAL8) + 15);
         XLAL_ERROR_NULL(XLAL_ENOMEM);
      }
      alignedinput = (void*)(((UINT8)allocinput+15) & ~15);
      memcpy(alignedinput, input->data, sizeof(REAL8)*2*roundedvectorlength);
      arr1 = (__m128d*)(void*)alignedinput;
   }
   
   //Allocate memory for aligning output vector if necessary
   if ( output->data==(void*)(((UINT8)output->data+15) & ~15) ) {
      outputaligned = 1;
      result = (__m128d*)(void*)output->data;
   } else {
      allocoutput = (REAL8*)XLALMalloc(2*roundedvectorlength*sizeof(REAL8) + 15);
      if (allocoutput==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 2*roundedvectorlength*sizeof(REAL8) + 15);
         XLAL_ERROR_NULL(XLAL_ENOMEM);
      }
      alignedoutput = (void*)(((UINT8)allocoutput+15) & ~15);
      result = (__m128d*)(void*)alignedoutput;
   }
   
   //Add the value to the vector and put in the output
   for (ii=0; ii<roundedvectorlength; ii++) {
      *result = _mm_add_pd(*arr1, scalefactor);
      arr1++;
      result++;
   }
   
   //Copy output aligned memory to non-aligned memory if necessary
   if (!outputaligned) memcpy(output->data, alignedoutput, 2*roundedvectorlength*sizeof(REAL8));
   
   //Finish up the remaining part
   for (ii=2*roundedvectorlength; ii<(INT4)input->length; ii++) output->data[ii] = input->data[ii] + scalar;
   
   //Free memory if necessary
   if (!vecaligned) XLALFree(allocinput);
   if (!outputaligned) XLALFree(allocoutput);
   
   return output;
#else
   fprintf(stderr, "%s: Failed because SSE2 is not supported, possibly because -msse2 flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR_NULL(XLAL_EFAILED);
#endif
   
}


//Scale a REAL4Vector with a single scale factor using SSE
REAL4Vector * sseScaleREAL4Vector(REAL4Vector *output, REAL4Vector *input, REAL4 scale)
{
   
#ifdef __SSE__
   INT4 roundedvectorlength = (INT4)input->length / 4;
   INT4 vecaligned = 0, outputaligned = 0, ii = 0;
   
   REAL4 *allocinput = NULL, *allocoutput = NULL, *alignedinput = NULL, *alignedoutput = NULL;
   __m128 *arr1, *result;
   
   __m128 scalefactor = _mm_set1_ps(scale);
   
   //Allocate memory for aligning input vector 1 if necessary
   if ( input->data==(void*)(((UINT8)input->data+15) & ~15) ) {
      vecaligned = 1;
      arr1 = (__m128*)(void*)input->data;
   } else {
      allocinput = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15);
      if (allocinput==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 4*roundedvectorlength*sizeof(REAL4) + 15);
         XLAL_ERROR_NULL(XLAL_ENOMEM);
      }
      alignedinput = (void*)(((UINT8)allocinput+15) & ~15);
      memcpy(alignedinput, input->data, sizeof(REAL4)*4*roundedvectorlength);
      arr1 = (__m128*)(void*)alignedinput;
   }
   
   //Allocate memory for aligning output vector if necessary
   if ( output->data==(void*)(((UINT8)output->data+15) & ~15) ) {
      outputaligned = 1;
      result = (__m128*)(void*)output->data;
   } else {
      allocoutput = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15);
      if (allocoutput==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 4*roundedvectorlength*sizeof(REAL4) + 15);
         XLAL_ERROR_NULL(XLAL_ENOMEM);
      }
      alignedoutput = (void*)(((UINT8)allocoutput+15) & ~15);
      result = (__m128*)(void*)alignedoutput;
   }
   
   //multiply the vector into the output
   for (ii=0; ii<roundedvectorlength; ii++) {
      *result = _mm_mul_ps(*arr1, scalefactor);
      arr1++;
      result++;
   }
   
   //Copy output aligned memory to non-aligned memory if necessary
   if (!outputaligned) memcpy(output->data, alignedoutput, 4*roundedvectorlength*sizeof(REAL4));
   
   //Finish up the remaining part
   for (ii=4*roundedvectorlength; ii<(INT4)input->length; ii++) output->data[ii] = input->data[ii] * scale;
   
   //Free memory if necessary
   if (!vecaligned) XLALFree(allocinput);
   if (!outputaligned) XLALFree(allocoutput);
   
   return output;
#else
   fprintf(stderr, "%s: Failed because SSE is not supported, possibly because -msse flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR_NULL(XLAL_EFAILED);
#endif
   
}


//Scale a REAL8Vector with a single scale factor using SSE
REAL8Vector * sseScaleREAL8Vector(REAL8Vector *output, REAL8Vector *input, REAL8 scale)
{
   
#ifdef __SSE2__
   INT4 roundedvectorlength = (INT4)input->length / 2;
   INT4 vecaligned = 0, outputaligned = 0, ii = 0;
   
   REAL8 *allocinput = NULL, *allocoutput = NULL, *alignedinput = NULL, *alignedoutput = NULL;
   __m128d *arr1, *result;
   
   __m128d scalefactor = _mm_set1_pd(scale);
   
   //Allocate memory for aligning input vector 1 if necessary
   if ( input->data==(void*)(((UINT8)input->data+15) & ~15) ) {
      vecaligned = 1;
      arr1 = (__m128d*)(void*)input->data;
   } else {
      allocinput = (REAL8*)XLALMalloc(2*roundedvectorlength*sizeof(REAL8) + 15);
      if (allocinput==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 2*roundedvectorlength*sizeof(REAL8) + 15);
         XLAL_ERROR_NULL(XLAL_ENOMEM);
      }
      alignedinput = (void*)(((UINT8)allocinput+15) & ~15);
      memcpy(alignedinput, input->data, sizeof(REAL8)*2*roundedvectorlength);
      arr1 = (__m128d*)(void*)alignedinput;
   }
   
   //Allocate memory for aligning output vector if necessary
   if ( output->data==(void*)(((UINT8)output->data+15) & ~15) ) {
      outputaligned = 1;
      result = (__m128d*)(void*)output->data;
   } else {
      allocoutput = (REAL8*)XLALMalloc(2*roundedvectorlength*sizeof(REAL8) + 15);
      if (allocoutput==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 2*roundedvectorlength*sizeof(REAL8) + 15);
         XLAL_ERROR_NULL(XLAL_ENOMEM);
      }
      alignedoutput = (void*)(((UINT8)allocoutput+15) & ~15);
      result = (__m128d*)(void*)alignedoutput;
   }
   
   //multiply the vector into the output
   for (ii=0; ii<roundedvectorlength; ii++) {
      *result = _mm_mul_pd(*arr1, scalefactor);
      arr1++;
      result++;
   }
   
   //Copy output aligned memory to non-aligned memory if necessary
   if (!outputaligned) memcpy(output->data, alignedoutput, 2*roundedvectorlength*sizeof(REAL8));
   
   //Finish up the remaining part
   for (ii=2*roundedvectorlength; ii<(INT4)input->length; ii++) output->data[ii] = input->data[ii] * scale;
   
   //Free memory if necessary
   if (!vecaligned) XLALFree(allocinput);
   if (!outputaligned) XLALFree(allocoutput);
   
   return output;
#else
   fprintf(stderr, "%s: Failed because SSE2 is not supported, possibly because -msse2 flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR_NULL(XLAL_EFAILED);
#endif
   
}



//Using SSE, sum up a sequence of values from two vectors in a REAL4VectorSequence
//vectorpos1 = vector number of first vector
//vectorpos2 = vector number of second vector
//outputvectorpos = vector number to load the sum into
//numvectors = number of times to repeat the process, incrementing the vector positions by 1 each time
void sseSSVectorSequenceSum(REAL4VectorSequence *output, REAL4VectorSequence *input1, REAL4VectorSequence *input2, INT4 vectorpos1, INT4 vectorpos2, INT4 outputvectorpos, INT4 numvectors)
{
   
#ifdef __SSE__
   INT4 roundedvectorlength = (INT4)input1->vectorLength / 4;
   
   REAL4* allocinput1 = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15);
   REAL4* allocinput2 = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15);
   REAL4* allocoutput = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15);
   if (allocinput1==NULL) {
      fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 4*roundedvectorlength*sizeof(REAL4) + 15);
      XLAL_ERROR_VOID(XLAL_ENOMEM);
   } else if (allocinput2==NULL) {
      fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 4*roundedvectorlength*sizeof(REAL4) + 15);
      XLAL_ERROR_VOID(XLAL_ENOMEM);
   } else if (allocoutput==NULL) {
      fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 4*roundedvectorlength*sizeof(REAL4) + 15);
      XLAL_ERROR_VOID(XLAL_ENOMEM);
   }
   REAL4* alignedinput1 = (void*)(((UINT8)allocinput1+15) & ~15);
   REAL4* alignedinput2 = (void*)(((UINT8)allocinput2+15) & ~15);
   REAL4* alignedoutput = (void*)(((UINT8)allocoutput+15) & ~15);
   
   INT4 ii, jj;
   for (ii=0; ii<numvectors; ii++) {
      INT4 vec1 = (vectorpos1+ii)*input1->vectorLength, vec2 = (vectorpos2+ii)*input2->vectorLength, outvec = (outputvectorpos+ii)*output->vectorLength;
      
      INT4 outputvecaligned = 0;
      __m128 *arr1, *arr2, *result;
      if ( &(input1->data[vec1])==(void*)(((UINT8)&(input1->data[vec1])+15) & ~15) ) {
         arr1 = (__m128*)(void*)&(input1->data[vec1]);
      } else {
         memcpy(alignedinput1, &(input1->data[vec1]), sizeof(REAL4)*4*roundedvectorlength);
         arr1 = (__m128*)(void*)alignedinput1;
      }
      if ( &(input2->data[vec2])==(void*)(((UINT8)&(input2->data[vec2])+15) & ~15) ) {
         arr2 = (__m128*)(void*)&(input2->data[vec2]);
      } else {
         memcpy(alignedinput2, &(input2->data[vec2]), sizeof(REAL4)*4*roundedvectorlength);
         arr2 = (__m128*)(void*)alignedinput2;
      }
      if ( &(output->data[outvec])==(void*)(((UINT8)&(output->data[outvec])+15) & ~15) ) {
         outputvecaligned = 1;
         result = (__m128*)(void*)&(output->data[outvec]);
      } else result = (__m128*)(void*)alignedoutput;
      
      for (jj=0; jj<roundedvectorlength; jj++) {
         *result = _mm_add_ps(*arr1, *arr2);
         arr1++;
         arr2++;
         result++;
      }
      
      if (!outputvecaligned) memcpy(&(output->data[outvec]), alignedoutput, sizeof(REAL4)*4*roundedvectorlength);
      
      REAL4 *a = &(input1->data[vec1+4*roundedvectorlength]);
      REAL4 *b = &(input2->data[vec2+4*roundedvectorlength]);
      REAL4 *c = &(output->data[outvec+4*roundedvectorlength]);
      INT4 n = output->vectorLength-4*roundedvectorlength;
      while (n-- > 0) {
         *c = (*a)+(*b);
         a++;
         b++;
         c++;
      }
   }
   
   XLALFree(allocinput1);
   XLALFree(allocinput2);
   XLALFree(allocoutput);
#else
   fprintf(stderr, "%s: Failed because SSE is not supported, possibly because -msse flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR_VOID(XLAL_EFAILED);
#endif
   
}


//Does a fast subtraction of 1 vector for a specific vector in a vector sequence (labeled by vectorpos1) using SSE
//Output is a single (REAL4) vector
void sseSSVectorSequenceSubtract(REAL4Vector *output, REAL4VectorSequence *input1, REAL4Vector *input2, INT4 vectorpos1)
{
   
#ifdef __SSE__
   INT4 roundedvectorlength = (INT4)input1->vectorLength / 4;
   INT4 vec1 = vectorpos1*input1->vectorLength, ii = 0;
   INT4 vec1aligned = 0, vec2aligned = 0, outputaligned = 0;
   REAL4 *allocinput1 = NULL, *allocinput2 = NULL, *allocoutput = NULL, *alignedinput1 = NULL, *alignedinput2 = NULL, *alignedoutput = NULL;
   __m128 *arr1, *arr2, *result;
   
   //Allocate memory for aligning input vector 1 if necessary
   if ( &(input1->data[vec1])==(void*)(((UINT8)&(input1->data[vec1])+15) & ~15) ) {
      vec1aligned = 1;
      arr1 = (__m128*)(void*)&(input1->data[vec1]);
   } else {
      allocinput1 = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15);
      if (allocinput1==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 4*roundedvectorlength*sizeof(REAL4) + 15);
         XLAL_ERROR_VOID(XLAL_ENOMEM);
      }
      alignedinput1 = (void*)(((UINT8)allocinput1+15) & ~15);
      memcpy(alignedinput1, &(input1->data[vec1]), sizeof(REAL4)*4*roundedvectorlength);
      arr1 = (__m128*)(void*)alignedinput1;
   }
   
   //Allocate memory for aligning input vector 2 if necessary
   if ( input2->data==(void*)(((UINT8)input2->data+15) & ~15) ) {
      vec2aligned = 1;
      arr2 = (__m128*)(void*)input2->data;
   } else {
      allocinput2 = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15);
      if (allocinput2==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 4*roundedvectorlength*sizeof(REAL4) + 15);
         XLAL_ERROR_VOID(XLAL_ENOMEM);
      }
      alignedinput2 = (void*)(((UINT8)allocinput2+15) & ~15);
      memcpy(alignedinput2, input2->data, sizeof(REAL4)*4*roundedvectorlength);
      arr2 = (__m128*)(void*)alignedinput2;
   }
   
   //Allocate memory for aligning output vector if necessary
   if ( output->data==(void*)(((UINT8)output->data+15) & ~15) ) {
      outputaligned = 1;
      result = (__m128*)(void*)output->data;
   } else {
      allocoutput = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15);
      if (allocoutput==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 4*roundedvectorlength*sizeof(REAL4) + 15);
         XLAL_ERROR_VOID(XLAL_ENOMEM);
      }
      alignedoutput = (void*)(((UINT8)allocoutput+15) & ~15);
      result = (__m128*)(void*)alignedoutput;
   }
   
   for (ii=0; ii<roundedvectorlength; ii++) {
      *result = _mm_sub_ps(*arr1, *arr2);
      arr1++;
      arr2++;
      result++;
   }
   
   if (!outputaligned) memcpy(output->data, alignedoutput, sizeof(REAL4)*4*roundedvectorlength);
   
   REAL4 *a = &(input1->data[vec1+4*roundedvectorlength]);
   REAL4 *b = &(input2->data[4*roundedvectorlength]);
   REAL4 *c = &(output->data[4*roundedvectorlength]);
   INT4 n = output->length-4*roundedvectorlength;
   while (n-- > 0) {
      *c = (*a)-(*b);
      a++;
      b++;
      c++;
   }
   
   if (!vec1aligned) XLALFree(allocinput1);
   if (!vec2aligned) XLALFree(allocinput2);
   if (!outputaligned) XLALFree(allocoutput);
#else
   fprintf(stderr, "%s: Failed because SSE is not supported, possibly because -msse flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR_VOID(XLAL_EFAILED);
#endif
   
}


//Compute from a look up table, the sin and cos of a vector of x values using SSE
//Can't use this for x values less than 0 or greater than 2.147483648e9
INT4 sse_sin_cos_2PI_LUT_REAL8Vector(REAL8Vector *sin2pix_vector, REAL8Vector *cos2pix_vector, REAL8Vector *x)
{
   
#ifdef __SSE2__
   INT4 roundedvectorlength = (INT4)x->length / 2;
   INT4 ii;
   
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
   
   REAL8 *allocinput = NULL, *allocoutput1 = NULL, *allocoutput2 = NULL, *alignedinput = NULL, *alignedoutput1 = NULL, *alignedoutput2 = NULL;
   INT4 vecaligned = 0, outputaligned1 = 0, outputaligned2 = 0;
   __m128d *arr1, *sinresult, *cosresult;
   
   //Allocate memory for aligning input vector if necessary
   if ( x->data==(void*)(((UINT8)x->data+15) & ~15) ) {
      vecaligned = 1;
      arr1 = (__m128d*)(void*)x->data;
   } else {
      allocinput = (REAL8*)XLALMalloc(2*roundedvectorlength*sizeof(REAL8) + 15);
      if (allocinput==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 2*roundedvectorlength*sizeof(REAL8) + 15);
         XLAL_ERROR(XLAL_ENOMEM);
      }
      alignedinput = (void*)(((UINT8)allocinput+15) & ~15);
      memcpy(alignedinput, x->data, sizeof(REAL8)*2*roundedvectorlength);
      arr1 = (__m128d*)(void*)alignedinput;
   }
   
   //Allocate memory for aligning output vector 1 if necessary
   if ( sin2pix_vector->data==(void*)(((UINT8)sin2pix_vector->data+15) & ~15) ) {
      outputaligned1 = 1;
      sinresult = (__m128d*)(void*)sin2pix_vector->data;
   } else {
      allocoutput1 = (REAL8*)XLALMalloc(2*roundedvectorlength*sizeof(REAL8) + 15);
      if (allocoutput1==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 2*roundedvectorlength*sizeof(REAL8) + 15);
         XLAL_ERROR(XLAL_ENOMEM);
      }
      alignedoutput1 = (void*)(((UINT8)allocoutput1+15) & ~15);
      sinresult = (__m128d*)(void*)alignedoutput1;
   }
   
   //Allocate memory for aligning output vector 2 if necessary
   if ( cos2pix_vector->data==(void*)(((UINT8)cos2pix_vector->data+15) & ~15) ) {
      outputaligned2 = 1;
      cosresult = (__m128d*)(void*)cos2pix_vector->data;
   } else {
      allocoutput2 = (REAL8*)XLALMalloc(2*roundedvectorlength*sizeof(REAL8) + 15);
      if (allocoutput2==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 2*roundedvectorlength*sizeof(REAL8) + 15);
         XLAL_ERROR(XLAL_ENOMEM);
      }
      alignedoutput2 = (void*)(((UINT8)allocoutput2+15) & ~15);
      cosresult = (__m128d*)(void*)alignedoutput2;
   }
   
   INT4 *intgerVect = (INT4*)XLALMalloc(2*sizeof(INT4)+15);
   if (intgerVect==NULL) {
      fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 2*sizeof(INT4)+15);
      XLAL_ERROR(XLAL_ENOMEM);
   }
   INT4 *I0 = (void*)(((UINT8)intgerVect+15) & ~15);
   memset(I0, 0, sizeof(INT4)*2);
   
   __m128d lutresf = _mm_set1_pd(LUT_RES_F);
   __m128d onehalf = _mm_set1_pd(0.5);
   __m128d oolutres = _mm_set1_pd(OO_LUT_RES);
   __m128d twopi = _mm_set1_pd(LAL_TWOPI);
   
   for (ii=0; ii<roundedvectorlength; ii++) {
      
      //Fractional part of x [0,1)
      __m128i xfloor = _mm_cvttpd_epi32(*arr1);
      __m128d xfloord = _mm_cvtepi32_pd(xfloor);
      __m128d xt = _mm_sub_pd(*arr1, xfloord);
      
      //i0 in [0, LUT_RES]
      __m128d d0 = _mm_mul_pd(xt, lutresf);
      __m128d d1 = _mm_add_pd(d0, onehalf);
      __m128i i0 = _mm_cvttpd_epi32(d1);
      __m128d i0d = _mm_cvtepi32_pd(i0);
      //__m128i *i_0 = &i0;
      __m128i i0copy = i0;
      
      //d and d2
      __m128d d_0 = _mm_mul_pd(oolutres, i0d);
      __m128d d_1 = _mm_sub_pd(xt, d_0);
      __m128d d = _mm_mul_pd(d_1, twopi);
      __m128d d2 = _mm_mul_pd(d, d);
      d2 = _mm_mul_pd(d2, onehalf);
      
      //I0 = (INT4*)i_0;
      I0[0] = _mm_cvtsi128_si32(_mm_srli_si128(i0copy, 0));
      I0[1] = _mm_cvtsi128_si32(_mm_srli_si128(i0copy, 4));
      __m128d ts = _mm_setr_pd(sinVal[I0[0]], sinVal[I0[1]]);
      __m128d tc = _mm_setr_pd(cosVal[I0[0]], cosVal[I0[1]]);
      
      __m128d dtimestc = _mm_mul_pd(d, tc);
      __m128d d2timests = _mm_mul_pd(d2, ts);
      __m128d dtimests = _mm_mul_pd(d, ts);
      __m128d d2timestc = _mm_mul_pd(d2, tc);
      __m128d tsplusdtimestc = _mm_add_pd(ts, dtimestc);
      __m128d tcminusdtimests = _mm_sub_pd(tc, dtimests);
      *sinresult = _mm_sub_pd(tsplusdtimestc, d2timests);
      *cosresult = _mm_sub_pd(tcminusdtimests, d2timestc);
      
      arr1++;
      sinresult++;
      cosresult++;
   }
   
   //Copy output aligned memory to non-aligned memory if necessary
   if (!outputaligned1) memcpy(sin2pix_vector->data, alignedoutput1, 2*roundedvectorlength*sizeof(REAL8));
   if (!outputaligned2) memcpy(cos2pix_vector->data, alignedoutput2, 2*roundedvectorlength*sizeof(REAL8));
   
   //Finish up the remaining part
   for (ii=2*roundedvectorlength; ii<(INT4)x->length; ii++) {
      twospect_sin_cos_2PI_LUT(&(sin2pix_vector->data[ii]), &(cos2pix_vector->data[ii]), x->data[ii]);
   }
   
   //Free memory if necessary
   if (!vecaligned) XLALFree(allocinput);
   if (!outputaligned1) XLALFree(allocoutput1);
   if (!outputaligned2) XLALFree(allocoutput2);
   
   //Free memory
   XLALFree(intgerVect);
   
   return XLAL_SUCCESS;
#else
   fprintf(stderr, "%s: Failed because SSE2 is not supported, possibly because -msse2 flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR(XLAL_EFAILED);
#endif
   
} /* sse_sin_cos_2PI_LUT_REAL8Vector() */


//Compute from a look up table, the sin and cos of a vector of x values using SSE
//Can't use this for x values less than 0 or greater than 2.147483647e9
INT4 sse_sin_cos_2PI_LUT_REAL4Vector(REAL4Vector *sin2pix_vector, REAL4Vector *cos2pix_vector, REAL4Vector *x)
{
   
#ifdef __SSE2__
   INT4 roundedvectorlength = (INT4)x->length / 4;
   INT4 ii;
   
   static BOOLEAN firstCall = TRUE;
   static REAL4 sinVal[LUT_RES+1], cosVal[LUT_RES+1];
   
   /* the first time we get called, we set up the lookup-table */
   if ( firstCall ) {
      UINT4 k;
      for (k=0; k <= LUT_RES; k++) {
         sinVal[k] = (REAL4)sinf( (REAL4)(LAL_TWOPI * k * OO_LUT_RES) );
         cosVal[k] = (REAL4)cosf( (REAL4)(LAL_TWOPI * k * OO_LUT_RES) );
      }
      firstCall = FALSE;
   }
   
   REAL4 *allocinput = NULL, *allocoutput1 = NULL, *allocoutput2 = NULL, *alignedinput = NULL, *alignedoutput1 = NULL, *alignedoutput2 = NULL;
   INT4 vecaligned = 0, outputaligned1 = 0, outputaligned2 = 0;
   __m128 *arr1, *sinresult, *cosresult;
   
   //Allocate memory for aligning input vector if necessary
   if ( x->data==(void*)(((UINT8)x->data+15) & ~15) ) {
      vecaligned = 1;
      arr1 = (__m128*)(void*)x->data;
   } else {
      allocinput = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15);
      if (allocinput==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 4*roundedvectorlength*sizeof(REAL4) + 15);
         XLAL_ERROR(XLAL_ENOMEM);
      }
      alignedinput = (void*)(((UINT8)allocinput+15) & ~15);
      memcpy(alignedinput, x->data, sizeof(REAL4)*4*roundedvectorlength);
      arr1 = (__m128*)(void*)alignedinput;
   }
   
   //Allocate memory for aligning output vector 1 if necessary
   if ( sin2pix_vector->data==(void*)(((UINT8)sin2pix_vector->data+15) & ~15) ) {
      outputaligned1 = 1;
      sinresult = (__m128*)(void*)sin2pix_vector->data;
   } else {
      allocoutput1 = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15);
      if (allocoutput1==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 4*roundedvectorlength*sizeof(REAL4) + 15);
         XLAL_ERROR(XLAL_ENOMEM);
      }
      alignedoutput1 = (void*)(((UINT8)allocoutput1+15) & ~15);
      sinresult = (__m128*)(void*)alignedoutput1;
   }
   
   //Allocate memory for aligning output vector 2 if necessary
   if ( cos2pix_vector->data==(void*)(((UINT8)cos2pix_vector->data+15) & ~15) ) {
      outputaligned2 = 1;
      cosresult = (__m128*)(void*)cos2pix_vector->data;
   } else {
      allocoutput2 = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15);
      if (allocoutput2==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 4*roundedvectorlength*sizeof(REAL4) + 15);
         XLAL_ERROR(XLAL_ENOMEM);
      }
      alignedoutput2 = (void*)(((UINT8)allocoutput2+15) & ~15);
      cosresult = (__m128*)(void*)alignedoutput2;
   }
   
   INT4 *intgerVect = (INT4*)XLALMalloc(4*sizeof(INT4)+15);
   if (intgerVect==NULL) {
      fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 4*sizeof(INT4)+15);
      XLAL_ERROR(XLAL_ENOMEM);
   }
   INT4 *I0 = (void*)(((UINT8)intgerVect+15) & ~15);
   memset(I0, 0, sizeof(INT4)*4);
   
   __m128 lutresf = _mm_set1_ps((REAL4)LUT_RES_F);
   __m128 onehalf = _mm_set1_ps(0.5f);
   __m128 oolutres = _mm_set1_ps((REAL4)OO_LUT_RES);
   __m128 twopi = _mm_set1_ps((REAL4)LAL_TWOPI);
   
   for (ii=0; ii<roundedvectorlength; ii++) {
      
      //Fractional part of x [0,1)
      __m128i xfloor = _mm_cvttps_epi32(*arr1);
      __m128 xfloors = _mm_cvtepi32_ps(xfloor);
      __m128 xt = _mm_sub_ps(*arr1, xfloors);
      
      //i0 in [0, LUT_RES]
      __m128 d0 = _mm_mul_ps(xt, lutresf);
      __m128 d1 = _mm_add_ps(d0, onehalf);
      __m128i i0 = _mm_cvttps_epi32(d1);
      __m128 i0s = _mm_cvtepi32_ps(i0);
      __m128i i0copy = i0;
      
      //d and d2
      __m128 d_0 = _mm_mul_ps(oolutres, i0s);
      __m128 d_1 = _mm_sub_ps(xt, d_0);
      __m128 d = _mm_mul_ps(d_1, twopi);
      __m128 d2 = _mm_mul_ps(d, d);
      d2 = _mm_mul_ps(d2, onehalf);
      
      //I0 = (INT4*)i_0;
      I0[0] = _mm_cvtsi128_si32(_mm_srli_si128(i0copy, 0));
      I0[1] = _mm_cvtsi128_si32(_mm_srli_si128(i0copy, 4));
      I0[2] = _mm_cvtsi128_si32(_mm_srli_si128(i0copy, 8));
      I0[3] = _mm_cvtsi128_si32(_mm_srli_si128(i0copy, 12));
      __m128 ts = _mm_setr_ps(sinVal[I0[0]], sinVal[I0[1]], sinVal[I0[2]], sinVal[I0[3]]);
      __m128 tc = _mm_setr_ps(cosVal[I0[0]], cosVal[I0[1]], cosVal[I0[2]], cosVal[I0[3]]);
      
      __m128 dtimestc = _mm_mul_ps(d, tc);
      __m128 d2timests = _mm_mul_ps(d2, ts);
      __m128 dtimests = _mm_mul_ps(d, ts);
      __m128 d2timestc = _mm_mul_ps(d2, tc);
      __m128 tsplusdtimestc = _mm_add_ps(ts, dtimestc);
      __m128 tcminusdtimests = _mm_sub_ps(tc, dtimests);
      *sinresult = _mm_sub_ps(tsplusdtimestc, d2timests);
      *cosresult = _mm_sub_ps(tcminusdtimests, d2timestc);
      
      arr1++;
      sinresult++;
      cosresult++;
   }
   
   //Copy output aligned memory to non-aligned memory if necessary
   if (!outputaligned1) memcpy(sin2pix_vector->data, alignedoutput1, 4*roundedvectorlength*sizeof(REAL4));
   if (!outputaligned2) memcpy(cos2pix_vector->data, alignedoutput2, 4*roundedvectorlength*sizeof(REAL4));
   
   //Finish up the remaining part
   REAL8 sinval = 0.0, cosval = 0.0;
   for (ii=4*roundedvectorlength; ii<(INT4)x->length; ii++) {
      twospect_sin_cos_2PI_LUT(&sinval, &cosval, x->data[ii]);
      sin2pix_vector->data[ii] = (REAL4)sinval;
      cos2pix_vector->data[ii] = (REAL4)cosval;
   }
   
   //Free memory if necessary
   if (!vecaligned) XLALFree(allocinput);
   if (!outputaligned1) XLALFree(allocoutput1);
   if (!outputaligned2) XLALFree(allocoutput2);
   
   //Free memory
   XLALFree(intgerVect);
   
   return XLAL_SUCCESS;
#else
   fprintf(stderr, "%s: Failed because SSE2 is not supported, possibly because -msse2 flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR(XLAL_EFAILED);
#endif
   
} /* sse_sin_cos_2PI_LUT_REAL4Vector() */



//Exponential of input vector is computed using SSE
//Cephes library based
void sse_exp_REAL8Vector(REAL8Vector *output, REAL8Vector *input)
{
   
#ifdef __SSE2__
   INT4 roundedvectorlength = (INT4)input->length / 2;
   INT4 ii;
   
   REAL8 *allocinput = NULL, *allocoutput = NULL, *alignedinput = NULL, *alignedoutput = NULL;
   INT4 vecaligned = 0, outputaligned = 0;
   
   __m128d *x, *result;
   
   //Allocate memory for aligning input vector if necessary
   if ( input->data==(void*)(((UINT8)input->data+15) & ~15) ) {
      vecaligned = 1;
      x = (__m128d*)(void*)input->data;
   } else {
      allocinput = (REAL8*)XLALMalloc(2*roundedvectorlength*sizeof(REAL8) + 15);
      if (allocinput==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 2*roundedvectorlength*sizeof(REAL8) + 15);
         XLAL_ERROR_VOID(XLAL_ENOMEM);
      }
      alignedinput = (void*)(((UINT8)allocinput+15) & ~15);
      memcpy(alignedinput, input->data, sizeof(REAL8)*2*roundedvectorlength);
      x = (__m128d*)(void*)alignedinput;
   }
   
   //Allocate memory for aligning output vector 1 if necessary
   if ( output->data==(void*)(((UINT8)output->data+15) & ~15) ) {
      outputaligned = 1;
      result = (__m128d*)(void*)output->data;
   } else {
      allocoutput = (REAL8*)XLALMalloc(2*roundedvectorlength*sizeof(REAL8) + 15);
      if (allocoutput==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 2*roundedvectorlength*sizeof(REAL8) + 15);
         XLAL_ERROR_VOID(XLAL_ENOMEM);
      }
      alignedoutput = (void*)(((UINT8)allocoutput+15) & ~15);
      result = (__m128d*)(void*)alignedoutput;
   }
   
   __m128i expoffset = _mm_set_epi32(0, 1023, 0, 1023); //__m128i expoffset = _mm_set1_epi64x(1023);      //Exponent mask for double precision
   __m128i maskupper32bits = _mm_set_epi32(0xffffffff, 0x00000000, 0xffffffff, 0x00000000); //__m128i maskupper32bits = _mm_set1_epi64x(0xffffffff00000000);    //mask for upper 32 bits
   __m128i masklower32bits = _mm_set_epi32(0x00000000, 0xffffffff, 0x00000000, 0xffffffff); //__m128i masklower32bits = _mm_set1_epi64x(0x00000000ffffffff);    //mask for lower 32 bits
   __m128d log2e = _mm_set1_pd(1.442695040888963);                   //ln(2)
   __m128d onehalf = _mm_set1_pd(0.5);             //0.5
   __m128d one = _mm_set1_pd(1.0);                 //1.0
   __m128d two = _mm_set1_pd(2.0);                 //2.0
   __m128d maxexp = _mm_set1_pd(7.007827128933840e+02);
   __m128d minexp = _mm_set1_pd(-7.007827128933840e+02);
   
   __m128d cephes_c1 = _mm_set1_pd(6.93145751953125e-1);
   __m128d cephes_c2 = _mm_set1_pd(1.42860682030941723212e-6);
   __m128d cephes_p0 = _mm_set1_pd(1.26177193074810590878e-4);
   __m128d cephes_p1 = _mm_set1_pd(3.02994407707441961300e-2);
   __m128d cephes_p2 = _mm_set1_pd(9.99999999999999999910e-1);
   __m128d cephes_q0 = _mm_set1_pd(3.00198505138664455042e-6);
   __m128d cephes_q1 = _mm_set1_pd(2.52448340349684104192e-3);
   __m128d cephes_q2 = _mm_set1_pd(2.27265548208155028766e-1);
   __m128d cephes_q3 = _mm_set1_pd(2.00000000000000000009e0);
   
   for (ii=0; ii<roundedvectorlength; ii++) {
      __m128d y = *x;
      y = _mm_max_pd(y, minexp);
      y = _mm_min_pd(y, maxexp);
      
      //Compute x*ln(2)
      __m128d log2etimesx = _mm_mul_pd(log2e, y);
      
      //Round
      //__m128d log2etimesxsubhalf = _mm_sub_pd(log2etimesx, onehalf);
      //__m128i log2etimesx_rounded_i = _mm_cvttpd_epi32(log2etimesxsubhalf);
      //__m128d log2etimesx_rounded = _mm_cvtepi32_pd(log2etimesx_rounded_i);
      __m128d log2etimesxplushalf = _mm_add_pd(log2etimesx, onehalf);
      __m128i log2etimesx_rounded_i = _mm_cvttpd_epi32(log2etimesxplushalf);
      __m128d log2etimesx_rounded = _mm_cvtepi32_pd(log2etimesx_rounded_i);
      __m128d mask = _mm_cmpgt_pd(log2etimesx_rounded, log2etimesxplushalf);
      mask = _mm_and_pd(mask, one);
      log2etimesx_rounded = _mm_sub_pd(log2etimesx_rounded, mask);
      log2etimesx_rounded_i = _mm_cvttpd_epi32(log2etimesx_rounded);
      
      //multiply and subtract as in the cephes code
      __m128d log2etimesx_rounded_times_c1 = _mm_mul_pd(log2etimesx_rounded, cephes_c1);
      __m128d log2etimesx_rounded_times_c2 = _mm_mul_pd(log2etimesx_rounded, cephes_c2);
      y = _mm_sub_pd(y, log2etimesx_rounded_times_c1);
      y = _mm_sub_pd(y, log2etimesx_rounded_times_c2);
      
      //x**2
      __m128d xsq = _mm_mul_pd(y, y);
      
      //Pade approximation with polynomials
      //Now the polynomial part 1
      __m128d polevlresult = cephes_p0;
      polevlresult = _mm_mul_pd(polevlresult, xsq);
      polevlresult = _mm_add_pd(polevlresult, cephes_p1);
      polevlresult = _mm_mul_pd(polevlresult, xsq);
      polevlresult = _mm_add_pd(polevlresult, cephes_p2);
      __m128d xtimespolresult = _mm_mul_pd(y, polevlresult);
      
      //And polynomial part 2
      polevlresult = cephes_q0;
      polevlresult = _mm_mul_pd(polevlresult, xsq);
      polevlresult = _mm_add_pd(polevlresult, cephes_q1);
      polevlresult = _mm_mul_pd(polevlresult, xsq);
      polevlresult = _mm_add_pd(polevlresult, cephes_q2);
      polevlresult = _mm_mul_pd(polevlresult, xsq);
      polevlresult = _mm_add_pd(polevlresult, cephes_q3);
      __m128d polevlresult_minus_xtimespolresult = _mm_sub_pd(polevlresult, xtimespolresult);
      y = _mm_div_pd(xtimespolresult, polevlresult_minus_xtimespolresult);
      
      //Finishing the pade approximation
      y = _mm_mul_pd(y, two);
      y = _mm_add_pd(y, one);
      
      //Need to get the integer to a 64 bit format for below even though we don't need that much
      //Note that this is still limited to an INT4 size
      //It's hard because there is no _mm_srai_epi64 command, so we do it ourselves
      __m128i log2etimesx_rounded_i_64bit = _mm_unpacklo_epi32(log2etimesx_rounded_i, log2etimesx_rounded_i);  //Interleve with itself
      __m128i maskedlower32 = _mm_and_si128(log2etimesx_rounded_i_64bit, masklower32bits);      //Save lower 32 bits
      log2etimesx_rounded_i_64bit = _mm_srai_epi32(log2etimesx_rounded_i_64bit, 32);            //Shift with the sign bit
      log2etimesx_rounded_i_64bit = _mm_and_si128(log2etimesx_rounded_i_64bit, maskupper32bits);   //Discard now useless lower 32 bits
      log2etimesx_rounded_i_64bit = _mm_xor_si128(log2etimesx_rounded_i_64bit, maskedlower32);  //Restore original lower 32 bits
      
      //Now construct 2**n
      __m128i log2etimesx_rounded_i_with_offset = _mm_add_epi64(log2etimesx_rounded_i_64bit, expoffset);
      log2etimesx_rounded_i_with_offset = _mm_slli_epi64(log2etimesx_rounded_i_with_offset, 52);
      __m128d pow2n = _mm_castsi128_pd(log2etimesx_rounded_i_with_offset);
      
      //And multiply
      *result = _mm_mul_pd(y, pow2n);
      
      x++;
      result++;
      
   }
   
   //Copy output aligned memory to non-aligned memory if necessary
   if (!outputaligned) memcpy(output->data, alignedoutput, 2*roundedvectorlength*sizeof(REAL8));
   
   //Finish up the remaining part
   for (ii=2*roundedvectorlength; ii<(INT4)input->length; ii++) {
      output->data[ii] = exp(input->data[ii]);
   }
   
   //Free memory if necessary
   if (!vecaligned) XLALFree(allocinput);
   if (!outputaligned) XLALFree(allocoutput);
#else
   fprintf(stderr, "%s: Failed because SSE2 is not supported, possibly because -msse2 flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR_VOID(XLAL_EFAILED);
#endif
   
}


//Exponential of input vector is computed using SSE
//Cephes library based
void sse_exp_REAL4Vector(REAL4Vector *output, REAL4Vector *input)
{
   
#ifdef __SSE2__
   INT4 roundedvectorlength = (INT4)input->length / 4;
   INT4 ii;
   
   REAL4 *allocinput = NULL, *allocoutput = NULL, *alignedinput = NULL, *alignedoutput = NULL;
   INT4 vecaligned = 0, outputaligned = 0;
   
   __m128 *x, *result;
   
   //Allocate memory for aligning input vector if necessary
   if ( input->data==(void*)(((UINT8)input->data+15) & ~15) ) {
      vecaligned = 1;
      x = (__m128*)(void*)input->data;
   } else {
      allocinput = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15);
      if (allocinput==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 4*roundedvectorlength*sizeof(REAL4) + 15);
         XLAL_ERROR_VOID(XLAL_ENOMEM);
      }
      alignedinput = (void*)(((UINT8)allocinput+15) & ~15);
      memcpy(alignedinput, input->data, sizeof(REAL4)*4*roundedvectorlength);
      x = (__m128*)(void*)alignedinput;
   }
   
   //Allocate memory for aligning output vector 1 if necessary
   if ( output->data==(void*)(((UINT8)output->data+15) & ~15) ) {
      outputaligned = 1;
      result = (__m128*)(void*)output->data;
   } else {
      allocoutput = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15);
      if (allocoutput==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 4*roundedvectorlength*sizeof(REAL4) + 15);
         XLAL_ERROR_VOID(XLAL_ENOMEM);
      }
      alignedoutput = (void*)(((UINT8)allocoutput+15) & ~15);
      result = (__m128*)(void*)alignedoutput;
   }
   
   __m128i expoffset = _mm_set1_epi32(127);      //Exponent mask for single precision
   __m128 log2e = _mm_set1_ps(1.442695040888963f);                   //ln(2)
   __m128 onehalf = _mm_set1_ps(0.5f);             //0.5
   __m128 one = _mm_set1_ps(1.0f);                 //1.0
   __m128 two = _mm_set1_ps(2.0f);                 //2.0
   __m128 maxexp = _mm_set1_ps(88.0f);
   __m128 minexp = _mm_set1_ps(-88.0f);
   
   __m128 cephes_c1 = _mm_set1_ps(6.93145751953125e-1f);
   __m128 cephes_c2 = _mm_set1_ps(1.42860682030941723212e-6f);
   __m128 cephes_p0 = _mm_set1_ps(1.26177193074810590878e-4f);
   __m128 cephes_p1 = _mm_set1_ps(3.02994407707441961300e-2f);
   __m128 cephes_p2 = _mm_set1_ps(9.99999999999999999910e-1f);
   __m128 cephes_q0 = _mm_set1_ps(3.00198505138664455042e-6f);
   __m128 cephes_q1 = _mm_set1_ps(2.52448340349684104192e-3f);
   __m128 cephes_q2 = _mm_set1_ps(2.27265548208155028766e-1f);
   __m128 cephes_q3 = _mm_set1_ps(2.00000000000000000009e0f);
   
   for (ii=0; ii<roundedvectorlength; ii++) {
      __m128 y = *x;
      y = _mm_max_ps(y, minexp);
      y = _mm_min_ps(y, maxexp);
      
      //Compute x*ln(2)
      __m128 log2etimesx = _mm_mul_ps(log2e, y);
      
      //Round
      //__m128 log2etimesxsubhalf = _mm_sub_ps(log2etimesx, onehalf);
      //__m128i log2etimesx_rounded_i = _mm_cvttps_epi32(log2etimesxsubhalf);
      //__m128 log2etimesx_rounded = _mm_cvtepi32_ps(log2etimesx_rounded_i);
      __m128 log2etimesxplushalf = _mm_add_ps(log2etimesx, onehalf);
      __m128i log2etimesx_rounded_i = _mm_cvttps_epi32(log2etimesxplushalf);
      __m128 log2etimesx_rounded = _mm_cvtepi32_ps(log2etimesx_rounded_i);
      __m128 mask = _mm_cmpgt_ps(log2etimesx_rounded, log2etimesxplushalf);
      mask = _mm_and_ps(mask, one);
      log2etimesx_rounded = _mm_sub_ps(log2etimesx_rounded, mask);
      log2etimesx_rounded_i = _mm_cvttps_epi32(log2etimesx_rounded);
      
      //multiply and subtract as in the cephes code
      __m128 log2etimesx_rounded_times_c1 = _mm_mul_ps(log2etimesx_rounded, cephes_c1);
      __m128 log2etimesx_rounded_times_c2 = _mm_mul_ps(log2etimesx_rounded, cephes_c2);
      y = _mm_sub_ps(y, log2etimesx_rounded_times_c1);
      y = _mm_sub_ps(y, log2etimesx_rounded_times_c2);
      
      //x**2
      __m128 xsq = _mm_mul_ps(y, y);
      
      //Pade approximation with polynomials
      //Now the polynomial part 1
      __m128 polevlresult = cephes_p0;
      polevlresult = _mm_mul_ps(polevlresult, xsq);
      polevlresult = _mm_add_ps(polevlresult, cephes_p1);
      polevlresult = _mm_mul_ps(polevlresult, xsq);
      polevlresult = _mm_add_ps(polevlresult, cephes_p2);
      __m128 xtimespolresult = _mm_mul_ps(y, polevlresult);
      
      //And polynomial part 2
      polevlresult = cephes_q0;
      polevlresult = _mm_mul_ps(polevlresult, xsq);
      polevlresult = _mm_add_ps(polevlresult, cephes_q1);
      polevlresult = _mm_mul_ps(polevlresult, xsq);
      polevlresult = _mm_add_ps(polevlresult, cephes_q2);
      polevlresult = _mm_mul_ps(polevlresult, xsq);
      polevlresult = _mm_add_ps(polevlresult, cephes_q3);
      __m128 polevlresult_minus_xtimespolresult = _mm_sub_ps(polevlresult, xtimespolresult);
      y = _mm_div_ps(xtimespolresult, polevlresult_minus_xtimespolresult);
      
      //Finishing the pade approximation
      y = _mm_mul_ps(y, two);
      y = _mm_add_ps(y, one);
      
      //Now construct 2**n
      __m128i log2etimesx_rounded_i_with_offset = _mm_add_epi32(log2etimesx_rounded_i, expoffset);
      log2etimesx_rounded_i_with_offset = _mm_slli_epi32(log2etimesx_rounded_i_with_offset, 23);
      __m128 pow2n = _mm_castsi128_ps(log2etimesx_rounded_i_with_offset);
      
      //And multiply
      *result = _mm_mul_ps(y, pow2n);
      
      x++;
      result++;
      
   }
   
   //Copy output aligned memory to non-aligned memory if necessary
   if (!outputaligned) memcpy(output->data, alignedoutput, 4*roundedvectorlength*sizeof(REAL4));
   
   //Finish up the remaining part
   for (ii=4*roundedvectorlength; ii<(INT4)input->length; ii++) {
      output->data[ii] = expf(input->data[ii]);
   }
   
   /* FILE *EXPVALS = fopen("./output/expvals.dat","w");
   for (ii=0; ii<(INT4)output->length; ii++) {
      fprintf(EXPVALS, "%f\n", output->data[ii]);
   }
   fclose(EXPVALS); */
   
   //Free memory if necessary
   if (!vecaligned) XLALFree(allocinput);
   if (!outputaligned) XLALFree(allocoutput);
#else
   fprintf(stderr, "%s: Failed because SSE2 is not supported, possibly because -msse2 flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR_VOID(XLAL_EFAILED);
#endif
   
}




//Arctangent of input vector is computed using SSE
//Cephes library based
// !!!!!!!! NOT READY TO BE USED!!!!!!!!!
/* void sse_atan_REAL8Vector(REAL8Vector *output, REAL8Vector *input)
{
   
#ifdef __SSE2__
   INT4 roundedvectorlength = (INT4)input->length / 2;
   INT4 ii;
   
   REAL8 *allocinput = NULL, *allocoutput = NULL, *alignedinput = NULL, *alignedoutput = NULL;
   INT4 vecaligned = 0, outputaligned = 0;
   
   __m128d *x, *result;
   
   //Allocate memory for aligning input vector if necessary
   if ( input->data==(void*)(((UINT8)input->data+15) & ~15) ) {
      vecaligned = 1;
      x = (__m128d*)(void*)input->data;
   } else {
      allocinput = (REAL8*)XLALMalloc(2*roundedvectorlength*sizeof(REAL8) + 15);
      if (allocinput==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 2*roundedvectorlength*sizeof(REAL8) + 15);
         XLAL_ERROR_VOID(XLAL_ENOMEM);
      }
      alignedinput = (void*)(((UINT8)allocinput+15) & ~15);
      memcpy(alignedinput, input->data, sizeof(REAL8)*2*roundedvectorlength);
      x = (__m128d*)(void*)alignedinput;
   }
   
   //Allocate memory for aligning output vector 1 if necessary
   if ( output->data==(void*)(((UINT8)output->data+15) & ~15) ) {
      outputaligned = 1;
      result = (__m128d*)(void*)output->data;
   } else {
      allocoutput = (REAL8*)XLALMalloc(2*roundedvectorlength*sizeof(REAL8) + 15);
      if (allocoutput==NULL) {
         fprintf(stderr, "%s: XLALMalloc(%zu) failed.\n", __func__, 2*roundedvectorlength*sizeof(REAL8) + 15);
         XLAL_ERROR_VOID(XLAL_ENOMEM);
      }
      alignedoutput = (void*)(((UINT8)allocoutput+15) & ~15);
      result = (__m128d*)(void*)alignedoutput;
   }
   
   //Define values
   __m128d cephes_p0 = _mm_set1_pd(-8.750608600031904122785e-1);
   __m128d cephes_p1 = _mm_set1_pd(-1.615753718733365076637e1);
   __m128d cephes_p2 = _mm_set1_pd(-7.500855792314704667340e1);
   __m128d cephes_p3 = _mm_set1_pd(-1.228866684490136173410e2);
   __m128d cephes_p4 = _mm_set1_pd(-6.485021904942025371773e1);
   __m128d cephes_q0 = _mm_set1_pd(2.485846490142306297962e1);
   __m128d cephes_q1 = _mm_set1_pd(1.650270098316988542046e2);
   __m128d cephes_q2 = _mm_set1_pd(4.328810604912902668951e2);
   __m128d cephes_q3 = _mm_set1_pd(4.853903996359136964868e2);
   __m128d cephes_q4 = _mm_set1_pd(1.945506571482613964425e2);
   __m128d tan3pO8 = _mm_set1_pd(2.41421356237309504880);
   __m128d almostTwoThirds = _mm_set1_pd(0.66);
   __m128d piOverTwo = _mm_set1_pd(LAL_PI_2);
   __m128d piOverFour = _mm_set1_pd(LAL_PI_4);
   __m128d zero = _mm_set1_pd(0.0);
   __m128d minusOne = _mm_set1_pd(-1.0);
   __m128d small = _mm_set1_pd(1.0e-50);
   __m128d morebits = _mm_set1_pd(6.123233995736765886130e-17);
   __m128d halfmorebits = _mm_set1_pd(3.061616997868382943065e-17);
   
   for (ii=0; ii<roundedvectorlength; ii++) {
      
      //Save sign values
      //__m128d signs = _mm_and_pd(*x, signbit);
      
      //Make values positive
      //__m128d intx = _mm_xor_pd(*x, signs);
      
      //Assume positive
      __m128d intx = *x;
      
      //Reduce range
      __m128d greaterThanAlmostTwoThirds = _mm_cmpgt_pd(intx, almostTwoThirds);
      __m128d greaterThanTan3pO8 = _mm_cmpgt_pd(intx, tan3pO8);
      __m128d inBetween = _mm_and_pd( _mm_cmple_pd(intx, tan3pO8), greaterThanAlmostTwoThirds);
      __m128d lessThanAlmostTwoThirds = _mm_cmple_pd(intx, almostTwoThirds);
      __m128d y = zero;
      __m128d x1 = _mm_and_pd( _mm_div_pd(minusOne, _mm_add_pd(intx, small)), greaterThanTan3pO8);
      __m128d y_1 = _mm_and_pd(piOverTwo, greaterThanTan3pO8);
      __m128d x2 = _mm_and_pd( _mm_div_pd( _mm_add_pd(intx, minusOne), _mm_sub_pd(intx, minusOne)), inBetween);
      __m128d y_2 = _mm_and_pd(piOverFour, inBetween);
      y = _mm_add_pd(y, y_1);
      y = _mm_add_pd(y, y_2);
      intx = _mm_and_pd(intx, lessThanAlmostTwoThirds);
      intx = _mm_add_pd(intx, x1);
      intx = _mm_add_pd(intx, x2);
      
      __m128d z = _mm_mul_pd(intx, intx);
      __m128d polevlresult = cephes_p0;
      polevlresult = _mm_mul_pd(polevlresult, z);
      polevlresult = _mm_add_pd(polevlresult, cephes_p1);
      polevlresult = _mm_mul_pd(polevlresult, z);
      polevlresult = _mm_add_pd(polevlresult, cephes_p2);
      polevlresult = _mm_mul_pd(polevlresult, z);
      polevlresult = _mm_add_pd(polevlresult, cephes_p3);
      polevlresult = _mm_mul_pd(polevlresult, z);
      polevlresult = _mm_add_pd(polevlresult, cephes_p4);
      __m128d zTimesPolEvlResult = _mm_mul_pd(z, polevlresult);
      
      polevlresult = z;
      polevlresult = _mm_add_pd(polevlresult, cephes_q0);
      polevlresult = _mm_mul_pd(polevlresult, z);
      polevlresult = _mm_add_pd(polevlresult, cephes_q1);
      polevlresult = _mm_mul_pd(polevlresult, z);
      polevlresult = _mm_add_pd(polevlresult, cephes_q2);
      polevlresult = _mm_mul_pd(polevlresult, z);
      polevlresult = _mm_add_pd(polevlresult, cephes_q3);
      polevlresult = _mm_mul_pd(polevlresult, z);
      polevlresult = _mm_add_pd(polevlresult, cephes_q4);
      z = _mm_div_pd(zTimesPolEvlResult, polevlresult);
      
      z = _mm_add_pd( _mm_mul_pd(intx, z), intx);
      
      __m128d addextrabits = zero;
      addextrabits = _mm_add_pd(addextrabits, _mm_and_pd(morebits, greaterThanTan3pO8) );
      addextrabits = _mm_add_pd(addextrabits, _mm_and_pd(halfmorebits, inBetween) );
      z = _mm_add_pd(z, addextrabits);
      
      y = _mm_add_pd(y, z);
      
      // *result = _mm_xor_pd(y, signs);
      *result = y;
      
      x++;
      result++;
      
   }
   
   //Copy output aligned memory to non-aligned memory if necessary
   if (!outputaligned) memcpy(output->data, alignedoutput, 2*roundedvectorlength*sizeof(REAL8));
   
   //Finish up the remaining part
   for (ii=2*roundedvectorlength; ii<(INT4)input->length; ii++) {
      output->data[ii] = atan(input->data[ii]);
   }
   
   //Free memory if necessary
   if (!vecaligned) XLALFree(allocinput);
   if (!outputaligned) XLALFree(allocoutput);
#else
   fprintf(stderr, "%s: Failed because SSE2 is not supported, possibly because -msse2 flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR_VOID(XLAL_EFAILED);
#endif
   
} */




