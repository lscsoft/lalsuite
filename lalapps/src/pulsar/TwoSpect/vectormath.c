/*
 *  Copyright (C) 2011, 2012, 2014, 2015 Evan Goetz
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
#ifdef __AVX__
#include <immintrin.h>
#endif

#include <math.h>

#include <lal/LALConstants.h>
#include <lal/SinCosLUT.h>

#include "vectormath.h"

#ifdef HAVE_STDINT_H
#include <stdint.h>
#else
typedef size_t uintptr_t;
#endif

#define OOTWOPI         (1.0 / LAL_TWOPI)
#define LUT_RES         1024                 /* resolution of lookup-table */
#define LUT_RES_F       (1.0 * LUT_RES)
#define OO_LUT_RES      (1.0 / LUT_RES)
#define X_TO_IND        (1.0 * LUT_RES * OOTWOPI )
#define IND_TO_X        (LAL_TWOPI * OO_LUT_RES)
#define TRUE            (1==1)
#define FALSE           (1==0)


alignedREAL8Vector * createAlignedREAL8Vector(UINT4 length, const size_t align)
{
   alignedREAL8Vector *vector;
   XLAL_CHECK_NULL( (vector = XLALMalloc(sizeof(*vector))) != NULL, XLAL_ENOMEM );
   vector->length = length;
   UINT4 paddedLength = length + align - 1;
   XLAL_CHECK_NULL( (vector->data0 = XLALMalloc(paddedLength*sizeof(REAL8))) != NULL, XLAL_ENOMEM );
   size_t remBytes = ((size_t)vector->data0) % align;
   size_t offsetBytes = (align - remBytes) % align;
   vector->data = (void*)(((char*)vector->data0) + offsetBytes);
   return vector;
}
void destroyAlignedREAL8Vector(alignedREAL8Vector *vector)
{
   if (!vector) return;
   if (vector->data0) XLALFree(vector->data0);
   XLALFree(vector);
}
alignedREAL8VectorArray * createAlignedREAL8VectorArray(const UINT4 length, const UINT4 vectorLength, const size_t align)
{
   alignedREAL8VectorArray *array = NULL;
   XLAL_CHECK_NULL( (array = XLALMalloc(sizeof(*array))) != NULL, XLAL_ENOMEM );
   array->length = length;
   XLAL_CHECK_NULL( (array->data = XLALMalloc(sizeof(*(array->data))*array->length)) != NULL, XLAL_ENOMEM );
   for (UINT4 ii=0; ii<length; ii++) {
      XLAL_CHECK_NULL( (array->data[ii] = createAlignedREAL8Vector(vectorLength, align)) != NULL, XLAL_EFUNC );
   }
   return array;
}
void destroyAlignedREAL8VectorArray(alignedREAL8VectorArray *array)
{
   if (!array) return;
   for (UINT4 ii=0; ii<array->length; ii++) {
      destroyAlignedREAL8Vector(array->data[ii]);
   }
   XLALFree(array->data);
   XLALFree(array);
}
REAL4VectorAlignedArray * createREAL4VectorAlignedArray(const UINT4 length, const UINT4 vectorLength, const size_t align)
{
   REAL4VectorAlignedArray *array = NULL;
   XLAL_CHECK_NULL( (array = XLALMalloc(sizeof(*array))) != NULL, XLAL_ENOMEM );
   array->length = length;
   XLAL_CHECK_NULL( (array->data = XLALMalloc(sizeof(*(array->data))*array->length)) != NULL, XLAL_ENOMEM );
   for (UINT4 ii=0; ii<length; ii++) {
      XLAL_CHECK_NULL( (array->data[ii] = XLALCreateREAL4VectorAligned(vectorLength, align)) != NULL, XLAL_EFUNC );
   }
   return array;
}
void destroyREAL4VectorAlignedArray(REAL4VectorAlignedArray *array)
{
   if (!array) return;
   for (UINT4 ii=0; ii<array->length; ii++) {
      XLALDestroyREAL4VectorAligned(array->data[ii]);
   }
   XLALFree(array->data);
   XLALFree(array);
}

INT4 DirichletRatioVector(COMPLEX8Vector *output, alignedREAL8Vector *delta0, alignedREAL8Vector *delta1, alignedREAL8Vector *scaling, const UserInput_t *params)
{

   XLAL_CHECK( output!=NULL && delta0!=NULL && delta1!=NULL && scaling!=NULL && params!=NULL, XLAL_EFUNC );

   REAL4VectorAligned *delta0_int = NULL, *delta1_int = NULL, *roundedDelta0_int = NULL, *sinPiDelta0 = NULL, *sinPiDelta1 = NULL, *PiDelta = NULL, *cosPiDelta0 = NULL, *cosPiDelta1 = NULL, *realTerms = NULL, *realTerms2 = NULL, *imagTerms = NULL, *imagTerms2 = NULL;
   XLAL_CHECK( (delta0_int = XLALCreateREAL4VectorAligned(delta0->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (delta1_int = XLALCreateREAL4VectorAligned(delta0->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (roundedDelta0_int = XLALCreateREAL4VectorAligned(delta0->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (sinPiDelta0 = XLALCreateREAL4VectorAligned(delta0->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (sinPiDelta1 = XLALCreateREAL4VectorAligned(delta0->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (cosPiDelta0 = XLALCreateREAL4VectorAligned(delta0->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (cosPiDelta1 = XLALCreateREAL4VectorAligned(delta0->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (PiDelta = XLALCreateREAL4VectorAligned(delta0->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (realTerms = XLALCreateREAL4VectorAligned(delta0->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (realTerms2 = XLALCreateREAL4VectorAligned(delta0->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (imagTerms = XLALCreateREAL4VectorAligned(delta0->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (imagTerms2 = XLALCreateREAL4VectorAligned(delta0->length, 32)) != NULL, XLAL_EFUNC );

   for (UINT4 ii=0; ii<delta0_int->length; ii++) {
      delta0_int->data[ii] = (REAL4)(delta0->data[ii]);
      delta1_int->data[ii] = (REAL4)(delta1->data[ii]);
   }

   XLAL_CHECK( VectorRoundREAL4(roundedDelta0_int, delta0_int, params->vectorMath) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( XLALVectorScaleREAL4(PiDelta->data, (REAL4)LAL_PI, delta0_int->data, delta0_int->length) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( XLALVectorSinCosREAL4(sinPiDelta0->data, cosPiDelta0->data, PiDelta->data, PiDelta->length) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( XLALVectorScaleREAL4(PiDelta->data, (REAL4)LAL_PI, delta1_int->data, delta1_int->length) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( XLALVectorSinCosREAL4(sinPiDelta1->data, cosPiDelta1->data, PiDelta->data, PiDelta->length) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( XLALVectorMultiplyREAL4(realTerms->data, cosPiDelta1->data, cosPiDelta0->data, cosPiDelta1->length) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( XLALVectorMultiplyREAL4(realTerms2->data, sinPiDelta1->data, sinPiDelta0->data, sinPiDelta1->length) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( XLALVectorAddREAL4(realTerms->data, realTerms->data, realTerms2->data, realTerms->length) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( XLALVectorMultiplyREAL4(imagTerms->data, sinPiDelta1->data, cosPiDelta0->data, sinPiDelta1->length) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( XLALVectorMultiplyREAL4(imagTerms2->data, cosPiDelta1->data, sinPiDelta0->data, cosPiDelta1->length) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( XLALVectorScaleREAL4(imagTerms2->data, -1.0, imagTerms2->data, imagTerms2->length) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( XLALVectorAddREAL4(imagTerms->data, imagTerms->data, imagTerms2->data, imagTerms->length) == XLAL_SUCCESS, XLAL_EFUNC );

   for (UINT4 ii=0; ii<output->length; ii++) {
      if (fabsf(delta1_int->data[ii])<(REAL4)1.0e-6) {
         if (fabsf(delta0_int->data[ii])<(REAL4)1.0e-6) output->data[ii] = 1.0;
         else if (fabsf((REAL4)(delta0_int->data[ii]*delta0_int->data[ii]-1.0))<(REAL4)1.0e-6) output->data[ii] = -2.0;
         else if (fabsf((REAL4)(delta0_int->data[ii]-roundedDelta0_int->data[ii]))<(REAL4)1.0e-6) output->data[ii] = 0.0;
	 else output->data[ii] = -LAL_PI*crectf(cosPiDelta0->data[ii], -sinPiDelta0->data[ii])*delta0_int->data[ii]*(delta0_int->data[ii]*delta0_int->data[ii] - 1.0)/sinPiDelta0->data[ii];
      } else if (fabsf((REAL4)(delta1_int->data[ii]*delta1_int->data[ii]-1.0))<(REAL4)1.0e-6) {
         if (fabsf(delta0_int->data[ii])<(REAL4)1.0e-6) output->data[ii] = -0.5;
         else if (fabsf((REAL4)(delta0_int->data[ii]*delta0_int->data[ii]-1.0))<(REAL4)1.0e-6) output->data[ii] = 1.0;
         else if (fabsf((REAL4)(delta0_int->data[ii]-roundedDelta0_int->data[ii]))<(REAL4)1.0e-6) output->data[ii] = 0.0;
	 else output->data[ii] = -LAL_PI_2*crectf(-cosPiDelta0->data[ii], sinPiDelta0->data[ii])*delta0_int->data[ii]*(delta0_int->data[ii]*delta0_int->data[ii] - 1.0)/sinPiDelta0->data[ii];
      } else if (fabsf(delta0_int->data[ii])<(REAL4)1.0e-6) output->data[ii] = -LAL_1_PI*crectf(cosPiDelta1->data[ii], sinPiDelta1->data[ii])*sinPiDelta1->data[ii]/(delta1_int->data[ii]*(delta1_int->data[ii]*delta1_int->data[ii]-1.0));
      else if (fabsf((REAL4)(delta0_int->data[ii]*delta0_int->data[ii] - 1.0))<(REAL4)1.0e-6) output->data[ii] = LAL_2_PI*crectf(cosPiDelta1->data[ii], sinPiDelta1->data[ii])*sinPiDelta1->data[ii]/(delta1_int->data[ii]*(delta1_int->data[ii]*delta1_int->data[ii]-1.0));
      else if (fabsf((REAL4)(delta0_int->data[ii]-roundedDelta0_int->data[ii]))<(REAL4)1.0e-6) output->data[ii] = 0.0;
      else output->data[ii] = scaling->data[ii]*sinPiDelta1->data[ii]/sinPiDelta0->data[ii]*crectf(realTerms->data[ii], imagTerms->data[ii]);
   }

   XLALDestroyREAL4VectorAligned(delta0_int);
   XLALDestroyREAL4VectorAligned(delta1_int);
   XLALDestroyREAL4VectorAligned(roundedDelta0_int);
   XLALDestroyREAL4VectorAligned(sinPiDelta0);
   XLALDestroyREAL4VectorAligned(sinPiDelta1);
   XLALDestroyREAL4VectorAligned(cosPiDelta0);
   XLALDestroyREAL4VectorAligned(cosPiDelta1);
   XLALDestroyREAL4VectorAligned(PiDelta);
   XLALDestroyREAL4VectorAligned(realTerms);
   XLALDestroyREAL4VectorAligned(realTerms2);
   XLALDestroyREAL4VectorAligned(imagTerms);
   XLALDestroyREAL4VectorAligned(imagTerms2);

   return XLAL_SUCCESS;

}

/**
 * \brief Computes a multiplication of two vectors with a stride and initial offset
 *
 * Be sure you know what you are doing or this could go wrong (no error checking for speed!)
 * \param [out] output  Pointer to REAL4VectorAligned output
 * \param [in]  input1  Pointer to first REAL4VectorAligned input
 * \param [in]  input2  Pointer to second REAL4VectorAligned input
 * \param [in]  stride1 Skip stride1 number of elements in input1
 * \param [in]  stride2 Skip stride2 number of elements in input2
 * \param [in]  offset1 Start at offset1 number of elements from the beginning of input1
 * \param [in]  offset2 Start at offset2 number of elements from the beginning of input2
 * \return Status value
 */
INT4 fastSSVectorMultiply_with_stride_and_offset(REAL4VectorAligned *output, const REAL4VectorAligned *input1, const REAL4VectorAligned *input2, const INT4 stride1, const INT4 stride2, const INT4 offset1, const INT4 offset2)
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

   return XLAL_SUCCESS;

} /* SSVectorMultiply_with_stride_and_offset() */

INT4 VectorSubtractREAL4(REAL4VectorAligned *output, REAL4VectorAligned *input1, REAL4VectorAligned *input2, INT4 vectorMath)
{
   XLAL_CHECK( output!=NULL && input1!=NULL && input2!=NULL, XLAL_EINVAL );
   if (vectorMath==1) XLAL_CHECK( sseSSVectorSubtract(output, input1, input2) == XLAL_SUCCESS, XLAL_EFUNC );
   else if (vectorMath==2) XLAL_CHECK( avxSSVectorSubtract(output, input1, input2) == XLAL_SUCCESS, XLAL_EFUNC );
   else for (UINT4 ii=0; ii<input1->length; ii++) output->data[ii] = input1->data[ii] - input2->data[ii];
   return XLAL_SUCCESS;
}

INT4 VectorScaleREAL8(alignedREAL8Vector *output, alignedREAL8Vector *input, REAL8 scale, INT4 vectorMath)
{
   XLAL_CHECK( output!=NULL && input!=NULL, XLAL_EINVAL );
   if (vectorMath==1) XLAL_CHECK( sseScaleREAL8Vector(output, input, scale) == XLAL_SUCCESS, XLAL_EFUNC );
   else if (vectorMath==2) XLAL_CHECK( avxScaleREAL8Vector(output, input, scale) == XLAL_SUCCESS, XLAL_EFUNC );
   else for (UINT4 ii=0; ii<input->length; ii++) output->data[ii] = input->data[ii] * scale;
   return XLAL_SUCCESS;
}

INT4 VectorShiftREAL8(alignedREAL8Vector *output, alignedREAL8Vector *input, REAL8 shift, INT4 vectorMath)
{
   XLAL_CHECK( output!=NULL && input!=NULL, XLAL_EINVAL );
   if (vectorMath==1) XLAL_CHECK( sseAddScalarToREAL8Vector(output, input, shift) == XLAL_SUCCESS, XLAL_EFUNC );
   else if (vectorMath==2) XLAL_CHECK( avxAddScalarToREAL8Vector(output, input, shift) == XLAL_SUCCESS, XLAL_EFUNC );
   else for (UINT4 ii=0; ii<input->length; ii++) output->data[ii] = input->data[ii] + shift;
   return XLAL_SUCCESS;
}

INT4 VectorAddREAL8(alignedREAL8Vector *output, alignedREAL8Vector *input1, alignedREAL8Vector *input2, INT4 vectorMath)
{
   XLAL_CHECK( output!=NULL && input1!=NULL && input2!=NULL, XLAL_EINVAL );
   if (vectorMath==1) XLAL_CHECK( sseDDVectorSum(output, input1, input2) == XLAL_SUCCESS, XLAL_EFUNC );
   else if (vectorMath==2) XLAL_CHECK( avxDDVectorSum(output, input1, input2) == XLAL_SUCCESS, XLAL_EFUNC );
   else for (UINT4 ii=0; ii<input1->length; ii++) output->data[ii] = input1->data[ii] + input2->data[ii];
   return XLAL_SUCCESS;
}

INT4 VectorSubtractREAL8(alignedREAL8Vector *output, alignedREAL8Vector *input1, alignedREAL8Vector *input2, INT4 vectorMath)
{
   XLAL_CHECK( output!=NULL && input1!=NULL && input2!=NULL, XLAL_EINVAL );
   if (vectorMath==1) XLAL_CHECK( sseDDVectorSubtract(output, input1, input2) == XLAL_SUCCESS, XLAL_EFUNC );
   else if (vectorMath==2) XLAL_CHECK( avxDDVectorSubtract(output, input1, input2) == XLAL_SUCCESS, XLAL_EFUNC );
   else for (UINT4 ii=0; ii<input1->length; ii++) output->data[ii] = input1->data[ii] - input2->data[ii];
   return XLAL_SUCCESS;
}

INT4 VectorMultiplyREAL8(alignedREAL8Vector *output, alignedREAL8Vector *input1, alignedREAL8Vector *input2, INT4 vectorMath)
{
   XLAL_CHECK( output!=NULL && input1!=NULL && input2!=NULL, XLAL_EINVAL );
   if (vectorMath==1) XLAL_CHECK( sseDDVectorMultiply(output, input1, input2) == XLAL_SUCCESS, XLAL_EFUNC );
   else if (vectorMath==2) XLAL_CHECK( avxDDVectorMultiply(output, input1, input2) == XLAL_SUCCESS, XLAL_EFUNC );
   else for (UINT4 ii=0; ii<input1->length; ii++) output->data[ii] = input1->data[ii] * input2->data[ii];
   return XLAL_SUCCESS;
}

INT4 VectorInvertREAL8(alignedREAL8Vector *output, alignedREAL8Vector *input, INT4 vectorMath)
{
   XLAL_CHECK( output!=NULL && input!=NULL, XLAL_EINVAL );
   if (vectorMath==1) XLAL_CHECK( sseInvertREAL8Vector(output, input) == XLAL_SUCCESS, XLAL_EFUNC );
   else if (vectorMath==2) XLAL_CHECK( avxInvertREAL8Vector(output, input) == XLAL_SUCCESS, XLAL_EFUNC );
   else for (UINT4 ii=0; ii<input->length; ii++) output->data[ii] = 1.0/input->data[ii];
   return XLAL_SUCCESS;
}

INT4 VectorFloorREAL8(alignedREAL8Vector *output, alignedREAL8Vector *input, INT4 vectorMath)
{
   XLAL_CHECK( output!=NULL && input!=NULL, XLAL_EINVAL );
   if (vectorMath==2) {
#ifdef __AVX__
      _mm256_zeroupper();
      INT4 roundedvectorlength = (INT4)input->length / 4;
      __m256d *arr, *result;
      arr = (__m256d*)(void*)input->data;
      result = (__m256d*)(void*)output->data;
      for (INT4 ii=0; ii<roundedvectorlength; ii++) {
         *result = _mm256_floor_pd(*arr);
         arr++;
         result++;
      }
      _mm256_zeroupper();
      for (UINT4 ii=4*roundedvectorlength; ii<input->length; ii++) output->data[ii] = floor(input->data[ii]);
#else
      (void)output;
      (void)input;
      fprintf(stderr, "%s: Failed because AVX is not supported, possibly because -mavx flag wasn't used for compiling.\n", __func__);
      XLAL_ERROR(XLAL_EFAILED);
#endif
   } else for (UINT4 ii=0; ii<input->length; ii++) output->data[ii] = floor(input->data[ii]);
   return XLAL_SUCCESS;
}

INT4 VectorRoundREAL4(REAL4VectorAligned *output, REAL4VectorAligned *input, INT4 vectorMath)
{
   XLAL_CHECK( output!=NULL && input!=NULL, XLAL_EINVAL );
   if (vectorMath==2) {
#ifdef __AVX__
      _mm256_zeroupper();
      INT4 roundedvectorlength = (INT4)input->length / 8;
      __m256 *arr, *result;
      arr = (__m256*)(void*)input->data;
      result = (__m256*)(void*)output->data;
      for (INT4 ii=0; ii<roundedvectorlength; ii++) {
         *result = _mm256_round_ps(*arr, _MM_FROUND_TO_NEAREST_INT);
         arr++;
         result++;
      }
      _mm256_zeroupper();
      for (UINT4 ii=8*roundedvectorlength; ii<input->length; ii++) output->data[ii] = roundf(input->data[ii]);
#else
      (void)output;
      (void)input;
      fprintf(stderr, "%s: Failed because AVX is not supported, possibly because -mavx flag wasn't used for compiling.\n", __func__);
      XLAL_ERROR(XLAL_EFAILED);
#endif
   } else for (UINT4 ii=0; ii<input->length; ii++) output->data[ii] = round(input->data[ii]);
   return XLAL_SUCCESS;
}

INT4 VectorRoundREAL8(alignedREAL8Vector *output, alignedREAL8Vector *input, INT4 vectorMath)
{
   XLAL_CHECK( output!=NULL && input!=NULL, XLAL_EINVAL );
   if (vectorMath==2) {
#ifdef __AVX__
      _mm256_zeroupper();
      INT4 roundedvectorlength = (INT4)input->length / 4;
      __m256d *arr, *result;
      arr = (__m256d*)(void*)input->data;
      result = (__m256d*)(void*)output->data;
      for (INT4 ii=0; ii<roundedvectorlength; ii++) {
         *result = _mm256_round_pd(*arr, _MM_FROUND_TO_NEAREST_INT);
         arr++;
         result++;
      }
      _mm256_zeroupper();
      for (UINT4 ii=4*roundedvectorlength; ii<input->length; ii++) output->data[ii] = round(input->data[ii]);
#else
      (void)output;
      (void)input;
      fprintf(stderr, "%s: Failed because AVX is not supported, possibly because -mavx flag wasn't used for compiling.\n", __func__);
      XLAL_ERROR(XLAL_EFAILED);
#endif
   } else for (UINT4 ii=0; ii<input->length; ii++) output->data[ii] = round(input->data[ii]);
   return XLAL_SUCCESS;
}

INT4 VectorAbsREAL4(REAL4VectorAligned *output, REAL4VectorAligned *input, INT4 vectorMath)
{
   XLAL_CHECK( output!=NULL && input!=NULL, XLAL_EINVAL );
   if (vectorMath==2) {
#ifdef __AVX__
      _mm256_zeroupper();
      INT4 roundedvectorlength = (INT4)input->length / 8;
      __m256 *arr, *result;
      __m256i mask = _mm256_set1_epi32(~0x80000000);
      arr = (__m256*)(void*)input->data;
      result = (__m256*)(void*)output->data;
      for (INT4 ii=0; ii<roundedvectorlength; ii++) {
         *result = _mm256_and_ps(*arr, (__m256)mask);
         arr++;
         result++;
      }
      _mm256_zeroupper();
      for (UINT4 ii=8*roundedvectorlength; ii<input->length; ii++) output->data[ii] = fabsf(input->data[ii]);
#else
      (void)output;
      (void)input;
      fprintf(stderr, "%s: Failed because AVX is not supported, possibly because -mavx flag wasn't used for compiling.\n", __func__);
      XLAL_ERROR(XLAL_EFAILED);
#endif
   } else if (vectorMath==1) {
#ifdef __SSE2__
      INT4 roundedvectorlength = (INT4)input->length / 4;
      __m128 *arr, *result;
      __m128i mask = _mm_set1_epi32(~0x80000000);
      arr = (__m128*)(void*)input->data;
      result = (__m128*)(void*)output->data;
      for (INT4 ii=0; ii<roundedvectorlength; ii++) {
         *result = _mm_and_ps(*arr, (__m128)mask);
         arr++;
         result++;
      }
      for (UINT4 ii=4*roundedvectorlength; ii<input->length; ii++) output->data[ii] = fabsf(input->data[ii]);
#else
      (void)output;
      (void)input;
      fprintf(stderr, "%s: Failed because SSE2 is not supported, possibly because -msse2 flag wasn't used for compiling.\n", __func__);
      XLAL_ERROR(XLAL_EFAILED);
#endif
   } else for (UINT4 ii=0; ii<input->length; ii++) output->data[ii] = fabsf(input->data[ii]);
   return XLAL_SUCCESS;
}

INT4 VectorAbsREAL8(alignedREAL8Vector *output, alignedREAL8Vector *input, INT4 vectorMath)
{
   XLAL_CHECK( output!=NULL && input!=NULL, XLAL_EINVAL );
   if (vectorMath==2) {
#ifdef __AVX__
      _mm256_zeroupper();
      INT4 roundedvectorlength = (INT4)input->length / 4;
      __m256d *arr, *result;
      __m256i mask = _mm256_set1_epi64x(~0x8000000000000000);
      arr = (__m256d*)(void*)input->data;
      result = (__m256d*)(void*)output->data;
      for (INT4 ii=0; ii<roundedvectorlength; ii++) {
         *result = _mm256_and_pd(*arr, (__m256d)mask);
         arr++;
         result++;
      }
      _mm256_zeroupper();
      for (UINT4 ii=4*roundedvectorlength; ii<input->length; ii++) output->data[ii] = fabs(input->data[ii]);
#else
      (void)output;
      (void)input;
      fprintf(stderr, "%s: Failed because AVX is not supported, possibly because -mavx flag wasn't used for compiling.\n", __func__);
      XLAL_ERROR(XLAL_EFAILED);
#endif
   } else if (vectorMath==1) {
#ifdef __SSE2__
      INT4 roundedvectorlength = (INT4)input->length / 2;
      __m128d *arr, *result;
      __m128i mask = _mm_set1_epi64x(~0x8000000000000000);
      arr = (__m128d*)(void*)input->data;
      result = (__m128d*)(void*)output->data;
      for (INT4 ii=0; ii<roundedvectorlength; ii++) {
         *result = _mm_and_pd(*arr, (__m128d)mask);
         arr++;
         result++;
      }
      for (UINT4 ii=2*roundedvectorlength; ii<input->length; ii++) output->data[ii] = fabs(input->data[ii]);
#else
      (void)output;
      (void)input;
      fprintf(stderr, "%s: Failed because SSE2 is not supported, possibly because -msse2 flag wasn't used for compiling.\n", __func__);
      XLAL_ERROR(XLAL_EFAILED);
#endif
   } else for (UINT4 ii=0; ii<input->length; ii++) output->data[ii] = fabs(input->data[ii]);
   return XLAL_SUCCESS;
}

INT4 VectorCabsfCOMPLEX8(REAL4VectorAligned *output, COMPLEX8Vector *input, INT4 vectorMath)
{
   XLAL_CHECK( output!=NULL && input!=NULL, XLAL_EINVAL );
   if (vectorMath==2) {
#ifdef __AVX__
      _mm256_zeroupper();
      INT4 roundedvectorlength = (INT4)input->length / 8;
      __m256 *result = (__m256*)(void*)output->data;
      INT4 position = 0;
      for (INT4 ii=0; ii<roundedvectorlength; ii++) {
	 __m256 input_re = _mm256_set_ps((float)crealf(input->data[position+7]), (float)crealf(input->data[position+6]), (float)crealf(input->data[position+5]), (float)crealf(input->data[position+4]), (float)crealf(input->data[position+3]), (float)crealf(input->data[position+2]), (float)crealf(input->data[position+1]), (float)crealf(input->data[position]));
	 __m256 input_im = _mm256_set_ps((float)cimagf(input->data[position+7]), (float)cimagf(input->data[position+6]), (float)cimagf(input->data[position+5]), (float)cimagf(input->data[position+4]), (float)cimagf(input->data[position+3]), (float)cimagf(input->data[position+2]), (float)cimagf(input->data[position+1]), (float)cimagf(input->data[position]));
	 __m256 input_re_sq = _mm256_mul_ps(input_re, input_re);
	 __m256 input_im_sq = _mm256_mul_ps(input_im, input_im);
         __m256 sqsum = _mm256_add_ps(input_re_sq, input_im_sq);
	 *result = _mm256_sqrt_ps(sqsum);
	 position += 8;
         result++;
      }
      _mm256_zeroupper();
      for (UINT4 ii=8*roundedvectorlength; ii<input->length; ii++) output->data[ii] = cabsf(input->data[ii]);
#else
      (void)output;
      (void)input;
      fprintf(stderr, "%s: Failed because AVX is not supported, possibly because -mavx flag wasn't used for compiling.\n", __func__);
      XLAL_ERROR(XLAL_EFAILED);
#endif
   } else if (vectorMath==1) {
#ifdef __SSE2__
      INT4 roundedvectorlength = (INT4)input->length / 4;
      __m128 *result = (__m128*)(void*)output->data;
      INT4 position = 0;
      for (INT4 ii=0; ii<roundedvectorlength; ii++) {
	 __m128 input_re = _mm_set_ps((float)crealf(input->data[position+3]), (float)crealf(input->data[position+2]), (float)crealf(input->data[position+1]), (float)crealf(input->data[position]));
	 __m128 input_im = _mm_set_ps((float)cimagf(input->data[position+3]), (float)cimagf(input->data[position+2]), (float)cimagf(input->data[position+1]), (float)cimagf(input->data[position]));
	 __m128 input_re_sq = _mm_mul_ps(input_re, input_re);
	 __m128 input_im_sq = _mm_mul_ps(input_im, input_im);
         __m128 sqsum = _mm_add_ps(input_re_sq, input_im_sq);
	 *result = _mm_sqrt_ps(sqsum);
	 position += 4;
         result++;
      }
      for (UINT4 ii=4*roundedvectorlength; ii<input->length; ii++) output->data[ii] = cabsf(input->data[ii]);
#else
      (void)output;
      (void)input;
      fprintf(stderr, "%s: Failed because SSE2 is not supported, possibly because -msse2 flag wasn't used for compiling.\n", __func__);
      XLAL_ERROR(XLAL_EFAILED);
#endif
   } else for (UINT4 ii=0; ii<input->length; ii++) output->data[ii] = cabsf(input->data[ii]);
   return XLAL_SUCCESS;
}
INT4 VectorCabsCOMPLEX8(alignedREAL8Vector *output, COMPLEX8Vector *input, INT4 vectorMath)
{
   XLAL_CHECK( output!=NULL && input!=NULL, XLAL_EINVAL );
   if (vectorMath==2) {
#ifdef __AVX__
      _mm256_zeroupper();
      INT4 roundedvectorlength = (INT4)input->length / 4;
      __m256d *result = (__m256d*)(void*)output->data;
      INT4 position = 0;
      for (INT4 ii=0; ii<roundedvectorlength; ii++) {
	 __m256d input_re = _mm256_set_pd((double)crealf(input->data[position+3]), (double)crealf(input->data[position+2]), (double)crealf(input->data[position+1]), (double)crealf(input->data[position]));
	 __m256d input_im = _mm256_set_pd((double)cimagf(input->data[position+3]), (double)cimagf(input->data[position+2]), (double)cimagf(input->data[position+1]), (double)cimagf(input->data[position]));
	 __m256d input_re_sq = _mm256_mul_pd(input_re, input_re);
	 __m256d input_im_sq = _mm256_mul_pd(input_im, input_im);
         __m256d sqsum = _mm256_add_pd(input_re_sq, input_im_sq);
	 *result = _mm256_sqrt_pd(sqsum);
	 position += 4;
         result++;
      }
      _mm256_zeroupper();
      for (UINT4 ii=4*roundedvectorlength; ii<input->length; ii++) output->data[ii] = cabs(input->data[ii]);
#else
      (void)output;
      (void)input;
      fprintf(stderr, "%s: Failed because AVX is not supported, possibly because -mavx flag wasn't used for compiling.\n", __func__);
      XLAL_ERROR(XLAL_EFAILED);
#endif
   } else if (vectorMath==1) {
#ifdef __SSE2__
      INT4 roundedvectorlength = (INT4)input->length / 2;
      __m128d *result = (__m128d*)(void*)output->data;
      INT4 position = 0;
      for (INT4 ii=0; ii<roundedvectorlength; ii++) {
	 __m128d input_re = _mm_set_pd((double)crealf(input->data[position+1]), (double)crealf(input->data[position]));
	 __m128d input_im = _mm_set_pd((double)cimagf(input->data[position+1]), (double)cimagf(input->data[position]));
	 __m128d input_re_sq = _mm_mul_pd(input_re, input_re);
	 __m128d input_im_sq = _mm_mul_pd(input_im, input_im);
         __m128d sqsum = _mm_add_pd(input_re_sq, input_im_sq);
	 *result = _mm_sqrt_pd(sqsum);
	 position += 2;
         result++;
      }
      for (UINT4 ii=2*roundedvectorlength; ii<input->length; ii++) output->data[ii] = cabs(input->data[ii]);
#else
      (void)output;
      (void)input;
      fprintf(stderr, "%s: Failed because SSE2 is not supported, possibly because -msse2 flag wasn't used for compiling.\n", __func__);
      XLAL_ERROR(XLAL_EFAILED);
#endif
   } else for (UINT4 ii=0; ii<input->length; ii++) output->data[ii] = cabs(input->data[ii]);
   return XLAL_SUCCESS;
}

/**
 * Sum two alignedREAL8Vector using SSE2
 * \param [out] output Pointer to a alignedREAL8Vector
 * \param [in]  input1 Pointer to a alignedREAL8Vector
 * \param [in]  input2 Pointer to a alignedREAL8Vector
 * \return Status value
 */
INT4 sseDDVectorSum(alignedREAL8Vector *output, alignedREAL8Vector *input1, alignedREAL8Vector *input2)
{

#ifdef __SSE2__
   INT4 roundedvectorlength = (INT4)input1->length / 2;

   __m128d *arr1, *arr2, *result;
   arr1 = (__m128d*)(void*)input1->data;
   arr2 = (__m128d*)(void*)input2->data;
   result = (__m128d*)(void*)output->data;

   //add the two vectors into the output
   for (INT4 ii=0; ii<roundedvectorlength; ii++) {
      *result = _mm_add_pd(*arr1, *arr2);
      arr1++;
      arr2++;
      result++;
   }

   //Finish up the remaining part
   for (UINT4 ii=2*roundedvectorlength; ii<input1->length; ii++) output->data[ii] = input1->data[ii] + input2->data[ii];

   return XLAL_SUCCESS;
#else
   (void)output;
   (void)input1;
   (void)input2;
   fprintf(stderr, "%s: Failed because SSE2 is not supported, possibly because -msse2 flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR(XLAL_EFAILED);
#endif

}


/**
 * Sum two alignedREAL8Vector using AVX
 * \param [out] output Pointer to a alignedREAL8Vector
 * \param [in]  input1 Pointer to a alignedREAL8Vector
 * \param [in]  input2 Pointer to a alignedREAL8Vector
 * \return Status value
 */
INT4 avxDDVectorSum(alignedREAL8Vector *output, alignedREAL8Vector *input1, alignedREAL8Vector *input2)
{

#ifdef __AVX__
   //Need to zero the upper 128 bits in case of an SSE function call after this AVX function
   _mm256_zeroupper();
   
   INT4 roundedvectorlength = (INT4)input1->length / 4;

   __m256d *arr1, *arr2, *result;
   arr1 = (__m256d*)(void*)input1->data;
   arr2 = (__m256d*)(void*)input2->data;
   result = (__m256d*)(void*)output->data;

   //add the two vectors into the output
   for (INT4 ii=0; ii<roundedvectorlength; ii++) {
      *result = _mm256_add_pd(*arr1, *arr2);
      arr1++;
      arr2++;
      result++;
   }

   //Finish up the remaining part
   for (UINT4 ii=4*roundedvectorlength; ii<input1->length; ii++) output->data[ii] = input1->data[ii] + input2->data[ii];

   //Need to zero the upper 128 bits in case of an SSE function call after this AVX function
   _mm256_zeroupper();

   return XLAL_SUCCESS;
#else
   (void)output;
   (void)input1;
   (void)input2;
   fprintf(stderr, "%s: Failed because AVX is not supported, possibly because -mavx flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR(XLAL_EFAILED);
#endif

}

/**
 * Subtract two REAL4VectorAligned using SSE
 * \param [out] output Pointer to a REAL4VectorAligned
 * \param [in]  input1 Pointer to a REAL4VectorAligned
 * \param [in]  input2 Pointer to a REAL4VectorAligned
 * \return Status value
 */
INT4 sseSSVectorSubtract(REAL4VectorAligned *output, REAL4VectorAligned *input1, REAL4VectorAligned *input2)
{

#ifdef __SSE__
   INT4 roundedvectorlength = (INT4)input1->length / 4;

   __m128 *arr1, *arr2, *result;
   arr1 = (__m128*)(void*)input1->data;
   arr2 = (__m128*)(void*)input2->data;
   result = (__m128*)(void*)output->data;

   //subtract the two vectors into the output
   for (INT4 ii=0; ii<roundedvectorlength; ii++) {
      *result = _mm_sub_ps(*arr1, *arr2);
      arr1++;
      arr2++;
      result++;
   }

   //Finish up the remaining part
   for (UINT4 ii=4*roundedvectorlength; ii<input1->length; ii++) output->data[ii] = input1->data[ii] - input2->data[ii];

   return XLAL_SUCCESS;
#else
   (void)output;
   (void)input1;
   (void)input2;
   fprintf(stderr, "%s: Failed because SSE is not supported, possibly because -msse flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR(XLAL_EFAILED);
#endif

}

/**
 * Subtract two REAL4VectorAligned using AVX
 * \param [out] output Pointer to a REAL4VectorAligned
 * \param [in]  input1 Pointer to a REAL4VectorAligned
 * \param [in]  input2 Pointer to a REAL4VectorAligned
 * \return Status value
 */
INT4 avxSSVectorSubtract(REAL4VectorAligned *output, REAL4VectorAligned *input1, REAL4VectorAligned *input2)
{

#ifdef __AVX__
   //Need to zero the upper 128 bits in case of an SSE function call after this AVX function
   _mm256_zeroupper();

   INT4 roundedvectorlength = (INT4)input1->length / 8;

   __m256 *arr1, *arr2, *result;
   arr1 = (__m256*)(void*)input1->data;
   arr2 = (__m256*)(void*)input2->data;
   result = (__m256*)(void*)output->data;

   //subtract the two vectors into the output
   for (INT4 ii=0; ii<roundedvectorlength; ii++) {
      *result = _mm256_sub_ps(*arr1, *arr2);
      arr1++;
      arr2++;
      result++;
   }

   //Finish up the remaining part
   for (UINT4 ii=8*roundedvectorlength; ii<input1->length; ii++) output->data[ii] = input1->data[ii] - input2->data[ii];

   //Need to zero the upper 128 bits in case of an SSE function call after this AVX function
   _mm256_zeroupper();

   return XLAL_SUCCESS;
#else
   (void)output;
   (void)input1;
   (void)input2;
   fprintf(stderr, "%s: Failed because AVX is not supported, possibly because -mavx flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR(XLAL_EFAILED);
#endif

}

/**
 * Subtract two alignedREAL8Vector using SSE
 * \param [out] output Pointer to a alignedREAL8Vector
 * \param [in]  input1 Pointer to a alignedREAL8Vector
 * \param [in]  input2 Pointer to a alignedREAL8Vector
 * \return Status value
 */
INT4 sseDDVectorSubtract(alignedREAL8Vector *output, alignedREAL8Vector *input1, alignedREAL8Vector *input2)
{

#ifdef __SSE2__
   INT4 roundedvectorlength = (INT4)input1->length / 2;

   __m128d *arr1, *arr2, *result;
   arr1 = (__m128d*)(void*)input1->data;
   arr2 = (__m128d*)(void*)input2->data;
   result = (__m128d*)(void*)output->data;

   //Subtract the two vectors into the output
   for (INT4 ii=0; ii<roundedvectorlength; ii++) {
      *result = _mm_sub_pd(*arr1, *arr2);
      arr1++;
      arr2++;
      result++;
   }

   //Finish up the remaining part
   for (UINT4 ii=2*roundedvectorlength; ii<input1->length; ii++) output->data[ii] = input1->data[ii] - input2->data[ii];

   return XLAL_SUCCESS;
#else
   (void)output;
   (void)input1;
   (void)input2;
   fprintf(stderr, "%s: Failed because SSE2 is not supported, possibly because -msse2 flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR(XLAL_EFAILED);
#endif

}

/**
 * Subtract two alignedREAL8Vector using AVX
 * \param [out] output Pointer to a alignedREAL8Vector
 * \param [in]  input1 Pointer to a alignedREAL8Vector
 * \param [in]  input2 Pointer to a alignedREAL8Vector
 * \return Status value
 */
INT4 avxDDVectorSubtract(alignedREAL8Vector *output, alignedREAL8Vector *input1, alignedREAL8Vector *input2)
{

#ifdef __AVX__
   //Need to zero the upper 128 bits in case of an SSE function call after this AVX function
   _mm256_zeroupper();

   INT4 roundedvectorlength = (INT4)input1->length / 4;

   __m256d *arr1, *arr2, *result;
   arr1 = (__m256d*)(void*)input1->data;
   arr2 = (__m256d*)(void*)input2->data;
   result = (__m256d*)(void*)output->data;

   //Subtract the two vectors into the output
   for (INT4 ii=0; ii<roundedvectorlength; ii++) {
      *result = _mm256_sub_pd(*arr1, *arr2);
      arr1++;
      arr2++;
      result++;
   }

   //Finish up the remaining part
   for (UINT4 ii=4*roundedvectorlength; ii<input1->length; ii++) output->data[ii] = input1->data[ii] - input2->data[ii];

   //Need to zero the upper 128 bits in case of an SSE function call after this AVX function
   _mm256_zeroupper();

   return XLAL_SUCCESS;
#else
   (void)output;
   (void)input1;
   (void)input2;
   fprintf(stderr, "%s: Failed because AVX is not supported, possibly because -mavx flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR(XLAL_EFAILED);
#endif

}

/**
 * Multiply two alignedREAL8Vector using SSE
 * \param [out] output Pointer to a alignedREAL8Vector
 * \param [in]  input1 Pointer to a alignedREAL8Vector
 * \param [in]  input2 Pointer to a alignedREAL8Vector
 * \return Status value
 */
INT4 sseDDVectorMultiply(alignedREAL8Vector *output, alignedREAL8Vector *input1, alignedREAL8Vector *input2)
{

#ifdef __SSE2__
   INT4 roundedvectorlength = (INT4)input1->length / 2;

   __m128d *arr1, *arr2, *result;
   arr1 = (__m128d*)(void*)input1->data;
   arr2 = (__m128d*)(void*)input2->data;
   result = (__m128d*)(void*)output->data;

   //multiply the two vectors into the output
   for (INT4 ii=0; ii<roundedvectorlength; ii++) {
      *result = _mm_mul_pd(*arr1, *arr2);
      arr1++;
      arr2++;
      result++;
   }

   //Finish up the remaining part
   for (INT4 ii=2*roundedvectorlength; ii<(INT4)input1->length; ii++) output->data[ii] = input1->data[ii] * input2->data[ii];

   return XLAL_SUCCESS;
#else
   (void)output;
   (void)input1;
   (void)input2;
   fprintf(stderr, "%s: Failed because SSE2 is not supported, possibly because -msse2 flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR(XLAL_EFAILED);
#endif

}


/**
 * Multiply two alignedREAL8Vector using AVX
 * \param [out] output Pointer to a alignedREAL8Vector
 * \param [in]  input1 Pointer to a alignedREAL8Vector
 * \param [in]  input2 Pointer to a alignedREAL8Vector
 * \return Status value
 */
INT4 avxDDVectorMultiply(alignedREAL8Vector *output, alignedREAL8Vector *input1, alignedREAL8Vector *input2)
{

#ifdef __AVX__
   //Need to zero the upper 128 bits in case of an SSE function call after this AVX function
   _mm256_zeroupper();

   INT4 roundedvectorlength = (INT4)input1->length / 4;

   __m256d *arr1, *arr2, *result;
   arr1 = (__m256d*)(void*)input1->data;
   arr2 = (__m256d*)(void*)input2->data;
   result = (__m256d*)(void*)output->data;

   //multiply the two vectors into the output
   for (INT4 ii=0; ii<roundedvectorlength; ii++) {
      *result = _mm256_mul_pd(*arr1, *arr2);
      arr1++;
      arr2++;
      result++;
   }

   //Finish up the remaining part
   for (INT4 ii=4*roundedvectorlength; ii<(INT4)input1->length; ii++) output->data[ii] = input1->data[ii] * input2->data[ii];

   //Need to zero the upper 128 bits in case of an SSE function call after this AVX function
   _mm256_zeroupper();

   return XLAL_SUCCESS;
#else
   (void)output;
   (void)input1;
   (void)input2;
   fprintf(stderr, "%s: Failed because AVX is not supported, possibly because -mavx flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR(XLAL_EFAILED);
#endif

}

/**
 * Invert a alignedREAL8Vector using SSE
 * \param [out] output Pointer to a alignedREAL8Vector
 * \param [in]  input1 Pointer to a alignedREAL8Vector
 * \return Status value
 */
INT4 sseInvertREAL8Vector(alignedREAL8Vector *output, alignedREAL8Vector *input1)
{

#ifdef __SSE__
   INT4 roundedvectorlength = (INT4)input1->length / 2;

   __m128d *arr1, *result;
   arr1 = (__m128d*)(void*)input1->data;
   result = (__m128d*)(void*)output->data;

   __m128d one = _mm_set1_pd(1.0);

   //Invert the vector
   for (INT4 ii=0; ii<roundedvectorlength; ii++) {
      *result = _mm_div_pd(one, *arr1);
      arr1++;
      result++;
   }

   //Finish up the remaining part
   for (INT4 ii=2*roundedvectorlength; ii<(INT4)input1->length; ii++) output->data[ii] = 1.0/input1->data[ii];

   return XLAL_SUCCESS;
#else
   (void)output;
   (void)input1;
   fprintf(stderr, "%s: Failed because SSE is not supported, possibly because -msse flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR(XLAL_EFAILED);
#endif

}


/**
 * Invert a alignedREAL8Vector using AVX
 * \param [out] output Pointer to a alignedREAL8Vector
 * \param [in]  input1 Pointer to a alignedREAL8Vector
 * \return Status value
 */
INT4 avxInvertREAL8Vector(alignedREAL8Vector *output, alignedREAL8Vector *input1)
{

#ifdef __AVX__
   //Need to zero the upper 128 bits in case of an SSE function call after this AVX function
   _mm256_zeroupper();

   INT4 roundedvectorlength = (INT4)input1->length / 4;

   __m256d *arr1, *result;
   arr1 = (__m256d*)(void*)input1->data;
   result = (__m256d*)(void*)output->data;

   __m256d one = _mm256_set1_pd(1.0);

   //Invert the vector
   for (INT4 ii=0; ii<roundedvectorlength; ii++) {
      *result = _mm256_div_pd(one, *arr1);
      arr1++;
      result++;
   }

   //Finish up the remaining part
   for (INT4 ii=4*roundedvectorlength; ii<(INT4)input1->length; ii++) output->data[ii] = 1.0/input1->data[ii];

   //Need to zero the upper 128 bits in case of an SSE function call after this AVX function
   _mm256_zeroupper();

   return XLAL_SUCCESS;
#else
   (void)output;
   (void)input1;
   fprintf(stderr, "%s: Failed because AVX is not supported, possibly because -mavx flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR(XLAL_EFAILED);
#endif

}

/**
 * Add a REAL8 scalar value to the elements of a alignedREAL8Vector using SSE
 * \param [out] output Pointer to a alignedREAL8Vector
 * \param [in]  input  Pointer to a alignedREAL8Vector
 * \param [in]  scalar Value to add to the elements of input
 * \return Status value
 */
INT4 sseAddScalarToREAL8Vector(alignedREAL8Vector *output, alignedREAL8Vector *input, REAL8 scalar)
{

#ifdef __SSE2__
   INT4 roundedvectorlength = (INT4)input->length / 2;

   __m128d *arr1, *result;
   arr1 = (__m128d*)(void*)input->data;
   result = (__m128d*)(void*)output->data;

   __m128d scalefactor = _mm_set1_pd(scalar);

   //Add the value to the vector and put in the output
   for (INT4 ii=0; ii<roundedvectorlength; ii++) {
      *result = _mm_add_pd(*arr1, scalefactor);
      arr1++;
      result++;
   }

   //Finish up the remaining part
   for (INT4 ii=2*roundedvectorlength; ii<(INT4)input->length; ii++) output->data[ii] = input->data[ii] + scalar;

   return XLAL_SUCCESS;
#else
   (void)output;
   (void)input;
   (void)scalar;
   fprintf(stderr, "%s: Failed because SSE2 is not supported, possibly because -msse2 flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR(XLAL_EFAILED);
#endif

}


/**
 * Add a REAL8 scalar value to the elements of a alignedREAL8Vector using AVX
 * \param [out] output Pointer to a alignedREAL8Vector
 * \param [in]  input  Pointer to a alignedREAL8Vector
 * \param [in]  scalar Value to add to the elements of input
 * \return Status value
 */
INT4 avxAddScalarToREAL8Vector(alignedREAL8Vector *output, alignedREAL8Vector *input, REAL8 scalar)
{

#ifdef __AVX__
   //Need to zero the upper 128 bits in case of an SSE function call after this AVX function
   _mm256_zeroupper();
   
   INT4 roundedvectorlength = (INT4)input->length / 4;

   __m256d *arr1, *result;
   arr1 = (__m256d*)(void*)input->data;
   result = (__m256d*)(void*)output->data;

   __m256d scalefactor = _mm256_set1_pd(scalar);

   //Add the value to the vector and put in the output
   for (INT4 ii=0; ii<roundedvectorlength; ii++) {
      *result = _mm256_add_pd(*arr1, scalefactor);
      arr1++;
      result++;
   }

   //Finish up the remaining part
   for (INT4 ii=4*roundedvectorlength; ii<(INT4)input->length; ii++) output->data[ii] = input->data[ii] + scalar;

   //Need to zero the upper 128 bits in case of an SSE function call after this AVX function
   _mm256_zeroupper();

   return XLAL_SUCCESS;
#else
   (void)output;
   (void)input;
   (void)scalar;
   fprintf(stderr, "%s: Failed because AVX is not supported, possibly because -mavx flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR(XLAL_EFAILED);
#endif

}

/**
 * Scale the elements of a alignedREAL8Vector by a REAL8 value using SSE
 * \param [out] output Pointer to a alignedREAL8Vector
 * \param [in]  input  Pointer to a alignedREAL8Vector
 * \param [in]  scale  Value to scale the elements of input
 * \return Status value
 */
INT4 sseScaleREAL8Vector(alignedREAL8Vector *output, alignedREAL8Vector *input, REAL8 scale)
{

#ifdef __SSE2__
   INT4 roundedvectorlength = (INT4)input->length / 2;

   __m128d *arr1, *result;
   arr1 = (__m128d*)(void*)input->data;
   result = (__m128d*)(void*)output->data;

   __m128d scalefactor = _mm_set1_pd(scale);

   //multiply the vector into the output
   for (INT4 ii=0; ii<roundedvectorlength; ii++) {
      *result = _mm_mul_pd(*arr1, scalefactor);
      arr1++;
      result++;
   }

   //Finish up the remaining part
   for (UINT4 ii=2*roundedvectorlength; ii<input->length; ii++) output->data[ii] = input->data[ii] * scale;

   return XLAL_SUCCESS;
#else
   (void)output;
   (void)input;
   (void)scale;
   fprintf(stderr, "%s: Failed because SSE2 is not supported, possibly because -msse2 flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR(XLAL_EFAILED);
#endif

}


/**
 * Scale the elements of a alignedREAL8Vector by a REAL8 value using AVX
 * \param [out] output Pointer to a alignedREAL8Vector
 * \param [in]  input  Pointer to a alignedREAL8Vector
 * \param [in]  scale  Value to scale the elements of input
 * \return Status value
 */
INT4 avxScaleREAL8Vector(alignedREAL8Vector *output, alignedREAL8Vector *input, REAL8 scale)
{

#ifdef __AVX__
   //Need to zero the upper 128 bits in case of an SSE function call after this AVX function
   _mm256_zeroupper();
   
   INT4 roundedvectorlength = (INT4)input->length / 4;

   __m256d *arr1, *result;
   arr1 = (__m256d*)(void*)input->data;
   result = (__m256d*)(void*)output->data;

   __m256d scalefactor = _mm256_set1_pd(scale);

   //multiply the vector into the output
   for (INT4 ii=0; ii<roundedvectorlength; ii++) {
      *result = _mm256_mul_pd(*arr1, scalefactor);
      arr1++;
      result++;
   }

   //Finish up the remaining part
   for (UINT4 ii=4*roundedvectorlength; ii<input->length; ii++) output->data[ii] = input->data[ii] * scale;

   //Need to zero the upper 128 bits in case of an SSE function call after this AVX function
   _mm256_zeroupper();

   return XLAL_SUCCESS;
#else
   (void)output;
   (void)input;
   (void)scale;
   fprintf(stderr, "%s: Failed because AVX is not supported, possibly because -mavx flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR(XLAL_EFAILED);
#endif

}

/**
 * Sum vectors from REAL4VectorAlignedArrays into an output REAL4VectorAlignedArray using SIMD
 * \param [out] output          Pointer to REAL4VectorAlignedArray
 * \param [in]  input1          Pointer to REAL4VectorAlignedArray
 * \param [in]  input2          Pointer to REAL4VectorAlignedArray
 * \param [in]  vectorpos1      Starting vector index for input1
 * \param [in]  vectorpos2      Starting vector index for input2
 * \param [in]  outputvectorpos Starting vector index for output
 * \param [in]  numvectors      Number of vectors to sum, incrementing vectorpos1, vectorpos2, and outputvectorpos by 1 each time
 * \return Status value
 */
INT4 VectorArraySum(REAL4VectorAlignedArray *output, REAL4VectorAlignedArray *input1, REAL4VectorAlignedArray *input2, INT4 vectorpos1, INT4 vectorpos2, INT4 outputvectorpos, INT4 numvectors)
{
   INT4 vec1 = vectorpos1, vec2 = vectorpos2, outvec = outputvectorpos;
   for (INT4 ii=0; ii<numvectors; ii++) {
      XLAL_CHECK( XLALVectorAddREAL4(output->data[outvec]->data, input1->data[vec1]->data, input2->data[vec2]->data, input1->data[vec1]->length) == XLAL_SUCCESS, XLAL_EFUNC );
      vec1++;
      vec2++;
      outvec++;
   }
   return XLAL_SUCCESS;
}

/*
 * (Deprecated) Compute from a look up table, the sin and cos of a vector of x values using SSE
 * Cannot use this for x values less than 0 or greater than 2.147483648e9
 * \param [out] sin2pix_vector Pointer to REAL8Vector of sin(2*pi*x)
 * \param [out] cos2pix_vector Pointer to REAL8Vector of cos(2*pi*x)
 * \param [in]  x              Pointer to REAL8Vector
 * \return Status value
 */
/* INT4 sse_sin_cos_2PI_LUT_REAL8Vector(REAL8Vector *sin2pix_vector, REAL8Vector *cos2pix_vector, REAL8Vector *x)
{

#ifdef __SSE2__
   INT4 roundedvectorlength = (INT4)x->length / 2;

   static BOOLEAN firstCall = TRUE;
   static REAL8 sinVal[LUT_RES+1], cosVal[LUT_RES+1];

   // the first time we get called, we set up the lookup-table
   if ( firstCall ) {
      for (UINT4 k=0; k <= LUT_RES; k++) {
         sinVal[k] = sin( LAL_TWOPI * k * OO_LUT_RES );
         cosVal[k] = cos( LAL_TWOPI * k * OO_LUT_RES );
      }
      firstCall = FALSE;
   }

   REAL8 *allocinput = NULL, *allocoutput1 = NULL, *allocoutput2 = NULL, *alignedinput = NULL, *alignedoutput1 = NULL, *alignedoutput2 = NULL;
   INT4 vecaligned = 0, outputaligned1 = 0, outputaligned2 = 0;
   __m128d *arr1, *sinresult, *cosresult;

   //Allocate memory for aligning input vector if necessary
   if ( x->data==(void*)(((uintptr_t)x->data+15) & ~15) ) {
      vecaligned = 1;
      arr1 = (__m128d*)(void*)x->data;
   } else {
      XLAL_CHECK( (allocinput = (REAL8*)XLALMalloc(2*roundedvectorlength*sizeof(REAL8) + 15)) != NULL, XLAL_ENOMEM );
      alignedinput = (void*)(((uintptr_t)allocinput+15) & ~15);
      memcpy(alignedinput, x->data, sizeof(REAL8)*2*roundedvectorlength);
      arr1 = (__m128d*)(void*)alignedinput;
   }

   //Allocate memory for aligning output vector 1 if necessary
   if ( sin2pix_vector->data==(void*)(((uintptr_t)sin2pix_vector->data+15) & ~15) ) {
      outputaligned1 = 1;
      sinresult = (__m128d*)(void*)sin2pix_vector->data;
   } else {
      XLAL_CHECK( (allocoutput1 = (REAL8*)XLALMalloc(2*roundedvectorlength*sizeof(REAL8) + 15)) != NULL, XLAL_ENOMEM );
      alignedoutput1 = (void*)(((uintptr_t)allocoutput1+15) & ~15);
      sinresult = (__m128d*)(void*)alignedoutput1;
   }

   //Allocate memory for aligning output vector 2 if necessary
   if ( cos2pix_vector->data==(void*)(((uintptr_t)cos2pix_vector->data+15) & ~15) ) {
      outputaligned2 = 1;
      cosresult = (__m128d*)(void*)cos2pix_vector->data;
   } else {
      XLAL_CHECK( (allocoutput2 = (REAL8*)XLALMalloc(2*roundedvectorlength*sizeof(REAL8) + 15)) != NULL, XLAL_ENOMEM );
      alignedoutput2 = (void*)(((uintptr_t)allocoutput2+15) & ~15);
      cosresult = (__m128d*)(void*)alignedoutput2;
   }

   INT4 *integerVect = NULL;
   XLAL_CHECK( (integerVect = (INT4*)XLALMalloc(2*sizeof(INT4)+15)) != NULL, XLAL_ENOMEM );
   INT4 *I0 = (void*)(((uintptr_t)integerVect+15) & ~15);
   memset(I0, 0, sizeof(INT4)*2);

   __m128d lutresf = _mm_set1_pd(LUT_RES_F);
   __m128d onehalf = _mm_set1_pd(0.5);
   __m128d oolutres = _mm_set1_pd(OO_LUT_RES);
   __m128d twopi = _mm_set1_pd(LAL_TWOPI);

   for (INT4 ii=0; ii<roundedvectorlength; ii++) {

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
   REAL4 sin2pix = 0.0, cos2pix = 0.0;
   for (INT4 ii=2*roundedvectorlength; ii<(INT4)x->length; ii++) {
      XLAL_CHECK( XLALSinCos2PiLUT(&sin2pix, &cos2pix, x->data[ii]) == XLAL_SUCCESS, XLAL_EFUNC );
      if (sin2pix>1.0) sin2pix = 1.0;
      else if (sin2pix<-1.0) sin2pix = -1.0;
      if (cos2pix>1.0) cos2pix = 1.0;
      else if (cos2pix<-1.0) cos2pix = -1.0;
      sin2pix_vector->data[ii] = (REAL8)sin2pix;
      cos2pix_vector->data[ii] = (REAL8)cos2pix;
   }

   //Free memory if necessary
   if (!vecaligned) XLALFree(allocinput);
   if (!outputaligned1) XLALFree(allocoutput1);
   if (!outputaligned2) XLALFree(allocoutput2);

   //Free memory
   XLALFree(integerVect);

   return XLAL_SUCCESS;
#else
   (void)sin2pix_vector;
   (void)cos2pix_vector;
   (void)x;
   fprintf(stderr, "%s: Failed because SSE2 is not supported, possibly because -msse2 flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR(XLAL_EFAILED);
#endif

} */ // sse_sin_cos_2PI_LUT_REAL8Vector()


/*
 * (Deprecated) Compute from a look up table, the sin and cos of a vector of x values using SSE
 * Cannot use this for x values less than 0 or greater than 2.147483648e9
 * \param [out] sin2pix_vector Pointer to REAL4Vector of sin(2*pi*x)
 * \param [out] cos2pix_vector Pointer to REAL4Vector of cos(2*pi*x)
 * \param [in]  x              Pointer to REAL4Vector
 * \return Status value
 */
/* INT4 sse_sin_cos_2PI_LUT_REAL4Vector(REAL4Vector *sin2pix_vector, REAL4Vector *cos2pix_vector, REAL4Vector *x)
{

#ifdef __SSE2__
   INT4 roundedvectorlength = (INT4)x->length / 4;

   static BOOLEAN firstCall = TRUE;
   static REAL4 sinVal[LUT_RES+1], cosVal[LUT_RES+1];

   // the first time we get called, we set up the lookup-table
   if ( firstCall ) {
      for (UINT4 k=0; k <= LUT_RES; k++) {
         sinVal[k] = (REAL4)sinf( (REAL4)(LAL_TWOPI * k * OO_LUT_RES) );
         cosVal[k] = (REAL4)cosf( (REAL4)(LAL_TWOPI * k * OO_LUT_RES) );
      }
      firstCall = FALSE;
   }

   REAL4 *allocinput = NULL, *allocoutput1 = NULL, *allocoutput2 = NULL, *alignedinput = NULL, *alignedoutput1 = NULL, *alignedoutput2 = NULL;
   INT4 vecaligned = 0, outputaligned1 = 0, outputaligned2 = 0;
   __m128 *arr1, *sinresult, *cosresult;

   //Allocate memory for aligning input vector if necessary
   if ( x->data==(void*)(((uintptr_t)x->data+15) & ~15) ) {
      vecaligned = 1;
      arr1 = (__m128*)(void*)x->data;
   } else {
      XLAL_CHECK( (allocinput = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15)) != NULL, XLAL_ENOMEM );
      alignedinput = (void*)(((uintptr_t)allocinput+15) & ~15);
      memcpy(alignedinput, x->data, sizeof(REAL4)*4*roundedvectorlength);
      arr1 = (__m128*)(void*)alignedinput;
   }

   //Allocate memory for aligning output vector 1 if necessary
   if ( sin2pix_vector->data==(void*)(((uintptr_t)sin2pix_vector->data+15) & ~15) ) {
      outputaligned1 = 1;
      sinresult = (__m128*)(void*)sin2pix_vector->data;
   } else {
      XLAL_CHECK( (allocoutput1 = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15)) != NULL, XLAL_ENOMEM );
      alignedoutput1 = (void*)(((uintptr_t)allocoutput1+15) & ~15);
      sinresult = (__m128*)(void*)alignedoutput1;
   }

   //Allocate memory for aligning output vector 2 if necessary
   if ( cos2pix_vector->data==(void*)(((uintptr_t)cos2pix_vector->data+15) & ~15) ) {
      outputaligned2 = 1;
      cosresult = (__m128*)(void*)cos2pix_vector->data;
   } else {
      XLAL_CHECK( (allocoutput2 = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15)) != NULL, XLAL_ENOMEM );
      alignedoutput2 = (void*)(((uintptr_t)allocoutput2+15) & ~15);
      cosresult = (__m128*)(void*)alignedoutput2;
   }

   INT4 *integerVect = NULL;
   XLAL_CHECK( (integerVect = (INT4*)XLALMalloc(4*sizeof(INT4)+15)) != NULL, XLAL_ENOMEM );
   INT4 *I0 = (void*)(((uintptr_t)integerVect+15) & ~15);
   memset(I0, 0, sizeof(INT4)*4);

   __m128 lutresf = _mm_set1_ps((REAL4)LUT_RES_F);
   __m128 onehalf = _mm_set1_ps(0.5f);
   __m128 oolutres = _mm_set1_ps((REAL4)OO_LUT_RES);
   __m128 twopi = _mm_set1_ps((REAL4)LAL_TWOPI);

   for (INT4 ii=0; ii<roundedvectorlength; ii++) {

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
   REAL4 sin2pix = 0.0, cos2pix = 0.0;
   for (INT4 ii=4*roundedvectorlength; ii<(INT4)x->length; ii++) {
      XLAL_CHECK( XLALSinCos2PiLUT(&(sin2pix), &(cos2pix), x->data[ii]) == XLAL_SUCCESS, XLAL_EFUNC );
      if (sin2pix>1.0) sin2pix = 1.0;
      else if (sin2pix<-1.0) sin2pix = -1.0;
      if (cos2pix>1.0) cos2pix = 1.0;
      else if (cos2pix<-1.0) cos2pix = -1.0;
      sin2pix_vector->data[ii] = sin2pix;
      cos2pix_vector->data[ii] = cos2pix;
   }

   //Free memory if necessary
   if (!vecaligned) XLALFree(allocinput);
   if (!outputaligned1) XLALFree(allocoutput1);
   if (!outputaligned2) XLALFree(allocoutput2);

   //Free memory
   XLALFree(integerVect);

   return XLAL_SUCCESS;
#else
   (void)sin2pix_vector;
   (void)cos2pix_vector;
   (void)x;
   fprintf(stderr, "%s: Failed because SSE2 is not supported, possibly because -msse2 flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR(XLAL_EFAILED);
#endif

} */ // sse_sin_cos_2PI_LUT_REAL4Vector()


/**
 * Exponential of input vector is computed using SSE, based on the Cephes library
 * \param [out] output Pointer to alignedREAL8Vector
 * \param [in]  input  Pointer to alignedREAL8Vector
 * \return Status value
 */
INT4 sse_exp_REAL8Vector(alignedREAL8Vector *output, alignedREAL8Vector *input)
{

#ifdef __SSE2__
   INT4 roundedvectorlength = (INT4)input->length / 2;

   __m128d *x, *result;
   x = (__m128d*)(void*)input->data;
   result = (__m128d*)(void*)output->data;

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

   for (INT4 ii=0; ii<roundedvectorlength; ii++) {
      __m128d y = *x;
      y = _mm_max_pd(y, minexp);
      y = _mm_min_pd(y, maxexp);

      //Compute x*ln(2)
      __m128d log2etimesx = _mm_mul_pd(log2e, y);

      //Round
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

   //Finish up the remaining part
   for (INT4 ii=2*roundedvectorlength; ii<(INT4)input->length; ii++) {
      output->data[ii] = exp(input->data[ii]);
   }

   return XLAL_SUCCESS;
#else
   (void)output;
   (void)input;
   fprintf(stderr, "%s: Failed because SSE2 is not supported, possibly because -msse2 flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR(XLAL_EFAILED);
#endif

}


/*
 * (Deprecated) Exponential of input vector is computed using SSE, based on the Cephes library
 * \param [out] output Pointer to REAL4Vector
 * \param [in]  input  Pointer to REAL4Vector
 * \return Status value
 */
/* INT4 sse_exp_REAL4Vector(REAL4Vector *output, REAL4Vector *input)
{

#ifdef __SSE2__
   INT4 roundedvectorlength = (INT4)input->length / 4;

   REAL4 *allocinput = NULL, *allocoutput = NULL, *alignedinput = NULL, *alignedoutput = NULL;
   INT4 vecaligned = 0, outputaligned = 0;

   __m128 *x, *result;

   //Allocate memory for aligning input vector if necessary
   if ( input->data==(void*)(((size_t)input->data+15) & ~15) ) {
      vecaligned = 1;
      x = (__m128*)(void*)input->data;
   } else {
      XLAL_CHECK( (allocinput = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15)) != NULL, XLAL_ENOMEM );
      alignedinput = (void*)(((size_t)allocinput+15) & ~15);
      memcpy(alignedinput, input->data, sizeof(REAL4)*4*roundedvectorlength);
      x = (__m128*)(void*)alignedinput;
   }

   //Allocate memory for aligning output vector 1 if necessary
   if ( output->data==(void*)(((size_t)output->data+15) & ~15) ) {
      outputaligned = 1;
      result = (__m128*)(void*)output->data;
   } else {
      XLAL_CHECK( (allocoutput = (REAL4*)XLALMalloc(4*roundedvectorlength*sizeof(REAL4) + 15)) != NULL, XLAL_ENOMEM );
      alignedoutput = (void*)(((size_t)allocoutput+15) & ~15);
      result = (__m128*)(void*)alignedoutput;
   }

   __m128i expoffset = _mm_set1_epi32(127);      //Exponent mask for single precision
   __m128 log2e = _mm_set1_ps(1.442695040888963f);                   //ln(2)
   __m128 onehalf = _mm_set1_ps(0.5f);             //0.5
   __m128 one = _mm_set1_ps(1.0f);                 //1.0
   __m128 two = _mm_set1_ps(2.0f);                 //2.0
   __m128 maxexp = _mm_set1_ps(88.3762626647949f);
   __m128 minexp = _mm_set1_ps(-88.3762626647949f);

   __m128 cephes_c1 = _mm_set1_ps(6.93145751953125e-1f);
   __m128 cephes_c2 = _mm_set1_ps(1.42860682030941723212e-6f);
   __m128 cephes_p0 = _mm_set1_ps(1.26177193074810590878e-4f);
   __m128 cephes_p1 = _mm_set1_ps(3.02994407707441961300e-2f);
   __m128 cephes_p2 = _mm_set1_ps(9.99999999999999999910e-1f);
   __m128 cephes_q0 = _mm_set1_ps(3.00198505138664455042e-6f);
   __m128 cephes_q1 = _mm_set1_ps(2.52448340349684104192e-3f);
   __m128 cephes_q2 = _mm_set1_ps(2.27265548208155028766e-1f);
   __m128 cephes_q3 = _mm_set1_ps(2.00000000000000000009e0f);

   for (INT4 ii=0; ii<roundedvectorlength; ii++) {
      __m128 y = *x;
      y = _mm_max_ps(y, minexp);
      y = _mm_min_ps(y, maxexp);

      //Compute x*ln(2)
      __m128 log2etimesx = _mm_mul_ps(log2e, y);

      //Round
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

   //Alternative method below. Similar output and errors w.r.t. libm exp()
   //__m128 cephes_c1 = _mm_set1_ps(0.693359375);
   //__m128 cephes_c2 = _mm_set1_ps(-2.12194440e-4);
   //__m128 cephes_p0 = _mm_set1_ps(1.9875691500E-4);
   //__m128 cephes_p1 = _mm_set1_ps(1.3981999507E-3);
   //__m128 cephes_p2 = _mm_set1_ps(8.3334519073E-3);
   //__m128 cephes_p3 = _mm_set1_ps(4.1665795894E-2);
   //__m128 cephes_p4 = _mm_set1_ps(1.6666665459E-1);
   //__m128 cephes_p5 = _mm_set1_ps(5.0000001201E-1);
   //for (INT4 ii=0; ii<roundedvectorlength; ii++) {
      //__m128 y = *x;
      //y = _mm_max_ps(y, minexp);
      //y = _mm_min_ps(y, maxexp);

      //Compute x*ln(2)
      //__m128 log2etimesx = _mm_mul_ps(log2e, y);

      //Round
      //__m128 log2etimesxplushalf = _mm_add_ps(log2etimesx, onehalf);
      //__m128i log2etimesx_rounded_i = _mm_cvttps_epi32(log2etimesxplushalf);
      //__m128 tmpval = _mm_cvtepi32_ps(log2etimesx_rounded_i);
      //__m128 mask = _mm_cmpgt_ps(tmpval, log2etimesxplushalf);
      //mask = _mm_and_ps(mask, one);
      //__m128 log2etimesx_rounded = _mm_sub_ps(tmpval, mask);

      //multiply and subtract as in the cephes code
      //tmpval = _mm_mul_ps(log2etimesx_rounded, cephes_c1);
      //__m128 tmpval2 = _mm_mul_ps(log2etimesx_rounded, cephes_c2);
      //y = _mm_sub_ps(y, tmpval);
      //y = _mm_sub_ps(y, tmpval2);

      //tmpval2 = _mm_mul_ps(y, y);

      //Pade approximation with polynomials
      //Now the polynomial part 1
      //__m128 polevlresult = cephes_p0;
      //polevlresult = _mm_mul_ps(polevlresult, y);
      //polevlresult = _mm_add_ps(polevlresult, cephes_p1);
      //polevlresult = _mm_mul_ps(polevlresult, y);
      //polevlresult = _mm_add_ps(polevlresult, cephes_p2);
      //polevlresult = _mm_mul_ps(polevlresult, y);
      //polevlresult = _mm_add_ps(polevlresult, cephes_p3);
      //polevlresult = _mm_mul_ps(polevlresult, y);
      //polevlresult = _mm_add_ps(polevlresult, cephes_p4);
      //polevlresult = _mm_mul_ps(polevlresult, y);
      //polevlresult = _mm_add_ps(polevlresult, cephes_p5);
      //polevlresult = _mm_mul_ps(polevlresult, tmpval2);
      //polevlresult = _mm_add_ps(polevlresult, y);
      //polevlresult = _mm_add_ps(polevlresult, one);

      //Now construct 2**n
      //log2etimesx_rounded_i = _mm_cvttps_epi32(log2etimesx_rounded);
      //__m128i log2etimesx_rounded_i_with_offset = _mm_add_epi32(log2etimesx_rounded_i, expoffset);
      //log2etimesx_rounded_i_with_offset = _mm_slli_epi32(log2etimesx_rounded_i_with_offset, 23);
      //__m128 pow2n = _mm_castsi128_ps(log2etimesx_rounded_i_with_offset);

      //And multiply
      // *result = _mm_mul_ps(polevlresult, pow2n);

      //x++;
      //result++;
   //}

   //Copy output aligned memory to non-aligned memory if necessary
   if (!outputaligned) memcpy(output->data, alignedoutput, 4*roundedvectorlength*sizeof(REAL4));

   //Finish up the remaining part
   for (INT4 ii=4*roundedvectorlength; ii<(INT4)input->length; ii++) {
      output->data[ii] = expf(input->data[ii]);
   }

   //FILE *EXPVALS = fopen("./output/expvals.dat","w");
   //for (ii=0; ii<(INT4)output->length; ii++) {
   //   fprintf(EXPVALS, "%f\n", output->data[ii]);
   //}
   //fclose(EXPVALS);

   //Free memory if necessary
   if (!vecaligned) XLALFree(allocinput);
   if (!outputaligned) XLALFree(allocoutput);

   return XLAL_SUCCESS;
#else
   (void)output;
   (void)input;
   fprintf(stderr, "%s: Failed because SSE2 is not supported, possibly because -msse2 flag wasn't used for compiling.\n", __func__);
   XLAL_ERROR(XLAL_EFAILED);
#endif

} */
