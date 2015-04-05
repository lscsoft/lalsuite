/*
 *  Copyright (C) 2011, 2015 Evan Goetz
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

#ifndef __VECTORMATH_H__
#define __VECTORMATH_H__

#include <lal/AVFactories.h>
#include <lal/VectorOps.h>
#include <lal/VectorMath.h>
#include "TwoSpectTypes.h"

alignedREAL8Vector * createAlignedREAL8Vector(UINT4 length, const size_t align);
void destroyAlignedREAL8Vector(alignedREAL8Vector *vector);
alignedREAL8VectorArray * createAlignedREAL8VectorArray(const UINT4 length, const UINT4 vectorLength, const size_t align);
void destroyAlignedREAL8VectorArray(alignedREAL8VectorArray *array);
alignedREAL4VectorArray * createAlignedREAL4VectorArray(const UINT4 length, const UINT4 vectorLength, const size_t align);
void destroyAlignedREAL4VectorArray(alignedREAL4VectorArray *array);

INT4 sseSSVectorArraySum(alignedREAL4VectorArray *output, alignedREAL4VectorArray *input1, alignedREAL4VectorArray *input2, INT4 vectorpos1, INT4 vectorpos2, INT4 outputvectorpos, INT4 numvectors);
INT4 avxSSVectorArraySum(alignedREAL4VectorArray *output, alignedREAL4VectorArray *input1, alignedREAL4VectorArray *input2, INT4 vectorpos1, INT4 vectorpos2, INT4 outputvectorpos, INT4 numvectors);

INT4 sse_exp_REAL8Vector(alignedREAL8Vector *output, alignedREAL8Vector *input);
//INT4 sse_exp_REAL4Vector(REAL4Vector *output, REAL4Vector *input);

INT4 fastSSVectorMultiply_with_stride_and_offset(REAL4Vector *output, REAL4Vector *input1, REAL4Vector *input2, INT4 stride1, INT4 stride2, INT4 offset1, INT4 offset2);
INT4 sseSSVectorSum(REAL4VectorAligned *output, REAL4VectorAligned *input1, REAL4VectorAligned *input2);
INT4 avxSSVectorSum(REAL4VectorAligned *output, REAL4VectorAligned *input1, REAL4VectorAligned *input2);
INT4 sseSSVectorSubtract(REAL4VectorAligned *output, REAL4VectorAligned *input1, REAL4VectorAligned *input2);
INT4 avxSSVectorSubtract(REAL4VectorAligned *output, REAL4VectorAligned *input1, REAL4VectorAligned *input2);
INT4 sseSSVectorMultiply(REAL4VectorAligned *output, REAL4VectorAligned *input1, REAL4VectorAligned *input2);
INT4 avxSSVectorMultiply(REAL4VectorAligned *output, REAL4VectorAligned *input1, REAL4VectorAligned *input2);
INT4 sseAddScalarToREAL4Vector(REAL4VectorAligned *output, REAL4VectorAligned *input, REAL4 scalar);
INT4 avxAddScalarToREAL4Vector(REAL4VectorAligned *output, REAL4VectorAligned *input, REAL4 scalar);
INT4 sseScaleREAL4Vector(REAL4VectorAligned *output, REAL4VectorAligned *input, REAL4 scale);
INT4 avxScaleREAL4Vector(REAL4VectorAligned *output, REAL4VectorAligned *input, REAL4 scale);

INT4 sseAddScalarToREAL8Vector(alignedREAL8Vector *output, alignedREAL8Vector *input, REAL8 scalar);
INT4 avxAddScalarToREAL8Vector(alignedREAL8Vector *output, alignedREAL8Vector *input, REAL8 scalar);
INT4 sseScaleREAL8Vector(alignedREAL8Vector *output, alignedREAL8Vector *input, REAL8 scale);
INT4 avxScaleREAL8Vector(alignedREAL8Vector *output, alignedREAL8Vector *input, REAL8 scale);
INT4 sseDDVectorSum(alignedREAL8Vector *output, alignedREAL8Vector *input1, alignedREAL8Vector *input2);
INT4 avxDDVectorSum(alignedREAL8Vector *output, alignedREAL8Vector *input1, alignedREAL8Vector *input2);
INT4 sseDDVectorSubtract(alignedREAL8Vector *output, alignedREAL8Vector *input1, alignedREAL8Vector *input2);
INT4 avxDDVectorSubtract(alignedREAL8Vector *output, alignedREAL8Vector *input1, alignedREAL8Vector *input2);
INT4 sseDDVectorMultiply(alignedREAL8Vector *output, alignedREAL8Vector *input1, alignedREAL8Vector *input2);
INT4 avxDDVectorMultiply(alignedREAL8Vector *output, alignedREAL8Vector *input1, alignedREAL8Vector *input2);
INT4 sseInvertREAL8Vector(alignedREAL8Vector *output, alignedREAL8Vector *input);
INT4 avxInvertREAL8Vector(alignedREAL8Vector *output, alignedREAL8Vector *input);

//INT4 sse_sin_cos_2PI_LUT_REAL8Vector(REAL8Vector *sin2pix_vector, REAL8Vector *cos2pix_vector, REAL8Vector *x);
//INT4 sse_sin_cos_2PI_LUT_REAL4Vector(REAL4Vector *sin2pix_vector, REAL4Vector *cos2pix_vector, REAL4Vector *x);


#endif
