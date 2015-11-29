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
REAL4VectorAlignedArray * createREAL4VectorAlignedArray(const UINT4 length, const UINT4 vectorLength, const size_t align);
void destroyREAL4VectorAlignedArray(REAL4VectorAlignedArray *array);

INT4 DirichletRatioVector(COMPLEX8Vector *output, alignedREAL8Vector *delta0, alignedREAL8Vector *delta1, alignedREAL8Vector *scaling, const UserInput_t *params);

INT4 VectorArraySum(REAL4VectorAlignedArray *output, REAL4VectorAlignedArray *input1, REAL4VectorAlignedArray *input2, INT4 vectorpos1, INT4 vectorpos2, INT4 outputvectorpos, INT4 numvectors);

INT4 sse_exp_REAL8Vector(alignedREAL8Vector *output, alignedREAL8Vector *input);
//INT4 sse_exp_REAL4Vector(REAL4Vector *output, REAL4Vector *input);

INT4 fastSSVectorMultiply_with_stride_and_offset(REAL4VectorAligned *output, const REAL4VectorAligned *input1, const REAL4VectorAligned *input2, const INT4 stride1, const INT4 stride2, const INT4 offset1, const INT4 offset2);

INT4 VectorSubtractREAL4(REAL4VectorAligned *output, REAL4VectorAligned *input1, REAL4VectorAligned *input2, INT4 vectorMath);
INT4 VectorScaleREAL8(alignedREAL8Vector *output, alignedREAL8Vector *input, REAL8 scale, INT4 vectorMath);
INT4 VectorShiftREAL8(alignedREAL8Vector *output, alignedREAL8Vector *input, REAL8 shift, INT4 vectorMath);
INT4 VectorAddREAL8(alignedREAL8Vector *output, alignedREAL8Vector *input1, alignedREAL8Vector *input2, INT4 vectorMath);
INT4 VectorSubtractREAL8(alignedREAL8Vector *output, alignedREAL8Vector *input1, alignedREAL8Vector *input2, INT4 vectorMath);
INT4 VectorMultiplyREAL8(alignedREAL8Vector *output, alignedREAL8Vector *input1, alignedREAL8Vector *input2, INT4 vectorMath);
INT4 VectorInvertREAL8(alignedREAL8Vector *output, alignedREAL8Vector *input, INT4 vectorMath);

INT4 VectorFloorREAL8(alignedREAL8Vector *output, alignedREAL8Vector *input, INT4 vectorMath);
INT4 VectorRoundREAL4(REAL4VectorAligned *output, REAL4VectorAligned *input, INT4 vectorMath);
INT4 VectorRoundREAL8(alignedREAL8Vector *output, alignedREAL8Vector *input, INT4 vectorMath);
INT4 VectorAbsREAL4(REAL4VectorAligned *output, REAL4VectorAligned *input, INT4 vectorMath);
INT4 VectorAbsREAL8(alignedREAL8Vector *output, alignedREAL8Vector *input, INT4 vectorMath);
INT4 VectorCabsfCOMPLEX8(REAL4VectorAligned *output, COMPLEX8Vector *input, INT4 vectorMath);
INT4 VectorCabsCOMPLEX8(alignedREAL8Vector *output, COMPLEX8Vector *input, INT4 vectorMath);

INT4 sseSSVectorSubtract(REAL4VectorAligned *output, REAL4VectorAligned *input1, REAL4VectorAligned *input2);
INT4 avxSSVectorSubtract(REAL4VectorAligned *output, REAL4VectorAligned *input1, REAL4VectorAligned *input2);

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
