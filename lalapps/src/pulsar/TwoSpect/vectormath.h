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

#ifndef __VECTORMATH_H__
#define __VECTORMATH_H__

#include <lal/AVFactories.h>

INT4 fastSSVectorSequenceSum(REAL4VectorSequence *output, REAL4VectorSequence *input1, REAL4VectorSequence *input2, INT4 vectorpos1, INT4 vectorpos2, INT4 outputvectorpos);
INT4 sseSSVectorSequenceSum(REAL4VectorSequence *output, REAL4VectorSequence *input1, REAL4VectorSequence *input2, INT4 vectorpos1, INT4 vectorpos2, INT4 outputvectorpos, INT4 numvectors);
INT4 avxSSVectorSequenceSum(REAL4VectorSequence *output, REAL4VectorSequence *input1, REAL4VectorSequence *input2, INT4 vectorpos1, INT4 vectorpos2, INT4 outputvectorpos, INT4 numvectors);
INT4 fastSSVectorSequenceSubtract(REAL4Vector *output, REAL4VectorSequence *input1, REAL4Vector *input2, INT4 vectorpos1);
INT4 sseSSVectorSequenceSubtract(REAL4Vector *output, REAL4VectorSequence *input1, REAL4Vector *input2, INT4 vectorpos1);
INT4 avxSSVectorSequenceSubtract(REAL4Vector *output, REAL4VectorSequence *input1, REAL4Vector *input2, INT4 vectorpos1);
INT4 sse_exp_REAL8Vector(REAL8Vector *output, REAL8Vector *input);
INT4 sse_exp_REAL4Vector(REAL4Vector *output, REAL4Vector *input);

INT4 fastSSVectorMultiply_with_stride_and_offset(REAL4Vector *output, REAL4Vector *input1, REAL4Vector *input2, INT4 stride1, INT4 stride2, INT4 offset1, INT4 offset2);
INT4 sseSSVectorSum(REAL4Vector *output, REAL4Vector *input1, REAL4Vector *input2);
INT4 avxSSVectorSum(REAL4Vector *output, REAL4Vector *input1, REAL4Vector *input2);
INT4 sseSSVectorMultiply(REAL4Vector *output, REAL4Vector *input1, REAL4Vector *input2);
INT4 avxSSVectorMultiply(REAL4Vector *output, REAL4Vector *input1, REAL4Vector *input2);
INT4 sseAddScalarToREAL4Vector(REAL4Vector *output, REAL4Vector *input, REAL4 scalar);
INT4 avxAddScalarToREAL4Vector(REAL4Vector *output, REAL4Vector *input, REAL4 scalar);
INT4 sseScaleREAL4Vector(REAL4Vector *output, REAL4Vector *input, REAL4 scale);
INT4 avxScaleREAL4Vector(REAL4Vector *output, REAL4Vector *input, REAL4 scale);

INT4 sseAddScalarToREAL8Vector(REAL8Vector *output, REAL8Vector *input, REAL8 scalar);
INT4 avxAddScalarToREAL8Vector(REAL8Vector *output, REAL8Vector *input, REAL8 scalar);
INT4 sseScaleREAL8Vector(REAL8Vector *output, REAL8Vector *input, REAL8 scale);
INT4 avxScaleREAL8Vector(REAL8Vector *output, REAL8Vector *input, REAL8 scale);

INT4 sse_sin_cos_2PI_LUT_REAL8Vector(REAL8Vector *sin2pix_vector, REAL8Vector *cos2pix_vector, REAL8Vector *x);
INT4 sse_sin_cos_2PI_LUT_REAL4Vector(REAL4Vector *sin2pix_vector, REAL4Vector *cos2pix_vector, REAL4Vector *x);


#endif
