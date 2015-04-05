/*
*  Copyright (C) 2010, 2011 Evan Goetz
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

#ifndef __TEMPLATES_H__
#define __TEMPLATES_H__

#include <lal/RealFFT.h>
#include "TwoSpectTypes.h"

TwoSpectTemplate * new_TwoSpectTemplate(INT4 length);
void resetTwoSpectTemplate(TwoSpectTemplate *template);
void free_TwoSpectTemplate(TwoSpectTemplate *template);
TwoSpectTemplateVector * new_TwoSpectTemplateVector(UINT4 numTemplates, UINT4 templateLength);
void free_TwoSpectTemplateVector(TwoSpectTemplateVector *vector);
TwoSpectTemplateVector * generateTwoSpectTemplateVector(REAL8 Pmin, REAL8 Pmax, REAL8 dfmin, REAL8 dfmax, REAL8 Tsft, REAL8 SFToverlap, REAL8 Tobs, UINT4 maxvectorlength, UINT4 minTemplateLength, UINT4 maxTemplateLength, UINT4 vectormathflag, BOOLEAN exactflag);
INT4 writeTwoSpectTemplateVector(TwoSpectTemplateVector *vector, CHAR *filename);
TwoSpectTemplateVector * readTwoSpectTemplateVector(CHAR *filename);
INT4 convertTemplateForSpecificFbin(TwoSpectTemplate *output, TwoSpectTemplate *input, REAL8 freq, UserInput_t *params);

INT4 makeTemplateGaussians(TwoSpectTemplate *output, candidate input, UserInput_t *params);
INT4 makeTemplateGaussians2(TwoSpectTemplate *output, REAL8 offset, REAL8 P, REAL8 deltaf, REAL8 Tsft, REAL8 SFToverlap, REAL8 Tobs, UINT4 minTemplateLength, UINT4 vectormathflag);
INT4 makeTemplate(TwoSpectTemplate *output, candidate intput, UserInput_t *params, INT4Vector *sftexist, REAL4FFTPlan *plan);
INT4 makeTemplate2(TwoSpectTemplate *output, REAL8 offset, REAL8 P, REAL8 deltaf, REAL8 Tsft, REAL8 SFToverlap, REAL8 Tobs, UINT4 minTemplateLength, UINT4 vectormathflag, REAL4FFTPlan *plan);
void insertionSort_template(TwoSpectTemplate *output, REAL4 weight, INT4 pixelloc);

REAL8 sincxoverxsqminusone(REAL8 overage);
REAL8 sqsincxoverxsqminusone(REAL8 x);

INT4 truncated_sseScaleREAL4Vector(REAL4VectorAligned *output, REAL4VectorAligned *input, REAL4 scale, UINT4 length);
INT4 truncated_sseAddScalarToREAL4Vector(REAL4VectorAligned *output, REAL4VectorAligned *input, REAL4 scalar, UINT4 length);
INT4 truncated_sseSSVectorMultiply(REAL4VectorAligned *output, REAL4VectorAligned *input1, REAL4VectorAligned *input2, UINT4 length);
INT4 truncated_avxScaleREAL4Vector(REAL4VectorAligned *output, REAL4VectorAligned *input, REAL4 scale, UINT4 length);
INT4 truncated_avxAddScalarToREAL4Vector(REAL4VectorAligned *output, REAL4VectorAligned *input, REAL4 scalar, UINT4 length);
INT4 truncated_avxSSVectorMultiply(REAL4VectorAligned *output, REAL4VectorAligned *input1, REAL4VectorAligned *input2, UINT4 length);

#endif

