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

#ifndef __IHS_H__
#define __IHS_H__

#include "TwoSpectTypes.h"

ihsMaximaStruct *new_ihsMaxima(INT4 fbins, INT4 rows);
void free_ihsMaxima(ihsMaximaStruct *data);
INT4 runIHS(ihsMaximaStruct *output, ffdataStruct *input, ihsfarStruct *ihsfarinput, inputParamsStruct *params, INT4 rows, REAL4Vector *aveNoise, REAL4Vector *FbinMean);

ihsVals * new_ihsVals(void);
void free_ihsVals(ihsVals *ihsvals);
INT4 incHarmSum(ihsVals *output, REAL4Vector *input, INT4 ihsfactor);
INT4 incHarmSumVector(REAL4Vector *output, REAL4Vector *input, INT4 ihsfactor);
INT4 incHarmSumVectorWeighted(REAL4Vector *output, REAL4Vector *input, REAL4Vector *aveNoise, INT4 ihsfactor);

ihsfarStruct * new_ihsfarStruct(INT4 rows, inputParamsStruct *params);
void free_ihsfarStruct(ihsfarStruct *ihsfarstruct);
INT4 genIhsFar(ihsfarStruct *output, inputParamsStruct *params, INT4 rows, REAL4Vector *aveNoise);

INT4 sumIHSSequenceFAR(ihsfarStruct *outputfar, REAL4VectorSequence *ihsvectorsequence, INT4 rows, REAL4Vector *FbinMean, inputParamsStruct *params);
INT4 sumIHSSequence(ihsMaximaStruct *output, ihsfarStruct *inputfar, REAL4VectorSequence *ihsvectorsequence, INT4 rows, REAL4Vector *FbinMean, inputParamsStruct *params);

void SSVectorSequenceSum(REAL4VectorSequence *output, REAL4VectorSequence *input1, REAL4VectorSequence *input2, INT4 vectorpos1, INT4 vectorpos2, INT4 outputvectorpos);

INT4 findIHScandidates(candidateVector **candlist, ihsfarStruct *ihsfarstruct, inputParamsStruct *params, ffdataStruct *ffdata, ihsMaximaStruct *ihsmaxima, REAL4Vector *fbinavgs, REAL4VectorSequence *trackedlines);

REAL4 ihsFOM(INT4Vector *locs, INT4 fomnorm);

REAL8 ihs2h0(REAL8 ihsval, inputParamsStruct *params);

#endif



