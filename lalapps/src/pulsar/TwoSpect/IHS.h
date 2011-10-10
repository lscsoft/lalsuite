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
void runIHS(ihsMaximaStruct *output, ffdataStruct *input, inputParamsStruct *params, INT4 rows, REAL4Vector *aveNoise, REAL4Vector *FbinMean);

ihsVals * new_ihsVals(void);
void free_ihsVals(ihsVals *ihsvals);
void incHarmSum(ihsVals *output, REAL4Vector *input, INT4 ihsfactor);
void incHarmSumVector(REAL4Vector *output, REAL4Vector *input, INT4 ihsfactor);

ihsfarStruct * new_ihsfarStruct(INT4 rows);
void free_ihsfarStruct(ihsfarStruct *ihsfarstruct);
void genIhsFar(ihsfarStruct *output, inputParamsStruct *params, INT4 rows, REAL4Vector *aveNoise);

void ihsSums(ihsMaximaStruct *output, REAL4Vector *ihss, INT4Vector *locs, INT4 rows, REAL4Vector *FbinMean, INT4 locationnormfactor);
void ihsSums2_withFAR(ihsMaximaStruct *output, ihsfarStruct *outputfar, REAL4VectorSequence *ihsvectorsequence, REAL4Vector *ihss, INT4Vector *locs, REAL4Vector *aveNoise, INT4 rows, REAL4Vector *FbinMean, INT4 locationnormfactor, inputParamsStruct *params, INT4 calcPInvVals);
void ihsSums2_withFAR_withnoise(ihsMaximaStruct *output, ihsfarStruct *outputfar, REAL4VectorSequence *ihsvectorsequence, REAL4Vector *ihss, INT4Vector *locs, REAL4Vector *aveNoise, INT4 rows, REAL4Vector *FbinMean, INT4 locationnormfactor, inputParamsStruct *params, INT4 calcPInvVals);
REAL4VectorSequence * ihsVectorSums(REAL4VectorSequence *input, INT4 rows);

void findIHScandidates(candidateVector *candlist, ihsfarStruct *ihsfarstruct, inputParamsStruct *params, ffdataStruct *ffdata, ihsMaximaStruct *ihsmaxima, REAL4Vector *aveNoise, REAL4Vector *fbinavgs);

REAL4 ihsFOM(REAL4Vector *ihss, INT4Vector *locs, REAL4Vector *sigma, INT4 locationnormfactor);
REAL4 ihsLoc(REAL4Vector *ihss, INT4Vector *locs, REAL4Vector *sigma);

REAL8 ihs2h0_withNoiseSubtraction(REAL8 ihsval, INT4 location, INT4 lowestfrequencybin, INT4 rows, inputParamsStruct *params, REAL4Vector *aveNoise, REAL4Vector *fbinavgs);
//REAL8 ihs2h0(REAL8 ihsval, inputParamsStruct *params);
REAL8 ihs2h0(REAL8 ihsval, inputParamsStruct *params);
REAL8 significance_of_IHSval(REAL8 ihsval, INT4 location, INT4 lowestfrequencybin, INT4 rows, inputParamsStruct *params, REAL4Vector *aveNoise, REAL4Vector *fbinavgs);

#endif



