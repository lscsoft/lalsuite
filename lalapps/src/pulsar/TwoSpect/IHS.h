/*
*  Copyright (C) 2010, 2011, 2015 Evan Goetz
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

ihsMaximaStruct *createihsMaxima(const UINT4 fbins, const UINT4 rows);
void destroyihsMaxima(ihsMaximaStruct *data);
INT4 runIHS(ihsMaximaStruct *output, const ffdataStruct *input, const ihsfarStruct *ihsfarinput, const UserInput_t *params, UINT4 rows, const REAL4VectorAligned *aveNoise, const REAL4VectorAligned *FbinMean);

ihsVals * createihsVals(void);
void destroyihsVals(ihsVals *ihsvals);
INT4 incHarmSum(ihsVals *output, const REAL4VectorAligned *input, const UINT4 ihsfactor);
INT4 incHarmSumVector(REAL4VectorAligned *output, const REAL4VectorAligned *input, const UINT4 ihsfactor);
INT4 incHarmSumVectorWeighted(REAL4VectorAligned *output, const REAL4VectorAligned *input, const REAL4VectorAligned *aveNoise, const UINT4 ihsfactor);

ihsfarStruct * createihsfarStruct(const UINT4 rows, const UserInput_t *params);
void destroyihsfarStruct(ihsfarStruct *ihsfarstruct);
INT4 genIhsFar(ihsfarStruct *output, const UserInput_t *params, UINT4 rows, const REAL4VectorAligned *aveNoise, const gsl_rng *rng);

INT4 sumIHSarrayFAR(ihsfarStruct *outputfar, REAL4VectorAlignedArray *ihsvectorarray, const UINT4 rows, const REAL4VectorAligned *FbinMean, const UserInput_t *params, const gsl_rng *rng);
INT4 sumIHSarray(ihsMaximaStruct *output, const ihsfarStruct *inputfar, REAL4VectorAlignedArray *ihsvectorarray, const UINT4 rows, const REAL4VectorAligned *FbinMean, const UserInput_t *params);

INT4 findIHScandidates(candidateVector **candlist, const ihsfarStruct *ihsfarstruct, const UserInput_t *params, const ffdataStruct *ffdata, const ihsMaximaStruct *ihsmaxima, const REAL4VectorAligned *fbinavgs, const REAL4VectorSequence *trackedlines);

REAL4 ihsFOM(const INT4Vector *locs, const INT4 fomnorm);

REAL8 ihs2h0(const REAL8 ihsval, const UserInput_t *params);

#endif



