/*
*  Copyright (C) 2010 Evan Goetz
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

#ifndef __TWOSPECT_H__
#define __TWOSPECT_H__

#include <lal/RealFFT.h>

#include <gsl/gsl_rng.h>
#include "TwoSpectTypes.h"


inputParamsStruct * new_inputParams(void);
void free_inputParams(inputParamsStruct *input);

ffdataStruct * new_ffdata(inputParamsStruct *input);
void free_ffdata(ffdataStruct *data);

REAL4Vector * readInSFTs(inputParamsStruct *input, REAL8 *normalization);

void slideTFdata(REAL4Vector *output, inputParamsStruct *input, REAL4Vector *tfdata, INT4Vector *binshifts);
void tfWeightMeanSubtract(REAL4Vector *output, REAL4Vector *tfdata, REAL4Vector *rngMeans, REAL4Vector *antPatternWeights, inputParamsStruct *input);
void tfRngMeans(REAL4Vector *output, REAL4Vector *tfdata, INT4 numffts, INT4 numfbins, INT4 blksize);
void makeSecondFFT(REAL4Vector *output, REAL8 normalization, REAL4Vector *tfdata, inputParamsStruct *input, REAL4FFTPlan *plan);
void ffPlaneNoise(REAL4Vector *aveNoise, inputParamsStruct *input, REAL4Vector *backgrnd, REAL4Vector *antweights, REAL8 *normalization);

REAL8 expRandNum(REAL8 mu, gsl_rng *ptrToGenerator);
REAL8 maxModDepth(REAL8 period, REAL8 cohtime);
REAL8 minPeriod(REAL8 moddepth, REAL8 cohtime);
REAL8 calculateR(REAL4Vector *ffdata, templateStruct *templatestruct, REAL4Vector *noise, REAL4Vector *fbinaveratios);
REAL8 calcMeanD(REAL8Vector *vector);
REAL8 calcStddevD(REAL8Vector *vector);

REAL4 avgTFdataBand(REAL4Vector *backgrnd, INT4 numfbins, INT4 numffts, INT4 binmin, INT4 binmax);
REAL4 rmsTFdataBand(REAL4Vector *backgrnd, INT4 numfbins, INT4 numffts, INT4 binmin, INT4 binmax);
REAL4 calcMean(REAL4Vector *vector);
REAL4 calcStddev(REAL4Vector *vector);
REAL4 calcRms(REAL4Vector *vector);

#endif



