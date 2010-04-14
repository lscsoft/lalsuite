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

#include <lal/LALStdlib.h>
#include <lal/DetResponse.h>
#include <lal/AVFactories.h>
#include <gsl/gsl_rng.h>
#include "TwoSpectTypes.h"


inputParamsStruct * new_inputParams(void);
void free_inputParams(inputParamsStruct *input);

ffdataStruct * new_ffdata(inputParamsStruct *param, INT4 mode);
void free_ffdata(ffdataStruct *data);

REAL8Vector * readInSFTs(inputParamsStruct *input);
REAL8Vector * slideTFdata(inputParamsStruct *input, REAL8Vector *tfdata, INT4Vector *binshifts);
REAL8Vector * ffPlaneNoise(inputParamsStruct *param, REAL8Vector *rngMeans, REAL8Vector *antPatternWeights);
REAL8Vector * tfWeightMeanSubtract(REAL8Vector *tfdata, REAL8Vector *rngMeans, REAL8Vector *antPatternWeights, inputParamsStruct *params);
REAL8Vector * tfRngMeans(REAL8Vector *tfdata, INT4 numffts, INT4 numfbins, INT4 blksize);
REAL8Vector * makeSecondFFT(REAL8Vector *tfdata, inputParamsStruct *params);

REAL8 calculateR(REAL8Vector *ffdata, templateStruct *templatestruct, REAL8Vector *noise);
REAL8 expRandNum(REAL8 mu, gsl_rng *ptrToGenerator);
REAL8 maxModDepth(REAL8 period, REAL8 cohtime);
REAL8 minPeriod(REAL8 moddepth, REAL8 cohtime);
REAL8 calcMean(REAL8Vector *vector);
REAL8 calcStddev(REAL8Vector *vector);
REAL8 calcRms(REAL8Vector *vector);
REAL8 minValue(REAL8Vector *vector);
REAL8 maxValue(REAL8Vector *vector);


#endif



