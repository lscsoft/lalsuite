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
#include "templates.h"


inputParamsStruct * new_inputParams(void);
void free_inputParams(inputParamsStruct *input);

ffdataStruct * new_ffdata(inputParamsStruct *param, INT4 mode);
void free_ffdata(ffdataStruct *data);
//void makeFakeBinarySigFF(ffdataStruct *ffdata, inputParamsStruct *input);


REAL4Vector * gaussRandNumVector(REAL4 sigma, UINT4 length, gsl_rng *ptrToGenerator);
REAL4Vector * expRandNumVector(REAL4 mu, UINT4 length, gsl_rng *ptrToGenerator);
REAL4Vector * readInSFTs(inputParamsStruct *input);
REAL4Vector * slideTFdata(inputParamsStruct *input, REAL4Vector *tfdata, INT4Vector *binshifts);
REAL4Vector * slideBackgroundData(inputParamsStruct *input, REAL4Vector *background, INT4Vector *binshifts);
REAL4Vector * ffPlaneNoise(inputParamsStruct *param, REAL4Vector *rngMeans, REAL4Vector *antPatternWeights);
REAL4Vector * tfWeightMeanSubtract(REAL4Vector *tfdata, REAL4Vector *rngMeans, REAL4Vector *antPatternWeights, inputParamsStruct *params);
REAL4Vector * tfRngMeans(REAL4Vector *tfdata, INT4 numffts, INT4 numfbins, INT4 blksize);
//REAL4Vector * makeFakeBinarySigTF(signalParamsStruct *in);
REAL4Vector * makeSecondFFT(REAL4Vector *tfdata, inputParamsStruct *params);

REAL4 calculateR(REAL4Vector *ffdata, templateStruct *templatestruct, REAL4Vector *noise);
REAL4 gaussRandNum(REAL4 sigma, gsl_rng *ptrToGenerator);
REAL4 expRandNum(REAL4 mu, gsl_rng *ptrToGenerator);
REAL4 maxModDepth(REAL4 period, REAL4 cohtime);
REAL4 minPeriod(REAL4 moddepth, REAL4 cohtime);
REAL4 calcMean(REAL4Vector *vector);
REAL4 calcStddev(REAL4Vector *vector);
REAL4 calcRms(REAL4Vector *vector);



#endif



