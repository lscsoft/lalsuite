/*
*  Copyright (C) 2010 -- 2014 Evan Goetz
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

#include "TwoSpectTypes.h"

ffdataStruct * createffdata(const UserInput_t *params);
void destroyffdata(ffdataStruct *data);

INT4Vector * detectLines_simple(const REAL4VectorAligned *TFdata, const ffdataStruct *ffdata, const UserInput_t *params);
REAL4VectorSequence * trackLines(const INT4Vector *lines, const INT4Vector *binshifts, const REAL4 minfbin, const REAL4 df);
INT4 cleanLines(REAL4VectorAligned *TFdata, const REAL4VectorAligned *background, const INT4Vector *lines, const UserInput_t *params, const gsl_rng *rng);
INT4 makeSecondFFT(ffdataStruct *ffdata, REAL4VectorAligned *tfdata, const REAL4FFTPlan *plan);
INT4 ffPlaneNoise(REAL4VectorAligned *aveNoise, const UserInput_t *params, const INT4Vector *sftexist, const REAL4VectorAligned *aveNoiseInTime, const REAL4VectorAligned *antweights, const REAL4VectorAligned *backgroundScaling, const REAL4FFTPlan *plan, const REAL4VectorAligned *expDistVals, const gsl_rng *rng, REAL8 *normalization);

REAL4 avgTFdataBand(const REAL4VectorAligned *backgrnd, UINT4 numfbins, UINT4 numffts, UINT4 binmin, UINT4 binmax);
REAL4 rmsTFdataBand(const REAL4VectorAligned *backgrnd, UINT4 numfbins, UINT4 numffts, UINT4 binmin, UINT4 binmax);
REAL4VectorAligned * calcAveTFnoisePerFbinRatio(const REAL4VectorAligned *background, const REAL4VectorAligned *backgroundScaling, const UINT4 numffts);
INT4 medianBackgroundBandInTime(REAL4VectorAligned *aveNoiseInTime, const REAL4VectorAligned *backgrnd, const INT4Vector *sftexist);

MultiLALDetector * setupMultiLALDetector(LALStringVector *IFO);

INT4 readTwoSpectInputParams(UserInput_t *uvar, int argc, char *argv[]);
INT4 printREAL4Vector2File(const REAL4Vector *vector, const CHAR *directory, const CHAR *filename);

#endif



