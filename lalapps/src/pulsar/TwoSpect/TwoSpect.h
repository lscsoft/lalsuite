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

ffdataStruct * new_ffdata(UserInput_t *params);
void free_ffdata(ffdataStruct *data);

INT4Vector * detectLines_simple(REAL4Vector *TFdata, ffdataStruct *ffdata, UserInput_t *params);
REAL4VectorSequence * trackLines(INT4Vector *lines, INT4Vector *binshifts, REAL4 minfbin, REAL4 df);
INT4 cleanLines(REAL4Vector *TFdata, REAL4Vector *background, INT4Vector *lines, UserInput_t *params, gsl_rng *rng);
INT4 makeSecondFFT(ffdataStruct *ffdata, REAL4Vector *tfdata, REAL4FFTPlan *plan);
INT4 ffPlaneNoise(REAL4VectorAligned *aveNoise, UserInput_t *params, INT4Vector *sftexist, REAL4Vector *backgrnd, REAL4Vector *antweights, REAL4FFTPlan *plan, REAL8 *normalization, gsl_rng *rng);

REAL4 avgTFdataBand(REAL4Vector *backgrnd, INT4 numfbins, INT4 numffts, INT4 binmin, INT4 binmax);
REAL4 rmsTFdataBand(REAL4Vector *backgrnd, INT4 numfbins, INT4 numffts, INT4 binmin, INT4 binmax);

INT4 readTwoSpectInputParams(UserInput_t *uvar, int argc, char *argv[]);
INT4 printREAL4Vector2File(REAL4Vector *vector, CHAR *filename);

#endif



