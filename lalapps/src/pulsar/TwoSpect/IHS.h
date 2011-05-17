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

ihsMaximaStruct *new_ihsMaxima(INT4 fbins, INT4 columns);
void free_ihsMaxima(ihsMaximaStruct *data);
//void runIHS(ihsMaximaStruct *output, ffdataStruct *input, inputParamsStruct *params, INT4 columns);
void runIHS(ihsMaximaStruct *output, ffdataStruct *input, inputParamsStruct *params, INT4 columns, REAL4Vector *FbinMean);

ihsVals * new_ihsVals(void);
void free_ihsVals(ihsVals *ihsvals);
void incHarmSum(ihsVals *output, REAL4Vector *input);

ihsfarStruct * new_ihsfarStruct(INT4 columns);
void free_ihsfarStruct(ihsfarStruct *ihsfarstruct);
void genIhsFar(ihsfarStruct *output, inputParamsStruct *params, INT4 columns, REAL4Vector *aveNoise);

//void ihsSums(ihsMaximaStruct *output, REAL4Vector *ihss, INT4Vector *locs, INT4 cols);
void ihsSums(ihsMaximaStruct *output, REAL4Vector *ihss, INT4Vector *locs, INT4 cols, REAL4Vector *FbinMean);

void findIHScandidates(candidateVector *candlist, ihsfarStruct *ihsfarstruct, inputParamsStruct *params, ffdataStruct *ffdata, ihsMaximaStruct *ihsmaxima, REAL4Vector *fbinavgs);

REAL4 ihsFOM(REAL4Vector *ihss, INT4Vector *locs, REAL4Vector *sigma);
REAL4 ihsLoc(REAL4Vector *ihss, INT4Vector *locs, REAL4Vector *sigma);
//REAL4 ihsFOM(REAL4Vector *ihss, INT4Vector *locs, REAL4 sigma);
//REAL4 ihsLoc(REAL4Vector *ihss, INT4Vector *locs, REAL4 sigma);

#endif



