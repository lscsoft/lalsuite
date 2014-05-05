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

#ifndef __CANDIDATES_H__
#define __CANDIDATES_H__

#include "TwoSpectTypes.h"

candidateVector * new_candidateVector(UINT4 length);
candidateVector * resize_candidateVector(candidateVector *vector, UINT4 length);
candidateVector * keepMostSignificantCandidates(candidateVector *input, inputParamsStruct *params);
void free_candidateVector(candidateVector *vector);

void loadCandidateData(candidate* output, REAL8 fsig, REAL8 period, REAL8 moddepth, REAL4 ra, REAL4 dec, REAL8 statval, REAL8 h0, REAL8 prob, INT4 proberrcode, REAL8 normalization);

INT4 clusterCandidates(candidateVector **output, candidateVector *input, ffdataStruct *ffdata, inputParamsStruct *params, REAL4Vector *ffplanenoise, REAL4Vector *fbinaveratios, INT4Vector *sftexist, INT4 option);
INT4 testIHScandidates(candidateVector **output, candidateVector *ihsCandidates, ffdataStruct *ffdata, REAL4Vector *aveNoise, REAL4Vector *aveTFnoisePerFbinRatio, REAL4 alpha, REAL4 delta, inputParamsStruct *params);

REAL8 maxModDepth(REAL8 period, REAL8 cohtime);
REAL8 minPeriod(REAL8 moddepth, REAL8 cohtime);
REAL8 calculateR(REAL4Vector *ffdata, templateStruct *templatestruct, REAL4Vector *noise, REAL4Vector *fbinaveratios);

#endif


