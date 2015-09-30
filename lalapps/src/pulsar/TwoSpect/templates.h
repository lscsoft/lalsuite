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

#ifndef __TEMPLATES_H__
#define __TEMPLATES_H__

#include <lal/RealFFT.h>
#include "TwoSpectTypes.h"

TwoSpectTemplate * createTwoSpectTemplate(const UINT4 length);
void resetTwoSpectTemplate(TwoSpectTemplate *template);
void destroyTwoSpectTemplate(TwoSpectTemplate *template);
TwoSpectTemplateVector * createTwoSpectTemplateVector(const UINT4 numTemplates, const UINT4 templateLength);
void destroyTwoSpectTemplateVector(TwoSpectTemplateVector *vector);
TwoSpectTemplateVector * generateTwoSpectTemplateVector(const REAL8 Pmin, const REAL8 Pmax, const REAL8 dfmin, const REAL8 dfmax, const REAL8 Tsft, const REAL8 SFToverlap, const REAL8 Tobs, const UINT4 maxvectorlength, const UINT4 minTemplateLength, const UINT4 maxTemplateLength, const UINT4 vectormathflag, const BOOLEAN exactflag);
INT4 writeTwoSpectTemplateVector(const TwoSpectTemplateVector *vector, const CHAR *filename);
TwoSpectTemplateVector * readTwoSpectTemplateVector(const CHAR *filename);
INT4 convertTemplateForSpecificFbin(TwoSpectTemplate *output, const TwoSpectTemplate *input, const REAL8 freq, const UserInput_t *params);

INT4 makeTemplateGaussians(TwoSpectTemplate *output, const candidate input, const UserInput_t *params);
INT4 makeTemplateGaussians2(TwoSpectTemplate *output, const REAL8 offset, const REAL8 P, const REAL8 deltaf, const REAL8 Tsft, const REAL8 SFToverlap, const REAL8 Tobs, const UINT4 minTemplateLength, const UINT4 vectormathflag);
INT4 makeTemplate(TwoSpectTemplate *output, const candidate intput, const UserInput_t *params, const INT4Vector *sftexist, const REAL4FFTPlan *plan);
INT4 makeTemplate2(TwoSpectTemplate *output, const REAL8 offset, const REAL8 P, const REAL8 deltaf, const REAL8 Tsft, const REAL8 SFToverlap, const REAL8 Tobs, const UINT4 minTemplateLength, const UINT4 vectormathflag, const REAL4FFTPlan *plan);
void insertionSort_template(TwoSpectTemplate *output, const REAL4 weight, const INT4 pixelloc);

REAL8 sincxoverxsqminusone(const REAL8 x);
REAL8 sqsincxoverxsqminusone(const REAL8 x);

#endif

