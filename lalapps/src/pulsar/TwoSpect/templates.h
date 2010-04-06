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

#ifndef __TEMPLATES_H__
#define __TEMPLATES_H__

#include <lal/RealFFT.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include "candidates.h"
#include "TwoSpect.h"
#include "TwoSpectTypes.h"


farStruct * new_farStruct(void);
void free_farStruct(farStruct *farstruct);
//void estimateFAR(farStruct *out, REAL4Vector *weights, topbinsStruct *topbinsstruct, REAL4 thresh, REAL4Vector *ffplanenoise);
void estimateFAR(farStruct *out, templateStruct *templatestruct, INT4 trials, REAL8 thresh, REAL8Vector *ffplanenoise);

templateStruct * new_templateStruct(INT4 length);
void free_templateStruct(templateStruct *nameoftemplate);
//void makeTemplateGaussians(ffdataStruct *out, candidate *in);
//void makeTemplate(ffdataStruct *out, candidate *in, REAL4FFTPlan *plan);
void makeTemplateGaussians(templateStruct *out, candidate *in, inputParamsStruct *params);
void makeTemplate(templateStruct *out, candidate *in, inputParamsStruct *params, REAL8FFTPlan *plan);

REAL8 sincxoverxsqminusone(REAL8 overage);

#endif

