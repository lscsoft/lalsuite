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

#include "TwoSpect.h"

struct gsl_probR_pars {
   templateStruct *templatestruct;
   REAL4Vector *ffplanenoise;
   REAL4Vector *fbinaveratios;
   REAL4 threshold;
   inputParamsStruct *inputParams;
   INT4 errcode;
};

farStruct * new_farStruct(void);
void free_farStruct(farStruct *farstruct);
void estimateFAR(farStruct *output, templateStruct *templatestruct, INT4 trials, REAL8 thresh, REAL4Vector *ffplanenoise, REAL4Vector *fbinaveratios);
INT4 numericFAR(farStruct *output, templateStruct *templatestruct, REAL8 thresh, REAL4Vector *ffplanenoise, REAL4Vector *fbinaveratios, inputParamsStruct *inputParams, INT4 method);
REAL8 gsl_probR(REAL8 R, void *pars);
REAL8 gsl_dprobRdR(REAL8 R, void *pars);
void gsl_probRandDprobRdR(REAL8 R, void *pars, REAL8 *probabilityR, REAL8 *dprobRdR);

templateStruct * new_templateStruct(INT4 length);
void resetTemplateStruct(templateStruct *templatestruct);
void free_templateStruct(templateStruct *nameoftemplate);

INT4 makeTemplateGaussians(templateStruct *output, candidate input, inputParamsStruct *params, INT4 numfbins, INT4 numfprbins);
INT4 makeTemplate(templateStruct *output, candidate intput, inputParamsStruct *params, INT4Vector *sftexist, REAL4FFTPlan *plan);
void insertionSort_template(templateStruct *output, REAL4 weight, INT4 pixelloc, INT4 firstfftfreq, INT4 secfftfreq);
INT4 analyzeOneTemplate(candidate *output, candidate *input, ffdataStruct *ffdata, REAL4Vector *aveNoise, REAL4Vector *aveTFnoisePerFbinRatio, inputParamsStruct *params, INT4Vector *sftexist, REAL4FFTPlan *plan);
INT4 bruteForceTemplateSearch(candidate *output, candidate input, REAL8 fminimum, REAL8 fmaximum, INT4 numfsteps, INT4 numperiodslonger, INT4 numperiodsshorter, REAL4 periodSpacingFactor, REAL8 dfmin, REAL8 dfmax, INT4 numdfsteps, inputParamsStruct *params, REAL4Vector *ffdata, INT4Vector *sftexist, REAL4Vector *aveNoise, REAL4Vector *aveTFnoisePerFbinRatio, REAL4FFTPlan *secondFFTplan, INT4 useExactTemplates);
INT4 bruteForceTemplateTest(candidateVector **output, candidate input, REAL8 fminimum, REAL8 fmaximum, INT4 numfsteps, INT4 numperiodslonger, INT4 numperiodsshorter, REAL4 periodSpacingFactor, REAL8 dfmin, REAL8 dfmax, INT4 numdfsteps, inputParamsStruct *params, REAL4Vector *ffdata, INT4Vector *sftexist, REAL4Vector *aveNoise, REAL4Vector *aveTFnoisePerFbinRatio, REAL4FFTPlan *secondFFTplan, INT4 useExactTemplates);
INT4 templateSearch_scox1Style(candidateVector **output, REAL8 fminimum, REAL8 fspan, REAL8 period, REAL8 asini, inputParamsStruct *params, REAL4Vector *ffdata, INT4Vector *sftexist, REAL4Vector *aveNoise, REAL4Vector *aveTFnoisePerFbinRatio, REAL4FFTPlan *secondFFTplan, INT4 useExactTemplates);

void efficientTemplateSearch(candidate *output, candidate input, REAL8 fminimum, REAL8 fmaximum, REAL8 minfstep, INT4 numperiods, REAL8 dfmin, REAL8 dfmax, REAL8 minDfstep, inputParamsStruct *params, REAL4Vector *ffdata, INT4Vector *sftexist, REAL4Vector *aveNoise, REAL4Vector *aveTFnoisePerFbinRatio, REAL4FFTPlan *secondFFTplan, INT4 useExactTemplates);

REAL8 probR(templateStruct *templatestruct, REAL4Vector *ffplanenoise, REAL4Vector *fbinaveratios, REAL8 R, inputParamsStruct *params, INT4 *errcode);
REAL8 sincxoverxsqminusone(REAL8 overage);
REAL8 sqsincxoverxsqminusone(REAL8 x);

INT4 twospect_sin_cos_2PI_LUT(REAL8 *sin2pix, REAL8 *cos2pix, REAL8 x);
INT4 twospect_sin_cos_LUT(REAL8 *sinx, REAL8 *cosx, REAL8 x);

#endif

