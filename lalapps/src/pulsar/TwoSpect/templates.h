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

#include "TwoSpectTypes.h"
#include "TwoSpect.h"

struct gsl_probR_pars {
   templateStruct *templatestruct;
   REAL4Vector *ffplanenoise;
   REAL4Vector *fbinaveratios;
   REAL4 threshold;
   INT4 errcode;
};

farStruct * new_farStruct(void);
void free_farStruct(farStruct *farstruct);
void estimateFAR(farStruct *output, templateStruct *templatestruct, INT4 trials, REAL8 thresh, REAL4Vector *ffplanenoise, REAL4Vector *fbinaveratios);
void numericFAR(farStruct *output, templateStruct *templatestruct, REAL8 thresh, REAL4Vector *ffplanenoise, REAL4Vector *fbinaveratios, INT4 method);
REAL8 gsl_probR(REAL8 R, void *pars);
REAL8 gsl_dprobRdR(REAL8 R, void *pars);
void gsl_probRandDprobRdR(REAL8 R, void *pars, REAL8 *probabilityR, REAL8 *dprobRdR);

templateStruct * new_templateStruct(INT4 length);
void free_templateStruct(templateStruct *nameoftemplate);
void makeTemplateGaussians(templateStruct *output, candidate input, inputParamsStruct *params);
void makeTemplate(templateStruct *output, candidate intput, inputParamsStruct *params, INT4Vector *sftexist, REAL4FFTPlan *plan);
void insertionSort_template(templateStruct *output, REAL4 weight, INT4 pixelloc, INT4 firstfftfreq, INT4 secfftfreq);

REAL8 probR(templateStruct *templatestruct, REAL4Vector *ffplanenoise, REAL4Vector *fbinaveratios, REAL8 R, INT4 *errcode);
REAL8 sincxoverxsqminusone(REAL8 overage);



#endif

