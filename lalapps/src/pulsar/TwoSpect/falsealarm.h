/*
*  Copyright (C) 2010, 2011, 2014 Evan Goetz
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

#ifndef __FALSEALARM_H__
#define __FALSEALARM_H__

#include "TwoSpectTypes.h"

struct gsl_probR_pars {
   TwoSpectTemplate *template;
   REAL4Vector *ffplanenoise;
   REAL4Vector *fbinaveratios;
   REAL4 threshold;
   UserInput_t *inputParams;
   gsl_rng *rng;
   INT4 errcode;
};

farStruct * new_farStruct(void);
void free_farStruct(farStruct *farstruct);
INT4 estimateFAR(farStruct *output, TwoSpectTemplate *template, INT4 trials, REAL8 thresh, REAL4Vector *ffplanenoise, REAL4Vector *fbinaveratios);
INT4 numericFAR(farStruct *output, TwoSpectTemplate *template, REAL8 thresh, REAL4Vector *ffplanenoise, REAL4Vector *fbinaveratios, UserInput_t *inputParams, gsl_rng *rng, INT4 method);
REAL8 probR(TwoSpectTemplate *template, REAL4Vector *ffplanenoise, REAL4Vector *fbinaveratios, REAL8 R, UserInput_t *params, gsl_rng *rng, INT4 *errcode);
REAL8 gsl_probR(REAL8 R, void *pars);
REAL8 gsl_dprobRdR(REAL8 R, void *pars);
void gsl_probRandDprobRdR(REAL8 R, void *pars, REAL8 *probabilityR, REAL8 *dprobRdR);

#endif
