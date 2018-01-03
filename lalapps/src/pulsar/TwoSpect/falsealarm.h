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
   const TwoSpectTemplate *template;
   const REAL4VectorAligned *ffplanenoise;
   const REAL4VectorAligned *fbinaveratios;
   const REAL4 threshold;
   const UserInput_t *inputParams;
   const gsl_rng *rng;
   INT4 errcode;
};

farStruct * createfarStruct(void);
void destroyfarStruct(farStruct *farstruct);
INT4 estimateFAR(farStruct *output, const TwoSpectTemplate *template, const UINT4 trials, const REAL8 thresh, const REAL4VectorAligned *ffplanenoise, const REAL4VectorAligned *fbinaveratios);
INT4 numericFAR(farStruct *output, const TwoSpectTemplate *template, const REAL8 thresh, const REAL4VectorAligned *ffplanenoise, const REAL4VectorAligned *fbinaveratios, const UserInput_t *inputParams, const gsl_rng *rng, const INT4 method);
REAL8 probR(const TwoSpectTemplate *template, const REAL4VectorAligned *ffplanenoise, const REAL4VectorAligned *fbinaveratios, const REAL8 R, const UserInput_t *params, const gsl_rng *rng, INT4 *errcode);
REAL8 gsl_probR(const REAL8 R, void *pars);
REAL8 gsl_dprobRdR(const REAL8 R, void *pars);
void gsl_probRandDprobRdR(const REAL8 R, void *pars, REAL8 *probabilityR, REAL8 *dprobRdR);

#endif
