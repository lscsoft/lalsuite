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

#ifndef __UPPERLIMITS_H__
#define __UPPERLIMITS_H__

#include "TwoSpectTypes.h"

struct ncx2cdf_solver_params
{
   REAL8 val;
   REAL8 dof;
   REAL8 ULpercent;
};

UpperLimitVector * new_UpperLimitVector(UINT4 length);
UpperLimitVector * resize_UpperLimitVector(UpperLimitVector *vector, UINT4 length);
void reset_UpperLimitStruct(UpperLimit *ul);
void free_UpperLimitVector(UpperLimitVector *vector);
void free_UpperLimitStruct(UpperLimit *ul);


void skypoint95UL(UpperLimit *ul, inputParamsStruct *params, ffdataStruct *ffdata, ihsMaximaStruct *ihsmaxima, ihsfarStruct *ihsfar, REAL4Vector *fbinavgs);
REAL8 gsl_ncx2cdf_solver(REAL8 x, void *p);
REAL8 gsl_ncx2cdf_float_solver(REAL8 x, void *p);
REAL8 gsl_ncx2cdf_withouttinyprob_solver(REAL8 x, void *p);
REAL8 gsl_ncx2cdf_float_withouttinyprob_solver(REAL8 x, void *p);
REAL8 ncx2cdf_withouttinyprob_withmatlabchi2cdf_solver(REAL8 x, void *p);
REAL8 ncx2cdf_float_withouttinyprob_withmatlabchi2cdf_solver(REAL8 x, void *p);

void outputUpperLimitToFile(FILE *outputfile, UpperLimit ul, INT4 printAllULvalues);

#endif
