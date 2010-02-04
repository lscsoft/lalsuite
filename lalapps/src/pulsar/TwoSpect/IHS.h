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

#ifndef __IHS_H__
#define __IHS_H__

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include "TwoSpect.h"

typedef struct
{
   REAL4Vector *maxima;
   INT4Vector *locations;
   INT4 columns;
} ihsMaximaStruct;

typedef struct
{
   REAL4 ihs;
   INT4 loc;
} ihsVals;

typedef struct
{
   REAL4Vector *ihsfar;
   REAL4Vector *ihsdistMean;
   REAL4Vector *ihsdistSigma;
} ihsfarStruct;


ihsMaximaStruct *new_ihsMaxima(ffdataStruct *ffdata, INT4 columns);
void free_ihsMaxima(ihsMaximaStruct *data);
void runIHS(ihsMaximaStruct *out, ffdataStruct *in, INT4 columns);

ihsVals * new_ihsVals(void);
void free_ihsVals(ihsVals *ihsvals);
void incHarmSum(ihsVals *out, REAL4Vector *in);

ihsfarStruct * new_ihsfarStruct(INT4 columns);
void free_ihsfarStruct(ihsfarStruct *ihsfarstruct);
void genIhsFar(ihsfarStruct *out, ffdataStruct *ffdata, INT4 columns, REAL4 threshold);

REAL4Vector * ihsSums(REAL4Vector *ihss, INT4 cols);

REAL4 ihsFOM(REAL4Vector *ihss, INT4Vector *locs, REAL4Vector *expect);
REAL4 ihsLoc(REAL4Vector *ihss, INT4Vector *locs, REAL4Vector *expect);


#endif



