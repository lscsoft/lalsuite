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

#ifndef __CANDIDATES_H__
#define __CANDIDATES_H__

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include "TwoSpect.h"

typedef struct
{
   REAL4 fsig; /* 0 value means candidate not valid */
   REAL4 period;
   REAL4 moddepth;
   REAL4 Tobs;
   REAL4 Tcoh;
   REAL4 fmin;
   REAL4 fspan;
   REAL4 stat;
   REAL4 snr;
} candidate;


candidate * new_candidate(void);
void free_candidate(candidate *cand);
void loadCandidateData(candidate *out, REAL4 fsig, REAL4 period, REAL4 moddepth, REAL4 Tobs, REAL4 Tcoh, REAL4 fmin, REAL4 fspan, REAL4 stat, REAL4 snr);
void clusterCandidates(candidate *out[], candidate *in[], ffdataStruct *ffdata, inputParamsStruct *params, REAL4Vector *ffplanenoise, INT4 numofcandidates, INT4 option);





#endif


