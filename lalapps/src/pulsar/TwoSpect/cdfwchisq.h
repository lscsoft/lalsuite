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

#ifndef __CDFWCHISQ_H__
#define __CDFWCHISQ_H__

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>

typedef struct
{
   REAL8 sigsq;      //The value in front of the additional normally distributed random variable Q = Sum(w_i X_i) + sigma X_0
   REAL8 wnmax;      //Maximum value of weights*noise
   REAL8 wnmin;      //Minimum values of weights*noise
   REAL8 wnmean;     //The 'mean' value (not really the true mean)
   REAL8 c;          //The threshold value for Prob(Q < c)
   REAL8 integrationValue;       //Integration value
   REAL8 integrationError;       //Error of integration
   INT4 count;       //Count for number of times entering specific functions
   INT4 lim;         //Limit to number of integration terms
   INT4 arrayNotSorted;      //Array has not been sorted
   INT4 fail;        //Fail flag if integration failes
   INT4 useSSE;      //Flag to specify use SSE integration function
   INT4 useAVX;      //Flag to specify use AVX integration function
   INT4Vector *dofs;     //Array to hold values of the d.o.f. for each chi-squared variable
   INT4Vector *sorting;      //Array to hold the sorted element values for weights*noise
   REAL8Vector *weights;         //Array of weights in front of each chi-squared variable to sum (in my case, weight*noise/2.0)
   REAL8Vector *noncentrality;   //Array of non-centrality parameters for each chi-squared variable in sum
} qfvars;

REAL8 cdfwchisq(qfvars *vars, REAL8 sigma, REAL8 acc, INT4 *ifault);
REAL8 cdfwchisq_twospect(qfvars *vars, REAL8 sigma, REAL8 acc, INT4 *ifault);

void order(qfvars *vars);
void findu(qfvars *vars, REAL8* utx, REAL8 accx);
void findu_twospect(qfvars *vars, REAL8* utx, REAL8 accx);
void integrate(qfvars *vars, INT4 nterm, REAL8 interv, REAL8 tausq, INT4 mainx);
void integrate_eg(qfvars *vars, INT4 nterm, REAL8 interv, REAL8 tausq, INT4 mainx);
void integrate_twospect(qfvars *vars, INT4 nterm, REAL8 interv, REAL8 tausq, INT4 mainx);
void integrate_twospect2(qfvars *vars, INT4 nterm, REAL8 interv, REAL8 tausq, INT4 mainx);
INT4 fast_integrate_twospect2(qfvars *vars, INT4 nterm, REAL8 interv, REAL8 tausq, INT4 mainx);

REAL8 exp1(REAL8 x);
REAL8 twospect_log_1plusx(REAL8 x);
REAL8 twospect_log_1plusx_mx(REAL8 x);
REAL8 errbound(qfvars *vars, REAL8 u, REAL8* cx);
REAL8 errbound_twospect(qfvars *vars, REAL8 u, REAL8* cx);
REAL8 cutoff(qfvars *vars, REAL8 accx, REAL8* upn);
REAL8 cutoff_twospect(qfvars *vars, REAL8 accx, REAL8* upn);
REAL8 truncation(qfvars *vars, REAL8 u, REAL8 tausq);
REAL8 truncation_twospect(qfvars *vars, REAL8 u, REAL8 tausq);
REAL8 coeff(qfvars *vars, REAL8 x);
REAL8 coeff_twospect(qfvars *vars, REAL8 x);

int compar(void *p, const void *a, const void *b);

int compar(void *p, const void *a, const void *b);

#endif
