/*
 *  Copyright (C) 2011 Evan Goetz
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


#ifndef __TWOSPECTSPECFUNC_H__
#define __TWOSPECTSPECFUNC_H__

#include <lal/LALStdlib.h>

struct cheb_series_struct {
   double * c;   /* coefficients                */
   int order;    /* order of expansion          */
   double a;     /* lower interval point        */
   double b;     /* upper interval point        */
   int order_sp; /* effective single precision order */
};
typedef struct cheb_series_struct cheb_series;

COMPLEX16 DirichletKernelLargeN(const REAL8 delta);
COMPLEX16 DirichletKernelLargeNHann(const REAL8 delta);
INT4 DirichletKernalLargeNHannRatio(COMPLEX8 *ratio, const REAL4 delta0, const REAL4 delta1, const REAL4 scaling);

REAL8 twospect_small(REAL8 q);
REAL8 twospect_intermediate(REAL8 r);
REAL8 twospect_tail(REAL8 r);
REAL8 rat_eval(const REAL8 a[], const size_t na, const REAL8 b[], const size_t nb, const REAL8 x);
REAL8 ran_gamma_pdf(REAL8 x, REAL8 a, REAL8 b);
INT4 sf_gamma_inc_P(REAL8 *out, REAL8 a, REAL8 x);
INT4 sf_gamma_inc_Q(REAL8 *out, REAL8 a, REAL8 x);
REAL8 matlab_gamma_inc(REAL8 x, REAL8 a, INT4 upper);
INT4 gamma_inc_P_series(REAL8 *out, REAL8 a, REAL8 x);
INT4 gamma_inc_Q_series(REAL8 *out, REAL8 a, REAL8 x);
INT4 gamma_inc_D(REAL8 *out, REAL8 a, REAL8 x);
INT4 twospect_sf_gammastar(REAL8 *out, REAL8 x);
REAL8 twospect_cheb_eval(const cheb_series * cs, REAL8 x);
REAL8 gammastar_ser(REAL8 x);
REAL8 sf_exprel_n_CF(REAL8 N, REAL8 x);
REAL8 gamma_inc_Q_asymp_unif(REAL8 a, REAL8 x);
INT4 gamma_inc_Q_CF(REAL8 *out, REAL8 a, REAL8 x);
INT4 gamma_inc_F_CF(REAL8 *out, REAL8 a, REAL8 x);
REAL8 gamma_inc_Q_large_x(REAL8 a, REAL8 x);
REAL8 epsval(REAL8 val);
REAL8 binodeviance(REAL8 x, REAL8 np);

REAL4 epsval_float(REAL4 val);

void sumseries(REAL8 *computedprob, REAL8 P, REAL8 C, REAL8 E, INT8 counter, REAL8 x, REAL8 dof, REAL8 halfdelta, REAL8 err, INT4 countdown);
void sumseries_eg(REAL8 *computedprob, REAL8 P, REAL8 C, REAL8 E, INT8 counter, REAL8 x, REAL8 dof, REAL8 halfdelta, REAL8 err, INT4 countdown);

#endif
