/*
 *  Copyright (C) 2014 Evan Goetz
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

#ifndef __CDFDIST_H__
#define __CDFDIST_H__

#include <lal/LALStdlib.h>

REAL8 cdf_chisq_Pinv(REAL8 P, REAL8 nu);
REAL8 cdf_chisq_Qinv(REAL8 Q, REAL8 nu);
REAL8 cdf_gamma_Pinv(REAL8 P, REAL8 a, REAL8 b);
INT4 cdf_gamma_Qinv(REAL8 *out, REAL8 Q, REAL8 a, REAL8 b);
REAL8 cdf_ugaussian_Pinv(REAL8 P);
INT4 cdf_ugaussian_Qinv(REAL8 *out, REAL8 Q);
REAL8 cdf_gamma_P(REAL8 x, REAL8 a, REAL8 b);
REAL8 cdf_gamma_P_usingmatlab(REAL8 x, REAL8 a, REAL8 b);
INT4 cdf_gamma_Q(REAL8 *out, REAL8 x, REAL8 a, REAL8 b);
REAL8 cdf_gamma_Q_usingmatlab(REAL8 x, REAL8 a, REAL8 b);
REAL8 ncx2cdf(REAL8 x, REAL8 dof, REAL8 delta);
REAL8 ncx2cdf_withouttinyprob(REAL8 x, REAL8 dof, REAL8 delta);
REAL8 ncx2cdf_withouttinyprob_withmatlabchi2cdf(REAL8 x, REAL8 dof, REAL8 delta);
REAL8 ncx2pdf(REAL8 x, REAL8 dof, REAL8 delta);
REAL8 ncx2inv(REAL8 p, REAL8 dof, REAL8 delta);
REAL4 ncx2inv_float(REAL8 p, REAL8 dof, REAL8 delta);
REAL8 norminv(REAL8 p, REAL8 mu, REAL8 sigma);
REAL8 twospect_cdf_chisq_P(REAL8 x, REAL8 nu);
REAL8 matlab_cdf_chisq_P(REAL8 x, REAL8 nu);
REAL8 unitGaussianSNR(REAL8 value, REAL8 dof);
REAL4 ncx2cdf_float(REAL4 x, REAL4 dof, REAL4 delta);
REAL4 ncx2cdf_float_withouttinyprob(REAL4 x, REAL4 dof, REAL4 delta);
REAL4 ncx2cdf_float_withouttinyprob_withmatlabchi2cdf(REAL4 x, REAL4 dof, REAL4 delta);

#endif
