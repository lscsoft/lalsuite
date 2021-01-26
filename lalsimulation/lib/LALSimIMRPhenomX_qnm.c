/*
 * Copyright (C) 2018 Sascha Husa, Geraint Pratten
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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */


/*
 * \author Sascha Husa
 * \author Geraint Pratten
 *
 */

#include "LALSimIMRPhenomX_qnm.h"

#if QNMfits == 1
/**
 *  evaluate fit QNMData_fring_22
 */
static double evaluate_QNMfit_fring22(double finalDimlessSpin){

double return_val;

if (fabs(finalDimlessSpin) > 1.0) {
       XLAL_ERROR(XLAL_EDOM, "PhenomX evaluate_QNMfit_fring22 \
function: |finalDimlessSpin| > 1.0 not supported");
}

double x2= finalDimlessSpin*finalDimlessSpin;
double x3= x2*finalDimlessSpin;
double x4= x2*x2;
double x5= x3*x2;
double x6= x3*x3;
double x7= x4*x3;

return_val = (0.05947169566573468 - \
0.14989771215394762*finalDimlessSpin + 0.09535606290986028*x2 + \
0.02260924869042963*x3 - 0.02501704155363241*x4 - \
0.005852438240997211*x5 + 0.0027489038393367993*x6 + \
0.0005821983163192694*x7)/(1 - 2.8570126619966296*finalDimlessSpin + \
2.373335413978394*x2 - 0.6036964688511505*x4 + \
0.0873798215084077*x6);
return return_val;
}


/**
 *  evaluate fit QNMData_fdamp_22
 */
static double evaluate_QNMfit_fdamp22(double finalDimlessSpin){

double return_val;

if (fabs(finalDimlessSpin) > 1.0) {
       XLAL_ERROR(XLAL_EDOM, "PhenomX evaluate_QNMfit_fdamp22 \
function: |finalDimlessSpin| > 1.0 not supported");
}

double x2= finalDimlessSpin*finalDimlessSpin;
double x3= x2*finalDimlessSpin;
double x4= x2*x2;
double x5= x3*x2;
double x6= x3*x3;

return_val = (0.014158792290965177 - \
0.036989395871554566*finalDimlessSpin + 0.026822526296575368*x2 + \
0.0008490933750566702*x3 - 0.004843996907020524*x4 - \
0.00014745235759327472*x5 + 0.0001504546201236794*x6)/(1 - \
2.5900842798681376*finalDimlessSpin + 1.8952576220623967*x2 - \
0.31416610693042507*x4 + 0.009002719412204133*x6);
return return_val;
}


#else
 /**
  *  evaluate interpolated data set QNMData_fring_22, see Berti et al, CQG, 26, 163001, (2009)
  */
static double interpolateQNMData_fring_22(double finalDimlessSpin) {
   double return_val;

   if (fabs(finalDimlessSpin) > 1.0) {
 	   XLAL_ERROR(XLAL_EDOM, "PhenomX interpolateQNMData_fring_22 function: |finalDimlessSpin| > 1.0 not supported");
   }
   gsl_interp_accel *acc = gsl_interp_accel_alloc();
   gsl_spline *iData = gsl_spline_alloc(gsl_interp_cspline, QNMData_fring_22_length);
   gsl_spline_init(iData, QNMData_a, QNMData_fring_22, QNMData_fring_22_length);

   return_val = gsl_spline_eval(iData, finalDimlessSpin, acc);

   gsl_spline_free(iData);
   gsl_interp_accel_free(acc);
   return return_val;
 }


 /**
  *  evaluate interpolated data set QNMData_fdamp_22, see Berti et al, CQG, 26, 163001, (2009)
  */
static double interpolateQNMData_fdamp_22(double finalDimlessSpin) {
   double return_val;

   if (fabs(finalDimlessSpin) > 1.0) {
 	XLAL_ERROR(XLAL_EDOM, "PhenomX interpolateQNMData_fdamp_22 function: |finalDimlessSpin| > 1.0 not supported");
   }
   gsl_interp_accel *acc = gsl_interp_accel_alloc();
   gsl_spline *iData = gsl_spline_alloc(gsl_interp_cspline, QNMData_fdamp_22_length);
   gsl_spline_init(iData, QNMData_a, QNMData_fdamp_22, QNMData_fdamp_22_length);

   return_val = gsl_spline_eval(iData, finalDimlessSpin, acc);

   gsl_spline_free(iData);
   gsl_interp_accel_free(acc);
   return return_val;
 }

 #endif
