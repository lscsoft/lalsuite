/*
 * Copyright (C) 2019  Marta Colleoni, Cecilio García Quirós
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
 *
 */

#ifndef _LALSIM_IMR_PHENOMXHM_QNM_H
#define _LALSIM_IMR_PHENOMXHM_QNM_H

#ifdef __cplusplus
extern "C" {
#endif

// Declaration of the functions for the ringdown and damping frequency fits for each mode.
static double evaluate_QNMfit_fring21(double finalDimlessSpin);
static double evaluate_QNMfit_fdamp21(double finalDimlessSpin);
static double evaluate_QNMfit_fring33(double finalDimlessSpin);
static double evaluate_QNMfit_fdamp33(double finalDimlessSpin);
static double evaluate_QNMfit_fring32(double finalDimlessSpin);
static double evaluate_QNMfit_fdamp32(double finalDimlessSpin);
static double evaluate_QNMfit_fring44(double finalDimlessSpin);
static double evaluate_QNMfit_fdamp44(double finalDimlessSpin);

// Declaration of the functions for the mixing coefficients (only 32 mode)
static double evaluate_QNMfit_re_l2m2lp2(double finalDimlessSpin);
static double evaluate_QNMfit_im_l2m2lp2(double finalDimlessSpin);
static double evaluate_QNMfit_re_l3m2lp2(double finalDimlessSpin);
static double evaluate_QNMfit_im_l3m2lp2(double finalDimlessSpin);
static double evaluate_QNMfit_re_l2m2lp3(double finalDimlessSpin);
static double evaluate_QNMfit_im_l2m2lp3(double finalDimlessSpin);
static double evaluate_QNMfit_re_l3m2lp3(double finalDimlessSpin);
static double evaluate_QNMfit_im_l3m2lp3(double finalDimlessSpin);

#ifdef __cplusplus
}
#endif

#endif
