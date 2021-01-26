
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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */

 /*
    This file contains the functions for the fits of the ringdown and damping frequencies for each mode
    and the fits for the real and imaginary part of the mixing coefficients (only for the 32 mode).
    The fits build on data taken from https://pages.jh.edu/~eberti2/ringdown/, see also https://arxiv.org/abs/gr-qc/0512160 and https://arxiv.org/abs/1408.1860. Explicit expressions are given in the supplementary material, in IMRPhenomXHM_PhaseFits.nb, as well as in the dedicated folders QNMs and MixingCoefficients, see https://dcc.ligo.org/LIGO-P2000011 and https://arxiv.org/abs/2001.10914
 */

#include "LALSimIMRPhenomXHM_qnm.h"


static double evaluate_QNMfit_fring21(double finalDimlessSpin){

  double return_val;

  if (fabs(finalDimlessSpin) > 1.0) {
    XLAL_ERROR(XLAL_EDOM, "PhenomXHM evaluate_QNMfit_fring21 function: |finalDimlessSpin| > 1.0 not supported");
  }

  double x2= finalDimlessSpin*finalDimlessSpin;
  double x3= x2*finalDimlessSpin;
  double x4= x2*x2;
  double x5= x3*x2;

  return_val = (0.059471695665734674 - 0.07585416297991414*finalDimlessSpin + 0.021967909664591865*x2 - 0.0018964744613388146*x3 + 0.001164879406179587*x4 - 0.0003387374454044957*x5)/(1 - 1.4437415542456158*finalDimlessSpin + 0.49246920313191234*x2);
  return return_val;
}

static double evaluate_QNMfit_fdamp21(double finalDimlessSpin){

  double return_val;

  if (fabs(finalDimlessSpin) > 1.0) {
    XLAL_ERROR(XLAL_EDOM, "PhenomXHM evaluate_QNMfit_fdamp21 function: |finalDimlessSpin| > 1.0 not supported");
  }

  double x2= finalDimlessSpin*finalDimlessSpin;
  double x3= x2*finalDimlessSpin;
  double x4= x2*x2;
  double x5= x3*x2;

  return_val = (2.0696914454467294 - 3.1358071947583093*finalDimlessSpin + 0.14456081596393977*x2 + 1.2194717985037946*x3 - 0.2947372598589144*x4 + 0.002943057145913646*x5)/(146.1779212636481 - 219.81790388304876*finalDimlessSpin + 17.7141194900164*x2 + 75.90115083917898*x3 - 18.975287709794745*x4);
  return return_val;
}

static double evaluate_QNMfit_fring33(double finalDimlessSpin){

  double return_val;

  if (fabs(finalDimlessSpin) > 1.0) {
    XLAL_ERROR(XLAL_EDOM, "PhenomXHM evaluate_QNMfit_fring33 function: |finalDimlessSpin| > 1.0 not supported");
  }

  double x2= finalDimlessSpin*finalDimlessSpin;
  double x3= x2*finalDimlessSpin;
  double x4= x2*x2;
  double x5= x3*x2;
  double x6= x3*x3;

  return_val = (0.09540436245212061 - 0.22799517865876945*finalDimlessSpin + 0.13402916709362475*x2 + 0.03343753057911253*x3 - 0.030848060170259615*x4 - 0.006756504382964637*x5 + 0.0027301732074159835*x6)/(1 - 2.7265947806178334*finalDimlessSpin + 2.144070539525238*x2 - 0.4706873667569393*x4 + 0.05321818246993958*x6);
  return return_val;
}

static double evaluate_QNMfit_fdamp33(double finalDimlessSpin){

  double return_val;

  if (fabs(finalDimlessSpin) > 1.0) {
    XLAL_ERROR(XLAL_EDOM, "PhenomXHM evaluate_QNMfit_fdamp33 function: |finalDimlessSpin| > 1.0 not supported");
  }

  double x2= finalDimlessSpin*finalDimlessSpin;
  double x3= x2*finalDimlessSpin;
  double x4= x2*x2;
  double x5= x3*x2;

  return_val = (0.014754148319335946 - 0.03124423610028678*finalDimlessSpin + 0.017192623913708124*x2 + 0.001034954865629645*x3 - 0.0015925124814622795*x4 - 0.0001414350555699256*x5)/(1 - 2.0963684630756894*finalDimlessSpin + 1.196809702382645*x2 - 0.09874113387889819*x4);
  return return_val;
}

static double evaluate_QNMfit_fring32(double finalDimlessSpin){

  double return_val;

  if (fabs(finalDimlessSpin) > 1.0) {
    XLAL_ERROR(XLAL_EDOM, "PhenomXHM evaluate_QNMfit_fring32 function: |finalDimlessSpin| > 1.0 not supported");
  }

  double x2= finalDimlessSpin*finalDimlessSpin;
  double x3= x2*finalDimlessSpin;
  double x4= x2*x2;
  double x5= x3*x2;
  double x6= x3*x3;

  return_val = (0.09540436245212061 - 0.13628306966373951*finalDimlessSpin + 0.030099881830507727*x2 - 0.000673589757007597*x3 + 0.0118277880067919*x4 + 0.0020533816327907334*x5 - 0.0015206141948469621*x6)/(1 - 1.6531854335715193*finalDimlessSpin + 0.5634705514193629*x2 + 0.12256204148002939*x4 - 0.027297817699401976*x6);
  return return_val;
}

static double evaluate_QNMfit_fdamp32(double finalDimlessSpin){

  double return_val;

  if (fabs(finalDimlessSpin) > 1.0) {
    XLAL_ERROR(XLAL_EDOM, "PhenomXHM evaluate_QNMfit_fdamp32 function: |finalDimlessSpin| > 1.0 not supported");
  }

  double x2= finalDimlessSpin*finalDimlessSpin;
  double x3= x2*finalDimlessSpin;
  double x4= x2*x2;

  return_val = (0.014754148319335946 - 0.03445752346074498*finalDimlessSpin + 0.02168855041940869*x2 + 0.0014945908223317514*x3 - 0.0034761714223258693*x4)/(1 - 2.320722660848874*finalDimlessSpin + 1.5096146036915865*x2 - 0.18791187563554512*x4);
  return return_val;
}

static double evaluate_QNMfit_fring44(double finalDimlessSpin){

  double return_val;

  if (fabs(finalDimlessSpin) > 1.0) {
    XLAL_ERROR(XLAL_EDOM, "PhenomXHM evaluate_QNMfit_fring44 function: |finalDimlessSpin| > 1.0 not supported");
  }

  double x2= finalDimlessSpin*finalDimlessSpin;
  double x3= x2*finalDimlessSpin;
  double x4= x2*x2;
  double x5= x3*x2;
  double x6= x3*x3;

  return_val = (0.1287821193485683 - 0.21224284094693793*finalDimlessSpin + 0.0710926778043916*x2 + 0.015487322972031054*x3 - 0.002795401084713644*x4 + 0.000045483523029172406*x5 + 0.00034775290179000503*x6)/(1 - 1.9931645124693607*finalDimlessSpin + 1.0593147376898773*x2 - 0.06378640753152783*x4);
  return return_val;
}

static double evaluate_QNMfit_fdamp44(double finalDimlessSpin){

  double return_val;

  if (fabs(finalDimlessSpin) > 1.0) {
    XLAL_ERROR(XLAL_EDOM, "PhenomXHM evaluate_QNMfit_fdamp44 function: |finalDimlessSpin| > 1.0 not supported");
  }

  double x2= finalDimlessSpin*finalDimlessSpin;
  double x3= x2*finalDimlessSpin;
  double x4= x2*x2;
  double x5= x3*x2;
  double x6= x3*x3;

  return_val = (0.014986847152355699 - 0.01722587715950451*finalDimlessSpin - 0.0016734788189065538*x2 + 0.0002837322846047305*x3 + 0.002510528746148588*x4 + 0.00031983835498725354*x5 + 0.000812185411753066*x6)/(1 - 1.1350205970682399*finalDimlessSpin - 0.0500827971270845*x2 + 0.13983808071522857*x4 + 0.051876225199833995*x6);
  return return_val;
}


static double evaluate_QNMfit_re_l2m2lp2(double finalDimlessSpin){

  double return_val;

  if (fabs(finalDimlessSpin) > 1.0) {
    XLAL_ERROR(XLAL_EDOM, "PhenomXHM evaluate_QNMfit_re_l2m2lp2 function: |finalDimlessSpin| > 1.0 not supported");
  }

  double x2= finalDimlessSpin*finalDimlessSpin;
  double x3= x2*finalDimlessSpin;
  double x4= x2*x2;
  double x5= x3*x2;
  double x6= x3*x3;

  return_val = (1 - 2.2956993576253635*finalDimlessSpin + 1.461988775298876*x2 + 0.0043296365593147035*x3 - 0.1695667458204109*x4 - 0.0006267849034466508*x5)/(1 - 2.2956977727459043*finalDimlessSpin + 1.4646339137818438*x2 - 0.16843226886562457*x4 - 0.00007150540890128118*x6);
  return return_val;
}

static double evaluate_QNMfit_im_l2m2lp2(double finalDimlessSpin){

  double return_val;

  if (fabs(finalDimlessSpin) > 1.0) {
    XLAL_ERROR(XLAL_EDOM, "PhenomXHM evaluate_QNMfit_im_l2m2lp2 function: |finalDimlessSpin| > 1.0 not supported");
  }

  double x2= finalDimlessSpin*finalDimlessSpin;
  double x3= x2*finalDimlessSpin;
  double x4= x2*x2;
  double x5= x3*x2;
  double x6= x3*x3;

  return_val = (finalDimlessSpin*(0.3826673013161342 - 0.47531267226013896*finalDimlessSpin - 0.05898102880105067*x2 + 0.0724525431346487*x3 + 0.054714637311702986*x4 + 0.024544862718252784*x5))/(-38.70835035062785 + 69.82140084545878*finalDimlessSpin - 27.99036444363243*x2 - 4.152310472191899*x4 + 1.*x6);
  return return_val;
}

static double evaluate_QNMfit_re_l3m2lp2(double finalDimlessSpin){

  double return_val;

  if (fabs(finalDimlessSpin) > 1.0) {
    XLAL_ERROR(XLAL_EDOM, "PhenomXHM evaluate_QNMfit_re_l3m2lp2 function: |finalDimlessSpin| > 1.0 not supported");
  }

  double x2= finalDimlessSpin*finalDimlessSpin;
  double x3= x2*finalDimlessSpin;
  double x4= x2*x2;
  double x5= x3*x2;

  return_val = (finalDimlessSpin*(0.47513455283841244 - 0.9016636384605536*finalDimlessSpin + 0.3844811236426182*x2 + 0.0855565148647794*x3 - 0.03620067426672167*x4 - 0.006557249133752502*x5))/(-6.76894063440646 + 15.170831931186493*finalDimlessSpin - 9.406169787571082*x2 + 1.*x4);
  return return_val;
}

static double evaluate_QNMfit_im_l3m2lp2(double finalDimlessSpin){

  double return_val;

  if (fabs(finalDimlessSpin) > 1.0) {
    XLAL_ERROR(XLAL_EDOM, "PhenomXHM evaluate_QNMfit_im_l3m2lp2 function: |finalDimlessSpin| > 1.0 not supported");
  }

  double x2= finalDimlessSpin*finalDimlessSpin;
  double x3= x2*finalDimlessSpin;
  double x4= x2*x2;
  double x5= x3*x2;
  double x6= x3*x3;

  return_val = (finalDimlessSpin*(-2.8704762147145533 + 4.436434016918535*finalDimlessSpin - 1.0115343326360486*x2 - 0.08965314412106505*x3 - 0.4236810894599512*x4 - 0.041787576033810676*x5))/(-171.80908957903395 + 272.362882450877*finalDimlessSpin - 76.68544453077854*x2 - 25.14197656531123*x4 + 1.*x6);
  return return_val;
}

static double evaluate_QNMfit_re_l2m2lp3(double finalDimlessSpin){

  double return_val;

  if (fabs(finalDimlessSpin) > 1.0) {
    XLAL_ERROR(XLAL_EDOM, "PhenomXHM evaluate_QNMfit_re_l2m2lp3 function: |finalDimlessSpin| > 1.0 not supported");
  }

  double x2= finalDimlessSpin*finalDimlessSpin;
  double x3= x2*finalDimlessSpin;
  double x4= x2*x2;
  double x5= x3*x2;
  double x6= x3*x3;

  return_val = (finalDimlessSpin*(18.522563276099167 - 37.978140351289014*finalDimlessSpin + 19.030390708998894*x2 + 3.0355668591803386*x3 - 2.210028290847915*x4 - 0.37117112862247975*x5))/(164.52480238697507 - 377.9093045285145*finalDimlessSpin + 243.3353695550844*x2 - 30.79738566181734*x4 + 1.*x6);
  return return_val;
}

static double evaluate_QNMfit_im_l2m2lp3(double finalDimlessSpin){

  double return_val;

  if (fabs(finalDimlessSpin) > 1.0) {
    XLAL_ERROR(XLAL_EDOM, "PhenomXHM evaluate_QNMfit_im_l2m2lp3 function: |finalDimlessSpin| > 1.0 not supported");
  }

  double x2= finalDimlessSpin*finalDimlessSpin;
  double x3= x2*finalDimlessSpin;
  double x4= x2*x2;
  double x5= x3*x2;
  double x6= x3*x3;

  return_val = (finalDimlessSpin*(-49.7688437256778 + 120.43773704442333*finalDimlessSpin - 82.95323455645332*x2 + 1.721453011852496*x3 + 11.540237244397877*x4 - 0.9819458637589314*x5))/(2858.5790831181725 - 6305.619505422591*finalDimlessSpin + 3825.6742092829054*x2 - 377.7822297815406*x4 + 1.*x6);
  return return_val;
}

static double evaluate_QNMfit_re_l3m2lp3(double finalDimlessSpin){

  double return_val;

  if (fabs(finalDimlessSpin) > 1.0) {
    XLAL_ERROR(XLAL_EDOM, "PhenomXHM evaluate_QNMfit_re_l3m2lp3 function: |finalDimlessSpin| > 1.0 not supported");
  }

  double x2= finalDimlessSpin*finalDimlessSpin;
  double x3= x2*finalDimlessSpin;
  double x4= x2*x2;
  double x5= x3*x2;
  double x6= x3*x3;

  return_val = (1 - 2.107852425643677*finalDimlessSpin + 1.1906393634562715*x2 + 0.02244848864087732*x3 - 0.09593447799423722*x4 - 0.0021343381708933025*x5 - 0.005319515989331159*x6)/(1 - 2.1078515887706324*finalDimlessSpin + 1.2043484690080966*x2 - 0.08910191596778137*x4 - 0.005471749827809503*x6);
  return return_val;
}

static double evaluate_QNMfit_im_l3m2lp3(double finalDimlessSpin){

  double return_val;

  if (fabs(finalDimlessSpin) > 1.0) {
    XLAL_ERROR(XLAL_EDOM, "PhenomXHM evaluate_QNMfit_im_l3m2lp3 function: |finalDimlessSpin| > 1.0 not supported");
  }

  double x2= finalDimlessSpin*finalDimlessSpin;
  double x3= x2*finalDimlessSpin;
  double x4= x2*x2;
  double x5= x3*x2;
  double x6= x3*x3;

  return_val = (finalDimlessSpin*(12.45701482868677 - 29.398484595717147*finalDimlessSpin + 18.26221675782779*x2 + 1.9308599142669403*x3 - 3.159763242921214*x4 - 0.0910871567367674*x5))/(345.52914639836257 - 815.4349339779621*finalDimlessSpin + 538.3888932415709*x2 - 69.3840921447381*x4 + 1.*x6);
  return return_val;
}
