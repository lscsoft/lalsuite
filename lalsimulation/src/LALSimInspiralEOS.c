/*
 *  Copyright (C) 2012 Walter Del Pozzo, Tjonnie Li, Michalis Agathos
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
 #include <stdlib.h> 
 #include <stdio.h>
 #include <string.h>
 #include <lal/LALSimInspiralEOS.h>
 #include <lal/LALSimInspiral.h>

 LALEquationOfState XLALSimEOSfromString(char eos_name[])
 {
   LALEquationOfState eos;
   if (!strcmp("MS1",eos_name)) eos = LAL_SIM_INSPIRAL_EOS_MS1;
   else if (!strcmp("H4",eos_name)) eos = LAL_SIM_INSPIRAL_EOS_H4;
   else if (!strcmp("SQM3",eos_name)) eos = LAL_SIM_INSPIRAL_EOS_SQM3;
   else if (!strcmp("MPA1",eos_name)) eos = LAL_SIM_INSPIRAL_EOS_MPA1;
   else if (!strcmp("GNH3",eos_name)) eos = LAL_SIM_INSPIRAL_EOS_GNH3;
   else if (!strcmp("A",eos_name)) eos = LAL_SIM_INSPIRAL_EOS_A;
   else if (!strcmp("AU",eos_name)) eos = LAL_SIM_INSPIRAL_EOS_AU;
   else if (!strcmp("FPS",eos_name)) eos = LAL_SIM_INSPIRAL_EOS_FPS;
   else if (!strcmp("APR",eos_name)) eos = LAL_SIM_INSPIRAL_EOS_APR;
   else if (!strcmp("UU",eos_name)) eos = LAL_SIM_INSPIRAL_EOS_UU;
   else if (!strcmp("L",eos_name)) eos = LAL_SIM_INSPIRAL_EOS_L;
   else if (!strcmp("PP",eos_name)) eos = LAL_SIM_INSPIRAL_EOS_NONE;
   else
   {
    XLALPrintError( "XLAL Error - %s: Equation of state %s not recognized.", __func__, eos_name);
    XLAL_ERROR(XLAL_EINVAL);
   }
   return eos;
 }

REAL8 XLALSimInspiralEOSLambda(LALEquationOfState eos_type, REAL8 m_intr_msun){/** this must be fed the INTRINSIC mass */

    /* this is fed the intrinsic masses and then computes the value of \Lambda(m) See Hinderer et al ( http://arxiv.org/abs/0911.3535 ) for details of the EOSes*/
    /* \Lambda(m) is in units of s^5 */
    REAL8 lambda=0.;
//  printf("EOS number: %d\n", eos_type);
//  printf("mass: %e\n", m_intr_msun);
    switch (eos_type)
    {
        case LAL_SIM_INSPIRAL_EOS_NONE:
            lambda = 0.0;
        break;
    // MS1
        case LAL_SIM_INSPIRAL_EOS_MS1:
           // printf("Using EOS MS1\n");
            lambda = 2.755956E-24*(2.19296 + 20.0273*m_intr_msun - 17.9443*m_intr_msun*m_intr_msun 
            + 5.75129*m_intr_msun*m_intr_msun*m_intr_msun - 0.699095*m_intr_msun*m_intr_msun*m_intr_msun*m_intr_msun);
        break;
    // H4
        case LAL_SIM_INSPIRAL_EOS_H4:
            lambda = 2.755956E-24*(0.743351 + 15.8917*m_intr_msun - 14.7348*m_intr_msun*m_intr_msun 
            + 5.32863*m_intr_msun*m_intr_msun*m_intr_msun - 0.942625*m_intr_msun*m_intr_msun*m_intr_msun*m_intr_msun);
        break; 
    // SQM3
        case LAL_SIM_INSPIRAL_EOS_SQM3:
            lambda = 2.755956E-24*(-5.55858 + 20.8977*m_intr_msun - 20.5583*m_intr_msun*m_intr_msun 
            + 9.55465*m_intr_msun*m_intr_msun*m_intr_msun - 1.84933*m_intr_msun*m_intr_msun*m_intr_msun*m_intr_msun);
        break;
    // MPA1
    case LAL_SIM_INSPIRAL_EOS_MPA1:
        lambda = 2.755956E-24*(0.276761 + 7.26925*m_intr_msun - 5.72102*m_intr_msun*m_intr_msun
        + 1.51347*m_intr_msun*m_intr_msun*m_intr_msun - 0.152181*m_intr_msun*m_intr_msun*m_intr_msun*m_intr_msun);
        break;
    // GNH3
    case LAL_SIM_INSPIRAL_EOS_GNH3:
        lambda = 2.755956E-24*(7.80715 + 0.683549*m_intr_msun + 1.21351*m_intr_msun*m_intr_msun
        - 3.50234*m_intr_msun*m_intr_msun*m_intr_msun + 0.894662*m_intr_msun*m_intr_msun*m_intr_msun*m_intr_msun);
        break;
    default:
        lambda = 0.0;
        break;
    }
//  printf("calculated love number: %e\n", lambda);
    if (lambda<0.0) return 0.0;
    else return lambda;
}

REAL8 XLALLambdaQuadratic(REAL8 c0, REAL8 c1, REAL8 c2, REAL8 mass) {
    mass = mass*LAL_MTSUN_SI;
    // [LAMBDA0] = SEC^5; [LAMBDA1] = SEC^4; [LAMBDA2] = SEC^3
    REAL8 lambda = 1.0E-23*c0 + 1.0E-18*(mass-1.4*LAL_MTSUN_SI)*c1 + 1.0E-13*(mass-1.4*LAL_MTSUN_SI)*(mass-1.4*LAL_MTSUN_SI)*c2;
    lambda = (lambda > 0.0) ? lambda : 0.0;
    return lambda;
}
 

REAL8 XLALSimInspiralEOSQfromLambda(REAL8 lambda) {
    /* Quadrupole-monopole parameter calculated from love number;
       see http://arxiv.org/abs/1303.1528 */
    REAL8 q, loglam;
    REAL8 tolerance = 1E-15;
    if(lambda<tolerance) { //printf("Love number is (nearly) zero; cannot compute QM parameter. Setting to 1.0 (BH value).\n");
                      q = 1.0; } 
    else {
    loglam = log(lambda);
    q =  0.194 + 0.0936*loglam + 0.0474*loglam*loglam;
    q -= 0.00421*loglam*loglam*loglam;
    q += 0.000123*loglam*loglam*loglam*loglam;
    q = exp(q);
    }

//  printf("%e %e\n", l, q); // Testing numerical results from these functions.

    return q;

}

REAL8 XLALSimInspiralEOSqmparameter(LALEquationOfState eos_type, REAL8 m_intr_msun){
  
  REAL8 q = 0.0 ;
  REAL8 m = m_intr_msun ;
  REAL8 m2 = m*m ;
  REAL8 m3 = m2*m ;
  
  switch (eos_type) {
  /*  */
  case LAL_SIM_INSPIRAL_EOS_A:
    q = -6.41414141*m3 + 30.70779221*m2 - 53.37417027*m + 35.62253247 ;
    break;
  /*  */
  case LAL_SIM_INSPIRAL_EOS_AU:
    q = -6.18686869*m3 + 30.15909091*m2 - 52.87806638*m + 35.86616883 ;
    break;
  /*  */
  case LAL_SIM_INSPIRAL_EOS_FPS:
    q = -3.86363636*m3 + 21.03030303*m2 - 42.19448052*m + 32.83722944 ;
    break;
  /*  */
  case LAL_SIM_INSPIRAL_EOS_APR:
    q = -10.55555556*m3 + 49.52380952*m2 - 82.77063492*m + 53.02428571 ;
    break;
  /*  */
  case LAL_SIM_INSPIRAL_EOS_UU:
    q = -8.03030303*m3 + 37.61363636*m2 - 63.48733766*m + 41.75080087 ;
    break;
  /*  */
  case LAL_SIM_INSPIRAL_EOS_L:
    q = -6.59090909*m3 + 33.67424242*m2 - 63.77034632*m + 48.98073593 ;
    break;
  case LAL_SIM_INSPIRAL_EOS_NONE:
    q = 1.0 ;
    break;

  default:
    q = 1.0 ;
    break ;
  }
  
  if (q < 1.0) {
    q = 1.0;
  }
  
  return q ;
}

/**
 * This function estimates the radius of a NS of a given mass and
 * tidal deformability parameter, based on the "I-Love-Q forever" relation
 * of Maselli et al, arXiv:1304.2052v1.
 * To be used for masses within [1.2,2]M_sun, and preferably not for strange
 * quark stars (since the relation is calibrated for this mass range and for 
 * the EoS APR4, MS1, H4).
 * For a BH, (lambda=0) it returns the Schwarzschild radius.
 * The arguments are:
 * m_intr_msun              the intrinsic mass in solar masses
 * barlambda                the dimensionless tidal deformability (lambda/m^5)
 * The return value is the radius in meters.
 */

REAL8 XLALSimInspiralNSRadiusOfLambdaM(REAL8 m_intr_msun, REAL8 barlambda){

  REAL8 loglambda;
  REAL8 compactness, radius ;
  REAL8 tolerance = 1E-15;

  /* Check for sign of lambda */
  if ( barlambda <= tolerance && barlambda >= 0.0 ) {
  /* This is a black hole */
    compactness = 0.5;
  }
  else if ( barlambda > tolerance ) {
  loglambda = log(barlambda);
  /* Calculate compactness according to arXiv:1304.2052v1 */
  compactness = 0.371 - 0.0391*loglambda + 0.001056*loglambda*loglambda;
  }
  else {
    XLALPrintError( "XLAL Error - %s: Tidal deformability is negative. Only positive values are acceptable.", __func__);
    XLAL_ERROR_REAL8(XLAL_EDOM);
  }

  /* Check that radius is larger than Schwarzschild radius */
  if ( compactness > 0.5 ) {
    XLALPrintWarning( "XLAL Warning - %s: Neutron Star is calculated to have compactness larger than a black hole (C = %f, lambda = %f, m = %f).\n Setting C=0.5 ...", __func__, compactness, barlambda, m_intr_msun);
    compactness = 0.5;
  }

  if ( compactness < 0.0 ) {
    XLALPrintError( "XLAL Error - %s: Neutron Star is calculated to have negative compactness (C = %f, lambda = %f, m = %f).", __func__, compactness, barlambda, m_intr_msun);
    XLAL_ERROR_REAL8(XLAL_ERANGE);
  }

  radius = LAL_MRSUN_SI * m_intr_msun / compactness;

 return radius;

}


/**
 * This function estimates the radius for a binary of given masses and
 * tidal deformability parameters.
 * It uses XLALSimInspiralNSRadiusOfLambdaM() to calculate radii (see above).
 * The arguments are:
 * m1_intr, m2_intr                      the intrinsic masses in solar masses
 * barlambda1, barlambda2                the dimensionless tidal deformabilities (lambda_i/m_i^5)
 * The return value is the GW contact frequency in Hz.
 */

REAL8 XLALSimInspiralContactFrequency(REAL8 m1_intr, REAL8 barlambda1, REAL8 m2_intr, REAL8 barlambda2){

  REAL8 r1, r2, rtot, mtot, f_gw_contact;

  /* Calculate radii for the two components */
  r1 = XLALSimInspiralNSRadiusOfLambdaM(m1_intr, barlambda1);
  r2 = XLALSimInspiralNSRadiusOfLambdaM(m2_intr, barlambda2);

  rtot = (r1 + r2)/LAL_C_SI;                             // Orbital distance in seconds
  mtot = (m1_intr + m2_intr)*LAL_MTSUN_SI;               // Total mass in seconds

  /* Calculate the GW contact frequency */
  f_gw_contact = sqrt(mtot/(rtot*rtot*rtot))/LAL_PI;
  if ( f_gw_contact < 0.0 ) {
    XLALPrintError( "XLAL Error - %s: Contact frequency is calculated to be negative  (fcontact = %f)", __func__, f_gw_contact);
    XLAL_ERROR_REAL8(XLAL_ERANGE);
  }

  return f_gw_contact;

}
