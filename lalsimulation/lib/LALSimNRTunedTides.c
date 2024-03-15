/*
 * Copyright (C) 2017 Tim Dietrich, Sebastiano Bernuzzi, Nathan Johnson-McDaniel,
 * Shasvath J Kapadia, Francesco Pannarale and Sebastian Khan, Michael Puerrer, Adrian Abac.
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/LALSimIMR.h>

#include "LALSimNRTunedTides.h"
#include "LALSimUniversalRelations.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/**
 * Planck taper window
 */
static REAL8 PlanckTaper(const REAL8 t, const REAL8 t1, const REAL8 t2) {
  REAL8 taper;
  if (t <= t1)
    taper = 0.;
  else if (t >= t2)
    taper = 1.;
  else
    taper = 1. / (exp((t2 - t1)/(t - t1) + (t2 - t1)/(t - t2)) + 1.);

  return taper;
}

/**
 * function to swap masses and lambda to enforce m1 >= m2
 */
static int EnforcePrimaryMassIsm1(REAL8 *m1, REAL8 *m2, REAL8 *lambda1, REAL8 *lambda2){
  if ((*m1 == *m2) && (*lambda1 != *lambda2))
    XLALPrintWarning("m1 == m2 (%g), but lambda1 != lambda2 (%g, %g).\n", *m1, *lambda1, *lambda2);

  REAL8 lambda1_tmp, lambda2_tmp, m1_tmp, m2_tmp;
  if (*m1>=*m2) {
    lambda1_tmp = *lambda1;
    lambda2_tmp = *lambda2;
    m1_tmp   = *m1;
    m2_tmp   = *m2;
  } else { 
    lambda1_tmp = *lambda2;
    lambda2_tmp = *lambda1;
    m1_tmp   = *m2;
    m2_tmp   = *m1;
  }
  *m1 = m1_tmp;
  *m2 = m2_tmp;
  *lambda1 = lambda1_tmp;
  *lambda2 = lambda2_tmp;

  if (*m1 < *m2)
    XLAL_ERROR(XLAL_EDOM, "XLAL_ERROR in EnforcePrimaryMassIsm1. When trying\
 to enfore that m1 should be the larger mass.\
 After trying to enforce this m1 = %f and m2 = %f\n", *m1, *m2);

  return XLAL_SUCCESS;
}

/**
 * function to swap masses, spins, and lambda to enforce m1 >= m2, This is mainly for NRTidalv3, 
 * which has a merger frequency fit that is dependent on the aligned spin component. See Eq. (41)
 * in https://arxiv.org/pdf/2311.07456.pdf.
 */
static int EnforcePrimaryMassIsm1_v3(REAL8 *m1, REAL8 *m2, REAL8 *lambda1, REAL8 *lambda2, REAL8 *chi1_AS, REAL8 *chi2_AS){
  if ((*m1 == *m2) && (*lambda1 != *lambda2))
    XLALPrintWarning("m1 == m2 (%g), but lambda1 != lambda2 (%g, %g).\n", *m1, *lambda1, *lambda2);

  REAL8 lambda1_tmp, lambda2_tmp, m1_tmp, m2_tmp, chi1_tmp, chi2_tmp;
  if (*m1>=*m2) {
    lambda1_tmp = *lambda1;
    lambda2_tmp = *lambda2;
    m1_tmp   = *m1;
    m2_tmp   = *m2;
    chi1_tmp = *chi1_AS;
    chi2_tmp = *chi2_AS;
  } else { /* swap spins, dimensionless tidal deformabilities, and masses */
    lambda1_tmp = *lambda2;
    lambda2_tmp = *lambda1;
    m1_tmp   = *m2;
    m2_tmp   = *m1;
    chi1_tmp = *chi2_AS;
    chi2_tmp = *chi1_AS;
  }
  *m1 = m1_tmp;
  *m2 = m2_tmp;
  *lambda1 = lambda1_tmp;
  *lambda2 = lambda2_tmp;
  *chi1_AS = chi1_tmp;
  *chi2_AS = chi2_tmp;

  if (*m1 < *m2)
    XLAL_ERROR(XLAL_EDOM, "XLAL_ERROR in EnforcePrimaryMassIsm1. When trying\
 to enfore that m1 should be the larger mass.\
 After trying to enforce this m1 = %f and m2 = %f\n", *m1, *m2);

  return XLAL_SUCCESS;
}


/**
 * convenient function to compute tidal coupling constant. Eq. 2 in arXiv:1706.02969
 * given masses and lambda numbers
 */
double XLALSimNRTunedTidesComputeKappa2T(
					 REAL8 m1_SI, /**< Mass of companion 1 (kg) */
					 REAL8 m2_SI, /**< Mass of companion 2 (kg) */
					 REAL8 lambda1, /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
					 REAL8 lambda2 /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
					 )
{
  int errcode = EnforcePrimaryMassIsm1(&m1_SI, &m2_SI, &lambda1, &lambda2);
  XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "EnforcePrimaryMassIsm1 failed");

  const REAL8 m1 = m1_SI / LAL_MSUN_SI;
  const REAL8 m2 = m2_SI / LAL_MSUN_SI;
  const REAL8 mtot = m1 + m2;

  /* Xa and Xb are the masses normalised for a total mass = 1 */
  /* not the masses appear symmetrically so we don't need to switch them. */
  const REAL8 Xa = m1 / mtot;
  const REAL8 Xb = m2 / mtot;

  /* tidal coupling constant. This is the
    kappa^T_eff = 2/13 [  (1 + 12 X_B/X_A) (X_A/C_A)^5 k^A_2 +  [A <- -> B]  ]
    from Tim Dietrich */

  /* Note that 2*k_2^A/c_A^5 = 3*lambda1 */
  const REAL8 term1 = ( 1.0 + 12.0*Xb/Xa ) * pow(Xa, 5.0) * lambda1;
  const REAL8 term2 = ( 1.0 + 12.0*Xa/Xb ) * pow(Xb, 5.0) * lambda2;
  const REAL8 kappa2T = (3.0/13.0) * ( term1 + term2 );

  return kappa2T;
}

/**
 * compute the merger frequency of a BNS system.
 * Tim's new fit that incorporates mass-ratio and asymptotes to zero for large kappa2T.
 */
double XLALSimNRTunedTidesMergerFrequency(
					  const REAL8 mtot_MSUN, /**< total mass of system (solar masses) */
					  const REAL8 kappa2T,   /**< tidal coupling constant. Eq. 2 in arXiv:1706.02969 */
					  const REAL8 q          /**< mass-ratio q >= 1 */
					  )
{
  if (q < 1.0) {
    XLALPrintError("XLAL Error - %s: q (%f) should be greater or equal to unity!\n",
		   __func__, q);
    XLAL_ERROR( XLAL_EDOM );
  }

  const REAL8 a_0 = 0.3586;
  const REAL8 n_1 = 3.35411203e-2;
  const REAL8 n_2 = 4.31460284e-5;
  const REAL8 d_1 = 7.54224145e-2;
  const REAL8 d_2 = 2.23626859e-4;

  const REAL8 kappa2T2 = kappa2T * kappa2T;

  const REAL8 num = 1.0 + n_1*kappa2T + n_2*kappa2T2;
  const REAL8 den = 1.0 + d_1*kappa2T + d_2*kappa2T2;
  const REAL8 Q_0 = a_0 / sqrt(q);

  /* dimensionless angular frequency of merger */
  const REAL8 Momega_merger = Q_0 * (num / den);

  /* convert from angular frequency to frequency (divide by 2*pi)
   * and then convert from
   * dimensionless frequency to Hz (divide by mtot * LAL_MTSUN_SI)
   */
  const REAL8 fHz_merger = Momega_merger / (LAL_TWOPI) / (mtot_MSUN * LAL_MTSUN_SI);

  return fHz_merger;
}


/**
 * compute the merger frequency of a BNS system for NRTidalv3 (https://arxiv.org/pdf/2311.07456.pdf).
 * This uses a new fit from Gonzalez, et. al (2022); Eq. (23) of https://arxiv.org/abs/2210.16366.
 */
double XLALSimNRTunedTidesMergerFrequency_v3(
            const REAL8 mtot_MSUN,   /**< total mass of system (solar masses) */
					  REAL8 lambda1, /**<dimensionless tidal deformability of star 1 */
            REAL8 lambda2, /**<dimensionless tidal deformability of star 2 */
            const REAL8 q, /**<mass ratio q>=1 */
            REAL8 chi1_AS, /**<aligned spin component of star 1 */
            REAL8 chi2_AS /**<aligned spin component of star 2 */
					  )
{
  if (q < 1.0) {
    XLALPrintError("XLAL Error - %s: q (%f) should be greater or equal to unity!\n",
		   __func__, q);
    XLAL_ERROR( XLAL_EDOM );
  }

  REAL8 Xa = q/(1.0+q);
  REAL8 Xb = 1.0 - Xa;

  REAL8 m1 = Xa*mtot_MSUN;
  REAL8 m2 = Xb*mtot_MSUN;

  int errcode = EnforcePrimaryMassIsm1_v3(&m1, &m2, &lambda1, &lambda2, &chi1_AS, &chi2_AS);
  XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "EnforcePrimaryMassIsm1_v3 failed");

  REAL8 nu = Xa*Xb;

  REAL8 Xa_2 = Xa*Xa;
  REAL8 Xa_3 = Xa_2*Xa;
  REAL8 Xb_2 = Xb*Xb;
  REAL8 Xb_3 = Xb_2*Xb;
  
  REAL8 kappa2eff = 3.0*nu*(Xa_3*lambda1 + Xb_3*lambda2);
  
  const REAL8 a_0 = 0.22;
  
  const REAL8 a_1M = 0.80;
  const REAL8 a_1S = 0.25;
  const REAL8 b_1S = -1.99; /**< This was not written in 2311.07456, but can be found in Table 2 of 2210.16366.*/

  const REAL8 a_1T = 0.0485;
  const REAL8 a_2T = 0.00000586;
  const REAL8 a_3T = 0.10;
  const REAL8 a_4T = 0.000186;

  const REAL8 b_1T = 1.80;
  const REAL8 b_2T = 599.99;
  const REAL8 b_3T = 7.80;
  const REAL8 b_4T = 84.76;

  const REAL8 Xval = 1.0 - 4.0*nu;

  const REAL8 p_1S = a_1S*(1.0 + b_1S*Xval);
  const REAL8 p_1T = a_1T*(1.0 + b_1T*Xval);
  const REAL8 p_2T = a_2T*(1.0 + b_2T*Xval);
  const REAL8 p_3T = a_3T*(1.0 + b_3T*Xval);
  const REAL8 p_4T = a_4T*(1.0 + b_4T*Xval);

  const REAL8 kappa2eff2 = kappa2eff*kappa2eff;
  
  const REAL8 Sval = Xa_2*chi1_AS + Xb_2*chi2_AS;

  const REAL8 QM = 1.0 + a_1M*Xval;
  const REAL8 QS = 1.0 + p_1S*Sval;

  const REAL8 num = 1.0 + p_1T*kappa2eff + p_2T*kappa2eff2;
  const REAL8 den = 1.0 + p_3T*kappa2eff + p_4T*kappa2eff2;
  const REAL8 QT = num / den;

  /* dimensionless angular frequency of merger */

  const REAL8 Qfit = a_0*QM*QS*QT;

  const REAL8 Momega_merger = nu*Qfit*(LAL_TWOPI); 


  /* convert from angular frequency to frequency (divide by 2*pi)
   * and then convert from
   * dimensionless frequency to Hz (divide by mtot * LAL_MTSUN_SI)
   */
  const REAL8 fHz_merger = Momega_merger / (LAL_TWOPI) / (mtot_MSUN * LAL_MTSUN_SI);

  return fHz_merger;
}


/**
 * Internal function only.
 * Function to call the frequency domain tidal correction.
 * Equation (7) in arXiv:1706.02969
 */
static double SimNRTunedTidesFDTidalPhase(
					  const REAL8 fHz, /**< Gravitational wave frequency (Hz) */
					  const REAL8 Xa, /**< Mass of companion 1 divided by total mass */
					  const REAL8 Xb, /**< Mass of companion 2 divided by total mass */
					  const REAL8 mtot, /**< total mass (Msun) */
					  const REAL8 kappa2T /**< tidal coupling constant. Eq. 2 in arXiv:1706.02969 */
					  )
{
  /* NRTunedTidesFDTidalPhase is Eq 7 in arXiv:1706.02969
   * and is a function of x = angular_orb_freq^(2/3)
   */
  REAL8 M_omega = LAL_PI * fHz * (mtot * LAL_MTSUN_SI); //dimensionless angular GW frequency

  REAL8 PN_x = pow(M_omega, 2.0/3.0);
  REAL8 PN_x_2 = PN_x * PN_x;
  REAL8 PN_x_3over2 = pow(PN_x, 3.0/2.0);
  REAL8 PN_x_5over2 = pow(PN_x, 5.0/2.0);

  /* model parameters */
  const REAL8 c_Newt = 2.4375; /* 39.0 / 16.0 */

  const REAL8 n_1 = -17.428;
  const REAL8 n_3over2 = 31.867;
  const REAL8 n_2 = -26.414;
  const REAL8 n_5over2 = 62.362;

  const REAL8 d_1 = n_1 - 2.496; /* 3115.0/1248.0 */
  const REAL8 d_3over2 = 36.089;

  REAL8 tidal_phase = - kappa2T * c_Newt / (Xa * Xb) * PN_x_5over2;

  REAL8 num = 1.0 + (n_1 * PN_x) + (n_3over2 * PN_x_3over2) + (n_2 * PN_x_2) + (n_5over2 * PN_x_5over2) ;
  REAL8 den = 1.0 + (d_1 * PN_x) + (d_3over2 * PN_x_3over2) ;

  REAL8 ratio = num / den;

  tidal_phase *= ratio;

  return tidal_phase;
}

/** 
 * Tidal amplitude corrections; only available for NRTidalv2;
 * Eq 24 of arxiv:1905.06011
 */
static REAL8 SimNRTunedTidesFDTidalAmplitude(
					     const REAL8 fHz, /**< Gravitational wave frequency (Hz) */
					     const REAL8 mtot, /**< Total mass in solar masses */
					     const REAL8 kappa2T /**< tidal coupling constant. Eq. 2 in arXiv:1706.02969 */
					     )
{
  const REAL8 M_sec   = (mtot * LAL_MTSUN_SI);

  REAL8 prefac = 0.0;
  prefac = 9.0*kappa2T;

  REAL8 x = pow(LAL_PI*M_sec*fHz, 2.0/3.0);
  REAL8 ampT = 0.0;
  REAL8 poly = 1.0;
  const REAL8 n1   = 4.157407407407407;
  const REAL8 n289 = 2519.111111111111;
  const REAL8 d    = 13477.8073677; 

  poly = (1.0 + n1*x + n289*pow(x, 2.89))/(1+d*pow(x,4.));
  ampT = - prefac*pow(x,3.25)*poly;

  return ampT;
}

/**
 * Set the NRTidalv2 phase coefficients in an array for use here and in the IMRPhenomX*_NRTidalv2 implementation
 */
int XLALSimNRTunedTidesSetFDTidalPhase_v2_Coeffs(REAL8 *NRTidalv2_coeffs)
{
  NRTidalv2_coeffs[0] =   2.4375; // c_Newt
  NRTidalv2_coeffs[1] = -12.615214237993088; // n_1
  NRTidalv2_coeffs[2] =  19.0537346970349; // n_3over2
  NRTidalv2_coeffs[3] = -21.166863146081035; // n_2
  NRTidalv2_coeffs[4] =  90.55082156324926; // n_5over2
  NRTidalv2_coeffs[5] = -60.25357801943598; // n_3
  NRTidalv2_coeffs[6] = -15.111207827736678; // d_1
  NRTidalv2_coeffs[7] =  22.195327350624694; // d_3over2
  NRTidalv2_coeffs[8] =   8.064109635305156; // d_2

  return XLAL_SUCCESS;
}

/** 
 * NRTunedTidesFDTidalPhase is Eq 22 of https://arxiv.org/abs/1905.06011 
 * and is a function of x = angular_orb_freq^(2/3)
 */
static double SimNRTunedTidesFDTidalPhase_v2(
					     const REAL8 fHz, /**< Gravitational wave frequency (Hz) */
					     const REAL8 Xa, /**< Mass of companion 1 divided by total mass */
					     const REAL8 Xb, /**< Mass of companion 2 divided by total mass */
					     const REAL8 mtot, /**< total mass (Msun) */
					     const REAL8 kappa2T /**< tidal coupling constant. Eq. 2 in arXiv:1706.02969 */
					     )
{

  REAL8 M_omega = LAL_PI * fHz * (mtot * LAL_MTSUN_SI); //dimensionless angular GW frequency
  REAL8 PN_x = pow(M_omega, 2.0/3.0);
  REAL8 PN_x_2 = PN_x * PN_x;
  REAL8 PN_x_3 = PN_x * PN_x_2;
  REAL8 PN_x_3over2 = pow(PN_x, 3.0/2.0);
  REAL8 PN_x_5over2 = pow(PN_x, 5.0/2.0);
  /* model parameters */
  REAL8 NRTidalv2_coeffs[9];

  int errcode;
  errcode = XLALSimNRTunedTidesSetFDTidalPhase_v2_Coeffs(NRTidalv2_coeffs);
  XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "Setting NRTidalv2 coefficients failed.\n");

  const REAL8 c_Newt   = NRTidalv2_coeffs[0];
  const REAL8 n_1      = NRTidalv2_coeffs[1];
  const REAL8 n_3over2 = NRTidalv2_coeffs[2];
  const REAL8 n_2      = NRTidalv2_coeffs[3];
  const REAL8 n_5over2 = NRTidalv2_coeffs[4];
  const REAL8 n_3      = NRTidalv2_coeffs[5];
  const REAL8 d_1      = NRTidalv2_coeffs[6];
  const REAL8 d_3over2 = NRTidalv2_coeffs[7];
  const REAL8 d_2      = NRTidalv2_coeffs[8];
  REAL8 tidal_phase = - kappa2T * c_Newt / (Xa * Xb) * PN_x_5over2;
  REAL8 num = 1.0 + (n_1 * PN_x) + (n_3over2 * PN_x_3over2) + (n_2 * PN_x_2) + (n_5over2 * PN_x_5over2) + (n_3 * PN_x_3);
  REAL8 den = 1.0 + (d_1 * PN_x) + (d_3over2 * PN_x_3over2) + (d_2 * PN_x_2) ;
  REAL8 ratio = num / den;
  tidal_phase *= ratio;
  return tidal_phase;
}

/** 
 * Coefficients or the PN tidal phase correction, at 7.5PN, to connect with NRTidalv3 Phase post-merger,
 * see Eq. (45) of https://arxiv.org/pdf/2311.07456.pdf
 */
int XLALSimNRTunedTidesSetFDTidalPhase_PN_Coeffs(REAL8 *PN_coeffs,/**<PN coefficients*/
               const REAL8 Xa /**< Mass of companion 1 divided by total mass*/)
{
  /* Powers of Xa and Xb */
  REAL8 Xb = 1.0 - Xa; //secondary mass divided by total mass
  REAL8 Xa_2 = Xa*Xa;
  REAL8 Xa_3 = Xa_2*Xa;
  REAL8 Xa_4 = Xa_3*Xa;
  REAL8 Xa_5 = Xa_4*Xa;

  REAL8 Xb_2 = Xb*Xb;
  REAL8 Xb_3 = Xb_2*Xb;
  REAL8 Xb_4 = Xb_3*Xb;
  REAL8 Xb_5 = Xb_4*Xb;

  REAL8 den_a = 11.0*Xa-12.0;
  REAL8 den_b = 11.0*Xb-12.0;

  /* 7.5PN Coefficients */
  PN_coeffs[0] = -3.0*den_a/(16*Xa*Xb_2);//(3.0)*(12.0 -11.0*Xa)/(16*Xa*Xb_2);
  PN_coeffs[1] = (-1300.0*Xa_3 + 11430.0*Xa_2 + 4595.0*Xa -15895.0)/(672.0*den_a);//-5.0*(260.0*Xa_3 - 2286.0*Xa_2 - 919.0*Xa + 3179.0)/(672.0*(11.0*Xa-12.0));
  PN_coeffs[2] = -LAL_PI;
  PN_coeffs[3] = (22861440.0*Xa_5 - 102135600.0*Xa_4 + 791891100.0*Xa_3 + 874828080.0*Xa_2 + 216234195.0*Xa -1939869350.0)/(27433728.0*den_a);//(5.0*(4572288.0*Xa_5 - 20427120.0*Xa_4 + 158378220.0*Xa_3 +174965616.0*Xa_2 + 43246839.0*Xa -387973870.0))/(27433728.0*(11.0*Xa - 12.0));
  PN_coeffs[4] = -LAL_PI*(10520.0*Xa_3 -7598.0*Xa_2 +22415.0*Xa - 27719.0)/(672.0*den_a);

  PN_coeffs[5] = -3.0*den_b/(16*Xb*Xa_2);//(3.0)*(12.0 -11.0*Xb)/(16*Xb*Xa_2);
  PN_coeffs[6] = (-1300.0*Xb_3 + 11430.0*Xb_2 + 4595.0*Xb -15895.0)/(672.0*den_b);//-5.0*(260.0*Xb_3 - 2286.0*Xb_2 - 919.0*Xb + 3179.0)/(672.0*(11.0*Xb-12.0));
  PN_coeffs[7] = -LAL_PI;
  PN_coeffs[8] = (22861440.0*Xb_5 - 102135600.0*Xb_4 + 791891100.0*Xb_3 + 874828080.0*Xb_2 + 216234195.0*Xb -1939869350.0)/(27433728.0*den_b);//(5.0*(4572288.0*Xb_5 - 20427120.0*Xb_4 + 158378220.0*Xb_3 +174965616.0*Xb_2 + 43246839.0*Xb -387973870.0))/(27433728.0*(11.0*Xb - 12.0));
  PN_coeffs[9] = -LAL_PI*(10520.0*Xb_3 -7598.0*Xb_2 +22415.0*Xb - 27719.0)/(672.0*den_b);

  return XLAL_SUCCESS;
}

/**
 * Set the NRTidalv3 effective love number and phase coefficients in an array for use here and in the IMRPhenomX*_NRTidalv3 implementation
 */
int XLALSimNRTunedTidesSetFDTidalPhase_v3_Coeffs(REAL8 *NRTidalv3_coeffs, /**< output of NRTidalv3 parameters to be used in computing the tidal phase corrections*/
               const REAL8 Xa, /**< Mass of companion 1 divided by total mass*/
               const REAL8 mtot, /**< total mass (Msun) */
               const REAL8 lambda1, /**< dimensionless tidal deformability of companion 1*/
               const REAL8 lambda2, /**< dimensionless tidal deformability of companion 2*/
               const REAL8 PN_coeffs[10] /**< 7.5 PN coefficients to be used for constraints*/
               )
{ 
  REAL8 Xb = 1.0 - Xa; //secondary mass divided by total mass
  REAL8 q = Xa/Xb;

  //Coefficients for the effective enhancement factor:
  const REAL8 s10 =   1.273000423; //s10
  const REAL8 s11 =   3.64169971e-03; //s11
  const REAL8 s12 =   1.76144380e-03; //s12
  
  const REAL8 s20 =   2.78793291e+01; //s20
  const REAL8 s21 =   1.18175396e-02; //s21
  const REAL8 s22 =   -5.39996790e-03; //s22 
  
  const REAL8 s30 =   1.42449682e-01; //s30
  const REAL8 s31 =   -1.70505852e-05; //s31
  const REAL8 s32 =   3.38040594e-05; //s32 

  REAL8 m1_SI = Xa * mtot * LAL_MSUN_SI;
  REAL8 m2_SI = Xb * mtot * LAL_MSUN_SI;

  REAL8 kappa2T = XLALSimNRTunedTidesComputeKappa2T(m1_SI, m2_SI, lambda1, lambda2);
  
  NRTidalv3_coeffs[0] = s10 + s11*kappa2T + s12*q*kappa2T; // s1
  NRTidalv3_coeffs[1] = s20 + s21*kappa2T + s22*q*kappa2T; // s2
  NRTidalv3_coeffs[2] = s30 + s31*kappa2T + s32*q*kappa2T; // s3

  REAL8 s2s3 = NRTidalv3_coeffs[1]*NRTidalv3_coeffs[2];
  NRTidalv3_coeffs[3] = cosh(s2s3) + sinh(s2s3);

  NRTidalv3_coeffs[4] = 3.0*Xb*Xa*Xa*Xa*Xa*lambda1; //kappaA
  NRTidalv3_coeffs[5] = 3.0*Xa*Xb*Xb*Xb*Xb*lambda2; //kappaB

//Exponent parameters:
  const REAL8 alpha =   -8.08155404e-03; //alpha
  const REAL8 beta =  -1.13695919e+00; //beta

  REAL8 kappaA_alpha = pow(NRTidalv3_coeffs[4] + 1, alpha); //kappaA_alpha
  REAL8 kappaB_alpha = pow(NRTidalv3_coeffs[5] + 1, alpha); //kappaB_alpha

  REAL8 Xa_beta = pow(Xa, beta); //Xa_beta
  REAL8 Xb_beta = pow(Xb, beta); //Xb_beta

//Pade approximant coefficients:
  const REAL8 n_5over20 =  -9.40654388e+02; //n_5over20
  const REAL8 n_5over21 =  6.26517157e+02; //n_5over21
  const REAL8 n_5over22 =  5.53629706e+02; //n_5over22
  const REAL8 n_5over23 =  8.84823087e+01; //n_5over23
  
  const REAL8 n_30 =  4.05483848e+02; //n_30
  const REAL8 n_31 =  -4.25525054e+02; //n_31
  const REAL8 n_32 = -1.92004957e+02; //n_32
  const REAL8 n_33 =  -5.10967553e+01; //n_33
  
  const REAL8 d_10 =  3.80343306e+00; //d_10
  const REAL8 d_11 =  -2.52026996e+01; //d_11
  const REAL8 d_12 =  -3.08054443e+00; //d_12

  NRTidalv3_coeffs[6] = n_5over20 + n_5over21*Xa + n_5over22*kappaA_alpha + n_5over23*Xa_beta; //n_5over2A
  NRTidalv3_coeffs[7] = n_30 + n_31*Xa + n_32*kappaA_alpha + n_33*Xa_beta; //n_3A
  NRTidalv3_coeffs[8] = d_10 + d_11*Xa + d_12*Xa_beta; //d_1A
  
  NRTidalv3_coeffs[9] = n_5over20  + n_5over21*Xb  + n_5over22*kappaB_alpha + n_5over23*Xb_beta; //n_5over2B
  NRTidalv3_coeffs[10] = n_30 + n_31*Xb + n_32*kappaB_alpha + n_33*Xb_beta; //n_3B
  NRTidalv3_coeffs[11] = d_10 + d_11*Xb + d_12*Xb_beta; //d_1B

  /* 7.5PN Coefficients */
  REAL8 c_1A = PN_coeffs[1];
  REAL8 c_3over2A = PN_coeffs[2];
  REAL8 c_2A = PN_coeffs[3];
  REAL8 c_5over2A = PN_coeffs[4];

  REAL8 c_1B = PN_coeffs[6];
  REAL8 c_3over2B = PN_coeffs[7];
  REAL8 c_2B = PN_coeffs[8];
  REAL8 c_5over2B = PN_coeffs[9];

  /* Pade Coefficients constrained with PN */
  REAL8 inv_c1_A = 1.0/c_1A;
  NRTidalv3_coeffs[12] = c_1A + NRTidalv3_coeffs[8]; //n_1A
  NRTidalv3_coeffs[13] = ((c_1A*c_3over2A) - c_5over2A - (c_3over2A)*NRTidalv3_coeffs[8] +  NRTidalv3_coeffs[6])*inv_c1_A; //n_3over2A
  NRTidalv3_coeffs[14] = c_2A + c_1A*NRTidalv3_coeffs[8];// + d_2A; //n_2A
  NRTidalv3_coeffs[15] = -(c_5over2A + c_3over2A*NRTidalv3_coeffs[8] -  NRTidalv3_coeffs[6])*inv_c1_A; //d_3over2A

  REAL8 inv_c1_B = 1.0/c_1B;
  NRTidalv3_coeffs[16] = c_1B + NRTidalv3_coeffs[11];//n_1B
  NRTidalv3_coeffs[17] = ((c_1B*c_3over2B) - c_5over2B - (c_3over2B)*NRTidalv3_coeffs[11] + NRTidalv3_coeffs[9])*inv_c1_B; //n_3over2B
  NRTidalv3_coeffs[18] = c_2B + c_1B*NRTidalv3_coeffs[11];// + d_2B; //n_2B
  NRTidalv3_coeffs[19] = -(c_5over2B + c_3over2B*NRTidalv3_coeffs[11] - NRTidalv3_coeffs[9])*inv_c1_B; //d_3over2B

  return XLAL_SUCCESS;
}

/** 
 * Tidal phase correction for NRTidalv3, Eq. (27,30), from Abac, et. al. (2023) (https://arxiv.org/pdf/2311.07456.pdf)
 * and is a function of x = angular_orb_freq^(2./3.)
 */
static double SimNRTunedTidesFDTidalPhase_v3(
               const REAL8 fHz, /**< Gravitational wave frequency (Hz) */
               const REAL8 mtot, /**< total mass (Msun) */
               const REAL8 NRTidalv3_coeffs[20], /**< NRTidalv3 coefficients */
               const REAL8 PN_coeffs[10] /**< 7.5 PN coefficients to be used as constraints*/
               )
{
  REAL8 M_omega = LAL_PI * fHz * (mtot * LAL_MTSUN_SI); //dimensionless angular GW frequency

  REAL8 s1 = NRTidalv3_coeffs[0]; 
  REAL8 s2 = NRTidalv3_coeffs[1]; 

  REAL8 exps2s3 = NRTidalv3_coeffs[3];

  REAL8 s2Mom = -s2*M_omega*2.0;
  REAL8 exps2Mom = cosh(s2Mom) + sinh(s2Mom);

  /* expression for the effective love number enhancement factor, see Eq. (27) of https://arxiv.org/pdf/2311.07456.pdf.*/
  REAL8 dynk2bar = 1.0 + ((s1) - 1)*(1.0/(1.0 + exps2Mom*exps2s3)) - ((s1-1.0)/(1.0 + exps2s3)) - 2.0*M_omega*((s1) - 1)*s2*exps2s3/((1.0 + exps2s3)*(1.0 + exps2s3));

  REAL8 PN_x = M_omega / cbrt(M_omega);               // pow(M_omega, 2.0/3.0)
  REAL8 PN_x_2 = PN_x * PN_x;                         // pow(PN_x, 2)
  REAL8 PN_x_3 = PN_x * PN_x_2;
  REAL8 PN_x_3over2 = PN_x * sqrt(PN_x);              // pow(PN_x, 3.0/2.0)
  REAL8 PN_x_5over2 = PN_x_3over2 * PN_x;      // pow(PN_x, 5.0/2.0)

  REAL8 kappaA = NRTidalv3_coeffs[4];
  REAL8 kappaB = NRTidalv3_coeffs[5];

  REAL8 dynkappaA = kappaA*dynk2bar;
  REAL8 dynkappaB = kappaB*dynk2bar;

  /* Pade Coefficients */
  REAL8 n_5over2A = NRTidalv3_coeffs[6];
  REAL8 n_3A = NRTidalv3_coeffs[7];
  REAL8 d_1A = NRTidalv3_coeffs[8];
  
  REAL8 n_5over2B = NRTidalv3_coeffs[9];
  REAL8 n_3B = NRTidalv3_coeffs[10];
  REAL8 d_1B = NRTidalv3_coeffs[11];

  /* 7.5PN Coefficients */
  REAL8 c_NewtA = PN_coeffs[0];
  REAL8 c_NewtB = PN_coeffs[5];

  /* Pade Coefficients constrained with PN */
  REAL8 n_1A = NRTidalv3_coeffs[12];
  REAL8 n_3over2A = NRTidalv3_coeffs[13];
  REAL8 n_2A = NRTidalv3_coeffs[14];
  REAL8 d_3over2A = NRTidalv3_coeffs[15];

  REAL8 n_1B = NRTidalv3_coeffs[16];
  REAL8 n_3over2B = NRTidalv3_coeffs[17];
  REAL8 n_2B = NRTidalv3_coeffs[18];
  REAL8 d_3over2B = NRTidalv3_coeffs[19];

  REAL8 factorA = -c_NewtA*PN_x_5over2*dynkappaA;
  REAL8 factorB = -c_NewtB*PN_x_5over2*dynkappaB;

  /* Pade approximant, see Eq. (32) of https://arxiv.org/pdf/2311.07456.pdf. */
  REAL8 numA = 1.0 + (n_1A*PN_x) + (n_3over2A*PN_x_3over2) + (n_2A*PN_x_2) + (n_5over2A*PN_x_5over2) + (n_3A*PN_x_3);
  REAL8 denA = 1.0 + (d_1A*PN_x) + (d_3over2A*PN_x_3over2);// + (d_2A*PN_x_2);
  
  REAL8 numB = 1.0 + (n_1B*PN_x) + (n_3over2B*PN_x_3over2) + (n_2B*PN_x_2) + (n_5over2B*PN_x_5over2) + (n_3B*PN_x_3);
  REAL8 denB = 1.0 + (d_1B*PN_x) + (d_3over2B*PN_x_3over2);// + (d_2B*PN_x_2);

  REAL8 ratioA = numA/denA;
  REAL8 ratioB = numB/denB;
  
  REAL8 tidal_phaseA = factorA*ratioA;
  REAL8 tidal_phaseB = factorB*ratioB;
  
  REAL8 tidal_phase = tidal_phaseA + tidal_phaseB;

  return tidal_phase;
}

/** 
 * PN tidal phase correction, at 7.5PN, to connect with NRTidalv3 Phase post-merger,
 * see Eq. (22) and (45) of https://arxiv.org/pdf/2311.07456.pdf
 * and is a function of x = angular_orb_freq^(2./3.)
 */
static double SimNRTunedTidesFDTidalPhase_PN(
               const REAL8 fHz, /**< Gravitational wave frequency (Hz) */
               const REAL8 Xa, /**< Mass of companion 1 divided by total mass */
               const REAL8 mtot, /**< total mass (Msun) */
               const REAL8 lambda1, /**< dimensionless tidal deformability of companion 1*/
               const REAL8 lambda2, /**< dimensionless tidal deformability of companion 2*/
               const REAL8 PN_coeffs[10] /**< 7.5 PN coefficients*/
               )
{
  REAL8 M_omega = LAL_PI * fHz * (mtot * LAL_MTSUN_SI); //dimensionless angular GW frequency

  REAL8 Xb = 1.0 - Xa;

  REAL8 PN_x = M_omega / cbrt(M_omega);               // pow(M_omega, 2.0/3.0)
  REAL8 PN_x_2 = PN_x * PN_x;                         // pow(PN_x, 2)
  REAL8 PN_x_3over2 = PN_x * sqrt(PN_x);              // pow(PN_x, 3.0/2.0)
  REAL8 PN_x_5over2 = PN_x_3over2 * PN_x;      // pow(PN_x, 5.0/2.0)

  REAL8 kappaA = 3.0*Xb*Xa*Xa*Xa*Xa*lambda1;
  REAL8 kappaB = 3.0*Xa*Xb*Xb*Xb*Xb*lambda2;

  /* 7.5PN Coefficients */
  REAL8 c_NewtA = PN_coeffs[0];
  REAL8 c_1A = PN_coeffs[1];
  REAL8 c_3over2A = PN_coeffs[2];
  REAL8 c_2A = PN_coeffs[3];
  REAL8 c_5over2A = PN_coeffs[4];

  REAL8 c_NewtB = PN_coeffs[5];
  REAL8 c_1B = PN_coeffs[6];
  REAL8 c_3over2B = PN_coeffs[7];
  REAL8 c_2B = PN_coeffs[8];
  REAL8 c_5over2B = PN_coeffs[9];

  REAL8 factorA = -c_NewtA*PN_x_5over2*kappaA;
  REAL8 factorB = -c_NewtB*PN_x_5over2*kappaB;

  REAL8 tidal_phasePNA = factorA*(1.0 + (c_1A*PN_x) + (c_3over2A*PN_x_3over2) + (c_2A*PN_x_2) + (c_5over2A*PN_x_5over2));
  REAL8 tidal_phasePNB = factorB*(1.0 + (c_1B*PN_x) + (c_3over2B*PN_x_3over2) + (c_2B*PN_x_2) + (c_5over2B*PN_x_5over2));
    
  REAL8 tidal_phasePN = tidal_phasePNA + tidal_phasePNB;

  return tidal_phasePN;
}

/** Function to call amplitude tidal series only; 
 * done for convenience to use for PhenomD_NRTidalv2 and 
 * SEOBNRv4_ROM_NRTidalv2
 */

int XLALSimNRTunedTidesFDTidalAmplitudeFrequencySeries(
						       const REAL8Sequence *amp_tidal, /**< [out] tidal amplitude frequency series */
						       const REAL8Sequence *fHz, /**< list of input Gravitational wave Frequency [Hz or dimensionless] */
						       REAL8 m1, /**< Mass of companion 1 in solar masses */
						       REAL8 m2, /**< Mass of companion 2 in solar masses */
						       REAL8 lambda1, /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
						       REAL8 lambda2 /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
						       )
{
  REAL8 m1_SI = m1 * LAL_MSUN_SI;
  REAL8 m2_SI = m2 * LAL_MSUN_SI;
  REAL8 f_dim_to_Hz;
  int errcode = EnforcePrimaryMassIsm1(&m1_SI, &m2_SI, &lambda1, &lambda2);
  XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "EnforcePrimaryMassIsm1 failed");
  
  if( lambda1 < 0 || lambda2 < 0 )
  XLAL_ERROR(XLAL_EFUNC, "lambda1 = %f, lambda2 = %f. Both should be greater than zero for NRTidal models", lambda1, lambda2);


  const REAL8 mtot = m1 + m2;
  /* SEOBNRv4ROM_NRTidalv2 and IMRPhenomD_NRTidalv2 deal with dimensionless freqs and freq in Hz;
   * If the value corresponding to the last index is above 1, we are safe to assume a frequency given in Hz,
   * otherwise a dimensionless frequency
   */

  if ((*fHz).data[(*fHz).length - 1] > 1.)
    f_dim_to_Hz = 1.;
  else 
    f_dim_to_Hz = mtot*LAL_MTSUN_SI;

  /* tidal coupling constant.*/
  const REAL8 kappa2T = XLALSimNRTunedTidesComputeKappa2T(m1_SI, m2_SI, lambda1, lambda2);

  for(UINT4 i = 0; i < (*fHz).length; i++)
    (*amp_tidal).data[i] = SimNRTunedTidesFDTidalAmplitude((*fHz).data[i]/f_dim_to_Hz, mtot, kappa2T);

  return XLAL_SUCCESS;
}

/**
 * Function to call the frequency domain tidal correction
 * over an array of input frequencies. This is
 * Equation (7) in arXiv:1706.02969 when NRTidal_version is NRTidal_V, 
 * or Equations (17)-(21) (for phasing) and Equation (24) (for amplitude) 
 * in arXiv:1905.06011 when NRTidal_version is NRTidalv2_V, 
 * or Equations (17)-(21) in arXiv:1905.06011 when NRTidal_version is NRTidalv2NoAmpCorr_V.
 * NoNRT_V specifies NO tidal phasing or amplitude is being added.
 * Note internally we use m1>=m2 - this is enforced in the code.
 * So any can be supplied
 *
 * The model for the tidal phase correction in NRTidal_V/NRTidalv2_V was calibrated
 * up to mass-ratio q=1.5 and kappa2T in [40, 5000].
 * The upper kappa2T limit is reached roughly for a
 * 1.4+1.4 BNS with lambda  = 2700 on both NSs.
 * In the high mass-ratio limit, the BNS merger frequency defined in
 * XLALSimNRTunedTidesMergerFrequency() asymptotes to zero. The waveform
 * amplitude should be tapered away starting at this frequency.
 * Therefore, no explicit limits are enforced.
 *
 * We also include here NRTidalv3_V, which was calibrated
 * up to mass ratio q = 2.0 and kappa2T in [40, 5000].
 * The NRTidalv3 tidal phase is connected to the 7.5PN tidal phase to minimize the presence of asymptotes post-merger.
 *
 */

int XLALSimNRTunedTidesFDTidalPhaseFrequencySeries(
						   const REAL8Sequence *phi_tidal, /**< [out] tidal phase frequency series */
						   const REAL8Sequence *amp_tidal, /**< [out] tidal amplitude frequency series */
						   const REAL8Sequence *planck_taper, /**< [out] planck tapering to be applied on overall signal */
						   const REAL8Sequence *fHz, /**< list of input Gravitational wave Frequency in Hz to evaluate */
						   REAL8 m1_SI, /**< Mass of companion 1 (kg) */
						   REAL8 m2_SI, /**< Mass of companion 2 (kg) */
						   REAL8 lambda1, /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
						   REAL8 lambda2, /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
               REAL8 chi1_AS, /**< aligned-spin component of mass 1*/
               REAL8 chi2_AS, /**< aligned-spin component of mass 2*/
						   NRTidal_version_type NRTidal_version /** < one of NRTidal_V, NRTidalv2_V, NRTidalv3_V, NRTidalv3NoAmpCorr_V, NRTidalv2NoAmpCorr_V or NoNRT_V */
						   )
{
  /* NOTE: internally m1 >= m2
   * This is enforced in the code below and we swap the lambda's
   * accordingly. For NRTidalv3, we also swap the aligned spin components chi1_AS and chi2_AS, on which the new merger frequency fit is dependent.
   */
  

  int errcode = EnforcePrimaryMassIsm1_v3(&m1_SI, &m2_SI, &lambda1, &lambda2, &chi1_AS, &chi2_AS);
  XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "EnforcePrimaryMassIsm1_v3 failed");
  
  if( lambda1 < 0 || lambda2 < 0 )
  XLAL_ERROR(XLAL_EFUNC, "lambda1 = %f, lambda2 = %f. Both should be greater than zero for NRTidal models", lambda1, lambda2);

  const REAL8 m1 = m1_SI / LAL_MSUN_SI;
  const REAL8 m2 = m2_SI / LAL_MSUN_SI;
  const REAL8 mtot = m1 + m2;
  const REAL8 q = m1 / m2;

  /* Xa and Xb are the masses normalised for a total mass = 1 */
  const REAL8 Xa = m1 / mtot;
  const REAL8 Xb = m2 / mtot;

  /* tidal coupling constant.*/
  const REAL8 kappa2T = XLALSimNRTunedTidesComputeKappa2T(m1_SI, m2_SI, lambda1, lambda2);

  /* Prepare tapering of amplitude beyond merger frequency */
  const REAL8 fHz_mrg = XLALSimNRTunedTidesMergerFrequency(mtot, kappa2T, q);
  const REAL8 fHz_mrg_v3 = XLALSimNRTunedTidesMergerFrequency_v3(mtot, lambda1, lambda2, q, chi1_AS, chi2_AS); //for NRTidalv3
  const REAL8 fHz_end_taper = 1.2*fHz_mrg;
  const REAL8 fHz_end_taper_v3 = 1.2*fHz_mrg_v3; //for NRTidalv3

  if (NRTidal_version == NRTidal_V) {
    for(UINT4 i = 0; i < (*fHz).length; i++){
      (*phi_tidal).data[i] = SimNRTunedTidesFDTidalPhase((*fHz).data[i], Xa, Xb, mtot, kappa2T);
      (*planck_taper).data[i] = 1.0 - PlanckTaper((*fHz).data[i], fHz_mrg, fHz_end_taper);
    }
  }
  else if (NRTidal_version == NRTidalv2_V) {
    for(UINT4 i = 0; i < (*fHz).length; i++) {
      (*phi_tidal).data[i] = SimNRTunedTidesFDTidalPhase_v2((*fHz).data[i], Xa, Xb, mtot, kappa2T);
      (*amp_tidal).data[i] = SimNRTunedTidesFDTidalAmplitude((*fHz).data[i], mtot, kappa2T);
      (*planck_taper).data[i] = 1.0 - PlanckTaper((*fHz).data[i], fHz_mrg, fHz_end_taper);
    }
  }
  else if (NRTidal_version == NRTidalv3_V) {
    int indexmin = -1; //initiate an invalid index where one first finds a minimum
    REAL8 PN_coeffs[10];
    XLALSimNRTunedTidesSetFDTidalPhase_PN_Coeffs(PN_coeffs, Xa);
    REAL8 NRTidalv3_coeffs[20];
    XLALSimNRTunedTidesSetFDTidalPhase_v3_Coeffs(NRTidalv3_coeffs, Xa, mtot, lambda1, lambda2, PN_coeffs);
    REAL8 fHzmrgcheck = 0.9 * fHz_mrg_v3; // start checking of minimum; 
    for (UINT4 i = 1; i < (*fHz).length; i++) {
        (*phi_tidal).data[i] = SimNRTunedTidesFDTidalPhase_v3((*fHz).data[i], mtot, NRTidalv3_coeffs, PN_coeffs);
        if ((*fHz).data[i] >= fHzmrgcheck && (*phi_tidal).data[i] >= (*phi_tidal).data[i-1]) {
            indexmin = i - 1;
            break;
        }
    }
    if (indexmin != -1) { //If a minimum is found, the tidal phase will be constant at that minimum value
        REAL8 tidal_min_value = (*phi_tidal).data[indexmin];
        for (UINT4 i = indexmin + 1; i < (*fHz).length; i++) {
            (*phi_tidal).data[i] = tidal_min_value;
        }
    }
    for(UINT4 i = 0; i < (*fHz).length; i++) {
      REAL8 planck_func = PlanckTaper((*fHz).data[i], 1.15*fHz_mrg_v3, 1.35*fHz_mrg_v3);
      /* We employ here the smooth connection between NRTidal and PN post-merger, Eq. (45) of https://arxiv.org/pdf/2311.07456.pdf*/
      (*phi_tidal).data[i] = (*phi_tidal).data[i]*(1.0 - planck_func) + SimNRTunedTidesFDTidalPhase_PN((*fHz).data[i], Xa, mtot, lambda1, lambda2, PN_coeffs)*planck_func;
      (*amp_tidal).data[i] = SimNRTunedTidesFDTidalAmplitude((*fHz).data[i], mtot, kappa2T);
      (*planck_taper).data[i] = 1.0 - PlanckTaper((*fHz).data[i], fHz_mrg_v3, fHz_end_taper_v3);
    }
  }
  else if (NRTidal_version == NRTidalv2NSBH_V) {
    for(UINT4 i = 0; i < (*fHz).length; i++) {
      (*phi_tidal).data[i] = SimNRTunedTidesFDTidalPhase_v2((*fHz).data[i], Xa, Xb, mtot, kappa2T);
      (*planck_taper).data[i] = 1.0;
    }
  }
  else if (NRTidal_version == NRTidalv2NoAmpCorr_V) {
    for(UINT4 i = 0; i < (*fHz).length; i++) {
      (*phi_tidal).data[i] = SimNRTunedTidesFDTidalPhase_v2((*fHz).data[i], Xa, Xb, mtot, kappa2T);
      (*planck_taper).data[i] = 1.0 - PlanckTaper((*fHz).data[i], fHz_mrg, fHz_end_taper);
    }
  }
  else if (NRTidal_version == NRTidalv3NoAmpCorr_V) {
    int indexmin = -1; //initiate an invalid index where one first finds a minimum
    REAL8 PN_coeffs[10];
    XLALSimNRTunedTidesSetFDTidalPhase_PN_Coeffs(PN_coeffs, Xa);
    REAL8 NRTidalv3_coeffs[20];
    XLALSimNRTunedTidesSetFDTidalPhase_v3_Coeffs(NRTidalv3_coeffs, Xa, mtot, lambda1, lambda2, PN_coeffs);
    REAL8 fHzmrgcheck = 0.9 * fHz_mrg_v3; // start checking of minimum; if a minimum is found, the tidal phase will be constant at that minimum value
    for (UINT4 i = 1; i < (*fHz).length; i++) {
        (*phi_tidal).data[i] = SimNRTunedTidesFDTidalPhase_v3((*fHz).data[i], mtot, NRTidalv3_coeffs, PN_coeffs);
        if ((*fHz).data[i] >= fHzmrgcheck && (*phi_tidal).data[i] >= (*phi_tidal).data[i-1]) {
            indexmin = i - 1;
            break;
        }
    }
    if (indexmin != -1) { //If a minimum is found, the tidal phase will be constant at that minimum value
        REAL8 tidal_min_value = (*phi_tidal).data[indexmin];
        for (UINT4 i = indexmin + 1; i < (*fHz).length; i++) {
            (*phi_tidal).data[i] = tidal_min_value;
        }
    }
    for(UINT4 i = 0; i < (*fHz).length; i++) {
      REAL8 planck_func = PlanckTaper((*fHz).data[i], 1.15*fHz_mrg_v3, 1.35*fHz_mrg_v3);
      /* We employ here the smooth connection between NRTidal and PN post-merger, Eq. (45) of https://arxiv.org/pdf/2311.07456.pdf*/
      (*phi_tidal).data[i] = (*phi_tidal).data[i]*(1.0 - planck_func) + SimNRTunedTidesFDTidalPhase_PN((*fHz).data[i], Xa, mtot, lambda1, lambda2, PN_coeffs)*planck_func;
      (*planck_taper).data[i] = 1.0 - PlanckTaper((*fHz).data[i], fHz_mrg_v3, fHz_end_taper_v3);
    }
  }
  else if (NRTidal_version == NoNRT_V)
    XLAL_ERROR( XLAL_EINVAL, "Trying to add NRTides to a BBH waveform!" );
  else
    XLAL_ERROR( XLAL_EINVAL, "Unknown version of NRTidal being used! At present, NRTidal_V, NRTidalv2_V, NRTidalv2NSBH_V, NRTidalv2NoAmpCorr_V and NoNRT_V are the only known ones!" );
  
  return XLAL_SUCCESS;
}

/**
 * Function to add 3.5PN spin-squared and 3.5PN spin-cubed terms. 
 * The spin-squared terms occur with the spin-induced quadrupole moment terms
 * while the spin-cubed terms occur with both spin-induced quadrupole as well as 
 * octupole moments. The terms are computed in arXiv:1806.01772 and are 
 * explicitly written out in Eq 27 of arXiv:1905.06011. The following terms
 * are specifically meant for BNS systems, and are added to the NRTidalv2
 * extensions of the approximants IMRPhenomPv2, IMRPhenomD and SEOBNRv4_ROM. 
 */

void XLALSimInspiralGetHOSpinTerms(
				   REAL8 *SS_3p5PN, /**< 3.5PN spin-spin tail term containing spin-induced quadrupole moment */
				   REAL8 *SSS_3p5PN, /**< 3.5 PN spin cubed term containing spin-induced octupole moment */
				   REAL8 X_A, /**< Mass fraction m_1/M for first component of binary */
				   REAL8 X_B, /**< Mass fraction m_2/M for second component of binary */
				   REAL8 chi1, /**< Aligned component of spin vector of first component of binary */
				   REAL8 chi2, /**< Aligned component of spin vector of second component of binary */
				   REAL8 quadparam1, /**< Spin-induced quadrupole moment parameter for component 1 */
				   REAL8 quadparam2 /**< Spin-induced quadrupole moment parameter for component 2 */
				   )
{
  REAL8 chi1_sq = 0., chi2_sq = 0.;
  REAL8 X_Asq = 0., X_Bsq = 0.;
  REAL8 octparam1 = 0, octparam2 = 0.;

  X_Asq = X_A*X_A;
  X_Bsq = X_B*X_B;

  chi1_sq = chi1*chi1;
  chi2_sq = chi2*chi2;

  /* Remove -1 to account for BBH baseline*/
  octparam1 = XLALSimUniversalRelationSpinInducedOctupoleVSSpinInducedQuadrupole(quadparam1)-1.;
  octparam2 = XLALSimUniversalRelationSpinInducedOctupoleVSSpinInducedQuadrupole(quadparam2)-1.;

  *SS_3p5PN = - 400.*LAL_PI*(quadparam1-1.)*chi1_sq*X_Asq - 400.*LAL_PI*(quadparam2-1.)*chi2_sq*X_Bsq;
  *SSS_3p5PN = 10.*((X_Asq+308./3.*X_A)*chi1+(X_Bsq-89./3.*X_B)*chi2)*(quadparam1-1.)*X_Asq*chi1_sq
    + 10.*((X_Bsq+308./3.*X_B)*chi2+(X_Asq-89./3.*X_A)*chi1)*(quadparam2-1.)*X_Bsq*chi2_sq
    - 440.*octparam1*X_A*X_Asq*chi1_sq*chi1 - 440.*octparam2*X_B*X_Bsq*chi2_sq*chi2;
}
