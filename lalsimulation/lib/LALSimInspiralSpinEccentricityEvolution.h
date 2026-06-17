/*
 * Copyright (C) 2024 Khun Sang Phukon, Nathan Johnson-McDaniel, Amitesh Singh, Anuradha Gupta
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


#ifndef _LALSIMINSPIRALSPINECCENTRICITYEVOLUTION_H
#define _LALSIMINSPIRALSPINECCENTRICITYEVOLUTION_H

#define LAL_NUM_SPIN_ECC_VARIABLES 17 /*number of equations or variables used for prcessing eccentric binaries */

#define LAL_GAMMA_1_3_MULT_GAMMA_2_3 3.6275987284684357011881565152843114645681324961854811511397698707 /*Multiplication of Gamma(1/3) and Gamma (2/3) */


#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif



#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include <lal/LALConstants.h>

/* Structure for static intrinsic parameters */
typedef struct
tagexpnCoeffsSpinEcc {
   /* mt: total mass in seconds, eta: symmetric mass ratio  */
   REAL8 mt,eta;

   // mt will be used in computation of GW frequency from orbital velocity and vice versa

   /* m: total mass, m1, m2: component masses (in solar masses), q: mass ratio (m2/m1) */
   REAL8 m, m1,m2, q;

   /* Squared masses, used in computing dimensionful spins (in total mass = 1 units) */
   REAL8 norm1, norm2;

   /* chi1,chi2: dimensionless spins  */
   REAL8 chi1, chi2;

   /* Cosines of twice the angle between projections of spins or spin combinations into
    * the orbital plane and the periastron line
    * cos2psi1 and cos2psi2 are for the individual spins
    * cos2psi_ReducedTotalSpin is for 'm1*\vec{chi1} * m2*\vec{chi2}'
    * */
   REAL8 cos2psi1, cos2psi2, cos2psi_ReducedTotalSpin;

   /* start and end frequencies of evolution */
   REAL8 fstart;
   REAL8 fend;

   /* variable to store values dOmega/dt, which is to be used in
    * computing the double derivative of Omega in one of the stopping tests
    *  */
   REAL8 prev_domega;

   /*enhancementFunc:  eccentricity enhancement function to be used from the following
    *                  1. EF_Arun_etal,                corresponding value 0
    *                  2. EF_LoutrelYunes_SuperAsym,   corresponding value 1
    *                  3. EF_LoutrelYunes_HyperAsym    corresponding value 2
    * */
   EnhancementFunction enhancementFunc;

   /* EccOrder: Desired order of eccentricity in the polynomials of enhancement functions
    * Eccentricity order values are: -1 (all orders), 20, 18, 16, 14,
    *                                12, 10, 8, 6, 4, 2
    * */
   LALSimInspiralSpinEccPolyEccOrder EccOrder;
   /* SpinOrder: Desired PN order of spin contributions in evolution equations
    * Spin order values with corresponding strings are:
    *                  1. LAL_SIM_INSPIRAL_SPIN_ORDER_ALL  -1
    *                  2. LAL_SIM_INSPIRAL_SPIN_ORDER_3PN   6
    *                  3. LAL_SIM_INSPIRAL_SPIN_ORDER_25PN  5
    *                  4. LAL_SIM_INSPIRAL_SPIN_ORDER_2PN   4
    *                  5. LAL_SIM_INSPIRAL_SPIN_ORDER_15PN  3
    *                  6. LAL_SIM_INSPIRAL_SPIN_ORDER_1PN   2
    *                  7. LAL_SIM_INSPIRAL_SPIN_ORDER_05PN  1
    *                  8. LAL_SIM_INSPIRAL_SPIN_ORDER_0PN   0
    *In this routine, maximum spin order is LAL_SIM_INSPIRAL_SPIN_ORDER_2PN
    * */
   LALSimInspiralSpinOrder SpinOrder;
   /* PhaseOrder: 2 times of desired PN order in evolution equation
    * Phase order values are: -1 (all orders), 6, 5, 4, 3, 2, 1, 0
    * */
   INT4 PhaseOrder;

   /* Flag for including the 2PN spin terms in the eccentricity evolution */
   EccEvol2PNSpinFlag et2PNSpinFlag;

   /* Newtonian orbital angular momentum */
   REAL8 LNmag;
   /* Semi-major and Semi-minor axes */
   REAL8 SemiMajor, SemiMinor;
   /* Dot and cross  products between spin and orbital angular momentum unit vectors */
   REAL8 s1dotl, s2dotl;
   REAL8 lcrosss1[3], lcrosss2[3], s1crosss2[3];
   /*Cross products between spin, orbital angular momentum and periastron line */
   REAL8 s1crossPhat[3], s2crossPhat[3], lcrossPhat[3];

   /* gamma1: Output of XLALSimInspiralSpinEccGamma1 */
   REAL8 gamma1;
   /* Fe: XLALSimInspiralSpinEccEnhancementFlux */
   REAL8 Fe;
   /* beta43: XLALSimInspiralSpinEccBeta(4,3, ...) */
   REAL8 beta43;

}expnCoeffsSpinEcc;

typedef struct tagXLALSimInspiralSpinEccPNEvolveOrbitParams
{
    REAL8 (*funcomega)(REAL8 omega, REAL8 etSq, REAL8 l, REAL8 lambda,  REAL8 s1x, REAL8 s1y, REAL8 s1z,  REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lx, REAL8 ly, REAL8 lz, REAL8 phatx, REAL8 phaty, REAL8 phatz, REAL8 k,  expnCoeffsSpinEcc *ak);
    REAL8 (*funcetSq)(REAL8 omega, REAL8 etSq, REAL8 l, REAL8 lambda,  REAL8 s1x, REAL8 s1y, REAL8 s1z,  REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lx, REAL8 ly, REAL8 lz, REAL8 phatx, REAL8 phaty, REAL8 phatz, REAL8 k, expnCoeffsSpinEcc *ak);
    REAL8 (*funcl)(REAL8 omega, REAL8 etSq, REAL8 l, REAL8 lambda,  REAL8 s1x, REAL8 s1y, REAL8 s1z,  REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lx, REAL8 ly, REAL8 lz, REAL8 phatx, REAL8 phaty, REAL8 phatz, REAL8 k, expnCoeffsSpinEcc *ak);
    REAL8 (*funclamb)(REAL8 omega, REAL8 etSq, REAL8 l, REAL8 lambda,  REAL8 s1x, REAL8 s1y, REAL8 s1z,  REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lx, REAL8 ly, REAL8 lz, REAL8 phatx, REAL8 phaty, REAL8 phatz, REAL8 k, expnCoeffsSpinEcc *ak);
    REAL8 (*funcs1x)(REAL8 omega, REAL8 etSq, REAL8 l, REAL8 lambda,  REAL8 s1x, REAL8 s1y, REAL8 s1z,  REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lx, REAL8 ly, REAL8 lz, REAL8 phatx, REAL8 phaty, REAL8 phatz, REAL8 k, expnCoeffsSpinEcc *ak);
    REAL8 (*funcs1y)(REAL8 omega, REAL8 etSq, REAL8 l, REAL8 lambda,  REAL8 s1x, REAL8 s1y, REAL8 s1z,  REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lx, REAL8 ly, REAL8 lz, REAL8 phatx, REAL8 phaty, REAL8 phatz, REAL8 k, expnCoeffsSpinEcc *ak);
    REAL8 (*funcs1z)(REAL8 omega, REAL8 etSq, REAL8 l, REAL8 lambda,  REAL8 s1x, REAL8 s1y, REAL8 s1z,  REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lx, REAL8 ly, REAL8 lz, REAL8 phatx, REAL8 phaty, REAL8 phatz, REAL8 k, expnCoeffsSpinEcc *ak);
    REAL8 (*funcs2x)(REAL8 omega, REAL8 etSq, REAL8 l, REAL8 lambda,  REAL8 s1x, REAL8 s1y, REAL8 s1z,  REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lx, REAL8 ly, REAL8 lz, REAL8 phatx, REAL8 phaty, REAL8 phatz, REAL8 k, expnCoeffsSpinEcc *ak);
    REAL8 (*funcs2y)(REAL8 omega, REAL8 etSq, REAL8 l, REAL8 lambda,  REAL8 s1x, REAL8 s1y, REAL8 s1z,  REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lx, REAL8 ly, REAL8 lz, REAL8 phatx, REAL8 phaty, REAL8 phatz, REAL8 k, expnCoeffsSpinEcc *ak);
    REAL8 (*funcs2z)(REAL8 omega, REAL8 etSq, REAL8 l, REAL8 lambda,  REAL8 s1x, REAL8 s1y, REAL8 s1z,  REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lx, REAL8 ly, REAL8 lz, REAL8 phatx, REAL8 phaty, REAL8 phatz, REAL8 k, expnCoeffsSpinEcc *ak);
    REAL8 (*funclx)(REAL8 omega, REAL8 etSq, REAL8 l, REAL8 lambda,  REAL8 s1x, REAL8 s1y, REAL8 s1z,  REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lx, REAL8 ly, REAL8 lz, REAL8 phatx, REAL8 phaty, REAL8 phatz, REAL8 k, expnCoeffsSpinEcc *ak);
    REAL8 (*funcly)(REAL8 omega, REAL8 etSq, REAL8 l, REAL8 lambda,  REAL8 s1x, REAL8 s1y, REAL8 s1z,  REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lx, REAL8 ly, REAL8 lz, REAL8 phatx, REAL8 phaty, REAL8 phatz, REAL8 k, expnCoeffsSpinEcc *ak);
    REAL8 (*funclz)(REAL8 omega, REAL8 etSq, REAL8 l, REAL8 lambda,  REAL8 s1x, REAL8 s1y, REAL8 s1z,  REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lx, REAL8 ly, REAL8 lz, REAL8 phatx, REAL8 phaty, REAL8 phatz, REAL8 k, expnCoeffsSpinEcc *ak);
     REAL8 (*funcpx)(REAL8 omega, REAL8 etSq, REAL8 l, REAL8 lambda,  REAL8 s1x, REAL8 s1y, REAL8 s1z,  REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lx, REAL8 ly, REAL8 lz, REAL8 phatx, REAL8 phaty, REAL8 phatz, REAL8 k, expnCoeffsSpinEcc *ak);
    REAL8 (*funcpy)(REAL8 omega, REAL8 etSq, REAL8 l, REAL8 lambda,  REAL8 s1x, REAL8 s1y, REAL8 s1z,  REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lx, REAL8 ly, REAL8 lz, REAL8 phatx, REAL8 phaty, REAL8 phatz, REAL8 k, expnCoeffsSpinEcc *ak);
    REAL8 (*funcpz)(REAL8 omega, REAL8 etSq, REAL8 l, REAL8 lambda,  REAL8 s1x, REAL8 s1y, REAL8 s1z,  REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lx, REAL8 ly, REAL8 lz, REAL8 phatx, REAL8 phaty, REAL8 phatz, REAL8 k, expnCoeffsSpinEcc *ak);
    REAL8 (*funck)(REAL8 omega, REAL8 etSq, REAL8 l, REAL8 lambda,  REAL8 s1x, REAL8 s1y, REAL8 s1z,  REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lx, REAL8 ly, REAL8 lz, REAL8 phatx, REAL8 phaty, REAL8 phatz, REAL8 k, expnCoeffsSpinEcc *ak);
    expnCoeffsSpinEcc  ak;

}XLALSimInspiralSpinEccPNEvolveOrbitParams;

typedef REAL8 (SimInspiralSpinEccEvolutionEquations)(
	REAL8 omega,		/* orbital frequency */
	REAL8 etSq, 		/* square of orbital eccentricity */
    REAL8 l,            /* mean anomaly */
    REAL8 lambda,       /* mean orbital phase */
    REAL8 s1x,          /* x-component of vector s1 */
    REAL8 s1y,          /* y-component of vector s1 */
    REAL8 s1z,          /* z-component of vector s1 */
    REAL8 s2x,          /* x-component of vector s2 */
    REAL8 s2y,          /* y-component of vector s2 */
    REAL8 s2z,          /* z-component of vector s2 */
    REAL8 lhatx,        /* x-component of unit vector of orbital angular momentum */
    REAL8 lhaty,        /* y-component of unit vector of orbital angular momentum */
    REAL8 lhatz,        /* z-component of unit vector of orbital angular momentum */
    REAL8 phatx,        /* x-component of unit vector of periastron line */
    REAL8 phaty,        /* y-component of unit vector of periastron line */
    REAL8 phatz,        /* z-component of unit vector of periastron line */
    REAL8 k,            /* periastron precession */
    expnCoeffsSpinEcc *ak
);

typedef struct
tagexpnFuncSpinEcc
{
    SimInspiralSpinEccEvolutionEquations *orbfreq;
    SimInspiralSpinEccEvolutionEquations *orbeccSq;
    SimInspiralSpinEccEvolutionEquations *meananom;
    SimInspiralSpinEccEvolutionEquations *meanorbitphase;
    SimInspiralSpinEccEvolutionEquations *compspin1x;
    SimInspiralSpinEccEvolutionEquations *compspin1y;
    SimInspiralSpinEccEvolutionEquations *compspin1z;
    SimInspiralSpinEccEvolutionEquations *compspin2x;
    SimInspiralSpinEccEvolutionEquations *compspin2y;
    SimInspiralSpinEccEvolutionEquations *compspin2z;
    SimInspiralSpinEccEvolutionEquations *orbangLx;
    SimInspiralSpinEccEvolutionEquations *orbangLy;
    SimInspiralSpinEccEvolutionEquations *orbangLz;
    SimInspiralSpinEccEvolutionEquations *periastronlinex;
    SimInspiralSpinEccEvolutionEquations *periastronliney;
    SimInspiralSpinEccEvolutionEquations *periastronlinez;
    SimInspiralSpinEccEvolutionEquations *periastronadvance;
} expnFuncSpinEcc;

// Function declaration
int XLALAdaptiveRungeKutta4HermiteOnlyFinalCheckStopping(LALAdaptiveRungeKuttaIntegrator *integrator,
                                                         void *params,
                                                         REAL8 *yinit,
                                                         REAL8 tinit,
                                                         REAL8 tend_in,
                                                         REAL8 deltat
                                                         );
int XLALSimInspiralSpinEccODESetup(
    expnCoeffsSpinEcc *ak,                          /* coefficients for TaylorT4 evolution [modified] */
    expnFuncSpinEcc *f,                             /* functions for TaylorT4 evolution [modified] */
    REAL8 m1_kg,                                    /* mass of companion 1 in kg */
    REAL8 m2_kg,                                    /* mass of companion 2 in kg*/
    REAL8 chi1,                                     /* dimensionless spin 1 */
    REAL8 chi2,                                     /* dimensionless spin 2 */
    REAL8 fstart,                                   /* start frequency of evolution   */
    REAL8 fend,                                     /* end frequency of evolution */
    EnhancementFunction enhancementFunc,            /* Enhancement functions */
    LALSimInspiralSpinEccPolyEccOrder EccOrder,     /* Order of Eccentricity in the Enhancement functions*/
    LALSimInspiralSpinOrder SpinOrder,              /* Twice PN order of spin effects */
    INT4 PhaseOrder,                                /* Twice PN phase order */
    EccEvol2PNSpinFlag et2PNSpinFlag                /* Whether to include the 2PN spin terms in the eccentricity evolution */
);
REAL8 XLALSimInspiralSpinEccLmagN( REAL8 omega, REAL8 etSq, expnCoeffsSpinEcc *ak  );
REAL8 XLALSimInspiralSpinEccBeta(  REAL8 a, REAL8 b,
        REAL8 s1x, REAL8 s1y, REAL8 s1z,
        REAL8 s2x, REAL8 s2y, REAL8 s2z,
        REAL8 lhatx, REAL8 lhaty, REAL8 lhatz,
        expnCoeffsSpinEcc *ak );
REAL8 XLALSimInspiralSpinEccGamma1(
        REAL8 s1x, REAL8 s1y, REAL8 s1z,
        REAL8 s2x, REAL8 s2y, REAL8 s2z,
        REAL8 lhatx, REAL8 lhaty, REAL8 lhatz,
        expnCoeffsSpinEcc *ak );
REAL8 XLALSimInspiralSpinEccPhi( REAL8 etSq );
REAL8 XLALSimInspiralSpinEccTildePhi( REAL8 etSq );
REAL8 XLALSimInspiralSpinEccPhiE_mult_etpow2( REAL8 etSq, EnhancementFunction enhancementFunc, LALSimInspiralSpinEccPolyEccOrder EccOrder );
REAL8 XLALSimInspiralSpinEccPsi( REAL8 etSq );
REAL8 XLALSimInspiralSpinEccTildePsi( REAL8 etSq );
REAL8 XLALSimInspiralSpinEccPsiE_mult_etpow2( REAL8 etSq, EnhancementFunction enhancementFunc, LALSimInspiralSpinEccPolyEccOrder EccOrder  );
REAL8 XLALSimInspiralSpinEccZeta( REAL8 etSq );
REAL8 XLALSimInspiralSpinEccTildeZeta( REAL8 etSq );
REAL8 XLALSimInspiralSpinEccTildeZetaSuperAsyLY( REAL8 etSq);
REAL8 XLALSimInspiralSpinEccTildeZetaHyperAsyLY( REAL8 etSq, LALSimInspiralSpinEccPolyEccOrder EccOrder);
REAL8 XLALSimInspiralSpinEccZetaE_mult_etpow2( REAL8 etSq, EnhancementFunction enhancementFunc, LALSimInspiralSpinEccPolyEccOrder EccOrder );
REAL8 XLALSimInspiralSpinEccTildeThetaSuperAsyLY( REAL8 etSq );
REAL8 XLALSimInspiralSpinEccEnhancementFlux( REAL8 etSq );
REAL8 XLALSimInspiralSpinEccEnhancementFluxE( REAL8 etSq );
REAL8 XLALSimInspiralSpinEccEnhancementKappa( REAL8 etSq );
REAL8 XLALSimInspiralSpinEccTildeEnhancementFlux( REAL8 etSq );
REAL8 XLALSimInspiralSpinEccEnhancementTildeKappa( REAL8 etSq );
REAL8 XLALSimInspiralSpinEccEnhancementKappaE_mult_etpow2( REAL8 etSq, EnhancementFunction enhancementFunc, LALSimInspiralSpinEccPolyEccOrder EccOrder );
REAL8 XLALSimInspiralSpinEccBetaSuperAsyLY( REAL8 etSq );
REAL8 XLALSimInspiralSpinEccPhiSuperAsyLY( REAL8 etSq);
REAL8 XLALSimInspiralSpinEccPhiHyperAsyLY( REAL8 etSq, LALSimInspiralSpinEccPolyEccOrder EccOrder);
REAL8 XLALSimInspiralSpinEccZetaSuperAsyLY( REAL8 etSq );
REAL8 XLALSimInspiralSpinEccEnhancementChiSuperAsyLY( REAL8 etSq );
REAL8 XLALSimInspiralSpinEccEnhancementChiHyperAsyLY( REAL8 etSq, LALSimInspiralSpinEccPolyEccOrder EccOrder );
REAL8 XLALSimInspiralSpinEccThetaSuperAsyLY( REAL8 etSq  );
REAL8 XLALSimInspiralSpinEccGammaSuperAsyLY( REAL8 etSq );
REAL8 XLALSimInspiralSpinEccTildeGammaSuperAsyLY( REAL8 etSq );
REAL8 XLALSimInspiralSpinEccZetaHyperAsyLY( REAL8 etSq, LALSimInspiralSpinEccPolyEccOrder EccOrder );
REAL8 XLALSimInspiralSpinEccAlphaSuperAsyLY( REAL8 etSq );
REAL8 XLALSimInspiralSpinEccTildeAlphaSuperAsyLY( REAL8 etSq );
REAL8 XLALSimInspiralSpinEccTildePsiSuperAsyLY( REAL8 etSq );
REAL8 XLALSimInspiralSpinEccTildePsiHyperAsyLY( REAL8 etSq, LALSimInspiralSpinEccPolyEccOrder EccOrder );
REAL8 XLALSimInspiralSpinEccPadePolyA( REAL8 etSq );
REAL8 XLALSimInspiralSpinEccPadePolyTildeA( REAL8 etSq );
REAL8 XLALSimInspiralSpinEccPsiHyperAsyLY( REAL8 etSq, LALSimInspiralSpinEccPolyEccOrder EccOrder );
REAL8 XLALSimInspiralSpinEccEnhancementKappaSuperAsyLY( REAL8 etSq );
REAL8 XLALSimInspiralSpinEccEnhancementKappaHyperAsyLY( REAL8 etSq, LALSimInspiralSpinEccPolyEccOrder EccOrder );
REAL8 XLALSimInspiralSpinEccTildePhiSuperAsyLY( REAL8 etSq );
REAL8 XLALSimInspiralSpinEccTildePhiHyperAsyLY( REAL8 etSq, LALSimInspiralSpinEccPolyEccOrder EccOrder );
REAL8 XLALSimInspiralSpinEccPsiSuperAsyLY( REAL8 etSq );
REAL8 XLALSimInspiralSpinEccPhik ( REAL8 etSq, EnhancementFunction enhancementFunc, LALSimInspiralSpinEccPolyEccOrder EccOrder );
REAL8 XLALSimInspiralSpinEccSigmaK ( REAL8 a, REAL8 b, REAL8 c, REAL8 a1a2, REAL8 b1b2, REAL8 c1c2,
                                    REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z,
                                    REAL8 lhatx, REAL8 lhaty, REAL8 lhatz, expnCoeffsSpinEcc *ak);
REAL8 XLALSimInspiralSpinEccEnhancementTildeChiSuperAsyLY( REAL8 etSq );
REAL8 XLALSimInspiralSpinEccEnhancementTildeChiHyperAsyLY( REAL8 etSq, LALSimInspiralSpinEccPolyEccOrder EccOrder );
REAL8 XLALSimInspiralSpinEccEnhancementTildeKappaSuperAsyLY( REAL8 etSq );
REAL8 XLALSimInspiralSpinEccEnhancementTildeKappaHyperAsyLY( REAL8 etSq, LALSimInspiralSpinEccPolyEccOrder EccOrder );
REAL8 XLALSimInspiralSpinEccPsiOmega( REAL8 etSq, EnhancementFunction enhancementFunc, LALSimInspiralSpinEccPolyEccOrder EccOrder );
REAL8 XLALSimInspiralSpinEccZetaOmega( REAL8 etSq, EnhancementFunction enhancementFunc, LALSimInspiralSpinEccPolyEccOrder EccOrder );
REAL8 XLALSimInspiralSpinEccTildeBetaSuperAsyLY ( REAL8 etSq );
REAL8 XLALSimInspiralSpinEccSemiMajorAxis( REAL8 omega, expnCoeffsSpinEcc *ak);
REAL8 XLALSimInspiralOrbitalEccentricitySqEvolution( REAL8 omega, REAL8 etSq, REAL8 UNUSED l, REAL8 UNUSED  lambda,
                                                   REAL8 s1x, REAL8 s1y, REAL8 s1z,
                                                   REAL8 s2x, REAL8 s2y, REAL8 s2z,
                                                   REAL8 lhatx, REAL8 lhaty, REAL8 lhatz,
                                                   REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
                                                   REAL8 UNUSED k, expnCoeffsSpinEcc *ak);
REAL8 XLALSimInspiralMeanAnomalyEvolution( REAL8 omega, REAL8 UNUSED etSq, REAL8 UNUSED l, REAL8 UNUSED lambda,
                                           REAL8 UNUSED s1x, REAL8 UNUSED s1y, REAL8 UNUSED s1z,
                                           REAL8 UNUSED s2x, REAL8 UNUSED s2y, REAL8 UNUSED s2z,
                                           REAL8 UNUSED lhatx, REAL8 UNUSED lhaty, REAL8 UNUSED lhatz,
                                           REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
                                           REAL8 k, expnCoeffsSpinEcc UNUSED *ak);
REAL8 XLALSimInspiralMeanOrbitalPhaseEvolution( REAL8 omega, REAL8 etSq, REAL8 UNUSED l, REAL8 UNUSED lambda,
                                                REAL8 UNUSED s1x, REAL8 UNUSED s1y, REAL8 UNUSED s1z,
                                                REAL8 UNUSED s2x, REAL8 UNUSED s2y, REAL8 UNUSED s2z,
                                                REAL8 UNUSED lhatx, REAL8 UNUSED lhaty, REAL8 UNUSED lhatz,
                                                REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
                                                REAL8 UNUSED k, expnCoeffsSpinEcc UNUSED *ak);
REAL8 XLALSimInspiralPeriastronPrecessionEvolution( REAL8 omega, REAL8 etSq, REAL8 UNUSED l, REAL8 UNUSED lambda,
                                                REAL8 UNUSED s1x, REAL8 UNUSED s1y, REAL8 UNUSED s1z,
                                                REAL8 UNUSED s2x, REAL8 UNUSED s2y, REAL8 UNUSED s2z,
                                                REAL8 UNUSED lhatx, REAL8 UNUSED lhaty, REAL8 UNUSED lhatz,
                                                REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
                                                REAL8 UNUSED k, expnCoeffsSpinEcc *ak);
REAL8 XLALSimInspiralOmegaEvolution( REAL8 omega, REAL8 etSq, REAL8 UNUSED l, REAL8 UNUSED lambda,
                                                REAL8 s1x, REAL8 s1y, REAL8 s1z,
                                                REAL8 s2x, REAL8 s2y, REAL8 s2z,
                                                REAL8 lhatx, REAL8 lhaty, REAL8 lhatz,
                                                REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
                                                REAL8 UNUSED k, expnCoeffsSpinEcc *ak);
REAL8 XLALSimInspiralS1xEvolution( REAL8 UNUSED omega, REAL8 UNUSED etSq, REAL8 UNUSED l, REAL8 UNUSED lambda,
                                                REAL8 UNUSED s1x, REAL8 UNUSED s1y, REAL8 UNUSED s1z,
                                                REAL8 UNUSED s2x, REAL8 UNUSED s2y, REAL8 UNUSED s2z,
                                                REAL8 UNUSED lhatx, REAL8 UNUSED lhaty, REAL8 UNUSED lhatz,
                                                REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
                                                REAL8 UNUSED k, expnCoeffsSpinEcc *ak);
REAL8 XLALSimInspiralS1yEvolution( REAL8 UNUSED omega, REAL8 UNUSED etSq, REAL8 UNUSED l, REAL8 UNUSED lambda,
                                                REAL8 UNUSED s1x, REAL8 UNUSED s1y, REAL8 UNUSED s1z,
                                                REAL8 UNUSED s2x, REAL8 UNUSED s2y, REAL8 UNUSED s2z,
                                                REAL8 UNUSED lhatx, REAL8 UNUSED lhaty, REAL8 UNUSED lhatz,
                                                REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
                                                REAL8 UNUSED k, expnCoeffsSpinEcc *ak);
REAL8 XLALSimInspiralS1zEvolution( REAL8 UNUSED omega, REAL8 UNUSED etSq, REAL8 UNUSED l, REAL8 UNUSED lambda,
                                                REAL8 UNUSED s1x, REAL8 UNUSED s1y, REAL8 UNUSED s1z,
                                                REAL8 UNUSED s2x, REAL8 UNUSED s2y, REAL8 UNUSED s2z,
                                                REAL8 UNUSED lhatx, REAL8 UNUSED lhaty, REAL8 UNUSED lhatz,
                                                REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
                                                REAL8 UNUSED k, expnCoeffsSpinEcc *ak);
REAL8 XLALSimInspiralS2xEvolution( REAL8 UNUSED omega, REAL8 UNUSED etSq, REAL8 UNUSED l, REAL8 UNUSED lambda,
                                                REAL8 UNUSED s1x, REAL8 UNUSED s1y, REAL8 UNUSED s1z,
                                                REAL8 UNUSED s2x, REAL8 UNUSED s2y, REAL8 UNUSED s2z,
                                                REAL8 UNUSED lhatx, REAL8 UNUSED lhaty, REAL8 UNUSED lhatz,
                                                REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
                                                REAL8 UNUSED k, expnCoeffsSpinEcc *ak);
REAL8 XLALSimInspiralS2yEvolution( REAL8 UNUSED omega, REAL8 UNUSED etSq, REAL8 UNUSED l, REAL8 UNUSED lambda,
                                                REAL8 UNUSED s1x, REAL8 UNUSED s1y, REAL8 UNUSED s1z,
                                                REAL8 UNUSED s2x, REAL8 UNUSED s2y, REAL8 UNUSED s2z,
                                                REAL8 UNUSED lhatx, REAL8 UNUSED lhaty, REAL8 UNUSED lhatz,
                                                REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
                                                REAL8 UNUSED k, expnCoeffsSpinEcc *ak);
REAL8 XLALSimInspiralS2zEvolution( REAL8 UNUSED omega, REAL8 UNUSED etSq, REAL8 UNUSED l, REAL8 UNUSED lambda,
                                                REAL8 UNUSED s1x, REAL8 UNUSED s1y, REAL8 UNUSED s1z,
                                                REAL8 UNUSED s2x, REAL8 UNUSED s2y, REAL8 UNUSED s2z,
                                                REAL8 UNUSED lhatx, REAL8 UNUSED lhaty, REAL8 UNUSED lhatz,
                                                REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
                                                REAL8 UNUSED k, expnCoeffsSpinEcc *ak);
REAL8 XLALSimInspiralLhatxEvolution( REAL8 omega, REAL8 etSq, REAL8 UNUSED l, REAL8 UNUSED lambda,
                                                REAL8 UNUSED s1x, REAL8 UNUSED s1y, REAL8 UNUSED s1z,
                                                REAL8 UNUSED s2x, REAL8 UNUSED s2y, REAL8 UNUSED s2z,
                                                REAL8 UNUSED lhatx, REAL8 UNUSED lhaty, REAL8 UNUSED lhatz,
                                                REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
                                                REAL8 UNUSED k, expnCoeffsSpinEcc *ak);
REAL8 XLALSimInspiralLhatyEvolution( REAL8 omega, REAL8 etSq, REAL8 UNUSED l, REAL8 UNUSED lambda,
                                                REAL8 UNUSED s1x, REAL8 UNUSED s1y, REAL8 UNUSED s1z,
                                                REAL8 UNUSED s2x, REAL8 UNUSED s2y, REAL8 UNUSED s2z,
                                                REAL8 UNUSED lhatx, REAL8 UNUSED lhaty, REAL8 UNUSED lhatz,
                                                REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
                                                REAL8 UNUSED k, expnCoeffsSpinEcc *ak);
REAL8 XLALSimInspiralLhatzEvolution( REAL8 omega, REAL8 etSq, REAL8 UNUSED l, REAL8 UNUSED lambda,
                                                REAL8 UNUSED s1x, REAL8 UNUSED s1y, REAL8 UNUSED s1z,
                                                REAL8 UNUSED s2x, REAL8 UNUSED s2y, REAL8 UNUSED s2z,
                                                REAL8 UNUSED lhatx, REAL8 UNUSED lhaty, REAL8 UNUSED lhatz,
                                                REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
                                                REAL8 UNUSED k, expnCoeffsSpinEcc *ak);
REAL8 XLALSimInspiralPeriastronLinexEvolution( REAL8 omega, REAL8 UNUSED etSq, REAL8 UNUSED l, REAL8 UNUSED lambda,
                                                REAL8 UNUSED s1x, REAL8 UNUSED s1y, REAL8 UNUSED s1z,
                                                REAL8 UNUSED s2x, REAL8 UNUSED s2y, REAL8 UNUSED s2z,
                                                REAL8 UNUSED lhatx, REAL8 UNUSED lhaty, REAL8 UNUSED lhatz,
                                                REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
                                                REAL8 k, expnCoeffsSpinEcc *ak);
REAL8 XLALSimInspiralPeriastronLineyEvolution( REAL8 omega, REAL8 UNUSED etSq, REAL8 UNUSED l, REAL8 UNUSED lambda,
                                                REAL8 UNUSED s1x, REAL8 UNUSED s1y, REAL8 UNUSED  s1z,
                                                REAL8 UNUSED s2x, REAL8 UNUSED s2y, REAL8 UNUSED s2z,
                                                REAL8 UNUSED lhatx, REAL8 UNUSED lhaty, REAL8 UNUSED lhatz,
                                                REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
                                                REAL8 k, expnCoeffsSpinEcc *ak);
REAL8 XLALSimInspiralPeriastronLinezEvolution( REAL8 omega, REAL8 UNUSED etSq, REAL8 UNUSED l, REAL8 UNUSED lambda,
                                                REAL8 UNUSED s1x, REAL8 UNUSED s1y, REAL8 UNUSED s1z,
                                                REAL8 UNUSED s2x, REAL8 UNUSED s2y, REAL8 UNUSED s2z,
                                                REAL8 UNUSED lhatx, REAL8 UNUSED lhaty, REAL8 UNUSED lhatz,
                                                REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
                                                REAL8 k, expnCoeffsSpinEcc *ak);
REAL8 XLALSimInspiralPeriastronPrecession( REAL8 xbar, REAL8 etSq, expnCoeffsSpinEcc *ak,
					        LALSimInspiralSpinOrder spinO, INT4 phaseO);

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif


#endif /* _LALSIMINSPIRALSPINECCENTRICITYEVOLUTION_H */
