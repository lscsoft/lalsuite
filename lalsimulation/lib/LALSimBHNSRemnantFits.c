#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/LALSimIMR.h>

#include "LALSimBHNSRemnantFits.h"

/**
 * @author Andrew Matas, Jonathan Thompson, Edward Fauchon-Jones, Sebastian Khan
 * @addtogroup LALSimBHNSRemnantFits_c Module LALSimBHNSRemnantFits.c
 * @ingroup lalsimulation_general
 * @brief Provides routines for NSBH remnant fits.
 *
 * This module provides an implementation of the NSBH remnant fits described by
 * \cite Zappa:2019ntl and the underlying BBH remnant fits described by
 * \cite Jimenez-Forteza:2016oae. For a standalone python implementation please
 * see https://git.tpi.uni-jena.de/core/bhnsremnant.
 * @{
 *
 * @name NSBH remnant routines
 * @{
 */

/**
 * Compute final black hole mass for aligned black hole spin
 *
 * See Table I (p7) in https://arxiv.org/abs/1903.11622 for the definition of
 * the coefficients of this model.
 */
REAL8 XLALBHNS_mass_aligned(
    const REAL8 m1,     /**< Mass of BH (companion 1) (solar masses) */
    const REAL8 m2,     /**< Mass of NS (companion 2) (solar masses) */
    const REAL8 chi1,   /**< Dimensionless spin of BH (companion 1) */
    const REAL8 lam     /**< Dimensionless tidal deformability of NS (companion 2) */
){
    //TODO: Insert error checkking
    REAL8 mtot=m1+m2;
    REAL8 nu=m1*m2/pow(mtot,2);
    REAL8 massc[]={-1.83417425e-03,  2.39226041e-03,  4.29407902e-03,  9.79775571e-03,
              2.33868869e-07, -8.28090025e-07, -1.64315549e-06, 8.08340931e-06,
             -2.00726981e-02,  1.31986011e-01,  6.50754064e-02, -1.42749961e-01};
    REAL8 model=model3a(nu,chi1,lam,massc);
    if(chi1<0 && nu<0.188){
        model=1;
    }
    if(chi1<-0.5){
        model=1;
    }
    if(model>1){
        model=1;
    }
    return XLALbbh_final_mass_non_precessing_UIB2016(m1,m2,chi1,0.) * model;
}

/**
 * Compute final black hole spin for aligned black hole spin
 *
 * See Table I (p7) in https://arxiv.org/abs/1903.11622 for the definition of
 * the coefficients of this model.
 */
REAL8 XLALBHNS_spin_aligned(
    const REAL8 m1,    /**< Mass of BH (companion 1) (solar masses) */
    const REAL8 m2,    /**< Mass of NS (companion 2) (solar masses) */
    const REAL8 chi1,  /**< Dimensionless spin of BH (companion 1) */
    const REAL8 lam    /**< Dimensionless tidal deformability of NS (companion 2) */
){
    //TODO: Insert error checking
    REAL8 mtot=m1+m2;
    REAL8 nu=m1*m2/pow(mtot,2);
    REAL8 spinc[]={-5.44187381e-03,  7.91165608e-03,  2.33362046e-02,  2.47764497e-02,
                 -8.56844797e-07, -2.81727682e-06,  6.61290966e-06,  4.28979016e-05,
                 -3.04174272e-02,  2.54889050e-01,  1.47549350e-01, -4.27905832e-01};

    REAL8 model=model3a(nu,chi1,lam,spinc);
    if(chi1<0 && nu<0.188){
        model=1;
    }
    if(chi1<-0.5){
        model=1;
    }
    if(model>1){
        model=1;
    }
    return XLALbbh_final_spin_non_precessing_UIB2016(m1,m2,chi1,0.) * model;
}

/**
 * Template 3D model for an NSBH remnant property
 *
 * See Eq. (4) (p7) in https://arxiv.org/abs/1903.11622 for the definition of
 * this function.
 */
static REAL8 model3a(
    const REAL8 nu,    /**< Symmetric mass ratio */
    const REAL8 ai,    /**< Dimensionless spin of BH (companion 1) */
    const REAL8 lam,   /**< Dimensionless tidal deformability of NS (companion 2) */
    const REAL8 *par   /**< Best fit parameters for particular instance of remnant fit */
){
    REAL8 p[]={0,0,0};
    pijk_to_pk(p,nu,ai,par);
    return model1a(lam,p);
}

/**
 * Calculate parameters for 1D model
 *
 * See Eq. (5,6) (p7) in https://arxiv.org/abs/1903.11622 for the definition of
 * this function.
 */
static void pijk_to_pk(
    REAL8 *p,          /**< Output: Coefficients for 1D model */
    const REAL8 nu,    /**< Symmetric mass ratio */
    const REAL8 ai,    /**< Dimensionless spin of BH (companion 1) */
    const REAL8 *par   /**< Best fit parameters for particular instance of remnant fit */
){
    //REAL8 p[] = {0, 0, 0};
    REAL8 p110 = par[0];
    REAL8 p111 = par[1];
    REAL8 p120 = par[2];
    REAL8 p121 = par[3];
    REAL8 p210 = par[4];
    REAL8 p211 = par[5];
    REAL8 p220 = par[6];
    REAL8 p221 = par[7];
    REAL8 p310 = par[8];
    REAL8 p311 = par[9];
    REAL8 p320 = par[10];
    REAL8 p321 = par[11];
    REAL8 p11 = p110*ai + p111;
    REAL8 p12 = p120*ai + p121;
    REAL8 p21 = p210*ai + p211;
    REAL8 p22 = p220*ai + p221;
    REAL8 p31 = p310*ai + p311;
    REAL8 p32 = p320*ai + p321;
    p[0] = (p11 + p12*nu)*nu;
    p[1] = (p21 + p22*nu)*nu;
    p[2] = (p31 + p32*nu)*nu;
    //return p;
}

/**
 * Template 1D model for an NSBH remnant property
 *
 * See Eq. (4) (p7) in https://arxiv.org/abs/1903.11622 for the definition of
 * this function.
 */
static double model1a(
    const REAL8 x,   /**< Input parameter for 1D model */
    const REAL8 *p   /**< Coefficients for 1D model */
){
    return (1 + x * p[0] + pow(x,2) * p[1]) / pow(1 + x * pow(p[2],2),2);
}

/**
 * @}
 *
 * @name BBH remnant routines
 * @{
 */

/**
 * Calculate the final mass with the aligned-spin NR fit
 *
 * Calculate the final mass with the aligned-spin NR fit by Xisco Jimenez
 * Forteza, David Keitel, Sascha Husa et al. [LIGO-P1600270]
 * [https://arxiv.org/abs/1611.00332] versions v1 and v2 use the same ansatz,
 * with v2 calibrated to additional SXS and RIT data
 */
REAL8 XLALbbh_final_mass_non_precessing_UIB2016(
    const REAL8 m1,      /**< component mass of first BH (solar masses) */
    const REAL8 m2,      /**< component mass of second BH (solar masses) */
    const REAL8 chi1,    /**< dimensionless spin of first BH */
    const REAL8 chi2     /**< dimensionless spin of second BH */
){
    // binary masses
    REAL8 m    = m1+m2;
    REAL8 msq  = m*m;
    REAL8 m1sq = m1*m1;
    REAL8 m2sq = m2*m2;
    // symmetric mass ratio
    REAL8 eta  = m1*m2/msq;
    if (eta>0.25){
      printf("Truncating eta from above to 0.25. This should only be necessary in some rounding corner cases, but better check your m1 and m2 inputs...");
      eta = 0.25;
    }
    if (eta<0.0){
      printf("Truncating negative eta to 0.0. This should only be necessary in some rounding corner cases, but better check your m1 and m2 inputs...");
      eta = 0.0;
    }
    REAL8 eta2 = eta*eta;
    REAL8 eta3 = eta2*eta;
    REAL8 eta4 = eta2*eta2;
    // spin variables (in m = 1 units)
    //REAL8 S1    = chi1*m1sq/msq; // spin angular momentum 1
    //REAL8 S2    = chi2*m2sq/msq; // spin angular momentum 2
    //REAL8 Stot  = S1+S2;         // total spin
    REAL8 Shat  = (chi1*m1sq+chi2*m2sq)/(m1sq+m2sq); // effective spin, = msq*Stot/(m1sq+m2sq)
    REAL8 Shat2 = Shat*Shat;
    REAL8 Shat3 = Shat2*Shat;
    //REAL8 Shat4 = Shat2*Shat2;
    // spin difference, assuming m1>m2
    REAL8 chidiff  = chi1 - chi2;
    if (m2>m1){ // fit assumes m1>m2
      chidiff = -chidiff;
     }
    REAL8 chidiff2 = chidiff*chidiff;
    // typical squareroots and functions of eta
    REAL8 sqrt2 = pow(2.,0.5);
    //REAL8 sqrt3 = pow(3.,0.5);
    REAL8 sqrt1m4eta = pow((1. - 4.*eta),0.5);


    // rational-function Pade coefficients (exact) from Eq. (22) of 1611.00332v2
    REAL8 b10 = 0.346;
    REAL8 b20 = 0.211;
    REAL8 b30 = 0.128;
    REAL8 b50 = -0.212;
    // fit coefficients from Tables VII-X of 1611.00332v2
    // values at increased numerical precision copied from
    // https://git.ligo.org/uib-papers/finalstate2016/blob/master/LALInference/EradUIB2016v2_pyform_coeffs.txt
    // git commit f490774d3593adff5bb09ae26b7efc6deab76a42
    REAL8 a2 = 0.5609904135313374;
    REAL8 a3 = -0.84667563764404;
    REAL8 a4 = 3.145145224278187;
    REAL8 b1 = -0.2091189048177395;
    REAL8 b2 = -0.19709136361080587;
    REAL8 b3 = -0.1588185739358418;
    REAL8 b5 = 2.9852925538232014;
    REAL8 f20 = 4.271313308472851;
    REAL8 f30 = 31.08987570280556;
    REAL8 f50 = 1.5673498395263061;
    REAL8 f10 = 1.8083565298668276;
    REAL8 f21 = 0.;
    REAL8 d10 = -0.09803730445895877;
    REAL8 d11 = -3.2283713377939134;
    REAL8 d20 = 0.01118530335431078;
    REAL8 d30 = -0.01978238971523653;
    REAL8 d31 = -4.91667749015812;
    REAL8 f11 = 15.738082204419655;
    REAL8 f31 = -243.6299258830685;
    REAL8 f51 = -0.5808669012986468;
    // Calculate the radiated-energy fit from Eq. (27) of 1611.00332
    REAL8 Erad = (((1. + -2.0/3.0*sqrt2)*eta + a2*eta2 + a3*eta3 + a4*eta4)*(1. + b10*b1*Shat*(f10 + f11*eta + (16. - 16.*f10 - 4.*f11)*eta2) + b20*b2*Shat2*(f20 + f21*eta + (16. - 16.*f20 - 4.*f21)*eta2) + b30*b3*Shat3*(f30 + f31*eta + (16. - 16.*f30 - 4.*f31)*eta2)))/(1. + b50*b5*Shat*(f50 + f51*eta + (16. - 16.*f50 - 4.*f51)*eta2)) + d10*sqrt1m4eta*eta2*(1. + d11*eta)*chidiff + d30*Shat*sqrt1m4eta*eta*(1. + d31*eta)*chidiff + d20*eta3*chidiff2;
    // Convert to actual final mass
    REAL8 Mf = m*(1.-Erad);
    return Mf;
}

/**
 * Calculate the final spin with the aligned-spin NR fit
 *
 * Calculate the final spin with the aligned-spin NR fit by Xisco Jimenez
 * Forteza, David Keitel, Sascha Husa et al. [LIGO-P1600270]
 * [https://arxiv.org/abs/1611.00332] versions v1 and v2 use the same ansatz,
 * with v2 calibrated to additional SXS and RIT data
 */
REAL8 XLALbbh_final_spin_non_precessing_UIB2016(
    const REAL8 m1,      /**< component mass of first BH (solar masses) */
    const REAL8 m2,      /**< component mass of second BH (solar masses) */
    const REAL8 chi1,    /**< dimensionless spin of first BH */
    const REAL8 chi2     /**< dimensionless spin of second BH */
){
    // binary masses
    REAL8 m    = m1+m2;
    REAL8 msq  = m*m;
    REAL8 m1sq = m1*m1;
    REAL8 m2sq = m2*m2;
    // symmetric mass ratio
    REAL8 eta  = m1*m2/msq;
    if (eta>0.25){
      printf("Truncating eta from above to 0.25. This should only be necessary in some rounding corner cases, but better check your m1 and m2 inputs...");
      eta = 0.25;
    }
    if (eta<0.0){
      printf("Truncating negative eta to 0.0. This should only be necessary in some rounding corner cases, but better check your m1 and m2 inputs...");
      eta = 0.0;
    }
    REAL8 eta2 = eta*eta;
    REAL8 eta3 = eta2*eta;
    //REAL8 eta4 = eta2*eta2;
    // spin variables (in m = 1 units)
    REAL8 S1    = chi1*m1sq/msq; // spin angular momentum 1
    REAL8 S2    = chi2*m2sq/msq; // spin angular momentum 2
    REAL8 Stot  = S1+S2;         // total spin
    REAL8 Shat  = (chi1*m1sq+chi2*m2sq)/(m1sq+m2sq); // effective spin, = msq*Stot/(m1sq+m2sq)
    REAL8 Shat2 = Shat*Shat;
    REAL8 Shat3 = Shat2*Shat;
    //REAL8 Shat4 = Shat2*Shat2;
    // spin difference, assuming m1>m2
    REAL8 chidiff  = chi1 - chi2;
    if (m2>m1){ // fit assumes m1>m2
      chidiff = -chidiff;
     }
    REAL8 chidiff2 = chidiff*chidiff;
    // typical squareroots and functions of eta
    //REAL8 sqrt2 = pow(2.,0.5);
    REAL8 sqrt3 = pow(3.,0.5);
    REAL8 sqrt1m4eta = pow((1. - 4.*eta),0.5);


    REAL8 a20 = 5.24;
    REAL8 a30 = 1.3;
    REAL8 a50 = 2.88;
    REAL8 b10 = -0.194;
    REAL8 b20 = 0.0851;
    REAL8 b30 = 0.00954;
    REAL8 b50 = -0.579;
    // fit coefficients from Tables I-IV of 1611.00332v2
    // values at increased numerical precision copied from
    // https://git.ligo.org/uib-papers/finalstate2016/blob/master/LALInference/FinalSpinUIB2016v2_pyform_coeffs.txt
    // git commit f490774d3593adff5bb09ae26b7efc6deab76a42
    REAL8 a2 = 3.8326341618708577;
    REAL8 a3 = -9.487364155598392;
    REAL8 a5 = 2.5134875145648374;
    REAL8 b1 = 1.0009563702914628;
    REAL8 b2 = 0.7877509372255369;
    REAL8 b3 = 0.6540138407185817;
    REAL8 b5 = 0.8396665722805308;
    REAL8 f21 = 8.77367320110712;
    REAL8 f31 = 22.830033250479833;
    REAL8 f50 = 1.8804718791591157;
    REAL8 f11 = 4.409160174224525;
    REAL8 f52 = 0.;
    REAL8 d10 = 0.3223660562764661;
    REAL8 d11 = 9.332575956437443;
    REAL8 d20 = -0.059808322561702126;
    REAL8 d30 = 2.3170397514509933;
    REAL8 d31 = -3.2624649875884852;
    REAL8 f12 = 0.5118334706832706;
    REAL8 f22 = -32.060648277652994;
    REAL8 f32 = -153.83722669033995;
    REAL8 f51 = -4.770246856212403;

    // Calculate the fit for the Lorb' quantity from Eq. (16) of 1611.00332
    REAL8 Lorb = (2.*sqrt3*eta + a20*a2*eta2 + a30*a3*eta3)/(1. + a50*a5*eta) + (b10*b1*Shat*(f11*eta + f12*eta2 + (64. - 16.*f11 - 4.*f12)*eta3) + b20*b2*Shat2*(f21*eta + f22*eta2 + (64. - 16.*f21 - 4.*f22)*eta3) + b30*b3*Shat3*(f31*eta + f32*eta2 + (64. - 16.*f31 - 4.*f32)*eta3))/(1. + b50*b5*Shat*(f50 + f51*eta + f52*eta2 + (64. - 64.*f50 - 16.*f51 - 4.*f52)*eta3)) + d10*sqrt1m4eta*eta2*(1. + d11*eta)*chidiff + d30*Shat*sqrt1m4eta*eta3*(1. + d31*eta)*chidiff + d20*eta3*chidiff2;
    // Convert to actual final spin
    REAL8 chif = Lorb + Stot;
    return chif;

}

/** @} */

/** @} */
