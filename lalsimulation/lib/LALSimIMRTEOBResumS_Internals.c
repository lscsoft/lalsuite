
/*
 * This file is part of TEOBResumS
 *
 * Copyright (C) 2017-2018  Alessandro Nagar, Sebastiano Bernuzzi,
 * Sarp Ackay, Gregorio Carullo, Walter Del Pozzo, Ka Wa Tsang, Michalis Agathos
 * LALSimulation implementation by Michalis Agathos
 *
 * TEOBResumS is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * TEOBResumS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see http://www.gnu.org/licenses/.
 *
 */

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include <complex.h>
#include <math.h>

#include <gsl/gsl_const.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv.h>

#include "LALSimTEOBResumS.h"


#define TEOB_NV (21) /* temp arrays */

/* hlmNewt coefficients for amplitude */
static const REAL8 ChlmNewt_ampli[35] = {2.1137745587232057, 6.341323676169617, 0.1412325034218127, 1.7864655618418102, 4.9229202032627635, 0.023872650234580958, 0.2250735048768909, 1.7053495827316825, 4.763908164911493, 0.001122165903318321, 0.06333806741197714, 0.2945348827200268, 1.755276012972272, 5.0817902739730565, 0.00014954736544380072, 0.005296595280255638, 0.10342548105284892, 0.3713362832603404, 1.8983258440274462, 5.727111757630886, 5.184622059790144e-06, 0.0012191691413436815, 0.011783593824950922, 0.14639129995388936, 0.4653654097044924, 2.125638973894669, 6.685178460621457, 5.54485779151375621e-7, 0.0000763473331250837455, 0.00353250998285463003,
    0.0204988821800401766, 0.19579402814926015, 0.584571015778149663, 2.44207899966355693, 7.9955401278763745};

/* hlmNewt additive coefficients for phase */
static const REAL8 ChlmNewt_phase[35] = {4.71238898038469, 3.141592653589793, 4.71238898038469, 0.0, 1.5707963267948966, 1.5707963267948966, 0.0, 4.71238898038469, 3.141592653589793, 1.5707963267948966, 3.141592653589793, 4.71238898038469, 0.0, 1.5707963267948966, 4.71238898038469, 3.141592653589793, 1.5707963267948966, 0.0, 4.71238898038469, 3.141592653589793, 4.71238898038469, 0.0, 1.5707963267948966, 3.141592653589793, 4.71238898038469, 0.0, 1.5707963267948966, 4.7123889803846898577, 0.0, 1.5707963267948966192, 3.1415926535897932385,
    4.7123889803846898577, 0.0, 1.5707963267948966192, 3.1415926535897932385};

/* Coefficients for Newtonian flux */
static const double CNlm[35] = {
    0.17777777777777778, 6.4,
    0.0007936507936507937, 0.5079365079365079, 8.678571428571429,
    2.2675736961451248e-05, 0.008062484252960444, 1.0414285714285714, 14.447971781305114,
    5.010421677088344e-08, 0.0006384836014465644, 0.03106534090909091, 1.9614216236438458, 25.688197074915823,
    8.898517797026696e-10, 4.464920289836115e-06, 0.003830520129221428, 0.08778390483440988, 3.584607999290539, 46.98226573426573,
    1.0695333890657086e-12, 2.3656367312710969e-07, 4.972309783123969e-05, 0.013643024456211269, 0.21542115380351798, 6.472046810332524, 87.1329124642076,
    1.2233225038333268e-14, 9.27700678929842e-10, 4.468579053461083e-06, 0.0002675102834551229, 0.03813282951314997, 0.4894825318738884, 11.627213264559293, 162.79300083906728
};

/* EOB nonspinning Hamiltonian */
void eob_ham(REAL8 nu, REAL8 r, REAL8 pphi, REAL8 prstar, REAL8 A, REAL8 dA,
             REAL8 *H, /* real EOB Hamiltonian divided by mu=m1m2/(m1+m2) */
             REAL8 *Heff, /* effective EOB Hamiltonian (divided by mu) */
             REAL8 *dHeff_dr, /* drvt Heff,r */
             REAL8 *dHeff_dprstar, /* drvt Heff,prstar */
             REAL8 *dHeff_dpphi /* drvt Heff,pphi */
)
{
    const REAL8 z3 = 2.0*nu*(4.0-3.0*nu);
    const REAL8 pphi2 = SQ(pphi);
    const REAL8 u  = 1./r;
    const REAL8 u2 = SQ(u);
    const REAL8 u3 = u2*u;
    const REAL8 prstar2 = SQ(prstar);
    const REAL8 prstar3 = prstar2*prstar;
    const REAL8 prstar4 = prstar2*prstar2;

    *Heff = sqrt(A*(1.0 + pphi2*u2) + prstar2 + z3*A*u2*prstar4);
    *H    = sqrt( 1.0 + 2.0*nu*(*Heff - 1) )/nu;

    if (dHeff_dr != NULL)      *dHeff_dr      = 0.5*(dA + (pphi2 + z3*prstar4)*(dA*u2 - 2*A*u3))/(*Heff);
    if (dHeff_dprstar != NULL) *dHeff_dprstar = (prstar + z3*2.0*A*u2*prstar3)/(*Heff);
    if (dHeff_dpphi != NULL)   *dHeff_dpphi   = A*pphi*u2/(*Heff);
    return;
}

/* EOB spinning Hamiltonian */
void eob_ham_s(REAL8 nu,
               REAL8 r,
               REAL8 rc,
               REAL8 drc_dr,
               REAL8 pphi,
               REAL8 prstar,
               REAL8 S,
               REAL8 Sstar,
               REAL8 chi1,
               REAL8 chi2,
               REAL8 X1,
               REAL8 X2,
               REAL8 aK2,
               REAL8 c3,
               REAL8 A,
               REAL8 dA,
               REAL8 *H,             /* real EOB Hamiltonian divided by mu=m1m2/(m1+m2) */
               REAL8 *Heff,          /* effective EOB Hamiltonian (divided by mu) */
               REAL8 *Heff_orb,
               REAL8 *dHeff_dr,      /* drvt Heff,r */
               REAL8 *dHeff_dprstar, /* drvt Heff,prstar */
               REAL8 *dHeff_dpphi,    /* drvt Heff,pphi */
               REAL8 *d2Heff_dprstar20
               )
{
    /* Shorthands */
    const REAL8 z3      = 2.0*nu*(4.0-3.0*nu);
    const REAL8 pphi2    = SQ(pphi);
    const REAL8 prstar2 = SQ(prstar);
    const REAL8 UNUSED prstar3 = prstar2*prstar;
    const REAL8 prstar4 = prstar2*prstar2;
    const REAL8 uc  = 1./rc;
    const REAL8 uc2 = uc*uc;
    const REAL8 uc3 = uc2*uc;

    /* Compute spin-related functions*/
    REAL8 ggm[14];
    eob_dyn_s_GS(r, rc, drc_dr, aK2, prstar, pphi, nu, chi1, chi2, X1, X2, c3, ggm);
    const REAL8 GS              = ggm[2];
    const REAL8 GSs             = ggm[3];
    const REAL8 dGS_dprstar     = ggm[4];
    const REAL8 dGSs_dprstar    = ggm[5];
    const REAL8 dGS_dr          = ggm[6];
    const REAL8 dGSs_dr         = ggm[7];
    const REAL8 dGSs_dpphi      = ggm[9];
    const REAL8 d2GS_dprstar20  = ggm[12];
    const REAL8 d2GSs_dprstar20 = ggm[13];

    /* Compute Hamiltonian and its derivatives */
    *Heff_orb         = sqrt( prstar2+A*(1. + pphi2*uc2 +  z3*prstar4*uc2) );
    *Heff             = *Heff_orb + (GS*S + GSs*Sstar)*pphi;
    *H                = sqrt( 1. + 2.*nu*(*Heff - 1.) )/nu;
    if (dHeff_dr != NULL)         *dHeff_dr         = pphi*(dGS_dr*S + dGSs_dr*Sstar) + 1./(2.*(*Heff_orb))*( dA*(1. + pphi2*uc2 + z3*prstar4*uc2) - 2.*A*uc3*drc_dr*(pphi2 + z3*prstar4) );
    if (dHeff_dprstar != NULL)    *dHeff_dprstar    = pphi*(dGS_dprstar*S + dGSs_dprstar*Sstar) + (prstar/(*Heff_orb))*(1. + 2.*A*uc2*z3*prstar2);
    if (d2Heff_dprstar20 != NULL) *d2Heff_dprstar20 = pphi*(d2GS_dprstar20*S + d2GSs_dprstar20*Sstar) +  (1./(*Heff_orb))*(1. + 2.*A*uc2*z3*prstar2); /* second derivative of Heff wrt to pr_star neglecting all pr_star^2 terms */
    if (dHeff_dpphi != NULL)      *dHeff_dpphi      = GS*S + (GSs + pphi*dGSs_dpphi)*Sstar + pphi*A*uc2/(*Heff_orb);

    return;
}


/* Function for root finder: Derivative of the effective Hamiltonian */
struct DHeff0_tmp_params
{
    REAL8 rorb, A, dA, rc, drc_dr, ak2, S, Ss, nu, chi1, chi2, X1, X2, c3;
};


REAL8 eob_dyn_DHeff0(REAL8 x, void *params)
{

    struct DHeff0_tmp_params *p
    = (struct DHeff0_tmp_params *) params;

    REAL8 rorb   = p->rorb;
    REAL8 A      = p->A;
    REAL8 dA     = p->dA;
    REAL8 rc     = p->rc;
    REAL8 drc_dr = p->drc_dr;
    REAL8 ak2    = p->ak2;
    REAL8 S      = p->S;
    REAL8 Ss     = p->Ss;
    REAL8 nu     = p->nu;
    REAL8 chi1   = p->chi1;
    REAL8 chi2   = p->chi2;
    REAL8 X1     = p->X1;
    REAL8 X2     = p->X2;
    REAL8 c3     = p->c3;

    REAL8 ggm0[14];
    eob_dyn_s_GS(rorb, rc, drc_dr, ak2, 0., x, nu, chi1, chi2, X1, X2, c3, ggm0);
    REAL8 dGS_dr  = ggm0[6];
    REAL8 dGSs_dr = ggm0[7];

    REAL8 x2  = SQ(x);
    REAL8 uc  = 1./rc;
    REAL8 uc2 = SQ(uc);
    REAL8 uc3 = uc2*uc;

    /* Orbital circular effective Hamiltonian */
    REAL8 Horbeff0 = sqrt(A*(1. + x2*uc2));
    REAL8 dHeff_dr = x*(dGS_dr*S + dGSs_dr*Ss) + 1./(2.*Horbeff0)*( dA*(1. + x2*uc2) - 2.*A*uc3*drc_dr*x2);
    return dHeff_dr;
}


/* Root finder: Compute minimum of Heff0 */
REAL8 eob_dyn_bisecHeff0_s(REAL8 nu,
                           REAL8 chi1,
                           REAL8 chi2,
                           REAL8 X1,
                           REAL8 X2,
                           REAL8 c3,
                           REAL8 pph,
                           REAL8 rorb,
                           REAL8 A,
                           REAL8 dA,
                           REAL8 rc,
                           REAL8 drc_dr,
                           REAL8 ak2,
                           REAL8 S,
                           REAL8 Ss)
{

    int status;
    INT4 iter = 0;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;

    double r;
    double x_lo =  0.5*(double) pph, x_hi = 1.5*(double) pph;

    gsl_function F;
    struct DHeff0_tmp_params p = {rorb,A,dA,rc,drc_dr,ak2,S,Ss,nu,chi1,chi2,X1,X2,c3};
    F.function = &eob_dyn_DHeff0;

    F.params = &p;
    T = gsl_root_fsolver_bisection;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, x_lo, x_hi);
    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        r      = gsl_root_fsolver_root (s);
        x_lo   = gsl_root_fsolver_x_lower (s);
        x_hi   = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (x_lo, x_hi, 0, TOLERANCE);
    }
    while (status == GSL_CONTINUE && iter < MAX_ITER);
    gsl_root_fsolver_free (s);

    return (REAL8) r;
}

/* Initial radius from initial frequency using Kepler's law */
REAL8 eob_dyn_r0_Kepler (REAL8 f0)
{
    const REAL8 omg_orb0 = LAL_PI*f0; // =2*LAL_PI*(f_orb0/2)
    return pow(omg_orb0, -2./3.);
}

/* Initial radius from initial frequency using EOB circular dynamics */
REAL8 eob_dyn_r0_eob (REAL8 f0, LALTEOBResumSDynamics *dyn)
{
    const REAL8 omg_orb0 = LAL_PI*f0;
    const REAL8 r0_kepl  = eob_dyn_r0_Kepler(f0);
    return eob_dyn_bisecOmegaorb0(dyn, omg_orb0, r0_kepl);
}

/* Function for root finder: omega = omega_circ */
struct Omegaorb0_tmp_params {
    REAL8 omg_orb0;
    LALTEOBResumSDynamics *dyn;
};

REAL8 eob_dyn_Omegaorb0(REAL8 r, void *params)
{

    /* Unpack parameters */
    struct Omegaorb0_tmp_params *p = (struct Omegaorb0_tmp_params *) params;
    REAL8   omg_orb0 = p->omg_orb0;
    LALTEOBResumSDynamics *dyn = p->dyn;

    const REAL8 nu    = dyn->nu;
    const REAL8 X1    = dyn->X1;
    const REAL8 X2    = dyn->X2;
    const REAL8 chi1  = dyn->chi1;
    const REAL8 chi2  = dyn->chi2;
    const REAL8 a1    = dyn->a1;
    const REAL8 a2    = dyn->a2;
    const REAL8 aK2   = dyn->aK2;
    const REAL8 S     = dyn->S;
    const REAL8 Sstar = dyn->Sstar;
    const REAL8 c3    = dyn->cN3LO;
    const REAL8 C_Q1  = dyn->C_Q1;
    const REAL8 C_Q2  = dyn->C_Q2;
    const REAL8 C_Oct1 = dyn->C_Oct1;
    const REAL8 C_Oct2 = dyn->C_Oct2;
    const REAL8 C_Hex1 = dyn->C_Hex1;
    const REAL8 C_Hex2 = dyn->C_Hex2;

    const int usetidal = dyn->use_tidal;
    const int usespins = dyn->use_spins;

    REAL8 A,B,dA,rc,drc_dr,G,dG_dr,uc,uc2,dAuc2_dr,j02,j0,H,Heff,Heff_orb,dHeff_dj0,omg_orb;
    REAL8 pl_hold,a_coeff,b_coeff,c_coeff,Delta,sol_p,sol_m;
    REAL8 ggm[14];

    /* Computing metric, centrifugal radius and ggm functions*/
    if(usespins)
    {
        eob_metric_s(r,dyn, &A, &B, &dA, &pl_hold, &pl_hold);
        dyn->eob_dyn_s_get_rc(r, nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, &rc, &drc_dr, &pl_hold);
        eob_dyn_s_GS(r, rc, drc_dr, aK2, 0.0, 0.0, nu, chi1, chi2, X1, X2, c3, ggm);
        G     = ggm[2]*S + ggm[3]*Sstar;    // tildeG = GS*S+GSs*Ss
        dG_dr = ggm[6]*S + ggm[7]*Sstar;
    }
    else
    {
        eob_metric(r ,dyn, &A, &B, &dA, &pl_hold, &pl_hold);
        rc     = r;   //Nonspinning case: rc = r; G = 0;
        drc_dr = 1;
        G      = 0.0;
        dG_dr  = 0.0;
    }

    /* Auxiliary variables*/
    uc       = 1./rc;
    uc2      = uc*uc;
    dAuc2_dr = uc2*(dA-2*A*uc*drc_dr);

    /* Circular angular momentum */
    if (usespins)
    {

        /* Quadratic equation \f$a*x^2+b*x+c=0\f$ */
        a_coeff = SQ(dAuc2_dr)  - 4*A*uc2*SQ(dG_dr);
        b_coeff = 2*dA*dAuc2_dr - 4*A*SQ(dG_dr);
        c_coeff = SQ(dA);

        Delta = SQ(b_coeff) - 4*a_coeff*c_coeff;

        if (S==0 && Sstar==0)
            Delta=0;             // dG_dr=0 -> Set Delta=0 to avoid num. errors

        sol_p   = (-b_coeff + sqrt(Delta))/(2*a_coeff);
        sol_m   = (-b_coeff - sqrt(Delta))/(2*a_coeff);

        if (dG_dr > 0)
            j02 = sol_p;
        else
            j02 = sol_m;

    } else {
        /* Linear equation \f$a*x+b=0\f$ */
        a_coeff = dAuc2_dr;
        b_coeff = dA;
        j02 = -b_coeff/a_coeff;
    }

    j0 = sqrt(j02);

    /* Circular Hamiltonians */
    Heff_orb = sqrt(A*(1+j02*uc2));
    Heff     = Heff_orb + j0*G;
    H        = sqrt(1+2*nu*(Heff-1))/nu;

    /* Circular orbital frequency */
    dHeff_dj0 = G + A*j0*uc2/Heff_orb;
    omg_orb   = dHeff_dj0/nu/H;

    /* Subtraction of initial evolution frequency */
    return (omg_orb - omg_orb0);
}

/* Root finder: Compute r0 such that omg_orb = omg_orb0 */
REAL8 eob_dyn_bisecOmegaorb0(LALTEOBResumSDynamics *dyn, REAL8 omg_orb0, REAL8 r0_kepl)
{
    int status;
    INT4 iter = 0;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    REAL8 r0;
    double x_lo = 0.5*r0_kepl, x_hi = 1.5*r0_kepl;
    gsl_function F;

    struct  Omegaorb0_tmp_params p = {omg_orb0,dyn};

    F.function = &eob_dyn_Omegaorb0;
    F.params = &p;
    T = gsl_root_fsolver_bisection;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, x_lo, x_hi);
    do {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        r0     = (REAL8) gsl_root_fsolver_root (s);
        x_lo   = gsl_root_fsolver_x_lower (s);
        x_hi   = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (x_lo, x_hi, 0, TOLERANCE);
    }
    while (status == GSL_CONTINUE && iter < MAX_ITER);
    gsl_root_fsolver_free (s);

    return r0;
}

/* Initial conditions calculation for non-spinning systems */
/* Post-post-circular initial data at separation r0
 r0       => relative separation
 pph0     => post-post-circular angular momentum
 pr       => post-circular radial momentum
 pprstar0 => post-circular r*-conjugate radial momentum
 j0       => circular angular momentum

 Three-step procedure
 1. Compute j0                         =>           circular ID, j!=0, pr =0
 2. From j0, compute pr*               =>      post circular ID, j!=0, pr !=0
 3. From pr* and j0, re-compute pph0   => post-post-circular ID, pph0!=j!=0, pr !=0
 */
void eob_dyn_ic(REAL8 r0, LALTEOBResumSDynamics *dyn, REAL8 y_init[])
{
    const REAL8 nu = dyn->nu;
    const REAL8 z3 = 2.0*nu*(4.0-3.0*nu);

    /* Build a small grid */
#define TEOB_IC_N (6)
    const REAL8 dr = 1e-10;

    REAL8 r[2*TEOB_IC_N], dA[2*TEOB_IC_N], j[2*TEOB_IC_N], j2[2*TEOB_IC_N], djdr[2*TEOB_IC_N]; /* j:angular momentum */
    REAL8 E0[2*TEOB_IC_N], Omega_j[2*TEOB_IC_N];
    REAL8 Fphi[2*TEOB_IC_N], Ctmp[2*TEOB_IC_N], prstar[2*TEOB_IC_N], pr[2*TEOB_IC_N], pph[2*TEOB_IC_N], dprstardr[2*TEOB_IC_N];

    REAL8 A, B, d2A, dB;
    REAL8 r2, r3, j3;
    REAL8 H0eff, H0, psi, r_omega, v_phi, jhat, x, dprstardt;

    for (int i = 0; i < 2*TEOB_IC_N; i++) {
        r[i] = r0+(i-TEOB_IC_N+1)*dr;
        r2   = SQ(r[i]);
        r3   = r2*r[i];

        /* Compute metric  */
        eob_metric(r[i], dyn, &A, &B, &dA[i], &d2A, &dB);

        /* Angular momentum for circular orbit: circular ID  */
        j2[i]   =  r3*dA[i]/(2.*A-r[i]*dA[i]);
        j[i]    =  sqrt(j2[i]);
        j3      =  j2[i]*j[i];
        djdr[i] = -j3/r3*( 2.0 - 3.0*A/(r[i]*dA[i]) - A*d2A/(dA[i]*dA[i]) );

        /* For circular orbit at r0=r(TEOB_IC_N)  */
        H0eff      = sqrt(A*(1.0 + j2[i]/r2));                     /* effective Hamiltonian H_0^eff  */
        E0[i]      = sqrt(1.0 + 2.0*nu*(H0eff - 1.0) );            /* real Hamiltonian      H_0  */
        H0         = E0[i]/nu;                                     /* H_0/nu  */
        Omega_j[i] = A*j[i]/(nu*r2*H0*H0eff);                      /* Orbital frequency (from Hamilton's equation)  */
        psi        = 2.*(1.0 + 2.0*nu*(H0eff - 1.0))/(r2*dA[i]);   /* correction factor to the radius  */
        r_omega    = r[i]*cbrt(psi);                               /* EOB-corrected radius  */
        v_phi      = Omega_j[i]*r_omega;                           /* "corrected" azimuthal velocity such that Kepler's law is satisfied, r_omg^3 Omg_i^2 = 1  */
        x          = v_phi * v_phi;
        jhat       = j[i]/(r_omega*v_phi);                         /* Newton-normalized angular momentum  */

        Fphi[i] = eob_flx_Flux(x,Omega_j[i],r_omega,E0[i],H0eff,jhat,r[i], 0,0,dyn);

        /* Radial momentum conjugate to r*: post-circular ID  */
        Ctmp[i]   = sqrt(B/A)*nu*H0*H0eff;
        prstar[i] = Ctmp[i]*Fphi[i]/djdr[i];

        /* Radial momentum conjugate to r  */
        pr[i] = prstar[i]*sqrt(B/A);

    }

    /* prstar by finite diff. */
    D0(prstar, dr, 2*TEOB_IC_N, dprstardr);

    int i = TEOB_IC_N-1;
    dprstardt = dprstardr[i] * Fphi[i]/djdr[i];
    pph[i] = j[i]*sqrt(1. + 2.*Ctmp[i]/dA[i]*dprstardt - z3*gsl_pow_int(prstar[i],4)/j2[i]);

    y_init[TEOB_ID_RAD]    = r[TEOB_IC_N-1];
    y_init[TEOB_ID_PHI]    = 0.;
    y_init[TEOB_ID_PPHI]   = pph[TEOB_IC_N-1];
    y_init[TEOB_ID_PRSTAR] = prstar[TEOB_IC_N-1];
    y_init[TEOB_ID_PR]     = pr[TEOB_IC_N-1];
    y_init[TEOB_ID_J]      = j[TEOB_IC_N-1];
    y_init[TEOB_ID_E0]     = E0[TEOB_IC_N-1];
    y_init[TEOB_ID_OMGJ]   = Omega_j[TEOB_IC_N-1];

}


/* Initial conditions calculation for spinning systems */
void eob_dyn_ic_s(REAL8 r0, LALTEOBResumSDynamics *dyn, REAL8 y_init[])
{
    const REAL8 nu   = dyn->nu;
    const REAL8 chi1 = dyn->chi1;
    const REAL8 chi2 = dyn->chi2;
    const REAL8 S1   = dyn->S1;
    const REAL8 S2   = dyn->S2;
    const REAL8 c3   = dyn->cN3LO;
    const REAL8 X1   = dyn->X1;
    const REAL8 X2   = dyn->X2;
    const REAL8 a1   = dyn->a1;
    const REAL8 a2   = dyn->a2;
    const REAL8 aK2  = dyn->aK2;
    const REAL8 C_Q1 = dyn->C_Q1;
    const REAL8 C_Q2 = dyn->C_Q2;
    const REAL8 C_Oct1 = dyn->C_Oct1;
    const REAL8 C_Oct2 = dyn->C_Oct2;
    const REAL8 C_Hex1 = dyn->C_Hex1;
    const REAL8 C_Hex2 = dyn->C_Hex2;

    const REAL8 S  = S1 + S2;
    const REAL8 Ss = X2*a1 + X1*a2;
    const REAL8 z3 = 2.0*nu*(4.0-3.0*nu);

    /* Build a small grid */
#define TEOB_IC_S_N (6)
    const REAL8 dr = 1e-4; /* do not change this */

    REAL8 r[2*TEOB_IC_S_N], dA[2*TEOB_IC_S_N], j[2*TEOB_IC_S_N]; /* j:angular momentum */
    REAL8 E0[2*TEOB_IC_S_N], Omega_j[2*TEOB_IC_S_N];
    REAL8 Fphi[2*TEOB_IC_S_N], UNUSED Ctmp[2*TEOB_IC_S_N], prstar[2*TEOB_IC_S_N], pr[2*TEOB_IC_S_N], pph[2*TEOB_IC_S_N];
    REAL8 rc[2*TEOB_IC_S_N], drc_dr[2*TEOB_IC_S_N], d2rc_dr2[2*TEOB_IC_S_N]; //, drc[2*TEOB_IC_S_N];
    REAL8 A[2*TEOB_IC_S_N],B[2*TEOB_IC_S_N],d2A[2*TEOB_IC_S_N],dB, sqrtAbyB;
    REAL8 pphorb, uc, uc2, psic, r_omg, v_phi, jhat, x, Omg;
    REAL8 UNUSED H0eff, H0, Horbeff0, Heff0, one_H0, dHeff_dprstarbyprstar, dHeff_dpph, Heff, H, Horbeff;
    REAL8 ggm0[14], GS_0, GSs_0, dGS_dr_0, dGSs_dr_0, dGSs_dpph_0, dGS_dprstarbyprstar_0, dGSs_dprstarbyprstar_0, GS, GSs, dGS_dr, dGSs_dr;
    REAL8 C0;
    REAL8 Gtilde, dGtilde_dr, duc_dr;

    int i;
    for (i = 0; i < 2*TEOB_IC_S_N; i++) {
        r[i] = r0+(i-TEOB_IC_S_N+1)*dr;

        /* Compute metric  */
        eob_metric_s(r[i], dyn, &A[i],&B[i],&dA[i],&d2A[i],&dB);

        /* Compute minimum of Heff0 using bisection method */
        pphorb = r[i]/sqrt(r[i]-3.);
        dyn->eob_dyn_s_get_rc(r[i], nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1,
                         C_Hex2, dyn->use_tidal, &rc[i], &drc_dr[i], &d2rc_dr2[i]);
        pph[i] = eob_dyn_bisecHeff0_s(nu,chi1,chi2,X1,X2,c3, pphorb, r[i], A[i],
                                      dA[i], rc[i], drc_dr[i], aK2, S, Ss);

    }

    /* Post-circular initial conditions */

    /* pph by finite diff. */
    REAL8 dpph_dr[2*TEOB_IC_S_N];
    D0(pph, dr, 2*TEOB_IC_S_N, dpph_dr);

    for (i = 0; i < 2*TEOB_IC_S_N; i++) {

        sqrtAbyB = sqrt(A[i]/B[i]);
        uc  = 1./rc[i];
        uc2 = uc*uc;

        /* Orbital effective Hamiltonian */
        Horbeff0 = sqrt(A[i]*(1. + SQ(pph[i])*uc2));

        /* Compute gyro-gravitomagnetic coupling functions */
        eob_dyn_s_GS(r[i], rc[i], drc_dr[i], aK2, 0, pph[i], nu, chi1, chi2, X1, X2, c3, ggm0);
        GS_0                   = ggm0[2];
        GSs_0                  = ggm0[3];
        dGS_dr_0               = ggm0[6];
        dGSs_dr_0              = ggm0[7];
        dGSs_dpph_0            = ggm0[9];
        dGS_dprstarbyprstar_0  = ggm0[10];
        dGSs_dprstarbyprstar_0 = ggm0[11];

        /* Final effective Hamiltonian */
        Heff0 = (GS_0*S + GSs_0*Ss)*pph[i] + Horbeff0;

        /* Real Hamiltonian: beware that this is NOT divided by nu */
        H0     = sqrt( 1. + 2.*nu*(Heff0 - 1.));
        one_H0 = 1./H0;

        /* Get gyro-gravitomagnetic (derivative) functions */
        dHeff_dprstarbyprstar = pph[i]*(dGS_dprstarbyprstar_0*S + dGSs_dprstarbyprstar_0*Ss) + 1./Horbeff0;

        C0 = sqrtAbyB*one_H0*dHeff_dprstarbyprstar;

        /* Orbital frequency for circular orbit */
        dHeff_dpph = GS_0*S + (GSs_0 + pph[i]*dGSs_dpph_0)*Ss + pph[i]*A[i]*uc2/Horbeff0;
        Omg        = one_H0*dHeff_dpph;

        /* Flux */
        Gtilde     =  GS_0*S     + GSs_0*Ss;
        dGtilde_dr =  dGS_dr_0*S + dGSs_dr_0*Ss;
        duc_dr     = -uc2*drc_dr[i];
        psic       = (duc_dr + dGtilde_dr*rc[i]*sqrt(A[i]/(SQ(pph[i])) + A[i]*uc2)/A[i])/(-0.5*dA[i]);
        r_omg      =  pow((pow(rc[i]*rc[i]*rc[i]*psic,-1./2)+Gtilde)*one_H0,-2./3.);
        v_phi      =  r_omg*Omg;
        x          =  v_phi*v_phi;
        jhat       =  pph[i]/(r_omg*v_phi);  /* Newton-normalized angular momentum */

        Fphi[i]    = eob_flx_Flux_s(x, Omg, r_omg, H0, Heff0, jhat, r[i], 0., 0., dyn);
        prstar[i]  = Fphi[i]/(dpph_dr[i]*C0);
        pr[i]      = prstar[i]/sqrtAbyB;

        j[i]       = pph[i];
        E0[i]      = H0;
        Omega_j[i] = Omg;

    }

#if (POSTPOSTCIRCULAR)

    /* Post-post-circular initial data */

    REAL8 dpi1bydj, dprstardr[2*TEOB_IC_S_N],djdr[2*TEOB_IC_S_N];
    D0(prstar, dr, 2*TEOB_IC_S_N, dprstardr);
    D0(j, dr, 2*TEOB_IC_S_N, djdr);

    REAL8 dpi1dt, prstar4, a,b,c;

    //for (i = 0; i < 2*TEOB_IC_S_N; i++) { //No need of the loop here
    i = TEOB_IC_S_N-1;

    sqrtAbyB = sqrt(A[i]/B[i]);
    uc  = 1./rc[i];
    uc2 = uc*uc;

    dpi1bydj = dprstardr[i]/djdr[i];
    dpi1dt   = dpi1bydj*Fphi[i];
    prstar4  = prstar[i]*prstar[i]*prstar[i]*prstar[i];

    /* Still circular, no pr* dependence here */
    Horbeff  = sqrt(A[i]*(1. + SQ(pph[i])*uc2));

    eob_dyn_s_GS(r[i], rc[i], drc_dr[i], aK2, 0, pph[i], nu, chi1, chi2, X1, X2, c3, ggm0);
    GS      = ggm0[2];
    GSs     = ggm0[3];
    dGS_dr  = ggm0[6];
    dGSs_dr = ggm0[7];

    /* Effective EOB energy */
    Heff     = (GS*S + GSs*Ss)*pph[i] + Horbeff;

    /* Total EOB energy */
    H        = sqrt( 1. + 2.*nu*(Heff - 1.));

    /* Setting up second order equation for the orbital angular momentum */
    a = -sqrtAbyB*uc2/(2.*H*Horbeff)*(dA[i]  - 2.*A[i]*uc*drc_dr[i]);
    b = -sqrtAbyB/H*(dGS_dr*S + dGSs_dr*Ss);
    c = -dpi1dt - sqrtAbyB/(2.*H*Horbeff)*(dA[i] + z3*prstar4*uc2*(dA[i] - 2.*A[i]*uc*drc_dr[i]));

    /* Fill out the array of the post-circular angular momentum */
    pph[i] = 0.5*(-b + sqrt(b*b-4*a*c))/a;

    //  }

#endif

    y_init[TEOB_ID_RAD]    = r[TEOB_IC_S_N-1];
    y_init[TEOB_ID_PHI]    = 0.;
    y_init[TEOB_ID_PPHI]   = pph[TEOB_IC_S_N-1];
    y_init[TEOB_ID_PRSTAR] = prstar[TEOB_IC_S_N-1];
    y_init[TEOB_ID_PR]     = pr[TEOB_IC_S_N-1];
    y_init[TEOB_ID_J]      = j[TEOB_IC_S_N-1];
    y_init[TEOB_ID_E0]     = E0[TEOB_IC_S_N-1];
    y_init[TEOB_ID_OMGJ]   = Omega_j[TEOB_IC_S_N-1];

    return;
}


/* r.h.s. of EOB Hamiltonian dynamics, no spins version */
int eob_dyn_rhs(REAL8 t, const REAL8 y[], REAL8 dy[], void *d)
{

    (void)(t); /* avoid unused parameter warning */
    LALTEOBResumSDynamics *dyn = d;

    const REAL8 nu = dyn->nu;
    const REAL8 z3 = 2.0*nu*(4.0-3.0*nu);

    /* Unpack y */
    const REAL8 UNUSED phi    = y[TEOB_EVOLVE_PHI];
    const REAL8 r      = y[TEOB_EVOLVE_RAD];
    const REAL8 pphi   = y[TEOB_EVOLVE_PPHI];
    const REAL8 prstar = y[TEOB_EVOLVE_PRSTAR];

    /* Compute EOB Metric */
    REAL8 A, B, dA, d2A, dB;
    eob_metric(r, d, &A, &B, &dA, &d2A, &dB);

    /* Compute Hamiltonian */
    REAL8 H, Heff, dHeff_dr,dHeff_dprstar;
    eob_ham(nu, r,pphi,prstar,A,dA, &H,&Heff,&dHeff_dr,&dHeff_dprstar,NULL);
    REAL8 E = nu*H;

    /* Shorthands */
    const REAL8 u  = 1./r;
    const REAL8 u2 = u*u;
    const REAL8 u3 = u2*u;
    const REAL8 pphi2    = SQ(pphi);
    const REAL8 prstar2  = prstar*prstar;
    const REAL8 prstar3  = prstar2*prstar;
    const REAL8 prstar4  = prstar3*prstar;
    const REAL8 sqrtAbyB = sqrt(A/B);
    const REAL8 divHE    = 1./(Heff*E);
    const REAL8 Omega    = A*pphi*u2*divHE;

    /* \f$d\phi/dt\f$ */
    dy[TEOB_EVOLVE_PHI] = Omega;

    /* \f$dr/dt\f$ (conservative part of) */
    dy[TEOB_EVOLVE_RAD] = sqrtAbyB*(prstar+4.0*nu*(4.0-3.0*nu)*A*u2*prstar3)*divHE;

    /* \f$dp_{r*}/dt\f$ (conservative part of) */
    dy[TEOB_EVOLVE_PRSTAR] = - 0.5*sqrtAbyB*( pphi2*u2*(dA-2.0*A*u) + dA + 2.0*nu*(4.0-3.0*nu)*(dA*u2 - 2.0*A*u3)*prstar4 )*divHE;

    /* Compute flux */
    const REAL8 sqrtW = sqrt(A*(1. + pphi2*u2));
    const REAL8 psi   = 2.*(1.0 + 2.0*nu*(sqrtW - 1.0))/(SQ(r)*dA);
    /*const REAL8 psi = 2.*(1.0 + 2.0*nu*(Heff - 1.0))/(r2*dA); */
    const REAL8 r_omega = r*cbrt(psi);
    const REAL8 v_phi   = r_omega*Omega;
    const REAL8 x       = v_phi * v_phi;
    const REAL8 jhat    = pphi/(r_omega*v_phi);
    const REAL8 tmpE    = 1./Heff+nu/(E*E);
    const REAL8 dprstar_dt    = dy[TEOB_EVOLVE_PRSTAR];
    const REAL8 dr_dt         = dy[TEOB_EVOLVE_RAD];
    const REAL8 ddotr_dr      = sqrtAbyB*( (prstar + z3*2.*A*u2*prstar3)*(0.5*(dA/A-dB/B)-dHeff_dr*tmpE)+ 2.0*z3*(dA*u2 - 2.*A*u3)*prstar3)*divHE;
    const REAL8 ddotr_dprstar = sqrtAbyB*( 1.+z3*6.*A*u2*prstar2-(prstar + z3*2.*A*u2*prstar3)*dHeff_dprstar*tmpE)*divHE;

    /* Approximate ddot(r) without Flux */
    const REAL8 ddotr = dprstar_dt*ddotr_dprstar + dr_dt*ddotr_dr;

    /* Compute flux and \f$dp_{\phi}/dt\f$ */
    if (dyn->noflux) dy[TEOB_EVOLVE_PPHI] = 0.;
    else            dy[TEOB_EVOLVE_PPHI] = eob_flx_Flux(x,Omega,r_omega,E,Heff,jhat,r, prstar,ddotr,dyn);

    if(dyn->store) {
        /* Store values */
        dyn->t = t;
        dyn->r = r;
        dyn->phi = y[TEOB_EVOLVE_PHI];
        dyn->pphi = pphi;
        dyn->prstar = prstar;
        dyn->Omg = Omega;
        dyn->Omg_orb = Omega;
        dyn->H = H;
        dyn->E = E;
        dyn->Heff = Heff;
        dyn->A = A;
        dyn->dA = dA;
        dyn->d2A = d2A;
        dyn->B = B;
        dyn->dB = dB;
        dyn->psi = psi;
        dyn->r_omega = r_omega;
        dyn->v_phi = v_phi;
        dyn->jhat = jhat;
        dyn->ddotr = ddotr;
    }

    return GSL_SUCCESS;

}


/* r.h.s. of EOB Hamiltonian dynamics, spins version */
int eob_dyn_rhs_s(REAL8 t, const REAL8 y[], REAL8 dy[], void *d)
{

    (void)(t); /* avoid unused parameter warning */
    LALTEOBResumSDynamics *dyn = d;

    /* Unpack values */
    const REAL8 nu    = dyn->nu;
    const REAL8 S     = dyn->S;
    const REAL8 Sstar = dyn->Sstar;
    const REAL8 chi1  = dyn->chi1;
    const REAL8 chi2  = dyn->chi2;
    const REAL8 X1    = dyn->X1;
    const REAL8 X2    = dyn->X2;
    const REAL8 c3    = dyn->cN3LO;
    const REAL8 aK2   = dyn->aK2;
    const REAL8 a1    = dyn->a1;
    const REAL8 a2    = dyn->a2;
    const REAL8 C_Q1  = dyn->C_Q1;
    const REAL8 C_Q2  = dyn->C_Q2;
    const REAL8 C_Oct1 = dyn->C_Oct1;
    const REAL8 C_Oct2 = dyn->C_Oct2;
    const REAL8 C_Hex1 = dyn->C_Hex1;
    const REAL8 C_Hex2 = dyn->C_Hex2;
    const int usetidal = dyn->use_tidal;
    const int UNUSED usespins = dyn->use_spins;

    /* Shorthands */
    const REAL8 r      = y[TEOB_EVOLVE_RAD];
    const REAL8 prstar = y[TEOB_EVOLVE_PRSTAR];
    const REAL8 pphi   = y[TEOB_EVOLVE_PPHI];
    const REAL8 pphi2  = pphi*pphi;

    /* Compute Metric */
    REAL8 A, B, dA, d2A, dB;
    eob_metric_s(r, d, &A, &B, &dA, &d2A, &dB);

    /* Compute centrifugal radius */
    REAL8 rc, drc_dr, d2rc_dr;
    dyn->eob_dyn_s_get_rc(r, nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, &rc, &drc_dr, &d2rc_dr);
    const REAL8 uc     = 1./rc;
    const REAL8 uc2    = uc*uc;
    const REAL8 UNUSED uc3    = uc2*uc;

    /* Compute Hamiltonian */
    REAL8 Heff_orb, Heff, H, dHeff_dr, dHeff_dprstar, d2Heff_dprstar20, dHeff_dpphi;
    eob_ham_s(nu, r, rc, drc_dr, pphi, prstar, S, Sstar, chi1, chi2, X1, X2, aK2, c3, A, dA,
              &H, &Heff, &Heff_orb, &dHeff_dr, &dHeff_dprstar, &dHeff_dpphi, &d2Heff_dprstar20);

    /* H follows the same convention of Heff, i.e. it is the energy per unit mass,
     while E is the real energy.*/
    REAL8 E = nu*H;
    const REAL8 ooH = 1./E;

    const REAL8 sqrtAbyB       = sqrt(A/B);
    const REAL8 dp_rstar_dt_0  = - sqrtAbyB*dHeff_dr*ooH;
    const REAL8 ddotr_dp_rstar = sqrtAbyB*d2Heff_dprstar20*ooH;
    const REAL8 Omg            = dHeff_dpphi*ooH;
    const REAL8 ddotr          = dp_rstar_dt_0*ddotr_dp_rstar; /* approximate ddot(r)_0 without Fphi, order pr_star^2 neglected */

    /* r evol eqn rhs */
    dy[TEOB_EVOLVE_RAD] = sqrtAbyB*dHeff_dprstar*ooH;

    /* phi evol eqn rhs */
    dy[TEOB_EVOLVE_PHI] = Omg;

    /* dp_{r*}/dt */
    dy[TEOB_EVOLVE_PRSTAR] = -sqrtAbyB*dHeff_dr*ooH;

    /* Compute here the new r_omg radius
     Compute same quantities with prstar=0. This to obtain psi.
     Procedure consistent with the nonspinning case. */
    REAL8 ggm0[14];
    eob_dyn_s_GS(r, rc, drc_dr, aK2, 0., pphi, nu, chi1, chi2, X1, X2, c3, ggm0);

    const REAL8 GS_0       = ggm0[2];
    const REAL8 GSs_0      = ggm0[3];
    const REAL8 dGS_dr_0   = ggm0[6];
    const REAL8 dGSs_dr_0  = ggm0[7];
    const REAL8 Heff_orb_0 = sqrt(A*(1.0 + pphi2*uc2));    /* effective Hamiltonian H_0^eff */
    const REAL8 Heff_0     = Heff_orb_0 + (GS_0*S + GSs_0*Sstar)*pphi;
    const REAL8 H0         = sqrt(1.0 + 2.0*nu*(Heff_0 - 1.0) );
    const REAL8 ooH0       = 1./H0;
    const REAL8 Gtilde     = GS_0*S     + GSs_0*Sstar;
    const REAL8 dGtilde_dr = dGS_dr_0*S + dGSs_dr_0*Sstar;
    const REAL8 duc_dr     = -uc2*drc_dr;
    const REAL8 psic       = (duc_dr + dGtilde_dr*rc*sqrt(A/pphi2 + A*uc2)/A)/(-0.5*dA);
    const REAL8 r_omg      = pow( ((1./sqrt(rc*rc*rc*psic))+Gtilde)*ooH0, -2./3. );
    const REAL8 v_phi      = r_omg*Omg;
    const REAL8 x          = v_phi*v_phi;
    const REAL8 jhat       = pphi/(r_omg*v_phi);

    /* Compute flux and \f$dp_{\phi}/dt\f$ */
    if (dyn->noflux) dy[TEOB_EVOLVE_PPHI] = 0.;
    else            dy[TEOB_EVOLVE_PPHI] = eob_flx_Flux_s(x,Omg,r_omg,E,Heff,jhat,r,prstar,ddotr,dyn);

    if (dyn->store) {
        /* Store values */
        dyn->t = t;
        dyn->r = r;
        dyn->phi = y[TEOB_EVOLVE_PHI];
        dyn->pphi = pphi;
        dyn->prstar = prstar;
        dyn->Omg = Omg;
        dyn->Omg_orb = ooH*pphi*A*uc2/Heff_orb;
        dyn->H = H;
        dyn->E = E;
        dyn->Heff = Heff;
        dyn->A = A;
        dyn->dA = dA;
        dyn->d2A = d2A;
        dyn->B = B;
        dyn->dB = dB;
        //    dyn->psi = psi;
        dyn->r_omega = r_omg;
        dyn->v_phi = v_phi;
        dyn->jhat = jhat;
        dyn->ddotr = ddotr;
    }

    return GSL_SUCCESS;
}


/* Computes the gyro-gravitomagnetic functions GS and GS*, that are called GS and GSs.
 r      => BL radius
 aK2    => squared Kerr parameter
 prstar => r* conjugate momentum
 nu     => symmetric mass ratio
 the CN3LO parameter is hard-coded in this routine
 ggm is the output structure. */

void eob_dyn_s_GS(REAL8 r,
                  REAL8 rc,
                  REAL8 drc_dr,
                  REAL8 UNUSED aK2,
                  REAL8 prstar,
                  REAL8 UNUSED pph,
                  REAL8 nu,
                  REAL8 UNUSED chi1,
                  REAL8 UNUSED chi2,
                  REAL8 UNUSED X1,
                  REAL8 UNUSED X2,
                  REAL8 cN3LO,
                  REAL8 *ggm)
{
//    static REAL8 c10,c20,c30,c02,c12,c04;
//    static REAL8 cs10,cs20,cs30,cs40,cs02,cs12,cs04;
    REAL8 c10,c20,c30,c02,c12,c04;
    REAL8 cs10,cs20,cs30,cs40,cs02,cs12,cs04;

    /* Compute the nu-dep. coefficient at first call only */
    // MA: FIXME: is this dangerous for lalinference?
//    static int firstcall = 1;
//    if (firstcall) {
//        firstcall = 0;
        REAL8 nu2   = nu*nu;
        /* coefficients of hat{GS} */
        c10 =  5./16.*nu;
        c20 =  51./8.*nu + 41./256.*nu2;
        c30 =  nu*cN3LO;
        c02 =  27./16.*nu;
        c12 =  12.*nu - 49./128.*nu2;
        c04 = -5./16.*nu + 169./256.*nu2;
        /* coefficients of hat{GS*} */
        cs10 = 3./4.   + nu/2.;
        cs20 = 27./16. + 29./4.*nu + 3./8.*nu2;
        cs02 = 5./4.   + 3./2.*nu;
        cs12 = 4.   + 11.*nu     - 7./8.*nu2;
        cs04 = 5./48.  + 25./12.*nu + 3./8.*nu2;
        cs30 = nu*cN3LO + 135./32.;
        cs40 = 2835./256.;
//    }

    REAL8 u   = 1./r;
    REAL8 u2  = u*u;

    REAL8 uc      = 1./rc;
    REAL8 uc2     = uc*uc;
    REAL8 uc3     = uc2*uc;
    REAL8 uc4     = uc3*uc;
    REAL8 prstar2 = prstar*prstar;
    REAL8 prstar4 = prstar2*prstar2;

    REAL8 GS0       = 2.*u*uc2;
    REAL8 dGS0_duc  = 2.*u2/drc_dr + 4.*u*uc;

    REAL8 GSs0          =  3./2.*uc3;
    REAL8 dGSs0_duc     =  9./2.*uc2;
    REAL8 dGSs0_dprstar =  0.0;
    REAL8 dGSs0_dpph    =  0.0;

    REAL8 hGS  = 1./(1.  + c10*uc + c20*uc2 + c30*uc3 + c02*prstar2 + c12*uc*prstar2 + c04*prstar4);
    REAL8 hGSs = 1./(1.  + cs10*uc + cs20*uc2  + cs30*uc3 + cs40*uc4 + cs02*prstar2 + cs12*uc*prstar2 + cs04*prstar4);

    /* complete gyro-gravitomagnetic functions */
    REAL8 GS  =  GS0*hGS;
    REAL8 GSs = GSs0*hGSs;

    /* Get derivatives of gyro-gravitomagnetic functions */
    REAL8 dhGS_dprstar  = -2.*prstar*hGS*hGS *( c02 +  c12*uc +  2.*c04*prstar2);
    REAL8 dhGSs_dprstar = -2.*prstar*hGSs*hGSs*(cs02 + cs12*uc + 2.*cs04*prstar2);

    REAL8 dGS_dprstar  = GS0 *dhGS_dprstar;
    REAL8 dGSs_dprstar = GSs0*dhGSs_dprstar + dGSs0_dprstar*hGSs;

    /* derivatives of hat{G} with respect to uc */
    REAL8 dhGS_duc  = -hGS*hGS*(c10 + 2.*c20*uc  + 3.*c30*uc2 + c12*prstar2);
    REAL8 dhGSs_duc = -hGSs*hGSs*(cs10 + 2.*cs20*uc + 3.*cs30*uc2 + 4.*cs40*uc3 + cs12*prstar2);

    /* derivatives of G with respect to uc */
    REAL8 dGS_duc  =  dGS0_duc*hGS  +  GS0*dhGS_duc;
    REAL8 dGSs_duc = dGSs0_duc*hGSs + GSs0*dhGSs_duc;

    /* derivatives of (G,G*) with respect to r */
    REAL8 dGS_dr  = -drc_dr*uc2*dGS_duc;
    REAL8 dGSs_dr = -drc_dr*uc2*dGSs_duc;

    /* derivatives of (G,G*) with respect to pph */
    REAL8 dGS_dpph  = 0.;
    REAL8 dGSs_dpph = dGSs0_dpph*hGSs;

    /* For initial data: compute the two ratios of ggm.dG_dprstar/prstar for GS and GSs */
    const REAL8 dGS_dprstarbyprstar  = -2.*GS0*hGS*hGS *( c02  +  c12*uc +  2.*c04*prstar2);
    const REAL8 dGSs_dprstarbyprstar = -2.*GSs0*hGSs*hGSs*(cs02 + cs12*uc + 2.*cs04*prstar2);

    /* For NQC: Second derivatives neglecting all pr_star^2 terms */
    const REAL8 d2GS_dprstar20  =  GS0*(-2.*hGS*hGS *( c02 +  c12*uc +  2.*c04*prstar2));
    const REAL8 d2GSs_dprstar20 =  GSs0*(-2.*hGSs*hGSs*(cs02 + cs12*uc + 2.*cs04*prstar2));

    ggm[0]=hGS;
    ggm[1]=hGSs;
    ggm[2]=GS;
    ggm[3]=GSs;
    ggm[4]=dGS_dprstar;
    ggm[5]=dGSs_dprstar;
    ggm[6]=dGS_dr;
    ggm[7]=dGSs_dr;
    ggm[8]=dGS_dpph;
    ggm[9]=dGSs_dpph;
    ggm[10]=dGS_dprstarbyprstar;
    ggm[11]=dGSs_dprstarbyprstar;
    ggm[12]=d2GS_dprstar20;
    ggm[13]=d2GSs_dprstar20;

    return;
}


/* Define radius rc that includes of LO spin-square coupling.  */
/*
 The S1*S2 term coincides with the BBH one, no effect of structure.
 The self-spin couplings, S1*S1 and S2*S2 get a EOS-dependent coefficient, CQ, that describe the quadrupole
 deformation due to spin. Notation of Levi-Steinhoff, JCAP 1412 (2014), no.12, 003. Notation analogous to
 the parameter a of Poisson, PRD 57, (1998) 5287-5290 or C_ES^2 in Porto & Rothstein, PRD 78 (2008), 044013

 The implementation uses the I-Love-Q fits of Table I of Yunes-Yagi
 paper, PRD 88, 023009, the \f$bar{Q}(bar{\lambda)^{tid})\f$ relation, line 3 of the table.
 The dimensionless \f$bar{\lambda}\f$ love number is related to our apsidal constant as \f$lambda = 2/3 k2/(C^5)\f$ so that both quantities have to appear here.
 */
void eob_dyn_s_get_rc_LO(REAL8 r,
                         REAL8 nu,
                         REAL8 at1,
                         REAL8 at2,
                         REAL8 aK2,
                         REAL8 C_Q1,
                         REAL8 C_Q2,
                         REAL8 UNUSED C_Oct1,
                         REAL8 UNUSED C_Oct2,
                         REAL8 UNUSED C_Hex1,
                         REAL8 UNUSED C_Hex2,
                         int usetidal,
                         REAL8 *rc,
                         REAL8 *drc_dr,
                         REAL8 *d2rc_dr2)
{

    REAL8 u   = 1./r;
    REAL8 u2  = u*u;
    REAL8 u3  = u*u2;
    REAL8 r2  = r*r;

    if (usetidal) {
#if (RC_EXCLUDESPINSPINTIDES)
        /* Switch off spin-spin-tidal couplings */
        /* See also: eob_wav_flm_s() */
        REAL8 rc2 = r2;
        *rc = r;
        *drc_dr = 1;
        *d2rc_dr2 = 0;
        /* Above code switch off everything,
         Alt. one can set C_Q1=C_Q2=0, but keep centrifugal radius */
        /*
         REAL8 a02  = 2.*at1*at2;
         REAL8 rc2  = r2 + a02*(1.+2.*u);
         *rc         = sqrt(rc2);
         *drc_dr     = r/(*rc)*(1.-a02*u3);
         *d2rc_dr2   = 1./(*rc)*(1.-(*drc_dr)*r/(*rc)*(1.-a02*u3)+2.*a02*u3);
         */
#else
        /* BNS effective spin parameter */
        REAL8 a02  = C_Q1*at1*at1 + 2.*at1*at2 + C_Q2*at2*at2;
        REAL8 rc2  = r2 + a02*(1.+2.*u); /* tidally-modified centrifugal radius */
        *rc         = sqrt(rc2);
        *drc_dr     = r/(*rc)*(1.-a02*u3);
        *d2rc_dr2   = 1./(*rc)*(1.-(*drc_dr)*r/(*rc)*(1.-a02*u3)+2.*a02*u3);
#endif
    } else {
        /*
         REAL8 X12      = sqrt(1.-4.*nu);
         REAL8 alphanu2 = 1. + 0.5/aK2*(- at2*at2*(5./4. + 5./4.*X12 + nu/2.) - at1*at1*(5./4. - 5./4.*X12 +nu/2.) + at1*at2*(-2.+nu));
         REAL8 rc2 = r2 + aK2*(1. + 2.*alphanu2/r);
         *rc         = sqrt(rc2);
         *drc_dr     = r/(*rc)*(1.+aK2*(-alphanu2*u3 ));
         *d2rc_dr2   = 1./(*rc)*(1.-(*drc_dr)*r/(*rc)*(1.-alphanu2*aK2*u3)+ 2.*alphanu2*aK2*u3);
         */
        /* Following implementation is regular (avoids 1/aK2) */
        REAL8 X12 = sqrt(1.-4.*nu);
        REAL8 c_ss_nlo = (- at2*at2*(1.25 + 1.25*X12 + 0.5*nu) - at1*at1*(1.25 - 1.25*X12 + 0.5*nu) + at1*at2*(-2.+nu));
        REAL8 rc2   = r2 + aK2*(1. + 2.*u) + u*c_ss_nlo;
        *rc          = sqrt(rc2);
        REAL8 divrc = 1.0/(*rc);
        *drc_dr      = r*divrc*(1-(aK2 + 0.5*c_ss_nlo)*u3);
        *d2rc_dr2    = divrc*(1.-(*drc_dr)*r*divrc*(1.-(aK2+0.5*c_ss_nlo)*u3)+ (2.*aK2 + c_ss_nlo)*u3);
    }

    return;
}

/* tidal rc with NLO coefficient that depends on C_Qi */
void eob_dyn_s_get_rc_NLO(REAL8 r,
                          REAL8 nu,
                          REAL8 at1,
                          REAL8 at2,
                          REAL8 aK2,
                          REAL8 C_Q1,
                          REAL8 C_Q2,
                          REAL8 UNUSED C_Oct1,
                          REAL8 UNUSED C_Oct2,
                          REAL8 UNUSED C_Hex1,
                          REAL8 UNUSED C_Hex2,
                          int usetidal,
                          REAL8 *rc,
                          REAL8 *drc_dr,
                          REAL8 *d2rc_dr2)
{

    REAL8 u   = 1./r;
    REAL8 u2  = u*u;
    REAL8 u3  = u*u2;
    REAL8 r2  = r*r;
    REAL8 X12 = sqrt(1.-4.*nu);

    if (usetidal) {

        /* BNS effective spin parameter */
        REAL8 a02      = C_Q1*at1*at1 + 2.*at1*at2 + C_Q2*at2*at2;

        REAL8 delta_a2 = X12*(at1*at1*(C_Q1+0.25) - at2*at2*(C_Q2+0.25))
        + at1*at1*(-17./4.+3.*C_Q1-0.5*nu)
        + at2*at2*(-17./4.+3.*C_Q2-0.5*nu)
        + at1*at2*(nu-2.0);

        REAL8 rc2 = r2 + a02*(1. + 2.*u) + delta_a2*u;
        *rc         = sqrt(rc2);
        REAL8 divrc = 1.0/(*rc);
        *drc_dr     = divrc*(r - (a02 + 0.5*delta_a2)*u2);
        *d2rc_dr2   = divrc*(1 + (2.*a02 + delta_a2)*u3 - (*drc_dr)*(*drc_dr));

    }
    else {

        REAL8 c_ss_nlo = (- at2*at2*(1.25 + 1.25*X12 + 0.5*nu) - at1*at1*(1.25 - 1.25*X12 + 0.5*nu) + at1*at2*(-2.+nu));
        REAL8 rc2   = r2 + aK2*(1. + 2.*u) + u*c_ss_nlo;
        *rc          = sqrt(rc2);
        REAL8 divrc = 1.0/(*rc);
        *drc_dr      = r*divrc*(1-(aK2 + 0.5*c_ss_nlo)*u3);
        *d2rc_dr2    = divrc*(1.-(*drc_dr)*r*divrc*(1.-(aK2+0.5*c_ss_nlo)*u3)+ (2.*aK2 + c_ss_nlo)*u3);

    }

}

/*  tidal rc with NNLO coefficient that depends on C_Qi */
void eob_dyn_s_get_rc_NNLO(REAL8 r,
                           REAL8 nu,
                           REAL8 at1,
                           REAL8 at2,
                           REAL8 aK2,
                           REAL8 C_Q1,
                           REAL8 C_Q2,
                           REAL8 UNUSED C_Oct1,
                           REAL8 UNUSED C_Oct2,
                           REAL8 UNUSED C_Hex1,
                           REAL8 UNUSED C_Hex2,
                           int usetidal,
                           REAL8 *rc,
                           REAL8 *drc_dr,
                           REAL8 *d2rc_dr2)
{

    REAL8 u   = 1./r;
    REAL8 u2  = u*u;
    REAL8 u3  = u*u2;
    REAL8 u4  = u*u3;
    REAL8 u5  = u*u4;
    REAL8 r2  = r*r;
    REAL8 X12 = sqrt(1.-4.*nu);

    if (usetidal)
    {

        /* BNS effective spin parameter */
        REAL8 a02      = C_Q1*at1*at1 + 2.*at1*at2 + C_Q2*at2*at2;

        REAL8 delta_a2 = X12*(at1*at1*(C_Q1+0.25) - at2*at2*(C_Q2+0.25))
        + at1*at1*(-17./4.+3.*C_Q1-0.5*nu)
        + at2*at2*(-17./4.+3.*C_Q2-0.5*nu)
        + at1*at2*(nu-2.0);

        REAL8 delta_a2_nnlo  =
        (  387./28.  - 207./28.*nu              )     *a02
        + (-2171./212. - 269./28.*nu + 0.375*nu*nu)     *(at1*at1+at2*at2)
        + (- 281./7    - 187./56.*nu - 0.75 *nu*nu)     *at1*at2
        +    163./28.                               *X12*(C_Q1*at1*at1-C_Q2*at2*at2)
        + (  -29./112. - 2.625   *nu              ) *X12*(at1*at1-at2*at2);

        REAL8 rc2   =  r2 + a02*(1. + 2.*u) + delta_a2*u + delta_a2_nnlo*u2;
        *rc          = sqrt(rc2);
        REAL8 divrc = 1.0/(*rc);
        *drc_dr      = divrc*(r - (a02 + 0.5*delta_a2)*u2 - delta_a2_nnlo*u3);
        *d2rc_dr2    = divrc*(1 + (2.*a02 + delta_a2)*u3
                              + 3*delta_a2_nnlo*u4 - (*drc_dr)*(*drc_dr));

    }
    else
    {

        REAL8 a0  = at1 + at2;
        REAL8 a12 = at1 - at2;

        REAL8 c_ss_nlo = -1.125*a0*a0 -(0.125+0.5+nu)*a12*a12 + 1.25*X12*a0*a12;

        REAL8 c_ss_nnlo = - (189./32. + 417.32*nu) *a0 *a0
        + ( 11./32. - 127.32*nu + 0.375*nu*nu) *a12*a12
        + ( 87.16   -  2.625*nu )*X12 *a0 *a12;


        REAL8 rc2   = r2 + aK2*(1. + 2.*u) + u*c_ss_nlo + u2*c_ss_nnlo;
        *rc          = sqrt(rc2);
        REAL8 divrc = 1.0/(*rc);
        *drc_dr      = r*divrc*(1-(aK2 + 0.5*c_ss_nlo)*u3 - 0.5*u4*c_ss_nnlo);
        *d2rc_dr2    = 1./r*(*drc_dr) + r*divrc*((3.*aK2+c_ss_nlo)*u4 + 2.*c_ss_nnlo*u5);

    }

}

/*  tidal rc @ NNLO with the addition of the LO spin^4 coefficient that depends on C_Q, C_Oct and C_Hex */
void eob_dyn_s_get_rc_NNLO_S4(REAL8 r,
                              REAL8 nu,
                              REAL8 at1,
                              REAL8 at2,
                              REAL8 aK2,
                              REAL8 C_Q1,
                              REAL8 C_Q2,
                              REAL8 C_Oct1,
                              REAL8 C_Oct2,
                              REAL8 C_Hex1,
                              REAL8 C_Hex2,
                              int usetidal,
                              REAL8 *rc,
                              REAL8 *drc_dr,
                              REAL8 *d2rc_dr2)
{

    REAL8 u   = 1./r;
    REAL8 u2  = u*u;
    REAL8 u3  = u*u2;
    REAL8 u4  = u*u3;
    REAL8 u5  = u*u4;
    REAL8 r2  = r*r;
    REAL8 X12 = sqrt(1.-4.*nu);

    if (usetidal)
    {

        /* BNS effective spin parameter */
        REAL8 a02      = C_Q1*at1*at1 + 2.*at1*at2 + C_Q2*at2*at2;

        REAL8 delta_a2 = X12*(at1*at1*(C_Q1+0.25) - at2*at2*(C_Q2+0.25))
        + at1*at1*(-17./4.+3.*C_Q1-0.5*nu)
        + at2*at2*(-17./4.+3.*C_Q2-0.5*nu)
        + at1*at2*(nu-2.0);

        REAL8 delta_a2_nnlo  =
        (  387./28.  - 207./28.*nu              )     *a02
        + (-2171./212. - 269./28.*nu + 0.375*nu*nu)     *(at1*at1+at2*at2)
        + (- 281./7    - 187./56.*nu - 0.75 *nu*nu)     *at1*at2
        +    163./28.                               *X12*(C_Q1*at1*at1-C_Q2*at2*at2)
        + (  -29./112. - 2.625   *nu              ) *X12*(at1*at1-at2*at2);

        REAL8 delta_a4_lo = 0.75*(C_Hex1 - C_Q1*C_Q1)*at1*at1*at1*at1
        + 3.*(C_Oct1 - C_Q1)     *at1*at1*at1*at2
        + 3.*(C_Q1*C_Q2 - 1)     *at1*at1*at2*at2
        + 3.*(C_Oct2 - C_Q2)     *at1*at2*at2*at2
        + 0.75*(C_Hex2 - C_Q2*C_Q2)*at2*at2*at2*at2;

        REAL8 rc2   =  r2 + a02*(1. + 2.*u) + delta_a2*u + (delta_a2_nnlo+delta_a4_lo)*u2;
        *rc          = sqrt(rc2);
        REAL8 divrc = 1.0/(*rc);
        *drc_dr      = divrc*(r - (a02 + 0.5*delta_a2)*u2 - (delta_a2_nnlo+delta_a4_lo)*u3);
        *d2rc_dr2    = divrc*(1 + (2.*a02 + delta_a2)*u3
                              + 3*(delta_a2_nnlo+delta_a4_lo)*u4 - (*drc_dr)*(*drc_dr));

    }
    else
    {

        REAL8 a0  = at1 + at2;
        REAL8 a12 = at1 - at2;

        REAL8 c_ss_nlo = -1.125*a0*a0 -(0.125+0.5+nu)*a12*a12 + 1.25*X12*a0*a12;

        REAL8 c_ss_nnlo = - (189./32. + 417.32*nu              )    *a0 *a0
        + ( 11./32. - 127.32*nu + 0.375*nu*nu)    *a12*a12
        + ( 87.16   -  2.625*nu              )*X12*a0 *a12;


        REAL8 rc2   = r2 + aK2*(1. + 2.*u) + u*c_ss_nlo + u2*c_ss_nnlo;
        *rc          = sqrt(rc2);
        REAL8 divrc = 1.0/(*rc);
        *drc_dr      = r*divrc*(1-(aK2 + 0.5*c_ss_nlo)*u3 - 0.5*u4*c_ss_nnlo);
        *d2rc_dr2    = 1./r*(*drc_dr) + r*divrc*((3.*aK2+c_ss_nlo)*u4 + 2.*c_ss_nnlo*u5);

    }

}

/* Non-spinning case -- rc = r */
void eob_dyn_s_get_rc_NOSPIN(REAL8 r,
                             REAL8 UNUSED nu,
                             REAL8 UNUSED at1,
                             REAL8 UNUSED at2,
                             REAL8 UNUSED aK2,
                             REAL8 UNUSED C_Q1,
                             REAL8 UNUSED C_Q2,
                             REAL8 UNUSED C_Oct1,
                             REAL8 UNUSED C_Oct2,
                             REAL8 UNUSED C_Hex1,
                             REAL8 UNUSED C_Hex2,
                             int UNUSED usetidal,
                             REAL8 UNUSED *rc,
                             REAL8 UNUSED *drc_dr,
                             REAL8 UNUSED *d2rc_dr2)
{
    *rc = r;
    *drc_dr = 1;
    *d2rc_dr2 = 0;
}

/* LO case with C_Q1 = 0 for tidal part*/
void eob_dyn_s_get_rc_NOTIDES(REAL8 r,
                              REAL8 nu,
                              REAL8 at1,
                              REAL8 at2,
                              REAL8 aK2,
                              REAL8 UNUSED C_Q1,
                              REAL8 UNUSED C_Q2,
                              REAL8 UNUSED C_Oct1,
                              REAL8 UNUSED C_Oct2,
                              REAL8 UNUSED C_Hex1,
                              REAL8 UNUSED C_Hex2,
                              int usetidal,
                              REAL8 *rc,
                              REAL8 *drc_dr,
                              REAL8 *d2rc_dr2)
{

    REAL8 u   = 1./r;
    REAL8 u2  = u*u;
    REAL8 u3  = u*u2;
    REAL8 r2  = r*r;

    if (usetidal) {
        /*  We set C_Q1=C_Q2=0, but keep centrifugal radius */

        REAL8 a02  = 2.*at1*at2;
        REAL8 rc2  = r2 + a02*(1.+2.*u);
        *rc         = sqrt(rc2);
        *drc_dr     = r/(*rc)*(1.-a02*u3);
        *d2rc_dr2   = 1./(*rc)*(1.-(*drc_dr)*r/(*rc)*(1.-a02*u3)+2.*a02*u3);

    } else {

        REAL8 X12 = sqrt(1.-4.*nu);
        REAL8 c_ss_nlo = (- at2*at2*(1.25 + 1.25*X12 + 0.5*nu) - at1*at1*(1.25 - 1.25*X12 + 0.5*nu) + at1*at2*(-2.+nu));
        REAL8 rc2   = r2 + aK2*(1. + 2.*u) + u*c_ss_nlo;
        *rc          = sqrt(rc2);
        REAL8 divrc = 1.0/(*rc);
        *drc_dr      = r*divrc*(1-(aK2 + 0.5*c_ss_nlo)*u3);
        *d2rc_dr2    = divrc*(1.-(*drc_dr)*r*divrc*(1.-(aK2+0.5*c_ss_nlo)*u3)+ (2.*aK2 + c_ss_nlo)*u3);
    }

}


/* EOB Metric A function 5PN log
 This function computes the Pade' (1,5) resummed A function (with its
 derivatives) starting from the 5PN-expanded version of the A function
 including 4PN and 5PN log terms.
 This represents the current, stable, most accurate implementation of
 the EOB effective potential

 Coefficients a5 and a6 are the nonlog contributions to the 4PN and 5PN terms.
 In practice, a5 is fixed to its GSF value computed in Akcay et al,

 \f$a5 \equiv a5_GSF = +23.50190(5) \approx +23.5\f$

 and \f$a6 \equiv a6(nu) = (-110.5 + 129*(1-4*nu)).*(1-1.5e-5/((0.26-nu)^2)\f$
 as obtained from comparison with the Caltech-Cornell-CITA numerical data.
 These values are used as default. */
void eob_metric_A5PNlog(REAL8 r, REAL8 nu, REAL8 *A, REAL8 *dA, REAL8 *d2A)
{

    /* shortcuts */
    REAL8 nu2 = nu*nu;
    REAL8 pi2 = SQ(LAL_PI);
    REAL8 pi4 = pi2*pi2;
    REAL8 u    = 1./r;
    REAL8 u2   = u*u;
    REAL8 u3   = u*u2;
    REAL8 u4   = u2*u2;
    REAL8 u5   = u4*u;
    REAL8 u6   = u5*u;
    REAL8 UNUSED u7   = u6*u;
    REAL8 UNUSED u10  = u5*u5;
    REAL8 u8   = u5*u3;
    REAL8 UNUSED u9   = u8*u;
    REAL8 logu = log(u);

    REAL8 a5c0 = -4237./60. + 2275./512.*pi2 + 256./5.*Log2 + 128./5.*LAL_GAMMA;
    REAL8 a5c1 = -221./6.   + 41./32.*pi2;
    REAL8 a5   =  a5c0 + nu*a5c1;
    REAL8 a6   =  3097.3*nu2 - 1330.6*nu + 81.38;

    /* 4PN and 5PN coefficients including all known log terms */
    REAL8 a5tot  = a5  + 64./5.*logu;
    REAL8 a6tot  = a6  + (-7004./105. - 144./5.*nu)*logu;
    REAL8 a5tot2 = a5tot*a5tot;

    /* Coefficients of the Padeed function */
    REAL8 N1 = (-3*(-512 - 32*nu2 + nu*(3520 + 32*a5tot + 8*a6tot - 123*pi2)))/(-768 + nu*(3584 + 24*a5tot - 123*pi2));
    REAL8 D1 = (nu*(-3392 - 48*a5tot - 24*a6tot + 96*nu + 123*pi2))/(-768 + nu*(3584 + 24*a5tot - 123*pi2));
    REAL8 D2 = (2*nu*(-3392 - 48*a5tot - 24*a6tot + 96*nu + 123*pi2))/(-768 + nu*(3584 + 24*a5tot - 123*pi2));
    REAL8 D3 = (-2*nu*(6016 + 48*a6tot + 3392*nu + 24*a5tot*(4 + nu) - 246*pi2 - 123*nu*pi2))/(-768 + nu*(3584 + 24*a5tot - 123*pi2));
    REAL8 D4 = -(nu*(-4608*a6tot*(-4 + nu) + a5tot*(36864 + nu*(72192 - 2952*pi2)) + nu*(2048*(5582 + 9*nu) - 834432*pi2 + 15129*pi4)))/(96.*(-768 + nu*(3584 + 24*a5tot - 123*pi2)));
    REAL8 D5 = (nu*(-24*a6tot*(1536 + nu*(-3776 + 123*pi2)) + nu*(-2304*a5tot2 + 96*a5tot*(-3392 + 123*pi2) - (-3776 + 123*pi2)*(-3008 - 96*nu + 123*pi2))))/(96.*(-768 + nu*(3584 + 24*a5tot - 123*pi2)));

    /* First derivatives */
    REAL8 dN1 = (160*nu*(-828672 - 32256*nu2 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*pi2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 174045*pi2)))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),2)*u);
    REAL8 dD1 = (160*nu*(-828672 - 32256*nu2 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*pi2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 174045*pi2)))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),2)*u);
    REAL8 dD2 = (320*nu*(-828672 - 32256*nu2 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*pi2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 174045*pi2)))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),2)*u);
    REAL8 dD3 = (640*nu*(-828672 - 32256*nu2 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*pi2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 174045*pi2)))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),2)*u);
    REAL8 dD4 = (-320*(-4 + nu)*nu*(-828672 - 32256*nu2 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*pi2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 174045*pi2)))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),2)*u);
    REAL8 dD5 = (nu*(-8400*nu*(-24*(a6 - (4*logu*(1751 + 756*nu))/105.)*(1536 + nu*(-3776 + 123*pi2)) + nu*(-2304*gsl_pow_int(a5 + (64*logu)/5.,2) + 96*(a5 + (64*logu)/5.)*(-3392 + 123*pi2) - (-3776 + 123*pi2)*(-32*(94 + 3*nu) + 123*pi2))) - (1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)))*(4128768*logu*nu + 5*(-2689536 + nu*(11170624 + 64512*a5 - 380685*pi2) - 756*nu*(1536 + nu*(-3776 + 123*pi2))))))/(2625.*gsl_pow_int(-768 + nu*(3584 + 24*(a5 + (64*logu)/5.) - 123*pi2),2)*u);

    /* Numerator and denominator of the Pade */
    REAL8 Num = 1 + N1*u;
    REAL8 Den = 1 + D1*u + D2*u2 + D3*u3 + D4*u4 + D5*u5;
    *A = Num/Den;

    /* First derivative */
    REAL8 dNum  = dN1*u + N1;
    REAL8 dDen  = D1 + u*(dD1 + 2*D2) + u2*(dD2 + 3*D3) + u3*(dD3 + 4*D4) + u4*(dD4 + 5*D5) + dD5*u5;

    /* Derivative of A function with respect to u */
    REAL8 prefactor = (*A)/(Num*Den);
    REAL8 dA_u      = prefactor*(dNum*Den - dDen*Num);

    /* Derivative of A with respect to r */
    /* *dA = -u2*dA_u; */

    *dA = dA_u;

    if (d2A != NULL) {

        /* Second derivatives of Pade coefficients */
        REAL8 d2N1 = (160*nu*(-3840 + 1536*logu*nu + nu*(20992 + 120*a5 - 615*pi2))*(828672 + nu*(-42024*a5 - 8064*a6 + 3584*(-1397 + 9*nu) + 174045*pi2) + 756*nu*(768 + nu*(-3584 - 24*a5 + 123*pi2))))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),3)*u2);
        REAL8 d2D1 = (160*nu*(-3840 + 1536*logu*nu + nu*(20992 + 120*a5 - 615*pi2))*(828672 + nu*(-42024*a5 - 8064*a6 + 3584*(-1397 + 9*nu) + 174045*pi2) + 756*nu*(768 + nu*(-3584 - 24*a5 + 123*pi2))))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),3)*u2);
        REAL8 d2D2 = (320*nu*(-3840 + 1536*logu*nu + nu*(20992 + 120*a5 - 615*pi2))*(828672 + nu*(-42024*a5 - 8064*a6 + 3584*(-1397 + 9*nu) + 174045*pi2) + 756*nu*(768 + nu*(-3584 - 24*a5 + 123*pi2))))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),3)*u2);
        REAL8 d2D3 = (640*nu*(-3840 + 1536*logu*nu + nu*(20992 + 120*a5 - 615*pi2))*(828672 + nu*(-42024*a5 - 8064*a6 + 3584*(-1397 + 9*nu) + 174045*pi2) + 756*nu*(768 + nu*(-3584 - 24*a5 + 123*pi2))))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),3)*u2);
        REAL8 d2D4 = (320*(-4 + nu)*nu*(-828672 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*pi2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 32256*nu - 174045*pi2))*(-3840 + 1536*logu*nu + nu*(20992 + 120*a5 - 615*pi2)))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),3)*u2);
        REAL8 d2D5 = (nu*(gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),2)*(4128768*logu*nu - 7680*(1751 + 756*nu) + nu*(64*(808193 + 5040*a5 + 223020*nu) - 615*(3095 + 756*nu)*pi2)) + 3072*nu*(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)))*(4128768*logu*nu - 7680*(1751 + 756*nu) + 5*nu*(64*(174541 + 1008*a5 + 44604*nu) - 123*(3095 + 756*nu)*pi2)) + 25804800*nu2*(-24*(a6 - (4*logu*(1751 + 756*nu))/105.)*(1536 + nu*(-3776 + 123*pi2)) + nu*(-2304*gsl_pow_int(a5 + (64*logu)/5.,2) + 96*(a5 + (64*logu)/5.)*(-3392 + 123*pi2) - (-3776 + 123*pi2)*(-32*(94 + 3*nu) + 123*pi2))) + 42000*nu*(-768 + nu*(3584 + 24*(a5 + (64*logu)/5.) - 123*pi2))*(-24*(a6 - (4*logu*(1751 + 756*nu))/105.)*(1536 + nu*(-3776 + 123*pi2)) + nu*(-2304*gsl_pow_int(a5 + (64*logu)/5.,2) + 96*(a5 + (64*logu)/5.)*(-3392 + 123*pi2) - (-3776 + 123*pi2)*(-32*(94 + 3*nu) + 123*pi2)))))/(13125.*gsl_pow_int(-768 + nu*(3584 + 24*(a5 + (64*logu)/5.) - 123*pi2),3)*u2);

        /* Second derivative of numerator and denominator */
        REAL8 d2Num = 2.*dN1 + d2N1*u;
        REAL8 d2Den = 2.*(D2 + dD1) + u*(6.*D3 + 4.*dD2 + d2D1) + u2*(12.*D4 + 6.*dD3 + d2D2) + u3*(20.*D5 + 8.*dD4 + d2D3) + u4*(10.*dD5 + d2D4) + u5*d2D5;

        /* Second derivative with respect of u */
        REAL8 d2A_u = prefactor*(2.*dDen*dDen*(*A) - 2.*dNum*dDen + Den*d2Num - d2Den*Num);

        *d2A = d2A_u;

        /* Second derivative with respect of r */
        /* *d2A = u4*d2A_u + 2.*u3*dA_u; */

    }
}

/* Tidal potential, two version implemented:
 1. TEOB NNLO, Bernuzzi+ 1205.3403
 2. TEOBResum, Bini&Damour, 1409.6933, Bernuzzi+ 1412.4553 */
void eob_metric_Atidal(REAL8 r,
                       LALTEOBResumSDynamics *dyn,
                       REAL8 *AT,
                       REAL8 *dAT,
                       REAL8 *d2AT)
{

    REAL8 A, dA_u, d2A_u, UNUSED dA, UNUSED d2A;
    A = dA_u = d2A_u = 0.;

    const REAL8 elsix = 1.833333333333333333333;  // 11/6
    const REAL8 eightthird = 2.6666666666666666667; // 8/3
    const REAL8 nu    = dyn->nu;
    const REAL8 rLR   = dyn->rLR_tidal;
    const REAL8 XA    = dyn->X1;
    const REAL8 XB    = dyn->X2;
    const REAL8 kapA2 = dyn->kapA2;
    const REAL8 kapB2 = dyn->kapB2;
    const REAL8 kapA3 = dyn->kapA3;
    const REAL8 kapB3 = dyn->kapB3;
    const REAL8 kapT2 = dyn->kapT2;
    const REAL8 kapT3 = dyn->kapT3;
    const REAL8 kapT4 = dyn->kapT4;
    const REAL8 kapA2j = dyn->kapA2j;
    const REAL8 kapB2j = dyn->kapB2j;
    const REAL8 kapT2j = dyn->kapT2j;

    /* Definition of the conservative tidal coefficients \f$\bar{\alpha}_n^{(\ell)}\f$,
     Eq.(37) of Damour&Nagar, PRD 81, 084016 (2010) */
    const REAL8 bar_alph2_1 = dyn->bar_alph2_1;
    const REAL8 bar_alph2_2 = dyn->bar_alph2_2;
    const REAL8 bar_alph3_1 = dyn->bar_alph3_1;
    const REAL8 bar_alph3_2 = dyn->bar_alph3_2;
    const REAL8 bar_alph2j_1 = dyn->bar_alph2j_1;


    const REAL8 p = dyn->pGSF_tidal;

    /* shortcuts */
    REAL8 UNUSED nu2  = nu*nu;
    REAL8 pi2  = LAL_PI*LAL_PI;
    REAL8 UNUSED pi4  = pi2*pi2;
    REAL8 u    = 1./r;
    REAL8 u2   = u*u;
    REAL8 u3   = u*u2;
    REAL8 u4   = u2*u2;
    REAL8 u5   = u4*u;
    REAL8 u6   = u5*u;
    REAL8 u7   = u6*u;
    REAL8 u10  = u5*u5;
    REAL8 u8   = u5*u3;
    REAL8 u9   = u8*u;
    REAL8 UNUSED logu = log(u);
    REAL8 oom3u  = 1./(1.-rLR*u);


    if (dyn->use_tidal==LAL_SIM_INSPIRAL_GETIDES_NNLO) {

        A    = -(kapT4*u10) - kapT2*u6*(1. + bar_alph2_1*u + bar_alph2_2*u2) - kapT3*u8*(1. + bar_alph3_1*u + bar_alph3_2*u2);
        dA_u = -10.*kapT4*u9 - kapT2*u6*(bar_alph2_1 + 2.*bar_alph2_2*u) - kapT3*u8*(bar_alph3_1 + 2.*bar_alph3_2*u)
        - 6.*kapT2*u5*(1. + bar_alph2_1*u + bar_alph2_2*u2) - 8.*kapT3*u7*(1. + bar_alph3_1*u + bar_alph3_2*u2);

        if (d2AT != NULL) {
            d2A_u = -90.*kapT4*u8
            - kapT2*(2*bar_alph2_2*u6 + 12.*u5*(bar_alph2_1 + 2*bar_alph2_2*u)
                     + 30.*u4*(1 + bar_alph2_1*u + bar_alph2_2*u2))
            - kapT3*(2.*bar_alph3_2*u8 + 16.*u7*(bar_alph3_1 + 2*bar_alph3_2*u) + 56.*u6*(1 + bar_alph3_1*u + bar_alph3_2*u2));
        }

    } else if (dyn->use_tidal==LAL_SIM_INSPIRAL_GETIDES_GSF2) {

        const REAL8 c1  =  8.533515908;      // OLD value 8.53353;
        const REAL8 c2  = 3.043093411;    // OLD value 3.04309;
        const REAL8 n1  =  0.8400636422;     // OLD value 0.840058;
        const REAL8 d2  =  17.7324036;    // OLD value 17.73239

        REAL8 Acub   = 5./2.* u * (1. -  (c1+c2)*u +   c1*c2*u2);
        REAL8 dAcub  = 5./2.*     (1. -2*(c1+c2)*u + 3*c1*c2*u2);
        REAL8 d2Acub = 5    *     (   -  (c1+c2)   + 3*c1*c2*u);
        REAL8 Den    = 1./(1. + d2*u2);
        REAL8 f23    = (1. + n1*u)*Den;
        REAL8 df23   = (n1 - 2*d2*u - n1*d2*u2)*(Den*Den);
        REAL8 A1SF   = Acub*f23;
        REAL8 dA1SF  = dAcub*f23 + Acub*df23;
        REAL8 A2SF   = 337./28.*u2;
        REAL8 dA2SF  = 674./28.*u;

        REAL8 f0     = 1 + 3*u2*oom3u;
        REAL8 f1     = A1SF *pow(oom3u,7./2.);
        REAL8 f2     = A2SF *pow(oom3u,p);

        REAL8 df0    = 3*u*(2.-rLR*u)*(oom3u*oom3u);
        REAL8 df1    = 0.5*(7*rLR*A1SF + 2*(1.-rLR*u)*dA1SF)*pow(oom3u,9./2.);
        REAL8 df2    = (rLR*p*A2SF + (1.-rLR*u)*dA2SF)*pow(oom3u,p+1);

        /* Gravito-electric tides for el = 2, 3, 4 */
        REAL8 AT2    = - kapA2*u6*( f0 + XA*f1 + XA*XA*f2 ) - kapB2*u6*( f0 + XB*f1 + XB*XB*f2 );
        REAL8 AT3    = - kapT3*u8*(1. + bar_alph3_1*u + bar_alph3_2*u2);
        REAL8 AT4    = - kapT4*u10;

        REAL8 dAT2  = - kapA2*6.*u5*( f0 + XA*f1 + XA*XA*f2 ) - kapB2*6.*u5*( f0 + XB*f1 + XB*XB*f2 ) - kapA2*u6*( df0 + XA*df1 + XA*XA*df2 ) - kapB2*u6*( df0 + XB*df1 + XB*XB*df2 );
        REAL8 dAT3  = - kapT3*(8.*u7 + 9*bar_alph3_1*u8 + 10*bar_alph3_2*u9);
        REAL8 dAT4  = - kapT4*10.*u9;

        A     = AT2   + AT3   + AT4;
        dA_u  = dAT2  + dAT3  + dAT4;

        if (d2AT != NULL) {
            REAL8 d2f23  = 2*d2*(-1 + 3*d2*u2 + n1*u*(-3+d2*u2))*(Den*Den*Den);
            REAL8 d2A1SF = d2Acub*f23 + 2*dAcub*df23 + Acub*d2f23;
            REAL8 d2A2SF = 674./28.;
            REAL8 d2f0   = 6*(oom3u*oom3u*oom3u);
            REAL8 d2f1   = 0.25*(63*(rLR*rLR)*A1SF + 4*(-1+rLR*u)*(-7*rLR*dA1SF + (-1+rLR*u)*d2A1SF))*pow(oom3u,11./2.);
            REAL8 d2f2   = (  rLR*p*(1+p)*rLR*A2SF +(-1+rLR*u)*( -2.*p*rLR*dA2SF +(-1.+rLR*u)*d2A2SF )  )*pow(oom3u,p+2);

            REAL8 d2AT2  = - kapA2*30*u4*( f0 + XA*f1 + XA*XA*f2 ) - kapB2*30*u4*( f0 + XB*f1 + XB*XB*f2 ) - 2*kapA2*6*u5*( df0 + XA*df1 + XA*XA*df2 ) - 2*kapB2*6*u5*( df0 + XB*df1 + XB*XB*df2 ) - kapA2*u6*( d2f0 + XA*d2f1 + XA*XA*d2f2 ) - kapB2*u6*( d2f0 + XB*d2f1 + XB*XB*d2f2 );
            REAL8 d2AT3  = - kapT3*(56*u6 + 72*bar_alph3_1*u7 + 90*bar_alph3_2*u8);
            REAL8 d2AT4  = - kapT4*90*u8;

            d2A_u = d2AT2 + d2AT3 + d2AT4;

        }
    } else if (dyn->use_tidal==LAL_SIM_INSPIRAL_GETIDES_GSF23) {

        const REAL8 c1  =  8.533515908;
        const REAL8 c2  = 3.043093411;
        const REAL8 n1  =  0.8400636422;
        const REAL8 d2  =  17.7324036;

        REAL8 Acub   = 5./2.* u * (1. -  (c1+c2)*u +   c1*c2*u2);
        REAL8 dAcub  = 5./2.*     (1. -2*(c1+c2)*u + 3*c1*c2*u2);
        REAL8 d2Acub = 5    *     (   -  (c1+c2)   + 3*c1*c2*u);
        REAL8 Den    = 1./(1. + d2*u2);
        REAL8 f23    = (1. + n1*u)*Den;
        REAL8 df23   = (n1 - 2*d2*u - n1*d2*u2)*(Den*Den);
        REAL8 A1SF   = Acub*f23;
        REAL8 dA1SF  = dAcub*f23 + Acub*df23;
        REAL8 A2SF   = 337./28.*u2;
        REAL8 dA2SF  = 674./28.*u;

        REAL8 f0     = 1 + 3*u2*oom3u;
        REAL8 f1     = A1SF *pow(oom3u,7./2.);
        REAL8 f2     = A2SF *pow(oom3u,p);

        REAL8 df0    = 3*u*(2.-rLR*u)*(oom3u*oom3u);
        REAL8 df1    = 0.5*(7*rLR*A1SF + 2*(1.-rLR*u)*dA1SF)*pow(oom3u,9./2.);
        REAL8 df2    = (rLR*p*A2SF + (1.-rLR*u)*dA2SF)*pow(oom3u,p+1);

        /* Gravito-electric tides for el = 2, 4; el = 3 added below as a GSF series */
        REAL8 AT2    = - kapA2*u6*( f0 + XA*f1 + XA*XA*f2 ) - kapB2*u6*( f0 + XB*f1 + XB*XB*f2 );
        REAL8 AT4    = - kapT4*u10;

        REAL8 dAT2  = - kapA2*6.*u5*( f0 + XA*f1 + XA*XA*f2 ) - kapB2*6.*u5*( f0 + XB*f1 + XB*XB*f2 ) - kapA2*u6*( df0 + XA*df1 + XA*XA*df2 ) - kapB2*u6*( df0 + XB*df1 + XB*XB*df2 );
        REAL8 dAT4  = - kapT4*10.*u9;

        /* el = 3+, i.e.,  even parity tidal potential **/

        /* 1GSF fitting parameters */
        const REAL8 C1 = -3.6820949997216643;
        const REAL8 C2 = 5.171003322924513;
        const REAL8 C3 = -7.639164165720986;
        const REAL8 C4 = -8.63278143009751;
        const REAL8 C5 = 12.319646912775516;
        const REAL8 C6 = 16.36009385150114;

        /* 0SF -- el = 3+, i.e.,  even parity terms */
        REAL8 A3hat_Sch    = (1.0 - 2.0*u)*( 1.0 + eightthird*u2*oom3u );
        REAL8 dA3hat_Sch   = (1.0 - 2.0*u)*( eightthird*rLR*u2*oom3u*oom3u + 2.0*eightthird*u*oom3u ) - 2.0*( 1.0 + eightthird*u2*oom3u );
        REAL8 d2A3hat_Sch  = (1.0 - 2.0*u)*( 2.0*eightthird*rLR*rLR*u2*oom3u*oom3u*oom3u + 4.0*eightthird*rLR*u*oom3u*oom3u + 2.0*eightthird*oom3u ) - 4.0*( eightthird*rLR*u2*oom3u*oom3u + 2.0*eightthird*u*oom3u );

        /* 1SF -- el = 3+, i.e.,  even parity terms */
        REAL8 Denom3    = 1./(1. + C5*u2);
        REAL8 A3tilde   = 7.5*u*( 1 + C1*u + C2*u2 + C3*u3 )*( 1 + C4*u + C6*u2 )*Denom3;
        REAL8 dA3tilde  = 7.5*( 1 + 3*C2*u2 + 3*C6*u2 + 4*C3*u3 + 5*C2*C6*u4 + 6*C3*C6*u5 + C1*u*(2 + 3*C4*u + 4*C6*u2) + C4*u*(2 + 4*C2*u2 + 5*C3*u3) )*Denom3 + ( -15.*C5*u2*(1. + C4*u + C6*u2)*(1. + C1*u + C2*u2 + C3*u3) )*Denom3*Denom3;
        REAL8 d2A3tilde = 15.*( C1*(1 + 3*C4*u - 3*C5*pow(u,2) + 6*C6*pow(u,2) - C4*C5*pow(u,3) + 3*C5*C6*pow(u,4) + pow(C5,2)*C6*pow(u,6)) + C4*(1 - 3*C5*pow(u,2) + 10*C3*pow(u,3) +  9*C3*C5*pow(u,5) + 3*C3*pow(C5,2)*pow(u,7) +  C2*pow(u,2)*(6 + 3*C5*pow(u,2) + pow(C5,2)*pow(u,4))) + u*(3*(C6 + 2*C3*u + 5*C3*C6*pow(u,3)) + C5*(-3 - C6*pow(u,2) + 3*C3*pow(u,3) +17*C3*C6*pow(u,5)) + pow(C5,2)*(pow(u,2) + C3*pow(u,5) + 6*C3*C6*pow(u,7)) + C2*(3 + 10*C6*pow(u,2) + 3*pow(C5,2)*C6*pow(u,6) + C5*pow(u,2)*(-1 + 9*C6*pow(u,2)))) )*Denom3*Denom3*Denom3;
        REAL8 A3hat1GSFfit = A3tilde*pow(oom3u, 3.5);
        REAL8 dA3hat1GSFfit = 3.5*rLR*A3tilde*pow(oom3u, 4.5) + dA3tilde*pow(oom3u, 3.5);
        REAL8 d2A3hat1GSFfit = 15.75*rLR*rLR*A3tilde*pow(oom3u, 5.5) + 7.0*rLR*dA3tilde*pow(oom3u, 4.5) + d2A3tilde*pow(oom3u, 3.5);

        /* 2SF -- el = 3+, i.e.,  even parity terms */
        REAL8 A3hat2GSF     =  36.666666666666666667*u2*pow(oom3u,p);
        REAL8 dA3hat2GSF    =  36.666666666666666667*u*( 2. + (p - 2.)*rLR*u ) * pow(oom3u, p+1);
        REAL8 d2A3hat2GSF   =  36.666666666666666667*( 2. + 4.*(p - 1.)*rLR*u + (2. - 3.*p + 1.*p*p)*rLR*rLR*u2 ) * pow(oom3u, p+2);

        /* Hatted el = 3+ potential as a GSF series */
        REAL8 A3hatA   = A3hat_Sch + XA*A3hat1GSFfit + XA*XA*A3hat2GSF;
        REAL8 dA3hatA  = dA3hat_Sch + XA*dA3hat1GSFfit + XA*XA*dA3hat2GSF;
        REAL8 A3hatB   = A3hat_Sch + XB*A3hat1GSFfit + XB*XB*A3hat2GSF;
        REAL8 dA3hatB  = dA3hat_Sch + XB*dA3hat1GSFfit + XB*XB*dA3hat2GSF;

        /* Total el = 3+ tidal potential */
        REAL8 AT3      = -1.*kapA3*u8*( A3hatA ) - 1.*kapB3*u8*( A3hatB );
        REAL8 dAT3     = -1.*kapA3*u7*( 8.*A3hatA + 1.*u*dA3hatA ) - 1.*kapB3*u7*( 8.*A3hatB + 1.*u*dA3hatB );

        A     = AT2   + AT3   + AT4;
        dA_u  = dAT2  + dAT3  + dAT4;

        if (d2AT != NULL) {
            REAL8 d2f23  = 2*d2*(-1 + 3*d2*u2 + n1*u*(-3+d2*u2))*(Den*Den*Den);
            REAL8 d2A1SF = d2Acub*f23 + 2*dAcub*df23 + Acub*d2f23;
            REAL8 d2A2SF = 674./28.;
            REAL8 d2f0   = 6*(oom3u*oom3u*oom3u);
            REAL8 d2f1   = 0.25*(63*(rLR*rLR)*A1SF + 4*(-1+rLR*u)*(-7*rLR*dA1SF + (-1+rLR*u)*d2A1SF))*pow(oom3u,11./2.);
            REAL8 d2f2   = (  rLR*p*(1+p)*rLR*A2SF +(-1+rLR*u)*( -2.*p*rLR*dA2SF +(-1.+rLR*u)*d2A2SF )  )*pow(oom3u,p+2);

            REAL8 d2AT2  = - kapA2*30*u4*( f0 + XA*f1 + XA*XA*f2 ) - kapB2*30*u4*( f0 + XB*f1 + XB*XB*f2 ) - 2*kapA2*6*u5*( df0 + XA*df1 + XA*XA*df2 ) - 2*kapB2*6*u5*( df0 + XB*df1 + XB*XB*df2 ) - kapA2*u6*( d2f0 + XA*d2f1 + XA*XA*d2f2 ) - kapB2*u6*( d2f0 + XB*d2f1 + XB*XB*d2f2 );
            REAL8 d2AT4  = - kapT4*90*u8;

            REAL8 d2A3hatA = d2A3hat_Sch + XA*d2A3hat1GSFfit + XA*XA*d2A3hat2GSF;
            REAL8 d2A3hatB = d2A3hat_Sch + XB*d2A3hat1GSFfit + XB*XB*d2A3hat2GSF;
            REAL8 d2AT3 = -1.*kapA3 * ( 56.*u6*A3hatA + 16.*u7*dA3hatA + 1.*u8*d2A3hatA ) - 1.*kapB3 * ( 56.*u6*A3hatB + 16.*u7*dA3hatB + 1.*u8*d2A3hatB );

            d2A_u += d2AT2  + d2AT3  + d2AT4;
        }
    }


    if (dyn->use_tidal_gravitomagnetic==LAL_SIM_INSPIRAL_GMTIDES_PN) {

        /* PN series for the (2-) tidal potential */
        A    +=-kapT2j*u7*(1. +  bar_alph2j_1*u);
        dA_u += -kapT2j*u7*bar_alph2j_1 - 7.*kapT2j*u6*(1. +  bar_alph2j_1*u);

        if (d2AT != NULL) {
            d2A_u += - 14.*kapT2j*u5*(3. + 4.*bar_alph2j_1*u);
        }

    } else if (dyn->use_tidal_gravitomagnetic==LAL_SIM_INSPIRAL_GMTIDES_GSF) {

        /* GSF series for the (2-) tidal potential */
        const REAL8 a1j =  0.728591192;
        const REAL8 a2j =  3.100367557;
        const REAL8 n1j = -15.04421708;
        const REAL8 d2j =  12.55229698;
        // Schwarzschild gravito-magnetic term
        REAL8 Ahat_Schj     =  (1.-2.*u)*oom3u;
        REAL8 dAhat_Schj    =  (rLR-2.)*oom3u*oom3u;
        REAL8 d2Ahat_Schj   =  2.*rLR*(rLR-2.)*pow(oom3u, 3.);

        /* 1SF -- el = 2 gravitomagnetic terms */
        REAL8 Denomj = 1./(1. + d2j*u2);
        REAL8 Ahat1GSFfitj = elsix*u*(1. - a1j*u)*(1. - a2j*u)*(1. + n1j*u)*Denomj*pow(oom3u, 3.5);
        REAL8 dAhat1GSFfitj = 0.5*elsix * Denomj * Denomj * (2 + 4*n1j*u + 5*rLR*u - 2*d2j*pow(u,2) + 3*n1j*rLR*pow(u,2) + 9*d2j*rLR*pow(u,3) + 7*d2j*n1j*rLR*pow(u,4) -    a2j*u*(4 + rLR*u*(3 + 7*d2j*pow(u,2)) + n1j*u*(6 + 2*d2j*pow(u,2) + rLR*(u + 5*d2j*pow(u,3)))) + a1j*u*(-4 - 3*rLR*u - 7*d2j*rLR*pow(u,3) - n1j*u*(6 + rLR*u + 2*d2j*pow(u,2) + 5*d2j*rLR*pow(u,3)) + a2j*u*(6 + rLR*u + 2*d2j*pow(u,2) + 5*d2j*rLR*pow(u,3) + n1j*u*(8 - rLR*u + 4*d2j*pow(u,2) + 3*d2j*rLR*pow(u,3)))) ) * pow(oom3u, 4.5);
        REAL8 d2Ahat1GSFfitj = 0.25*elsix * Denomj * Denomj * Denomj * ( 8*(1 + n1j*u)*pow(-1 + rLR*u,2)*pow(1 + d2j*pow(u,2),2)*(-a2j + a1j*(-1 + 3*a2j*u)) +    4*(1 - rLR*u)*(1 + d2j*pow(u,2))*(1 - 2*a2j*u + a1j*u*(-2 + 3*a2j*u))*(-4*d2j*u + rLR*(7 + 11*d2j*pow(u,2)) + n1j*(2 - 2*d2j*pow(u,2) + rLR*u*(5 + 9*d2j*pow(u,2)))) + u*(-1 + a1j*u)*(-1 + a2j*u)*(7*rLR*(9*rLR + n1j*(4 + 5*rLR*u)) + 2*d2j*(-4 - 20*rLR*u + 87*pow(rLR,2)*pow(u,2) +3*n1j*u*(-4 + 8*rLR*u + 17*pow(rLR,2)*pow(u,2))) + pow(d2j,2)*pow(u,2)*(24 - 104*rLR*u + 143*pow(rLR,2)*pow(u,2) + n1j*u*(8 - 44*rLR*u + 99*pow(rLR,2)*pow(u,2)))) ) * pow(oom3u, 5.5);

        /* 2SF -- el = 2 gravitomagnetic terms */
        REAL8 Ahat2GSFj    =  u*pow(oom3u,p);
        REAL8 dAhat2GSFj    =  ( 1.+ (p-1.)*rLR*u ) * pow(oom3u, p+1);
        REAL8 d2Ahat2GSFj   =  p*rLR * ( 2.+ (p-1.)*rLR*u ) * pow(oom3u, p+2);

        /* Total el = 2 gravitomagnetic potential as a GSF series */
        REAL8 AhatjA   = Ahat_Schj + XA*Ahat1GSFfitj + XA*XA*Ahat2GSFj;
        REAL8 dAhatjA  = dAhat_Schj + XA*dAhat1GSFfitj + XA*XA*dAhat2GSFj;
        REAL8 AhatjB   = Ahat_Schj + XB*Ahat1GSFfitj + XB*XB*Ahat2GSFj;
        REAL8 dAhatjB  = dAhat_Schj + XB*dAhat1GSFfitj + XB*XB*dAhat2GSFj;

        /* el = 2 gravitomagnetic total contribution */
        REAL8 ATj_2      = -1.*kapA2j*u7*( AhatjA ) - 1.*kapB2j*u7*( AhatjB );
        REAL8 dATj_2     = -1.*kapA2j * ( 7.*u6*AhatjA + u7*dAhatjA ) - 1.*kapB2j * ( 7.*u6*AhatjB + u7*dAhatjB );

        A    += ATj_2;
        dA_u += dATj_2;

        if (d2AT != NULL) {
            REAL8 d2AhatjA = d2Ahat_Schj + XA*d2Ahat1GSFfitj + XA*XA*d2Ahat2GSFj;
            REAL8 d2AhatjB = d2Ahat_Schj + XB*d2Ahat1GSFfitj + XB*XB*d2Ahat2GSFj;
            REAL8 d2ATj_2    = -1.*kapA2j * ( 42.*u5*AhatjA + 14.*u6*dAhatjA + u7*d2AhatjA ) - 1.*kapB2j * ( 42.*u5*AhatjB + 14.*u6*dAhatjB + u7*d2AhatjB );

            d2A_u += d2ATj_2;
        }

    }

    *AT   = A;
    *dAT  = dA_u;
    if (d2AT != NULL) *d2AT = d2A_u;

    return;
}

/* EOB Metric potentials A(r), B(r), and their derivatives, no spin version */
void eob_metric(REAL8 r,
                LALTEOBResumSDynamics *dyn,
                REAL8 *A,
                REAL8 *B,
                REAL8 *dA,
                REAL8 *d2A,
                REAL8 *dB)
{
    const REAL8 nu    = dyn->nu;
    const REAL8 u     = 1./r;
    const REAL8 u2    = u*u;
    const REAL8 u3    = u2*u;
    const REAL8 u4    = u2*u2;
    const REAL8 UNUSED u6    = u2*u4;

    REAL8 Atmp=0., dAtmp_u=0., d2Atmp_u=0.;
    REAL8 Btmp=0., dBtmp_r=0.;

    /* A potential and derivative with respect to u */
    eob_metric_A5PNlog(r, nu, &Atmp, &dAtmp_u, &d2Atmp_u);

    /* Add here tides if needed */
    if (dyn->use_tidal) {
        REAL8 AT, dAT_u, d2AT_u;
        REAL8 UNUSED BT, UNUSED dBT;
        eob_metric_Atidal(r, dyn, &AT, &dAT_u, &d2AT_u);
        Atmp     += AT;
        dAtmp_u  += dAT_u;
        d2Atmp_u += d2AT_u;
#if (USEBTIDALPOTENTIAL)
        /* Vines, Flanagan 1PN term in B */
        REAL8 kT2 = dyn->kapT2;
        BT  = kT2*(8. - 15.*nu)*u6;
        dBT = -kT2*6.*(8. - 15.*nu)*u4*u3;
        Btmp    += BT;
        dBtmp_r += dBT;
#endif
    }

    /* A potential and derivative with respect to r */
    *A   = Atmp;
    *dA  = -dAtmp_u*u2;
    *d2A = 2.*dAtmp_u*u3 + d2Atmp_u*u4;

    /* D potential and derivative with respect to r */
    const REAL8 Dp  = 1.0 + 6.*nu*u2 - 2.*(3.0*nu-26.0)*nu*u3; // Pade' resummation of D
    const REAL8 D   = 1./Dp;
    const REAL8 dD  = 6.*u2*(2.*nu*u-(3.*nu-26.)*nu*u2)*D*D;

    /* B potential and derivative with respect to r */
    Btmp    += D/(Atmp);
    dBtmp_r += (dD*(Atmp) - D*(*dA))/((Atmp)*(Atmp));

    *B  = Btmp;
    *dB = dBtmp_r;

}

/* EOB Metric potentials A(r), B(r), and their derivatives, spin version */
void eob_metric_s(REAL8 r, LALTEOBResumSDynamics *dyn, REAL8 *A, REAL8 *B, REAL8 *dA, REAL8 *d2A, REAL8 *dB)
{

    const REAL8 nu    = dyn->nu;
    const REAL8 a1    = dyn->a1;
    const REAL8 a2    = dyn->a2;
    const REAL8 aK2   = dyn->aK2;
    const REAL8 C_Q1  = dyn->C_Q1;
    const REAL8 C_Q2  = dyn->C_Q2;
    const REAL8 C_Oct1 = dyn->C_Oct1;
    const REAL8 C_Oct2 = dyn->C_Oct2;
    const REAL8 C_Hex1 = dyn->C_Hex1;
    const REAL8 C_Hex2 = dyn->C_Hex2;
    const int usetidal = dyn->use_tidal;

    const REAL8 u   = 1./r;
    const REAL8 u2  = u*u;
    const REAL8 u3  = u2*u;
    const REAL8 u4  = u2*u2;

    REAL8 rc, drc, d2rc;
    dyn->eob_dyn_s_get_rc(r, nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, &rc, &drc, &d2rc);

    /* A potential and derivative with respect to u */
    REAL8 Aorb, dAorb_u, d2Aorb_u;
    eob_metric_A5PNlog(rc, nu, &Aorb, &dAorb_u, &d2Aorb_u);

    /* Add here tides if needed */
    if (usetidal) {
        REAL8 AT, dAT_u, d2AT_u;
        eob_metric_Atidal(rc, dyn, &AT, &dAT_u, &d2AT_u);
        Aorb     += AT;
        dAorb_u  += dAT_u;
        d2Aorb_u += d2AT_u;
    }

    /* A potential and derivative with respect to r */
    REAL8 uc  = 1./rc;
    REAL8 uc2 = uc*uc;
    REAL8 uc3 = uc2*uc;
    REAL8 uc4 = uc2*uc2;

    REAL8 dAorb  = -dAorb_u*uc2;
    REAL8 d2Aorb = 2.*dAorb_u*uc3 + d2Aorb_u*uc4;

    /* Correct A for spin */
    REAL8 AKerr_Multipole = (1.+2.*uc)/(1.+2.*u);
    REAL8 fss = 1.;

    *A   = Aorb*AKerr_Multipole*fss;
    *dA  = dAorb*drc*(1.+2.*uc)/(1.+2.*u) - 2.*Aorb*drc*uc2/(1.+2.*u) + 2.*Aorb*(1.+2.*uc)*u2/((1.+2.*u)*(1.+2.*u));
    *d2A = d2Aorb*(1.+2.*uc)/(1.+2.*u) + 4.*dAorb*( u2*(1.+2.*uc)/((1.+2.*u)*(1.+2.*u)) - uc2/(1.+2.*u)*drc) + Aorb*(-4.*u3*(1.+2.*uc)/((1.+2.*u)*(1.+2.*u)) + 8.*u4*(1.+2.*uc)/((1.+2.*u)*(1.+2.*u)*(1.+2.*u))+4.*uc3*(1.+2.*u)*drc*drc - 2.*uc2/(1.+2.*u)*d2rc);

    /* D potential and derivative with respect to r */
    REAL8 Dp = 1.0 + 6.*nu*uc2 - 2.*(3.0*nu-26.0)*nu*uc3; // Pade' resummation of D
    REAL8 D  = 1./Dp;
    REAL8 dD = 6.*uc2*(2.*nu*uc-(3.*nu-26.)*nu*uc2)*D*D;

    /* B potential and derivative with respect to r */
    *B   = r*r*uc2*D/(*A);
    *dB  = (dD*(*A) - D*(*dA))/((*A)*(*A));

    /* Add here tides if needed */

}

/* Newtonian partial fluxes */
void eob_flx_FlmNewt(REAL8 x, REAL8 nu, REAL8 *Nlm)
{

    /* Shorthands*/
    const REAL8 nu2 = nu*nu;
    const REAL8 nu3 = nu2*nu;
    const REAL8 x5  = x*x*x*x*x;
    const REAL8 x6  = x*x5;
    const REAL8 x7  = x*x6;
    const REAL8 x8  = x*x7;
    const REAL8 x9  = x*x8;
    const REAL8 x10 = x*x9;
    const REAL8 x11 = x*x10;
    const REAL8 x12 = x*x11;

    const REAL8 sp2 = 1.-4.*nu;
    const REAL8 sp4 = (1-4*nu)*SQ((1-2*nu));
    const REAL8 sp3 = (1.-3.*nu)*(1.-3.*nu);
    const REAL8 sp5 = (1.-5.*nu+5.*nu2)*(1.-5.*nu+5.*nu2);
    const REAL8 sp6 = (1-4*nu)*(3*nu2-4*nu +1)*(3*nu2-4*nu +1);
    const REAL8 sp7 = (1 - 7*nu + 14*nu2 - 7*nu3)*(1 - 7*nu + 14*nu2 - 7*nu3);
    const REAL8 sp8 = (1 - 4*nu)*(1 - 6*nu + 10*nu2 - 4*nu3)*(1 - 6*nu + 10*nu2 - 4*nu3);

    REAL8 spx[] = {
        sp2 * x6, x5,
        sp2 * x6, sp3 * x7, sp2 * x6,
        sp4 * x8, sp3 * x7, sp4 * x8, sp3 * x7,
        sp4 * x8, sp5 * x9, sp4 * x8, sp5 * x9, sp4 * x8,
        sp6 * x10, sp5 * x9, sp6 * x10, sp5 * x9, sp6 * x10, sp5 * x9,
        sp6 * x10, sp7 * x11, sp6 * x10, sp7 * x11, sp6 * x10, sp7 * x11, sp6 * x10,
        sp8 * x12, sp7 * x11, sp8 * x12, sp7 * x11, sp8 * x12, sp7 * x11, sp8 * x12, (7*nu3-14*nu2+7*nu-1)*(7*nu3-14*nu2+7*nu-1) * x11
    };

    /* Newtonian partial fluxes*/
    for (int k = 0; k < KMAX; k++) {
        Nlm[k] = CNlm[k] * spx[k];
    }

    return;
}

/* Factorials evaluated for the tail term */
static const REAL8 f14[] =
{1.,         1.,          2.,
    6.,         24.,         120.,
    720.,       5040.,       40320.,
    362880.,    3628800.,    39916800.,
    479001600., 6227020800., 87178291200.};


/* Tail term (modulus) */
void eob_flx_Tlm(const REAL8 w, REAL8 *MTlm)
{
    REAL8 hhatk, x2, y, prod;
    INT4 k, j;
    for (k = 0; k < KMAX; k++ ) {
        hhatk = TEOB_MINDEX[k] * w;
        x2    = 4.*hhatk*hhatk;
        prod  = 1.;
        for (j=1; j <= TEOB_LINDEX[k]; j++ ) {
            prod *= ( j*j + x2 );
        }
        y  = 4.*LAL_PI*hhatk;
        y /= ( 1. - exp(-y) );
        MTlm[k] = sqrt( 1./(f14[TEOB_LINDEX[k]]*f14[TEOB_LINDEX[k]]) * y * prod );
    }
    return;
}


/* Compute horizon-absorbed fluxes. no spin case.
 * Nagar & Akcay, PRD 85, 044025 (2012)
 * Bernuzzi, Nagar & Zenginoglu, PRD 86, 104038 (2012)
 */
REAL8 eob_flx_HorizonFlux(REAL8 x, REAL8 Heff, REAL8 jhat, REAL8 nu)
{
    REAL8 rhoHlm[2]; /* only 21,22 multipoles -> k=0,1 */
    REAL8 FlmHLO[2];
    REAL8 FlmH[2];

    /* Shorthands */
    REAL8 nu2 = nu*nu;
    REAL8 nu3 = nu*nu2;
    REAL8 x2  = x*x;
    REAL8 x3  = x*x2;
    REAL8 x4  = x*x3;
    REAL8 x5  = x*x4;
    REAL8 x9  = x4*x5;
    REAL8 x10 = x*x9;

    /* The Newtonian asymptotic contribution */
    const REAL8 FNewt22 = 32./5.*x5;

    /* Compute leading-order part (nu-dependent) */
    FlmHLO[1] = 32./5.*(1-4*nu+2*nu2)*x9;
    FlmHLO[0] = 32./5.*(1-4*nu+2*nu2)*x10;

    /* Compute rho_lm */
    REAL8 c1[2];
    REAL8 c2[2];
    REAL8 c3[2];
    REAL8 c4[2];

    c1[0] = 0.58121;
    c2[0] = 1.01059;
    c3[0] = 7.955729;
    c4[0] = 1.650228;

    c1[1] = (4.-21.*nu + 27.*nu2 - 8.*nu3)/(4.*(1.-4.*nu+2.*nu2));
    c2[1] =  4.78752;
    c3[1] = 26.760136;
    c4[1] = 43.861478;

    rhoHlm[1] = 1. + c1[1]*x + c2[1]*x2 + c3[1]*x3 + c4[1]*x4;
    rhoHlm[0] = 1. + c1[0]*x + c2[0]*x2 + c3[0]*x3 + c4[0]*x4;

    /* Compute horizon multipolar flux (only l=2) */
    const REAL8 Heff2 = Heff*Heff;
    const REAL8 jhat2 = jhat*jhat;

    FlmH[1] = FlmHLO[1] * Heff2 * gsl_pow_int(rhoHlm[1],4);
    FlmH[0] = FlmHLO[0] * jhat2 * gsl_pow_int(rhoHlm[0],4);

    /* Sum over multipoles and normalize to the 22 Newtonian multipole */
    REAL8 hatFH = (FlmH[0]+FlmH[1])/FNewt22;

    return hatFH;
}


/* Compute horizon-absorbed fluxes. spin case. */
REAL8 eob_flx_HorizonFlux_s(REAL8 x,
                            REAL8 UNUSED Heff,
                            REAL8 UNUSED jhat,
                            REAL8 UNUSED nu,
                            REAL8 X1,
                            REAL8 X2,
                            REAL8 chi1,
                            REAL8 chi2)
{

    REAL8 x2 = x*x;
    REAL8 x3 = x2*x;
    REAL8 x4 = x3*x;
    REAL8 x5 = x4*x;
    REAL8 v5 = sqrt(x5);

    REAL8 cv5[2];
    REAL8 cv8[2];

    /* Coefficients of the v^5 term (Alvi leading order) */
    cv5[0] = -1./4.*chi1*(1.+3.*chi1*chi1)*X1*X1*X1;
    cv5[1] = -1./4.*chi2*(1.+3.*chi2*chi2)*X2*X2*X2;

    /* Coefficients of the v^8=x^4 term */
    cv8[0] = 0.5*(1.+sqrt(1.-chi1*chi1))*(1.+3.*chi1*chi1)*X1*X1*X1*X1;
    cv8[1] = 0.5*(1.+sqrt(1.-chi2*chi2))*(1.+3.*chi2*chi2)*X2*X2*X2*X2;

    REAL8 FH22_S = (cv5[0]+cv5[1])*v5;
    REAL8 FH22   = (cv8[0]+cv8[1])*x4;
    REAL8 FH21   =  0.0;

    /* Newton-normalized horizon flux: use only l=2 fluxes */
    REAL8 hatFH  = FH22_S + FH22 + FH21;

    return hatFH;
}


/* Flux calculation for spinning systems */
// TODO: NQC are not applied in spin case!
REAL8 eob_flx_Flux_s(REAL8 x,
                     REAL8 Omega,
                     REAL8 r_omega,
                     REAL8 E,
                     REAL8 Heff,
                     REAL8 jhat,
                     REAL8 r,
                     REAL8 pr_star,
                     REAL8 ddotr,
                     LALTEOBResumSDynamics *dyn)
{
    const REAL8 nu = dyn->nu;
    const REAL8 chi1 = dyn->chi1;
    const REAL8 chi2 = dyn->chi2;
    const REAL8 X1 = dyn->X1;
    const REAL8 X2 = dyn->X2;
    const REAL8 a1 = dyn->a1;
    const REAL8 a2 = dyn->a2;
    const REAL8 C_Q1 = dyn->C_Q1;
    const REAL8 C_Q2 = dyn->C_Q2;
    const REAL8 X12  = X1-X2; /* sqrt(1-4nu) */
    const REAL8 UNUSED X12sq = SQ(X12); /* (1-4nu) */

    const INT4 usetidal = dyn->use_tidal;
    const INT4 usespins = dyn->use_spins;

    REAL8 prefact[] = {
        jhat, Heff,
        Heff, jhat, Heff,
        jhat, Heff, jhat, Heff,
        Heff, jhat, Heff, jhat, Heff,
        jhat, Heff, jhat, Heff, jhat, Heff,
        Heff, jhat, Heff, jhat, Heff, jhat, Heff,
        jhat, Heff, jhat, Heff, jhat, Heff, jhat, Heff};

    REAL8 FNewt22, sum_k=0.;
    REAL8 rholm[KMAX], flm[KMAX], FNewtlm[KMAX], MTlm[KMAX], hlmTidal[KMAX], hlmNQC[KMAX];
    REAL8 Modhhatlm[KMAX];

    /* Newtonian flux */
    eob_flx_FlmNewt(x, nu, FNewtlm);

    /* Correct amplitudes for specific multipoles and cases */
    if (usespins) {
        /* Correct (2,1), (3,1) and (3,3) ( sp2 = 1 ) */
        REAL8 x6 = gsl_pow_int(x, 6);
        FNewtlm[0] = CNlm[0] * x6; /* (2,1) */
        FNewtlm[2] = CNlm[2] * x6; /* (3,1) */
        FNewtlm[4] = CNlm[4] * x6; /* (3,3) */
        /* Correct (4,1), (4,3)  ( sp4 = (1-2nu)^2 ) */
        REAL8 sp4x8 = SQ((1-2*nu)) * gsl_pow_int(x, 8);
        FNewtlm[5] = CNlm[5] * sp4x8; /* (4,1) */
        FNewtlm[7] = CNlm[7] * sp4x8; /* (4,3) */
    } else {
        if (usetidal) {
            /* Correct (2,1), (3,1) and (3,3) ( sp2 = 1 ) */
            REAL8 x6 = gsl_pow_int(x, 6);
            FNewtlm[0] = CNlm[0] * x6; /* (2,1) */
            FNewtlm[2] = CNlm[2] * x6; /* (3,1) */
            FNewtlm[4] = CNlm[4] * x6; /* (3,3) */
        }
    }

    /* Tail term */
    eob_flx_Tlm(E*Omega, MTlm);

    /* Amplitudes */
    if (usespins) {
        //dyn->eob_wav_flm_s_old(x,nu, X1,X2,chi1,chi2,a1,a2,C_Q1,C_Q2, usetidal, rholm, flm);
        dyn->eob_wav_flm_s(x, nu, X1, X2, chi1, chi2, a1, a2, C_Q1, C_Q2, dyn->clm, usetidal, rholm, flm);
    } else {
        //eob_wav_flm_old(x,nu, rholm, flm);
        eob_wav_flm(x, nu, dyn->clm, rholm, flm);
    }

    FNewt22 = FNewtlm[1];

    /* NQC correction to the modulus of the (l,m) waveform */
    for (int k = 0; k < KMAX; k++) hlmNQC[k] = 1.; /* no NQC */

    if (dyn->nqc_coeffs_flx != NQC_OFF)
    {
        LALTEOBResumSWaveformModeSingleTime hNQC;
        /* eob_wav_hlmNQC_nospin201602(nu,r,pr_star,Omega,ddotr, &hNQC); */
        eob_wav_hlmNQC(nu, r, pr_star, Omega, ddotr, dyn->NQC->flx, &hNQC);
        /*
         const INT4 UNUSED maxk = MIN(KMAX, NQC->hlm->maxk+1);
         for (int k = 0; k < maxk; k++) {
         if (NQC->hlm->activemode[k]) {
         hlmNQC[k] = hNQC.ampli[k];
         }
         }
         */
        /* Use only the 22:  */
        hlmNQC[1] = hNQC.ampli[1];
    }

    /* Compute modulus of hhat_lm (with NQC) */
    for (int k = 0; k < KMAX; k++) {
        Modhhatlm[k] = prefact[k] * MTlm[k] * flm[k] * hlmNQC[k];
    }

    if (usetidal) {
        /* Tidal amplitudes */
        eob_wav_hlmTidal(x, dyn, hlmTidal);
        if (!(usespins)) {
            /* Correct normalization of (2,1) (3,1), (3,3) point-mass amplitudes */
            Modhhatlm[0] *= X12;
            Modhhatlm[2] *= X12;
            Modhhatlm[4] *= X12;
        }
        /* Add tidal amplitudes */
        for (int k = 0; k < KMAX; k++) {
            Modhhatlm[k] += MTlm[k] * hlmTidal[k];
        }
    }

    /* Total multipolar flux */
    for (int k = KMAX; k--;) sum_k += SQ(Modhhatlm[k]) * FNewtlm[k];

    /* Normalize to the 22 Newtonian multipole */
    REAL8 hatf = sum_k/(FNewt22);

    /* Horizon flux */
    if (!(usetidal)) {
        REAL8 hatFH;
        if (usespins) {
            hatFH = eob_flx_HorizonFlux_s(x, Heff, jhat, nu, X1, X2, chi1, chi2);
        } else {
            hatFH = eob_flx_HorizonFlux(x,Heff,jhat,nu);
        }
        hatf += hatFH;
    }

    /* return Fphi */
    return (-32./5. * nu * gsl_pow_int(r_omega,4) * gsl_pow_int(Omega,5) * hatf);
}


/* Flux calculation for Newton-Normalized energy flux
 Use the DIN resummation procedure.
 Add non-QC and non-K corrections to (2,2) partial flux. */
REAL8 eob_flx_Flux(REAL8 x,
                   REAL8 Omega,
                   REAL8 r_omega,
                   REAL8 E,
                   REAL8 Heff,
                   REAL8 jhat,
                   REAL8 r,
                   REAL8 pr_star,
                   REAL8 ddotr,
                   LALTEOBResumSDynamics *dyn)
{
    return eob_flx_Flux_s(x, Omega, r_omega, E, Heff, jhat, r, pr_star, ddotr,dyn);
}


/* Alternative implementation of the phase of the tail factor */
void eob_wav_speedyTail(REAL8 Omega, REAL8 Hreal, REAL8 bphys, LALTEOBResumSWaveformModeSingleTime *tlm)
{
    REAL8 x;
    REAL8 x2;
    REAL8 x3;
    REAL8 x4;
    REAL8 x5;
    REAL8 tlm_ang;
    REAL8 num_ang;

    /* Fit coefficients*/
    const REAL8 b1[] = {
        0.1113090643348557, 0.1112593821157397, 0.0424759238428813, 0.0424489015884926, 0.0424717446800903, 0.0215953972500844, 0.0215873812155663, 0.0215776183122621, 0.0216017621863542, 0.0128123696874894, 0.0128097056242375, 0.0128038943888768, 0.0128025242617949, 0.0128202485907368, 0.0083762045692408, 0.0083751913886140, 0.0083724067460769, 0.0083694435961860, 0.0083710364141552, 0.0083834483913443, 0.0058540393221396, 0.0058536069384738, 0.0058522594457692, 0.0058502436535615, 0.0058491157293566, 0.0058514875071582, 0.0058602498033381, 0.0042956812356573, 0.0042954784390887, 0.0042947951664056, 0.0042935886137697, 0.0042923691461384, 0.0042922256848799, 0.0042945927126022, 0.0043009106861259};

    const REAL8 b2[] = {
        0.0004643273300862, 0.0009375605440004, 0.0000597134489198, 0.0002551406918111, 0.0001741036904709, 0.0000124649041611, 0.0000685496215625, 0.0001131160409390, 0.0000419907542591, 0.0000035218982282, 0.0000219211271097, 0.0000473186962874, 0.0000524142634057, 0.0000106823372552, 0.0000012237574387, 0.0000081742188269, 0.0000201940563214, 0.0000295722761753, 0.0000260539631956, 0.0000018994753518, 0.0000004932942990, 0.0000034477210351, 0.0000092294406360, 0.0000155143073237, 0.0000183386499818, 0.0000137922469695, -0.0000007075155453, 0.0000002223410995, 0.0000016045317657, 0.0000045260028113, 0.0000082655700107, 0.0000112393599417, 0.0000115758243113, 0.0000076838709956, -0.0000014020591745};

    const REAL8 b3[] = {
        -0.0221835462237291, -0.0235386333304348, -0.0042911639711832, -0.0047431560217121, -0.0046577314472149, -0.0013089557502947, -0.0014343968205390, -0.0014978542575474, -0.0014329302934532, -0.0005167994164556, -0.0005573939123058, -0.0005921030407223, -0.0005978284714483, -0.0005673965369076, -0.0002409269302708, -0.0002561516055118, -0.0002723768586352, -0.0002815958312453, -0.0002792078156272, -0.0002646630240693, -0.0001261183503407, -0.0001325622938779, -0.0001403198638518, -0.0001464084186977, -0.0001485971591029, -0.0001459023931717, -0.0001384829633836, -0.0000719062974278, -0.0000749128468013, -0.0000788187384314, -0.0000824202283094, -0.0000846673495936, -0.0000849054394951, -0.0000829269749240, -0.0000788883333858};

    const REAL8 b4[] = {
        0.0058366730167965, 0.0070452306758401, 0.0006914295465364, 0.0010322294603561, 0.0010057563135650, 0.0001394203795507, 0.0002309706405978, 0.0002596611624417, 0.0002409588083156, 0.0000386949167221, 0.0000679154947896, 0.0000830199015202, 0.0000850120755064, 0.0000780125513602, 0.0000133034384660, 0.0000241813441339, 0.0000311573885555, 0.0000340233089866, 0.0000335167900637, 0.0000307571022927, 0.0000053305073331, 0.0000099143129290, 0.0000132296989826, 0.0000150959309402, 0.0000156304390748, 0.0000151274875147, 0.0000139320508803, 0.0000023959090314, 0.0000045285807761, 0.0000061918979830, 0.0000072894226381, 0.0000078251853305, 0.0000078772667984, 0.0000075606242809, 0.0000069956215270
    };

    REAL8 Tlm_real[KMAX];
    eob_flx_Tlm(Omega*Hreal, Tlm_real);

    /* Pre-computed psi */
    const REAL8 psi[] = {0.9227843350984671394, 0.9227843350984671394,
        1.256117668431800473, 1.256117668431800473, 1.256117668431800473,
        1.506117668431800473, 1.506117668431800473, 1.506117668431800473, 1.506117668431800473,
        1.706117668431800473, 1.706117668431800473, 1.706117668431800473, 1.706117668431800473, 1.706117668431800473,
        1.872784335098467139, 1.872784335098467139, 1.872784335098467139, 1.872784335098467139, 1.872784335098467139, 1.872784335098467139,
        2.015641477955609997, 2.015641477955609997, 2.015641477955609997, 2.015641477955609997, 2.015641477955609997, 2.015641477955609997, 2.015641477955609997,
        2.140641477955609997, 2.140641477955609997, 2.140641477955609997, 2.140641477955609997, 2.140641477955609997, 2.140641477955609997, 2.140641477955609997, 2.140641477955609997};

    REAL8 k;
    int i;
    // TODO: [optimization] vectorize?
    for (i=0; i<KMAX; i++) {
        k  = TEOB_MINDEX[i] * Omega;
        x  = k * Hreal; /* hathatk */
        x2 = x * x;
        x3 = x2 * x;
        x4 = x3 * x;
        x5 = x4 * x;
        num_ang   = 1. + b1[i]*x2 + b2[i]*x3 + b3[i]*x4 + b4[i]*x5;
        tlm_ang   = (- 2. * psi[i] * x * num_ang) + 2.*x* log(2. * k * bphys);
        tlm->ampli[i] = Tlm_real[i];
        tlm->phase[i] = tlm_ang;
    }

}

/* Tail contribution to the resummed wave.
 Ref. Damour, Iyer & Nagar, PRD 79, 064004 (2009) */
void eob_wav_hhatlmTail(REAL8 Omega, REAL8 Hreal, REAL8 bphys, LALTEOBResumSWaveformModeSingleTime *tlm)
{
    REAL8 k;
    REAL8 hhatk;

    gsl_sf_result num_rad;
    gsl_sf_result num_phase;
    gsl_sf_result denom_rad;
    gsl_sf_result denom_phase;

    REAL8 ratio_rad;
    REAL8 ratio_ang;
    REAL8 tlm_rad;
    REAL8 tlm_phase;

    int i;
    for (i = 0; i < KMAX; i++) {
        k     = TEOB_MINDEX[i] * Omega;
        hhatk = k * Hreal;

        gsl_sf_lngamma_complex_e((double) TEOB_LINDEX[i] + 1., (double) (-2.*hhatk), &num_rad, &num_phase);
        gsl_sf_lngamma_complex_e((double) TEOB_LINDEX[i] + 1., 0., &denom_rad, &denom_phase);

        ratio_rad     = (REAL8) (num_rad.val - denom_rad.val);
        ratio_ang     = (REAL8) num_phase.val - 0.;

        tlm_rad       = ratio_rad + LAL_PI * hhatk;
        tlm_phase     = ratio_ang + 2.*hhatk*log(2.*k*bphys);

        tlm->ampli[i] = exp(tlm_rad);
        tlm->phase[i] = tlm_phase;
    }

}


/* Resummed amplitudes for the spin case.
 This function computes the residual amplitude corrections flm's as
 introduced in Damour, Iyer & Nagar, PRD 79, 064004 (2008).
 The orbital part is taken at the usual 3^{+2} PN order, i.e. 3PN terms
 are integrated by the 4PN and 5PN test-particle terms, with the higher
 modes obtained by Fujita & Iyer.
 It only includes spin-spin interaction at LO for the (2,2) mode.
 Note that the variables called here (a1,a2)
 are what we usually cal \f$tilde{a}_1\f$ and \f$tilde{a}_2\f$ and are defined as
 a1 = X1*chi1, a2=X2*chi2 and are passed here as parameters. Special
 combinations of these quantities are used here to write the spin-dependent
 part of the waveform in particularly compact form, so that the (spinning)
 test-particle limit is recovered just by visual inspection of the equations */
void eob_wav_flm_s_SSLO(REAL8 x,
                        REAL8 nu,
                        REAL8 X1,
                        REAL8 X2,
                        REAL8 UNUSED chi1,
                        REAL8 UNUSED chi2,
                        REAL8 a1,
                        REAL8 a2,
                        REAL8 C_Q1,
                        REAL8 C_Q2,
                        REAL8 clm[KMAX][6],
                        int usetidal,
                        REAL8 *rholm,
                        REAL8 *flm)
{

    /* Orbital part */
    eob_wav_flm(x, nu, clm, rholm, flm);

    /* Spin corrections */
    REAL8 rho22S;
    REAL8 rho32S;
    REAL8 rho44S;
    REAL8 rho42S;
    REAL8 f21S;
    REAL8 f33S;
    REAL8 f31S;
    REAL8 f43S;
    REAL8 f41S;

    const REAL8 a0      = a1+a2;
    const REAL8 a12     = a1-a2;
    const REAL8 X12     = X1-X2;
    const REAL8 a0X12   = a0*X12;
    const REAL8 a12X12  = a12*X12;

    const REAL8 v  = sqrt(x);
    const REAL8 v2 = x;
    const REAL8 v3 = v*v2;
    const REAL8 v4 = v3*v;
    const REAL8 v5 = v4*v;

    /* l=m=2 multipole */
    /* spin-orbit */
    const REAL8 cSO_lo    = (-0.5*a0 - a12X12/6.);
    const REAL8 cSO_nlo   = (-52./63.-19./504.*nu)*a0 - (50./63.+209./504.*nu)*a12X12;

    /* SPIN-SPIN contribution */
    REAL8 cSS_lo;
    if (usetidal) {
#if (RC_EXCLUDESPINSPINTIDES)
        /* Switch off spin-spin-tidal couplings */
        /* See also: eob_dyn_s_get_rc() */
        cSS_lo = 0.;
        /* Above code switch off everything,
         Alt. one can set C_Q1=C_Q2=0, but keep the term: */
        /*
         cSS_lo = a1*a2;
         */
#else
        cSS_lo = 0.5*(C_Q1*a1*a1 + 2.*a1*a2 + C_Q2*a2*a2);
#endif
    } else {
        cSS_lo = 0.5*a0*a0;
    }

    /* rho_22^S: Eq. (80) of Damour & Nagar, PRD 90, 044018 (2014) */
    rho22S = cSO_lo*v3 + cSS_lo*v4 + cSO_nlo*v5 ;

    /* l>=3, m=even: multipoles rewritten in compact and self-explanatory form */
    rho32S = (a0-a12X12)/(3.*(1.-3.*nu))*v;
    rho44S = (-19./30.*a0 -  (1.-21.*nu)/(30.-90.*nu)*a12X12)*v3;
    rho42S = ( -1./30.*a0 - (19.-39.*nu)/(30.-90.*nu)*a12X12)*v3;

    /* l>=2, m=odd: multipoles rewritten in compact and self-explanatory form */
    f21S = -1.5*a12*v + ((110./21. + 79./84.*nu)*a12 - 13./84.*a0X12)*v3;
    f33S = ((-0.25 + 2.5*nu)*a12 - 1.75*a0X12)*v3;
    f31S = ((-2.25 + 6.5*nu)*a12 + 0.25*a0X12)*v3;
    f43S = (( 5. -10.*nu)*a12 - 5.*a0X12)/(-4.+8.*nu)*v;
    f41S = f43S;

    /* Amplitudes (correct with spin terms) */
    flm[0] = SQ(rholm[0]);
    flm[0] = (X12*flm[0] + f21S);

    flm[1] = SQ(rholm[1]+ rho22S);

    flm[2] = gsl_pow_int(rholm[2], 3);
    flm[2] = (X12*flm[2] + f31S);

    flm[3] = gsl_pow_int(rholm[3]+ rho32S, 3);

    flm[4] = gsl_pow_int(rholm[4], 3);
    flm[4] = (X12*flm[4] + f33S);

    flm[5] = gsl_pow_int(rholm[5], 4);
    flm[5] = (X12*flm[5] + f41S);

    flm[6] = gsl_pow_int(rholm[6] + rho42S, 4);

    flm[7] = gsl_pow_int(rholm[7], 4);
    flm[7] = (X12*flm[7] + f43S);

    flm[8] = gsl_pow_int(rholm[8] + rho44S, 4);

}

/*
 Resummed amplitudes for the spin case.
 This function computes the residual amplitude corrections flm's as
 introduced in Damour, Iyer & Nagar, PRD 79, 064004 (2008).
 The orbital part is taken at the usual 3^{+2} PN order, i.e. 3PN terms
 are integrated by the 4PN and 5PN test-particle terms, with the higher
 modes obtained by Fujita & Iyer.
 The function includes spin-spin interaction at NLO for the (2,2) mode
 and at LO for the (2,1),(3,1) and (3,3) modes.
 Note that the variables called here (a1,a2)
 are what we usually cal \f$tilde{a}_1\f$ and \f$tilde{a}_2\f$ and are defined as
 a1 = X1*chi1, a2=X2*chi2 and are passed here as parameters. */
void eob_wav_flm_s_SSNLO(REAL8 x,
                         REAL8 nu,
                         REAL8 X1,
                         REAL8 X2,
                         REAL8 UNUSED chi1,
                         REAL8 UNUSED chi2,
                         REAL8 a1,
                         REAL8 a2,
                         REAL8 C_Q1,
                         REAL8 C_Q2,
                         REAL8 clm[KMAX][6],
                         int usetidal,
                         REAL8 *rholm,
                         REAL8 *flm)
{

    /* Orbital part */
    eob_wav_flm(x, nu, clm, rholm, flm);

    /* Spin corrections */
    REAL8 rho22S;
    REAL8 rho32S;
    REAL8 rho44S;
    REAL8 rho42S;
    REAL8 f21S;
    REAL8 f33S;
    REAL8 f31S;
    REAL8 f43S;
    REAL8 f41S;

    const REAL8 a0      = a1+a2;
    const REAL8 a12     = a1-a2;
    const REAL8 X12     = X1-X2;
    const REAL8 a0X12   = a0*X12;
    const REAL8 a12X12  = a12*X12;

    const REAL8 v  = sqrt(x);
    const REAL8 v2 = x;
    const REAL8 v3 = v2*v;
    const REAL8 v4 = v3*v;
    const REAL8 v5 = v4*v;
    const REAL8 v6 = v5*v;
    const REAL8 UNUSED v7 = v6*v;

    /* l=m=2 multipole */
    /* spin-orbit */
    const REAL8 cSO_lo    = (-0.5*a0 - a12X12/6.);
    const REAL8 cSO_nlo   = (-52./63.-19./504.*nu)*a0 - (50./63.+209./504.*nu)*a12X12;
    const REAL8 UNUSED cSO_nnlo  = (32873./21168 + 477563./42336.*nu + 147421./84672.*nu*nu)*a0 - (23687./63504 - 171791./127008.*nu + 50803./254016.*nu*nu)*a12X12;

    /* SPIN-SPIN contribution */
    REAL8 cSS_lo = 0.;
    REAL8 cSS_nlo = 0.;
    if (usetidal) {
#if (RC_EXCLUDESPINSPINTIDES)
        /* Switch off spin-spin-tidal couplings */
        /* See also: eob_dyn_s_get_rc() */
        cSS_lo  = 0.;
        cSS_nlo = 0.;
        /* Above code switch off everything,
         Alt. one can set C_Q1=C_Q2=0, but keep the term: */
        /*
         cSS_lo = a1*a2;
         */
#else
        cSS_lo  = 0.5*(C_Q1*a1*a1 + 2.*a1*a2 + C_Q2*a2*a2);
        cSS_nlo = (-85./63. + 383./252.*nu)*a1*a2 + (-2./3. - 5./18.*nu)*(a1*a1 + a2*a2) + (1./7. + 27./56.*nu)*(C_Q1*a1*a1 + C_Q2*a2*a2) + 2./9.*X12*(a1*a1 - a2*a2) + 55./84.*X12*(C_Q1*a1*a1 - C_Q2*a2*a2);
#endif
    } else {
        cSS_lo  = 0.5*a0*a0;
        cSS_nlo = 1./504.*(2.*(19. - 70.*nu)*a12*a12 + (-302. + 243.*nu)*a0*a0 + 442.*X12*a0*a12);
    }

    /* rho_22^S: Eq. (80) of Damour & Nagar, PRD 90, 044018 (2014) */
    rho22S = cSO_lo*v3 + cSS_lo*v4 + cSO_nlo*v5;

    // Adding NLO SS term w.r.t. eob_wav_flm_s_SSLO
    rho22S += cSS_nlo*v6;

    /* l>=3, m=even: multipoles rewritten in compact and self-explanatory form */
    rho32S = (a0-a12X12)/(3.*(1.-3.*nu))*v;
    rho44S = (-19./30.*a0 -  (1.-21.*nu)/(30.-90.*nu)*a12X12)*v3;
    rho42S = ( -1./30.*a0 - (19.-39.*nu)/(30.-90.*nu)*a12X12)*v3;

    /* l>=2, m=odd*/
    /* spin-orbit */
    f21S = -1.5*a12*v + ((110./21. + 79./84.*nu)*a12 - 13./84.*a0X12)*v3;
    f33S = ((-0.25 + 2.5*nu)*a12 - 1.75*a0X12)*v3;
    f31S = ((-2.25 + 6.5*nu)*a12 + 0.25*a0X12)*v3;
    f43S = (( 5. -10.*nu)*a12 - 5.*a0X12)/(-4.+8.*nu)*v;
    f41S = f43S;

    /* SPIN-SPIN contribution */
    REAL8 c21SS_lo;
    REAL8 c33SS_lo;
    REAL8 c31SS_lo;
    if (usetidal) {
#if (RC_EXCLUDESPINSPINTIDES)
        /* Switch off spin-spin-tidal couplings */
        /* See also: eob_dyn_s_get_rc() */
        c21SS_lo  = 0.;
        c33SS_lo  = 0.;
        c31SS_lo  = 0.;
        /* Above code switch off everything,
         Alt. one can set C_Q1=C_Q2=0, but keep the term: */
#else
        c21SS_lo  = -19./8.*(a1*a1 - a2*a2) - (C_Q1*a1*a1 - C_Q2*a2*a2) + 1./8.*(-9.*a1*a1 + 10*a1*a2 -9.*a2*a2 + 12.*(C_Q1*a1*a1 + C_Q2*a2*a2))*X12;
        c33SS_lo  = 3.*(a1*a2 + 0.5*(C_Q1*a1*a1 + C_Q2*a2*a2))*X12;
        c31SS_lo  = -4.*(C_Q1*a1*a1 - C_Q2*a2*a2) + 3.*(a1*a2 + 0.5*(C_Q1*a1*a1 + C_Q2*a2*a2))*X12;
#endif
    } else {
        c21SS_lo  = 1./8.*(-27.*(a1*a1 - a2*a2) + (3.*a1*a1 + 10.*a1*a2 + 3.*a2*a2)*X12);
        c33SS_lo  = 3./2.*a0*a0*X12;
        c31SS_lo  = -4.*(a1*a1 - a2*a2) + 3./2.*a0*a0*X12;
    }

    /* Adding LO SS term w.r.t. eob_wav_flm_s_SSLO */
    f21S += c21SS_lo*v4;
    f33S += c33SS_lo*v4;
    f31S += c31SS_lo*v4;

    /* Amplitudes (correct with spin terms) */
    flm[0] = SQ(rholm[0]);
    flm[0] = (X12*flm[0] + f21S);

    flm[1] = SQ(rholm[1]+ rho22S);

    flm[2] = gsl_pow_int(rholm[2], 3);
    flm[2] = (X12*flm[2] + f31S);

    flm[3] = gsl_pow_int(rholm[3]+ rho32S, 3);

    flm[4] = gsl_pow_int(rholm[4], 3);
    flm[4] = (X12*flm[4] + f33S);

    flm[5] = gsl_pow_int(rholm[5], 4);
    flm[5] = (X12*flm[5] + f41S);

    flm[6] = gsl_pow_int(rholm[6] + rho42S, 4);

    flm[7] = gsl_pow_int(rholm[7], 4);
    flm[7] = (X12*flm[7] + f43S);

    flm[8] = gsl_pow_int(rholm[8] + rho44S, 4);

}

/* Pre-calculate coefficients for resummed amplitudes.
 *  Refs:
 *  . Damour, Iyer & Nagar, PRD 79, 064004 (2009)     [theory]
 *  . Fujita & Iyer, PRD 82, 044051 (2010)            [test-mass 5.5PN]
 *  . Damour, Nagar & Bernuzzi, PRD 87, 084035 (2013) [complete information]
 */
void eob_wav_flm_coeffs(REAL8 nu, REAL8 clm[KMAX][6])
{

    const REAL8 nu2 = nu*nu;
    const REAL8 nu3 = nu*nu2;
    const REAL8 nu4 = nu*nu3;

    for (int k=0; k<KMAX; k++) clm[k][0] = 1.;
    for (int k=0; k<KMAX; k++) for (int n=1; n<6; n++) clm[k][n] = 0.;

    /* (2,1) */
    clm[0][1] = (-1.0535714285714286 + 0.27380952380952384 *nu);
    clm[0][2] = (-0.8327841553287982 - 0.7789824263038548  *nu + 0.13116496598639457*nu2);
    /* clm[0][3] = (2.9192806270460925  - 1.019047619047619   *el1); */
    /* clm[0][4] = (-1.28235780892213   + 1.073639455782313   *el1); */
    /* clm[0][5] = (-3.8466571723355227 + 0.8486467106683944  *el1)*PMTERMS_eps; */

    /* (2,2) */
    clm[1][1] = (-1.0238095238095237 + 0.6547619047619048*nu);
    clm[1][2] = (-1.94208238851096   - 1.5601379440665155*nu + 0.4625614134542706*nu2);
    /* clm[1][3] = (12.736034731834051  - 2.902228713904598 *nu - 1.9301558466099282*nu2 + 0.2715020968103451*nu3 - 4.076190476190476*el2); */
    /* clm[1][4] = (-2.4172313935587004 + 4.173242630385488 *el2); */
    /* clm[1][5] = (-30.14143102836864  + 7.916297736025627 *el2); */

    /* (3,1) */
    clm[2][1] = (-0.7222222222222222 - 0.2222222222222222*nu);
    clm[2][2] = (0.014169472502805836 - 0.9455667789001122*nu - 0.46520763187429853*nu2);
    /* clm[2][3] = (1.9098284139598072 - 0.4126984126984127*el1); */
    /* clm[2][4] = (0.5368150316615179 + 0.2980599647266314*el1); */
    /* clm[2][5] = (1.4497991763035063 - 0.0058477188106817735*el1)*PMTERMS_eps; */

    /* (3,2) */
    clm[3][1] = (0.003703703703703704*(328. - 1115.*nu + 320.*nu2))/(-1. + 3.*nu);
    clm[3][2] = (6.235191420376606e-7*(-1.444528e6 + 8.050045e6*nu - 4.725605e6*nu2 - 2.033896e7*nu3 + 3.08564e6*nu4))/((-1. + 3.*nu)*(-1. + 3.*nu));
    /* clm[3][3] = (6.220997955214429 - 1.6507936507936507*el2); */
    /* clm[3][4] = (-3.4527288879001268 + 2.005408583186361*el2)*PMTERMS_eps; */

    /* (3,3) */
    clm[4][1] =  (-1.1666666666666667 + 0.6666666666666666*nu);
    clm[4][2] = (-1.6967171717171716 - 1.8797979797979798*nu + 0.45151515151515154*nu2);
    /* clm[4][3] = (14.10891386831863 - 3.7142857142857144*el3); */
    /* clm[4][4] = (-6.723375314944128 + 4.333333333333333*el3); */
    /* clm[4][5] = (-29.568699895427518 + 6.302092352092352*el3)*PMTERMS_eps; */

    /* (4,1) */
    clm[5][1] = (0.001893939393939394*(602. - 1385.*nu + 288.*nu2))/(-1. + 2.*nu);
    clm[5][2] = (- 0.36778992787515513);
    /* clm[5][3] = (0.6981550175535535 - 0.2266955266955267*el1); */
    /* clm[5][4] = (-0.7931524512893319 + 0.2584672482399755*el1)*PMTERMS_eps; */

    /* (4,2) */
    clm[6][1] = (0.0007575757575757576*(1146. - 3530.*nu + 285.*nu2))/(-1. + 3.*nu);
    clm[6][2] = - (3.1534122443213353e-9*(1.14859044e8 - 2.95834536e8*nu - 1.204388696e9*nu2 + 3.04798116e9*nu3 + 3.79526805e8*nu4))/((-1. + 3.*nu)*(-1. + 3.*nu));
    /* clm[6][3] = 4.550378418934105e-12*(8.48238724511e11 - 1.9927619712e11*el2); */
    /* clm[6][4] = (-0.6621921297263365 + 0.787251738160829*el2)*PMTERMS_eps; */

    /* (4,3) */
    clm[7][1] = (0.005681818181818182*(222. - 547.*nu + 160.*nu2))/(-1. + 2.*nu);
    clm[7][2] = (- 0.9783218202252293);
    /* clm[7][3] = (8.519456157072423 - 2.0402597402597404*el3)*PMTERMS_eps; */
    /* clm[7][4] = (-5.353216984886716 + 2.5735094451003544*el3)*PMTERMS_eps; */

    /* (4,4) */
    clm[8][1] = (0.0007575757575757576*(1614. - 5870.*nu + 2625.*nu2))/(-1. + 3.*nu);
    clm[8][2] = (3.1534122443213353e-9*(-5.11573572e8 + 2.338945704e9*nu - 3.13857376e8*nu2 - 6.733146e9*nu3 + 1.252563795e9*nu4))/((-1. + 3.*nu)*(-1. + 3.*nu));
    /* clm[8][3] = (15.108111214795123 - 3.627128427128427*el4); */
    /* clm[8][4] = (-8.857121657199649 + 4.434988849534304*el4)*PMTERMS_eps; */

    /* (5,1) */
    clm[9][1] = (0.002564102564102564*(319. - 626.*nu + 8.*nu2))/(-1. + 2.*nu);
    clm[9][2] = (- 0.1047896120973044);
    /* clm[9][3] = (0.642701885362399 - 0.14414918414918415*el1)*PMTERMS_eps; */
    /* clm[9][4] = (-0.07651588046467575 + 0.11790664036817883*el1)*PMTERMS_eps; */

    /* (5,2) */
    clm[10][1] = (0.00007326007326007326*(-15828. + 84679.*nu - 104930.*nu2 + 21980.*nu3))/(1. - 5.*nu + 5.*nu2);
    clm[10][2] = (- 0.4629337197600934)*PMTERMS_eps;
    /* clm[10][3] = (2.354458371550237 - 0.5765967365967366*el2)*PMTERMS_eps; */

    /* (5,3) */
    clm[11][1] = (0.002564102564102564*(375. - 850.*nu + 176.*nu2))/(-1. + 2.*nu);
    clm[11][2] = (- 0.5788010707241477);
    /* clm[11][3] = (5.733973288504755 - 1.2973426573426574*el3)*PMTERMS_eps; */
    /* clm[11][4] = (-1.9573287625526001 + 1.2474448628294783*el3)*PMTERMS_eps; */

    /* (5,4) */
    clm[12][1] = (0.00007326007326007326*(-17448. + 96019.*nu - 127610.*nu2 + 33320.*nu3))/(1. - 5.*nu + 5.*nu2);
    clm[12][2] = (- 1.0442142414362194)*PMTERMS_eps;
    /* clm[12][3] = (10.252052781721588 - 2.3063869463869464*el4)*PMTERMS_eps; */

    /* (5,5) */
    clm[13][1] = (0.002564102564102564*(487. - 1298.*nu + 512.*nu2))/(-1. + 2.*nu);
    clm[13][2] = (- 1.5749727622804546);
    /* clm[13][3] = (15.939827047208668 - 3.6037296037296036*el5)*PMTERMS_eps; */
    /* clm[13][4] = (-10.272578060123237 + 4.500041838503377*el5)*PMTERMS_eps; */

    /* (6,1) */
    clm[14][1] = (0.006944444444444444*(-161. + 694.*nu - 670.*nu2 + 124.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[14][2] = (- 0.29175486850885135)*PMTERMS_eps;
    /* clm[14][3] = (0.21653486654395454 - 0.10001110001110002*el1)*PMTERMS_eps; */

    /* (6,2) */
    clm[15][1] = (0.011904761904761904*(-74. + 378.*nu - 413.*nu2 + 49.*nu3))/(1. - 5.*nu + 5.*nu2);
    clm[15][2] = ( - 0.24797525070634313)*PMTERMS_eps;
    /* clm[15][3] = (1.7942694138754138 - 0.40004440004440006*el2)*PMTERMS_eps; */

    /* (6,3) */
    clm[16][1] = (0.006944444444444444*(-169. + 742.*nu - 750.*nu2 + 156.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[16][2] = (- 0.5605554442947213)*PMTERMS_eps;
    /* clm[16][3] = (4.002558222882566 - 0.9000999000999002*el3)*PMTERMS_eps; */

    /* (6,4) */
    clm[17][1] = (0.011904761904761904*(-86. + 462.*nu - 581.*nu2 + 133.*nu3))/(1. - 5.*nu + 5.*nu2);
    clm[17][2] = (- 0.7228451986855349)*PMTERMS_eps;
    /* clm[17][3] = (7.359388663371044 - 1.6001776001776002*el4)*PMTERMS_eps; */

    /* (6,5) */
    clm[18][1] = (0.006944444444444444*(-185. + 838.*nu - 910.*nu2 + 220.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[18][2] = (- 1.0973940686333457)*PMTERMS_eps;
    /* clm[18][3] = (11.623366217471297 - 2.5002775002775004*el5)*PMTERMS_eps; */

    /* (6,6) */
    clm[19][1] = (0.011904761904761904*(-106. + 602.*nu - 861.*nu2 + 273.*nu3))/(1. - 5.*nu + 5.*nu2);
    clm[19][2] = (- 1.5543111183867486)*PMTERMS_eps;
    /* clm[19][3] = (16.645950799433503 - 3.6003996003996006*el6)*PMTERMS_eps; */

    /* (7,1) */
    clm[20][1] = (0.0014005602240896359*(-618. + 2518.*nu - 2083.*nu2 + 228.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[20][2] = ( - 0.1508235111143767)*PMTERMS_eps;
    /* clm[20][3] = (0.2581280702019663 - 0.07355557607658449*el1)*PMTERMS_eps; */

    /* (7,2) */
    clm[21][1] = (0.00006669334400426837*(16832. - 123489.*nu + 273924.*nu2 - 190239.*nu3 + 32760.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[21][2] = (- 0.351319484450667)*PMTERMS_eps;

    /* (7,3) */
    clm[22][1] = (0.0014005602240896359*(-666. + 2806.*nu - 2563.*nu2 + 420.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[22][2] = (- 0.37187416047628863)*PMTERMS_eps;
    /* clm[22][3] = (3.0835293524055283 - 0.6620001846892604*el3)*PMTERMS_eps; */

    /* (7,4) */
    clm[23][1] = (0.00006669334400426837*(17756. - 131805.*nu + 298872.*nu2 - 217959.*nu3 + 41076.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[23][2] = (- 0.6473746896670599)*PMTERMS_eps;

    /* (7,5) */
    clm[24][1] = (0.0014005602240896359*(-762. + 3382.*nu - 3523.*nu2 + 804.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[24][2] = (- 0.8269193364414116)*PMTERMS_eps;
    /* clm[24][3] = (8.750589067052443 - 1.838889401914612*el5)*PMTERMS_eps; */

    /* (7,6) */
    clm[25][1] = (0.0006002400960384153*(2144. - 16185.*nu + 37828.*nu2 - 29351.*nu3 + 6104.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[25][2] = (- 1.1403265020692532)*PMTERMS_eps;

    /* (7,7) */
    clm[26][1] = (0.0014005602240896359*(-906. + 4246.*nu - 4963.*nu2 + 1380.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[26][2] = (- 1.5418467934923434)*PMTERMS_eps;
    /* clm[26][3] = (17.255875091408523 - 3.6042232277526396*el7)*PMTERMS_eps; */

    /* (8,1) */
    clm[27][1] = (0.00005482456140350877*(20022. - 126451.*nu + 236922.*nu2 - 138430.*nu3 + 21640.*nu4))/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
    clm[27][2] = (- 0.26842133517043704)*PMTERMS_eps;

    /* (8,2) */
    clm[28][1] = (0.0003654970760233918*(2462. - 17598.*nu + 37119.*nu2 - 22845.*nu3 + 3063.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[28][2] = (- 0.2261796441029474)*PMTERMS_eps;

    /* (8,3) */
    clm[29][1] = (0.00005482456140350877*(20598. - 131059.*nu + 249018.*nu2 - 149950.*nu3 + 24520.*nu4))/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
    clm[29][2] = (- 0.4196774909106648)*PMTERMS_eps;

    /* (8,4) */
    clm[30][1] = (0.0003654970760233918*(2666. - 19434.*nu + 42627.*nu2 - 28965.*nu3 + 4899.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[30][2] = (- 0.47652059150068155)*PMTERMS_eps;

    /* (8,5) */
    clm[31][1] = (0.00027412280701754384*(4350. - 28055.*nu + 54642.*nu2 - 34598.*nu3 + 6056.*nu4))/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
    clm[31][2] = (- 0.7220789990670207)*PMTERMS_eps;

    /* (8,6) */
    clm[32][1] = (0.0010964912280701754*(1002. - 7498.*nu + 17269.*nu2 - 13055.*nu3 + 2653.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[32][2] = (- 0.9061610303170207)*PMTERMS_eps;

    /* (8,7) */
    clm[33][1] = (0.00005482456140350877*(23478. - 154099.*nu + 309498.*nu2 - 207550.*nu3 + 38920.*nu4))/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
    clm[33][2] = (- 1.175404252991305)*PMTERMS_eps;

    /* (8,8) */
    clm[34][1] = (0.0003654970760233918*(3482. - 26778.*nu + 64659.*nu2 - 53445.*nu3 + 12243.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[34][2] = (- 1.5337092502821381)*PMTERMS_eps;

    return;
}


/* Resummed amplitudes in the general nu-dependent case.
 *  Refs:
 *  . Damour, Iyer & Nagar, PRD 79, 064004 (2009)     [theory]
 *  . Fujita & Iyer, PRD 82, 044051 (2010)            [test-mass 5.5PN]
 *  . Damour, Nagar & Bernuzzi, PRD 87, 084035 (2013) [complete information]
 */
void eob_wav_flm(REAL8 x, REAL8 nu, REAL8 clm[KMAX][6], REAL8 *rholm, REAL8 *flm)
{

    /* Coefficients */
    const REAL8 nu2 = nu*nu;
    const REAL8 nu3 = nu*nu2;
    const REAL8 Pi2 = SQ(LAL_PI);


    /* Compute EulerLogs */
    const REAL8 el1 = Eulerlog(x,1);
    const REAL8 el2 = Eulerlog(x,2);
    const REAL8 el3 = Eulerlog(x,3);
    const REAL8 el4 = Eulerlog(x,4);
    const REAL8 el5 = Eulerlog(x,5);
    const REAL8 el6 = Eulerlog(x,6);
    const REAL8 el7 = Eulerlog(x,7);

    /* Coefs with Eulerlogs */
    clm[0][3] = (2.9192806270460925  - 1.019047619047619   *el1);
    clm[0][4] = (-1.28235780892213   + 1.073639455782313   *el1);
    clm[0][5] = (-3.8466571723355227 + 0.8486467106683944  *el1)*PMTERMS_eps;

    clm[1][3] = (12.736034731834051  - 2.902228713904598 *nu - 1.9301558466099282*nu2 + 0.2715020968103451*nu3 - 4.076190476190476*el2);
    clm[1][4] = (-2.4172313935587004 + 4.173242630385488 *el2);
    clm[1][5] = (-30.14143102836864  + 7.916297736025627 *el2);

    clm[2][3] = (1.9098284139598072 - 0.4126984126984127*el1+ (-4.646868015386534 + (0.21354166666666666)*Pi2)*nu + 2.3020866307903347*nu2 - 0.5813492634480288*nu3);
    clm[2][4] = (0.5368150316615179 + 0.2980599647266314*el1);
    clm[2][5] = (1.4497991763035063 - 0.0058477188106817735*el1)*PMTERMS_eps;

    clm[3][3] = (6.220997955214429 - 1.6507936507936507*el2);
    clm[3][4] = (-3.4527288879001268 + 2.005408583186361*el2)*PMTERMS_eps;

    clm[4][3] = (14.10891386831863 - 3.7142857142857144*el3 + (-5.031429681429682 + (0.21354166666666666)*Pi2)*nu - 1.7781727531727531*nu2 + 0.25923767590434255*nu3);
    clm[4][4] = (-6.723375314944128 + 4.333333333333333*el3);
    clm[4][5] = (-29.568699895427518 + 6.302092352092352*el3)*PMTERMS_eps;

    clm[5][3] = (0.6981550175535535 - 0.2266955266955267*el1);
    clm[5][4] = (-0.7931524512893319 + 0.2584672482399755*el1)*PMTERMS_eps;

    clm[6][3] = 4.550378418934105e-12*(8.48238724511e11 - 1.9927619712e11*el2);
    clm[6][4] = (-0.6621921297263365 + 0.787251738160829*el2)*PMTERMS_eps;

    clm[7][3] = (8.519456157072423 - 2.0402597402597404*el3)*PMTERMS_eps;
    clm[7][4] = (-5.353216984886716 + 2.5735094451003544*el3)*PMTERMS_eps;

    clm[8][3] = (15.108111214795123 - 3.627128427128427*el4);
    clm[8][4] = (-8.857121657199649 + 4.434988849534304*el4)*PMTERMS_eps;

    clm[9][3] = (0.642701885362399 - 0.14414918414918415*el1)*PMTERMS_eps;
    clm[9][4] = (-0.07651588046467575 + 0.11790664036817883*el1)*PMTERMS_eps;

    clm[10][3] = (2.354458371550237 - 0.5765967365967366*el2)*PMTERMS_eps;

    clm[11][3] = (5.733973288504755 - 1.2973426573426574*el3)*PMTERMS_eps;
    clm[11][4] = (-1.9573287625526001 + 1.2474448628294783*el3)*PMTERMS_eps;

    clm[12][3] = (10.252052781721588 - 2.3063869463869464*el4)*PMTERMS_eps;

    clm[13][3] = (15.939827047208668 - 3.6037296037296036*el5)*PMTERMS_eps;
    clm[13][4] = (-10.272578060123237 + 4.500041838503377*el5)*PMTERMS_eps;

    clm[14][3] = (0.21653486654395454 - 0.10001110001110002*el1)*PMTERMS_eps;

    clm[15][3] = (1.7942694138754138 - 0.40004440004440006*el2)*PMTERMS_eps;

    clm[16][3] = (4.002558222882566 - 0.9000999000999002*el3)*PMTERMS_eps;

    clm[17][3] = (7.359388663371044 - 1.6001776001776002*el4)*PMTERMS_eps;

    clm[18][3] = (11.623366217471297 - 2.5002775002775004*el5)*PMTERMS_eps;

    clm[19][3] = (16.645950799433503 - 3.6003996003996006*el6)*PMTERMS_eps;

    clm[20][3] = (0.2581280702019663 - 0.07355557607658449*el1)*PMTERMS_eps;

    clm[22][3] = (3.0835293524055283 - 0.6620001846892604*el3)*PMTERMS_eps;

    clm[24][3] = (8.750589067052443 - 1.838889401914612*el5)*PMTERMS_eps;

    clm[26][3] = (17.255875091408523 - 3.6042232277526396*el7)*PMTERMS_eps;

    /* rho_lm */
    const REAL8 x2  = x*x;
    const REAL8 x3  = x*x2;
    const REAL8 x4  = x*x3;
    const REAL8 x5  = x*x4;
    const REAL8 xn[] = {1.,x,x2,x3,x4,x5};

    for (int k=0; k<KMAX; k++) {
        /* Note: the two sums give different result */
#if (1)
        rholm[k] = clm[k][0];
        for (int n=1; n<6; n++) {
            rholm[k] += clm[k][n]*xn[n];
        }
#else
        rholm[k] = x5*clm[k][5];
        for (int n=5; n-- >1; ) { // 4,3,2,1 //
            rholm[k] += clm[k][n]*xn[n];
        }
        rholm[k] += clm[k][0];
#endif
    }

//#pragma omp parallel
//#pragma omp for
    /* Amplitudes */
    for (int k = 0; k < KMAX; k++) {
        flm[k] = gsl_pow_int(rholm[k], TEOB_LINDEX[k]);
    }

    return;
}


/* Residual phase corrections \f$delta_{lm}\f$ up to \f$l=m=5\f$.
 Reference(s)
 Damour, Iyer & Nagar, PRD 79, 064004 (2008)
 Fujita & Iyer, PRD 82 044051 (2010)
 Faye et al., Class. Q. Grav. 29 175004 (2012)
 Damour, Nagar & Bernuzzi, PRD 87, 084035 (2013) */
//TODO: [optimization] this routine can be optimized: precompute coefficients c(nu)
void eob_wav_deltalm(REAL8 Hreal, REAL8 Omega, REAL8 nu, REAL8 *dlm)
{

    /* Useful shorthands*/
    const REAL8 Pi2 = SQ(LAL_PI);
    REAL8 nu2    = SQ(nu);
    REAL8 y      = cbrt(Hreal*Omega*Hreal*Omega);
    REAL8 sqrt_y = sqrt(y);
    REAL8 y3     = y*y*y;
    REAL8 y32    = Hreal*Omega;

    /* Leading order contributions*/
    REAL8 delta22LO = 7./3.   * y32;
    REAL8 delta21LO = 2./3.   * y32;
    REAL8 delta33LO = 13./10. * y32;
    REAL8 delta31LO = 13./30. * y32;

    /* Init phase */
    for (int k = 0; k < KMAX; k++) {
        dlm[k] = 0.;
    }

    /* Residual phases in Pade-resummed form when possible */
    REAL8 num;
    REAL8 den;

    /* l=2 */
    /* Pade(1,2) approximant */
    num        = 69020.*nu + 5992.*LAL_PI*sqrt_y;
    den        = 5992.*LAL_PI*sqrt_y + 2456.*nu*(28.+493.*nu* y);
    dlm[0] = delta21LO*num/den;
    /* Pade(2,2) approximant */
    num        = (808920.*nu*LAL_PI*sqrt(y) + 137388.*Pi2*y + 35.*nu2*(136080. + (154975. - 1359276.*nu)*y));
    den        = (808920.*nu*LAL_PI*sqrt(y) + 137388.*Pi2*y + 35.*nu2*(136080. + (154975. + 40404.*nu)*y));
    dlm[1] = delta22LO*num/den;

    /* l=3 */
    /* Pade(1,2) approximant */
    num        = 4641.*nu + 1690.*LAL_PI*sqrt_y;
    den        = num + 18207.*nu2*y;
    dlm[2] = delta31LO*num/den;
    /* Taylor-expanded form */
    num        = 1.  + 94770.*LAL_PI/(566279.*nu)*sqrt_y;
    den        = num + 80897.* nu/3159.*y;
    dlm[3] = (10.+33.*nu)/(15.*(1.-3.*nu)) * y32 + 52./21.*LAL_PI*y3;
    /* Pade(1,2) approximant */
    dlm[4] = delta33LO*num/den;

    /* l=4 */
    dlm[5] =   (2.+507.*nu)/(10.*(1.-2.*nu))*y32   + 1571./3465.*LAL_PI*y3;
    dlm[6] =  7.*(1.+6.*nu)/(15.*(1.-3.*nu))*y32   + 6284./3465.*LAL_PI*y3;
    dlm[7] = (486.+4961.*nu)/(810.*(1.-2.*nu))*y32 + 1571./385.*LAL_PI*y3;
    dlm[8] =  (112.+219.*nu)/(120.*(1.-3.*nu))*y32 + 25136./3465.*LAL_PI*y3;

    /* l=5 */
    dlm[9] = (96875. + 857528.*nu)/(131250.*(1.-2.*nu))*y32;

    return;
}


/* Leading-order (Newtonian) prefactor  of the multipolar resummed waveform.
 Reference: Damour, Iyer & Nagar, PRD 79, 064004 (2009) */
void eob_wav_hlmNewt(REAL8 r,
                     REAL8 Omega,
                     REAL8 phi,
                     REAL8 nu,
                     LALTEOBResumSWaveformModeSingleTime *hlmNewt)
{
    /* Shorthands */
    REAL8 nu2   = nu*nu;
    REAL8 nu3   = nu*nu2;

    REAL8 vphi  = r*Omega;
    REAL8 vphi2 = vphi*vphi;
    REAL8 vphi3 = vphi*vphi2;
    REAL8 vphi4 = vphi*vphi3;
    REAL8 vphi5 = vphi*vphi4;
    REAL8 vphi6 = vphi*vphi5;
    REAL8 vphi7 = vphi*vphi6;
    REAL8 vphi8 = vphi*vphi7;
    REAL8 vphi9 = vphi*vphi8;

    /* Polynomials in nu */
    const REAL8 p1 = 1.;
    const REAL8 p2 = sqrt(1.-4.*nu);
    const REAL8 p3 = (3.*nu-1.);
    const REAL8 p4 = (2.*nu-1.)*sqrt(1.-4.*nu);
    const REAL8 p5 = 1.-5.*nu+5.*nu2;
    const REAL8 p6 = (1.-4.*nu+3.*nu2)*sqrt(1.-4.*nu);
    const REAL8 p7 = 7.*nu3 - 14.*nu2 + 7.*nu -1.;
    const REAL8 p8 = (4.*nu3 - 10.*nu2 + 6.*nu -1.)*sqrt(1.-4.*nu);

    const REAL8 phix2 = 2. * phi;
    const REAL8 phix3 = 3. * phi;
    const REAL8 phix4 = 4. * phi;
    const REAL8 phix5 = 5. * phi;
    const REAL8 phix6 = 6. * phi;
    const REAL8 phix7 = 7. * phi;

    const REAL8 pv23 = p2 * vphi3;
    const REAL8 pv34 = p3 * vphi4;
    const REAL8 pv45 = p4 * vphi5;
    const REAL8 pv56 = p5 * vphi6;
    const REAL8 pv67 = p6 * vphi7;
    const REAL8 pv78 = p7 * vphi8;
    const REAL8 pv89 = p8 * vphi9;

    REAL8 phim[KMAX] = {
        phi,phix2,
        phi,phix2,phix3,
        phi,phix2,phix3,phix4,
        phi,phix2,phix3,phix4,phix5,
        phi,phix2,phix3,phix4,phix5,phix6,
        phi,phix2,phix3,phix4,phix5,phix6,phix7,
        phi,phix2,phix3,phix4,phix5,phix6,phix7,8.*phi
    };

    REAL8 Alm[KMAX] = {
        pv23, p1 * vphi2,
        pv23, pv34, pv23,
        pv45, pv34, pv45, pv34,
        pv45, pv56, pv45, pv56, pv45,
        pv67, pv56, pv67, pv56, pv67, pv56,
        pv67, pv78, pv67, pv78, pv67, pv78, pv67,
        pv89, pv78, pv89, pv78, pv89, pv78, pv89, pv78
    };

    /* Compute hlmNewt (without phase factor) in complex Polar coords */
    for (int k = 0; k < KMAX; k++) {
        hlmNewt->phase[k] = - phim[k] + ChlmNewt_phase[k];
        hlmNewt->ampli[k] = ChlmNewt_ampli[k] * Alm[k];
    }

}


/* Calculate tidal correction to multipolar waveform amplitude
 Ref. Damour, Nagar & Villain, Phys.Rev. D85 (2012) 123007 */
void eob_wav_hlmTidal(REAL8 x, LALTEOBResumSDynamics *dyn, REAL8 *hTidallm)
{
    const REAL8 nu       = dyn->nu;
    const REAL8 XA       = dyn->X1;
    const REAL8 XB       = dyn->X2;
    const REAL8 khatA_2  = dyn->khatA2;
    const REAL8 khatB_2  = dyn->khatB2;
    const REAL8 kapA2j   = dyn->kapA2j;
    const REAL8 kapB2j   = dyn->kapB2j;
    const REAL8 kapT2j   = dyn->kapT2j;

    const REAL8 x5 = (REAL8) gsl_pow_int((double) x, 5);
    const REAL8 x6 = (REAL8) gsl_pow_int((double) x, 6);

    REAL8 hA[KMAX], hB[KMAX], betaA1[KMAX], betaB1[KMAX];

    memset(hTidallm, 0., KMAX*sizeof(REAL8));
    memset(hA, 0., KMAX*sizeof(REAL8));
    memset(hB, 0., KMAX*sizeof(REAL8));
    memset(betaA1, 0., KMAX*sizeof(REAL8));
    memset(betaB1, 0., KMAX*sizeof(REAL8));

    /* l=2 */
    hA[1]     = 2 * khatA_2 *(XA/XB+3);
    hB[1]     = 2 * khatB_2 *(XB/XA+3);

    betaA1[1] = (-202. + 560*XA - 340*XA*XA + 45*XA*XA*XA)/(42*(3-2*XA));
    betaB1[1] = (-202. + 560*XB - 340*XB*XB + 45*XB*XB*XB)/(42*(3-2*XB));

    hA[0]     = 3 * khatA_2 * (3-4*XA);
    hB[0]     = 3 * khatB_2 * (3-4*XB);

    /* l=3 */
    hA[2] = 12 * khatA_2 * XB;
    hB[2] = 12 * khatB_2 * XA;

    betaA1[2] = (-6. -5.*XA +131.*XA*XA -130.*XA*XA*XA)/(36.*(1.-XA));
    betaB1[2] = (-6. -5.*XB +131.*XB*XB -130.*XB*XB*XB)/(36.*(1.-XB));

    hA[4] = hA[2];
    hB[4] = hB[2];

    betaA1[4] = ( (XA-3.)*(10.*XA*XA - 25.*XA+ 14.) )/(12.*(1.-XA));
    betaB1[4] = ( (XB-3.)*(10.*XB*XB - 25.*XB+ 14.) )/(12.*(1.-XB));

    /* l=2 */
    /* (2,1) */
    hTidallm[0] = ( -hA[0] + hB[0] )*x5;
    /* (2,2) */
    hTidallm[1] = ( hA[1]*(1. + betaA1[1]*x) + hB[1]*(1. + betaB1[1]*x) )*x5;

    /* l=3 */
    /* (3,1) */
    hTidallm[2] = ( -hA[2]*(1. + betaA1[2]*x) + hB[2]*(1. + betaB1[2]*x) )*x5;
    /* (3,2) */
    hTidallm[3] = 8.*( khatA_2*(1. -2.*XB + 3.*XB*XB) +khatB_2*(1. -2.*XA + 3.*XA*XA) )*x5/(1.-3.*nu);
    /* (3,3) */
    hTidallm[4] = ( -hA[4]*(1. + betaA1[4]*x) + hB[4]*(1. + betaB1[4]*x) )*x5;
    // MA: TODO: replace


    if ( (dyn->use_tidal_gravitomagnetic==LAL_SIM_INSPIRAL_GMTIDES_GSF) || (dyn->use_tidal_gravitomagnetic==LAL_SIM_INSPIRAL_GMTIDES_PN) ) {
        const double fourtnine= 1.5555555555555555556;  // 14/9 = 112/(3*24)
        const double fourthird = 1.3333333333333333333; // 32/24 = 4/3
        hTidallm[0] += 0.5*( -1.*kapA2j/XB + kapB2j/XA )*x5;
        hTidallm[1] += fourtnine*kapT2j*x6;
        hTidallm[2] += 0.5*( kapA2j*(4. - 17.*XB) - kapB2j*(4. - 17.*XA) )*x6;
        hTidallm[3] += fourthird*kapT2j*x5/(1.-3.*nu);
        hTidallm[4] += 0.5*( kapA2j*(4. - 9.*XB) - kapB2j*(4. - 9.*XA) )*x6;
    }

}

/* NQC corrections to the RWZ multipolar waveform
 Nagar, Damour, Reisswig, Pollney http://arxiv.org/abs/1506.08457
 Nonspinning case, Current fits: 9/02/2016 */
void eob_wav_hlmNQC_nospin201602(REAL8  nu,
                                 REAL8  r,
                                 REAL8  prstar,
                                 REAL8  Omega,
                                 REAL8  ddotr,
                                 LALTEOBResumSWaveformModeSingleTime *hlmnqc)
{
    const REAL8 xnu  = 1-4*nu;
    const REAL8 xnu2 = SQ(xnu);

    REAL8 a1[KMAX], a2[KMAX], a3[KMAX];
    REAL8 b1[KMAX], b2[KMAX], b3[KMAX];
    REAL8 n[KMAX][6];

    const INT4 k21 = 0;
    const INT4 k22 = 1;
    const INT4 k33 = 4;
    int k;

    /* NR fits */
    for (int ki = 0; ki < KMAX; ki++) {
        a1[ki] = 0.;
        a2[ki] = 0.;
        a3[ki] = 0.;
        b1[ki] = 0.;
        b2[ki] = 0.;
        b3[ki] = 0.;
    }

    /* (2,1) */
    a1[k21] =  0.0162387198*(7.32653082*xnu2 + 1.19616248*xnu + 0.73496656);
    a2[k21] =                -1.80492460*xnu2 + 1.78172686*xnu + 0.30865284;
    a3[k21] =                                                           0.0;

    b1[k21] =  -0.0647955017*(3.59934444*xnu2 - 4.08628784*xnu + 1.37890907);
    b2[k21] =   1.3410693180*(0.38491989*xnu2 + 0.10969453*xnu + 0.97513971);
    b3[k21] =                                                            0.0;

    /* (2,2) */
    a1[k22]   = -0.0805236959*( 1 - 2.00332326*xnu2)/( 1 + 3.08595088*xnu2);
    a2[k22]   =  1.5299534255*( 1 + 1.16438929*xnu2)/( 1 + 1.92033923*xnu2);
    a3[k22]   =  0.0;

    b1[k22]   = 0.146768094955*( 0.07417121*xnu + 1.01691256);
    b2[k22]   = 0.896911234248*(-0.61072011*xnu + 0.94295129);
    b3[k22]   = 0.0;

    /* (3,3) */
    a1[k33]   = -0.0377680000*(1 - 14.61548907*xnu2)/( 1 + 2.44559263*xnu2);
    a2[k33]   =  1.9898000000*(1 + 2.09750346 *xnu2)/( 1 + 2.57489466*xnu2);
    a3[k33]   =  0.0;

    b1[k33]   = 0.1418400000*(1.07430512 - 1.23906804*xnu + 4.44910652*xnu2);
    b2[k33]   = 0.6191300000*(0.80672432 + 4.07432829*xnu - 7.47270977*xnu2);
    b3[k33]   = 0.0;

    /* NQC corrections to the modulus and phase */
    for (int ki = 0; ki < KMAX; ki++) {
        for (int j = 0; j < 6; j++) {
            n[ki][j] = 0.;
        }
    }

    k = k21;
    n[k][0] = (prstar/(r*Omega))*(prstar/(r*Omega));
    n[k][1] = ddotr/(r*Omega*Omega);
    n[k][2] = n[k][0]*prstar*prstar;
    n[k][3] = prstar/(r*Omega);
    n[k][4] = n[k][3]*cbrt(Omega*Omega);
    n[k][5] = n[k][4]*prstar*prstar;

    k = k22;
    n[k][0] = (prstar/(r*Omega))*(prstar/(r*Omega));
    n[k][1] = ddotr/(r*Omega*Omega);
    n[k][2] = n[k][0]*prstar*prstar;
    n[k][3] = prstar/(r*Omega);
    /* n[k][4] = n[k][3]*cbrt(Omega*Omega); */
    n[k][4] = n[k][3]*(r*Omega)*(r*Omega);
    n[k][5] = n[k][4]*prstar*prstar;

    k = k33;
    n[k][0] = (prstar/(r*Omega))*(prstar/(r*Omega));
    n[k][1] = ddotr/(r*Omega*Omega);
    n[k][2] = n[k][0]*prstar*prstar;
    n[k][3] = prstar/(r*Omega);
    n[k][4] = n[k][3]*cbrt(Omega*Omega);
    n[k][5] = n[k][4]*prstar*prstar;

    /* NQC factor */
    for (int ki = 0; ki < KMAX; ki++) {
        hlmnqc->ampli[ki] = 1.;
        hlmnqc->phase[ki] = 0.;
    }

    k = k21; /* (2,1) */
    hlmnqc->ampli[k] = 1. + a1[k]*n[k][0] + a2[k]*n[k][1] + a3[k]*n[k][2];
    hlmnqc->phase[k] =      b1[k]*n[k][3] + b2[k]*n[k][4] + b3[k]*n[k][5];

    k = k22; /* (2,2) */
    hlmnqc->ampli[k] = 1. + a1[k]*n[k][0] + a2[k]*n[k][1] + a3[k]*n[k][2];
    hlmnqc->phase[k] =      b1[k]*n[k][3] + b2[k]*n[k][4] + b3[k]*n[k][5];

    k = k33; /* (3,3) */
    hlmnqc->ampli[k] = 1. + a1[k]*n[k][0] + a2[k]*n[k][1] + a3[k]*n[k][2];
    hlmnqc->phase[k] =      b1[k]*n[k][3] + b2[k]*n[k][4] + b3[k]*n[k][5];

    return;
}

/* Generic routine for NQC */
void eob_wav_hlmNQC(REAL8 UNUSED nu, REAL8 r, REAL8 prstar, REAL8 Omega, REAL8 ddotr, NQCcoefs *nqc,
                    LALTEOBResumSWaveformModeSingleTime *hlmnqc)
{
    const INT4 maxk = MIN(KMAX, nqc->maxk+1);

    /* Multipoles with special treatment */
    const INT4 k22 = 1;

    /* Shorthand */
    const REAL8 n0 = (prstar/(r*Omega))*(prstar/(r*Omega));
    const REAL8 n1 = ddotr/(r*Omega*Omega);
    const REAL8 n2 = n0*SQ(prstar);
    const REAL8 n3 = prstar/(r*Omega);
    const REAL8 n4 = n3*cbrt(Omega*Omega);
    const REAL8 n5 = n4*SQ(prstar);

    const REAL8 n4_k = n3*SQ((r*Omega));
    const REAL8 n5_k = n4*SQ(prstar);

    /* n functions */
    for (int k = 0; k < maxk; k++) {
        if (nqc->activemode[k]) {
            nqc->n[k][0] = n0;
            nqc->n[k][1] = n1;
            nqc->n[k][2] = n2;
            nqc->n[k][3] = n3;
            nqc->n[k][4] = n4;
            nqc->n[k][5] = n5;
        }
    }

    /* Change special multipoles */
    INT4 k0 = k22;
    nqc->n[k0][4] = n4_k;
    nqc->n[k0][5] = n5_k;

    /* NQC wave factor */
    for (int k = 0; k < KMAX; k++) {
        hlmnqc->ampli[k] = 1.;
        hlmnqc->phase[k] = 0.;
    }

    for (int k = 0; k < maxk; k++) {
        if (nqc->activemode[k]) {
            hlmnqc->ampli[k] += nqc->a1[k]*nqc->n[k][0] + nqc->a2[k]*nqc->n[k][1] + nqc->a3[k]*nqc->n[k][2];
            hlmnqc->phase[k] += nqc->b1[k]*nqc->n[k][3] + nqc->b2[k]*nqc->n[k][4] + nqc->b3[k]*nqc->n[k][5];
        }
    }

}

/* Computes the factors and the coefficients that build the
 NQC corrections to the waveform in the spinning case */
void eob_wav_hlmNQC_find_a1a2a3(LALTEOBResumSDynamics *dyn, SphHarmPolarTimeSeries *h, SphHarmPolarTimeSeries *hnqc)
{

    REAL8 A_tmp, dA_tmp, omg_tmp, domg_tmp;

    const REAL8 nu   = dyn->nu;
    const REAL8 X1   = dyn->X1;
    const REAL8 X2   = dyn->X2;
    const REAL8 chi1 = dyn->chi1;
    const REAL8 chi2 = dyn->chi2;
    const REAL8 aK   = dyn->a1 + dyn->a2;

    const REAL8 nu2  = SQ(nu);
    const REAL8 nu3  = nu2*nu;
    const REAL8 X12  = X1 - X2;
    const REAL8 aK2  = SQ(aK);
    const REAL8 aK3  = aK2*aK;
    const REAL8 aK4  = aK2*aK2;
    const REAL8 a12  = X1*chi1 - X2*chi2;
    const REAL8 aeff     = aK + 1./3.*a12*X12;
    const REAL8 aeff_omg = aK + a12*X12;

    REAL8 *t       = h->tdata->data;
    REAL8 *r       = dyn->data[TEOB_RAD];
    REAL8 *w       = dyn->data[TEOB_MOMG]; /* Omega */
    REAL8 UNUSED *pph     = dyn->data[TEOB_PPHI];
    REAL8 *pr_star = dyn->data[TEOB_PRSTAR];
    REAL8 *Omg_orb = dyn->data[TEOB_OMGORB]; /* Omega orbital */
    REAL8 *ddotr   = dyn->data[TEOB_DDOTR];

    REAL8 c_p1,     c_p2,     c_p3,   c_p4;
    REAL8 c_pdA1,   c_pdA2,   c_pdA3, c_pdA4;
    REAL8 c_pdomg1, c_pdomg2;
    REAL8 n0, d1;
    REAL8 a0_omg_tmp, a1_omg_tmp, a2_omg_tmp, b0_omg_tmp, b1_omg_tmp, b2_omg_tmp, a0_domg_tmp, a1_domg_tmp, a2_domg_tmp, b0_domg_tmp, b1_domg_tmp, b2_domg_tmp, a0_A_tmp, a1_A_tmp , a2_A_tmp, b0_A_tmp, b1_A_tmp, b2_A_tmp, a0_dA_tmp, a1_dA_tmp, a2_dA_tmp, b0_dA_tmp, b1_dA_tmp, b2_dA_tmp, omg_tmp_nu, omg_tmp_equal, domg_tmp_nu, domg_tmp_equal,  A_tmp_scale_nu, A_tmp_scale_equal, dA_tmp_scale_nu, dA_tmp_scale_equal ;

    REAL8 P[2], M[4], p1[2], p2[2], p3[2], p4[2], pA[5], pdA[5];
    REAL8 pomg[5], pdomg[5], pn0[2], pd1[2], ppdomg1[2], ppdomg2[2], pdA1[2],pdA2[2],pdA3[2],pdA4[2];
    REAL8 max_A[KMAX], max_dA[KMAX], UNUSED d2max[KMAX], UNUSED d3max[KMAX], max_omg[KMAX], max_domg[KMAX], UNUSED maxd2omg[KMAX], UNUSED DeltaT[KMAX];
    REAL8 ai[KMAX][2];
    REAL8 bi[KMAX][2];

    const int size = (int) h->tdata->length;
    for (int i = 0; i < size; i++) {
        hnqc->tdata->data[i] = t[i];
    }

    REAL8 *omg[KMAX], *domg[KMAX];
    REAL8 *n1,*n2,*n3,*n4,*n5,*n6, *d_n4,*d_n5, *d2_n4,*d2_n5;
    REAL8 *m11[KMAX], *m12[KMAX], *m13[KMAX], *m21[KMAX], *m22[KMAX];
    REAL8 *p1tmp[KMAX], *p2tmp[KMAX]; /* RWZ amplitude and derivative */

    for (int k=0; k<KMAX; k++) {
        omg[k]  = (REAL8*) calloc (size,sizeof(REAL8));
        domg[k] = (REAL8*) calloc (size,sizeof(REAL8));
        m11[k] = (REAL8*) calloc (size,sizeof(REAL8));
        m12[k] = (REAL8*) calloc (size,sizeof(REAL8));
        m13[k] = (REAL8*) calloc (size,sizeof(REAL8));
        m21[k] = (REAL8*) calloc (size,sizeof(REAL8));
        m22[k] = (REAL8*) calloc (size,sizeof(REAL8));
        p1tmp[k] = (REAL8*) calloc (size,sizeof(REAL8));
        p2tmp[k] = (REAL8*) calloc (size,sizeof(REAL8));
    }

    n1 = (REAL8*) calloc (size,sizeof(REAL8));
    n2 = (REAL8*) calloc (size,sizeof(REAL8));
    n3 = (REAL8*) calloc (size,sizeof(REAL8));
    n4 = (REAL8*) calloc (size,sizeof(REAL8));
    n5 = (REAL8*) calloc (size,sizeof(REAL8));
    n6 = (REAL8*) calloc (size,sizeof(REAL8));
    d_n4 = (REAL8*) calloc (size,sizeof(REAL8));
    d_n5 = (REAL8*) calloc (size,sizeof(REAL8));
    d2_n4 = (REAL8*) calloc (size,sizeof(REAL8));
    d2_n5 = (REAL8*) calloc (size,sizeof(REAL8));

    /* omega derivatives */
    const REAL8 dt = t[1]-t[0];
    SphHarmPolarTimeSeries *this;
    this = h;
    for (int k=0; k<KMAX; k++) {
        D0(this->phase->data->data, dt, size, omg[k]);
        D0(omg[k], dt, size, domg[k]);
        this = this->next;
    }

    /* NR fits */
    if (DEQUAL(nu,0.25,1e-9)) {

        pA[0]    =  0.00178195;
        pA[1]    =  0.00435589;
        pA[2]    =  0.00344489;
        pA[3]    = -0.00076165;
        pA[4]    =  0.31973334;
        A_tmp    =  pA[0]*aK4    + pA[1]*aK3   + pA[2]*aK2    + pA[3]*aK     + pA[4];

        pdA[0]   =  0.00000927;
        pdA[1]   = -0.00024550;
        pdA[2]   =  0.00012469;
        pdA[3]   =  0.00123845;
        pdA[4]   = -0.00195014;
        dA_tmp   =  pdA[0]*aK4   + pdA[1]*aK3   + pdA[2]*aK2   + pdA[3]*aK   + pdA[4];

        pomg[0]  =  0.00603482;
        pomg[1]  =  0.01604555;
        pomg[2]  =  0.02290799;
        pomg[3]  =  0.07084587;
        pomg[4]  =  0.38321834;
        omg_tmp  =  pomg[0]*aK4  + pomg[1]*aK3  + pomg[2]*aK2  + pomg[3]*aK  + pomg[4];

        pdomg[0] =  0.00024066;
        pdomg[1] =  0.00038123;
        pdomg[2] = -0.00049714;
        pdomg[3] =  0.00041219;
        pdomg[4] =  0.01190548;
        domg_tmp =  pdomg[0]*aK4 + pdomg[1]*aK3 + pdomg[2]*aK2 + pdomg[3]*aK + pdomg[4];

    }  else if( nu > 0.16) {

        p1[0]      =  0.04680896;
        p1[1]      = -0.00632114;
        p2[0]      =  0.06586192;
        p2[1]      = -0.01180039;
        p3[0]      = -0.11617413;
        p3[1]      =  0.02704959;
        p4[0]      =  0.15597465;
        p4[1]      =  0.28034978;
        c_p1       =  p1[0]*nu + p1[1];
        c_p2       =  p2[0]*nu + p2[1];
        c_p3       =  p3[0]*nu + p3[1];
        c_p4       =  p4[0]*nu + p4[1];
        A_tmp      =  c_p1*aK3 + c_p2*aK2 + c_p3*aK + c_p4;

        pdA1[0]    = -0.00130824;
        pdA1[1]    =  0.00006202;
        pdA2[0]    =  0.00199855;
        pdA2[1]    = -0.00027474;
        pdA3[0]    =  0.00218838;
        pdA3[1]    =  0.00071540;
        pdA4[0]    = -0.00362779;
        pdA4[1]    = -0.00105397;
        c_pdA1     =  pdA1[0]*nu + pdA1[1];
        c_pdA2     =  pdA2[0]*nu + pdA2[1];
        c_pdA3     =  pdA3[0]*nu + pdA3[1];
        c_pdA4     =  pdA4[0]*nu + pdA4[1];
        dA_tmp     =  c_pdA1*aK3   + c_pdA2*aK2 + c_pdA3*aK+ c_pdA4;

        pn0[0]     =  0.46908067;
        pn0[1]     =  0.27022141;
        pd1[0]     =  0.64131115;
        pd1[1]     = -0.37878384;
        n0         =  pn0[0]*nu + pn0[1];
        d1         =  pd1[0]*nu + pd1[1];
        omg_tmp    =  n0/(1 + d1*aK);

        ppdomg1[0] =  0.00061175;
        ppdomg1[1] =  0.00074001;
        ppdomg2[0] =  0.02504442;
        ppdomg2[1] =  0.00548217;
        c_pdomg1   =  ppdomg1[0]*nu + ppdomg1[1];
        c_pdomg2   =  ppdomg2[0]*nu + ppdomg2[1];
        domg_tmp   =  c_pdomg1*aK   + c_pdomg2;

    }  else {

        /* Fit by G.Riemanschneider incorporating the test-particle NQC point
         obtained from the most-recent Teukolsky waveforms done by
         M. Colleoni using the 6PN-accurare iResum-radiation reaction.
         These points assure a smooth connection between merger and
         ringdown also outside the "calibration" domain, notably for
         large-mass ratios (though q<=20) and large (negative) spins
         Updated, 28/09/2017 */

        a0_omg_tmp    = -0.1460961247;
        a1_omg_tmp    =  0.0998056;
        a2_omg_tmp    = -0.118098;
        b0_omg_tmp    = -0.3430184009;
        b1_omg_tmp    =  0.0921551;
        b2_omg_tmp    = -0.0740285;
        omg_tmp_nu    = +0.5427169903*nu2 +0.2512395608*nu +0.2863992248;
        omg_tmp_equal =((a2_omg_tmp*X12*X12 + a1_omg_tmp*X12 + a0_omg_tmp)*aeff_omg+1)/((b2_omg_tmp*X12*X12 +b1_omg_tmp*X12 + b0_omg_tmp)*aeff_omg+1);
        omg_tmp       = omg_tmp_nu*omg_tmp_equal;

        a0_domg_tmp    = +0.0604556289;
        b0_domg_tmp    = -0.0299583285;
        a1_domg_tmp    = 0.0711715;
        a2_domg_tmp    = -0.0500886;
        b1_domg_tmp    = 0.0461239;
        b2_domg_tmp    = -0.0153068;

        domg_tmp_nu    = ( +0.0045213831*nu +0.0064934920)/( -1.4466409969*nu+1);
        domg_tmp_equal = (a2_domg_tmp*X12*X12 +a1_domg_tmp*X12 +b0_domg_tmp)*aeff_omg*aeff_omg +(b2_domg_tmp*X12*X12 +b1_domg_tmp*X12+a0_domg_tmp)*aeff_omg+1;
        domg_tmp       = domg_tmp_nu*domg_tmp_equal;

        a0_A_tmp     = -0.2750516062;
        b0_A_tmp     = -0.4693776065;
        a1_A_tmp     =  0.143066;
        a2_A_tmp     = -0.0425947;
        b1_A_tmp     =  0.176955;
        b2_A_tmp     = -0.111902;

        A_tmp_scale_nu    = -0.9862040409*nu3 +0.8167558040*nu2 -0.0427442282*nu+0.2948879452;
        A_tmp_scale_equal = ((a2_A_tmp*X12*X12 + a1_A_tmp*X12 +a0_A_tmp)*aeff+1)/((b2_A_tmp*X12*X12 + b1_A_tmp*X12 +b0_A_tmp)*aeff+1);
        A_tmp             = A_tmp_scale_nu*A_tmp_scale_equal*(1-0.5*omg_tmp*aeff);

        a0_dA_tmp     = +0.0037461628;
        b0_dA_tmp     = +0.0636082543;
        a1_dA_tmp     =  0.00129393;
        a2_dA_tmp     = -0.00239069;
        b1_dA_tmp     = -0.0534209;
        b2_dA_tmp     = -0.186101;

        dA_tmp_scale_nu    = ( -0.0847947167*nu -0.0042142765)/( +16.1559461812*nu+1);
        dA_tmp_scale_equal = ((a2_dA_tmp*X12*X12 + a1_dA_tmp*X12+ a0_dA_tmp)*aeff)/((b2_dA_tmp*X12*X12 + b1_dA_tmp*X12 + b0_dA_tmp)*aeff+1);
        dA_tmp             = (dA_tmp_scale_nu +dA_tmp_scale_equal)*omg_tmp;

    }

    /* Switch on the 22 values (only) */
    for (int k=0; k<KMAX; k++) {
        max_A[k]    = 0.;
        max_dA[k]   = 0.;
        max_omg[k]  = 0.;
        max_domg[k] = 0.;
    }
    max_A[1]    = A_tmp;
    max_dA[1]   = dA_tmp;
    max_omg[1]  = omg_tmp;
    max_domg[1] = domg_tmp;

    //    if (TEOB_VERBOSE) {
    //        printf("NR values for NQC determination:\n");
    //        PRFORMd("A22_mrg",max_A[1]);
    //        PRFORMd("dA22_mrg",max_dA[1]);
    //        PRFORMd("omg22_mrg",max_omg[1]);
    //        PRFORMd("domg22_mrg",max_domg[1]);
    //    }

    /* NQC corrections to AMPLITUDE (n1,n2,n3) and PHASE (n4,n5,n6)
     * NQC basis for (2,2) waveform : AMPLITUDE
     * note: n3 and n6 are not used
     */
    REAL8 pr_star2, r2, w2;
    for (int j=0; j<size; j++) {
        pr_star2 = SQ(pr_star[j]);
        r2       = SQ(r[j]);
        w2       = SQ(w[j]); //CHECKME: Omg or Omg_orbital ?
        n1[j]  = pr_star2/(r2*w2);         /* [pr*\/(r Omg)]^2 */
        n2[j]  = ddotr[j]/(r[j]*w2);       /* [ddot{r}/(r Omg^2)] */
        n4[j]  = pr_star[j]/(r[j]*w[j]);   /* pr*\/(r Omg) */
        n5[j]  = n4[j]*r2*w2;              /* (pr*)*(r Omg) */
    }

    //#if (TEOB_DEBUG)
    //    FILE* fp = fopen("nqc_nfunc.txt", "w");
    //    for (int j=0; j<size; j++) {
    //        fprintf(fp, "%20.12f\t%.16e\t%.16e\t%.16e\t%.16e\n", t[j], n1[j], n2[j], n4[j], n5[j]);
    //    }
    //    fclose(fp);
    //#endif

    /* Derivatives for the phase */
    D0(n4,dt,size, d_n4);
    D0(n5,dt,size, d_n5);
    D0(d_n4,dt,size, d2_n4);
    D0(d_n5,dt,size, d2_n5);

    //#if (TEOB_DEBUG)
    //    fp = fopen("nqc_dfunc.txt", "w");
    //    for (int j=0; j<size; j++) {
    //        fprintf(fp, "%f\t%.16e\t%.16e\t%.16e\t%.16e\n", t[j], d_n4[j], d_n5[j], d2_n4[j], d2_n5[j]);
    //    }
    //    fclose(fp);
    //#endif

    /* Find max Omg */
    //TODO: search backwards!
    int Omgmax_index = 0;
    REAL8 Omg_max   = Omg_orb[0];
    for (int j=0; j<size; j++) {
        if (Omg_orb[j] > Omg_max) {
            Omg_max = Omg_orb[j];
            Omgmax_index = j;
        }
    }

    /* Time */
    REAL8 tOmgOrb_pk = t[Omgmax_index];
    REAL8 DeltaT_nqc = eob_nqc_timeshift(nu, chi1);
    REAL8 tNQC = tOmgOrb_pk - DeltaT_nqc;

    //    if (TEOB_VERBOSE) {
    //        printf("NQC info:\n");
    //        PRFORMd("DeltaT_tNQC",DeltaT_nqc);
    //        PRFORMd("tNQC[bare]",tNQC);
    //    }

    /* Find jmax: t[jmax] <= tNQC */
    int jmax = 0;
    for (int j=0; j<size; j++) {
        if(t[j] > tNQC) {
            jmax = j-2;
            break;
        }
    }

    /* Solve the linear systems */

    /* Regge-Wheeler-Zerilli normalized amplitude.
     The ringdown coefficient refer to this normalization.
     Nagar & Rezzolla, CQG 22 (2005) R167 */
    this = h;
    for (int k=0; k<KMAX; k++) {
        REAL8 nlm = 1./(sqrt( (this->l+2)*(this->l+1)*this->l*(this->l-1) ) );
        for (int j=0; j<size; j++) {
            p1tmp[k][j] = this->ampl->data->data[j] * nlm;
        }
        this = this->next;
    }

    /* Matrix elements: waveform amplitude at all points */
    for (int k=0; k<KMAX; k++) {
        for (int j=0; j<size; j++) {
            m11[k][j] = n1[j] * p1tmp[k][j];
            m12[k][j] = n2[j] * p1tmp[k][j];
        }
    }

    /* Take FD derivatives */
    for (int k=0; k<KMAX; k++) {
        D0(m11[k],dt,size, m21[k]);
        D0(m12[k],dt,size, m22[k]);
        D0(p1tmp[k],dt,size, p2tmp[k]);
    }

    //#if (TEOB_DEBUG)
    //    fp = fopen("nqc_amp_func.txt", "w");
    //    for (int j=0; j<size; j++) {
    //        fprintf(fp, "%e\t%e\t%e\n", t[j], p1tmp[1][j], p2tmp[1][j]);
    //    }
    //    fclose(fp);
    //#endif

    //    REAL8 detM = 1.;
    REAL8 oodetM = 1.;
    for (int k=0; k<KMAX; k++) {

        ai[k][0] = ai[k][1] = 0.;
        bi[k][0] = bi[k][1] = 0.;

        /* Computation of ai coefficients at Omega peak */
        P[0]     = max_A[k]  - p1tmp[k][jmax];
        P[1]     = max_dA[k] - p2tmp[k][jmax];

        M[0]     = m11[k][jmax];
        M[1]     = m12[k][jmax];
        M[2]     = m21[k][jmax];
        M[3]     = m22[k][jmax];

        /* detM     = M[0]*M[3]-M[1]*M[2];
         ai[k][0] = (M[3]*P[0] - M[1]*P[1])/detM;
         ai[k][1] = (M[0]*P[1] - M[2]*P[0])/detM; */
        /* safe version (amplitude can be zero) */
        oodetM   = 1.0/(M[0]*M[3]-M[1]*M[2]);
        if (isfinite((double) oodetM)) {
            ai[k][0] = (M[3]*P[0] - M[1]*P[1])*oodetM;
            ai[k][1] = (M[0]*P[1] - M[2]*P[0])*oodetM;
        }

        /* Computation of bi coefficients at Omega peak */
        P[0]     = omg[k][jmax]   - max_omg[k];
        P[1]     = domg[k][jmax]  - max_domg[k];

        M[0]     = d_n4[jmax];
        M[1]     = d_n5[jmax];
        M[2]     = d2_n4[jmax];
        M[3]     = d2_n5[jmax];

        /* detM     =  M[0]*M[3] - M[1]*M[2];
         bi[k][0] = (M[3]*P[0] - M[1]*P[1])/detM;
         bi[k][1] = (M[0]*P[1] - M[2]*P[0])/detM; */
        /* safe version (phase can be zero) */
        oodetM   = 1.0/(M[0]*M[3]-M[1]*M[2]);
        if (isfinite((double) oodetM)) {
            bi[k][0] = (M[3]*P[0] - M[1]*P[1])*oodetM;
            bi[k][1] = (M[0]*P[1] - M[2]*P[0])*oodetM;
        }

    }

    //    if (TEOB_VERBOSE){
    //        printf("NQC coefficients for 22 mode:\n");
    //        PRFORMd("a1",ai[1][0]);
    //        PRFORMd("a2",ai[1][1]);
    //        PRFORMd("b1",bi[1][0]);
    //        PRFORMd("b2",bi[1][1]);
    //    }

    /* Set amplitude and phase */
    // TODO: [optimization] vectorize!
    this = hnqc;
    for (int k=0; k<KMAX; k++) {
        for (int j=0; j<size; j++) {
            this->ampl->data->data[j]  = 1. + ai[k][0]*n1[j] + ai[k][1]*n2[j];
            this->phase->data->data[j] =      bi[k][0]*n4[j] + bi[k][1]*n5[j];
        }
        this = this->next;
    }

    this = h;
    SphHarmPolarTimeSeries *that = hnqc;
    /* Multiply waveform to NQC */
    // TODO: [optimization] vectorize!
    while (this && that) {
        for (int j=0; j<size; j++) {
            this->ampl->data->data[j] *= that->ampl->data->data[j];
            this->phase->data->data[j] -= that->phase->data->data[j];
        }
        this = this->next;
        that = that->next;
    }

    /* Free mem */
    for (int k=0; k<KMAX; k++) {
        XLALFree(omg[k]);
        XLALFree(domg[k]);
        XLALFree(m11[k]);
        XLALFree(m12[k]);
        XLALFree(m13[k]);
        XLALFree(m21[k]);
        XLALFree(m22[k]);
        XLALFree(p1tmp[k]);
        XLALFree(p2tmp[k]);
    }
    XLALFree(n1);
    XLALFree(n2);
    XLALFree(n3);
    XLALFree(n4);
    XLALFree(n5);
    XLALFree(n6);
    XLALFree(d_n4);
    XLALFree(d_n5);
    XLALFree(d2_n4);
    XLALFree(d2_n5);

}

/* Computes the factors and the coefficients that build the
 NQC corrections to the waveform in the spinning case.
 This routine works around merger with dyn and h and
 then add everything also to hlm */
void eob_wav_hlmNQC_find_a1a2a3_mrg(LALTEOBResumSDynamics *dyn_mrg, SphHarmPolarTimeSeries *hlm_mrg, SphHarmPolarTimeSeries *hnqc, LALTEOBResumSDynamics *dyn, SphHarmPolarTimeSeries *hlm)
{

    REAL8 A_tmp, dA_tmp, omg_tmp, domg_tmp;

    const REAL8 nu   = dyn->nu;
    const REAL8 X1   = dyn->X1;
    const REAL8 X2   = dyn->X2;
    const REAL8 chi1 = dyn->chi1;
    const REAL8 chi2 = dyn->chi2;
    const REAL8 aK   = dyn->a1 + dyn->a2;

    const REAL8 nu2  = SQ(nu);
    const REAL8 nu3  = nu2*nu;
    const REAL8 X12  = X1 - X2;
    const REAL8 aK2  = SQ(aK);
    const REAL8 aK3  = aK2*aK;
    const REAL8 aK4  = aK2*aK2;
    const REAL8 a12  = X1*chi1 - X2*chi2;
    const REAL8 aeff     = aK + 1./3.*a12*X12;
    const REAL8 aeff_omg = aK + a12*X12;

    REAL8 *t       = hlm_mrg->tdata->data;
    REAL8 *r       = dyn_mrg->data[TEOB_RAD];
    REAL8 *w       = dyn_mrg->data[TEOB_MOMG]; /* Omega */
    REAL8 UNUSED *pph     = dyn_mrg->data[TEOB_PPHI];
    REAL8 *pr_star = dyn_mrg->data[TEOB_PRSTAR];
    REAL8 *Omg_orb = dyn_mrg->data[TEOB_OMGORB]; /* Omega orbital */
    REAL8 *ddotr   = dyn_mrg->data[TEOB_DDOTR];

    REAL8 c_p1,     c_p2,     c_p3,   c_p4;
    REAL8 c_pdA1,   c_pdA2,   c_pdA3, c_pdA4;
    REAL8 c_pdomg1, c_pdomg2;
    REAL8 n0, d1;
    REAL8 a0_omg_tmp, a1_omg_tmp, a2_omg_tmp, b0_omg_tmp, b1_omg_tmp, b2_omg_tmp, a0_domg_tmp, a1_domg_tmp, a2_domg_tmp, b0_domg_tmp, b1_domg_tmp, b2_domg_tmp, a0_A_tmp, a1_A_tmp , a2_A_tmp, b0_A_tmp, b1_A_tmp, b2_A_tmp, a0_dA_tmp, a1_dA_tmp, a2_dA_tmp, b0_dA_tmp, b1_dA_tmp, b2_dA_tmp, omg_tmp_nu, omg_tmp_equal, domg_tmp_nu, domg_tmp_equal,  A_tmp_scale_nu, A_tmp_scale_equal, dA_tmp_scale_nu, dA_tmp_scale_equal ;

    REAL8 P[2], M[4], p1[2], p2[2], p3[2], p4[2], pA[5], pdA[5];
    REAL8 pomg[5], pdomg[5], pn0[2], pd1[2], ppdomg1[2], ppdomg2[2], pdA1[2],pdA2[2],pdA3[2],pdA4[2];
    REAL8 max_A[KMAX], max_dA[KMAX], UNUSED d2max[KMAX], UNUSED d3max[KMAX], max_omg[KMAX], max_domg[KMAX], UNUSED maxd2omg[KMAX], UNUSED DeltaT[KMAX];
    REAL8 ai[KMAX][2];
    REAL8 bi[KMAX][2];

    /* SphHarmSetTData has already taken care of pointing to time sequence */
    const int size = (int) hlm_mrg->tdata->length;

    REAL8 *omg[KMAX], *domg[KMAX];
    REAL8 *n1,*n2,*n3,*n4,*n5,*n6, *d_n4,*d_n5, *d2_n4,*d2_n5;
    REAL8 *m11[KMAX], *m12[KMAX], *m13[KMAX], *m21[KMAX], *m22[KMAX];
    REAL8 *p1tmp[KMAX], *p2tmp[KMAX]; /* RWZ amplitude and derivative */

    for (int k=0; k<KMAX; k++)
    {
        omg[k]  = (REAL8*) calloc (size,sizeof(REAL8));
        domg[k] = (REAL8*) calloc (size,sizeof(REAL8));
        m11[k]  = (REAL8*) calloc (size,sizeof(REAL8));
        m12[k]  = (REAL8*) calloc (size,sizeof(REAL8));
        m13[k]  = (REAL8*) calloc (size,sizeof(REAL8));
        m21[k]  = (REAL8*) calloc (size,sizeof(REAL8));
        m22[k]  = (REAL8*) calloc (size,sizeof(REAL8));
        p1tmp[k] = (REAL8*) calloc (size,sizeof(REAL8));
        p2tmp[k] = (REAL8*) calloc (size,sizeof(REAL8));
    }

    n1 = (REAL8*) calloc (size,sizeof(REAL8));
    n2 = (REAL8*) calloc (size,sizeof(REAL8));
    n3 = (REAL8*) calloc (size,sizeof(REAL8));
    n4 = (REAL8*) calloc (size,sizeof(REAL8));
    n5 = (REAL8*) calloc (size,sizeof(REAL8));
    n6 = (REAL8*) calloc (size,sizeof(REAL8));
    d_n4  = (REAL8*) calloc (size,sizeof(REAL8));
    d_n5  = (REAL8*) calloc (size,sizeof(REAL8));
    d2_n4 = (REAL8*) calloc (size,sizeof(REAL8));
    d2_n5 = (REAL8*) calloc (size,sizeof(REAL8));

    /* omega derivatives */
    const REAL8 dt = t[1]-t[0];
    SphHarmPolarTimeSeries *this;
    this = hlm_mrg;
    for (int k=0; k<KMAX; k++)
    {
        D0(this->phase->data->data, dt, size, omg[k]);
        D0(omg[k], dt, size, domg[k]);
        this = this->next;
    }

    /* NR fits */
    if (DEQUAL(nu,0.25,1e-9))
    {

        pA[0]    =  0.00178195;
        pA[1]    =  0.00435589;
        pA[2]    =  0.00344489;
        pA[3]    = -0.00076165;
        pA[4]    =  0.31973334;
        A_tmp    =  pA[0]*aK4    + pA[1]*aK3   + pA[2]*aK2    + pA[3]*aK     + pA[4];

        pdA[0]   =  0.00000927;
        pdA[1]   = -0.00024550;
        pdA[2]   =  0.00012469;
        pdA[3]   =  0.00123845;
        pdA[4]   = -0.00195014;
        dA_tmp   =  pdA[0]*aK4   + pdA[1]*aK3   + pdA[2]*aK2   + pdA[3]*aK   + pdA[4];

        pomg[0]  =  0.00603482;
        pomg[1]  =  0.01604555;
        pomg[2]  =  0.02290799;
        pomg[3]  =  0.07084587;
        pomg[4]  =  0.38321834;
        omg_tmp  =  pomg[0]*aK4  + pomg[1]*aK3  + pomg[2]*aK2  + pomg[3]*aK  + pomg[4];

        pdomg[0] =  0.00024066;
        pdomg[1] =  0.00038123;
        pdomg[2] = -0.00049714;
        pdomg[3] =  0.00041219;
        pdomg[4] =  0.01190548;
        domg_tmp =  pdomg[0]*aK4 + pdomg[1]*aK3 + pdomg[2]*aK2 + pdomg[3]*aK + pdomg[4];

    }  else if( nu > 0.16) {

        p1[0]      =  0.04680896;
        p1[1]      = -0.00632114;
        p2[0]      =  0.06586192;
        p2[1]      = -0.01180039;
        p3[0]      = -0.11617413;
        p3[1]      =  0.02704959;
        p4[0]      =  0.15597465;
        p4[1]      =  0.28034978;
        c_p1       =  p1[0]*nu + p1[1];
        c_p2       =  p2[0]*nu + p2[1];
        c_p3       =  p3[0]*nu + p3[1];
        c_p4       =  p4[0]*nu + p4[1];
        A_tmp      =  c_p1*aK3 + c_p2*aK2 + c_p3*aK + c_p4;

        pdA1[0]    = -0.00130824;
        pdA1[1]    =  0.00006202;
        pdA2[0]    =  0.00199855;
        pdA2[1]    = -0.00027474;
        pdA3[0]    =  0.00218838;
        pdA3[1]    =  0.00071540;
        pdA4[0]    = -0.00362779;
        pdA4[1]    = -0.00105397;
        c_pdA1     =  pdA1[0]*nu + pdA1[1];
        c_pdA2     =  pdA2[0]*nu + pdA2[1];
        c_pdA3     =  pdA3[0]*nu + pdA3[1];
        c_pdA4     =  pdA4[0]*nu + pdA4[1];
        dA_tmp     =  c_pdA1*aK3   + c_pdA2*aK2 + c_pdA3*aK+ c_pdA4;

        pn0[0]     =  0.46908067;
        pn0[1]     =  0.27022141;
        pd1[0]     =  0.64131115;
        pd1[1]     = -0.37878384;
        n0         =  pn0[0]*nu + pn0[1];
        d1         =  pd1[0]*nu + pd1[1];
        omg_tmp    =  n0/(1 + d1*aK);

        ppdomg1[0] =  0.00061175;
        ppdomg1[1] =  0.00074001;
        ppdomg2[0] =  0.02504442;
        ppdomg2[1] =  0.00548217;
        c_pdomg1   =  ppdomg1[0]*nu + ppdomg1[1];
        c_pdomg2   =  ppdomg2[0]*nu + ppdomg2[1];
        domg_tmp   =  c_pdomg1*aK   + c_pdomg2;

    }  else {

        /* Fit by G.Riemanschneider incorporating the test-particle NQC point
         obtained from the most-recent Teukolsky waveforms done by
         M. Colleoni using the 6PN-accurare iResum-radiation reaction.
         These points assure a smooth connection between merger and
         ringdown also outside the "calibration" domain, notably for
         large-mass ratios (though q<=20) and large (negative) spins
         Updated, 28/09/2017 */

        a0_omg_tmp    = -0.1460961247;
        a1_omg_tmp    =  0.0998056;
        a2_omg_tmp    = -0.118098;
        b0_omg_tmp    = -0.3430184009;
        b1_omg_tmp    =  0.0921551;
        b2_omg_tmp    = -0.0740285;
        omg_tmp_nu    = +0.5427169903*nu2 +0.2512395608*nu +0.2863992248;
        omg_tmp_equal =((a2_omg_tmp*X12*X12 + a1_omg_tmp*X12 + a0_omg_tmp)*aeff_omg+1)/((b2_omg_tmp*X12*X12 +b1_omg_tmp*X12 + b0_omg_tmp)*aeff_omg+1);
        omg_tmp       = omg_tmp_nu*omg_tmp_equal;

        a0_domg_tmp    = +0.0604556289;
        b0_domg_tmp    = -0.0299583285;
        a1_domg_tmp    = 0.0711715;
        a2_domg_tmp    = -0.0500886;
        b1_domg_tmp    = 0.0461239;
        b2_domg_tmp    = -0.0153068;

        domg_tmp_nu    = ( +0.0045213831*nu +0.0064934920)/( -1.4466409969*nu+1);
        domg_tmp_equal = (a2_domg_tmp*X12*X12 +a1_domg_tmp*X12 +b0_domg_tmp)*aeff_omg*aeff_omg +(b2_domg_tmp*X12*X12 +b1_domg_tmp*X12+a0_domg_tmp)*aeff_omg+1;
        domg_tmp       = domg_tmp_nu*domg_tmp_equal;

        a0_A_tmp     = -0.2750516062;
        b0_A_tmp     = -0.4693776065;
        a1_A_tmp     =  0.143066;
        a2_A_tmp     = -0.0425947;
        b1_A_tmp     =  0.176955;
        b2_A_tmp     = -0.111902;

        A_tmp_scale_nu    = -0.9862040409*nu3 +0.8167558040*nu2 -0.0427442282*nu+0.2948879452;
        A_tmp_scale_equal = ((a2_A_tmp*X12*X12 + a1_A_tmp*X12 +a0_A_tmp)*aeff+1)/((b2_A_tmp*X12*X12 + b1_A_tmp*X12 +b0_A_tmp)*aeff+1);
        A_tmp             = A_tmp_scale_nu*A_tmp_scale_equal*(1-0.5*omg_tmp*aeff);

        a0_dA_tmp     = +0.0037461628;
        b0_dA_tmp     = +0.0636082543;
        a1_dA_tmp     =  0.00129393;
        a2_dA_tmp     = -0.00239069;
        b1_dA_tmp     = -0.0534209;
        b2_dA_tmp     = -0.186101;

        dA_tmp_scale_nu    = ( -0.0847947167*nu -0.0042142765)/( +16.1559461812*nu+1);
        dA_tmp_scale_equal = ((a2_dA_tmp*X12*X12 + a1_dA_tmp*X12+ a0_dA_tmp)*aeff)/((b2_dA_tmp*X12*X12 + b1_dA_tmp*X12 + b0_dA_tmp)*aeff+1);
        dA_tmp             = (dA_tmp_scale_nu +dA_tmp_scale_equal)*omg_tmp;

    }

    /* Switch on the 22 values (only) */
    for (int k=0; k<KMAX; k++) {
        max_A[k]    = 0.;
        max_dA[k]   = 0.;
        max_omg[k]  = 0.;
        max_domg[k] = 0.;
    }
    max_A[1]    = A_tmp;
    max_dA[1]   = dA_tmp;
    max_omg[1]  = omg_tmp;
    max_domg[1] = domg_tmp;

    //    if (TEOB_VERBOSE) {
    //        printf("NR values for NQC determination:\n");
    //        PRFORMd("A22_mrg",max_A[1]);
    //        PRFORMd("dA22_mrg",max_dA[1]);
    //        PRFORMd("omg22_mrg",max_omg[1]);
    //        PRFORMd("domg22_mrg",max_domg[1]);
    //    }

    /* NQC corrections to AMPLITUDE (n1,n2,n3) and PHASE (n4,n5,n6)
     * NQC basis for (2,2) waveform : AMPLITUDE
     * note: n3 and n6 are not used
     */
    REAL8 pr_star2, r2, w2;
    for (int j=0; j<size; j++)
    {
        pr_star2 = SQ(pr_star[j]);
        r2       = SQ(r[j]);
        w2       = SQ(w[j]);               // FIXME: Omg or Omg_orbital ?
        n1[j]  = pr_star2/(r2*w2);         /* [pr*\/(r Omg)]^2 */
        n2[j]  = ddotr[j]/(r[j]*w2);       /* [ddot{r}/(r Omg^2)] */
        n4[j]  = pr_star[j]/(r[j]*w[j]);   /* pr*\/(r Omg) */
        n5[j]  = n4[j]*r2*w2;              /* (pr*)*(r Omg) */
    }

//#if (TEOB_DEBUG)
//    fp = fopen("nqc_nfunc.txt", "w");
//    for (int j=0; j<size; j++) {
//        fprintf(fp, "%20.12f\t%.16e\t%.16e\t%.16e\t%.16e\n", t[j], n1[j], n2[j], n4[j], n5[j]);
//    }
//    fclose(fp);
//#endif

    /* Derivatives for the phase */
    D0(n4,dt,size, d_n4);
    D0(n5,dt,size, d_n5);
    D0(d_n4,dt,size, d2_n4);
    D0(d_n5,dt,size, d2_n5);

//#if (TEOB_DEBUG)
//    fp = fopen("nqc_dfunc.txt", "w");
//    for (int j=0; j<size; j++) {
//        fprintf(fp, "%f\t%.16e\t%.16e\t%.16e\t%.16e\n", t[j], d_n4[j], d_n5[j], d2_n4[j], d2_n5[j]);
//    }
//    fclose(fp);
//#endif

    /* Find max Omg */
    INT4 Omgmax_index = 0;
    REAL8 Omg_max   = Omg_orb[0];
    for (int j=0; j<size; j++)
    {
        if (Omg_orb[j] > Omg_max)
        {
            Omg_max = Omg_orb[j];
            Omgmax_index = j;
        }
    }
    //TODO:
    //Test the search backwards
    /*
     int Omgmax_index = size-1;
     REAL8 Omg_max   = Omg_orb[Omgmax_index];
     for (int j=(size-2); j--; ) {
     if (Omg_orb[j] < Omg_max)
     break;
     Omg_max = Omg_orb[j];
     Omgmax_index = j;
     }
     */

    /* Time */
    REAL8 tOmgOrb_pk = t[Omgmax_index];
    REAL8 DeltaT_nqc = eob_nqc_timeshift(nu, chi1);
    REAL8 tNQC = tOmgOrb_pk - DeltaT_nqc;

//    if (TEOB_VERBOSE) {
//        printf("NQC info:\n");
//        PRFORMd("DeltaT_tNQC",DeltaT_nqc);
//        PRFORMd("tNQC[bare]",tNQC);
//    }

    /* Find jmax: t[jmax] <= tNQC */
    INT4 jmax = 0;
    for (INT4 j=0; j<size; j++)
    {
        if(t[j] > tNQC)
        {
            jmax = j-2;
            break;
        }
    }

    /* Solve the linear systems */

    /* Regge-Wheeler-Zerilli normalized amplitude.
     The ringdown coefficient refer to this normalization.
     Nagar & Rezzolla, CQG 22 (2005) R167 */
    this = hlm_mrg;
    for (int k=0; k<KMAX; k++)
    {
        REAL8 nlm = 1./(sqrt( (this->l+2)*(this->l+1)*this->l*(this->l-1) ) );
        for (int j=0; j<size; j++)
        {
            p1tmp[k][j] = this->ampl->data->data[j] * nlm;
        }
        this = this->next;
    }

    /* Matrix elements: waveform amplitude at all points */
    for (int k=0; k<KMAX; k++)
    {
        for (int j=0; j<size; j++)
        {
            m11[k][j] = n1[j] * p1tmp[k][j];
            m12[k][j] = n2[j] * p1tmp[k][j];
        }
    }

    /* Take FD derivatives */
    for (int k=0; k<KMAX; k++)
    {
        D0(m11[k],dt,size, m21[k]);
        D0(m12[k],dt,size, m22[k]);
        D0(p1tmp[k],dt,size, p2tmp[k]);
    }

//#if (TEOB_DEBUG)
//    fp = fopen("nqc_amp_func.txt", "w");
//    for (int j=0; j<size; j++) {
//        fprintf(fp, "%e\t%e\t%e\n", t[j], p1tmp[1][j], p2tmp[1][j]);
//    }
//    fclose(fp);
//#endif

//    REAL8 detM = 1.;
    REAL8 oodetM = 1.;
    for (int k=0; k<KMAX; k++) {

        ai[k][0] = ai[k][1] = 0.;
        bi[k][0] = bi[k][1] = 0.;

        /* Computation of ai coefficients at Omega peak */
        P[0]     = max_A[k]  - p1tmp[k][jmax];
        P[1]     = max_dA[k] - p2tmp[k][jmax];

        M[0]     = m11[k][jmax];
        M[1]     = m12[k][jmax];
        M[2]     = m21[k][jmax];
        M[3]     = m22[k][jmax];

        /* detM     = M[0]*M[3]-M[1]*M[2];
         ai[k][0] = (M[3]*P[0] - M[1]*P[1])/detM;
         ai[k][1] = (M[0]*P[1] - M[2]*P[0])/detM; */
        /* safe version (amplitude can be zero) */
        oodetM   = 1.0/(M[0]*M[3]-M[1]*M[2]);
        if (isfinite(oodetM)) {
            ai[k][0] = (M[3]*P[0] - M[1]*P[1])*oodetM;
            ai[k][1] = (M[0]*P[1] - M[2]*P[0])*oodetM;
        }

        /* Computation of bi coefficients at Omega peak */
        P[0]     = omg[k][jmax]   - max_omg[k];
        P[1]     = domg[k][jmax]  - max_domg[k];

        M[0]     = d_n4[jmax];
        M[1]     = d_n5[jmax];
        M[2]     = d2_n4[jmax];
        M[3]     = d2_n5[jmax];

        /* detM     =  M[0]*M[3] - M[1]*M[2];
         bi[k][0] = (M[3]*P[0] - M[1]*P[1])/detM;
         bi[k][1] = (M[0]*P[1] - M[2]*P[0])/detM; */
        /* safe version (phase can be zero) */
        oodetM   = 1.0/(M[0]*M[3]-M[1]*M[2]);
        if (isfinite(oodetM)) {
            bi[k][0] = (M[3]*P[0] - M[1]*P[1])*oodetM;
            bi[k][1] = (M[0]*P[1] - M[2]*P[0])*oodetM;
        }

    }

//    if (TEOB_VERBOSE){
//        printf("NQC coefficients for 22 mode:\n");
//        PRFORMd("a1",ai[1][0]);
//        PRFORMd("a2",ai[1][1]);
//        PRFORMd("b1",bi[1][0]);
//        PRFORMd("b2",bi[1][1]);
//    }

    /* Set amplitude and phase */
    this = hnqc;
    for (int k=0; k<KMAX; k++) {
        for (int j=0; j<size; j++) {
            this->ampl->data->data[j]  = 1. + ai[k][0]*n1[j] + ai[k][1]*n2[j];
            this->phase->data->data[j] =      bi[k][0]*n4[j] + bi[k][1]*n5[j];
        }
        this = this->next;
    }

    /* Multiply merger waveform to NQC */
    this = hlm_mrg;
    SphHarmPolarTimeSeries *that = hnqc;
    while (this && that) {
        for (int j=0; j<size; j++) {
            this->ampl->data->data[j] *= that->ampl->data->data[j];
            this->phase->data->data[j] -= that->phase->data->data[j];
            if (DEBUG)
                XLAL_CHECK_VOID((this->l == that->l) && (this->m == that->m), XLAL_EFAULT, "Spherical harmonic mode numbers do not agree.");
        }
        this = this->next;
        that = that->next;
    }
    if (DEBUG)
        XLAL_CHECK_VOID(!this && !that, XLAL_EFAULT, "Spherical harmonics link lists do not have the same length.");

    /* Multiply full waveform to NQC */

    r       = dyn->data[TEOB_RAD];
    w       = dyn->data[TEOB_MOMG]; /* Omega */
    pph     = dyn->data[TEOB_PPHI];
    pr_star = dyn->data[TEOB_PRSTAR];
    Omg_orb = dyn->data[TEOB_OMGORB]; /* Omega orbital */
    ddotr   = dyn->data[TEOB_DDOTR];

    free(n1);
    free(n2);
    free(n4);
    free(n5);

    const INT4 fullsize = hlm->tdata->length;

    n1 = (REAL8*) calloc (fullsize,sizeof(REAL8));
    n2 = (REAL8*) calloc (fullsize,sizeof(REAL8));
    n4 = (REAL8*) calloc (fullsize,sizeof(REAL8));
    n5 = (REAL8*) calloc (fullsize,sizeof(REAL8));

    for (int j=0; j<fullsize; j++) {
        pr_star2 = SQ(pr_star[j]);
        r2       = SQ(r[j]);
        w2       = SQ(w[j]); //CHECKME: Omg or Omg_orbital ?
        n1[j]  = pr_star2/(r2*w2);         /* [pr*\/(r Omg)]^2 */
        n2[j]  = ddotr[j]/(r[j]*w2);       /* [ddot{r}/(r Omg^2)] */
        n4[j]  = pr_star[j]/(r[j]*w[j]);   /* pr*\/(r Omg) */
        n5[j]  = n4[j]*r2*w2;              /* (pr*)*(r Omg) */
    }

    this = hlm;
    for (int k=0; k<KMAX; k++) {
        for (int j=0; j<fullsize; j++) {
            this->ampl->data->data[j] *= (1. + ai[k][0]*n1[j] + ai[k][1]*n2[j]);
            this->phase->data->data[j] -= (bi[k][0]*n4[j] + bi[k][1]*n5[j]);
        }
        this = this->next;
    }

    /* Free mem */
    for (int k=0; k<KMAX; k++)
    {
        free(omg[k]);
        free(domg[k]);
        free(m11[k]);
        free(m12[k]);
        free(m13[k]);
        free(m21[k]);
        free(m22[k]);
        free(p1tmp[k]);
        free(p2tmp[k]);
    }
    free(n1);
    free(n2);
    free(n3);
    free(n4);
    free(n5);
    free(n6);
    free(d_n4);
    free(d_n5);
    free(d2_n4);
    free(d2_n5);
}


/* Main routine for factorized EOB waveform */
void eob_wav_hlm(LALTEOBResumSDynamics *dyn, LALTEOBResumSWaveformModeSingleTime *hlm)
{

    const REAL8 nu = dyn->nu;
    const REAL8 chi1 = dyn->chi1;
    const REAL8 chi2 = dyn->chi2;
    const REAL8 a1 = dyn->a1;
    const REAL8 a2 = dyn->a2;
    const REAL8 X1 = dyn->X1;
    const REAL8 X2 = dyn->X2;
    const REAL8 C_Q1 = dyn->C_Q1;
    const REAL8 C_Q2 = dyn->C_Q2;
    const int usetidal = dyn->use_tidal;
    const int usespins = dyn->use_spins;
    const int usespeedytail = 1;
    const REAL8 X12 = X1-X2; /* sqrt(1-4nu) */

    const REAL8 t   = dyn->t;
    const REAL8 phi = dyn->phi;
    const REAL8 r   = dyn->r;
    const REAL8 UNUSED pph = dyn->pphi;
    const REAL8 prstar = dyn->prstar;
    const REAL8 Omega  = dyn->Omg;
    const REAL8 ddotr  = dyn->ddotr;
    const REAL8 H      = dyn->H;
    const REAL8 Heff   = dyn->Heff;
    const REAL8 jhat   = dyn->jhat;
    const REAL8 rw     = dyn->r_omega;

    hlm->time = t;

    /* Source term */
    REAL8 source[] = {
        jhat,Heff,
        Heff,jhat,Heff,
        jhat,Heff,jhat,Heff,
        Heff,jhat,Heff,jhat,Heff,
        jhat,Heff,jhat,Heff,jhat,Heff,
        Heff,jhat,Heff,jhat,Heff,jhat,Heff,
        jhat,Heff,jhat,Heff,jhat,Heff,jhat,Heff
    };

    /* Newtonian waveform */
    LALTEOBResumSWaveformModeSingleTime hNewt;
    eob_wav_hlmNewt(rw,Omega,phi,nu, &hNewt);

    if (usetidal) {
        /* Need to correct some of the m=odd modes.
         The Newtonian factor has a different normalization when entering the point-mass
         and the tidal term. The factor X12 = sqrt*1-4nu) is re-introduced in the point-mass term
         in eob_wav_hlm() */
        REAL8 vphi3 = gsl_pow_int(rw*Omega,3);
        hNewt.ampli[0] = ChlmNewt_ampli[0] * vphi3;
        hNewt.ampli[2] = ChlmNewt_ampli[2] * vphi3;
        hNewt.ampli[4] = ChlmNewt_ampli[4] * vphi3;
        REAL8 p4_vphi5 = (2.*nu-1) * gsl_pow_int(rw*Omega,5);
        hNewt.ampli[5]  = ChlmNewt_ampli[5]  * p4_vphi5;
        hNewt.ampli[7]  = ChlmNewt_ampli[7]  * p4_vphi5;
        hNewt.ampli[9]  = ChlmNewt_ampli[9]  * p4_vphi5;
        hNewt.ampli[11] = ChlmNewt_ampli[11] * p4_vphi5;
        hNewt.ampli[13] = ChlmNewt_ampli[13] * p4_vphi5;
    }

    if (usespins) {
        /* Special treatment when spin is on because of the singularity in the sqrt(1-4*nu)
         for m=odd mode and nu=1/4. See discussion in
         Damour & Nagar, PRD 90, 044018, Sec. 4, Eq.(89).
         This is not done for multipoles l>4 because no spinning information is included there */

        REAL8 vphi3 = gsl_pow_int(rw*Omega,3);
        hNewt.ampli[0] = ChlmNewt_ampli[0] * vphi3; /* (2,1) */
        hNewt.ampli[2] = ChlmNewt_ampli[2] * vphi3; /* (3,1) */
        hNewt.ampli[4] = ChlmNewt_ampli[4] * vphi3; /* (3,3) */

        REAL8 p4_vphi5 = (2.*nu-1) * gsl_pow_int(rw*Omega,5);
        hNewt.ampli[5]  = ChlmNewt_ampli[5]  * p4_vphi5; /* (4,1) */
        hNewt.ampli[7]  = ChlmNewt_ampli[7]  * p4_vphi5; /* (4,3) */
        hNewt.ampli[9]  = ChlmNewt_ampli[9]  * p4_vphi5 * X12; /* (5,1) */
        hNewt.ampli[11] = ChlmNewt_ampli[11] * p4_vphi5 * X12; /* (5,3) */
        hNewt.ampli[13] = ChlmNewt_ampli[13] * p4_vphi5 * X12; /* (5,5) */
    }

    /* Compute corrections */
    REAL8 rholm[KMAX], flm[KMAX];
    REAL8 x = SQ(rw*Omega);
    if (usespins){
        /* eob_wav_flm_s_old(x, nu, X1,X2,chi1,chi2,a1,a2,C_Q1,C_Q2,usetidal,rholm,flm); */
        dyn->eob_wav_flm_s(x, nu, X1, X2, chi1, chi2, a1, a2, C_Q1, C_Q2, dyn->clm, usetidal, rholm, flm);
    } else {
        /* eob_wav_flm_old(x, nu, rholm,flm); */
        eob_wav_flm(x, nu, dyn->clm, rholm, flm);
    }

    /* Computing the tail */
#define RTAIL (1.213061319425267e+00)
    const REAL8 Hreal = H * nu;
    LALTEOBResumSWaveformModeSingleTime tlm;
    if (usespeedytail) {
        eob_wav_speedyTail(Omega, Hreal, RTAIL, &tlm);
    } else {
        eob_wav_hhatlmTail(Omega, Hreal, RTAIL, &tlm);
    }

    /* Residual phase corrections delta_{lm} */
    REAL8 dlm[KMAX];
    eob_wav_deltalm(Hreal, Omega, nu, dlm);

    /* Point-mass h_lm */
    for (int k = 0; k < KMAX; k++) {
        hlm->ampli[k] =  hNewt.ampli[k] * flm[k] * source[k] * tlm.ampli[k];
        hlm->phase[k] = -( hNewt.phase[k] + tlm.phase[k] + dlm[k]); /* Minus sign by convention */
    }

    /* NQC */
    if (!(dyn->nqc_coeffs_hlm == NQC_OFF) &&
        !(dyn->nqc_coeffs_hlm == NQC_COMPUTE))
    {
        LALTEOBResumSWaveformModeSingleTime hNQC;
        /* eob_wav_hlmNQC_nospin2016(nu,r,prstar,Omega,ddotr, &hNQC); */
        eob_wav_hlmNQC(nu, r, prstar, Omega, ddotr, dyn->NQC->hlm, &hNQC);

        const INT4 maxk = MIN(KMAX, dyn->NQC->hlm->maxk+1);
        /* Add NQC correction where necessary */
        for (INT4 k = 0; k < maxk; k++) {
            if (dyn->NQC->hlm->activemode[k]) {
                hlm->ampli[k] *= hNQC.ampli[k];
                hlm->phase[k] -= hNQC.phase[k];
            }
        }
    }

    if (usetidal) {
        /* Tidal contribution */
        REAL8 hlmtidal[KMAX];
        eob_wav_hlmTidal(x, dyn, hlmtidal);
        if( !(usespins) ) {
            /* Correct normalization of point-mass wave for some of the m=odd modes */
            hlm->ampli[0] *= X12;
            hlm->ampli[2] *= X12;
            hlm->ampli[4] *= X12;
            hlm->ampli[5] *= X12;
            hlm->ampli[7] *= X12;
            hlm->ampli[9] *= X12;
            hlm->ampli[11] *= X12;
            hlm->ampli[13] *= X12;
        }
        /* Add tidal contribution to waveform */
        for (int k = 0; k < KMAX; k++) {
            hlm->ampli[k] += (hNewt.ampli[k] * tlm.ampli[k] * hlmtidal[k]);
        }

    }

    return;
}

/* Fit of c3, TEOBResumS paper Nagar et al. (2018)
 Note: c3 = 0 with tides*/
REAL8 eob_c3_fit_global(REAL8 nu, REAL8 UNUSED chi1, REAL8 UNUSED chi2, REAL8 UNUSED X1, REAL8 UNUSED X2, REAL8 a1, REAL8 a2)
{
    const REAL8 nu2 = nu*nu;
    const REAL8 nu3 = nu2*nu;
    const REAL8 X12 = sqrt(1.-4.*nu);
    const REAL8 a12 = a1+a2;

    /* Equal-mass, equal-spin coefficients */
    const REAL8 c0 =  43.371638;
    const REAL8 n1 =  -1.174839;
    const REAL8 n2 =   0.354064;
    const REAL8 d1 =  -0.151961;

    const REAL8 c3_eq = c0*(1. + n1*a12 + n2*a12*a12)/(1.+d1*a12);

    /* Coefficients 10/05/2018 */
    const REAL8 cnu    =  929.579;
    const REAL8 cnu2   = -9178.87;
    const REAL8 cnu3   =  23632.3;
    const REAL8 ca1_a2 = -104.891;

    const REAL8 c3_uneq = cnu*a12*nu*X12 + cnu2*a12*nu2*X12 + cnu3*a12*nu3*X12 + ca1_a2*(a1-a2)*nu2;

    return c3_eq + c3_uneq;
}

/* Time-shift for NQC */
REAL8 eob_nqc_timeshift(REAL8 nu, REAL8 chi1)
{

    REAL8 DeltaT_nqc = 1.0;

    /* Additional time-shift only needed ONLY for large, negative, spins.
     This change from 1.0 to 4.0 eliminates unphysical features in the
     frequency related to the imperfect behavior of the NQC functions */

    /* Old Delta_T NQC
     if ((chi1 <-0.85) && (nu <= 14./225.)) {
     DeltaT_nqc = 4.0;
     } else {
     DeltaT_nqc = 1.0; // standard choice inspired by test-particle results
     }
     */

    /* New Delta_T NQC determined in TEOBResumS paper (arXiv:1806.01772) */
    if (((chi1 < -0.9) && (nu < 8./81.)) || ((chi1 < -0.8) && (nu < 11./144.))) {
        DeltaT_nqc = 4.0;
    } else {
        DeltaT_nqc = 1.0; // standard choice inspired by test-particle results
    }

    return DeltaT_nqc;
}

/* Set NQC coefficients */
void eob_nqc_setcoefs(LALTEOBResumSDynamics *dyn)
{
    NQCdata *nqc;
    nqc = dyn->NQC;
    REAL8 nu = dyn->nu;

    nqc->flx->add = 1;
    nqc->hlm->add = 1;

    if (dyn->nqc_coeffs_hlm == NQC_OFF)
        nqc->hlm->add = 0;
    if (dyn->nqc_coeffs_flx == NQC_OFF)
        nqc->flx->add = 0;

    /* Init NQC coefs to zero */
    for (int k = 0; k < KMAX; k++) {
        for (int j = 0; j < 6; j++) {
            nqc->flx->n[k][j] = 0.;
            nqc->hlm->n[k][j] = 0.;
        }
        nqc->flx->a1[k] = 0.;
        nqc->flx->a2[k] = 0.;
        nqc->flx->a3[k] = 0.;
        nqc->flx->b1[k] = 0.;
        nqc->flx->b2[k] = 0.;
        nqc->flx->b3[k] = 0.;
        nqc->flx->activemode[k] = 0;
        nqc->hlm->a1[k] = 0.;
        nqc->hlm->a2[k] = 0.;
        nqc->hlm->a3[k] = 0.;
        nqc->hlm->b1[k] = 0.;
        nqc->hlm->b2[k] = 0.;
        nqc->hlm->b3[k] = 0.;
        nqc->hlm->activemode[k] = 0;
    }
    nqc->flx->maxk = -1;
    nqc->hlm->maxk = -1;

    if (nqc->flx->add + nqc->hlm->add == 0)
        return;

    /*
     * NOTE: Option coefs from file is not used;
     * NR implies nrfit_nospin201602
     */

    // if (STREQUAL(par_get_s("nqc_coefs_flx"),"nrfit_nospin201602"))
    // TODO: ADD HERE YOUR LATEST FITS
    //else if (STREQUAL(par_get_s("nqc_coefs_hlm"),"nrfit_spin_202001"))
    if (dyn->nqc_coeffs_flx == NQC_NR_NOSPIN)
        eob_nqc_setcoefs_nospin201602(nu, nqc->flx);


    //    if (STREQUAL(par_get_s("nqc_coefs_hlm"),"nrfit_nospin201602"))
    // TODO: ADD HERE YOUR LATEST FITS
    //else if (STREQUAL(par_get_s("nqc_coefs_hlm"),"nrfit_spin_202001"))
    if (dyn->nqc_coeffs_hlm == NQC_NR_NOSPIN)
        eob_nqc_setcoefs_nospin201602(nu, nqc->hlm);

}

/* Set NQC coefficients
 NR fits for nonspinning case 2016/02/09
 Hardcoded in eob_wav_hlmNQC_nospin201602() */
void eob_nqc_setcoefs_nospin201602(REAL8 nu, NQCcoefs *nqc)
{

    const REAL8 xnu  = 1-4*nu;
    const REAL8 xnu2 = SQ(xnu);

    const INT4 k21 = 0;
    const INT4 k22 = 1;
    const INT4 k33 = 4;

    nqc->activemode[k21]=1;
    nqc->activemode[k22]=1;
    nqc->activemode[k33]=1;

    /* (2,1) */
    nqc->a1[k21] =  0.0162387198*(7.32653082*xnu2 + 1.19616248*xnu + 0.73496656);
    nqc->a2[k21] = -1.80492460*xnu2 + 1.78172686*xnu + 0.30865284;
    nqc->a3[k21] =  0.0;

    nqc->b1[k21] = -0.0647955017*(3.59934444*xnu2 - 4.08628784*xnu + 1.37890907);
    nqc->b2[k21] =  1.3410693180*(0.38491989*xnu2 + 0.10969453*xnu + 0.97513971);
    nqc->b3[k21] =  0.0;

    /* (2,2) */
    nqc->a1[k22] = -0.0805236959*( 1 - 2.00332326*xnu2)/( 1 + 3.08595088*xnu2);
    nqc->a2[k22] =  1.5299534255*( 1 + 1.16438929*xnu2)/( 1 + 1.92033923*xnu2);
    nqc->a3[k22] =  0.0;

    nqc->b1[k22] = 0.146768094955*( 0.07417121*xnu + 1.01691256);
    nqc->b2[k22] = 0.896911234248*(-0.61072011*xnu + 0.94295129);
    nqc->b3[k22] = 0.0;

    /* (3,3) */
    nqc->a1[k33] = -0.0377680000*(1 - 14.61548907*xnu2)/( 1 + 2.44559263*xnu2);
    nqc->a2[k33] =  1.9898000000*(1 + 2.09750346 *xnu2)/( 1 + 2.57489466*xnu2);
    nqc->a3[k33] =  0.0;

    nqc->b1[k33] = 0.1418400000*(1.07430512 - 1.23906804*xnu + 4.44910652*xnu2);
    nqc->b2[k33] = 0.6191300000*(0.80672432 + 4.07432829*xnu - 7.47270977*xnu2);
    nqc->b3[k33] = 0.0;

    nqc->add = 1;
    nqc->maxk = k33;

}


/* QNM fits for the 22 mode for spinning systems */
void QNMHybridFitCab(REAL8 nu,
                     REAL8 X1,
                     REAL8 X2,
                     REAL8 chi1,
                     REAL8 chi2,
                     REAL8 aK,
                     REAL8 Mbh,
                     REAL8 abh,
                     REAL8 *a1,
                     REAL8 *a2,
                     REAL8 *a3,
                     REAL8 *a4,
                     REAL8 *b1,
                     REAL8 *b2,
                     REAL8 *b3,
                     REAL8 *b4,
                     REAL8 *sigmar,
                     REAL8 *sigmai,
                     int usespins)
{

    const REAL8 a12        = X1*chi1 - X2*chi2;
    const REAL8 X12        = X1 - X2;
    const REAL8 aeff       = aK + 1./3.*a12*X12;
    const REAL8 aeff_omg   = aK + a12*X12;
    const REAL8 af         = abh;
    const REAL8 nu2        = SQ(nu);
    const REAL8 nu3        = nu2*nu;
    const REAL8 aeff2      = SQ(aeff);
    const REAL8 aeff3      = aeff2*aeff;
    const REAL8 af2        = SQ(af);
    const REAL8 af3        = af2*af;
    const REAL8 aeff_omg2  = SQ(aeff_omg);
    const REAL8 aeff_omg3  = aeff_omg2*aeff_omg;
    const REAL8 aeff_omg4  = SQ(aeff_omg2);
    const REAL8 X12_2      = SQ(X12);

    REAL8 alpha21[KMAX], alpha1[KMAX], omega1[KMAX], c3A[KMAX], c3phi[KMAX], c4phi[KMAX], Domg[KMAX], Amrg[KMAX], c2A[KMAX];


    INT4 modeon[KMAX];
    const INT4 k21 = 0;
    const INT4 k22 = 1;
    const INT4 k33 = 4;
    const INT4 UNUSED k44 = 8;

    for (int k=0; k<KMAX; k++) {
        modeon[k] = 0; /* off */
        sigmar[k] = sigmai[k] = 0.;
        a1[k] = a2[k] = a3[k] = a4[k] = 0.;
        b1[k] = b2[k] = b3[k] = b4[k] = 0.;
    }

    if (!(usespins)) {

        modeon[k21]=1;
        modeon[k22]=1;
        modeon[k33]=1;

        /* Last updates: 05/09/2017 from CoM extrapolated SXS data */

        // l=2 -------------------------------------------------------------------

        /* (l=2, m=2)*/
        alpha21[k22] = -0.3025985041156393 *nu2 +  0.0032794155172817 *nu +  0.1828276903682022;
        alpha1[k22]  = -0.1615300454109702 *nu2 +  0.0147030662812516 *nu +  0.0878204175700328;
        c3A[k22]     =  0.8118901739129283 *nu  -  0.5584875090785957;
        c3phi[k22]   =  0.7156419884962878 *nu  +  3.8436474282409803;
        c4phi[k22]   =  2.2336960710670901 *nu  +  1.4736119175780844;
        Domg[k22]    =  0.8846304360111242 *nu2 +  0.0872792137250448 *nu +  0.1058414813686749;
        Amrg[k22]     = 1.4935750287318139 *nu2 +  0.2157497669089671 *nu +  1.4292027468283439;

        /* (l=2, m=1)*/
        alpha21[k21] = -0.2741607253846813 *nu2 +  0.0079342900879431 *nu +  0.1835522430667348;
        alpha1[k21]  = -0.1277546304610336 *nu2 +  0.0093615534859368 *nu +  0.0882855170502398;
        c3A[k21]     = -0.9431151070942140 *nu  +  0.2569989171628133;
        c3phi[k21]   = -3.4479482376671666 *nu  +  2.4755856452648359;
        c4phi[k21]   = -3.4024504071619841 *nu  +  1.0650118588151427;
        Domg[k21]    =  0.2660644668923829 *nu2 +  0.2276854484140649 *nu +  0.0884880283627388;
        Amrg[k21]    = -5.7236432632743952 *nu2 +  0.0390010969627653 *nu +  0.4291847351869338;

        // l=3 ------------------------------------------------------------------
        /* (l=3,m=3)*/
        alpha21[k33] = -0.3620553934265325 *nu2 +  0.0171973908686402 *nu +  0.1865364041200878;
        alpha1[k33]  = -0.1821867653548689 *nu2 +  0.0134440240947561 *nu +  0.0916720214797975;
        c3A[k33]     =  2.7565431398030675 *nu  -  0.5506682334306747;
        c3phi[k33]   = -0.2497526471104979 *nu  +  2.3737675006958683;
        c4phi[k33]   = -2.9538823110315420 *nu  +  1.4483501341373066;
        Domg[k33]    =  1.3341439550896721 *nu2 -  0.1717105341058959 *nu +  0.1694617455660599;
        Amrg[k33]    = -9.3034388918614841 *nu2 +  1.0189351143222705 *nu +  0.4533252110436300;

        // l=4 ------------------------------------------------------------------
        /* (l=4,m=4)*/
        /* alpha21[k44] = -0.3991680748908423 *nu2 +   0.0287698202159666 *nu +  0.1880112530796091; */
        /* alpha1[k44]  = -0.2003781755488581 *nu2 +   0.0171888841352427 *nu +  0.0930836242032652; */
        /* c3A[k44]     =  3.1899853343683140 *nu  +  -0.4131730594856833; */
        /* c3phi[k44]   = 31.5753575286023747 *nu  +  -1.0375600524681363; */
        /* c4phi[k44]   = 25.4170586178559716 *nu  +  -0.4151371540505313; */
        /* Domg[k44]    = -1.5342842283421341 *nu2 +   1.5224173843877831 *nu +  0.0897013049238634; */
        /* Amrg[k44]    =  0.9438333992719329 *nu2 +  -1.0464153920266663 *nu +  0.2897769169572948; */

        //sigma[k][0] = -0.208936*nu3 - 0.028103*nu2 - 0.005383*nu + 0.08896;
        //sigma[k][1] =  0.733477*nu3 + 0.188359*nu2 + 0.220659*nu + 0.37367;
        sigmar[k21] = -0.208936*nu3 - 0.028103*nu2 - 0.005383*nu + 0.08896;
        sigmai[k21] =  0.733477*nu3 + 0.188359*nu2 + 0.220659*nu + 0.37367;

        sigmar[k22] = -0.364177*nu3 + 0.010951*nu2 - 0.010591*nu + 0.08896;
        sigmai[k22] =  2.392808*nu3 + 0.051309*nu2 + 0.449425*nu + 0.37365;

        sigmar[k33] = -0.319703*nu3 - 0.030076*nu2-0.009034*nu + 0.09270;
        sigmai[k33] =  2.957425*nu3 + 0.178146*nu2 + 0.709560*nu + 0.59944;

    } else {

        modeon[k22]=1;

        /* Setting up coefficients from the phenomenological description of the ringdown.
         For notation: Damour&Nagar, PRD 90 (2015), 024054 and Del Pozzo & Nagar, PRD 95 (2017), 124034
         Current global fits are new. See Nagar+ 2017 (in preparation) for a global performance
         and Riemenschneider& Nagar (2017) in preparation for the description of the fits */

        /* omg1 - imaginary part of the fundamental mode */
        REAL8 omega1_c    = -0.0598837831 * af3 + 0.8082136788 * af2 - 1.7408467418 * af + 1;
        REAL8 omega1_d    = -0.2358960279 * af3 + 1.3152369374 * af2 - 2.0764065380 * af + 1;
        omega1[k22]        =  0.3736716844 * (omega1_c/omega1_d);

        /* alpha1 - real part (damping time) of the fundamental mode */
        REAL8 alpha1_c    =  0.1211263886 * af3 + 0.7015835813 * af2 - 1.8226060896 * af + 1;
        REAL8 alpha1_d    =  0.0811633377 * af3 + 0.7201166020 * af2 - 1.8002031358 * af + 1;
        alpha1[k22]        =  0.0889623157 * (alpha1_c/alpha1_d);

        /* alpha2 - alpha1 */
        REAL8 alpha21_c   =  0.4764196512 * af3 - 0.0593165805 * af2 - 1.4168096833 * af + 1;
        REAL8 alpha21_d   =  0.4385578151 * af3 - 0.0763529088 * af2 - 1.3595491146 * af + 1;
        alpha21[k22]       =  0.1849525596 * (alpha21_c/alpha21_d);

        /* c3A */
        REAL8 a_c3A     =  0.0169543;
        REAL8 b_c3A     = -0.0799343;
        REAL8 c_c3A     = -0.115928;
        REAL8 c3A_nu       =  0.8298678603 * nu - 0.5615838975;
        REAL8 c3A_eq       =  (c_c3A * X12 + 0.0907476903) * aeff3 + (b_c3A * X12 + 0.0227344099) * aeff2 + (a_c3A * X12 - 0.1994944332)*aeff;
        c3A[k22]            =  c3A_nu + c3A_eq;

        /* c3_phi */
        REAL8 a_c3phi      = -0.462321;
        REAL8 b_c3phi      = -0.904512;
        REAL8 c_c3phi      =  0.437747;
        REAL8 d_c3phi      =  1.8275;
        REAL8 c3phi_nu     =  0.4558467286 * nu + 3.8883812141;
        REAL8 c3phi_equal  =  (d_c3phi*X12-2.0575868122) * aeff_omg4 +(c_c3phi*X12-0.5051534498)*aeff_omg3 +(b_c3phi*X12+2.5742292762)*aeff_omg2 +(a_c3phi*X12+2.5599640181)*aeff_omg;
        c3phi[k22]          = c3phi_nu + c3phi_equal;

        /* c4_phi */
        REAL8 a_c4phi      = -0.449976;
        REAL8 b_c4phi      = -0.980913;
        REAL8 c4phi_nu     =  2.0822327682 * nu + 1.4996868401;
        REAL8 c4phi_equal  =  (b_c4phi*X12+3.5695199109) * aeff_omg2 + (a_c4phi * X12 + 4.1312404030) * aeff_omg;
        c4phi[k22]          =  c4phi_nu + c4phi_equal;

        /* omg_mrg: the "merger frequency", i.e. the frequency at the peak of |h22| */
        /* Special scaling and independent variables used for the fit. AN&GR 2017 */
        REAL8 a2_omgmx     = -0.122735;
        REAL8 a1_omgmx     =  0.0857478;
        REAL8 b2_omgmx     = -0.0760023;
        REAL8 b1_omgmx     =  0.0826514;
        REAL8 omgmx_eq_c   =  (a2_omgmx*X12_2 +a1_omgmx*X12 -0.1416002395) * aeff_omg + 1;
        REAL8 omgmx_eq_d   =  (b2_omgmx*X12_2 +b1_omgmx*X12 -0.3484804901) * aeff_omg + 1;
        REAL8 omgmx_eq     =  omgmx_eq_c/omgmx_eq_d;
        REAL8 omgmx        =  (0.481958619443355 * nu2 + 0.223976694441952 * nu + 0.273813064427363) * omgmx_eq;

        /* the peak of the h22 metric (strain) waveform.*/
        /* Special scaling and independent variables used for the fit. AN& GR 2017*/
        REAL8 a2_A_scaled = -0.0820894;
        REAL8 a1_A_scaled = 0.176126;
        REAL8 b2_A_scaled = -0.150239;
        REAL8 b1_A_scaled = 0.20491;
        REAL8 A_scaled_eq = ((a2_A_scaled*X12*X12 + a1_A_scaled*X12 -0.2935238329)*aeff + 1)/((b2_A_scaled*X12*X12 + b1_A_scaled*X12 -0.4728707630)*aeff + 1);
        REAL8 A_scaled    = (+1.826573640739664*nu2 +0.100709438291872*nu +1.438424467327531)*A_scaled_eq;

        Amrg[k22]      = A_scaled*(1-0.5*omgmx*aeff);
        Domg[k22]      = omega1[k22] - Mbh*omgmx;

        /* renaming real & imaginary part of the QNM complex frequency sigma */
        //sigma[k22][0] = alpha1[k22];
        //sigma[k22][1] = omega1[k22];
        sigmar[k22] = alpha1[k22];
        sigmai[k22] = omega1[k22];

    }

    for (int k=0; k<KMAX; k++) {
        if (modeon[k]) {
            c2A[k] = 0.5*alpha21[k];
            REAL8 cosh_c3A = cosh(c3A[k]);
            a1[k] = Amrg[k] * alpha1[k] * cosh_c3A * cosh_c3A / c2A[k];
            a2[k] = c2A[k];
            a3[k] = c3A[k];
            a4[k] = Amrg[k] - a1[k] * tanh(c3A[k]);
            b2[k] = alpha21[k];
            b3[k] = c3phi[k];
            b4[k] = c4phi[k];
            b1[k] = Domg[k] * (1+c3phi[k]+c4phi[k]) / (b2[k]*(c3phi[k] + 2.*c4phi[k]));
        }
    }

    return;
}


/* Ringdown waveform template */
void eob_wav_ringdown_template(REAL8 x,
                               REAL8 a1,
                               REAL8 a2,
                               REAL8 a3,
                               REAL8 a4,
                               REAL8 b1,
                               REAL8 b2,
                               REAL8 b3,
                               REAL8 b4,
                               REAL8 sigmar,
                               REAL8 sigmai,
                               REAL8 *psi)
{
    REAL8 amp   = ( a1 * tanh(a2*x +a3) + a4 ) ;
    REAL8 phase = -b1*log((1. + b3*exp(-b2*x) + b4*exp(-2.*b2*x))/(1.+b3+b4));
    psi[0] = amp * exp(-sigmar*x); /* amplitude */
    psi[1] = - (phase - sigmai*x); /* phase, minus sign in front by convention */
}


/* Ringdown calculation and match to the dynamics */
void eob_wav_ringdown(LALTEOBResumSDynamics *dyn, SphHarmPolarTimeSeries *hlm)
{
    const REAL8 Mbh   = dyn->Mbhf;
    const REAL8 abh   = dyn->abhf;
    const REAL8 nu    = dyn->nu;
    const REAL8 chi1  = dyn->chi1;
    const REAL8 chi2  = dyn->chi2;
    const REAL8 X1    = dyn->X1;
    const REAL8 X2    = dyn->X2;
    const REAL8 aK    = dyn->a1+dyn->a2;

    const REAL8 xnu   = (1.-4.*nu);
    const REAL8 ooMbh = 1./Mbh;
    const REAL8 dt = dyn->dt;

    SphHarmPolarTimeSeries *this_hlm = hlm;

    //REAL8 *Omega = dyn->data[TEOB_MOMG];
    REAL8 *Omega = dyn->data[TEOB_OMGORB]; /* use this for spin */

    /* Note:
     dynsize < size , since the wf has been extended
     but the two time arrays agree up to dynsize */
    const UINT4 dynsize = dyn->size;
    const INT4 size = hlm->tdata->length;
    REAL8 *t = hlm->tdata->data;

    const INT4 k21 = 0;
    const INT4 k22 = 1;
    const INT4 k33 = 4;

    /* Find peak of Omega */
    /* Assume a monotonically increasing function, then the peak
     start from after the peak */
    UINT4 index_pk = dynsize-1;
    REAL8 Omega_pk = Omega[index_pk];
    for (INT4 j = dynsize-2; j-- ; ) {
        if (Omega[j] < Omega_pk)
            break;
        index_pk = j;
        Omega_pk = Omega[j];
    }

#if (1)

  /* This is a hard-fix that always guarantee the 7 points */
  /* Make sure to comment the following line in main:
     dt_merger_interp = MIN(dt_merger_interp, (dyn->time[size-1] - dyn->tMOmgpeak)/4 );
     and uncomment:
     dt_merger_interp = MIN(dt_merger_interp, dyn->dt);
  */
  REAL8 *Omega_ptr = &Omega[index_pk-3];
  REAL8 tOmg_pk; /* New interpolated value of the Omega peak */
  const INT4 n = 7; /* USE 7, it seems we need at least 7 points to determine t_Omega_peak properly */
  REAL8 tmax = dyn->time[index_pk];

  if ( (index_pk + (n-1)/2) > (dynsize-1) ) {
    /* Here there are not enough points after the Omega peak
       We always need 3; we compute what we need by linear extrapolation */
    REAL8 Omega_pk_grid[7]; /* Temporary buffer for the 7-point interp */
    const INT4 ni = (index_pk + (n-1)/2) - (dynsize-1) ; /* Pts to extrap, 0 <  ni <= 3 */

    /* Copy the pts we have */
    for (int j = 0; j < (7-ni); j++)
      Omega_pk_grid[j] = Omega_ptr[j];
    /* Extrapolate the others */
    if (ni==1) {
      Omega_pk_grid[6] = 2.*Omega_pk_grid[5]-Omega_pk_grid[4];
      //Omega_pk_grid[6] =3.*Omega_pk_grid[5]-3.*Omega_pk_grid[4]+Omega_pk_grid[3];//quadratic, PLEASE CHECK
    } else if (ni==2) {
      Omega_pk_grid[5] = 2.*Omega_pk_grid[4]-Omega_pk_grid[3];
      Omega_pk_grid[6] = 2.*Omega_pk_grid[5]-Omega_pk_grid[4];
    } else if (ni==3) {
      Omega_pk_grid[4] = 2.*Omega_pk_grid[3]-Omega_pk_grid[2];
      Omega_pk_grid[5] = 2.*Omega_pk_grid[4]-Omega_pk_grid[3];
      Omega_pk_grid[6] = 2.*Omega_pk_grid[5]-Omega_pk_grid[4];
    } else XLAL_ERROR_VOID(XLAL_EBADLEN, "Wrong counting (ni).\n");

    /* Now we have 7 */
    tOmg_pk = find_max(n, dt, tmax, Omega_pk_grid, NULL);
  } else {
    /* Here everything is good */
    tOmg_pk = find_max(n, dt, tmax, Omega_ptr, NULL);
  }

  /* Scale peak value by BH mass */
  tOmg_pk *= ooMbh;

#else

    const INT4 n = 7; /* USE 7, it seems we need at least 7 points to determine t_Omega_peak properly */
    XLAL_CHECK_VOID( ( (index_pk + (n-1)/2) <= (dynsize-1) ), XLAL_EBADLEN, "Not enough points to interpolate.\n");

    REAL8 tmax = dyn->time[index_pk];
    REAL8 *Omega_ptr = &Omega[index_pk-3];
    REAL8 tOmg_pk = find_max(n, dt, tmax, Omega_ptr, NULL);
    tOmg_pk *= ooMbh;

#endif

    /* Merger time t_max(A22) */
    REAL8 DeltaT_nqc = eob_nqc_timeshift(nu, chi1);
    REAL8 tmrg[KMAX], tmatch[KMAX], dtmrg[KMAX];

    /* nonspinning case */
    /* tmrg[k22]  = tOmg_pk-3./Mbh; */ /* OLD */
    REAL8 tmrgA22 = tOmg_pk-(DeltaT_nqc + 2.)/Mbh;

    for (int k=0; k<KMAX; k++) {
        tmrg[k] = tmrgA22;
    }

    /* The following values are the difference between the time of the peak of
     the 22 waveform and the 21 and 33. These specific values refer to the
     nonspinning case. They are different in the spinning case, which
     is however not implemented. These are here only as placeholder */
    dtmrg[k21] = 5.70364338 + 1.85804796*xnu  + 4.0332262*xnu*xnu; //k21
    dtmrg[k33] = 4.29550934 - 0.85938*xnu;                         //k33
    tmrg[k21]  = tmrg[k22] + dtmrg[k21]/Mbh;     // t_max(A21) => peak of 21 mode
    tmrg[k33]  = tmrg[k22] + dtmrg[k33]/Mbh;     // t_max(A33) => peak of 33 mode

    /* Postmerger-Ringdown matching time */
    for (int k=0; k<KMAX; k++) {
        tmatch[k] = 2.*ooMbh + tmrg[k];
    }

    /* Compute QNM */
    REAL8 sigma[2][KMAX];
    REAL8 a1[KMAX], a2[KMAX], a3[KMAX], a4[KMAX];
    REAL8 b1[KMAX], b2[KMAX], b3[KMAX], b4[KMAX];
    QNMHybridFitCab(nu, X1, X2, chi1, chi2, aK,  Mbh, abh,
                    a1, a2, a3, a4, b1, b2, b3, b4,
                    sigma[0], sigma[1], dyn->use_spins);

    /* Define a time vector for each multipole, scale by mass
     Ringdown of each multipole has its own starting time */
    REAL8 *t_lm[KMAX];
    for (int k=0; k<KMAX; k++) {
        t_lm[k] =  malloc ( size * sizeof(REAL8) );
        for (int j = 0; j < size; j++ ) {
            t_lm[k][j] = t[j] * ooMbh;
        }
    }

    /* Find attachment index */
    int idx[KMAX];
    for (int k = 0; k < KMAX; k++) {
        for (int j = size-1; j-- ; ) {
            if (t_lm[k][j] < tmatch[k]) {
                idx[k] = j - 1;
                break;
            }
        }
    }

    /* Compute Ringdown waveform for t>=tmatch */
    REAL8 t0, tm, psi[2];
    REAL8 Deltaphi[KMAX];

    for (int k = 0; k < KMAX; k++) {
        XLAL_CHECK_VOID(this_hlm, XLAL_EBADLEN, "Mode does not exist. Reached NULL pointer instead.\n");
        /* Calculate Deltaphi */
        t0 = t_lm[k][idx[k]] - tmrg[k];
        eob_wav_ringdown_template(t0, a1[k], a2[k], a3[k], a4[k], b1[k], b2[k], b3[k], b4[k], sigma[0][k], sigma[1][k], psi);
        Deltaphi[k] = psi[1] - this_hlm->phase->data->data[idx[k]];
        /* Compute and attach ringdown */
        for (int j = idx[k]; j < size ; j++ ) {
            tm = t_lm[k][j] - tmrg[k];
            eob_wav_ringdown_template(tm, a1[k], a2[k], a3[k], a4[k], b1[k], b2[k], b3[k], b4[k], sigma[0][k], sigma[1][k], psi);
            this_hlm->ampl->data->data[j] = psi[0];
            this_hlm->phase->data->data[j] = psi[1] - Deltaphi[k];
        }
        this_hlm = this_hlm->next;
    }
    XLAL_CHECK_VOID( !(this_hlm), XLAL_EBADLEN, "More modes\n");

    /* Free mem. */
    for (int k=0; k<KMAX; k++)
    {
        if(t_lm[k]) free(t_lm[k]);
    }

}


/* Root function to compute light-ring */
//TODO: THIS IS FOR NOSPIN
REAL8 eob_dyn_fLR(REAL8 r, void  *params)
{
    LALTEOBResumSDynamics *dyn = params;
    REAL8 A,B,dA,d2A,dB;
    //if (dyn->use_spins) eob_metric_s(r, dyn, &A,&B,&dA,&d2A,&dB);
    //else
    eob_metric  (r, dyn, &A,&B,&dA,&d2A,&dB);
    REAL8 u = 1./r;
    REAL8 dA_u = (-dA)*SQ(r);
    return A + 0.5 * u * dA_u;
}

/* Root finder for adiabatic light-ring */
int eob_dyn_adiabLR(LALTEOBResumSDynamics *dyn, REAL8 *rLR, INT4 tidesFlag)
{
    int status;
    INT4 iter = 0, max_iter = MAX_ITER;
    const double epsabs = 0.; /* if converges, precision is |r-r*| = epsabs + epsrel r*  */
    const double epsrel = 1e-10;
    const gsl_root_fsolver_type *T;
    double x, x_lo, x_hi;

    /* Set interval to search root */
    if (dyn->use_tidal) {
        /* Tides are always temporarily set as = NNLO to compute LR,
         But we may want to define different searches intervals */
        switch (tidesFlag) {
            case LAL_SIM_INSPIRAL_GETIDES_NNLO:
            case LAL_SIM_INSPIRAL_GETIDES_GSF2:
            case LAL_SIM_INSPIRAL_GETIDES_GSF23:
                /* BNS */
                x_lo = 2.1; // nu~1/4 kappaT2 ~ 12
                x_hi = 5.9; // nu~1/4 kappaT2 ~ 600
                break;
            default:
                XLAL_ERROR_REAL8(XLAL_EINVAL, "Invalid tides flag.\n");
        }
        if (dyn->bhns > 0) {
            /* BHNS */
            x_lo = 1.8;
            x_hi = 5.6; // nu~1/4 kappaT2 ~ 600
        }
    }
    else {
        /* BBH */
        x_lo = 1.8; // 1.818461553848201e+00 nu = 1/4
        x_hi = 3.1; // 3. nu = 0
        /* x_lo = 0.9*eob_approxLR(dyn->nu);
         x_hi = 1.1*eob_approxLR(dyn->nu); */
    }

    gsl_root_fsolver *s;
    gsl_function F;
    F.function = &eob_dyn_fLR;
    F.params = dyn;
    //T = gsl_root_fsolver_bisection;
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, x_lo, x_hi);

    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        x      = gsl_root_fsolver_root (s);
        x_lo   = gsl_root_fsolver_x_lower (s);
        x_hi   = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (x_lo, x_hi, epsabs, epsrel);
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free (s);

    *rLR = 0.;
    if (isfinite(x)) *rLR = (REAL8) x;

    //if (status == ???) {
    //  return ROOT_ERRORS_BRACKET;
    //}
    if (status == GSL_SUCCESS) {
        return XLAL_SUCCESS;
    }
    if (iter >= max_iter) {
        return XLAL_EMAXITER;
    }
    if (status != GSL_SUCCESS) {
        return XLAL_EFUNC;
    }

    return status;
}


/* Root function to compute LSO */
//TODO: THIS IS FOR NOSPIN
//REAL8 eob_dyn_fLSO(REAL8 r, void  *params)
//{
//    LALTEOBResumSDynamics *dyn = params;
//    REAL8 A,B,dA,d2A,dB;
//    //if (dyn->use_spins) eob_metric_s(r, dyn, &A,&B,&dA,&d2A,&dB);
//    //else
//    eob_metric  (r, dyn, &A,&B,&dA,&d2A,&dB);
//    REAL8 u = 1./r;
//    REAL8 u2  = SQ(u);
//    REAL8 dA_u = (-dA)*SQ(r);
//    REAL8 d2A_u = d2A*SQ(r)*SQ(r) + 2*dA*SQ(r)*r;
//    dB = u2*dA_u + 2.*A*u;
//    REAL8 d2B = d2A_u*u2 + 4.*u*dA_u + 2*A;
//    return ( dA_u*d2B - d2A_u*(dB) );
//}

///* Root finder for adiabatic LSO */
//int eob_dyn_adiabLSO(LALTEOBResumSDynamics *dyn, REAL8 *rLSO)
//{
//    int status;
//    INT4 iter = 0, max_iter = MAX_ITER;
//    const double epsabs = 0.; /* if converges, precision is |r-r*| = epsabs + epsrel r*  */
//    const double epsrel = 1e-10;
//    const gsl_root_fsolver_type *T;
//    double x;
//    double x_lo = 4.5; // 4.532648e+00 nu= 1/4
//    double x_hi = 6.2; // 6 nu=0
//    if (dyn->use_tidal) x_hi = 36.;
//
//    gsl_root_fsolver *s;
//    gsl_function F;
//    F.function = &eob_dyn_fLSO;
//    F.params = dyn;
//    //T = gsl_root_fsolver_bisection;
//    T = gsl_root_fsolver_brent;
//    s = gsl_root_fsolver_alloc (T);
//    gsl_root_fsolver_set (s, &F, x_lo, x_hi);
//
//    do
//    {
//        iter++;
//        status = gsl_root_fsolver_iterate (s);
//        x      = gsl_root_fsolver_root (s);
//        x_lo   = gsl_root_fsolver_x_lower (s);
//        x_hi   = gsl_root_fsolver_x_upper (s);
//        status = gsl_root_test_interval (x_lo, x_hi, epsabs, epsrel);
//    }
//    while (status == GSL_CONTINUE && iter < max_iter);
//    gsl_root_fsolver_free (s);
//
//    *rLSO = 0.;
//    if (isfinite(x)) *rLSO = (REAL8) x;
//
//    //if (status == ???) {
//    //  return ROOT_ERRORS_BRACKET;
//    //}
//    if (status == GSL_SUCCESS) {
//        return ROOT_ERRORS_NO;
//    }
//    if (iter >= max_iter) {
//        return ROOT_ERRORS_MAXITS;
//    }
//    if (status != GSL_SUCCESS) {
//        return ROOT_ERRORS_NOSUCC;
//    }
//
//    return status;
//}


/* POST-ADIABATIC MODULE */

/* Post-adiabatic dynamics */
int eob_dyn_Npostadiabatic(LALTEOBResumSDynamics *dyn, const REAL8 r0)
{
    /* Unpack values */
    const REAL8 nu    = dyn->nu;
    const REAL8 S     = dyn->S;
    const REAL8 Sstar = dyn->Sstar;
    const REAL8 chi1  = dyn->chi1;
    const REAL8 chi2  = dyn->chi2;
    const REAL8 X1    = dyn->X1;
    const REAL8 X2    = dyn->X2;
    const REAL8 c3    = dyn->cN3LO;
    const REAL8 aK2   = dyn->aK2;
    const REAL8 a1    = dyn->a1;
    const REAL8 a2    = dyn->a2;
    const REAL8 C_Q1  = dyn->C_Q1;
    const REAL8 C_Q2  = dyn->C_Q2;
    const REAL8 C_Oct1 = dyn->C_Oct1;
    const REAL8 C_Oct2 = dyn->C_Oct2;
    const REAL8 C_Hex1 = dyn->C_Hex1;
    const REAL8 C_Hex2 = dyn->C_Hex2;
    const REAL8 z3    = 2.0*nu*(4.0-3.0*nu);
    const INT4 usetidal = dyn->use_tidal;
    const INT4 usespins = dyn->use_spins;
    // MA: probably
    const INT4 size = dyn->size;

    /* Parameters for post adiabatic dynamics */
    //    const REAL8 rmin = POSTADIABATIC_RMIN;
    //    const REAL8 dr = (r0 - rmin)/(size-1); /* Uniform grid spacing */
    const REAL8 dr = POSTADIABATIC_DR;

    /* Additional memory */
    REAL8 *buffer[TEOB_NV];
    for (int v=0; v < TEOB_NV; v++)
        buffer[v] = malloc(size * sizeof (REAL8));

    REAL8 *A_vec                  = buffer[0];
    REAL8 *dA_vec                 = buffer[1];
    REAL8 *B_vec                  = buffer[2];
    REAL8 *sqrtAbyB_vec           = buffer[3];
    REAL8 *rc_vec                 = buffer[4];
    REAL8 *drc_dr_vec             = buffer[5];
    REAL8 *uc2_vec                = buffer[6];
    REAL8 *duc_dr_vec             = buffer[7];
    REAL8 *dAuc2_dr_vec           = buffer[8];
    REAL8 *dG_dr_vec              = buffer[9];
    REAL8 *dG_dprstar_vec         = buffer[10];
    REAL8 *dG_dprstarbyprstar_vec = buffer[11];
    REAL8 *G0_vec                 = buffer[12];
    REAL8 *dG_dr0_vec             = buffer[13];
    REAL8 *E_vec                  = buffer[14];
    REAL8 *Heff_vec               = buffer[15];
    REAL8 *Heff_orb_vec           = buffer[16];
    REAL8 *dpphi_dr_vec           = buffer[17];
    REAL8 *dprstar_dr_vec         = buffer[18];
    REAL8 *dphi_dr_vec            = buffer[19];
    REAL8 *dt_dr_vec              = buffer[20];

    REAL8 ggm[14];
    REAL8 a_coeff, b_coeff, c_coeff, Delta, sol_p, sol_m, j02, uc, u2, prstar2, dHeff_dpphi, dHeff_dprstar, dHeff_dr, dHeff_dprstarbyprstar, d2Heff_dprstar20, H, G, pl_hold, x, jhat, psi, r_omg, v_phi, Fphi, dr_dtbyprstar, prstar4, Heff_orb_f, Heff_f, E_f;

    /*
     * Compute circular dynamics
     */

    for (int i = 0; i < size; i++)
    {
        /* Current radius */
        dyn->r = r0 - i*dr;

        /* Computing metric functions and centrifugal radius */
        if(usespins)
        {
            eob_metric_s(dyn->r,dyn, &A_vec[i], &B_vec[i], &dA_vec[i], &pl_hold, &pl_hold);
            dyn->eob_dyn_s_get_rc(dyn->r, nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1,
                                  C_Oct2, C_Hex1, C_Hex2, usetidal, &rc_vec[i],
                                  &drc_dr_vec[i], &pl_hold);
            eob_dyn_s_GS(dyn->r, rc_vec[i], drc_dr_vec[i], aK2, 0.0, 0.0, nu, chi1, chi2, X1, X2, c3, ggm);

            G                         = ggm[2] *S+ggm[3] *Sstar;    // tildeG = GS*S+GSs*Ss
            dG_dr_vec[i]              = ggm[6] *S+ggm[7] *Sstar;
            dG_dprstar_vec[i]         = ggm[4] *S+ggm[5] *Sstar;
            dG_dprstarbyprstar_vec[i] = ggm[10]*S+ggm[11]*Sstar;
        }
        else
        {
            eob_metric(dyn->r ,dyn, &A_vec[i], &B_vec[i], &dA_vec[i], &pl_hold, &pl_hold);

            rc_vec[i]                 = dyn->r; //Nonspinning case: rc = r
            drc_dr_vec[i]             = 1;

            G                         = 0.0;
            dG_dr_vec[i]              = 0.0;
            dG_dprstar_vec[i]         = 0.0;
            dG_dprstarbyprstar_vec[i] = 0.0;
        }

        /* Defining circular quantities for the flux calculation.
         Must not be overwritten in successive iterations, thus
         we define separate quantities with the subscripts 0. */
        G0_vec[i]        = G;
        dG_dr0_vec[i]    = dG_dr_vec[i];

        /* Auxiliary variables */
        sqrtAbyB_vec[i] = sqrt(A_vec[i]/B_vec[i]);
        uc              = 1./rc_vec[i];
        uc2_vec[i]      = uc*uc;
        duc_dr_vec[i]   = -uc2_vec[i]*drc_dr_vec[i];
        dAuc2_dr_vec[i] = uc2_vec[i]*(dA_vec[i]-2*A_vec[i]*uc*drc_dr_vec[i]);

        /* Computing the circular angular momentum by solving eq. (A15) of TEOBResumS paper
         (which is equivalent to solve eq.(4)=0 of arXiv:1805.03891).
         */

        if (usespins) {
            a_coeff = SQ(dAuc2_dr_vec[i]) - 4*A_vec[i]*uc2_vec[i]*SQ(dG_dr_vec[i]);  /* First coefficient of the quadratic equation a*x^2+b*x+c=0 */
            b_coeff = 2*dA_vec[i]*dAuc2_dr_vec[i] - 4*A_vec[i]*SQ(dG_dr_vec[i]);     /* Second coefficient of the quadratic equation */
            c_coeff = SQ(dA_vec[i]);                                                 /* Third coefficient of the quadratic equation */

            Delta   = SQ(b_coeff) - 4*a_coeff*c_coeff ;                              /* Delta of the quadratic equation */
            if (Delta<0)
            /* If the spins are very small,
               numerical fluctuations sometimes make Delta negative (e.g. -1e-30).
               Setting it to 0 by hand */
                Delta=0.;

            sol_p   = (-b_coeff + sqrt(Delta))/(2*a_coeff); /* Plus  solution of the quadratic equation */
            sol_m   = (-b_coeff - sqrt(Delta))/(2*a_coeff); /* Minus solution of the quadratic equation */

            /* dGdr sign determines choice of solution */
            if (dG_dr0_vec[i] > 0) j02 = sol_p;
            else       j02 = sol_m;

        }
        else
        {

            a_coeff = dAuc2_dr_vec[i];
            b_coeff = dA_vec[i];

            j02 = -b_coeff/a_coeff;

        }

        /* Define momenta in the circular orbit approximation */
        dyn->pphi                = sqrt(j02);
        dyn->prstar              = 0.0;
        dprstar_dr_vec[i]        = 0.0;

        /* Circular Hamiltonians, ref: arXiv: 1406.6913 */
        if(usespins) {

            eob_ham_s(nu,dyn->r,rc_vec[i],drc_dr_vec[i],dyn->pphi,dyn->prstar,S,
                      Sstar,chi1,chi2,X1,X2,aK2,c3,A_vec[i],dA_vec[i],
                      &H,               /* real EOB Hamiltonian divided by mu=m1m2/(m1+m2) */
                      &Heff_vec[i],     /* effective EOB Hamiltonian (divided by mu)       */
                      &Heff_orb_vec[i],
                      &dHeff_dr,        /* drvt Heff,r      */
                      NULL,             /* drvt Heff,prstar */
                      &dHeff_dpphi,     /* drvt Heff,pphi   */
                      &d2Heff_dprstar20);

            E_vec[i] = nu*H;

        } else {

            eob_ham(nu, dyn->r, dyn->pphi, dyn->prstar, A_vec[i], dA_vec[i],
                    &H,               /* real EOB Hamiltonian divided by mu=m1m2/(m1+m2) */
                    &Heff_orb_vec[i], /* effective EOB Hamiltonian (divided by mu). */
                    &dHeff_dr,        /* drvt Heff,r      */
                    NULL,             /* drvt Heff,prstar */
                    &dHeff_dpphi);    /* drvt Heff,pphi   */

            d2Heff_dprstar20 = 1/Heff_orb_vec[i];

            Heff_vec[i] = Heff_orb_vec[i]; /* Heff coincides with Heff_orb for the non-spinning case */
            E_vec[i] = nu*H;

        }

        /* Circular orbital frequency */
        dyn->Omg     = dHeff_dpphi/E_vec[i];

        /* Circular real orbital frequency */
        dyn->Omg_orb = (dyn->pphi*A_vec[i]*uc2_vec[i])/(E_vec[i]*Heff_orb_vec[i]);

        /* ddotr */
        dyn->ddotr   = -A_vec[i]/B_vec[i]*dHeff_dr*d2Heff_dprstar20;

        dyn->data[TEOB_RAD][i]    = dyn->r;
        dyn->data[TEOB_PPHI][i]   = dyn->pphi;
        dyn->data[TEOB_PRSTAR][i] = dyn->prstar;
        dyn->data[TEOB_DDOTR][i]  = dyn->ddotr;
        dyn->data[TEOB_MOMG][i]   = dyn->Omg;
        dyn->data[TEOB_OMGORB][i] = dyn->Omg_orb;

    } // END r-GRID FOR

    /* Computing angular momentum derivative */
    D0(dyn->data[TEOB_PPHI],-dr, size, dpphi_dr_vec); /* dJ0/dr */

    /*
     * Post-Adiabatic dynamics
     */

    int parity = 1; /* parity of the post-adiab iteration */

    /* For on PA orders */
    for (int n = 1; n <= POSTADIABATIC_N; n++) {

        /* Separating even and odd orders */
        if (n%2==0) parity = 0;
        else        parity = 1;

        /* For on r-grid */
        for (int i = 0; i < size; i++) {

            /* Setting loop variables to help reader */
            dyn->r       = dyn->data[TEOB_RAD][i];
            dyn->phi     = dyn->data[TEOB_PHI][i];
            dyn->pphi    = dyn->data[TEOB_PPHI][i];
            dyn->Omg     = dyn->data[TEOB_MOMG][i];
            dyn->ddotr   = dyn->data[TEOB_DDOTR][i];
            dyn->prstar  = dyn->data[TEOB_PRSTAR][i];
            dyn->Omg_orb = dyn->data[TEOB_OMGORB][i];

            if (parity)  {

                /* ***********************************
                 * Odd PA orders : prstar corrections
                 * ********************************** */

                /* Calculating the flux Fphi */
                //FIXME USE C-routines, jhat etc. are already present inside dynamics

                //FIXME Non-spinning routine gives 1e-2 difference between PA and full EOB waveform. Tested cases: bbh q 1 f 0.001 and q 5 f 0.006.
                if (usespins) {

                    /* Variables for which Kepler's law is still valid */
                    Heff_orb_f = sqrt(A_vec[i]*(1.0 + SQ(dyn->pphi)*uc2_vec[i]));
                    Heff_f     = G0_vec[i]*dyn->pphi + Heff_orb_f;
                    E_f        = sqrt(1 + 2*nu*(Heff_f - 1));
                    psi        = (duc_dr_vec[i] + dG_dr0_vec[i]*rc_vec[i]*sqrt(A_vec[i]/(SQ(dyn->pphi)) + A_vec[i]*uc2_vec[i])/A_vec[i])/(-0.5*dA_vec[i]);
                    r_omg      = 1.0/cbrt(SQ(((1./sqrt(rc_vec[i]*rc_vec[i]*rc_vec[i]*psi))+G0_vec[i])/(E_f)));
                    v_phi      = r_omg*dyn->Omg;
                    x          = SQ(v_phi);
                    jhat       = dyn->pphi/(r_omg*v_phi);

                    Fphi = eob_flx_Flux_s(x,dyn->Omg,r_omg, E_vec[i], Heff_vec[i],jhat,dyn->r,dyn->prstar, dyn->ddotr, dyn);

                } else {

                    Heff_orb_f = sqrt(A_vec[i]*(1.0 + SQ(dyn->pphi)*uc2_vec[i]));
                    Heff_f     = Heff_orb_f;
                    psi        = 2.*(1.0 + 2.0*nu*(Heff_orb_f - 1.0))/(SQ(dyn->r)*dA_vec[i]);
                    r_omg      = dyn->r*cbrt(psi);
                    v_phi      = r_omg*dyn->Omg;
                    x          = SQ(v_phi);
                    jhat       = dyn->pphi/(r_omg*v_phi);

                    Fphi = eob_flx_Flux(x,dyn->Omg,r_omg, E_vec[i], Heff_vec[i],jhat,dyn->r,dyn->prstar, dyn->ddotr, dyn);

                }

                /* Calculating prstar */
                dHeff_dprstarbyprstar = dyn->pphi*dG_dprstarbyprstar_vec[i] + (1+2*z3*A_vec[i]*uc2_vec[i]*SQ(dyn->prstar))/Heff_orb_vec[i];
                dr_dtbyprstar         = sqrtAbyB_vec[i]/(E_vec[i])*dHeff_dprstarbyprstar;
                dyn->prstar           = Fphi/dpphi_dr_vec[i]/dr_dtbyprstar;

                /* Note: p_phi does not change at odd orders
                 Computing first PA using the approximation detailed above A19 of TEOBResumS paper and Hamilton's equations.
                 */

                /* New GGM functions */
                eob_dyn_s_GS(dyn->r, rc_vec[i], drc_dr_vec[i], aK2, dyn->prstar, 0.0, nu, chi1, chi2, X1, X2, c3, ggm);

                dG_dr_vec[i]              = ggm[6] *S+ggm[7] *Sstar;
                dG_dprstar_vec[i]         = ggm[4] *S+ggm[5] *Sstar;
                dG_dprstarbyprstar_vec[i] = ggm[10]*S+ggm[11]*Sstar;

            } else {

                /* ***********************************
                 * Even PA orders : pphi corrections
                 * ********************************** */

                prstar4 = SQ(SQ(dyn->prstar));
                a_coeff = dAuc2_dr_vec[i];                   /* coefficients of the quadratic equation a*x^2+b*x+c=0 */
                b_coeff = 2*Heff_orb_vec[i]*(dG_dr_vec[i] + dG_dprstar_vec[i]*dprstar_dr_vec[i]);
                c_coeff = dA_vec[i] + 2*dyn->prstar*dprstar_dr_vec[i]*(1+2*z3*A_vec[i]*uc2_vec[i]*SQ(dyn->prstar)) + z3*dAuc2_dr_vec[i]*prstar4;
                Delta   = SQ(b_coeff) - 4*a_coeff*c_coeff;   /* Delta of the quadratic equation */

                /* sol_p = (-b_coeff + sqrt(Delta))/(2*a_coeff); */  /* Plus solution - Unphysical */
                sol_m = (-b_coeff - sqrt(Delta))/(2*a_coeff);  /* Minus solution of the quadratic equation */
                dyn->pphi = sol_m;

                /* Note: prstar and G functions do not change at even orders
                 G does not change because of the chosen gauge,
                 which eliminates the dependence of G from pphi).
                 */

            } //END IF-ELSE parity

            /* New Hamiltonians */
            if(usespins) {

                eob_ham_s(nu, dyn->r, rc_vec[i], drc_dr_vec[i], dyn->pphi, dyn->prstar, S, Sstar, chi1, chi2, X1, X2, aK2, c3, A_vec[i], dA_vec[i],
                          &H,               /* real EOB Hamiltonian divided by mu=m1m2/(m1+m2) */
                          &Heff_vec[i],     /* effective EOB Hamiltonian (divided by mu). Heff coincides with Heff_orb for the non-spinning case */
                          &Heff_orb_vec[i],
                          &dHeff_dr,        /* drvt Heff,r      */
                          &dHeff_dprstar,   /* drvt Heff,prstar */
                          &dHeff_dpphi,     /* drvt Heff,pphi   */
                          &d2Heff_dprstar20);

                E_vec[i] = nu*H;

            } else {

                eob_ham(nu, dyn->r, dyn->pphi, dyn->prstar, A_vec[i], dA_vec[i],
                        &H,               /* real EOB Hamiltonian divided by mu=m1m2/(m1+m2) */
                        &Heff_orb_vec[i], /* effective EOB Hamiltonian (divided by mu). Heff coincides with Heff_orb for the non-spinning case */
                        &dHeff_dr,        /* drvt Heff,r      */
                        &dHeff_dprstar,   /* drvt Heff,prstar */
                        &dHeff_dpphi);    /* drvt Heff,pphi   */

                u2      = 1./((dyn->r)*(dyn->r));
                prstar2 = (dyn->prstar)*(dyn->prstar);
                d2Heff_dprstar20 = (1. + 2.*A_vec[i]*u2*z3*prstar2)/Heff_orb_vec[i];

                Heff_vec[i] = Heff_orb_vec[i]; /* Heff coincides with Heff_orb for the non-spinning case */
                E_vec[i] = nu*H;

            }

            /* Orbital Frequency */
            dyn->Omg = dHeff_dpphi/E_vec[i];

            /* Real Orbital Frequency */
            dyn->Omg_orb = (dyn->pphi*A_vec[i]*uc2_vec[i])/(E_vec[i]*Heff_orb_vec[i]);

            /* ddotr */
            dyn->ddotr = -A_vec[i]/B_vec[i]*dHeff_dr*d2Heff_dprstar20;

            /* Time and phase radial derivatives */
            dt_dr_vec[i]   = E_vec[i]/(sqrtAbyB_vec[i]*dHeff_dprstar); /* dt_dr = 1/dr_dt */
            dphi_dr_vec[i] = dyn->Omg*dt_dr_vec[i];                    /* d(phi)_dr = d(phi)_dt*dt_dr */

            /* Re-assigning quantities to array elements */
            dyn->data[TEOB_PHI][i]    = dyn->phi;
            dyn->data[TEOB_PPHI][i]   = dyn->pphi;
            dyn->data[TEOB_MOMG][i]   = dyn->Omg;
            dyn->data[TEOB_DDOTR][i]  = dyn->ddotr;
            dyn->data[TEOB_PRSTAR][i] = dyn->prstar;
            dyn->data[TEOB_OMGORB][i] = dyn->Omg_orb;

        } // END R-GRID FOR

        /* Computing derivatives of the momenta */
        if (parity) D0(dyn->data[TEOB_PRSTAR],-dr, size, dprstar_dr_vec);
        else        D0(dyn->data[TEOB_PPHI],-dr, size, dpphi_dr_vec);

    } // END PA-CORRECTIONS FOR

    /*
     * Computing integrals for time and phase
     */

    /* Compute time */
    cumint3(dt_dr_vec, dyn->data[TEOB_RAD], size, dyn->time);

    /* Set last value for evolution */
    dyn->t = dyn->time[size-1];

    /* Compute orbital phase */
    cumint3(dphi_dr_vec, dyn->data[TEOB_RAD], size, dyn->data[TEOB_PHI]);


    /* Free memory */
    for (int v=0; v < TEOB_NV; v++)
        LALFree(buffer[v]);

    /*
     XLALFree(A_vec);
     XLALFree(dA_vec);
     XLALFree(B_vec);
     XLALFree(sqrtAbyB_vec);
     XLALFree(rc_vec);
     XLALFree(drc_dr_vec);
     XLALFree(uc2_vec);
     XLALFree(duc_dr_vec);
     XLALFree(dAuc2_dr_vec);
     XLALFree(dG_dr_vec);
     XLALFree(dG_dprstar_vec);
     XLALFree(dG_dprstarbyprstar_vec);
     XLALFree(G0_vec);
     XLALFree(dG_dr0_vec);
     XLALFree(E_vec);
     XLALFree(Heff_vec);
     XLALFree(Heff_orb_vec);
     XLALFree(dpphi_dr_vec);
     XLALFree(dprstar_dr_vec);
     XLALFree(dphi_dr_vec);
     XLALFree(dt_dr_vec);
     */

    return XLAL_SUCCESS;
}
