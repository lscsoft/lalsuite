#ifndef _LALSIM_IMR_PHENOMX_PRECESSION_H
#define _LALSIM_IMR_PHENOMX_PRECESSION_H
/*
 * Copyright (C) 2018/2019 Geraint Pratten
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


/**
 * \author Geraint Pratten
 *
 */

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include "LALSimIMRPhenomX.h"
#include "LALSimIMRPhenomX_internals.h"
#include "LALSimIMRPhenomXUtilities.h"

#include "LALSimIMRPhenomXHM.h"
#include "LALSimIMRPhenomXHM_internals.h"

#include <lal/LALStdlib.h>
#include <lal/LALSimIMR.h>
#include <lal/LALConstants.h>

#include <lal/FrequencySeries.h>
#include <lal/LALSimInspiral.h>


/**
 * Tolerance used below which numbers are treated as zero for the calculation of atan2
 */
#define MAX_TOL_ATAN 1.0e-15

/* Cartesian 3-vector */
typedef struct tagvector
{
    double x; // x-component
    double y; // y-component
    double z; // z-component
} vector;

/* Spherical Polar 3-vector */
typedef struct tagsphpolvector
{
    double r;     // radial-component
    double theta; // theta-component
    double phi;   // phi-component
} sphpolvector;

/* Struct to store all precession relevant variables */
typedef struct tagIMRPhenomXPrecessionStruct
{
        /* Flag to define the version of IMRPhenomXP called */
        INT4 IMRPhenomXPrecVersion; /**< Flag to set version of Euler angles used. */

        /* Debug flag */
        INT4 debug_prec; /**< Set debugging level. */

        /* Mass and spin weightings */
        REAL8 A1;        /**< Mass weighted pre-factor, see Eq. 3.2 of Schmidt et al, arXiv:1408.1810 */
        REAL8 A2;        /**< Mass weighted pre-factor, see Eq. 3.2 of Schmidt et al, arXiv:1408.1810 */
        REAL8 ASp1;      /**< \f$A1 * S_{1 \perp}\f$, see Eq. 3.3 of Schmidt et al, arXiv:1408.1810 */
        REAL8 ASp2;      /**< \f$A2 * S_{2 \perp}\f$, see Eq. 3.3 of Schmidt et al, arXiv:1408.1810 */
        REAL8 chi1_perp; /**< \f$ \chi_{1 \perp} \f$ */
        REAL8 chi2_perp; /**< \f$ \chi_{2 \perp} \f$ */
        REAL8 S1_perp;   /**< \f$ S_{1 \perp} \f$ */
        REAL8 S2_perp;   /**< \f$ S_{2 \perp} \f$ */
        REAL8 SL;        /**< \f$ \chi_{1 L} m^1_2 + \chi_{2 L} m^2_2 \f$ */
        REAL8 Sperp;     /**< \f$ \chi_{p} m^1_2 \f$ */
        REAL8 STot_perp; /**< \f$ S_{1 \perp} + S_{2 \perp} \f$ */
        REAL8 chi_perp;
        REAL8 chiTot_perp; /**< \f$ S_{1 \perp} + S_{2 \perp} \f$ */

        /* Effective precessing spin parameter: Schmidt et al, Phys. Rev. D 91, 024043 (2015), arXiv:1408.1810 */
        REAL8 chi_p;   /**< \f$ \chi_{p} = S_p / (A_1 \chi^2_1) \f$, Eq. 3.4 of Schmidt et al, arXiv:1408.1810 */

        /* Dimensionless aligned spin components on companion 1 and 2 respectively */
        REAL8 chi1L;   /**< \f$ \chi_{1L} = \chi_{1} \cdot \hat{L} \f$ */
        REAL8 chi2L;   /**< \f$ \chi_{2L} = \chi_{2} \cdot \hat{L} \f$ */

        /* Angle between J0 and line of sight (z-direction) */
        REAL8 thetaJN; /**< Angle between J0 and line of sight (z-direction) */

        /* The initial phase that we pass to the underlying aligned spin IMR waveform model */
        REAL8 phi0_aligned; /**< Initial phase to feed the underlying aligned-spin model */

        /* Angle to rotate the polarization */
        REAL8 zeta_polarization; /**< Angle to rotate the polarizations */

        /* Post-Newtonian Euler angles */
        REAL8 alpha0;       /**< Coefficient of \f$\alpha\f$ */
        REAL8 alpha1;       /**< Coefficient of \f$\alpha\f$ */
        REAL8 alpha2;       /**< Coefficient of \f$\alpha\f$ */
        REAL8 alpha3;       /**< Coefficient of \f$\alpha\f$ */
        REAL8 alpha4L;      /**< Coefficient of \f$\alpha\f$ */
        REAL8 alpha5;       /**< Coefficient of \f$\alpha\f$ */
        REAL8 alpha6;       /**< Coefficient of \f$\alpha\f$ */
        REAL8 alphaNNLO;    /**< Post Newtonian \f$\alpha\f$ at NNLO. */
        REAL8 alpha_offset; /**< Offset for \f$\alpha\f$ */

        REAL8 epsilon0;  /**< Coefficient of \f$\epsilon \f$ */
        REAL8 epsilon1;  /**< Coefficient of \f$\epsilon \f$ */
        REAL8 epsilon2;  /**< Coefficient of \f$\epsilon \f$ */
        REAL8 epsilon3;  /**< Coefficient of \f$\epsilon \f$ */
        REAL8 epsilon4L; /**< Coefficient of \f$\epsilon \f$ */
        REAL8 epsilon5;  /**< Coefficient of \f$\epsilon \f$ */
        REAL8 epsilon6;  /**< Coefficient of \f$\epsilon \f$ */
        REAL8 epsilonNNLO;    /**< Post Newtonian \f$\epsilon \f$ at NNLO. */
        REAL8 epsilon_offset; /**< Offset for \f$\epsilon \f$ */

        /* Alpha and epsilon offset for mprime !=2. alpha_offset corresponds to mprime=2 */
        REAL8 alpha_offset_1; /**< \f$\alpha\f$ offset passed to \f$m = 1\f$ modes. */
        REAL8 alpha_offset_3; /**< \f$\alpha\f$ offset passed to \f$m = 3\f$ modes. */
        REAL8 alpha_offset_4; /**< \f$\alpha\f$ offset passed to \f$m = 4\f$ modes. */
        REAL8 epsilon_offset_1; /**< \f$\epsilon\f$ offset passed to \f$m = 1\f$ modes. */
        REAL8 epsilon_offset_3; /**< \f$\epsilon\f$ offset passed to \f$m = 3\f$ modes. */
        REAL8 epsilon_offset_4; /**< \f$\epsilon\f$ offset passed to \f$m = 4\f$ modes. */

        /* Complex exponential of the Euler angles */
        COMPLEX16 cexp_i_alpha; /**< \f$e^{i \alpha}\f$ */
        COMPLEX16 cexp_i_epsilon; /**< \f$e^{i \epsilon}\f$ */
        COMPLEX16 cexp_i_betah;

        /* Multibanding applied to Euler angles */
        INT4 MBandPrecVersion; /**< Flag to control multibanding for Euler angles.  */

        /* Source Frame Variables */
        REAL8 J0x_Sf; /**< \f$ J_{0,x}\f$ in L frame. */
        REAL8 J0y_Sf; /**< \f$ J_{0,y}\f$ in L frame. */
        REAL8 J0z_Sf; /**< \f$ J_{0,z}\f$ in L frame. */
        REAL8 J0; /**< \f$ J_{0}\f$ in L frame. */
        REAL8 thetaJ_Sf; /**< Angle between \f$J_0\f$ and \f$ L_{\rm{N}} \f$ (z-direction) */
        REAL8 phiJ_Sf; /**< Azimuthal angle of \f$J_0\f$ in the L frame */
        REAL8 Nx_Sf; /**< Line-of-sight vector component \f$ N_{x}\f$ in L frame. */
        REAL8 Ny_Sf; /**< Line-of-sight vector component \f$ N_{y}\f$ in L frame. */
        REAL8 Nz_Sf; /**< Line-of-sight vector component \f$ N_{z}\f$ in L frame. */
        REAL8 Xx_Sf; /**< Component of triad basis vector X in L frame. */
        REAL8 Xy_Sf; /**< Component of triad basis vector X in L frame. */
        REAL8 Xz_Sf; /**< Component of triad basis vector X in L frame. */
        REAL8 kappa; /**< Eq. C12 of arXiv:XXXX.YYYY */

        /* J-frame variables */
        REAL8 Nx_Jf; /**< Line-of-sight vector component \f$ N_{x}\f$ in J frame. */
        REAL8 Ny_Jf; /**< Line-of-sight vector component \f$ N_{x}\f$ in J frame. */
        REAL8 Nz_Jf; /**< Line-of-sight vector component \f$ N_{x}\f$ in LJ frame. */
        REAL8 PArunx_Jf; /**< Component of triad basis vector P in J frame, arXiv:0810.5336. */
        REAL8 PAruny_Jf; /**< Component of triad basis vector P in J frame, arXiv:0810.5336. */
        REAL8 PArunz_Jf; /**< Component of triad basis vector P in J frame, arXiv:0810.5336. */
        REAL8 QArunx_Jf; /**< Component of triad basis vector Q in J frame, arXiv:0810.5336. */
        REAL8 QAruny_Jf; /**< Component of triad basis vector Q in J frame, arXiv:0810.5336. */
        REAL8 QArunz_Jf; /**< Component of triad basis vector Q in J frame, arXiv:0810.5336. */
        REAL8 XdotPArun; /**< \f$ X \cdot P \f$ */
        REAL8 XdotQArun; /**< \f$ X \cdot Q \f$ */

        /* Orbital angular momentum */
        REAL8 L0, L1, L2, L3, L4, L5, L6, L7, L8, L8L, LN, LOrb, LRef,LInit;

        /* Reference frequencies */
        REAL8 omega_ref, logomega_ref, omega_ref_cbrt, omega_ref_cbrt2;

        /* Spin weighted spherical harmonics */
        /* l = 2 */
        COMPLEX16 Y2m2, Y2m1, Y20, Y21, Y22;

        /* l = 3 */
        COMPLEX16 Y3m3, Y3m2, Y3m1, Y30, Y31, Y32, Y33;

        /* l = 4 */
        COMPLEX16 Y4m4, Y4m3, Y4m2, Y4m1, Y40, Y41, Y42, Y43, Y44;

        /* Useful sqare roots */
        REAL8 sqrt2, sqrt5, sqrt6, sqrt7, sqrt10, sqrt14, sqrt15, sqrt70, sqrt30, sqrt2p5;

        /* Variables for MSA precession angles of Chatziioannou et al, arXiv:1703.03967 */
        /* Lhat = {0,0,1} */
        REAL8 Lhat_theta, Lhat_phi, Lhat_norm, Lhat_cos_theta;

        /* Cartesian Dimensionful Spins */
        REAL8 S1x; /**< \f$ S_{1,x} \f$ in L frame */
        REAL8 S1y; /**< \f$ S_{1,y} \f$ in L frame */
        REAL8 S1z; /**< \f$ S_{1,z} \f$ in L frame */
        REAL8 S2x; /**< \f$ S_{2,x} \f$ in L frame */
        REAL8 S2y; /**< \f$ S_{2,y} \f$ in L frame */
        REAL8 S2z; /**< \f$ S_{2,z} \f$ in L frame */

        /* Spherical Polar Dimensionful Spins */
        REAL8 S1_norm; /**< \f$ \left| S_{1} \right| \f$ */
        REAL8 S1_theta; /**< Spherical polar component \f$ S_{1,\theta} \f$ in L frame */
        REAL8 S1_phi; /**< Spherical polar component \f$ S_{1,\phi} \f$ in L frame */
        REAL8 S1_cos_theta; /**< Spherical polar component \f$ \cos S_{1,\theta} \f$ in L frame */

        REAL8 S2_norm; /**< \f$ \left| S_{2} \right| \f$ */
        REAL8 S2_theta; /**< Spherical polar component \f$ S_{2,\theta} \f$ in L frame */
        REAL8 S2_phi; /**< Spherical polar component \f$ S_{2,\phi} \f$ in L frame */
        REAL8 S2_cos_theta; /**< Spherical polar component \f$ \cos S_{2,\theta} \f$ in L frame */

        /* Cartesian Dimensionless Spin Variables */
        REAL8 chi1x; /**< \f$ \chi_{1,x} \f$ in L frame */
        REAL8 chi1y; /**< \f$ \chi_{1,y} \f$ in L frame */
        REAL8 chi1z; /**< \f$ \chi_{1,z} \f$ in L frame */

        REAL8 chi2x; /**< \f$ \chi_{2,x} \f$ in L frame */
        REAL8 chi2y; /**< \f$ \chi_{2,y} \f$ in L frame */
        REAL8 chi2z; /**< \f$ \chi_{2,z} \f$ in L frame */

        /* Spherical Polar Dimensionless Spins */
        REAL8 chi1_theta; /**< Spherical polar component \f$ \chi_{1,\theta} \f$ in L frame */
        REAL8 chi1_phi; /**< Spherical polar component \f$ \chi_{1,\phi} \f$ in L frame */
        REAL8 chi1_norm; /**< \f$ \left| \chi_{1} \right| \f$ */
        REAL8 chi1_cos_theta; /**< Spherical polar component \f$ \cos \chi_{1,\theta} \f$ in L frame */

        REAL8 chi2_theta; /**< Spherical polar component \f$ \chi_{2,\theta} \f$ in L frame */
        REAL8 chi2_phi; /**< Spherical polar component \f$ \chi_{2,\phi} \f$ in L frame */
        REAL8 chi2_norm; /**< \f$ \left| \chi_{2} \right| \f$ */
        REAL8 chi2_cos_theta; /**< Spherical polar component \f$ \cos \chi_{2,\theta} \f$ in L frame */

        INT4 ExpansionOrder; /**< Flag to control expansion order of MSA system of equations. */

        REAL8 twopiGM, piGM;

        REAL8 L_Norm_N; /**< Norm of Newtonian orbital angular momentum \f$ \left| L_N \right| \f$ */
        REAL8 L_Norm_3PN; /**< Norm of orbital angular momentum at 3PN \f$ \left| L_{3 \rm{PN}} \right| \f$ */
        REAL8 J_Norm_N; /**< Norm of total angular momentum using Newtonian orbital angular momentum \f$ \left| J_{N} \right| \f$ */
        REAL8 J_Norm_3PN; /**< Norm of total angular momentum at 3PN \f$ \left| J_{3 \rm{PN}} \right| \f$ */

        REAL8 dotS1L; /**< \f$ S_1 \cdot \hat{L} \f$ */
        REAL8 dotS1Ln; /**< \f$ \hat{S}_1 \cdot \hat{L} \f$ */
        REAL8 dotS2L; /**< \f$ S_2 \cdot \hat{L} \f$ */
        REAL8 dotS2Ln; /**< \f$ \hat{S}_1 \cdot \hat{L} \f$ */
        REAL8 dotS1S2; /**< \f$ S_1 \cdot S_2 \f$ */
        REAL8 Seff; /**< \f$ S_{\rm{eff}} = (1 + q^{-1}) S_1 \cdot \hat{L} + (1 + q) S_2 \cdot \hat{L} \f$, Eq. 7 of arXiv:1703.03967. Note convention for q. */
        REAL8 Seff2; /**< \f$ S_{\rm{eff}}^2 \f$ */

        vector S1_0; /**< Initial value for \f$ S_{1} \f$ */
        vector S2_0; /**< Initial value for \f$ S_{2} \f$ */
        vector L_0; /**< Initial value for \f$ L \f$ */
        vector Lhat_0; /**< Initial value for \f$ \hat{L} \f$ */
        vector S_0; /**< Initial value for \f$ S_1 + S_2 \f$ */
        vector J_0; /**< Initial value for \f$ J \f$ */

        REAL8 S_0_norm, S_0_norm_2;
        REAL8 J_0_norm, J_0_norm_2;
        REAL8 L_0_norm, L_0_norm_2;

        REAL8 deltam_over_M; /**< \f$ (m_1 - m_2) / (m_1 + m_2) \f$ */

        //REAL8 phiz_0, phiz_1, phiz_2, phiz_3, phiz_4, phiz_5;
        REAL8 phiz_0_coeff, phiz_1_coeff, phiz_2_coeff, phiz_3_coeff, phiz_4_coeff, phiz_5_coeff;

        // Omegaz_i coefficients in D10 - D15
        REAL8 Omegaz0_coeff, Omegaz1_coeff, Omegaz2_coeff, Omegaz3_coeff, Omegaz4_coeff, Omegaz5_coeff;

        // Omegaz terms in D16 - D21
        REAL8 Omegaz0, Omegaz1, Omegaz2, Omegaz3, Omegaz4, Omegaz5;

        REAL8 Omegazeta0_coeff, Omegazeta1_coeff, Omegazeta2_coeff, Omegazeta3_coeff, Omegazeta4_coeff, Omegazeta5_coeff;
        REAL8 Omegazeta0, Omegazeta1, Omegazeta2, Omegazeta3, Omegazeta4, Omegazeta5;

        // MSA-SUA Euler Angles
        REAL8 phiz; /**< Azimuthal angle of L around J */
        REAL8 zeta; /**< Angle to describe L w.r.t. J */
        REAL8 cos_theta_L; /**< Cosine of polar angle between L and J */

        // First order MSA corrections
        REAL8 zeta_0_MSA; /**< First MSA correction \f$ \zeta_0 \f$, Eq. F19 of arXiv:1703.03967 */
        REAL8 phiz_0_MSA; /**< First MSA correction \f$ \phi_{z,0} \f$, Eq. 67 of arXiv:1703.03967 */

        // Initial values
        REAL8 zeta_0; /**< Initial value of \f$ \zeta \f$ */
        REAL8 phiz_0; /**< Initial value of \f$ \phi_{z,0} \f$ */

        // Orbital velocity, v and v^2
        REAL8 v; /**< Orbital velocity, \f$ v \f$ */
        REAL8 v2; /**< Orbital velocity squared, \f$ v^2 \f$ */

        // Reference orbital velocity, v and v^2
        REAL8 v_0; /**< Orbital velocity at reference frequency, \f$ v_{\rm{ref}} \f$ */
        REAL8 v_0_2; /**< Orbital velocity at reference frequency squared, \f$ v_{\rm{ref}}^2 \f$ */

        REAL8 Delta; /**< Eq. C3 of arXiv:1703.03967 */

        REAL8 D2;
        REAL8 D3;

        // Precession averaged total spin, Eq. 45
        REAL8 SAv; /**< \f$ S_{\rm{av}} \f$ as defined in Eq. 45 of arXiv:1703.03967 */
        REAL8 SAv2; /**< \f$ S_{\rm{av}}^2 \f$ */
        REAL8 invSAv2; /**< \f$ 1 / S_{\rm{av}}^2 \f$ */
        REAL8 invSAv; /**< \f$ 1 / S_{\rm{av}} \f$ */

        // Eq. C1, C2 for Eq. 51
        REAL8 psi1; /**< \f$ \psi_1 \f$ as defined in Eq. C1 of arXiv:1703.03967 */
        REAL8 psi2; /**< \f$ \psi_2 \f$ as defined in Eq. C2 of arXiv:1703.03967 */

        // Integration constant in Eq. 51
        REAL8 psi0; /**< \f$ \psi_0 \f$ as per Eq. 51 of arXiv:1703.03967 */

        // Eq. 51 and Eq. 24
        REAL8 psi; /**< \f$ \psi \f$ as defined by Eq. 51 of arXiv:1703.03967 */
        REAL8 psi_dot; /**< \f$ \dot{\psi} \f$ as per Eq. 50 of arXiv:1703.03967 */


        REAL8 Cphi; /**< \f$ C_{\phi} \f$ as defined by Eq. B14 of arXiv:1703.03967 */
        REAL8 Dphi; /**< \f$ C_{\phi} \f$ as defined by Eq. B15 of arXiv:1703.03967 */

        //REAL8 phiz_0_MSA_Cphi_term, phiz_0_MSA_Dphi_term;

        // PN Coefficients in Appendix A of Chatziioannou et al, PRD, 88, 063011, (2013)
        REAL8 a0, a1, a2, a3, a4, a5, a6, a7, b6, a0_2, a0_3, a2_2;
        REAL8 c0, c1, c2, c4, c12, c1_over_eta;

        // Eq. B9 - B11
        REAL8 d0, d2, d4;

        // Eq. A1 - A8
        REAL8 g0, g2, g3, g4, g5, g6, g7, g6L;

        // Spin-Orbit couplings
        REAL8 beta3, beta5, beta6, beta7, sigma4;

        // Total normed spin of primary and secondary BH squared
        REAL8 S1_norm_2, S2_norm_2;

        // Precession averaged spin couplings in A9 - A14
        REAL8 S1L_pav; /**< Precession averaged coupling \f$ \langle S_1 \cdot \hat{L} \rangle_{\rm{pr}} \f$, Eq. A9 of arXiv:1703.03967 */
        REAL8 S2L_pav; /**< Precession averaged coupling \f$ \langle S_2 \cdot \hat{L} \rangle_{\rm{pr}} \f$, Eq. A10 of arXiv:1703.03967 */
        REAL8 S1S2_pav; /**< Precession averaged coupling \f$ \langle S_1 \cdot S_2 \rangle_{\rm{pr}} \f$, Eq. A11 of arXiv:1703.03967 */
        REAL8 S1Lsq_pav; /**< Precession averaged coupling \f$ \langle (S_1 \cdot \hat{L})^2 \rangle_{\rm{pr}} \f$, Eq. A12 of arXiv:1703.03967 */
        REAL8 S2Lsq_pav; /**< Precession averaged coupling \f$ \langle (S_2 \cdot \hat{L})^2 \rangle_{\rm{pr}} \f$, Eq. A13 of arXiv:1703.03967 */
        REAL8 S1LS2L_pav; /**< Precession averaged coupling \f$ \langle (S_1 \cdot \hat{L}) (S_2 \cdot \hat{L}) \rangle_{\rm{pr}} \f$, Eq. A14 of arXiv:1703.03967 */

        // Total spin in Eq. 23 of Chatziioannou et al PRD, 95, 104004, (2017)
        REAL8 S_norm, S_norm_2;

        REAL8 Spl2; /**< Largest root of polynomial \f$ S^2_+ \f$, Eq. 22 of arXiv:1703.03967 */
        REAL8 Smi2; /**< Smallest root of polynomial \f$ S^2_- \f$, Eq. 22 of arXiv:1703.03967 */
        REAL8 S32; /**< Third root of polynomial \f$ S^2_3 \f$, Eq. 22 of arXiv:1703.03967 */
        REAL8 Spl; /**< \f$ S_+ \f$ */
        REAL8 Smi; /**< \f$ S_- \f$ */
        REAL8 S3; /**< \f$ S_3 \f$ */
        REAL8 Spl2mSmi2; /**< \f$ S^2_+ - S^2_- \f$ */
        REAL8 Spl2pSmi2; /**< \f$ S^2_+ + S^2_- \f$ */

        REAL8 A_coeff, B_coeff, C_coeff, D_coeff;

        REAL8 qq, invqq, eta, eta2, eta3, eta4;
        REAL8 delta_qq, delta2_qq, delta3_qq, delta4_qq;
        REAL8 inveta, inveta2, inveta3, inveta4, sqrt_inveta;

        REAL8 SConst;

        REAL8 LPN_coefficients[6];

        REAL8 constants_L[5];

        INT4 MSA_ERROR; /**< Flag to track errors in initialization of MSA system. */

} IMRPhenomXPrecessionStruct;

double IMRPhenomX_L_norm_3PN_of_v(const double v, const double v2, const double L_norm, IMRPhenomXPrecessionStruct *pPrec);

/* ~~~~~~~~~~ Perform transformations & initialize struct ~~~~~~~~~~ */
int IMRPhenomXGetAndSetPrecessionVariables(IMRPhenomXWaveformStruct *pWF, IMRPhenomXPrecessionStruct *pPrec,
	REAL8 m1_SI, REAL8 m2_SI,
	REAL8 S1x, REAL8 S1y, REAL8 S1z,
	REAL8 S2x, REAL8 S2y, REAL8 S2z,
	LALDict *lalParams, INT4 debug_flag
);

/* ~~~~~~~~~~ Post Newtonian Orbital Angular Momentum ~~~~~~~~~~ */
/* Functions to generate post-Newtonian Orbital Angular Momenta */
REAL8 XLALSimIMRPhenomXL2PNNS(REAL8 v, REAL8 eta);
REAL8 XLALSimIMRPhenomXL3PNAS(REAL8 v, REAL8 eta, REAL8 chi1L, REAL8 chi2L, REAL8 delta);
REAL8 XLALSimIMRPhenomXL4PNAS(REAL8 v, REAL8 eta, REAL8 chi1L, REAL8 chi2L, REAL8 delta);
REAL8 XLALSimIMRPhenomXL4PNLOSIAS(REAL8 v, REAL8 eta, REAL8 chi1L, REAL8 chi2L, REAL8 delta);

REAL8 XLALSimIMRPhenomXLPNAnsatz(REAL8 v, REAL8 LN, REAL8 L0, REAL8 L1, REAL8 L2, REAL8 L3, REAL8 L4, REAL8 L5, REAL8 L6, REAL8 L7, REAL8 L8, REAL8 L8L);

/* Function to check the maximum opening angle */
int IMRPhenomXPCheckMaxOpeningAngle(
	IMRPhenomXWaveformStruct *pWF,      /**< IMRPhenomX Waveform Struct */
  IMRPhenomXPrecessionStruct *pPrec   /**< IMRPhenomX Precession Struct */
);

/* ~~~~~~~~~~ Wigner Functions ~~~~~~~~~~ */
/* Wigner d Functions */
int IMRPhenomXWignerdCoefficients(
  REAL8 *cos_beta_half,   						/**< [out] cos(beta/2) */
  REAL8 *sin_beta_half,   						/**< [out] sin(beta/2) */
  const REAL8 v,          						/**< Cubic root of Pi * Frequency in geometric units */
	IMRPhenomXWaveformStruct *pWF,		  /**< IMRPhenomX Waveform Struct   */
	IMRPhenomXPrecessionStruct *pPrec 	/**< IMRPhenomX Precession Struct */
);

/* Wigner d Functions */
int IMRPhenomXWignerdCoefficients_cosbeta(
  REAL8 *cos_beta_half,   						/**< [out] cos(beta/2) */
  REAL8 *sin_beta_half,   						/**< [out] sin(beta/2) */
	REAL8 cos_beta                      /**< cos(beta) */
);


/* ~~~~~~~~~~ Routine to twist up aligned spin model ~~~~~~~~~~ */
/* Generic routine for twisting up 22-only models */
int IMRPhenomXPTwistUp22(
  const REAL8 Mf,                            /**< Frequency (Hz) */
  const COMPLEX16 hAS,                       /**< Underlying aligned-spin IMRPhenomXAS waveform*/
  IMRPhenomXWaveformStruct *pWF,             /**< IMRPhenomX Waveform Struct */
  IMRPhenomXPrecessionStruct *pPrec,         /**< IMRPhenomX Precession Struct */
  COMPLEX16 *hp,                             /**< [out] h_+ polarization \f$\tilde h_+\f$ */
  COMPLEX16 *hc                              /**< [out] h_x polarization \f$\tilde h_x\f$ */
);


/* ~~~~~~~~~~ NNLO post-Newtonian Euler Angles ~~~~~~~~~~ */
REAL8 XLALSimIMRPhenomXPNEuleralphaNNLO(
  REAL8 f,          /**< Geometric frequency                                                          */
  REAL8 eta,        /**< Symmetric mass rato                                                          */
  REAL8 chi1L,      /**< Dimensionless aligned spin of larger BH                                      */
  REAL8 chi2L,      /**< Dimensionless aligned spin of smaller BH                                     */
  REAL8 chip,       /**< Effective precession parameter: Schmidt, Ohme, Hannam, PRD, 91,024043 (2015) */
  REAL8 epsilon0    /**< Euler angle at reference Frequency, defines a constant offset                */
);

REAL8 XLALSimIMRPhenomXPNEulerepsilonNNLO(
  REAL8 f,          /**< Geometric frequency                                                          */
  REAL8 eta,        /**< Symmetric mass rato                                                          */
  REAL8 chi1L,      /**< Dimensionless aligned spin of larger BH                                      */
  REAL8 chi2L,      /**< Dimensionless aligned spin of smaller BH                                     */
  REAL8 chip,       /**< Effective precession parameter: Schmidt, Ohme, Hannam, PRD, 91,024043 (2015) */
  REAL8 epsilon0    /**< Euler angle at reference Frequency, defines a constant offset                */
);


/* Get alpha and epsilon offset depending of the mprime (second index of the non-precessing mode) */
void Get_alphaepsilon_atfref(
	REAL8 *alpha_offset,                 /**< [out] alpha offset for mprime */
	REAL8 *epsilon_offset, 							 /**< [out] epsilon offset for mprime */
	UINT4 mprime,                        /**< Second index non-precessing mode */
	IMRPhenomXPrecessionStruct *pPrec,   /**< IMRPhenomX Precession Struct */
	IMRPhenomXWaveformStruct *pWF        /**< IMRPhenomX Waveform Struct */
);

double IMRPhenomX_PN_Euler_alpha_NNLO(
  IMRPhenomXPrecessionStruct *pPrec,        /**< IMRPhenomX Precession Struct */
  const double omega,                       /**< Orbital frequency */
  const double omega_cbrt2,                 /**< Orbital frequency */
  const double omega_cbrt,                  /**< Cubic root of orbital frequency */
  const double logomega                     /**< Natural logarithm of orbital frequency */
);

double IMRPhenomX_PN_Euler_epsilon_NNLO(
  IMRPhenomXPrecessionStruct *pPrec,        /**< IMRPhenomX Precession Struct */
  const double omega,                       /**< Orbital frequency */
  const double omega_cbrt2,                 /**< Orbital frequency */
  const double omega_cbrt,                  /**< Cubic root of orbital frequency */
  const double logomega                     /**< Natural logarithm of orbital frequency */
);


/* ~~~~~~~~~~ MSA-SUA Euler Angle Functions ~~~~~~~~~~ */
double IMRPhenomX_Spin_Evolution_Equation_MSA(IMRPhenomXWaveformStruct *pWF, IMRPhenomXPrecessionStruct *pPrec);
vector IMRPhenomX_Return_Spin_Evolution_Coefficients_MSA(double LNorm, double JNorm, const IMRPhenomXPrecessionStruct *pPrec);
vector IMRPhenomX_Return_Constants_c_MSA(double v, const double JNorm, const IMRPhenomXPrecessionStruct *pPrec);
vector IMRPhenomX_Return_Constants_d_MSA(const double LNorm, const double JNorm, const IMRPhenomXPrecessionStruct *pPrec);
double IMRPhenomX_costhetaLJ(const double J_norm, const double L_norm, const double S_norm);
double IMRPhenomX_Return_Psi_MSA(double v, double v2, const IMRPhenomXPrecessionStruct *pPrec);
double IMRPhenomX_Return_Psi_dot_MSA(double v, const IMRPhenomXPrecessionStruct *pPrec);
double IMRPhenomX_Return_SNorm_MSA(const double v, IMRPhenomXPrecessionStruct *pPrec);
double IMRPhenomX_JNorm_MSA(const double LNorm, IMRPhenomXPrecessionStruct *pPrec);
int IMRPhenomX_Get_MSA_Euler_Angles(REAL8 *alpha, REAL8 *beta, REAL8 *mprime_epsilon, REAL8 fHz, INT4 mprime, const REAL8 twopi_Msec, IMRPhenomXPrecessionStruct *pPrec);
vector IMRPhenomX_Return_Roots_MSA(const double LNorm, const double JNorm, const IMRPhenomXPrecessionStruct *pPrec);
vector IMRPhenomX_Return_phi_zeta_costhetaL_MSA(const double xi, IMRPhenomXWaveformStruct *pWF, IMRPhenomXPrecessionStruct *pPrec);
int IMRPhenomX_Initialize_MSA_System(IMRPhenomXWaveformStruct *pWF, IMRPhenomXPrecessionStruct *pPrec, const int ExpansionOrder);
vector IMRPhenomX_Return_MSA_Corrections_MSA(double v, const double LNorm, const double JNorm, const IMRPhenomXPrecessionStruct *pPrec);
double IMRPhenomX_Return_zeta_MSA(const double v, const IMRPhenomXPrecessionStruct *pPrec);
double IMRPhenomX_Return_phiz_MSA(const double v, const double JNorm, const IMRPhenomXPrecessionStruct *pPrec);
double IMRPhenomX_psiofv(const double v, const double v2, const double psi0, const double psi1, const double psi2, const IMRPhenomXPrecessionStruct *pPrec);
double IMRPhenomX_Spin_Evolution_Equation_MSA(IMRPhenomXWaveformStruct *pWF, IMRPhenomXPrecessionStruct *pPrec);


/* ~~~~~~~~~~ Spin Couplings for post-Newtonian Orbital Angular Momentum ~~~~~~~~~~ */
double IMRPhenomX_Get_PN_beta(const double a, const double b, const IMRPhenomXPrecessionStruct *pPrec);
double IMRPhenomX_Get_PN_sigma(const double a, const double b, const IMRPhenomXPrecessionStruct *pPrec);
double IMRPhenomX_Get_PN_tau(const double a, const double b, const IMRPhenomXPrecessionStruct *pPrec);




/* ~~~~~~~~~~ Vector Utility Functions ~~~~~~~~~~ */
/* Function to calculate dot product: vout = v1 . v2 */
double IMRPhenomX_vector_dot_product(const vector v1, const vector v2);
/* Function to calculate cross product: vout = v1 x v2 */
vector IMRPhenomX_vector_cross_product(const vector v1, const vector v2);
/* Function to calculate L2 norm: vout = |v1| */
double IMRPhenomX_vector_L2_norm(const vector v1);
/* Function to perform scalar multiplications: vout = a * v1 */
vector IMRPhenomX_vector_scalar(const vector v1, const double a);
/* Vector addition: vout = v1 + v2 */
vector IMRPhenomX_vector_sum(const vector v1, const vector v2);
/* Vector subtraction: vout = v1 - v2 */
vector IMRPhenomX_vector_diff(const vector v1, const vector v2);

/* Function to transform polar to Cartesian */
vector IMRPhenomX_vector_PolarToCartesian(const sphpolvector v1);
vector IMRPhenomX_vector_PolarToCartesian_components(const REAL8 mag, const REAL8 theta, const REAL8 phi);

/* Function to transform Cartesian to polar */
sphpolvector IMRPhenomX_vector_CartesianToPolar(const vector v1);
/* Function to rotate vector about z axis by given angle */
vector IMRPhenomX_vector_rotate_z(const REAL8 angle, const vector v1);
/* Function to rotate vector about y axis by given angle */
vector IMRPhenomX_vector_rotate_y(const REAL8 angle, const vector v1);
/* Function to rotate vector about z axis by given angle */
void IMRPhenomX_rotate_z(const REAL8 angle, REAL8 *vx, REAL8 *vy, REAL8 *vz);
/* Function to rotate vector about y axis by given angle */
void IMRPhenomX_rotate_y(const REAL8 angle, REAL8 *vx, REAL8 *vy, REAL8 *vz);

/* Vector Utility Functions */
REAL8 IMRPhenomX_Cartesian_to_SphericalPolar_theta(const double x, const double y, const double z);
REAL8 IMRPhenomX_Cartesian_to_SphericalPolar_phi(const double x, const double y, const double z);
void IMRPhenomX_CartesianToPolar(REAL8 *polar,REAL8 *azimuthal,REAL8 *magnitude,REAL8 x,REAL8 y,REAL8 z);
vector IMRPhenomXCreateSphere(const double r, const double th, const double ph);



#ifdef __cplusplus
}
#endif

#endif /* _LALSIM_IMR_PHENOMX_PRECESSION_H */
