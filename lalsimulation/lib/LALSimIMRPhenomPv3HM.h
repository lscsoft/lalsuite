#ifndef _LALSIM_IMR_PHENOMPv3HM_H
#define _LALSIM_IMR_PHENOMPv3HM_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <lal/LALSimIMR.h>
#include <lal/LALDict.h>

/* IMRPhenomPv3 - uses the angles from arXiv 1703.03967*/
#include "LALSimInspiralFDPrecAngles_internals.c"

/* default and constant value places lhat = (0,0,1) */
#define LHAT_COS_THETA 1.0 /* Cosine of Polar angle of orbital angular momentum */
#define LHAT_PHI 0.0       /* Azimuthal angle of orbital angular momentum */

/* structs */

/**
 * Structure storing initial and derived variables for IMRPhenomPv3HM
 */
typedef struct tagPhenomPv3HMStorage
{
    INT4 PRECESSING;   /**< integer to signify if system is precessing, 1 for false (not precessing), 0 for true (precessing) */
    /* input parameters */
    REAL8 m1_SI;       /**< mass of primary in SI (kg) */
    REAL8 m2_SI;       /**< mass of secondary in SI (kg) */
    REAL8 chi1x;       /**< x-component of dimensionless spin on primary w.r.t. Lhat = (0,0,1) */
    REAL8 chi1y;       /**< y-component of dimensionless spin on primary w.r.t. Lhat = (0,0,1) */
    REAL8 chi1z;       /**< z-component of dimensionless spin on primary w.r.t. Lhat = (0,0,1) */
    REAL8 chi2x;       /**< x-component of dimensionless spin on secondary w.r.t. Lhat = (0,0,1) */
    REAL8 chi2y;       /**< y-component of dimensionless spin on secondary w.r.t. Lhat = (0,0,1) */
    REAL8 chi2z;       /**< z-component of dimensionless spin on secondary w.r.t. Lhat = (0,0,1) */
    REAL8 distance_SI; /**< distance to source in SI (m) */
    REAL8 inclination; /**< inclination - used to compute the angle thetaJN (rad) */
    REAL8 phiRef;      /**< */
    REAL8 deltaF;      /**< frequency spacing (Hz) */
    REAL8 f_min;       /**< starting GW frequency (Hz) */
    REAL8 f_max;       /**< ending GW frequency (Hz) */
    REAL8 f_ref;       /**< reference GW frequency (Hz) */
    /* derived parameters */
    REAL8 m1_Msun;      /**< mass of primary in solar masses */
    REAL8 m2_Msun;      /**< mass of secondary in solar masses */
    REAL8 Mtot_SI;      /**< total mass in SI (kg) */
    REAL8 Mtot_Msun;    /**< total mass in solar masses */
    REAL8 eta;          /**< Symmetric mass ratio*/
    REAL8 q;            /* with m1>=m2 so q>=1 */
    REAL8 Msec;         /**< Total mass in seconds */
    REAL8 f_ref_Orb_Hz; /**< Reference orbital frequency (Hz) [It's the reference GW frequency converted to orbital frequency] */
    REAL8 twopi_Msec;   /**< LAL_TWOPI * Msec */
    REAL8 amp0;         /**< frequency domain physical scaling */
    /* variables used when rotating input parameters (LAL frame) into PhenomP intrinsic parameters  */
    REAL8 chip;         /**< effective precessing parameter */
    REAL8 thetaJN;      /**< Angle between J0 and line of sight (z-direction) */
    REAL8 alpha0;       /**< Initial value of alpha angle (azimuthal precession angle) */
    REAL8 phi_aligned;  /**< Initial phase to feed the underlying aligned-spin model */
    REAL8 zeta_polariz; /**< Angle to rotate the polarizations */
    /* compute spins in polar coordinates */
    REAL8 chi1_mag;   /**< dimensionless spin magnitude on primary */
    REAL8 chi1_theta; /**< polar angle w.r.t. Lhat = (0,0,1) on primary */
    REAL8 chi1_phi;   /**< azimuthal angle w.r.t. Lhat = (0,0,1) on primary */
    REAL8 chi2_mag;   /**< dimensionless spin magnitude on secondary */
    REAL8 chi2_theta; /**< polar angle w.r.t. Lhat = (0,0,1) on secondary */
    REAL8 chi2_phi;   /**< azimuthal angle w.r.t. Lhat = (0,0,1) on secondary */
    /* Precession angles at reference frequency */
    REAL8 alphaRef;   /**< azimuthal precession angle at f_ref */
    REAL8 epsilonRef; /**< epsilon precession angle at f_ref */
    REAL8 betaRef;    /**< beta (opening angle) precession angle at f_ref */
    // REAL8 t_corr; /**< time shift for peak */
    // REAL8 finspin; /**< final spin */
} PhenomPv3HMStorage;

/* function prototypes */

static LALDict *IMRPhenomPv3HM_setup_mode_array(LALDict *extraParams);

static int init_PhenomPv3HM_Storage(PhenomPv3HMStorage *p, sysq *q, REAL8 m1_SI, REAL8 m2_SI, REAL8 S1x, REAL8 S1y, REAL8 S1z, REAL8 S2x, REAL8 S2y, REAL8 S2z, const REAL8 distance, const REAL8 inclination, const REAL8 phiRef, const REAL8 deltaF, const REAL8 f_min, const REAL8 f_max, const REAL8 f_ref);

static int IMRPhenomPv3HM_check_mode_array(LALValue *ModeArray);

static int IMRPhenomPv3HM_Compute_Mode(COMPLEX16FrequencySeries **hptilde, COMPLEX16FrequencySeries **hctilde, UINT4 ell, INT4 mprime, const REAL8 Mtot_Msun, PhenomPv3HMStorage *pv3HM, SphHarmFrequencySeries **hlmsD, sysq *pAngles, REAL8Sequence *freqs_seq);

static int IMRPhenomPv3HM_Compute_a_b_e(REAL8 *alpha, REAL8 *beta, REAL8 *mprime_epsilon, REAL8 fHz, INT4 mprime, const REAL8 twopi_Msec, PhenomPv3HMStorage *pv3HM, sysq *pAngles);

static int IMRPhenomPv3HM_wigner_loop(COMPLEX16 *Term1, COMPLEX16 *Term2, INT4 ell, INT4 mprime, IMRPhenomPv3HMYlmStruct *ylms, IMRPhenomPv3HMAlphaStruct *als, IMRPhenomPv3HMWignderStruct *wigs);

#ifdef __cplusplus
}
#endif

#endif /* _LALSIM_IMR_PHENOMPv3HM_H */
