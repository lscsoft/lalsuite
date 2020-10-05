
/**
 * This file is part of TEOBResumS
 *
 * Copyright (C) 2017-2018 Alessandro Nagar, Sebastiano Bernuzzi,
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

#ifndef _LALSIM_TEOB_RESUMS_H
#define _LALSIM_TEOB_RESUMS_H

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * @file LALSimTEOBResumS.h
 * @brief Header file of the TEOBResumS C code
 *
 * This file contains all the macros, typdef, and routine prototype.
 * Doxygen documentation should go here.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include <lal/Date.h>
#include <lal/LALString.h>
#include <lal/LALDatatypes.h>
#include <lal/Units.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/SphericalHarmonics.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/Sequence.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeFreqFFT.h>
#include "LALSimSphHarmSeries.h"
#include "LALSimSphHarmMode.h"
#include "LALSimIMR.h"


#ifndef POSTPOSTCIRCULAR
#define POSTPOSTCIRCULAR 1 /* use post-post-circular initial conditions by default */
#endif

#ifndef RC_EXCLUDESPINSPINTIDES
#define RC_EXCLUDESPINSPINTIDES 0 /* use tidally deformed centr. radius with self-spin and tides by default */
#endif

#ifndef USEBTIDALPOTENTIAL
#define USEBTIDALPOTENTIAL 1 /* add B LO tidal potential */
#endif

#ifndef ODESTOPAFTERNDT
#define ODESTOPAFTERNDT (4) /* UNUSED: stop ODE integration after this many timesteps beyond the peak */
#endif

#define MAX_ITER (200)
#define TOLERANCE (1e-14)

/* Macros for postadiabatic stage and ID */
#define POSTADIABATIC_RMIN_BNS (14)
#define POSTADIABATIC_RMIN_BBH (14)
#define POSTADIABATIC_DR (0.1)
#define POSTADIABATIC_NSTEP_MIN (10)
#define POSTADIABATIC_N (8)
#define TEOB_R0_THRESHOLD (14)

#define TEOB_LAMBDA_TOL (1.0)

#define TEOB_GE_TIDES_DEFAULT LAL_SIM_INSPIRAL_GETIDES_GSF23
#define TEOB_GM_TIDES_DEFAULT LAL_SIM_INSPIRAL_GMTIDES_PN

#define TEOB_MODES_BNS_DEFAULT TEOB_MODES_22
#define TEOB_MODES_BBH_DEFAULT TEOB_MODES_22
#define TEOB_MODES_NSBH_DEFAULT TEOB_MODES_22
#define TEOB_MODES_BHNS_DEFAULT TEOB_MODES_22

#define RINGDOWN_EXTEND_ARRAY (1500)
// #define DT_MERGER_INTERP (0.5)

#define INTERP_UNIFORM_GRID_DEFAULT INTERP_UNIFORM_GRID_HLM

/** Macros */
#ifndef DEBUG
#define DEBUG 0
#endif

#define TEOB_STRLEN (128) /** Standard string length */
#define TEOBResumS_Info "TEOBResumS code (C) 2017\n"
#define TEOBResumS_Usage(x) {printf("%sUSAGE:\t%s <parfile>\n", TEOBResumS_Info, x);}
#define SIGN(x,y) ((y) >= 0.0 ? fabs(x) : -fabs(x))
#define typeof __typeof__
#define MAX(a,b)                \
  ({ typeof (a) _a = (a);       \
    typeof (b) _b = (b);        \
    _a > _b ? _a : _b; })
#define MIN(a,b)                \
  ({ typeof (a) _a = (a);       \
    typeof (b) _b = (b);        \
    _a < _b ? _a : _b; })
#define MAX3(a,b,c) (((a) > (b)) ? MAX(a,c) : MAX(b,c))
#define MIN3(a,b,c) (((a) < (b)) ? MIN(a,c) : MIN(b,c))
#define SQ(a) ((a)*(a))
#define DEQUAL(a,b,eps) (fabs((a)-(b))<(eps)) /** double compare */
#define DUNEQUAL(a,b,eps) (fabs((a)-(b))>(eps))
#define STREQUAL(s,t) ((strcmp((s),(t))==0)) /** string compare */
#define SWAPTRS(a,b)   \
  ({		    \
    typeof(a) temp; \
    temp = a;	    \
    a = b;	    \
    b = temp;	    \
  })

/* Useful constants */
// LAL_SQRT2
// #define Sqrt2 (1.41421356237309504880168872420969808)
#define Sqrt3 (1.73205080756887729352744634150587237)
// LAL_SQRT1_2
//#define ooSqrt2 (0.707106781186547524400844362104849039284836)
#define Log1  (0.)
// LAL_LN2
#define Log2  (0.693147180559945309417232)
#define Log3  (1.09861228866810969139525)
#define Log4  (1.38629436111989061883446)
#define Log5  (1.60943791243410037460076)
#define Log6  (1.79175946922805500081248)
#define Log7  (1.94591014905531330510535)
//#define EulerGamma_Log2 (1.27036284546147817002374) /** EulerGamma + Log2 */

/** Index list of EOB evolved variables */
enum{
    TEOB_EVOLVE_RAD,
    TEOB_EVOLVE_PHI,
    TEOB_EVOLVE_PRSTAR,
    TEOB_EVOLVE_PPHI,
    TEOB_EVOLVE_NVARS
};
// static const char* LALTEOBResumSEvolutionVariables[] = {"r","phi","Prstar","Pphi"};

/** Index list of EOB variables for initial data */
enum{
    TEOB_ID_RAD,
    TEOB_ID_PHI,
    TEOB_ID_PPHI,
    TEOB_ID_PRSTAR,
    TEOB_ID_PR,
    TEOB_ID_J,
    TEOB_ID_E0,
    TEOB_ID_OMGJ,
    TEOB_ID_NVARS
};
// static const char* LALTEOBResumSVariables[] = {"r","phi","Pphi","Prstar","Pr","j","E0","Omega"};

/** Index list of EOB dynamical variables (to be stored in arrays) */
enum{
    TEOB_RAD,
    TEOB_PHI,
    TEOB_PPHI,
    TEOB_MOMG,
    TEOB_DDOTR,
    TEOB_PRSTAR,
    TEOB_OMGORB,
    TEOB_E0,
    TEOB_DYNAMICS_NVARS
};

/** Index list of TEOB mode-list options */
enum{
    TEOB_MODES_ALL,
    TEOB_MODES_22,
    TEOB_MODES_NOPT
};

// static const char* LALTEOBResumSDynamicalVariables[] = {"r","phi","Pphi","MOmega","ddor","Prstar","MOmega_orb","E"};

#define KMAX (35) /** Multipolar linear index, max value */
#define PMTERMS_eps (1) /** Switch on Fujita-Iyer point-mass terms. This is hard-coded here */

#define OUTPUT_LM_LEN (35)


/** List of options for centrifugal radius */
enum{
    RC_LO,
    RC_NLO,
    RC_NNLO,
    RC_NNLOS4,
    RC_NOSPIN,
    RC_NOTIDES,
    RC_NOPT
};

enum{
    SS_LO,
    SS_NLO,
    SS_NOPT
};

/** List of options for ODE timestepping */
enum{
    ODE_TSTEP_UNIFORM,
    ODE_TSTEP_ADAPTIVE,
    ODE_TSTEP_ADAPTIVE_UNIFORM_AFTER_LSO,
    ODE_TSTEP_NOPT
};
static const char* const LALTEOBResumSODETimeStep[] = {"uniform","adaptive","adaptive+uniform_after_LSO","undefined"};

/** List of options for postadiabatic inspiral */
enum{
    TEOB_PA_OFF,
    TEOB_PA_ON,
    TEOB_PA_ONLY,
    TEOB_PA_NOPT
};

/** List of options for interp_uniform_grid */
enum{
    INTERP_UNIFORM_GRID_OFF,
    INTERP_UNIFORM_GRID_HPC,
    INTERP_UNIFORM_GRID_HLM,
    INTERP_UNIFORM_GRID_NOPT
};

/** List of options for NQC waveform nqc_coeffs_hlm and nqc_coeffs_flx */
enum{
    NQC_OFF,
    NQC_COMPUTE,
    NQC_NR_NOSPIN,
    NQC_FILE,
    NQC_NOPT
};

/** Error handler for root finders */
enum{
    ROOT_ERRORS_NO,
    ROOT_ERRORS_BRACKET,
    ROOT_ERRORS_MAXITS,
    ROOT_ERRORS_NOSUCC,
    ROOT_ERRORS
};
static const char* const LALTEOBResumSRootErrors[] = {"none","root is not bracketed.","root finder did not converged.","root finder failed."};
#define ROOTFINDER(i, x) {if ( ((i) = (x)) && ((i)>ROOT_ERRORS_NO) )  { exit(LALTEOBResumSRootErrors[(i)]); }} //TODO: CHECK THIS MACRO (LOGIC INVOLVED)

/** Maps between linear index and the corresponding (l, m) multipole indices */
extern const INT4 TEOB_LINDEX[KMAX]; /* defined in LALSimIMRTEOBResumS.c */
extern const INT4 TEOB_MINDEX[KMAX]; /* defined in LALSimIMRTEOBResumS.c */


/** Multipolar waveform at given time, comes at handy */
typedef struct tagLALTEOBResumSWaveformModeSingleTime
{
    REAL8 time;
    REAL8 ampli[KMAX]; /* amplitude */
    REAL8 phase[KMAX]; /* phase */
    /* kmask is currently unused in the LAL version */
    // INT4 kmask[KMAX]; /* mask for multipoles */
}  LALTEOBResumSWaveformModeSingleTime;

/** Func pointer types */
typedef void (*EOBWavFlmSFunc)(REAL8, REAL8, REAL8, REAL8, REAL8, REAL8, REAL8, REAL8, REAL8, REAL8, REAL8 [KMAX][6], int, REAL8 *, REAL8 *);
typedef void (*EOBDynSGetRCFunc)(REAL8, REAL8, REAL8, REAL8, REAL8, REAL8, REAL8, REAL8, REAL8, REAL8, REAL8, INT4, REAL8*, REAL8*, REAL8*);

/** Multipolar coefficients for NQC waveform */
typedef struct tagNQCcoefs
{
    double a1[KMAX];
    double a2[KMAX];
    double a3[KMAX];
    double b1[KMAX];
    double b2[KMAX];
    double b3[KMAX];
    double n[KMAX][6];
    int activemode[KMAX]; /* mask for modes with nonzero NQC */
    int maxk; /* index of maximum active mode */
    int add; /* flag */
} NQCcoefs;

/** NQC data for flux and waveform */
typedef struct tagNQCdata
{
    NQCcoefs *flx;
    NQCcoefs *hlm;
} NQCdata;

/** Dynamics data type */
typedef struct tagLALTEOBResumSDynamics
{
    char name[TEOB_STRLEN];
    /* various pointwise variables */
    INT4 store; /* store following values? */
    INT4 noflux; /* compute rhs without flux */
    REAL8 t, r, phi, pphi, prstar, ddotr, Omg, Omg_orb;
    REAL8 H, Heff, Heff_orb, E, jhat, r_omega, psi, v_phi;
    REAL8 A, dA, d2A, B, dB;
    REAL8 MOmg, MOmg_prev, tMOmgpeak;
    /* stuff for ODE solver */
    REAL8 y[TEOB_EVOLVE_NVARS]; /* rhs storage */
    REAL8 dy[TEOB_EVOLVE_NVARS];
    REAL8 y0[TEOB_ID_NVARS]; /* ID storage */
    REAL8 dt, t_stop, ti;
    REAL8 r_stop;
    INT4 ode_timestep;
    bool ode_stop, ode_stop_MOmgpeak, ode_stop_radius;
    INT4 nqc_coeffs_hlm, nqc_coeffs_flx;
    NQCdata *NQC;
    REAL8 clm[KMAX][6];
    /* arrays */
    INT4 size;
    REAL8 *time;
    REAL8 *data[TEOB_DYNAMICS_NVARS];
    /* func pointers */
    EOBDynSGetRCFunc eob_dyn_s_get_rc;
    EOBWavFlmSFunc eob_wav_flm_s;
    /* key parameters for quick access */
    REAL8 M, nu, q, X1, X2;
    REAL8 chi1, chi2, S1, S2, S, Sstar, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, cN3LO;
    REAL8 rLR, rLSO;
    REAL8 kapA2,kapA3,kapA4, kapB2,kapB3,kapB4, kapT2,kapT3,kapT4, khatA2,khatB2;
    REAL8 bar_alph2_1, bar_alph2_2, bar_alph3_1, bar_alph3_2, bar_alph2j_1;
    REAL8 kapA2j, kapB2j, kapT2j;
    REAL8 rLR_tidal, pGSF_tidal;
    REAL8 Mbhf, abhf; /* final BH */
    INT4 use_tidal, use_spins, use_tidal_gravitomagnetic;
    INT4 bhns; /* 0:BNS, 1:BHNS, 2:NSBH */
} LALTEOBResumSDynamics;

/* Function protoypes grouped based on file */

/** Should not be needed in LAL (all parameters passed to interface function) */
/* TEOBResumSPars.c */
//void par_db_init (void);
//void par_db_free (void);
//void par_db_default (void);
//void par_file_parse (const char *fname);
//void par_file_parse_merge (const char *fname);
//void par_db_write_file (const char *fname);
//void par_db_screen (void);
//void par_set_i(const char *key, int val);
//void par_set_b(const char *key, int val);
//void par_set_d(const char *key, double val);
//void par_set_s(const char *key, const char *val);
//void par_set_arrayi (const char *key, int *array, int n);
//void par_set_arrayd (const char *key, double *array, int n);
//int par_get_i(const char *key);
//int par_get_b(const char *key);
//double par_get_d(const char *key);
//const char * par_get_s(const char *key);
//int * par_get_arrayi(const char *key, int *n);
//double * par_get_arrayd(const char *key, int *n);
//void eob_set_params(char *s, int n);
//void eob_free_params(void);

INT4 XLALSetup_TEOB_mode_array(LALValue *ModeArray, INT4 modeType);
INT4 XLALCheck_TEOB_mode_array(LALValue *ModeArray, UINT4 use_tidal);

/* TEOBResumSUtil.c */
REAL8 q_to_nu(const REAL8 q);
REAL8 nu_to_X1(const REAL8 nu);
REAL8 Eulerlog(const REAL8 x,const INT4 m);
void interp_spline(REAL8 *t, REAL8 *y, INT4 n, REAL8 *ti, INT4 ni, REAL8 *yi);
int find_point_bisection(REAL8 x, INT4 n, REAL8 *xp, INT4 o);
REAL8 baryc_f(REAL8 xx, INT4 n, REAL8 *f, REAL8 *x);
void baryc_weights(INT4 n, REAL8 *x, REAL8 *omega);
REAL8 baryc_f_weights(REAL8 xx, INT4 n, REAL8 *f, REAL8 *x, REAL8 *omega);
REAL8 interp1d (const INT4 order, REAL8 xx, INT4 nx, REAL8 *f, REAL8 *x);
REAL8 find_max (const INT4 n, REAL8 dx, REAL8 x0, REAL8 *f, REAL8 *fmax);
//double fact(int n);
/** XLALWignerdMatrix, XLALWignerDMatrix in <lal/SphericalHarmonics.h> */
// double wigner_d_function(int l, int m, int s, double i);
//int spinsphericalharm(double *rY, double *iY, int s, int l, int m, double phi, double i);

//void compute_hpc(SphHarmPolarTimeSeries *hlm, double nu, double M, double distance, double amplitude_prefactor, double psi, double iota, COMPLEX16TimeSeries *hpc);

int D0(REAL8 *f, REAL8 dx, INT4 n, REAL8 *df);
int D2(REAL8 *f, REAL8 dx, INT4 n, REAL8 *d2f);
//int D0_x(double *f, double *x, int n, double *df);

//double cumtrapz(double *f, double *x, const int n, double *sum);
REAL8 cumint3(REAL8 *f, REAL8 *x, const INT4 n, REAL8 *sum);

void unwrap(REAL8 *p, const INT4 size);
void unwrap_proxy(REAL8 *p, REAL8 *r, const INT4 size, const INT4 shift0);

int get_uniform_size(const REAL8 tf, const REAL8 t0, const REAL8 dt);

/** Replaced by LAL COMPLEX16TimeSeries functions */
//void Waveform_alloc (Waveform **wav, int size, const char *name);
//void Waveform_push (Waveform **wav, int size);
//void Waveform_output (Waveform *wav);
//void Waveform_free (Waveform *wav);
//void Waveform_lm_alloc (LALTEOBResumSWaveformMode **wav, int size, const char *name);
//void Waveform_lm_push (LALTEOBResumSWaveformMode **wav, int size);
//void Waveform_lm_output (LALTEOBResumSWaveformMode *wav);
//void Waveform_lm_output_reim (LALTEOBResumSWaveformMode *wav);
//void Waveform_lm_free (LALTEOBResumSWaveformMode *wav);
void Waveform_lm_t_alloc (LALTEOBResumSWaveformModeSingleTime **wav);
void Waveform_lm_t_free (LALTEOBResumSWaveformModeSingleTime *wav);

void XLALTEOBDynamicsInit (LALTEOBResumSDynamics **dyn, INT4 size, const CHAR *name);
void XLALTEOBDynamicsPush (LALTEOBResumSDynamics **dyn, INT4 size);
void XLALFreeTEOBDynamics (LALTEOBResumSDynamics *dyn);
void XLALTEOBDynamicsSetParams(LALTEOBResumSDynamics *dyn,
                               LALValue **pModeArray,
                               const REAL8 m1,
                               const REAL8 m2,
                               const REAL8 S1x,
                               const REAL8 S1y,
                               const REAL8 S1z,
                               const REAL8 S2x,
                               const REAL8 S2y,
                               const REAL8 S2z,
                               const REAL8 deltaT,
                               const REAL8 lambda1,
                               const REAL8 lambda2,
                               const REAL8 lambda1oct,
                               const REAL8 lambda2oct,
                               const REAL8 lambda1hex,
                               const REAL8 lambda2hex,
                               LALDict *LALparams);

void NQCdata_alloc (NQCdata **nqc);
void NQCdata_free (NQCdata *nqc);

REAL8 time_units_factor(REAL8 M);
REAL8 time_units_conversion(REAL8 M, REAL8 t);
REAL8 radius0(REAL8 M, REAL8 fHz);
REAL8 get_mrg_timestep(REAL8 q, REAL8 chi1, REAL8 chi2);

void XLALSimIMRComputePolarisations(REAL8Sequence *hplus_out, REAL8Sequence *hcross_out, SphHarmTimeSeries *hlm, LALValue *modeArray, REAL8 amplitude_prefactor, REAL8 theta, REAL8 phi);

void XLALSimIMRComputePolarisationsPolar(REAL8Sequence *hplus_out, REAL8Sequence *hcross_out, SphHarmPolarTimeSeries *hlm, LALValue *modeArray, REAL8 amplitude_prefactor, REAL8 theta, REAL8 phi);

/* TEOBResumSFits.c */
//double eob_c3_fit_global(double nu, double chi1, double chi2, double X1, double X2, double a1, double a2);
//double eob_nqc_dtfit(const double chi, const double chi0);
REAL8 eob_nqc_timeshift(REAL8 nu, REAL8 chi1);
REAL8 logQ(REAL8 x);
REAL8 Yagi13_fit_barlamdel(REAL8 barlam2, INT4 ell);
REAL8 Yagi13_fit_barsigmalambda(REAL8 barlam2);
REAL8 Yagi14_fit_Coct(REAL8 C_Q);
REAL8 Yagi14_fit_Chex(REAL8 C_Q);
REAL8 JFAPG_fit_Sigma_Static(REAL8 barlam2);
REAL8 JFAPG_fit_Sigma_Irrotational(REAL8 barlam2);
void HealyBBHFitRemnant(REAL8 chi1, REAL8 chi2, REAL8 q, REAL8 *mass, REAL8 *spin);
REAL8 JimenezFortezaRemnantSpin(REAL8 nu, REAL8 X1, REAL8 X2, REAL8 chi1, REAL8 chi2);

void QNMHybridFitCab(REAL8 nu, REAL8 X1, REAL8 X2, REAL8 chi1, REAL8 chi2, REAL8 aK, REAL8 Mbh, REAL8 abh, REAL8 *a1, REAL8 *a2, REAL8 *a3, REAL8 *a4, REAL8 *b1, REAL8 *b2, REAL8 *b3, REAL8 *b4, REAL8 *sigmar, REAL8 *sigmai, INT4 usespins);
//double eob_approxLR(const double nu);

/* TEOBResumSDynamics.c */
int eob_dyn_rhs(REAL8 t, const REAL8 y[], REAL8 dy[], void *params);
void eob_ham(REAL8 nu, REAL8 r, REAL8 pph, REAL8 prstar, REAL8 A, REAL8 dA,
	     REAL8 *H, REAL8 *Heff, REAL8 *dHeff_dr, REAL8 *dHeff_dprstar, REAL8 *dHeff_dpphi);
int eob_dyn_rhs_s(REAL8 t, const REAL8 y[], REAL8 dy[], void *params);
void eob_ham_s(REAL8 nu, REAL8 r, REAL8 rc, REAL8 drc_dr, REAL8 pphi, REAL8 prstar, REAL8 S, REAL8 Sstar, REAL8 chi1, REAL8 chi2, REAL8 X1, REAL8 X2, REAL8 aK2, REAL8 c3, REAL8 A, REAL8 dA, REAL8 *H, REAL8 *Heff, REAL8 *Heff_orb, REAL8 *dHeff_dr, REAL8 *dHeff_dprstar, REAL8 *dHeff_dpphi, REAL8 *d2Heff_dprstar20);
void eob_dyn_s_GS(REAL8 r, REAL8 rc, REAL8 drc_dr, REAL8 aK2, REAL8 prstar, REAL8 pph, REAL8 nu, REAL8 chi1, REAL8 chi2, REAL8 X1, REAL8 X2, REAL8 cN3LO, REAL8 *ggm);
//void eob_dyn_s_get_rc(REAL8 r, REAL8 nu, REAL8 at1,REAL8 at2, REAL8 aK2, REAL8 C_Q1, REAL8 C_Q2, INT4 usetidal, REAL8 *rc, REAL8 *drc_dr, REAL8 *d2rc_dr2);

void eob_dyn_s_get_rc_LO(REAL8 r, REAL8 nu, REAL8 at1, REAL8 at2, REAL8 aK2, REAL8 C_Q1, REAL8 C_Q2, REAL8 UNUSED C_Oct1, REAL8 UNUSED C_Oct2, REAL8 UNUSED C_Hex1, REAL8 UNUSED C_Hex2, INT4 usetidal, REAL8 *rc, REAL8 *drc_dr, REAL8 *d2rc_dr2);
void eob_dyn_s_get_rc_NLO(REAL8 r, REAL8 nu, REAL8 at1, REAL8 at2, REAL8 aK2, REAL8 C_Q1, REAL8 C_Q2, REAL8 UNUSED C_Oct1, REAL8 UNUSED C_Oct2, REAL8 UNUSED C_Hex1, REAL8 UNUSED C_Hex2, INT4 usetidal, REAL8 *rc, REAL8 *drc_dr, REAL8 *d2rc_dr2);
void eob_dyn_s_get_rc_NNLO(REAL8 r, REAL8 nu, REAL8 at1, REAL8 at2, REAL8 aK2, REAL8 C_Q1, REAL8 C_Q2, REAL8 UNUSED C_Oct1, REAL8 UNUSED C_Oct2, REAL8 UNUSED C_Hex1, REAL8 UNUSED C_Hex2, INT4 usetidal, REAL8 *rc, REAL8 *drc_dr, REAL8 *d2rc_dr2);
void eob_dyn_s_get_rc_NNLO_S4(REAL8 r, REAL8 nu, REAL8 at1,REAL8 at2, REAL8 aK2, REAL8 C_Q1, REAL8 C_Q2, REAL8 C_Oct1, REAL8 C_Oct2, REAL8 C_Hex1, REAL8 C_Hex2, INT4 usetidal, REAL8 *rc, REAL8 *drc_dr, REAL8 *d2rc_dr2);
void eob_dyn_s_get_rc_NOSPIN(REAL8 r, REAL8 UNUSED nu, REAL8 UNUSED at1, REAL8 UNUSED at2, REAL8 UNUSED aK2, REAL8 UNUSED C_Q1, REAL8 UNUSED C_Q2, REAL8 UNUSED C_Oct1, REAL8 UNUSED C_Oct2, REAL8 UNUSED C_Hex1, REAL8 UNUSED C_Hex2, INT4 UNUSED usetidal, REAL8 UNUSED *rc, REAL8 UNUSED *drc_dr, REAL8 UNUSED *d2rc_dr2);
void eob_dyn_s_get_rc_NOTIDES(REAL8 r, REAL8 nu, REAL8 at1,REAL8 at2, REAL8 aK2, REAL8 C_Q1, REAL8 C_Q2, REAL8 UNUSED C_Oct1, REAL8 UNUSED C_Oct2, REAL8 UNUSED C_Hex1, REAL8 UNUSED C_Hex2, INT4 usetidal, REAL8 *rc, REAL8 *drc_dr, REAL8 *d2rc_dr2);
REAL8 eob_dyn_fLR(REAL8 r, void * params);
int eob_dyn_adiabLR(LALTEOBResumSDynamics *dyn, REAL8 *rLR, INT4 tidesFlag);
//double eob_dyn_fLSO(double r, void * params);
//int eob_dyn_adiabLSO(LALTEOBResumSDynamics *dyn, REAL8 *rLSO);

/* TEOBResumSPostAdiabatic.c */
int eob_dyn_Npostadiabatic(LALTEOBResumSDynamics *dyn, REAL8 r0);

/* TEOBResumSInitialCondition.c */
void eob_dyn_ic(REAL8 r0, LALTEOBResumSDynamics *dyn, REAL8 y_init[]);
void eob_dyn_ic_s(REAL8 r0, LALTEOBResumSDynamics *dyn, REAL8 y_init[]);
REAL8 eob_dyn_bisecHeff0_s(REAL8 nu, REAL8 chi1, REAL8 chi2, REAL8 X1, REAL8 X2, REAL8 c3, REAL8 pph, REAL8 rorb, REAL8 A, REAL8 dA, REAL8 rc, REAL8 drc_dr, REAL8 ak2, REAL8 S, REAL8 Ss);
REAL8 eob_dyn_DHeff0(REAL8 x, void *params);

/* TEOBResumSMetric.c */
void eob_metric_A5PNlog(REAL8 r, REAL8 nu, REAL8 *A, REAL8 *dA, REAL8 *d2A);
void eob_metric_Atidal(REAL8 r, LALTEOBResumSDynamics *dyn, REAL8 *AT, REAL8 *dAT, REAL8 *d2AT);
void eob_metric(REAL8 r, LALTEOBResumSDynamics *dyn, REAL8 *A, REAL8 *B, REAL8 *dA, REAL8 *d2A, REAL8 *dB);
void eob_metric_s(REAL8 r, LALTEOBResumSDynamics *dyn, REAL8 *A, REAL8 *B, REAL8 *dA, REAL8 *d2A, REAL8 *dB);

/* TEOBResumSFlux.c */
REAL8 eob_flx_Flux(REAL8 x, REAL8 Omega, REAL8 r_omega, REAL8 E, REAL8 Heff, REAL8 jhat, REAL8 r, REAL8 pr_star, REAL8 ddotr, LALTEOBResumSDynamics *dyn);
REAL8 eob_flx_Flux_s(REAL8 x, REAL8 Omega, REAL8 r_omega, REAL8 E, REAL8 Heff, REAL8 jhat, REAL8 r, REAL8 pr_star, REAL8 ddotr, LALTEOBResumSDynamics *dyn);
void eob_flx_Tlm(REAL8 w, REAL8 *MTlm);
void eob_flx_FlmNewt(REAL8 x, REAL8 nu, REAL8 *Nlm);
REAL8 eob_flx_HorizonFlux(REAL8 x, REAL8 Heff, REAL8 jhat, REAL8 nu);
REAL8 eob_flx_HorizonFlux_s(REAL8 x, REAL8 Heff, REAL8 jhat, REAL8 nu, REAL8 X1, REAL8 X2, REAL8 chi1, REAL8 chi2);

/* TEOBResumSWaveform.c */
void eob_wav_hlm(LALTEOBResumSDynamics *dyn, LALTEOBResumSWaveformModeSingleTime *hlm_t);

void eob_wav_deltalm(REAL8 Hreal,REAL8 Omega,REAL8 nu, REAL8 *dlm);
void eob_wav_hhatlmTail(REAL8 Omega,REAL8 Hreal,REAL8 bphys, LALTEOBResumSWaveformModeSingleTime *tlm);
void eob_wav_speedyTail(REAL8 Omega, REAL8 Hreal, REAL8 bphys, LALTEOBResumSWaveformModeSingleTime *tlm);
void eob_wav_hlmNewt(REAL8 r, REAL8 Omega, REAL8 phi, REAL8 nu, LALTEOBResumSWaveformModeSingleTime *hNewt);
void eob_wav_hlmTidal(REAL8 x, LALTEOBResumSDynamics *dyn, REAL8 *hTidallm);
void eob_wav_flm_coeffs(REAL8 nu, REAL8 clm[KMAX][6]);
void eob_wav_flm(REAL8 x,REAL8 nu, REAL8 clm[KMAX][6], REAL8 *rholm, REAL8 *flm);
//void eob_wav_flm_old(double x,double nu, double *rholm, double *flm);

void eob_wav_flm_s_SSLO(REAL8 x, REAL8 nu, REAL8 X1, REAL8 X2, REAL8 chi1, REAL8 chi2, REAL8 a1, REAL8 a2, REAL8 C_Q1, REAL8 C_Q2, REAL8 clm[KMAX][6], int usetidal, REAL8 *rholm, REAL8 *flm);
void eob_wav_flm_s_SSNLO(REAL8 x, REAL8 nu, REAL8 X1, REAL8 X2, REAL8 chi1, REAL8 chi2, REAL8 a1, REAL8 a2, REAL8 C_Q1, REAL8 C_Q2, REAL8 clm[KMAX][6], int usetidal, REAL8 *rholm, REAL8 *flm);
//void eob_wav_flm_s_old(double x, double nu, double X1, double X2, double chi1, double chi2, double a1, double a2, double C_Q1, double C_Q2, int usetidal, double *rholm, double *flm);
void eob_wav_hlmNQC_find_a1a2a3(LALTEOBResumSDynamics *dyn, SphHarmPolarTimeSeries *h, SphHarmPolarTimeSeries *hnqc);
void eob_wav_hlmNQC_find_a1a2a3_mrg(LALTEOBResumSDynamics *dyn_mrg, SphHarmPolarTimeSeries *hlm_mrg, SphHarmPolarTimeSeries *hnqc, LALTEOBResumSDynamics *dyn, SphHarmPolarTimeSeries *hlm);
void eob_wav_hlmNQC(REAL8  nu, REAL8  r, REAL8  prstar, REAL8  Omega, REAL8  ddotr, NQCcoefs *nqc, LALTEOBResumSWaveformModeSingleTime *psilmnqc);
void eob_wav_ringdown_template(REAL8 x, REAL8 a1, REAL8 a2, REAL8 a3, REAL8 a4, REAL8 b1, REAL8 b2, REAL8 b3, REAL8 b4, REAL8 sigmar, REAL8 sigmai, REAL8 *psi);
void eob_wav_ringdown(LALTEOBResumSDynamics *dyn, SphHarmPolarTimeSeries *hlm);

REAL8 eob_dyn_bisecOmegaorb0(LALTEOBResumSDynamics *dyn, REAL8 omg_orb0, REAL8 r0_kepl);
REAL8 eob_dyn_r0_Kepler (REAL8 f0);
REAL8 eob_dyn_r0_eob (REAL8 f0, LALTEOBResumSDynamics *dyn);
REAL8 eob_dyn_Omegaorb0(REAL8 r, void *params);
void unwrap(REAL8 *p, const INT4 size);
INT4 get_uniform_size(const REAL8 tN, const REAL8 t0, const REAL8 dt);
void XLALTEOBDynamicsInterp (LALTEOBResumSDynamics *dyn, const INT4 size, const REAL8 t0, const REAL8 dt, const char *name);
void XLALTEOBDynamicsExtract (LALTEOBResumSDynamics *dyna, const REAL8 to, const REAL8 tn, LALTEOBResumSDynamics **dynb, const char *name);
void XLALTEOBDynamicsJoin (LALTEOBResumSDynamics *dyna, LALTEOBResumSDynamics *dynb, REAL8 to);
void eob_wav_hlmNQC_nospin201602(REAL8  nu,
                                 REAL8  r,
                                 REAL8  prstar,
                                 REAL8  Omega,
                                 REAL8  ddotr,
                                 LALTEOBResumSWaveformModeSingleTime *hlmnqc);
void eob_nqc_setcoefs(LALTEOBResumSDynamics *dyn);
void eob_nqc_setcoefs_nospin201602(REAL8 nu, NQCcoefs *nqc);
REAL8 eob_c3_fit_global(REAL8 nu, REAL8 chi1, REAL8 chi2, REAL8 X1, REAL8 X2, REAL8 a1, REAL8 a2);

void XLALSphHarmPolarJoin (SphHarmPolarTimeSeries *hlma, SphHarmPolarTimeSeries *hlmb, REAL8 to);


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif // of #ifndef _LALSIM_TEOB_RESUMS_H
