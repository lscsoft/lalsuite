#ifndef _LALSIM_IMR_PHENOMHM_H
#define _LALSIM_IMR_PHENOMHM_H

#include <lal/LALDatatypes.h>
#include <lal/Sequence.h>
#include <lal/LALDict.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/FrequencySeries.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __GNUC__
#define UNUSED __attribute__((unused))
#else
#define UNUSED
#endif

/**
 * Dimensionless frequency (Mf) at which define the end of the waveform
 */
#define PHENOMHM_DEFAULT_MF_MAX 0.5

/**
 * eta is the symmetric mass-ratio.
 * This number corresponds to a mass-ratio of 20
 * The underlying PhenomD model was calibrated to mass-ratio 18
 * simulations. We choose mass-ratio 20 as a conservative
 * limit on where the model should be accurate.
 */
#define MAX_ALLOWED_ETA 0.045351

/**
 * Maximum number of (l,m) mode paris PhenomHM models.
 * Only count positive 'm' modes.
 * Used to set default mode array
 */
#define NMODES_MAX 6

/**
 * Highest ell multipole PhenomHM models + 1.
 * Used to set array sizes
 */
#define L_MAX_PLUS_1 5

/* Activates amplitude part of the model */
#define AmpFlagTrue 1
#define AmpFlagFalse 0

LALDict *IMRPhenomHM_setup_mode_array(
    LALDict *extraParams);

static int IMRPhenomHM_check_mode_array(LALValue *ModeArray);

/**
 * useful powers in GW waveforms: 1/6, 1/3, 2/3, 4/3, 5/3, 2, 7/3, 8/3, -1, -1/6, -7/6, -1/3, -2/3, -5/3
 * calculated using only one invocation of 'pow', the rest are just multiplications and divisions
 */
typedef struct tagPhenomHMUsefulPowers
{
    REAL8 third;
    REAL8 two_thirds;
    REAL8 four_thirds;
    REAL8 five_thirds;
    REAL8 two;
    REAL8 seven_thirds;
    REAL8 eight_thirds;
    REAL8 inv;
    REAL8 m_seven_sixths;
    REAL8 m_third;
    REAL8 m_two_thirds;
    REAL8 m_five_thirds;
} PhenomHMUsefulPowers;

/**
  * Useful powers of Mf: 1/6, 1/3, 2/3, 4/3, 5/3, 2, 7/3, 8/3, -7/6, -5/6, -1/2, -1/6, 1/2
  * calculated using only one invocation of 'pow' and one of 'sqrt'.
  * The rest are just multiplications and divisions.  Also including Mf itself in here.
  */
typedef struct tagPhenomHMUsefulMfPowers
{
    REAL8 itself;
    REAL8 sixth;
    REAL8 third;
    REAL8 two_thirds;
    REAL8 four_thirds;
    REAL8 five_thirds;
    REAL8 two;
    REAL8 seven_thirds;
    REAL8 eight_thirds;
    REAL8 m_seven_sixths;
    REAL8 m_five_sixths;
    REAL8 m_sqrt;
    REAL8 m_sixth;
    REAL8 sqrt;
} PhenomHMUsefulMfPowers;

/**
 * must be called before the first usage of *p
 */
int PhenomHM_init_useful_mf_powers(PhenomHMUsefulMfPowers *p, REAL8 number);

/**
 * must be called before the first usage of *p
 */
int PhenomHM_init_useful_powers(PhenomHMUsefulPowers *p, REAL8 number);

/**
  * Structure holding Higher Mode Phase pre-computations
  */
typedef struct tagHMPhasePreComp
{
    double ai;
    double bi;
    double am;
    double bm;
    double ar;
    double br;
    double fi;
    double fr;
    double PhDBconst;
    double PhDCconst;
    double PhDBAterm;
} HMPhasePreComp;

/**
 * Structure storing pre-determined quantities
 * that describe the frequency array
 * and tells you over which indices will contain non-zero values.
 */
typedef struct tagPhenomHMFrequencyBoundsStorage
{
    REAL8 deltaF;
    REAL8 f_min;
    REAL8 f_max;
    REAL8 f_ref;
    UINT4 freq_is_uniform; /**< If = 1 then assume uniform spaced, If = 0 then assume arbitrarily spaced. */
    size_t npts;           /**< number of points in waveform array */
    size_t ind_min;        /**< start index containing non-zero values */
    size_t ind_max;        /**< end index containing non-zero values */
} PhenomHMFrequencyBoundsStorage;

int init_IMRPhenomHMGet_FrequencyBounds_storage(
    PhenomHMFrequencyBoundsStorage *p,
    REAL8Sequence *freqs,
    REAL8 Mtot,
    REAL8 deltaF,
    REAL8 f_ref_in);

UINT4 IMRPhenomHM_is_freq_uniform(
    REAL8Sequence *freqs,
    REAL8 deltaF);

/**
 * Structure storing pre-determined quantities
 * complying to the conventions of the PhenomHM model.
 * convensions such as m1>=m2
 */
typedef struct tagPhenomHMStorage
{
    REAL8 m1;    /**< mass of larger body in solar masses */
    REAL8 m2;    /**< mass of lighter body in solar masses */
    REAL8 m1_SI; /**< mass of larger body in kg */
    REAL8 m2_SI; /**< mass of lighter body in kg */
    REAL8 Mtot;  /**< total mass in solar masses */
    REAL8 eta;   /**< symmetric mass-ratio */
    REAL8 chi1x; /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 chi1y; /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 chi1z; /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 chi2x; /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    REAL8 chi2y; /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    REAL8 chi2z; /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    REAL8Sequence *freqs;
    REAL8 deltaF;
    REAL8 f_min;
    REAL8 f_max;
    REAL8 f_ref;
    REAL8 Mf_ref; /**< reference frequnecy in geometric units */
    REAL8 phiRef;
    UINT4 freq_is_uniform; /**< If = 1 then assume uniform spaced, If = 0 then assume arbitrarily spaced. */
    size_t npts;           /**< number of points in waveform array */
    size_t ind_min;        /**< start index containing non-zero values */
    size_t ind_max;        /**< end index containing non-zero values */
    REAL8 finmass;
    REAL8 finspin;
    REAL8 Mf_RD_22;
    REAL8 Mf_DM_22;
    REAL8 PhenomHMfring[L_MAX_PLUS_1][L_MAX_PLUS_1];
    REAL8 PhenomHMfdamp[L_MAX_PLUS_1][L_MAX_PLUS_1];
    REAL8 Rholm[L_MAX_PLUS_1][L_MAX_PLUS_1]; /**< ratio of (2,2) mode to (l,m) mode ringdown frequency */
    REAL8 Taulm[L_MAX_PLUS_1][L_MAX_PLUS_1]; /**< ratio of (l,m) mode to (2,2) mode damping time */
} PhenomHMStorage;

static int init_PhenomHM_Storage(
    PhenomHMStorage *p, /**< [out] PhenomHMStorage struct */
    const REAL8 m1_SI,  /**< mass of companion 1 (kg) */
    const REAL8 m2_SI,  /**< mass of companion 2 (kg) */
    const REAL8 chi1x,  /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    const REAL8 chi1y,  /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    const REAL8 chi1z,  /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    const REAL8 chi2x,  /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    const REAL8 chi2y,  /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    const REAL8 chi2z,  /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    REAL8Sequence *freqs, /**< Input frequency sequence (Hz) */
    const REAL8 deltaF,   /**< frequency spacing (Hz) */
    const REAL8 f_ref,    /**< reference GW frequency (hz) */
    const REAL8 phiRef    /**< orbital phase at f_ref */
);

double IMRPhenomHMTrd(
    REAL8 Mf,
    REAL8 Mf_RD_22,
    REAL8 Mf_RD_lm,
    const INT4 AmpFlag,
    const INT4 ell,
    const INT4 mm,
    PhenomHMStorage *pHM);

double IMRPhenomHMTi(
    REAL8 Mf,
    const INT4 mm);

int IMRPhenomHMSlopeAmAndBm(
    double *Am,
    double *Bm,
    const INT4 mm,
    REAL8 fi,
    REAL8 fr,
    REAL8 Mf_RD_22,
    REAL8 Mf_RD_lm,
    const INT4 AmpFlag,
    const INT4 ell,
    PhenomHMStorage *pHM);

int IMRPhenomHMMapParams(
    REAL8 *a,
    REAL8 *b,
    REAL8 flm,
    REAL8 fi,
    REAL8 fr,
    REAL8 Ai,
    REAL8 Bi,
    REAL8 Am,
    REAL8 Bm,
    REAL8 Ar,
    REAL8 Br);

int IMRPhenomHMFreqDomainMapParams(
    REAL8 *a,
    REAL8 *b,
    REAL8 *fi,
    REAL8 *fr,
    REAL8 *f1,
    const REAL8 flm,
    const INT4 ell,
    const INT4 mm,
    PhenomHMStorage *pHM,
    const int AmpFlag);

double IMRPhenomHMFreqDomainMap(
    REAL8 Mflm,
    const INT4 ell,
    const INT4 mm,
    PhenomHMStorage *pHM,
    const int AmpFlag);

int IMRPhenomHMPhasePreComp(
    HMPhasePreComp *q,
    const INT4 ell,
    const INT4 mm,
    PhenomHMStorage *pHM,
    LALDict *extraParams);

COMPLEX16 IMRPhenomHMOnePointFiveSpinPN(
    REAL8 fM,
    INT4 l,
    INT4 m,
    REAL8 M1,
    REAL8 M2,
    REAL8 X1z,
    REAL8 X2z);

int IMRPhenomHMCore(
    COMPLEX16FrequencySeries **hptilde,
    COMPLEX16FrequencySeries **hctilde,
    REAL8Sequence *freqs,
    REAL8 m1_SI,
    REAL8 m2_SI,
    REAL8 chi1z,
    REAL8 chi2z,
    const REAL8 distance,
    const REAL8 inclination,
    const REAL8 phiRef,
    const REAL8 deltaF,
    REAL8 f_ref,
    LALDict *extraParams);

int IMRPhenomHMEvaluateOnehlmMode(
    COMPLEX16FrequencySeries **hlm,
    REAL8Sequence *amps,
    REAL8Sequence *phases,
    REAL8Sequence *freqs_geom,
    PhenomHMStorage *pHM,
    UINT4 ell,
    INT4 mm,
    REAL8 phi0,
    LALDict *extraParams);

int IMRPhenomHMAmplitude(
    REAL8Sequence *amps,
    REAL8Sequence *freqs_geom,
    PhenomHMStorage *pHM,
    UINT4 ell,
    INT4 mm,
    LALDict *extraParams);

int IMRPhenomHMPhase(
    REAL8Sequence *phases,
    REAL8Sequence *freqs_geom,
    PhenomHMStorage *pHM,
    UINT4 ell,
    INT4 mm,
    LALDict *extraParams);

int IMRPhenomHMGetRingdownFrequency(
    REAL8 *fringdown,
    REAL8 *fdamp,
    UINT4 ell,
    INT4 mm,
    REAL8 finalmass,
    REAL8 finalspin);

#ifdef __cplusplus
}
#endif

#endif /* _LALSIM_IMR_PHENOMHM_H */