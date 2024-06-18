#include <lal/LALStdlib.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALDict.h>
#include "LALSimInspiralGenerator_private.h"
#include <lal/LALSimIMR.h>
#include "check_series_macros.h"
#include "check_waveform_macros.h"
#include "LALSimUniversalRelations.h"
#include "rotation_macros.h"
#include "fix_reference_frequency_macro.h"

#include <lal/FrequencySeries.h>
#include <lal/AVFactories.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif


static int XLALSimInspiralChooseTDWaveform_legacy(
    REAL8TimeSeries **hplus,
    REAL8TimeSeries **hcross,
    REAL8 m1,
    REAL8 m2,
    REAL8 S1x,
    REAL8 S1y,
    REAL8 S1z,
    REAL8 S2x,
    REAL8 S2y,
    REAL8 S2z,
    REAL8 distance,
    REAL8 inclination,
    REAL8 phiRef,
    REAL8 longAscNodes,
    REAL8 eccentricity,
    REAL8 meanPerAno,
    REAL8 deltaT,
    REAL8 f_min,
    REAL8 f_ref,
    LALDict *params,
    Approximant approximant
);

static int XLALSimInspiralChooseFDWaveform_legacy(
    COMPLEX16FrequencySeries ** hptilde,
    COMPLEX16FrequencySeries ** hctilde,
    REAL8 m1,
    REAL8 m2,
    REAL8 S1x,
    REAL8 S1y,
    REAL8 S1z,
    REAL8 S2x,
    REAL8 S2y,
    REAL8 S2z,
    REAL8 distance,
    REAL8 inclination,
    REAL8 phiRef,
    REAL8 longAscNodes,
    REAL8 eccentricity,
    REAL8 meanPerAno,
    REAL8 deltaF,
    REAL8 f_min,
    REAL8 f_max,
    REAL8 f_ref,
    LALDict * params,
    Approximant approximant
);

static SphHarmTimeSeries *XLALSimInspiralChooseTDModes_legacy(
    REAL8 phiRef,
    REAL8 deltaT,
    REAL8 m1,
    REAL8 m2,
    REAL8 S1x,
    REAL8 S1y,
    REAL8 S1z,
    REAL8 S2x,
    REAL8 S2y,
    REAL8 S2z,
    REAL8 f_min,
    REAL8 f_ref,
    REAL8 r,
    LALDict * LALpars,
    int lmax,
    Approximant approximant
);

static SphHarmFrequencySeries *XLALSimInspiralChooseFDModes_legacy(
    REAL8 m1,
    REAL8 m2,
    REAL8 S1x,
    REAL8 S1y,
    REAL8 S1z,
    REAL8 S2x,
    REAL8 S2y,
    REAL8 S2z,
    REAL8 deltaF,
    REAL8 f_min,
    REAL8 f_max,
    REAL8 f_ref,
    REAL8 phiRef,
    REAL8 distance,
    REAL8 inclination,
    LALDict * params,
    Approximant approximant
);

/** 
 * Method to initialize generator. This activates the conditioning for proper Fourier transform in the option `condition` is on in the LALDict
 */

static int initialize(LALSimInspiralGenerator * myself, LALDict *params)
{
    if (params && XLALDictContains(params, "condition")) {
        int condition;
        int errnum;
        XLAL_TRY(condition = XLALDictLookupINT4Value(params, "condition"), errnum);
        if (errnum)
            XLAL_ERROR(errnum);
        if (condition) {
            Approximant approximant = *(Approximant *)myself->internal_data;
            return XLALSimInspiralGeneratorAddConditioningForApproximant(myself, approximant);
        }
    }
    return 0;
}

/**
 * Define waveform generator methods to generate polarizations or modes in time or Fourier domain for legacy approximants
 */
 
/** Fourier domain modes */
static int generate_fd_modes(
    SphHarmFrequencySeries **hlm,
    LALDict *params,
    LALSimInspiralGenerator *myself
)
{
    REAL8 m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, deltaF, f_min, f_max, f_ref, phiRef, distance, inclination;
    Approximant approximant;

    /* approximant for this generator */
    approximant = *(Approximant *)myself->internal_data;

    XLALSimInspiralParseDictionaryToChooseFDModes(&m1, &m2, &S1x, &S1y, &S1z, &S2x, &S2y, &S2z, &deltaF, &f_min, &f_max, &f_ref, &phiRef, &distance, &inclination, params);

    *hlm = XLALSimInspiralChooseFDModes_legacy(m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, deltaF, f_min, f_max, f_ref, phiRef, distance, inclination, params, approximant);

    return 0;
}

/** Fourier domain polarizations */
static int generate_fd_waveform(
    COMPLEX16FrequencySeries **hplus,
    COMPLEX16FrequencySeries **hcross,
    LALDict *params,
    LALSimInspiralGenerator *myself
)
{
    REAL8 m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, deltaF, f_min, f_max, f_ref; 
    Approximant approximant;

    /* approximant for this generator */
    approximant = *(Approximant *)myself->internal_data;

    XLALSimInspiralParseDictionaryToChooseFDWaveform(&m1, &m2, &S1x, &S1y, &S1z, &S2x, &S2y, &S2z, &distance, &inclination, &phiRef, &longAscNodes, &eccentricity, &meanPerAno, &deltaF, &f_min, &f_max, &f_ref, params);

    return XLALSimInspiralChooseFDWaveform_legacy(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, deltaF, f_min, f_max, f_ref, params, approximant);
}

/** Time domain modes */
static int generate_td_modes(
    SphHarmTimeSeries **hlm,
    LALDict *params,
    LALSimInspiralGenerator *myself
)
{
    REAL8 phiRef, deltaT, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref, distance;
    INT4 lmax;
    Approximant approximant;

    /* approximant for this generator */
    approximant = *(Approximant *)myself->internal_data;

    XLALSimInspiralParseDictionaryToChooseTDModes(&phiRef, &deltaT, &m1, &m2, &S1x, &S1y, &S1z, &S2x, &S2y, &S2z, &f_min, &f_ref, &distance, &lmax, params);

    *hlm = XLALSimInspiralChooseTDModes_legacy(phiRef, deltaT, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref, distance, params, lmax, approximant);

    return 0;
}

/** Time domain polarizations */
static int generate_td_waveform(
    REAL8TimeSeries **hplus,
    REAL8TimeSeries **hcross,
    LALDict *params,
    LALSimInspiralGenerator *myself
)
{
    REAL8 m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, deltaT, f_min, f_ref;
    Approximant approximant;

    /* approximant for this generator */
    approximant = *(Approximant *)myself->internal_data;

    /* read parameters needed for ChooseTDWaveform from lal dictionary */
    XLALSimInspiralParseDictionaryToChooseTDWaveform(&m1, &m2, &S1x, &S1y, &S1z, &S2x, &S2y, &S2z, &distance, &inclination, &phiRef, &longAscNodes, &eccentricity, &meanPerAno, &deltaT, &f_min, &f_ref, params);

    return XLALSimInspiralChooseTDWaveform_legacy(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, deltaT, f_min, f_ref, params, approximant);
}

/** 
 * Define which methods are supported by every legacy approximant
 */

#define DEFINE_GENERATOR_TEMPLATE(approx, fd_modes, fd_waveform, td_modes, td_waveform) \
    static Approximant _lal ## approx ## GeneratorInternalData = approx; \
    const LALSimInspiralGenerator lal ## approx ## GeneratorTemplate = { \
        .name = #approx,  \
        .initialize = initialize, \
        .finalize = NULL, \
        .generate_fd_modes = fd_modes, \
        .generate_fd_waveform = fd_waveform, \
        .generate_td_modes = td_modes, \
        .generate_td_waveform = td_waveform, \
        .internal_data = &_lal ## approx ## GeneratorInternalData \
    };

/* TD POLARIZATIONS ONLY */
DEFINE_GENERATOR_TEMPLATE(EccentricTD, NULL, NULL, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(HGimri, NULL, NULL, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(IMRPhenomT, NULL, NULL, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(IMRPhenomTP, NULL, NULL, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(NR_hdf5, NULL, NULL, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(PhenSpinTaylor, NULL, NULL, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(PhenSpinTaylorRD, NULL, NULL, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(SEOBNRv1, NULL, NULL, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(SEOBNRv2, NULL, NULL, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(SEOBNRv2T, NULL, NULL, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(SEOBNRv2_opt, NULL, NULL, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(SEOBNRv3, NULL, NULL, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(SEOBNRv3_opt, NULL, NULL, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(SEOBNRv3_opt_rk4, NULL, NULL, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(SEOBNRv3_pert, NULL, NULL, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(SEOBNRv4, NULL, NULL, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(SEOBNRv4HM, NULL, NULL, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(SEOBNRv4T, NULL, NULL, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(SEOBNRv4_opt, NULL, NULL, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(SpinDominatedWf, NULL, NULL, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(TEOBResumS, NULL, NULL, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(TEOBResum_ROM, NULL, NULL, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(TaylorEt, NULL, NULL, NULL, generate_td_waveform)

/* TD POLARIZATIONS AND MODES ONLY */
DEFINE_GENERATOR_TEMPLATE(EOBNRv2, NULL, NULL, generate_td_modes, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(EOBNRv2HM, NULL, NULL, generate_td_modes, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(IMRPhenomTHM, NULL, NULL, generate_td_modes, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(IMRPhenomTPHM, NULL, NULL, generate_td_modes, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(NRHybSur3dq8, NULL, NULL, generate_td_modes, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(NRSur7dq2, NULL, NULL, generate_td_modes, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(NRSur7dq4, NULL, NULL, generate_td_modes, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(SEOBNRv4P, NULL, NULL, generate_td_modes, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(SEOBNRv4PHM, NULL, NULL, generate_td_modes, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(SEOBNRv4HM_PA, NULL, NULL, generate_td_modes, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(pSEOBNRv4HM_PA, NULL, NULL, generate_td_modes, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(SpinTaylorT1, NULL, NULL, generate_td_modes, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(SpinTaylorT4, NULL, NULL, generate_td_modes, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(SpinTaylorT5, NULL, NULL, generate_td_modes, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(TaylorT1, NULL, NULL, generate_td_modes, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(TaylorT2, NULL, NULL, generate_td_modes, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(TaylorT3, NULL, NULL, generate_td_modes, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(TaylorT4, NULL, NULL, generate_td_modes, generate_td_waveform)

/* FD POLARIZATIONS ONLY */
DEFINE_GENERATOR_TEMPLATE(EOBNRv2HM_ROM, NULL, generate_fd_waveform, NULL, NULL)
DEFINE_GENERATOR_TEMPLATE(EOBNRv2_ROM, NULL, generate_fd_waveform, NULL, NULL)
DEFINE_GENERATOR_TEMPLATE(EccentricFD, NULL, generate_fd_waveform, NULL, NULL)
DEFINE_GENERATOR_TEMPLATE(IMRPhenomD_NRTidal, NULL, generate_fd_waveform, NULL, NULL)
DEFINE_GENERATOR_TEMPLATE(IMRPhenomP, NULL, generate_fd_waveform, NULL, NULL)
DEFINE_GENERATOR_TEMPLATE(Lackey_Tidal_2013_SEOBNRv2_ROM, NULL, generate_fd_waveform, NULL, NULL)
DEFINE_GENERATOR_TEMPLATE(NRSur4d2s, NULL, generate_fd_waveform, NULL, NULL)
DEFINE_GENERATOR_TEMPLATE(SEOBNRv1_ROM_DoubleSpin, NULL, generate_fd_waveform, NULL, NULL)
DEFINE_GENERATOR_TEMPLATE(SEOBNRv1_ROM_EffectiveSpin, NULL, generate_fd_waveform, NULL, NULL)
DEFINE_GENERATOR_TEMPLATE(SEOBNRv2_ROM_DoubleSpin, NULL, generate_fd_waveform, NULL, NULL)
DEFINE_GENERATOR_TEMPLATE(SEOBNRv2_ROM_DoubleSpin_HI, NULL, generate_fd_waveform, NULL, NULL)
DEFINE_GENERATOR_TEMPLATE(SEOBNRv2_ROM_EffectiveSpin, NULL, generate_fd_waveform, NULL, NULL)
DEFINE_GENERATOR_TEMPLATE(SEOBNRv4T_surrogate, NULL, generate_fd_waveform, NULL, NULL)
DEFINE_GENERATOR_TEMPLATE(SEOBNRv4_ROM, NULL, generate_fd_waveform, NULL, NULL)
DEFINE_GENERATOR_TEMPLATE(SEOBNRv4_ROM_NRTidal, NULL, generate_fd_waveform, NULL, NULL)
DEFINE_GENERATOR_TEMPLATE(SpinTaylorF2, NULL, generate_fd_waveform, NULL, NULL)
DEFINE_GENERATOR_TEMPLATE(SpinTaylorT4Fourier, NULL, generate_fd_waveform, NULL, NULL)
DEFINE_GENERATOR_TEMPLATE(SpinTaylorT5Fourier, NULL, generate_fd_waveform, NULL, NULL)
DEFINE_GENERATOR_TEMPLATE(TaylorF2, NULL, generate_fd_waveform, NULL, NULL)
DEFINE_GENERATOR_TEMPLATE(TaylorF2Ecc, NULL, generate_fd_waveform, NULL, NULL)
DEFINE_GENERATOR_TEMPLATE(TaylorF2NLTides, NULL, generate_fd_waveform, NULL, NULL)
DEFINE_GENERATOR_TEMPLATE(TaylorF2RedSpin, NULL, generate_fd_waveform, NULL, NULL)
DEFINE_GENERATOR_TEMPLATE(TaylorF2RedSpinTidal, NULL, generate_fd_waveform, NULL, NULL)
DEFINE_GENERATOR_TEMPLATE(TaylorR2F4, NULL, generate_fd_waveform, NULL, NULL)

/* FD POLARIZATIONS AND MODES ONLY */
DEFINE_GENERATOR_TEMPLATE(SEOBNRv4HM_ROM, generate_fd_modes, generate_fd_waveform, NULL, NULL)
DEFINE_GENERATOR_TEMPLATE(SEOBNRv5_ROM, generate_fd_modes, generate_fd_waveform, NULL, NULL)

/* TD AND FD POLARIZATIONS ONLY */
DEFINE_GENERATOR_TEMPLATE(IMRPhenomA, NULL, generate_fd_waveform, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(IMRPhenomB, NULL, generate_fd_waveform, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(IMRPhenomC, NULL, generate_fd_waveform, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(IMRPhenomD, NULL, generate_fd_waveform, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(IMRPhenomD_NRTidalv2, NULL, generate_fd_waveform, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(IMRPhenomNSBH, NULL, generate_fd_waveform, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(IMRPhenomPv2, NULL, generate_fd_waveform, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(IMRPhenomPv2_NRTidal, NULL, generate_fd_waveform, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(IMRPhenomPv2_NRTidalv2, NULL, generate_fd_waveform, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(IMRPhenomPv3, NULL, generate_fd_waveform, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(IMRPhenomPv3HM, NULL, generate_fd_waveform, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(IMRPhenomXAS, NULL, generate_fd_waveform, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(IMRPhenomXP, NULL, generate_fd_waveform, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(IMRPhenomXAS_NRTidalv2, NULL, generate_fd_waveform, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(IMRPhenomXP_NRTidalv2, NULL, generate_fd_waveform, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(IMRPhenomXAS_NRTidalv3, NULL, generate_fd_waveform, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(IMRPhenomXP_NRTidalv3, NULL, generate_fd_waveform, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(SEOBNRv4_ROM_NRTidalv2, NULL, generate_fd_waveform, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(SEOBNRv4_ROM_NRTidalv2_NSBH, NULL, generate_fd_waveform, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(SEOBNRv5_ROM_NRTidalv3, NULL, generate_fd_waveform, NULL, generate_td_waveform)

/* TD POLARIZATIONS AND FD POLARIZATIONS AND MODES ONLY */
DEFINE_GENERATOR_TEMPLATE(IMRPhenomHM, generate_fd_modes, generate_fd_waveform, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(IMRPhenomXHM, generate_fd_modes, generate_fd_waveform, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(IMRPhenomXPHM, generate_fd_modes, generate_fd_waveform, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(IMRPhenomXO4a, generate_fd_modes, generate_fd_waveform, NULL, generate_td_waveform)
DEFINE_GENERATOR_TEMPLATE(SEOBNRv5HM_ROM, generate_fd_modes, generate_fd_waveform, NULL, generate_td_waveform)


/**
 * Copy of the old code of XLALSimInspiralChooseTDWaveform(). The new version of XLALSimInspiralChooseTDWaveform() is just a wrapper over XLALSimInspiralGenerateTDWaveform().
 * XLALSimInspiralGenerateTDWaveform() internally calls this function for legacy approximants.
 */
static int XLALSimInspiralChooseTDWaveform_legacy(
  REAL8TimeSeries **hplus,              /* +-polarization waveform */
  REAL8TimeSeries **hcross,             /* x-polarization waveform */
  REAL8 m1,                             /* mass of companion 1 (kg) */
  REAL8 m2,                             /* mass of companion 2 (kg) */
  REAL8 S1x,                            /* x-component of the dimensionless spin of object 1 */
  REAL8 S1y,                            /* y-component of the dimensionless spin of object 1 */
  REAL8 S1z,                            /* z-component of the dimensionless spin of object 1 */
  REAL8 S2x,                            /* x-component of the dimensionless spin of object 2 */
  REAL8 S2y,                            /* y-component of the dimensionless spin of object 2 */
  REAL8 S2z,                            /* z-component of the dimensionless spin of object 2 */
  REAL8 distance,                       /* distance of source (m) */
  REAL8 inclination,                    /* inclination of source (rad) */
  REAL8 phiRef,                         /* reference orbital phase (rad) */
  REAL8 longAscNodes,                   /* longitude of ascending nodes, degenerate with the polarization angle, Omega in documentation */
  REAL8 eccentricity,                   /* eccentrocity at reference epoch */
  REAL8 meanPerAno,                     /* mean anomaly of periastron */
  REAL8 deltaT,                         /* sampling interval (s) */
  REAL8 f_min,                          /* starting GW frequency (Hz) */
  REAL8 f_ref,                          /* reference GW frequency (Hz) */
  LALDict *params,                      /* LAL dictionary containing accessory parameters */
  Approximant approximant               /* post-Newtonian approximant to use for waveform production */

)
{
    REAL8 LNhatx, LNhaty, LNhatz, E1x, E1y, E1z;
    //REAL8 tmp1, tmp2;
    int ret;
    /* N.B. the quadrupole of a spinning compact body labeled by A is
     * Q_A = - quadparam_A chi_A^2 m_A^3 (see gr-qc/9709032)
     * where quadparam = 1 for BH ~= 4-8 for NS.
     * This affects the quadrupole-monopole interaction.
     */
    REAL8 v0 = 1.;
    /* Note: approximant SEOBNRv2T/v4T will by default compute dQuadMon1, dQuadMon2 */
    /* from TidalLambda1, TidalLambda2 using universal relations, */
    /* or use the input value if it is present in the dictionary params */
    REAL8 quadparam1 = 1. + XLALSimInspiralWaveformParamsLookupdQuadMon1(params);
    REAL8 quadparam2 = 1. + XLALSimInspiralWaveformParamsLookupdQuadMon2(params);
    REAL8 lambda1 = XLALSimInspiralWaveformParamsLookupTidalLambda1(params);
    REAL8 lambda2 = XLALSimInspiralWaveformParamsLookupTidalLambda2(params);
    int amplitudeO = XLALSimInspiralWaveformParamsLookupPNAmplitudeOrder(params);
    int phaseO = XLALSimInspiralWaveformParamsLookupPNPhaseOrder(params);
    /* Tidal parameters to be computed, if required, by universal relations */
    REAL8 lambda3A_UR = 0.;
    REAL8 omega2TidalA_UR = 0.;
    REAL8 omega3TidalA_UR = 0.;
    REAL8 lambda3B_UR = 0.;
    REAL8 omega2TidalB_UR = 0.;
    REAL8 omega3TidalB_UR = 0.;
    REAL8 quadparam1_UR = 0.;
    REAL8 quadparam2_UR = 0.;

    /* General sanity checks that will abort
     *
     * If non-GR approximants are added, include them in
     * XLALSimInspiralApproximantAcceptTestGRParams()
     */
    if (!XLALSimInspiralWaveformParamsNonGRAreDefault(params) && XLALSimInspiralApproximantAcceptTestGRParams(approximant) != LAL_SIM_INSPIRAL_TESTGR_PARAMS) {
        XLALPrintError("XLAL Error - %s: Passed in non-NULL pointer to LALSimInspiralTestGRParam for an approximant that does not use LALSimInspiralTestGRParam\n", __func__);
        XLAL_ERROR(XLAL_EINVAL);
    }
    /* Support variables for precessing wfs */
    REAL8 incl;


    /* SEOBNR flag for spin aligned model version. 1 for SEOBNRv1, 2 for SEOBNRv2 */
    UINT4 SpinAlignedEOBversion;
    REAL8 spin1x, spin1y, spin1z;
    REAL8 spin2x, spin2y, spin2z;
    REAL8 polariz = longAscNodes;

    /* SEOBNR flag for precessing model version. 3 for SEOBNRv3, 300 for SEOBNRv3_opt, 401 for SEOBNRv4P, 402 for SEOBNRv4PHM */
    UINT4 PrecEOBversion;
    REAL8 spin1[3], spin2[3];

    REAL8 maxamp = 0;
    INT4 loopi = 0;
    INT4 maxind = 0;

    //LIGOTimeGPS epoch = LIGOTIMEGPSZERO;

    /* General sanity check the input parameters - only give warnings! */
    if (deltaT > 1.)
        XLALPrintWarning("XLAL Warning - %s: Large value of deltaT = %e requested.\nPerhaps sample rate and time step size were swapped?\n", __func__, deltaT);
    if (deltaT < 1. / 16385.)
        XLALPrintWarning("XLAL Warning - %s: Small value of deltaT = %e requested.\nCheck for errors, this could create very large time series.\n", __func__, deltaT);
    if (m1 < 0.09 * LAL_MSUN_SI)
        XLALPrintWarning("XLAL Warning - %s: Small value of m1 = %e (kg) = %e (Msun) requested.\nPerhaps you have a unit conversion error?\n", __func__, m1, m1 / LAL_MSUN_SI);
    if (m2 < 0.09 * LAL_MSUN_SI)
        XLALPrintWarning("XLAL Warning - %s: Small value of m2 = %e (kg) = %e (Msun) requested.\nPerhaps you have a unit conversion error?\n", __func__, m2, m2 / LAL_MSUN_SI);
    if (m1 + m2 > 1000. * LAL_MSUN_SI)
        XLALPrintWarning("XLAL Warning - %s: Large value of total mass m1+m2 = %e (kg) = %e (Msun) requested.\nSignal not likely to be in band of ground-based detectors.\n", __func__, m1 + m2, (m1 + m2) / LAL_MSUN_SI);
    if (S1x * S1x + S1y * S1y + S1z * S1z > 1.000001)
        XLALPrintWarning("XLAL Warning - %s: S1 = (%e,%e,%e) with norm > 1 requested.\nAre you sure you want to violate the Kerr bound?\n", __func__, S1x, S1y, S1z);
    if (S2x * S2x + S2y * S2y + S2z * S2z > 1.000001)
        XLALPrintWarning("XLAL Warning - %s: S2 = (%e,%e,%e) with norm > 1 requested.\nAre you sure you want to violate the Kerr bound?\n", __func__, S2x, S2y, S2z);
    if (f_min < 1.)
        XLALPrintWarning("XLAL Warning - %s: Small value of fmin = %e requested.\nCheck for errors, this could create a very long waveform.\n", __func__, f_min);
    if (f_min > 40.000001)
        XLALPrintWarning("XLAL Warning - %s: Large value of fmin = %e requested.\nCheck for errors, the signal will start in band.\n", __func__, f_min);

    /* adjust the reference frequency for certain precessing approximants:
     * if that approximate interprets f_ref==0 to be f_min, set f_ref=f_min;
     * otherwise do nothing */
    FIX_REFERENCE_FREQUENCY(f_ref, f_min, approximant);

    switch (approximant) {
        /* non-spinning inspiral-only models */
    case TaylorEt:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero spins were given, but this is a non-spinning approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        if (f_ref != 0.)
            XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);
        /* Call the waveform driver routine */
        ret = XLALSimInspiralTaylorEtPNGenerator(hplus, hcross, phiRef, v0, deltaT, m1, m2, f_min, distance, inclination, amplitudeO, phaseO);
        break;

    case TaylorT1:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFrameAxisIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
        if (!XLALSimInspiralWaveformParamsModesChoiceIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        if (!XLALSimInspiralWaveformParamsPNSpinOrderIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralSpinOrder provided, but this approximant does not use that flag.");
        if (!checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero spins were given, but this is a non-spinning approximant.");
        /* Call the waveform driver routine */
        ret = XLALSimInspiralTaylorT1PNGenerator(hplus, hcross, phiRef, v0, deltaT, m1, m2, f_min, f_ref, distance, inclination, lambda1, lambda2, XLALSimInspiralWaveformParamsLookupPNTidalOrder(params), amplitudeO, phaseO);
        break;

    case TaylorT2:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFrameAxisIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
        if (!XLALSimInspiralWaveformParamsModesChoiceIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        if (!XLALSimInspiralWaveformParamsPNSpinOrderIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralSpinOrder provided, but this approximant does not use that flag.");
        if (!checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero spins were given, but this is a non-spinning approximant.");
        /* Call the waveform driver routine */
        ret = XLALSimInspiralTaylorT2PNGenerator(hplus, hcross, phiRef, v0, deltaT, m1, m2, f_min, f_ref, distance, inclination, lambda1, lambda2, XLALSimInspiralWaveformParamsLookupPNTidalOrder(params), amplitudeO, phaseO);
        break;

    case TaylorT3:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFrameAxisIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
        if (!XLALSimInspiralWaveformParamsModesChoiceIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        if (!XLALSimInspiralWaveformParamsPNSpinOrderIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralSpinOrder provided, but this approximant does not use that flag.");
        if (!checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero spins were given, but this is a non-spinning approximant.");
        /* Call the waveform driver routine */
        ret = XLALSimInspiralTaylorT3PNGenerator(hplus, hcross, phiRef, v0, deltaT, m1, m2, f_min, f_ref, distance, inclination, lambda1, lambda2, XLALSimInspiralWaveformParamsLookupPNTidalOrder(params), amplitudeO, phaseO);
        break;

    case TaylorT4:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFrameAxisIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
        if (!XLALSimInspiralWaveformParamsModesChoiceIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        if (!XLALSimInspiralWaveformParamsPNSpinOrderIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralSpinOrder provided, but this approximant does not use that flag.");
        if (!checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero spins were given, but this is a non-spinning approximant.");
        /* Call the waveform driver routine */
        ret = XLALSimInspiralTaylorT4PNGenerator(hplus, hcross, phiRef, v0, deltaT, m1, m2, f_min, f_ref, distance, inclination, lambda1, lambda2, XLALSimInspiralWaveformParamsLookupPNTidalOrder(params), amplitudeO, phaseO);
        break;

    case TEOBResum_ROM:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFrameAxisIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
        if (!XLALSimInspiralWaveformParamsModesChoiceIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        if (!XLALSimInspiralWaveformParamsPNSpinOrderIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralSpinOrder provided, but this approximant does not use that flag.");
        if (!checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero spins were given, but this is a non-spinning approximant.");
        /* Call the waveform driver routine */
        ret = XLALSimInspiralTEOBResumROM(hplus, hcross, phiRef, deltaT, f_min, f_ref, distance, inclination, m1, m2, lambda1, lambda2);
        break;

    case TEOBResumS:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFrameAxisIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
        if (!XLALSimInspiralWaveformParamsPNSpinOrderIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralSpinOrder provided, but this approximant does not use that flag.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (f_ref != 0.)
            XLALPrintWarning("XLAL Warning - %s: This approximant does not use f_ref. The reference phase will be defined at coalescence.\n", __func__);
        /* Comply with master convention on orientation angles */
        polariz += LAL_PI_2;
        /* Call the waveform driver routine */
        //  MA: TODO if( something about modes choice )

        /* Make sure params exists (otherwise segfault) */
        if (!params)
            params = XLALCreateDict();
        ret = XLALSimIMRTEOBResumS(hplus, hcross, phiRef, deltaT, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, lambda1, lambda2, distance, inclination, longAscNodes, params, eccentricity, meanPerAno, f_min, f_ref);
        break;


    case EccentricTD:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFrameAxisIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
        if (!XLALSimInspiralWaveformParamsModesChoiceIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        if (!XLALSimInspiralWaveformParamsPNSpinOrderIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralSpinOrder provided, but this approximant does not use that flag.");
        if (!checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero spins were given, but this is a non-spinning approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        /* Call the waveform driver routine */
        ret = XLALSimInspiralEccentricTDPNGenerator(hplus, hcross, phiRef, deltaT, m1, m2, f_min, f_ref, distance, inclination, eccentricity, amplitudeO, phaseO);
        if (ret == XLAL_FAILURE)
            XLAL_ERROR(XLAL_EFUNC);
        break;

        /* non-spinning inspiral-merger-ringdown models */
    case IMRPhenomA:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero spins were given, but this is a non-spinning approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        if (f_ref != 0.)
            XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);
        /* Call the waveform driver routine */
        // NB: f_max = 0 will generate up to the ringdown cut-off frequency
        ret = XLALSimIMRPhenomAGenerateTD(hplus, hcross, phiRef, deltaT, m1, m2, f_min, 0., distance, inclination);
        break;

    case EOBNRv2HM:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero spins were given, but this is a non-spinning approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        if (f_ref != 0.)
            XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);
        /* Call the waveform driver routine */
        // FIXME: need to create a function to take in different modes or produce an error if all modes not given
        ret = XLALSimIMREOBNRv2AllModes(hplus, hcross, phiRef, deltaT, m1, m2, f_min, distance, inclination);
        break;

    case EOBNRv2:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero spins were given, but this is a non-spinning approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        if (f_ref != 0.)
            XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);
        /* Call the waveform driver routine */
        ret = XLALSimIMREOBNRv2DominantMode(hplus, hcross, phiRef, deltaT, m1, m2, f_min, distance, inclination);
        break;

        /* spinning inspiral-only models */
    case SpinTaylorT5:
        /* Waveform-specific sanity checks */
        //if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) && (XLALSimInspiralWaveformParamsLookupPNSpinOrder(params)>5) )
        //XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given with spinOrder %d, but this spins dynamics is not implemented to this order for precessing spins.", XLALSimInspiralWaveformParamsLookupPNSpinOrder);
        XLALSimInspiralInitialConditionsPrecessingApproxs(&incl, &spin1x, &spin1y, &spin1z, &spin2x, &spin2y, &spin2z, inclination, S1x, S1y, S1z, S2x, S2y, S2z, m1, m2, f_ref, phiRef, XLALSimInspiralWaveformParamsLookupFrameAxis(params));
        LNhatx = sin(incl);
        LNhaty = 0.;
        LNhatz = cos(incl);
        E1x = 0.;
        E1y = 1.;
        E1z = 0.;
        polariz += LAL_PI / 2.;
        /* Call the waveform driver routine */
        ret = XLALSimInspiralSpinTaylorT5(hplus, hcross, phiRef, deltaT, m1, m2, f_min, f_ref, distance, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, LNhatx, LNhaty, LNhatz, E1x, E1y, E1z, params);
        break;

        // need to make a consistent choice for SpinTaylorT4 and PSpinInspiralRD waveform inputs
        // proposal: TotalJ frame of PSpinInspiralRD
        // inclination denotes the angle between the view direction
        // and J (J is constant during the evolution, J//z, both N and initial
        // L are in the x-z plane) and the spin coordinates are given wrt
        // initial ** L **.
    case SpinTaylorT4:
        /* Waveform-specific sanity checks */
        //if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) && (XLALSimInspiralWaveformParamsLookupPNSpinOrder(params)>5) )
        //  XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given with spinOrder %d, but this spins dynamics is not implemented to this order for precessing spins.", XLALSimInspiralWaveformParamsLookupPNSpinOrder);
        XLALSimInspiralInitialConditionsPrecessingApproxs(&incl, &spin1x, &spin1y, &spin1z, &spin2x, &spin2y, &spin2z, inclination, S1x, S1y, S1z, S2x, S2y, S2z, m1, m2, f_ref, phiRef, XLALSimInspiralWaveformParamsLookupFrameAxis(params));
        LNhatx = sin(incl);
        LNhaty = 0.;
        LNhatz = cos(incl);
        E1x = 0.;
        E1y = 1.;
        E1z = 0.;
        polariz += LAL_PI / 2.;
        /* Call the waveform driver routine */
        ret = XLALSimInspiralSpinTaylorT4(hplus, hcross, phiRef, deltaT, m1, m2, f_min, f_ref, distance, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, LNhatx, LNhaty, LNhatz, E1x, E1y, E1z, params);
        break;

    case SpinTaylorT1:
        /* Waveform-specific sanity checks */
        /* Waveform-specific sanity checks */
        //if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) && (XLALSimInspiralWaveformParamsLookupPNSpinOrder(params)>5) )
        //  XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given with spinOrder %d, but this spins dynamics is not implemented to this order for precessing spins.", XLALSimInspiralWaveformParamsLookupPNSpinOrder);
        XLALSimInspiralInitialConditionsPrecessingApproxs(&incl, &spin1x, &spin1y, &spin1z, &spin2x, &spin2y, &spin2z, inclination, S1x, S1y, S1z, S2x, S2y, S2z, m1, m2, f_ref, phiRef, XLALSimInspiralWaveformParamsLookupFrameAxis(params));
        LNhatx = sin(incl);
        LNhaty = 0.;
        LNhatz = cos(incl);
        E1x = 0.;
        E1y = 1.;
        E1z = 0.;
        polariz += LAL_PI / 2.;
        /* Call the waveform driver routine */
        ret = XLALSimInspiralSpinTaylorT1(hplus, hcross, phiRef, deltaT, m1, m2, f_min, f_ref, distance, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, LNhatx, LNhaty, LNhatz, E1x, E1y, E1z, params);
        break;

    case SpinDominatedWf:
        // waveform specific sanity checks
        if (S2x != 0. || S2y != 0. || S2z != 0.) {
            XLALPrintError("XLAL Error : The spindominatedwf approximant is only for 1 spin case.\n");
            XLAL_ERROR(XLAL_EDOM);
        }
        /*Maximal PN amplitude order is 1.5, maximal phase order is 2 PN */
        if (amplitudeO > 3) {
            XLALPrintError("XLAL Error : Foe the spindominatedwf approximant maximal amplitude correction is 1.5 PN\n");
            XLAL_ERROR(XLAL_EDOM);
        }
        if (phaseO > 4) {
            XLALPrintError("XLAL Error : For the spindominatedwf approximant maximal phase correction is 2 PN\n");
            XLAL_ERROR(XLAL_EDOM);
        }
        incl = inclination;
        LNhatx = 0.;
        LNhaty = 0.;
        LNhatz = 1.;
        /* Call the waveform driver routine */
        ret = XLALSimInspiralSpinDominatedWaveformInterfaceTD(hplus, hcross, deltaT, m1, m2, f_min, f_ref, distance, S1x, S1y, S1z, LNhatx, LNhaty, LNhatz, incl, phaseO, amplitudeO, phiRef);
        break;

        /* spin aligned inspiral-merger-ringdown models */
    case IMRPhenomB:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        if (f_ref != 0.)
            XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);
        /* Call the waveform driver routine */
        // NB: f_max = 0 will generate up to the ringdown cut-off frequency
        ret = XLALSimIMRPhenomBGenerateTD(hplus, hcross, phiRef, deltaT, m1, m2, XLALSimIMRPhenomBComputeChi(m1, m2, S1z, S2z), f_min, 0., distance, inclination);
        break;

    case PhenSpinTaylor:
        /* Waveform-specific sanity checks */

        /* Call the waveform driver routine */
        XLALSimInspiralInitialConditionsPrecessingApproxs(&incl, &spin1x, &spin1y, &spin1z, &spin2x, &spin2y, &spin2z, inclination, S1x, S1y, S1z, S2x, S2y, S2z, m1, m2, f_ref, phiRef, XLALSimInspiralWaveformParamsLookupFrameAxis(params));
        polariz += LAL_PI / 2.;

        ret = XLALSimSpinInspiralGenerator(hplus, hcross, phiRef, deltaT, m1, m2, f_min, f_ref, distance, incl, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, phaseO, amplitudeO, lambda1, lambda2, quadparam1, quadparam2, params);
        break;

    case IMRPhenomC:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        if (f_ref != 0.)
            XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);
        /* Call the waveform driver routine */
        // NB: f_max = 0 will generate up to the ringdown cut-off frequency
        ret = XLALSimIMRPhenomCGenerateTD(hplus, hcross, phiRef, deltaT, m1, m2, XLALSimIMRPhenomBComputeChi(m1, m2, S1z, S2z), f_min, 0., distance, inclination, params);
        break;

    case IMRPhenomD:
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        // generate TD waveforms with zero inclincation so that amplitude can be
        // calculated from hplus and hcross, apply inclination-dependent factors
        // in loop below
        ret = XLALSimInspiralTDFromFD(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, 0., phiRef, longAscNodes, eccentricity, meanPerAno, deltaT, f_min, f_ref, params, approximant);
        maxamp = 0;
        REAL8TimeSeries *hp = *hplus;
        REAL8TimeSeries *hc = *hcross;
        maxind = hp->data->length - 1;
        const REAL8 cfac = cos(inclination);
        const REAL8 pfac = 0.5 * (1. + cfac * cfac);
        for (loopi = hp->data->length - 1; loopi > -1; loopi--) {
            REAL8 ampsqr = (hp->data->data[loopi]) * (hp->data->data[loopi]) + (hc->data->data[loopi]) * (hc->data->data[loopi]);
            if (ampsqr > maxamp) {
                maxind = loopi;
                maxamp = ampsqr;
            }
            hp->data->data[loopi] *= pfac;
            hc->data->data[loopi] *= cfac;
        }
        XLALGPSSetREAL8(&(hp->epoch), (-1.) * deltaT * maxind);
        XLALGPSSetREAL8(&(hc->epoch), (-1.) * deltaT * maxind);
        break;

    case IMRPhenomHM:
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");

        ret = XLALSimInspiralTDFromFD(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, deltaT, f_min, f_ref, params, approximant);
        /*
         * NOTE: We enforce that the hp**2 + hx**2 peaks at t=0
         * see the wiki page from phenomHM review
         * https://git.ligo.org/waveforms/reviews/phenomhm/wikis/time-domain-behaviour
         */
        maxamp = 0;
        maxind = (*hplus)->data->length - 1;
        for (loopi = (*hplus)->data->length - 1; loopi > -1; loopi--) {
            REAL8 ampsqr = ((*hplus)->data->data[loopi]) * ((*hplus)->data->data[loopi]) + ((*hcross)->data->data[loopi]) * ((*hcross)->data->data[loopi]);
            if (ampsqr > maxamp) {
                maxind = loopi;
                maxamp = ampsqr;
            }
        }
        XLALGPSSetREAL8(&((*hplus)->epoch), (-1.) * deltaT * maxind);
        XLALGPSSetREAL8(&((*hcross)->epoch), (-1.) * deltaT * maxind);
        break;

    case IMRPhenomPv2:
        ret = XLALSimInspiralTDFromFD(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, deltaT, f_min, f_ref, params, approximant);
        break;

    case IMRPhenomPv3:
        ret = XLALSimInspiralTDFromFD(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, deltaT, f_min, f_ref, params, approximant);
        break;

    case IMRPhenomPv3HM:
        ret = XLALSimInspiralTDFromFD(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, deltaT, f_min, f_ref, params, approximant);
        break;

    case IMRPhenomPv2_NRTidal:
        ret = XLALSimInspiralTDFromFD(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, deltaT, f_min, f_ref, params, approximant);
        break;

    case IMRPhenomPv2_NRTidalv2:
        ret = XLALSimInspiralTDFromFD(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, deltaT, f_min, f_ref, params, approximant);
        break;

    case IMRPhenomD_NRTidalv2:
        ret = XLALSimInspiralTDFromFD(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, deltaT, f_min, f_ref, params, approximant);
        break;

    case SEOBNRv4_ROM_NRTidalv2:
        ret = XLALSimInspiralTDFromFD(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, deltaT, f_min, f_ref, params, approximant);
        break;

    case SEOBNRv5_ROM_NRTidalv3:
        ret = XLALSimInspiralTDFromFD(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, deltaT, f_min, f_ref, params, approximant);
        break;

    case IMRPhenomNSBH:
        ret = XLALSimInspiralTDFromFD(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, deltaT, f_min, f_ref, params, approximant);
        break;

    case SEOBNRv4_ROM_NRTidalv2_NSBH:
        ret = XLALSimInspiralTDFromFD(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, deltaT, f_min, f_ref, params, approximant);
        break;

    case SEOBNRv5HM_ROM:
        ret = XLALSimInspiralTDFromFD(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, deltaT, f_min, f_ref, params, approximant);
        break;

    case PhenSpinTaylorRD:
        /* Waveform-specific sanity checks */
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        if (f_ref != 0.)
            XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at the start.\n", __func__);
        /* Call the waveform driver routine */
        XLALSimInspiralInitialConditionsPrecessingApproxs(&incl, &spin1x, &spin1y, &spin1z, &spin2x, &spin2y, &spin2z, inclination, S1x, S1y, S1z, S2x, S2y, S2z, m1, m2, f_ref, phiRef, XLALSimInspiralWaveformParamsLookupFrameAxis(params));
        polariz += LAL_PI / 2.;
        ret = XLALSimIMRPhenSpinInspiralRDGenerator(hplus, hcross, phiRef, deltaT, m1, m2, f_min, f_ref, distance, incl, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, phaseO, amplitudeO, lambda1, lambda2, quadparam1, quadparam2, params);
        break;

    case SEOBNRv1:
    case SEOBNRv2_opt:
    case SEOBNRv2:
    case SEOBNRv4_opt:
    case SEOBNRv4:
    case SEOBNRv4HM:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        if (f_ref != 0.)
            XLALPrintWarning("XLAL Warning - %s: This approximant does not use f_ref. The reference phase will be defined at coalescence.\n", __func__);
        /* Call the waveform driver routine */
        polariz += -LAL_PI / 2.;
        //R.C. this rotation of -pi/2 is needed to go from the EOB wave frame to the LAL wave frame, see slide 9 of https://git.ligo.org/waveforms/reviews/SEOBNRv4HM/blob/master/tests/conventions/conventions.pdf
        if (approximant == SEOBNRv1)
            SpinAlignedEOBversion = 1;
        if (approximant == SEOBNRv2)
            SpinAlignedEOBversion = 2;
        if (approximant == SEOBNRv2_opt)
            SpinAlignedEOBversion = 200;
        if (approximant == SEOBNRv4)
            SpinAlignedEOBversion = 4;
        if (approximant == SEOBNRv4_opt)
            SpinAlignedEOBversion = 400;
        if (approximant == SEOBNRv4HM)
            SpinAlignedEOBversion = 41;
        ret = XLALSimIMRSpinAlignedEOBWaveform(hplus, hcross, phiRef, deltaT, m1, m2, f_min, distance, inclination, S1z, S2z, SpinAlignedEOBversion, params);
        break;

    case SEOBNRv4HM_PA:
    case pSEOBNRv4HM_PA:
        if( !XLALSimInspiralWaveformParamsFlagsAreDefault(params) )
          XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
          XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if( !checkTidesZero(lambda1, lambda2) )
          XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant does not have tidal corrections.");
        if( f_ref != 0.)
          XLALPrintWarning("XLAL Warning - %s: This approximant does not use f_ref. The reference phase will be defined at coalescence.\n", __func__);
        /* Call the waveform driver routine */
        polariz+=-LAL_PI/2.;
        //R.C. this rotation of -pi/2 is needed to go from the EOB wave frame to the LAL wave frame, see slide 9 of https://git.ligo.org/waveforms/reviews/SEOBNRv4HM/blob/master/tests/conventions/con \ventions.pdf                                                                                                                                                                                                
        if(approximant == SEOBNRv4HM_PA) SpinAlignedEOBversion = 4111;
        if(approximant == pSEOBNRv4HM_PA) SpinAlignedEOBversion = 4112;

        ret = XLALSimIMRSpinAlignedEOBWaveform(
          hplus, hcross,
          phiRef, deltaT,
          m1, m2,
          f_min,
          distance, inclination,
          S1z, S2z,
          SpinAlignedEOBversion, params
        );
        break;

    case SEOBNRv3_opt_rk4:
    case SEOBNRv3_opt:
    case SEOBNRv3_pert:
    case SEOBNRv3:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        if (f_ref != 0.)
            XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);
        /* Call the waveform driver routine */
        //REAL8 spin1[3], spin2[3];
        spin1[0] = S1x;
        spin1[1] = S1y;
        spin1[2] = S1z;
        spin2[0] = S2x;
        spin2[1] = S2y;
        spin2[2] = S2z;
        polariz += -LAL_PI / 2.;
        PrecEOBversion = 3;
        if (approximant == SEOBNRv3_pert) {
            const double m1pert = m1 * (1.0 + 1e-15);
            ret = XLALSimIMRSpinEOBWaveform(hplus, hcross, /*&epoch, */ phiRef,
                                            deltaT, m1pert, m2, f_min, distance, inclination, spin1, spin2, PrecEOBversion);

        } else {
            if (approximant == SEOBNRv3_opt)
                PrecEOBversion = 300;
            if (approximant == SEOBNRv3_opt_rk4)
                PrecEOBversion = 304;
            ret = XLALSimIMRSpinEOBWaveform(hplus, hcross, /*&epoch, */ phiRef,
                                            deltaT, m1, m2, f_min, distance, inclination, spin1, spin2, PrecEOBversion);
        }
        break;

    case SEOBNRv4P:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        if (f_ref != 0.)
            XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);

        spin1[0] = S1x;
        spin1[1] = S1y;
        spin1[2] = S1z;
        spin2[0] = S2x;
        spin2[1] = S2y;
        spin2[2] = S2z;
        polariz += -LAL_PI / 2.;
        PrecEOBversion = 401;
        ret = XLALSimIMRSpinPrecEOBWaveform(hplus, hcross, phiRef, deltaT, m1, m2, f_min, distance, inclination, spin1, spin2, PrecEOBversion, params);
        break;
    case SEOBNRv4PHM:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        if (f_ref != 0.)
            XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);

        spin1[0] = S1x;
        spin1[1] = S1y;
        spin1[2] = S1z;
        spin2[0] = S2x;
        spin2[1] = S2y;
        spin2[2] = S2z;
        polariz += -LAL_PI / 2.;
        PrecEOBversion = 402;
        ret = XLALSimIMRSpinPrecEOBWaveform(hplus, hcross, phiRef, deltaT, m1, m2, f_min, distance, inclination, spin1, spin2, PrecEOBversion, params);
        break;

    case SEOBNRv2T:
    case SEOBNRv4T:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (f_ref != 0.)
            XLALPrintWarning("XLAL Warning - %s: This approximant does not use f_ref. The reference phase will be defined at coalescence.\n", __func__);
        /* If tides-related parameter was not input by the user, use universal relations to compute it from quadrupolar lambda (or from octupolar lambda, itself either input or computed, for omega03) - else use the input value given by the user */
        if (!XLALDictContains(params, "TidalOctupolarLambda1")) {
            lambda3A_UR = XLALSimUniversalRelationlambda3TidalVSlambda2Tidal(lambda1);
            XLALSimInspiralWaveformParamsInsertTidalOctupolarLambda1(params, lambda3A_UR);
        }
        if (!XLALDictContains(params, "TidalOctupolarLambda2")) {
            lambda3B_UR = XLALSimUniversalRelationlambda3TidalVSlambda2Tidal(lambda2);
            XLALSimInspiralWaveformParamsInsertTidalOctupolarLambda2(params, lambda3B_UR);
        }
        if (!XLALDictContains(params, "TidalQuadrupolarFMode1")) {
            omega2TidalA_UR = XLALSimUniversalRelationomega02TidalVSlambda2Tidal(lambda1);
            XLALSimInspiralWaveformParamsInsertTidalQuadrupolarFMode1(params, omega2TidalA_UR);
        }
        if (!XLALDictContains(params, "TidalQuadrupolarFMode2")) {
            omega2TidalB_UR = XLALSimUniversalRelationomega02TidalVSlambda2Tidal(lambda2);
            XLALSimInspiralWaveformParamsInsertTidalQuadrupolarFMode2(params, omega2TidalB_UR);
        }
        if (!XLALDictContains(params, "TidalOctupolarFMode1")) {
            omega3TidalA_UR = XLALSimUniversalRelationomega03TidalVSlambda3Tidal(lambda3A_UR);
            XLALSimInspiralWaveformParamsInsertTidalOctupolarFMode1(params, omega3TidalA_UR);
        }
        if (!XLALDictContains(params, "TidalOctupolarFMode2")) {
            omega3TidalB_UR = XLALSimUniversalRelationomega03TidalVSlambda3Tidal(lambda3B_UR);
            XLALSimInspiralWaveformParamsInsertTidalOctupolarFMode2(params, omega3TidalB_UR);
        }
        if (!XLALDictContains(params, "dQuadMon1")) {
            quadparam1_UR = XLALSimUniversalRelationQuadMonVSlambda2Tidal(lambda1);
            XLALSimInspiralWaveformParamsInsertdQuadMon1(params, quadparam1_UR - 1.);
        }
        if (!XLALDictContains(params, "dQuadMon2")) {
            quadparam2_UR = XLALSimUniversalRelationQuadMonVSlambda2Tidal(lambda2);
            XLALSimInspiralWaveformParamsInsertdQuadMon2(params, quadparam2_UR - 1.);
        }
        /* Call the waveform driver routine */
        if (approximant == SEOBNRv2T)
            SpinAlignedEOBversion = 201;
        if (approximant == SEOBNRv4T)
            SpinAlignedEOBversion = 401;
        ret = XLALSimIMRSpinAlignedEOBWaveform(hplus, hcross, phiRef, deltaT, m1, m2, f_min, distance, inclination, S1z, S2z, SpinAlignedEOBversion, params);
        break;

    case HGimri:
        /* Waveform-specific sanity checks */
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        if (!checkCOSpinZero(S2x, S2y, S2z))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero CO spin given, but this approximant does not support this case.");
        /* Call the waveform driver */
        ret = XLALHGimriGenerator(hplus, hcross, phiRef, deltaT, m1, m2, f_min, distance, inclination, S1z);
        break;

    case NR_hdf5:
        /* Waveform-specific sanity checks */
        /* Call the waveform driver routine */
        ret = XLALSimInspiralNRWaveformGetHplusHcross(hplus, hcross,
                                                      phiRef, inclination, deltaT, m1, m2, distance, f_min, f_ref, S1x, S1y, S1z,
                                                      S2x, S2y, S2z, XLALSimInspiralWaveformParamsLookupNumRelData(params), XLALSimInspiralWaveformParamsLookupModeArray(params));
        break;

    case NRSur7dq2:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        /* Call the waveform driver routine */
        ret = XLALSimInspiralPrecessingNRSurPolarizations(hplus, hcross, phiRef, inclination, deltaT, m1, m2, distance, f_min, f_ref, S1x, S1y, S1z, S2x, S2y, S2z, params, approximant);
        break;

    case NRSur7dq4:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        /* Call the waveform driver routine */
        ret = XLALSimInspiralPrecessingNRSurPolarizations(hplus, hcross, phiRef, inclination, deltaT, m1, m2, distance, f_min, f_ref, S1x, S1y, S1z, S2x, S2y, S2z, params, approximant);
        break;

    case NRHybSur3dq8:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");

        /* Call the waveform driver routine */
        ret = XLALSimIMRNRHybSur3dq8Polarizations(hplus, hcross, phiRef, inclination, deltaT, m1, m2, distance, f_min, f_ref, S1z, S2z, params);
        break;

    case IMRPhenomXAS:
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        // generate TD waveforms with zero inclincation so that amplitude can be
        // calculated from hplus and hcross, apply inclination-dependent factors
        // in loop below
        polariz = 0;

        ret = XLALSimInspiralTDFromFD(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, deltaT, f_min, f_ref, params, approximant);

        if (ret == XLAL_FAILURE)
            XLAL_ERROR(XLAL_EFUNC);

        // The Fourier domain model is built such that the TD transformation peakds approximately at zero.
        // Here we force an exact alignment at zero by computing the maximum of hp^2 + hc^2.
        maxamp = 0;
        maxind = (*hplus)->data->length - 1;
        for (loopi = (*hplus)->data->length - 1; loopi > -1; loopi--) {
            REAL8 ampsqr = ((*hplus)->data->data[loopi]) * ((*hplus)->data->data[loopi]) + ((*hcross)->data->data[loopi]) * ((*hcross)->data->data[loopi]);
            if (ampsqr > maxamp) {
                maxind = loopi;
                maxamp = ampsqr;
            }
        }
        // Shift peak to t=0.
        XLALGPSSetREAL8(&((*hplus)->epoch), (-1.) * deltaT * maxind);
        XLALGPSSetREAL8(&((*hcross)->epoch), (-1.) * deltaT * maxind);
        break;

    case IMRPhenomXAS_NRTidalv2:
        if( !XLALSimInspiralWaveformParamsFlagsAreDefault(params) )
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
            
        // generate TD waveforms with zero inclincation so that amplitude can be
        // calculated from hplus and hcross, apply inclination-dependent factors
        // in loop below
            polariz = 0;
            
        ret = XLALSimInspiralTDFromFD(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef,
                longAscNodes, eccentricity, meanPerAno, deltaT, f_min, f_ref, params, approximant);
                
            if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
            
            break;

    case IMRPhenomXAS_NRTidalv3:
        if( !XLALSimInspiralWaveformParamsFlagsAreDefault(params) )
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
            
        // generate TD waveforms with zero inclincation so that amplitude can be
        // calculated from hplus and hcross, apply inclination-dependent factors
        // in loop below
            polariz = 0;
            
        ret = XLALSimInspiralTDFromFD(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef,
                longAscNodes, eccentricity, meanPerAno, deltaT, f_min, f_ref, params, approximant);
                
            if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
            
            break;       

    case IMRPhenomXHM:
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        polariz = 0;

        ret = XLALSimInspiralTDFromFD(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, deltaT, f_min, f_ref, params, approximant);

        if (ret == XLAL_FAILURE)
            XLAL_ERROR(XLAL_EFUNC);

        // The Fourier domain model is built such that the TD transformation peakds approximately at zero.
        // Here we force an exact alignment at zero by computing the maximum of hp^2 + hc^2.
        maxamp = 0;
        maxind = (*hplus)->data->length - 1;
        for (loopi = (*hplus)->data->length - 1; loopi > -1; loopi--) {
            REAL8 ampsqr = ((*hplus)->data->data[loopi]) * ((*hplus)->data->data[loopi]) + ((*hcross)->data->data[loopi]) * ((*hcross)->data->data[loopi]);
            if (ampsqr > maxamp) {
                maxind = loopi;
                maxamp = ampsqr;
            }
        }
        // Shift peak to t=0.
        XLALGPSSetREAL8(&((*hplus)->epoch), (-1.) * deltaT * maxind);
        XLALGPSSetREAL8(&((*hcross)->epoch), (-1.) * deltaT * maxind);
        break;

    case IMRPhenomXP:
        polariz = 0;
        ret = XLALSimInspiralTDFromFD(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, deltaT, f_min, f_ref, params, approximant);
        break;
    case IMRPhenomXP_NRTidalv2:
        polariz = 0;
        ret = XLALSimInspiralTDFromFD(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, deltaT, f_min, f_ref, params, approximant);
        break;

    case IMRPhenomXP_NRTidalv3:
        polariz = 0;
        ret = XLALSimInspiralTDFromFD(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, deltaT, f_min, f_ref, params, approximant);
        break;

    case IMRPhenomXPHM:
        polariz = 0;
        ret = XLALSimInspiralTDFromFD(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, deltaT, f_min, f_ref, params, approximant);
        break;

    case IMRPhenomT:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");

        ret = XLALSimIMRPhenomT(hplus, hcross, m1, m2, S1z, S2z, distance, inclination, deltaT, f_min, f_ref, phiRef, params);
        break;


    case IMRPhenomTHM:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");

        ret = XLALSimIMRPhenomTHM(hplus, hcross, m1, m2, S1z, S2z, distance, inclination, deltaT, f_min, f_ref, phiRef, params);
        break;

    case IMRPhenomTP:
        /* Waveform-specific sanity checks */
        /* FIXME: CHECK ADDITIONAL CHECKS OF XP */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");

        ret = XLALSimIMRPhenomTP(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, deltaT, f_min, f_ref, phiRef, params);
        break;

    case IMRPhenomTPHM:
        /* Waveform-specific sanity checks */
        /* FIXME: CHECK ADDITIONAL CHECKS OF XPHM */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");

        ret = XLALSimIMRPhenomTPHM(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, deltaT, f_min, f_ref, phiRef, params);
        break;

    case IMRPhenomXO4a:
    	 		polariz = 0;
    			ret = XLALSimInspiralTDFromFD(hplus, hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, longAscNodes, eccentricity, meanPerAno, deltaT, f_min, f_ref, params, approximant);
    			break;

    default:
        XLALPrintError("TD version of approximant not implemented in lalsimulation\n");
        XLAL_ERROR(XLAL_EINVAL);
    }
    //R.C.: here's the reference explaining why we perform this rotation https://dcc.ligo.org/LIGO-G1900275
    if (polariz && (*hplus) && (*hcross)) {
        REAL8 tmpP, tmpC;
        REAL8 cp = cos(2. * polariz);
        REAL8 sp = sin(2. * polariz);
        for (UINT4 idx = 0; idx < (*hplus)->data->length; idx++) {
            tmpP = (*hplus)->data->data[idx];
            tmpC = (*hcross)->data->data[idx];
            (*hplus)->data->data[idx] = cp * tmpP + sp * tmpC;
            (*hcross)->data->data[idx] = cp * tmpC - sp * tmpP;
        }
    }

    if (ret == XLAL_FAILURE)
        XLAL_ERROR(XLAL_EFUNC);

    return ret;
}

/**
 * Copy of the old code of XLALSimInspiralChooseFDWaveform(). The new version of XLALSimInspiralChooseFDWaveform() is just a wrapper over XLALSimInspiralGenerateFDWaveform.
 * XLALSimInspiralGenerateFDWaveform() internally calls this function for legacy approximants.
 */
static int XLALSimInspiralChooseFDWaveform_legacy(
  COMPLEX16FrequencySeries **hptilde,     /* FD plus polarization */
  COMPLEX16FrequencySeries **hctilde,     /* FD cross polarization */
  REAL8 m1,                         /* mass of companion 1 (kg) */
  REAL8 m2,                         /* mass of companion 2 (kg) */
  REAL8 S1x,                        /* x-component of the dimensionless spin of object 1 */
  REAL8 S1y,                        /* y-component of the dimensionless spin of object 1 */
  REAL8 S1z,                        /* z-component of the dimensionless spin of object 1 */
  REAL8 S2x,                        /* x-component of the dimensionless spin of object 2 */
  REAL8 S2y,                        /* y-component of the dimensionless spin of object 2 */
  REAL8 S2z,                        /* z-component of the dimensionless spin of object 2 */
  REAL8 distance,                   /* distance of source (m) */
  REAL8 inclination,                /* inclination of source (rad) */
  REAL8 phiRef,                     /* reference orbital phase (rad) */
  REAL8 longAscNodes,               /* longitude of ascending nodes, degenerate with the polarization angle, Omega in documentation */
  REAL8 eccentricity,               /* eccentricity at reference epoch */
  REAL8 UNUSED meanPerAno,          /* mean anomaly of periastron */
  // frequency sampling parameters, no default value
  REAL8 deltaF,                     /* sampling interval (Hz) */
  REAL8 f_min,                      /* starting GW frequency (Hz) */
  REAL8 f_max,                      /* ending GW frequency (Hz) */
  REAL8 f_ref,                      /* Reference frequency (Hz) */
  LALDict *params,                  /* LAL dictionary containing accessory parameters */
  Approximant approximant           /* post-Newtonian approximant to use for waveform production */
)
{
    REAL8 LNhatx, LNhaty, LNhatz;
    REAL8 tmp1, tmp2;
    REAL8 E1x, E1y, E1z;
    REAL8 kMax;
    REAL8 v0, fStart;
    int ret;
    unsigned int j;
    REAL8 pfac, cfac;
    INT4 phiRefAtEnd;
    int amplitudeO = XLALSimInspiralWaveformParamsLookupPNAmplitudeOrder(params);
    int phaseO = XLALSimInspiralWaveformParamsLookupPNPhaseOrder(params);

    REAL8 quadparam1 = 1. + XLALSimInspiralWaveformParamsLookupdQuadMon1(params);
    REAL8 quadparam2 = 1. + XLALSimInspiralWaveformParamsLookupdQuadMon2(params);
    REAL8 lambda1 = XLALSimInspiralWaveformParamsLookupTidalLambda1(params);
    REAL8 lambda2 = XLALSimInspiralWaveformParamsLookupTidalLambda2(params);

    /* Support variables for precessing wfs */
    REAL8 spin1x, spin1y, spin1z;
    REAL8 spin2x, spin2y, spin2z;

    /* Variables for IMRPhenomP and IMRPhenomPv2 */
    REAL8 chi1_l, chi2_l, chip, thetaJN, alpha0, phi_aligned, zeta_polariz;
    COMPLEX16 PhPpolp, PhPpolc;

    /* General sanity checks that will abort
     *
     * If non-GR approximants are added, include them in
     * XLALSimInspiralApproximantAcceptTestGRParams()
     */
    if (!XLALSimInspiralWaveformParamsNonGRAreDefault(params) && XLALSimInspiralApproximantAcceptTestGRParams(approximant) != LAL_SIM_INSPIRAL_TESTGR_PARAMS) {
        XLALPrintError("XLAL Error - %s: Passed in non-NULL pointer to LALSimInspiralTestGRParam for an approximant that does not use LALSimInspiralTestGRParam\n", __func__);
        XLAL_ERROR(XLAL_EINVAL);
    }

    /* General sanity check the input parameters - only give warnings! */
    if (deltaF > 1.)
        XLALPrintWarning("XLAL Warning - %s: Large value of deltaF = %e requested...This corresponds to a very short TD signal (with padding). Consider a smaller value.\n", __func__, deltaF);
    if (deltaF < 1. / 4096.)
        XLALPrintWarning("XLAL Warning - %s: Small value of deltaF = %e requested...This corresponds to a very long TD signal. Consider a larger value.\n", __func__, deltaF);
    if (m1 < 0.09 * LAL_MSUN_SI)
        XLALPrintWarning("XLAL Warning - %s: Small value of m1 = %e (kg) = %e (Msun) requested...Perhaps you have a unit conversion error?\n", __func__, m1, m1 / LAL_MSUN_SI);
    if (m2 < 0.09 * LAL_MSUN_SI)
        XLALPrintWarning("XLAL Warning - %s: Small value of m2 = %e (kg) = %e (Msun) requested...Perhaps you have a unit conversion error?\n", __func__, m2, m2 / LAL_MSUN_SI);
    if (m1 + m2 > 1000. * LAL_MSUN_SI)
        XLALPrintWarning("XLAL Warning - %s: Large value of total mass m1+m2 = %e (kg) = %e (Msun) requested...Signal not likely to be in band of ground-based detectors.\n", __func__, m1 + m2, (m1 + m2) / LAL_MSUN_SI);
    if (S1x * S1x + S1y * S1y + S1z * S1z > 1.000001)
        XLALPrintWarning("XLAL Warning - %s: S1 = (%e,%e,%e) with norm > 1 requested...Are you sure you want to violate the Kerr bound?\n", __func__, S1x, S1y, S1z);
    if (S2x * S2x + S2y * S2y + S2z * S2z > 1.000001)
        XLALPrintWarning("XLAL Warning - %s: S2 = (%e,%e,%e) with norm > 1 requested...Are you sure you want to violate the Kerr bound?\n", __func__, S2x, S2y, S2z);
    if (f_min < 1.)
        XLALPrintWarning("XLAL Warning - %s: Small value of fmin = %e requested...Check for errors, this could create a very long waveform.\n", __func__, f_min);
    if (f_min > 40.000001)
        XLALPrintWarning("XLAL Warning - %s: Large value of fmin = %e requested...Check for errors, the signal will start in band.\n", __func__, f_min);

    /* adjust the reference frequency for certain precessing approximants:
     * if that approximate interprets f_ref==0 to be f_min, set f_ref=f_min;
     * otherwise do nothing */
    FIX_REFERENCE_FREQUENCY(f_ref, f_min, approximant);

    /* The non-precessing waveforms return h(f) for optimal orientation
     * (i=0, Fp=1, Fc=0; Lhat pointed toward the observer)
     * To get generic polarizations we multiply by inclination dependence
     * and note hc(f) \propto -I * hp(f)
     * Non-precessing waveforms multiply hp by pfac, hc by -I*cfac
     */
    cfac = cos(inclination);
    pfac = 0.5 * (1. + cfac * cfac);

    switch (approximant) {
        /* inspiral-only models */
    case EccentricFD:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFrameAxisIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
        if (!XLALSimInspiralWaveformParamsModesChoiceIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        /* Call the waveform driver routine */
        /* Note that for generic inclined eccentric waveforms
         *             * it is not possible to decompose hc(f) \propto I * hp(f)
         *                         * we call both polarizations independently
         *                                     */
        /*ret = XLALSimInspiralEFD(hptilde, hctilde, phiRef, deltaF, m1, m2,
         *             f_min, f_max, i, r, lambda1, lambda2, phaseO);*/
        // IMPORTANT CHECK: verify that inclination_azimuth is the longitude of ascending nodes
        ret = XLALSimInspiralEFD(hptilde, hctilde, phiRef, deltaF, m1, m2, f_min, f_max, inclination, distance, longAscNodes, eccentricity, phaseO);
        if (ret == XLAL_FAILURE)
            XLAL_ERROR(XLAL_EFUNC);
        break;

    case TaylorF2:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFrameAxisIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
        if (!XLALSimInspiralWaveformParamsModesChoiceIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");

        /* Call the waveform driver routine */
        ret = XLALSimInspiralTaylorF2(hptilde, phiRef, deltaF, m1, m2, S1z, S2z, f_min, f_max, f_ref, distance, params);
        if (ret == XLAL_FAILURE)
            XLAL_ERROR(XLAL_EFUNC);
        /* Produce both polarizations */
        *hctilde = XLALCreateCOMPLEX16FrequencySeries("FD hcross", &((*hptilde)->epoch), (*hptilde)->f0, (*hptilde)->deltaF, &((*hptilde)->sampleUnits), (*hptilde)->data->length);
        for (j = 0; j < (*hptilde)->data->length; j++) {
            (*hctilde)->data->data[j] = -I * cfac * (*hptilde)->data->data[j];
            (*hptilde)->data->data[j] *= pfac;
        }
        break;

    case TaylorF2Ecc:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFrameAxisIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
        if (!XLALSimInspiralWaveformParamsModesChoiceIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        REAL8 f_ecc = XLALSimInspiralWaveformParamsLookupEccentricityFreq(params);     /** get f_ecc */
        if (eccentricity > 0.0 && eccentricity < 1.0 && f_ecc < 0.0) {
            /* we set f_ecc to be f_ref for correct eccentricity but not specifying f_ecc. */
            f_ecc = f_ref;
            if (f_ecc == 0)
                f_ecc = f_min;
            XLALSimInspiralWaveformParamsInsertEccentricityFreq(params, f_ecc);
            XLAL_PRINT_WARNING("Warning... The reference frequency for eccentricity was set as default value(%f). This might be not optimal case for you.\n", f_ecc);
        }

        /* Call the waveform driver routine */
        ret = XLALSimInspiralTaylorF2Ecc(hptilde, phiRef, deltaF, m1, m2, S1z, S2z, f_min, f_max, f_ref, distance, eccentricity, params);
        if (ret == XLAL_FAILURE)
            XLAL_ERROR(XLAL_EFUNC);
        /* Produce both polarizations */
        *hctilde = XLALCreateCOMPLEX16FrequencySeries("FD hcross", &((*hptilde)->epoch), (*hptilde)->f0, (*hptilde)->deltaF, &((*hptilde)->sampleUnits), (*hptilde)->data->length);
        for (j = 0; j < (*hptilde)->data->length; j++) {
            (*hctilde)->data->data[j] = -I * cfac * (*hptilde)->data->data[j];
            (*hptilde)->data->data[j] *= pfac;
        }
        break;

    case TaylorF2NLTides:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFrameAxisIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
        if (!XLALSimInspiralWaveformParamsModesChoiceIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");

        // FIXME : add checks for NL tidal parameters?

        /* Call the waveform driver routine */
        ret = XLALSimInspiralTaylorF2NLTides(hptilde, phiRef, deltaF, m1, m2, S1z, S2z, f_min, f_max, f_ref, distance, params);
        if (ret == XLAL_FAILURE)
            XLAL_ERROR(XLAL_EFUNC);
        /* Produce both polarizations */
        *hctilde = XLALCreateCOMPLEX16FrequencySeries("FD hcross", &((*hptilde)->epoch), (*hptilde)->f0, (*hptilde)->deltaF, &((*hptilde)->sampleUnits), (*hptilde)->data->length);
        for (j = 0; j < (*hptilde)->data->length; j++) {
            (*hctilde)->data->data[j] = -I * cfac * (*hptilde)->data->data[j];
            (*hptilde)->data->data[j] *= pfac;
        }
        break;

        /* non-spinning inspiral-merger-ringdown models */
    case IMRPhenomA:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero spins were given, but this is a non-spinning approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        /* Call the waveform driver routine */
        ret = XLALSimIMRPhenomAGenerateFD(hptilde, phiRef, deltaF, m1, m2, f_min, f_max, distance);
        if (ret == XLAL_FAILURE)
            XLAL_ERROR(XLAL_EFUNC);
        /* Produce both polarizations */
        *hctilde = XLALCreateCOMPLEX16FrequencySeries("FD hcross", &((*hptilde)->epoch), (*hptilde)->f0, (*hptilde)->deltaF, &((*hptilde)->sampleUnits), (*hptilde)->data->length);
        for (j = 0; j < (*hptilde)->data->length; j++) {
            (*hctilde)->data->data[j] = -I * cfac * (*hptilde)->data->data[j];
            (*hptilde)->data->data[j] *= pfac;
        }
        break;

        /* spinning inspiral-only models */
    case SpinTaylorF2:
        /* Waveform-specific sanity checks */
        /* Sanity check unused fields of params */
        if (!XLALSimInspiralWaveformParamsFrameAxisIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
        if (!XLALSimInspiralWaveformParamsModesChoiceIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        if (!checkCOSpinZero(S2x, S2y, S2z))    // This is a single-spin model
            XLAL_ERROR(XLAL_EINVAL, "Non-zero CO spin given, but this approximant does not support this case.");
        spin1x = S1x;
        spin1y = S1y;
        spin1z = S1z;
        ROTATEY(inclination, spin1x, spin1y, spin1z);
        LNhatx = sin(inclination);
        LNhaty = 0.;
        LNhatz = cos(inclination);
        /* Maximum PN amplitude order for precessing waveforms is
         * MAX_PRECESSING_AMP_PN_ORDER */
        amplitudeO = 0;         /* amplitudeO <= MAX_PRECESSING_AMP_PN_ORDER ?
                                   amplitudeO : MAX_PRECESSING_AMP_PN_ORDER */ ;
        /* Call the waveform driver routine */
        ret = XLALSimInspiralSpinTaylorF2(hptilde, hctilde, phiRef, deltaF, m1, m2, spin1x, spin1y, spin1z, LNhatx, LNhaty, LNhatz, f_min, f_max, f_ref, distance, params, phaseO, amplitudeO);
        if (ret == XLAL_FAILURE)
            XLAL_ERROR(XLAL_EFUNC);
        break;

        /* FIXME: Comment out this case, as I don't have its source code */
        //case TaylorR2F4:
        //    /* Waveform-specific sanity checks */
        //    if( !XLALSimInspiralWaveformParamsFlagsAreDefault(params) )
        //        XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        //    if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
        //        XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        //    /* Call the waveform driver routine */
        //    ret = XLALSimInspiralTaylorR2F4(hptilde, phiRef, deltaF, m1, m2,
        //            S1z, S2z, f_min, r, phaseO, amplitudeO);
        //    break;

    case TaylorF2RedSpin:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        /* Call the waveform driver routine */
        ret = XLALSimInspiralTaylorF2ReducedSpin(hptilde, phiRef, deltaF, m1, m2, XLALSimInspiralTaylorF2ReducedSpinComputeChi(m1, m2, S1z, S2z), f_min, f_max, distance, phaseO, amplitudeO);
        if (ret == XLAL_FAILURE)
            XLAL_ERROR(XLAL_EFUNC);
        /* Produce both polarizations */
        *hctilde = XLALCreateCOMPLEX16FrequencySeries("FD hcross", &((*hptilde)->epoch), (*hptilde)->f0, (*hptilde)->deltaF, &((*hptilde)->sampleUnits), (*hptilde)->data->length);
        for (j = 0; j < (*hptilde)->data->length; j++) {
            (*hctilde)->data->data[j] = -I * cfac * (*hptilde)->data->data[j];
            (*hptilde)->data->data[j] *= pfac;
        }
        break;

    case TaylorF2RedSpinTidal:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        /* Call the waveform driver routine */
        ret = XLALSimInspiralTaylorF2ReducedSpinTidal(hptilde, phiRef, deltaF, m1, m2, XLALSimIMRPhenomBComputeChi(m1, m2, S1z, S2z), lambda1, lambda2, f_min, f_max, distance, phaseO, amplitudeO);
        if (ret == XLAL_FAILURE)
            XLAL_ERROR(XLAL_EFUNC);
        /* Produce both polarizations */
        *hctilde = XLALCreateCOMPLEX16FrequencySeries("FD hcross", &((*hptilde)->epoch), (*hptilde)->f0, (*hptilde)->deltaF, &((*hptilde)->sampleUnits), (*hptilde)->data->length);
        for (j = 0; j < (*hptilde)->data->length; j++) {
            (*hctilde)->data->data[j] = -I * cfac * (*hptilde)->data->data[j];
            (*hptilde)->data->data[j] *= pfac;
        }
        break;

        /* spinning inspiral-merger-ringdown models */
    case IMRPhenomB:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        /* Call the waveform driver routine */
        ret = XLALSimIMRPhenomBGenerateFD(hptilde, phiRef, deltaF, m1, m2, XLALSimIMRPhenomBComputeChi(m1, m2, S1z, S2z), f_min, f_max, distance);
        if (ret == XLAL_FAILURE)
            XLAL_ERROR(XLAL_EFUNC);
        /* Produce both polarizations */
        *hctilde = XLALCreateCOMPLEX16FrequencySeries("FD hcross", &((*hptilde)->epoch), (*hptilde)->f0, (*hptilde)->deltaF, &((*hptilde)->sampleUnits), (*hptilde)->data->length);
        for (j = 0; j < (*hptilde)->data->length; j++) {
            (*hctilde)->data->data[j] = -I * cfac * (*hptilde)->data->data[j];
            (*hptilde)->data->data[j] *= pfac;
        }
        break;

    case IMRPhenomC:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        /* Call the waveform driver routine */
        ret = XLALSimIMRPhenomCGenerateFD(hptilde, phiRef, deltaF, m1, m2, XLALSimIMRPhenomBComputeChi(m1, m2, S1z, S2z), f_min, f_max, distance, params);
        if (ret == XLAL_FAILURE)
            XLAL_ERROR(XLAL_EFUNC);
        /* Produce both polarizations */
        *hctilde = XLALCreateCOMPLEX16FrequencySeries("FD hcross", &((*hptilde)->epoch), (*hptilde)->f0, (*hptilde)->deltaF, &((*hptilde)->sampleUnits), (*hptilde)->data->length);
        for (j = 0; j < (*hptilde)->data->length; j++) {
            (*hctilde)->data->data[j] = -I * cfac * (*hptilde)->data->data[j];
            (*hptilde)->data->data[j] *= pfac;
        }
        break;

    case IMRPhenomD:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        /* Call the waveform driver routine */

        ret = XLALSimIMRPhenomDGenerateFD(hptilde, phiRef, f_ref, deltaF, m1, m2, S1z, S2z, f_min, f_max, distance, params, NoNRT_V);
        if (ret == XLAL_FAILURE)
            XLAL_ERROR(XLAL_EFUNC);
        /* Produce both polarizations */
        *hctilde = XLALCreateCOMPLEX16FrequencySeries("FD hcross", &((*hptilde)->epoch), (*hptilde)->f0, (*hptilde)->deltaF, &((*hptilde)->sampleUnits), (*hptilde)->data->length);
        for (j = 0; j < (*hptilde)->data->length; j++) {
            (*hctilde)->data->data[j] = -I * cfac * (*hptilde)->data->data[j];
            (*hptilde)->data->data[j] *= pfac;
        }
        break;

    case IMRPhenomD_NRTidal:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (lambda1 < 0 || lambda2 < 0)
            XLAL_ERROR(XLAL_EFUNC, "lambda1 = %f, lambda2 = %f. Both should be greater than zero for IMRPhenomD_NRTidal", lambda1, lambda2);
        /* Call the waveform driver routine */
        ret = XLALSimIMRPhenomDNRTidal(hptilde, phiRef, deltaF, f_min, f_max, f_ref, distance, m1, m2, S1z, S2z, lambda1, lambda2, params, NRTidal_V);
        if (ret == XLAL_FAILURE)
            XLAL_ERROR(XLAL_EFUNC);
        /* Produce both polarizations */
        *hctilde = XLALCreateCOMPLEX16FrequencySeries("FD hcross", &((*hptilde)->epoch), (*hptilde)->f0, (*hptilde)->deltaF, &((*hptilde)->sampleUnits), (*hptilde)->data->length);
        for (j = 0; j < (*hptilde)->data->length; j++) {
            (*hctilde)->data->data[j] = -I * cfac * (*hptilde)->data->data[j];
            (*hptilde)->data->data[j] *= pfac;
        }
        break;

    case IMRPhenomD_NRTidalv2:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (lambda1 < 0 || lambda2 < 0)
            XLAL_ERROR(XLAL_EFUNC, "lambda1 = %f, lambda2 = %f. Both should be greater than zero for IMRPhenomD_NRTidalv2", lambda1, lambda2);
        ret = XLALSimInspiralSetQuadMonParamsFromLambdas(params);
        XLAL_CHECK(XLAL_SUCCESS == ret, ret, "Failed to set QuadMon from Lambdas for IMRPhenomD_NRTidalv2");
        /* Call the waveform driver routine */
        ret = XLALSimIMRPhenomDNRTidal(hptilde, phiRef, deltaF, f_min, f_max, f_ref, distance, m1, m2, S1z, S2z, lambda1, lambda2, params, NRTidalv2_V);
        if (ret == XLAL_FAILURE)
            XLAL_ERROR(XLAL_EFUNC);
        /* Produce both polarizations */
        *hctilde = XLALCreateCOMPLEX16FrequencySeries("FD hcross", &((*hptilde)->epoch), (*hptilde)->f0, (*hptilde)->deltaF, &((*hptilde)->sampleUnits), (*hptilde)->data->length);
        for (j = 0; j < (*hptilde)->data->length; j++) {
            (*hctilde)->data->data[j] = -I * cfac * (*hptilde)->data->data[j];
            (*hptilde)->data->data[j] *= pfac;
        }
        break;

    case IMRPhenomNSBH:
        /* Waveform-specific sanity checks */
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (lambda1 != 0 || lambda2 < 0)
            XLAL_ERROR(XLAL_EDOM, "lambda1 = %f, lambda2 = %f. lambda1 should be equal to zero and lambda2 should be greater than or equal to zero for IMRPhenomNSBH", lambda1, lambda2);
        /* Call the waveform driver routine */
        ret = XLALSimIMRPhenomNSBH(hptilde, phiRef, deltaF, f_min, f_max, f_ref, distance, m1, m2, S1z, S2z, params);
        if (ret == XLAL_FAILURE)
            XLAL_ERROR(XLAL_EFUNC);
        /* Produce both polarizations */
        *hctilde = XLALCreateCOMPLEX16FrequencySeries("FD hcross", &((*hptilde)->epoch), (*hptilde)->f0, (*hptilde)->deltaF, &((*hptilde)->sampleUnits), (*hptilde)->data->length);
        for (j = 0; j < (*hptilde)->data->length; j++) {
            (*hctilde)->data->data[j] = -I * cfac * (*hptilde)->data->data[j];
            (*hptilde)->data->data[j] *= pfac;
        }
        break;

    case IMRPhenomHM:
        /* Waveform-specific sanity checks */
        // if( !XLALSimInspiralWaveformParamsFlagsAreDefault(params) )
        // XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        /* Call the waveform driver routine */

        REAL8Sequence *freqs = XLALCreateREAL8Sequence(2);
        freqs->data[0] = f_min;
        freqs->data[1] = f_max;
        ret = XLALSimIMRPhenomHM(hptilde, hctilde, freqs, m1, m2, S1z, S2z, distance, inclination, phiRef, deltaF, f_ref, params);
        if (ret == XLAL_FAILURE)
            XLAL_ERROR(XLAL_EFUNC);
        XLALDestroyREAL8Sequence(freqs);
        break;

    case EOBNRv2_ROM:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero spins were given, but this is a non-spinning approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");

        ret = XLALSimIMREOBNRv2HMROM(hptilde, hctilde, phiRef, deltaF, f_min, f_max, f_ref, distance, inclination, m1, m2, 0);
        break;

    case EOBNRv2HM_ROM:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero spins were given, but this is a non-spinning approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");

        ret = XLALSimIMREOBNRv2HMROM(hptilde, hctilde, phiRef, deltaF, f_min, f_max, f_ref, distance, inclination, m1, m2, 1);
        break;

    case SEOBNRv1_ROM_EffectiveSpin:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        if (!checkAlignedSpinsEqual(S1z, S2z)) {
            XLALPrintError("XLAL Error - %s: SEOBNRv1ROM Effective Spin model called with unequal aligned spins: %lf, %lf.\n", __func__, S1z, S2z);
            XLAL_ERROR(XLAL_EINVAL);
        }

        ret = XLALSimIMRSEOBNRv1ROMEffectiveSpin(hptilde, hctilde, phiRef, deltaF, f_min, f_max, f_ref, distance, inclination, m1, m2, XLALSimIMRPhenomBComputeChi(m1, m2, S1z, S2z));
        break;

    case SEOBNRv1_ROM_DoubleSpin:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");

        ret = XLALSimIMRSEOBNRv1ROMDoubleSpin(hptilde, hctilde, phiRef, deltaF, f_min, f_max, f_ref, distance, inclination, m1, m2, S1z, S2z);
        break;

    case SEOBNRv2_ROM_EffectiveSpin:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        if (!checkAlignedSpinsEqual(S1z, S2z)) {
            XLALPrintError("XLAL Error - %s: SEOBNRv2ROM Effective Spin model called with unequal aligned spins: %lf, %lf.\n", __func__, S1z, S2z);
            XLAL_ERROR(XLAL_EINVAL);
        }

        ret = XLALSimIMRSEOBNRv2ROMEffectiveSpin(hptilde, hctilde, phiRef, deltaF, f_min, f_max, f_ref, distance, inclination, m1, m2, XLALSimIMRPhenomBComputeChi(m1, m2, S1z, S2z));
        break;

    case SEOBNRv2_ROM_DoubleSpin:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");

        ret = XLALSimIMRSEOBNRv2ROMDoubleSpin(hptilde, hctilde, phiRef, deltaF, f_min, f_max, f_ref, distance, inclination, m1, m2, S1z, S2z);
        break;

    case SEOBNRv2_ROM_DoubleSpin_HI:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");

        ret = XLALSimIMRSEOBNRv2ROMDoubleSpinHI(hptilde, hctilde, phiRef, deltaF, f_min, f_max, f_ref, distance, inclination, m1, m2, S1z, S2z, -1);
        break;

    case SEOBNRv4_ROM:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");

        ret = XLALSimIMRSEOBNRv4ROM(hptilde, hctilde, phiRef, deltaF, f_min, f_max, f_ref, distance, inclination, m1, m2, S1z, S2z, -1, params, NoNRT_V);
        break;

    case SEOBNRv4HM_ROM:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");

        ret = XLALSimIMRSEOBNRv4HMROM(hptilde, hctilde, phiRef, deltaF, f_min, f_max, f_ref, distance, inclination, m1, m2, S1z, S2z, -1, 5, 1, params);
        break;

    case SEOBNRv5_ROM:
        /* Waveform-specific sanity checks */
        if( !XLALSimInspiralWaveformParamsFlagsAreDefault(params) )
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if( !checkTidesZero(lambda1, lambda2) )
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");

        ret = XLALSimIMRSEOBNRv5HMROM(hptilde, hctilde,
                phiRef, deltaF, f_min, f_max, f_ref, distance, inclination, m1, m2, S1z, S2z, -1, 1, true, params, NoNRT_V);
        break;

    case SEOBNRv5HM_ROM:
        /* Waveform-specific sanity checks */
        if( !XLALSimInspiralWaveformParamsFlagsAreDefault(params) )
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if( !checkTidesZero(lambda1, lambda2) )
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");

        LALValue *mode_arr = NULL;
        LALDict *pars_aux;
        INT2Sequence *modeseq = NULL;
        UINT4 nmodes = 7;
        UINT2 eobmodesv5hm = 7;
    
        if(params == NULL){
            pars_aux = XLALCreateDict();
        }
        else{
            pars_aux = XLALDictDuplicate(params);
        }
        mode_arr = XLALSimInspiralWaveformParamsLookupModeArray(pars_aux);

        if(mode_arr != NULL)
        {
            modeseq = XLALSimInspiralModeArrayReadModes(mode_arr);
            nmodes = modeseq->length/2;
            if (
                nmodes == 2 &&
                modeseq->data[0] == 2 && abs(modeseq->data[1]) == 2 &&
                modeseq->data[2] == 3 && abs(modeseq->data[3]) == 3)
            {   
                eobmodesv5hm = 2;
            }
            if (
                nmodes == 3 &&
                modeseq->data[0] == 2 && abs(modeseq->data[1]) == 2 &&
                modeseq->data[4] == 3 && abs(modeseq->data[5]) == 3 && 
                modeseq->data[2] == 2 && abs(modeseq->data[3]) == 1)
            {   
                eobmodesv5hm = 3;
            }
            if (
                nmodes == 4 &&
                modeseq->data[0] == 2 && abs(modeseq->data[1]) == 2 &&
                modeseq->data[4] == 3 && abs(modeseq->data[5]) == 3 && 
                modeseq->data[2] == 2 && abs(modeseq->data[3]) == 1 && 
                modeseq->data[6] == 4 && abs(modeseq->data[7]) == 4)
            {   
                eobmodesv5hm = 4;
            }
            if (
                nmodes == 5 &&
                modeseq->data[0] == 2 && abs(modeseq->data[1]) == 2 &&
                modeseq->data[4] == 3 && abs(modeseq->data[5]) == 3 && 
                modeseq->data[2] == 2 && abs(modeseq->data[3]) == 1 && 
                modeseq->data[6] == 4 && abs(modeseq->data[7]) == 4 && 
                modeseq->data[8] == 5 && abs(modeseq->data[9]) == 5)
            {
                eobmodesv5hm = 5;
            }
            if (
                nmodes == 6 &&
                modeseq->data[0] == 2 && abs(modeseq->data[1]) == 2 &&
                modeseq->data[4] == 3 && abs(modeseq->data[5]) == 3 && 
                modeseq->data[2] == 2 && abs(modeseq->data[3]) == 1 && 
                modeseq->data[8] == 4 && abs(modeseq->data[9]) == 4 && 
                modeseq->data[10] == 5 && abs(modeseq->data[11]) == 5 && 
                modeseq->data[6] == 3 && abs(modeseq->data[7]) == 2)
            {
                eobmodesv5hm = 6;
            }
        }
        ret = XLALSimIMRSEOBNRv5HMROM(hptilde, hctilde,
            phiRef, deltaF, f_min, f_max, f_ref, distance, inclination, m1, m2, S1z, S2z, -1, eobmodesv5hm, true, params, NoNRT_V);
        XLALDestroyINT2Sequence(modeseq);
        XLALDestroyValue(mode_arr);
        XLALDestroyDict(pars_aux);
        break;

    case SEOBNRv4_ROM_NRTidal:

        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (lambda1 < 0 || lambda2 < 0)
            XLAL_ERROR(XLAL_EFUNC, "lambda1 = %f, lambda2 = %f. Both should be greater than zero for SEOBNRv4_ROM_NRTidal", lambda1, lambda2);
        ret = XLALSimInspiralSetQuadMonParamsFromLambdas(params);
        XLAL_CHECK(XLAL_SUCCESS == ret, ret, "Failed to set QuadMon from Lambdas for SEOBNRv4_ROM_NRTidal");
        ret = XLALSimIMRSEOBNRv4ROMNRTidal(hptilde, hctilde, phiRef, deltaF, f_min, f_max, f_ref, distance, inclination, m1, m2, S1z, S2z, lambda1, lambda2, params, NRTidal_V);
        break;

    case SEOBNRv4_ROM_NRTidalv2:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (lambda1 < 0 || lambda2 < 0)
            XLAL_ERROR(XLAL_EFUNC, "lambda1 = %f, lambda2 = %f. Both should be greater than zero for SEOBNRv4_ROM_NRTidal", lambda1, lambda2);
        ret = XLALSimInspiralSetQuadMonParamsFromLambdas(params);
        XLAL_CHECK(XLAL_SUCCESS == ret, ret, "Failed to set QuadMon from Lambdas for SEOBNRv4_ROM_NRTidalv2");
        ret = XLALSimIMRSEOBNRv4ROMNRTidal(hptilde, hctilde, phiRef, deltaF, f_min, f_max, f_ref, distance, inclination, m1, m2, S1z, S2z, lambda1, lambda2, params, NRTidalv2_V);
        break;

    case SEOBNRv4_ROM_NRTidalv2_NSBH:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (m1 < m2)
            XLAL_ERROR(XLAL_EDOM, "m1 = %e, m2=%e. m1 should be greater than or equal to m2 for SEOBNRv4_ROM_NRTidalv2_NSBH", m1, m2);
        if (lambda1 != 0)
            XLAL_ERROR(XLAL_EDOM, "lambda1 = %f. lambda1 should be zero for SEOBNRv4_ROM_NRTidalv2_NSBH", lambda1);
        if (lambda2 < 0)
            XLAL_ERROR(XLAL_EDOM, "lambda2 = %f. lambda2 should be nonnegative for SEOBNRv4_ROM_NRTidalv2_NSBH", lambda2);
        if (lambda2 > 5000)
            XLAL_ERROR(XLAL_EDOM, "lambda2 = %f. lambda2 should be < 5000", lambda2);
        if (S2z != 0)
            XLAL_PRINT_WARNING("WARNING: S2z = %f. SEOBNRv4_ROM_NRTidalv2_NSBH is calibrated to NR data for which the NS spin is zero.", S2z);
        if (m2 < 1 * LAL_MSUN_SI)
            XLAL_PRINT_WARNING("WARNING: m2=%e MSun. SEOBNRv4_ROM_NRTidalv2_NSBH is calibrated to NR data for which the NS mass is >=1 solar mass.", m2 / LAL_MSUN_SI);
        if (m2 > 3 * LAL_MSUN_SI)
            XLAL_ERROR(XLAL_EDOM, "m2=%e Msun. NS Mass should be <=3 solar masses", m2 / LAL_MSUN_SI);
        if (m1 / m2 > 100)
            XLAL_ERROR(XLAL_EDOM, "m1/m2=%e mass ratio should be < 100", m1 / m2);

        ret = XLALSimIMRSEOBNRv4ROMNRTidal(hptilde, hctilde, phiRef, deltaF, f_min, f_max, f_ref, distance, inclination, m1, m2, S1z, S2z, lambda1, lambda2, params, NRTidalv2NSBH_V);
        break;

    case SEOBNRv4T_surrogate:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");

        ret = XLALSimIMRSEOBNRv4TSurrogate(hptilde, hctilde, phiRef, deltaF, f_min, f_max, f_ref, distance, inclination, m1, m2, S1z, S2z, lambda1, lambda2, SEOBNRv4TSurrogate_CUBIC);
        break;

    case SEOBNRv5_ROM_NRTidalv3:

        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (lambda1 < 0 || lambda2 < 0)
            XLAL_ERROR(XLAL_EFUNC, "lambda1 = %f, lambda2 = %f. Both should be greater than zero for SEOBNRv5_ROM_NRTidal", lambda1, lambda2);
        ret = XLALSimInspiralSetQuadMonParamsFromLambdas(params);
        XLAL_CHECK(XLAL_SUCCESS == ret, ret, "Failed to set QuadMon from Lambdas for SEOBNRv5_ROM_NRTidal");
        ret = XLALSimIMRSEOBNRv5ROMNRTidal(hptilde, hctilde, phiRef, deltaF, f_min, f_max, f_ref, distance, inclination, m1, m2, S1z, S2z, lambda1, lambda2, params, NRTidalv3_V);
        break;

    case Lackey_Tidal_2013_SEOBNRv2_ROM:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        ret = XLALSimIMRLackeyTidal2013(hptilde, hctilde, phiRef, deltaF, f_min, f_max, f_ref, distance, inclination, m1, m2, S1z, lambda2);
        break;

    case IMRPhenomP:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFrameAxisIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");      /* Default is LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L : z-axis along direction of orbital angular momentum. */
        if (!XLALSimInspiralWaveformParamsModesChoiceIsDefault( /* Default is (2,2) or l=2 modes. */ params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        /* Tranform to model parameters */
        if (f_ref == 0.0)
            f_ref = f_min;      /* Default reference frequency is minimum frequency */
        XLALSimIMRPhenomPCalculateModelParametersFromSourceFrame(&chi1_l, &chi2_l, &chip, &thetaJN, &alpha0, &phi_aligned, &zeta_polariz, m1, m2, f_ref, phiRef, inclination, S1x, S1y, S1z, S2x, S2y, S2z, IMRPhenomPv1_V);
        /* Call the waveform driver routine */
        ret = XLALSimIMRPhenomP(hptilde, hctilde, chi1_l, chi2_l, chip, thetaJN, m1, m2, distance, alpha0, phi_aligned, deltaF, f_min, f_max, f_ref, IMRPhenomPv1_V, NoNRT_V, params);
        if (ret == XLAL_FAILURE)
            XLAL_ERROR(XLAL_EFUNC);
        for (UINT4 idx = 0; idx < (*hptilde)->data->length; idx++) {
            PhPpolp = (*hptilde)->data->data[idx];
            PhPpolc = (*hctilde)->data->data[idx];
            (*hptilde)->data->data[idx] = cos(2. * zeta_polariz) * PhPpolp + sin(2. * zeta_polariz) * PhPpolc;
            (*hctilde)->data->data[idx] = cos(2. * zeta_polariz) * PhPpolc - sin(2. * zeta_polariz) * PhPpolp;
        }
        break;

    case IMRPhenomPv2:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFrameAxisIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");      /* Default is LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L : z-axis along direction of orbital angular momentum. */
        if (!XLALSimInspiralWaveformParamsModesChoiceIsDefault( /* Default is (2,2) or l=2 modes. */ params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        /* Tranform to model parameters */
        if (f_ref == 0.0)
            f_ref = f_min;      /* Default reference frequency is minimum frequency */
        XLALSimIMRPhenomPCalculateModelParametersFromSourceFrame(&chi1_l, &chi2_l, &chip, &thetaJN, &alpha0, &phi_aligned, &zeta_polariz, m1, m2, f_ref, phiRef, inclination, S1x, S1y, S1z, S2x, S2y, S2z, IMRPhenomPv2_V);
        /* Call the waveform driver routine */
        ret = XLALSimIMRPhenomP(hptilde, hctilde, chi1_l, chi2_l, chip, thetaJN, m1, m2, distance, alpha0, phi_aligned, deltaF, f_min, f_max, f_ref, IMRPhenomPv2_V, NoNRT_V, params);
        if (ret == XLAL_FAILURE)
            XLAL_ERROR(XLAL_EFUNC);
        for (UINT4 idx = 0; idx < (*hptilde)->data->length; idx++) {
            PhPpolp = (*hptilde)->data->data[idx];
            PhPpolc = (*hctilde)->data->data[idx];
            (*hptilde)->data->data[idx] = cos(2. * zeta_polariz) * PhPpolp + sin(2. * zeta_polariz) * PhPpolc;
            (*hctilde)->data->data[idx] = cos(2. * zeta_polariz) * PhPpolc - sin(2. * zeta_polariz) * PhPpolp;
        }
        break;

    case IMRPhenomPv2_NRTidal:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFrameAxisIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");      /* Default is LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L : z-axis along direction of orbital angular momentum. */
        if (!XLALSimInspiralWaveformParamsModesChoiceIsDefault( /* Default is (2,2) or l=2 modes. */ params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        /* Tranform to model parameters */
        if (f_ref == 0.0)
            f_ref = f_min;      /* Default reference frequency is minimum frequency */
        XLALSimIMRPhenomPCalculateModelParametersFromSourceFrame(&chi1_l, &chi2_l, &chip, &thetaJN, &alpha0, &phi_aligned, &zeta_polariz, m1, m2, f_ref, phiRef, inclination, S1x, S1y, S1z, S2x, S2y, S2z, IMRPhenomPv2NRTidal_V);
        /* Call the waveform driver routine */
        ret = XLALSimIMRPhenomP(hptilde, hctilde, chi1_l, chi2_l, chip, thetaJN, m1, m2, distance, alpha0, phi_aligned, deltaF, f_min, f_max, f_ref, IMRPhenomPv2NRTidal_V, NRTidal_V, params);
        if (ret == XLAL_FAILURE)
            XLAL_ERROR(XLAL_EFUNC);
        for (UINT4 idx = 0; idx < (*hptilde)->data->length; idx++) {
            PhPpolp = (*hptilde)->data->data[idx];
            PhPpolc = (*hctilde)->data->data[idx];
            (*hptilde)->data->data[idx] = cos(2. * zeta_polariz) * PhPpolp + sin(2. * zeta_polariz) * PhPpolc;
            (*hctilde)->data->data[idx] = cos(2. * zeta_polariz) * PhPpolc - sin(2. * zeta_polariz) * PhPpolp;
        }
        break;
    case IMRPhenomPv2_NRTidalv2:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFrameAxisIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");      /* Default is LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L : z-axis along direction of orbital angular momentum. */
        if (!XLALSimInspiralWaveformParamsModesChoiceIsDefault( /* Default is (2,2) or l=2 modes. */ params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        /* Tranform to model parameters */
        if (f_ref == 0.0)
            f_ref = f_min;      /* Default reference frequency is minimum frequency */
        XLALSimIMRPhenomPCalculateModelParametersFromSourceFrame(&chi1_l, &chi2_l, &chip, &thetaJN, &alpha0, &phi_aligned, &zeta_polariz, m1, m2, f_ref, phiRef, inclination, S1x, S1y, S1z, S2x, S2y, S2z, IMRPhenomPv2NRTidal_V);
        /* Call the waveform driver routine */
        ret = XLALSimIMRPhenomP(hptilde, hctilde, chi1_l, chi2_l, chip, thetaJN, m1, m2, distance, alpha0, phi_aligned, deltaF, f_min, f_max, f_ref, IMRPhenomPv2NRTidal_V, NRTidalv2_V, params);
        if (ret == XLAL_FAILURE)
            XLAL_ERROR(XLAL_EFUNC);
        for (UINT4 idx = 0; idx < (*hptilde)->data->length; idx++) {
            PhPpolp = (*hptilde)->data->data[idx];
            PhPpolc = (*hctilde)->data->data[idx];
            (*hptilde)->data->data[idx] = cos(2. * zeta_polariz) * PhPpolp + sin(2. * zeta_polariz) * PhPpolc;
            (*hctilde)->data->data[idx] = cos(2. * zeta_polariz) * PhPpolc - sin(2. * zeta_polariz) * PhPpolp;
        }
        break;

    case IMRPhenomPv3:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFrameAxisIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");      /* Default is LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L : z-axis along direction of orbital angular momentum. */
        if (!XLALSimInspiralWaveformParamsModesChoiceIsDefault( /* Default is (2,2) or l=2 modes. */ params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        /* Tranform to model parameters */
        if (f_ref == 0.0)
            f_ref = f_min;      /* Default reference frequency is minimum frequency */
        REAL8Sequence *freqspv3 = XLALCreateREAL8Sequence(2);
        freqspv3->data[0] = f_min;
        freqspv3->data[1] = f_max;
        ret = XLALSimIMRPhenomPv3(hptilde, hctilde, freqspv3, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, deltaF, f_ref, params);
        if (ret == XLAL_FAILURE)
            XLAL_ERROR(XLAL_EFUNC);
        XLALDestroyREAL8Sequence(freqspv3);
        break;

    case IMRPhenomPv3HM:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFrameAxisIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");      /* Default is LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L : z-axis along direction of orbital angular momentum. */
        if (!XLALSimInspiralWaveformParamsModesChoiceIsDefault( /* Default is (2,2) or l=2 modes. */ params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        /* Tranform to model parameters */
        if (f_ref == 0.0)
            f_ref = f_min;      /* Default reference frequency is minimum frequency */

        REAL8Sequence *freqspv3hm = XLALCreateREAL8Sequence(2);
        freqspv3hm->data[0] = f_min;
        freqspv3hm->data[1] = f_max;
        ret = XLALSimIMRPhenomPv3HMGetHplusHcross(hptilde, hctilde, freqspv3hm, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, deltaF, f_ref, params);

        if (ret == XLAL_FAILURE)
            XLAL_ERROR(XLAL_EFUNC);
        XLALDestroyREAL8Sequence(freqspv3hm);

        break;

    case SpinTaylorT4Fourier:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFrameAxisIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
        if (!XLALSimInspiralWaveformParamsModesChoiceIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        spin1x = S1x;
        spin1y = S1y;
        spin1z = S1z;
        spin2x = S2x;
        spin2y = S2y;
        spin2z = S2z;
        ROTATEY(inclination, spin1x, spin1y, spin1z);
        ROTATEY(inclination, spin2x, spin2y, spin2z);
        LNhatx = sin(inclination);
        LNhaty = 0.;
        LNhatz = cos(inclination);
        E1x = 0.;
        E1y = 1.;
        E1z = 0.;
        // default kMax = 3
        kMax = 3;
        // default v0 = 1
        v0 = 1.;
        // default fStart = 0.9*fMin
        fStart = 0.9 * f_min;
        phiRefAtEnd = 0;
        // if f_ref = 0, set it to f_min, and tell the driver routine that we came from there
        if (f_ref == 0) {
            f_ref = f_min;
            phiRefAtEnd = 1;
        }
        // default quadparams are for black holes. Replace by ~2-12 for neutron stars
        /* Call the waveform driver routine */
        ret = XLALSimInspiralSpinTaylorT4Fourier(hptilde, hctilde,
                                                 f_min, f_max, deltaF, kMax, phiRef, v0, m1, m2, fStart, f_ref, distance, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, LNhatx, LNhaty, LNhatz, E1x, E1y, E1z, lambda1, lambda2, quadparam1, quadparam2,
                                                 params, phaseO, amplitudeO, phiRefAtEnd);
        if (ret == XLAL_FAILURE)
            XLAL_ERROR(XLAL_EFUNC);
        break;

    case SpinTaylorT5Fourier:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFrameAxisIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
        if (!XLALSimInspiralWaveformParamsModesChoiceIsDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        spin1x = S1x;
        spin1y = S1y;
        spin1z = S1z;
        spin2x = S2x;
        spin2y = S2y;
        spin2z = S2z;
        ROTATEY(inclination, spin1x, spin1y, spin1z);
        ROTATEY(inclination, spin2x, spin2y, spin2z);
        LNhatx = sin(inclination);
        LNhaty = 0.;
        LNhatz = cos(inclination);
        E1x = 0.;
        E1y = 1.;
        E1z = 0.;
        // default kMax = 3
        kMax = 3;
        // default v0 = 1
        v0 = 1.;
        // default fStart = 0.9*fMin
        fStart = 0.9 * f_min;
        phiRefAtEnd = 0;
        // if f_ref = 0, set it to f_min, and tell the driver routine that we came from there
        if (f_ref == 0) {
            f_ref = f_min;
            phiRefAtEnd = 1;
        }
        // default quadparams are for black holes. Replace by ~2-12 for neutron stars
        /* Call the waveform driver routine */
        ret = XLALSimInspiralSpinTaylorT5Fourier(hptilde, hctilde,
                                                 f_min, f_max, deltaF, kMax, phiRef, v0, m1, m2, fStart, f_ref, distance, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, LNhatx, LNhaty, LNhatz, E1x, E1y, E1z, lambda1, lambda2, quadparam1, quadparam2,
                                                 params, phaseO, amplitudeO, phiRefAtEnd);
        if (ret == XLAL_FAILURE)
            XLAL_ERROR(XLAL_EFUNC);
        break;


    case NRSur4d2s:

        ret = XLALSimNRSur4d2s(hptilde, hctilde, phiRef, deltaF, f_min, f_max, distance, inclination, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z);
        if (ret == XLAL_FAILURE)
            XLAL_ERROR(XLAL_EFUNC);
        break;

    case IMRPhenomXAS:
      {
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");

        /* This is the factor that comes from Y_22star + (-1)^l * Y_2-2 without the dependence in inclination, that is included in pfac and cfac */
        /* Ylm(inclination, beta), with beta = PI/2 - phiRef. phiRef is included in the individual mode */
        COMPLEX16 Ylmfactor = 2.0 * sqrt(5.0 / (64.0 * LAL_PI)) * cexp(-I * 2 * (LAL_PI_2));
        /* The factor for hc is the same but opposite sign */

        /* Call the waveform driver routine. */
        /* It returns h_2-2(f) for positive frequencies. h_2-2 is zero for negative frequencies. */
        /* h_22(f) is zero for positive frequencies. For negatives frequencies h_22(f) = Conjugate[h_2-2(-f)] */
        /* We returns h2_-2 because it is the mode with the positive frequencies,
           and we need this mode because XLALSimInspiralTDFromFD assumes that the input array is in the positive frequencies regime. */
        ret = XLALSimIMRPhenomXASGenerateFD(hptilde, m1, m2, S1z, S2z, distance, f_min, f_max, deltaF, phiRef, f_ref, params);
        if (ret == XLAL_FAILURE)
            XLAL_ERROR(XLAL_EFUNC);

        /* Produce both polarizations for positive frequencies */
        *hctilde = XLALCreateCOMPLEX16FrequencySeries("FD hcross", &((*hptilde)->epoch), (*hptilde)->f0, (*hptilde)->deltaF, &((*hptilde)->sampleUnits), (*hptilde)->data->length);
        for (j = 0; j < (*hptilde)->data->length; j++) {
            (*hctilde)->data->data[j] = -I * cfac * (*hptilde)->data->data[j] * Ylmfactor;
            (*hptilde)->data->data[j] *= pfac * Ylmfactor;
        }
        break;
      }

    case IMRPhenomXAS_NRTidalv2:
        {
        /* Waveform-specific sanity checks */
        if( !XLALSimInspiralWaveformParamsFlagsAreDefault(params) )
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        
        LALDict *LALparams_aux;
                if (params == NULL)
                            LALparams_aux = XLALCreateDict();
                else
                            LALparams_aux = XLALDictDuplicate(params);
                        
        /* PhenomXTidalFlag maps to the NRTidalv* version used: 2-> NRTidalv2 */
        XLALSimInspiralWaveformParamsInsertPhenomXTidalFlag(LALparams_aux,2);
        
        /* to employ multibanding, we need to call XHM: hence, we activate a reduced mode array */
        LALValue *ModeArray = XLALSimInspiralCreateModeArray();
        XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 2);
        XLALSimInspiralModeArrayActivateMode(ModeArray, 2, -2);
        XLALSimInspiralWaveformParamsInsertModeArray(LALparams_aux, ModeArray);
                                    
        ret = XLALSimIMRPhenomXHM(hptilde, hctilde, m1, m2, S1z, S2z, f_min, f_max, deltaF, distance, inclination, phiRef, f_ref, LALparams_aux);
                            
        XLALDestroyValue(ModeArray);
        
        if (ret == XLAL_FAILURE)
                        XLAL_ERROR(XLAL_EFUNC);
                    
        XLALDestroyDict(LALparams_aux);
                
            break;
	}

        case IMRPhenomXAS_NRTidalv3:
        {
        /* Waveform-specific sanity checks */
        if( !XLALSimInspiralWaveformParamsFlagsAreDefault(params) )
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        
        LALDict *LALparams_aux;
                if (params == NULL)
                            LALparams_aux = XLALCreateDict();
                else
                            LALparams_aux = XLALDictDuplicate(params);
                        
        /* PhenomXTidalFlag maps to the NRTidalv* version used: 3-> NRTidalv3 */
        XLALSimInspiralWaveformParamsInsertPhenomXTidalFlag(LALparams_aux,3);
        
        /* to employ multibanding, we need to call XHM: hence, we activate a reduced mode array */
        LALValue *ModeArray = XLALSimInspiralCreateModeArray();
        XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 2);
        XLALSimInspiralModeArrayActivateMode(ModeArray, 2, -2);
        XLALSimInspiralWaveformParamsInsertModeArray(LALparams_aux, ModeArray);
                                    
        ret = XLALSimIMRPhenomXHM(hptilde, hctilde, m1, m2, S1z, S2z, f_min, f_max, deltaF, distance, inclination, phiRef, f_ref, LALparams_aux);
                            
        XLALDestroyValue(ModeArray);
        
        if (ret == XLAL_FAILURE)
                        XLAL_ERROR(XLAL_EFUNC);
                    
        XLALDestroyDict(LALparams_aux);
                
            break;
	}

    case IMRPhenomXHM:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");

        /* Activate or not the debug info. */
#ifndef PHENOMXHMDEBUG          // Defined at compilation time. ./configure --prefix=... CFLAGS='-g -D PHENOMXHMDEBUG'
#define DEBUG 0
#else
#define DEBUG 1                 //print debugging info
#endif

        /* Return hp and hc for positive frequencies. Only negatives modes contribute to positive frequencies.  */
        /* The negative frquencies contribution is the complex conjugate of the positive one. */

        /* Take input/default value for the threshold of the Multibanding. If = 0 then do not use Multibanding. */
        REAL8 resTest = XLALSimInspiralWaveformParamsLookupPhenomXHMThresholdMband(params);

        /* If total mass is very high (>500 solar masses), we only have a few points in the ringdown, interpolation is not efficient, do not use Multibanding */
        REAL8 Mtot = (m1 + m2) / LAL_MSUN_SI;
        if (resTest != 0 && Mtot > 500) {
            resTest = 0.;
        }

        if (resTest == 0.) {    //Do not use multibanding
            ret = XLALSimIMRPhenomXHM2(hptilde, hctilde, m1, m2, S1z, S2z, f_min, f_max, deltaF, distance, inclination, phiRef, f_ref, params);
        } else {                // Use multibanding
            ret = XLALSimIMRPhenomXHM(hptilde, hctilde, m1, m2, S1z, S2z, f_min, f_max, deltaF, distance, inclination, phiRef, f_ref, params);
        }

        if (ret == XLAL_FAILURE)
            XLAL_ERROR(XLAL_EFUNC);

#if DEBUG == 1
        printf("\n\n**********Leaving ChooseFDWaveform *********************\n\n");
#endif

        break;

    case IMRPhenomXP:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFrameAxisIsDefault(params)) {
            /* Default is LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L : z-axis along direction of orbital angular momentum. */
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
        }
        if (!XLALSimInspiralWaveformParamsModesChoiceIsDefault(params)) {
            /* Default is (2,2) or l=2 modes. */
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        }
        if (!checkTidesZero(lambda1, lambda2)) {
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        }
        if (f_ref == 0.0) {
            /* Default reference frequency is minimum frequency */
            f_ref = f_min;
        }

        /* Call the main waveform driver. Note that we pass the full spin vectors
           with XLALSimIMRPhenomXPCalculateModelParametersFromSourceFrame being
           effectively called in the initialization of the pPrec struct
         */
        ret = XLALSimIMRPhenomXPGenerateFD(hptilde, hctilde, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, f_min, f_max, deltaF, f_ref, params);
        if (ret == XLAL_FAILURE) {
            XLAL_ERROR(XLAL_EFUNC);
        }

        break;


    case IMRPhenomXP_NRTidalv2:
        
        {
            
        /* Waveform-specific sanity checks */
        if( !XLALSimInspiralWaveformParamsFrameAxisIsDefault(params) )
        {
            /* Default is LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L : z-axis along direction of orbital angular momentum. */
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
        }
        if(!XLALSimInspiralWaveformParamsModesChoiceIsDefault(params))
        {
            /* Default is (2,2) or l=2 modes. */
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        }
            
        LALDict *LALparams_aux;
        if (params == NULL) LALparams_aux = XLALCreateDict();
        else LALparams_aux = XLALDictDuplicate(params);
                            
        /* PhenomXTidalFlag maps to the NRTidalv* version used: 2-> NRTidalv2 */
        XLALSimInspiralWaveformParamsInsertPhenomXTidalFlag(LALparams_aux,2);
        
        /* to employ multibanding, we need to call XPHM: hence, we activate a reduced mode array */
        LALValue *ModeArray = XLALSimInspiralCreateModeArray();
        XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 2);
        XLALSimInspiralModeArrayActivateMode(ModeArray, 2, -2);
        XLALSimInspiralWaveformParamsInsertModeArray(LALparams_aux, ModeArray);
        
	    // Insert quadrupole parameters for use in twisting up
        ret = XLALSimInspiralSetQuadMonParamsFromLambdas(LALparams_aux);
        XLAL_CHECK(XLAL_SUCCESS == ret, ret, "Failed to set quadrupole parameters from lambdas for IMRPhenomXP_NRTidalv2");
            
        if(f_ref==0.0)
        {
            /* Default reference frequency is minimum frequency */
            f_ref = f_min;
        }
            
        ret = XLALSimIMRPhenomXPHM(hptilde, hctilde, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination,phiRef, f_min, f_max, deltaF, f_ref, LALparams_aux);
        
        XLALDestroyValue(ModeArray);
    
        if (ret == XLAL_FAILURE)    XLAL_ERROR(XLAL_EFUNC);
                        
        XLALDestroyDict(LALparams_aux);
            
        }
        break;

    case IMRPhenomXP_NRTidalv3:
        
        {
            
        /* Waveform-specific sanity checks */
        if( !XLALSimInspiralWaveformParamsFrameAxisIsDefault(params) )
        {
            /* Default is LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L : z-axis along direction of orbital angular momentum. */
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
        }
        if(!XLALSimInspiralWaveformParamsModesChoiceIsDefault(params))
        {
            /* Default is (2,2) or l=2 modes. */
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        }
            
        LALDict *LALparams_aux;
        if (params == NULL) LALparams_aux = XLALCreateDict();
        else LALparams_aux = XLALDictDuplicate(params);
                            
        /* PhenomXTidalFlag maps to the NRTidalv* version used: 3-> NRTidalv3 */
        XLALSimInspiralWaveformParamsInsertPhenomXTidalFlag(LALparams_aux,3);
        
        /* to employ multibanding, we need to call XPHM: hence, we activate a reduced mode array */
        LALValue *ModeArray = XLALSimInspiralCreateModeArray();
        XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 2);
        XLALSimInspiralModeArrayActivateMode(ModeArray, 2, -2);
        XLALSimInspiralWaveformParamsInsertModeArray(LALparams_aux, ModeArray);
        
	    // Insert quadrupole parameters for use in twisting up
        ret = XLALSimInspiralSetQuadMonParamsFromLambdas(LALparams_aux);
        XLAL_CHECK(XLAL_SUCCESS == ret, ret, "Failed to set quadrupole parameters from lambdas for IMRPhenomXP_NRTidalv3");
            
        if(f_ref==0.0)
        {
            /* Default reference frequency is minimum frequency */
            f_ref = f_min;
        }
            
        ret = XLALSimIMRPhenomXPHM(hptilde, hctilde, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination,phiRef, f_min, f_max, deltaF, f_ref, LALparams_aux);
        
        XLALDestroyValue(ModeArray);
    
        if (ret == XLAL_FAILURE)    XLAL_ERROR(XLAL_EFUNC);
                        
        XLALDestroyDict(LALparams_aux);
            
        }
        break;

    case IMRPhenomXPHM:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFrameAxisIsDefault(params)) {
            /* Default is LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L : z-axis along direction of orbital angular momentum. */
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
        }
        if (!XLALSimInspiralWaveformParamsModesChoiceIsDefault(params)) {
            /* Default is (2,2) or l=2 modes. */
            XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        }
        if (!checkTidesZero(lambda1, lambda2)) {
            XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        }
        if (f_ref == 0.0) {
            /* Default reference frequency is minimum frequency */
            f_ref = f_min;
        }

        /* Call the main waveform driver. Note that we pass the full spin vectors
           with XLALSimIMRPhenomXPCalculateModelParametersFromSourceFrame being
           effectively called in the initialization of the pPrec struct
         */
        INT4 usemodes = XLALSimInspiralWaveformParamsLookupPhenomXPHMUseModes(params);

        if (usemodes == 0) {
            ret = XLALSimIMRPhenomXPHM(hptilde, hctilde, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, f_min, f_max, deltaF, f_ref, params);
        } else {
            ret = XLALSimIMRPhenomXPHMFromModes(hptilde, hctilde, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, distance, inclination, phiRef, f_min, f_max, deltaF, f_ref, params);
        }

        if (ret == XLAL_FAILURE) {
            XLAL_ERROR(XLAL_EFUNC);
        }

        break;
        
    case IMRPhenomXO4a:
    {
        LALDict *params_aux;
        if (params == NULL){
            params_aux = XLALCreateDict();
        }
        else{
            params_aux = XLALDictDuplicate(params);
        }

        /* XO4 uses previous version of XHM */
        XLALSimInspiralWaveformParamsInsertPhenomXHMReleaseVersion(params_aux, 122019);
	
    		/* Waveform-specific sanity checks */
    		if( !XLALSimInspiralWaveformParamsFrameAxisIsDefault(params_aux) )
    		{
    			/* Default is LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L : z-axis along direction of orbital angular momentum. */
    			XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
    		}
    		if(!XLALSimInspiralWaveformParamsModesChoiceIsDefault(params_aux))
    		{
    			/* Default is (2,2) or l=2 modes. */
    			XLAL_ERROR(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
    		}
    		if( !checkTidesZero(lambda1, lambda2) )
    		{
    			XLAL_ERROR(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
    		}
    		if(f_ref==0.0)
    		{
    			/* Default reference frequency is minimum frequency */
    			f_ref = f_min;
    		}

    		/* Call the main waveform driver. Note that we pass the full spin vectors
    			 with XLALSimIMRPhenomXPCalculateModelParametersFromSourceFrame being
    			 effectively called in the initialization of the pPrec struct
    		*/

        /* Toggle on PNR angles */
        if(!XLALDictContains(params_aux, "PNRUseTunedAngles")){
            XLALSimInspiralWaveformParamsInsertPhenomXPNRUseTunedAngles(params_aux, 1);
        }

        /* Toggle on tuned coprecessing strain */
        if(!XLALDictContains(params_aux, "PNRUseTunedCoprec")){
            XLALSimInspiralWaveformParamsInsertPhenomXPNRUseTunedCoprec(params_aux, 1);
        }
        if(!XLALDictContains(params_aux, "PNRForceXHMAlignment")){
            XLALSimInspiralWaveformParamsInsertPhenomXPNRForceXHMAlignment(params_aux, 0);
        }

        /* Toggle on antisymmetric contributions */
        if(!XLALDictContains(params_aux, "AntisymmetricWaveform")){
            XLALSimInspiralWaveformParamsInsertPhenomXAntisymmetricWaveform(params_aux, 1);
        }

        if(XLALSimInspiralWaveformParamsLookupPhenomXAntisymmetricWaveform(params_aux))
        {
            if(!XLALSimInspiralWaveformParamsLookupPhenomXPNRUseTunedAngles(params_aux))
            {
                XLAL_ERROR(XLAL_EFUNC,"Error: Antisymmetric waveform generation not supported without PNR angles, please turn on PNR angles to produce waveform with asymmetries in the (2,2) and (2,-2) modes \n");
            }
        }

        /* Toggle on reviewed PrecVersion and FinalSpinMod */
        if(!XLALDictContains(params_aux, "PrecVersion")){
            XLALSimInspiralWaveformParamsInsertPhenomXPrecVersion(params_aux, 300);
        }

    		usemodes = XLALSimInspiralWaveformParamsLookupPhenomXPHMUseModes(params_aux);

    		if(usemodes == 0){
    			ret = XLALSimIMRPhenomXPHM(
    				hptilde, hctilde,
    				m1, m2,
    				S1x, S1y, S1z,
    				S2x, S2y, S2z,
    				distance, inclination,
    				phiRef, f_min, f_max, deltaF, f_ref, params_aux
    			);
    		}
    		else{
    			ret = XLALSimIMRPhenomXPHMFromModes(
    				hptilde, hctilde,
    				m1, m2,
    				S1x, S1y, S1z,
    				S2x, S2y, S2z,
    				distance, inclination,
    				phiRef, f_min, f_max, deltaF, f_ref, params_aux
    			);
    		}

        XLALDestroyDict(params_aux);


    		if (ret == XLAL_FAILURE)
    		{
    			XLAL_ERROR(XLAL_EFUNC);
    		}

    		break;
    }
    default:
        XLALPrintError("FD version of approximant not implemented in lalsimulation\n");
        XLAL_ERROR(XLAL_EINVAL);
    }

    REAL8 polariz = longAscNodes;
    if (polariz) {
        COMPLEX16 tmpP, tmpC;
        for (UINT4 idx = 0; idx < (*hptilde)->data->length; idx++) {
            tmpP = (*hptilde)->data->data[idx];
            tmpC = (*hctilde)->data->data[idx];
            (*hptilde)->data->data[idx] = cos(2. * polariz) * tmpP + sin(2. * polariz) * tmpC;
            (*hctilde)->data->data[idx] = cos(2. * polariz) * tmpC - sin(2. * polariz) * tmpP;
        }
    }

    if (ret == XLAL_FAILURE)
        XLAL_ERROR(XLAL_EFUNC);
    if (XLALSimInspiralWaveformParamsLookupEnableLIV(params))
        ret = XLALSimLorentzInvarianceViolationTerm(hptilde, hctilde, m1 / LAL_MSUN_SI, m2 / LAL_MSUN_SI, distance, params);
    if (ret == XLAL_FAILURE)
        XLAL_ERROR(XLAL_EFUNC);

    return ret;
}

/**
 * Copy of the old code of XLALSimInspiralChooseTDModes(). The new version of XLALSimInspiralChooseTDModes() is just a wrapper over XLALSimInspiralGenerateTDModes().
 * XLALSimInspiralGenerateTDModes() internally calls this function for legacy approximants.
 */
static SphHarmTimeSeries *XLALSimInspiralChooseTDModes_legacy(
  REAL8 phiRef,                        /* reference orbital phase (rad). This variable is not used and only kept here for backwards compatibility */
  REAL8 deltaT,                               /* sampling interval (s) */
  REAL8 m1,                                   /* mass of companion 1 (kg) */
  REAL8 m2,                                   /* mass of companion 2 (kg) */
  REAL8 S1x,                                  /* x-component of the dimensionless spin of object 1 */
  REAL8 S1y,                                  /* y-component of the dimensionless spin of object 1 */
  REAL8 S1z,                                  /* z-component of the dimensionless spin of object 1 */
  REAL8 S2x,                                  /* x-component of the dimensionless spin of object 2 */
  REAL8 S2y,                                  /* y-component of the dimensionless spin of object 2 */
  REAL8 S2z,                                  /* z-component of the dimensionless spin of object 2 */
  REAL8 f_min,                                /* starting GW frequency (Hz) */
  REAL8 f_ref,                                /* reference GW frequency (Hz) */
  REAL8 r,                                    /* distance of source (m) */
  LALDict *LALpars,                           /* LAL dictionary containing accessory parameters */
  int lmax,                                   /* generate all modes with l <= lmax */
  Approximant approximant                     /* post-Newtonian approximant to use for waveform production */
)
{
    XLALPrintWarning("WARNING: The phiRef argument in XLALSimInspiralChooseTDModes will be removed in the future and is currently not used. \n");
    REAL8 v0 = 1.;
    SphHarmTimeSeries *hlm = NULL;
    INT4 errCode = 0;

    /* SEOBNR flag for precessing model version. 3 for SEOBNRv3, 300 for SEOBNRv3_opt, 401 for SEOBNRv4P, 402 for SEOBNRv4PHM */
    UINT4 PrecEOBversion;
    REAL8 spin1[3], spin2[3];

    /* General sanity checks that will abort */
    /*
     * If non-GR approximants are added, change the below to
     * if( nonGRparams && approximant != nonGR1 && approximant != nonGR2 )
     */
    if (!XLALSimInspiralWaveformParamsNonGRAreDefault(LALpars)) {
        XLALPrintError("XLAL Error - %s: Passed in non-NULL pointer to LALSimInspiralTestGRParam for an approximant that does not use LALSimInspiralTestGRParam\n", __func__);
        XLAL_ERROR_NULL(XLAL_EINVAL);
    }

    /* General sanity check the input parameters - only give warnings! */
    if (deltaT > 1.)
        XLALPrintWarning("XLAL Warning - %s: Large value of deltaT = %e requested.\nPerhaps sample rate and time step size were swapped?\n", __func__, deltaT);
    if (deltaT < 1. / 16385.)
        XLALPrintWarning("XLAL Warning - %s: Small value of deltaT = %e requested.\nCheck for errors, this could create very large time series.\n", __func__, deltaT);
    if (m1 < 0.09 * LAL_MSUN_SI)
        XLALPrintWarning("XLAL Warning - %s: Small value of m1 = %e (kg) = %e (Msun) requested.\nPerhaps you have a unit conversion error?\n", __func__, m1, m1 / LAL_MSUN_SI);
    if (m2 < 0.09 * LAL_MSUN_SI)
        XLALPrintWarning("XLAL Warning - %s: Small value of m2 = %e (kg) = %e (Msun) requested.\nPerhaps you have a unit conversion error?\n", __func__, m2, m2 / LAL_MSUN_SI);
    if (m1 + m2 > 1000. * LAL_MSUN_SI)
        XLALPrintWarning("XLAL Warning - %s: Large value of total mass m1+m2 = %e (kg) = %e (Msun) requested.\nSignal not likely to be in band of ground-based detectors.\n", __func__, m1 + m2, (m1 + m2) / LAL_MSUN_SI);
    if (S1x * S1x + S1y * S1y + S1z * S1z > 1.000001)
        XLALPrintWarning("XLAL Warning - %s: S1 = (%e,%e,%e) with norm > 1 requested.\nAre you sure you want to violate the Kerr bound?\n", __func__, S1x, S1y, S1z);
    if (S2x * S2x + S2y * S2y + S2z * S2z > 1.000001)
        XLALPrintWarning("XLAL Warning - %s: S2 = (%e,%e,%e) with norm > 1 requested.\nAre you sure you want to violate the Kerr bound?\n", __func__, S2x, S2y, S2z);
    if (f_min < 1.)
        XLALPrintWarning("XLAL Warning - %s: Small value of fmin = %e requested.\nCheck for errors, this could create a very long waveform.\n", __func__, f_min);
    if (f_min > 40.000001)
        XLALPrintWarning("XLAL Warning - %s: Large value of fmin = %e requested.\nCheck for errors, the signal will start in band.\n", __func__, f_min);

    REAL8 lambda1 = XLALSimInspiralWaveformParamsLookupTidalLambda1(LALpars);
    REAL8 lambda2 = XLALSimInspiralWaveformParamsLookupTidalLambda2(LALpars);
    int amplitudeO = XLALSimInspiralWaveformParamsLookupPNAmplitudeOrder(LALpars);
    int phaseO = XLALSimInspiralWaveformParamsLookupPNPhaseOrder(LALpars);
    UINT4 l;

    switch (approximant) {
    case TaylorT1:
        /* Waveform-specific sanity checks */
        if (!checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-zero spins were given, but this is a non-spinning approximant.");
        if (!XLALSimInspiralWaveformParamsFrameAxisIsDefault(LALpars))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
        if (!XLALSimInspiralWaveformParamsModesChoiceIsDefault(LALpars))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        /* Call the waveform driver routine */
        hlm = XLALSimInspiralTaylorT1PNModes(v0, deltaT, m1, m2, f_min, f_ref, r, lambda1, lambda2, XLALSimInspiralWaveformParamsLookupPNTidalOrder(LALpars), amplitudeO, phaseO, lmax);
        break;
    case TaylorT2:
        /* Waveform-specific sanity checks */
        if (!checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-zero spins were given, but this is a non-spinning approximant.");
        if (!XLALSimInspiralWaveformParamsFrameAxisIsDefault(LALpars))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
        if (!XLALSimInspiralWaveformParamsModesChoiceIsDefault(LALpars))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        /* Call the waveform driver routine */
        hlm = XLALSimInspiralTaylorT2PNModes(v0, deltaT, m1, m2, f_min, f_ref, r, lambda1, lambda2, XLALSimInspiralWaveformParamsLookupPNTidalOrder(LALpars), amplitudeO, phaseO, lmax);
        break;
    case TaylorT3:
        /* Waveform-specific sanity checks */
        if (!checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-zero spins were given, but this is a non-spinning approximant.");
        if (!XLALSimInspiralWaveformParamsFrameAxisIsDefault(LALpars))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
        if (!XLALSimInspiralWaveformParamsModesChoiceIsDefault(LALpars))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        /* Call the waveform driver routine */
        hlm = XLALSimInspiralTaylorT3PNModes(v0, deltaT, m1, m2, f_min, f_ref, r, lambda1, lambda2, 0, amplitudeO, phaseO, lmax);
        break;
    case TaylorT4:
        /* Waveform-specific sanity checks */
        if (!checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-zero spins were given, but this is a non-spinning approximant.");
        if (!XLALSimInspiralWaveformParamsFrameAxisIsDefault(LALpars))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
        if (!XLALSimInspiralWaveformParamsModesChoiceIsDefault(LALpars))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        /* Call the waveform driver routine */
        hlm = XLALSimInspiralTaylorT4PNModes(v0, deltaT, m1, m2, f_min, f_ref, r, lambda1, lambda2, 0, amplitudeO, phaseO, lmax);
        break;
    case EOBNRv2:
    case EOBNRv2HM:
        /* Waveform-specific sanity checks */
        if (!checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-zero spins were given, but this is a non-spinning approximant.");
        if (!XLALSimInspiralWaveformParamsFrameAxisIsDefault(LALpars))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.");
        if (!XLALSimInspiralWaveformParamsModesChoiceIsDefault(LALpars))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.");
        /* Call the waveform driver routine */
        hlm = XLALSimIMREOBNRv2Modes(deltaT, m1, m2, f_min, r);
        // EOB driver only outputs modes with m>0, add m<0 modes by symmetry
        size_t j;
        int m;
        for (l = 2; l <= XLALSphHarmTimeSeriesGetMaxL(hlm); l++) {
            for (m = -l; m < 0; m++) {
                COMPLEX16TimeSeries *inmode = XLALSphHarmTimeSeriesGetMode(hlm, l, -m);
                if (!inmode)
                    continue;
                COMPLEX16TimeSeries *tmpmode = XLALCutCOMPLEX16TimeSeries(inmode, 0, inmode->data->length);
                for (j = 0; j < tmpmode->data->length; j++) {
                    tmpmode->data->data[j] = cpow(-1, l)
                        * conj(tmpmode->data->data[j]);
                }
                hlm = XLALSphHarmTimeSeriesAddMode(hlm, tmpmode, l, m);
                XLALDestroyCOMPLEX16TimeSeries(tmpmode);
            }
        }
        break;

    case NRSur7dq2:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(LALpars))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        /* Call the waveform driver routine */
        hlm = XLALSimInspiralPrecessingNRSurModes(deltaT, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref, r, LALpars, approximant);
        break;

    case NRSur7dq4:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(LALpars))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        /* Call the waveform driver routine */
        hlm = XLALSimInspiralPrecessingNRSurModes(deltaT, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref, r, LALpars, approximant);
        break;

    case NRHybSur3dq8:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(LALpars))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");

        /* Call the waveform driver routine */
        hlm = XLALSimIMRNRHybSur3dq8Modes(deltaT, m1, m2, S1z, S2z, f_min, f_ref, r, LALpars);
        break;

    case IMRPhenomTHM:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(LALpars))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");

        /* Call the waveform driver routine */
        hlm = XLALSimIMRPhenomTHM_Modes(m1, m2, S1z, S2z, r, deltaT, f_min, f_ref, phiRef, LALpars);

        break;

    case IMRPhenomTPHM:
        /* Waveform-specific sanity checks */
        /* FIXME: CHECK XPHM CHECKS */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(LALpars))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        /* Call the waveform driver routine. */
        hlm = XLALSimIMRPhenomTPHM_ChooseTDModes(m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, r, deltaT, f_min, f_ref, LALpars);
        break;
    
    case SEOBNRv4HM_PA:
    case pSEOBNRv4HM_PA:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(LALpars))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        if (f_ref != 0.)
            XLALPrintWarning("XLAL Warning - %s: This approximant does not use f_ref. The reference phase will be defined at coalescence.\n", __func__);

        UINT4 SpinAlignedEOBversion;
        if (approximant == SEOBNRv4HM_PA) SpinAlignedEOBversion = 4111;
        if (approximant == pSEOBNRv4HM_PA) SpinAlignedEOBversion = 4112;

        REAL8Vector *dynamics = NULL;
        REAL8Vector *dynamicsHi = NULL;

        REAL8 lambda2Tidal1 = 0.0;
        REAL8 lambda2Tidal2 = 0.0;
        REAL8 omega02Tidal1 = 0.0;
        REAL8 omega02Tidal2 = 0.0;
        REAL8 lambda3Tidal1 = 0.0;
        REAL8 lambda3Tidal2 = 0.0;
        REAL8 omega03Tidal1 = 0.0;
        REAL8 omega03Tidal2 = 0.0;
        REAL8 quadparam1 = 1.0;
        REAL8 quadparam2 = 1.0;
        REAL8Vector *nqcCoeffsInput = NULL;
        INT4 nqcFlag = 0;

        LALDict *PAParams = XLALCreateDict();
        XLALDictInsertUINT4Value(PAParams, "PAFlag", 1);
        XLALDictInsertUINT4Value(PAParams, "PAOrder", 8);
        XLALDictInsertREAL8Value(PAParams, "rFinal", 1.8);
        XLALDictInsertREAL8Value(PAParams, "rSwitch", 1.8);
        XLALDictInsertUINT2Value(PAParams, "analyticFlag", 1);

        LALDict *TGRParams = XLALCreateDict();
        REAL8 domega220 = 0.0;
        REAL8 dtau220 = 0.0;
        REAL8 domega210 = 0.0;
        REAL8 dtau210 = 0.0;
        REAL8 domega330 = 0.0;
        REAL8 dtau330 = 0.0;
        REAL8 domega440 = 0.0;
        REAL8 dtau440 = 0.0;
        REAL8 domega550 = 0.0;
        REAL8 dtau550 = 0.0;

        domega220 = XLALSimInspiralWaveformParamsLookupDOmega220(LALpars);
        dtau220 = XLALSimInspiralWaveformParamsLookupDTau220(LALpars);
        domega210 = XLALSimInspiralWaveformParamsLookupDOmega210(LALpars);
        dtau210 = XLALSimInspiralWaveformParamsLookupDTau210(LALpars);
        domega330 = XLALSimInspiralWaveformParamsLookupDOmega330(LALpars);
        dtau330 = XLALSimInspiralWaveformParamsLookupDTau330(LALpars);
        domega440 = XLALSimInspiralWaveformParamsLookupDOmega440(LALpars);
        dtau440 = XLALSimInspiralWaveformParamsLookupDTau440(LALpars);
        domega550 = XLALSimInspiralWaveformParamsLookupDOmega550(LALpars);
        dtau550 = XLALSimInspiralWaveformParamsLookupDTau550(LALpars);

        UINT2 TGRflag = 0;
        if (approximant == pSEOBNRv4HM_PA) TGRflag = 1;
        
        XLALSimInspiralWaveformParamsInsertDOmega220(TGRParams, domega220);
        XLALSimInspiralWaveformParamsInsertDTau220(TGRParams, dtau220);
        XLALSimInspiralWaveformParamsInsertDOmega210(TGRParams, domega210);
        XLALSimInspiralWaveformParamsInsertDTau210(TGRParams, dtau210);
        XLALSimInspiralWaveformParamsInsertDOmega330(TGRParams, domega330);
        XLALSimInspiralWaveformParamsInsertDTau330(TGRParams, dtau330);
        XLALSimInspiralWaveformParamsInsertDOmega440(TGRParams, domega440);
        XLALSimInspiralWaveformParamsInsertDTau440(TGRParams, dtau440);
        XLALSimInspiralWaveformParamsInsertDOmega550(TGRParams, domega550);
        XLALSimInspiralWaveformParamsInsertDTau550(TGRParams, dtau550);

        XLALDictInsertUINT2Value(TGRParams, "TGRflag", TGRflag);

        if(XLALSimIMRSpinAlignedEOBModes (
            &hlm,
            &dynamics, &dynamicsHi,
            deltaT,
            m1, m2,
            f_min,
            r,
            S1z, S2z,
            SpinAlignedEOBversion,
            lambda2Tidal1, lambda2Tidal2,
            omega02Tidal1, omega02Tidal2,
            lambda3Tidal1, lambda3Tidal2,
            omega03Tidal1, omega03Tidal2,
            quadparam1, quadparam2,
            nqcCoeffsInput, nqcFlag,
            PAParams,
            TGRParams) == XLAL_FAILURE
        ){
            XLAL_ERROR_NULL (XLAL_EFUNC);
        };

        if(dynamics) XLALDestroyREAL8Vector(dynamics);
        if(dynamicsHi) XLALDestroyREAL8Vector(dynamicsHi);
        XLALDestroyDict(PAParams);
        XLALDestroyDict(TGRParams);
        
        UINT4 i;
        UINT4 modeArrayCreated = 0;

        LALValue *modeArray = XLALSimInspiralWaveformParamsLookupModeArray(
            LALpars
        );
        
            if (modeArray == NULL) {
                modeArray = XLALSimInspiralCreateModeArray();
                modeArrayCreated = 1;

                XLALSimInspiralModeArrayActivateMode(modeArray, 2, 2);
                XLALSimInspiralModeArrayActivateMode(modeArray, 2, 1);
                XLALSimInspiralModeArrayActivateMode(modeArray, 3, 3);
                XLALSimInspiralModeArrayActivateMode(modeArray, 4, 4);
                XLALSimInspiralModeArrayActivateMode(modeArray, 5, 5);
        }

        SphHarmTimeSeries *modes = hlm;
        COMPLEX16TimeSeries *tmpMode = NULL;
        char modeName[40];

        while (modes) {
            if (XLALSimInspiralModeArrayIsModeActive(
                modeArray, modes->l, -modes->m
            ) == 1) {
                sprintf(modeName, "h%dm%d", modes->l, modes->m);

                tmpMode = XLALCreateCOMPLEX16TimeSeries(
                    modeName,
                    &modes->mode->epoch,
                    0,
                    deltaT,
                    &lalStrainUnit,
                    modes->mode->data->length
                );

                for (i = 0; i < modes->mode->data->length; i++) {
                    tmpMode->data->data[i] = pow(-1, modes->l) * conj(- modes->mode->data->data[i]);
                }

                hlm = XLALSphHarmTimeSeriesAddMode(hlm, tmpMode, modes->l, -(modes->m));
            }

            if (XLALSimInspiralModeArrayIsModeActive(
                modeArray, modes->l, modes->m
            ) == 1) {
                
                for (i = 0; i < modes->mode->data->length; i++)
                    modes->mode->data->data[i] *= -1;
            }
            
            modes = modes->next;
        }

        if (modeArrayCreated) {
            XLALDestroyValue(modeArray);
        }
        
        XLALDestroyCOMPLEX16TimeSeries(tmpMode);

        break;

    case SEOBNRv4P:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(LALpars))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        if (f_ref != 0.)
            XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);

        spin1[0] = S1x;
        spin1[1] = S1y;
        spin1[2] = S1z;
        spin2[0] = S2x;
        spin2[1] = S2y;
        spin2[2] = S2z;
        PrecEOBversion = 401;
        hlm = XLALSimIMRSpinPrecEOBModes(deltaT, m1, m2, f_min, r, spin1, spin2, PrecEOBversion, LALpars);
        break;

    case SEOBNRv4PHM:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(LALpars))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");
        if (f_ref != 0.)
            XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);

        spin1[0] = S1x;
        spin1[1] = S1y;
        spin1[2] = S1z;
        spin2[0] = S2x;
        spin2[1] = S2y;
        spin2[2] = S2z;
        PrecEOBversion = 402;
        hlm = XLALSimIMRSpinPrecEOBModes(deltaT, m1, m2, f_min, r, spin1, spin2, PrecEOBversion, LALpars);
        break;

    case SpinTaylorT1:
    case SpinTaylorT5:
    case SpinTaylorT4:
        if (lmax > 4)
            XLALPrintError("XLAL ERROR - %s: maximum l implemented for SpinTaylors is 4, = %d requested.\n", __func__, lmax);

        REAL8TimeSeries *V = NULL;
        REAL8TimeSeries *Phi = NULL;
        REAL8TimeSeries *Spin1x = NULL;
        REAL8TimeSeries *Spin1y = NULL;
        REAL8TimeSeries *Spin1z = NULL;
        REAL8TimeSeries *Spin2x = NULL;
        REAL8TimeSeries *Spin2y = NULL;
        REAL8TimeSeries *Spin2z = NULL;
        REAL8TimeSeries *LNhx = NULL;
        REAL8TimeSeries *LNhy = NULL;
        REAL8TimeSeries *LNhz = NULL;
        REAL8TimeSeries *E1x = NULL;
        REAL8TimeSeries *E1y = NULL;
        REAL8TimeSeries *E1z = NULL;

        /* Here we start dynamics with L//z and e1//x
         * which is not the standard case for SpinTaylor
         */
        REAL8 lnhx = 0.;
        REAL8 lnhy = 0.;
        REAL8 lnhz = 1.;
        REAL8 e1x = 1.;
        REAL8 e1y = 0.;
        REAL8 e1z = 0.;
        //phi_ref is added later
        errCode +=
            XLALSimInspiralSpinTaylorDriver(NULL, NULL, &V, &Phi, &Spin1x, &Spin1y, &Spin1z, &Spin2x, &Spin2y, &Spin2z, &LNhx, &LNhy, &LNhz, &E1x, &E1y, &E1z, 0., deltaT, m1, m2, f_min, f_ref, r, S1x, S1y, S1z, S2x, S2y, S2z, lnhx, lnhy, lnhz, e1x, e1y,
                                            e1z, LALpars, approximant);
        INT4 ma_needs_destroy = 0;
        LALValue *modearray = XLALSimInspiralWaveformParamsLookupModeArray(LALpars);
        if (modearray == NULL) {
            modearray = XLALSimInspiralCreateModeArray();
            ma_needs_destroy = 1;
            for (l = 2; l <= (UINT4) lmax; l++)
                XLALSimInspiralModeArrayActivateAllModesAtL(modearray, l);
        }
        errCode += XLALSimInspiralSpinTaylorHlmModesFromOrbit(&hlm, V, Phi, LNhx, LNhy, LNhz, E1x, E1y, E1z, Spin1x, Spin1y, Spin1z, Spin2x, Spin2y, Spin2z, m1, m2, r, XLALSimInspiralWaveformParamsLookupPNAmplitudeOrder(LALpars), modearray);

        XLALDestroyREAL8TimeSeries(V);
        XLALDestroyREAL8TimeSeries(Phi);
        XLALDestroyREAL8TimeSeries(Spin1x);
        XLALDestroyREAL8TimeSeries(Spin1y);
        XLALDestroyREAL8TimeSeries(Spin1z);
        XLALDestroyREAL8TimeSeries(Spin2x);
        XLALDestroyREAL8TimeSeries(Spin2y);
        XLALDestroyREAL8TimeSeries(Spin2z);
        XLALDestroyREAL8TimeSeries(LNhx);
        XLALDestroyREAL8TimeSeries(LNhy);
        XLALDestroyREAL8TimeSeries(LNhz);
        XLALDestroyREAL8TimeSeries(E1x);
        XLALDestroyREAL8TimeSeries(E1y);
        XLALDestroyREAL8TimeSeries(E1z);
        if (ma_needs_destroy)
            XLALDestroyValue(modearray);
        break;

    default:
        XLALPrintError("Cannot generate modes for this approximant\n");
        XLAL_ERROR_NULL(XLAL_EINVAL);
    }
    if (errCode || !(hlm))
        XLAL_ERROR_NULL(XLAL_EFUNC);

    return hlm;
}

/**
 * Copy of the old code of XLALSimInspiralChooseFDModes(). The new version of XLALSimInspiralChooseFDModes() is just a wrapper over XLALSimInspiralGenerateFDModes().
 * XLALSimInspiralGeenrateFDModes() internally calls this function for legacy approximants.
 */
static SphHarmFrequencySeries *XLALSimInspiralChooseFDModes_legacy(
    REAL8 m1,                                   /* mass of companion 1 (kg) */
    REAL8 m2,                                   /* mass of companion 2 (kg) */
    REAL8 S1x,                                  /* x-component of the dimensionless spin of object 1 */
    REAL8 S1y,                                  /* y-component of the dimensionless spin of object 1 */
    REAL8 S1z,                                  /* z-component of the dimensionless spin of object 1 */
    REAL8 S2x,                                  /* x-component of the dimensionless spin of object 2 */
    REAL8 S2y,                                  /* y-component of the dimensionless spin of object 2 */
    REAL8 S2z,                                  /* z-component of the dimensionless spin of object 2 */
    REAL8 deltaF,                               /* sampling interval (s) */
    REAL8 f_min,                                /* starting GW frequency (Hz) */
    REAL8 f_max,                                /* ending GW frequency (Hz) */
    REAL8 f_ref,                                /* reference GW frequency (Hz) */
    REAL8 phiRef,                               /* reference phase (rad) */
    REAL8 distance,                             /* distance of source (m) */
    REAL8 inclination,                          /* inclination of source (rad) */
    LALDict *params,                            /* LAL dictionary containing accessory parameters (optional mode array) */
    Approximant approximant                     /* approximant to use for waveform production */
)
{

    REAL8 lambda1 = XLALSimInspiralWaveformParamsLookupTidalLambda1(params);
    REAL8 lambda2 = XLALSimInspiralWaveformParamsLookupTidalLambda2(params);

    /* General sanity checks that will abort
     *
     * If non-GR approximants are added, include them in
     * XLALSimInspiralApproximantAcceptTestGRParams()
     */
    if (!XLALSimInspiralWaveformParamsNonGRAreDefault(params) && XLALSimInspiralApproximantAcceptTestGRParams(approximant) != LAL_SIM_INSPIRAL_TESTGR_PARAMS) {
        XLALPrintError("XLAL Error - %s: Passed in non-NULL pointer to LALSimInspiralTestGRParam for an approximant that does not use LALSimInspiralTestGRParam\n", __func__);
        XLAL_ERROR_NULL(XLAL_EINVAL);
    }

    /* General sanity check the input parameters - only give warnings! */
    if (deltaF > 1.)
        XLALPrintWarning("XLAL Warning - %s: Large value of deltaF = %e requested...This corresponds to a very short TD signal (with padding). Consider a smaller value.\n", __func__, deltaF);
    if (deltaF < 1. / 4096.)
        XLALPrintWarning("XLAL Warning - %s: Small value of deltaF = %e requested...This corresponds to a very long TD signal. Consider a larger value.\n", __func__, deltaF);
    if (m1 < 0.09 * LAL_MSUN_SI)
        XLALPrintWarning("XLAL Warning - %s: Small value of m1 = %e (kg) = %e (Msun) requested...Perhaps you have a unit conversion error?\n", __func__, m1, m1 / LAL_MSUN_SI);
    if (m2 < 0.09 * LAL_MSUN_SI)
        XLALPrintWarning("XLAL Warning - %s: Small value of m2 = %e (kg) = %e (Msun) requested...Perhaps you have a unit conversion error?\n", __func__, m2, m2 / LAL_MSUN_SI);
    if (m1 + m2 > 1000. * LAL_MSUN_SI)
        XLALPrintWarning("XLAL Warning - %s: Large value of total mass m1+m2 = %e (kg) = %e (Msun) requested...Signal not likely to be in band of ground-based detectors.\n", __func__, m1 + m2, (m1 + m2) / LAL_MSUN_SI);
    if (S1x * S1x + S1y * S1y + S1z * S1z > 1.000001)
        XLALPrintWarning("XLAL Warning - %s: S1 = (%e,%e,%e) with norm > 1 requested...Are you sure you want to violate the Kerr bound?\n", __func__, S1x, S1y, S1z);
    if (S2x * S2x + S2y * S2y + S2z * S2z > 1.000001)
        XLALPrintWarning("XLAL Warning - %s: S2 = (%e,%e,%e) with norm > 1 requested...Are you sure you want to violate the Kerr bound?\n", __func__, S2x, S2y, S2z);
    if (f_min < 1.)
        XLALPrintWarning("XLAL Warning - %s: Small value of fmin = %e requested...Check for errors, this could create a very long waveform.\n", __func__, f_min);
    if (f_min > 40.000001)
        XLALPrintWarning("XLAL Warning - %s: Large value of fmin = %e requested...Check for errors, the signal will start in band.\n", __func__, f_min);

    /* Adjust the reference frequency for certain precessing approximants:
     * if that approximate interprets f_ref==0 to be f_min, set f_ref=f_min;
     * otherwise do nothing */
    FIX_REFERENCE_FREQUENCY(f_ref, f_min, approximant);

    /* Output object, structure with the individual modes required.
       The values of each mode are returned both for positive and negative frequencies to be consistent with the precessing models. */
    SphHarmFrequencySeries *hlms = NULL;

    /* Frequency array of each mode. It will have both positive and negative values. */
    REAL8Sequence *freqsSphH = NULL;


    /* The following variables are only used for PhenomHM and SEOBNRv4HM_ROM since some extra operations are needed for them. */

    /* Input ModeArray. If not specified in the LAL dictionary, it will return all the available modes in the model. */
    LALValue *ModeArray = NULL;
    LALDict *params_aux;
    /* This is an auxiliar, easy to read list with the modes in the ModeArray option.
       E.g. if (2,-2), (3,-3) are activated, the it would be (2, -2, 3, -3). */
    INT2Sequence *modeseq;
    /* Variable for the number of modes in the ModeArray */
    UINT4 nmodes;
    /* Variable for the length of individual modes in half the frequency spectrum. */
    INT4 length;
    /* Auxiliar variable to store the individual modes computed from the internal functions of each model which later we will
       apply some operations to be consistent with LAL conventions. */
    SphHarmFrequencySeries **hlms_tmp = XLALMalloc(sizeof(SphHarmFrequencySeries));
    *hlms_tmp = NULL;

    INT4 retcode;

    switch (approximant) {
    case IMRPhenomXHM:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");

        /* Compute individual modes from IMRPhenomXHM */
        XLALSimIMRPhenomXHMModes(&hlms, m1, m2, S1z, S2z, deltaF, f_min, f_max, f_ref, phiRef, distance, params);
        break;

    case IMRPhenomXPHM:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");

        /* Compute individual modes in the J-frame from IMRPhenomXPHM */
        XLALSimIMRPhenomXPHMModes(&hlms, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, deltaF, f_min, f_max, f_ref, phiRef, distance, inclination, params);
        break;

    case IMRPhenomXO4a:
        if (params == NULL){
            params_aux = XLALCreateDict();
        }
        else{
            params_aux = XLALDictDuplicate(params);
        }

        /* XO4 uses previous version of XHM */
        XLALSimInspiralWaveformParamsInsertPhenomXHMReleaseVersion(params_aux, 122019);
	
  			/* Waveform-specific sanity checks */
  			if( !XLALSimInspiralWaveformParamsFlagsAreDefault(params_aux) )
  					XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
  			if( !checkTidesZero(lambda1, lambda2) )
  					XLAL_ERROR_NULL(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");

        /* Toggle on PNR angles */
        if(!XLALDictContains(params_aux, "PNRUseTunedAngles")){
            XLALSimInspiralWaveformParamsInsertPhenomXPNRUseTunedAngles(params_aux, 1);
        }

        /* Toggle on tuned coprecessing strain */
        if(!XLALDictContains(params_aux, "PNRUseTunedCoprec")){
            XLALSimInspiralWaveformParamsInsertPhenomXPNRUseTunedCoprec(params_aux, 1);
        }

        /* Ensure that 33 tuning is set to preferred value */
        if(!XLALDictContains(params_aux, "PNRUseTunedCoprec33")){
            XLALSimInspiralWaveformParamsInsertPhenomXPNRUseTunedCoprec33(params_aux, 0);
        }
        if(!XLALDictContains(params_aux, "PNRForceXHMAlignment")){
            XLALSimInspiralWaveformParamsInsertPhenomXPNRForceXHMAlignment(params_aux, 0);
        }

        /* Toggle on antisymmetric contributions */
        if(!XLALDictContains(params_aux, "AntisymmetricWaveform")){
            XLALSimInspiralWaveformParamsInsertPhenomXAntisymmetricWaveform(params_aux, 1);
        }

        if(XLALSimInspiralWaveformParamsLookupPhenomXAntisymmetricWaveform(params_aux))
        {
            if(!XLALSimInspiralWaveformParamsLookupPhenomXPNRUseTunedAngles(params_aux))
            {
                XLAL_ERROR_NULL(XLAL_EFUNC,"Error: Antisymmetric waveform generation not supported without PNR angles, please turn on PNR angles to produce waveform with asymmetries in the (2,2) and (2,-2) modes \n");
            }
        }

        /* Toggle on reviewed PrecVersion and FinalSpinMod */
        if(!XLALDictContains(params_aux, "PrecVersion")){
            XLALSimInspiralWaveformParamsInsertPhenomXPrecVersion(params_aux, 300);
        }

	      /* Compute individual modes in the J-frame from IMRPhenomXPHM */
        XLALSimIMRPhenomXPHMModes(&hlms, m1, m2, S1x, S1y, S1z,	S2x, S2y, S2z, deltaF, f_min, f_max, f_ref, phiRef, distance, inclination, params_aux);

        XLALDestroyDict(params_aux);

        break;


    case SEOBNRv4HM_ROM:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");

        /* First we define the mode array of the output SphHarmFrequencySeries.
           Although the user can choose this array, the model computes internally all the modes
           and then we just pick those specified by the user.
           The only exception is when only the 2,-2 mode is required, in such case SEOBNRv4_ROM is called.
         */
        if (params == NULL) {
            params_aux = XLALCreateDict();
        } else {
            params_aux = XLALDictDuplicate(params);
        }
        ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(params_aux);
        if (ModeArray == NULL) {
            /* If not specified, fill array with default modes of IMRPhenomHM */
            ModeArray = XLALSimInspiralCreateModeArray();
            XLALSimInspiralModeArrayActivateMode(ModeArray, 2, -2);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 2, -1);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 3, -3);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 4, -4);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 5, -5);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 2);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 1);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 3, 3);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 4, 4);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 5, 5);

            modeseq = XLALSimInspiralModeArrayReadModes(ModeArray);

            XLALDestroyValue(ModeArray);
            nmodes = modeseq->length / 2;
        } else                  // This is just to avoid killing the kernel when you ask for a mode that is not available.
        {
            modeseq = XLALSimInspiralModeArrayReadModes(ModeArray);
            XLALDestroyValue(ModeArray);
            nmodes = modeseq->length / 2;

            /* Check that there are not unavailable modes. */
            LALValue *DefaultModeArray = XLALSimInspiralCreateModeArray();
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 2, -2);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 2, -1);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 3, -3);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 4, -4);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 5, -5);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 2, 2);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 2, 1);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 3, 3);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 4, 4);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 5, 5);

            for (UINT4 i = 0; i < nmodes; i++) {
                INT2 l, m;
                l = modeseq->data[2 * i];
                m = modeseq->data[2 * i + 1];
                if (XLALSimInspiralModeArrayIsModeActive(DefaultModeArray, l, m) == 0) {
                    XLALDestroyValue(DefaultModeArray);
                    XLALDestroyINT2Sequence(modeseq);
                    XLALFree(hlms_tmp);
                    XLAL_ERROR_NULL(XLAL_EINVAL, "Mode (%i,%i) is not available in SEOBNRv4HM_ROM.\n", l, m);
                }
            }
            XLALDestroyValue(DefaultModeArray);
        }
        XLALDestroyDict(params_aux);

        UINT2 eobmodes = 5;
        if (nmodes == 1 && modeseq->data[0] == 2 && abs(modeseq->data[0]) == 2) {
            eobmodes = 1;       // This will  internally call SEOBNRv4_ROM instead of all the modes, therefore saving time.
        }

        /* Compute individual modes of SEOBNRv4HM_ROM */
        retcode = XLALSimIMRSEOBNRv4HMROM_Modes(hlms_tmp, phiRef, deltaF, f_min, f_max, f_ref, distance, m1, m2, S1z, S2z, -1, eobmodes, 1);
        if (retcode != XLAL_SUCCESS) {
            XLALFree(hlms_tmp);
            XLAL_ERROR_NULL(XLAL_EFUNC);
        }


        /* This is the length of half of the frequency spectrum.
           Later we will resize series to add the negative frequency regime. */
        length = (*hlms_tmp)->mode->data->length - 1;


        /* Loop over modes in the SphHarmFrequencySeries. Resize each mode. */
        for (UINT4 i = 0; i < nmodes; i++) {
            INT2 l, m;
            l = modeseq->data[2 * i];
            m = modeseq->data[2 * i + 1];

            COMPLEX16FrequencySeries *hlm = XLALSphHarmFrequencySeriesGetMode(*hlms_tmp, l, -abs(m));


            if (m < 0) {
                /* Resize series to add the negative frequency regime */
                hlm = XLALResizeCOMPLEX16FrequencySeries(hlm, -length, 2 * length + 1);
            } else {
                /* Use equatorial symmetry to transform negative to positive mode. */
                INT4 minus1l = -1;
                if (l % 2 == 0) {
                    minus1l = 1;
                }
                hlm = XLALResizeCOMPLEX16FrequencySeries(hlm, 0, 2 * length + 1);
                for (INT4 j = 0; j < length; j++) {
                    hlm->data->data[j] = minus1l * conj(hlm->data->data[hlm->data->length - 1 - j]);
                    hlm->data->data[hlm->data->length - 1 - j] = 0.;
                }
            }

            hlms = XLALSphHarmFrequencySeriesAddMode(hlms, hlm, l, m);
        }
        XLALDestroyINT2Sequence(modeseq);
        XLALDestroySphHarmFrequencySeries(*hlms_tmp);

        /* Add frequency array to SphHarmFrequencySeries */
        freqsSphH = XLALCreateREAL8Sequence(2 * length + 1);
        for (INT4 i = -length; i <= length; i++) {
            freqsSphH->data[i + length] = i * deltaF;
        }
        XLALSphHarmFrequencySeriesSetFData(hlms, freqsSphH);
        break;

    case SEOBNRv5_ROM:
        /* Waveform-specific sanity checks */
        if( !XLALSimInspiralWaveformParamsFlagsAreDefault(params) )
                XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                XLAL_ERROR_NULL(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if( !checkTidesZero(lambda1, lambda2) )
                XLAL_ERROR_NULL(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");

        if(params == NULL){
            params_aux = XLALCreateDict();
        }
        else{
            params_aux = XLALDictDuplicate(params);
        }
        ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(params_aux);
        if(ModeArray == NULL)
        {
            /* If not specified, fill array with default modes of SEOBNRv5_ROM */
            ModeArray = XLALSimInspiralCreateModeArray();
            XLALSimInspiralModeArrayActivateMode(ModeArray, 2, -2);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 2);

            modeseq = XLALSimInspiralModeArrayReadModes(ModeArray);

            XLALDestroyValue(ModeArray);
            nmodes = modeseq->length/2;
        }
        else // This is just to avoid killing the kernel when you ask for a mode that is not available.
        {
            modeseq = XLALSimInspiralModeArrayReadModes(ModeArray);
            XLALDestroyValue(ModeArray);
            nmodes = modeseq->length/2;

            /* Check that there are not unavailable modes. */
            LALValue *DefaultModeArray = XLALSimInspiralCreateModeArray();
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 2, -2);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 2, 2);

            for(UINT4 i=0; i<nmodes; i++)
            {
                INT2 l, m;
                l = modeseq->data[2*i];
                m = modeseq->data[2*i+1];
                if(XLALSimInspiralModeArrayIsModeActive(DefaultModeArray, l, m) == 0){
                    XLALDestroyValue(DefaultModeArray);
                    XLALDestroyINT2Sequence(modeseq);
                    XLALFree(hlms_tmp);
                    XLAL_ERROR_NULL(XLAL_EINVAL, "Mode (%i,%i) is not available in SEOBNRv5_ROM.\n", l, m);
                }
            }
            XLALDestroyValue(DefaultModeArray);
        }
        XLALDestroyDict(params_aux);

        UINT2 eobmodesv5 = 1;

        /* Compute individual modes of SEOBNRv5_ROM */
        retcode = XLALSimIMRSEOBNRv5HMROM_Modes(hlms_tmp, phiRef, deltaF, f_min, f_max, f_ref, distance, m1, m2, S1z, S2z, -1, eobmodesv5, true, params, NoNRT_V);
        if( retcode != XLAL_SUCCESS){
            XLALFree(hlms_tmp);
            XLAL_ERROR_NULL(XLAL_EFUNC);
        }


        /* This is the length of half of the frequency spectrum.
           Later we will resize series to add the negative frequency regime. */
        length = (*hlms_tmp)->mode->data->length -1;


        /* Loop over modes in the SphHarmFrequencySeries. Resize each mode. */
        for(UINT4 i=0; i<nmodes; i++)
        {
            INT2 l, m;
            l = modeseq->data[2*i];
            m = modeseq->data[2*i+1];

            COMPLEX16FrequencySeries *hlm = XLALSphHarmFrequencySeriesGetMode(*hlms_tmp, l, -abs(m));


            if(m<0){
                /* Resize series to add the negative frequency regime */
                hlm = XLALResizeCOMPLEX16FrequencySeries(hlm, -length, 2*length+1);
            }
            else{
                /* Use equatorial symmetry to transform negative to positive mode. */
                INT4 minus1l = -1;
                if (l%2 == 0){
                    minus1l = 1;
                }
                hlm = XLALResizeCOMPLEX16FrequencySeries(hlm, 0, 2*length+1);
                for(INT4 j=0; j<length; j++)
                {
                    hlm->data->data[j] = minus1l * conj(hlm->data->data[hlm->data->length -1 - j]);
                    hlm->data->data[hlm->data->length -1 - j] = 0.;
                }
            }

            hlms = XLALSphHarmFrequencySeriesAddMode(hlms, hlm, l, m);
        }
        XLALDestroyINT2Sequence(modeseq);
        XLALDestroySphHarmFrequencySeries(*hlms_tmp);

        /* Add frequency array to SphHarmFrequencySeries */
         freqsSphH = XLALCreateREAL8Sequence(2*length+1);
        for (INT4 i = -length; i<=length; i++)
         {
             freqsSphH->data[i+length] = i*deltaF;
         }
         XLALSphHarmFrequencySeriesSetFData(hlms, freqsSphH);
        break;

    case SEOBNRv5HM_ROM:
        /* Waveform-specific sanity checks */
        if( !XLALSimInspiralWaveformParamsFlagsAreDefault(params) )
                XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                XLAL_ERROR_NULL(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if( !checkTidesZero(lambda1, lambda2) )
                XLAL_ERROR_NULL(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");

        if(params == NULL){
            params_aux = XLALCreateDict();
        }
        else{
            params_aux = XLALDictDuplicate(params);
        }
        ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(params_aux);
        if(ModeArray == NULL)
        {
            /* If not specified, fill array with default modes of SEOBNRv5HM_ROM */
            ModeArray = XLALSimInspiralCreateModeArray();
            XLALSimInspiralModeArrayActivateMode(ModeArray, 2, -2);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 2, -1);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 3, -3);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 4, -4);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 5, -5);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 3, -2);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 4, -3);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 2);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 1);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 3, 3);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 4, 4);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 5, 5);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 3, 2);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 4, 3);


            modeseq = XLALSimInspiralModeArrayReadModes(ModeArray);

            XLALDestroyValue(ModeArray);
            nmodes = modeseq->length/2;
        }
        else // This is just to avoid killing the kernel when you ask for a mode that is not available.
        {
            modeseq = XLALSimInspiralModeArrayReadModes(ModeArray);
            XLALDestroyValue(ModeArray);
            nmodes = modeseq->length/2;

            /* Check that there are not unavailable modes. */
            LALValue *DefaultModeArray = XLALSimInspiralCreateModeArray();
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 2, -2);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 2, -1);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 3, -3);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 4, -4);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 5, -5);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 3, -2);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 4, -3);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 2, 2);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 2, 1);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 3, 3);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 4, 4);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 5, 5);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 3, 2);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 4, 3);

            for(UINT4 i=0; i<nmodes; i++)
            {
                INT2 l, m;
                l = modeseq->data[2*i];
                m = modeseq->data[2*i+1];
                if(XLALSimInspiralModeArrayIsModeActive(DefaultModeArray, l, m) == 0){
                    XLALDestroyValue(DefaultModeArray);
                    XLALDestroyINT2Sequence(modeseq);
                    XLALFree(hlms_tmp);
                    XLAL_ERROR_NULL(XLAL_EINVAL, "Mode (%i,%i) is not available in SEOBNRv5HM_ROM.\n", l, m);
                }
            }
            XLALDestroyValue(DefaultModeArray);
        }
        XLALDestroyDict(params_aux);

        UINT2 eobmodesv5hm = 7;
        if(nmodes == 1 && modeseq->data[0]==2 && abs(modeseq->data[1])==2)
        {
            eobmodesv5hm = 1; // This will  internally call SEOBNRv5_ROM instead of all the modes, therefore saving time.
        }

        /* If asking for a subset of modes in the proper order only those are computed, and not all 7 */
        if (
            nmodes == 2 &&
            modeseq->data[0] == 2 && abs(modeseq->data[1]) == 2 &&
            modeseq->data[2] == 3 && abs(modeseq->data[3]) == 3)
        {   
            eobmodesv5hm = 2;
        }
        if (
            nmodes == 3 &&
            modeseq->data[0] == 2 && abs(modeseq->data[1]) == 2 &&
            modeseq->data[4] == 3 && abs(modeseq->data[5]) == 3 && 
            modeseq->data[2] == 2 && abs(modeseq->data[3]) == 1)
        {   
            eobmodesv5hm = 3;
        }
        if (
            nmodes == 4 &&
            modeseq->data[0] == 2 && abs(modeseq->data[1]) == 2 &&
            modeseq->data[4] == 3 && abs(modeseq->data[5]) == 3 && 
            modeseq->data[2] == 2 && abs(modeseq->data[3]) == 1 && 
            modeseq->data[6] == 4 && abs(modeseq->data[7]) == 4)
        {   
            eobmodesv5hm = 4;
        }
        if (
            nmodes == 5 &&
            modeseq->data[0] == 2 && abs(modeseq->data[1]) == 2 &&
            modeseq->data[4] == 3 && abs(modeseq->data[5]) == 3 && 
            modeseq->data[2] == 2 && abs(modeseq->data[3]) == 1 && 
            modeseq->data[6] == 4 && abs(modeseq->data[7]) == 4 && 
            modeseq->data[8] == 5 && abs(modeseq->data[9]) == 5)
        {
            eobmodesv5hm = 5;
        }
        if (
            nmodes == 6 &&
            modeseq->data[0] == 2 && abs(modeseq->data[1]) == 2 &&
            modeseq->data[4] == 3 && abs(modeseq->data[5]) == 3 && 
            modeseq->data[2] == 2 && abs(modeseq->data[3]) == 1 && 
            modeseq->data[8] == 4 && abs(modeseq->data[9]) == 4 && 
            modeseq->data[10] == 5 && abs(modeseq->data[11]) == 5 && 
            modeseq->data[6] == 3 && abs(modeseq->data[7]) == 2)
        {
            eobmodesv5hm = 6;
        }
        /* Compute individual modes of SEOBNRv5HM_ROM */
        retcode = XLALSimIMRSEOBNRv5HMROM_Modes(hlms_tmp, phiRef, deltaF, f_min, f_max, f_ref, distance, m1, m2, S1z, S2z, -1, eobmodesv5hm, true, params, NoNRT_V);
        if( retcode != XLAL_SUCCESS){
            XLALFree(hlms_tmp);
            XLAL_ERROR_NULL(XLAL_EFUNC);
        }


        /* This is the length of half of the frequency spectrum.
           Later we will resize series to add the negative frequency regime. */
        length = (*hlms_tmp)->mode->data->length -1;


        /* Loop over modes in the SphHarmFrequencySeries. Resize each mode. */
        for(UINT4 i=0; i<nmodes; i++)
        {
            INT2 l, m;
            l = modeseq->data[2*i];
            m = modeseq->data[2*i+1];

            COMPLEX16FrequencySeries *hlm = XLALSphHarmFrequencySeriesGetMode(*hlms_tmp, l, -abs(m));


            if(m<0){
                /* Resize series to add the negative frequency regime */
                hlm = XLALResizeCOMPLEX16FrequencySeries(hlm, -length, 2*length+1);
            }
            else{
                /* Use equatorial symmetry to transform negative to positive mode. */
                INT4 minus1l = -1;
                if (l%2 == 0){
                    minus1l = 1;
                }
                hlm = XLALResizeCOMPLEX16FrequencySeries(hlm, 0, 2*length+1);
                for(INT4 j=0; j<length; j++)
                {
                    hlm->data->data[j] = minus1l * conj(hlm->data->data[hlm->data->length -1 - j]);
                    hlm->data->data[hlm->data->length -1 - j] = 0.;
                }
            }

            hlms = XLALSphHarmFrequencySeriesAddMode(hlms, hlm, l, m);
        }
        XLALDestroyINT2Sequence(modeseq);
        XLALDestroySphHarmFrequencySeries(*hlms_tmp);

        /* Add frequency array to SphHarmFrequencySeries */
         freqsSphH = XLALCreateREAL8Sequence(2*length+1);
        for (INT4 i = -length; i<=length; i++)
         {
             freqsSphH->data[i+length] = i*deltaF;
         }
         XLALSphHarmFrequencySeriesSetFData(hlms, freqsSphH);
        break;

    case IMRPhenomHM:
        /* Waveform-specific sanity checks */
        if (!XLALSimInspiralWaveformParamsFlagsAreDefault(params))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-default flags given, but this approximant does not support this case.");
        if (!checkTransverseSpinsZero(S1x, S1y, S2x, S2y))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-zero transverse spins were given, but this is a non-precessing approximant.");
        if (!checkTidesZero(lambda1, lambda2))
            XLAL_ERROR_NULL(XLAL_EINVAL, "Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.");


        /* First we define the mode array of the output SphHarmFrequencySeries.
           PhenomHM only computes those modes specified in this array.
           We use an auxiliary LALDictionary params_I */
        if (params == NULL) {
            params_aux = XLALCreateDict();
        } else {
            params_aux = XLALDictDuplicate(params);
        }
        ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(params_aux);
        if (ModeArray == NULL) {
            /* If not specified, fill array with default modes of IMRPhenomHM */
            ModeArray = XLALSimInspiralCreateModeArray();
            XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 2);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 1);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 3, 3);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 3, 2);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 4, 4);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 4, 3);

            XLALSimInspiralWaveformParamsInsertModeArray(params_aux, ModeArray);

            XLALSimInspiralModeArrayActivateMode(ModeArray, 2, -2);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 2, -1);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 3, -3);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 3, -2);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 4, -4);
            XLALSimInspiralModeArrayActivateMode(ModeArray, 4, -3);

            modeseq = XLALSimInspiralModeArrayReadModes(ModeArray);
            nmodes = modeseq->length / 2;
        } else                  // This is to avoid killing the kernel when you ask for a mode that is not available.
        {
            modeseq = XLALSimInspiralModeArrayReadModes(ModeArray);
            nmodes = modeseq->length / 2;

            /* Modes supported by IMRPhenomHM */
            LALValue *DefaultModeArray = XLALSimInspiralCreateModeArray();
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 2, 2);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 2, 1);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 3, 3);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 3, 2);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 4, 4);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 4, 3);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 2, -2);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 2, -1);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 3, -3);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 3, -2);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 4, -4);
            XLALSimInspiralModeArrayActivateMode(DefaultModeArray, 4, -3);

            /* Check if there is any mode not supported */
            for (UINT4 i = 0; i < nmodes; i++) {
                INT2 l, m;
                l = modeseq->data[2 * i];
                m = modeseq->data[2 * i + 1];

                if (XLALSimInspiralModeArrayIsModeActive(DefaultModeArray, l, m) == 0) {
                    XLALDestroyValue(ModeArray);
                    XLALDestroyValue(DefaultModeArray);
                    XLALDestroyINT2Sequence(modeseq);
                    XLALFree(hlms_tmp);
                    XLAL_ERROR_NULL(XLAL_EINVAL, "Mode (%i,%i) is not available in IMRPhenomHM.\n", l, m);
                }
                /* For the internal function of IMRPhenomHM we must pass an array only with positive modes */
                if (m < 0) {
                    XLALSimInspiralModeArrayDeactivateMode(ModeArray, l, m);
                    XLALSimInspiralModeArrayActivateMode(ModeArray, l, abs(m));
                }
            }
            XLALSimInspiralWaveformParamsInsertModeArray(params_aux, ModeArray);
            XLALDestroyValue(DefaultModeArray);
        }

        /* Build structure for minimum and maximum frequencies */
        REAL8Sequence *freqs = XLALCreateREAL8Sequence(2);
        freqs->data[0] = f_min;
        freqs->data[1] = f_max;


        /* Call individual modes of PhenomHM */
        retcode = XLALSimIMRPhenomHMGethlmModes(hlms_tmp, freqs, m1, m2, 0., 0., S1z, 0., 0., S2z, phiRef, deltaF, f_ref, params_aux);
        XLALDestroyREAL8Sequence(freqs);
        if (retcode != XLAL_SUCCESS) {
            XLALFree(hlms_tmp);
            XLAL_ERROR_NULL(XLAL_EFUNC);
        }


        /* This is the length of half of the frequency spectrum.
           Later we will resize series to add the negative frequency regime. */
        length = (*hlms_tmp)->mode->data->length - 1;


        /* PhenomHM returns the modes in geometrical units, we need to multiply them by amp0 to obtain physical units. */
        const REAL8 Mtot_Msun = (m1 + m2) / LAL_MSUN_SI;
        const REAL8 amp0 = Mtot_Msun * LAL_MRSUN_SI * Mtot_Msun * LAL_MTSUN_SI / distance;

        /* PhenomHM neglects the LAL convention that the azimuthal angle of the spherical harmonics Ylm is PI/2 - phiRef.
           Here we compesate by this factor so it is consistent with the polarizations construction. */
        COMPLEX16 extra_phase = cexp(-I * (LAL_PI_2 - phiRef));

        /* Loop over modes in the SphHarmFrequencySeries.
           We add the previous factors and resize the series. */
        for (UINT4 i = 0; i < nmodes; i++) {
            INT2 l, m;          // Indexes of mode
            l = modeseq->data[2 * i];
            m = modeseq->data[2 * i + 1];

            /* Get one individual mode.
               Either if m is positive or negative we read the same mode and transform accordingly later. */
            COMPLEX16FrequencySeries *hlm = XLALSphHarmFrequencySeriesGetMode(*hlms_tmp, l, abs(m));

            INT4 minus1l = -1;
            if (l % 2 == 0)
                minus1l = 1;

            /* Incorporate correct units and */
            COMPLEX16 extra_factor_lm = minus1l * amp0 * cpow(extra_phase, m);

            if (m < 0) {
                for (UINT4 j = 0; j < hlm->data->length; j++) {
                    hlm->data->data[j] = hlm->data->data[j] * extra_factor_lm;
                }
                hlm = XLALResizeCOMPLEX16FrequencySeries(hlm, -length, 2 * length + 1);
            } else {
                if (XLALSimInspiralModeArrayIsModeActive(ModeArray, l, -m) == 1) {
                    extra_factor_lm = minus1l;
                } else {
                    extra_factor_lm = minus1l * extra_factor_lm;
                }
                hlm = XLALResizeCOMPLEX16FrequencySeries(hlm, 0, 2 * length + 1);
                for (INT4 j = 0; j < length; j++) {
                    hlm->data->data[j] = conj(hlm->data->data[hlm->data->length - 1 - j]) * extra_factor_lm;
                    hlm->data->data[hlm->data->length - 1 - j] = 0.;
                }
            }

            /* Add the mode to the SphHarmFrequencySeries */
            hlms = XLALSphHarmFrequencySeriesAddMode(hlms, hlm, l, m);

        }
        XLALDestroyINT2Sequence(modeseq);
        XLALDestroySphHarmFrequencySeries(*hlms_tmp);
        XLALDestroyValue(ModeArray);
        XLALDestroyDict(params_aux);

        /* Add frequency array to SphHarmFrequencySeries */
        /* Here we build the whole frequency regime (negative and positive). */
        freqsSphH = XLALCreateREAL8Sequence(hlms->mode->data->length);
        for (INT4 i = -length; i <= length; i++) {
            freqsSphH->data[i + length] = i * deltaF;
        }
        XLALSphHarmFrequencySeriesSetFData(hlms, freqsSphH);

        break;
    default:
        XLALPrintError("XLAL ERROR - %s approximant not supported  by ChooseFDModes.\n", XLALSimInspiralGetStringFromApproximant(approximant));
        XLAL_ERROR_NULL(XLAL_EINVAL);
    }
    XLALFree(hlms_tmp);

    if (!(hlms))
        XLAL_ERROR_NULL(XLAL_EFUNC);

    return hlms;
}
